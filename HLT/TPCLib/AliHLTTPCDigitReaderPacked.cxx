// $Id$

//*************************************************************************
// This file is property of and copyright by the ALICE HLT Project        * 
// ALICE Experiment at CERN, All rights reserved.                         *
//                                                                        *
// Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
//                  Timm Steinbeck <timm@kip.uni-heidelberg.de>           *
//                  Jochen Thaeder <thaeder@kip.uni-heidelberg.de>        *
//                  for The ALICE HLT Project.                            *
//                                                                        *
// Permission to use, copy, modify and distribute this software and its   *
// documentation strictly for non-commercial purposes is hereby granted   *
// without fee, provided that the above copyright notice appears in all   *
// copies and that both the copyright notice and this permission notice   *
// appear in the supporting documentation. The authors make no claims     *
// about the suitability of this software for any purpose. It is          *
// provided "as is" without express or implied warranty.                  *
//*************************************************************************/

/** @file   AliHLTTPCDigitReaderPacked.cxx
    @author Timm Steinbeck, Jochen Thaeder, Matthias Richter, Kenneth Aamodt
    @date   
    @brief  A digit reader implementation for raw data, using the offline
            AliAltroRawStream/AliTPCRawStream.
*/

#if __GNUC__>= 3
using namespace std;
#endif

#include "AliHLTTPCDigitReaderPacked.h"

#include "AliTPCRawStream.h"
#include "AliRawReaderMemory.h"
#include "AliRawDataHeader.h"

#include "AliHLTTPCTransform.h"
#include "AliHLTStdIncludes.h"

ClassImp(AliHLTTPCDigitReaderPacked)

AliHLTTPCDigitReaderPacked::AliHLTTPCDigitReaderPacked()
  :
  fTPCRawStream(NULL),
  fCurrentRow(0),
  fCurrentPad(0),
  fCurrentBin(-1),
  fRowOffset(0),
  fNRows(0),
  fNPads(0),
  fData(NULL),
  fUnsorted(kFALSE),
  fDataBunch(),
  fCurrentChannel(-1),
  fbHaveData(false),
  fCurrentPatch(0)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt


  // Matthias Sep 2008: the pad sorting functionality needs a deep
  // revision of the code
  // I just stumbled over a few awkward implementation
  // - each instance allocates the buffer for sorted data, this is
  //   approx. 8 MByte each
  // - each instance loops to extract the buffer size
  // - the AliTPCRawStream reads the Altro mapping for each instance
  //
  // There is now an instance handling of the buffer and the buffer size
  // is only calculated once The sorting of pads is going to be a common
  // functionality of the DigitReader base class in the future.
  //
  // Same applies to the instance of the AliRawReaderMemory and
  // AliTPCRawStream. Since processing of multiple instances is always
  // sequential, one instance can be used for all readers.

  // get max number of rows
  if (fNMaxRows<0) {
  for (Int_t ii=0; ii < 6; ii++)
      if (AliHLTTPCTransform::GetNRows(ii) > fNMaxRows) 
	  fNMaxRows = AliHLTTPCTransform::GetNRows(ii);
  }

  // get max number of pads
  if (fNMaxPads<0) {
  for (Int_t ii=0; ii < AliHLTTPCTransform::GetNRows();ii++ )
      if (AliHLTTPCTransform::GetNPads(ii) > fNMaxPads) 
	  fNMaxPads = AliHLTTPCTransform::GetNPads(ii);
  }

  // get max number of bins
  if (fNTimeBins<0) {
  fNTimeBins = AliHLTTPCTransform::GetNTimeBins();
  }

  fgObjectCount++;
}

Int_t AliHLTTPCDigitReaderPacked::fNMaxRows=-1;
Int_t AliHLTTPCDigitReaderPacked::fNMaxPads=-1;
Int_t AliHLTTPCDigitReaderPacked::fNTimeBins=-1;
Int_t* AliHLTTPCDigitReaderPacked::fgpFreeBufferInstance=NULL;
Int_t* AliHLTTPCDigitReaderPacked::fgpIssuedBufferInstance=NULL;
AliHLTTPCDigitReaderPacked::AliHLTTPCRawStream* AliHLTTPCDigitReaderPacked::fgpFreeStreamInstance=NULL;
AliHLTTPCDigitReaderPacked::AliHLTTPCRawStream* AliHLTTPCDigitReaderPacked::fgpIssuedStreamInstance=NULL;
Int_t AliHLTTPCDigitReaderPacked::fgObjectCount=0;

AliHLTTPCDigitReaderPacked::~AliHLTTPCDigitReaderPacked()
{
  // see header file for class documentation
  if (fData) ReleaseBufferInstance(fData);
  fData=NULL;

  if (fTPCRawStream) ReleaseRawStreamInstance(fTPCRawStream);
  fTPCRawStream=NULL;

  if (--fgObjectCount==0) {
    if (fgpFreeBufferInstance) delete fgpFreeBufferInstance;
    fgpFreeBufferInstance=NULL;
    if (fgpIssuedBufferInstance) delete fgpIssuedBufferInstance;
    fgpIssuedBufferInstance=NULL;
    if (fgpFreeStreamInstance) delete fgpFreeStreamInstance;
    fgpFreeStreamInstance=NULL;
    if (fgpIssuedStreamInstance) delete fgpIssuedStreamInstance;
    fgpIssuedStreamInstance=NULL;
  }
}

Int_t AliHLTTPCDigitReaderPacked::InitBlock(void* ptr,ULong_t size, Int_t patch, Int_t slice)
{
  // see header file for class documentation
  fTPCRawStream=GetRawStreamInstance();
  if (!fTPCRawStream) return -ENODEV;

  fCurrentPatch=patch;
  
  //get DDL ID in order to tell the memory reader which slice/patch to use
  Int_t DDLid= 0;
  if (patch < 2)
    DDLid = 768 + 2*slice + patch;
  else 
    DDLid = 840 + 4*slice + patch-2;

  fTPCRawStream->SetMemory(DDLid, reinterpret_cast<UChar_t*>( ptr ), size );

  // fCurrentRow always is the row number within a partition, whereas the TPCRawStream
  // counts rows within the inner and outer sector, i.e. 0 to 62 for the inner and
  // 63 to 158 for the outer sector
  // fRowOffset is the offset of the first row of the current partition with respect
  // to inner or outer sector.
  fCurrentRow = 0;
  fCurrentPad = 0;
  fCurrentBin = -1;

  Int_t firstrow=AliHLTTPCTransform::GetFirstRow(patch);
  Int_t lastrow=AliHLTTPCTransform::GetLastRow(patch);
  fNRows = lastrow - firstrow + 1;
  for (Int_t ii=firstrow; ii <= lastrow;ii++ ) {
    if (AliHLTTPCTransform::GetNPads(ii) > fNPads) 
      fNPads = AliHLTTPCTransform::GetNPads(ii);
  }

  Int_t offset=0;
  if (patch > 1) offset =  AliHLTTPCTransform::GetFirstRow( 2 );

  fRowOffset = firstrow - offset;
  firstrow -= offset;
  lastrow  -= offset;

  fbHaveData=false;

  if(!fUnsorted){

  // get the global instance of the array
  fData=GetBufferInstance();
  if (!fData) return -ENOMEM;

  // Init array with -1
  memset( fData, 0xFF, sizeof(Int_t)*(fNMaxRows*fNMaxPads*fNTimeBins) );

  // read data and fill in array
  while( fTPCRawStream->Next()){

      Int_t row = fTPCRawStream->GetRow();
      Int_t pad = fTPCRawStream->GetPad();
      Int_t bin = fTPCRawStream->GetTime();

      if ( row < firstrow || row > lastrow || pad > AliHLTTPCTransform::GetNPads(row + offset) || bin > fNTimeBins){
	HLTFatal("Index out of Range Probably wrong patch! %d - %d", slice, patch);
	if ( row < firstrow || row > lastrow ) 
	  HLTFatal("Row out of Range %d < %d < %d",firstrow, row, lastrow);
	if ( pad > AliHLTTPCTransform::GetNPads(row + offset) ) 
	  HLTFatal("Pad out of Range %d < %d < %d",pad, AliHLTTPCTransform::GetNPads(row + offset));
	if ( bin > fNTimeBins )
	  HLTFatal("Bin out of Range %d < %d < %d",bin, fNTimeBins);
      }
      else {  
	  if ((row-fRowOffset)*fNMaxPads*fNTimeBins+ pad*fNTimeBins + bin >=  fNMaxRows*fNMaxPads*fNTimeBins ) {
	      HLTFatal("Index out of array range PAD=%d ||| ROW=%d ||| BIN=%d ||| OFFSET=%d ||| ROWOFFSET=%d", pad, row, bin, offset, fRowOffset);
	      continue;
	  }
	  else {
	      fData[ (row-fRowOffset)*fNMaxPads*fNTimeBins+ pad*fNTimeBins + bin ] = fTPCRawStream->GetSignal() ;
	  }
      }
  }
  }
  return 0;
}

int AliHLTTPCDigitReaderPacked::Reset()
{
  // see header file for class documentation
  if (fData) ReleaseBufferInstance(fData);
  fData=NULL;
  if (fTPCRawStream) ReleaseRawStreamInstance(fTPCRawStream);
  fTPCRawStream=NULL;
  return 0;
}

Bool_t AliHLTTPCDigitReaderPacked::NextSignal()
{
  // see header file for class documentation
  Bool_t readvalue = kTRUE;

  if(!fUnsorted) {
    if (!fData) return false;
    while (1) {
      fCurrentBin++;
      if (fCurrentBin >= fNTimeBins){
	  fCurrentBin = 0;
	  fCurrentPad++;
     
	  if (fCurrentPad >=fNPads){
	      fCurrentPad = 0;
	      fCurrentRow++;
	      
	      if (fCurrentRow >= fNRows){
		  readvalue = kFALSE;
		  break;
	      }
	  }
      }

      if (fCurrentRow*fNMaxPads*fNTimeBins+ fCurrentPad*fNTimeBins + fCurrentBin >=  fNMaxRows*fNMaxPads*fNTimeBins ) {
	  HLTFatal("Overflow: row=%d pad=%d bin=%d", fCurrentRow, fCurrentPad, fCurrentBin);
	  readvalue = kFALSE;
	  break;
      }

      if (fData[ fCurrentRow*fNMaxPads*fNTimeBins + fCurrentPad*fNTimeBins + fCurrentBin  ] != -1) break;
    }
  } else{
    if ((readvalue = fTPCRawStream->Next())) {
      fCurrentBin=fTPCRawStream->GetTime();
    }
  }

  fbHaveData=readvalue;
  return readvalue;
}

Int_t AliHLTTPCDigitReaderPacked::GetRow()
{
  // see header file for class documentation
  if(!fUnsorted){
  return (fCurrentRow + fRowOffset);
  }
  else{
    return (Int_t) fTPCRawStream->GetRow()-fRowOffset;
  }
}

int AliHLTTPCDigitReaderPacked::GetPad()
{
  // see header file for class documentation
  if(!fUnsorted){
    return fCurrentPad;
  }
  else{
    return fTPCRawStream->GetPad();
  }
}

AliHLTUInt32_t AliHLTTPCDigitReaderPacked::GetAltroBlockHWaddr() const
{
  // see header file for class documentation
  return fTPCRawStream->GetHWAddress();
}

int AliHLTTPCDigitReaderPacked::GetRCUTrailerSize()
{
  // see header file for class documentation
  if(fTPCRawStream){
    return fTPCRawStream->GetRCUTrailerSize();
  }
  return 0;
}

bool AliHLTTPCDigitReaderPacked::GetRCUTrailerData(UChar_t*& trData)
{
  // see header file for class documentation
  if(fTPCRawStream){
    return fTPCRawStream->GetRCUTrailerData(trData);
  }
  return false;
}

Int_t AliHLTTPCDigitReaderPacked::GetSignal()
{ 
  // see header file for class documentation
  if(!fUnsorted){
    // check for validity of fData is in NextSignal, no check at here
    return fData[ fCurrentRow*fNMaxPads*fNTimeBins+ fCurrentPad*fNTimeBins + fCurrentBin ];
  }
  else{
    return fTPCRawStream->GetSignal();
  }
}

Int_t AliHLTTPCDigitReaderPacked::GetTime()
{
  // see header file for class documentation
  return fCurrentBin;
}

Int_t AliHLTTPCDigitReaderPacked::GetTimeOfUnsortedSignal(){
  return fTPCRawStream->GetTime();
}

bool AliHLTTPCDigitReaderPacked::NextChannel()
{
  // see header file for class documentation

  // return true if a channel is available. This is true if
  // 1. the current row position >=0 : a signal has already been read in the stream
  //    but it was not part of the previous bunch
  // 2. the current row position <0 : this is the first invocation at all, read
  //    signal and send result according to the availability
  if(fbHaveData || // data available from the last NextSignal call?
     NextSignal()) { // there is data
    fCurrentChannel=GetAltroBlockHWaddr();
    return true;
  }
  return false;
}

int AliHLTTPCDigitReaderPacked::NextBunch()
{  
  // see header file for class documentation
  if (fCurrentChannel<0) return 0;

  fDataBunch.clear();
  //adding the first signal (will always be the leftover from either NextChannel call or Previous bunch)
  fDataBunch.push_back(GetSignal());

  Int_t prevTime=GetTimeOfUnsortedSignal();
  do{
    if(NextSignal()){
      if((int)GetAltroBlockHWaddr()==fCurrentChannel){//check if there is a change in channel(new row and pad)
	if(prevTime==GetTimeOfUnsortedSignal()+1){//if true means that we have consecutive signals
	  prevTime=GetTimeOfUnsortedSignal();
	  fDataBunch.insert(fDataBunch.begin(), GetSignal());
	}
	else{//end of bunch but not of channel
	  break;
	}
      }
      else{
	fCurrentChannel=-1;
	break;
      }
    }
    else{
      // end of data, but there is one bunch to be completed
      fCurrentChannel=-1;
      break;
    }
  }while(1);

  fCurrentBin=prevTime;
  return fDataBunch.size();

}

int AliHLTTPCDigitReaderPacked::GetBunchSize(){
  // see header file for class documentation
  return fDataBunch.size();
}

const UInt_t* AliHLTTPCDigitReaderPacked::GetSignals()
{
  // see header file for class documentation
  return &fDataBunch[0];
}

Int_t* AliHLTTPCDigitReaderPacked::GetBufferInstance()
{
  // see header file for class documentation

  // for the moment only a singleton of the buffer is foreseen
  // could be extended but very unlikly to be worth the effort
  // because pad sorting is just a debug feature.

  // This is just a poor man's solution, no synchronization for the
  // moment
  AliHLTLogging log;
  if (fgpIssuedBufferInstance) {
    log.LoggingVarargs(kHLTLogError, "AliHLTTPCDigitReaderPacked", "GetBufferInstance" , __FILE__ , __LINE__ ,
		       "instance of sorted buffer has not been released or multiple instances requested. Only available as global singleton for DigitReaderPacked");
    return NULL;
  }

  if (!fgpFreeBufferInstance) {
    if (fNMaxRows<0 || fNMaxPads<0 || fNTimeBins<0) {
      log.LoggingVarargs(kHLTLogError, "AliHLTTPCDigitReaderPacked", "GetBufferInstance" , __FILE__ , __LINE__ ,
			 "can not determine size of buffer for sorted data");
      return NULL;
    }
    fgpFreeBufferInstance=new Int_t[ fNMaxRows*fNMaxPads*fNTimeBins ];
    log.LoggingVarargs(kHLTLogDebug, "AliHLTTPCDigitReaderPacked", "GetBufferInstance" , __FILE__ , __LINE__ , 
		       "Array Borders  ||| MAXPAD=%d ||| MAXROW=%d ||| MAXBIN=%d ||| MAXMUL=%d", 
		       fNMaxPads, fNMaxRows, fNTimeBins, fNTimeBins*fNMaxRows*fNMaxPads);
  }

  fgpIssuedBufferInstance=fgpFreeBufferInstance;
  fgpFreeBufferInstance=NULL;
  return fgpIssuedBufferInstance;
}

void AliHLTTPCDigitReaderPacked::ReleaseBufferInstance(Int_t* pInstance)
{
  // see header file for class documentation
  if (!pInstance) return;
  if (pInstance!=fgpIssuedBufferInstance) {
    AliHLTLogging log;
    log.LoggingVarargs(kHLTLogError, "AliHLTTPCDigitReaderPacked", "ReleaseBufferInstance" , __FILE__ , __LINE__ ,
		       "wrong instance %p, expecting %p", pInstance, fgpIssuedBufferInstance);
    return;
  }
  fgpFreeBufferInstance=fgpIssuedBufferInstance;
  fgpIssuedBufferInstance=NULL;
}

AliHLTTPCDigitReaderPacked::AliHLTTPCRawStream* AliHLTTPCDigitReaderPacked::GetRawStreamInstance()
{
  // see header file for class documentation
  if (fgpIssuedStreamInstance) {
    HLTError("instance of TPCRawStream has not been released or multiple instances requested. Only available as global singleton for DigitReaderPacked");
    return NULL;
  }

  if (!fgpFreeStreamInstance) {
    fgpFreeStreamInstance=new AliHLTTPCDigitReaderPacked::AliHLTTPCRawStream;
  }

  fgpIssuedStreamInstance=fgpFreeStreamInstance;
  fgpFreeStreamInstance=NULL;
  return fgpIssuedStreamInstance;
}

void AliHLTTPCDigitReaderPacked::ReleaseRawStreamInstance(AliHLTTPCDigitReaderPacked::AliHLTTPCRawStream* pInstance)
{
  // see header file for class documentation
  if (!pInstance) return;
  if (pInstance!=fgpIssuedStreamInstance) {
    HLTError("wrong instance %p, expecting %p", pInstance, fgpIssuedStreamInstance);
    return;
  }
  fgpFreeStreamInstance=fgpIssuedStreamInstance;
  fgpIssuedStreamInstance=NULL;
}

AliHLTTPCDigitReaderPacked::AliHLTTPCRawStream::AliHLTTPCRawStream()
  :
  fRawMemoryReader(new AliRawReaderMemory),
  fTPCRawStream(new AliTPCRawStream(fRawMemoryReader))
{
  // see header file for class documentation
}

AliHLTTPCDigitReaderPacked::AliHLTTPCRawStream::~AliHLTTPCRawStream()
{
  // see header file for class documentation
  if (fRawMemoryReader) delete fRawMemoryReader;
  fRawMemoryReader=NULL;
  if (fTPCRawStream) delete fTPCRawStream;
  fTPCRawStream=NULL;
}

Bool_t AliHLTTPCDigitReaderPacked::AliHLTTPCRawStream::SetMemory(Int_t ddlId, UChar_t* memory, ULong_t size )
{
  // see header file for class documentation
  if (!fRawMemoryReader) return false;
  Bool_t result=fRawMemoryReader->SetMemory(memory, size );
  fRawMemoryReader->SetEquipmentID(ddlId);
  fRawMemoryReader->RewindEvents();
  fRawMemoryReader->NextEvent();
  return result;
}

bool AliHLTTPCDigitReaderPacked::AliHLTTPCRawStream::Next()
{
  // see header file for class documentation
  if (!fTPCRawStream) return -1;
  return fTPCRawStream->Next();
}

Int_t AliHLTTPCDigitReaderPacked::AliHLTTPCRawStream::GetRow() const
{
  // see header file for class documentation
  if (!fTPCRawStream) return -1;
  return fTPCRawStream->GetRow();
}

Int_t AliHLTTPCDigitReaderPacked::AliHLTTPCRawStream::GetPad() const
{
  // see header file for class documentation
  if (!fTPCRawStream) return -1;
  return fTPCRawStream->GetPad();
}

Int_t AliHLTTPCDigitReaderPacked::AliHLTTPCRawStream::GetTime() const
{
  // see header file for class documentation
  if (!fTPCRawStream) return -1;
  return fTPCRawStream->GetTime();
}

Int_t AliHLTTPCDigitReaderPacked::AliHLTTPCRawStream::GetSignal() const
{
  // see header file for class documentation
  if (!fTPCRawStream) return -1;
  return fTPCRawStream->GetSignal();
}

Int_t AliHLTTPCDigitReaderPacked::AliHLTTPCRawStream::GetHWAddress() const
{
  // see header file for class documentation
  if (!fTPCRawStream) return -1;
  return fTPCRawStream->GetHWAddress();
}

Bool_t  AliHLTTPCDigitReaderPacked::AliHLTTPCRawStream::GetRCUTrailerData(UChar_t*& data) const
{
  // see header file for class documentation
  if (!fTPCRawStream) return -1;
  return fTPCRawStream->GetRCUTrailerData(data);
}

Int_t   AliHLTTPCDigitReaderPacked::AliHLTTPCRawStream::GetRCUTrailerSize() const
{
  // see header file for class documentation
  if (!fTPCRawStream) return -1;
  return fTPCRawStream->GetRCUTrailerSize();
}
