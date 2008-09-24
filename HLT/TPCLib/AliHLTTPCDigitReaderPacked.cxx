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
    @brief  A digit reader implementation for simulated, packed TPC 'raw' data.
*/

#if __GNUC__>= 3
using namespace std;
#endif

#include "AliHLTTPCDigitReaderPacked.h"

#include "AliTPCRawStream.h"
#include "AliRawReaderMemory.h"
#include "AliRawDataHeader.h"

//#if ENABLE_PAD_SORTING
#include "AliHLTTPCTransform.h"
//#endif // ENABLE_PAD_SORTING
#include "AliHLTStdIncludes.h"

ClassImp(AliHLTTPCDigitReaderPacked)

AliHLTTPCDigitReaderPacked::AliHLTTPCDigitReaderPacked()
  :
  fRawMemoryReader(NULL),
  fTPCRawStream(NULL),
  //#if ENABLE_PAD_SORTING
  fCurrentRow(0),
  fCurrentPad(0),
  fCurrentBin(-1),
  fRowOffset(0),
  fNRows(0),
  fData(NULL),
  //#endif // ENABLE_PAD_SORTING  
  fUnsorted(kFALSE),
  fDataBunch(),
  fNextChannelFlag(kFALSE),
  fCurrentPatch(0)
{
  fRawMemoryReader = new AliRawReaderMemory;
  
  fTPCRawStream = new AliTPCRawStream( fRawMemoryReader );

  //#if ENABLE_PAD_SORTING

  // Matthias Sep 2008: the pad sorting functionality needs a deep
  // revision of the code
  // I just stumbled over a few awkward realizations
  // - each instance allocates the buffer for sorted data, this is
  //   approx. 8 MByte each
  // - each instance to loops to extract the buffer size
  //
  // quick reaction: there is now an instance handling of the buffer
  // and the buffer size is only calculated once
  // The sorting of pads is going to be a common functionality of the
  // DigitReader base class

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

  //#endif // ENABLE_PAD_SORTING
}

Int_t AliHLTTPCDigitReaderPacked::fNMaxRows=-1;
Int_t AliHLTTPCDigitReaderPacked::fNMaxPads=-1;
Int_t AliHLTTPCDigitReaderPacked::fNTimeBins=-1;
Int_t* AliHLTTPCDigitReaderPacked::fgpFreeInstance=NULL;
Int_t* AliHLTTPCDigitReaderPacked::fgpIssuedInstance=NULL;

AliHLTTPCDigitReaderPacked::~AliHLTTPCDigitReaderPacked()
{
  if (fData)
    ReleaseBufferInstance(fData);

  if ( fRawMemoryReader )
    delete fRawMemoryReader;
  fRawMemoryReader = NULL;
  if ( fTPCRawStream )
      delete fTPCRawStream;
  fTPCRawStream = NULL;
}

Int_t AliHLTTPCDigitReaderPacked::InitBlock(void* ptr,ULong_t size, Int_t patch, Int_t slice)
{

  fRawMemoryReader->SetMemory( reinterpret_cast<UChar_t*>( ptr ), size );

  fCurrentPatch=patch;
  
  //get DDL ID in order to tell the memory reader which slice/patch to use
  Int_t DDLid= 0;
  if (patch < 2)
    DDLid = 768 + 2*slice + patch;
  else 
    DDLid = 840 + 4*slice + patch-2;

  fRawMemoryReader->SetEquipmentID(DDLid);
  //fRawMemoryReader->SetEquipmentID(1);
  fRawMemoryReader->RewindEvents();
  fRawMemoryReader->NextEvent();

  if(!fUnsorted){
  //#if ENABLE_PAD_SORTING

  fCurrentRow = 0;
  fCurrentPad = 0;
  fCurrentBin = -1;

  Int_t firstrow=AliHLTTPCTransform::GetFirstRow(patch);
  Int_t lastrow=AliHLTTPCTransform::GetLastRow(patch);
  fNRows = lastrow - firstrow + 1;

  Int_t offset=0;
  if (patch > 1) offset =  AliHLTTPCTransform::GetFirstRow( 2 );

  fRowOffset = firstrow - offset;
  firstrow -= offset;
  lastrow  -= offset;

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
  //#endif // ENABLE_PAD_SORTING
  }
  return 0;
}

int AliHLTTPCDigitReaderPacked::Reset()
{
  // see header file for class documentation
  if (fData) ReleaseBufferInstance(fData);
  fData=NULL;
  return 0;
}

Bool_t AliHLTTPCDigitReaderPacked::NextSignal(){
  Bool_t readvalue = kTRUE;

  if(!fUnsorted){//added for test
    if (!fData) return false;
    //#if ENABLE_PAD_SORTING
    while (1) {
      fCurrentBin++;
      if (fCurrentBin >= fNTimeBins){
	  fCurrentBin = 0;
	  fCurrentPad++;
     
	  if (fCurrentPad >=fNMaxPads){
	      fCurrentPad = 0;
	      fCurrentRow++;
	      
	      if (fCurrentRow >= fNMaxRows){
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
  }// added for test
  else{//added for test
    //#else // !ENABLE_PAD_SORTING
    readvalue = fTPCRawStream->Next();
  }//added for test
  //#endif // ENABLE_PAD_SORTING

  return readvalue;
}

Int_t AliHLTTPCDigitReaderPacked::GetRow(){
  /*#if ENABLE_PAD_SORTING
  return (fCurrentRow + fRowOffset);
#else // !ENABLE_PAD_SORTING
  return (Int_t) fTPCRawStream->GetRow();
#endif // ENABLE_PAD_SORTING
  */
  if(!fUnsorted){
  return (fCurrentRow + fRowOffset);
  }
  else{
    if(fCurrentPatch>1){
      return (Int_t) fTPCRawStream->GetRow()-AliHLTTPCTransform::GetFirstRow(fCurrentPatch)+AliHLTTPCTransform::GetFirstRow(2);
    }
    else{
      return (Int_t) fTPCRawStream->GetRow()-AliHLTTPCTransform::GetFirstRow(fCurrentPatch);
    }
  }
}

int AliHLTTPCDigitReaderPacked::GetPad(){
  /*#if ENABLE_PAD_SORTING
    return fCurrentPad;
    #else // !ENABLE_PAD_SORTING
    return fTPCRawStream->GetPad();
    #endif // ENABLE_PAD_SORTING
  */
  if(!fUnsorted){
    return fCurrentPad;
  }
  else{
    return fTPCRawStream->GetPad();
  }
}

AliHLTUInt32_t AliHLTTPCDigitReaderPacked::GetAltroBlockHWaddr() const
{
  return fTPCRawStream->GetHWAddress();
}

Int_t AliHLTTPCDigitReaderPacked::GetSignal(){ 
  /*
    #if ENABLE_PAD_SORTING
    return fData[ fCurrentRow*fNMaxPads*fNTimeBins+ fCurrentPad*fNTimeBins + fCurrentBin ];
    #else // !ENABLE_PAD_SORTING
    return fTPCRawStream->GetSignal();
    #endif // ENABLE_PAD_SORTING
  */
  if(!fUnsorted){
    // check for validity of fData is in NextSignal, no check at here
    return fData[ fCurrentRow*fNMaxPads*fNTimeBins+ fCurrentPad*fNTimeBins + fCurrentBin ];
  }
  else{
    return fTPCRawStream->GetSignal();
  }
}

Int_t AliHLTTPCDigitReaderPacked::GetTime(){
  /*
    #if ENABLE_PAD_SORTING
    return fCurrentBin;
    #else // !ENABLE_PAD_SORTING
    return fTPCRawStream->GetTime();
    #endif // ENABLE_PAD_SORTING
  */
  if(!fUnsorted){
    return fCurrentBin;
  }
  else{
    if((Int_t)(fTPCRawStream->GetTime()-fDataBunch.size()+1)>0 &&(Int_t)(fTPCRawStream->GetTime()-fDataBunch.size()+1)<=AliHLTTPCTransform::GetNTimeBins()){
      return fTPCRawStream->GetTime()-fDataBunch.size()+1;
    }
    else{
      HLTDebug("Timebin is out of range: %d",fTPCRawStream->GetTime()-fDataBunch.size()+1);
      return 0;
    }
  }
}

Int_t AliHLTTPCDigitReaderPacked::GetTimeOfUnsortedSignal(){
  return fTPCRawStream->GetTime();
}

bool AliHLTTPCDigitReaderPacked::NextChannel(){
  bool iResult=false;
  if(fNextChannelFlag==kFALSE){
    if(!NextSignal()){//if there are no more signals
      iResult=false;
    }
    else{
      iResult=true;
    }
  }
  else{
    iResult=true;
  }
  return iResult;
}

int AliHLTTPCDigitReaderPacked::NextBunch(){
  
  fDataBunch.clear();
  //adding the first signal (will always be the leftover from either NextChannel call or Previous bunch)
  fDataBunch.push_back(GetSignal());

  int iResult=1;
  Bool_t continueLoop=kTRUE;
  AliHLTUInt32_t prevHWAddress=GetAltroBlockHWaddr();
  Int_t prevTime=GetTimeOfUnsortedSignal();
  do{
    if(NextSignal()){
      if(GetAltroBlockHWaddr()==prevHWAddress){//check if there is a change in channel(new row and pad)
	if(prevTime==GetTimeOfUnsortedSignal()+1){//if true means that we have consecutive signals
	  prevTime=GetTimeOfUnsortedSignal();
	  fDataBunch.push_back(GetSignal());
	}
	else{//end of bunch but not of channel
	  continueLoop=kFALSE;
	}
      }
      else{
	iResult=0;//end of bunch
	continueLoop=kFALSE;
	fNextChannelFlag=kTRUE;
      }
    }
    else{
      continueLoop=kFALSE;
      fNextChannelFlag=kFALSE;
      if(fDataBunch.size()>0){//we reached end of data in total, but we still have a bunch
	iResult = 0;
      }
    }
  }while(continueLoop);

  return iResult;

}

int AliHLTTPCDigitReaderPacked::GetBunchSize(){
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
  if (fgpIssuedInstance) {
    log.LoggingVarargs(kHLTLogError, "AliHLTTPCDigitReaderPacked", "GetBufferInstance" , __FILE__ , __LINE__ ,
		       "instance of sorted buffer has not been released or multiple instances requested. Only available as global singleton for DigitReaderPacked");
    return NULL;
  }

  if (!fgpFreeInstance) {
    if (fNMaxRows<0 || fNMaxPads<0 || fNTimeBins<0) {
      log.LoggingVarargs(kHLTLogError, "AliHLTTPCDigitReaderPacked", "GetBufferInstance" , __FILE__ , __LINE__ ,
			 "can not determine size of buffer for sorted data");
      return NULL;
    }
    fgpFreeInstance=new Int_t[ fNMaxRows*fNMaxPads*fNTimeBins ];
    log.LoggingVarargs(kHLTLogDebug, "AliHLTTPCDigitReaderPacked", "GetBufferInstance" , __FILE__ , __LINE__ , 
		       "Array Borders  ||| MAXPAD=%d ||| MAXROW=%d ||| MAXBIN=%d ||| MAXMUL=%d", 
		       fNMaxPads, fNMaxRows, fNTimeBins, fNTimeBins*fNMaxRows*fNMaxPads);
  }

  fgpIssuedInstance=fgpFreeInstance;
  fgpFreeInstance=NULL;
  return fgpIssuedInstance;
}

void AliHLTTPCDigitReaderPacked::ReleaseBufferInstance(Int_t* pInstance)
{
  // see header file for class documentation
  if (!pInstance) return;
  if (pInstance!=fgpIssuedInstance) {
    AliHLTLogging log;
    log.LoggingVarargs(kHLTLogError, "AliHLTTPCDigitReaderPacked", "ReleaseBufferInstance" , __FILE__ , __LINE__ ,
		       "wrong instance %p, expecting %p", pInstance, fgpIssuedInstance);
    return;
  }
  fgpFreeInstance=fgpIssuedInstance;
  fgpIssuedInstance=NULL;
}
