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

/** @file   AliHLTTPCDigitReaderUnpacked.cxx
    @author Timm Steinbeck, Jochen Thaeder, Matthias Richter
    @date   
    @brief  A digit reader implementation for unpacked TPC data.
*/

#if __GNUC__== 3
using namespace std;
#endif

#include <cassert>
#include "AliHLTTPCDigitReaderUnpacked.h"
#include "AliHLTTPCDigitData.h"
#include "AliHLTTPCTransform.h"
#include "AliHLTStdIncludes.h"
#include "AliHLTTPCMapping.h"

ClassImp(AliHLTTPCDigitReaderUnpacked)

AliHLTTPCDigitReaderUnpacked::AliHLTTPCDigitReaderUnpacked()
  :
  fDigitRowData(NULL),
  fActRowData(NULL),
  fData(NULL),
  fPtr(NULL),
  fSize(0),
  fBin(0),
  fRow(0),
  fFirstRow(0),
  fLastRow(0),
  fUnsorted(kFALSE),
  fDataBunch(),
  fTrackIDs(),
  fTrackIDCounts(),
  fEndOfDataReached(kFALSE),
  fEndOfChannelReached(kFALSE),
  fPrevTime(0),
  fEndTimeBinOfBunch(0),
  fPrevSignal(0),
  fPrevPad(-1),
  fPrevRow(-1),
  fNextChannelIsAlreadyConfirmed(kFALSE),
  fMapping(NULL),
  fDigitsVector(),
  fBinRowPositionSorted(),
  fPatch(0)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTTPCDigitReaderUnpacked::~AliHLTTPCDigitReaderUnpacked(){
  // see header file for class documentation
  if(fMapping){
    delete fMapping;
  }
  fMapping=NULL;
}

int AliHLTTPCDigitReaderUnpacked::InitBlock(void* ptr,unsigned long size, Int_t patch, Int_t slice){
  // see header file for class documentation
  int iResult=0;
  AliHLTTPCUnpackedRawData *tmpptr=NULL;
  fPtr = ptr;
  fSize = size;

  fPatch=patch;
  fEndOfDataReached=kFALSE;
  fEndOfChannelReached=kFALSE;
  fPrevTime=0;
  fEndTimeBinOfBunch=0;
  fPrevSignal=0;
  fPrevPad=-1;
  fPrevRow=-1;
  fNextChannelIsAlreadyConfirmed=kFALSE;
  fDigitsVector.clear();

  tmpptr = reinterpret_cast<AliHLTTPCUnpackedRawData*>(fPtr);
  fDigitRowData = (AliHLTTPCDigitRowData*) tmpptr->fDigits;
  fActRowData = fDigitRowData;

  if(!fMapping){
    fMapping = new AliHLTTPCMapping(patch);
  }

  while (fActRowData && ((iResult=GetNextRowData(fActRowData))>=0)) {/* empty body */};

  if (iResult>=0) {
  fActRowData = fDigitRowData;
  fBin = -1;

  int dummy=0;
  AliHLTTPCTransform::Slice2Sector(slice, AliHLTTPCTransform::GetFirstRow(patch), dummy, fFirstRow);
  AliHLTTPCTransform::Slice2Sector(slice, AliHLTTPCTransform::GetLastRow(patch), dummy, fLastRow);

  fRow = fFirstRow; 

  if ((Int_t)fActRowData->fRow != fRow){
      HLTWarning("Row number should match! fActRowData->fRow=%d fRow=%d", fActRowData->fRow, fRow);
  }
  } else {
    fActRowData=NULL;
  }
  return iResult;
}

int AliHLTTPCDigitReaderUnpacked::GetNextRowData(AliHLTTPCDigitRowData*& pRow) const
{
  // get new row data from the current row data
  int iResult=0;
  AliHLTTPCDigitRowData* pCurrent=pRow;
  assert(pCurrent);
  pRow=NULL;
  Byte_t *tmp = (Byte_t*) pCurrent;
  Int_t size = sizeof(AliHLTTPCDigitRowData) + pCurrent->fNDigit*sizeof(AliHLTTPCDigitData);
  tmp += size;
  pRow = reinterpret_cast<AliHLTTPCDigitRowData*>(tmp);
  
  // check if the new pointer is within the range
  if (((Byte_t*)fPtr) + fSize <= tmp){
    if (((Byte_t*)fPtr) + fSize < tmp) {
      // if the increment does not match exactly there is a format error
      HLTError("input format not recognized: buffer %p %d, current row data %p, %d digits", fPtr, fSize, pCurrent, pCurrent->fNDigit);
      iResult=-EBADF;
    }
    pRow=NULL;
  } else {
    // check if the current row structure has the right format
    size = sizeof(AliHLTTPCDigitRowData) + pRow->fNDigit*sizeof(AliHLTTPCDigitData);
    tmp += size;
    if (((Byte_t*)fPtr) + fSize < tmp){
      HLTError("Current row data not recognized %p (buffer %p %d) %d digits", pRow, fPtr, fSize, pRow->fNDigit);
      pRow=NULL;
      iResult=-EBADF;
    }
  }
    
  return iResult;
}
void AliHLTTPCDigitReaderUnpacked::SortBunchBinVector(){
  fBinRowPositionSorted.clear();

  fBinRowPositionSorted.push_back(0);
  Int_t nAdded=0;
  Int_t beginningOfPadIndex=0;
  Int_t totalAdded=0;
  for(Int_t i=1;i<(Int_t)fActRowData->fNDigit;i++){
     if(fData[i-1].fPad == fData[i].fPad){// means that these sinals belong to the same pad
      if(fData[i-1].fTime+1 == fData[i].fTime){ //means that the signal belong to the same bunch
	nAdded++;
	totalAdded++;
	fBinRowPositionSorted.insert(fBinRowPositionSorted.begin()+beginningOfPadIndex+nAdded,i);
      }
      else{//we have a new bunch on this pad, put it in fornt of the previous bunch
	totalAdded++;
	nAdded=0;
	fBinRowPositionSorted.insert(fBinRowPositionSorted.begin()+beginningOfPadIndex,i);
      }
    }
    else{
      totalAdded++;
      beginningOfPadIndex=totalAdded;
      fBinRowPositionSorted.push_back(i);
      nAdded=0;
    }
  }
}

bool AliHLTTPCDigitReaderUnpacked::NextSignal(){
  // see header file for class documentation
  if (fActRowData==NULL) return false;

  bool rreadvalue = true;

  fBin++;

  while ( fBin >= (Int_t)fActRowData->fNDigit ){
    fRow++;
    if ((fRow >= fFirstRow) && (fRow <= fLastRow)){
      
      //new row 
      if (GetNextRowData(fActRowData)<0) {
	rreadvalue = false;
	return rreadvalue;
      }
      fBin = 0;
    }
    else {
      rreadvalue = false;
      return rreadvalue;
    }
    if(!fActRowData){
      return false;
    }
    if ((Int_t)fActRowData->fRow != fRow){
      HLTWarning("Row number should match! fActRowData->fRow=%d fRow=%d", fActRowData->fRow, fRow);
    }
  }
  
  fData = fActRowData->fDigitData;

  if(fBin==0){
    SortBunchBinVector(); 
  }

  if(fPrevPad == -1){
    fPrevPad = GetSortedPad();
  }
 
  return rreadvalue;
}

int AliHLTTPCDigitReaderUnpacked::GetRow(){
  // see header file for class documentation
  int rrow;

  if(fUnsorted == kFALSE){
    rrow = fRow;
  }
  else{
   rrow = fRow-AliHLTTPCTransform::GetFirstRow(fPatch);
   if(fPatch>1){
    rrow += AliHLTTPCTransform::GetFirstRow(2);
   }
  }

  return rrow;
}

int AliHLTTPCDigitReaderUnpacked::GetPad(){
  // see header file for class documentation
  int rpad;
  if(fUnsorted == kFALSE){
    rpad = GetSortedPad();
  }
  else{
    rpad = fPrevPad;
  }
  return rpad   ;
}

int AliHLTTPCDigitReaderUnpacked::GetSignal(){ 
  // see header file for class documentation
  int rsignal;
  if(fUnsorted == kFALSE){
    rsignal = GetSortedSignal();
  }
  else{
    rsignal = fPrevSignal;
  }
  return rsignal;
}

int AliHLTTPCDigitReaderUnpacked::GetTime(){
  // see header file for class documentation
  int rtime;
  if(fUnsorted==kFALSE){
    rtime = GetSortedTime();
  }
  else{
    rtime = fPrevTime+1-fDataBunch.size();
  }
  return rtime;
}

AliHLTUInt32_t AliHLTTPCDigitReaderUnpacked::GetAltroBlockHWaddr() const
{
  // see header file for class documentation
  return (AliHLTUInt32_t)(fMapping->GetHwAddress((UInt_t)GetSortedRow(),(UInt_t)GetSortedPad()));//fTPCRawStream->GetHWAddress();
}

Int_t AliHLTTPCDigitReaderUnpacked::GetSortedTime(){
  // see header file for class documentation
  assert(fData);
  if (!fData) return -1;
  int rtime;
  rtime = (int)fData[fBinRowPositionSorted.at(fBin)].fTime;
  if(fDataBunch.size()>1){
    fEndTimeBinOfBunch=rtime;
  }
  return rtime;
}

Int_t AliHLTTPCDigitReaderUnpacked::GetSortedSignal(){
  // see header file for class documentation
  assert(fData);
  if (!fData) return -1;
  int rsignal;
  rsignal = (int)fData[fBinRowPositionSorted.at(fBin)].fCharge;
  return rsignal;

}

Int_t AliHLTTPCDigitReaderUnpacked::GetSortedPad() const{
  // see header file for class documentation
  assert(fData);
  if (!fData) return -1;
  int rpad;
  rpad = (int)fData[fBinRowPositionSorted.at(fBin)].fPad;
  return rpad   ;
}

int AliHLTTPCDigitReaderUnpacked::GetSortedRow() const {
  // see header file for class documentation
  int rrow;
  rrow = fRow-AliHLTTPCTransform::GetFirstRow(fPatch);
  if(fPatch>1){
    rrow += AliHLTTPCTransform::GetFirstRow(2);
  }
  return rrow;
}

bool AliHLTTPCDigitReaderUnpacked::NextChannel()
{
  // see header file for class documentation

  // If the next channel is already confirmed by the next bunch function
  // or there are more signals (this will only be for the first signal)
  
  if(fEndOfDataReached == kTRUE){
    return false;
  }
  if(fNextChannelIsAlreadyConfirmed){
    fNextChannelIsAlreadyConfirmed = kFALSE;
    fPrevTime=GetSortedTime();
    fPrevSignal=GetSortedSignal();
    fPrevPad = GetSortedPad();
    fPrevRow = GetSortedRow();

    if(GetAltroBlockHWaddr() == (UInt_t)-1){
      fPrevPad = -1;
      return NextChannel();
    }

    return true;
  }
  else if(NextSignal()) { // there is data

    if(GetAltroBlockHWaddr() == (UInt_t)-1){
      fPrevPad = -1;
      return NextChannel();
    }
    return true;
  }
  return false;
}

int AliHLTTPCDigitReaderUnpacked::NextBunch()
{  
  // see header file for class documentation

  if(fEndOfDataReached == kTRUE || fEndOfChannelReached == kTRUE){
    // sets fEndOfChannelReached back to false, for the next channel
    // and returns 0 to tell stop the NextBunch calls for this channel.
    fEndOfChannelReached=kFALSE;
    return 0;
  }


  fDataBunch.clear();
  fDigitsVector.clear();

  //adding the first signal (will always be the leftover from either NextChannel call or previous bunch) 
  fPrevTime=GetSortedTime();
  fPrevSignal=GetSortedSignal();
  fPrevPad = GetSortedPad();
  fPrevRow = GetSortedRow();
  fDataBunch.push_back(GetSortedSignal());
  fDigitsVector.push_back(fData[fBin]);

  do{
    if(NextSignal()){
	if((fPrevPad == GetSortedPad()) && (fPrevRow == GetSortedRow())){//check if there is a change in channel(new pad or row)
	  if(fPrevTime+1 == GetSortedTime()){//if true means that we have consecutive signals
	    fPrevTime = GetSortedTime();
	    //fDataBunch.insert(fDataBunch.begin(), GetSortedSignal());// add the signal to the beginning of the buffer	    
	    fDataBunch.push_back(GetSortedSignal());// add the signal to the beginning of the buffer
	    fDigitsVector.push_back(fData[fBin]);
	  }
	  else{//end of bunch but not of channel
	    break;
	  }
	}
	else{
	  // end of channel, last bunch will be completed
	  fEndOfChannelReached = kTRUE;
	  // the next channel is already confirmed since the next channel returned true
	  fNextChannelIsAlreadyConfirmed = kTRUE;
	  break;
	}
    }
    else{
      // end of data, but there is one bunch to be completed.
      fEndOfDataReached = kTRUE;
      break;
    }
  }while(1);

  return fDataBunch.size();
}

int AliHLTTPCDigitReaderUnpacked::GetBunchSize(){
  // see header file for class documentation
  return fDataBunch.size();
}

const UInt_t* AliHLTTPCDigitReaderUnpacked::GetSignals()
{
  // see header file for class documentation
  return &fDataBunch[0];
}

const AliHLTTPCDigitData* AliHLTTPCDigitReaderUnpacked::GetBunchDigits()
{
  // see header file for class documentation
  return &fDigitsVector[0];
}
int AliHLTTPCDigitReaderUnpacked::GetRowOffset() const
{
  // see header file for class documentation
  return AliHLTTPCTransform::GetFirstRow(fPatch);
}
