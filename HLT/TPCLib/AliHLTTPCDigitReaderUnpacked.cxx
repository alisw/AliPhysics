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
  fLastRow(0)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTTPCDigitReaderUnpacked::~AliHLTTPCDigitReaderUnpacked(){
  // see header file for class documentation
}

int AliHLTTPCDigitReaderUnpacked::InitBlock(void* ptr,unsigned long size, Int_t patch, Int_t slice){
  // see header file for class documentation
  int iResult=0;
  AliHLTTPCUnpackedRawData *tmpptr=NULL;
  fPtr = ptr;
  fSize = size;

  tmpptr = reinterpret_cast<AliHLTTPCUnpackedRawData*>(fPtr);
  fDigitRowData = (AliHLTTPCDigitRowData*) tmpptr->fDigits;
  fActRowData = fDigitRowData;

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

bool AliHLTTPCDigitReaderUnpacked::NextSignal(){
  // see header file for class documentation
  if (fActRowData==NULL) return false;

  bool rreadvalue = true;

  fBin++;

  if ( fBin >= (Int_t)fActRowData->fNDigit ){
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
 
  return rreadvalue;
}

int AliHLTTPCDigitReaderUnpacked::GetRow(){
  // see header file for class documentation
  int rrow;
  rrow = fRow;
  return rrow;
}

int AliHLTTPCDigitReaderUnpacked::GetPad(){
  // see header file for class documentation
  assert(fData);
  if (!fData) return -1;
  int rpad;
  rpad = (int)fData[fBin].fPad;
  return rpad   ;
}

int AliHLTTPCDigitReaderUnpacked::GetSignal(){ 
  // see header file for class documentation
  assert(fData);
  if (!fData) return -1;
  int rsignal;
  rsignal = (int)fData[fBin].fCharge;
  return rsignal;
}

int AliHLTTPCDigitReaderUnpacked::GetTime(){
  // see header file for class documentation
  assert(fData);
  if (!fData) return -1;
  int rtime;
  rtime = (int)fData[fBin].fTime;
  return rtime;
}
