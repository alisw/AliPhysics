// $Id$

/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 *                                                                        *
 * Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
 *                  Timm Steinbeck <timm@kip.uni-heidelberg.de>           *
 *                  Jochen Thaeder <thaeder@kip.uni-heidelberg.de>        *
 *                  for The ALICE HLT Project.                            *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/** @file   AliHLTTPCDigitReaderUnpacked.cxx
    @author Timm Steinbeck, Jochen Thaeder, Matthias Richter
    @date   
    @brief  A digit reader implementation for unpacked TPC data.
*/

#if __GNUC__== 3
using namespace std;
#endif

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
  AliHLTTPCUnpackedRawData *tmpptr=NULL;
  fPtr = ptr;
  fSize = size;

  tmpptr = reinterpret_cast<AliHLTTPCUnpackedRawData*>(fPtr);
  fDigitRowData = (AliHLTTPCDigitRowData*) tmpptr->fDigits;
  fActRowData = (AliHLTTPCDigitRowData*) fDigitRowData;

  fBin = -1;

  int dummy=0;
  AliHLTTPCTransform::Slice2Sector(slice, AliHLTTPCTransform::GetFirstRow(patch), dummy, fFirstRow);
  AliHLTTPCTransform::Slice2Sector(slice, AliHLTTPCTransform::GetLastRow(patch), dummy, fLastRow);

  fRow = fFirstRow; 

  if ((Int_t)fActRowData->fRow != fRow){
      HLTWarning("Row number should match! fActRowData->fRow=%d fRow=%d", fActRowData->fRow, fRow);
  }
  return 0;
}

bool AliHLTTPCDigitReaderUnpacked::Next(){
  // see header file for class documentation
  bool rreadvalue = true;

  fBin++;

  if ( fBin >= (Int_t)fActRowData->fNDigit ){

    fRow++;

    if ((fRow >= fFirstRow) && (fRow <= fLastRow)){

      //new row 
      Byte_t *tmp = (Byte_t*) fActRowData;
      Int_t size = sizeof(AliHLTTPCDigitRowData) + fActRowData->fNDigit*sizeof(AliHLTTPCDigitData);
      tmp += size;
      fActRowData = (AliHLTTPCDigitRowData*) tmp;
     
      if (((Byte_t*)fPtr) + fSize <= tmp){
	rreadvalue = false;
	return rreadvalue;
      }

      fBin = 0;
    }
    else {
      rreadvalue = false;
      return rreadvalue;
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
  int rpad;
  rpad = (int)fData[fBin].fPad;
  return rpad   ;
}

int AliHLTTPCDigitReaderUnpacked::GetSignal(){ 
  // see header file for class documentation
  int rsignal;
  rsignal = (int)fData[fBin].fCharge;
  return rsignal;
}

int AliHLTTPCDigitReaderUnpacked::GetTime(){
  // see header file for class documentation
  int rtime;
  rtime = (int)fData[fBin].fTime;
  return rtime;
}
