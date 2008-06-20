// $Id$

//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
//*                  Timm Steinbeck <timm@kip.uni-heidelberg.de>           *
//*                  Jochen Thaeder <thaeder@kip.uni-heidelberg.de>        *
//*                  for The ALICE HLT Project.                            *
//*                                                                        *
//* Permission to use, copy, modify and distribute this software and its   *
//* documentation strictly for non-commercial purposes is hereby granted   *
//* without fee, provided that the above copyright notice appears in all   *
//* copies and that both the copyright notice and this permission notice   *
//* appear in the supporting documentation. The authors make no claims     *
//* about the suitability of this software for any purpose. It is          *
//* provided "as is" without express or implied warranty.                  *
//**************************************************************************

/** @file   AliHLTTPCDigitReaderDecoder.cxx
    @author Kenneth Aamodt, Matthias Richter
    @date   
    @brief  DigitReader implementation for the fast ALTRO Decoder
*/

#if __GNUC__>= 3
using namespace std;
#endif

#include <cassert>
#include "AliHLTTPCDigitReaderDecoder.h"
#include "AliHLTTPCMapping.h"
#include "AliAltroDecoder.h"
#include "AliAltroData.h"
#include "AliAltroBunch.h"
#include "AliHLTTPCTransform.h"

ClassImp(AliHLTTPCDigitReaderDecoder)

AliHLTTPCDigitReaderDecoder::AliHLTTPCDigitReaderDecoder()
  :
  AliHLTTPCDigitReader(),
  fAltroDecoder(NULL),
  fAltroData(),
  fAltroBunch(NULL),
  fMapping(NULL),
  // initialization due to the logic in NextSignals
  fNextCounter(-1),
  fNextSignalMethodUsed(kFALSE)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTTPCDigitReaderDecoder::~AliHLTTPCDigitReaderDecoder()
{
  // see header file for class documentation
  if(fAltroDecoder){
    delete fAltroDecoder;
  }
  if(fAltroBunch){
    delete fAltroBunch;
  }
  if(fMapping){
    delete fMapping;
  }
}

int AliHLTTPCDigitReaderDecoder::InitBlock(void* ptr,unsigned long size, Int_t patch, Int_t /*slice*/)
{
  // see header file for class documentation
  //  HLTDebug("Initializing block in decoder");
  if(!fMapping){
    fMapping = new AliHLTTPCMapping(patch);
  }
  if(!fAltroDecoder){
    fAltroDecoder = new AliAltroDecoder();
  }
  if(!fAltroBunch){
    fAltroBunch = new AliAltroBunch();
  }
  fAltroDecoder->SetMemory((UChar_t*)ptr, size);
  fAltroDecoder->Decode();
  return 0;
}

void AliHLTTPCDigitReaderDecoder::SetUnsorted(bool unsorted)
{
  // see header file for class documentation

  // The DigitReaderDecoder does not support sorted data, forward to
  // default if sorted data requested
  if (!unsorted) AliHLTTPCDigitReader::SetUnsorted(unsorted);
}

bool AliHLTTPCDigitReaderDecoder::NextChannel()
{
  // see header file for class documentation
  Bool_t result=fAltroDecoder->NextChannel(&fAltroData);
  if(result && !fMapping->IsValidHWAddress(fAltroData.GetHadd())){
    result = fAltroDecoder->NextChannel(&fAltroData);
  }
  return result;
}

int AliHLTTPCDigitReaderDecoder::NextBunch()
{
  // see header file for class documentation
  return fAltroData.NextBunch(fAltroBunch);
}

bool AliHLTTPCDigitReaderDecoder::NextSignal()
{
  // see header file for class documentation
  fNextSignalMethodUsed=kTRUE;
  do {
    if (fNextCounter>0) {
      // there is data available in the current bunch
      fNextCounter--;
      return true;
    }

    // there is no data left in the current bunch, search for the next one
    while (NextBunch()) if (GetBunchSize()>0) {
      fNextCounter=GetBunchSize()-1;
      return true;
    }

    fNextCounter=-1;
    // there is no bunch left, go to the next channel
  } while (NextChannel());
  
  return false;
}

const UInt_t* AliHLTTPCDigitReaderDecoder::GetSignals()
{
  // see header file for class documentation
  return fAltroBunch->GetData();
}

int AliHLTTPCDigitReaderDecoder::GetRow()
{
  // see header file for class documentation
  return fMapping->GetRow(fAltroData.GetHadd());
}

int AliHLTTPCDigitReaderDecoder::GetPad()
{
  // see header file for class documentation
  return fMapping->GetPad(fAltroData.GetHadd());
}

int AliHLTTPCDigitReaderDecoder::GetSignal()
{
  // see header file for class documentation
  if (fNextSignalMethodUsed) {
    const  UInt_t* pData=GetSignals();
    if (pData && fNextCounter>=0) {
      assert(fNextCounter<GetBunchSize());
      return pData[fNextCounter];
    }
  }
  return 0;
}

int AliHLTTPCDigitReaderDecoder::GetTime()
{
  // see header file for class documentation
  int iResult=0;
  if(!fNextSignalMethodUsed){// this is true if the bunch approach is used
    
    iResult= fAltroBunch->GetStartTimeBin();
  }
  else{
    assert(fNextCounter>=0);
    iResult = fAltroBunch->GetStartTimeBin()+fNextCounter;
  }
  if(iResult<0 || iResult>AliHLTTPCTransform::GetNTimeBins()){
    iResult=0;
  }
  return iResult;
}

int AliHLTTPCDigitReaderDecoder::GetBunchSize()
{
  // see header file for class documentation
  return fAltroBunch->GetBunchSize();
}

int AliHLTTPCDigitReaderDecoder::GetRowOffset() const
{
  return fMapping->GetRowOffset();
}
AliHLTUInt32_t AliHLTTPCDigitReaderDecoder::GetAltroBlockHWaddr() const
{
  // see header file for class documentation
  return (AliHLTUInt32_t)fAltroData.GetHadd();
}
AliHLTUInt32_t AliHLTTPCDigitReaderDecoder::GetAltroBlockHWaddr(Int_t row, Int_t pad) const
{
  // see header file for class documentation
  if(fMapping){
    return fMapping->GetHwAddress(row,pad);
  }
  else{
    return 0;
  }
}
