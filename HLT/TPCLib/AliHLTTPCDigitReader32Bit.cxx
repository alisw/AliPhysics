
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

/** @file   AliHLTTPCDigitReader32Bit.cxx
    @author Kenneth Aamodt
    @date   
    @brief  DigitReader implementation for the 32 bit offline decoder
*/

#if __GNUC__>= 3
using namespace std;
#endif

#include <cassert>
#include "AliHLTTPCDigitReader32Bit.h"
#include "AliHLTTPCMapping.h"
#include "AliRawReader.h"
#include "AliRawReaderMemory.h"
#ifndef HAVE_NOT_ALTRORAWSTREAMV3
#include "AliAltroRawStreamV3.h"
#endif //HAVE_NOT_ALTRORAWSTREAMV3
#include "AliHLTTPCTransform.h"

ClassImp(AliHLTTPCDigitReader32Bit)

AliHLTTPCDigitReader32Bit::AliHLTTPCDigitReader32Bit()
  :
  AliHLTTPCDigitReader(),
  fRawReader(NULL),
  fRawReaderMemory(NULL),
  fAltroRawStreamV3(NULL),
  fMapping(NULL)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
  
  // initlialized here to get more accurate comparison with the 
  // digitReaderDecoder when using SimpleComponentWrapper
  if(fRawReaderMemory ==NULL){
    fRawReaderMemory = new AliRawReaderMemory();
  }
}

AliHLTTPCDigitReader32Bit::~AliHLTTPCDigitReader32Bit()
{
  // see header file for class documentation
  if (fRawReaderMemory){
    delete fRawReaderMemory;
    fRawReaderMemory=NULL;
  }

#ifndef HAVE_NOT_ALTRORAWSTREAMV3
  if (fAltroRawStreamV3){
    delete fAltroRawStreamV3;
    fAltroRawStreamV3 = NULL;
  }
#endif //HAVE_NOT_ALTRORAWSTREAMV3

  if(fMapping){
    delete fMapping;
  }
}

int AliHLTTPCDigitReader32Bit::InitBlock(void* ptr,unsigned long size, Int_t patch, Int_t slice)
{
  // see header file for class documentation
  
  Int_t ddlno=768;
  if (patch>1) ddlno+=72+4*slice+(patch-2);
  else ddlno+=2*slice+patch;

  if(fRawReaderMemory == NULL){
    fRawReaderMemory = new AliRawReaderMemory();
  }
  if(!fRawReaderMemory){
    return -ENODEV;
  }
  fRawReaderMemory->SetMemory(reinterpret_cast<UChar_t*>(ptr), ULong_t(size));
  fRawReaderMemory->SetEquipmentID(ddlno);
  fRawReaderMemory->Reset();
  fRawReaderMemory->NextEvent();

#ifndef HAVE_NOT_ALTRORAWSTREAMV3
  if(fAltroRawStreamV3 != NULL){
    delete fAltroRawStreamV3;
    fAltroRawStreamV3=NULL;
  }
  fAltroRawStreamV3= new AliAltroRawStreamV3(fRawReaderMemory);

  if (!fAltroRawStreamV3){
    return -ENODEV;
  }
  fAltroRawStreamV3->NextDDL();
#else
  HLTError("AltroRawStreamV3 is not available in this AliRoot version");
#endif //HAVE_NOT_ALTRORAWSTREAMV3

  if(!fMapping){
    fMapping = new AliHLTTPCMapping(patch);
    if(!fMapping){
      return -ENODEV;
    }
  }
  return 0;
}

int AliHLTTPCDigitReader32Bit::Reset()
{
  // see header file for class documentation
  fRawReaderMemory->ClearBuffers();
  return 0;
}

void AliHLTTPCDigitReader32Bit::SetUnsorted(bool unsorted)
{
  // see header file for class documentation

  // The DigitReaderDecoder does not support sorted data, forward to
  // default if sorted data requested
  if (!unsorted) AliHLTTPCDigitReader::SetUnsorted(unsorted);
}

bool AliHLTTPCDigitReader32Bit::NextChannel()
{
  // see header file for class documentation
#ifndef HAVE_NOT_ALTRORAWSTREAMV3
  return fAltroRawStreamV3->NextChannel(); 
#else
  return false;
#endif //HAVE_NOT_ALTRORAWSTREAMV3

}

int AliHLTTPCDigitReader32Bit::NextBunch()
{
  // see header file for class documentation
#ifndef HAVE_NOT_ALTRORAWSTREAMV3
  return fAltroRawStreamV3->NextBunch();
#else
  return false;
#endif //HAVE_NOT_ALTRORAWSTREAMV3
}

bool AliHLTTPCDigitReader32Bit::NextSignal()
{
  // see header file for class documentation
  return false;
}

const UInt_t* AliHLTTPCDigitReader32Bit::GetSignals()
{
  // see header file for class documentation
  HLTError("AliHLTTPCDigitReader32Bit does not support the UInt_t* format, use GetSignalsShort instead");
  return 0;
}

const UShort_t* AliHLTTPCDigitReader32Bit::GetSignalsShort()
{
  // see header file for class documentation
#ifndef HAVE_NOT_ALTRORAWSTREAMV3
  return fAltroRawStreamV3->GetSignals();
#else
  return false;
#endif //HAVE_NOT_ALTRORAWSTREAMV3
}

int AliHLTTPCDigitReader32Bit::GetRow()
{
  // see header file for class documentation
#ifndef HAVE_NOT_ALTRORAWSTREAMV3
  return fMapping->GetRow(fAltroRawStreamV3->GetHWAddress());
#else
  return -1;
#endif //HAVE_NOT_ALTRORAWSTREAMV3
}

int AliHLTTPCDigitReader32Bit::GetPad()
{
  // see header file for class documentation
#ifndef HAVE_NOT_ALTRORAWSTREAMV3
  return fMapping->GetPad(fAltroRawStreamV3->GetHWAddress());
#else
  return -1;
#endif //HAVE_NOT_ALTRORAWSTREAMV3
}

int AliHLTTPCDigitReader32Bit::GetSignal()
{
  // see header file for class documentation
  return 0;
}

int AliHLTTPCDigitReader32Bit::GetTime()
{
  // see header file for class documentation
  int iResult=-1;
#ifndef HAVE_NOT_ALTRORAWSTREAMV3
  iResult=fAltroRawStreamV3->GetStartTimeBin()-fAltroRawStreamV3->GetBunchLength()+1;
#endif //HAVE_NOT_ALTRORAWSTREAMV3
  return iResult;
}

int AliHLTTPCDigitReader32Bit::GetBunchSize()
{
  // see header file for class documentation
#ifndef HAVE_NOT_ALTRORAWSTREAMV3
  return fAltroRawStreamV3->GetBunchLength();
#else
  return -1;
#endif //HAVE_NOT_ALTRORAWSTREAMV3
}

int AliHLTTPCDigitReader32Bit::GetRowOffset() const
{
  return fMapping->GetRowOffset();
}

AliHLTUInt32_t AliHLTTPCDigitReader32Bit::GetAltroBlockHWaddr() const
{
  // see header file for class documentation
#ifndef HAVE_NOT_ALTRORAWSTREAMV3
  return (AliHLTUInt32_t)fAltroRawStreamV3->GetHWAddress();
#else
  return 0;
#endif //HAVE_NOT_ALTRORAWSTREAMV3
}

AliHLTUInt32_t AliHLTTPCDigitReader32Bit::GetAltroBlockHWaddr(Int_t row, Int_t pad) const
{
  // see header file for class documentation
  if(fMapping){
    return fMapping->GetHwAddress(row,pad);
  }
  else{
    return 0;
  }
}

int AliHLTTPCDigitReader32Bit::GetRCUTrailerSize()
{
  // see header file for class documentation
  return 0;
}

bool AliHLTTPCDigitReader32Bit::GetRCUTrailerData(UChar_t*& /*trData*/)
{
  // see header file for class documentation
  return false;
}
