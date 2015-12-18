// $Id$

//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
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

/// @file   AliHLTDataInflaterHuffman.cxx
/// @author Matthias Richter
/// @date   2011-09-01
/// @brief  Data inflater implementation for huffman encoded data
/// @note   

#include "AliHLTDataInflaterHuffman.h"
#include "AliHLTHuffman.h"
#include "TList.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTDataInflaterHuffman)

AliHLTDataInflaterHuffman::AliHLTDataInflaterHuffman()
  : AliHLTDataInflater()
  , fHuffmanCoders()
  , fHuffmanCoderList(NULL)
  , fCurrentParameter(-1)
  , fLegacyMode(-1)
  , fInput(0)
  , fInputLength(0)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTDataInflaterHuffman::~AliHLTDataInflaterHuffman()
{
  // destructor
  if (fHuffmanCoderList) delete fHuffmanCoderList;
  fHuffmanCoderList=NULL;
}

int AliHLTDataInflaterHuffman::AddParameterDefinition(const char* name, unsigned bitLength)
{
  /// search a parameter definition in the decoder configuration, and set the index
  /// array, return reference id
  /// TODO: this code is a copy of AliHLTDataDeflaterHuffman::AddParameterDefinition
  /// make a common base class
  if (!name) return -EINVAL;
  if (!fHuffmanCoderList) return -ENODEV;
  TObject* pObj=fHuffmanCoderList->FindObject(name);
  if (!pObj) {
    HLTError("can not find decoder of id '%s'", name);
    return -ENOENT;
  }
  AliHLTHuffman* pHuffman=dynamic_cast<AliHLTHuffman*>(pObj);
  if (!pHuffman) {
    HLTError("object %s has wrong type, expected AliHLTHuffman", name);
    return -EBADF;
  }
  if (pHuffman->GetMaxBits()!=bitLength) {
    HLTError("mismatch in bitlengt: can not use decoder %s of length %d for encoding of %d bits", pHuffman->GetName(), pHuffman->GetMaxBits(), bitLength);
    return -EPERM;
  }

  fHuffmanCoders.push_back(pHuffman);
  return fHuffmanCoders.size()-1;
}

int AliHLTDataInflaterHuffman::InitDecoders(TList* decoderlist)
{
  /// init list of decoders
  /// expects to be an external pointer, valid throughout the livetime of
  /// the instance
  /// TODO: this code is a copy of AliHLTDataDeflaterHuffman::InitDecoders
  /// make a common base class
  if (!decoderlist) return -EINVAL;
  if (!fHuffmanCoderList) {
    fHuffmanCoderList=new TList;
  } else {
    if (fHuffmanCoderList->GetEntries()>0 && fHuffmanCoderList->IsOwner()) {
      HLTWarning("list of decoders owns already %d object(s), but disabling ownership now because of new external pointers");
    }
  }
  if (!fHuffmanCoderList) return -ENOMEM;
  fHuffmanCoderList->SetOwner(kFALSE);
  TIter next(decoderlist);
  TObject* pObj=NULL;
  while ((pObj=next())!=NULL) {
    AliHLTHuffman* coder=NULL;
    if ((coder=dynamic_cast<AliHLTHuffman*>(pObj))==NULL) continue;
    if (fHuffmanCoderList->FindObject(pObj->GetName())) {
      HLTError("duplicate entry of name '%s'", pObj->GetName());
      return -EEXIST;
    }
    fHuffmanCoderList->Add(pObj);
    coder->InitMaxCodeLength();
  }

  return fHuffmanCoderList->GetEntries();
}

bool AliHLTDataInflaterHuffman::NextValue(AliHLTUInt64_t& value, AliHLTUInt32_t& length)
{
  /// overloaded from AliHLTDataInflater
  /// functions reads the sequence of parameters as defined by the decoder
  /// list, than it starts at the first parameter again
  value=0;
  length=0;
  if (fLegacyMode!=0) {
  if ((++fCurrentParameter)>=(int)fHuffmanCoders.size()) fCurrentParameter=0;
  fLegacyMode=1;
  }
  if (fHuffmanCoders.size()==0 || fCurrentParameter<0) return false;
  if (fInputLength<fHuffmanCoders[fCurrentParameter]->GetMaxCodeLength() ||
      fHuffmanCoders[fCurrentParameter]->GetMaxCodeLength()==0)
  {
    AliHLTUInt64_t input=0;
    AliHLTUInt32_t inputLength=64-fInputLength;
    if (GetRemainingBitDataSizeBytes()<=sizeof(AliHLTUInt64_t)) {
      inputLength=8*GetRemainingBitDataSizeBytes();
      inputLength-=(7-GetCurrentBitInputPosition());
    }
    if (64-fInputLength<inputLength) inputLength=64-fInputLength;
    if (!InputBits(input, inputLength)) return false;
    input<<=(64-inputLength);
    input>>=fInputLength;
    fInput|=input;
    fInputLength+=inputLength;
  }
  AliHLTUInt32_t codeLength=0;
  if (!fHuffmanCoders[fCurrentParameter]->DecodeMSB(fInput, value, length, codeLength)) return false;
  if (fInputLength<codeLength) {
    HLTError("huffman decoder '%s' pretends to have %d bit(s) decoded, but only %d available",
	     fHuffmanCoders[fCurrentParameter]->GetName(), codeLength, fInputLength);
    return false;
  }
  fInput<<=codeLength;
  fInputLength-=codeLength;

  return true;
}

bool AliHLTDataInflaterHuffman::InputBit( AliHLTUInt8_t & value )
{
  /// special overload of InputBit method to consider the
  /// internal register
  if (fInputLength > 0) {
    const int shiftval=sizeof(fInput)*8 - 1;
    value = (fInput>>shiftval) & 0x1;
    fInput<<=1;
    fInputLength-=1;
    return true;
  }

  return AliHLTDataInflater::InputBit(value);
}

void AliHLTDataInflaterHuffman::Print(Option_t* option) const
{
  /// Print info
  for (vector<AliHLTHuffman*>::const_iterator coder=fHuffmanCoders.begin();
       coder!=fHuffmanCoders.end(); coder++) {
    if (!*coder) continue;
    (*coder)->Print(option);
  }
}

void AliHLTDataInflaterHuffman::Clear(Option_t * option)
{
  /// clear the object
  fCurrentParameter=-1;
  fInput=0;
  fInputLength=0;

  if (strcmp(option, "all")==0) {
    fHuffmanCoders.clear();
    if (fHuffmanCoderList) delete fHuffmanCoderList;
    fHuffmanCoderList=NULL;
  }

  AliHLTDataInflater::Clear(option);
}
