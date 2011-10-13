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
    if (dynamic_cast<AliHLTHuffman*>(pObj)==NULL) continue;
    if (fHuffmanCoderList->FindObject(pObj->GetName())) {
      HLTError("duplicate entry of name '%s'", pObj->GetName());
      return -EEXIST;
    }
    fHuffmanCoderList->Add(pObj);
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
  AliHLTUInt64_t input=0;
  AliHLTUInt32_t inputLength=64;
  if (GetRemainingBitDataSizeBytes()<=sizeof(AliHLTUInt64_t)) {
    inputLength=8*GetRemainingBitDataSizeBytes();
    inputLength-=(7-GetCurrentBitInputPosition());
  }
  if (!InputBits(input, inputLength)) return false;
  if (fLegacyMode!=0) {
  if ((++fCurrentParameter)>=(int)fHuffmanCoders.size()) fCurrentParameter=0;
  fLegacyMode=1;
  }
  if (fHuffmanCoders.size()==0 || fCurrentParameter<0) return false;
  // the huffman code is decoded from bit 0 corresponding to the top node and then to
  // the left. The bitstream stores the reversed huffman code from MSB to LSB to ensure
  // a continous bit stream, that's why then input word needs to be reversed before
  // decoding.
  // TODO: introducing DecodeReverse into AliHLTHuffman can speed up the reading
  std::bitset<64> bits;
  for (unsigned bit=0; bit<inputLength; bit++) {bits[inputLength-1-bit]=input&0x1;input>>=1;}
  AliHLTUInt32_t codeLength=0;
  if (!fHuffmanCoders[fCurrentParameter]->Decode(bits, value, length, codeLength)) return false;
  if (inputLength<codeLength) {
    HLTError("huffman decoder '%s' pretends to have %d bit(s) decoded, but only %d available",
	     fHuffmanCoders[fCurrentParameter]->GetName(), codeLength, inputLength);
    return false;
  }
  inputLength-=codeLength;
  RewindBitPosition(inputLength);

  return true;
}
