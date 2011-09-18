// $Id$

//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Thorsten Kollegger <kollegge@ikf.uni-frankfurt.de>    *
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

/// @file   AliHLTDataDeflaterHuffman.cxx
/// @author Thorsten Kollegger, Matthias Richter
/// @date   2011-08-10
/// @brief  Deflater implementation using huffman code

#include "AliHLTDataDeflaterHuffman.h"
#include "AliHLTHuffman.h"
#include "TObjArray.h"
#include "TList.h"
#include "TString.h"
#include "TFile.h"
#include <memory>
#include <algorithm>
#include <iostream>

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTDataDeflaterHuffman)

AliHLTDataDeflaterHuffman::AliHLTDataDeflaterHuffman(bool bTrainingMode)
  : AliHLTDataDeflater()
  , fHuffmanCoders()
  , fHuffmanCoderList(NULL)
  , fTrainingMode(bTrainingMode)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
  if (bTrainingMode) {
    HLTInfo("using DataDeflaterHuffman in training mode");
  }
}

AliHLTDataDeflaterHuffman::~AliHLTDataDeflaterHuffman()
{
  // destructor
  if (fHuffmanCoderList) delete fHuffmanCoderList;
  fHuffmanCoderList=NULL;

  Clear();
}

int AliHLTDataDeflaterHuffman::AddParameterDefinition(const char* name, unsigned bitLength)
{
  /// search a parameter definition in the decoder configuration, and set the index
  /// array, return reference id
  if (IsTrainingMode())
    return AddTrainingParameter(name, bitLength);

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

int AliHLTDataDeflaterHuffman::InitDecoders(TList* decoderlist)
{
  /// init list of decoders
  /// expects to be an external pointer, valid throughout the livetime of
  /// the instance
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

bool AliHLTDataDeflaterHuffman::OutputParameterBits( int memberId, AliHLTUInt64_t const & value )
{
  // write huffman encoded bit pattern of a member to the current byte and position
  if (IsTrainingMode())
    return AddTrainingValue(memberId, value);

  if (memberId>=(int)fHuffmanCoders.size()) {
    return false;
  }

  AliHLTUInt64_t length = 0;
  const std::bitset<64>& v=fHuffmanCoders[memberId]->Encode((value>fHuffmanCoders[memberId]->GetMaxValue())?fHuffmanCoders[memberId]->GetMaxValue():value, length);
  if (length>0) {
    return OutputBits(v, length);
  }
  return false;
}

int AliHLTDataDeflaterHuffman::AddTrainingParameter(const char* name, unsigned bitLength)
{
  /// add a parameter definition for huffman training

  /// returns index in the array
  if (!fHuffmanCoderList) {
    fHuffmanCoderList=new TList;
    if (!fHuffmanCoderList) return -ENOMEM;
    // always set ownership for the new list since it is supposed to
    // contain only internal pointers
    fHuffmanCoderList->SetOwner();
  } else if (!fHuffmanCoderList->IsOwner()) {
    // not sure about the pointers which are already in the list
    if (fHuffmanCoderList->GetEntries()>0) {
      HLTWarning("skip setting ownership because list contains already %d object(s), possible memory leak at cleanup");
    } else {
      fHuffmanCoderList->SetOwner();
    }
  }
  AliHLTHuffman* pHuffman=new AliHLTHuffman(name, bitLength);
  if (!pHuffman) return -ENOMEM;
  fHuffmanCoderList->Add(pHuffman);

  fHuffmanCoders.push_back(pHuffman);
  return fHuffmanCoders.size()-1;
}

bool AliHLTDataDeflaterHuffman::AddTrainingValue( int memberId, AliHLTUInt64_t const & value )
{
  /// add a training value for the specified parameter
  if (memberId>=(int)fHuffmanCoders.size()) {
    return false;
  }
  return fHuffmanCoders[memberId]->AddTrainingValue(value);
}

const TList* AliHLTDataDeflaterHuffman::GenerateHuffmanTree()
{
  /// generate the huffman tree
  for (unsigned i=0; i<fHuffmanCoders.size(); i++) {
    if (!fHuffmanCoders[i]->GenerateHuffmanTree()) {
      HLTError("failed to generate huffman tree for parameter %s", fHuffmanCoders[i]->GetName());
    }
  }
  return fHuffmanCoderList;
}

void AliHLTDataDeflaterHuffman::Clear(Option_t * option)
{
  // internal cleanup

  AliHLTDataDeflater::Clear(option);
}

void AliHLTDataDeflaterHuffman::Print(Option_t *option) const
{
  // print info
  Print(cout, option);
}

void AliHLTDataDeflaterHuffman::Print(ostream& out, Option_t * /*option*/) const
{
  // print to stream
  out << "AliHLTDataDeflaterHuffman:" << endl;
  for (vector<AliHLTHuffman*>::const_iterator it=fHuffmanCoders.begin();
       it!=fHuffmanCoders.end(); it++) {
    (*it)->Print("short"); cout << endl;
  }
  if (fHuffmanCoders.size()==0 && fHuffmanCoderList && fHuffmanCoderList->GetEntries()>0) {
    TIter next(fHuffmanCoderList);
    TObject* pObj=NULL;
    while ((pObj=next())) {
      pObj->Print("short"); cout << endl;
    }
  }
}

TObject *AliHLTDataDeflaterHuffman::FindObject(const char *name) const
{
  /// find object: 'DeflaterConfiguration'
  if (strcmp(name, "DeflaterConfiguration")==0) {
    for (unsigned i=0; i<fHuffmanCoders.size(); i++) {
      if (!fHuffmanCoders[i]->GenerateHuffmanTree()) {
	HLTError("generation of huffman tree for parameter '%s' failed", fHuffmanCoders[i]->GetName());
      }
    }
    
    return fHuffmanCoderList;
  }

  return NULL;
}

void AliHLTDataDeflaterHuffman::SaveAs(const char *filename, Option_t *option) const
{
  /// save data according to option
  TString remainingOptions;
  bool bWriteHuffmanConf=false; // write the huffman configuration
  TString strOption(option);
  std::auto_ptr<TObjArray> tokens(strOption.Tokenize(" "));
  if (tokens.get()) {
    for (int i=0; i<tokens->GetEntriesFast(); i++) {
      if (!tokens->At(i)) continue;
      const char* key="";
      TString arg=tokens->At(i)->GetName();

      key="huffmanconf";
      if (arg.BeginsWith(key)) {
	bWriteHuffmanConf=true;
	continue;
      }

      if (!remainingOptions.IsNull()) remainingOptions+=" ";
      remainingOptions+=arg;
    }
  }

  if (bWriteHuffmanConf) {
    std::auto_ptr<TFile> output(TFile::Open(filename, "RECREATE"));
    if (!output.get() || output->IsZombie()) {
      HLTError("can not open file %s from writing", filename);
      return;
    }

    if (!fHuffmanCoderList || fHuffmanCoderList->GetEntries()==0) {
      HLTError("no huffman instances available for writing");
      return;
    }

    for (unsigned i=0; i<fHuffmanCoders.size(); i++) {
      if (!fHuffmanCoders[i]->GenerateHuffmanTree()) {
	HLTError("generation of huffman tree for parameter '%s' failed", fHuffmanCoders[i]->GetName());
      }
    }

    output->cd();
    fHuffmanCoderList->Write("DeflaterConfiguration", TObject::kSingleKey);
    output->Close();
  }

  return AliHLTDataDeflater::SaveAs(remainingOptions);
}

ostream& operator<<(ostream &out, const AliHLTDataDeflaterHuffman& me)
{
  me.Print(out);
  return out;
}

