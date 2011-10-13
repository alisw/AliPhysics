//-*- Mode: C++ -*-
// $Id$
#ifndef ALIHLTDATAINFLATERHUFFMAN_H
#define ALIHLTDATAINFLATERHUFFMAN_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/// @file   AliHLTDataInflaterHuffman.h
/// @author Matthias Richter
/// @date   2011-09-01
/// @brief  Data inflater implementation for huffman encoded data
/// @note   

#include "AliHLTDataInflater.h"

class AliHLTHuffman;
class TList;

class AliHLTDataInflaterHuffman : public AliHLTDataInflater
{
public:
  /// standard constructor
  AliHLTDataInflaterHuffman();
  /// destructor
  ~AliHLTDataInflaterHuffman();

  /// add a parameter definition to the configuration, return reference id
  int AddParameterDefinition(const char* name, unsigned bitLength);

  /// init list of decoders
  int InitDecoders(TList* decoderlist);

  /// overloaded from AliHLTDataInflater
  virtual bool NextValue(AliHLTUInt64_t& value, AliHLTUInt32_t& length);
  /// switch to next parameter
  virtual int NextParameter() {
    if (fHuffmanCoders.size()==0) return -1;
    if (fLegacyMode>0) return fCurrentParameter;
    fLegacyMode=0;
    if ((++fCurrentParameter)>=(int)fHuffmanCoders.size()) fCurrentParameter=0;
    return fCurrentParameter;
  }

protected:
private:
  /** copy constructor prohibited */
  AliHLTDataInflaterHuffman(const AliHLTDataInflaterHuffman&);
  /** assignment operator prohibited */
  AliHLTDataInflaterHuffman& operator=(const AliHLTDataInflaterHuffman&);

  /// index of the decoders in the decoder list
  vector<AliHLTHuffman*> fHuffmanCoders; //! index of decoders

  /// list of huffman coders identified by parameter name
  TList* fHuffmanCoderList; //! list of huffman coders

  /// current parameter during reading
  int fCurrentParameter; //!
  /// legacy mode to handle code not using NextParameter()
  int fLegacyMode;

  ClassDef(AliHLTDataInflaterHuffman, 0)
};

#endif //ALIHLTDATAINFLATERHUFFMAN_H
