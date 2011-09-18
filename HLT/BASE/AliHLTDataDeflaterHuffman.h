//-*- Mode: C++ -*-
// $Id$
#ifndef ALIHLTDATADEFLATERHUFFMAN_H
#define ALIHLTDATADEFLATERHUFFMAN_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/// @file   AliHLTDataDeflaterHuffman.h
/// @author Thorsten Kollegger, Matthias Richter
/// @date   2011-08-10
/// @brief  Data deflater class using huffman coding

#include "AliHLTDataDeflater.h"
#include "AliHLTHuffman.h"
#include <vector>
#include <string>
#include <sstream>

/**
 * @class AliHLTDataDeflaterHuffman
 * Deflater implementation using standard huffman code.
 *
 * @ingroup alihlt_base
 */
class AliHLTDataDeflaterHuffman : public AliHLTDataDeflater
{
public:
  /// standard constructor
  AliHLTDataDeflaterHuffman(bool bTrainingMode=false);
  /// destructor
  ~AliHLTDataDeflaterHuffman();

  /// add a parameter definition to the configuration, return reference id
  int AddParameterDefinition(const char* name, unsigned bitLength);

  /// init list of decoders
  int InitDecoders(TList* decoderlist);

  /// inherited from AliHLTDataDeflater: write bit pattern according to configuration
  virtual bool OutputParameterBits( int parameterId, AliHLTUInt64_t const & value );

  /// add a parameter definition for huffman training
  int AddTrainingParameter(const char* name, unsigned bitLength);

  /// add a training value for the specified parameter
  bool AddTrainingValue( int memberId, AliHLTUInt64_t const & value );

  /// generate huffman trees for all parameters and return list
  const TList* GenerateHuffmanTree();

  const TList* GetList() const {return fHuffmanCoderList;}

  /// clear the object and reset pointer references
  virtual void Clear(Option_t * /*option*/ ="");

  /// print info
  virtual void Print(Option_t *option="") const;

  /// print info
  virtual void Print(ostream& out, Option_t *option="") const;

  /// find object: 'DeflaterConfiguration'
  virtual TObject *FindObject(const char *name) const;

  /// save data according to option
  virtual void SaveAs(const char *filename="",Option_t *option="") const;

  /// DataDeflaterHuffman version (ID) is 2
  virtual int GetDeflaterVersion() const {return 2;}

  /// check if in training mode
  bool IsTrainingMode() const {return fTrainingMode;}

protected:
private:
  /// copy constructor prohibited
  AliHLTDataDeflaterHuffman(const AliHLTDataDeflaterHuffman&);
  /// assigment operator prohibited
  AliHLTDataDeflaterHuffman& operator=(const AliHLTDataDeflaterHuffman&);

  /// index of the decoders in the decoder list
  vector<AliHLTHuffman*> fHuffmanCoders; //! index of decoders

  /// list of huffman coders identified by parameter name
  TList* fHuffmanCoderList; //! list of huffman coders

  bool fTrainingMode; //! indicate training mode

  ClassDef(AliHLTDataDeflaterHuffman, 0)
};

ostream& operator<<(ostream &out, const AliHLTDataDeflaterHuffman& me);

#endif
