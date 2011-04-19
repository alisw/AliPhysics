//-*- Mode: C++ -*-
// $Id$
#ifndef ALIHLTCOMPHUFFMANCODEDATA_H
#define ALIHLTCOMPHUFFMANCODEDATA_H

///**************************************************************************
///* This file is property of and copyright by the ALICE HLT Project        * 
///* All rights reserved.                                                   *
///*                                                                        *
///* Primary Author: Jenny Wagner  (jwagner@cern.ch)                        *
///*                                                                        *
///* Permission to use, copy, modify and distribute this software and its   *
///* documentation strictly for non-commercial purposes is hereby granted   *
///* without fee, provided that the above copyright notice appears in all   *
///* copies and that both the copyright notice and this permission notice   *
///* appear in the supporting documentation. The authors make no claims     *
///* about the suitability of this software for any purpose. It is          * 
///* provided "as is" without express or implied warranty.                  *
///**************************************************************************

/// @file   AliHLTCOMPHuffmanCodeData.h
/// @author Jenny Wagner
/// @date   29-08-2007
/// @brief  Data class containing the Huffman code table of ADV-values
///

#include "AliHLTLogging.h"
#include "AliHLTDataTypes.h"

// typedef for the final HuffmanCodeArray
/** @class   AliHLTCOMPHuffmanCodeData
    @author Jenny Wagner
    @date   29-08-2007
    @brief  The Huffman Code Data containing the struture of the Huffman code table 
*/
class AliHLTCOMPHuffmanCodeData : public TObject, public AliHLTLogging
{
public:

  /** typedef for the structs in Huffman code table */
  typedef struct AliHLTCOMPHuffmanCodeStruct
  {
    Int_t famplitude; // ADC-value
    AliHLTUInt64_t fhuffmancode; // maximal possible codelength: 64 bits
    UInt_t fvalidcodelength; // validcodelength needed as code is of variable length!!
  } AliHLTCOMPHuffmanCodeStruct;
  
  /** standard constructor */
  AliHLTCOMPHuffmanCodeData();

  /** destructor */
  virtual ~AliHLTCOMPHuffmanCodeData();

  /** convert Huffman code struct into class instance of HuffmanCodeData */
  void SetHuffmanCodeData(const AliHLTCOMPHuffmanCodeStruct& codetableentry);

  /** return Huffman code struct */
  AliHLTCOMPHuffmanCodeStruct* GetHuffmanCodeData(AliHLTCOMPHuffmanCodeStruct* codetableentry) const;

private:

  /** copy constructor prohibited */
  AliHLTCOMPHuffmanCodeData(const AliHLTCOMPHuffmanCodeData&);

  /** assignment operator prohibited */
  AliHLTCOMPHuffmanCodeData& operator=(const AliHLTCOMPHuffmanCodeData&);

  Int_t famplitude;       // 10-bit ADC-value
  AliHLTUInt64_t fhuffmancode;   // respective Huffman code with maximal possible codelength: 64 bits
  Int_t fvalidcodelength; // variable to store the respective valid codelength
  
  ClassDef(AliHLTCOMPHuffmanCodeData, 1)
    
    };
#endif

