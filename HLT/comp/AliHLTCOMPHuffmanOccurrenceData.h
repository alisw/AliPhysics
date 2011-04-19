//-*- Mode: C++ -*-
// $Id$
#ifndef ALIHLTCOMPHUFFMANOCCURRENCEDATA_H
#define ALIHLTCOMPHUFFMANOCCURRENCEDATA_H

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

/// @file   AliHLTCOMPHuffmanOccurrenceData.h
/// @author Jenny Wagner
/// @date   29-08-2007
/// @brief  Data class containing the occurrence table of ADC-values
///

#include "AliHLTLogging.h"
#include "AliHLTDataTypes.h"

// type definitions needed for the Huffman compression

/** @class   AliHLTCOMPHuffmanOccurrenceData
    @author Jenny Wagner
    @date   29-08-2007
    @brief  The Huffman Occurrence Data containing the struture of the ADC-value occurrence table 
*/
class AliHLTCOMPHuffmanOccurrenceData : public TObject, public AliHLTLogging
{
public:

  /** typedef for the array data of occurrence data */
  typedef struct
  {
    Int_t famplitude; // ADC-value
    Int_t fabundance; // occurrence = abundance of ADV-value
    Int_t fcode;       // internal boolean used to decide whether 0 or 1 is appended to Huffman code
  } AliHLTCOMPHuffmanDataStruct;
  
  /** standard constructor */
  AliHLTCOMPHuffmanOccurrenceData();

  /** destructor */
  virtual ~AliHLTCOMPHuffmanOccurrenceData();

  /** convert one entry of occurrence table into a class instance of HuffmanOccurrenceData */
  void SetHuffmanOccurrenceData(const AliHLTCOMPHuffmanDataStruct& occurrencetableentry);

  /** return one entry of occurrence table */
  AliHLTCOMPHuffmanDataStruct* GetHuffmanOccurrenceData(AliHLTCOMPHuffmanDataStruct* occurrencetableentry) const;

  int GetAmplitude() const {return famplitude;}
  int GetAbundance() const {return fabundance;}

private:

   /** copy constructor prohibited */
  AliHLTCOMPHuffmanOccurrenceData(const AliHLTCOMPHuffmanOccurrenceData&);

  /** assignment operator prohibited */
  AliHLTCOMPHuffmanOccurrenceData& operator=(const AliHLTCOMPHuffmanOccurrenceData&);

  /** 10-bit ADC value used for conversion from struct to class */
  Int_t famplitude;  // 10-bit ADC-value
  /** occurrence = abundance used for conversion from struct to class */
  Int_t fabundance; // occurrence of one 10-bit ADC-value
  Int_t fcode;       // internal variable used for sorting the binary tree (nothing to do with Huffman code!)
  
  ClassDef(AliHLTCOMPHuffmanOccurrenceData, 1)
    
    };
#endif

