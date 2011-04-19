//-*- Mode: C++ -*-
// $Id$
#ifndef ALIHLTCOMPHUFFMANDATA_H
#define ALIHLTCOMPHUFFMANDATA_H

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

/// @file   AliHLTCOMPHuffmanData.h
/// @author Jenny Wagner
/// @date   29-08-2007
/// @brief  The Huffman Data containing the Huffman code table and the occurrence table of ADV-values
///

#include "AliHLTLogging.h"
#include "AliHLTCompDefinitions.h"
#include "AliHLTCOMPHuffmanOccurrenceData.h"
#include "AliHLTCOMPHuffmanCodeData.h"

#define TIMEBINS 1024

// type definitions needed for the Huffman compression
/** @class   AliHLTCOMPHuffmanData
    @author Jenny Wagner
    @date   29-08-2007
    @brief  The Huffman Data containing the Huffman code table and the occurrence table of ADC-values
            it uses the classes @ref AliHLTCOMPHuffmanOccurrenceData and @ref AliHLTCOMPHuffmanCodeData
            to convert all structs into classes (and instances) such that this information can be written
            to a ROOT file
*/


class AliHLTCOMPHuffmanData : public TObject, public AliHLTLogging
{
public:

  /** typedef for the nodes in the Huffman tree */
  typedef struct AliHLTCOMPHuffmanTreeDataStruct  // struct containing all information for Huffman tree
  {
    AliHLTCOMPHuffmanOccurrenceData::AliHLTCOMPHuffmanDataStruct fleafcontents;  // contains all info about one leaf (node)
    AliHLTCOMPHuffmanTreeDataStruct* fleft;      // pointer to left children
    AliHLTCOMPHuffmanTreeDataStruct* fright;     // pointer to right children
    AliHLTCOMPHuffmanTreeDataStruct* fparent;    // pointer to parent node
    
    AliHLTCOMPHuffmanTreeDataStruct* fnext;      // pointer to next element in remaining list to build tree from
    AliHLTCOMPHuffmanTreeDataStruct* fprevious;  // pointer to previous element in remaining list to build tree from
   
  } AliHLTCOMPHuffmanTreeDataStruct;

  /** standard constructor */
  AliHLTCOMPHuffmanData();

  /** destructor */
  virtual ~AliHLTCOMPHuffmanData();

  /** init internal occurrencetable and codetable
   * @param occurrencetable   pointer to occurrence table
   * @param codetable         pointer to Huffman code table
   */
  int InitHuffmanData(const AliHLTCOMPHuffmanOccurrenceData::AliHLTCOMPHuffmanDataStruct* occurrencetable, const AliHLTCOMPHuffmanCodeData::AliHLTCOMPHuffmanCodeStruct* codetable);

  /** copy occurrencetable to target
   * @param occurrencetable table  pointer to occurrence table
   */
  AliHLTCOMPHuffmanOccurrenceData::AliHLTCOMPHuffmanDataStruct* GetOccurrenceTable(AliHLTCOMPHuffmanOccurrenceData::AliHLTCOMPHuffmanDataStruct* occurrencetable) const;

  /** copy codetable to target
   * @param codetable   pointer to Huffman code table
  */
  AliHLTCOMPHuffmanCodeData::AliHLTCOMPHuffmanCodeStruct* GetCodeTable(AliHLTCOMPHuffmanCodeData::AliHLTCOMPHuffmanCodeStruct* codetable) const;

  /** set specifications for Huffman data to be recognised again in the Preprocessor before storing the data in the OCDB 
   * @param origin string to detector origin
   * @param dataspec integer to define data specification
   * @return 0 upon success
   **/
  Int_t SetOCDBSpecifications(const TString& origin, Int_t dataspec);

  /** get detector origin (used in Preprocessor of OCDB)
   * @return fOrigin (TString)
   **/
  TString GetOrigin() const {return fOrigin;}

  /** get dataspec (i.e. table number from 1 to 6) (used in Preprocessor of OCDB)
   * @return dataspec (Int_t)
   **/
  Int_t GetDataSpec() const {return fDataSpec;}

private:

  /** copy constructor prohibited */
  AliHLTCOMPHuffmanData(const AliHLTCOMPHuffmanData&);

  /** assignment operator prohibited */
  AliHLTCOMPHuffmanData& operator=(const AliHLTCOMPHuffmanData&);
  
  /** array of instances of HuffmanOccurrenceData that contains occurrence table */
  AliHLTCOMPHuffmanOccurrenceData fOccurrenceTable[TIMEBINS]; // occurrence table for all ADC-values
  /** array of instances of HuffmanCodeData thtat contains complete Huffman code table */
  AliHLTCOMPHuffmanCodeData fCodeTable[TIMEBINS];             // Huffman translation table for all ADC-values

 /** define detector where this data comes from **/
  TString fOrigin;

  /** define dataspecification where this data comes from **/
  Int_t fDataSpec;
    
    
  ClassDef(AliHLTCOMPHuffmanData, 1)
    
    };
#endif

