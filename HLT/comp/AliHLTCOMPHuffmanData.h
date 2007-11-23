//-*- Mode: C++ -*-

#ifndef ALIHLTCOMPHUFFMANDATA_H
#define ALIHLTCOMPHUFFMANDATA_H

/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * All rights reserved.                                                   *
 *                                                                        *
 * Primary Author: Jenny Wagner  (jwagner@cern.ch)                        *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          * 
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/** @file   AliHLTCOMPHuffmanData.h
    @author Jenny Wagner
    @date   29-08-2007
    @brief  The Huffman Data containing the Huffman code table and the occurrence table of ADV-values
*/

#include "AliHLTLogging.h"
#include "AliHLTDataTypes.h"
#include "AliHLTCOMPHuffmanOccurrenceData.h"
#include "AliHLTCOMPHuffmanCodeData.h"

// type definitions needed for the Huffman compression
/** @class   AliHLTCOMPHuffmanData
    @author Jenny Wagner
    @date   29-08-2007
    @brief  The Huffman Data containing the Huffman code table and the occurrence table of ADC-values
*/

#define TIMEBINS 1024

class AliHLTCOMPHuffmanData : public TObject, public AliHLTLogging
{
public:

  /** typedef for the nodes in the Huffman tree */
  typedef struct AliHLTCOMPHuffmanTreeData_t
  {
    AliHLTCOMPHuffmanOccurrenceData::AliHLTCOMPHuffmanData_t leafcontents;
    AliHLTCOMPHuffmanTreeData_t* left;
    AliHLTCOMPHuffmanTreeData_t* right;
    AliHLTCOMPHuffmanTreeData_t* parent;
    
    AliHLTCOMPHuffmanTreeData_t* next;
    AliHLTCOMPHuffmanTreeData_t* previous;
   
  } AliHLTCOMPHuffmanTreeData_t;
  
  
  /** standard constructor */
  AliHLTCOMPHuffmanData();

  /** destructor */
  virtual ~AliHLTCOMPHuffmanData();

  /** get data from OCDB (currently from ROOT-file) and write into instance of HuffmanData 
   * @param occurrencetable   pointer to occurrence table
   * @param codetable         pointer to Huffman code table
   */
  void InitHuffmanData(AliHLTCOMPHuffmanOccurrenceData::AliHLTCOMPHuffmanData_t* occurrencetable, AliHLTCOMPHuffmanCodeData::AliHLTCOMPHuffmanCode_t* codetable);

  /** return occurrence table to be used/written somewhere else (intialisation of Huffman Compressor tables)
   * @param occurrence table  pointer to occurrence table
   */
  AliHLTCOMPHuffmanOccurrenceData::AliHLTCOMPHuffmanData_t* GetOccurrenceTable(AliHLTCOMPHuffmanOccurrenceData::AliHLTCOMPHuffmanData_t* occurrencetable);
  
  /** return code table to be used/ written somewhere else (initialisation of the Huffman Compressor tables)
   * @param codetable   pointer to Huffman code table
  */
  AliHLTCOMPHuffmanCodeData::AliHLTCOMPHuffmanCode_t* GetCodeTable(AliHLTCOMPHuffmanCodeData::AliHLTCOMPHuffmanCode_t* codetable);

private:

  /** copy constructor prohibited */
  AliHLTCOMPHuffmanData(const AliHLTCOMPHuffmanData&);

  /** assignment operator prohibited */
  AliHLTCOMPHuffmanData& operator=(const AliHLTCOMPHuffmanData&);
  
  /** array of instances of HuffmanOccurrenceData that contains occurrence table */
  AliHLTCOMPHuffmanOccurrenceData fOccurrenceTable[TIMEBINS]; // occurrence table for all ADC-values
  /** array of instances of HuffmanCodeData thtat contains complete Huffman code table */
  AliHLTCOMPHuffmanCodeData fCodeTable[TIMEBINS];             // Huffman translation table for all ADC-values
    
  ClassDef(AliHLTCOMPHuffmanData, 1)
    
    };
#endif

