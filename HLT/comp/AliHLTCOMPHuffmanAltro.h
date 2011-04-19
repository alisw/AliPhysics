// $Id$

#ifndef ALIHLTCOMPHUFFMANALTRO_H
#define ALIHLTCOMPHUFFMANALTRO_H

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

/** @file   AliHLTCOMPHuffmanAltro.h
    @author Jenny Wagner
    @date   20-11-2007
    @brief  The Huffman compressor
*/

#include "AliHLTLogging.h"
#include "AliHLTCOMPHuffmanData.h"

class AliRawReaderMemory;
class AliAltroRawStreamV3;

/** @class   AliHLTCOMPHuffmanAltro
    @author Jenny Wagner
    @date   20-11-2007
    @brief  The Huffman Compressor with functions training (for Calibration), compress and decompress, calculate entropy  
*/
class AliHLTCOMPHuffmanAltro : public AliHLTLogging
{
 public:

  /** constructor for test use only */
  AliHLTCOMPHuffmanAltro();

  /** standard constructor
   *
   * @param compressionswitch   decides whether to compress or decompress (TRUE for calibration/training)
   * @param trainingmode        decides whether to create a new Huffman table 
   * @param translationtable    pointer to lookup Huffman table for compression and decompression
   * @param nrcutrailerwords    number of NRCU trailer words (ranging from 1 to 3)
   */
  AliHLTCOMPHuffmanAltro(Bool_t compressionswitch, Bool_t trainingmode, AliHLTCOMPHuffmanCodeData::AliHLTCOMPHuffmanCodeStruct* translationtable, Int_t nrcutrailerwords);
  /** destructor */
  virtual ~AliHLTCOMPHuffmanAltro();

  /** SetInputData 
   * Data buffer has to be valid all the time, no internal copy.
   *  @param memory       pointer to input data
   *  @param size         size of input data
   *  @param equipmentId  equipment (ddl) id
   */
  int AddInputData(UChar_t* memory, ULong_t size, Int_t equipmentId);
  
  /** SetOutputData
   *  @param outputdata  pointer to output data
   *  @param outputsize  size of output data
   */
  int SetOutputData(AliHLTUInt8_t* outputdata, unsigned long outputsize);

  /**
   * Reset and prepare for new event
   */
  int Reset();

  /** GetOutputDataSize (which is unknown at first and has to be calculated
   * @return output data size
  */
  unsigned long GetOutputDataSize();

  /** initialise training table */
  int InitNewTrainingTable();

  /** write out new HuffmanData
   * @param huffmandata   pointer to Huffman data 
  */
  int SetHuffmanData(AliHLTCOMPHuffmanData* huffmandata) const;

  /** get training table from HuffmanData instance
   * @param huffmandata  pointer to Huffman data
  */
  int GetTrainingTable(const AliHLTCOMPHuffmanData* huffmandata);

  /** initialise the translation table
   * @param huffmandata  pointer to Huffman data
  */
  int GetTranslationTable(const AliHLTCOMPHuffmanData* huffmandata);

  /** function to compress or decompress data */
  void ProcessData();

 /** function to create Huffman code table by means of private functions
  * @return zero upon success
 */
  Int_t CreateCodeTable();

  /** calculate the entropy of the input data
   * @return calculated entropy (in double precision)
  */
  Double_t GetEntropy();

  /** function to read 10 bit data in and write same 10 bit data out (performance test) */
  Int_t CopyData();

  /** print content */
  virtual void Print(Option_t* option = "") const;

 private:

  /** copy constructor prohibited */
  AliHLTCOMPHuffmanAltro(const AliHLTCOMPHuffmanAltro&);
  /** assignment operator prohibited */
  AliHLTCOMPHuffmanAltro& operator=(const AliHLTCOMPHuffmanAltro&);

  /** function to calculate the entropy of the incoming data */
  /// FIXME: the size of the array is not known, everything assumes size TIMEBINS
  Int_t CalcEntropy(const AliHLTCOMPHuffmanOccurrenceData::AliHLTCOMPHuffmanDataStruct* occurrencetable=NULL);

  /** function for merge sorting the array data 
   * @param unsortedarray   unsorted array of data from occurrence table
   * @param n               number of entries in the unsorted array
   * @return pointer to first element of sorted list
   */
  AliHLTCOMPHuffmanOccurrenceData::AliHLTCOMPHuffmanDataStruct* Mergesort(AliHLTCOMPHuffmanOccurrenceData::AliHLTCOMPHuffmanDataStruct* unsortedarray, Int_t n);

  /** function to create the HuffmanCode and write it in a sorted array 
   * @param treeroot           pointer to root of Huffman tree
   * @param HuffmanCodearray   pointer to final Huffman code table (array)
   * @return zero upon success
   */ 
 Int_t HuffmanCode(AliHLTCOMPHuffmanData::AliHLTCOMPHuffmanTreeDataStruct* treeroot, AliHLTCOMPHuffmanCodeData::AliHLTCOMPHuffmanCodeStruct* HuffmanCodearray);
  
  /** function to create the Huffmantree
   * @param listroot   pointer to first element of sorted list to build Huffman tree from
   * @param listend    pointer to last element of sorted list to build Huffman tree from
   * @param n          number of entries in the list
   * @return pointer to root of Huffman tree
   */
 AliHLTCOMPHuffmanData::AliHLTCOMPHuffmanTreeDataStruct* CreateHuffmanTree(AliHLTCOMPHuffmanData::AliHLTCOMPHuffmanTreeDataStruct* listroot,  AliHLTCOMPHuffmanData::AliHLTCOMPHuffmanTreeDataStruct* listend, Int_t n);
  
  /** entropy encoding function: read in data, table lookup to encode, write encoded data to output array and set fOutputDataSize
   * @return zero upon success
   */
 Int_t EntropyCompression();

  /** merge sorting function for translation table
   * @param unsortedarray   Huffman code table which is not sorted for decompression yet
   * @param n               number of entries in the unsorted array
   * @return pointer to first element of sorted list
   */
 AliHLTCOMPHuffmanCodeData::AliHLTCOMPHuffmanCodeStruct* TTMergesort(AliHLTCOMPHuffmanCodeData::AliHLTCOMPHuffmanCodeStruct* unsortedarray, Int_t n);
    
  /** entropy decoding function
   * @return zero upon success
   */
  Int_t EntropyDecompression();
  
  /** training function for the translation table
   * @return zero upon success
   */
  Int_t TrainingData(); 

  /// raw reader instance
  AliRawReaderMemory* fpRawReader;                    //! transient
  /// raw stream decoder
  AliAltroRawStreamV3* fpAltroRawStream;          //! transient

  /** boolean to decide whether to process data (FALSE) or create a new code table (TRUE) */
  Bool_t fTrainingMode;                         // choice if new codetable is created or not
  /** boolean to decide whether to compress (TRUE) or decompress (FALSE) incoming data (automatically TRUE for code creation) */
  Bool_t fCompressionSwitch;                    // mode choice (training, compress, decompress)
  /** pointer to input data */
  AliHLTUInt8_t* fPointer2InData;               // pointer to input data (uncompressed raw data)
  /** pointer to output data */
  AliHLTUInt64_t* fPointer2OutData;             // pointer to output data (compressed data)
  /** input data size */
  unsigned long fInputDataSize;                 // input data size
  /** output data size */
  UInt_t fOutputDataSize;                        // output data size
  /** number of NRCU trailer words */
  UInt_t fNrcuTrailerwords;                      // number of RCU trailerwords
  /** calculated entropy of input data */
  Double_t fEntropy;                            // entropy of the file

  int fVerbosity;                               //! verbosity

  //AliHLTUInt8_t fSlice;                         // transient variables to specify
  //AliHLTUInt8_t fPatch;                         // where output root file comes from
 
  /** pointer to occurrence table */
  AliHLTCOMPHuffmanOccurrenceData::AliHLTCOMPHuffmanDataStruct* fTrainingTable;       // training table with amplitudes and resp. abundances
  /** pointer to Huffman code table */
  AliHLTCOMPHuffmanCodeData::AliHLTCOMPHuffmanCodeStruct* fTranslationTable;    // Huffman Entropy Code Table


  ClassDef(AliHLTCOMPHuffmanAltro, 0)
	};
#endif
