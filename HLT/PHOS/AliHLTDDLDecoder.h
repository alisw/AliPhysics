#ifndef ALIHLTDDLDECODER_H
#define ALIHLTDDLDECODER_H

/**************************************************************************
 * This file is property of and copyright by the Experimental Nuclear     *
 * Physics Group, Dep. of Physics                                         *
 * University of Oslo, Norway, 2007                                       *
 *                                                                        *
 * Authors: Øystein Djuvsland <oystein.djuvsland@gmail.com>               * 
 *      &   Per Thomas Hille <perthi@fys.uio.no>                          *
 * for the ALICE HLT Project.                                             *
 * Contributors are mentioned in the code where appropriate.              *
 * Please report bugs to perthi@fys.uio.no                                * 
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#define DDL_32BLOCK_SIZE 5

#include "Rtypes.h"
#include <iostream>
#include "AliHLTPHOSConstants.h"

using  std::cout;
using  std::endl;
using namespace PhosHLTConst;


class AliHLTAltroData;


class AliHLTDDLDecoder
{
 public:
 
  /*
   *Default constructor
   */
  AliHLTDDLDecoder();

  /*
   *Default destructor
   */
  virtual ~AliHLTDDLDecoder();

  /*
   *Check wether or not there is consistency between the number of 40 bit altro words given by
   *the RCU payload and the number of 40 bit words calculated from the size of the RCU payload.
   */
  bool CheckPayload();

  /*
   *Decode the RCU/DDL payload 
   */
  bool Decode();

  /*
   *Reads the next altro channels 
   */
  bool NextChannel(AliHLTAltroData *altroDataPtr);

  /* 
   * DONT use !
   * For debugging purphoses only, will be removed in near future
   */
  template<typename T> 
  void  DumpData(T *array, int N, int nPerLine)
  {
    cout <<   "DumpData N=  " << N <<endl;
    for(int i= 0; i< N; i++)
      {
	if((i%nPerLine == 0)  &&  (i != 0))
	  {
	    printf("\n");
	  }

	cout << array[i]<< "\t";

      }
  }


  void SetMemory(UChar_t  *dtaPtr, UInt_t size);
  
  void PrintInfo(AliHLTAltroData &altrodata, int n = 0, int nPerLine = 4);

  /*
   * Prints to stdout the percent of altroblocks that
   * is missing the 2aaa trailer.
   */
  float GetFailureRate();

 private:
  /*
   *Decode one 160 bit DDL block into 16 x 16 bit integers (only least significant 10 bits are filled)
   */
  void DecodeDDLBlock();
 
  /*
   *Decode one 160 bit DDL block into 16 integers, 
   *if the las blaock does not align with 160 bits then first pad with zeroes 
  */
  void DecodeLastDDLBlock();

  /*
   *Use for simulated data only.
   *Patch for incorrectly simulated data. Counts the number of 
   *aaaa word in the trailer of the payload and tries to figure out
   *the correct number of 40 bit altro words in the RCU pauload
   */
  int countAAApaddings();

  UInt_t  *f32DtaPtr;                     /**<Pointer to dat of the input buffer in entities of 32 bit words (the RCU/DDL block) */
  UChar_t *f8DtaPtr;                      /**<Pointer to dat of the input buffer in entities of 8 bit words (the RCU/DDL block) */
  const long int fN32HeaderWords;         /**<Number of 32 bit words in the common data header*/
  const long int fN32RcuTrailerWords;     /**<Number of 32 bit words in the RCU trailer*/
  int  fN40AltroWords;                    /**<Number of 40 bit altro words contained in the RCU payload as calculated form the payload size*/
  int  fN40RcuAltroWords;                 /**<Number of 40 bit altro words contained in the RCU payload as given by the RCU trailer*/            
  int fNDDLBlocks;                        /**<Number of DDL blocks in the payload (the last blocj might/ight not be 160 bits )*/
  int  f32LastDDLBlockSize;               /**<Size of the last DDL block*/ 
  UInt_t fDDLBlockDummy[DDL_BLOCK_SIZE];  /**<buffer to contain the las DDL block, if the block is not aligned with 160 bitm the remaining fileds are padded with zeroes*/
  UInt_t f32PayloadSize;                  /**<The size of the payload in entities of 32 bit words (after subtraction of the RCU header and the RCU trailer words)*/
  long int fOutBufferIndex;               /**<current buffer position of the buffer for the decoded data (10 bit words represnted as int's)*/
  UInt_t  fSize;                          /**<The size of the input RCU/DDL payload in entities of bytes, inluding the RCU header and trailer */  
  UInt_t fOutBuffer[N_FEECS*N_BRANCHES*8*N_ALTROCHANNELS*(ALTRO_MAX_SAMPLES + ALTRO_MAX_TRALER_SIZE)];  /**<Buffer to hold the decoded data*/
  UInt_t fNAltro10bitWords;               /**<The total number of 10 bit altro words in the RCU payload, including trailers (disregardin that the altro trialer is not aligned with 10 bit)*/
  int fComplete;                          /**<Number of altro channels that is only partially read out  (0x2aaa pattern missing in trailer)*/
  int fInComplete;                        /**<Number of altro channels that is read out properly*/       
  bool fDecodeIfCorruptedTrailer;         /**<Wether or not to try to decode the data if the RCU trailer is incorrect (will succseed in most cases)*/ 
  bool fIsDecoded;                        /**<Wether or not the buffer set last by the "SetMemory()" function has been decoded*/
};

#endif
