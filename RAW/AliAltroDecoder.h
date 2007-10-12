//#-*- Mode: c++ -*-

#ifndef ALIALTRODECODER_H
#define ALIALTRODECODER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <TObject.h>
 
#define DDL_32BLOCK_SIZE         5
#define MAX_BRANCHES             2
#define MAX_FEE_PER_BRANCH       16
#define MAX_ALTROS_PER_FEE       8
#define CHANNELS_PER_ALTRO       16
#define MAX_SAMPLES_PER_CHANNEL  1024
#define ALTRO_TRAILER_SIZE       4
#define MAX_TRAILER_WORDS        3

class AliAltroData;

class AliAltroDecoder: public TObject {
 public:
 
  /*
   *Default constructor
   **/
  AliAltroDecoder();

  /*
   *Default destructor
   **/
  virtual ~AliAltroDecoder();

  /*
   *Decode the RCU/DDL payload 
   **/
  Bool_t Decode();

  /*
   *Reads the next altro channels 
   **/
  Bool_t NextChannel(AliAltroData *altroDataPtr);

  /* 
   * DONT use !
   * For debugging purphoses only, will be removed in near future
   **/
  template<typename T> 
  void  DumpData(T *array, Int_t N, Int_t nPerLine)
  {
    cout <<   "DumpData N=  " << N <<endl;
    for(Int_t i= 0; i< N; i++)
      {
	if((i%nPerLine == 0)  &&  (i != 0))
	  {
	    printf("\n");
	  }

	cout << array[i]<< "\t";

      }
  }

  int SetMemory(UChar_t  *dtaPtr, UInt_t size);

  void PrintInfo(AliAltroData &altrodata, Int_t n = 0, Int_t nPerLine = 4);

  /*
   *Prints to stdout the percent of altroblocks that
   *is missing the 2aaa trailer.
   **/
  Float_t GetFailureRate();

 private:

  AliAltroDecoder& operator = (const AliAltroDecoder& decoder);
  AliAltroDecoder(const AliAltroDecoder& decoder);

  /*
   *Check wether or not there is consistency between the number of 40 bit altro words given by
   *the RCU payload and the number of 40 bit words calculated from the size of the RCU payload.
   **/
  Bool_t CheckPayloadTrailer();

  /*
   *Decode one 160 bit DDL block into 16 x 16 bit integers (only least significant 10 bits are filled)
   **/
  void DecodeDDLBlock();
 
  /*
   *Decode one 160 bit DDL block into 16 integers. 
   *In order to use the same decoding function (DecodeDDLBlock()) 
   *a copy of the the last DDL block is made and  
   *if the las block does not align with 160 bits then it is padded with zeroes 
   **/
  void DecodeLastDDLBlock();

  /*
   *Use for simulated data only.
   *Patch for incorrectly simulated data. Counts the number of 
   *2aaa word in the trailer of the payload and tries to figure out
   *the correct number of 40 bit altro words in the RCU pauload
   **/
  Int_t countAAApaddings();

  UInt_t  *f32DtaPtr;                        // Pointer to dat of the input buffer in entities of 32 bit words (the RCU/DDL block)
  UChar_t *f8DtaPtr;                         // Pointer to dat of the input buffer in entities of 8 bit words (the RCU/DDL block)
  const Long_t fN32HeaderWords;              // Number of 32 bit words in the common data header
  Int_t    fN40AltroWords;                   // Number of 40 bit altro words contained in the RCU payload as calculated form the payload size
  Int_t    fN40RcuAltroWords;                // Number of 40 bit altro words contained in the RCU payload as given by the RCU trailer        
  Int_t    fNDDLBlocks;                      // Number of DDL blocks in the payload (the last blocj might/ight not be 160 bits )
  Int_t    f32LastDDLBlockSize;              // Size of the last DDL block
  UInt_t   fDDLBlockDummy[DDL_32BLOCK_SIZE]; // buffer to contain the las DDL block, if the block is not aligned with 160 bitm the remaining fileds are padded with zeroes
  UInt_t   f32PayloadSize;                   // The size of the payload in entities of 32 bit words (after subtraction of the RCU header and the RCU trailer words)
  Long_t   fOutBufferIndex;                  // current buffer position of the buffer for the decoded data (10 bit words represnted as int's)
  UInt_t   fSize;                            // The size of the input RCU/DDL payload in entities of bytes, inluding the RCU header and trailer
  UInt_t   fOutBuffer[MAX_FEE_PER_BRANCH*MAX_BRANCHES*MAX_ALTROS_PER_FEE*CHANNELS_PER_ALTRO*(MAX_SAMPLES_PER_CHANNEL + ALTRO_TRAILER_SIZE)]; // Buffer to hold the decoded data
  UInt_t   fNAltro10bitWords;                // The total number of 10 bit altro words in the RCU payload, including trailers (disregardin that the altro trialer is not aligned with 10 bit)
  Int_t    fComplete;                        // Number of altro channels that is only partially read out  (0x2aaa pattern missing in trailer)
  Int_t    fInComplete;                      // Number of altro channels that is read out properly
  Bool_t   fDecodeIfCorruptedTrailer;        // Wether or not to try to decode the data if the RCU trailer is incorrect (will succseed in most cases)
  Bool_t   fIsDecoded;                       // Wether or not the buffer set last by the "SetMemory()" function has been decoded

  Bool_t  fIsFatalCorruptedTrailer;         //If trailer is fataly corrupted, not possible in any way to recover, then it is not allowed to decode the DDL payload.  

  ClassDef(AliAltroDecoder, 0)  // class for decoding Altro payload
};

#endif
