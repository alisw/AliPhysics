//#-*- Mode: c++ -*-
// $Id$

#ifndef ALIALTRODECODER_H
#define ALIALTRODECODER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/**
   @file AliAltroDocoder.h
   @author Per Thomas Hille, Oystein Djuvsland
   @date   
   @brief High performance decoder class for the RCU/Altro data format
*/

///////////////////////////////////////////////////////////////////////////////
///
/// This is the class for fast decoding of TPC/PHOS/EMCAL raw data
//  see .cxx file for more detailed comments.
///
///////////////////////////////////////////////////////////////////////////////


#define DECODERERROR -3

#include <TObject.h>

#include <iostream>
using namespace std;
 
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
  AliAltroDecoder();
  virtual ~AliAltroDecoder();
  Bool_t Decode();
  Bool_t NextChannel(AliAltroData *altroDataPtr);

  /**
   * Copy the original 10/40 bit encecoded data of the current channel.
   * The funtions copies the data to the end of the provided buffer.
   * @param pBuffer    target buffer
   * @param bufferSize size of target buffer
   * @return number of copied bytes, neg. error code if failed
   */
  Int_t CopyBackward(Byte_t* pBuffer, Int_t bufferSize);

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
  Float_t GetFailureRate();

  /**
   * Provide a pointer to RCU trailer.
   * The type of the parameter might not be optimal, but the function has
   * been chosen that way to be similar to the counterpart in
   * AliAltroRawStream.
   * @return kTRUE if trailer available;
   */
  Bool_t  GetRCUTrailerData(UChar_t*& data) const;
  Int_t   GetRCUTrailerSize() const;

 private:

  AliAltroDecoder& operator = (const AliAltroDecoder& decoder);
  AliAltroDecoder(const AliAltroDecoder& decoder);
  Bool_t CheckPayloadTrailer() const;
  void DecodeDDLBlock();
  void DecodeLastDDLBlock();
  Int_t CountAAApaddings() const;
  UInt_t  *f32DtaPtr;                        // Pointer to dat of the input buffer in entities of 32 bit words (the RCU/DDL block)
  UChar_t *f8DtaPtr;                         // Pointer to dat of the input buffer in entities of 8 bit words (the RCU/DDL block)
  const Long_t fkN32HeaderWords;              // Number of 32 bit words in the common data header
  Int_t    fN40AltroWords;                   // Number of 40 bit altro words contained in the RCU payload as calculated form the payload size
  Int_t    fN40RcuAltroWords;                // Number of 40 bit altro words contained in the RCU payload as given by the RCU trailer        
  Int_t    fNDDLBlocks;                      // Number of DDL blocks in the payload (the last blocj might/ight not be 160 bits )
  Int_t    f32LastDDLBlockSize;              // Size of the last DDL block
  UInt_t   fDDLBlockDummy[DDL_32BLOCK_SIZE]; // buffer to contain the las DDL block, if the block is not aligned with 160 bitm the remaining fileds are padded with zeroes
  UInt_t   f8PayloadSize;                    // The size of the payload in bytes (after subtraction of the RCU header and the RCU trailer words)
  Long_t   fOutBufferIndex;                  // current buffer position of the buffer for the decoded data (10 bit words represnted as int's)
  UInt_t   fSize;                            // The size of the input RCU/DDL payload in entities of bytes, inluding the RCU header and trailer
  UInt_t   fOutBuffer[MAX_FEE_PER_BRANCH*MAX_BRANCHES*MAX_ALTROS_PER_FEE*CHANNELS_PER_ALTRO*(MAX_SAMPLES_PER_CHANNEL + ALTRO_TRAILER_SIZE)]; // Buffer to hold the decoded data
  UInt_t   fNAltro10bitWords;                // The total number of 10 bit altro words in the RCU payload, including trailers (disregardin that the altro trialer is not aligned with 10 bit)
  Int_t    fComplete;                        // Number of altro channels that is only partially read out  (0x2aaa pattern missing in trailer)
  Int_t    fInComplete;                      // Number of altro channels that is read out properly
  Bool_t   fDecodeIfCorruptedTrailer;        // Wether or not to try to decode the data if the RCU trailer is incorrect (will succseed in most cases)
  Bool_t   fIsDecoded;                       // Wether or not the buffer set last by the "SetMemory()" function has been decoded
  Bool_t   fIsFatalCorruptedTrailer;          // If trailer is fataly corrupted, not possible in any way to recover, then it is not allowed to decode the DDL payload.  

  //  Bool_t   fIsFirstChannelOfPayload;


  ClassDef(AliAltroDecoder, 0)  // class for decoding Altro payload
};

#endif
