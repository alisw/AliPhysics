// XEmacs -*-C++-*-
// $Id$

#ifndef ALIHLTTPCCOMPMODELINFLATER_H
#define ALIHLTTPCCOMPMODELINFLATER_H
/* TPCCompModelInflaterright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full TPCCompModelInflaterright notice                               */

/** @file   AliHLTTPCCompModelInflater.h
    @author Timm Steinbeck
    @date   
    @brief  Declaration of an uncompressor class for TPC track data. */

#include "AliHLTTPCTrackArray.h"
#include "AliHLTTPCTrackletDataFormat.h"
#include "AliHLTTPCClusterDataFormat.h"
#include "AliHLTLogging.h"
#include "AliHLTTPCModels.h"

/**
 * @class AliHLTTPCCompModelInflater
 * @brief A HLT TPC track compressor. 
 *
 * A class that can compress HLT TPC track model data
 */
class AliHLTTPCCompModelInflater: public AliHLTLogging
    {
    public:
      /** standard constructor */
      AliHLTTPCCompModelInflater();
      /** standard destructor */
      virtual ~AliHLTTPCCompModelInflater();

      /** function to decompress tracks 
       * @param inData     AliHLTUInt8_t* pointer to input
       * @param inputSize  UInt_t const& input data size
       * @param output     AliHLTUInt8_t* pointer to output
       * @param outputSize UInt_t& output data size
       * @return zero upon success
       */
	int DecompressTracks( AliHLTUInt8_t* inData, UInt_t const& inputSize, AliHLTUInt8_t* output, UInt_t& outputSize );

      /** function to decompress remaining clusters
       * @param inData     AliHLTUInt8_t* pointer to input
       * @param inputSize  UInt_t const& input data size
       * @param output     AliHLTUInt8_t* pointer to output
       * @param outputSize UInt_t& output data size
       * @return zero upon success
       */
	int DecompressRemainingClusters( AliHLTUInt8_t* inData, UInt_t const& inputSize, AliHLTUInt8_t* output, UInt_t& outputSize );

    protected:

      /** function to initialise the bit data input
       * @param input  AliHLTUInt8_t* pointer to input data
       * @param inputSize UInt_t input data size
       */
      void InitBitDataInput( AliHLTUInt8_t* input, UInt_t inputSize )
      {
	fBitDataCurrentWord = 0;
	fBitDataCurrentPosInWord = 7;
	fBitDataCurrentInput = fBitDataCurrentInputStart = input;
	fBitDataCurrentInputEnd = input+inputSize;
	fBitDataCurrentWord = *fBitDataCurrentInput;
      }

      /** function to close bit data input */	
      void CloseBitDataInput()
      {
      }
      
      /** function to get current byte input position 
       * @return unsigned long value of current byte input position
       */
      unsigned long GetCurrentByteInputPosition() const
      {
	return (unsigned long)( fBitDataCurrentInput - fBitDataCurrentInputStart );
      }

      /** function to get current bit input position
       * @return unsigned value of current bit input position
       */
      unsigned GetCurrentBitInputPosition() const
      {
	return fBitDataCurrentPosInWord;
      }

      /** function to get current input byte
       * @return AliHLTUInt8_t value of current input byte
       */
      AliHLTUInt8_t GetCurrentInputByte() const
      {
		return fBitDataCurrentWord;
      }

      /** function to determine end of bit input
       * @return boolean if end is reached or not
       */
      bool EndOfBitInput() const
      {
	return (fBitDataCurrentInput>=fBitDataCurrentInputEnd);
      }
      
      /** function to get bit data input size bytes
       * @return UInt_t value of bit data input size bytes
       */
      UInt_t GetBitDataInputSizeBytes() const
      {
	return fBitDataCurrentInput-fBitDataCurrentInputStart;
      }

      /** function to determine input bit
       * @return boolean (if bit is 1 or 0)
       */
      bool InputBit( AliHLTUInt8_t & value )
      {
	if ( fBitDataCurrentInput>=fBitDataCurrentInputEnd )
	  return false;
	value = (fBitDataCurrentWord >> fBitDataCurrentPosInWord) & 1;
	if ( fBitDataCurrentPosInWord )
	  fBitDataCurrentPosInWord--;
	else
	  {
	    fBitDataCurrentInput++;
	    if ( fBitDataCurrentInput<fBitDataCurrentInputEnd )
	      {
		fBitDataCurrentWord = *fBitDataCurrentInput;
		fBitDataCurrentPosInWord = 7;
	      }
	  }
		return true;
      }

      /** function to determine input bits below 8 bits
       * @param value    AliHLTUInt8_t &
       * @param bitCount UInt_t const &
       * @return boolean 
       */
      bool InputBits( AliHLTUInt8_t & value, UInt_t const & bitCount )
      {
	if ( bitCount>8 )
	  {
	    HLTFatal( "Internal error: Attempt to write more than 32 bits (%u)", (unsigned)bitCount );
	    return false;
	  }
	AliHLTUInt64_t temp;
	if ( !InputBits( temp, bitCount ) )
	  return false;
	value = (AliHLTUInt8_t)( temp & (AliHLTUInt64_t)0xFFFFFFFFULL );
	return true;
      }

      /** function to determine input bits between 8 and 16 bits
       * @param value    AliHLTUInt16_t &
       * @param bitCount UInt_t const &
       * @return boolean 
       */
      bool InputBits( AliHLTUInt16_t & value, UInt_t const & bitCount )
      {
	if ( bitCount>16 )
	  {
	    HLTFatal( "Internal error: Attempt to write more than 32 bits (%u)", (unsigned)bitCount );
	    return false;
	  }
	AliHLTUInt64_t temp;
	if ( !InputBits( temp, bitCount ) )
	  return false;
	value = (AliHLTUInt16_t)( temp & (AliHLTUInt64_t)0xFFFFFFFFULL );
	return true;
      }

      /** function to determine input bits between 16 and 32 bits
       * @param value    AliHLTUInt32_t &
       * @param bitCount UInt_t const &
       * @return boolean 
       */
      bool InputBits( AliHLTUInt32_t & value, UInt_t const & bitCount )
      {
	if ( bitCount>32 )
	  {
	    HLTFatal( "Internal error: Attempt to write more than 32 bits (%u)", (unsigned)bitCount );
	    return false;
	  }
	AliHLTUInt64_t temp;
	if ( !InputBits( temp, bitCount ) )
	  return false;
	value = (AliHLTUInt32_t)( temp & (AliHLTUInt64_t)0xFFFFFFFFULL );
	return true;
      }

      /** function to determine input bits between 16 and 32 bits
       * @param value    Int_t_t &
       * @param bitCount UInt_t const &
       * @return boolean 
       */
      bool InputBits( Int_t & value, UInt_t const & bitCount )
      {
	if ( bitCount>32 )
	  {
	    HLTFatal( "Internal error: Attempt to write more than 32 bits (%u)", (unsigned)bitCount );
	    return false;
	  }
	AliHLTUInt64_t temp;
	if ( !InputBits( temp, bitCount ) )
	  return false;
	value = (Int_t)( temp & (AliHLTUInt64_t)0xFFFFFFFFULL );
	return true;
      }

      /** function to determine input bits between 32 and 64 bits
       * @param value    AliHLTUInt64_t &
       * @param bitCount UInt_t const &
       * @return boolean 
       */
      bool InputBits( AliHLTUInt64_t & value, UInt_t const & bitCount )
      {
	if ( bitCount>64 )
	  {
	    HLTFatal( "Internal error: Attempt to write more than 64 bits (%u)", (unsigned)bitCount );
	    return false;
	  }
	UInt_t bitsToRead=bitCount;
	UInt_t curBitCount;
	value = 0;
	while ( bitsToRead>0 )
	  {
	    if ( fBitDataCurrentInput>=fBitDataCurrentInputEnd )
	      return false;
	    if ( bitsToRead >= fBitDataCurrentPosInWord+1 )
	      curBitCount = fBitDataCurrentPosInWord+1;
	    else
	      curBitCount = bitsToRead;
	    value = (value << curBitCount) | ( (fBitDataCurrentWord >> (fBitDataCurrentPosInWord-curBitCount+1)) & ((1 << curBitCount)-1) );
	    if ( fBitDataCurrentPosInWord < curBitCount )
	      {
		fBitDataCurrentInput++;
		if ( fBitDataCurrentInput<fBitDataCurrentInputEnd )
		  {
		    fBitDataCurrentWord = *fBitDataCurrentInput;
		    fBitDataCurrentPosInWord = 7;
		  }
	      }
	    else
	      fBitDataCurrentPosInWord -= curBitCount;
	    bitsToRead -= curBitCount;
	  }
	return true;
      }

      /** function pad 8 bits */
      void Pad8Bits()
      {
	if ( fBitDataCurrentPosInWord == 7 )
	  return;
	fBitDataCurrentInput++;
	if ( fBitDataCurrentInput<fBitDataCurrentInputEnd )
	  {
	    fBitDataCurrentWord = *fBitDataCurrentInput;
	    fBitDataCurrentPosInWord = 7;
	  }
      }

      /** function to determine input bytes
       * @param data       AliHLTUInt8_t* pointer to input data
       * @param byteCount  UInt_t const &
       * @return boolean
       */
      bool InputBytes( AliHLTUInt8_t* data, UInt_t const & byteCount )
      {
	Pad8Bits();
	if ( fBitDataCurrentInput+byteCount>fBitDataCurrentInputEnd )
	  return false;
	memcpy( data, fBitDataCurrentInput, byteCount );
	fBitDataCurrentInput += byteCount;
	if ( fBitDataCurrentInput<fBitDataCurrentInputEnd )
	  {
	    fBitDataCurrentWord = *fBitDataCurrentInput;
	    fBitDataCurrentPosInWord = 7;
	  }
	return true;
      }
      
    private:
      /** copy constructor prohibited */
      AliHLTTPCCompModelInflater(const AliHLTTPCCompModelInflater&);
      /** assignment operator prohibited */
      AliHLTTPCCompModelInflater& operator=(const AliHLTTPCCompModelInflater&);

      /** member variable for bit data current word */
      AliHLTUInt8_t fBitDataCurrentWord;
      /** member variable for bit data current position in word */
      UInt_t fBitDataCurrentPosInWord;
      /** member variable for bit data current input */
      AliHLTUInt8_t *fBitDataCurrentInput;
      /** member variable for bit data current input start */
      AliHLTUInt8_t *fBitDataCurrentInputStart;
      /** member variable for bit data current input end */
      AliHLTUInt8_t *fBitDataCurrentInputEnd;
      
      ClassDef(AliHLTTPCCompModelInflater, 0);

	    

    };
#endif
