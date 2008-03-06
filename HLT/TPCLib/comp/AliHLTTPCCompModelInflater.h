// XEmacs -*-C++-*-
// $Id$

#ifndef ALIHLTTPCCOMPMODELINFLATER_H
#define ALIHLTTPCCOMPMODELINFLATER_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

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
 * The model inflater is the counterpart of the deflater component
 * which decompresses the data after a conversion and compression.
 * The deflater is followed by a deconverter in order to get the
 * original standard HLT cluster track format back
 *
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
      void InitBitDataInput( AliHLTUInt8_t* input, UInt_t inputSize );

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
      bool InputBit( AliHLTUInt8_t & value );

      /** function to determine input bits below 8 bits
       * @param value    AliHLTUInt8_t &
       * @param bitCount UInt_t const &
       * @return boolean 
       */
      bool InputBits( AliHLTUInt8_t & value, UInt_t const & bitCount );

      /** function to determine input bits between 8 and 16 bits
       * @param value    AliHLTUInt16_t &
       * @param bitCount UInt_t const &
       * @return boolean 
       */
      bool InputBits( AliHLTUInt16_t & value, UInt_t const & bitCount );

      /** function to determine input bits between 16 and 32 bits
       * @param value    AliHLTUInt32_t &
       * @param bitCount UInt_t const &
       * @return boolean 
       */
      bool InputBits( AliHLTUInt32_t & value, UInt_t const & bitCount );

      /** function to determine input bits between 16 and 32 bits
       * @param value    Int_t_t &
       * @param bitCount UInt_t const &
       * @return boolean 
       */
      bool InputBits( Int_t & value, UInt_t const & bitCount );

      /** function to determine input bits between 32 and 64 bits
       * @param value    AliHLTUInt64_t &
       * @param bitCount UInt_t const &
       * @return boolean 
       */
      bool InputBits( AliHLTUInt64_t & value, UInt_t const & bitCount );
 
      /** function pad 8 bits */
      void Pad8Bits();

      /** function to determine input bytes
       * @param data       AliHLTUInt8_t* pointer to input data
       * @param byteCount  UInt_t const &
       * @return boolean
       */
      bool InputBytes( AliHLTUInt8_t* data, UInt_t const & byteCount );
      
    private:
      /** copy constructor prohibited */
      AliHLTTPCCompModelInflater(const AliHLTTPCCompModelInflater&);
      /** assignment operator prohibited */
      AliHLTTPCCompModelInflater& operator=(const AliHLTTPCCompModelInflater&);

      /** member variable for bit data current word */
      AliHLTUInt8_t fBitDataCurrentWord; // member variable for bit data current word
      /** member variable for bit data current position in word */
      UInt_t fBitDataCurrentPosInWord;// member variable for bit data current position in word
      /** member variable for bit data current input */
      AliHLTUInt8_t *fBitDataCurrentInput; // member variable for bit data current input
      /** member variable for bit data current input start */
      AliHLTUInt8_t *fBitDataCurrentInputStart; // member variable for bit data current input star
      /** member variable for bit data current input end */
      AliHLTUInt8_t *fBitDataCurrentInputEnd; // member variable for bit data current input end
      
      ClassDef(AliHLTTPCCompModelInflater, 0);

	    

    };
#endif
