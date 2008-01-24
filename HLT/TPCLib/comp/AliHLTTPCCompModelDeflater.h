// XEmacs -*-C++-*-
// $Id: AliHLTTPCCompModelDeflater.h,v 1.2 2006/08/10 09:46:51 richterm Exp $

#ifndef ALIHLTTPCCOMPMODELDEFLATER_H
#define ALIHLTTPCCOMPMODELDEFLATER_H
/* TPCCompModelDeflaterright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full TPCCompModelDeflaterright notice                               */

/** @file   AliHLTTPCCompModelDeflater.h
    @author Timm Steinbeck
    @date   
    @brief  Declaration of a compressor class for TPC track data. */

#include "AliHLTTPCTrackArray.h"
#include "AliHLTTPCTrackletDataFormat.h"
#include "AliHLTTPCClusterDataFormat.h"
#include "AliHLTLogging.h"
#include "AliHLTTPCModels.h"

/**
 * @class AliHLTTPCCompModelDeflater
 * @brief A HLT TPC track compressor. 
 *
 * A class that can compress HLT TPC track model data
 */
class AliHLTTPCCompModelDeflater: public AliHLTLogging
    {
    public:
      /** standard constructor */
      AliHLTTPCCompModelDeflater();
      /** standard destructor */
      virtual ~AliHLTTPCCompModelDeflater();

      /** function to write shape */
      void WriteShape( bool write=true )
      {
	fWriteShape = write;
      }

      /** function to compress tracks
       * @param inData     AliHLTUInt8_t* pointer to input data
       * @param inputSize  UInt_t const& input size
       * @param output     AliHLTUInt8_t* pointer to output data
       * @param outputSize UInt_t& output size
       * @return zero upon success
      */
	int CompressTracks( AliHLTUInt8_t* inData, UInt_t const& inputSize, AliHLTUInt8_t* output, UInt_t& outputSize );

      /** function to compress remaining clusters
       * @param inData     AliHLTUInt8_t* pointer to input data
       * @param inputSize  UInt_t const& input size
       * @param output     AliHLTUInt8_t* pointer to output data
       * @param outputSize UInt_t& output size
       * @return zero upon success
       */
	int CompressRemainingClusters( AliHLTUInt8_t* inData, UInt_t const& inputSize, AliHLTUInt8_t* output, UInt_t& outputSize );

    protected:

      /** member variable to write shape */
      bool fWriteShape;
      
      /** function to initialise bit data output
       * @param output AliHLTUInt8_t* pointer to output data
       * @param outputSize UInt_t output size
       */
      void InitBitDataOutput( AliHLTUInt8_t* output, UInt_t outputSize )
      {
	fBitDataCurrentWord = 0;
	fBitDataCurrentPosInWord = 7;
	fBitDataCurrentOutput = fBitDataCurrentOutputStart = output;
	fBitDataCurrentOutputEnd = output+outputSize;
      }
      
      /** function to close bit data output */
      void CloseBitDataOutput()
      {
	Pad8Bits();
      }
      
      /** function to get current byte output position
       * @return unsigned long value for current byte output position
       */
      unsigned long GetCurrentByteOutputPosition() const
      {
	return (unsigned long)( fBitDataCurrentOutput - fBitDataCurrentOutputStart );
      }

      /** function to get current bit output position
       * @return unsigned long value for current bit output position
       */
      unsigned GetCurrentBitOutputPosition() const
      {
	return fBitDataCurrentPosInWord;
      }

      /** function to get current output byte
       * @param offset Int_t (set to zero if not specified explicitly)
       * @return AliHLTUInt8_t value for current output byte
       */
      AliHLTUInt8_t GetCurrentOutputByte( Int_t offset=0 ) const
      {
	if ( !offset )
	  return fBitDataCurrentWord;
	else
	  return *(fBitDataCurrentOutput+offset);
      }
      
      /** function to get bit data output size bytes
       * @return UInt_t value of bit data output size bytes
       */
      UInt_t GetBitDataOutputSizeBytes() const
      {
	return fBitDataCurrentOutput-fBitDataCurrentOutputStart;
      }
      
      /** function for output bit
       * @param value  AliHLTUInt32_t const & input
       * @return boolean (output bit)
       */
      bool OutputBit( AliHLTUInt32_t const & value )
      {
	if ( fBitDataCurrentOutput>=fBitDataCurrentOutputEnd )
	  return false;
	fBitDataCurrentWord |= (value & 1) << fBitDataCurrentPosInWord;
	if ( fBitDataCurrentPosInWord )
	  fBitDataCurrentPosInWord--;
	else
	  {
	    *fBitDataCurrentOutput = fBitDataCurrentWord;
	    fBitDataCurrentPosInWord = 7;
	    fBitDataCurrentOutput++;
	    fBitDataCurrentWord = 0;
	  }
	return true;
      }

      /** function to output bits 
       * @param value     AliHLTUInt64_t const &
       * @param bitCount  UInt_t const &
       * @return zero upon success
       */
      bool OutputBits( AliHLTUInt64_t const & value, UInt_t const & bitCount )
      {
	if ( bitCount>64 )
	  {
	    HLTFatal( "Internal error: Attempt to write more than 64 bits (%u)", (unsigned)bitCount );
	    return false;
	  }
	UInt_t bitsToWrite=bitCount;
	UInt_t curBitCount;
	while ( bitsToWrite>0 )
	  {
	    if ( fBitDataCurrentOutput>=fBitDataCurrentOutputEnd )
	      return false;
#if 1
	    if ( bitsToWrite >= fBitDataCurrentPosInWord+1 )
	      curBitCount = fBitDataCurrentPosInWord+1;
	    else
	      curBitCount = bitsToWrite;
	    fBitDataCurrentWord |= ( (value >> (bitsToWrite-curBitCount)) & ((1<<curBitCount)-1) ) << (fBitDataCurrentPosInWord+1-curBitCount);
	    if ( fBitDataCurrentPosInWord < curBitCount )
	      {
		*fBitDataCurrentOutput = fBitDataCurrentWord;
		fBitDataCurrentPosInWord = 7;
		fBitDataCurrentOutput++;
		fBitDataCurrentWord = 0;
	      }
	    else
	      fBitDataCurrentPosInWord -= curBitCount;
	    bitsToWrite -= curBitCount;
	    
#else
	    AliHLTUInt8_t curValue;
	    if ( bitsToWrite>=8 )
	      {
		curBitCount=8;
		curValue = (value >> bitsToWrite-8) & 0xFF;
		bitsToWrite -= 8;
	      }
	    else
	      {
		curBitCount=bitsToWrite;
		curValue = value & ( (1<<bitsToWrite)-1 );
		bitsToWrite = 0;
	      }
	    if ( fBitDataCurrentPosInWord+1>curBitCount )
	      {
		fBitDataCurrentWord |= curValue << (fBitDataCurrentPosInWord-curBitCount+1);
		fBitDataCurrentPosInWord -= curBitCount;
	      }
	    else if ( fBitDataCurrentPosInWord+1==curBitCount )
	      {
		fBitDataCurrentWord |= curValue;
		*fBitDataCurrentOutput = fBitDataCurrentWord;
		fBitDataCurrentPosInWord = 7;
		fBitDataCurrentOutput++;
		fBitDataCurrentWord = 0;
	      }
	    else
	      {
		const UInt_t first = fBitDataCurrentPosInWord+1; // Number of bits for first block
		const UInt_t second = curBitCount-first; // Number of bits for second block
		fBitDataCurrentWord |= ( curValue >> second ) & ((1<<first)-1);
		*fBitDataCurrentOutput = fBitDataCurrentWord;
		fBitDataCurrentOutput++;
		if ( fBitDataCurrentOutput>=fBitDataCurrentOutputEnd )
		  return false;
		fBitDataCurrentWord = curValue & ((1<<second)-1) << (8-second);
		fBitDataCurrentPosInWord = 7-second;
	      }
#endif
	  }
	return true;
      }

      /** function pad 8 bits */
      void Pad8Bits()
      {
	if ( fBitDataCurrentPosInWord==7 )
	  return;
	*fBitDataCurrentOutput = fBitDataCurrentWord;
	fBitDataCurrentPosInWord = 7;
	fBitDataCurrentOutput++;
	fBitDataCurrentWord = 0;
      }

      /** function to output bytes
       * @param data  AliHLTUInt8_t const *
       * @param byteCount UInt_t const &
       * @return boolean (output bytes)
       */
      bool OutputBytes( AliHLTUInt8_t const * data, UInt_t const & byteCount )
      {
	Pad8Bits();
	if ( fBitDataCurrentOutput+byteCount>fBitDataCurrentOutputEnd )
	  return false;
	memcpy( fBitDataCurrentOutput, data, byteCount );
	fBitDataCurrentOutput += byteCount;
	return true;
      }
      
    private:
      /** copy constructor prohibited */
      AliHLTTPCCompModelDeflater(const AliHLTTPCCompModelDeflater&);
      /** assignment operator prohibited */
      AliHLTTPCCompModelDeflater& operator=(const AliHLTTPCCompModelDeflater&);
      
      /** member variable for bit data current word */
      AliHLTUInt8_t fBitDataCurrentWord;
      /** member variable for bit data current position in word */
      UInt_t fBitDataCurrentPosInWord;
      /** member variable for bit data current output */
      AliHLTUInt8_t *fBitDataCurrentOutput;
      /** member variable for bit data current output start */
      AliHLTUInt8_t *fBitDataCurrentOutputStart;
      /** member variable for bit data current output end */
      AliHLTUInt8_t *fBitDataCurrentOutputEnd;
      
      ClassDef(AliHLTTPCCompModelDeflater, 0);	    

    };
#endif
