// XEmacs -*-C++-*-
// $Id$

#ifndef ALIHLTTPCCOMPMODELDEFLATER_H
#define ALIHLTTPCCOMPMODELDEFLATER_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

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
      bool fWriteShape; // member variable to write shape 
      
      /** function to initialise bit data output
       * @param output AliHLTUInt8_t* pointer to output data
       * @param outputSize UInt_t output size
       */
      void InitBitDataOutput( AliHLTUInt8_t* output, UInt_t outputSize );
      
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
      AliHLTUInt8_t GetCurrentOutputByte( Int_t offset=0 ) const;
      
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
      bool OutputBit( AliHLTUInt32_t const & value );

      /** function to output bits 
       * @param value     AliHLTUInt64_t const &
       * @param bitCount  UInt_t const &
       * @return zero upon success
       */
      bool OutputBits( AliHLTUInt64_t const & value, UInt_t const & bitCount );

      /** function pad 8 bits */
      void Pad8Bits();

      /** function to output bytes
       * @param data  AliHLTUInt8_t const *
       * @param byteCount UInt_t const &
       * @return boolean (output bytes)
       */
      bool OutputBytes( AliHLTUInt8_t const * data, UInt_t const & byteCount );      
    private:
      /** copy constructor prohibited */
      AliHLTTPCCompModelDeflater(const AliHLTTPCCompModelDeflater&);
      /** assignment operator prohibited */
      AliHLTTPCCompModelDeflater& operator=(const AliHLTTPCCompModelDeflater&);
      
      /** member variable for bit data current word */
      AliHLTUInt8_t fBitDataCurrentWord; // member variable for bit data current word
      /** member variable for bit data current position in word */
      UInt_t fBitDataCurrentPosInWord; // member variable for bit data current position in word
      /** member variable for bit data current output */
      AliHLTUInt8_t *fBitDataCurrentOutput; // member variable for bit data current output
      /** member variable for bit data current output start */
      AliHLTUInt8_t *fBitDataCurrentOutputStart; // member variable for bit data current output start
      /** member variable for bit data current output end */
      AliHLTUInt8_t *fBitDataCurrentOutputEnd; // member variable for bit data current output end 
      
      ClassDef(AliHLTTPCCompModelDeflater, 0);	    

    };
#endif
