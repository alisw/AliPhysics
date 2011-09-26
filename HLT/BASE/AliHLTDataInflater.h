//-*- Mode: C++ -*-
// $Id$
#ifndef ALIHLTDATAINFLATER_H
#define ALIHLTDATAINFLATER_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/// @file   AliHLTDataInflater.h
/// @author Matthias Richter, Timm Steinbeck
/// @date   2011-08-10
/// @brief  Data inflater reading the bitstream from the AliHLTDataDeflater
/// @note   Code original from AliHLTTPCCompModelInflater

#include "AliHLTLogging.h"
#include "AliHLTDataTypes.h"
#include "AliHLTStdIncludes.h"

/**
 * @class AliHLTDataInflater
 * Data inflating interface to a bitstream encoded by AliHLTDataDeflater
 *
 * Setting up the data inflater:
 * <pre>
  AliHLTDataInflater inflater;
  if (inflater.InitBitDataInput(pData, dataSize)<0) {
    // can not initialize
  }
 * </pre>
 *
 * If the layout in the bitstream is known the parts of known bit length
 * can be read by
 * - InputBit(value)
 * - InputBits(value, length)
 *
 * For inflating huffman encoded streams where the symbol length is unknown
 * one has to read a fixed value, retrieve the symbol and rewind the stream.
 * Example how to read a fixed number of bits from the stream and rewind
 * after the length of the symbol has been determined.
 * <pre>
  AliHLTUInt64_t inputWord=0;
  int inputLength=sizeof(inputWord)*8;
  while (outClusterCnt<nofClusters && inflater.InputBits(inputWord, inputLength)) {
    // check how many of the bits belong to the symbol and rewind the
    // input stream
    int symbolLength=...;
    inputWord>>=(inputLength-symbolLength);
    if (!inflater.RewindBitPosition(inputLength-symbolLength)) {
      // some error
      break;
    }

    // do something with data in input word

    // check if there is less remaining data than the full input word bit length
    UInt_t bytes=inflater.GetRemainingBitDataSizeBytes();
    if (bytes>0 && bytes<=sizeof(inputWord)) {
      // reading the last bytes
      // available bits are determined by cuurent position+1 (because
      // position 0 means 1 bit still available) and the bits in the
      // available bytes
      inputLength=inflater.GetCurrentBitInputPosition()+1+(bytes-1)*8;
    }
  }
 * </pre>
 * 
 *
 * @ingroup alihlt_base
 */
class AliHLTDataInflater : public AliHLTLogging
{
public:
  /// standard constructor
  AliHLTDataInflater();
  /// destructor
  ~AliHLTDataInflater();

  /** init inflater for reading
   * @param input  AliHLTUInt8_t* pointer to input data
   * @param inputSize UInt_t input data size
   */
  int InitBitDataInput(const AliHLTUInt8_t* input, UInt_t inputSize );

  /** close inflater for reading */	
  void CloseBitDataInput();
      
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

  /** get number of remaining input bytes
   * the last byte might be already partially read, use
   * GetCurrentBitInputPosition()
   */
  UInt_t GetRemainingBitDataSizeBytes() const
  {
    return fBitDataCurrentInputEnd-fBitDataCurrentInput;
  }

  /** function to determine input bit
   * @return boolean (if bit is 1 or 0)
   */
  bool InputBit( AliHLTUInt8_t & value );

  /** function to read bits from bitstream
   * @param value
   * @param bitCount
   * @return boolean 
   */
  template<typename T>
  bool InputBits( T & value, UInt_t const & bitCount );

  /**
   * Rewind the current bit position by the given number of bits.
   */
  bool RewindBitPosition(UInt_t const & bitCount);
 
  /** function pad 8 bits */
  void Pad8Bits();

  /** function to determine input bytes
   * @param data       AliHLTUInt8_t* pointer to input data
   * @param byteCount  UInt_t const &
   * @return boolean
   */
  bool InputBytes( AliHLTUInt8_t* data, UInt_t const & byteCount );

  /**
   * Read the next value.
   * Data read function for inflaters for different formats
   */
  virtual bool NextValue(AliHLTUInt64_t& /*value*/, AliHLTUInt32_t& /*length*/) {return false;}

  /// clear the object and reset pointer references
  virtual void Clear(Option_t * /*option*/ ="");

  /// print info
  virtual void Print(Option_t *option="") const;

  /// print info
  virtual void Print(ostream& out, Option_t *option="") const;

protected:
private:
  /** copy constructor prohibited */
  AliHLTDataInflater(const AliHLTDataInflater&);
  /** assignment operator prohibited */
  AliHLTDataInflater& operator=(const AliHLTDataInflater&);

  /** member variable for bit data current word */
  AliHLTUInt8_t fBitDataCurrentWord; // member variable for bit data current word
  /** member variable for bit data current position in word */
  UInt_t fBitDataCurrentPosInWord;// member variable for bit data current position in word
  /** member variable for bit data current input */
  const AliHLTUInt8_t *fBitDataCurrentInput; // member variable for bit data current input
  /** member variable for bit data current input start */
  const AliHLTUInt8_t *fBitDataCurrentInputStart; // member variable for bit data current input star
  /** member variable for bit data current input end */
  const AliHLTUInt8_t *fBitDataCurrentInputEnd; // member variable for bit data current input end

  ClassDef(AliHLTDataInflater, 0)
};

template<typename T>
bool AliHLTDataInflater::InputBits( T & value, UInt_t const & bitCount )
{
  // read bits from the input stream into variable
  if (bitCount>sizeof(T)*8) {
    HLTFatal( "variable of type size %u too small to read %u bits", sizeof(T)*8, (unsigned)bitCount);
    return false;
  }
  UInt_t bitsToRead=bitCount;
  UInt_t curBitCount;
  value = 0;
  while ( bitsToRead>0 ) {
    if ( fBitDataCurrentInput>=fBitDataCurrentInputEnd )
      return false;
    if ( bitsToRead >= fBitDataCurrentPosInWord+1 )
      curBitCount = fBitDataCurrentPosInWord+1;
    else
      curBitCount = bitsToRead;
    value = (value << curBitCount) | ( (fBitDataCurrentWord >> (fBitDataCurrentPosInWord-curBitCount+1)) & ((1 << curBitCount)-1) );
    if ( fBitDataCurrentPosInWord < curBitCount ) {
      fBitDataCurrentInput++;
      fBitDataCurrentPosInWord = 7;
      if ( fBitDataCurrentInput<fBitDataCurrentInputEnd ) {
	fBitDataCurrentWord = *fBitDataCurrentInput;
      }
    }
    else
      fBitDataCurrentPosInWord -= curBitCount;
    bitsToRead -= curBitCount;
  }
  return true;
}

ostream& operator<<(ostream &out, const AliHLTDataInflater& me);

#endif
