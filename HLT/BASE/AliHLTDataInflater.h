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
 * Data inflating interface to a bitstream encoded by AliHLTDataInflater
 * 
 * 
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
 
  /** function pad 8 bits */
  void Pad8Bits();

  /** function to determine input bytes
   * @param data       AliHLTUInt8_t* pointer to input data
   * @param byteCount  UInt_t const &
   * @return boolean
   */
  bool InputBytes( AliHLTUInt8_t* data, UInt_t const & byteCount );

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
      if ( fBitDataCurrentInput<fBitDataCurrentInputEnd ) {
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

ostream& operator<<(ostream &out, const AliHLTDataInflater& me);

#endif
