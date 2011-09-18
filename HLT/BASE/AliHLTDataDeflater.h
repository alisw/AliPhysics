//-*- Mode: C++ -*-
// $Id$
#ifndef ALIHLTDATADEFLATER_H
#define ALIHLTDATADEFLATER_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/// @file   AliHLTDataDeflater.h
/// @author Matthias Richter, Timm Steinbeck
/// @date   2011-08-10
/// @brief  Data deflater class storing only necessary bits
/// @note   Code original from AliHLTTPCCompModelDeflater

#include "AliHLTLogging.h"
#include "AliHLTDataTypes.h"
#include "AliHLTStdIncludes.h"
#include <bitset>

/**
 * @class AliHLTDataDeflater
 * Data deflater class to write a bitstream into a buffer. The necessary
 * number of bits for each value is written to a stream without gaps.
 * The buffer can be padded to fill full bytes and continue writing at the
 * next byte.
 *
 * @ingroup alihlt_base
 */
class AliHLTDataDeflater : public AliHLTLogging
{
public:
  /// standard constructor
  AliHLTDataDeflater();
  /// destructor
  ~AliHLTDataDeflater();

  /** function to initialise bit data output
   * @param output AliHLTUInt8_t* pointer to output data
   * @param outputSize UInt_t output size
   */
  int InitBitDataOutput( AliHLTUInt8_t* output, UInt_t outputSize );

  /** function to close bit data output */
  void CloseBitDataOutput();

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

  /** function to output bits from a bitset
   * @param value     AliHLTUInt64_t const &
   * @param bitCount  UInt_t const &
   * @return zero upon success
   */
  bool OutputBits( std::bitset<64> const & value, UInt_t const & bitCount );

  /* function pad 8 bits */
  void Pad8Bits();

  /** function to output bytes
   * @param data  AliHLTUInt8_t const *
   * @param byteCount UInt_t const &
   * @return boolean (output bytes)
   */
  bool OutputBytes( AliHLTUInt8_t const * data, UInt_t const & byteCount );      

  /// clear the object and reset pointer references
  virtual void Clear(Option_t * /*option*/ ="");

  /// print info
  virtual void Print(Option_t *option="") const;

  /// print info
  virtual void Print(ostream& out, Option_t *option="") const;

  /// find object
  virtual TObject *FindObject(const char */*name*/) const {return NULL;}

  /// save data according to option
  virtual void SaveAs(const char */*filename*/="",Option_t */*option*/="") const {}

  /// write bit pattern of a parameter to the current byte and position
  virtual bool OutputParameterBits( int parameterId, AliHLTUInt64_t const & value );

  /// return unique version of the deflater, base class has version 0
  virtual int GetDeflaterVersion() const {return 0;}

protected:

private:
  /// copy constructor prohibited
  AliHLTDataDeflater(const AliHLTDataDeflater&);
  /// assignment operator prohibited
  AliHLTDataDeflater& operator=(const AliHLTDataDeflater&);

  /// bit data current word
  AliHLTUInt8_t fBitDataCurrentWord; //! bit data current word
  /// bit data current position in word
  UInt_t fBitDataCurrentPosInWord; //! data current position in word
  /// bit data current output
  AliHLTUInt8_t *fBitDataCurrentOutput; //! bit data current output
  /// bit data current output start
  AliHLTUInt8_t *fBitDataCurrentOutputStart; //! bit data current output start
  /// bit data current output end
  AliHLTUInt8_t *fBitDataCurrentOutputEnd; //! bit data current output end 

  ClassDef(AliHLTDataDeflater, 0)
};

ostream& operator<<(ostream &out, const AliHLTDataDeflater& me);

#endif
