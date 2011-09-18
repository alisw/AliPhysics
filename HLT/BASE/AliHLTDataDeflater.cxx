// $Id$

//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
//*                  for The ALICE HLT Project.                            *
//*                                                                        *
//* Permission to use, copy, modify and distribute this software and its   *
//* documentation strictly for non-commercial purposes is hereby granted   *
//* without fee, provided that the above copyright notice appears in all   *
//* copies and that both the copyright notice and this permission notice   *
//* appear in the supporting documentation. The authors make no claims     *
//* about the suitability of this software for any purpose. It is          *
//* provided "as is" without express or implied warranty.                  *
//**************************************************************************

/// @file   AliHLTDataDeflater.cxx
/// @author Matthias Richter, Timm Steinbeck
/// @date   2011-08-10
/// @brief  Data deflater class storing only necessary bits
/// @note   Code original from AliHLTTPCCompModelDeflater

#include "AliHLTDataDeflater.h"
#include "AliHLTErrorGuard.h"
#include <memory>
#include <algorithm>
#include <iostream>

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTDataDeflater)

AliHLTDataDeflater::AliHLTDataDeflater()
  : AliHLTLogging()
  , fBitDataCurrentWord(0)
  , fBitDataCurrentPosInWord(0)
  , fBitDataCurrentOutput(NULL)
  , fBitDataCurrentOutputStart(NULL)
  , fBitDataCurrentOutputEnd(NULL)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTDataDeflater::~AliHLTDataDeflater()
{
  // destructor
  Clear();
}

int AliHLTDataDeflater::InitBitDataOutput( AliHLTUInt8_t* output, UInt_t outputSize)
{
  // init the target buffer
  fBitDataCurrentWord = 0;
  fBitDataCurrentPosInWord = 7;
  fBitDataCurrentOutput = fBitDataCurrentOutputStart = output;
  fBitDataCurrentOutputEnd = output+outputSize;

  return 0;
}

void AliHLTDataDeflater::CloseBitDataOutput()
{
  // pad to full byte and clear internal pointer references
  Pad8Bits();
  fBitDataCurrentWord=0;
  fBitDataCurrentPosInWord=0;
  fBitDataCurrentOutput=NULL;
  fBitDataCurrentOutputStart=NULL;
  fBitDataCurrentOutputEnd=NULL;
}

AliHLTUInt8_t AliHLTDataDeflater::GetCurrentOutputByte( Int_t offset ) const
{
  // get the current byte
  if ( !offset )
    return fBitDataCurrentWord;
  else
    return *(fBitDataCurrentOutput+offset);
}

bool AliHLTDataDeflater::OutputBit( AliHLTUInt32_t const & value )
{
  // write one bit to the current byte and position
  if ( fBitDataCurrentOutput>=fBitDataCurrentOutputEnd )
    return false;
  fBitDataCurrentWord |= (value & 1) << fBitDataCurrentPosInWord;
  if ( fBitDataCurrentPosInWord )
    fBitDataCurrentPosInWord--;
  else {
    *fBitDataCurrentOutput = fBitDataCurrentWord;
    fBitDataCurrentPosInWord = 7;
    fBitDataCurrentOutput++;
    fBitDataCurrentWord = 0;
  }
  return true;
}

bool AliHLTDataDeflater::OutputBits( AliHLTUInt64_t const & value, UInt_t const & bitCount )
{
  // write bit pattern to the current byte and position
  if ( bitCount>64 ) {
    HLTFatal( "Internal error: Attempt to write more than 64 bits (%u)", (unsigned)bitCount );
    return false;
  }
  UInt_t bitsToWrite=bitCount;
  UInt_t curBitCount;
  while ( bitsToWrite>0 ) {
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

bool AliHLTDataDeflater::OutputBits( std::bitset<64> const & value, UInt_t const & bitCount )
{
  // write bit pattern to the current byte and position
  if ( bitCount>64 ) {
    HLTFatal( "Internal error: Attempt to write more than 64 bits (%u)", (unsigned)bitCount );
    return false;
  }
  static const std::bitset<64> mask8bit(255ul);
  UInt_t bitsToWrite=bitCount;
  UInt_t curBitCount;
  while ( bitsToWrite>0 ) {
    if ( fBitDataCurrentOutput>=fBitDataCurrentOutputEnd )
      return false;
    if ( bitsToWrite >= fBitDataCurrentPosInWord+1 )
      curBitCount = fBitDataCurrentPosInWord+1;
    else
      curBitCount = bitsToWrite;
    std::bitset<64> valwrite=(value >> (bitsToWrite-curBitCount)) & mask8bit;
    fBitDataCurrentWord |= ( valwrite.to_ulong() & ((1<<curBitCount)-1) ) << (fBitDataCurrentPosInWord+1-curBitCount);
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
  }
  return true;
}

void AliHLTDataDeflater::Pad8Bits()
{
  // finish the current word
  if ( fBitDataCurrentPosInWord==7 )
    return;
  *fBitDataCurrentOutput = fBitDataCurrentWord;
  fBitDataCurrentPosInWord = 7;
  fBitDataCurrentOutput++;
  fBitDataCurrentWord = 0;
}

bool AliHLTDataDeflater::OutputBytes( AliHLTUInt8_t const * data, UInt_t const & byteCount )
{
  // write sequence of bytes
  Pad8Bits();
  if ( fBitDataCurrentOutput+byteCount>fBitDataCurrentOutputEnd )
    return false;
  memcpy( fBitDataCurrentOutput, data, byteCount );
  fBitDataCurrentOutput += byteCount;
  return true;
}

bool AliHLTDataDeflater::OutputParameterBits( int /*(parameterId*/, AliHLTUInt64_t const & /*value*/ )
{
  // write bit pattern of a member to the current byte and position
  ALIHLTERRORGUARD(1,"method needs to be implemented in child class");
  return false;
}

void AliHLTDataDeflater::Clear(Option_t * /*option*/)
{
  // internal cleanup
}

void AliHLTDataDeflater::Print(Option_t *option) const
{
  // print info
  Print(cout, option);
}

void AliHLTDataDeflater::Print(ostream& out, Option_t */*option*/) const
{
  // print to stream
  out << "AliHLTDataDeflater: " << endl;
}

ostream& operator<<(ostream &out, const AliHLTDataDeflater& me)
{
  me.Print(out);
  return out;
}
