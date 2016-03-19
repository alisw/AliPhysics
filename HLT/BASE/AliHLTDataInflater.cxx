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

/// @file   AliHLTDataInflater.cxx
/// @author Matthias Richter, Timm Steinbeck
/// @date   2011-08-10
/// @brief  Data inflater reading the bitstream from the AliHLTDataDeflater
/// @note   Code original from AliHLTTPCCompModelInflater

#include "AliHLTDataInflater.h"
#include "AliHLTErrorGuard.h"
#include <memory>
#include <algorithm>
#include <iostream>

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTDataInflater)

AliHLTDataInflater::AliHLTDataInflater()
  : AliHLTLogging()
  , fBitDataCurrentWord(0)
  , fBitDataCurrentPosInWord(0)
  , fBitDataCurrentInput(NULL)
  , fBitDataCurrentInputStart(NULL)
  , fBitDataCurrentInputEnd(NULL)
{
  // constructor, see header file for class documentation
}

AliHLTDataInflater::~AliHLTDataInflater()
{
  // destructor
  Clear();
}

int AliHLTDataInflater::InitBitDataInput(const AliHLTUInt8_t* input, UInt_t inputSize )
{
  // init inflater for reading
  fBitDataCurrentWord = 0;
  fBitDataCurrentPosInWord = 7;
  fBitDataCurrentInput = fBitDataCurrentInputStart = input;
  fBitDataCurrentInputEnd = input+inputSize;
  fBitDataCurrentWord = *fBitDataCurrentInput;
  return 0;
}

void AliHLTDataInflater::CloseBitDataInput()
{
  // close inflater for reading
  fBitDataCurrentWord=0;
  fBitDataCurrentPosInWord=0;
  fBitDataCurrentInput=NULL;
  fBitDataCurrentInputStart=NULL;
  fBitDataCurrentInputEnd=NULL;
}

bool AliHLTDataInflater::InputBit( AliHLTUInt8_t & value )
{
  // see header file for class documenation
  if ( fBitDataCurrentInput>=fBitDataCurrentInputEnd )
    return false;
  value = (fBitDataCurrentWord >> fBitDataCurrentPosInWord) & 1;
  if ( fBitDataCurrentPosInWord )
    fBitDataCurrentPosInWord--;
  else {
    fBitDataCurrentInput++;
    if ( fBitDataCurrentInput<fBitDataCurrentInputEnd ) {
      fBitDataCurrentWord = *fBitDataCurrentInput;
      fBitDataCurrentPosInWord = 7;
    }
  }
  HLTDebug("   code 0x%08x  length 1", value);
  return true;
}

bool AliHLTDataInflater::RewindBitPosition(UInt_t const & bitCount)
{
  // Reverse the current bit position by the given number of bits.
  UInt_t bitDataCurrentPosInWord=fBitDataCurrentPosInWord+bitCount;
  if ( bitDataCurrentPosInWord > 7) {
    UInt_t byteShift=bitDataCurrentPosInWord/8;
    if (fBitDataCurrentInputStart+byteShift>fBitDataCurrentInput) {
      return false;
    }
    fBitDataCurrentInput-=byteShift;
    fBitDataCurrentWord = *fBitDataCurrentInput;
  }
  fBitDataCurrentPosInWord = bitDataCurrentPosInWord%8;
  return true;
}

void AliHLTDataInflater::Pad8Bits()
{
  // see header file for class documenation
  if ( fBitDataCurrentPosInWord == 7 )
    return;
  fBitDataCurrentInput++;
  if ( fBitDataCurrentInput<fBitDataCurrentInputEnd ) {
    fBitDataCurrentWord = *fBitDataCurrentInput;
    fBitDataCurrentPosInWord = 7;
  }
}

bool AliHLTDataInflater::InputBytes( AliHLTUInt8_t* data, UInt_t const & byteCount )
{
  // see header file for class documenation
  Pad8Bits();
  if ( fBitDataCurrentInput+byteCount>fBitDataCurrentInputEnd )
    return false;
  memcpy( data, fBitDataCurrentInput, byteCount );
  fBitDataCurrentInput += byteCount;
  if ( fBitDataCurrentInput<fBitDataCurrentInputEnd ) {
    fBitDataCurrentWord = *fBitDataCurrentInput;
    fBitDataCurrentPosInWord = 7;
  }
  return true;
}

void AliHLTDataInflater::Clear(Option_t * /*option*/)
{
  // internal cleanup
}

void AliHLTDataInflater::Print(Option_t *option) const
{
  // print info
  Print(cout, option);
}

void AliHLTDataInflater::Print(ostream& out, Option_t */*option*/) const
{
  // print to stream
  out << "AliHLTDataInflater: " << endl;
}

ostream& operator<<(ostream &out, const AliHLTDataInflater& me)
{
  me.Print(out);
  return out;
}
