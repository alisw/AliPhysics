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

/// @file   AliHLTDataDeflaterSimple.cxx
/// @author Matthias Richter
/// @date   2011-08-10
/// @brief  Simple deflater implementation storing frequent values below a
///         maximum value with a reduced bit number and others with the full
///         number of bits.

#include "AliHLTDataDeflaterSimple.h"
#include <memory>
#include <algorithm>
#include <iostream>

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTDataDeflaterSimple)

AliHLTDataDeflaterSimple::AliHLTDataDeflaterSimple()
  : AliHLTDataDeflater()
  , fParameterDefinitions()
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTDataDeflaterSimple::~AliHLTDataDeflaterSimple()
{
  // destructor
  Clear();
}

bool AliHLTDataDeflaterSimple::OutputParameterBits( int memberId, AliHLTUInt64_t const & value )
{
  // write bit pattern of a member to the current byte and position
  if (memberId>=(int)fParameterDefinitions.size()) return false;

  AliHLTUInt32_t switchBit=fParameterDefinitions[memberId].SwitchBit(value); // 0 -> reduced, 1 -> full
  AliHLTUInt64_t v=fParameterDefinitions[memberId].Value(value);
  AliHLTUInt32_t length=fParameterDefinitions[memberId].ValueLength(value);
  fParameterDefinitions[memberId].IncrementBitCount(value);

  if (!OutputBit(switchBit)) return false;
  return OutputBits(v, length);
}

void AliHLTDataDeflaterSimple::Clear(Option_t * /*option*/)
{
  // internal cleanup
}

void AliHLTDataDeflaterSimple::Print(Option_t *option) const
{
  // print info
  Print(cout, option);
}

void AliHLTDataDeflaterSimple::Print(ostream& out, Option_t *option) const
{
  // print to stream
  out << "AliHLTDataDeflaterSimple:" << endl;
  AliHLTUInt64_t bitCount=0;
  AliHLTUInt64_t fullSize=0;
  for (vector<AliHLTDataDeflaterParameter>::const_iterator m=fParameterDefinitions.begin();
       m!=fParameterDefinitions.end(); m++) {
    cout << "   "; m->Print(option);
    bitCount+=m->GetBitCount();
    fullSize+=m->GetValueCount()*m->GetBitLength();
  }
  out << " total: " << bitCount << "/" << fullSize << " " << (fullSize>0?float(bitCount)/fullSize:0.0) << endl;
}

void AliHLTDataDeflaterSimple::AliHLTDataDeflaterParameter::Print(const char* /*option*/) const
{
  // print info
  cout << fName << " (" << fFullBitLength << "," << fReducedBitLength << "): "
       << fValueCount << " entries  "
       << fBitCount << "/" << fFullBitLength*fValueCount;
  if (fFullBitLength && fValueCount) {
    cout << " " << float(fBitCount)/(fValueCount*fFullBitLength);
  }
  cout << endl;
}

ostream& operator<<(ostream &out, const AliHLTDataDeflaterSimple& me)
{
  me.Print(out);
  return out;
}
