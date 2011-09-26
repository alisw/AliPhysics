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

/// @file   AliHLTDataInflaterSimple.cxx
/// @author Matthias Richter
/// @date   2011-09-01
/// @brief  Data inflater implementation for format of AliHLTDataDeflaterSimple
/// @note   

#include "AliHLTDataInflaterSimple.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTDataInflaterSimple)

AliHLTDataInflaterSimple::AliHLTDataInflaterSimple()
  : AliHLTDataInflater()
  , fParameterDefinitions()
  , fCurrentParameter(-1)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTDataInflaterSimple::~AliHLTDataInflaterSimple()
{
  // destructor
}

int AliHLTDataInflaterSimple::AddParameterDefinition(const char* name, int bitLength, int reducedBitLength)
{
  /// add a parameter definition to the configuration, return reference id
  fParameterDefinitions.push_back(AliHLTDataDeflaterSimple::AliHLTDataDeflaterParameter(name, bitLength, reducedBitLength));
  return fParameterDefinitions.size()-1;
}

bool AliHLTDataInflaterSimple::NextValue(AliHLTUInt64_t& value, AliHLTUInt32_t& length)
{
  /// overloaded from AliHLTDataInflater
  /// functions reads the sequence of parameters as defined by the decoder
  /// list, than it starts at the first parameter again
  value=0;
  length=0;
  if (fParameterDefinitions.size()==0) return false;
  if ((++fCurrentParameter)>=(int)fParameterDefinitions.size()) fCurrentParameter=0;
  const AliHLTDataDeflaterSimple::AliHLTDataDeflaterParameter& parameter
    =fParameterDefinitions[fCurrentParameter];

  AliHLTUInt8_t switchBit=0;
  if (!InputBit(switchBit))
    return false;
  int readlength=switchBit?parameter.GetBitLength():parameter.GetReducedBitLength();
  if (!InputBits(value, readlength))
    return false;
  length=parameter.GetBitLength();

  return true;
}
