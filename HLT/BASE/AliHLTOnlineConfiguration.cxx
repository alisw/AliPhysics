// $Id$

//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors:                                                       *
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

/// @file   
/// @author 
/// @date   
/// @brief  
///

#include "AliHLTOnlineConfiguration.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTOnlineConfiguration)

AliHLTOnlineConfiguration::AliHLTOnlineConfiguration()
  : TObject()
  , fXMLBuffer()
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

}

AliHLTOnlineConfiguration::~AliHLTOnlineConfiguration()
{
  // destructor
}

int AliHLTOnlineConfiguration::LoadConfiguration(const char* /*filename*/)
{
  /// load configuration from file

  return 0;
}

int AliHLTOnlineConfiguration::Compress()
{
  /// compress the xml buffer

  return 0;
}

int AliHLTOnlineConfiguration::Uncompress()
{
  /// compress the xml buffer

  return 0;
}

void AliHLTOnlineConfiguration::Print(const char* options) const
{
  /// overloaded from TObject, print info

  TObject::Print(options);
}

void AliHLTOnlineConfiguration::Dump() const
{
  /// overloaded from TObject, more crude data dump

  TObject::Dump();
}

void AliHLTOnlineConfiguration::Clear(Option_t * option)
{
  /// overloaded from TObject, clear object

  TObject::Clear(option);
}

TObject * AliHLTOnlineConfiguration::Clone(const char *newname) const
{
  /// overloaded from TObject, clone object

  return TObject::Clone(newname);
}

void AliHLTOnlineConfiguration::Copy(TObject &object) const
{
  /// overloaded from TObject, copy object

  TObject::Copy(object);
}

void AliHLTOnlineConfiguration::Draw(Option_t */*option*/)
{
  /// overloaded from TObject, draw graph of the configuration
}
