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

#include <cerrno>
#include <fstream>
#include <iostream>

#include "TString.h"

#include "AliHLTOnlineConfiguration.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTOnlineConfiguration)

AliHLTOnlineConfiguration::AliHLTOnlineConfiguration()
  : TObject()
  , fXMLBuffer()
  , fXMLSize(0)
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

int AliHLTOnlineConfiguration::LoadConfiguration(const char* filename)
{
  /// load configuration from file
  
  ifstream in;
  in.open(filename);
  if (!in.is_open())
    return -EIO;

  size_t filesize = 0;
  in.seekg(0, std::ios::end );
  filesize = in.tellg();
  in.seekg(0, std::ios::beg);

  char * content = new char[filesize];
  in.read(content, filesize);
  in.close();

  fXMLBuffer.Adopt(filesize, content);
  fXMLSize = filesize;
  SetBit(kLoaded);

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
  const UInt_t defaultSampleSize = 200;

  TObject::Print(options);
  printf("Configuration loaded: %s\n", (TestBit(kLoaded) ? "YES" : "NO"));
  TString opt = options;
  opt.ToLower();
  Bool_t full = opt.Contains("full");

  if (TestBit(kLoaded)) {
    if (full) {
      char configuration[fXMLSize + 1];
      strncpy(configuration, fXMLBuffer.GetArray(), fXMLSize);
      printf("%s\n\n", configuration);
    } else {
      Int_t sampleSize = (defaultSampleSize <= fXMLSize) ?
	defaultSampleSize : fXMLSize;
      char sample[sampleSize];
      for (int i = 0; i < sampleSize - 1; i++)
	sample[i] = fXMLBuffer.At(i);
      sample[sampleSize - 1] = '\0';
      printf("%s...\n\n", sample);
    }
  }

  printf("XML size (uncompressed): %d\n", fXMLSize);
  printf("Configuration compressed: %s\n", (TestBit(kCompressed) ? "YES" : "NO"));
  printf("Configuration parsed: %s\n", (TestBit(kParsed) ? "YES" : "NO"));
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
  fXMLBuffer.Reset();
  fXMLSize = 0;
  ResetBit(kLoaded);
  ResetBit(kCompressed);
  ResetBit(kParsed);
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
