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

/// @file   AliHLTObjectCompressionComponent.cxx
/// @author Matthias Richter
/// @date   2009-11-09
/// @brief  Component for compression adaption of TObjects
///

//#include <cstdlib>
//#include <cassert>
#include "AliHLTObjectCompressionComponent.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTObjectCompressionComponent)

AliHLTObjectCompressionComponent::AliHLTObjectCompressionComponent()
  : AliHLTProcessor()
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTObjectCompressionComponent::~AliHLTObjectCompressionComponent()
{
  // see header file for class documentation
}

void AliHLTObjectCompressionComponent::GetInputDataTypes(AliHLTComponentDataTypeList& list)
{
  // see header file for class documentation
  list.clear();
  list.push_back(kAliHLTAnyDataType);
}

AliHLTComponentDataType AliHLTObjectCompressionComponent::GetOutputDataType()
{
  // see header file for class documentation
  return kAliHLTAnyDataType;
}

void AliHLTObjectCompressionComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier )
{
  // see header file for class documentation
  constBase=0;
  inputMultiplier=2.0; // have to adjust an upper limit according to compression level
}

int AliHLTObjectCompressionComponent::DoInit( int argc, const char** argv )
{
  // see header file for class documentation
  int iResult=0;

  iResult=ConfigureFromArgumentString(argc, argv);

  return iResult;
}

int AliHLTObjectCompressionComponent::ScanConfigurationArgument(int argc, const char** argv)
{
  // see header file for class documentation
  if (argc<=0) return 0;
  int i=0;
  TString argument=argv[i];

  // -verbose
  if (argument.CompareTo("-verbose")==0) {
    // 
    return 1;
  }    

  return 0;
}

int AliHLTObjectCompressionComponent::DoDeinit()
{
  // see header file for class documentation
  int iResult=0;
  return iResult;
}

int AliHLTObjectCompressionComponent::DoEvent(const AliHLTComponentEventData& /*evtData*/,
				   AliHLTComponentTriggerData& /*trigData*/)
{
  // see header file for class documentation
  int iResult=0;
  for (const TObject* pObject=GetFirstInputObject();
       pObject!=NULL; 
       pObject=GetNextInputObject()) {
    PushBack(pObject, GetDataType(), GetSpecification());
  }

  return iResult;
}
