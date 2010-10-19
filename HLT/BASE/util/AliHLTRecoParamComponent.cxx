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

/// @file   AliHLTRecoParamComponent.cxx
/// @author Matthias Richter
/// @date   2010-10-18
/// @brief  Online HLT RecoParam generator component
///

#include "AliHLTRecoParamComponent.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTRecoParamComponent)

AliHLTRecoParamComponent::AliHLTRecoParamComponent()
  : AliHLTCalibrationProcessor()
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

}

AliHLTRecoParamComponent::~AliHLTRecoParamComponent()
{
  // see header file for class documentation
}

void AliHLTRecoParamComponent::GetInputDataTypes(AliHLTComponentDataTypeList& list)
{
  // see header file for class documentation
  list.push_back(kAliHLTAnyDataType);
}

AliHLTComponentDataType AliHLTRecoParamComponent::GetOutputDataType()
{
  // see header file for class documentation
  return kAliHLTDataTypeFXSCalib;
}

void AliHLTRecoParamComponent::GetOutputDataSize(unsigned long& constBase, double& inputMultiplier)
{
  // see header file for class documentation

  // this is nothing more than an assumption, in fact it's very difficult to predict how
  // much output the component produces
  constBase=100*1024;
  inputMultiplier=1;
}

void AliHLTRecoParamComponent::GetOCDBObjectDescription( TMap* const /*targetArray*/)
{
  // see header file for class documentation
}

int AliHLTRecoParamComponent::InitCalibration()
{
  // see header file for class documentation

  int iResult=0;

  return iResult;
}

int AliHLTRecoParamComponent::DeinitCalibration()
{
  // see header file for class documentation

  int iResult=0;

  return iResult;
}

int AliHLTRecoParamComponent::ProcessCalibration( const AliHLTComponentEventData& /*evtData*/,
							    AliHLTComponentTriggerData& /*trigData*/ )
{
  // see header file for class documentation
  int iResult=0;

  return iResult;
}

int AliHLTRecoParamComponent::ShipDataToFXS( const AliHLTComponentEventData& /*evtData*/,
						       AliHLTComponentTriggerData& /*trigData*/)
{
  // see header file for class documentation

  return 0;
}

int AliHLTRecoParamComponent::ScanConfigurationArgument(int argc, const char** argv)
{
  // see header file for class documentation
  int iResult=0;
  if (argc<=0) return 0;
  int i=0;
  TString argument=argv[i];

  // -firstarg
  if (argument.Contains("-firstarg")) {
    return 1;
  }

  return iResult;
}
