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

#include <cstring>

#include "AliHLTRecoParamComponent.h"
#include "AliHLTReadoutList.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTRecoParamComponent)

AliHLTRecoParamComponent::AliHLTRecoParamComponent()
  : AliHLTCalibrationProcessor()
  , fOnlineConfig()
  , fOutputSize(0)
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
  const UInt_t streamerInfoEstSize = 1024; // Estimated size of streamer info
  // total size: FXSHeader + StreamerInfo + XML configuration
  constBase = AliHLTCalibrationProcessor::fgkFXSProtocolHeaderSize +
    streamerInfoEstSize + fOutputSize;
  inputMultiplier = 0;
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

  AliHLTReadoutList rdList(AliHLTReadoutList::kHLT);
  PushToFXS(&fOnlineConfig, "HLT", "OnlineRecoParam", &rdList);
  return 0;
}

int AliHLTRecoParamComponent::ScanConfigurationArgument(int argc, const char** argv)
{
  // see header file for class documentation
  int iResult=0;
  int result=0;
  char* configFile;
  if (argc == 1) {
    int argLen = strlen(argv[0]);
    char argument[argLen+1];
    strcpy(argument, argv[0]);
    argument[argLen] = '\0';
    if (strstr(argument, "-configfile")) {
      strtok(argument, "=");
      configFile = strtok(0, "=");
      if (configFile)
        result = fOnlineConfig.LoadConfiguration(configFile);
      if (result > 0) {
        fOutputSize = result; // configuration file was successfully read
	iResult = 1;
      }
    }
  }
  if (result == 0) {
    HLTError("Missing argument -configfile");
    iResult = -EPROTO;
  }
  else if (result < 0) {
    HLTError("Could not read configuration file %s", configFile);
    iResult = -ENOENT;
  }
  return iResult;
}
