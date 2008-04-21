// $Id$

//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
//*                                                                        *
//* Permission to use, copy, modify and distribute this software and its   *
//* documentation strictly for non-commercial purposes is hereby granted   *
//* without fee, provided that the above copyright notice appears in all   *
//* copies and that both the copyright notice and this permission notice   *
//* appear in the supporting documentation. The authors make no claims     *
//* about the suitability of this software for any purpose. It is          *
//* provided "as is" without express or implied warranty.                  *
//**************************************************************************

/** @file   AliHLTTPCOfflineTrackerComponent.cxx
    @author Matthias Richter
    @date   
    @brief  Wrapper component to the TPC offline tracker
*/

#include "AliHLTTPCOfflineTrackerComponent.h"
#include "TString.h"
#include "TObjArray.h"
#include "TObjString.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTTPCOfflineTrackerComponent)

AliHLTTPCOfflineTrackerComponent::AliHLTTPCOfflineTrackerComponent()
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTTPCOfflineTrackerComponent::~AliHLTTPCOfflineTrackerComponent()
{
  // see header file for class documentation
}

const char* AliHLTTPCOfflineTrackerComponent::GetComponentID()
{
  // see header file for class documentation
  return "TPCOfflineTracker";
}

void AliHLTTPCOfflineTrackerComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list)
{
  // see header file for class documentation
  list.push_back(kAliHLTDataTypeAliTreeR|kAliHLTDataOriginTPC);
}

AliHLTComponentDataType AliHLTTPCOfflineTrackerComponent::GetOutputDataType()
{
  // see header file for class documentation
  return kAliHLTDataTypeESDTree|kAliHLTDataOriginTPC;
}

void AliHLTTPCOfflineTrackerComponent::GetOutputDataSize(unsigned long& constBase, double& inputMultiplier)
{
  // see header file for class documentation
  constBase = 0;inputMultiplier = 1;
}

AliHLTComponent* AliHLTTPCOfflineTrackerComponent::Spawn()
{
  // see header file for class documentation
  return new AliHLTTPCOfflineTrackerComponent;
}

int AliHLTTPCOfflineTrackerComponent::DoInit( int argc, const char** argv )
{
  // see header file for class documentation
  int iResult=0;

  TString argument="";
  TString configuration=""; 
  int bMissingParam=0;
  for (int i=0; i<argc && iResult>=0; i++) {
    argument=argv[i];
    if (argument.IsNull()) continue;

  }
  if (bMissingParam) {
    HLTError("missing parameter for argument %s", argument.Data());
    iResult=-EINVAL;
  }

  if (iResult>=0 && !configuration.IsNull()) {
    iResult=Configure(configuration.Data());
  } else {
    iResult=Reconfigure(NULL, NULL);
  }

  return iResult;
}

int AliHLTTPCOfflineTrackerComponent::DoDeinit()
{
  // see header file for class documentation
  return 0;
}

int AliHLTTPCOfflineTrackerComponent::DoEvent( const AliHLTComponentEventData& /*evtData*/, AliHLTComponentTriggerData& /*trigData*/)
{
  // see header file for class documentation
  HLTInfo("processing data");

  return 0;
}

int AliHLTTPCOfflineTrackerComponent::Configure(const char* arguments)
{
  // see header file for class documentation
  int iResult=0;
  if (!arguments) return iResult;

  TString allArgs=arguments;
  TString argument;
  int bMissingParam=0;

  TObjArray* pTokens=allArgs.Tokenize(" ");
  if (pTokens) {
    for (int i=0; i<pTokens->GetEntries() && iResult>=0; i++) {
      argument=((TObjString*)pTokens->At(i))->GetString();
      if (argument.IsNull()) continue;

      if (argument.CompareTo("-something")==0) {
	if ((bMissingParam=(++i>=pTokens->GetEntries()))) break;

      } else {
	HLTError("unknown argument %s", argument.Data());
	iResult=-EINVAL;
	break;
      }
    }
    delete pTokens;
  }
  if (bMissingParam) {
    HLTError("missing parameter for argument %s", argument.Data());
    iResult=-EINVAL;
  }
  return iResult;
}

int AliHLTTPCOfflineTrackerComponent::Reconfigure(const char* /*cdbEntry*/, const char* /*chainId*/)
{
  // see header file for class documentation
  int iResult=0;
  return iResult;
}
