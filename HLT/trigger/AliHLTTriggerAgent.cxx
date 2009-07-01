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

/** @file   AliHLTTriggerAgent.cxx
    @author Matthias Richter
    @date   
    @brief  Agent of the libAliHLTTrigger library
*/

#include <cassert>
#include "AliHLTTriggerAgent.h"
#include "TObjString.h"
#include "TObjArray.h"

// header files of library components
#include "AliHLTEventSummaryProducerComponent.h"
#include "AliHLTRunSummaryProducerComponent.h"
#include "AliHLTTriggerBarrelMultiplicity.h"
#include "AliHLTTriggerBarrelCosmic.h"
#include "AliHLTGlobalTriggerComponent.h"

/** global instance for agent registration */
AliHLTTriggerAgent gAliHLTTriggerAgent;

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTTriggerAgent)

AliHLTTriggerAgent::AliHLTTriggerAgent()
  :
  AliHLTModuleAgent("Trigger")
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTTriggerAgent::~AliHLTTriggerAgent()
{
  // see header file for class documentation
}

int AliHLTTriggerAgent::RegisterComponents(AliHLTComponentHandler* pHandler) const
{
  // see header file for class documentation
  assert(pHandler);
  if (!pHandler) return -EINVAL;
  pHandler->AddComponent(new AliHLTGlobalTriggerComponent);
  pHandler->AddComponent(new AliHLTTriggerBarrelMultiplicity);
  pHandler->AddComponent(new AliHLTTriggerBarrelCosmic);
  return 0;
}

int AliHLTTriggerAgent::CreateConfigurations(AliHLTConfigurationHandler* pHandler,
					    AliRawReader* /*rawReader*/,
					    AliRunLoader* /*runloader*/) const
{
  // see header file for class documentation
  if (!pHandler) return -EINVAL;

  TString triggerInputs;
  TString triggerOutputs;
  TString configurationId;
  /////////////////////////////////////////////////////////////////////////////////////
  //
  // a central barrel charged particle multiplicity trigger
  configurationId="TRIGGER-Barrel-Multiplicity";

  // define the inputsfor the BarrelMultiplicityTrigger
  triggerInputs="GLOBAL-esd-converter";

  // check for the availibility
  TObjArray* pTokens=triggerInputs.Tokenize(" ");
  triggerInputs="";
  if (pTokens) {
    for (int n=0; n<pTokens->GetEntriesFast(); n++) {
      TString module=((TObjString*)pTokens->At(n))->GetString();
      if (pHandler->FindConfiguration(module.Data())) {
	triggerInputs+=module;
	triggerInputs+=" ";
      }
    }
    delete pTokens;
  }

  if (triggerInputs.Length()>0) {
    HLTInfo("Configuring inputs for %s: %s", configurationId.Data(), triggerInputs.Data());
    pHandler->CreateConfiguration(configurationId.Data(), "BarrelMultiplicityTrigger", triggerInputs.Data(), "");
    if (triggerOutputs.Length()>0) triggerOutputs+=" ";
    triggerOutputs+=configurationId;
  } else {
    HLTWarning("No inputs for %s found, skipping component", configurationId.Data());
  }

  /////////////////////////////////////////////////////////////////////////////////////
  //
  // the global trigger component
  configurationId="GLOBAL-Trigger";
  HLTInfo("setting inputs for %s: %s", configurationId.Data(), triggerOutputs.IsNull()?"none":triggerOutputs.Data());
  pHandler->CreateConfiguration(configurationId.Data(), "HLTGlobalTrigger", triggerOutputs.Data(), "");
  
  return 0;
}

const char* AliHLTTriggerAgent::GetReconstructionChains(AliRawReader* /*rawReader*/,
						    AliRunLoader* runloader) const
{
  // see header file for class documentation
  if (runloader) {
    // reconstruction chains for AliRoot simulation
    // Note: run loader is only available while running embedded into
    // AliRoot simulation

    // currently disabled due to a problem compiling the runtime trigger library
    //return "GLOBAL-Trigger";
  }
  return NULL;
}

const char* AliHLTTriggerAgent::GetRequiredComponentLibraries() const
{
  // see header file for class documentation

  return "libAliHLTUtil.so libAliHLTRCU.so libAliHLTTPC.so libAliHLTITS.so libAliHLTGlobal.so";
}

int AliHLTTriggerAgent::GetHandlerDescription(AliHLTComponentDataType dt,
					   AliHLTUInt32_t /*spec*/,
					  AliHLTOUTHandlerDesc& desc) const
{
  // see header file for class documentation

  // handler for the HLT readou list and trigger data data blocks {'HLTRDLST':'HLT '}
  if (dt==AliHLTComponentDataTypeInitializer("HLTRDLST", kAliHLTDataOriginOut) ||
      dt==AliHLTComponentDataTypeInitializer("HLTTRGDT", kAliHLTDataOriginOut)) {
      desc=AliHLTOUTHandlerDesc(kProprietary, dt, GetModuleId());
      return 1;
  }

  return 0;
}

AliHLTOUTHandler* AliHLTTriggerAgent::GetOutputHandler(AliHLTComponentDataType dt,
						   AliHLTUInt32_t /*spec*/)
{
  // see header file for class documentation

  // handler for the HLT readou list and trigger data data blocks {'HLTRDLST':'HLT '}
  if (dt==AliHLTComponentDataTypeInitializer("HLTRDLST", kAliHLTDataOriginOut) ||
      dt==AliHLTComponentDataTypeInitializer("HLTTRGDT", kAliHLTDataOriginOut)) {
    return NULL;
  }

  return NULL;
}

int AliHLTTriggerAgent::DeleteOutputHandler(AliHLTOUTHandler* pInstance)
{
  // see header file for class documentation
  if (pInstance==NULL) return -EINVAL;

  return 0;
}
