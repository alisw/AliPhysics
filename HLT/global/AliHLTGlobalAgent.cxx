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

/** @file   AliHLTGlobalAgent.cxx
    @author Matthias Richter
    @date   
    @brief  Agent of the libAliHLTGlobal library
*/

#include <cassert>
#include "AliHLTGlobalAgent.h"
#include "AliHLTConfigurationHandler.h"
#include "TObjString.h"
#include "TObjArray.h"

// header files of library components
#include "AliHLTGlobalTrackMergerComponent.h"
#include "AliHLTGlobalEsdConverterComponent.h"
#include "AliHLTV0HistoComponent.h"

/** global instance for agent registration */
AliHLTGlobalAgent gAliHLTGlobalAgent;

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTGlobalAgent)

AliHLTGlobalAgent::AliHLTGlobalAgent()
  :
  AliHLTModuleAgent("Global")
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTGlobalAgent::~AliHLTGlobalAgent()
{
  // see header file for class documentation
}

int AliHLTGlobalAgent::RegisterComponents(AliHLTComponentHandler* pHandler) const
{
  // see header file for class documentation
  assert(pHandler);
  if (!pHandler) return -EINVAL;
  pHandler->AddComponent(new AliHLTGlobalTrackMergerComponent);
  pHandler->AddComponent(new AliHLTGlobalEsdConverterComponent);
  pHandler->AddComponent(new AliHLTV0HistoComponent );
  return 0;
}

int AliHLTGlobalAgent::CreateConfigurations(AliHLTConfigurationHandler* pHandler,
					    AliRawReader* /*rawReader*/,
					    AliRunLoader* /*runloader*/) const
{
  // see header file for class documentation
  if (!pHandler) return -EINVAL;

  /////////////////////////////////////////////////////////////////////////////////////
  //
  // assembly of the global ESD

  // define the inputs to the global ESD
  TString esdInputs="TPC-globalmerger TPC-mcTrackMarker";

  // check for the availibility
  TObjArray* pTokens=esdInputs.Tokenize(" ");
  esdInputs="";
  if (pTokens) {
    for (int n=0; n<pTokens->GetEntriesFast(); n++) {
      TString module=((TObjString*)pTokens->At(n))->GetString();
      if (pHandler->FindConfiguration(module.Data())) {
	esdInputs+=module;
	esdInputs+=" ";
      }
    }
    delete pTokens;
  }

  if (esdInputs.Length()>0) {
    HLTInfo("Configuring inputs to global HLT ESD: %s", esdInputs.Data());
  } else {
    HLTWarning("No inputs to global HLT ESD found");
  }

  pHandler->CreateConfiguration("GLOBAL-esd-converter", "GlobalEsdConverter", esdInputs.Data(), "");
  
  return 0;
}

const char* AliHLTGlobalAgent::GetReconstructionChains(AliRawReader* /*rawReader*/,
						    AliRunLoader* runloader) const
{
  // see header file for class documentation
  if (runloader) {
    // reconstruction chains for AliRoot simulation
    // Note: run loader is only available while running embedded into
    // AliRoot simulation
    return "GLOBAL-esd-converter";
  }
  return NULL;
}

const char* AliHLTGlobalAgent::GetRequiredComponentLibraries() const
{
  // see header file for class documentation

  return "libAliHLTUtil.so libAliHLTRCU.so libAliHLTTPC.so libAliHLTITS.so";
}

int AliHLTGlobalAgent::GetHandlerDescription(AliHLTComponentDataType /*dt*/,
					     AliHLTUInt32_t /*spec*/,
					     AliHLTOUTHandlerDesc& /*desc*/) const
{
  // see header file for class documentation

  return 0;
}

AliHLTOUTHandler* AliHLTGlobalAgent::GetOutputHandler(AliHLTComponentDataType /*dt*/,
						      AliHLTUInt32_t /*spec*/)
{
  // see header file for class documentation

  return NULL;
}

int AliHLTGlobalAgent::DeleteOutputHandler(AliHLTOUTHandler* pInstance)
{
  // see header file for class documentation
  if (pInstance==NULL) return -EINVAL;

  return 0;
}
