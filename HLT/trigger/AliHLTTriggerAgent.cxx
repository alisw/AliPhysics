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
#include "AliHLTTriggerDecision.h"
#include "AliHLTOUT.h"
#include "AliHLTMessage.h"
#include "AliESDEvent.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "TArrayC.h"

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
  : AliHLTModuleAgent("Trigger")
  , fTriggerDecisionHandler(NULL)
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
    return "GLOBAL-Trigger";
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

  // handler of the trigger decisions {'ROOTTOBJ':'HLT '}
  // currently stored as a TObject with the common data type and origin
  // HLTOUT. However we might need a separate data type in order to
  // avoid interference with other handlers
  // the handler produces an ESD object in order to be merged to the
  // hltEsd afterwards
  if (dt==(kAliHLTDataTypeTObject|kAliHLTDataOriginOut)) {
    desc=AliHLTOUTHandlerDesc(AliHLTModuleAgent::kEsd, dt, GetModuleId());
    return 1;
  }

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

  // raw data blocks to be fed into offline reconstruction
  if (dt==(kAliHLTDataTypeTObject|kAliHLTDataOriginOut)) {
    if (!fTriggerDecisionHandler) {
      fTriggerDecisionHandler=new AliHLTTriggerAgent::AliHLTTriggerDecisionHandler;
    }
    return fTriggerDecisionHandler;
  }

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

  if (pInstance==fTriggerDecisionHandler) {
    delete fTriggerDecisionHandler;
    fTriggerDecisionHandler=NULL;
  }

  return 0;
}

AliHLTTriggerAgent::AliHLTTriggerDecisionHandler::AliHLTTriggerDecisionHandler()
  : AliHLTOUTHandler() 
  , fESD(NULL)
  , fpData(NULL)
  , fSize(0)
{
  // see header file for class documentation
}

AliHLTTriggerAgent::AliHLTTriggerDecisionHandler::~AliHLTTriggerDecisionHandler()
{
  // see header file for class documentation
  if (fESD) delete fESD;
  fESD=NULL;

  if (fpData) delete fpData;
  fpData=NULL;
  fSize=0;
}

int AliHLTTriggerAgent::AliHLTTriggerDecisionHandler::ProcessData(AliHLTOUT* pData)
{
  // see header file for class documentation
  if (!pData) return -EINVAL;
  pData->SelectFirstDataBlock();
  AliHLTComponentDataType dt=kAliHLTVoidDataType;
  AliHLTUInt32_t spec=kAliHLTVoidDataSpec;
  int iResult=pData->GetDataBlockDescription(dt, spec);
  if (iResult>=0) {
    TObject* pObject=pData->GetDataObject();
    if (pObject) {
      AliHLTTriggerDecision* pDecision=dynamic_cast<AliHLTTriggerDecision*>(pObject);
      if (pDecision) {
	//pDecision->Print();
	HLTDebug("extracted %s", pDecision->GetName());
	if (!fESD) {
	  // create the ESD container, but without std content
	  fESD = new AliESDEvent;
	}
	if (!fpData) fpData=new TArrayC;
	if (fESD && fpData) {
	  fESD->Reset();
	  TObject* pESDObject=fESD->FindListObject("HLTGlobalTrigger");
	  if (pESDObject) {
	    // copy the content to the already existing object
	    pObject->Copy(*pESDObject);
	  } else {
	    // add a new object
	    fESD->AddObject(pObject->Clone());
	  }
	  AliHLTMessage* pMsg=AliHLTMessage::Stream(fESD);
	  if (pMsg) {
	    if (!pMsg->CompBuffer()) {
	      fSize=pMsg->Length();
	      fpData->Set(fSize, pMsg->Buffer());
	    } else {
	      fSize=pMsg->CompLength();
	      fpData->Set(fSize, pMsg->CompBuffer());
	    }
	  } else {
	    HLTError("streaming of objects failed");
	  }
	} else {
	  HLTError("memory allocation failed");
	  iResult=-ENOMEM;
	}
      } else {
	HLTError("object %s is not an AliHLTTriggerDecision", pObject->GetName());
	iResult=-ENODATA;
      }
      pData->ReleaseDataObject(pObject);
      pObject=NULL;
    } else {
      HLTError("can not get TObject from HLTOUT buffer");
      iResult=-ENODATA;
    }
  }
  if (iResult>=0) {
    if (pData->SelectNextDataBlock()>=0) {
      HLTWarning("current implementation of trigger decision handler can only handle one block");
    }
    return fSize;
  }
  fSize=0;
  return iResult;
}

int AliHLTTriggerAgent::AliHLTTriggerDecisionHandler::GetProcessedData(const AliHLTUInt8_t* &pData)
{
  // see header file for class documentation
  if (!fpData) {
    pData=NULL;
    return 0;
  }

  pData=reinterpret_cast<AliHLTUInt8_t*>(fpData->GetArray());
  return fSize;
}

int AliHLTTriggerAgent::AliHLTTriggerDecisionHandler::ReleaseProcessedData(const AliHLTUInt8_t* pData, int size)
{
  // see header file for class documentation
  int iResult=0;
  if (!fpData || size != fSize ||
      const_cast<AliHLTUInt8_t*>(pData) != reinterpret_cast<AliHLTUInt8_t*>(fpData->GetArray())) {
    HLTError("attempt to release to wrong data buffer %p size %d, expected %p size %d", pData, size, fpData?fpData->GetArray():NULL, fSize);
  }
  fSize=0;
  return iResult;
}
