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
#include "AliHLTGlobalTriggerDecision.h"
#include "AliHLTOUT.h"
#include "AliHLTMessage.h"
#include "AliESDEvent.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "TArrayC.h"
#include "TFile.h"
#include "TTree.h"

// header files of library components
#include "AliHLTTriggerBarrelMultiplicity.h"
#include "AliHLTD0Trigger.h"
#include "AliHLTTriggerITSMultiplicity.h"
#include "AliHLTTriggerBarrelGeomMultiplicity.h"
#include "AliHLTGlobalTriggerComponent.h"
#include "AliHLTTriggerPhosClusterEnergy.h"
#include "AliHLTTriggerEmcalClusterEnergy.h"
#include "AliHLTTriggerPhosMip.h"
#include "AliHLTTriggerTrdClusterMultiplicity.h"
#include "AliHLTTriggerGammaConversion.h"
#include "AliHLTMuonSpectroTriggerComponent.h"
#include "AliHLTUpcTriggerComponent.h"
#include "AliHLTTriggerCosmics.h"
#include "AliHLTTriggerCounterComponent.h"
#include "AliHLTTriggerEmcalElectron.h"
//#include "AliHLTTriggerFastJet.h"
#include "AliHLTFastJetMonitorComponent.h"
#include "AliHLTEmcalElectronMonitorComponent.h"

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
  pHandler->AddComponent(new AliHLTTriggerITSMultiplicity);
  pHandler->AddComponent(new AliHLTD0Trigger);
  pHandler->AddComponent(new AliHLTTriggerBarrelGeomMultiplicity);
  pHandler->AddComponent(new AliHLTTriggerPhosClusterEnergy); 
  pHandler->AddComponent(new AliHLTTriggerEmcalClusterEnergy); 
  pHandler->AddComponent(new AliHLTTriggerPhosMip); 
  pHandler->AddComponent(new AliHLTTriggerTrdClusterMultiplicity);
  pHandler->AddComponent(new AliHLTTriggerGammaConversion);
  pHandler->AddComponent(new AliHLTMuonSpectroTriggerComponent);
  pHandler->AddComponent(new AliHLTUpcTriggerComponent);
  pHandler->AddComponent(new AliHLTTriggerCosmics);
  pHandler->AddComponent(new AliHLTTriggerCounterComponent);
  pHandler->AddComponent(new AliHLTTriggerEmcalElectron);
  //pHandler->AddComponent(new AliHLTTriggerFastJet);
  pHandler->AddComponent(new AliHLTFastJetMonitorComponent);
  pHandler->AddComponent(new AliHLTEmcalElectronMonitorComponent);
 return 0;
}

int AliHLTTriggerAgent::CreateConfigurations(AliHLTConfigurationHandler* pHandler,
					    AliRawReader* rawReader,
					    AliRunLoader* runloader) const
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

  // define the inputs for the BarrelMultiplicityTrigger
  triggerInputs="GLOBAL-esd-converter";

  // check for the availibility
  TObjArray* pTokens=triggerInputs.Tokenize(" ");
  triggerInputs="";
  if (pTokens) {
    for (int n=0; n<pTokens->GetEntriesFast(); n++) {
      if (!pTokens->At(n)) continue;
      TString module=pTokens->At(n)->GetName();
      if (pHandler->FindConfiguration(module.Data())) {
	triggerInputs+=module;
	triggerInputs+=" ";
      }
    }
    delete pTokens;
  }

  TString arg;
  if (triggerInputs.Length()>0) {
    // define multiple instances of the BarrelMultiplicityTrigger with different settings
    HLTInfo("Configuring inputs for %s: %s", configurationId.Data(), triggerInputs.Data());
    pHandler->CreateConfiguration(configurationId.Data(), "BarrelMultiplicityTrigger", triggerInputs.Data(), "");
    if (triggerOutputs.Length()>0) triggerOutputs+=" ";
    triggerOutputs+=configurationId;

    configurationId="TRIGGER-Barrel-HighMultiplicity";
    arg="-triggername BarrelHighMultiplicity";
    pHandler->CreateConfiguration(configurationId.Data(), "BarrelMultiplicityTrigger", triggerInputs.Data(), arg.Data());
    if (triggerOutputs.Length()>0) triggerOutputs+=" ";
    triggerOutputs+=configurationId;

    configurationId="TRIGGER-Barrel-Pt_v01";
    arg="-triggername BarrelPt_v01";
    pHandler->CreateConfiguration(configurationId.Data(), "BarrelMultiplicityTrigger", triggerInputs.Data(), arg.Data());
    if (triggerOutputs.Length()>0) triggerOutputs+=" ";
    triggerOutputs+=configurationId;

    configurationId="TRIGGER-Barrel-Pt_v02";
    arg="-triggername BarrelPt_v02";
    pHandler->CreateConfiguration(configurationId.Data(), "BarrelMultiplicityTrigger", triggerInputs.Data(), arg.Data());
    if (triggerOutputs.Length()>0) triggerOutputs+=" ";
    triggerOutputs+=configurationId;

    configurationId="TRIGGER-Barrel-Pt_v03";
    arg="-triggername BarrelPt_v03";
    pHandler->CreateConfiguration(configurationId.Data(), "BarrelMultiplicityTrigger", triggerInputs.Data(), arg.Data());
    if (triggerOutputs.Length()>0) triggerOutputs+=" ";
    triggerOutputs+=configurationId;
  } else {
    HLTWarning("No inputs for %s found, skipping component", configurationId.Data());
  }
  
  /////////////////////////////////////////////////////////////////////////////////////
  // The muon spectrometer trigger
  configurationId = "TRIGGER-Muon-Spectrometer";

  // define the inputsfor the muon spectrometer trigger.
  if (pHandler->FindConfiguration("dHLT-sim-fromRaw")) {
    triggerInputs = "dHLT-sim-fromRaw";
  }
  else if (pHandler->FindConfiguration("dHLT-sim")) {
    triggerInputs = "dHLT-sim";
  }
  else if (pHandler->FindConfiguration("dHLT-sim-fromMC")) {
    triggerInputs = "dHLT-sim-fromMC";
  }

  if (triggerInputs.Length() > 0) {
    HLTInfo("Configuring inputs for %s: %s", configurationId.Data(), triggerInputs.Data());
    pHandler->CreateConfiguration(configurationId.Data(), "MuonSpectroTrigger", triggerInputs.Data(), "-makestats");
    if (triggerOutputs.Length() > 0) triggerOutputs += " ";
    triggerOutputs += configurationId;
  } else {
    HLTWarning("No inputs for %s found, skipping component.", configurationId.Data());
  }

  /////////////////////////////////////////////////////////////////////////////////////
  // D0 trigger
  configurationId = "TRIGGER-D0";
  if(runloader && !rawReader){
    // simulation without simulated raw data
    // use ESD as input
    triggerInputs="GLOBAL-esd-converter ";
  }
  else{
    // simulation with simulated raw data, or raw data reconstruction
    // use input from ITS tracker and vertexer directly
    triggerInputs="ITS-tracker GLOBAL-vertexer";
  }

  // check for the availibility of inputs
  pTokens=triggerInputs.Tokenize(" ");
  triggerInputs="";
  if (pTokens) {
    for (int n=0; n<pTokens->GetEntriesFast(); n++) {
      if (!pTokens->At(n)) continue;
      TString module=pTokens->At(n)->GetName();
      if (pHandler->FindConfiguration(module.Data())) {
	triggerInputs+=module;
	triggerInputs+=" ";
      }
    }
    delete pTokens;
  }

  TString argD0 = "";
  if (triggerInputs.Length()>0) {
    HLTInfo("Configuring inputs for %s: %s", configurationId.Data(), triggerInputs.Data());
    pHandler->CreateConfiguration(configurationId.Data(), "D0Trigger", triggerInputs.Data(), argD0.Data());
    // FIXME: due to a rare segfault for reconstruction of PbPb data the
    // component is temporarily excluded
    // https://savannah.cern.ch/bugs/?72590
    //if (triggerOutputs.Length()>0) triggerOutputs+=" ";
    //triggerOutputs+=configurationId;
  } else {
    HLTWarning("No inputs for %s found, skipping component", configurationId.Data());
  }

  /////////////////////////////////////////////////////////////////////////////////////
  //
  // the global trigger component
  configurationId="GLOBAL-Trigger";
  HLTInfo("setting inputs for %s: %s", configurationId.Data(), triggerOutputs.IsNull()?"none":triggerOutputs.Data());
  pHandler->CreateConfiguration(configurationId.Data(), "HLTGlobalTrigger", triggerOutputs.Data(), "-skipctp");
  
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

  return "";
}

int AliHLTTriggerAgent::GetHandlerDescription(AliHLTComponentDataType dt,
					   AliHLTUInt32_t /*spec*/,
					  AliHLTOUTHandlerDesc& desc) const
{
  // get description of HLTOUT data blocks handled by this agent

  // handler of the trigger decisions {'ROOTTOBJ':'HLT '}
  // currently stored as a TObject with the common data type and origin
  // HLTOUT. However we might need a separate data type in order to
  // avoid interference with other handlers
  // the handler produces an ESD object in order to be merged to the
  // hltEsd afterwards
  // 2009-11-17 adding the data types for (global) trigger decisions
  // {GLOBTRIG:HLT } and {TRIG_DEC:HLT }
  // the TObject data types stays for a while in order to preserve
  // backward compatibility
  if (dt==(kAliHLTDataTypeTObject|kAliHLTDataOriginOut) ||
      dt==kAliHLTDataTypeTriggerDecision ||
      dt==kAliHLTDataTypeGlobalTrigger) {
    desc=AliHLTOUTHandlerDesc(AliHLTModuleAgent::kEsd, dt, GetModuleId());
    return 1;
  }

  // handler for the HLT readout list and trigger data data blocks {'HLTRDLST':'HLT '}
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
  if ((dt==(kAliHLTDataTypeTObject|kAliHLTDataOriginOut) ||
       (dt==kAliHLTDataTypeTriggerDecision) ||
       (dt==kAliHLTDataTypeGlobalTrigger))) {
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
  , fpESDfile(NULL)
  , fpESDtree(NULL)
{
  // see header file for class documentation
}

AliHLTTriggerAgent::AliHLTTriggerDecisionHandler::~AliHLTTriggerDecisionHandler()
{
  // see header file for class documentation
  if (fpESDtree) {
    fpESDtree->GetUserInfo()->Clear();
    delete fpESDtree;
  }
  fpESDtree=NULL;

  if (fpESDfile) {
    fpESDfile->Close();
    delete fpESDfile;
  }
  fpESDfile=NULL;

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
  int iResult=0;
  AliHLTGlobalTriggerDecision* pGlobalDecision=NULL;
  TObjArray triggerDecisions;
  triggerDecisions.SetOwner(kTRUE);
  for (iResult=pData->SelectFirstDataBlock(); iResult>=0; iResult=pData->SelectNextDataBlock()) {
    AliHLTComponentDataType dt=kAliHLTVoidDataType;
    AliHLTUInt32_t spec=kAliHLTVoidDataSpec;
    if ((iResult=pData->GetDataBlockDescription(dt, spec))<0) break;
    TObject* pObject=pData->GetDataObject();
    if (pObject) {
      if(dt==kAliHLTDataTypeGlobalTrigger) {
	if (!pGlobalDecision) {
	  if ((pGlobalDecision=dynamic_cast<AliHLTGlobalTriggerDecision*>(pObject))==NULL ||
	      (pGlobalDecision=dynamic_cast<AliHLTGlobalTriggerDecision*>(pGlobalDecision->Clone()))==NULL) {
	    HLTFatal("can not convert object of name %s (%s) to HLTGlobalTriggerDecsion according to data type", pObject->GetName(), pObject->Class()->GetName());
	  }
	} else {
	  HLTWarning("multiple HLT GlobalTrigger decision objects, ignoring all but the first one");
	}
      } else if (dt==kAliHLTDataTypeTriggerDecision) {
	if (pObject->IsA() == AliHLTTriggerDecision::Class() &&
	    !(pObject->IsA() == AliHLTGlobalTriggerDecision::Class())) {
	  AliHLTTriggerDecision* pDecision=dynamic_cast<AliHLTTriggerDecision*>(pObject);
	  if (pDecision) {
	    if (pGlobalDecision) {
	      // add directly
	      pGlobalDecision->AddTriggerInput(*pDecision);
	    } else {
	      // schedule
	      triggerDecisions.Add(pDecision->Clone());
	    }
	  } else {
	    HLTFatal("can not convert object of name %s (%s) to HLT TriggerDecsion according to data type", pObject->GetName(), pObject->Class()->GetName());
	  }
	}
      } else if (dt==(kAliHLTDataTypeTObject|kAliHLTDataOriginOut)){
	// this is the branch for keeping compatibility
	// the first version of the trigger framework was using the kAliHLTDataTypeTObject
	// data type instead of the specific data types for HLT triggers
	  // this effects the cosmic data taken Sep to Oct 2009
	if (pObject->IsA() == AliHLTGlobalTriggerDecision::Class()) {
	  if (!pGlobalDecision) {
	    if ((pGlobalDecision=dynamic_cast<AliHLTGlobalTriggerDecision*>(pObject))==NULL ||
		(pGlobalDecision=dynamic_cast<AliHLTGlobalTriggerDecision*>(pGlobalDecision->Clone()))==NULL) {
	      HLTFatal("can not convert object of name %s (%s) to HLTGlobalTriggerDecsion according to data type", pObject->GetName(), pObject->Class()->GetName());
	    }
	  } else {
	    HLTWarning("multiple HLT GlobalTrigger decision objects, ignoring all but the first one");
	  }
	} else if (pObject->IsA() == AliHLTTriggerDecision::Class()) {
	  AliHLTTriggerDecision* pDecision=dynamic_cast<AliHLTTriggerDecision*>(pObject);
	  if (pDecision) {
	    if (pGlobalDecision) {
	      // add directly
	      pGlobalDecision->AddTriggerInput(*pDecision);
	    } else {
	      // schedule
	      triggerDecisions.Add(pDecision->Clone());
	    }
	  } else {
	    HLTFatal("can not convert object of name %s (%s) to HLT TriggerDecsion according to data type", pObject->GetName(), pObject->Class()->GetName());
	  }
	}
      }
      pData->ReleaseDataObject(pObject);
      pObject=NULL;
    } else {
      HLTError("can not get TObject from HLTOUT buffer");
      iResult=-ENODATA;
    }
  }
  // -ENOENT just signals that there  are no more entries
  if (iResult==-ENOENT) iResult=0;

  if (pGlobalDecision) {
    for (int i=0; i<triggerDecisions.GetEntriesFast(); i++) {
      if (triggerDecisions[i]) {
	pGlobalDecision->AddTriggerInput(*((AliHLTTriggerDecision*)triggerDecisions[i]));
      }
    }
    triggerDecisions.Delete();
    AliHLTTriggerDecision* pDecision=pGlobalDecision;
    {
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
	    pDecision->Copy(*pESDObject);
	  } else {
	    // add a new object
	    fESD->AddObject(pDecision->Clone());
	  }
	  WriteESD();
	  AliHLTMessage* pMsg=AliHLTMessage::Stream(fESD);
	  if (pMsg) {
	    if (!pMsg->CompBuffer()) {
	      fSize=pMsg->Length();
	      fpData->Set(fSize, pMsg->Buffer());
	    } else {
	      fSize=pMsg->CompLength();
	      fpData->Set(fSize, pMsg->CompBuffer());
	    }
	    delete pMsg;
	    pMsg=NULL;
	  } else {
	    HLTError("streaming of objects failed");
	  }
	} else {
	  HLTError("memory allocation failed");
	  iResult=-ENOMEM;
	}
      }
    }
    delete pGlobalDecision;
    pGlobalDecision=NULL;
  } else {
    HLTError("no global trigger found in data collection");
  }

  if (iResult>=0) {
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

int AliHLTTriggerAgent::AliHLTTriggerDecisionHandler::WriteESD()
{
  // see header file for class documentation
  int iResult=0;
  if (!fESD) return 0;
  if (!fpESDfile) {
    fpESDfile=new TFile("HLTdecision.root", "RECREATE");
  }
  if (!fpESDtree) {
    fpESDtree=new TTree("HLTesdTree", "Tree with HLT ESD containing HLT decision");
    if (fpESDtree) {
      fESD->WriteToTree(fpESDtree);
      fpESDtree->GetUserInfo()->Add(fESD);
    }
  }
  if (!fpESDfile || !fpESDtree) return -ENOMEM;

  fpESDtree->Fill();
  fpESDfile->cd();
  fpESDtree->Write(fpESDtree->GetName(),TObject::kOverwrite);

  return iResult;
}
