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

//  @file   AliHLTReconstructor.cxx
//  @author Matthias Richter
//  @date   
//  @brief  Binding class for HLT reconstruction in AliRoot
//          Implements bot the interface to run HLT chains embedded into
//          AliReconstruction and the unpacking and treatment of HLTOUT

#include <TSystem.h>
#include <TObjString.h>
#include "TFile.h"
#include "TTree.h"
#include "TObject.h"
#include "TObjArray.h"
#include "TClass.h"
#include "TStreamerInfo.h"
#include "AliHLTReconstructor.h"
#include "AliLog.h"
#include "AliRawReader.h"
#include "AliESDEvent.h"
#include "AliHLTSystem.h"
#include "AliHLTOUTRawReader.h"
#include "AliHLTOUTDigitReader.h"
#include "AliHLTEsdManager.h"
#include "AliHLTPluginBase.h"
#include "AliHLTMisc.h"
#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliHLTMessage.h"
#include "AliCentralTrigger.h"
#include "AliTriggerConfiguration.h"
#include "AliTriggerClass.h"
#include "AliTriggerCluster.h"
#include "AliDAQ.h"
#include "AliRunLoader.h"

class AliCDBEntry;

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTReconstructor)

AliHLTReconstructor::AliHLTReconstructor()
  : AliReconstructor()
  , fpEsdManager(NULL)
  , fpPluginBase(new AliHLTPluginBase)
  , fFlags(0)
{ 
  //constructor
}

AliHLTReconstructor::AliHLTReconstructor(const char* options)
  : AliReconstructor()
  , fpEsdManager(NULL)
  , fpPluginBase(new AliHLTPluginBase)
  , fFlags(0)
{ 
  //constructor
  if (options) Init(options);
}

AliHLTReconstructor::~AliHLTReconstructor()
{ 
  //destructor

  if (fpEsdManager) AliHLTEsdManager::Delete(fpEsdManager);
  fpEsdManager=NULL;

  if (fpPluginBase) {
  AliHLTSystem* pSystem=fpPluginBase->GetInstance();
  if (pSystem) {
    AliDebug(0, Form("terminate HLT system: status %#x", pSystem->GetStatusFlags()));
    if (pSystem->CheckStatus(AliHLTSystem::kStarted)) {
      // send specific 'event' to execute the stop sequence
      pSystem->Reconstruct(0, NULL, NULL);
    }
  }
  delete fpPluginBase;
  }
  fpPluginBase=NULL;

}

void AliHLTReconstructor::Init(const char* options)
{
  // init the reconstructor
  SetOption(options);
  Init();
}

void AliHLTReconstructor::Init()
{
  // init the reconstructor
  if (!fpPluginBase) {
    AliError("internal memory error: can not get AliHLTSystem instance from plugin");
    return;
  }

  AliHLTSystem* pSystem=fpPluginBase->GetInstance();
  if (!pSystem) {
    AliError("can not create AliHLTSystem object");
    return;
  }
  if (pSystem->CheckStatus(AliHLTSystem::kError)) {
    AliError("HLT system in error state");
    return;
  }

  TString esdManagerOptions;

  // the options scan has been moved to AliHLTSystem, the old code
  // here is kept to be able to run an older version of the HLT code
  // with newer AliRoot versions.
  TString option = GetOption();
  TObjArray* pTokens=option.Tokenize(" ");
  option="";
  if (pTokens) {
    int iEntries=pTokens->GetEntries();
    for (int i=0; i<iEntries; i++) {
      TString token=(((TObjString*)pTokens->At(i))->GetString());
      if (token.Contains("loglevel=")) {
	TString param=token.ReplaceAll("loglevel=", "");
	if (param.IsDigit()) {
	  pSystem->SetGlobalLoggingLevel((AliHLTComponentLogSeverity)param.Atoi());
	} else if (param.BeginsWith("0x") &&
		   param.Replace(0,2,"",0).IsHex()) {
	  int severity=0;
	  sscanf(param.Data(),"%x", &severity);
	  pSystem->SetGlobalLoggingLevel((AliHLTComponentLogSeverity)severity);
	} else {
	  AliWarning("wrong parameter for option \'loglevel=\', (hex) number expected");
	}
      } else if (token.Contains("alilog=off")) {
	pSystem->SwitchAliLog(0);
      } else if (token.CompareTo("ignore-hltout")==0) {
	fFlags|=kAliHLTReconstructorIgnoreHLTOUT;
      } else if (token.CompareTo("run-online-config")==0) {
        fFlags|=kAliHLTReconstructorIgnoreHLTOUT;
	if (option.Length()>0) option+=" ";
	option+=token;
      } else if (token.CompareTo("ignore-ctp")==0) {
	fFlags|=kAliHLTReconstructorIgnoreCTP;
      } else if (token.Contains("esdmanager=")) {
	token.ReplaceAll("esdmanager=", "");
	token.ReplaceAll(","," ");
	token.ReplaceAll("'","");
	esdManagerOptions=token;
      } else {
	if (option.Length()>0) option+=" ";
	option+=token;
      }
    }
    delete pTokens;
  }

  TString ecsParam;
  TString ctpParam;
  if ((fFlags&kAliHLTReconstructorIgnoreCTP)==0 &&
      BuildCTPTriggerClassString(ctpParam)>=0) {
    if (!ecsParam.IsNull()) ecsParam+=";";
    ecsParam+="CTP_TRIGGER_CLASS=";
    ecsParam+=ctpParam;
  }

  if (!ecsParam.IsNull()) {
    option+=" ECS=";
    option+=ecsParam;
  }

  if (!pSystem->CheckStatus(AliHLTSystem::kReady)) {
    if (pSystem->ScanOptions(option.Data())<0) {
      AliError("error setting options for HLT system");
      return;
    }
    if ((pSystem->Configure())<0) {
      AliError("error during HLT system configuration");
      return;
    }
  }

  fpEsdManager=AliHLTEsdManager::New();
  if (fpEsdManager) {
    fpEsdManager->SetOption(esdManagerOptions.Data());
  }

  AliHLTMisc::Instance().InitStreamerInfos(fgkCalibStreamerInfoEntry);
}

void AliHLTReconstructor::Terminate() const
{
  AliInfo("terminating");
  if (fpPluginBase) {
  AliHLTSystem* pSystem=fpPluginBase->GetInstance();
  if (pSystem) {
    AliDebug(0, Form("terminate HLT system: status %#x", pSystem->GetStatusFlags()));
    if (pSystem->CheckStatus(AliHLTSystem::kStarted)) {
      // send specific 'event' to execute the stop sequence
      pSystem->Reconstruct(0, NULL, NULL);
    }
  }
  }
}

const char* AliHLTReconstructor::fgkCalibStreamerInfoEntry="HLT/Calib/StreamerInfo";

void AliHLTReconstructor::Reconstruct(AliRawReader* rawReader, TTree* /*clustersTree*/) const 
{
  // reconstruction of real data without writing of ESD
  // For each event, HLT reconstruction chains can be executed and
  // added to the existing HLTOUT data
  // The HLTOUT data is finally processed in FillESD

  if (!fpPluginBase) {
    AliError("internal memory error: can not get AliHLTSystem instance from plugin");
    return;
  }

  int iResult=0;
  AliHLTSystem* pSystem=fpPluginBase->GetInstance();

  if (pSystem) {
    AliHLTOUT* pHLTOUT=NULL;
    pSystem->InvalidateHLTOUT(&pHLTOUT);
    if (pHLTOUT) {
      delete pHLTOUT;
      pHLTOUT=NULL;
    }
    if (pSystem->CheckStatus(AliHLTSystem::kError)) {
      AliError("HLT system in error state");
      return;
    }
    if (!pSystem->CheckStatus(AliHLTSystem::kReady)) {
      AliError("HLT system in wrong state");
      return;
    }

    // init the HLTOUT instance for the current event
    // not nice. Have to query the global run loader to get the current event no.
    Int_t eventNo=-1;
    AliRunLoader* runloader = AliRunLoader::Instance();
    if (runloader) {
      eventNo=runloader->GetEventNumber();
    }
    if (eventNo>=0) {
      AliRawReader* input=NULL;
      if ((fFlags&kAliHLTReconstructorIgnoreHLTOUT) == 0 ) {
	input=rawReader;
      }
      pHLTOUT=new AliHLTOUTRawReader(input, eventNo, fpEsdManager);
      if (pHLTOUT) {
	if (pHLTOUT->Init()>=0) {
	  pSystem->InitHLTOUT(pHLTOUT);
	} else {
	  AliError("error : initialization of HLTOUT handler failed");
	}
      } else {
	AliError("memory allocation failed: can not create AliHLTOUT object");
      }
    } else {
      AliError("can not get event number");
    }

    if ((iResult=pSystem->Reconstruct(1, NULL, rawReader))>=0) {
    }
  }
}

void AliHLTReconstructor::FillESD(AliRawReader* rawReader, TTree* /*clustersTree*/, 
				  AliESDEvent* esd) const
{
  // reconstruct real data and fill ESD
  if (!rawReader || !esd) {
    AliError("missing raw reader or esd object");
    return;
  }

  if (!fpPluginBase) {
    AliError("internal memory error: can not get AliHLTSystem instance from plugin");
    return;
  }

  AliHLTSystem* pSystem=fpPluginBase->GetInstance();

  if (pSystem) {
    if (pSystem->CheckStatus(AliHLTSystem::kError)) {
      AliError("HLT system in error state");
      return;
    }
    if (!pSystem->CheckStatus(AliHLTSystem::kReady)) {
      AliError("HLT system in wrong state");
      return;
    }
    pSystem->FillESD(-1, NULL, esd);

    // the HLTOUT handler has either been created in the AliHLTReconstructor::Reconstruct
    // step of this event or is created now. In either case the instance is deleted after
    // the processing
    AliHLTOUT* pHLTOUT=NULL;
    pSystem->InvalidateHLTOUT(&pHLTOUT);
    if (!pHLTOUT) {
      AliRawReader* input=NULL;
      if ((fFlags&kAliHLTReconstructorIgnoreHLTOUT) == 0 ) {
	input=rawReader;
      }
      pHLTOUT=new AliHLTOUTRawReader(input, esd->GetEventNumberInFile(), fpEsdManager);
    }
    if (pHLTOUT) {
      ProcessHLTOUT(pHLTOUT, esd, (pSystem->GetGlobalLoggingLevel()&kHLTLogDebug)!=0);
      delete pHLTOUT;
    } else {
      AliError("error creating HLTOUT handler");
    }
  }
}

void AliHLTReconstructor::Reconstruct(TTree* /*digitsTree*/, TTree* /*clustersTree*/) const
{
  // reconstruct simulated data

  AliHLTSystem* pSystem=fpPluginBase->GetInstance();

  if (pSystem) {
    // create the HLTOUT instance in order to be available for other detector reconstruction
    // first cleanup any existing instance
    AliHLTOUT* pHLTOUT=NULL;
    pSystem->InvalidateHLTOUT(&pHLTOUT);
    if (pHLTOUT) {
      delete pHLTOUT;
      pHLTOUT=NULL;
    }

    // not nice. Have to query the global run loader to get the current event no.
    // This is related to the missing AliLoader for HLT.
    // Since AliReconstruction can not provide a digits tree, the file needs to be accessed
    // explicitely, and the corresponding event needs to be selected.
    Int_t eventNo=-1;
    AliRunLoader* runloader = AliRunLoader::Instance();
    if (runloader) {
      eventNo=runloader->GetEventNumber();
    }
    if (eventNo>=0) {
      const char* digitfile=NULL;
      if ((fFlags&kAliHLTReconstructorIgnoreHLTOUT) == 0 ) {
	digitfile="HLT.Digits.root";
      }

      pHLTOUT=new AliHLTOUTDigitReader(eventNo, fpEsdManager, digitfile);
      if (pHLTOUT) {
	if (pHLTOUT->Init()>=0) {
	  pSystem->InitHLTOUT(pHLTOUT);
	} else {
	  AliError("error : initialization of HLTOUT handler failed");
	}
      } else {
	AliError("memory allocation failed: can not create AliHLTOUT object");
      }
    } else {
      AliError("can not get event number");
    }

    // all data processing happens in FillESD
  }
}

void AliHLTReconstructor::FillESD(TTree* /*digitsTree*/, TTree* /*clustersTree*/, AliESDEvent* esd) const
{
  // reconstruct simulated data and fill ESD

  // later this is the place to extract the simulated HLT data
  // for now it's only an user failure condition as he tries to run HLT reconstruction
  // on simulated data 
  TString option = GetOption();
  if (!option.IsNull() && 
      (option.Contains("config=") || option.Contains("chains="))) {
    AliWarning(Form("You are trying to run a custom HLT chain on digits data.\n\n"
		    "HLT reconstruction can be run embedded into AliReconstruction from\n"
		    "raw data (real or simulated)). Reconstruction of digit data takes\n"
		    "place in AliSimulation, appropriate input conversion is needed to\n"
		    "feed data from the detector digits into the HLT chain.\n"
		    "Consider running embedded into AliSimulation.\n"
		    "        /***  run macro *****************************************/\n"
		    "        AliSimulation sim;\n"
		    "        sim.SetRunHLT(\"%s\");\n"
		    "        sim.SetRunGeneration(kFALSE);\n"
		    "        sim.SetMakeDigits(\"\");\n"
		    "        sim.SetMakeSDigits(\"\");\n"
		    "        sim.SetMakeDigitsFromHits(\"\");\n"
		    "        sim.Run();\n"
		    "        /*********************************************************/\n\n",
		    option.Data()));
  }
  if (!fpPluginBase) {
    AliError("internal memory error: can not get AliHLTSystem instance from plugin");
    return;
  }

  AliHLTSystem* pSystem=fpPluginBase->GetInstance();
  if (pSystem) {
    if (pSystem->CheckStatus(AliHLTSystem::kError)) {
      AliError("HLT system in error state");
      return;
    }
    if (!pSystem->CheckStatus(AliHLTSystem::kReady)) {
      AliError("HLT system in wrong state");
      return;
    }

    // the HLTOUT handler has either been created in the AliHLTReconstructor::Reconstruct
    // step of this event or is created now. In either case the instance is deleted after
    // the processing
    AliHLTOUT* pHLTOUT=NULL;
    pSystem->InvalidateHLTOUT(&pHLTOUT);
    if (!pHLTOUT) {
      const char* digitfile=NULL;
      if ((fFlags&kAliHLTReconstructorIgnoreHLTOUT) == 0 ) {
	digitfile="HLT.Digits.root";
      }
      pHLTOUT=new AliHLTOUTDigitReader(esd->GetEventNumberInFile(), fpEsdManager, digitfile);
    }

    if (pHLTOUT) {
      ProcessHLTOUT(pHLTOUT, esd, (pSystem->GetGlobalLoggingLevel()&kHLTLogDebug)!=0);
      delete pHLTOUT;
    } else {
      AliError("error creating HLTOUT handler");
    }
  }
}

void AliHLTReconstructor::ProcessHLTOUT(AliHLTOUT* pHLTOUT, AliESDEvent* esd, bool bVerbose) const
{
  // treatment of simulated or real HLTOUT data
  if (!pHLTOUT) return;
  if (!fpPluginBase) {
    AliError("internal memory error: can not get AliHLTSystem instance from plugin");
    return;
  }

  AliHLTSystem* pSystem=fpPluginBase->GetInstance();
  if (!pSystem) {
    AliError("error getting HLT system instance");
    return;
  }

  if (pHLTOUT->Init()<0) {
    AliError("error : initialization of HLTOUT handler failed");
    return;
  }

  if (bVerbose)
    PrintHLTOUTContent(pHLTOUT);

  int blockindex=pHLTOUT->SelectFirstDataBlock(kAliHLTDataTypeStreamerInfo);
  if (blockindex>=0) {
    const AliHLTUInt8_t* pBuffer=NULL;
    AliHLTUInt32_t size=0;
    if (pHLTOUT->GetDataBuffer(pBuffer, size)>=0) {
      TObject* pObject=AliHLTMessage::Extract(pBuffer, size);
      if (pObject) {
	TObjArray* pArray=dynamic_cast<TObjArray*>(pObject);
	if (pArray) {
	  AliHLTMisc::Instance().InitStreamerInfos(pArray);
	} else {
	  AliError(Form("wrong class type of streamer info list: expected TObjArray, but object is of type %s", pObject->Class()->GetName()));
	}
      } else {
	AliError(Form("failed to extract object from data block of type %s", AliHLTComponent::DataType2Text(kAliHLTDataTypeStreamerInfo).c_str()));
      }
    } else {
      AliError(Form("failed to get data buffer for block of type %s", AliHLTComponent::DataType2Text(kAliHLTDataTypeStreamerInfo).c_str()));
    }
  }

  if (pSystem->ProcessHLTOUT(pHLTOUT, esd)<0) {
    AliError("error processing HLTOUT");
  }

  if (bVerbose && esd) {
    AliInfo("HLT ESD content:");
    esd->Print();
  }
  pHLTOUT->Reset();
}

void AliHLTReconstructor::ProcessHLTOUT(const char* digitFile, AliESDEvent* pEsd) const
{
  // debugging/helper function to examine simulated data
  if (!digitFile) return;

  // read the number of events
  TFile f(digitFile);
  if (f.IsZombie()) return;
  TTree* pTree=NULL;
  f.GetObject("rawhltout", pTree);
  if (!pTree) {
    AliWarning(Form("can not find tree rawhltout in file %s", digitFile));
    return ;
  }
  int nofEvents=pTree->GetEntries();
  f.Close();
  //delete pTree; OF COURSE NOT! its an object in the file
  pTree=NULL;

  for (int event=0; event<nofEvents; event++) {
    AliHLTOUTDigitReader* pHLTOUT=new AliHLTOUTDigitReader(event, fpEsdManager, digitFile);
    if (pHLTOUT) {
      AliInfo(Form("event %d", event));
      ProcessHLTOUT(pHLTOUT, pEsd, true);
      delete pHLTOUT;
    } else {
      AliError("error creating HLTOUT handler");
    }
  }
}

void AliHLTReconstructor::ProcessHLTOUT(AliRawReader* pRawReader, AliESDEvent* pEsd) const
{
  // debugging/helper function to examine simulated or real HLTOUT data
  if (!pRawReader) return;

  pRawReader->RewindEvents();
  for (int event=0; pRawReader->NextEvent(); event++) {
    AliHLTOUTRawReader* pHLTOUT=new AliHLTOUTRawReader(pRawReader, event, fpEsdManager);
    if (pHLTOUT) {
      AliInfo(Form("event %d", event));
      // the two event fields contain: period - orbit - bunch crossing counter
      //        id[0]               id[1]
      // |32                0|32                0|
      //
      // |      28 bit    |       24 bit     | 12|
      //        period          orbit         bcc
      AliHLTUInt64_t eventId=0;
      const UInt_t* rawreaderEventId=pRawReader->GetEventId();
      if (rawreaderEventId) {
	eventId=rawreaderEventId[0];
	eventId=eventId<<32;
	eventId|=rawreaderEventId[1];
      }
      AliInfo(Form("Event Id from rawreader:\t 0x%016llx", eventId));
      ProcessHLTOUT(pHLTOUT, pEsd, true);
      delete pHLTOUT;
    } else {
      AliError("error creating HLTOUT handler");
    }
  }
}

void AliHLTReconstructor::PrintHLTOUTContent(AliHLTOUT* pHLTOUT) const
{
  // print the block specifications of the HLTOUT data blocks
  if (!pHLTOUT) return;
  int iResult=0;

  AliInfo(Form("Event Id from hltout:\t 0x%016llx", pHLTOUT->EventId()));
  for (iResult=pHLTOUT->SelectFirstDataBlock();
       iResult>=0;
       iResult=pHLTOUT->SelectNextDataBlock()) {
    AliHLTComponentDataType dt=kAliHLTVoidDataType;
    AliHLTUInt32_t spec=kAliHLTVoidDataSpec;
    pHLTOUT->GetDataBlockDescription(dt, spec);
    const AliHLTUInt8_t* pBuffer=NULL;
    AliHLTUInt32_t size=0;
    if (pHLTOUT->GetDataBuffer(pBuffer, size)>=0) {
      pHLTOUT->ReleaseDataBuffer(pBuffer);
      pBuffer=NULL; // just a dummy
    }
    AliInfo(Form("   %s  0x%x: size %d", AliHLTComponent::DataType2Text(dt).c_str(), spec, size));
  }
}

int AliHLTReconstructor::BuildCTPTriggerClassString(TString& triggerclasses) const
{
  // build the CTP trigger class string from the OCDB entry of the CTP trigger
  int iResult=0;
  
  triggerclasses.Clear();
  AliCentralTrigger* pCTP = new AliCentralTrigger();
  AliTriggerConfiguration *config=NULL;
  TString configstr("");
  if (pCTP->LoadConfiguration(configstr) && 
      (config = pCTP->GetConfiguration())!=NULL) {
    const TObjArray& classesArray = config->GetClasses();
    int nclasses = classesArray.GetEntriesFast();
    for( int iclass=0; iclass < nclasses; iclass++ ) {
      AliTriggerClass* trclass = NULL;
      if (classesArray.At(iclass) && (trclass=dynamic_cast<AliTriggerClass*>(classesArray.At(iclass)))!=NULL) {
	TString entry;
	int trindex = TMath::Nint(TMath::Log2(trclass->GetMask()));
	entry.Form("%02d:%s:", trindex, trclass->GetName());
	AliTriggerCluster* cluster=NULL;
	TObject* clusterobj=config->GetClusters().FindObject(trclass->GetCluster());
	if (clusterobj && (cluster=dynamic_cast<AliTriggerCluster*>(clusterobj))!=NULL) {
	  TString detectors=cluster->GetDetectorsInCluster();
	  TObjArray* pTokens=detectors.Tokenize(" ");
	  if (pTokens) {
	    for (int dix=0; dix<pTokens->GetEntriesFast(); dix++) {
	      int id=AliDAQ::DetectorID(((TObjString*)pTokens->At(dix))->GetString());
	      if (id>=0) {
		TString detstr; detstr.Form("%s%02d", dix>0?"-":"", id);
		entry+=detstr;
	      } else {
		AliError(Form("invalid detector name extracted from trigger cluster: %s (%s)", ((TObjString*)pTokens->At(dix))->GetString().Data(), detectors.Data()));
		iResult=-EPROTO;
		break;
	      }
	    }
	    delete pTokens;
	  }
	} else {
	  AliError(Form("can not find trigger cluster %s in config", trclass->GetCluster()?trclass->GetCluster()->GetName():"NULL"));
	  iResult=-EPROTO;
	  break;
	}
	if (!triggerclasses.IsNull()) triggerclasses+=",";
	triggerclasses+=entry;
      }
    }
  }

  return iResult;
}
