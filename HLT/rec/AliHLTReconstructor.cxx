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

/** @file   AliHLTReconstructor.cxx
    @author Matthias Richter
    @date   
    @brief  Binding class for HLT reconstruction in AliRoot. */

#include <TSystem.h>
#include <TObjString.h>
#include "TFile.h"
#include "TTree.h"
#include "AliHLTReconstructor.h"
#include "AliLog.h"
#include "AliRawReader.h"
#include "AliESDEvent.h"
#include "AliHLTSystem.h"
#include "AliHLTOUTRawReader.h"
#include "AliHLTOUTDigitReader.h"
#include "AliHLTEsdManager.h"

ClassImp(AliHLTReconstructor)

AliHLTReconstructor::AliHLTReconstructor()
  : 
  AliReconstructor(),
  AliHLTReconstructorBase(),
  fFctProcessHLTOUT(NULL),
  fpEsdManager(NULL)
{ 
  //constructor
}

AliHLTReconstructor::AliHLTReconstructor(const char* options)
  : 
  AliReconstructor(),
  AliHLTReconstructorBase(),
  fFctProcessHLTOUT(NULL),
  fpEsdManager(NULL)
{ 
  //constructor
  if (options) Init(options);
}

AliHLTReconstructor::~AliHLTReconstructor()
{ 
  //destructor

  AliHLTSystem* pSystem=GetInstance();
  if (pSystem) {
    AliDebug(0, Form("delete HLT system: status %#x", pSystem->GetStatusFlags()));
    if (pSystem->CheckStatus(AliHLTSystem::kStarted)) {
      // send specific 'event' to execute the stop sequence
      pSystem->Reconstruct(0, NULL, NULL);
    }
  }

  if (fpEsdManager) AliHLTEsdManager::Delete(fpEsdManager);
  fpEsdManager=NULL;
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
  AliHLTSystem* pSystem=GetInstance();
  if (!pSystem) {
    AliError("can not create AliHLTSystem object");
    return;
  }
  if (pSystem->CheckStatus(AliHLTSystem::kError)) {
    AliError("HLT system in error state");
    return;
  }

  // the options scan has been moved to AliHLTSystem, the old code
  // here is kept to be able to run an older version of the HLT code
  // with newer AliRoot versions.
  TString libs("");
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
      } else if (token.BeginsWith("lib") && token.EndsWith(".so")) {
	libs+=token;
	libs+=" ";
      } else {
	if (option.Length()>0) option+=" ";
	option+=token;
      }
    }
    delete pTokens;
  }

  if (!libs.IsNull() &&
      (!pSystem->CheckStatus(AliHLTSystem::kLibrariesLoaded)) &&
      (pSystem->LoadComponentLibraries(libs.Data())<0)) {
    AliError("error while loading HLT libraries");
    return;
  }

  if (!pSystem->CheckStatus(AliHLTSystem::kReady)) {
    typedef int (*AliHLTSystemSetOptions)(AliHLTSystem* pInstance, const char* options);
    gSystem->Load("libHLTinterface.so");
    AliHLTSystemSetOptions pFunc=(AliHLTSystemSetOptions)(gSystem->DynFindSymbol("libHLTinterface.so", "AliHLTSystemSetOptions"));
    if (pFunc) {
      if ((pFunc)(pSystem, option.Data())<0) {
      AliError("error setting options for HLT system");
      return;	
      }
    } else if (option.Length()>0) {
      AliError(Form("version of HLT system does not support the options \'%s\'", option.Data()));
      return;
    }
    if ((pSystem->Configure())<0) {
      AliError("error during HLT system configuration");
      return;
    }
  }

  gSystem->Load("libHLTinterface.so");
  fFctProcessHLTOUT=(void (*)())gSystem->DynFindSymbol("libHLTinterface.so", "AliHLTSystemProcessHLTOUT");

  fpEsdManager=AliHLTEsdManager::New();
}

void AliHLTReconstructor::Reconstruct(AliRawReader* rawReader, TTree* /*clustersTree*/) const 
{
  // reconstruction of real data without writing of ESD
  // For each event, HLT reconstruction chains can be executed and
  // added to the existing HLTOUT data
  // The HLTOUT data is finally processed in FillESD
  int iResult=0;
  AliHLTSystem* pSystem=GetInstance();

  if (pSystem) {
    if (pSystem->CheckStatus(AliHLTSystem::kError)) {
      AliError("HLT system in error state");
      return;
    }
    if (!pSystem->CheckStatus(AliHLTSystem::kReady)) {
      AliError("HLT system in wrong state");
      return;
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

  AliHLTSystem* pSystem=GetInstance();

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

    AliHLTOUTRawReader* pHLTOUT=new AliHLTOUTRawReader(rawReader, esd->GetEventNumberInFile(), fpEsdManager);
    if (pHLTOUT) {
      ProcessHLTOUT(pHLTOUT, esd);
      delete pHLTOUT;
    } else {
      AliError("error creating HLTOUT handler");
    }
  }
}

void AliHLTReconstructor::Reconstruct(TTree* /*digitsTree*/, TTree* /*clustersTree*/) const
{
  // reconstruct simulated data

  // all reconstruction has been moved to FillESD
  //AliReconstructor::Reconstruct(digitsTree,clustersTree);
  AliInfo("running digit data reconstruction");
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
    AliWarning(Form("HLT reconstruction can be run embedded into Alireconstruction from\n"
		    "raw data (real or simulated)). Reconstruction of of digit data takes\n"
		    "place in AliSimulation, appropriate input conversion is needed.\n"
		    "Consider running embedded into AliSimulation."
		    "        /***  run macro *****************************************/\n"
		    "        AliSimulation sim;\n"
		    "        sim.SetRunHLT(\"%s\");\n"
		    "        sim.SetRunGeneration(kFALSE);\n"
		    "        sim.SetMakeDigits(\"\");\n"
		    "        sim.SetMakeSDigits(\"\");\n"
		    "        sim.SetMakeDigitsFromHits(\"\");\n"
		    "        sim.Run();\n"
		    "        /*********************************************************/", option.Data()));
  }
  AliHLTSystem* pSystem=GetInstance();
  if (pSystem) {
    if (pSystem->CheckStatus(AliHLTSystem::kError)) {
      AliError("HLT system in error state");
      return;
    }
    if (!pSystem->CheckStatus(AliHLTSystem::kReady)) {
      AliError("HLT system in wrong state");
      return;
    }

    AliHLTOUTDigitReader* pHLTOUT=new AliHLTOUTDigitReader(esd->GetEventNumberInFile(), fpEsdManager);
    if (pHLTOUT) {
      ProcessHLTOUT(pHLTOUT, esd);
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
  AliHLTSystem* pSystem=GetInstance();
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

  if (fFctProcessHLTOUT) {
    typedef int (*AliHLTSystemProcessHLTOUT)(AliHLTSystem* pInstance, AliHLTOUT* pHLTOUT, AliESDEvent* esd);
    AliHLTSystemProcessHLTOUT pFunc=(AliHLTSystemProcessHLTOUT)fFctProcessHLTOUT;
    if ((pFunc)(pSystem, pHLTOUT, esd)<0) {
      AliError("error processing HLTOUT");
    }
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
