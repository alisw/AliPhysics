// $Id$

//**************************************************************************
//* This file is property of and copyright by the                          * 
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
//**************************************************************************/

/// @file   AliHLTSampleMonitoringComponent.cxx
/// @author Matthias Richter
/// @date   
/// @brief  A sample monitoring component for the HLT.
///

#include "AliHLTSampleMonitoringComponent.h"
#include "TH1F.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "AliCDBEntry.h"
#include "AliCDBManager.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTSampleMonitoringComponent)

AliHLTSampleMonitoringComponent::AliHLTSampleMonitoringComponent()
  :
  fPushHistograms(false),
  fPushTTree(false),
  fPushTObjArray(false),
  fHpx(NULL),
  fHpy(NULL)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTSampleMonitoringComponent::~AliHLTSampleMonitoringComponent()
{
  // see header file for class documentation
}

const char* AliHLTSampleMonitoringComponent::GetComponentID()
{
  // see header file for class documentation
  return "Sample-MonitoringComponent";
}

void AliHLTSampleMonitoringComponent::GetInputDataTypes( AliHLTComponentDataTypeList& list)
{
  // see header file for class documentation
  list.push_back(kAliHLTAnyDataType);
}

AliHLTComponentDataType AliHLTSampleMonitoringComponent::GetOutputDataType()
{
  // see header file for class documentation
  return kAliHLTVoidDataType;
}

void AliHLTSampleMonitoringComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier )
{
  // see header file for class documentation
  constBase = 100000;
  inputMultiplier = 0;
}

AliHLTComponent* AliHLTSampleMonitoringComponent::Spawn()
{
  // see header file for class documentation
  return new AliHLTSampleMonitoringComponent;
}

int AliHLTSampleMonitoringComponent::DoInit( int argc, const char** argv )
{
  // see header file for class documentation
  int iResult=0;

  TString argument="";
  TString configuration=""; 
  int bMissingParam=0;

  for (int i=0; i<argc && iResult>=0; i++) {
    argument=argv[i];
    if (argument.IsNull()) continue;

    // -push-histograms
    if (argument.CompareTo("-push-histograms")==0) {
      fPushHistograms=true;

    // -push-ttree
    } else if (argument.CompareTo("-push-ttree")==0) {
      fPushTTree=true;

    // -push-array
    } else if (argument.CompareTo("-push-array")==0) {
      fPushTObjArray=true;

    } else {
      // the remaining arguments are treated as configuration
      if (!configuration.IsNull()) configuration+=" ";
      configuration+=argument;
    }
  }
  if (bMissingParam) {
    HLTError("missing parameter for argument %s", argument.Data());
    iResult=-EINVAL;
  }

  // choose fPushTTree as default if none is set
  if (!(fPushTTree || fPushTObjArray || fPushHistograms)) fPushTTree=true;

  // strictly speaking I would prefer to use local or dynamic variables
  // locally in DoEvent, but there is a ROOT bug or feature (related to
  // garbage collection) which causes seg faults after a while. 
  fHpx = new TH1F("hpx","px distribution",100,-4,4);
  fHpy = new TH1F("hpy","py distribution",100,-10,10);

  return iResult;
}

int AliHLTSampleMonitoringComponent::DoDeinit()
{
  // see header file for class documentation
  if (fHpx) delete fHpx;
  fHpx=NULL;
  if (fHpy) delete fHpy;
  fHpy=NULL;
  return 0;
}

int AliHLTSampleMonitoringComponent::DoEvent( const AliHLTComponentEventData& /*evtData*/, AliHLTComponentTriggerData& /*trigData*/)
{
  // see header file for class documentation

  // the function ignores all input blocks and fakes some monitoring histogram
  int iResult=0;

  fHpx->Reset();
  fHpx->FillRandom("gaus",100*(GetEventCount()+1));

  fHpy->Reset();
  fHpy->FillRandom("gaus",500*(GetEventCount()+1));

  if (fPushHistograms) {
    PushBack(fHpx, "ROOTTH1F", "EXPL", 0);
    PushBack(fHpy, "ROOTTH1F", "EXPL", 1);
  }

  if (fPushTTree) {
    TString event;
    TTree *pTree = new TTree("T","A Root Tree");
    if (pTree) {
      pTree->SetDirectory(0);
      event.Form("event_%d_hpx", GetEventCount());
      pTree->Branch(event, "TH1F", &fHpx, 32000, 0);
      event.Form("event_%d_hpy", GetEventCount());
      pTree->Branch(event, "TH1F", &fHpy, 32000, 0);

      PushBack(pTree, "ROOTTREE", "EXPL");
      delete pTree;
    } else {
      iResult=-ENOMEM;
    }
  }

  if (fPushTObjArray) {
    TObjArray* pArray=new TObjArray;
    if (pArray) {
      pArray->Add(fHpx);
      pArray->Add(fHpy);
      
      PushBack(pArray, "ROOTOBJA", "EXPL");
      delete pArray;
    }
  }

  return iResult;
}

int AliHLTSampleMonitoringComponent::Configure(const char* arguments)
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

      // -reset
      if (argument.CompareTo("-reset")==0) {

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

int AliHLTSampleMonitoringComponent::Reconfigure(const char* cdbEntry, const char* /*chainId*/)
{
  // see header file for class documentation
  int iResult=0;
  const char* path=NULL;
#ifdef __DEBUG
  const char* defaultNotify="";
#endif
  if (cdbEntry) {
    path=cdbEntry;
#ifdef __DEBUG
    defaultNotify=" (default)";
#endif
  }
  if (path) {
    AliCDBEntry *pEntry = AliCDBManager::Instance()->Get(path/*,GetRunNo()*/);
    if (pEntry) {
      TObjString* pString=dynamic_cast<TObjString*>(pEntry->GetObject());
      if (pString) {
	iResult=Configure(pString->GetString().Data());
      } else {
	HLTError("configuration object \"%s\" has wrong type, required TObjString", path);
      }
    } else {
      HLTError("can not fetch object \"%s\" from CDB", path);
    }
  }
  return iResult;
}
