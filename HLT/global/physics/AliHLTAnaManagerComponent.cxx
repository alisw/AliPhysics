//-*- Mode: C++ -*-
// $Id: AliHLTAnaManagerComponent.cxx $
/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 *                                                                        *
 * Primary Authors: David Rohr, Jens Wiechula, C. Zampolli                *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/** @file    AliHLTAnaManagerComponent.cxx
    @author  David Rohr, Jens Wiechula, C. Zampolli
    @brief   Component for Testing Analysis Manager inside HLT component
*/

#include "TMap.h"
#include "TSystem.h"
#include "TTimeStamp.h"
#include "TObjString.h"
#include "TH1F.h"
#include "TList.h"
#include "AliESDtrackCuts.h"
#include "AliHLTErrorGuard.h"
#include "AliHLTDataTypes.h"
#include "AliHLTAnaManagerComponent.h"
#include "AliHLTITSClusterDataFormat.h"
#include "AliAnalysisManager.h"
#include "AliHLTTestInputHandler.h"
#include "AliAnalysisTaskPt.h"
#include "AliAnalysisDataContainer.h"
#include "TTree.h"

using namespace std;

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTAnaManagerComponent)

/*
 * ---------------------------------------------------------------------------------
 *                            Constructor / Destructor
 * ---------------------------------------------------------------------------------
 */

// #################################################################################
AliHLTAnaManagerComponent::AliHLTAnaManagerComponent() :
  AliHLTProcessor(),
  fUID(0),
  fAnalysisManager(NULL),
  fInputHandler(NULL) {
  // an example component which implements the ALICE HLT processor
  // interface and does some analysis on the input raw data
  //
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
  //
  // NOTE: all helper classes should be instantiated in DoInit()
}

// #################################################################################
AliHLTAnaManagerComponent::~AliHLTAnaManagerComponent() {
  // see header file for class documentation
}

/*
 * ---------------------------------------------------------------------------------
 * Public functions to implement AliHLTComponent's interface.
 * These functions are required for the registration process
 * ---------------------------------------------------------------------------------
 */

// #################################################################################
const Char_t* AliHLTAnaManagerComponent::GetComponentID() { 
  // see header file for class documentation
  return "AnaManagerComponent";
}

// #################################################################################
void AliHLTAnaManagerComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list) {
  // see header file for class documentation
  list.push_back(kAliHLTDataTypeESDObject|kAliHLTDataOriginAny);
  list.push_back(kAliHLTDataTypeClusters|kAliHLTDataOriginITSSPD);
  list.push_back(kAliHLTDataTypeESDContent|kAliHLTDataOriginVZERO);
}

// #################################################################################
AliHLTComponentDataType AliHLTAnaManagerComponent::GetOutputDataType() {
  // see header file for class documentation
  return kAliHLTDataTypeTObject|kAliHLTDataOriginHLT;
}

// #################################################################################
void AliHLTAnaManagerComponent::GetOutputDataSize( ULong_t& constBase, Double_t& inputMultiplier ) {
  // see header file for class documentation
  constBase = 100000;
  inputMultiplier = 0.5;
}

// #################################################################################
void AliHLTAnaManagerComponent::GetOCDBObjectDescription( TMap* const targetMap) {
  // see header file for class documentation

  if (!targetMap) return;
  targetMap->Add(new TObjString("HLT/ConfigGlobal/MultiplicityCorrelations"),
		 new TObjString("configuration object"));
  targetMap->Add(new TObjString("HLT/ConfigGlobal/MultiplicityCorrelationsCentrality"),
		 new TObjString("centrality configuration object"));

  return;
}

// #################################################################################
AliHLTComponent* AliHLTAnaManagerComponent::Spawn() {
  // see header file for class documentation
  return new AliHLTAnaManagerComponent;
}

/*
 * ---------------------------------------------------------------------------------
 * Protected functions to implement AliHLTComponent's interface.
 * These functions provide initialization as well as the actual processing
 * capabilities of the component. 
 * ---------------------------------------------------------------------------------
 */

// #################################################################################
Int_t AliHLTAnaManagerComponent::DoInit( Int_t argc, const Char_t** argv ) {
  // see header file for class documentation
  printf("AliHLTAnaManagerComponent::DoInit\n");

  Int_t iResult=0;

  fAnalysisManager = new AliAnalysisManager;
  fInputHandler    = new AliHLTTestInputHandler;
  fAnalysisManager->SetInputEventHandler(fInputHandler);
  fAnalysisManager->SetExternalLoop(kTRUE); 

  AliAnalysisTaskPt *task = new AliAnalysisTaskPt("TaskPt");
  fAnalysisManager->AddTask(task);
  AliAnalysisDataContainer *cinput  = fAnalysisManager->GetCommonInputContainer();
  Printf("Defining output file");
  AliAnalysisDataContainer *coutput1 = fAnalysisManager->CreateContainer("pt", TH1::Class(),
      AliAnalysisManager::kOutputContainer, "Pt.ESD.root");

  //           connect containers
  fAnalysisManager->ConnectInput  (task,  0, cinput );
  Printf("---> Connecting output...");
  fAnalysisManager->ConnectOutput (task,  0, coutput1);
  Printf("---> ...connected.");

  fAnalysisManager->InitAnalysis();
  fAnalysisManager->StartAnalysis("local", (TTree*)new TTree);

  return iResult;
}



// #################################################################################
Int_t AliHLTAnaManagerComponent::DoDeinit() {
  // see header file for class documentation

  fUID = 0;
  fAnalysisManager->SetSkipTerminate(kTRUE);
  fAnalysisManager->Terminate();

  delete fAnalysisManager;

  return 0;
}

// #################################################################################
Int_t AliHLTAnaManagerComponent::DoEvent(const AliHLTComponentEventData& evtData,
					AliHLTComponentTriggerData& /*trigData*/) {
  // see header file for class documentation

  printf("AliHLTAnaManagerComponent::DoEvent\n");
  Int_t iResult=0;

  // -- Only use data event
  if (!IsDataEvent()) 
    return 0;
  
  if( fUID == 0 ){
    TTimeStamp t;
    fUID = ( gSystem->GetPid() + t.GetNanoSec())*10 + evtData.fEventID;
  }
  
  // -- Get ESD object
  // -------------------
  AliESDEvent *esdEvent = NULL;
  for ( const TObject *iter = GetFirstInputObject(kAliHLTDataTypeESDObject); iter != NULL; iter = GetNextInputObject() ) {
    esdEvent = dynamic_cast<AliESDEvent*>(const_cast<TObject*>( iter ) );
    if( !esdEvent ){ 
      HLTWarning("Wrong ESDEvent object received");
      iResult = -1;
      continue;
    }
    esdEvent->GetStdContent();
  }
  printf("ESDEvent: %p\n",esdEvent);

  fInputHandler->SetEvent(esdEvent);
  fInputHandler->BeginEvent(0);
  fAnalysisManager->ExecAnalysis();
  fInputHandler->FinishEvent();


  // -- Send histlist
//  PushBack(dynamic_cast<TObject*>(fCorrObj->GetHistList()),
//	   kAliHLTDataTypeTObject|kAliHLTDataOriginHLT,fUID);
 
  return iResult;
}

// #################################################################################
Int_t AliHLTAnaManagerComponent::Reconfigure(const Char_t* cdbEntry, const Char_t* chainId) {
  // see header file for class documentation

  Int_t iResult=0;
  TString cdbPath;
  if (cdbEntry) {
    cdbPath=cdbEntry;
  } else {
    cdbPath="HLT/ConfigGlobal/";
    cdbPath+=GetComponentID();
  }

  AliInfoClass(Form("reconfigure '%s' from entry %s%s", chainId, cdbPath.Data(), cdbEntry?"":" (default)"));
  iResult=ConfigureFromCDBTObjString(cdbPath);

  return iResult;
}

// #################################################################################
Int_t AliHLTAnaManagerComponent::ReadPreprocessorValues(const Char_t* /*modules*/) {
  // see header file for class documentation
  ALIHLTERRORGUARD(5, "ReadPreProcessorValues not implemented for this component");
  return 0;
}

