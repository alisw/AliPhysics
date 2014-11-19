//-*- Mode: C++ -*-
// $Id: AliHLTFlatAnaManagerComponent.cxx $
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

/** @file    AliHLTFlatAnaManagerComponent.cxx
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
#include "AliFlatESDEvent.h"
#include "AliESDfriend.h"
#include "AliHLTErrorGuard.h"
#include "AliHLTDataTypes.h"
#include "AliHLTFlatAnaManagerComponent.h"
#include "AliHLTITSClusterDataFormat.h"
#include "AliAnalysisManager.h"
#include "AliHLTVEventInputHandler.h"
#include "AliAnalysisTaskPt.h"
#include "AliAnalysisDataContainer.h"
#include "TTree.h"

using namespace std;

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTFlatAnaManagerComponent)

/*
 * ---------------------------------------------------------------------------------
 *                            Constructor / Destructor
 * ---------------------------------------------------------------------------------
 */

// #################################################################################
AliHLTFlatAnaManagerComponent::AliHLTFlatAnaManagerComponent() :
  AliHLTProcessor(),
  fUID(0),
  fAnalysisManager(NULL),
  fInputHandler(NULL){
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
AliHLTFlatAnaManagerComponent::~AliHLTFlatAnaManagerComponent() {
  // see header file for class documentation
}

/*
 * ---------------------------------------------------------------------------------
 * Public functions to implement AliHLTComponent's interface.
 * These functions are required for the registration process
 * ---------------------------------------------------------------------------------
 */

// #################################################################################
const Char_t* AliHLTFlatAnaManagerComponent::GetComponentID() { 
  // see header file for class documentation
  return "FlatAnaManagerComponent";
}

// #################################################################################
void AliHLTFlatAnaManagerComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list) {
  // see header file for class documentation
  list.push_back(kAliHLTDataTypeFlatESD|kAliHLTDataOriginAny);
  //  list.push_back(kAliHLTDataTypeESDfriendObject|kAliHLTDataOriginAny);
}

// #################################################################################
AliHLTComponentDataType AliHLTFlatAnaManagerComponent::GetOutputDataType() {
  // see header file for class documentation
  return kAliHLTDataTypeTObject|kAliHLTDataOriginHLT;
}

// #################################################################################
void AliHLTFlatAnaManagerComponent::GetOutputDataSize( ULong_t& constBase, Double_t& inputMultiplier ) {
  // see header file for class documentation
  constBase = 100000;
  inputMultiplier = 0.5;
}

// #################################################################################
void AliHLTFlatAnaManagerComponent::GetOCDBObjectDescription( TMap* const targetMap) {
  // see header file for class documentation

  if (!targetMap) return;
  /*  targetMap->Add(new TObjString("HLT/ConfigGlobal/MultiplicityCorrelations"),
		 new TObjString("configuration object"));
  targetMap->Add(new TObjString("HLT/ConfigGlobal/MultiplicityCorrelationsCentrality"),
		 new TObjString("centrality configuration object"));
  */
  return;
}

// #################################################################################
AliHLTComponent* AliHLTFlatAnaManagerComponent::Spawn() {
  // see header file for class documentation
  return new AliHLTFlatAnaManagerComponent;
}

/*
 * ---------------------------------------------------------------------------------
 * Protected functions to implement AliHLTComponent's interface.
 * These functions provide initialization as well as the actual processing
 * capabilities of the component. 
 * ---------------------------------------------------------------------------------
 */

// #################################################################################
Int_t AliHLTFlatAnaManagerComponent::DoInit( Int_t /*argc*/, const Char_t** /*argv*/ ) {
  // see header file for class documentation
  printf("AliHLTFlatAnaManagerComponent::DoInit\n");

  Int_t iResult=0;

  Printf("----> AliHLTFlatAnaManagerComponent::DoInit"); 
  fAnalysisManager = new AliAnalysisManager;
  fInputHandler    = new AliHLTVEventInputHandler;
  fAnalysisManager->SetInputEventHandler(fInputHandler);
  fAnalysisManager->SetExternalLoop(kTRUE); 

  AliAnalysisTaskPt *task = new AliAnalysisTaskPt("TaskPt");
  task->SetUseFriends(kFALSE);
  fAnalysisManager->AddTask(task);
  AliAnalysisDataContainer *cinput  = fAnalysisManager->GetCommonInputContainer();
  Printf("Defining output file");
  AliAnalysisDataContainer *coutput1 = fAnalysisManager->CreateContainer("pt", TList::Class(),
      AliAnalysisManager::kOutputContainer, "Pt.ESD.root");

  //           connect containers
  Printf("---> Connecting input...");
  fAnalysisManager->ConnectInput  (task,  0, cinput ); 
  Printf("---> ...connected.");
  Printf("---> Connecting output...");
  fAnalysisManager->ConnectOutput (task,  0, coutput1);
  Printf("---> ...connected.");

  Printf("----> Calling InitAnalysis");
  fAnalysisManager->InitAnalysis();
  Printf("----> Done.");
  Printf("----> Calling StartAnalysis");
  fAnalysisManager->StartAnalysis("local", (TTree*)new TTree);
  //fAnalysisManager->StartAnalysis("local", (TTree*)NULL);
  Printf("----> Done.");

  return iResult;
}



// #################################################################################
Int_t AliHLTFlatAnaManagerComponent::DoDeinit() {
  // see header file for class documentation

  fUID = 0;
  fAnalysisManager->SetSkipTerminate(kTRUE);
  fAnalysisManager->Terminate();

  delete fAnalysisManager;

  return 0;
}

// #################################################################################
Int_t AliHLTFlatAnaManagerComponent::DoEvent(const AliHLTComponentEventData& evtData,
					AliHLTComponentTriggerData& /*trigData*/) {
  // see header file for class documentation

  printf("AliHLTFlatAnaManagerComponent::DoEvent\n");
  Int_t iResult=0;

  // -- Only use data event
  if (!IsDataEvent()) {
    Printf("-------> no data event");
    return 0;
  }
  
  if( fUID == 0 ){
    TTimeStamp t;
    fUID = ( gSystem->GetPid() + t.GetNanoSec())*10 + evtData.fEventID;
  }
  
  Bool_t isESDfound = kFALSE;
  // -- Get ESD object
  // -------------------
  AliFlatESDEvent *esdEvent = NULL;
  AliESDfriend *esdFriend = NULL;
  Int_t dataBlockIndex = 0;
  for ( const AliHLTComponentBlockData *iter = GetFirstInputBlock(kAliHLTDataTypeFlatESD); iter != NULL; iter = GetNextInputBlock() ) {
    if(iter->fDataType != kAliHLTDataTypeFlatESD) {
      Printf("Data type of block is not flatESD");
      continue;
    }
    dataBlockIndex++;
    Printf("----> dataBlockIndex = %d", dataBlockIndex);
    //iter->Print();
    esdEvent = (AliFlatESDEvent*)(iter->fPtr);
    if( !esdEvent ){ 
      HLTWarning("Wrong ESDEvent object received");
      Printf("Current object is not ESD event");
      iResult = -1;
      continue;
    }
    else isESDfound = kTRUE;
    //    esdEvent->GetStdContent();
  }
  if (!isESDfound) {
    Printf("-------> no ESD found!!");
    return 0;
  }
  printf("----> ESDEvent %p has %d tracks: \n", esdEvent, esdEvent->GetNumberOfTracks());

  /* // no friends for teh time being... 
  for ( const TObject *iter = GetFirstInputObject(kAliHLTDataTypeESDfriendObject); iter != NULL; iter = GetNextInputObject() ) {
    esdFriend = dynamic_cast<AliESDfriend*>(const_cast<TObject*>( iter ) );
    if( !esdFriend ){ 
      HLTWarning("Wrong ESDFriend object received");
      iResult = -1;
      continue;
    }
  }
  printf("----> ESDFriend %p has %d tracks: \n", esdFriend, esdFriend->GetNumberOfTracks());
  */
  fAnalysisManager->InitInputData(esdEvent, esdFriend);
  //  fInputHandler->BeginEvent(0);
  fAnalysisManager->ExecAnalysis();
  fInputHandler->FinishEvent();


  // -- Send histlist
//  PushBack(dynamic_cast<TObject*>(fCorrObj->GetHistList()),
//	   kAliHLTDataTypeTObject|kAliHLTDataOriginHLT,fUID);
 
  return iResult;
}

// #################################################################################
Int_t AliHLTFlatAnaManagerComponent::Reconfigure(const Char_t* cdbEntry, const Char_t* chainId) {
  // see header file for class documentation

  Int_t iResult=0;
  TString cdbPath;
  if (cdbEntry) {
    cdbPath=cdbEntry;
  } else {
    cdbPath="HLT/ConfigGlobal/";
    cdbPath+=GetComponentID();
  }

  TString strcdb = "";
  if (!cdbEntry) strcdb = "default";
  Printf("reconfigure '%s' from entry %s%s", chainId, cdbPath.Data(), strcdb.Data());
  //AliInfoClass(Form("reconfigure '%s' from entry %s%s", chainId, cdbPath.Data(), cdbEntry?"":" (default)"));
  iResult=ConfigureFromCDBTObjString(cdbPath);

  return iResult;
}

// #################################################################################
Int_t AliHLTFlatAnaManagerComponent::ReadPreprocessorValues(const Char_t* /*modules*/) {
  // see header file for class documentation
  ALIHLTERRORGUARD(5, "ReadPreProcessorValues not implemented for this component");
  return 0;
}

