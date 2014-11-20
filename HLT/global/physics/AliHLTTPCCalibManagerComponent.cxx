//-*- Mode: C++ -*-
// $Id: AliHLTTPCCalibManagerComponent.cxx $
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

/** @file    AliHLTTPCCalibManagerComponent.cxx
    @author  David Rohr, Jens Wiechula, C. Zampolli, I. Vorobyev
    @brief   Component for Testing TPC Calibration inside HLT component
*/

#include "TMap.h"
#include "TSystem.h"
#include "TTimeStamp.h"
#include "TObjString.h"
#include "TH1F.h"
#include "TList.h"
//#include "AliESDtrackCuts.h"
#include "AliESDEvent.h"
#include "AliHLTErrorGuard.h"
#include "AliHLTDataTypes.h"
#include "AliHLTTPCCalibManagerComponent.h"
#include "AliHLTITSClusterDataFormat.h"
#include "AliAnalysisManager.h"
#include "AliHLTVEventInputHandler.h"
#include "AliTPCAnalysisTaskcalib.h"
#include "AliAnalysisDataContainer.h"
#include "TTree.h"
#include "TChain.h"
#include "AliFlatESDEvent.h"
#include "AliFlatESDFriend.h"
#include "AliVEvent.h"
#include "AliVfriendEvent.h"

#include "AliTPCcalibAlign.h"
#include "AliTPCcalibBase.h"

#include "TROOT.h"

using namespace std;

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTTPCCalibManagerComponent)

/*
 * ---------------------------------------------------------------------------------
 *                            Constructor / Destructor
 * ---------------------------------------------------------------------------------
 */

// #################################################################################
AliHLTTPCCalibManagerComponent::AliHLTTPCCalibManagerComponent() :
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
AliHLTTPCCalibManagerComponent::~AliHLTTPCCalibManagerComponent() {
  // see header file for class documentation
}

/*
 * ---------------------------------------------------------------------------------
 * Public functions to implement AliHLTComponent's interface.
 * These functions are required for the registration process
 * ---------------------------------------------------------------------------------
 */

// #################################################################################
const Char_t* AliHLTTPCCalibManagerComponent::GetComponentID() {
  // see header file for class documentation
  return "TPCCalibManagerComponent";
}

// #################################################################################
void AliHLTTPCCalibManagerComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list) {
  // see header file for class documentation
  list.push_back(kAliHLTDataTypeESDObject|kAliHLTDataOriginAny);
  list.push_back(kAliHLTDataTypeClusters|kAliHLTDataOriginITSSPD);
  list.push_back(kAliHLTDataTypeESDContent|kAliHLTDataOriginVZERO);
  list.push_back(kAliHLTDataTypeESDfriendObject|kAliHLTDataOriginAny);
  list.push_back(kAliHLTDataTypeFlatESD|kAliHLTDataOriginOut);
  list.push_back(kAliHLTDataTypeFlatESDFriend|kAliHLTDataOriginOut);
}

// #################################################################################
AliHLTComponentDataType AliHLTTPCCalibManagerComponent::GetOutputDataType() {
  // see header file for class documentation
  return kAliHLTDataTypeTObject|kAliHLTDataOriginHLT;
}

// #################################################################################
void AliHLTTPCCalibManagerComponent::GetOutputDataSize( ULong_t& constBase, Double_t& inputMultiplier ) {
  // see header file for class documentation
  constBase = 100000;
  inputMultiplier = 0.5;
}

// #################################################################################
void AliHLTTPCCalibManagerComponent::GetOCDBObjectDescription( TMap* const targetMap) {
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
AliHLTComponent* AliHLTTPCCalibManagerComponent::Spawn() {
  // see header file for class documentation
  return new AliHLTTPCCalibManagerComponent;
}

/*
 * ---------------------------------------------------------------------------------
 * Protected functions to implement AliHLTComponent's interface.
 * These functions provide initialization as well as the actual processing
 * capabilities of the component. 
 * ---------------------------------------------------------------------------------
 */

// #################################################################################
Int_t AliHLTTPCCalibManagerComponent::DoInit( Int_t /*argc*/, const Char_t** /*argv*/ ) {
  // see header file for class documentation
  printf("AliHLTTPCCalibManagerComponent::DoInit\n");

  gROOT->Macro("$ALICE_ROOT/PWGPP/CalibMacros/CPass0/LoadLibraries.C");

  Int_t iResult=0;

  Printf("----> AliHLTTPCCalibManagerComponent::DoInit");
  fAnalysisManager = new AliAnalysisManager;
  fInputHandler    = new AliHLTVEventInputHandler;
  fAnalysisManager->SetInputEventHandler(fInputHandler);
  fAnalysisManager->SetExternalLoop(kTRUE);

  //gROOT->Macro("$ALICE_ROOT/PWGPP/CalibMacros/CPass0/AddTaskTPCCalib.C");

  // setup task TPCCalib
  AliTPCAnalysisTaskcalib *task1=new AliTPCAnalysisTaskcalib("CalibObjectsTrain1");

  AliTPCcalibAlign *calibAlign = new AliTPCcalibAlign("alignTPC","Alignment of the TPC sectors");
  calibAlign->SetDebugLevel(0);
  calibAlign->SetStreamLevel(0);
  calibAlign->SetTriggerMask(-1,-1,kTRUE);        //accept everything
  task1->AddJob(calibAlign);

  fAnalysisManager->AddTask(task1);
  AliAnalysisDataContainer *cinput1 = fAnalysisManager->GetCommonInputContainer();
  if (!cinput1) cinput1 = fAnalysisManager->CreateContainer("cchain",TChain::Class(),
                                      AliAnalysisManager::kInputContainer);
  for (Int_t i=0; i<task1->GetJobs()->GetEntries(); i++) {
    if (task1->GetJobs()->At(i)) {
      AliAnalysisDataContainer* coutput = fAnalysisManager->CreateContainer(task1->GetJobs()->At(i)->GetName(),
                                                               AliTPCcalibBase::Class(),
                                                               AliAnalysisManager::kOutputContainer,
                                                               "AliESDfriends_v1.root:TPCAlign");
      fAnalysisManager->ConnectOutput(task1,i,coutput);
    }
  }
  fAnalysisManager->ConnectInput(task1,0,cinput1);

  Printf("----> Calling InitAnalysis");
  fAnalysisManager->InitAnalysis();
  Printf("----> Done.");
  Printf("----> Calling StartAnalysis");
  fAnalysisManager->StartAnalysis("local", (TTree*)new TTree);
  Printf("----> Done.");

  return iResult;
}



// #################################################################################
Int_t AliHLTTPCCalibManagerComponent::DoDeinit() {
  // see header file for class documentation

  fUID = 0;
  fAnalysisManager->SetSkipTerminate(kTRUE);
  fAnalysisManager->Terminate();

  delete fAnalysisManager;

  return 0;
}

// #################################################################################
Int_t AliHLTTPCCalibManagerComponent::DoEvent(const AliHLTComponentEventData& evtData,
    AliHLTComponentTriggerData& /*trigData*/) {
  // see header file for class documentation

  printf("AliHLTTPCCalibManagerComponent::DoEvent\n");
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
  AliVEvent* vEvent=NULL;
  AliVfriendEvent* vFriend=NULL;

  for ( const TObject *iter = GetFirstInputObject(kAliHLTDataTypeESDObject); iter != NULL; iter = GetNextInputObject() ) {
    vEvent = dynamic_cast<AliESDEvent*>(const_cast<TObject*>( iter ) );
    if( !vEvent ){ 
      HLTWarning("Wrong ESDEvent object received");
      iResult = -1;
      continue;
    }
    vEvent->GetStdContent();
  }
  if (vEvent) {printf("----> ESDEvent %p has %d tracks: \n", vEvent, vEvent->GetNumberOfTracks());}
  for ( const TObject *iter = GetFirstInputObject(kAliHLTDataTypeESDfriendObject); iter != NULL; iter = GetNextInputObject() ) {
    vFriend = dynamic_cast<AliESDfriend*>(const_cast<TObject*>( iter ) );
    if( !vFriend ){ 
      HLTWarning("Wrong ESDFriend object received");
      iResult = -1;
      continue;
    }
  }
  if(vFriend) {printf("----> ESDFriend %p has %d tracks: \n", vFriend, vFriend->GetNumberOfTracks());}

  if (!vEvent){
    for (const AliHLTComponentBlockData* pBlock=GetFirstInputBlock(kAliHLTDataTypeFlatESD|kAliHLTDataOriginOut);
        pBlock!=NULL; pBlock=GetNextInputBlock()) {
      AliFlatESDEvent* tmpFlatEvent=reinterpret_cast<AliFlatESDEvent*>( pBlock->fPtr );
      if (tmpFlatEvent->GetSize()==pBlock->fSize ){
        tmpFlatEvent->Reinitialize();
      } else {
        tmpFlatEvent = NULL;
        HLTWarning("data mismatch in block %s (0x%08x): size %d -> ignoring flatESD information",
            DataType2Text(pBlock->fDataType).c_str(), pBlock->fSpecification, pBlock->fSize);
      }
      vEvent = tmpFlatEvent;
      break;
    }

    if( vEvent ){
      for (const AliHLTComponentBlockData* pBlock=GetFirstInputBlock(kAliHLTDataTypeFlatESDFriend|kAliHLTDataOriginOut);
          pBlock!=NULL; pBlock=GetNextInputBlock()) {
        AliFlatESDFriend* tmpFlatFriend = reinterpret_cast<AliFlatESDFriend*>( pBlock->fPtr );
        if (tmpFlatFriend->GetSize()==pBlock->fSize ){
          tmpFlatFriend->Reinitialize();
        } else {
          tmpFlatFriend = NULL;
          HLTWarning("data mismatch in block %s (0x%08x): size %d -> ignoring flatESDfriend information", 
              DataType2Text(pBlock->fDataType).c_str(), pBlock->fSpecification, pBlock->fSize);
        }
        break;
        vFriend = tmpFlatFriend;
      }   
    }
  }

  fAnalysisManager->InitInputData(vEvent, vFriend);
  Printf("The friend is: %p",vEvent->FindFriend());


  //  fInputHandler->BeginEvent(0);
  Printf("----> ExecAnalysis");
  fAnalysisManager->ExecAnalysis();
  fInputHandler->FinishEvent();


  // -- Send histlist
  //  PushBack(dynamic_cast<TObject*>(fCorrObj->GetHistList()),
  //	   kAliHLTDataTypeTObject|kAliHLTDataOriginHLT,fUID);

  return iResult;
}

// #################################################################################
Int_t AliHLTTPCCalibManagerComponent::Reconfigure(const Char_t* cdbEntry, const Char_t* chainId) {
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
Int_t AliHLTTPCCalibManagerComponent::ReadPreprocessorValues(const Char_t* /*modules*/) {
  // see header file for class documentation
  ALIHLTERRORGUARD(5, "ReadPreProcessorValues not implemented for this component");
  return 0;
}

