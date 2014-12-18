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
#include "AliHLTAnalysisManager.h"
#include "AliHLTVEventInputHandler.h"
#include "AliTPCAnalysisTaskcalib.h"
#include "AliAnalysisDataContainer.h"
#include "TTree.h"
#include "TChain.h"
#include "AliFlatESDEvent.h"
#include "AliFlatESDFriend.h"
#include "AliVEvent.h"
#include "AliVfriendEvent.h"
#include "AliSysInfo.h"

#include "AliTPCcalibBase.h"

#include "AliTPCcalibAlign.h"
#include "AliTPCcalibLaser.h"
#include "AliTPCcalibCosmic.h"

#include "AliTPCcalibCalib.h"
#include "AliTPCcalibTimeGain.h"
#include "AliTPCcalibGainMult.h"
#include "AliTPCcalibTime.h"

#include "AliTPCcalibTracks.h"

#include "AliTPCParam.h"
#include "AliMagF.h"
#include "TGeoGlobalMagField.h"
#include "AliTPCTransform.h"
#include "AliCDBEntry.h"
#include "AliTPCRecoParam.h"
#include "AliCDBManager.h"
#include "AliGRPObject.h"
#include "AliTPCcalibDB.h"
#include "AliTPCClusterParam.h"
#include "AliTPCcalibTracksCuts.h"

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
  HLTInfo("AliHLTTPCCalibManagerComponent::DoInit\n");

  //shut up AliSysInfo globally, prevents syswatch.log files from being created
  //dont open any debug files
  AliSysInfo::SetDisabled(kTRUE);
  TTreeSRedirector::SetDisabled(kTRUE);

  HLTInfo("----> AliHLTTPCCalibManagerComponent::DoInit");
  fAnalysisManager = new AliHLTAnalysisManager();
  fInputHandler    = new AliHLTVEventInputHandler("HLTinputHandler","HLT input handler");
  fAnalysisManager->SetInputEventHandler(fInputHandler);

  AddCalibTasks();

  fAnalysisManager->InitAnalysis();

  return 0;
}

// #################################################################################
Int_t AliHLTTPCCalibManagerComponent::DoDeinit() {
  // see header file for class documentation

  fUID = 0;
  fAnalysisManager->WriteAnalysisToFile();
  delete fAnalysisManager;
  return 0;
}

// #################################################################################
Int_t AliHLTTPCCalibManagerComponent::DoEvent(const AliHLTComponentEventData& evtData,
    AliHLTComponentTriggerData& /*trigData*/) {
  // see header file for class documentation

  HLTInfo("AliHLTTPCCalibManagerComponent::DoEvent\n");
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

  ReadInput(vEvent,vFriend);

  if (vEvent) {HLTInfo("----> event %p has %d tracks: \n", vEvent, vEvent->GetNumberOfTracks());}
  if(vFriend) {HLTInfo("----> friend %p has %d tracks: \n", vFriend, vFriend->GetNumberOfTracks());}

  fInputHandler->InitTaskInputData(vEvent, vFriend, fAnalysisManager->GetTasks());
  fAnalysisManager->ExecAnalysis();
  fInputHandler->FinishEvent();

  //PushBack(dynamic_cast<TObject*>(fAnalysisManager->GetOutputs()), 
  //    kAliHLTDataTypeTObject|kAliHLTDataOriginHLT,fUID);

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

// #################################################################################
Int_t AliHLTTPCCalibManagerComponent::AddCalibTasks()
{
  //add tasks to the manager
  // Set run time ranges (time stamps)
  // doesn't really work at the moment, to be checked

  gROOT->Macro("$ALICE_ROOT/PWGPP/CalibMacros/CPass0/LoadLibraries.C");

  AliCDBEntry* entry = AliCDBManager::Instance()->Get("GRP/GRP/Data");
  if(!entry) {
    ::Error("AddCalibTimeGain","Cannot get AliCDBEntry");
    //return;
  }
  const AliGRPObject* grpData = dynamic_cast<AliGRPObject*>(entry->GetObject());
  if(!grpData) {
    ::Error("AddCalibTimeGain","Cannot get AliGRPObject");
    //return;
  }

  time_t sTime = grpData->GetTimeStart();
  time_t eTime = grpData->GetTimeEnd();
  TTimeStamp startRunTime(sTime);
  TTimeStamp stopRunTime(eTime);
  UInt_t year;
  startRunTime.GetDate(kTRUE,0,&year);
  TTimeStamp startTime(year,1,1,0,0,0);
  TTimeStamp stopTime(year,12,31,23,59,59);
  Bool_t useQmax = (grpData->GetBeamType()).Contains("Pb-Pb");

  // setup task TPCCalib-------------------------------------------------------------
  AliTPCAnalysisTaskcalib *task1=new AliTPCAnalysisTaskcalib("CalibObjectsTrain1");

  AliTPCcalibCalib *calibCalib = new AliTPCcalibCalib("calibTPC","calibTPC");
  calibCalib->SetDebugLevel(0);
  calibCalib->SetStreamLevel(0);
  calibCalib->SetTriggerMask(-1,-1,kFALSE);        //accept everything
  task1->AddJob(calibCalib);

  AliTPCcalibTimeGain *calibTimeGain = new AliTPCcalibTimeGain("calibTimeGain","calibTimeGain", startTime.GetSec(), stopTime.GetSec(), 10*60);
  calibTimeGain->SetIsCosmic(kFALSE);
  calibTimeGain->SetUseCookAnalytical(kTRUE);
  calibTimeGain->SetUseMax(useQmax);
  calibTimeGain->SetDebugLevel(0);
  calibTimeGain->SetStreamLevel(0);
  calibTimeGain->SetTriggerMask(-1,-1,kTRUE);        //reject laser
  calibTimeGain->SetLowerTrunc(0.02);
  calibTimeGain->SetUpperTrunc(0.6);
  task1->AddJob(calibTimeGain);

  AliTPCcalibGainMult *calibGainMult = new AliTPCcalibGainMult("calibGainMult","calibGainMult");
  calibGainMult->SetUseMax(useQmax);
  calibGainMult->SetDebugLevel(0);
  calibGainMult->SetStreamLevel(0);
  calibGainMult->SetTriggerMask(-1,-1,kTRUE);        //reject laser
  calibGainMult->SetLowerTrunc(0.02);
  calibGainMult->SetUpperTrunc(0.6);
  task1->AddJob(calibGainMult);

  //caliTime crashes at the moment, more investigation is needed
  AliTPCcalibTime *calibTime = new AliTPCcalibTime("calibTime","calibTime",  startTime.GetSec(), stopTime.GetSec(), 10*60, 2);
  calibTime->SetDebugLevel(0);
  calibTime->SetStreamLevel(0);
  calibTime->SetTriggerMask(-1,-1,kFALSE);        //accept everything
  calibTime->SetCutTracks(15000);              // max 15000 tracks per event
  task1->AddJob(calibTime);

  fAnalysisManager->AddTask(task1);
  //AliAnalysisDataContainer *cinput1 = fAnalysisManager->GetCommonInputContainer();
  //if (!cinput1) cinput1 = fAnalysisManager->CreateContainer("cchain",TChain::Class(),
  //                                    AliAnalysisManager::kInputContainer);
  for (Int_t i=0; i<task1->GetJobs()->GetEntries(); i++) {
    if (task1->GetJobs()->At(i)) {
      AliAnalysisDataContainer* coutput = fAnalysisManager->CreateContainer(task1->GetJobs()->At(i)->GetName(),
                                                               AliTPCcalibBase::Class(),
                                                               AliAnalysisManager::kOutputContainer,
                                                               "AliESDfriends_v1.root:TPCCalib");
      fAnalysisManager->ConnectOutput(task1,i,coutput);
    }
  }
  //fAnalysisManager->ConnectInput(task1,0,cinput1);

  // setup task TPCAlign--------------------------------------------------------------------------
  //AliTPCAnalysisTaskcalib *taskAlign=new AliTPCAnalysisTaskcalib("CalibObjectsTrain1");

  //AliTPCcalibAlign *calibAlign = new AliTPCcalibAlign("alignTPC","Alignment of the TPC sectors"); 
  //calibAlign->SetDebugLevel(0);
  //calibAlign->SetStreamLevel(0);
  //calibAlign->SetTriggerMask(-1,-1,kTRUE);        //accept everything
  //taskAlign->AddJob(calibAlign);

  //Laser crashes, more investigation is needed
  /*AliTPCcalibLaser *calibLaser = new AliTPCcalibLaser("laserTPC","laserTPC");
  calibLaser->SetDebugLevel(0);
  calibLaser->SetStreamLevel(0);
  calibLaser->SetTriggerMask(-1,-1,kFALSE);        //accept everything
  taskAlign->AddJob(calibLaser);*/

  //AliTPCcalibCosmic *calibCosmic = new AliTPCcalibCosmic("cosmicTPC","cosmicTPC");
  //calibCosmic->SetDebugLevel(0);
  //calibCosmic->SetStreamLevel(1);
  //calibCosmic->SetTriggerMask(-1,-1,kTRUE);        //accept everything
  //taskAlign->AddJob(calibCosmic);

  //fAnalysisManager->AddTask(taskAlign);
  //for (Int_t i=0; i<taskAlign->GetJobs()->GetEntries(); i++) {
  //  if (taskAlign->GetJobs()->At(i)) {
  //    AliAnalysisDataContainer* coutput = fAnalysisManager->CreateContainer(taskAlign->GetJobs()->At(i)->GetName(),
  //                                                             AliTPCcalibBase::Class(),
  //                                                             AliAnalysisManager::kOutputContainer,
  //                                                             "AliESDfriends_v1.root:TPCAlign");
  //    fAnalysisManager->ConnectOutput(taskAlign,i,coutput);
  //  }
  //}
  //fAnalysisManager->ConnectInput(taskAlign,0,cinput1);

  // setup task TPCCluster----------------------------------------------------------------
  //AliTPCAnalysisTaskcalib *taskCluster=new AliTPCAnalysisTaskcalib("CalibObjectsTrain1");

  //AliTPCClusterParam * clusterParam = AliTPCcalibDB::Instance()->GetClusterParam();
  //AliTPCcalibTracksCuts *cuts = new AliTPCcalibTracksCuts(30, 0.4, 5, 0.13, 0.018);
  //
  //AliTPCcalibTracks *calibTracks =  new AliTPCcalibTracks("calibTracks", "Resolution calibration object for tracks", clusterParam, cuts);
  //calibTracks->SetDebugLevel(0);
  //calibTracks->SetStreamLevel(0);
  //calibTracks->SetTriggerMask(-1,-1,kTRUE);
  //taskCluster->AddJob(calibTracks);

  //fAnalysisManager->AddTask(taskCluster);
  //for (Int_t i=0; i<taskCluster->GetJobs()->GetEntries(); i++) {
  //  if (taskCluster->GetJobs()->At(i)) {
  //    AliAnalysisDataContainer* coutput = fAnalysisManager->CreateContainer(taskCluster->GetJobs()->At(i)->GetName(),
  //                                                             AliTPCcalibBase::Class(),
  //                                                             AliAnalysisManager::kOutputContainer,
  //                                                             "AliESDfriends_v1.root:TPCCluster");
  //    fAnalysisManager->ConnectOutput(taskCluster,i,coutput);
  //  }
  //}
  //fAnalysisManager->ConnectInput(taskCluster,0,cinput1);
  //--------------------------------------------------------------------------------------
  return 0;
}

// #################################################################################
Int_t AliHLTTPCCalibManagerComponent::ReadInput(AliVEvent*& vEvent, AliVfriendEvent*& vFriend)
{
  Int_t iResult=0;
  for ( const TObject *iter = GetFirstInputObject(kAliHLTDataTypeESDObject); iter != NULL; iter = GetNextInputObject() ) {
    vEvent = dynamic_cast<AliESDEvent*>(const_cast<TObject*>( iter ) );
    if( !vEvent ){ 
      HLTWarning("Wrong ESDEvent object received");
      iResult = -1;
      continue;
    }
    vEvent->GetStdContent();
  }
  for ( const TObject *iter = GetFirstInputObject(kAliHLTDataTypeESDfriendObject); iter != NULL; iter = GetNextInputObject() ) {
    vFriend = dynamic_cast<AliESDfriend*>(const_cast<TObject*>( iter ) );
    if( !vFriend ){ 
      HLTWarning("Wrong ESDFriend object received");
      iResult = -1;
      continue;
    }
  }

  if (!vEvent){ 
    {
      const AliHLTComponentBlockData* pBlock=GetFirstInputBlock(kAliHLTDataTypeFlatESD|kAliHLTDataOriginOut);
      AliFlatESDEvent* tmpFlatEvent=reinterpret_cast<AliFlatESDEvent*>( pBlock->fPtr );
      if (tmpFlatEvent->GetSize()==pBlock->fSize ){
        tmpFlatEvent->Reinitialize();
      } else {
        tmpFlatEvent = NULL;
        HLTWarning("data mismatch in block %s (0x%08x): size %d -> ignoring flatESD information",
            DataType2Text(pBlock->fDataType).c_str(), pBlock->fSpecification, pBlock->fSize);
      }
      vEvent = tmpFlatEvent;
    }

    if( vEvent ){
      const AliHLTComponentBlockData* pBlock=GetFirstInputBlock(kAliHLTDataTypeFlatESDFriend|kAliHLTDataOriginOut);
      AliFlatESDFriend* tmpFlatFriend = reinterpret_cast<AliFlatESDFriend*>( pBlock->fPtr );
      if (tmpFlatFriend->GetSize()==pBlock->fSize ){
        tmpFlatFriend->Reinitialize();
      } else {
        tmpFlatFriend = NULL;
        HLTWarning("data mismatch in block %s (0x%08x): size %d -> ignoring flatESDfriend information", 
            DataType2Text(pBlock->fDataType).c_str(), pBlock->fSpecification, pBlock->fSize);
      }
      vFriend = tmpFlatFriend;
    }
  }
  return iResult;
}
