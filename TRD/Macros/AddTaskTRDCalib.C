//=============================================================================
//
// *** AddTaskTRDCalib
//
// This macros setup the TRD calibration task
//
//=============================================================================

const AliTRDCalDet *GetCalDet(Int_t runNumber);

AliAnalysisTask  *AddTaskTRDCalib(Int_t runNumber)
{
  gSystem->Load("libTRDcalib");
  // pointer to the analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskTRDCalib", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // check the input handler
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTask", "This task requires an input event handler");
    return NULL;
  }  

  /////////////////////////
  // The TRD calib Task
  /////////////////////////
  AliTRDCalibTask *calibTask = new AliTRDCalibTask();
  calibTask->SetHisto2d(kTRUE);
  calibTask->SetVector2d(kFALSE);
  calibTask->SetVdriftLinear(kTRUE);
  calibTask->SetNz(0,0);
  calibTask->SetNrphi(0,0);
  calibTask->SetNz(0,1);
  calibTask->SetNrphi(0,1);
  calibTask->SetNz(0,2);
  calibTask->SetNrphi(0,2);
  calibTask->SetLow(0);
  calibTask->SetHigh(30);
  calibTask->SetFillZero(kTRUE);
  calibTask->AddSelectedTriggerClass("CINT1B-ABCE-NOPF-ALL");
  calibTask->SetReject(kFALSE);
  //calibTask->SetDebug(3);
  calibTask->SetNbTimeBins(30);
  //calibTask->SetMaxEvent(20);
  calibTask->SetRequirePrimaryVertex(kTRUE);
  calibTask->SetMinNbOfContributors(1);
  calibTask->SetMaxCluster(100.0);
  calibTask->SetNbMaxCluster(2);
  calibTask->SetCalDetGain(GetCalDet(runNumber));

  /////////////////////////////
  // Track cuts
  /////////////////////////////
  AliESDtrackCuts *trackCuts = new AliESDtrackCuts("trackcuts","trackcuts");
  trackCuts->SetMinNClustersTPC(50);
  trackCuts->SetMaxChi2PerClusterTPC(3.5);
  //trackCuts->SetMaxCovDiagonalElements(2,2,0.5,0.5,2);
  trackCuts->SetRequireTPCRefit(kTRUE);
  //trackCuts->SetRequireITSRefit(kTRUE);
  //trackCuts->SetMinNsigmaToVertex(10);
  trackCuts->SetRequireSigmaToVertex(kFALSE);
  trackCuts->SetAcceptKinkDaughters(kFALSE);
  trackCuts->SetMaxDCAToVertexZ(30.0);
  trackCuts->SetMaxDCAToVertexXY(3.0);
  trackCuts->SetDCAToVertex2D(kFALSE);

  calibTask->SetESDtrackCuts(trackCuts);

  mgr->AddTask(calibTask);

  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  
  if (!cinput) cinput = mgr->CreateContainer("cchain",TChain::Class(), 
                                      AliAnalysisManager::kInputContainer);

  AliAnalysisDataContainer *coutput =mgr->CreateContainer("TRDCalib",TList::Class(), AliAnalysisManager::kOutputContainer, "AliESDfriends_v1.root");  


  mgr->ConnectInput(calibTask,0,cinput);
  mgr->ConnectOutput(calibTask,1,coutput);
  return calibTask;

}

const AliTRDCalDet *GetCalDet(Int_t runNumber){
  // 
  // Get Cal Det for gas gain
  //
  
  AliCDBEntry *entry = AliCDBManager::Instance()->Get("TRD/Calib/ChamberGainFactor",runNumber);
  if(!entry) {
    printf("Found no entry\n");
    return 0x0;
  }
  const AliTRDCalDet* calDet = (AliTRDCalDet *)entry->GetObject();

  return calDet;

}

