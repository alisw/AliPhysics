AliTOFAnalysisTaskCalibTree *AddTOFAnalysisTaskCalibTree(Bool_t savecoords = kFALSE) {

  /* check analysis manager */
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddAnalysisTaskEventTime", "cannot get analysis manager");
    return NULL;
  }

  /* check input event handler */
  if (!mgr->GetInputEventHandler()) {
    Error("AddAnalysisTaskEventTime", "cannot get input event handler");
    return NULL;
  }
  
  /* check input data type */
  TString str = mgr->GetInputEventHandler()->GetDataType();
  if (str.CompareTo("ESD")) {
    Error("AddAnalysisTaskEventTime", "input data type is not \"ESD\"");
    return NULL;
  }

  /* get common input data container */
  AliAnalysisDataContainer *inputc = mgr->GetCommonInputContainer();
  if (!inputc) {
    Error("AddAnalysisTaskEventTime", "cannot get common input container");
    return NULL;
  }
  
  /* create output data container */
  // setup output event handler
  AliAnalysisDataContainer *coutput   = mgr->CreateContainer(Form("aodTree"), TTree::Class(), AliAnalysisManager::kOutputContainer, "TOFcalibTree.root"); // tree
  if (!coutput) {
    Error("AddTOFAnalysisTaskCalibTree", "cannot create output container");
    return NULL;
  }

  /*  create task and connect input/output */
  AliTOFAnalysisTaskCalibTree *task = new AliTOFAnalysisTaskCalibTree();
  Printf("After initializing the TOF task: task = %p", task);
  // adding the task
  mgr->AddTask(task);

  mgr->ConnectInput(task, 0, inputc);
  mgr->ConnectOutput(task, 1, coutput);

  // setup task 
  task->SetEventSelectionFlag(kFALSE);
  task->SetVertexSelectionFlag(kTRUE);
  task->SetVertexCut(50.0);
  task->SetDiscardPileupEventFlag(kFALSE);
  task->SetPrimaryDCASelectionFlag(kFALSE);
  task->SetCalibrateTOFsignal(kFALSE);
  task->SetComputeT0TOF(kTRUE);
  task->SetUseT0TOF(kFALSE);
  task->SetUseLHCClockPhase(kFALSE);
  task->SetSaveCoordinates(savecoords);
 

  //  task->SetSpecificStorageParOffline("alien://?folder=/alice/cern.ch/user/r/rpreghen/OCDB");
  //  task->SetSpecificStorageRunParams("alien://?folder=/alice/cern.ch/user/r/rpreghen/OCDB");

  // setup event cuts 
  task->GetEventCuts()->SetAnalyzeMC(kFALSE);

  // setup TOF calib 
  task->GetTOFcalib()->SetRemoveMeanT0(kFALSE);
  task->GetTOFcalib()->SetCalibrateTOFsignal(kTRUE);
  task->GetTOFcalib()->SetCorrectTExp(kFALSE);

  //setup resolution 
  Double_t timeReso = 100.;

  // setup TOF response 
  //task->GetESDpid()->GetTOFResponse().SetTimeResolution(timeReso);

  // setup TOF-T0 maker 
  task->GetTOFT0maker()->SetTimeResolution(timeReso);

  // setup track cuts 
  AliESDtrackCuts *trackCuts = task->GetTrackCuts();
  trackCuts->SetPtRange(0.15, 10.);
  trackCuts->SetEtaRange(-1.0, 1.0);
  trackCuts->SetRequireITSRefit(kTRUE);
  trackCuts->SetMinNClustersITS(1);
  //  trackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);
  trackCuts->SetRequireTPCRefit(kTRUE);
  trackCuts->SetMinNClustersTPC(70);
  trackCuts->SetMaxChi2PerClusterTPC(4.);
  trackCuts->SetAcceptKinkDaughters(kFALSE);
  trackCuts->SetMaxDCAToVertexZ(3.2);
  trackCuts->SetMaxDCAToVertexXY(2.4);
  trackCuts->SetDCAToVertex2D(kTRUE);

  /* return task */
  return task;

}
