AliAnalysisTask *AddTaskTPCcalibResidualPID(Bool_t producePIDqa = kTRUE, Bool_t useTPCCutMIGeo = kTRUE,
                                            Bool_t correctdEdxEtaDependence = kFALSE, 
                                            Bool_t correctdEdxMultiplicityDependence = kFALSE,
                                            Bool_t useMCinfo = kTRUE){
  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_statsQA_TPCresPID", "No analysis manager found.");
    return 0;
  }

  TString trainConfig=gSystem->Getenv("CONFIG_FILE");

  TString list=gSystem->Getenv("LIST");
  Bool_t isLHC11h = list.Contains("LHC11h");


  //========= Add task to the ANALYSIS manager =====
  AliTPCcalibResidualPID *task=new AliTPCcalibResidualPID("TPCresPID");
  task->SelectCollisionCandidates(AliVEvent::kMB | AliVEvent::kINT7);

  // which THnSparse should be produced
  task->SetProducePIDqa(producePIDqa);
  
  AliESDtrackCuts* esdTrackCuts = 0x0;
  AliESDtrackCuts* esdTrackCutsV0 = 0x0;
  
  //TODO LHC12 + check if cuts are correct in this way!
  TString listOfFiles = gSystem->Getenv("LIST");
  if (listOfFiles.Contains("LHC11") || listOfFiles.Contains("LHC12") || listOfFiles.Contains("LHC13")) {
    esdTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kTRUE);
    esdTrackCutsV0 = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE);
    task->SetESDtrackCuts(esdTrackCuts);
    //task->SetESDtrackCutsV0(esdTrackCutsV0);
    printf("Using standard ITS-TPC track cuts 2011.\n");
  }
  else if (listOfFiles.Contains("LHC10")) {
    esdTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kTRUE);
    esdTrackCutsV0 = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kFALSE);
    task->SetESDtrackCuts(esdTrackCuts);
    //task->SetESDtrackCutsV0(esdTrackCutsV0);
    printf("Using standard ITS-TPC track cuts 2010.\n");
  }
  else  {
    esdTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kTRUE);
    esdTrackCutsV0 = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE);
    task->SetESDtrackCuts(esdTrackCuts);
    //task->SetESDtrackCutsV0(esdTrackCutsV0);
    
    printf("WARNING: Cuts not configured for this period!!! Using standard ITS-TPC track cuts 2011\n");
  }
  
  // Test whether we have pPb or Pbp
  if (listOfFiles.Contains("pPb") || listOfFiles.Contains("Pbp")) {
    task->SetIsPbpOrpPb(kTRUE);
    printf("pPb/Pbp detected -> Adapting vertex cuts!\n");
  }
  else  {
    task->SetIsPbpOrpPb(kFALSE);
    printf("Collision type different from pPb/Pbp detected -> Using standard vertex cuts!\n");
  }

  task->SetZvtxCutEvent(10.0);
  printf("Cut on z position of vertex: %.2f cm\n", task->GetZvtxCutEvent());
  
  task->SetUseTPCCutMIGeo(useTPCCutMIGeo);
  printf("UseTPCCutMIGeo: %d\n", task->GetUseTPCCutMIGeo());
  
  task->SetCorrectdEdxEtaDependence(correctdEdxEtaDependence);
  
  task->SetCorrectdEdxMultiplicityDependence(correctdEdxMultiplicityDependence);
  
  task->SetUseMCinfo(useMCinfo);
  
  printf("Eta correction: %s for this task\n", 
         task->GetCorrectdEdxEtaDependence() ? "enabled (only works if enabled in PIDresponse!)" : "explicitly disabled");
  
  printf("Multiplicity correction: %s for this task\n", 
         task->GetCorrectdEdxMultiplicityDependence() ? "enabled (only works if enabled in PIDresponse!)" : "explicitly disabled");

  printf("Use MC info: %d\n", task->GetUseMCinfo());
  
  
  mgr->AddTask(task);

  //================================================
  //              data containers
  //================================================
  //            find input container
  
  
  //            define output containers, please use 'username'_'somename'
  AliAnalysisDataContainer *coutput1 =
    mgr->CreateContainer("TPCresPID", TObjArray::Class(),
                         AliAnalysisManager::kOutputContainer,"TPCresidualPID.root");
    
  AliAnalysisDataContainer *coutput2 =
    mgr->CreateContainer("TPCresPIDQA", TObjArray::Class(),
                         AliAnalysisManager::kOutputContainer,"TPCresidualPIDQA.root");
  
  //           connect containers
  mgr->ConnectInput  (task,  0, mgr->GetCommonInputContainer() );
  mgr->ConnectOutput (task,  1, coutput1);
  mgr->ConnectOutput (task,  2, coutput2);
  
  return task;
}
