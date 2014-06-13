AliAnalysisTask *AddTaskTPCcalibResidualPID(TString period = "", Bool_t isPbpOrpPb = kFALSE,
                                            Bool_t producePIDqa = kTRUE, Bool_t useTPCCutMIGeo = kTRUE,
                                            Bool_t correctdEdxEtaDependence = kFALSE, 
                                            Bool_t correctdEdxMultiplicityDependence = kFALSE,
                                            Bool_t useMCinfo = kTRUE){
  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_statsQA_TPCresPID", "No analysis manager found.");
    return 0;
  }

  //========= Add task to the ANALYSIS manager =====
  AliTPCcalibResidualPID *task=new AliTPCcalibResidualPID("TPCresPID");
  task->SelectCollisionCandidates(AliVEvent::kMB | AliVEvent::kINT7);

  // which THnSparse should be produced
  task->SetProducePIDqa(producePIDqa);
  
  AliESDtrackCuts* esdTrackCuts = 0x0;
  AliESDtrackCuts* esdTrackCutsV0 = 0x0;
  
  if (period.Contains("LHC11") || period.Contains("LHC12") || period.Contains("LHC13")) {
    esdTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kTRUE);
    esdTrackCutsV0 = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE);
    task->SetESDtrackCuts(esdTrackCuts);
    //task->SetESDtrackCutsV0(esdTrackCutsV0);
    printf("Using standard ITS-TPC track cuts 2011.\n");
  }
  else if (period.Contains("LHC10")) {
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
  
  task->SetIsPbpOrpPb(isPbpOrpPb);
  if (task->GetIsPbpOrpPb()) {
    printf("Collision type pPb/Pbp set -> Adapting vertex cuts!\n");
  }
  else  {
    printf("Collision type different from pPb/Pbp -> Using standard vertex cuts!\n");
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
