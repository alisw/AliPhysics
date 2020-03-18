AliAnalysisTaskSE *AddTaskTPCcalibResidualPID_Modified(TString period = "", Bool_t isPbpOrpPb = kFALSE,
                                            Bool_t producePIDqa = kTRUE, Bool_t useTPCCutMIGeo = kFALSE,
                                            Bool_t writeAdditionalQA = kFALSE,
                                            Bool_t cutOnProdRadiusForV0el = kTRUE,
                                            Bool_t correctdEdxEtaDependence = kFALSE,
                                            Bool_t correctdEdxMultiplicityDependence = kFALSE,
                                            Bool_t useMCinfo = kTRUE,
                                            Bool_t runInQA = kFALSE
                                           ){
  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_statsQA_TPCresPID", "No analysis manager found.");
    return 0;
  }

  //========= Add task to the ANALYSIS manager =====
  AliTPCcalibResidualPID_Modified *task=new AliTPCcalibResidualPID_Modified("TPCresPID");
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

  task->SetWriteAdditionalOutput(writeAdditionalQA);

  task->SetCutOnProdRadiusForV0el(cutOnProdRadiusForV0el);

  task->SetCorrectdEdxEtaDependence(correctdEdxEtaDependence);

  task->SetCorrectdEdxMultiplicityDependence(correctdEdxMultiplicityDependence);

  task->SetUseMCinfo(useMCinfo);

  printf("Eta correction: %s for this task\n",
         task->GetCorrectdEdxEtaDependence() ? "enabled (only works if enabled in PIDresponse!)" : "explicitly disabled");

  printf("Multiplicity correction: %s for this task\n",
         task->GetCorrectdEdxMultiplicityDependence() ? "enabled (only works if enabled in PIDresponse!)" : "explicitly disabled");

  printf("Use MC info: %d\n", task->GetUseMCinfo());

  printf("WriteAdditionalOutput: %d\n", task->GetWriteAdditionalOutput());

  printf("CutOnProdRadiusForV0el: %d\n", task->GetCutOnProdRadiusForV0el());


  mgr->AddTask(task);

  //================================================
  //              data containers
  //================================================
  //            find input container


  //            define output containers
  TString outputFile = "TPCresidualPID.root";
  if (runInQA) outputFile = Form("%s:TPCresPID", AliAnalysisManager::GetCommonFileName());
  AliAnalysisDataContainer *coutput1 =
    mgr->CreateContainer("TPCresPID", TObjArray::Class(),
                         AliAnalysisManager::kOutputContainer,outputFile);

  //           connect containers
  mgr->ConnectInput  (task,  0, mgr->GetCommonInputContainer() );
  mgr->ConnectOutput (task,  1, coutput1);

  if (task->GetWriteAdditionalOutput()) {
    AliAnalysisDataContainer *coutput2 =
      mgr->CreateContainer("TPCresPIDQA", TObjArray::Class(),
                          AliAnalysisManager::kOutputContainer,"TPCresidualPIDQA.root");

    AliAnalysisDataContainer *coutput3 =
      mgr->CreateContainer("TPCresPIDTreeEl", TTree::Class(),
                          AliAnalysisManager::kOutputContainer,"TPCresidualPIDTree_el.root");

    AliAnalysisDataContainer *coutput4 =
      mgr->CreateContainer("TPCresPIDTreePi", TTree::Class(),
                          AliAnalysisManager::kOutputContainer,"TPCresidualPIDTree_pi.root");

    AliAnalysisDataContainer *coutput5 =
      mgr->CreateContainer("TPCresPIDTreePr", TTree::Class(),
                          AliAnalysisManager::kOutputContainer,"TPCresidualPIDTree_pr.root");


    mgr->ConnectOutput (task,  2, coutput2);
    mgr->ConnectOutput (task,  3, coutput3);
    mgr->ConnectOutput (task,  4, coutput4);
    mgr->ConnectOutput (task,  5, coutput5);

  }

  //Store Reduced Tree (Alberto Caliva)
  AliAnalysisDataContainer *coutput6 =
  mgr->CreateContainer("ReducedTree_TPCSplines", TTree::Class(),
                      AliAnalysisManager::kOutputContainer,"TPCsplined_ReducedTree.root");

  mgr->ConnectOutput (task,  6, coutput6);//(Alberto Caliva)

  return task;
}
