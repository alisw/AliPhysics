AliAnalysisTask *AddTaskTPCPIDEtaQA(TString period = "", Bool_t isPbpOrpPb = kFALSE,
                                    Int_t tpcCutType = AliTPCPIDBase::kTPCCutMIGeo /*AliTPCPIDBase::kTPCnclCut*/,
                                    Bool_t usePhiCut = kFALSE,
                                    Double_t ptThresholdForPhiCut = 0.0,
                                    TString centralityEstimator = ""/*"ITSTPCtracklets" or "ppMultV0M" or ""*/){
  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_bhess_PIDetaAdv", "No analysis manager found.");
    return 0;
  }
  
  //========= Add task to the ANALYSIS manager =====
  AliTPCPIDEtaQA *task = new AliTPCPIDEtaQA("TPCPIDEtaQA");
  task->SelectCollisionCandidates(AliVEvent::kMB | AliVEvent::kINT7);
  
  //
  // Add track filters
  //
  AliAnalysisFilter* trackFilter = new AliAnalysisFilter("trackFilter");
  AliESDtrackCuts* esdTrackCutsL = 0x0;
  
  printf("\nSettings:\n");
  if (period.Contains("LHC11") || period.Contains("LHC12") || period.Contains("LHC13")) {
    esdTrackCutsL = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kTRUE);
    printf("Using standard ITS-TPC track cuts 2011\n");
  }
  else if (period.Contains("LHC10")) {
    esdTrackCutsL = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kTRUE);
    printf("Using standard ITS-TPC track cuts 2010\n");
  }
  else  {
    esdTrackCutsL = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kTRUE);
    printf("WARNING: Cuts not configured for this period!!! Using standard ITS-TPC track cuts 2011\n");
  }
  
/*
  esdTrackCutsL->SetMinNCrossedRowsTPC(120);
  esdTrackCutsL->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
  esdTrackCutsL->SetMaxChi2PerClusterITS(36);
  esdTrackCutsL->SetMaxFractionSharedTPCClusters(0.4);
  esdTrackCutsL->SetMaxChi2TPCConstrainedGlobal(36);
*/  

  task->SetIsPbpOrpPb(isPbpOrpPb);
  if (task->GetIsPbpOrpPb()) {
    printf("Collision type pPb/Pbp set -> Adapting vertex cuts!\n");
  }
  else  {
    printf("Collision type different from pPb/Pbp -> Using standard vertex cuts!\n");
  }
  
  trackFilter->AddCuts(esdTrackCutsL);
  task->SetTrackFilter(trackFilter);
  
  task->SetEtaCut(0.9);
  task->SetUsePhiCut(usePhiCut);
  task->SetPtThresholdForPhiCut(ptThresholdForPhiCut);
  task->SetTPCcutType(tpcCutType);
  task->SetCentralityEstimator(centralityEstimator);
  
  printf("Eta cut: %f\n", task->GetEtaCut());
  printf("UsePhiCut: %d\n", task->GetUsePhiCut());
  if (task->GetUsePhiCut())
    printf("PtThresholdForPhiCut: %f\n", task->GetPtThresholdForPhiCut());
  printf("UseTPCCutMIGeo: %d\n", task->GetUseTPCCutMIGeo());
  printf("UseTPCnclCut: %d\n", task->GetUseTPCnclCut());
  printf("Centrality estimator: \"%s\"\n", task->GetCentralityEstimator().Data());
  
  
  
  
  
  
  task->SetZvtxCutEvent(10.0);
  printf("Cut on z position of vertex: %.2f cm\n", task->GetZvtxCutEvent());
  
  printf("UsePhiCut: %d\nPtThresholdForPhiCut: %.3f GeV/c\n\n", task->GetUsePhiCut(), task->GetPtThresholdForPhiCut());
  mgr->AddTask(task);


  //================================================
  //              data containers
  //================================================
  //            find input container
  //below the trunk version
  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();

  //dumm output container
  AliAnalysisDataContainer *coutput0 =
      mgr->CreateContainer("TPCPIDEtaQA_tree",
                           TTree::Class(),
                           AliAnalysisManager::kExchangeContainer,
                           "TPCPIDEtaQA_default");

  //define output containers, please use 'username'_'somename'
  AliAnalysisDataContainer *coutput1 = 
      mgr->CreateContainer("TPCPIDEtaQA", TObjArray::Class(),
                           AliAnalysisManager::kOutputContainer,"TPCPIDEtaQA.root");

  //connect containers
  mgr->ConnectInput  (task,  0, cinput );
  mgr->ConnectOutput (task,  0, coutput0);
  mgr->ConnectOutput (task,  1, coutput1);

  return task;
}
