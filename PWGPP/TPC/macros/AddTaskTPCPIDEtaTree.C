AliAnalysisTask *AddTaskTPCPIDEtaTree(TString period = "", Bool_t isPbpOrpPb = kFALSE, Bool_t storeMultiplicity = kTRUE,
                                      Bool_t correctdEdxEtaDependence = kFALSE, Bool_t correctdEdxMultiplicityDependence = kFALSE,
                                      Bool_t setDoAdditionalQA = kFALSE,
                                      AliTPCPIDBase::TPCcutType tpcCutType = AliTPCPIDBase::kTPCCutMIGeo /*AliTPCPIDBase::kTPCnclCut*/,
                                      Bool_t usePhiCut = kFALSE,
                                      Bool_t useFilteredTreeCuts = kFALSE){
  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskTPCPIDetaTree", "No analysis manager found.");
    return 0;
  }
  
  //========= Add task to the ANALYSIS manager =====
  AliTPCPIDEtaTree *task = new AliTPCPIDEtaTree("TPCPIDEtaTree");
  task->SelectCollisionCandidates(AliVEvent::kMB | AliVEvent::kINT7 | AliVEvent::kINT8);
  
  //
  // Add track filters
  //
  AliAnalysisFilter* trackFilter = new AliAnalysisFilter("trackFilter");
  AliESDtrackCuts* esdTrackCutsL = 0x0;
  
  printf("\nSettings:\n");
  if (period.Contains("LHC11") || period.Contains("LHC12") || period.Contains("LHC13")) {
    esdTrackCutsL = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kTRUE);
    printf("Using standard ITS-TPC track cuts 2011.\n");
  } else if (period.Contains("LHC10")) {
    esdTrackCutsL = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kTRUE);
    printf("Using standard ITS-TPC track cuts 2010.\n");
  } else  {
    esdTrackCutsL = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kTRUE);
    printf("WARNING: Cuts not configured for this period!!! Using standard ITS-TPC track cuts 2011\n");
  }
  
  task->SetStoreMultiplicity(storeMultiplicity);
  if (task->GetStoreMultiplicity()) {
    printf("Storing multiplicity in tree!\n");
  } else {
    printf("NOT storing multiplicity in tree!\n");
  }
  
  if (useFilteredTreeCuts == kTRUE) {
    esdTrackCutsL->SetMinNClustersTPC(50);
    esdTrackCutsL->SetMinNCrossedRowsTPC(70);
    esdTrackCutsL->SetMaxDCAToVertexXY(3.0);
    esdTrackCutsL->SetMaxDCAToVertexZ(3.0);
  }

  task->SetIsPbpOrpPb(isPbpOrpPb);
  if (task->GetIsPbpOrpPb()) {
    printf("Collision type pPb/Pbp set -> Adapting vertex cuts!\n");
  } else {
    printf("Collision type different from pPb/Pbp -> Using standard vertex cuts!\n");
  }
  
  trackFilter->AddCuts(esdTrackCutsL);
  task->SetTrackFilter(trackFilter);
  task->SetUsePhiCut(usePhiCut);
  task->SetTPCcutType(tpcCutType);
  
  printf("UsePhiCut: %d\n", task->GetUsePhiCut());
  printf("UseTPCCutMIGeo: %d\n", task->GetUseTPCCutMIGeo());
  printf("UseTPCnclCut: %d\n", task->GetUseTPCnclCut());
  
  task->SetDoAdditionalQA(setDoAdditionalQA);
  
  if (task->GetDoAdditionalQA())
    printf("Storing histos for additional QA!\n");
  else
    printf("NOT storing histos for additional QA!\n");
  
  task->SetZvtxCutEvent(10.0);
  printf("Cut on z position of vertex: %.2f cm\n", task->GetZvtxCutEvent());
  
  task->SetEtaCut(0.9);
  printf("EtaCut: %.2f\n", task->GetEtaCut());
  
  task->SetPtpcPionCut(0.6);
  printf("P_TPC_Pion cut: %.2f\n", task->GetPtpcPionCut());

  task->SetUseFilteredTreeCuts(useFilteredTreeCuts);
  if (useFilteredTreeCuts == kTRUE){
  printf("Using additional cuts from filtered tree spline creation!\n");
  } else {
  printf("Using old configurations\n");
  }

  task->SetStoreNumOfSubthresholdclusters(kTRUE);
  printf("Store num subthreshold clusters: %d\n", task->GetStoreNumOfSubthresholdclusters());
  
  task->SetStoreNumClustersInActiveVolume(kTRUE);
  printf("Store num clusters in active volume: %d\n", task->GetStoreNumClustersInActiveVolume());
  
  task->SetCorrectdEdxEtaDependence(correctdEdxEtaDependence);  
  task->SetCorrectdEdxMultiplicityDependence(correctdEdxMultiplicityDependence);
  
  printf("Eta correction: %s for this task\n", 
         task->GetCorrectdEdxEtaDependence() ? "enabled (only works if enabled in PIDresponse!)" : "explicitly disabled");
  printf("Multiplicity correction: %s for this task\n\n", 
         task->GetCorrectdEdxMultiplicityDependence() ? "enabled (only works if enabled in PIDresponse!)" : "explicitly disabled");
  
  mgr->AddTask(task);


  //================================================
  //              data containers
  //================================================
  //            find input container
  //below the trunk version
  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
  mgr->ConnectInput(task,  0, cinput );

  //dummy output container
  AliAnalysisDataContainer *coutput0 =
      mgr->CreateContainer("TPCPIDEtaTree_tree",
                           TTree::Class(),
                           AliAnalysisManager::kExchangeContainer,
                           "TPCPIDEtaTree_default");
  mgr->ConnectOutput(task,  0, coutput0);

  //define output containers, please use 'username'_'somename'
  AliAnalysisDataContainer *coutput1 = 
      mgr->CreateContainer("TPCPIDEtaTree", TTree::Class(),
                           AliAnalysisManager::kOutputContainer,"TPCPIDEtaTree.root");
  mgr->ConnectOutput(task,  1, coutput1);
  
  AliAnalysisDataContainer *coutput2 = 
      mgr->CreateContainer("TPCPIDEtaTreePions", TTree::Class(),
                           AliAnalysisManager::kOutputContainer,"TPCPIDEtaTreePions.root");
  mgr->ConnectOutput(task, 2, coutput2); 
      
  if (task->GetDoAdditionalQA()) {
      AliAnalysisDataContainer *coutput3 = 
         mgr->CreateContainer("TPCPIDEtaTreeAdditionalQA", TObjArray::Class(),
                              AliAnalysisManager::kOutputContainer,"TPCPIDEtaTreeAdditonalQA.root");
      mgr->ConnectOutput(task, 3, coutput3); 
  }
  
  return task;
}
