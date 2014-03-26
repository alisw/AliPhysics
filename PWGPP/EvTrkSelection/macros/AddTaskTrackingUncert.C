AliAnalysisTask *AddTaskTrackingUncert() {
  //
  // add task of tracking uncertainty
  //
  //
  //get the current analysis manager
  //
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskTrackingUncert", "No analysis manager found.");
    return 0;
  }
  //  
  //========= Add task for standard analysis to the ANALYSIS manager ====
  //
  AliAnalysisTrackingUncertainties *task    = new AliAnalysisTrackingUncertainties("trackingUncertainty");
  //
  task->SelectCollisionCandidates(AliVEvent::kMB|AliVEvent::kINT7);
  mgr->AddTask(task);
  //  
  //  
  //======================================================================
  //              data containers
  //======================================================================
  //            find input container
  //below the trunk version
  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();

  //dummy output container
  AliAnalysisDataContainer *coutput0 = mgr->CreateContainer("dummyTreeUncert",TTree::Class(),AliAnalysisManager::kExchangeContainer,"defaultTreeUncert");

  //define output containers
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("trackingUncert", TList::Class(),AliAnalysisManager::kOutputContainer,"AnalysisResults.root");

  //connect containers
  mgr->ConnectInput  (task, 0, cinput );
  mgr->ConnectOutput (task, 0, coutput0);
  mgr->ConnectOutput (task, 1, coutput1);
  //
  //
  //
  return task;

}

