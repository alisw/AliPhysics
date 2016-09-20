AliAnalysisTask *AddTask_HypTritEventTree(UInt_t triggerMask = AliVEvent::kINT7, Bool_t pidQa = kTRUE) {
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_lkreis_HypTritEventTree", "No analysis manager found.");
    return 0;
  }
  AliAnalysisTaskHypTritEventTree *task = new AliAnalysisTaskHypTritEventTree("lkreisTaskHypTritEventTree");
  task->SelectCollisionCandidates(triggerMask);
  task->SetPidQa(pidQa);
  //Data Containers
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput1 =
    mgr->CreateContainer("histograms", TList::Class(),AliAnalysisManager::kOutputContainer,mgr->GetCommonFileName());
  AliAnalysisDataContainer *coutput2 =
    mgr->CreateContainer("tree", TTree::Class(),AliAnalysisManager::kOutputContainer,mgr->GetCommonFileName());
  AliAnalysisDataContainer *coutput3 =
    mgr->CreateContainer("tree_mc", TTree::Class(),AliAnalysisManager::kOutputContainer,mgr->GetCommonFileName());
  mgr->ConnectInput(task, 0, cinput);
  mgr->ConnectOutput(task, 1, coutput1);
  mgr->ConnectOutput(task, 2, coutput2);
  mgr->ConnectOutput(task, 3, coutput3);
  return task;
}
