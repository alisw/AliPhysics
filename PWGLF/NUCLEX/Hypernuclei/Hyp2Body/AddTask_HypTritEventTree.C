AliAnalysisTask *AddTask_HypTritEventTree() {

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_lkreis_HypTritTree", "No analysis manager found.");
    return 0;
  }

  AliAnalysisTaskHypTritEventTree *task = new AliAnalysisTaskHypTritEventTree("taskHypTritEventTree");
  task->SelectCollisionCandidates(AliVEvent::kMB+AliVEvent::kCentral+AliVEvent::kSemiCentral);

  //Data Containers
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();

  AliAnalysisDataContainer *coutput1 =
      mgr->CreateContainer("lkreis_HypTrit", TList::Class(),
                           AliAnalysisManager::kOutputContainer,
                           mgr->GetCommonFileName());

  AliAnalysisDataContainer *coutput2 =
      mgr->CreateContainer("lkreis_HypTrit_Tree", TTree::Class(),
                           AliAnalysisManager::kOutputContainer,
                           mgr->GetCommonFileName());

  mgr->ConnectInput(task, 0, cinput);
  mgr->ConnectOutput(task, 1, coutput1);
  mgr->ConnectOutput(task, 2, coutput2);

  return task;

}

