AliAnalysisTask *AddTask_HypTritEventTree() {

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_lkreis_HypTritEventTree", "No analysis manager found.");
    return 0;
  }

  AliAnalysisTaskHypTritEventTree *task = new AliAnalysisTaskHypTritEventTree("lkreisTaskHypTritEventTree");
  // task->SelectCollisionCandidates(AliVEvent::kMB);


  //Data Containers
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
//  AliAnalysisDataContainer *coutput0 =
//      mgr->CreateContainer("lkreis_tree",
//                           TTree::Class(),
//                           AliAnalysisManager::kExchangeContainer,
//                           "lkreis_default");

  AliAnalysisDataContainer *coutput1 =
      mgr->CreateContainer("lkreis_HypTrit", TList::Class(),
                           AliAnalysisManager::kOutputContainer,
                           mgr->GetCommonFileName());

  AliAnalysisDataContainer *coutput2 =
      mgr->CreateContainer("lkreis_HypTrit_Tree", TTree::Class(),
                           AliAnalysisManager::kOutputContainer,
                           mgr->GetCommonFileName());

  AliAnalysisDataContainer *coutput3 =
      mgr->CreateContainer("lkreis_HypTrit_TreeGen", TTree::Class(),
                           AliAnalysisManager::kOutputContainer,
                           mgr->GetCommonFileName());

  mgr->ConnectInput(task, 0, cinput);
 // mgr->ConnectOutput(task, 0, coutput0);
  mgr->ConnectOutput(task, 1, coutput1);
  mgr->ConnectOutput(task, 2, coutput2);
  mgr->ConnectOutput(task, 3, coutput3);

  return task;

}
