AliAnalysisTask *AddTask_HypTritTestAOD() {
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_HypTritTest", "No analysis manager found.");
    return 0;
  }
  AliAnalysisTaskHypTritTestAOD *task = new AliAnalysisTaskHypTritTestAOD("mhartungTaskHypTritTestAOD");
  //task->SelectCollisionCandidates(triggerMask);
  //task->SetTriggerMask(triggerMask);
  //task->SetPeriod(period);
  //task->SetBetheSplines(betheSplines);

  //Data Containers
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput1 =
    mgr->CreateContainer("histogramsKF", TList::Class(),AliAnalysisManager::kOutputContainer,mgr->GetCommonFileName());
  AliAnalysisDataContainer *coutput2 =
    mgr->CreateContainer("treeKF", TTree::Class(),AliAnalysisManager::kOutputContainer,mgr->GetCommonFileName());
  AliAnalysisDataContainer *coutput3 =
    mgr->CreateContainer("treeKF_mc", TTree::Class(),AliAnalysisManager::kOutputContainer,mgr->GetCommonFileName());
  mgr->ConnectInput(task, 0, cinput);
  mgr->ConnectOutput(task, 1, coutput1);
  mgr->ConnectOutput(task, 2, coutput2);
  mgr->ConnectOutput(task, 3, coutput3);
  return task;
}
