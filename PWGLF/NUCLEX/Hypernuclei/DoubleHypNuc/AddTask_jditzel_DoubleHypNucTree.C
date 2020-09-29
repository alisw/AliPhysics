AliAnalysisTask *AddTask_jditzel_DoubleHypNucTree(UInt_t triggerMask = AliVEvent::kAny, Bool_t betheSplines = kFALSE, Bool_t pidch = kFALSE) {
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_jditzel_DoubleHypNucTree", "No analysis manager found.");
    return 0;
  }
  
  AliAnalysisTaskDoubleHypNucTree *task = new AliAnalysisTaskDoubleHypNucTree("jditzelTaskDoubleHypNucTree");
  
  task->SelectCollisionCandidates(triggerMask);
  task->SetTriggerMask(triggerMask);
  task->SetBetheSplines(betheSplines);
  task->SelectPIDcheckOnly(pidch);
  
  //Data Containers
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("histograms", TList::Class(),AliAnalysisManager::kOutputContainer,mgr->GetCommonFileName());
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("tree4LH", TTree::Class(),AliAnalysisManager::kOutputContainer,mgr->GetCommonFileName());
  AliAnalysisDataContainer *coutput3 = mgr->CreateContainer("tree4Li", TTree::Class(),AliAnalysisManager::kOutputContainer,mgr->GetCommonFileName());
  AliAnalysisDataContainer *coutput4 = mgr->CreateContainer("tree5LHe", TTree::Class(),AliAnalysisManager::kOutputContainer,mgr->GetCommonFileName());
  AliAnalysisDataContainer *coutput5 = mgr->CreateContainer("tree5LHe2", TTree::Class(),AliAnalysisManager::kOutputContainer,mgr->GetCommonFileName());
  AliAnalysisDataContainer *coutput6 = mgr->CreateContainer("tree4LHe", TTree::Class(),AliAnalysisManager::kOutputContainer,mgr->GetCommonFileName());
  AliAnalysisDataContainer *coutput7 = mgr->CreateContainer("tree4LH3B", TTree::Class(),AliAnalysisManager::kOutputContainer,mgr->GetCommonFileName());
  AliAnalysisDataContainer *coutput8 = mgr->CreateContainer("tree4LLH", TTree::Class(),AliAnalysisManager::kOutputContainer,mgr->GetCommonFileName());
  mgr->ConnectInput(task, 0, cinput);
  mgr->ConnectOutput(task, 1, coutput1);
  mgr->ConnectOutput(task, 2, coutput2);
  mgr->ConnectOutput(task, 3, coutput3);
  mgr->ConnectOutput(task, 4, coutput4);
  mgr->ConnectOutput(task, 5, coutput5);
  mgr->ConnectOutput(task, 6, coutput6);
  mgr->ConnectOutput(task, 7, coutput7);
  mgr->ConnectOutput(task, 8, coutput8);
  return task;
}
