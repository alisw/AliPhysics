AliAnalysisTask *AddTask_jditzel_DoubleHypNucTree(UInt_t triggerMask = (AliVEvent::kAny), Bool_t pidch = kFALSE) {
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_jditzel_DoubleHypNucTree", "No analysis manager found.");
    return 0;
  }
  if (!mgr->GetInputEventHandler()) {
    return 0x0;
  }
  
  AliAnalysisTaskDoubleHypNucTree *task = new AliAnalysisTaskDoubleHypNucTree("jditzelTaskDoubleHypNucTree");
  if(!task) return 0x0;
  
  task->SetTriggerMask(triggerMask);
  task->SelectPIDcheckOnly(pidch);
  
  //Data Containers
  AliAnalysisDataContainer *cinput    = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput1  = mgr->CreateContainer("histograms", TList::Class(),AliAnalysisManager::kOutputContainer,mgr->GetCommonFileName());
  AliAnalysisDataContainer *coutput2  = mgr->CreateContainer("fTree",      TTree::Class(),AliAnalysisManager::kOutputContainer,mgr->GetCommonFileName());
  AliAnalysisDataContainer *coutput3  = mgr->CreateContainer("fTreeKF",    TTree::Class(),AliAnalysisManager::kOutputContainer,mgr->GetCommonFileName());
  AliAnalysisDataContainer *coutput4  = mgr->CreateContainer("gTree",      TTree::Class(),AliAnalysisManager::kOutputContainer,mgr->GetCommonFileName());
  AliAnalysisDataContainer *coutput5  = mgr->CreateContainer("gTreeKF",    TTree::Class(),AliAnalysisManager::kOutputContainer,mgr->GetCommonFileName());
  AliAnalysisDataContainer *coutput6  = mgr->CreateContainer("hTree",      TTree::Class(),AliAnalysisManager::kOutputContainer,mgr->GetCommonFileName());
  AliAnalysisDataContainer *coutput7  = mgr->CreateContainer("hTreeKF",    TTree::Class(),AliAnalysisManager::kOutputContainer,mgr->GetCommonFileName());
  AliAnalysisDataContainer *coutput8  = mgr->CreateContainer("iTree",      TTree::Class(),AliAnalysisManager::kOutputContainer,mgr->GetCommonFileName());
  AliAnalysisDataContainer *coutput9  = mgr->CreateContainer("iTreeKF",    TTree::Class(),AliAnalysisManager::kOutputContainer,mgr->GetCommonFileName());
  AliAnalysisDataContainer *coutput10 = mgr->CreateContainer("gTreeGen",   TTree::Class(),AliAnalysisManager::kOutputContainer,mgr->GetCommonFileName());
  AliAnalysisDataContainer *coutput11 = mgr->CreateContainer("fTreeGen",   TTree::Class(),AliAnalysisManager::kOutputContainer,mgr->GetCommonFileName());
  mgr->ConnectInput(task,  0,  cinput);
  mgr->ConnectOutput(task, 1,  coutput1);
  mgr->ConnectOutput(task, 2,  coutput2);
  mgr->ConnectOutput(task, 3,  coutput3);
  mgr->ConnectOutput(task, 4,  coutput4);
  mgr->ConnectOutput(task, 5,  coutput5);
  mgr->ConnectOutput(task, 6,  coutput6);
  mgr->ConnectOutput(task, 7,  coutput7);
  mgr->ConnectOutput(task, 8,  coutput8);
  mgr->ConnectOutput(task, 9,  coutput9);
  mgr->ConnectOutput(task, 10, coutput10);
  mgr->ConnectOutput(task, 11, coutput11);
  return task;
}
