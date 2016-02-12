// -*- C++ -*-
// $Id$

AliAnalysisTaskADChargeMonitoring* AddTaskADChargeMonitoring() {  
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_ADChargeMonitoring", "No analysis manager found.");
    return MNULL;
  }
  
  if (!mgr->GetInputEventHandler()) {
    Error("AddTask_ADChargeMonitoring", "This task requires an input event handler");
    return 0;
  }
      
  AliAnalysisTaskADChargeMonitoring *task = new AliAnalysisTaskADChargeMonitoring();
  mgr->AddTask(task);
  
  // Create containers for input/output
  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput = mgr->CreateContainer("TE",
							   TTree::Class(),
							   AliAnalysisManager::kOutputContainer,
							   AliAnalysisManager::GetCommonFileName());

  // Connect input/output
  mgr->ConnectInput (task, 0, cinput);
  mgr->ConnectOutput(task, 1, coutput);
  
  return task;
}
