// -*- C++ -*-
// $Id$

AliAnalysisTaskADCent* AddTaskADCent() {
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_ADCent", "No analysis manager found.");
    return NULL;
  }

  if (!mgr->GetInputEventHandler()) {
    Error("AddTask_ADCent", "This task requires an input event handler");
    return NULL;
  }

  AliAnalysisTaskADCent *task = new AliAnalysisTaskADCent();
  mgr->AddTask(task);

  task->SelectCollisionCandidates(AliVEvent::kINT7);

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
