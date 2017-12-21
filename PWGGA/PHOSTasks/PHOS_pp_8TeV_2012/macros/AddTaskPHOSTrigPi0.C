AliAnalysisTaskPHOSTrigPi0* AddTaskPHOSTrigPi0(const char* name = "PHOSTrigPi0",
					       const char* options = "",
					       const char* kPeriod = "LHC12h",
					       UInt_t offlineTriggerMask = AliVEvent::kINT7)
{
  
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskPHOSTrigPi0", "No analysis manager to connect to");
    return NULL;
  }
  
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskPHOSTrigPi0", "This task requires an input event handler");
    return NULL;
  }
  
  AliAnalysisTaskPHOSTrigPi0* task = new AliAnalysisTaskPHOSTrigPi0(Form("%sTask", name));  
  task->SelectCollisionCandidates(offlineTriggerMask);
  mgr->AddTask(task);
  
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer() );
  
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(Form("list_%s_Cluster", name),
							    TList::Class(), 
							    AliAnalysisManager::kOutputContainer, 
							    AliAnalysisManager::GetCommonFileName());
  mgr->ConnectOutput(task, 1, coutput1);
  
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer(Form("list_%s_Pi0", name),
							    TList::Class(), 
							    AliAnalysisManager::kOutputContainer,
							    AliAnalysisManager::GetCommonFileName());
  mgr->ConnectOutput(task, 2, coutput2);
  
  AliAnalysisDataContainer *coutput3 = mgr->CreateContainer(Form("list_%s_Event", name),
							    TList::Class(), 
							    AliAnalysisManager::kOutputContainer,
							    AliAnalysisManager::GetCommonFileName());
  mgr->ConnectOutput(task, 3, coutput3);

  return task;
}
