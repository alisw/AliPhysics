AliAnalysisTaskSE* AddTaskNanoAODFilter(Int_t iMC, Bool_t savecuts = 0) 
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskNanoAODFilter", "No analysis manager to connect to.");
    return NULL;
  }  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskNanoAODFilter", "This task requires an input event handler");
    return NULL;
  }
    
  // Configure analysis
  //===========================================================================
  AliAnalysisTaskNanoAODFilter* task = new AliAnalysisTaskNanoAODFilter("NanoAODFilter", savecuts);
  
  mgr->AddTask(task);
  task->SetMCMode(iMC);

  mgr->ConnectInput (task, 0, mgr->GetCommonInputContainer());

  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("normalisation", AliNanoFilterNormalisation::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s", AliAnalysisManager::GetCommonFileName(), "NanoAODFilter"));
  mgr->ConnectOutput(task, 1, coutput1);
  
  if (savecuts) {
    AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("qa", TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s", AliAnalysisManager::GetCommonFileName(), "NanoAODFilter"));
    mgr->ConnectOutput(task, 2, coutput2);
  }

  return task;
}
