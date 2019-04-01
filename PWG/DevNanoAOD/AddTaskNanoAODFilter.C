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
  
  if(savecuts) {
    AliAnalysisDataContainer *coutputpt1 = mgr->CreateContainer("evtcuts", AliAnalysisCuts::Class(),  AliAnalysisManager::kOutputContainer,"cuts.root");
    AliAnalysisDataContainer *coutputpt2 = mgr->CreateContainer("trkcuts", AliAnalysisCuts::Class(),  AliAnalysisManager::kOutputContainer,"cuts.root");
    AliAnalysisDataContainer *coutputpt3 = mgr->CreateContainer("v0cuts", AliAnalysisCuts::Class(),  AliAnalysisManager::kOutputContainer,"cuts.root");
    mgr->ConnectOutput(task, 1, coutputpt1);
    mgr->ConnectOutput(task, 2, coutputpt2);
    mgr->ConnectOutput(task, 3, coutputpt3);
  }

  return task;
}
