AliAnalysisTaskTrigChEff *AddTaskMTRchamberEfficiency(Bool_t useGhosts = kFALSE) 
{
  //
  // Task for the determination of the MUON trigger chamber efficiency
  //
  // stocco@subatech.in2p3.fr
  //


  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskMTRchamberEfficiency", "No analysis manager to connect to.");
    return NULL;
  }   

  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AliAnalysisTaskTrigChEff", "This task requires an input event handler");
    return NULL;
  }

  // Create the task
  AliAnalysisTaskTrigChEff* taskTrigChEff = new AliAnalysisTaskTrigChEff("TriggerChamberEfficiency");
  taskTrigChEff->SetUseGhostTracks(useGhosts);
  // Add to the manager
  mgr->AddTask(taskTrigChEff);

  //
  // Create containers for input/output
  AliAnalysisDataContainer *cOutputTrigChEff = mgr->CreateContainer("triggerChamberEff", TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:MUON.TriggerEfficiencyMap", mgr->GetCommonFileName()));

  // Attach input
  mgr->ConnectInput(taskTrigChEff,0,mgr->GetCommonInputContainer());
  // Attach output
  mgr->ConnectOutput(taskTrigChEff,1,cOutputTrigChEff);
  
  return taskTrigChEff;
}
