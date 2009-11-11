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
    ::Error("AddTask", "No analysis manager to connect to.");
    return NULL;
  }   

  // Create the task
  AliAnalysisTaskTrigChEff* taskTrigChEff = new AliAnalysisTaskTrigChEff("TriggerChamberEfficiency");
  taskTrigChEff->SetUseGhostTracks(useGhosts);
  // Add to the manager
  mgr->AddTask(taskTrigChEff);

  //
  // Create containers for input/output
  AliAnalysisDataContainer *cOutputTrigChEff = mgr->CreateContainer("triggerChamberEff", TList::Class(), AliAnalysisManager::kOutputContainer, "MUON.TriggerEfficiencyMap.root");

  // Attach input
  mgr->ConnectInput(taskTrigChEff,0,mgr->GetCommonInputContainer());
  // Attach output
  mgr->ConnectOutput(taskTrigChEff,0,cOutputTrigChEff);
  
  return taskTrigChEff;
}
