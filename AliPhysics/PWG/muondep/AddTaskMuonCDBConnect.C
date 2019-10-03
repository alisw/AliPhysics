AliAnalysisTaskMuonCDBConnect* AddTaskMuonCDBConnect()
{
  /// Add AliAnalysisTaskMuonCDConnect to the train (Philippe Pillot)
  
  // Get the pointer to the existing analysis manager via the static access method.
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if(!mgr) {
    Error("AddTaskMuonCDBConnect","AliAnalysisManager not set!");
    return NULL;
  }
  
  // Create and configure task
  AliAnalysisTaskMuonCDBConnect *task = new AliAnalysisTaskMuonCDBConnect("MuonCDBConnect");
  if (!task) {
    Error("AddTaskMuonCDBConnect", "Muon CDB connection task cannot be created!");
    return NULL;
  }
  
  // Add task to analysis manager
  mgr->AddTask(task);
  
  // Connect input container
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  
  return task;
  
}

