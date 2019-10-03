AliAnalysisTaskMuonRefitVtx* AddTaskMuonRefitVtx(Bool_t useMeanVtxSPD = kFALSE, Bool_t useMCVtx = kFALSE,
                                                 Bool_t useTrackVtx = kFALSE)
{
  /// Add AliAnalysisTaskMuonRefitVtx to the train (Philippe Pillot)
  
  // Get the pointer to the existing analysis manager via the static access method.
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if(!mgr) { 
    Error("AddTaskMuonRefitVtx","AliAnalysisManager not set!");
    return NULL;
  }
  
  // This task run on ESDs
  TString type = mgr->GetInputEventHandler()->GetDataType();
  if (!type.Contains("ESD")) {
    Error("AddTaskMuonRefitVtx", "ESD input handler needed!");
    return NULL;
  }
  
  // Create and configure task
  AliAnalysisTaskMuonRefitVtx *task = new AliAnalysisTaskMuonRefitVtx("MuonRefitVtx");
  if (!task) {
    Error("AddTaskMuonRefitVtx", "Muon refit vtx task cannot be created!");
    return NULL;
  }
  task->UseMeanVtxSPD(useMeanVtxSPD);
  task->UseMCVtx(useMCVtx);
  task->UseTrackVtx(useTrackVtx);
  
  // Add task to analysis manager
  mgr->AddTask(task);
  
  // Connect input container
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  
  return task;
}

