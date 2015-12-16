///
/// This task will create, from an AOD (whatever kind, as long as it contains muon tracks)
/// an AliAODDimuon object and attach it to the AliAODEvent, so it can be used
/// in subsequent tasks. It is not intended to be made persistent.
///
/// \author L. Aphecetche

AliAnalysisTask* AddTaskOnTheFlyAliAODDimuon()
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskMuMu", "No analysis manager to connect to.");
    return NULL;
  }

  AliAnalysisTask* task = new AliAnalysisTaskOnTheFlyAliAODDimuon;
  mgr->AddTask(task);
  
  // Create containers for input/output
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  
  // Connect input/output
  mgr->ConnectInput(task, 0, cinput);
  
  return task;
}