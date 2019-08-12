AliAnalysisTaskCheckGenKine* AddTaskCheckGenKine(TString suffix = "")
{
  // Get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskCheckGenKine", "No analysis manager found.");
    return 0x0;
  }
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskCheckGenKine", "This task requires an input event handler");
    return 0x0;
  }
  AliMCEventHandler* handler = (AliMCEventHandler*)mgr->GetMCtruthEventHandler();
  if (!handler) {
    ::Error("AddTaskCheckGenKine","MC handler not present");
    return 0x0;
  }
  
  AliAnalysisTaskCheckGenKine *taskgen = new AliAnalysisTaskCheckGenKine();
  mgr->AddTask(taskgen);
  
  TString outputFileName = AliAnalysisManager::GetCommonFileName();
  outputFileName += ":TrackGenKine";

  TString listname="listGenKine";
  listname+=suffix.Data();

  AliAnalysisDataContainer *coutput = mgr->CreateContainer(listname,
							   TList::Class(),
							   AliAnalysisManager::kOutputContainer,
							   outputFileName.Data() );

  mgr->ConnectInput  (taskgen,  0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput (taskgen,  1, coutput);
  return taskgen;
}
