AliAnalysisTaskNucleiSkimAOD* AddTaskNucleiSkimAOD (TString tskname = "NucleiSkimAOD", TString suffix = "") {

  // Get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskNucleiSkimAOD", "No analysis manager found.");
    return 0x0;
  }

  // Check the analysis type using the event handlers connected to the analysis manager.
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskNucleiSkimAOD", "This task requires an input event handler");
    return 0x0;
  }

  tskname.Append(suffix.Data());

  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(Form("%s_summary", tskname.Data()), TList::Class(),
                                                            AliAnalysisManager::kOutputContainer,
                                                            "AnalysisResults.root");


  mgr->ConnectInput  (task,  0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput (task,  1, coutput1);
  return task;
}
