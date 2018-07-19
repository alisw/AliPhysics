AliAnalysisTaskStrangenessLifetimes *AddTaskStrangenessLifetimes(
    TString tskname = "LifetimesFiltering", TString suffix = "") {
  // Get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskStrangenessLifetimes", "No analysis manager found.");
    return 0x0;
  }

  // Check the analysis type using the event handlers connected to the analysis
  // manager.
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskStrangenessLifetimes",
            "This task requires an input event handler");
    return 0x0;
  }

  tskname.Append(suffix.Data());
  AliAnalysisTaskStrangenessLifetimes *task =
      new AliAnalysisTaskStrangenessLifetimes(tskname.Data());

  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(
      Form("%s_summary", tskname.Data()), TList::Class(),
      AliAnalysisManager::kOutputContainer, "AnalysisResults.root");

  AliAnalysisDataContainer *coutput2 =
      mgr->CreateContainer(Form("V0tree%s", suffix.Data()), TTree::Class(),
                           AliAnalysisManager::kOutputContainer, "V0tree.root");
  coutput2->SetSpecialOutput();

  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, coutput1);
  mgr->ConnectOutput(task, 2, coutput2);
  return task;
}
