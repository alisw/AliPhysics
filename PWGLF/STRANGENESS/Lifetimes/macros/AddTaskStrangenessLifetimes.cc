AliAnalysisTaskStrangenessLifetimes *AddTaskStrangenessLifetimes(bool isMC = false,bool V0s=false,
    TString tskname = "LifetimesFiltering", TString suffix = "std") {
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
      new AliAnalysisTaskStrangenessLifetimes(isMC, tskname.Data(),1,true,V0s);

  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(
      Form("%s_summary", tskname.Data()), TList::Class(),
      AliAnalysisManager::kOutputContainer, "AnalysisResults.root");

  AliAnalysisDataContainer *coutput2 =
      mgr->CreateContainer(Form("V0tree%s", suffix.Data()), TTree::Class(),
                           AliAnalysisManager::kOutputContainer, Form("V0tree.root:%s",suffix.Data()));
  coutput2->SetSpecialOutput();

  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, coutput1);
  mgr->ConnectOutput(task, 2, coutput2);
  return task;
}
