AliAnalysisCODEXtask* AddCODEXtask (TString tskname = "AliCODEXfiltering", TString suffix = "") {

  // Get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskNucleiYield", "No analysis manager found.");
    return 0x0;
  }

  // Check the analysis type using the event handlers connected to the analysis manager.
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskNucleiYield", "This task requires an input event handler");
    return 0x0;
  }

  tskname.Append(suffix.Data());
  AliAnalysisCODEXtask* task = new AliAnalysisCODEXtask(tskname.Data());
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(Form("AliCODEX%s",suffix.Data()), TTree::Class(),
                                                            AliAnalysisManager::kOutputContainer,
                                                            "AliCODEX.root");
  coutput1->SetSpecialOutput();

  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer(Form("%s_summary", tskname.Data()), TList::Class(),
							    AliAnalysisManager::kOutputContainer,
							    "AnalysisResults.root");


  mgr->ConnectInput  (task,  0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput (task,  1, coutput1);
  mgr->ConnectOutput (task,  2, coutput2);
  return task;
}
