AliAnalysisCODEXtask* AddCODEXtask (TString tskname = "AliCODEX", TString suffix = "") {

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
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer(tskname.Data(), TTree::Class(),
                                                            AliAnalysisManager::kOutputContainer,
                                                            Form("%s.root",tskname.Data()));
  coutput2->SetSpecialOutput();

  AliAnalysisDataContainer *coutput = mgr->CreateContainer(Form("%s_summary%s", tskname.Data(), suffix.Data()), TList::Class(),
							    AliAnalysisManager::kOutputContainer,
							    "AnalysisResults.root");


  mgr->ConnectInput  (task,  0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput (task,  1, coutput2);
  mgr->ConnectOutput (task,  2, coutput);
  return task;
}
