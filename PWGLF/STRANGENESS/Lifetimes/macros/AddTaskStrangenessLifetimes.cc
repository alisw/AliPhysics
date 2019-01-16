AliAnalysisTaskStrangenessLifetimes *AddTaskStrangenessLifetimes(bool isMC = false, bool V0s = true,bool Hypertritons=true,
TString tskname = "LifetimesFiltering", TString suffix = "" ) {
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
      new AliAnalysisTaskStrangenessLifetimes(isMC, tskname.Data(), 1, Hypertritons, V0s);

  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(
      Form("%s_summary", tskname.Data()), TList::Class(),
      AliAnalysisManager::kOutputContainer, "AnalysisResults.root");

  AliAnalysisDataContainer *coutput2 =
      mgr->CreateContainer(Form("V0tree%s", suffix.Data()), TTree::Class(),
                           AliAnalysisManager::kOutputContainer, "V0tree.root");
  coutput2->SetSpecialOutput();
                    

    AliAnalysisDataContainer *coutput3 =
      mgr->CreateContainer(Form("V0treeHypertriton%s", suffix.Data()), TTree::Class(),
                           AliAnalysisManager::kOutputContainer, "V0treeHypertriton.root");
    coutput3->SetSpecialOutput();


    AliAnalysisDataContainer *coutput4 =
      mgr->CreateContainer(Form("V0treeMC%s", suffix.Data()), TTree::Class(),
                           AliAnalysisManager::kOutputContainer, "V0treeMC.root");
  coutput4->SetSpecialOutput();
    
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, coutput1);
  mgr->ConnectOutput(task, 2, coutput2); 
  mgr->ConnectOutput(task, 3, coutput3);
  mgr->ConnectOutput(task, 4, coutput4);
    

    
  return task;
}
