AliAnalysisTaskFilterHe3 *AddTask_FilterHe3(){


  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_FilterHe3", "No analysis manager found.");
    return 0;
  }

  // Check the analysis type using the event handlers connected to the analysis manager.
   //==============================================================================
   if (!mgr->GetInputEventHandler()) {
      ::Error("AddTask_FilterHe3", "This task requires an input event handler");
      return NULL;
   }  


  //========= Add task to the ANALYSIS manager =====
  AliAnalysisTaskFilterHe3* task = new AliAnalysisTaskFilterHe3("akalweitTaskFilterHe3");
  //task->SelectCollisionCandidates(AliVEvent::kINT7|AliVEvent::kCentral|AliVEvent::kSemiCentral);

   
  //================================================
  //              data containers
  //================================================


  //dumm output container
  AliAnalysisDataContainer *coutput0 =
      mgr->CreateContainer("akalweit_treeFilterHe3",
                           TTree::Class(),
                           AliAnalysisManager::kExchangeContainer,
                           "akalweit_default");

// Create containers for input/output
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("akalweit_filterHe3_hist", TList::Class(),AliAnalysisManager::kOutputContainer, "AnalysisResults.root:FilterHe3");
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("akalweit_filterHe3_names", TTree::Class(),AliAnalysisManager::kOutputContainer, "AnalysisResults.root:FilterHe3_List");


  // Connect input/output
  mgr->ConnectInput(task, 0, cinput);
  mgr->ConnectOutput(task, 1, coutput1);
  mgr->ConnectOutput(task, 2, coutput2);


 
  return task;
}
