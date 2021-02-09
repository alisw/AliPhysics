AliAnalysisTaskFilterHe3 *AddTask_FilterHe3(TString name = "standard")
{
  
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
  AliAnalysisTaskFilterHe3 *task = new AliAnalysisTaskFilterHe3(Form("akalweitTaskFilter_%s", name.Data()));
  //task->SetParticleType(ParticleType);
  //task->SelectCollisionCandidates(AliVEvent::kINT7|AliVEvent::kCentral|AliVEvent::kSemiCentral);
  
  //================================================
  //              data containers
  //================================================

  //dumm output container
  AliAnalysisDataContainer *coutput0 =
    mgr->CreateContainer(Form("akalweit_treeFilter_%s", name.Data()),
			 TTree::Class(),
			 AliAnalysisManager::kExchangeContainer,
			 "akalweit_default");
  
  // Create containers for input/output
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(Form("akalweit_filter_%s_hist", name.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, Form("AnalysisResults.root:Filter_%s", name.Data()));
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer(Form("akalweit_filter_%s_names", name.Data()), TTree::Class(), AliAnalysisManager::kOutputContainer, Form("AnalysisResults.root:Filter_List_%s", name.Data()));
  
  // Connect input/output
  mgr->ConnectInput(task, 0, cinput);
  mgr->ConnectOutput(task, 1, coutput1);
  mgr->ConnectOutput(task, 2, coutput2);


 
  return task;
}
