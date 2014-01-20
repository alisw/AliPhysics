AliAnalysisTaskUpcK0sK0s *AddTaskUpcK0sK0s(Bool_t runTree = kTRUE,Bool_t runHist = kTRUE){

  
  //--- get the current analysis manager ---//
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
      Error("AddTask_UpcPsi2s", "No analysis manager found.");
      return 0;
   }
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    Error("AddTask_UpcPsi2s", "This task requires an input event handler");
    return 0;
  }
	
  TString inputDataType = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  
  // Create tasks
  AliAnalysisTaskUpcK0sK0s *task = new AliAnalysisTaskUpcK0sK0s(inputDataType.Data());
  task->SetRunTree(runTree);
  task->SetRunHist(runHist);
  mgr->AddTask(task);


   // Create containers for input/output
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput = mgr->CreateContainer("K0sTree", TTree::Class(), AliAnalysisManager::kOutputContainer,Form("%s:K0sK0s", AliAnalysisManager::GetCommonFileName()));
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("Counter", TH1I::Class(), AliAnalysisManager::kOutputContainer, Form("%s:K0sK0s", AliAnalysisManager::GetCommonFileName()));
  AliAnalysisDataContainer *coutput3 = mgr->CreateContainer("ListHist", TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:K0sK0s", AliAnalysisManager::GetCommonFileName()));

  // Connect input/output
  mgr->ConnectInput(task, 0, cinput);
  mgr->ConnectOutput(task, 1, coutput);
  mgr->ConnectOutput(task, 2, coutput2);
  mgr->ConnectOutput(task, 3, coutput3);

return task;
}
