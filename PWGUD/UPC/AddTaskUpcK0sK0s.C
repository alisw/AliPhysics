AliAnalysisTaskUpcK0sK0s *AddTaskUpcK0sK0s(){

  
  //--- get the current analysis manager ---//
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
      Error("AddTask_UpcK0sK0s", "No analysis manager found.");
      return 0;
   }
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    Error("AddTask_UpcK0sK0s", "This task requires an input event handler");
    return 0;
  }
	
  TString inputDataType = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  
  // Create tasks
  AliAnalysisTaskUpcK0sK0s *task = new AliAnalysisTaskUpcK0sK0s(inputDataType.Data());
  mgr->AddTask(task);


   // Create containers for input/output
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput = mgr->CreateContainer("K0sTree", TTree::Class(), AliAnalysisManager::kOutputContainer,Form("%s:K0sK0s", AliAnalysisManager::GetCommonFileName()));
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("ListTrig", TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:K0sK0s", AliAnalysisManager::GetCommonFileName()));

  // Connect input/output
  mgr->ConnectInput(task, 0, cinput);
  mgr->ConnectOutput(task, 1, coutput);
  mgr->ConnectOutput(task, 2, coutput2);

return task;
}
