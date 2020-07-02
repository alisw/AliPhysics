AliAnalysisTaskSatellites *AddTaskSatellites(){

  //--- get the current analysis manager ---//
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
      Error("AddTaskSatellites", "No analysis manager found.");
      return 0;
   }

  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    Error("AddTaskSatellites", "This task requires an input event handler");
    return 0;
  }

  TString inputDataType = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"  

  // Create tasks
  AliAnalysisTaskSatellites *task = new AliAnalysisTaskSatellites(inputDataType.Data());
  mgr->AddTask(task);
    
  // Create containers for input/output
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutList = mgr->CreateContainer("ListHist", TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:Satellites", AliAnalysisManager::GetCommonFileName()));
  AliAnalysisDataContainer *coutTree = mgr->CreateContainer("Tree", TTree::Class(), AliAnalysisManager::kOutputContainer, Form("%s:Satellites", AliAnalysisManager::GetCommonFileName()));

  // Connect input/output
  mgr->ConnectInput(task, 0, cinput);
  mgr->ConnectOutput(task, 1, coutList);
  mgr->ConnectOutput(task, 2, coutTree);

return task;
}
