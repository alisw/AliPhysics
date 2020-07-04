AliAnalysisTaskUpc4Prongs *AddTaskUpc4Prongs()
{
  //--- get the current analysis manager ---//
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    Error("AddTask_Upc4Prongs", "No analysis manager found.");
    return 0;
  }

  // Check the analysis type using the event handlers connected to the analysis
  // manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler())
  {
    Error("AddTask_Upc4Prongs", "This task requires an input event handler");
    return 0;
  }

  TString inputDataType =
      mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"

  // Create tasks
  AliAnalysisTaskUpc4Prongs *task =
      new AliAnalysisTaskUpc4Prongs(inputDataType.Data());
  mgr->AddTask(task);

  // Create containers for input/output
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(
      "DataTree", TTree::Class(), AliAnalysisManager::kOutputContainer,
      Form("%s:4Prongs", AliAnalysisManager::GetCommonFileName()));

  // Connect input/output
  mgr->ConnectInput(task, 0, cinput);
  mgr->ConnectOutput(task, 1, coutput1);
  // mgr->ConnectOutput(task, 2, coutput2);

  return task;
}
