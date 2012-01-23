AliAnalysisTaskSEMonitNorm *AddTaskMonitNorm() 
{
  //
  
  //


  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTask", "No analysis manager to connect to.");
    return NULL;
  }   

  // Create the task
  AliAnalysisTaskSEMonitNorm *taskMonitNorm = new AliAnalysisTaskSEMonitNorm("Monit_Norm");
  
  AliLog::SetClassDebugLevel("AliAnalysisTaskSEMonitNorm",10);
  // Add to the manager
  mgr->AddTask(taskMonitNorm);

  //
  // Create containers for input/output
  AliAnalysisDataContainer *cInputVtxESD = mgr->CreateContainer("cInputMonitNorm",TChain::Class(),AliAnalysisManager::kInputContainer);

  AliAnalysisDataContainer *cOutputList = mgr->CreateContainer("cOutputList", TList::Class(),AliAnalysisManager::kOutputContainer, "AnalysisResults.root");

  AliAnalysisDataContainer *cOutputMonitNorm1 = mgr->CreateContainer("cOutputMonitNorm1",AliCounterCollection::Class(),AliAnalysisManager::kOutputContainer, "AnalysisResults.root");

  AliAnalysisDataContainer *cOutputMonitNorm2 = mgr->CreateContainer("cOutputMonitNorm2",AliCounterCollection::Class(),AliAnalysisManager::kOutputContainer, "AnalysisResults.root");

  AliAnalysisDataContainer *cOutputMonitNorm3 = mgr->CreateContainer("cOutputMonitNorm3",AliCounterCollection::Class(),AliAnalysisManager::kOutputContainer, "AnalysisResults.root");

  AliAnalysisDataContainer *cOutputMonitNorm4 = mgr->CreateContainer("cOutputMonitNorm4",AliCounterCollection::Class(),AliAnalysisManager::kOutputContainer, "AnalysisResults.root");

  // Attach input
  mgr->ConnectInput(taskMonitNorm,0,mgr->GetCommonInputContainer());
  // Attach output
  mgr->ConnectOutput(taskMonitNorm,1,cOutputList);
  mgr->ConnectOutput(taskMonitNorm,2,cOutputMonitNorm1);
  mgr->ConnectOutput(taskMonitNorm,3,cOutputMonitNorm2);
  mgr->ConnectOutput(taskMonitNorm,4,cOutputMonitNorm3);
  mgr->ConnectOutput(taskMonitNorm,5,cOutputMonitNorm4);

  return taskMonitNorm;
}
