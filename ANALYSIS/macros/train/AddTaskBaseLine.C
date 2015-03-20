AliAnalysisTask* AddTaskBaseLine() 
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskPhysicsSelection", "No analysis manager to connect to.");
    return NULL;
  }  
  
  AliAnalysisTask *task = new AliAnalysisTaskBaseLine("baseline");
  mgr->AddTask(task);

  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());

  // TODO adding the output container removes some warnings. however, this needs also changes in the baseline task to post an output
//   AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("cont", TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s", AliAnalysisManager::GetCommonFileName(), "baselinefolder"));
//   mgr->ConnectOutput (task, 0, coutput1 );
  
  return task;
}   
