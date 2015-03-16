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
  
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("cont", TList::Class(),AliAnalysisManager::kOutputContainer, "out");
  
  mgr->ConnectOutput (task, 0, coutput1 );
  
  
  return task;
}   
