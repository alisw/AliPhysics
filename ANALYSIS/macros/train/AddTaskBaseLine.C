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
  
  return task;
}   
