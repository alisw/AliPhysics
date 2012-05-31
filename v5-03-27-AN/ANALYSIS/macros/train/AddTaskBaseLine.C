AliAnalysisTask* AddTaskBaseLine() 
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskPhysicsSelection", "No analysis manager to connect to.");
    return NULL;
  }  
  
  AliAnalysisTask *task = new AliAnalysisTaskBaseLine("baseline");
  mgr->AddTask(task);
  
  return task;
}   
