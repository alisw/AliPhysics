AliTaskCDBconnect* AddTaskCDBconnect(const char *path/*="raw://"*/, Int_t run=0) 
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskCDBconnect", "No analysis manager to connect to.");
    return NULL;
  }   
  AliTaskCDBconnect *task= new AliTaskCDBconnect("CDBconnect", path, run);
  mgr->AddTask(task);
  AliAnalysisDataContainer *cinput1 = mgr->GetCommonInputContainer();    
  mgr->ConnectInput(task,  0, cinput1);
  return task;
}   

AliTaskCDBconnect* AddTaskCDBconnect()
{
  return AddTaskCDBconnect("cvmfs://");
}
