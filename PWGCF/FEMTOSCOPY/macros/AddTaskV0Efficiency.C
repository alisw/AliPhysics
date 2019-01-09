AliAnalysisV0Efficiency *AddTaskV0Efficiency(const char* aName, bool aIgnoreInjectedV0s)
{
  // Adds the Lambda-Lambda femtoscopy task to the manger

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskV0Efficiency", "No analysis manager to connect to.");
    return NULL;
  }
  
  AliAnalysisV0Efficiency *myTask = new AliAnalysisV0Efficiency(aName, aIgnoreInjectedV0s);
  if(!myTask) return NULL;
  mgr->AddTask(myTask);

  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  if (!cinput) {
    ::Error("AddTaskV0Efficiency", "No common input container.");
    return NULL;
  }

  TString containerName = "MyListStudy";
  
  TString outputFileName = AliAnalysisManager::GetCommonFileName();
  outputFileName += ":Results";

  
  AliAnalysisDataContainer *coutput = mgr->CreateContainer(containerName, TList::Class(), AliAnalysisManager::kOutputContainer, outputFileName.Data());

  mgr->ConnectInput(myTask, 0, cinput);	
  mgr->ConnectOutput(myTask, 1, coutput);	
  
  return myTask;
}
