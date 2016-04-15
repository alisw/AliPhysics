AliAnalysisV0Lam *AddTaskLamLam(Int_t sysStudyType, Int_t varCutType, Bool_t flattenCent)
{
  // Adds the Lambda-Lambda femtoscopy task to the manger

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskLamLam", "No analysis manager to connect to.");
    return NULL;
  }
  
  AliAnalysisV0Lam *myTask = new AliAnalysisV0Lam("LamLamTask", sysStudyType, varCutType, flattenCent);
  if(!myTask) return NULL;
  mgr->AddTask(myTask);

  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  if (!cinput) {
    ::Error("AddTaskLamLam", "No common input container.");
    return NULL;
  }

  TString containerName = "MyListStudy";
  containerName += sysStudyType;
  containerName += "Var";
  containerName += varCutType;
  if(!flattenCent) containerName += "NotFlat";
  
  TString outputFileName = AliAnalysisManager::GetCommonFileName();
  outputFileName += ":Results";

  
  AliAnalysisDataContainer *coutput = mgr->CreateContainer(containerName, TList::Class(), AliAnalysisManager::kOutputContainer, outputFileName.Data());

  mgr->ConnectInput(myTask, 0, cinput);	
  mgr->ConnectOutput(myTask, 1, coutput);	
  
  return myTask;
}
