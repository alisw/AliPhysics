AliFemtoK0Analysis *AddTaskK0Femto(){

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskK0Femto", "No analysis manager to connect to.");
    return NULL;
  }

  AliFemtoK0Analysis *K0Task = new AliFemtoK0Analysis("K0Task");
  if(!K0Task) return NULL;
  mgr->AddTask(K0Task);

  TString outputFileName = AliAnalysisManager::GetCommonFileName();
  outputFileName += ":PWGCF.outputK0Analysis.root";
  AliAnalysisDataContainer *coutputK0 = mgr->CreateContainer("MyList", TList::Class(), AliAnalysisManager::kOutputContainer, outputFileName.Data());

  mgr->ConnectInput(K0Task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(K0Task, 1, coutputK0);

  return K0Task;
}


  
