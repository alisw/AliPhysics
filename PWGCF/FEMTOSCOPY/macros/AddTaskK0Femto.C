AliFemtoK0Analysis *AddTaskK0Femto(bool FieldPositive = kTRUE){

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskBF", "No analysis manager to connect to.");
    return NULL;
  }

  AliFemtoK0Analysis *K0Task = new AliFemtoK0Analysis("K0Task", FieldPositive);
  if(!K0Task) exit(-1);
  mgr->AddTask(K0Task);

  TString outputFileName = AliAnalysisManager::GetCommonFileName();
  if(FieldPositive) outputFileName += ":PWGCF.outputK0Analysis_FieldPos.root";
  else outputFileName += ":PWGCF.outputK0Analysis_FieldNeg.root";
  AliAnalysisDataContainer *coutputK0 = mgr->CreateContainer("MyList", TList::Class(), AliAnalysisManager::kOutputContainer, outputFileName.Data());

  mgr->ConnectInput(K0Task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(K0Task, 1, coutputK0);

  return K0Task;
}


  
