AliAnalysisTaskNucleiKine* AddTaskNucleiKine() {

  AliAnalysisTaskNucleiKine *task = new AliAnalysisTaskNucleiKine();

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    cout<<"AliAnalysisTaskNucleiKine","No analysis manager to connect to."<<endl;
    return NULL;
  }

  mgr->AddTask(task);

  // Create containers for input/output
  TString outputFileName = AliAnalysisManager::GetCommonFileName();

  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput = mgr->CreateContainer("fOutputList", TList::Class(),  AliAnalysisManager::kOutputContainer,outputFileName);

  //connect containers
  mgr->ConnectInput(task, 0, cinput);
  mgr->ConnectOutput(task, 1, coutput);

  return task;
}

