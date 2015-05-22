AliAnalysisTaskPythiaNuclei* AddTaskPythiaNuclei(TString suffix = ""){
  
  AliAnalysisTaskPythiaNuclei *task = new AliAnalysisTaskPythiaNuclei(Form("nuclei%s",suffix.Data()));

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  
  if (!mgr) {
    cout<<"AliAnalysisTaskPythiaNuclei","No analysis manager to connect to."<<endl;
    return NULL;
  }

  mgr->AddTask(task); 

  // Create containers for input/output
  TString outputFileName = AliAnalysisManager::GetCommonFileName();
  
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();      
  AliAnalysisDataContainer *coutput = mgr->CreateContainer(Form("fOutputList%s",suffix.Data()), TList::Class(),  
    AliAnalysisManager::kOutputContainer,
    outputFileName);
  
  //connect containers
  mgr->ConnectInput(task, 0, cinput);
  mgr->ConnectOutput(task, 1, coutput);

  return task;
}
