AliAnalysisTaskPythiaMpi* AddTaskPythiaMpi(){
  
  AliAnalysisTaskPythiaMpi *task = new AliAnalysisTaskPythiaMpi("");

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) 
    {
      cout<<"AliAnalysisTaskPythiaMpi","No analysis manager to connect to."<<endl;
      return NULL;
    }

  mgr->AddTask(task);
   

  // Create containers for input/output
  TString outputFileName = AliAnalysisManager::GetCommonFileName();
  cout<<"-------------- outputFileName:  "<<outputFileName<<endl;
  
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();      
  AliAnalysisDataContainer *coutput = mgr->CreateContainer("fOutputList", TList::Class(),  AliAnalysisManager::kOutputContainer,outputFileName);
  
  //connect containers
  mgr->ConnectInput(task, 0, cinput);
  mgr->ConnectOutput(task, 1, coutput);

  return task;
}
