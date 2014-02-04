AliAnalysisTaskPhiEffMc* AddTaskPhiEffMc(Bool_t mc=kFALSE){
  
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) 
    {
      ::Error("AddAliAnalysisTaskPhiEffMc", "No analysis manager to connect to.");
      return NULL;
    }   
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) 
    {
      ::Error("AddTaskITSsaTracks", "This task requires an input event handler");
      return NULL;
    }   
  
  TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  if(type.Contains("ESD"))
    {
      ::Error("AddTaskITSsaTracks", "This task requires to run on AOD");
      return NULL;
    }
  
  AliAnalysisTaskPhiEffMc *task = new AliAnalysisTaskPhiEffMc("PhiEffMc");
  
  TString outputFileName = AliAnalysisManager::GetCommonFileName();
  cout<<"-------------- outputFileName:  "<<outputFileName<<endl;
  
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();      
  AliAnalysisDataContainer *coutputptA = mgr->CreateContainer("fOutput", TList::Class(), AliAnalysisManager::kOutputContainer,outputFileName);
  AliAnalysisDataContainer *coutputpt1 = mgr->CreateContainer("cpidpt", AliHelperPID::Class(), AliAnalysisManager::kOutputContainer, outputFileName);
  AliAnalysisDataContainer *coutputpt2 = mgr->CreateContainer("cvcutpt", AliSpectraAODEventCuts::Class(), AliAnalysisManager::kOutputContainer, outputFileName);
  AliAnalysisDataContainer *coutputpt3 = mgr->CreateContainer("ctcutpt", AliSpectraAODTrackCuts::Class(), AliAnalysisManager::kOutputContainer, outputFileName);
  
  //connect containers
  mgr->ConnectInput(task, 0, cinput);
  mgr->ConnectOutput(task, 1, coutputptA);
  mgr->ConnectOutput(task, 2, coutputpt1);
  mgr->ConnectOutput(task, 3, coutputpt2);
  mgr->ConnectOutput(task, 4, coutputpt3);
  mgr->AddTask(task);
  return task;
}
