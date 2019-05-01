AliAnalysisTaskUpcEtaC *AddTaskUpcEtaC(Bool_t runTree = kTRUE,Bool_t runHist = kTRUE,Bool_t runSyst = kFALSE,Int_t tracking = 0){

  
  //--- get the current analysis manager ---//
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
      Error("AddTask_UpcEtaC", "No analysis manager found.");
      return 0;
   }
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    Error("AddTask_UpcEtaC", "This task requires an input event handler");
    return 0;
  }
        
  TString inputDataType = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  Bool_t isMC;
  if(mgr->GetMCtruthEventHandler()) isMC = kTRUE;
  
  // Create tasks
  AliAnalysisTaskUpcEtaC *task = new AliAnalysisTaskUpcEtaC(inputDataType.Data());
  task->SetRunTree(runTree);
  task->SetRunHist(runHist);
  task->SetIsMC(isMC);
  task->SetRunSyst(runSyst);
  task->SetTracking(tracking);
  mgr->AddTask(task);


  // Create containers for input/output
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput = mgr->CreateContainer("EtaCK0sChannelTree", TTree::Class(), AliAnalysisManager::kOutputContainer,Form("%s:EtaCUpc", AliAnalysisManager::GetCommonFileName()));
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("EtaCTree", TTree::Class(), AliAnalysisManager::kOutputContainer, Form("%s:EtaCUpc", AliAnalysisManager::GetCommonFileName()));
  AliAnalysisDataContainer *coutput3 = mgr->CreateContainer("ListTrig", TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:EtaCUpc", AliAnalysisManager::GetCommonFileName()));
  AliAnalysisDataContainer *coutput4 = mgr->CreateContainer("ListHist", TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:EtaCUpc", AliAnalysisManager::GetCommonFileName()));  
  AliAnalysisDataContainer *coutput5 = mgr->CreateContainer("ListHistKstar", TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:EtaCUpc", AliAnalysisManager::GetCommonFileName()));  
  AliAnalysisDataContainer *coutput6 = mgr->CreateContainer("ListHist2Rho4Pion", TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:EtaCUpc", AliAnalysisManager::GetCommonFileName()));  
  AliAnalysisDataContainer *coutput7 = mgr->CreateContainer("ListHistK0s3PiPi4K", TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:EtaCUpc", AliAnalysisManager::GetCommonFileName()));  
  AliAnalysisDataContainer *coutput8 = mgr->CreateContainer("ListHistZDC", TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:EtaCUpc", AliAnalysisManager::GetCommonFileName()));  

  // Connect input/output
  mgr->ConnectInput(task, 0, cinput);
  mgr->ConnectOutput(task, 1, coutput);
  mgr->ConnectOutput(task, 2, coutput2);
  mgr->ConnectOutput(task, 3, coutput3);
  mgr->ConnectOutput(task, 4, coutput4);
  mgr->ConnectOutput(task, 5, coutput5);
  mgr->ConnectOutput(task, 6, coutput6);
  mgr->ConnectOutput(task, 7, coutput7);
  mgr->ConnectOutput(task, 8, coutput8);

return task;
}
