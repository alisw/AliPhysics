AliAnalysisTaskUpcPsi2s *AddTaskUpcPsi2s(Bool_t runTree = kTRUE,Bool_t runHist = kTRUE,Bool_t runSyst = kFALSE,Int_t tracking = 0){

  
  //--- get the current analysis manager ---//
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
      Error("AddTask_UpcPsi2s", "No analysis manager found.");
      return 0;
   }
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    Error("AddTask_UpcPsi2s", "This task requires an input event handler");
    return 0;
  }
	
  TString inputDataType = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  Bool_t isMC;
  if(mgr->GetMCtruthEventHandler()) isMC = kTRUE;
  
  // Create tasks
  AliAnalysisTaskUpcPsi2s *task = new AliAnalysisTaskUpcPsi2s(inputDataType.Data());
  task->SetRunTree(runTree);
  task->SetRunHist(runHist);
  task->SetIsMC(isMC);
  task->SetRunSyst(runSyst);
  task->SetTracking(tracking);
  mgr->AddTask(task);


   // Create containers for input/output
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput = mgr->CreateContainer("JPsiTree", TTree::Class(), AliAnalysisManager::kOutputContainer,Form("%s:Psi2sCentral", AliAnalysisManager::GetCommonFileName()));
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("Psi2sTree", TTree::Class(), AliAnalysisManager::kOutputContainer, Form("%s:Psi2sCentral", AliAnalysisManager::GetCommonFileName()));
  AliAnalysisDataContainer *coutput3 = mgr->CreateContainer("ListTrig", TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:Psi2sCentral", AliAnalysisManager::GetCommonFileName()));
  AliAnalysisDataContainer *coutput4 = mgr->CreateContainer("Psi2sListHist", TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:Psi2sCentral", AliAnalysisManager::GetCommonFileName()));  

  // Connect input/output
  mgr->ConnectInput(task, 0, cinput);
  mgr->ConnectOutput(task, 1, coutput);
  mgr->ConnectOutput(task, 2, coutput2);
  mgr->ConnectOutput(task, 3, coutput3);
  mgr->ConnectOutput(task, 4, coutput4);

return task;
}
