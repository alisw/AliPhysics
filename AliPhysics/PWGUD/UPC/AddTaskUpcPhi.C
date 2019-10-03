AliAnalysisTaskUpcPhi *AddTaskUpcPhi(Bool_t runTree = kTRUE,Bool_t runHist = kTRUE,Bool_t runSyst = kFALSE){

  
  //--- get the current analysis manager ---//
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
      Error("AddTask_UpcPhi", "No analysis manager found.");
      return 0;
   }
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    Error("AddTask_UpcPhi", "This task requires an input event handler");
    return 0;
  }
	
  TString inputDataType = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  Bool_t isMC;
  if(mgr->GetMCtruthEventHandler()) isMC = kTRUE;
  
  // Create tasks
  AliAnalysisTaskUpcPhi *task = new AliAnalysisTaskUpcPhi(inputDataType.Data());
  task->SetRunTree(runTree);
  task->SetRunHist(runHist);
  task->SetIsMC(isMC);
  task->SetRunSyst(runSyst);
  mgr->AddTask(task);


   // Create containers for input/output
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("ITSTree", TTree::Class(), AliAnalysisManager::kOutputContainer, Form("%s:Phi", AliAnalysisManager::GetCommonFileName()));
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("TPCTree", TTree::Class(), AliAnalysisManager::kOutputContainer, Form("%s:Phi", AliAnalysisManager::GetCommonFileName()));
  AliAnalysisDataContainer *coutput3 = mgr->CreateContainer("ListTrig", TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:Phi", AliAnalysisManager::GetCommonFileName()));
  AliAnalysisDataContainer *coutput4 = mgr->CreateContainer("PhiListHist", TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:Phi", AliAnalysisManager::GetCommonFileName()));  

  // Connect input/output
  mgr->ConnectInput(task, 0, cinput);
  mgr->ConnectOutput(task, 1, coutput1);
  mgr->ConnectOutput(task, 2, coutput2);
  mgr->ConnectOutput(task, 3, coutput3);
  mgr->ConnectOutput(task, 4, coutput4);

return task;
}
