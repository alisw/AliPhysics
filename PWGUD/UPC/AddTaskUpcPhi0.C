/* Modified by Abdennacer Hamdi */
AliAnalysisTaskUpcPhi0 *AddTaskUpcPhi0(){

  //--- get the current analysis manager ---//
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
      Error("AddTask_UpcPhi0", "No analysis manager found.");
      return 0;
   }
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    Error("AddTask_UpcPhi0", "This task requires an input event handler");
    return 0;
  }
	
  TString inputDataType = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  Bool_t isMC = kFALSE;
  if(mgr->GetMCtruthEventHandler()) isMC = kTRUE;
  
  // Create tasks
  AliAnalysisTaskUpcPhi0 *task = new AliAnalysisTaskUpcPhi0(inputDataType.Data(), isMC);
  mgr->AddTask(task);
  task->SetIsMC(isMC);

   // Create containers for input/output
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("Tree", TTree::Class(), AliAnalysisManager::kOutputContainer,Form("%s:Central", AliAnalysisManager::GetCommonFileName()));
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("Histograms", TList::Class(), AliAnalysisManager::kOutputContainer,Form("%s:Central", AliAnalysisManager::GetCommonFileName()));
  AliAnalysisDataContainer *coutput3 = NULL;
  if (isMC) coutput3 = mgr->CreateContainer("MCTree", TTree::Class(), AliAnalysisManager::kOutputContainer,Form("%s:Central", AliAnalysisManager::GetCommonFileName()));

  // Connect input/output
  mgr->ConnectInput(task, 0, cinput);
  mgr->ConnectOutput(task, 1, coutput1);
  mgr->ConnectOutput(task, 2, coutput2);
  if (isMC) mgr->ConnectOutput(task, 3, coutput3);
  

return task;
}
