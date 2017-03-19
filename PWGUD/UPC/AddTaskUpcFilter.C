AliAnalysisTaskUpcFilter *AddTaskUpcFilter() {
 
  //--- get the current analysis manager ---//
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
      Error("AddTask_UpcFilter", "No analysis manager found.");
      return 0;
   }
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    Error("AddTask_UpcFilter", "This task requires an input event handler");
    return 0;
  }

  Bool_t isESD = kFALSE;
  TString inputDataType = mgr->GetInputEventHandler()->GetDataType();
  if( inputDataType.Contains("ESD") ) isESD = kTRUE;
  Bool_t isMC = kFALSE;
  if( mgr->GetMCtruthEventHandler() ) isMC = kTRUE;

  // Create tasks
  AliAnalysisTaskUpcFilter *task = new AliAnalysisTaskUpcFilter();
  task->SetIsESD( isESD );
  task->SetIsMC( isMC );
  task->SetAllTrg(kTRUE);
  //task->SetTrgClass(15, kTRUE);
  //task->SetTrgClass(16, kTRUE);
  mgr->AddTask(task);

  // Create containers for input/output
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput = mgr->CreateContainer("UPCTree", TTree::Class(), AliAnalysisManager::kOutputContainer, Form("%s:UpcFilter", AliAnalysisManager::GetCommonFileName()));
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("HistList", TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:UpcFilter", AliAnalysisManager::GetCommonFileName()));

  // Connect input/output
  mgr->ConnectInput(task, 0, cinput);
  mgr->ConnectOutput(task, 1, coutput);
  mgr->ConnectOutput(task, 2, coutput2);

return task;

}

