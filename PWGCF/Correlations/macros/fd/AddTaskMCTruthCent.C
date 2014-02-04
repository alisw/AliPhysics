AliMCTruthCent* AddTaskMCTruthCent(Bool_t doHists = kFALSE)
{
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskMCTruthCent", "No analysis manager to connect to.");
    return NULL;
  }
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler())
  {
    ::Error("AddTaskMCTruthCent", "This task requires an input event handler");
    return NULL;
  }
  
  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------
  TString name("AliMCTruthCent");
  AliMCTruthCent *eTask = new AliMCTruthCent(name);
  if (doHists) {
    eTask->SetFillHistos();
  }
  
  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------
  mgr->AddTask(eTask);
  
  // Create containers for input/output
  AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer();
  mgr->ConnectInput  (eTask, 0,  cinput1 );
  
  if (doHists) {
    AliAnalysisDataContainer *coutput = mgr->CreateContainer("MCTruthCentStat",
                                                             TList::Class(),
                                                             AliAnalysisManager::kOutputContainer,
                                                             "MCTruthCent.root");
    mgr->ConnectOutput(eTask,1,coutput);
  }
  
  return eTask;
}
