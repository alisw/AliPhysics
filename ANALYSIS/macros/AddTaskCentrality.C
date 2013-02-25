AliCentralitySelectionTask *AddTaskCentrality(Bool_t fillHistos=kTRUE, Bool_t aod=kFALSE)
{
// Macro to connect a centrality selection task to an existing analysis manager.
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskCentrality", "No analysis manager to connect to.");
    return NULL;
  }    
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskCentrality", "This task requires an input event handler");
    return NULL;
  }
  TString inputDataType = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  if (!aod && (inputDataType != "ESD")) {
    ::Error("AddTaskCentrality", "This task works only on ESD analysis");
    return NULL;
  }
  AliCentralitySelectionTask *centralityTask = new AliCentralitySelectionTask("CentralitySelection");
  centralityTask->SetInput(inputDataType);
  centralityTask->SelectCollisionCandidates(AliVEvent::kAny);
  mgr->AddTask(centralityTask);
  
  mgr->ConnectInput(centralityTask, 0, mgr->GetCommonInputContainer());
  if (fillHistos) {
    centralityTask->SetFillHistos();
    AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("CentralityStat",
                                                              TList::Class(), 
                                                              AliAnalysisManager::kOutputContainer,
                                                              "EventStat_temp.root");
    mgr->ConnectOutput(centralityTask,1,coutput1);
  }

  return centralityTask;
}   
