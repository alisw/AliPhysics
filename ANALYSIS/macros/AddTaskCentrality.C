AliCentralitySelectionTask *AddTaskCentrality(Int_t passNumber = 2)
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
  if (inputDataType != "ESD") {
    ::Error("AddTaskCentrality", "This task works only on ESD analysis");
    return NULL;
  }
  AliCentralitySelectionTask *centralityTask = new AliCentralitySelectionTask("CentralitySelection");
  centralityTask->SetPass(passNumber);
  centralityTask->SelectCollisionCandidates(AliVEvent::kAny);
  mgr->AddTask(centralityTask);
  
  AliAnalysisDataContainer *cinput0 = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("CentralityStat",
                TList::Class(), AliAnalysisManager::kOutputContainer,
                "EventStat_temp.root");
  
  mgr->ConnectInput(centralityTask, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(centralityTask,1,coutput1);

  return centralityTask;
}   
