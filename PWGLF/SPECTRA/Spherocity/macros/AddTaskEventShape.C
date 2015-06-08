AliEventShapeSelectionTask *AddTaskEventShape(Bool_t fillHistos=kTRUE, Bool_t aod=kFALSE)
{
// Macro to connect a centrality selection task to an existing analysis manager.
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskEventShape", "No analysis manager to connect to.");
    return NULL;
  }    
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskEventShape", "This task requires an input event handler");
    return NULL;
  }
  TString inputDataType = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  if (!aod && (inputDataType != "ESD")) {
    ::Error("AddTaskEventShape", "This task works only on ESD analysis");
    return NULL;
  }
  //
  AliInputEventHandler* hdl = (AliInputEventHandler*)mgr->GetInputEventHandler();
  if (hdl) hdl->SetNeedField(); 
  //
  AliEventShapeSelectionTask *eventshapeTask = new AliEventShapeSelectionTask("EventShapeSelection");
  eventshapeTask->SetInput(inputDataType);
  eventshapeTask->SelectCollisionCandidates(AliVEvent::kAny);
  mgr->AddTask(eventshapeTask);
  
  mgr->ConnectInput(eventshapeTask, 0, mgr->GetCommonInputContainer());
  if (fillHistos) {
    centralityTask->SetFillHistos();
    AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("EventShapeStat",
                                                              TList::Class(), 
                                                              AliAnalysisManager::kOutputContainer,
                                                              "EventStat_temp.root");
    mgr->ConnectOutput(eventshapeTask,1,coutput1);
  }

  return centralityTask;
}   
