AlidNdPtAnalysisPbPbAOD *AddTask_dNdPt_PbPbAOD( UInt_t uTriggerMask = AliVEvent::kMB | AliVEvent::kCentral | AliVEvent::kSemiCentral , 
						Double_t dNCrossedRowsTPC = 120, 
						char *contName = "dNdPtPbPbAOD")
{
  // Creates, configures and attaches to the train a cascades check task.
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTask_dNdPt_PbPbAOD", "No analysis manager to connect to.");
    return NULL;
  }   
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTask_dNdPt_PbPbAOD", "This task requires an input event handler");
    return NULL;
  }   
  TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  
  // Create and configure the task
  AlidNdPtAnalysisPbPbAOD *task = new AlidNdPtAnalysisPbPbAOD("dNdPtPbPbAOD");
  //   UInt_t triggerMask = AliVEvent::kMB;
  //   triggerMask |= AliVEvent::kCentral;
  //   triggerMask |= AliVEvent::kSemiCentral;
  
  task->SelectCollisionCandidates(uTriggerMask);
  
  task->SetCutMinNCrossedRowsTPC(dNCrossedRowsTPC);
  
  ::Info("AddTask_dNdPt_PbPbAOD",Form("CrossedRowCut set to %.0f", task->GetCutMinNCrossedRowsTPC()));
  
  mgr->AddTask(task);
  
  AliAnalysisDataContainer *coutput = mgr->CreateContainer(Form("%s", contName), TList::Class(),  AliAnalysisManager::kOutputContainer, Form("%s:dNdPtHistos", mgr->GetCommonFileName()));
  
  mgr->ConnectInput( task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, coutput);
  
  return task;
}   
