AlidNdPtAnalysisPbPbAOD *AddTask_dNdPt_PbPbAOD( UInt_t uTriggerMask = AliVEvent::kMB | AliVEvent::kCentral | AliVEvent::kSemiCentral , 
						Double_t dNCrossedRowsTPC = 100, 
						char *contName = "dNdPtPbPbAOD")
{
  Printf("===============BAUM================");
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
  task->SetCutMinNClustersTPC(0);
  //task->SetCutLengthInTPCPtDependent();
  //task->SetPrefactorLengthInTPCPtDependent(0.85);
//   task->EnableRelativeCuts();
//   task->SetCutPercMinNClustersTPC(1.0);
//   task->SetCutPercMinNCrossedRowsTPC(1.0);
  
  Double_t binsPtCheck[] = {0., 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 3.0, 4.0, 5.0, 10.0, 13.0, 15.0, 20.0, 25.0, 30.0, 40.0, 50.0, 70.0, 100.0, 150.0, 200.0}; 
  Int_t nBinPtCheck = sizeof(binsPtCheck)/sizeof(Double_t);
  task->SetBinsPtCheck(nBinPtCheck, binsPtCheck);
  ::Info("AddTask_dNdPt_PbPbAOD",Form("CrossedRowCut set to %.0f", task->GetCutMinNCrossedRowsTPC()));
  
  mgr->AddTask(task);
  
  AliAnalysisDataContainer *coutput = mgr->CreateContainer(Form("%s", contName), TList::Class(),  AliAnalysisManager::kOutputContainer, Form("%s:dNdPtHistos", mgr->GetCommonFileName()));
  
  mgr->ConnectInput( task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, coutput);
  
  return task;
}   
