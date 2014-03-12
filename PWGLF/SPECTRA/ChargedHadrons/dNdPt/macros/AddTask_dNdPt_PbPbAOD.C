AlidNdPtAnalysisPbPbAOD *AddTask_dNdPt_PbPbAOD( UInt_t uTriggerMask = AliVEvent::kMB | AliVEvent::kCentral | AliVEvent::kSemiCentral , 
						Double_t dNCrossedRowsTPC = 100,
						Int_t iFilterBit = AliAODTrack::kTrkGlobal,
						char *contName = "dNdPtPbPbAOD",
						Double_t dNClustersTPC = 0,
						Bool_t bDoCutTPCLength = kTRUE
											  )
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
  task->SetCutMinNClustersTPC(dNClustersTPC);
  task->SetCutLengthInTPCPtDependent(bDoCutTPCLength);
  //task->SetCutLengthInTPCPtDependent();
  //task->SetPrefactorLengthInTPCPtDependent(0.85);
//   task->EnableRelativeCuts();
//   task->SetCutPercMinNClustersTPC(1.0);
//   task->SetCutPercMinNCrossedRowsTPC(1.0);
  
  Double_t binsPtCheck[] = {0., 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 3.0, 4.0, 5.0, 10.0, 13.0, 15.0, 20.0, 25.0, 30.0, 40.0, 50.0, 70.0, 100.0, 150.0, 200.0}; 
  Int_t nBinPtCheck = sizeof(binsPtCheck)/sizeof(Double_t);
  task->SetBinsPtCheck(nBinPtCheck, binsPtCheck);
  
  Double_t binsPhi[] = {0 ,0.10472 ,0.20944 ,0.314159 ,0.418879 ,0.523599 ,0.628319 ,0.733038 ,0.837758 ,0.942478 ,1.0472 ,1.15192 ,1.25664 ,1.36136 ,1.46608 ,1.5708 ,1.67552 ,1.78024 ,1.88496 ,1.98968 ,2.0944 ,2.19911 ,2.30383 ,2.40855 ,2.51327 ,2.61799 ,2.72271 ,2.82743 ,2.93215 ,3.03687 ,3.14159 ,3.24631 ,3.35103 ,3.45575 ,3.56047 ,3.66519 ,3.76991 ,3.87463 ,3.97935 ,4.08407 ,4.18879 ,4.29351 ,4.39823 ,4.50295 ,4.60767 ,4.71239 ,4.81711 ,4.92183 ,5.02655 ,5.13127 ,5.23599 ,5.34071 ,5.44543 ,5.55015 ,5.65487 ,5.75959 ,5.86431 ,5.96903 ,6.07375 ,6.17847 ,6.28319};
  Int_t nBinPhi = sizeof(binsPhi)/sizeof(Double_t);
  task->SetBinsPhi(nBinPhi, binsPhi);
    
  task->SetFilterBit(iFilterBit);
  
  ::Info("AddTask_dNdPt_PbPbAOD",Form("CrossedRowCut set to %.0f", task->GetCutMinNCrossedRowsTPC()));
  ::Info("AddTask_dNdPt_PbPbAOD",Form("FilterBit set to %d", task->GetFilterBit()));
  
  mgr->AddTask(task);
  
  AliAnalysisDataContainer *coutput = mgr->CreateContainer(Form("%s", contName), TList::Class(),  AliAnalysisManager::kOutputContainer, Form("%s:dNdPtHistos", mgr->GetCommonFileName()));
  
  mgr->ConnectInput( task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, coutput);
  
  return task;
}   
