AliAnalysisTaskMuonPerformance *AddTaskMuonPerformance(Bool_t correctClusterResForSystematics = kTRUE,
						       Bool_t fitClusterResiduals = kTRUE)
{
  /// Add AliAnalysisTaskMuonPerformance to the train
  
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskMuonPerformance", "No analysis manager to connect to.");
    return NULL;
  }   
  
  // This task requires an ESD input handler.
  // Check this using the analysis manager.
  //===============================================================================
  TString type = mgr->GetInputEventHandler()->GetDataType();
  if (!type.Contains("ESD")) {
    ::Error("AddTaskMuonPerformance", "MuonPerformance task needs the manager to have an ESD input handler.");
    return NULL;
  }
  
  TString baseOutName = "muonPerformance.root";
  TString outputfile = mgr->GetCommonFileName();
  if ( ! outputfile.IsNull() ) outputfile += ":MUON_Performances";
  else outputfile = baseOutName;
  
  AliAnalysisDataContainer *coutput1  = mgr->CreateContainer("EffContainer",AliCFContainer::Class(),AliAnalysisManager::kOutputContainer,outputfile);
  AliAnalysisDataContainer *coutput2  = mgr->CreateContainer("Efficiency",TObjArray::Class(),AliAnalysisManager::kParamContainer,outputfile);
  AliAnalysisDataContainer *coutput3  = mgr->CreateContainer("TriggerResolution",TObjArray::Class(),AliAnalysisManager::kOutputContainer,outputfile);
  AliAnalysisDataContainer *coutput4  = mgr->CreateContainer("TrackerResolution",TObjArray::Class(),AliAnalysisManager::kOutputContainer,outputfile);
  AliAnalysisDataContainer *coutput5  = mgr->CreateContainer("MomentumAtVtx",TObjArray::Class(),AliAnalysisManager::kParamContainer,outputfile);
  AliAnalysisDataContainer *coutput6  = mgr->CreateContainer("SlopeAtVtx",TObjArray::Class(),AliAnalysisManager::kParamContainer,outputfile);
  AliAnalysisDataContainer *coutput7  = mgr->CreateContainer("EtaAtVtx",TObjArray::Class(),AliAnalysisManager::kParamContainer,outputfile);
  AliAnalysisDataContainer *coutput8  = mgr->CreateContainer("PhiAtVtx",TObjArray::Class(),AliAnalysisManager::kParamContainer,outputfile);
  AliAnalysisDataContainer *coutput9  = mgr->CreateContainer("MomentumAtFirstCl",TObjArray::Class(),AliAnalysisManager::kParamContainer,outputfile);
  AliAnalysisDataContainer *coutput10 = mgr->CreateContainer("SlopeAtFirstCl",TObjArray::Class(),AliAnalysisManager::kParamContainer,outputfile);
  AliAnalysisDataContainer *coutput11 = mgr->CreateContainer("DCA",TObjArray::Class(),AliAnalysisManager::kParamContainer,outputfile);
  AliAnalysisDataContainer *coutput12 = mgr->CreateContainer("Clusters",TObjArray::Class(),AliAnalysisManager::kParamContainer,outputfile);
  
  // Create the task, add it to the manager and configure it.
  //===========================================================================   
  AliAnalysisTaskMuonPerformance *muonPerformanceTask = new AliAnalysisTaskMuonPerformance("muonPerformanceTask");
  mgr->AddTask(muonPerformanceTask);
  muonPerformanceTask->CorrectClusterResForSystematics(correctClusterResForSystematics);
  muonPerformanceTask->FitClusterResiduals(fitClusterResiduals);
  
  // Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================
  mgr->ConnectInput  (muonPerformanceTask, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput (muonPerformanceTask, 1, coutput1);
  mgr->ConnectOutput (muonPerformanceTask, 2, coutput2);
  mgr->ConnectOutput (muonPerformanceTask, 3, coutput3);
  mgr->ConnectOutput (muonPerformanceTask, 4, coutput4);
  mgr->ConnectOutput (muonPerformanceTask, 5, coutput5);
  mgr->ConnectOutput (muonPerformanceTask, 6, coutput6);
  mgr->ConnectOutput (muonPerformanceTask, 7, coutput7);
  mgr->ConnectOutput (muonPerformanceTask, 8, coutput8);
  mgr->ConnectOutput (muonPerformanceTask, 9, coutput9);
  mgr->ConnectOutput (muonPerformanceTask, 10, coutput10);
  mgr->ConnectOutput (muonPerformanceTask, 11, coutput11);
  mgr->ConnectOutput (muonPerformanceTask, 12, coutput12);
  
  return muonPerformanceTask;
}   
