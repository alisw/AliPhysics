AliAnalysisTaskMuonTrackingEff *AddTaskMUONTrackingEfficiency(Bool_t setDefaultTrackCuts, Bool_t isMC, TString extension = "")
{
  //
  // Task for the determination of the MUON tracking chamber efficiency
  //
  // lenhardt@subatech.in2p3.fr
  // lardeux@subatech.in2p3.fr
  //
  
  // Get the pointer to the existing analysis manager via the static access method
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskMUONTrackingEfficiency", "No analysis manager to connect to.");
    return NULL;
  }   
  
  // This task run on ESDs
  TString type = mgr->GetInputEventHandler()->GetDataType();
  if (!type.Contains("ESD")) {
    ::Error("AddTaskMUONTrackingEfficiency", "ESD input handler needed!");
    return NULL;
  }
  
  // Define output file directory
  TString fileName = AliAnalysisManager::GetCommonFileName();
  if (fileName.IsNull()) {
    ::Error("AddTaskMUONTrackingEfficiency", "Common output file is not defined!");
    return NULL;
  }
  fileName += ":MUON_Efficiency";
  TString suffix = (!extension.IsNull()) ? Form("_%s",extension.Data()) : "";
  fileName += suffix;
  
  
  // Create and configure task
  TString name = Form("MuonTrackingEfficiency%s",extension.Data());
  AliAnalysisTaskMuonTrackingEff* taskMuonTrackingEff = new AliAnalysisTaskMuonTrackingEff(name.Data());
  if (setDefaultTrackCuts) taskMuonTrackingEff->SetDefaultMuonTrackCuts(isMC);
  
  // Add to the manager
  mgr->AddTask(taskMuonTrackingEff);
  
  // Connect input container
  mgr->ConnectInput (taskMuonTrackingEff, 0, mgr->GetCommonInputContainer());
  
  // Create and connect output containers
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(Form("ClustersCounters%s",suffix.Data()), AliCounterCollection::Class(), AliAnalysisManager::kOutputContainer, fileName);
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer(Form("EventsCounters%s",suffix.Data()), AliCounterCollection::Class(), AliAnalysisManager::kOutputContainer, fileName);
  AliAnalysisDataContainer *coutput3 = mgr->CreateContainer(Form("TracksDetectedPerChamber%s",suffix.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, fileName);
  AliAnalysisDataContainer *coutput4 = mgr->CreateContainer(Form("TotalTracksPerChamber%s",suffix.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, fileName);
  AliAnalysisDataContainer *coutput5 = mgr->CreateContainer(Form("SingleDetectedPerChamber%s",suffix.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, fileName);
  AliAnalysisDataContainer *coutput6 = mgr->CreateContainer(Form("ExtraHistos%s",suffix.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, fileName);
  mgr->ConnectOutput(taskMuonTrackingEff, 1, coutput1);
  mgr->ConnectOutput(taskMuonTrackingEff, 2, coutput2);
  mgr->ConnectOutput(taskMuonTrackingEff, 3, coutput3);
  mgr->ConnectOutput(taskMuonTrackingEff, 4, coutput4);
  mgr->ConnectOutput(taskMuonTrackingEff, 5, coutput5);
  mgr->ConnectOutput(taskMuonTrackingEff, 6, coutput6);
  
  return taskMuonTrackingEff;
}

