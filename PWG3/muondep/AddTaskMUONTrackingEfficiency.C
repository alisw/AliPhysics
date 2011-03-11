AliAnalysisTaskMuonTrackingEff *AddTaskMUONTrackingEfficiency(Bool_t matchTrig = kTRUE, Bool_t applyAccCut = kTRUE) 
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
  
  
  // Create and configure task
  AliAnalysisTaskMuonTrackingEff* taskMuonTrackingEff = new AliAnalysisTaskMuonTrackingEff("MuonTrackingEfficiency");
  taskMuonTrackingEff->MatchTrigger(matchTrig);
  taskMuonTrackingEff->ApplyAccCut(applyAccCut);
  
  // Add to the manager
  mgr->AddTask(taskMuonTrackingEff);
  
  // Connect input container
  mgr->ConnectInput (taskMuonTrackingEff, 0, mgr->GetCommonInputContainer());
  
  // Create and connect output containers
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("TracksDetectedPerDE", TList::Class(),AliAnalysisManager::kOutputContainer,  fileName);
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("TotalTracksPerDE", TList::Class(),AliAnalysisManager::kOutputContainer,  fileName);
  AliAnalysisDataContainer *coutput3 = mgr->CreateContainer("TracksDetectedPerChamber", TList::Class(),AliAnalysisManager::kOutputContainer,  fileName);
  AliAnalysisDataContainer *coutput4 = mgr->CreateContainer("TotalTracksPerChamber", TList::Class(),AliAnalysisManager::kOutputContainer,  fileName);
  mgr->ConnectOutput(taskMuonTrackingEff, 1, coutput1);
  mgr->ConnectOutput(taskMuonTrackingEff, 2, coutput2);
  mgr->ConnectOutput(taskMuonTrackingEff, 3, coutput3);
  mgr->ConnectOutput(taskMuonTrackingEff, 4, coutput4);
  
  return taskMuonTrackingEff;
}
