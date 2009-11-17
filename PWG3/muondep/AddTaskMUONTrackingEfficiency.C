AliAnalysisTaskMuonTrackingEff *AddTaskMUONTrackingEfficiency(Bool_t isCosmicData = kFALSE) 
{
  //
  // Task for the determination of the MUON trigger chamber efficiency
  //
  // lenhardt@subatech.in2p3.fr
  //


  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTask", "No analysis manager to connect to.");
    return NULL;
  }   

  
  const Char_t* fileName = "MUON.TrackerEfficiency.root";


  // Load the geometry
  AliMUONGeometryTransformer transformer = new AliMUONGeometryTransformer();
  transformer->LoadGeometryData();


  // Create the task
  AliAnalysisTaskMuonTrackingEff* taskMuonTrackingEff = new AliAnalysisTaskMuonTrackingEff("MuonTrackingEfficiency", transformer, isCosmicData);
  // Add to the manager
  mgr->AddTask(taskMuonTrackingEff);


  //
  // Create containers for input/output
    AliAnalysisDataContainer* cinput0  =	mgr->GetCommonInputContainer();

    AliAnalysisDataContainer *coutput0 =
	mgr->CreateContainer("TracksDetectedPerDE", TClonesArray::Class(),AliAnalysisManager::kOutputContainer,  fileName);
    AliAnalysisDataContainer *coutput1 =
	mgr->CreateContainer("TotalTracksPerDE", TClonesArray::Class(),AliAnalysisManager::kOutputContainer,  fileName);
    AliAnalysisDataContainer *coutput2 =
	mgr->CreateContainer("EfficiencyPerDE", TClonesArray::Class(),AliAnalysisManager::kOutputContainer,  fileName);
    AliAnalysisDataContainer *coutput3 =
	mgr->CreateContainer("TracksDetectedPerChamber", TClonesArray::Class(),AliAnalysisManager::kOutputContainer,  fileName);
    AliAnalysisDataContainer *coutput4 =
	mgr->CreateContainer("TotalTracksPerChamber", TClonesArray::Class(),AliAnalysisManager::kOutputContainer,  fileName);
    AliAnalysisDataContainer *coutput5 =
	mgr->CreateContainer("EfficiencyPerChamber", TClonesArray::Class(),AliAnalysisManager::kOutputContainer,  fileName);

  // Attach input
    mgr->ConnectInput (taskMuonTrackingEff, 0, cinput0 );
  // Attach output
    mgr->ConnectOutput(taskMuonTrackingEff, 0, coutput0);
    mgr->ConnectOutput(taskMuonTrackingEff, 1, coutput1);
    mgr->ConnectOutput(taskMuonTrackingEff, 2, coutput2);
    mgr->ConnectOutput(taskMuonTrackingEff, 3, coutput3);
    mgr->ConnectOutput(taskMuonTrackingEff, 4, coutput4);
    mgr->ConnectOutput(taskMuonTrackingEff, 5, coutput5);
  
    delete transformer;
    return taskMuonTrackingEff;
}
