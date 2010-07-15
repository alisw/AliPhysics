AliAnalysisTaskMuonTrackingEff *AddTaskMUONTrackingEfficiency() 
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


  // Create the task
  AliAnalysisTaskMuonTrackingEff* taskMuonTrackingEff = new AliAnalysisTaskMuonTrackingEff("MuonTrackingEfficiency");
  // Add to the manager
  mgr->AddTask(taskMuonTrackingEff);

  AliESDInputHandler* esdH = new AliESDInputHandler();
  esdH->SetReadFriends(kFALSE);
  esdH->SetReadTags();
  mgr->SetInputEventHandler(esdH);


  //
  // Create containers for input/output
    AliAnalysisDataContainer* cinput0  =	mgr->GetCommonInputContainer();

    AliAnalysisDataContainer *coutput0 =
	mgr->CreateContainer("TracksDetectedPerDE", TList::Class(),AliAnalysisManager::kOutputContainer,  fileName);
    AliAnalysisDataContainer *coutput1 =
	mgr->CreateContainer("TotalTracksPerDE", TList::Class(),AliAnalysisManager::kOutputContainer,  fileName);
    AliAnalysisDataContainer *coutput2 =
	mgr->CreateContainer("TracksDetectedPerChamber", TList::Class(),AliAnalysisManager::kOutputContainer,  fileName);
    AliAnalysisDataContainer *coutput3 =
	mgr->CreateContainer("TotalTracksPerChamber", TList::Class(),AliAnalysisManager::kOutputContainer,  fileName);

  // Attach input
    mgr->ConnectInput (taskMuonTrackingEff, 0, cinput0 );
  // Attach output
    mgr->ConnectOutput(taskMuonTrackingEff, 0, coutput0);
    mgr->ConnectOutput(taskMuonTrackingEff, 1, coutput1);
    mgr->ConnectOutput(taskMuonTrackingEff, 2, coutput2);
    mgr->ConnectOutput(taskMuonTrackingEff, 3, coutput3);
  
    return taskMuonTrackingEff;
}
