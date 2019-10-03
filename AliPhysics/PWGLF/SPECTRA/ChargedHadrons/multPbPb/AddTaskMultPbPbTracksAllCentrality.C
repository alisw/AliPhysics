AliAnalysisTaskMultPbTracks ** AddTaskMultPbPbTracksAllCentrality(TString outfilename, AliESDtrackCuts * esdTrackCuts = 0, AliAnalysisMultPbCentralitySelector * centr, Int_t ncentr, Float_t * minCentr, Float_t *maxCentr)
{
  // TODO: add some parameters to set the centrality for this task, and maybe the name of the task
  // TODO: shall I use the same file and different dirs for the different centralities?

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskPhysicsSelection", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskPhysicsSelection", "This task requires an input event handler");
    return NULL;
  }
  TString inputDataType = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  
  if (inputDataType != "ESD") {
    Printf("ERROR! This task can only run on ESDs!");
  }

  // Configure analysis
  //===========================================================================
    
    
  cout << "Booking " << ncentr << " Tasks" << endl;
  AliAnalysisTaskMultPbTracks ** tasks = new AliAnalysisTaskMultPbTracks*[ncentr];
  
  for(Int_t icentr = 0; icentr < ncentr; icentr++){  

    tasks[icentr] = new AliAnalysisTaskMultPbTracks(Form("TaskMultPbTracks_%d",icentr));
    mgr->AddTask(tasks[icentr]);
  
    // Set Cuts
    if (!esdTrackCuts)
      {
	printf("ERROR: esdTrackCuts could not be created\n");
	return;
      }  
    tasks[icentr]->SetTrackCuts(esdTrackCuts);

    // set centrality
    AliAnalysisMultPbCentralitySelector * centrBin = (AliAnalysisMultPbCentralitySelector*) centr->Clone();
    centrBin->SetMultRange(minCentr[icentr],maxCentr[icentr]);
    tasks[icentr]->SetCentralitySelector(centrBin);

    // TODO:
    // IO into folders in a file?

    // Set I/O
    TString outfilenameCentr = outfilename;
    outfilenameCentr.ReplaceAll(".root",Form("_%2.2d.root",icentr));
    AliAnalysisDataContainer *cinput0 = mgr->GetCommonInputContainer();
    AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(Form("cmultPbTracksOutHM_%d",icentr),
							      AliAnalysisMultPbTrackHistoManager::Class(),
							      AliAnalysisManager::kOutputContainer,
							      outfilenameCentr.Data());
    AliAnalysisDataContainer *coutput2 = mgr->CreateContainer(Form("cmultPbTracksOutCT_%d",icentr),
							      AliESDtrackCuts::Class(),
							      AliAnalysisManager::kOutputContainer,
							      outfilenameCentr.Data());
    AliAnalysisDataContainer *coutput3 = mgr->CreateContainer(Form("cmultPbTracksOutCM_%d",icentr),
							      AliAnalysisMultPbCentralitySelector::Class(),
							      AliAnalysisManager::kOutputContainer,
							      outfilenameCentr.Data());

    mgr->ConnectInput(tasks[icentr], 0, mgr->GetCommonInputContainer());
    mgr->ConnectOutput(tasks[icentr],1,coutput1);
    mgr->ConnectOutput(tasks[icentr],2,coutput2);
    mgr->ConnectOutput(tasks[icentr],3,coutput3);

  }
  return tasks;
}   
