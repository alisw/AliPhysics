AliAnalysisTaskMultPbTracks * AddTaskMultPbPbTracks(const char * outfilename, AliESDtrackCuts * esdTrackCuts = 0, AliAnalysisMultPbCentralitySelector * centr)
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
    
    

  AliAnalysisTaskMultPbTracks *task = new AliAnalysisTaskMultPbTracks("TaskMultPbTracks");
  mgr->AddTask(task);
  
  // Set Cuts
  if (!esdTrackCuts)
    {
      printf("ERROR: esdTrackCuts could not be created\n");
      return;
    }  
  task->SetTrackCuts(esdTrackCuts);

  // set centrality
  task->SetCentralitySelector(centr);

  // TODO:
  // IO into folders in a file?

  // Set I/O
  AliAnalysisDataContainer *cinput0 = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("cmultPbTracksOutHM",
							    AliAnalysisMultPbTrackHistoManager::Class(),
							    AliAnalysisManager::kOutputContainer,
							    outfilename);
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("cmultPbTracksOutCT",
							    AliESDtrackCuts::Class(),
							    AliAnalysisManager::kOutputContainer,
							    outfilename);
  AliAnalysisDataContainer *coutput3 = mgr->CreateContainer("cmultPbTracksOutCM",
  							    AliAnalysisMultPbCentralitySelector::Class(),
  							    AliAnalysisManager::kOutputContainer,
  							    outfilename);

  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task,1,coutput1);
  mgr->ConnectOutput(task,2,coutput2);
  mgr->ConnectOutput(task,3,coutput3);

  return task;
}   
