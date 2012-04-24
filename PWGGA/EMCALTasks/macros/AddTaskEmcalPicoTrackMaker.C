AliEmcalPicoTrackMaker* AddTaskEmcalPicoTrackMaker(
						       const char *name       = "PicoTracks",
						       const char *inname     = "tracks",
						       const char *runPeriod  = "LHC11h",
						       AliESDtrackCuts *cuts  = 0
                                                       )
{  
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskEmcalPicoTrackMaker", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler())
  {
    ::Error("AddTaskEmcalPicoTrackMaker", "This task requires an input event handler");
    return NULL;
  }
  
  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------

  // Add aod track filter task.

  AliEmcalPicoTrackMaker *eTask = new AliEmcalPicoTrackMaker();
  eTask->SetTracksOutName(name);
  eTask->SetTracksInName(inname);
  if (!strcmp(runPeriod, "LHC11h")) {
    eTask->SetAODfilterBits(256,512,1024); // hybrid tracks for LHC11h
  }
  else {
    AliWarning(Form("Run period %s not known. AOD filter bit not set.", renPeriod));
  }
  eTask->SetESDtrackCuts(cuts);

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------

  mgr->AddTask(eTask);
  
  // Create containers for input/output
  AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer();
  
  mgr->ConnectInput  (eTask, 0,  cinput1 );
  
  return eTask;
  
}
