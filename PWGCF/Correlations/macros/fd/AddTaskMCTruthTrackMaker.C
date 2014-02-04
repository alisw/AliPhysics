AliMCTruthTrackMaker* AddTaskMCTruthTrackMaker(
  const char *outname    = "MCTruthTracks",
  Bool_t      mumothers  = kFALSE,
  Bool_t      ch         = kFALSE,
  Double_t    etamax     = 1.0
)
{  
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskMCTrackSelector", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler())
  {
    ::Error("AddTaskMCTrackSelector", "This task requires an input event handler");
    return NULL;
  }
  
  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------
  TString name("AliMCTruthTrackMaker");
  name += outname;
  AliMCTruthTrackMaker *eTask = new AliMCTruthTrackMaker(name);
  eTask->SetFillMuMothers(mumothers);
  eTask->SetTracksOutName(outname);
  eTask->SetChargedMC(ch);
  eTask->SetEtaMax(etamax);

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------
  mgr->AddTask(eTask);
  
  // Create containers for input/output
  AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer();
  mgr->ConnectInput  (eTask, 0,  cinput1 );
  
  return eTask;
}
