// $Id$

AliEmcalPicoTrackFromJetMaker* AddTaskEmcalPicoTrackFromJetMaker(
  const char *name         = "PicoTracks",
  const char *inname       = "jets",
  const char *taskName     = "AliEmcalPicoTrackFromJetMaker"
)
{  
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskEmcalPicoTrackFromJetMaker", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler())
  {
    ::Error("AddTaskEmcalPicoTrackFromJetMaker", "This task requires an input event handler");
    return NULL;
  }
  
  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------

  AliEmcalPicoTrackFromJetMaker *eTask = new AliEmcalPicoTrackFromJetMaker(taskName);
  eTask->SetTracksOutName(name);
  eTask->SetJetsInName(inname);

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------
  mgr->AddTask(eTask);
  
  // Create containers for input/output
  AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer();
  mgr->ConnectInput(eTask, 0, cinput1 );
  
  return eTask;
}
