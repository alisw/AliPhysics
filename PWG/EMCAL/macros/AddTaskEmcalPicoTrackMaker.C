// $Id$

AliEmcalPicoTrackMaker* AddTaskEmcalPicoTrackMaker(
  const char *name         = "PicoTracks",
  const char *inname       = "tracks",
  const char *runperiod    = "",
  Bool_t      includeNoITS = kTRUE,
  AliESDtrackCuts *cuts    = 0
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

  AliEmcalPicoTrackMaker *eTask = new AliEmcalPicoTrackMaker();
  eTask->SetTracksOutName(name);
  eTask->SetTracksInName(inname);
  eTask->SetIncludeNoITS(includeNoITS);
  TString runPeriod(runperiod);
  runPeriod.ToLower();
  if (runPeriod == "lhc11h") {
    eTask->SetAODfilterBits(256,512); // hybrid tracks for LHC11h
    eTask->SetMC(kFALSE);
  }
  else if (runPeriod == "lhc11a") {
    eTask->SetAODfilterBits(256,16); // hybrid tracks for LHC11a
    eTask->SetMC(kFALSE);
  }
  else if (runPeriod == "lhc12a15a" || runPeriod == "lhc12a15e") {
    eTask->SetAODfilterBits(256,16); // hybrid tracks for LHC12a15a and LHC12a15e
    eTask->SetMC(kTRUE);
  }
  else {
    if (runPeriod.IsNull())
      ::Warning("Run period %s not known. It will use IsHybridGlobalConstrainedGlobal.");
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
