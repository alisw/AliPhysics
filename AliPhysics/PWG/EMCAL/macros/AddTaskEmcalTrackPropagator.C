AliEmcalTrackPropagatorTask* AddTaskEmcalTrackPropagator(const char* nTracks       = "tracks",
                                                         const Double_t d          = 440,
                                                         const Bool_t onlyIfNotSet = kTRUE,
                                                         const Bool_t onlyIfEmcal  = kFALSE)
{  
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskEmcalTrackPropagator", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  AliVEventHandler *evhand = mgr->GetInputEventHandler();
  if (!evhand) {
    ::Error("AddTaskEmcalTrackPropagator", "This task requires an input event handler");
    return NULL;
  }

  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------
  TString tname(Form("AliEmcalTrackPropagatorTask_%s", nTracks));
  AliEmcalTrackPropagatorTask* propagator = new AliEmcalTrackPropagatorTask(tname);
  AliParticleContainer *trackCont = propagator->AddParticleContainer(nTracks);
  if (d > 0) propagator->SetDist(d);
  propagator->SetOnlyIfNotSet(onlyIfNotSet);
  propagator->SetOnlyIfEmcal(onlyIfEmcal);

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------
  mgr->AddTask(propagator);
  
  // Create containers for input/output
  AliAnalysisDataContainer *cinput1 = mgr->GetCommonInputContainer();
  mgr->ConnectInput(propagator, 0, cinput1 );
  
  return propagator;
}
