// $Id$

AliEmcalTrackPropagatorTask* AddTaskEmcalTrackPropagator(
  const char *name         = "Tracks",
  const Double_t d         = -1,
  const Double_t pt        = -1
)
{  
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskEmcalTrackPropagator", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler())
  {
    ::Error("AddTaskEmcalTrackPropagator", "This task requires an input event handler");
    return NULL;
  }
  
  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------
  AliEmcalTrackPropagatorTask* propagator = new AliEmcalTrackPropagatorTask();
  propagator->SetTracksName(name);
  if (d > -1) propagator->SetDist(d);
  if (pt > -1) propagator->SetMinPt(pt);

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------
  mgr->AddTask(propagator);
  
  // Create containers for input/output
  AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer()  ;
  mgr->ConnectInput  (propagator, 0,  cinput1 );
  
  return propagator;
}
