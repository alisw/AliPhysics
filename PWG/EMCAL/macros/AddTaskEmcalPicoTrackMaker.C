// $Id$

AliEmcalPicoTrackMaker* AddTaskEmcalPicoTrackMaker(
  const char *name         = "PicoTracks",
  const char *inname       = "FilterTracks",
  Double_t ptmin           = 0,
  Double_t ptmax           = 1000,
  Double_t etamin          = -10,
  Double_t etamax          = +10,
  Double_t phimin          = -10,
  Double_t phimax          = +10,
  Double_t trackeff        = 1.0,
  const char *taskName     = "AliEmcalPicoTrackMaker"
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

  AliEmcalPicoTrackMaker *eTask = new AliEmcalPicoTrackMaker(taskName);
  eTask->SetTracksOutName(name);
  eTask->SetTracksInName(inname);
  eTask->SetTrackPtLimits(ptmin, ptmax);
  eTask->SetTrackEtaLimits(etamin, etamax);
  eTask->SetTrackPhiLimits(phimin, phimax);
  eTask->SetTrackEfficiency(trackeff);

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------
  mgr->AddTask(eTask);
  
  // Create containers for input/output
  AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer();
  mgr->ConnectInput(eTask, 0, cinput1 );
  
  return eTask;
}
