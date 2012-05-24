// $Id$

AliAnalysisTaskRhoAverage* AddTaskRhoAverage(
   const char *nJets          = "Jets",
   const char *nTracks        = "PicoTracks",   
   const char *nClusters      = "CaloClustersCorr",  
   const char *nRho           = "Rho",
   const Double_t minPhi      = 0,
   const Double_t maxPhi      = 2 * TMath::Pi(),
   const Double_t minEta      = -0.9,
   const Double_t maxEta      = 0.9,
   const Double_t minPt       = 0.15
)
{  
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskRhoAverage", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler())
  {
    ::Error("AddTaskRhoAverage", "This task requires an input event handler");
    return NULL;
  }
  
  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------

  TString name(Form("RhoAverage_%s", nJets));
  AliAnalysisTaskRhoAverage *rhotask = new AliAnalysisTaskRhoAverage(name);
  rhotask->SetJetsName(nJets);
  rhotask->SetTracksName(nTracks);
  rhotask->SetClustersName(nClusters);
  rhotask->SetRhoName(nRho);
  rhotask->SetPhiLimits(minPhi,maxPhi);
  rhotask->SetEtaLimits(minEta,maxEta);
  rhotask->SetPtMin(minPt);

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------

  mgr->AddTask(rhotask);

  // Create containers for input/output
  mgr->ConnectInput (rhotask, 0, mgr->GetCommonInputContainer() );

  return rhotask;
}
