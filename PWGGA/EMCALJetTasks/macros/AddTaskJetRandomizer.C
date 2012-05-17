// $Id$

AliJetRandomizerTask* AddTaskJetRandomizer(
  const char     *tracksName   = "Tracks",
  const char     *clusName     = "CaloClustersCorr",
  const char     *taskName     = "JetEmbeddingTask",
  const Double_t  minEta       = -1,
  const Double_t  maxEta       = 1,
  const Double_t  minPhi       = 0,
  const Double_t  maxPhi       = TMath::Pi() * 2,
  const Int_t     nTracks      = 1,
  const Int_t     nClus        = 0,
  const Bool_t    copyArray    = kFALSE
)
{  
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskJetRandomizer", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler())
  {
    ::Error("AddTaskJetRandomizer", "This task requires an input event handler");
    return NULL;
  }
  
  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------

  AliJetRandomizerTask *jetRand = new AliJetRandomizerTask(taskName);
  jetRand->SetTracksName(tracksName);
  jetRand->SetClusName(clusName);
  jetRand->SetEtaRange(minEta, maxEta);
  jetRand->SetPhiRange(minPhi, maxPhi);
  jetRand->SetCopyArray(copyArray);
  jetRand->SetNClusters(nClus);
  jetRand->SetNTracks(nTracks);

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------

  mgr->AddTask(jetRand);
    
  // Create containers for input/output
  mgr->ConnectInput (jetRand, 0, mgr->GetCommonInputContainer() );

  return jetRand;
}
