// $Id$

AliJetEmbedding* AddTaskJetEmbedding(
  const char     *tracksName   = "Tracks",
  const char     *clusName     = "CaloClustersCorr",
  const Float_t   minPt        = 50,
  const Float_t   maxPt        = 60,
  const Float_t   minEta       = -1,
  const Float_t   maxEta       = 1,
  const Float_t   minPhi       = 0,
  const Float_t   maxPhi       = TMath::Pi() * 2,
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
    ::Error("AddTaskJetEmbedding", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler())
  {
    ::Error("AddTaskJetEmbedding", "This task requires an input event handler");
    return NULL;
  }
  
  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------

  AliJetEmbedding *jetEmb = new AliJetEmbedding("AliJetEmbedding");
  jetEmb->SetTracksName(tracksName);
  jetEmb->SetClusName(clusName);
  jetEmb->SetEtaRange(minEta, maxEta);
  jetEmb->SetPhiRange(minPhi, maxPhi);
  jetEmb->SetPtRange(minPt, maxPt);
  jetEmb->SetCopyArray(copyArray);
  jetEmb->SetNClusters(nClus);
  jetEmb->SetNTracks(nTracks);

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------

  mgr->AddTask(jetEmb);
    
  // Create containers for input/output
  mgr->ConnectInput (jetEmb, 0, mgr->GetCommonInputContainer() );

  return jetEmb;
}
