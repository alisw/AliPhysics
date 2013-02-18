// $Id: AddTaskJetEmbeddingFromPYTHIA.C  $

AliJetEmbeddingFromPYTHIATask* AddTaskJetEmbeddingFromPYTHIA(
  const char     *tracksName    = "Tracks",
  const char     *clusName      = "",
  const char     *cellsName     = "EMCALCells",
  const char     *MCPartName    = "",
  const char     *simpath       = "/alice/sim/2012/LHC12a15e",
  Int_t           nPtHard       = 11,
  Double_t       *ptHardScaling = 0,
  const char     *aodTreeName   = "aodTree",
  const char     *aodTracksName = "tracks",
  const char     *aodClusName   = "",
  const char     *aodCellsName  = "emcalCells",
  const char     *aodMCPartName = "mcparticles",
  const char     *runperiod     = "lhc12a15e",
  Bool_t          includeNoITS  = kTRUE,
  Double_t        minCent       = -1,
  Double_t        maxCent       = -1,
  UInt_t          mask          = AliVEvent::kAny,
  const Int_t     nTracks       = 1234567890,
  const Int_t     nClus         = 0,
  const Int_t     nCells        = 1234567890,
  const Bool_t    copyArray     = kTRUE,
  const Int_t     nFiles        = 1234567890,
  const Double_t  minPt         = 0,
  const Double_t  maxPt         = 1000,
  const Double_t  minEta        = -0.9,
  const Double_t  maxEta        = 0.9,
  const Double_t  minPhi        = 0,
  const Double_t  maxPhi        = TMath::Pi() * 2,
  const char     *taskName      = "JetEmbeddingFromPYTHIATask"
)
{  
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskJetEmbeddingFromPYTHIA", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler())
  {
    ::Error("AddTaskJetEmbeddingFromPYTHIA", "This task requires an input event handler");
    return NULL;
  }
  
  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------

  AliJetEmbeddingFromPYTHIATask *jetEmb = new AliJetEmbeddingFromPYTHIATask(taskName);
  jetEmb->SetTracksName(tracksName);
  jetEmb->SetClusName(clusName);
  jetEmb->SetCellsName(cellsName);
  jetEmb->SetMCParticlesName(MCPartName);
  jetEmb->SetAODTreeName(aodTreeName);
  jetEmb->SetAODTracksName(aodTracksName);
  jetEmb->SetAODClusName(aodClusName);
  jetEmb->SetAODCellsName(aodCellsName);
  jetEmb->SetAODMCParticlesName(aodMCPartName);
  jetEmb->SetCentralityRange(minCent, maxCent);
  jetEmb->SetTriggerMask(mask);
  jetEmb->SetNCells(nCells);
  jetEmb->SetNClusters(nClus);
  jetEmb->SetNTracks(nTracks);
  jetEmb->SetCopyArray(copyArray);
  jetEmb->SetEtaRange(minEta, maxEta);
  jetEmb->SetPhiRange(minPhi, maxPhi);
  jetEmb->SetPtRange(minPt, maxPt);

  jetEmb->SetIncludeNoITS(includeNoITS);
  TString runPeriod(runperiod);
  runPeriod.ToLower();
  if (runPeriod == "lhc12a15a" || runPeriod == "lhc12a15e") {
    jetEmb->SetAODfilterBits(256,16);
  }
  else {
    if (runPeriod.IsNull())
      ::Warning("Run period %s not known. It will use IsHybridGlobalConstrainedGlobal.");
  }

  jetEmb->SetPYTHIAPath(simpath);

  if (nPtHard > 0) {
    if (ptHardScaling==0) {
      ptHardScaling = new Double_t[nPtHard];
      for (Int_t i = 0; i < nPtHard; i++)
	ptHardScaling[i] = 1;
    }
    jetEmb->SetPtHardBinScaling(nPtHard, ptHardScaling);
  }

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------

  mgr->AddTask(jetEmb);
    
  // Create containers for input/output
  mgr->ConnectInput(jetEmb, 0, mgr->GetCommonInputContainer());

  return jetEmb;
}
