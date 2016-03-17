AliHadCorrTask* AddTaskHadCorr(
  const char *nTracks        = "EmcalTracks",
  const char *nClusters      = "EmcalClusters",
  const char *outClusName    = "CaloClustersCorr",
  const Double_t hadcorr     = 1,
  const Double_t minPt       = 0.15,
  const Double_t phiMatch    = 0.050,
  const Double_t etaMatch    = 0.025,
  const Double_t Eexcl       = 0,
  const Bool_t trackClus     = kTRUE,
  const Bool_t   histo       = kFALSE,
  const char *outputname     = "AnalysisResults.root"
)
{  
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskHadCorr", "No analysis manager to connect to.");
    return 0;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  AliVEventHandler* handler = mgr->GetInputEventHandler();
  if (!handler)
  {
    ::Error("AddTaskHadCorr", "This task requires an input event handler");
    return 0;
  }

  enum EDataType_t {
    kUnknown,
    kESD,
    kAOD
  };

  EDataType_t dataType = kUnknown;

  if (handler->InheritsFrom("AliESDInputHandler")) {
    dataType = kESD;
  }
  else if (handler->InheritsFrom("AliAODInputHandler")) {
    dataType = kAOD;
  }
  
  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------

  TString trackName(nTracks);
  TString clusName(nClusters);

  if (trackName == "usedefault") {
    if (dataType == kESD) {
      trackName = "Tracks";
    }
    else if (dataType == kAOD) {
      trackName = "tracks";
    }
    else {
      trackName = "";
    }
  }

  if (clusName == "usedefault") {
    if (dataType == kESD) {
      clusName = "CaloClusters";
    }
    else if (dataType == kAOD) {
      clusName = "caloClusters";
    }
    else {
      clusName = "";
    }
  }

  TString name(Form("HadCorr_%s_%s", trackName.Data(), clusName.Data()));
  if (strcmp(outClusName, "") != 0) name += Form("_%s", outClusName);
  AliHadCorrTask *hcor = new AliHadCorrTask(name, histo);
  hcor->SetOutClusName(outClusName);
  hcor->SetPhiMatch(phiMatch);
  hcor->SetEtaMatch(etaMatch);
  hcor->SetHadCorr(hadcorr);
  hcor->SetEexcl(Eexcl);
  hcor->SetTrackClus(trackClus);

  AliParticleContainer* trackCont = 0;
  if (trackName == "Tracks" || trackName == "tracks") {
    trackCont = hcor->AddTrackContainer(trackName);
  }
  else {
    trackCont = hcor->AddParticleContainer(trackName);
  }
  trackCont->SetParticlePtCut(minPt);

  AliClusterContainer *clusCont = hcor->AddClusterContainer(clusName);
  if (clusCont) clusCont->SetClusNonLinCorrEnergyCut(minPt);

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------

  mgr->AddTask(hcor);
    
  // Create containers for input/output
  mgr->ConnectInput (hcor, 0, mgr->GetCommonInputContainer() );

  if (histo) {
    AliAnalysisDataContainer *cohcor = mgr->CreateContainer(name,
							    TList::Class(),
							    AliAnalysisManager::kOutputContainer,
							    outputname);
    mgr->ConnectOutput(hcor,1,cohcor);
  }
    
  return hcor;
}
