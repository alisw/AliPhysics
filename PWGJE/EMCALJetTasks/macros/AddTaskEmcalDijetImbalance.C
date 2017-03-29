AliAnalysisTaskEmcalDijetImbalance* AddTaskEmcalDijetImbalance(
  const char *ntracks            = "usedefault",
  const char *nclusters          = "usedefault",
  const Double_t deltaPhiMin     = 2*TMath::Pi()/3,   // Minimum delta phi between di-jets
  const Bool_t doGeomMatching    = kFALSE,            // Set whether to enable constituent study with geometrical matching
  const Double_t minTrPtHardCore = 3.0,               // Minimum track pT in high-threshold track container (for hard-core jets)
  const Double_t minClPtHardCore = 3.0,               // Minimum cluster E in standard cluster container (for hard-core jets)
  const Double_t jetR            = 0.2,               // jet R (for hard-core jets)
  const Bool_t includePHOS       = kTRUE,             // Set whether to include PHOS clusters (if enabled, must also include phos clusters in jet finder)
  const Double_t minTrPt         = 0.15,              // Minimum track pT in standard track container
  const Double_t minClPt         = 0.30,              // Minimum cluster E in standard cluster container
  const char *suffix             = ""
)
{
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskEmcalDijetImbalance", "No analysis manager to connect to.");
    return 0;
  }

  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  AliVEventHandler* handler = mgr->GetInputEventHandler();
  if (!handler)
  {
    ::Error("AddTaskEmcalDijetImbalance", "This task requires an input event handler");
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

  TString trackName(ntracks);
  TString clusName(nclusters);

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

  TString name("AliAnalysisTaskEmcalDijetImbalance");
  if (!trackName.IsNull()) {
    name += "_";
    name += trackName;
  }
  if (!clusName.IsNull()) {
    name += "_";
    name += clusName;
  }
  if (strcmp(suffix,"") != 0) {
    name += "_";
    name += suffix;
  }

  /////////////////////////////////////////////////////////////
  // Configure di-jet task
  AliAnalysisTaskEmcalDijetImbalance* dijetTask = new AliAnalysisTaskEmcalDijetImbalance(name);
  dijetTask->SetDeltaPhiCut(deltaPhiMin);
  if (doGeomMatching) dijetTask->SetDoGeometricalMatching(doGeomMatching, jetR, minTrPtHardCore, minClPtHardCore);
  
  /////////////////////////////////////////////////////////////
  // Create track and cluster containers with the standard cuts
  
  AliParticleContainer* partCont = 0;
  if (trackName == "mcparticles") {
    AliMCParticleContainer* mcpartCont = new AliMCParticleContainer(trackName);
    partCont = mcpartCont;
  }
  else if (trackName == "tracks" || trackName == "Tracks") {
    AliTrackContainer* trackCont = new AliTrackContainer(trackName);
    partCont = trackCont;
  }
  if (partCont) partCont->SetParticlePtCut(minTrPt);
  if (partCont) dijetTask->AdoptParticleContainer(partCont);
  
  AliClusterContainer* clusCont = 0;
  if (!clusName.IsNull()) {
    clusCont = new AliClusterContainer(clusName);
    clusCont->SetClusECut(0.);
    clusCont->SetClusPtCut(0.);
    clusCont->SetClusNonLinCorrEnergyCut(0.);
    clusCont->SetClusHadCorrEnergyCut(minClPt);
    clusCont->SetDefaultClusterEnergy(AliVCluster::kHadCorr);
    if (includePHOS) {
      clusCont->SetIncludePHOS(kTRUE);
      clusCont->SetPhosMinNcells(3);
      clusCont->SetPhosMinM02(0.2);
    }
  }
  if (clusCont) dijetTask->AdoptClusterContainer(clusCont);
  
  /////////////////////////////////////////////////////////////
  // Create track and cluster containers for constituent study with geometrical matching (if enabled)
  
  if (doGeomMatching) {
    AliParticleContainer* partContThresh = 0;
    if (trackName == "mcparticles") {
      AliMCParticleContainer* mcpartCont = new AliMCParticleContainer(trackName);
      partContThresh = mcpartCont;
    }
    else if (trackName == "tracks" || trackName == "Tracks") {
      AliTrackContainer* trackCont = new AliTrackContainer(trackName);
      partContThresh = trackCont;
    }
    if (partContThresh) {
      partContThresh->SetParticlePtCut(minTrPtHardCore);
      partContThresh->SetName("tracksThresh");
      dijetTask->AdoptParticleContainer(partContThresh);
    }
  
    AliClusterContainer* clusContThresh = 0;
    if (!clusName.IsNull()) {
      clusContThresh = new AliClusterContainer(clusName);
      clusContThresh->SetName("caloClustersThresh");
      clusContThresh->SetClusECut(0.);
      clusContThresh->SetClusPtCut(0.);
      clusContThresh->SetClusNonLinCorrEnergyCut(0.);
      clusContThresh->SetClusHadCorrEnergyCut(minClPtHardCore);
      clusContThresh->SetDefaultClusterEnergy(AliVCluster::kHadCorr);
      if (includePHOS) {
        clusContThresh->SetIncludePHOS(kTRUE);
        clusContThresh->SetPhosMinNcells(3);
        clusContThresh->SetPhosMinM02(0.2);
      }
    }
    if (clusContThresh) dijetTask->AdoptClusterContainer(clusContThresh);
  }

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------

  mgr->AddTask(dijetTask);

  // Create containers for input/output
  AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer()  ;
  TString contname(name);
  contname += "_histos";
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(contname.Data(),
      TList::Class(),AliAnalysisManager::kOutputContainer,
      Form("%s", AliAnalysisManager::GetCommonFileName()));
  mgr->ConnectInput  (dijetTask, 0,  cinput1 );
  mgr->ConnectOutput (dijetTask, 1, coutput1 );

  return dijetTask;
}
