AliAnalysisTaskRho* AddTaskRhoNew(
   const char    *nTracks                        = "usedefault",
   const char    *nClusters                      = "usedefault",
   const char    *nRho                           = "Rho",
   Double_t       jetradius                      = 0.2,
   AliJetContainer::JetAcceptanceType acceptance = AliJetContainer::kTPCfid,
   AliJetContainer::EJetType_t jetType           = AliJetContainer::kChargedJet,
   const Bool_t   histo                          = kFALSE,
   AliJetContainer::ERecoScheme_t rscheme = AliJetContainer::pt_scheme,
   const char    *suffix                         = ""
)
{  
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskRho", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  AliVEventHandler* handler = mgr->GetInputEventHandler();
  if (!handler)
  {
    ::Error("AddTaskEmcalJetSpectraQA", "This task requires an input event handler");
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

  TString name("AliAnalysisTaskRho");
  if (strcmp(suffix,"") != 0) {
    name += "_";
    name += suffix;
  }

  AliAnalysisTaskRho* mgrTask = mgr->GetTask(name.Data());
  if (mgrTask) return mgrTask;

  AliAnalysisTaskRho *rhotask = new AliAnalysisTaskRho(name, histo);
  rhotask->SetOutRhoName(nRho);

  if (trackName == "mcparticles") {
    AliMCParticleContainer* mcpartCont = rhotask->AddMCParticleContainer(trackName);
    mcpartCont->SelectPhysicalPrimaries(kTRUE);
  }
  else if (trackName == "tracks" || trackName == "Tracks") {
    AliTrackContainer* trackCont = rhotask->AddTrackContainer(trackName);
  }
  else if (!trackName.IsNull()) {
    rhotask->AddParticleContainer(trackName);
  }
  AliParticleContainer* partCont = rhotask->GetParticleContainer(0);

  AliClusterContainer *clusterCont = rhotask->AddClusterContainer(clusName);
  if (clusterCont) {
    clusterCont->SetClusECut(0.);
    clusterCont->SetClusPtCut(0.);
    clusterCont->SetClusHadCorrEnergyCut(0.3);
    clusterCont->SetDefaultClusterEnergy(AliVCluster::kHadCorr);
  }

  AliJetContainer *jetCont = rhotask->AddJetContainer(jetType, AliJetContainer::kt_algorithm, rscheme, jetradius, acceptance, partCont, clusterCont);
  if (jetCont) jetCont->SetJetPtCut(0);

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------

  mgr->AddTask(rhotask);

  // Create containers for input/output
  mgr->ConnectInput(rhotask, 0, mgr->GetCommonInputContainer());
  if (histo) {
    TString contname(name);
    contname += "_histos";
    AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(contname.Data(), 
        TList::Class(), AliAnalysisManager::kOutputContainer,
        Form("%s", AliAnalysisManager::GetCommonFileName()));
    mgr->ConnectOutput(rhotask, 1, coutput1);
  }

  return rhotask;
}
