// AddTaskEmcalJetSpectraQA.C

AliAnalysisHFjetTagHFE* AddTaskHFjetTagHFE(
  const char *ntracks            = "usedefault",
  const char *nclusters          = "usedefault",
  const char *njets              = "Jets",
  const char *nrho               = "Rho",
  Double_t    jetradius          = 0.2,
  Double_t    jetptcut           = 1,
  Double_t    jetareacut         = 0.557,
  const char *cutType            = "TPCfid",
  Int_t       leadhadtype        = 0,
  const char *suffix             = ""
)
{  
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskHFjetTagHFE", "No analysis manager to connect to.");
    return 0;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  AliVEventHandler* handler = mgr->GetInputEventHandler();
  if (!handler)
  {
    ::Error("AddTaskHFjetTagHFE", "This task requires an input event handler");
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

  TString name("AliAnalysisHFjetTagHFE");
  if (strcmp(njets,"")) {
    name += "_";
    name += njets;
  }
  if (strcmp(nrho,"")) {
    name += "_";
    name += nrho;
  }
  name += "_";
  name += cutType;

  if (strcmp(suffix,"")) {
    name += "_";
    name += suffix;
  }

  AliAnalysisHFjetTagHFE* jetTask = new AliAnalysisHFjetTagHFE(name);
  jetTask->SetVzRange(-10,10);
  jetTask->SetNeedEmcalGeom(kFALSE);

  /*
  AliParticleContainer *trackCont = jetTask->AddTrackContainer(trackName);
  AliClusterContainer *clusterCont = jetTask->AddClusterContainer(clusName);
  */

  if (trackName == "mcparticles") {
    AliMCParticleContainer* mcpartCont = jetTask->AddMCParticleContainer(trackName);
    mcpartCont->SelectPhysicalPrimaries(kTRUE);
  }
  else if (trackName == "tracks" || trackName == "Tracks") {
    AliTrackContainer* trackCont = jetTask->AddTrackContainer(trackName);
    trackCont->SetFilterHybridTracks(kTRUE);
  }
  else if (!trackName.IsNull()) {
    jetTask->AddParticleContainer(trackName);
  }

 /*
  AliParticleContainer *partCont = jetTask->GetParticleContainer(0);
  if (partCont) {
    partCont->SetParticlePtCut(trackPtCut);
  }
*/ 

  AliClusterContainer *clusterCont = jetTask->AddClusterContainer(clusName);
  if (clusterCont) {
    clusterCont->SetClusECut(0.);
    clusterCont->SetClusPtCut(0.);
    clusterCont->SetClusHadCorrEnergyCut(clusECut);
    clusterCont->SetDefaultClusterEnergy(AliVCluster::kHadCorr);
  }

  AliJetContainer *jetCont = jetTask->AddJetContainer(njets, cutType, jetradius);
  if (jetCont) {
    jetCont->SetRhoName(nrho);
    jetCont->SetPercAreaCut(jetareacut);
    jetCont->SetJetPtCut(jetptcut);
    jetCont->ConnectParticleContainer(trackCont);
    jetCont->ConnectClusterContainer(clusterCont);
    jetCont->SetLeadingHadronType(leadhadtype);
    jetCont->SetMaxTrackPt(1000);
    jetCont->SetZLeadingCut(0.98,0.98);
  }

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------
  
  mgr->AddTask(jetTask);
  
  // Create containers for input/output
  AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer()  ;
  TString contname(name);
  contname += "_histos";
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(contname.Data(), 
							    TList::Class(),AliAnalysisManager::kOutputContainer,
							    Form("%s", AliAnalysisManager::GetCommonFileName()));
  mgr->ConnectInput  (jetTask, 0,  cinput1 );
  mgr->ConnectOutput (jetTask, 1, coutput1 );
  
  return jetTask;
}
