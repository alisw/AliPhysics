// AddTaskEmcalJetTree.C

AliAnalysisTaskEmcalJetTreeBase* AddTaskEmcalJetTree(
    const char *ntracks            = "usedefault",
    const char *nclusters          = "usedefault",
    Double_t    trackPtCut         = 0.15,
    Double_t    clusECut           = 0.30,
    AliAnalysisTaskEmcalJetTreeBase::EAnalisysType_t type = AliAnalysisTaskEmcalJetTreeBase::kJetPP,
    const char *suffix             = ""
)
{  
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskEmcalJetTree", "No analysis manager to connect to.");
    return 0;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  AliVEventHandler* handler = mgr->GetInputEventHandler();
  if (!handler)
  {
    ::Error("AddTaskEmcalJetTree", "This task requires an input event handler");
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

  TString name("AliAnalysisTaskEmcalJetTree");
  if (strcmp(suffix,"")) {
    name += "_";
    name += suffix;
  }

  AliAnalysisTaskEmcalJetTreeBase* jetTask = AliAnalysisTaskEmcalJetTreeBase::CreateInstance(name, type);
  jetTask->SetVzRange(-10,10);
  jetTask->SetNeedEmcalGeom(kFALSE);

  if (trackName == "mcparticles") {
    AliMCParticleContainer* mcpartCont = jetTask->AddMCParticleContainer(trackName);
    mcpartCont->SelectPhysicalPrimaries(kTRUE);
  }
  else if (trackName == "tracks" || trackName == "Tracks") {
    AliTrackContainer* trackCont = jetTask->AddTrackContainer(trackName);
  }
  else if (!trackName.IsNull()) {
    jetTask->AddParticleContainer(trackName);
  }

  AliParticleContainer *partCont = jetTask->GetParticleContainer(0);
  if (partCont) {
    partCont->SetParticlePtCut(trackPtCut);
  }

  AliClusterContainer *clusterCont = jetTask->AddClusterContainer(clusName);
  if (clusterCont) {
    clusterCont->SetClusECut(0.);
    clusterCont->SetClusPtCut(0.);
    clusterCont->SetClusHadCorrEnergyCut(clusECut);
    clusterCont->SetDefaultClusterEnergy(AliVCluster::kHadCorr);
  }

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------
  mgr->AddTask(jetTask);
  
  // Create containers for input/output
  AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer()  ;
  TString contname0(name);
  contname0 += "_jets";
  TString contname1(name);
  contname1 += "_histos";
  AliAnalysisDataContainer *coutput0 = mgr->CreateContainer(contname0.Data(),
                  TTree::Class(),AliAnalysisManager::kOutputContainer,
                  Form("%s", AliAnalysisManager::GetCommonFileName()));
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(contname1.Data(),
							    TList::Class(),AliAnalysisManager::kOutputContainer,
							    Form("%s", AliAnalysisManager::GetCommonFileName()));
  mgr->ConnectInput  (jetTask, 0,  cinput1 );
  mgr->ConnectOutput (jetTask, 1, coutput1 );
  mgr->ConnectOutput (jetTask, 2, coutput0 );

  return jetTask;
}
