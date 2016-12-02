// AddTaskPWGJEQA.C

AliAnalysisTaskPWGJEQA* AddTaskPWGJEQA(
                                       const char* ntracks            = "usedefault",
                                       const char* nclusters          = "usedefault",
                                       const char* ncells             = "usedefault",
                                       const char *nGenLev            = "mcparticles",
                                       Bool_t      doTrackQA          = kTRUE,
                                       Bool_t      doEmcalQA          = kTRUE,
                                       Bool_t      doJetQA            = kTRUE,
                                       Bool_t      doEventQA          = kTRUE,
                                       Double_t    trackPtCut         = 0.15,
                                       Double_t    clusECut           = 0.30,
                                       const char* suffix             = ""
                                       )
{
  
  // Get the pointer to the existing analysis manager via the static access method
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskPWGJEQA", "No analysis manager to connect to.");
    return NULL;
  }
  
  // Check the analysis type using the event handlers connected to the analysis manager
  AliVEventHandler* handler = mgr->GetInputEventHandler();
  if (!handler) {
    ::Error("AddTaskPWGJEQA", "This task requires an input event handler");
    return NULL;
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
  
  // Init the task and do settings
  TString trackName(ntracks);
  TString clusName(nclusters);
  TString cellName(ncells);
  
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
  
  if (ncells == "usedefault") {
    if (dataType == kESD) {
      cellName = "EMCALCells";
    }
    else if (dataType == kAOD) {
      cellName = "emcalCells";
    }
    else {
      cellName = "";
    }
  }
  
  TString name("AliAnalysisTaskPWGJEQA");
  if (!trackName.IsNull()) {
    name += "_";
    name += trackName;
  }
  if (!clusName.IsNull()) {
    name += "_";
    name += clusName;
  }
  
  if (!cellName.IsNull()) {
    name += "_";
    name += cellName;
  }
   
  if (strcmp(suffix,"")) {
    name += "_";
    name += suffix;
  }
  
  AliAnalysisTaskPWGJEQA* qaTask = new AliAnalysisTaskPWGJEQA(name);
  qaTask->SetVzRange(-10,10);
  qaTask->SetNeedEmcalGeom(kFALSE);
  qaTask->SetCaloCellsName(cellName);
  qaTask->SetDetectorLevelName(trackName);
  if (nGenLev && strcmp(nGenLev,"")!=0) qaTask->SetGeneratorLevelName(nGenLev);
  qaTask->SetDoTrackQA(doTrackQA);
  qaTask->SetDoEmcalQA(doEmcalQA);
  qaTask->SetDoJetQA(doJetQA);
  qaTask->SetDoEventQA(doEventQA);
  
  // Add the detector-level track container
  if (trackName == "tracks" || trackName == "Tracks") {
    AliTrackContainer* trackCont = qaTask->AddTrackContainer(trackName);
    trackCont->SetFilterHybridTracks(kTRUE);
  }
  else if (!trackName.IsNull()) {
    qaTask->AddParticleContainer(trackName);
  }
  AliParticleContainer *partCont = qaTask->GetParticleContainer(trackName);
  if (partCont) {
    partCont->SetParticlePtCut(trackPtCut);
  }
  
  // Add the generator-level container, if specified
  if (nGenLev && strcmp(nGenLev,"")!=0) {
    AliMCParticleContainer* mcpartCont = qaTask->AddMCParticleContainer(nGenLev);
    mcpartCont->SelectPhysicalPrimaries(kTRUE);
    mcpartCont->SetParticlePtCut(0);
  }

  // Add the cluster container
  AliClusterContainer *clusterCont = qaTask->AddClusterContainer(clusName);
  if (clusterCont) {
    clusterCont->SetClusECut(0.);
    clusterCont->SetClusPtCut(0.);
    clusterCont->SetClusHadCorrEnergyCut(clusECut);
    clusterCont->SetDefaultClusterEnergy(AliVCluster::kHadCorr);
  }

  // Final settings, pass to manager and set the containers
  mgr->AddTask(qaTask);
  
  // Create containers for input/output
  AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer()  ;
  
  TString contName = TString::Format("%s_histos", name.Data());
  TString commonoutput;
  commonoutput = mgr->GetCommonFileName();
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(contName.Data(),
                                                            TList::Class(),AliAnalysisManager::kOutputContainer,
                                                            commonoutput);
  mgr->ConnectInput  (qaTask, 0,  cinput1 );
  mgr->ConnectOutput (qaTask, 1, coutput1 );
  
  return qaTask;
}
