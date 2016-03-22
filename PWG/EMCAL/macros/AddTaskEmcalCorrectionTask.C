AliEmcalCorrectionTask* AddTaskEmcalCorrectionTask(
  const char *nCells                         = "usedefault",
  const char *nClusters                      = "usedefault",
  const char *nTracks                        = "usedefault",
  const Double_t minTrPt                     = 0.15,
  const Double_t minClPt                     = 0.30
)
{  
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskEmcalCorrectionTask", "No analysis manager to connect to.");
    return 0;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  AliVEventHandler* handler = mgr->GetInputEventHandler();
  if (!handler)
  {
    ::Error("AddTaskEmcalCorrectionTask", "This task requires an input event handler");
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

  TString cellsName(nCells);
  TString clusName(nClusters);
  TString trackName(nTracks);

  if (cellsName == "usedefault") {
    if (dataType == kESD) {
      cellsName = "EMCALCells";
    }
    else if (dataType == kAOD) {
      cellsName = "emcalCells";
    }
    else {
      cellsName = "";
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

  AliParticleContainer* partCont = 0;
  if (trackName == "mcparticles") {
    AliMCParticleContainer* mcpartCont = new AliMCParticleContainer(trackName);
    mcpartCont->SelectPhysicalPrimaries(kTRUE);
    partCont = mcpartCont;
  }
  else if (trackName == "tracks" || trackName == "Tracks") {
    AliTrackContainer* trackCont = new AliTrackContainer(trackName);
    trackCont->SetFilterHybridTracks(kTRUE);
    partCont = trackCont;
  }
  else if (!trackName.IsNull()) {
    partCont = new AliParticleContainer(trackName);
  }
  if (partCont) partCont->SetParticlePtCut(minTrPt);

  AliClusterContainer* clusCont = 0;
  if (!clusName.IsNull()) {
    clusCont = new AliClusterContainer(clusName);
    clusCont->SetClusECut(0.);
    clusCont->SetClusPtCut(0.);
    //clusCont->SetClusHadCorrEnergyCut(minClPt);
    //clusCont->SetDefaultClusterEnergy(AliVCluster::kHadCorr);
  }

  TString name = "AliEmcalCorrectionTask";
  
  Printf("name: %s", name.Data());

  AliEmcalCorrectionTask* mgrTask = static_cast<AliEmcalCorrectionTask *>(mgr->GetTask(name.Data()));
  if (mgrTask) return mgrTask;
  
  Printf("name past return: %s", name.Data());

  // Create the task that manages
  AliEmcalCorrectionTask* correctionTask = new AliEmcalCorrectionTask(name.Data());
  if (!cellsName.IsNull()) correctionTask->SetCaloCellsName(cellsName);
  if (clusCont) correctionTask->AdoptClusterContainer(clusCont);
  if (partCont) correctionTask->AdoptParticleContainer(partCont);
  correctionTask->SetUserConfigurationFilename("");
  correctionTask->SetDefaultConfigurationFilename("AliEmcalCorrectionConfiguration.yaml");
  //correctionTask->InitializeConfiguration();

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------

  mgr->AddTask(correctionTask);

  Printf("Creating i/o containers");
  // Create containers for input/output
  AliAnalysisDataContainer* cInput = mgr->GetCommonInputContainer();
  
  TString outputContainerName(name);
  outputContainerName += "_histos";
  
  AliAnalysisDataContainer * cOutput = mgr->CreateContainer(outputContainerName.Data(),
                                                            TList::Class(),
                                                            AliAnalysisManager::kOutputContainer,
                                                            Form("%s", AliAnalysisManager::GetCommonFileName()));
  
  mgr->ConnectInput(correctionTask, 0, cInput);
  mgr->ConnectOutput(correctionTask, 1, cOutput);
  
  //TObjArray* cnt = mgr->GetContainers();

  return correctionTask;
}
