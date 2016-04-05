AliAnalysisTaskEmcalJetHMEC* AddTaskEmcalJetHMEC(
   const char *nTracks        = "usedefault",
   const char *nCaloClusters  = "usedefault",
   // Jet options
   const Double_t trackBias     = 5,
   const Double_t clusterBias   = 5,
   const Double_t minJetArea     = 0.4,
   // Mixed event options
   const Int_t nTracksMixedEvent = 0,  // Additionally acts as a switch for enabling mixed events
   const Int_t minNTracksMixedEvent = 5000,
   const Int_t minNEventsMixedEvent = 5,
   const UInt_t nCentBinsMixedEvent = 10,
   // Triggers
   UInt_t trigEvent           = AliVEvent::kAny,
   UInt_t mixEvent            = AliVEvent::kAny,
   // Options
   const char *suffix          = "biased",
   const char *CentEst         = "V0M",
   const Int_t nCentBins       = 5,
   const Double_t trackEta     = 0.9,
   const Bool_t lessSparseAxes = 0,
   const Bool_t widerTrackBin  = 0,
   // Corrections
   const Int_t doEffCorrSW     = 0,
   const Bool_t embeddingCorrection = kFALSE,
   const char * embeddingCorrectionFilename = "alien:///alice/cern.ch/user/r/rehlersi/embeddingCorrection.root",
   const char * embeddingCorrectionHistName = "embeddingCorrection",
   // Beam type
   const Short_t beamType      = AliAnalysisTaskEmcal::kAA, 
)
{  
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskEmcalJetHMEC", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  AliVEventHandler* handler = mgr->GetInputEventHandler();
  if (!handler)
  {
    ::Error("AddTaskEmcalJetHMEC", "This task requires an input event handler");
    return NULL;
  }
  
  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------

  // Determine data type
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
  
  // Determine cluster and track names
  TString trackName(nTracks);
  TString clusName(nCaloClusters);

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

  TString name("AliAnalysisTaskJetH");
  if (!trackName.IsNull()) {
    name += TString::Format("_%s", trackName.Data());
  }
  if (!clusName.IsNull()) {
    name += TString::Format("_%s", clusName.Data());
  }
  if (strcmp(suffix, "") != 0) {
    name += TString::Format("_%s", suffix);
  }

  AliAnalysisTaskEmcalJetHMEC *correlationTask = new AliAnalysisTaskEmcalJetHMEC(name);
  // Set jet bias
  correlationTask->SetTrackBias(trackBias);
  correlationTask->SetClusterBias(clusterBias);
  // Mixed events
  correlationTask->SetEventMixing(static_cast<Bool_t>(nTracksMixedEvent));
  correlationTask->SetNumberOfMixingTracks(nTracksMixedEvent);
  correlationTask->SetMinNTracksForMixedEvents(minNTracksMixedEvent);
  correlationTask->SetMinNEventsForMixedEvents(minNEventsMixedEvent);
  correlationTask->SetNCentBinsMixedEvent(nCentBinsMixedEvent);
  // Triggers
  correlationTask->SetTriggerType(trigEvent);
  correlationTask->SetMixedEventTriggerType(mixEvent);
  // Options
  correlationTask->SetCentralityEstimator(CentEst);
  correlationTask->SetNCentBins(nCentBins);
  correlationTask->SetForceBeamType(beamType);
  correlationTask->SetVzRange(-10,10);
  correlationTask->SetDoLessSparseAxes(lessSparseAxes);
  correlationTask->SetDoWiderTrackBin(widerTrackBin);
  // Corrections
  correlationTask->SetDoEffCorr(doEffCorrSW);
  if (embeddingCorrection == kTRUE)
  {
    // Open file containing the correction
    TFile * embeddingCorrectionFile = TFile::Open(embeddingCorrectionFilename);
    if (!embeddingCorrectionFile || embeddingCorrectionFile->IsZombie()) {
        ::Error("AddTaskEmcalJetHMEC", Form("Could not open embedding correction file %s", embeddingCorrectionFilename));
        return NULL;
    }

    // Retrieve the histogram containing the correction and save add it to the task.
    TH2F * embeddingCorrectionHist = dynamic_cast<TH2F*>(file->Get(embeddingCorrectionHistName));
    if (embeddingCorrectionHist) {
      ::Info("AddTaskEmcalJetHMEC", Form("Embedding correction %s loaded from file %s.", embeddingCorrectionHistName, embeddingCorrectionFilename));
    }
    else {
      ::Error("AddTaskEmcalJetHMEC", Form("Embedding correction %s not found in file %s.", embeddingCorrectionHistName, embeddingCorrectionFilename));
      return NULL;
    }

    correlationTask->SetEmbeddingCorrectionHist(embeddingCorrectionHist);
  }

  // Jet parameters determined by how we ran the jet finder
  Double_t jetRadius = 0.2;
  Double_t minClusterPt = 3;
  Double_t minTrackPt = 3;

  // Add Containers
  // Clusters
  AliClusterContainer * clusterContainer = correlationTask->AddClusterContainer(clusName);
  clusterContainer->SetMinE(minClusterPt);

  // Tracks
  // For jet finding
  AliTrackContainer * tracksForJets = new AliTrackContainer(trackName);
  tracksForJets->SetName("tracksForJets");
  tracksForJets->SetMinPt(minTrackPt);
  tracksForJets->SetEtaLimits(-1.0*trackEta, trackEta);
  // Adopt the container
  correlationTask->AdoptParticleContainer(tracksForJets);
  // For correlations
  AliTrackContainer * tracksForCorrelations = new AliTrackContainer(trackName);
  tracksForCorrelations->SetName("tracksForCorrelations");
  tracksForCorrelations->SetMinPt(0.15);
  tracksForCorrelations->SetEtaLimits(-1.0*trackEta, trackEta);
  // Adopt the container
  correlationTask->AdoptParticleContainer(tracksForCorrelations);

  // Jets
  AliJetContainer * jetContainer = correlationTask->AddJetContainer(AliJetContainer::kFullJet,
                                   AliJetContainer::antikt_algorithm,
                                   AliJetContainer::pt_scheme,
                                   jetRadius,
                                   AliJetContainer::kEMCALfid,
                                   tracksForJets,
                                   clusterContainer);
  jetContainer->SetJetAreaCut(minJetArea);
  jetContainer->SetMaxTrackPt(100);
  jetContainer->SetJetPtCut(0.1);

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------

  mgr->AddTask(correlationTask);

  // Create containers for input/output
  mgr->ConnectInput (correlationTask, 0, mgr->GetCommonInputContainer() );
  AliAnalysisDataContainer * cojeth = mgr->CreateContainer(name,
                TList::Class(),
                AliAnalysisManager::kOutputContainer,
							    Form("%s", AliAnalysisManager::GetCommonFileName()));
  mgr->ConnectOutput(correlationTask, 1, cojeth);

  return correlationTask;
}
