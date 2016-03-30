// $Id$

AliAnalysisTaskEmcalJetHMEC* AddTaskEmcalJetHMEC(
   const char *nTracks        = "PicoTracks",
   const char *nCaloClusters  = "CaloClustersCorr",
   const Double_t minArea     = 0.4,
   const Int_t EvtMix         = 0, 
   const Double_t TrkBias     = 5,
   const Double_t ClusBias    = 5,
   const Double_t TrkEta      = 0.9,
   const Int_t nmixingTR      = 5000,
   const Int_t nmixingEV      = 5,
   UInt_t trigevent           = AliVEvent::kAny,
   UInt_t mixevent            = AliVEvent::kAny,
   Bool_t lessSparseAxes      = 0,
   Bool_t widertrackbin       = 0,
   UInt_t centbinsize         = 1,
   const Int_t doEffcorrSW    = 0,
   const char *suffix         = "biased",
   const char *CentEst        = "V0M",
   const Short_t beamType     = AliAnalysisTaskEmcal::kAA, 
   Bool_t embeddingCorrection = kFALSE,
   const char * embeddingCorrectionFilename = "alien:///alice/cern.ch/user/r/rehlersi/embeddingCorrection.root",
   const char * embeddingCorrectionHistName = "embeddingCorrection"
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

  AliAnalysisTaskEmcalJetHMEC *correlationtask = new AliAnalysisTaskEmcalJetHMEC(name);
  if(EvtMix>0){
    correlationtask->SetMixingTracks(EvtMix);
    correlationtask->SetEventMixing(1);
    correlationtask->SetNMixedTracks(nmixingTR);
    correlationtask->SetNMixedEvents(nmixingEV);
  }else{
    correlationtask->SetEventMixing(EvtMix);
  }
  correlationtask->SetTrkBias(TrkBias);
  correlationtask->SetClusBias(ClusBias);
  correlationtask->SetTrkEta(TrkEta);
  correlationtask->SetTrigType(trigevent);
  correlationtask->SetMixType(mixevent);
  correlationtask->SetDoLessSparseAxes(lessSparseAxes);
  correlationtask->SetDoWiderTrackBin(widertrackbin);
  correlationtask->SetCentBinSize(centbinsize);
  correlationtask->SetDoEffCorr(doEffcorrSW);
  correlationtask->SetCentralityEstimator(CentEst);
  correlationtask->SetForceBeamType(beamType);
  correlationtask->SetVzRange(-10,10);

  // Determined by how we run the jet finder
  Double_t jetRadius = 0.2;
  Double_t minClusterPt = 3;
  Double_t minTrackPt = 3;
  // Add Containers
  // Clusters
  AliClusterContainer * clusterContainer = correlationtask->AddClusterContainer(clusName);
  clusterContainer->SetMinE(minClusterPt);
  // Tracks
  // For jet finding
  AliTrackContainer * tracksForJets = new AliTrackContainer(trackName);
  tracksForJets->SetName("tracksForJets");
  tracksForJets->SetMinPt(minTrackPt);
  tracksForJets->SetEtaLimits(-1.0*TrkEta, TrkEta);
  // Adopt the container
  correlationtask->AdoptTrackContainer(tracksForJets);
    
  // For correlations
  AliTrackContainer * tracksForCorrelations = new AliTrackContainer(trackName);
  tracksForCorrelations->SetName("tracksForCorrelations");
  tracksForCorrelations->SetMinPt(0.15);
  tracksForCorrelations->SetEtaLimits(-1.0*TrkEta, TrkEta);
  // Adopt the container
  correlationtask->AdoptTrackContainer(tracksForCorrelations);

  // Jets
  AliJetContainer * jetContainer = correlationtask->AddJetContainer(AliJetContainer::kFullJet,
                                   AliJetContainer::antikt_algorithm,
                                   AliJetContainer::pt_scheme,
                                   jetRadius,
                                   AliJetContainer::kEMCALfid,
                                   tracksForJets,
                                   clusterContainer);
  jetContainer->SetJetAreaCut(minArea);
  jetContainer->SetMaxTrackPt(100);
  jetContainer->SetJetPtCut(0.1);

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

    correlationtask->SetEmbeddingCorrectionHist(embeddingCorrectionHist);
  }

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------

  mgr->AddTask(correlationtask);

  // Create containers for input/output
  mgr->ConnectInput (correlationtask, 0, mgr->GetCommonInputContainer() );
  AliAnalysisDataContainer *cojeth = mgr->CreateContainer(name,
                TList::Class(),
                AliAnalysisManager::kOutputContainer,
							    Form("%s", AliAnalysisManager::GetCommonFileName()));
  mgr->ConnectOutput(correlationtask,1,cojeth);

  return correlationtask;
}
