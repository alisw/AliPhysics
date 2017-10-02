AliAnalysisTaskEmcalJetHCorrelations* AddTaskEmcalJetHCorrelations(
   const char *nTracks              = "usedefault",
   const char *nCaloClusters        = "usedefault",
   // Jet options
   const Double_t trackBias         = 5,
   const Double_t clusterBias       = 5,
   const Double_t minJetArea        = 0.4,
   // Mixed event options
   const Int_t nTracksMixedEvent    = 0,  // Additionally acts as a switch for enabling mixed events
   const Int_t minNTracksMixedEvent = 5000,
   const Int_t minNEventsMixedEvent = 5,
   const UInt_t nCentBinsMixedEvent = 10,
   // Triggers
   UInt_t trigEvent                 = AliVEvent::kAny,
   UInt_t mixEvent                  = AliVEvent::kAny,
   // Options
   const char *CentEst              = "V0M",
   const Int_t nCentBins            = 5,
   const Double_t trackEta          = 0.9,
   const Bool_t lessSparseAxes      = 0,
   const Bool_t widerTrackBin       = 0,
   // Corrections
   const Int_t doEffCorrSW          = 0,
   const Bool_t embeddingCorrection = kFALSE,
   const char * embeddingCorrectionFilename = "alien:///alice/cern.ch/user/r/rehlersi/embeddingCorrection.root",
   const char * embeddingCorrectionHistName = "embeddingCorrection",
   const char *suffix               = "biased"
)
{  
  AliAnalysisTaskEmcalJetHCorrelations * task = AliAnalysisTaskEmcalJetHCorrelations::AddTaskEmcalJetHCorrelations(
                          nTracks, nCaloClusters,
                          trackBias, clusterBias, minJetArea,
                          nTracksMixedEvent, minNTracksMixedEvent, minNEventsMixedEvent,
                          nCentBinsMixedEvent,
                          trigEvent, mixEvent,
                          CentEst, nCentBins,
                          trackEta,
                          lessSparseAxes,
                          widerTrackBin,
                          doEffCorrSW,
                          embeddingCorrection,
                          embeddingCorrectionFilename, embeddingCorrectionHistName,
                          suffix
                          );
  return task;
}
