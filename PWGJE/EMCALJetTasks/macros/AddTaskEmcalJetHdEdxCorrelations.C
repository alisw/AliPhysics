PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalJetHdEdxCorrelations *AddTaskEmcalJetHdEdxCorrelations(
    const char *nTracks = "usedefault",
    const char *nCaloClusters = "usedefault",
    // Jet options
    const Double_t trackBias = 5,
    const Double_t clusterBias = 5,
    // Mixed event options
    const Int_t nTracksMixedEvent = 0, // Additionally acts as a switch for enabling mixed events
    const Int_t minNTracksMixedEvent = 5000,
    const Int_t minNEventsMixedEvent = 5,
    const UInt_t nCentBinsMixedEvent = 10,
    // Triggers
    UInt_t trigEvent = AliVEvent::kAny,
    UInt_t mixEvent = AliVEvent::kAny,
    // Options
    const Bool_t lessSparseAxes = kFALSE,
    const Bool_t widerTrackBin = kFALSE,
    // Corrections
    const Bool_t embeddingCorrection = kFALSE,
    const char *embeddingCorrectionFilename = "alien:///alice/cern.ch/user/r/rehlersi/embeddingCorrection.root",
    const char *embeddingCorrectionHistName = "embeddingCorrection",
    Bool_t makeEPcorrectionsFile = kFALSE,
    Bool_t doEPshiftCorrection = kFALSE,
    const char *epCorrectionsFilename = "",
    const char *suffix = "biased")
{
  PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalJetHdEdxCorrelations *task = PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalJetHdEdxCorrelations::AddTaskEmcalJetHdEdxCorrelations(
      nTracks, nCaloClusters,
      trackBias, clusterBias,
      nTracksMixedEvent, minNTracksMixedEvent, minNEventsMixedEvent,
      nCentBinsMixedEvent,
      trigEvent, mixEvent,
      lessSparseAxes,
      widerTrackBin,
      embeddingCorrection,
      embeddingCorrectionFilename,
      embeddingCorrectionHistName,
      makeEPcorrectionsFile,
      doEPshiftCorrection,
      epCorrectionsFilename,
          suffix);
  return task;
}
