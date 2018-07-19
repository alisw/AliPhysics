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
  return AliAnalysisTaskEmcalDijetImbalance::AddTaskEmcalDijetImbalance(ntracks,
                                                                        nclusters,
                                                                        deltaPhiMin,
                                                                        doGeomMatching,
                                                                        minTrPtHardCore,
                                                                        minClPtHardCore,
                                                                        jetR,
                                                                        includePHOS,
                                                                        minTrPt,
                                                                        minClPt,
                                                                        suffix);
}
