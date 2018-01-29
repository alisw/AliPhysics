AliAnalysisTaskEmcalJetPerformance* AddTaskEmcalJetPerformance(
  const char *ntracks            = "usedefault",
  const char *nclusters          = "usedefault",
  const char *nGenLev            = "mcparticles",
  const Double_t minTrPt         = 0.15,              // Minimum track pT in standard track container
  const Double_t minClPt         = 0.30,              // Minimum cluster E in standard cluster container
  const char *suffix             = "")
{
  return AliAnalysisTaskEmcalJetPerformance::AddTaskEmcalJetPerformance(ntracks,
                                                                        nclusters,
                                                                        nGenLev,
                                                                        minTrPt,
                                                                        minClPt,
                                                                        suffix);
}
