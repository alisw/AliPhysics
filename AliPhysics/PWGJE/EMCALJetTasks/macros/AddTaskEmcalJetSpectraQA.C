// AddTaskEmcalJetSpectraQA.C

AliAnalysisTaskEmcalJetSpectraQA* AddTaskEmcalJetSpectraQA(
  const char *ntracks            = "usedefault",
  const char *nclusters          = "usedefault",
  Double_t    trackPtCut         = 0.15,
  Double_t    clusECut           = 0.30,
  const char *suffix             = ""
)
{
  return AliAnalysisTaskEmcalJetSpectraQA::AddTaskEmcalJetSpectraQA(ntracks, nclusters, trackPtCut, clusECut, suffix);
}
