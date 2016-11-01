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
  return AliAnalysisTaskEmcalJetTreeBase::AddTaskEmcalJetTree(ntracks, nclusters, trackPtCut, clusECut, type, suffix);
}
