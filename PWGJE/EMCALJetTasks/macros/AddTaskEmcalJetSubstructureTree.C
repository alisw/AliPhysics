EmcalTriggerJets::AliAnalysisTaskEmcalJetSubstructureTree *AddTaskEmcalJetSubstructureTree(Bool_t isMC, Bool_t isData, Double_t jetradius, const char *trigger) {
  return EmcalTriggerJets::AliAnalysisTaskEmcalJetSubstructureTree::AddEmcalJetSubstructureTreeMaker(isMC, isData, jetradius, trigger);
}
