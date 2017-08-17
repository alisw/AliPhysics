EmcalTriggerJets::AliAnalysisTaskEmcalJetSubstructureTree *AddTaskEmcalJetSubstructureTree(Bool_t isMC, Bool_t isData, Double_t jetradius, AliJetContainer::ERecoScheme_t recoscheme, const char *trigger) {
  return EmcalTriggerJets::AliAnalysisTaskEmcalJetSubstructureTree::AddEmcalJetSubstructureTreeMaker(isMC, isData, jetradius, recoscheme, trigger);
}
