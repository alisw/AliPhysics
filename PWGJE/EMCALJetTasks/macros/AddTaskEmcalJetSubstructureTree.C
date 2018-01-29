EmcalTriggerJets::AliAnalysisTaskEmcalJetSubstructureTree *AddTaskEmcalJetSubstructureTree(Bool_t isMC, Bool_t isData, Double_t jetradius, AliJetContainer::EJetType_t jettype, AliJetContainer::ERecoScheme_t recoscheme, const char *trigger) {
  return EmcalTriggerJets::AliAnalysisTaskEmcalJetSubstructureTree::AddEmcalJetSubstructureTreeMaker(isMC, isData, jetradius, jettype, recoscheme, trigger);
}
