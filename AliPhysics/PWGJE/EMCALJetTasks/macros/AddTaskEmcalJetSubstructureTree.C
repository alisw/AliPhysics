PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalJetSubstructureTree *AddTaskEmcalJetSubstructureTree(Bool_t isMC, Bool_t isData, Double_t jetradius, AliJetContainer::EJetType_t jettype, AliJetContainer::ERecoScheme_t recoscheme, Bool_t useDCAL, const char *trigger) {
  return PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalJetSubstructureTree::AddEmcalJetSubstructureTreeMaker(isMC, isData, jetradius, jettype, recoscheme, useDCAL, trigger);
}
