PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalJetEnergyScale *AddTaskEmcalJetEnergyScale(AliJetContainer::EJetType_t jettype, AliJetContainer::ERecoScheme_t recoscheme, double jetradius, Bool_t useDCAL, const char *namepartcont, const char *trigger, const char *suffix){
  return PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalJetEnergyScale::AddTaskJetEnergyScale(jettype, recoscheme, jetradius, useDCAL, namepartcont, trigger, suffix);
}
