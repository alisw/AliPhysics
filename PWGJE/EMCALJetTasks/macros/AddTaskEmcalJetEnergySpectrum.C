PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalJetEnergySpectrum *AddTaskEmcalJetEnergySpectrum(Bool_t isMC, AliJetContainer::EJetType_t jettype, AliJetContainer::ERecoScheme_t recoscheme, Double_t radius, const char *namepartcont, const char *trigger, const char *suffix){
  return PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalJetEnergySpectrum::AddTaskJetEnergySpectrum(isMC, jettype, recoscheme, radius, namepartcont, trigger, suffix);
}
