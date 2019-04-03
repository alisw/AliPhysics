EmcalTriggerJets::AliAnalysisTaskEmcalJetEnergySpectrum *AddTaskEmcalJetEnergySpectrum(Bool_t isMC, AliJetContainer::EJetType_t jettype, AliJetContainer::ERecoScheme_t recoscheme, Double_t radius, const char *namepartcont, const char *trigger){
  return EmcalTriggerJets::AliAnalysisTaskEmcalJetEnergySpectrum::AddTaskJetEnergySpectrum(isMC, jettype, recoscheme, radius, namepartcont, trigger);
}
