EmcalTriggerJets::AliAnalysisTaskEmcalJetEnergySpectrum *AddTaskJetEnergySpectrum(Bool_t isMC, AliJetContainer::EJetType_t jettype, Double_t radius, const char *trigger, const char *suffix){
  return EmcalTriggerJets::AliAnalysisTaskEmcalJetEnergySpectrum::AddTaskEmcalJetEnergySpectrum(isMC, jettype, radius, trigger, suffix);
}