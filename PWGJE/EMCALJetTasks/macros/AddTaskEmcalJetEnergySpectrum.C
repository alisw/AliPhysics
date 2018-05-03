EmcalTriggerJets::AliAnalysisTaskEmcalJetEnergySpectrum *AddTaskEmcalJetEnergySpectrum(Bool_t isMC, AliJetContainer::EJetType_t jettype, Double_t radius, const char *trigger, const char *suffix){
  return EmcalTriggerJets::AliAnalysisTaskEmcalJetEnergySpectrum::AddTaskJetEnergySpectrum(isMC, jettype, radius, trigger, suffix);
}
