EmcalTriggerJets::AliAnalysisTaskEmcalJetEnergyScale *AddTaskEmcalJetEnergyScale(AliJetContainer::EJetType_t jettype, double jetradius, Bool_t useDCAL, const char *namepartcont, const char *trigger){
  return EmcalTriggerJets::AliAnalysisTaskEmcalJetEnergyScale::AddTaskJetEnergyScale(jettype, jetradius, useDCAL, namepartcont, trigger);
}
