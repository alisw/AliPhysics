EmcalTriggerJets::AliAnalysisTaskEmcalJetEnergyScale *AddTaskEmcalJetEnergyScale(AliJetContainer::EJetType_t jettype, double jetradius, Bool_t useDCAL, const char *trigger){
  return EmcalTriggerJets::AliAnalysisTaskEmcalJetEnergyScale::AddTaskJetEnergyScale(jettype, jetradius, useDCAL, trigger);
}
