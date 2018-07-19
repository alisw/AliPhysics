EmcalTriggerJets::AliAnalysisTaskEmcalJetConstituentQA *AddTaskEmcalJetConstituentQA(AliJetContainer::EJetType_t jettype, bool part, const char *trigger) {
  return EmcalTriggerJets::AliAnalysisTaskEmcalJetConstituentQA::AddTaskEmcalJetConstituentQA(trigger, jettype, part);
}
