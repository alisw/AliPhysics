PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalJetConstituentQA *AddTaskEmcalJetConstituentQA(AliJetContainer::EJetType_t jettype, bool part, const char *trigger) {
  return PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalJetConstituentQA::AddTaskEmcalJetConstituentQA(trigger, jettype, part);
}
