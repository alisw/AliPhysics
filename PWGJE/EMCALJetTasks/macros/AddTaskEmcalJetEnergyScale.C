PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalJetEnergyScale *AddTaskEmcalJetEnergyScale(AliJetContainer::EJetType_t jettype, AliJetContainer::ERecoScheme_t recoscheme, AliVCluster::VCluUserDefEnergy_t energydef, double jetradius, Bool_t useDCAL, const char *namepartcont, const char *trigger, const char *suffix){
  return PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalJetEnergyScale::AddTaskJetEnergyScale(jettype, recoscheme, energydef, jetradius, useDCAL, namepartcont, trigger, suffix);
}
