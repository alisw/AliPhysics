PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalJetEnergyScale *AddTaskEmcalJetEnergyScaleBkgSub(AliJetContainer::EJetType_t jettype, AliJetContainer::ERecoScheme_t recoscheme, AliVCluster::VCluUserDefEnergy_t energydef, double jetradius, Bool_t useDCAL, const char *namepartcont, const char *nRho, const char *nRhoMC, const char *trigger, const char *suffix){
  return PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalJetEnergyScale::AddTaskJetEnergyScaleBkgSub(jettype, recoscheme, energydef, jetradius, useDCAL, namepartcont, nRho, nRhoMC, trigger, suffix);
}

PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalJetEnergyScale *AddTaskEmcalJetEnergyScaleBkgSub(AliJetContainer::EJetType_t jettype, AliJetContainer::ERecoScheme_t recoscheme, AliVCluster::VCluUserDefEnergy_t energydef, double jetradius, Bool_t useDCAL, const char *namepartcont, const char *nRho, const char *nRhoMC, const char *trigger, const char *nametrackcont, const char *nameclustercont, const char *suffix){
  return PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalJetEnergyScale::AddTaskJetEnergyScaleBkgSub(jettype, recoscheme, energydef, jetradius, useDCAL, namepartcont, nRho, nRhoMC, trigger, nametrackcont, nameclustercont, suffix);
}
