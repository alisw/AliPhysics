PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalJetEnergySpectrum *AddTaskEmcalJetEnergySpectrum(Bool_t isMC, AliJetContainer::EJetType_t jettype, AliJetContainer::ERecoScheme_t recoscheme, AliVCluster::VCluUserDefEnergy_t energydef, Double_t radius, const char *namepartcont, const char *trigger, const char *suffix){
  return PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalJetEnergySpectrum::AddTaskJetEnergySpectrum(isMC, jettype, recoscheme, energydef, radius, namepartcont, trigger, suffix);
}

PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalJetEnergySpectrum *AddTaskEmcalJetEnergySpectrum(Bool_t isMC, AliJetContainer::EJetType_t jettype, AliJetContainer::ERecoScheme_t recoscheme, AliVCluster::VCluUserDefEnergy_t energydef, Double_t radius, const char *namepartcont, const char *trigger, const char *nametrackcont, const char *nameclustercont, const char *suffix){
  return PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalJetEnergySpectrum::AddTaskJetEnergySpectrum(isMC, jettype, recoscheme, energydef, radius, namepartcont, trigger, nametrackcont, nameclustercont, suffix);
}
