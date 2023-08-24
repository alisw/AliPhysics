PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalJetEnergySpectrum *AddTaskEmcalJetEnergySpectrumBkgSub(Bool_t isMC, AliJetContainer::EJetType_t jettype, AliJetContainer::ERecoScheme_t recoscheme, AliVCluster::VCluUserDefEnergy_t energydef, Double_t radius, const char *namepartcont, const char *trigger, const char *nrho, const char *suffix){
  return PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalJetEnergySpectrum::AddTaskJetEnergySpectrumBkgSub(isMC, jettype, recoscheme, energydef, radius, namepartcont, trigger, nrho, suffix);
}

PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalJetEnergySpectrum *AddTaskEmcalJetEnergySpectrumBkgSub(Bool_t isMC, AliJetContainer::EJetType_t jettype, AliJetContainer::ERecoScheme_t recoscheme, AliVCluster::VCluUserDefEnergy_t energydef, Double_t radius, const char *namepartcont, const char *trigger, const char *nrho, const char *nametrackcont, const char *nameclustercont, const char *suffix){
  return PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalJetEnergySpectrum::AddTaskJetEnergySpectrumBkgSub(isMC, jettype, recoscheme, energydef, radius, namepartcont, trigger, nrho, nametrackcont, nameclustercont, suffix);
}
