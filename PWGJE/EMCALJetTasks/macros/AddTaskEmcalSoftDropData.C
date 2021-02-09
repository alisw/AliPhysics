PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalSoftDropData * AddTaskEmcalSoftDropData(Double_t jetradius, AliJetContainer::EJetType_t jettype, AliJetContainer::ERecoScheme_t recombinationScheme, AliVCluster::VCluUserDefEnergy_t energydef, const char *trigger) {
    return PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalSoftDropData::AddTaskEmcalSoftDropData(jetradius, jettype, recombinationScheme, energydef, trigger);
}
