PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalTriggerEmulation *AddTaskEmcalTriggerEmulation(
    TString nametrackcont,
    TString nameclustercont,
    TString namemcparticlecont,
    const char *uff
    )
{
  PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalTriggerEmulation *trgtask = new PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalTriggerEmulation(Form("trgtask%s", uff));
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  mgr->AddTask(trgtask);

  TString fname = mgr->GetCommonFileName();
  fname += TString::Format("TriggerEmulation_%s", uff);

  // Create the containers
  if(namemcparticlecont.Length()){
    AliMCParticleContainer *mccont = trgtask->AddMCParticleContainer(namemcparticlecont.Data());
    trgtask->SetNameMCParticles(namemcparticlecont);
    mccont->SetCharge(AliParticleContainer::kCharged);
    mccont->SetEtaLimits(-0.5, 0.5);
    mccont->SetPhiLimits(1.41, 3.1);
  }

  if(nameclustercont.Length()){
    AliClusterContainer *clustercont = trgtask->AddClusterContainer(nameclustercont.Data());
    trgtask->SetNameClusters(nameclustercont);
  }

  if(nametrackcont.Length()){
    AliTrackContainer *trackcont = trgtask->AddTrackContainer(nametrackcont.Data());
    trgtask->SetNameTracks(nametrackcont.Data());
    trackcont->SetTrackFilterType(AliTrackContainer::k)
  }

  mgr->ConnectInput(trgtask, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(trgtask, 1, mgr->CreateContainer(Form("TriggerEmulation_%s", uff), TList::Class(), AliAnalysisManager::kOutputContainer, fname.Data()));

  return trgtask;
}
