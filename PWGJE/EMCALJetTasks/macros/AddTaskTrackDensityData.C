PWGJE::EMCALJetTasks::AliAnalysisTaskTrackDensityData *AddTaskTrackDensityData(
    const char *jetcontainername,
    const char *trackcontainername,
    const char *period,
    const char *suffix
) {
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  Bool_t isAOD = mgr->GetInputEventHandler()->IsA() == AliAODInputHandler::Class();

  PWGJE::EMCALJetTasks::AliAnalysisTaskTrackDensityData *densitytask = new PWGJE::EMCALJetTasks::AliAnalysisTaskTrackDensityData("densitytaskdata");
  mgr->AddTask(densitytask);
  mgr->ConnectInput(densitytask, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(densitytask, 1,
      mgr->CreateContainer(Form("histos_%s", suffix),
          TList::Class(), AliAnalysisManager::kOutputContainer,
          Form("%s:TrackDensityData_%s", mgr->GetCommonFileName(), suffix)));

  densitytask->SetOffTrigger(AliVEvent::kINT7);
  densitytask->SetUseAliAnaUtils(true, true);
  densitytask->SetEmcalTrackSelection(PWGJE::EMCALJetTasks::AliEmcalAnalysisFactory::TrackCutsFactory("standard", isAOD));

  AliTrackContainer *trackcont = densitytask->AddTrackContainer(trackcontainername);
  trackcont->SetName("trackcontainer");
  trackcont->SetEtaLimits(-0.8, 0.8);
  trackcont->SetTrackCutsPeriod(period);
  trackcont->SetTrackFilterType(AliEmcalTrackSelection::kHybridTracks);
  densitytask->SetNameTrackContainer("trackcontainer");

  AliJetContainer *jetcont = densitytask->AddJetContainer(jetcontainername, AliEmcalJet::kTPC, 0.4);
  jetcont->SetName("jetcontainer");
  jetcont->SetJetPtCut(20.);
  jetcont->SetJetEtaLimits(-0.5, 0.5);
  jetcont->ConnectParticleContainer(trackcont);
  densitytask->SetNameJetContainer("jetcontainer");

  return densitytask;
}
