EMCalTriggerPtAnalysis::AliAnalysisTaskTrackDensity *AddTaskTrackDensity(
    const char *mcjetcontainername,
    const char *suffix
) {
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  EMCalTriggerPtAnalysis::AliAnalysisTaskTrackDensity *densitytask = new EMCalTriggerPtAnalysis::AliAnalysisTaskTrackDensity("densitytask");
  mgr->AddTask(densitytask);
  mgr->ConnectInput(densitytask, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1,
      mgr->CreateContainer(Form("histos_%s", suffix),
          TList::Class(), AliAnalysisManager::kOutputContainer,
          Form("%s:TrackDensity_%s", mgr->GetCommonFileName(), suffix)));

  densitytask->SetIsPythia(true);

  AliMCParticleContainer *partcont = densitytask->AddMCParticleContainer("MCParticlesSelected");
  densitytask->SetMCParticleContainer("MCParticlesSelected");

  AliJetContainer *jetcont = densitytask->AddJetContainer(mcjetcontainername, AliJetContainer::kTPC, 0.4);
  jetcont->SetName("MCJetContainer");
  jetcont->SetJetPtCut(20.);
  jetcont->SetJetEtaLimits(-0.5, 0.5);
  densitytask->SetMCJetContainer("MCJetContainer");

  return densitytask;
}
