PWGJE::EMCALJetTasks::AliAnalysisTaskTrackDensity *AddTaskTrackDensity(
    const char *mcparticlecontainername,
    const char *mcjetcontainername,
    const char *suffix
) {
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  Bool_t isAOD = mgr->GetInputEventHandler()->IsA() == AliAODInputHandler::Class();

  PWGJE::EMCALJetTasks::AliAnalysisTaskTrackDensity *densitytask = new PWGJE::EMCALJetTasks::AliAnalysisTaskTrackDensity("densitytask");
  mgr->AddTask(densitytask);
  mgr->ConnectInput(densitytask, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(densitytask, 1,
      mgr->CreateContainer(Form("histos_%s", suffix),
          TList::Class(), AliAnalysisManager::kOutputContainer,
          Form("%s:TrackDensity_%s", mgr->GetCommonFileName(), suffix)));

  densitytask->SetOffTrigger(AliVEvent::kINT7);
  densitytask->SetIsPythia(true);
  densitytask->SetUseAliAnaUtils(true, true);
  densitytask->SetTrackSelection(PWGJE::EMCALJetTasks::AliEmcalAnalysisFactory::TrackCutsFactory("standard", isAOD));

  AliMCParticleContainer *partcont = densitytask->AddMCParticleContainer(mcparticlecontainername);
  partcont->SetName("MCParticles");
  densitytask->SetMCParticleContainer("MCParticles");

  AliJetContainer *jetcont = densitytask->AddJetContainer(mcjetcontainername, AliEmcalJet::kTPC, 0.4);
  jetcont->SetName("MCJetContainer");
  jetcont->SetJetPtCut(20.);
  jetcont->SetJetEtaLimits(-0.5, 0.5);
  jetcont->ConnectParticleContainer(partcont);
  densitytask->SetMCJetContainer("MCJetContainer");

  return densitytask;
}
