PWGJE::EMCALJetTasks::AliAnalysisTaskChargedParticlesMCTriggerMimic * AddTaskChargedParticlesMCTriggerMimic(const char *dummy = "", const char *suffix= ""){
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  PWGJE::EMCALJetTasks::AliAnalysisTaskChargedParticlesMCTriggerMimic * triggertask = new PWGJE::EMCALJetTasks::AliAnalysisTaskChargedParticlesMCTriggerMimic(Form("chargedParticlesMCtrigger_%s", suffix));
  mgr->AddTask(triggertask);
  triggertask->SetJetPtFactor(4.);
  triggertask->SetTrackPtFactor(1.5);
  triggertask->SetEmcalTrackSelection(
      PWGJE::EMCALJetTasks::AliEmcalAnalysisFactory::TrackCutsFactory(
          "standard",
          mgr->GetInputEventHandler()->IsA() == AliAODInputHandler::Class()
      )
  );


  TString outputfile = mgr->GetCommonFileName();
  outputfile += TString::Format(":ChargedParticleResults_%s", suffix);

  mgr->ConnectInput(triggertask, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(triggertask, 1, mgr->CreateContainer(Form("TrackResults_%s", suffix), AliEmcalList::Class(), AliAnalysisManager::kOutputContainer, outputfile.Data()));

  return triggertask;
}
