EMCalTriggerPtAnalysis::AliAnalysisTaskChargedParticlesMCTriggerMimic * AddTaskChargedParticlesMCTriggerMimic(const char *dummy = "", const char *suffix= ""){
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  EMCalTriggerPtAnalysis::AliAnalysisTaskChargedParticlesMCTriggerMimic * triggertask = new EMCalTriggerPtAnalysis::AliAnalysisTaskChargedParticlesMCTriggerMimic(Form("chargedParticlesMCtrigger_%s", suffix));
  mgr->AddTask(triggertask);
  triggertask->SetJetPtFactor(4.);
  triggertask->SetTrackPtFactor(1.5);
  triggertask->SetTrackSelection(
      EMCalTriggerPtAnalysis::AliEmcalAnalysisFactory::TrackCutsFactory(
          "standard",
          mgr->GetInputEventHandler()->IsA() == AliAODInputHandler::Class()
      )
  );


  TString outputfile = mgr->GetCommonFileName();
  outputfile += TString::Format(":ChargedParticleResults_%s", suffix);

  mgr->ConnectInput(triggertask, 1, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(triggertask, 1, mgr->CreateContainer(Form("TrackResults_%s", suffix), AliEmcalList::Class(), AliAnalysisManager::kOutputContainer, outputfile.Data()));

  return triggertask;
}
