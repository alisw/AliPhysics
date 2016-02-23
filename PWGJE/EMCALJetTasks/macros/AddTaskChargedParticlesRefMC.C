EMCalTriggerPtAnalysis::AliAnalysisTaskChargedParticlesRefMC *AddTaskChargedParticlesRefMC(const char *cutname = "standard"){

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  EMCalTriggerPtAnalysis::AliAnalysisTaskChargedParticlesRefMC *task = new EMCalTriggerPtAnalysis::AliAnalysisTaskChargedParticlesRefMC(Form("chargedParticleMCQA_%s", suffix));
  task->SetOutlierCut(-1);
  // Set Energy thresholds for additional patch selection:
  // These are events with offline patches of a given type where the trigger reached already the plateau
  // These numers are determined as:
  // EMC7: 5 GeV
  // EG1:  14 GeV
  // EG2:  8 GeV
  // EJ1:  22 GeV
  // EJ2:  12 GeV
  mgr->AddTask(task);
  task->SetOfflineTriggerSelection(
      EMCalTriggerPtAnalysis::AliEmcalAnalysisFactory::TriggerSelectionFactory(5, 14, 8, 22, 12)
  );
  task->SetTrackSelection(
      EMCalTriggerPtAnalysis::AliEmcalAnalysisFactory::TrackCutsFactory(
          cutname,
          mgr->GetInputEventHandler()->IsA() == AliAODInputHandler::Class()
      )
  );

  TString outfile(mgr->GetCommonFileName());
  outfile += ":ChargedParticleQA";

  task->ConnectInput(0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, mgr->CreateContainer("TrackResults", TList::Class(), AliAnalysisManager::kOutputContainer, outfile.Data()));

  return task;
}
