EMCalTriggerPtAnalysis::AliAnalysisTaskChargedParticlesRef *AddTaskChargedParticlesRef(const char *cutname = "standard"){

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  EMCalTriggerPtAnalysis::AliAnalysisTaskChargedParticlesRef *task = new EMCalTriggerPtAnalysis::AliAnalysisTaskChargedParticlesRef("chargedParticleQA");
  mgr->AddTask(task);

  // Set Energy thresholds for additional patch selection:
  // These are events with offline patches of a given type where the trigger reached already the plateau
  // These numers are determined as:
  // EMC7: 5 GeV
  // EG1:  14 GeV
  // EG2:  8 GeV
  // EJ1:  22 GeV
  // EJ2:  12 GeV
  task->SetOfflineEnergyThreshold(EMCalTriggerPtAnalysis::AliAnalysisTaskChargedParticlesRef::kCPREL0, 5);
  task->SetOfflineEnergyThreshold(EMCalTriggerPtAnalysis::AliAnalysisTaskChargedParticlesRef::kCPREG1, 14);
  task->SetOfflineEnergyThreshold(EMCalTriggerPtAnalysis::AliAnalysisTaskChargedParticlesRef::kCPREG2, 8);
  task->SetOfflineEnergyThreshold(EMCalTriggerPtAnalysis::AliAnalysisTaskChargedParticlesRef::kCPREJ1, 22);
  task->SetOfflineEnergyThreshold(EMCalTriggerPtAnalysis::AliAnalysisTaskChargedParticlesRef::kCPREJ2, 12);
  mgr->AddTask(task);
  task->SetTrackSelection(
      EMCalTriggerPtAnalysis::AliAnalysisTaskChargedParticlesRef::TrackCutsFactory(
          cutname,
          mgr->GetInputEventHandler()->IsA() == AliAODInputHandler::Class()
      )
  );

  TString outfile(mgr->GetCommonFileName());
  outfile += TString::Format(":ChargedParticleQA_%s", cutname);

  task->ConnectInput(0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, mgr->CreateContainer(Form("TrackResults_%s", cutname), TList::Class(), AliAnalysisManager::kOutputContainer, outfile.Data()));

  return task;
}
