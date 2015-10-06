EMCalTriggerPtAnalysis::AliAnalysisTaskChargedParticlesRef *AddTaskChargedParticlesRef(){

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

  TString outfile(mgr->GetCommonFileName());
  outfile += ":ChargedParticleQA";

  task->ConnectInput(0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, mgr->CreateContainer("TrackResults", TList::Class(), AliAnalysisManager::kOutputContainer, outfile.Data()));

  return task;
}
