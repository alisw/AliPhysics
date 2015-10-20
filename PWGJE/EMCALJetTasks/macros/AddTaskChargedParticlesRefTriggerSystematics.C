EMCalTriggerPtAnalysis::AliAnalysisTaskChargedParticlesRef *AddTaskChargedParticlesRefTriggerSystematics(
    double enEMC7,
    double enEG1,
    double enEG2,
    double enEJ1,
    double enEJ2
){

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  TString selectionstring = TString::Format("EMC7%dEG1%dEG2%dEJ1%dEJ2%d", int(enEMC7), int(enEG1), int(enEG2), int(enEJ1), int(enEJ2)),
      taskname = "ChargedParticleQA_" + selectionstring,
      listname = "TrackResults_" + selectionstring;

  EMCalTriggerPtAnalysis::AliAnalysisTaskChargedParticlesRef *task = new EMCalTriggerPtAnalysis::AliAnalysisTaskChargedParticlesRef(taskname.Data());
  mgr->AddTask(task);

  // Set Energy thresholds for additional patch selection:
  // These are events with offline patches of a given type where the trigger reached already the plateau
  // Tunable in order to study effects of the trigger systematics
  task->SetOfflineEnergyThreshold(EMCalTriggerPtAnalysis::AliAnalysisTaskChargedParticlesRef::kCPREL0, enEMC7);
  task->SetOfflineEnergyThreshold(EMCalTriggerPtAnalysis::AliAnalysisTaskChargedParticlesRef::kCPREG1, enEG1);
  task->SetOfflineEnergyThreshold(EMCalTriggerPtAnalysis::AliAnalysisTaskChargedParticlesRef::kCPREG2, enEG2);
  task->SetOfflineEnergyThreshold(EMCalTriggerPtAnalysis::AliAnalysisTaskChargedParticlesRef::kCPREJ1, enEJ1);
  task->SetOfflineEnergyThreshold(EMCalTriggerPtAnalysis::AliAnalysisTaskChargedParticlesRef::kCPREJ2, enEJ2);

  TString outfile(mgr->GetCommonFileName());
  outfile += ":" + taskname.Data();

  task->ConnectInput(0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, mgr->CreateContainer(listname.Data(), TList::Class(), AliAnalysisManager::kOutputContainer, outfile.Data()));

  return task;
}
