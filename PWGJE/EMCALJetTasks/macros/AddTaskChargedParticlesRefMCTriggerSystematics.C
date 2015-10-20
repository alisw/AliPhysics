EMCalTriggerPtAnalysis::AliAnalysisTaskChargedParticlesRefMC *AddTaskChargedParticlesRefMCTriggerSystematics(
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

  EMCalTriggerPtAnalysis::AliAnalysisTaskChargedParticlesRefMC *task = new EMCalTriggerPtAnalysis::AliAnalysisTaskChargedParticlesRefMC(taskname.Data());
  task->SetOutlierCut(-1);
  // Set Energy thresholds for additional patch selection:
  // These are events with offline patches of a given type where the trigger reached already the plateau
  // These numers are determined as:
  // EMC7: 5 GeV
  // EG1:  14 GeV
  // EG2:  8 GeV
  // EJ1:  22 GeV
  // EJ2:  12 GeV
  task->SetOfflineEnergyThreshold(EMCalTriggerPtAnalysis::AliAnalysisTaskChargedParticlesRefMC::kCPREL0, enEMC7);
  task->SetOfflineEnergyThreshold(EMCalTriggerPtAnalysis::AliAnalysisTaskChargedParticlesRefMC::kCPREG1, enEG1);
  task->SetOfflineEnergyThreshold(EMCalTriggerPtAnalysis::AliAnalysisTaskChargedParticlesRefMC::kCPREG2, enEG2);
  task->SetOfflineEnergyThreshold(EMCalTriggerPtAnalysis::AliAnalysisTaskChargedParticlesRefMC::kCPREJ1, enEJ1);
  task->SetOfflineEnergyThreshold(EMCalTriggerPtAnalysis::AliAnalysisTaskChargedParticlesRefMC::kCPREJ2, enEJ2);
  mgr->AddTask(task);

  TString outfile(mgr->GetCommonFileName());
  outfile += ":" + taskname;

  task->ConnectInput(0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, mgr->CreateContainer(listname.Data(), TList::Class(), AliAnalysisManager::kOutputContainer, outfile.Data()));

  return task;
}
