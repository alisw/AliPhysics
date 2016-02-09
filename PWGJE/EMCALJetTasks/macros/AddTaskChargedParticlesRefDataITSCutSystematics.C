EMCalTriggerPtAnalysis::AliAnalysisTaskChargedParticlesRefMC *AddTaskChargedParticlesRefMCITSCutSystematics(int itscut = 0, const char *suffix = ""){

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  TString cutname;

  EMCalTriggerPtAnalysis::AliAnalysisTaskChargedParticlesRefMC *task = new EMCalTriggerPtAnalysis::AliAnalysisTaskChargedParticlesRefMC("chargedParticleMCQA");
  // Set Energy thresholds for additional patch selection:
  // These are events with offline patches of a given type where the trigger reached already the plateau
  // These numers are determined as:
  // EMC7: 5 GeV
  // EG1:  14 GeV
  // EG2:  8 GeV
  // EJ1:  22 GeV
  // EJ2:  12 GeV
  task->SetOfflineEnergyThreshold(EMCalTriggerPtAnalysis::AliAnalysisTaskChargedParticlesRefMC::kCPREL0, 5);
  task->SetOfflineEnergyThreshold(EMCalTriggerPtAnalysis::AliAnalysisTaskChargedParticlesRefMC::kCPREG1, 14);
  task->SetOfflineEnergyThreshold(EMCalTriggerPtAnalysis::AliAnalysisTaskChargedParticlesRefMC::kCPREG2, 8);
  task->SetOfflineEnergyThreshold(EMCalTriggerPtAnalysis::AliAnalysisTaskChargedParticlesRefMC::kCPREJ1, 22);
  task->SetOfflineEnergyThreshold(EMCalTriggerPtAnalysis::AliAnalysisTaskChargedParticlesRefMC::kCPREJ2, 12);
  switch(itscut){
  case 0:
    cutname = "Standard";
    break;
  case 1:
    task->SwitchoffSPDCut();
    cutname = "NoSPDcut";
    break;
  case 2:
    task->SwitchoffITSCut();
    cutname = "NoITScut";
    break;
  }
  mgr->AddTask(task);

  TString outfile(mgr->GetCommonFileName());
  outfile += ":ChargedParticleQA_ITScut" + cutname;

  task->ConnectInput(0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, mgr->CreateContainer(Form("TrackResults_ITS%s", cutname.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, outfile.Data()));

  return task;
}
