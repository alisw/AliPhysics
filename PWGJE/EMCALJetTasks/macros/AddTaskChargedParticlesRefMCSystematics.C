EMCalTriggerPtAnalysis::AliAnalysisTaskChargedParticlesRefMC *AddTaskChargedParticlesRefMCSystematics(const char *dummy, const char *suffix){

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  TString taskname = "chargedParticleMCQA_" + TString(suffix);

  EMCalTriggerPtAnalysis::AliAnalysisTaskChargedParticlesRefMC *task = new EMCalTriggerPtAnalysis::AliAnalysisTaskChargedParticlesRefMC(taskname.Data());
  task->SetJetPtFactor(4.);
  task->SetTrackPtFactor(1.5);
  mgr->AddTask(task);

  TString outfile(mgr->GetCommonFileName());
  outfile += ":ChargedParticleQA_" + TString(suffix);

  task->ConnectInput(0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, mgr->CreateContainer(Form("TrackResults_%s", suffix), TList::Class(), AliAnalysisManager::kOutputContainer, outfile.Data()));

  return task;
}
