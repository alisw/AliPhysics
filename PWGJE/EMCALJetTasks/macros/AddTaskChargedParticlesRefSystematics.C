EMCalTriggerPtAnalysis::AliAnalysisTaskChargedParticlesRef *AddTaskChargedParticlesRefSystematics(const char *dummy, const char *suffix){

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  TString taskname = "chargedParticleQA_" + TString(suffix);

  EMCalTriggerPtAnalysis::AliAnalysisTaskChargedParticlesRef *task = new EMCalTriggerPtAnalysis::AliAnalysisTaskChargedParticlesRef(taskname);
  mgr->AddTask(task);

  TString outfile(mgr->GetCommonFileName());
  outfile += ":ChargedParticleQA_" + TString(suffix);

  task->ConnectInput(0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, mgr->CreateContainer("TrackResults", TList::Class(), AliAnalysisManager::kOutputContainer, outfile.Data()));

  return task;
}
