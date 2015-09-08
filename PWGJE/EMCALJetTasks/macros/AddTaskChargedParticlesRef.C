EMCalTriggerPtAnalysis::AliAnalysisTaskChargedParticlesRef *AddTaskChargedParticlesRef(){

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  EMCalTriggerPtAnalysis::AliAnalysisTaskChargedParticlesRef *task = new EMCalTriggerPtAnalysis::AliAnalysisTaskChargedParticlesRef("chargedParticleQA");
  mgr->AddTask(task);

  TString outfile(mgr->GetCommonFileName());
  outfile += ":ChargedParticleQA";

  task->ConnectInput(0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, mgr->CreateContainer("TrackResults", TList::Class(), AliAnalysisManager::kOutputContainer, outfile.Data()));

  return task;
}
