EMCalTriggerPtAnalysis::AliAnalysisTaskChargedParticlesRefMC *AddTaskChargedParticlesRefMC(){

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  EMCalTriggerPtAnalysis::AliAnalysisTaskChargedParticlesRefMC *task = new EMCalTriggerPtAnalysis::AliAnalysisTaskChargedParticlesRefMC("chargedParticleMCQA");
  mgr->AddTask(task);

  TString outfile(mgr->GetCommonFileName());
  outfile += ":ChargedParticleQA";

  task->ConnectInput(0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, mgr->CreateContainer("Results", TList::Class(), AliAnalysisManager::kOutputContainer, outfile.Data()));

  return task;
}
