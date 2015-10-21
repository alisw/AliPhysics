EMCalTriggerPtAnalysis::AliAnalysisTaskChargedParticlesRefMC *AddTaskChargedParticlesRefMCSystematics(const char *suffix){

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  TString taskname = "chargedParticleMCQA_" + TString(suffix);

  EMCalTriggerPtAnalysis::AliAnalysisTaskChargedParticlesRefMC *task = new EMCalTriggerPtAnalysis::AliAnalysisTaskChargedParticlesRefMC(taskname.Data());
  task->SetOutlierCut(-1);

  mgr->AddTask(task);

  TString outfile(mgr->GetCommonFileName());
  outfile += ":ChargedParticleQA_" + TString(suffix);

  task->ConnectInput(0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, mgr->CreateContainer(Form("TrackResults_outlier%d", int(outliercut * 10)), TList::Class(), AliAnalysisManager::kOutputContainer, outfile.Data()));

  return task;
}
