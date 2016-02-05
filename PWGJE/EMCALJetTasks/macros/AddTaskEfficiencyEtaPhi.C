EMCalTriggerPtAnalysis::AliAnalysisTaskEtaPhiEfficiency *AddTaskEfficiencyEtaPhi(){
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  EMCalTriggerPtAnalysis::AliAnalysisTaskEtaPhiEfficiency *efftask = new EMCalTriggerPtAnalysis::AliAnalysisTaskEtaPhiEfficiency("efficiencyTask");
  mgr->AddTask(efftask);

  TString outfile = mgr->GetCommonFileName();
  outfile += ":EfficiencyEtaPhi";

  mgr->ConnectInput(efftask, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(efftask, 1, mgr->CreateContainer("efficiencyHistos", TList::Class(), AliAnalysisManager::kOutputContainer, outfile.Data()));
  return efftask;
}
