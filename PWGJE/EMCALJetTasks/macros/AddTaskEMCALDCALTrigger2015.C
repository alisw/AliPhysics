EMCalTriggerPtAnalysis::AliAnalysisTaskEMCALDCALTrigger2015 AddTaskEMCALDCALTrigger2015(TString clustercontainer){

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  EMCalTriggerPtAnalysis::AliAnalysisTaskEMCALDCALTrigger2015 *task = new EMCalTriggerPtAnalysis::AliAnalysisTaskEMCALDCALTrigger2015("EMCALDCALL0");
  mgr->AddTask(task);
  task->SetClusterContainerName(clustercontainer.Data());

  TString outfile(mgr->GetCommonFileName());
  outfile += ":EMCALDCALL0";

  task->ConnectInput(0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, mgr->CreateContainer("EMCALDCALL0Results", TList::Class(), AliAnalysisManager::kOutputContainer, outfile.Data()));

  return task;
}
