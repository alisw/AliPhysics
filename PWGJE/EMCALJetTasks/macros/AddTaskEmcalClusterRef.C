EMCalTriggerPtAnalysis::AliAnalysisTaskEmcalClustersRef *AddTaskEmcalClusterRef(TString clustercontname){

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  EMCalTriggerPtAnalysis::AliAnalysisTaskEmcalClustersRef *task = new EMCalTriggerPtAnalysis::AliAnalysisTaskEmcalClustersRef("emcalClusterQA");
  task->SetClusterContainer(clustercontname);
  mgr->AddTask(task);

  // Set Energy thresholds for additional patch selection:
  // These are events with offline patches of a given type where the trigger reached already the plateau
  // These numers are determined as:
  // EMC7: 3.5 GeV
  // EG1:  14 GeV
  // EG2:  8 GeV
  // EJ1:  22 GeV
  // EJ2:  12 GeV
  task->SetOfflineTriggerSelection(
      EMCalTriggerPtAnalysis::AliEmcalAnalysisFactory::TriggerSelectionFactory(5, 14, 8, 22, 12)
  );

  TString outfile(mgr->GetCommonFileName());
  outfile += ":ClusterQA";

  task->ConnectInput(0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, mgr->CreateContainer("ClusterResults", TList::Class(), AliAnalysisManager::kOutputContainer, outfile.Data()));

  return task;
}
