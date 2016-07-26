EMCalTriggerPtAnalysis::AliAnalysisTaskEmcalClustersRef *AddTaskEmcalClusterRefSystematics(TString clustercontname, const char *suffix){

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  TString taskname = "emcalClusterQA_" + TString(suffix);

  EMCalTriggerPtAnalysis::AliAnalysisTaskEmcalClustersRef *task = new EMCalTriggerPtAnalysis::AliAnalysisTaskEmcalClustersRef(taskname.Data());
  task->SetClusterContainer(clustercontname);
  mgr->AddTask(task);

  TString outfile(mgr->GetCommonFileName());
  outfile += ":ClusterQA_" + TString(suffix);
  TString containername = "ClusterResults_" + TString(suffix);
  printf("Outfile: %s, container: %s\n", outfile.Data(), containername.Data());

  task->ConnectInput(0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, mgr->CreateContainer(containername.Data(), TList::Class(), AliAnalysisManager::kOutputContainer, outfile.Data()));

  return task;
}
