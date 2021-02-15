PWGJE::EMCALJetTasks::AliAnalysisTaskEventSelectionRef *AddTaskEventSelectionRef(
    const char *clustercontainername,
    const char *dummy
    )
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  TString taskname = "EventSelectionQA_" + dummy,
      listname = "EventSelectionResults_" + dummy;

  PWGJE::EMCALJetTasks::AliAnalysisTaskEventSelectionRef *task = new PWGJE::EMCALJetTasks::AliAnalysisTaskEventSelectionRef(taskname.Data());
  mgr->AddTask(task);
  task->SetClusterContainer(clustercontainername);

  TString outfile(mgr->GetCommonFileName());
  outfile += ":" + taskname;

  task->ConnectInput(0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, mgr->CreateContainer(listname.Data(), TList::Class(), AliAnalysisManager::kOutputContainer, outfile.Data()));

  return task;

}
