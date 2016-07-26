EMCalTriggerPtAnalysis::AliAnalysisTaskEventSelectionRef *AddTaskEventSelectionRef(
    const char *clustercontainername,
    Double_t enEMC7 = 5,
    Double_t enEG1 = 14,
    Double_t enEG2 = 8,
    Double_t enEJ1 = 22,
    Double_t enEJ2 = 12
    )
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  TString selectionstring = TString::Format("EMC7%dEG1%dEG2%dEJ1%dEJ2%d", int(enEMC7), int(enEG1), int(enEG2), int(enEJ1), int(enEJ2)),
      taskname = "EventSelectionQA_" + selectionstring,
      listname = "EventSelectionResults_" + selectionstring;

  EMCalTriggerPtAnalysis::AliAnalysisTaskEventSelectionRef *task = new EMCalTriggerPtAnalysis::AliAnalysisTaskEventSelectionRef(taskname.Data());
  mgr->AddTask(task);

  // Set Energy thresholds for additional patch selection
  // Threholds are configurable to test the influence of stronger event selections
  task->SetOfflineTriggerSelection(
      EMCalTriggerPtAnalysis::AliEmcalAnalysisFactory::TriggerSelectionFactory(enEMC7, enEG1, enEG2, enEJ1, enEJ2)
  );
  task->SetClusterContainer(clustercontainername);

  TString outfile(mgr->GetCommonFileName());
  outfile += ":" + taskname;

  task->ConnectInput(0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, mgr->CreateContainer(listname.Data(), TList::Class(), AliAnalysisManager::kOutputContainer, outfile.Data()));

  return task;

}
