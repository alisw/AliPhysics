EMCalTriggerPtAnalysis::AliAnalysisTaskEventSelectionRef *AddTaskEventSelectionRef(
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
  task->SetOfflineEnergyThreshold(EMCalTriggerPtAnalysis::AliAnalysisTaskEventSelectionRef::kCPREL0, enEMC7);
  task->SetOfflineEnergyThreshold(EMCalTriggerPtAnalysis::AliAnalysisTaskEventSelectionRef::kCPREG1, enEG1);
  task->SetOfflineEnergyThreshold(EMCalTriggerPtAnalysis::AliAnalysisTaskEventSelectionRef::kCPREG2, enEG2);
  task->SetOfflineEnergyThreshold(EMCalTriggerPtAnalysis::AliAnalysisTaskEventSelectionRef::kCPREJ1, enEJ1);
  task->SetOfflineEnergyThreshold(EMCalTriggerPtAnalysis::AliAnalysisTaskEventSelectionRef::kCPREJ2, enEJ2);

  TString outfile(mgr->GetCommonFileName());
  outfile += ":" + taskname;

  task->ConnectInput(0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, mgr->CreateContainer(listname.Data(), TList::Class(), AliAnalysisManager::kOutputContainer, outfile.Data()));

  return task;

}
