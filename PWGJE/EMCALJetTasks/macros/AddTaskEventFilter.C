void AddTaskEventFilter() {
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  EMCalTriggerPtAnalysis::AliAnalysisTaskEventFilter *task = new EMCalTriggerPtAnalysis::AliAnalysisTaskEventFilter("eventFilterTask");
  mgr->AddTask(task);

  TString outputcont = TString(mgr->GetCommonFileName()) + ":eventfilter";
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, mgr->CreateContainer("eventfilterhists", TList::Class(), AliAnalysisManager::kOutputContainer, outputcont.Data()));
}
