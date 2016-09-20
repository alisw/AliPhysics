EMCalTriggerPtAnalysis::AliAnalysisTaskEGAMonitor *AddTaskEGAMonitor() {
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  EMCalTriggerPtAnalysis::AliAnalysisTaskEGAMonitor *egamonitor = new EMCalTriggerPtAnalysis::AliAnalysisTaskEGAMonitor("EGAMonitor");
  mgr->AddTask(egamonitor);

  TString outputcont = mgr->GetCommonFileName();
  outputcont += ":EGAmonitor";

  mgr->ConnectInput(egamonitor, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(egamonitor, 1, mgr->CreateContainer("EGAmonitorHists", TList::Class(), AliAnalysisManager::kOutputContainer, outputcont.Data()));

  return egamonitor;
}
