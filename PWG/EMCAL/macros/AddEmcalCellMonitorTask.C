AliEmcalCellMonitorTask *AddEmcalCellMonitorTask() {
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  AliEmcalCellMonitorTask *mon = new AliEmcalCellMonitorTask("cellMonitor");
  mgr->AddTask(mon);

  TString contname(mgr->GetCommonFileName());
  contname += ":cellMonitor";

  mgr->ConnectInput(mon, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(mon, 1, mgr->CreateContainer("cellMonitorHists", TList::Class(), AliAnalysisManager::kOutputContainer, contname.Data()));

  return mon;
}
