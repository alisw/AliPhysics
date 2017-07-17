PWG::EMCAL::AliEmcalFastOrMonitorTask *AddEmcalFastOrMonitor(const char *dummy = "", const char *subwagonname = "fastOrMonitor"){
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  PWG::EMCAL::AliEmcalFastOrMonitorTask *mon = new AliEmcalFastOrMonitorTask(subwagonname);
  mgr->AddTask(mon);

  TString contname(mgr->GetCommonFileName());
  contname += Form(":%s", subwagonname);

  mgr->ConnectInput(mon, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(mon, 1, mgr->CreateContainer(Form("%sHists", subwagonname), TList::Class(), AliAnalysisManager::kOutputContainer, contname.Data()));

  return mon;
}
