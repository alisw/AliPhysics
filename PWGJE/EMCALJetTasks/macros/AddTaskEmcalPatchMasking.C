PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalPatchMasking *AddTaskEmcalPatchMasking(const char *dummy, const char *suffix){
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalPatchMasking *maskingtask = new PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalPatchMasking(Form("MaskingTask%s", suffix));
  mgr->AddTask(maskingtask);

  TString outputcont = mgr->GetCommonFileName();
  outputcont += TString::Format(":PatchMasking%s", suffix);

  mgr->ConnectInput(maskingtask, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(maskingtask, 1, mgr->CreateContainer(Form("HistosPatchMasking%s", suffix), AliEmcalList::Class(), AliAnalysisManager::kOutputContainer, outputcont.Data()));

  return maskingtask;
}
