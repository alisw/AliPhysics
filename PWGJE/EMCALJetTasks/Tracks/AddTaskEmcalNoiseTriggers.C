EMCalTriggerPtAnalysis::AliAnalysisTaskEmcalNoiseTriggers *AddTaskEmcalNoiseTriggers(const char *dummy, const char *suffix){
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  EMCalTriggerPtAnalysis::AliAnalysisTaskEmcalNoiseTriggers *task = new EMCalTriggerPtAnalysis::AliAnalysisTaskEmcalNoiseTriggers(Form("NoiseTriggerStudies%s", suffix));
  mgr->AddTask(task);

  TString outputfile(mgr->GetCommonFileName());
  outputfile += TString::Format(":NoiseTriggerStudies%s", suffix);

  mgr->ConnectInput(task, 1, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, mgr->CreateContainer(Form("NoiseTriggerHistograms%s", suffix), AliEmcalList::Class(), AliAnalysisManager::kOutputContainer, outputfile.Data()));

  return task;
}
