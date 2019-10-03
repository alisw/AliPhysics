EMCalTriggerPtAnalysis::AliAnalysisTaskEmcalOfflinePatchesRef *AddTaskEmcalOfflinePatchesRef(){
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  EMCalTriggerPtAnalysis::AliAnalysisTaskEmcalOfflinePatchesRef *task = new EMCalTriggerPtAnalysis::AliAnalysisTaskEmcalOfflinePatchesRef("emcalOnlinePatchQA");
  mgr->AddTask(task);

  TString outfile(mgr->GetCommonFileName());
  outfile += ":OfflinePatchQA";

  task->ConnectInput(0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, mgr->CreateContainer("OfflinePatchResults", TList::Class(), AliAnalysisManager::kOutputContainer, outfile.Data()));

  return task;
}
