EMCalTriggerPtAnalysis::AliAnalysisTaskEmcalOnlinePatchesRef *AddTaskEmcalOnlinePatchesRef(){
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  EMCalTriggerPtAnalysis::AliAnalysisTaskEmcalOnlinePatchesRef *task = new EMCalTriggerPtAnalysis::AliAnalysisTaskEmcalOnlinePatchesRef("emcalOnlinePatchQA");
  mgr->AddTask(task);

  TString outfile(mgr->GetCommonFileName());
  outfile += ":OnlinePatchQA";

  task->ConnectInput(0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, mgr->CreateContainer("OnlinePatchResults", TList::Class(), AliAnalysisManager::kOutputContainer, outfile.Data()));

  return task;
}
