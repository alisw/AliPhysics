EMCalTriggerPtAnalysis::AliAnalysisTaskEmcalPatchesRef *AddTaskEmcalPatchesRefSystematics(const char *suffix){

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  TString taskname = "emcalPatchQA_" + TString(suffix);

  EMCalTriggerPtAnalysis::AliAnalysisTaskEmcalPatchesRef *task = new EMCalTriggerPtAnalysis::AliAnalysisTaskEmcalPatchesRef(taskname.Data());

  mgr->AddTask(task);

  TString outfile(mgr->GetCommonFileName());
  outfile += ":PatchQA_" + TString(suffix);

  task->ConnectInput(0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, mgr->CreateContainer("PatchResults", TList::Class(), AliAnalysisManager::kOutputContainer, outfile.Data()));

  return task;
}
