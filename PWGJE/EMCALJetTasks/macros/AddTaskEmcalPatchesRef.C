EMCalTriggerPtAnalysis::AliAnalysisTaskEmcalPatchesRef *AddTaskEmcalPatchesRef(TString clustercontname){

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  EMCalTriggerPtAnalysis::AliAnalysisTaskEmcalPatchesRef *task = new EMCalTriggerPtAnalysis::AliAnalysisTaskEmcalPatchesRef("emcalPatchQA");
  mgr->AddTask(task);

  TString outfile(mgr->GetCommonFileName());
  outfile += ":ClusterQA";

  task->ConnectInput(0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, mgr->CreateContainer("PatchResults", TList::Class(), AliAnalysisManager::kOutputContainer, outfile.Data()));

  return task;
}
