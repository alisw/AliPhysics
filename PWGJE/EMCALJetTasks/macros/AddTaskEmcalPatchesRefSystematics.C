PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalPatchesRef *AddTaskEmcalPatchesRefSystematics(const char *dummy, const char *suffix){

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  TString taskname = "emcalPatchQA_" + TString(suffix);

  PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalPatchesRef *task = new PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalPatchesRef(taskname.Data());

  mgr->AddTask(task);

  TString outfile(mgr->GetCommonFileName());
  outfile += ":PatchQA_" + TString(suffix);
  TString containername = "PatchResults_" + TString(suffix);
  printf("Outfile: %s, container: %s\n", outfile.Data(), containername.Data());

  task->ConnectInput(0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, mgr->CreateContainer(containername.Data(), TList::Class(), AliAnalysisManager::kOutputContainer, outfile.Data()));

  return task;
}
