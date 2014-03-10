AliAnalysisTask* AddTaskTRDtriggerCheck(const char *name = "trd_trgcheck")
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    cerr << "No Analysis manager available" << endl;
    return 0x0;
  }

  AliAnalysisTaskTRDtriggerCheck *task = new AliAnalysisTaskTRDtriggerCheck("TRDtriggerCheck");
  // task->SetDebugLevel(2);
  mgr->AddTask(task);

  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput = mgr->GetCommonOutputContainer();

  AliAnalysisDataContainer *chist =
    mgr->CreateContainer(Form("hist_%s", name), TList::Class(), AliAnalysisManager::kOutputContainer,
                         Form("%s:PWGTRD_trgcheck", AliAnalysisManager::GetCommonFileName()));

  if (!chist) {
    ::Error("AddTaskTRDtriggerCheck", "no output container created");
    return 0x0;
  }

  mgr->ConnectInput(task, 0, cinput);

  if (coutput)
    mgr->ConnectOutput(task, 0, coutput);
  mgr->ConnectOutput(task, 1, chist);

  return task;
}
