AliAnalysisTask* AddTaskTRDtriggerCheck()
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    cerr << "No Analysis manager available" << endl;
    return 0x0;
  }

  AliAnalysisTaskTRDtriggerCheck *task = new AliAnalysisTaskTRDtriggerCheck("TRDtriggerCheck");
  mgr->AddTask(task);

  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput = mgr->GetCommonOutputContainer();

  AliAnalysisDataContainer *ctrklqa =
    mgr->CreateContainer("TRDtriggerCheck", TList::Class(),
                         AliAnalysisManager::kOutputContainer,
			 "trg_check.root");

  mgr->ConnectInput(task, 0, cinput);

  if (coutput)
    mgr->ConnectOutput(task, 0, coutput);
  mgr->ConnectOutput(task, 1, ctrklqa);

  return task;
}
