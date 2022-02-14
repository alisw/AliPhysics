AliAnalysisTaskSEMultDQdev *AddTaskMultDQdev(const Bool_t bOutputList = kTRUE)
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  if (!mgr) {
    ::Error("AddTaskMuonDistributions", "No analysis manager to connect to.");
     return nullptr;
  }
//=============================================================================

  auto task(new AliAnalysisTaskSEMultDQdev("AliAnalysisTaskSEMultDQdev"));
  mgr->AddTask(task);
  mgr->ConnectInput( task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 0, mgr->GetCommonOutputContainer());
  if (bOutputList) mgr->ConnectOutput(task, 1, mgr->CreateContainer("listEvts", TList::Class(),
                                                                    AliAnalysisManager::kOutputContainer,
                                                                    AliAnalysisManager::GetCommonFileName()));
//=============================================================================

  return task;
}
