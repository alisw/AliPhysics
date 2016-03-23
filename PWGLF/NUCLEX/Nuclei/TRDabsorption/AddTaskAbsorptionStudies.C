AliAnalysisTaskAbsorptionStudies *AddTaskAbsorptionStudies(Bool_t usePhysicsSelection = true) {

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskAbsorptionStudies", "No analysis manager to connect to.");
    return NULL;
  }

  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("clist", TList::Class(), AliAnalysisManager::kOutputContainer, mgr->GetCommonFileName());
  AliAnalysisTaskAbsorptionStudies *simplePtTask = new AliAnalysisTaskAbsorptionStudies("AliAnalysisTaskAbsorptionStudies");
  mgr->AddTask(simplePtTask);

  mgr->ConnectInput(simplePtTask, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(simplePtTask, 1, coutput1);

  return simplePtTask;
}
