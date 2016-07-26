AliAnalysisTaskAbsorptionStudies *AddTaskAbsorptionStudies(Bool_t usePhysicsSelection = true, TString suffix = "") {
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskAbsorptionStudies", "No analysis manager to connect to.");
    return NULL;
  }

  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(Form("absorption%s",suffix.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, mgr->GetCommonFileName());
  AliAnalysisTaskAbsorptionStudies *simplePtTask = new AliAnalysisTaskAbsorptionStudies(Form("AbsorptionStudies%s",suffix.Data()));
  mgr->AddTask(simplePtTask);

  mgr->ConnectInput(simplePtTask, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(simplePtTask, 1, coutput1);

  return simplePtTask;
}
