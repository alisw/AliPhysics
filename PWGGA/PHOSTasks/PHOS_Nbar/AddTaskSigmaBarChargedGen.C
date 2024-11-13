AliAnalysisSigmaBarChargedGen *AddTaskSigmaBarChargedGen() {

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AliAnalysisSigmaBarChargedGen",
            "No analysis manager to connect to.");
    return NULL;
  }

  if (!mgr->GetInputEventHandler())
    printf("Caution: No input handler detected!");

  AliAnalysisSigmaBarChargedGen *taskMCPred =
      new AliAnalysisSigmaBarChargedGen("taskMCPred");

  mgr->AddTask(taskMCPred);
  TString outputFileName = AliAnalysisManager::GetCommonFileName();
  outputFileName += ":Sigma_MCpred";
  Printf("Set OutputFileName : \n %s\n", outputFileName.Data());

  AliAnalysisDataContainer *coutputList = mgr->CreateContainer(
      "fHistos_misc", TList::Class(), AliAnalysisManager::kOutputContainer,
      outputFileName);

  mgr->ConnectInput(taskMCPred, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(taskMCPred, 1, coutputList);

  return taskMCPred;
}
