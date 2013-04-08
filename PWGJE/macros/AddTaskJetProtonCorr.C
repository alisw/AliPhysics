AliAnalysisTaskJetProtonCorr* AddTaskJetProtonCorr(const char *name = "jet_prot_corr_01", const char *jetBranchName = "")
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  if(!mgr){
    ::Error("AddTaskJetProtonCorr", "No analysis manager to connect to.");
    return 0x0;
  }
  if(!mgr->GetInputEventHandler()){
    ::Error("AddTaskJetProtonCorr", "This task requires an input event handler.");
    return 0x0;
  }

  AliAnalysisTaskJetProtonCorr *task = new AliAnalysisTaskJetProtonCorr(name);

  if (strlen(jetBranchName) > 0) {
    task->SetJetBranchName(jetBranchName);
  }
  else {
    // branch existing in AOD115
    task->SetJetBranchName("clustersAOD_ANTIKT04_B1_Filter00768_Cut00150_Skip00");
  }

  AliAnalysisDataContainer *coutput =
    mgr->CreateContainer(Form("hist_%s", name), TList::Class(), AliAnalysisManager::kOutputContainer,
                         Form("%s:PWGJE_jet_prot_corr", AliAnalysisManager::GetCommonFileName()));

  if (!coutput) {
    ::Error("AddTaskJetProtonCorr", "no output container created");
    return 0x0;
  }

  mgr->AddTask(task);

  mgr->ConnectInput (task, 0, mgr->GetCommonInputContainer());
  // mgr->ConnectOutput(task, 0, mgr->GetCommonOutputContainer());
  mgr->ConnectOutput(task, 1, coutput);

  return task;
}
