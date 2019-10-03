AliAnalysisTaskJetsTriggerTRD* AddTaskJetsTriggerTRD(const char *name = "jets_trg_trd", const char *jetBranchName = "")
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if(!mgr){
    ::Error("AddTaskJetsTriggerTRD", "No analysis manager to connect to.");
    return NULL;
  }
  if(!mgr->GetInputEventHandler()){
    ::Error("AddTaskJetsTriggerTRD", "This task requires an input event handler.");
    return NULL;
  }

  Bool_t isMC = (mgr->GetMCtruthEventHandler() != 0x0);
  TString inputDataType = mgr->GetInputEventHandler()->GetDataType();

  AliAnalysisTaskJetsTriggerTRD *task = new AliAnalysisTaskJetsTriggerTRD(name);

  // if no jet branch is specified, set default values depending on data type
  if (strlen(jetBranchName) > 0) {
    task->SetJetBranchName(jetBranchName);
  }
  else {
    // for ESDs (MC and real)  we use a preceding jet finder
    if (isMC) {
      task->SetJetBranchName("clustersAODMC_ANTIKT04_B0_Filter00272_Cut00150_Skip00");
    }
    else if (inputDataType.Contains("ESD")) {
      task->SetJetBranchName("clustersAOD_ANTIKT04_B0_Filter00272_Cut00150_Skip00");
    }
    // for AODs we use an existing jet branch
    else if (inputDataType.Contains("AOD")) {
      task->SetJetBranchName("clustersAOD_ANTIKT02_B0_Filter00272_Cut00150_Skip00");
    }
    else {
      printf("unknown input data type\n");
    }
  }

  mgr->AddTask(task);

  AliAnalysisDataContainer *coutput =
    mgr->CreateContainer(Form("hist_%s", name), TList::Class(), AliAnalysisManager::kOutputContainer,
                         Form("%s:PWGJE_jets_trg_trd", AliAnalysisManager::GetCommonFileName()));

  mgr->ConnectInput (task, 0, mgr->GetCommonInputContainer());
  if (mgr->GetCommonOutputContainer())
    mgr->ConnectOutput(task, 0, mgr->GetCommonOutputContainer());
  mgr->ConnectOutput(task, 1, coutput);

  return task;
}
