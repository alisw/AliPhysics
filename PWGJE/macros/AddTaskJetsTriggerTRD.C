AliAnalysisTaskJetsTriggerTRD* AddTaskJetsTriggerTRD(const char *name = "jets_trg_trd", const char *jetBranchName = "")
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if(!mgr){
    ::Error("AddTaskTriggerJets", "No analysis manager to connect to.");
    return NULL;
  }
  if(!mgr->GetInputEventHandler()){
    ::Error("AddTaskTriggerJets", "This task requires an input event handler.");
    return NULL;
  }

  AliAnalysisTaskJetsTriggerTRD *task = new AliAnalysisTaskJetsTriggerTRD(name);

  if (strlen(jetBranchName) > 0) {
    task->SetJetBranchName(jetBranchName);
  }
  else {
    TString inputDataType = mgr->GetInputEventHandler()->GetDataType();

    if (inputDataType.Contains("ESD")) {
      // for ESDs we use a preceding jet finder
      task->SetJetBranchName("clustersAOD_ANTIKT04_B0_Filter00272_Cut00150_Skip00");
    }
    else if (inputDataType.Contains("AOD")) {
      // for AODs we use an existing jet branch
      task->SetJetBranchName("jetsAOD_UA104_B0_Filter00768_Cut00150");
    }
    else {
      printf("unknown input data type\n");
    }
  }

  mgr->AddTask(task);

  AliAnalysisDataContainer *coutput =
    mgr->CreateContainer("histograms", TList::Class(), AliAnalysisManager::kOutputContainer,
                         Form("%s:PWGJE_trg_trd", AliAnalysisManager::GetCommonFileName()));

  mgr->ConnectInput (task, 0, mgr->GetCommonInputContainer());
  // mgr->ConnectOutput(task, 0, mgr->GetCommonOutputContainer());
  mgr->ConnectOutput(task, 1, coutput);

  return task;
}
