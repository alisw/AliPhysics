EMCalTriggerPtAnalysis::AliAnalysisTaskEmcalClustersRef *AddTaskEmcalClusterRefSystematics(const char * nClusters, const char *suffix){

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  enum EDataType_t {
    kUnknown,
    kESD,
    kAOD
  };

  AliVEventHandler* handler = mgr->GetInputEventHandler();
  EDataType_t dataType = kUnknown;

  if (handler->InheritsFrom("AliESDInputHandler")) {
    dataType = kESD;
  }
  else if (handler->InheritsFrom("AliAODInputHandler")) {
    dataType = kAOD;
  }

  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------

  TString clusName(nClusters);

  if (clusName == "usedefault") {
    if (dataType == kESD) {
      clusName = "CaloClusters";
    }
    else if (dataType == kAOD) {
      clusName = "caloClusters";
    }
    else {
      clusName = "";
    }
  }

  TString taskname = "emcalClusterQA_" + TString(suffix);

  EMCalTriggerPtAnalysis::AliAnalysisTaskEmcalClustersRef *task = new EMCalTriggerPtAnalysis::AliAnalysisTaskEmcalClustersRef(taskname.Data());
  task->AddClusterContainer(clusName);
  task->SetClusterContainer(clusName);
  mgr->AddTask(task);

  TString outfile(mgr->GetCommonFileName());
  outfile += ":ClusterQA_" + TString(suffix);
  TString containername = "ClusterResults_" + TString(suffix);
  printf("Outfile: %s, container: %s\n", outfile.Data(), containername.Data());

  task->ConnectInput(0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, mgr->CreateContainer(containername.Data(), TList::Class(), AliAnalysisManager::kOutputContainer, outfile.Data()));

  return task;
}
