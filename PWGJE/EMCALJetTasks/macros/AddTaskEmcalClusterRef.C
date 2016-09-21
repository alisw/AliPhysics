EMCalTriggerPtAnalysis::AliAnalysisTaskEmcalClustersRef *AddTaskEmcalClusterRef(TString nClusters = "usedefault"){

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


  EMCalTriggerPtAnalysis::AliAnalysisTaskEmcalClustersRef *task = new EMCalTriggerPtAnalysis::AliAnalysisTaskEmcalClustersRef("emcalClusterQA");
  task->AddClusterContainer(clusName.Data());
  task->SetClusterContainer(clusName);
  mgr->AddTask(task);

  // Set Energy thresholds for additional patch selection:
  // These are events with offline patches of a given type where the trigger reached already the plateau
  // These numers are determined as:
  // EMC7: 3.5 GeV
  // EG1:  14 GeV
  // EG2:  8 GeV
  // EJ1:  22 GeV
  // EJ2:  12 GeV
  task->SetOfflineTriggerSelection(
      EMCalTriggerPtAnalysis::AliEmcalAnalysisFactory::TriggerSelectionFactory(5, 14, 8, 22, 12)
  );

  TString outfile(mgr->GetCommonFileName());
  outfile += ":ClusterQA";

  task->ConnectInput(0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, mgr->CreateContainer("ClusterResults", AliEmcalList::Class(), AliAnalysisManager::kOutputContainer, outfile.Data()));

  return task;
}
