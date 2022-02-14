AliAnalysisTaskJpsiJetFilter* AddTaskJpsiJetFilter_pp(Bool_t storeLS = kFALSE, Bool_t storeTR = kFALSE){
  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AliAnalysisTaskJpsiJetFilter", "No analysis manager found.");
    exit(1);
  }
  
  //check for output aod handler
  if (!mgr->GetOutputEventHandler()||mgr->GetOutputEventHandler()->IsA()!=AliAODHandler::Class()) {
    Error("AliAnalysisTaskJpsiJetFilter", "No AOD output handler available. Not adding the task!");
    exit(1);
  }

  //Do we have an MC handler?
  Bool_t hasMC=(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()!=0x0);
  //Do we run on AOD?
  Bool_t isAOD=mgr->GetInputEventHandler()->IsA()==AliAODInputHandler::Class();
  if(!isAOD || hasMC){
    Error("AliAnalysisTaskJpsiJetFilter", "This task ONLY supports AOD event");
    exit(1);
  }

  //add options to AliAODHandler to duplicate input event
  AliAODHandler *aodHandler = (AliAODHandler*)mgr->GetOutputEventHandler();
  aodHandler->SetCreateNonStandardAOD();
  aodHandler->SetNeedsHeaderReplication();
  aodHandler->SetNeedsTOFHeaderReplication();
  aodHandler->SetNeedsVZEROReplication();
  
  //Create task and add it to the analysis manager
  AliAnalysisTaskJpsiJetFilter *task=new AliAnalysisTaskJpsiJetFilter("jpsi2ee_EMCalFilter");
  task->SetTriggerMask(AliVEvent::kEMCEGA); 
  task->UsePhysicsSelection();
  task->SetRejectPileup(kTRUE);
  task->SetToMerge(kTRUE);
  task->SetStoreLikeSignCandidates(storeLS);
  task->SetStoreRotatedPairs(storeTR);
  //Add event filter
  AliDielectronEventCuts *eventCuts = new AliDielectronEventCuts("eventCuts", "Vertex Track && |vtxZ|<10 && ncontrib>0");
  eventCuts->SetVertexType(AliDielectronEventCuts::kVtxAny);
  eventCuts->SetRequireVertex();
  eventCuts->SetMinVtxContributors(1);
  eventCuts->SetVertexZ(-10., 10.);
  task->SetEventFilter(eventCuts);
  // Add AliDielectron - if not set, apply internal config
  task->InitDielectron();
  task->InitHistogramsForDielectron("DieFilterHistos");

  mgr->AddTask(task);

  //----------------------
  //create data containers
  //----------------------
  
  
  TString containerName = mgr->GetCommonFileName();
  containerName += ":PWGDQ_dielectronFilter";
  
  //create output container
  
  AliAnalysisDataContainer *cOutputHist1 =
    mgr->CreateContainer("jpsi_FilterQA",
                         THashList::Class(),
                         AliAnalysisManager::kOutputContainer,
                         containerName.Data());
  
  AliAnalysisDataContainer *cOutputHist2 =
    mgr->CreateContainer("jpsi_FilterEventStat",
                         TH1D::Class(),
                         AliAnalysisManager::kOutputContainer,
                         containerName.Data());
  
  
  mgr->ConnectInput(task,  0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, cOutputHist1);
  mgr->ConnectOutput(task, 2, cOutputHist2);
  
  return task;
}
