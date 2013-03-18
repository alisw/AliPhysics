AliAnalysisTask *AddTaskJPSIFilter(Bool_t storeLS = kTRUE, Bool_t hasMC_aod = kFALSE){
  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskJPSIFilter", "No analysis manager found.");
    return 0;
  }
  
  //check for output aod handler
  if (!mgr->GetOutputEventHandler()||mgr->GetOutputEventHandler()->IsA()!=AliAODHandler::Class()) {
    Warning("AddTaskJPSIFilter","No AOD output handler available. Not adding the task!");
    return 0;
  }

  //Do we have an MC handler?
  Bool_t hasMC=(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()!=0x0)||hasMC_aod;
  
  //Do we run on AOD?
  Bool_t isAOD=mgr->GetInputEventHandler()->IsA()==AliAODInputHandler::Class();

  //gROOT->LoadMacro("$ALICE_ROOT/PWGDQ/dielectron/macros/ConfigBJpsi_ff_PbPbFilter.C");
  gROOT->LoadMacro("$ALICE_ROOT/PWGDQ/dielectron/macrosJPSI/ConfigBJpsi_ff_PbPbFilter.C");
  AliDielectron *jpsi=ConfigBJpsi_ff_PbPbFilter(0,hasMC);
  
  if(isAOD) {
    //add options to AliAODHandler to duplicate input event
    AliAODHandler *aodHandler = (AliAODHandler*)mgr->GetOutputEventHandler();
    aodHandler->SetCreateNonStandardAOD();
    aodHandler->SetNeedsHeaderReplication();
    aodHandler->SetNeedsTOFHeaderReplication();
    aodHandler->SetNeedsVZEROReplication();
    /*aodHandler->SetNeedsTracksBranchReplication();
    aodHandler->SetNeedsCaloClustersBranchReplication();
    aodHandler->SetNeedsVerticesBranchReplication();
    aodHandler->SetNeedsV0sBranchReplication();
    aodHandler->SetNeedsCascadesBranchReplication();
    aodHandler->SetNeedsTrackletsBranchReplication();
    aodHandler->SetNeedsPMDClustersBranchReplication();
    aodHandler->SetNeedsJetsBranchReplication();
    aodHandler->SetNeedsFMDClustersBranchReplication();
    //aodHandler->SetNeedsMCParticlesBranchReplication();
    aodHandler->SetNeedsDimuonsBranchReplication();*/
    if(hasMC) aodHandler->SetNeedsMCParticlesBranchReplication();
    jpsi->SetHasMC(hasMC);
  }
  
  //Create task and add it to the analysis manager
  AliAnalysisTaskDielectronFilter *task=new AliAnalysisTaskDielectronFilter("jpsi_DielectronFilter");
  task->SetTriggerMask(AliVEvent::kMB+AliVEvent::kCentral+AliVEvent::kSemiCentral);  
  if (!hasMC) task->UsePhysicsSelection();

  //Add event filter
  AliDielectronEventCuts *eventCuts=new AliDielectronEventCuts("eventCuts","Vertex Track && |vtxZ|<10 && ncontrib>0");
  if(!hasMC) eventCuts->SetRequireVertex();
  if (isAOD) eventCuts->SetVertexType(AliDielectronEventCuts::kVtxAny);
  eventCuts->SetMinVtxContributors(1);
  eventCuts->SetVertexZ(-10.,10.);
  eventCuts->SetCentralityRange(0.0,80.0);
  task->SetEventFilter(eventCuts);

  task->SetDielectron(jpsi);
  if(storeLS) task->SetStoreLikeSignCandidates(storeLS);
  task->SetCreateNanoAODs(kTRUE);
  task->SetStoreEventsWithSingleTracks(kTRUE);
  //task->SetStoreHeader(kTRUE);
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
