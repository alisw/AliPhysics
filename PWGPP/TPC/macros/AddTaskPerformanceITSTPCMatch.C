//____________________________________________
AliPerformanceTask* AddTaskPerformanceITSTPCMatch(Bool_t bUseMCInfo=kFALSE, Bool_t bUseESDfriend=kTRUE, const char *triggerClass=0)
{
  //
  // Add AliPerformanceTask with TPC performance components
  //
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if(!mgr) { 
    Error("AddTaskPerformanceTPC","AliAnalysisManager not set!");
    return NULL;
  }
  
  TString type = mgr->GetInputEventHandler()->GetDataType();
  if (!type.Contains("ESD")) {
    Error("AddTaskPerformanceTPC", "ESD input handler needed!");
    return NULL;
  }
  
  AliMCEventHandler *mcH = (AliMCEventHandler*)mgr->GetMCtruthEventHandler();
  if (!mcH && bUseMCInfo) {
    Error("AddTaskPerformanceTPC", "MC input handler needed!");
    return NULL;
  }

  //
  // Create task
  //
  AliPerformanceTask *task = new AliPerformanceTask("Performance","TPC Performance");
  if (!task) {
    Error("AddTaskPerformanceTPC", "TPC performance task cannot be created!");
    return NULL;
  }
  task->SetUseMCInfo(bUseMCInfo);
  task->SetUseESDfriend(bUseESDfriend);
  task->SelectCollisionCandidates();

  //
  // Add task to analysis manager
  //
  mgr->AddTask(task);


  //
  // Create TPC-ESD track reconstruction cuts
  // MATCH tracks
  //
  AliRecInfoCuts *pRecInfoCutsMATCH = new AliRecInfoCuts("pRecInfoCutsMATCH"); 
  if(pRecInfoCutsMATCH) {
    pRecInfoCutsMATCH->SetMaxDCAToVertexXY(3.0);
    pRecInfoCutsMATCH->SetMaxDCAToVertexZ(3..0);
    pRecInfoCutsMATCH->SetRequireSigmaToVertex(kFALSE);
    pRecInfoCutsMATCH->SetRequireTPCRefit(kFALSE);
    pRecInfoCutsMATCH->SetAcceptKinkDaughters(kTRUE);
    pRecInfoCutsMATCH->SetMinNClustersTPC(70);
    pRecInfoCutsMATCH->SetMaxChi2PerClusterTPC(1000000.);
    pRecInfoCutsMATCH->SetDCAToVertex2D(kFALSE);
    pRecInfoCutsMATCH->SetTPCITSMatchingRadius(70); 
    pRecInfoCutsMATCH->SetTPCTRDMatchingRadius(260); 
    pRecInfoCutsMATCH->SetMinNClustersITS(3);

    pRecInfoCutsMATCH->SetHistogramsOn(kFALSE); 
  } 
  else {
    Error("AddTaskPerformanceTPC", "AliRecInfoCutsTPC cannot be created!");
    return NULL;
  }

  //
  // Create TPC-MC track reconstruction cuts
  //
  AliMCInfoCuts  *pMCInfoCuts = new AliMCInfoCuts("pMCInfoCuts");
  if(pMCInfoCuts) {
    pMCInfoCuts->SetMinTrackLength(70);
  } 
  else {
    Error("AddTaskPerformanceTPC", "AliMCInfoCuts cannot be created!");
    return NULL;
  }

  //
  // Create performance objects for TPC and set cuts 
  //
  // TPC+ITS matching performance
  //
  AliPerformanceMatch *pCompMatch0 = new AliPerformanceMatch("AliPerformanceMatchTPCITS","AliPerformanceMatchTPCITS",0,kFALSE); 
  if(!pCompMatch0) {
    Error("AddTaskPerformanceMatch", "Cannot create AliPerformanceMatchTPCITS");
  }
  pCompMatch0->SetAliRecInfoCuts(pRecInfoCutsMATCH);
  pCompMatch0->SetAliMCInfoCuts(pMCInfoCuts);
  //
  // TPC+TRD matching performance
  //
  AliPerformanceMatch *pCompMatch1 = new AliPerformanceMatch("AliPerformanceMatchTPCTRD","AliPerformanceMatchTPCTRD",1,kFALSE); 
  if(!pCompMatch1) {
    Error("AddTaskPerformanceMatch", "Cannot create AliPerformanceMatchTPCTRD");
  }
  pCompMatch1->SetAliRecInfoCuts(pRecInfoCutsMATCH);
  pCompMatch1->SetAliMCInfoCuts(pMCInfoCuts);
 
  AliPerformanceMatch *pCompMatch2 = new AliPerformanceMatch("AliPerformanceMatchTPCEFF","AliPerformanceMatchTPCEFF",2,kFALSE); 
  if(!pCompMatch2) {
    Error("AddTaskPerformanceMatch", "Cannot create AliPerformanceMatchTPCEFF");
  }
  pCompMatch2->SetAliRecInfoCuts(pRecInfoCutsMATCH);
  pCompMatch2->SetAliMCInfoCuts(pMCInfoCuts);

  //
  // Add components to the performance task
  //
  task->AddPerformanceObject( pCompMatch0 );
  //task->AddPerformanceObject( pCompMatch1 );
  task->AddPerformanceObject( pCompMatch2 );

  //
  if(!bUseMCInfo &&  triggerClass) {
    pCompDEdx3->SetTriggerClass(triggerClass);
    pCompDCA0->SetTriggerClass(triggerClass);
    pCompTPC0->SetTriggerClass(triggerClass);
    pCompMatch0->SetTriggerClass(triggerClass);
    pCompMatch1->SetTriggerClass(triggerClass);
    pCompMatch2->SetTriggerClass(triggerClass);
  }

  //
  // Create containers for input
  //
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());

  //
  // Create containers for output
  //
  AliAnalysisDataContainer *coutput_tpc = mgr->CreateContainer("ITSTPCMatch", TList::Class(), AliAnalysisManager::kOutputContainer, Form("TPC.%sDet.root", task->GetName()));
  mgr->ConnectOutput(task, 1, coutput_tpc);

  AliAnalysisDataContainer *coutput2_tpc = mgr->CreateContainer("TPCTPCMatchSummary", TTree::Class(), AliAnalysisManager::kParamContainer, "trending.root:SummaryTPCQA"); 
  mgr->ConnectOutput(task, 2, coutput2_tpc);

return task;  
}
