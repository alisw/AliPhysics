AliAnalysisTaskHIMultCorr *AddTaskHIMultCorr(){

#define USE_STREAMER 0

  //===============================================
  //-- Get the current analysis manager

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_jthaeder_HIMultCorr", "No analysis manager found.");
    return 0;
  }

  //===============================================
  // -- Check for MC

  Bool_t isMC=(mgr->GetMCtruthEventHandler() != NULL);
  if (isMC){
    Info("AddTask_jthaeder_HIMultCorr", "This task has MC.");
  }

  //===============================================
  // -- Create ESD cuts

  AliESDtrackCuts* esdCuts = new AliESDtrackCuts("AliESDtrackCuts");
  AliESDtrackCuts* esdCuts2 = new AliESDtrackCuts("AliESDtrackCuts2");

  // -- dNdPt track cuts 23
  esdCuts->SetRequireSigmaToVertex(kFALSE);
  esdCuts->SetRequireITSRefit(kFALSE);
  esdCuts->SetRequireTPCRefit(kFALSE);
  esdCuts->SetRequireTPCStandAlone(kTRUE);
  esdCuts->SetAcceptKinkDaughters(kFALSE);
  esdCuts->SetMinNClustersTPC(70);

  esdCuts->SetMaxChi2PerClusterTPC(4.0);

  esdCuts->SetMaxDCAToVertexXY(2.4); // cm
  esdCuts->SetMaxDCAToVertexZ(3.2);  // cm
  esdCuts->SetDCAToVertex2D(kTRUE);

  // -- dNdPt track cuts 72
  esdCuts2->SetRequireSigmaToVertex(kFALSE);
  esdCuts2->SetRequireTPCRefit(kTRUE);
  esdCuts2->SetRequireITSRefit(kTRUE);
  esdCuts2->SetAcceptKinkDaughters(kFALSE);
  esdCuts2->SetMinNClustersTPC(70);

  esdCuts2->SetMaxChi2PerClusterTPC(4.0);
  esdCuts2->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);

  esdCuts2->SetMaxDCAToVertexZ(2.0);
  esdCuts2->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01"); // 7*(0.0026+0.0050/pt^1.01)
  esdCuts2->SetDCAToVertex2D(kFALSE);

  // -- acceptance
  esdCuts->SetEtaRange(-0.8,0.8);
  esdCuts2->SetEtaRange(-0.8,0.8);

  //===============================================
  // -- Add task to the ANALYSIS manager

  AliAnalysisTaskHIMultCorr *task =  new AliAnalysisTaskHIMultCorr("AliAnalysisTaskHIMultCorrVZERO");
  task->SelectCollisionCandidates(AliVEvent::kMB);
  task->SetUseCentrality(1);  // V0
  task->SetESDCuts(esdCuts);
  task->SetESDCuts2(esdCuts2);
  task->SetMaxVertexZ(10.);    // cm
  if (isMC) task->SetIsMC();
  mgr->AddTask(task);
  
#ifdef USE_SPD
  AliAnalysisTaskHIMultCorr *taskSPD = new AliAnalysisTaskHIMultCorr("AliAnalysisTaskHIMultCorrSPD");
  taskSPD->SelectCollisionCandidates(AliVEvent::kMB);
  taskSPD->SetUseCentrality(2);  // SPD
  taskSPD->SetESDCuts(esdCuts);
  taskSPD->SetESDCuts2(esdCuts2);
  taskSPD->SetMaxVertexZ(10.);   // cm
  if (isMC) task->SetIsMC();
  mgr->AddTask(taskSPD);
#endif

  //================================================
  // -- Output for Streamer
#if USE_STREAMER  
  TString addoutput = gSystem->Getenv("ADD_OUTPUT_FILES");
  if (addoutput.Length()) addoutput+=",";
  addoutput+="eventInfoCoirr.root";
  gSystem->Setenv("ADD_OUTPUT_FILES",addoutput.Data());
#endif
  //================================================
  // -- data containers - input

  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();

  //================================================
  // -- data containers - VZERO
  
  TString outputFileName = AliAnalysisManager::GetCommonFileName();
  AliAnalysisDataContainer *coutput = 
    mgr->CreateContainer("HIMultCorr_VZERO", TList::Class(),
			 AliAnalysisManager::kOutputContainer,
			 outputFileName);

  mgr->ConnectInput  (task,  0, cinput );
  mgr->ConnectOutput (task,  1, coutput);

  //================================================
  // -- data containers - SPD
  
#ifdef USE_SPD
  AliAnalysisDataContainer *coutputSPD = 
    mgr->CreateContainer("HIMultCorr_SPD", TList::Class(),
                           AliAnalysisManager::kOutputContainer,
			                  outputFileName);

  mgr->ConnectInput  (taskSPD,  0, cinput );
  mgr->ConnectOutput (taskSPD,  1, coutputSPD);
#endif

  return task;
}
