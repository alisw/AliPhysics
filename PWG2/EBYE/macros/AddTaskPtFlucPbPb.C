AliAnalysisTask *AddTaskPtFlucPbPb(){
  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_sheckel_PtFlucPbPb", "No analysis manager found.");
    return 0;
  }


  // --- Check for MC ---
  AliMCEventHandler  *mcH = static_cast<AliMCEventHandler*>(mgr->GetMCtruthEventHandler());


  //============= Set Task Name ===================
  TString taskName=("AliAnalysisTaskPtFlucPbPb.cxx+");
  //===============================================
  //            Load the task
  gROOT->LoadMacro(taskName.Data());
  if (gProof){
    TString taskSO=gSystem->pwd();
    taskSO+="/";
    taskSO+=taskName(0,taskName.First('.'))+"_cxx.so";
    gProof->Exec(Form("gSystem->Load(\"%s\")",taskSO.Data()),kTRUE);
  }


  //========= Add task to the ANALYSIS manager =====
  AliAnalysisTaskPtFlucPbPb *task = new AliAnalysisTaskPtFlucPbPb("sheckelTaskPtFlucPbPb");


  // --- Set ESD track Cuts ---

  AliESDtrackCuts *esdTrackCuts = new AliESDtrackCuts("esdTrackCuts"); // My pp track cuts (ITS+TPC)
  AliESDtrackCuts *esdCuts = new AliESDtrackCuts("esdCuts"); // TPC-only track cuts

  // -- My pp track cuts (ITS+TPC)
  // TPC
  esdTrackCuts->SetRequireTPCStandAlone(kFALSE);
  esdTrackCuts->SetMinNClustersTPC(70);
  esdTrackCuts->SetMaxChi2PerClusterTPC(4.0);
  esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
  esdTrackCuts->SetRequireTPCRefit(kTRUE);
  // ITS
  esdTrackCuts->SetRequireITSRefit(kTRUE);
  esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kAny);
  // DCA to vertex
  esdTrackCuts->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01"); // 7*(0.0026+0.0050/pt^1.01)
  esdTrackCuts->SetMaxDCAToVertexZ(2.);
  esdTrackCuts->SetDCAToVertex2D(kFALSE);
  esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
  // Pt range and acceptance in eta
  esdTrackCuts->SetPtRange(0.15,2.);
  esdTrackCuts->SetEtaRange(-0.8,0.8);


  // -- TPC-only track cuts
  // TPC
  esdCuts->SetRequireTPCStandAlone(kTRUE);
  esdCuts->SetMinNClustersTPC(70);
  esdCuts->SetMaxChi2PerClusterTPC(4.0);
  esdCuts->SetAcceptKinkDaughters(kFALSE);
  esdCuts->SetRequireTPCRefit(kFALSE);
  // ITS
  esdCuts->SetRequireITSRefit(kFALSE);
  // DCA to vertex
  esdCuts->SetMaxDCAToVertexXY(2.4); // cm
  esdCuts->SetMaxDCAToVertexZ(3.2);  // cm
  esdCuts->SetDCAToVertex2D(kTRUE);
  esdCuts->SetRequireSigmaToVertex(kFALSE);
  // Pt range and acceptance in eta
  esdCuts->SetPtRange(0.15,2.);
  esdCuts->SetEtaRange(-0.8,0.8);

  // --- End track Cuts ---


//   task->SetAliESDtrackCuts(esdTrackCuts);
  task->SetAliESDtrackCuts(esdCuts);

  task->SelectCollisionCandidates(AliVEvent::kMB);
  task->SetMaxVertexZ(10.);    // cm
  task->SetNContributors(1);
  task->SetUseCentrality(1);   // 0=off, 1=VZERO, 2=SPD(not implemented)

  if (mcH) task->SetMC(kTRUE);


  mgr->AddTask(task);


  //================================================
  //              data containers
  //================================================

  // input container
  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();

  // output container
  AliAnalysisDataContainer *coutput =
      mgr->CreateContainer("sheckel_PtFlucPbPb", TList::Class(),
                           AliAnalysisManager::kOutputContainer,"sheckel_PtFlucPbPb.root");

  // connect containers
  mgr->ConnectInput  (task,  0, cinput );
  mgr->ConnectOutput (task,  1, coutput );

  return task;

}
