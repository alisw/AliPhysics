AliAnalysisTaskMeanPtRaw *AddTask_MeanPtRaw(char *contName = "MeanPtRaw", char *suffix = ""){
  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_MeanPtRaw", "No analysis manager found.");
    return 0;
  }

// Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTask_MeanPtRaw", "This task requires an input event handler");
    return NULL;
  }   
  TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"

  ::Info("AddTask_MeanPtRaw", Form("data type is: %s", type.Data()));

  TString combinedName;
  combinedName.Form("%s%s",contName, suffix);
  
  //========= Add task to the ANALYSIS manager =====
  AliAnalysisTaskMeanPtRaw *task = new AliAnalysisTaskMeanPtRaw("pluettigMeanPt");

  // --- Set ESD track Cuts ---
  AliESDtrackCuts* esdTrackCutsTPCOnly = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
//   AliESDtrackCuts* esdTrackCutsTPCITSDefault = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011();
//   AliESDtrackCuts* esdTrackCutsTPCITS120Rows = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011();
//   esdTrackCutsTPCITS120Rows->SetMinNCrossedRowsTPC(120);
  

 // --- End track Cuts ---

  task->AddAliESDtrackCut(esdTrackCutsTPCOnly, "TPC-only tracking");
//   task->AddAliESDtrackCut(esdTrackCutsTPCITSDefault, "TPC+ITS combined tracking, default");
//   task->AddAliESDtrackCut(esdTrackCutsTPCITS120Rows, "TPC+ITS combined tracking, +120 crossedRows");
  task->SelectCollisionCandidates(AliVEvent::kMB);
  task->SetMaxVertexZ(10.);    // cm


  mgr->AddTask(task);
  
  
  AliAnalysisDataContainer *coutput = mgr->CreateContainer(Form("%s", combinedName.Data()), TList::Class(),  AliAnalysisManager::kOutputContainer, Form("%s:MeanPtRaw", mgr->GetCommonFileName()));
  
  mgr->ConnectInput( task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, coutput);
  

  return task;
}

