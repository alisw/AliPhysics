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

  gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/CreatedNdPtTrackCuts.C");
  AliESDtrackCuts *esdTrackCuts0 = CreatedNdPtTrackCuts(200);
  AliESDtrackCuts *esdTrackCuts1 = CreatedNdPtTrackCuts(42);
  AliESDtrackCuts *esdTrackCuts2 = CreatedNdPtTrackCuts(72);
  

 // --- End track Cuts ---

  task->AddAliESDtrackCut(esdTrackCuts0, "TPC-only tracking, cut number 42");
  task->AddAliESDtrackCut(esdTrackCuts1, "TPC+ITS combine tracking + DCAr(pt) (2011)");
  task->AddAliESDtrackCut(esdTrackCuts2, "TPC+ITS combine tracking + DCAr(pt) (2010)");
  task->SelectCollisionCandidates(AliTrigger::kMB);
  task->SetMaxVertexZ(10.);    // cm


  mgr->AddTask(task);
  
  
  AliAnalysisDataContainer *coutput = mgr->CreateContainer(Form("%s", combinedName.Data()), TList::Class(),  AliAnalysisManager::kOutputContainer, Form("%s:MeanPtRaw", mgr->GetCommonFileName()));
  
  mgr->ConnectInput( task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, coutput);
  

  return task;
}


