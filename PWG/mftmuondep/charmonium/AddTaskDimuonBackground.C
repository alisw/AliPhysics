//====================================================================================================================================================

AliAnalysisTaskDimuonBackground* AddTaskDimuonBackground() {
  
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskDimuonBackground", "No analysis manager to connect to.");
    return NULL;
  }  
  
  AliAnalysisTaskDimuonBackground *task = new AliAnalysisTaskDimuonBackground("AliAnalysisTaskDimuonBackground");

  // in cm, taken from Fig. 7.4 of the ITS Upgrade TDR, in the hypothesis of ~80 tracks participating to the vtx
  // task -> SetVtxResolutionITS(5.e-4, 5.e-4, 4.e-4);

  task -> SetVertexMode(AliAnalysisTaskDimuonBackground::kReconstructed);

  task -> SetMinTriggerMatch(1);
  task -> SetSingleMuonMinPt(1.0);
  task -> SetSingleMuonMinEta(-3.6);
  task -> SetSingleMuonMaxEta(-2.5);
  task -> SetSingleMuonMaxChi2(9999.);
  
  // create output container(s)
  AliAnalysisDataContainer *histogramList = mgr->CreateContainer("list", TList::Class(), AliAnalysisManager::kOutputContainer, "AnalysisMFT_Output.root");
  
  mgr->AddTask(task);

  // connect input and output
  mgr->ConnectInput(task,  0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, histogramList);

  return task;

}

//====================================================================================================================================================
