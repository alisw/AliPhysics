AliAnalysisTaskPHOSCluster* AddTaskPHOSCluster (Float_t fZVertex  = 10.0,
																Int_t   fMinCells = 3,
																Bool_t  fHitmapFilling = kTRUE,
															   TString tskname   = "AliPHOSCluster") {
  // Get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskPHOSCluster", "No analysis manager found.");
    return 0x0;
  }

  // Set up the task
  AliAnalysisTaskPHOSCluster* task = new AliAnalysisTaskPHOSCluster(tskname.Data());
  task ->SetZVertexCut(fZVertex);
  task ->SetMinCellsPerCluster(fMinCells);
  task ->SetHitmapFillingCellByCell(fHitmapFilling);

  // Add task to the manager
  mgr->AddTask(task);

  // Create output container
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(Form("%s_summary", tskname.Data()), TList::Class(),
                                                            AliAnalysisManager::kOutputContainer,
                                                            "AnalysisResults.root");


  // Connect input and output
  mgr->ConnectInput  (task,  0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput (task,  1, coutput1);
  return task;
}
