AliTrackMatchingTPCITSCosmics *AddTaskTrackMatchingTPCITS() 
{
  //
  // Task for the extraction of the ITS alignment data (AliTrackPoints)
  //
  // andrea.dainese@pd.infn.it
  //


  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTask", "No analysis manager to connect to.");
    return NULL;
  }   

  // Create the task
  AliTrackMatchingTPCITSCosmics *taskMatch = new AliTrackMatchingTPCITSCosmics("matching");
  AliLog::SetClassDebugLevel("AliTrackMatchingTPCITSCosmics",10);
  //taskMatch->SetOnlySPDFO();
  taskMatch->SetGeometryFileName("alien:///alice/cern.ch/user/d/dainesea/geometry.root");
  // Add to the manager
  taskMatch->SetReadHLTESD(kTRUE);
  mgr->AddTask(taskMatch);

  //
  // Create containers for input/output
  AliAnalysisDataContainer *cInput = mgr->CreateContainer("cInput",TChain::Class(),AliAnalysisManager::kInputContainer);

  AliAnalysisDataContainer *cOutput = mgr->CreateContainer("cOutput",TList::Class(),AliAnalysisManager::kOutputContainer,"TPCITSMatching.root");


  // Attach input
  mgr->ConnectInput(taskMatch,0,mgr->GetCommonInputContainer());
  // Attach output
  mgr->ConnectOutput(taskMatch,0,cOutput);
  
  return taskMatch;
}
