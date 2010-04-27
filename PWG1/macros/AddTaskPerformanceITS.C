AliAnalysisTaskITSTrackingCheck *AddTaskPerformanceITS(Bool_t readMC=kFALSE,
						       Bool_t readRP=kFALSE,
						       Bool_t fillNtuples=kFALSE) 
{
  //
  // Task for check of ITS tracking
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
  AliAnalysisTaskITSTrackingCheck *taskITS = new AliAnalysisTaskITSTrackingCheck("ITStracking");
  taskITS->SetReadMC(readMC);
  taskITS->SetReadRPLabels(readRP);
  taskITS->SetFillNtuples(fillNtuples);
  taskITS->SetUseITSSAforNtuples(kFALSE);
  AliLog::SetClassDebugLevel("AliAnalysisTaskITSTrackingCheck",10);
  // Add to the manager
  mgr->AddTask(taskITS);

  //
  // Create containers for input/output
  AliAnalysisDataContainer *cOutputITS = mgr->CreateContainer("cOutputITS",TList::Class(),AliAnalysisManager::kOutputContainer,Form("%s:ITS_Performance",mgr->GetCommonFileName()));


  // Attach input
  mgr->ConnectInput(taskITS,0,mgr->GetCommonInputContainer());
  // Attach output
  mgr->ConnectOutput(taskITS, 1,cOutputITS);
  
  return taskITS;
}
