AliAnalysisTaskITSTrackingCheck *AddTaskPerformanceITS(Bool_t readMC=kFALSE,
						       Bool_t readRP=kFALSE,
						       Bool_t fillNtuples=kFALSE,
						       Int_t minmult=0,
						       Int_t maxmult=1000000,
						       Int_t checkSDDIsIn=1,
						       TString suffix="") 
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
  taskITS->SetMultiplicityRange(minmult,maxmult);
  taskITS->SetReadMC(readMC);
  taskITS->SetReadRPLabels(readRP);
  taskITS->SetFillNtuples(fillNtuples);
  taskITS->SetUseITSSAforNtuples(kFALSE);
  taskITS->SetCheckSDDIsIn(checkSDDIsIn);
  //taskITS->SetOCDBPath("alien://folder=/alice/data/2011/OCDB"); // to be commented for the QAtrain
  AliLog::SetClassDebugLevel("AliAnalysisTaskITSTrackingCheck",10);
  // Add to the manager
  mgr->AddTask(taskITS);

  //
  // Create containers for input/output
  TString cname="cOutputITS";
  cname+=suffix.Data();
  if(maxmult<1000000) {
    cname.Append("_"); cname+=minmult; 
    cname.Append("_"); cname+=maxmult;
  } 


  AliAnalysisDataContainer *cOutputITS = mgr->CreateContainer(cname.Data(),TList::Class(),AliAnalysisManager::kOutputContainer,Form("%s:ITS_Performance",mgr->GetCommonFileName()));


  // Attach input
  mgr->ConnectInput(taskITS,0,mgr->GetCommonInputContainer());
  // Attach output
  mgr->ConnectOutput(taskITS, 1,cOutputITS);
  
  return taskITS;
}
