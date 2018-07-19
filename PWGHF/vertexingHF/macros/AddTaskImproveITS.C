AliAnalysisTaskSEImproveITS *AddTaskImproveITS(Bool_t isRunInVertexing=kFALSE, // set to kTRUE to run during AODvertexingHF creation
					       const char *period="LHC15o",
                           const char *systematic="central",
					       Int_t ndebug=0) {
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddAliAnalysisTaskSEImproveITS", "No analysis manager to connect to.");
    return 0;
  }
                
  AliAnalysisTaskSEImproveITS *task
    =new AliAnalysisTaskSEImproveITS("ITSImprover",
                                     period,
                                     systematic,
                                     isRunInVertexing,
                                     ndebug);
  mgr->AddTask(task);
                
  TString outputFileName=AliAnalysisManager::GetCommonFileName();
  outputFileName+=":ITSImprover";
  AliAnalysisDataContainer *coutput
     =mgr->CreateContainer("debug",
                           TNtuple::Class(),
                           AliAnalysisManager::kOutputContainer,
                           outputFileName);
  
  mgr->ConnectInput (task,0,mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task,1,coutput);

  return task;
}

