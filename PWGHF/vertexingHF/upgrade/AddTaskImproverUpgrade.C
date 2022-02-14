AliAnalysisTaskSEImproveITS3 *AddTaskImproverUpgrade(Bool_t isRunInVertexing=kFALSE, // set to kTRUE to run during AODvertexingHF creation
                                               const char *resfileCurURI="ITSgraphs_NewAll-X0.3-Res4um.root",
                                               const char *resfileUpgURI="cylIB_bp17_B05Tnew.root",
					       Bool_t isImproveDeuteron=kFALSE,
					       Int_t ndebug=0) {
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddAliAnalysisTaskSEImproveITS3", "No analysis manager to connect to.");
    return 0;
  }
            
  AliAnalysisTaskSEImproveITS3 *task
    =new AliAnalysisTaskSEImproveITS3("ITSImprover",
                                     resfileCurURI,
                                     resfileUpgURI,
                                     isRunInVertexing,
                                     isImproveDeuteron,
                                     ndebug);

  mgr->AddTask(task);

                
 TString outputFileName=AliAnalysisManager::GetCommonFileName();
 outputFileName+=":ITSImprover";
  

  AliAnalysisDataContainer *coutput
     =mgr->CreateContainer("debug",
                           TList::Class(),
                           AliAnalysisManager::kOutputContainer,
                           outputFileName.Data());
 
  mgr->ConnectInput (task,0,mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task,1,coutput);


  return task;
}

