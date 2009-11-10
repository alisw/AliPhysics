AliAnalysisTaskSEBkgLikeSignJPSI *AddTaskBkgLikeSignJPSI() 
{
  //
  // Test macro for the AliAnalysisTaskSEBkgLikeSignJPSI
  // starting from AliAOD.root file with HF + Like Sign candidates.
  // C.Di Giglio, carmelo.digiglio@ba.infn.it
  //


  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskBkgLikeSignJPSI", "No analysis manager to connect to.");
    return NULL;
  }   


  // Like-sign background analysis task    
  AliAnalysisTaskSEBkgLikeSignJPSI *lsTask = new AliAnalysisTaskSEBkgLikeSignJPSI("CmpLikeSignAnalysis");
  lsTask->SetDebugLevel(0);

  mgr->AddTask(lsTask);

  //
  // Create containers for input/output
  AliAnalysisDataContainer *cinputLS = mgr->CreateContainer("cinputLikeSignJPSI",TChain::Class(), 
							  AliAnalysisManager::kInputContainer);
  TString outputfile = AliAnalysisManager::GetCommonFileName();
  outputfile += ":PWG3_D2H_CmpLikesignJPSI";
  AliAnalysisDataContainer *coutputLS = mgr->CreateContainer("coutputLikeSignJPSI",TList::Class(),
                                                           AliAnalysisManager::kOutputContainer,
							     outputfile.Data());

  mgr->ConnectInput(lsTask,0,mgr->GetCommonInputContainer());
  mgr->ConnectOutput(lsTask,1,coutputLS);


  return lsTask;
}
