AliAnalysisTaskSEBkgLikeSignJPSI *AddTaskBkgLikeSign() 
{
  //
  // Test macro for the AliAnalysisTaskSEBtoJPSItoEle 
  // starting from AliAOD.root file with HF + Like Sign candidates.
  // C.Di Giglio, carmelo.digiglio@ba.infn.it
  //


  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskCompareHF", "No analysis manager to connect to.");
    return NULL;
  }   

  // Like-sign background analysis task    
  AliAnalysisTaskSEBkgLikeSignJPSI *lsTask = new AliAnalysisTaskSEBkgLikeSignJPSI("CmpLikeSignAnalysis");
  lsTask->SetDebugLevel(2);

  mgr->AddTask(lsTask);

  //
  // Create containers for input/output
  AliAnalysisDataContainer *cinputLS = mgr->CreateContainer("cinputLikeSign",TChain::Class(), 
							  AliAnalysisManager::kInputContainer);
  AliAnalysisDataContainer *coutputLS = mgr->CreateContainer("coutputLikeSign",TList::Class(),
                                                           AliAnalysisManager::kOutputContainer,
                                                           "CmpLikeSignJPSI.root");

  mgr->ConnectInput(lsTask,0,mgr->GetCommonInputContainer());
  mgr->ConnectOutput(lsTask,1,coutputLS);


  return lsTask;
}
