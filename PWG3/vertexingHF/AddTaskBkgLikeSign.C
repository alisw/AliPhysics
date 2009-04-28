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
  AliAnalysisTaskSEBkgLikeSignJPSI *hfTask = new AliAnalysisTaskSEBkgLikeSignJPSI("CmpLikeSignAnalysis");
  hfTask->SetDebugLevel(2);

  mgr->AddTask(hfTask);

  //
  // Create containers for input/output
  mgr->ConnectInput(hfTask,0,mgr->GetCommonInputContainer());

  AliAnalysisDataContainer *coutput = mgr->CreateContainer("coutputLikeSign",TList::Class(),
                                                           AliAnalysisManager::kOutputContainer,
                                                           "CmpLikeSignJPSI.root");

  mgr->ConnectOutput(hfTask,1,coutput);


  return hfTask;
}
