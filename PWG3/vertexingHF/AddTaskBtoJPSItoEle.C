AliAnalysisTaskSEBtoJPSItoEle *AddTaskBtoJPSItoEle() 
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
    ::Error("AddTaskBtoJPSItoEle", "No analysis manager to connect to.");
    return NULL;
  }   

  // Cdf unbinned log-likelihood fit analysis task    
  AliAnalysisTaskSEBtoJPSItoEle *hfTask = new AliAnalysisTaskSEBtoJPSItoEle("CdfFitAnalysis");
  hfTask->SetDebugLevel(0);

  mgr->AddTask(hfTask);

  //
  // Create containers for input/output
  mgr->ConnectInput(hfTask,0,mgr->GetCommonInputContainer());

  AliAnalysisDataContainer *coutput = mgr->CreateContainer("coutputCdfFit",TList::Class(),
                                                           AliAnalysisManager::kOutputContainer,       
							   "CdfFit.root");
  mgr->ConnectOutput(hfTask,1,coutput);

  return hfTask;
}
