AliAnalysisTaskSEBkgLikeSignD0 *AddTaskBkgLikeSignD0() 
{
  //
  // Test macro for the AliAnalysisTaskSEBkgLikeSignD0 
  // starting from AliAOD.root file with HF + Like Sign candidates.
  // C.Di Giglio, carmelo.digiglio@ba.infn.it
  //


  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskBkgLikeSignD0", "No analysis manager to connect to.");
    return NULL;
  }   

  // Like-sign background analysis task    
  AliAnalysisTaskSEBkgLikeSignD0 *lsD0Task = new AliAnalysisTaskSEBkgLikeSignD0("CmpLikeSignD0Analysis");
  lsD0Task->SetDebugLevel(0);

  mgr->AddTask(lsD0Task);

  //
  // Create containers for input/output
  AliAnalysisDataContainer *cinputLSD0 = mgr->CreateContainer("cinputLikeSignD0",TChain::Class(), 
							  AliAnalysisManager::kInputContainer);
  TString outputfile = AliAnalysisManager::GetCommonFileName();
  outputfile += ":PWG3_D2H_CmpLikeSignD0";
  AliAnalysisDataContainer *coutput1LSD0 = mgr->CreateContainer("coutput1LikeSignD0",TList::Class(),
                                                           AliAnalysisManager::kOutputContainer,
								outputfile.Data());
  AliAnalysisDataContainer *coutput2LSD0 = mgr->CreateContainer("coutput2LikeSignD0",TH1F::Class(),
                                                           AliAnalysisManager::kOutputContainer,
								outputfile.Data());

  mgr->ConnectInput(lsD0Task,0,mgr->GetCommonInputContainer());
  mgr->ConnectOutput(lsD0Task,1,coutput1LSD0);
  mgr->ConnectOutput(lsD0Task,2,coutput2LSD0);

  return lsD0Task;
}
