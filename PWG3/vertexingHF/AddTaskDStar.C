AliAnalysisTaskSEDStar *AddTaskDStar(Bool_t readMC=kTRUE)
{
  //
  // Test macro for the AliAnalysisTaskSE for D*+ candidates
  // invariant mass histogram and association with MC truth 
  // (using MC info in AOD)
  // Yifei Wang, yifei@pi0.physi.uni-heidelberg.de
  //


  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskDStar", "No analysis manager to connect to.");
    return NULL;
  }   
  TString filename = AliAnalysisManager::GetCommonFileName();
  filename += ":PWG3_D2H_DStar";
  
  // Aanalysis task    
  AliAnalysisTaskSEDStar *DStarTask = new AliAnalysisTaskSEDStar("DStarAnalysis");
  DStarTask->SetDebugLevel(0);
  DStarTask->SetReadMC(readMC);
  mgr->AddTask(DStarTask);
  
  //
  // Create containers for input/output
  AliAnalysisDataContainer *cinputDStar = mgr->CreateContainer("cinputDStar",TChain::Class(), 
							  AliAnalysisManager::kInputContainer);
  AliAnalysisDataContainer *coutputDStar = mgr->CreateContainer("coutputDStar",TList::Class(),
							   AliAnalysisManager::kOutputContainer, 
							   filename.Data());

  mgr->ConnectInput(DStarTask,0,mgr->GetCommonInputContainer());

  mgr->ConnectOutput(DStarTask,1,coutputDStar);

  return DStarTask;
}
