AliAnalysisTaskSED0Mass *AddTaskD0Mass()
{
  //
  // Test macro for the AliAnalysisTaskSE for D0 candidates
  // invariant mass histogram and association with MC truth 
  // (using MC info in AOD)
  // C.Bianchin  chiara.bianchin@pd.infn.it
  //


  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskD0Mass", "No analysis manager to connect to.");
    return NULL;
  }   

  
  // Aanalysis task    
  AliAnalysisTaskSED0Mass *massD0Task = new AliAnalysisTaskSED0Mass("D0MassAnalysis");
  massD0Task->SetDebugLevel(2);
  mgr->AddTask(massD0Task);
  
  //
  // Create containers for input/output
  AliAnalysisDataContainer *cinputmassD0 = mgr->CreateContainer("cinputmassD0",TChain::Class(), 
							  AliAnalysisManager::kInputContainer);
  AliAnalysisDataContainer *coutputmassD0 = mgr->CreateContainer("coutputmassD0",TList::Class(),
							   AliAnalysisManager::kOutputContainer, 
							   "D0InvMass.root");
  mgr->ConnectInput(massD0Task,0,mgr->GetCommonInputContainer());

  mgr->ConnectOutput(massD0Task,1,coutputmassD0);

  return massD0Task;
}
