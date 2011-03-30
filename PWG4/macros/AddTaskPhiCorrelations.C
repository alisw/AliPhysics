AliAnalysisTaskPhiCorrelations *AddTaskPhiCorrelations(Int_t analysisMode = 0)
{
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskPhiCorrelations", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Create the task and configure it.
  //===========================================================================
  AliAnalysisTaskPhiCorrelations* ana = new  AliAnalysisTaskPhiCorrelations("PhiCorrelations");
  ana->SetMode(analysisMode);// data or corrections mode
  
//  if (analysisMode == 0) // data
//    ana->SelectCollisionCandidates(AliVEvent::kMB);

  // common config,
  ana->SetDebugLevel(0); 
  //  ana->SetFilterBit(16);  
  //ana->SetFilterBit(64+32);  
  ana->SetFilterBit(1);  
  ana->SetTrackEtaCut(0.8);
  ana->SetPtMin(0.15);
  //ana->SetEventSelectionBit(AliAnalysisHelperJetTasks::kIsPileUp);
  ana->SetReduceMemoryFootprint(kTRUE);
  //ana->SetSelectCharge(2);
  
  if (1)
  {
    Printf("\n\n\n+++++++++++++++ Configuring for p+p! +++++++++++++++++\n\n\n");
    ana->SetCentralityMethod(""); // only for pp
  }    
  
  if (0)
  {
    Printf("\n\n\n++++++++++ Using SPD centrality selection ++++++++++++\n\n\n");
    ana->SetCentralityMethod("CL1");
  }    
  
  mgr->AddTask(ana);
  
  // Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("histosPhiCorrelations", TList::Class(),AliAnalysisManager::kOutputContainer,Form("%s:PWG4_PhiCorrelations", AliAnalysisManager::GetCommonFileName()));
  
  mgr->ConnectInput  (ana, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput (ana, 0, coutput1 );
   
  return ana;
}
