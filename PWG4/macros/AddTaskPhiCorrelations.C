AliAnalysisTaskPhiCorrelations *AddTaskPhiCorrelations(Int_t analysisMode = 0, Bool_t ppRun = kFALSE, const char* outputFileName = 0, Bool_t eventMixing = kTRUE)
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
  
//   Int_t bit = 1;
  Int_t bit = 128;
//   Int_t bit = 256;
//   Int_t bit = 512;
//   Int_t bit = 1024;
  ana->SetFilterBit(bit);  
  
  Printf("AddTaskPhiCorrelations:\n\n\n++++++++++ Using bit %d ++++++++++++\n\n\n", bit);
  
  ana->SetTrackEtaCut(1.0);
//   ana->SetTrackEtaCut(0.8);
//   ana->SetPtMin(0.15);
  ana->SetPtMin(1.0);
  //ana->SetEventSelectionBit(AliAnalysisHelperJetTasks::kIsPileUp);
  ana->SetReduceMemoryFootprint(kTRUE);
  //ana->SetSelectCharge(2);
  
  ana->SetEventMixing(eventMixing);
  
//   ana->SetCompareCentralities(kTRUE);
//   ana->SetTwoTrackEfficiencyStudy(kTRUE);
  
//   ana->SetFillpT(kTRUE);
  
  if (0)
  {
    Printf("AddTaskPhiCorrelations:\n\n\n++++++++++ Using SPD centrality selection ++++++++++++\n\n\n");
    ana->SetCentralityMethod("CL1");
  }    
  
  if (0)
  {
    Printf("AddTaskPhiCorrelations:\n\n\n++++++++++ Using ZDC centrality selection ++++++++++++\n\n\n");
    ana->SetCentralityMethod("ZEMvsZDC");
  }    
  

  if (ppRun)
  {
    Printf("AddTaskPhiCorrelations:\n\n\n+++++++++++++++ Configuring for p+p! +++++++++++++++++\n\n\n");
    ana->SetCentralityMethod(""); // only for pp
  }    
  
//   gSystem->Sleep(3000);
  
  mgr->AddTask(ana);
  
  // Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================
  if (!outputFileName)
    outputFileName = AliAnalysisManager::GetCommonFileName();
  
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("histosPhiCorrelations", TList::Class(),AliAnalysisManager::kOutputContainer,Form("%s:PWG4_PhiCorrelations", outputFileName));
  
  mgr->ConnectInput  (ana, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput (ana, 0, coutput1 );
   
  return ana;
}
