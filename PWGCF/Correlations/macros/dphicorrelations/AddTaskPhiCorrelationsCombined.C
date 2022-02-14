AliAnalysisTaskPhiCorrelations *AddTaskPhiCorrelationsCombined(Int_t analysisMode = 0, Bool_t ppRun = kFALSE, const char* outputFileName = 0, Bool_t eventMixing = kTRUE, Int_t zVtxAxis = 0, const char* containerName = "histosPhiCorrelations", const char* folderName = "PWG4_PhiCorrelations", const char* suffix = "")
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
  TString combinedName;
  if(suffix!="")
    combinedName.Form("%s_%s", containerName, suffix);
  else
    combinedName=containerName;
  
  AliAnalysisTaskPhiCorrelations* ana = new  AliAnalysisTaskPhiCorrelations(combinedName);
  ana->SetMode(analysisMode);// data or corrections mode
  
  ana->SetDebugLevel(0); 

  Int_t bit = 128;
  ana->SetFilterBit(bit);  
  
  Printf("AddTaskPhiCorrelations:\n\n\n++++++++++ Using bit %d ++++++++++++\n\n\n", bit);
  
  ana->SetTrackEtaCut(0.9);

  ana->SetPtMin(0.5);

  //ana->SetReduceMemoryFootprint(kTRUE);
  
  ana->SetEventMixing(eventMixing);
  ana->SetUseVtxAxis(zVtxAxis);
  
  ana->SetZVertex(10);
  
  if (ppRun)
  {
    Printf("AddTaskPhiCorrelations:\n\n\n+++++++++++++++ Configuring for p+p! +++++++++++++++++\n\n\n");
    ana->SetCentralityMethod(""); // only for pp
  }    
  
  mgr->AddTask(ana);
  
  // Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================
  if (!outputFileName)
    outputFileName = AliAnalysisManager::GetCommonFileName();
  
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(combinedName, TList::Class(),AliAnalysisManager::kOutputContainer,Form("%s:%s", outputFileName, folderName));
  
  mgr->ConnectInput  (ana, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput (ana, 0, coutput1 );
   
  return ana;
}
