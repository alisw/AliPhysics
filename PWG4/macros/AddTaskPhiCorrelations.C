AliAnalysisTaskPhiCorrelations *AddTaskPhiCorrelations(Int_t analysisMode = 0)
{
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskPhiCorrelations", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskPhiCorrelations", "This task requires an input event handler");
    return NULL;
  }
  
  // Create the task and configure it.
  //===========================================================================
  AliAnalysisTaskPhiCorrelations* ana = new  AliAnalysisTaskPhiCorrelations("PhiCorrelations");
  ana->SetMode(analysisMode);// data or corrections mode
  // common config,
  ana->SetDebugLevel(0); 
  //  ana->SetFilterBit(16);  
  //ana->SetFilterBit(64+32);  
  ana->SetFilterBit(1);  
  ana->SetTrackEtaCut(0.8);
  ana->SetPtMin(2.0);
  //ana->SetEventSelectionBit(AliAnalysisHelperJetTasks::kIsPileUp);
  ana->SetReduceMemoryFootprint(kTRUE);
  
  if (0)
  {
    file = TFile::Open("$ALICE_ROOT/PWG4/JetTasks/inputFiles/ue_trackingefficiency.root");
    trackingEff = (TH1D*) file->Get("trackingefficiency");
    ana->SetTrackingEfficiency(trackingEff);
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
