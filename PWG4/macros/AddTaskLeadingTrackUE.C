

void ConfigTaskUE(AliAnalysisTaskLeadingTrackUE * ueana );          // common config, extend with different cases

AliAnalysisTaskLeadingTrackUE *AddTaskLeadingTrackUE(Int_t analysisMode = 0)
{
// Creates a jet fider task, configures it and adds it to the analysis manager.

   // Get the pointer to the existing analysis manager via the static access method.
   //==============================================================================
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      ::Error("AddTaskLeadingTrackUE", "No analysis manager to connect to.");
      return NULL;
   }  
   
   // Check the analysis type using the event handlers connected to the analysis manager.
   //==============================================================================
   if (!mgr->GetInputEventHandler()) {
      ::Error("AddTaskLeadingTrackUE", "This task requires an input event handler");
      return NULL;
   }

   // Create the task and configure it.
   //===========================================================================
   
   AliAnalysisTaskLeadingTrackUE* ueana = new  AliAnalysisTaskLeadingTrackUE("UEAnalysis_LeadingTrack");
   ueana->SetMode(analysisMode);// data or corrections mode
   ConfigTaskUE(ueana);

   mgr->AddTask(ueana);
   
   // Create ONLY the output containers for the data produced by the task.
   // Get and connect other common input/output containers via the manager as below
   //==============================================================================
   AliAnalysisDataContainer *coutput1_UE = 0;
   coutput1_UE = mgr->CreateContainer("histosLeadingTrackUE", TList::Class(),AliAnalysisManager::kOutputContainer,Form("%s:PWG4_LeadingTrackUE",AliAnalysisManager::GetCommonFileName()));
   
   mgr->ConnectInput  (ueana, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput (ueana,     0, coutput1_UE );
   
   return ueana;
}

void ConfigTaskUE(AliAnalysisTaskLeadingTrackUE * ueana){
  // common config,
  ueana->SetDebugLevel(0); 
  ueana->SetPtRangeInHist(100, 0., 100.);
  //  ueana->SetFilterBit(16);  
  ueana->SetFilterBit(64+32);  
  ueana->SetTrackEtaCut(0.8);
  ueana->SetLeadingTrackEtaCut(0.8);
  ueana->SetEventSelectionBit(AliAnalysisHelperJetTasks::kIsPileUp);
  ueana->SetReduceMemoryFootprint(kTRUE);
 
  if (1)
  {
    file = TFile::Open("$ALICE_ROOT/PWG4/JetTasks/inputFiles/ue_trackingefficiency.root");
    trackingEff = (TH1D*) file->Get("trackingefficiency");
    ueana->SetTrackingEfficiency(trackingEff);
  }
}
