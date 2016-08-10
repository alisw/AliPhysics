AliAnalysisTaskPSHFE* AddTaskPSHFE(const char* taskname, Bool_t trkCutsStrong=kTRUE){
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskEmcalJetSample", "No analysis manager to connect to.");
    return NULL;
  }

  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler())
  {
    ::Error("AddTaskEmcalJetSample", "This task requires an input event handler");
    return NULL;
  }
  
  AliAnalysisTaskPSHFE* PSHFEtask = new AliAnalysisTaskPSHFE(taskname);
    
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/macros/CreateTrackCutsPWGJE.C");

  AliESDtrackCuts* globaltrackCuts = 0x0;
  AliESDtrackCuts* comptrackCuts = 0x0;
    
  globaltrackCuts = CreateTrackCutsPWGJE(10001006);
  comptrackCuts = CreateTrackCutsPWGJE(10011008);
    
  PSHFEtask->SetTrackCuts(globaltrackCuts, comptrackCuts);
  PSHFEtask->SetElectronTrackCuts(trkCutsStrong);
  mgr->AddTask(PSHFEtask);  
    
  TString outfilename = "PSHFE_histos";
  
  // create containers for input/output
    AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
    
    AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("Min-Bias", TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s", AliAnalysisManager::GetCommonFileName()));
    
    AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("EMCal7", TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s", AliAnalysisManager::GetCommonFileName()));
    
    AliAnalysisDataContainer *coutput3 = mgr->CreateContainer("EMCal8", TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s", AliAnalysisManager::GetCommonFileName()));
    
    AliAnalysisDataContainer *coutput4 = mgr->CreateContainer("EMCalJet", TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s", AliAnalysisManager::GetCommonFileName()));
    
    
    // connect input/output
    mgr->ConnectInput(PSHFEtask, 0, cinput);
    mgr->ConnectOutput(PSHFEtask, 1, coutput1);
    mgr->ConnectOutput(PSHFEtask, 2, coutput2);
    mgr->ConnectOutput(PSHFEtask, 3, coutput3);
    mgr->ConnectOutput(PSHFEtask, 4, coutput4);
    
return PSHFEtask;
}