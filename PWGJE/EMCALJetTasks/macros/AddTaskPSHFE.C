AliAnalysisTaskPSHFE* AddTaskPSHFE(const char* taskname, Bool_t trkCutsStrong=kFALSE, Bool_t SSCuts=kFALSE)
{
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
  PSHFEtask->SetSSCutBool(SSCuts);
  mgr->AddTask(PSHFEtask);  
    
  TString contname = "";
  
    if(trkCutsStrong){contname+="_Strong";}else{contname+="_Weak";}
    if(SSCuts){contname+="_SS";}else{contname+="_NoSS";}
  // create containers for input/output
    AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
    
    AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(Form("Min-Bias%s",contname), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s", AliAnalysisManager::GetCommonFileName()));
    
    AliAnalysisDataContainer *coutput2 = mgr->CreateContainer(Form("EMCal7%s",contname), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s", AliAnalysisManager::GetCommonFileName()));
    
    AliAnalysisDataContainer *coutput3 = mgr->CreateContainer(Form("EMCalJet%s",contname), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s", AliAnalysisManager::GetCommonFileName()));
    
    
    
    // connect input/output
    mgr->ConnectInput(PSHFEtask, 0, cinput);
    mgr->ConnectOutput(PSHFEtask, 1, coutput1);
    mgr->ConnectOutput(PSHFEtask, 2, coutput2);
    mgr->ConnectOutput(PSHFEtask, 3, coutput3);
    
return PSHFEtask;
}
