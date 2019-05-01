void AddTask_dNdPtPbPb_CutAnalysis(Int_t cutMode =222 , char *controlString ="default", char* eventTrigger="kINT7"){

  TString stEventTrigger(eventTrigger);
  TString stControlString(controlString);
  
  gSystem->Load("libPWG0base");
  gSystem->Load("libPWG0dep");
  gSystem->Load("libPWG0selectors");

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {Error("AddTask_dNdPtPbPb_CutAnalysis", "No analysis manager found.");return 0;}

  Bool_t hasMC=(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()!=0x0);

  // Switch off all AliInfo (too much output!!!)
  AliLog::SetGlobalLogLevel(AliLog::kError);
  mgr->SetDebugLevel(0);


  /// Create event cuts
  Float_t zvWindow = 10.;

  AlidNdPtEventCuts *evtCuts = new AlidNdPtEventCuts("AlidNdPtEventCuts","Event cuts");
  evtCuts->SetZvRange(-zvWindow,zvWindow);
  evtCuts->SetMeanXYZv(0.0,0.0,0.0);
  evtCuts->SetSigmaMeanXYZv(1.0,1.0,10.0);
  evtCuts->SetTriggerRequired(kTRUE);

  // Create geom. acceptance cuts
  Float_t etaWindow = 0.8;
  Float_t ptMin = 0.0;    
  AlidNdPtAcceptanceCuts *accCuts = new AlidNdPtAcceptanceCuts("AlidNdPtAcceptanceCuts","Geom. acceptance cuts");
  accCuts->SetEtaRange(-etaWindow,etaWindow);
  accCuts->SetPtRange(ptMin,1.e10);

  // Create standard esd track cuts
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/CreatedNdPtTrackCuts.C");
  AliESDtrackCuts* esdTrackCuts = CreatedNdPtTrackCuts(cutMode,stControlString);
  if (!esdTrackCuts) { printf("ERROR: esdTrackCuts could not be created\n"); return; }
  esdTrackCuts->SetHistogramsOn(kFALSE);

  // Create task
  //
  AliAnalysisTaskCutTest *task = new AliAnalysisTaskCutTest("cuttestANDresolution");
  task->SetUseMCInfo(hasMC);

  // trigger selection: MB
  if(stEventTrigger.Contains("kINT7")) task->SelectCollisionCandidates(AliVEvent::kINT7);
  else if(stEventTrigger.Contains("kMB")) task->SelectCollisionCandidates(AliVEvent::kMB);

  task->SetTrackCuts(esdTrackCuts);
  task->SetAcceptanceCuts(accCuts);
  task->SetEventCuts(evtCuts);

 mgr->AddTask(task);


  // Create containers for input
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  mgr->ConnectInput(task, 0, cinput);

  TString stContainerName;
  stContainerName = Form("dNdPt_PbPb_CutTest_%d",cutMode);
  if(!stControlString.Contains("default")) {
    stContainerName = stContainerName + "_" + stControlString;
  }

  AliAnalysisDataContainer *coutput = mgr->CreateContainer(stContainerName,TList::Class(),AliAnalysisManager::kOutputContainer, mgr->GetCommonFileName());

  mgr->ConnectOutput(task, 1, coutput);
}

