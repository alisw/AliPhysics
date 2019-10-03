void AddTask_jacek_dNdPtTrackDumpTaskPbPb_TPCITS()
{
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libTender");
  gSystem->Load("libCORRFW");
  gSystem->Load("libPWGUDbase");
  gSystem->Load("libTPCcalib");
  gSystem->Load("libPWGPP");
  gSystem->Load("libPWGLFspectra");


  gRandom->SetSeed(0);

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  if (!mgr) {
    Error("AddTask_dNdPtTrackDumpTaskPbPb_TPCITS", "No analysis manager found.");
    return 0;
  }

  // Switch off all AliInfo (too much output!!!)
  AliLog::SetGlobalLogLevel(AliLog::kError);
  mgr->SetDebugLevel(0);

  

  //
  // Create physics trigger selection class
  //
  AliPhysicsSelection *physTrigSel =  new AliPhysicsSelection();

  //
  // Create event cuts
  //
  Float_t zvWindow = 30. ;

  AlidNdPtEventCuts *evtCuts = new AlidNdPtEventCuts("AlidNdPtEventCuts","Event cuts");
  evtCuts->SetZvRange(-zvWindow,zvWindow);
  evtCuts->SetMeanXYZv(0.0,0.0,0.0);
  evtCuts->SetSigmaMeanXYZv(1.0,1.0,10.0);
  evtCuts->SetTriggerRequired(kFALSE);
  //evtCuts->SetTriggerRequired(kTRUE);

  //
  // Create geom. acceptance cuts
  //
  Float_t etaWindow = 1.0 ;
  Float_t ptMin = 0.15 ;

  AlidNdPtAcceptanceCuts *accCuts = new AlidNdPtAcceptanceCuts("AlidNdPtAcceptanceCuts","Geom. acceptance cuts");
  accCuts->SetEtaRange(-etaWindow,etaWindow);
  accCuts->SetPtRange(ptMin,1.e10);
  accCuts->SetMaxDCAr(3.0);
  accCuts->SetMaxDCAz(30.0);

  //
  // Create standard esd track cuts
  //
  Int_t cutMode = 154;
  //Int_t cutMode = 156;
  //Int_t cutMode = 200;

  gROOT->LoadMacro("$ALICE_ROOT/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/CreatedNdPtTrackCuts.C");
  AliESDtrackCuts* esdTrackCuts = CreatedNdPtTrackCuts(cutMode);
  if (!esdTrackCuts) {
    printf("ERROR: esdTrackCuts could not be created\n");
    return;
  } else {
    esdTrackCuts->SetHistogramsOn(kTRUE);
    esdTrackCuts->SetRequireITSRefit(kFALSE); // no ITS refit
  }

  Bool_t hasMC=(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()!=0x0);

  //
  // Create task
  //
  AlidNdPtTrackDumpTask *task = new AlidNdPtTrackDumpTask("AlidNdPtTrackDumpTask");
  task->SetUseMCInfo(hasMC);
  //task->SetLowPtTrackDownscaligF(1.e4);
  //task->SetLowPtV0DownscaligF(1.e2);
  task->SetLowPtTrackDownscaligF(1.e7);
  task->SetLowPtV0DownscaligF(1.e4);
  task->SetProcessAll(kTRUE);
  task->SetProcessCosmics(kTRUE);
  //task->SetProcessAll(kFALSE);

  // trigger
  //task->SelectCollisionCandidates(AliVEvent::kMB); 

  //
  // set analysis options from the Helper here !!!
  //
  // AlidNdPtHelper::OutputObject outputObject = AlidNdPtHelper::kCutAnalysisPbPb;
  // AlidNdPtHelper::ParticleMode particleMode = AlidNdPtHelper::kAllPart ;
  
  AlidNdPtHelper::AnalysisMode analysisMode = AlidNdPtHelper::kTPCITS;

  task->SetUseMCInfo(hasMC);
  task->SetEventCuts(evtCuts);
  task->SetAcceptanceCuts(accCuts);
  task->SetTrackCuts(esdTrackCuts);
  task->SetAnalysisMode(analysisMode); 
  task->SetCentralityEstimator("V0M");
    
  // Add task
  mgr->AddTask(task);

  // Create containers for input
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  mgr->ConnectInput(task, 0, cinput);

  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("filtered1", TTree::Class(), AliAnalysisManager::kOutputContainer, "FilterEvents_Trees.root");
  mgr->ConnectOutput(task, 1, coutput1);
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("filtered2", TTree::Class(), AliAnalysisManager::kOutputContainer, "FilterEvents_Trees.root");
  mgr->ConnectOutput(task, 2, coutput2);
  AliAnalysisDataContainer *coutput3 = mgr->CreateContainer("filtered3", TTree::Class(), AliAnalysisManager::kOutputContainer, "FilterEvents_Trees.root");
  mgr->ConnectOutput(task, 3, coutput3);
  AliAnalysisDataContainer *coutput4 = mgr->CreateContainer("filtered4", TTree::Class(), AliAnalysisManager::kOutputContainer, "FilterEvents_Trees.root");
  mgr->ConnectOutput(task, 4, coutput4);
  AliAnalysisDataContainer *coutput5 = mgr->CreateContainer("filtered5", TTree::Class(), AliAnalysisManager::kOutputContainer, "FilterEvents_Trees.root");
  mgr->ConnectOutput(task, 5, coutput5);
  AliAnalysisDataContainer *coutput6 = mgr->CreateContainer("filtered6", TTree::Class(), AliAnalysisManager::kOutputContainer, "FilterEvents_Trees.root");
  mgr->ConnectOutput(task, 6, coutput6);
}

