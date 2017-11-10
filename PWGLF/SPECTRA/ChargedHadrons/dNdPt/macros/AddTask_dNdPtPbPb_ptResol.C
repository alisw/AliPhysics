void AddTask_dNdPtPbPb_ptResol(Int_t cutMode =222 , Double_t smearing = 0.000, char* eventTrigger="kINT7"){

  TString stEventTrigger(eventTrigger);

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {Error("AddTask_dNdPtPbPb_ptResol", "No analysis manager found.");return 0;}

  Bool_t hasMC=(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()!=0x0);

  // Switch off all AliInfo (too much output!!!)
  AliLog::SetGlobalLogLevel(AliLog::kError);
  mgr->SetDebugLevel(0);


  /// Create event cuts
  Float_t zvWindow = 30. ;

  AlidNdPtEventCuts *evtCuts = new AlidNdPtEventCuts("AlidNdPtEventCuts","Event cuts");
  evtCuts->SetZvRange(-zvWindow,zvWindow);
  evtCuts->SetMeanXYZv(0.0,0.0,0.0);
  evtCuts->SetSigmaMeanXYZv(1.0,1.0,10.0);
  evtCuts->SetTriggerRequired(kTRUE);

  // Create geom. acceptance cuts
  Float_t etaWindow = 0.8 ;
  Float_t ptMin = 0.15;

  AlidNdPtAcceptanceCuts *accCuts = new AlidNdPtAcceptanceCuts("AlidNdPtAcceptanceCuts","Geom. acceptance cuts");
  accCuts->SetEtaRange(-etaWindow,etaWindow);
  accCuts->SetPtRange(ptMin,1.e10);
  accCuts->SetMaxDCAr(3.0);
  accCuts->SetMaxDCAz(30.0);

  // Create standard esd track cuts
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/CreatedNdPtTrackCuts.C");
  AliESDtrackCuts* esdTrackCuts = CreatedNdPtTrackCuts(cutMode);
  if (!esdTrackCuts) { printf("ERROR: esdTrackCuts could not be created\n"); return; }
  esdTrackCuts->SetHistogramsOn(kTRUE);

  // Create task
  AlidNdPtTask *task = new AlidNdPtTask("AlidNdPtTask");
  task->SetUseMCInfo(hasMC);

  // trigger selection: MB
  if(stEventTrigger.Contains("kINT7")) task->SelectCollisionCandidates(AliVEvent::kINT7);
  else if(stEventTrigger.Contains("kMB")) task->SelectCollisionCandidates(AliVEvent::kMB);
  // Create analysis object

  AliPtResolAnalysisPbPb *fdNdPtAnalysis = new AliPtResolAnalysisPbPb("PtResolAnalysis","pT Resolution Analysis");
  fdNdPtAnalysis->SetEventCuts(evtCuts);
  fdNdPtAnalysis->SetAcceptanceCuts(accCuts);
  fdNdPtAnalysis->SetTrackCuts(esdTrackCuts);
  fdNdPtAnalysis->SetAnalysisMode(AlidNdPtHelper::kTPCITS); 
  fdNdPtAnalysis->SetParticleMode(AlidNdPtHelper::kAllPart); 

  if(stEventTrigger.Contains("kINT7")) fdNdPtAnalysis->SetTriggerMask(AliVEvent::kINT7);
  else if(stEventTrigger.Contains("kMB")) fdNdPtAnalysis->SetTriggerMask(AliVEvent::kMB);

  fdNdPtAnalysis->SetSigmaScale(smearing);
  fdNdPtAnalysis->SetUseMCInfo(hasMC);
  fdNdPtAnalysis->SetCentralityEstimator("V0M");
 
  // Add analysis object
  task->AddAnalysisObject( fdNdPtAnalysis );
  // Add task
  mgr->AddTask(task);

  // Create containers for input
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  mgr->ConnectInput(task, 0, cinput);

  TString stContainerName;
  stContainerName = Form("dNdPt_PbPb_ptResol_%1.0f",smearing*10000);

  AliAnalysisDataContainer *coutput = mgr->CreateContainer(stContainerName,TList::Class(),AliAnalysisManager::kOutputContainer, mgr->GetCommonFileName());

  mgr->ConnectOutput(task, 1, coutput);
}

