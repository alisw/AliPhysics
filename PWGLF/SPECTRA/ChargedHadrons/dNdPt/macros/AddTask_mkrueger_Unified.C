AlidNdPtUnifiedAnalysisTask* AddTask_mkrueger_Unified(
  Int_t nBinsMultiplicity = 100, 
  Float_t etaCut = 0.8, 
  Float_t upperPtCut = 10., 
  Bool_t is2015Data = kTRUE, 
  Bool_t oldTrigger = kFALSE, 
  Bool_t isPbPbAnalysis = kFALSE, 
  Bool_t includeSigmas = kTRUE)
{

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  if (!mgr) {
    Error("AddTask_mkrueger_Unified", "No analysis manager found.");
    return 0;
  }

  // Switch off all AliInfo (too much output!!!)
  AliLog::SetGlobalLogLevel(AliLog::kError);
  mgr->SetDebugLevel(0);


  TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  Bool_t hasMC=(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()!=0x0);

  //AliESDtrackCuts* esdTrackCuts = new AliESDtrackCuts("AliESDtrackCuts");

  AlidNdPtUnifiedAnalysisTask *task = new AlidNdPtUnifiedAnalysisTask("AlidNdPtUnifiedAnalysisTask_mkrueger");
  task->SetUseMC(hasMC);
  if(type.Contains("ESD")) task->SetUseESD();
  else task->SetUseAOD();
  task->SetUseMultiplicity();
  task->SetUseCountedMult();
  task->SetIncludeSigmas(includeSigmas);	// to cross-check effects of particle composition correction

  if(isPbPbAnalysis) task->SetCentralityCut(60., 80.);
  
  task->SelectCollisionCandidates(AliVEvent::kINT7);
  task->SetTriggerMask(AliVEvent::kINT7 );
  
  if(oldTrigger){
    task->SelectCollisionCandidates(AliVEvent::kMB);
    task->SetTriggerMask(AliVEvent::kMB);    
  }
    

  const Int_t multNbins = nBinsMultiplicity;
  Double_t binsMult[multNbins+1];
  for (Int_t ii = 0; ii <= multNbins; ii++){binsMult[ii] = ii-0.5;}
  task->SetBinsMultCent(multNbins,binsMult);

  /// Acceptance cuts for tracks

  task->SetMinEta(-etaCut);
  task->SetMaxEta(etaCut);
  task->SetMinPt(0.15);
  task->SetMaxPt(upperPtCut);

  task->Set2013pA(kFALSE);
  task->Set2015data(is2015Data);

  ///TOF pileup, kTRUE only for Matching efficiency studies
  task->SetTOFbunchCrossing(kFALSE);

  task->SetMeanXYZv(0.0,0.0,0.0);
  task->SetSigmaMeanXYZv(1.0,1.0,10.0);
  task->SetZvtx(10.); //30
  task->SetEventTriggerRequired(kTRUE);

  // Quality cuts for tracks
  task->SetTPCRefit(kTRUE);
  task->SetITSRefit(kTRUE);
  task->SetKinkDaughters(kFALSE);
  task->SetRatioCrossedRowsOverFindableClustersTPC(0.8);
  task->SetFractionSharedClustersTPC(0.4);
  task->SetMaxchi2perTPCclu(4.);
  task->SetClusterReqITS(kTRUE);
  task->SetMaxchi2perITSclu(36.);
  task->SetDCAtoVertex2D(kFALSE);
  task->SetSigmaToVertex(kFALSE);
  task->SetDCAtoVertexZ(2.0);
  task->SetDCAtoVertexXYPtDep("0.0182+0.0350/pt^1.01");
  // task->SetMaxChi2TPCConstrained(36.);
  task->SetMinLenghtInActiveZoneTPC(0);
  task->SetGeometricalCut(kTRUE,3,130,1.5,0.85,0.7); ///if kTRUE comment CrossedRowsTPC cut
  // task->SetMinCrossedRowsTPC(120);

  mgr->AddTask(task);

  char containerName[60] = ""; 
  sprintf(containerName, "mkrueger_dNdPt_mult_%d_eta_%.2f_ptMax_%.2f", nBinsMultiplicity, etaCut, upperPtCut);
  
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput = mgr->CreateContainer(containerName,
							   TList::Class(),
							   AliAnalysisManager::kOutputContainer,
							   "AnalysisResults.root");



  mgr->ConnectInput(task, 0, cinput);
  mgr->ConnectOutput(task, 1, coutput);

  return task;

}
