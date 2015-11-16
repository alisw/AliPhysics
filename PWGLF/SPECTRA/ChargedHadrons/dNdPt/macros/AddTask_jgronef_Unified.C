AlidNdPtUnifiedAnalysisTask* AddTask_jgronef_Unified()
{

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  if (!mgr) {
    Error("AddTask_jgronef_Unified", "No analysis manager found.");
    return 0;
  }

  // Switch off all AliInfo (too much output!!!)
  AliLog::SetGlobalLogLevel(AliLog::kError);
  mgr->SetDebugLevel(0);


  TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  Bool_t hasMC=(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()!=0x0);

  //AliESDtrackCuts* esdTrackCuts = new AliESDtrackCuts("AliESDtrackCuts");

  AlidNdPtUnifiedAnalysisTask *task = new AlidNdPtUnifiedAnalysisTask("AlidNdPtUnifiedAnalysisTask_jgronef");  
  task->SetUseMC(hasMC);
  if(type.Contains("ESD")) task->SetUseESD();
  else task->SetUseAOD();
  task->SetUseMultiplicity();
  task->SelectCollisionCandidates(AliVEvent::kMB);
  if(hasMC) task->SetMCParticleType(AlidNdPtUnifiedAnalysisTask::ParticleType::kPrimary); //only if MC, particle dependent MC analysis
  //task->SelectCollisionCandidates(AliVEvent::kINT7);
  
  Int_t multNbins = 252;  
  Double_t binsMult[253];
  for (int i=0; i<=multNbins; i++) { binsMult[i] = -0.5 + i; }
  binsMult[252] = 1000.;  
  task->SetBinsMultCent(multNbins,binsMult);
  
  

  /// Acceptance cuts for tracks
  
  task->SetMinEta(-0.8);
  task->SetMaxEta(0.8);
  task->SetMinPt(0.10);
  
  task->SetMeanXYZv(0.0,0.0,0.0);
  task->SetSigmaMeanXYZv(1.0,1.0,10.0);
  task->SetZvtx(10.);
  task->SetEventTriggerRequired(kTRUE);

  // Quality cuts for tracks
  task->SetTPCRefit(kTRUE);
  task->SetITSRefit(kTRUE);
  task->SetKinkDaughters(kFALSE);
  task->SetMinCrossedRowsTPC(120);
  task->SetRatioCrossedRowsOverFindableClustersTPC(0.8);
  task->SetFractionSharedClustersTPC(0.4);
  task->SetMaxchi2perTPCclu(4.);
  task->SetMaxchi2perITSclu(36.);
  task->SetDCAtoVertex2D(kFALSE);
  task->SetSigmaToVertex(kFALSE);
  task->SetDCAtoVertexZ(2.0);
//   task->SetMaxChi2TPCConstrained(36.);
  task->SetMinLenghtInActiveZoneTPC(0);

  mgr->AddTask(task);

  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput = mgr->CreateContainer("dNdPt",
							   TList::Class(),
							   AliAnalysisManager::kOutputContainer,
							   Form("%s:dNdPtHistos", mgr->GetCommonFileName())); 

  
  
  mgr->ConnectInput(task, 0, cinput);
  mgr->ConnectOutput(task, 1, coutput);
  
  return task; 

}

