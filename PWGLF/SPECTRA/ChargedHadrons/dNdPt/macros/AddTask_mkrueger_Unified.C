AlidNdPtUnifiedAnalysisTask* AddTask_mkrueger_Unified(TString controlstring, Int_t cutModeLow = 100, Int_t cutModeHigh = 120)
{
  // INFO: for MC use 4 different trains (cutModeLow, cutModeHigh) = (100,105) (105,111), (111,116), (116,120)

  // settings:
  UInt_t offlineTriggerMask = AliVEvent::kINT7;
  Bool_t is2015Data = kFALSE;
  Bool_t isPbPbAnalysis = kFALSE; 
  
  Int_t nBinsMultiplicity = 100; 
  Float_t etaCut = 0.8;
  if(controlstring.Contains("eta03")) etaCut = 0.3;  
  
  Float_t upperPtCut = 10.; 
  Bool_t includeSigmas = kTRUE;
  string colsys = "pp";

  
  if(controlstring.Contains("pp")){    
    if(controlstring.Contains("5TeV")) is2015Data = kTRUE;
    if(controlstring.Contains("7TeV")) offlineTriggerMask = AliVEvent::kMB;
  }
  if(controlstring.Contains("pPb"))  {nBinsMultiplicity = 150; colsys = "pPb";}
  if(controlstring.Contains("PbPb")) {nBinsMultiplicity = 500; isPbPbAnalysis = kTRUE; is2015Data = kTRUE; colsys = "PbPb";}
  if(controlstring.Contains("excludeSigmas")) includeSigmas = kFALSE;

  // Binning in Multiplicity
  const Int_t multNbins = nBinsMultiplicity;
  Double_t binsMult[multNbins+1];
  for (Int_t ii = 0; ii <= multNbins; ii++){binsMult[ii] = ii-0.5;}


  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  if (!mgr) {
    Error("AddTask_mkrueger_Unified", "No analysis manager found.");
    return 0;
  }

  // Switch off all AliInfo (too much output!!!)
  AliLog::SetGlobalLogLevel(AliLog::kError);
  mgr->SetDebugLevel(0);

  TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  Bool_t hasMC = (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler() != 0x0);


  AlidNdPtUnifiedAnalysisTask* mainTask = NULL;

  char taskName[100] = ""; 

  
  for(Int_t cutMode = cutModeLow; cutMode < cutModeHigh; cutMode++){
    sprintf(taskName, "mkrueger_%s_mult_%d_eta_%.2f_cutMode_%d", colsys.c_str(), multNbins, etaCut, cutMode);
    
    AlidNdPtUnifiedAnalysisTask *task = new AlidNdPtUnifiedAnalysisTask(taskName);

    // general cuts
    task->SelectCollisionCandidates(offlineTriggerMask);
    task->SetTriggerMask(offlineTriggerMask);
    if(isPbPbAnalysis) task->SetCentralityCut(60., 80.);

    task->SetUseMC(hasMC);
    if(type.Contains("ESD")) task->SetUseESD();
    else task->SetUseAOD();
    task->SetUseMultiplicity();
    task->SetUseCountedMult();
    task->SetIncludeSigmas(includeSigmas);		// to cross-check effects of particle composition correction
    task->SetBinsMultCent(multNbins, binsMult);

    // nominal cut-setting:

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
    task->SetZvtx(10.);
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
    task->SetMaxChi2TPCConstrained(36.); // was excluded for some reason?
    task->SetMinLenghtInActiveZoneTPC(0);
    task->SetGeometricalCut(kTRUE,3,130,1.5,0.85,0.7); ///if kTRUE comment CrossedRowsTPC cut
    // task->SetMinCrossedRowsTPC(120);
    
//MC1--    
      if(cutMode == 100) mainTask = task; // has no particular use... just to return task with nominal cut setting in macro
    
    // cut-variation:
      if(cutMode == 101) {task->SetMaxchi2perITSclu(25.);}
      if(cutMode == 102) {task->SetMaxchi2perITSclu(49.);}

      if(cutMode == 103) {task->SetMaxchi2perTPCclu(3); }
      if(cutMode == 104) {task->SetMaxchi2perTPCclu(5); }
//MC2--	
      if(cutMode == 105) {task->SetRatioCrossedRowsOverFindableClustersTPC(0.7);}
      if(cutMode == 106) {task->SetRatioCrossedRowsOverFindableClustersTPC(0.9);}
	
      if(cutMode == 107) {task->SetFractionSharedClustersTPC(0.2);}
      if(cutMode == 108) {task->SetFractionSharedClustersTPC(1.0);}
	
      if(cutMode == 109) {task->SetMaxChi2TPCConstrained(25.);}
      if(cutMode == 110) {task->SetMaxChi2TPCConstrained(49.);}
//MC3--	
      if(cutMode == 111) {task->SetDCAtoVertexXYPtDep("0.0104+0.0200/pt^1.01");}
      if(cutMode == 112) {task->SetDCAtoVertexXYPtDep("0.0260+0.0500/pt^1.01");}
	
      if(cutMode == 113) {task->SetDCAtoVertexZ(1.0);}
      if(cutMode == 114) {task->SetDCAtoVertexZ(5.0);}

      if(cutMode == 115) {task->SetClusterReqITS(kFALSE);}
//MC4--		  
      if(cutMode == 116) {task->SetGeometricalCut(kTRUE,3,120,1.5,0.85,0.7);}	
      if(cutMode == 117) {task->SetGeometricalCut(kTRUE,3,140,1.5,0.85,0.7);}	

      if(cutMode == 118) {task->SetGeometricalCut(kTRUE,4,130,1.5,0.85,0.7);}
      if(cutMode == 119) {task->SetGeometricalCut(kTRUE,2,130,1.5,0.85,0.7);}
      


  // hang task in train
    
    mgr->AddTask(task);

    AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
    AliAnalysisDataContainer *coutput = mgr->CreateContainer(taskName,
							    TList::Class(),
							    AliAnalysisManager::kOutputContainer,
							    "AnalysisResults.root");



    mgr->ConnectInput(task, 0, cinput);
    mgr->ConnectOutput(task, 1, coutput);

  
  }  
    
  
  return mainTask;

}
