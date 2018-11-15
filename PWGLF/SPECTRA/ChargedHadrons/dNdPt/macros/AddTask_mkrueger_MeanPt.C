AliMeanPtAnalysisTask* AddTask_mkrueger_MeanPt(TString controlstring, Int_t cutModeLow = 100, Int_t cutModeHigh = 101)
{
  // INFO: for MC use 4 different trains (cutModeLow, cutModeHigh) = (100,105), (105,109), (109,113), (113,116), (116,120)

  // settings:
  Bool_t includeCrosscheckHistos = kFALSE;
  Bool_t is2013pA = kFALSE;
  UInt_t offlineTriggerMask = AliVEvent::kINT7;
  Bool_t is2015Data = kFALSE;
  Bool_t isPbPbAnalysis = kFALSE;

  Int_t maxMultiplicity = 100;
  Float_t etaCut = 0.8;
  if(controlstring.Contains("eta03")) etaCut = 0.3;

  Float_t lowerPtCut = 0.15;
  Float_t upperPtCut = 10.;


  Bool_t includeSigmas = kTRUE;
  string colsys = "pp";

  const Int_t nMultSteps = 3;
  Int_t multSteps[] = {100, 0, 0};
  Int_t multBinWidth[] = {1,1,1};

  Int_t nBinsCent = 1;
  Double_t centBinEdgesDummy[2] = {0., 100.};
  Double_t centBinEdgesUsed[9] = {0., 5., 10., 20., 40., 60., 80., 90., 100.};
  Double_t* centBinEdges = centBinEdgesDummy;

  if(controlstring.Contains("pp")){
    if(controlstring.Contains("performanceHistos")) includeCrosscheckHistos = kTRUE;
    if(controlstring.Contains("5TeV")) is2015Data = kTRUE;
    if(controlstring.Contains("7TeV")) offlineTriggerMask = AliVEvent::kMB;
    if(controlstring.Contains("2TeV")) offlineTriggerMask = AliVEvent::kMB;
    if(controlstring.Contains("09TeV")) offlineTriggerMask = AliVEvent::kMB;
  }
  if(controlstring.Contains("XeXe"))  {
    colsys = "XeXe";
    multSteps[0] = 3500;  multBinWidth[0] = 1;
    multSteps[1] = 0;    multBinWidth[1] = 1;
    multSteps[2] = 0;    multBinWidth[2] = 1;
    nBinsCent = 8;
    centBinEdges = centBinEdgesUsed;
  }
  if(controlstring.Contains("pPb"))  {
    colsys = "pPb";
    multSteps[0] = 300;   multBinWidth[0] = 1;
    multSteps[1] = 0;     multBinWidth[1] = 1;
    multSteps[2] = 0;     multBinWidth[2] = 1;
    is2013pA = kTRUE;
  }
  if(controlstring.Contains("PbPb")) {
    isPbPbAnalysis = kTRUE;
    is2015Data = kTRUE;
    colsys = "PbPb";
    multSteps[0] = 4500;   multBinWidth[0] = 1;
    multSteps[1] = 0;    multBinWidth[1] = 1;
    multSteps[2] = 0;    multBinWidth[2] = 1;
    nBinsCent = 8;
    centBinEdges = centBinEdgesUsed;
    if(controlstring.Contains("2TeV")){offlineTriggerMask = AliVEvent::kMB; is2015Data = kFALSE;}
  }
  if(controlstring.Contains("excludeSigmas")) includeSigmas = kFALSE;
  if(controlstring.Contains("oldTrigger")){offlineTriggerMask = AliVEvent::kMB;}

  // Binning in Multiplicity
  const Int_t nBinsMult = multSteps[0]+ multSteps[1] + multSteps[2] + 1;
  Double_t multBinEdges[nBinsMult+1];
  multBinEdges[0] = -0.5;
  multBinEdges[1] = 0.5;
  Int_t startBin = 1;
  Int_t endBin = 1;
  for(Int_t multStep = 0; multStep < nMultSteps; multStep++){
    endBin += multSteps[multStep];
    for(Int_t multBin = startBin; multBin < endBin; multBin++)  {
      multBinEdges[multBin+1] = multBinEdges[multBin] + multBinWidth[multStep];
    }
    startBin = endBin;
  }


  AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();

  if (!mgr) {
    Error("AddTask_mkrueger_MeanPt", "No analysis manager found.");
    return 0;
  }

  // Switch off all AliInfo (too much output!!!)
  AliLog::SetGlobalLogLevel(AliLog::kError);
  mgr->SetDebugLevel(0);

  TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  Bool_t hasMC = (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler() != 0x0);

  if(controlstring.Contains("PbPb") && hasMC) offlineTriggerMask = AliVEvent::kMB;

  AliMeanPtAnalysisTask* mainTask = NULL;

  char taskName[100] = "";


  for(Int_t cutMode = cutModeLow; cutMode < cutModeHigh; cutMode++){
    sprintf(taskName, "mkrueger_%s_eta_%.2f_cutMode_%d", colsys.c_str(), etaCut, cutMode);

    AliMeanPtAnalysisTask* task = new AliMeanPtAnalysisTask(taskName);

    if(controlstring.Contains("fullPt")){
      upperPtCut = 50;
      lowerPtCut = 0.;
//      Double_t pTBinEdgesLarge[69] = {0.,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.2,2.4,2.6,2.8,3.0,3.2,3.4,3.6,3.8,4.0,4.5,5.0,5.5,6.0,6.5,7.0,8.0,9.0,10.0,11.0,12.0,13.0,14.0,15.0,16.0,18.0,20.0,22.0,24.0,26.0,28.0,30.0,32.0,34.0,36.0,40.0,45.0,50.0};
      Double_t pTBinEdgesLarge[54] = {0.,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.2,2.4,2.6,2.8,3.0,3.2,3.4,3.6,3.8,4.0,4.5,5.0,5.5,6.0,6.5,7.0,8.0,9.0,10.0,20.0,30.0,40.0,50.0};
      task->SetBinsPt(53, pTBinEdgesLarge);
    }

    if(controlstring.Contains("logBinsPt")){

      Double_t minPt = 0.1;
      Double_t maxPt = 100.;
      Int_t nbinsPt = 30;

      Double_t logminPt = TMath::Log10(minPt);
      Double_t logmaxPt = TMath::Log10(maxPt);
      Double_t binwidth = (logmaxPt-logminPt)/nbinsPt;
      Double_t *binsPt =  new Double_t[nbinsPt+1];
      binsPt[0] = minPt;
      for (Int_t i = 1; i <= nbinsPt; i++) {
        binsPt[i] = minPt + TMath::Power(10, logminPt + i*binwidth);
      }
    }

    task->SetIncludeCrosscheckHistos(includeCrosscheckHistos);
    task->Set2013pA(is2013pA);

    // general cuts
    task->SelectCollisionCandidates(offlineTriggerMask);
    task->SetTriggerMask(offlineTriggerMask);

    task->SetUseMC(hasMC);
    if(type.Contains("ESD")) task->SetUseESD();
    else task->SetUseAOD();
    task->SetBinsMult(nBinsMult, multBinEdges);
    task->SetBinsCent(nBinsCent, centBinEdges);

    // nominal cut-setting:
    task->SetMinEta(-etaCut);
    task->SetMaxEta(etaCut);
    task->SetMinPt(lowerPtCut);
    task->SetMaxPt(upperPtCut);

    task->Set2013pA(kFALSE);
    task->Set2015data(is2015Data);

    task->SetMeanXYZv(0.0,0.0,0.0);
    task->SetSigmaMeanXYZv(1.0,1.0,10.0);
    task->SetZvtx(10.);

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
//MC3--
      if(cutMode == 109) {task->SetMaxChi2TPCConstrained(25.);}
      if(cutMode == 110) {task->SetMaxChi2TPCConstrained(49.);}

      if(cutMode == 111) {task->SetDCAtoVertexXYPtDep("0.0104+0.0200/pt^1.01");}
      if(cutMode == 112) {task->SetDCAtoVertexXYPtDep("0.0260+0.0500/pt^1.01");}
//MC4--
      if(cutMode == 113) {task->SetDCAtoVertexZ(1.0);}
      if(cutMode == 114) {task->SetDCAtoVertexZ(5.0);}

      if(cutMode == 115) {task->SetClusterReqITS(kFALSE);}
//MC5--
      if(cutMode == 116) {task->SetGeometricalCut(kTRUE,3,120,1.5,0.85,0.7);}
      if(cutMode == 117) {task->SetGeometricalCut(kTRUE,3,140,1.5,0.85,0.7);}

      if(cutMode == 118) {task->SetGeometricalCut(kTRUE,4,130,1.5,0.85,0.7);}
      if(cutMode == 119) {task->SetGeometricalCut(kTRUE,2,130,1.5,0.85,0.7);}

      // event cut varaitions
      if(cutMode == 120) task->SetZvtx(5.);
      if(cutMode == 121) task->SetZvtx(20.);


  // hang task in train

    mgr->AddTask(task);

    AliAnalysisDataContainer* cinput = mgr->GetCommonInputContainer();
    AliAnalysisDataContainer* coutput = mgr->CreateContainer(taskName,
							    TList::Class(),
							    AliAnalysisManager::kOutputContainer,
							    "AnalysisResults.root");



    mgr->ConnectInput(task, 0, cinput);
    mgr->ConnectOutput(task, 1, coutput);


  }


  return mainTask;

}
