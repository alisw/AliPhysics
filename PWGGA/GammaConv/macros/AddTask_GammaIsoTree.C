void AddTask_GammaIsoTree(
  Int_t     trainConfig                   = 1,
  Int_t     isMC                          = 0,
  Int_t     IsHeavyIon                    = 0,
  TString   photonCutNumberV0Reader      = "00200009327000008250400000",
  TString   periodNameV0Reader            = "",
  Bool_t    kHistograms                   = kTRUE,
  TString   corrTaskSetting = "",
  Int_t     enableExtMatchAndQA           = 0,
  Bool_t    enableTriggerOverlapRej       = kTRUE,
  TString   settingMaxFacPtHard           = "3.",       // maximum factor between hardest jet and ptHard generated
  Bool_t    makeAdditionalHistos          = kFALSE,
  Bool_t    storeTracks                   = kTRUE,
  Bool_t    storeEMCalCluster             = kTRUE,
  Bool_t    storePHOSCluster              = kTRUE,
  Bool_t    storeConversions              = kTRUE,
  Bool_t    doIsolation                   = kTRUE,
  Bool_t    doOwnTrackMatching            = kFALSE
  ){

  //
  // ─── SET CONFIG ─────────────────────────────────────────────────────────────────
  //

  // Default
  TString   TaskEventCutnumber                = "00010113";
  TString   TaskClusterCutnumberEMC           = "111110001f022700000";
  TString   TaskClusterCutnumberBackgroundEMC = "111110001f022700000";
  TString   TaskClusterCutnumberPHOS          = "2444411044013300000";
  TString   TaskConvCutnumber                 = "0dm0000922700000dge0404000";

  Float_t trackIsoR[2] = {0.2,0.4};
  Float_t neutralIsoR[2] = {0.2,0.4};
  Int_t trackMatcherRunningMode = 0; // CaloTrackMatcher running mode
  Bool_t backgroundTrackMatching = kTRUE;
  Bool_t doNeutralIso            = kTRUE;
  Bool_t doChargedIso            = kTRUE;
  Bool_t doCellIso               = kTRUE;
  Bool_t doTagging               = kTRUE;
  
  // pp 8 TeV
  // ────────────────────────────────────────────────────────────────────────────────
  if(trainConfig == 1){ 
      TaskEventCutnumber                = "00010113";
      TaskClusterCutnumberEMC           = "1111132060032230000";
      TaskClusterCutnumberBackgroundEMC = "1111100060022700000";
      TaskClusterCutnumberPHOS          = "2444411044013300000";
      TaskConvCutnumber                 = "0dm00009f9730000dge0404000";
  } else if(trainConfig == 2){ 
      TaskEventCutnumber                = "00052113";
      TaskClusterCutnumberEMC           = "1111132060032230000";
      TaskClusterCutnumberBackgroundEMC = "1111100060022700000";
      TaskClusterCutnumberPHOS          = "2444411044013300000";
      TaskConvCutnumber                 = "0dm00009f9730000dge0404000";
  } else if(trainConfig == 3){ 
      TaskEventCutnumber                = "00081113";
      TaskClusterCutnumberEMC           = "1111132010032230000";
      TaskClusterCutnumberBackgroundEMC = "1111100060022700000";
      TaskClusterCutnumberPHOS          = "2444411044013300000";
      TaskConvCutnumber                 = "0dm00009f9730000dge0404000";

  } else if(trainConfig == 4){  // min bias loose cluster cuts
      TaskEventCutnumber                = "00010113";
      TaskClusterCutnumberEMC           = "111113200f002000000";
      TaskClusterCutnumberBackgroundEMC = "111113206f022700000";
      TaskClusterCutnumberPHOS          = "2444411044013300000";
      TaskConvCutnumber                 = "0dm00009f9730000dge0404000";

      backgroundTrackMatching = kTRUE;
      doNeutralIso = kTRUE;
      doChargedIso = kTRUE;
      doTagging = kTRUE;
      doCellIso = kTRUE;

  // pPb 8 TeV
  // ────────────────────────────────────────────────────────────────────────────────
  } else if(trainConfig == 10){
      TaskEventCutnumber                     = "80010123";
      TaskClusterCutnumberEMC                = "1111111060032230000";
      TaskClusterCutnumberBackgroundEMC      = "1111111060032230000";
      TaskClusterCutnumberPHOS               = "2444411044013300000";
      TaskConvCutnumber                      = "0dm00009f9730000dge0404000";
  }

  

  // ================== GetAnalysisManager ===============================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    return ;
  }

  // ================== GetInputEventHandler =============================
  AliVEventHandler *inputHandler=mgr->GetInputEventHandler();

  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();

//=========  Set Cutnumber for V0Reader ================================
  TString cutnumberPhoton     = photonCutNumberV0Reader.Data();
  TString cutnumberEvent      = "00000003";
  if(IsHeavyIon==1)
    cutnumberEvent = "10000003";
  else if(IsHeavyIon==2)
    cutnumberEvent = "80000003";

//========= Check V0 Reader in  ANALYSIS manager  =====
  TString V0ReaderName        = Form("V0ReaderV1_%s_%s",cutnumberEvent.Data(),cutnumberPhoton.Data());
  AliV0ReaderV1 *fV0ReaderV1  =  NULL;
  if( !(AliV0ReaderV1*)mgr->GetTask(V0ReaderName.Data()) ){
    cout << "V0Reader: " << V0ReaderName.Data() << " not found!!"<< endl;
    return;
  } else {
    cout << "V0Reader: " << V0ReaderName.Data() << " found!!"<< endl;
  }

  AliAnalysisTaskGammaIsoTree *fQA = new AliAnalysisTaskGammaIsoTree("GammaIsoTree");
  
  fQA->SetV0ReaderName(V0ReaderName);
  fQA->SetIsMC(isMC);
  fQA->SetYCutMC(0.9);
  
  // fQA->SetSaveClusterCells(doSaveClusterCells);
  // fQA->SetSaveEventProperties(doSaveEventProp);
  // fQA->SetDoAdditionalHistos(makeAdditionalHistos);

  TObjArray *rmaxFacPtHardSetting = settingMaxFacPtHard.Tokenize("_");
  if(rmaxFacPtHardSetting->GetEntries()<1){cout << "ERROR: AddTask_ClusterQA during parsing of settingMaxFacPtHard String '" << settingMaxFacPtHard.Data() << "'" << endl; return;}
  Bool_t fMinPtHardSet        = kFALSE;
  Double_t minFacPtHard       = -1;
  Bool_t fMaxPtHardSet        = kFALSE;
  Double_t maxFacPtHard       = 100;
  Bool_t fSingleMaxPtHardSet  = kFALSE;
  Double_t maxFacPtHardSingle = 100;
  for(Int_t i = 0; i<rmaxFacPtHardSetting->GetEntries() ; i++){
    TObjString* tempObjStrPtHardSetting     = (TObjString*) rmaxFacPtHardSetting->At(i);
    TString strTempSetting                  = tempObjStrPtHardSetting->GetString();
    if(strTempSetting.BeginsWith("MINPTHFAC:")){
      strTempSetting.Replace(0,10,"");
      minFacPtHard               = strTempSetting.Atof();
      cout << "running with min pT hard jet fraction of: " << minFacPtHard << endl;
      fMinPtHardSet        = kTRUE;
    } else if(strTempSetting.BeginsWith("MAXPTHFAC:")){
      strTempSetting.Replace(0,10,"");
      maxFacPtHard               = strTempSetting.Atof();
      cout << "running with max pT hard jet fraction of: " << maxFacPtHard << endl;
      fMaxPtHardSet        = kTRUE;
    } else if(strTempSetting.BeginsWith("MAXPTHFACSINGLE:")){
      strTempSetting.Replace(0,16,"");
      maxFacPtHardSingle         = strTempSetting.Atof();
      cout << "running with max single particle pT hard fraction of: " << maxFacPtHardSingle << endl;
      fSingleMaxPtHardSet        = kTRUE;
    } else if(rmaxFacPtHardSetting->GetEntries()==1 && strTempSetting.Atof()>0){
      maxFacPtHard               = strTempSetting.Atof();
      cout << "running with max pT hard jet fraction of: " << maxFacPtHard << endl;
      fMaxPtHardSet        = kTRUE;
    }
  }

  // add track matcher if do own trackmatching is enabled

  TString TrackMatcherNameSignal = Form("CaloTrackMatcher_Signal_%s_%i",TaskClusterCutnumberEMC.Data(),trackMatcherRunningMode);
  TString TrackMatcherNameBackground = Form("CaloTrackMatcher_Background_%s_%i",TaskClusterCutnumberBackgroundEMC.Data(),trackMatcherRunningMode);
  
  if(!doOwnTrackMatching){
    
    // matching for signal clusters
    if(corrTaskSetting.CompareTo("")){
      TrackMatcherNameSignal = TrackMatcherNameSignal+"_"+corrTaskSetting.Data();
      cout << "Using separate track matcher for correction framework setting: " << TrackMatcherNameSignal.Data() << endl;
    }
    if( !(AliCaloTrackMatcher*)mgr->GetTask(TrackMatcherNameSignal.Data()) ){
      AliCaloTrackMatcher* fTrackMatcherSignal = new AliCaloTrackMatcher(TrackMatcherNameSignal.Data(),1,trackMatcherRunningMode);
      fTrackMatcherSignal->SetV0ReaderName(V0ReaderName);
      fTrackMatcherSignal->SetCorrectionTaskSetting(corrTaskSetting);
      mgr->AddTask(fTrackMatcherSignal);
      mgr->ConnectInput(fTrackMatcherSignal,0,cinput);
    }

    // matching for background clusters
    if(corrTaskSetting.CompareTo("")){
      TrackMatcherNameBackground = TrackMatcherNameBackground+"_"+corrTaskSetting.Data();
      cout << "Using separate track matcher for correction framework setting: " << TrackMatcherNameBackground.Data() << endl;
    }
    if( !(AliCaloTrackMatcher*)mgr->GetTask(TrackMatcherNameBackground.Data()) ){
      AliCaloTrackMatcher* fTrackMatcherBackground = new AliCaloTrackMatcher(TrackMatcherNameBackground.Data(),1,trackMatcherRunningMode);
      fTrackMatcherBackground->SetV0ReaderName(V0ReaderName);
      fTrackMatcherBackground->SetCorrectionTaskSetting(corrTaskSetting);
      mgr->AddTask(fTrackMatcherBackground);
      mgr->ConnectInput(fTrackMatcherBackground,0,cinput);
    }
  }

  // Create Event Cuts
  AliConvEventCuts *analysisEventCuts = new AliConvEventCuts();
  analysisEventCuts->SetV0ReaderName(V0ReaderName);

  analysisEventCuts->SetTriggerOverlapRejecion(enableTriggerOverlapRej);
  if(fMinPtHardSet)
    analysisEventCuts->SetMinFacPtHard(minFacPtHard);
  if(fMaxPtHardSet)
    analysisEventCuts->SetMaxFacPtHard(maxFacPtHard);
  if(fSingleMaxPtHardSet)
    analysisEventCuts->SetMaxFacPtHardSingleParticle(maxFacPtHardSingle);
  analysisEventCuts->SetCorrectionTaskSetting(corrTaskSetting);
  if (periodNameV0Reader.CompareTo("") != 0) analysisEventCuts->SetPeriodEnum(periodNameV0Reader);
  analysisEventCuts->InitializeCutsFromCutString(TaskEventCutnumber.Data());
  analysisEventCuts->SetFillCutHistograms("",kFALSE);

  // EMC signal cluster cuts (used to store in tree)
  AliCaloPhotonCuts *analysisClusterCutsEMC = new AliCaloPhotonCuts(isMC,"analysisClusterCutsEMC","analysisClusterCutsEMC");
  analysisClusterCutsEMC->SetV0ReaderName(V0ReaderName);
  analysisClusterCutsEMC->SetCorrectionTaskSetting(corrTaskSetting);
  analysisClusterCutsEMC->SetCaloTrackMatcherName(TrackMatcherNameSignal);
  analysisClusterCutsEMC->SetExtendedMatchAndQA(enableExtMatchAndQA);
  analysisClusterCutsEMC->InitializeCutsFromCutString(TaskClusterCutnumberEMC.Data());
  analysisClusterCutsEMC->SetFillCutHistograms("");


  // EMC background cluster cuts (used to calculate iso and tagging)
  AliCaloPhotonCuts *analysisClusterCutsBackgroundEMC = new AliCaloPhotonCuts(isMC,"analysisClusterCutsBackgroundEMC","analysisClusterCutsBackgroundEMC");
  analysisClusterCutsBackgroundEMC->SetV0ReaderName(V0ReaderName);
  analysisClusterCutsBackgroundEMC->SetCorrectionTaskSetting(corrTaskSetting);
  analysisClusterCutsBackgroundEMC->SetCaloTrackMatcherName(TrackMatcherNameBackground);
  analysisClusterCutsBackgroundEMC->SetExtendedMatchAndQA(enableExtMatchAndQA);
  analysisClusterCutsBackgroundEMC->InitializeCutsFromCutString(TaskClusterCutnumberBackgroundEMC.Data());
  analysisClusterCutsBackgroundEMC->SetFillCutHistograms("");

  // PHOS cluster cuts
  AliCaloPhotonCuts *analysisClusterCutsPHOS = new AliCaloPhotonCuts(isMC,"analysisClusterCutsPHOS","analysisClusterCutsPHOS");
  analysisClusterCutsPHOS->SetV0ReaderName(V0ReaderName);
  // analysisClusterCutsPHOS->SetCaloTrackMatcherName(TrackMatcherNamePHOS);
  analysisClusterCutsPHOS->SetExtendedMatchAndQA(enableExtMatchAndQA);
  analysisClusterCutsPHOS->InitializeCutsFromCutString(TaskClusterCutnumberPHOS.Data());
  analysisClusterCutsPHOS->SetFillCutHistograms("");
  
  AliConversionPhotonCuts* analysisConvCuts = new AliConversionPhotonCuts();
  analysisConvCuts->SetV0ReaderName(V0ReaderName);
  analysisConvCuts->InitializeCutsFromCutString(TaskConvCutnumber.Data());
  analysisConvCuts->SetFillCutHistograms("");
  
  fQA->SetEventCuts(analysisEventCuts,IsHeavyIon);
  fQA->SetClusterCutsEMC(analysisClusterCutsEMC,IsHeavyIon);
  fQA->SetClusterCutsBackgroundEMC(analysisClusterCutsBackgroundEMC,IsHeavyIon);
  fQA->SetClusterCutsPHOS(analysisClusterCutsPHOS,IsHeavyIon);
  fQA->SetConvCuts(analysisConvCuts,IsHeavyIon);
  fQA->SetDoTrackIso(kTRUE);
  fQA->SetDoBackgroundTrackMatching(backgroundTrackMatching);
  fQA->SetDoTrackIso(doChargedIso);
  fQA->SetDoNeutralIso(doNeutralIso);
  fQA->SetDoTagging(doTagging);
  fQA->SetDoCellIso(doCellIso);
  fQA->SetTrackIsoR(trackIsoR[0],trackIsoR[1]);
  fQA->SetNeutralIsoR(trackIsoR[0],trackIsoR[1]);
  fQA->SetCorrectionTaskSetting(corrTaskSetting);
  fQA->SetSaveConversions(storeConversions);
  fQA->SetSaveEMCClusters(storeEMCalCluster);
  fQA->SetSavePHOSClusters(storePHOSCluster);
  fQA->SetSaveTracks(storeTracks);
  fQA->SetBuffSize(60*1024*1024);
  fQA->SetTrackMatcherRunningMode(trackMatcherRunningMode);
  fQA->SetDoOwnTrackMatching(doOwnTrackMatching);
  
  mgr->AddTask(fQA);

  mgr->ConnectInput(fQA, 0,  cinput );
  AliAnalysisDataContainer *coutput = mgr->CreateContainer( Form("GammaIsoTree_%d",trainConfig),
                                                            TTree::Class(),
                                                            AliAnalysisManager::kOutputContainer,
                                                            Form("GammaIsoTree_%d.root",trainConfig));
  AliAnalysisDataContainer *histos= mgr->CreateContainer( Form("GammaIsoTree_histos_%d",trainConfig),
                                                            TList::Class(),
                                                            AliAnalysisManager::kOutputContainer,
                                                            Form("GammaIsoTree_histos_%d.root",trainConfig));
  mgr->ConnectOutput (fQA, 1, histos );
  mgr->ConnectOutput (fQA, 2, coutput );

  // mgr->ConnectOutput (fQA, 2, "GammaIsoTreeHistos", TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:TreeAnalysisHistos", mgr->GetCommonFileName())) );
  return;
}
