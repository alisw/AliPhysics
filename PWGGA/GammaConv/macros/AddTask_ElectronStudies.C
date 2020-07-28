void AddTask_ElectronStudies(
  Int_t     trainConfig                   = 1,
  Int_t     isMC                          = 0,
  Int_t     IsHeavyIon                    = 0,
  TString   photonCutNumberV0Reader       = "00200009327000008250400000",
  TString   periodNameV0Reader            = "",
  Bool_t    useHistograms                 = kFALSE, // if activated, analysis will be performed hist based instead of cut based
  TString   corrTaskSetting = "",
  Int_t     enableExtMatchAndQA           = 0,
  Bool_t    enableTriggerOverlapRej       = kTRUE,
  Int_t     enableTriggerMimicking        = 0,        // enable trigger mimicking
  TString   settingMaxFacPtHard           = "3.",       // maximum factor between hardest jet and ptHard generated
  Bool_t    makeAdditionalHistos          = kFALSE,
  // subwagon config
  TString   additionalTrainConfig         = "0"       // additional counter for trainconfig
  ){

  //
  // ─── SET CONFIG ─────────────────────────────────────────────────────────────────
  //

  // Default
  TString   TaskEventCutnumber                = "00010113";
  TString   TaskClusterCutnumberEMC           = "111110001f022700000";
  TString   TaskConvCutnumber                 = "0dm0000922700000dge0404000";


  Int_t trackMatcherRunningMode = 0; // CaloTrackMatcher running mode

  if (additionalTrainConfig.Atoi() > 0){
    trainConfig = trainConfig + additionalTrainConfig.Atoi();
    cout << "INFO: running additionalTrainConfig '" << additionalTrainConfig.Atoi() << "', train config: '" << trainConfig << "'" << endl;
  }
  
  // pp 8 TeV
  // ────────────────────────────────────────────────────────────────────────────────
  Int_t                       fMinClsTPC = 70;  
  Double_t                    fChi2PerClsTPC = 5;   
  Int_t                       fMinClsITS = 0;  
  Double_t                    fEtaCut = 0.9;  
  Double_t                    fPtCut= 0.5;  
  Double_t                    fYMCCut = 9999;  
  if(trainConfig == 1){  // min bias (cuts from PCMEMC 84 + loose iso)
      TaskEventCutnumber                = "00010113";
      TaskClusterCutnumberEMC           = "111113206f032000000";
      TaskConvCutnumber                 = "0dm00009f9730000dge0404000";
  } else if(trainConfig == 2){  // trigger
      TaskEventCutnumber                = "00052113";
      TaskClusterCutnumberEMC           = "111113206f032000000";
      TaskConvCutnumber                 = "0dm00009f9730000dge0404000";

  } else if(trainConfig == 3){  // trigger
      TaskEventCutnumber                = "00081113";
      TaskClusterCutnumberEMC           = "111113206f032000000";
      TaskConvCutnumber                 = "0dm00009f9730000dge0404000";
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

  AliAnalysisTaskElectronStudies *fQA = new AliAnalysisTaskElectronStudies("ElectronStudies");
  
  fQA->SetV0ReaderName(V0ReaderName);
  fQA->SetIsMC(isMC);
  
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

  TString TrackMatcherNameSignal = Form("CaloTrackMatcher_Signal_%s_%i",TaskClusterCutnumberEMC.Data(),trackMatcherRunningMode);


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
  analysisEventCuts->SetTriggerMimicking(enableTriggerMimicking);
  analysisEventCuts->SetFillCutHistograms("",kFALSE);

  // EMC signal cluster cuts (used to store in tree)
  AliCaloPhotonCuts *analysisClusterCutsEMC = new AliCaloPhotonCuts(isMC,"analysisClusterCutsEMC","analysisClusterCutsEMC");
  analysisClusterCutsEMC->SetV0ReaderName(V0ReaderName);
  analysisClusterCutsEMC->SetCorrectionTaskSetting(corrTaskSetting);
  analysisClusterCutsEMC->SetCaloTrackMatcherName(TrackMatcherNameSignal);
  analysisClusterCutsEMC->SetExtendedMatchAndQA(enableExtMatchAndQA);
  analysisClusterCutsEMC->InitializeCutsFromCutString(TaskClusterCutnumberEMC.Data());
  analysisClusterCutsEMC->SetFillCutHistograms("");
  
  AliConversionPhotonCuts* analysisConvCuts = new AliConversionPhotonCuts();
  analysisConvCuts->SetV0ReaderName(V0ReaderName);
  analysisConvCuts->InitializeCutsFromCutString(TaskConvCutnumber.Data());
  analysisConvCuts->SetFillCutHistograms("");
  
  fQA->SetEventCuts(analysisEventCuts,IsHeavyIon);
  fQA->SetClusterCutsEMC(analysisClusterCutsEMC,IsHeavyIon);
  fQA->SetConvCuts(analysisConvCuts,IsHeavyIon);
  fQA->SetCorrectionTaskSetting(corrTaskSetting);
  fQA->SetTrackMatcherRunningMode(trackMatcherRunningMode);  
  fQA->SetMinClsTPC(fMinClsTPC);
  fQA->SetMinClsITS(fMinClsITS);
  fQA->SetChi2PerClsTPC(fChi2PerClsTPC);
  fQA->SetEtaCut(fEtaCut);
  fQA->SetMinPtCut(fPtCut);

  mgr->AddTask(fQA);

  mgr->ConnectInput(fQA, 0,  cinput );
  AliAnalysisDataContainer *coutput = mgr->CreateContainer( Form("ElectronStudies_%d",trainConfig),
                                                            TTree::Class(),
                                                            AliAnalysisManager::kOutputContainer,
                                                            Form("ElectronStudies_%d.root",trainConfig));
  AliAnalysisDataContainer *histos= mgr->CreateContainer( Form("ElectronStudies_histos_%d",trainConfig),
                                                            TList::Class(),
                                                            AliAnalysisManager::kOutputContainer,
                                                            Form("ElectronStudies_histos_%d.root",trainConfig));
  mgr->ConnectOutput (fQA, 1, histos );
  mgr->ConnectOutput (fQA, 2, coutput );

  // mgr->ConnectOutput (fQA, 2, "ElectronStudiesHistos", TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:TreeAnalysisHistos", mgr->GetCommonFileName())) );
  return;
}
