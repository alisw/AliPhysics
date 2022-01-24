void AddTask_ElectronStudies(
  Int_t     trainConfig                   = 1,
  Int_t     isMC                          = 0,
  Int_t     IsHeavyIon                    = 0,
  TString   photonCutNumberV0Reader       = "00200009327000008250400000",
  TString   periodNameV0Reader            = "",
  Bool_t    useHistograms                 = kFALSE, // if activated, analysis will be performed hist based instead of cut based
  Int_t     enableExtMatchAndQA           = 0,
  Bool_t    enableTriggerOverlapRej       = kTRUE,
  Int_t     enableTriggerMimicking        = 0,        // enable trigger mimicking
  TString   settingMaxFacPtHard           = "3.",       // maximum factor between hardest jet and ptHard generated
  Bool_t    makeAdditionalHistos          = kFALSE,
  TString   fileNameExternalInputs        = "",
  // subwagon config
  TString   additionalTrainConfig         = "0"       // additional counter for trainconfig
  ){

  //
  // ─── SET CONFIG ─────────────────────────────────────────────────────────────────
  //
  AliCutHandlerPCM cuts(13); // only for tokenize
  TString addTaskName = "AddTask_ElectronStudies";
  // Default
  TString   TaskEventCutnumber                = "00010113";
  TString   TaskClusterCutnumberEMC           = "1111100010022700000";
  TString   TaskTMCut                         = "111110000f000000000";
  TString   TaskConvCutnumber                 = "0dm0000922700000dge0404000";



  Int_t trackMatcherRunningMode = 7; // CaloTrackMatcher running mode
  TString sAdditionalTrainConfig      = cuts.GetSpecialSettingFromAddConfig(additionalTrainConfig, "", "", addTaskName);
  if (sAdditionalTrainConfig.Atoi() > 0){
    trainConfig = trainConfig + sAdditionalTrainConfig.Atoi();
    cout << "INFO: running additionalTrainConfig '" << sAdditionalTrainConfig.Atoi() << "', train config: '" << trainConfig << "'" << endl;
  }

  TString fileNamePtWeights           = cuts.GetSpecialFileNameFromString (fileNameExternalInputs, "FPTW:");
  TString fileNameMultWeights         = cuts.GetSpecialFileNameFromString (fileNameExternalInputs, "FMUW:");
  TString fileNameCustomTriggerMimicOADB   = cuts.GetSpecialFileNameFromString (fileNameExternalInputs, "FTRM:");

  TString corrTaskSetting             = cuts.GetSpecialSettingFromAddConfig(additionalTrainConfig, "CF", "", addTaskName);
  if(corrTaskSetting.CompareTo(""))
    cout << "corrTaskSetting: " << corrTaskSetting.Data() << endl;
  
  // pp 8 TeV
  // ────────────────────────────────────────────────────────────────────────────────
  Int_t                       fMinClsTPC = 70;  
  Double_t                    fChi2PerClsTPC = 5;   
  Int_t                       fMinClsITS = 0;  
  Double_t                    fEtaCut = 0.9;  
  Double_t                    fPtCut= 0.5;  
  Double_t                    fYMCCut = 9999;  
  Double_t                    fMaxDCAxy = 9999;  
  Double_t                    fMaxDCAz = 9999;  
  Double_t                    fMinFracClsTPC = 0;

  Double_t               fEtaMatch[3]={0.,0.,0.};
  Double_t                 fPhiMatch[3]={0.,0.,0.};
  Bool_t  fUseRTrackMatching = kTRUE;
  Double_t fRTrackMatching   = 0.01; 

  Float_t fIsoRadius = 0.2;

  // no M02 cut, no exotics cut
  // add NCells

  if(trainConfig == 1){  // min bias 
      // TM pt dep by default to reduce file size
      TaskEventCutnumber                = "00010113";
      TaskClusterCutnumberEMC           = "4117900060l30000000";
                                         //411792106fe32220000 latest and greates
      TaskConvCutnumber                 = "0dm00009f9730000dge0404000";
      TaskTMCut                         = "4117900066l30000000";
  } else if(trainConfig == 2){  // trigger
      TaskEventCutnumber                = "0008e113";
      TaskClusterCutnumberEMC           = "4117900060l30000000";
      TaskConvCutnumber                 = "0dm00009f9730000dge0404000";
      TaskTMCut                         = "4117900066l30000000"; // only used for track mathing

  } else if(trainConfig == 3){  // trigger
      TaskEventCutnumber                = "0008d113";
      TaskClusterCutnumberEMC           = "4117900060l30000000";
      TaskConvCutnumber                 = "0dm00009f9730000dge0404000";
      TaskTMCut                         = "4117900066l30000000";
  // Same cluster cuts as default
  // but trying to replicate track cuts used for electrons as close as possible
  } else if(trainConfig == 4){  // mb
      // Track cuts
      fMinClsTPC = 80;  
      fChi2PerClsTPC = 4;   
      fMinClsITS = 3;  
      fEtaCut = 0.8;  
      fPtCut= 0.5;  
      fYMCCut = 9999;  
      fMaxDCAxy = 1;  
      fMaxDCAz = 2;  
      fMinFracClsTPC = 0.6;
      // cluster cuts
      TaskEventCutnumber                = "00010113";
      TaskClusterCutnumberEMC           = "4117921060e32000000";
                                         //411792106fe32220000 latest and greates
      TaskConvCutnumber                 = "0dm00009f9730000dge0404000";
      TaskTMCut                         = "4117921062e32000000";
  } else if(trainConfig == 5){  // trigger
      // Track cuts
      fMinClsTPC = 80;  
      fChi2PerClsTPC = 4;   
      fMinClsITS = 3;  
      fEtaCut = 0.8;  
      fPtCut= 0.5;  
      fYMCCut = 9999;  
      fMaxDCAxy = 1;  
      fMaxDCAz = 2;  
      fMinFracClsTPC = 0.6;
      // cluster cuts
      TaskEventCutnumber                = "0008e113";
      TaskClusterCutnumberEMC           = "4117921060e32000000";
      TaskConvCutnumber                 = "0dm00009f9730000dge0404000";
      TaskTMCut                         = "4117921062e32000000"; // only used for track mathing

  } else if(trainConfig == 6){  // trigger
      // Track cuts
      fMinClsTPC = 80;  
      fChi2PerClsTPC = 4;   
      fMinClsITS = 3;  
      fEtaCut = 0.8;  
      fPtCut= 0.5;  
      fYMCCut = 9999;  
      fMaxDCAxy = 1;  
      fMaxDCAz = 2;  
      fMinFracClsTPC = 0.6;
      // cluster cuts
      TaskEventCutnumber                = "0008d113";
      TaskClusterCutnumberEMC           = "4117921060e32000000";
      TaskConvCutnumber                 = "0dm00009f9730000dge0404000";
      TaskTMCut                         = "4117921062e32000000";
  } else if(trainConfig == 100){  // no event cuts (to be used for particle gun)
      TaskEventCutnumber                = "00000000";
      TaskClusterCutnumberEMC           = "4117900060l30000000";
                                         //411792106fe32220000 latest and greates
      TaskConvCutnumber                 = "0dm00009f9730000dge0404000";
      TaskTMCut                         = "4117900062l30000000";

  //
  // ─── OLD NL ─────────────────────────────────────────────────────────────────────
  //

  } else if(trainConfig == 10){  // min bias 
      TaskEventCutnumber                = "00010113";
      TaskClusterCutnumberEMC           = "4117905060e32000000";
                                         //411790506fe32220000 latest and greates
      TaskConvCutnumber                 = "0dm00009f9730000dge0404000";
      TaskTMCut                         = "4117905062e32000000";
  } else if(trainConfig == 11){  // trigger
      TaskEventCutnumber                = "0008e113";
      TaskClusterCutnumberEMC           = "4117905060e32000000";
      TaskConvCutnumber                 = "0dm00009f9730000dge0404000";
      TaskTMCut                         = "4117905062e32000000"; // only used for track mathing

  } else if(trainConfig == 12){  // trigger
      TaskEventCutnumber                = "0008d113";
      TaskClusterCutnumberEMC           = "4117905060e32000000";
      TaskConvCutnumber                 = "0dm00009f9730000dge0404000";
      TaskTMCut                         = "4117905062e32000000";
  //
  // ─── Variation new Nonlin1 ─────────────────────────────────────────────────────────────────────
  //

  } else if(trainConfig == 13){  // min bias 
      TaskEventCutnumber                = "00010113";
      TaskClusterCutnumberEMC           = "4117939060e32000000";
                                         //411793906fe32220000 latest and greates
      TaskConvCutnumber                 = "0dm00009f9730000dge0404000";
      TaskTMCut                         = "4117939062e32000000";
  } else if(trainConfig == 14){  // trigger
      TaskEventCutnumber                = "0008e113";
      TaskClusterCutnumberEMC           = "4117939060e32000000";
      TaskConvCutnumber                 = "0dm00009f9730000dge0404000";
      TaskTMCut                         = "4117939062e32000000"; // only used for track mathing

  } else if(trainConfig == 15){  // trigger
      TaskEventCutnumber                = "0008d113";
      TaskClusterCutnumberEMC           = "4117939060e32000000";
      TaskConvCutnumber                 = "0dm00009f9730000dge0404000";
      TaskTMCut                         = "4117939062e32000000";
  //
  // ─── Variation new Nonlin2 ─────────────────────────────────────────────────────────────────────
  //
  } else if(trainConfig == 16){  // min bias 
      TaskEventCutnumber                = "00010113";
      TaskClusterCutnumberEMC           = "4117938060e32000000";
                                         //411793806fe32220000 latest and greates
      TaskConvCutnumber                 = "0dm00009f9730000dge0404000";
      TaskTMCut                         = "4117938062e32000000";
  } else if(trainConfig == 17){  // trigger
      TaskEventCutnumber                = "0008e113";
      TaskClusterCutnumberEMC           = "4117938060e32000000";
      TaskConvCutnumber                 = "0dm00009f9730000dge0404000";
      TaskTMCut                         = "4117938062e32000000"; // only used for track mathing

  } else if(trainConfig == 18){  // trigger
      TaskEventCutnumber                = "0008d113";
      TaskClusterCutnumberEMC           = "4117938060e32000000";
      TaskConvCutnumber                 = "0dm00009f9730000dge0404000";
      TaskTMCut                         = "4117938062e32000000";
  }
  

  TString clusterTypeString(TaskTMCut(0,1));
  Int_t clusterType = clusterTypeString.Atoi();

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
  Bool_t fJetFinderUsage      = kFALSE;
  Bool_t fUsePtHardFromFile      = kFALSE;
  Bool_t fUseAddOutlierRej      = kFALSE;
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
    } else if(strTempSetting.BeginsWith("USEJETFINDER:")){
      strTempSetting.Replace(0,13,"");
      if(strTempSetting.Atoi()==1){
        cout << "using MC jet finder for outlier removal" << endl;
        fJetFinderUsage        = kTRUE;
      }
    } else if(strTempSetting.BeginsWith("PTHFROMFILE:")){
      strTempSetting.Replace(0,12,"");
      if(strTempSetting.Atoi()==1){
        cout << "using MC jet finder for outlier removal" << endl;
        fUsePtHardFromFile        = kTRUE;
      }
    } else if(strTempSetting.BeginsWith("ADDOUTLIERREJ:")){
      strTempSetting.Replace(0,14,"");
      if(strTempSetting.Atoi()==1){
        cout << "using path based outlier removal" << endl;
        fUseAddOutlierRej        = kTRUE;
      }
    } else if(rmaxFacPtHardSetting->GetEntries()==1 && strTempSetting.Atof()>0){
      maxFacPtHard               = strTempSetting.Atof();
      cout << "running with max pT hard jet fraction of: " << maxFacPtHard << endl;
      fMaxPtHardSet        = kTRUE;
    }
  }

  TString TrackMatcherNameSignal = Form("CaloTrackMatcher_Signal_Elec_%s_%i",TaskTMCut.Data(),trackMatcherRunningMode);


  // matching for signal clusters
  if(corrTaskSetting.CompareTo("")){
      TrackMatcherNameSignal = TrackMatcherNameSignal+"_"+corrTaskSetting.Data();
      cout << "Using separate track matcher for correction framework setting: " << TrackMatcherNameSignal.Data() << endl;
  }
  if( !(AliCaloTrackMatcher*)mgr->GetTask(TrackMatcherNameSignal.Data()) ){
      AliCaloTrackMatcher* fTrackMatcherSignal = new AliCaloTrackMatcher(TrackMatcherNameSignal.Data(),clusterType,trackMatcherRunningMode);
      fTrackMatcherSignal->SetV0ReaderName(V0ReaderName);
      fTrackMatcherSignal->SetCorrectionTaskSetting(corrTaskSetting);
      if(TrackMatcherNameSignal.Contains("Elec")){
        fTrackMatcherSignal->SetMassHypothesis(0.000510999);
      }
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
  if(fJetFinderUsage)
      analysisEventCuts->SetUseJetFinderForOutliers(kTRUE);
  if(fUsePtHardFromFile)
    analysisEventCuts->SetUsePtHardBinFromFile(kTRUE);
  if(fUseAddOutlierRej)
    analysisEventCuts->SetUseAdditionalOutlierRejection(kTRUE);
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

  // EMC signal cluster cuts (used to store in tree)
  AliCaloPhotonCuts *analysisClusterCutsEMCTrackMatching = new AliCaloPhotonCuts(isMC,"analysisClusterCutsEMCTrackMatching","analysisClusterCutsEMCTrackMatching");
  analysisClusterCutsEMCTrackMatching->SetV0ReaderName(V0ReaderName);
  analysisClusterCutsEMCTrackMatching->SetCorrectionTaskSetting(corrTaskSetting);
  analysisClusterCutsEMCTrackMatching->SetCaloTrackMatcherName(TrackMatcherNameSignal);
  analysisClusterCutsEMCTrackMatching->SetExtendedMatchAndQA(enableExtMatchAndQA);
  analysisClusterCutsEMCTrackMatching->InitializeCutsFromCutString(TaskTMCut.Data());
  analysisClusterCutsEMCTrackMatching->SetFillCutHistograms("");
  
  AliConversionPhotonCuts* analysisConvCuts = new AliConversionPhotonCuts();
  analysisConvCuts->SetV0ReaderName(V0ReaderName);
  analysisConvCuts->InitializeCutsFromCutString(TaskConvCutnumber.Data());
  analysisConvCuts->SetFillCutHistograms("");
  
  fQA->SetEventCuts(analysisEventCuts,IsHeavyIon);
  fQA->SetClusterCutsEMC(analysisClusterCutsEMC,IsHeavyIon);
  fQA->SetTMCuts(analysisClusterCutsEMCTrackMatching,IsHeavyIon);
  fQA->SetConvCuts(analysisConvCuts,IsHeavyIon);
  fQA->SetCorrectionTaskSetting(corrTaskSetting);
  fQA->SetTrackMatcherRunningMode(trackMatcherRunningMode);  
  fQA->SetMinClsTPC(fMinClsTPC);
  fQA->SetMinClsITS(fMinClsITS);
  fQA->SetChi2PerClsTPC(fChi2PerClsTPC);
  fQA->SetEtaCut(fEtaCut);
  fQA->SetMinPtCut(fPtCut);
  fQA->SetTrackMatcherName(TrackMatcherNameSignal);
  fQA->SetMaxDCA(fMaxDCAxy,fMaxDCAz);
  fQA->SetMinFracClsTPC(fMinFracClsTPC);
  fQA->SetMaxIsoRadius(fIsoRadius);

  fQA->SetEtaMatching(fEtaMatch[0],fEtaMatch[1],fEtaMatch[2]);
  fQA->SetPhiMatching(fPhiMatch[0],fPhiMatch[1],fPhiMatch[2]);
  fQA->SetUseRTrackMatching(fUseRTrackMatching);
  fQA->SetRTrackMatching(fRTrackMatching);


  mgr->AddTask(fQA);

  mgr->ConnectInput(fQA, 0,  cinput );
    AliAnalysisDataContainer *coutput = NULL;
  AliAnalysisDataContainer *histos = NULL;

  if(corrTaskSetting.CompareTo("")){
    coutput =mgr->CreateContainer( Form("ElectronStudies_%d_%s",trainConfig,corrTaskSetting.Data()),
                                                              TTree::Class(),
                                                              AliAnalysisManager::kOutputContainer,
                                                              Form("ElectronStudies_%d.root",trainConfig));
    histos = mgr->CreateContainer( Form("ElectronStudies_histos_%d_%s",trainConfig,corrTaskSetting.Data()),
                                                              TList::Class(),
                                                              AliAnalysisManager::kOutputContainer,
                                                              Form("ElectronStudies_histos_%d.root",trainConfig));
  } else{
    coutput =mgr->CreateContainer( Form("ElectronStudies_%d",trainConfig),
                                                              TTree::Class(),
                                                              AliAnalysisManager::kOutputContainer,
                                                              Form("ElectronStudies_%d.root",trainConfig));
    histos = mgr->CreateContainer( Form("ElectronStudies_histos_%d",trainConfig),
                                                              TList::Class(),
                                                              AliAnalysisManager::kOutputContainer,
                                                              Form("ElectronStudies_histos_%d.root",trainConfig));
   
  }
  mgr->ConnectOutput (fQA, 1, histos );
  mgr->ConnectOutput (fQA, 2, coutput );

  // mgr->ConnectOutput (fQA, 2, "ElectronStudiesHistos", TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:TreeAnalysisHistos", mgr->GetCommonFileName())) );
  return;
}
