void AddTask_ClusterQA(
  TString   photonCutNumberV0Reader       = "00200009327000008250400000",
  TString   TaskEventCutnumber            = "00010113",
  TString   TaskClusterCutnumber          = "1111100010022700000",
  // TString   TaskMesonCutnumber            = "0163300000000000",
  Int_t     minNLM                        = 1,
  Int_t     maxNLM                        = 1,
  Int_t     isMC                          = 0,
  Int_t     IsHeavyIon                    = 0,
  Bool_t    kHistograms                   = kTRUE,
  Double_t  kTree                         = 1.0,  // 0. / 0 / kFALSE for no, 1. / 1 / kTRUE for yes,  x > 1.0 will use only 1/x of the event statistics for the tree
  TString   V0ReaderCutNumberAODBranch    = "0000000060084001001500000",
  Bool_t    doEtaShiftV0Reader            = kFALSE,
  Bool_t    enableV0findingEffi           = kFALSE,
  TString   periodNameV0Reader            = "",
  TString   corrTaskSetting = "",
  Int_t     enableExtMatchAndQA           = 0,
  Bool_t    doSaveSurroundingTracks       = 1,
  Bool_t    doSaveSurroundingCells        = 1,
  Float_t   minClusterEnergy              = 2.0,
  Float_t   maxConeRadius                 = 0.6,
  Float_t   minTrackMomentum              = 0.3,
  Bool_t    doSaveClusterCells            = 1,
  Bool_t    doSaveEventProp               = 1,
  Bool_t    doSaveEventwiseClusters       = 0,
  Bool_t    enableTriggerOverlapRej       = kTRUE,
  TString   settingMaxFacPtHard           = "3.",       // maximum factor between hardest jet and ptHard generated
  Bool_t    makeAdditionalHistos          = kFALSE
  ){




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

  // TString caloCutPosEMC = TaskClusterCutnumberEMC;
  // caloCutPosEMC.Resize(1);
  // TString TrackMatcherNameEMC = Form("CaloTrackMatcher_%s",caloCutPosEMC.Data());
  // if( !(AliCaloTrackMatcher*)mgr->GetTask(TrackMatcherNameEMC.Data()) ){
  //   AliCaloTrackMatcher* fTrackMatcherEMC = new AliCaloTrackMatcher(TrackMatcherNameEMC.Data(),caloCutPosEMC.Atoi());
  //   fTrackMatcherEMC->SetV0ReaderName(V0ReaderName);
  //   //fTrackMatcherEMC->SetCorrectionTaskSetting(corrTaskSetting);
  //   mgr->AddTask(fTrackMatcherEMC);
  //   mgr->ConnectInput(fTrackMatcherEMC,0,cinput);
  // }
  // TString caloCutPosDMC = TaskClusterCutnumberDMC;
  // caloCutPosDMC.Resize(1);
  // TString TrackMatcherNameDMC = Form("CaloTrackMatcher_%s",caloCutPosDMC.Data());
  // if( !(AliCaloTrackMatcher*)mgr->GetTask(TrackMatcherNameDMC.Data()) ){
  //   AliCaloTrackMatcher* fTrackMatcherDMC = new AliCaloTrackMatcher(TrackMatcherNameDMC.Data(),caloCutPosDMC.Atoi());
  //   fTrackMatcherDMC->SetV0ReaderName(V0ReaderName);
  //   //fTrackMatcherDMC->SetCorrectionTaskSetting(corrTaskSetting);
  //   mgr->AddTask(fTrackMatcherDMC);
  //   mgr->ConnectInput(fTrackMatcherDMC,0,cinput);
  // }

  AliCaloPhotonCuts *analysisClusterCutsEMC = new AliCaloPhotonCuts();
  analysisClusterCutsEMC->SetV0ReaderName(V0ReaderName);
  analysisClusterCutsEMC->SetCorrectionTaskSetting(corrTaskSetting);
  // analysisClusterCutsEMC->SetCaloTrackMatcherName(TrackMatcherNameEMC);
  analysisClusterCutsEMC->SetExtendedMatchAndQA(enableExtMatchAndQA);
  analysisClusterCutsEMC->InitializeCutsFromCutString(TaskClusterCutnumber.Data());
  analysisClusterCutsEMC->SetFillCutHistograms("");

  // AliCaloPhotonCuts *analysisClusterCutsDMC = new AliCaloPhotonCuts();
  // analysisClusterCutsDMC->SetV0ReaderName(V0ReaderName);
  // //analysisClusterCutsDMC->SetCorrectionTaskSetting(corrTaskSetting);
  // // analysisClusterCutsDMC->SetCaloTrackMatcherName(TrackMatcherNameEMC);
  // analysisClusterCutsDMC->SetExtendedMatchAndQA(enableExtMatchAndQA);
  // analysisClusterCutsDMC->InitializeCutsFromCutString(TaskClusterCutnumberDMC.Data());
  // analysisClusterCutsDMC->SetFillCutHistograms("");

  // AliConversionMesonCuts *analysisMesonCuts    = new AliConversionMesonCuts();
  // analysisMesonCuts->SetEnableOpeningAngleCut(kFALSE);
  // analysisMesonCuts->SetIsMergedClusterCut(1);
  // analysisMesonCuts->InitializeCutsFromCutString(TaskMesonCutnumber.Data());
  // analysisMesonCuts->SetFillCutHistograms("");


  AliAnalysisTaskClusterQA *fQA = new AliAnalysisTaskClusterQA(Form("%s_%s_QA",TaskEventCutnumber.Data(),TaskClusterCutnumber.Data()));
  fQA->SetEventCuts(analysisEventCuts,IsHeavyIon);
  fQA->SetClusterCutsEMC(analysisClusterCutsEMC,IsHeavyIon);
  // fQA->SetClusterCutsDMC(analysisClusterCutsDMC,IsHeavyIon);
  // fQA->SetMesonCuts(analysisMesonCuts,IsHeavyIon);
  fQA->SetCorrectionTaskSetting(corrTaskSetting);
  fQA->SetMinMaxNLMCut(minNLM,maxNLM);
  fQA->FillType(kTree,kHistograms);
  fQA->SetIsMC(isMC);
  if(isMC)
    fQA->SetSaveMCInformation(kTRUE);
  fQA->SetSaveSurroundingCells(doSaveSurroundingCells);
  fQA->SetSaveClusterCells(doSaveClusterCells);
  fQA->SetSaveSurroundingTracks(doSaveSurroundingTracks);
  fQA->SetMaxConeRadius(maxConeRadius);
  fQA->SetMinTrackPt(minTrackMomentum);
  fQA->SetMinClusterEnergy(minClusterEnergy);
  fQA->SetSaveEventProperties(doSaveEventProp);
  fQA->SetV0ReaderName(V0ReaderName);
  fQA->SetDoAdditionalHistos(makeAdditionalHistos);
  fQA->SetEventwiseClusterOutput(doSaveEventProp);
  mgr->AddTask(fQA);

  mgr->ConnectInput  (fQA, 0,  cinput );
  mgr->ConnectOutput (fQA, 1, mgr->CreateContainer(!(corrTaskSetting.CompareTo("")) ?  Form("GammaCaloQA_%s_%s", TaskEventCutnumber.Data(), TaskClusterCutnumber.Data()) : Form("GammaCaloQA_%s_%s_%s", TaskEventCutnumber.Data(), TaskClusterCutnumber.Data(),corrTaskSetting.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:GammaCaloQA_%s_%s", mgr->GetCommonFileName(), TaskEventCutnumber.Data(), TaskClusterCutnumber.Data())) );
  mgr->ConnectOutput (fQA, 2, mgr->CreateContainer(!(corrTaskSetting.CompareTo("")) ?  Form("ClusterQA_%s_%s", TaskEventCutnumber.Data(), TaskClusterCutnumber.Data()) : Form("ClusterQA_%s_%s_%s", TaskEventCutnumber.Data(), TaskClusterCutnumber.Data(),corrTaskSetting.Data()), TTree::Class(), AliAnalysisManager::kOutputContainer, mgr->GetCommonFileName()) );

  return;
}
