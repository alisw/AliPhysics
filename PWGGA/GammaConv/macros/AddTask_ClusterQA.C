void AddTask_ClusterQA(
  TString   photonCutNumberV0Reader       = "00200009327000008250400000",
  TString   TaskEventCutnumber            = "00010113",
  TString   TaskClusterCutnumberEMC       = "1111100010022700000",
  TString   TaskClusterCutnumberDMC       = "3885500010022700000",
  TString   TaskMesonCutnumber            = "0163300000000000",
  Int_t     minNLM                        = 1,
  Int_t     maxNLM                        = 1,
  Bool_t    isMC                          = kFALSE,
  Int_t     IsHeavyIon                    = 0,
  Bool_t    kHistograms                   = kTRUE,
  Double_t  kTree                         = 1.0,  // 0. / 0 / kFALSE for no, 1. / 1 / kTRUE for yes,  x > 1.0 will use only 1/x of the event statistics for the tree
  TString   V0ReaderCutNumberAODBranch    = "0000000060084001001500000",
  Bool_t    doEtaShiftV0Reader            = kFALSE,
  Bool_t    enableV0findingEffi           = kFALSE,
  TString   periodNameV0Reader            = "",
  TString   corrTaskSetting = "",
  Int_t     enableExtMatchAndQA           = 5,
  Bool_t    doSaveSurroundingTracks       = 1,
  Bool_t    doSaveSurroundingCells        = 1,
  Float_t   minClusterEnergy              = 2.0,
  Float_t   maxConeRadius                 = 0.6,
  Float_t   minTrackMomentum              = 0.3,
  Bool_t    doSaveClusterCells            = 1,
  Bool_t    doSaveEventProp               = 1,
  Bool_t    enableTriggerOverlapRej       = kTRUE,
  Float_t   maxFacPtHard                  = 3.
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

  AliConvEventCuts *analysisEventCuts = new AliConvEventCuts();
  analysisEventCuts->SetV0ReaderName(V0ReaderName);

  analysisEventCuts->SetTriggerOverlapRejecion(enableTriggerOverlapRej);
  analysisEventCuts->SetMaxFacPtHard(maxFacPtHard);
  analysisEventCuts->SetCorrectionTaskSetting(corrTaskSetting);
  if (periodNameV0Reader.CompareTo("") != 0) analysisEventCuts->SetPeriodEnum(periodNameV0Reader);  
  analysisEventCuts->InitializeCutsFromCutString(TaskEventCutnumber.Data());
  analysisEventCuts->SetFillCutHistograms("",kFALSE);

  TString caloCutPosEMC = TaskClusterCutnumberEMC;
  caloCutPosEMC.Resize(1);
  TString TrackMatcherNameEMC = Form("CaloTrackMatcher_%s",caloCutPosEMC.Data());
  if( !(AliCaloTrackMatcher*)mgr->GetTask(TrackMatcherNameEMC.Data()) ){
    AliCaloTrackMatcher* fTrackMatcherEMC = new AliCaloTrackMatcher(TrackMatcherNameEMC.Data(),caloCutPosEMC.Atoi());
    fTrackMatcherEMC->SetV0ReaderName(V0ReaderName);
    //fTrackMatcherEMC->SetCorrectionTaskSetting(corrTaskSetting);
    mgr->AddTask(fTrackMatcherEMC);
    mgr->ConnectInput(fTrackMatcherEMC,0,cinput);
  }
  TString caloCutPosDMC = TaskClusterCutnumberDMC;
  caloCutPosDMC.Resize(1);
  TString TrackMatcherNameDMC = Form("CaloTrackMatcher_%s",caloCutPosDMC.Data());
  if( !(AliCaloTrackMatcher*)mgr->GetTask(TrackMatcherNameDMC.Data()) ){
    AliCaloTrackMatcher* fTrackMatcherDMC = new AliCaloTrackMatcher(TrackMatcherNameDMC.Data(),caloCutPosDMC.Atoi());
    fTrackMatcherDMC->SetV0ReaderName(V0ReaderName);
    //fTrackMatcherDMC->SetCorrectionTaskSetting(corrTaskSetting);
    mgr->AddTask(fTrackMatcherDMC);
    mgr->ConnectInput(fTrackMatcherDMC,0,cinput);
  }
  
  AliCaloPhotonCuts *analysisClusterCutsEMC = new AliCaloPhotonCuts();
  analysisClusterCutsEMC->SetV0ReaderName(V0ReaderName);
  //analysisClusterCutsEMC->SetCorrectionTaskSetting(corrTaskSetting);
  analysisClusterCutsEMC->SetCaloTrackMatcherName(TrackMatcherNameEMC);
  analysisClusterCutsEMC->SetExtendedMatchAndQA(enableExtMatchAndQA);
  analysisClusterCutsEMC->InitializeCutsFromCutString(TaskClusterCutnumberEMC.Data());
  analysisClusterCutsEMC->SetFillCutHistograms("");
  
  AliCaloPhotonCuts *analysisClusterCutsDMC = new AliCaloPhotonCuts();
  analysisClusterCutsDMC->SetV0ReaderName(V0ReaderName);
  //analysisClusterCutsDMC->SetCorrectionTaskSetting(corrTaskSetting);
  analysisClusterCutsDMC->SetCaloTrackMatcherName(TrackMatcherNameEMC);
  analysisClusterCutsDMC->SetExtendedMatchAndQA(enableExtMatchAndQA);
  analysisClusterCutsDMC->InitializeCutsFromCutString(TaskClusterCutnumberDMC.Data());
  analysisClusterCutsDMC->SetFillCutHistograms("");

  AliConversionMesonCuts *analysisMesonCuts    = new AliConversionMesonCuts();
  analysisMesonCuts->SetEnableOpeningAngleCut(kFALSE);
  analysisMesonCuts->SetIsMergedClusterCut(1);
  analysisMesonCuts->InitializeCutsFromCutString(TaskMesonCutnumber.Data());
  analysisMesonCuts->SetFillCutHistograms("");


  AliAnalysisTaskClusterQA *fQA = new AliAnalysisTaskClusterQA(Form("%s_%s_%s_QA",TaskEventCutnumber.Data(),TaskClusterCutnumberEMC.Data(),TaskClusterCutnumberEMC.Data()));
  fQA->SetEventCuts(analysisEventCuts,IsHeavyIon);
  fQA->SetClusterCutsEMC(analysisClusterCutsEMC,IsHeavyIon);
  fQA->SetClusterCutsDMC(analysisClusterCutsDMC,IsHeavyIon);
  fQA->SetMesonCuts(analysisMesonCuts,IsHeavyIon);
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
  mgr->AddTask(fQA);


  mgr->ConnectInput  (fQA, 0,  cinput );
  mgr->ConnectOutput (fQA, 1, mgr->CreateContainer(Form("GammaCaloQA_%s_%s_%s", TaskEventCutnumber.Data(), TaskClusterCutnumberEMC.Data(),TaskClusterCutnumberEMC.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:GammaCaloQA_%s_%s_%s", mgr->GetCommonFileName(), TaskEventCutnumber.Data(), TaskClusterCutnumberEMC.Data(),TaskClusterCutnumberEMC.Data())) );
  mgr->ConnectOutput (fQA, 2, mgr->CreateContainer("ClusterTree", TTree::Class(), AliAnalysisManager::kOutputContainer, mgr->GetCommonFileName()) );

  return;
}

