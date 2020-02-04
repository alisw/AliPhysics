/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                       *
 * Author: Friederike Bock, Daniel Mühlheim                     *
 * Version 1.0                                 *
 *                                       *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its    *
 * documentation strictly for non-commercial purposes is hereby granted    *
 * without fee, provided that the above copyright notice appears in all    *
 * copies and that both the copyright notice and this permission notice    *
 * appear in the supporting documentation. The authors make no claims    *
 * about the suitability of this software for any purpose. It is      *
 * provided "as is" without express or implied warranty.               *
 **************************************************************************/

//***************************************************************************************
//This AddTask is supposed to set up the main task
//($ALIPHYSICS/PWGGA/GammaConv/AliAnalysisTaskGammaConvCalo.cxx) for
//pPb together with all supporting classes
//***************************************************************************************

//***************************************************************************************
//main function
//***************************************************************************************
void AddTask_GammaConvCalo_pPb(
  Int_t     trainConfig                   = 1,        // change different set of cuts
  Int_t     isMC                          = 0,        // run MC
  TString   photonCutNumberV0Reader       = "",       // 00000008400000000100000000 nom. B, 00000088400000000100000000 low B
  TString   periodNameV0Reader            = "",
  // general setting for task
  Int_t     enableQAMesonTask             = 0,        // enable QA in AliAnalysisTaskGammaConvV1
  Int_t     enableQAPhotonTask            = 0,        // enable additional QA task
  Int_t     enableExtMatchAndQA           = 0,        // disabled (0), extMatch (1), extQA_noCellQA (2), extMatch+extQA_noCellQA (3), extQA+cellQA (4), extMatch+extQA+cellQA (5)
  Int_t     enableLightOutput             = 0,        // switch to run light output (only essential histograms for afterburner)
  Bool_t    enableTHnSparse               = kFALSE,   // switch on THNsparse
  Int_t     enableTriggerMimicking        = 0,        // enable trigger mimicking
  Bool_t    enableTriggerOverlapRej       = kFALSE,   // enable trigger overlap rejection
  TString   settingMaxFacPtHard           = "3.",       // maximum factor between hardest jet and ptHard generated
  Int_t     debugLevel                    = 0,        // introducing debug levels for grid running
  // settings for weights
  // FPTW:fileNamePtWeights, FMUW:fileNameMultWeights, FMAW:fileNameMatBudWeights, separate with ;
 // Material Budget Weights file for Run 2
  // FMAW:alien:///alice/cern.ch/user/a/amarin//MBW/MCInputFileMaterialBudgetWeightsLHC16_Pythia_00010103_0d000009266300008850404000_date181214.root
  TString   fileNameExternalInputs        = "",
  Int_t     doWeightingPart               = 0,        // enable Weighting
  TString   generatorName                 = "DPMJET", // generator Name
  Bool_t    enableMultiplicityWeighting   = kFALSE,   //
  Int_t     enableMatBudWeightsPi0        = 0,        // 1 = three radial bins, 2 = 10 radial bins (2 is the default when using weights)
  Bool_t    enableElecDeDxPostCalibration = kFALSE,
  TString   periodNameAnchor              = "",       //
  // special settings
  Bool_t    enableSortingMCLabels         = kTRUE,    // enable sorting for MC cluster labels
  Bool_t    doPrimaryTrackMatching        = kTRUE,                  // enable basic track matching for all primary tracks to cluster
  // subwagon config
  TString   additionalTrainConfig         = "0"       // additional counter for trainconfig
  ) {

  AliCutHandlerPCM cuts;

  TString fileNamePtWeights           = cuts.GetSpecialFileNameFromString (fileNameExternalInputs, "FPTW:");
  TString fileNameMultWeights         = cuts.GetSpecialFileNameFromString (fileNameExternalInputs, "FMUW:");
  TString fileNameMatBudWeights       = cuts.GetSpecialFileNameFromString (fileNameExternalInputs, "FMAW:");
  TString fileNamedEdxPostCalib       = cuts.GetSpecialFileNameFromString (fileNameExternalInputs, "FEPC:");
  TString fileNameCustomTriggerMimicOADB   = cuts.GetSpecialFileNameFromString (fileNameExternalInputs, "FTRM:");

  TString addTaskName                 = "AddTask_GammaConvCalo_pPb";
  TString sAdditionalTrainConfig      = cuts.GetSpecialSettingFromAddConfig(additionalTrainConfig, "", "", addTaskName);
  if (sAdditionalTrainConfig.Atoi() > 0){
    trainConfig = trainConfig + sAdditionalTrainConfig.Atoi();
    cout << "INFO: " << addTaskName.Data() << " running additionalTrainConfig '" << sAdditionalTrainConfig.Atoi() << "', train config: '" << trainConfig << "'" << endl;
  }
  TString corrTaskSetting             = cuts.GetSpecialSettingFromAddConfig(additionalTrainConfig, "CF", "", addTaskName);
  if(corrTaskSetting.CompareTo(""))
    cout << "corrTaskSetting: " << corrTaskSetting.Data() << endl;

  if(additionalTrainConfig.Contains("MaterialBudgetWeights"))
    fileNameMatBudWeights         = cuts.GetSpecialSettingFromAddConfig(additionalTrainConfig, "MaterialBudgetWeights",fileNameMatBudWeights, addTaskName);


  Int_t trackMatcherRunningMode = 0; // CaloTrackMatcher running mode
  TString strTrackMatcherRunningMode  = cuts.GetSpecialSettingFromAddConfig(additionalTrainConfig, "TM", "", addTaskName);
  if(additionalTrainConfig.Contains("TM"))
    trackMatcherRunningMode = strTrackMatcherRunningMode.Atoi();

  Bool_t doTreeClusterShowerShape = kFALSE; // switch to produce EOverP tree
  TString strdoTreeClusterShowerShape             = cuts.GetSpecialSettingFromAddConfig(additionalTrainConfig, "INVMASSCLUSTree", "", addTaskName);
  if(strdoTreeClusterShowerShape.Atoi()==1)
    doTreeClusterShowerShape = kTRUE;

  TH1S* histoAcc = 0x0;         // histo for modified acceptance
  TString strModifiedAcc              = cuts.GetSpecialSettingFromAddConfig(additionalTrainConfig, "MODIFYACC", "", addTaskName);
  if(strModifiedAcc.Contains("MODIFYACC")){
    TString tempType = strModifiedAcc;
    tempType.Replace(0,9,"");
    cout << "INFO: connecting to alien..." << endl;
    TGrid::Connect("alien://");
    cout << "done!" << endl;
    TFile *w = TFile::Open(fileNamePtWeights.Data());
    if(!w){cout << "ERROR: Could not open file: " << fileNamePtWeights.Data() << endl;return;}
    histoAcc = (TH1S*) w->Get(tempType.Data());
    if(!histoAcc) {cout << "ERROR: Could not find histo: " << tempType.Data() << endl;return;}
    cout << "found: " << histoAcc << endl;
  }

  Int_t localDebugFlag = 0;
  TString strLocalDebugFlag              = cuts.GetSpecialSettingFromAddConfig(additionalTrainConfig, "LOCALDEBUGFLAG", "", addTaskName);
  if(strLocalDebugFlag.Atoi()>0)
    localDebugFlag = strLocalDebugFlag.Atoi();

  TObjArray *rmaxFacPtHardSetting = settingMaxFacPtHard.Tokenize("_");
  if(rmaxFacPtHardSetting->GetEntries()<1){cout << "ERROR: AddTask_GammaConvCalo_pPb during parsing of settingMaxFacPtHard String '" << settingMaxFacPtHard.Data() << "'" << endl; return;}
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

  Int_t isHeavyIon = 2;

  // ================== GetAnalysisManager ===============================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error(Form("%s_%i", addTaskName.Data(),  trainConfig), "No analysis manager found.");
    return ;
  }

  // ================== GetInputEventHandler =============================
  AliVEventHandler *inputHandler=mgr->GetInputEventHandler();

  //=========  Set Cutnumber for V0Reader ================================
  TString cutnumberPhoton = photonCutNumberV0Reader.Data();
  TString cutnumberEvent = "80000003";
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();

    //========= Add V0 Reader to  ANALYSIS manager if not yet existent =====
  TString V0ReaderName        = Form("V0ReaderV1_%s_%s",cutnumberEvent.Data(),cutnumberPhoton.Data());
  AliV0ReaderV1 *fV0ReaderV1  =  NULL;
  if( !(AliV0ReaderV1*)mgr->GetTask(V0ReaderName.Data()) ){
    cout << "V0Reader: " << V0ReaderName.Data() << " not found!!"<< endl;
    return;
  } else {
    cout << "V0Reader: " << V0ReaderName.Data() << " found!!"<< endl;
  }

  //================================================
  //========= Add task to the ANALYSIS manager =====
  //================================================
  AliAnalysisTaskGammaConvCalo *task=NULL;
  task= new AliAnalysisTaskGammaConvCalo(Form("GammaConvCalo_%i",trainConfig));
  task->SetIsHeavyIon(isHeavyIon);
  task->SetIsMC(isMC);
  task->SetV0ReaderName(V0ReaderName);
  task->SetCorrectionTaskSetting(corrTaskSetting);
  if (enableLightOutput > 1) task->SetLightOutput(kTRUE);
  task->SetDoPrimaryTrackMatching(doPrimaryTrackMatching);
  task->SetTrackMatcherRunningMode(trackMatcherRunningMode);
  if(trainConfig >= 1520 && trainConfig < 1530) task->SetDoHBTHistoOutput(kTRUE);

  // cluster cuts
  // 0 "ClusterType",  1 "EtaMin", 2 "EtaMax", 3 "PhiMin", 4 "PhiMax", 5 "DistanceToBadChannel", 6 "Timing", 7 "TrackMatching", 8 "ExoticCell",
  // 9 "MinEnergy", 10 "MinNCells", 11 "MinM02", 12 "MaxM02", 13 "MinMaxM20", 14 "RecConv", 15 "MaximumDispersion", 16 "NLM"

  //************************************************ PCM- EDC analysis 5 TeV pPb *********************************************
  if (trainConfig == 1){ // EMC  INT7 run1 & run2
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","411790105f032230000","0h63103100000010"); // 0-100% without NL
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","111110105f032230000","0h63103100000010"); // 0-100% without NL just EMC
  } else if (trainConfig == 2){ // EMC  INT7 run1 & run2
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","411793105f032230000","0h63103100000010"); // 0-100% PCM NL
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","111113105f032230000","0h63103100000010"); // 0-100% PCM NL just EMC
  } else if (trainConfig == 3){ // EMC EMC triggers
    cuts.AddCutPCMCalo("80052113","00200009f9730000dge0400000","111110105f032230000","0h63103100000010"); // 0-100% without NL just EMC, EMC7
    cuts.AddCutPCMCalo("80052113","00200009f9730000dge0400000","111113105f032230000","0h63103100000010"); // 0-100% PCM NL just EMC, EMC7
  } else if (trainConfig == 4){ // EMC EMC triggers
    cuts.AddCutPCMCalo("80083113","00200009f9730000dge0400000","111110105f032230000","0h63103100000010"); // 0-100% without NL just EMC, EG1
    cuts.AddCutPCMCalo("80083113","00200009f9730000dge0400000","111113105f032230000","0h63103100000010"); // 0-100% PCM NL just EMC, EG1
  } else if (trainConfig == 5){ // EMC EMC triggers
    cuts.AddCutPCMCalo("80085113","00200009f9730000dge0400000","111110105f032230000","0h63103100000010"); // 0-100% without NL just EMC, EG2
    cuts.AddCutPCMCalo("80085113","00200009f9730000dge0400000","111113105f032230000","0h63103100000010"); // 0-100% PCM NL just EMC, EG2
  } else if (trainConfig == 6){ // EMC EMC triggers
    cuts.AddCutPCMCalo("80083123","00200009f9730000dge0400000","111110105f032230000","0h63103100000010"); // 0-100% without NL just EMC, EG1
    cuts.AddCutPCMCalo("80083123","00200009f9730000dge0400000","111113105f032230000","0h63103100000010"); // 0-100% PCM NL just EMC, EG1
  } else if (trainConfig == 7){ // EMC EMC triggers
    cuts.AddCutPCMCalo("80085123","00200009f9730000dge0400000","111110105f032230000","0h63103100000010"); // 0-100% without NL just EMC, EG2
    cuts.AddCutPCMCalo("80085123","00200009f9730000dge0400000","111113105f032230000","0h63103100000010"); // 0-100% PCM NL just EMC, EG2
  } else if (trainConfig == 8){ // EMC  INT7 run1 & run2
    cuts.AddCutPCMCalo("80010123","00200009f9730000dge0400000","411793105f032230000","0h63103100000010"); // 0-100% PCM NL
    cuts.AddCutPCMCalo("80010123","00200009f9730000dge0400000","111113105f032230000","0h63103100000010"); // 0-100% PCM NL just EMC

  //************************************************ PCM- EDC analysis 5 TeV pPb INT7 sys *********************************
  } else if (trainConfig == 10) { // PCM variations
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","411793105f032230000","0h63103100000010"); //New standard cut for 8TeV analysis for RpA
    cuts.AddCutPCMCalo("80010013","00200009f9730000dge0400000","411793105f032230000","0h63103100000010"); // no SPD pileup cut
    cuts.AddCutPCMCalo("80010113","00100009f9730000dge0400000","411793105f032230000","0h63103100000010"); // R cut 2.8 -180 cm
    cuts.AddCutPCMCalo("80010113","00500009f9730000dge0400000","411793105f032230000","0h63103100000010"); // R cut 10. -180 cm
  } else if (trainConfig == 11) {
    cuts.AddCutPCMCalo("80010113","00200069f9730000dge0400000","411793105f032230000","0h63103100000010"); // min pT 40 MeV
    cuts.AddCutPCMCalo("80010113","00200049f9730000dge0400000","411793105f032230000","0h63103100000010"); // min pT 75 MeV
    cuts.AddCutPCMCalo("80010113","00200019f9730000dge0400000","411793105f032230000","0h63103100000010"); // min pT 100MeV
  } else if (trainConfig == 12) {
    cuts.AddCutPCMCalo("80010113","00200068f9730000dge0400000","411793105f032230000","0h63103100000010"); // TPC cluster 35%
    cuts.AddCutPCMCalo("80010113","00200066f9730000dge0400000","411793105f032230000","0h63103100000010"); // TPC cluster 70%
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0600000","411793105f032230000","0h63103100000010"); // cosPA 0.9
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0300000","411793105f032230000","0h63103100000010"); // cosPA 0.75
  } else if (trainConfig == 13) {
    cuts.AddCutPCMCalo("80010113","0020000939730000dge0400000","411793105f032230000","0h63103100000010"); // nsig electron   -4,5
    cuts.AddCutPCMCalo("80010113","0020000969730000dge0400000","411793105f032230000","0h63103100000010"); // nsig electron -2.5,4
    cuts.AddCutPCMCalo("80010113","00200009f5730000dge0400000","411793105f032230000","0h63103100000010"); // nsig pion 2,-10
    cuts.AddCutPCMCalo("80010113","00200009f1730000dge0400000","411793105f032230000","0h63103100000010"); // nsig pion 0,-10
  } else if (trainConfig == 14) {
    cuts.AddCutPCMCalo("80010113","00200009f9030000dge0400000","411793105f032230000","0h63103100000010"); // pion nsig min mom 0.50 GeV/c
    cuts.AddCutPCMCalo("80010113","00200009f9630000dge0400000","411793105f032230000","0h63103100000010"); // pion nsig min mom 0.25 GeV/c
    cuts.AddCutPCMCalo("80010113","00200009f9760000dge0400000","411793105f032230000","0h63103100000010"); // pion nsig max mom 2.00 GeV/c
    cuts.AddCutPCMCalo("80010113","00200009f9710000dge0400000","411793105f032230000","0h63103100000010"); // pion nsig max mom 5.00 GeV/c
  } else if (trainConfig == 15) {
    cuts.AddCutPCMCalo("80010113","00200009f9730000age0400000","411793105f032230000","0h63103100000010"); // qT max 0.040, qt pt max 0.11
    cuts.AddCutPCMCalo("80010113","00200009f9730000ege0400000","411793105f032230000","0h63103100000010"); // qT max 0.060, qt pt max 0.14
    cuts.AddCutPCMCalo("80010113","00200009f9730000fge0400000","411793105f032230000","0h63103100000010"); // qT max 0.070, qt pt max 0.16
  } else if (trainConfig == 16) {
    cuts.AddCutPCMCalo("80010113","00200009f9730000d1e0400000","411793105f032230000","0h63103100000010"); // chi2 50 no chi2 dep.
    cuts.AddCutPCMCalo("80010113","00200009f9730000dfe0400000","411793105f032230000","0h63103100000010"); // chi2 50 chi2 dep -0.065
    cuts.AddCutPCMCalo("80010113","00200009f9730000dhe0400000","411793105f032230000","0h63103100000010"); // chi2 50 chi2 dep -0.050
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0404000","411793105f032230000","0h63103100000010"); // reject close v0
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0406000","411793105f032230000","0h63103100000010"); // double count with open angle 0.04
  } else if (trainConfig == 17) {
    cuts.AddCutPCMCalo("80010113","00200009f9730000dgd0400000","411793105f032230000","0h63103100000010"); // Psi pair 0.15 dep
    cuts.AddCutPCMCalo("80010113","00200009f9730000dgf0400000","411793105f032230000","0h63103100000010"); // Psi pair 0.20 dep
    cuts.AddCutPCMCalo("80010113","00200009f9730000dgg0400000","411793105f032230000","0h63103100000010"); // Psi pair 0.30 dep
  } else if (trainConfig == 18) {
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","411793105f032230000","0263103100000010"); // variation BG scheme track mult
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","411793105f032230000","0h63107100000010"); // alpha meson 0.85
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","411793105f032230000","0h63105100000010"); // alpha meson 0.75
    cuts.AddCutPCMCalo("80010113","00200009227300008250404000","411793105f032230000","0h63103100000010"); // old cuts (run1)
  } else if (trainConfig == 19) { // CALO variations
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","411793105f022230000","0h63103100000010"); // 0.6 GeV/c
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","411793105f042230000","0h63103100000010"); // 0.8 GeV/c
  } else if (trainConfig == 20) {
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","411793105f031230000","0h63103100000010"); // n cells >= 1
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","411793105f033230000","0h63103100000010"); // n cells >= 3
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","411793105f032200000","0h63103100000010"); // no max M02 cut
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","411793105f032250000","0h63103100000010"); // M02 < 0.3
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","411793105f0322k0000","0h63103100000010"); // M02, pT-dep
  } else if (trainConfig == 21) {
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","411793105e032230000","0h63103100000010"); // TM var e/p 2.0
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","411793105g032230000","0h63103100000010"); // TM var e/p 1.5
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","411793105h032230000","0h63103100000010"); // TM var e/p 1.25
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","4117931057032230000","0h63103100000010"); // TM var no veto
  } else if (trainConfig == 22) {
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","411793205f032230000","0h63103100000010"); // NL 32
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","411793305f032230000","0h63103100000010"); // NL 33
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","411793405f032230000","0h63103100000010"); // NL 34
  } else if (trainConfig == 23) {
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","411793106f032230000","0h63103100000010"); // 30/35ns
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","411793104f032230000","0h63103100000010"); // 100ns
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","411793107f032230000","0h63103100000010"); // 30ns

  //************************************************ PCM- EDC analysis 5 TeV pPb EG2 sys *********************************
  } else if (trainConfig == 40) { // PCM variations
    cuts.AddCutPCMCalo("80085113","00200009f9730000dge0400000","411793105f032230000","0h63103100000010"); //New standard cut for 8TeV analysis for RpA
    cuts.AddCutPCMCalo("80085013","00200009f9730000dge0400000","411793105f032230000","0h63103100000010"); // no SPD pileup cut
    cuts.AddCutPCMCalo("80085113","00100009f9730000dge0400000","411793105f032230000","0h63103100000010"); // R cut 2.8 -180 cm
    cuts.AddCutPCMCalo("80085113","00500009f9730000dge0400000","411793105f032230000","0h63103100000010"); // R cut 10. -180 cm
  } else if (trainConfig == 41) {
    cuts.AddCutPCMCalo("80085113","00200069f9730000dge0400000","411793105f032230000","0h63103100000010"); // min pT 40 MeV
    cuts.AddCutPCMCalo("80085113","00200049f9730000dge0400000","411793105f032230000","0h63103100000010"); // min pT 75 MeV
    cuts.AddCutPCMCalo("80085113","00200019f9730000dge0400000","411793105f032230000","0h63103100000010"); // min pT 100MeV
  } else if (trainConfig == 42) {
    cuts.AddCutPCMCalo("80085113","00200068f9730000dge0400000","411793105f032230000","0h63103100000010"); // TPC cluster 35%
    cuts.AddCutPCMCalo("80085113","00200066f9730000dge0400000","411793105f032230000","0h63103100000010"); // TPC cluster 70%
    cuts.AddCutPCMCalo("80085113","00200009f9730000dge0600000","411793105f032230000","0h63103100000010"); // cosPA 0.9
    cuts.AddCutPCMCalo("80085113","00200009f9730000dge0300000","411793105f032230000","0h63103100000010"); // cosPA 0.75
  } else if (trainConfig == 43) {
    cuts.AddCutPCMCalo("80085113","0020000939730000dge0400000","411793105f032230000","0h63103100000010"); // nsig electron   -4,5
    cuts.AddCutPCMCalo("80085113","0020000969730000dge0400000","411793105f032230000","0h63103100000010"); // nsig electron -2.5,4
    cuts.AddCutPCMCalo("80085113","00200009f5730000dge0400000","411793105f032230000","0h63103100000010"); // nsig pion 2,-10
    cuts.AddCutPCMCalo("80085113","00200009f1730000dge0400000","411793105f032230000","0h63103100000010"); // nsig pion 0,-10
  } else if (trainConfig == 44) {
    cuts.AddCutPCMCalo("80085113","00200009f9030000dge0400000","411793105f032230000","0h63103100000010"); // pion nsig min mom 0.50 GeV/c
    cuts.AddCutPCMCalo("80085113","00200009f9630000dge0400000","411793105f032230000","0h63103100000010"); // pion nsig min mom 0.25 GeV/c
    cuts.AddCutPCMCalo("80085113","00200009f9760000dge0400000","411793105f032230000","0h63103100000010"); // pion nsig max mom 2.00 GeV/c
    cuts.AddCutPCMCalo("80085113","00200009f9710000dge0400000","411793105f032230000","0h63103100000010"); // pion nsig max mom 5.00 GeV/c
  } else if (trainConfig == 45) {
    cuts.AddCutPCMCalo("80085113","00200009f9730000age0400000","411793105f032230000","0h63103100000010"); // qT max 0.040, qt pt max 0.11
    cuts.AddCutPCMCalo("80085113","00200009f9730000ege0400000","411793105f032230000","0h63103100000010"); // qT max 0.060, qt pt max 0.14
    cuts.AddCutPCMCalo("80085113","00200009f9730000fge0400000","411793105f032230000","0h63103100000010"); // qT max 0.070, qt pt max 0.16
  } else if (trainConfig == 46) {
    cuts.AddCutPCMCalo("80085113","00200009f9730000d1e0400000","411793105f032230000","0h63103100000010"); // chi2 50 no chi2 dep.
    cuts.AddCutPCMCalo("80085113","00200009f9730000dfe0400000","411793105f032230000","0h63103100000010"); // chi2 50 chi2 dep -0.065
    cuts.AddCutPCMCalo("80085113","00200009f9730000dhe0400000","411793105f032230000","0h63103100000010"); // chi2 50 chi2 dep -0.050
    cuts.AddCutPCMCalo("80085113","00200009f9730000dge0404000","411793105f032230000","0h63103100000010"); // reject close v0
    cuts.AddCutPCMCalo("80085113","00200009f9730000dge0406000","411793105f032230000","0h63103100000010"); // double count with open angle 0.04
  } else if (trainConfig == 47) {
    cuts.AddCutPCMCalo("80085113","00200009f9730000dgd0400000","411793105f032230000","0h63103100000010"); // Psi pair 0.15 dep
    cuts.AddCutPCMCalo("80085113","00200009f9730000dgf0400000","411793105f032230000","0h63103100000010"); // Psi pair 0.20 dep
    cuts.AddCutPCMCalo("80085113","00200009f9730000dgg0400000","411793105f032230000","0h63103100000010"); // Psi pair 0.30 dep
  } else if (trainConfig == 48) {
    cuts.AddCutPCMCalo("80085113","00200009f9730000dge0400000","411793105f032230000","0263103100000010"); // variation BG scheme track mult
    cuts.AddCutPCMCalo("80085113","00200009f9730000dge0400000","411793105f032230000","0h63107100000010"); // alpha meson 0.85
    cuts.AddCutPCMCalo("80085113","00200009f9730000dge0400000","411793105f032230000","0h63105100000010"); // alpha meson 0.75
    cuts.AddCutPCMCalo("80085113","00200009227300008250404000","411793105f032230000","0h63103100000010"); // old cuts (run1)
  } else if (trainConfig == 49) { // CALO variations
    cuts.AddCutPCMCalo("80085113","00200009f9730000dge0400000","411793105f022230000","0h63103100000010"); // 0.6 GeV/c
    cuts.AddCutPCMCalo("80085113","00200009f9730000dge0400000","411793105f042230000","0h63103100000010"); // 0.8 GeV/c
  } else if (trainConfig == 50) {
    cuts.AddCutPCMCalo("80085113","00200009f9730000dge0400000","411793105f031230000","0h63103100000010"); // n cells >= 1
    cuts.AddCutPCMCalo("80085113","00200009f9730000dge0400000","411793105f033230000","0h63103100000010"); // n cells >= 3
    cuts.AddCutPCMCalo("80085113","00200009f9730000dge0400000","411793105f032200000","0h63103100000010"); // no max M02 cut
    cuts.AddCutPCMCalo("80085113","00200009f9730000dge0400000","411793105f032250000","0h63103100000010"); // M02 < 0.3
    cuts.AddCutPCMCalo("80085113","00200009f9730000dge0400000","411793105f0322k0000","0h63103100000010"); // M02, pT-dep
  } else if (trainConfig == 51) {
    cuts.AddCutPCMCalo("80085113","00200009f9730000dge0400000","411793105e032230000","0h63103100000010"); // TM var e/p 2.0
    cuts.AddCutPCMCalo("80085113","00200009f9730000dge0400000","411793105g032230000","0h63103100000010"); // TM var e/p 1.5
    cuts.AddCutPCMCalo("80085113","00200009f9730000dge0400000","411793105h032230000","0h63103100000010"); // TM var e/p 1.25
    cuts.AddCutPCMCalo("80085113","00200009f9730000dge0400000","4117931057032230000","0h63103100000010"); // TM var no veto
  } else if (trainConfig == 52) {
    cuts.AddCutPCMCalo("80085113","00200009f9730000dge0400000","411793205f032230000","0h63103100000010"); // NL 32
    cuts.AddCutPCMCalo("80085113","00200009f9730000dge0400000","411793305f032230000","0h63103100000010"); // NL 33
    cuts.AddCutPCMCalo("80085113","00200009f9730000dge0400000","411793405f032230000","0h63103100000010"); // NL 34
  } else if (trainConfig == 53) {
    cuts.AddCutPCMCalo("80085113","00200009f9730000dge0400000","411793106f032230000","0h63103100000010"); // 30/35ns
    cuts.AddCutPCMCalo("80085113","00200009f9730000dge0400000","411793104f032230000","0h63103100000010"); // 100ns
    cuts.AddCutPCMCalo("80085113","00200009f9730000dge0400000","411793107f032230000","0h63103100000010"); // 30ns

    //************************************************ PCM- EDC analysis 5 TeV pPb EG1 sys *********************************
  } else if (trainConfig == 70) { // PCM variations
    cuts.AddCutPCMCalo("80083113","00200009f9730000dge0400000","411793105f032230000","0h63103100000010"); //New standard cut for 8TeV analysis for RpA
    cuts.AddCutPCMCalo("80083013","00200009f9730000dge0400000","411793105f032230000","0h63103100000010"); // no SPD pileup cut
    cuts.AddCutPCMCalo("80083113","00100009f9730000dge0400000","411793105f032230000","0h63103100000010"); // R cut 2.8 -180 cm
    cuts.AddCutPCMCalo("80083113","00500009f9730000dge0400000","411793105f032230000","0h63103100000010"); // R cut 10. -180 cm
  } else if (trainConfig == 71) {
    cuts.AddCutPCMCalo("80083113","00200069f9730000dge0400000","411793105f032230000","0h63103100000010"); // min pT 40 MeV
    cuts.AddCutPCMCalo("80083113","00200049f9730000dge0400000","411793105f032230000","0h63103100000010"); // min pT 75 MeV
    cuts.AddCutPCMCalo("80083113","00200019f9730000dge0400000","411793105f032230000","0h63103100000010"); // min pT 100MeV
  } else if (trainConfig == 72) {
    cuts.AddCutPCMCalo("80083113","00200068f9730000dge0400000","411793105f032230000","0h63103100000010"); // TPC cluster 35%
    cuts.AddCutPCMCalo("80083113","00200066f9730000dge0400000","411793105f032230000","0h63103100000010"); // TPC cluster 70%
    cuts.AddCutPCMCalo("80083113","00200009f9730000dge0600000","411793105f032230000","0h63103100000010"); // cosPA 0.9
    cuts.AddCutPCMCalo("80083113","00200009f9730000dge0300000","411793105f032230000","0h63103100000010"); // cosPA 0.75
  } else if (trainConfig == 73) {
    cuts.AddCutPCMCalo("80083113","0020000939730000dge0400000","411793105f032230000","0h63103100000010"); // nsig electron   -4,5
    cuts.AddCutPCMCalo("80083113","0020000969730000dge0400000","411793105f032230000","0h63103100000010"); // nsig electron -2.5,4
    cuts.AddCutPCMCalo("80083113","00200009f5730000dge0400000","411793105f032230000","0h63103100000010"); // nsig pion 2,-10
    cuts.AddCutPCMCalo("80083113","00200009f1730000dge0400000","411793105f032230000","0h63103100000010"); // nsig pion 0,-10
  } else if (trainConfig == 74) {
    cuts.AddCutPCMCalo("80083113","00200009f9030000dge0400000","411793105f032230000","0h63103100000010"); // pion nsig min mom 0.50 GeV/c
    cuts.AddCutPCMCalo("80083113","00200009f9630000dge0400000","411793105f032230000","0h63103100000010"); // pion nsig min mom 0.25 GeV/c
    cuts.AddCutPCMCalo("80083113","00200009f9760000dge0400000","411793105f032230000","0h63103100000010"); // pion nsig max mom 2.00 GeV/c
    cuts.AddCutPCMCalo("80083113","00200009f9710000dge0400000","411793105f032230000","0h63103100000010"); // pion nsig max mom 5.00 GeV/c
  } else if (trainConfig == 75) {
    cuts.AddCutPCMCalo("80083113","00200009f9730000age0400000","411793105f032230000","0h63103100000010"); // qT max 0.040, qt pt max 0.11
    cuts.AddCutPCMCalo("80083113","00200009f9730000ege0400000","411793105f032230000","0h63103100000010"); // qT max 0.060, qt pt max 0.14
    cuts.AddCutPCMCalo("80083113","00200009f9730000fge0400000","411793105f032230000","0h63103100000010"); // qT max 0.070, qt pt max 0.16
  } else if (trainConfig == 76) {
    cuts.AddCutPCMCalo("80083113","00200009f9730000d1e0400000","411793105f032230000","0h63103100000010"); // chi2 50 no chi2 dep.
    cuts.AddCutPCMCalo("80083113","00200009f9730000dfe0400000","411793105f032230000","0h63103100000010"); // chi2 50 chi2 dep -0.065
    cuts.AddCutPCMCalo("80083113","00200009f9730000dhe0400000","411793105f032230000","0h63103100000010"); // chi2 50 chi2 dep -0.050
    cuts.AddCutPCMCalo("80083113","00200009f9730000dge0404000","411793105f032230000","0h63103100000010"); // reject close v0
    cuts.AddCutPCMCalo("80083113","00200009f9730000dge0406000","411793105f032230000","0h63103100000010"); // double count with open angle 0.04
  } else if (trainConfig == 77) {
    cuts.AddCutPCMCalo("80083113","00200009f9730000dgd0400000","411793105f032230000","0h63103100000010"); // Psi pair 0.15 dep
    cuts.AddCutPCMCalo("80083113","00200009f9730000dgf0400000","411793105f032230000","0h63103100000010"); // Psi pair 0.20 dep
    cuts.AddCutPCMCalo("80083113","00200009f9730000dgg0400000","411793105f032230000","0h63103100000010"); // Psi pair 0.30 dep
  } else if (trainConfig == 78) {
    cuts.AddCutPCMCalo("80083113","00200009f9730000dge0400000","411793105f032230000","0263103100000010"); // variation BG scheme track mult
    cuts.AddCutPCMCalo("80083113","00200009f9730000dge0400000","411793105f032230000","0h63107100000010"); // alpha meson 0.85
    cuts.AddCutPCMCalo("80083113","00200009f9730000dge0400000","411793105f032230000","0h63105100000010"); // alpha meson 0.75
    cuts.AddCutPCMCalo("80083113","00200009227300008250404000","411793105f032230000","0h63103100000010"); // old cuts (run1)
  } else if (trainConfig == 79) { // CALO variations
    cuts.AddCutPCMCalo("80083113","00200009f9730000dge0400000","411793105f022230000","0h63103100000010"); // 0.6 GeV/c
    cuts.AddCutPCMCalo("80083113","00200009f9730000dge0400000","411793105f042230000","0h63103100000010"); // 0.8 GeV/c
  } else if (trainConfig == 80) {
    cuts.AddCutPCMCalo("80083113","00200009f9730000dge0400000","411793105f031230000","0h63103100000010"); // n cells >= 1
    cuts.AddCutPCMCalo("80083113","00200009f9730000dge0400000","411793105f033230000","0h63103100000010"); // n cells >= 3
    cuts.AddCutPCMCalo("80083113","00200009f9730000dge0400000","411793105f032200000","0h63103100000010"); // no max M02 cut
    cuts.AddCutPCMCalo("80083113","00200009f9730000dge0400000","411793105f032250000","0h63103100000010"); // M02 < 0.3
    cuts.AddCutPCMCalo("80083113","00200009f9730000dge0400000","411793105f0322k0000","0h63103100000010"); // M02, pT-dep
  } else if (trainConfig == 81) {
    cuts.AddCutPCMCalo("80083113","00200009f9730000dge0400000","411793105e032230000","0h63103100000010"); // TM var e/p 2.0
    cuts.AddCutPCMCalo("80083113","00200009f9730000dge0400000","411793105g032230000","0h63103100000010"); // TM var e/p 1.5
    cuts.AddCutPCMCalo("80083113","00200009f9730000dge0400000","411793105h032230000","0h63103100000010"); // TM var e/p 1.25
    cuts.AddCutPCMCalo("80083113","00200009f9730000dge0400000","4117931057032230000","0h63103100000010"); // TM var no veto
  } else if (trainConfig == 82) {
    cuts.AddCutPCMCalo("80083113","00200009f9730000dge0400000","411793205f032230000","0h63103100000010"); // NL 32
    cuts.AddCutPCMCalo("80083113","00200009f9730000dge0400000","411793305f032230000","0h63103100000010"); // NL 33
    cuts.AddCutPCMCalo("80083113","00200009f9730000dge0400000","411793405f032230000","0h63103100000010"); // NL 34
  } else if (trainConfig == 83) {
    cuts.AddCutPCMCalo("80083113","00200009f9730000dge0400000","411793106f032230000","0h63103100000010"); // 30/35ns
    cuts.AddCutPCMCalo("80083113","00200009f9730000dge0400000","411793104f032230000","0h63103100000010"); // 100ns
    cuts.AddCutPCMCalo("80083113","00200009f9730000dge0400000","411793107f032230000","0h63103100000010"); // 30ns


  //************************************************ PCM- EDC analysis 5 TeV pPb cent dep ************************************
  } else if (trainConfig == 100){ // EMC  INT7 run1 & run2
    cuts.AddCutPCMCalo("a0110113","00200009f9730000dge0400000","411793105f032230000","0h63103100000010"); // 0-5% PCM NL
    cuts.AddCutPCMCalo("a1210113","00200009f9730000dge0400000","411793105f032230000","0h63103100000010"); // 5-10% PCM NL
    cuts.AddCutPCMCalo("81210113","00200009f9730000dge0400000","411793105f032230000","0h63103100000010"); // 10-20% PCM NL
  } else if (trainConfig == 101){ // EMC  INT7 run1 & run2
    cuts.AddCutPCMCalo("a0110113","00200009f9730000dge0400000","111113105f032230000","0h63103100000010"); // 0-5% PCM NL
    cuts.AddCutPCMCalo("a1210113","00200009f9730000dge0400000","111113105f032230000","0h63103100000010"); // 5-10% PCM NL
    cuts.AddCutPCMCalo("81210113","00200009f9730000dge0400000","111113105f032230000","0h63103100000010"); // 10-20% PCM NL
  } else if (trainConfig == 102){ // EMC  EMC7 run1
    cuts.AddCutPCMCalo("a0152113","00200009f9730000dge0400000","111113105f032230000","0h63103100000010"); // 0-5% PCM NL
    cuts.AddCutPCMCalo("a1252113","00200009f9730000dge0400000","111113105f032230000","0h63103100000010"); // 5-10% PCM NL
    cuts.AddCutPCMCalo("81252113","00200009f9730000dge0400000","111113105f032230000","0h63103100000010"); // 10-20% PCM NL
  } else if (trainConfig == 103){ // EMC  EMC7 run1
    cuts.AddCutPCMCalo("a0185113","00200009f9730000dge0400000","111113105f032230000","0h63103100000010"); // 0-5% PCM NL
    cuts.AddCutPCMCalo("a1285113","00200009f9730000dge0400000","111113105f032230000","0h63103100000010"); // 5-10% PCM NL
    cuts.AddCutPCMCalo("81285113","00200009f9730000dge0400000","111113105f032230000","0h63103100000010"); // 10-20% PCM NL
  } else if (trainConfig == 104){ // EMC  EMC7 run1
    cuts.AddCutPCMCalo("a0183113","00200009f9730000dge0400000","111113105f032230000","0h63103100000010"); // 0-5% PCM NL
    cuts.AddCutPCMCalo("a1283113","00200009f9730000dge0400000","111113105f032230000","0h63103100000010"); // 5-10% PCM NL
    cuts.AddCutPCMCalo("81283113","00200009f9730000dge0400000","111113105f032230000","0h63103100000010"); // 10-20% PCM NL
  } else if (trainConfig == 105){ // EMC  INT7 run1 & run2
    cuts.AddCutPCMCalo("82410113","00200009f9730000dge0400000","411793105f032230000","0h63103100000010"); // 20-40% PCM NL
    cuts.AddCutPCMCalo("84610113","00200009f9730000dge0400000","411793105f032230000","0h63103100000010"); // 40-60% PCM NL
    cuts.AddCutPCMCalo("86810113","00200009f9730000dge0400000","411793105f032230000","0h63103100000010"); // 60-80% PCM NL
    cuts.AddCutPCMCalo("88010113","00200009f9730000dge0400000","411793105f032230000","0h63103100000010"); // 80-100% PCM NL
  } else if (trainConfig == 106){ // EMC  INT7 run1 & run2
    cuts.AddCutPCMCalo("82410113","00200009f9730000dge0400000","111113105f032230000","0h63103100000010"); // 20-40% PCM NL
    cuts.AddCutPCMCalo("84610113","00200009f9730000dge0400000","111113105f032230000","0h63103100000010"); // 40-60% PCM NL
    cuts.AddCutPCMCalo("86810113","00200009f9730000dge0400000","111113105f032230000","0h63103100000010"); // 60-80% PCM NL
    cuts.AddCutPCMCalo("88010113","00200009f9730000dge0400000","111113105f032230000","0h63103100000010"); // 80-100% PCM NL
  } else if (trainConfig == 107){ // EMC  EMC7 run1
    cuts.AddCutPCMCalo("82452113","00200009f9730000dge0400000","111113105f032230000","0h63103100000010"); // 20-40% PCM NL
    cuts.AddCutPCMCalo("84652113","00200009f9730000dge0400000","111113105f032230000","0h63103100000010"); // 40-60% PCM NL
    cuts.AddCutPCMCalo("86852113","00200009f9730000dge0400000","111113105f032230000","0h63103100000010"); // 60-80% PCM NL
    cuts.AddCutPCMCalo("88052113","00200009f9730000dge0400000","111113105f032230000","0h63103100000010"); // 80-100% PCM NL
  } else if (trainConfig == 108){ // EMC  EMC7 run1
    cuts.AddCutPCMCalo("82485113","00200009f9730000dge0400000","111113105f032230000","0h63103100000010"); // 20-40% PCM NL
    cuts.AddCutPCMCalo("84685113","00200009f9730000dge0400000","111113105f032230000","0h63103100000010"); // 40-60% PCM NL
    cuts.AddCutPCMCalo("86885113","00200009f9730000dge0400000","111113105f032230000","0h63103100000010"); // 60-80% PCM NL
    cuts.AddCutPCMCalo("88085113","00200009f9730000dge0400000","111113105f032230000","0h63103100000010"); // 80-100% PCM NL
  } else if (trainConfig == 109){ // EMC  EMC7 run1
    cuts.AddCutPCMCalo("82483113","00200009f9730000dge0400000","111113105f032230000","0h63103100000010"); // 20-40% PCM NL
    cuts.AddCutPCMCalo("84683113","00200009f9730000dge0400000","111113105f032230000","0h63103100000010"); // 40-60% PCM NL
    cuts.AddCutPCMCalo("86883113","00200009f9730000dge0400000","111113105f032230000","0h63103100000010"); // 60-80% PCM NL
    cuts.AddCutPCMCalo("88083113","00200009f9730000dge0400000","111113105f032230000","0h63103100000010"); // 80-100% PCM NL

  //************************************************ PCM- EDC analysis 5 TeV pPb INT7 sys *********************************
  } else if (trainConfig == 110) { // PCM variations
    cuts.AddCutPCMCalo("80010123","00200009f9730000dge0400000","411793105f032230000","0h63103100000010"); //New standard cut for 8TeV analysis for RpA
    cuts.AddCutPCMCalo("80010023","00200009f9730000dge0400000","411793105f032230000","0h63103100000010"); // no SPD pileup cut
    cuts.AddCutPCMCalo("80010123","00100009f9730000dge0400000","411793105f032230000","0h63103100000010"); // R cut 2.8 -180 cm
    cuts.AddCutPCMCalo("80010123","00500009f9730000dge0400000","411793105f032230000","0h63103100000010"); // R cut 10. -180 cm
  } else if (trainConfig == 111) {
    cuts.AddCutPCMCalo("80010123","00200069f9730000dge0400000","411793105f032230000","0h63103100000010"); // min pT 40 MeV
    cuts.AddCutPCMCalo("80010123","00200049f9730000dge0400000","411793105f032230000","0h63103100000010"); // min pT 75 MeV
    cuts.AddCutPCMCalo("80010123","00200019f9730000dge0400000","411793105f032230000","0h63103100000010"); // min pT 100MeV
  } else if (trainConfig == 112) {
    cuts.AddCutPCMCalo("80010123","00200068f9730000dge0400000","411793105f032230000","0h63103100000010"); // TPC cluster 35%
    cuts.AddCutPCMCalo("80010123","00200066f9730000dge0400000","411793105f032230000","0h63103100000010"); // TPC cluster 70%
    cuts.AddCutPCMCalo("80010123","00200009f9730000dge0600000","411793105f032230000","0h63103100000010"); // cosPA 0.9
    cuts.AddCutPCMCalo("80010123","00200009f9730000dge0300000","411793105f032230000","0h63103100000010"); // cosPA 0.75
  } else if (trainConfig == 113) {
    cuts.AddCutPCMCalo("80010123","0020000939730000dge0400000","411793105f032230000","0h63103100000010"); // nsig electron   -4,5
    cuts.AddCutPCMCalo("80010123","0020000969730000dge0400000","411793105f032230000","0h63103100000010"); // nsig electron -2.5,4
    cuts.AddCutPCMCalo("80010123","00200009f5730000dge0400000","411793105f032230000","0h63103100000010"); // nsig pion 2,-10
    cuts.AddCutPCMCalo("80010123","00200009f1730000dge0400000","411793105f032230000","0h63103100000010"); // nsig pion 0,-10
  } else if (trainConfig == 114) {
    cuts.AddCutPCMCalo("80010123","00200009f9030000dge0400000","411793105f032230000","0h63103100000010"); // pion nsig min mom 0.50 GeV/c
    cuts.AddCutPCMCalo("80010123","00200009f9630000dge0400000","411793105f032230000","0h63103100000010"); // pion nsig min mom 0.25 GeV/c
    cuts.AddCutPCMCalo("80010123","00200009f9760000dge0400000","411793105f032230000","0h63103100000010"); // pion nsig max mom 2.00 GeV/c
    cuts.AddCutPCMCalo("80010123","00200009f9710000dge0400000","411793105f032230000","0h63103100000010"); // pion nsig max mom 5.00 GeV/c
  } else if (trainConfig == 115) {
    cuts.AddCutPCMCalo("80010123","00200009f9730000age0400000","411793105f032230000","0h63103100000010"); // qT max 0.040, qt pt max 0.11
    cuts.AddCutPCMCalo("80010123","00200009f9730000ege0400000","411793105f032230000","0h63103100000010"); // qT max 0.060, qt pt max 0.14
    cuts.AddCutPCMCalo("80010123","00200009f9730000fge0400000","411793105f032230000","0h63103100000010"); // qT max 0.070, qt pt max 0.16
  } else if (trainConfig == 116) {
    cuts.AddCutPCMCalo("80010123","00200009f9730000d1e0400000","411793105f032230000","0h63103100000010"); // chi2 50 no chi2 dep.
    cuts.AddCutPCMCalo("80010123","00200009f9730000dfe0400000","411793105f032230000","0h63103100000010"); // chi2 50 chi2 dep -0.065
    cuts.AddCutPCMCalo("80010123","00200009f9730000dhe0400000","411793105f032230000","0h63103100000010"); // chi2 50 chi2 dep -0.050
    cuts.AddCutPCMCalo("80010123","00200009f9730000dge0404000","411793105f032230000","0h63103100000010"); // reject close v0
    cuts.AddCutPCMCalo("80010123","00200009f9730000dge0406000","411793105f032230000","0h63103100000010"); // double count with open angle 0.04
  } else if (trainConfig == 117) {
    cuts.AddCutPCMCalo("80010123","00200009f9730000dgd0400000","411793105f032230000","0h63103100000010"); // Psi pair 0.15 dep
    cuts.AddCutPCMCalo("80010123","00200009f9730000dgf0400000","411793105f032230000","0h63103100000010"); // Psi pair 0.20 dep
    cuts.AddCutPCMCalo("80010123","00200009f9730000dgg0400000","411793105f032230000","0h63103100000010"); // Psi pair 0.30 dep
  } else if (trainConfig == 118) {
    cuts.AddCutPCMCalo("80010123","00200009f9730000dge0400000","411793105f032230000","0263103100000010"); // variation BG scheme track mult
    cuts.AddCutPCMCalo("80010123","00200009f9730000dge0400000","411793105f032230000","0h63107100000010"); // alpha meson 0.85
    cuts.AddCutPCMCalo("80010123","00200009f9730000dge0400000","411793105f032230000","0h63105100000010"); // alpha meson 0.75
    cuts.AddCutPCMCalo("80010123","00200009227300008250404000","411793105f032230000","0h63103100000010"); // old cuts (run1)
  } else if (trainConfig == 119) { // CALO variations
    cuts.AddCutPCMCalo("80010123","00200009f9730000dge0400000","411793105f022230000","0h63103100000010"); // 0.6 GeV/c
    cuts.AddCutPCMCalo("80010123","00200009f9730000dge0400000","411793105f042230000","0h63103100000010"); // 0.8 GeV/c
  } else if (trainConfig == 120) {
    cuts.AddCutPCMCalo("80010123","00200009f9730000dge0400000","411793105f031230000","0h63103100000010"); // n cells >= 1
    cuts.AddCutPCMCalo("80010123","00200009f9730000dge0400000","411793105f033230000","0h63103100000010"); // n cells >= 3
    cuts.AddCutPCMCalo("80010123","00200009f9730000dge0400000","411793105f032200000","0h63103100000010"); // no max M02 cut
    cuts.AddCutPCMCalo("80010123","00200009f9730000dge0400000","411793105f032250000","0h63103100000010"); // M02 < 0.3
    cuts.AddCutPCMCalo("80010123","00200009f9730000dge0400000","411793105f0322k0000","0h63103100000010"); // M02, pT-dep
  } else if (trainConfig == 121) {
    cuts.AddCutPCMCalo("80010123","00200009f9730000dge0400000","411793105e032230000","0h63103100000010"); // TM var e/p 2.0
    cuts.AddCutPCMCalo("80010123","00200009f9730000dge0400000","411793105g032230000","0h63103100000010"); // TM var e/p 1.5
    cuts.AddCutPCMCalo("80010123","00200009f9730000dge0400000","411793105h032230000","0h63103100000010"); // TM var e/p 1.25
    cuts.AddCutPCMCalo("80010123","00200009f9730000dge0400000","4117931057032230000","0h63103100000010"); // TM var no veto
  } else if (trainConfig == 122) {
    cuts.AddCutPCMCalo("80010123","00200009f9730000dge0400000","411793205f032230000","0h63103100000010"); // NL 32
    cuts.AddCutPCMCalo("80010123","00200009f9730000dge0400000","411793305f032230000","0h63103100000010"); // NL 33
    cuts.AddCutPCMCalo("80010123","00200009f9730000dge0400000","411793405f032230000","0h63103100000010"); // NL 34
  } else if (trainConfig == 123) {
    cuts.AddCutPCMCalo("80010123","00200009f9730000dge0400000","411793106f032230000","0h63103100000010"); // 30/35ns
    cuts.AddCutPCMCalo("80010123","00200009f9730000dge0400000","411793104f032230000","0h63103100000010"); // 100ns
    cuts.AddCutPCMCalo("80010123","00200009f9730000dge0400000","411793107f032230000","0h63103100000010"); // 30ns

  //************************************************ PCM- EDC analysis 5 TeV pPb EG2 sys *********************************
  } else if (trainConfig == 140) { // PCM variations
    cuts.AddCutPCMCalo("80085123","00200009f9730000dge0400000","411793105f032230000","0h63103100000010"); //New standard cut for 8TeV analysis for RpA
    cuts.AddCutPCMCalo("80085023","00200009f9730000dge0400000","411793105f032230000","0h63103100000010"); // no SPD pileup cut
    cuts.AddCutPCMCalo("80085123","00100009f9730000dge0400000","411793105f032230000","0h63103100000010"); // R cut 2.8 -180 cm
    cuts.AddCutPCMCalo("80085123","00500009f9730000dge0400000","411793105f032230000","0h63103100000010"); // R cut 10. -180 cm
  } else if (trainConfig == 141) {
    cuts.AddCutPCMCalo("80085123","00200069f9730000dge0400000","411793105f032230000","0h63103100000010"); // min pT 40 MeV
    cuts.AddCutPCMCalo("80085123","00200049f9730000dge0400000","411793105f032230000","0h63103100000010"); // min pT 75 MeV
    cuts.AddCutPCMCalo("80085123","00200019f9730000dge0400000","411793105f032230000","0h63103100000010"); // min pT 100MeV
  } else if (trainConfig == 142) {
    cuts.AddCutPCMCalo("80085123","00200068f9730000dge0400000","411793105f032230000","0h63103100000010"); // TPC cluster 35%
    cuts.AddCutPCMCalo("80085123","00200066f9730000dge0400000","411793105f032230000","0h63103100000010"); // TPC cluster 70%
    cuts.AddCutPCMCalo("80085123","00200009f9730000dge0600000","411793105f032230000","0h63103100000010"); // cosPA 0.9
    cuts.AddCutPCMCalo("80085123","00200009f9730000dge0300000","411793105f032230000","0h63103100000010"); // cosPA 0.75
  } else if (trainConfig == 143) {
    cuts.AddCutPCMCalo("80085123","0020000939730000dge0400000","411793105f032230000","0h63103100000010"); // nsig electron   -4,5
    cuts.AddCutPCMCalo("80085123","0020000969730000dge0400000","411793105f032230000","0h63103100000010"); // nsig electron -2.5,4
    cuts.AddCutPCMCalo("80085123","00200009f5730000dge0400000","411793105f032230000","0h63103100000010"); // nsig pion 2,-10
    cuts.AddCutPCMCalo("80085123","00200009f1730000dge0400000","411793105f032230000","0h63103100000010"); // nsig pion 0,-10
  } else if (trainConfig == 144) {
    cuts.AddCutPCMCalo("80085123","00200009f9030000dge0400000","411793105f032230000","0h63103100000010"); // pion nsig min mom 0.50 GeV/c
    cuts.AddCutPCMCalo("80085123","00200009f9630000dge0400000","411793105f032230000","0h63103100000010"); // pion nsig min mom 0.25 GeV/c
    cuts.AddCutPCMCalo("80085123","00200009f9760000dge0400000","411793105f032230000","0h63103100000010"); // pion nsig max mom 2.00 GeV/c
    cuts.AddCutPCMCalo("80085123","00200009f9710000dge0400000","411793105f032230000","0h63103100000010"); // pion nsig max mom 5.00 GeV/c
  } else if (trainConfig == 145) {
    cuts.AddCutPCMCalo("80085123","00200009f9730000age0400000","411793105f032230000","0h63103100000010"); // qT max 0.040, qt pt max 0.11
    cuts.AddCutPCMCalo("80085123","00200009f9730000ege0400000","411793105f032230000","0h63103100000010"); // qT max 0.060, qt pt max 0.14
    cuts.AddCutPCMCalo("80085123","00200009f9730000fge0400000","411793105f032230000","0h63103100000010"); // qT max 0.070, qt pt max 0.16
  } else if (trainConfig == 146) {
    cuts.AddCutPCMCalo("80085123","00200009f9730000d1e0400000","411793105f032230000","0h63103100000010"); // chi2 50 no chi2 dep.
    cuts.AddCutPCMCalo("80085123","00200009f9730000dfe0400000","411793105f032230000","0h63103100000010"); // chi2 50 chi2 dep -0.065
    cuts.AddCutPCMCalo("80085123","00200009f9730000dhe0400000","411793105f032230000","0h63103100000010"); // chi2 50 chi2 dep -0.050
    cuts.AddCutPCMCalo("80085123","00200009f9730000dge0404000","411793105f032230000","0h63103100000010"); // reject close v0
    cuts.AddCutPCMCalo("80085123","00200009f9730000dge0406000","411793105f032230000","0h63103100000010"); // double count with open angle 0.04
  } else if (trainConfig == 147) {
    cuts.AddCutPCMCalo("80085123","00200009f9730000dgd0400000","411793105f032230000","0h63103100000010"); // Psi pair 0.15 dep
    cuts.AddCutPCMCalo("80085123","00200009f9730000dgf0400000","411793105f032230000","0h63103100000010"); // Psi pair 0.20 dep
    cuts.AddCutPCMCalo("80085123","00200009f9730000dgg0400000","411793105f032230000","0h63103100000010"); // Psi pair 0.30 dep
  } else if (trainConfig == 148) {
    cuts.AddCutPCMCalo("80085123","00200009f9730000dge0400000","411793105f032230000","0263103100000010"); // variation BG scheme track mult
    cuts.AddCutPCMCalo("80085123","00200009f9730000dge0400000","411793105f032230000","0h63107100000010"); // alpha meson 0.85
    cuts.AddCutPCMCalo("80085123","00200009f9730000dge0400000","411793105f032230000","0h63105100000010"); // alpha meson 0.75
    cuts.AddCutPCMCalo("80085123","00200009227300008250404000","411793105f032230000","0h63103100000010"); // old cuts (run1)
  } else if (trainConfig == 149) { // CALO variations
    cuts.AddCutPCMCalo("80085123","00200009f9730000dge0400000","411793105f022230000","0h63103100000010"); // 0.6 GeV/c
    cuts.AddCutPCMCalo("80085123","00200009f9730000dge0400000","411793105f042230000","0h63103100000010"); // 0.8 GeV/c
  } else if (trainConfig == 150) {
    cuts.AddCutPCMCalo("80085123","00200009f9730000dge0400000","411793105f031230000","0h63103100000010"); // n cells >= 1
    cuts.AddCutPCMCalo("80085123","00200009f9730000dge0400000","411793105f033230000","0h63103100000010"); // n cells >= 3
    cuts.AddCutPCMCalo("80085123","00200009f9730000dge0400000","411793105f032200000","0h63103100000010"); // no max M02 cut
    cuts.AddCutPCMCalo("80085123","00200009f9730000dge0400000","411793105f032250000","0h63103100000010"); // M02 < 0.3
    cuts.AddCutPCMCalo("80085123","00200009f9730000dge0400000","411793105f0322k0000","0h63103100000010"); // M02, pT-dep
  } else if (trainConfig == 151) {
    cuts.AddCutPCMCalo("80085123","00200009f9730000dge0400000","411793105e032230000","0h63103100000010"); // TM var e/p 2.0
    cuts.AddCutPCMCalo("80085123","00200009f9730000dge0400000","411793105g032230000","0h63103100000010"); // TM var e/p 1.5
    cuts.AddCutPCMCalo("80085123","00200009f9730000dge0400000","411793105h032230000","0h63103100000010"); // TM var e/p 1.25
    cuts.AddCutPCMCalo("80085123","00200009f9730000dge0400000","4117931057032230000","0h63103100000010"); // TM var no veto
  } else if (trainConfig == 152) {
    cuts.AddCutPCMCalo("80085123","00200009f9730000dge0400000","411793205f032230000","0h63103100000010"); // NL 32
    cuts.AddCutPCMCalo("80085123","00200009f9730000dge0400000","411793305f032230000","0h63103100000010"); // NL 33
    cuts.AddCutPCMCalo("80085123","00200009f9730000dge0400000","411793405f032230000","0h63103100000010"); // NL 34
  } else if (trainConfig == 153) {
    cuts.AddCutPCMCalo("80085123","00200009f9730000dge0400000","411793106f032230000","0h63103100000010"); // 30/35ns
    cuts.AddCutPCMCalo("80085123","00200009f9730000dge0400000","411793104f032230000","0h63103100000010"); // 100ns
    cuts.AddCutPCMCalo("80085123","00200009f9730000dge0400000","411793107f032230000","0h63103100000010"); // 30ns

    //************************************************ PCM- EDC analysis 5 TeV pPb EG1 sys *********************************
  } else if (trainConfig == 170) { // PCM variations
    cuts.AddCutPCMCalo("80083123","00200009f9730000dge0400000","411793105f032230000","0h63103100000010"); //New standard cut for 8TeV analysis for RpA
    cuts.AddCutPCMCalo("80083023","00200009f9730000dge0400000","411793105f032230000","0h63103100000010"); // no SPD pileup cut
    cuts.AddCutPCMCalo("80083123","00100009f9730000dge0400000","411793105f032230000","0h63103100000010"); // R cut 2.8 -180 cm
    cuts.AddCutPCMCalo("80083123","00500009f9730000dge0400000","411793105f032230000","0h63103100000010"); // R cut 10. -180 cm
  } else if (trainConfig == 171) {
    cuts.AddCutPCMCalo("80083123","00200069f9730000dge0400000","411793105f032230000","0h63103100000010"); // min pT 40 MeV
    cuts.AddCutPCMCalo("80083123","00200049f9730000dge0400000","411793105f032230000","0h63103100000010"); // min pT 75 MeV
    cuts.AddCutPCMCalo("80083123","00200019f9730000dge0400000","411793105f032230000","0h63103100000010"); // min pT 100MeV
  } else if (trainConfig == 172) {
    cuts.AddCutPCMCalo("80083123","00200068f9730000dge0400000","411793105f032230000","0h63103100000010"); // TPC cluster 35%
    cuts.AddCutPCMCalo("80083123","00200066f9730000dge0400000","411793105f032230000","0h63103100000010"); // TPC cluster 70%
    cuts.AddCutPCMCalo("80083123","00200009f9730000dge0600000","411793105f032230000","0h63103100000010"); // cosPA 0.9
    cuts.AddCutPCMCalo("80083123","00200009f9730000dge0300000","411793105f032230000","0h63103100000010"); // cosPA 0.75
  } else if (trainConfig == 173) {
    cuts.AddCutPCMCalo("80083123","0020000939730000dge0400000","411793105f032230000","0h63103100000010"); // nsig electron   -4,5
    cuts.AddCutPCMCalo("80083123","0020000969730000dge0400000","411793105f032230000","0h63103100000010"); // nsig electron -2.5,4
    cuts.AddCutPCMCalo("80083123","00200009f5730000dge0400000","411793105f032230000","0h63103100000010"); // nsig pion 2,-10
    cuts.AddCutPCMCalo("80083123","00200009f1730000dge0400000","411793105f032230000","0h63103100000010"); // nsig pion 0,-10
  } else if (trainConfig == 174) {
    cuts.AddCutPCMCalo("80083123","00200009f9030000dge0400000","411793105f032230000","0h63103100000010"); // pion nsig min mom 0.50 GeV/c
    cuts.AddCutPCMCalo("80083123","00200009f9630000dge0400000","411793105f032230000","0h63103100000010"); // pion nsig min mom 0.25 GeV/c
    cuts.AddCutPCMCalo("80083123","00200009f9760000dge0400000","411793105f032230000","0h63103100000010"); // pion nsig max mom 2.00 GeV/c
    cuts.AddCutPCMCalo("80083123","00200009f9710000dge0400000","411793105f032230000","0h63103100000010"); // pion nsig max mom 5.00 GeV/c
  } else if (trainConfig == 175) {
    cuts.AddCutPCMCalo("80083123","00200009f9730000age0400000","411793105f032230000","0h63103100000010"); // qT max 0.040, qt pt max 0.11
    cuts.AddCutPCMCalo("80083123","00200009f9730000ege0400000","411793105f032230000","0h63103100000010"); // qT max 0.060, qt pt max 0.14
    cuts.AddCutPCMCalo("80083123","00200009f9730000fge0400000","411793105f032230000","0h63103100000010"); // qT max 0.070, qt pt max 0.16
  } else if (trainConfig == 176) {
    cuts.AddCutPCMCalo("80083123","00200009f9730000d1e0400000","411793105f032230000","0h63103100000010"); // chi2 50 no chi2 dep.
    cuts.AddCutPCMCalo("80083123","00200009f9730000dfe0400000","411793105f032230000","0h63103100000010"); // chi2 50 chi2 dep -0.065
    cuts.AddCutPCMCalo("80083123","00200009f9730000dhe0400000","411793105f032230000","0h63103100000010"); // chi2 50 chi2 dep -0.050
    cuts.AddCutPCMCalo("80083123","00200009f9730000dge0404000","411793105f032230000","0h63103100000010"); // reject close v0
    cuts.AddCutPCMCalo("80083123","00200009f9730000dge0406000","411793105f032230000","0h63103100000010"); // double count with open angle 0.04
  } else if (trainConfig == 177) {
    cuts.AddCutPCMCalo("80083123","00200009f9730000dgd0400000","411793105f032230000","0h63103100000010"); // Psi pair 0.15 dep
    cuts.AddCutPCMCalo("80083123","00200009f9730000dgf0400000","411793105f032230000","0h63103100000010"); // Psi pair 0.20 dep
    cuts.AddCutPCMCalo("80083123","00200009f9730000dgg0400000","411793105f032230000","0h63103100000010"); // Psi pair 0.30 dep
  } else if (trainConfig == 178) {
    cuts.AddCutPCMCalo("80083123","00200009f9730000dge0400000","411793105f032230000","0263103100000010"); // variation BG scheme track mult
    cuts.AddCutPCMCalo("80083123","00200009f9730000dge0400000","411793105f032230000","0h63107100000010"); // alpha meson 0.85
    cuts.AddCutPCMCalo("80083123","00200009f9730000dge0400000","411793105f032230000","0h63105100000010"); // alpha meson 0.75
    cuts.AddCutPCMCalo("80083123","00200009227300008250404000","411793105f032230000","0h63103100000010"); // old cuts (run1)
  } else if (trainConfig == 179) { // CALO variations
    cuts.AddCutPCMCalo("80083123","00200009f9730000dge0400000","411793105f022230000","0h63103100000010"); // 0.6 GeV/c
    cuts.AddCutPCMCalo("80083123","00200009f9730000dge0400000","411793105f042230000","0h63103100000010"); // 0.8 GeV/c
  } else if (trainConfig == 180) {
    cuts.AddCutPCMCalo("80083123","00200009f9730000dge0400000","411793105f031230000","0h63103100000010"); // n cells >= 1
    cuts.AddCutPCMCalo("80083123","00200009f9730000dge0400000","411793105f033230000","0h63103100000010"); // n cells >= 3
    cuts.AddCutPCMCalo("80083123","00200009f9730000dge0400000","411793105f032200000","0h63103100000010"); // no max M02 cut
    cuts.AddCutPCMCalo("80083123","00200009f9730000dge0400000","411793105f032250000","0h63103100000010"); // M02 < 0.3
    cuts.AddCutPCMCalo("80083123","00200009f9730000dge0400000","411793105f0322k0000","0h63103100000010"); // M02, pT-dep
  } else if (trainConfig == 181) {
    cuts.AddCutPCMCalo("80083123","00200009f9730000dge0400000","411793105e032230000","0h63103100000010"); // TM var e/p 2.0
    cuts.AddCutPCMCalo("80083123","00200009f9730000dge0400000","411793105g032230000","0h63103100000010"); // TM var e/p 1.5
    cuts.AddCutPCMCalo("80083123","00200009f9730000dge0400000","411793105h032230000","0h63103100000010"); // TM var e/p 1.25
    cuts.AddCutPCMCalo("80083123","00200009f9730000dge0400000","4117931057032230000","0h63103100000010"); // TM var no veto
  } else if (trainConfig == 182) {
    cuts.AddCutPCMCalo("80083123","00200009f9730000dge0400000","411793205f032230000","0h63103100000010"); // NL 32
    cuts.AddCutPCMCalo("80083123","00200009f9730000dge0400000","411793305f032230000","0h63103100000010"); // NL 33
    cuts.AddCutPCMCalo("80083123","00200009f9730000dge0400000","411793405f032230000","0h63103100000010"); // NL 34
  } else if (trainConfig == 183) {
    cuts.AddCutPCMCalo("80083123","00200009f9730000dge0400000","411793106f032230000","0h63103100000010"); // 30/35ns
    cuts.AddCutPCMCalo("80083123","00200009f9730000dge0400000","411793104f032230000","0h63103100000010"); // 100ns
    cuts.AddCutPCMCalo("80083123","00200009f9730000dge0400000","411793107f032230000","0h63103100000010"); // 30ns
    
  //************************************************ PCM- PHOS analysis 5 TeV pPb ********************************************
  } else if (trainConfig == 1000){ // PHOS  INT7 run1
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","244440004a013200000","0h63103100000010"); // 0-100% without NL
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","244445304a013200000","0h63103100000010"); // 0-100% with NL 1
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","244445404a013200000","0h63103100000010"); // 0-100% with NL 2
  } else if (trainConfig == 1001){ // PHOS  PHI7 run1
    cuts.AddCutPCMCalo("80062113","00200009f9730000dge0400000","244440004a013200000","0h63103100000010"); // 0-100% without NL
    cuts.AddCutPCMCalo("80062113","00200009f9730000dge0400000","244445304a013200000","0h63103100000010"); // 0-100% with NL 1
    cuts.AddCutPCMCalo("80062113","00200009f9730000dge0400000","244445404a013200000","0h63103100000010"); // 0-100% with NL 2
  } else if (trainConfig == 1002) {  // PHOS  INT7 run2
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","24466000ha012200000","0h63103100000010"); // 0-100% without NL
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","24466530ha012200000","0h63103100000010"); // 0-100% without NL
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","24466540ha012200000","0h63103100000010"); // 0-100% without NL

  } else if (trainConfig == 1003){ // PHOS  INT7 run1
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","244445104a013200000","0h63103100000010"); // 0-100% PCM NL
  } else if (trainConfig == 1004){ // PHOS  PHI7 run1
    cuts.AddCutPCMCalo("80062113","00200009f9730000dge0400000","244445104a013200000","0h63103100000010"); // 0-100% PCM NL
  } else if (trainConfig == 1005) {  // PHOS  INT7 run2
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","24466410ha012200000","0h63103100000010"); // 0-100% PCM NL
  
  } else if (trainConfig == 1006) {  // PHOS  INT7 run2, new PHOS run2 default cuts (M02, cluster energy)
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","24466000ha01cc00000","0h63103100000010"); // 0-100% without NL
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","24466530ha01cc00000","0h63103100000010"); // 0-100% with NL 1
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","24466540ha01cc00000","0h63103100000010"); // 0-100% with NL 2

  //************************************************ PCM- PHOS analysis 5 TeV pPb cent dep ************************************
  } else if (trainConfig == 1100){ // centrality dependent and with latest TM run 1
    cuts.AddCutPCMCalo("a0110113","00200009f9730000dge0400000","244445104a013200000","0h63103100000010"); // 0-5% with PCM NL
    cuts.AddCutPCMCalo("a1210113","00200009f9730000dge0400000","244445104a013200000","0h63103100000010"); // 5-10% with PCM NL
    cuts.AddCutPCMCalo("81210113","00200009f9730000dge0400000","244445104a013200000","0h63103100000010"); // 10-20% with PCM NL
  } else if (trainConfig == 1101){ // centrality dependent and with latest TM run 1
    cuts.AddCutPCMCalo("a0162113","00200009f9730000dge0400000","244445104a013200000","0h63103100000010"); // 0-5% with PCM NL
    cuts.AddCutPCMCalo("a1262113","00200009f9730000dge0400000","244445104a013200000","0h63103100000010"); // 5-10% with PCM NL
    cuts.AddCutPCMCalo("81262113","00200009f9730000dge0400000","244445104a013200000","0h63103100000010"); // 10-20% with PCM NL
  } else if (trainConfig == 1103) {  // PHOS  INT7 with cents
    cuts.AddCutPCMCalo("a0110113","00200009f9730000dge0400000","24466410ha012200000","0h63103100000010"); // non lin 0-5%
    cuts.AddCutPCMCalo("a1210113","00200009f9730000dge0400000","24466410ha012200000","0h63103100000010"); // non lin 5-10%
    cuts.AddCutPCMCalo("81210113","00200009f9730000dge0400000","24466410ha012200000","0h63103100000010"); // non lin 10-20%
  } else if (trainConfig == 1104){ // centrality dependent and with latest TM run 1
    cuts.AddCutPCMCalo("82410113","00200009f9730000dge0400000","244445104a013200000","0h63103100000010"); // 20-40% with PCM NL
    cuts.AddCutPCMCalo("84610113","00200009f9730000dge0400000","244445104a013200000","0h63103100000010"); // 40-60% with PCM NL
    cuts.AddCutPCMCalo("86810113","00200009f9730000dge0400000","244445104a013200000","0h63103100000010"); // 60-80% with PCM NL
    cuts.AddCutPCMCalo("88010113","00200009f9730000dge0400000","244445104a013200000","0h63103100000010"); // 80-100% with PCM NL
  } else if (trainConfig == 1105){ // centrality dependent and with latest TM run 1
    cuts.AddCutPCMCalo("82462113","00200009f9730000dge0400000","244445104a013200000","0h63103100000010"); // 20-40% with PCM NL
    cuts.AddCutPCMCalo("84662113","00200009f9730000dge0400000","244445104a013200000","0h63103100000010"); // 40-60% with PCM NL
    cuts.AddCutPCMCalo("86862113","00200009f9730000dge0400000","244445104a013200000","0h63103100000010"); // 60-80% with PCM NL
    cuts.AddCutPCMCalo("88062113","00200009f9730000dge0400000","244445104a013200000","0h63103100000010"); // 80-100% with PCM NL
  } else if (trainConfig == 1106) {  // PHOS  INT7 with cents
    cuts.AddCutPCMCalo("82410113","00200009f9730000dge0400000","24466410ha012200000","0h63103100000010"); // non lin 20-40%
    cuts.AddCutPCMCalo("84610113","00200009f9730000dge0400000","24466410ha012200000","0h63103100000010"); // non lin 40-60%
    cuts.AddCutPCMCalo("86810113","00200009f9730000dge0400000","24466410ha012200000","0h63103100000010"); // non lin 60-80%
    cuts.AddCutPCMCalo("88010113","00200009f9730000dge0400000","24466410ha012200000","0h63103100000010"); // non lin 80-100%

  //************************************************ PCM- PHOS analysis 5 TeV pPb HBT ************************************
  } else if (trainConfig == 1520) {  // PHOS  INT7
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0404000","24466410ha012200000","0h63103100000010"); //
  } else if (trainConfig == 1521) {  // PHOS  INT7 with cents
    cuts.AddCutPCMCalo("80110113","00200009f9730000dge0404000","24466410ha012200000","0h63103100000010"); // non lin 0-10%
    cuts.AddCutPCMCalo("81210113","00200009f9730000dge0404000","24466410ha012200000","0h63103100000010"); // non lin 10-20%
    cuts.AddCutPCMCalo("82410113","00200009f9730000dge0404000","24466410ha012200000","0h63103100000010"); // non lin 20-40%
  } else if (trainConfig == 1522) {  // PHOS  INT7 with cents
    cuts.AddCutPCMCalo("84610113","00200009f9730000dge0404000","24466410ha012200000","0h63103100000010"); // non lin 40-60%
    cuts.AddCutPCMCalo("86810113","00200009f9730000dge0404000","24466410ha012200000","0h63103100000010"); // non lin 60-80%
    cuts.AddCutPCMCalo("88010113","00200009f9730000dge0404000","24466410ha012200000","0h63103100000010"); // non lin 80-100%


  //************************************************ PCM - EDC analysis 8 TeV pPb *********************************************
  // 8 TeV pPb variations with new PCM cut
  } else if (trainConfig == 2020) { // standard
    cuts.AddCutPCMCalo("80010123","0dm00009f9730000dge0404000","411793105f032230000","0h63103100000010"); //New standard cut for 8TeV analysis for RpA
    cuts.AddCutPCMCalo("8008e123","0dm00009f9730000dge0404000","411793105f032230000","0h63103100000010"); //New standard cut for 8TeV analysis for RpA
    cuts.AddCutPCMCalo("8008d123","0dm00009f9730000dge0404000","411793105f032230000","0h63103100000010"); //New standard cut for 8TeV analysis for RpA
  } else if (trainConfig == 2021) { // standard copy of 2020 for timing cut studies on clusterizer
    cuts.AddCutPCMCalo("80010123","0dm00009f9730000dge0404000","411793105f032230000","0h63103100000010");
    cuts.AddCutPCMCalo("8008e123","0dm00009f9730000dge0404000","411793105f032230000","0h63103100000010");
    cuts.AddCutPCMCalo("8008d123","0dm00009f9730000dge0404000","411793105f032230000","0h63103100000010");
  } else if (trainConfig == 2022) { // standard copy of 2020 for timing cut studies on clusterizer
    cuts.AddCutPCMCalo("80010123","0dm00009f9730000dge0404000","411793105f032230000","0h63103100000010");
    cuts.AddCutPCMCalo("8008e123","0dm00009f9730000dge0404000","411793105f032230000","0h63103100000010");
    cuts.AddCutPCMCalo("8008d123","0dm00009f9730000dge0404000","411793105f032230000","0h63103100000010");
  } else if (trainConfig == 2023) { // standard copy of 2020 for timing cut studies on clusterizer
    cuts.AddCutPCMCalo("80010123","0dm00009f9730000dge0404000","411793105f032230000","0h63103100000010");
    cuts.AddCutPCMCalo("8008e123","0dm00009f9730000dge0404000","411793105f032230000","0h63103100000010");
    cuts.AddCutPCMCalo("8008d123","0dm00009f9730000dge0404000","411793105f032230000","0h63103100000010");
  } else if (trainConfig == 2024) { // standard copy of 2020 for timing cut studies on clusterizer
    cuts.AddCutPCMCalo("80010123","0dm00009f9730000dge0404000","411793105f032230000","0h63103100000010");
    cuts.AddCutPCMCalo("8008e123","0dm00009f9730000dge0404000","411793105f032230000","0h63103100000010");
    cuts.AddCutPCMCalo("8008d123","0dm00009f9730000dge0404000","411793105f032230000","0h63103100000010");
  } else if (trainConfig == 2025) { // including eta<0.8, DC, R region rej.
    cuts.AddCutPCMCalo("80010123","0dm00009f9730000dge0404000","411793105f032230000","0h63103100000010");
    cuts.AddCutPCMCalo("8008e123","0dm00009f9730000dge0404000","411793105f032230000","0h63103100000010");
    cuts.AddCutPCMCalo("8008d123","0dm00009f9730000dge0404000","411793105f032230000","0h63103100000010");
  } else if (trainConfig == 2026) { // including eta<0.8, DC, R region rej.
    cuts.AddCutPCMCalo("80010123","0dm00009f9730000dge0404000","411793805f032230000","0h63103100000010");
    cuts.AddCutPCMCalo("8008e123","0dm00009f9730000dge0404000","411793805f032230000","0h63103100000010");
    cuts.AddCutPCMCalo("8008d123","0dm00009f9730000dge0404000","411793805f032230000","0h63103100000010");
  } else if (trainConfig == 2027) { // including eta<0.8, DC, R region rej.
    cuts.AddCutPCMCalo("80010123","0dm00009f9730000dge0404000","411793905f032230000","0h63103100000010");
    cuts.AddCutPCMCalo("8008e123","0dm00009f9730000dge0404000","411793905f032230000","0h63103100000010");
    cuts.AddCutPCMCalo("8008d123","0dm00009f9730000dge0404000","411793905f032230000","0h63103100000010");
  } else if (trainConfig == 2028) { // including eta<0.8, DC, R region rej.
    cuts.AddCutPCMCalo("80010623","0dm00009f9730000dge0404000","411793105f032230000","0h63103100000010");
    cuts.AddCutPCMCalo("8008e623","0dm00009f9730000dge0404000","411793105f032230000","0h63103100000010");
    cuts.AddCutPCMCalo("8008d623","0dm00009f9730000dge0404000","411793105f032230000","0h63103100000010");
  } else if (trainConfig == 2029) { // including eta<0.8, DC, R region rej.
    cuts.AddCutPCMCalo("80010723","0dm00009f9730000dge0404000","411793105f032230000","0h63103100000010");
    cuts.AddCutPCMCalo("8008e723","0dm00009f9730000dge0404000","411793105f032230000","0h63103100000010");
    cuts.AddCutPCMCalo("8008d723","0dm00009f9730000dge0404000","411793105f032230000","0h63103100000010");

  // Nonlin testing configs (TB only)
  } else if (trainConfig == 2030) { // NL 01 -> 100 MeV aggregation
    cuts.AddCutPCMCalo("80010123","00200009f9730000dge0400000","411790105f032230000","0h63103100000010");
  } else if (trainConfig == 2031) { // NL 01 -> 100 MeV aggregation
    cuts.AddCutPCMCalo("8008e123","00200009f9730000dge0400000","411790105f032230000","0h63103100000010");
    cuts.AddCutPCMCalo("8008d123","00200009f9730000dge0400000","411790105f032230000","0h63103100000010");
  } else if (trainConfig == 2032) { // NL 02 -> 50 MeV aggregation
    cuts.AddCutPCMCalo("80010123","00200009f9730000dge0400000","411790205f032230000","0h63103100000010");
  } else if (trainConfig == 2033) { // NL 02 -> 50 MeV aggregation
    cuts.AddCutPCMCalo("8008e123","00200009f9730000dge0400000","411790205f032230000","0h63103100000010");
    cuts.AddCutPCMCalo("8008d123","00200009f9730000dge0400000","411790205f032230000","0h63103100000010");
  } else if (trainConfig == 2034) { // NL 03 -> 150 MeV aggregation
    cuts.AddCutPCMCalo("80010123","00200009f9730000dge0400000","411790305f032230000","0h63103100000010");
  } else if (trainConfig == 2035) { // NL 03 -> 150 MeV aggregation
    cuts.AddCutPCMCalo("8008e123","00200009f9730000dge0400000","411790305f032230000","0h63103100000010");
    cuts.AddCutPCMCalo("8008d123","00200009f9730000dge0400000","411790305f032230000","0h63103100000010");
  } else if (trainConfig == 2036) { // NL 04 -> 300 MeV aggregation
    cuts.AddCutPCMCalo("80010123","00200009f9730000dge0400000","411790405f032230000","0h63103100000010");
  } else if (trainConfig == 2037) { // NL 04 -> 300 MeV aggregation
    cuts.AddCutPCMCalo("8008e123","00200009f9730000dge0400000","411790405f032230000","0h63103100000010");
    cuts.AddCutPCMCalo("8008d123","00200009f9730000dge0400000","411790405f032230000","0h63103100000010");

  // Nonlin testing configs (TB + finetuning)
  } else if (trainConfig == 2040) { // NL 01 -> 100 MeV aggregation
    cuts.AddCutPCMCalo("80010123","00200009f9730000dge0400000","411793105f032230000","0h63103100000010");
  } else if (trainConfig == 2041) { // NL 01 -> 100 MeV aggregation
    cuts.AddCutPCMCalo("8008e123","00200009f9730000dge0400000","411793105f032230000","0h63103100000010");
    cuts.AddCutPCMCalo("8008d123","00200009f9730000dge0400000","411793105f032230000","0h63103100000010");
  } else if (trainConfig == 2042) { // NL 02 -> 50 MeV aggregation
    cuts.AddCutPCMCalo("80010123","00200009f9730000dge0400000","411793205f032230000","0h63103100000010");
  } else if (trainConfig == 2043) { // NL 02 -> 50 MeV aggregation
    cuts.AddCutPCMCalo("8008e123","00200009f9730000dge0400000","411793205f032230000","0h63103100000010");
    cuts.AddCutPCMCalo("8008d123","00200009f9730000dge0400000","411793205f032230000","0h63103100000010");
  } else if (trainConfig == 2044) { // NL 03 -> 150 MeV aggregation
    cuts.AddCutPCMCalo("80010123","00200009f9730000dge0400000","411793305f032230000","0h63103100000010");
  } else if (trainConfig == 2045) { // NL 03 -> 150 MeV aggregation
    cuts.AddCutPCMCalo("8008e123","00200009f9730000dge0400000","411793305f032230000","0h63103100000010");
    cuts.AddCutPCMCalo("8008d123","00200009f9730000dge0400000","411793305f032230000","0h63103100000010");
  } else if (trainConfig == 2046) { // NL 04 -> 300 MeV aggregation
    cuts.AddCutPCMCalo("80010123","00200009f9730000dge0400000","411793405f032230000","0h63103100000010");
  } else if (trainConfig == 2047) { // NL 04 -> 300 MeV aggregation
    cuts.AddCutPCMCalo("8008e123","00200009f9730000dge0400000","411793405f032230000","0h63103100000010");
    cuts.AddCutPCMCalo("8008d123","00200009f9730000dge0400000","411793405f032230000","0h63103100000010");

  // Nonlin testing configs (TB + finetuning)
  } else if (trainConfig == 2050) { // NL 01 -> 100 MeV aggregation
    cuts.AddCutPCMCalo("80010123","00200009f9730000dge0400000","411796105f032230000","0h63103100000010");
  } else if (trainConfig == 2051) { // NL 01 -> 100 MeV aggregation
    cuts.AddCutPCMCalo("8008e123","00200009f9730000dge0400000","411796105f032230000","0h63103100000010");
    cuts.AddCutPCMCalo("8008d123","00200009f9730000dge0400000","411796105f032230000","0h63103100000010");
  } else if (trainConfig == 2052) { // NL 02 -> 50 MeV aggregation
    cuts.AddCutPCMCalo("80010123","00200009f9730000dge0400000","411796205f032230000","0h63103100000010");
  } else if (trainConfig == 2053) { // NL 02 -> 50 MeV aggregation
    cuts.AddCutPCMCalo("8008e123","00200009f9730000dge0400000","411796205f032230000","0h63103100000010");
    cuts.AddCutPCMCalo("8008d123","00200009f9730000dge0400000","411796205f032230000","0h63103100000010");
  } else if (trainConfig == 2054) { // NL 03 -> 150 MeV aggregation
    cuts.AddCutPCMCalo("80010123","00200009f9730000dge0400000","411796305f032230000","0h63103100000010");
  } else if (trainConfig == 2055) { // NL 03 -> 150 MeV aggregation
    cuts.AddCutPCMCalo("8008e123","00200009f9730000dge0400000","411796305f032230000","0h63103100000010");
    cuts.AddCutPCMCalo("8008d123","00200009f9730000dge0400000","411796305f032230000","0h63103100000010");
  } else if (trainConfig == 2056) { // NL 04 -> 300 MeV aggregation
    cuts.AddCutPCMCalo("80010123","00200009f9730000dge0400000","411796405f032230000","0h63103100000010");
  } else if (trainConfig == 2057) { // NL 04 -> 300 MeV aggregation
    cuts.AddCutPCMCalo("8008e123","00200009f9730000dge0400000","411796405f032230000","0h63103100000010");
    cuts.AddCutPCMCalo("8008d123","00200009f9730000dge0400000","411796405f032230000","0h63103100000010");
  } else if (trainConfig == 2058) { // NL 04 -> 300 MeV aggregation
    cuts.AddCutPCMCalo("80010123","00200009f9730000dge0400000","411796905f032230000","0h63103100000010");
  } else if (trainConfig == 2059) { // NL 04 -> 300 MeV aggregation
    cuts.AddCutPCMCalo("8008e123","00200009f9730000dge0400000","411796905f032230000","0h63103100000010");
    cuts.AddCutPCMCalo("8008d123","00200009f9730000dge0400000","411796905f032230000","0h63103100000010");

  } else if (trainConfig == 2060) { // no NCell cut
    cuts.AddCutPCMCalo("80010103","0dm00009f9730000dge0404000","411793105f030230000","0h63103100000010");
    cuts.AddCutPCMCalo("8008e103","0dm00009f9730000dge0404000","411793105f030230000","0h63103100000010");
    cuts.AddCutPCMCalo("8008d103","0dm00009f9730000dge0404000","411793105f030230000","0h63103100000010");

  } else if (trainConfig == 2070) { // TOF timing requirement one leg
    cuts.AddCutPCMCalo("80010123","0dm00009f9730600dge0404000","411793105f032230000","0h63103100000010");
    cuts.AddCutPCMCalo("8008e123","0dm00009f9730600dge0404000","411793105f032230000","0h63103100000010");
    cuts.AddCutPCMCalo("8008d123","0dm00009f9730600dge0404000","411793105f032230000","0h63103100000010");
  } else if (trainConfig == 2071) { // TOF timing requirement both legs
    cuts.AddCutPCMCalo("80010123","0dm00009f9730700dge0404000","411793105f032230000","0h63103100000010");
    cuts.AddCutPCMCalo("8008e123","0dm00009f9730700dge0404000","411793105f032230000","0h63103100000010");
    cuts.AddCutPCMCalo("8008d123","0dm00009f9730700dge0404000","411793105f032230000","0h63103100000010");
  } else if (trainConfig == 2072) { // TOF timing requirement both legs
    cuts.AddCutPCMCalo("80010123","0dm00009f9730800dge0404000","411793105f032230000","0h63103100000010");
    cuts.AddCutPCMCalo("8008e123","0dm00009f9730800dge0404000","411793105f032230000","0h63103100000010");
    cuts.AddCutPCMCalo("8008d123","0dm00009f9730800dge0404000","411793105f032230000","0h63103100000010");
  } else if (trainConfig == 2073) { // TOF timing requirement both legs
    cuts.AddCutPCMCalo("80010123","0dm00009f9730900dge0404000","411793105f032230000","0h63103100000010");
    cuts.AddCutPCMCalo("8008e123","0dm00009f9730900dge0404000","411793105f032230000","0h63103100000010");
    cuts.AddCutPCMCalo("8008d123","0dm00009f9730900dge0404000","411793105f032230000","0h63103100000010");

  } else if (trainConfig == 2080) { // EMCal-only config for decay gamma MC
    cuts.AddCutPCMCalo("80010103","0dm00009f9730000dge0404000","111113105f032230000","0h63103100000010");
    cuts.AddCutPCMCalo("80085103","0dm00009f9730000dge0404000","111113105f032230000","0h63103100000010");
    cuts.AddCutPCMCalo("80083103","0dm00009f9730000dge0404000","111113105f032230000","0h63103100000010");
  } else if (trainConfig == 2081) { // DCal-only config for decay gamma MC
    cuts.AddCutPCMCalo("80010103","0dm00009f9730000dge0404000","388553105f032230000","0h63103100000010");
    cuts.AddCutPCMCalo("8008b103","0dm00009f9730000dge0404000","388553105f032230000","0h63103100000010");
    cuts.AddCutPCMCalo("80089103","0dm00009f9730000dge0404000","388553105f032230000","0h63103100000010");
  } else if (trainConfig == 2082) { // EMCal+DCal config for decay gamma MC
    cuts.AddCutPCMCalo("80010103","0dm00009f9730000dge0404000","411793105f032230000","0h63103100000010");
    cuts.AddCutPCMCalo("8008e103","0dm00009f9730000dge0404000","411793105f032230000","0h63103100000010");
    cuts.AddCutPCMCalo("8008d103","0dm00009f9730000dge0404000","411793105f032230000","0h63103100000010");

  } else if (trainConfig == 2083) { // EMCal+DCal config for decay gamma MC, T0 based cuts
    cuts.AddCutPCMCalo("80011103","0dm00009f9730000dge0404000","411793105f032230000","0h63103100000010");
    cuts.AddCutPCMCalo("8008g103","0dm00009f9730000dge0404000","411793105f032230000","0h63103100000010");
    cuts.AddCutPCMCalo("8008f103","0dm00009f9730000dge0404000","411793105f032230000","0h63103100000010");

  } else if (trainConfig == 2090) { // EMCal-only config for decay gamma MC with smearing
    cuts.AddCutPCMCalo("80010103","0dm00009f9730000dge0404000","111113105f032230000","0h63103100b00010");
    cuts.AddCutPCMCalo("80085103","0dm00009f9730000dge0404000","111113105f032230000","0h63103100b00010");
    cuts.AddCutPCMCalo("80083103","0dm00009f9730000dge0404000","111113105f032230000","0h63103100b00010");
  } else if (trainConfig == 2091) { // DCal-only config for decay gamma MC with smearing
    cuts.AddCutPCMCalo("80010103","0dm00009f9730000dge0404000","388553105f032230000","0h63103100b00010");
    cuts.AddCutPCMCalo("8008b103","0dm00009f9730000dge0404000","388553105f032230000","0h63103100b00010");
    cuts.AddCutPCMCalo("80089103","0dm00009f9730000dge0404000","388553105f032230000","0h63103100b00010");
  } else if (trainConfig == 2092) { // EMCal+DCal config for decay gamma MC with smearing
    cuts.AddCutPCMCalo("80010103","0dm00009f9730000dge0404000","411793105f032230000","0h63103100b00010");
    cuts.AddCutPCMCalo("8008e103","0dm00009f9730000dge0404000","411793105f032230000","0h63103100b00010");
    cuts.AddCutPCMCalo("8008d103","0dm00009f9730000dge0404000","411793105f032230000","0h63103100b00010");

  } else if (trainConfig == 2730) { // PCM variations
    cuts.AddCutPCMCalo("80010103","00200009f9730000dge0400000","411793105f032230000","0h63103100000010"); //New standard cut for 8TeV analysis for RpA
    cuts.AddCutPCMCalo("80010003","00200009f9730000dge0400000","411793105f032230000","0h63103100000010"); // no SPD pileup cut
    cuts.AddCutPCMCalo("80010103","00100009f9730000dge0400000","411793105f032230000","0h63103100000010"); // R cut 2.8 -180 cm
    cuts.AddCutPCMCalo("80010103","00500009f9730000dge0400000","411793105f032230000","0h63103100000010"); // R cut 10. -180 cm
  } else if (trainConfig == 2731) {
    cuts.AddCutPCMCalo("80010103","00200069f9730000dge0400000","411793105f032230000","0h63103100000010"); // min pT 40 MeV
    cuts.AddCutPCMCalo("80010103","00200049f9730000dge0400000","411793105f032230000","0h63103100000010"); // min pT 75 MeV
    cuts.AddCutPCMCalo("80010103","00200019f9730000dge0400000","411793105f032230000","0h63103100000010"); // min pT 100MeV
  } else if (trainConfig == 2732) {
    cuts.AddCutPCMCalo("80010103","00200068f9730000dge0400000","411793105f032230000","0h63103100000010"); // TPC cluster 35%
    cuts.AddCutPCMCalo("80010103","00200066f9730000dge0400000","411793105f032230000","0h63103100000010"); // TPC cluster 70%
    cuts.AddCutPCMCalo("80010103","00200009f9730000dge0600000","411793105f032230000","0h63103100000010"); // cosPA 0.9
    cuts.AddCutPCMCalo("80010103","00200009f9730000dge0300000","411793105f032230000","0h63103100000010"); // cosPA 0.75
  } else if (trainConfig == 2733) {
    cuts.AddCutPCMCalo("80010103","0020000939730000dge0400000","411793105f032230000","0h63103100000010"); // nsig electron   -4,5
    cuts.AddCutPCMCalo("80010103","0020000969730000dge0400000","411793105f032230000","0h63103100000010"); // nsig electron -2.5,4
    cuts.AddCutPCMCalo("80010103","00200009f5730000dge0400000","411793105f032230000","0h63103100000010"); // nsig pion 2,-10
    cuts.AddCutPCMCalo("80010103","00200009f1730000dge0400000","411793105f032230000","0h63103100000010"); // nsig pion 0,-10
  } else if (trainConfig == 2734) {
    cuts.AddCutPCMCalo("80010103","00200009f9030000dge0400000","411793105f032230000","0h63103100000010"); // pion nsig min mom 0.50 GeV/c
    cuts.AddCutPCMCalo("80010103","00200009f9630000dge0400000","411793105f032230000","0h63103100000010"); // pion nsig min mom 0.25 GeV/c
    cuts.AddCutPCMCalo("80010103","00200009f9760000dge0400000","411793105f032230000","0h63103100000010"); // pion nsig max mom 2.00 GeV/c
    cuts.AddCutPCMCalo("80010103","00200009f9710000dge0400000","411793105f032230000","0h63103100000010"); // pion nsig max mom 5.00 GeV/c
  } else if (trainConfig == 2735) {
    cuts.AddCutPCMCalo("80010103","00200009f9730000age0400000","411793105f032230000","0h63103100000010"); // qT max 0.040, qt pt max 0.11
    cuts.AddCutPCMCalo("80010103","00200009f9730000ege0400000","411793105f032230000","0h63103100000010"); // qT max 0.060, qt pt max 0.14
    cuts.AddCutPCMCalo("80010103","00200009f9730000fge0400000","411793105f032230000","0h63103100000010"); // qT max 0.070, qt pt max 0.16
  } else if (trainConfig == 2736) {
    cuts.AddCutPCMCalo("80010103","00200009f9730000d1e0400000","411793105f032230000","0h63103100000010"); // chi2 50 no chi2 dep.
    cuts.AddCutPCMCalo("80010103","00200009f9730000dfe0400000","411793105f032230000","0h63103100000010"); // chi2 50 chi2 dep -0.065
    cuts.AddCutPCMCalo("80010103","00200009f9730000dhe0400000","411793105f032230000","0h63103100000010"); // chi2 50 chi2 dep -0.050
    cuts.AddCutPCMCalo("80010103","00200009f9730000dge0404000","411793105f032230000","0h63103100000010"); // reject close v0
    cuts.AddCutPCMCalo("80010103","00200009f9730000dge0406000","411793105f032230000","0h63103100000010"); // double count with open angle 0.04
  } else if (trainConfig == 2737) {
    cuts.AddCutPCMCalo("80010103","00200009f9730000dgd0400000","411793105f032230000","0h63103100000010"); // Psi pair 0.15 dep
    cuts.AddCutPCMCalo("80010103","00200009f9730000dgf0400000","411793105f032230000","0h63103100000010"); // Psi pair 0.20 dep
    cuts.AddCutPCMCalo("80010103","00200009f9730000dgg0400000","411793105f032230000","0h63103100000010"); // Psi pair 0.30 dep
  } else if (trainConfig == 2738) {
    cuts.AddCutPCMCalo("80010103","00200009f9730000dge0400000","411793105f032230000","0263103100000010"); // variation BG scheme track mult
    cuts.AddCutPCMCalo("80010103","00200009f9730000dge0400000","411793105f032230000","0h63107100000010"); // alpha meson 0.85
    cuts.AddCutPCMCalo("80010103","00200009f9730000dge0400000","411793105f032230000","0h63105100000010"); // alpha meson 0.75
    cuts.AddCutPCMCalo("80010103","00200009227300008250404000","411793105f032230000","0h63103100000010"); // old cuts (run1)
  } else if (trainConfig == 2739) {
    cuts.AddCutPCMCalo("80010203","00200009f9730000dge0400000","411793105f032230000","0h63103100000010"); //same as std + maximum past future rejection
    cuts.AddCutPCMCalo("80010503","00200009f9730000dge0400000","411793105f032230000","0h63103100000010"); //same as std + medium past future rejection
  } else if (trainConfig == 2740) { // CALO variations
    cuts.AddCutPCMCalo("80010103","00200009f9730000dge0400000","411793105f022230000","0h63103100000010"); // 0.6 GeV/c
    cuts.AddCutPCMCalo("80010103","00200009f9730000dge0400000","411793105f042230000","0h63103100000010"); // 0.8 GeV/c
  } else if (trainConfig == 2741) {
    cuts.AddCutPCMCalo("80010103","00200009f9730000dge0400000","411793105f031230000","0h63103100000010"); // n cells >= 1
    cuts.AddCutPCMCalo("80010103","00200009f9730000dge0400000","411793105f033230000","0h63103100000010"); // n cells >= 3
    cuts.AddCutPCMCalo("80010103","00200009f9730000dge0400000","411793105f032200000","0h63103100000010"); // no max M02 cut
    cuts.AddCutPCMCalo("80010103","00200009f9730000dge0400000","411793105f032250000","0h63103100000010"); // M02 < 0.3
    cuts.AddCutPCMCalo("80010103","00200009f9730000dge0400000","411793105f0322k0000","0h63103100000010"); // M02, pT-dep
  } else if (trainConfig == 2742) {
    cuts.AddCutPCMCalo("80010103","00200009f9730000dge0400000","411793105e032230000","0h63103100000010"); // TM var e/p 2.0
    cuts.AddCutPCMCalo("80010103","00200009f9730000dge0400000","411793105g032230000","0h63103100000010"); // TM var e/p 1.5
    cuts.AddCutPCMCalo("80010103","00200009f9730000dge0400000","411793105h032230000","0h63103100000010"); // TM var e/p 1.25
    cuts.AddCutPCMCalo("80010103","00200009f9730000dge0400000","4117931057032230000","0h63103100000010"); // TM var no veto
  } else if (trainConfig == 2743) {
    cuts.AddCutPCMCalo("80010103","00200009f9730000dge0400000","411793205f032230000","0h63103100000010"); // NL 32
    cuts.AddCutPCMCalo("80010103","00200009f9730000dge0400000","411793305f032230000","0h63103100000010"); // NL 33
    cuts.AddCutPCMCalo("80010103","00200009f9730000dge0400000","411793405f032230000","0h63103100000010"); // NL 34
  } else if (trainConfig == 2744) {
    cuts.AddCutPCMCalo("80010103","00200009f9730000dge0400000","411793106f032230000","0h63103100000010"); // 30/35ns
    cuts.AddCutPCMCalo("80010103","00200009f9730000dge0400000","411793104f032230000","0h63103100000010"); // 100ns
    cuts.AddCutPCMCalo("80010103","00200009f9730000dge0400000","411793107f032230000","0h63103100000010"); // 30ns

  } else if (trainConfig == 2750) { // PCM variations
    cuts.AddCutPCMCalo("8008e103","00200009f9730000dge0400000","411793105f032230000","0h63103100000010"); //New standard cut for 8TeV analysis for RpA
    cuts.AddCutPCMCalo("8008e003","00200009f9730000dge0400000","411793105f032230000","0h63103100000010"); // no SPD pileup cut
    cuts.AddCutPCMCalo("8008e103","00100009f9730000dge0400000","411793105f032230000","0h63103100000010"); // R cut 2.8 -180 cm
    cuts.AddCutPCMCalo("8008e103","00500009f9730000dge0400000","411793105f032230000","0h63103100000010"); // R cut 10. -180 cm
  } else if (trainConfig == 2751) {
    cuts.AddCutPCMCalo("8008e103","00200069f9730000dge0400000","411793105f032230000","0h63103100000010"); // min pT 40 MeV
    cuts.AddCutPCMCalo("8008e103","00200049f9730000dge0400000","411793105f032230000","0h63103100000010"); // min pT 75 MeV
    cuts.AddCutPCMCalo("8008e103","00200019f9730000dge0400000","411793105f032230000","0h63103100000010"); // min pT 100MeV
  } else if (trainConfig == 2752) {
    cuts.AddCutPCMCalo("8008e103","00200068f9730000dge0400000","411793105f032230000","0h63103100000010"); // TPC cluster 35%
    cuts.AddCutPCMCalo("8008e103","00200066f9730000dge0400000","411793105f032230000","0h63103100000010"); // TPC cluster 70%
    cuts.AddCutPCMCalo("8008e103","00200009f9730000dge0600000","411793105f032230000","0h63103100000010"); // cosPA 0.9
    cuts.AddCutPCMCalo("8008e103","00200009f9730000dge0300000","411793105f032230000","0h63103100000010"); // cosPA 0.75
  } else if (trainConfig == 2753) {
    cuts.AddCutPCMCalo("8008e103","0020000939730000dge0400000","411793105f032230000","0h63103100000010"); // nsig electron   -4,5
    cuts.AddCutPCMCalo("8008e103","0020000969730000dge0400000","411793105f032230000","0h63103100000010"); // nsig electron -2.5,4
    cuts.AddCutPCMCalo("8008e103","00200009f5730000dge0400000","411793105f032230000","0h63103100000010"); // nsig pion 2,-10
    cuts.AddCutPCMCalo("8008e103","00200009f1730000dge0400000","411793105f032230000","0h63103100000010"); // nsig pion 0,-10
  } else if (trainConfig == 2754) {
    cuts.AddCutPCMCalo("8008e103","00200009f9030000dge0400000","411793105f032230000","0h63103100000010"); // pion nsig min mom 0.50 GeV/c
    cuts.AddCutPCMCalo("8008e103","00200009f9630000dge0400000","411793105f032230000","0h63103100000010"); // pion nsig min mom 0.25 GeV/c
    cuts.AddCutPCMCalo("8008e103","00200009f9760000dge0400000","411793105f032230000","0h63103100000010"); // pion nsig max mom 2.00 GeV/c
    cuts.AddCutPCMCalo("8008e103","00200009f9710000dge0400000","411793105f032230000","0h63103100000010"); // pion nsig max mom 5.00 GeV/c
  } else if (trainConfig == 2755) {
    cuts.AddCutPCMCalo("8008e103","00200009f9730000age0400000","411793105f032230000","0h63103100000010"); // qT max 0.040, qt pt max 0.11
    cuts.AddCutPCMCalo("8008e103","00200009f9730000ege0400000","411793105f032230000","0h63103100000010"); // qT max 0.060, qt pt max 0.14
    cuts.AddCutPCMCalo("8008e103","00200009f9730000fge0400000","411793105f032230000","0h63103100000010"); // qT max 0.070, qt pt max 0.16
  } else if (trainConfig == 2756) {
    cuts.AddCutPCMCalo("8008e103","00200009f9730000d1e0400000","411793105f032230000","0h63103100000010"); // chi2 50 no chi2 dep.
    cuts.AddCutPCMCalo("8008e103","00200009f9730000dfe0400000","411793105f032230000","0h63103100000010"); // chi2 50 chi2 dep -0.065
    cuts.AddCutPCMCalo("8008e103","00200009f9730000dhe0400000","411793105f032230000","0h63103100000010"); // chi2 50 chi2 dep -0.050
    cuts.AddCutPCMCalo("8008e103","00200009f9730000dge0404000","411793105f032230000","0h63103100000010"); // reject close v0
    cuts.AddCutPCMCalo("8008e103","00200009f9730000dge0406000","411793105f032230000","0h63103100000010"); // double count with open angle 0.04
  } else if (trainConfig == 2757) {
    cuts.AddCutPCMCalo("8008e103","00200009f9730000dgd0400000","411793105f032230000","0h63103100000010"); // Psi pair 0.15 dep
    cuts.AddCutPCMCalo("8008e103","00200009f9730000dgf0400000","411793105f032230000","0h63103100000010"); // Psi pair 0.20 dep
    cuts.AddCutPCMCalo("8008e103","00200009f9730000dgg0400000","411793105f032230000","0h63103100000010"); // Psi pair 0.30 dep
  } else if (trainConfig == 2758) {
    cuts.AddCutPCMCalo("8008e103","00200009f9730000dge0400000","411793105f032230000","0263103100000010"); // variation BG scheme track mult
    cuts.AddCutPCMCalo("8008e103","00200009f9730000dge0400000","411793105f032230000","0h63107100000010"); // alpha meson 0.85
    cuts.AddCutPCMCalo("8008e103","00200009f9730000dge0400000","411793105f032230000","0h63105100000010"); // alpha meson 0.75
    cuts.AddCutPCMCalo("8008e103","00200009227300008250404000","411793105f032230000","0h63103100000010"); // old cuts (run1)
  } else if (trainConfig == 2759) {
    cuts.AddCutPCMCalo("8008e203","00200009f9730000dge0400000","411793105f032230000","0h63103100000010"); //same as std + maximum past future rejection
    cuts.AddCutPCMCalo("8008e503","00200009f9730000dge0400000","411793105f032230000","0h63103100000010"); //same as std + medium past future rejection
  } else if (trainConfig == 2760) { // CALO variations
    cuts.AddCutPCMCalo("8008e103","00200009f9730000dge0400000","411793105f022230000","0h63103100000010"); // 0.6 GeV/c
    cuts.AddCutPCMCalo("8008e103","00200009f9730000dge0400000","411793105f042230000","0h63103100000010"); // 0.8 GeV/c
  } else if (trainConfig == 2761) {
    cuts.AddCutPCMCalo("8008e103","00200009f9730000dge0400000","411793105f031230000","0h63103100000010"); // n cells >= 1
    cuts.AddCutPCMCalo("8008e103","00200009f9730000dge0400000","411793105f033230000","0h63103100000010"); // n cells >= 3
    cuts.AddCutPCMCalo("8008e103","00200009f9730000dge0400000","411793105f032200000","0h63103100000010"); // no max M02 cut
    cuts.AddCutPCMCalo("8008e103","00200009f9730000dge0400000","411793105f032250000","0h63103100000010"); // M02 < 0.3
    cuts.AddCutPCMCalo("8008e103","00200009f9730000dge0400000","411793105f0322k0000","0h63103100000010"); // M02, pT-dep
  } else if (trainConfig == 2762) {
    cuts.AddCutPCMCalo("8008e103","00200009f9730000dge0400000","411793105e032230000","0h63103100000010"); // TM var e/p 2.0
    cuts.AddCutPCMCalo("8008e103","00200009f9730000dge0400000","411793105g032230000","0h63103100000010"); // TM var e/p 1.5
    cuts.AddCutPCMCalo("8008e103","00200009f9730000dge0400000","411793105h032230000","0h63103100000010"); // TM var e/p 1.25
    cuts.AddCutPCMCalo("8008e103","00200009f9730000dge0400000","4117931057032230000","0h63103100000010"); // TM var no veto
  } else if (trainConfig == 2763) {
    cuts.AddCutPCMCalo("8008e103","00200009f9730000dge0400000","411793205f032230000","0h63103100000010"); // NL 32
    cuts.AddCutPCMCalo("8008e103","00200009f9730000dge0400000","411793305f032230000","0h63103100000010"); // NL 33
    cuts.AddCutPCMCalo("8008e103","00200009f9730000dge0400000","411793405f032230000","0h63103100000010"); // NL 34
  } else if (trainConfig == 2764) {
    cuts.AddCutPCMCalo("8008e103","00200009f9730000dge0400000","411793106f032230000","0h63103100000010"); // 30/35ns
    cuts.AddCutPCMCalo("8008e103","00200009f9730000dge0400000","411793104f032230000","0h63103100000010"); // 100ns
    cuts.AddCutPCMCalo("8008e103","00200009f9730000dge0400000","411793107f032230000","0h63103100000010"); // 30ns

  } else if (trainConfig == 2770) { // PCM variations
    cuts.AddCutPCMCalo("8008d103","00200009f9730000dge0400000","411793105f032230000","0h63103100000010"); //New standard cut for 8TeV analysis for RpA
    cuts.AddCutPCMCalo("8008d003","00200009f9730000dge0400000","411793105f032230000","0h63103100000010"); // no SPD pileup cut
    cuts.AddCutPCMCalo("8008d103","00100009f9730000dge0400000","411793105f032230000","0h63103100000010"); // R cut 2.8 -180 cm
    cuts.AddCutPCMCalo("8008d103","00500009f9730000dge0400000","411793105f032230000","0h63103100000010"); // R cut 10. -180 cm
  } else if (trainConfig == 2771) {
    cuts.AddCutPCMCalo("8008d103","00200069f9730000dge0400000","411793105f032230000","0h63103100000010"); // min pT 40 MeV
    cuts.AddCutPCMCalo("8008d103","00200049f9730000dge0400000","411793105f032230000","0h63103100000010"); // min pT 75 MeV
    cuts.AddCutPCMCalo("8008d103","00200019f9730000dge0400000","411793105f032230000","0h63103100000010"); // min pT 100MeV
  } else if (trainConfig == 2772) {
    cuts.AddCutPCMCalo("8008d103","00200068f9730000dge0400000","411793105f032230000","0h63103100000010"); // TPC cluster 35%
    cuts.AddCutPCMCalo("8008d103","00200066f9730000dge0400000","411793105f032230000","0h63103100000010"); // TPC cluster 70%
    cuts.AddCutPCMCalo("8008d103","00200009f9730000dge0600000","411793105f032230000","0h63103100000010"); // cosPA 0.9
    cuts.AddCutPCMCalo("8008d103","00200009f9730000dge0300000","411793105f032230000","0h63103100000010"); // cosPA 0.75
  } else if (trainConfig == 2773) {
    cuts.AddCutPCMCalo("8008d103","0020000939730000dge0400000","411793105f032230000","0h63103100000010"); // nsig electron   -4,5
    cuts.AddCutPCMCalo("8008d103","0020000969730000dge0400000","411793105f032230000","0h63103100000010"); // nsig electron -2.5,4
    cuts.AddCutPCMCalo("8008d103","00200009f5730000dge0400000","411793105f032230000","0h63103100000010"); // nsig pion 2,-10
    cuts.AddCutPCMCalo("8008d103","00200009f1730000dge0400000","411793105f032230000","0h63103100000010"); // nsig pion 0,-10
  } else if (trainConfig == 2774) {
    cuts.AddCutPCMCalo("8008d103","00200009f9030000dge0400000","411793105f032230000","0h63103100000010"); // pion nsig min mom 0.50 GeV/c
    cuts.AddCutPCMCalo("8008d103","00200009f9630000dge0400000","411793105f032230000","0h63103100000010"); // pion nsig min mom 0.25 GeV/c
    cuts.AddCutPCMCalo("8008d103","00200009f9760000dge0400000","411793105f032230000","0h63103100000010"); // pion nsig max mom 2.00 GeV/c
    cuts.AddCutPCMCalo("8008d103","00200009f9710000dge0400000","411793105f032230000","0h63103100000010"); // pion nsig max mom 5.00 GeV/c
  } else if (trainConfig == 2775) {
    cuts.AddCutPCMCalo("8008d103","00200009f9730000age0400000","411793105f032230000","0h63103100000010"); // qT max 0.040, qt pt max 0.11
    cuts.AddCutPCMCalo("8008d103","00200009f9730000ege0400000","411793105f032230000","0h63103100000010"); // qT max 0.060, qt pt max 0.14
    cuts.AddCutPCMCalo("8008d103","00200009f9730000fge0400000","411793105f032230000","0h63103100000010"); // qT max 0.070, qt pt max 0.16
  } else if (trainConfig == 2776) {
    cuts.AddCutPCMCalo("8008d103","00200009f9730000d1e0400000","411793105f032230000","0h63103100000010"); // chi2 50 no chi2 dep.
    cuts.AddCutPCMCalo("8008d103","00200009f9730000dfe0400000","411793105f032230000","0h63103100000010"); // chi2 50 chi2 dep -0.065
    cuts.AddCutPCMCalo("8008d103","00200009f9730000dhe0400000","411793105f032230000","0h63103100000010"); // chi2 50 chi2 dep -0.050
    cuts.AddCutPCMCalo("8008d103","00200009f9730000dge0404000","411793105f032230000","0h63103100000010"); // reject close v0
    cuts.AddCutPCMCalo("8008d103","00200009f9730000dge0406000","411793105f032230000","0h63103100000010"); // double count with open angle 0.04
  } else if (trainConfig == 2777) {
    cuts.AddCutPCMCalo("8008d103","00200009f9730000dgd0400000","411793105f032230000","0h63103100000010"); // Psi pair 0.15 dep
    cuts.AddCutPCMCalo("8008d103","00200009f9730000dgf0400000","411793105f032230000","0h63103100000010"); // Psi pair 0.20 dep
    cuts.AddCutPCMCalo("8008d103","00200009f9730000dgg0400000","411793105f032230000","0h63103100000010"); // Psi pair 0.30 dep
  } else if (trainConfig == 2778) {
    cuts.AddCutPCMCalo("8008d103","00200009f9730000dge0400000","411793105f032230000","0263103100000010"); // variation BG scheme track mult
    cuts.AddCutPCMCalo("8008d103","00200009f9730000dge0400000","411793105f032230000","0h63107100000010"); // alpha meson 0.85
    cuts.AddCutPCMCalo("8008d103","00200009f9730000dge0400000","411793105f032230000","0h63105100000010"); // alpha meson 0.75
    cuts.AddCutPCMCalo("8008d103","00200009227300008250404000","411793105f032230000","0h63103100000010"); // old cuts (run1)
  } else if (trainConfig == 2779) {
    cuts.AddCutPCMCalo("8008d203","00200009f9730000dge0400000","411793105f032230000","0h63103100000010"); //same as std + maximum past future rejection
    cuts.AddCutPCMCalo("8008d503","00200009f9730000dge0400000","411793105f032230000","0h63103100000010"); //same as std + medium past future rejection
  } else if (trainConfig == 2780) { // CALO variations
    cuts.AddCutPCMCalo("8008d103","00200009f9730000dge0400000","411793105f022230000","0h63103100000010"); // 0.6 GeV/c
    cuts.AddCutPCMCalo("8008d103","00200009f9730000dge0400000","411793105f042230000","0h63103100000010"); // 0.8 GeV/c
  } else if (trainConfig == 2781) {
    cuts.AddCutPCMCalo("8008d103","00200009f9730000dge0400000","411793105f031230000","0h63103100000010"); // n cells >= 1
    cuts.AddCutPCMCalo("8008d103","00200009f9730000dge0400000","411793105f033230000","0h63103100000010"); // n cells >= 3
    cuts.AddCutPCMCalo("8008d103","00200009f9730000dge0400000","411793105f032200000","0h63103100000010"); // no max M02 cut
    cuts.AddCutPCMCalo("8008d103","00200009f9730000dge0400000","411793105f032250000","0h63103100000010"); // M02 < 0.3
    cuts.AddCutPCMCalo("8008d103","00200009f9730000dge0400000","411793105f0322k0000","0h63103100000010"); // M02, pT-dep
  } else if (trainConfig == 2782) {
    cuts.AddCutPCMCalo("8008d103","00200009f9730000dge0400000","411793105e032230000","0h63103100000010"); // TM var e/p 2.0
    cuts.AddCutPCMCalo("8008d103","00200009f9730000dge0400000","411793105g032230000","0h63103100000010"); // TM var e/p 1.5
    cuts.AddCutPCMCalo("8008d103","00200009f9730000dge0400000","411793105h032230000","0h63103100000010"); // TM var e/p 1.25
    cuts.AddCutPCMCalo("8008d103","00200009f9730000dge0400000","4117931057032230000","0h63103100000010"); // TM var no veto
  } else if (trainConfig == 2783) {
    cuts.AddCutPCMCalo("8008d103","00200009f9730000dge0400000","411793205f032230000","0h63103100000010"); // NL 32
    cuts.AddCutPCMCalo("8008d103","00200009f9730000dge0400000","411793305f032230000","0h63103100000010"); // NL 33
    cuts.AddCutPCMCalo("8008d103","00200009f9730000dge0400000","411793405f032230000","0h63103100000010"); // NL 34
  } else if (trainConfig == 2784) {
    cuts.AddCutPCMCalo("8008d103","00200009f9730000dge0400000","411793106f032230000","0h63103100000010"); // 30/35ns
    cuts.AddCutPCMCalo("8008d103","00200009f9730000dge0400000","411793104f032230000","0h63103100000010"); // 100ns
    cuts.AddCutPCMCalo("8008d103","00200009f9730000dge0400000","411793107f032230000","0h63103100000010"); // 30ns


  // pPb 8 TeV PHOS new default with timing effi
  } else if (trainConfig == 3000){ // PHOS  INT7
    cuts.AddCutPCMCalo("80010103","0dm00009f9730000dge0404000","24466000ha012200000","0h63103100000010"); // 0-100% without NL
  } else if (trainConfig == 3001){ // PHOS  PHI7
    cuts.AddCutPCMCalo("80062103","0dm00009f9730000dge0404000","24466000ha012200000","0h63103100000010"); // 0-100% without NL
  } else if (trainConfig == 3002){
    cuts.AddCutPCMCalo("80010103","0dm00009f9730000dge0404000","24466590ha012200000","0h63103100000010"); // 59 NL
  } else if (trainConfig == 3003){
    cuts.AddCutPCMCalo("80062103","0dm00009f9730000dge0404000","24466590ha012200000","0h63103100000010"); // 59 NL
  } else if (trainConfig == 3004){
    cuts.AddCutPCMCalo("80010103","0dm00009f9730000dge0404000","24466690ha012200000","0h63103100000010"); // 69 NL
  } else if (trainConfig == 3005){
    cuts.AddCutPCMCalo("80062103","0dm00009f9730000dge0404000","24466690ha012200000","0h63103100000010"); // 69 NL

  } else {
    Error(Form("GammaConvCalo_%i",trainConfig), "wrong trainConfig variable no cuts have been specified for the configuration");
    return;
  }

  if(!cuts.AreValid()){
    cout << "\n\n****************************************************" << endl;
    cout << "ERROR: No valid cuts stored in CutHandlerConvCalo! Returning..." << endl;
    cout << "****************************************************\n\n" << endl;
    return;
  }

  TList *EventCutList   = new TList();
  TList *ConvCutList    = new TList();
  TList *ClusterCutList = new TList();
  TList *MesonCutList   = new TList();

  Int_t numberOfCuts = cuts.GetNCuts();

  TList *HeaderList     = new TList();
  if (doWeightingPart==1) {
    TObjString *Header1 = new TObjString("pi0_1");
    HeaderList->Add(Header1);
  }
  if (doWeightingPart==2){
    TObjString *Header3 = new TObjString("eta_2");
    HeaderList->Add(Header3);
  }
  if (doWeightingPart==3) {
    TObjString *Header1 = new TObjString("pi0_1");
    HeaderList->Add(Header1);
    TObjString *Header3 = new TObjString("eta_2");
    HeaderList->Add(Header3);
  }
  if (periodNameV0Reader.Contains("LHC18b9")||periodNameV0Reader.Contains("LHC17g8")){
    TObjString *HeaderPMB = new TObjString("EPOSLHC_0");
    TObjString *HeaderP8J = new TObjString("Pythia8Jets_1");
    if (doWeightingPart==4) { // all headers
      HeaderList->Add(HeaderPMB);
      HeaderList->Add(HeaderP8J);
    } else if (doWeightingPart==5) { // only MB header
      HeaderList->Add(HeaderPMB);
    } else { // only JJ header
      HeaderList->Add(HeaderP8J);
    }
  } else if (periodNameV0Reader.Contains("LHC17g6a2") || periodNameV0Reader.Contains("LHC17g6a3") ){
    TObjString *HeaderPMB = new TObjString("Dpmjet_0");
    TObjString *HeaderP8J = new TObjString("Pythia8JetsGammaTrg_1");
    if (doWeightingPart==4) { // all headers
      HeaderList->Add(HeaderPMB);
      HeaderList->Add(HeaderP8J);
    } else if (doWeightingPart==5) { // only MB header
      HeaderList->Add(HeaderPMB);
    } else { // only JJ header
      HeaderList->Add(HeaderP8J);
    }
  }

  EventCutList->SetOwner(kTRUE);
  AliConvEventCuts **analysisEventCuts        = new AliConvEventCuts*[numberOfCuts];
  ConvCutList->SetOwner(kTRUE);
  AliConversionPhotonCuts **analysisCuts      = new AliConversionPhotonCuts*[numberOfCuts];
  ClusterCutList->SetOwner(kTRUE);
  AliCaloPhotonCuts **analysisClusterCuts     = new AliCaloPhotonCuts*[numberOfCuts];
  MesonCutList->SetOwner(kTRUE);
  AliConversionMesonCuts **analysisMesonCuts  = new AliConversionMesonCuts*[numberOfCuts];

  Bool_t initializedMatBudWeigths_existing    = kFALSE;

  for(Int_t i = 0; i<numberOfCuts; i++){
    //create AliCaloTrackMatcher instance, if there is none present
    TString caloCutPos = cuts.GetClusterCut(i);
    caloCutPos.Resize(1);
    TString TrackMatcherName = Form("CaloTrackMatcher_%s_%i",caloCutPos.Data(),trackMatcherRunningMode);
    if(corrTaskSetting.CompareTo("")){
      TrackMatcherName = TrackMatcherName+"_"+corrTaskSetting.Data();
      cout << "Using separate track matcher for correction framework setting: " << TrackMatcherName.Data() << endl;
    }
    if( !(AliCaloTrackMatcher*)mgr->GetTask(TrackMatcherName.Data()) ){
      AliCaloTrackMatcher* fTrackMatcher = new AliCaloTrackMatcher(TrackMatcherName.Data(),caloCutPos.Atoi(),trackMatcherRunningMode);
      fTrackMatcher->SetV0ReaderName(V0ReaderName);
      fTrackMatcher->SetCorrectionTaskSetting(corrTaskSetting);
      mgr->AddTask(fTrackMatcher);
      mgr->ConnectInput(fTrackMatcher,0,cinput);
    }

    analysisEventCuts[i] = new AliConvEventCuts();

    TString dataInputMultHisto  = "";
    TString mcInputMultHisto    = "";
    TString triggerString       = cuts.GetEventCut(i);
    triggerString               = triggerString(3,2);

    dataInputMultHisto          = Form("%s_%s", periodNameAnchor.Data(), triggerString.Data());
    mcInputMultHisto            = Form("%s_%s", periodNameV0Reader.Data(), triggerString.Data());

    if (enableMultiplicityWeighting){
      cout << "enabling mult weighting" << endl;
      analysisEventCuts[i]->SetUseWeightMultiplicityFromFile( kTRUE, fileNameMultWeights, dataInputMultHisto, mcInputMultHisto );
    }

    analysisEventCuts[i]->SetTriggerMimicking(enableTriggerMimicking);
    if(fileNameCustomTriggerMimicOADB.CompareTo("") != 0)
      analysisEventCuts[i]->SetCustomTriggerMimicOADBFile(fileNameCustomTriggerMimicOADB);
    analysisEventCuts[i]->SetTriggerOverlapRejecion(enableTriggerOverlapRej);
    if(fMinPtHardSet)
      analysisEventCuts[i]->SetMinFacPtHard(minFacPtHard);
    if(fMaxPtHardSet)
      analysisEventCuts[i]->SetMaxFacPtHard(maxFacPtHard);
    if(fSingleMaxPtHardSet)
      analysisEventCuts[i]->SetMaxFacPtHardSingleParticle(maxFacPtHardSingle);
    if(fJetFinderUsage)
      analysisEventCuts[i]->SetUseJetFinderForOutliers(kTRUE);
    if(fUsePtHardFromFile)
      analysisEventCuts[i]->SetUsePtHardBinFromFile(kTRUE);
    if(fUseAddOutlierRej)
      analysisEventCuts[i]->SetUseAdditionalOutlierRejection(kTRUE);
    analysisEventCuts[i]->SetV0ReaderName(V0ReaderName);
    analysisEventCuts[i]->SetCorrectionTaskSetting(corrTaskSetting);
    if (periodNameV0Reader.CompareTo("") != 0) analysisEventCuts[i]->SetPeriodEnum(periodNameV0Reader);
    if (enableLightOutput > 0) analysisEventCuts[i]->SetLightOutput(kTRUE);
    analysisEventCuts[i]->InitializeCutsFromCutString((cuts.GetEventCut(i)).Data());
    EventCutList->Add(analysisEventCuts[i]);


    analysisEventCuts[i]->SetFillCutHistograms("",kFALSE);

    analysisCuts[i] = new AliConversionPhotonCuts();

    if (enableMatBudWeightsPi0 > 0){
      if (isMC > 0){
        if (analysisCuts[i]->InitializeMaterialBudgetWeights(enableMatBudWeightsPi0,fileNameMatBudWeights)){
          initializedMatBudWeigths_existing = kTRUE;}
        else {cout << "ERROR The initialization of the materialBudgetWeights did not work out." << endl;}
      }
      else {cout << "ERROR 'enableMatBudWeightsPi0'-flag was set > 0 even though this is not a MC task. It was automatically reset to 0." << endl;}
    }

    analysisCuts[i]->SetV0ReaderName(V0ReaderName);
    if (enableElecDeDxPostCalibration){
      if (isMC == 0){
        if(fileNamedEdxPostCalib.CompareTo("") != 0){
          analysisCuts[i]->SetElecDeDxPostCalibrationCustomFile(fileNamedEdxPostCalib);
          cout << "Setting custom dEdx recalibration file: " << fileNamedEdxPostCalib.Data() << endl;
        }
        analysisCuts[i]->SetDoElecDeDxPostCalibration(enableElecDeDxPostCalibration);
        cout << "Enabled TPC dEdx recalibration." << endl;
      } else{
        cout << "ERROR enableElecDeDxPostCalibration set to True even if MC file. Automatically reset to 0"<< endl;
        enableElecDeDxPostCalibration=kFALSE;
        analysisCuts[i]->SetDoElecDeDxPostCalibration(kFALSE);
      }
    }
    if (enableLightOutput > 0) analysisCuts[i]->SetLightOutput(kTRUE);
    analysisCuts[i]->InitializeCutsFromCutString((cuts.GetPhotonCut(i)).Data());
    analysisCuts[i]->SetIsHeavyIon(isHeavyIon);
    ConvCutList->Add(analysisCuts[i]);
    analysisCuts[i]->SetFillCutHistograms("",kFALSE);

    analysisClusterCuts[i] = new AliCaloPhotonCuts(isMC);
    analysisClusterCuts[i]->SetHistoToModifyAcceptance(histoAcc);
    analysisClusterCuts[i]->SetV0ReaderName(V0ReaderName);
    analysisClusterCuts[i]->SetCorrectionTaskSetting(corrTaskSetting);
    analysisClusterCuts[i]->SetCaloTrackMatcherName(TrackMatcherName);
    if (enableLightOutput > 0) analysisClusterCuts[i]->SetLightOutput(kTRUE);
    analysisClusterCuts[i]->InitializeCutsFromCutString((cuts.GetClusterCut(i)).Data());
    ClusterCutList->Add(analysisClusterCuts[i]);
    analysisClusterCuts[i]->SetExtendedMatchAndQA(enableExtMatchAndQA);
    analysisClusterCuts[i]->SetFillCutHistograms("");

    analysisMesonCuts[i] = new AliConversionMesonCuts();
    if (enableLightOutput > 0) analysisMesonCuts[i]->SetLightOutput(kTRUE);
    analysisMesonCuts[i]->SetRunningMode(2);
    analysisMesonCuts[i]->InitializeCutsFromCutString((cuts.GetMesonCut(i)).Data());
    MesonCutList->Add(analysisMesonCuts[i]);
    analysisMesonCuts[i]->SetFillCutHistograms("");
    analysisEventCuts[i]->SetAcceptedHeader(HeaderList);
  }

  task->SetEventCutList(numberOfCuts,EventCutList);
  task->SetConversionCutList(numberOfCuts,ConvCutList);
  task->SetCaloCutList(numberOfCuts,ClusterCutList);
  task->SetMesonCutList(numberOfCuts,MesonCutList);
  task->SetMoveParticleAccordingToVertex(kTRUE);
  task->SetCorrectionTaskSetting(corrTaskSetting);
  task->SetDoMesonAnalysis(kTRUE);
  task->SetDoMesonQA(enableQAMesonTask); //Attention new switch for Pi0 QA
  task->SetDoPhotonQA(enableQAPhotonTask);  //Attention new switch small for Photon QA
  task->SetDoClusterQA(1);  //Attention new switch small for Cluster QA
  task->SetUseTHnSparse(enableTHnSparse);
  task->SetDoTreeInvMassShowerShape(doTreeClusterShowerShape);
  task->SetEnableSortingOfMCClusLabels(enableSortingMCLabels);
  if(enableExtMatchAndQA > 1){ task->SetPlotHistsExtQA(kTRUE);}

  if (initializedMatBudWeigths_existing) {
      task->SetDoMaterialBudgetWeightingOfGammasForTrueMesons(kTRUE);
  }


  //connect containers
  AliAnalysisDataContainer *coutput =
    mgr->CreateContainer(!(corrTaskSetting.CompareTo("")) ? Form("GammaConvCalo_%i",trainConfig) : Form("GammaConvCalo_%i_%s",trainConfig,corrTaskSetting.Data()), TList::Class(),
              AliAnalysisManager::kOutputContainer,Form("GammaConvCalo_%i.root",trainConfig) );

  mgr->AddTask(task);
  mgr->ConnectInput(task,0,cinput);
  mgr->ConnectOutput(task,1,coutput);
  Int_t nContainer = 2;
  for(Int_t i = 0; i<numberOfCuts; i++){
    if(enableQAPhotonTask>1){
      mgr->ConnectOutput(task,nContainer,mgr->CreateContainer(Form("%s_%s_%s_%s Photon DCA tree",(cuts.GetEventCut(i)).Data(),(cuts.GetPhotonCut(i)).Data(),(cuts.GetClusterCut(i)).Data(),(cuts.GetMesonCut(i)).Data()), TTree::Class(), AliAnalysisManager::kOutputContainer, Form("GammaConvCalo_%i.root",trainConfig)) );
      nContainer++;
    }
    if(enableQAMesonTask>1){
	    mgr->ConnectOutput(task,nContainer,mgr->CreateContainer(Form("%s_%s_%s_%s Meson DCA tree",(cuts.GetEventCut(i)).Data(),(cuts.GetPhotonCut(i)).Data(),(cuts.GetClusterCut(i)).Data(),(cuts.GetMesonCut(i)).Data()), TTree::Class(), AliAnalysisManager::kOutputContainer, Form("GammaConvCalo_%i.root",trainConfig)) );
      nContainer++;
    }
  }

  return;

}
