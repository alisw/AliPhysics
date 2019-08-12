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
  Int_t     enableExtMatchAndQA           = 0,                            // disabled (0), extMatch (1), extQA_noCellQA (2), extMatch+extQA_noCellQA (3), extQA+cellQA (4), extMatch+extQA+cellQA (5)
  Int_t     enableLightOutput             = 0,   // switch to run light output (only essential histograms for afterburner)
  Bool_t    enableTHnSparse               = kFALSE,   // switch on THNsparse
  Bool_t    enableTriggerMimicking        = kFALSE,   // enable trigger mimicking
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
  if(trainConfig >= 520 && trainConfig < 530) task->SetDoHBTHistoOutput(kTRUE);

  // cluster cuts
  // 0 "ClusterType",  1 "EtaMin", 2 "EtaMax", 3 "PhiMin", 4 "PhiMax", 5 "DistanceToBadChannel", 6 "Timing", 7 "TrackMatching", 8 "ExoticCell",
  // 9 "MinEnergy", 10 "MinNCells", 11 "MinM02", 12 "MaxM02", 13 "MinMaxM20", 14 "RecConv", 15 "MaximumDispersion", 16 "NLM"

  //************************************************ EMCAL clusters **********************************************************
  //---------------------------------------------------------------------------------------------
  // no non linearity cuts
  //---------------------------------------------------------------------------------------------
  if (trainConfig == 1){ // min energy = 0.3 GeV/c
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","1111100057022230000","0h63103100000010"); //standart cut, kINT7
  } else if (trainConfig == 2){  // min energy = 0.3 GeV/c
    cuts.AddCutPCMCalo("80052113","00200009f9730000dge0400000","1111100057022230000","0h63103100000010"); //standard cut, kEMC7
    cuts.AddCutPCMCalo("80083113","00200009f9730000dge0400000","1111100057022230000","0h63103100000010"); //standard cut, kEMCEG1 based on INT7
    cuts.AddCutPCMCalo("80085113","00200009f9730000dge0400000","1111100057022230000","0h63103100000010"); //standard cut, kEMCEG2 based on INT7
  } else if (trainConfig == 3){ // min energy = 0.4 GeV/c
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","1111100057032230000","0h63103100000010"); //standart cut, kINT7
  } else if (trainConfig == 4){ // min energy = 0.4 GeV/
    cuts.AddCutPCMCalo("80052113","00200009f9730000dge0400000","1111100057032230000","0h63103100000010"); //standard cut, kEMC7
    cuts.AddCutPCMCalo("80083113","00200009f9730000dge0400000","1111100057032230000","0h63103100000010"); //standard cut, kEMCEG1 based on INT7
    cuts.AddCutPCMCalo("80085113","00200009f9730000dge0400000","1111100057032230000","0h63103100000010"); //standard cut, kEMCEG2 based on INT7
  } else if (trainConfig == 5){ // min energy = 0.4 GeV/c w/o time cut
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","1111100007032230000","0h63103100000010"); //standart cut, kINT7
  } else if (trainConfig == 6){ // min energy = 0.4 GeV/ w/o time cut
    cuts.AddCutPCMCalo("80052113","00200009f9730000dge0400000","1111100007032230000","0h63103100000010"); //standard cut, kEMC7
    cuts.AddCutPCMCalo("80083113","00200009f9730000dge0400000","1111100007032230000","0h63103100000010"); //standard cut, kEMCEG1 based on INT7
    cuts.AddCutPCMCalo("80085113","00200009f9730000dge0400000","1111100007032230000","0h63103100000010"); //standard cut, kEMCEG2 based on INT7
  } else if (trainConfig == 7){ // introduce fast future protection
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","1111100007032230000","0h63103100000010"); // no PF
    cuts.AddCutPCMCalo("80010313","00200009f9730000dge0400000","1111100007032230000","0h63103100000010"); // 0.1 \mus protected
    cuts.AddCutPCMCalo("80010413","00200009f9730000dge0400000","1111100007032230000","0h63103100000010"); // 0.25 \mus protected
    cuts.AddCutPCMCalo("80010513","00200009f9730000dge0400000","1111100007032230000","0h63103100000010"); // 1.075 \mus protected
  } else if (trainConfig == 8){ // PF for cent var, 1.075 \mus protected
    cuts.AddCutPCMCalo("80210513","00200009f9730000dge0400000","1111100007032230000","0h63103100000010"); // 0-20
    cuts.AddCutPCMCalo("82410513","00200009f9730000dge0400000","1111100007032230000","0h63103100000010"); // 20-40
    cuts.AddCutPCMCalo("84610513","00200009f9730000dge0400000","1111100007032230000","0h63103100000010"); // 40-60
    cuts.AddCutPCMCalo("86010513","00200009f9730000dge0400000","1111100007032230000","0h63103100000010"); // 60-100
  //---------------------------------------------------------------------------------------------
  // standard cuts
  //---------------------------------------------------------------------------------------------
  } else if (trainConfig == 9){ // EMCAL clusters standard cut MB + cent dependent
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","1111141057032230000","0i63103100000010"); // 0-100
    cuts.AddCutPCMCalo("80210113","00200009f9730000dge0400000","1111141057032230000","0i63103100000010"); // 0-20
    cuts.AddCutPCMCalo("82410113","00200009f9730000dge0400000","1111141057032230000","0i63103100000010"); // 20-40
    cuts.AddCutPCMCalo("84610113","00200009f9730000dge0400000","1111141057032230000","0i63103100000010"); // 40-60
    cuts.AddCutPCMCalo("86010113","00200009f9730000dge0400000","1111141057032230000","0i63103100000010"); // 60-100
  } else if (trainConfig == 10){ // introduce fast future protection
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","1111141007032230000","0h63103100000010"); // no PF
    cuts.AddCutPCMCalo("80010313","00200009f9730000dge0400000","1111141007032230000","0h63103100000010"); // 0.1 \mus protected
    cuts.AddCutPCMCalo("80010413","00200009f9730000dge0400000","1111141007032230000","0h63103100000010"); // 0.25 \mus protected
    cuts.AddCutPCMCalo("80010513","00200009f9730000dge0400000","1111141007032230000","0h63103100000010"); // 1.075 \mus protected
  } else if (trainConfig == 11) {
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","2444451044013200000","0h63103100000010"); // standart cut, kINT7
    cuts.AddCutPCMCalo("80062113","00200009f9730000dge0400000","2444451044013200000","0h63103100000010"); // standard cut, kPHI7

  } else if (trainConfig == 12){ // EMCAL clusters standard cuts cent dependent
    cuts.AddCutPCMCalo("80210113","00200009f9730000dge0400000","1111141057032230000","0h63103100000010"); // 0-20
  } else if (trainConfig == 13){ // EMCAL clusters standard cuts cent dependent
    cuts.AddCutPCMCalo("82410113","00200009f9730000dge0400000","1111141057032230000","0h63103100000010"); // 20-40
  } else if (trainConfig == 14){ // EMCAL clusters standard cuts cent dependent
    cuts.AddCutPCMCalo("84610113","00200009f9730000dge0400000","1111141057032230000","0h63103100000010"); // 40-60
  } else if (trainConfig == 15){ // EMCAL clusters standard cuts cent dependent
    cuts.AddCutPCMCalo("86010113","00200009f9730000dge0400000","1111141057032230000","0h63103100000010"); // 60-100
  } else if (trainConfig == 16){ // EMCAL clusters standard cut MB
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","1111141057032230000","0h63103100000010"); // INT7
  } else if (trainConfig == 17){ // EMCAL clusters standard cut MB + cent dependent
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","1111141057032230000","0h63103100000010"); // 0-100
    cuts.AddCutPCMCalo("80210113","00200009f9730000dge0400000","1111141057032230000","0h63103100000010"); // 0-20
    cuts.AddCutPCMCalo("82410113","00200009f9730000dge0400000","1111141057032230000","0h63103100000010"); // 20-40
    cuts.AddCutPCMCalo("84610113","00200009f9730000dge0400000","1111141057032230000","0h63103100000010"); // 40-60
    cuts.AddCutPCMCalo("86010113","00200009f9730000dge0400000","1111141057032230000","0h63103100000010"); // 60-100
  } else if (trainConfig == 18){ // EMCAL clusters standard cut MB + cent dependent CL1 est
    cuts.AddCutPCMCalo("90010113","00200009f9730000dge0400000","1111141057032230000","0h63103100000010"); // 0-20
    cuts.AddCutPCMCalo("90210113","00200009f9730000dge0400000","1111141057032230000","0h63103100000010"); // 0-20
    cuts.AddCutPCMCalo("92410113","00200009f9730000dge0400000","1111141057032230000","0h63103100000010"); // 20-40
    cuts.AddCutPCMCalo("94610113","00200009f9730000dge0400000","1111141057032230000","0h63103100000010"); // 40-60
    cuts.AddCutPCMCalo("96010113","00200009f9730000dge0400000","1111141057032230000","0h63103100000010"); // 60-100
  } else if (trainConfig == 19){ // EMCAL clusters standard cut MB + cent dependent CL1 est
    cuts.AddCutPCMCalo("e0010113","00200009f9730000dge0400000","1111141057032230000","0h63103100000010"); // 0-20
    cuts.AddCutPCMCalo("e0210113","00200009f9730000dge0400000","1111141057032230000","0h63103100000010"); // 0-20
    cuts.AddCutPCMCalo("e2410113","00200009f9730000dge0400000","1111141057032230000","0h63103100000010"); // 20-40
    cuts.AddCutPCMCalo("e4610113","00200009f9730000dge0400000","1111141057032230000","0h63103100000010"); // 40-60
    cuts.AddCutPCMCalo("e6010113","00200009f9730000dge0400000","1111141057032230000","0h63103100000010"); // 60-100

  } else if (trainConfig == 20){ // EMCAL clusters standard cuts
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","1111141057032230000","0h63103100000010"); // INT7
    cuts.AddCutPCMCalo("80052113","00200009f9730000dge0400000","1111141057032230000","0h63103100000010"); // EMC7
    cuts.AddCutPCMCalo("80085113","00200009f9730000dge0400000","1111141057032230000","0h63103100000010"); // EG2
    cuts.AddCutPCMCalo("80083113","00200009f9730000dge0400000","1111141057032230000","0h63103100000010"); // EG1

  //---------------------------------------------------------------------------------------------
  // minimum bias variations
  //---------------------------------------------------------------------------------------------
  } else if (trainConfig == 21){ //EMCAL minEnergy variation
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","1111141057032230000","0h63103100000010"); // 0.7 GeV/c default
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","1111141057042230000","0h63103100000010"); // 0.8 GeV/c
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","1111141057052230000","0h63103100000010"); // 0.9 GeV/c
  } else if (trainConfig == 22){ //EMCAL minNCells variation
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","1111141057031230000","0h63103100000010"); // n cells >= 1
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","1111141057033230000","0h63103100000010"); // n cells >= 3
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","1111141057032200000","0h63103100000010"); // no max M02 cut
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","1111141057032250000","0h63103100000010"); // M02 < 0.3
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","1111141057032260000","0h63103100000010"); // M02 < 0.27
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","1112141057032230000","0h63103100000010"); // only modules with TRD infront
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","1111341057032230000","0h63103100000010"); // no modules with TRD infront
  } else if (trainConfig == 23){ // EMCAL track matching variations
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","1111141053032230000","0h63103100000010"); // fixed TM
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","1111141056032230000","0h63103100000010"); // pt dep var 1
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","1111141058032230000","0h63103100000010"); // pt dep var 1
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","1111141059032230000","0h63103100000010"); // pt dep var 1
  } else if (trainConfig == 24){ // PCM variations
    cuts.AddCutPCMCalo("80010113","00200009227000008250400000","1111141057032230000","0h63103100000010"); // dEdx e -3, 5
    cuts.AddCutPCMCalo("80010113","00200009127000008250400000","1111141057032230000","0h63103100000010"); // dEdx e -5, 5
    cuts.AddCutPCMCalo("80010113","00200009357000008250400000","1111141057032230000","0h63103100000010"); // dEdx pi 2
    cuts.AddCutPCMCalo("80010113","00200009317000008250400000","1111141057032230000","0h63103100000010"); // dEdx pi 0
    cuts.AddCutPCMCalo("80010113","00200009387300008250400000","1111141057032230000","0h63103100000010"); // dEdx pi 2 high 1 (> 3.5 GeV)
  } else if (trainConfig == 25){ // PCM variations
    cuts.AddCutPCMCalo("80010113","00200009327000009250400000","1111141057032230000","0h63103100000010"); // qt 2D 0.03
    cuts.AddCutPCMCalo("80010113","00200009327000003250400000","1111141057032230000","0h63103100000010"); // qt 1D 0.05
    cuts.AddCutPCMCalo("80010113","00200009327000002250400000","1111141057032230000","0h63103100000010"); // qt 1D 0.07
    cuts.AddCutPCMCalo("80010113","00200049327000008250400000","1111141057032230000","0h63103100000010"); // single pt > 0.075
    cuts.AddCutPCMCalo("80010113","00200019327000008250400000","1111141057032230000","0h63103100000010"); // single pt > 0.1
  } else if (trainConfig == 26){ // PCM variations
    cuts.AddCutPCMCalo("80010113","00200009327000008850400000","1111141057032230000","0h63103100000010"); // 2D psi pair chi2 var
    cuts.AddCutPCMCalo("80010113","00200009327000008260400000","1111141057032230000","0h63103100000010"); // 2D psi pair chi2 var
    cuts.AddCutPCMCalo("80010113","00200009327000008860400000","1111141057032230000","0h63103100000010"); // 2D psi pair chi2 var
    cuts.AddCutPCMCalo("80010113","00200009327000008280400000","1111141057032230000","0h63103100000010"); // 2D psi pair chi2 var
    cuts.AddCutPCMCalo("80010113","00200009327000008880400000","1111141057032230000","0h63103100000010"); // 2D psi pair chi2 var
  } else if (trainConfig == 27){ // PCM variations
    cuts.AddCutPCMCalo("80010113","00200006327000008250400000","1111141057032230000","0h63103100000010"); // min TPC cl > 0.7
    cuts.AddCutPCMCalo("80010113","00200008327000008250400000","1111141057032230000","0h63103100000010"); // min TPC cl > 0.35
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","1111141057032230000","0h63106100000010"); // alpha < 0.8
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","1111141057032230000","0h63105100000010"); // alpha < 0.75
  } else if (trainConfig == 28){ // PCM variations
    cuts.AddCutPCMCalo("80010113","00202209327000008250400000","1111141057032230000","0h63103100000010"); // restrict acceptance to EMCAL loose
    cuts.AddCutPCMCalo("80010113","00204409327000008250400000","1111141057032230000","0h63103100000010"); // restrict acceptance to EMCAL tight
  } else if (trainConfig == 29){ // PCM variations pi dEdx
    cuts.AddCutPCMCalo("80010113","00200009317300008250400000","1111141057032230000","0h63103100000010"); // dEdx pi: 0: 0.4-3.5, -10: 3.5 ->
    cuts.AddCutPCMCalo("80010113","00200009327300008250400000","1111141057032230000","0h63103100000010"); // dEdx pi: 1: 0.4-3.5, -10: 3.5 ->
    cuts.AddCutPCMCalo("80010113","00200009325000008250400000","1111141057032230000","0h63103100000010"); // dEdx pi: 1: 0.3 ->
    cuts.AddCutPCMCalo("80010113","00200009320000008250400000","1111141057032230000","0h63103100000010"); // dEdx pi: 1: 0.5 ->
  } else if (trainConfig == 30){ // PCM variations pi dEdx
    cuts.AddCutPCMCalo("80010113","00200009327600008250400000","1111141057032230000","0h63103100000010"); // dEdx pi: 1: 0.4-2, -10: 2. ->
    cuts.AddCutPCMCalo("80010113","00200009327400008250400000","1111141057032230000","0h63103100000010"); // dEdx pi: 1: 0.4-3, -10: 3. ->
    cuts.AddCutPCMCalo("80010113","00200009315600008250400000","1111141057032230000","0h63103100000010"); // dEdx pi: 0: 0.3-2, -10: 2. ->
    cuts.AddCutPCMCalo("80010113","00200009367400008250400000","1111141057032230000","0h63103100000010"); // dEdx pi: 2: 0.4-3, 0.5: 3. ->
    cuts.AddCutPCMCalo("80010113","00200009347400008250400000","1111141057032230000","0h63103100000010"); // dEdx pi: 3: 0.4-3, 1: 3. ->
  } else if (trainConfig == 31){ // PCM variations to close V0s
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","1111141057032230000","0h63103100000010");
    cuts.AddCutPCMCalo("80010113","00200009327000008250401000","1111141057032230000","0h63103100000010");
    cuts.AddCutPCMCalo("80010113","00200009327000008250402000","1111141057032230000","0h63103100000010");
    cuts.AddCutPCMCalo("80010113","00200009327000008250403000","1111141057032230000","0h63103100000010");
  } else if (trainConfig == 32){ // EMCal cluster, non lin variations
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","1111141057032230000","0h63103100000010");
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","1111142057032230000","0h63103100000010");
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","1111151057032230000","0h63103100000010");
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","1111152057032230000","0h63103100000010");
  } else if (trainConfig == 33){  // EMCal cluster, non lin variations
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","1111100057032230000","0h63103100000010");
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","1111101057032230000","0h63103100000010");
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","1111102057032230000","0h63103100000010");
  } else if (trainConfig == 34){ // EMCal cluster, non lin variations
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","1111143057032230000","0h63103100000010");
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","1111144057032230000","0h63103100000010");
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","1111153057032230000","0h63103100000010");
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","1111154057032230000","0h63103100000010");

  } else if (trainConfig == 35){ // EMCAL clusters standard cuts TESTING CONV CALO MIXING
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","1111141057032230000","0163103100000010"); // INT7
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","1111141057032230000","0h63103100000010"); // INT7

  } else if (trainConfig == 36){ // EMCAL clusters standard cut MB, new TB Nico
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","1111100057032230000","0h63103100000010"); // INT7
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","1111141057032230000","0h63103100000010"); // INT7
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","1111142057032230000","0h63103100000010"); // INT7
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","1111165057032230000","0h63103100000010"); // INT7
  } else if (trainConfig == 37){ // EMCAL clusters standard cut MB, new TB Nico
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","1111100057032230000","0h63103100000010"); // INT7
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","3885500057032230000","0h63103100000010"); // INT7
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","4117900057032230000","0h63103100000010"); // INT7

  //---------------------------------------------------------------------------------------------
  // 0-20 % variations
  //---------------------------------------------------------------------------------------------
  } else if (trainConfig == 60){ // EMCAL clusters standard cuts
    cuts.AddCutPCMCalo("80210113","00200009f9730000dge0400000","1111141057032230000","0h63103100000010"); // INT7
    cuts.AddCutPCMCalo("80252113","00200009f9730000dge0400000","1111141057032230000","0h63103100000010"); // EMC7
    cuts.AddCutPCMCalo("80285113","00200009f9730000dge0400000","1111141057022230000","0h63103100000010"); // EG2
    cuts.AddCutPCMCalo("80283113","00200009f9730000dge0400000","1111141057022230000","0h63103100000010"); // EG1
  } else if (trainConfig == 61){ //EMCAL minEnergy variation
    cuts.AddCutPCMCalo("80210113","00200009f9730000dge0400000","1111141057032230000","0h63103100000010"); // 0.6 GeV/c default
    cuts.AddCutPCMCalo("80210113","00200009f9730000dge0400000","1111141057042230000","0h63103100000010"); // 0.7 GeV/c
    cuts.AddCutPCMCalo("80210113","00200009f9730000dge0400000","1111141057052230000","0h63103100000010"); // 0.9 GeV/c
  } else if (trainConfig == 62){ //EMCAL minNCells variation
    cuts.AddCutPCMCalo("80210113","00200009f9730000dge0400000","1111141057031230000","0h63103100000010"); // n cells >= 1
    cuts.AddCutPCMCalo("80210113","00200009f9730000dge0400000","1111141057033230000","0h63103100000010"); // n cells >= 3
    cuts.AddCutPCMCalo("80210113","00200009f9730000dge0400000","1111141057032200000","0h63103100000010"); // no max M02 cut
    cuts.AddCutPCMCalo("80210113","00200009f9730000dge0400000","1111141057032250000","0h63103100000010"); // M02 < 0.3
    cuts.AddCutPCMCalo("80210113","00200009f9730000dge0400000","1111141057032260000","0h63103100000010"); // M02 < 0.27
    cuts.AddCutPCMCalo("80210113","00200009f9730000dge0400000","1112141057032230000","0h63103100000010"); // only modules with TRD infront
    cuts.AddCutPCMCalo("80210113","00200009f9730000dge0400000","1111341057032230000","0h63103100000010"); // no modules with TRD infront
  } else if (trainConfig == 63){ // EMCAL track matching variations
    cuts.AddCutPCMCalo("80210113","00200009f9730000dge0400000","1111141053032230000","0h63103100000010"); //
    cuts.AddCutPCMCalo("80210113","00200009f9730000dge0400000","1111141056032230000","0h63103100000010"); //
    cuts.AddCutPCMCalo("80210113","00200009f9730000dge0400000","1111141058032230000","0h63103100000010"); //
    cuts.AddCutPCMCalo("80210113","00200009f9730000dge0400000","1111141059032230000","0h63103100000010"); //
  } else if (trainConfig == 64){ // PCM variations
    cuts.AddCutPCMCalo("80210113","00200009227000008250400000","1111141057032230000","0h63103100000010"); // dEdx e -3, 5
    cuts.AddCutPCMCalo("80210113","00200009127000008250400000","1111141057032230000","0h63103100000010"); // dEdx e -5, 5
    cuts.AddCutPCMCalo("80210113","00200009357000008250400000","1111141057032230000","0h63103100000010"); // dEdx pi 2
    cuts.AddCutPCMCalo("80210113","00200009317000008250400000","1111141057032230000","0h63103100000010"); // dEdx pi 0
    cuts.AddCutPCMCalo("80210113","00200009387300008250400000","1111141057032230000","0h63103100000010"); // dEdx pi 2 high 1 (> 3.5 GeV)
  } else if (trainConfig == 65){ // PCM variations
    cuts.AddCutPCMCalo("80210113","00200009327000009250400000","1111141057032230000","0h63103100000010"); // qt 2D 0.03
    cuts.AddCutPCMCalo("80210113","00200009327000003250400000","1111141057032230000","0h63103100000010"); // qt 1D 0.05
    cuts.AddCutPCMCalo("80210113","00200009327000002250400000","1111141057032230000","0h63103100000010"); // qt 1D 0.07
    cuts.AddCutPCMCalo("80210113","00200049327000008250400000","1111141057032230000","0h63103100000010"); // single pt > 0.075
    cuts.AddCutPCMCalo("80210113","00200019327000008250400000","1111141057032230000","0h63103100000010"); // single pt > 0.1
  } else if (trainConfig == 66){ // PCM variations
    cuts.AddCutPCMCalo("80210113","00200009327000008850400000","1111141057032230000","0h63103100000010"); // 2D psi pair chi2 var
    cuts.AddCutPCMCalo("80210113","00200009327000008260400000","1111141057032230000","0h63103100000010"); // 2D psi pair chi2 var
    cuts.AddCutPCMCalo("80210113","00200009327000008860400000","1111141057032230000","0h63103100000010"); // 2D psi pair chi2 var
    cuts.AddCutPCMCalo("80210113","00200009327000008280400000","1111141057032230000","0h63103100000010"); // 2D psi pair chi2 var
    cuts.AddCutPCMCalo("80210113","00200009327000008880400000","1111141057032230000","0h63103100000010"); // 2D psi pair chi2 var
  } else if (trainConfig == 67){ // PCM variations
    cuts.AddCutPCMCalo("80210113","00200006327000008250400000","1111141057032230000","0h63103100000010"); // min TPC cl > 0.7
    cuts.AddCutPCMCalo("80210113","00200008327000008250400000","1111141057032230000","0h63103100000010"); // min TPC cl > 0.35
    cuts.AddCutPCMCalo("80210113","00200009f9730000dge0400000","1111141057032230000","0h63106100000010"); // alpha < 0.8
    cuts.AddCutPCMCalo("80210113","00200009f9730000dge0400000","1111141057032230000","0h63105100000010"); // alpha < 0.75
  } else if (trainConfig == 68){ // PCM variations
    cuts.AddCutPCMCalo("80210113","00202209327000008250400000","1111141057032230000","0h63103100000010"); // restrict acceptance to EMCAL loose
    cuts.AddCutPCMCalo("80210113","00204409327000008250400000","1111141057032230000","0h63103100000010"); // restrict acceptance to EMCAL tight
  } else if (trainConfig == 69){ // PCM variations pi dEdx
    cuts.AddCutPCMCalo("80210113","00200009317300008250400000","1111141057032230000","0h63103100000010"); // dEdx pi: 0: 0.4-3.5, -10: 3.5 ->
    cuts.AddCutPCMCalo("80210113","00200009327300008250400000","1111141057032230000","0h63103100000010"); // dEdx pi: 1: 0.4-3.5, -10: 3.5 ->
    cuts.AddCutPCMCalo("80210113","00200009325000008250400000","1111141057032230000","0h63103100000010"); // dEdx pi: 1: 0.3 ->
    cuts.AddCutPCMCalo("80210113","00200009320000008250400000","1111141057032230000","0h63103100000010"); // dEdx pi: 1: 0.5 ->
  } else if (trainConfig == 70){ // PCM variations pi dEdx
    cuts.AddCutPCMCalo("80210113","00200009327600008250400000","1111141057032230000","0h63103100000010"); // dEdx pi: 1: 0.4-2, -10: 2. ->
    cuts.AddCutPCMCalo("80210113","00200009327400008250400000","1111141057032230000","0h63103100000010"); // dEdx pi: 1: 0.4-3, -10: 3. ->
    cuts.AddCutPCMCalo("80210113","00200009315600008250400000","1111141057032230000","0h63103100000010"); // dEdx pi: 0: 0.3-2, -10: 2. ->
    cuts.AddCutPCMCalo("80210113","00200009367400008250400000","1111141057032230000","0h63103100000010"); // dEdx pi: 2: 0.4-3, 0.5: 3. ->
    cuts.AddCutPCMCalo("80210113","00200009347400008250400000","1111141057032230000","0h63103100000010"); // dEdx pi: 3: 0.4-3, 1: 3. ->
  } else if (trainConfig == 71){ // PCM variations to close V0s
    cuts.AddCutPCMCalo("80210113","00200009f9730000dge0400000","1111141057032230000","0h63103100000010");
    cuts.AddCutPCMCalo("80210113","00200009327000008250401000","1111141057032230000","0h63103100000010");
    cuts.AddCutPCMCalo("80210113","00200009327000008250402000","1111141057032230000","0h63103100000010");
    cuts.AddCutPCMCalo("80210113","00200009327000008250403000","1111141057032230000","0h63103100000010");
  } else if (trainConfig == 72){ // EMCal cluster, non lin variations
    cuts.AddCutPCMCalo("80210113","00200009f9730000dge0400000","1111142057032230000","0h63103100000010");
    cuts.AddCutPCMCalo("80210113","00200009f9730000dge0400000","1111143057032230000","0h63103100000010");
    cuts.AddCutPCMCalo("80210113","00200009f9730000dge0400000","1111151057032230000","0h63103100000010");
    cuts.AddCutPCMCalo("80210113","00200009f9730000dge0400000","1111152057032230000","0h63103100000010");
  } else if (trainConfig == 73){  // EMCal cluster, non lin variations
    cuts.AddCutPCMCalo("80210113","00200009f9730000dge0400000","1111100057032230000","0h63103100000010");
    cuts.AddCutPCMCalo("80210113","00200009f9730000dge0400000","1111101057032230000","0h63103100000010");
    cuts.AddCutPCMCalo("80210113","00200009f9730000dge0400000","1111102057032230000","0h63103100000010");
    cuts.AddCutPCMCalo("80210113","00200009f9730000dge0400000","1111144057032230000","0h63103100000010");

  //---------------------------------------------------------------------------------------------
  // 20-40 % variations
  //---------------------------------------------------------------------------------------------
  } else if (trainConfig == 80){ // EMCAL clusters standard cuts
    cuts.AddCutPCMCalo("82410113","00200009f9730000dge0400000","1111141057032230000","0h63103100000010"); // INT7
    cuts.AddCutPCMCalo("82452113","00200009f9730000dge0400000","1111141057032230000","0h63103100000010"); // EMC7
    cuts.AddCutPCMCalo("82485113","00200009f9730000dge0400000","1111141057022230000","0h63103100000010"); // EG2
    cuts.AddCutPCMCalo("82483113","00200009f9730000dge0400000","1111141057022230000","0h63103100000010"); // EG1
  } else if (trainConfig == 81){ //EMCAL minEnergy variation
    cuts.AddCutPCMCalo("82410113","00200009f9730000dge0400000","1111141057032230000","0h63103100000010"); // 0.7 GeV/c default
    cuts.AddCutPCMCalo("82410113","00200009f9730000dge0400000","1111141057042230000","0h63103100000010"); // 0.8 GeV/c
    cuts.AddCutPCMCalo("82410113","00200009f9730000dge0400000","1111141057052230000","0h63103100000010"); // 0.9 GeV/c
  } else if (trainConfig == 82){ //EMCAL minNCells variation
    cuts.AddCutPCMCalo("82410113","00200009f9730000dge0400000","1111141057031230000","0h63103100000010"); // n cells >= 1
    cuts.AddCutPCMCalo("82410113","00200009f9730000dge0400000","1111141057033230000","0h63103100000010"); // n cells >= 3
    cuts.AddCutPCMCalo("82410113","00200009f9730000dge0400000","1111141057032200000","0h63103100000010"); // no max M02 cut
    cuts.AddCutPCMCalo("82410113","00200009f9730000dge0400000","1111141057032250000","0h63103100000010"); // M02 < 0.3
    cuts.AddCutPCMCalo("82410113","00200009f9730000dge0400000","1111141057032260000","0h63103100000010"); // M02 < 0.27
    cuts.AddCutPCMCalo("82410113","00200009f9730000dge0400000","1112141057032230000","0h63103100000010"); // only modules with TRD infront
    cuts.AddCutPCMCalo("82410113","00200009f9730000dge0400000","1111341057032230000","0h63103100000010"); // no modules with TRD infront
  } else if (trainConfig == 83){ // EMCAL track matching variations
    cuts.AddCutPCMCalo("82410113","00200009f9730000dge0400000","1111141053032230000","0h63103100000010"); //
    cuts.AddCutPCMCalo("82410113","00200009f9730000dge0400000","1111141056032230000","0h63103100000010"); //
    cuts.AddCutPCMCalo("82410113","00200009f9730000dge0400000","1111141058032230000","0h63103100000010"); //
    cuts.AddCutPCMCalo("82410113","00200009f9730000dge0400000","1111141059032230000","0h63103100000010"); //
  } else if (trainConfig == 84){ // PCM variations
    cuts.AddCutPCMCalo("82410113","00200009227000008250400000","1111141057032230000","0h63103100000010"); // dEdx e -3, 5
    cuts.AddCutPCMCalo("82410113","00200009127000008250400000","1111141057032230000","0h63103100000010"); // dEdx e -5, 5
    cuts.AddCutPCMCalo("82410113","00200009357000008250400000","1111141057032230000","0h63103100000010"); // dEdx pi 2
    cuts.AddCutPCMCalo("82410113","00200009317000008250400000","1111141057032230000","0h63103100000010"); // dEdx pi 0
    cuts.AddCutPCMCalo("82410113","00200009387300008250400000","1111141057032230000","0h63103100000010"); // dEdx pi 2 high 1 (> 3.5 GeV)
  } else if (trainConfig == 85){ // PCM variations
    cuts.AddCutPCMCalo("82410113","00200009327000009250400000","1111141057032230000","0h63103100000010"); // qt 2D 0.03
    cuts.AddCutPCMCalo("82410113","00200009327000003250400000","1111141057032230000","0h63103100000010"); // qt 1D 0.05
    cuts.AddCutPCMCalo("82410113","00200009327000002250400000","1111141057032230000","0h63103100000010"); // qt 1D 0.07
    cuts.AddCutPCMCalo("82410113","00200049327000008250400000","1111141057032230000","0h63103100000010"); // single pt > 0.075
    cuts.AddCutPCMCalo("82410113","00200019327000008250400000","1111141057032230000","0h63103100000010"); // single pt > 0.1
  } else if (trainConfig == 86){ // PCM variations
    cuts.AddCutPCMCalo("82410113","00200009327000008850400000","1111141057032230000","0h63103100000010"); // 2D psi pair chi2 var
    cuts.AddCutPCMCalo("82410113","00200009327000008260400000","1111141057032230000","0h63103100000010"); // 2D psi pair chi2 var
    cuts.AddCutPCMCalo("82410113","00200009327000008860400000","1111141057032230000","0h63103100000010"); // 2D psi pair chi2 var
    cuts.AddCutPCMCalo("82410113","00200009327000008280400000","1111141057032230000","0h63103100000010"); // 2D psi pair chi2 var
    cuts.AddCutPCMCalo("82410113","00200009327000008880400000","1111141057032230000","0h63103100000010"); // 2D psi pair chi2 var
  } else if (trainConfig == 87){ // PCM variations
    cuts.AddCutPCMCalo("82410113","00200006327000008250400000","1111141057032230000","0h63103100000010"); // min TPC cl > 0.7
    cuts.AddCutPCMCalo("82410113","00200008327000008250400000","1111141057032230000","0h63103100000010"); // min TPC cl > 0.35
    cuts.AddCutPCMCalo("82410113","00200009f9730000dge0400000","1111141057032230000","0h63106100000010"); // alpha < 0.8
    cuts.AddCutPCMCalo("82410113","00200009f9730000dge0400000","1111141057032230000","0h63105100000010"); // alpha < 0.75
  } else if (trainConfig == 88){ // PCM variations
    cuts.AddCutPCMCalo("82410113","00202209327000008250400000","1111141057032230000","0h63103100000010"); // restrict acceptance to EMCAL loose
    cuts.AddCutPCMCalo("82410113","00204409327000008250400000","1111141057032230000","0h63103100000010"); // restrict acceptance to EMCAL tight
  } else if (trainConfig == 89){ // PCM variations pi dEdx
    cuts.AddCutPCMCalo("82410113","00200009317300008250400000","1111141057032230000","0h63103100000010"); // dEdx pi: 0: 0.4-3.5, -10: 3.5 ->
    cuts.AddCutPCMCalo("82410113","00200009327300008250400000","1111141057032230000","0h63103100000010"); // dEdx pi: 1: 0.4-3.5, -10: 3.5 ->
    cuts.AddCutPCMCalo("82410113","00200009325000008250400000","1111141057032230000","0h63103100000010"); // dEdx pi: 1: 0.3 ->
    cuts.AddCutPCMCalo("82410113","00200009320000008250400000","1111141057032230000","0h63103100000010"); // dEdx pi: 1: 0.5 ->
  } else if (trainConfig == 90){ // PCM variations pi dEdx
    cuts.AddCutPCMCalo("82410113","00200009327600008250400000","1111141057032230000","0h63103100000010"); // dEdx pi: 1: 0.4-2, -10: 2. ->
    cuts.AddCutPCMCalo("82410113","00200009327400008250400000","1111141057032230000","0h63103100000010"); // dEdx pi: 1: 0.4-3, -10: 3. ->
    cuts.AddCutPCMCalo("82410113","00200009315600008250400000","1111141057032230000","0h63103100000010"); // dEdx pi: 0: 0.3-2, -10: 2. ->
    cuts.AddCutPCMCalo("82410113","00200009367400008250400000","1111141057032230000","0h63103100000010"); // dEdx pi: 2: 0.4-3, 0.5: 3. ->
    cuts.AddCutPCMCalo("82410113","00200009347400008250400000","1111141057032230000","0h63103100000010"); // dEdx pi: 3: 0.4-3, 1: 3. ->
  } else if (trainConfig == 91){ // PCM variations to close V0s
    cuts.AddCutPCMCalo("82410113","00200009f9730000dge0400000","1111141057032230000","0h63103100000010");
    cuts.AddCutPCMCalo("82410113","00200009327000008250401000","1111141057032230000","0h63103100000010");
    cuts.AddCutPCMCalo("82410113","00200009327000008250402000","1111141057032230000","0h63103100000010");
    cuts.AddCutPCMCalo("82410113","00200009327000008250403000","1111141057032230000","0h63103100000010");
  } else if (trainConfig == 92){ // EMCal cluster, non lin variations
    cuts.AddCutPCMCalo("82410113","00200009f9730000dge0400000","1111142057032230000","0h63103100000010");
    cuts.AddCutPCMCalo("82410113","00200009f9730000dge0400000","1111143057032230000","0h63103100000010");
    cuts.AddCutPCMCalo("82410113","00200009f9730000dge0400000","1111151057032230000","0h63103100000010");
    cuts.AddCutPCMCalo("82410113","00200009f9730000dge0400000","1111152057032230000","0h63103100000010");
  } else if (trainConfig == 93){  // EMCal cluster, non lin variations
    cuts.AddCutPCMCalo("82410113","00200009f9730000dge0400000","1111100057032230000","0h63103100000010");
    cuts.AddCutPCMCalo("82410113","00200009f9730000dge0400000","1111101057032230000","0h63103100000010");
    cuts.AddCutPCMCalo("82410113","00200009f9730000dge0400000","1111102057032230000","0h63103100000010");
    cuts.AddCutPCMCalo("82410113","00200009f9730000dge0400000","1111144057032230000","0h63103100000010");

  //---------------------------------------------------------------------------------------------
  // 40-60 % variations
  //---------------------------------------------------------------------------------------------
  } else if (trainConfig == 100){ // EMCAL clusters standard cuts
    cuts.AddCutPCMCalo("84610113","00200009f9730000dge0400000","1111141057032230000","0h63103100000010"); // INT7
    cuts.AddCutPCMCalo("84652113","00200009f9730000dge0400000","1111141057032230000","0h63103100000010"); // EMC7
    cuts.AddCutPCMCalo("84685113","00200009f9730000dge0400000","1111141057022230000","0h63103100000010"); // EG2
    cuts.AddCutPCMCalo("84683113","00200009f9730000dge0400000","1111141057022230000","0h63103100000010"); // EG1
  } else if (trainConfig == 101){ //EMCAL minEnergy variation
    cuts.AddCutPCMCalo("84610113","00200009f9730000dge0400000","1111141057032230000","0h63103100000010"); // 0.7 GeV/c default
    cuts.AddCutPCMCalo("84610113","00200009f9730000dge0400000","1111141057042230000","0h63103100000010"); // 0.8 GeV/c
    cuts.AddCutPCMCalo("84610113","00200009f9730000dge0400000","1111141057052230000","0h63103100000010"); // 0.9 GeV/c
  } else if (trainConfig == 102){ //EMCAL minNCells variation
    cuts.AddCutPCMCalo("84610113","00200009f9730000dge0400000","1111141057031230000","0h63103100000010"); // n cells >= 1
    cuts.AddCutPCMCalo("84610113","00200009f9730000dge0400000","1111141057033230000","0h63103100000010"); // n cells >= 3
    cuts.AddCutPCMCalo("84610113","00200009f9730000dge0400000","1111141057032200000","0h63103100000010"); // no max M02 cut
    cuts.AddCutPCMCalo("84610113","00200009f9730000dge0400000","1111141057032250000","0h63103100000010"); // M02 < 0.3
    cuts.AddCutPCMCalo("84610113","00200009f9730000dge0400000","1111141057032260000","0h63103100000010"); // M02 < 0.27
    cuts.AddCutPCMCalo("84610113","00200009f9730000dge0400000","1112141057032230000","0h63103100000010"); // only modules with TRD infront
    cuts.AddCutPCMCalo("84610113","00200009f9730000dge0400000","1111341057032230000","0h63103100000010"); // no modules with TRD infront
  } else if (trainConfig == 103){ // EMCAL track matching variations
    cuts.AddCutPCMCalo("84610113","00200009f9730000dge0400000","1111141053032230000","0h63103100000010"); //
    cuts.AddCutPCMCalo("84610113","00200009f9730000dge0400000","1111141056032230000","0h63103100000010"); //
    cuts.AddCutPCMCalo("84610113","00200009f9730000dge0400000","1111141058032230000","0h63103100000010"); //
    cuts.AddCutPCMCalo("84610113","00200009f9730000dge0400000","1111141059032230000","0h63103100000010"); //
  } else if (trainConfig == 104){ // PCM variations
    cuts.AddCutPCMCalo("84610113","00200009227000008250400000","1111141057032230000","0h63103100000010"); // dEdx e -3, 5
    cuts.AddCutPCMCalo("84610113","00200009127000008250400000","1111141057032230000","0h63103100000010"); // dEdx e -5, 5
    cuts.AddCutPCMCalo("84610113","00200009357000008250400000","1111141057032230000","0h63103100000010"); // dEdx pi 2
    cuts.AddCutPCMCalo("84610113","00200009317000008250400000","1111141057032230000","0h63103100000010"); // dEdx pi 0
    cuts.AddCutPCMCalo("84610113","00200009387300008250400000","1111141057032230000","0h63103100000010"); // dEdx pi 2 high 1 (> 3.5 GeV)
  } else if (trainConfig == 105){ // PCM variations
    cuts.AddCutPCMCalo("84610113","00200009327000009250400000","1111141057032230000","0h63103100000010"); // qt 2D 0.03
    cuts.AddCutPCMCalo("84610113","00200009327000003250400000","1111141057032230000","0h63103100000010"); // qt 1D 0.05
    cuts.AddCutPCMCalo("84610113","00200009327000002250400000","1111141057032230000","0h63103100000010"); // qt 1D 0.07
    cuts.AddCutPCMCalo("84610113","00200049327000008250400000","1111141057032230000","0h63103100000010"); // single pt > 0.075
    cuts.AddCutPCMCalo("84610113","00200019327000008250400000","1111141057032230000","0h63103100000010"); // single pt > 0.1
  } else if (trainConfig == 106){ // PCM variations
    cuts.AddCutPCMCalo("84610113","00200009327000008850400000","1111141057032230000","0h63103100000010"); // 2D psi pair chi2 var
    cuts.AddCutPCMCalo("84610113","00200009327000008260400000","1111141057032230000","0h63103100000010"); // 2D psi pair chi2 var
    cuts.AddCutPCMCalo("84610113","00200009327000008860400000","1111141057032230000","0h63103100000010"); // 2D psi pair chi2 var
    cuts.AddCutPCMCalo("84610113","00200009327000008280400000","1111141057032230000","0h63103100000010"); // 2D psi pair chi2 var
    cuts.AddCutPCMCalo("84610113","00200009327000008880400000","1111141057032230000","0h63103100000010"); // 2D psi pair chi2 var
  } else if (trainConfig == 107){ // PCM variations
    cuts.AddCutPCMCalo("84610113","00200006327000008250400000","1111141057032230000","0h63103100000010"); // min TPC cl > 0.7
    cuts.AddCutPCMCalo("84610113","00200008327000008250400000","1111141057032230000","0h63103100000010"); // min TPC cl > 0.35
    cuts.AddCutPCMCalo("84610113","00200009f9730000dge0400000","1111141057032230000","0h63106100000010"); // alpha < 0.8
    cuts.AddCutPCMCalo("84610113","00200009f9730000dge0400000","1111141057032230000","0h63105100000010"); // alpha < 0.75
  } else if (trainConfig == 108){ // PCM variations
    cuts.AddCutPCMCalo("84610113","00202209327000008250400000","1111141057032230000","0h63103100000010"); // restrict acceptance to EMCAL loose
    cuts.AddCutPCMCalo("84610113","00204409327000008250400000","1111141057032230000","0h63103100000010"); // restrict acceptance to EMCAL tight
  } else if (trainConfig == 109){ // PCM variations pi dEdx
    cuts.AddCutPCMCalo("84610113","00200009317300008250400000","1111141057032230000","0h63103100000010"); // dEdx pi: 0: 0.4-3.5, -10: 3.5 ->
    cuts.AddCutPCMCalo("84610113","00200009327300008250400000","1111141057032230000","0h63103100000010"); // dEdx pi: 1: 0.4-3.5, -10: 3.5 ->
    cuts.AddCutPCMCalo("84610113","00200009325000008250400000","1111141057032230000","0h63103100000010"); // dEdx pi: 1: 0.3 ->
    cuts.AddCutPCMCalo("84610113","00200009320000008250400000","1111141057032230000","0h63103100000010"); // dEdx pi: 1: 0.5 ->
  } else if (trainConfig == 110){ // PCM variations pi dEdx
    cuts.AddCutPCMCalo("84610113","00200009327600008250400000","1111141057032230000","0h63103100000010"); // dEdx pi: 1: 0.4-2, -10: 2. ->
    cuts.AddCutPCMCalo("84610113","00200009327400008250400000","1111141057032230000","0h63103100000010"); // dEdx pi: 1: 0.4-3, -10: 3. ->
    cuts.AddCutPCMCalo("84610113","00200009315600008250400000","1111141057032230000","0h63103100000010"); // dEdx pi: 0: 0.3-2, -10: 2. ->
    cuts.AddCutPCMCalo("84610113","00200009367400008250400000","1111141057032230000","0h63103100000010"); // dEdx pi: 2: 0.4-3, 0.5: 3. ->
    cuts.AddCutPCMCalo("84610113","00200009347400008250400000","1111141057032230000","0h63103100000010"); // dEdx pi: 3: 0.4-3, 1: 3. ->
  } else if (trainConfig == 111){ // PCM variations to close V0s
    cuts.AddCutPCMCalo("84610113","00200009f9730000dge0400000","1111141057032230000","0h63103100000010");
    cuts.AddCutPCMCalo("84610113","00200009327000008250401000","1111141057032230000","0h63103100000010");
    cuts.AddCutPCMCalo("84610113","00200009327000008250402000","1111141057032230000","0h63103100000010");
    cuts.AddCutPCMCalo("84610113","00200009327000008250403000","1111141057032230000","0h63103100000010");
  } else if (trainConfig == 112){ // EMCal cluster, non lin variations
    cuts.AddCutPCMCalo("84610113","00200009f9730000dge0400000","1111142057032230000","0h63103100000010");
    cuts.AddCutPCMCalo("84610113","00200009f9730000dge0400000","1111143057032230000","0h63103100000010");
    cuts.AddCutPCMCalo("84610113","00200009f9730000dge0400000","1111151057032230000","0h63103100000010");
    cuts.AddCutPCMCalo("84610113","00200009f9730000dge0400000","1111152057032230000","0h63103100000010");
  } else if (trainConfig == 113){  // EMCal cluster, non lin variations
    cuts.AddCutPCMCalo("84610113","00200009f9730000dge0400000","1111100057032230000","0h63103100000010");
    cuts.AddCutPCMCalo("84610113","00200009f9730000dge0400000","1111101057032230000","0h63103100000010");
    cuts.AddCutPCMCalo("84610113","00200009f9730000dge0400000","1111102057032230000","0h63103100000010");
    cuts.AddCutPCMCalo("84610113","00200009f9730000dge0400000","1111144057032230000","0h63103100000010");

  //---------------------------------------------------------------------------------------------
  // 60-100 % variations
  //---------------------------------------------------------------------------------------------
  } else if (trainConfig == 120){ // EMCAL clusters standard cuts
    cuts.AddCutPCMCalo("86010113","00200009f9730000dge0400000","1111141057032230000","0h63103100000010"); // INT7
    cuts.AddCutPCMCalo("86052113","00200009f9730000dge0400000","1111141057032230000","0h63103100000010"); // EMC7
    cuts.AddCutPCMCalo("86085113","00200009f9730000dge0400000","1111141057022230000","0h63103100000010"); // EG2
    cuts.AddCutPCMCalo("86083113","00200009f9730000dge0400000","1111141057022230000","0h63103100000010"); // EG1
  } else if (trainConfig == 121){ //EMCAL minEnergy variation
    cuts.AddCutPCMCalo("86010113","00200009f9730000dge0400000","1111141057032230000","0h63103100000010"); // 0.7 GeV/c default
    cuts.AddCutPCMCalo("86010113","00200009f9730000dge0400000","1111141057042230000","0h63103100000010"); // 0.8 GeV/c
    cuts.AddCutPCMCalo("86010113","00200009f9730000dge0400000","1111141057052230000","0h63103100000010"); // 0.9 GeV/c
  } else if (trainConfig == 122){ //EMCAL minNCells variation
    cuts.AddCutPCMCalo("86010113","00200009f9730000dge0400000","1111141057031230000","0h63103100000010"); // n cells >= 1
    cuts.AddCutPCMCalo("86010113","00200009f9730000dge0400000","1111141057033230000","0h63103100000010"); // n cells >= 3
    cuts.AddCutPCMCalo("86010113","00200009f9730000dge0400000","1111141057032200000","0h63103100000010"); // no max M02 cut
    cuts.AddCutPCMCalo("86010113","00200009f9730000dge0400000","1111141057032250000","0h63103100000010"); // M02 < 0.3
    cuts.AddCutPCMCalo("86010113","00200009f9730000dge0400000","1111141057032260000","0h63103100000010"); // M02 < 0.27
    cuts.AddCutPCMCalo("86010113","00200009f9730000dge0400000","1112141057032230000","0h63103100000010"); // only modules with TRD infront
    cuts.AddCutPCMCalo("86010113","00200009f9730000dge0400000","1111341057032230000","0h63103100000010"); // no modules with TRD infront
  } else if (trainConfig == 123){ // EMCAL track matching variations
    cuts.AddCutPCMCalo("86010113","00200009f9730000dge0400000","1111141053032230000","0h63103100000010"); //
    cuts.AddCutPCMCalo("86010113","00200009f9730000dge0400000","1111141056032230000","0h63103100000010"); //
    cuts.AddCutPCMCalo("86010113","00200009f9730000dge0400000","1111141058032230000","0h63103100000010"); //
    cuts.AddCutPCMCalo("86010113","00200009f9730000dge0400000","1111141059032230000","0h63103100000010"); //
  } else if (trainConfig == 124){ // PCM variations
    cuts.AddCutPCMCalo("86010113","00200009227000008250400000","1111141057032230000","0h63103100000010"); // dEdx e -3, 5
    cuts.AddCutPCMCalo("86010113","00200009127000008250400000","1111141057032230000","0h63103100000010"); // dEdx e -5, 5
    cuts.AddCutPCMCalo("86010113","00200009357000008250400000","1111141057032230000","0h63103100000010"); // dEdx pi 2
    cuts.AddCutPCMCalo("86010113","00200009317000008250400000","1111141057032230000","0h63103100000010"); // dEdx pi 0
    cuts.AddCutPCMCalo("86010113","00200009387300008250400000","1111141057032230000","0h63103100000010"); // dEdx pi 2 high 1 (> 3.5 GeV)
  } else if (trainConfig == 125){ // PCM variations
    cuts.AddCutPCMCalo("86010113","00200009327000009250400000","1111141057032230000","0h63103100000010"); // qt 2D 0.03
    cuts.AddCutPCMCalo("86010113","00200009327000003250400000","1111141057032230000","0h63103100000010"); // qt 1D 0.05
    cuts.AddCutPCMCalo("86010113","00200009327000002250400000","1111141057032230000","0h63103100000010"); // qt 1D 0.07
    cuts.AddCutPCMCalo("86010113","00200049327000008250400000","1111141057032230000","0h63103100000010"); // single pt > 0.075
    cuts.AddCutPCMCalo("86010113","00200019327000008250400000","1111141057032230000","0h63103100000010"); // single pt > 0.1
  } else if (trainConfig == 126){ // PCM variations
    cuts.AddCutPCMCalo("86010113","00200009327000008850400000","1111141057032230000","0h63103100000010"); // 2D psi pair chi2 var
    cuts.AddCutPCMCalo("86010113","00200009327000008260400000","1111141057032230000","0h63103100000010"); // 2D psi pair chi2 var
    cuts.AddCutPCMCalo("86010113","00200009327000008860400000","1111141057032230000","0h63103100000010"); // 2D psi pair chi2 var
    cuts.AddCutPCMCalo("86010113","00200009327000008280400000","1111141057032230000","0h63103100000010"); // 2D psi pair chi2 var
    cuts.AddCutPCMCalo("86010113","00200009327000008880400000","1111141057032230000","0h63103100000010"); // 2D psi pair chi2 var
  } else if (trainConfig == 127){ // PCM variations
    cuts.AddCutPCMCalo("86010113","00200006327000008250400000","1111141057032230000","0h63103100000010"); // min TPC cl > 0.7
    cuts.AddCutPCMCalo("86010113","00200008327000008250400000","1111141057032230000","0h63103100000010"); // min TPC cl > 0.35
    cuts.AddCutPCMCalo("86010113","00200009f9730000dge0400000","1111141057032230000","0163106100000010"); // alpha < 0.8
    cuts.AddCutPCMCalo("86010113","00200009f9730000dge0400000","1111141057032230000","0163105100000010"); // alpha < 0.75
  } else if (trainConfig == 128){ // PCM variations
    cuts.AddCutPCMCalo("86010113","00202209327000008250400000","1111141057032230000","0h63103100000010"); // restrict acceptance to EMCAL loose
    cuts.AddCutPCMCalo("86010113","00204409327000008250400000","1111141057032230000","0h63103100000010"); // restrict acceptance to EMCAL tight
  } else if (trainConfig == 129){ // PCM variations pi dEdx
    cuts.AddCutPCMCalo("86010113","00200009317300008250400000","1111141057032230000","0h63103100000010"); // dEdx pi: 0: 0.4-3.5, -10: 3.5 ->
    cuts.AddCutPCMCalo("86010113","00200009327300008250400000","1111141057032230000","0h63103100000010"); // dEdx pi: 1: 0.4-3.5, -10: 3.5 ->
    cuts.AddCutPCMCalo("86010113","00200009325000008250400000","1111141057032230000","0h63103100000010"); // dEdx pi: 1: 0.3 ->
    cuts.AddCutPCMCalo("86010113","00200009320000008250400000","1111141057032230000","0h63103100000010"); // dEdx pi: 1: 0.5 ->
  } else if (trainConfig == 130){ // PCM variations pi dEdx
    cuts.AddCutPCMCalo("86010113","00200009327600008250400000","1111141057032230000","0h63103100000010"); // dEdx pi: 1: 0.4-2, -10: 2. ->
    cuts.AddCutPCMCalo("86010113","00200009327400008250400000","1111141057032230000","0h63103100000010"); // dEdx pi: 1: 0.4-3, -10: 3. ->
    cuts.AddCutPCMCalo("86010113","00200009315600008250400000","1111141057032230000","0h63103100000010"); // dEdx pi: 0: 0.3-2, -10: 2. ->
    cuts.AddCutPCMCalo("86010113","00200009367400008250400000","1111141057032230000","0h63103100000010"); // dEdx pi: 2: 0.4-3, 0.5: 3. ->
    cuts.AddCutPCMCalo("86010113","00200009347400008250400000","1111141057032230000","0h63103100000010"); // dEdx pi: 3: 0.4-3, 1: 3. ->
  } else if (trainConfig == 131){ // PCM variations to close V0s
    cuts.AddCutPCMCalo("86010113","00200009f9730000dge0400000","1111141057032230000","0h63103100000010");
    cuts.AddCutPCMCalo("86010113","00200009327000008250401000","1111141057032230000","0h63103100000010");
    cuts.AddCutPCMCalo("86010113","00200009327000008250402000","1111141057032230000","0h63103100000010");
    cuts.AddCutPCMCalo("86010113","00200009327000008250403000","1111141057032230000","0h63103100000010");
  } else if (trainConfig == 132){ // EMCal cluster, non lin variations
    cuts.AddCutPCMCalo("86010113","00200009f9730000dge0400000","1111142057032230000","0h63103100000010");
    cuts.AddCutPCMCalo("86010113","00200009f9730000dge0400000","1111143057032230000","0h63103100000010");
    cuts.AddCutPCMCalo("86010113","00200009f9730000dge0400000","1111151057032230000","0h63103100000010");
    cuts.AddCutPCMCalo("86010113","00200009f9730000dge0400000","1111152057032230000","0h63103100000010");
  } else if (trainConfig == 133){  // EMCal cluster, non lin variations
    cuts.AddCutPCMCalo("86010113","00200009f9730000dge0400000","1111100057032230000","0h63103100000010");
    cuts.AddCutPCMCalo("86010113","00200009f9730000dge0400000","1111101057032230000","0h63103100000010");
    cuts.AddCutPCMCalo("86010113","00200009f9730000dge0400000","1111102057032230000","0h63103100000010");
    cuts.AddCutPCMCalo("86010113","00200009f9730000dge0400000","1111144057032230000","0h63103100000010");

  // ===============================================================================================
  // Run 2 data EMC clusters pPb 5TeV
  // ===============================================================================================
  } else if (trainConfig == 200){ // EMCAL clusters standard cuts, triggers, no nonlin, open timing
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","1111100057032230000","0h63103100000010"); // INT7
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","1111100017032230000","0h63103100000010"); // INT7
  } else if (trainConfig == 201){
    cuts.AddCutPCMCalo("80052113","00200009f9730000dge0400000","1111100057032230000","0h63103100000010"); // EMC7
    cuts.AddCutPCMCalo("80085113","00200009f9730000dge0400000","1111100057032230000","0h63103100000010"); // EG2
    cuts.AddCutPCMCalo("80083113","00200009f9730000dge0400000","1111100057032230000","0h63103100000010"); // EG1
  } else if (trainConfig == 202){
    cuts.AddCutPCMCalo("80052113","00200009f9730000dge0400000","1111100017032230000","0h63103100000010"); // EMC7
    cuts.AddCutPCMCalo("80085113","00200009f9730000dge0400000","1111100017032230000","0h63103100000010"); // EG2
    cuts.AddCutPCMCalo("80083113","00200009f9730000dge0400000","1111100017032230000","0h63103100000010"); // EG1
  } else if (trainConfig == 203){ // EMCAL clusters standard cuts, triggers, no nonlin, open timing
    cuts.AddCutPCMCalo("80110113","00200009f9730000dge0400000","1111100057032230000","0h63103100000010"); // 0-10
    cuts.AddCutPCMCalo("81210113","00200009f9730000dge0400000","1111100057032230000","0h63103100000010"); // 10-20
    cuts.AddCutPCMCalo("82410113","00200009f9730000dge0400000","1111100057032230000","0h63103100000010"); // 20-40
    cuts.AddCutPCMCalo("84610113","00200009f9730000dge0400000","1111100057032230000","0h63103100000010"); // 40-60
    cuts.AddCutPCMCalo("86810113","00200009f9730000dge0400000","1111100057032230000","0h63103100000010"); // 60-80
    cuts.AddCutPCMCalo("88010113","00200009f9730000dge0400000","1111100057032230000","0h63103100000010"); // 80-100
  } else if (trainConfig == 204){ // EMCAL clusters standard cuts, triggers, no nonlin, open timing
    cuts.AddCutPCMCalo("80210113","00200009f9730000dge0400000","1111100057032230000","0h63103100000010"); // 0-20
    cuts.AddCutPCMCalo("86010113","00200009f9730000dge0400000","1111100057032230000","0h63103100000010"); // 60-100
    cuts.AddCutPCMCalo("a0110113","00200009f9730000dge0400000","1111100057032230000","0h63103100000010"); // 0-5
    cuts.AddCutPCMCalo("a1210113","00200009f9730000dge0400000","1111100057032230000","0h63103100000010"); // 5-10
  } else if (trainConfig == 205){ // EMCAL clusters standard cuts, triggers, no nonlin, open timing
    cuts.AddCutPCMCalo("80110113","00200009f9730000dge0400000","1111100017032230000","0h63103100000010"); // 0-10
    cuts.AddCutPCMCalo("81210113","00200009f9730000dge0400000","1111100017032230000","0h63103100000010"); // 10-20
    cuts.AddCutPCMCalo("82410113","00200009f9730000dge0400000","1111100017032230000","0h63103100000010"); // 20-40
    cuts.AddCutPCMCalo("84610113","00200009f9730000dge0400000","1111100017032230000","0h63103100000010"); // 40-60
    cuts.AddCutPCMCalo("86810113","00200009f9730000dge0400000","1111100017032230000","0h63103100000010"); // 60-80
    cuts.AddCutPCMCalo("88010113","00200009f9730000dge0400000","1111100017032230000","0h63103100000010"); // 80-100
  } else if (trainConfig == 206){ // EMCAL clusters standard cuts, triggers, no nonlin, open timing
    cuts.AddCutPCMCalo("80210113","00200009f9730000dge0400000","1111100017032230000","0h63103100000010"); // 0-20
    cuts.AddCutPCMCalo("86010113","00200009f9730000dge0400000","1111100017032230000","0h63103100000010"); // 60-100
    cuts.AddCutPCMCalo("a0110113","00200009f9730000dge0400000","1111100017032230000","0h63103100000010"); // 0-5
    cuts.AddCutPCMCalo("a1210113","00200009f9730000dge0400000","1111100017032230000","0h63103100000010"); // 5-10

  } else if (trainConfig == 207){ // EMCAL clusters standard cuts, nonlin variations
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","1111141057032230000","0h63103100000010"); // 0-100
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","1111151057032230000","0h63103100000010"); // 0-100
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","1111142057032230000","0h63103100000010"); // 0-100
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","1111152057032230000","0h63103100000010"); // 0-100
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","1111155057032230000","0h63103100000010"); // 0-100
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","1111156057032230000","0h63103100000010"); // 0-100
  } else if (trainConfig == 208){ // EMCAL clusters standard cuts, nonlin variations
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","1111100057032230000","0h63103100000010"); // 0-100
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","1111102057032230000","0h63103100000010"); // 0-100
  } else if (trainConfig == 209){ // EMCAL clusters standard cuts, nonlin variations
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","1111141017032230000","0h63103100000010"); // 0-100
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","1111151017032230000","0h63103100000010"); // 0-100
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","1111142017032230000","0h63103100000010"); // 0-100
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","1111152017032230000","0h63103100000010"); // 0-100
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","1111155017032230000","0h63103100000010"); // 0-100
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","1111156017032230000","0h63103100000010"); // 0-100
  } else if (trainConfig == 210){ // EMCAL clusters standard cuts, nonlin variations
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","1111100017032230000","0h63103100000010"); // 0-100
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","1111102017032230000","0h63103100000010"); // 0-100

  // non lin variations with cent
  } else if (trainConfig == 211){ // EMCAL clusters standard cuts, triggers, no nonlin, open timing
    cuts.AddCutPCMCalo("80110113","00200009f9730000dge0400000","1111141057032230000","0h63103100000010"); // 0-10
    cuts.AddCutPCMCalo("81210113","00200009f9730000dge0400000","1111141057032230000","0h63103100000010"); // 10-20
    cuts.AddCutPCMCalo("82410113","00200009f9730000dge0400000","1111141057032230000","0h63103100000010"); // 20-40
    cuts.AddCutPCMCalo("84610113","00200009f9730000dge0400000","1111141057032230000","0h63103100000010"); // 40-60
    cuts.AddCutPCMCalo("86810113","00200009f9730000dge0400000","1111141057032230000","0h63103100000010"); // 60-80
    cuts.AddCutPCMCalo("88010113","00200009f9730000dge0400000","1111141057032230000","0h63103100000010"); // 80-100
  } else if (trainConfig == 212){ // EMCAL clusters standard cuts, triggers, no nonlin, open timing
    cuts.AddCutPCMCalo("80210113","00200009f9730000dge0400000","1111141057032230000","0h63103100000010"); // 0-20
    cuts.AddCutPCMCalo("86010113","00200009f9730000dge0400000","1111141057032230000","0h63103100000010"); // 60-100
    cuts.AddCutPCMCalo("a0110113","00200009f9730000dge0400000","1111141057032230000","0h63103100000010"); // 0-5
    cuts.AddCutPCMCalo("a1210113","00200009f9730000dge0400000","1111141057032230000","0h63103100000010"); // 5-10
  } else if (trainConfig == 213){ // EMCAL clusters standard cuts, triggers, no nonlin, open timing
    cuts.AddCutPCMCalo("80110113","00200009f9730000dge0400000","1111151057032230000","0h63103100000010"); // 0-10
    cuts.AddCutPCMCalo("81210113","00200009f9730000dge0400000","1111151057032230000","0h63103100000010"); // 10-20
    cuts.AddCutPCMCalo("82410113","00200009f9730000dge0400000","1111151057032230000","0h63103100000010"); // 20-40
    cuts.AddCutPCMCalo("84610113","00200009f9730000dge0400000","1111151057032230000","0h63103100000010"); // 40-60
    cuts.AddCutPCMCalo("86810113","00200009f9730000dge0400000","1111151057032230000","0h63103100000010"); // 60-80
    cuts.AddCutPCMCalo("88010113","00200009f9730000dge0400000","1111151057032230000","0h63103100000010"); // 80-100
  } else if (trainConfig == 214){ // EMCAL clusters standard cuts, triggers, no nonlin, open timing
    cuts.AddCutPCMCalo("80210113","00200009f9730000dge0400000","1111151057032230000","0h63103100000010"); // 0-20
    cuts.AddCutPCMCalo("86010113","00200009f9730000dge0400000","1111151057032230000","0h63103100000010"); // 60-100
    cuts.AddCutPCMCalo("a0110113","00200009f9730000dge0400000","1111151057032230000","0h63103100000010"); // 0-5
    cuts.AddCutPCMCalo("a1210113","00200009f9730000dge0400000","1111151057032230000","0h63103100000010"); // 5-10

  // non lin variations with diff minE cut
  } else if (trainConfig == 215){ // EMCAL clusters nonlin variations, min E=0.6 GeV
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","1111100057022230000","0h63103100000010"); // 0-100
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","1111141057022230000","0h63103100000010"); // 0-100
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","1111151057022230000","0h63103100000010"); // 0-100
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","1111142057022230000","0h63103100000010"); // 0-100
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","1111152057022230000","0h63103100000010"); // 0-100
  } else if (trainConfig == 216){ // EMCAL clusters nonlin variations, min E=0.65 GeV
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","11111000570b2230000","0h63103100000010"); // 0-100
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","11111410570b2230000","0h63103100000010"); // 0-100
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","11111510570b2230000","0h63103100000010"); // 0-100
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","11111420570b2230000","0h63103100000010"); // 0-100
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","11111520570b2230000","0h63103100000010"); // 0-100
  } else if (trainConfig == 217){ // EMCAL clusters nonlin variations, min E=0.675 GeV
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","11111000570c2230000","0h63103100000010"); // 0-100
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","11111410570c2230000","0h63103100000010"); // 0-100
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","11111510570c2230000","0h63103100000010"); // 0-100
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","11111420570c2230000","0h63103100000010"); // 0-100
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","11111520570c2230000","0h63103100000010"); // 0-100
  } else if (trainConfig == 218){ // EMCAL clusters nonlin variations, min E=0.625 GeV
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","11111000570d2230000","0h63103100000010"); // 0-100
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","11111410570d2230000","0h63103100000010"); // 0-100
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","11111510570d2230000","0h63103100000010"); // 0-100
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","11111420570d2230000","0h63103100000010"); // 0-100
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","11111520570d2230000","0h63103100000010"); // 0-100

  // narrow cent variations
  } else if (trainConfig == 219){
    cuts.AddCutPCMCalo("c0110113","00200009f9730000dge0400000","1111151057032230000","0h63103100000010"); // 0-1
    cuts.AddCutPCMCalo("c0210113","00200009f9730000dge0400000","1111151057032230000","0h63103100000010"); // 0-2
  } else if (trainConfig == 220){
    cuts.AddCutPCMCalo("c0110113","00200009f9730000dge0400000","1111141057032230000","0h63103100000010"); // 0-1
    cuts.AddCutPCMCalo("c0210113","00200009f9730000dge0400000","1111141057032230000","0h63103100000010"); // 0-2
  // AOD validation
  } else if (trainConfig == 221){
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","1111151057032230000","0h63103100000010"); // 0-100

  // ===============================================================================================
  // Run 1 data PHOS clusters
  // ===============================================================================================
  } else if (trainConfig == 301) {  // min energy = 0.3 GeV/c
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","2444400041013200000","0h63103100000010"); // standart cut, kINT7
    cuts.AddCutPCMCalo("80062113","00200009f9730000dge0400000","2444400041013200000","0h63103100000010"); // standard cut, kPHI7
  } else if (trainConfig == 302) { //PHOS
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","2444400041013200000","0h63103100000010");
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","2444400042013200000","0h63103100000010");
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","2444400043013200000","0h63103100000010");
  } else if (trainConfig == 303) { //PHOS
    cuts.AddCutPCMCalo("80010023","00200009f9730000dge0400000","2444400041013200000","0h63103100000010");
    cuts.AddCutPCMCalo("80010023","00200009f9730000dge0400000","2444400042013200000","0h63103100000010");
    cuts.AddCutPCMCalo("80010023","00200009f9730000dge0400000","2444400043013200000","0h63103100000010");
  } else if (trainConfig == 305) { //PHOS timing cut variations
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","2444400011013200000","0h63103100000010"); // 1000ns
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","2444400031013200000","0h63103100000010"); // 200ns
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","2444400051013200000","0h63103100000010"); // 50ns
  } else if (trainConfig == 306) {  // Non lin variations PHOS INT7
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","2444401041013200000","0h63103100000010"); // PHOS group default
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","2444451041013200000","0h63103100000010"); // CCMF PHOS
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","2444452041013200000","0h63103100000010"); // CMF PHOS
  } else if (trainConfig == 307) {  // Non lin variations PHOS PHI7
    cuts.AddCutPCMCalo("80062113","00200009f9730000dge0400000","2444401041013200000","0h63103100000010"); // PHOS group default
    cuts.AddCutPCMCalo("80062113","00200009f9730000dge0400000","2444451041013200000","0h63103100000010"); // CCMF PHOS
    cuts.AddCutPCMCalo("80062113","00200009f9730000dge0400000","2444452041013200000","0h63103100000010"); // CMF PHOS
  } else if (trainConfig == 308) {  // CCMF cent dependet
    cuts.AddCutPCMCalo("80210113","00200009f9730000dge0400000","2444451041013200000","0h63103100000010"); // 0-20
    cuts.AddCutPCMCalo("82410113","00200009f9730000dge0400000","2444451041013200000","0h63103100000010"); // 20-40
    cuts.AddCutPCMCalo("84610113","00200009f9730000dge0400000","2444451041013200000","0h63103100000010"); // 40-60
    cuts.AddCutPCMCalo("86010113","00200009f9730000dge0400000","2444451041013200000","0h63103100000010"); // 60-100
  } else if (trainConfig == 309) {  // PHOS default cent dependet
    cuts.AddCutPCMCalo("80210113","00200009f9730000dge0400000","2444401041013200000","0h63103100000010"); // 0-20
    cuts.AddCutPCMCalo("82410113","00200009f9730000dge0400000","2444401041013200000","0h63103100000010"); // 20-40
    cuts.AddCutPCMCalo("84610113","00200009f9730000dge0400000","2444401041013200000","0h63103100000010"); // 40-60
    cuts.AddCutPCMCalo("86010113","00200009f9730000dge0400000","2444401041013200000","0h63103100000010"); // 60-100

  // testing past future protection
  } else if (trainConfig == 310) {  // PHOS INT7 PF test
    cuts.AddCutPCMCalo("80010213","00200009f9730000dge0400000","2444451041013200000","0h63103100000010"); // PF 2.25 \mu
    cuts.AddCutPCMCalo("80010513","00200009f9730000dge0400000","2444451041013200000","0h63103100000010"); // PF 1.075 \mu
    cuts.AddCutPCMCalo("80010413","00200009f9730000dge0400000","2444451041013200000","0h63103100000010"); // PF 0.25 \mu
  } else if (trainConfig == 311) {  // PHOS INT7 PF test, open cluster time
    cuts.AddCutPCMCalo("80010213","00200009f9730000dge0400000","2444451011013200000","0h63103100000010"); // PF 2.25 \mu
    cuts.AddCutPCMCalo("80010513","00200009f9730000dge0400000","2444451011013200000","0h63103100000010"); // PF 1.075 \mu
    cuts.AddCutPCMCalo("80010413","00200009f9730000dge0400000","2444451011013200000","0h63103100000010"); // PF 0.25 \mu
  } else if (trainConfig == 312) {  // PHOS PHI7 PF test
    cuts.AddCutPCMCalo("80062213","00200009f9730000dge0400000","2444451041013200000","0h63103100000010"); // PF 2.25 \mu
    cuts.AddCutPCMCalo("80062513","00200009f9730000dge0400000","2444451041013200000","0h63103100000010"); // PF 1.075 \mu
    cuts.AddCutPCMCalo("80062413","00200009f9730000dge0400000","2444451041013200000","0h63103100000010"); // PF 0.25 \mu
  } else if (trainConfig == 313) {  // PHOS PHI7 PF test, open cluster time
    cuts.AddCutPCMCalo("80062213","00200009f9730000dge0400000","2444451011013200000","0h63103100000010"); // PF 2.25 \mu
    cuts.AddCutPCMCalo("80062513","00200009f9730000dge0400000","2444451011013200000","0h63103100000010"); // PF 1.075 \mu
    cuts.AddCutPCMCalo("80062413","00200009f9730000dge0400000","2444451011013200000","0h63103100000010"); // PF 0.25 \mu

    //PCM-PHOS systematics
  } else if(trainConfig == 320){//dEdx e-line variation and dE/dx pi-line variation
    cuts.AddCutPCMCalo("80010113","00200009127000008250400000","2444451041013200000","0h63103100000010"); //-5 < sigma < 5
    cuts.AddCutPCMCalo("80010113","00200009227000008250400000","2444451041013200000","0h63103100000010"); //-3 < sigma < 5
    cuts.AddCutPCMCalo("80010113","00200009327400008250400000","2444451041013200000","0h63103100000010"); //1, -10, 0.4, 3
    cuts.AddCutPCMCalo("80010113","00200009367400008250400000","2444451041013200000","0h63103100000010"); //2, 0.5, 0.4, 3
  } else if(trainConfig == 321){//dE/dx pi-line variation and single pT variation and 2D triangular chi2 and psi pair
    cuts.AddCutPCMCalo("80010113","00200009317400008250400000","2444451041013200000","0h63103100000010"); //0, -10, 0.4, 3
    cuts.AddCutPCMCalo("80010113","00200049327000008250400000","2444451041013200000","0h63103100000010"); //single pT 0.075 GeV/c
    cuts.AddCutPCMCalo("80010113","00200019327000008250400000","2444451041013200000","0h63103100000010"); //single pT 0.1   GeV/c
    cuts.AddCutPCMCalo("80010113","00200009327000008850400000","2444451041013200000","0h63103100000010"); //20 & 0.1
  } else if(trainConfig == 322){//2D triangular chi2 and psi pair
    cuts.AddCutPCMCalo("80010113","00200009327000008260400000","2444451041013200000","0h63103100000010"); //30 & 0.05
    cuts.AddCutPCMCalo("80010113","00200009327000008860400000","2444451041013200000","0h63103100000010"); //20 & 0.05
    cuts.AddCutPCMCalo("80010113","00200009327000008280400000","2444451041013200000","0h63103100000010"); //30 & 0.2
    cuts.AddCutPCMCalo("80010113","00200009327000008880400000","2444451041013200000","0h63103100000010"); //20 & 0.2
  } else if(trainConfig == 323){//min TPC clusters and //elliptic cut Qt and alpha
    cuts.AddCutPCMCalo("80010113","00200000327000008250400000","2444451041013200000","0h63103100000010"); //0
    cuts.AddCutPCMCalo("80010113","00200008327000008250400000","2444451041013200000","0h63103100000010"); //0.35
    cuts.AddCutPCMCalo("80010113","00200009327000009250400000","2444451041013200000","0h63103100000010"); // qT 0.03 no quadratic
    cuts.AddCutPCMCalo("80010113","00200009327000003250400000","2444451041013200000","0h63103100000010"); // qT 0.05 y  quadratic
  } else if(trainConfig == 324){// and min phi max phi and NL variations
    cuts.AddCutPCMCalo("80010113","00209909327000008250400000","2444451041013200000","0h63103100000010"); //4.54 > phi > 5.59
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","2444452041013200000","0h63103100000010"); // PHOS calo NL
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","2444401041013200000","0h63103100000010"); // PHOS people NL

    //Cluster Cuts varation
  } else if(trainConfig == 325){ // first set of variations CLUSTER
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","2444451041073200000","0h63103100000010"); // min energy 0.2 GeV
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","2444451041083200000","0h63103100000010"); // min energy 0.4 GeV
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","2444451041023200000","0h63103100000010"); // min energy 0.5 GeV
  } else if(trainConfig == 326){ // second set of variations CLUSTER
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","2444451041013270000","0h63103100000010"); // min/max M02  0.1<M<1
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","2444451041013280000","0h63103100000010"); // min/max M02  0.1<M<1
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","2444451041012200000","0h63103100000010"); // min number 2 cells
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","2444451041014200000","0h63103100000010"); // min number 4 cells
  } else if(trainConfig == 327){ // MESON
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","2444451041013200000","0h63403100000010"); // rapidity variation  y<0.5
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","2444451041013200000","0h63803100000010"); // rapidity variation  y<0.25
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","2444451041013200000","0h63106100000010"); // alpha meson variation 1 0<alpha<0.8
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","2444451041013200000","0h63105100000010"); // alpha meson variation 2 0<alpha<0.75
  } else if(trainConfig == 328){ // fourth set of variations
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","2444451044013200000","0h63103100000010"); // tm variation
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","2444451045013200000","0h63103100000010"); // tm variation
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","2444451046013200000","0h63103100000010"); // tm variation
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","2444451041013200000","0h63103100000000"); // min opening angle 0    -> open
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","2444451041013200000","0h63103100000030"); // min opening angle 0.01 -> 2 cell diag
  } else if(trainConfig == 329){ // centrality dependent and with latest TM
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","2444451044013200000","0h63103100000010"); // 0-100% with PCM NL
    cuts.AddCutPCMCalo("80210113","00200009f9730000dge0400000","2444451044013200000","0h63103100000010"); // 0-20% with PCM NL
    cuts.AddCutPCMCalo("82410113","00200009f9730000dge0400000","2444451044013200000","0h63103100000010"); // 20-40% with PCM NL
    cuts.AddCutPCMCalo("84610113","00200009f9730000dge0400000","2444451044013200000","0h63103100000010"); // 40-60% with PCM NL
    cuts.AddCutPCMCalo("86010113","00200009f9730000dge0400000","2444451044013200000","0h63103100000010"); // 60-100% with PCM NL

    //-----------------------------------------------------------------------------------------------
    // PCM-PHOS systematics run 1
    //-----------------------------------------------------------------------------------------------
    //0 - 20
  } else if(trainConfig == 330){//dEdx e-line variation and dE/dx pi-line variation
    cuts.AddCutPCMCalo("80210113","00200009127000008250400000","2444451041013200000","0h63103100000010"); //-5 < sigma < 5
    cuts.AddCutPCMCalo("80210113","00200009227000008250400000","2444451041013200000","0h63103100000010"); //-3 < sigma < 5
    cuts.AddCutPCMCalo("80210113","00200009327400008250400000","2444451041013200000","0h63103100000010"); //1, -10, 0.4, 3
    cuts.AddCutPCMCalo("80210113","00200009367400008250400000","2444451041013200000","0h63103100000010"); //2, 0.5, 0.4, 3
  } else if(trainConfig == 331){//dE/dx pi-line variation and single pT variation and 2D triangular chi2 and psi pair
    cuts.AddCutPCMCalo("80210113","00200009317400008250400000","2444451041013200000","0h63103100000010"); //0, -10, 0.4, 3
    cuts.AddCutPCMCalo("80210113","00200049327000008250400000","2444451041013200000","0h63103100000010"); //single pT 0.075 GeV/c
    cuts.AddCutPCMCalo("80210113","00200019327000008250400000","2444451041013200000","0h63103100000010"); //single pT 0.1   GeV/c
    cuts.AddCutPCMCalo("80210113","00200009327000008850400000","2444451041013200000","0h63103100000010"); //20 & 0.1
  } else if(trainConfig == 332){//2D triangular chi2 and psi pair
    cuts.AddCutPCMCalo("80210113","00200009327000008260400000","2444451041013200000","0h63103100000010"); //30 & 0.05
    cuts.AddCutPCMCalo("80210113","00200009327000008860400000","2444451041013200000","0h63103100000010"); //20 & 0.05
    cuts.AddCutPCMCalo("80210113","00200009327000008280400000","2444451041013200000","0h63103100000010"); //30 & 0.2
    cuts.AddCutPCMCalo("80210113","00200009327000008880400000","2444451041013200000","0h63103100000010"); //20 & 0.2
  } else if(trainConfig == 333){//min TPC clusters and //elliptic cut Qt and alpha
    cuts.AddCutPCMCalo("80210113","00200000327000008250400000","2444451041013200000","0h63103100000010"); //0
    cuts.AddCutPCMCalo("80210113","00200008327000008250400000","2444451041013200000","0h63103100000010"); //0.35
    cuts.AddCutPCMCalo("80210113","00200009327000009250400000","2444451041013200000","0h63103100000010"); // qT 0.03 no quadratic
    cuts.AddCutPCMCalo("80210113","00200009327000003250400000","2444451041013200000","0h63103100000010"); // qT 0.05 y  quadratic
  } else if(trainConfig == 333){// and min phi max phi and NL variations
    cuts.AddCutPCMCalo("80210113","00209909327000008250400000","2444451041013200000","0h63103100000010"); //4.54 > phi > 5.59
    cuts.AddCutPCMCalo("80210113","00200009f9730000dge0400000","2444452041013200000","0h63103100000010"); // PHOS calo NL
    cuts.AddCutPCMCalo("80210113","00200009f9730000dge0400000","2444401041013200000","0h63103100000010"); // PHOS people NL

    //Cluster Cuts varation
  } else if(trainConfig == 334){ // first set of variations CLUSTER
    cuts.AddCutPCMCalo("80210113","00200009f9730000dge0400000","2444451041073200000","0h63103100000010"); // min energy 0.2 GeV
    cuts.AddCutPCMCalo("80210113","00200009f9730000dge0400000","2444451041083200000","0h63103100000010"); // min energy 0.4 GeV
    cuts.AddCutPCMCalo("80210113","00200009f9730000dge0400000","2444451041023200000","0h63103100000010"); // min energy 0.5 GeV
  } else if(trainConfig == 335){ // second set of variations CLUSTER
    cuts.AddCutPCMCalo("80210113","00200009f9730000dge0400000","2444451041013270000","0h63103100000010"); // min/max M02  0.1<M<1
    cuts.AddCutPCMCalo("80210113","00200009f9730000dge0400000","2444451041013280000","0h63103100000010"); // min/max M02  0.1<M<1
    cuts.AddCutPCMCalo("80210113","00200009f9730000dge0400000","2444451041012200000","0h63103100000010"); // min number 2 cells
    cuts.AddCutPCMCalo("80210113","00200009f9730000dge0400000","2444451041014200000","0h63103100000010"); // min number 4 cells
  } else if(trainConfig == 336){ // MESON
    cuts.AddCutPCMCalo("80210113","00200009f9730000dge0400000","2444451041013200000","0h63403100000010"); // rapidity variation  y<0.5
    cuts.AddCutPCMCalo("80210113","00200009f9730000dge0400000","2444451041013200000","0h63803100000010"); // rapidity variation  y<0.25
    cuts.AddCutPCMCalo("80210113","00200009f9730000dge0400000","2444451041013200000","0h63106100000010"); // alpha meson variation 1 0<alpha<0.8
    cuts.AddCutPCMCalo("80210113","00200009f9730000dge0400000","2444451041013200000","0h63105100000010"); // alpha meson variation 2 0<alpha<0.75
  } else if(trainConfig == 337){ // fourth set of variations
    cuts.AddCutPCMCalo("80210113","00200009f9730000dge0400000","2444451044013200000","0h63103100000010"); // tm variation
    cuts.AddCutPCMCalo("80210113","00200009f9730000dge0400000","2444451045013200000","0h63103100000010"); // tm variation
    cuts.AddCutPCMCalo("80210113","00200009f9730000dge0400000","2444451046013200000","0h63103100000010"); // tm variation
    cuts.AddCutPCMCalo("80210113","00200009f9730000dge0400000","2444451041013200000","0h63103100000000"); // min opening angle 0    -> open
    cuts.AddCutPCMCalo("80210113","00200009f9730000dge0400000","2444451041013200000","0h63103100000030"); // min opening angle 0.01 -> 2 cell diag
    //20 - 40
  } else if(trainConfig == 340){//dEdx e-line variation and dE/dx pi-line variation
    cuts.AddCutPCMCalo("82410113","00200009127000008250400000","2444451041013200000","0h63103100000010"); //-5 < sigma < 5
    cuts.AddCutPCMCalo("82410113","00200009227000008250400000","2444451041013200000","0h63103100000010"); //-3 < sigma < 5
    cuts.AddCutPCMCalo("82410113","00200009327400008250400000","2444451041013200000","0h63103100000010"); //1, -10, 0.4, 3
    cuts.AddCutPCMCalo("82410113","00200009367400008250400000","2444451041013200000","0h63103100000010"); //2, 0.5, 0.4, 3
  } else if(trainConfig == 341){//dE/dx pi-line variation and single pT variation and 2D triangular chi2 and psi pair
    cuts.AddCutPCMCalo("82410113","00200009317400008250400000","2444451041013200000","0h63103100000010"); //0, -10, 0.4, 3
    cuts.AddCutPCMCalo("82410113","00200049327000008250400000","2444451041013200000","0h63103100000010"); //single pT 0.075 GeV/c
    cuts.AddCutPCMCalo("82410113","00200019327000008250400000","2444451041013200000","0h63103100000010"); //single pT 0.1   GeV/c
    cuts.AddCutPCMCalo("82410113","00200009327000008850400000","2444451041013200000","0h63103100000010"); //20 & 0.1
  } else if(trainConfig == 342){//2D triangular chi2 and psi pair
    cuts.AddCutPCMCalo("82410113","00200009327000008260400000","2444451041013200000","0h63103100000010"); //30 & 0.05
    cuts.AddCutPCMCalo("82410113","00200009327000008860400000","2444451041013200000","0h63103100000010"); //20 & 0.05
    cuts.AddCutPCMCalo("82410113","00200009327000008280400000","2444451041013200000","0h63103100000010"); //30 & 0.2
    cuts.AddCutPCMCalo("82410113","00200009327000008880400000","2444451041013200000","0h63103100000010"); //20 & 0.2
  } else if(trainConfig == 343){//min TPC clusters and //elliptic cut Qt and alpha
    cuts.AddCutPCMCalo("82410113","00200000327000008250400000","2444451041013200000","0h63103100000010"); //0
    cuts.AddCutPCMCalo("82410113","00200008327000008250400000","2444451041013200000","0h63103100000010"); //0.35
    cuts.AddCutPCMCalo("82410113","00200009327000009250400000","2444451041013200000","0h63103100000010"); // qT 0.03 no quadratic
    cuts.AddCutPCMCalo("82410113","00200009327000003250400000","2444451041013200000","0h63103100000010"); // qT 0.05 y  quadratic
  } else if(trainConfig == 344){// and min phi max phi and NL variations
    cuts.AddCutPCMCalo("82410113","00209909327000008250400000","2444451041013200000","0h63103100000010"); //4.54 > phi > 5.59
    cuts.AddCutPCMCalo("82410113","00200009f9730000dge0400000","2444452041013200000","0h63103100000010"); // PHOS calo NL
    cuts.AddCutPCMCalo("82410113","00200009f9730000dge0400000","2444401041013200000","0h63103100000010"); // PHOS people NL

    //Cluster Cuts varation
  } else if(trainConfig == 345){ // first set of variations CLUSTER
    cuts.AddCutPCMCalo("82410113","00200009f9730000dge0400000","2444451041073200000","0h63103100000010"); // min energy 0.2 GeV
    cuts.AddCutPCMCalo("82410113","00200009f9730000dge0400000","2444451041083200000","0h63103100000010"); // min energy 0.4 GeV
    cuts.AddCutPCMCalo("82410113","00200009f9730000dge0400000","2444451041023200000","0h63103100000010"); // min energy 0.5 GeV
  } else if(trainConfig == 346){ // second set of variations CLUSTER
    cuts.AddCutPCMCalo("82410113","00200009f9730000dge0400000","2444451041013270000","0h63103100000010"); // min/max M02  0.1<M<1
    cuts.AddCutPCMCalo("82410113","00200009f9730000dge0400000","2444451041013280000","0h63103100000010"); // min/max M02  0.1<M<1
    cuts.AddCutPCMCalo("82410113","00200009f9730000dge0400000","2444451041012200000","0h63103100000010"); // min number 2 cells
    cuts.AddCutPCMCalo("82410113","00200009f9730000dge0400000","2444451041014200000","0h63103100000010"); // min number 4 cells
  } else if(trainConfig == 347){ // MESON
    cuts.AddCutPCMCalo("82410113","00200009f9730000dge0400000","2444451041013200000","0163403100000010"); // rapidity variation  y<0.5
    cuts.AddCutPCMCalo("82410113","00200009f9730000dge0400000","2444451041013200000","0163803100000010"); // rapidity variation  y<0.25
    cuts.AddCutPCMCalo("82410113","00200009f9730000dge0400000","2444451041013200000","0163106100000010"); // alpha meson variation 1 0<alpha<0.8
    cuts.AddCutPCMCalo("82410113","00200009f9730000dge0400000","2444451041013200000","0163105100000010"); // alpha meson variation 2 0<alpha<0.75
  } else if(trainConfig == 348){ // fourth set of variations
    cuts.AddCutPCMCalo("82410113","00200009f9730000dge0400000","2444451044013200000","0h63103100000010"); // tm variation
    cuts.AddCutPCMCalo("82410113","00200009f9730000dge0400000","2444451045013200000","0h63103100000010"); // tm variation
    cuts.AddCutPCMCalo("82410113","00200009f9730000dge0400000","2444451046013200000","0h63103100000010"); // tm variation
    cuts.AddCutPCMCalo("82410113","00200009f9730000dge0400000","2444451041013200000","0h63103100000000"); // min opening angle 0    -> open
    cuts.AddCutPCMCalo("82410113","00200009f9730000dge0400000","2444451041013200000","0h63103100000030"); // min opening angle 0.01 -> 2 cell diag
    //40 - 60
  } else if(trainConfig == 350){//dEdx e-line variation and dE/dx pi-line variation
    cuts.AddCutPCMCalo("84610113","00200009127000008250400000","2444451041013200000","0h63103100000010"); //-5 < sigma < 5
    cuts.AddCutPCMCalo("84610113","00200009227000008250400000","2444451041013200000","0h63103100000010"); //-3 < sigma < 5
    cuts.AddCutPCMCalo("84610113","00200009327400008250400000","2444451041013200000","0h63103100000010"); //1, -10, 0.4, 3
    cuts.AddCutPCMCalo("84610113","00200009367400008250400000","2444451041013200000","0h63103100000010"); //2, 0.5, 0.4, 3
  } else if(trainConfig == 351){//dE/dx pi-line variation and single pT variation and 2D triangular chi2 and psi pair
    cuts.AddCutPCMCalo("84610113","00200009317400008250400000","2444451041013200000","0h63103100000010"); //0, -10, 0.4, 3
    cuts.AddCutPCMCalo("84610113","00200049327000008250400000","2444451041013200000","0h63103100000010"); //single pT 0.075 GeV/c
    cuts.AddCutPCMCalo("84610113","00200019327000008250400000","2444451041013200000","0h63103100000010"); //single pT 0.1   GeV/c
    cuts.AddCutPCMCalo("84610113","00200009327000008850400000","2444451041013200000","0h63103100000010"); //20 & 0.1
  } else if(trainConfig == 352){//2D triangular chi2 and psi pair
    cuts.AddCutPCMCalo("84610113","00200009327000008260400000","2444451041013200000","0h63103100000010"); //30 & 0.05
    cuts.AddCutPCMCalo("84610113","00200009327000008860400000","2444451041013200000","0h63103100000010"); //20 & 0.05
    cuts.AddCutPCMCalo("84610113","00200009327000008280400000","2444451041013200000","0h63103100000010"); //30 & 0.2
    cuts.AddCutPCMCalo("84610113","00200009327000008880400000","2444451041013200000","0h63103100000010"); //20 & 0.2
  } else if(trainConfig == 353){//min TPC clusters and //elliptic cut Qt and alpha
    cuts.AddCutPCMCalo("84610113","00200000327000008250400000","2444451041013200000","0h63103100000010"); //0
    cuts.AddCutPCMCalo("84610113","00200008327000008250400000","2444451041013200000","0h63103100000010"); //0.35
    cuts.AddCutPCMCalo("84610113","00200009327000009250400000","2444451041013200000","0h63103100000010"); // qT 0.03 no quadratic
    cuts.AddCutPCMCalo("84610113","00200009327000003250400000","2444451041013200000","0h63103100000010"); // qT 0.05 y  quadratic
  } else if(trainConfig == 354){// and min phi max phi and NL variations
    cuts.AddCutPCMCalo("84610113","00209909327000008250400000","2444451041013200000","0h63103100000010"); //4.54 > phi > 5.59
    cuts.AddCutPCMCalo("84610113","00200009f9730000dge0400000","2444452041013200000","0h63103100000010"); // PHOS calo NL
    cuts.AddCutPCMCalo("84610113","00200009f9730000dge0400000","2444401041013200000","0h63103100000010"); // PHOS people NL

    //Cluster Cuts varation
  } else if(trainConfig == 355){ // first set of variations CLUSTER
    cuts.AddCutPCMCalo("84610113","00200009f9730000dge0400000","2444451041073200000","0h63103100000010"); // min energy 0.2 GeV
    cuts.AddCutPCMCalo("84610113","00200009f9730000dge0400000","2444451041083200000","0h63103100000010"); // min energy 0.4 GeV
    cuts.AddCutPCMCalo("84610113","00200009f9730000dge0400000","2444451041023200000","0h63103100000010"); // min energy 0.5 GeV
  } else if(trainConfig == 356){ // second set of variations CLUSTER
    cuts.AddCutPCMCalo("84610113","00200009f9730000dge0400000","2444451041013270000","0h63103100000010"); // min/max M02  0.1<M<1
    cuts.AddCutPCMCalo("84610113","00200009f9730000dge0400000","2444451041013280000","0h63103100000010"); // min/max M02  0.1<M<1
    cuts.AddCutPCMCalo("84610113","00200009f9730000dge0400000","2444451041012200000","0h63103100000010"); // min number 2 cells
    cuts.AddCutPCMCalo("84610113","00200009f9730000dge0400000","2444451041014200000","0h63103100000010"); // min number 4 cells
  } else if(trainConfig == 357){ // MESON
    cuts.AddCutPCMCalo("84610113","00200009f9730000dge0400000","2444451041013200000","0h63403100000010"); // rapidity variation  y<0.5
    cuts.AddCutPCMCalo("84610113","00200009f9730000dge0400000","2444451041013200000","0h63803100000010"); // rapidity variation  y<0.25
    cuts.AddCutPCMCalo("84610113","00200009f9730000dge0400000","2444451041013200000","0h63106100000010"); // alpha meson variation 1 0<alpha<0.8
    cuts.AddCutPCMCalo("84610113","00200009f9730000dge0400000","2444451041013200000","0h63105100000010"); // alpha meson variation 2 0<alpha<0.75
  } else if(trainConfig == 358){ // fourth set of variations
    cuts.AddCutPCMCalo("84610113","00200009f9730000dge0400000","2444451044013200000","0h63103100000010"); // tm variation
    cuts.AddCutPCMCalo("84610113","00200009f9730000dge0400000","2444451045013200000","0h63103100000010"); // tm variation
    cuts.AddCutPCMCalo("84610113","00200009f9730000dge0400000","2444451046013200000","0h63103100000010"); // tm variation
    cuts.AddCutPCMCalo("84610113","00200009f9730000dge0400000","2444451041013200000","0h63103100000000"); // min opening angle 0    -> open
    cuts.AddCutPCMCalo("84610113","00200009f9730000dge0400000","2444451041013200000","0h63103100000030"); // min opening angle 0.01 -> 2 cell diag
    //60 - 100
  } else if(trainConfig == 360){//dEdx e-line variation and dE/dx pi-line variation
    cuts.AddCutPCMCalo("86010113","00200009127000008250400000","2444451041013200000","0h63103100000010"); //-5 < sigma < 5
    cuts.AddCutPCMCalo("86010113","00200009227000008250400000","2444451041013200000","0h63103100000010"); //-3 < sigma < 5
    cuts.AddCutPCMCalo("86010113","00200009327400008250400000","2444451041013200000","0h63103100000010"); //1, -10, 0.4, 3
    cuts.AddCutPCMCalo("86010113","00200009367400008250400000","2444451041013200000","0h63103100000010"); //2, 0.5, 0.4, 3
  } else if(trainConfig == 361){//dE/dx pi-line variation and single pT variation and 2D triangular chi2 and psi pair
    cuts.AddCutPCMCalo("86010113","00200009317400008250400000","2444451041013200000","0h63103100000010"); //0, -10, 0.4, 3
    cuts.AddCutPCMCalo("86010113","00200049327000008250400000","2444451041013200000","0h63103100000010"); //single pT 0.075 GeV/c
    cuts.AddCutPCMCalo("86010113","00200019327000008250400000","2444451041013200000","0h63103100000010"); //single pT 0.1   GeV/c
    cuts.AddCutPCMCalo("86010113","00200009327000008850400000","2444451041013200000","0h63103100000010"); //20 & 0.1
  } else if(trainConfig == 362){//2D triangular chi2 and psi pair
    cuts.AddCutPCMCalo("86010113","00200009327000008260400000","2444451041013200000","0h63103100000010"); //30 & 0.05
    cuts.AddCutPCMCalo("86010113","00200009327000008860400000","2444451041013200000","0h63103100000010"); //20 & 0.05
    cuts.AddCutPCMCalo("86010113","00200009327000008280400000","2444451041013200000","0h63103100000010"); //30 & 0.2
    cuts.AddCutPCMCalo("86010113","00200009327000008880400000","2444451041013200000","0h63103100000010"); //20 & 0.2
  } else if(trainConfig == 363){//min TPC clusters and //elliptic cut Qt and alpha
    cuts.AddCutPCMCalo("86010113","00200000327000008250400000","2444451041013200000","0h63103100000010"); //0
    cuts.AddCutPCMCalo("86010113","00200008327000008250400000","2444451041013200000","0h63103100000010"); //0.35
    cuts.AddCutPCMCalo("86010113","00200009327000009250400000","2444451041013200000","0h63103100000010"); // qT 0.03 no quadratic
    cuts.AddCutPCMCalo("86010113","00200009327000003250400000","2444451041013200000","0h63103100000010"); // qT 0.05 y  quadratic
  } else if(trainConfig == 364){// and min phi max phi and NL variations
    cuts.AddCutPCMCalo("86010113","00209909327000008250400000","2444451041013200000","0h63103100000010"); //4.54 > phi > 5.59
    cuts.AddCutPCMCalo("86010113","00200009f9730000dge0400000","2444452041013200000","0h63103100000010"); // PHOS calo NL
    cuts.AddCutPCMCalo("86010113","00200009f9730000dge0400000","2444401041013200000","0h63103100000010"); // PHOS people NL

    //Cluster Cuts varation
  } else if(trainConfig == 365){ // first set of variations CLUSTER
    cuts.AddCutPCMCalo("86010113","00200009f9730000dge0400000","2444451041073200000","0h63103100000010"); // min energy 0.2 GeV
    cuts.AddCutPCMCalo("86010113","00200009f9730000dge0400000","2444451041083200000","0h63103100000010"); // min energy 0.4 GeV
    cuts.AddCutPCMCalo("86010113","00200009f9730000dge0400000","2444451041023200000","0h63103100000010"); // min energy 0.5 GeV
  } else if(trainConfig == 366){ // second set of variations CLUSTER
    cuts.AddCutPCMCalo("86010113","00200009f9730000dge0400000","2444451041013270000","0h63103100000010"); // min/max M02  0.1<M<1
    cuts.AddCutPCMCalo("86010113","00200009f9730000dge0400000","2444451041013280000","0h63103100000010"); // min/max M02  0.1<M<1
    cuts.AddCutPCMCalo("86010113","00200009f9730000dge0400000","2444451041012200000","0h63103100000010"); // min number 2 cells
    cuts.AddCutPCMCalo("86010113","00200009f9730000dge0400000","2444451041014200000","0h63103100000010"); // min number 4 cells
  } else if(trainConfig == 367){ // MESON
    cuts.AddCutPCMCalo("86010113","00200009f9730000dge0400000","2444451041013200000","0h63403100000010"); // rapidity variation  y<0.5
    cuts.AddCutPCMCalo("86010113","00200009f9730000dge0400000","2444451041013200000","0h63803100000010"); // rapidity variation  y<0.25
    cuts.AddCutPCMCalo("86010113","00200009f9730000dge0400000","2444451041013200000","0h63106100000010"); // alpha meson variation 1 0<alpha<0.8
    cuts.AddCutPCMCalo("86010113","00200009f9730000dge0400000","2444451041013200000","0h63105100000010"); // alpha meson variation 2 0<alpha<0.75
  } else if(trainConfig == 368){ // fourth set of variations
    cuts.AddCutPCMCalo("86010113","00200009f9730000dge0400000","2444451044013200000","0h63103100000010"); // tm variation
    cuts.AddCutPCMCalo("86010113","00200009f9730000dge0400000","2444451045013200000","0h63103100000010"); // tm variation
    cuts.AddCutPCMCalo("86010113","00200009f9730000dge0400000","2444451046013200000","0h63103100000010"); // tm variation
    cuts.AddCutPCMCalo("86010113","00200009f9730000dge0400000","2444451041013200000","0h63103100000000"); // min opening angle 0    -> open
    cuts.AddCutPCMCalo("86010113","00200009f9730000dge0400000","2444451041013200000","0h63103100000030"); // min opening angle 0.01 -> 2 cell diag

  } else if(trainConfig == 371){ // PCM-PHOS other cent estimators
    cuts.AddCutPCMCalo("90010113","00200009f9730000dge0400000","2444451044013200000","0h63103100000010"); // 0-100% with PCM NL
    cuts.AddCutPCMCalo("90210113","00200009f9730000dge0400000","2444451044013200000","0h63103100000010"); // 0-20% with PCM NL
    cuts.AddCutPCMCalo("92410113","00200009f9730000dge0400000","2444451044013200000","0h63103100000010"); // 20-40% with PCM NL
    cuts.AddCutPCMCalo("94610113","00200009f9730000dge0400000","2444451044013200000","0h63103100000010"); // 40-60% with PCM NL
    cuts.AddCutPCMCalo("96010113","00200009f9730000dge0400000","2444451044013200000","0h63103100000010"); // 60-100% with PCM NL
  } else if(trainConfig == 372){ // PCM-PHOS other cent estimators
    cuts.AddCutPCMCalo("e0010113","00200009f9730000dge0400000","2444451044013200000","0h63103100000010"); // 0-100% with PCM NL
    cuts.AddCutPCMCalo("e0210113","00200009f9730000dge0400000","2444451044013200000","0h63103100000010"); // 0-20% with PCM NL
    cuts.AddCutPCMCalo("e2410113","00200009f9730000dge0400000","2444451044013200000","0h63103100000010"); // 20-40% with PCM NL
    cuts.AddCutPCMCalo("e4610113","00200009f9730000dge0400000","2444451044013200000","0h63103100000010"); // 40-60% with PCM NL
    cuts.AddCutPCMCalo("e6010113","00200009f9730000dge0400000","2444451044013200000","0h63103100000010"); // 60-100% with PCM NL

  // ===============================================================================================
  // Run 2 data EMC clusters pPb 8TeV
  // ===============================================================================================
  } else if (trainConfig == 400){  // EMCal clusters standard cuts, triggers, no nonlin, open timing
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","1111100017032230000","0h63103100000010"); // INT7
  } else if (trainConfig == 401){  // EMCal clusters standard cuts, triggers, no nonlin, open timing
    cuts.AddCutPCMCalo("80085113","00200009f9730000dge0400000","1111100017032230000","0h63103100000010"); // EG2
    cuts.AddCutPCMCalo("80083113","00200009f9730000dge0400000","1111100017032230000","0h63103100000010"); // EG1
  } else if (trainConfig == 402){  // EMCal clusters standard cuts, triggers, no nonlin, open timing
    cuts.AddCutPCMCalo("80095113","00200009f9730000dge0400000","1111100017032230000","0h63103100000010"); // EJ2
    cuts.AddCutPCMCalo("80093113","00200009f9730000dge0400000","1111100017032230000","0h63103100000010"); // EJ1

  } else if (trainConfig == 404){  // EMCal clusters standard cuts, triggers, no nonlin, +-50ns timing (5)
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","1111100057032230000","0h63103100000010"); // INT7
  } else if (trainConfig == 405){  // EMCal clusters standard cuts, triggers, no nonlin, +-50ns timing (5)
    cuts.AddCutPCMCalo("80085113","00200009f9730000dge0400000","1111100057032230000","0h63103100000010"); // EG2
    cuts.AddCutPCMCalo("80083113","00200009f9730000dge0400000","1111100057032230000","0h63103100000010"); // EG1
  } else if (trainConfig == 406){  // EMCal clusters standard cuts, triggers, no nonlin, +-50ns timing (5)
    cuts.AddCutPCMCalo("80095113","00200009f9730000dge0400000","1111100057032230000","0h63103100000010"); // EJ2
    cuts.AddCutPCMCalo("80093113","00200009f9730000dge0400000","1111100057032230000","0h63103100000010"); // EJ1

  } else if (trainConfig == 407){  // EMCal clusters standard cuts, triggers, Nico TB NL, +-50ns timing (5)
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","1111165057032230000","0h63103100000010"); // INT7
  } else if (trainConfig == 408){  // EMCal clusters standard cuts, triggers, Nico TB NL, +-50ns timing (5)
    cuts.AddCutPCMCalo("80085113","00200009f9730000dge0400000","1111165057032230000","0h63103100000010"); // EG2
    cuts.AddCutPCMCalo("80083113","00200009f9730000dge0400000","1111165057032230000","0h63103100000010"); // EG1
  } else if (trainConfig == 409){  // EMCal clusters standard cuts, triggers, Nico TB NL, +-50ns timing (5)
    cuts.AddCutPCMCalo("80095113","00200009f9730000dge0400000","1111165057032230000","0h63103100000010"); // EJ2
    cuts.AddCutPCMCalo("80093113","00200009f9730000dge0400000","1111165057032230000","0h63103100000010"); // EJ1

  } else if (trainConfig == 410){  // EMCal clusters standard cuts, triggers, Nico TB NL, +-50ns timing (5)
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","1111147057032230000","0h63103100000010"); // INT7
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","1111148057032230000","0h63103100000010"); // INT7
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","1111157057032230000","0h63103100000010"); // INT7
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","1111158057032230000","0h63103100000010"); // INT7
  } else if (trainConfig == 411){  // EMCal clusters standard cuts, triggers, Nico TB NL, +-50ns timing (5)
    cuts.AddCutPCMCalo("80085113","00200009f9730000dge0400000","1111147057032230000","0h63103100000010"); // EG2
    cuts.AddCutPCMCalo("80085113","00200009f9730000dge0400000","1111148057032230000","0h63103100000010"); // EG2
    cuts.AddCutPCMCalo("80085113","00200009f9730000dge0400000","1111157057032230000","0h63103100000010"); // EG2
    cuts.AddCutPCMCalo("80085113","00200009f9730000dge0400000","1111158057032230000","0h63103100000010"); // EG2
  } else if (trainConfig == 412){  // EMCal clusters standard cuts, triggers, Nico TB NL, +-50ns timing (5)
    cuts.AddCutPCMCalo("80083113","00200009f9730000dge0400000","1111147057032230000","0h63103100000010"); // EG1
    cuts.AddCutPCMCalo("80083113","00200009f9730000dge0400000","1111148057032230000","0h63103100000010"); // EG1
    cuts.AddCutPCMCalo("80083113","00200009f9730000dge0400000","1111157057032230000","0h63103100000010"); // EG1
    cuts.AddCutPCMCalo("80083113","00200009f9730000dge0400000","1111158057032230000","0h63103100000010"); // EG1

  } else if (trainConfig == 413){  // EMCal clusters standard cuts, triggers, Nico TB NL, +-50ns timing (5), TM E/p 1.75 (f)
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","111114705f032230000","0h63103100000010"); // INT7
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","111114805f032230000","0h63103100000010"); // INT7
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","111115705f032230000","0h63103100000010"); // INT7
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","111115805f032230000","0h63103100000010"); // INT7
  } else if (trainConfig == 414){  // EMCal clusters standard cuts, triggers, Nico TB NL, +-50ns timing (5), TM E/p 1.75 (f)
    cuts.AddCutPCMCalo("80085113","00200009f9730000dge0400000","111114705f032230000","0h63103100000010"); // EG2
    cuts.AddCutPCMCalo("80085113","00200009f9730000dge0400000","111114805f032230000","0h63103100000010"); // EG2
    cuts.AddCutPCMCalo("80085113","00200009f9730000dge0400000","111115705f032230000","0h63103100000010"); // EG2
    cuts.AddCutPCMCalo("80085113","00200009f9730000dge0400000","111115805f032230000","0h63103100000010"); // EG2
  } else if (trainConfig == 415){  // EMCal clusters standard cuts, triggers, Nico TB NL, +-50ns timing (5), TM E/p 1.75 (f)
    cuts.AddCutPCMCalo("80083113","00200009f9730000dge0400000","111114705f032230000","0h63103100000010"); // EG1
    cuts.AddCutPCMCalo("80083113","00200009f9730000dge0400000","111114805f032230000","0h63103100000010"); // EG1
    cuts.AddCutPCMCalo("80083113","00200009f9730000dge0400000","111115705f032230000","0h63103100000010"); // EG1
    cuts.AddCutPCMCalo("80083113","00200009f9730000dge0400000","111115805f032230000","0h63103100000010"); // EG1

  // Systematics pPb 8 TeV Min Bias
  } else if (trainConfig == 420){ //EMCAL minEnergy variation
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","1111151057032230000","0h63103100000010"); // 0.7 GeV/c default
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","1111151057042230000","0h63103100000010"); // 0.8 GeV/c
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","1111151057052230000","0h63103100000010"); // 0.9 GeV/c
  } else if (trainConfig == 421){ //EMCAL minNCells variation
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","1111151057031230000","0h63103100000010"); // n cells >= 1
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","1111151057033230000","0h63103100000010"); // n cells >= 3
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","1111151057032200000","0h63103100000010"); // no max M02 cut
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","1111151057032250000","0h63103100000010"); // M02 < 0.3
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","1111151057032260000","0h63103100000010"); // M02 < 0.27
  } else if (trainConfig == 422){ // EMCAL track matching variations
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","1111151053032230000","0h63103100000010"); //
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","1111151056032230000","0h63103100000010"); //
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","1111151058032230000","0h63103100000010"); //
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","1111151059032230000","0h63103100000010"); //
  } else if (trainConfig == 423){ // PCM variations
    cuts.AddCutPCMCalo("80010113","00200009227000008250400000","1111151057032230000","0h63103100000010"); // dEdx e -3, 5
    cuts.AddCutPCMCalo("80010113","00200009127000008250400000","1111151057032230000","0h63103100000010"); // dEdx e -5, 5
    cuts.AddCutPCMCalo("80010113","00200009357000008250400000","1111151057032230000","0h63103100000010"); // dEdx pi 2
    cuts.AddCutPCMCalo("80010113","00200009317000008250400000","1111151057032230000","0h63103100000010"); // dEdx pi 0
    cuts.AddCutPCMCalo("80010113","00200009387300008250400000","1111151057032230000","0h63103100000010"); // dEdx pi 2 high 1 (> 3.5 GeV)
  } else if (trainConfig == 424){ // PCM variations
    cuts.AddCutPCMCalo("80010113","00200009327000009250400000","1111151057032230000","0h63103100000010"); // qt 2D 0.03
    cuts.AddCutPCMCalo("80010113","00200009327000003250400000","1111151057032230000","0h63103100000010"); // qt 1D 0.05
    cuts.AddCutPCMCalo("80010113","00200009327000002250400000","1111151057032230000","0h63103100000010"); // qt 1D 0.07
    cuts.AddCutPCMCalo("80010113","00200049327000008250400000","1111151057032230000","0h63103100000010"); // single pt > 0.075
    cuts.AddCutPCMCalo("80010113","00200019327000008250400000","1111151057032230000","0h63103100000010"); // single pt > 0.1
  } else if (trainConfig == 425){ // PCM variations
    cuts.AddCutPCMCalo("80010113","00200009327000008850400000","1111151057032230000","0h63103100000010"); // 2D psi pair chi2 var
    cuts.AddCutPCMCalo("80010113","00200009327000008260400000","1111151057032230000","0h63103100000010"); // 2D psi pair chi2 var
    cuts.AddCutPCMCalo("80010113","00200009327000008860400000","1111151057032230000","0h63103100000010"); // 2D psi pair chi2 var
    cuts.AddCutPCMCalo("80010113","00200009327000008280400000","1111151057032230000","0h63103100000010"); // 2D psi pair chi2 var
    cuts.AddCutPCMCalo("80010113","00200009327000008880400000","1111151057032230000","0h63103100000010"); // 2D psi pair chi2 var
  } else if (trainConfig == 426){ // PCM variations
    cuts.AddCutPCMCalo("80010113","00200006327000008250400000","1111151057032230000","0h63103100000010"); // min TPC cl > 0.7
    cuts.AddCutPCMCalo("80010113","00200008327000008250400000","1111151057032230000","0h63103100000010"); // min TPC cl > 0.35
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","1111151057032230000","0163106100000010"); // alpha < 0.9
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","1111151057032230000","0163105100000010"); // alpha < 0.75
  } else if (trainConfig == 427){ // PCM variations
    cuts.AddCutPCMCalo("80010113","00202209327000008250400000","1111151057032230000","0h63103100000010"); // restrict acceptance to EMCAL loose
    cuts.AddCutPCMCalo("80010113","00204409327000008250400000","1111151057032230000","0h63103100000010"); // restrict acceptance to EMCAL tight
  } else if (trainConfig == 428){ // PCM variations pi dEdx
    cuts.AddCutPCMCalo("80010113","00200009317300008250400000","1111151057032230000","0h63103100000010"); // dEdx pi: 0: 0.4-3.5, -10: 3.5 ->
    cuts.AddCutPCMCalo("80010113","00200009327300008250400000","1111151057032230000","0h63103100000010"); // dEdx pi: 1: 0.4-3.5, -10: 3.5 ->
    cuts.AddCutPCMCalo("80010113","00200009325000008250400000","1111151057032230000","0h63103100000010"); // dEdx pi: 1: 0.3 ->
    cuts.AddCutPCMCalo("80010113","00200009320000008250400000","1111151057032230000","0h63103100000010"); // dEdx pi: 1: 0.5 ->
  } else if (trainConfig == 429){ // PCM variations pi dEdx
    cuts.AddCutPCMCalo("80010113","00200009327600008250400000","1111151057032230000","0h63103100000010"); // dEdx pi: 1: 0.4-2, -10: 2. ->
    cuts.AddCutPCMCalo("80010113","00200009327400008250400000","1111151057032230000","0h63103100000010"); // dEdx pi: 1: 0.4-3, -10: 3. ->
    cuts.AddCutPCMCalo("80010113","00200009315600008250400000","1111151057032230000","0h63103100000010"); // dEdx pi: 0: 0.3-2, -10: 2. ->
    cuts.AddCutPCMCalo("80010113","00200009367400008250400000","1111151057032230000","0h63103100000010"); // dEdx pi: 2: 0.4-3, 0.5: 3. ->
    cuts.AddCutPCMCalo("80010113","00200009347400008250400000","1111151057032230000","0h63103100000010"); // dEdx pi: 3: 0.4-3, 1: 3. ->
  } else if (trainConfig == 430){ // PCM variations to close V0s
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","1111151057032230000","0h63103100000010");
    cuts.AddCutPCMCalo("80010113","00200009327000008250401000","1111151057032230000","0h63103100000010");
    cuts.AddCutPCMCalo("80010113","00200009327000008250402000","1111151057032230000","0h63103100000010");
    cuts.AddCutPCMCalo("80010113","00200009327000008250403000","1111151057032230000","0h63103100000010");

  } else if ( trainConfig == 480){ // EMCAL standard cut but CALO+CALOFAST - NO TM
    cuts.AddCutPCMCalo("800a0113","00200009f9730000dge0400000","1111100050032230000","0h63103100000010"); // std INT7
    cuts.AddCutPCMCalo("800a1113","00200009f9730000dge0400000","1111100050032230000","0h63103100000010"); // std EMC7
    cuts.AddCutPCMCalo("800a2113","00200009f9730000dge0400000","1111100050032230000","0h63103100000010"); // std EG2
    cuts.AddCutPCMCalo("800a3113","00200009f9730000dge0400000","1111100050032230000","0h63103100000010"); // std EG1
  } else if ( trainConfig == 481){ // EMCAL standard cut but standard readout - NO TM
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","1111100050032230000","0h63103100000010"); // std INT7
    cuts.AddCutPCMCalo("80052113","00200009f9730000dge0400000","1111100050032230000","0h63103100000010"); // std EMC7
    cuts.AddCutPCMCalo("80085113","00200009f9730000dge0400000","1111100050032230000","0h63103100000010"); // std EG2
    cuts.AddCutPCMCalo("80083113","00200009f9730000dge0400000","1111100050032230000","0h63103100000010"); // std EG1

  // ===============================================================================================
  // Run 2 data PHOS clusters 5TeV
  // ===============================================================================================
  // INT7 triggers
  } else if (trainConfig == 500) {  // PHOS  INT7
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","2446600051012200000","0h63103100000010"); // PHOS group standard +-50ns
  } else if (trainConfig == 501) {  // PHOS  INT7
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","2446600041012200000","0h63103100000010"); // PHOS group standard
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","2446600011012200000","0h63103100000010"); // PHOS group standard, 1000 \mus
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","2446600061012200000","0h63103100000010"); // PHOS group standard, -30, 50ns
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","24466000a1012200000","0h63103100000010"); // PHOS group standard, -12.5, 13ns

  } else if (trainConfig == 502) {  // PHOS  INT7 with cents
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","2446600051012200000","0h63103100000010"); // non lin 0-100%
    cuts.AddCutPCMCalo("80110113","00200009f9730000dge0400000","2446600051012200000","0h63103100000010"); // non lin 0-10%
    cuts.AddCutPCMCalo("81210113","00200009f9730000dge0400000","2446600051012200000","0h63103100000010"); // non lin 10-20%
    cuts.AddCutPCMCalo("82410113","00200009f9730000dge0400000","2446600051012200000","0h63103100000010"); // non lin 20-40%
    cuts.AddCutPCMCalo("84610113","00200009f9730000dge0400000","2446600051012200000","0h63103100000010"); // non lin 40-60%
    cuts.AddCutPCMCalo("86810113","00200009f9730000dge0400000","2446600051012200000","0h63103100000010"); // non lin 60-80%
    cuts.AddCutPCMCalo("88010113","00200009f9730000dge0400000","2446600051012200000","0h63103100000010"); // non lin 80-100%
  } else if (trainConfig == 503) {  // PHOS  INT7 with cents
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","2446600051012200000","0h63103100000010"); // non lin 0-100%
    cuts.AddCutPCMCalo("80210113","00200009f9730000dge0400000","2446600051012200000","0h63103100000010"); // non lin 0-10%
    cuts.AddCutPCMCalo("86010113","00200009f9730000dge0400000","2446600051012200000","0h63103100000010"); // non lin 10-20%
    cuts.AddCutPCMCalo("a0110113","00200009f9730000dge0400000","2446600051012200000","0h63103100000010"); // non lin 0-5%
    cuts.AddCutPCMCalo("a1210113","00200009f9730000dge0400000","2446600051012200000","0h63103100000010"); // non lin 5-10%
  } else if (trainConfig == 504) {  // PHOS  INT7
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","2446600051012200000","0h63103100000010"); // non lin none
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","2446641051012200000","0h63103100000010"); // non lin 0-100%
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","2446651051012200000","0h63103100000010"); // non lin 0-100%
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","2446601051012200000","0h63103100000010"); // non lin PHOS 0-100%
  } else if (trainConfig == 505) {  // PHOS  INT7 with cents
    cuts.AddCutPCMCalo("80110113","00200009f9730000dge0400000","2446601051012200000","0h63103100000010"); // non lin 0-10%
    cuts.AddCutPCMCalo("81210113","00200009f9730000dge0400000","2446601051012200000","0h63103100000010"); // non lin 10-20%
    cuts.AddCutPCMCalo("82410113","00200009f9730000dge0400000","2446601051012200000","0h63103100000010"); // non lin 20-40%
    cuts.AddCutPCMCalo("84610113","00200009f9730000dge0400000","2446601051012200000","0h63103100000010"); // non lin 40-60%
    cuts.AddCutPCMCalo("86810113","00200009f9730000dge0400000","2446601051012200000","0h63103100000010"); // non lin 60-80%
    cuts.AddCutPCMCalo("88010113","00200009f9730000dge0400000","2446601051012200000","0h63103100000010"); // non lin 80-100%
  } else if (trainConfig == 506) {  // PHOS  INT7 with cents
    cuts.AddCutPCMCalo("80210113","00200009f9730000dge0400000","2446601051012200000","0h63103100000010"); // non lin 0-10%
    cuts.AddCutPCMCalo("86010113","00200009f9730000dge0400000","2446601051012200000","0h63103100000010"); // non lin 10-20%
    cuts.AddCutPCMCalo("a0110113","00200009f9730000dge0400000","2446601051012200000","0h63103100000010"); // non lin 0-5%
    cuts.AddCutPCMCalo("a1210113","00200009f9730000dge0400000","2446601051012200000","0h63103100000010"); // non lin 5-10%
  } else if (trainConfig == 507) {  // PHOS  INT7 with narrow cents
    cuts.AddCutPCMCalo("c0110113","00200009f9730000dge0400000","2446641051012200000","0h63103100000010"); // non lin 0-1%
    cuts.AddCutPCMCalo("c0210113","00200009f9730000dge0400000","2446641051012200000","0h63103100000010"); // non lin 0-2%
  } else if (trainConfig == 508) {  // AOD validation
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","2446641051012200000","0h63103100000010"); // non lin 0-100%
  } else if (trainConfig == 509) {  // AOD validation PHOS NL
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","2446601051012200000","0h63103100000010"); // non lin 0-100%
  } else if (trainConfig == 510) {  // AOD validation no NL
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","2446600051012200000","0h63103100000010"); // non lin 0-100%
  } else if (trainConfig == 511) {  // JJ AOD validation no NL
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","2446600051012200000","0h63103100000010"); // no non lin 0-100%
    cuts.AddCutPCMCalo("80010123","00200009f9730000dge0400000","2446600051012200000","0h63103100000010"); // no non lin 0-100%

  } else if (trainConfig == 512) {  // PHOS  INT7
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","2446600050012200000","0h63103100000010"); // TM studies
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","2446600051012200000","0h63103100000010"); // TM studies
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","2446600054012200000","0h63103100000010"); // TM studies
  } else if (trainConfig == 513) {  // PHOS  INT7
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","2446600051012200000","0h63103100000010"); // non-lin variations
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","2446642051012200000","0h63103100000010"); // non-lin variations
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","2446652051012200000","0h63103100000010"); // non-lin variations
  } else if (trainConfig == 514) {  // PHOS  INT7
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","2446641054012200000","0h63103100000010"); //
  } else if (trainConfig == 515) {  // PHOS  INT7 with cents
    cuts.AddCutPCMCalo("80110113","00200009f9730000dge0400000","2446641054012200000","0h63103100000010"); // non lin 0-10%
    cuts.AddCutPCMCalo("81210113","00200009f9730000dge0400000","2446641054012200000","0h63103100000010"); // non lin 10-20%
    cuts.AddCutPCMCalo("82410113","00200009f9730000dge0400000","2446641054012200000","0h63103100000010"); // non lin 20-40%
    cuts.AddCutPCMCalo("84610113","00200009f9730000dge0400000","2446641054012200000","0h63103100000010"); // non lin 40-60%
    cuts.AddCutPCMCalo("86810113","00200009f9730000dge0400000","2446641054012200000","0h63103100000010"); // non lin 60-80%
    cuts.AddCutPCMCalo("88010113","00200009f9730000dge0400000","2446641054012200000","0h63103100000010"); // non lin 80-100%
  } else if (trainConfig == 516) {  // PHOS  INT7
    cuts.AddCutPCMCalo("80210113","00200009f9730000dge0400000","2446641054012200000","0h63103100000010"); //
    cuts.AddCutPCMCalo("86010113","00200009f9730000dge0400000","2446641054012200000","0h63103100000010"); //
    cuts.AddCutPCMCalo("a0110113","00200009f9730000dge0400000","2446641054012200000","0h63103100000010"); //
    cuts.AddCutPCMCalo("a1210113","00200009f9730000dge0400000","2446641054012200000","0h63103100000010"); //
  // PCM-PHOS run2 HBT studies
  } else if (trainConfig == 520) {  // PHOS  INT7
    cuts.AddCutPCMCalo("80010113","00600009a27000006250800000","244664105a012200000","0h63103100000010"); //
  } else if (trainConfig == 521) {  // PHOS  INT7 with cents
    cuts.AddCutPCMCalo("80110113","00600009a27000006250800000","244664105a012200000","0h63103100000010"); // non lin 0-10%
    cuts.AddCutPCMCalo("81210113","00600009a27000006250800000","244664105a012200000","0h63103100000010"); // non lin 10-20%
    cuts.AddCutPCMCalo("82410113","00600009a27000006250800000","244664105a012200000","0h63103100000010"); // non lin 20-40%
    cuts.AddCutPCMCalo("84610113","00600009a27000006250800000","244664105a012200000","0h63103100000010"); // non lin 40-60%
    cuts.AddCutPCMCalo("86810113","00600009a27000006250800000","244664105a012200000","0h63103100000010"); // non lin 60-80%
    cuts.AddCutPCMCalo("88010113","00600009a27000006250800000","244664105a012200000","0h63103100000010"); // non lin 80-100%

  //PCM-PHOS run 2 MB systematics
  } else if(trainConfig == 530){//dEdx e-line variation and dE/dx pi-line variation
    cuts.AddCutPCMCalo("80010113","00200009127000008250400000","2446641054012200000","0h63103100000010"); //-5 < sigma < 5
    cuts.AddCutPCMCalo("80010113","00200009227000008250400000","2446641054012200000","0h63103100000010"); //-3 < sigma < 5
    cuts.AddCutPCMCalo("80010113","00200009327400008250400000","2446641054012200000","0h63103100000010"); //1, -10, 0.4, 3
    cuts.AddCutPCMCalo("80010113","00200009367400008250400000","2446641054012200000","0h63103100000010"); //2, 0.5, 0.4, 3
  } else if(trainConfig == 531){//dE/dx pi-line variation and single pT variation and 2D triangular chi2 and psi pair
    cuts.AddCutPCMCalo("80010113","00200009317400008250400000","2446641054012200000","0h63103100000010"); //0, -10, 0.4, 3
    cuts.AddCutPCMCalo("80010113","00200049327000008250400000","2446641054012200000","0h63103100000010"); //single pT 0.075 GeV/c
    cuts.AddCutPCMCalo("80010113","00200019327000008250400000","2446641054012200000","0h63103100000010"); //single pT 0.1   GeV/c
    cuts.AddCutPCMCalo("80010113","00200009327000008850400000","2446641054012200000","0h63103100000010"); //20 & 0.1
  } else if(trainConfig == 532){//2D triangular chi2 and psi pair
    cuts.AddCutPCMCalo("80010113","00200009327000008260400000","2446641054012200000","0h63103100000010"); //30 & 0.05
    cuts.AddCutPCMCalo("80010113","00200009327000008860400000","2446641054012200000","0h63103100000010"); //20 & 0.05
    cuts.AddCutPCMCalo("80010113","00200009327000008280400000","2446641054012200000","0h63103100000010"); //30 & 0.2
    cuts.AddCutPCMCalo("80010113","00200009327000008880400000","2446641054012200000","0h63103100000010"); //20 & 0.2
  } else if(trainConfig == 533){//min TPC clusters and //elliptic cut Qt and alpha
    cuts.AddCutPCMCalo("80010113","00200000327000008250400000","2446641054012200000","0h63103100000010"); //0
    cuts.AddCutPCMCalo("80010113","00200008327000008250400000","2446641054012200000","0h63103100000010"); //0.35
    cuts.AddCutPCMCalo("80010113","00200009327000009250400000","2446641054012200000","0h63103100000010"); // qT 0.03 no quadratic
    cuts.AddCutPCMCalo("80010113","00200009327000003250400000","2446641054012200000","0h63103100000010"); // qT 0.05 y  quadratic
  } else if(trainConfig == 534){// and min phi max phi and NL variations
    cuts.AddCutPCMCalo("80010113","00209909327000008250400000","2446641054012200000","0h63103100000010"); //4.54 > phi > 5.59
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","2446642054012200000","0h63103100000010"); // PHOS calo NL
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","2446652054012200000","0h63103100000010"); // PHOS calo NL
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","2446601054012200000","0h63103100000010"); // PHOS people NL
  } else if(trainConfig == 535){ // first set of variations CLUSTER
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","2446641054072200000","0h63103100000010"); // min energy 0.2 GeV
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","2446641054082200000","0h63103100000010"); // min energy 0.4 GeV
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","2446641054022200000","0h63103100000010"); // min energy 0.5 GeV
  } else if(trainConfig == 536){ // second set of variations CLUSTER
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","2446641054012270000","0h63103100000010"); // min/max M02  0.1<M<1
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","2446641054012280000","0h63103100000010"); // min/max M02  0.1<M<1
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","2446641054013200000","0h63103100000010"); // min number 2 cells
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","2446641054014200000","0h63103100000010"); // min number 4 cells
  } else if(trainConfig == 537){ // MESON
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","2446641054012200000","0h63403100000010"); // rapidity variation  y<0.5
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","2446641054012200000","0h63803100000010"); // rapidity variation  y<0.25
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","2446641054012200000","0h63106100000010"); // alpha meson variation 1 0<alpha<0.8
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","2446641054012200000","0h63105100000010"); // alpha meson variation 2 0<alpha<0.75
  } else if(trainConfig == 538){ // fourth set of variations
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","2446641051012200000","0h63103100000010"); // tm variation
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","2446641055012200000","0h63103100000010"); // tm variation
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","2446641056012200000","0h63103100000010"); // tm variation
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","2446641054012200000","0h63103100000000"); // min opening angle 0    -> open
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","2446641054012200000","0h63103100000030"); // min opening angle 0.01 -> 2 cell diag

  // ===============================================================================================
  // Run 2 data PHOS clusters 8TeV
  // ===============================================================================================
  } else if (trainConfig == 600) {  // PHOS PHI7
    cuts.AddCutPCMCalo("80001113","00200009f9730000dge0400000","2446600051012200000","0h63103100000010");// PHI7 triggers
    cuts.AddCutPCMCalo("80001113","00200009f9730000dge0400000","2446600011012200000","0h63103100000010");// PHI7 triggers
  } else if (trainConfig == 601) {  // PHOS PHI7
    cuts.AddCutPCMCalo("80062113","00200009f9730000dge0400000","2446600051012200000","0h63103100000010");// PHI7 triggers
    cuts.AddCutPCMCalo("80062113","00200009f9730000dge0400000","2446600011012200000","0h63103100000010");// PHI7 triggers
  } else if (trainConfig == 602){ // No non-lin corr
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","244660007a012200000","0h63103100000010"); // No NL
    cuts.AddCutPCMCalo("80062113","00200009f9730000dge0400000","244660007a012200000","0h63103100000010"); // No NL

  // ===============================================================================================
  // Run 2 data DMC clusters pPb 5TeV
  // ===============================================================================================
  } else if (trainConfig == 700){ // DCal clusters standard cuts, triggers, no nonlin, open timing
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","3885500057032230000","0h63103100000010"); // INT7
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","3885500017032230000","0h63103100000010"); // INT7
  } else if (trainConfig == 701){  // DCal clusters standard cuts, triggers, no nonlin
    cuts.AddCutPCMCalo("00055113","00200009f9730000dge0400000","3885500057032230000","0h63103100000010"); // EMC7
    cuts.AddCutPCMCalo("00089113","00200009f9730000dge0400000","3885500057032230000","0h63103100000010"); // EG2
    cuts.AddCutPCMCalo("0008b113","00200009f9730000dge0400000","3885500057032230000","0h63103100000010"); // EG1
  } else if (trainConfig == 702){  // DCal clusters standard cuts, triggers, no nonlin, open timing
    cuts.AddCutPCMCalo("00055113","00200009f9730000dge0400000","3885500017032230000","0h63103100000010"); // EMC7
    cuts.AddCutPCMCalo("00089113","00200009f9730000dge0400000","3885500017032230000","0h63103100000010"); // EG2
    cuts.AddCutPCMCalo("0008b113","00200009f9730000dge0400000","3885500017032230000","0h63103100000010"); // EG1
  } else if (trainConfig == 703){ // EMCAL clusters standard cuts, triggers, no nonlin, open timing
    cuts.AddCutPCMCalo("80110113","00200009f9730000dge0400000","3885500057032230000","0h63103100000010"); // 0-10
    cuts.AddCutPCMCalo("81210113","00200009f9730000dge0400000","3885500057032230000","0h63103100000010"); // 10-20
    cuts.AddCutPCMCalo("82410113","00200009f9730000dge0400000","3885500057032230000","0h63103100000010"); // 20-40
    cuts.AddCutPCMCalo("84610113","00200009f9730000dge0400000","3885500057032230000","0h63103100000010"); // 40-60
    cuts.AddCutPCMCalo("86810113","00200009f9730000dge0400000","3885500057032230000","0h63103100000010"); // 60-80
    cuts.AddCutPCMCalo("88010113","00200009f9730000dge0400000","3885500057032230000","0h63103100000010"); // 80-100
  } else if (trainConfig == 704){ // EMCAL clusters standard cuts, triggers, no nonlin, open timing
    cuts.AddCutPCMCalo("80210113","00200009f9730000dge0400000","3885500057032230000","0h63103100000010"); // 0-20
    cuts.AddCutPCMCalo("86010113","00200009f9730000dge0400000","3885500057032230000","0h63103100000010"); // 60-100
    cuts.AddCutPCMCalo("a0110113","00200009f9730000dge0400000","3885500057032230000","0h63103100000010"); // 0-5
    cuts.AddCutPCMCalo("a1210113","00200009f9730000dge0400000","3885500057032230000","0h63103100000010"); // 5-10
  } else if (trainConfig == 705){ // EMCAL clusters standard cuts, triggers, no nonlin, open timing
    cuts.AddCutPCMCalo("80110113","00200009f9730000dge0400000","3885500017032230000","0h63103100000010"); // 0-10
    cuts.AddCutPCMCalo("81210113","00200009f9730000dge0400000","3885500017032230000","0h63103100000010"); // 10-20
    cuts.AddCutPCMCalo("82410113","00200009f9730000dge0400000","3885500017032230000","0h63103100000010"); // 20-40
    cuts.AddCutPCMCalo("84610113","00200009f9730000dge0400000","3885500017032230000","0h63103100000010"); // 40-60
    cuts.AddCutPCMCalo("86810113","00200009f9730000dge0400000","3885500017032230000","0h63103100000010"); // 60-80
    cuts.AddCutPCMCalo("88010113","00200009f9730000dge0400000","3885500017032230000","0h63103100000010"); // 80-100
  } else if (trainConfig == 706){ // EMCAL clusters standard cuts, triggers, no nonlin, open timing
    cuts.AddCutPCMCalo("80210113","00200009f9730000dge0400000","3885500017032230000","0h63103100000010"); // 0-20
    cuts.AddCutPCMCalo("86010113","00200009f9730000dge0400000","3885500017032230000","0h63103100000010"); // 60-100
    cuts.AddCutPCMCalo("a0110113","00200009f9730000dge0400000","3885500017032230000","0h63103100000010"); // 0-5
    cuts.AddCutPCMCalo("a1210113","00200009f9730000dge0400000","3885500017032230000","0h63103100000010"); // 5-10
  } else if (trainConfig == 707){ // AOD validation
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","3885500057032230000","0h63103100000010"); // INT7

  } else if (trainConfig == 780){ // DCal clusters standard cuts, CALO+CALOFAST triggers, no nonlin, open timing, no TM
    cuts.AddCutPCMCalo("800a0113","00200009f9730000dge0400000","3885500050032230000","0h63103100000010"); // INT7
    cuts.AddCutPCMCalo("800a6113","00200009f9730000dge0400000","3885500050032230000","0h63103100000010"); // EMC7
    cuts.AddCutPCMCalo("800a7113","00200009f9730000dge0400000","3885500050032230000","0h63103100000010"); // EG2
    cuts.AddCutPCMCalo("800a8113","00200009f9730000dge0400000","3885500050032230000","0h63103100000010"); // EG1
  } else if (trainConfig == 781){ // DCal clusters standard cuts, triggers, no nonlin, open timing, no TM
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","3885500050032230000","0h63103100000010"); // INT7
    cuts.AddCutPCMCalo("80055113","00200009f9730000dge0400000","3885500050032230000","0h63103100000010"); // EMC7
    cuts.AddCutPCMCalo("80089113","00200009f9730000dge0400000","3885500050032230000","0h63103100000010"); // EG2
    cuts.AddCutPCMCalo("8008b113","00200009f9730000dge0400000","3885500050032230000","0h63103100000010"); // EG1
  // ===============================================================================================
  // Run 2 data DMC clusters pPb 8TeV
  // ===============================================================================================
  } else if (trainConfig == 800){  // DCal clusters standard cuts, triggers, no nonlin, open timing
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","3885500017032230000","0h63103100000010"); // INT7
  } else if (trainConfig == 801){  // DCal clusters standard cuts, triggers, no nonlin, open timing
    cuts.AddCutPCMCalo("80089113","00200009f9730000dge0400000","3885500017032230000","0h63103100000010"); // DG2
    cuts.AddCutPCMCalo("8008b113","00200009f9730000dge0400000","3885500017032230000","0h63103100000010"); // DG1
  } else if (trainConfig == 802){  // DCal clusters standard cuts, triggers, no nonlin, open timing
    cuts.AddCutPCMCalo("80099113","00200009f9730000dge0400000","3885500017032230000","0h63103100000010"); // DJ2
    cuts.AddCutPCMCalo("80097113","00200009f9730000dge0400000","3885500017032230000","0h63103100000010"); // DJ1

  } else if (trainConfig == 804){  // DCal clusters standard cuts, triggers, no nonlin, +-50ns timing (5)
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","3885500057032230000","0h63103100000010"); // INT7
  } else if (trainConfig == 805){  // DCal clusters standard cuts, triggers, no nonlin, +-50ns timing (5)
    cuts.AddCutPCMCalo("80089113","00200009f9730000dge0400000","3885500057032230000","0h63103100000010"); // DG2
    cuts.AddCutPCMCalo("8008b113","00200009f9730000dge0400000","3885500057032230000","0h63103100000010"); // DG1
  } else if (trainConfig == 806){  // DCal clusters standard cuts, triggers, no nonlin, +-50ns timing (5)
    cuts.AddCutPCMCalo("80099113","00200009f9730000dge0400000","3885500057032230000","0h63103100000010"); // DJ2
    cuts.AddCutPCMCalo("80097113","00200009f9730000dge0400000","3885500057032230000","0h63103100000010"); // DJ1

  } else if (trainConfig == 807){  // DCal clusters standard cuts, triggers, Nico TB NL, +-50ns timing (5)
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","3885565057032230000","0h63103100000010"); // INT7
  } else if (trainConfig == 808){  // DCal clusters standard cuts, triggers, Nico TB NL, +-50ns timing (5)
    cuts.AddCutPCMCalo("80089113","00200009f9730000dge0400000","3885565057032230000","0h63103100000010"); // DG2
    cuts.AddCutPCMCalo("8008b113","00200009f9730000dge0400000","3885565057032230000","0h63103100000010"); // DG1
  } else if (trainConfig == 809){  // DCal clusters standard cuts, triggers, Nico TB NL, +-50ns timing (5)
    cuts.AddCutPCMCalo("80099113","00200009f9730000dge0400000","3885565057032230000","0h63103100000010"); // DJ2
    cuts.AddCutPCMCalo("80097113","00200009f9730000dge0400000","3885565057032230000","0h63103100000010"); // DJ1

  } else if (trainConfig == 810){  // DCal clusters standard cuts, triggers, Nico TB NL, +-50ns timing (5)
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","3885547057032230000","0h63103100000010"); // INT7
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","3885548057032230000","0h63103100000010"); // INT7
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","3885557057032230000","0h63103100000010"); // INT7
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","3885558057032230000","0h63103100000010"); // INT7
  } else if (trainConfig == 811){  // DCal clusters standard cuts, triggers, Nico TB NL, +-50ns timing (5)
    cuts.AddCutPCMCalo("80089113","00200009f9730000dge0400000","3885547057032230000","0h63103100000010"); // DG2
    cuts.AddCutPCMCalo("80089113","00200009f9730000dge0400000","3885548057032230000","0h63103100000010"); // DG2
    cuts.AddCutPCMCalo("80089113","00200009f9730000dge0400000","3885557057032230000","0h63103100000010"); // DG2
    cuts.AddCutPCMCalo("80089113","00200009f9730000dge0400000","3885558057032230000","0h63103100000010"); // DG2
  } else if (trainConfig == 812){  // DCal clusters standard cuts, triggers, Nico TB NL, +-50ns timing (5)
    cuts.AddCutPCMCalo("8008b113","00200009f9730000dge0400000","3885547057032230000","0h63103100000010"); // DG1
    cuts.AddCutPCMCalo("8008b113","00200009f9730000dge0400000","3885548057032230000","0h63103100000010"); // DG1
    cuts.AddCutPCMCalo("8008b113","00200009f9730000dge0400000","3885557057032230000","0h63103100000010"); // DG1
    cuts.AddCutPCMCalo("8008b113","00200009f9730000dge0400000","3885558057032230000","0h63103100000010"); // DG1

  } else if (trainConfig == 813){  // DCal clusters standard cuts, triggers, Nico TB NL, +-50ns timing (5), TM E/p 1.75 (f)
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","388554705f032230000","0h63103100000010"); // INT7
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","388554805f032230000","0h63103100000010"); // INT7
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","388555705f032230000","0h63103100000010"); // INT7
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","388555805f032230000","0h63103100000010"); // INT7
  } else if (trainConfig == 814){  // DCal clusters standard cuts, triggers, Nico TB NL, +-50ns timing (5), TM E/p 1.75 (f)
    cuts.AddCutPCMCalo("80089113","00200009f9730000dge0400000","388554705f032230000","0h63103100000010"); // DG2
    cuts.AddCutPCMCalo("80089113","00200009f9730000dge0400000","388554805f032230000","0h63103100000010"); // DG2
    cuts.AddCutPCMCalo("80089113","00200009f9730000dge0400000","388555705f032230000","0h63103100000010"); // DG2
    cuts.AddCutPCMCalo("80089113","00200009f9730000dge0400000","388555805f032230000","0h63103100000010"); // DG2
  } else if (trainConfig == 815){  // DCal clusters standard cuts, triggers, Nico TB NL, +-50ns timing (5), TM E/p 1.75 (f)
    cuts.AddCutPCMCalo("8008b113","00200009f9730000dge0400000","388554705f032230000","0h63103100000010"); // DG1
    cuts.AddCutPCMCalo("8008b113","00200009f9730000dge0400000","388554805f032230000","0h63103100000010"); // DG1
    cuts.AddCutPCMCalo("8008b113","00200009f9730000dge0400000","388555705f032230000","0h63103100000010"); // DG1
    cuts.AddCutPCMCalo("8008b113","00200009f9730000dge0400000","388555805f032230000","0h63103100000010"); // DG1

    // Systematics pPb 8 TeV Min Bias
  } else if (trainConfig == 820){ //EMCAL minEnergy variation
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","3885552057032230000","0h63103100000010"); // 0.7 GeV/c default
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","3885552057042230000","0h63103100000010"); // 0.8 GeV/c
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","3885552057052230000","0h63103100000010"); // 0.9 GeV/c
  } else if (trainConfig == 821){ //EMCAL minNCells variation
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","3885552057031230000","0h63103100000010"); // n cells >= 1
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","3885552057033230000","0h63103100000010"); // n cells >= 3
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","3885552057032200000","0h63103100000010"); // no max M02 cut
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","3885552057032250000","0h63103100000010"); // M02 < 0.3
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","3885552057032260000","0h63103100000010"); // M02 < 0.27
  } else if (trainConfig == 822){ // EMCAL track matching variations
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","3885552053032230000","0h63103100000010"); //
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","3885552056032230000","0h63103100000010"); //
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","3885552058032230000","0h63103100000010"); //
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","3885552059032230000","0h63103100000010"); //
  } else if (trainConfig == 823){ // PCM variations
    cuts.AddCutPCMCalo("80010113","00200009227000008250400000","3885552057032230000","0h63103100000010"); // dEdx e -3, 5
    cuts.AddCutPCMCalo("80010113","00200009127000008250400000","3885552057032230000","0h63103100000010"); // dEdx e -5, 5
    cuts.AddCutPCMCalo("80010113","00200009357000008250400000","3885552057032230000","0h63103100000010"); // dEdx pi 2
    cuts.AddCutPCMCalo("80010113","00200009317000008250400000","3885552057032230000","0h63103100000010"); // dEdx pi 0
    cuts.AddCutPCMCalo("80010113","00200009387300008250400000","3885552057032230000","0h63103100000010"); // dEdx pi 2 high 1 (> 3.5 GeV)
  } else if (trainConfig == 824){ // PCM variations
    cuts.AddCutPCMCalo("80010113","00200009327000009250400000","3885552057032230000","0h63103100000010"); // qt 2D 0.03
    cuts.AddCutPCMCalo("80010113","00200009327000003250400000","3885552057032230000","0h63103100000010"); // qt 1D 0.05
    cuts.AddCutPCMCalo("80010113","00200009327000002250400000","3885552057032230000","0h63103100000010"); // qt 1D 0.07
    cuts.AddCutPCMCalo("80010113","00200049327000008250400000","3885552057032230000","0h63103100000010"); // single pt > 0.075
    cuts.AddCutPCMCalo("80010113","00200019327000008250400000","3885552057032230000","0h63103100000010"); // single pt > 0.1
  } else if (trainConfig == 825){ // PCM variations
    cuts.AddCutPCMCalo("80010113","00200009327000008850400000","3885552057032230000","0h63103100000010"); // 2D psi pair chi2 var
    cuts.AddCutPCMCalo("80010113","00200009327000008260400000","3885552057032230000","0h63103100000010"); // 2D psi pair chi2 var
    cuts.AddCutPCMCalo("80010113","00200009327000008860400000","3885552057032230000","0h63103100000010"); // 2D psi pair chi2 var
    cuts.AddCutPCMCalo("80010113","00200009327000008280400000","3885552057032230000","0h63103100000010"); // 2D psi pair chi2 var
    cuts.AddCutPCMCalo("80010113","00200009327000008880400000","3885552057032230000","0h63103100000010"); // 2D psi pair chi2 var
  } else if (trainConfig == 826){ // PCM variations
    cuts.AddCutPCMCalo("80010113","00200006327000008250400000","3885552057032230000","0h63103100000010"); // min TPC cl > 0.7
    cuts.AddCutPCMCalo("80010113","00200008327000008250400000","3885552057032230000","0h63103100000010"); // min TPC cl > 0.35
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","3885552057032230000","0163106100000010"); // alpha < 0.9
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","3885552057032230000","0163105100000010"); // alpha < 0.75
  } else if (trainConfig == 827){ // PCM variations
    cuts.AddCutPCMCalo("80010113","00202209327000008250400000","3885552057032230000","0h63103100000010"); // restrict acceptance to EMCAL loose
    cuts.AddCutPCMCalo("80010113","00204409327000008250400000","3885552057032230000","0h63103100000010"); // restrict acceptance to EMCAL tight
  } else if (trainConfig == 828){ // PCM variations pi dEdx
    cuts.AddCutPCMCalo("80010113","00200009317300008250400000","3885552057032230000","0h63103100000010"); // dEdx pi: 0: 0.4-3.5, -10: 3.5 ->
    cuts.AddCutPCMCalo("80010113","00200009327300008250400000","3885552057032230000","0h63103100000010"); // dEdx pi: 1: 0.4-3.5, -10: 3.5 ->
    cuts.AddCutPCMCalo("80010113","00200009325000008250400000","3885552057032230000","0h63103100000010"); // dEdx pi: 1: 0.3 ->
    cuts.AddCutPCMCalo("80010113","00200009320000008250400000","3885552057032230000","0h63103100000010"); // dEdx pi: 1: 0.5 ->
  } else if (trainConfig == 829){ // PCM variations pi dEdx
    cuts.AddCutPCMCalo("80010113","00200009327600008250400000","3885552057032230000","0h63103100000010"); // dEdx pi: 1: 0.4-2, -10: 2. ->
    cuts.AddCutPCMCalo("80010113","00200009327400008250400000","3885552057032230000","0h63103100000010"); // dEdx pi: 1: 0.4-3, -10: 3. ->
    cuts.AddCutPCMCalo("80010113","00200009315600008250400000","3885552057032230000","0h63103100000010"); // dEdx pi: 0: 0.3-2, -10: 2. ->
    cuts.AddCutPCMCalo("80010113","00200009367400008250400000","3885552057032230000","0h63103100000010"); // dEdx pi: 2: 0.4-3, 0.5: 3. ->
    cuts.AddCutPCMCalo("80010113","00200009347400008250400000","3885552057032230000","0h63103100000010"); // dEdx pi: 3: 0.4-3, 1: 3. ->
  } else if (trainConfig == 830){ // PCM variations to close V0s
    cuts.AddCutPCMCalo("80010113","00200009f9730000dge0400000","3885552057032230000","0h63103100000010");
    cuts.AddCutPCMCalo("80010113","00200009327000008250401000","3885552057032230000","0h63103100000010");
    cuts.AddCutPCMCalo("80010113","00200009327000008250402000","3885552057032230000","0h63103100000010");
    cuts.AddCutPCMCalo("80010113","00200009327000008250403000","3885552057032230000","0h63103100000010");
  // ===============================================================================================
  // Run 1 data EMC clusters pPb 5TeV triggers only
  // ===============================================================================================
  // EMC 7
  } else if (trainConfig == 1000){ // QA setups
    cuts.AddCutPCMCalo("80052113","00200009f9730000dge0400000","1111100057022230000","0h63103100000010");
    cuts.AddCutPCMCalo("80052113","00200009f9730000dge0400000","1111100007032230000","0h63103100000010");
  } else if (trainConfig == 1001){ // default cut
    cuts.AddCutPCMCalo("80052113","00200009f9730000dge0400000","1111141057032230000","0h63103100000010");
    cuts.AddCutPCMCalo("80252113","00200009f9730000dge0400000","1111141057032230000","0h63103100000010");
    cuts.AddCutPCMCalo("82452113","00200009f9730000dge0400000","1111141057032230000","0h63103100000010");
    cuts.AddCutPCMCalo("84652113","00200009f9730000dge0400000","1111141057032230000","0h63103100000010");
    cuts.AddCutPCMCalo("86052113","00200009f9730000dge0400000","1111141057032230000","0h63103100000010");
  } else if (trainConfig == 1002){ //EMCAL minEnergy variation
    cuts.AddCutPCMCalo("80052113","00200009f9730000dge0400000","1111141057032230000","0h63103100000010"); // 0.7 GeV/c default
    cuts.AddCutPCMCalo("80052113","00200009f9730000dge0400000","1111141057042230000","0h63103100000010"); // 0.8 GeV/c
    cuts.AddCutPCMCalo("80052113","00200009f9730000dge0400000","1111141057052230000","0h63103100000010"); // 0.9 GeV/c
  } else if (trainConfig == 1003){ //EMCAL minNCells variation
    cuts.AddCutPCMCalo("80052113","00200009f9730000dge0400000","1111141057031230000","0h63103100000010"); // n cells >= 1
    cuts.AddCutPCMCalo("80052113","00200009f9730000dge0400000","1111141057033230000","0h63103100000010"); // n cells >= 3
    cuts.AddCutPCMCalo("80052113","00200009f9730000dge0400000","1111141057032200000","0h63103100000010"); // no max M02 cut
    cuts.AddCutPCMCalo("80052113","00200009f9730000dge0400000","1111141057032250000","0h63103100000010"); // M02 < 0.3
    cuts.AddCutPCMCalo("80052113","00200009f9730000dge0400000","1111141057032260000","0h63103100000010"); // M02 < 0.27
    cuts.AddCutPCMCalo("80052113","00200009f9730000dge0400000","1112141057032230000","0h63103100000010"); // only modules with TRD infront
    cuts.AddCutPCMCalo("80052113","00200009f9730000dge0400000","1111341057032230000","0h63103100000010"); // no modules with TRD infront
  } else if (trainConfig == 1004){ // EMCAL track matching variations
    cuts.AddCutPCMCalo("80052113","00200009f9730000dge0400000","1111141053032230000","0h63103100000010"); // fixed TM
    cuts.AddCutPCMCalo("80052113","00200009f9730000dge0400000","1111141056032230000","0h63103100000010"); // pt dep var 1
    cuts.AddCutPCMCalo("80052113","00200009f9730000dge0400000","1111141058032230000","0h63103100000010"); // pt dep var 1
    cuts.AddCutPCMCalo("80052113","00200009f9730000dge0400000","1111141059032230000","0h63103100000010"); // pt dep var 1
  } else if (trainConfig == 1005){ // PCM variations
    cuts.AddCutPCMCalo("80052113","00200009227000008250400000","1111141057032230000","0h63103100000010"); // dEdx e -3, 5
    cuts.AddCutPCMCalo("80052113","00200009127000008250400000","1111141057032230000","0h63103100000010"); // dEdx e -5, 5
    cuts.AddCutPCMCalo("80052113","00200009357000008250400000","1111141057032230000","0h63103100000010"); // dEdx pi 2
    cuts.AddCutPCMCalo("80052113","00200009317000008250400000","1111141057032230000","0h63103100000010"); // dEdx pi 0
    cuts.AddCutPCMCalo("80052113","00200009387300008250400000","1111141057032230000","0h63103100000010"); // dEdx pi 2 high 1 (> 3.5 GeV)
  } else if (trainConfig == 1006){ // PCM variations
    cuts.AddCutPCMCalo("80052113","00200009327000009250400000","1111141057032230000","0h63103100000010"); // qt 2D 0.03
    cuts.AddCutPCMCalo("80052113","00200009327000003250400000","1111141057032230000","0h63103100000010"); // qt 1D 0.05
    cuts.AddCutPCMCalo("80052113","00200009327000002250400000","1111141057032230000","0h63103100000010"); // qt 1D 0.07
    cuts.AddCutPCMCalo("80052113","00200049327000008250400000","1111141057032230000","0h63103100000010"); // single pt > 0.075
    cuts.AddCutPCMCalo("80052113","00200019327000008250400000","1111141057032230000","0h63103100000010"); // single pt > 0.1
  } else if (trainConfig == 1007){ // PCM variations
    cuts.AddCutPCMCalo("80052113","00200009327000008850400000","1111141057032230000","0h63103100000010"); // 2D psi pair chi2 var
    cuts.AddCutPCMCalo("80052113","00200009327000008260400000","1111141057032230000","0h63103100000010"); // 2D psi pair chi2 var
    cuts.AddCutPCMCalo("80052113","00200009327000008860400000","1111141057032230000","0h63103100000010"); // 2D psi pair chi2 var
    cuts.AddCutPCMCalo("80052113","00200009327000008280400000","1111141057032230000","0h63103100000010"); // 2D psi pair chi2 var
    cuts.AddCutPCMCalo("80052113","00200009327000008880400000","1111141057032230000","0h63103100000010"); // 2D psi pair chi2 var
  } else if (trainConfig == 1008){ // PCM variations
    cuts.AddCutPCMCalo("80052113","00200006327000008250400000","1111141057032230000","0h63103100000010"); // min TPC cl > 0.7
    cuts.AddCutPCMCalo("80052113","00200008327000008250400000","1111141057032230000","0h63103100000010"); // min TPC cl > 0.35
    cuts.AddCutPCMCalo("80052113","00200009f9730000dge0400000","1111141057032230000","0h63106100000010"); // alpha < 0.8
    cuts.AddCutPCMCalo("80052113","00200009f9730000dge0400000","1111141057032230000","0h63105100000010"); // alpha < 0.75
  } else if (trainConfig == 1009){ // PCM variations
    cuts.AddCutPCMCalo("80052113","00202209327000008250400000","1111141057032230000","0h63103100000010"); // restrict acceptance to EMCAL loose
    cuts.AddCutPCMCalo("80052113","00204409327000008250400000","1111141057032230000","0h63103100000010"); // restrict acceptance to EMCAL tight
  } else if (trainConfig == 1010){ // PCM variations pi dEdx
    cuts.AddCutPCMCalo("80052113","00200009317300008250400000","1111141057032230000","0h63103100000010"); // dEdx pi: 0: 0.4-3.5, -10: 3.5 ->
    cuts.AddCutPCMCalo("80052113","00200009327300008250400000","1111141057032230000","0h63103100000010"); // dEdx pi: 1: 0.4-3.5, -10: 3.5 ->
    cuts.AddCutPCMCalo("80052113","00200009325000008250400000","1111141057032230000","0h63103100000010"); // dEdx pi: 1: 0.3 ->
    cuts.AddCutPCMCalo("80052113","00200009320000008250400000","1111141057032230000","0h63103100000010"); // dEdx pi: 1: 0.5 ->
  } else if (trainConfig == 1011){ // PCM variations pi dEdx
    cuts.AddCutPCMCalo("80052113","00200009327600008250400000","1111141057032230000","0h63103100000010"); // dEdx pi: 1: 0.4-2, -10: 2. ->
    cuts.AddCutPCMCalo("80052113","00200009327400008250400000","1111141057032230000","0h63103100000010"); // dEdx pi: 1: 0.4-3, -10: 3. ->
    cuts.AddCutPCMCalo("80052113","00200009315600008250400000","1111141057032230000","0h63103100000010"); // dEdx pi: 0: 0.3-2, -10: 2. ->
    cuts.AddCutPCMCalo("80052113","00200009367400008250400000","1111141057032230000","0h63103100000010"); // dEdx pi: 2: 0.4-3, 0.5: 3. ->
    cuts.AddCutPCMCalo("80052113","00200009347400008250400000","1111141057032230000","0h63103100000010"); // dEdx pi: 3: 0.4-3, 1: 3. ->
  } else if (trainConfig == 1012){ // PCM variations to close V0s
    cuts.AddCutPCMCalo("80052113","00200009f9730000dge0400000","1111141057032230000","0h63103100000010");
    cuts.AddCutPCMCalo("80052113","00200009327000008250401000","1111141057032230000","0h63103100000010");
    cuts.AddCutPCMCalo("80052113","00200009327000008250402000","1111141057032230000","0h63103100000010");
    cuts.AddCutPCMCalo("80052113","00200009327000008250403000","1111141057032230000","0h63103100000010");
  } else if (trainConfig == 1013){ // EMCal cluster, non lin variations
    cuts.AddCutPCMCalo("80052113","00200009f9730000dge0400000","1111141057032230000","0h63103100000010");
    cuts.AddCutPCMCalo("80052113","00200009f9730000dge0400000","1111142057032230000","0h63103100000010");
    cuts.AddCutPCMCalo("80052113","00200009f9730000dge0400000","1111151057032230000","0h63103100000010");
    cuts.AddCutPCMCalo("80052113","00200009f9730000dge0400000","1111152057032230000","0h63103100000010");
  } else if (trainConfig == 1014){  // EMCal cluster, non lin variations
    cuts.AddCutPCMCalo("80052113","00200009f9730000dge0400000","1111100057032230000","0h63103100000010");
    cuts.AddCutPCMCalo("80052113","00200009f9730000dge0400000","1111101057032230000","0h63103100000010");
    cuts.AddCutPCMCalo("80052113","00200009f9730000dge0400000","1111102057032230000","0h63103100000010");
  } else if (trainConfig == 1015){ // EMCal cluster, non lin variations
    cuts.AddCutPCMCalo("80052113","00200009f9730000dge0400000","1111143057032230000","0h63103100000010");
    cuts.AddCutPCMCalo("80052113","00200009f9730000dge0400000","1111144057032230000","0h63103100000010");
    cuts.AddCutPCMCalo("80052113","00200009f9730000dge0400000","1111153057032230000","0h63103100000010");
    cuts.AddCutPCMCalo("80052113","00200009f9730000dge0400000","1111154057032230000","0h63103100000010");
  } else if (trainConfig == 1016){ // default cut
    cuts.AddCutPCMCalo("80052113","00200009f9730000dge0400000","1111141057032230000","0i63103100000010");
    cuts.AddCutPCMCalo("80252113","00200009f9730000dge0400000","1111141057032230000","0i63103100000010");
    cuts.AddCutPCMCalo("82452113","00200009f9730000dge0400000","1111141057032230000","0i63103100000010");
    cuts.AddCutPCMCalo("84652113","00200009f9730000dge0400000","1111141057032230000","0i63103100000010");
    cuts.AddCutPCMCalo("86052113","00200009f9730000dge0400000","1111141057032230000","0i63103100000010");

    // EG2 triggers
  } else if (trainConfig == 1030){ // QA setups
    cuts.AddCutPCMCalo("80085113","00200009f9730000dge0400000","1111100057022230000","0h63103100000010");
    cuts.AddCutPCMCalo("80085113","00200009f9730000dge0400000","1111100007032230000","0h63103100000010");
  } else if (trainConfig == 1031){ // default cut
    cuts.AddCutPCMCalo("80085113","00200009f9730000dge0400000","1111141057032230000","0h63103100000010");
    cuts.AddCutPCMCalo("80285113","00200009f9730000dge0400000","1111141057032230000","0h63103100000010");
    cuts.AddCutPCMCalo("82485113","00200009f9730000dge0400000","1111141057032230000","0h63103100000010");
    cuts.AddCutPCMCalo("84685113","00200009f9730000dge0400000","1111141057032230000","0h63103100000010");
    cuts.AddCutPCMCalo("86085113","00200009f9730000dge0400000","1111141057032230000","0h63103100000010");
  } else if (trainConfig == 1032){ //EMCAL minEnergy variation
    cuts.AddCutPCMCalo("80085113","00200009f9730000dge0400000","1111141057032230000","0h63103100000010"); // 0.7 GeV/c default
    cuts.AddCutPCMCalo("80085113","00200009f9730000dge0400000","1111141057042230000","0h63103100000010"); // 0.8 GeV/c
    cuts.AddCutPCMCalo("80085113","00200009f9730000dge0400000","1111141057052230000","0h63103100000010"); // 0.9 GeV/c
  } else if (trainConfig == 1033){ //EMCAL minNCells variation
    cuts.AddCutPCMCalo("80085113","00200009f9730000dge0400000","1111141057031230000","0h63103100000010"); // n cells >= 1
    cuts.AddCutPCMCalo("80085113","00200009f9730000dge0400000","1111141057033230000","0h63103100000010"); // n cells >= 3
    cuts.AddCutPCMCalo("80085113","00200009f9730000dge0400000","1111141057032200000","0h63103100000010"); // no max M02 cut
    cuts.AddCutPCMCalo("80085113","00200009f9730000dge0400000","1111141057032250000","0h63103100000010"); // M02 < 0.3
    cuts.AddCutPCMCalo("80085113","00200009f9730000dge0400000","1111141057032260000","0h63103100000010"); // M02 < 0.27
    cuts.AddCutPCMCalo("80085113","00200009f9730000dge0400000","1112141057032230000","0h63103100000010"); // only modules with TRD infront
    cuts.AddCutPCMCalo("80085113","00200009f9730000dge0400000","1111341057032230000","0h63103100000010"); // no modules with TRD infront
  } else if (trainConfig == 1034){ // EMCAL track matching variations
    cuts.AddCutPCMCalo("80085113","00200009f9730000dge0400000","1111141053032230000","0h63103100000010"); // fixed TM
    cuts.AddCutPCMCalo("80085113","00200009f9730000dge0400000","1111141056032230000","0h63103100000010"); // pt dep var 1
    cuts.AddCutPCMCalo("80085113","00200009f9730000dge0400000","1111141058032230000","0h63103100000010"); // pt dep var 1
    cuts.AddCutPCMCalo("80085113","00200009f9730000dge0400000","1111141059032230000","0h63103100000010"); // pt dep var 1
  } else if (trainConfig == 1035){ // PCM variations
    cuts.AddCutPCMCalo("80085113","00200009227000008250400000","1111141057032230000","0h63103100000010"); // dEdx e -3, 5
    cuts.AddCutPCMCalo("80085113","00200009127000008250400000","1111141057032230000","0h63103100000010"); // dEdx e -5, 5
    cuts.AddCutPCMCalo("80085113","00200009357000008250400000","1111141057032230000","0h63103100000010"); // dEdx pi 2
    cuts.AddCutPCMCalo("80085113","00200009317000008250400000","1111141057032230000","0h63103100000010"); // dEdx pi 0
    cuts.AddCutPCMCalo("80085113","00200009387300008250400000","1111141057032230000","0h63103100000010"); // dEdx pi 2 high 1 (> 3.5 GeV)
  } else if (trainConfig == 1036){ // PCM variations
    cuts.AddCutPCMCalo("80085113","00200009327000009250400000","1111141057032230000","0h63103100000010"); // qt 2D 0.03
    cuts.AddCutPCMCalo("80085113","00200009327000003250400000","1111141057032230000","0h63103100000010"); // qt 1D 0.05
    cuts.AddCutPCMCalo("80085113","00200009327000002250400000","1111141057032230000","0h63103100000010"); // qt 1D 0.07
    cuts.AddCutPCMCalo("80085113","00200049327000008250400000","1111141057032230000","0h63103100000010"); // single pt > 0.075
    cuts.AddCutPCMCalo("80085113","00200019327000008250400000","1111141057032230000","0h63103100000010"); // single pt > 0.1
  } else if (trainConfig == 1037){ // PCM variations
    cuts.AddCutPCMCalo("80085113","00200009327000008850400000","1111141057032230000","0h63103100000010"); // 2D psi pair chi2 var
    cuts.AddCutPCMCalo("80085113","00200009327000008260400000","1111141057032230000","0h63103100000010"); // 2D psi pair chi2 var
    cuts.AddCutPCMCalo("80085113","00200009327000008860400000","1111141057032230000","0h63103100000010"); // 2D psi pair chi2 var
    cuts.AddCutPCMCalo("80085113","00200009327000008280400000","1111141057032230000","0h63103100000010"); // 2D psi pair chi2 var
    cuts.AddCutPCMCalo("80085113","00200009327000008880400000","1111141057032230000","0h63103100000010"); // 2D psi pair chi2 var
  } else if (trainConfig == 1038){ // PCM variations
    cuts.AddCutPCMCalo("80085113","00200006327000008250400000","1111141057032230000","0h63103100000010"); // min TPC cl > 0.7
    cuts.AddCutPCMCalo("80085113","00200008327000008250400000","1111141057032230000","0h63103100000010"); // min TPC cl > 0.35
    cuts.AddCutPCMCalo("80085113","00200009f9730000dge0400000","1111141057032230000","0h63106100000010"); // alpha < 0.8
    cuts.AddCutPCMCalo("80085113","00200009f9730000dge0400000","1111141057032230000","0h63105100000010"); // alpha < 0.75
  } else if (trainConfig == 1039){ // PCM variations
    cuts.AddCutPCMCalo("80085113","00202209327000008250400000","1111141057032230000","0h63103100000010"); // restrict acceptance to EMCAL loose
    cuts.AddCutPCMCalo("80085113","00204409327000008250400000","1111141057032230000","0h63103100000010"); // restrict acceptance to EMCAL tight
  } else if (trainConfig == 1040){ // PCM variations pi dEdx
    cuts.AddCutPCMCalo("80085113","00200009317300008250400000","1111141057032230000","0h63103100000010"); // dEdx pi: 0: 0.4-3.5, -10: 3.5 ->
    cuts.AddCutPCMCalo("80085113","00200009327300008250400000","1111141057032230000","0h63103100000010"); // dEdx pi: 1: 0.4-3.5, -10: 3.5 ->
    cuts.AddCutPCMCalo("80085113","00200009325000008250400000","1111141057032230000","0h63103100000010"); // dEdx pi: 1: 0.3 ->
    cuts.AddCutPCMCalo("80085113","00200009320000008250400000","1111141057032230000","0h63103100000010"); // dEdx pi: 1: 0.5 ->
  } else if (trainConfig == 1041){ // PCM variations pi dEdx
    cuts.AddCutPCMCalo("80085113","00200009327600008250400000","1111141057032230000","0h63103100000010"); // dEdx pi: 1: 0.4-2, -10: 2. ->
    cuts.AddCutPCMCalo("80085113","00200009327400008250400000","1111141057032230000","0h63103100000010"); // dEdx pi: 1: 0.4-3, -10: 3. ->
    cuts.AddCutPCMCalo("80085113","00200009315600008250400000","1111141057032230000","0h63103100000010"); // dEdx pi: 0: 0.3-2, -10: 2. ->
    cuts.AddCutPCMCalo("80085113","00200009367400008250400000","1111141057032230000","0h63103100000010"); // dEdx pi: 2: 0.4-3, 0.5: 3. ->
    cuts.AddCutPCMCalo("80085113","00200009347400008250400000","1111141057032230000","0h63103100000010"); // dEdx pi: 3: 0.4-3, 1: 3. ->
  } else if (trainConfig == 1042){ // PCM variations to close V0s
    cuts.AddCutPCMCalo("80085113","00200009f9730000dge0400000","1111141057032230000","0h63103100000010");
    cuts.AddCutPCMCalo("80085113","00200009327000008250401000","1111141057032230000","0h63103100000010");
    cuts.AddCutPCMCalo("80085113","00200009327000008250402000","1111141057032230000","0h63103100000010");
    cuts.AddCutPCMCalo("80085113","00200009327000008250403000","1111141057032230000","0h63103100000010");
  } else if (trainConfig == 1043){ // EMCal cluster, non lin variations
    cuts.AddCutPCMCalo("80085113","00200009f9730000dge0400000","1111141057032230000","0h63103100000010");
    cuts.AddCutPCMCalo("80085113","00200009f9730000dge0400000","1111142057032230000","0h63103100000010");
    cuts.AddCutPCMCalo("80085113","00200009f9730000dge0400000","1111151057032230000","0h63103100000010");
    cuts.AddCutPCMCalo("80085113","00200009f9730000dge0400000","1111152057032230000","0h63103100000010");
  } else if (trainConfig == 1044){  // EMCal cluster, non lin variations
    cuts.AddCutPCMCalo("80085113","00200009f9730000dge0400000","1111100057032230000","0h63103100000010");
    cuts.AddCutPCMCalo("80085113","00200009f9730000dge0400000","1111101057032230000","0h63103100000010");
    cuts.AddCutPCMCalo("80085113","00200009f9730000dge0400000","1111102057032230000","0h63103100000010");
  } else if (trainConfig == 1045){ // EMCal cluster, non lin variations
    cuts.AddCutPCMCalo("80085113","00200009f9730000dge0400000","1111143057032230000","0h63103100000010");
    cuts.AddCutPCMCalo("80085113","00200009f9730000dge0400000","1111144057032230000","0h63103100000010");
    cuts.AddCutPCMCalo("80085113","00200009f9730000dge0400000","1111153057032230000","0h63103100000010");
    cuts.AddCutPCMCalo("80085113","00200009f9730000dge0400000","1111154057032230000","0h63103100000010");
  } else if (trainConfig == 1046){ // default cut with sector mixing
    cuts.AddCutPCMCalo("80085113","00200009f9730000dge0400000","1111141057032230000","0i63103100000010");
    cuts.AddCutPCMCalo("80285113","00200009f9730000dge0400000","1111141057032230000","0i63103100000010");
    cuts.AddCutPCMCalo("82485113","00200009f9730000dge0400000","1111141057032230000","0i63103100000010");
    cuts.AddCutPCMCalo("84685113","00200009f9730000dge0400000","1111141057032230000","0i63103100000010");
    cuts.AddCutPCMCalo("86085113","00200009f9730000dge0400000","1111141057032230000","0i63103100000010");

    // EG1 triggers
  } else if (trainConfig == 1060){ // QA setups
    cuts.AddCutPCMCalo("80083113","00200009f9730000dge0400000","1111100057022230000","0h63103100000010");
    cuts.AddCutPCMCalo("80083113","00200009f9730000dge0400000","1111100007032230000","0h63103100000010");
  } else if (trainConfig == 1061){ // default cut
    cuts.AddCutPCMCalo("80083113","00200009f9730000dge0400000","1111141057032230000","0h63103100000010");
    cuts.AddCutPCMCalo("80283113","00200009f9730000dge0400000","1111141057032230000","0h63103100000010");
    cuts.AddCutPCMCalo("82483113","00200009f9730000dge0400000","1111141057032230000","0h63103100000010");
    cuts.AddCutPCMCalo("84683113","00200009f9730000dge0400000","1111141057032230000","0h63103100000010");
    cuts.AddCutPCMCalo("86083113","00200009f9730000dge0400000","1111141057032230000","0h63103100000010");
  } else if (trainConfig == 1062){ //EMCAL minEnergy variation
    cuts.AddCutPCMCalo("80083113","00200009f9730000dge0400000","1111141057032230000","0h63103100000010"); // 0.7 GeV/c default
    cuts.AddCutPCMCalo("80083113","00200009f9730000dge0400000","1111141057042230000","0h63103100000010"); // 0.8 GeV/c
    cuts.AddCutPCMCalo("80083113","00200009f9730000dge0400000","1111141057052230000","0h63103100000010"); // 0.9 GeV/c
  } else if (trainConfig == 1063){ //EMCAL minNCells variation
    cuts.AddCutPCMCalo("80083113","00200009f9730000dge0400000","1111141057031230000","0h63103100000010"); // n cells >= 1
    cuts.AddCutPCMCalo("80083113","00200009f9730000dge0400000","1111141057033230000","0h63103100000010"); // n cells >= 3
    cuts.AddCutPCMCalo("80083113","00200009f9730000dge0400000","1111141057032200000","0h63103100000010"); // no max M02 cut
    cuts.AddCutPCMCalo("80083113","00200009f9730000dge0400000","1111141057032250000","0h63103100000010"); // M02 < 0.3
    cuts.AddCutPCMCalo("80083113","00200009f9730000dge0400000","1111141057032260000","0h63103100000010"); // M02 < 0.27
    cuts.AddCutPCMCalo("80083113","00200009f9730000dge0400000","1112141057032230000","0h63103100000010"); // only modules with TRD infront
    cuts.AddCutPCMCalo("80083113","00200009f9730000dge0400000","1111341057032230000","0h63103100000010"); // no modules with TRD infront
  } else if (trainConfig == 1064){ // EMCAL track matching variations
    cuts.AddCutPCMCalo("80083113","00200009f9730000dge0400000","1111141053032230000","0h63103100000010"); // fixed TM
    cuts.AddCutPCMCalo("80083113","00200009f9730000dge0400000","1111141056032230000","0h63103100000010"); // pt dep var 1
    cuts.AddCutPCMCalo("80083113","00200009f9730000dge0400000","1111141058032230000","0h63103100000010"); // pt dep var 1
    cuts.AddCutPCMCalo("80083113","00200009f9730000dge0400000","1111141059032230000","0h63103100000010"); // pt dep var 1
  } else if (trainConfig == 1065){ // PCM variations
    cuts.AddCutPCMCalo("80083113","00200009227000008250400000","1111141057032230000","0h63103100000010"); // dEdx e -3, 5
    cuts.AddCutPCMCalo("80083113","00200009127000008250400000","1111141057032230000","0h63103100000010"); // dEdx e -5, 5
    cuts.AddCutPCMCalo("80083113","00200009357000008250400000","1111141057032230000","0h63103100000010"); // dEdx pi 2
    cuts.AddCutPCMCalo("80083113","00200009317000008250400000","1111141057032230000","0h63103100000010"); // dEdx pi 0
    cuts.AddCutPCMCalo("80083113","00200009387300008250400000","1111141057032230000","0h63103100000010"); // dEdx pi 2 high 1 (> 3.5 GeV)
  } else if (trainConfig == 1066){ // PCM variations
    cuts.AddCutPCMCalo("80083113","00200009327000009250400000","1111141057032230000","0h63103100000010"); // qt 2D 0.03
    cuts.AddCutPCMCalo("80083113","00200009327000003250400000","1111141057032230000","0h63103100000010"); // qt 1D 0.05
    cuts.AddCutPCMCalo("80083113","00200009327000002250400000","1111141057032230000","0h63103100000010"); // qt 1D 0.07
    cuts.AddCutPCMCalo("80083113","00200049327000008250400000","1111141057032230000","0h63103100000010"); // single pt > 0.075
    cuts.AddCutPCMCalo("80083113","00200019327000008250400000","1111141057032230000","0h63103100000010"); // single pt > 0.1
  } else if (trainConfig == 1067){ // PCM variations
    cuts.AddCutPCMCalo("80083113","00200009327000008850400000","1111141057032230000","0h63103100000010"); // 2D psi pair chi2 var
    cuts.AddCutPCMCalo("80083113","00200009327000008260400000","1111141057032230000","0h63103100000010"); // 2D psi pair chi2 var
    cuts.AddCutPCMCalo("80083113","00200009327000008860400000","1111141057032230000","0h63103100000010"); // 2D psi pair chi2 var
    cuts.AddCutPCMCalo("80083113","00200009327000008280400000","1111141057032230000","0h63103100000010"); // 2D psi pair chi2 var
    cuts.AddCutPCMCalo("80083113","00200009327000008880400000","1111141057032230000","0h63103100000010"); // 2D psi pair chi2 var
  } else if (trainConfig == 1068){ // PCM variations
    cuts.AddCutPCMCalo("80083113","00200006327000008250400000","1111141057032230000","0h63103100000010"); // min TPC cl > 0.7
    cuts.AddCutPCMCalo("80083113","00200008327000008250400000","1111141057032230000","0h63103100000010"); // min TPC cl > 0.35
    cuts.AddCutPCMCalo("80083113","00200009f9730000dge0400000","1111141057032230000","0h63106100000010"); // alpha < 0.8
    cuts.AddCutPCMCalo("80083113","00200009f9730000dge0400000","1111141057032230000","0h63105100000010"); // alpha < 0.75
  } else if (trainConfig == 1069){ // PCM variations
    cuts.AddCutPCMCalo("80083113","00202209327000008250400000","1111141057032230000","0h63103100000010"); // restrict acceptance to EMCAL loose
    cuts.AddCutPCMCalo("80083113","00204409327000008250400000","1111141057032230000","0h63103100000010"); // restrict acceptance to EMCAL tight
  } else if (trainConfig == 1070){ // PCM variations pi dEdx
    cuts.AddCutPCMCalo("80083113","00200009317300008250400000","1111141057032230000","0h63103100000010"); // dEdx pi: 0: 0.4-3.5, -10: 3.5 ->
    cuts.AddCutPCMCalo("80083113","00200009327300008250400000","1111141057032230000","0h63103100000010"); // dEdx pi: 1: 0.4-3.5, -10: 3.5 ->
    cuts.AddCutPCMCalo("80083113","00200009325000008250400000","1111141057032230000","0h63103100000010"); // dEdx pi: 1: 0.3 ->
    cuts.AddCutPCMCalo("80083113","00200009320000008250400000","1111141057032230000","0h63103100000010"); // dEdx pi: 1: 0.5 ->
  } else if (trainConfig == 1071){ // PCM variations pi dEdx
    cuts.AddCutPCMCalo("80083113","00200009327600008250400000","1111141057032230000","0h63103100000010"); // dEdx pi: 1: 0.4-2, -10: 2. ->
    cuts.AddCutPCMCalo("80083113","00200009327400008250400000","1111141057032230000","0h63103100000010"); // dEdx pi: 1: 0.4-3, -10: 3. ->
    cuts.AddCutPCMCalo("80083113","00200009315600008250400000","1111141057032230000","0h63103100000010"); // dEdx pi: 0: 0.3-2, -10: 2. ->
    cuts.AddCutPCMCalo("80083113","00200009367400008250400000","1111141057032230000","0h63103100000010"); // dEdx pi: 2: 0.4-3, 0.5: 3. ->
    cuts.AddCutPCMCalo("80083113","00200009347400008250400000","1111141057032230000","0h63103100000010"); // dEdx pi: 3: 0.4-3, 1: 3. ->
  } else if (trainConfig == 1072){ // PCM variations to close V0s
    cuts.AddCutPCMCalo("80083113","00200009f9730000dge0400000","1111141057032230000","0h63103100000010");
    cuts.AddCutPCMCalo("80083113","00200009327000008250401000","1111141057032230000","0h63103100000010");
    cuts.AddCutPCMCalo("80083113","00200009327000008250402000","1111141057032230000","0h63103100000010");
    cuts.AddCutPCMCalo("80083113","00200009327000008250403000","1111141057032230000","0h63103100000010");
  } else if (trainConfig == 1073){ // EMCal cluster, non lin variations
    cuts.AddCutPCMCalo("80083113","00200009f9730000dge0400000","1111141057032230000","0h63103100000010");
    cuts.AddCutPCMCalo("80083113","00200009f9730000dge0400000","1111142057032230000","0h63103100000010");
    cuts.AddCutPCMCalo("80083113","00200009f9730000dge0400000","1111151057032230000","0h63103100000010");
    cuts.AddCutPCMCalo("80083113","00200009f9730000dge0400000","1111152057032230000","0h63103100000010");
  } else if (trainConfig == 1074){  // EMCal cluster, non lin variations
    cuts.AddCutPCMCalo("80083113","00200009f9730000dge0400000","1111100057032230000","0h63103100000010");
    cuts.AddCutPCMCalo("80083113","00200009f9730000dge0400000","1111101057032230000","0h63103100000010");
    cuts.AddCutPCMCalo("80083113","00200009f9730000dge0400000","1111102057032230000","0h63103100000010");
  } else if (trainConfig == 1075){ // EMCal cluster, non lin variations
    cuts.AddCutPCMCalo("80083113","00200009f9730000dge0400000","1111143057032230000","0h63103100000010");
    cuts.AddCutPCMCalo("80083113","00200009f9730000dge0400000","1111144057032230000","0h63103100000010");
    cuts.AddCutPCMCalo("80083113","00200009f9730000dge0400000","1111153057032230000","0h63103100000010");
    cuts.AddCutPCMCalo("80083113","00200009f9730000dge0400000","1111154057032230000","0h63103100000010");
  } else if (trainConfig == 1076){ // default cut with sector mixing
    cuts.AddCutPCMCalo("80083113","00200009f9730000dge0400000","1111141057032230000","0i63103100000010");
    cuts.AddCutPCMCalo("80283113","00200009f9730000dge0400000","1111141057032230000","0i63103100000010");
    cuts.AddCutPCMCalo("82483113","00200009f9730000dge0400000","1111141057032230000","0i63103100000010");
    cuts.AddCutPCMCalo("84683113","00200009f9730000dge0400000","1111141057032230000","0i63103100000010");
    cuts.AddCutPCMCalo("86083113","00200009f9730000dge0400000","1111141057032230000","0i63103100000010");


  // configurations for pPb 5 and 8 TeV Run2 with EMCAL + DCAL
  } else if (trainConfig == 2000){  // EMCal+DCAL clusters standard cuts, triggers, Nico TB NL, +-50ns timing (5)
    cuts.AddCutPCMCalo("80010123","00200009f9730000dge0400000","4117947057032230000","0h63103100000010"); // INT7
    cuts.AddCutPCMCalo("80010123","00200009f9730000dge0400000","4117948057032230000","0h63103100000010"); // INT7
    cuts.AddCutPCMCalo("80010123","00200009f9730000dge0400000","4117957057032230000","0h63103100000010"); // INT7
    cuts.AddCutPCMCalo("80010123","00200009f9730000dge0400000","4117958057032230000","0h63103100000010"); // INT7
  } else if (trainConfig == 2001){  // EMCal+DCAL clusters standard cuts, triggers, Nico TB NL, +-50ns timing (5)
    cuts.AddCutPCMCalo("8008e123","00200009f9730000dge0400000","4117947057032230000","0h63103100000010"); // EG2+DG2
    cuts.AddCutPCMCalo("8008e123","00200009f9730000dge0400000","4117948057032230000","0h63103100000010"); // EG2+DG2
    cuts.AddCutPCMCalo("8008e123","00200009f9730000dge0400000","4117957057032230000","0h63103100000010"); // EG2+DG2
    cuts.AddCutPCMCalo("8008e123","00200009f9730000dge0400000","4117958057032230000","0h63103100000010"); // EG2+DG2
  } else if (trainConfig == 2002){  // EMCal+DCAL clusters standard cuts, triggers, Nico TB NL, +-50ns timing (5)
    cuts.AddCutPCMCalo("8008d123","00200009f9730000dge0400000","4117947057032230000","0h63103100000010"); // EG1+DG1
    cuts.AddCutPCMCalo("8008d123","00200009f9730000dge0400000","4117948057032230000","0h63103100000010"); // EG1+DG1
    cuts.AddCutPCMCalo("8008d123","00200009f9730000dge0400000","4117957057032230000","0h63103100000010"); // EG1+DG1
    cuts.AddCutPCMCalo("8008d123","00200009f9730000dge0400000","4117958057032230000","0h63103100000010"); // EG1+DG1
  } else if (trainConfig == 2003){  // EMCal+DCal clusters standard cuts, triggers, Nico TB NL, +-50ns timing (5)
    cuts.AddCutPCMCalo("80010123","00200009f9730000dge0400000","4117958057032230000","0h63103100000010"); // INT7
    cuts.AddCutPCMCalo("8008e123","00200009f9730000dge0400000","4117958057032230000","0h63103100000010"); // EG2+DG2
    cuts.AddCutPCMCalo("8008d123","00200009f9730000dge0400000","4117958057032230000","0h63103100000010"); // EG1+DG1
  } else if (trainConfig == 2004){  // EMCal clusters standard cuts, triggers, Nico TB NL, +-50ns timing (5)
    cuts.AddCutPCMCalo("80010123","00200009f9730000dge0400000","1111158057032230000","0h63103100000010"); // INT7
    cuts.AddCutPCMCalo("80085123","00200009f9730000dge0400000","1111158057032230000","0h63103100000010"); // EG2
    cuts.AddCutPCMCalo("80083123","00200009f9730000dge0400000","1111158057032230000","0h63103100000010"); // EG1
  } else if (trainConfig == 2005){  // DCal clusters standard cuts, triggers, Nico TB NL, +-50ns timing (5)
    cuts.AddCutPCMCalo("80010123","00200009f9730000dge0400000","3885558057032230000","0h63103100000010"); // INT7
    cuts.AddCutPCMCalo("8008b123","00200009f9730000dge0400000","3885558057032230000","0h63103100000010"); // DG2
    cuts.AddCutPCMCalo("80089123","00200009f9730000dge0400000","3885558057032230000","0h63103100000010"); // DG1
  } else if (trainConfig == 2006){  // EMCal+DCAL clusters standard cuts, triggers, Nico TB NL and no NL, +-50ns timing (5)
    cuts.AddCutPCMCalo("80010123","00200009f9730000dge0400000","4117900057032230000","0h63103100000010"); // INT7
    cuts.AddCutPCMCalo("80010123","00200009f9730000dge0400000","4117906057032230000","0h63103100000010"); // INT7
    cuts.AddCutPCMCalo("80010123","00200009f9730000dge0400000","4117965057032230000","0h63103100000010"); // INT7
  } else if (trainConfig == 2007){  // EMCal+DCAL clusters standard cuts, triggers, Nico TB NL, +-50ns timing (5)
    cuts.AddCutPCMCalo("8008e123","00200009f9730000dge0400000","4117900057032230000","0h63103100000010"); // EG2+DG2
    cuts.AddCutPCMCalo("8008e123","00200009f9730000dge0400000","4117906057032230000","0h63103100000010"); // EG2+DG2
    cuts.AddCutPCMCalo("8008e123","00200009f9730000dge0400000","4117965057032230000","0h63103100000010"); // EG2+DG2
  } else if (trainConfig == 2008){  // EMCal+DCAL clusters standard cuts, triggers, Nico TB NL, +-50ns timing (5)
    cuts.AddCutPCMCalo("8008d123","00200009f9730000dge0400000","4117900057032230000","0h63103100000010"); // EG1+DG1
    cuts.AddCutPCMCalo("8008d123","00200009f9730000dge0400000","4117906057032230000","0h63103100000010"); // EG1+DG1
    cuts.AddCutPCMCalo("8008d123","00200009f9730000dge0400000","4117965057032230000","0h63103100000010"); // EG1+DG1
  } else if (trainConfig == 2010){  // EMCal+DCAL clusters standard cuts, triggers, Nico TB NL and no NL, +-50ns timing (5)
    cuts.AddCutPCMCalo("80010123","00200009f9730000dge0400000","4117900057032230000","0h63103100000010"); // INT7
    cuts.AddCutPCMCalo("80010123","00200009f9730000dge0400000","4117906057032230000","0h63103100000010"); // INT7
    cuts.AddCutPCMCalo("80010123","00200009f9730000dge0400000","4117965057032230000","0h63103100000010"); // INT7
  } else if (trainConfig == 2011){  // EMCal+DCAL clusters standard cuts, triggers, Nico TB NL, +-50ns timing (5)
    cuts.AddCutPCMCalo("8008e123","00200009f9730000dge0400000","4117900057032230000","0h63103100000010"); // EG2+DG2
    cuts.AddCutPCMCalo("8008e123","00200009f9730000dge0400000","4117906057032230000","0h63103100000010"); // EG2+DG2
    cuts.AddCutPCMCalo("8008e123","00200009f9730000dge0400000","4117965057032230000","0h63103100000010"); // EG2+DG2
  } else if (trainConfig == 2012){  // EMCal+DCAL clusters standard cuts, triggers, Nico TB NL, +-50ns timing (5)
    cuts.AddCutPCMCalo("8008d123","00200009f9730000dge0400000","4117900057032230000","0h63103100000010"); // EG1+DG1
    cuts.AddCutPCMCalo("8008d123","00200009f9730000dge0400000","4117906057032230000","0h63103100000010"); // EG1+DG1
    cuts.AddCutPCMCalo("8008d123","00200009f9730000dge0400000","4117965057032230000","0h63103100000010"); // EG1+DG1
  } else if (trainConfig == 2013){  // EMCal+DCAL clusters standard cuts, triggers, Nico TB NL + FineTuning 57, +-50ns timing (5)
    cuts.AddCutPCMCalo("80010123","00200009f9730000dge0400000","4117900057032230000","0h63103100000010"); // INT7
    cuts.AddCutPCMCalo("8008e123","00200009f9730000dge0400000","4117900057032230000","0h63103100000010"); // EG2+DG2
    cuts.AddCutPCMCalo("8008d123","00200009f9730000dge0400000","4117900057032230000","0h63103100000010"); // EG1+DG1
  } else if (trainConfig == 2014){  // EMCal+DCAL clusters standard cuts, triggers, Nico TB NL + FineTuning 57, +-50ns timing (5)
    cuts.AddCutPCMCalo("80010123","00200009f9730000dge0400000","411790005f032230000","0h63103100000010"); // INT7
    cuts.AddCutPCMCalo("8008e123","00200009f9730000dge0400000","411790005f032230000","0h63103100000010"); // EG2+DG2
    cuts.AddCutPCMCalo("8008d123","00200009f9730000dge0400000","411790005f032230000","0h63103100000010"); // EG1+DG1
  } else if (trainConfig == 2015){  // EMCal+DCAL clusters standard cuts, triggers, Nico TB NL + FineTuning 57, +-50ns timing (5)
    cuts.AddCutPCMCalo("80010123","00200009f9730000dge0400000","411796505f032230000","0h63103100000010"); // INT7
    cuts.AddCutPCMCalo("8008e123","00200009f9730000dge0400000","411796505f032230000","0h63103100000010"); // EG2+DG2
    cuts.AddCutPCMCalo("8008d123","00200009f9730000dge0400000","411796505f032230000","0h63103100000010"); // EG1+DG1
  } else if (trainConfig == 2016){  // EMCal+DCAL clusters standard cuts, triggers, Nico TB NL + FineTuning 57, +-50ns timing (5)
    cuts.AddCutPCMCalo("80010123","00200009f9730000dge0400000","4117957057032230000","0h63103100000010"); // INT7
    cuts.AddCutPCMCalo("8008e123","00200009f9730000dge0400000","4117957057032230000","0h63103100000010"); // EG2+DG2
    cuts.AddCutPCMCalo("8008d123","00200009f9730000dge0400000","4117957057032230000","0h63103100000010"); // EG1+DG1
  } else if (trainConfig == 2017){  // EMCal+DCAL clusters standard cuts, triggers, Nico TB NL + FineTuning 57, +-50ns timing (5)
    cuts.AddCutPCMCalo("80010123","00200009f9730000dge0400000","411795705f032230000","0h63103100000010"); // INT7
    cuts.AddCutPCMCalo("8008e123","00200009f9730000dge0400000","411795705f032230000","0h63103100000010"); // EG2+DG2
    cuts.AddCutPCMCalo("8008d123","00200009f9730000dge0400000","411795705f032230000","0h63103100000010"); // EG1+DG1

  } else if (trainConfig == 2018){  // EMCal+DCAL clusters standard cuts, triggers, Nico TB NL + FineTuning 57, +-50ns timing (5)
    cuts.AddCutPCMCalo("80010123","00200009f9730000dge0400000","411790105f032230000","0h63103100000010"); // INT7
    cuts.AddCutPCMCalo("8008e123","00200009f9730000dge0400000","411790105f032230000","0h63103100000010"); // EG2+DG2
    cuts.AddCutPCMCalo("8008d123","00200009f9730000dge0400000","411790105f032230000","0h63103100000010"); // EG1+DG1
  } else if (trainConfig == 2019){  // EMCal+DCAL clusters standard cuts, triggers, Nico TB NL + FineTuning 57, +-50ns timing (5)
    cuts.AddCutPCMCalo("80010123","00200009f9730000dge0400000","411790105f022230000","0h63103100000010"); // INT7
    cuts.AddCutPCMCalo("8008e123","00200009f9730000dge0400000","411790105f022230000","0h63103100000010"); // EG2+DG2
    cuts.AddCutPCMCalo("8008d123","00200009f9730000dge0400000","411790105f022230000","0h63103100000010"); // EG1+DG1

  } else if (trainConfig == 2020){  // EMCal+DCAL clusters standard cuts, triggers, Nico TB NL + FineTuning 57, +-50ns timing (5)
    cuts.AddCutPCMCalo("80010123","00200009f9730000dge0400000","411796105f032230000","0h63103100000010"); // INT7
    cuts.AddCutPCMCalo("8008e123","00200009f9730000dge0400000","411796105f032230000","0h63103100000010"); // EG2+DG2
    cuts.AddCutPCMCalo("8008d123","00200009f9730000dge0400000","411796105f032230000","0h63103100000010"); // EG1+DG1
  } else if (trainConfig == 2021){  // EMCal+DCAL clusters standard cuts, triggers, Nico TB NL + FineTuning 57, +-50ns timing (5)
    cuts.AddCutPCMCalo("80010123","00200009f9730000dge0400000","411796205f032230000","0h63103100000010"); // INT7
    cuts.AddCutPCMCalo("8008e123","00200009f9730000dge0400000","411796205f032230000","0h63103100000010"); // EG2+DG2
    cuts.AddCutPCMCalo("8008d123","00200009f9730000dge0400000","411796205f032230000","0h63103100000010"); // EG1+DG1
  } else if (trainConfig == 2022){  // EMCal+DCAL clusters standard cuts, triggers, Nico TB NL + FineTuning 57, +-50ns timing (5)
    cuts.AddCutPCMCalo("80010123","00200009f9730000dge0400000","411796305f032230000","0h63103100000010"); // INT7
    cuts.AddCutPCMCalo("8008e123","00200009f9730000dge0400000","411796305f032230000","0h63103100000010"); // EG2+DG2
    cuts.AddCutPCMCalo("8008d123","00200009f9730000dge0400000","411796305f032230000","0h63103100000010"); // EG1+DG1
  } else if (trainConfig == 2023){  // EMCal+DCAL clusters standard cuts, triggers, Nico TB NL + FineTuning 57, +-50ns timing (5)
    cuts.AddCutPCMCalo("80010123","00200009f9730000dge0400000","411796405f032230000","0h63103100000010"); // INT7
    cuts.AddCutPCMCalo("8008e123","00200009f9730000dge0400000","411796405f032230000","0h63103100000010"); // EG2+DG2
    cuts.AddCutPCMCalo("8008d123","00200009f9730000dge0400000","411796405f032230000","0h63103100000010"); // EG1+DG1
  } else if (trainConfig == 2024){  // EMCal+DCAL clusters standard cuts, triggers, Nico TB NL + FineTuning 57, +-50ns timing (5)
    cuts.AddCutPCMCalo("80010123","00200009f9730000dge0400000","411796505f032230000","0h63103100000010"); // INT7
    cuts.AddCutPCMCalo("8008e123","00200009f9730000dge0400000","411796505f032230000","0h63103100000010"); // EG2+DG2
    cuts.AddCutPCMCalo("8008d123","00200009f9730000dge0400000","411796505f032230000","0h63103100000010"); // EG1+DG1
  } else if (trainConfig == 2025){  // EMCal+DCAL clusters standard cuts, triggers, Nico TB NL + FineTuning 57, +-50ns timing (5)
    cuts.AddCutPCMCalo("80010123","00200009f9730000dge0400000","411790005f032230000","0h63103100000010"); // INT7
    cuts.AddCutPCMCalo("8008e123","00200009f9730000dge0400000","411790005f032230000","0h63103100000010"); // EG2+DG2
    cuts.AddCutPCMCalo("8008d123","00200009f9730000dge0400000","411790005f032230000","0h63103100000010"); // EG1+DG1
  } else if (trainConfig == 2026){  // EMCal+DCAL new PCM cut
    cuts.AddCutPCMCalo("80010123","00200009f9730000dge0400000","411790105f032230000","0h63103100000010"); // INT7
    cuts.AddCutPCMCalo("8008e123","00200009f9730000dge0400000","411790105f032230000","0h63103100000010"); // EG2+DG2
    cuts.AddCutPCMCalo("8008d123","00200009f9730000dge0400000","411790105f032230000","0h63103100000010"); // EG1+DG1

  // standard cut configs with no NL for SM-wise correction
  } else if (trainConfig == 2100){  // EMCal+DCAL clusters standard cuts, triggers, no NL, +-50ns timing (5)
    cuts.AddCutPCMCalo("80010123","00200009f9730000dge0400000","4117900057032230000","0h63103100000010"); // INT7
    cuts.AddCutPCMCalo("8008e123","00200009f9730000dge0400000","4117900057032230000","0h63103100000010"); // EG2+DG2
    cuts.AddCutPCMCalo("8008d123","00200009f9730000dge0400000","4117900057032230000","0h63103100000010"); // EG1+DG1
  } else if (trainConfig == 2101){  // EMCal+DCAL clusters standard cuts, triggers, no NL, +-50ns timing (5)
    cuts.AddCutPCMCalo("80010123","00200009f9730000dge0400000","4117900057032230000","0h63103100000010"); // INT7
    cuts.AddCutPCMCalo("8008e123","00200009f9730000dge0400000","4117900057032230000","0h63103100000010"); // EG2+DG2
    cuts.AddCutPCMCalo("8008d123","00200009f9730000dge0400000","4117900057032230000","0h63103100000010"); // EG1+DG1
  } else if (trainConfig == 2102){  // EMCal+DCAL clusters standard cuts, triggers, no NL, +-50ns timing (5)
    cuts.AddCutPCMCalo("80010123","00200009f9730000dge0400000","4117900057032230000","0h63103100000010"); // INT7
    cuts.AddCutPCMCalo("8008e123","00200009f9730000dge0400000","4117900057032230000","0h63103100000010"); // EG2+DG2
    cuts.AddCutPCMCalo("8008d123","00200009f9730000dge0400000","4117900057032230000","0h63103100000010"); // EG1+DG1
  } else if (trainConfig == 2103){  // EMCal+DCAL clusters standard cuts, triggers, no NL, +-50ns timing (5)
    cuts.AddCutPCMCalo("80010123","00200009f9730000dge0400000","4117900057032230000","0h63103100000010"); // INT7
    cuts.AddCutPCMCalo("8008e123","00200009f9730000dge0400000","4117900057032230000","0h63103100000010"); // EG2+DG2
    cuts.AddCutPCMCalo("8008d123","00200009f9730000dge0400000","4117900057032230000","0h63103100000010"); // EG1+DG1
  } else if (trainConfig == 2104){  // EMCal+DCAL clusters standard cuts, triggers, no NL, +-50ns timing (5)
    cuts.AddCutPCMCalo("80010123","00200009f9730000dge0400000","4117900057032230000","0h63103100000010"); // INT7
    cuts.AddCutPCMCalo("8008e123","00200009f9730000dge0400000","4117900057032230000","0h63103100000010"); // EG2+DG2
    cuts.AddCutPCMCalo("8008d123","00200009f9730000dge0400000","4117900057032230000","0h63103100000010"); // EG1+DG1
  } else if (trainConfig == 2105){  // EMCal+DCAL clusters standard cuts, triggers, no NL, +-50ns timing (5)
    cuts.AddCutPCMCalo("80010123","00200009f9730000dge0400000","4117900057032230000","0h63103100000010"); // INT7
    cuts.AddCutPCMCalo("8008e123","00200009f9730000dge0400000","4117900057032230000","0h63103100000010"); // EG2+DG2
    cuts.AddCutPCMCalo("8008d123","00200009f9730000dge0400000","4117900057032230000","0h63103100000010"); // EG1+DG1
  } else if (trainConfig == 2106){  // EMCal+DCAL clusters standard cuts, triggers, no NL, +-50ns timing (5)
    cuts.AddCutPCMCalo("80010123","00200009f9730000dge0400000","4117900057032230000","0h63103100000010"); // INT7
    cuts.AddCutPCMCalo("8008e123","00200009f9730000dge0400000","4117900057032230000","0h63103100000010"); // EG2+DG2
    cuts.AddCutPCMCalo("8008d123","00200009f9730000dge0400000","4117900057032230000","0h63103100000010"); // EG1+DG1
  } else if (trainConfig == 2107){  // EMCal+DCAL clusters standard cuts, triggers, no NL, +-50ns timing (5)
    cuts.AddCutPCMCalo("80010123","00200009f9730000dge0400000","4117900057032230000","0h63103100000010"); // INT7
    cuts.AddCutPCMCalo("8008e123","00200009f9730000dge0400000","4117900057032230000","0h63103100000010"); // EG2+DG2
    cuts.AddCutPCMCalo("8008d123","00200009f9730000dge0400000","4117900057032230000","0h63103100000010"); // EG1+DG1
  } else if (trainConfig == 2108){  // EMCal+DCAL clusters standard cuts, triggers, no NL, +-50ns timing (5)
    cuts.AddCutPCMCalo("80010123","00200009f9730000dge0400000","4117900057032230000","0h63103100000010"); // INT7
    cuts.AddCutPCMCalo("8008e123","00200009f9730000dge0400000","4117900057032230000","0h63103100000010"); // EG2+DG2
    cuts.AddCutPCMCalo("8008d123","00200009f9730000dge0400000","4117900057032230000","0h63103100000010"); // EG1+DG1
  } else if (trainConfig == 2109){  // EMCal+DCAL clusters standard cuts, triggers, no NL, +-50ns timing (5)
    cuts.AddCutPCMCalo("80010123","00200009f9730000dge0400000","4117900057032230000","0h63103100000010"); // INT7
    cuts.AddCutPCMCalo("8008e123","00200009f9730000dge0400000","4117900057032230000","0h63103100000010"); // EG2+DG2
    cuts.AddCutPCMCalo("8008d123","00200009f9730000dge0400000","4117900057032230000","0h63103100000010"); // EG1+DG1
  } else if (trainConfig == 2110){  // EMCal+DCAL clusters standard cuts, triggers, no NL, +-50ns timing (5)
    cuts.AddCutPCMCalo("80010123","00200009f9730000dge0400000","4117900057032230000","0h63103100000010"); // INT7
    cuts.AddCutPCMCalo("8008e123","00200009f9730000dge0400000","4117900057032230000","0h63103100000010"); // EG2+DG2
    cuts.AddCutPCMCalo("8008d123","00200009f9730000dge0400000","4117900057032230000","0h63103100000010"); // EG1+DG1
  } else if (trainConfig == 2111){  // EMCal+DCAL clusters standard cuts, triggers, no NL, +-50ns timing (5)
    cuts.AddCutPCMCalo("80010123","00200009f9730000dge0400000","4117900057032230000","0h63103100000010"); // INT7
    cuts.AddCutPCMCalo("8008e123","00200009f9730000dge0400000","4117900057032230000","0h63103100000010"); // EG2+DG2
    cuts.AddCutPCMCalo("8008d123","00200009f9730000dge0400000","4117900057032230000","0h63103100000010"); // EG1+DG1
  } else if (trainConfig == 2112){  // EMCal+DCAL clusters standard cuts, triggers, no NL, +-50ns timing (5)
    cuts.AddCutPCMCalo("80010123","00200009f9730000dge0400000","4117900057032230000","0h63103100000010"); // INT7
    cuts.AddCutPCMCalo("8008e123","00200009f9730000dge0400000","4117900057032230000","0h63103100000010"); // EG2+DG2
    cuts.AddCutPCMCalo("8008d123","00200009f9730000dge0400000","4117900057032230000","0h63103100000010"); // EG1+DG1
  } else if (trainConfig == 2113){  // EMCal+DCAL clusters standard cuts, triggers, no NL, +-50ns timing (5)
    cuts.AddCutPCMCalo("80010123","00200009f9730000dge0400000","4117900057032230000","0h63103100000010"); // INT7
    cuts.AddCutPCMCalo("8008e123","00200009f9730000dge0400000","4117900057032230000","0h63103100000010"); // EG2+DG2
    cuts.AddCutPCMCalo("8008d123","00200009f9730000dge0400000","4117900057032230000","0h63103100000010"); // EG1+DG1
  } else if (trainConfig == 2114){  // EMCal+DCAL clusters standard cuts, triggers, no NL, +-50ns timing (5)
    cuts.AddCutPCMCalo("80010123","00200009f9730000dge0400000","4117900057032230000","0h63103100000010"); // INT7
    cuts.AddCutPCMCalo("8008e123","00200009f9730000dge0400000","4117900057032230000","0h63103100000010"); // EG2+DG2
    cuts.AddCutPCMCalo("8008d123","00200009f9730000dge0400000","4117900057032230000","0h63103100000010"); // EG1+DG1
  } else if (trainConfig == 2115){  // EMCal+DCAL clusters standard cuts, triggers, no NL, +-50ns timing (5)
    cuts.AddCutPCMCalo("80010123","00200009f9730000dge0400000","4117900057032230000","0h63103100000010"); // INT7
    cuts.AddCutPCMCalo("8008e123","00200009f9730000dge0400000","4117900057032230000","0h63103100000010"); // EG2+DG2
    cuts.AddCutPCMCalo("8008d123","00200009f9730000dge0400000","4117900057032230000","0h63103100000010"); // EG1+DG1
  } else if (trainConfig == 2116){  // EMCal+DCAL clusters standard cuts, triggers, no NL, +-50ns timing (5)
    cuts.AddCutPCMCalo("80010123","00200009f9730000dge0400000","4117900057032230000","0h63103100000010"); // INT7
    cuts.AddCutPCMCalo("8008e123","00200009f9730000dge0400000","4117900057032230000","0h63103100000010"); // EG2+DG2
    cuts.AddCutPCMCalo("8008d123","00200009f9730000dge0400000","4117900057032230000","0h63103100000010"); // EG1+DG1
  } else if (trainConfig == 2117){  // EMCal+DCAL clusters standard cuts, triggers, no NL, +-50ns timing (5)
    cuts.AddCutPCMCalo("80010123","00200009f9730000dge0400000","4117900057032230000","0h63103100000010"); // INT7
    cuts.AddCutPCMCalo("8008e123","00200009f9730000dge0400000","4117900057032230000","0h63103100000010"); // EG2+DG2
    cuts.AddCutPCMCalo("8008d123","00200009f9730000dge0400000","4117900057032230000","0h63103100000010"); // EG1+DG1
  } else if (trainConfig == 2118){  // EMCal+DCAL clusters standard cuts, triggers, no NL, +-50ns timing (5)
    cuts.AddCutPCMCalo("80010123","00200009f9730000dge0400000","4117900057032230000","0h63103100000010"); // INT7
    cuts.AddCutPCMCalo("8008e123","00200009f9730000dge0400000","4117900057032230000","0h63103100000010"); // EG2+DG2
    cuts.AddCutPCMCalo("8008d123","00200009f9730000dge0400000","4117900057032230000","0h63103100000010"); // EG1+DG1
  } else if (trainConfig == 2119){  // EMCal+DCAL clusters standard cuts, triggers, no NL, +-50ns timing (5)
    cuts.AddCutPCMCalo("80010123","00200009f9730000dge0400000","4117900057032230000","0h63103100000010"); // INT7
    cuts.AddCutPCMCalo("8008e123","00200009f9730000dge0400000","4117900057032230000","0h63103100000010"); // EG2+DG2
    cuts.AddCutPCMCalo("8008d123","00200009f9730000dge0400000","4117900057032230000","0h63103100000010"); // EG1+DG1
  } else if (trainConfig == 2120){  // EMCal+DCAL clusters standard cuts, triggers, no NL, +-50ns timing (5)
    cuts.AddCutPCMCalo("80010123","00200009f9730000dge0400000","4117900057032230000","0h63103100000010"); // INT7
    cuts.AddCutPCMCalo("8008e123","00200009f9730000dge0400000","4117900057032230000","0h63103100000010"); // EG2+DG2
    cuts.AddCutPCMCalo("8008d123","00200009f9730000dge0400000","4117900057032230000","0h63103100000010"); // EG1+DG1
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
    analysisEventCuts[i]->SetTriggerOverlapRejecion(enableTriggerOverlapRej);
    if(fMinPtHardSet)
      analysisEventCuts[i]->SetMinFacPtHard(minFacPtHard);
    if(fMaxPtHardSet)
      analysisEventCuts[i]->SetMaxFacPtHard(maxFacPtHard);
    if(fSingleMaxPtHardSet)
      analysisEventCuts[i]->SetMaxFacPtHardSingleParticle(maxFacPtHardSingle);    analysisEventCuts[i]->SetV0ReaderName(V0ReaderName);
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
  if(enableQAPhotonTask>1){
    for(Int_t i = 0; i<numberOfCuts; i++){
      mgr->ConnectOutput(task,2+i,mgr->CreateContainer(Form("%s_%s_%s_%s Photon DCA tree",(cuts.GetEventCut(i)).Data(),(cuts.GetPhotonCut(i)).Data(),cuts.GetClusterCut(i).Data(),(cuts.GetMesonCut(i)).Data()), TTree::Class(), AliAnalysisManager::kOutputContainer, Form("GammaConvCalo_%i.root",trainConfig)) );
    }
  }

  return;

}
