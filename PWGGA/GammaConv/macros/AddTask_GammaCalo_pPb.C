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
//($ALIPHYSICS/PWGGA/GammaConv/AliAnalysisTaskGammaCalo.cxx) for
//pPb together with all supporting classes
//***************************************************************************************

//***************************************************************************************
//main function
//***************************************************************************************
void AddTask_GammaCalo_pPb(
  Int_t     trainConfig                   = 1,        // change different set of cuts
  Int_t     isMC                          = 0,        // run MC
  TString   photonCutNumberV0Reader       = "",       // 00000008400000000100000000 nom. B, 00000088400000000100000000 low B
  TString   periodNameV0Reader            = "",
  // general setting for task
  Int_t     enableQAMesonTask             = 0,        // enable QA in AliAnalysisTaskGammaConvV1
  Int_t     enableQAClusterTask           = 0,        // enable additional QA task
  Int_t     enableExtMatchAndQA           = 0,        // disabled (0), extMatch (1), extQA_noCellQA (2), extMatch+extQA_noCellQA (3), extQA+cellQA (4), extMatch+extQA+cellQA (5)
  Bool_t    enableLightOutput             = kFALSE,   // switch to run light output (only essential histograms for afterburner)
  Bool_t    enableTHnSparse               = kFALSE,   // switch on THNsparse
  Bool_t    enableTriggerMimicking        = kFALSE,   // enable trigger mimicking
  Bool_t    enableTriggerOverlapRej       = kFALSE,   // enable trigger overlap rejection
  TString   settingMaxFacPtHard           = "3.",     // maximum factor between hardest jet and ptHard generated
  Int_t     debugLevel                    = 0,        // introducing debug levels for grid running
  // settings for weights
  // FPTW:fileNamePtWeights, FMUW:fileNameMultWeights, separate with ;
  TString   fileNameExternalInputs        = "",
  Int_t     doWeightingPart               = 0,        // enable Weighting
  TString   generatorName                 = "DPMJET", // generator Name
  Bool_t    enableMultiplicityWeighting   = kFALSE,   //
  TString   periodNameAnchor              = "",       //
  // special settings
  Bool_t    enableSortingMCLabels         = kTRUE,    // enable sorting for MC cluster labels
  // subwagon config
  TString   additionalTrainConfig         = "0"       // additional counter for trainconfig
                           ) {

  AliCutHandlerPCM cuts;

  TString fileNamePtWeights     = cuts.GetSpecialFileNameFromString (fileNameExternalInputs, "FPTW:");
  TString fileNameMultWeights   = cuts.GetSpecialFileNameFromString (fileNameExternalInputs, "FMUW:");

  TString addTaskName                 = "AddTask_GammaCalo_pPb";
  TString sAdditionalTrainConfig      = cuts.GetSpecialSettingFromAddConfig(additionalTrainConfig, "", "", addTaskName);
  if (sAdditionalTrainConfig.Atoi() > 0){
    trainConfig = trainConfig + sAdditionalTrainConfig.Atoi();
    cout << "INFO: " << addTaskName.Data() << " running additionalTrainConfig '" << sAdditionalTrainConfig.Atoi() << "', train config: '" << trainConfig << "'" << endl;
  }

  TString corrTaskSetting             = cuts.GetSpecialSettingFromAddConfig(additionalTrainConfig, "CF", "", addTaskName);
  if(corrTaskSetting.CompareTo(""))
    cout << "corrTaskSetting: " << corrTaskSetting.Data() << endl;

  Int_t trackMatcherRunningMode = 0; // CaloTrackMatcher running mode
  TString strTrackMatcherRunningMode  = cuts.GetSpecialSettingFromAddConfig(additionalTrainConfig, "TM", "", addTaskName);
  if(additionalTrainConfig.Contains("TM"))
    trackMatcherRunningMode = strTrackMatcherRunningMode.Atoi();

  Bool_t doTreeEOverP = kFALSE; // switch to produce EOverP tree
  TString strdoTreeEOverP             = cuts.GetSpecialSettingFromAddConfig(additionalTrainConfig, "EPCLUSTree", "", addTaskName);
  if(strdoTreeEOverP.Atoi()==1)
    doTreeEOverP = kTRUE;

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
  if(rmaxFacPtHardSetting->GetEntries()<1){cout << "ERROR: AddTask_GammaCalo_pPb during parsing of settingMaxFacPtHard String '" << settingMaxFacPtHard.Data() << "'" << endl; return;}
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
  AliAnalysisTaskGammaCalo *task=NULL;
  task= new AliAnalysisTaskGammaCalo(Form("GammaCalo_%i",trainConfig));
  task->SetIsHeavyIon(isHeavyIon);
  task->SetIsMC(isMC);
  task->SetV0ReaderName(V0ReaderName);
  task->SetCorrectionTaskSetting(corrTaskSetting);
  task->SetLightOutput(enableLightOutput);
  task->SetTrackMatcherRunningMode(trackMatcherRunningMode);

  // cluster cuts
  // 0 "ClusterType",  1 "EtaMin", 2 "EtaMax", 3 "PhiMin", 4 "PhiMax", 5 "DistanceToBadChannel", 6 "Timing", 7 "TrackMatching", 8 "ExoticCell",
  // 9 "MinEnergy", 10 "MinNCells", 11 "MinM02", 12 "MaxM02", 13 "MinMaxM20", 14 "RecConv", 15 "MaximumDispersion", 16 "NLM"

  // ===============================================================================================
  // Run 1 data EMC clusters pPb 5TeV
  // ===============================================================================================
  if (trainConfig == 1){ // no non lin
    cuts.AddCutCalo("80010113","1111100057032230000","01631031000000d0");
    cuts.AddCutCalo("80052113","1111100057032230000","01631031000000d0");
    cuts.AddCutCalo("80085113","1111100057032230000","01631031000000d0");
    cuts.AddCutCalo("80083113","1111100057032230000","01631031000000d0");
  } else if (trainConfig == 2){ // no non lin
    cuts.AddCutCalo("80010113","1111100057032230000","01631031000000d0");
  } else if (trainConfig == 3){ // no non lin
    cuts.AddCutCalo("80210113","1111100057032230000","01631031000000d0"); // 0-20
    cuts.AddCutCalo("82410113","1111100057032230000","01631031000000d0"); // 20-40
    cuts.AddCutCalo("84610113","1111100057032230000","01631031000000d0"); // 40-60
    cuts.AddCutCalo("86010113","1111100057032230000","01631031000000d0"); // 60-80
  } else if (trainConfig == 4){ // no non lin, no time cut
    cuts.AddCutCalo("80010113","1111100007032230000","01631031000000d0");
    cuts.AddCutCalo("80052113","1111100007032230000","01631031000000d0");
    cuts.AddCutCalo("80085113","1111100007032230000","01631031000000d0");
    cuts.AddCutCalo("80083113","1111100007032230000","01631031000000d0");
  } else if (trainConfig == 5){ // no non lin, no time
    cuts.AddCutCalo("80010113","1111100007032230000","01631031000000d0");
  } else if (trainConfig == 6){ // no non lin, no time cut, cent dep
    cuts.AddCutCalo("80210113","1111100007032230000","01631031000000d0"); // 0-20
    cuts.AddCutCalo("82410113","1111100007032230000","01631031000000d0"); // 20-40
    cuts.AddCutCalo("84610113","1111100007032230000","01631031000000d0"); // 40-60
    cuts.AddCutCalo("86010113","1111100007032230000","01631031000000d0"); // 60-80

  } else if (trainConfig == 7){
    cuts.AddCutCalo("90010113","1111141057032230000","01631031000000d0"); // 0-100
    cuts.AddCutCalo("90210113","1111141057032230000","01631031000000d0"); // 0-20
    cuts.AddCutCalo("92410113","1111141057032230000","01631031000000d0"); // 20-40
    cuts.AddCutCalo("94610113","1111141057032230000","01631031000000d0"); // 40-60
    cuts.AddCutCalo("96010113","1111141057032230000","01631031000000d0"); // 60-100
  } else if (trainConfig == 8){
    cuts.AddCutCalo("e0010113","1111141057032230000","01631031000000d0"); // 0-100
    cuts.AddCutCalo("e0210113","1111141057032230000","01631031000000d0"); // 0-20
    cuts.AddCutCalo("e2410113","1111141057032230000","01631031000000d0"); // 20-40
    cuts.AddCutCalo("e4610113","1111141057032230000","01631031000000d0"); // 40-60
    cuts.AddCutCalo("e6010113","1111141057032230000","01631031000000d0"); // 60-100

  } else if (trainConfig == 9){ // new TB Nico
    cuts.AddCutCalo("80010113","1111100057032230000","01631031000000d0"); // 0-100
    cuts.AddCutCalo("80010113","1111141057032230000","01631031000000d0"); // 0-100
    cuts.AddCutCalo("80010113","1111142057032230000","01631031000000d0"); // 0-100
    cuts.AddCutCalo("80010113","1111165057032230000","01631031000000d0"); // 0-100
    //-----------------------------------------------------------------------------------------------
  // Standard cuts
  //-----------------------------------------------------------------------------------------------
  } else if (trainConfig == 10){ // default cutstring smaller openangle 1 cell
    cuts.AddCutCalo("80010113","1111141057032230000","01631031000000d0"); // default with 1 cell
    cuts.AddCutCalo("80010113","1111141057032230000","0163103100000060"); // default with 17mrad
  } else if (trainConfig == 11){ // new default cut
    cuts.AddCutCalo("80010113","1111141057032230000","01631031000000d0"); // default with 1cell
    cuts.AddCutCalo("80010113","1111141057032230000","0i631031000000d0"); // default with 1cell sector mixing
  } else if (trainConfig == 12){ //all default triggers
    cuts.AddCutCalo("80010113","1111141057032230000","01631031000000d0"); // default MB
    cuts.AddCutCalo("80052113","1111141057032230000","01631031000000d0"); // default EMC7
    cuts.AddCutCalo("80083113","1111141057032230000","01631031000000d0"); // default EG1
    cuts.AddCutCalo("80085113","1111141057032230000","01631031000000d0"); // default EG2
  } else if (trainConfig == 13){
    cuts.AddCutCalo("80010113","11111410570322l0000","01631031000000d0"); // default MB
    cuts.AddCutCalo("80052113","11111410570322l0000","01631031000000d0"); // default EMC7
    cuts.AddCutCalo("80083113","11111410570322l0000","01631031000000d0"); // default EG1
    cuts.AddCutCalo("80085113","11111410570322l0000","01631031000000d0"); // default EG2
  //-----------------------------------------------------------------------------------------------
  // Standard cuts cent dependent
  //-----------------------------------------------------------------------------------------------
  } else if (trainConfig == 14){
    cuts.AddCutCalo("80010113","1111141057032230000","01631031000000d0"); // 0-100
    cuts.AddCutCalo("80210113","1111141057032230000","01631031000000d0"); // 0-20
    cuts.AddCutCalo("82410113","1111141057032230000","01631031000000d0"); // 20-40
    cuts.AddCutCalo("84610113","1111141057032230000","01631031000000d0"); // 40-60
    cuts.AddCutCalo("86010113","1111141057032230000","01631031000000d0"); // 60-80
  } else if (trainConfig == 15){ // direct photons sector mixing
    cuts.AddCutCalo("80010113","1111141057032230000","0i631031000000d0"); // 0-100
    cuts.AddCutCalo("80210113","1111141057032230000","0i631031000000d0"); // 0-20
    cuts.AddCutCalo("82410113","1111141057032230000","0i631031000000d0"); // 20-40
    cuts.AddCutCalo("84610113","1111141057032230000","0i631031000000d0"); // 40-60
    cuts.AddCutCalo("86010113","1111141057032230000","0i631031000000d0"); // 60-80
  } else if (trainConfig == 17){ // direct photons sector mixing
    cuts.AddCutCalo("80010113","11111410570322l0000","0i631031000000d0");
    cuts.AddCutCalo("80210113","11111410570322l0000","0i631031000000d0");
    cuts.AddCutCalo("82410113","11111410570322l0000","0i631031000000d0");
    cuts.AddCutCalo("84610113","11111410570322l0000","0i631031000000d0");
    cuts.AddCutCalo("86010113","11111410570322l0000","0i631031000000d0");
  } else if (trainConfig == 18){ // direct photons
    cuts.AddCutCalo("80010113","11111410570322l0000","01631031000000d0");
    cuts.AddCutCalo("80210113","11111410570322l0000","01631031000000d0");
    cuts.AddCutCalo("82410113","11111410570322l0000","01631031000000d0");
    cuts.AddCutCalo("84610113","11111410570322l0000","01631031000000d0");
    cuts.AddCutCalo("86010113","11111410570322l0000","01631031000000d0");
  } else if (trainConfig == 19){ // direct photons CL1
    cuts.AddCutCalo("90010113","11111410570322l0000","01631031000000d0");
    cuts.AddCutCalo("90210113","11111410570322l0000","01631031000000d0");
    cuts.AddCutCalo("92410113","11111410570322l0000","01631031000000d0");
    cuts.AddCutCalo("94610113","11111410570322l0000","01631031000000d0");
    cuts.AddCutCalo("96010113","11111410570322l0000","01631031000000d0");

  } else if (trainConfig == 21){ // default cutstring, 1cell distance lead cell
    cuts.AddCutCalo("80010113","1111141057032230000","01631031000000a0"); // 1 cell lead cell
    cuts.AddCutCalo("80010113","1111141057032230000","01631031000000d0"); // 1 cell lead cell, 17mrad open
    cuts.AddCutCalo("80010113","1111141057032230000","01631031000000b0"); // 1 cell lead cell, 15mrad open
  } else if (trainConfig == 22){ // default cutstring, M02 variations
    cuts.AddCutCalo("80010113","1111141057032250000","01631031000000d0"); // 0.3
    cuts.AddCutCalo("80010113","1111141057032260000","01631031000000d0"); // 0.27
    cuts.AddCutCalo("80010113","1111141057032240000","01631031000000d0"); // 0.4
  } else if (trainConfig == 23){ // default cutstring, M02 variations
    cuts.AddCutCalo("80010113","1111141057032290000","01631031000000d0"); // 0.35
    cuts.AddCutCalo("80010113","11111410570322a0000","01631031000000d0"); // 0.33
    cuts.AddCutCalo("80010113","11111410570322b0000","01631031000000d0"); // 0.28
    cuts.AddCutCalo("80010113","11111410570322c0000","01631031000000d0"); // 0.32
  } else if (trainConfig == 24){ // default cutstring, cluster energy variations, decreased tender thresholds
    cuts.AddCutCalo("80010113","1111141057012230000","01631031000000d0"); // E cluster > 0.5
    cuts.AddCutCalo("80010113","1111141057032230000","01631031000000d0"); // E cluster > 0.7
    cuts.AddCutCalo("80010113","1111141057052230000","01631031000000d0"); // E cluster > 0.9
    cuts.AddCutCalo("80010113","11111410570b2230000","01631031000000d0"); // E cluster > 1.0
    cuts.AddCutCalo("80010113","11111410570a2230000","01631031000000d0"); // E cluster > 1.5
  } else if (trainConfig == 25){ // default cutstring, cluster energy variations, same tender thresholds
    cuts.AddCutCalo("80010113","1111141057032230000","01631031000000d0"); // E cluster > 0.7
    cuts.AddCutCalo("80010113","1111141057052230000","01631031000000d0"); // E cluster > 0.9
    cuts.AddCutCalo("80010113","11111410570b2230000","01631031000000d0"); // E cluster > 1.0
    cuts.AddCutCalo("80010113","11111410570a2230000","01631031000000d0"); // E cluster > 1.5
  } else if (trainConfig == 26){ // default cutstring, cluster energy variations, increased tender thresholds
    cuts.AddCutCalo("80010113","1111141057052230000","01631031000000d0"); // E cluster > 0.9
    cuts.AddCutCalo("80010113","11111410570b2230000","01631031000000d0"); // E cluster > 1.0
    cuts.AddCutCalo("80010113","11111410570a2230000","01631031000000d0"); // E cluster > 1.5
  } else if (trainConfig == 27){ // direct photons ZNA
    cuts.AddCutCalo("e0010113","11111410570322l0000","01631031000000d0");
    cuts.AddCutCalo("e0210113","11111410570322l0000","01631031000000d0");
    cuts.AddCutCalo("e2410113","11111410570322l0000","01631031000000d0");
    cuts.AddCutCalo("e4610113","11111410570322l0000","01631031000000d0");
    cuts.AddCutCalo("e6010113","11111410570322l0000","01631031000000d0");
    //-----------------------------------------------------------------------------------------------
  // Systematics variations MB
  //-----------------------------------------------------------------------------------------------
  } else if (trainConfig == 30){ // nonlinearity variations
    cuts.AddCutCalo("80010113","1111142057032230000","01631031000000d0"); // CRF
    cuts.AddCutCalo("80010113","1111151057032230000","01631031000000d0"); // CCMF
    cuts.AddCutCalo("80010113","1111152057032230000","01631031000000d0"); // CMF
  } else if (trainConfig == 31){ // second set of variations CLUSTER
    cuts.AddCutCalo("80010113","1111141057022230000","01631031000000d0"); // min energy cluster variation 1  600 MeV
    cuts.AddCutCalo("80010113","1111141057042230000","01631031000000d0"); // min energy cluster variation 2  800 MeV
    cuts.AddCutCalo("80010113","1111141057052230000","01631031000000d0"); // min energy cluster variation 3  900 MeV
    cuts.AddCutCalo("80010113","1111141057032230000","01631031000000d0"); // min/max M02  0.1<M<0.5
    cuts.AddCutCalo("80010113","1111141057032200000","01631031000000d0"); // min/max M02  0.1<M<100
  } else if (trainConfig == 32){ // third set of variations CLUSTER
    cuts.AddCutCalo("80010113","1111141057032250000","01631031000000d0"); // min/max M02  0.1<M<0.3
    cuts.AddCutCalo("80010113","1111141057032260000","01631031000000d0"); // min/max M02  0.1<M<0.27
    cuts.AddCutCalo("80010113","1111141057031230000","01631031000000d0"); // min number of cells variation 1  1 cell
    cuts.AddCutCalo("80010113","1112141057032230000","01631031000000d0"); // only modules with TRD infront
    cuts.AddCutCalo("80010113","1111341057032230000","01631031000000d0"); // no modules with TRD infront
  } else if (trainConfig == 33){ // third set of variations MESON
    cuts.AddCutCalo("80010113","1111141057032230000","01633031000000d0"); // rapidity variation  y<0.6
    cuts.AddCutCalo("80010113","1111141057032230000","01634031000000d0"); // rapidity variation  y<0.5
    cuts.AddCutCalo("80010113","1111141057032230000","01631061000000d0"); // alpha meson variation 1   0<alpha<0.8
    cuts.AddCutCalo("80010113","1111141057032230000","01631051000000d0"); // alpha meson variation 2  0<alpha<0.75
  } else if (trainConfig == 34){ // opening angle variations
    cuts.AddCutCalo("80010113","1111141057032230000","0163103100000040"); // min opening angle 0.0152
    cuts.AddCutCalo("80010113","1111141057032230000","0163103100000060"); // min opening angle 0.017
    cuts.AddCutCalo("80010113","1111141057032230000","0163103100000070"); // min opening angle 0.016
    cuts.AddCutCalo("80010113","1111141057032230000","0163103100000080"); // min opening angle 0.018
    cuts.AddCutCalo("80010113","1111141057032230000","0163103100000090"); // min opening angle 0.018
  } else if (trainConfig == 35){ // TM variations
    cuts.AddCutCalo("80010113","1111141053032230000","01631031000000d0"); // fixed window
    cuts.AddCutCalo("80010113","1111141056032230000","01631031000000d0"); // tm pt dependent var 1
    cuts.AddCutCalo("80010113","1111141058032230000","01631031000000d0"); // tm pt dependent var 2
    cuts.AddCutCalo("80010113","1111141059032230000","01631031000000d0"); // tm pt dependent var 3
  } else if (trainConfig == 36){ // TM variations
    cuts.AddCutCalo("80010113","1111153057032230000","01631031000000d0"); // CCMF nonlin+testbeam
    cuts.AddCutCalo("80010113","1111154057032230000","01631031000000d0"); // CMF nonlin+testbeam
    cuts.AddCutCalo("80010113","1111143057032230000","01631031000000d0"); // CCRF testbeam nonlin
    cuts.AddCutCalo("80010113","1111144057032230000","01631031000000d0"); // CRF testbeam nonlin
    cuts.AddCutCalo("80010113","1111102057032230000","01631031000000d0"); // testbeam nonlin

  } else if (trainConfig == 37){
    cuts.AddCutCalo("80010113","11111410570322l0000","01631031000000d0"); // M02 pt dep with new std: 0.32, 0.0072, 0.5
    cuts.AddCutCalo("80010113","11111410570322d0000","01631031000000d0"); // M02, pt dep with  0.27, 0.0072, 0.4
    cuts.AddCutCalo("80010113","11111410570322e0000","01631031000000d0"); // M02, pt dep with  0.31, 0.0072, 0.5
    cuts.AddCutCalo("80010113","11111410570322f0000","01631031000000d0"); // M02, pt dep with  0.36, 0.0072, 0.7
  } else if (trainConfig == 38){
    cuts.AddCutCalo("80010113","11111410570322m0000","01631031000000d0"); // M02, pt dep with  0.32, 0.0152, 0.5
    cuts.AddCutCalo("80010113","11111410570322g0000","01631031000000d0"); // M02, pt dep with  0.37, 0.0072, 0.7
    cuts.AddCutCalo("80010113","11111410570322h0000","01631031000000d0"); // M02, pT-dep with  0.30, 0.0072, 0.5
    cuts.AddCutCalo("80010113","11111410570322i0000","01631031000000d0"); // M02, pT-dep with  0.35, 0.0072, 0.7
  } else if (trainConfig == 39){
    cuts.AddCutCalo("80010113","11111410570322j0000","01631031000000d0"); // M02, pT-dep with  0.25, 0.0072, 0.39
    cuts.AddCutCalo("80010113","11111410570322r0000","01631031000000d0"); // M02, pT-dep with  0.25, 0.0072, 0.5
    cuts.AddCutCalo("80010113","11111410570322s0000","01631031000000d0"); // M02, pT-dep with  0.32, 0.0238, 0.7
    cuts.AddCutCalo("80010113","11111410570322n0000","01631031000000d0"); // M02, pT-dep with  0.32, 0.0238, 0.7

  //-----------------------------------------------------------------------------------------------
  // Systematics variations 0-20
  //-----------------------------------------------------------------------------------------------
  } else if (trainConfig == 40){ // nonlinearity variations
    cuts.AddCutCalo("80210113","1111142057032230000","01631031000000d0"); // CRF
    cuts.AddCutCalo("80210113","1111151057032230000","01631031000000d0"); // CCMF
    cuts.AddCutCalo("80210113","1111152057032230000","01631031000000d0"); // CMF
  } else if (trainConfig == 41){ // second set of variations CLUSTER
    cuts.AddCutCalo("80210113","1111141057022230000","01631031000000d0"); // min energy cluster variation 1  600 MeV
    cuts.AddCutCalo("80210113","1111141057042230000","01631031000000d0"); // min energy cluster variation 2  800 MeV
    cuts.AddCutCalo("80210113","1111141057052230000","01631031000000d0"); // min energy cluster variation 3  900 MeV
    cuts.AddCutCalo("80210113","1111141057032230000","01631031000000d0"); // min/max M02  0.1<M<0.5
    cuts.AddCutCalo("80210113","1111141057032200000","01631031000000d0"); // min/max M02  0.1<M<100
  } else if (trainConfig == 42){ // third set of variations CLUSTER
    cuts.AddCutCalo("80210113","1111141057032250000","01631031000000d0"); // min/max M02  0.1<M<0.3
    cuts.AddCutCalo("80210113","1111141057032260000","01631031000000d0"); // min/max M02  0.1<M<0.27
    cuts.AddCutCalo("80210113","1111141057031230000","01631031000000d0"); // min number of cells variation 1  1 cell
    cuts.AddCutCalo("80210113","1112141057032230000","01631031000000d0"); // only modules with TRD infront
    cuts.AddCutCalo("80210113","1111341057032230000","01631031000000d0"); // no modules with TRD infront
  } else if (trainConfig == 43){ // third set of variations MESON
    cuts.AddCutCalo("80210113","1111141057032230000","01633031000000d0"); // rapidity variation  y<0.6
    cuts.AddCutCalo("80210113","1111141057032230000","01634031000000d0"); // rapidity variation  y<0.5
    cuts.AddCutCalo("80210113","1111141057032230000","01631061000000d0"); // alpha meson variation 1   0<alpha<0.8
    cuts.AddCutCalo("80210113","1111141057032230000","01631051000000d0"); // alpha meson variation 2  0<alpha<0.75
  } else if (trainConfig == 44){ // opening angle variations
    cuts.AddCutCalo("80210113","1111141057032230000","0163103100000040"); // min opening angle 0.0152
    cuts.AddCutCalo("80210113","1111141057032230000","0163103100000060"); // min opening angle 0.017
    cuts.AddCutCalo("80210113","1111141057032230000","0163103100000070"); // min opening angle 0.016
    cuts.AddCutCalo("80210113","1111141057032230000","0163103100000080"); // min opening angle 0.018
    cuts.AddCutCalo("80210113","1111141057032230000","0163103100000090"); // min opening angle 0.018
  } else if (trainConfig == 45){ // TM variations
    cuts.AddCutCalo("80210113","1111141053032230000","01631031000000d0"); // fixed window
    cuts.AddCutCalo("80210113","1111141056032230000","01631031000000d0"); // tm pt dependent var 1
    cuts.AddCutCalo("80210113","1111141058032230000","01631031000000d0"); // tm pt dependent var 2
    cuts.AddCutCalo("80210113","1111141059032230000","01631031000000d0"); // tm pt dependent var 3
  } else if (trainConfig == 46){ // TM variations
    cuts.AddCutCalo("80210113","1111153057032230000","01631031000000d0"); // CCMF nonlin+testbeam
    cuts.AddCutCalo("80210113","1111154057032230000","01631031000000d0"); // CMF nonlin+testbeam
    cuts.AddCutCalo("80210113","1111143057032230000","01631031000000d0"); // CCRF testbeam nonlin
    cuts.AddCutCalo("80210113","1111144057032230000","01631031000000d0"); // CRF testbeam nonlin
    cuts.AddCutCalo("80210113","1111102057032230000","01631031000000d0"); // testbeam nonlin
  } else if (trainConfig == 47){ //all default triggers
    cuts.AddCutCalo("80210113","1111141057032230000","01631031000000d0"); // default MB
    cuts.AddCutCalo("80252113","1111141057032230000","01631031000000d0"); // default EMC7
    cuts.AddCutCalo("80283113","1111141057032230000","01631031000000d0"); // default EG1
    cuts.AddCutCalo("80285113","1111141057032230000","01631031000000d0"); // default EG2

  //-----------------------------------------------------------------------------------------------
  // Systematics variations 20-40
  //-----------------------------------------------------------------------------------------------
  } else if (trainConfig == 50){ // nonlinearity variations
    cuts.AddCutCalo("82410113","1111142057032230000","01631031000000d0"); // CRF
    cuts.AddCutCalo("82410113","1111151057032230000","01631031000000d0"); // CCMF
    cuts.AddCutCalo("82410113","1111152057032230000","01631031000000d0"); // CMF
  } else if (trainConfig == 51){ // second set of variations CLUSTER
    cuts.AddCutCalo("82410113","1111141057022230000","01631031000000d0"); // min energy cluster variation 1  600 MeV
    cuts.AddCutCalo("82410113","1111141057042230000","01631031000000d0"); // min energy cluster variation 2  800 MeV
    cuts.AddCutCalo("82410113","1111141057052230000","01631031000000d0"); // min energy cluster variation 3  900 MeV
    cuts.AddCutCalo("82410113","1111141057032230000","01631031000000d0"); // min/max M02  0.1<M<0.5
    cuts.AddCutCalo("82410113","1111141057032200000","01631031000000d0"); // min/max M02  0.1<M<100
  } else if (trainConfig == 52){ // third set of variations CLUSTER
    cuts.AddCutCalo("82410113","1111141057032250000","01631031000000d0"); // min/max M02  0.1<M<0.3
    cuts.AddCutCalo("82410113","1111141057032260000","01631031000000d0"); // min/max M02  0.1<M<0.27
    cuts.AddCutCalo("82410113","1111141057031230000","01631031000000d0"); // min number of cells variation 1  1 cell
    cuts.AddCutCalo("82410113","1112141057032230000","01631031000000d0"); // only modules with TRD infront
    cuts.AddCutCalo("82410113","1111341057032230000","01631031000000d0"); // no modules with TRD infront
  } else if (trainConfig == 53){ // third set of variations MESON
    cuts.AddCutCalo("82410113","1111141057032230000","01633031000000d0"); // rapidity variation  y<0.6
    cuts.AddCutCalo("82410113","1111141057032230000","01634031000000d0"); // rapidity variation  y<0.5
    cuts.AddCutCalo("82410113","1111141057032230000","01631061000000d0"); // alpha meson variation 1   0<alpha<0.8
    cuts.AddCutCalo("82410113","1111141057032230000","01631051000000d0"); // alpha meson variation 2  0<alpha<0.75
  } else if (trainConfig == 54){ // opening angle variations
    cuts.AddCutCalo("82410113","1111141057032230000","0163103100000040"); // min opening angle 0.0152
    cuts.AddCutCalo("82410113","1111141057032230000","0163103100000060"); // min opening angle 0.017
    cuts.AddCutCalo("82410113","1111141057032230000","0163103100000070"); // min opening angle 0.016
    cuts.AddCutCalo("82410113","1111141057032230000","0163103100000080"); // min opening angle 0.018
    cuts.AddCutCalo("82410113","1111141057032230000","0163103100000090"); // min opening angle 0.018
  } else if (trainConfig == 55){ // TM variations
    cuts.AddCutCalo("82410113","1111141053032230000","01631031000000d0"); // fixed window
    cuts.AddCutCalo("82410113","1111141056032230000","01631031000000d0"); // tm pt dependent var 1
    cuts.AddCutCalo("82410113","1111141058032230000","01631031000000d0"); // tm pt dependent var 2
    cuts.AddCutCalo("82410113","1111141059032230000","01631031000000d0"); // tm pt dependent var 3
  } else if (trainConfig == 56){ // TM variations
    cuts.AddCutCalo("82410113","1111153057032230000","01631031000000d0"); // CCMF nonlin+testbeam
    cuts.AddCutCalo("82410113","1111154057032230000","01631031000000d0"); // CMF nonlin+testbeam
    cuts.AddCutCalo("82410113","1111143057032230000","01631031000000d0"); // CCRF testbeam nonlin
    cuts.AddCutCalo("82410113","1111144057032230000","01631031000000d0"); // CRF testbeam nonlin
    cuts.AddCutCalo("82410113","1111102057032230000","01631031000000d0"); // testbeam nonlin
  } else if (trainConfig == 57){ //all default triggers
    cuts.AddCutCalo("82410113","1111141057032230000","01631031000000d0"); // default MB
    cuts.AddCutCalo("82452113","1111141057032230000","01631031000000d0"); // default EMC7
    cuts.AddCutCalo("82483113","1111141057032230000","01631031000000d0"); // default EG1
    cuts.AddCutCalo("82485113","1111141057032230000","01631031000000d0"); // default EG2

  //-----------------------------------------------------------------------------------------------
  // Systematics variations 40-60
  //-----------------------------------------------------------------------------------------------
  } else if (trainConfig == 60){ // nonlinearity variations
    cuts.AddCutCalo("84610113","1111142057032230000","01631031000000d0"); // CRF
    cuts.AddCutCalo("84610113","1111151057032230000","01631031000000d0"); // CCMF
    cuts.AddCutCalo("84610113","1111152057032230000","01631031000000d0"); // CMF
  } else if (trainConfig == 61){ // second set of variations CLUSTER
    cuts.AddCutCalo("84610113","1111141057022230000","01631031000000d0"); // min energy cluster variation 1  600 MeV
    cuts.AddCutCalo("84610113","1111141057042230000","01631031000000d0"); // min energy cluster variation 2  800 MeV
    cuts.AddCutCalo("84610113","1111141057052230000","01631031000000d0"); // min energy cluster variation 3  900 MeV
    cuts.AddCutCalo("84610113","1111141057032230000","01631031000000d0"); // min/max M02  0.1<M<0.5
    cuts.AddCutCalo("84610113","1111141057032200000","01631031000000d0"); // min/max M02  0.1<M<100
  } else if (trainConfig == 62){ // third set of variations CLUSTER
    cuts.AddCutCalo("84610113","1111141057032250000","01631031000000d0"); // min/max M02  0.1<M<0.3
    cuts.AddCutCalo("84610113","1111141057032260000","01631031000000d0"); // min/max M02  0.1<M<0.27
    cuts.AddCutCalo("84610113","1111141057031230000","01631031000000d0"); // min number of cells variation 1  1 cell
    cuts.AddCutCalo("84610113","1112141057032230000","01631031000000d0"); // only modules with TRD infront
    cuts.AddCutCalo("84610113","1111341057032230000","01631031000000d0"); // no modules with TRD infront
  } else if (trainConfig == 63){ // third set of variations MESON
    cuts.AddCutCalo("84610113","1111141057032230000","01633031000000d0"); // rapidity variation  y<0.6
    cuts.AddCutCalo("84610113","1111141057032230000","01634031000000d0"); // rapidity variation  y<0.5
    cuts.AddCutCalo("84610113","1111141057032230000","01631061000000d0"); // alpha meson variation 1   0<alpha<0.8
    cuts.AddCutCalo("84610113","1111141057032230000","01631051000000d0"); // alpha meson variation 2  0<alpha<0.75
  } else if (trainConfig == 64){ // opening angle variations
    cuts.AddCutCalo("84610113","1111141057032230000","0163103100000040"); // min opening angle 0.0152
    cuts.AddCutCalo("84610113","1111141057032230000","0163103100000060"); // min opening angle 0.017
    cuts.AddCutCalo("84610113","1111141057032230000","0163103100000070"); // min opening angle 0.016
    cuts.AddCutCalo("84610113","1111141057032230000","0163103100000080"); // min opening angle 0.018
    cuts.AddCutCalo("84610113","1111141057032230000","0163103100000090"); // min opening angle 0.018
  } else if (trainConfig == 65){ // TM variations
    cuts.AddCutCalo("84610113","1111141053032230000","01631031000000d0"); // fixed window
    cuts.AddCutCalo("84610113","1111141056032230000","01631031000000d0"); // tm pt dependent var 1
    cuts.AddCutCalo("84610113","1111141058032230000","01631031000000d0"); // tm pt dependent var 2
    cuts.AddCutCalo("84610113","1111141059032230000","01631031000000d0"); // tm pt dependent var 3
  } else if (trainConfig == 66){ // TM variations
    cuts.AddCutCalo("84610113","1111153057032230000","01631031000000d0"); // CCMF nonlin+testbeam
    cuts.AddCutCalo("84610113","1111154057032230000","01631031000000d0"); // CMF nonlin+testbeam
    cuts.AddCutCalo("84610113","1111143057032230000","01631031000000d0"); // CCRF testbeam nonlin
    cuts.AddCutCalo("84610113","1111144057032230000","01631031000000d0"); // CRF testbeam nonlin
    cuts.AddCutCalo("84610113","1111102057032230000","01631031000000d0"); // testbeam nonlin
  } else if (trainConfig == 67){ //all default triggers
    cuts.AddCutCalo("84610113","1111141057032230000","01631031000000d0"); // default MB
    cuts.AddCutCalo("84652113","1111141057032230000","01631031000000d0"); // default EMC7
    cuts.AddCutCalo("84683113","1111141057032230000","01631031000000d0"); // default EG1
    cuts.AddCutCalo("84685113","1111141057032230000","01631031000000d0"); // default EG2

  //-----------------------------------------------------------------------------------------------
  // Systematics variations 60-100
  //-----------------------------------------------------------------------------------------------
  } else if (trainConfig == 70){ // nonlinearity variations
    cuts.AddCutCalo("86010113","1111142057032230000","01631031000000d0"); // CRF
    cuts.AddCutCalo("86010113","1111151057032230000","01631031000000d0"); // CCMF
    cuts.AddCutCalo("86010113","1111152057032230000","01631031000000d0"); // CMF
  } else if (trainConfig == 71){ // second set of variations CLUSTER
    cuts.AddCutCalo("86010113","1111141057022230000","01631031000000d0"); // min energy cluster variation 1  600 MeV
    cuts.AddCutCalo("86010113","1111141057042230000","01631031000000d0"); // min energy cluster variation 2  800 MeV
    cuts.AddCutCalo("86010113","1111141057052230000","01631031000000d0"); // min energy cluster variation 3  900 MeV
    cuts.AddCutCalo("86010113","1111141057032230000","01631031000000d0"); // min/max M02  0.1<M<0.5
    cuts.AddCutCalo("86010113","1111141057032200000","01631031000000d0"); // min/max M02  0.1<M<100
  } else if (trainConfig == 72){ // third set of variations CLUSTER
    cuts.AddCutCalo("86010113","1111141057032250000","01631031000000d0"); // min/max M02  0.1<M<0.3
    cuts.AddCutCalo("86010113","1111141057032260000","01631031000000d0"); // min/max M02  0.1<M<0.27
    cuts.AddCutCalo("86010113","1111141057031230000","01631031000000d0"); // min number of cells variation 1  1 cell
    cuts.AddCutCalo("86010113","1112141057032230000","01631031000000d0"); // only modules with TRD infront
    cuts.AddCutCalo("86010113","1111341057032230000","01631031000000d0"); // no modules with TRD infront
  } else if (trainConfig == 73){ // third set of variations MESON
    cuts.AddCutCalo("86010113","1111141057032230000","01633031000000d0"); // rapidity variation  y<0.6
    cuts.AddCutCalo("86010113","1111141057032230000","01634031000000d0"); // rapidity variation  y<0.5
    cuts.AddCutCalo("86010113","1111141057032230000","01631061000000d0"); // alpha meson variation 1   0<alpha<0.8
    cuts.AddCutCalo("86010113","1111141057032230000","01631051000000d0"); // alpha meson variation 2  0<alpha<0.75
  } else if (trainConfig == 74){ // opening angle variations
    cuts.AddCutCalo("86010113","1111141057032230000","0163103100000040"); // min opening angle 0.0152
    cuts.AddCutCalo("86010113","1111141057032230000","0163103100000060"); // min opening angle 0.017
    cuts.AddCutCalo("86010113","1111141057032230000","0163103100000070"); // min opening angle 0.016
    cuts.AddCutCalo("86010113","1111141057032230000","0163103100000080"); // min opening angle 0.018
    cuts.AddCutCalo("86010113","1111141057032230000","0163103100000090"); // min opening angle 0.018
  } else if (trainConfig == 75){ // TM variations
    cuts.AddCutCalo("86010113","1111141053032230000","01631031000000d0"); // fixed window
    cuts.AddCutCalo("86010113","1111141056032230000","01631031000000d0"); // tm pt dependent var 1
    cuts.AddCutCalo("86010113","1111141058032230000","01631031000000d0"); // tm pt dependent var 2
    cuts.AddCutCalo("86010113","1111141059032230000","01631031000000d0"); // tm pt dependent var 3
  } else if (trainConfig == 76){ // TM variations
    cuts.AddCutCalo("86010113","1111153057032230000","01631031000000d0"); // CCMF nonlin+testbeam
    cuts.AddCutCalo("86010113","1111154057032230000","01631031000000d0"); // CMF nonlin+testbeam
    cuts.AddCutCalo("86010113","1111143057032230000","01631031000000d0"); // CCRF testbeam nonlin
    cuts.AddCutCalo("86010113","1111144057032230000","01631031000000d0"); // CRF testbeam nonlin
    cuts.AddCutCalo("86010113","1111102057032230000","01631031000000d0"); // testbeam nonlin
  } else if (trainConfig == 77){ //all default triggers
    cuts.AddCutCalo("86010113","1111141057032230000","01631031000000d0"); // default MB
    cuts.AddCutCalo("86052113","1111141057032230000","01631031000000d0"); // default EMC7
    cuts.AddCutCalo("86083113","1111141057032230000","01631031000000d0"); // default EG1
    cuts.AddCutCalo("86085113","1111141057032230000","01631031000000d0"); // default EG2

  //-----------------------------------------------------------------------------------------------
  //testing different event mix methods
  //-----------------------------------------------------------------------------------------------
  } else if (trainConfig == 81){
    cuts.AddCutCalo("80010113","1111141057032230000","01631031000000d0"); // default (V0 mult)
    cuts.AddCutCalo("80010113","1111141057032230000","02631031000000d0"); // using track mult
    cuts.AddCutCalo("80010113","1111141057032230000","09631031000000d0"); // using PtMax method

  //-----------------------------------------------------------------------------------------------
  // Systematics variations 0-20
  //-----------------------------------------------------------------------------------------------
  } else if (trainConfig == 90){ // nonlinearity variations
    cuts.AddCutCalo("80210113","11111420570322l0000","01631031000000d0"); // CRF
    cuts.AddCutCalo("80210113","11111510570322l0000","01631031000000d0"); // CCMF
    cuts.AddCutCalo("80210113","11111520570322l0000","01631031000000d0"); // CMF
  } else if (trainConfig == 91){ // second set of variations CLUSTER
    cuts.AddCutCalo("80210113","11111410570222l0000","01631031000000d0"); // min energy cluster variation 1  600 MeV
    cuts.AddCutCalo("80210113","11111410570422l0000","01631031000000d0"); // min energy cluster variation 2  800 MeV
    cuts.AddCutCalo("80210113","11111410570522l0000","01631031000000d0"); // min energy cluster variation 3  900 MeV
    cuts.AddCutCalo("80210113","11111410570322l0000","01631031000000d0"); // min/max M02  0.1<M<0.5
  } else if (trainConfig == 92){ // third set of variations CLUSTER
    cuts.AddCutCalo("80210113","1111141057032250000","01631031000000d0"); // min/max M02  0.1<M<0.3
    cuts.AddCutCalo("80210113","1111141057032260000","01631031000000d0"); // min/max M02  0.1<M<0.27
    cuts.AddCutCalo("80210113","11111410570312l0000","01631031000000d0"); // min number of cells variation 1  1 cell
    cuts.AddCutCalo("80210113","11121410570322l0000","01631031000000d0"); // only modules with TRD infront
    cuts.AddCutCalo("80210113","11113410570322l0000","01631031000000d0"); // no modules with TRD infront
  } else if (trainConfig == 93){ // third set of variations MESON
    cuts.AddCutCalo("80210113","11111410570322l0000","01633031000000d0"); // rapidity variation  y<0.6
    cuts.AddCutCalo("80210113","11111410570322l0000","01634031000000d0"); // rapidity variation  y<0.5
    cuts.AddCutCalo("80210113","11111410570322l0000","01631061000000d0"); // alpha meson variation 1   0<alpha<0.8
    cuts.AddCutCalo("80210113","11111410570322l0000","01631051000000d0"); // alpha meson variation 2  0<alpha<0.75
  } else if (trainConfig == 94){ // opening angle variations
    cuts.AddCutCalo("80210113","11111410570322l0000","0163103100000040"); // min opening angle 0.0152
    cuts.AddCutCalo("80210113","11111410570322l0000","0163103100000060"); // min opening angle 0.017
    cuts.AddCutCalo("80210113","11111410570322l0000","0163103100000070"); // min opening angle 0.016
    cuts.AddCutCalo("80210113","11111410570322l0000","0163103100000080"); // min opening angle 0.018
    cuts.AddCutCalo("80210113","11111410570322l0000","0163103100000090"); // min opening angle 0.018
  } else if (trainConfig == 95){ // TM variations
    cuts.AddCutCalo("80210113","11111410530322l0000","01631031000000d0"); // fixed window
    cuts.AddCutCalo("80210113","11111410560322l0000","01631031000000d0"); // tm pt dependent var 1
    cuts.AddCutCalo("80210113","11111410580322l0000","01631031000000d0"); // tm pt dependent var 2
    cuts.AddCutCalo("80210113","11111410590322l0000","01631031000000d0"); // tm pt dependent var 3
  } else if (trainConfig == 96){ // TM variations
    cuts.AddCutCalo("80210113","11111530570322l0000","01631031000000d0"); // CCMF nonlin+testbeam
    cuts.AddCutCalo("80210113","11111540570322l0000","01631031000000d0"); // CMF nonlin+testbeam
    cuts.AddCutCalo("80210113","11111430570322l0000","01631031000000d0"); // CCRF testbeam nonlin
    cuts.AddCutCalo("80210113","11111440570322l0000","01631031000000d0"); // CRF testbeam nonlin
    cuts.AddCutCalo("80210113","11111020570322l0000","01631031000000d0"); // testbeam nonlin
  } else if (trainConfig == 97){
    cuts.AddCutCalo("80210113","1111141057032230000","01631031000000d0"); // M02 pt dep with new std: 0.32, 0.0072, 0.5
    cuts.AddCutCalo("80210113","11111410570322d0000","01631031000000d0"); // M02, pt dep with  0.27, 0.0072, 0.4
    cuts.AddCutCalo("80210113","11111410570322e0000","01631031000000d0"); // M02, pt dep with  0.31, 0.0072, 0.5
    cuts.AddCutCalo("80210113","11111410570322f0000","01631031000000d0"); // M02, pt dep with  0.36, 0.0072, 0.7
  } else if (trainConfig == 98){
    cuts.AddCutCalo("80210113","11111410570322m0000","01631031000000d0"); // M02, pt dep with  0.32, 0.0152, 0.5
    cuts.AddCutCalo("80210113","11111410570322g0000","01631031000000d0"); // M02, pt dep with  0.37, 0.0072, 0.7
    cuts.AddCutCalo("80210113","11111410570322h0000","01631031000000d0"); // M02, pT-dep with  0.30, 0.0072, 0.5
    cuts.AddCutCalo("80210113","11111410570322i0000","01631031000000d0"); // M02, pT-dep with  0.35, 0.0072, 0.7
  } else if (trainConfig == 99){
    cuts.AddCutCalo("80210113","11111410570322j0000","01631031000000d0"); // M02, pT-dep with  0.25, 0.0072, 0.39
    cuts.AddCutCalo("80210113","11111410570322r0000","01631031000000d0"); // M02, pT-dep with  0.25, 0.0072, 0.5
    cuts.AddCutCalo("80210113","11111410570322s0000","01631031000000d0"); // M02, pT-dep with  0.32, 0.0238, 0.7
    cuts.AddCutCalo("80210113","11111410570322n0000","01631031000000d0"); // M02, pT-dep with  0.32, 0.0238, 0.7

  //-----------------------------------------------------------------------------------------------
  // Systematics variations 20-40
  //-----------------------------------------------------------------------------------------------
  } else if (trainConfig == 100){ // nonlinearity variations
    cuts.AddCutCalo("82410113","11111420570322l0000","01631031000000d0"); // CRF
    cuts.AddCutCalo("82410113","11111510570322l0000","01631031000000d0"); // CCMF
    cuts.AddCutCalo("82410113","11111520570322l0000","01631031000000d0"); // CMF
  } else if (trainConfig == 101){ // second set of variations CLUSTER
    cuts.AddCutCalo("82410113","11111410570222l0000","01631031000000d0"); // min energy cluster variation 1  600 MeV
    cuts.AddCutCalo("82410113","11111410570422l0000","01631031000000d0"); // min energy cluster variation 2  800 MeV
    cuts.AddCutCalo("82410113","11111410570522l0000","01631031000000d0"); // min energy cluster variation 3  900 MeV
    cuts.AddCutCalo("82410113","11111410570322l0000","01631031000000d0"); // min/max M02  0.1<M<0.5
  } else if (trainConfig == 102){ // third set of variations CLUSTER
    cuts.AddCutCalo("82410113","1111141057032250000","01631031000000d0"); // min/max M02  0.1<M<0.3
    cuts.AddCutCalo("82410113","1111141057032260000","01631031000000d0"); // min/max M02  0.1<M<0.27
    cuts.AddCutCalo("82410113","11111410570312l0000","01631031000000d0"); // min number of cells variation 1  1 cell
    cuts.AddCutCalo("82410113","11121410570322l0000","01631031000000d0"); // only modules with TRD infront
    cuts.AddCutCalo("82410113","11113410570322l0000","01631031000000d0"); // no modules with TRD infront
  } else if (trainConfig == 103){ // third set of variations MESON
    cuts.AddCutCalo("82410113","11111410570322l0000","01633031000000d0"); // rapidity variation  y<0.6
    cuts.AddCutCalo("82410113","11111410570322l0000","01634031000000d0"); // rapidity variation  y<0.5
    cuts.AddCutCalo("82410113","11111410570322l0000","01631061000000d0"); // alpha meson variation 1   0<alpha<0.8
    cuts.AddCutCalo("82410113","11111410570322l0000","01631051000000d0"); // alpha meson variation 2  0<alpha<0.75
  } else if (trainConfig == 104){ // opening angle variations
    cuts.AddCutCalo("82410113","11111410570322l0000","0163103100000040"); // min opening angle 0.0152
    cuts.AddCutCalo("82410113","11111410570322l0000","0163103100000060"); // min opening angle 0.017
    cuts.AddCutCalo("82410113","11111410570322l0000","0163103100000070"); // min opening angle 0.016
    cuts.AddCutCalo("82410113","11111410570322l0000","0163103100000080"); // min opening angle 0.018
    cuts.AddCutCalo("82410113","11111410570322l0000","0163103100000090"); // min opening angle 0.018
  } else if (trainConfig == 105){ // TM variations
    cuts.AddCutCalo("82410113","11111410530322l0000","01631031000000d0"); // fixed window
    cuts.AddCutCalo("82410113","11111410560322l0000","01631031000000d0"); // tm pt dependent var 1
    cuts.AddCutCalo("82410113","11111410580322l0000","01631031000000d0"); // tm pt dependent var 2
    cuts.AddCutCalo("82410113","11111410590322l0000","01631031000000d0"); // tm pt dependent var 3
  } else if (trainConfig == 106){ // TM variations
    cuts.AddCutCalo("82410113","11111530570322l0000","01631031000000d0"); // CCMF nonlin+testbeam
    cuts.AddCutCalo("82410113","11111540570322l0000","01631031000000d0"); // CMF nonlin+testbeam
    cuts.AddCutCalo("82410113","11111430570322l0000","01631031000000d0"); // CCRF testbeam nonlin
    cuts.AddCutCalo("82410113","11111440570322l0000","01631031000000d0"); // CRF testbeam nonlin
    cuts.AddCutCalo("82410113","11111020570322l0000","01631031000000d0"); // testbeam nonlin
  } else if (trainConfig == 107){
    cuts.AddCutCalo("82410113","1111141057032230000","01631031000000d0"); // M02 pt dep with new std: 0.32, 0.0072, 0.5
    cuts.AddCutCalo("82410113","11111410570322d0000","01631031000000d0"); // M02, pt dep with  0.27, 0.0072, 0.4
    cuts.AddCutCalo("82410113","11111410570322e0000","01631031000000d0"); // M02, pt dep with  0.31, 0.0072, 0.5
    cuts.AddCutCalo("82410113","11111410570322f0000","01631031000000d0"); // M02, pt dep with  0.36, 0.0072, 0.7
  } else if (trainConfig == 108){
    cuts.AddCutCalo("82410113","11111410570322m0000","01631031000000d0"); // M02, pt dep with  0.32, 0.0152, 0.5
    cuts.AddCutCalo("82410113","11111410570322g0000","01631031000000d0"); // M02, pt dep with  0.37, 0.0072, 0.7
    cuts.AddCutCalo("82410113","11111410570322h0000","01631031000000d0"); // M02, pT-dep with  0.30, 0.0072, 0.5
    cuts.AddCutCalo("82410113","11111410570322i0000","01631031000000d0"); // M02, pT-dep with  0.35, 0.0072, 0.7
  } else if (trainConfig == 109){
    cuts.AddCutCalo("82410113","11111410570322j0000","01631031000000d0"); // M02, pT-dep with  0.25, 0.0072, 0.39
    cuts.AddCutCalo("82410113","11111410570322r0000","01631031000000d0"); // M02, pT-dep with  0.25, 0.0072, 0.5
    cuts.AddCutCalo("82410113","11111410570322s0000","01631031000000d0"); // M02, pT-dep with  0.32, 0.0238, 0.7
    cuts.AddCutCalo("82410113","11111410570322n0000","01631031000000d0"); // M02, pT-dep with  0.32, 0.0238, 0.7

  //-----------------------------------------------------------------------------------------------
  // Systematics variations 40-60
  //-----------------------------------------------------------------------------------------------
  } else if (trainConfig == 110){ // nonlinearity variations
    cuts.AddCutCalo("84610113","11111420570322l0000","01631031000000d0"); // CRF
    cuts.AddCutCalo("84610113","11111510570322l0000","01631031000000d0"); // CCMF
    cuts.AddCutCalo("84610113","11111520570322l0000","01631031000000d0"); // CMF
  } else if (trainConfig == 111){ // second set of variations CLUSTER
    cuts.AddCutCalo("84610113","11111410570222l0000","01631031000000d0"); // min energy cluster variation 1  600 MeV
    cuts.AddCutCalo("84610113","11111410570422l0000","01631031000000d0"); // min energy cluster variation 2  800 MeV
    cuts.AddCutCalo("84610113","11111410570522l0000","01631031000000d0"); // min energy cluster variation 3  900 MeV
    cuts.AddCutCalo("84610113","11111410570322l0000","01631031000000d0"); // min/max M02  0.1<M<0.5
  } else if (trainConfig == 112){ // third set of variations CLUSTER
    cuts.AddCutCalo("84610113","1111141057032250000","01631031000000d0"); // min/max M02  0.1<M<0.3
    cuts.AddCutCalo("84610113","1111141057032260000","01631031000000d0"); // min/max M02  0.1<M<0.27
    cuts.AddCutCalo("84610113","11111410570312l0000","01631031000000d0"); // min number of cells variation 1  1 cell
    cuts.AddCutCalo("84610113","11121410570322l0000","01631031000000d0"); // only modules with TRD infront
    cuts.AddCutCalo("84610113","11113410570322l0000","01631031000000d0"); // no modules with TRD infront
  } else if (trainConfig == 113){ // third set of variations MESON
    cuts.AddCutCalo("84610113","11111410570322l0000","01633031000000d0"); // rapidity variation  y<0.6
    cuts.AddCutCalo("84610113","11111410570322l0000","01634031000000d0"); // rapidity variation  y<0.5
    cuts.AddCutCalo("84610113","11111410570322l0000","01631061000000d0"); // alpha meson variation 1   0<alpha<0.8
    cuts.AddCutCalo("84610113","11111410570322l0000","01631051000000d0"); // alpha meson variation 2  0<alpha<0.75
  } else if (trainConfig == 114){ // opening angle variations
    cuts.AddCutCalo("84610113","11111410570322l0000","0163103100000040"); // min opening angle 0.0152
    cuts.AddCutCalo("84610113","11111410570322l0000","0163103100000060"); // min opening angle 0.017
    cuts.AddCutCalo("84610113","11111410570322l0000","0163103100000070"); // min opening angle 0.016
    cuts.AddCutCalo("84610113","11111410570322l0000","0163103100000080"); // min opening angle 0.018
    cuts.AddCutCalo("84610113","11111410570322l0000","0163103100000090"); // min opening angle 0.018
  } else if (trainConfig == 115){ // TM variations
    cuts.AddCutCalo("84610113","11111410530322l0000","01631031000000d0"); // fixed window
    cuts.AddCutCalo("84610113","11111410560322l0000","01631031000000d0"); // tm pt dependent var 1
    cuts.AddCutCalo("84610113","11111410580322l0000","01631031000000d0"); // tm pt dependent var 2
    cuts.AddCutCalo("84610113","11111410590322l0000","01631031000000d0"); // tm pt dependent var 3
  } else if (trainConfig == 116){ // TM variations
    cuts.AddCutCalo("84610113","11111530570322l0000","01631031000000d0"); // CCMF nonlin+testbeam
    cuts.AddCutCalo("84610113","11111540570322l0000","01631031000000d0"); // CMF nonlin+testbeam
    cuts.AddCutCalo("84610113","11111430570322l0000","01631031000000d0"); // CCRF testbeam nonlin
    cuts.AddCutCalo("84610113","11111440570322l0000","01631031000000d0"); // CRF testbeam nonlin
    cuts.AddCutCalo("84610113","11111020570322l0000","01631031000000d0"); // testbeam nonlin
  } else if (trainConfig == 117){
    cuts.AddCutCalo("84610113","1111141057032230000","01631031000000d0"); // M02 pt dep with new std: 0.32, 0.0072, 0.5
    cuts.AddCutCalo("84610113","11111410570322d0000","01631031000000d0"); // M02, pt dep with  0.27, 0.0072, 0.4
    cuts.AddCutCalo("84610113","11111410570322e0000","01631031000000d0"); // M02, pt dep with  0.31, 0.0072, 0.5
    cuts.AddCutCalo("84610113","11111410570322f0000","01631031000000d0"); // M02, pt dep with  0.36, 0.0072, 0.7
  } else if (trainConfig == 118){
    cuts.AddCutCalo("84610113","11111410570322m0000","01631031000000d0"); // M02, pt dep with  0.32, 0.0152, 0.5
    cuts.AddCutCalo("84610113","11111410570322g0000","01631031000000d0"); // M02, pt dep with  0.37, 0.0072, 0.7
    cuts.AddCutCalo("84610113","11111410570322h0000","01631031000000d0"); // M02, pT-dep with  0.30, 0.0072, 0.5
    cuts.AddCutCalo("84610113","11111410570322i0000","01631031000000d0"); // M02, pT-dep with  0.35, 0.0072, 0.7
  } else if (trainConfig == 119){
    cuts.AddCutCalo("84610113","11111410570322j0000","01631031000000d0"); // M02, pT-dep with  0.25, 0.0072, 0.39
    cuts.AddCutCalo("84610113","11111410570322r0000","01631031000000d0"); // M02, pT-dep with  0.25, 0.0072, 0.5
    cuts.AddCutCalo("84610113","11111410570322s0000","01631031000000d0"); // M02, pT-dep with  0.32, 0.0238, 0.7
    cuts.AddCutCalo("84610113","11111410570322n0000","01631031000000d0"); // M02, pT-dep with  0.32, 0.0238, 0.7

  //-----------------------------------------------------------------------------------------------
  // Systematics variations 60-100
  //-----------------------------------------------------------------------------------------------
  } else if (trainConfig == 120){ // nonlinearity variations
    cuts.AddCutCalo("86010113","11111420570322l0000","01631031000000d0"); // CRF
    cuts.AddCutCalo("86010113","11111510570322l0000","01631031000000d0"); // CCMF
    cuts.AddCutCalo("86010113","11111520570322l0000","01631031000000d0"); // CMF
  } else if (trainConfig == 121){ // second set of variations CLUSTER
    cuts.AddCutCalo("86010113","11111410570222l0000","01631031000000d0"); // min energy cluster variation 1  600 MeV
    cuts.AddCutCalo("86010113","11111410570422l0000","01631031000000d0"); // min energy cluster variation 2  800 MeV
    cuts.AddCutCalo("86010113","11111410570522l0000","01631031000000d0"); // min energy cluster variation 3  900 MeV
    cuts.AddCutCalo("86010113","11111410570322l0000","01631031000000d0"); // min/max M02  0.1<M<0.5
  } else if (trainConfig == 122){ // third set of variations CLUSTER
    cuts.AddCutCalo("86010113","1111141057032250000","01631031000000d0"); // min/max M02  0.1<M<0.3
    cuts.AddCutCalo("86010113","1111141057032260000","01631031000000d0"); // min/max M02  0.1<M<0.27
    cuts.AddCutCalo("86010113","11111410570312l0000","01631031000000d0"); // min number of cells variation 1  1 cell
    cuts.AddCutCalo("86010113","11121410570322l0000","01631031000000d0"); // only modules with TRD infront
    cuts.AddCutCalo("86010113","11113410570322l0000","01631031000000d0"); // no modules with TRD infront
  } else if (trainConfig == 123){ // third set of variations MESON
    cuts.AddCutCalo("86010113","11111410570322l0000","01633031000000d0"); // rapidity variation  y<0.6
    cuts.AddCutCalo("86010113","11111410570322l0000","01634031000000d0"); // rapidity variation  y<0.5
    cuts.AddCutCalo("86010113","11111410570322l0000","01631061000000d0"); // alpha meson variation 1   0<alpha<0.8
    cuts.AddCutCalo("86010113","11111410570322l0000","01631051000000d0"); // alpha meson variation 2  0<alpha<0.75
  } else if (trainConfig == 124){ // opening angle variations
    cuts.AddCutCalo("86010113","11111410570322l0000","0163103100000040"); // min opening angle 0.0152
    cuts.AddCutCalo("86010113","11111410570322l0000","0163103100000060"); // min opening angle 0.017
    cuts.AddCutCalo("86010113","11111410570322l0000","0163103100000070"); // min opening angle 0.016
    cuts.AddCutCalo("86010113","11111410570322l0000","0163103100000080"); // min opening angle 0.018
    cuts.AddCutCalo("86010113","11111410570322l0000","0163103100000090"); // min opening angle 0.018
  } else if (trainConfig == 125){ // TM variations
    cuts.AddCutCalo("86010113","11111410530322l0000","01631031000000d0"); // fixed window
    cuts.AddCutCalo("86010113","11111410560322l0000","01631031000000d0"); // tm pt dependent var 1
    cuts.AddCutCalo("86010113","11111410580322l0000","01631031000000d0"); // tm pt dependent var 2
    cuts.AddCutCalo("86010113","11111410590322l0000","01631031000000d0"); // tm pt dependent var 3
  } else if (trainConfig == 126){ // TM variations
    cuts.AddCutCalo("86010113","11111530570322l0000","01631031000000d0"); // CCMF nonlin+testbeam
    cuts.AddCutCalo("86010113","11111540570322l0000","01631031000000d0"); // CMF nonlin+testbeam
    cuts.AddCutCalo("86010113","11111430570322l0000","01631031000000d0"); // CCRF testbeam nonlin
    cuts.AddCutCalo("86010113","11111440570322l0000","01631031000000d0"); // CRF testbeam nonlin
    cuts.AddCutCalo("86010113","11111020570322l0000","01631031000000d0"); // testbeam nonlin
  } else if (trainConfig == 127){
    cuts.AddCutCalo("86010113","1111141057032230000","01631031000000d0"); // M02 pt dep with new std: 0.32, 0.0072, 0.5
    cuts.AddCutCalo("86010113","11111410570322d0000","01631031000000d0"); // M02, pt dep with  0.27, 0.0072, 0.4
    cuts.AddCutCalo("86010113","11111410570322e0000","01631031000000d0"); // M02, pt dep with  0.31, 0.0072, 0.5
    cuts.AddCutCalo("86010113","11111410570322f0000","01631031000000d0"); // M02, pt dep with  0.36, 0.0072, 0.7
  } else if (trainConfig == 128){
    cuts.AddCutCalo("86010113","11111410570322m0000","01631031000000d0"); // M02, pt dep with  0.32, 0.0152, 0.5
    cuts.AddCutCalo("86010113","11111410570322g0000","01631031000000d0"); // M02, pt dep with  0.37, 0.0072, 0.7
    cuts.AddCutCalo("86010113","11111410570322h0000","01631031000000d0"); // M02, pT-dep with  0.30, 0.0072, 0.5
    cuts.AddCutCalo("86010113","11111410570322i0000","01631031000000d0"); // M02, pT-dep with  0.35, 0.0072, 0.7
  } else if (trainConfig == 129){
    cuts.AddCutCalo("86010113","11111410570322j0000","01631031000000d0"); // M02, pT-dep with  0.25, 0.0072, 0.39
    cuts.AddCutCalo("86010113","11111410570322r0000","01631031000000d0"); // M02, pT-dep with  0.25, 0.0072, 0.5
    cuts.AddCutCalo("86010113","11111410570322s0000","01631031000000d0"); // M02, pT-dep with  0.32, 0.0238, 0.7
    cuts.AddCutCalo("86010113","11111410570322n0000","01631031000000d0"); // M02, pT-dep with  0.32, 0.0238, 0.7

  // ===============================================================================================
  // Run 2 data EMC clusters pPb 5TeV
  // ===============================================================================================
  } else if (trainConfig == 200){ // EMCAL clusters standard cuts,
    cuts.AddCutCalo("80010113","1111100057032230000","01631031000000d0"); // 0-100
  } else if (trainConfig == 201){ // EMCAL clusters standard cuts, no nonlin, open timing
    cuts.AddCutCalo("80010113","1111100017032230000","01631031000000d0"); // INT7
  } else if (trainConfig == 202){ // EMCAL clusters standard cuts, no nonlin, +-50ns, trigger
    cuts.AddCutCalo("80052113","1111100057032230000","01631031000000d0"); // EMC7
    cuts.AddCutCalo("80085113","1111100057032230000","01631031000000d0"); // EG2
    cuts.AddCutCalo("80083113","1111100057032230000","01631031000000d0"); // EG1
  } else if (trainConfig == 203){ // EMCAL clusters standard cuts, no nonlin, open timing, trigger
    cuts.AddCutCalo("80052113","1111100017032230000","01631031000000d0"); // EMC7
    cuts.AddCutCalo("80085113","1111100017032230000","01631031000000d0"); // EG2
    cuts.AddCutCalo("80083113","1111100017032230000","01631031000000d0"); // EG1
  } else if (trainConfig == 204){ // EMCAL clusters standard cuts, triggers, no nonlin, +-50ns
    cuts.AddCutCalo("80110113","1111100057032230000","01631031000000d0"); // 0-10
    cuts.AddCutCalo("81210113","1111100057032230000","01631031000000d0"); // 10-20
    cuts.AddCutCalo("82410113","1111100057032230000","01631031000000d0"); // 20-40
    cuts.AddCutCalo("84610113","1111100057032230000","01631031000000d0"); // 40-60
    cuts.AddCutCalo("86810113","1111100057032230000","01631031000000d0"); // 60-80
    cuts.AddCutCalo("88010113","1111100057032230000","01631031000000d0"); // 80-100
  } else if (trainConfig == 205){ // EMCAL clusters standard cuts, triggers, no nonlin, open timing
    cuts.AddCutCalo("80110113","1111100017032230000","01631031000000d0"); // 0-10
    cuts.AddCutCalo("81210113","1111100017032230000","01631031000000d0"); // 10-20
    cuts.AddCutCalo("82410113","1111100017032230000","01631031000000d0"); // 20-40
    cuts.AddCutCalo("84610113","1111100017032230000","01631031000000d0"); // 40-60
    cuts.AddCutCalo("86810113","1111100017032230000","01631031000000d0"); // 60-80
    cuts.AddCutCalo("88010113","1111100017032230000","01631031000000d0"); // 80-100
  } else if (trainConfig == 206){ // EMCAL clusters standard cuts, no nonlin, +-50ns
    cuts.AddCutCalo("80210113","1111100057032230000","01631031000000d0"); // 0-20
    cuts.AddCutCalo("a0110113","1111100057032230000","01631031000000d0"); // 0-5
    cuts.AddCutCalo("a1210113","1111100057032230000","01631031000000d0"); // 5-10
    cuts.AddCutCalo("86010113","1111100057032230000","01631031000000d0"); // 60-100
  } else if (trainConfig == 207){ // EMCAL clusters standard cuts, no nonlin, open timing
    cuts.AddCutCalo("80210113","1111100017032230000","01631031000000d0"); // 0-20
    cuts.AddCutCalo("a0110113","1111100017032230000","01631031000000d0"); // 0-5
    cuts.AddCutCalo("a1210113","1111100017032230000","01631031000000d0"); // 5-10
    cuts.AddCutCalo("86010113","1111100017032230000","01631031000000d0"); // 60-100

  // non lin variations
  } else if (trainConfig == 208){ // EMCAL clusters standard cuts, non lin variations
    cuts.AddCutCalo("80010113","1111141057032230000","01631031000000d0"); // 0-100
    cuts.AddCutCalo("80010113","1111142057032230000","01631031000000d0"); // 0-100
    cuts.AddCutCalo("80010113","1111151057032230000","01631031000000d0"); // 0-100
    cuts.AddCutCalo("80010113","1111152057032230000","01631031000000d0"); // 0-100
    cuts.AddCutCalo("80010113","1111155057032230000","01631031000000d0"); // 0-100
    cuts.AddCutCalo("80010113","1111156057032230000","01631031000000d0"); // 0-100
  } else if (trainConfig == 209){ // EMCAL clusters standard cuts, non lin variations
    cuts.AddCutCalo("80010113","1111102057032230000","01631031000000d0"); // 0-100
    cuts.AddCutCalo("80010113","1111100057032230000","01631031000000d0"); // 0-100
  } else if (trainConfig == 210){ // EMCAL clusters standard cuts, non lin variations
    cuts.AddCutCalo("80010113","1111141017032230000","01631031000000d0"); // 0-100
    cuts.AddCutCalo("80010113","1111142017032230000","01631031000000d0"); // 0-100
    cuts.AddCutCalo("80010113","1111151017032230000","01631031000000d0"); // 0-100
    cuts.AddCutCalo("80010113","1111152017032230000","01631031000000d0"); // 0-100
  } else if (trainConfig == 211){ // EMCAL clusters standard cuts, non lin variations
    cuts.AddCutCalo("80010113","1111102017032230000","01631031000000d0"); // 0-100
    cuts.AddCutCalo("80010113","1111100017032230000","01631031000000d0"); // 0-100

  // non lin variations
  } else if (trainConfig == 212){ // EMCAL clusters standard cuts, cent vars, non lin variations
    cuts.AddCutCalo("80110113","1111141057032230000","01631031000000d0"); // 0-10
    cuts.AddCutCalo("81210113","1111141057032230000","01631031000000d0"); // 10-20
    cuts.AddCutCalo("82410113","1111141057032230000","01631031000000d0"); // 20-40
    cuts.AddCutCalo("84610113","1111141057032230000","01631031000000d0"); // 40-60
    cuts.AddCutCalo("86810113","1111141057032230000","01631031000000d0"); // 60-80
    cuts.AddCutCalo("88010113","1111141057032230000","01631031000000d0"); // 80-100
  } else if (trainConfig == 213){ // EMCAL clusters standard cuts, cent vars, non lin variations
    cuts.AddCutCalo("80210113","1111141057032230000","01631031000000d0"); // 0-20
    cuts.AddCutCalo("86010113","1111141057032230000","01631031000000d0"); // 60-100
    cuts.AddCutCalo("a0110113","1111141057032230000","01631031000000d0"); // 0-5
    cuts.AddCutCalo("a1210113","1111141057032230000","01631031000000d0"); // 5-10
  } else if (trainConfig == 214){ // EMCAL clusters standard cuts, cent vars, non lin variations
    cuts.AddCutCalo("80110113","1111151057032230000","01631031000000d0"); // 0-10
    cuts.AddCutCalo("81210113","1111151057032230000","01631031000000d0"); // 10-20
    cuts.AddCutCalo("82410113","1111151057032230000","01631031000000d0"); // 20-40
    cuts.AddCutCalo("84610113","1111151057032230000","01631031000000d0"); // 40-60
    cuts.AddCutCalo("86810113","1111151057032230000","01631031000000d0"); // 60-80
    cuts.AddCutCalo("88010113","1111151057032230000","01631031000000d0"); // 80-100
  } else if (trainConfig == 215){ // EMCAL clusters standard cuts, cent vars, non lin variations
    cuts.AddCutCalo("80210113","1111151057032230000","01631031000000d0"); // 0-20
    cuts.AddCutCalo("86010113","1111151057032230000","01631031000000d0"); // 60-100
    cuts.AddCutCalo("a0110113","1111151057032230000","01631031000000d0"); // 0-5
    cuts.AddCutCalo("a1210113","1111151057032230000","01631031000000d0"); // 5-10
  } else if (trainConfig == 216){ // EMCAL clusters standard cuts, dir gamma cut
    cuts.AddCutCalo("80010113","11111410570322l0000","01631031000000d0"); // 0-100
    cuts.AddCutCalo("80010113","11111420570322l0000","01631031000000d0"); // 0-100
    cuts.AddCutCalo("80010113","11111510570322l0000","01631031000000d0"); // 0-100
    cuts.AddCutCalo("80010113","11111520570322l0000","01631031000000d0"); // 0-100
    cuts.AddCutCalo("80010113","11111020570322l0000","01631031000000d0"); // 0-100
  } else if (trainConfig == 217){ // EMCAL clusters standard cuts, dir gamma cut
    cuts.AddCutCalo("80110113","11111410570322l0000","01631031000000d0");
    cuts.AddCutCalo("81210113","11111410570322l0000","01631031000000d0");
    cuts.AddCutCalo("82410113","11111410570322l0000","01631031000000d0");
    cuts.AddCutCalo("84610113","11111410570322l0000","01631031000000d0");
    cuts.AddCutCalo("86810113","11111410570322l0000","01631031000000d0");
    cuts.AddCutCalo("88010113","11111410570322l0000","01631031000000d0");
  } else if (trainConfig == 218){ // EMCAL clusters standard cuts, dir gamma cut
    cuts.AddCutCalo("80210113","11111410570322l0000","01631031000000d0"); // 0-20
    cuts.AddCutCalo("86010113","11111410570322l0000","01631031000000d0"); // 60-100
    cuts.AddCutCalo("a0110113","11111410570322l0000","01631031000000d0"); // 0-5
    cuts.AddCutCalo("a1210113","11111410570322l0000","01631031000000d0"); // 5-10
  } else if (trainConfig == 219){ // EMCAL clusters standard cuts, dir gamma cut
    cuts.AddCutCalo("80110113","11111510570322l0000","01631031000000d0");
    cuts.AddCutCalo("81210113","11111510570322l0000","01631031000000d0");
    cuts.AddCutCalo("82410113","11111510570322l0000","01631031000000d0");
    cuts.AddCutCalo("84610113","11111510570322l0000","01631031000000d0");
    cuts.AddCutCalo("86810113","11111510570322l0000","01631031000000d0");
    cuts.AddCutCalo("88010113","11111510570322l0000","01631031000000d0");
  } else if (trainConfig == 220){ // EMCAL clusters standard cuts, dir gamma cut
    cuts.AddCutCalo("80210113","11111510570322l0000","01631031000000d0"); // 0-20
    cuts.AddCutCalo("86010113","11111510570322l0000","01631031000000d0"); // 60-100
    cuts.AddCutCalo("a0110113","11111510570322l0000","01631031000000d0"); // 0-5
    cuts.AddCutCalo("a1210113","11111510570322l0000","01631031000000d0"); // 5-10

  // Nonlin variations with different min E cut
  } else if (trainConfig == 221){ // EMCAL clusters, non lin variations, min E = 0.6
    cuts.AddCutCalo("80010113","1111100017022230000","01631031000000d0"); // 0-100
    cuts.AddCutCalo("80010113","1111141017022230000","01631031000000d0"); // 0-100
    cuts.AddCutCalo("80010113","1111142017022230000","01631031000000d0"); // 0-100
    cuts.AddCutCalo("80010113","1111151017022230000","01631031000000d0"); // 0-100
    cuts.AddCutCalo("80010113","1111152017022230000","01631031000000d0"); // 0-100
  } else if (trainConfig == 222){ // EMCAL clusters, non lin variations, min E = 0.65
    cuts.AddCutCalo("80010113","11111000170b2230000","01631031000000d0"); // 0-100
    cuts.AddCutCalo("80010113","11111410170b2230000","01631031000000d0"); // 0-100
    cuts.AddCutCalo("80010113","11111420170b2230000","01631031000000d0"); // 0-100
    cuts.AddCutCalo("80010113","11111510170b2230000","01631031000000d0"); // 0-100
    cuts.AddCutCalo("80010113","11111520170b2230000","01631031000000d0"); // 0-100
  } else if (trainConfig == 223){ // EMCAL clusters, non lin variations, min E = 0.675
    cuts.AddCutCalo("80010113","11111000170c2230000","01631031000000d0"); // 0-100
    cuts.AddCutCalo("80010113","11111410170c2230000","01631031000000d0"); // 0-100
    cuts.AddCutCalo("80010113","11111420170c2230000","01631031000000d0"); // 0-100
    cuts.AddCutCalo("80010113","11111510170c2230000","01631031000000d0"); // 0-100
    cuts.AddCutCalo("80010113","11111520170c2230000","01631031000000d0"); // 0-100
  } else if (trainConfig == 224){ // EMCAL clusters, non lin variations, min E = 0.625
    cuts.AddCutCalo("80010113","11111000170d2230000","01631031000000d0"); // 0-100
    cuts.AddCutCalo("80010113","11111410170d2230000","01631031000000d0"); // 0-100
    cuts.AddCutCalo("80010113","11111420170d2230000","01631031000000d0"); // 0-100
    cuts.AddCutCalo("80010113","11111510170d2230000","01631031000000d0"); // 0-100
    cuts.AddCutCalo("80010113","11111520170d2230000","01631031000000d0"); // 0-100

  // Small centrality intervals
  } else if (trainConfig == 225){
    cuts.AddCutCalo("c0110113","1111141057032230000","01631031000000d0"); // 0-1
    cuts.AddCutCalo("c0210113","1111141057032230000","01631031000000d0"); // 0-2
    cuts.AddCutCalo("c0310113","1111141057032230000","01631031000000d0"); // 0-3
    cuts.AddCutCalo("c0410113","1111141057032230000","01631031000000d0"); // 0-4
  // AOD validation
  } else if (trainConfig == 226){
    cuts.AddCutCalo("80010113","1111141057032230000","01631031000000d0"); // 0-100
  // pPb 8 TeV EPOS+PythiaJets JJ simulation QA
  } else if (trainConfig == 227){ // variations when using special JJ case (Jets+MB)
    cuts.AddCutCalo("80010143","1111100017032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("80083143","1111100017032230000","01631031000000d0"); // EG1
    cuts.AddCutCalo("80085143","1111100017032230000","01631031000000d0"); // EG2
  } else if (trainConfig == 228){ // variations for JJ MC (all, MB only and special headers)
    cuts.AddCutCalo("80010103","1111100017032230000","01631031000000d0"); // all headers
    cuts.AddCutCalo("80010113","1111100017032230000","01631031000000d0"); // MB only
    cuts.AddCutCalo("80010123","1111100017032230000","01631031000000d0"); // special header
  } else if (trainConfig == 229){ // variations when using special JJ case (Jets only)
    cuts.AddCutCalo("80010123","1111100017032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("80083123","1111100017032230000","01631031000000d0"); // EG1
    cuts.AddCutCalo("80085123","1111100017032230000","01631031000000d0"); // EG2
  //-----------------------------------------------------------------------------------------------
  // Systematics variations MB run 2 EMC pPb std cut: cuts.AddCutCalo("80010113","1111141057032230000","01631031000000d0"); // 0-100
  //-----------------------------------------------------------------------------------------------
  } else if (trainConfig == 230){ // set of variations CLUSTER
    cuts.AddCutCalo("80010113","1111141057012230000","01631031000000d0"); // min energy cluster variation 1  500 MeV
    cuts.AddCutCalo("80010113","1111141057022230000","01631031000000d0"); // min energy cluster variation 1  600 MeV
    cuts.AddCutCalo("80010113","1111141057042230000","01631031000000d0"); // min energy cluster variation 2  800 MeV
    cuts.AddCutCalo("80010113","1111141057052230000","01631031000000d0"); // min energy cluster variation 3  900 MeV
  } else if (trainConfig == 231){ // set of variations CLUSTER
    cuts.AddCutCalo("80010113","1111141057032200000","01631031000000d0"); // min/max M02  0.1<M<100
    cuts.AddCutCalo("80010113","1111141057032210000","01631031000000d0"); // min/max M02  0.1<M<1
    cuts.AddCutCalo("80010113","1111141057032220000","01631031000000d0"); // min/max M02  0.1<M<0.7
    cuts.AddCutCalo("80010113","1111141057032240000","01631031000000d0"); // min/max M02  0.1<M<0.4
    cuts.AddCutCalo("80010113","1111141057032250000","01631031000000d0"); // min/max M02  0.1<M<0.3
    cuts.AddCutCalo("80010113","1111141057032260000","01631031000000d0"); // min/max M02  0.1<M<0.27
  } else if (trainConfig == 232){ // set of variations CLUSTER
    cuts.AddCutCalo("80010113","1111141057031230000","01631031000000d0"); // min number of cells variation 1  1 cell
    cuts.AddCutCalo("80010113","1112141057032230000","01631031000000d0"); // only modules with TRD infront
    cuts.AddCutCalo("80010113","1111341057032230000","01631031000000d0"); // no modules with TRD infront
  } else if (trainConfig == 233){ // set of variations MESON
    cuts.AddCutCalo("80010113","1111141057032230000","01633031000000d0"); // rapidity variation  y<0.6
    cuts.AddCutCalo("80010113","1111141057032230000","01634031000000d0"); // rapidity variation  y<0.5
    cuts.AddCutCalo("80010113","1111141057032230000","01631061000000d0"); // alpha meson variation 1   0<alpha<0.8
    cuts.AddCutCalo("80010113","1111141057032230000","01631051000000d0"); // alpha meson variation 2  0<alpha<0.75
  } else if (trainConfig == 234){ // opening angle variations
    cuts.AddCutCalo("80010113","1111141057032230000","0163103100000040"); // min opening angle 0.0152
    cuts.AddCutCalo("80010113","1111141057032230000","0163103100000060"); // min opening angle 0.017
    cuts.AddCutCalo("80010113","1111141057032230000","0163103100000070"); // min opening angle 0.016
    cuts.AddCutCalo("80010113","1111141057032230000","0163103100000080"); // min opening angle 0.018
    cuts.AddCutCalo("80010113","1111141057032230000","0163103100000090"); // min opening angle 0.018
  } else if (trainConfig == 235){ // TM variations
    cuts.AddCutCalo("80010113","1111141053032230000","01631031000000d0"); // fixed window
    cuts.AddCutCalo("80010113","1111141056032230000","01631031000000d0"); // tm pt dependent var 1
    cuts.AddCutCalo("80010113","1111141058032230000","01631031000000d0"); // tm pt dependent var 2
    cuts.AddCutCalo("80010113","1111141059032230000","01631031000000d0"); // tm pt dependent var 3
  } else if (trainConfig == 236){ // EMCAL clusters standard cuts, non lin variations
    cuts.AddCutCalo("80010113","1111100057032230000","01631031000000d0"); // 0-100
    cuts.AddCutCalo("80010113","1111142057032230000","01631031000000d0"); // 0-100
    cuts.AddCutCalo("80010113","1111151057032230000","01631031000000d0"); // 0-100
    cuts.AddCutCalo("80010113","1111152057032230000","01631031000000d0"); // 0-100

  //-----------------------------------------------------------------------------------------------
  // Systematics variations 0-10% run 2 EMC pPb std cut: cuts.AddCutCalo("80010113","1111141057032230000","01631031000000d0"); // 0-10
  //-----------------------------------------------------------------------------------------------
  } else if (trainConfig == 240){ // set of variations CLUSTER
    cuts.AddCutCalo("80110113","1111141057012230000","01631031000000d0"); // min energy cluster variation 1  500 MeV
    cuts.AddCutCalo("80110113","1111141057022230000","01631031000000d0"); // min energy cluster variation 1  600 MeV
    cuts.AddCutCalo("80110113","1111141057042230000","01631031000000d0"); // min energy cluster variation 2  800 MeV
    cuts.AddCutCalo("80110113","1111141057052230000","01631031000000d0"); // min energy cluster variation 3  900 MeV
  } else if (trainConfig == 241){ // set of variations CLUSTER
    cuts.AddCutCalo("80110113","1111141057032200000","01631031000000d0"); // min/max M02  0.1<M<100
    cuts.AddCutCalo("80110113","1111141057032210000","01631031000000d0"); // min/max M02  0.1<M<1
    cuts.AddCutCalo("80110113","1111141057032220000","01631031000000d0"); // min/max M02  0.1<M<0.7
    cuts.AddCutCalo("80110113","1111141057032240000","01631031000000d0"); // min/max M02  0.1<M<0.4
    cuts.AddCutCalo("80110113","1111141057032250000","01631031000000d0"); // min/max M02  0.1<M<0.3
    cuts.AddCutCalo("80110113","1111141057032260000","01631031000000d0"); // min/max M02  0.1<M<0.27
  } else if (trainConfig == 242){ // set of variations CLUSTER
    cuts.AddCutCalo("80110113","1111141057031230000","01631031000000d0"); // min number of cells variation 1  1 cell
    cuts.AddCutCalo("80110113","1112141057032230000","01631031000000d0"); // only modules with TRD infront
    cuts.AddCutCalo("80110113","1111341057032230000","01631031000000d0"); // no modules with TRD infront
  } else if (trainConfig == 243){ // set of variations MESON
    cuts.AddCutCalo("80110113","1111141057032230000","01633031000000d0"); // rapidity variation  y<0.6
    cuts.AddCutCalo("80110113","1111141057032230000","01634031000000d0"); // rapidity variation  y<0.5
    cuts.AddCutCalo("80110113","1111141057032230000","01631061000000d0"); // alpha meson variation 1   0<alpha<0.8
    cuts.AddCutCalo("80110113","1111141057032230000","01631051000000d0"); // alpha meson variation 2  0<alpha<0.75
  } else if (trainConfig == 244){ // opening angle variations
    cuts.AddCutCalo("80110113","1111141057032230000","0163103100000040"); // min opening angle 0.0152
    cuts.AddCutCalo("80110113","1111141057032230000","0163103100000060"); // min opening angle 0.017
    cuts.AddCutCalo("80110113","1111141057032230000","0163103100000070"); // min opening angle 0.016
    cuts.AddCutCalo("80110113","1111141057032230000","0163103100000080"); // min opening angle 0.018
    cuts.AddCutCalo("80110113","1111141057032230000","0163103100000090"); // min opening angle 0.018
  } else if (trainConfig == 245){ // TM variations
    cuts.AddCutCalo("80110113","1111141053032230000","01631031000000d0"); // fixed window
    cuts.AddCutCalo("80110113","1111141056032230000","01631031000000d0"); // tm pt dependent var 1
    cuts.AddCutCalo("80110113","1111141058032230000","01631031000000d0"); // tm pt dependent var 2
    cuts.AddCutCalo("80110113","1111141059032230000","01631031000000d0"); // tm pt dependent var 3
  } else if (trainConfig == 246){ // EMCAL clusters standard cuts, non lin variations
    cuts.AddCutCalo("80110113","1111142057032230000","01631031000000d0"); // 0-100
    cuts.AddCutCalo("80110113","1111151057032230000","01631031000000d0"); // 0-100
    cuts.AddCutCalo("80110113","1111152057032230000","01631031000000d0"); // 0-100
    cuts.AddCutCalo("80110113","1111155057032230000","01631031000000d0"); // 0-100
    cuts.AddCutCalo("80110113","1111156057032230000","01631031000000d0"); // 0-100

  //-----------------------------------------------------------------------------------------------
  // Systematics variations 80-100% run 2 EMC pPb std cut: cuts.AddCutCalo("80010113","1111141057032230000","01631031000000d0"); // 80-100
  //-----------------------------------------------------------------------------------------------
  } else if (trainConfig == 250){ // set of variations CLUSTER
    cuts.AddCutCalo("88010113","1111141057012230000","01631031000000d0"); // min energy cluster variation 1  500 MeV
    cuts.AddCutCalo("88010113","1111141057022230000","01631031000000d0"); // min energy cluster variation 1  600 MeV
    cuts.AddCutCalo("88010113","1111141057042230000","01631031000000d0"); // min energy cluster variation 2  800 MeV
    cuts.AddCutCalo("88010113","1111141057052230000","01631031000000d0"); // min energy cluster variation 3  900 MeV
  } else if (trainConfig == 251){ // set of variations CLUSTER
    cuts.AddCutCalo("88010113","1111141057032200000","01631031000000d0"); // min/max M02  0.1<M<100
    cuts.AddCutCalo("88010113","1111141057032210000","01631031000000d0"); // min/max M02  0.1<M<1
    cuts.AddCutCalo("88010113","1111141057032220000","01631031000000d0"); // min/max M02  0.1<M<0.7
    cuts.AddCutCalo("88010113","1111141057032240000","01631031000000d0"); // min/max M02  0.1<M<0.4
    cuts.AddCutCalo("88010113","1111141057032250000","01631031000000d0"); // min/max M02  0.1<M<0.3
    cuts.AddCutCalo("88010113","1111141057032260000","01631031000000d0"); // min/max M02  0.1<M<0.27
  } else if (trainConfig == 252){ // set of variations CLUSTER
    cuts.AddCutCalo("88010113","1111141057031230000","01631031000000d0"); // min number of cells variation 1  1 cell
    cuts.AddCutCalo("88010113","1112141057032230000","01631031000000d0"); // only modules with TRD infront
    cuts.AddCutCalo("88010113","1111341057032230000","01631031000000d0"); // no modules with TRD infront
  } else if (trainConfig == 253){ // set of variations MESON
    cuts.AddCutCalo("88010113","1111141057032230000","01633031000000d0"); // rapidity variation  y<0.6
    cuts.AddCutCalo("88010113","1111141057032230000","01634031000000d0"); // rapidity variation  y<0.5
    cuts.AddCutCalo("88010113","1111141057032230000","01631061000000d0"); // alpha meson variation 1   0<alpha<0.8
    cuts.AddCutCalo("88010113","1111141057032230000","01631051000000d0"); // alpha meson variation 2  0<alpha<0.75
  } else if (trainConfig == 254){ // opening angle variations
    cuts.AddCutCalo("88010113","1111141057032230000","0163103100000040"); // min opening angle 0.0152
    cuts.AddCutCalo("88010113","1111141057032230000","0163103100000060"); // min opening angle 0.017
    cuts.AddCutCalo("88010113","1111141057032230000","0163103100000070"); // min opening angle 0.016
    cuts.AddCutCalo("88010113","1111141057032230000","0163103100000080"); // min opening angle 0.018
    cuts.AddCutCalo("88010113","1111141057032230000","0163103100000090"); // min opening angle 0.018
  } else if (trainConfig == 255){ // TM variations
    cuts.AddCutCalo("88010113","1111141053032230000","01631031000000d0"); // fixed window
    cuts.AddCutCalo("88010113","1111141056032230000","01631031000000d0"); // tm pt dependent var 1
    cuts.AddCutCalo("88010113","1111141058032230000","01631031000000d0"); // tm pt dependent var 2
    cuts.AddCutCalo("88010113","1111141059032230000","01631031000000d0"); // tm pt dependent var 3
  } else if (trainConfig == 256){ // EMCAL clusters standard cuts, non lin variations
    cuts.AddCutCalo("88010113","1111142057032230000","01631031000000d0"); // 0-100
    cuts.AddCutCalo("88010113","1111151057032230000","01631031000000d0"); // 0-100
    cuts.AddCutCalo("88010113","1111152057032230000","01631031000000d0"); // 0-100
    cuts.AddCutCalo("88010113","1111155057032230000","01631031000000d0"); // 0-100
    cuts.AddCutCalo("88010113","1111156057032230000","01631031000000d0"); // 0-100

  // EMC, DCal, EMC+DCal calibs
  } else if (trainConfig == 260){
    cuts.AddCutCalo("80010113","1111100057032230000","01631031000000d0"); // EMC
    cuts.AddCutCalo("80010113","3885500057032230000","01631031000000d0"); // DCal
    cuts.AddCutCalo("80010113","4117900057032230000","01631031000000d0"); // EDC
  } else if (trainConfig == 261){
    cuts.AddCutCalo("80010113","4117900057032230000","01631031000000d0"); // EDC
  } else if (trainConfig == 262){
    cuts.AddCutCalo("80010113","4117900007032230000","01631031000000d0"); // EDC
    cuts.AddCutCalo("80052113","4117900007032230000","01631031000000d0"); // EDC
    cuts.AddCutCalo("80085113","4117900007032230000","01631031000000d0"); // EDC
    cuts.AddCutCalo("80081113","4117900007032230000","01631031000000d0"); // EDC

  } else if (trainConfig == 263){
    cuts.AddCutCalo("80010113","4117941057032230000","01631031000000d0"); // 0-100
  } else if (trainConfig == 264){ // EMCAL clusters standard cuts, cent vars, non lin variations
    cuts.AddCutCalo("80110113","4117941057032230000","01631031000000d0"); // 0-10
    cuts.AddCutCalo("81210113","4117941057032230000","01631031000000d0"); // 10-20
    cuts.AddCutCalo("82410113","4117941057032230000","01631031000000d0"); // 20-40
    cuts.AddCutCalo("84610113","4117941057032230000","01631031000000d0"); // 40-60
    cuts.AddCutCalo("86810113","4117941057032230000","01631031000000d0"); // 60-80
    cuts.AddCutCalo("88010113","4117941057032230000","01631031000000d0"); // 80-100
  } else if (trainConfig == 265){
    cuts.AddCutCalo("80010113","4117900057032230000","01631031000000d0"); // 0-100
  // ===============================================================================================
  // Run 1 data PHOS clusters pPb 5TeV
  // ===============================================================================================
  } else if (trainConfig == 301) {  // min energy = 0.3 GeV/c
    cuts.AddCutCalo("80010113","2444400041013200000","0163103100000010"); //standart cut, kINT7 // PHOS clusters
    cuts.AddCutCalo("80062113","2444400041013200000","0163103100000010"); //standard cut, kPHI7  // PHOS clusters
  } else if (trainConfig == 302){ // Validation PHOS
    cuts.AddCutCalo("80010113","2444400041013200000","0163103100000010");
  } else if (trainConfig == 303){ // Validation PHOS, only added signals
    cuts.AddCutCalo("80010023","2444400041013200000","0163103100000010");
  } else if (trainConfig == 304){ // min energy = 0.3 GeV/c
    cuts.AddCutCalo("80010113","2444400041013200000","0163103100000000"); // kINT7 // PHOS clusters no open angle cut
    cuts.AddCutCalo("80062113","2444400041013200000","0163103100000000"); // kPHI7 // PHOS clusters no open angle cut
    cuts.AddCutCalo("80010113","2444400041013200000","0163103100000030"); // kINT7 // PHOS clusters open angle cut 0.1
    cuts.AddCutCalo("80062113","2444400041013200000","0163103100000030"); // kPHI7 // PHOS clusters  open angle cut 0.1
  } else if (trainConfig == 305){ // timing cut variations
    cuts.AddCutCalo("80010113","2444400011013200000","0163103100000010"); // 1000ns
    cuts.AddCutCalo("80010113","2444400031013200000","0163103100000010"); // 200ns
    cuts.AddCutCalo("80010113","2444400051013200000","0163103100000010"); // 50ns
  } else if (trainConfig == 306) {  // PHOS non lin var INT7
    cuts.AddCutCalo("80010113","2444401041013200000","0163103100000010"); // PHOS group standard
    cuts.AddCutCalo("80010113","2444451041013200000","0163103100000010"); // CCMF PHOS
    cuts.AddCutCalo("80010113","2444452041013200000","0163103100000010"); // CMF PHOS
  } else if (trainConfig == 307) {  // PHOS non lin var PHI7
    cuts.AddCutCalo("80062113","2444401041013200000","0163103100000010"); // PHOS group standard
    cuts.AddCutCalo("80062113","2444451041013200000","0163103100000010"); // CCMF PHOS
    cuts.AddCutCalo("80062113","2444452041013200000","0163103100000010"); // CMF PHOS
  } else if (trainConfig == 308) {  // PHOS CCMF cent dep
    cuts.AddCutCalo("80210113","2444451041013200000","0163103100000010"); // 0-20
    cuts.AddCutCalo("82410113","2444451041013200000","0163103100000010"); // 20-40
    cuts.AddCutCalo("84610113","2444451041013200000","0163103100000010"); // 40-60
    cuts.AddCutCalo("86010113","2444451041013200000","0163103100000010"); // 60-100
  } else if (trainConfig == 309) {  // PHOS default cent dep
    cuts.AddCutCalo("80210113","2444401041013200000","0163103100000010"); // 0-20
    cuts.AddCutCalo("82410113","2444401041013200000","0163103100000010"); // 20-40
    cuts.AddCutCalo("84610113","2444401041013200000","0163103100000010"); // 40-60
    cuts.AddCutCalo("86010113","2444401041013200000","0163103100000010"); // 60-100

  } else if(trainConfig == 310){ // first set of variations CLUSTER
    cuts.AddCutCalo("80010113","2444451041073200000","0163103100000010"); // min energy 0.2 GeV
    cuts.AddCutCalo("80010113","2444451041083200000","0163103100000010"); // min energy 0.4 GeV
    cuts.AddCutCalo("80010113","2444451041023200000","0163103100000010"); // min energy 0.5 GeV
  } else if(trainConfig == 311){ // second set of variations CLUSTER
    cuts.AddCutCalo("80010113","2444451041013270000","0163103100000010"); // min/max M02  0.1<M<1
    cuts.AddCutCalo("80010113","2444451041013280000","0163103100000010"); // min/max M02  0.1<M<1
    cuts.AddCutCalo("80010113","2444451041012200000","0163103100000010"); // min number 2 cells
    cuts.AddCutCalo("80010113","2444451041014200000","0163103100000010"); // min number 4 cells
  } else if(trainConfig == 312){ // MESON
    cuts.AddCutCalo("80010113","2444451041013200000","0163403100000010"); // rapidity variation  y<0.5
    cuts.AddCutCalo("80010113","2444451041013200000","0163803100000010"); // rapidity variation  y<0.25
    cuts.AddCutCalo("80010113","2444451041013200000","0163106100000010"); // alpha meson variation 1 0<alpha<0.8
    cuts.AddCutCalo("80010113","2444451041013200000","0163105100000010"); // alpha meson variation 2 0<alpha<0.75
  } else if(trainConfig == 313){ // fourth set of variations
    cuts.AddCutCalo("80010113","2444451044013200000","0163103100000010"); // tm variation
    cuts.AddCutCalo("80010113","2444451045013200000","0163103100000010"); // tm variation
    cuts.AddCutCalo("80010113","2444451046013200000","0163103100000010"); // tm variation
    cuts.AddCutCalo("80010113","2444451041013200000","0163103100000000"); // min opening angle 0    -> open
    cuts.AddCutCalo("80010113","2444451041013200000","0163103100000030"); // min opening angle 0.01 -> 2 cell diag

  } else if(trainConfig == 320){ // reproducing Dmitri's results pi0, eta
    cuts.AddCutCalo("80010113","2444400040013300000","0163103100000010"); // dmitri default pi0/eta w/ opening angle
    cuts.AddCutCalo("80010113","2444400040013300000","0163103100000000"); // dmitri default pi0/eta w/o opening angle
  } else if(trainConfig == 321){ // reproducing Dmitri's results pi0, eta
    cuts.AddCutCalo("80210113","2444400040013300000","0163103100000010"); // dmitri default pi0/eta w/ opening angle
    cuts.AddCutCalo("82410113","2444400040013300000","0163103100000010"); // dmitri default pi0/eta w/ opening angle
    cuts.AddCutCalo("84610113","2444400040013300000","0163103100000010"); // dmitri default pi0/eta w/ opening angle
    cuts.AddCutCalo("86010113","2444400040013300000","0163103100000010"); // dmitri default pi0/eta w/ opening angle
  } else if(trainConfig == 322){ // reproducing Dmitri's results pi0, eta
    cuts.AddCutCalo("80210113","2444400040013300000","0163103100000000"); // dmitri default pi0/eta w/o opening angle
    cuts.AddCutCalo("82410113","2444400040013300000","0163103100000000"); // dmitri default pi0/eta w/o opening angle
    cuts.AddCutCalo("84610113","2444400040013300000","0163103100000000"); // dmitri default pi0/eta w/o opening angle
    cuts.AddCutCalo("86010113","2444400040013300000","0163103100000000"); // dmitri default pi0/eta w/o opening angle
  } else if(trainConfig == 323){ // reproducing Dmitri's results pi0, eta, PHI7 trigg
    cuts.AddCutCalo("80062113","2444400040013300000","0163103100000010"); // dmitri default pi0/eta w/ opening angle
    cuts.AddCutCalo("80062113","2444400040013300000","0163103100000000"); // dmitri default pi0/eta w/ opening angle
  } else if(trainConfig == 324){ // reproducing Dmitri's results pi0, eta, PHI7 trigg
    cuts.AddCutCalo("80262113","2444400040013300000","0163103100000010"); // dmitri default pi0/eta w/ opening angle
    cuts.AddCutCalo("82462113","2444400040013300000","0163103100000010"); // dmitri default pi0/eta w/ opening angle
    cuts.AddCutCalo("84662113","2444400040013300000","0163103100000010"); // dmitri default pi0/eta w/ opening angle
    cuts.AddCutCalo("86062113","2444400040013300000","0163103100000010"); // dmitri default pi0/eta w/ opening angle
  } else if(trainConfig == 325){ // reproducing Dmitri's results pi0, eta, PHI7 trigg
    cuts.AddCutCalo("80262113","2444400040013300000","0163103100000000"); // dmitri default pi0/eta w/o opening angle
    cuts.AddCutCalo("82462113","2444400040013300000","0163103100000000"); // dmitri default pi0/eta w/o opening angle
    cuts.AddCutCalo("84662113","2444400040013300000","0163103100000000"); // dmitri default pi0/eta w/o opening angle
    cuts.AddCutCalo("86062113","2444400040013300000","0163103100000000"); // dmitri default pi0/eta w/o opening angle
  } else if(trainConfig == 326){ // reproducing Dmitri's results pi0, eta
    cuts.AddCutCalo("80010113","2444451040013300000","0163103100000010"); // dmitri default pi0/eta w/ opening angle
    cuts.AddCutCalo("80010113","2444451040013300000","0163103100000000"); // dmitri default pi0/eta w/o opening angle
  } else if(trainConfig == 327){ // reproducing Dmitri's results pi0, eta
    cuts.AddCutCalo("80210113","2444451040013300000","0163103100000010"); // dmitri default pi0/eta w/ opening angle PCM-PHOS NL
    cuts.AddCutCalo("82410113","2444451040013300000","0163103100000010"); // dmitri default pi0/eta w/ opening angle PCM-PHOS NL
    cuts.AddCutCalo("84610113","2444451040013300000","0163103100000010"); // dmitri default pi0/eta w/ opening angle PCM-PHOS NL
    cuts.AddCutCalo("86010113","2444451040013300000","0163103100000010"); // dmitri default pi0/eta w/ opening angle PCM-PHOS NL
  } else if(trainConfig == 328){ // reproducing Dmitri's results pi0, eta, PHI7 trigg
    cuts.AddCutCalo("80062113","2444451040013300000","0163103100000010"); // dmitri default pi0/eta w/ opening angle PCM-PHOS NL
    cuts.AddCutCalo("80062113","2444451040013300000","0163103100000000"); // dmitri default pi0/eta w/ opening angle PCM-PHOS NL
  } else if(trainConfig == 329){ // reproducing Dmitri's results pi0, eta, PHI7 trigg
    cuts.AddCutCalo("80262113","2444451040013300000","0163103100000010"); // dmitri default pi0/eta w/ opening angle PCM-PHOS NL
    cuts.AddCutCalo("82462113","2444451040013300000","0163103100000010"); // dmitri default pi0/eta w/ opening angle PCM-PHOS NL
    cuts.AddCutCalo("84662113","2444451040013300000","0163103100000010"); // dmitri default pi0/eta w/ opening angle PCM-PHOS NL
    cuts.AddCutCalo("86062113","2444451040013300000","0163103100000010"); // dmitri default pi0/eta w/ opening angle PCM-PHOS NL

  } else if(trainConfig == 340){ // reproducing Dmitri's results gamma
    cuts.AddCutCalo("80010113","2444400048013300020","0163103100000010"); // dmitri default gamma w/ opening angle
    cuts.AddCutCalo("80010113","2444400048013300020","0163103100000000"); // dmitri default gamma w/o opening angle
  } else if(trainConfig == 341){ // reproducing Dmitri's results gamma
    cuts.AddCutCalo("80210113","2444400048013300020","0163103100000010"); // dmitri default gamma w/ opening angle
    cuts.AddCutCalo("82410113","2444400048013300020","0163103100000010"); // dmitri default gamma w/ opening angle
    cuts.AddCutCalo("84610113","2444400048013300020","0163103100000010"); // dmitri default gamma w/ opening angle
    cuts.AddCutCalo("86010113","2444400048013300020","0163103100000010"); // dmitri default gamma w/ opening angle
  } else if(trainConfig == 342){ // reproducing Dmitri's results gamma
    cuts.AddCutCalo("80210113","2444400048013300020","0163103100000000"); // dmitri default gamma w/o opening angle
    cuts.AddCutCalo("82410113","2444400048013300020","0163103100000000"); // dmitri default gamma w/o opening angle
    cuts.AddCutCalo("84610113","2444400048013300020","0163103100000000"); // dmitri default gamma w/o opening angle
    cuts.AddCutCalo("86010113","2444400048013300020","0163103100000000"); // dmitri default gamma w/o opening angle
  } else if(trainConfig == 343){ // reproducing Dmitri's results gamma PHI7 trigg
    cuts.AddCutCalo("80062113","2444400048013300020","0163103100000010"); // dmitri default gamma w/ opening angle
    cuts.AddCutCalo("80062113","2444400048013300020","0163103100000000"); // dmitri default gamma w/o opening angle
  } else if(trainConfig == 344){ // reproducing Dmitri's results gamma PHI7 trigg
    cuts.AddCutCalo("80262113","2444400048013300020","0163103100000010"); // dmitri default gamma w/ opening angle
    cuts.AddCutCalo("82462113","2444400048013300020","0163103100000010"); // dmitri default gamma w/ opening angle
    cuts.AddCutCalo("84662113","2444400048013300020","0163103100000010"); // dmitri default gamma w/ opening angle
    cuts.AddCutCalo("86062113","2444400048013300020","0163103100000010"); // dmitri default gamma w/ opening angle
  } else if(trainConfig == 345){ // reproducing Dmitri's results gamma PHI7 trigg
    cuts.AddCutCalo("80262113","2444400048013300020","0163103100000000"); // dmitri default gamma w/o opening angle
    cuts.AddCutCalo("82462113","2444400048013300020","0163103100000000"); // dmitri default gamma w/o opening angle
    cuts.AddCutCalo("84662113","2444400048013300020","0163103100000000"); // dmitri default gamma w/o opening angle
    cuts.AddCutCalo("86062113","2444400048013300020","0163103100000000"); // dmitri default gamma w/o opening angle
  } else if(trainConfig == 346){ // reproducing Dmitri's results pi0, eta
    cuts.AddCutCalo("80010113","2444451048013300020","0163103100000010"); // dmitri default gamma w/ opening angle
    cuts.AddCutCalo("80010113","2444451048013300020","0163103100000000"); // dmitri default gamma w/o opening angle
  } else if(trainConfig == 347){ // reproducing Dmitri's results pi0, eta
    cuts.AddCutCalo("80210113","2444451048013300020","0163103100000010"); // dmitri default gamma w/ opening angle PCM-PHOS NL
    cuts.AddCutCalo("82410113","2444451048013300020","0163103100000010"); // dmitri default gamma w/ opening angle PCM-PHOS NL
    cuts.AddCutCalo("84610113","2444451048013300020","0163103100000010"); // dmitri default gamma w/ opening angle PCM-PHOS NL
    cuts.AddCutCalo("86010113","2444451048013300020","0163103100000010"); // dmitri default gamma w/ opening angle PCM-PHOS NL
  } else if(trainConfig == 348){ // reproducing Dmitri's results pi0, eta, PHI7 trigg
    cuts.AddCutCalo("80062113","2444451048013300020","0163103100000010"); // dmitri default gamma w/ opening angle PCM-PHOS NL
    cuts.AddCutCalo("80062113","2444451048013300020","0163103100000000"); // dmitri default gamma w/ opening angle PCM-PHOS NL
  } else if(trainConfig == 349){ // reproducing Dmitri's results pi0, eta, PHI7 trigg
    cuts.AddCutCalo("80262113","2444451048013300020","0163103100000010"); // dmitri default gamma w/ opening angle PCM-PHOS NL
    cuts.AddCutCalo("82462113","2444451048013300020","0163103100000010"); // dmitri default gamma w/ opening angle PCM-PHOS NL
    cuts.AddCutCalo("84662113","2444451048013300020","0163103100000010"); // dmitri default gamma w/ opening angle PCM-PHOS NL
    cuts.AddCutCalo("86062113","2444451048013300020","0163103100000010"); // dmitri default gamma w/ opening angle PCM-PHOS NL

  // ===============================================================================================
  // Run 2 data EMC clusters pPb 8TeV
  // ===============================================================================================
  // pPb 8TeV variations for QA
  } else if (trainConfig == 400){ // EMCAL clusters standard cuts, triggers, no nonlin, open timing (same as 201 but split)
    cuts.AddCutCalo("80010113","1111100017032230000","01631031000000d0"); // INT7
  } else if (trainConfig == 401){ // EMCAL clusters standard cuts, triggers, no nonlin, open timing (same as 201 but split)
    cuts.AddCutCalo("80085113","1111100017032230000","01631031000000d0"); // EG2
    cuts.AddCutCalo("80083113","1111100017032230000","01631031000000d0"); // EG1
  } else if (trainConfig == 402){ // EMCAL clusters standard cuts, triggers, no nonlin, open timing (same as 201 but split)
    cuts.AddCutCalo("80095113","1111100017032230000","01631031000000d0"); // EJ2
    cuts.AddCutCalo("80093113","1111100017032230000","01631031000000d0"); // EJ1

  } else if (trainConfig == 404){ // EMCAL clusters standard cuts, triggers, no nonlin, +-50ns timing (5)
    cuts.AddCutCalo("80010113","1111100057032230000","01631031000000d0"); // INT7
  } else if (trainConfig == 405){ // EMCAL clusters standard cuts, triggers, no nonlin, +-50ns timing (5)
    cuts.AddCutCalo("80085113","1111100057032230000","01631031000000d0"); // EG2
    cuts.AddCutCalo("80083113","1111100057032230000","01631031000000d0"); // EG1
  } else if (trainConfig == 406){ // EMCAL clusters standard cuts, triggers, no nonlin, +-50ns timing (5)
    cuts.AddCutCalo("80095113","1111100057032230000","01631031000000d0"); // EJ2
    cuts.AddCutCalo("80093113","1111100057032230000","01631031000000d0"); // EJ1

  } else if (trainConfig == 407){ // EMCAL clusters standard cuts, triggers, no nonlin, +-50ns timing (5), Nico TB NL
    cuts.AddCutCalo("80010113","1111165057032230000","01631031000000d0"); // INT7
  } else if (trainConfig == 408){ // EMCAL clusters standard cuts, triggers, no nonlin, +-50ns timing (5), Nico TB NL
    cuts.AddCutCalo("80085113","1111165057032230000","01631031000000d0"); // EG2
    cuts.AddCutCalo("80083113","1111165057032230000","01631031000000d0"); // EG1
  } else if (trainConfig == 409){ // EMCAL clusters standard cuts, triggers, no nonlin, +-50ns timing (5), Nico TB NL
    cuts.AddCutCalo("80095113","1111165057032230000","01631031000000d0"); // EJ2
    cuts.AddCutCalo("80093113","1111165057032230000","01631031000000d0"); // EJ1

  } else if (trainConfig == 410){ // EMCAL clusters standard cuts, triggers, NL vars
    cuts.AddCutCalo("80010113","1111147057032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("80010113","1111148057032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("80010113","1111157057032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("80010113","1111158057032230000","01631031000000d0"); // INT7
  } else if (trainConfig == 411){ // EMCAL clusters standard cuts, triggers, NL vars
    cuts.AddCutCalo("80085113","1111147057032230000","01631031000000d0"); // EG2
    cuts.AddCutCalo("80085113","1111148057032230000","01631031000000d0"); // EG2
    cuts.AddCutCalo("80085113","1111157057032230000","01631031000000d0"); // EG2
    cuts.AddCutCalo("80085113","1111158057032230000","01631031000000d0"); // EG2
  } else if (trainConfig == 412){ // EMCAL clusters standard cuts, triggers, NL vars
    cuts.AddCutCalo("80083113","1111147057032230000","01631031000000d0"); // EG1
    cuts.AddCutCalo("80083113","1111148057032230000","01631031000000d0"); // EG1
    cuts.AddCutCalo("80083113","1111157057032230000","01631031000000d0"); // EG1
    cuts.AddCutCalo("80083113","1111158057032230000","01631031000000d0"); // EG1

  } else if (trainConfig == 413){ // EMCAL clusters standard cuts, triggers, NL vars, TM E/p 1.75 (f)
    cuts.AddCutCalo("80010113","111114705f032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("80010113","111114805f032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("80010113","111115705f032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("80010113","111115805f032230000","01631031000000d0"); // INT7
  } else if (trainConfig == 414){ // EMCAL clusters standard cuts, triggers, NL vars, TM E/p 1.75 (f)
    cuts.AddCutCalo("80085113","111114705f032230000","01631031000000d0"); // EG2
    cuts.AddCutCalo("80085113","111114805f032230000","01631031000000d0"); // EG2
    cuts.AddCutCalo("80085113","111115705f032230000","01631031000000d0"); // EG2
    cuts.AddCutCalo("80085113","111115805f032230000","01631031000000d0"); // EG2
  } else if (trainConfig == 415){ // EMCAL clusters standard cuts, triggers, NL vars, TM E/p 1.75 (f)
    cuts.AddCutCalo("80083113","111114705f032230000","01631031000000d0"); // EG1
    cuts.AddCutCalo("80083113","111114805f032230000","01631031000000d0"); // EG1
    cuts.AddCutCalo("80083113","111115705f032230000","01631031000000d0"); // EG1
    cuts.AddCutCalo("80083113","111115805f032230000","01631031000000d0"); // EG1

  } else if (trainConfig == 419){ // EMCAL clusters standard cuts, triggers, TB NL
    cuts.AddCutCalo("80010113","1111106057032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("80085113","1111106057032230000","01631031000000d0"); // EG2
    cuts.AddCutCalo("80083113","1111106057032230000","01631031000000d0"); // EG1

  } else if (trainConfig == 420){ // nonlinearity variations
    cuts.AddCutCalo("80010113","1111141057032230000","01631031000000d0"); // CRF
    cuts.AddCutCalo("80010113","1111142057032230000","01631031000000d0"); // CRF
    cuts.AddCutCalo("80010113","1111151057032230000","01631031000000d0"); // CCMF
    cuts.AddCutCalo("80010113","1111152057032230000","01631031000000d0"); // CMF
  } else if (trainConfig == 421){ // second set of variations CLUSTER
    cuts.AddCutCalo("80010113","1111142057022230000","01631031000000d0"); // min energy cluster variation 1  600 MeV
    cuts.AddCutCalo("80010113","1111142057042230000","01631031000000d0"); // min energy cluster variation 2  800 MeV
    cuts.AddCutCalo("80010113","1111142057052230000","01631031000000d0"); // min energy cluster variation 3  900 MeV
    cuts.AddCutCalo("80010113","1111142057032230000","01631031000000d0"); // min/max M02  0.1<M<0.5
  } else if (trainConfig == 422){ // third set of variations CLUSTER
    cuts.AddCutCalo("80010113","1111142057032250000","01631031000000d0"); // min/max M02  0.1<M<0.3
    cuts.AddCutCalo("80010113","1111142057032260000","01631031000000d0"); // min/max M02  0.1<M<0.27
    cuts.AddCutCalo("80010113","1111142057031230000","01631031000000d0"); // min number of cells variation 1  1 cell
    cuts.AddCutCalo("80010113","1112141057032230000","01631031000000d0"); // only modules with TRD infront
    cuts.AddCutCalo("80010113","1111341057032230000","01631031000000d0"); // no modules with TRD infront
  } else if (trainConfig == 423){ // third set of variations MESON
    cuts.AddCutCalo("80010113","1111142057032230000","01633031000000d0"); // rapidity variation  y<0.6
    cuts.AddCutCalo("80010113","1111142057032230000","01634031000000d0"); // rapidity variation  y<0.5
    cuts.AddCutCalo("80010113","1111142057032230000","01631061000000d0"); // alpha meson variation 1   0<alpha<0.8
    cuts.AddCutCalo("80010113","1111142057032230000","01631051000000d0"); // alpha meson variation 2  0<alpha<0.75
  } else if (trainConfig == 424){ // opening angle variations
    cuts.AddCutCalo("80010113","1111142057032230000","0163103100000040"); // min opening angle 0.0152
    cuts.AddCutCalo("80010113","1111142057032230000","0163103100000060"); // min opening angle 0.017
    cuts.AddCutCalo("80010113","1111142057032230000","0163103100000070"); // min opening angle 0.016
    cuts.AddCutCalo("80010113","1111142057032230000","0163103100000080"); // min opening angle 0.018
    cuts.AddCutCalo("80010113","1111142057032230000","0163103100000090"); // min opening angle 0.018
  } else if (trainConfig == 425){ // TM variations
    cuts.AddCutCalo("80010113","1111142053032230000","01631031000000d0"); // fixed window
    cuts.AddCutCalo("80010113","1111142056032230000","01631031000000d0"); // tm pt dependent var 1
    cuts.AddCutCalo("80010113","1111142058032230000","01631031000000d0"); // tm pt dependent var 2
    cuts.AddCutCalo("80010113","1111142059032230000","01631031000000d0"); // tm pt dependent var 3
  } else if (trainConfig == 426){
    cuts.AddCutCalo("80010113","11111420570322d0000","01631031000000d0"); // M02, pt dep with  0.27, 0.0072, 0.4
    cuts.AddCutCalo("80010113","11111420570322e0000","01631031000000d0"); // M02, pt dep with  0.31, 0.0072, 0.5
    cuts.AddCutCalo("80010113","11111420570322f0000","01631031000000d0"); // M02, pt dep with  0.36, 0.0072, 0.7
  } else if (trainConfig == 427){
    cuts.AddCutCalo("80010113","11111420570322m0000","01631031000000d0"); // M02, pt dep with  0.32, 0.0152, 0.5
    cuts.AddCutCalo("80010113","11111420570322g0000","01631031000000d0"); // M02, pt dep with  0.37, 0.0072, 0.7
    cuts.AddCutCalo("80010113","11111420570322h0000","01631031000000d0"); // M02, pT-dep with  0.30, 0.0072, 0.5
    cuts.AddCutCalo("80010113","11111420570322i0000","01631031000000d0"); // M02, pT-dep with  0.35, 0.0072, 0.7

  } else if (trainConfig == 480){ // EMCAL clusters standard cuts, triggers, no nonlin, open timing CALO+CALOFAST readout  - NO TM
    cuts.AddCutCalo("800a0113","1111100010032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("800a1113","1111100010032230000","01631031000000d0"); // EMC7
    cuts.AddCutCalo("800a2113","1111100010032230000","01631031000000d0"); // EG2
    cuts.AddCutCalo("800a3113","1111100010032230000","01631031000000d0"); // EG1
  } else if (trainConfig == 481){ // EMCAL clusters standard cuts, triggers, no nonlin, open timing  - NO TM
    cuts.AddCutCalo("80010113","1111100010032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("80052113","1111100010032230000","01631031000000d0"); // EMC7
    cuts.AddCutCalo("80085113","1111100010032230000","01631031000000d0"); // EG2
    cuts.AddCutCalo("80083113","1111100010032230000","01631031000000d0"); // EG1
  // ===============================================================================================
  // Run 2 data PHOS clusters pPb 5TeV
  // ===============================================================================================
  // INT7 triggers
  } else if (trainConfig == 500) {  // PHOS  INT7
    cuts.AddCutCalo("80010113","24466420ha012200000","0163103100000010"); // standard
    cuts.AddCutCalo("80010113","24466000ha012200000","0163103100000010"); // standard without non-lin
    cuts.AddCutCalo("80010113","244664205a012200000","0163103100000010"); // standard
    cuts.AddCutCalo("80010113","244660005ha012200000","0163103100000010"); // standard without non-lin
  } else if (trainConfig == 501) {  // PHOS  INT7
    cuts.AddCutCalo("80010113","2446600041012200000","0163103100000010"); // no non lin
    cuts.AddCutCalo("80010113","2446600011012200000","0163103100000010"); // no non lin 1000 \mus
    cuts.AddCutCalo("80010113","2446600061012200000","0163103100000010"); // no non lin, -30, 50ns
    cuts.AddCutCalo("80010113","24466000a1012200000","0163103100000010"); // no non lin, -12.5, 13ns
  } else if (trainConfig == 502) {  // PHOS  INT7 non lin vars
    cuts.AddCutCalo("80010113","2446600051012200000","0163103100000010"); //
    cuts.AddCutCalo("80010113","2446641051012200000","0163103100000010"); //
    cuts.AddCutCalo("80010113","2446642051012200000","0163103100000010"); //
    cuts.AddCutCalo("80010113","2446651051012200000","0163103100000010"); //
    cuts.AddCutCalo("80010113","2446652051012200000","0163103100000010"); //
  } else if (trainConfig == 503) {  // PHOS  INT7 with cents
    cuts.AddCutCalo("80010113","24466420ha012200000","0163103100000010"); // non lin 0-100%
    cuts.AddCutCalo("80110113","24466420ha012200000","0163103100000010"); // non lin 0-10%
    cuts.AddCutCalo("81210113","24466420ha012200000","0163103100000010"); // non lin 10-20%
    cuts.AddCutCalo("82410113","24466420ha012200000","0163103100000010"); // non lin 20-40%
    cuts.AddCutCalo("84610113","24466420ha012200000","0163103100000010"); // non lin 40-60%
    cuts.AddCutCalo("86810113","24466420ha012200000","0163103100000010"); // non lin 60-80%
    cuts.AddCutCalo("88010113","24466420ha012200000","0163103100000010"); // non lin 80-100%
  } else if (trainConfig == 504) {  // PHOS  INT7 with cents
    cuts.AddCutCalo("80010113","24466420ha012200000","0163103100000010"); // non lin 0-100%
    cuts.AddCutCalo("80210113","24466420ha012200000","0163103100000010"); // non lin 0-20%
    cuts.AddCutCalo("86010113","24466420ha012200000","0163103100000010"); // non lin 60-100%
    cuts.AddCutCalo("a0110113","24466420ha012200000","0163103100000010"); // non lin 0-5%
    cuts.AddCutCalo("a1210113","24466420ha012200000","0163103100000010"); // non lin 5-10%
  } else if (trainConfig == 505) {  // PHOS  INT7 with cents
    cuts.AddCutCalo("80110113","24466420ha012200000","0163103100000010"); // non lin 0-10%
    cuts.AddCutCalo("81210113","24466420ha012200000","0163103100000010"); // non lin 10-20%
    cuts.AddCutCalo("82410113","24466420ha012200000","0163103100000010"); // non lin 20-40%
    cuts.AddCutCalo("84610113","24466420ha012200000","0163103100000010"); // non lin 40-60%
    cuts.AddCutCalo("86810113","24466420ha012200000","0163103100000010"); // non lin 60-80%
    cuts.AddCutCalo("88010113","24466420ha012200000","0163103100000010"); // non lin 80-100%
  } else if (trainConfig == 506) {  // PHOS  INT7 with cents
    cuts.AddCutCalo("80210113","24466420ha012200000","0163103100000010"); // no non lin 0-10%
    cuts.AddCutCalo("86010113","24466420ha012200000","0163103100000010"); // no non lin 60-100%
    cuts.AddCutCalo("a0110113","24466420ha012200000","0163103100000010"); // no non lin 0-5%
    cuts.AddCutCalo("a1210113","24466420ha012200000","0163103100000010"); // no non lin 5-10%
  } else if (trainConfig == 507){
    cuts.AddCutCalo("80010113","2446641051012200000","0163103100000010"); // TM on
    cuts.AddCutCalo("80010113","2446641050012200000","0163103100000010"); // TM off
  } else if (trainConfig == 508){ // JJ MC AOD validation
    cuts.AddCutCalo("80010123","2446641051012200000","0163103100000010"); // TM on
    cuts.AddCutCalo("80010123","2446641050012200000","0163103100000010"); // TM off
  } else if (trainConfig == 509){ // JJ MC AOD validation PHOS NL
    cuts.AddCutCalo("80010123","2446601051012200000","0163103100000010"); // TM on
    cuts.AddCutCalo("80010123","2446601050012200000","0163103100000010"); // TM off
  } else if (trainConfig == 510){ // JJ MC AOD validation w/o NL
    cuts.AddCutCalo("80010123","2446600051012200000","0163103100000010"); // TM on
    cuts.AddCutCalo("80010123","2446600050012200000","0163103100000010"); // TM off
  } else if (trainConfig == 511){ // cut variation: dist to bad channel, Ncells, min energy
    cuts.AddCutCalo("80010113","2446641051013200000","0163103100000010"); // no dist to bad channel cut, Ncells: 3, min energy 300 MeV
    cuts.AddCutCalo("80010113","2446641051083200000","0163103100000010"); // no dist to bad channel cut, Ncells: 3, min energy 400 MeV
    cuts.AddCutCalo("80010113","2446641051012200000","0163103100000010"); // no dist to bad channel cut, Ncells: 2, min energy 300 MeV
    cuts.AddCutCalo("80010113","2446641051082200000","0163103100000010"); // no dist to bad channel cut, Ncells: 2, min energy 400 MeV
  } else if (trainConfig == 512){ // cut variation: dist to bad channel, Ncells, min energy
    cuts.AddCutCalo("80010113","2446641151013200000","0163103100000010"); // dist to bad channel: 1, Ncells: 3, min energy 300 MeV
    cuts.AddCutCalo("80010113","2446641151083200000","0163103100000010"); // dist to bad channel: 1, Ncells: 3, min energy 400 MeV
    cuts.AddCutCalo("80010113","2446641151012200000","0163103100000010"); // dist to bad channel: 1, Ncells: 2, min energy 300 MeV
    cuts.AddCutCalo("80010113","2446641151082200000","0163103100000010"); // dist to bad channel: 1, Ncells: 2, min energy 400 MeV
  } else if (trainConfig == 513){ // cut variation: dist to bad channel
    cuts.AddCutCalo("80010113","2446642051013200000","0163103100000010"); // dist to bad channel: 0
    cuts.AddCutCalo("80010113","2446642151083200000","0163103100000010"); // dist to bad channel: 1
    cuts.AddCutCalo("80010113","2446642251012200000","0163103100000010"); // dist to bad channel: 2
    cuts.AddCutCalo("80010113","2446642351082200000","0163103100000010"); // dist to bad channel: 3
  } else if (trainConfig == 514){
    cuts.AddCutCalo("80010113","2446642151012200000","0163103100000010"); // standard
    cuts.AddCutCalo("80010113","2446642101012200000","0163103100000010"); // standard without timing cut

  } else if (trainConfig == 520) {  // JJ MC
    cuts.AddCutCalo("80010123","24466420ha012200000","0163103100000010"); // standard
  } else if (trainConfig == 521) {
    cuts.AddCutCalo("80010113","2446600004012200000","0163103100000010"); // no non lin 1000 \mus
    cuts.AddCutCalo("80062113","2446600004012200000","0163103100000010"); // no non lin 1000 \mus
  } else if (trainConfig == 522) {  // JJ MC with cents
    cuts.AddCutCalo("80010123","24466420ha012200000","0163103100000010"); // non lin 0-100%
    cuts.AddCutCalo("80110123","24466420ha012200000","0163103100000010"); // non lin 0-10%
    cuts.AddCutCalo("81210123","24466420ha012200000","0163103100000010"); // non lin 10-20%
    cuts.AddCutCalo("82410123","24466420ha012200000","0163103100000010"); // non lin 20-40%
    cuts.AddCutCalo("84610123","24466420ha012200000","0163103100000010"); // non lin 40-60%
    cuts.AddCutCalo("86810123","24466420ha012200000","0163103100000010"); // non lin 60-80%
    cuts.AddCutCalo("88010123","24466420ha012200000","0163103100000010"); // non lin 80-100%
  } else if (trainConfig == 523) {  // JJ MC with cents
    cuts.AddCutCalo("80010123","24466420ha012200000","0163103100000010"); // non lin 0-100%
    cuts.AddCutCalo("80110123","24466420ha012200000","0163103100000010"); // non lin 0-10%
    cuts.AddCutCalo("81210123","24466420ha012200000","0163103100000010"); // non lin 10-20%
    cuts.AddCutCalo("82410123","24466420ha012200000","0163103100000010"); // non lin 20-40%
    cuts.AddCutCalo("84610123","24466420ha012200000","0163103100000010"); // non lin 40-60%
    cuts.AddCutCalo("86810123","24466420ha012200000","0163103100000010"); // non lin 60-80%
    cuts.AddCutCalo("88010123","24466420ha012200000","0163103100000010"); // non lin 80-100%

  // Variations for systematics
  } else if (trainConfig == 530) { // NL variations (standard: 42 PHOS ML)
    cuts.AddCutCalo("80010113","2446641051012200000","0163103100000010"); // PCM-PHOS NL
    cuts.AddCutCalo("80010113","2446652051012200000","0163103100000010"); // PHOS NL 2
    cuts.AddCutCalo("80010113","2446601051012200000","0163103100000010"); // PHOS group NL
  } else if (trainConfig == 531) { // distance to bad channel variations (standard: no cut)
    cuts.AddCutCalo("80010113","2446642051012200000","0163103100000010"); // dist. to bad channel = 0
    cuts.AddCutCalo("80010113","2446642151012200000","0163103100000010"); // dist. to bad channel = 1
    cuts.AddCutCalo("80010113","2446642251012200000","0163103100000010"); // dist. to bad channel = 2
    cuts.AddCutCalo("80010113","2446642351012200000","0163103100000010"); // dist. to bad channel = 3
  } else if (trainConfig == 532) { // timing variations - bunch spacing: 100ns (standard: 50ns)
    cuts.AddCutCalo("80010113","2446642071012200000","0163103100000010"); // 30ns
    cuts.AddCutCalo("80010113","2446642081012200000","0163103100000010"); // -20ns/25ns
  } else if (trainConfig == 533) { // track matching variations
    cuts.AddCutCalo("80010113","2446642050012200000","0163103100000010"); // without track matching
    cuts.AddCutCalo("80010113","2446642055012200000","0163103100000010"); // tm variation
    cuts.AddCutCalo("80010113","2446642056012200000","0163103100000010"); // tm variation
  } else if (trainConfig == 534) { // min. cluster energy variations (standard: 0.5 MeV)
    cuts.AddCutCalo("80010113","2446642051092200000","0163103100000010"); // 0.1 MeV
    cuts.AddCutCalo("80010113","2446642051072200000","0163103100000010"); // 0.2 MeV
    cuts.AddCutCalo("80010113","2446642051012200000","0163103100000010"); // 0.3 MeV
    cuts.AddCutCalo("80010113","2446642051082200000","0163103100000010"); // 0.4 MeV
    cuts.AddCutCalo("80010113","2446642051022200000","0163103100000010"); // 0.5 MeV
    cuts.AddCutCalo("80010113","2446642051032200000","0163103100000010"); // 0.6 MeV
  } else if (trainConfig == 535) { // min. number of cells per cluster variations (standard: 2)
    cuts.AddCutCalo("80010113","2446642051013200000","0163103100000010"); // min. 3 cells
    cuts.AddCutCalo("80010113","2446642051014200000","0163103100000010"); // min. 4 cells
  } else if (trainConfig == 536) { // M02/dispersion variations (standard: M02 > 0.1)
    cuts.AddCutCalo("80010113","2446642051012300000","0163103100000010"); // M02 > 0.2
    cuts.AddCutCalo("80010113","2446642051012000010","0163103100000010"); // dispersion < 2
    cuts.AddCutCalo("80010113","2446642051012000030","0163103100000010"); // dispersion < 2*2

    // Variations for systematics - JJ MC
  } else if (trainConfig == 540) { // NL variations (standard: 42 PHOS ML)
    cuts.AddCutCalo("80010123","2446641051012200000","0163103100000010"); // PCM-PHOS NL
    cuts.AddCutCalo("80010123","2446652051012200000","0163103100000010"); // PHOS NL 2
    cuts.AddCutCalo("80010123","2446601051012200000","0163103100000010"); // PHOS group NL
  } else if (trainConfig == 541) { // distance to bad channel variations (standard: no cut)
    cuts.AddCutCalo("80010123","2446642051012200000","0163103100000010"); // dist. to bad channel = 0
    cuts.AddCutCalo("80010123","2446642151012200000","0163103100000010"); // dist. to bad channel = 1
    cuts.AddCutCalo("80010123","2446642251012200000","0163103100000010"); // dist. to bad channel = 2
    cuts.AddCutCalo("80010123","2446642351012200000","0163103100000010"); // dist. to bad channel = 3
  } else if (trainConfig == 542) { // timing variations - bunch spacing: 100ns (standard: 50ns)
    cuts.AddCutCalo("80010123","2446642071012200000","0163103100000010"); // 30ns
    cuts.AddCutCalo("80010123","2446642081012200000","0163103100000010"); // -20ns/25ns
  } else if (trainConfig == 543) { // track matching variations
    cuts.AddCutCalo("80010123","2446642050012200000","0163103100000010"); // without track matching
    cuts.AddCutCalo("80010123","2446642055012200000","0163103100000010"); // tm variation
    cuts.AddCutCalo("80010123","2446642056012200000","0163103100000010"); // tm variation
  } else if (trainConfig == 544) { // min. cluster energy variations (standard: 0.5 MeV)
    cuts.AddCutCalo("80010123","2446642051092200000","0163103100000010"); // 0.1 MeV
    cuts.AddCutCalo("80010123","2446642051072200000","0163103100000010"); // 0.2 MeV
    cuts.AddCutCalo("80010123","2446642051012200000","0163103100000010"); // 0.3 MeV
    cuts.AddCutCalo("80010123","2446642051082200000","0163103100000010"); // 0.4 MeV
    cuts.AddCutCalo("80010123","2446642051022200000","0163103100000010"); // 0.5 MeV
    cuts.AddCutCalo("80010123","2446642051032200000","0163103100000010"); // 0.6 MeV
  } else if (trainConfig == 545) { // min. number of cells per cluster variations (standard: 2)
    cuts.AddCutCalo("80010123","2446642051013200000","0163103100000010"); // min. 3 cells
    cuts.AddCutCalo("80010123","2446642051014200000","0163103100000010"); // min. 4 cells
  } else if (trainConfig == 546) { // M02/dispersion variations (standard: M02 > 0.1)
    cuts.AddCutCalo("80010123","2446642051012300000","0163103100000010"); // M02 > 0.2
    cuts.AddCutCalo("80010123","2446642051012000010","0163103100000010"); // dispersion < 2
    cuts.AddCutCalo("80010123","2446642051012000030","0163103100000010"); // dispersion < 2*2

  // ===============================================================================================
  // Run 2 data PHOS clusters pPb 8TeV
  // ===============================================================================================
  // pPb 8TeV variations for QA
  } else if (trainConfig == 600){ // PHOS clusters standard cuts, triggers, no nonlin, +-50ns
    cuts.AddCutCalo("80010113","2446600051012200000","0163103100000010"); // INT7
  } else if (trainConfig == 601){ // PHOS clusters standard cuts, triggers, no nonlin, +-50ns
    cuts.AddCutCalo("80062113","2446600051012200000","0163103100000010"); // PHI7
  } else if( trainConfig == 602){ // No non-lin corr
    cuts.AddCutCalo("80010113","244660007a012200000","0163103100000010"); // no NL
    cuts.AddCutCalo("80010113","2446600070012200000","0163103100000010"); // no NL
    cuts.AddCutCalo("800ap113","2446600070012200000","0163103100000010"); // PHI7 CALOFAST
    cuts.AddCutCalo("80062113","244660007a012200000","0163103100000010"); // PHI7

  // ===============================================================================================
  // Run 2 data DMC clusters pPb 5TeV
  // ===============================================================================================
  // pPb 8TeV variations for QA
  } else if (trainConfig == 700){ // DCal clusters standard cuts, triggers, no nonlin, open timing
    cuts.AddCutCalo("80010113","3885500017032230000","01631031000000d0"); // INT7
  } else if (trainConfig == 701){ // DCal clusters standard cuts, triggers, no nonlin, open timing
    cuts.AddCutCalo("80010113","3885500057032230000","01631031000000d0"); // EMC7
  } else if (trainConfig == 702){ // DCal clusters standard cuts, triggers, no nonlin, +-50ns
    cuts.AddCutCalo("80110113","3885500057032230000","01631031000000d0"); // 0-10
    cuts.AddCutCalo("81210113","3885500057032230000","01631031000000d0"); // 10-20
    cuts.AddCutCalo("82410113","3885500057032230000","01631031000000d0"); // 20-40
    cuts.AddCutCalo("84610113","3885500057032230000","01631031000000d0"); // 40-60
    cuts.AddCutCalo("86810113","3885500057032230000","01631031000000d0"); // 60-80
    cuts.AddCutCalo("88010113","3885500057032230000","01631031000000d0"); // 80-100
  } else if (trainConfig == 703){ // DCal clusters standard cuts, triggers, no nonlin, open timing
    cuts.AddCutCalo("80110113","3885500017032230000","01631031000000d0"); // 0-10
    cuts.AddCutCalo("81210113","3885500017032230000","01631031000000d0"); // 10-20
    cuts.AddCutCalo("82410113","3885500017032230000","01631031000000d0"); // 20-40
    cuts.AddCutCalo("84610113","3885500017032230000","01631031000000d0"); // 40-60
    cuts.AddCutCalo("86810113","3885500017032230000","01631031000000d0"); // 60-80
    cuts.AddCutCalo("88010113","3885500017032230000","01631031000000d0"); // 80-100
  } else if (trainConfig == 704){ // DCal clusters standard cuts, no nonlin, +-50ns
    cuts.AddCutCalo("80210113","3885500057032230000","01631031000000d0"); // 0-20
    cuts.AddCutCalo("a0110113","3885500057032230000","01631031000000d0"); // 0-5
    cuts.AddCutCalo("a1210113","3885500057032230000","01631031000000d0"); // 5-10
    cuts.AddCutCalo("86010113","3885500057032230000","01631031000000d0"); // 60-100
  } else if (trainConfig == 705){ // DCal clusters standard cuts, no nonlin, open timing
    cuts.AddCutCalo("80210113","3885500017032230000","01631031000000d0"); // 0-20
    cuts.AddCutCalo("a0110113","3885500017032230000","01631031000000d0"); // 0-5
    cuts.AddCutCalo("a1210113","3885500017032230000","01631031000000d0"); // 5-10
    cuts.AddCutCalo("86010113","3885500017032230000","01631031000000d0"); // 60-100

  } else if (trainConfig == 780){ // DCal clusters standard cuts, triggers, no nonlin, open timing CALO+CALOFAST
    cuts.AddCutCalo("800a0113","3885500010032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("800a6113","3885500010032230000","01631031000000d0"); // DMC7
    cuts.AddCutCalo("800a7113","3885500010032230000","01631031000000d0"); // DG2
    cuts.AddCutCalo("800a8113","3885500010032230000","01631031000000d0"); // DG1
  } else if (trainConfig == 781){ // DCal clusters standard cuts, triggers, no nonlin, open timing CALO+CALOFAST
    cuts.AddCutCalo("80010113","3885500010032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("80055113","3885500010032230000","01631031000000d0"); // DMC7
    cuts.AddCutCalo("80089113","3885500010032230000","01631031000000d0"); // DG2
    cuts.AddCutCalo("8008b113","3885500010032230000","01631031000000d0"); // DG1
  // ===============================================================================================
  // Run 2 data DMC clusters pPb 8TeV
  // ===============================================================================================
  // pPb 8TeV variations for QA
  } else if (trainConfig == 800){ // DCal clusters standard cuts, triggers, no nonlin, open timing
    cuts.AddCutCalo("80010113","3885500017032230000","01631031000000d0"); // INT7
  } else if (trainConfig == 801){ // DCal clusters standard cuts, triggers, no nonlin, open timing
    cuts.AddCutCalo("00089113","3885500017032230000","01631031000000d0"); // DG2
    cuts.AddCutCalo("0008b113","3885500017032230000","01631031000000d0"); // DG1
  } else if (trainConfig == 802){ // DCal clusters standard cuts, triggers, no nonlin, open timing
    cuts.AddCutCalo("00099113","3885500017032230000","01631031000000d0"); // DJ2
    cuts.AddCutCalo("00097113","3885500017032230000","01631031000000d0"); // DJ1
  } else if (trainConfig == 804){ // DCal clusters standard cuts, triggers, no nonlin, +-50ns timing (5)
    cuts.AddCutCalo("80010113","3885500057032230000","01631031000000d0"); // INT7
  } else if (trainConfig == 805){ // DCal clusters standard cuts, triggers, no nonlin
    cuts.AddCutCalo("00089113","3885500057032230000","01631031000000d0"); // DG2
    cuts.AddCutCalo("0008b113","3885500057032230000","01631031000000d0"); // DG1
  } else if (trainConfig == 806){ // DCal clusters standard cuts, triggers, no nonlin
    cuts.AddCutCalo("00099113","3885500057032230000","01631031000000d0"); // DJ2
    cuts.AddCutCalo("00097113","3885500057032230000","01631031000000d0"); // DJ1
  } else if (trainConfig == 807){ // DCal clusters standard cuts, triggers, Nico TB NL
    cuts.AddCutCalo("80010113","3885565057032230000","01631031000000d0"); // INT7
  } else if (trainConfig == 808){ // DCal clusters standard cuts, triggers, Nico TB NL
    cuts.AddCutCalo("00089113","3885565057032230000","01631031000000d0"); // DG2
    cuts.AddCutCalo("0008b113","3885565057032230000","01631031000000d0"); // DG1
  } else if (trainConfig == 809){ // DCal clusters standard cuts, triggers, Nico TB NL
    cuts.AddCutCalo("00099113","3885565057032230000","01631031000000d0"); // DJ2
    cuts.AddCutCalo("00097113","3885565057032230000","01631031000000d0"); // DJ1

  } else if (trainConfig == 810){ // DCal clusters standard cuts, triggers, NL variations
    cuts.AddCutCalo("80010113","3885547057032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("80010113","3885548057032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("80010113","3885557057032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("80010113","3885558057032230000","01631031000000d0"); // INT7
  } else if (trainConfig == 811){ // DCal clusters standard cuts, triggers, kSDM NL DMC
    cuts.AddCutCalo("00089113","3885547057032230000","01631031000000d0"); // DG2
    cuts.AddCutCalo("00089113","3885548057032230000","01631031000000d0"); // DG2
    cuts.AddCutCalo("00089113","3885557057032230000","01631031000000d0"); // DG2
    cuts.AddCutCalo("00089113","3885558057032230000","01631031000000d0"); // DG2
  } else if (trainConfig == 812){ // DCal clusters standard cuts, triggers, kSDM NL DMC
    cuts.AddCutCalo("0008b113","3885547057032230000","01631031000000d0"); // DG1
    cuts.AddCutCalo("0008b113","3885548057032230000","01631031000000d0"); // DG1
    cuts.AddCutCalo("0008b113","3885557057032230000","01631031000000d0"); // DG1
    cuts.AddCutCalo("0008b113","3885558057032230000","01631031000000d0"); // DG1

  } else if (trainConfig == 813){ // DCal clusters standard cuts, triggers, NL variations, TM E/p 1.75 (f)
    cuts.AddCutCalo("80010113","388554705f032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("80010113","388554805f032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("80010113","388555705f032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("80010113","388555805f032230000","01631031000000d0"); // INT7
  } else if (trainConfig == 814){ // DCal clusters standard cuts, triggers, kSDM NL DMC, TM E/p 1.75 (f)
    cuts.AddCutCalo("00089113","388554705f032230000","01631031000000d0"); // DG2
    cuts.AddCutCalo("00089113","388554805f032230000","01631031000000d0"); // DG2
    cuts.AddCutCalo("00089113","388555705f032230000","01631031000000d0"); // DG2
    cuts.AddCutCalo("00089113","388555805f032230000","01631031000000d0"); // DG2
  } else if (trainConfig == 815){ // DCal clusters standard cuts, triggers, kSDM NL DMC, TM E/p 1.75 (f)
    cuts.AddCutCalo("0008b113","388554705f032230000","01631031000000d0"); // DG1
    cuts.AddCutCalo("0008b113","388554805f032230000","01631031000000d0"); // DG1
    cuts.AddCutCalo("0008b113","388555705f032230000","01631031000000d0"); // DG1
    cuts.AddCutCalo("0008b113","388555805f032230000","01631031000000d0"); // DG1

  } else if (trainConfig == 820){ // nonlinearity variations
    cuts.AddCutCalo("80010113","3885541057032230000","01631031000000d0"); // CRF
    cuts.AddCutCalo("80010113","3885552057032230000","01631031000000d0"); // CRF
    cuts.AddCutCalo("80010113","3885551057032230000","01631031000000d0"); // CCMF
    cuts.AddCutCalo("80010113","3885552057032230000","01631031000000d0"); // CMF
  } else if (trainConfig == 821){ // second set of variations CLUSTER
    cuts.AddCutCalo("80010113","3885552057022230000","01631031000000d0"); // min energy cluster variation 1  600 MeV
    cuts.AddCutCalo("80010113","3885552057042230000","01631031000000d0"); // min energy cluster variation 2  800 MeV
    cuts.AddCutCalo("80010113","3885552057052230000","01631031000000d0"); // min energy cluster variation 3  900 MeV
    cuts.AddCutCalo("80010113","3885552057032230000","01631031000000d0"); // min/max M02  0.1<M<0.5
  } else if (trainConfig == 822){ // third set of variations CLUSTER
    cuts.AddCutCalo("80010113","3885552057032250000","01631031000000d0"); // min/max M02  0.1<M<0.3
    cuts.AddCutCalo("80010113","3885552057032260000","01631031000000d0"); // min/max M02  0.1<M<0.27
    cuts.AddCutCalo("80010113","3885552057031230000","01631031000000d0"); // min number of cells variation 1  1 cell
  } else if (trainConfig == 823){ // third set of variations MESON
    cuts.AddCutCalo("80010113","3885552057032230000","01633031000000d0"); // rapidity variation  y<0.6
    cuts.AddCutCalo("80010113","3885552057032230000","01634031000000d0"); // rapidity variation  y<0.5
    cuts.AddCutCalo("80010113","3885552057032230000","01631061000000d0"); // alpha meson variation 1   0<alpha<0.8
    cuts.AddCutCalo("80010113","3885552057032230000","01631051000000d0"); // alpha meson variation 2  0<alpha<0.75
  } else if (trainConfig == 824){ // opening angle variations
    cuts.AddCutCalo("80010113","3885552057032230000","0163103100000040"); // min opening angle 0.0152
    cuts.AddCutCalo("80010113","3885552057032230000","0163103100000060"); // min opening angle 0.017
    cuts.AddCutCalo("80010113","3885552057032230000","0163103100000070"); // min opening angle 0.016
    cuts.AddCutCalo("80010113","3885552057032230000","0163103100000080"); // min opening angle 0.018
    cuts.AddCutCalo("80010113","3885552057032230000","0163103100000090"); // min opening angle 0.018
  } else if (trainConfig == 825){ // TM variations
    cuts.AddCutCalo("80010113","3885552053032230000","01631031000000d0"); // fixed window
    cuts.AddCutCalo("80010113","3885552056032230000","01631031000000d0"); // tm pt dependent var 1
    cuts.AddCutCalo("80010113","3885552058032230000","01631031000000d0"); // tm pt dependent var 2
    cuts.AddCutCalo("80010113","3885552059032230000","01631031000000d0"); // tm pt dependent var 3
  } else if (trainConfig == 826){
    cuts.AddCutCalo("80010113","38855520570322d0000","01631031000000d0"); // M02, pt dep with  0.27, 0.0072, 0.4
    cuts.AddCutCalo("80010113","38855520570322e0000","01631031000000d0"); // M02, pt dep with  0.31, 0.0072, 0.5
    cuts.AddCutCalo("80010113","38855520570322f0000","01631031000000d0"); // M02, pt dep with  0.36, 0.0072, 0.7
  } else if (trainConfig == 827){
    cuts.AddCutCalo("80010113","38855520570322m0000","01631031000000d0"); // M02, pt dep with  0.32, 0.0152, 0.5
    cuts.AddCutCalo("80010113","38855520570322g0000","01631031000000d0"); // M02, pt dep with  0.37, 0.0072, 0.7
    cuts.AddCutCalo("80010113","38855520570322h0000","01631031000000d0"); // M02, pT-dep with  0.30, 0.0072, 0.5
    cuts.AddCutCalo("80010113","38855520570322i0000","01631031000000d0"); // M02, pT-dep with  0.35, 0.0072, 0.7

  // ===============================================================================================
  // Run 1 data EMC triggers only
  // ===============================================================================================
  // EMC7 triggers
  } else if (trainConfig == 1000){ // QA setups
    cuts.AddCutCalo("80052113","1111100007032230000","01631031000000d0");
    cuts.AddCutCalo("80052113","1111100057032230000","01631031000000d0");
  } else if (trainConfig == 1001){ // new default cuts
    cuts.AddCutCalo("80052113","1111141057032230000","01631031000000d0"); // meson
    cuts.AddCutCalo("80052113","11111410570322l0000","01631031000000d0"); // dir gamma
  } else if (trainConfig == 1002){
    cuts.AddCutCalo("80052113","1111142057032230000","01631031000000d0"); // CRF
    cuts.AddCutCalo("80052113","1111151057032230000","01631031000000d0"); // CCMF
    cuts.AddCutCalo("80052113","1111152057032230000","01631031000000d0"); // CMF
  } else if (trainConfig == 1003){ // second set of variations CLUSTER
    cuts.AddCutCalo("80052113","1111141057022230000","01631031000000d0"); // min energy cluster variation 1  600 MeV
    cuts.AddCutCalo("80052113","1111141057042230000","01631031000000d0"); // min energy cluster variation 2  800 MeV
    cuts.AddCutCalo("80052113","1111141057052230000","01631031000000d0"); // min energy cluster variation 3  900 MeV
    cuts.AddCutCalo("80052113","1111141057032230000","01631031000000d0"); // min/max M02  0.1<M<0.5
    cuts.AddCutCalo("80052113","1111141057032200000","01631031000000d0"); // min/max M02  0.1<M<100
  } else if (trainConfig == 1004){ // third set of variations CLUSTER
    cuts.AddCutCalo("80052113","1111141057032250000","01631031000000d0"); // min/max M02  0.1<M<0.3
    cuts.AddCutCalo("80052113","1111141057032260000","01631031000000d0"); // min/max M02  0.1<M<0.27
    cuts.AddCutCalo("80052113","1111141057031230000","01631031000000d0"); // min number of cells variation 1  1 cell
    cuts.AddCutCalo("80052113","1112141057032230000","01631031000000d0"); // only modules with TRD infront
    cuts.AddCutCalo("80052113","1111341057032230000","01631031000000d0"); // no modules with TRD infront
  } else if (trainConfig == 1005){ // third set of variations MESON
    cuts.AddCutCalo("80052113","1111141057032230000","01633031000000d0"); // rapidity variation  y<0.6
    cuts.AddCutCalo("80052113","1111141057032230000","01634031000000d0"); // rapidity variation  y<0.5
    cuts.AddCutCalo("80052113","1111141057032230000","01631061000000d0"); // alpha meson variation 1   0<alpha<0.8
    cuts.AddCutCalo("80052113","1111141057032230000","01631051000000d0"); // alpha meson variation 2  0<alpha<0.75
  } else if (trainConfig == 1006){ // opening angle variations
    cuts.AddCutCalo("80052113","1111141057032230000","0163103100000040"); // min opening angle 0.0152
    cuts.AddCutCalo("80052113","1111141057032230000","0163103100000060"); // min opening angle 0.017
    cuts.AddCutCalo("80052113","1111141057032230000","0163103100000070"); // min opening angle 0.016
    cuts.AddCutCalo("80052113","1111141057032230000","0163103100000080"); // min opening angle 0.018
    cuts.AddCutCalo("80052113","1111141057032230000","0163103100000090"); // min opening angle 0.018
  } else if (trainConfig == 1007){ // TM variations
    cuts.AddCutCalo("80052113","1111141053032230000","01631031000000d0"); // fixed window
    cuts.AddCutCalo("80052113","1111141056032230000","01631031000000d0"); // tm pt dependent var 1
    cuts.AddCutCalo("80052113","1111141058032230000","01631031000000d0"); // tm pt dependent var 2
    cuts.AddCutCalo("80052113","1111141059032230000","01631031000000d0"); // tm pt dependent var 3
  } else if (trainConfig == 1008){ // TM variations
    cuts.AddCutCalo("80052113","1111153057032230000","01631031000000d0"); // CCMF nonlin+testbeam
    cuts.AddCutCalo("80052113","1111154057032230000","01631031000000d0"); // CMF nonlin+testbeam
    cuts.AddCutCalo("80052113","1111143057032230000","01631031000000d0"); // CCRF testbeam nonlin
    cuts.AddCutCalo("80052113","1111144057032230000","01631031000000d0"); // CRF testbeam nonlin
    cuts.AddCutCalo("80052113","1111102057032230000","01631031000000d0"); // testbeam nonlin
  } else if (trainConfig == 1009){
    cuts.AddCutCalo("80052113","11111410570322l0000","01631031000000d0"); // M02 pt dep with new std: 0.32, 0.0072, 0.5
    cuts.AddCutCalo("80052113","11111410570322d0000","01631031000000d0"); // M02, pt dep with  0.27, 0.0072, 0.4
    cuts.AddCutCalo("80052113","11111410570322e0000","01631031000000d0"); // M02, pt dep with  0.31, 0.0072, 0.5
    cuts.AddCutCalo("80052113","11111410570322f0000","01631031000000d0"); // M02, pt dep with  0.36, 0.0072, 0.7
  } else if (trainConfig == 1010){
    cuts.AddCutCalo("80052113","11111410570322m0000","01631031000000d0"); // M02, pt dep with  0.32, 0.0152, 0.5
    cuts.AddCutCalo("80052113","11111410570322g0000","01631031000000d0"); // M02, pt dep with  0.37, 0.0072, 0.7
    cuts.AddCutCalo("80052113","11111410570322h0000","01631031000000d0"); // M02, pT-dep with  0.30, 0.0072, 0.5
    cuts.AddCutCalo("80052113","11111410570322i0000","01631031000000d0"); // M02, pT-dep with  0.35, 0.0072, 0.7
  } else if (trainConfig == 1011){
    cuts.AddCutCalo("80052113","11111410570322j0000","01631031000000d0"); // M02, pT-dep with  0.25, 0.0072, 0.39
    cuts.AddCutCalo("80052113","11111410570322r0000","01631031000000d0"); // M02, pT-dep with  0.25, 0.0072, 0.5
    cuts.AddCutCalo("80052113","11111410570322s0000","01631031000000d0"); // M02, pT-dep with  0.32, 0.0238, 0.7
    cuts.AddCutCalo("80052113","11111410570322n0000","01631031000000d0"); // M02, pT-dep with  0.32, 0.0238, 0.7

  } else if (trainConfig == 1012){
    cuts.AddCutCalo("80252113","1111141057032230000","01631031000000d0"); // 0-20 EMC7
    cuts.AddCutCalo("82452113","1111141057032230000","01631031000000d0"); // 20-40 EMC7
    cuts.AddCutCalo("84652113","1111141057032230000","01631031000000d0"); // 40-60 EMC7
    cuts.AddCutCalo("86052113","1111141057032230000","01631031000000d0"); // 60-80 EMC7
  } else if (trainConfig == 1013){
    cuts.AddCutCalo("80252113","11111410570322l0000","01631031000000d0"); // 0-20 EMC7
    cuts.AddCutCalo("82452113","11111410570322l0000","01631031000000d0"); // 20-40 EMC7
    cuts.AddCutCalo("84652113","11111410570322l0000","01631031000000d0"); // 40-60 EMC7
    cuts.AddCutCalo("86052113","11111410570322l0000","01631031000000d0"); // 60-80 EMC7
  } else if (trainConfig == 1014){ // new default cuts with sector mixing
    cuts.AddCutCalo("80052113","1111141057032230000","0i631031000000d0"); // meson
    cuts.AddCutCalo("80052113","11111410570322l0000","0i631031000000d0"); // dir gamma
  } else if (trainConfig == 1015){
    cuts.AddCutCalo("80252113","1111141057032230000","0i631031000000d0"); // 0-20 EMC7
    cuts.AddCutCalo("82452113","1111141057032230000","0i631031000000d0"); // 20-40 EMC7
    cuts.AddCutCalo("84652113","1111141057032230000","0i631031000000d0"); // 40-60 EMC7
    cuts.AddCutCalo("86052113","1111141057032230000","0i631031000000d0"); // 60-80 EMC7
  } else if (trainConfig == 1016){
    cuts.AddCutCalo("80252113","11111410570322l0000","0i631031000000d0"); // 0-20 EMC7
    cuts.AddCutCalo("82452113","11111410570322l0000","0i631031000000d0"); // 20-40 EMC7
    cuts.AddCutCalo("84652113","11111410570322l0000","0i631031000000d0"); // 40-60 EMC7
    cuts.AddCutCalo("86052113","11111410570322l0000","0i631031000000d0"); // 60-80 EMC7


  // EG2 triggers
  } else if (trainConfig == 1020){  // QA setups
    cuts.AddCutCalo("80085113","1111100007032230000","01631031000000d0");
    cuts.AddCutCalo("80085113","1111100057032230000","01631031000000d0");
  } else if (trainConfig == 1021){ // new default cuts
    cuts.AddCutCalo("80085113","1111141057032230000","01631031000000d0"); // meson
    cuts.AddCutCalo("80085113","11111410570322l0000","01631031000000d0"); // dir gamma
  } else if (trainConfig == 1022){
    cuts.AddCutCalo("80085113","1111142057032230000","01631031000000d0"); // CRF
    cuts.AddCutCalo("80085113","1111151057032230000","01631031000000d0"); // CCMF
    cuts.AddCutCalo("80085113","1111152057032230000","01631031000000d0"); // CMF
  } else if (trainConfig == 1023){ // second set of variations CLUSTER
    cuts.AddCutCalo("80085113","1111141057022230000","01631031000000d0"); // min energy cluster variation 1  600 MeV
    cuts.AddCutCalo("80085113","1111141057042230000","01631031000000d0"); // min energy cluster variation 2  800 MeV
    cuts.AddCutCalo("80085113","1111141057052230000","01631031000000d0"); // min energy cluster variation 3  900 MeV
    cuts.AddCutCalo("80085113","1111141057032230000","01631031000000d0"); // min/max M02  0.1<M<0.5
    cuts.AddCutCalo("80085113","1111141057032200000","01631031000000d0"); // min/max M02  0.1<M<100
  } else if (trainConfig == 1024){ // third set of variations CLUSTER
    cuts.AddCutCalo("80085113","1111141057032250000","01631031000000d0"); // min/max M02  0.1<M<0.3
    cuts.AddCutCalo("80085113","1111141057032260000","01631031000000d0"); // min/max M02  0.1<M<0.27
    cuts.AddCutCalo("80085113","1111141057031230000","01631031000000d0"); // min number of cells variation 1  1 cell
    cuts.AddCutCalo("80085113","1112141057032230000","01631031000000d0"); // only modules with TRD infront
    cuts.AddCutCalo("80085113","1111341057032230000","01631031000000d0"); // no modules with TRD infront
  } else if (trainConfig == 1025){ // third set of variations MESON
    cuts.AddCutCalo("80085113","1111141057032230000","01633031000000d0"); // rapidity variation  y<0.6
    cuts.AddCutCalo("80085113","1111141057032230000","01634031000000d0"); // rapidity variation  y<0.5
    cuts.AddCutCalo("80085113","1111141057032230000","01631061000000d0"); // alpha meson variation 1   0<alpha<0.8
    cuts.AddCutCalo("80085113","1111141057032230000","01631051000000d0"); // alpha meson variation 2  0<alpha<0.75
  } else if (trainConfig == 1026){ // opening angle variations
    cuts.AddCutCalo("80085113","1111141057032230000","0163103100000040"); // min opening angle 0.0152
    cuts.AddCutCalo("80085113","1111141057032230000","0163103100000060"); // min opening angle 0.017
    cuts.AddCutCalo("80085113","1111141057032230000","0163103100000070"); // min opening angle 0.016
    cuts.AddCutCalo("80085113","1111141057032230000","0163103100000080"); // min opening angle 0.018
    cuts.AddCutCalo("80085113","1111141057032230000","0163103100000090"); // min opening angle 0.018
  } else if (trainConfig == 1027){ // TM variations
    cuts.AddCutCalo("80085113","1111141053032230000","01631031000000d0"); // fixed window
    cuts.AddCutCalo("80085113","1111141056032230000","01631031000000d0"); // tm pt dependent var 1
    cuts.AddCutCalo("80085113","1111141058032230000","01631031000000d0"); // tm pt dependent var 2
    cuts.AddCutCalo("80085113","1111141059032230000","01631031000000d0"); // tm pt dependent var 3
  } else if (trainConfig == 1028){ // TM variations
    cuts.AddCutCalo("80085113","1111153057032230000","01631031000000d0"); // CCMF nonlin+testbeam
    cuts.AddCutCalo("80085113","1111154057032230000","01631031000000d0"); // CMF nonlin+testbeam
    cuts.AddCutCalo("80085113","1111143057032230000","01631031000000d0"); // CCRF testbeam nonlin
    cuts.AddCutCalo("80085113","1111144057032230000","01631031000000d0"); // CRF testbeam nonlin
    cuts.AddCutCalo("80085113","1111102057032230000","01631031000000d0"); // testbeam nonlin
  } else if (trainConfig == 1029){
    cuts.AddCutCalo("80085113","11111410570322l0000","01631031000000d0"); // M02 pt dep with new std: 0.32, 0.0072, 0.5
    cuts.AddCutCalo("80085113","11111410570322d0000","01631031000000d0"); // M02, pt dep with  0.27, 0.0072, 0.4
    cuts.AddCutCalo("80085113","11111410570322e0000","01631031000000d0"); // M02, pt dep with  0.31, 0.0072, 0.5
    cuts.AddCutCalo("80085113","11111410570322f0000","01631031000000d0"); // M02, pt dep with  0.36, 0.0072, 0.7
  } else if (trainConfig == 1030){
    cuts.AddCutCalo("80085113","11111410570322m0000","01631031000000d0"); // M02, pt dep with  0.32, 0.0152, 0.5
    cuts.AddCutCalo("80085113","11111410570322g0000","01631031000000d0"); // M02, pt dep with  0.37, 0.0072, 0.7
    cuts.AddCutCalo("80085113","11111410570322h0000","01631031000000d0"); // M02, pT-dep with  0.30, 0.0072, 0.5
    cuts.AddCutCalo("80085113","11111410570322i0000","01631031000000d0"); // M02, pT-dep with  0.35, 0.0072, 0.7
  } else if (trainConfig == 1031){
    cuts.AddCutCalo("80085113","11111410570322j0000","01631031000000d0"); // M02, pT-dep with  0.25, 0.0072, 0.39
    cuts.AddCutCalo("80085113","11111410570322r0000","01631031000000d0"); // M02, pT-dep with  0.25, 0.0072, 0.5
    cuts.AddCutCalo("80085113","11111410570322s0000","01631031000000d0"); // M02, pT-dep with  0.32, 0.0238, 0.7
    cuts.AddCutCalo("80085113","11111410570322n0000","01631031000000d0"); // M02, pT-dep with  0.32, 0.0238, 0.7

  } else if (trainConfig == 1032){
    cuts.AddCutCalo("80285113","1111141057032230000","01631031000000d0"); // 0-20 EG2
    cuts.AddCutCalo("82485113","1111141057032230000","01631031000000d0"); // 20-40 EG2
    cuts.AddCutCalo("84685113","1111141057032230000","01631031000000d0"); // 40-60 EG2
    cuts.AddCutCalo("86085113","1111141057032230000","01631031000000d0"); // 60-80 EG2
  } else if (trainConfig == 1033){
    cuts.AddCutCalo("80285113","11111410570322l0000","01631031000000d0"); // 0-20 EG2
    cuts.AddCutCalo("82485113","11111410570322l0000","01631031000000d0"); // 20-40 EG2
    cuts.AddCutCalo("84685113","11111410570322l0000","01631031000000d0"); // 40-60 EG2
    cuts.AddCutCalo("86085113","11111410570322l0000","01631031000000d0"); // 60-80 EG2
  } else if (trainConfig == 1034){ // new default cuts with sector mixing
    cuts.AddCutCalo("80085113","1111141057032230000","0i631031000000d0"); // meson
    cuts.AddCutCalo("80085113","11111410570322l0000","0i631031000000d0"); // dir gamma
  } else if (trainConfig == 1035){
    cuts.AddCutCalo("80285113","1111141057032230000","0i631031000000d0"); // 0-20 EG2
    cuts.AddCutCalo("82485113","1111141057032230000","0i631031000000d0"); // 20-40 EG2
    cuts.AddCutCalo("84685113","1111141057032230000","0i631031000000d0"); // 40-60 EG2
    cuts.AddCutCalo("86085113","1111141057032230000","0i631031000000d0"); // 60-80 EG2
  } else if (trainConfig == 1036){
    cuts.AddCutCalo("80285113","11111410570322l0000","0i631031000000d0"); // 0-20 EG2
    cuts.AddCutCalo("82485113","11111410570322l0000","0i631031000000d0"); // 20-40 EG2
    cuts.AddCutCalo("84685113","11111410570322l0000","0i631031000000d0"); // 40-60 EG2
    cuts.AddCutCalo("86085113","11111410570322l0000","0i631031000000d0"); // 60-80 EG2


  // EG1 triggers
  } else if (trainConfig == 1040){  // QA setups
    cuts.AddCutCalo("80083113","1111100007032230000","01631031000000d0");
    cuts.AddCutCalo("80083113","1111100057032230000","01631031000000d0");
  } else if (trainConfig == 1041){ // new default cuts
    cuts.AddCutCalo("80083113","1111141057032230000","01631031000000d0"); // meson
    cuts.AddCutCalo("80083113","11111410570322l0000","01631031000000d0"); // dir gamma
  } else if (trainConfig == 1042){
    cuts.AddCutCalo("80083113","1111142057032230000","01631031000000d0"); // CRF
    cuts.AddCutCalo("80083113","1111151057032230000","01631031000000d0"); // CCMF
    cuts.AddCutCalo("80083113","1111152057032230000","01631031000000d0"); // CMF
  } else if (trainConfig == 1043){ // second set of variations CLUSTER
    cuts.AddCutCalo("80083113","1111141057022230000","01631031000000d0"); // min energy cluster variation 1  600 MeV
    cuts.AddCutCalo("80083113","1111141057042230000","01631031000000d0"); // min energy cluster variation 2  800 MeV
    cuts.AddCutCalo("80083113","1111141057052230000","01631031000000d0"); // min energy cluster variation 3  900 MeV
    cuts.AddCutCalo("80083113","1111141057032230000","01631031000000d0"); // min/max M02  0.1<M<0.5
    cuts.AddCutCalo("80083113","1111141057032200000","01631031000000d0"); // min/max M02  0.1<M<100
  } else if (trainConfig == 1044){ // third set of variations CLUSTER
    cuts.AddCutCalo("80083113","1111141057032250000","01631031000000d0"); // min/max M02  0.1<M<0.3
    cuts.AddCutCalo("80083113","1111141057032260000","01631031000000d0"); // min/max M02  0.1<M<0.27
    cuts.AddCutCalo("80083113","1111141057031230000","01631031000000d0"); // min number of cells variation 1  1 cell
    cuts.AddCutCalo("80083113","1112141057032230000","01631031000000d0"); // only modules with TRD infront
    cuts.AddCutCalo("80083113","1111341057032230000","01631031000000d0"); // no modules with TRD infront
  } else if (trainConfig == 1045){ // third set of variations MESON
    cuts.AddCutCalo("80083113","1111141057032230000","01633031000000d0"); // rapidity variation  y<0.6
    cuts.AddCutCalo("80083113","1111141057032230000","01634031000000d0"); // rapidity variation  y<0.5
    cuts.AddCutCalo("80083113","1111141057032230000","01631061000000d0"); // alpha meson variation 1   0<alpha<0.8
    cuts.AddCutCalo("80083113","1111141057032230000","01631051000000d0"); // alpha meson variation 2  0<alpha<0.75
  } else if (trainConfig == 1046){ // opening angle variations
    cuts.AddCutCalo("80083113","1111141057032230000","0163103100000040"); // min opening angle 0.0152
    cuts.AddCutCalo("80083113","1111141057032230000","0163103100000060"); // min opening angle 0.017
    cuts.AddCutCalo("80083113","1111141057032230000","0163103100000070"); // min opening angle 0.016
    cuts.AddCutCalo("80083113","1111141057032230000","0163103100000080"); // min opening angle 0.018
    cuts.AddCutCalo("80083113","1111141057032230000","0163103100000090"); // min opening angle 0.018
  } else if (trainConfig == 1047){ // TM variations
    cuts.AddCutCalo("80083113","1111141053032230000","01631031000000d0"); // fixed window
    cuts.AddCutCalo("80083113","1111141056032230000","01631031000000d0"); // tm pt dependent var 1
    cuts.AddCutCalo("80083113","1111141058032230000","01631031000000d0"); // tm pt dependent var 2
    cuts.AddCutCalo("80083113","1111141059032230000","01631031000000d0"); // tm pt dependent var 3
  } else if (trainConfig == 1048){ // TM variations
    cuts.AddCutCalo("80083113","1111153057032230000","01631031000000d0"); // CCMF nonlin+testbeam
    cuts.AddCutCalo("80083113","1111154057032230000","01631031000000d0"); // CMF nonlin+testbeam
    cuts.AddCutCalo("80083113","1111143057032230000","01631031000000d0"); // CCRF testbeam nonlin
    cuts.AddCutCalo("80083113","1111144057032230000","01631031000000d0"); // CRF testbeam nonlin
    cuts.AddCutCalo("80083113","1111102057032230000","01631031000000d0"); // testbeam nonlin
  } else if (trainConfig == 1049){
    cuts.AddCutCalo("80083113","11111410570322l0000","01631031000000d0"); // M02 pt dep with new std: 0.32, 0.0072, 0.5
    cuts.AddCutCalo("80083113","11111410570322d0000","01631031000000d0"); // M02, pt dep with  0.27, 0.0072, 0.4
    cuts.AddCutCalo("80083113","11111410570322e0000","01631031000000d0"); // M02, pt dep with  0.31, 0.0072, 0.5
    cuts.AddCutCalo("80083113","11111410570322f0000","01631031000000d0"); // M02, pt dep with  0.36, 0.0072, 0.7
  } else if (trainConfig == 1050){
    cuts.AddCutCalo("80083113","11111410570322m0000","01631031000000d0"); // M02, pt dep with  0.32, 0.0152, 0.5
    cuts.AddCutCalo("80083113","11111410570322g0000","01631031000000d0"); // M02, pt dep with  0.37, 0.0072, 0.7
    cuts.AddCutCalo("80083113","11111410570322h0000","01631031000000d0"); // M02, pT-dep with  0.30, 0.0072, 0.5
    cuts.AddCutCalo("80083113","11111410570322i0000","01631031000000d0"); // M02, pT-dep with  0.35, 0.0072, 0.7
  } else if (trainConfig == 1051){
    cuts.AddCutCalo("80083113","11111410570322j0000","01631031000000d0"); // M02, pT-dep with  0.25, 0.0072, 0.39
    cuts.AddCutCalo("80083113","11111410570322r0000","01631031000000d0"); // M02, pT-dep with  0.25, 0.0072, 0.5
    cuts.AddCutCalo("80083113","11111410570322s0000","01631031000000d0"); // M02, pT-dep with  0.32, 0.0238, 0.7
    cuts.AddCutCalo("80083113","11111410570322n0000","01631031000000d0"); // M02, pT-dep with  0.32, 0.0238, 0.7

  } else if (trainConfig == 1052){
    cuts.AddCutCalo("80283113","1111141057032230000","01631031000000d0"); // 0-20 EG1
    cuts.AddCutCalo("82483113","1111141057032230000","01631031000000d0"); // 20-40 EG1
    cuts.AddCutCalo("84683113","1111141057032230000","01631031000000d0"); // 40-60 EG1
    cuts.AddCutCalo("86083113","1111141057032230000","01631031000000d0"); // 60-80 EG1
  } else if (trainConfig == 1053){
    cuts.AddCutCalo("80283113","11111410570322l0000","01631031000000d0"); // 0-20 EG1
    cuts.AddCutCalo("82483113","11111410570322l0000","01631031000000d0"); // 20-40 EG1
    cuts.AddCutCalo("84683113","11111410570322l0000","01631031000000d0"); // 40-60 EG1
    cuts.AddCutCalo("86083113","11111410570322l0000","01631031000000d0"); // 60-80 EG1
  } else if (trainConfig == 1054){ // new default cuts with sector mixing
    cuts.AddCutCalo("80083113","1111141057032230000","0i631031000000d0"); // meson
    cuts.AddCutCalo("80083113","11111410570322l0000","0i631031000000d0"); // dir gamma
  } else if (trainConfig == 1055){
    cuts.AddCutCalo("80283113","1111141057032230000","0i631031000000d0"); // 0-20 EG1
    cuts.AddCutCalo("82483113","1111141057032230000","0i631031000000d0"); // 20-40 EG1
    cuts.AddCutCalo("84683113","1111141057032230000","0i631031000000d0"); // 40-60 EG1
    cuts.AddCutCalo("86083113","1111141057032230000","0i631031000000d0"); // 60-80 EG1
  } else if (trainConfig == 1056){
    cuts.AddCutCalo("80283113","11111410570322l0000","0i631031000000d0"); // 0-20 EG1
    cuts.AddCutCalo("82483113","11111410570322l0000","0i631031000000d0"); // 20-40 EG1
    cuts.AddCutCalo("84683113","11111410570322l0000","0i631031000000d0"); // 40-60 EG1
    cuts.AddCutCalo("86083113","11111410570322l0000","0i631031000000d0"); // 60-80 EG1

  // configurations for pPb 5 and 8 TeV Run2 with EMCAL + DCAL
  } else if (trainConfig == 2000){ // EMCAL+DCAL clusters standard cuts, triggers, NL vars
    cuts.AddCutCalo("80010123","4117947057032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("80010123","4117948057032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("80010123","4117957057032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("80010123","4117958057032230000","01631031000000d0"); // INT7
  } else if (trainConfig == 2001){ // EMCAL+DCAL clusters standard cuts, triggers, NL vars
    cuts.AddCutCalo("8008e123","4117947057032230000","01631031000000d0"); // EG2
    cuts.AddCutCalo("8008e123","4117948057032230000","01631031000000d0"); // EG2
    cuts.AddCutCalo("8008e123","4117957057032230000","01631031000000d0"); // EG2
    cuts.AddCutCalo("8008e123","4117958057032230000","01631031000000d0"); // EG2
  } else if (trainConfig == 2002){ // EMCAL+DCAL clusters standard cuts, triggers, NL vars
    cuts.AddCutCalo("8008d123","4117947057032230000","01631031000000d0"); // EG1
    cuts.AddCutCalo("8008d123","4117948057032230000","01631031000000d0"); // EG1
    cuts.AddCutCalo("8008d123","4117957057032230000","01631031000000d0"); // EG1
    cuts.AddCutCalo("8008d123","4117958057032230000","01631031000000d0"); // EG1
  } else if (trainConfig == 2003){ // EMCAL+DCAL clusters standard cuts, triggers, NL vars
    cuts.AddCutCalo("80010123","4117958057032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("8008e123","4117958057032230000","01631031000000d0"); // EG2+DG2
    cuts.AddCutCalo("8008d123","4117958057032230000","01631031000000d0"); // EG1+DG1
  } else if (trainConfig == 2004){ // EMCAL clusters standard cuts, triggers, NL vars
    cuts.AddCutCalo("80010123","1111158057032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("80085123","1111158057032230000","01631031000000d0"); // EG2
    cuts.AddCutCalo("80083123","1111158057032230000","01631031000000d0"); // EG1
  } else if (trainConfig == 2005){ // DCAL clusters standard cuts, triggers, NL vars
    cuts.AddCutCalo("80010123","3885558057032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("8008b123","3885558057032230000","01631031000000d0"); // DG2
    cuts.AddCutCalo("80089123","3885558057032230000","01631031000000d0"); // DG1
  } else if (trainConfig == 2006){ // EMCAL+DCAL clusters standard cuts, triggers, NL vars
    cuts.AddCutCalo("80010123","4117900057032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("80010123","4117906057032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("80010123","4117965057032230000","01631031000000d0"); // INT7
  } else if (trainConfig == 2007){ // EMCAL+DCAL clusters standard cuts, triggers, NL vars
    cuts.AddCutCalo("8008e123","4117900057032230000","01631031000000d0"); // EG2
    cuts.AddCutCalo("8008e123","4117906057032230000","01631031000000d0"); // EG2
    cuts.AddCutCalo("8008e123","4117965057032230000","01631031000000d0"); // EG2
  } else if (trainConfig == 2008){ // EMCAL+DCAL clusters standard cuts, triggers, NL vars
    cuts.AddCutCalo("8008d123","4117900057032230000","01631031000000d0"); // EG1
    cuts.AddCutCalo("8008d123","4117906057032230000","01631031000000d0"); // EG1
    cuts.AddCutCalo("8008d123","4117965057032230000","01631031000000d0"); // EG1
  // includes both stripes EMCal and DCal
  } else if (trainConfig == 2010){ // EMCAL+DCAL clusters standard cuts, triggers, NL vars
    cuts.AddCutCalo("80010123","4117900057032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("80010123","4117906057032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("80010123","4117965057032230000","01631031000000d0"); // INT7
  } else if (trainConfig == 2011){ // EMCAL+DCAL clusters standard cuts, triggers, NL vars
    cuts.AddCutCalo("8008e123","4117900057032230000","01631031000000d0"); // EG2
    cuts.AddCutCalo("8008e123","4117906057032230000","01631031000000d0"); // EG2
    cuts.AddCutCalo("8008e123","4117965057032230000","01631031000000d0"); // EG2
  } else if (trainConfig == 2012){ // EMCAL+DCAL clusters standard cuts, triggers, NL vars
    cuts.AddCutCalo("8008d123","4117900057032230000","01631031000000d0"); // EG1
    cuts.AddCutCalo("8008d123","4117906057032230000","01631031000000d0"); // EG1
    cuts.AddCutCalo("8008d123","4117965057032230000","01631031000000d0"); // EG1

  } else if (trainConfig == 2013){ // EMCAL+DCAL clusters standard cuts, triggers, NL vars
    cuts.AddCutCalo("80010123","4117900057032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("8008e123","4117900057032230000","01631031000000d0"); // EG2
    cuts.AddCutCalo("8008d123","4117900057032230000","01631031000000d0"); // EG1
  } else if (trainConfig == 2014){ // EMCAL+DCAL clusters standard cuts, triggers, NL vars
    cuts.AddCutCalo("80010123","411790005f032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("8008e123","411790005f032230000","01631031000000d0"); // EG2
    cuts.AddCutCalo("8008d123","411790005f032230000","01631031000000d0"); // EG1
  } else if (trainConfig == 2015){ // EMCAL+DCAL clusters standard cuts, triggers, NL vars
    cuts.AddCutCalo("80010123","411796505f032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("8008e123","411796505f032230000","01631031000000d0"); // EG2
    cuts.AddCutCalo("8008d123","411796505f032230000","01631031000000d0"); // EG1
  } else if (trainConfig == 2016){ // EMCAL+DCAL clusters standard cuts, triggers, NL vars
    cuts.AddCutCalo("80010123","4117957057032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("8008e123","4117957057032230000","01631031000000d0"); // EG2
    cuts.AddCutCalo("8008d123","4117957057032230000","01631031000000d0"); // EG1
  } else if (trainConfig == 2017){ // EMCAL+DCAL clusters standard cuts, triggers, NL vars
    cuts.AddCutCalo("80010123","411795705f032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("8008e123","411795705f032230000","01631031000000d0"); // EG2
    cuts.AddCutCalo("8008d123","411795705f032230000","01631031000000d0"); // EG1
  } else if (trainConfig == 2018){ // EMCAL+DCAL clusters standard cuts, triggers, NL vars
    cuts.AddCutCalo("80010123","411790105f032230000","01631031000000d0"); // INT7
  } else if (trainConfig == 2019){ // EMCAL+DCAL clusters standard cuts, triggers, NL vars
    cuts.AddCutCalo("8008e123","411790105f032230000","01631031000000d0"); // EG2
    cuts.AddCutCalo("8008d123","411790105f032230000","01631031000000d0"); // EG1

  } else if (trainConfig == 2020){ // EMCAL+DCAL clusters standard cuts, triggers, NL vars
    cuts.AddCutCalo("80010123","411793105f032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("8008e123","411793105f032230000","01631031000000d0"); // EG2
    cuts.AddCutCalo("8008d123","411793105f032230000","01631031000000d0"); // EG1
  } else if (trainConfig == 2021){ // EMCAL+DCAL clusters standard cuts, triggers, NL vars
    cuts.AddCutCalo("80010123","411793205f032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("8008e123","411793205f032230000","01631031000000d0"); // EG2
    cuts.AddCutCalo("8008d123","411793205f032230000","01631031000000d0"); // EG1
  } else if (trainConfig == 2022){ // EMCAL+DCAL clusters standard cuts, triggers, NL vars
    cuts.AddCutCalo("80010123","411793305f032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("8008e123","411793305f032230000","01631031000000d0"); // EG2
    cuts.AddCutCalo("8008d123","411793305f032230000","01631031000000d0"); // EG1
  } else if (trainConfig == 2023){ // EMCAL+DCAL clusters standard cuts, triggers, NL vars
    cuts.AddCutCalo("80010123","411793405f032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("8008e123","411793405f032230000","01631031000000d0"); // EG2
    cuts.AddCutCalo("8008d123","411793405f032230000","01631031000000d0"); // EG1
  } else if (trainConfig == 2024){ // open timing, TB NL
    cuts.AddCutCalo("80010123","411796500f032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("80057123","411796500f032230000","01631031000000d0"); // L0
  } else if (trainConfig == 2025){ // open timing, TB NL
    cuts.AddCutCalo("8008e123","411796500f032230000","01631031000000d0"); // L1 low
    cuts.AddCutCalo("8008d123","411796500f032230000","01631031000000d0"); // L1 high
  } else if (trainConfig == 2026){ // open timing, no NL
    cuts.AddCutCalo("80010123","411790000f032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("80057123","411790000f032230000","01631031000000d0"); // L0
  } else if (trainConfig == 2027){ // open timing, TB NL
    cuts.AddCutCalo("8008e123","411790000f032230000","01631031000000d0"); // L1 low
    cuts.AddCutCalo("8008d123","411790000f032230000","01631031000000d0"); // L1 high

  // no TM
  } else if (trainConfig == 2030){ // EMCAL+DCAL clusters standard cuts, triggers, NL vars
    cuts.AddCutCalo("80010123","4117947050032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("8008e123","4117947050032230000","01631031000000d0"); // EG2
    cuts.AddCutCalo("8008d123","4117947050032230000","01631031000000d0"); // EG1
  } else if (trainConfig == 2031){ // EMCAL+DCAL clusters standard cuts, triggers, NL vars
    cuts.AddCutCalo("80010123","4117948050032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("8008e123","4117948050032230000","01631031000000d0"); // EG2
    cuts.AddCutCalo("8008d123","4117948050032230000","01631031000000d0"); // EG1
  } else if (trainConfig == 2032){ // EMCAL+DCAL clusters standard cuts, triggers, NL vars
    cuts.AddCutCalo("80010123","4117957050032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("8008e123","4117957050032230000","01631031000000d0"); // EG2
    cuts.AddCutCalo("8008d123","4117957050032230000","01631031000000d0"); // EG1
  } else if (trainConfig == 2033){ // EMCAL+DCAL clusters standard cuts, triggers, NL vars
    cuts.AddCutCalo("80010123","4117958050032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("8008e123","4117958050032230000","01631031000000d0"); // EG2
    cuts.AddCutCalo("8008d123","4117958050032230000","01631031000000d0"); // EG1
  } else if (trainConfig == 2034){ // EMCAL+DCAL clusters standard cuts, triggers, NL vars
    cuts.AddCutCalo("80010123","4117965050032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("8008e123","4117965050032230000","01631031000000d0"); // EG2
    cuts.AddCutCalo("8008d123","4117965050032230000","01631031000000d0"); // EG1
  } else if (trainConfig == 2035){ // EMCAL+DCAL clusters standard cuts, triggers, NL vars
    cuts.AddCutCalo("80010123","4117900050032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("8008e123","4117900050032230000","01631031000000d0"); // EG2
    cuts.AddCutCalo("8008d123","4117900050032230000","01631031000000d0"); // EG1

  } else if (trainConfig == 2050){ // TB NL tests 100 MeV aggregation
    cuts.AddCutCalo("80010123","411790105f032230000","01631031000000d0"); // INT7
  } else if (trainConfig == 2051){ // TB NL tests 100 MeV aggregation
    cuts.AddCutCalo("8008e123","411790105f032230000","01631031000000d0"); // EG2
    cuts.AddCutCalo("8008d123","411790105f032230000","01631031000000d0"); // EG1
  } else if (trainConfig == 2052){ // TB NL tests 50 MeV aggregation
    cuts.AddCutCalo("80010123","411790205f032230000","01631031000000d0"); // INT7
  } else if (trainConfig == 2053){ // TB NL tests 50 MeV aggregation
    cuts.AddCutCalo("8008e123","411790205f032230000","01631031000000d0"); // EG2
    cuts.AddCutCalo("8008d123","411790205f032230000","01631031000000d0"); // EG1
  } else if (trainConfig == 2054){ // TB NL tests 150 MeV aggregation
    cuts.AddCutCalo("80010123","411790305f032230000","01631031000000d0"); // INT7
  } else if (trainConfig == 2055){ // TB NL tests 150 MeV aggregation
    cuts.AddCutCalo("8008e123","411790305f032230000","01631031000000d0"); // EG2
    cuts.AddCutCalo("8008d123","411790305f032230000","01631031000000d0"); // EG1
  } else if (trainConfig == 2056){ // TB NL tests 300 MeV aggregation
    cuts.AddCutCalo("80010123","411790405f032230000","01631031000000d0"); // INT7
  } else if (trainConfig == 2057){ // TB NL tests 300 MeV aggregation
    cuts.AddCutCalo("8008e123","411790405f032230000","01631031000000d0"); // EG2
    cuts.AddCutCalo("8008d123","411790405f032230000","01631031000000d0"); // EG1

  // standard cut configs with no NL for SM-wise correction
  } else if (trainConfig == 2100){ // EMCAL+DCAL clusters standard cuts, triggers, no NL
    cuts.AddCutCalo("80010123","4117900057032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("8008e123","4117900057032230000","01631031000000d0"); // EG2+DG2
    cuts.AddCutCalo("8008d123","4117900057032230000","01631031000000d0"); // EG1+DG1
  } else if (trainConfig == 2101){ // EMCAL+DCAL clusters standard cuts, triggers, no NL
    cuts.AddCutCalo("80010123","4117900057032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("8008e123","4117900057032230000","01631031000000d0"); // EG2+DG2
    cuts.AddCutCalo("8008d123","4117900057032230000","01631031000000d0"); // EG1+DG1
  } else if (trainConfig == 2102){ // EMCAL+DCAL clusters standard cuts, triggers, no NL
    cuts.AddCutCalo("80010123","4117900057032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("8008e123","4117900057032230000","01631031000000d0"); // EG2+DG2
    cuts.AddCutCalo("8008d123","4117900057032230000","01631031000000d0"); // EG1+DG1
  } else if (trainConfig == 2103){ // EMCAL+DCAL clusters standard cuts, triggers, no NL
    cuts.AddCutCalo("80010123","4117900057032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("8008e123","4117900057032230000","01631031000000d0"); // EG2+DG2
    cuts.AddCutCalo("8008d123","4117900057032230000","01631031000000d0"); // EG1+DG1
  } else if (trainConfig == 2104){ // EMCAL+DCAL clusters standard cuts, triggers, no NL
    cuts.AddCutCalo("80010123","4117900057032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("8008e123","4117900057032230000","01631031000000d0"); // EG2+DG2
    cuts.AddCutCalo("8008d123","4117900057032230000","01631031000000d0"); // EG1+DG1
  } else if (trainConfig == 2105){ // EMCAL+DCAL clusters standard cuts, triggers, no NL
    cuts.AddCutCalo("80010123","4117900057032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("8008e123","4117900057032230000","01631031000000d0"); // EG2+DG2
    cuts.AddCutCalo("8008d123","4117900057032230000","01631031000000d0"); // EG1+DG1
  } else if (trainConfig == 2106){ // EMCAL+DCAL clusters standard cuts, triggers, no NL
    cuts.AddCutCalo("80010123","4117900057032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("8008e123","4117900057032230000","01631031000000d0"); // EG2+DG2
    cuts.AddCutCalo("8008d123","4117900057032230000","01631031000000d0"); // EG1+DG1
  } else if (trainConfig == 2107){ // EMCAL+DCAL clusters standard cuts, triggers, no NL
    cuts.AddCutCalo("80010123","4117900057032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("8008e123","4117900057032230000","01631031000000d0"); // EG2+DG2
    cuts.AddCutCalo("8008d123","4117900057032230000","01631031000000d0"); // EG1+DG1
  } else if (trainConfig == 2108){ // EMCAL+DCAL clusters standard cuts, triggers, no NL
    cuts.AddCutCalo("80010123","4117900057032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("8008e123","4117900057032230000","01631031000000d0"); // EG2+DG2
    cuts.AddCutCalo("8008d123","4117900057032230000","01631031000000d0"); // EG1+DG1
  } else if (trainConfig == 2109){ // EMCAL+DCAL clusters standard cuts, triggers, no NL
    cuts.AddCutCalo("80010123","4117900057032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("8008e123","4117900057032230000","01631031000000d0"); // EG2+DG2
    cuts.AddCutCalo("8008d123","4117900057032230000","01631031000000d0"); // EG1+DG1
  } else if (trainConfig == 2110){ // EMCAL+DCAL clusters standard cuts, triggers, no NL
    cuts.AddCutCalo("80010123","4117900057032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("8008e123","4117900057032230000","01631031000000d0"); // EG2+DG2
    cuts.AddCutCalo("8008d123","4117900057032230000","01631031000000d0"); // EG1+DG1
  } else if (trainConfig == 2111){ // EMCAL+DCAL clusters standard cuts, triggers, no NL
    cuts.AddCutCalo("80010123","4117900057032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("8008e123","4117900057032230000","01631031000000d0"); // EG2+DG2
    cuts.AddCutCalo("8008d123","4117900057032230000","01631031000000d0"); // EG1+DG1
  } else if (trainConfig == 2112){ // EMCAL+DCAL clusters standard cuts, triggers, no NL
    cuts.AddCutCalo("80010123","4117900057032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("8008e123","4117900057032230000","01631031000000d0"); // EG2+DG2
    cuts.AddCutCalo("8008d123","4117900057032230000","01631031000000d0"); // EG1+DG1
  } else if (trainConfig == 2113){ // EMCAL+DCAL clusters standard cuts, triggers, no NL
    cuts.AddCutCalo("80010123","4117900057032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("8008e123","4117900057032230000","01631031000000d0"); // EG2+DG2
    cuts.AddCutCalo("8008d123","4117900057032230000","01631031000000d0"); // EG1+DG1
  } else if (trainConfig == 2114){ // EMCAL+DCAL clusters standard cuts, triggers, no NL
    cuts.AddCutCalo("80010123","4117900057032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("8008e123","4117900057032230000","01631031000000d0"); // EG2+DG2
    cuts.AddCutCalo("8008d123","4117900057032230000","01631031000000d0"); // EG1+DG1
  } else if (trainConfig == 2115){ // EMCAL+DCAL clusters standard cuts, triggers, no NL
    cuts.AddCutCalo("80010123","4117900057032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("8008e123","4117900057032230000","01631031000000d0"); // EG2+DG2
    cuts.AddCutCalo("8008d123","4117900057032230000","01631031000000d0"); // EG1+DG1
  } else if (trainConfig == 2116){ // EMCAL+DCAL clusters standard cuts, triggers, no NL
    cuts.AddCutCalo("80010123","4117900057032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("8008e123","4117900057032230000","01631031000000d0"); // EG2+DG2
    cuts.AddCutCalo("8008d123","4117900057032230000","01631031000000d0"); // EG1+DG1
  } else if (trainConfig == 2117){ // EMCAL+DCAL clusters standard cuts, triggers, no NL
    cuts.AddCutCalo("80010123","4117900057032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("8008e123","4117900057032230000","01631031000000d0"); // EG2+DG2
    cuts.AddCutCalo("8008d123","4117900057032230000","01631031000000d0"); // EG1+DG1
  } else if (trainConfig == 2118){ // EMCAL+DCAL clusters standard cuts, triggers, no NL
    cuts.AddCutCalo("80010123","4117900057032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("8008e123","4117900057032230000","01631031000000d0"); // EG2+DG2
    cuts.AddCutCalo("8008d123","4117900057032230000","01631031000000d0"); // EG1+DG1
  } else if (trainConfig == 2119){ // EMCAL+DCAL clusters standard cuts, triggers, no NL
    cuts.AddCutCalo("80010123","4117900057032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("8008e123","4117900057032230000","01631031000000d0"); // EG2+DG2
    cuts.AddCutCalo("8008d123","4117900057032230000","01631031000000d0"); // EG1+DG1

  // systematics for pPb8TeV PRL
  } else if (trainConfig == 2200) { // CALO variations
    cuts.AddCutCalo("80010123","411793105f032230000","01631031000000d0"); // std
    cuts.AddCutCalo("80010123","411793105f022230000","01631031000000d0"); // 0.6 GeV/c
    cuts.AddCutCalo("80010123","411793105f042230000","01631031000000d0"); // 0.8 GeV/c
  } else if (trainConfig == 2201) {
    cuts.AddCutCalo("80010123","411793105f031230000","01631031000000d0"); // n cells >= 1
    cuts.AddCutCalo("80010123","411793105f033230000","01631031000000d0"); // n cells >= 3
    cuts.AddCutCalo("80010123","411793105f032200000","01631031000000d0"); // no max M02 cut
    cuts.AddCutCalo("80010123","411793105f032250000","01631031000000d0"); // M02 < 0.3
    cuts.AddCutCalo("80010123","411793105f0322k0000","01631031000000d0"); // M02, pT-dep
  } else if (trainConfig == 2202) {
    cuts.AddCutCalo("80010123","411793105e032230000","01631031000000d0"); // TM var e/p 2.0
    cuts.AddCutCalo("80010123","411793105g032230000","01631031000000d0"); // TM var e/p 1.5
    cuts.AddCutCalo("80010123","411793105h032230000","01631031000000d0"); // TM var e/p 1.25
    cuts.AddCutCalo("80010123","4117931057032230000","01631031000000d0"); // TM var no veto
  } else if (trainConfig == 2203) {
    cuts.AddCutCalo("80010123","411793205f032230000","01631031000000d0"); // NL 32
    cuts.AddCutCalo("80010123","411793305f032230000","01631031000000d0"); // NL 33
    cuts.AddCutCalo("80010123","411793405f032230000","01631031000000d0"); // NL 34
  } else if (trainConfig == 2204) {
    cuts.AddCutCalo("80010123","411793106f032230000","01631031000000d0"); // 30/35ns
    cuts.AddCutCalo("80010123","411793104f032230000","01631031000000d0"); // 100ns
    cuts.AddCutCalo("80010123","411793107f032230000","01631031000000d0"); // 30ns

  } else if (trainConfig == 2210) { // CALO variations
    cuts.AddCutCalo("8008e123","411793105f032230000","01631031000000d0"); // std
    cuts.AddCutCalo("8008e123","411793105f022230000","01631031000000d0"); // 0.6 GeV/c
    cuts.AddCutCalo("8008e123","411793105f042230000","01631031000000d0"); // 0.8 GeV/c
  } else if (trainConfig == 2211) {
    cuts.AddCutCalo("8008e123","411793105f031230000","01631031000000d0"); // n cells >= 1
    cuts.AddCutCalo("8008e123","411793105f033230000","01631031000000d0"); // n cells >= 3
    cuts.AddCutCalo("8008e123","411793105f032200000","01631031000000d0"); // no max M02 cut
    cuts.AddCutCalo("8008e123","411793105f032250000","01631031000000d0"); // M02 < 0.3
    cuts.AddCutCalo("8008e123","411793105f0322k0000","01631031000000d0"); // M02, pT-dep
  } else if (trainConfig == 2212) {
    cuts.AddCutCalo("8008e123","411793105e032230000","01631031000000d0"); // TM var e/p 2.0
    cuts.AddCutCalo("8008e123","411793105g032230000","01631031000000d0"); // TM var e/p 1.5
    cuts.AddCutCalo("8008e123","411793105h032230000","01631031000000d0"); // TM var e/p 1.25
    cuts.AddCutCalo("8008e123","4117931057032230000","01631031000000d0"); // TM var no veto
  } else if (trainConfig == 2213) {
    cuts.AddCutCalo("8008e123","411793205f032230000","01631031000000d0"); // NL 32
    cuts.AddCutCalo("8008e123","411793305f032230000","01631031000000d0"); // NL 33
    cuts.AddCutCalo("8008e123","411793405f032230000","01631031000000d0"); // NL 34
  } else if (trainConfig == 2214) {
    cuts.AddCutCalo("8008e123","411793106f032230000","01631031000000d0"); // 30/35ns
    cuts.AddCutCalo("8008e123","411793104f032230000","01631031000000d0"); // 100ns
    cuts.AddCutCalo("8008e123","411793107f032230000","01631031000000d0"); // 30ns

  } else if (trainConfig == 2220) { // CALO variations
    cuts.AddCutCalo("8008d123","411793105f032230000","01631031000000d0"); // std
    cuts.AddCutCalo("8008d123","411793105f022230000","01631031000000d0"); // 0.6 GeV/c
    cuts.AddCutCalo("8008d123","411793105f042230000","01631031000000d0"); // 0.8 GeV/c
  } else if (trainConfig == 2221) {
    cuts.AddCutCalo("8008d123","411793105f031230000","01631031000000d0"); // n cells >= 1
    cuts.AddCutCalo("8008d123","411793105f033230000","01631031000000d0"); // n cells >= 3
    cuts.AddCutCalo("8008d123","411793105f032200000","01631031000000d0"); // no max M02 cut
    cuts.AddCutCalo("8008d123","411793105f032250000","01631031000000d0"); // M02 < 0.3
    cuts.AddCutCalo("8008d123","411793105f0322k0000","01631031000000d0"); // M02, pT-dep
  } else if (trainConfig == 2222) {
    cuts.AddCutCalo("8008d123","411793105e032230000","01631031000000d0"); // TM var e/p 2.0
    cuts.AddCutCalo("8008d123","411793105g032230000","01631031000000d0"); // TM var e/p 1.5
    cuts.AddCutCalo("8008d123","411793105h032230000","01631031000000d0"); // TM var e/p 1.25
    cuts.AddCutCalo("8008d123","4117931057032230000","01631031000000d0"); // TM var no veto
  } else if (trainConfig == 2223) {
    cuts.AddCutCalo("8008d123","411793205f032230000","01631031000000d0"); // NL 32
    cuts.AddCutCalo("8008d123","411793305f032230000","01631031000000d0"); // NL 33
    cuts.AddCutCalo("8008d123","411793405f032230000","01631031000000d0"); // NL 34
  } else if (trainConfig == 2224) {
    cuts.AddCutCalo("8008d123","411793106f032230000","01631031000000d0"); // 30/35ns
    cuts.AddCutCalo("8008d123","411793104f032230000","01631031000000d0"); // 100ns
    cuts.AddCutCalo("8008d123","411793107f032230000","01631031000000d0"); // 30ns


  } else {
    Error(Form("GammaCalo_%i",trainConfig), "wrong trainConfig variable no cuts have been specified for the configuration");
    return;
  }

  if(!cuts.AreValid()){
    cout << "\n\n****************************************************" << endl;
    cout << "ERROR: No valid cuts stored in CutHandlerCalo! Returning..." << endl;
    cout << "****************************************************\n\n" << endl;
    return;
  }

  Int_t numberOfCuts = cuts.GetNCuts();

  TList *EventCutList     = new TList();
  TList *ClusterCutList   = new TList();
  TList *MesonCutList     = new TList();

  TList *HeaderList = new TList();
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
  ClusterCutList->SetOwner(kTRUE);
  AliCaloPhotonCuts **analysisClusterCuts     = new AliCaloPhotonCuts*[numberOfCuts];
  MesonCutList->SetOwner(kTRUE);
  AliConversionMesonCuts **analysisMesonCuts  = new AliConversionMesonCuts*[numberOfCuts];

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

    if (enableMultiplicityWeighting) analysisEventCuts[i]->SetUseWeightMultiplicityFromFile( kTRUE, fileNameMultWeights, dataInputMultHisto, mcInputMultHisto );

    analysisEventCuts[i]->SetTriggerMimicking(enableTriggerMimicking);
    analysisEventCuts[i]->SetTriggerOverlapRejecion(enableTriggerOverlapRej);
    if(fMinPtHardSet)
      analysisEventCuts[i]->SetMinFacPtHard(minFacPtHard);
    if(fMaxPtHardSet)
      analysisEventCuts[i]->SetMaxFacPtHard(maxFacPtHard);
    if(fSingleMaxPtHardSet)
      analysisEventCuts[i]->SetMaxFacPtHardSingleParticle(maxFacPtHardSingle);
    analysisEventCuts[i]->SetV0ReaderName(V0ReaderName);
    analysisEventCuts[i]->SetCorrectionTaskSetting(corrTaskSetting);
    analysisEventCuts[i]->SetLightOutput(enableLightOutput);
    if (periodNameV0Reader.CompareTo("") != 0) analysisEventCuts[i]->SetPeriodEnum(periodNameV0Reader);
    analysisEventCuts[i]->InitializeCutsFromCutString((cuts.GetEventCut(i)).Data());
    EventCutList->Add(analysisEventCuts[i]);
    analysisEventCuts[i]->SetFillCutHistograms("",kFALSE);

    analysisClusterCuts[i] = new AliCaloPhotonCuts(isMC);
    analysisClusterCuts[i]->SetHistoToModifyAcceptance(histoAcc);
    analysisClusterCuts[i]->SetV0ReaderName(V0ReaderName);
    analysisClusterCuts[i]->SetCorrectionTaskSetting(corrTaskSetting);
    analysisClusterCuts[i]->SetCaloTrackMatcherName(TrackMatcherName);
    analysisClusterCuts[i]->SetLightOutput(enableLightOutput);
    analysisClusterCuts[i]->InitializeCutsFromCutString((cuts.GetClusterCut(i)).Data());
    ClusterCutList->Add(analysisClusterCuts[i]);
    analysisClusterCuts[i]->SetExtendedMatchAndQA(enableExtMatchAndQA);
    analysisClusterCuts[i]->SetFillCutHistograms("");

    analysisMesonCuts[i] = new AliConversionMesonCuts();
    analysisMesonCuts[i]->SetLightOutput(enableLightOutput);
    analysisMesonCuts[i]->InitializeCutsFromCutString((cuts.GetMesonCut(i)).Data());
    analysisMesonCuts[i]->SetIsMergedClusterCut(2);
    analysisMesonCuts[i]->SetCaloMesonCutsObject(analysisClusterCuts[i]);
    MesonCutList->Add(analysisMesonCuts[i]);
    analysisMesonCuts[i]->SetFillCutHistograms("");
    analysisEventCuts[i]->SetAcceptedHeader(HeaderList);
  }
  task->SetEventCutList(numberOfCuts,EventCutList);
  task->SetCaloCutList(numberOfCuts,ClusterCutList);
  task->SetMesonCutList(numberOfCuts,MesonCutList);
  task->SetCorrectionTaskSetting(corrTaskSetting);
  task->SetDoMesonAnalysis(kTRUE);
  task->SetDoMesonQA(enableQAMesonTask); //Attention new switch for Pi0 QA
  task->SetDoClusterQA(enableQAClusterTask);  //Attention new switch small for Cluster QA
  task->SetDoTHnSparse(enableTHnSparse);
  task->SetProduceTreeEOverP(doTreeEOverP);
  task->SetEnableSortingOfMCClusLabels(enableSortingMCLabels);
  if(enableExtMatchAndQA > 1){ task->SetPlotHistsExtQA(kTRUE);}

  //connect containers
  AliAnalysisDataContainer *coutput =
    mgr->CreateContainer(!(corrTaskSetting.CompareTo("")) ? Form("GammaCalo_%i",trainConfig) : Form("GammaCalo_%i_%s",trainConfig,corrTaskSetting.Data()), TList::Class(),
              AliAnalysisManager::kOutputContainer, Form("GammaCalo_%i.root",trainConfig) );

  mgr->AddTask(task);
  mgr->ConnectInput(task,0,cinput);
  mgr->ConnectOutput(task,1,coutput);

  Int_t nContainer = 2;
  for(Int_t i = 0; i<numberOfCuts; i++){
      if(enableQAMesonTask==5){
          mgr->ConnectOutput(task,nContainer,mgr->CreateContainer(Form("%s_%s_%s ClusterTimingEff",(cuts.GetEventCut(i)).Data(),(cuts.GetClusterCut(i)).Data(),(cuts.GetMesonCut(i)).Data()), TTree::Class(), AliAnalysisManager::kOutputContainer, Form("GammaCalo_%i.root",trainConfig)) );
          nContainer++;
      }
  }

  return;

}
