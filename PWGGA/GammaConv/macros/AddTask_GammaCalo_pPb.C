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
  Int_t     enableTriggerMimicking        = 0,        // enable trigger mimicking
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
  TString fileNameCustomTriggerMimicOADB   = cuts.GetSpecialFileNameFromString (fileNameExternalInputs, "FTRM:");

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
  Bool_t fJetFinderUsage  = kFALSE;
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
  // EMC +(DMC) clusters pPb 5TeV
  // ===============================================================================================
  if (trainConfig == 1){
    cuts.AddCutCalo("80010113","411790105f032230000","01631031000000d0"); // 0-100
    cuts.AddCutCalo("80010113","111110105f032230000","01631031000000d0"); // 0-100
  } else if (trainConfig == 2){
    cuts.AddCutCalo("80010113","411793105f032230000","01631031000000d0"); // 0-100
    cuts.AddCutCalo("80010113","111113105f032230000","01631031000000d0"); // 0-100
  } else if (trainConfig == 3){ // no time cut, TB NL
    cuts.AddCutCalo("80010113","411793100f032230000","01631031000000d0"); // 0-100
    cuts.AddCutCalo("80010113","111113100f032230000","01631031000000d0"); // 0-100
  } else if (trainConfig == 4){
    cuts.AddCutCalo("80010123","411793105f032230000","01631031000000d0"); // 0-100
    cuts.AddCutCalo("80010123","111113105f032230000","01631031000000d0"); // 0-100
  } else if (trainConfig == 5){
    cuts.AddCutCalo("80010113","411793105f032230000","01631031000000d0"); // 0-100
    cuts.AddCutCalo("80010113","411793105f032230000","01631031000000i0"); // 0-100 only pair in same calo
    cuts.AddCutCalo("80010113","411793105f032230300","01631031000000i0"); // 0-100 remove 1 leg of conv cand Mgg < 0.03
    cuts.AddCutCalo("80010113","411793105f032230200","01631031000000i0"); // 0-100 remove 1 leg of conv cand Mgg < 0.025
    cuts.AddCutCalo("80010113","411793105f032230400","01631031000000i0"); // 0-100 remove 1 leg of conv cand Mgg < 0.045
  } else if (trainConfig == 6){
    cuts.AddCutCalo("80010113","111113105f032230000","01631031000000d0"); // 0-100
    cuts.AddCutCalo("80010113","111113105f032230000","01631031000000i0"); // 0-100 only pair in same calo
    cuts.AddCutCalo("80010113","111113105f032230300","01631031000000i0"); // 0-100 remove 1 leg of conv cand Mgg < 0.03
    cuts.AddCutCalo("80010113","111113105f032230200","01631031000000i0"); // 0-100 remove 1 leg of conv cand Mgg < 0.025
    cuts.AddCutCalo("80010113","111113105f032230400","01631031000000i0"); // 0-100 remove 1 leg of conv cand Mgg < 0.035
  } else if (trainConfig == 7){
    cuts.AddCutCalo("80010123","111113105f032230000","01631031000000d0"); // 0-100
    cuts.AddCutCalo("80010123","111113105f032230000","01631031000000i0"); // 0-100 only pair in same calo
    cuts.AddCutCalo("80010123","111113105f032230300","01631031000000i0"); // 0-100 remove 1 leg of conv cand Mgg < 0.03
    cuts.AddCutCalo("80010123","111113105f032230200","01631031000000i0"); // 0-100 remove 1 leg of conv cand Mgg < 0.025
    cuts.AddCutCalo("80010123","111113105f032230400","01631031000000i0"); // 0-100 remove 1 leg of conv cand Mgg < 0.035


  // CONFIGS for sys MB MCs & JJs
  } else if (trainConfig == 10) { // CALO variations
    cuts.AddCutCalo("80010113","411793105f032230000","01631031000000d0"); // std
    cuts.AddCutCalo("80010113","411793105f022230000","01631031000000d0"); // 0.6 GeV/c
    cuts.AddCutCalo("80010113","411793105f042230000","01631031000000d0"); // 0.8 GeV/c
  } else if (trainConfig == 11) {
    cuts.AddCutCalo("80010113","411793105f031230000","01631031000000d0"); // n cells >= 1
    cuts.AddCutCalo("80010113","411793105f033230000","01631031000000d0"); // n cells >= 3
    cuts.AddCutCalo("80010113","411793105f032200000","01631031000000d0"); // no max M02 cut
    cuts.AddCutCalo("80010113","411793105f032250000","01631031000000d0"); // M02 < 0.3
    cuts.AddCutCalo("80010113","411793105f0322k0000","01631031000000d0"); // M02, pT-dep
  } else if (trainConfig == 12) {
    cuts.AddCutCalo("80010113","411793105e032230000","01631031000000d0"); // TM var e/p 2.0
    cuts.AddCutCalo("80010113","411793105g032230000","01631031000000d0"); // TM var e/p 1.5
    cuts.AddCutCalo("80010113","411793105h032230000","01631031000000d0"); // TM var e/p 1.25
    cuts.AddCutCalo("80010113","4117931057032230000","01631031000000d0"); // TM var no veto
  } else if (trainConfig == 13) {
    cuts.AddCutCalo("80010113","411793205f032230000","01631031000000d0"); // NL 32
    cuts.AddCutCalo("80010113","411793305f032230000","01631031000000d0"); // NL 33
    cuts.AddCutCalo("80010113","411793405f032230000","01631031000000d0"); // NL 34
    cuts.AddCutCalo("80010113","411793805f032230000","01631031000000d0"); // NL 38
    cuts.AddCutCalo("80010113","411793905f032230000","01631031000000d0"); // NL 39
  } else if (trainConfig == 14) {
    cuts.AddCutCalo("80010113","411793106f032230000","01631031000000d0"); // 30/35ns
    cuts.AddCutCalo("80010113","411793104f032230000","01631031000000d0"); // 100ns
    cuts.AddCutCalo("80010113","411793107f032230000","01631031000000d0"); // 30ns
  } else if (trainConfig == 15) { // CALO variations
    cuts.AddCutCalo("80010113","111113105f032230000","01631031000000d0"); // std
    cuts.AddCutCalo("80010113","111113105f022230000","01631031000000d0"); // 0.6 GeV/c
    cuts.AddCutCalo("80010113","111113105f042230000","01631031000000d0"); // 0.8 GeV/c
  } else if (trainConfig == 16) {
    cuts.AddCutCalo("80010113","111113105f031230000","01631031000000d0"); // n cells >= 1
    cuts.AddCutCalo("80010113","111113105f033230000","01631031000000d0"); // n cells >= 3
    cuts.AddCutCalo("80010113","111113105f032200000","01631031000000d0"); // no max M02 cut
    cuts.AddCutCalo("80010113","111113105f032250000","01631031000000d0"); // M02 < 0.3
    cuts.AddCutCalo("80010113","111113105f0322k0000","01631031000000d0"); // M02, pT-dep
  } else if (trainConfig == 17) {
    cuts.AddCutCalo("80010113","111113105e032230000","01631031000000d0"); // TM var e/p 2.0
    cuts.AddCutCalo("80010113","111113105g032230000","01631031000000d0"); // TM var e/p 1.5
    cuts.AddCutCalo("80010113","111113105h032230000","01631031000000d0"); // TM var e/p 1.25
    cuts.AddCutCalo("80010113","1111131057032230000","01631031000000d0"); // TM var no veto
  } else if (trainConfig == 18) {
    cuts.AddCutCalo("80010113","111113205f032230000","01631031000000d0"); // NL 32
    cuts.AddCutCalo("80010113","111113305f032230000","01631031000000d0"); // NL 33
    cuts.AddCutCalo("80010113","111113405f032230000","01631031000000d0"); // NL 34
    cuts.AddCutCalo("80010113","111113805f032230000","01631031000000d0"); // NL 34
    cuts.AddCutCalo("80010113","111113905f032230000","01631031000000d0"); // NL 34
  } else if (trainConfig == 19) {
    cuts.AddCutCalo("80010113","111113106f032230000","01631031000000d0"); // 30/35ns
    cuts.AddCutCalo("80010113","111113104f032230000","01631031000000d0"); // 100ns
    cuts.AddCutCalo("80010113","111113107f032230000","01631031000000d0"); // 30ns

  // CONFIGS for sys JJ MC with MB BG + JJ
  } else if (trainConfig == 20) { // CALO variations
    cuts.AddCutCalo("80010123","111113105f032230000","01631031000000d0"); // std
    cuts.AddCutCalo("80010123","111113105f022230000","01631031000000d0"); // 0.6 GeV/c
    cuts.AddCutCalo("80010123","111113105f042230000","01631031000000d0"); // 0.8 GeV/c
  } else if (trainConfig == 21) {
    cuts.AddCutCalo("80010123","111113105f031230000","01631031000000d0"); // n cells >= 1
    cuts.AddCutCalo("80010123","111113105f033230000","01631031000000d0"); // n cells >= 3
    cuts.AddCutCalo("80010123","111113105f032200000","01631031000000d0"); // no max M02 cut
    cuts.AddCutCalo("80010123","111113105f032250000","01631031000000d0"); // M02 < 0.3
    cuts.AddCutCalo("80010123","111113105f0322k0000","01631031000000d0"); // M02, pT-dep
  } else if (trainConfig == 22) {
    cuts.AddCutCalo("80010123","111113105e032230000","01631031000000d0"); // TM var e/p 2.0
    cuts.AddCutCalo("80010123","111113105g032230000","01631031000000d0"); // TM var e/p 1.5
    cuts.AddCutCalo("80010123","111113105h032230000","01631031000000d0"); // TM var e/p 1.25
    cuts.AddCutCalo("80010123","1111131057032230000","01631031000000d0"); // TM var no veto
  } else if (trainConfig == 23) {
    cuts.AddCutCalo("80010123","111113205f032230000","01631031000000d0"); // NL 32
    cuts.AddCutCalo("80010123","111113305f032230000","01631031000000d0"); // NL 33
    cuts.AddCutCalo("80010123","111113405f032230000","01631031000000d0"); // NL 34
    cuts.AddCutCalo("80010123","111113805f032230000","01631031000000d0"); // NL 34
    cuts.AddCutCalo("80010123","111113905f032230000","01631031000000d0"); // NL 34
  } else if (trainConfig == 24) {
    cuts.AddCutCalo("80010123","111113106f032230000","01631031000000d0"); // 30/35ns
    cuts.AddCutCalo("80010123","111113104f032230000","01631031000000d0"); // 100ns
    cuts.AddCutCalo("80010123","111113107f032230000","01631031000000d0"); // 30ns
  } else if (trainConfig == 25){ // opening angle variations
    cuts.AddCutCalo("80010113","111113105f032230000","01631031000000b0"); // min opening angle 0.0152
    cuts.AddCutCalo("80010113","111113105f032230000","01631031000000c0"); // min opening angle 0.016
    cuts.AddCutCalo("80010113","111113105f032230000","01631031000000e0"); // min opening angle 0.018
    cuts.AddCutCalo("80010113","111113105f032230000","01631031000000f0"); // min opening angle 0.019
  } else if (trainConfig == 26){
    cuts.AddCutCalo("80010113","111113105f0322l0000","01631031000000d0"); // M02 pt dep with new std: 0.32, 0.0072, 0.5
    cuts.AddCutCalo("80010113","111113105f0322d0000","01631031000000d0"); // M02, pt dep with  0.27, 0.0072, 0.4
    cuts.AddCutCalo("80010113","111113105f0322e0000","01631031000000d0"); // M02, pt dep with  0.31, 0.0072, 0.5
    cuts.AddCutCalo("80010113","111113105f0322f0000","01631031000000d0"); // M02, pt dep with  0.36, 0.0072, 0.7
    cuts.AddCutCalo("80010113","111113105f0322m0000","01631031000000d0"); // M02, pt dep with  0.32, 0.0152, 0.5
  } else if (trainConfig == 27){ // opening angle variations
    cuts.AddCutCalo("80010123","111113105f032230000","01631031000000b0"); // min opening angle 0.0152
    cuts.AddCutCalo("80010123","111113105f032230000","01631031000000c0"); // min opening angle 0.016
    cuts.AddCutCalo("80010123","111113105f032230000","01631031000000e0"); // min opening angle 0.018
    cuts.AddCutCalo("80010123","111113105f032230000","01631031000000f0"); // min opening angle 0.019
  } else if (trainConfig == 28){
    cuts.AddCutCalo("80010123","111113105f0322l0000","01631031000000d0"); // M02 pt dep with new std: 0.32, 0.0072, 0.5
    cuts.AddCutCalo("80010123","111113105f0322d0000","01631031000000d0"); // M02, pt dep with  0.27, 0.0072, 0.4
    cuts.AddCutCalo("80010123","111113105f0322e0000","01631031000000d0"); // M02, pt dep with  0.31, 0.0072, 0.5
    cuts.AddCutCalo("80010123","111113105f0322f0000","01631031000000d0"); // M02, pt dep with  0.36, 0.0072, 0.7
    cuts.AddCutCalo("80010123","111113105f0322m0000","01631031000000d0"); // M02, pt dep with  0.32, 0.0152, 0.5
  } else if (trainConfig == 29){
    cuts.AddCutCalo("80010113","111113105f032230000","01633031000000d0"); // rapidity variation  y<0.6
    cuts.AddCutCalo("80010113","111113105f032230000","01634031000000d0"); // rapidity variation  y<0.5
    cuts.AddCutCalo("80010113","111113105f032230000","01631061000000d0"); // alpha meson variation 1   0<alpha<0.8
    cuts.AddCutCalo("80010113","111113105f032230000","01631051000000d0"); // alpha meson variation 2  0<alpha<0.75
  } else if (trainConfig == 30){
    cuts.AddCutCalo("80010123","111113105f032230000","01633031000000d0"); // rapidity variation  y<0.6
    cuts.AddCutCalo("80010123","111113105f032230000","01634031000000d0"); // rapidity variation  y<0.5
    cuts.AddCutCalo("80010123","111113105f032230000","01631061000000d0"); // alpha meson variation 1   0<alpha<0.8
    cuts.AddCutCalo("80010123","111113105f032230000","01631051000000d0"); // alpha meson variation 2  0<alpha<0.75
  } else if (trainConfig == 31){ // opening angle variations
    cuts.AddCutCalo("80010113","411793105f032230000","01631031000000b0"); // min opening angle 0.0152
    cuts.AddCutCalo("80010113","411793105f032230000","01631031000000c0"); // min opening angle 0.016
    cuts.AddCutCalo("80010113","411793105f032230000","01631031000000e0"); // min opening angle 0.018
    cuts.AddCutCalo("80010113","411793105f032230000","01631031000000f0"); // min opening angle 0.019
  } else if (trainConfig == 32){
    cuts.AddCutCalo("80010113","411793105f0322l0000","01631031000000d0"); // M02 pt dep with new std: 0.32, 0.0072, 0.5
    cuts.AddCutCalo("80010113","411793105f0322d0000","01631031000000d0"); // M02, pt dep with  0.27, 0.0072, 0.4
    cuts.AddCutCalo("80010113","411793105f0322e0000","01631031000000d0"); // M02, pt dep with  0.31, 0.0072, 0.5
    cuts.AddCutCalo("80010113","411793105f0322f0000","01631031000000d0"); // M02, pt dep with  0.36, 0.0072, 0.7
    cuts.AddCutCalo("80010113","411793105f0322m0000","01631031000000d0"); // M02, pt dep with  0.32, 0.0152, 0.5
  } else if (trainConfig == 33){ // opening angle variations
    cuts.AddCutCalo("80010123","411793105f032230000","01631031000000b0"); // min opening angle 0.0152
    cuts.AddCutCalo("80010123","411793105f032230000","01631031000000c0"); // min opening angle 0.016
    cuts.AddCutCalo("80010123","411793105f032230000","01631031000000e0"); // min opening angle 0.018
    cuts.AddCutCalo("80010123","411793105f032230000","01631031000000f0"); // min opening angle 0.019
  } else if (trainConfig == 34){
    cuts.AddCutCalo("80010123","411793105f0322l0000","01631031000000d0"); // M02 pt dep with new std: 0.32, 0.0072, 0.5
    cuts.AddCutCalo("80010123","411793105f0322d0000","01631031000000d0"); // M02, pt dep with  0.27, 0.0072, 0.4
    cuts.AddCutCalo("80010123","411793105f0322e0000","01631031000000d0"); // M02, pt dep with  0.31, 0.0072, 0.5
    cuts.AddCutCalo("80010123","411793105f0322f0000","01631031000000d0"); // M02, pt dep with  0.36, 0.0072, 0.7
    cuts.AddCutCalo("80010123","411793105f022230000","01631031000000d0"); // M02, pt dep with  0.32, 0.0152, 0.5
  } else if (trainConfig == 35){
    cuts.AddCutCalo("80010113","411793105f032230000","01633031000000d0"); // rapidity variation  y<0.6
    cuts.AddCutCalo("80010113","411793105f032230000","01634031000000d0"); // rapidity variation  y<0.5
    cuts.AddCutCalo("80010113","411793105f032230000","01631061000000d0"); // alpha meson variation 1   0<alpha<0.8
    cuts.AddCutCalo("80010113","411793105f032230000","01631051000000d0"); // alpha meson variation 2  0<alpha<0.75
  } else if (trainConfig == 36){
    cuts.AddCutCalo("80010123","411793105f032230000","01633031000000d0"); // rapidity variation  y<0.6
    cuts.AddCutCalo("80010123","411793105f032230000","01634031000000d0"); // rapidity variation  y<0.5
    cuts.AddCutCalo("80010123","411793105f032230000","01631061000000d0"); // alpha meson variation 1   0<alpha<0.8
    cuts.AddCutCalo("80010123","411793105f032230000","01631051000000d0"); // alpha meson variation 2  0<alpha<0.75  
    
    
  } else if (trainConfig == 100){ // EMCAL clusters standard cuts, cent vars
    cuts.AddCutCalo("a0110113","111113105f032230000","01631031000000d0"); // 0-5
    cuts.AddCutCalo("a1210113","111113105f032230000","01631031000000d0"); // 5-10
    cuts.AddCutCalo("81210113","111113105f032230000","01631031000000d0"); // 10-20
  } else if (trainConfig == 101){ // EMCAL+DCAL clusters standard cuts, cent vars
    cuts.AddCutCalo("a0110113","411793105f032230000","01631031000000d0"); // 0-5
    cuts.AddCutCalo("a1210113","411793105f032230000","01631031000000d0"); // 5-10
    cuts.AddCutCalo("81210113","411793105f032230000","01631031000000d0"); // 10-20
  } else if (trainConfig == 102){ // EMCAL clusters standard cuts, cent vars
    cuts.AddCutCalo("82410113","111113105f032230000","01631031000000d0"); // 20-40
    cuts.AddCutCalo("84610113","111113105f032230000","01631031000000d0"); // 40-60
    cuts.AddCutCalo("86810113","111113105f032230000","01631031000000d0"); // 60-80
    cuts.AddCutCalo("88010113","111113105f032230000","01631031000000d0"); // 80-100
  } else if (trainConfig == 103){ // EMCAL+DCAL clusters standard cuts, cent vars
    cuts.AddCutCalo("82410113","411793105f032230000","01631031000000d0"); // 20-40
    cuts.AddCutCalo("84610113","411793105f032230000","01631031000000d0"); // 40-60
    cuts.AddCutCalo("86810113","411793105f032230000","01631031000000d0"); // 60-80
    cuts.AddCutCalo("88010113","411793105f032230000","01631031000000d0"); // 80-100

  } else if (trainConfig == 190){ // EDC JETS - EDC
    cuts.AddCutCalo("80010113","411793105f032230000","2l631031000000d0"); // Standard EDC Jets
    cuts.AddCutCalo("80010113","411793105f032230000","01631031000000d0"); // Standard EDC MB
  } else if (trainConfig == 191){ // EDC JETS - EDC - enabling jet QA
    cuts.AddCutCalo("80010113","411793105f032230000","3l631031000000d0"); // Standard EDC Jets

  // ===============================================================================================
  // EMC +(DMC) clusters pPb 5TeV triggrs
  // ===============================================================================================

  } else if (trainConfig == 200){ // EMCAL clusters standard cuts
    cuts.AddCutCalo("80052113","111113105f032230000","01631031000000d0"); // 0-100
  } else if (trainConfig == 201){ // no time cut
    cuts.AddCutCalo("80052113","111113100f032230000","01631031000000d0"); // 0-100
  } else if (trainConfig == 202){ // EMCAL clusters standard cuts, cent vars, non lin variations
    cuts.AddCutCalo("a0152113","111113105f032230000","01631031000000d0"); // 0-5
    cuts.AddCutCalo("a1252113","111113105f032230000","01631031000000d0"); // 5-10
    cuts.AddCutCalo("81252113","111113105f032230000","01631031000000d0"); // 10-20
  } else if (trainConfig == 203){ // EMCAL clusters standard cuts, cent vars, non lin variations
    cuts.AddCutCalo("82452113","111113105f032230000","01631031000000d0"); // 20-40
    cuts.AddCutCalo("84652113","111113105f032230000","01631031000000d0"); // 40-60
    cuts.AddCutCalo("86852113","111113105f032230000","01631031000000d0"); // 60-80
    cuts.AddCutCalo("88052113","111113105f032230000","01631031000000d0"); // 80-100

  } else if (trainConfig == 204){
    cuts.AddCutCalo("80085113","111113105f032230000","01631031000000d0"); // 0-100
    cuts.AddCutCalo("80085113","111113105f032230000","01631031000000i0"); // 0-100 only pair in same calo
    cuts.AddCutCalo("80085113","111113105f032230300","01631031000000i0"); // 0-100 remove 1 leg of conv cand Mgg < 0.03
    cuts.AddCutCalo("80085113","111113105f032230200","01631031000000i0"); // 0-100 remove 1 leg of conv cand Mgg < 0.025
    cuts.AddCutCalo("80085113","111113105f032230400","01631031000000i0"); // 0-100 remove 1 leg of conv cand Mgg < 0.035
  } else if (trainConfig == 205){
    cuts.AddCutCalo("80085123","111113105f032230000","01631031000000d0"); // 0-100
    cuts.AddCutCalo("80085123","111113105f032230000","01631031000000i0"); // 0-100 only pair in same calo
    cuts.AddCutCalo("80085123","111113105f032230300","01631031000000i0"); // 0-100 remove 1 leg of conv cand Mgg < 0.03
    cuts.AddCutCalo("80085123","111113105f032230200","01631031000000i0"); // 0-100 remove 1 leg of conv cand Mgg < 0.025
    cuts.AddCutCalo("80085123","111113105f032230400","01631031000000i0"); // 0-100 remove 1 leg of conv cand Mgg < 0.035
  } else if (trainConfig == 206){
    cuts.AddCutCalo("80083113","111113105f032230000","01631031000000d0"); // 0-100
    cuts.AddCutCalo("80083113","111113105f032230000","01631031000000i0"); // 0-100 only pair in same calo
    cuts.AddCutCalo("80083113","111113105f032230300","01631031000000i0"); // 0-100 remove 1 leg of conv cand Mgg < 0.03
    cuts.AddCutCalo("80083113","111113105f032230200","01631031000000i0"); // 0-100 remove 1 leg of conv cand Mgg < 0.025
    cuts.AddCutCalo("80083113","111113105f032230400","01631031000000i0"); // 0-100 remove 1 leg of conv cand Mgg < 0.035
  } else if (trainConfig == 207){
    cuts.AddCutCalo("80083123","111113105f032230000","01631031000000d0"); // 0-100
    cuts.AddCutCalo("80083123","111113105f032230000","01631031000000i0"); // 0-100 only pair in same calo
    cuts.AddCutCalo("80083123","111113105f032230300","01631031000000i0"); // 0-100 remove 1 leg of conv cand Mgg < 0.03
    cuts.AddCutCalo("80083123","111113105f032230200","01631031000000i0"); // 0-100 remove 1 leg of conv cand Mgg < 0.025
    cuts.AddCutCalo("80083123","111113105f032230400","01631031000000i0"); // 0-100 remove 1 leg of conv cand Mgg < 0.035

  } else if (trainConfig == 220){ // EMCAL clusters standard cuts
    cuts.AddCutCalo("80085113","111113105f032230000","01631031000000d0"); // 0-100
  } else if (trainConfig == 221){ // no time cut
    cuts.AddCutCalo("80085113","111113100f032230000","01631031000000d0"); // 0-100
  } else if (trainConfig == 222){ // EMCAL clusters standard cuts, cent vars, non lin variations
    cuts.AddCutCalo("a0185113","111113105f032230000","01631031000000d0"); // 0-5
    cuts.AddCutCalo("a1285113","111113105f032230000","01631031000000d0"); // 5-10
    cuts.AddCutCalo("81285113","111113105f032230000","01631031000000d0"); // 10-20
  } else if (trainConfig == 223){ // EMCAL clusters standard cuts, cent vars, non lin variations
    cuts.AddCutCalo("82485113","111113105f032230000","01631031000000d0"); // 20-40
    cuts.AddCutCalo("84685113","111113105f032230000","01631031000000d0"); // 40-60
    cuts.AddCutCalo("86885113","111113105f032230000","01631031000000d0"); // 60-80
    cuts.AddCutCalo("88085113","111113105f032230000","01631031000000d0"); // 80-100
  } else if (trainConfig == 224){
    cuts.AddCutCalo("80085123","111113105f032230000","01631031000000d0"); // 0-100
    // CONFIGS for sys JJ MC or MB MC
  } else if (trainConfig == 225) { // CALO variations
    cuts.AddCutCalo("80085113","111113105f032230000","01631031000000d0"); // std
    cuts.AddCutCalo("80085113","111113105f022230000","01631031000000d0"); // 0.6 GeV/c
    cuts.AddCutCalo("80085113","111113105f042230000","01631031000000d0"); // 0.8 GeV/c
  } else if (trainConfig == 226) {
    cuts.AddCutCalo("80085113","111113105f031230000","01631031000000d0"); // n cells >= 1
    cuts.AddCutCalo("80085113","111113105f033230000","01631031000000d0"); // n cells >= 3
    cuts.AddCutCalo("80085113","111113105f032200000","01631031000000d0"); // no max M02 cut
    cuts.AddCutCalo("80085113","111113105f032250000","01631031000000d0"); // M02 < 0.3
    cuts.AddCutCalo("80085113","111113105f0322k0000","01631031000000d0"); // M02, pT-dep
  } else if (trainConfig == 227) {
    cuts.AddCutCalo("80085113","111113105e032230000","01631031000000d0"); // TM var e/p 2.0
    cuts.AddCutCalo("80085113","111113105g032230000","01631031000000d0"); // TM var e/p 1.5
    cuts.AddCutCalo("80085113","111113105h032230000","01631031000000d0"); // TM var e/p 1.25
    cuts.AddCutCalo("80085113","1111131057032230000","01631031000000d0"); // TM var no veto
    cuts.AddCutCalo("80085113","111113105f032230000","01631061000000d0"); // alpha meson variation 1   0<alpha<0.8
    cuts.AddCutCalo("80085113","111113105f032230000","01631051000000d0"); // alpha meson variation 2  0<alpha<0.75
  } else if (trainConfig == 228) {
    cuts.AddCutCalo("80085113","111113205f032230000","01631031000000d0"); // NL 32
    cuts.AddCutCalo("80085113","111113305f032230000","01631031000000d0"); // NL 33
    cuts.AddCutCalo("80085113","111113405f032230000","01631031000000d0"); // NL 34
    cuts.AddCutCalo("80085113","111113805f032230000","01631031000000d0"); // NL 38
    cuts.AddCutCalo("80085113","111113905f032230000","01631031000000d0"); // NL 39
  } else if (trainConfig == 229) {
    cuts.AddCutCalo("80085113","111113106f032230000","01631031000000d0"); // 30/35ns
    cuts.AddCutCalo("80085113","111113104f032230000","01631031000000d0"); // 100ns
    cuts.AddCutCalo("80085113","111113107f032230000","01631031000000d0"); // 30ns
  // CONFIGS for sys JJ MC with MB BG + JJ
  } else if (trainConfig == 230) { // CALO variations
    cuts.AddCutCalo("80085123","111113105f032230000","01631031000000d0"); // std
    cuts.AddCutCalo("80085123","111113105f022230000","01631031000000d0"); // 0.6 GeV/c
    cuts.AddCutCalo("80085123","111113105f042230000","01631031000000d0"); // 0.8 GeV/c
  } else if (trainConfig == 231) {
    cuts.AddCutCalo("80085123","111113105f031230000","01631031000000d0"); // n cells >= 1
    cuts.AddCutCalo("80085123","111113105f033230000","01631031000000d0"); // n cells >= 3
    cuts.AddCutCalo("80085123","111113105f032200000","01631031000000d0"); // no max M02 cut
    cuts.AddCutCalo("80085123","111113105f032250000","01631031000000d0"); // M02 < 0.3
    cuts.AddCutCalo("80085123","111113105f0322k0000","01631031000000d0"); // M02, pT-dep
  } else if (trainConfig == 232) {
    cuts.AddCutCalo("80085123","111113105e032230000","01631031000000d0"); // TM var e/p 2.0
    cuts.AddCutCalo("80085123","111113105g032230000","01631031000000d0"); // TM var e/p 1.5
    cuts.AddCutCalo("80085123","111113105h032230000","01631031000000d0"); // TM var e/p 1.25
    cuts.AddCutCalo("80085123","1111131057032230000","01631031000000d0"); // TM var no veto
    cuts.AddCutCalo("80085123","111113105f032230000","01631061000000d0"); // alpha meson variation 1   0<alpha<0.8
    cuts.AddCutCalo("80085123","111113105f032230000","01631051000000d0"); // alpha meson variation 2  0<alpha<0.75
  } else if (trainConfig == 233) {
    cuts.AddCutCalo("80085123","111113205f032230000","01631031000000d0"); // NL 32
    cuts.AddCutCalo("80085123","111113305f032230000","01631031000000d0"); // NL 33
    cuts.AddCutCalo("80085123","111113405f032230000","01631031000000d0"); // NL 34
    cuts.AddCutCalo("80085123","111113805f032230000","01631031000000d0"); // NL 38
    cuts.AddCutCalo("80085123","111113905f032230000","01631031000000d0"); // NL 39
  } else if (trainConfig == 234) {
    cuts.AddCutCalo("80085123","111113106f032230000","01631031000000d0"); // 30/35ns
    cuts.AddCutCalo("80085123","111113104f032230000","01631031000000d0"); // 100ns
    cuts.AddCutCalo("80085123","111113107f032230000","01631031000000d0"); // 30ns
  } else if (trainConfig == 235){ // opening angle variations
    cuts.AddCutCalo("80085113","111113105f032230000","01631031000000b0"); // min opening angle 0.0152
    cuts.AddCutCalo("80085113","111113105f032230000","01631031000000c0"); // min opening angle 0.016
    cuts.AddCutCalo("80085113","111113105f032230000","01631031000000e0"); // min opening angle 0.018
    cuts.AddCutCalo("80085113","111113105f032230000","01631031000000f0"); // min opening angle 0.019
    cuts.AddCutCalo("80085113","111113105f032230000","01633031000000d0"); // rapidity variation  y<0.6
    cuts.AddCutCalo("80085113","111113105f032230000","01634031000000d0"); // rapidity variation  y<0.5
  } else if (trainConfig == 236){
    cuts.AddCutCalo("80085113","111113105f0322l0000","01631031000000d0"); // M02 pt dep with new std: 0.32, 0.0072, 0.5
    cuts.AddCutCalo("80085113","111113105f0322d0000","01631031000000d0"); // M02, pt dep with  0.27, 0.0072, 0.4
    cuts.AddCutCalo("80085113","111113105f0322e0000","01631031000000d0"); // M02, pt dep with  0.31, 0.0072, 0.5
    cuts.AddCutCalo("80085113","111113105f0322f0000","01631031000000d0"); // M02, pt dep with  0.36, 0.0072, 0.7
    cuts.AddCutCalo("80085113","111113105f0322m0000","01631031000000d0"); // M02, pt dep with  0.32, 0.0152, 0.5
  } else if (trainConfig == 237){ // opening angle variations
    cuts.AddCutCalo("80085123","111113105f032230000","01631031000000b0"); // min opening angle 0.0152
    cuts.AddCutCalo("80085123","111113105f032230000","01631031000000c0"); // min opening angle 0.016
    cuts.AddCutCalo("80085123","111113105f032230000","01631031000000e0"); // min opening angle 0.018
    cuts.AddCutCalo("80085123","111113105f032230000","01631031000000f0"); // min opening angle 0.019
    cuts.AddCutCalo("80085123","111113105f032230000","01633031000000d0"); // rapidity variation  y<0.6
    cuts.AddCutCalo("80085123","111113105f032230000","01634031000000d0"); // rapidity variation  y<0.5
  } else if (trainConfig == 238){
    cuts.AddCutCalo("80085123","111113105f0322l0000","01631031000000d0"); // M02 pt dep with new std: 0.32, 0.0072, 0.5
    cuts.AddCutCalo("80085123","111113105f0322d0000","01631031000000d0"); // M02, pt dep with  0.27, 0.0072, 0.4
    cuts.AddCutCalo("80085123","111113105f0322e0000","01631031000000d0"); // M02, pt dep with  0.31, 0.0072, 0.5
    cuts.AddCutCalo("80085123","111113105f0322f0000","01631031000000d0"); // M02, pt dep with  0.36, 0.0072, 0.7
    cuts.AddCutCalo("80085123","111113105f0322m0000","01631031000000d0"); // M02, pt dep with  0.32, 0.0152, 0.5



  } else if (trainConfig == 240){ // EMCAL clusters standard cuts
    cuts.AddCutCalo("80083113","111113105f032230000","01631031000000d0"); // 0-100
  } else if (trainConfig == 241){ // no time cut
    cuts.AddCutCalo("80083113","111113100f032230000","01631031000000d0"); // 0-100
  } else if (trainConfig == 242){ // EMCAL clusters standard cuts, cent vars, non lin variations
    cuts.AddCutCalo("a0183113","111113105f032230000","01631031000000d0"); // 0-5
    cuts.AddCutCalo("a1283113","111113105f032230000","01631031000000d0"); // 5-10
    cuts.AddCutCalo("81283113","111113105f032230000","01631031000000d0"); // 10-20
  } else if (trainConfig == 243){ // EMCAL clusters standard cuts, cent vars, non lin variations
    cuts.AddCutCalo("82483113","111113105f032230000","01631031000000d0"); // 20-40
    cuts.AddCutCalo("84683113","111113105f032230000","01631031000000d0"); // 40-60
    cuts.AddCutCalo("86883113","111113105f032230000","01631031000000d0"); // 60-80
    cuts.AddCutCalo("88083113","111113105f032230000","01631031000000d0"); // 80-100
  } else if (trainConfig == 244){
    cuts.AddCutCalo("80083123","111113105f032230000","01631031000000d0"); // 0-100
    // CONFIGS for sys JJ MC or MB MC
  } else if (trainConfig == 245) { // CALO variations
    cuts.AddCutCalo("80083113","111113105f032230000","01631031000000d0"); // std
    cuts.AddCutCalo("80083113","111113105f022230000","01631031000000d0"); // 0.6 GeV/c
    cuts.AddCutCalo("80083113","111113105f042230000","01631031000000d0"); // 0.8 GeV/c
  } else if (trainConfig == 246) {
    cuts.AddCutCalo("80083113","111113105f031230000","01631031000000d0"); // n cells >= 1
    cuts.AddCutCalo("80083113","111113105f033230000","01631031000000d0"); // n cells >= 3
    cuts.AddCutCalo("80083113","111113105f032200000","01631031000000d0"); // no max M02 cut
    cuts.AddCutCalo("80083113","111113105f032250000","01631031000000d0"); // M02 < 0.3
    cuts.AddCutCalo("80083113","111113105f0322k0000","01631031000000d0"); // M02, pT-dep
  } else if (trainConfig == 247) {
    cuts.AddCutCalo("80083113","111113105e032230000","01631031000000d0"); // TM var e/p 2.0
    cuts.AddCutCalo("80083113","111113105g032230000","01631031000000d0"); // TM var e/p 1.5
    cuts.AddCutCalo("80083113","111113105h032230000","01631031000000d0"); // TM var e/p 1.25
    cuts.AddCutCalo("80083113","1111131057032230000","01631031000000d0"); // TM var no veto
    cuts.AddCutCalo("80083113","111113105f032230000","01631061000000d0"); // alpha meson variation 1   0<alpha<0.8
    cuts.AddCutCalo("80083113","111113105f032230000","01631051000000d0"); // alpha meson variation 2  0<alpha<0.75
  } else if (trainConfig == 248) {
    cuts.AddCutCalo("80083113","111113205f032230000","01631031000000d0"); // NL 32
    cuts.AddCutCalo("80083113","111113305f032230000","01631031000000d0"); // NL 33
    cuts.AddCutCalo("80083113","111113405f032230000","01631031000000d0"); // NL 34
    cuts.AddCutCalo("80083113","111113805f032230000","01631031000000d0"); // NL 38
    cuts.AddCutCalo("80083113","111113905f032230000","01631031000000d0"); // NL 39
  } else if (trainConfig == 249) {
    cuts.AddCutCalo("80083113","111113106f032230000","01631031000000d0"); // 30/35ns
    cuts.AddCutCalo("80083113","111113104f032230000","01631031000000d0"); // 100ns
    cuts.AddCutCalo("80083113","111113107f032230000","01631031000000d0"); // 30ns
  // CONFIGS for sys JJ MC with MB BG + JJ
  } else if (trainConfig == 250) { // CALO variations
    cuts.AddCutCalo("80083123","111113105f032230000","01631031000000d0"); // std
    cuts.AddCutCalo("80083123","111113105f022230000","01631031000000d0"); // 0.6 GeV/c
    cuts.AddCutCalo("80083123","111113105f042230000","01631031000000d0"); // 0.8 GeV/c
  } else if (trainConfig == 251) {
    cuts.AddCutCalo("80083123","111113105f031230000","01631031000000d0"); // n cells >= 1
    cuts.AddCutCalo("80083123","111113105f033230000","01631031000000d0"); // n cells >= 3
    cuts.AddCutCalo("80083123","111113105f032200000","01631031000000d0"); // no max M02 cut
    cuts.AddCutCalo("80083123","111113105f032250000","01631031000000d0"); // M02 < 0.3
    cuts.AddCutCalo("80083123","111113105f0322k0000","01631031000000d0"); // M02, pT-dep
  } else if (trainConfig == 252) {
    cuts.AddCutCalo("80083123","111113105e032230000","01631031000000d0"); // TM var e/p 2.0
    cuts.AddCutCalo("80083123","111113105g032230000","01631031000000d0"); // TM var e/p 1.5
    cuts.AddCutCalo("80083123","111113105h032230000","01631031000000d0"); // TM var e/p 1.25
    cuts.AddCutCalo("80083123","1111131057032230000","01631031000000d0"); // TM var no veto
    cuts.AddCutCalo("80083123","111113105f032230000","01631061000000d0"); // alpha meson variation 1   0<alpha<0.8
    cuts.AddCutCalo("80083123","111113105f032230000","01631051000000d0"); // alpha meson variation 2  0<alpha<0.75
  } else if (trainConfig == 253) {
    cuts.AddCutCalo("80083123","111113205f032230000","01631031000000d0"); // NL 32
    cuts.AddCutCalo("80083123","111113305f032230000","01631031000000d0"); // NL 33
    cuts.AddCutCalo("80083123","111113405f032230000","01631031000000d0"); // NL 34
    cuts.AddCutCalo("80083123","111113805f032230000","01631031000000d0"); // NL 38
    cuts.AddCutCalo("80083123","111113905f032230000","01631031000000d0"); // NL 39
  } else if (trainConfig == 254) {
    cuts.AddCutCalo("80083123","111113106f032230000","01631031000000d0"); // 30/35ns
    cuts.AddCutCalo("80083123","111113104f032230000","01631031000000d0"); // 100ns
    cuts.AddCutCalo("80083123","111113107f032230000","01631031000000d0"); // 30ns
  } else if (trainConfig == 255){ // opening angle variations
    cuts.AddCutCalo("80083113","111113105f032230000","01631031000000b0"); // min opening angle 0.0152
    cuts.AddCutCalo("80083113","111113105f032230000","01631031000000c0"); // min opening angle 0.016
    cuts.AddCutCalo("80083113","111113105f032230000","01631031000000e0"); // min opening angle 0.018
    cuts.AddCutCalo("80083113","111113105f032230000","01631031000000f0"); // min opening angle 0.019
    cuts.AddCutCalo("80083113","111113105f032230000","01633031000000d0"); // rapidity variation  y<0.6
    cuts.AddCutCalo("80083113","111113105f032230000","01634031000000d0"); // rapidity variation  y<0.5
  } else if (trainConfig == 256){
    cuts.AddCutCalo("80083113","111113105f0322l0000","01631031000000d0"); // M02 pt dep with new std: 0.32, 0.0072, 0.5
    cuts.AddCutCalo("80083113","111113105f0322d0000","01631031000000d0"); // M02, pt dep with  0.27, 0.0072, 0.4
    cuts.AddCutCalo("80083113","111113105f0322e0000","01631031000000d0"); // M02, pt dep with  0.31, 0.0072, 0.5
    cuts.AddCutCalo("80083113","111113105f0322f0000","01631031000000d0"); // M02, pt dep with  0.36, 0.0072, 0.7
    cuts.AddCutCalo("80083113","111113105f0322m0000","01631031000000d0"); // M02, pt dep with  0.32, 0.0152, 0.5
  } else if (trainConfig == 257){ // opening angle variations
    cuts.AddCutCalo("80083123","111113105f032230000","01631031000000b0"); // min opening angle 0.0152
    cuts.AddCutCalo("80083123","111113105f032230000","01631031000000c0"); // min opening angle 0.016
    cuts.AddCutCalo("80083123","111113105f032230000","01631031000000e0"); // min opening angle 0.018
    cuts.AddCutCalo("80083123","111113105f032230000","01631031000000f0"); // min opening angle 0.019
    cuts.AddCutCalo("80083123","111113105f032230000","01633031000000d0"); // rapidity variation  y<0.6
    cuts.AddCutCalo("80083123","111113105f032230000","01634031000000d0"); // rapidity variation  y<0.5
  } else if (trainConfig == 258){
    cuts.AddCutCalo("80083123","111113105f0322l0000","01631031000000d0"); // M02 pt dep with new std: 0.32, 0.0072, 0.5
    cuts.AddCutCalo("80083123","111113105f0322d0000","01631031000000d0"); // M02, pt dep with  0.27, 0.0072, 0.4
    cuts.AddCutCalo("80083123","111113105f0322e0000","01631031000000d0"); // M02, pt dep with  0.31, 0.0072, 0.5
    cuts.AddCutCalo("80083123","111113105f0322f0000","01631031000000d0"); // M02, pt dep with  0.36, 0.0072, 0.7
    cuts.AddCutCalo("80083123","111113105f0322m0000","01631031000000d0"); // M02, pt dep with  0.32, 0.0152, 0.5

  } else if (trainConfig == 270) { // CALO variations
    cuts.AddCutCalo("80010113","411793205f03h230000","01631031000000d0"); // NCell cut effi MC (param 1)
    cuts.AddCutCalo("80010113","411793205f03i230000","01631031000000d0"); // NCell cut effi MC (param 2)
    cuts.AddCutCalo("80010113","411793205f03j230000","01631031000000d0"); // NCell cut effi data (param 1)
    cuts.AddCutCalo("80010113","411793205f03k230000","01631031000000d0"); // NCell cut effi data (param 2)
  } else if (trainConfig == 271) { // CALO variations
    cuts.AddCutCalo("80010113","111113205f03h230000","01631031000000d0"); // NCell cut effi MC (param 1)
    cuts.AddCutCalo("80010113","111113205f03i230000","01631031000000d0"); // NCell cut effi MC (param 2)
    cuts.AddCutCalo("80010113","111113205f03j230000","01631031000000d0"); // NCell cut effi data (param 1)
    cuts.AddCutCalo("80010113","111113205f03k230000","01631031000000d0"); // NCell cut effi data (param 2)
  } else if (trainConfig == 272) { // CALO variations
    cuts.AddCutCalo("80010123","111113205f03h230000","01631031000000d0"); // NCell cut effi MC (param 1)
    cuts.AddCutCalo("80010123","111113205f03i230000","01631031000000d0"); // NCell cut effi MC (param 2)
    cuts.AddCutCalo("80010123","111113205f03j230000","01631031000000d0"); // NCell cut effi data (param 1)
    cuts.AddCutCalo("80010123","111113205f03k230000","01631031000000d0"); // NCell cut effi data (param 2)
  } else if (trainConfig == 273) { // CALO variations
    cuts.AddCutCalo("80085113","111113205f03h230000","01631031000000d0"); // NCell cut effi MC (param 1)
    cuts.AddCutCalo("80085113","111113205f03i230000","01631031000000d0"); // NCell cut effi MC (param 2)
    cuts.AddCutCalo("80085113","111113205f03j230000","01631031000000d0"); // NCell cut effi data (param 1)
    cuts.AddCutCalo("80085113","111113205f03k230000","01631031000000d0"); // NCell cut effi data (param 2)
  } else if (trainConfig == 274) { // CALO variations
    cuts.AddCutCalo("80085123","111113205f03h230000","01631031000000d0"); // NCell cut effi MC (param 1)
    cuts.AddCutCalo("80085123","111113205f03i230000","01631031000000d0"); // NCell cut effi MC (param 2)
    cuts.AddCutCalo("80085123","111113205f03j230000","01631031000000d0"); // NCell cut effi data (param 1)
    cuts.AddCutCalo("80085123","111113205f03k230000","01631031000000d0"); // NCell cut effi data (param 2)
  } else if (trainConfig == 275) { // CALO variations
    cuts.AddCutCalo("80083113","111113205f03h230000","01631031000000d0"); // NCell cut effi MC (param 1)
    cuts.AddCutCalo("80083113","111113205f03i230000","01631031000000d0"); // NCell cut effi MC (param 2)
    cuts.AddCutCalo("80083113","111113205f03j230000","01631031000000d0"); // NCell cut effi data (param 1)
    cuts.AddCutCalo("80083113","111113205f03k230000","01631031000000d0"); // NCell cut effi data (param 2)
  } else if (trainConfig == 276) { // CALO variations
    cuts.AddCutCalo("80083123","111113205f03h230000","01631031000000d0"); // NCell cut effi MC (param 1)
    cuts.AddCutCalo("80083123","111113205f03i230000","01631031000000d0"); // NCell cut effi MC (param 2)
    cuts.AddCutCalo("80083123","111113205f03j230000","01631031000000d0"); // NCell cut effi data (param 1)
    cuts.AddCutCalo("80083123","111113205f03k230000","01631031000000d0"); // NCell cut effi data (param 2)
    
  } else if (trainConfig == 277) { // CALO variations
    cuts.AddCutCalo("80010113","411793205f032230000","0r631031000000d0"); // 
    cuts.AddCutCalo("80010113","411793205f032230000","0s631031000000d0"); // 
    cuts.AddCutCalo("80010113","411793205f032230000","0t631031000000d0"); // 
  } else if (trainConfig == 278) { // CALO variations
    cuts.AddCutCalo("80010113","411793205f032230000","0u631031000000d0"); // 
    cuts.AddCutCalo("80010113","411793205f032230000","0v631031000000d0"); // 
    cuts.AddCutCalo("80010113","411793205f032230000","0w631031000000d0"); // 
    cuts.AddCutCalo("80010113","411793205f032230000","0x631031000000d0"); // 
  } else if (trainConfig == 279) { // CALO variations
    cuts.AddCutCalo("80010113","111113205f032230000","0r631031000000d0"); // 
    cuts.AddCutCalo("80010113","111113205f032230000","0s631031000000d0"); // 
    cuts.AddCutCalo("80010113","111113205f032230000","0t631031000000d0"); // 
  } else if (trainConfig == 280) { // CALO variations
    cuts.AddCutCalo("80010113","111113205f032230000","0u631031000000d0"); // 
    cuts.AddCutCalo("80010113","111113205f032230000","0v631031000000d0"); // 
    cuts.AddCutCalo("80010113","111113205f032230000","0w631031000000d0"); // 
    cuts.AddCutCalo("80010113","111113205f032230000","0x631031000000d0"); // 
  } else if (trainConfig == 281) { // CALO variations
    cuts.AddCutCalo("80010123","111113205f032230000","0r631031000000d0"); // 
    cuts.AddCutCalo("80010123","111113205f032230000","0s631031000000d0"); // 
    cuts.AddCutCalo("80010123","111113205f032230000","0t631031000000d0"); // 
  } else if (trainConfig == 282) { // CALO variations
    cuts.AddCutCalo("80010123","111113205f032230000","0u631031000000d0"); // 
    cuts.AddCutCalo("80010123","111113205f032230000","0v631031000000d0"); // 
    cuts.AddCutCalo("80010123","111113205f032230000","0w631031000000d0"); // 
    cuts.AddCutCalo("80010123","111113205f032230000","0x631031000000d0"); // 
  } else if (trainConfig == 283) { // CALO variations
    cuts.AddCutCalo("80085113","111113205f032230000","0r631031000000d0"); // 
    cuts.AddCutCalo("80085113","111113205f032230000","0s631031000000d0"); // 
    cuts.AddCutCalo("80085113","111113205f032230000","0t631031000000d0"); // 
  } else if (trainConfig == 284) { // CALO variations  
    cuts.AddCutCalo("80085113","111113205f032230000","0u631031000000d0"); // 
    cuts.AddCutCalo("80085113","111113205f032230000","0v631031000000d0"); // 
    cuts.AddCutCalo("80085113","111113205f032230000","0w631031000000d0"); // 
    cuts.AddCutCalo("80085113","111113205f032230000","0x631031000000d0"); // 
  } else if (trainConfig == 285) { // CALO variations
    cuts.AddCutCalo("80083113","111113205f032230000","0r631031000000d0"); // 
    cuts.AddCutCalo("80083113","111113205f032230000","0s631031000000d0"); // 
    cuts.AddCutCalo("80083113","111113205f032230000","0t631031000000d0"); // 
  } else if (trainConfig == 286) { // CALO variations  
    cuts.AddCutCalo("80083113","111113205f032230000","0u631031000000d0"); // 
    cuts.AddCutCalo("80083113","111113205f032230000","0v631031000000d0"); // 
    cuts.AddCutCalo("80083113","111113205f032230000","0w631031000000d0"); // 
    cuts.AddCutCalo("80083113","111113205f032230000","0x631031000000d0"); // 
  } else if (trainConfig == 287) { // CALO variations
    cuts.AddCutCalo("80085123","111113205f032230000","0r631031000000d0"); // 
    cuts.AddCutCalo("80085123","111113205f032230000","0s631031000000d0"); // 
    cuts.AddCutCalo("80085123","111113205f032230000","0t631031000000d0"); // 
  } else if (trainConfig == 288) { // CALO variations  
    cuts.AddCutCalo("80085123","111113205f032230000","0u631031000000d0"); // 
    cuts.AddCutCalo("80085123","111113205f032230000","0v631031000000d0"); // 
    cuts.AddCutCalo("80085123","111113205f032230000","0w631031000000d0"); // 
    cuts.AddCutCalo("80085123","111113205f032230000","0x631031000000d0"); // 
  } else if (trainConfig == 289) { // CALO variations  
    cuts.AddCutCalo("80083123","111113205f032230000","0r631031000000d0"); // 
    cuts.AddCutCalo("80083123","111113205f032230000","0s631031000000d0"); // 
    cuts.AddCutCalo("80083123","111113205f032230000","0t631031000000d0"); // 
  } else if (trainConfig == 290) { // CALO variations  
    cuts.AddCutCalo("80083123","111113205f032230000","0u631031000000d0"); // 
    cuts.AddCutCalo("80083123","111113205f032230000","0v631031000000d0"); // 
    cuts.AddCutCalo("80083123","111113205f032230000","0w631031000000d0"); // 
    cuts.AddCutCalo("80083123","111113205f032230000","0x631031000000d0"); // 
    
  } else if (trainConfig == 291) {
    cuts.AddCutCalo("00010113","388553205f032230000","01631031000000d0"); 
    cuts.AddCutCalo("00010113","388553205f032230000","0r631031000000d0"); 
    cuts.AddCutCalo("00010113","388553205f032230000","0s631031000000d0"); 
    cuts.AddCutCalo("00010113","388553205f032230000","0t631031000000d0"); 
    cuts.AddCutCalo("00010113","388553205f032230000","0u631031000000d0"); 
    
  // ===============================================================================================
  // Run 1 data PHOS clusters pPb 5TeV
  // ===============================================================================================
  } else if (trainConfig == 301) {  // min energy = 0.3 GeV/c
    cuts.AddCutCalo("80010113","2444453041013200000","0163103100000010"); //standart cut, kINT7 // PHOS clusters
    cuts.AddCutCalo("80062113","2444453041013200000","0163103100000010"); //standard cut, kPHI7  // PHOS clusters
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
    cuts.AddCutCalo("80010113","2444453041033200000","0163103100000010"); // min energy 0.1 GeV
    cuts.AddCutCalo("80010113","2444453041093200000","0163103100000010"); // min energy 0.1 GeV
    cuts.AddCutCalo("80010113","2444453041073200000","0163103100000010"); // min energy 0.2 GeV
    cuts.AddCutCalo("80010113","2444453041083200000","0163103100000010"); // min energy 0.4 GeV
    cuts.AddCutCalo("80010113","2444453041023200000","0163103100000010"); // min energy 0.5 GeV
  } else if(trainConfig == 311){ // second set of variations CLUSTER
    cuts.AddCutCalo("80010113","2444453041013270000","0163103100000010"); // min/max M02  0.1<M<1
    cuts.AddCutCalo("80010113","2444453041013280000","0163103100000010"); // min/max M02  0.1<M<1
    cuts.AddCutCalo("80010113","2444453041012200000","0163103100000010"); // min number 2 cells
    cuts.AddCutCalo("80010113","2444453041014200000","0163103100000010"); // min number 4 cells
  } else if(trainConfig == 312){ // MESON
    cuts.AddCutCalo("80010113","2444453041013200000","0163403100000010"); // rapidity variation  y<0.5
    cuts.AddCutCalo("80010113","2444453041013200000","0163803100000010"); // rapidity variation  y<0.25
    cuts.AddCutCalo("80010113","2444453041013200000","0163106100000010"); // alpha meson variation 1 0<alpha<0.8
    cuts.AddCutCalo("80010113","2444453041013200000","0163105100000010"); // alpha meson variation 2 0<alpha<0.75
  } else if(trainConfig == 313){ // fourth set of variations
    cuts.AddCutCalo("80010113","2444453044013200000","0163103100000010"); // tm variation
    cuts.AddCutCalo("80010113","2444453040013200000","0163103100000010"); // tm variation
    cuts.AddCutCalo("80010113","2444453045013200000","0163103100000010"); // tm variation
    cuts.AddCutCalo("80010113","2444453046013200000","0163103100000010"); // tm variation
    cuts.AddCutCalo("80010113","2444453041013200000","0163103100000000"); // min opening angle 0    -> open
    cuts.AddCutCalo("80010113","2444453041013200000","0163103100000030"); // min opening angle 0.01 -> 2 cell diag
  } else if(trainConfig == 314){ 
    cuts.AddCutCalo("80010113","2444454041033200000","0163103100000010"); // 
    cuts.AddCutCalo("80010113","2444451041033200000","0163103100000010"); // 
    cuts.AddCutCalo("80010113","2444452041033200000","0163103100000010"); // 
  } else if(trainConfig == 315){ 
    cuts.AddCutCalo("80010113","2444453141033200000","0163103100000010"); // 
    cuts.AddCutCalo("80010113","2444453241033200000","0163103100000010"); // 
    cuts.AddCutCalo("80010113","2444453341033200000","0163103100000010"); // 
  } else if( trainConfig == 316){ 
    cuts.AddCutCalo("80010113","2444453041033200010","0163103100000010"); // 
    cuts.AddCutCalo("80010113","2444453041033200020","0163103100000010"); // 
    cuts.AddCutCalo("80010113","2444453041033200000","0263103100000010"); // 
  } else if( trainConfig == 317){ 
    cuts.AddCutCalo("80010113","2444453041033200300","0163103100000010"); // 
    cuts.AddCutCalo("80010113","2444453041033200200","0163103100000010"); // 
    cuts.AddCutCalo("80010113","2444453041033200400","0163103100000010"); // 
    
    
  } else if(trainConfig == 320){ // first set of variations CLUSTER
    cuts.AddCutCalo("80010123","2444453041033200000","0163103100000010"); // min energy 0.1 GeV
    cuts.AddCutCalo("80010123","2444453041093200000","0163103100000010"); // min energy 0.1 GeV
    cuts.AddCutCalo("80010123","2444453041073200000","0163103100000010"); // min energy 0.2 GeV
    cuts.AddCutCalo("80010123","2444453041083200000","0163103100000010"); // min energy 0.4 GeV
    cuts.AddCutCalo("80010123","2444453041023200000","0163103100000010"); // min energy 0.5 GeV
  } else if(trainConfig == 321){ // second set of variations CLUSTER
    cuts.AddCutCalo("80010123","2444453041013270000","0163103100000010"); // min/max M02  0.1<M<1
    cuts.AddCutCalo("80010123","2444453041013280000","0163103100000010"); // min/max M02  0.1<M<1
    cuts.AddCutCalo("80010123","2444453041012200000","0163103100000010"); // min number 2 cells
    cuts.AddCutCalo("80010123","2444453041014200000","0163103100000010"); // min number 4 cells
  } else if(trainConfig == 322){ // MESON
    cuts.AddCutCalo("80010123","2444453041013200000","0163403100000010"); // rapidity variation  y<0.5
    cuts.AddCutCalo("80010123","2444453041013200000","0163803100000010"); // rapidity variation  y<0.25
    cuts.AddCutCalo("80010123","2444453041013200000","0163106100000010"); // alpha meson variation 1 0<alpha<0.8
    cuts.AddCutCalo("80010123","2444453041013200000","0163105100000010"); // alpha meson variation 2 0<alpha<0.75
  } else if(trainConfig == 323){ // fourth set of variations
    cuts.AddCutCalo("80010123","2444453044013200000","0163103100000010"); // tm variation
    cuts.AddCutCalo("80010123","2444453040013200000","0163103100000010"); // tm variation
    cuts.AddCutCalo("80010123","2444453045013200000","0163103100000010"); // tm variation
    cuts.AddCutCalo("80010123","2444453046013200000","0163103100000010"); // tm variation
    cuts.AddCutCalo("80010123","2444453041013200000","0163103100000000"); // min opening angle 0    -> open
    cuts.AddCutCalo("80010123","2444453041013200000","0163103100000030"); // min opening angle 0.01 -> 2 cell diag
  } else if(trainConfig == 324){ 
    cuts.AddCutCalo("80010123","2444454041033200000","0163103100000010"); // 
    cuts.AddCutCalo("80010123","2444451041033200000","0163103100000010"); // 
    cuts.AddCutCalo("80010123","2444452041033200000","0163103100000010"); 
  } else if(trainConfig == 325){ 
    cuts.AddCutCalo("80010123","2444453141033200000","0163103100000010"); 
    cuts.AddCutCalo("80010123","2444453241033200000","0163103100000010"); 
    cuts.AddCutCalo("80010123","2444453341033200000","0163103100000010"); 
  } else if(trainConfig == 326){ 
    cuts.AddCutCalo("80010123","2444453041033200010","0163103100000010"); 
    cuts.AddCutCalo("80010123","2444453041033200020","0163103100000010"); 
    cuts.AddCutCalo("80010123","2444453041033200000","0263103100000010"); 
  } else if( trainConfig == 327){
    cuts.AddCutCalo("80010123","2444453041033200300","0163103100000010"); 
    cuts.AddCutCalo("80010123","2444453041033200200","0163103100000010"); 
    cuts.AddCutCalo("80010123","2444453041033200400","0163103100000010"); 
    
  } else if(trainConfig == 330){ // first set of variations CLUSTER
    cuts.AddCutCalo("80062113","2444453041033200000","0163103100000010"); // min energy 0.1 GeV
    cuts.AddCutCalo("80062113","2444453041093200000","0163103100000010"); // min energy 0.1 GeV
    cuts.AddCutCalo("80062113","2444453041073200000","0163103100000010"); // min energy 0.2 GeV
    cuts.AddCutCalo("80062113","2444453041083200000","0163103100000010"); // min energy 0.4 GeV
    cuts.AddCutCalo("80062113","2444453041023200000","0163103100000010"); // min energy 0.5 GeV
  } else if(trainConfig == 331){ // second set of variations CLUSTER
    cuts.AddCutCalo("80062113","2444453041013270000","0163103100000010"); // min/max M02  0.1<M<1
    cuts.AddCutCalo("80062113","2444453041013280000","0163103100000010"); // min/max M02  0.1<M<1
    cuts.AddCutCalo("80062113","2444453041012200000","0163103100000010"); // min number 2 cells
    cuts.AddCutCalo("80062113","2444453041014200000","0163103100000010"); // min number 4 cells
  } else if(trainConfig == 332){ // MESON
    cuts.AddCutCalo("80062113","2444453041013200000","0163403100000010"); // rapidity variation  y<0.5
    cuts.AddCutCalo("80062113","2444453041013200000","0163803100000010"); // rapidity variation  y<0.25
    cuts.AddCutCalo("80062113","2444453041013200000","0163106100000010"); // alpha meson variation 1 0<alpha<0.8
    cuts.AddCutCalo("80062113","2444453041013200000","0163105100000010"); // alpha meson variation 2 0<alpha<0.75
  } else if(trainConfig == 333){ // fourth set of variations
    cuts.AddCutCalo("80062113","2444453044013200000","0163103100000010"); // tm variation
    cuts.AddCutCalo("80062113","2444453040013200000","0163103100000010"); // tm variation
    cuts.AddCutCalo("80062113","2444453045013200000","0163103100000010"); // tm variation
    cuts.AddCutCalo("80062113","2444453046013200000","0163103100000010"); // tm variation
    cuts.AddCutCalo("80062113","2444453041013200000","0163103100000000"); // min opening angle 0    -> open
    cuts.AddCutCalo("80062113","2444453041013200000","0163103100000030"); // min opening angle 0.01 -> 2 cell diag
  } else if(trainConfig == 334){ // fourth set of variations
    cuts.AddCutCalo("80062113","2444454041033200000","0163103100000010"); // 
    cuts.AddCutCalo("80062113","2444451041033200000","0163103100000010"); // 
    cuts.AddCutCalo("80062113","2444452041033200000","0163103100000010"); // 
  } else if(trainConfig == 335){ // fourth set of variations
    cuts.AddCutCalo("80062113","2444453141033200000","0163103100000010"); // 
    cuts.AddCutCalo("80062113","2444453241033200000","0163103100000010"); // 
    cuts.AddCutCalo("80062113","2444453341033200000","0163103100000010"); // 
  } else if(trainConfig == 336){ // fourth set of variations
    cuts.AddCutCalo("80062113","2444453041033200010","0163103100000010"); // 
    cuts.AddCutCalo("80062113","2444453041033200020","0163103100000010"); // 
    cuts.AddCutCalo("80062113","2444453041033200000","0263103100000010"); // 
  } else if( trainConfig == 337){ // fourth set of variations
    cuts.AddCutCalo("80062113","2444453041033200300","0163103100000010"); // 
    cuts.AddCutCalo("80062113","2444453041033200200","0163103100000010"); // 
    cuts.AddCutCalo("80062113","2444453041033200400","0163103100000010"); // 
    
  } else if(trainConfig == 340){ // first set of variations CLUSTER
    cuts.AddCutCalo("80062123","2444453041033200000","0163103100000010"); // min energy 0.1 GeV
    cuts.AddCutCalo("80062123","2444453041093200000","0163103100000010"); // min energy 0.1 GeV
    cuts.AddCutCalo("80062123","2444453041073200000","0163103100000010"); // min energy 0.2 GeV
    cuts.AddCutCalo("80062123","2444453041083200000","0163103100000010"); // min energy 0.4 GeV
    cuts.AddCutCalo("80062123","2444453041023200000","0163103100000010"); // min energy 0.5 GeV
  } else if(trainConfig == 341){ // second set of variations CLUSTER
    cuts.AddCutCalo("80062123","2444453041013270000","0163103100000010"); // min/max M02  0.1<M<1
    cuts.AddCutCalo("80062123","2444453041013280000","0163103100000010"); // min/max M02  0.1<M<1
    cuts.AddCutCalo("80062123","2444453041012200000","0163103100000010"); // min number 2 cells
    cuts.AddCutCalo("80062123","2444453041014200000","0163103100000010"); // min number 4 cells
  } else if(trainConfig == 342){ // MESON
    cuts.AddCutCalo("80062123","2444453041013200000","0163403100000010"); // rapidity variation  y<0.5
    cuts.AddCutCalo("80062123","2444453041013200000","0163803100000010"); // rapidity variation  y<0.25
    cuts.AddCutCalo("80062123","2444453041013200000","0163106100000010"); // alpha meson variation 1 0<alpha<0.8
    cuts.AddCutCalo("80062123","2444453041013200000","0163105100000010"); // alpha meson variation 2 0<alpha<0.75
  } else if(trainConfig == 343){ // fourth set of variations
    cuts.AddCutCalo("80062123","2444453044013200000","0163103100000010"); // tm variation
    cuts.AddCutCalo("80062123","2444453040013200000","0163103100000010"); // tm variation
    cuts.AddCutCalo("80062123","2444453045013200000","0163103100000010"); // tm variation
    cuts.AddCutCalo("80062123","2444453046013200000","0163103100000010"); // tm variation
    cuts.AddCutCalo("80062123","2444453041013200000","0163103100000000"); // min opening angle 0    -> open
    cuts.AddCutCalo("80062123","2444453041013200000","0163103100000030"); // min opening angle 0.01 -> 2 cell diag
  } else if(trainConfig == 344){ 
    cuts.AddCutCalo("80062123","2444454041033200000","0163103100000010"); // 
    cuts.AddCutCalo("80062123","2444451041033200000","0163103100000010"); // 
    cuts.AddCutCalo("80062123","2444452041033200000","0163103100000010"); 
  } else if(trainConfig == 345){ 
    cuts.AddCutCalo("80062123","2444453141033200000","0163103100000010"); 
    cuts.AddCutCalo("80062123","2444453241033200000","0163103100000010"); 
    cuts.AddCutCalo("80062123","2444453341033200000","0163103100000010"); 
  } else if(trainConfig == 346){ 
    cuts.AddCutCalo("80062123","2444453041033200010","0163103100000010"); 
    cuts.AddCutCalo("80062123","2444453041033200020","0163103100000010"); 
    cuts.AddCutCalo("80062123","2444453041033200000","0263103100000010"); 
  } else if( trainConfig == 347){ 
    cuts.AddCutCalo("80062123","2444453041033200300","0163103100000010"); 
    cuts.AddCutCalo("80062123","2444453041033200200","0163103100000010"); 
    cuts.AddCutCalo("80062123","2444453041033200400","0163103100000010"); 

  // JJ MC cut strings
  } else if (trainConfig == 360) {  // min energy = 0.3 GeV/c
    cuts.AddCutCalo("80010123","2444400041013200000","0163103100000010"); //standart cut, kINT7 // PHOS clusters
    cuts.AddCutCalo("80062123","2444400041013200000","0163103100000010"); //standard cut, kPHI7  // PHOS clusters
  } else if (trainConfig == 361) {  // PHOS default cent dep
    cuts.AddCutCalo("80210123","2444401041013200000","0163103100000010"); // 0-20
    cuts.AddCutCalo("82410123","2444401041013200000","0163103100000010"); // 20-40
    cuts.AddCutCalo("84610123","2444401041013200000","0163103100000010"); // 40-60
    cuts.AddCutCalo("86010123","2444401041013200000","0163103100000010"); // 60-100

  // ===============================================================================================
  // Run 2 data PHOS clusters pPb 5TeV
  // ===============================================================================================
  // INT7 triggers
  } else if (trainConfig == 500) {  // PHOS  INT7
    cuts.AddCutCalo("80010113","24466410ha012200000","0163103100000010"); // standard without non-lin
    cuts.AddCutCalo("80010113","244664105a012200000","0163103100000010"); // standard without non-lin
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
    cuts.AddCutCalo("80010113","24466410ha012200000","0163103100000010"); // non lin 0-100%
    cuts.AddCutCalo("80110113","24466410ha012200000","0163103100000010"); // non lin 0-10%
    cuts.AddCutCalo("81210113","24466410ha012200000","0163103100000010"); // non lin 10-20%
    cuts.AddCutCalo("82410113","24466410ha012200000","0163103100000010"); // non lin 20-40%
    cuts.AddCutCalo("84610113","24466410ha012200000","0163103100000010"); // non lin 40-60%
    cuts.AddCutCalo("86810113","24466410ha012200000","0163103100000010"); // non lin 60-80%
    cuts.AddCutCalo("88010113","24466410ha012200000","0163103100000010"); // non lin 80-100%
  } else if (trainConfig == 504) {  // PHOS  INT7 with cents
    cuts.AddCutCalo("80010113","24466410ha012200000","0163103100000010"); // non lin 0-100%
    cuts.AddCutCalo("80210113","24466410ha012200000","0163103100000010"); // non lin 0-20%
    cuts.AddCutCalo("86010113","24466410ha012200000","0163103100000010"); // non lin 60-100%
    cuts.AddCutCalo("a0110113","24466410ha012200000","0163103100000010"); // non lin 0-5%
    cuts.AddCutCalo("a1210113","24466410ha012200000","0163103100000010"); // non lin 5-10%
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
  } else if (trainConfig == 515) {  // PHOS  INT7
    cuts.AddCutCalo("80010113","24466000ha012200000","0163103100000010"); // standard without non-lin
    cuts.AddCutCalo("80010113","244660005a012200000","0163103100000010"); // standard without non-lin
  } else if (trainConfig == 516) {  // PHOS  INT7
    cuts.AddCutCalo("80010113","24466530ha012200000","0163103100000010"); // non-lin: shifting to pi0 mass using PHOS-PHOS
    cuts.AddCutCalo("80010113","24466540ha012200000","0163103100000010"); // non-lin: shifting to pi0 mass using PCM-PHOS
    cuts.AddCutCalo("80010113","24466000ha012200000","0163103100000010"); // without non-lin
  } else if (trainConfig == 517) { // new PHOS default cut for RUN2: Use M02 and Ncell cut only for clusters with a minimum energy of 1 GeV
    cuts.AddCutCalo("80010113","24466530ha01c200000","0163103100000010"); // use Ncell cut only for clusters with a minimum energy of 1 GeV
    cuts.AddCutCalo("80010113","24466530ha012c00000","0163103100000010"); // use M02 cut only for clusters with a minimum energy of 1 GeV
    cuts.AddCutCalo("80010113","24466530ha01cc00000","0163103100000010"); // use Ncell & M02 cut only for clusters with a minimum energy of 1 GeV
  } else if (trainConfig == 518) { // new PHOS default cut for RUN2, different non-lins
    cuts.AddCutCalo("80010113","24466000ha01cc00000","0163103100000010"); // without non-lin
    cuts.AddCutCalo("80010113","24466530ha01cc00000","0163103100000010"); // non-lin 1
    cuts.AddCutCalo("80010113","24466540ha01cc00000","0163103100000010"); // non-lin 2

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
  } else if (trainConfig == 524) {  // JJ MC different non lins (shiftig to pi0 mass)
    cuts.AddCutCalo("80010123","24466000ha012200000","0163103100000010"); // 0-100% without non-lin
    cuts.AddCutCalo("80010123","24466530ha012200000","0163103100000010"); // 0-100% non-lin 1
    cuts.AddCutCalo("80010123","24466540ha012200000","0163103100000010"); // 0-100% non-lin 2
  } else if (trainConfig == 525) {  // JJ MC new PHOS default cuts with different non lins (shiftig to pi0 mass)
    cuts.AddCutCalo("80010123","24466000ha01cc00000","0163103100000010"); // 0-100% without non-lin
    cuts.AddCutCalo("80010123","24466530ha01cc00000","0163103100000010"); // 0-100% non-lin 1
    cuts.AddCutCalo("80010123","24466540ha01cc00000","0163103100000010"); // 0-100% non-lin 2

  // Variations for systematics
  } else if ( trainConfig == 530) { // NL variations (standard: 42 PHOS ML)
    cuts.AddCutCalo("80010113","24466530ha01cc00000","0163103100000010");
    cuts.AddCutCalo("80010113","24466540ha01cc00000","0163103100000010");
    cuts.AddCutCalo("80010113","24466010ha01cc00000","0163103100000010"); // PHOS group NL
    cuts.AddCutCalo("80010113","24466000ha01cc00000","0163103100000010"); // PHOS group NL
  } else if ( trainConfig == 531) { // distance to bad channel variations (standard: no cut)
    cuts.AddCutCalo("80010113","24466531ha01cc00000","0163103100000010"); // dist. to bad channel = 1
    cuts.AddCutCalo("80010113","24466532ha01cc00000","0163103100000010"); // dist. to bad channel = 2
    cuts.AddCutCalo("80010113","24466533ha01cc00000","0163103100000010"); // dist. to bad channel = 3
  } else if ( trainConfig == 532) {  // timing variations - bunch spacing: 100ns (standard: 50ns)
                                    // min. number of cells per cluster variations (standard: 2)
    cuts.AddCutCalo("80010113","24466530ia01cc00000","0163103100000010"); // 30ns
    cuts.AddCutCalo("80010113","24466530ga01cc00000","0163103100000010"); // -20ns/25ns
    cuts.AddCutCalo("80010113","24466530ha01dc00000","0163103100000010"); // nCells > 3 above 1 GeV
  } else if ( trainConfig == 533) { // track matching variations
    cuts.AddCutCalo("80010113","24466530h001cc00000","0163103100000010"); // without track matching
    cuts.AddCutCalo("80010113","24466530h501cc00000","0163103100000010"); // tm variation
    cuts.AddCutCalo("80010113","24466530h601cc00000","0163103100000010"); // tm variation
    cuts.AddCutCalo("80010113","24466530ha01cc00000","0263303100000010"); // BG
  } else if ( trainConfig == 534) { // min. cluster energy variations (standard: 0.5 MeV)
    cuts.AddCutCalo("80010113","24466530ha07cc00000","0163103100000010"); // 0.2 MeV
    cuts.AddCutCalo("80010113","24466530ha09cc00000","0163103100000010"); // 0.1 MeV
    cuts.AddCutCalo("80010113","24466530ha08cc00000","0163103100000010"); // 0.4 MeV
    cuts.AddCutCalo("80010113","24466530ha02cc00000","0163103100000010"); // 0.5 MeV
  } else if ( trainConfig == 535) {
    cuts.AddCutCalo("80010113","24466530ha01cb00000","0163103100000010"); // M02 min 0.002, E > 1 GeV
    cuts.AddCutCalo("80010113","24466530ha01cd00000","0163103100000010"); // M02 min 0.2, E > 1 GeV
    cuts.AddCutCalo("80010113","24466530ha01cc00010","0163103100000010"); // dispersion < 2
    cuts.AddCutCalo("80010113","24466530ha01cc00030","0163103100000010"); // dispersion < 2
  } else if ( trainConfig == 536){ // MESON
    cuts.AddCutCalo("80010113","24466530ha01cc00000","0163303100000010"); // rapidity variation  y<0.5
    cuts.AddCutCalo("80010113","24466530ha01cc00000","0163803100000010"); // rapidity variation  y<0.5
    cuts.AddCutCalo("80010113","24466530ha01cc00000","0163106100000010"); // alpha meson variation 1 0<alpha<0.8
    cuts.AddCutCalo("80010113","24466530ha01cc00000","0163105100000010"); // alpha meson variation 2 0<alpha<0.75
  } else if( trainConfig == 537){ // fourth set of variations
    cuts.AddCutCalo("80010113","24466530ha01cc00300","0163103100000000"); // conv rec 0.03
    cuts.AddCutCalo("80010113","24466530ha01cc00200","0163103100000030"); // conv rec 0.025
    cuts.AddCutCalo("80010113","24466530ha01cc00400","0163103100000010"); // conv rec 0.035
    
// Variations for systematics - JJ MC
  } else if (trainConfig == 540) { // NL variations (standard: 42 PHOS ML)
    cuts.AddCutCalo("80010123","24466530ha01cc00000","0163103100000010");
    cuts.AddCutCalo("80010123","24466540ha01cc00000","0163103100000010");
    cuts.AddCutCalo("80010123","24466010ha01cc00000","0163103100000010"); // PHOS group NL
    cuts.AddCutCalo("80010123","24466000ha01cc00000","0163103100000010"); // PHOS group NL
  } else if (trainConfig == 541) { // distance to bad channel variations (standard: no cut)
    cuts.AddCutCalo("80010123","24466531ha01cc00000","0163103100000010"); // dist. to bad channel = 1
    cuts.AddCutCalo("80010123","24466532ha01cc00000","0163103100000010"); // dist. to bad channel = 2
    cuts.AddCutCalo("80010123","24466533ha01cc00000","0163103100000010"); // dist. to bad channel = 3
  } else if (trainConfig == 542) {  // timing variations - bunch spacing: 100ns (standard: 50ns)
                                    // min. number of cells per cluster variations (standard: 2)
    cuts.AddCutCalo("80010123","24466530ia01cc00000","0163103100000010"); // 30ns
    cuts.AddCutCalo("80010123","24466530ga01cc00000","0163103100000010"); // -20ns/25ns
    cuts.AddCutCalo("80010123","24466530ha01dc00000","0163103100000010"); // nCells > 3 above 1 GeV
  } else if (trainConfig == 543) { // track matching variations
    cuts.AddCutCalo("80010123","24466530h001cc00000","0163103100000010"); // without track matching
    cuts.AddCutCalo("80010123","24466530h501cc00000","0163103100000010"); // tm variation
    cuts.AddCutCalo("80010123","24466530h601cc00000","0163103100000010"); // tm variation
    cuts.AddCutCalo("80010123","24466530ha01cc00000","0263303100000010"); // BG
  } else if (trainConfig == 544) { // min. cluster energy variations (standard: 0.5 MeV)
    cuts.AddCutCalo("80010123","24466530ha07cc00000","0163103100000010"); // 0.2 MeV
    cuts.AddCutCalo("80010123","24466530ha09cc00000","0163103100000010"); // 0.1 MeV
    cuts.AddCutCalo("80010123","24466530ha08cc00000","0163103100000010"); // 0.4 MeV
    cuts.AddCutCalo("80010123","24466530ha02cc00000","0163103100000010"); // 0.5 MeV
  } else if (trainConfig == 545) {
    cuts.AddCutCalo("80010123","24466530ha01cb00000","0163103100000010"); // M02 min 0.002, E > 1 GeV
    cuts.AddCutCalo("80010123","24466530ha01cd00000","0163103100000010"); // M02 min 0.2, E > 1 GeV
    cuts.AddCutCalo("80010123","24466530ha01cc00010","0163103100000010"); // dispersion < 2
    cuts.AddCutCalo("80010123","24466530ha01cc00030","0163103100000010"); // dispersion < 2
  } else if ( trainConfig == 536){ // MESON
    cuts.AddCutCalo("80010123","24466530ha01cc00000","0163303100000010"); // rapidity variation  y<0.5
    cuts.AddCutCalo("80010123","24466530ha01cc00000","0163803100000010"); // rapidity variation  y<0.5
    cuts.AddCutCalo("80010123","24466530ha01cc00000","0163106100000010"); // alpha meson variation 1 0<alpha<0.8
    cuts.AddCutCalo("80010123","24466530ha01cc00000","0163105100000010"); // alpha meson variation 2 0<alpha<0.75
  } else if( trainConfig == 567){ // fourth set of variations
    cuts.AddCutCalo("80010123","24466530ha01cc00300","0163103100000000"); // conv rec 0.03
    cuts.AddCutCalo("80010123","24466530ha01cc00200","0163103100000030"); // conv rec 0.025
    cuts.AddCutCalo("80010123","24466530ha01cc00400","0163103100000010"); // conv rec 0.035
    

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

  // PHOS clusters standard cuts, triggers, no nonlin with TM
  } else if (trainConfig == 700){
    cuts.AddCutCalo("80010123","24466000ha012200000","0163103100000010"); // INT7
  } else if (trainConfig == 701){
    cuts.AddCutCalo("80062123","24466000ha012200000","0163103100000010"); // PHI7
  } else if (trainConfig == 702){
    cuts.AddCutCalo("80010123","24466590ha012200000","0163103100000010"); // INT7
  } else if (trainConfig == 703){
    cuts.AddCutCalo("80062123","24466590ha012200000","0163103100000010"); // PHI7
  } else if (trainConfig == 704){
    cuts.AddCutCalo("80010123","24466690ha012200000","0163103100000010"); // INT7
  } else if (trainConfig == 705){
    cuts.AddCutCalo("80062123","24466690ha012200000","0163103100000010"); // PHI7

  // PHOS clusters standard cuts, triggers, no nonlin, no TM
  } else if (trainConfig == 710){
    cuts.AddCutCalo("80010103","24466000h0012200000","0163103100000010"); // INT7
  } else if (trainConfig == 711){
    cuts.AddCutCalo("80062103","24466000h0012200000","0163103100000010"); // PHI7
  } else if (trainConfig == 712){
    cuts.AddCutCalo("80010103","24466590h0012200000","0163103100000010"); // INT7
  } else if (trainConfig == 713){
    cuts.AddCutCalo("80062103","24466590h0012200000","0163103100000010"); // PHI7
  } else if (trainConfig == 714){
    cuts.AddCutCalo("80010103","24466690h0012200000","0163103100000010"); // INT7
  } else if (trainConfig == 715){
    cuts.AddCutCalo("80062103","24466690h0012200000","0163103100000010"); // PHI7

  // ===============================================================================================
  // Run 1 data EMC triggers only
  // ===============================================================================================

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
    cuts.AddCutCalo("80010103","411793105f032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("8008e103","411793105f032230000","01631031000000d0"); // EG2
    cuts.AddCutCalo("8008d103","411793105f032230000","01631031000000d0"); // EG1
  } else if (trainConfig == 2021){ // EMCAL+DCAL clusters standard cuts, triggers, NL vars
    cuts.AddCutCalo("80010103","411793205f032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("8008e103","411793205f032230000","01631031000000d0"); // EG2
    cuts.AddCutCalo("8008d103","411793205f032230000","01631031000000d0"); // EG1
  } else if (trainConfig == 2022){ // EMCAL+DCAL clusters standard cuts, triggers, NL vars
    cuts.AddCutCalo("80010103","411793305f032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("8008e103","411793305f032230000","01631031000000d0"); // EG2
    cuts.AddCutCalo("8008d103","411793305f032230000","01631031000000d0"); // EG1
  } else if (trainConfig == 2023){ // EMCAL+DCAL clusters standard cuts, triggers, NL vars
    cuts.AddCutCalo("80010103","411793405f032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("8008e103","411793405f032230000","01631031000000d0"); // EG2
    cuts.AddCutCalo("8008d103","411793405f032230000","01631031000000d0"); // EG1
  } else if (trainConfig == 2024){ // open timing, TB NL
    cuts.AddCutCalo("80010103","411796500f032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("80057103","411796500f032230000","01631031000000d0"); // L0
  } else if (trainConfig == 2025){ // open timing, TB NL
    cuts.AddCutCalo("8008e103","411796500f032230000","01631031000000d0"); // L1 low
    cuts.AddCutCalo("8008d103","411796500f032230000","01631031000000d0"); // L1 high
  } else if (trainConfig == 2026){ // open timing, no NL
    cuts.AddCutCalo("80010103","411790000f032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("80057103","411790000f032230000","01631031000000d0"); // L0
  } else if (trainConfig == 2027){ // open timing, TB NL
    cuts.AddCutCalo("8008e103","411790000f032230000","01631031000000d0"); // L1 low
    cuts.AddCutCalo("8008d103","411790000f032230000","01631031000000d0"); // L1 high
  } else if (trainConfig == 2028){ // EMCAL+DCAL clusters standard cuts, triggers, NL via lead cell
    cuts.AddCutCalo("80010103","411793805f032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("8008e103","411793805f032230000","01631031000000d0"); // EG2
    cuts.AddCutCalo("8008d103","411793805f032230000","01631031000000d0"); // EG1
  } else if (trainConfig == 2029){ // EMCAL+DCAL clusters standard cuts, triggers, NL via indiv cells
    cuts.AddCutCalo("80010103","411793905f032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("8008e103","411793905f032230000","01631031000000d0"); // EG2
    cuts.AddCutCalo("8008d103","411793905f032230000","01631031000000d0"); // EG1
  // no TM
  } else if (trainConfig == 2030){ // EMCAL+DCAL clusters standard cuts, triggers, NL vars
    cuts.AddCutCalo("80010103","4117931050032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("8008e103","4117931050032230000","01631031000000d0"); // EG2
    cuts.AddCutCalo("8008d103","4117931050032230000","01631031000000d0"); // EG1
  } else if (trainConfig == 2031){ // EMCAL+DCAL clusters standard cuts, triggers, NL vars
    cuts.AddCutCalo("80010103","4117932050032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("8008e103","4117932050032230000","01631031000000d0"); // EG2
    cuts.AddCutCalo("8008d103","4117932050032230000","01631031000000d0"); // EG1
  } else if (trainConfig == 2032){ // EMCAL+DCAL clusters standard cuts, triggers, NL vars
    cuts.AddCutCalo("80010103","4117933050032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("8008e103","4117933050032230000","01631031000000d0"); // EG2
    cuts.AddCutCalo("8008d103","4117933050032230000","01631031000000d0"); // EG1
  } else if (trainConfig == 2033){ // EMCAL+DCAL clusters standard cuts, triggers, NL vars
    cuts.AddCutCalo("80010103","4117938050032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("8008e103","4117938050032230000","01631031000000d0"); // EG2
    cuts.AddCutCalo("8008d103","4117938050032230000","01631031000000d0"); // EG1
  } else if (trainConfig == 2034){ // EMCAL+DCAL clusters standard cuts, triggers, NL vars
    cuts.AddCutCalo("80010103","4117965050032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("8008e103","4117965050032230000","01631031000000d0"); // EG2
    cuts.AddCutCalo("8008d103","4117965050032230000","01631031000000d0"); // EG1
  } else if (trainConfig == 2035){ // EMCAL+DCAL clusters standard cuts, triggers, NL vars
    cuts.AddCutCalo("80010103","4117900050032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("8008e103","4117900050032230000","01631031000000d0"); // EG2
    cuts.AddCutCalo("8008d103","4117900050032230000","01631031000000d0"); // EG1
  } else if (trainConfig == 2038){ // EMCAL+DCAL clusters standard cuts, triggers, NL via lead cell
    cuts.AddCutCalo("80010103","4117938050032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("8008e103","4117938050032230000","01631031000000d0"); // EG2
    cuts.AddCutCalo("8008d103","4117938050032230000","01631031000000d0"); // EG1
  } else if (trainConfig == 2039){ // EMCAL+DCAL clusters standard cuts, triggers, NL via indiv cells
    cuts.AddCutCalo("80010103","4117939050032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("8008e103","4117939050032230000","01631031000000d0"); // EG2
    cuts.AddCutCalo("8008d103","4117939050032230000","01631031000000d0"); // EG1
  } else if (trainConfig == 2040){ // copy of 2020 for clusterizer cell timing studies
    cuts.AddCutCalo("80010103","411793105f032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("8008e103","411793105f032230000","01631031000000d0"); // EG2
    cuts.AddCutCalo("8008d103","411793105f032230000","01631031000000d0"); // EG1
  } else if (trainConfig == 2041){ // copy of 2020 for clusterizer cell timing studies
    cuts.AddCutCalo("80010103","411793105f032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("8008e103","411793105f032230000","01631031000000d0"); // EG2
    cuts.AddCutCalo("8008d103","411793105f032230000","01631031000000d0"); // EG1
  } else if (trainConfig == 2042){ // copy of 2020 for clusterizer cell timing studies
    cuts.AddCutCalo("80010103","411793105f032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("8008e103","411793105f032230000","01631031000000d0"); // EG2
    cuts.AddCutCalo("8008d103","411793105f032230000","01631031000000d0"); // EG1
  } else if (trainConfig == 2043){ // copy of 2020 for clusterizer cell timing studies
    cuts.AddCutCalo("80010103","411793105f032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("8008e103","411793105f032230000","01631031000000d0"); // EG2
    cuts.AddCutCalo("8008d103","411793105f032230000","01631031000000d0"); // EG1

  } else if (trainConfig == 2044){ // T0 based triggers
    cuts.AddCutCalo("80011103","411793105f032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("8008g103","411793105f032230000","01631031000000d0"); // EG2
    cuts.AddCutCalo("8008f103","411793105f032230000","01631031000000d0"); // EG1
  } else if (trainConfig == 2045){ // T0 based triggers no TM
    cuts.AddCutCalo("80011103","4117931050032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("8008g103","4117931050032230000","01631031000000d0"); // EG2
    cuts.AddCutCalo("8008f103","4117931050032230000","01631031000000d0"); // EG1

  } else if (trainConfig == 2046){ // copy of 2030 for clusterizer cell timing studies NO TRACK MATCHING
    cuts.AddCutCalo("80010123","4117931050032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("8008e123","4117931050032230000","01631031000000d0"); // EG2
    cuts.AddCutCalo("8008d123","4117931050032230000","01631031000000d0"); // EG1
  } else if (trainConfig == 2047){ // copy of 2030 for clusterizer cell timing studies NO TRACK MATCHING
    cuts.AddCutCalo("80010123","4117931050032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("8008e123","4117931050032230000","01631031000000d0"); // EG2
    cuts.AddCutCalo("8008d123","4117931050032230000","01631031000000d0"); // EG1
  } else if (trainConfig == 2048){ // copy of 2030 for clusterizer cell timing studies NO TRACK MATCHING
    cuts.AddCutCalo("80010123","4117931050032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("8008e123","4117931050032230000","01631031000000d0"); // EG2
    cuts.AddCutCalo("8008d123","4117931050032230000","01631031000000d0"); // EG1
  } else if (trainConfig == 2049){ // copy of 2030 for clusterizer cell timing studies NO TRACK MATCHING
    cuts.AddCutCalo("80010123","4117931050032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("8008e123","4117931050032230000","01631031000000d0"); // EG2
    cuts.AddCutCalo("8008d123","4117931050032230000","01631031000000d0"); // EG1

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

  } else if (trainConfig == 2060){ // TB NL tests 100 MeV aggregation
    cuts.AddCutCalo("80010123","411796105f032230000","01631031000000d0"); // INT7
  } else if (trainConfig == 2061){ // TB NL tests 100 MeV aggregation
    cuts.AddCutCalo("8008e123","411796105f032230000","01631031000000d0"); // EG2
    cuts.AddCutCalo("8008d123","411796105f032230000","01631031000000d0"); // EG1
  } else if (trainConfig == 2062){ // TB NL tests 50 MeV aggregation
    cuts.AddCutCalo("80010123","411796205f032230000","01631031000000d0"); // INT7
  } else if (trainConfig == 2063){ // TB NL tests 50 MeV aggregation
    cuts.AddCutCalo("8008e123","411796205f032230000","01631031000000d0"); // EG2
    cuts.AddCutCalo("8008d123","411796205f032230000","01631031000000d0"); // EG1
  } else if (trainConfig == 2064){ // TB NL tests 150 MeV aggregation
    cuts.AddCutCalo("80010123","411796305f032230000","01631031000000d0"); // INT7
  } else if (trainConfig == 2065){ // TB NL tests 150 MeV aggregation
    cuts.AddCutCalo("8008e123","411796305f032230000","01631031000000d0"); // EG2
    cuts.AddCutCalo("8008d123","411796305f032230000","01631031000000d0"); // EG1
  } else if (trainConfig == 2066){ // TB NL tests 300 MeV aggregation
    cuts.AddCutCalo("80010123","411796405f032230000","01631031000000d0"); // INT7
  } else if (trainConfig == 2067){ // TB NL tests 300 MeV aggregation
    cuts.AddCutCalo("8008e123","411796405f032230000","01631031000000d0"); // EG2
    cuts.AddCutCalo("8008d123","411796405f032230000","01631031000000d0"); // EG1
  } else if (trainConfig == 2068){ // TB NL tests 300 MeV aggregation
    cuts.AddCutCalo("80010123","411796905f032230000","01631031000000d0"); // INT7
  } else if (trainConfig == 2069){ // TB NL tests 300 MeV aggregation
    cuts.AddCutCalo("8008e123","411796905f032230000","01631031000000d0"); // EG2
    cuts.AddCutCalo("8008d123","411796905f032230000","01631031000000d0"); // EG1


  } else if (trainConfig == 2070){ // EMC only
    cuts.AddCutCalo("80010103","411793405f032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("8008e103","411793405f032230000","01631031000000d0"); // EG2
    cuts.AddCutCalo("8008d103","411793405f032230000","01631031000000d0"); // EG1
  } else if (trainConfig == 2071){ // EMC only
    cuts.AddCutCalo("80010103","411793205f032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("8008e103","411793205f032230000","01631031000000d0"); // EG2
    cuts.AddCutCalo("8008d103","411793205f032230000","01631031000000d0"); // EG1
  } else if (trainConfig == 2072){ // EMC only
    cuts.AddCutCalo("80010103","411793305f032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("8008e103","411793305f032230000","01631031000000d0"); // EG2
    cuts.AddCutCalo("8008d103","411793305f032230000","01631031000000d0"); // EG1
  } else if (trainConfig == 2073){ // EMC only
    cuts.AddCutCalo("80010103","411793805f032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("8008e103","411793805f032230000","01631031000000d0"); // EG2
    cuts.AddCutCalo("8008d103","411793805f032230000","01631031000000d0"); // EG1
  } else if (trainConfig == 2074){ // EMC only
    cuts.AddCutCalo("80010103","411793905f032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("8008e103","411793905f032230000","01631031000000d0"); // EG2
    cuts.AddCutCalo("8008d103","411793905f032230000","01631031000000d0"); // EG1

  } else if (trainConfig == 2075){ // EMC only
    cuts.AddCutCalo("80010103","111113405f032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("80085103","111113405f032230000","01631031000000d0"); // EG2
    cuts.AddCutCalo("80083103","111113405f032230000","01631031000000d0"); // EG1
  } else if (trainConfig == 2076){ // EMC only
    cuts.AddCutCalo("80010103","111113205f032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("80085103","111113205f032230000","01631031000000d0"); // EG2
    cuts.AddCutCalo("80083103","111113205f032230000","01631031000000d0"); // EG1
  } else if (trainConfig == 2077){ // EMC only
    cuts.AddCutCalo("80010103","111113305f032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("80085103","111113305f032230000","01631031000000d0"); // EG2
    cuts.AddCutCalo("80083103","111113305f032230000","01631031000000d0"); // EG1
  } else if (trainConfig == 2078){ // EMC only
    cuts.AddCutCalo("80010103","111113805f032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("80085103","111113805f032230000","01631031000000d0"); // EG2
    cuts.AddCutCalo("80083103","111113805f032230000","01631031000000d0"); // EG1
  } else if (trainConfig == 2079){ // EMC only
    cuts.AddCutCalo("80010103","111113905f032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("80085103","111113905f032230000","01631031000000d0"); // EG2
    cuts.AddCutCalo("80083103","111113905f032230000","01631031000000d0"); // EG1

  } else if (trainConfig == 2080){ // EMC only
    cuts.AddCutCalo("80010103","4117934050032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("8008e103","4117934050032230000","01631031000000d0"); // EG2
    cuts.AddCutCalo("8008d103","4117934050032230000","01631031000000d0"); // EG1
  } else if (trainConfig == 2081){ // EMC only
    cuts.AddCutCalo("80010103","4117932050032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("8008e103","4117932050032230000","01631031000000d0"); // EG2
    cuts.AddCutCalo("8008d103","4117932050032230000","01631031000000d0"); // EG1
  } else if (trainConfig == 2082){ // EMC only
    cuts.AddCutCalo("80010103","4117933050032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("8008e103","4117933050032230000","01631031000000d0"); // EG2
    cuts.AddCutCalo("8008d103","4117933050032230000","01631031000000d0"); // EG1
  } else if (trainConfig == 2083){ // EMC only
    cuts.AddCutCalo("80010103","4117938050032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("8008e103","4117938050032230000","01631031000000d0"); // EG2
    cuts.AddCutCalo("8008d103","4117938050032230000","01631031000000d0"); // EG1
  } else if (trainConfig == 2084){ // EMC only
    cuts.AddCutCalo("80010103","4117939050032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("8008e103","4117939050032230000","01631031000000d0"); // EG2
    cuts.AddCutCalo("8008d103","4117939050032230000","01631031000000d0"); // EG1

  } else if (trainConfig == 2085){ // EMC only
    cuts.AddCutCalo("80010103","1111131050032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("80085103","1111131050032230000","01631031000000d0"); // EG2
    cuts.AddCutCalo("80083103","1111131050032230000","01631031000000d0"); // EG1
  } else if (trainConfig == 2086){ // EMC only
    cuts.AddCutCalo("80010103","1111132050032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("80085103","1111132050032230000","01631031000000d0"); // EG2
    cuts.AddCutCalo("80083103","1111132050032230000","01631031000000d0"); // EG1
  } else if (trainConfig == 2087){ // EMC only
    cuts.AddCutCalo("80010103","1111133050032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("80085103","1111133050032230000","01631031000000d0"); // EG2
    cuts.AddCutCalo("80083103","1111133050032230000","01631031000000d0"); // EG1
  } else if (trainConfig == 2088){ // EMC only
    cuts.AddCutCalo("80010103","1111138050032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("80085103","1111138050032230000","01631031000000d0"); // EG2
    cuts.AddCutCalo("80083103","1111138050032230000","01631031000000d0"); // EG1
  } else if (trainConfig == 2089){ // EMC only
    cuts.AddCutCalo("80010103","1111139050032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("80085103","1111139050032230000","01631031000000d0"); // EG2
    cuts.AddCutCalo("80083103","1111139050032230000","01631031000000d0"); // EG1


  } else if (trainConfig == 2090){ // EMC only
    cuts.AddCutCalo("80010103","1111131050032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("80085103","1111131050032230000","01631031000000d0"); // EG2
    cuts.AddCutCalo("80083103","1111131050032230000","01631031000000d0"); // EG1
  } else if (trainConfig == 2091){ // DMC only
    cuts.AddCutCalo("80010103","3885531050032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("8008b103","3885531050032230000","01631031000000d0"); // EG2
    cuts.AddCutCalo("80089103","3885531050032230000","01631031000000d0"); // EG1
  } else if (trainConfig == 2092){ // EDC
    cuts.AddCutCalo("80010103","4117931050032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("8008e103","4117931050032230000","01631031000000d0"); // EG2
    cuts.AddCutCalo("8008d103","4117931050032230000","01631031000000d0"); // EG1
  } else if (trainConfig == 2093){ // EMC only
    cuts.AddCutCalo("80010103","1111132050032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("80085103","1111132050032230000","01631031000000d0"); // EG2
    cuts.AddCutCalo("80083103","1111132050032230000","01631031000000d0"); // EG1
  } else if (trainConfig == 2094){ // EMC only
    cuts.AddCutCalo("80010103","111113205f032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("80085103","111113205f032230000","01631031000000d0"); // EG2
    cuts.AddCutCalo("80083103","111113205f032230000","01631031000000d0"); // EG1
  } else if (trainConfig == 2095){ // EMC only
    cuts.AddCutCalo("80010103","111113805f032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("80085103","111113805f032230000","01631031000000d0"); // EG2
    cuts.AddCutCalo("80083103","111113805f032230000","01631031000000d0"); // EG1
  } else if (trainConfig == 2096){ // EMC only
    cuts.AddCutCalo("80010103","111113905f032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("80085103","111113905f032230000","01631031000000d0"); // EG2
    cuts.AddCutCalo("80083103","111113905f032230000","01631031000000d0"); // EG1

  } else if (trainConfig == 2099){ // no NCell cut
    cuts.AddCutCalo("80010103","4117931050030230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("8008e103","4117931050030230000","01631031000000d0"); // EG2
    cuts.AddCutCalo("8008d103","4117931050030230000","01631031000000d0"); // EG1


  } else if (trainConfig == 2100) { // CALO variations
    cuts.AddCutCalo("80010103","111113205f03h230000","01631031000000d0"); // NCell cut effi MC (param 1)
    cuts.AddCutCalo("80010103","111113205f03i230000","01631031000000d0"); // NCell cut effi MC (param 2)
    cuts.AddCutCalo("80010103","111113205f03j230000","01631031000000d0"); // NCell cut effi data (param 1)
    cuts.AddCutCalo("80010103","111113205f03k230000","01631031000000d0"); // NCell cut effi data (param 2)
  } else if (trainConfig == 2101) { // CALO variations
    cuts.AddCutCalo("80085103","111113205f03h230000","01631031000000d0"); // NCell cut effi MC (param 1)
    cuts.AddCutCalo("80085103","111113205f03i230000","01631031000000d0"); // NCell cut effi MC (param 2)
    cuts.AddCutCalo("80085103","111113205f03j230000","01631031000000d0"); // NCell cut effi data (param 1)
    cuts.AddCutCalo("80085103","111113205f03k230000","01631031000000d0"); // NCell cut effi data (param 2)
  } else if (trainConfig == 2102) { // CALO variations
    cuts.AddCutCalo("80083103","111113205f03h230000","01631031000000d0"); // NCell cut effi MC (param 1)
    cuts.AddCutCalo("80083103","111113205f03i230000","01631031000000d0"); // NCell cut effi MC (param 2)
    cuts.AddCutCalo("80083103","111113205f03j230000","01631031000000d0"); // NCell cut effi data (param 1)
    cuts.AddCutCalo("80083103","111113205f03k230000","01631031000000d0"); // NCell cut effi data (param 2)

  // systematics for pPb8TeV PRL
  } else if (trainConfig == 2200) { // CALO variations
    cuts.AddCutCalo("80010103","411793205f032230000","01631031000000d0"); // std
    cuts.AddCutCalo("80010103","411793205f022230000","01631031000000d0"); // 0.6 GeV/c
    cuts.AddCutCalo("80010103","411793205f042230000","01631031000000d0"); // 0.8 GeV/c
  } else if (trainConfig == 2201) {
    cuts.AddCutCalo("80010103","411793205f031230000","01631031000000d0"); // n cells >= 1
    cuts.AddCutCalo("80010103","411793205f033230000","01631031000000d0"); // n cells >= 3
    cuts.AddCutCalo("80010103","411793205f032200000","01631031000000d0"); // no max M02 cut
    cuts.AddCutCalo("80010103","411793205f032250000","01631031000000d0"); // M02 < 0.3
    cuts.AddCutCalo("80010103","411793205f0322k0000","01631031000000d0"); // M02, pT-dep
  } else if (trainConfig == 2202) {
    cuts.AddCutCalo("80010103","411793205e032230000","01631031000000d0"); // TM var e/p 2.0
    cuts.AddCutCalo("80010103","411793205g032230000","01631031000000d0"); // TM var e/p 1.5
    cuts.AddCutCalo("80010103","411793205h032230000","01631031000000d0"); // TM var e/p 1.25
    cuts.AddCutCalo("80010103","4117932057032230000","01631031000000d0"); // TM var no veto
  } else if (trainConfig == 2203) {
    cuts.AddCutCalo("80010103","411793305f032230000","01631031000000d0"); // NL 33
    cuts.AddCutCalo("80010103","411793405f032230000","01631031000000d0"); // NL 34
    cuts.AddCutCalo("80010103","411793805f032230000","01631031000000d0"); // NL 34
    cuts.AddCutCalo("80010103","411793905f032230000","01631031000000d0"); // NL 34
  } else if (trainConfig == 2204) {
    cuts.AddCutCalo("80010103","411793206f032230000","01631031000000d0"); // 30/35ns
    cuts.AddCutCalo("80010103","411793204f032230000","01631031000000d0"); // 100ns
    cuts.AddCutCalo("80010103","411793207f032230000","01631031000000d0"); // 30ns
  } else if (trainConfig == 2205) {
    cuts.AddCutCalo("80010103","411793205f032230000","01631031000000b0"); // min opening angle 0.0152
    cuts.AddCutCalo("80010103","411793205f032230000","01631031000000c0"); // min opening angle 0.016
    cuts.AddCutCalo("80010103","411793205f032230000","01631031000000e0"); // min opening angle 0.018
    cuts.AddCutCalo("80010103","411793205f032230000","01631031000000f0"); // min opening angle 0.019
  } else if (trainConfig == 2206) {
    cuts.AddCutCalo("80010103","411793205f032230000","01633031000000d0"); // rapidity variation  y<0.6
    cuts.AddCutCalo("80010103","411793205f032230000","01634031000000d0"); // rapidity variation  y<0.5
    cuts.AddCutCalo("80010103","411793205f032230000","01631061000000d0"); // alpha meson variation 1   0<alpha<0.8
    cuts.AddCutCalo("80010103","411793205f032230000","01631051000000d0"); // alpha meson variation 2  0<alpha<0.75
  } else if (trainConfig == 2207) { // CALO variations
    cuts.AddCutCalo("80010103","411793205f03h230000","01631031000000d0"); // NCell cut effi MC (param 1)
    cuts.AddCutCalo("80010103","411793205f03i230000","01631031000000d0"); // NCell cut effi MC (param 2)
    cuts.AddCutCalo("80010103","411793205f03j230000","01631031000000d0"); // NCell cut effi data (param 1)
    cuts.AddCutCalo("80010103","411793205f03k230000","01631031000000d0"); // NCell cut effi data (param 2)


  } else if (trainConfig == 2210) { // CALO variations
    cuts.AddCutCalo("8008e103","411793205f032230000","01631031000000d0"); // std
    cuts.AddCutCalo("8008e103","411793205f022230000","01631031000000d0"); // 0.6 GeV/c
    cuts.AddCutCalo("8008e103","411793205f042230000","01631031000000d0"); // 0.8 GeV/c
  } else if (trainConfig == 2211) {
    cuts.AddCutCalo("8008e103","411793205f031230000","01631031000000d0"); // n cells >= 1
    cuts.AddCutCalo("8008e103","411793205f033230000","01631031000000d0"); // n cells >= 3
    cuts.AddCutCalo("8008e103","411793205f032200000","01631031000000d0"); // no max M02 cut
    cuts.AddCutCalo("8008e103","411793205f032250000","01631031000000d0"); // M02 < 0.3
    cuts.AddCutCalo("8008e103","411793205f0322k0000","01631031000000d0"); // M02, pT-dep
  } else if (trainConfig == 2212) {
    cuts.AddCutCalo("8008e103","411793205e032230000","01631031000000d0"); // TM var e/p 2.0
    cuts.AddCutCalo("8008e103","411793205g032230000","01631031000000d0"); // TM var e/p 1.5
    cuts.AddCutCalo("8008e103","411793205h032230000","01631031000000d0"); // TM var e/p 1.25
    cuts.AddCutCalo("8008e103","4117932057032230000","01631031000000d0"); // TM var no veto
  } else if (trainConfig == 2213) {
    cuts.AddCutCalo("8008e103","411793305f032230000","01631031000000d0"); // NL 33
    cuts.AddCutCalo("8008e103","411793405f032230000","01631031000000d0"); // NL 34
    cuts.AddCutCalo("8008e103","411793805f032230000","01631031000000d0"); // NL 38
    cuts.AddCutCalo("8008e103","411793905f032230000","01631031000000d0"); // NL 39
  } else if (trainConfig == 2214) {
    cuts.AddCutCalo("8008e103","411793206f032230000","01631031000000d0"); // 30/35ns
    cuts.AddCutCalo("8008e103","411793204f032230000","01631031000000d0"); // 100ns
    cuts.AddCutCalo("8008e103","411793207f032230000","01631031000000d0"); // 30ns
  } else if (trainConfig == 2215) {
    cuts.AddCutCalo("8008e103","411793205f032230000","01631031000000b0"); // min opening angle 0.0152
    cuts.AddCutCalo("8008e103","411793205f032230000","01631031000000c0"); // min opening angle 0.016
    cuts.AddCutCalo("8008e103","411793205f032230000","01631031000000e0"); // min opening angle 0.018
    cuts.AddCutCalo("8008e103","411793205f032230000","01631031000000f0"); // min opening angle 0.019
  } else if (trainConfig == 2216) {
    cuts.AddCutCalo("8008e103","411793205f032230000","01633031000000d0"); // rapidity variation  y<0.6
    cuts.AddCutCalo("8008e103","411793205f032230000","01634031000000d0"); // rapidity variation  y<0.5
    cuts.AddCutCalo("8008e103","411793205f032230000","01631061000000d0"); // alpha meson variation 1   0<alpha<0.8
    cuts.AddCutCalo("8008e103","411793205f032230000","01631051000000d0"); // alpha meson variation 2  0<alpha<0.75

  } else if (trainConfig == 2220) { // CALO variations
    cuts.AddCutCalo("8008d103","411793205f032230000","01631031000000d0"); // std
    cuts.AddCutCalo("8008d103","411793205f022230000","01631031000000d0"); // 0.6 GeV/c
    cuts.AddCutCalo("8008d103","411793205f042230000","01631031000000d0"); // 0.8 GeV/c
  } else if (trainConfig == 2221) {
    cuts.AddCutCalo("8008d103","411793205f031230000","01631031000000d0"); // n cells >= 1
    cuts.AddCutCalo("8008d103","411793205f033230000","01631031000000d0"); // n cells >= 3
    cuts.AddCutCalo("8008d103","411793205f032200000","01631031000000d0"); // no max M02 cut
    cuts.AddCutCalo("8008d103","411793205f032250000","01631031000000d0"); // M02 < 0.3
    cuts.AddCutCalo("8008d103","411793205f0322k0000","01631031000000d0"); // M02, pT-dep
  } else if (trainConfig == 2222) {
    cuts.AddCutCalo("8008d103","411793205e032230000","01631031000000d0"); // TM var e/p 2.0
    cuts.AddCutCalo("8008d103","411793205g032230000","01631031000000d0"); // TM var e/p 1.5
    cuts.AddCutCalo("8008d103","411793205h032230000","01631031000000d0"); // TM var e/p 1.25
    cuts.AddCutCalo("8008d103","4117932057032230000","01631031000000d0"); // TM var no veto
  } else if (trainConfig == 2223) {
    cuts.AddCutCalo("8008d103","411793305f032230000","01631031000000d0"); // NL 33
    cuts.AddCutCalo("8008d103","411793405f032230000","01631031000000d0"); // NL 34
    cuts.AddCutCalo("8008d103","411793805f032230000","01631031000000d0"); // NL 38
    cuts.AddCutCalo("8008d103","411793905f032230000","01631031000000d0"); // NL 39
  } else if (trainConfig == 2224) {
    cuts.AddCutCalo("8008d103","411793206f032230000","01631031000000d0"); // 30/35ns
    cuts.AddCutCalo("8008d103","411793204f032230000","01631031000000d0"); // 100ns
    cuts.AddCutCalo("8008d103","411793207f032230000","01631031000000d0"); // 30ns
  } else if (trainConfig == 2225) {
    cuts.AddCutCalo("8008d103","411793205f032230000","01631031000000b0"); // min opening angle 0.0152
    cuts.AddCutCalo("8008d103","411793205f032230000","01631031000000c0"); // min opening angle 0.016
    cuts.AddCutCalo("8008d103","411793205f032230000","01631031000000e0"); // min opening angle 0.018
    cuts.AddCutCalo("8008d103","411793205f032230000","01631031000000f0"); // min opening angle 0.019
  } else if (trainConfig == 2226) {
    cuts.AddCutCalo("8008d103","411793205f032230000","01633031000000d0"); // rapidity variation  y<0.6
    cuts.AddCutCalo("8008d103","411793205f032230000","01634031000000d0"); // rapidity variation  y<0.5
    cuts.AddCutCalo("8008d103","411793205f032230000","01631061000000d0"); // alpha meson variation 1   0<alpha<0.8
    cuts.AddCutCalo("8008d103","411793205f032230000","01631051000000d0"); // alpha meson variation 2  0<alpha<0.75


  // systematics for pPb8TeV PRL without TM
  } else if (trainConfig == 2300) { // CALO variations
    cuts.AddCutCalo("80010103","4117932050032230000","01631031000000d0"); // std
    cuts.AddCutCalo("80010103","4117932050022230000","01631031000000d0"); // 0.6 GeV/c
    cuts.AddCutCalo("80010103","4117932050042230000","01631031000000d0"); // 0.8 GeV/c
  } else if (trainConfig == 2301) {
    cuts.AddCutCalo("80010103","4117932050031230000","01631031000000d0"); // n cells >= 1
    cuts.AddCutCalo("80010103","4117932050033230000","01631031000000d0"); // n cells >= 3
    cuts.AddCutCalo("80010103","4117932050032200000","01631031000000d0"); // no max M02 cut
    cuts.AddCutCalo("80010103","4117932050032250000","01631031000000d0"); // M02 < 0.3
    cuts.AddCutCalo("80010103","41179320500322k0000","01631031000000d0"); // M02, pT-dep
  } else if (trainConfig == 2302) {
    cuts.AddCutCalo("80010103","4117933050032230000","01631031000000d0"); // NL 33
    cuts.AddCutCalo("80010103","4117934050032230000","01631031000000d0"); // NL 34
    cuts.AddCutCalo("80010103","4117938050032230000","01631031000000d0"); // NL 38
    cuts.AddCutCalo("80010103","4117939050032230000","01631031000000d0"); // NL 39
  } else if (trainConfig == 2303) {
    cuts.AddCutCalo("80010103","4117932060032230000","01631031000000d0"); // 30/35ns
    cuts.AddCutCalo("80010103","4117932040032230000","01631031000000d0"); // 100ns
    cuts.AddCutCalo("80010103","4117932070032230000","01631031000000d0"); // 30ns
  } else if (trainConfig == 2304) {
    cuts.AddCutCalo("80010103","4117932050032230000","01631031000000b0"); // min opening angle 0.0152
    cuts.AddCutCalo("80010103","4117932050032230000","01631031000000c0"); // min opening angle 0.016
    cuts.AddCutCalo("80010103","4117932050032230000","01631031000000e0"); // min opening angle 0.018
    cuts.AddCutCalo("80010103","4117932050032230000","01631031000000f0"); // min opening angle 0.019
  } else if (trainConfig == 2305) {
    cuts.AddCutCalo("80010103","4117932050032230000","01633031000000d0"); // rapidity variation  y<0.6
    cuts.AddCutCalo("80010103","4117932050032230000","01634031000000d0"); // rapidity variation  y<0.5
    cuts.AddCutCalo("80010103","4117932050032230000","01631061000000d0"); // alpha meson variation 1   0<alpha<0.8
    cuts.AddCutCalo("80010103","4117932050032230000","01631051000000d0"); // alpha meson variation 2  0<alpha<0.75


  } else if (trainConfig == 2310) { // CALO variations no TM
    cuts.AddCutCalo("8008e103","4117932050032230000","01631031000000d0"); // std
    cuts.AddCutCalo("8008e103","4117932050022230000","01631031000000d0"); // 0.6 GeV/c
    cuts.AddCutCalo("8008e103","4117932050042230000","01631031000000d0"); // 0.8 GeV/c
  } else if (trainConfig == 2311) {
    cuts.AddCutCalo("8008e103","4117932050031230000","01631031000000d0"); // n cells >= 1
    cuts.AddCutCalo("8008e103","4117932050033230000","01631031000000d0"); // n cells >= 3
    cuts.AddCutCalo("8008e103","4117932050032200000","01631031000000d0"); // no max M02 cut
    cuts.AddCutCalo("8008e103","4117932050032250000","01631031000000d0"); // M02 < 0.3
    cuts.AddCutCalo("8008e103","41179320500322k0000","01631031000000d0"); // M02, pT-dep
  } else if (trainConfig == 2312) {
    cuts.AddCutCalo("8008e103","4117933050032230000","01631031000000d0"); // NL 33
    cuts.AddCutCalo("8008e103","4117934050032230000","01631031000000d0"); // NL 34
    cuts.AddCutCalo("8008e103","4117938050032230000","01631031000000d0"); // NL 38
    cuts.AddCutCalo("8008e103","4117939050032230000","01631031000000d0"); // NL 39
  } else if (trainConfig == 2313) {
    cuts.AddCutCalo("8008e103","4117932060032230000","01631031000000d0"); // 30/35ns
    cuts.AddCutCalo("8008e103","4117932040032230000","01631031000000d0"); // 100ns
    cuts.AddCutCalo("8008e103","4117932070032230000","01631031000000d0"); // 30ns
  } else if (trainConfig == 2314) {
    cuts.AddCutCalo("8008e103","4117932050032230000","01631031000000b0"); // min opening angle 0.0152
    cuts.AddCutCalo("8008e103","4117932050032230000","01631031000000c0"); // min opening angle 0.016
    cuts.AddCutCalo("8008e103","4117932050032230000","01631031000000e0"); // min opening angle 0.018
    cuts.AddCutCalo("8008e103","4117932050032230000","01631031000000f0"); // min opening angle 0.019
  } else if (trainConfig == 2315) {
    cuts.AddCutCalo("8008e103","4117932050032230000","01633031000000d0"); // rapidity variation  y<0.6
    cuts.AddCutCalo("8008e103","4117932050032230000","01634031000000d0"); // rapidity variation  y<0.5
    cuts.AddCutCalo("8008e103","4117932050032230000","01631061000000d0"); // alpha meson variation 1   0<alpha<0.8
    cuts.AddCutCalo("8008e103","4117932050032230000","01631051000000d0"); // alpha meson variation 2  0<alpha<0.75

  } else if (trainConfig == 2320) { // CALO variations
    cuts.AddCutCalo("8008d103","4117932050032230000","01631031000000d0"); // std
    cuts.AddCutCalo("8008d103","4117932050022230000","01631031000000d0"); // 0.6 GeV/c
    cuts.AddCutCalo("8008d103","4117932050042230000","01631031000000d0"); // 0.8 GeV/c
  } else if (trainConfig == 2321) {
    cuts.AddCutCalo("8008d103","4117932050031230000","01631031000000d0"); // n cells >= 1
    cuts.AddCutCalo("8008d103","4117932050033230000","01631031000000d0"); // n cells >= 3
    cuts.AddCutCalo("8008d103","4117932050032200000","01631031000000d0"); // no max M02 cut
    cuts.AddCutCalo("8008d103","4117932050032250000","01631031000000d0"); // M02 < 0.3
    cuts.AddCutCalo("8008d103","41179320500322k0000","01631031000000d0"); // M02, pT-dep
  } else if (trainConfig == 2322) {
    cuts.AddCutCalo("8008d103","4117933050032230000","01631031000000d0"); // NL 33
    cuts.AddCutCalo("8008d103","4117934050032230000","01631031000000d0"); // NL 34
    cuts.AddCutCalo("8008d103","4117938050032230000","01631031000000d0"); // NL 38
    cuts.AddCutCalo("8008d103","4117939050032230000","01631031000000d0"); // NL 39
  } else if (trainConfig == 2323) {
    cuts.AddCutCalo("8008d103","4117932060032230000","01631031000000d0"); // 30/35ns
    cuts.AddCutCalo("8008d103","4117932040032230000","01631031000000d0"); // 100ns
    cuts.AddCutCalo("8008d103","4117932070032230000","01631031000000d0"); // 30ns
  } else if (trainConfig == 2324) {
    cuts.AddCutCalo("8008d103","4117932050032230000","01631031000000b0"); // min opening angle 0.0152
    cuts.AddCutCalo("8008d103","4117932050032230000","01631031000000c0"); // min opening angle 0.016
    cuts.AddCutCalo("8008d103","4117932050032230000","01631031000000e0"); // min opening angle 0.018
    cuts.AddCutCalo("8008d103","4117932050032230000","01631031000000f0"); // min opening angle 0.019
  } else if (trainConfig == 2325) {
    cuts.AddCutCalo("8008d103","4117932050032230000","01633031000000d0"); // rapidity variation  y<0.6
    cuts.AddCutCalo("8008d103","4117932050032230000","01634031000000d0"); // rapidity variation  y<0.5
    cuts.AddCutCalo("8008d103","4117932050032230000","01631061000000d0"); // alpha meson variation 1   0<alpha<0.8
    cuts.AddCutCalo("8008d103","4117932050032230000","01631051000000d0"); // alpha meson variation 2  0<alpha<0.75

  // PHOS pPb8TeV
  } else if (trainConfig == 3000){ // settings Jens
    cuts.AddCutCalo("80010103","24466190ua012300000","0163103100000010"); // INT7
    cuts.AddCutCalo("80062103","24466190ua012300000","0163103100000010"); // PHI7

  } else if (trainConfig == 3010){ // settings Dmitri
    cuts.AddCutCalo("80010103","24466640ya09dc00000","0163103100000010"); // INT7
  } else if (trainConfig == 3011){ // settings Dmitri
    cuts.AddCutCalo("80062103","24466640ya09dc00000","0163103100000010"); // PHI7
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
    analysisEventCuts[i]->SetLightOutput(enableLightOutput);
    if (periodNameV0Reader.CompareTo("") != 0) analysisEventCuts[i]->SetPeriodEnum(periodNameV0Reader);
    analysisEventCuts[i]->InitializeCutsFromCutString((cuts.GetEventCut(i)).Data());
    EventCutList->Add(analysisEventCuts[i]);
    analysisEventCuts[i]->SetFillCutHistograms("",kFALSE);
    analysisEventCuts[i]->SetDebugLevel(debugLevel);

    analysisClusterCuts[i] = new AliCaloPhotonCuts(isMC);
    analysisClusterCuts[i]->SetHistoToModifyAcceptance(histoAcc);
    analysisClusterCuts[i]->SetV0ReaderName(V0ReaderName);
    analysisClusterCuts[i]->SetCorrectionTaskSetting(corrTaskSetting);
    analysisClusterCuts[i]->SetCaloTrackMatcherName(TrackMatcherName);
    analysisClusterCuts[i]->SetLightOutput(enableLightOutput);
    analysisClusterCuts[i]->InitializeCutsFromCutString((cuts.GetClusterCut(i)).Data());
    ClusterCutList->Add(analysisClusterCuts[i]);
    analysisClusterCuts[i]->SetExtendedMatchAndQA(enableExtMatchAndQA);

    analysisMesonCuts[i] = new AliConversionMesonCuts();
    analysisMesonCuts[i]->SetLightOutput(enableLightOutput);
    analysisMesonCuts[i]->InitializeCutsFromCutString((cuts.GetMesonCut(i)).Data());
    analysisMesonCuts[i]->SetIsMergedClusterCut(2);
    analysisMesonCuts[i]->SetCaloMesonCutsObject(analysisClusterCuts[i]);
    MesonCutList->Add(analysisMesonCuts[i]);
    analysisMesonCuts[i]->SetFillCutHistograms("");
    analysisEventCuts[i]->SetAcceptedHeader(HeaderList);

    if(analysisMesonCuts[i]->DoGammaSwappForBg()) analysisClusterCuts[i]->SetUseEtaPhiMapForBackCand(kTRUE);
    analysisClusterCuts[i]->SetFillCutHistograms("");
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
