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
//pp together with all supporting classes
//***************************************************************************************

//***************************************************************************************
//main function
//***************************************************************************************
void AddTask_GammaCalo_pp(
  Int_t     trainConfig                   = 1,        // change different set of cuts
  Int_t     isMC                          = 0,        // run MC
  TString   photonCutNumberV0Reader       = "",       // 00000008400000000100000000 nom. B, 00000088400000000100000000 low B
  TString   periodNameV0Reader            = "",
  // general setting for task
  Int_t     enableQAMesonTask             = 0,        // enable QA in AliAnalysisTaskGammaConvV1
  Int_t     enableQAClusterTask           = 0,        // enable additional QA task
  Int_t     enableExtMatchAndQA           = 0,        // disabled (0), extMatch (1), extQA_noCellQA (2), extMatch+extQA_noCellQA (3), extQA+cellQA (4), extMatch+extQA+cellQA (5)
  Int_t     enableLightOutput             = kFALSE,   // switch to run light output (only essential histograms for afterburner)
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

  AliCutHandlerPCM cuts(13);


  TString addTaskName                 = "AddTask_GammaCalo_pp";
  TString sAdditionalTrainConfig      = cuts.GetSpecialSettingFromAddConfig(additionalTrainConfig, "", "", addTaskName);
  if (sAdditionalTrainConfig.Atoi() > 0){
    trainConfig = trainConfig + sAdditionalTrainConfig.Atoi();
    cout << "INFO: " << addTaskName.Data() << " running additionalTrainConfig '" << sAdditionalTrainConfig.Atoi() << "', train config: '" << trainConfig << "'" << endl;
  }

  TString fileNamePtWeights           = cuts.GetSpecialFileNameFromString (fileNameExternalInputs, "FPTW:");
  TString fileNameMultWeights         = cuts.GetSpecialFileNameFromString (fileNameExternalInputs, "FMUW:");
  TString fileNameCustomTriggerMimicOADB   = cuts.GetSpecialFileNameFromString (fileNameExternalInputs, "FTRM:");

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
  if(rmaxFacPtHardSetting->GetEntries()<1){cout << "ERROR: AddTask_GammaCalo_pp during parsing of settingMaxFacPtHard String '" << settingMaxFacPtHard.Data() << "'" << endl; return;}
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

  Int_t isHeavyIon = 0;

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
  TString cutnumberEvent = "00000003";
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
  if (enableLightOutput > 0) task->SetLightOutput(kTRUE);
  if (enableLightOutput == 5) task->SetECalibOutput(kTRUE);
  task->SetCorrectionTaskSetting(corrTaskSetting);
  task->SetTrackMatcherRunningMode(trackMatcherRunningMode);

  // here is the order of the cluster cut string
  // usually for EMCal we start with 11111: default values for            "ClusterType", "EtaMin", "EtaMax", "PhiMin", "PhiMax"
  // then two numbers for nonlinearity, e.g. 21: this is                  "NonLinearity1", "NonLinearity2"
  // Then some cuts on the clusters, e.g. 06003222: this is               "DistanceToBadChannel", "Timing", "TrackMatching", "ExoticCell", "MinEnergy", "MinNCells", "MinM02", "MaxM02"
  // finally some for now unused cuts, usually 0000: this is              "MinMaxM20", "RecConv", "MaximumDispersion", "NLM"


  // *****************************************************************************************************
  // ******************** pp 2.76 TeV cuts paper EMC *****************************************************
  // *****************************************************************************************************
  if (trainConfig == 1){ // pp 2.76 TeV LHC11a paper cuts
    cuts.AddCutCalo("00003013","1111121057032220000","0163103100000050"); // MB w/o pileup
    cuts.AddCutCalo("00003113","1111121057032220000","0163103100000050"); // MB
    cuts.AddCutCalo("00051013","1111121057032220000","0163103100000050"); // EMC1
  } else if (trainConfig == 2){  // pp 2.76 TeV LHC13g paper cuts
    cuts.AddCutCalo("00010113","1111121067032220000","0163103100000050");
    cuts.AddCutCalo("00010013","1111121067032220000","0163103100000050"); // without pile-up correction
    cuts.AddCutCalo("00052013","1111121067032220000","0163103100000050"); // EMC7
    cuts.AddCutCalo("00083013","1111121067032220000","0163103100000050"); // EMCEG1,
    cuts.AddCutCalo("00085013","1111121067032220000","0163103100000050"); // EMCEG2,

  } else if (trainConfig == 3){
    cuts.AddCutCalo("00003113","1111121057032250000","0163103100000050"); // MB
    cuts.AddCutCalo("00003113","1111121057032260000","0163103100000050"); // MB
    cuts.AddCutCalo("00003113","1111121057032240000","0163103100000050"); // MB
    cuts.AddCutCalo("00003113","1111121057032290000","0163103100000050"); // MB
  } else if (trainConfig == 4){
    cuts.AddCutCalo("00003113","11111210570322a0000","0163103100000050"); // MB
    cuts.AddCutCalo("00003113","11111210570322b0000","0163103100000050"); // MB
    cuts.AddCutCalo("00003113","11111210570322c0000","0163103100000050"); // MB
  } else if (trainConfig == 5){
    cuts.AddCutCalo("00003113","1111121057032250000","01631031000000d0"); // MB
    cuts.AddCutCalo("00003113","1111121057032260000","01631031000000d0"); // MB
    cuts.AddCutCalo("00003113","1111121057032240000","01631031000000d0"); // MB
    cuts.AddCutCalo("00003113","1111121057032290000","01631031000000d0"); // MB
  } else if (trainConfig == 6){
    cuts.AddCutCalo("00003113","11111210570322a0000","01631031000000d0"); // MB
    cuts.AddCutCalo("00003113","11111210570322b0000","01631031000000d0"); // MB
    cuts.AddCutCalo("00003113","11111210570322c0000","01631031000000d0"); // MB

  // *****************************************************************************************************
  // ************************************* Calibration configuration EMC *********************************
  // *****************************************************************************************************
  } else if (trainConfig == 40){ // EMCAL clusters 2.76 TeV LHC11a, with SDD (0), kEMC1 (1) with TM
    cuts.AddCutCalo("00003113","1111100053032220000","0163103100000050"); // MB
    cuts.AddCutCalo("00051013","1111100053032220000","0163103100000050"); // EMC1
  } else if (trainConfig == 41){ // EMCAL clusters 2.76 TeV LHC11a, with SDD (0), kEMC1 (1) without TM
    cuts.AddCutCalo("00003113","1111100050032220000","0163103100000050"); // 700 MeV cluster min energy
    cuts.AddCutCalo("00051013","1111100050032220000","0163103100000050"); // 700 MeV cluster min energy
  } else if (trainConfig == 42){  // EMCAL clusters 2.76TeV LHC13g with TM
    cuts.AddCutCalo("00010113","1111100063032220000","0163103100000050");
    cuts.AddCutCalo("00052013","1111100063032220000","0163103100000050"); // EMC7
    cuts.AddCutCalo("00085013","1111100063032220000","0163103100000050"); // EG2
    cuts.AddCutCalo("00083013","1111100063032220000","0163103100000050"); // EG1
  } else if (trainConfig == 43){  // EMCAL clusters 2.76TeV LHC13g without TM
    cuts.AddCutCalo("00010113","1111100060032220000","0163103100000050");
    cuts.AddCutCalo("00052013","1111100060032220000","0163103100000050"); // EMC7
    cuts.AddCutCalo("00085013","1111100060032220000","0163103100000050"); // EG2
    cuts.AddCutCalo("00083013","1111100060032220000","0163103100000050"); // EG1
  } else if (trainConfig == 44){   // EMCAL clusters 7TeV LHC10
    cuts.AddCutCalo("00000113","1111100010032220000","0163103100000050"); // wo TM
    cuts.AddCutCalo("00000113","1111100013032220000","0163103100000050"); // w TM
  } else if (trainConfig == 45){  // EMCAL clusters, 8TeV LHC12 with TM
    cuts.AddCutCalo("00010113","1111100067032230000","0163103100000060");
    cuts.AddCutCalo("00052013","1111100067032230000","0163103100000060"); // EMC7
    cuts.AddCutCalo("00081013","1111100067032230000","0163103100000060"); // EMCEGA
  } else if (trainConfig == 46){  // EMCAL clusters, 8TeV LHC12 without TM
    cuts.AddCutCalo("00010113","1111100060032230000","0163103100000060");
    cuts.AddCutCalo("00052013","1111100060032230000","0163103100000060"); // EMC7
    cuts.AddCutCalo("00081013","1111100060032230000","0163103100000060"); // EMCEGA
  // *****************************************************************************************************
  // ************************************* Direct Photon Configurations  *********************************
  // *****************************************************************************************************
  } else if (trainConfig == 60){
    cuts.AddCutCalo("00003113","11111210570322c0000","01631031000000d0"); // MB std
    cuts.AddCutCalo("00003113","1111121057032200000","01631031000000d0"); // no M02
    cuts.AddCutCalo("00003113","11111210570322k0000","01631031000000d0"); // M02, pT-dep with 0.27-0.5
    cuts.AddCutCalo("00003113","11111210570322l0000","01631031000000d0"); // M02, pT-dep with 0.32-0.5
    cuts.AddCutCalo("00003113","11111210570322m0000","01631031000000d0"); // M02, pT-dep with 0.32-0.7
  } else if (trainConfig == 61){
    cuts.AddCutCalo("00003113","11111210570322n0000","01631031000000d0"); // M02, pT-dep with 0.32-0.7
    cuts.AddCutCalo("00003113","11111210570322o0000","01631031000000d0"); // M02, pT-dep with 0.27-0.7
    cuts.AddCutCalo("00003113","11111210570322p0000","01631031000000d0"); // M02, pT-dep with 0.32-0.7
    cuts.AddCutCalo("00003113","11111210570322q0000","01631031000000d0"); // M02, pT-dep with 0.34-0.7
  } else if (trainConfig == 62){
    cuts.AddCutCalo("00051013","11111210570322c0000","01631031000000d0"); // MB std
    cuts.AddCutCalo("00051013","1111121057032200000","01631031000000d0"); // no M02
    cuts.AddCutCalo("00051013","11111210570322k0000","01631031000000d0"); // M02, pT-dep with 0.27-0.5
    cuts.AddCutCalo("00051013","11111210570322l0000","01631031000000d0"); // M02, pT-dep with 0.32-0.5
    cuts.AddCutCalo("00051013","11111210570322m0000","01631031000000d0"); // M02, pT-dep with 0.32-0.7
  } else if (trainConfig == 63){
    cuts.AddCutCalo("00051013","11111210570322n0000","01631031000000d0"); // M02, pT-dep with 0.32-0.7
    cuts.AddCutCalo("00051013","11111210570322o0000","01631031000000d0"); // M02, pT-dep with 0.27-0.7
    cuts.AddCutCalo("00051013","11111210570322p0000","01631031000000d0"); // M02, pT-dep with 0.32-0.7
    cuts.AddCutCalo("00051013","11111210570322q0000","01631031000000d0"); // M02, pT-dep with 0.34-0.7
  } else if (trainConfig == 64){
    cuts.AddCutCalo("00003113","11111210570322l0000","01631031000000d0"); // M02 pt dep with new std: 0.32, 0.0072, 0.5
    cuts.AddCutCalo("00003113","11111210570322d0000","01631031000000d0"); // M02, pt dep with  0.27, 0.0072, 0.4
    cuts.AddCutCalo("00003113","11111210570322e0000","01631031000000d0"); // M02, pt dep with  0.31, 0.0072, 0.5
    cuts.AddCutCalo("00003113","11111210570322f0000","01631031000000d0"); // M02, pt dep with  0.36, 0.0072, 0.7
  } else if (trainConfig == 65){
    cuts.AddCutCalo("00003113","11111210570322m0000","01631031000000d0"); // M02, pt dep with  0.32, 0.0152, 0.5
    cuts.AddCutCalo("00003113","11111210570322g0000","01631031000000d0"); // M02, pt dep with  0.37, 0.0072, 0.7
    cuts.AddCutCalo("00003113","11111210570322h0000","01631031000000d0"); // M02, pT-dep with  0.30, 0.0072, 0.5
    cuts.AddCutCalo("00003113","11111210570322i0000","01631031000000d0"); // M02, pT-dep with  0.35, 0.0072, 0.7
  } else if (trainConfig == 66){
    cuts.AddCutCalo("00003113","11111210570322j0000","01631031000000d0"); // M02, pT-dep with  0.25, 0.0072, 0.39
    cuts.AddCutCalo("00003113","11111210570322r0000","01631031000000d0"); // M02, pT-dep with  0.25, 0.0072, 0.5
    cuts.AddCutCalo("00003113","11111210570322s0000","01631031000000d0"); // M02, pT-dep with  0.32, 0.0238, 0.7
    cuts.AddCutCalo("00003113","11111210570322n0000","01631031000000d0"); // M02, pT-dep with  0.32, 0.0238, 0.7
  // *****************************************************************************************************
  // 8 TeV configs
  // *****************************************************************************************************
    //std, but no opening angle cut
  } else if (trainConfig == 99){ // EMCAL clusters pp 8 TeV
    cuts.AddCutCalo("00010113","1111111067032220000","0163103100000000"); // std
    cuts.AddCutCalo("00052113","1111111067032220000","0163103100000000"); // std
    cuts.AddCutCalo("00081113","1111111067032220000","0163103100000000"); // std
    //standard cuts
  } else if (trainConfig == 100){ // EMCAL clusters pp 8 TeV
    cuts.AddCutCalo("00010113","1111111067032220000","01631031000000d0"); // std
    cuts.AddCutCalo("00052113","1111111067032220000","01631031000000d0"); // std
    cuts.AddCutCalo("00081113","1111111067032220000","01631031000000d0"); // std
  } else if (trainConfig == 101){ // EMCAL clusters pp 8 TeV
    cuts.AddCutCalo("00010113","1111111067032220000","01631031000000d0"); // std
    cuts.AddCutCalo("00052113","1111111067032220000","01631031000000d0"); // std
    cuts.AddCutCalo("00081113","1111111067032220000","01631031000000d0"); // std

    // ATTENTION: adapted for 8 TeV dirGamma - ADJUSTED M02 -> l
    // 8 TeV variations
  } else if (trainConfig == 102){ // EMCAL clusters pp 8 TeV, timing+minEnergy variation
    cuts.AddCutCalo("00010113","11111110570322l0000","01631031000000d0"); // time -50ns_50ns
    cuts.AddCutCalo("00010113","11111110770322l0000","01631031000000d0"); // time -30ns_30ns
    cuts.AddCutCalo("00010113","11111110870322l0000","01631031000000d0"); // time -20ns_30ns
    cuts.AddCutCalo("00010113","11111110670122l0000","01631031000000d0"); //0.5 GeV/c
    cuts.AddCutCalo("00010113","11111110670222l0000","01631031000000d0"); //0.6 GeV/c
    cuts.AddCutCalo("00010113","11111110670422l0000","01631031000000d0"); //0.8 GeV/c
    cuts.AddCutCalo("00010113","11111110670522l0000","01631031000000d0"); //0.9 GeV/c
  } else if (trainConfig == 103){ //EMCAL M02 variation
    cuts.AddCutCalo("00010113","1111111067032200000","01631031000000d0"); //no max M02 cut
    cuts.AddCutCalo("00010113","11111110670322l0000","01631031000000d0"); // M02 pt dep with new std: 0.32, 0.0072, 0.5
    cuts.AddCutCalo("00010113","11111110670322d0000","01631031000000d0"); // M02, pt dep with  0.27, 0.0072, 0.4
    cuts.AddCutCalo("00010113","11111110670322e0000","01631031000000d0"); // M02, pt dep with  0.31, 0.0072, 0.5
    cuts.AddCutCalo("00010113","11111110670322f0000","01631031000000d0"); // M02, pt dep with  0.36, 0.0072, 0.7
    cuts.AddCutCalo("00010113","11111110670322m0000","01631031000000d0"); // M02, pt dep with  0.32, 0.0152, 0.5
    cuts.AddCutCalo("00010113","11111110670322g0000","01631031000000d0"); // M02, pt dep with  0.37, 0.0072, 0.7
    cuts.AddCutCalo("00010113","11111110670322h0000","01631031000000d0"); // M02, pT-dep with  0.30, 0.0072, 0.5
    cuts.AddCutCalo("00010113","11111110670322i0000","01631031000000d0"); // M02, pT-dep with  0.35, 0.0072, 0.7
    cuts.AddCutCalo("00010113","11111110670322j0000","01631031000000d0"); // M02, pT-dep with  0.25, 0.0072, 0.39
    cuts.AddCutCalo("00010113","11111110670322r0000","01631031000000d0"); // M02, pT-dep with  0.25, 0.0072, 0.5
    cuts.AddCutCalo("00010113","11111110670322s0000","01631031000000d0"); // M02, pT-dep with  0.32, 0.0238, 0.7
    cuts.AddCutCalo("00010113","11111110670322n0000","01631031000000d0"); // M02, pT-dep with  0.32, 0.0238, 0.7
  } else if (trainConfig == 104){ //EMCAL minNCells,with/without TRD variation
    cuts.AddCutCalo("00010113","11111110670312l0000","01631031000000d0"); //n cells >= 1
    cuts.AddCutCalo("00010113","11111110670332l0000","01631031000000d0"); //n cells >= 3
    cuts.AddCutCalo("00010113","11121110670322l0000","01631031000000d0"); //only modules with TRD infront
    cuts.AddCutCalo("00010113","11113110670322l0000","01631031000000d0"); //no modules with TRD infront
  } else if (trainConfig == 105){  // trackMatching variations
    cuts.AddCutCalo("00010113","11111110660322l0000","01631031000000d0"); //
    cuts.AddCutCalo("00010113","11111110680322l0000","01631031000000d0"); //
    cuts.AddCutCalo("00010113","11111110690322l0000","01631031000000d0"); //
    cuts.AddCutCalo("00010113","11111110630322l0000","01631031000000d0"); //
    cuts.AddCutCalo("00010113","11111110600322l0000","01631031000000d0"); //
  } else if (trainConfig == 106){ // EMCAL clusters pp 8 TeV, combining cluster within time window and without
    cuts.AddCutCalo("00010113","11111110070322l0000","01631031000000d0"); //
  } else if (trainConfig == 107){ // EMCAL clusters open angle variation
    cuts.AddCutCalo("00010113","11111110670322l0000","01631031000000a0"); // min open angle - 0. + 1 cell diag
    cuts.AddCutCalo("00010113","11111110670322l0000","01631031000000d0"); // min open angle - 0.017 + 1 cell diag
    cuts.AddCutCalo("00010113","11111110670322l0000","01631031000000b0"); // min open angle - 0.0152 + 1 cell diag
    cuts.AddCutCalo("00010113","11111110670322l0000","01631031000000c0"); // min open angle - 0.016 + 1 cell diag
    cuts.AddCutCalo("00010113","11111110670322l0000","01631031000000e0"); // min open angle - 0.018 + 1 cell diag
    cuts.AddCutCalo("00010113","11111110670322l0000","01631031000000f0"); // min open angle - 0.019 + 1 cell diag
  } else if (trainConfig == 108){ // EMCAL clusters pp 8 TeV, Different DistanceToBadChannels
    cuts.AddCutCalo("00010113","11111111670322l0000","01631031000000d0"); //
    cuts.AddCutCalo("00010113","11111112670322l0000","01631031000000d0"); //
    cuts.AddCutCalo("00010113","11111113670322l0000","01631031000000d0"); //
    cuts.AddCutCalo("00010113","11111115670322l0000","01631031000000d0"); //
    cuts.AddCutCalo("00010113","11111116670322l0000","01631031000000d0"); //
  } else if (trainConfig == 109){ // EMCAL clusters pp 8 TeV, Different NonLinearities
    cuts.AddCutCalo("00010113","11111010670322l0000","01631031000000d0"); // NonLinearity kSDMv5
    cuts.AddCutCalo("00010113","11111130670322l0000","01631031000000d0"); // NonLinearity kTestBeamv2 + LHC12 ConvCalo
    cuts.AddCutCalo("00010113","11111140670322l0000","01631031000000d0"); // NonLinearity kTestBeamv2 + LHC12 Calo
  } else if (trainConfig == 110){ // EMCAL clusters pp 8 TeV, Different NonLinearities
    cuts.AddCutCalo("00010113","11111110670322l0000","01631031000000d0"); // NonLinearity LHC12 ConvCalo
    cuts.AddCutCalo("00010113","11111120670322l0000","01631031000000d0"); // NonLinearity LHC12 Calo
    cuts.AddCutCalo("00010113","11111210670322l0000","01631031000000d0"); // NonLinearity LHC12 ConvCalo MassRatioFits
    cuts.AddCutCalo("00010113","11111220670322l0000","01631031000000d0"); // NonLinearity LHC12 Calo MassRatioFits
    cuts.AddCutCalo("00010113","11111000670322l0000","01631031000000d0"); // NonLinearity none
  } else if (trainConfig == 111){  // EMCAL clusters, different triggers no NonLinearity
    cuts.AddCutCalo("00010113","11111000670322l0000","01631031000000d0");
    cuts.AddCutCalo("00052113","11111000670322l0000","01631031000000d0"); // EMC7
    cuts.AddCutCalo("00081113","11111000670322l0000","01631031000000d0"); // EMCEG1,
  } else if (trainConfig == 112){ // EMCAL clusters, exotic cut var
    cuts.AddCutCalo("00010113","11111110672322l0000","01631031000000d0"); //
    cuts.AddCutCalo("00010113","11111110673322l0000","01631031000000d0"); //
    cuts.AddCutCalo("00010113","11111110675322l0000","01631031000000d0"); //
    cuts.AddCutCalo("00010113","11111110677322l0000","01631031000000d0"); //
    cuts.AddCutCalo("00010113","11111110679322l0000","01631031000000d0"); //
  } else if (trainConfig == 113){  // trackMatching variations
    cuts.AddCutCalo("00010113","11111110620322l0000","01631031000000d0"); //
    cuts.AddCutCalo("00010113","11111110630322l0000","01631031000000d0"); //
    cuts.AddCutCalo("00010113","11111110640322l0000","01631031000000d0"); //
    cuts.AddCutCalo("00010113","11111110650322l0000","01631031000000d0"); //
  } else if (trainConfig == 114){ // EMCAL clusters pp 8 TeV
    cuts.AddCutCalo("00010013","11111110670322l0000","01631031000000d0"); // std - no pileup cut
  } else if (trainConfig == 115){ // EMCAL clusters pp 8 TeV
    cuts.AddCutCalo("00010113","11111110670322l0000","01631051000000d0"); // alpha
    cuts.AddCutCalo("00010113","11111110670322l0000","01631061000000d0"); //
  } else if (trainConfig == 116){ // EMCAL clusters pp 8 TeV
    cuts.AddCutCalo("00010113","11111110670322l0000","01633031000000d0"); // rapidity
    cuts.AddCutCalo("00010113","11111110670322l0000","01634031000000d0"); //


  } else if (trainConfig == 118){ // EMCAL clusters pp 8 TeV - no SPD PileUp
    cuts.AddCutCalo("00010113","11111110670322l0000","01631031000000d0"); // std
    cuts.AddCutCalo("00010013","11111110670322l0000","01631031000000d0"); // std - no pileup cut
    cuts.AddCutCalo("00052113","11111110670322l0000","01631031000000d0"); // std
    cuts.AddCutCalo("00052013","11111110670322l0000","01631031000000d0"); // std - no pileup cut
    cuts.AddCutCalo("00081113","11111110670322l0000","01631031000000d0"); // std
    cuts.AddCutCalo("00081013","11111110670322l0000","01631031000000d0"); // std - no pileup cut

  // only std cuts
  } else if (trainConfig == 119){ // EMCAL clusters pp 8 TeV
    cuts.AddCutCalo("00010113","1111111067032220000","01631031000000d0"); // std
  } else if (trainConfig == 120){ // EMCAL clusters pp 8 TeV
    cuts.AddCutCalo("00010113","11111110670322l0000","01631031000000d0"); // std - dirGamma

    // ATTENTION: adapted for 8 TeV dirGamma  - ADJUSTED M02 Cut 2 -> p
    //8 TeV kEMC7 variations
  } else if (trainConfig == 121){ // EMCAL clusters pp 8 TeV, timing+minEnergy variation
    cuts.AddCutCalo("00052113","11111110570322l0000","01631031000000d0"); // time -50ns_50ns
    cuts.AddCutCalo("00052113","11111110770322l0000","01631031000000d0"); // time -30ns_30ns
    cuts.AddCutCalo("00052113","11111110870322l0000","01631031000000d0"); // time -20ns_30ns
    cuts.AddCutCalo("00052113","11111110670122l0000","01631031000000d0"); //0.5 GeV/c
    cuts.AddCutCalo("00052113","11111110670222l0000","01631031000000d0"); //0.6 GeV/c
    cuts.AddCutCalo("00052113","11111110670422l0000","01631031000000d0"); //0.8 GeV/c
    cuts.AddCutCalo("00052113","11111110670522l0000","01631031000000d0"); //0.9 GeV/c
  } else if (trainConfig == 122){ //EMCAL M02 variation
    cuts.AddCutCalo("00052113","1111111067032200000","01631031000000d0"); //no max M02 cut
    cuts.AddCutCalo("00052113","11111110670322l0000","01631031000000d0"); // M02 pt dep with new std: 0.32, 0.0072, 0.5
    cuts.AddCutCalo("00052113","11111110670322d0000","01631031000000d0"); // M02, pt dep with  0.27, 0.0072, 0.4
    cuts.AddCutCalo("00052113","11111110670322e0000","01631031000000d0"); // M02, pt dep with  0.31, 0.0072, 0.5
    cuts.AddCutCalo("00052113","11111110670322f0000","01631031000000d0"); // M02, pt dep with  0.36, 0.0072, 0.7
    cuts.AddCutCalo("00052113","11111110670322m0000","01631031000000d0"); // M02, pt dep with  0.32, 0.0152, 0.5
    cuts.AddCutCalo("00052113","11111110670322g0000","01631031000000d0"); // M02, pt dep with  0.37, 0.0072, 0.7
    cuts.AddCutCalo("00052113","11111110670322h0000","01631031000000d0"); // M02, pT-dep with  0.30, 0.0072, 0.5
    cuts.AddCutCalo("00052113","11111110670322i0000","01631031000000d0"); // M02, pT-dep with  0.35, 0.0072, 0.7
    cuts.AddCutCalo("00052113","11111110670322j0000","01631031000000d0"); // M02, pT-dep with  0.25, 0.0072, 0.39
    cuts.AddCutCalo("00052113","11111110670322r0000","01631031000000d0"); // M02, pT-dep with  0.25, 0.0072, 0.5
    cuts.AddCutCalo("00052113","11111110670322s0000","01631031000000d0"); // M02, pT-dep with  0.32, 0.0238, 0.7
    cuts.AddCutCalo("00052113","11111110670322n0000","01631031000000d0"); // M02, pT-dep with  0.32, 0.0238, 0.7
  } else if (trainConfig == 123){ //EMCAL minNCells, with/without TRD variation
    cuts.AddCutCalo("00052113","11111110670312l0000","01631031000000d0"); //n cells >= 1
    cuts.AddCutCalo("00052113","11111110670332l0000","01631031000000d0"); //n cells >= 3
    cuts.AddCutCalo("00052113","11121110670322l0000","01631031000000d0"); //only modules with TRD infront
    cuts.AddCutCalo("00052113","11113110670322l0000","01631031000000d0"); //no modules with TRD infront
  } else if (trainConfig == 124){  // trackMatching variations
    cuts.AddCutCalo("00052113","11111110660322l0000","01631031000000d0"); //
    cuts.AddCutCalo("00052113","11111110670322l0000","01631031000000d0"); //
    cuts.AddCutCalo("00052113","11111110680322l0000","01631031000000d0"); //
    cuts.AddCutCalo("00052113","11111110690322l0000","01631031000000d0"); //
    cuts.AddCutCalo("00052113","11111110630322l0000","01631031000000d0"); //
    cuts.AddCutCalo("00052113","11111110600322l0000","01631031000000d0"); //
  } else if (trainConfig == 125){ // EMCAL clusters pp 8 TeV, combining cluster within time window and without
    cuts.AddCutCalo("00052113","11111110070322l0000","01631031000000d0"); //
  } else if (trainConfig == 126){ // EMCAL clusters open angle variation
    cuts.AddCutCalo("00052113","11111110670322l0000","01631031000000a0"); // min open angle - 0. + 1 cell diag
    cuts.AddCutCalo("00052113","11111110670322l0000","01631031000000d0"); // min open angle - 0.017 + 1 cell diag
    cuts.AddCutCalo("00052113","11111110670322l0000","01631031000000b0"); // min open angle - 0.0152 + 1 cell diag
    cuts.AddCutCalo("00052113","11111110670322l0000","01631031000000c0"); // min open angle - 0.016 + 1 cell diag
    cuts.AddCutCalo("00052113","11111110670322l0000","01631031000000e0"); // min open angle - 0.018 + 1 cell diag
    cuts.AddCutCalo("00052113","11111110670322l0000","01631031000000f0"); // min open angle - 0.019 + 1 cell diag
  } else if (trainConfig == 127){ // EMCAL clusters pp 8 TeV, Different DistanceToBadChannels
    cuts.AddCutCalo("00052113","11111110670322l0000","01631031000000d0"); //
    cuts.AddCutCalo("00052113","11111111670322l0000","01631031000000d0"); //
    cuts.AddCutCalo("00052113","11111112670322l0000","01631031000000d0"); //
    cuts.AddCutCalo("00052113","11111113670322l0000","01631031000000d0"); //
    cuts.AddCutCalo("00052113","11111115670322l0000","01631031000000d0"); //
    cuts.AddCutCalo("00052113","11111116670322l0000","01631031000000d0"); //
  } else if (trainConfig == 128){ // EMCAL clusters pp 8 TeV, Different NonLinearities
    cuts.AddCutCalo("00052113","11111010670322l0000","01631031000000d0"); // NonLinearity kSDMv5
    cuts.AddCutCalo("00052113","11111130670322l0000","01631031000000d0"); // NonLinearity kTestBeamv2 + LHC12 ConvCalo
    cuts.AddCutCalo("00052113","11111140670322l0000","01631031000000d0"); // NonLinearity kTestBeamv2 + LHC12 Calo
  } else if (trainConfig == 129){ // EMCAL clusters pp 8 TeV, Different NonLinearities
    cuts.AddCutCalo("00052113","11111110670322l0000","01631031000000d0"); // NonLinearity LHC12 ConvCalo
    cuts.AddCutCalo("00052113","11111120670322l0000","01631031000000d0"); // NonLinearity LHC12 Calo
    cuts.AddCutCalo("00052113","11111210670322l0000","01631031000000d0"); // NonLinearity LHC12 ConvCalo MassRatioFits
    cuts.AddCutCalo("00052113","11111220670322l0000","01631031000000d0"); // NonLinearity LHC12 Calo MassRatioFits
    cuts.AddCutCalo("00052113","11111000670322l0000","01631031000000d0"); // NonLinearity none
  } else if (trainConfig == 130){ // EMCAL clusters, exotic cut var
    cuts.AddCutCalo("00052113","11111110670322l0000","01631031000000d0"); //
    cuts.AddCutCalo("00052113","11111110672322l0000","01631031000000d0"); //
    cuts.AddCutCalo("00052113","11111110673322l0000","01631031000000d0"); //
    cuts.AddCutCalo("00052113","11111110675322l0000","01631031000000d0"); //
    cuts.AddCutCalo("00052113","11111110677322l0000","01631031000000d0"); //
    cuts.AddCutCalo("00052113","11111110679322l0000","01631031000000d0"); //
  } else if (trainConfig == 131){  // trackMatching variations
    cuts.AddCutCalo("00052113","11111110670322l0000","01631031000000d0"); //
    cuts.AddCutCalo("00052113","11111110620322l0000","01631031000000d0"); //
    cuts.AddCutCalo("00052113","11111110630322l0000","01631031000000d0"); //
    cuts.AddCutCalo("00052113","11111110640322l0000","01631031000000d0"); //
  } else if (trainConfig == 132){ // EMCAL clusters pp 8 TeV
    cuts.AddCutCalo("00052013","11111110670322l0000","01631031000000d0"); // std - no pileup cut
  } else if (trainConfig == 133){ // EMCAL clusters pp 8 TeV
    cuts.AddCutCalo("00052113","11111110670322l0000","01631051000000d0"); //
    cuts.AddCutCalo("00052113","11111110670322l0000","01631061000000d0"); //
  } else if (trainConfig == 134){ // EMCAL clusters pp 8 TeV
    cuts.AddCutCalo("00052113","11111110670322l0000","01633031000000d0"); //
    cuts.AddCutCalo("00052113","11111110670322l0000","01634031000000d0"); //

  // only std cuts
  } else if (trainConfig == 139){ // EMCAL clusters pp 8 TeV
    cuts.AddCutCalo("00052113","1111111067032220000","01631031000000d0"); // std
  } else if (trainConfig == 140){ // EMCAL clusters pp 8 TeV
    cuts.AddCutCalo("00052113","11111110670322l0000","01631031000000d0"); // std - dirGAMMA

    // ATTENTION: adapted for 8 TeV dirGamma - ADJUSTED M02 Cut 2 -> p
    //8 TeV kEMCEGA variations
  } else if (trainConfig == 141){ // EMCAL clusters pp 8 TeV, timing+minEnergy variation
    cuts.AddCutCalo("00081113","11111110570322l0000","01631031000000d0"); // time -50ns_50ns
    cuts.AddCutCalo("00081113","11111110770322l0000","01631031000000d0"); // time -30ns_30ns
    cuts.AddCutCalo("00081113","11111110870322l0000","01631031000000d0"); // time -20ns_30ns
    cuts.AddCutCalo("00081113","11111110670122l0000","01631031000000d0"); //0.5 GeV/c
    cuts.AddCutCalo("00081113","11111110670222l0000","01631031000000d0"); //0.6 GeV/c
    cuts.AddCutCalo("00081113","11111110670422l0000","01631031000000d0"); //0.8 GeV/c
    cuts.AddCutCalo("00081113","11111110670522l0000","01631031000000d0"); //0.9 GeV/c
  } else if (trainConfig == 142){ //EMCAL M02 variation
    cuts.AddCutCalo("00081113","1111111067032200000","01631031000000d0"); //no max M02 cut
    cuts.AddCutCalo("00081113","11111110670322l0000","01631031000000d0"); //M02, pT-dep
    cuts.AddCutCalo("00081113","11111110670322d0000","01631031000000d0"); //M02, pT-dep
    cuts.AddCutCalo("00081113","11111110670322e0000","01631031000000d0"); //M02, pT-dep
    cuts.AddCutCalo("00081113","11111110670322f0000","01631031000000d0"); //M02, pT-dep
    cuts.AddCutCalo("00081113","11111110670322m0000","01631031000000d0"); //M02, pT-dep
    cuts.AddCutCalo("00081113","11111110670322g0000","01631031000000d0"); //M02, pT-dep
    cuts.AddCutCalo("00081113","11111110670322h0000","01631031000000d0"); //M02, pT-dep
    cuts.AddCutCalo("00081113","11111110670322i0000","01631031000000d0"); //M02, pT-dep
    cuts.AddCutCalo("00081113","11111110670322j0000","01631031000000d0"); //M02, pT-dep
    cuts.AddCutCalo("00081113","11111110670322r0000","01631031000000d0"); //M02, pT-dep
    cuts.AddCutCalo("00081113","11111110670322s0000","01631031000000d0"); //M02, pT-dep
    cuts.AddCutCalo("00081113","11111110670322n0000","01631031000000d0"); //M02, pT-dep
  } else if (trainConfig == 143){ //EMCAL minNCells, with/without TRD variation
    cuts.AddCutCalo("00081113","11111110670312l0000","01631031000000d0"); //n cells >= 1
    cuts.AddCutCalo("00081113","11111110670332l0000","01631031000000d0"); //n cells >= 3
    cuts.AddCutCalo("00081113","11121110670322l0000","01631031000000d0"); //only modules with TRD infront
    cuts.AddCutCalo("00081113","11113110670322l0000","01631031000000d0"); //no modules with TRD infront
  } else if (trainConfig == 144){  // trackMatching variations
    cuts.AddCutCalo("00081113","11111110660322l0000","01631031000000d0"); //
    cuts.AddCutCalo("00081113","11111110670322l0000","01631031000000d0"); //
    cuts.AddCutCalo("00081113","11111110680322l0000","01631031000000d0"); //
    cuts.AddCutCalo("00081113","11111110690322l0000","01631031000000d0"); //
    cuts.AddCutCalo("00081113","11111110630322l0000","01631031000000d0"); //
    cuts.AddCutCalo("00081113","11111110600322l0000","01631031000000d0"); //
  } else if (trainConfig == 145){ // EMCAL clusters pp 8 TeV, combining cluster within time window and without
    cuts.AddCutCalo("00081113","11111110070322l0000","01631031000000d0"); //
  } else if (trainConfig == 146){ // EMCAL clusters open angle variation
    cuts.AddCutCalo("00081113","11111110670322l0000","01631031000000a0"); // min open angle - 0. + 1 cell diag
    cuts.AddCutCalo("00081113","11111110670322l0000","01631031000000d0"); // min open angle - 0.017 + 1 cell diag
    cuts.AddCutCalo("00081113","11111110670322l0000","01631031000000b0"); // min open angle - 0.0152 + 1 cell diag
    cuts.AddCutCalo("00081113","11111110670322l0000","01631031000000c0"); // min open angle - 0.016 + 1 cell diag
    cuts.AddCutCalo("00081113","11111110670322l0000","01631031000000e0"); // min open angle - 0.018 + 1 cell diag
    cuts.AddCutCalo("00081113","11111110670322l0000","01631031000000f0"); // min open angle - 0.019 + 1 cell diag
  } else if (trainConfig == 147){ // EMCAL clusters pp 8 TeV, Different DistanceToBadChannels
    cuts.AddCutCalo("00081113","11111110670322l0000","01631031000000d0"); //
    cuts.AddCutCalo("00081113","11111111670322l0000","01631031000000d0"); //
    cuts.AddCutCalo("00081113","11111112670322l0000","01631031000000d0"); //
    cuts.AddCutCalo("00081113","11111113670322l0000","01631031000000d0"); //
    cuts.AddCutCalo("00081113","11111115670322l0000","01631031000000d0"); //
    cuts.AddCutCalo("00081113","11111116670322l0000","01631031000000d0"); //
  } else if (trainConfig == 148){ // EMCAL clusters pp 8 TeV, Different NonLinearities
    cuts.AddCutCalo("00081113","11111010670322l0000","01631031000000d0"); // NonLinearity kSDMv5
    cuts.AddCutCalo("00081113","11111130670322l0000","01631031000000d0"); // NonLinearity kTestBeamv2 + LHC12 ConvCalo
    cuts.AddCutCalo("00081113","11111140670322l0000","01631031000000d0"); // NonLinearity kTestBeamv2 + LHC12 Calo
  } else if (trainConfig == 149){ // EMCAL clusters pp 8 TeV, Different NonLinearities
    cuts.AddCutCalo("00081113","11111110670322l0000","01631031000000d0"); // NonLinearity LHC12 ConvCalo
    cuts.AddCutCalo("00081113","11111120670322l0000","01631031000000d0"); // NonLinearity LHC12 Calo
    cuts.AddCutCalo("00081113","11111210670322l0000","01631031000000d0"); // NonLinearity LHC12 ConvCalo MassRatioFits
    cuts.AddCutCalo("00081113","11111220670322l0000","01631031000000d0"); // NonLinearity LHC12 Calo MassRatioFits
    cuts.AddCutCalo("00081113","11111000670322l0000","01631031000000d0"); // NonLinearity none
  } else if (trainConfig == 150){ // EMCAL clusters, exotic cut var
    cuts.AddCutCalo("00081113","11111110670322l0000","01631031000000d0"); //
    cuts.AddCutCalo("00081113","11111110672322l0000","01631031000000d0"); //
    cuts.AddCutCalo("00081113","11111110673322l0000","01631031000000d0"); //
    cuts.AddCutCalo("00081113","11111110675322l0000","01631031000000d0"); //
    cuts.AddCutCalo("00081113","11111110677322l0000","01631031000000d0"); //
    cuts.AddCutCalo("00081113","11111110679322l0000","01631031000000d0"); //
  } else if (trainConfig == 151){  // trackMatching variations
    cuts.AddCutCalo("00081113","11111110670322l0000","01631031000000d0"); //
    cuts.AddCutCalo("00081113","11111110620322l0000","01631031000000d0"); //
    cuts.AddCutCalo("00081113","11111110630322l0000","01631031000000d0"); //
    cuts.AddCutCalo("00081113","11111110640322l0000","01631031000000d0"); //
  } else if (trainConfig == 152){ // EMCAL clusters pp 8 TeV
    cuts.AddCutCalo("00081013","11111110670322l0000","01631031000000d0"); // std - no pileup cut
  } else if (trainConfig == 153){ // EMCAL clusters pp 8 TeV
    cuts.AddCutCalo("00081113","11111110670322l0000","01631051000000d0"); //
    cuts.AddCutCalo("00081113","11111110670322l0000","01631061000000d0"); //
  } else if (trainConfig == 154){ // EMCAL clusters pp 8 TeV
    cuts.AddCutCalo("00081113","11111110670322l0000","01633031000000d0"); //
    cuts.AddCutCalo("00081113","11111110670322l0000","01634031000000d0"); //

  // only std cuts
  } else if (trainConfig == 159){ // EMCAL clusters pp 8 TeV
    cuts.AddCutCalo("00081113","1111111067032220000","01631031000000d0"); // std
  } else if (trainConfig == 160){ // EMCAL clusters pp 8 TeV
    cuts.AddCutCalo("00081113","11111110670322l0000","01631031000000d0"); // std - dirGAMMA

  // CutStudies for DirGamma
  } else if (trainConfig == 161){ //
    cuts.AddCutCalo("00010113","11111110670322l0000","01631031000000d0"); // std
    cuts.AddCutCalo("00052113","11111110670322l0000","01631031000000d0"); // std
    cuts.AddCutCalo("00081113","11111110670322l0000","01631031000000d0"); // std

  } else if (trainConfig == 162){ //
    cuts.AddCutCalo("00010113","11111110670322l0000","01631031000000d0"); // std
  } else if (trainConfig == 163){ //
    cuts.AddCutCalo("00052113","11111110670322l0000","01631031000000d0"); // std
  } else if (trainConfig == 164){ //
    cuts.AddCutCalo("00081113","11111110670322l0000","01631031000000d0"); // std

  } else if (trainConfig == 165){
    cuts.AddCutCalo("00010113","11111110670322l0000","01631031000000d0"); // M02 pt dep with new std: 0.32, 0.0072, 0.5
    cuts.AddCutCalo("00010113","11111110670322d0000","01631031000000d0"); // M02, pt dep with  0.27, 0.0072, 0.4
    cuts.AddCutCalo("00010113","11111110670322e0000","01631031000000d0"); // M02, pt dep with  0.31, 0.0072, 0.5
    cuts.AddCutCalo("00010113","11111110670322f0000","01631031000000d0"); // M02, pt dep with  0.36, 0.0072, 0.7
  } else if (trainConfig == 166){
    cuts.AddCutCalo("00010113","11111110670322m0000","01631031000000d0"); // M02, pt dep with  0.32, 0.0152, 0.5
    cuts.AddCutCalo("00010113","11111110670322g0000","01631031000000d0"); // M02, pt dep with  0.37, 0.0072, 0.7
    cuts.AddCutCalo("00010113","11111110670322h0000","01631031000000d0"); // M02, pT-dep with  0.30, 0.0072, 0.5
    cuts.AddCutCalo("00010113","11111110670322i0000","01631031000000d0"); // M02, pT-dep with  0.35, 0.0072, 0.7
  } else if (trainConfig == 167){
    cuts.AddCutCalo("00010113","11111110670322j0000","01631031000000d0"); // M02, pT-dep with  0.25, 0.0072, 0.39
    cuts.AddCutCalo("00010113","11111110670322r0000","01631031000000d0"); // M02, pT-dep with  0.25, 0.0072, 0.5
    cuts.AddCutCalo("00010113","11111110670322s0000","01631031000000d0"); // M02, pT-dep with  0.32, 0.0238, 0.7
    cuts.AddCutCalo("00010113","11111110670322n0000","01631031000000d0"); // M02, pT-dep with  0.32, 0.0238, 0.7

  } else if (trainConfig == 170){ // EMCAL clusters pp 8 TeV 100MeV aggregation
    cuts.AddCutCalo("00010113","111110106f032230000","01631031000000d0"); // std
    cuts.AddCutCalo("00052113","111110106f032230000","01631031000000d0"); // std
    cuts.AddCutCalo("00081113","111110106f032230000","01631031000000d0"); // std
  } else if (trainConfig == 171){ // EMCAL clusters pp 8 TeV 50MeV aggregation
    cuts.AddCutCalo("00010113","111110206f032230000","01631031000000d0"); // std
    cuts.AddCutCalo("00052113","111110206f032230000","01631031000000d0"); // std
    cuts.AddCutCalo("00081113","111110206f032230000","01631031000000d0"); // std
  } else if (trainConfig == 172){ // EMCAL clusters pp 8 TeV 150MeV aggregation
    cuts.AddCutCalo("00010113","111110306f032230000","01631031000000d0"); // std
    cuts.AddCutCalo("00052113","111110306f032230000","01631031000000d0"); // std
    cuts.AddCutCalo("00081113","111110306f032230000","01631031000000d0"); // std
  } else if (trainConfig == 173){ // EMCAL clusters pp 8 TeV 300MeV aggregation
    cuts.AddCutCalo("00010113","111110406f032230000","01631031000000d0"); // std
    cuts.AddCutCalo("00052113","111110406f032230000","01631031000000d0"); // std
    cuts.AddCutCalo("00081113","111110406f032230000","01631031000000d0"); // std

  } else if (trainConfig == 174){ // EMCAL clusters pp 8 TeV, TB+finetuning CCRF
    cuts.AddCutCalo("00010113","111113806f032230000","01631031000000d0"); // std
    cuts.AddCutCalo("00052113","111113806f032230000","01631031000000d0"); // std
    cuts.AddCutCalo("00081113","111113806f032230000","01631031000000d0"); // std
  } else if (trainConfig == 175){ // EMCAL clusters pp 8 TeV, TB+finetuning CCRF
    cuts.AddCutCalo("00010113","111113906f032230000","01631031000000d0"); // std
    cuts.AddCutCalo("00052113","111113906f032230000","01631031000000d0"); // std
    cuts.AddCutCalo("00081113","111113906f032230000","01631031000000d0"); // std
  } else if (trainConfig == 176){ // EMCAL clusters pp 8 TeV, TB+finetuning CCRF
    cuts.AddCutCalo("00010113","111113106f032230000","01631031000000d0"); // std
    cuts.AddCutCalo("00052113","111113106f032230000","01631031000000d0"); // std
    cuts.AddCutCalo("00081113","111113106f032230000","01631031000000d0"); // std
  } else if (trainConfig == 177){ // EMCAL clusters pp 8 TeV, TB+finetuning CRF
    cuts.AddCutCalo("00010113","111113206f032230000","01631031000000d0"); // std
    cuts.AddCutCalo("00052113","111113206f032230000","01631031000000d0"); // std
    cuts.AddCutCalo("00081113","111113206f032230000","01631031000000d0"); // std
  } else if (trainConfig == 178){ // EMCAL clusters pp 8 TeV, TB+finetuning CCMF
    cuts.AddCutCalo("00010113","111113306f032230000","01631031000000d0"); // std
    cuts.AddCutCalo("00052113","111113306f032230000","01631031000000d0"); // std
    cuts.AddCutCalo("00081113","111113306f032230000","01631031000000d0"); // std
  } else if (trainConfig == 179){ // EMCAL clusters pp 8 TeV, TB+finetuning CMF
    cuts.AddCutCalo("00010113","111113406f032230000","01631031000000d0"); // std
    cuts.AddCutCalo("00052113","111113406f032230000","01631031000000d0"); // std
    cuts.AddCutCalo("00081113","111113406f032230000","01631031000000d0"); // std

  } else if (trainConfig == 180){ // EMCAL clusters pp 8 TeV, TB+finetuning CCRF
    cuts.AddCutCalo("00010113","111116106f032230000","01631031000000d0"); // std
    cuts.AddCutCalo("00052113","111116106f032230000","01631031000000d0"); // std
    cuts.AddCutCalo("00081113","111116106f032230000","01631031000000d0"); // std
  } else if (trainConfig == 181){ // EMCAL clusters pp 8 TeV, TB+finetuning CRF
    cuts.AddCutCalo("00010113","111116206f032230000","01631031000000d0"); // std
    cuts.AddCutCalo("00052113","111116206f032230000","01631031000000d0"); // std
    cuts.AddCutCalo("00081113","111116206f032230000","01631031000000d0"); // std
  } else if (trainConfig == 182){ // EMCAL clusters pp 8 TeV, TB+finetuning CCMF
    cuts.AddCutCalo("00010113","111116306f032230000","01631031000000d0"); // std
    cuts.AddCutCalo("00052113","111116306f032230000","01631031000000d0"); // std
    cuts.AddCutCalo("00081113","111116306f032230000","01631031000000d0"); // std
  } else if (trainConfig == 183){ // EMCAL clusters pp 8 TeV, TB+finetuning CMF
    cuts.AddCutCalo("00010113","111116406f032230000","01631031000000d0"); // std
    cuts.AddCutCalo("00052113","111116406f032230000","01631031000000d0"); // std
    cuts.AddCutCalo("00081113","111116406f032230000","01631031000000d0"); // std
  } else if (trainConfig == 184){ // special MC fit
    cuts.AddCutCalo("00010113","111116906f032230000","01631031000000d0"); // std
    cuts.AddCutCalo("00052113","111116906f032230000","01631031000000d0"); // std
    cuts.AddCutCalo("00081113","111116906f032230000","01631031000000d0"); // std

  } else if (trainConfig == 185){ // no NCell cut
    cuts.AddCutCalo("00010113","111113106f030230000","01631031000000d0"); // std
    cuts.AddCutCalo("00052113","111113106f030230000","01631031000000d0"); // std
    cuts.AddCutCalo("00081113","111113106f030230000","01631031000000d0"); // std

  } else if (trainConfig == 186){ // EMCAL clusters pp 8 TeV, TB+finetuning CRF
    cuts.AddCutCalo("00010113","1111132060032230000","01631031000000d0"); // std
    cuts.AddCutCalo("00052113","1111132060032230000","01631031000000d0"); // std
    cuts.AddCutCalo("00081113","1111132060032230000","01631031000000d0"); // std

  //multiple std dirGAMMA cuts for different studies
  } else if (trainConfig == 190){ // EMCAL clusters pp 8 TeV
    cuts.AddCutCalo("00010113","11111110670322l0000","01631031000000d0"); // std
  } else if (trainConfig == 191){ // EMCAL clusters pp 8 TeV
    cuts.AddCutCalo("00010113","11111110670322l0000","01631031000000d0"); // std
  } else if (trainConfig == 192){ // EMCAL clusters pp 8 TeV
    cuts.AddCutCalo("00010113","11111110670322l0000","01631031000000d0"); // std

  // 7 TeV
  } else if (trainConfig == 200){ // EMCAL clusters pp 7 TeV, pT dep matching
    cuts.AddCutCalo("00000113","11111110b7032220000","01631031000000d0"); // std NL
  } else if (trainConfig == 201){ // EMCAL clusters pp 7 TeV, pT dep matching
    cuts.AddCutCalo("00000113","11111310b7032220000","01631031000000d0"); // std TB NL

    // ATTENTION: adapted for dirGamma - ADJUSTED M02 -> l
  } else if (trainConfig == 202){ // EMCAL clusters pp 7 TeV, timing+minEnergy variation
    cuts.AddCutCalo("00000113","11111110370322l0000","01631031000000d0"); // time
    cuts.AddCutCalo("00000113","11111110470322l0000","01631031000000d0"); // time
    cuts.AddCutCalo("00000113","11111110570322l0000","01631031000000d0"); // time
    cuts.AddCutCalo("00000113","11111110b70122l0000","01631031000000d0"); //0.5 GeV/c
    cuts.AddCutCalo("00000113","11111110b70222l0000","01631031000000d0"); //0.6 GeV/c
    cuts.AddCutCalo("00000113","11111110b70422l0000","01631031000000d0"); //0.8 GeV/c
    cuts.AddCutCalo("00000113","11111110b70522l0000","01631031000000d0"); //0.9 GeV/c
  } else if (trainConfig == 203){ //EMCAL M02 variation
    cuts.AddCutCalo("00000113","11111110b7032200000","01631031000000d0"); //no max M02 cut
    cuts.AddCutCalo("00000113","11111110b70322l0000","01631031000000d0"); // M02 pt dep with new std: 0.32, 0.0072, 0.5
    cuts.AddCutCalo("00000113","11111110b70322d0000","01631031000000d0"); // M02, pt dep with  0.27, 0.0072, 0.4
    cuts.AddCutCalo("00000113","11111110b70322e0000","01631031000000d0"); // M02, pt dep with  0.31, 0.0072, 0.5
    cuts.AddCutCalo("00000113","11111110b70322f0000","01631031000000d0"); // M02, pt dep with  0.36, 0.0072, 0.7
    cuts.AddCutCalo("00000113","11111110b70322m0000","01631031000000d0"); // M02, pt dep with  0.32, 0.0152, 0.5
    cuts.AddCutCalo("00000113","11111110b70322g0000","01631031000000d0"); // M02, pt dep with  0.37, 0.0072, 0.7
    cuts.AddCutCalo("00000113","11111110b70322h0000","01631031000000d0"); // M02, pT-dep with  0.30, 0.0072, 0.5
    cuts.AddCutCalo("00000113","11111110b70322i0000","01631031000000d0"); // M02, pT-dep with  0.35, 0.0072, 0.7
    cuts.AddCutCalo("00000113","11111110b70322j0000","01631031000000d0"); // M02, pT-dep with  0.25, 0.0072, 0.39
    cuts.AddCutCalo("00000113","11111110b70322r0000","01631031000000d0"); // M02, pT-dep with  0.25, 0.0072, 0.5
    cuts.AddCutCalo("00000113","11111110b70322s0000","01631031000000d0"); // M02, pT-dep with  0.32, 0.0238, 0.7
    cuts.AddCutCalo("00000113","11111110b70322n0000","01631031000000d0"); // M02, pT-dep with  0.32, 0.0238, 0.7
  } else if (trainConfig == 204){ //EMCAL minNCells,with/without TRD variation
    cuts.AddCutCalo("00000113","11111110b70312l0000","01631031000000d0"); //n cells >= 1
    cuts.AddCutCalo("00000113","11111110b70332l0000","01631031000000d0"); //n cells >= 3
    cuts.AddCutCalo("00000113","11121110b70322l0000","01631031000000d0"); //only modules with TRD infront
    cuts.AddCutCalo("00000113","11113110b70322l0000","01631031000000d0"); //no modules with TRD infront
  } else if (trainConfig == 205){  // trackMatching variations
    cuts.AddCutCalo("00000113","11111110b60322l0000","01631031000000d0"); //
    cuts.AddCutCalo("00000113","11111110b80322l0000","01631031000000d0"); //
    cuts.AddCutCalo("00000113","11111110b90322l0000","01631031000000d0"); //
    cuts.AddCutCalo("00000113","11111110b30322l0000","01631031000000d0"); //
    cuts.AddCutCalo("00000113","11111110b00322l0000","01631031000000d0"); //
  } else if (trainConfig == 207){ // EMCAL clusters open angle variation
    cuts.AddCutCalo("00000113","11111110b70322l0000","01631031000000a0"); // min open angle - 0. + 1 cell diag
    cuts.AddCutCalo("00000113","11111110b70322l0000","01631031000000d0"); // min open angle - 0.017 + 1 cell diag
    cuts.AddCutCalo("00000113","11111110b70322l0000","01631031000000b0"); // min open angle - 0.0152 + 1 cell diag
    cuts.AddCutCalo("00000113","11111110b70322l0000","01631031000000c0"); // min open angle - 0.016 + 1 cell diag
    cuts.AddCutCalo("00000113","11111110b70322l0000","01631031000000e0"); // min open angle - 0.018 + 1 cell diag
    cuts.AddCutCalo("00000113","11111110b70322l0000","01631031000000f0"); // min open angle - 0.019 + 1 cell diag
  } else if (trainConfig == 208){ // EMCAL clusters pp 7 TeV, Different DistanceToBadChannels
    cuts.AddCutCalo("00000113","11111111b70322l0000","01631031000000d0"); //
    cuts.AddCutCalo("00000113","11111112b70322l0000","01631031000000d0"); //
    cuts.AddCutCalo("00000113","11111113b70322l0000","01631031000000d0"); //
    cuts.AddCutCalo("00000113","11111115b70322l0000","01631031000000d0"); //
    cuts.AddCutCalo("00000113","11111116b70322l0000","01631031000000d0"); //
  } else if (trainConfig == 209){ // EMCAL clusters pp 7 TeV, Different NonLinearities
    cuts.AddCutCalo("00000113","11111010b70322l0000","01631031000000d0"); // NonLinearity kSDMv5
    cuts.AddCutCalo("00000113","11111130b70322l0000","01631031000000d0"); // NonLinearity kTestBeamv2 + LHC12 ConvCalo
    cuts.AddCutCalo("00000113","11111140b70322l0000","01631031000000d0"); // NonLinearity kTestBeamv2 + LHC12 Calo
  } else if (trainConfig == 210){ // EMCAL clusters pp 7 TeV, Different NonLinearities
    cuts.AddCutCalo("00000113","11111110b70322l0000","01631031000000d0"); // NonLinearity LHC12 ConvCalo
    cuts.AddCutCalo("00000113","11111120b70322l0000","01631031000000d0"); // NonLinearity LHC12 Calo
    cuts.AddCutCalo("00000113","11111210b70322l0000","01631031000000d0"); // NonLinearity LHC12 ConvCalo MassRatioFits
    cuts.AddCutCalo("00000113","11111220b70322l0000","01631031000000d0"); // NonLinearity LHC12 Calo MassRatioFits
    cuts.AddCutCalo("00000113","11111000b70322l0000","01631031000000d0"); // NonLinearity none
  } else if (trainConfig == 211){  // EMCAL clusters, different triggers no NonLinearity
    cuts.AddCutCalo("00000113","11111000b70322l0000","01631031000000d0");
    cuts.AddCutCalo("00052113","11111000b70322l0000","01631031000000d0"); // EMC7
    cuts.AddCutCalo("00081113","11111000b70322l0000","01631031000000d0"); // EMCEG1,
  } else if (trainConfig == 212){ // EMCAL clusters, exotic cut var
    cuts.AddCutCalo("00000113","11111110b72322l0000","01631031000000d0"); //
    cuts.AddCutCalo("00000113","11111110b73322l0000","01631031000000d0"); //
    cuts.AddCutCalo("00000113","11111110b75322l0000","01631031000000d0"); //
    cuts.AddCutCalo("00000113","11111110b77322l0000","01631031000000d0"); //
    cuts.AddCutCalo("00000113","11111110b79322l0000","01631031000000d0"); //
  } else if (trainConfig == 213){  // trackMatching variations
    cuts.AddCutCalo("00000113","11111110b20322l0000","01631031000000d0"); //
    cuts.AddCutCalo("00000113","11111110b30322l0000","01631031000000d0"); //
    cuts.AddCutCalo("00000113","11111110b40322l0000","01631031000000d0"); //
    cuts.AddCutCalo("00000113","11111110b50322l0000","01631031000000d0"); //
  } else if (trainConfig == 214){ // EMCAL clusters pp 7 TeV
    cuts.AddCutCalo("00010013","11111110b70322l0000","01631031000000d0"); // std - no pileup cut
  } else if (trainConfig == 215){ // EMCAL clusters pp 7 TeV
    cuts.AddCutCalo("00000113","11111110b70322l0000","01631051000000d0"); // alpha
    cuts.AddCutCalo("00000113","11111110b70322l0000","01631061000000d0"); //
  } else if (trainConfig == 216){ // EMCAL clusters pp 7 TeV
    cuts.AddCutCalo("00000113","11111110b70322l0000","01633031000000d0"); // rapidity
    cuts.AddCutCalo("00000113","11111110b70322l0000","01634031000000d0"); //


  // only std cuts
  } else if (trainConfig == 219){ // EMCAL clusters pp 7 TeV
    cuts.AddCutCalo("00000113","11111110b7032220000","01631031000000d0"); // std
  } else if (trainConfig == 220){ // EMCAL clusters pp 7 TeV
    cuts.AddCutCalo("00000113","11111110b70322l0000","01631031000000d0"); // std - dirGAMMA

    // std
    // ATTENTION: adapted for dirGamma - ADJUSTED M02 -> l
  } else if (trainConfig == 221){ // EMCAL clusters pp 7 TeV
    cuts.AddCutCalo("00000113","11111110b70322l0000","01631031000000d0"); // std NL
  } else if (trainConfig == 222){ // EMCAL clusters pp 7 TeV
    cuts.AddCutCalo("00000113","11111310b70322l0000","01631031000000d0"); // std TB NL

  } else if (trainConfig == 223){ // EMCAL clusters pp 7 TeV
    cuts.AddCutCalo("00000113","11111110b7032220000","01631031000000d0"); // std
  } else if (trainConfig == 224){ // EMCAL clusters pp 7 TeV
    cuts.AddCutCalo("00000113","11111110b7032220000","01631031000000d0"); // std
  } else if (trainConfig == 225){ // EMCAL clusters pp 7 TeV
    cuts.AddCutCalo("00000113","11111110b7032220000","01631031000000d0"); // std
  } else if (trainConfig == 226){ // EMCAL clusters pp 7 TeV
    cuts.AddCutCalo("00000113","11111110b7032220000","01631031000000d0"); // std

  // variations for omega analysis
  } else if (trainConfig == 230){  // standard
    cuts.AddCutCalo("00000113","1111111047032230000","01631031000000d0");
  } else if (trainConfig == 231){ // pileup
    cuts.AddCutCalo("00000013","1111111047032230000","01631031000000d0"); // rmeove pileup
  } else if (trainConfig == 232){ // nonlin
    cuts.AddCutCalo("00000113","1111112047032230000","01631031000000d0"); // NonLinearity ConvCalo - kTestBeamv3 + shifting MC
    cuts.AddCutCalo("00000113","1111113047032230000","01631031000000d0"); // NonLinearity pp Calo - only shifting MC - no timing cut
    cuts.AddCutCalo("00000113","1111121047032230000","01631031000000d0"); // NonLinearity pp ConvCalo - only shifting MC - no timing cut (Fits)
    cuts.AddCutCalo("00000113","1111122047032230000","01631031000000d0"); // NonLinearity pp ConvCalo - only shifting MC - no timing cut (Fits)
    cuts.AddCutCalo("00000113","1111123047032230000","01631031000000d0");  // NonLinearity ConvCalo - kTestBeamv3 + shifting MC
  } else if (trainConfig == 233){ // timing
    cuts.AddCutCalo("00000113","1111111037032230000","01631031000000d0"); // timing diff
    cuts.AddCutCalo("00000113","1111111057032230000","01631031000000d0"); // timing diff
    cuts.AddCutCalo("00000113","1111111077032230000","01631031000000d0"); // timing diff
    cuts.AddCutCalo("00000113","1111111097032230000","01631031000000d0"); // timing diff
  } else if (trainConfig == 234){ // TrackMatching
    cuts.AddCutCalo("00000113","1111111046032230000","01631031000000d0");
    cuts.AddCutCalo("00000113","1111111048032230000","01631031000000d0");
    cuts.AddCutCalo("00000113","1111111049032230000","01631031000000d0");
    cuts.AddCutCalo("00000113","111111104a032230000","01631031000000d0");
    cuts.AddCutCalo("00000113","111111104b032230000","01631031000000d0");
    cuts.AddCutCalo("00000113","1111111043032230000","01631031000000d0");
    cuts.AddCutCalo("00000113","111111104f032230000","01631031000000d0");
  } else if (trainConfig == 235){ // MinEnergy (of cluster)
    cuts.AddCutCalo("00000113","1111111047022230000","01631031000000d0"); // 0.6
    cuts.AddCutCalo("00000113","1111111047042230000","01631031000000d0"); // 0.8
    cuts.AddCutCalo("00000113","1111111047052230000","01631031000000d0"); // 0.9
  } else if (trainConfig == 236){ // MinNCells
    cuts.AddCutCalo("00000113","1111111047031230000","01631031000000d0"); // 1
    cuts.AddCutCalo("00000113","1111111047033230000","01631031000000d0"); // 3
  } else if (trainConfig == 237){ // MinMaxM02
    cuts.AddCutCalo("00000113","1111111047032330000","01631031000000d0"); // 0.2 - 0.5
    cuts.AddCutCalo("00000113","1111111047032130000","01631031000000d0"); // 0.002 - 0.5
    cuts.AddCutCalo("00000113","1111111047032240000","01631031000000d0"); // 0.1 - 0.4
    cuts.AddCutCalo("00000113","1111111047032220000","01631031000000d0"); // 0.1 - 0.7

  // LHC11cd configs V0OR and V0AND
  } else if (trainConfig == 250){  // EMCAL clusters 7 TeV LHC11 TM on, +-30ns, std TM, no NL
    cuts.AddCutCalo("00010113","1111100067032230000","01631031000000d0"); // VOAND
    cuts.AddCutCalo("00052113","1111100067032230000","01631031000000d0"); // EMC7
  } else if (trainConfig == 251){  // EMCAL clusters 7 TeV LHC10 MB only
    cuts.AddCutCalo("00000113","11111000b7032230000","01631031000000d0"); // VOAND
  } else if (trainConfig == 252){  // EMCAL clusters 7 TeV LHC10 MB only new 31 TB NL
    cuts.AddCutCalo("00000113","1111a3104f032230000","01631031000000d0"); // V0OR
  } else if (trainConfig == 254){  // QA for settings of omega analysis
    cuts.AddCutCalo("00000113","1111111047032230000","0163503800000000");
  // LHC11cd configs V0OR and V0AND with nonlinearity
  } else if (trainConfig == 260){  // EMCAL clusters 7 TeV LHC11 TM on, +-30ns, std TM, no NL
    cuts.AddCutCalo("00010113","1111111067032230000","01631031000000d0"); // VOAND
    cuts.AddCutCalo("00010113","1111121067032230000","01631031000000d0"); // VOAND
    cuts.AddCutCalo("00052113","1111111067032230000","01631031000000d0"); // EMC7
    cuts.AddCutCalo("00052113","1111121067032230000","01631031000000d0"); // EMC7
  } else if (trainConfig == 261){  // EMCAL clusters 7 TeV LHC10 MB only
    cuts.AddCutCalo("00000113","11111110b7032230000","01631031000000d0"); // VOAND
    cuts.AddCutCalo("00000113","11111210b7032230000","01631031000000d0"); // VOAND
  // LHC11cd configs V0OR and V0AND with nonlinearity
  } else if (trainConfig == 262){  // EMCAL clusters 7 TeV LHC11 TB NL TM on, +-30ns
    cuts.AddCutCalo("00010113","111110106f032230000","01631031000000d0"); // VOAND
    cuts.AddCutCalo("00052113","111110106f032230000","01631031000000d0"); // EMC7
  } else if (trainConfig == 263){  // EMCAL clusters 7 TeV LHC11 TB NL TM on, +-30ns
    cuts.AddCutCalo("00010c13","111113106f032230000","01631031000000d0"); // VOAND
    cuts.AddCutCalo("00052c13","111113106f032230000","01631031000000d0"); // EMC7
  } else if (trainConfig == 264){  // EMCAL clusters 7 TeV LHC11 TB NL TM on, +-30ns
    cuts.AddCutCalo("00010113","111113206f032230000","01631031000000d0"); // VOAND
    cuts.AddCutCalo("00052113","111113206f032230000","01631031000000d0"); // EMC7
  } else if (trainConfig == 265){  // EMCAL clusters 7 TeV LHC11 TB NL TM on, +-30ns
    cuts.AddCutCalo("00010113","111113306f032230000","01631031000000d0"); // VOAND
    cuts.AddCutCalo("00052113","111113306f032230000","01631031000000d0"); // EMC7
  } else if (trainConfig == 266){  // EMCAL clusters 7 TeV LHC11 TB NL TM on, +-30ns
    cuts.AddCutCalo("00010113","111113406f032230000","01631031000000d0"); // VOAND
    cuts.AddCutCalo("00052113","111113406f032230000","01631031000000d0"); // EMC7
  } else if (trainConfig == 267){  // EMCAL clusters 7 TeV LHC11 for Omega QA
    cuts.AddCutCalo("00010113","111111105f032230000","01631031000000d0"); // VOAND
    cuts.AddCutCalo("00052113","111111105f032230000","01631031000000d0"); // EMC7
  } else if (trainConfig == 268){  // EMCAL clusters 7 TeV LHC11 for Omega QA
    cuts.AddCutCalo("00010113","111113106f032230000","01631031000000d0"); // VOAND
    cuts.AddCutCalo("00052113","111113106f032230000","01631031000000d0"); // EMC7
  //multiple std dirGAMMA cuts for different studies
  } else if (trainConfig == 281){ // EMCAL clusters pp 7 TeV
    cuts.AddCutCalo("00000113","11111110b70322l0000","01631031000000d0"); // std
  } else if (trainConfig == 282){ // EMCAL clusters pp 7 TeV
    cuts.AddCutCalo("00000113","11111110b70322l0000","01631031000000d0"); // std
  } else if (trainConfig == 283){ // EMCAL clusters pp 7 TeV
    cuts.AddCutCalo("00000113","11111110b70322l0000","01631031000000d0"); // std
  } else if (trainConfig == 284){ // EMCAL clusters pp 7 TeV
    cuts.AddCutCalo("00000113","11111110b70322l0000","01631031000000d0"); // std


  // *****************************************************************************************************
  // ************************************* PHOS cuts ****************************************************
  // *****************************************************************************************************
  // pp 2.76 TeV
  } else if (trainConfig == 301) { //PHOS clusters
    cuts.AddCutCalo("00003113","2444400040033200000","0163803100000010"); //pp LHC11a with SDD, PHOS
    cuts.AddCutCalo("00010113","2444400040033200000","0163803100000010"); //pp LHC13g default MB
    cuts.AddCutCalo("00061113","2444400040033200000","0163803100000010"); //pp LHC11a PHI1
    cuts.AddCutCalo("00062113","2444400040033200000","0163803100000010"); //pp LHC11a PHI7
  } else if (trainConfig == 302){ // Validation PHOS
    cuts.AddCutCalo("00003113","2444400040033200000","0163803100000010");
  } else if (trainConfig == 303){ // PHOS clusters, without and with added signals
    cuts.AddCutCalo("00003113","2444400040033200000","0163803100000010");
    cuts.AddCutCalo("00003123","2444400040033200000","0163803100000010");

  // pp 7 TeV direct photon PHOS
  } else if (trainConfig == 351){
    cuts.AddCutCalo("00000113","2444400000013300000","0163803100000010"); // no nonlinearity
    cuts.AddCutCalo("00000113","2444401000013300000","0163803100000010"); // with PHOS nonlinearity
  } else if (trainConfig == 352){
    cuts.AddCutCalo("00000113","2444400040013300000","0163803100000010"); // 100ns timing cut, no track matching
    cuts.AddCutCalo("00000113","2444400043013300000","0163803100000010"); // 100ns timing cut
    cuts.AddCutCalo("00000113","2444400043013350000","0163803100000010"); // 100ns timing cut, M02<0.3
    cuts.AddCutCalo("00000113","2444400043013330000","0163803100000010"); // 100ns timing cut, M02<0.5
    cuts.AddCutCalo("00000113","2444400043013320000","0163803100000010"); // 100ns timing cut, M02<0.7
  } else if (trainConfig == 353){ // same as 352 but with PHOS nonlinearity
    cuts.AddCutCalo("00000113","2444401040013300000","0163803100000010"); // 100ns timing cut, no track matching
    cuts.AddCutCalo("00000113","2444401043013300000","0163803100000010"); // 100ns timing cut
    cuts.AddCutCalo("00000113","2444401043013350000","0163803100000010"); // 100ns timing cut, M02<0.3
    cuts.AddCutCalo("00000113","2444401043013330000","0163803100000010"); // 100ns timing cut, M02<0.5
    cuts.AddCutCalo("00000113","2444401043013320000","0163803100000010"); // 100ns timing cut, M02<0.7

  // pp 7 TeV PHOS
  } else if (trainConfig == 361){
    cuts.AddCutCalo("00000113","2444400000013300000","0163803100000010"); // QA
    cuts.AddCutCalo("00000113","2444400040013300000","0163803100000010"); // 100ns timing cut, no track matching
    cuts.AddCutCalo("00000113","2444400043013300000","0163803100000010"); // 100ns timing cut
  } else if (trainConfig == 362){
    cuts.AddCutCalo("00062113","2444400000013300000","0163803100000010"); // QA
    cuts.AddCutCalo("00062113","2444400040013300000","0163803100000010"); // 100ns timing cut, no track matching
    cuts.AddCutCalo("00062113","2444400043013300000","0163803100000010"); // 100ns timing cut
  } else if (trainConfig == 363){ // train config for bad channels and NonLin Variation
    cuts.AddCutCalo("00000113","2444400043012300000","0163803100000010"); // no NonLin
    cuts.AddCutCalo("00000113","2444401043012300000","0163803100000010"); // ext PHOS NonLin
    cuts.AddCutCalo("00000113","2444412043012300000","0163803100000010"); // own Calo NonLin
    cuts.AddCutCalo("00000113","2444411043012300000","0163803100000010"); // own ConvCalo NonLin (std)
    cuts.AddCutCalo("00000113","2444421043012300000","0163803100000010"); // own ConvCalo NonLin + ExtPHOS

  // pp 7 TeV PHOS for omega anlysis (ratio)
  } else if (trainConfig == 370){ // std
    cuts.AddCutCalo("00000113","2444411044013300000","0163803100000010");
  } else if (trainConfig == 371){ // remove pilup
    cuts.AddCutCalo("00000013","2444411044013300000","0163803100000010");
  } else if (trainConfig == 372){ // timing diff
    cuts.AddCutCalo("00000113","24444110b4013300000","0163803100000010");
    cuts.AddCutCalo("00000113","24444110c4013300000","0163803100000010");
    cuts.AddCutCalo("00000113","24444110d4013300000","0163803100000010");
    cuts.AddCutCalo("00000113","24444110e4013300000","0163803100000010");
    cuts.AddCutCalo("00000113","24444110f4013300000","0163803100000010");
  } else if (trainConfig == 373){ // non lin
    cuts.AddCutCalo("00000113","2444412044013300000","0163803100000010");
    cuts.AddCutCalo("00000113","2444421044013300000","0163803100000010");
  } else if (trainConfig == 374){ // track matching
    cuts.AddCutCalo("00000113","2444411041013300000","0163803100000010");
    cuts.AddCutCalo("00000113","2444411043013300000","0163803100000010");
    cuts.AddCutCalo("00000113","2444411045013300000","0163803100000010");
    cuts.AddCutCalo("00000113","2444411046013300000","0163803100000010");
    cuts.AddCutCalo("00000113","2444411047013300000","0163803100000010");
    cuts.AddCutCalo("00000113","2444411048013300000","0163803100000010");
    cuts.AddCutCalo("00000113","2444411049013300000","0163803100000010");
  } else if (trainConfig == 375){ // min energy
    cuts.AddCutCalo("00000113","2444411044023300000","0163803100000010");
    cuts.AddCutCalo("00000113","2444411044033300000","0163803100000010");
    cuts.AddCutCalo("00000113","2444411044043300000","0163803100000010");
    cuts.AddCutCalo("00000113","2444411044083300000","0163803100000010");
  } else if (trainConfig == 376){ // Min N of cells (std is 2)
    cuts.AddCutCalo("00000113","2444411044011300000","0163803100000010");
    cuts.AddCutCalo("00000113","2444411044013300000","0163803100000010");
    cuts.AddCutCalo("00000113","2444411044014300000","0163803100000010");
  } else if (trainConfig == 377){ // MinMaxM02 (std is >0.2)
    cuts.AddCutCalo("00000113","2444411044013200000","0163803100000010");
    cuts.AddCutCalo("00000113","2444411044013100000","0163803100000010");
  // PHOS @ 8 TeV
  } else if (trainConfig == 381){ // PHOS clusters
    cuts.AddCutCalo("00010113","2444400040013300000","0163103100000010");
  } else if (trainConfig == 382){ // PHOS clusters
    cuts.AddCutCalo("00062113","2444400040013300000","0163103100000010");

  // *********************************************************************************************************
  // 5 TeV  pp Run2 - EMC configurations
  // *********************************************************************************************************
  } else if (trainConfig == 401){ // EMCAL clusters
    cuts.AddCutCalo("00010113","1111100017032220000","01631031000000d0"); // 1000ns timing cut, no NL INT7
    cuts.AddCutCalo("00052013","1111100017032220000","01631031000000d0"); // 1000ns timing cut, no NL EMC7
    cuts.AddCutCalo("00085013","1111100017032220000","01631031000000d0"); // 1000ns timing cut, no NL EG2
    cuts.AddCutCalo("00083013","1111100017032220000","01631031000000d0"); // 1000ns timing cut, no NL EG1
  } else if (trainConfig == 402){ // EMCAL clusters
    cuts.AddCutCalo("00010113","1111100067032220000","01631031000000d0"); // -50ns, 30ns timing cut, no NL INT7
    cuts.AddCutCalo("00052013","1111100067032220000","01631031000000d0"); // -50ns, 30ns timing cut, no NL EMC7
    cuts.AddCutCalo("00085013","1111100067032220000","01631031000000d0"); // -50ns, 30ns timing cut, no NL EG2
    cuts.AddCutCalo("00083013","1111100067032220000","01631031000000d0"); // -50ns, 30ns timing cut, no NL EG1
  } else if (trainConfig == 403){ // EMCAL clusters - NonLin INT7
    cuts.AddCutCalo("00010113","1111111017032220000","01631031000000d0"); // 1000ns timing cut, NL kSDM PCMEMC
    cuts.AddCutCalo("00010113","1111112017032220000","01631031000000d0"); // 1000ns timing cut, NL kSDM EMC
    cuts.AddCutCalo("00010113","1111121017032220000","01631031000000d0"); // 1000ns timing cut, NL DExt PCMEMC
    cuts.AddCutCalo("00010113","1111122017032220000","01631031000000d0"); // 1000ns timing cut, NL DExt EMC
  } else if (trainConfig == 404){ // EMCAL clusters - NonLin INT7
    cuts.AddCutCalo("00010113","1111111067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL kSDM PCMEMC
    cuts.AddCutCalo("00010113","1111112067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL kSDM EMC
    cuts.AddCutCalo("00010113","1111121067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL DExt PCMEMC
    cuts.AddCutCalo("00010113","1111122067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL DExt EMC
  } else if (trainConfig == 405){ // EMCAL clusters - NonLin EMC7
    cuts.AddCutCalo("00052013","1111111017032220000","01631031000000d0"); // 1000ns timing cut, NL kSDM PCMEMC
    cuts.AddCutCalo("00052013","1111112017032220000","01631031000000d0"); // 1000ns timing cut, NL kSDM EMC
    cuts.AddCutCalo("00052013","1111121017032220000","01631031000000d0"); // 1000ns timing cut, NL DExt PCMEMC
    cuts.AddCutCalo("00052013","1111122017032220000","01631031000000d0"); // 1000ns timing cut, NL DExt EMC
  } else if (trainConfig == 406){ // EMCAL clusters - NonLin EMC7
    cuts.AddCutCalo("00052013","1111111067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL kSDM PCMEMC
    cuts.AddCutCalo("00052013","1111112067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL kSDM EMC
    cuts.AddCutCalo("00052013","1111121067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL DExt PCMEMC
    cuts.AddCutCalo("00052013","1111122067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL DExt EMC
  } else if (trainConfig == 407){ // EMCAL clusters - NonLin EG2
    cuts.AddCutCalo("00085013","1111111017032220000","01631031000000d0"); // 1000ns timing cut, NL kSDM PCMEMC
    cuts.AddCutCalo("00085013","1111112017032220000","01631031000000d0"); // 1000ns timing cut, NL kSDM EMC
    cuts.AddCutCalo("00085013","1111121017032220000","01631031000000d0"); // 1000ns timing cut, NL DExt PCMEMC
    cuts.AddCutCalo("00085013","1111122017032220000","01631031000000d0"); // 1000ns timing cut, NL DExt EMC
  } else if (trainConfig == 408){ // EMCAL clusters - NonLin EG2
    cuts.AddCutCalo("00085013","1111111067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL kSDM PCMEMC
    cuts.AddCutCalo("00085013","1111112067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL kSDM EMC
    cuts.AddCutCalo("00085013","1111121067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL DExt PCMEMC
    cuts.AddCutCalo("00085013","1111122067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL DExt EMC
  } else if (trainConfig == 409){ // EMCAL clusters - NonLin EG1
    cuts.AddCutCalo("00083013","1111111017032220000","01631031000000d0"); // 1000ns timing cut, NL kSDM PCMEMC
    cuts.AddCutCalo("00083013","1111112017032220000","01631031000000d0"); // 1000ns timing cut, NL kSDM EMC
    cuts.AddCutCalo("00083013","1111121017032220000","01631031000000d0"); // 1000ns timing cut, NL DExt PCMEMC
    cuts.AddCutCalo("00083013","1111122017032220000","01631031000000d0"); // 1000ns timing cut, NL DExt EMC
  } else if (trainConfig == 410){ // EMCAL clusters - NonLin EG1
    cuts.AddCutCalo("00083013","1111111067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL kSDM PCMEMC
    cuts.AddCutCalo("00083013","1111112067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL kSDM EMC
    cuts.AddCutCalo("00083013","1111121067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL DExt PCMEMC
    cuts.AddCutCalo("00083013","1111122067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL DExt EMC
  } else if (trainConfig == 411){ // EMCAL clusters - NonLin INT7
    cuts.AddCutCalo("00010113","1111113017032220000","01631031000000d0"); // 1000ns timing cut, NL kSDM PCMEMC + testbeam
    cuts.AddCutCalo("00010113","1111114017032220000","01631031000000d0"); // 1000ns timing cut, NL kSDM EMC + testbeam
    cuts.AddCutCalo("00010113","1111123017032220000","01631031000000d0"); // 1000ns timing cut, NL DExt PCMEMC + testbeam
    cuts.AddCutCalo("00010113","1111124017032220000","01631031000000d0"); // 1000ns timing cut, NL DExt EMC + testbeam
  } else if (trainConfig == 412){ // EMCAL clusters - NonLin INT7
    cuts.AddCutCalo("00010113","1111113067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL kSDM PCMEMC + testbeam
    cuts.AddCutCalo("00010113","1111114067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL kSDM EMC + testbeam
    cuts.AddCutCalo("00010113","1111123067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL DExt PCMEMC + testbeam
    cuts.AddCutCalo("00010113","1111124067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL DExt EMC + testbeam
  } else if (trainConfig == 413){ // EMCAL clusters - NonLin EMC7
    cuts.AddCutCalo("00052013","1111113017032220000","01631031000000d0"); // 1000ns timing cut, NL kSDM PCMEMC
    cuts.AddCutCalo("00052013","1111114017032220000","01631031000000d0"); // 1000ns timing cut, NL kSDM EMC
    cuts.AddCutCalo("00052013","1111123017032220000","01631031000000d0"); // 1000ns timing cut, NL DExt PCMEMC
    cuts.AddCutCalo("00052013","1111124017032220000","01631031000000d0"); // 1000ns timing cut, NL DExt EMC
  } else if (trainConfig == 414){ // EMCAL clusters - NonLin EMC7
    cuts.AddCutCalo("00052013","1111113067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL kSDM PCMEMC
    cuts.AddCutCalo("00052013","1111114067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL kSDM EMC
    cuts.AddCutCalo("00052013","1111123067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL DExt PCMEMC
    cuts.AddCutCalo("00052013","1111124067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL DExt EMC
  } else if (trainConfig == 415){ // EMCAL clusters - NonLin INT7
    cuts.AddCutCalo("00010113","1111102017032220000","01631031000000d0"); // 1000ns timing cut, NL BeamTest ONLY
    cuts.AddCutCalo("00010113","1111109017032220000","01631031000000d0"); // 1000ns timing cut, NL 8 TeV 11
    cuts.AddCutCalo("00010113","1111117017032220000","01631031000000d0"); // 1000ns timing cut, NL kSDM PCMEMC - fixed E-dependence
    cuts.AddCutCalo("00010113","1111127017032220000","01631031000000d0"); // 1000ns timing cut, NL DExt PCMEMC - fixed E-dependence
  } else if (trainConfig == 416){ // EMCAL clusters - NonLin INT7
    cuts.AddCutCalo("00010113","1111102067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL BeamTest ONLY
    cuts.AddCutCalo("00010113","1111109067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL 8 TeV 11
    cuts.AddCutCalo("00010113","1111117067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL kSDM PCMEMC - fixed E-dependence
    cuts.AddCutCalo("00010113","1111127067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL DExt PCMEMC - fixed E-dependence
  } else if (trainConfig == 417){ // EMCAL clusters - NonLin INT7 - with BeamTest
    cuts.AddCutCalo("00010113","1111113017032220000","01631031000000d0"); // 1000ns timing cut, NL kSDM PCMEMC+BTv3
    cuts.AddCutCalo("00010113","1111114017032220000","01631031000000d0"); // 1000ns timing cut, NL kSDM PCMEMC+BTv3
    cuts.AddCutCalo("00010113","1111123017032220000","01631031000000d0"); // 1000ns timing cut, NL DExt PCMEMC+BTv3
    cuts.AddCutCalo("00010113","1111124017032220000","01631031000000d0"); // 1000ns timing cut, NL DExt PCMEMC+BTv3
  } else if (trainConfig == 418){ // EMCAL clusters - NonLin INT7 - with BeamTest
    cuts.AddCutCalo("00010113","1111113067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL kSDM PCMEMC+BTv3
    cuts.AddCutCalo("00010113","1111114067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL kSDM PCMEMC+BTv3
    cuts.AddCutCalo("00010113","1111123067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL DExt PCMEMC+BTv3
    cuts.AddCutCalo("00010113","1111124067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL DExt PCMEMC+BTv3
  } else if (trainConfig == 420){ // EMCAL clusters no NL
    cuts.AddCutCalo("00010113","11111000a7032230000","01631031000000d0"); // -12.5ns, 13ns timing cut, no NL INT7
    cuts.AddCutCalo("00010113","11111060a7032230000","01631031000000d0"); // -12.5ns, 13ns timing cut, TB NL INT7
    cuts.AddCutCalo("00010113","11111650a7032230000","01631031000000d0"); // -12.5ns, 13ns timing cut, TB NL INT7
    cuts.AddCutCalo("00010113","11111660a7032230000","01631031000000d0"); // -12.5ns, 13ns timing cut, TB NL INT7
  } else if (trainConfig == 421){ // EMCAL clusters no NonLin
    cuts.AddCutCalo("00010113","11111610a7032230000","01631031000000d0"); // -12.5ns, 13ns timing cut, no NL INT7
    cuts.AddCutCalo("00010113","11111620a7032230000","01631031000000d0"); // -12.5ns, 13ns timing cut, TB NL INT7
    cuts.AddCutCalo("00010113","11111630a7032230000","01631031000000d0"); // -12.5ns, 13ns timing cut, TB NL INT7
    cuts.AddCutCalo("00010113","11111640a7032230000","01631031000000d0"); // -12.5ns, 13ns timing cut, TB NL INT7
  } else if (trainConfig == 422){ // EMCAL clusters no NonLin
    cuts.AddCutCalo("00010113","11111670a7032230000","01631031000000d0"); // -12.5ns, 13ns timing cut, TB NL INT7
    cuts.AddCutCalo("00010113","11111680a7032230000","01631031000000d0"); // -12.5ns, 13ns timing cut, TB NL INT7
  } else if (trainConfig == 423){ // EMCAL clusters no NonLin
    cuts.AddCutCalo("00010113","11111110a7032230000","01631031000000d0"); // -12.5ns, 13ns timing cut,
    cuts.AddCutCalo("00010113","11111120a7032230000","01631031000000d0"); // -12.5ns, 13ns timing cut,
    cuts.AddCutCalo("00010113","11111210a7032230000","01631031000000d0"); // -12.5ns, 13ns timing cut,
    cuts.AddCutCalo("00010113","11111220a7032230000","01631031000000d0"); // -12.5ns, 13ns timing cut,
    // 5 TeV JetJet configs without trackmatching
  } else if (trainConfig == 430){ // EMCAL clusters no NonLin
    cuts.AddCutCalo("00010113","1111100060032220000","01631031000000d0"); // -50ns, 30ns timing cut, no NL INT7
    cuts.AddCutCalo("00052013","1111100060032220000","01631031000000d0"); // -50ns, 30ns timing cut, no NL EMC7
    cuts.AddCutCalo("00085013","1111100060032220000","01631031000000d0"); // -50ns, 30ns timing cut, no NL EG2
    cuts.AddCutCalo("00083013","1111100060032220000","01631031000000d0"); // -50ns, 30ns timing cut, no NL EG1
  } else if (trainConfig == 431){ // EMCAL clusters - NonLin INT7 - no trackmatching for 5TeV JJ
    cuts.AddCutCalo("00010113","1111111060032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL kSDM PCMEMC
    cuts.AddCutCalo("00010113","1111112060032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL kSDM EMC
    cuts.AddCutCalo("00010113","1111121060032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL DExt PCMEMC
    cuts.AddCutCalo("00010113","1111122060032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL DExt EMC
  } else if (trainConfig == 432){ // EMCAL clusters - NonLin EMC7 - no trackmatching for 5TeV JJ
    cuts.AddCutCalo("00052013","1111111060032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL kSDM PCMEMC
    cuts.AddCutCalo("00052013","1111112060032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL kSDM EMC
    cuts.AddCutCalo("00052013","1111121060032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL DExt PCMEMC
    cuts.AddCutCalo("00052013","1111122060032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL DExt EMC
  } else if (trainConfig == 433){ // EMCAL clusters - NonLin EG2  - no trackmatching for 5TeV JJ
    cuts.AddCutCalo("00085013","1111111060032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL kSDM PCMEMC
    cuts.AddCutCalo("00085013","1111112060032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL kSDM EMC
    cuts.AddCutCalo("00085013","1111121060032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL DExt PCMEMC
    cuts.AddCutCalo("00085013","1111122060032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL DExt EMC
  } else if (trainConfig == 434){ // EMCAL clusters - NonLin EG1  - no trackmatching for 5TeV JJ
    cuts.AddCutCalo("00083013","1111111060032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL kSDM PCMEMC
    cuts.AddCutCalo("00083013","1111112060032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL kSDM EMC
    cuts.AddCutCalo("00083013","1111121060032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL DExt PCMEMC
    cuts.AddCutCalo("00083013","1111122060032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL DExt EMC

  } else if (trainConfig == 435){ // EMCAL clusters pp 5 TeV V0M mult cuts
    cuts.AddCutCalo("m0110113","4117911077032230000","01631031000000d0"); // std 0-1%
    cuts.AddCutCalo("m1510113","4117911077032230000","01631031000000d0"); // std 1-5%
    cuts.AddCutCalo("m5k10113","4117911077032230000","01631031000000d0"); // std 5-20%
    cuts.AddCutCalo("n2410113","4117911077032230000","01631031000000d0"); // std 20-40%
    cuts.AddCutCalo("n4710113","4117911077032230000","01631031000000d0"); // std 40-70%
    cuts.AddCutCalo("n7a10113","4117911077032230000","01631031000000d0"); // std 70-100%
  } else if (trainConfig == 436){ // EMCAL clusters pp 5 TeV SPD mult cuts
    cuts.AddCutCalo("o0110113","4117911077032230000","01631031000000d0"); // std 0-1%
    cuts.AddCutCalo("o0210113","4117911077032230000","01631031000000d0"); // std 0-2%
    cuts.AddCutCalo("o0510113","4117911077032230000","01631031000000d0"); // std 2-5%
    cuts.AddCutCalo("o5k10113","4117911077032230000","01631031000000d0"); // std 5-20%
    cuts.AddCutCalo("p2610113","4117911077032230000","01631031000000d0"); // std 20-60%
    cuts.AddCutCalo("p6a10113","4117911077032230000","01631031000000d0"); // std 60-100%
  } else if (trainConfig == 437){ // EMCAL clusters pp 5 TeV V0M mult cuts - EMC7
    cuts.AddCutCalo("m01a1113","4117911070032230000","01631031000000d0"); // std 0-1%
    cuts.AddCutCalo("m15a1113","4117911070032230000","01631031000000d0"); // std 1-5%
    cuts.AddCutCalo("m5ka1113","4117911070032230000","01631031000000d0"); // std 5-20%
    cuts.AddCutCalo("n24a1113","4117911070032230000","01631031000000d0"); // std 20-40%
    cuts.AddCutCalo("n47a1113","4117911070032230000","01631031000000d0"); // std 40-70%
    cuts.AddCutCalo("n7aa1113","4117911070032230000","01631031000000d0"); // std 70-100%
  } else if (trainConfig == 438){ // EMCAL clusters pp 5 TeV V0M mult cuts - EG2
    cuts.AddCutCalo("m01a2113","4117911070032230000","01631031000000d0"); // std 0-1%
    cuts.AddCutCalo("m15a2113","4117911070032230000","01631031000000d0"); // std 1-5%
    cuts.AddCutCalo("m5ka2113","4117911070032230000","01631031000000d0"); // std 5-20%
    cuts.AddCutCalo("n24a2113","4117911070032230000","01631031000000d0"); // std 20-40%
    cuts.AddCutCalo("n47a2113","4117911070032230000","01631031000000d0"); // std 40-70%
    cuts.AddCutCalo("n7aa2113","4117911070032230000","01631031000000d0"); // std 70-100%
  } else if (trainConfig == 439){ // EMCAL clusters pp 5 TeV V0M mult cuts
    cuts.AddCutCalo("m0110113","4117911070032230000","01631031000000d0"); // std 0-1%
    cuts.AddCutCalo("m1510113","4117911070032230000","01631031000000d0"); // std 1-5%
    cuts.AddCutCalo("m5k10113","4117911070032230000","01631031000000d0"); // std 5-20%
    cuts.AddCutCalo("n2410113","4117911070032230000","01631031000000d0"); // std 20-40%
    cuts.AddCutCalo("n4710113","4117911070032230000","01631031000000d0"); // std 40-70%
    cuts.AddCutCalo("n7a10113","4117911070032230000","01631031000000d0"); // std 70-100%

    //Sphericity Cuts
  } else if (trainConfig == 440){ // EMCAL clusters pp 5 TeV Sphericity Cuts
    cuts.AddCutCalo("00010113","11111110a7032230000","01631031000000d0"); // std
    cuts.AddCutCalo("h0510113","11111110a7032230000","01631031000000d0"); // std
    cuts.AddCutCalo("h5a10113","11111110a7032230000","01631031000000d0"); // std
    cuts.AddCutCalo("h0a10113","11111110a7032230000","01631031000000d0"); // std
    cuts.AddCutCalo("h0310113","11111110a7032230000","01631031000000d0"); // std
    cuts.AddCutCalo("h7a10113","11111110a7032230000","01631031000000d0"); // std
  } else if (trainConfig == 441){ // EMCAL clusters pp 5 TeV Sphericity Cuts
    cuts.AddCutCalo("i0a10113","11111110a7032230000","01631031000000d0"); // std
    cuts.AddCutCalo("i0510113","11111110a7032230000","01631031000000d0"); // std
    cuts.AddCutCalo("i5a10113","11111110a7032230000","01631031000000d0"); // std
    cuts.AddCutCalo("j0a10113","11111110a7032230000","01631031000000d0"); // std
    cuts.AddCutCalo("j0510113","11111110a7032230000","01631031000000d0"); // std
    cuts.AddCutCalo("j5a10113","11111110a7032230000","01631031000000d0"); // std
  } else if (trainConfig == 442){ // EMCAL clusters pp 5 TeV opening angle studies
    cuts.AddCutCalo("00010113","11111110a7032230000","01631031000000d0"); // std
    cuts.AddCutCalo("00010113","11111110a7032230000","01631031000000d3"); // std pi0 & eta with max opening angle cut
    cuts.AddCutCalo("00010113","11111110a7032230000","01631031000000h3"); // std eta excluding the pi0 & max cut
  } else if (trainConfig == 443){ // EMCAL clusters pp 5 TeV mixing studies
    cuts.AddCutCalo("00010113","11111110a7032230000","01631031000000d0"); // std
    cuts.AddCutCalo("00010113","11111110a7032230000","0i631031000000d0"); // std
  } else if (trainConfig == 444){ // EMCAL clusters pp 5 TeV studies for flat energy subtraction
    cuts.AddCutCalo("00010113","11111110a7032230000","01631031000000d0"); // std
    cuts.AddCutCalo("00010113","11111110a00f2230000","01631031000000d0"); // std
  } else if (trainConfig == 445){ // EMCAL clusters pp 5 TeV Sphericity Cuts in jet events
    cuts.AddCutCalo("00010113","4117911077032230000","21631031000000d0"); // std
    cuts.AddCutCalo("h0510113","4117911077032230000","21631031000000d0"); // std
    cuts.AddCutCalo("h5a10113","4117911077032230000","21631031000000d0"); // std
    cuts.AddCutCalo("h0a10113","4117911077032230000","21631031000000d0"); // std
    cuts.AddCutCalo("h0310113","4117911077032230000","21631031000000d0"); // std
    cuts.AddCutCalo("h7a10113","4117911077032230000","21631031000000d0"); // std
  } else if (trainConfig == 446){ // EMCAL clusters pp 5 TeV Sphericity Cuts in soft events
    cuts.AddCutCalo("00010113","11111110a7032230000","21631031000000d0"); // std
    cuts.AddCutCalo("h0510113","11111110a7032230000","21631031000000d0"); // std
    cuts.AddCutCalo("h5a10113","11111110a7032230000","21631031000000d0"); // std
    cuts.AddCutCalo("h0a10113","11111110a7032230000","21631031000000d0"); // std
    cuts.AddCutCalo("h0310113","11111110a7032230000","21631031000000d0"); // std
    cuts.AddCutCalo("h7a10113","11111110a7032230000","21631031000000d0"); // std
  } else if (trainConfig == 447){ // EMCAL clusters pp 5 TeV True Sphericity Cuts
    cuts.AddCutCalo("h0a10113","11111110a7032230000","01631031000000d0"); // std
  } else if (trainConfig == 448){ // EMCAL clusters pp 5 TeV Sphericity axis cuts
    cuts.AddCutCalo("k0310113","11111110a7032230000","01631031000000d0"); // std
    cuts.AddCutCalo("l0310113","11111110a7032230000","01631031000000d0"); // std

  } else if (trainConfig == 449){ // EMCAL standard cuts, different triggers - no TM
    cuts.AddCutCalo("00010113","11111110a0032230000","01631031000000d0"); // -50ns, 30ns timing cut, NL kSDM PCMEMC INT7
    cuts.AddCutCalo("00052013","11111110a0032230000","01631031000000d0"); // -50ns, 30ns timing cut, NL kSDM PCMEMC
    cuts.AddCutCalo("00085013","11111110a0032230000","01631031000000d0"); // -50ns, 30ns timing cut, NL kSDM PCMEMC EG2
    cuts.AddCutCalo("00083013","11111110a0032230000","01631031000000d0"); // -50ns, 30ns timing cut, NL kSDM PCMEMC EG1
  } else if (trainConfig == 450){ // EMCAL standard cuts, different triggers
    cuts.AddCutCalo("00010113","11111110a7032230000","01631031000000d0"); // -50ns, 30ns timing cut, NL kSDM PCMEMC INT7
    cuts.AddCutCalo("00052013","11111110a7032230000","01631031000000d0"); // -50ns, 30ns timing cut, NL kSDM PCMEMC
    cuts.AddCutCalo("00085013","11111110a7032230000","01631031000000d0"); // -50ns, 30ns timing cut, NL kSDM PCMEMC EG2
    cuts.AddCutCalo("00083013","11111110a7032230000","01631031000000d0"); // -50ns, 30ns timing cut, NL kSDM PCMEMC EG1
  } else if (trainConfig == 451){ // EMCAL syst 1/7
    cuts.AddCutCalo("00010113","11111110a7022230000","01631031000000d0"); // min energy cluster variation 1  600 MeV
    cuts.AddCutCalo("00010113","11111110a7042230000","01631031000000d0"); // min energy cluster variation 2  800 MeV
    cuts.AddCutCalo("00010113","11111110a7052230000","01631031000000d0"); // min energy cluster variation 3  900 MeV
    cuts.AddCutCalo("00010113","11111110a7032230000","01631061000000d0"); // alpha meson variation 1 0<alpha<0.8
    cuts.AddCutCalo("00010113","11111110a7032230000","01631051000000d0"); // alpha meson variation 2 0<alpha<0.75
  } else if (trainConfig == 452){ // EMCAL syst 2/7
    cuts.AddCutCalo("00010113","11111110a7032230000","01631031000000d0"); // min/max M02  0.1<M<0.5
    cuts.AddCutCalo("00010113","11111110a7032200000","01631031000000d0"); // min/max M02  0.1<M<100
    cuts.AddCutCalo("00010113","11111110a7032250000","01631031000000d0"); // min/max M02  0.1<M<0.3
    cuts.AddCutCalo("00010113","11111110a7032260000","01631031000000d0"); // min/max M02  0.1<M<0.27
  } else if (trainConfig == 453){ // EMCAL syst 3/7
    cuts.AddCutCalo("00010113","11121110a7032230000","01631031000000d0"); // only modules with TRD infront
    cuts.AddCutCalo("00010113","11113110a7032230000","01631031000000d0"); // no modules with TRD infront
    cuts.AddCutCalo("00010113","11111110a7032230000","01633031000000d0"); // rapidity variation  y<0.6
    cuts.AddCutCalo("00010113","11111110a7032230000","01634031000000d0"); // rapidity variation  y<0.5
  } else if (trainConfig == 454){ // EMCAL syst 4/7
    cuts.AddCutCalo("00010113","11111110a7032230000","0163103100000040"); // min opening angle 0.0152
    cuts.AddCutCalo("00010113","11111110a7032230000","0163103100000060"); // min opening angle 0.017
    cuts.AddCutCalo("00010113","11111110a7032230000","0163103100000070"); // min opening angle 0.016
    cuts.AddCutCalo("00010113","11111110a7032230000","0163103100000080"); // min opening angle 0.018
    cuts.AddCutCalo("00010113","11111110a7032230000","0163103100000090"); // min opening angle 0.019
  } else if (trainConfig == 455){ // EMCAL syst 5/7
    cuts.AddCutCalo("00010113","11111110a3032230000","01631031000000d0"); // fixed window
    cuts.AddCutCalo("00010113","11111110a6032230000","01631031000000d0"); // tm pt dependent var 1
    cuts.AddCutCalo("00010113","11111110a8032230000","01631031000000d0"); // tm pt dependent var 2
    cuts.AddCutCalo("00010113","11111110a9032230000","01631031000000d0"); // tm pt dependent var 3
  } else if (trainConfig == 456){ // EMCAL syst 6/7
    cuts.AddCutCalo("00010113","11111010a7032230000","01631031000000d0"); // NL variation, TB
    cuts.AddCutCalo("00010113","11111110a7032230000","01631031000000d0"); // NL variation, ConvCalo
    cuts.AddCutCalo("00010113","11111120a7032230000","01631031000000d0"); // NL variation, Calo
  } else if (trainConfig == 457){ // EMCAL syst 7/7
    cuts.AddCutCalo("00010113","1111111047032230000","01631031000000d0"); // cluster timing cut
    cuts.AddCutCalo("00010113","1111111057032230000","01631031000000d0"); // cluster timing cut
    cuts.AddCutCalo("00010113","1111111077032230000","01631031000000d0"); // cluster timing cut
    cuts.AddCutCalo("00010113","11111110a7032230000","01631031000000d0"); // cluster timing cut
  } else if (trainConfig == 458){ // DEFAULT 2018 oct 31 no NL
    cuts.AddCutCalo("00010113","11111000a7032230000","01631031000000d0"); // cluster timing cut
  } else if (trainConfig == 459){ // DEFAULT 2018 oct 31
    cuts.AddCutCalo("00010113","11111110a7032230000","01631031000000d0"); // cluster timing cut

  } else if (trainConfig == 460){ // INT7 EMCAL standard cut but with E/p TM veto
    cuts.AddCutCalo("00010113","1111111067032220000","01631031000000d0"); // std INT7
    cuts.AddCutCalo("00010113","111111106c032220000","01631031000000d0"); // fEOverPMax = 9e9
    cuts.AddCutCalo("00010113","111111106d032220000","01631031000000d0"); // fEOverPMax = 3.0
    cuts.AddCutCalo("00010113","111111106e032220000","01631031000000d0"); // fEOverPMax = 2.0
  } else if (trainConfig == 461){ // INT7 EMCAL standard cut but with E/p TM veto
    cuts.AddCutCalo("00010113","111111106f032220000","01631031000000d0"); // fEOverPMax = 1.75
    cuts.AddCutCalo("00010113","111111106g032220000","01631031000000d0"); // fEOverPMax = 1.5
    cuts.AddCutCalo("00010113","111111106h032220000","01631031000000d0"); // fEOverPMax = 1.25
  } else if (trainConfig == 462){ // EMC7 EMCAL standard cut but with E/p TM veto
    cuts.AddCutCalo("00052113","1111111067032220000","01631031000000d0"); // std EMC7
    cuts.AddCutCalo("00052113","111111106c032220000","01631031000000d0"); // fEOverPMax = 9e9
    cuts.AddCutCalo("00052113","111111106d032220000","01631031000000d0"); // fEOverPMax = 3.0
    cuts.AddCutCalo("00052113","111111106e032220000","01631031000000d0"); // fEOverPMax = 2.0
  } else if (trainConfig == 463){ // EMC7 EMCAL standard cut but with E/p TM veto
    cuts.AddCutCalo("00052113","111111106f032220000","01631031000000d0"); // fEOverPMax = 1.75
    cuts.AddCutCalo("00052113","111111106g032220000","01631031000000d0"); // fEOverPMax = 1.5
    cuts.AddCutCalo("00052113","111111106h032220000","01631031000000d0"); // fEOverPMax = 1.25
  } else if (trainConfig == 464){ // EG1 EMCAL standard cut but with E/p TM veto
    cuts.AddCutCalo("00083113","1111111067032220000","01631031000000d0"); // std EG1
    cuts.AddCutCalo("00083113","111111106c032220000","01631031000000d0"); // fEOverPMax = 9e9
    cuts.AddCutCalo("00083113","111111106d032220000","01631031000000d0"); // fEOverPMax = 3.0
    cuts.AddCutCalo("00083113","111111106e032220000","01631031000000d0"); // fEOverPMax = 2.0
  } else if (trainConfig == 465){ // EG1 EMCAL standard cut but with E/p TM veto
    cuts.AddCutCalo("00083113","111111106f032220000","01631031000000d0"); // fEOverPMax = 1.75
    cuts.AddCutCalo("00083113","111111106g032220000","01631031000000d0"); // fEOverPMax = 1.5
    cuts.AddCutCalo("00083113","111111106h032220000","01631031000000d0"); // fEOverPMax = 1.25
  } else if (trainConfig == 466){ // EG2 EMCAL standard cut but with E/p TM veto
    cuts.AddCutCalo("00085113","1111111067032220000","01631031000000d0"); // std EG2
    cuts.AddCutCalo("00085113","111111106c032220000","01631031000000d0"); // fEOverPMax = 9e9
    cuts.AddCutCalo("00085113","111111106d032220000","01631031000000d0"); // fEOverPMax = 3.0
    cuts.AddCutCalo("00085113","111111106e032220000","01631031000000d0"); // fEOverPMax = 2.0
  } else if (trainConfig == 467){ // EG2 EMCAL standard cut but with E/p TM veto
    cuts.AddCutCalo("00085113","111111106f032220000","01631031000000d0"); // fEOverPMax = 1.75
    cuts.AddCutCalo("00085113","111111106g032220000","01631031000000d0"); // fEOverPMax = 1.5
    cuts.AddCutCalo("00085113","111111106h032220000","01631031000000d0"); // fEOverPMax = 1.25

  } else if (trainConfig == 480){ // INT7 EMCAL standard cut with triggers - NO TM - CALO+CALOFAST readout triggers
    cuts.AddCutCalo("000a0113","11111110a0032230000","01631031000000d0"); // std INT7
    cuts.AddCutCalo("000a1113","11111110a0032230000","01631031000000d0"); // std EMC7
    cuts.AddCutCalo("000a2113","11111110a0032230000","01631031000000d0"); // std EG2
    cuts.AddCutCalo("000a3113","11111110a0032230000","01631031000000d0"); // std EG1
  } else if (trainConfig == 481){ // INT7 EMCAL standard cut with triggers - NO TM
    cuts.AddCutCalo("00010113","11111110a0032230000","01631031000000d0"); // std INT7
    cuts.AddCutCalo("00052113","11111110a0032230000","01631031000000d0"); // std EMC7
    cuts.AddCutCalo("00085113","11111110a0032230000","01631031000000d0"); // std EG2
    cuts.AddCutCalo("00083113","11111110a0032230000","01631031000000d0"); // std EG1
  } else if (trainConfig == 482){ // INT7 EMCAL standard cut with triggers - NO TM - CALO+CALOFAST readout triggers
    cuts.AddCutCalo("000a0113","4117932090032230000","01631031000000d0"); // std INT7
    cuts.AddCutCalo("000a1113","4117932090032230000","01631031000000d0"); // std EMC7
    cuts.AddCutCalo("000a2113","4117932090032230000","01631031000000d0"); // std EG2

    // *********************************************************************************************************
  // 13 TeV  pp Run2 - EMC configurations
  // *********************************************************************************************************
  } else if (trainConfig == 500){ // EMCAL clusters
    cuts.AddCutCalo("00010113","1111100017032220000","01631031000000d0"); // 1000ns timing cut, no NL INT7
    cuts.AddCutCalo("00010113","1111100067032220000","01631031000000d0"); // -30ns, 35ns timing cut, no NL INT7
  } else if (trainConfig == 501){ // EMCAL clusters
    cuts.AddCutCalo("00052113","1111100017032220000","01631031000000d0"); // 1000ns timing cut, no NL EMC7
    cuts.AddCutCalo("00085113","1111100017032220000","01631031000000d0"); // 1000ns timing cut, no NL EG2
    cuts.AddCutCalo("00083113","1111100017032220000","01631031000000d0"); // 1000ns timing cut, no NL EG1
  } else if (trainConfig == 502){ // EMCAL clusters
    cuts.AddCutCalo("00052113","1111100067032220000","01631031000000d0"); // -50ns, 30ns timing cut, no NL EMC7
    cuts.AddCutCalo("00085113","1111100067032220000","01631031000000d0"); // -50ns, 30ns timing cut, no NL EG2
    cuts.AddCutCalo("00083113","1111100067032220000","01631031000000d0"); // -50ns, 30ns timing cut, no NL EG1
  } else if (trainConfig == 503){ // EMCAL clusters
    cuts.AddCutCalo("00010113","1111112067032220000","01631031000000d0"); // -30ns, 35ns timing cut, NL1 INT7
    cuts.AddCutCalo("00010113","1111122067032220000","01631031000000d0"); // -30ns, 35ns timing cut, NL2 INT7
    cuts.AddCutCalo("00010113","1111111067032220000","01631031000000d0"); // -30ns, 35ns timing cut, NL1 INT7
    cuts.AddCutCalo("00010113","1111121067032220000","01631031000000d0"); // -30ns, 35ns timing cut, NL2 INT7
  } else if (trainConfig == 504){ // EMCAL clusters
    cuts.AddCutCalo("00052113","1111112067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL1 EMC7
    cuts.AddCutCalo("00052113","1111122067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL2 EMC7
    cuts.AddCutCalo("00052113","1111111067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL3 EMC7
    cuts.AddCutCalo("00052113","1111121067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL4 EMC7
  } else if (trainConfig == 505){ // EMCAL clusters
    cuts.AddCutCalo("00085113","1111112067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL1 EG2
    cuts.AddCutCalo("00085113","1111122067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL2 EG2
    cuts.AddCutCalo("00085113","1111111067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL3 EG2
    cuts.AddCutCalo("00085113","1111121067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL4 EG2
  } else if (trainConfig == 506){ // EMCAL clusters
    cuts.AddCutCalo("00083113","1111112067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL1 EG1
    cuts.AddCutCalo("00083113","1111122067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL2 EG1
    cuts.AddCutCalo("00083113","1111121067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL4 EG1
    cuts.AddCutCalo("00083113","1111111067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL3 EG1
  } else if (trainConfig == 507){ // EMCAL clusters
    cuts.AddCutCalo("00010113","1111102067032220000","01631031000000d0"); // -30ns, 35ns timing cut, NLtestbeam INT7
    cuts.AddCutCalo("00052113","1111102067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NLtestbeam EMC7
    cuts.AddCutCalo("00085113","1111102067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NLtestbeam EG2
    cuts.AddCutCalo("00083113","1111102067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NLtestbeam EG1
  } else if (trainConfig == 508){ // EMCAL  clusters standard cuts, INT7, NL , std TM
    cuts.AddCutCalo("00010113","111111206f032230000","01631031000000d0"); // INT7 NL 12 + TB
  } else if (trainConfig == 509){ // EMCAL  clusters standard cuts, EG2+EG1, NL , std TM
    cuts.AddCutCalo("00085113","111111206f032230000","01631031000000d0"); // EG2  NL 12 + TB
    cuts.AddCutCalo("00083113","111111206f032230000","01631031000000d0"); // EG1  NL 12 + TB

  // *********************************************************************************************************
  // 13 TeV  DMC configurations
  // *********************************************************************************************************
  } else if (trainConfig==510){ //DCAL
    cuts.AddCutCalo("00010113","3885500017032220000","01631031000000d0"); // 1000ns timing cut, no NL INT7
    cuts.AddCutCalo("00010113","3885500067032220000","01631031000000d0"); // -30ns, 35ns timing cut, no NL INT7
  } else if (trainConfig == 511){ // DCAL clusters
    cuts.AddCutCalo("00055113","3885500017032220000","01631031000000d0"); // 1000ns timing cut, no NL DMC7
    cuts.AddCutCalo("00089113","3885500017032220000","01631031000000d0"); // 1000ns timing cut, no NL DG2
    cuts.AddCutCalo("0008b113","3885500017032220000","01631031000000d0"); // 1000ns timing cut, no NL DG1
  } else if (trainConfig == 512){ // DCAL clusters
    cuts.AddCutCalo("00055113","3885500067032220000","01631031000000d0"); // -50ns, 30ns timing cut, no NL DMC7
    cuts.AddCutCalo("00089113","3885500067032220000","01631031000000d0"); // -50ns, 30ns timing cut, no NL DG2
    cuts.AddCutCalo("0008b113","3885500067032220000","01631031000000d0"); // -50ns, 30ns timing cut, no NL DG1
  } else if (trainConfig == 513){ // DCal  clusters standard cuts, INT7, NL , std TM
    cuts.AddCutCalo("00010113","388551206f032230000","01631031000000d0"); // INT7 NL 12 + TB
  } else if (trainConfig == 514){ // DCal  clusters standard cuts, EG2+EG1, NL , std TM
    cuts.AddCutCalo("00089113","388551206f032230000","01631031000000d0"); // DG2  NL 12 + TB
    cuts.AddCutCalo("0008b113","388551206f032230000","01631031000000d0"); // DG1  NL 12 + TB

  // *********************************************************************************************************
  // 13 TeV  pp Run2 - EMC configurations HM trigg
  // *********************************************************************************************************
  } else if (trainConfig == 520){ // EMCAL clusters
    cuts.AddCutCalo("00074113","1111100017032220000","01631031000000d0"); // 1000ns timing cut, no NL VOHM
    cuts.AddCutCalo("00076113","1111100017032220000","01631031000000d0"); // 1000ns timing cut, no NL VOHM with SPD
    cuts.AddCutCalo("00074113","1111100067032220000","01631031000000d0"); // -50ns, 30ns timing cut, no NL VOHM
    cuts.AddCutCalo("00076113","1111100067032220000","01631031000000d0"); // -50ns, 30ns timing cut, no NL VOHM with SPD
  } else if (trainConfig == 521){ // EMCAL clusters
    cuts.AddCutCalo("00074113","1111112067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL1 VOHM
    cuts.AddCutCalo("00074113","1111122067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL2 VOHM
    cuts.AddCutCalo("00074113","1111102067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NLtestbeam VOHM
    cuts.AddCutCalo("00074113","1111111067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL3 VOHM
    cuts.AddCutCalo("00074113","1111121067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL4 VOHM
  } else if (trainConfig == 522){ // EMCAL clusters
    cuts.AddCutCalo("00076113","1111121067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL4 VOHM with SPD
    cuts.AddCutCalo("00076113","1111111067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL3 VOHM with SPD
    cuts.AddCutCalo("00076113","1111102067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NLtestbeam VOHM with SPD
    cuts.AddCutCalo("00076113","1111122067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL2 VOHM with SPD
    cuts.AddCutCalo("00076113","1111112067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL1 VOHM with SPD

  // *********************************************************************************************************
  // 13 TeV  DMC configurations HM trigg
  // *********************************************************************************************************
  } else if (trainConfig==530){ //DCAL
    cuts.AddCutCalo("00074113","3885500017032220000","01631031000000d0"); // 1000ns timing cut, no NL VOHM
    cuts.AddCutCalo("00074113","3885500067032220000","01631031000000d0"); // -50ns, 30ns timing cut, no NL VOHM
    cuts.AddCutCalo("00076113","3885500017032220000","01631031000000d0"); // 1000ns timing cut, no NL VOHM with SPD
    cuts.AddCutCalo("00076113","3885500067032220000","01631031000000d0"); // -50ns, 30ns timing cut, no NL VOHM with SPD

  // *********************************************************************************************************
  // 13 TeV  EMC configurations Mult dependant
  // *********************************************************************************************************
  } else if (trainConfig == 540){ // EMCAL
    cuts.AddCutCalo("00110113","1111122067032220000","01631031000000d0"); // -30ns, 35ns timing cut, NL2 INT7
    cuts.AddCutCalo("01210113","1111122067032220000","01631031000000d0"); // -30ns, 35ns timing cut, NL2 INT7
    cuts.AddCutCalo("02310113","1111122067032220000","01631031000000d0"); // -30ns, 35ns timing cut, NL2 INT7
    cuts.AddCutCalo("03410113","1111122067032220000","01631031000000d0"); // -30ns, 35ns timing cut, NL2 INT7
    cuts.AddCutCalo("04510113","1111122067032220000","01631031000000d0"); // -30ns, 35ns timing cut, NL2 INT7
    cuts.AddCutCalo("05610113","1111122067032220000","01631031000000d0"); // -30ns, 35ns timing cut, NL2 INT7
    cuts.AddCutCalo("06710113","1111122067032220000","01631031000000d0"); // -30ns, 35ns timing cut, NL2 INT7
    cuts.AddCutCalo("07810113","1111122067032220000","01631031000000d0"); // -30ns, 35ns timing cut, NL2 INT7
  } else if (trainConfig == 541){ // EMCAL HM Trigger
    cuts.AddCutCalo("00174113","1111122067032220000","01631031000000d0"); // -30ns, 35ns timing cut, NL2 INT7, VOHM
    cuts.AddCutCalo("01274113","1111122067032220000","01631031000000d0"); // -30ns, 35ns timing cut, NL2 INT7, VOHM
    cuts.AddCutCalo("02374113","1111122067032220000","01631031000000d0"); // -30ns, 35ns timing cut, NL2 INT7, VOHM
    cuts.AddCutCalo("03474113","1111122067032220000","01631031000000d0"); // -30ns, 35ns timing cut, NL2 INT7, VOHM
    cuts.AddCutCalo("04574113","1111122067032220000","01631031000000d0"); // -30ns, 35ns timing cut, NL2 INT7, VOHM
    cuts.AddCutCalo("05674113","1111122067032220000","01631031000000d0"); // -30ns, 35ns timing cut, NL2 INT7, VOHM
    cuts.AddCutCalo("06774113","1111122067032220000","01631031000000d0"); // -30ns, 35ns timing cut, NL2 INT7, VOHM
    cuts.AddCutCalo("07874113","1111122067032220000","01631031000000d0"); // -30ns, 35ns timing cut, NL2 INT7, VOHM
  } else if (trainConfig == 542){ // EMCAL HM Trigger
    cuts.AddCutCalo("00176113","1111122067032220000","01631031000000d0"); // -30ns, 35ns timing cut, NL2 INT7, VOHM with SPD
    cuts.AddCutCalo("01276113","1111122067032220000","01631031000000d0"); // -30ns, 35ns timing cut, NL2 INT7, VOHM with SPD
    cuts.AddCutCalo("02376113","1111122067032220000","01631031000000d0"); // -30ns, 35ns timing cut, NL2 INT7, VOHM with SPD
    cuts.AddCutCalo("03476113","1111122067032220000","01631031000000d0"); // -30ns, 35ns timing cut, NL2 INT7, VOHM with SPD
    cuts.AddCutCalo("04576113","1111122067032220000","01631031000000d0"); // -30ns, 35ns timing cut, NL2 INT7, VOHM with SPD
    cuts.AddCutCalo("05676113","1111122067032220000","01631031000000d0"); // -30ns, 35ns timing cut, NL2 INT7, VOHM with SPD
    cuts.AddCutCalo("06776113","1111122067032220000","01631031000000d0"); // -30ns, 35ns timing cut, NL2 INT7, VOHM with SPD
    cuts.AddCutCalo("07876113","1111122067032220000","01631031000000d0"); // -30ns, 35ns timing cut, NL2 INT7, VOHM with SPD

  } else if (trainConfig == 550){  // EMCAL+DCAL clusters 13 TeV std. QA
    cuts.AddCutCalo("00010113","4117900007032220000","01631031000000d0"); // no timing cut, no NL INT7
    cuts.AddCutCalo("00010113","4117900067032220000","01631031000000d0"); // no timing cut, no NL INT7
  } else if (trainConfig == 551){  // EMCAL+DCAL clusters 13 TeV std. QA
    cuts.AddCutCalo("00010113","4117900067032220000","01631031000000d0"); // -30ns, 35ns timing cut, no NL INT7
    cuts.AddCutCalo("00010113","41179000a7032220000","01631031000000d0"); // -12.5ns, 12.5ns timing cut, no NL INT7
    cuts.AddCutCalo("00010113","411790006f032220000","01631031000000d0"); // -30ns, 35ns timing cut, no NL INT7
    cuts.AddCutCalo("00010113","41179000af032220000","01631031000000d0"); // -12.5ns, 12.5ns timing cut, no NL INT7
  } else if (trainConfig == 552){ // EMCAL+DCAL clusters standard cuts, triggers, no NL, std TM, tight timing
    cuts.AddCutCalo("00010113","41179000a7032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("0008e113","41179000a7032230000","01631031000000d0"); // EG2
    cuts.AddCutCalo("0008d113","41179000a7032230000","01631031000000d0"); // EG1
    cuts.AddCutCalo("0009b113","41179060a7032230000","01631031000000d0"); // EJ1
  } else if (trainConfig == 553){ // EMCAL+DCAL clusters standard cuts, triggers, no NL, E/p TM, tight timing
    cuts.AddCutCalo("00010113","41179000af032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("0008e113","41179000af032230000","01631031000000d0"); // EG2
    cuts.AddCutCalo("0008d113","41179000af032230000","01631031000000d0"); // EG1
    cuts.AddCutCalo("0009b113","41179060af032230000","01631031000000d0"); // EJ1
  } else if (trainConfig == 554){ // EMCAL+DCAL clusters standard cuts, triggers, no NL, std TM, -50ns, 30ns timing cut
    cuts.AddCutCalo("00010113","4117900067032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("0008e113","4117900067032230000","01631031000000d0"); // EG2
    cuts.AddCutCalo("0008d113","4117900067032230000","01631031000000d0"); // EG1
    cuts.AddCutCalo("0009b113","4117906067032230000","01631031000000d0"); // EJ1
  } else if (trainConfig == 555){ // EMCAL+DCAL clusters standard cuts, triggers, no NL, E/p TM, tight timing
    cuts.AddCutCalo("00010113","411790006f032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("0008e113","411790006f032230000","01631031000000d0"); // EG2
    cuts.AddCutCalo("0008d113","411790006f032230000","01631031000000d0"); // EG1
    cuts.AddCutCalo("0009b113","411790606f032230000","01631031000000d0"); // EJ1
  } else if (trainConfig == 556){ // EMCAL+DCAL clusters standard cuts, triggers, no NL, E/p TM, tight timing
    cuts.AddCutCalo("00010113","41179110af032220000","01631031000000d0"); // -12.5ns, 12.5ns timing cut
    cuts.AddCutCalo("00010113","41179120af032220000","01631031000000d0"); // -12.5ns, 12.5ns timing cut
    cuts.AddCutCalo("00010113","41179210af032220000","01631031000000d0"); // -12.5ns, 12.5ns timing cut
    cuts.AddCutCalo("00010113","41179220af032220000","01631031000000d0"); // -12.5ns, 12.5ns timing cut
  } else if (trainConfig == 570){ //EMC + DCal HM trigger
    cuts.AddCutCalo("00010113","4117900067032230000","01631031000000d0"); // -50ns, 30ns timing cut, MB trigger
    cuts.AddCutCalo("00010113","4117900007032230000","01631031000000d0"); // no timing cut, MB trigger
    cuts.AddCutCalo("00074113","4117900067032230000","01631031000000d0"); // -50ns, 30ns timing cut, no NL VOHM
    cuts.AddCutCalo("00076113","4117900067032230000","01631031000000d0"); // -50ns, 30ns timing cut, no NL VOHM with SPD
  // *********************************************************************************************************
  // 5 TeV 2015 pp Run2 - DMC configurations
  // *********************************************************************************************************
  //                 LHC17pq
  } else if (trainConfig == 600){
    cuts.AddCutCalo("00010113","3885512087032220000","01631031000000d0"); // std
    cuts.AddCutCalo("00010113","3885512017032220000","01631031000000d0"); // QA
  } else if (trainConfig == 601){ // DCAL clusters pp 5.02 TeV (Triggered)
    cuts.AddCutCalo("00055113","3885512087032220000","01631031000000d0"); // EMC7
    // Changed BGEvents to 20(80) for variations
  } else if (trainConfig == 602){ // timing Cut variation  std -20+50ns
    cuts.AddCutCalo("00010113","3885512017032220000","01631031000000d0"); //     -1000  +1000 ns
    cuts.AddCutCalo("00010113","3885512077032220000","01631031000000d0"); //     -30    +30   ns
    cuts.AddCutCalo("00010113","3885512097032220000","01631031000000d0"); //     -20    +25   ns
    cuts.AddCutCalo("00010113","38855120a7032220000","01631031000000d0"); //     -12.5  +13   ns
  } else if (trainConfig == 603){ // track matching variation
    cuts.AddCutCalo("00010113","3885512080032220000","01631031000000d0"); //
    cuts.AddCutCalo("00010113","3885512081032220000","01631031000000d0"); //
    cuts.AddCutCalo("00010113","3885512086032220000","01631031000000d0"); //
  } else if (trainConfig == 604){ // opening angle variation
    cuts.AddCutCalo("00010113","3885512087032220000","0163103100000000"); //
    cuts.AddCutCalo("00010113","3885512087032220000","0163103100000010"); //
    cuts.AddCutCalo("00010113","3885512087032220000","0163103100000030"); //
    cuts.AddCutCalo("00010113","3885512087032220000","0163103100000040"); //
    cuts.AddCutCalo("00010113","3885512087032220000","0163103100000050"); //
    cuts.AddCutCalo("00010113","3885512087032220000","0163103100000080"); //
  } else if (trainConfig == 605){ // min energy variation std 0.7 GeV/c
    cuts.AddCutCalo("00010113","3885512087002220000","01631031000000d0"); //     0.1 GeV/c
    cuts.AddCutCalo("00010113","3885512087012220000","01631031000000d0"); //     0.5 GeV/c
    cuts.AddCutCalo("00010113","3885512087022220000","01631031000000d0"); //     0.6 GeV/c
    cuts.AddCutCalo("00010113","3885512087042220000","01631031000000d0"); //     0.8 GeV/c
    cuts.AddCutCalo("00010113","3885512087052220000","01631031000000d0"); //     0.9 GeV/c
  } else if (trainConfig == 606){ // min nCells & M02 variation
    // std: min nCells = 1; M02 max=0.7, min=0.1
    cuts.AddCutCalo("00010113","3885512087031220000","01631031000000d0"); //   min nCells = 1
    cuts.AddCutCalo("00010113","3885512087033220000","01631031000000d0"); //   min nCells = 3
    cuts.AddCutCalo("00010113","3885512087032210000","01631031000000d0"); //   max M02    = 1
    cuts.AddCutCalo("00010113","3885512087032240000","01631031000000d0"); //   max M02    = 0.4
  } else if (trainConfig == 607){ // NonLin variation
    cuts.AddCutCalo("00010113","3885500087032220000","01631031000000d0"); // no NL
    cuts.AddCutCalo("00010113","3885511087032220000","01631031000000d0"); // PCM-DCal kSDM
    cuts.AddCutCalo("00010113","3885512087032220000","01631031000000d0"); // DCal kSDM
    cuts.AddCutCalo("00010113","3885521087032220000","01631031000000d0"); // PCM-DCal DExp/DPow
    cuts.AddCutCalo("00010113","3885522087032220000","01631031000000d0"); // DCal DExp/DPow
  } else if (trainConfig == 608){ // LHC15n
    cuts.AddCutCalo("00010113","3885522087032220000","01631031000000d0"); // std
    cuts.AddCutCalo("00010113","3885500017032220000","01631031000000d0"); // QA
  } else if (trainConfig == 610){ // no NL tight timing
    cuts.AddCutCalo("00010113","38855000a7032230000","01631031000000d0"); // -12.5 to 13 ns timing and M02 0.1-0.5
    cuts.AddCutCalo("00010113","38855060a7032230000","01631031000000d0"); // -12.5 to 13 ns timing and M02 0.1-0.5 TB NL
  } else if (trainConfig == 611){ // new NL variations
    cuts.AddCutCalo("00010113","38855610a7032230000","01631031000000d0"); // -12.5 to 13 ns timing and M02 0.1-0.5
    cuts.AddCutCalo("00010113","38855630a7032230000","01631031000000d0"); // -12.5 to 13 ns timing and M02 0.1-0.5
  } else if (trainConfig == 612){ // new NL variations
    cuts.AddCutCalo("00010113","38855620a7032230000","01631031000000d0"); // -12.5 to 13 ns timing and M02 0.1-0.5
    cuts.AddCutCalo("00010113","38855640a7032230000","01631031000000d0"); // -12.5 to 13 ns timing and M02 0.1-0.5

  } else if (trainConfig == 660){ // INT7 DCAL standard cut but with E/p TM veto
    cuts.AddCutCalo("00010113","3885511087041220000","01631031000000d0"); // std INT7
    cuts.AddCutCalo("00010113","388551108c041220000","01631031000000d0"); // fEOverPMax = 9e9
    cuts.AddCutCalo("00010113","388551108d041220000","01631031000000d0"); // fEOverPMax = 3.0
    cuts.AddCutCalo("00010113","388551108e041220000","01631031000000d0"); // fEOverPMax = 2.0
  } else if (trainConfig == 661){ // INT7 DCAL standard cut but with E/p TM veto
    cuts.AddCutCalo("00010113","388551108f041220000","01631031000000d0"); // fEOverPMax = 1.75
    cuts.AddCutCalo("00010113","388551108g041220000","01631031000000d0"); // fEOverPMax = 1.5
    cuts.AddCutCalo("00010113","388551108h041220000","01631031000000d0"); // fEOverPMax = 1.25
  } else if (trainConfig == 662){ // DMC7 DCAL standard cut but with E/p TM veto
    cuts.AddCutCalo("00055113","3885511087041220000","01631031000000d0"); // std DMC7
    cuts.AddCutCalo("00055113","388551108c041220000","01631031000000d0"); // fEOverPMax = 9e9
    cuts.AddCutCalo("00055113","388551108d041220000","01631031000000d0"); // fEOverPMax = 3.0
    cuts.AddCutCalo("00055113","388551108e041220000","01631031000000d0"); // fEOverPMax = 2.0
  } else if (trainConfig == 663){ // DMC7 DCAL standard cut but with E/p TM veto
    cuts.AddCutCalo("00055113","388551108f041220000","01631031000000d0"); // fEOverPMax = 1.75
    cuts.AddCutCalo("00055113","388551108g041220000","01631031000000d0"); // fEOverPMax = 1.5
    cuts.AddCutCalo("00055113","388551108h041220000","01631031000000d0"); // fEOverPMax = 1.25
  } else if (trainConfig == 664){ // DG1 DCAL standard cut but with E/p TM veto
    cuts.AddCutCalo("0008b113","3885511087041220000","01631031000000d0"); // std DG1
    cuts.AddCutCalo("0008b113","388551108c041220000","01631031000000d0"); // fEOverPMax = 9e9
    cuts.AddCutCalo("0008b113","388551108d041220000","01631031000000d0"); // fEOverPMax = 3.0
    cuts.AddCutCalo("0008b113","388551108e041220000","01631031000000d0"); // fEOverPMax = 2.0
  } else if (trainConfig == 665){ // DG1 DCAL standard cut but with E/p TM veto
    cuts.AddCutCalo("0008b113","388551108f041220000","01631031000000d0"); // fEOverPMax = 1.75
    cuts.AddCutCalo("0008b113","388551108g041220000","01631031000000d0"); // fEOverPMax = 1.5
    cuts.AddCutCalo("0008b113","388551108h041220000","01631031000000d0"); // fEOverPMax = 1.25
  } else if (trainConfig == 666){ // DG2 DCAL standard cut but with E/p TM veto
    cuts.AddCutCalo("00089113","3885511087041220000","01631031000000d0"); // std DG2
    cuts.AddCutCalo("00089113","388551108c041220000","01631031000000d0"); // fEOverPMax = 9e9
    cuts.AddCutCalo("00089113","388551108d041220000","01631031000000d0"); // fEOverPMax = 3.0
    cuts.AddCutCalo("00089113","388551108e041220000","01631031000000d0"); // fEOverPMax = 2.0
  } else if (trainConfig == 667){ // DG2 DCAL standard cut but with E/p TM veto
    cuts.AddCutCalo("00089113","388551108f041220000","01631031000000d0"); // fEOverPMax = 1.75
    cuts.AddCutCalo("00089113","388551108g041220000","01631031000000d0"); // fEOverPMax = 1.5
    cuts.AddCutCalo("00089113","388551108h041220000","01631031000000d0"); // fEOverPMax = 1.25

  } else if (trainConfig == 680){ // DCAL standard for CALO+CALOFAST readout triggers - NO TM
    cuts.AddCutCalo("000a0113","3885511080041220000","01631031000000d0"); // std INT7
    cuts.AddCutCalo("000a6113","3885511080041220000","01631031000000d0"); // std DMC7
    cuts.AddCutCalo("000a7113","3885511080041220000","01631031000000d0"); // std DG2
    cuts.AddCutCalo("000a8113","3885511080041220000","01631031000000d0"); // std DG1
  } else if (trainConfig == 681){ // DCAL standard for standard readout triggers - NO TM
    cuts.AddCutCalo("00010113","3885511080041220000","01631031000000d0"); // std INT7
    cuts.AddCutCalo("00055113","3885511080041220000","01631031000000d0"); // std DMC7
    cuts.AddCutCalo("00089113","3885511080041220000","01631031000000d0"); // std DG2
    cuts.AddCutCalo("0008b113","3885511080041220000","01631031000000d0"); // std DG1
  // *********************************************************************************************************
  // 5 TeV 2015 pp Run2 - PHOS configurations
  // *********************************************************************************************************
  } else if (trainConfig == 700){ // PHOS clusters with larger acceptance
    cuts.AddCutCalo("00010113","2446600040012300000","0163103100000010"); // INT7
    cuts.AddCutCalo("00062113","2446600040012300000","0163103100000010"); // PHI7
  } else if (trainConfig == 701){ // Default cut, No TM
    cuts.AddCutCalo("00010113","2446651040012300000","0163103100000010"); // INT7
    cuts.AddCutCalo("00062113","2446651040012300000","0163103100000010"); // PHI7
  } else if (trainConfig == 702){ // Default cut, with TM
    cuts.AddCutCalo("00010113","2446651044012300000","0163103100000010"); // INT7
    cuts.AddCutCalo("00062113","2446651044012300000","0163103100000010"); // PHI7
  } else if( trainConfig == 703){ // DEFAULT 2018 oct 31 no NL
    cuts.AddCutCalo("00010113","2446600044012300000","0163103100000010"); //
  } else if( trainConfig == 704){ // DEFAULT 2018 oct 31 NL variations
    cuts.AddCutCalo("00010113","2446600044012300000","0163103100000010"); //
    cuts.AddCutCalo("00010113","2446601044012300000","0163103100000010"); // PHOS people NL
    cuts.AddCutCalo("00010113","2446651044012300000","0163103100000010"); // INT7
    cuts.AddCutCalo("00010113","2446652044012300000","0163103100000010"); // PHOS calo NL
  } else if( trainConfig == 705){ // Timing cut efficiency studies
    cuts.AddCutCalo("00010113","244665100a012200000","0163103100000010"); //
    cuts.AddCutCalo("000ap113","2446651000012200000","0163103100000010"); // PHI7
  } else if( trainConfig == 706){ // DEFAULT 2019 july 18
    cuts.AddCutCalo("00010113","244665107a012200000","0163103100000010"); //
    cuts.AddCutCalo("00010113","2446651070012200000","0163103100000010"); //
    cuts.AddCutCalo("000ap113","2446651070012200000","0163103100000010"); // PHI7
  } else if( trainConfig == 707){ // No non-lin corr, use with Run2Tune / Run2TuneMC
    cuts.AddCutCalo("00010113","244660007a012200000","0163103100000010"); //
    cuts.AddCutCalo("00010113","2446600070012200000","0163103100000010"); //
    cuts.AddCutCalo("000ap113","2446600070012200000","0163103100000010"); // PHI7
  } else if( trainConfig == 708){ // Sphericity PHOS PHOS
    cuts.AddCutCalo("00010113","24466510ga012200000","0163103100000010"); //
    cuts.AddCutCalo("h0510113","24466510ga012200000","0163103100000010"); //
    cuts.AddCutCalo("h5a10113","24466510ga012200000","0163103100000010"); //
    cuts.AddCutCalo("h0a10113","24466510ga012200000","0163103100000010"); //
    cuts.AddCutCalo("h0310113","24466510ga012200000","0163103100000010"); //
    cuts.AddCutCalo("h7a10113","24466510ga012200000","0163103100000010"); //
  } else if( trainConfig == 709){ // V0M multiplicity cuts PHOS
    cuts.AddCutCalo("m0110113","24466510ga012200000","0163103100000010"); // 0-1%
    cuts.AddCutCalo("m1510113","24466510ga012200000","0163103100000010"); // 1-5%
    cuts.AddCutCalo("m5k10113","24466510ga012200000","0163103100000010"); // 5-20%
    cuts.AddCutCalo("n2410113","24466510ga012200000","0163103100000010"); // 20-40%
    cuts.AddCutCalo("n4710113","24466510ga012200000","0163103100000010"); // 40-70%
    cuts.AddCutCalo("n7a10113","24466510ga012200000","0163103100000010"); // 70-100%
  } else if( trainConfig == 710){ // SPD multiplicity cuts PHOS
    cuts.AddCutCalo("o0110113","24466510ga012200000","0163103100000010"); // 0-1%
    cuts.AddCutCalo("o0210113","24466510ga012200000","0163103100000010"); // 0-2%
    cuts.AddCutCalo("o0510113","24466510ga012200000","0163103100000010"); // 0-5%
    cuts.AddCutCalo("o5k10113","24466510ga012200000","0163103100000010"); // 5-20%
    cuts.AddCutCalo("p2610113","24466510ga012200000","0163103100000010"); // 20-60%
    cuts.AddCutCalo("p6a10113","24466510ga012200000","0163103100000010"); // 60-100%
  } else if( trainConfig == 711){ // New standard cut August 14 2019
    cuts.AddCutCalo("00010113","24466510ga012200000","0163103100000010"); //
    cuts.AddCutCalo("00010113","24466510g0012200000","0163103100000010"); //
    cuts.AddCutCalo("000ap113","24466510g0012200000","0163103100000010"); // PHI7
  } else if( trainConfig == 712){ // Timing cut efficiency studies
    cuts.AddCutCalo("00010113","244665107a012200000","0163103100000010"); // std timing +-30ns
    cuts.AddCutCalo("00010113","24466510aa012200000","0163103100000010"); // timing -12.5 + 13 ns
    cuts.AddCutCalo("00010113","24466510ga012200000","0163103100000010"); // std timing +-30ns    ,w MC Timing cut eff. corr.
    cuts.AddCutCalo("00010113","24466510ia012200000","0163103100000010"); // timing -12.5 + 13 ns ,w MC Timing cut eff. corr.
  } else if( trainConfig == 713){ // V0M multiplicity cuts PHOS -PHI7 trigger
    cuts.AddCutCalo("m01ap113","24466510g0012200000","0163103100000010"); // 0-1%
    cuts.AddCutCalo("m15ap113","24466510g0012200000","0163103100000010"); // 1-5%
    cuts.AddCutCalo("m5kap113","24466510g0012200000","0163103100000010"); // 5-20%
    cuts.AddCutCalo("n24ap113","24466510g0012200000","0163103100000010"); // 20-40%
    cuts.AddCutCalo("n47ap113","24466510g0012200000","0163103100000010"); // 40-70%
    cuts.AddCutCalo("n7aap113","24466510g0012200000","0163103100000010"); // 70-100%
  } else if( trainConfig == 714){ // V0M multiplicity cuts PHOS -PHI7 trigger
    cuts.AddCutCalo("m0110113","24466510g0012200000","0163103100000010"); // 0-1%
    cuts.AddCutCalo("m1510113","24466510g0012200000","0163103100000010"); // 1-5%
    cuts.AddCutCalo("m5k10113","24466510g0012200000","0163103100000010"); // 5-20%
    cuts.AddCutCalo("n2410113","24466510g0012200000","0163103100000010"); // 20-40%
    cuts.AddCutCalo("n4710113","24466510g0012200000","0163103100000010"); // 40-70%
    cuts.AddCutCalo("n7a10113","24466510g0012200000","0163103100000010"); // 70-100%

  // Variations for systematics
  } else if ( trainConfig == 730) { // NL variations (standard: 42 PHOS ML)
    cuts.AddCutCalo("00010113","24466530ga01cc00000","0163103100000010");
    cuts.AddCutCalo("00010113","24466540ga01cc00000","0163103100000010");
    cuts.AddCutCalo("00010113","24466010ga01cc00000","0163103100000010"); // PHOS group NL
    cuts.AddCutCalo("00010113","24466000ga01cc00000","0163103100000010"); // PHOS group NL
  } else if ( trainConfig == 731) { // distance to bad channel variations (standard: no cut)
    cuts.AddCutCalo("00010113","24466531ga01cc00000","0163103100000010"); // dist. to bad channel = 1
    cuts.AddCutCalo("00010113","24466532ga01cc00000","0163103100000010"); // dist. to bad channel = 2
    cuts.AddCutCalo("00010113","24466533ga01cc00000","0163103100000010"); // dist. to bad channel = 3
  } else if ( trainConfig == 732) {  // timing variations - bunch spacing: 100ns (standard: 50ns) 
                                    // min. number of cells per cluster variations (standard: 2)
    cuts.AddCutCalo("00010113","24466533ha01cc00000","0163103100000010"); // 30ns
    cuts.AddCutCalo("00010113","24466533ia01cc00000","0163103100000010"); // -20ns/25ns
    cuts.AddCutCalo("00010113","24466530ga01dc00000","0163103100000010"); // nCells > 3 above 1 GeV
  } else if ( trainConfig == 733) { // track matching variations
    cuts.AddCutCalo("00010113","24466533g001cc00000","0163103100000010"); // without track matching
    cuts.AddCutCalo("00010113","24466533g501cc00000","0163103100000010"); // tm variation
    cuts.AddCutCalo("00010113","24466533g601cc00000","0163103100000010"); // tm variation
  } else if ( trainConfig == 734) { // min. cluster energy variations (standard: 0.5 MeV)
    cuts.AddCutCalo("00010113","24466530ga07cc00000","0163103100000010"); // 0.2 MeV
    cuts.AddCutCalo("00010113","24466530ga09cc00000","0163103100000010"); // 0.1 MeV
    cuts.AddCutCalo("00010113","24466530ga08cc00000","0163103100000010"); // 0.4 MeV
    cuts.AddCutCalo("00010113","24466530ga02cc00000","0163103100000010"); // 0.5 MeV
  } else if ( trainConfig == 735) { 
    cuts.AddCutCalo("00010113","24466530ga01cb00000","0163103100000010"); // M02 min 0.002, E > 1 GeV
    cuts.AddCutCalo("00010113","24466530ga01cd00000","0163103100000010"); // M02 min 0.2, E > 1 GeV
    cuts.AddCutCalo("00010113","24466530ga01cc00010","0163103100000010"); // dispersion < 2
    cuts.AddCutCalo("00010113","24466530ga01cc00030","0163103100000010"); // dispersion < 2
  } else if ( trainConfig == 736){ // MESON
    cuts.AddCutCalo("00010113","24466530ga01cc00000","0163303100000010"); // rapidity variation  y<0.5
    cuts.AddCutCalo("00010113","24466530ga01cc00000","0163803100000010"); // rapidity variation  y<0.5
    cuts.AddCutCalo("00010113","24466530ga01cc00000","0163106100000010"); // alpha meson variation 1 0<alpha<0.8
    cuts.AddCutCalo("00010113","24466530ga01cc00000","0163105100000010"); // alpha meson variation 2 0<alpha<0.75
  } else if( trainConfig == 737){ // fourth set of variations
    cuts.AddCutCalo("00010113","24466530ga01cc00300","0163103100000000"); // conv rec 0.03
    cuts.AddCutCalo("00010113","24466530ga01cc00200","0163103100000030"); // conv rec 0.025
    cuts.AddCutCalo("00010113","24466530ga01cc00400","0163103100000010"); // conv rec 0.035
    cuts.AddCutCalo("00010113","24466530ga01cc00000","0263303100000010"); // BG

  // Variations for systematics wo TM
  } else if ( trainConfig == 740) { // NL variations (standard: 42 PHOS ML)
    cuts.AddCutCalo("00010113","24466530g001cc00000","0163103100000010");
    cuts.AddCutCalo("00010113","24466540g001cc00000","0163103100000010");
    cuts.AddCutCalo("00010113","24466010g001cc00000","0163103100000010"); // PHOS group NL
    cuts.AddCutCalo("00010113","24466000g001cc00000","0163103100000010"); // PHOS group NL
  } else if ( trainConfig == 741) { // distance to bad channel variations (standard: no cut)
    cuts.AddCutCalo("00010113","24466531g001cc00000","0163103100000010"); // dist. to bad channel = 1
    cuts.AddCutCalo("00010113","24466532g001cc00000","0163103100000010"); // dist. to bad channel = 2
    cuts.AddCutCalo("00010113","24466533g001cc00000","0163103100000010"); // dist. to bad channel = 3
  } else if ( trainConfig == 742) {  // timing variations - bunch spacing: 100ns (standard: 50ns) 
                                    // min. number of cells per cluster variations (standard: 2)
    cuts.AddCutCalo("00010113","24466533h001cc00000","0163103100000010"); // 30ns
    cuts.AddCutCalo("00010113","24466533i001cc00000","0163103100000010"); // -20ns/25ns
    cuts.AddCutCalo("00010113","24466530g001dc00000","0163103100000010"); // nCells > 3 above 1 GeV
  } else if ( trainConfig == 744) { // min. cluster energy variations (standard: 0.5 MeV)
    cuts.AddCutCalo("00010113","24466530g007cc00000","0163103100000010"); // 0.2 MeV
    cuts.AddCutCalo("00010113","24466530g009cc00000","0163103100000010"); // 0.1 MeV
    cuts.AddCutCalo("00010113","24466530g008cc00000","0163103100000010"); // 0.4 MeV
    cuts.AddCutCalo("00010113","24466530g002cc00000","0163103100000010"); // 0.5 MeV
  } else if ( trainConfig == 745) { 
    cuts.AddCutCalo("00010113","24466530g001cb00000","0163103100000010"); // M02 min 0.002, E > 1 GeV
    cuts.AddCutCalo("00010113","24466530g001cd00000","0163103100000010"); // M02 min 0.2, E > 1 GeV
    cuts.AddCutCalo("00010113","24466530g001cc00010","0163103100000010"); // dispersion < 2
    cuts.AddCutCalo("00010113","24466530g001cc00030","0163103100000010"); // dispersion < 2
  } else if ( trainConfig == 746){ // MESON
    cuts.AddCutCalo("00010113","24466530g001cc00000","0163303100000010"); // rapidity variation  y<0.5
    cuts.AddCutCalo("00010113","24466530g001cc00000","0163803100000010"); // rapidity variation  y<0.5
    cuts.AddCutCalo("00010113","24466530g001cc00000","0163106100000010"); // alpha meson variation 1 0<alpha<0.8
    cuts.AddCutCalo("00010113","24466530g001cc00000","0163105100000010"); // alpha meson variation 2 0<alpha<0.75
  } else if( trainConfig == 747){ // fourth set of variations
    cuts.AddCutCalo("00010113","24466530g001cc00300","0163103100000000"); // conv rec 0.03
    cuts.AddCutCalo("00010113","24466530g001cc00200","0163103100000030"); // conv rec 0.025
    cuts.AddCutCalo("00010113","24466530g001cc00400","0163103100000010"); // conv rec 0.035
    cuts.AddCutCalo("00010113","24466530g001cc00000","0263303100000010"); // BG

  // Variations for systematics
  } else if ( trainConfig == 750) { // NL variations (standard: 42 PHOS ML)
    cuts.AddCutCalo("00062113","24466530g001cc00000","0163103100000010");
    cuts.AddCutCalo("00062113","24466540g001cc00000","0163103100000010");
    cuts.AddCutCalo("00062113","24466010g001cc00000","0163103100000010"); // PHOS group NL
    cuts.AddCutCalo("00062113","24466000g001cc00000","0163103100000010"); // PHOS group NL
  } else if ( trainConfig == 751) { // distance to bad channel variations (standard: no cut)
    cuts.AddCutCalo("00062113","24466531g001cc00000","0163103100000010"); // dist. to bad channel = 1
    cuts.AddCutCalo("00062113","24466532g001cc00000","0163103100000010"); // dist. to bad channel = 2
    cuts.AddCutCalo("00062113","24466533g001cc00000","0163103100000010"); // dist. to bad channel = 3
  } else if ( trainConfig == 752) {  // timing variations - bunch spacing: 100ns (standard: 50ns) 
                                    // min. number of cells per cluster variations (standard: 2)
    cuts.AddCutCalo("00062113","244665330001cc00000","0163103100000010"); // 30ns
    cuts.AddCutCalo("00062113","24466533i001cc00000","0163103100000010"); // 30ns
    cuts.AddCutCalo("00062113","24466533h001cc00000","0163103100000010"); // -20ns/25ns
    cuts.AddCutCalo("00062113","24466530g001dc00000","0163103100000010"); // nCells > 3 above 1 GeV
  } else if ( trainConfig == 754) { // min. cluster energy variations (standard: 0.5 MeV)
    cuts.AddCutCalo("00062113","24466530g007cc00000","0163103100000010"); // 0.2 MeV
    cuts.AddCutCalo("00062113","24466530g009cc00000","0163103100000010"); // 0.1 MeV
    cuts.AddCutCalo("00062113","24466530g008cc00000","0163103100000010"); // 0.4 MeV
    cuts.AddCutCalo("00062113","24466530g002cc00000","0163103100000010"); // 0.5 MeV
  } else if ( trainConfig == 755) { 
    cuts.AddCutCalo("00062113","24466530g001cb00000","0163103100000010"); // M02 min 0.002, E > 1 GeV
    cuts.AddCutCalo("00062113","24466530g001cd00000","0163103100000010"); // M02 min 0.2, E > 1 GeV
    cuts.AddCutCalo("00062113","24466530g001cc00010","0163103100000010"); // dispersion < 2
    cuts.AddCutCalo("00062113","24466530g001cc00030","0163103100000010"); // dispersion < 2
  } else if ( trainConfig == 756){ // MESON
    cuts.AddCutCalo("00062113","24466530g001cc00000","0163303100000010"); // rapidity variation  y<0.5
    cuts.AddCutCalo("00062113","24466530g001cc00000","0163803100000010"); // rapidity variation  y<0.5
    cuts.AddCutCalo("00062113","24466530g001cc00000","0163106100000010"); // alpha meson variation 1 0<alpha<0.8
    cuts.AddCutCalo("00062113","24466530g001cc00000","0163105100000010"); // alpha meson variation 2 0<alpha<0.75
  } else if( trainConfig == 757){ // fourth set of variations
    cuts.AddCutCalo("00062113","24466530g001cc00300","0163103100000000"); // conv rec 0.03
    cuts.AddCutCalo("00062113","24466530g001cc00200","0163103100000030"); // conv rec 0.025
    cuts.AddCutCalo("00062113","24466530g001cc00400","0163103100000010"); // conv rec 0.035
    cuts.AddCutCalo("00062113","24466530g001cc00000","0263303100000010"); // BG
    
  // *********************************************************************************************************
  // 13 TeV 2015 pp Run2 - PHOS configurations
  // *********************************************************************************************************
  } else if (trainConfig == 800){ // PHOS INT7, 300MeV
    cuts.AddCutCalo("00010113","24466190pa01cc00000","0163103100000010"); //Int7 no Trigger
  } else if (trainConfig == 801){ // PHOS PHI7, 300MeV
    cuts.AddCutCalo("00062113","24466190pa01cc00000","0163103100000010"); //PHI7
  } else if (trainConfig == 802){ // PHOS clusters with larger acceptance w/ TM NCells 3
    cuts.AddCutCalo("00010113","24466190pa012300000","0163103100000010"); // INT7
    cuts.AddCutCalo("00062113","24466190pa012300000","0163103100000010"); // PHI7
  } else if (trainConfig == 803){ // PHOS clusters with larger acceptance w/ TM NCells 2
    cuts.AddCutCalo("00010113","24466190pa01cc00000","0163103100000010"); // INT7
    cuts.AddCutCalo("00062113","24466190pa01cc00000","0163103100000010"); // PHI7
  } else if (trainConfig == 804){ // QA
    cuts.AddCutCalo("00010113","24466190n0012300000","0163103100000010"); // INT7 NCells 3
    cuts.AddCutCalo("00010113","24466190pa012300000","0163103100000010"); // INT7 w/ TM NCells 3
  } else if (trainConfig == 805){ // NCell Cut Variations
    cuts.AddCutCalo("00010113","24466190pa012200000","0163103100000010"); // INT7
    cuts.AddCutCalo("00010113","24466190pa012c00000","0163103100000010"); // INT7
    cuts.AddCutCalo("00010113","24466190pa01c200000","0163103100000010"); // INT7
    cuts.AddCutCalo("00010113","24466190pa01cc00000","0163103100000010"); // INT7
    cuts.AddCutCalo("00010113","24466190pa01d200000","0163103100000010"); // INT7
    cuts.AddCutCalo("00010113","24466190pa01dc00000","0163103100000010"); // INT7
  } else if (trainConfig ==806){//Comparing CellQA Config from GammaConv
    cuts.AddCutCalo("00010113","24466190pa01cc00000","0163103100000010"); // INT7
  } else if (trainConfig ==807){//Non Lin Studies
    //cuts.AddCutCalo("00010113","24466000pa01cc00000","0163103100000010"); // No Nonlin MB
    //cuts.AddCutCalo("00062113","24466000pa01cc00000","0163103100000010"); // No Nonlin Triggered
    //-
    //cuts.AddCutCalo("00010113","24466510pa01cc00000","0163103100000010"); // 51 Nonlin MB
    //cuts.AddCutCalo("00062113","24466510pa01cc00000","0163103100000010"); // 51 Nonlin Triggered
    //-
    cuts.AddCutCalo("00010113","24466110pa01cc00000","0163103100000010"); // INT7 //case 11=> FunctionNL_kSDM MB PCMPHOS
    cuts.AddCutCalo("00062113","24466110pa01cc00000","0163103100000010"); // INT7 //case 11=> FunctionNL_kSDM Triggered PCMPHOS
    //-
    cuts.AddCutCalo("00010113","24466120pa01cc00000","0163103100000010"); // INT7 //case 12=> FunctionNL_kSDM MB PHOSPHOS
    cuts.AddCutCalo("00062113","24466120pa01cc00000","0163103100000010"); // INT7 //case 12=> FunctionNL_kSDM Triggered PHOSPHOS
    //
    cuts.AddCutCalo("00010113","24466190pa01cc00000","0163103100000010"); // INT7 //case 19=> FunctionNL_kSDM MB PCMPHOS, Tuned with PHOSPHOS
    cuts.AddCutCalo("00062113","24466190pa01cc00000","0163103100000010"); // INT7 //case 19=> FunctionNL_kSDM Triggered PCMPHOS, Tuned with PHOSPHOS
    //-
    //cuts.AddCutCalo("00010113","24466210pa01cc00000","0163103100000010"); // INT7 //case 21=> unctionNL_DPOW MB
    //cuts.AddCutCalo("00062113","24466210pa01cc00000","0163103100000010"); // INT7 //case 21=> unctionNL_DPOW Triggered
  } else if (trainConfig ==808){//PHOS MB and PHI7, 100MeV
    cuts.AddCutCalo("00010113","24466190pa09cc00000","0163103100000010"); //no Trigger
    cuts.AddCutCalo("00062113","24466190pa09cc00000","0163103100000010"); //PHI7
  } else if (trainConfig ==810){//PHOS Sphericity Check
    cuts.AddCutCalo("h0510113","24466190pa01cc00000","0163103100000010"); //  0.    - 0.5
    cuts.AddCutCalo("h5a10113","24466190pa01cc00000","0163103100000010"); //  0.5    - 1.
    cuts.AddCutCalo("h0a10113","24466190pa01cc00000","0163103100000010"); //  0.    - 1.
  } else if (trainConfig ==811){//PHOS Mult Check
    cuts.AddCutCalo("n0110113","24466190pa01cc00000","0163103100000010"); // INT7 0-10%
    cuts.AddCutCalo("n1210113","24466190pa01cc00000","0163103100000010"); // INT7 10-20%
    cuts.AddCutCalo("n2510113","24466190pa01cc00000","0163103100000010"); // INT7 20-50%
    cuts.AddCutCalo("n5a10113","24466190pa01cc00000","0163103100000010"); // INT7 50-100%
    cuts.AddCutCalo("m0110113","24466190pa01cc00000","0163103100000010"); // INT7
    cuts.AddCutCalo("m1510113","24466190pa01cc00000","0163103100000010"); // INT7
    cuts.AddCutCalo("m5a10113","24466190pa01cc00000","0163103100000010"); // INT7
  } else if (trainConfig ==812){//PHOS Triggers Timing Cut 0
    cuts.AddCutCalo("00010113","244661900a09cc00000","0163103100000010"); //no Trigger
    cuts.AddCutCalo("00062113","244661900a09cc00000","0163103100000010"); //PHI7
  } else if (trainConfig ==813){//PHOS Triggers Timing Cut Studies
    cuts.AddCutCalo("00010113","24466190ga01cc00000","0163103100000010"); //no Trigger, Mike's Timing
    cuts.AddCutCalo("00062113","24466190ga01cc00000","0163103100000010"); //PHI7, Mike's Timing
  } else if (trainConfig ==814){ //PHOS Triggers Timing Cut Studies without throwing out clusters
    cuts.AddCutCalo("00010113","244661907a01cc00000","0163103100000010"); //no Trigger, Mike's Timing
    cuts.AddCutCalo("00062113","244661907a01cc00000","0163103100000010"); //PHI7, Mike's Timing
  } else if (trainConfig == 815){ // NCell Cut Variations, without E>1GeV
    cuts.AddCutCalo("00010113","24466190pa012200000","0163103100000010"); // INT7 NCells 3
    cuts.AddCutCalo("00010113","24466190pa01cc00000","0163103100000010"); // INT7 NCells 2
  }  else if (trainConfig == 840){ // PHOS INT7, 100MeV, with Timing Efficiency
    cuts.AddCutCalo("00010113","24466000pa09cc00000","0163103100000010"); //Int7 no Trigger
  }  else if (trainConfig == 841){ // PHOS INT7, 100MeV, no Timing Efficiency
    cuts.AddCutCalo("00010113","244660000a09cc00000","0163103100000010"); //Int7 no Trigger
  } else if( trainConfig == 870){ // PHOS HM trigger
    cuts.AddCutCalo("00010113","2446600044012300000","0163103100000010"); // -50ns, 30ns timing cut, MB trigger
    cuts.AddCutCalo("00010113","2446600004012300000","0163103100000010"); // no timing, MB trigger
    cuts.AddCutCalo("00074113","2446600044012300000","0163103100000010"); // -50ns, 30ns timing cut, no NL VOHM
    cuts.AddCutCalo("00076113","2446600044012300000","0163103100000010"); // -50ns, 30ns timing cut, no NL VOHM with SPD

  // *********************************************************************************************************
  // 5 TeV 2017 pp - Jet configurations
  // *********************************************************************************************************

  } else if (trainConfig == 900){
    cuts.AddCutCalo("00010113","1111111077032230000","21631031000000d0"); // std
  } else if (trainConfig == 901){
    cuts.AddCutCalo("00010113","1111111077032230000","2i631031000000d0"); // MaxPtSector mixing
  } else if (trainConfig == 902){
    cuts.AddCutCalo("00010113","1111111077032230000","2j631031000000d0"); // JetSector mixing
  } else if (trainConfig == 903){
    cuts.AddCutCalo("00010113","111111107l032230000","2k631031000000d0"); // Jet mixing
  } else if (trainConfig == 904){
    cuts.AddCutCalo("00010113","111111107l032230000","2l631031000000d0"); // JetRotation mixing
  } else if (trainConfig == 905){
    cuts.AddCutCalo("00010113","111111107l032230000","2m631031000000d0"); // Jet mixing with Jet pt
  } else if (trainConfig == 906){
    cuts.AddCutCalo("00010113","111111107l032230000","2n631031000000d0"); // JetRotation mixing with Jet pt
  } else if (trainConfig == 907){
    cuts.AddCutCalo("00010113","111111107l032230000","21631031000000d0"); // Secondary TrackMatching
  } else if (trainConfig == 908){
    cuts.AddCutCalo("00010113","1111111077032230000","31631031000000d0"); // std QA
  } else if (trainConfig == 909){
    cuts.AddCutCalo("00010113","1111111077022230000","21631031000000d0"); // min energy cluster variation 1  600 MeV
    cuts.AddCutCalo("00010113","1111111077042230000","21631031000000d0"); // min energy cluster variation 2  800 MeV
    cuts.AddCutCalo("00010113","1111111077052230000","21631031000000d0"); // min energy cluster variation 2  900 MeV
    cuts.AddCutCalo("00010113","1111111077032230000","2163103100000040"); // min opening angle 0.0152
    cuts.AddCutCalo("00010113","1111111077032230000","2163103100000060"); // min opening angle 0.017
  } else if (trainConfig == 910){
    cuts.AddCutCalo("00010113","1111111077032230000","2163103100000070"); // min opening angle 0.016
    cuts.AddCutCalo("00010113","1111111070032230000","21631031000000d0"); // tm off
    cuts.AddCutCalo("00010113","1111111073032230000","21631031000000d0"); // fixed window
    cuts.AddCutCalo("00010113","1111111076032230000","21631031000000d0"); // tm pt dependent var 1
    cuts.AddCutCalo("00010113","1111111078032230000","21631031000000d0"); // tm pt dependent var 2
  } else if (trainConfig == 911){ //EMCAL+DCAL+JETS
    cuts.AddCutCalo("00010113","4117700077032230000","21631031000000d0"); // INT7 - NO NL
    cuts.AddCutCalo("00010113","4117706077032230000","21631031000000d0"); // INT7 - TB NL
    cuts.AddCutCalo("00010113","4117711077032230000","21631031000000d0"); // Standard EDC
 } else if (trainConfig == 912){ //EMCAL+DCAL+JETS + new mixing
    cuts.AddCutCalo("00010113","411790007l032230000","2l631031000000d0"); // INT7 - NO NL
    cuts.AddCutCalo("00010113","411790607l032230000","2l631031000000d0"); // INT7 - TB NL
    cuts.AddCutCalo("00010113","411791107l032230000","2l631031000000d0"); // Standard EDC
 } else if (trainConfig == 913){ //// Jet QA for EMCAL+DCAL
    cuts.AddCutCalo("00010113","411791107l032230000","3l631031000000d0"); // Standard EDC INT7
    cuts.AddCutCalo("0008d113","411791107l032230000","3l631031000000d0"); // Standard EDC EG1
    cuts.AddCutCalo("0008e113","411791107l032230000","3l631031000000d0"); // Standard EDC EG2
 } else if (trainConfig == 914){ //PHOS+JETS
    cuts.AddCutCalo("00010113","2446651044012300000","2163103100000010"); //
 } else if (trainConfig == 915){ //PHOS+JetQA
    cuts.AddCutCalo("00010113","2446651044012300000","3163103100000010"); // PHOS QA
 } else if (trainConfig == 916){ //PHOS+JETS
    cuts.AddCutCalo("00010113","2446600044012300000","2163103100000010"); // PHOS No NL
 } else if (trainConfig == 917){ //MB - EMCal+JETS
    cuts.AddCutCalo("00010113","411790607l032230000","0l631031000000d0"); // MB - INT7 - TB NL
 } else if (trainConfig == 950){ // EMCal+JETS clusters standard cuts triggered analysis
    cuts.AddCutCalo("00010113","411790607l032230000","2l631031000000d0"); // INT7 - TB NL
    cuts.AddCutCalo("0008d113","411790607l032230000","2l631031000000d0"); // EG1  - TB NL
    cuts.AddCutCalo("0008e113","411790607l032230000","2l631031000000d0"); // EG2  - TB NL
    cuts.AddCutCalo("0009b113","411790607l032230000","2l631031000000d0"); // EJ1  - TB NL
    cuts.AddCutCalo("0009c113","411790607l032230000","2l631031000000d0"); // EJ2  - TB NL
 } else if (trainConfig == 951){ // EMCal+JETS cut var. min cluster energy & NCells
    cuts.AddCutCalo("00010113","411790607l022230000","2l631031000000d0"); // min energy cluster variation 1 0.6 GeV
    cuts.AddCutCalo("00010113","411790607l042230000","2l631031000000d0"); // min energy cluster variation 2 0.8 GeV
    cuts.AddCutCalo("00010113","411790607l052230000","2l631031000000d0"); // min energy cluster variation 2 0.9 GeV
    cuts.AddCutCalo("00010113","411790607l031230000","2l631031000000d0"); // std, NCell = 1
    cuts.AddCutCalo("00010113","411790607l033230000","2l631031000000d0"); // std, NCell = 3
 } else if (trainConfig == 952){ // EMCal+JETS cut var. time
    cuts.AddCutCalo("00010113","411790605l032230000","2l631031000000d0"); // std, timing = 50ns
    cuts.AddCutCalo("00010113","411790609l032230000","2l631031000000d0"); // std, timing = -20~25ns
    cuts.AddCutCalo("00010113","411790608l032230000","2l631031000000d0"); // std, timing = -20~30ns
 } else if (trainConfig == 953){ // EMCal+JETS cut var. cluster shape, TM
    cuts.AddCutCalo("00010113","411790607l032220000","2l631031000000d0"); // std, M02 = 0.1~0,7
    cuts.AddCutCalo("00010113","411790607l032250000","2l631031000000d0"); // std, M02 = 0.1~0.3
    cuts.AddCutCalo("00010113","411790607l032230000","2l631031000000d0"); // std, M02 = 0.1~0.3
    cuts.AddCutCalo("00010113","411790607e032230000","2l631031000000d0"); // std, TM fEOverPMax = 2, 	not secondary
    cuts.AddCutCalo("00010113","411790607g032230000","2l631031000000d0"); // std, TM fEOverPMax = 1.5, 	not secondary
    cuts.AddCutCalo("00010113","411790607h032230000","2l631031000000d0"); // std, TM fEOverPMax = 1.25, not secondary
    cuts.AddCutCalo("00010113","4117906077032230000","2l631031000000d0"); // std, TM not secconday
 } else if (trainConfig == 954){ // EMCal+JETS cut var. opening anlge, alpha
    cuts.AddCutCalo("00010113","411790607l032230000","2l631031000000b0"); // std, opening angle min opening angle = 0.0152
    cuts.AddCutCalo("00010113","411790607l032230000","2l631031000000g0"); // std, opening angle min opening angle = 0.0202
    cuts.AddCutCalo("00010113","411790607l032230000","2l631031000000a0"); // std, opening angle min opening angle = 0
    cuts.AddCutCalo("00010113","411790607l032230000","2l631041000000d0"); // std, alpha = 0.65
    cuts.AddCutCalo("00010113","411790607l032230000","2l631051000000d0"); // std, alpha = 0.75
    cuts.AddCutCalo("00010113","411790607l032230000","2l631061000000d0"); // std, alpha = 0.8

  } else if (trainConfig == 2000){ // EMCAL+DCAL clusters standard cuts
    cuts.AddCutCalo("00010113","4117900077032230000","01631031000000d0"); // INT7 - NO NL
    cuts.AddCutCalo("00010113","4117906077032230000","01631031000000d0"); // INT7 - TB NL
    cuts.AddCutCalo("00010113","4117911077032230000","01631031000000d0"); // Standard EDC
  } else if (trainConfig == 2001){ // EMCAL+DCAL clusters standard cuts
    cuts.AddCutCalo("00010113","411790007f032230000","01631031000000d0"); // INT7 - NO NL
    cuts.AddCutCalo("00010113","411790607f032230000","01631031000000d0"); // INT7 - TB NL
    cuts.AddCutCalo("00010113","411791107f032230000","01631031000000d0"); // Standard EDC
  } else if (trainConfig == 2002){  // EMCAL+DCAL clusters 13 TeV std. QA
    cuts.AddCutCalo("00010113","4117900017032220000","01631031000000d0"); // no timing cut, no NL INT7
    cuts.AddCutCalo("00010113","4117900077032220000","01631031000000d0"); // -30ns, 35ns timing cut, no NL INT7
  } else if (trainConfig == 2003){ // EMCAL+DCAL clusters standard cuts - QA extra stripes
    cuts.AddCutCalo("00010113","4117911077032230000","01631031000000d0"); // Standard EDC
  } else if (trainConfig == 2004){ // EMCAL+DCAL clusters standard cuts - CALOFAST TRIGGERS
    cuts.AddCutCalo("00010113","4117911070032230000","01631031000000d0"); // Standard EDC
    cuts.AddCutCalo("000a0113","4117911070032230000","01631031000000d0"); // std INT7
    cuts.AddCutCalo("000a1113","4117911070032230000","01631031000000d0"); // std EMC7
    cuts.AddCutCalo("000a2113","4117911070032230000","01631031000000d0"); // std EG2
    cuts.AddCutCalo("000a3113","4117911070032230000","01631031000000d0"); // std EG1
  } else if (trainConfig == 2005){ // Timing cut efficiency studies
    cuts.AddCutCalo("00010113","4117911007032230000","01631031000000d0"); // Standard EDC - open timing cut
    cuts.AddCutCalo("00010113","4117911077032230000","01631031000000d0"); // Standard EDC - standard timing cut
    cuts.AddCutCalo("00010113","41179110j7032230000","01631031000000d0"); // Standard EDC - Timing efficiency cut
  } else if (trainConfig == 2006){ // EMCAL+DCAL clusters pp 5 TeV Sphericity Cuts
    cuts.AddCutCalo("00010113","4117911077032230000","01631031000000d0"); // std
    cuts.AddCutCalo("h0510113","4117911077032230000","01631031000000d0"); // std
    cuts.AddCutCalo("h5a10113","4117911077032230000","01631031000000d0"); // std
    cuts.AddCutCalo("h0a10113","4117911077032230000","01631031000000d0"); // std
    cuts.AddCutCalo("h0310113","4117911077032230000","01631031000000d0"); // std
    cuts.AddCutCalo("h7a10113","4117911077032230000","01631031000000d0"); // std
  } else if (trainConfig == 2007){ // EMCAL+DCAL clusters standard cuts pp 5TeV MB
    cuts.AddCutCalo("00010113","4117911077032230000","01631031000000d0"); // Standard EDC
    cuts.AddCutCalo("00010113","4117911070032230000","01631031000000d0"); // Standard EDC, no TM
  } else if (trainConfig == 2008){ // EMCAL+DCAL clusters standard cuts pp 5TeV triggered analysis
    cuts.AddCutCalo("000a1113","4117911070032230000","01631031000000d0"); // std EMC7
    cuts.AddCutCalo("000a2113","4117911070032230000","01631031000000d0"); // std EG2
  } else if (trainConfig == 2009){ // EMCAL+DCAL clusters standard cuts pp 5TeV triggered analysis
    cuts.AddCutCalo("00010113","4117911007032230000","01631031000000d0"); // Standard EDC
    cuts.AddCutCalo("00010113","4117911000032230000","01631031000000d0"); // Standard EDC, no TM
    cuts.AddCutCalo("000a1113","4117911000032230000","01631031000000d0"); // std EMC7
    cuts.AddCutCalo("000a2113","4117911000032230000","01631031000000d0"); // std EG2

  // includes both stripes EMCal and DCal
  } else if (trainConfig == 2010){ // EMCAL+DCAL clusters standard cuts, triggers, no NL, std TM, tight timing
    cuts.AddCutCalo("00010113","4117900077032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("0008e113","4117900077032230000","01631031000000d0"); // EG2
    cuts.AddCutCalo("0008d113","4117900077032230000","01631031000000d0"); // EG1
    cuts.AddCutCalo("0009b113","4117906077032230000","01631031000000d0"); // EJ1
  } else if (trainConfig == 2011){ // EMCAL+DCAL clusters standard cuts, triggers, no NL, E/p TM, tight timing
    cuts.AddCutCalo("00010113","411790007f032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("0008e113","411790007f032230000","01631031000000d0"); // EG2
    cuts.AddCutCalo("0008d113","411790007f032230000","01631031000000d0"); // EG1
    cuts.AddCutCalo("0009b113","411790607f032230000","01631031000000d0"); // EJ1
  } else if (trainConfig == 2012){ // EMCAL+DCAL clusters standard cuts, triggers, no NL, std TM, -50ns, 30ns timing cut
    cuts.AddCutCalo("00010113","4117900057032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("0008e113","4117900057032230000","01631031000000d0"); // EG2
    cuts.AddCutCalo("0008d113","4117900057032230000","01631031000000d0"); // EG1
    cuts.AddCutCalo("0009b113","4117906057032230000","01631031000000d0"); // EJ1
  } else if (trainConfig == 2013){ // EMCAL+DCAL clusters standard cuts, triggers, no NL, E/p TM, -50ns, 30ns timing cut
    cuts.AddCutCalo("00010113","411790005f032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("0008e113","411790005f032230000","01631031000000d0"); // EG2
    cuts.AddCutCalo("0008d113","411790005f032230000","01631031000000d0"); // EG1
    cuts.AddCutCalo("0009b113","411790605f032230000","01631031000000d0"); // EJ1
  } else if (trainConfig == 2014){ // EMCAL+DCAL clusters standard cuts, triggers, no NL, std TM, open timing
    cuts.AddCutCalo("00010113","4117900007032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("0008e113","4117900007032230000","01631031000000d0"); // EG2
    cuts.AddCutCalo("0008d113","4117900007032230000","01631031000000d0"); // EG1
    cuts.AddCutCalo("0009b113","4117906007032230000","01631031000000d0"); // EJ1
  } else if (trainConfig == 2015){ // EMCAL+DCAL clusters standard cuts, triggers, no NL, E/p TM, open timing
    cuts.AddCutCalo("00010113","411790000f032230000","01631031000000d0"); // INT7
    cuts.AddCutCalo("0008e113","411790000f032230000","01631031000000d0"); // EG2
    cuts.AddCutCalo("0008d113","411790000f032230000","01631031000000d0"); // EG1
    cuts.AddCutCalo("0009b113","411790600f032230000","01631031000000d0"); // EJ1
  } else if (trainConfig == 2016){ // EMCAL+DCAL clusters standard cuts
    cuts.AddCutCalo("00010113","4117900007032230000","01631031000000d0"); // no NL, std TM, open timing
    cuts.AddCutCalo("00010113","4117900057032230000","01631031000000d0"); // no NL, std TM, -50ns, 30ns timing cut
  } else if (trainConfig == 2017){ // EMCAL+DCAL clusters standard cuts
    cuts.AddCutCalo("00010113","4117900077032230000","01631031000000d0"); // no NL, std TM, tight timing
    cuts.AddCutCalo("00010113","411790007f032230000","01631031000000d0"); // no NL, E/p TM, tight timing
    cuts.AddCutCalo("00010113","411790005f032230000","01631031000000d0"); // no NL, E/p TM, -50ns, 30ns timing cut
    cuts.AddCutCalo("00010113","411790000f032230000","01631031000000d0"); // no NL, E/p TM, open timing
  } else if (trainConfig == 2018){ // EMCAL+DCAL clusters standard cuts
    cuts.AddCutCalo("00010113","411791105f032230000","01631031000000d0"); // std TM, -50ns, 30ns timing cut
    cuts.AddCutCalo("00010113","411791205f032230000","01631031000000d0"); // std TM, -50ns, 30ns timing cut
    cuts.AddCutCalo("00010113","411792105f032230000","01631031000000d0"); // std TM, -50ns, 30ns timing cut
    cuts.AddCutCalo("00010113","411792205f032230000","01631031000000d0"); // std TM, -50ns, 30ns timing cut
  } else if (trainConfig == 2019){ // EMCAL+DCAL clusters standard cuts
    cuts.AddCutCalo("00010113","4117900007032230000","01631031000000d0"); // no NL, std TM, open timing

// EDC 13 TeV 2016 & 2017 settings with MC fine tuning correction
  } else if (trainConfig == 2020){ // EMCAL+DCAL clusters standard cuts, INT7, NL , std TM
    cuts.AddCutCalo("00010113","411791206f032230000","01631031000000d0"); // INT7 NL 12 + TB
  } else if (trainConfig == 2021){ // EMCAL+DCAL clusters standard cuts, EG2, NL , std TM
    cuts.AddCutCalo("0008e113","411791206f032230000","01631031000000d0"); // EG2  NL 12 + TB
  } else if (trainConfig == 2022){ // EMCAL+DCAL clusters standard cuts, EG1, NL , std TM
    cuts.AddCutCalo("0008d113","411791206f032230000","01631031000000d0"); // EG1  NL 12 + TB
  } else if (trainConfig == 2023){ //EMCal + DCal INT7 cut var. NonLins
    cuts.AddCutCalo("00010113","411791106f032230000","01631031000000d0"); // INT7 NL11
    cuts.AddCutCalo("00010113","411792106f032230000","01631031000000d0"); // INT7 NL21
    cuts.AddCutCalo("00010113","411792206f032230000","01631031000000d0"); // INT7 NL22
    cuts.AddCutCalo("00010113","411793506f032230000","01631031000000d0"); // INT7 NL35
    cuts.AddCutCalo("00010113","411793606f032230000","01631031000000d0"); // INT7 NL36
  } else if (trainConfig == 2024){ //EMCal + DCal INT7 cut var. time
    cuts.AddCutCalo("00010113","411791205f032230000","01631031000000d0"); // INT7 time -50+50
    cuts.AddCutCalo("00010113","411791209f032230000","01631031000000d0"); // INT7 time -20+25
    cuts.AddCutCalo("00010113","411791208f032230000","01631031000000d0"); // INT7 time -20+30
  } else if (trainConfig == 2025){ //EMCal + DCal INT7 cut var. energy and NCell
    cuts.AddCutCalo("00010113","411791206f022230000","01631031000000d0"); // INT7 energy 0.6 GeV
    cuts.AddCutCalo("00010113","411791206f042230000","01631031000000d0"); // INT7 energy 0.8 GeV
    cuts.AddCutCalo("00010113","411791206f052230000","01631031000000d0"); // INT7 energy 0.9 GeV
    cuts.AddCutCalo("00010113","411791206f031230000","01631031000000d0"); // INT7 NCells 1
    cuts.AddCutCalo("00010113","411791206f033230000","01631031000000d0"); // INT7 NCells 3
  } else if (trainConfig == 2026){ //EMCal + DCal INT7 cut var. M02 and TM
    cuts.AddCutCalo("00010113","411791206f032220000","01631031000000d0"); // INT7 M02 0.7
    cuts.AddCutCalo("00010113","411791206f032250000","01631031000000d0"); // INT7 M02 0.3
    cuts.AddCutCalo("00010113","411791206f0322k0000","01631031000000d0"); // INT7 M02 E dep
    cuts.AddCutCalo("00010113","411791206e032230000","01631031000000d0"); // INT7 TM var
    cuts.AddCutCalo("00010113","411791206g032230000","01631031000000d0"); // INT7 TM var
    cuts.AddCutCalo("00010113","411791206h032230000","01631031000000d0"); // INT7 TM var
    cuts.AddCutCalo("00010113","4117912067032230000","01631031000000d0"); // INT7 TM var
  } else if (trainConfig == 2027){ //EMCal + DCal INT7 cut var. open. angle and alpha
    cuts.AddCutCalo("00010113","411791206f032230000","01631031000000b0"); // INT7 Op. Ang. var 1 cell dist + 0.0152
    cuts.AddCutCalo("00010113","411791206f032230000","01631031000000g0"); // INT7 Op. Ang. var 1 cell dist + 0.0202
    cuts.AddCutCalo("00010113","411791206f032230000","01631031000000a0"); // INT7 Op. Ang. var 1 cell dist + 0
    cuts.AddCutCalo("00010113","411791206f032230000","01631051000000d0"); // INT7 alpha cut 0-0.75
    cuts.AddCutCalo("00010113","411791206f032230000","01631081000000d0"); // INT7 alpha cut 0-0.6

  } else if (trainConfig == 2028){ //EMCal + DCal EG2 cut var. NonLins
    cuts.AddCutCalo("0008e113","411791106f032230000","01631031000000d0"); // EG2 NL11
    cuts.AddCutCalo("0008e113","411792106f032230000","01631031000000d0"); // EG2 NL21
    cuts.AddCutCalo("0008e113","411792206f032230000","01631031000000d0"); // EG2 NL22
    cuts.AddCutCalo("0008e113","411793506f032230000","01631031000000d0"); // EG2 NL35
    cuts.AddCutCalo("0008e113","411793606f032230000","01631031000000d0"); // EG2 NL36
  } else if (trainConfig == 2029){ //EMCal + DCal EG2 cut var. time
    cuts.AddCutCalo("0008e113","411791205f032230000","01631031000000d0"); // EG2 time -50+50
    cuts.AddCutCalo("0008e113","411791209f032230000","01631031000000d0"); // EG2 time -20+25
    cuts.AddCutCalo("0008e113","411791208f032230000","01631031000000d0"); // EG2 time -20+30
  } else if (trainConfig == 2030){ //EMCal + DCal EG2 cut var. energy and NCell
    cuts.AddCutCalo("0008e113","411791206f022230000","01631031000000d0"); // EG2 energy 0.6 GeV
    cuts.AddCutCalo("0008e113","411791206f042230000","01631031000000d0"); // EG2 energy 0.8 GeV
    cuts.AddCutCalo("0008e113","411791206f052230000","01631031000000d0"); // EG2 energy 0.9 GeV
    cuts.AddCutCalo("0008e113","411791206f031230000","01631031000000d0"); // EG2 NCells 1
    cuts.AddCutCalo("0008e113","411791206f033230000","01631031000000d0"); // EG2 NCells 3
  } else if (trainConfig == 2031){ //EMCal + DCal EG2 cut var. M02 and TM
    cuts.AddCutCalo("0008e113","411791206f032220000","01631031000000d0"); // EG2 M02 0.7
    cuts.AddCutCalo("0008e113","411791206f032250000","01631031000000d0"); // EG2 M02 0.3
    cuts.AddCutCalo("0008e113","411791206f0322k0000","01631031000000d0"); // EG2 M02 E dep
    cuts.AddCutCalo("0008e113","411791206e032230000","01631031000000d0"); // EG2 TM var
    cuts.AddCutCalo("0008e113","411791206g032230000","01631031000000d0"); // EG2 TM var
    cuts.AddCutCalo("0008e113","411791206h032230000","01631031000000d0"); // EG2 TM var
    cuts.AddCutCalo("0008e113","4117912067032230000","01631031000000d0"); // EG2 TM var
  } else if (trainConfig == 2032){ //EMCal + DCal EG2 cut var. open. angle and alpha
    cuts.AddCutCalo("0008e113","411791206f032230000","01631031000000b0"); // EG2 Op. Ang. var 1 cell dist + 0.0152
    cuts.AddCutCalo("0008e113","411791206f032230000","01631031000000g0"); // EG2 Op. Ang. var 1 cell dist + 0.0202
    cuts.AddCutCalo("0008e113","411791206f032230000","01631031000000a0"); // EG2 Op. Ang. var 1 cell dist + 0
    cuts.AddCutCalo("0008e113","411791206f032230000","01631051000000d0"); // EG2 alpha cut 0-0.75
    cuts.AddCutCalo("0008e113","411791206f032230000","01631081000000d0"); // EG2 alpha cut 0-0.6

  } else if (trainConfig == 2033){ //EMCal + DCal EG1 cut var. NonLins
    cuts.AddCutCalo("0008d113","411791106f032230000","01631031000000d0"); // EG1 NL11
    cuts.AddCutCalo("0008d113","411792106f032230000","01631031000000d0"); // EG1 NL21
    cuts.AddCutCalo("0008d113","411792206f032230000","01631031000000d0"); // EG1 NL22
    cuts.AddCutCalo("0008d113","411793506f032230000","01631031000000d0"); // EG1 NL35
    cuts.AddCutCalo("0008d113","411793606f032230000","01631031000000d0"); // EG1 NL36
  } else if (trainConfig == 2034){ //EMCal + DCal EG1 cut var. time
    cuts.AddCutCalo("0008d113","411791205f032230000","01631031000000d0"); // EG1 time -50+50
    cuts.AddCutCalo("0008d113","411791209f032230000","01631031000000d0"); // EG1 time -20+25
    cuts.AddCutCalo("0008d113","411791208f032230000","01631031000000d0"); // EG1 time -20+30
  } else if (trainConfig == 2035){ //EMCal + DCal EG1 cut var. energy and NCell
    cuts.AddCutCalo("0008d113","411791206f022230000","01631031000000d0"); // EG1 energy 0.6 GeV
    cuts.AddCutCalo("0008d113","411791206f042230000","01631031000000d0"); // EG1 energy 0.8 GeV
    cuts.AddCutCalo("0008d113","411791206f052230000","01631031000000d0"); // EG1 energy 0.9 GeV
    cuts.AddCutCalo("0008d113","411791206f031230000","01631031000000d0"); // EG1 NCells 1
    cuts.AddCutCalo("0008d113","411791206f033230000","01631031000000d0"); // EG1 NCells 3
  } else if (trainConfig == 2036){ //EMCal + DCal EG1 cut var. M02 and TM
    cuts.AddCutCalo("0008d113","411791206f032220000","01631031000000d0"); // EG1 M02 0.7
    cuts.AddCutCalo("0008d113","411791206f032250000","01631031000000d0"); // EG1 M02 0.3
    cuts.AddCutCalo("0008d113","411791206f0322k0000","01631031000000d0"); // EG1 M02 E dep
    cuts.AddCutCalo("0008d113","411791206e032230000","01631031000000d0"); // EG1 TM var
    cuts.AddCutCalo("0008d113","411791206g032230000","01631031000000d0"); // EG1 TM var
    cuts.AddCutCalo("0008d113","411791206h032230000","01631031000000d0"); // EG1 TM var
    cuts.AddCutCalo("0008d113","4117912067032230000","01631031000000d0"); // EG1 TM var
  } else if (trainConfig == 2037){ //EMCal + DCal EG1 cut var. open. angle and alpha
    cuts.AddCutCalo("0008d113","411791206f032230000","01631031000000b0"); // EG1 Op. Ang. var 1 cell dist + 0.0152
    cuts.AddCutCalo("0008d113","411791206f032230000","01631031000000g0"); // EG1 Op. Ang. var 1 cell dist + 0.0202
    cuts.AddCutCalo("0008d113","411791206f032230000","01631031000000a0"); // EG1 Op. Ang. var 1 cell dist + 0
    cuts.AddCutCalo("0008d113","411791206f032230000","01631051000000d0"); // EG1 alpha cut 0-0.75
    cuts.AddCutCalo("0008d113","411791206f032230000","01631081000000d0"); // EG1 alpha cut 0-0.6
  } else if (trainConfig == 2038){ // //EMCal + DCal EG1 cut var. trigger mimick
    cuts.AddCutCalo("0008e113","411791206f032230000","01631031000000d0"); // EG2  NL 12 + TB
    cuts.AddCutCalo("0008d113","411791206f032230000","01631031000000d0"); // EG1  NL 12 + TB


  } else if (trainConfig == 2040){ //EMCal + DCal EG1 mult. diff.
    cuts.AddCutCalo("m0110113","411791206f032230000","01631031000000d0"); // INT7, NL12, mult. dep 1.0% - 2.0%
    cuts.AddCutCalo("m1210113","411791206f032230000","01631031000000d0"); // INT7, NL12, mult. dep 1.0% - 2.0%
    cuts.AddCutCalo("m2310113","411791206f032230000","01631031000000d0"); // INT7, NL12, mult. dep 2.0% - 3.0%
    cuts.AddCutCalo("m3410113","411791206f032230000","01631031000000d0"); // INT7, NL12, mult. dep 3.0% - 4.0%
    cuts.AddCutCalo("m4510113","411791206f032230000","01631031000000d0"); // INT7, NL12, mult. dep 4.0% - 5.0%
    cuts.AddCutCalo("m5710113","411791206f032230000","01631031000000d0"); // INT7, NL12, mult. dep 5.0% - 7.0%
    cuts.AddCutCalo("m7a10113","411791206f032230000","01631031000000d0"); // INT7, NL12, mult. dep 7.0% - 10.0%
  } else if (trainConfig == 2041){
    cuts.AddCutCalo("n1210113","411791206f032230000","01631031000000d0"); // INT7, NL12, mult. dep 10.0% - 20.0%
    cuts.AddCutCalo("n2510113","411791206f032230000","01631031000000d0"); // INT7, NL12, mult. dep 20.0% - 50.0%
    cuts.AddCutCalo("n5a10113","411791206f032230000","01631031000000d0"); // INT7, NL12, mult. dep 50.0% - 100.0%

  } else if (trainConfig == 2042){  // high mult trigger
    cuts.AddCutCalo("00076113","411791206f032230000","01631031000000d0"); // INT7, NL12, high mult V0M,
  } else if (trainConfig == 2043){ //EMCal + DCal EG1 mult. diff.
    cuts.AddCutCalo("q0176113","411791206f032230000","01631031000000d0"); // INT7, NL12, high mult V0M, mult. dep 0% - 0.1%
    cuts.AddCutCalo("q1276113","411791206f032230000","01631031000000d0"); // INT7, NL12, high mult V0M, mult. dep 0.1% - 0.2%
    cuts.AddCutCalo("q2376113","411791206f032230000","01631031000000d0"); // INT7, NL12, high mult V0M, mult. dep 0.2% - 0.3%
    cuts.AddCutCalo("q3476113","411791206f032230000","01631031000000d0"); // INT7, NL12, high mult V0M, mult. dep 0.3% - 0.4%
    cuts.AddCutCalo("q4576113","411791206f032230000","01631031000000d0"); // INT7, NL12, high mult V0M, mult. dep 0.4% - 0.5%
    cuts.AddCutCalo("q5776113","411791206f032230000","01631031000000d0"); // INT7, NL12, high mult V0M, mult. dep 0.5% - 0.7%
    cuts.AddCutCalo("q7a76113","411791206f032230000","01631031000000d0"); // INT7, NL12, high mult V0M, mult. dep 0.7% - 1.0%
  } else if (trainConfig == 2044){
    cuts.AddCutCalo("m0176113","411791206f032230000","01631031000000d0"); // INT7, NL12, high mult V0M, mult. dep 0.0% - 1.0%
    cuts.AddCutCalo("m1276113","411791206f032230000","01631031000000d0"); // INT7, NL12, high mult V0M, mult. dep 1.0% - 2.0%
    cuts.AddCutCalo("m2376113","411791206f032230000","01631031000000d0"); // INT7, NL12, high mult V0M, mult. dep 2.0% - 3.0%
    cuts.AddCutCalo("m3476113","411791206f032230000","01631031000000d0"); // INT7, NL12, high mult V0M, mult. dep 3.0% - 4.0%
    cuts.AddCutCalo("m4576113","411791206f032230000","01631031000000d0"); // INT7, NL12, high mult V0M, mult. dep 4.0% - 5.0%
    cuts.AddCutCalo("m5776113","411791206f032230000","01631031000000d0"); // INT7, NL12, high mult V0M, mult. dep 5.0% - 7.0%
    cuts.AddCutCalo("m7a76113","411791206f032230000","01631031000000d0"); // INT7, NL12, high mult V0M, mult. dep 7.0% - 10.0%



  } else if (trainConfig == 2045){ //EMCal + DCal EG1 sphericity.
    cuts.AddCutCalo("h0510113","411791106f032230000","01631031000000d0"); // INT7, NL12, sphericity
    cuts.AddCutCalo("h5a10113","411791106f032230000","01631031000000d0"); // INT7, NL12, sphericity
  } else if (trainConfig == 2046){ //EMCal + DCal EG1 sphericity.
    cuts.AddCutCalo("h058e113","411791106f032230000","01631031000000d0"); // EG2, NL12, sphericity
    cuts.AddCutCalo("h5a8e113","411791106f032230000","01631031000000d0"); // EG2, NL12, sphericity
  } else if (trainConfig == 2047){ //EMCal + DCal EG1 sphericity.
    cuts.AddCutCalo("h058d113","411791106f032230000","01631031000000d0"); // EG1, NL12, sphericity
    cuts.AddCutCalo("h5a8d113","411791106f032230000","01631031000000d0"); // EG1, NL12, sphericity


  } else if (trainConfig == 2050){  // EMCAL+DCAL clusters standard cuts, INT7, NL , std TM
    cuts.AddCutCalo("00010113","411791206f032230000","01631031000000d0"); // INT7 NL 12 + TB dir. gamma
  } else if (trainConfig == 2051){ // EMCAL+DCAL clusters standard cuts, EG2, NL , std TM
    cuts.AddCutCalo("0008e113","411791206f032230000","01631031000000d0"); // EG2  NL 12 + TB dir. gamma
  } else if (trainConfig == 2052){ // EMCAL+DCAL clusters standard cuts, EG1, NL , std TM
    cuts.AddCutCalo("0008d113","411791206f032230000","01631031000000d0"); // EG1  NL 12 + TB dir. gamma


// EDC settings with TB correction
  } else if (trainConfig == 2100){ // 100 MeV aggregation
    cuts.AddCutCalo("00010113","411790106f032230000","01631031000000d0"); // INT7 test beam NL
  } else if (trainConfig == 2101){ // 100 MeV aggregation
    cuts.AddCutCalo("0008e113","411790106f032230000","01631031000000d0"); // EG2  test beam NL
    cuts.AddCutCalo("0008d113","411790106f032230000","01631031000000d0"); // EG1  test beam NL
  } else if (trainConfig == 2102){ // 50 MeV aggregation
    cuts.AddCutCalo("00010113","411790206f032230000","01631031000000d0"); // INT7 test beam NL
  } else if (trainConfig == 2103){ // 50 MeV aggregation
    cuts.AddCutCalo("0008e113","411790206f032230000","01631031000000d0"); // EG2  test beam NL
    cuts.AddCutCalo("0008d113","411790206f032230000","01631031000000d0"); // EG1  test beam NL
  } else if (trainConfig == 2104){ // 150 MeV aggregation
    cuts.AddCutCalo("00010113","411790306f032230000","01631031000000d0"); // INT7 test beam NL
  } else if (trainConfig == 2105){ // 150 MeV aggregation
    cuts.AddCutCalo("0008e113","411790306f032230000","01631031000000d0"); // EG2  test beam NL
    cuts.AddCutCalo("0008d113","411790306f032230000","01631031000000d0"); // EG1  test beam NL
  } else if (trainConfig == 2106){ // 300 MeV aggregation
    cuts.AddCutCalo("00010113","411790406f032230000","01631031000000d0"); // INT7 test beam NL
  } else if (trainConfig == 2107){ // 300 MeV aggregation
    cuts.AddCutCalo("0008e113","411790406f032230000","01631031000000d0"); // EG2  test beam NL
    cuts.AddCutCalo("0008d113","411790406f032230000","01631031000000d0"); // EG1  test beam NL
  } else if (trainConfig == 2108){ // any aggregation no NL
    cuts.AddCutCalo("00010113","411790006f032230000","01631031000000d0"); // INT7 test beam NL
  } else if (trainConfig == 2109){ // any aggregation no NL
    cuts.AddCutCalo("0008e113","411790006f032230000","01631031000000d0"); // EG2  test beam NL
    cuts.AddCutCalo("0008d113","411790006f032230000","01631031000000d0"); // EG1  test beam NL

  } else if (trainConfig == 2110){ // clusterizer timing cut studies (5TeV pp std cut)
    cuts.AddCutCalo("00010113","411793106f032230000","01631031000000d0"); // INT7 test beam NL
  } else if (trainConfig == 2111){ // clusterizer timing cut studies (5TeV pp std cut)
    cuts.AddCutCalo("00010113","411793106f032230000","01631031000000d0"); // INT7 test beam NL
  } else if (trainConfig == 2112){ // clusterizer timing cut studies (5TeV pp std cut)
    cuts.AddCutCalo("00010113","411793106f032230000","01631031000000d0"); // INT7 test beam NL
  } else if (trainConfig == 2113){ // clusterizer timing cut studies (5TeV pp std cut)
    cuts.AddCutCalo("00010113","411793106f032230000","01631031000000d0"); // INT7 test beam NL
  } else if (trainConfig == 2114){ // clusterizer timing cut studies (5TeV pp std cut)
    cuts.AddCutCalo("00010113","411793106f032230000","01631031000000d0"); // INT7 test beam NL

  } else if (trainConfig == 2150){ // EMCAL clusters pp 8 TeV 100MeV aggregation
    cuts.AddCutCalo("00010113","111110106f032230000","01631031000000d0"); // std
  } else if (trainConfig == 2151){ // EMCAL clusters pp 8 TeV 100MeV aggregation
    cuts.AddCutCalo("00052113","111110106f032230000","01631031000000d0"); // std
    cuts.AddCutCalo("00081113","111110106f032230000","01631031000000d0"); // std
  } else if (trainConfig == 2152){ // EMCAL clusters pp 8 TeV 50MeV aggregation
    cuts.AddCutCalo("00010113","111110206f032230000","01631031000000d0"); // std
  } else if (trainConfig == 2153){ // EMCAL clusters pp 8 TeV 50MeV aggregation
    cuts.AddCutCalo("00052113","111110206f032230000","01631031000000d0"); // std
    cuts.AddCutCalo("00081113","111110206f032230000","01631031000000d0"); // std
  } else if (trainConfig == 2154){ // EMCAL clusters pp 8 TeV 150MeV aggregation
    cuts.AddCutCalo("00010113","111110306f032230000","01631031000000d0"); // std
  } else if (trainConfig == 2155){ // EMCAL clusters pp 8 TeV 150MeV aggregation
    cuts.AddCutCalo("00052113","111110306f032230000","01631031000000d0"); // std
    cuts.AddCutCalo("00081113","111110306f032230000","01631031000000d0"); // std
  } else if (trainConfig == 2156){ // EMCAL clusters pp 8 TeV 300MeV aggregation
    cuts.AddCutCalo("00010113","111110406f032230000","01631031000000d0"); // std
  } else if (trainConfig == 2157){ // EMCAL clusters pp 8 TeV 300MeV aggregation
    cuts.AddCutCalo("00052113","111110406f032230000","01631031000000d0"); // std
    cuts.AddCutCalo("00081113","111110406f032230000","01631031000000d0"); // std

  // New standard cut PHOS (with timing effi)
  } else if( trainConfig == 2180){
    cuts.AddCutCalo("00010113","24444000ga012200000","0163103100000010"); //
  } else if( trainConfig == 2181){
    cuts.AddCutCalo("00062113","24444000ga012200000","0163103100000010"); //
  } else if( trainConfig == 2182){
    cuts.AddCutCalo("00010113","24444590ga012200000","0163103100000010"); //
  } else if( trainConfig == 2183){
    cuts.AddCutCalo("00062113","24444590ga012200000","0163103100000010"); //
  } else if( trainConfig == 2184){
    cuts.AddCutCalo("00010113","24444690ga012200000","0163103100000010"); //
  } else if( trainConfig == 2185){
    cuts.AddCutCalo("00062113","24444690ga012200000","0163103100000010");
  // New standard cut PHOS (with timing effi) no TM
  } else if( trainConfig == 2190){
    cuts.AddCutCalo("00010113","24444000g0012200000","0163103100000010"); //
  } else if( trainConfig == 2191){
    cuts.AddCutCalo("00062113","24444000g0012200000","0163103100000010"); //
  } else if( trainConfig == 2192){
    cuts.AddCutCalo("00010113","24444590g0012200000","0163103100000010"); //
  } else if( trainConfig == 2193){
    cuts.AddCutCalo("00062113","24444590g0012200000","0163103100000010"); //
  } else if( trainConfig == 2194){
    cuts.AddCutCalo("00010113","24444690g0012200000","0163103100000010"); //
  } else if( trainConfig == 2195){
    cuts.AddCutCalo("00062113","24444690g0012200000","0163103100000010");

  // emc variations pp 8 tev for ppb reference RpA
  } else if (trainConfig == 2200) { // CALO variations
    cuts.AddCutCalo("00010113","111113106f032230000","01631031000000d0"); // std
    cuts.AddCutCalo("00010113","111113106f022230000","01631031000000d0"); // 0.6 GeV/c
    cuts.AddCutCalo("00010113","111113106f042230000","01631031000000d0"); // 0.8 GeV/c
  } else if (trainConfig == 2201) {
    cuts.AddCutCalo("00010113","111113106f031230000","01631031000000d0"); // n cells >= 1
    cuts.AddCutCalo("00010113","111113106f033230000","01631031000000d0"); // n cells >= 3
    cuts.AddCutCalo("00010113","111113106f032200000","01631031000000d0"); // no max M02 cut
    cuts.AddCutCalo("00010113","111113106f032250000","01631031000000d0"); // M02 < 0.3
    cuts.AddCutCalo("00010113","111113106f0322k0000","01631031000000d0"); // M02, pT-dep
  } else if (trainConfig == 2202) {
    cuts.AddCutCalo("00010113","111113106e032230000","01631031000000d0"); // TM var e/p 2.0
    cuts.AddCutCalo("00010113","111113106g032230000","01631031000000d0"); // TM var e/p 1.5
    cuts.AddCutCalo("00010113","111113106h032230000","01631031000000d0"); // TM var e/p 1.25
    cuts.AddCutCalo("00010113","1111131067032230000","01631031000000d0"); // TM var no veto
  } else if (trainConfig == 2203) {
    cuts.AddCutCalo("00010113","111113206f032230000","01631031000000d0"); // NL 32
    cuts.AddCutCalo("00010113","111113306f032230000","01631031000000d0"); // NL 33
    cuts.AddCutCalo("00010113","111113806f032230000","01631031000000d0"); // NL 38
    cuts.AddCutCalo("00010113","111113906f032230000","01631031000000d0"); // NL 39
  } else if (trainConfig == 2204) {
    cuts.AddCutCalo("00010113","111113105f032230000","01631031000000d0"); // 50ns
    cuts.AddCutCalo("00010113","111113104f032230000","01631031000000d0"); // 100ns
    cuts.AddCutCalo("00010113","111113107f032230000","01631031000000d0"); // 30ns

  } else if (trainConfig == 2210) { // CALO variations
    cuts.AddCutCalo("00052113","111113106f032230000","01631031000000d0"); // std
    cuts.AddCutCalo("00052113","111113106f022230000","01631031000000d0"); // 0.6 GeV/c
    cuts.AddCutCalo("00052113","111113106f042230000","01631031000000d0"); // 0.8 GeV/c
  } else if (trainConfig == 2211) {
    cuts.AddCutCalo("00052113","111113106f031230000","01631031000000d0"); // n cells >= 1
    cuts.AddCutCalo("00052113","111113106f033230000","01631031000000d0"); // n cells >= 3
    cuts.AddCutCalo("00052113","111113106f032200000","01631031000000d0"); // no max M02 cut
    cuts.AddCutCalo("00052113","111113106f032250000","01631031000000d0"); // M02 < 0.3
    cuts.AddCutCalo("00052113","111113106f0322k0000","01631031000000d0"); // M02, pT-dep
  } else if (trainConfig == 2212) {
    cuts.AddCutCalo("00052113","111113106e032230000","01631031000000d0"); // TM var e/p 2.0
    cuts.AddCutCalo("00052113","111113106g032230000","01631031000000d0"); // TM var e/p 1.5
    cuts.AddCutCalo("00052113","111113106h032230000","01631031000000d0"); // TM var e/p 1.25
    cuts.AddCutCalo("00052113","1111131067032230000","01631031000000d0"); // TM var no veto
  } else if (trainConfig == 2213) {
    cuts.AddCutCalo("00052113","111113206f032230000","01631031000000d0"); // NL 32
    cuts.AddCutCalo("00052113","111113306f032230000","01631031000000d0"); // NL 33
    cuts.AddCutCalo("00052113","111113806f032230000","01631031000000d0"); // NL 38
    cuts.AddCutCalo("00052113","111113906f032230000","01631031000000d0"); // NL 39
  } else if (trainConfig == 2214) {
    cuts.AddCutCalo("00052113","111113105f032230000","01631031000000d0"); // 50ns
    cuts.AddCutCalo("00052113","111113104f032230000","01631031000000d0"); // 100ns
    cuts.AddCutCalo("00052113","111113107f032230000","01631031000000d0"); // 30ns

  } else if (trainConfig == 2220) { // CALO variations
    cuts.AddCutCalo("00081113","111113106f032230000","01631031000000d0"); // std
    cuts.AddCutCalo("00081113","111113106f022230000","01631031000000d0"); // 0.6 GeV/c
    cuts.AddCutCalo("00081113","111113106f042230000","01631031000000d0"); // 0.8 GeV/c
  } else if (trainConfig == 2221) {
    cuts.AddCutCalo("00081113","111113106f031230000","01631031000000d0"); // n cells >= 1
    cuts.AddCutCalo("00081113","111113106f033230000","01631031000000d0"); // n cells >= 3
    cuts.AddCutCalo("00081113","111113106f032200000","01631031000000d0"); // no max M02 cut
    cuts.AddCutCalo("00081113","111113106f032250000","01631031000000d0"); // M02 < 0.3
    cuts.AddCutCalo("00081113","111113106f0322k0000","01631031000000d0"); // M02, pT-dep
  } else if (trainConfig == 2222) {
    cuts.AddCutCalo("00081113","111113106e032230000","01631031000000d0"); // TM var e/p 2.0
    cuts.AddCutCalo("00081113","111113106g032230000","01631031000000d0"); // TM var e/p 1.5
    cuts.AddCutCalo("00081113","111113106h032230000","01631031000000d0"); // TM var e/p 1.25
    cuts.AddCutCalo("00081113","1111131067032230000","01631031000000d0"); // TM var no veto
  } else if (trainConfig == 2223) {
    cuts.AddCutCalo("00081113","111113206f032230000","01631031000000d0"); // NL 32
    cuts.AddCutCalo("00081113","111113306f032230000","01631031000000d0"); // NL 33
    cuts.AddCutCalo("00081113","111113806f032230000","01631031000000d0"); // NL 38
    cuts.AddCutCalo("00081113","111113906f032230000","01631031000000d0"); // NL 39
  } else if (trainConfig == 2224) {
    cuts.AddCutCalo("00081113","111113105f032230000","01631031000000d0"); // 50ns
    cuts.AddCutCalo("00081113","111113104f032230000","01631031000000d0"); // 100ns
    cuts.AddCutCalo("00081113","111113107f032230000","01631031000000d0"); // 30ns

  } else if (trainConfig == 2230){ // EMCAL clusters pp 8 TeV, TB+finetuning CCRF
    cuts.AddCutCalo("00010113","111113106f032230000","01631031000000d0"); // std
    cuts.AddCutCalo("00052113","111113106f032230000","01631031000000d0"); // std
    cuts.AddCutCalo("00081113","111113106f032230000","01631031000000d0"); // std
  } else if (trainConfig == 2231){ // EMCAL clusters pp 8 TeV, TB+finetuning CCRF
    cuts.AddCutCalo("00010113","111113806f032230000","01631031000000d0"); // std
    cuts.AddCutCalo("00052113","111113806f032230000","01631031000000d0"); // std
    cuts.AddCutCalo("00081113","111113806f032230000","01631031000000d0"); // std
  } else if (trainConfig == 2232){ // EMCAL clusters pp 8 TeV, TB+finetuning CCRF
    cuts.AddCutCalo("00010113","111113906f032230000","01631031000000d0"); // std
    cuts.AddCutCalo("00052113","111113906f032230000","01631031000000d0"); // std
    cuts.AddCutCalo("00081113","111113906f032230000","01631031000000d0"); // std

  } else if (trainConfig == 2235){ // T0-based triggers
    cuts.AddCutCalo("00011113","111113106f032230000","01631031000000d0"); // std
    cuts.AddCutCalo("00053113","111113106f032230000","01631031000000d0"); // std
    cuts.AddCutCalo("00082113","111113106f032230000","01631031000000d0"); // std


  } else if (trainConfig == 2240){ // EMCAL clusters pp 8 TeV, TB+finetuning CCRF
    cuts.AddCutCalo("00010113","1111131060032230000","01631031000000d0"); // std
    cuts.AddCutCalo("00052113","1111131060032230000","01631031000000d0"); // std
    cuts.AddCutCalo("00081113","1111131060032230000","01631031000000d0"); // std
  } else if (trainConfig == 2241){ // EMCAL clusters pp 8 TeV, TB+finetuning CCRF
    cuts.AddCutCalo("00010113","1111138060032230000","01631031000000d0"); // std
    cuts.AddCutCalo("00052113","1111138060032230000","01631031000000d0"); // std
    cuts.AddCutCalo("00081113","1111138060032230000","01631031000000d0"); // std
  } else if (trainConfig == 2242){ // EMCAL clusters pp 8 TeV, TB+finetuning CCRF
    cuts.AddCutCalo("00010113","1111139060032230000","01631031000000d0"); // std
    cuts.AddCutCalo("00052113","1111139060032230000","01631031000000d0"); // std
    cuts.AddCutCalo("00081113","1111139060032230000","01631031000000d0"); // std

  } else if (trainConfig == 2245){ // T0 based triggers
    cuts.AddCutCalo("00011113","1111131060032230000","01631031000000d0"); // std
    cuts.AddCutCalo("00053113","1111131060032230000","01631031000000d0"); // std
    cuts.AddCutCalo("00082113","1111131060032230000","01631031000000d0"); // std


  // CONFIGS for sys pp 5 TeV MB with TM
  } else if (trainConfig == 2300) { // CALO variations
    cuts.AddCutCalo("00010113","411793106f032230000","01631031000000d0"); // std
    cuts.AddCutCalo("00010113","411793106f022230000","01631031000000d0"); // 0.6 GeV/c
    cuts.AddCutCalo("00010113","411793106f042230000","01631031000000d0"); // 0.8 GeV/c
  } else if (trainConfig == 2301) {
    cuts.AddCutCalo("00010113","411793106f031230000","01631031000000d0"); // n cells >= 1
    cuts.AddCutCalo("00010113","411793106f033230000","01631031000000d0"); // n cells >= 3
    cuts.AddCutCalo("00010113","411793106f032200000","01631031000000d0"); // no max M02 cut
    cuts.AddCutCalo("00010113","411793106f032250000","01631031000000d0"); // M02 < 0.3
    cuts.AddCutCalo("00010113","411793106f0322k0000","01631031000000d0"); // M02, pT-dep
  } else if (trainConfig == 2302) {
    cuts.AddCutCalo("00010113","411793106e032230000","01631031000000d0"); // TM var e/p 2.0
    cuts.AddCutCalo("00010113","411793106g032230000","01631031000000d0"); // TM var e/p 1.5
    cuts.AddCutCalo("00010113","411793106h032230000","01631031000000d0"); // TM var e/p 1.25
    cuts.AddCutCalo("00010113","4117931067032230000","01631031000000d0"); // TM var no veto
  } else if (trainConfig == 2303) {
    cuts.AddCutCalo("00010113","411793206f032230000","01631031000000d0"); // NL 32
    cuts.AddCutCalo("00010113","411793306f032230000","01631031000000d0"); // NL 33
    cuts.AddCutCalo("00010113","411793806f032230000","01631031000000d0"); // NL 38
    cuts.AddCutCalo("00010113","411793906f032230000","01631031000000d0"); // NL 39
  } else if (trainConfig == 2304) {
    cuts.AddCutCalo("00010113","411793109f032230000","01631031000000d0"); // 20/25ns
    cuts.AddCutCalo("00010113","411793104f032230000","01631031000000d0"); // 100ns
    cuts.AddCutCalo("00010113","411793107f032230000","01631031000000d0"); // 30ns
  } else if (trainConfig == 2305){ // opening angle variations
    cuts.AddCutCalo("00010113","411793106f032230000","01631031000000b0"); // min opening angle 0.0152
    cuts.AddCutCalo("00010113","411793106f032230000","01631031000000c0"); // min opening angle 0.016
    cuts.AddCutCalo("00010113","411793106f032230000","01631031000000e0"); // min opening angle 0.018
    cuts.AddCutCalo("00010113","411793106f032230000","01631031000000f0"); // min opening angle 0.019
  } else if (trainConfig == 2306){
    cuts.AddCutCalo("00010113","411793106f0322l0000","01631031000000d0"); // M02 pt dep with new std: 0.32, 0.0072, 0.5
    cuts.AddCutCalo("00010113","411793106f0322d0000","01631031000000d0"); // M02, pt dep with  0.27, 0.0072, 0.4
    cuts.AddCutCalo("00010113","411793106f0322e0000","01631031000000d0"); // M02, pt dep with  0.31, 0.0072, 0.5
    cuts.AddCutCalo("00010113","411793106f0322f0000","01631031000000d0"); // M02, pt dep with  0.36, 0.0072, 0.7
    cuts.AddCutCalo("00010113","411793106r0322m0000","01631031000000d0"); // M02, pt dep with  0.32, 0.0152, 0.5
  } else if (trainConfig == 2307){
    cuts.AddCutCalo("00010113","411793106f032230000","01633031000000d0"); // rapidity variation  y<0.6
    cuts.AddCutCalo("00010113","411793106f032230000","01634031000000d0"); // rapidity variation  y<0.5
    cuts.AddCutCalo("00010113","411793106f032230000","01631061000000d0"); // alpha meson variation 1   0<alpha<0.8
    cuts.AddCutCalo("00010113","411793106f032230000","01631051000000d0"); // alpha meson variation 2  0<alpha<0.75

  // CONFIGS for sys pp 5 TeV MB without TM
  } else if (trainConfig == 2310) { // CALO variations
    cuts.AddCutCalo("00010113","4117931060032230000","01631031000000d0"); // std
    cuts.AddCutCalo("00010113","4117931060022230000","01631031000000d0"); // 0.6 GeV/c
    cuts.AddCutCalo("00010113","4117931060042230000","01631031000000d0"); // 0.8 GeV/c
  } else if (trainConfig == 2311) {
    cuts.AddCutCalo("00010113","4117931060031230000","01631031000000d0"); // n cells >= 1
    cuts.AddCutCalo("00010113","4117931060033230000","01631031000000d0"); // n cells >= 3
    cuts.AddCutCalo("00010113","4117931060032200000","01631031000000d0"); // no max M02 cut
    cuts.AddCutCalo("00010113","4117931060032250000","01631031000000d0"); // M02 < 0.3
    cuts.AddCutCalo("00010113","41179310600322k0000","01631031000000d0"); // M02, pT-dep
  } else if (trainConfig == 2313) {
    cuts.AddCutCalo("00010113","4117932060032230000","01631031000000d0"); // NL 32
    cuts.AddCutCalo("00010113","4117933060032230000","01631031000000d0"); // NL 33
    cuts.AddCutCalo("00010113","4117938060032230000","01631031000000d0"); // NL 34
    cuts.AddCutCalo("00010113","4117939060032230000","01631031000000d0"); // NL 34
  } else if (trainConfig == 2314) {
    cuts.AddCutCalo("00010113","4117931090032230000","01631031000000d0"); // 20/25ns
    cuts.AddCutCalo("00010113","4117931040032230000","01631031000000d0"); // 100ns
    cuts.AddCutCalo("00010113","4117931070032230000","01631031000000d0"); // 30ns
  } else if (trainConfig == 2315){ // opening angle variations
    cuts.AddCutCalo("00010113","4117931060032230000","01631031000000b0"); // min opening angle 0.0152
    cuts.AddCutCalo("00010113","4117931060032230000","01631031000000c0"); // min opening angle 0.016
    cuts.AddCutCalo("00010113","4117931060032230000","01631031000000e0"); // min opening angle 0.018
    cuts.AddCutCalo("00010113","4117931060032230000","01631031000000f0"); // min opening angle 0.019
  } else if (trainConfig == 2316){
    cuts.AddCutCalo("00010113","41179310600322l0000","01631031000000d0"); // M02 pt dep with new std: 0.32, 0.0072, 0.5
    cuts.AddCutCalo("00010113","41179310600322d0000","01631031000000d0"); // M02, pt dep with  0.27, 0.0072, 0.4
    cuts.AddCutCalo("00010113","41179310600322e0000","01631031000000d0"); // M02, pt dep with  0.31, 0.0072, 0.5
    cuts.AddCutCalo("00010113","41179310600322f0000","01631031000000d0"); // M02, pt dep with  0.36, 0.0072, 0.7
    cuts.AddCutCalo("00010113","411793106r0322m0000","01631031000000d0"); // M02, pt dep with  0.32, 0.0152, 0.5
  } else if (trainConfig == 2317){
    cuts.AddCutCalo("00010113","4117931060032230000","01633031000000d0"); // rapidity variation  y<0.6
    cuts.AddCutCalo("00010113","4117931060032230000","01634031000000d0"); // rapidity variation  y<0.5
    cuts.AddCutCalo("00010113","4117931060032230000","01631061000000d0"); // alpha meson variation 1   0<alpha<0.8
    cuts.AddCutCalo("00010113","4117931060032230000","01631051000000d0"); // alpha meson variation 2  0<alpha<0.75

  // EMC7 configs LHC15n triggered
  } else if (trainConfig == 2320) { // CALO variations
    cuts.AddCutCalo("000a1113","4117931060032230000","01631031000000d0"); // std
    cuts.AddCutCalo("000a1113","4117931060022230000","01631031000000d0"); // 0.6 GeV/c
    cuts.AddCutCalo("000a1113","4117931060042230000","01631031000000d0"); // 0.8 GeV/c
  } else if (trainConfig == 2321) {
    cuts.AddCutCalo("000a1113","4117931060031230000","01631031000000d0"); // n cells >= 1
    cuts.AddCutCalo("000a1113","4117931060033230000","01631031000000d0"); // n cells >= 3
    cuts.AddCutCalo("000a1113","4117931060032200000","01631031000000d0"); // no max M02 cut
    cuts.AddCutCalo("000a1113","4117931060032250000","01631031000000d0"); // M02 < 0.3
    cuts.AddCutCalo("000a1113","41179310600322k0000","01631031000000d0"); // M02, pT-dep
  } else if (trainConfig == 2323) {
    cuts.AddCutCalo("000a1113","4117932060032230000","01631031000000d0"); // NL 32
    cuts.AddCutCalo("000a1113","4117933060032230000","01631031000000d0"); // NL 33
    cuts.AddCutCalo("000a1113","4117938060032230000","01631031000000d0"); // NL 34
    cuts.AddCutCalo("000a1113","4117939060032230000","01631031000000d0"); // NL 34
  } else if (trainConfig == 2324) {
    cuts.AddCutCalo("000a1113","4117931090032230000","01631031000000d0"); // 20/25ns
    cuts.AddCutCalo("000a1113","4117931040032230000","01631031000000d0"); // 100ns
    cuts.AddCutCalo("000a1113","4117931070032230000","01631031000000d0"); // 30ns
  } else if (trainConfig == 2325){ // opening angle variations
    cuts.AddCutCalo("000a1113","4117931060032230000","01631031000000b0"); // min opening angle 0.0152
    cuts.AddCutCalo("000a1113","4117931060032230000","01631031000000c0"); // min opening angle 0.016
    cuts.AddCutCalo("000a1113","4117931060032230000","01631031000000e0"); // min opening angle 0.018
    cuts.AddCutCalo("000a1113","4117931060032230000","01631031000000f0"); // min opening angle 0.019
  } else if (trainConfig == 2326){
    cuts.AddCutCalo("000a1113","41179310600322l0000","01631031000000d0"); // M02 pt dep with new std: 0.32, 0.0072, 0.5
    cuts.AddCutCalo("000a1113","41179310600322d0000","01631031000000d0"); // M02, pt dep with  0.27, 0.0072, 0.4
    cuts.AddCutCalo("000a1113","41179310600322e0000","01631031000000d0"); // M02, pt dep with  0.31, 0.0072, 0.5
    cuts.AddCutCalo("000a1113","41179310600322f0000","01631031000000d0"); // M02, pt dep with  0.36, 0.0072, 0.7
    cuts.AddCutCalo("000a1113","411793106r0322m0000","01631031000000d0"); // M02, pt dep with  0.32, 0.0152, 0.5
  } else if (trainConfig == 2327){
    cuts.AddCutCalo("000a1113","4117931060032230000","01633031000000d0"); // rapidity variation  y<0.6
    cuts.AddCutCalo("000a1113","4117931060032230000","01634031000000d0"); // rapidity variation  y<0.5
    cuts.AddCutCalo("000a1113","4117931060032230000","01631061000000d0"); // alpha meson variation 1   0<alpha<0.8
    cuts.AddCutCalo("000a1113","4117931060032230000","01631051000000d0"); // alpha meson variation 2  0<alpha<0.75

  // EG2 configs LHC17pq triggered
  } else if (trainConfig == 2330) { // CALO variations
    cuts.AddCutCalo("000a2113","4117931060032230000","01631031000000d0"); // std
    cuts.AddCutCalo("000a2113","4117931060022230000","01631031000000d0"); // 0.6 GeV/c
    cuts.AddCutCalo("000a2113","4117931060042230000","01631031000000d0"); // 0.8 GeV/c
  } else if (trainConfig == 2331) {
    cuts.AddCutCalo("000a2113","4117931060031230000","01631031000000d0"); // n cells >= 1
    cuts.AddCutCalo("000a2113","4117931060033230000","01631031000000d0"); // n cells >= 3
    cuts.AddCutCalo("000a2113","4117931060032200000","01631031000000d0"); // no max M02 cut
    cuts.AddCutCalo("000a2113","4117931060032250000","01631031000000d0"); // M02 < 0.3
    cuts.AddCutCalo("000a2113","41179310600322k0000","01631031000000d0"); // M02, pT-dep
  } else if (trainConfig == 2333) {
    cuts.AddCutCalo("000a2113","4117932060032230000","01631031000000d0"); // NL 32
    cuts.AddCutCalo("000a2113","4117933060032230000","01631031000000d0"); // NL 33
    cuts.AddCutCalo("000a2113","4117938060032230000","01631031000000d0"); // NL 34
    cuts.AddCutCalo("000a2113","4117939060032230000","01631031000000d0"); // NL 34
  } else if (trainConfig == 2334) {
    cuts.AddCutCalo("000a2113","4117931090032230000","01631031000000d0"); // 20/25ns
    cuts.AddCutCalo("000a2113","4117931040032230000","01631031000000d0"); // 100ns
    cuts.AddCutCalo("000a2113","4117931070032230000","01631031000000d0"); // 30ns
  } else if (trainConfig == 2335){ // opening angle variations
    cuts.AddCutCalo("000a2113","4117931060032230000","01631031000000b0"); // min opening angle 0.0152
    cuts.AddCutCalo("000a2113","4117931060032230000","01631031000000c0"); // min opening angle 0.016
    cuts.AddCutCalo("000a2113","4117931060032230000","01631031000000e0"); // min opening angle 0.018
    cuts.AddCutCalo("000a2113","4117931060032230000","01631031000000f0"); // min opening angle 0.019
  } else if (trainConfig == 2336){
    cuts.AddCutCalo("000a2113","41179310600322l0000","01631031000000d0"); // M02 pt dep with new std: 0.32, 0.0072, 0.5
    cuts.AddCutCalo("000a2113","41179310600322d0000","01631031000000d0"); // M02, pt dep with  0.27, 0.0072, 0.4
    cuts.AddCutCalo("000a2113","41179310600322e0000","01631031000000d0"); // M02, pt dep with  0.31, 0.0072, 0.5
    cuts.AddCutCalo("000a2113","41179310600322f0000","01631031000000d0"); // M02, pt dep with  0.36, 0.0072, 0.7
    cuts.AddCutCalo("000a2113","411793106r0322m0000","01631031000000d0"); // M02, pt dep with  0.32, 0.0152, 0.5
  } else if (trainConfig == 2337){
    cuts.AddCutCalo("000a2113","4117931060032230000","01633031000000d0"); // rapidity variation  y<0.6
    cuts.AddCutCalo("000a2113","4117931060032230000","01634031000000d0"); // rapidity variation  y<0.5
    cuts.AddCutCalo("000a2113","4117931060032230000","01631061000000d0"); // alpha meson variation 1   0<alpha<0.8
    cuts.AddCutCalo("000a2113","4117931060032230000","01631051000000d0"); // alpha meson variation 2  0<alpha<0.75

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

  TList *EventCutList = new TList();
  TList *ClusterCutList = new TList();
  TList *MesonCutList = new TList();

  TList *HeaderList = new TList();
  if (generatorName.Contains("LHC12i3")){
    TObjString *Header2 = new TObjString("BOX");
    HeaderList->Add(Header2);
  } else if (generatorName.CompareTo("LHC14e2b")==0){
    TObjString *Header2 = new TObjString("pi0_1");
    HeaderList->Add(Header2);
    TObjString *Header3 = new TObjString("eta_2");
    HeaderList->Add(Header3);
  }

  TString energy = "";
  TString mcName = "";
  TString mcNameAdd = "";
  if (generatorName.Contains("WOSDD")){
    mcNameAdd = "_WOSDD";
  } else if (generatorName.Contains("WSDD")){
    mcNameAdd = "_WSDD";
  }
  if (generatorName.Contains("LHC12i3")){
    energy = "2760GeV";
    mcName = "Pythia8_LHC12i3";
  } else if (generatorName.Contains("LHC12f1a")){
    energy = "2760GeV";
    mcName = "Pythia8_LHC12f1a";
  } else if (generatorName.Contains("LHC12f1b")){
    energy = "2760GeV";
    mcName = "Phojet_LHC12f1b";
  } else if (generatorName.Contains("LHC14e2a")){
    energy = "8TeV";
    mcName = "Pythia8_LHC14e2a";
  } else if (generatorName.Contains("LHC14e2b")){
    energy = "8TeV";
    mcName = "Pythia8_LHC14e2b";
  } else if (generatorName.Contains("LHC14e2c")){
    energy = "8TeV";
    mcName = "Phojet_LHC14e2c";
  } else if (generatorName.Contains("LHC16c2")){
    energy            = "8TeV";
    mcName            = "LHC16c2";
  } else if (generatorName.Contains("LHC16h3")){
    energy            = "5TeV";
    mcName            = "PythiaJets_LHC16h3";
  } else if (generatorName.Contains("LHC18b8")){
    energy            = "5TeV";
    mcName            = "PythiaJets_LHC18b8";
  }

  EventCutList->SetOwner(kTRUE);
  AliConvEventCuts **analysisEventCuts = new AliConvEventCuts*[numberOfCuts];
  ClusterCutList->SetOwner(kTRUE);
  AliCaloPhotonCuts **analysisClusterCuts = new AliCaloPhotonCuts*[numberOfCuts];
  MesonCutList->SetOwner(kTRUE);
  AliConversionMesonCuts **analysisMesonCuts = new AliConversionMesonCuts*[numberOfCuts];

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

    // definition of weighting input
    TString fitNamePi0 = Form("Pi0_Fit_Data_%s",energy.Data());
    TString fitNameEta = Form("Eta_Fit_Data_%s",energy.Data());
    TString fAddedSignalString = cuts.GetEventCut(i);
    fAddedSignalString = fAddedSignalString(6,1);
    Bool_t fAddedSignal = kFALSE;
    if (fAddedSignalString.CompareTo("2") == 0) fAddedSignal = kTRUE;

    TString mcInputNamePi0 = "";
    TString mcInputNameEta = "";
    if (fAddedSignal && (generatorName.Contains("LHC12i3") || generatorName.CompareTo("LHC14e2b")==0)){
      mcInputNamePi0 = Form("Pi0_%s%s_addSig_%s", mcName.Data(), mcNameAdd.Data(), energy.Data() );
      mcInputNameEta = Form("Eta_%s%s_addSig_%s", mcName.Data(), mcNameAdd.Data(), energy.Data() );
    } else {
      mcInputNamePi0 = Form("Pi0_%s%s_%s", mcName.Data(), mcNameAdd.Data(), energy.Data() );
      mcInputNameEta = Form("Eta_%s%s_%s", mcName.Data(), mcNameAdd.Data(), energy.Data() );
    }

    if (doWeightingPart) analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE, fileNamePtWeights, mcInputNamePi0, mcInputNameEta, "",fitNamePi0,fitNameEta);

    TString dataInputMultHisto  = "";
    TString mcInputMultHisto    = "";
    TString triggerString   = cuts.GetEventCut(i);
    triggerString           = triggerString(3,2);
    if (triggerString.CompareTo("03")==0)
      triggerString         = "00";
    if (periodNameAnchor.CompareTo("LHC13g") == 0 && triggerString.CompareTo("10")== 0 )
      triggerString         = "00";


    dataInputMultHisto      = Form("%s_%s", periodNameAnchor.Data(), triggerString.Data());
    mcInputMultHisto        = Form("%s_%s", periodNameV0Reader.Data(), triggerString.Data());

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
    if (enableLightOutput > 0) analysisEventCuts[i]->SetLightOutput(kTRUE);
    if(trainConfig == 447) analysisEventCuts[i]->SetUseSphericityTrue(kTRUE);
    if (periodNameV0Reader.CompareTo("") != 0) analysisEventCuts[i]->SetPeriodEnum(periodNameV0Reader);
    analysisEventCuts[i]->InitializeCutsFromCutString((cuts.GetEventCut(i)).Data());
    EventCutList->Add(analysisEventCuts[i]);
    analysisEventCuts[i]->SetFillCutHistograms("",kFALSE);

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
  task->SetDoMesonAnalysis(kTRUE);
  task->SetCorrectionTaskSetting(corrTaskSetting);
  task->SetDoMesonQA(enableQAMesonTask); //Attention new switch for Pi0 QA
  task->SetDoClusterQA(enableQAClusterTask);  //Attention new switch small for Cluster QA
  task->SetDoTHnSparse(enableTHnSparse);
  task->SetProduceTreeEOverP(doTreeEOverP);
  task->SetEnableSortingOfMCClusLabels(enableSortingMCLabels);
  if(trainConfig == 2002 || (trainConfig >= 2020 && trainConfig <= 2045)) task->SetDoPi0Only(kTRUE);
  if(trainConfig == 446) task->SetSoftAnalysis(kTRUE);
  if(enableExtMatchAndQA > 1){ task->SetPlotHistsExtQA(kTRUE);}
  if(trainConfig == 106 || trainConfig == 125 || trainConfig == 145){
    task->SetInOutTimingCluster(-30e-9,35e-9);
  }
  task->SetLocalDebugFlag(localDebugFlag);

  //connect containers
  AliAnalysisDataContainer *coutput =
    mgr->CreateContainer(!(corrTaskSetting.CompareTo("")) ? Form("GammaCalo_%i",trainConfig) : Form("GammaCalo_%i_%s",trainConfig,corrTaskSetting.Data()), TList::Class(),
              AliAnalysisManager::kOutputContainer,Form("GammaCalo_%i.root",trainConfig));

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
