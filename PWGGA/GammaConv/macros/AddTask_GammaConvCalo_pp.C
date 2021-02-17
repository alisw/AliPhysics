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
//pp together with all supporting classes
//***************************************************************************************

//***************************************************************************************
//main function
//***************************************************************************************
void AddTask_GammaConvCalo_pp(
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
  // FPTW:fileNamePtWeights, FMUW:fileNameMultWeights,  FMAW:fileNameMatBudWeights,  separate with ;
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
  Bool_t    enableTreeConvGammaShape      = kFALSE,   // enable additional tree for conversion properties for clusters
  Bool_t    doSmear                       = kFALSE,   // switches to run user defined smearing
  Double_t  bremSmear                     = 1.,
  Double_t  smearPar                      = 0.,       // conv photon smearing params
  Double_t  smearParConst                 = 0.,       // conv photon smearing params
  Bool_t    doPrimaryTrackMatching        = kTRUE,    // enable basic track matching for all primary tracks to cluster
  // subwagon config
  TString   additionalTrainConfig         = "0"       // additional counter for trainconfig
  ) {

  AliCutHandlerPCM cuts;


  Int_t TriggerMimickingDDLEffiFlag = 2;
  if (enableTriggerMimicking >= 10) {
      if (enableTriggerMimicking >= 20) {
          TriggerMimickingDDLEffiFlag = 0;
          enableTriggerMimicking -= 20;
      } else {
          TriggerMimickingDDLEffiFlag = 1;
          enableTriggerMimicking -= 10;
      }
  }

  TString fileNamePtWeights     = cuts.GetSpecialFileNameFromString (fileNameExternalInputs, "FPTW:");
  TString fileNameMultWeights   = cuts.GetSpecialFileNameFromString (fileNameExternalInputs, "FMUW:");
  TString fileNameMatBudWeights = cuts.GetSpecialFileNameFromString (fileNameExternalInputs, "FMAW:");
  TString fileNamedEdxPostCalib = cuts.GetSpecialFileNameFromString (fileNameExternalInputs, "FEPC:");
  TString fileNameCustomTriggerMimicOADB   = cuts.GetSpecialFileNameFromString (fileNameExternalInputs, "FTRM:");

  TString addTaskName                 = "AddTask_GammaConvCalo_pp";
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

  Bool_t doTreeEOverP = kFALSE; // switch to produce EOverP tree
  TString strdoTreeEOverP             = cuts.GetSpecialSettingFromAddConfig(additionalTrainConfig, "EPCLUSTree", "", addTaskName);
  if(strdoTreeEOverP.Atoi()==1)
    doTreeEOverP = kTRUE;

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
  if(rmaxFacPtHardSetting->GetEntries()<1){cout << "ERROR: AddTask_GammaConvCalo_pp during parsing of settingMaxFacPtHard String '" << settingMaxFacPtHard.Data() << "'" << endl; return;}
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
  AliAnalysisTaskGammaConvCalo *task=NULL;
  task= new AliAnalysisTaskGammaConvCalo(Form("GammaConvCalo_%i",trainConfig));
  task->SetIsHeavyIon(isHeavyIon);
  task->SetIsMC(isMC);
  task->SetV0ReaderName(V0ReaderName);
  if(enableLightOutput >= 20){  // eta
    enableLightOutput -= 20;
    task->SetPi0EtaSwitch(2);
  } else if(enableLightOutput >= 10){  // pi0
    enableLightOutput -= 10;
    task->SetPi0EtaSwitch(1);
  } else {
    task->SetPi0EtaSwitch(0);
  }
  if (enableLightOutput > 1 && enableLightOutput != 4) task->SetLightOutput(1);
  else if (enableLightOutput == 4) task->SetLightOutput(2);
  if (enableLightOutput == 5) task->SetECalibOutput(kTRUE);
  task->SetDoPrimaryTrackMatching(doPrimaryTrackMatching);
  task->SetTrackMatcherRunningMode(trackMatcherRunningMode);
  if(trainConfig >= 850 && trainConfig < 860) task->SetDoHBTHistoOutput(kTRUE);

  // cluster cuts
  // 0 "ClusterType",  1 "EtaMin", 2 "EtaMax", 3 "PhiMin", 4 "PhiMax", 5 "DistanceToBadChannel", 6 "Timing", 7 "TrackMatching", 8 "ExoticCell",
  // 9 "MinEnergy", 10 "MinNCells", 11 "MinM02", 12 "MaxM02", 13 "MinMaxM20", 14 "RecConv", 15 "MaximumDispersion", 16 "NLM"

  // *****************************************************************************************************
  // pp 2.76 TeV setup, paper cuts
  // *****************************************************************************************************
  if (trainConfig == 1){ // EMCAL clusters 2.76 TeV LHC11a, with SDD (0), kEMC1 (1) final analysis cuts
    cuts.AddCutPCMCalo("00003113","00200009327000008250400000","1111121057032230000","0163103100000010");
    cuts.AddCutPCMCalo("00051013","00200009327000008250400000","1111121057032230000","0163103100000010");
  } else if ( trainConfig == 2){ // LHC11a no non linearity
    cuts.AddCutPCMCalo("00003113","00200009327000008250400000","1111100057032230000","0163103100000010");
    cuts.AddCutPCMCalo("00051013","00200009327000008250400000","1111100057032230000","0163103100000010");
  } else if ( trainConfig == 3){  // LHC13g final analysis cuts
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111121067032230000","0163103100000010"); // INT7
    cuts.AddCutPCMCalo("00052013","00200009327000008250400000","1111121067032230000","0163103100000010"); // EMC7
    cuts.AddCutPCMCalo("00083013","00200009327000008250400000","1111121067032230000","0163103100000010"); // EMCEG1,
    cuts.AddCutPCMCalo("00085013","00200009327000008250400000","1111121067032230000","0163103100000010"); // EMCEG2,
  } else if ( trainConfig == 4){  // EMCal, all triggers // LHC13g new conv calo non lienarity with pileup
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111121067032230000","0163103100000010"); // INT7
    cuts.AddCutPCMCalo("00052113","00200009327000008250400000","1111121067032230000","0163103100000010"); // EMC7
    cuts.AddCutPCMCalo("00083113","00200009327000008250400000","1111121067032230000","0163103100000010"); // EMCEG1,
    cuts.AddCutPCMCalo("00085113","00200009327000008250400000","1111121067032230000","0163103100000010"); // EMCEG2,
  } else if ( trainConfig == 5){  // EMCal, all triggers without non linearity
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111100067032230000","0163103100000010"); // INT7
    cuts.AddCutPCMCalo("00052013","00200009327000008250400000","1111100067032230000","0163103100000010"); // EMC7
    cuts.AddCutPCMCalo("00083013","00200009327000008250400000","1111100067032230000","0163103100000010"); // EMCEG1,
    cuts.AddCutPCMCalo("00085013","00200009327000008250400000","1111100067032230000","0163103100000010"); // EMCEG2,

  // INT1 variations
  } else if ( trainConfig == 10){ //EMCal acceptance variations
    cuts.AddCutPCMCalo("00003113","00200009327000008250400000","1113111057032230000","0163103100000010"); //only modules with TRD infront
    cuts.AddCutPCMCalo("00003113","00200009327000008250400000","1111211057032230000","0163103100000010"); //no modules with TRD infront
  } else if ( trainConfig == 11){ // NonLinearity variations
    cuts.AddCutPCMCalo("00003113","00200009327000008250400000","1111100057032230000","0163103100000010"); // NonLinearity none
    cuts.AddCutPCMCalo("00003113","00200009327000008250400000","1111101057032230000","0163103100000010"); // NonLinearity kSDMv5
    cuts.AddCutPCMCalo("00003113","00200009327000008250400000","1111122057032230000","0163103100000010"); // NonLinearity CMF
    cuts.AddCutPCMCalo("00003113","00200009327000008250400000","1111111057032230000","0163103100000010"); // NonLinearity CCRF
    cuts.AddCutPCMCalo("00003113","00200009327000008250400000","1111112057032230000","0163103100000010"); // NonLinearity CRF
  // EMC1 variations
  } else if ( trainConfig == 12){ //EMCal acceptance variations
    cuts.AddCutPCMCalo("00051013","00200009327000008250400000","1113111057032230000","0163103100000010"); //only modules with TRD infront
    cuts.AddCutPCMCalo("00051013","00200009327000008250400000","1111211057032230000","0163103100000010"); //no modules with TRD infront
  } else if ( trainConfig == 13){ // LHC11a NonLinearity variations
    cuts.AddCutPCMCalo("00051013","00200009327000008250400000","1111100057032230000","0163103100000010"); // NonLinearity none
    cuts.AddCutPCMCalo("00051013","00200009327000008250400000","1111101057032230000","0163103100000010"); // NonLinearity kSDMv5
    cuts.AddCutPCMCalo("00051013","00200009327000008250400000","1111122057032230000","0163103100000010"); // NonLinearity CMF
    cuts.AddCutPCMCalo("00051013","00200009327000008250400000","1111111057032230000","0163103100000010"); // NonLinearity CCRF
    cuts.AddCutPCMCalo("00051013","00200009327000008250400000","1111112057032230000","0163103100000010"); // NonLinearity CRF
  // INT7 variations
  } else if ( trainConfig == 14){ //EMCal acceptance variations
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1112111067032230000","0163103100000010"); //only modules with TRD infront
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111311067032230000","0163103100000010"); //no modules with TRD infront
  } else if ( trainConfig == 15){  //LHC11a NonLinearity variations
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111100067032230000","0163103100000010"); // NonLinearity none
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111101067032230000","0163103100000010"); // NonLinearity kSDMv5
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111122067032230000","0163103100000010"); // NonLinearity CMF
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111111067032230000","0163103100000010"); // NonLinearity CCRF
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111112067032230000","0163103100000010"); // NonLinearity CRF
  // EMC7 variations
  } else if ( trainConfig == 16){ //EMCal acceptance variations
    cuts.AddCutPCMCalo("00052013","00200009327000008250400000","1112111067032230000","0163103100000010"); //only modules with TRD infront
    cuts.AddCutPCMCalo("00052013","00200009327000008250400000","1111311067032230000","0163103100000010"); //no modules with TRD infront
  } else if ( trainConfig == 17){  //LHC11a NonLinearity variations
    cuts.AddCutPCMCalo("00052013","00200009327000008250400000","1111100067032230000","0163103100000010"); // NonLinearity none
    cuts.AddCutPCMCalo("00052013","00200009327000008250400000","1111101067032230000","0163103100000010"); // NonLinearity kSDMv5
    cuts.AddCutPCMCalo("00052013","00200009327000008250400000","1111122067032230000","0163103100000010"); // NonLinearity CMF
    cuts.AddCutPCMCalo("00052013","00200009327000008250400000","1111111067032230000","0163103100000010"); // NonLinearity CCRF
    cuts.AddCutPCMCalo("00052013","00200009327000008250400000","1111112067032230000","0163103100000010"); // NonLinearity CRF
  // EG2 variations
  } else if ( trainConfig == 18){ //EMCal acceptance variations
    cuts.AddCutPCMCalo("00085013","00200009327000008250400000","1112111067032230000","0163103100000010"); //only modules with TRD infront
    cuts.AddCutPCMCalo("00085013","00200009327000008250400000","1111311067032230000","0163103100000010"); //no modules with TRD infront
  } else if ( trainConfig == 19){  //LHC11a NonLinearity variations
    cuts.AddCutPCMCalo("00085013","00200009327000008250400000","1111100067032230000","0163103100000010"); // NonLinearity none
    cuts.AddCutPCMCalo("00085013","00200009327000008250400000","1111101067032230000","0163103100000010"); // NonLinearity kSDMv5
    cuts.AddCutPCMCalo("00085013","00200009327000008250400000","1111122067032230000","0163103100000010"); // NonLinearity CMF
    cuts.AddCutPCMCalo("00085013","00200009327000008250400000","1111111067032230000","0163103100000010"); // NonLinearity CCRF
    cuts.AddCutPCMCalo("00085013","00200009327000008250400000","1111112067032230000","0163103100000010"); // NonLinearity CRF
  // EG2 variations
  } else if ( trainConfig == 20){ //EMCal acceptance variations
    cuts.AddCutPCMCalo("00083013","00200009327000008250400000","1112111067032230000","0163103100000010"); //only modules with TRD infront
    cuts.AddCutPCMCalo("00083013","00200009327000008250400000","1111311067032230000","0163103100000010"); //no modules with TRD infront
  } else if ( trainConfig == 21){  //LHC11a NonLinearity variations
    cuts.AddCutPCMCalo("00083013","00200009327000008250400000","1111100067032230000","0163103100000010"); // NonLinearity none
    cuts.AddCutPCMCalo("00083013","00200009327000008250400000","1111101067032230000","0163103100000010"); // NonLinearity kSDMv5
    cuts.AddCutPCMCalo("00083013","00200009327000008250400000","1111122067032230000","0163103100000010"); // NonLinearity CMF
    cuts.AddCutPCMCalo("00083013","00200009327000008250400000","1111111067032230000","0163103100000010"); // NonLinearity CCRF
    cuts.AddCutPCMCalo("00083013","00200009327000008250400000","1111112067032230000","0163103100000010"); // NonLinearity CRF

  // Configurations without non lin
  } else if ( trainConfig == 31){  // LHC12 without non linearity
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111100067032230000","0163103100000010"); // INT7
    cuts.AddCutPCMCalo("00052113","00200009327000008250400000","1111100067032230000","0163103100000010"); // EMC7
    cuts.AddCutPCMCalo("00081113","00200009327000008250400000","1111100067032230000","0163103100000010"); // EMCEG1,
  } else if ( trainConfig == 32){  // LHC10 without non linearity
    cuts.AddCutPCMCalo("00000113","00200009327000008250400000","1111100017032230000","0163103100000010"); // MB


  // Multiplicity dependent cuts
  } else if ( trainConfig == 40){ // MB - with multiplicity bins
    cuts.AddCutPCMCalo("00103113","00200009327000008250400000","1111121053032230000","0163103100000010"); // 0 -2
    cuts.AddCutPCMCalo("01203113","00200009327000008250400000","1111121053032230000","0163103100000010"); // 2 -5
    cuts.AddCutPCMCalo("02303113","00200009327000008250400000","1111121053032230000","0163103100000010"); // 5 -10
    cuts.AddCutPCMCalo("03403113","00200009327000008250400000","1111121053032230000","0163103100000010"); // 10 -30
    cuts.AddCutPCMCalo("04503113","00200009327000008250400000","1111121053032230000","0163103100000010"); // 30 -100
  } else if ( trainConfig == 41){ // INT7 - with multiplicity bins
    cuts.AddCutPCMCalo("00110113","00200009327000008250400000","1111121063032230000","0163103100000010"); // 0 -2
    cuts.AddCutPCMCalo("01210113","00200009327000008250400000","1111121063032230000","0163103100000010"); // 2 -5
    cuts.AddCutPCMCalo("02310113","00200009327000008250400000","1111121063032230000","0163103100000010"); // 5 -10
    cuts.AddCutPCMCalo("03410113","00200009327000008250400000","1111121063032230000","0163103100000010"); // 10 -30
    cuts.AddCutPCMCalo("04510113","00200009327000008250400000","1111121063032230000","0163103100000010"); // 30 -100

  //*************************************************************************************************
  // 8 TeV EMC setup
  //*************************************************************************************************
  } else if ( trainConfig == 80){ // EMCAL clusters 8 TeV LHC12, TB 100 MeV
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0400000","111110106f032230000","0163103100000010"); // std
    cuts.AddCutPCMCalo("00052113","00200009f9730000dge0400000","111110106f032230000","0163103100000010"); // std
    cuts.AddCutPCMCalo("00081113","00200009f9730000dge0400000","111110106f032230000","0163103100000010"); // std
  } else if ( trainConfig == 81){ // EMCAL clusters 8 TeV LHC12, TB 50 MeV
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0400000","111110206f032230000","0163103100000010"); // std
    cuts.AddCutPCMCalo("00052113","00200009f9730000dge0400000","111110206f032230000","0163103100000010"); // std
    cuts.AddCutPCMCalo("00081113","00200009f9730000dge0400000","111110206f032230000","0163103100000010"); // std
  } else if ( trainConfig == 82){ // EMCAL clusters 8 TeV LHC12, TB 150 MeV
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0400000","111110306f032230000","0163103100000010"); // std
    cuts.AddCutPCMCalo("00052113","00200009f9730000dge0400000","111110306f032230000","0163103100000010"); // std
    cuts.AddCutPCMCalo("00081113","00200009f9730000dge0400000","111110306f032230000","0163103100000010"); // std
  } else if ( trainConfig == 83){ // EMCAL clusters 8 TeV LHC12, TB 300 MeV
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0400000","111110406f032230000","0163103100000010"); // std
    cuts.AddCutPCMCalo("00052113","00200009f9730000dge0400000","111110406f032230000","0163103100000010"); // std
    cuts.AddCutPCMCalo("00081113","00200009f9730000dge0400000","111110406f032230000","0163103100000010"); // std
  } else if ( trainConfig == 84){ // EMCAL clusters 8 TeV LHC12, smearing
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","111113206fg32230000","0163103100b00010"); // std
    cuts.AddCutPCMCalo("00052113","0dm00009f9730000dge0404000","111113206fg32230000","0163103100b00010"); // std
    cuts.AddCutPCMCalo("00081113","0dm00009f9730000dge0404000","111113206fg32230000","0163103100b00010"); // std
  } else if ( trainConfig == 85){ // EMCAL clusters 8 TeV LHC12, smearing
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","111113206f032230000","0163103100b00010"); // std
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","111113206f03h230000","0163103100b00010"); // std
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","111113206f03k230000","0163103100b00010"); // std


  } else if ( trainConfig == 86){ // T0-based triggers
    cuts.AddCutPCMCalo("00011113","0dm00009f9730000dge0404000","111113206fg32230000","0163103100000010"); // std
    cuts.AddCutPCMCalo("00053113","0dm00009f9730000dge0404000","111113206fg32230000","0163103100000010"); // std
    cuts.AddCutPCMCalo("00082113","0dm00009f9730000dge0404000","111113206fg32230000","0163103100000010"); // std
  } else if ( trainConfig == 87){ // EMCAL clusters 8 TeV LHC12, TB 300 MeV
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","111113206fg32230000","0163103100000010"); // std
    cuts.AddCutPCMCalo("00052113","0dm00009f9730000dge0404000","111113206fg32230000","0163103100000010"); // std
    cuts.AddCutPCMCalo("00081113","0dm00009f9730000dge0404000","111113206fg32230000","0163103100000010"); // std


  } else if ( trainConfig == 88){ // T0-based triggers
    cuts.AddCutPCMCalo("00011113","0dm00009f9730000dge0404000","111113206f032230000","0163103100000010"); // std
    cuts.AddCutPCMCalo("00053113","0dm00009f9730000dge0404000","111113206f032230000","0163103100000010"); // std
    cuts.AddCutPCMCalo("00082113","0dm00009f9730000dge0404000","111113206f032230000","0163103100000010"); // std

  } else if ( trainConfig == 90){ // EMCAL clusters 8 TeV LHC12, TB 100 MeV
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0400000","111116106f032230000","0163103100000010"); // std
    cuts.AddCutPCMCalo("00052113","00200009f9730000dge0400000","111116106f032230000","0163103100000010"); // std
    cuts.AddCutPCMCalo("00081113","00200009f9730000dge0400000","111116106f032230000","0163103100000010"); // std
  } else if ( trainConfig == 91){ // EMCAL clusters 8 TeV LHC12, TB 50 MeV
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0400000","111116206f032230000","0163103100000010"); // std
    cuts.AddCutPCMCalo("00052113","00200009f9730000dge0400000","111116206f032230000","0163103100000010"); // std
    cuts.AddCutPCMCalo("00081113","00200009f9730000dge0400000","111116206f032230000","0163103100000010"); // std
  } else if ( trainConfig == 92){ // EMCAL clusters 8 TeV LHC12, TB 150 MeV
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0400000","111116306f032230000","0163103100000010"); // std
    cuts.AddCutPCMCalo("00052113","00200009f9730000dge0400000","111116306f032230000","0163103100000010"); // std
    cuts.AddCutPCMCalo("00081113","00200009f9730000dge0400000","111116306f032230000","0163103100000010"); // std
  } else if ( trainConfig == 93){ // EMCAL clusters 8 TeV LHC12, TB 300 MeV
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0400000","111116406f032230000","0163103100000010"); // std
    cuts.AddCutPCMCalo("00052113","00200009f9730000dge0400000","111116406f032230000","0163103100000010"); // std
    cuts.AddCutPCMCalo("00081113","00200009f9730000dge0400000","111116406f032230000","0163103100000010"); // std
  } else if ( trainConfig == 94){ // EMCAL clusters 8 TeV LHC12, special TB 100 MeV
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0400000","111116906f032230000","0163103100000010"); // std
    cuts.AddCutPCMCalo("00052113","00200009f9730000dge0400000","111116906f032230000","0163103100000010"); // std
    cuts.AddCutPCMCalo("00081113","00200009f9730000dge0400000","111116906f032230000","0163103100000010"); // std

  } else if ( trainConfig == 100){ // EMCAL clusters 8 TeV LHC12 std cuts
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111111067032230000","0163103100000010"); // std
    cuts.AddCutPCMCalo("00052113","00200009327000008250400000","1111111067032230000","0163103100000010"); // std
    cuts.AddCutPCMCalo("00081113","00200009327000008250400000","1111111067032230000","0163103100000010"); // std
  } else if ( trainConfig == 101){ // EMCAL clusters 8 TeV LHC12
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111111067032230000","0163103100000010"); // std
    cuts.AddCutPCMCalo("00052113","00200009327000008250400000","1111111067032230000","0163103100000010"); // std
    cuts.AddCutPCMCalo("00081113","00200009327000008250400000","1111111067032230000","0163103100000010"); // std
  } else if ( trainConfig == 102){ //EMCAL minEnergy variation
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111111067022230000","0163103100000010"); //0.6 GeV/c
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111111067042230000","0163103100000010"); //0.8 GeV/c
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111111067052230000","0163103100000010"); //0.9 GeV/c
  } else if ( trainConfig == 103){ //EMCAL minNCells variation
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111111067031230000","0163103100000010"); //n cells >= 1
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111111067033230000","0163103100000010"); //n cells >= 3
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111111067032200000","0163103100000010"); //no max M02 cut
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111111067032250000","0163103100000010"); //M02 < 0.3
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111111067032260000","0163103100000010"); //M02 < 0.27
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","11111110670322k0000","0163103100000010"); //M02, pT-dep
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1112111067032230000","0163103100000010"); //only modules with TRD infront
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111311067032230000","0163103100000010"); //no modules with TRD infront
  } else if ( trainConfig == 104){ // EMCAL track matching variations
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111111066032230000","0163103100000010"); // track matching variations
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111111068032230000","0163103100000010"); //
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111111069032230000","0163103100000010"); //
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111111063032230000","0163103100000010"); //
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111111060032230000","0163103100000010"); //
  } else if ( trainConfig == 105){ // EMCAL clusters, timing variation
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111111057032230000","0163103100000010"); // time 50ns
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111111047032230000","0163103100000010"); // time 100ns
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111111037032230000","0163103100000010"); // time 200ns
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111111027032230000","0163103100000010"); // time 500ns
  } else if ( trainConfig == 106){ // EMCAL clusters, timing variation
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111111077032230000","0163103100000010"); // time
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111111087032230000","0163103100000010"); // time
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111111097032230000","0163103100000010"); // time
  } else if ( trainConfig == 107){ // EMCAL clusters, exotic cut var
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111111067232230000","0163103100000010"); //
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111111067332230000","0163103100000010"); //
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111111067532230000","0163103100000010"); //
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111111067732230000","0163103100000010"); //
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111111067932230000","0163103100000010"); //
  } else if ( trainConfig == 108){ // EMCAL track matching variations
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111111062032230000","0163103100000010"); //
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111111063032230000","0163103100000010"); //
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111111064032230000","0163103100000010"); //
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111111065032230000","0163103100000010"); //
  } else if ( trainConfig == 109){  //Different NonLinearities part2
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111101067032230000","0163103100000010"); // NonLinearity kSDMv5
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111113067032230000","0163103100000010"); // NonLinearity LHC12 kTestBeamv2+ConvCalo
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111114067032230000","0163103100000010"); // NonLinearity LHC12 kTestBeamv2+Calo
  } else if ( trainConfig == 110){ // Different NonLinearities
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111111067032230000","0163103100000010"); // NonLinearity LHC12 ConvCalo
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111112067032230000","0163103100000010"); // NonLinearity LHC12 Calo
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111121067032230000","0163103100000010"); // NonLinearity LHC12 ConvCalo MassRatioFits
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111122067032230000","0163103100000010"); // NonLinearity LHC12 Calo MassRatioFits
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111100067032230000","0163103100000010"); // NonLinearity none
  } else if ( trainConfig == 111){ // No NonLinearities
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111100067032230000","0163103100000010"); //
    cuts.AddCutPCMCalo("00052113","00200009327000008250400000","1111100067032230000","0163103100000010"); //
    cuts.AddCutPCMCalo("00081113","00200009327000008250400000","1111100067032230000","0163103100000010"); //
  } else if ( trainConfig == 112){ // Variations DistanceToBadChannel
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111111167032230000","0163103100000010"); //
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111111267032230000","0163103100000010"); //
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111111367032230000","0163103100000010"); //
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111111567032230000","0163103100000010"); //
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111111667032230000","0163103100000010"); //
  } else if ( trainConfig == 113){ // PCM variations
    cuts.AddCutPCMCalo("00010113","00200009227000008250400000","1111111067032230000","0163103100000010"); // dEdx e -3, 5
    cuts.AddCutPCMCalo("00010113","00200009127000008250400000","1111111067032230000","0163103100000010"); // dEdx e -5, 5
    cuts.AddCutPCMCalo("00010113","00200009357000008250400000","1111111067032230000","0163103100000010"); // dEdx pi 2
    //cuts.AddCutPCMCalo("00010113","00200009317000008250400000","1111111067032230000","0163103100000010"); // dEdx pi 0
    //cuts.AddCutPCMCalo("00010113","00200009387300008250400000","1111111067032230000","0163103100000010"); // dEdx pi 2 high 1 (> 3.5 GeV)
  } else if ( trainConfig == 114){ // PCM variations pi dEdx
    cuts.AddCutPCMCalo("00010113","00200009317300008250400000","1111111067032230000","0163103100000010"); // dEdx pi: 0: 0.4-3.5, -10: 3.5 ->
    cuts.AddCutPCMCalo("00010113","00200009327300008250400000","1111111067032230000","0163103100000010"); // dEdx pi: 1: 0.4-3.5, -10: 3.5 ->
    cuts.AddCutPCMCalo("00010113","00200009325000008250400000","1111111067032230000","0163103100000010"); // dEdx pi: 1: 0.3 ->
    cuts.AddCutPCMCalo("00010113","00200009320000008250400000","1111111067032230000","0163103100000010"); // dEdx pi: 1: 0.5 ->
  } else if ( trainConfig == 115){ // PCM variations pi dEdx
    cuts.AddCutPCMCalo("00010113","00200009327600008250400000","1111111067032230000","0163103100000010"); // dEdx pi: 1: 0.4-2, -10: 2. ->
    cuts.AddCutPCMCalo("00010113","00200009327400008250400000","1111111067032230000","0163103100000010"); // dEdx pi: 1: 0.4-3, -10: 3. ->
    cuts.AddCutPCMCalo("00010113","00200009367400008250400000","1111111067032230000","0163103100000010"); // dEdx pi: 2: 0.4-3, 0.5: 3. ->
    //cuts.AddCutPCMCalo("00010113","00200009347400008250400000","1111111067032230000","0163103100000010"); // dEdx pi: 3: 0.4-3, 1: 3. ->
    //cuts.AddCutPCMCalo("00010113","00200009315600008250400000","1111111067032230000","0163103100000010"); // dEdx pi: 0: 0.3-2, -10: 2. ->
  } else if ( trainConfig == 116){ // PCM variations
    cuts.AddCutPCMCalo("00010113","00200009327000009250400000","1111111067032230000","0163103100000010"); // qt 2D 0.03
    cuts.AddCutPCMCalo("00010113","00200009327000003250400000","1111111067032230000","0163103100000010"); // qt 1D 0.05
    cuts.AddCutPCMCalo("00010113","00200009327000002250400000","1111111067032230000","0163103100000010"); // qt 1D 0.07
    cuts.AddCutPCMCalo("00010113","00200049327000008250400000","1111111067032230000","0163103100000010"); // single pt > 0.075
    cuts.AddCutPCMCalo("00010113","00200019327000008250400000","1111111067032230000","0163103100000010"); // single pt > 0.1
  } else if ( trainConfig == 117){ // PCM variations
    cuts.AddCutPCMCalo("00010113","00200009327000008850400000","1111111067032230000","0163103100000010"); // 2D psi pair chi2 var
    cuts.AddCutPCMCalo("00010113","00200009327000008260400000","1111111067032230000","0163103100000010"); // 2D psi pair chi2 var
    cuts.AddCutPCMCalo("00010113","00200009327000008860400000","1111111067032230000","0163103100000010"); // 2D psi pair chi2 var
    cuts.AddCutPCMCalo("00010113","00200009327000008280400000","1111111067032230000","0163103100000010"); // 2D psi pair chi2 var
    cuts.AddCutPCMCalo("00010113","00200009327000008880400000","1111111067032230000","0163103100000010"); // 2D psi pair chi2 var
  } else if ( trainConfig == 118){ // PCM variations
    cuts.AddCutPCMCalo("00010113","00200006327000008250400000","1111111067032230000","0163103100000010"); // min TPC cl > 0.7
    cuts.AddCutPCMCalo("00010113","00200008327000008250400000","1111111067032230000","0163103100000010"); // min TPC cl > 0.35
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111111067032230000","0163106100000010"); // alpha < 0.8
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111111067032230000","0163105100000010"); // alpha < 0.75
  } else if ( trainConfig == 119){ // PCM variations
    cuts.AddCutPCMCalo("00010113","00202209327000008250400000","1111111067032230000","0163103100000010"); // restrict acceptance to EMCAL loose
    cuts.AddCutPCMCalo("00010113","00204409327000008250400000","1111111067032230000","0163103100000010"); // restrict acceptance to EMCAL tight
  } else if ( trainConfig == 120){ // PCM variations to close V0s
    cuts.AddCutPCMCalo("00010113","00200009327000008250401000","1111111067032230000","0163103100000010"); //
    cuts.AddCutPCMCalo("00010113","00200009327000008250402000","1111111067032230000","0163103100000010"); //
    cuts.AddCutPCMCalo("00010113","00200009327000008250403000","1111111067032230000","0163103100000010"); //

  // std cuts with pT dep matching
  } else if ( trainConfig == 121){ // EMCAL clusters 8 TeV LHC12 - no SPD PileUp
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111111067032230000","0163103100000010"); // std
    cuts.AddCutPCMCalo("00010013","00200009327000008250400000","1111111067032230000","0163103100000010"); // std - no pileup cut
  } else if ( trainConfig == 122){ // EMCAL clusters 8 TeV LHC12 - no SPD PileUp
    cuts.AddCutPCMCalo("00052113","00200009327000008250400000","1111111067032230000","0163103100000010"); // std
    cuts.AddCutPCMCalo("00052013","00200009327000008250400000","1111111067032230000","0163103100000010"); // std - no pileup cut
  } else if ( trainConfig == 123){ // EMCAL clusters 8 TeV LHC12 - no SPD PileUp
    cuts.AddCutPCMCalo("00081113","00200009327000008250400000","1111111067032230000","0163103100000010"); // std
    cuts.AddCutPCMCalo("00081013","00200009327000008250400000","1111111067032230000","0163103100000010"); // std - no pileup cut

  } else if ( trainConfig == 124){ // EMCAL clusters pp 8 TeV
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111111067032230000","0163105100000010"); // std
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111111067032230000","0163106100000010"); // std/

  // std cuts with fix track matching
  } else if ( trainConfig == 125){ // EMCAL clusters 8 TeV LHC12
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111111063032230000","0163103100000010"); //
    cuts.AddCutPCMCalo("00052113","00200009327000008250400000","1111111063032230000","0163103100000010"); //
    cuts.AddCutPCMCalo("00081113","00200009327000008250400000","1111111063032230000","0163103100000010"); //

  } else if ( trainConfig == 126){ // EMCAL clusters 8 TeV LHC12 std cuts
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0400000","111110106f032230000","0163103100000010"); // std
    cuts.AddCutPCMCalo("00052113","00200009f9730000dge0400000","111110106f032230000","0163103100000010"); // std
    cuts.AddCutPCMCalo("00081113","00200009f9730000dge0400000","111110106f032230000","0163103100000010"); // std

  } else if ( trainConfig == 127){ // EMCAL clusters 8 TeV LHC12, TB+finetuning CCRF
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0400000","111113106f032230000","0163103100000010"); // std
    cuts.AddCutPCMCalo("00052113","00200009f9730000dge0400000","111113106f032230000","0163103100000010"); // std
    cuts.AddCutPCMCalo("00081113","00200009f9730000dge0400000","111113106f032230000","0163103100000010"); // std
  } else if ( trainConfig == 128){ // EMCAL clusters 8 TeV LHC12, TB+finetuning CRF
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0400000","111113206f032230000","0163103100000010"); // std
    cuts.AddCutPCMCalo("00052113","00200009f9730000dge0400000","111113206f032230000","0163103100000010"); // std
    cuts.AddCutPCMCalo("00081113","00200009f9730000dge0400000","111113206f032230000","0163103100000010"); // std
  } else if ( trainConfig == 129){ // EMCAL clusters 8 TeV LHC12, TB+finetuning CCMF
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0400000","111113306f032230000","0163103100000010"); // std
    cuts.AddCutPCMCalo("00052113","00200009f9730000dge0400000","111113306f032230000","0163103100000010"); // std
    cuts.AddCutPCMCalo("00081113","00200009f9730000dge0400000","111113306f032230000","0163103100000010"); // std
  } else if ( trainConfig == 130){ // EMCAL clusters 8 TeV LHC12, TB+finetuning CMF
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0400000","111113406f032230000","0163103100000010"); // std
    cuts.AddCutPCMCalo("00052113","00200009f9730000dge0400000","111113406f032230000","0163103100000010"); // std
    cuts.AddCutPCMCalo("00081113","00200009f9730000dge0400000","111113406f032230000","0163103100000010"); // std

  // only std cuts
  } else if ( trainConfig == 131){ // EMCAL clusters 8 TeV LHC12
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111111067032230000","0163103100000010"); // std

  //kEMC7
  } else if ( trainConfig == 132){ //EMCAL minEnergy variation
    cuts.AddCutPCMCalo("00052113","00200009327000008250400000","1111111067022230000","0163103100000010"); //0.6 GeV/c
    cuts.AddCutPCMCalo("00052113","00200009327000008250400000","1111111067042230000","0163103100000010"); //0.8 GeV/c
    cuts.AddCutPCMCalo("00052113","00200009327000008250400000","1111111067052230000","0163103100000010"); //0.9 GeV/c
  } else if ( trainConfig == 133){ //EMCAL minNCells variation
    cuts.AddCutPCMCalo("00052113","00200009327000008250400000","1111111067031230000","0163103100000010"); //n cells >= 1
    cuts.AddCutPCMCalo("00052113","00200009327000008250400000","1111111067033230000","0163103100000010"); //n cells >= 3
    cuts.AddCutPCMCalo("00052113","00200009327000008250400000","1111111067032200000","0163103100000010"); //no max M02 cut
    cuts.AddCutPCMCalo("00052113","00200009327000008250400000","1111111067032250000","0163103100000010"); //M02 < 0.3
    cuts.AddCutPCMCalo("00052113","00200009327000008250400000","1111111067032260000","0163103100000010"); //M02 < 0.27
    cuts.AddCutPCMCalo("00052113","00200009327000008250400000","1112111067032230000","0163103100000010"); //only modules with TRD infront
    cuts.AddCutPCMCalo("00052113","00200009327000008250400000","1111311067032230000","0163103100000010"); //no modules with TRD infront
  } else if ( trainConfig == 134){ // EMCAL track matching variations
    cuts.AddCutPCMCalo("00052113","00200009327000008250400000","1111111066032230000","0163103100000010"); // track matching variations
    cuts.AddCutPCMCalo("00052113","00200009327000008250400000","1111111068032230000","0163103100000010"); //
    cuts.AddCutPCMCalo("00052113","00200009327000008250400000","1111111069032230000","0163103100000010"); //
    cuts.AddCutPCMCalo("00052113","00200009327000008250400000","1111111063032230000","0163103100000010"); //
    cuts.AddCutPCMCalo("00052113","00200009327000008250400000","1111111060032230000","0163103100000010"); //
  } else if ( trainConfig == 135){ // EMCAL clusters, timing variation
    cuts.AddCutPCMCalo("00052113","00200009327000008250400000","1111111057032230000","0163103100000010"); // time 50ns
    cuts.AddCutPCMCalo("00052113","00200009327000008250400000","1111111047032230000","0163103100000010"); // time 100ns
    cuts.AddCutPCMCalo("00052113","00200009327000008250400000","1111111037032230000","0163103100000010"); // time 200ns
    cuts.AddCutPCMCalo("00052113","00200009327000008250400000","1111111027032230000","0163103100000010"); // time 500ns
  } else if ( trainConfig == 136){ // EMCAL clusters, timing variation
    cuts.AddCutPCMCalo("00052113","00200009327000008250400000","1111111077032230000","0163103100000010"); // time
    cuts.AddCutPCMCalo("00052113","00200009327000008250400000","1111111087032230000","0163103100000010"); // time
    cuts.AddCutPCMCalo("00052113","00200009327000008250400000","1111111097032230000","0163103100000010"); // time
  } else if ( trainConfig == 137){ // EMCAL clusters, exotic cut var
    cuts.AddCutPCMCalo("00052113","00200009327000008250400000","1111111067232230000","0163103100000010"); //
    cuts.AddCutPCMCalo("00052113","00200009327000008250400000","1111111067332230000","0163103100000010"); //
    cuts.AddCutPCMCalo("00052113","00200009327000008250400000","1111111067532230000","0163103100000010"); //
    cuts.AddCutPCMCalo("00052113","00200009327000008250400000","1111111067732230000","0163103100000010"); //
    cuts.AddCutPCMCalo("00052113","00200009327000008250400000","1111111067932230000","0163103100000010"); //
  } else if ( trainConfig == 138){ // EMCAL track matching variations
    cuts.AddCutPCMCalo("00052113","00200009327000008250400000","1111111062032230000","0163103100000010"); //
    cuts.AddCutPCMCalo("00052113","00200009327000008250400000","1111111063032230000","0163103100000010"); //
    cuts.AddCutPCMCalo("00052113","00200009327000008250400000","1111111064032230000","0163103100000010"); //
  } else if ( trainConfig == 139){  //Different NonLinearities part2
    cuts.AddCutPCMCalo("00052113","00200009327000008250400000","1111101067032230000","0163103100000010"); // NonLinearity kSDMv5
    cuts.AddCutPCMCalo("00052113","00200009327000008250400000","1111113067032230000","0163103100000010"); // NonLinearity LHC12 kTestBeamv2+ConvCalo
    cuts.AddCutPCMCalo("00052113","00200009327000008250400000","1111114067032230000","0163103100000010"); // NonLinearity LHC12 kTestBeamv2+Calo
  } else if ( trainConfig == 140){ // Different NonLinearities
    cuts.AddCutPCMCalo("00052113","00200009327000008250400000","1111111067032230000","0163103100000010"); // NonLinearity LHC12 ConvCalo
    cuts.AddCutPCMCalo("00052113","00200009327000008250400000","1111112067032230000","0163103100000010"); // NonLinearity LHC12 Calo
    cuts.AddCutPCMCalo("00052113","00200009327000008250400000","1111121067032230000","0163103100000010"); // NonLinearity LHC12 ConvCalo MassRatioFits
    cuts.AddCutPCMCalo("00052113","00200009327000008250400000","1111122067032230000","0163103100000010"); // NonLinearity LHC12 Calo MassRatioFits
    cuts.AddCutPCMCalo("00052113","00200009327000008250400000","1111100067032230000","0163103100000010"); // NonLinearity none

  } else if ( trainConfig == 142){ // Variations DistanceToBadChannel
    cuts.AddCutPCMCalo("00052113","00200009327000008250400000","1111111167032230000","0163103100000010"); //
    cuts.AddCutPCMCalo("00052113","00200009327000008250400000","1111111267032230000","0163103100000010"); //
    cuts.AddCutPCMCalo("00052113","00200009327000008250400000","1111111367032230000","0163103100000010"); //
    cuts.AddCutPCMCalo("00052113","00200009327000008250400000","1111111567032230000","0163103100000010"); //
    cuts.AddCutPCMCalo("00052113","00200009327000008250400000","1111111667032230000","0163103100000010"); //
  } else if ( trainConfig == 143){ // PCM variations
    cuts.AddCutPCMCalo("00052113","00200009227000008250400000","1111111067032230000","0163103100000010"); // dEdx e -3, 5
    cuts.AddCutPCMCalo("00052113","00200009127000008250400000","1111111067032230000","0163103100000010"); // dEdx e -5, 5
    cuts.AddCutPCMCalo("00052113","00200009357000008250400000","1111111067032230000","0163103100000010"); // dEdx pi 2
    cuts.AddCutPCMCalo("00052113","00200009317000008250400000","1111111067032230000","0163103100000010"); // dEdx pi 0
    cuts.AddCutPCMCalo("00052113","00200009387300008250400000","1111111067032230000","0163103100000010"); // dEdx pi 2 high 1 (> 3.5 GeV)
  } else if ( trainConfig == 144){ // PCM variations pi dEdx
    cuts.AddCutPCMCalo("00052113","00200009317300008250400000","1111111067032230000","0163103100000010"); // dEdx pi: 0: 0.4-3.5, -10: 3.5 ->
    cuts.AddCutPCMCalo("00052113","00200009327300008250400000","1111111067032230000","0163103100000010"); // dEdx pi: 1: 0.4-3.5, -10: 3.5 ->
    cuts.AddCutPCMCalo("00052113","00200009325000008250400000","1111111067032230000","0163103100000010"); // dEdx pi: 1: 0.3 ->
    cuts.AddCutPCMCalo("00052113","00200009320000008250400000","1111111067032230000","0163103100000010"); // dEdx pi: 1: 0.5 ->
  } else if ( trainConfig == 145){ // PCM variations pi dEdx
    cuts.AddCutPCMCalo("00052113","00200009327600008250400000","1111111067032230000","0163103100000010"); // dEdx pi: 1: 0.4-2, -10: 2. ->
    cuts.AddCutPCMCalo("00052113","00200009327400008250400000","1111111067032230000","0163103100000010"); // dEdx pi: 1: 0.4-3, -10: 3. ->
    cuts.AddCutPCMCalo("00052113","00200009315600008250400000","1111111067032230000","0163103100000010"); // dEdx pi: 0: 0.3-2, -10: 2. ->
    cuts.AddCutPCMCalo("00052113","00200009367400008250400000","1111111067032230000","0163103100000010"); // dEdx pi: 2: 0.4-3, 0.5: 3. ->
    cuts.AddCutPCMCalo("00052113","00200009347400008250400000","1111111067032230000","0163103100000010"); // dEdx pi: 3: 0.4-3, 1: 3. ->
  } else if ( trainConfig == 146){ // PCM variations
    cuts.AddCutPCMCalo("00052113","00200009327000009250400000","1111111067032230000","0163103100000010"); // qt 2D 0.03
    cuts.AddCutPCMCalo("00052113","00200009327000003250400000","1111111067032230000","0163103100000010"); // qt 1D 0.05
    cuts.AddCutPCMCalo("00052113","00200009327000002250400000","1111111067032230000","0163103100000010"); // qt 1D 0.07
    cuts.AddCutPCMCalo("00052113","00200049327000008250400000","1111111067032230000","0163103100000010"); // single pt > 0.075
    cuts.AddCutPCMCalo("00052113","00200019327000008250400000","1111111067032230000","0163103100000010"); // single pt > 0.1
  } else if ( trainConfig == 147){ // PCM variations
    cuts.AddCutPCMCalo("00052113","00200009327000008850400000","1111111067032230000","0163103100000010"); // 2D psi pair chi2 var
    cuts.AddCutPCMCalo("00052113","00200009327000008260400000","1111111067032230000","0163103100000010"); // 2D psi pair chi2 var
    cuts.AddCutPCMCalo("00052113","00200009327000008860400000","1111111067032230000","0163103100000010"); // 2D psi pair chi2 var
    cuts.AddCutPCMCalo("00052113","00200009327000008280400000","1111111067032230000","0163103100000010"); // 2D psi pair chi2 var
    cuts.AddCutPCMCalo("00052113","00200009327000008880400000","1111111067032230000","0163103100000010"); // 2D psi pair chi2 var
  } else if ( trainConfig == 148){ // PCM variations
    cuts.AddCutPCMCalo("00052113","00200006327000008250400000","1111111067032230000","0163103100000010"); // min TPC cl > 0.7
    cuts.AddCutPCMCalo("00052113","00200008327000008250400000","1111111067032230000","0163103100000010"); // min TPC cl > 0.35
    cuts.AddCutPCMCalo("00052113","00200009327000008250400000","1111111067032230000","0163106100000010"); // alpha < 0.8
    cuts.AddCutPCMCalo("00052113","00200009327000008250400000","1111111067032230000","0163105100000010"); // alpha < 0.75
  } else if ( trainConfig == 149){ // PCM variations
    cuts.AddCutPCMCalo("00052113","00202209327000008250400000","1111111067032230000","0163103100000010"); // restrict acceptance to EMCAL loose
    cuts.AddCutPCMCalo("00052113","00204409327000008250400000","1111111067032230000","0163103100000010"); // restrict acceptance to EMCAL tight
  } else if ( trainConfig == 150){ // PCM variations to close V0s
    cuts.AddCutPCMCalo("00052113","00200009327000008250401000","1111111067032230000","0163103100000010"); //
    cuts.AddCutPCMCalo("00052113","00200009327000008250402000","1111111067032230000","0163103100000010"); //
    cuts.AddCutPCMCalo("00052113","00200009327000008250403000","1111111067032230000","0163103100000010"); //
  } else if ( trainConfig == 151){ // EMCAL clusters pp 8 TeV
    cuts.AddCutPCMCalo("00052113","00200009327000008250400000","1111111067032230000","0163105100000010"); // std
    cuts.AddCutPCMCalo("00052113","00200009327000008250400000","1111111067032230000","0163106100000010"); // std/

  // only std cuts
  } else if ( trainConfig == 158){ //std EMC7
    cuts.AddCutPCMCalo("00052113","00200009327000008250400000","1111111067032230000","0163103100000010"); // only EMC7
  } else if ( trainConfig == 159){ //std EMC7
    cuts.AddCutPCMCalo("00052113","00200009327000008250400000","1111111067032230000","0163103100000010"); // only EMC7

  //kEMCEGA
  } else if ( trainConfig == 162){ //EMCAL minEnergy variation
    cuts.AddCutPCMCalo("00081113","00200009327000008250400000","1111111067022230000","0163103100000010"); //0.6 GeV/c
    cuts.AddCutPCMCalo("00081113","00200009327000008250400000","1111111067042230000","0163103100000010"); //0.8 GeV/c
    cuts.AddCutPCMCalo("00081113","00200009327000008250400000","1111111067052230000","0163103100000010"); //0.9 GeV/c
  } else if ( trainConfig == 163){ //EMCAL minNCells variation
    cuts.AddCutPCMCalo("00081113","00200009327000008250400000","1111111067031230000","0163103100000010"); //n cells >= 1
    cuts.AddCutPCMCalo("00081113","00200009327000008250400000","1111111067033230000","0163103100000010"); //n cells >= 3
    cuts.AddCutPCMCalo("00081113","00200009327000008250400000","1111111067032200000","0163103100000010"); //no max M02 cut
    cuts.AddCutPCMCalo("00081113","00200009327000008250400000","1111111067032250000","0163103100000010"); //M02 < 0.3
    cuts.AddCutPCMCalo("00081113","00200009327000008250400000","1111111067032260000","0163103100000010"); //M02 < 0.27
    cuts.AddCutPCMCalo("00081113","00200009327000008250400000","1112111067032230000","0163103100000010"); //only modules with TRD infront
    cuts.AddCutPCMCalo("00081113","00200009327000008250400000","1111311067032230000","0163103100000010"); //no modules with TRD infront
  } else if ( trainConfig == 164){ // EMCAL track matching variations
    cuts.AddCutPCMCalo("00081113","00200009327000008250400000","1111111066032230000","0163103100000010"); // track matching variations
    cuts.AddCutPCMCalo("00081113","00200009327000008250400000","1111111068032230000","0163103100000010"); //
    cuts.AddCutPCMCalo("00081113","00200009327000008250400000","1111111069032230000","0163103100000010"); //
    cuts.AddCutPCMCalo("00081113","00200009327000008250400000","1111111063032230000","0163103100000010"); //
    cuts.AddCutPCMCalo("00081113","00200009327000008250400000","1111111060032230000","0163103100000010"); //
  } else if ( trainConfig == 165){ // EMCAL clusters, timing variation
    cuts.AddCutPCMCalo("00081113","00200009327000008250400000","1111111057032230000","0163103100000010"); // time 50ns
    cuts.AddCutPCMCalo("00081113","00200009327000008250400000","1111111047032230000","0163103100000010"); // time 100ns
    cuts.AddCutPCMCalo("00081113","00200009327000008250400000","1111111037032230000","0163103100000010"); // time 200ns
    cuts.AddCutPCMCalo("00081113","00200009327000008250400000","1111111027032230000","0163103100000010"); // time 500ns
  } else if ( trainConfig == 166){ // EMCAL clusters, timing variation
    cuts.AddCutPCMCalo("00081113","00200009327000008250400000","1111111077032230000","0163103100000010"); // time
    cuts.AddCutPCMCalo("00081113","00200009327000008250400000","1111111087032230000","0163103100000010"); // time
    cuts.AddCutPCMCalo("00081113","00200009327000008250400000","1111111097032230000","0163103100000010"); // time
  } else if ( trainConfig == 167){ // EMCAL clusters, exotic cut var
    cuts.AddCutPCMCalo("00081113","00200009327000008250400000","1111111067232230000","0163103100000010"); //
    cuts.AddCutPCMCalo("00081113","00200009327000008250400000","1111111067332230000","0163103100000010"); //
    cuts.AddCutPCMCalo("00081113","00200009327000008250400000","1111111067532230000","0163103100000010"); //
    cuts.AddCutPCMCalo("00081113","00200009327000008250400000","1111111067732230000","0163103100000010"); //
    cuts.AddCutPCMCalo("00081113","00200009327000008250400000","1111111067932230000","0163103100000010"); //
  } else if ( trainConfig == 168){ // EMCAL track matching variations
    cuts.AddCutPCMCalo("00081113","00200009327000008250400000","1111111062032230000","0163103100000010"); //
    cuts.AddCutPCMCalo("00081113","00200009327000008250400000","1111111063032230000","0163103100000010"); //
    cuts.AddCutPCMCalo("00081113","00200009327000008250400000","1111111064032230000","0163103100000010"); //
  } else if ( trainConfig == 169){  //Different NonLinearities part2
    cuts.AddCutPCMCalo("00081113","00200009327000008250400000","1111101067032230000","0163103100000010"); // NonLinearity kSDMv5
    cuts.AddCutPCMCalo("00081113","00200009327000008250400000","1111113067032230000","0163103100000010"); // NonLinearity LHC12 kTestBeamv2+ConvCalo
    cuts.AddCutPCMCalo("00081113","00200009327000008250400000","1111114067032230000","0163103100000010"); // NonLinearity LHC12 kTestBeamv2+Calo
  } else if ( trainConfig == 170){ // Different NonLinearities
    cuts.AddCutPCMCalo("00081113","00200009327000008250400000","1111111067032230000","0163103100000010"); // NonLinearity LHC12 ConvCalo
    cuts.AddCutPCMCalo("00081113","00200009327000008250400000","1111112067032230000","0163103100000010"); // NonLinearity LHC12 Calo
    cuts.AddCutPCMCalo("00081113","00200009327000008250400000","1111121067032230000","0163103100000010"); // NonLinearity LHC12 ConvCalo MassRatioFits
    cuts.AddCutPCMCalo("00081113","00200009327000008250400000","1111122067032230000","0163103100000010"); // NonLinearity LHC12 Calo MassRatioFits
    cuts.AddCutPCMCalo("00081113","00200009327000008250400000","1111100067032230000","0163103100000010"); // NonLinearity none
  } else if ( trainConfig == 171){ // EMCAL clusters pp 8 TeV
    cuts.AddCutPCMCalo("00081113","00200009327000008250400000","1111111067032230000","0163105100000010"); // std
    cuts.AddCutPCMCalo("00081113","00200009327000008250400000","1111111067032230000","0163106100000010"); // std/

  } else if ( trainConfig == 172){ // Variations DistanceToBadChannel
    cuts.AddCutPCMCalo("00081113","00200009327000008250400000","1111111167032230000","0163103100000010"); //
    cuts.AddCutPCMCalo("00081113","00200009327000008250400000","1111111267032230000","0163103100000010"); //
    cuts.AddCutPCMCalo("00081113","00200009327000008250400000","1111111367032230000","0163103100000010"); //
    cuts.AddCutPCMCalo("00081113","00200009327000008250400000","1111111567032230000","0163103100000010"); //
    cuts.AddCutPCMCalo("00081113","00200009327000008250400000","1111111667032230000","0163103100000010"); //
  } else if ( trainConfig == 173){ // PCM variations
    cuts.AddCutPCMCalo("00081113","00200009227000008250400000","1111111067032230000","0163103100000010"); // dEdx e -3, 5
    cuts.AddCutPCMCalo("00081113","00200009127000008250400000","1111111067032230000","0163103100000010"); // dEdx e -5, 5
    cuts.AddCutPCMCalo("00081113","00200009357000008250400000","1111111067032230000","0163103100000010"); // dEdx pi 2
    cuts.AddCutPCMCalo("00081113","00200009317000008250400000","1111111067032230000","0163103100000010"); // dEdx pi 0
    cuts.AddCutPCMCalo("00081113","00200009387300008250400000","1111111067032230000","0163103100000010"); // dEdx pi 2 high 1 (> 3.5 GeV)
  } else if ( trainConfig == 174){ // PCM variations pi dEdx
    cuts.AddCutPCMCalo("00081113","00200009317300008250400000","1111111067032230000","0163103100000010"); // dEdx pi: 0: 0.4-3.5, -10: 3.5 ->
    cuts.AddCutPCMCalo("00081113","00200009327300008250400000","1111111067032230000","0163103100000010"); // dEdx pi: 1: 0.4-3.5, -10: 3.5 ->
    cuts.AddCutPCMCalo("00081113","00200009325000008250400000","1111111067032230000","0163103100000010"); // dEdx pi: 1: 0.3 ->
    cuts.AddCutPCMCalo("00081113","00200009320000008250400000","1111111067032230000","0163103100000010"); // dEdx pi: 1: 0.5 ->
  } else if ( trainConfig == 175){ // PCM variations pi dEdx
    cuts.AddCutPCMCalo("00081113","00200009327600008250400000","1111111067032230000","0163103100000010"); // dEdx pi: 1: 0.4-2, -10: 2. ->
    cuts.AddCutPCMCalo("00081113","00200009327400008250400000","1111111067032230000","0163103100000010"); // dEdx pi: 1: 0.4-3, -10: 3. ->
    cuts.AddCutPCMCalo("00081113","00200009315600008250400000","1111111067032230000","0163103100000010"); // dEdx pi: 0: 0.3-2, -10: 2. ->
    cuts.AddCutPCMCalo("00081113","00200009367400008250400000","1111111067032230000","0163103100000010"); // dEdx pi: 2: 0.4-3, 0.5: 3. ->
    cuts.AddCutPCMCalo("00081113","00200009347400008250400000","1111111067032230000","0163103100000010"); // dEdx pi: 3: 0.4-3, 1: 3. ->
  } else if ( trainConfig == 176){ // PCM variations
    cuts.AddCutPCMCalo("00081113","00200009327000009250400000","1111111067032230000","0163103100000010"); // qt 2D 0.03
    cuts.AddCutPCMCalo("00081113","00200009327000003250400000","1111111067032230000","0163103100000010"); // qt 1D 0.05
    cuts.AddCutPCMCalo("00081113","00200009327000002250400000","1111111067032230000","0163103100000010"); // qt 1D 0.07
    cuts.AddCutPCMCalo("00081113","00200049327000008250400000","1111111067032230000","0163103100000010"); // single pt > 0.075
    cuts.AddCutPCMCalo("00081113","00200019327000008250400000","1111111067032230000","0163103100000010"); // single pt > 0.1
  } else if ( trainConfig == 177){ // PCM variations
    cuts.AddCutPCMCalo("00081113","00200009327000008850400000","1111111067032230000","0163103100000010"); // 2D psi pair chi2 var
    cuts.AddCutPCMCalo("00081113","00200009327000008260400000","1111111067032230000","0163103100000010"); // 2D psi pair chi2 var
    cuts.AddCutPCMCalo("00081113","00200009327000008860400000","1111111067032230000","0163103100000010"); // 2D psi pair chi2 var
    cuts.AddCutPCMCalo("00081113","00200009327000008280400000","1111111067032230000","0163103100000010"); // 2D psi pair chi2 var
    cuts.AddCutPCMCalo("00081113","00200009327000008880400000","1111111067032230000","0163103100000010"); // 2D psi pair chi2 var
  } else if ( trainConfig == 178){ // PCM variations
    cuts.AddCutPCMCalo("00081113","00200006327000008250400000","1111111067032230000","0163103100000010"); // min TPC cl > 0.7
    cuts.AddCutPCMCalo("00081113","00200008327000008250400000","1111111067032230000","0163103100000010"); // min TPC cl > 0.35
    cuts.AddCutPCMCalo("00081113","00200009327000008250400000","1111111067032230000","0163106100000010"); // alpha < 0.8
    cuts.AddCutPCMCalo("00081113","00200009327000008250400000","1111111067032230000","0163105100000010"); // alpha < 0.75
  } else if ( trainConfig == 179){ // PCM variations
    cuts.AddCutPCMCalo("00081113","00202209327000008250400000","1111111067032230000","0163103100000010"); // restrict acceptance to EMCAL loose
    cuts.AddCutPCMCalo("00081113","00204409327000008250400000","1111111067032230000","0163103100000010"); // restrict acceptance to EMCAL tight
  } else if ( trainConfig == 180){ // PCM variations to close V0s
    cuts.AddCutPCMCalo("00081113","00200009327000008250401000","1111111067032230000","0163103100000010"); //
    cuts.AddCutPCMCalo("00081113","00200009327000008250402000","1111111067032230000","0163103100000010"); //
    cuts.AddCutPCMCalo("00081113","00200009327000008250403000","1111111067032230000","0163103100000010"); //
  // only std cuts
  } else if ( trainConfig == 181){ //std EGA
    cuts.AddCutPCMCalo("00081113","00200009327000008250400000","1111111067032230000","0163103100000010"); // only EGA
  } else if ( trainConfig == 182){ //std EGA
    cuts.AddCutPCMCalo("00081113","00200009327000008250400000","1111111067032230000","0163103100000010"); // only EGA

  //multiple std cuts for different studies
  } else if ( trainConfig == 183){ // EMCAL clusters 8 TeV LHC12
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111111067032230000","0163103100000010"); // std
    cuts.AddCutPCMCalo("00052113","00200009327000008250400000","1111111067032230000","0163103100000010"); // std
    cuts.AddCutPCMCalo("00081113","00200009327000008250400000","1111111067032230000","0163103100000010"); // std
  } else if ( trainConfig == 184){ // EMCAL clusters 8 TeV LHC12
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111111067032230000","0163103100000010"); // std
    cuts.AddCutPCMCalo("00052113","00200009327000008250400000","1111111067032230000","0163103100000010"); // std
    cuts.AddCutPCMCalo("00081113","00200009327000008250400000","1111111067032230000","0163103100000010"); // std
  } else if ( trainConfig == 185){ // EMCAL clusters 8 TeV LHC12
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111111067032230000","0163103100000010"); // std
    cuts.AddCutPCMCalo("00052113","00200009327000008250400000","1111111067032230000","0163103100000010"); // std
    cuts.AddCutPCMCalo("00081113","00200009327000008250400000","1111111067032230000","0163103100000010"); // std
  } else if ( trainConfig == 186){ // EMCAL clusters 8 TeV LHC12
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111111067032230000","0163103100000010"); // std
    cuts.AddCutPCMCalo("00052113","00200009327000008250400000","1111111067032230000","0163103100000010"); // std
    cuts.AddCutPCMCalo("00081113","00200009327000008250400000","1111111067032230000","0163103100000010"); // std
  } else if ( trainConfig == 187){ // EMCAL clusters 8 TeV LHC12
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111111067032230000","0163103100000010"); // std
    cuts.AddCutPCMCalo("00052113","00200009327000008250400000","1111111067032230000","0163103100000010"); // std
    cuts.AddCutPCMCalo("00081113","00200009327000008250400000","1111111067032230000","0163103100000010"); // std
  } else if ( trainConfig == 188){ // EMCAL clusters 8 TeV LHC12
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111111067032230000","0163103100000010"); // std
    cuts.AddCutPCMCalo("00052113","00200009327000008250400000","1111111067032230000","0163103100000010"); // std
    cuts.AddCutPCMCalo("00081113","00200009327000008250400000","1111111067032230000","0163103100000010"); // std
  } else if ( trainConfig == 189){ // EMCAL clusters 8 TeV LHC12
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111111067032230000","0163103100000010"); // std
    cuts.AddCutPCMCalo("00052113","00200009327000008250400000","1111111067032230000","0163103100000010"); // std
    cuts.AddCutPCMCalo("00081113","00200009327000008250400000","1111111067032230000","0163103100000010"); // std
  } else if ( trainConfig == 190){ // EMCAL clusters 8 TeV LHC12
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111111067032230000","0163103100000010"); // std
    cuts.AddCutPCMCalo("00052113","00200009327000008250400000","1111111067032230000","0163103100000010"); // std
    cuts.AddCutPCMCalo("00081113","00200009327000008250400000","1111111067032230000","0163103100000010"); // std
  } else if ( trainConfig == 191){ // EMCAL clusters 8 TeV LHC12
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111111067032230000","0163103100000010"); // std
    cuts.AddCutPCMCalo("00052113","00200009327000008250400000","1111111067032230000","0163103100000010"); // std
    cuts.AddCutPCMCalo("00081113","00200009327000008250400000","1111111067032230000","0163103100000010"); // std
  } else if ( trainConfig == 192){ // EMCAL clusters 8 TeV LHC12
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111111067032230000","0163103100000010"); // std
    cuts.AddCutPCMCalo("00052113","00200009327000008250400000","1111111067032230000","0163103100000010"); // std
    cuts.AddCutPCMCalo("00081113","00200009327000008250400000","1111111067032230000","0163103100000010"); // std
  } else if ( trainConfig == 193){ // EMCAL clusters 8 TeV LHC12
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111111067032230000","0163103100000010"); // std
    cuts.AddCutPCMCalo("00052113","00200009327000008250400000","1111111067032230000","0163103100000010"); // std
    cuts.AddCutPCMCalo("00081113","00200009327000008250400000","1111111067032230000","0163103100000010"); // std

  //INT8
  } else if ( trainConfig == 194){ // EMCAL clusters, INT8
    cuts.AddCutPCMCalo("00011113","00200009327000008250400000","1111111067032230000","0163103100000010"); // INT8
    cuts.AddCutPCMCalo("00053113","00200009327000008250400000","1111111067032230000","0163103100000010"); // EMC8
    cuts.AddCutPCMCalo("00082113","00200009327000008250400000","1111111067032230000","0163103100000010"); // EGA+INT8

  // eta variations
  } else if ( trainConfig == 195){ // EMCAL clusters 8 TeV LHC12, |eta| < 0.7, y < 0.7
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1551111067032230000","0163203100000010"); //
    cuts.AddCutPCMCalo("00052113","00200009327000008250400000","1551111067032230000","0163203100000010"); //
    cuts.AddCutPCMCalo("00081113","00200009327000008250400000","1551111067032230000","0163203100000010"); //
  } else if ( trainConfig == 196){ // EMCAL clusters 8 TeV LHC12, |eta| < 0.3, y < 0.3
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1661111067032230000","0163703100000010"); //
    cuts.AddCutPCMCalo("00052113","00200009327000008250400000","1661111067032230000","0163703100000010"); //
    cuts.AddCutPCMCalo("00081113","00200009327000008250400000","1661111067032230000","0163703100000010"); //

  //*************************************************************************************************
  // 7 TeV EMC setup (LHC10x)
  //*************************************************************************************************
  } else if ( trainConfig == 200){ // EMCAL clusters pp 7 TeV, pT dep matching
    cuts.AddCutPCMCalo("00000113","00200009327000008250400000","11111110b7032230000","0163103100000010"); // std
  } else if ( trainConfig == 201){ // EMCAL clusters pp 7 TeV, pT dep matching
    cuts.AddCutPCMCalo("00000113","00200009327000008250400000","11111310b7032230000","0163103100000010"); // std TB NL
  } else if ( trainConfig == 202){ //EMCAL minEnergy variation
    cuts.AddCutPCMCalo("00000113","00200009327000008250400000","11111110b7022230000","0163103100000010"); //0.6 GeV/c
    cuts.AddCutPCMCalo("00000113","00200009327000008250400000","11111110b7032230000","0163103100000010"); //0.7 GeV/c default
    cuts.AddCutPCMCalo("00000113","00200009327000008250400000","11111110b7042230000","0163103100000010"); //0.8 GeV/c
    cuts.AddCutPCMCalo("00000113","00200009327000008250400000","11111110b7052230000","0163103100000010"); //0.9 GeV/c
  } else if ( trainConfig == 203){ //EMCAL minNCells variation
    cuts.AddCutPCMCalo("00000113","00200009327000008250400000","11111110b7031230000","0163103100000010"); //n cells >= 1
    cuts.AddCutPCMCalo("00000113","00200009327000008250400000","11111110b7033230000","0163103100000010"); //n cells >= 3
    cuts.AddCutPCMCalo("00000113","00200009327000008250400000","11111110b7032200000","0163103100000010"); //no max M02 cut
    cuts.AddCutPCMCalo("00000113","00200009327000008250400000","11111110b7032250000","0163103100000010"); //M02 < 0.3
    cuts.AddCutPCMCalo("00000113","00200009327000008250400000","11111110b7032260000","0163103100000010"); //M02 < 0.27
    cuts.AddCutPCMCalo("00000113","00200009327000008250400000","11131110b7032230000","0163103100000010"); //only modules with TRD infront
    cuts.AddCutPCMCalo("00000113","00200009327000008250400000","11112110b7032230000","0163103100000010"); //no modules with TRD infront
  } else if ( trainConfig == 204){ // EMCAL track matching variations
    cuts.AddCutPCMCalo("00000113","00200009327000008250400000","11111110b6032230000","0163103100000010"); // track matching variations
    cuts.AddCutPCMCalo("00000113","00200009327000008250400000","11111110b7032230000","0163103100000010"); //
    cuts.AddCutPCMCalo("00000113","00200009327000008250400000","11111110b8032230000","0163103100000010"); //
    cuts.AddCutPCMCalo("00000113","00200009327000008250400000","11111110b9032230000","0163103100000010"); //
    cuts.AddCutPCMCalo("00000113","00200009327000008250400000","11111110b3032230000","0163103100000010"); //
    cuts.AddCutPCMCalo("00000113","00200009327000008250400000","11111110b0032230000","0163103100000010"); //
  } else if ( trainConfig == 205){ // EMCAL clusters, timing variation
    cuts.AddCutPCMCalo("00000113","00200009327000008250400000","1111111057032230000","0163103100000010"); // time 50ns
    cuts.AddCutPCMCalo("00000113","00200009327000008250400000","1111111047032230000","0163103100000010"); // time 100ns
    cuts.AddCutPCMCalo("00000113","00200009327000008250400000","1111111037032230000","0163103100000010"); // time 200ns
    cuts.AddCutPCMCalo("00000113","00200009327000008250400000","1111111027032230000","0163103100000010"); // time 500ns
  } else if ( trainConfig == 206){ // EMCAL clusters, timing variation
    cuts.AddCutPCMCalo("00000113","00200009327000008250400000","1111111067032230000","0163103100000010"); // time
    cuts.AddCutPCMCalo("00000113","00200009327000008250400000","1111111077032230000","0163103100000010"); // time
    cuts.AddCutPCMCalo("00000113","00200009327000008250400000","1111111087032230000","0163103100000010"); // time
    cuts.AddCutPCMCalo("00000113","00200009327000008250400000","1111111097032230000","0163103100000010"); // time
  } else if ( trainConfig == 207){ // EMCAL track matching variations
    cuts.AddCutPCMCalo("00000113","00200009327000008250400000","11111110b7032230000","0163103100000010"); // track matching variations
    cuts.AddCutPCMCalo("00000113","00200009327000008250400000","11111110b2032230000","0163103100000010"); //
    cuts.AddCutPCMCalo("00000113","00200009327000008250400000","11111110b3032230000","0163103100000010"); //
    cuts.AddCutPCMCalo("00000113","00200009327000008250400000","11111110b4032230000","0163103100000010"); //
    cuts.AddCutPCMCalo("00000113","00200009327000008250400000","11111110b5032230000","0163103100000010"); //

  } else if ( trainConfig == 209){  //Different NonLinearities part2
    cuts.AddCutPCMCalo("00000113","00200009327000008250400000","11111010b7032230000","0163103100000010"); // NonLinearity kSDMv5
    cuts.AddCutPCMCalo("00000113","00200009327000008250400000","11111130b7032230000","0163103100000010"); // NonLinearity LHC12 kTestBeamv2+ConvCalo
    cuts.AddCutPCMCalo("00000113","00200009327000008250400000","11111140b7032230000","0163103100000010"); // NonLinearity LHC12 kTestBeamv2+Calo
  } else if ( trainConfig == 210){ // Different NonLinearities
    cuts.AddCutPCMCalo("00000113","00200009327000008250400000","11111110b7032230000","0163103100000010"); // NonLinearity LHC12 ConvCalo
    cuts.AddCutPCMCalo("00000113","00200009327000008250400000","11111120b7032230000","0163103100000010"); // NonLinearity LHC12 Calo
    cuts.AddCutPCMCalo("00000113","00200009327000008250400000","11111210b7032230000","0163103100000010"); // NonLinearity LHC12 ConvCalo MassRatioFits
    cuts.AddCutPCMCalo("00000113","00200009327000008250400000","11111220b7032230000","0163103100000010"); // NonLinearity LHC12 Calo MassRatioFits
    cuts.AddCutPCMCalo("00000113","00200009327000008250400000","11111000b7032230000","0163103100000010"); // NonLinearity none
  } else if ( trainConfig == 211){ // No NonLinearities
    cuts.AddCutPCMCalo("00000113","00200009327000008250400000","11111000b7032230000","0163103100000010"); //
  } else if ( trainConfig == 212){ // Variations DistanceToBadChannel
    cuts.AddCutPCMCalo("00000113","00200009327000008250400000","11111110b7032230000","0163103100000010"); //
    cuts.AddCutPCMCalo("00000113","00200009327000008250400000","11111111b7032230000","0163103100000010"); //
    cuts.AddCutPCMCalo("00000113","00200009327000008250400000","11111112b7032230000","0163103100000010"); //
    cuts.AddCutPCMCalo("00000113","00200009327000008250400000","11111113b7032230000","0163103100000010"); //
    cuts.AddCutPCMCalo("00000113","00200009327000008250400000","11111115b7032230000","0163103100000010"); //
    cuts.AddCutPCMCalo("00000113","00200009327000008250400000","11111116b7032230000","0163103100000010"); //
  } else if ( trainConfig == 213){ // PCM variations
    cuts.AddCutPCMCalo("00000113","00200009227000008250400000","11111110b7032230000","0163103100000010"); // dEdx e -3, 5
    cuts.AddCutPCMCalo("00000113","00200009127000008250400000","11111110b7032230000","0163103100000010"); // dEdx e -5, 5
    cuts.AddCutPCMCalo("00000113","00200009357000008250400000","11111110b7032230000","0163103100000010"); // dEdx pi 2
    cuts.AddCutPCMCalo("00000113","00200009317000008250400000","11111110b7032230000","0163103100000010"); // dEdx pi 0
    cuts.AddCutPCMCalo("00000113","00200009387300008250400000","11111110b7032230000","0163103100000010"); // dEdx pi 2 high 1 (> 3.5 GeV)
  } else if ( trainConfig == 214){ // PCM variations pi dEdx
    cuts.AddCutPCMCalo("00000113","00200009317300008250400000","11111110b7032230000","0163103100000010"); // dEdx pi: 0: 0.4-3.5, -10: 3.5 ->
    cuts.AddCutPCMCalo("00000113","00200009327300008250400000","11111110b7032230000","0163103100000010"); // dEdx pi: 1: 0.4-3.5, -10: 3.5 ->
    cuts.AddCutPCMCalo("00000113","00200009325000008250400000","11111110b7032230000","0163103100000010"); // dEdx pi: 1: 0.3 ->
    cuts.AddCutPCMCalo("00000113","00200009320000008250400000","11111110b7032230000","0163103100000010"); // dEdx pi: 1: 0.5 ->
  } else if ( trainConfig == 215){ // PCM variations pi dEdx
    cuts.AddCutPCMCalo("00000113","00200009327600008250400000","11111110b7032230000","0163103100000010"); // dEdx pi: 1: 0.4-2, -10: 2. ->
    cuts.AddCutPCMCalo("00000113","00200009327400008250400000","11111110b7032230000","0163103100000010"); // dEdx pi: 1: 0.4-3, -10: 3. ->
    cuts.AddCutPCMCalo("00000113","00200009315600008250400000","11111110b7032230000","0163103100000010"); // dEdx pi: 0: 0.3-2, -10: 2. ->
    cuts.AddCutPCMCalo("00000113","00200009367400008250400000","11111110b7032230000","0163103100000010"); // dEdx pi: 2: 0.4-3, 0.5: 3. ->
    cuts.AddCutPCMCalo("00000113","00200009347400008250400000","11111110b7032230000","0163103100000010"); // dEdx pi: 3: 0.4-3, 1: 3. ->
  } else if ( trainConfig == 216){ // PCM variations
    cuts.AddCutPCMCalo("00000113","00200009327000009250400000","11111110b7032230000","0163103100000010"); // qt 2D 0.03
    cuts.AddCutPCMCalo("00000113","00200009327000003250400000","11111110b7032230000","0163103100000010"); // qt 1D 0.05
    cuts.AddCutPCMCalo("00000113","00200009327000002250400000","11111110b7032230000","0163103100000010"); // qt 1D 0.07
    cuts.AddCutPCMCalo("00000113","00200049327000008250400000","11111110b7032230000","0163103100000010"); // single pt > 0.075
    cuts.AddCutPCMCalo("00000113","00200019327000008250400000","11111110b7032230000","0163103100000010"); // single pt > 0.1
  } else if ( trainConfig == 217){ // PCM variations
    cuts.AddCutPCMCalo("00000113","00200009327000008850400000","11111110b7032230000","0163103100000010"); // 2D psi pair chi2 var
    cuts.AddCutPCMCalo("00000113","00200009327000008260400000","11111110b7032230000","0163103100000010"); // 2D psi pair chi2 var
    cuts.AddCutPCMCalo("00000113","00200009327000008860400000","11111110b7032230000","0163103100000010"); // 2D psi pair chi2 var
    cuts.AddCutPCMCalo("00000113","00200009327000008280400000","11111110b7032230000","0163103100000010"); // 2D psi pair chi2 var
    cuts.AddCutPCMCalo("00000113","00200009327000008880400000","11111110b7032230000","0163103100000010"); // 2D psi pair chi2 var
  } else if ( trainConfig == 218){ // PCM variations
    cuts.AddCutPCMCalo("00000113","00200006327000008250400000","11111110b7032230000","0163103100000010"); // min TPC cl > 0.7
    cuts.AddCutPCMCalo("00000113","00200008327000008250400000","11111110b7032230000","0163103100000010"); // min TPC cl > 0.35
    cuts.AddCutPCMCalo("00000113","00200009327000008250400000","11111110b7032230000","0163106100000010"); // alpha < 0.8
    cuts.AddCutPCMCalo("00000113","00200009327000008250400000","11111110b7032230000","0163105100000010"); // alpha < 0.75
  } else if ( trainConfig == 219){ // PCM variations
    cuts.AddCutPCMCalo("00000113","00202209327000008250400000","11111110b7032230000","0163103100000010"); // restrict acceptance to EMCAL loose
    cuts.AddCutPCMCalo("00000113","00204409327000008250400000","11111110b7032230000","0163103100000010"); // restrict acceptance to EMCAL tight
  } else if ( trainConfig == 220){ // PCM variations to close V0s
    cuts.AddCutPCMCalo("00000113","00200009327000008250400000","11111110b7032230000","0163103100000010"); //
    cuts.AddCutPCMCalo("00000113","00200009327000008250401000","11111110b7032230000","0163103100000010"); //
    cuts.AddCutPCMCalo("00000113","00200009327000008250402000","11111110b7032230000","0163103100000010"); //
    cuts.AddCutPCMCalo("00000113","00200009327000008250403000","11111110b7032230000","0163103100000010"); //


  //*************************************************************************************************
  // 7 TeV EMC setup (LHC11x)
  //*************************************************************************************************
  } else if ( trainConfig == 221){ // EMCAL clusters pp 7 TeV, std matching
    cuts.AddCutPCMCalo("00000113","00200009327000008250400000","11111110b3032230000","0163103100000010"); // std

  } else if ( trainConfig == 222){ // EMCAL clusters pp 7 TeV, no SPD pileup
    cuts.AddCutPCMCalo("00000113","00200009327000008250400000","11111110b7032230000","0163103100000010"); // std
    cuts.AddCutPCMCalo("00000013","00200009327000008250400000","11111110b7032230000","0163103100000010"); // std - no SPD pileup
  } else if ( trainConfig == 223){ // EMCAL clusters pp 7 TeV
    cuts.AddCutPCMCalo("00000113","0dm0000922700000dge0404000","1111a3104f032230000","0163103100000010"); // std

  // LHC11cd configs V0OR and V0AND
  } else if ( trainConfig == 250){ // EMCAL clusters 7 TeV LHC11 TM on, +-30ns, std TM, no NL
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111100067032230000","0163103100000010"); // MBAND
    cuts.AddCutPCMCalo("00052113","00200009327000008250400000","1111100067032230000","0163103100000010"); // EMC7
  } else if ( trainConfig == 251){ // EMCAL clusters 7 TeV LHC10 MB only
    cuts.AddCutPCMCalo("00000113","00200009327000008250400000","11111000b7032230000","0163103100000010"); // MBAND
  // LHC11cd configs V0OR and V0AND with Nonlinearity
  } else if ( trainConfig == 260){ // EMCAL clusters 7 TeV LHC11 TM on, +-30ns, std TM, no NL
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111111067032230000","0163103100000010"); // MBAND
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111121067032230000","0163103100000010"); // MBAND
    cuts.AddCutPCMCalo("00052113","00200009327000008250400000","1111111067032230000","0163103100000010"); // EMC7
    cuts.AddCutPCMCalo("00052113","00200009327000008250400000","1111121067032230000","0163103100000010"); // EMC7
  } else if ( trainConfig == 261){ // EMCAL clusters 7 TeV LHC10 MB only
    cuts.AddCutPCMCalo("00000113","00200009327000008250400000","11111110b7032230000","0163103100000010"); // MBAND
    cuts.AddCutPCMCalo("00000113","00200009327000008250400000","11111210b7032230000","0163103100000010"); // MBAND

  } else if ( trainConfig == 262){ // EMCAL clusters 7 TeV LHC11 TM on, +-30ns, std TM, TB NL
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","111110106f032230000","0163103100000010"); // MBAND
    cuts.AddCutPCMCalo("00052113","00200009327000008250400000","111110106f032230000","0163103100000010"); // EMC7
  } else if ( trainConfig == 263){ // EMCAL clusters 7 TeV LHC11 TM on, +-30ns, std TM, TB NL
    cuts.AddCutPCMCalo("00010c13","00200009327000008250400000","111113106f032230000","0163103100000010"); // MBAND
    cuts.AddCutPCMCalo("00052c13","00200009327000008250400000","111113106f032230000","0163103100000010"); // EMC7
  } else if ( trainConfig == 264){ // EMCAL clusters 7 TeV LHC11 TM on, +-30ns, std TM, TB NL
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","111113206f032230000","0163103100000010"); // MBAND
    cuts.AddCutPCMCalo("00052113","00200009327000008250400000","111113206f032230000","0163103100000010"); // EMC7
  } else if ( trainConfig == 265){ // EMCAL clusters 7 TeV LHC11 TM on, +-30ns, std TM, TB NL
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","111113306f032230000","0163103100000010"); // MBAND
    cuts.AddCutPCMCalo("00052113","00200009327000008250400000","111113306f032230000","0163103100000010"); // EMC7
  } else if ( trainConfig == 266){ // EMCAL clusters 7 TeV LHC11 TM on, +-30ns, std TM, TB NL
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","111113406f032230000","0163103100000010"); // MBAND
    cuts.AddCutPCMCalo("00052113","00200009327000008250400000","111113406f032230000","0163103100000010"); // EMC7
  } else if ( trainConfig == 267){ // EMCAL clusters 7 TeV LHC11 new cuts for PCM + 31 TB NL
    cuts.AddCutPCMCalo("00010113","0dm0000922700000dge0404000","111113106f032230000","0163103100000010"); // MBAND
    cuts.AddCutPCMCalo("00052113","0dm0000922700000dge0404000","111113106f032230000","0163103100000010"); // EMC7
  } else if ( trainConfig == 268){ // EMCAL clusters 7 TeV LHC11 new cuts for PCM + 31 TB NL + SPD cut
    cuts.AddCutPCMCalo("00010c13","0dm0000922700000dge0404000","111113106f032230000","0163103100000010"); // MBAND
    cuts.AddCutPCMCalo("00052c13","0dm0000922700000dge0404000","111113106f032230000","0163103100000010"); // EMC7

  } else if ( trainConfig == 270){ // EMCAL clusters 8 TeV LHC12, TB 300 MeV
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","111113206f032230000","0163103100000010"); // 0.8
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","111113206f032230000","0163403100000010"); // 0.5
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","111113206f032230000","0163a03100000010"); // 0-0.8
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","111113206f032230000","0163b03100000010"); // -.8-0
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","111113206f032230000","0163e03100000010"); // 0-1
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","111113206f032230000","0163f03100000010"); // -1-0
  } else if ( trainConfig == 271){ // EMCAL clusters 8 TeV LHC12, TB 300 MeV
    cuts.AddCutPCMCalo("00052113","0dm00009f9730000dge0404000","111113206f032230000","0163103100000010"); // 0.8
    cuts.AddCutPCMCalo("00052113","0dm00009f9730000dge0404000","111113206f032230000","0163403100000010"); // 0.5
    cuts.AddCutPCMCalo("00052113","0dm00009f9730000dge0404000","111113206f032230000","0163a03100000010"); // 0-0.8
    cuts.AddCutPCMCalo("00052113","0dm00009f9730000dge0404000","111113206f032230000","0163b03100000010"); // -0.8-0
    cuts.AddCutPCMCalo("00052113","0dm00009f9730000dge0404000","111113206f032230000","0163e03100000010"); // 0-1
    cuts.AddCutPCMCalo("00052113","0dm00009f9730000dge0404000","111113206f032230000","0163f03100000010"); // -1-0
  } else if ( trainConfig == 272){ // EMCAL clusters 8 TeV LHC12, TB 300 MeV
    cuts.AddCutPCMCalo("00081113","0dm00009f9730000dge0404000","111113206f032230000","0163103100000010"); // 0.8
    cuts.AddCutPCMCalo("00081113","0dm00009f9730000dge0404000","111113206f032230000","0163403100000010"); // 0.5
    cuts.AddCutPCMCalo("00081113","0dm00009f9730000dge0404000","111113206f032230000","0163a03100000010"); // 0-0.8
    cuts.AddCutPCMCalo("00081113","0dm00009f9730000dge0404000","111113206f032230000","0163b03100000010"); // -0.8-0
    cuts.AddCutPCMCalo("00081113","0dm00009f9730000dge0404000","111113206f032230000","0163e03100000010"); // 0-1
    cuts.AddCutPCMCalo("00081113","0dm00009f9730000dge0404000","111113206f032230000","0163f03100000010"); // -1-0

  //multiple std dirGAMMA cuts for different studies
  } else if (trainConfig == 281){ // EMCAL clusters pp 7 TeV
    cuts.AddCutPCMCalo("00000113","00200009327000008250400000","11111110b7032230000","0163103100000010"); // std
  } else if (trainConfig == 282){ // EMCAL clusters pp 7 TeV
    cuts.AddCutPCMCalo("00000113","00200009327000008250400000","11111110b7032230000","0163103100000010"); // std
  } else if (trainConfig == 283){ // EMCAL clusters pp 7 TeV
    cuts.AddCutPCMCalo("00000113","00200009327000008250400000","11111110b7032230000","0163103100000010"); // std
  } else if (trainConfig == 284){ // EMCAL clusters pp 7 TeV
    cuts.AddCutPCMCalo("00000113","00200009327000008250400000","11111110b7032230000","0163103100000010"); // std

  //*************************************************************************************************
  // 2.76 TeV (LHC11a) PHOS setup
  //*************************************************************************************************
  } else if ( trainConfig == 300) { //PHOS clusters
    cuts.AddCutPCMCalo("00003113","00200009327000008250400000","2444400041033200000","0163103100000010");
    cuts.AddCutPCMCalo("00003113","00200009327000008250400000","2444400042033200000","0163103100000010");
    cuts.AddCutPCMCalo("00003113","00200009327000008250400000","2444400043033200000","0163103100000010");
    cuts.AddCutPCMCalo("00061113","00200009327000008250400000","2444400041033200000","0163103100000010");
    cuts.AddCutPCMCalo("00061113","00200009327000008250400000","2444400042033200000","0163103100000010");
    cuts.AddCutPCMCalo("00061113","00200009327000008250400000","2444400043033200000","0163103100000010");
  } else if ( trainConfig == 301) { //PHOS clusters without and with added signals
    cuts.AddCutPCMCalo("00003113","00200009327000008250400000","2444400042033200000","0163103100000010");
    cuts.AddCutPCMCalo("00003123","00200009327000008250400000","2444400042033200000","0163103100000010");

  //*************************************************************************************************
  // 2.76 TeV (LHC13g) PHOS setup
  //*************************************************************************************************
  } else if ( trainConfig == 310) { //PHOS clusters
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","2444400041033200000","0163103100000010");
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","2444400042033200000","0163103100000010");
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","2444400043033200000","0163103100000010");
    cuts.AddCutPCMCalo("00062113","00200009327000008250400000","2444400041033200000","0163103100000010");
    cuts.AddCutPCMCalo("00062113","00200009327000008250400000","2444400042033200000","0163103100000010");
    cuts.AddCutPCMCalo("00062113","00200009327000008250400000","2444400043033200000","0163103100000010");

  //*************************************************************************************************
  // 8 TeV PHOS setup
  //*************************************************************************************************
  } else if ( trainConfig == 330){ // PHOS clusters 8 TeV
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","2444400062023200000","0163103100000010"); // 600 MeV cluster min energy, -30,50ns timing
  } else if ( trainConfig == 331){ // PHOS clusters 8 TeV
    cuts.AddCutPCMCalo("00062113","00200009327000008250400000","2444400062023200000","0163103100000010"); // 600 MeV cluster min energy, -30,50ns timing
  } else if ( trainConfig == 332){ // With/without Added Signals
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","2444400062023200000","0163103100000010"); //
    cuts.AddCutPCMCalo("00010123","00200009327000008250400000","2444400062023200000","0163103100000010"); //
  } else if ( trainConfig == 333){ // INT7
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","2444400000013300000","0163103100000010"); // QA
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","2444400060013300000","0163103100000010"); // QA, -30,50ns timing
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","2444400063013300000","0163103100000010"); // QA, -30,50ns timing, TM on with default EMC params
  } else if ( trainConfig == 334){ // PHI7
    cuts.AddCutPCMCalo("00062113","00200009327000008250400000","2444400000013300000","0163103100000010"); // QA
    cuts.AddCutPCMCalo("00062113","00200009327000008250400000","2444400060013300000","0163103100000010"); // QA, -30,50ns timing
    cuts.AddCutPCMCalo("00062113","00200009327000008250400000","2444400063013300000","0163103100000010"); // QA, -30,50ns timing, TM on with default EMC params
  // PHOS new default with timing effi
  } else if ( trainConfig == 340){
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0400000","24444000ga012200000","0h63103100000010"); // No NL
  } else if ( trainConfig == 341){
    cuts.AddCutPCMCalo("00062113","00200009f9730000dge0400000","24444000ga012200000","0h63103100000010"); // No NL
  } else if ( trainConfig == 342){
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0400000","24444590ga012200000","0h63103100000010"); // 59 NL
  } else if ( trainConfig == 343){
    cuts.AddCutPCMCalo("00062113","00200009f9730000dge0400000","24444590ga012200000","0h63103100000010"); // 59 NL
  } else if ( trainConfig == 344){
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0400000","24444690ga012200000","0h63103100000010"); // 69 NL
  } else if ( trainConfig == 345){
    cuts.AddCutPCMCalo("00062113","00200009f9730000dge0400000","24444690ga012200000","0h63103100000010"); // 69 NL
  //*************************************************************************************************
  // 7 TeV PHOS setup
  //*************************************************************************************************
  } else if ( trainConfig == 350){
    cuts.AddCutPCMCalo("00000113","00200009327000008250400000","2444400000012300000","0163103100000010"); // QA
  } else if ( trainConfig == 351){
    cuts.AddCutPCMCalo("00000113","00200009327000008250400000","2444400000013300000","0163103100000010"); // QA
  } else if ( trainConfig == 352){
    cuts.AddCutPCMCalo("00000113","00200009327000008250400000","2444400040013300000","0163103100000010"); // 100ns timing cut, no track matching
    cuts.AddCutPCMCalo("00000113","00200009327000008250400000","2444400043013300000","0163103100000010"); // 100ns timing cut
    cuts.AddCutPCMCalo("00000113","00200009327000008250400000","2444400043013350000","0163103100000010"); // 100ns timing cut, M02<0.3
    cuts.AddCutPCMCalo("00000113","00200009327000008250400000","2444400043013330000","0163103100000010"); // 100ns timing cut, M02<0.5
    cuts.AddCutPCMCalo("00000113","00200009327000008250400000","2444400043013320000","0163103100000010"); // 100ns timing cut, M02<0.7
  } else if ( trainConfig == 353){
    cuts.AddCutPCMCalo("00000113","00200009327000008250400000","2444411000013300000","0163103100000010"); // own constant ConvCalo NonLin
    cuts.AddCutPCMCalo("00000113","00200009327000008250400000","2444421000013300000","0163103100000010"); // ext. PHOS NonLin * correction

  //*************************************************************************************************
  // 7 TeV EMC for omega analysis variations
  //*************************************************************************************************
  } else if ( trainConfig == 360){ // std
    cuts.AddCutPCMCalo("00000113","00200009227000008250400000","1111111047032230000","0163103100000010"); // std
  } else if ( trainConfig == 361){ // pileup
    cuts.AddCutPCMCalo("00000113","00200009227000008250400000","1111111047032230000","0163103100000010");
  } else if ( trainConfig == 362){ // singlept
    cuts.AddCutPCMCalo("00000113","00200019227000008250400000","1111111047032230000","0163103100000010"); // 0.100 GeV
    cuts.AddCutPCMCalo("00000113","00200049227000008250400000","1111111047032230000","0163103100000010"); // 0.075 GeV
    cuts.AddCutPCMCalo("00000113","00200069227000008250400000","1111111047032230000","0163103100000010"); // 0.04 GeV
    cuts.AddCutPCMCalo("00000113","00200059227000008250400000","1111111047032230000","0163103100000010"); // 0.125 GeV
  } else if ( trainConfig == 363){ // clstpc
    cuts.AddCutPCMCalo("00000113","00200008227000008250400000","1111111047032230000","0163103100000010"); // fMinClsTPCToF= 0.35;fUseCorrectedTPCClsInfo=1;
    cuts.AddCutPCMCalo("00000113","00200006227000008250400000","1111111047032230000","0163103100000010"); // fMinClsTPCToF= 0.70;fUseCorrectedTPCClsInfo=1;
    cuts.AddCutPCMCalo("00000113","00200001227000008250400000","1111111047032230000","0163103100000010"); // fMinClsTPCToF= 60;fUseCorrectedTPCClsInfo=0;
    cuts.AddCutPCMCalo("00000113","00200002227000008250400000","1111111047032230000","0163103100000010"); // fMinClsTPCToF= 80;fUseCorrectedTPCClsInfo=0;
    cuts.AddCutPCMCalo("00000113","00200003227000008250400000","1111111047032230000","0163103100000010"); // fMinClsTPCToF= 100;fUseCorrectedTPCClsInfo=0;
  } else if ( trainConfig == 364){ // TPCdEdxCutElectron
    cuts.AddCutPCMCalo("00000113","00200009327000008250400000","1111111047032230000","0163103100000010"); // -4,5
    cuts.AddCutPCMCalo("00000113","00200009627000008250400000","1111111047032230000","0163103100000010"); // -2.5,4
    cuts.AddCutPCMCalo("00000113","00200009427000008250400000","1111111047032230000","0163103100000010"); // -6,7
    cuts.AddCutPCMCalo("00000113","00200009527000008250400000","1111111047032230000","0163103100000010"); // -4,4
    cuts.AddCutPCMCalo("00000113","00200009627000008250400000","1111111047032230000","0163103100000010"); // -2.5,4
  } else if ( trainConfig == 365){ // TPCdEdxCutPion
    cuts.AddCutPCMCalo("00000113","00200009217000008250400000","1111111047032230000","0163103100000010"); // fPIDnSigmaAbovePionLine=0; fPIDnSigmaAbovePionLineHighPt=-10;
    cuts.AddCutPCMCalo("00000113","00200009237000008250400000","1111111047032230000","0163103100000010"); // fPIDnSigmaAbovePionLine=2.5; fPIDnSigmaAbovePionLineHighPt=-10;
    cuts.AddCutPCMCalo("00000113","00200009247000008250400000","1111111047032230000","0163103100000010"); // fPIDnSigmaAbovePionLine=0.5; fPIDnSigmaAbovePionLineHighPt=-10;
    cuts.AddCutPCMCalo("00000113","00200009257000008250400000","1111111047032230000","0163103100000010"); // fPIDnSigmaAbovePionLine=2; fPIDnSigmaAbovePionLineHighPt=-10;
  } else if ( trainConfig == 366){ // QtMaxCut
    cuts.AddCutPCMCalo("00000113","00200009227000003250400000","1111111047032230000","0163103100000010");
    cuts.AddCutPCMCalo("00000113","00200009227000009250400000","1111111047032230000","0163103100000010");
    cuts.AddCutPCMCalo("00000113","00200009227000002250400000","1111111047032230000","0163103100000010");
    cuts.AddCutPCMCalo("00000113","00200009227000006250400000","1111111047032230000","0163103100000010");
  } else if ( trainConfig == 367){ // Chi2GammaCut
    cuts.AddCutPCMCalo("00000113","00200009227000008150400000","1111111047032230000","0163103100000010");
    cuts.AddCutPCMCalo("00000113","00200009227000008850400000","1111111047032230000","0163103100000010");
    cuts.AddCutPCMCalo("00000113","00200009227000008a50400000","1111111047032230000","0163103100000010");
    cuts.AddCutPCMCalo("00000113","00200009227000008950400000","1111111047032230000","0163103100000010");
  } else if ( trainConfig == 368){ // PsiPair
    cuts.AddCutPCMCalo("00000113","00200009227000008260400000","1111111047032230000","0163103100000010");
    cuts.AddCutPCMCalo("00000113","00200009227000008280400000","1111111047032230000","0163103100000010");
    cuts.AddCutPCMCalo("00000113","00200009227000008210400000","1111111047032230000","0163103100000010");
    cuts.AddCutPCMCalo("00000113","00200009227000008220400000","1111111047032230000","0163103100000010");

  } else if ( trainConfig == 369){ // NonLin
    cuts.AddCutPCMCalo("00000113","00200009227000008250400000","1111112047032230000","0163103100000010");
    cuts.AddCutPCMCalo("00000113","00200009227000008250400000","1111113047032230000","0163103100000010");
    cuts.AddCutPCMCalo("00000113","00200009227000008250400000","1111121047032230000","0163103100000010");
    cuts.AddCutPCMCalo("00000113","00200009227000008250400000","1111122047032230000","0163103100000010");
    cuts.AddCutPCMCalo("00000113","00200009227000008250400000","1111123047032230000","0163103100000010");
  } else if ( trainConfig == 370){ // Timing diff
    cuts.AddCutPCMCalo("00000113","00200009227000008250400000","11111110b7032230000","0163103100000010");
    cuts.AddCutPCMCalo("00000113","00200009227000008250400000","11111110c7032230000","0163103100000010");
    cuts.AddCutPCMCalo("00000113","00200009227000008250400000","11111110d7032230000","0163103100000010");
    cuts.AddCutPCMCalo("00000113","00200009227000008250400000","11111110e7032230000","0163103100000010");
    cuts.AddCutPCMCalo("00000113","00200009227000008250400000","11111110f7032230000","0163103100000010");
  } else if ( trainConfig == 371){ // Track Matching
    cuts.AddCutPCMCalo("00000113","00200009227000008250400000","1111111046032230000","0163103100000010");
    cuts.AddCutPCMCalo("00000113","00200009227000008250400000","1111111048032230000","0163103100000010");
    cuts.AddCutPCMCalo("00000113","00200009227000008250400000","1111111049032230000","0163103100000010");
    cuts.AddCutPCMCalo("00000113","00200009227000008250400000","111111104a032230000","0163103100000010");
    cuts.AddCutPCMCalo("00000113","00200009227000008250400000","111111104b032230000","0163103100000010");
    cuts.AddCutPCMCalo("00000113","00200009227000008250400000","1111111043032230000","0163103100000010");
  } else if ( trainConfig == 372){ // MinEnergy (of cluster)
    cuts.AddCutPCMCalo("00000113","00200009227000008250400000","1111111047022230000","0163103100000010");
    cuts.AddCutPCMCalo("00000113","00200009227000008250400000","1111111047042230000","0163103100000010");
    cuts.AddCutPCMCalo("00000113","00200009227000008250400000","1111111047052230000","0163103100000010");
  } else if ( trainConfig == 373){ // MinNCells
    cuts.AddCutPCMCalo("00000113","00200009227000008250400000","1111111047032130000","0163103100000010");
    cuts.AddCutPCMCalo("00000113","00200009227000008250400000","1111111047032330000","0163103100000010");
  } else if ( trainConfig == 374){ // MinMaxM02
    cuts.AddCutPCMCalo("00000113","00200009227000008250400000","1111111047032330000","0163103100000010");
    cuts.AddCutPCMCalo("00000113","00200009227000008250400000","1111111047032130000","0163103100000010");
    cuts.AddCutPCMCalo("00000113","00200009227000008250400000","1111111047032240000","0163103100000010");
    cuts.AddCutPCMCalo("00000113","00200009227000008250400000","1111111047032220000","0163103100000010");

 //*************************************************************************************************
  // 7 TeV PHOS for omega analysis variations
  //*************************************************************************************************
  } else if ( trainConfig == 380){ // std
    cuts.AddCutPCMCalo("00000113","00200009227000008250400000","2444411044013300000","0163103100000010"); // std
  } else if ( trainConfig == 381){ // pileup
    cuts.AddCutPCMCalo("00000113","00200009227000008250400000","2444411044013300000","0163103100000010");
  } else if ( trainConfig == 382){ // singlept
    cuts.AddCutPCMCalo("00000113","00200019227000008250400000","2444411044013300000","0163103100000010"); // 0.100 GeV
    cuts.AddCutPCMCalo("00000113","00200049227000008250400000","2444411044013300000","0163103100000010"); // 0.075 GeV
    cuts.AddCutPCMCalo("00000113","00200069227000008250400000","2444411044013300000","0163103100000010"); // 0.04 GeV
    cuts.AddCutPCMCalo("00000113","00200059227000008250400000","2444411044013300000","0163103100000010"); // 0.125 GeV
  } else if ( trainConfig == 383){ // clstpc
    cuts.AddCutPCMCalo("00000113","00200008227000008250400000","2444411044013300000","0163103100000010"); // fMinClsTPCToF= 0.35;fUseCorrectedTPCClsInfo=1;
    cuts.AddCutPCMCalo("00000113","00200006227000008250400000","2444411044013300000","0163103100000010"); // fMinClsTPCToF= 0.70;fUseCorrectedTPCClsInfo=1;
    cuts.AddCutPCMCalo("00000113","00200001227000008250400000","2444411044013300000","0163103100000010"); // fMinClsTPCToF= 60;fUseCorrectedTPCClsInfo=0;
    cuts.AddCutPCMCalo("00000113","00200002227000008250400000","2444411044013300000","0163103100000010"); // fMinClsTPCToF= 80;fUseCorrectedTPCClsInfo=0;
    cuts.AddCutPCMCalo("00000113","00200003227000008250400000","2444411044013300000","0163103100000010"); // fMinClsTPCToF= 100;fUseCorrectedTPCClsInfo=0;
  } else if ( trainConfig == 384){ // TPCdEdxCutElectron
    cuts.AddCutPCMCalo("00000113","00200009327000008250400000","2444411044013300000","0163103100000010"); // -4,5
    cuts.AddCutPCMCalo("00000113","00200009627000008250400000","2444411044013300000","0163103100000010"); // -2.5,4
    cuts.AddCutPCMCalo("00000113","00200009427000008250400000","2444411044013300000","0163103100000010"); // -6,7
    cuts.AddCutPCMCalo("00000113","00200009527000008250400000","2444411044013300000","0163103100000010"); // -4,4
    cuts.AddCutPCMCalo("00000113","00200009627000008250400000","2444411044013300000","0163103100000010"); // -2.5,4
  } else if ( trainConfig == 385){ // TPCdEdxCutPion
    cuts.AddCutPCMCalo("00000113","00200009217000008250400000","2444411044013300000","0163103100000010"); // fPIDnSigmaAbovePionLine=0; fPIDnSigmaAbovePionLineHighPt=-10;
    cuts.AddCutPCMCalo("00000113","00200009237000008250400000","2444411044013300000","0163103100000010"); // fPIDnSigmaAbovePionLine=2.5; fPIDnSigmaAbovePionLineHighPt=-10;
    cuts.AddCutPCMCalo("00000113","00200009247000008250400000","2444411044013300000","0163103100000010"); // fPIDnSigmaAbovePionLine=0.5; fPIDnSigmaAbovePionLineHighPt=-10;
    cuts.AddCutPCMCalo("00000113","00200009257000008250400000","2444411044013300000","0163103100000010"); // fPIDnSigmaAbovePionLine=2; fPIDnSigmaAbovePionLineHighPt=-10;
  } else if ( trainConfig == 386){ // QtMaxCut
    cuts.AddCutPCMCalo("00000113","00200009227000003250400000","2444411044013300000","0163103100000010");
    cuts.AddCutPCMCalo("00000113","00200009227000009250400000","2444411044013300000","0163103100000010");
    cuts.AddCutPCMCalo("00000113","00200009227000002250400000","2444411044013300000","0163103100000010");
    cuts.AddCutPCMCalo("00000113","00200009227000006250400000","2444411044013300000","0163103100000010");
  } else if ( trainConfig == 387){ // Chi2GammaCut
    cuts.AddCutPCMCalo("00000113","00200009227000008150400000","2444411044013300000","0163103100000010");
    cuts.AddCutPCMCalo("00000113","00200009227000008850400000","2444411044013300000","0163103100000010");
    cuts.AddCutPCMCalo("00000113","00200009227000008a50400000","2444411044013300000","0163103100000010");
    cuts.AddCutPCMCalo("00000113","00200009227000008950400000","2444411044013300000","0163103100000010");
  } else if ( trainConfig == 388){ // PsiPair
    cuts.AddCutPCMCalo("00000113","00200009227000008260400000","2444411044013300000","0163103100000010");
    cuts.AddCutPCMCalo("00000113","00200009227000008280400000","2444411044013300000","0163103100000010");
    cuts.AddCutPCMCalo("00000113","00200009227000008210400000","2444411044013300000","0163103100000010");
    cuts.AddCutPCMCalo("00000113","00200009227000008220400000","2444411044013300000","0163103100000010");

  } else if ( trainConfig == 389){ // NonLin
    cuts.AddCutPCMCalo("00000113","00200009227000008250400000","2444421044012300000","0163103100000010");
    cuts.AddCutPCMCalo("00000113","00200009227000008250400000","2444412044012300000","0163103100000010");
  } else if ( trainConfig == 390){ // Timing diff
    cuts.AddCutPCMCalo("00000113","00200009227000008250400000","24444110b4013300000","0163103100000010");
    cuts.AddCutPCMCalo("00000113","00200009227000008250400000","24444110c4013300000","0163103100000010");
    cuts.AddCutPCMCalo("00000113","00200009227000008250400000","24444110d4013300000","0163103100000010");
    cuts.AddCutPCMCalo("00000113","00200009227000008250400000","24444110e4013300000","0163103100000010");
    cuts.AddCutPCMCalo("00000113","00200009227000008250400000","24444110f4013300000","0163103100000010");
  } else if ( trainConfig == 391){ // Track Matching
    cuts.AddCutPCMCalo("00000113","00200009227000008250400000","2444411041013300000","0163103100000010");
    cuts.AddCutPCMCalo("00000113","00200009227000008250400000","2444411043013300000","0163103100000010");
    cuts.AddCutPCMCalo("00000113","00200009227000008250400000","2444411045013300000","0163103100000010");
    cuts.AddCutPCMCalo("00000113","00200009227000008250400000","2444411046013300000","0163103100000010");
    cuts.AddCutPCMCalo("00000113","00200009227000008250400000","2444411047013300000","0163103100000010");
    cuts.AddCutPCMCalo("00000113","00200009227000008250400000","2444411048013300000","0163103100000010");
    cuts.AddCutPCMCalo("00000113","00200009227000008250400000","2444411049013300000","0163103100000010");
  } else if ( trainConfig == 392){ // MinEnergy (of cluster)
    cuts.AddCutPCMCalo("00000113","00200009227000008250400000","2444411044023300000","0163103100000010");
    cuts.AddCutPCMCalo("00000113","00200009227000008250400000","2444411044033300000","0163103100000010");
    cuts.AddCutPCMCalo("00000113","00200009227000008250400000","2444411044043300000","0163103100000010");
    cuts.AddCutPCMCalo("00000113","00200009227000008250400000","2444411044083300000","0163103100000010");
  } else if ( trainConfig == 393){ // MinNCells
    cuts.AddCutPCMCalo("00000113","00200009227000008250400000","2444411044011300000","0163103100000010");
    cuts.AddCutPCMCalo("00000113","00200009227000008250400000","2444411044012300000","0163103100000010");
    cuts.AddCutPCMCalo("00000113","00200009227000008250400000","2444411044014300000","0163103100000010");
  } else if ( trainConfig == 394){ // MinMaxM02
    cuts.AddCutPCMCalo("00000113","00200009227000008250400000","2444411044013200000","0163103100000010");
    cuts.AddCutPCMCalo("00000113","00200009227000008250400000","2444411044013100000","0163103100000010");

  //*************************************************************************************************
  // 5 TeV EDC setup
  //*************************************************************************************************
  } else if ( trainConfig == 400){ // with and without smearing
    cuts.AddCutPCMCalo("80010103","0dm00009f9730000dge0404000","411790106fg32230000","0163103100b00010");
    cuts.AddCutPCMCalo("80010103","0dm00009f9730000dge0404000","411790106fg32230000","0163103100000010");
  } else if ( trainConfig == 401){ // with and without smearing
    cuts.AddCutPCMCalo("80010103","0dm00009f9730000dge0404000","411790206fg32230000","0163103100b00010");
    cuts.AddCutPCMCalo("80010103","0dm00009f9730000dge0404000","411790206fg32230000","0163103100000010");

    //5 TeV JJ MC configs without trackmatching
  } else if ( trainConfig == 430){ // EMCAL clusters no NonLin
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111100060032230000","0163103100000010"); // -30ns, 35ns timing cut, no NL INT7
    cuts.AddCutPCMCalo("00052013","00200009327000008250400000","1111100060032230000","0163103100000010"); // -30ns, 35ns timing cut, no NL EMC7
    cuts.AddCutPCMCalo("00085013","00200009327000008250400000","1111100060032230000","0163103100000010"); // -30ns, 35ns timing cut, no NL EG2
    cuts.AddCutPCMCalo("00083013","00200009327000008250400000","1111100060032230000","0163103100000010"); // -30ns, 35ns timing cut, no NL EG1
  } else if ( trainConfig == 431){ // EMCAL clusters - NonLin INT7
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111111060032230000","0163103100000010"); // -30ns, 35ns timing cut, NL kSDM PCMEMC
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111112060032230000","0163103100000010"); // -30ns, 35ns timing cut, NL kSDM EMC
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111121060032230000","0163103100000010"); // -30ns, 35ns timing cut, NL DExt PCMEMC
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111122060032230000","0163103100000010"); // -30ns, 35ns timing cut, NL DExt EMC
    // 5 TeV PCM-EMC Systematics
  } else if ( trainConfig == 450){ // PCM based systematics
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111127067032230000","0163103100000010"); // std
    cuts.AddCutPCMCalo("00010013","00200009327000008250400000","1111127067032230000","0163103100000010"); // std - no pileup cut
    cuts.AddCutPCMCalo("00010013","00202209327000008250400000","1111127067032230000","0163103100000010"); // restrict acceptance to EMCAL loose
    cuts.AddCutPCMCalo("00010013","00204409327000008250400000","1111127067032230000","0163103100000010"); // restrict acceptance to EMCAL tight
  } else if ( trainConfig == 451){ // PCM based systematics
    cuts.AddCutPCMCalo("00010113","00200009327400008250400000","1111127067032230000","0163103100000010"); // dEdx pi: 1: 0.4-3, -10: 3. ->
    cuts.AddCutPCMCalo("00010113","00200009367400008250400000","1111127067032230000","0163103100000010"); // dEdx pi: 2: 0.4-3, 0.5: 3. ->
    cuts.AddCutPCMCalo("00010113","00200009227000008250400000","1111127067032230000","0163103100000010"); // dEdx e -3, 5
    cuts.AddCutPCMCalo("00010113","00200009127000008250400000","1111127067032230000","0163103100000010"); // dEdx e -5, 5
  } else if ( trainConfig == 452){ // PCM based systematics
    cuts.AddCutPCMCalo("00010113","00200049327000008250400000","1111127067032230000","0163103100000010"); // single pt > 0.075
    cuts.AddCutPCMCalo("00010113","00200019327000008250400000","1111127067032230000","0163103100000010"); // single pt > 0.1
    cuts.AddCutPCMCalo("00010113","00200009327000009250400000","1111127067032230000","0163103100000010"); // qt 2D 0.03
    cuts.AddCutPCMCalo("00010113","00200009327000003250400000","1111127067032230000","0163103100000010"); // qt 1D 0.05
    cuts.AddCutPCMCalo("00010113","00200009327000002250400000","1111127067032230000","0163103100000010"); // qt 1D 0.07
  } else if ( trainConfig == 453){ // PCM based systematics
    cuts.AddCutPCMCalo("00010113","00200009327000008850400000","1111127067032230000","0163103100000010"); // 2D psi pair chi2 var
    cuts.AddCutPCMCalo("00010113","00200009327000008260400000","1111127067032230000","0163103100000010"); // 2D psi pair chi2 var
    cuts.AddCutPCMCalo("00010113","00200009327000008860400000","1111127067032230000","0163103100000010"); // 2D psi pair chi2 var
    cuts.AddCutPCMCalo("00010113","00200009327000008280400000","1111127067032230000","0163103100000010"); // 2D psi pair chi2 var
    cuts.AddCutPCMCalo("00010113","00200009327000008880400000","1111127067032230000","0163103100000010"); // 2D psi pair chi2 var
  } else if ( trainConfig == 454){ // PCM based systematics
    cuts.AddCutPCMCalo("00010113","00200006327000008250400000","1111127067032230000","0163103100000010"); // min TPC cl > 0.7
    cuts.AddCutPCMCalo("00010113","00200008327000008250400000","1111127067032230000","0163103100000010"); // min TPC cl > 0.35
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111127067032230000","0163106100000010"); // alpha < 0.8
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111127067032230000","0163105100000010"); // alpha < 0.75

  } else if ( trainConfig == 460){ // EMC min E variations
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111127067022230000","0163103100000010"); //0.6 GeV/c
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111127067042230000","0163103100000010"); //0.8 GeV/c
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111127067052230000","0163103100000010"); //0.9 GeV/c
  } else if ( trainConfig == 461){ //EMCAL minNCells variation
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111127067031230000","0163103100000010"); //n cells >= 1
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111127067033230000","0163103100000010"); //n cells >= 3
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111127067032200000","0163103100000010"); //no max M02 cut
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111127067032250000","0163103100000010"); //M02 < 0.3
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","11111270670322k0000","0163103100000010"); //M02, pT-dep
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1112127067032230000","0163103100000010"); //only modules with TRD infront
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111327067032230000","0163103100000010"); //no modules with TRD infront
  } else if ( trainConfig == 462){ // EMCAL track matching variations
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111127066032230000","0163103100000010"); // track matching variations
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111127068032230000","0163103100000010"); //
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111127069032230000","0163103100000010"); //
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111127063032230000","0163103100000010"); //
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111127060032230000","0163103100000010"); //
  } else if ( trainConfig == 463){ // EMCAL clusters, timing variation
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111127057032230000","0163103100000010"); // time 50ns
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111127047032230000","0163103100000010"); // time 100ns
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111127037032230000","0163103100000010"); // time 200ns
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111127027032230000","0163103100000010"); // time 500ns
  } else if ( trainConfig == 464){ // EMCAL clusters, timing variation
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111127077032230000","0163103100000010"); // time
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111127087032230000","0163103100000010"); // time
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111127097032230000","0163103100000010"); // time

  } else if ( trainConfig == 470){ // PCM based systematics
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111111067032230000","0163103100000010"); // std
    cuts.AddCutPCMCalo("00010013","00200009327000008250400000","1111111067032230000","0163103100000010"); // std - no pileup cut
    cuts.AddCutPCMCalo("00010013","00202209327000008250400000","1111111067032230000","0163103100000010"); // restrict acceptance to EMCAL loose
    cuts.AddCutPCMCalo("00010013","00204409327000008250400000","1111111067032230000","0163103100000010"); // restrict acceptance to EMCAL tight
  } else if ( trainConfig == 471){ // PCM based systematics
    cuts.AddCutPCMCalo("00010113","00200009327400008250400000","1111111067032230000","0163103100000010"); // dEdx pi: 1: 0.4-3, -10: 3. ->
    cuts.AddCutPCMCalo("00010113","00200009367400008250400000","1111111067032230000","0163103100000010"); // dEdx pi: 2: 0.4-3, 0.5: 3. ->
    cuts.AddCutPCMCalo("00010113","00200009227000008250400000","1111111067032230000","0163103100000010"); // dEdx e -3, 5
    cuts.AddCutPCMCalo("00010113","00200009127000008250400000","1111111067032230000","0163103100000010"); // dEdx e -5, 5
  } else if ( trainConfig == 472){ // PCM based systematics
    cuts.AddCutPCMCalo("00010113","00200049327000008250400000","1111111067032230000","0163103100000010"); // single pt > 0.075
    cuts.AddCutPCMCalo("00010113","00200019327000008250400000","1111111067032230000","0163103100000010"); // single pt > 0.1
    cuts.AddCutPCMCalo("00010113","00200009327000009250400000","1111111067032230000","0163103100000010"); // qt 2D 0.03
    cuts.AddCutPCMCalo("00010113","00200009327000003250400000","1111111067032230000","0163103100000010"); // qt 1D 0.05
    cuts.AddCutPCMCalo("00010113","00200009327000002250400000","1111111067032230000","0163103100000010"); // qt 1D 0.07
  } else if ( trainConfig == 473){ // PCM based systematics
    cuts.AddCutPCMCalo("00010113","00200009327000008850400000","1111111067032230000","0163103100000010"); // 2D psi pair chi2 var
    cuts.AddCutPCMCalo("00010113","00200009327000008260400000","1111111067032230000","0163103100000010"); // 2D psi pair chi2 var
    cuts.AddCutPCMCalo("00010113","00200009327000008860400000","1111111067032230000","0163103100000010"); // 2D psi pair chi2 var
    cuts.AddCutPCMCalo("00010113","00200009327000008280400000","1111111067032230000","0163103100000010"); // 2D psi pair chi2 var
    cuts.AddCutPCMCalo("00010113","00200009327000008880400000","1111111067032230000","0163103100000010"); // 2D psi pair chi2 var
  } else if ( trainConfig == 474){ // PCM based systematics
    cuts.AddCutPCMCalo("00010113","00200006327000008250400000","1111111067032230000","0163103100000010"); // min TPC cl > 0.7
    cuts.AddCutPCMCalo("00010113","00200008327000008250400000","1111111067032230000","0163103100000010"); // min TPC cl > 0.35
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111111067032230000","0163106100000010"); // alpha < 0.8
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111111067032230000","0163105100000010"); // alpha < 0.75

  } else if ( trainConfig == 480){ // EMC min E variations
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111111067022230000","0163103100000010"); //0.6 GeV/c
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111111067042230000","0163103100000010"); //0.8 GeV/c
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111111067052230000","0163103100000010"); //0.9 GeV/c
  } else if ( trainConfig == 481){ //EMCAL minNCells variation
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111111067031230000","0163103100000010"); //n cells >= 1
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111111067033230000","0163103100000010"); //n cells >= 3
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111111067032200000","0163103100000010"); //no max M02 cut
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111111067032250000","0163103100000010"); //M02 < 0.3
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","11111110670322k0000","0163103100000010"); //M02, pT-dep
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1112111067032230000","0163103100000010"); //only modules with TRD infront
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111311067032230000","0163103100000010"); //no modules with TRD infront
  } else if ( trainConfig == 482){ // EMCAL track matching variations
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111111066032230000","0163103100000010"); // track matching variations
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111111068032230000","0163103100000010"); //
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111111069032230000","0163103100000010"); //
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111111063032230000","0163103100000010"); //
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111111060032230000","0163103100000010"); //
  } else if ( trainConfig == 483){ // EMCAL clusters, timing variation
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111111057032230000","0163103100000010"); // time 50ns
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111111047032230000","0163103100000010"); // time 100ns
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111111037032230000","0163103100000010"); // time 200ns
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111111027032230000","0163103100000010"); // time 500ns
  } else if ( trainConfig == 484){ // EMCAL clusters, timing variation
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111111077032230000","0163103100000010"); // time
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111111087032230000","0163103100000010"); // time
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111111097032230000","0163103100000010"); // time

  } else if ( trainConfig == 490){ // INT7 DCAL standard cut but with E/p TM veto
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111111067032230000","0163103100000010"); // std INT7
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","111111106c032230000","0163103100000010"); // fEOverPMax = 9e9
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","111111106d032230000","0163103100000010"); // fEOverPMax = 3.0
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","111111106e032230000","0163103100000010"); // fEOverPMax = 2.0
  } else if ( trainConfig == 491){ // INT7 DCAL standard cut but with E/p TM veto
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","111111106f032230000","0163103100000010"); // fEOverPMax = 1.75
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","111111106g032230000","0163103100000010"); // fEOverPMax = 1.5
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","111111106h032230000","0163103100000010"); // fEOverPMax = 1.25
  } else if ( trainConfig == 492){ // DMC7 DCAL standard cut but with E/p TM veto
    cuts.AddCutPCMCalo("00052113","00200009327000008250400000","1111111067032230000","0163103100000010"); // std DMC7
    cuts.AddCutPCMCalo("00052113","00200009327000008250400000","111111106c032230000","0163103100000010"); // fEOverPMax = 9e9
    cuts.AddCutPCMCalo("00052113","00200009327000008250400000","111111106d032230000","0163103100000010"); // fEOverPMax = 3.0
    cuts.AddCutPCMCalo("00052113","00200009327000008250400000","111111106e032230000","0163103100000010"); // fEOverPMax = 2.0
  } else if ( trainConfig == 493){ // DMC7 DCAL standard cut but with E/p TM veto
    cuts.AddCutPCMCalo("00052113","00200009327000008250400000","111111106f032230000","0163103100000010"); // fEOverPMax = 1.75
    cuts.AddCutPCMCalo("00052113","00200009327000008250400000","111111106g032230000","0163103100000010"); // fEOverPMax = 1.5
    cuts.AddCutPCMCalo("00052113","00200009327000008250400000","111111106h032230000","0163103100000010"); // fEOverPMax = 1.25
  } else if ( trainConfig == 494){ // DG1 DCAL standard cut but with E/p TM veto
    cuts.AddCutPCMCalo("00083113","00200009327000008250400000","1111111067032230000","0163103100000010"); // std DG1
    cuts.AddCutPCMCalo("00083113","00200009327000008250400000","111111106c032230000","0163103100000010"); // fEOverPMax = 9e9
    cuts.AddCutPCMCalo("00083113","00200009327000008250400000","111111106d032230000","0163103100000010"); // fEOverPMax = 3.0
    cuts.AddCutPCMCalo("00083113","00200009327000008250400000","111111106e032230000","0163103100000010"); // fEOverPMax = 2.0
  } else if ( trainConfig == 495){ // DG1 DCAL standard cut but with E/p TM veto
    cuts.AddCutPCMCalo("00083113","00200009327000008250400000","111111106f032230000","0163103100000010"); // fEOverPMax = 1.75
    cuts.AddCutPCMCalo("00083113","00200009327000008250400000","111111106g032230000","0163103100000010"); // fEOverPMax = 1.5
    cuts.AddCutPCMCalo("00083113","00200009327000008250400000","111111106h032230000","0163103100000010"); // fEOverPMax = 1.25
  } else if ( trainConfig == 496){ // DG2 DCAL standard cut but with E/p TM veto
    cuts.AddCutPCMCalo("00085113","00200009327000008250400000","1111111067032230000","0163103100000010"); // std DG2
    cuts.AddCutPCMCalo("00085113","00200009327000008250400000","111111106c032230000","0163103100000010"); // fEOverPMax = 9e9
    cuts.AddCutPCMCalo("00085113","00200009327000008250400000","111111106d032230000","0163103100000010"); // fEOverPMax = 3.0
    cuts.AddCutPCMCalo("00085113","00200009327000008250400000","111111106e032230000","0163103100000010"); // fEOverPMax = 2.0
  } else if ( trainConfig == 497){ // DG2 DCAL standard cut but with E/p TM veto
    cuts.AddCutPCMCalo("00085113","00200009327000008250400000","111111106f032230000","0163103100000010"); // fEOverPMax = 1.75
    cuts.AddCutPCMCalo("00085113","00200009327000008250400000","111111106g032230000","0163103100000010"); // fEOverPMax = 1.5
    cuts.AddCutPCMCalo("00085113","00200009327000008250400000","111111106h032230000","0163103100000010"); // fEOverPMax = 1.25

  } else if ( trainConfig == 498){ // EMCAL standard cut but standard readout - NO TM
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111111060032230000","0163103100000010"); // std INT7
    cuts.AddCutPCMCalo("00052113","00200009327000008250400000","1111111060032230000","0163103100000010"); // std EMC7
    cuts.AddCutPCMCalo("00085113","00200009327000008250400000","1111111060032230000","0163103100000010"); // std EG2
    cuts.AddCutPCMCalo("00083113","00200009327000008250400000","1111111060032230000","0163103100000010"); // std EG1
  } else if ( trainConfig == 499){ // EMCAL standard cut but CALO+CALOFAST - NO TM
    cuts.AddCutPCMCalo("000a0113","00200009327000008250400000","1111111060032230000","0163103100000010"); // std INT7
    cuts.AddCutPCMCalo("000a1113","00200009327000008250400000","1111111060032230000","0163103100000010"); // std EMC7
    cuts.AddCutPCMCalo("000a2113","00200009327000008250400000","1111111060032230000","0163103100000010"); // std EG2
    cuts.AddCutPCMCalo("000a3113","00200009327000008250400000","1111111060032230000","0163103100000010"); // std EG1
  //*************************************************************************************************
  // 13 TeV EMC setup
  //*************************************************************************************************
  } else if ( trainConfig == 500){ // EMCAL clusters
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111100017032230000","0163103100000010"); // 1000ns timing cut, no NL INT7
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111100067032230000","0163103100000010"); // -30ns, 35ns timing cut, no NL INT7
  } else if ( trainConfig == 501){ // EMCAL clusters
    cuts.AddCutPCMCalo("00052113","00200009327000008250400000","1111100017032230000","0163103100000010"); // 1000ns timing cut, no NL EMC7
    cuts.AddCutPCMCalo("00052113","00200009327000008250400000","1111100067032230000","0163103100000010"); // -30ns, 35ns timing cut, no NL EMC7
  } else if ( trainConfig == 502){ // EMCAL clusters
    cuts.AddCutPCMCalo("00085113","00200009327000008250400000","1111100067032230000","0163103100000010"); // -30ns, 35ns timing cut, no NL EG2
    cuts.AddCutPCMCalo("00083113","00200009327000008250400000","1111100067032230000","0163103100000010"); // -30ns, 35ns timing cut, no NL EG1
  } else if ( trainConfig == 503){ // EMCAL clusters NL vars
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111111067032230000","0163103100000010"); // NL
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111112067032230000","0163103100000010"); // NL
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111121067032230000","0163103100000010"); // NL
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111122067032230000","0163103100000010"); // NL
  } else if ( trainConfig == 504){ // EMCAL clusters
    cuts.AddCutPCMCalo("00052113","00200009327000008250400000","1111111067032230000","0163103100000010"); // -30ns, 35ns timing cut, NL1 EMC7
    cuts.AddCutPCMCalo("00052113","00200009327000008250400000","1111112067032230000","0163103100000010"); // -30ns, 35ns timing cut, NL2 EMC7
    cuts.AddCutPCMCalo("00052113","00200009327000008250400000","1111121067032230000","0163103100000010"); // -30ns, 35ns timing cut, NL3 EMC7
    cuts.AddCutPCMCalo("00052113","00200009327000008250400000","1111122067032230000","0163103100000010"); // -30ns, 35ns timing cut, NL4 EMC7
  } else if ( trainConfig == 505){ // EMCAL clusters
    cuts.AddCutPCMCalo("00085113","00200009327000008250400000","1111111067032230000","0163103100000010"); // -30ns, 35ns timing cut, NL1 EG2
    cuts.AddCutPCMCalo("00085113","00200009327000008250400000","1111112067032230000","0163103100000010"); // -30ns, 35ns timing cut, NL2 EG2
    cuts.AddCutPCMCalo("00085113","00200009327000008250400000","1111121067032230000","0163103100000010"); // -30ns, 35ns timing cut, NL3 EG2
    cuts.AddCutPCMCalo("00085113","00200009327000008250400000","1111122067032230000","0163103100000010"); // -30ns, 35ns timing cut, NL4 EG2
  } else if ( trainConfig == 506){ // EMCAL clusters
    cuts.AddCutPCMCalo("00083113","00200009327000008250400000","1111111067032230000","0163103100000010"); // -30ns, 35ns timing cut, NL1 EG1
    cuts.AddCutPCMCalo("00083113","00200009327000008250400000","1111121067032230000","0163103100000010"); // -30ns, 35ns timing cut, NL3 EG1
    cuts.AddCutPCMCalo("00083113","00200009327000008250400000","1111112067032230000","0163103100000010"); // -30ns, 35ns timing cut, NL2 EG1
    cuts.AddCutPCMCalo("00083113","00200009327000008250400000","1111122067032230000","0163103100000010"); // -30ns, 35ns timing cut, NL4 EG1
  } else if ( trainConfig == 507){ // EMCAL clusters
    cuts.AddCutPCMCalo("00085113","00200009327000008250400000","1111100017032230000","0163103100000010"); // 1000ns timing cut, no NL EG2
    cuts.AddCutPCMCalo("00083113","00200009327000008250400000","1111100017032230000","0163103100000010"); // 1000ns timing cut, no NL EG1

  //*****************************************************************************************************
  // pp 13 TeV DMC setup
  //*****************************************************************************************************
  } else if ( trainConfig==510){ //DCAL clusters
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","3885500017032230000","01631031000000d0");
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","3885500067032230000","0163103100000010"); // -30ns, 35ns timing cut, no NL INT7
  } else if ( trainConfig == 511){ // DCAL clusters
    cuts.AddCutPCMCalo("00055113","00200009327000008250400000","3885500017032230000","0163103100000010"); // 1000ns timing cut, no NL DMC7
    cuts.AddCutPCMCalo("00055113","00200009327000008250400000","3885500067032230000","0163103100000010"); // -30ns, 35ns timing cut, no NL DMC7
  } else if ( trainConfig == 512){ // DCAL clusters
    cuts.AddCutPCMCalo("00089113","00200009327000008250400000","3885500067032230000","0163103100000010"); // -30ns, 35ns timing cut, no NL DG2
    cuts.AddCutPCMCalo("0008b113","00200009327000008250400000","3885500067032230000","0163103100000010"); // -30ns, 35ns timing cut, no NL DG1
  } else if ( trainConfig == 513){ // DCAL clusters
    cuts.AddCutPCMCalo("00089113","00200009327000008250400000","3885500017032230000","0163103100000010"); // 1000ns timing cut, no NL DG2
    cuts.AddCutPCMCalo("0008b113","00200009327000008250400000","3885500017032230000","0163103100000010"); // 1000ns timing cut, no NL DG1

  //*************************************************************************************************
  // 13 TeV EMC setup high mult setups
  //*************************************************************************************************
  } else if ( trainConfig == 520 ) { // high mult triggers
    cuts.AddCutPCMCalo("00074113","00200009327000008250400000","1111100017032230000","0163103100000010"); // 1000ns timing cut, no NL V0HM
    cuts.AddCutPCMCalo("00074113","00200009327000008250400000","1111100067032230000","0163103100000010"); // -30ns, 35ns timing cut, no NL V0HM
    cuts.AddCutPCMCalo("00076113","00200009327000008250400000","1111100017032230000","0163103100000010"); // 1000ns timing cut, no NL V0HM + SPD1
    cuts.AddCutPCMCalo("00076113","00200009327000008250400000","1111100067032230000","0163103100000010"); // -30ns, 35ns timing cut, no NL V0HM + SPD1
  } else if ( trainConfig == 521 ) { // high mult triggers NL vars
    cuts.AddCutPCMCalo("00074113","00200009327000008250400000","1111111067032230000","0163103100000010"); // -30ns, 35ns timing cut, NL1 V0HM
    cuts.AddCutPCMCalo("00074113","00200009327000008250400000","1111112067032230000","0163103100000010"); // -30ns, 35ns timing cut, NL2 V0HM
    cuts.AddCutPCMCalo("00074113","00200009327000008250400000","1111121067032230000","0163103100000010"); // -30ns, 35ns timing cut, NL3 V0HM
    cuts.AddCutPCMCalo("00074113","00200009327000008250400000","1111122067032230000","0163103100000010"); // -30ns, 35ns timing cut, NL4 V0HM
  } else if ( trainConfig == 522 ) { // high mult triggers NL vars
    cuts.AddCutPCMCalo("00076113","00200009327000008250400000","1111111067032230000","0163103100000010"); // -30ns, 35ns timing cut, NL1 V0HM + SPD1
    cuts.AddCutPCMCalo("00076113","00200009327000008250400000","1111112067032230000","0163103100000010"); // -30ns, 35ns timing cut, NL2 V0HM + SPD1
    cuts.AddCutPCMCalo("00076113","00200009327000008250400000","1111121067032230000","0163103100000010"); // -30ns, 35ns timing cut, NL3 V0HM + SPD1
    cuts.AddCutPCMCalo("00076113","00200009327000008250400000","1111122067032230000","0163103100000010"); // -30ns, 35ns timing cut, NL4 V0HM + SPD1
  //*****************************************************************************************************
  // pp 13 TeV DMC setup high mult setups
  //*****************************************************************************************************
  } else if ( trainConfig == 530 ) { // high mult triggers
    cuts.AddCutPCMCalo("00074113","00200009327000008250400000","3885500017032230000","0163103100000010"); // 1000ns timing cut, no NL V0HM
    cuts.AddCutPCMCalo("00074113","00200009327000008250400000","3885500067032230000","0163103100000010"); // -30ns, 35ns timing cut, no NL V0HM
    cuts.AddCutPCMCalo("00076113","00200009327000008250400000","3885500017032230000","0163103100000010"); // 1000ns timing cut, no NL V0HM + SPD1
    cuts.AddCutPCMCalo("00076113","00200009327000008250400000","3885500067032230000","0163103100000010"); // -30ns, 35ns timing cut, no NL V0HM + SPD1
  //*************************************************************************************************
  // 13 TeV EDC setup low B periods
  //*************************************************************************************************
  } else if ( trainConfig == 540){ // EMCAL+DCal clusters
    cuts.AddCutPCMCalo("00010113","00200089327000008250400000","411790001f032230000","0163103100000010"); // 1000ns timing cut, no NL INT7
    cuts.AddCutPCMCalo("00010113","00200089327000008250400000","411790006f032230000","0163103100000010"); // -30ns, 35ns timing cut, no NL INT7
    cuts.AddCutPCMCalo("00010113","00200089327000008250400000","41179000af032230000","0163103100000010"); // -12.5ns, 12.5ns timing cut, no NL INT7
    cuts.AddCutPCMCalo("00010113","00200089327000008250400000","4117900060032230000","0163103100000010"); // -30ns, 35ns timing cut, no NL INT7
  } else if ( trainConfig == 541){ // EMCAL+DCal clusters
    cuts.AddCutPCMCalo("00010113","00200089327000008250400000","411791101f032230000","0163103100000010"); // 1000ns timing cut, INT7
    cuts.AddCutPCMCalo("00010113","00200089327000008250400000","411791106f032230000","0163103100000010"); // -30ns, 35ns timing cut, INT7
    cuts.AddCutPCMCalo("00010113","00200089327000008250400000","41179110af032230000","0163103100000010"); // -12.5ns, 12.5ns timing cut, INT7
    cuts.AddCutPCMCalo("00010113","00200089327000008250400000","4117911060032230000","0163103100000010"); // -30ns, 35ns timing cut, no NL INT7
  } else if ( trainConfig == 542){ // EMCAL+DCal clusters
    cuts.AddCutPCMCalo("00010113","00200089327000008250400000","411791201f032230000","0163103100000010"); // 1000ns timing cut, INT7
    cuts.AddCutPCMCalo("00010113","00200089327000008250400000","411791206f032230000","0163103100000010"); // -30ns, 35ns timing cut, INT7
    cuts.AddCutPCMCalo("00010113","00200089327000008250400000","41179120af032230000","0163103100000010"); // -12.5ns, 12.5ns timing cut, INT7
    cuts.AddCutPCMCalo("00010113","00200089327000008250400000","4117912060032230000","0163103100000010"); // -30ns, 35ns timing cut, INT7
  } else if ( trainConfig == 543){ // EMCAL+DCal clusters
    cuts.AddCutPCMCalo("00010113","00200089327000008250400000","411791106f032230000","0163103100000010"); // -30ns, 35ns timing cut TBNL + finetuning
    cuts.AddCutPCMCalo("00010113","00200089327000008250400000","411791206f032230000","0163103100000010"); // -30ns, 35ns timing cut TBNL + finetuning
    cuts.AddCutPCMCalo("00010113","00200089327000008250400000","411792106f032230000","0163103100000010"); // -30ns, 35ns timing cut TBNL + finetuning
    cuts.AddCutPCMCalo("00010113","00200089327000008250400000","411792206f032230000","0163103100000010"); // -30ns, 35ns timing cut TBNL + finetuning
  } else if ( trainConfig == 544){ // EMCAL+DCal clusters
    cuts.AddCutPCMCalo("00010113","00200089297000001280000000","411790006f032230000","0163103100000010"); // -30ns, 35ns timing cut, no NL INT7, open PCM cuts
    cuts.AddCutPCMCalo("00010113","00200089327000008250400000","411790006f032230000","0163103100000010"); // -30ns, 35ns timing cut, no NL INT7, std. PCM cuts
    cuts.AddCutPCMCalo("00010113","0020008932700000iih0400000","411790006f032230000","0163103100000010"); // -30ns, 35ns timing cut, no NL INT7, optimized PCM cuts
    cuts.AddCutPCMCalo("00010113","0020008932700000i280400000","411790006f032230000","0163103100000010"); // -30ns, 35ns timing cut, no NL INT7, partially optimized PCM cuts
    cuts.AddCutPCMCalo("00010113","00200089327000001ih0400000","411790006f032230000","0163103100000010"); // -30ns, 35ns timing cut, no NL INT7, partially optimized PCM cuts
  } else if ( trainConfig == 545){ // EMCAL+DCal clusters
    cuts.AddCutPCMCalo("00010113","0dm00089f9730000iih0404000","411791109fe30220000","0163103100000010"); // -30ns, 35ns timing cut TBNL
  } else if ( trainConfig == 546){  //   R Bins // weights 1
    cuts.AddCutPCMCalo("00010113", "0d200089f9730000iih0404000","411791109fe30220000", "0163103100000010"); // RBins    min = 5,      max = 180
    cuts.AddCutPCMCalo("00010113", "0dm00089f9730000iih0404000","411791109fe30220000", "0163103100000010"); // RBins    min = 5,      max = 180 without 55 -72
    cuts.AddCutPCMCalo("00010113", "0dh00089f9730000iih0404000","411791109fe30220000", "0163103100000010"); // RBins    min = 5,      max = 13
    cuts.AddCutPCMCalo("00010113", "0di00089f9730000iih0404000","411791109fe30220000", "0163103100000010"); // RBins    min = 13,     max = 33.5
  } else if ( trainConfig == 547){  //   R Bins // weights 1
    cuts.AddCutPCMCalo("00010113", "0dj00089f9730000iih0404000","411791109fe30220000", "0163103100000010"); // RBins    min = 33.5,   max = 55
    cuts.AddCutPCMCalo("00010113", "0dk00089f9730000iih0404000","411791109fe30220000", "0163103100000010"); // RBins    min = 55,     max = 72
    cuts.AddCutPCMCalo("00010113", "0dl00089f9730000iih0404000","411791109fe30220000", "0163103100000010"); // RBins    min = 72,     max = 95
    cuts.AddCutPCMCalo("00010113", "0dg00089f9730000iih0404000","411791109fe30220000", "0163103100000010"); // RBins    min = 95,     max = 180
  } else if ( trainConfig == 548){  //   R Bins // weights 2
    cuts.AddCutPCMCalo("00010113", "0d200089f9730000iih0404000","411791109fe30220000", "0163103100000010"); // RBins    min = 5,      max = 180
    cuts.AddCutPCMCalo("00010113", "0dm00089f9730000iih0404000","411791109fe30220000", "0163103100000010"); // RBins    min = 5,      max = 180 without 55 -72
    cuts.AddCutPCMCalo("00010113", "0dh00089f9730000iih0404000","411791109fe30220000", "0163103100000010"); // RBins    min = 5,      max = 13
    cuts.AddCutPCMCalo("00010113", "0di00089f9730000iih0404000","411791109fe30220000", "0163103100000010"); // RBins    min = 13,     max = 33.5
  } else if ( trainConfig == 549){  //   R Bins // weights 2
    cuts.AddCutPCMCalo("00010113", "0dj00089f9730000iih0404000","411791109fe30220000", "0163103100000010"); // RBins    min = 33.5,   max = 55
    cuts.AddCutPCMCalo("00010113", "0dk00089f9730000iih0404000","411791109fe30220000", "0163103100000010"); // RBins    min = 55,     max = 72
    cuts.AddCutPCMCalo("00010113", "0dl00089f9730000iih0404000","411791109fe30220000", "0163103100000010"); // RBins    min = 72,     max = 95
    cuts.AddCutPCMCalo("00010113", "0dg00089f9730000iih0404000","411791109fe30220000", "0163103100000010"); // RBins    min = 95,     max = 180
  } else if ( trainConfig == 550){  //   R Bins // weights 3
    cuts.AddCutPCMCalo("00010113", "0d200089f9730000iih0404000","411791109fe30220000", "0163103100000010"); // RBins    min = 5,      max = 180
    cuts.AddCutPCMCalo("00010113", "0dm00089f9730000iih0404000","411791109fe30220000", "0163103100000010"); // RBins    min = 5,      max = 180 without 55 -72
    cuts.AddCutPCMCalo("00010113", "0dh00089f9730000iih0404000","411791109fe30220000", "0163103100000010"); // RBins    min = 5,      max = 13
    cuts.AddCutPCMCalo("00010113", "0di00089f9730000iih0404000","411791109fe30220000", "0163103100000010"); // RBins    min = 13,     max = 33.5
  } else if ( trainConfig == 551){  //   R Bins // weights 3
    cuts.AddCutPCMCalo("00010113", "0dj00089f9730000iih0404000","411791109fe30220000", "0163103100000010"); // RBins    min = 33.5,   max = 55
    cuts.AddCutPCMCalo("00010113", "0dk00089f9730000iih0404000","411791109fe30220000", "0163103100000010"); // RBins    min = 55,     max = 72
    cuts.AddCutPCMCalo("00010113", "0dl00089f9730000iih0404000","411791109fe30220000", "0163103100000010"); // RBins    min = 72,     max = 95
    cuts.AddCutPCMCalo("00010113", "0dg00089f9730000iih0404000","411791109fe30220000", "0163103100000010"); // RBins    min = 95,     max = 180
  } else if ( trainConfig == 552){ // EMCAL clusters
    cuts.AddCutPCMCalo("00010113", "0dl00089f9730000iih0404000","411791109fe32220000","0163103100000010"); // OLD NCell >= 2
    cuts.AddCutPCMCalo("00010113", "0dl00089f9730000iih0404000","411791106fe30220000","0163103100000010"); // OLD Timing
    cuts.AddCutPCMCalo("00010113", "0dl00089f9730000iih0404000","411791109fe30230000","0163103100000010"); // OLD M02
    cuts.AddCutPCMCalo("00010113", "0dl00089f9730000iih0404000","4117911097e30220000","0163103100000010"); // OLD CPV

  //*************************************************************************************************
  // 13 TeV EDC setup
  //*************************************************************************************************
  } else if (trainConfig == 570){ //EMC + DCal HM trigger
    cuts.AddCutPCMCalo("00010113","00200089327000008250400000","4117900067032230000","01631031000000d0"); // -30ns, 35ns timing cut, MB trigger
    cuts.AddCutPCMCalo("00010113","00200089327000008250400000","4117900007032230000","01631031000000d0"); // no timing, MB trigger
    cuts.AddCutPCMCalo("00074113","00200089327000008250400000","4117900067032230000","01631031000000d0"); // -30ns, 35ns timing cut, no NL VOHM
    cuts.AddCutPCMCalo("00076113","00200089327000008250400000","4117900067032230000","01631031000000d0"); // -30ns, 35ns timing cut, no NL VOHM with SPD

  // R-bin variations with MBW
  } else if (trainConfig == 580) { // no smear - full R
    cuts.AddCutPCMCalo("00010113","0d200009f9730000dge0404000","111110106fg32230000","0163103100000010");
  } else if (trainConfig == 581) { // no smear - full R
    cuts.AddCutPCMCalo("00010113","0dh00009f9730000dge0404000","111110106fg32230000","0163103100000010");
    cuts.AddCutPCMCalo("00010113","0di00009f9730000dge0404000","111110106fg32230000","0163103100000010");
    cuts.AddCutPCMCalo("00010113","0dj00009f9730000dge0404000","111110106fg32230000","0163103100000010");
  } else if (trainConfig == 582) { // no smear - full R
    cuts.AddCutPCMCalo("00010113","0dk00009f9730000dge0404000","111110106fg32230000","0163103100000010");
    cuts.AddCutPCMCalo("00010113","0dl00009f9730000dge0404000","111110106fg32230000","0163103100000010");
    cuts.AddCutPCMCalo("00010113","0dg00009f9730000dge0404000","111110106fg32230000","0163103100000010");

  } else if (trainConfig == 583) { // no smear - full R
    cuts.AddCutPCMCalo("00052113","0d200009f9730000dge0404000","111110106fg32230000","0163103100000010");
  } else if (trainConfig == 584) { // no smear - full R
    cuts.AddCutPCMCalo("00052113","0dh00009f9730000dge0404000","111110106fg32230000","0163103100000010");
    cuts.AddCutPCMCalo("00052113","0di00009f9730000dge0404000","111110106fg32230000","0163103100000010");
    cuts.AddCutPCMCalo("00052113","0dj00009f9730000dge0404000","111110106fg32230000","0163103100000010");
  } else if (trainConfig == 585) { // no smear - full R
    cuts.AddCutPCMCalo("00052113","0dk00009f9730000dge0404000","111110106fg32230000","0163103100000010");
    cuts.AddCutPCMCalo("00052113","0dl00009f9730000dge0404000","111110106fg32230000","0163103100000010");
    cuts.AddCutPCMCalo("00052113","0dg00009f9730000dge0404000","111110106fg32230000","0163103100000010");

  } else if (trainConfig == 586) { // no smear - full R
    cuts.AddCutPCMCalo("00081113","0d200009f9730000dge0404000","111110106fg32230000","0163103100000010");
  } else if (trainConfig == 587) { // no smear - full R
    cuts.AddCutPCMCalo("00081113","0dh00009f9730000dge0404000","111110106fg32230000","0163103100000010");
    cuts.AddCutPCMCalo("00081113","0di00009f9730000dge0404000","111110106fg32230000","0163103100000010");
    cuts.AddCutPCMCalo("00081113","0dj00009f9730000dge0404000","111110106fg32230000","0163103100000010");
  } else if (trainConfig == 588) { // no smear - full R
    cuts.AddCutPCMCalo("00081113","0dk00009f9730000dge0404000","111110106fg32230000","0163103100000010");
    cuts.AddCutPCMCalo("00081113","0dl00009f9730000dge0404000","111110106fg32230000","0163103100000010");
    cuts.AddCutPCMCalo("00081113","0dg00009f9730000dge0404000","111110106fg32230000","0163103100000010");

  // R-bin variations with MBW and smearing
  } else if (trainConfig == 590) { // smearing - full R
    cuts.AddCutPCMCalo("00010113","0d200009f9730000dge0404000","111110106fg32230000","0163103100b00010");
  } else if (trainConfig == 591) { // smearing - full R
    cuts.AddCutPCMCalo("00010113","0dh00009f9730000dge0404000","111110106fg32230000","0163103100b00010");
    cuts.AddCutPCMCalo("00010113","0di00009f9730000dge0404000","111110106fg32230000","0163103100b00010");
    cuts.AddCutPCMCalo("00010113","0dj00009f9730000dge0404000","111110106fg32230000","0163103100b00010");
  } else if (trainConfig == 592) { // smearing - full R
    cuts.AddCutPCMCalo("00010113","0dk00009f9730000dge0404000","111110106fg32230000","0163103100b00010");
    cuts.AddCutPCMCalo("00010113","0dl00009f9730000dge0404000","111110106fg32230000","0163103100b00010");
    cuts.AddCutPCMCalo("00010113","0dg00009f9730000dge0404000","111110106fg32230000","0163103100b00010");

  } else if (trainConfig == 593) { // smearing - full R
    cuts.AddCutPCMCalo("00052113","0d200009f9730000dge0404000","111110106fg32230000","0163103100b00010");
  } else if (trainConfig == 594) { // smearing - full R
    cuts.AddCutPCMCalo("00052113","0dh00009f9730000dge0404000","111110106fg32230000","0163103100b00010");
    cuts.AddCutPCMCalo("00052113","0di00009f9730000dge0404000","111110106fg32230000","0163103100b00010");
    cuts.AddCutPCMCalo("00052113","0dj00009f9730000dge0404000","111110106fg32230000","0163103100b00010");
  } else if (trainConfig == 595) { // smearing - full R
    cuts.AddCutPCMCalo("00052113","0dk00009f9730000dge0404000","111110106fg32230000","0163103100b00010");
    cuts.AddCutPCMCalo("00052113","0dl00009f9730000dge0404000","111110106fg32230000","0163103100b00010");
    cuts.AddCutPCMCalo("00052113","0dg00009f9730000dge0404000","111110106fg32230000","0163103100b00010");

  } else if (trainConfig == 596) { // smearing - full R
    cuts.AddCutPCMCalo("00081113","0d200009f9730000dge0404000","111110106fg32230000","0163103100b00010");
  } else if (trainConfig == 597) { // smearing - full R
    cuts.AddCutPCMCalo("00081113","0dh00009f9730000dge0404000","111110106fg32230000","0163103100b00010");
    cuts.AddCutPCMCalo("00081113","0di00009f9730000dge0404000","111110106fg32230000","0163103100b00010");
    cuts.AddCutPCMCalo("00081113","0dj00009f9730000dge0404000","111110106fg32230000","0163103100b00010");
  } else if (trainConfig == 598) { // smearing - full R
    cuts.AddCutPCMCalo("00081113","0dk00009f9730000dge0404000","111110106fg32230000","0163103100b00010");
    cuts.AddCutPCMCalo("00081113","0dl00009f9730000dge0404000","111110106fg32230000","0163103100b00010");
    cuts.AddCutPCMCalo("00081113","0dg00009f9730000dge0404000","111110106fg32230000","0163103100b00010");
  //*****************************************************************************************************
  // pp 5 TeV DMC setup
  //*****************************************************************************************************
  } else if ( trainConfig == 600){ // DCAL clusters
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","3885500081041220000","0163103100000010"); //
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","3885500011041220000","0163103100000010"); //
  } else if ( trainConfig == 601){ // PCM-DCAL NonLinVariations
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","3885500081041220000","0163103100000010"); //
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","3885511081041220000","0163103100000010"); //
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","3885512081041220000","0163103100000010"); //
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","3885521081041220000","0163103100000010"); //
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","3885522081041220000","0163103100000010"); //
  } else if ( trainConfig == 610){ // DCAL clusters
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","38855000a7041230000","0163103100000010"); // -12.5 to 13 ns and M02 0.1-0.5
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","38855060a7041230000","0163103100000010"); // -12.5 to 13 ns and M02 0.1-0.5 TB NL
  } else if ( trainConfig == 611){ // DCAL clusters
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","38855610a7041230000","0163103100000010"); // -12.5 to 13 ns and M02 0.1-0.5
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","38855630a7041230000","0163103100000010"); // -12.5 to 13 ns and M02 0.1-0.5 TB NL
  } else if ( trainConfig == 612){ // DCAL clusters
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","38855620a7041230000","0163103100000010"); // -12.5 to 13 ns and M02 0.1-0.5
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","38855640a7041230000","0163103100000010"); // -12.5 to 13 ns and M02 0.1-0.5 TB NL
  } else if ( trainConfig == 660){ // INT7 DCAL standard cut but with E/p TM veto
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","3885511087041220000","0163103100000010"); // std INT7
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","388551108c041220000","0163103100000010"); // fEOverPMax = 9e9
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","388551108d041220000","0163103100000010"); // fEOverPMax = 3.0
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","388551108e041220000","0163103100000010"); // fEOverPMax = 2.0
  } else if ( trainConfig == 661){ // INT7 DCAL standard cut but with E/p TM veto
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","388551108f041220000","0163103100000010"); // fEOverPMax = 1.75
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","388551108g041220000","0163103100000010"); // fEOverPMax = 1.5
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","388551108h041220000","0163103100000010"); // fEOverPMax = 1.25
  } else if ( trainConfig == 662){ // DMC7 DCAL standard cut but with E/p TM veto
    cuts.AddCutPCMCalo("00055113","00200009327000008250400000","3885511087041220000","0163103100000010"); // std DMC7
    cuts.AddCutPCMCalo("00055113","00200009327000008250400000","388551108c041220000","0163103100000010"); // fEOverPMax = 9e9
    cuts.AddCutPCMCalo("00055113","00200009327000008250400000","388551108d041220000","0163103100000010"); // fEOverPMax = 3.0
    cuts.AddCutPCMCalo("00055113","00200009327000008250400000","388551108e041220000","0163103100000010"); // fEOverPMax = 2.0
  } else if ( trainConfig == 663){ // DMC7 DCAL standard cut but with E/p TM veto
    cuts.AddCutPCMCalo("00055113","00200009327000008250400000","388551108f041220000","0163103100000010"); // fEOverPMax = 1.75
    cuts.AddCutPCMCalo("00055113","00200009327000008250400000","388551108g041220000","0163103100000010"); // fEOverPMax = 1.5
    cuts.AddCutPCMCalo("00055113","00200009327000008250400000","388551108h041220000","0163103100000010"); // fEOverPMax = 1.25
  } else if ( trainConfig == 664){ // DG1 DCAL standard cut but with E/p TM veto
    cuts.AddCutPCMCalo("0008b113","00200009327000008250400000","3885511087041220000","0163103100000010"); // std DG1
    cuts.AddCutPCMCalo("0008b113","00200009327000008250400000","388551108c041220000","0163103100000010"); // fEOverPMax = 9e9
    cuts.AddCutPCMCalo("0008b113","00200009327000008250400000","388551108d041220000","0163103100000010"); // fEOverPMax = 3.0
    cuts.AddCutPCMCalo("0008b113","00200009327000008250400000","388551108e041220000","0163103100000010"); // fEOverPMax = 2.0
  } else if ( trainConfig == 665){ // DG1 DCAL standard cut but with E/p TM veto
    cuts.AddCutPCMCalo("0008b113","00200009327000008250400000","388551108f041220000","0163103100000010"); // fEOverPMax = 1.75
    cuts.AddCutPCMCalo("0008b113","00200009327000008250400000","388551108g041220000","0163103100000010"); // fEOverPMax = 1.5
    cuts.AddCutPCMCalo("0008b113","00200009327000008250400000","388551108h041220000","0163103100000010"); // fEOverPMax = 1.25
  } else if ( trainConfig == 666){ // DG2 DCAL standard cut but with E/p TM veto
    cuts.AddCutPCMCalo("00089113","00200009327000008250400000","3885511087041220000","0163103100000010"); // std DG2
    cuts.AddCutPCMCalo("00089113","00200009327000008250400000","388551108c041220000","0163103100000010"); // fEOverPMax = 9e9
    cuts.AddCutPCMCalo("00089113","00200009327000008250400000","388551108d041220000","0163103100000010"); // fEOverPMax = 3.0
    cuts.AddCutPCMCalo("00089113","00200009327000008250400000","388551108e041220000","0163103100000010"); // fEOverPMax = 2.0
  } else if ( trainConfig == 667){ // DG2 DCAL standard cut but with E/p TM veto
    cuts.AddCutPCMCalo("00089113","00200009327000008250400000","388551108f041220000","0163103100000010"); // fEOverPMax = 1.75
    cuts.AddCutPCMCalo("00089113","00200009327000008250400000","388551108g041220000","0163103100000010"); // fEOverPMax = 1.5
    cuts.AddCutPCMCalo("00089113","00200009327000008250400000","388551108h041220000","0163103100000010"); // fEOverPMax = 1.25

  } else if ( trainConfig == 680){ // INT7 DCAL standard cut with CALO_+CALOFAST readout
    cuts.AddCutPCMCalo("000a0113","00200009327000008250400000","3885511080041220000","0163103100000010"); // std INT7
    cuts.AddCutPCMCalo("000a6113","00200009327000008250400000","3885511080041220000","0163103100000010"); // std DMC7
    cuts.AddCutPCMCalo("000a7113","00200009327000008250400000","3885511080041220000","0163103100000010"); // std DG2
    cuts.AddCutPCMCalo("000a8113","00200009327000008250400000","3885511080041220000","0163103100000010"); // std DG1
  } else if ( trainConfig == 681){ // INT7 DCAL standard cut with standard readout
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","3885511080041220000","0163103100000010"); // std INT7
    cuts.AddCutPCMCalo("00055113","00200009327000008250400000","3885511080041220000","0163103100000010"); // std DMC7
    cuts.AddCutPCMCalo("00089113","00200009327000008250400000","3885511080041220000","0163103100000010"); // std DG2
    cuts.AddCutPCMCalo("0008b113","00200009327000008250400000","3885511080041220000","0163103100000010"); // std DG1
  //*************************************************************************************************
  // 2.76 TeV EMC setup variations
  //*************************************************************************************************
  } else if ( trainConfig == 700){ // 2011 variations
    cuts.AddCutPCMCalo("00003113","00200009327000008250400000","1111121057032230000","0163103100000010");
    cuts.AddCutPCMCalo("00003113","00200009327000008250400000","1111121067022230000","0163103100000010"); //0.6 GeV/c
    cuts.AddCutPCMCalo("00003113","00200009327000008250400000","1111121067032230000","0163103100000010"); //0.7 GeV/c default
  } else if ( trainConfig == 701){ // 2011 variations
    cuts.AddCutPCMCalo("00003113","00200009327000008250400000","1111121067042230000","0163103100000010"); //0.8 GeV/c
    cuts.AddCutPCMCalo("00003113","00200009327000008250400000","1111121067052230000","0163103100000010"); //0.9 GeV/c
    cuts.AddCutPCMCalo("00003113","00200009327000008250400000","1113121057032230000","0163103100000010"); //only modules with TRD infront
  } else if ( trainConfig == 702){ // 2011 variations
    cuts.AddCutPCMCalo("00003113","00200009327000008250400000","1111221057032230000","0163103100000010"); //no modules with TRD infront
    cuts.AddCutPCMCalo("00003113","00200009327000008250400000","1111100057032230000","0163103100000010"); // NonLinearity none
    cuts.AddCutPCMCalo("00003113","00200009327000008250400000","1111101057032230000","0163103100000010"); // NonLinearity kSDMv5
  } else if ( trainConfig == 703){ // 2011 variations
    cuts.AddCutPCMCalo("00003113","00200009327000008250400000","1111122057032230000","0163103100000010"); // NonLinearity CMF
    cuts.AddCutPCMCalo("00003113","00200009327000008250400000","1111111057032230000","0163103100000010"); // NonLinearity CCRF
    cuts.AddCutPCMCalo("00003113","00200009327000008250400000","1111112057032230000","0163103100000010"); // NonLinearity CRF
  } else if ( trainConfig == 704){ // 2011 variations
    cuts.AddCutPCMCalo("00003113","00200009327000008250400000","1111121057031230000","0163103100000010"); //n cells >= 1
    cuts.AddCutPCMCalo("00003113","00200009327000008250400000","1111121057033230000","0163103100000010"); //n cells >= 3
    cuts.AddCutPCMCalo("00003113","00200009327000008250400000","1111121057032200000","0163103100000010"); //no max M02 cut
  } else if ( trainConfig == 705){ // 2011 variations
    cuts.AddCutPCMCalo("00003113","00200009327000008250400000","1111121057032250000","0163103100000010"); //M02 < 0.3
    cuts.AddCutPCMCalo("00003113","00200009327000008250400000","1111121057032260000","0163103100000010"); //M02 < 0.27
    cuts.AddCutPCMCalo("00003113","00200009327000008250400000","1111121056032230000","0163103100000010"); // track matching variations
  } else if ( trainConfig == 706){ // 2011 variations
    cuts.AddCutPCMCalo("00003113","00200009327000008250400000","1111121057032230000","0163103100000010"); //
    cuts.AddCutPCMCalo("00003113","00200009327000008250400000","1111121058032230000","0163103100000010"); //
    cuts.AddCutPCMCalo("00003113","00200009327000008250400000","1111121059032230000","0163103100000010"); //
  } else if ( trainConfig == 707){ // 2011 variations
    cuts.AddCutPCMCalo("00003113","00200009327000008250400000","1111121053032230000","0163103100000010"); //
    cuts.AddCutPCMCalo("00003113","00200009327000008250400000","1111121050032230000","0163103100000010"); //
    cuts.AddCutPCMCalo("00003113","00200009327000008250400000","1111121067032230000","0163103100000010"); // time 50ns
  } else if ( trainConfig == 708){ // 2011 variations
    cuts.AddCutPCMCalo("00003113","00200009327000008250400000","1111121047032230000","0163103100000010"); // time 100ns
    cuts.AddCutPCMCalo("00003113","00200009327000008250400000","1111121027032230000","0163103100000010"); // time 500ns
    cuts.AddCutPCMCalo("00003113","00200009327000008250400000","1111121057232230000","0163103100000010"); //
  } else if ( trainConfig == 710){ // 2011 variations
    cuts.AddCutPCMCalo("00003113","00200009327000008250400000","1111123057032230000","0163103100000010"); // NonLinearity kSDMv5
    cuts.AddCutPCMCalo("00003113","00200009327000008250400000","1111113057032230000","0163103100000010"); // NonLinearity LHC12 kTestBeamv2+ConvCalo
    cuts.AddCutPCMCalo("00003113","00200009327000008250400000","1111114057032230000","0163103100000010"); // NonLinearity LHC12 kTestBeamv2+Calo
  } else if ( trainConfig == 711){ // 2011 variations
    cuts.AddCutPCMCalo("00003113","00200009327000008250400000","1111124057032230000","0163103100000010"); // NonLinearity LHC12 kTestBeamv2+Calo
    cuts.AddCutPCMCalo("00003113","00200009227000008250400000","1111121057032230000","0163103100000010"); // dEdx e -3, 5
    cuts.AddCutPCMCalo("00003113","00200009127000008250400000","1111121057032230000","0163103100000010"); // dEdx e -5, 5
  } else if ( trainConfig == 712){ // 2011 variations
    cuts.AddCutPCMCalo("00003113","00200009357000008250400000","1111121057032230000","0163103100000010"); // dEdx pi 2
    cuts.AddCutPCMCalo("00003113","00200009317000008250400000","1111121057032230000","0163103100000010"); // dEdx pi 0
    cuts.AddCutPCMCalo("00003113","00200009387300008250400000","1111121057032230000","0163103100000010"); // dEdx pi 2 high 1 (> 3.5 GeV)
  } else if ( trainConfig == 713){ // 2011 variations
    cuts.AddCutPCMCalo("00003113","00200009317300008250400000","1111121057032230000","0163103100000010"); // dEdx pi: 0: 0.4-3.5, -10: 3.5 ->
    cuts.AddCutPCMCalo("00003113","00200009327300008250400000","1111121057032230000","0163103100000010"); // dEdx pi: 1: 0.4-3.5, -10: 3.5 ->
    cuts.AddCutPCMCalo("00003113","00200009325000008250400000","1111121057032230000","0163103100000010"); // dEdx pi: 1: 0.3 ->
  } else if ( trainConfig == 714){ // 2011 variations
    cuts.AddCutPCMCalo("00003113","00200009320000008250400000","1111121057032230000","0163103100000010"); // dEdx pi: 1: 0.5 ->
    cuts.AddCutPCMCalo("00003113","00200009327600008250400000","1111121057032230000","0163103100000010"); // dEdx pi: 1: 0.4-2, -10: 2. ->
    cuts.AddCutPCMCalo("00003113","00200009327400008250400000","1111121057032230000","0163103100000010"); // dEdx pi: 1: 0.4-3, -10: 3. ->
  } else if ( trainConfig == 715){ // 2011 variations
    cuts.AddCutPCMCalo("00003113","00200009315600008250400000","1111121057032230000","0163103100000010"); // dEdx pi: 0: 0.3-2, -10: 2. ->
    cuts.AddCutPCMCalo("00003113","00200009367400008250400000","1111121057032230000","0163103100000010"); // dEdx pi: 2: 0.4-3, 0.5: 3. ->
    cuts.AddCutPCMCalo("00003113","00200009347400008250400000","1111121057032230000","0163103100000010"); // dEdx pi: 3: 0.4-3, 1: 3. ->
  } else if ( trainConfig == 716){ // 2011 variations
    cuts.AddCutPCMCalo("00003113","00200009327000009250400000","1111121057032230000","0163103100000010"); // qt 2D 0.03
    cuts.AddCutPCMCalo("00003113","00200009327000003250400000","1111121057032230000","0163103100000010"); // qt 1D 0.05
    cuts.AddCutPCMCalo("00003113","00200009327000002250400000","1111121057032230000","0163103100000010"); // qt 1D 0.07
  } else if ( trainConfig == 717){ // 2011 variations
    cuts.AddCutPCMCalo("00003113","00200049327000008250400000","1111121057032230000","0163103100000010"); // single pt > 0.075
    cuts.AddCutPCMCalo("00003113","00200019327000008250400000","1111121057032230000","0163103100000010"); // single pt > 0.1
    cuts.AddCutPCMCalo("00003113","00200009327000008850400000","1111121057032230000","0163103100000010"); // 2D psi pair chi2 var
  } else if ( trainConfig == 718){ // 2011 variations
    cuts.AddCutPCMCalo("00003113","00200009327000008260400000","1111121057032230000","0163103100000010"); // 2D psi pair chi2 var
    cuts.AddCutPCMCalo("00003113","00200009327000008860400000","1111121057032230000","0163103100000010"); // 2D psi pair chi2 var
    cuts.AddCutPCMCalo("00003113","00200009327000008280400000","1111121057032230000","0163103100000010"); // 2D psi pair chi2 var
  } else if ( trainConfig == 719){ // 2011 variations
    cuts.AddCutPCMCalo("00003113","00200009327000008880400000","1111121057032230000","0163103100000010"); // 2D psi pair chi2 var
    cuts.AddCutPCMCalo("00003113","00200006327000008250400000","1111121057032230000","0163103100000010"); // min TPC cl > 0.7
    cuts.AddCutPCMCalo("00003113","00200008327000008250400000","1111121057032230000","0163103100000010"); // min TPC cl > 0.35
  } else if ( trainConfig == 720){ // 2011 variations
    cuts.AddCutPCMCalo("00003113","00200009327000008250400000","1111121057032230000","0163106100000010"); // alpha < 0.8
    cuts.AddCutPCMCalo("00003113","00200009327000008250400000","1111121057032230000","0163105100000010"); // alpha < 0.75
    cuts.AddCutPCMCalo("00003113","00202209327000008250400000","1111121057032230000","0163103100000010"); // restrict acceptance to EMCAL loose
  } else if ( trainConfig == 721){ // 2011 variations
    cuts.AddCutPCMCalo("00003113","00204409327000008250400000","1111121057032230000","0163103100000010"); // restrict acceptance to EMCAL tight
    cuts.AddCutPCMCalo("00003113","00200009327000008250400000","1111121057032230000","0163103100000010"); //
    cuts.AddCutPCMCalo("00003113","00200009327000008250401000","1111121057032230000","0163103100000010"); //
  } else if ( trainConfig == 722){ // 2011 variations
    cuts.AddCutPCMCalo("00003113","00200009327000008250402000","1111121057032230000","0163103100000010"); //
    cuts.AddCutPCMCalo("00003113","00200009327000008250403000","1111121057032230000","0163103100000010"); //

  // 8 TeV pp variations with new PCM cut
  } else if (trainConfig == 730) { // PCM variations
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0400000","111113206f032230000","0163103100000010"); //New standard cut for 8TeV analysis for RpA
    cuts.AddCutPCMCalo("00010013","00200009f9730000dge0400000","111113206f032230000","0163103100000010"); // no SPD pileup cut
    cuts.AddCutPCMCalo("00010113","00100009f9730000dge0400000","111113206f032230000","0163103100000010"); // R cut 2.8 -180 cm
    cuts.AddCutPCMCalo("00010113","00500009f9730000dge0400000","111113206f032230000","0163103100000010"); // R cut 10. -180 cm
  } else if (trainConfig == 731) {
    cuts.AddCutPCMCalo("00010113","00200069f9730000dge0400000","111113206f032230000","0163103100000010"); // min pT 40 MeV
    cuts.AddCutPCMCalo("00010113","00200049f9730000dge0400000","111113206f032230000","0163103100000010"); // min pT 75 MeV
    cuts.AddCutPCMCalo("00010113","00200019f9730000dge0400000","111113206f032230000","0163103100000010"); // min pT 100MeV
  } else if (trainConfig == 732) {
    cuts.AddCutPCMCalo("00010113","00200068f9730000dge0400000","111113206f032230000","0163103100000010"); // TPC cluster 35%
    cuts.AddCutPCMCalo("00010113","00200066f9730000dge0400000","111113206f032230000","0163103100000010"); // TPC cluster 70%
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0600000","111113206f032230000","0163103100000010"); // cosPA 0.9
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0300000","111113206f032230000","0163103100000010"); // cosPA 0.75
  } else if (trainConfig == 733) {
    cuts.AddCutPCMCalo("00010113","0020000939730000dge0400000","111113206f032230000","0163103100000010"); // nsig electron   -4,5
    cuts.AddCutPCMCalo("00010113","0020000969730000dge0400000","111113206f032230000","0163103100000010"); // nsig electron -2.5,4
    cuts.AddCutPCMCalo("00010113","00200009f5730000dge0400000","111113206f032230000","0163103100000010"); // nsig pion 2,-10
    cuts.AddCutPCMCalo("00010113","00200009f1730000dge0400000","111113206f032230000","0163103100000010"); // nsig pion 0,-10
  } else if (trainConfig == 734) {
    cuts.AddCutPCMCalo("00010113","00200009f9030000dge0400000","111113206f032230000","0163103100000010"); // pion nsig min mom 0.50 GeV/c
    cuts.AddCutPCMCalo("00010113","00200009f9630000dge0400000","111113206f032230000","0163103100000010"); // pion nsig min mom 0.25 GeV/c
    cuts.AddCutPCMCalo("00010113","00200009f9760000dge0400000","111113206f032230000","0163103100000010"); // pion nsig max mom 2.00 GeV/c
    cuts.AddCutPCMCalo("00010113","00200009f9710000dge0400000","111113206f032230000","0163103100000010"); // pion nsig max mom 5.00 GeV/c
  } else if (trainConfig == 735) {
    cuts.AddCutPCMCalo("00010113","00200009f9730000age0400000","111113206f032230000","0163103100000010"); // qT max 0.040, qt pt max 0.11
    cuts.AddCutPCMCalo("00010113","00200009f9730000ege0400000","111113206f032230000","0163103100000010"); // qT max 0.060, qt pt max 0.14
    cuts.AddCutPCMCalo("00010113","00200009f9730000fge0400000","111113206f032230000","0163103100000010"); // qT max 0.070, qt pt max 0.16
  } else if (trainConfig == 736) {
    cuts.AddCutPCMCalo("00010113","00200009f9730000d1e0400000","111113206f032230000","0163103100000010"); // chi2 50 no chi2 dep.
    cuts.AddCutPCMCalo("00010113","00200009f9730000dfe0400000","111113206f032230000","0163103100000010"); // chi2 50 chi2 dep -0.065
    cuts.AddCutPCMCalo("00010113","00200009f9730000dhe0400000","111113206f032230000","0163103100000010"); // chi2 50 chi2 dep -0.050
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0404000","111113206f032230000","0163103100000010"); // reject close v0
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0406000","111113206f032230000","0163103100000010"); // double count with open angle 0.04
  } else if (trainConfig == 737) {
    cuts.AddCutPCMCalo("00010113","00200009f9730000dgd0400000","111113206f032230000","0163103100000010"); // Psi pair 0.15 dep
    cuts.AddCutPCMCalo("00010113","00200009f9730000dgf0400000","111113206f032230000","0163103100000010"); // Psi pair 0.20 dep
    cuts.AddCutPCMCalo("00010113","00200009f9730000dgg0400000","111113206f032230000","0163103100000010"); // Psi pair 0.30 dep
  } else if (trainConfig == 738) {
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0400000","111113206f032230000","0263103100000010"); // variation BG scheme track mult
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0400000","111113206f032230000","0163107100000010"); // alpha meson 0.85
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0400000","111113206f032230000","0163105100000010"); // alpha meson 0.75
    cuts.AddCutPCMCalo("00010113","00200009227300008250404000","111113206f032230000","0163103100000010"); // old cuts (run1)
  } else if (trainConfig == 739) {
    cuts.AddCutPCMCalo("00010213","00200009f9730000dge0400000","111113206f032230000","0163103100000010"); //same as std + maximum past future rejection
    cuts.AddCutPCMCalo("00010513","00200009f9730000dge0400000","111113206f032230000","0163103100000010"); //same as std + medium past future rejection
  } else if (trainConfig == 740) { // CALO variations
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0400000","111113206f022230000","0163103100000010"); // 0.6 GeV/c
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0400000","111113206f042230000","0163103100000010"); // 0.8 GeV/c
  } else if (trainConfig == 741) {
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0400000","111113206f031230000","0163103100000010"); // n cells >= 1
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0400000","111113206f033230000","0163103100000010"); // n cells >= 3
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0400000","111113206f032200000","0163103100000010"); // no max M02 cut
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0400000","111113206f032250000","0163103100000010"); // M02 < 0.3
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0400000","111113206f0322k0000","0163103100000010"); // M02, pT-dep
  } else if (trainConfig == 742) {
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0400000","111113206e032230000","0163103100000010"); // TM var e/p 2.0
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0400000","111113206g032230000","0163103100000010"); // TM var e/p 1.5
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0400000","111113206h032230000","0163103100000010"); // TM var e/p 1.25
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0400000","1111132067032230000","0163103100000010"); // TM var no veto
  } else if (trainConfig == 743) {
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0400000","111113306f032230000","0163103100000010"); // NL 32
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0400000","111113406f032230000","0163103100000010"); // NL 32
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0400000","111113806f032230000","0163103100000010"); // NL 33
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0400000","111113906f032230000","0163103100000010"); // NL 34
  } else if (trainConfig == 744) {
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0400000","111113205f032230000","0163103100000010"); // 50ns
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0400000","111113204f032230000","0163103100000010"); // 100ns
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0400000","111113207f032230000","0163103100000010"); // 30ns

  } else if (trainConfig == 750) { // PCM variations
    cuts.AddCutPCMCalo("00052113","00200009f9730000dge0400000","111113206f032230000","0163103100000010"); //New standard cut for 8TeV analysis for RpA
    cuts.AddCutPCMCalo("00052013","00200009f9730000dge0400000","111113206f032230000","0163103100000010"); // no SPD pileup cut
    cuts.AddCutPCMCalo("00052113","00100009f9730000dge0400000","111113206f032230000","0163103100000010"); // R cut 2.8 -180 cm
    cuts.AddCutPCMCalo("00052113","00500009f9730000dge0400000","111113206f032230000","0163103100000010"); // R cut 10. -180 cm
  } else if (trainConfig == 751) {
    cuts.AddCutPCMCalo("00052113","00200069f9730000dge0400000","111113206f032230000","0163103100000010"); // min pT 40 MeV
    cuts.AddCutPCMCalo("00052113","00200049f9730000dge0400000","111113206f032230000","0163103100000010"); // min pT 75 MeV
    cuts.AddCutPCMCalo("00052113","00200019f9730000dge0400000","111113206f032230000","0163103100000010"); // min pT 100MeV
  } else if (trainConfig == 752) {
    cuts.AddCutPCMCalo("00052113","00200068f9730000dge0400000","111113206f032230000","0163103100000010"); // TPC cluster 35%
    cuts.AddCutPCMCalo("00052113","00200066f9730000dge0400000","111113206f032230000","0163103100000010"); // TPC cluster 70%
    cuts.AddCutPCMCalo("00052113","00200009f9730000dge0600000","111113206f032230000","0163103100000010"); // cosPA 0.9
    cuts.AddCutPCMCalo("00052113","00200009f9730000dge0300000","111113206f032230000","0163103100000010"); // cosPA 0.75
  } else if (trainConfig == 753) {
    cuts.AddCutPCMCalo("00052113","0020000939730000dge0400000","111113206f032230000","0163103100000010"); // nsig electron   -4,5
    cuts.AddCutPCMCalo("00052113","0020000969730000dge0400000","111113206f032230000","0163103100000010"); // nsig electron -2.5,4
    cuts.AddCutPCMCalo("00052113","00200009f5730000dge0400000","111113206f032230000","0163103100000010"); // nsig pion 2,-10
    cuts.AddCutPCMCalo("00052113","00200009f1730000dge0400000","111113206f032230000","0163103100000010"); // nsig pion 0,-10
  } else if (trainConfig == 754) {
    cuts.AddCutPCMCalo("00052113","00200009f9030000dge0400000","111113206f032230000","0163103100000010"); // pion nsig min mom 0.50 GeV/c
    cuts.AddCutPCMCalo("00052113","00200009f9630000dge0400000","111113206f032230000","0163103100000010"); // pion nsig min mom 0.25 GeV/c
    cuts.AddCutPCMCalo("00052113","00200009f9760000dge0400000","111113206f032230000","0163103100000010"); // pion nsig max mom 2.00 GeV/c
    cuts.AddCutPCMCalo("00052113","00200009f9710000dge0400000","111113206f032230000","0163103100000010"); // pion nsig max mom 5.00 GeV/c
  } else if (trainConfig == 755) {
    cuts.AddCutPCMCalo("00052113","00200009f9730000age0400000","111113206f032230000","0163103100000010"); // qT max 0.040, qt pt max 0.11
    cuts.AddCutPCMCalo("00052113","00200009f9730000ege0400000","111113206f032230000","0163103100000010"); // qT max 0.060, qt pt max 0.14
    cuts.AddCutPCMCalo("00052113","00200009f9730000fge0400000","111113206f032230000","0163103100000010"); // qT max 0.070, qt pt max 0.16
  } else if (trainConfig == 756) {
    cuts.AddCutPCMCalo("00052113","00200009f9730000d1e0400000","111113206f032230000","0163103100000010"); // chi2 50 no chi2 dep.
    cuts.AddCutPCMCalo("00052113","00200009f9730000dfe0400000","111113206f032230000","0163103100000010"); // chi2 50 chi2 dep -0.065
    cuts.AddCutPCMCalo("00052113","00200009f9730000dhe0400000","111113206f032230000","0163103100000010"); // chi2 50 chi2 dep -0.050
    cuts.AddCutPCMCalo("00052113","00200009f9730000dge0404000","111113206f032230000","0163103100000010"); // reject close v0
    cuts.AddCutPCMCalo("00052113","00200009f9730000dge0406000","111113206f032230000","0163103100000010"); // double count with open angle 0.04
  } else if (trainConfig == 757) {
    cuts.AddCutPCMCalo("00052113","00200009f9730000dgd0400000","111113206f032230000","0163103100000010"); // Psi pair 0.15 dep
    cuts.AddCutPCMCalo("00052113","00200009f9730000dgf0400000","111113206f032230000","0163103100000010"); // Psi pair 0.20 dep
    cuts.AddCutPCMCalo("00052113","00200009f9730000dgg0400000","111113206f032230000","0163103100000010"); // Psi pair 0.30 dep
  } else if (trainConfig == 758) {
    cuts.AddCutPCMCalo("00052113","00200009f9730000dge0400000","111113206f032230000","0263103100000010"); // variation BG scheme track mult
    cuts.AddCutPCMCalo("00052113","00200009f9730000dge0400000","111113206f032230000","0163107100000010"); // alpha meson 0.85
    cuts.AddCutPCMCalo("00052113","00200009f9730000dge0400000","111113206f032230000","0163105100000010"); // alpha meson 0.75
    cuts.AddCutPCMCalo("00052113","00200009227300008250404000","111113206f032230000","0163103100000010"); // old cuts (run1)
  } else if (trainConfig == 759) {
    cuts.AddCutPCMCalo("00052213","00200009f9730000dge0400000","111113206f032230000","0163103100000010"); //same as std + maximum past future rejection
    cuts.AddCutPCMCalo("00052513","00200009f9730000dge0400000","111113206f032230000","0163103100000010"); //same as std + medium past future rejection
  } else if (trainConfig == 760) { // CALO variations
    cuts.AddCutPCMCalo("00052113","00200009f9730000dge0400000","111113206f022230000","0163103100000010"); // 0.6 GeV/c
    cuts.AddCutPCMCalo("00052113","00200009f9730000dge0400000","111113206f042230000","0163103100000010"); // 0.8 GeV/c
  } else if (trainConfig == 761) {
    cuts.AddCutPCMCalo("00052113","00200009f9730000dge0400000","111113206f031230000","0163103100000010"); // n cells >= 1
    cuts.AddCutPCMCalo("00052113","00200009f9730000dge0400000","111113206f033230000","0163103100000010"); // n cells >= 3
    cuts.AddCutPCMCalo("00052113","00200009f9730000dge0400000","111113206f032200000","0163103100000010"); // no max M02 cut
    cuts.AddCutPCMCalo("00052113","00200009f9730000dge0400000","111113206f032250000","0163103100000010"); // M02 < 0.3
    cuts.AddCutPCMCalo("00052113","00200009f9730000dge0400000","111113206f0322k0000","0163103100000010"); // M02, pT-dep
  } else if (trainConfig == 762) {
    cuts.AddCutPCMCalo("00052113","00200009f9730000dge0400000","111113206e032230000","0163103100000010"); // TM var e/p 2.0
    cuts.AddCutPCMCalo("00052113","00200009f9730000dge0400000","111113206g032230000","0163103100000010"); // TM var e/p 1.5
    cuts.AddCutPCMCalo("00052113","00200009f9730000dge0400000","111113206h032230000","0163103100000010"); // TM var e/p 1.25
    cuts.AddCutPCMCalo("00052113","00200009f9730000dge0400000","1111132067032230000","0163103100000010"); // TM var no veto
  } else if (trainConfig == 763) {
    cuts.AddCutPCMCalo("00052113","00200009f9730000dge0400000","111113306f032230000","0163103100000010"); // NL 32
    cuts.AddCutPCMCalo("00052113","00200009f9730000dge0400000","111113406f032230000","0163103100000010"); // NL 33
    cuts.AddCutPCMCalo("00052113","00200009f9730000dge0400000","111113806f032230000","0163103100000010"); // NL 34
    cuts.AddCutPCMCalo("00052113","00200009f9730000dge0400000","111113906f032230000","0163103100000010"); // NL 34
  } else if (trainConfig == 764) {
    cuts.AddCutPCMCalo("00052113","00200009f9730000dge0400000","111113205f032230000","0163103100000010"); // 50ns
    cuts.AddCutPCMCalo("00052113","00200009f9730000dge0400000","111113204f032230000","0163103100000010"); // 100ns
    cuts.AddCutPCMCalo("00052113","00200009f9730000dge0400000","111113207f032230000","0163103100000010"); // 30ns

  } else if (trainConfig == 770) { // PCM variations
    cuts.AddCutPCMCalo("00081113","00200009f9730000dge0400000","111113206f032230000","0163103100000010"); //New standard cut for 8TeV analysis for RpA
    cuts.AddCutPCMCalo("00081013","00200009f9730000dge0400000","111113206f032230000","0163103100000010"); // no SPD pileup cut
    cuts.AddCutPCMCalo("00081113","00100009f9730000dge0400000","111113206f032230000","0163103100000010"); // R cut 2.8 -180 cm
    cuts.AddCutPCMCalo("00081113","00500009f9730000dge0400000","111113206f032230000","0163103100000010"); // R cut 10. -180 cm
  } else if (trainConfig == 771) {
    cuts.AddCutPCMCalo("00081113","00200069f9730000dge0400000","111113206f032230000","0163103100000010"); // min pT 40 MeV
    cuts.AddCutPCMCalo("00081113","00200049f9730000dge0400000","111113206f032230000","0163103100000010"); // min pT 75 MeV
    cuts.AddCutPCMCalo("00081113","00200019f9730000dge0400000","111113206f032230000","0163103100000010"); // min pT 100MeV
  } else if (trainConfig == 772) {
    cuts.AddCutPCMCalo("00081113","00200068f9730000dge0400000","111113206f032230000","0163103100000010"); // TPC cluster 35%
    cuts.AddCutPCMCalo("00081113","00200066f9730000dge0400000","111113206f032230000","0163103100000010"); // TPC cluster 70%
    cuts.AddCutPCMCalo("00081113","00200009f9730000dge0600000","111113206f032230000","0163103100000010"); // cosPA 0.9
    cuts.AddCutPCMCalo("00081113","00200009f9730000dge0300000","111113206f032230000","0163103100000010"); // cosPA 0.75
  } else if (trainConfig == 773) {
    cuts.AddCutPCMCalo("00081113","0020000939730000dge0400000","111113206f032230000","0163103100000010"); // nsig electron   -4,5
    cuts.AddCutPCMCalo("00081113","0020000969730000dge0400000","111113206f032230000","0163103100000010"); // nsig electron -2.5,4
    cuts.AddCutPCMCalo("00081113","00200009f5730000dge0400000","111113206f032230000","0163103100000010"); // nsig pion 2,-10
    cuts.AddCutPCMCalo("00081113","00200009f1730000dge0400000","111113206f032230000","0163103100000010"); // nsig pion 0,-10
  } else if (trainConfig == 774) {
    cuts.AddCutPCMCalo("00081113","00200009f9030000dge0400000","111113206f032230000","0163103100000010"); // pion nsig min mom 0.50 GeV/c
    cuts.AddCutPCMCalo("00081113","00200009f9630000dge0400000","111113206f032230000","0163103100000010"); // pion nsig min mom 0.25 GeV/c
    cuts.AddCutPCMCalo("00081113","00200009f9760000dge0400000","111113206f032230000","0163103100000010"); // pion nsig max mom 2.00 GeV/c
    cuts.AddCutPCMCalo("00081113","00200009f9710000dge0400000","111113206f032230000","0163103100000010"); // pion nsig max mom 5.00 GeV/c
  } else if (trainConfig == 775) {
    cuts.AddCutPCMCalo("00081113","00200009f9730000age0400000","111113206f032230000","0163103100000010"); // qT max 0.040, qt pt max 0.11
    cuts.AddCutPCMCalo("00081113","00200009f9730000ege0400000","111113206f032230000","0163103100000010"); // qT max 0.060, qt pt max 0.14
    cuts.AddCutPCMCalo("00081113","00200009f9730000fge0400000","111113206f032230000","0163103100000010"); // qT max 0.070, qt pt max 0.16
  } else if (trainConfig == 776) {
    cuts.AddCutPCMCalo("00081113","00200009f9730000d1e0400000","111113206f032230000","0163103100000010"); // chi2 50 no chi2 dep.
    cuts.AddCutPCMCalo("00081113","00200009f9730000dfe0400000","111113206f032230000","0163103100000010"); // chi2 50 chi2 dep -0.065
    cuts.AddCutPCMCalo("00081113","00200009f9730000dhe0400000","111113206f032230000","0163103100000010"); // chi2 50 chi2 dep -0.050
    cuts.AddCutPCMCalo("00081113","00200009f9730000dge0404000","111113206f032230000","0163103100000010"); // reject close v0
    cuts.AddCutPCMCalo("00081113","00200009f9730000dge0406000","111113206f032230000","0163103100000010"); // double count with open angle 0.04
  } else if (trainConfig == 777) {
    cuts.AddCutPCMCalo("00081113","00200009f9730000dgd0400000","111113206f032230000","0163103100000010"); // Psi pair 0.15 dep
    cuts.AddCutPCMCalo("00081113","00200009f9730000dgf0400000","111113206f032230000","0163103100000010"); // Psi pair 0.20 dep
    cuts.AddCutPCMCalo("00081113","00200009f9730000dgg0400000","111113206f032230000","0163103100000010"); // Psi pair 0.30 dep
  } else if (trainConfig == 778) {
    cuts.AddCutPCMCalo("00081113","00200009f9730000dge0400000","111113206f032230000","0263103100000010"); // variation BG scheme track mult
    cuts.AddCutPCMCalo("00081113","00200009f9730000dge0400000","111113206f032230000","0163107100000010"); // alpha meson 0.85
    cuts.AddCutPCMCalo("00081113","00200009f9730000dge0400000","111113206f032230000","0163105100000010"); // alpha meson 0.75
    cuts.AddCutPCMCalo("00081113","00200009227300008250404000","111113206f032230000","0163103100000010"); // old cuts (run1)
  } else if (trainConfig == 779) {
    cuts.AddCutPCMCalo("00081213","00200009f9730000dge0400000","111113206f032230000","0163103100000010"); //same as std + maximum past future rejection
    cuts.AddCutPCMCalo("00081513","00200009f9730000dge0400000","111113206f032230000","0163103100000010"); //same as std + medium past future rejection
  } else if (trainConfig == 780) { // CALO variations
    cuts.AddCutPCMCalo("00081113","00200009f9730000dge0400000","111113206f022230000","0163103100000010"); // 0.6 GeV/c
    cuts.AddCutPCMCalo("00081113","00200009f9730000dge0400000","111113206f042230000","0163103100000010"); // 0.8 GeV/c
  } else if (trainConfig == 781) {
    cuts.AddCutPCMCalo("00081113","00200009f9730000dge0400000","111113206f031230000","0163103100000010"); // n cells >= 1
    cuts.AddCutPCMCalo("00081113","00200009f9730000dge0400000","111113206f033230000","0163103100000010"); // n cells >= 3
    cuts.AddCutPCMCalo("00081113","00200009f9730000dge0400000","111113206f032200000","0163103100000010"); // no max M02 cut
    cuts.AddCutPCMCalo("00081113","00200009f9730000dge0400000","111113206f032250000","0163103100000010"); // M02 < 0.3
    cuts.AddCutPCMCalo("00081113","00200009f9730000dge0400000","111113206f0322k0000","0163103100000010"); // M02, pT-dep
  } else if (trainConfig == 782) {
    cuts.AddCutPCMCalo("00081113","00200009f9730000dge0400000","111113206e032230000","0163103100000010"); // TM var e/p 2.0
    cuts.AddCutPCMCalo("00081113","00200009f9730000dge0400000","111113206g032230000","0163103100000010"); // TM var e/p 1.5
    cuts.AddCutPCMCalo("00081113","00200009f9730000dge0400000","111113206h032230000","0163103100000010"); // TM var e/p 1.25
    cuts.AddCutPCMCalo("00081113","00200009f9730000dge0400000","1111132067032230000","0163103100000010"); // TM var no veto
  } else if (trainConfig == 783) {
    cuts.AddCutPCMCalo("00081113","00200009f9730000dge0400000","111113306f032230000","0163103100000010"); // NL 32
    cuts.AddCutPCMCalo("00081113","00200009f9730000dge0400000","111113406f032230000","0163103100000010"); // NL 33
    cuts.AddCutPCMCalo("00081113","00200009f9730000dge0400000","111113806f032230000","0163103100000010"); // NL 34
    cuts.AddCutPCMCalo("00081113","00200009f9730000dge0400000","111113906f032230000","0163103100000010"); // NL 34
  } else if (trainConfig == 784) {
    cuts.AddCutPCMCalo("00081113","00200009f9730000dge0400000","111113205f032230000","0163103100000010"); // 50ns
    cuts.AddCutPCMCalo("00081113","00200009f9730000dge0400000","111113204f032230000","0163103100000010"); // 100ns
    cuts.AddCutPCMCalo("00081113","00200009f9730000dge0400000","111113207f032230000","0163103100000010"); // 30ns


  } else if (trainConfig == 789){ // same as 790 but with smearing
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","111113106f032230000","0163103100b00010"); // std
    cuts.AddCutPCMCalo("00052113","0dm00009f9730000dge0404000","111113106f032230000","0163103100b00010"); // std
    cuts.AddCutPCMCalo("00081113","0dm00009f9730000dge0404000","111113106f032230000","0163103100b00010"); // std
  } else if (trainConfig == 790){ // EMCAL clusters 8 TeV LHC12, TB+finetuning CCRF
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","111113106f032230000","0163103100000010"); // std
    cuts.AddCutPCMCalo("00052113","0dm00009f9730000dge0404000","111113106f032230000","0163103100000010"); // std
    cuts.AddCutPCMCalo("00081113","0dm00009f9730000dge0404000","111113106f032230000","0163103100000010"); // std
  } else if (trainConfig == 791){ // EMCAL clusters 8 TeV LHC12, TB+finetuning CCRF
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","111113806f032230000","0163103100000010"); // std
    cuts.AddCutPCMCalo("00052113","0dm00009f9730000dge0404000","111113806f032230000","0163103100000010"); // std
    cuts.AddCutPCMCalo("00081113","0dm00009f9730000dge0404000","111113806f032230000","0163103100000010"); // std
  } else if (trainConfig == 792){ // EMCAL clusters 8 TeV LHC12, TB+finetuning CCRF
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","111113906f032230000","0163103100000010"); // std
    cuts.AddCutPCMCalo("00052113","0dm00009f9730000dge0404000","111113906f032230000","0163103100000010"); // std
    cuts.AddCutPCMCalo("00081113","0dm00009f9730000dge0404000","111113906f032230000","0163103100000010"); // std
  } else if (trainConfig == 793){ // EMCAL clusters 8 TeV LHC12, TB+finetuning CCRF
    cuts.AddCutPCMCalo("00010613","0dm00009f9730000dge0404000","111113106f032230000","0163103100000010"); // std
    cuts.AddCutPCMCalo("00052613","0dm00009f9730000dge0404000","111113106f032230000","0163103100000010"); // std
    cuts.AddCutPCMCalo("00081613","0dm00009f9730000dge0404000","111113106f032230000","0163103100000010"); // std
  } else if (trainConfig == 794){ // EMCAL clusters 8 TeV LHC12, TB+finetuning CCRF
    cuts.AddCutPCMCalo("00010713","0dm00009f9730000dge0404000","111113106f032230000","0163103100000010"); // std
    cuts.AddCutPCMCalo("00052713","0dm00009f9730000dge0404000","111113106f032230000","0163103100000010"); // std
    cuts.AddCutPCMCalo("00081713","0dm00009f9730000dge0404000","111113106f032230000","0163103100000010"); // std
  } else if (trainConfig == 795){ // TOF single leg requirement
    cuts.AddCutPCMCalo("00010113","0dm00009f9730600dge0404000","111113106f032230000","0163103100000010"); // std
    cuts.AddCutPCMCalo("00052113","0dm00009f9730600dge0404000","111113106f032230000","0163103100000010"); // std
    cuts.AddCutPCMCalo("00081113","0dm00009f9730600dge0404000","111113106f032230000","0163103100000010"); // std
  } else if (trainConfig == 796){ // TOF both leg requirement
    cuts.AddCutPCMCalo("00010113","0dm00009f9730700dge0404000","111113106f032230000","0163103100000010"); // std
    cuts.AddCutPCMCalo("00052113","0dm00009f9730700dge0404000","111113106f032230000","0163103100000010"); // std
    cuts.AddCutPCMCalo("00081113","0dm00009f9730700dge0404000","111113106f032230000","0163103100000010"); // std
  } else if (trainConfig == 797){ // TOF single leg requirement
    cuts.AddCutPCMCalo("00010113","0dm00009f9730800dge0404000","111113106f032230000","0163103100000010"); // std
    cuts.AddCutPCMCalo("00052113","0dm00009f9730800dge0404000","111113106f032230000","0163103100000010"); // std
    cuts.AddCutPCMCalo("00081113","0dm00009f9730800dge0404000","111113106f032230000","0163103100000010"); // std
  } else if (trainConfig == 798){ // TOF both leg requirement
    cuts.AddCutPCMCalo("00010113","0dm00009f9730900dge0404000","111113106f032230000","0163103100000010"); // std
    cuts.AddCutPCMCalo("00052113","0dm00009f9730900dge0404000","111113106f032230000","0163103100000010"); // std
    cuts.AddCutPCMCalo("00081113","0dm00009f9730900dge0404000","111113106f032230000","0163103100000010"); // std
  } else if (trainConfig == 799){ // no NCell cut
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","111113106f030230000","0163103100000010"); // std
    cuts.AddCutPCMCalo("00052113","0dm00009f9730000dge0404000","111113106f030230000","0163103100000010"); // std
    cuts.AddCutPCMCalo("00081113","0dm00009f9730000dge0404000","111113106f030230000","0163103100000010"); // std

  //*************************************************************************************************
  // 5 TeV PHOS - setup
  //*************************************************************************************************
  } else if ( trainConfig == 800){ // QA
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","2446600040012300000","0163103100000010"); // INT7
    cuts.AddCutPCMCalo("00062113","00200009327000008250400000","2446600040012300000","0163103100000010"); // PHI7
  } else if ( trainConfig == 801){ // Default cut, No TM
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","2446651040012300000","0163103100000010"); // INT7
    cuts.AddCutPCMCalo("00062113","00200009327000008250400000","2446651040012300000","0163103100000010"); // PHI7
  } else if ( trainConfig == 802){ // Default cut, with TM
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","2446651044012300000","0163103100000010"); // INT7
    cuts.AddCutPCMCalo("00062113","00200009327000008250400000","2446651044012300000","0163103100000010"); // PHI7
  } else if ( trainConfig == 803){ // NL variations
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","2446600044012300000","0163103100000010"); // No NL
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","2446601044012300000","0163103100000010"); // PHOS people NL
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","2446651044012300000","0163103100000010"); // INT7
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","2446652044012300000","0163103100000010"); // PHOS calo NL
  } else if ( trainConfig == 804){ //dEdx e-line variation and dE/dx pi-line variation
    cuts.AddCutPCMCalo("00010113","00200009127000008250400000","2446651044012300000","0163103100000010"); //-5 < sigma < 5
    cuts.AddCutPCMCalo("00010113","00200009227000008250400000","2446651044012300000","0163103100000010"); //-3 < sigma < 5
    cuts.AddCutPCMCalo("00010113","00200009327400008250400000","2446651044012300000","0163103100000010"); //1, -10, 0.4, 3
    cuts.AddCutPCMCalo("00010113","00200009367400008250400000","2446651044012300000","0163103100000010"); //2, 0.5, 0.4, 3
  } else if ( trainConfig == 805){ //dE/dx pi-line variation and single pT variation and 2D triangular chi2 and psi pair
    cuts.AddCutPCMCalo("00010113","00200009317400008250400000","2446651044012300000","0163103100000010"); //0, -10, 0.4, 3
    cuts.AddCutPCMCalo("00010113","00200049327000008250400000","2446651044012300000","0163103100000010"); //single pT 0.075 GeV/c
    cuts.AddCutPCMCalo("00010113","00200019327000008250400000","2446651044012300000","0163103100000010"); //single pT 0.1   GeV/c
    cuts.AddCutPCMCalo("00010113","00200009327000008850400000","2446651044012300000","0163103100000010"); //20 & 0.1
  } else if ( trainConfig == 806){ //2D triangular chi2 and psi pair
    cuts.AddCutPCMCalo("00010113","00200009327000008260400000","2446651044012300000","0163103100000010"); //30 & 0.05
    cuts.AddCutPCMCalo("00010113","00200009327000008860400000","2446651044012300000","0163103100000010"); //20 & 0.05
    cuts.AddCutPCMCalo("00010113","00200009327000008280400000","2446651044012300000","0163103100000010"); //30 & 0.2
    cuts.AddCutPCMCalo("00010113","00200009327000008880400000","2446651044012300000","0163103100000010"); //20 & 0.2
  } else if ( trainConfig == 807){ //min TPC clusters and //elliptic cut Qt and alpha
    cuts.AddCutPCMCalo("00010113","00200000327000008250400000","2446651044012300000","0163103100000010"); //0
    cuts.AddCutPCMCalo("00010113","00200008327000008250400000","2446651044012300000","0163103100000010"); //0.35
    cuts.AddCutPCMCalo("00010113","00200009327000009250400000","2446651044012300000","0163103100000010"); // qT 0.03 no quadratic
    cuts.AddCutPCMCalo("00010113","00200009327000003250400000","2446651044012300000","0163103100000010"); // qT 0.05 y  quadratic
  } else if ( trainConfig == 808){ // first set of variations CLUSTER
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","2446651044072300000","0163103100000010"); // min energy 0.2 GeV
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","2446651044082300000","0163103100000010"); // min energy 0.4 GeV
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","2446651044022300000","0163103100000010"); // min energy 0.5 GeV
  } else if ( trainConfig == 809){ // second set of variations CLUSTER
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","2446651044012370000","0163103100000010"); // min/max M02  0.2<M<1
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","2446651044012380000","0163103100000010"); // min/max M02  0.2<M<1
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","2446651044013300000","0163103100000010"); // min number 3 cells
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","2446651044014300000","0163103100000010"); // min number 4 cells
  } else if ( trainConfig == 810){ // MESON variations
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","2446651044012300000","0163403100000010"); // rapidity variation  y<0.5
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","2446651044012300000","0163803100000010"); // rapidity variation  y<0.25
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","2446651044012300000","0163106100000010"); // alpha meson variation 1 0<alpha<0.8
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","2446651044012300000","0163105100000010"); // alpha meson variation 2 0<alpha<0.75
  } else if ( trainConfig == 811){ // TM & opening angle
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","2446651045012300000","0163103100000010"); // tm variation
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","2446651046012300000","0163103100000010"); // tm variation
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","2446651044012300000","0163103100000000"); // min opening angle 0    -> open
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","2446651044012300000","0163103100000030"); // min opening angle 0.01 -> 2 cell diag
  } else if ( trainConfig == 812){ // NL variations
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","2446600044012300000","0163103100000010"); // No NL
  } else if ( trainConfig == 813){ // new default 2019 july 18
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0400000","244665107a012200000","0h63103100000010"); // No NL
  } else if ( trainConfig == 814){ // No non-lin corr, use with Run2Tune / Run2TuneMC
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0400000","244660007a012200000","0h63103100000010"); // No NL
  } else if ( trainConfig == 815){ // Sphericity PCMPHOS
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0400000","24466510ga012200000","0h63103100000010"); // No NL
    cuts.AddCutPCMCalo("h0510113","00200009f9730000dge0400000","24466510ga012200000","0h63103100000010"); // No NL
    cuts.AddCutPCMCalo("h5a10113","00200009f9730000dge0400000","24466510ga012200000","0h63103100000010"); // No NL
    cuts.AddCutPCMCalo("h0a10113","00200009f9730000dge0400000","24466510ga012200000","0h63103100000010"); // No NL
    cuts.AddCutPCMCalo("h0310113","00200009f9730000dge0400000","24466510ga012200000","0h63103100000010"); // No NL
    cuts.AddCutPCMCalo("h7a10113","00200009f9730000dge0400000","24466510ga012200000","0h63103100000010"); // No NL
  } else if ( trainConfig == 816){ // V0M mult selection PCMPHOS
    cuts.AddCutPCMCalo("m0110113","00200009f9730000dge0400000","24466510ga012200000","0h63103100000010"); // 0-1%
    cuts.AddCutPCMCalo("m1510113","00200009f9730000dge0400000","24466510ga012200000","0h63103100000010"); // 1-5%
    cuts.AddCutPCMCalo("m5k10113","00200009f9730000dge0400000","24466510ga012200000","0h63103100000010"); // 5-20%
    cuts.AddCutPCMCalo("n2410113","00200009f9730000dge0400000","24466510ga012200000","0h63103100000010"); // 20-40%
    cuts.AddCutPCMCalo("n4710113","00200009f9730000dge0400000","24466510ga012200000","0h63103100000010"); // 40-70%
    cuts.AddCutPCMCalo("n7a10113","00200009f9730000dge0400000","24466510ga012200000","0h63103100000010"); // 70-100%
  } else if ( trainConfig == 817){ // SPD mult selection PCMPHOS
    cuts.AddCutPCMCalo("o0110113","00200009f9730000dge0400000","24466510ga012200000","0h63103100000010"); // 0-1%
    cuts.AddCutPCMCalo("o0210113","00200009f9730000dge0400000","24466510ga012200000","0h63103100000010"); // 0-2%
    cuts.AddCutPCMCalo("o0510113","00200009f9730000dge0400000","24466510ga012200000","0h63103100000010"); // 0-5%
    cuts.AddCutPCMCalo("o5k10113","00200009f9730000dge0400000","24466510ga012200000","0h63103100000010"); // 5-20%
    cuts.AddCutPCMCalo("p2610113","00200009f9730000dge0400000","24466510ga012200000","0h63103100000010"); // 20-60%
    cuts.AddCutPCMCalo("p6a10113","00200009f9730000dge0400000","24466510ga012200000","0h63103100000010"); // 60-100%
  } else if ( trainConfig == 818){ // new default 2019 august 14
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0400000","24466510ga012200000","0h63103100000010"); // No NL

    //Normal B Option
  } else if ( trainConfig == 820){ // Default cut, with TM   with eta<0.8
    cuts.AddCutPCMCalo("00010113","0d200009327000008250400000","2446651044012300000","0163103100000010"); // INT7
  } else if ( trainConfig == 821){ // Default cut, with TM   with eta<0.8
    cuts.AddCutPCMCalo("00010113","0da00009327000008250400000","2446651044012300000","0163103100000010"); // INT7   RBins
    cuts.AddCutPCMCalo("00010113","0db00009327000008250400000","2446651044012300000","0163103100000010"); // INT7
    cuts.AddCutPCMCalo("00010113","0dc00009327000008250400000","2446651044012300000","0163103100000010"); // INT7
  } else if ( trainConfig == 822){ // Default cut, with TM   with eta<0.8
    cuts.AddCutPCMCalo("00010113","0dh00009327000008250400000","2446651044012300000","0163103100000010"); // INT7   RBins
    cuts.AddCutPCMCalo("00010113","0di00009327000008250400000","2446651044012300000","0163103100000010"); // INT7

      //Low B Option
  } else if ( trainConfig == 823){ // Default cut, with TM   with eta<0.8
    cuts.AddCutPCMCalo("00010113","0d200089327000008250400000","2446651044012300000","0163103100000010"); // INT7
  } else if ( trainConfig == 824){ // Default cut, with TM   with eta<0.8
    cuts.AddCutPCMCalo("00010113","0da00089327000008250400000","2446651044012300000","0163103100000010"); // INT7   RBins
    cuts.AddCutPCMCalo("00010113","0db00089327000008250400000","2446651044012300000","0163103100000010"); // INT7
    cuts.AddCutPCMCalo("00010113","0dc00089327000008250400000","2446651044012300000","0163103100000010"); // INT7
  } else if ( trainConfig == 825){ // Default cut, with TM   with eta<0.8
    cuts.AddCutPCMCalo("00010113","0dh00089327000008250400000","2446651044012300000","0163103100000010"); // INT7   RBins
    cuts.AddCutPCMCalo("00010113","0di00089327000008250400000","2446651044012300000","0163103100000010"); // INT7

      //Normal B Option
  } else if ( trainConfig == 830){ // Default cut, with TM   with eta<0.8    //To be used with MBW
    cuts.AddCutPCMCalo("00010113","0d200009327000008250400000","2446651044012300000","0163103100000010"); // INT7
  } else if ( trainConfig == 831){ // Default cut, with TM   with eta<0.8
    cuts.AddCutPCMCalo("00010113","0da00009327000008250400000","2446651044012300000","0163103100000010"); // INT7   RBins
    cuts.AddCutPCMCalo("00010113","0db00009327000008250400000","2446651044012300000","0163103100000010"); // INT7
    cuts.AddCutPCMCalo("00010113","0dc00009327000008250400000","2446651044012300000","0163103100000010"); // INT7
  } else if ( trainConfig == 832){ // Default cut, with TM   with eta<0.8
    cuts.AddCutPCMCalo("00010113","0dh00009327000008250400000","2446651044012300000","0163103100000010"); // INT7   RBins
    cuts.AddCutPCMCalo("00010113","0di00009327000008250400000","2446651044012300000","0163103100000010"); // INT7

      //Low B Option
  } else if ( trainConfig == 833){ // Default cut, with TM   with eta<0.8    //To be used with MBW
    cuts.AddCutPCMCalo("00010113","0d200089327000008250400000","2446651044012300000","0163103100000010"); // INT7
  } else if ( trainConfig == 834){ // Default cut, with TM   with eta<0.8
    cuts.AddCutPCMCalo("00010113","0da00089327000008250400000","2446651044012300000","0163103100000010"); // INT7   RBins
    cuts.AddCutPCMCalo("00010113","0db00089327000008250400000","2446651044012300000","0163103100000010"); // INT7
    cuts.AddCutPCMCalo("00010113","0dc00089327000008250400000","2446651044012300000","0163103100000010"); // INT7
  } else if ( trainConfig == 835){ // Default cut, with TM   with eta<0.8
    cuts.AddCutPCMCalo("00010113","0dh00089327000008250400000","2446651044012300000","0163103100000010"); // INT7   RBins
    cuts.AddCutPCMCalo("00010113","0di00089327000008250400000","2446651044012300000","0163103100000010"); // INT7

  //PCM-PHOS pp HBT studies <- reserved 850 to 860
  } else if ( trainConfig == 850){ // Default cut, with TM
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0404000","24466510ga012200000","0h63103100000010"); // INT7

  //*************************************************************************************************
  // 13 TeV PHOS - setup
  //*************************************************************************************************
  } else if ( trainConfig == 900){ // INT7, 300MeV
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163103100000010"); // INT7 no Trigger
  } else if ( trainConfig == 901){ // PHI7, 300MeV
    cuts.AddCutPCMCalo("00062113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163103100000010"); //PHI7
  } else if ( trainConfig == 902){ // INT7 Track Matching Off
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","24466190s001cc00000","0163103100000010"); // INT7 no Trigger
  } else if ( trainConfig == 903){ // PHI7  Track Matching Off
    cuts.AddCutPCMCalo("00062113","0dm00009f9730000dge0404000","24466190s001cc00000","0163103100000010"); //PHI7
  } else if ( trainConfig == 904){ //PHI7 Gamma Cut
    cuts.AddCutPCMCalo("00062113","0dm00009f9730000dge0404000","24466190sa01cc00000","01631g3100000010"); //PHI7  Gamma Energy cut 6 GeV on one Cluster
    cuts.AddCutPCMCalo("00062113","0dm00009f9730000dge0404000","24466190sa01cc00000","01631j3100000010"); //PHI7  Gamma Energy cut 2 GeV on one Cluster
    cuts.AddCutPCMCalo("00062113","0dm00009f9730000dge0404000","24466190sa01cc00000","01631k3100000010"); //PHI7  Gamma Energy cut 4 GeV on one Cluster
  } else if ( trainConfig == 905){ // NCell Cuts Variations
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa012200000","0163103100000010"); // INT7
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa012c00000","0163103100000010"); // INT7
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa01c200000","0163103100000010"); // INT7
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163103100000010"); // INT7
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa01d200000","0163103100000010"); // INT7
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa01dc00000","0163103100000010"); // INT7
  } else if ( trainConfig == 906){ // StdCut // reconstructed conversion
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa01cc00300","0163103100000010"); //no Trigger
  } else if ( trainConfig == 907){ // NonLinStudies
    //cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","24466000sa01cc00000","0163103100000010"); // QA
    //cuts.AddCutPCMCalo("00062113","0dm00009f9730000dge0404000","24466000sa01cc00000","0163103100000010"); // QA
    //-
    //cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","24466510sa01cc00000","0163103100000010"); // 51 Nonlin MB
    //cuts.AddCutPCMCalo("00062113","0dm00009f9730000dge0404000","24466510sa01cc00000","0163103100000010"); // 51 Nonlin Triggerd
    //-
    //cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","24466110sa01cc00000","0163103100000010"); // QA //case 11=> FunctionNL_kSDM MB PCMPHOS
    //cuts.AddCutPCMCalo("00062113","0dm00009f9730000dge0404000","24466110sa01cc00000","0163103100000010"); // QA //case 11=> FunctionNL_kSDM Triggerd PCMPHOS
    //-
    //cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","24466120sa01cc00000","0163103100000010"); // QA //case 12=> FunctionNL_kSDM MB PHOSPHOS
    //cuts.AddCutPCMCalo("00062113","0dm00009f9730000dge0404000","24466120sa01cc00000","0163103100000010"); // QA //case 12=> FunctionNL_kSDM Triggerd PHOSPHOS
    //-
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163103100000010"); // QA //case 19=> FunctionNL_kSDM MB PCMPHOS, Tuned with PHOSPHOS
    cuts.AddCutPCMCalo("00062113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163103100000010"); // QA //case 19=> FunctionNL_kSDM Triggerd PCMPHOS, Tuned with PHOSPHOS
    //-
    //cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","24466210sa01cc00000","0163103100000010"); // QA //case 21=> unctionNL_DPOW MB
    //cuts.AddCutPCMCalo("00062113","0dm00009f9730000dge0404000","24466210sa01cc00000","0163103100000010"); // QA //case 21=> unctionNL_DPOW Triggerd
  } else if ( trainConfig == 908){  //PHOS PHI7 // reconstructed conversion
    cuts.AddCutPCMCalo("00062113","0dm00009f9730000dge0404000","24466190sa01cc00300","0163103100000010"); //PHI7
  } else if ( trainConfig == 909){  //PHOS Smearing Check
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163103100b00010"); //b
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163103100o00010"); //o  new approach
  } else if ( trainConfig == 910){  //PHOS Sphericiry Check
    cuts.AddCutPCMCalo("h0510113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163103100000010"); //  0.    - 0.5
    cuts.AddCutPCMCalo("h5a10113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163103100000010"); //  0.5    - 1.
    cuts.AddCutPCMCalo("h0a10113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163103100000010"); //  0.    - 1.
  } else if ( trainConfig == 911){  //PHOS Mult Check
    cuts.AddCutPCMCalo("n0110113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163103100000010"); // INT7 0-10%
    cuts.AddCutPCMCalo("n1210113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163103100000010"); // INT7 10-20%
    cuts.AddCutPCMCalo("n2510113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163103100000010"); // INT7 20-50%
    cuts.AddCutPCMCalo("n5a10113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163103100000010"); // INT7 50-100%
    cuts.AddCutPCMCalo("m0110113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163103100000010"); // INT7
    cuts.AddCutPCMCalo("m1510113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163103100000010"); // INT7
    cuts.AddCutPCMCalo("m5a10113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163103100000010"); // INT7
  } else if ( trainConfig == 912){  //PHOS MB and PHI7 Timing Cut 0; 100MeV
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","244661900a09cc00000","0163103100000010"); //no Trigger
    cuts.AddCutPCMCalo("00062113","0dm00009f9730000dge0404000","244661900a09cc00000","0163103100000010"); //PHI7
  } else if ( trainConfig == 913){ //PHOS MB and PHI7 Timing Cut 0; 300MeV
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","244661900a01cc00000","0163103100000010"); //no Trigger
    cuts.AddCutPCMCalo("00062113","0dm00009f9730000dge0404000","244661900a01cc00000","0163103100000010"); //PHI7
  } else if ( trainConfig == 914){  //PHOS Triggers Timing Cut Studies
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","244661907a01cc00000","0163103100000010"); //no Trigger, Mike's Timing
    cuts.AddCutPCMCalo("00062113","0dm00009f9730000dge0404000","244661907a01cc00000","0163103100000010"); //PHI7, Mike's Timing
  } else if ( trainConfig == 915){ // PHI7, 300MeV,  Maybe Bad DDLs for Trigger thrown out
    cuts.AddCutPCMCalo("00062113","0dm00009f9730000dge0404000","2446b190sa01cc00000","0163103100000010"); //PHI7,  Maybe Bad DDLs for Trigger thrown out
  } else if ( trainConfig == 916){ // INT7, 300MeV; Trigger Helper only uses events with (INT7 & PHI7) flag
    cuts.AddCutPCMCalo("00014113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163103100000010"); // INT7 no Trigger
  } else if ( trainConfig == 917){ // TimingEff; 2GeV<ETag<5.5GeV, |TimingTag|<30ns, |TimingProbe|<1000ns, SignalExtraction, LowPt from Trigger, HighPt const
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","24466190va01cc00000","0163103100000010"); //INT7
    cuts.AddCutPCMCalo("00062113","0dm00009f9730000dge0404000","24466190va01cc00000","0163103100000010"); //PHI7
  } else if ( trainConfig == 918){ // TimingEff; 2GeV<ETag<5.5GeV, |TimingTag|<30ns, |TimingProbe|<1000ns, SignalExtraction, LowPt from MB
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","24466190wa01cc00000","0163103100000010"); //INT7
    cuts.AddCutPCMCalo("00062113","0dm00009f9730000dge0404000","24466190wa01cc00000","0163103100000010"); //PHI7
  } else if ( trainConfig == 919){ // TimingEff; 2GeV<ETag<5.5GeV, |TimingTag|<30ns, |TimingProbe|<1000ns, SignalExtraction, LowPt from Trigger
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","24466190xa01cc00000","0163103100000010"); //INT7
    cuts.AddCutPCMCalo("00062113","0dm00009f9730000dge0404000","24466190xa01cc00000","0163103100000010"); //PHI7
    //Normal B Option
 } else if ( trainConfig == 920){ // Default cut, with TM   with eta<0.8
    cuts.AddCutPCMCalo("00010113","0d200009f9730000dge0404000","24466190sa01cc00000","0163103100000010"); // INT7
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163103100000010"); // INT7
 } else if ( trainConfig == 921){ // Default cut, with TM   with eta<0.8
    cuts.AddCutPCMCalo("00010113","0da00009f9730000dge0404000","24466190sa01cc00000","0163103100000010"); // INT7   RBins
    cuts.AddCutPCMCalo("00010113","0db00009f9730000dge0404000","24466190sa01cc00000","0163103100000010"); // INT7
    cuts.AddCutPCMCalo("00010113","0dc00009f9730000dge0404000","24466190sa01cc00000","0163103100000010"); // INT7
 } else if ( trainConfig == 922){ // Default cut, with TM   with eta<0.8
    cuts.AddCutPCMCalo("00010113","0dh00009f9730000dge0404000","24466190sa01cc00000","0163103100000010"); // INT7   RBins
    cuts.AddCutPCMCalo("00010113","0di00009f9730000dge0404000","24466190sa01cc00000","0163103100000010"); // INT7
    cuts.AddCutPCMCalo("00010113","0dj00009f9730000dge0404000","24466190sa01cc00000","0163103100000010"); // INT7
    cuts.AddCutPCMCalo("00010113","0dk00009f9730000dge0404000","24466190sa01cc00000","0163103100000010"); // INT7
    cuts.AddCutPCMCalo("00010113","0dl00009f9730000dge0404000","24466190sa01cc00000","0163103100000010"); // INT7
    cuts.AddCutPCMCalo("00010113","0dg00009f9730000dge0404000","24466190sa01cc00000","0163103100000010"); // INT7

    //Mimic var Option
  } else if ( trainConfig == 923){ // INT7, 300MeV mimic var
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163103100000010"); // INT7 no Trigger
  } else if ( trainConfig == 924){ // INT7, 300MeV mimic var
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163103100000010"); // INT7 no Trigger
  } else if ( trainConfig == 925){ // PHI7, 300MeV mimic var
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163103100000010"); // INT7 no Trigger
  } else if ( trainConfig == 926){ // INT7, 300MeV mimic var
    cuts.AddCutPCMCalo("00062113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163103100000010"); //PHI7
  } else if ( trainConfig == 927){ // INT7, 300MeV mimic var
    cuts.AddCutPCMCalo("00062113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163103100000010"); //PHI7
  } else if ( trainConfig == 928){ // PHI7, 300MeV mimic var
    cuts.AddCutPCMCalo("00062113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163103100000010"); //PHI7
  } else if ( trainConfig == 929){ // PHI7, 300MeV mimic var
    cuts.AddCutPCMCalo("00062113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163103100000010"); //PHI7

    //Normal B Option
 } else if ( trainConfig == 930){ // Default cut, with TM   with eta<0.8    //To be used with MBW
    cuts.AddCutPCMCalo("00010113","0d200009f9730000dge0404000","24466190sa01cc00000","0163103100000010"); // INT7
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163103100000010"); // INT7
 } else if ( trainConfig == 931){ // Default cut, with TM   with eta<0.8
    cuts.AddCutPCMCalo("00010113","0da00009f9730000dge0404000","24466190sa01cc00000","0163103100000010"); // INT7   RBins
    cuts.AddCutPCMCalo("00010113","0db00009f9730000dge0404000","24466190sa01cc00000","0163103100000010"); // INT7
    cuts.AddCutPCMCalo("00010113","0dc00009f9730000dge0404000","24466190sa01cc00000","0163103100000010"); // INT7
    cuts.AddCutPCMCalo("00010113","0di00009f9730000dge0404000","24466190sa01cc00000","0163103100000010"); // INT7
 } else if ( trainConfig == 932){ // Default cut, with TM   with eta<0.8
    cuts.AddCutPCMCalo("00010113","0dh00009f9730000dge0404000","24466190sa01cc00000","0163103100000010"); // INT7   RBins
    cuts.AddCutPCMCalo("00010113","0dj00009f9730000dge0404000","24466190sa01cc00000","0163103100000010"); // INT7
    cuts.AddCutPCMCalo("00010113","0dk00009f9730000dge0404000","24466190sa01cc00000","0163103100000010"); // INT7
    cuts.AddCutPCMCalo("00010113","0dl00009f9730000dge0404000","24466190sa01cc00000","0163103100000010"); // INT7
    cuts.AddCutPCMCalo("00010113","0dg00009f9730000dge0404000","24466190sa01cc00000","0163103100000010"); // INT7

    //Low B Option
 } else if ( trainConfig == 933){ // Default cut, with TM   with eta<0.8    //To be used with MBW
    cuts.AddCutPCMCalo("00010113","0d200089327000008250404000","24466190sa01cc00000","0163103100000010"); // INT7
 } else if ( trainConfig == 934){ // Default cut, with TM   with eta<0.8
    cuts.AddCutPCMCalo("00010113","0da00089327000008250404000","24466190sa01cc00000","0163103100000010"); // INT7   RBins
    cuts.AddCutPCMCalo("00010113","0db00089327000008250404000","24466190sa01cc00000","0163103100000010"); // INT7
    cuts.AddCutPCMCalo("00010113","0dc00089327000008250404000","24466190sa01cc00000","0163103100000010"); // INT7
 } else if ( trainConfig == 935){ // Default cut, with TM   with eta<0.8
    cuts.AddCutPCMCalo("00010113","0dh00089327000008250404000","24466190sa01cc00000","0163103100000010"); // INT7   RBins
    cuts.AddCutPCMCalo("00010113","0di00089327000008250404000","24466190sa01cc00000","0163103100000010"); // INT7

  } else if ( trainConfig == 940){ // INT7, 100MeV, with Timing Efficiency
      cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","24466000sa09cc00000","0163103100000010"); // INT7 no Trigger
  } else if ( trainConfig == 941){ // INT7, 100MeV, noTimingEff
      cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","244660000a09cc00000","0163103100000010"); // INT7 no Trigger
  } else if ( trainConfig == 942){ // INT7, 100MeV, with Timing Efficiency, NCells 2
      cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","24466000sa092200000","0163103100000010"); // INT7 no Trigger
  } else if ( trainConfig == 943){ // INT7, 100MeV, noTimingEff, NCells 2
      cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","244660000a092200000","0163103100000010"); // INT7 no Trigger
  } else if ( trainConfig == 944){ // INT7, 100MeV, with Timing Efficiency, NCells 3
      cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","24466000sa093200000","0163103100000010"); // INT7 no Trigger
  } else if ( trainConfig == 945){ // INT7, 100MeV, noTimingEff, NCells 3
      cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","244660000a093200000","0163103100000010"); // INT7 no Trigger

  } else if ( trainConfig == 946){ // INT7, 300MeV, with Timing Efficiency, NCells 2
      cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa092200000","0163103100000010"); // INT7 no Trigger
  } else if ( trainConfig == 947){ // INT7, 300MeV, noTimingEff, NCells 2
      cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","244661900a092200000","0163103100000010"); // INT7 no Trigger
  } else if ( trainConfig == 948){ // INT7, 300MeV, with Timing Efficiency, NCells 3
      cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa093200000","0163103100000010"); // INT7 no Trigger
  } else if ( trainConfig == 949){ // INT7, 300MeV, noTimingEff, NCells 3
      cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","244661900a093200000","0163103100000010"); // INT7 no Trigger

  //*************************************************************************************************
  // 13 TeV PHOS low B - setup
  //*************************************************************************************************
  } else if ( trainConfig == 950){ // INT7
    cuts.AddCutPCMCalo("00010113","00200089327000008250400000","2446600000013300000","0163103100000010"); // QA
    cuts.AddCutPCMCalo("00062113","00200089327000008250400000","2446600000013300000","0163103100000010"); // QA
  } else if ( trainConfig == 951){ // NCells >= 2
    cuts.AddCutPCMCalo("00010113","00200089327000008250400000","2446600060012300000","0163103100000010"); // QA, -30,50ns timing INT7
    cuts.AddCutPCMCalo("00062113","00200089327000008250400000","2446600060012300000","0163103100000010"); // QA, -30,50ns timing PHI7
  } else if ( trainConfig == 952){ // NCells >= 2
    cuts.AddCutPCMCalo("00010113","00200089327000008250400000","2446611060012300000","0163103100000010"); // QA, -30,50ns timing INT7 w/NL
    cuts.AddCutPCMCalo("00010113","00200089327000008250400000","2446621060012300000","0163103100000010"); // QA, -30,50ns timing INT7 w/NL
  } else if ( trainConfig == 953){ // NCells >= 2
    cuts.AddCutPCMCalo("00010113","00200089327000008250400000","2446600044012300000","0163103100000010"); // QA, -30,50ns
    cuts.AddCutPCMCalo("00010113","00200089327000008250400000","2446600064012300000","0163103100000010"); // QA, -30,50ns
  } else if ( trainConfig == 954){ // NCells >= 2
    cuts.AddCutPCMCalo("00010113","00200089297000001280000000","2446600067012300000","0163103100000010"); // QA, -30,50ns
    cuts.AddCutPCMCalo("00010113","00200089327000008250400000","2446600067012300000","0163103100000010"); // QA, -30,50ns
    cuts.AddCutPCMCalo("00010113","0020008932700000iih0400000","2446600067012300000","0163103100000010"); // QA, -30,50ns
    cuts.AddCutPCMCalo("00010113","0020008932700000i280400000","2446600067012300000","0163103100000010"); // QA, -30,50ns
    cuts.AddCutPCMCalo("00010113","00200089327000001ih0400000","2446600067012300000","0163103100000010"); // QA, -30,50ns

  //*************************************************************************************************
  // 13 TeV PHOS HM trigger
  //*************************************************************************************************
  } else if( trainConfig == 970){ // PHOS HM trigger
    cuts.AddCutPCMCalo("00010113","00200089327000008250400000","2446600044012300000","0163103100000010"); // -30ns, 35ns timing cut, MB trigger
    cuts.AddCutPCMCalo("00010113","00200089327000008250400000","2446600004012300000","0163103100000010"); // no timing, MB trigger
    cuts.AddCutPCMCalo("00074113","00200089327000008250400000","2446600044012300000","0163103100000010"); // -30ns, 35ns timing cut, no NL VOHM
    cuts.AddCutPCMCalo("00076113","00200089327000008250400000","2446600044012300000","0163103100000010"); // -30ns, 35ns timing cut, no NL VOHM with SPD

  //*************************************************************************************************
  // 13 TeV EDC (EMCal + DCal)
  //*************************************************************************************************
  } else if ( trainConfig == 2000){ // EMCAL clusters
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe30220000","0163103100b00010"); // STD old was 4117901067032230000
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe30220000","0163103100000010"); // no smearing // eta
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411790000fe30220000","0163103100000010"); // no timing cut, no NL INT7, no smearing   // for QA
  } else if ( trainConfig == 2001){ // EMCAL clusters  Pi0
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe30220000","0163103100b00010"); // STD
  } else if ( trainConfig == 2002){ // EMCAL clusters  Eta
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe30220000","0163103100000010"); // STD // no smearing // eta
  } else if ( trainConfig == 2003){ // EMCAL clusters       wo trigger mimicing
    cuts.AddCutPCMCalo("0008e113","0dm00009f9730000dge0404000","411790109fe30220000","0163103100b00010"); // EG2+DG2
    cuts.AddCutPCMCalo("0008e113","0dm00009f9730000dge0404000","411790109fe30220000","0163103100000010"); // EG2+DG2
    cuts.AddCutPCMCalo("0008d113","0dm00009f9730000dge0404000","411790109fe30220000","0163103100b00010"); // EG1+DG1
    cuts.AddCutPCMCalo("0008d113","0dm00009f9730000dge0404000","411790109fe30220000","0163103100000010"); // EG1+DG1
  } else if ( trainConfig == 2004){ // EMCAL clusters          with trigger mimicing
    cuts.AddCutPCMCalo("0008e113","0dm00009f9730000dge0404000","411790109fe30220000","0163103100b00010"); // EG2+DG2
    cuts.AddCutPCMCalo("0008e113","0dm00009f9730000dge0404000","411790109fe30220000","0163103100000010"); // EG2+DG2
    cuts.AddCutPCMCalo("0008d113","0dm00009f9730000dge0404000","411790109fe30220000","0163103100b00010"); // EG1+DG1
    cuts.AddCutPCMCalo("0008d113","0dm00009f9730000dge0404000","411790109fe30220000","0163103100000010"); // EG1+DG1
  } else if (trainConfig == 2005){  // EMCal+DCAL clusters standard cuts, Sphericity
    cuts.AddCutPCMCalo("h0310113","0dm00009f9730000dge0404000","411790109fe30220000","0163103100b00010"); //  0.    - 0.3
    cuts.AddCutPCMCalo("h0510113","0dm00009f9730000dge0404000","411790109fe30220000","0163103100b00010"); //  0.    - 0.5
    cuts.AddCutPCMCalo("h3710113","0dm00009f9730000dge0404000","411790109fe30220000","0163103100b00010"); //  0.3    - 0.7
    cuts.AddCutPCMCalo("h5a10113","0dm00009f9730000dge0404000","411790109fe30220000","0163103100b00010"); //  0.5   - 1.
    cuts.AddCutPCMCalo("h7a10113","0dm00009f9730000dge0404000","411790109fe30220000","0163103100b00010"); //  0.7   - 1.
    cuts.AddCutPCMCalo("h0a10113","0dm00009f9730000dge0404000","411790109fe30220000","0163103100b00010"); //  0.    - 1.
  } else if (trainConfig == 2006){  // EMCal+DCAL clusters standard cuts, V0M mult selections
    cuts.AddCutPCMCalo("n0a10113","0dm00009f9730000dge0404000","411790109fe30220000","0163103100b00010"); // INT7 0-100%
    cuts.AddCutPCMCalo("m0110113","0dm00009f9730000dge0404000","411790109fe30220000","0163103100b00010"); // INT7 0-1%
    cuts.AddCutPCMCalo("m1510113","0dm00009f9730000dge0404000","411790109fe30220000","0163103100b00010"); // INT7 1-5%
    cuts.AddCutPCMCalo("m5a10113","0dm00009f9730000dge0404000","411790109fe30220000","0163103100b00010"); // INT7 5-10%
    cuts.AddCutPCMCalo("maf10113","0dm00009f9730000dge0404000","411790109fe30220000","0163103100b00010"); // INT7 10-15%
    cuts.AddCutPCMCalo("mfk10113","0dm00009f9730000dge0404000","411790109fe30220000","0163103100b00010"); // INT7 15-20%
    cuts.AddCutPCMCalo("n2310113","0dm00009f9730000dge0404000","411790109fe30220000","0163103100b00010"); // INT7 20-30%
    cuts.AddCutPCMCalo("n3410113","0dm00009f9730000dge0404000","411790109fe30220000","0163103100b00010"); // INT7 30-40%
    cuts.AddCutPCMCalo("n4510113","0dm00009f9730000dge0404000","411790109fe30220000","0163103100b00010"); // INT7 40-50%
    cuts.AddCutPCMCalo("n5a10113","0dm00009f9730000dge0404000","411790109fe30220000","0163103100b00010"); // INT7 50-100%
  } else if (trainConfig == 2007) {  //   (fPSigSmearing, fPSigSmearingCte)
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe30220000","0163103100b00010"); // smearing (0.025,  0.030)
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe30220000","0163103100o00010"); // smearing new approach
  } else if ( trainConfig == 2008){ // EMCAL clusters
    cuts.AddCutPCMCalo("0008e113","0dm00009f9730000dge0404000","411790109fe30220000","0163103100000010"); // EDC
    cuts.AddCutPCMCalo("0008d113","0dm00009f9730000dge0404000","411790109fe30220000","0163103100000010"); // EDC
    cuts.AddCutPCMCalo("00083113","0dm00009f9730000dge0404000","111110109fe30220000","0163103100000010"); // EMCal only EG1
    cuts.AddCutPCMCalo("00085113","0dm00009f9730000dge0404000","111110109fe30220000","0163103100000010"); // EMCal only EG2
    cuts.AddCutPCMCalo("00089113","0dm00009f9730000dge0404000","388550109fe30220000","0163103100000010"); // DCal only DG1
    cuts.AddCutPCMCalo("0008b113","0dm00009f9730000dge0404000","388550109fe30220000","0163103100000010"); // DCal only DG2
  } else if (trainConfig == 2009){ //EMCal + DCal EG1 mult. diff.
    cuts.AddCutPCMCalo("q0176113","0dm00009f9730000dge0404000","411792106f032220000","01631031000000d0"); // INT7, NL12, high mult V0M, mult. dep 0% - 0.1%
    cuts.AddCutPCMCalo("q1276113","0dm00009f9730000dge0404000","411792106f032220000","01631031000000d0"); // INT7, NL12, high mult V0M, mult. dep 0.1% - 0.2%
    cuts.AddCutPCMCalo("q2376113","0dm00009f9730000dge0404000","411792106f032220000","01631031000000d0"); // INT7, NL12, high mult V0M, mult. dep 0.2% - 0.3%
    cuts.AddCutPCMCalo("q3476113","0dm00009f9730000dge0404000","411792106f032220000","01631031000000d0"); // INT7, NL12, high mult V0M, mult. dep 0.3% - 0.4%
    cuts.AddCutPCMCalo("q4576113","0dm00009f9730000dge0404000","411792106f032220000","01631031000000d0"); // INT7, NL12, high mult V0M, mult. dep 0.4% - 0.5%
    cuts.AddCutPCMCalo("q5776113","0dm00009f9730000dge0404000","411792106f032220000","01631031000000d0"); // INT7, NL12, high mult V0M, mult. dep 0.5% - 0.7%
    cuts.AddCutPCMCalo("q7a76113","0dm00009f9730000dge0404000","411792106f032220000","01631031000000d0"); // INT7, NL12, high mult V0M, mult. dep 0.7% - 1.0%
  } else if ( trainConfig == 2010){ // EMCAL clusters
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","111110109fe30220000","0163103100000010"); // EMCal only
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","388550109fe30220000","0163103100000010"); // DCal only
  } else if ( trainConfig == 2011){ // EMCAL clusters
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe30220000","0r63103100b00010"); // 90 degree rotation wo evt. weight
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe30220000","0t63103100b00010"); // random angle with evt. weight
  } else if ( trainConfig == 2012){ // EMCAL clusters
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe37220000","0163103100b00010"); //   min nCells = 2, E > 3
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe38220000","0163103100b00010"); //   min nCells = 2, E > 2
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe39220000","0163103100b00010"); //   min nCells = 2, E > 1
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe3a220000","0163103100b00010"); //   min nCells = 2, E > 5
  } else if ( trainConfig == 2013){ // EMCAL clusters
    cuts.AddCutPCMCalo("0008e113","0dm00009f9730000dge0404000","411790109fe32220000","0163103100000010"); // OLD NCell >= 2
    cuts.AddCutPCMCalo("0008e113","0dm00009f9730000dge0404000","411790106fe30220000","0163103100000010"); // OLD Timing
    cuts.AddCutPCMCalo("0008e113","0dm00009f9730000dge0404000","411790109fe30230000","0163103100000010"); // OLD M02
    cuts.AddCutPCMCalo("0008e113","0dm00009f9730000dge0404000","4117901097e30220000","0163103100000010"); // OLD CPV
  } else if ( trainConfig == 2014){ // EMCAL clusters
    cuts.AddCutPCMCalo("0008d113","0dm00009f9730000dge0404000","411790109fe32220000","0163103100000010"); // OLD NCell >= 2
    cuts.AddCutPCMCalo("0008d113","0dm00009f9730000dge0404000","411790106fe30220000","0163103100000010"); // OLD Timing
    cuts.AddCutPCMCalo("0008d113","0dm00009f9730000dge0404000","411790109fe30230000","0163103100000010"); // OLD M02
    cuts.AddCutPCMCalo("0008d113","0dm00009f9730000dge0404000","4117901097e30220000","0163103100000010"); // OLD CPV

  } else if ( trainConfig == 2015){ // min bias // EDC  // no nl
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411790009fe30220000","0163103100b00010"); // No NonLin
  } else if ( trainConfig == 2016){ // min bias // EDC  // NonLin settings
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe30220000","0163103100b00010"); // NonLin TB

  } else if (trainConfig == 2017){ // INT7, NL , STD without smearing
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411791109fe30220000","0163103100000010"); // INT7 NL11
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411791209fe30220000","0163103100000010"); // INT7 NL12
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411792109fe30220000","0163103100000010"); // INT7 NL21
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411792209fe30220000","0163103100000010"); // INT7 NL22
  } else if (trainConfig == 2018){ // EG2, NL , STD without smearing
    cuts.AddCutPCMCalo("0008e113","0dm00009f9730000dge0404000","411791109fe30220000","0163103100000010"); // EG2  NL11
    cuts.AddCutPCMCalo("0008e113","0dm00009f9730000dge0404000","411791209fe30220000","0163103100000010"); // EG2  NL12
    cuts.AddCutPCMCalo("0008e113","0dm00009f9730000dge0404000","411792109fe30220000","0163103100000010"); // EG2  NL21
    cuts.AddCutPCMCalo("0008e113","0dm00009f9730000dge0404000","411792209fe30220000","0163103100000010"); // EG2  NL22
  } else if (trainConfig == 2019){ // EG1, NL , STD without smearing
    cuts.AddCutPCMCalo("0008d113","0dm00009f9730000dge0404000","411791109fe30220000","0163103100000010"); // EG1  NL11
    cuts.AddCutPCMCalo("0008d113","0dm00009f9730000dge0404000","411791209fe30220000","0163103100000010"); // EG1  NL12
    cuts.AddCutPCMCalo("0008d113","0dm00009f9730000dge0404000","411792109fe30220000","0163103100000010"); // EG1  NL21
    cuts.AddCutPCMCalo("0008d113","0dm00009f9730000dge0404000","411792209fe30220000","0163103100000010"); // EG1  NL22

  } else if (trainConfig == 2020){ // INT7, NL , STD with smearing
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411791109fe30220000","0163103100b00010"); // INT7 NL11
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411791209fe30220000","0163103100b00010"); // INT7 NL12
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411792109fe30220000","0163103100b00010"); // INT7 NL21
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411792209fe30220000","0163103100b00010"); // INT7 NL22
  } else if (trainConfig == 2021){ // EG2, NL , STD with smearing
    cuts.AddCutPCMCalo("0008e113","0dm00009f9730000dge0404000","411791109fe30220000","0163103100b00010"); // EG2  NL11
    cuts.AddCutPCMCalo("0008e113","0dm00009f9730000dge0404000","411791209fe30220000","0163103100b00010"); // EG2  NL12
    cuts.AddCutPCMCalo("0008e113","0dm00009f9730000dge0404000","411792109fe30220000","0163103100b00010"); // EG2  NL21
    cuts.AddCutPCMCalo("0008e113","0dm00009f9730000dge0404000","411792209fe30220000","0163103100b00010"); // EG2  NL22
  } else if (trainConfig == 2022){ // EG1, NL , STD with smearing
    cuts.AddCutPCMCalo("0008d113","0dm00009f9730000dge0404000","411791109fe30220000","0163103100b00010"); // EG1  NL11
    cuts.AddCutPCMCalo("0008d113","0dm00009f9730000dge0404000","411791209fe30220000","0163103100b00010"); // EG1  NL12
    cuts.AddCutPCMCalo("0008d113","0dm00009f9730000dge0404000","411792109fe30220000","0163103100b00010"); // EG1  NL21
    cuts.AddCutPCMCalo("0008d113","0dm00009f9730000dge0404000","411792209fe30220000","0163103100b00010"); // EG1  NL22

    // Variations for systematics
    // Variations of EDC Part
  } else if (trainConfig == 2023){ // timing Cut variation  std -20+25ns
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411790106fe30220000","0163103100b00010"); //     -30    +35   ns
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","41179010afe30220000","0163103100b00010"); //     -12.5  +13   ns
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","41179010lfe30220000","0163103100b00010"); //     -12.5  +13   ns + timing effi
  } else if (trainConfig == 2024){ // track matching variation
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","4117901090e30220000","0163103100b00010"); //
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","4117901091e30220000","0163103100b00010"); //
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","4117901096e30220000","0163103100b00010"); //
  } else if (trainConfig == 2025){ // min nCells  // std: min nCells no cut
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe30220000","0163103100b00010"); //   no Cut
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe31220000","0163103100b00010"); //   min nCells = 1
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe32220000","0163103100b00010"); //   min nCells = 2
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe33220000","0163103100b00010"); //   min nCells = 3
  } else if (trainConfig == 2026){ // min energy variation std 0.7 GeV/c
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe00220000","0163103100b00010"); //     0.1 GeV/c
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe10220000","0163103100b00010"); //     0.5 GeV/c
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe20220000","0163103100b00010"); //     0.6 GeV/c
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe40220000","0163103100b00010"); //     0.8 GeV/c
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe50220000","0163103100b00010"); //     0.9 GeV/c
  } else if (trainConfig == 2027){ // M02 variation  // std: M02 max=0.7, min=0.1
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe30000000","0163103100b00010"); //   min max M02    = 0, 0
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe30210000","0163103100b00010"); //   min max M02    = 0.1 1
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe30230000","0163103100b00010"); //   min max M02    = 0.1 0.5
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe30250000","0163103100b00010"); //   min max M02    = 0.1 0.3
  } else if (trainConfig == 2028){ // exotic cluster  // STD: ExC 97 + TCard > 50
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411790109f230220000","0163103100b00010"); // ExC 99
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411790109f530220000","0163103100b00010"); // ExC 97
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411790109f830220000","0163103100b00010"); // ExC 95
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fb30220000","0163103100b00010"); // ExC 95 + TCard > 50
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411790109f030220000","0163103100b00010"); // no cut

  // Variations of PCM Part  // std 0dm00009f9730000dge0404000
  } else if (trainConfig == 2029) {
    cuts.AddCutPCMCalo("00010113","0dm00069f9730000dge0404000","411790109fe30220000","0163103100b00010"); // min pT 40 MeV
    cuts.AddCutPCMCalo("00010113","0dm00049f9730000dge0404000","411790109fe30220000","0163103100b00010"); // min pT 75 MeV
    cuts.AddCutPCMCalo("00010113","0dm00019f9730000dge0404000","411790109fe30220000","0163103100b00010"); // min pT 100MeV
  } else if (trainConfig == 2030) {
    cuts.AddCutPCMCalo("00010113","0dm00008f9730000dge0404000","411790109fe30220000","0163103100b00010"); // TPC cluster 35%
    cuts.AddCutPCMCalo("00010113","0dm00006f9730000dge0404000","411790109fe30220000","0163103100b00010"); // TPC cluster 70%
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0604000","411790109fe30220000","0163103100b00010"); // cosPA 0.9
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0304000","411790109fe30220000","0163103100b00010"); // cosPA 0.75
  } else if (trainConfig == 2031) {
    cuts.AddCutPCMCalo("00010113","0dm0000939730000dge0404000","411790109fe30220000","0163103100b00010"); // nsig electron   -4,5
    cuts.AddCutPCMCalo("00010113","0dm0000969730000dge0404000","411790109fe30220000","0163103100b00010"); // nsig electron -2.5,4
    cuts.AddCutPCMCalo("00010113","0dm00009f5730000dge0404000","411790109fe30220000","0163103100b00010"); // nsig pion 2,-10
    cuts.AddCutPCMCalo("00010113","0dm00009f1730000dge0404000","411790109fe30220000","0163103100b00010"); // nsig pion 0,-10
  } else if (trainConfig == 2032) {
    cuts.AddCutPCMCalo("00010113","0dm00009f9030000dge0404000","411790109fe30220000","0163103100b00010"); // pion nsig min mom 0.50 GeV/c
    cuts.AddCutPCMCalo("00010113","0dm00009f9630000dge0404000","411790109fe30220000","0163103100b00010"); // pion nsig min mom 0.25 GeV/c
    cuts.AddCutPCMCalo("00010113","0dm00009f9760000dge0404000","411790109fe30220000","0163103100b00010"); // pion nsig max mom 2.00 GeV/c
    cuts.AddCutPCMCalo("00010113","0dm00009f9710000dge0404000","411790109fe30220000","0163103100b00010"); // pion nsig max mom 5.00 GeV/c
  } else if (trainConfig == 2033){
    cuts.AddCutPCMCalo("00010113","0dm00009f97300008ge0404000","411790109fe30220000","0163103100b00010"); // qT max 0.05 1D
    cuts.AddCutPCMCalo("00010113","0dm00009f97300003ge0404000","411790109fe30220000","0163103100b00010"); // qT max 0.05 1D
    cuts.AddCutPCMCalo("00010113","0dm00009f97300002ge0404000","411790109fe30220000","0163103100b00010"); // qT max 0.06 2D
    cuts.AddCutPCMCalo("00010113","0dm00009f97300009ge0404000","411790109fe30220000","0163103100b00010"); // qT max 0.03 2D
  } else if (trainConfig == 2034) {
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dg50404000","411790109fe30220000","0163103100b00010"); // Psi pair 0.1  1D
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dg10404000","411790109fe30220000","0163103100b00010"); // Psi pair 0.1  1D
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dg60404000","411790109fe30220000","0163103100b00010"); // Psi pair 0.05 2D
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dg80404000","411790109fe30220000","0163103100b00010"); // Psi pair 0.2  2D
  } else if ( trainConfig == 2036){ // qT 2D pT dep
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe30220000","0163103100b00010"); // qT<0.110pT (2D) alpha<0.99
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000c259404000","411790109fe30220000","0163103100b00010"); // qT<0.110pT (2D) alpha<0.99
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000a259404000","411790109fe30220000","0163103100b00010"); // qT<0.125pT (2D) alpha<0.99
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000e259404000","411790109fe30220000","0163103100b00010"); // qT<0.130pT (2D) alpha<0.99
  } else if ( trainConfig == 2037){ // chi2, PsiPair cont.
    cuts.AddCutPCMCalo("00010113","0dm00009f97300008fd0404000","411790109fe30220000","0163103100b00010"); // PsiPair<0.15exp(-0.065chi2)
    cuts.AddCutPCMCalo("00010113","0dm00009f97300008ge0404000","411790109fe30220000","0163103100b00010"); // PsiPair<0.18exp(-0.055chi2)
    cuts.AddCutPCMCalo("00010113","0dm00009f97300008hf0404000","411790109fe30220000","0163103100b00010"); // PsiPair<0.20exp(-0.050chi2)
  } else if ( trainConfig == 2038){//alpha, std 3 == <=1.0
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe30220000","0163106100b00010"); //6: alpha meson variation 1 0<alpha<0.8
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe30220000","0163105100b00010"); //5: alpha meson variation 2 0<alpha<0.75
  } else if ( trainConfig == 2039){ //opening angle, std 1 == >0.5, 1 cell diagonal
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe30220000","0163103100b00000"); //0: min opening angle 0    -> open
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe30220000","0163103100b00030"); //3: min opening angle 0.01 -> 2 cell diag

  } else if ( trainConfig == 2040){ // EDC  ///   R Bins // without weights 1
    cuts.AddCutPCMCalo("00010113","0d200009f9730000dge0404000","411790109fe30220000","0163103100b00010"); // RBins    min = 5,      max = 180
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe30220000","0163103100b00010"); // RBins    min = 5,      max = 180 without 55 -72
    cuts.AddCutPCMCalo("00010113","0dh00009f9730000dge0404000","411790109fe30220000","0163103100b00010"); // RBins    min = 5,      max = 13
    cuts.AddCutPCMCalo("00010113","0di00009f9730000dge0404000","411790109fe30220000","0163103100b00010"); // RBins    min = 13,     max = 33.5
  } else if ( trainConfig == 2041){ // EDC  ///   R Bins // without weights 2
    cuts.AddCutPCMCalo("00010113","0dj00009f9730000dge0404000","411790109fe30220000","0163103100b00010"); // RBins    min = 33.5,   max = 55
    cuts.AddCutPCMCalo("00010113","0dk00009f9730000dge0404000","411790109fe30220000","0163103100b00010"); // RBins    min = 55,     max = 72
    cuts.AddCutPCMCalo("00010113","0dl00009f9730000dge0404000","411790109fe30220000","0163103100b00010"); // RBins    min = 72,     max = 95
    cuts.AddCutPCMCalo("00010113","0dg00009f9730000dge0404000","411790109fe30220000","0163103100b00010"); // RBins    min = 95,     max = 180
  } else if ( trainConfig == 2042){ // EDC  ///   R Bins // with weights 1
    cuts.AddCutPCMCalo("00010113","0d200009f9730000dge0404000","411790109fe30220000","0163103100b00010"); // RBins    min = 5,      max = 180
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe30220000","0163103100b00010"); // RBins    min = 5,      max = 180 without 55 -72
    cuts.AddCutPCMCalo("00010113","0dh00009f9730000dge0404000","411790109fe30220000","0163103100b00010"); // RBins    min = 5,      max = 13
    cuts.AddCutPCMCalo("00010113","0di00009f9730000dge0404000","411790109fe30220000","0163103100b00010"); // RBins    min = 13,     max = 33.5
  } else if ( trainConfig == 2043){ // EDC  ///   R Bins // with weights 2
    cuts.AddCutPCMCalo("00010113","0dj00009f9730000dge0404000","411790109fe30220000","0163103100b00010"); // RBins    min = 33.5,   max = 55
    cuts.AddCutPCMCalo("00010113","0dk00009f9730000dge0404000","411790109fe30220000","0163103100b00010"); // RBins    min = 55,     max = 72
    cuts.AddCutPCMCalo("00010113","0dl00009f9730000dge0404000","411790109fe30220000","0163103100b00010"); // RBins    min = 72,     max = 95
    cuts.AddCutPCMCalo("00010113","0dg00009f9730000dge0404000","411790109fe30220000","0163103100b00010"); // RBins    min = 95,     max = 180
  } else if ( trainConfig == 2044){ // EMCAL only //   R Bins // without weights 1
    cuts.AddCutPCMCalo("00010113","0d200009f9730000dge0404000","111110109fe30220000","0163103100b00010"); // RBins    min = 5,      max = 180
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","111110109fe30220000","0163103100b00010"); // RBins    min = 5,      max = 180 without 55 -72
    cuts.AddCutPCMCalo("00010113","0dh00009f9730000dge0404000","111110109fe30220000","0163103100b00010"); // RBins    min = 5,      max = 13
    cuts.AddCutPCMCalo("00010113","0di00009f9730000dge0404000","111110109fe30220000","0163103100b00010"); // RBins    min = 13,     max = 33.5
  } else if ( trainConfig == 2045){ // EMCAL only //   R Bins // without weights 2
    cuts.AddCutPCMCalo("00010113","0dj00009f9730000dge0404000","111110109fe30220000","0163103100b00010"); // RBins    min = 33.5,   max = 55
    cuts.AddCutPCMCalo("00010113","0dk00009f9730000dge0404000","111110109fe30220000","0163103100b00010"); // RBins    min = 55,     max = 72
    cuts.AddCutPCMCalo("00010113","0dl00009f9730000dge0404000","111110109fe30220000","0163103100b00010"); // RBins    min = 72,     max = 95
    cuts.AddCutPCMCalo("00010113","0dg00009f9730000dge0404000","111110109fe30220000","0163103100b00010"); // RBins    min = 95,     max = 180
  } else if ( trainConfig == 2046){ // EMCAL only //   R Bins // with weights 1
    cuts.AddCutPCMCalo("00010113","0d200009f9730000dge0404000","111110109fe30220000","0163103100b00010"); // RBins    min = 5,      max = 180
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","111110109fe30220000","0163103100b00010"); // RBins    min = 5,      max = 180 without 55 -72
    cuts.AddCutPCMCalo("00010113","0dh00009f9730000dge0404000","111110109fe30220000","0163103100b00010"); // RBins    min = 5,      max = 13
    cuts.AddCutPCMCalo("00010113","0di00009f9730000dge0404000","111110109fe30220000","0163103100b00010"); // RBins    min = 13,     max = 33.5
  } else if ( trainConfig == 2047){ // EMCAL only //   R Bins // with weights 2
    cuts.AddCutPCMCalo("00010113","0dj00009f9730000dge0404000","111110109fe30220000","0163103100b00010"); // RBins    min = 33.5,   max = 55
    cuts.AddCutPCMCalo("00010113","0dk00009f9730000dge0404000","111110109fe30220000","0163103100b00010"); // RBins    min = 55,     max = 72
    cuts.AddCutPCMCalo("00010113","0dl00009f9730000dge0404000","111110109fe30220000","0163103100b00010"); // RBins    min = 72,     max = 95
    cuts.AddCutPCMCalo("00010113","0dg00009f9730000dge0404000","111110109fe30220000","0163103100b00010"); // RBins    min = 95,     max = 180
  } else if ( trainConfig == 2048){ // DCAL only //   R Bins // without weights 1
    cuts.AddCutPCMCalo("00010113","0d200009f9730000dge0404000","388550109fe30220000","0163103100b00010"); // RBins    min = 5,      max = 180
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","388550109fe30220000","0163103100b00010"); // RBins    min = 5,      max = 180 without 55 -72
    cuts.AddCutPCMCalo("00010113","0dh00009f9730000dge0404000","388550109fe30220000","0163103100b00010"); // RBins    min = 5,      max = 13
    cuts.AddCutPCMCalo("00010113","0di00009f9730000dge0404000","388550109fe30220000","0163103100b00010"); // RBins    min = 13,     max = 33.5
  } else if ( trainConfig == 2049){ // DCAL only //   R Bins // without weights 2
    cuts.AddCutPCMCalo("00010113","0dj00009f9730000dge0404000","388550109fe30220000","0163103100b00010"); // RBins    min = 33.5,   max = 55
    cuts.AddCutPCMCalo("00010113","0dk00009f9730000dge0404000","388550109fe30220000","0163103100b00010"); // RBins    min = 55,     max = 72
    cuts.AddCutPCMCalo("00010113","0dl00009f9730000dge0404000","388550109fe30220000","0163103100b00010"); // RBins    min = 72,     max = 95
    cuts.AddCutPCMCalo("00010113","0dg00009f9730000dge0404000","388550109fe30220000","0163103100b00010"); // RBins    min = 95,     max = 180
  } else if ( trainConfig == 2050){ // DCAL only //   R Bins // with weights 1
    cuts.AddCutPCMCalo("00010113","0d200009f9730000dge0404000","388550109fe30220000","0163103100b00010"); // RBins    min = 5,      max = 180
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","388550109fe30220000","0163103100b00010"); // RBins    min = 5,      max = 180 without 55 -72
    cuts.AddCutPCMCalo("00010113","0dh00009f9730000dge0404000","388550109fe30220000","0163103100b00010"); // RBins    min = 5,      max = 13
    cuts.AddCutPCMCalo("00010113","0di00009f9730000dge0404000","388550109fe30220000","0163103100b00010"); // RBins    min = 13,     max = 33.5
  } else if ( trainConfig == 2051){ // DCAL only //   R Bins // with weights 2
    cuts.AddCutPCMCalo("00010113","0dj00009f9730000dge0404000","388550109fe30220000","0163103100b00010"); // RBins    min = 33.5,   max = 55
    cuts.AddCutPCMCalo("00010113","0dk00009f9730000dge0404000","388550109fe30220000","0163103100b00010"); // RBins    min = 55,     max = 72
    cuts.AddCutPCMCalo("00010113","0dl00009f9730000dge0404000","388550109fe30220000","0163103100b00010"); // RBins    min = 72,     max = 95
    cuts.AddCutPCMCalo("00010113","0dg00009f9730000dge0404000","388550109fe30220000","0163103100b00010"); // RBins    min = 95,     max = 180

  } else if (trainConfig == 2060){ // cluster swapping for background calculation
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe30220000","0r63103100b00010"); // INT7 TBNL
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe30220000","0s63103100b00010"); // INT7 TBNL
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe30220000","0t63103100b00010"); // INT7 TBNL
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe30220000","0u63103100b00010"); // INT7 TBNL

  } else if (trainConfig == 2070){  // EMCal+DCAL MB  TB NL tests, only TBNL (no MC finetuning)
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe30220000","0163103100000010"); // INT7
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411793709fe30220000","0163103100000010"); // INT7 NL37
  } else if (trainConfig == 2071){  // EMCal+DCAL EG2 TB NL tests, only TBNL (no MC finetuning)
    cuts.AddCutPCMCalo("0008e113","0dm00009f9730000dge0404000","411790109fe30220000","0163103100000010"); // EG2
    cuts.AddCutPCMCalo("0008e113","0dm00009f9730000dge0404000","411793709fe30220000","0163103100000010"); // EG2 NL37
  } else if (trainConfig == 2072){ // EMCal+DCAL EG1 TB NL tests, only TBNL (no MC finetuning)
    cuts.AddCutPCMCalo("0008d113","0dm00009f9730000dge0404000","411790109fe30220000","0163103100000010"); // EG1
    cuts.AddCutPCMCalo("0008d113","0dm00009f9730000dge0404000","411793709fe30220000","0163103100000010"); // EG1 NL37

  // configs without smearing
  } else if (trainConfig == 2080){ // EMCAL+DCAL clusters standard cuts, INT7, NL , std TM,
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe30220000","0163103100000010"); // INT7 no smearing
  } else if (trainConfig == 2081){ // EMCAL+DCAL clusters standard cuts, INT7, NL , std TM,
    cuts.AddCutPCMCalo("0008e113","0dm00009f9730000dge0404000","411790109fe30220000","0163103100000010"); // EG2 no smearing
  } else if (trainConfig == 2082){ // EMCAL+DCAL clusters standard cuts, INT7, NL , std TM,
    cuts.AddCutPCMCalo("0008d113","0dm00009f9730000dge0404000","411790109fe30220000","0163103100000010"); // EG1 no smearing

  // configs with smearing
  } else if (trainConfig == 2083){ // EMCAL+DCAL clusters standard cuts, INT7, NL , std TM,
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe30220000","0163103100b00010"); // INT7 smearing
  } else if (trainConfig == 2084){ // EMCAL+DCAL clusters standard cuts, INT7, NL , std TM,
    cuts.AddCutPCMCalo("0008e113","0dm00009f9730000dge0404000","411790109fe30220000","0163103100b00010"); // EG2 smearing
  } else if (trainConfig == 2085){ // EMCAL+DCAL clusters standard cuts, INT7, NL , std TM,
    cuts.AddCutPCMCalo("0008d113","0dm00009f9730000dge0404000","411790109fe30220000","0163103100b00010"); // EG1 smearing

  // configs with and without smearing
  } else if (trainConfig == 2086){ // EMCAL+DCAL clusters standard cuts, INT7, NL , std TM,
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe30220000","0163103100b00010"); // INT7 smearing
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe30220000","0163103100000010"); // INT7 no smearing
  } else if (trainConfig == 2087){ // EMCAL+DCAL clusters standard cuts, INT7, NL , std TM,
    cuts.AddCutPCMCalo("0008e113","0dm00009f9730000dge0404000","411790109fe30220000","0163103100b00010"); // EG2 smearing
    cuts.AddCutPCMCalo("0008e113","0dm00009f9730000dge0404000","411790109fe30220000","0163103100000010"); // EG2 no smearing
  } else if (trainConfig == 2088){ // EMCAL+DCAL clusters standard cuts, INT7, NL , std TM,
    cuts.AddCutPCMCalo("0008d113","0dm00009f9730000dge0404000","411790109fe30220000","0163103100b00010"); // EG1 smearing
    cuts.AddCutPCMCalo("0008d113","0dm00009f9730000dge0404000","411790109fe30220000","0163103100000010"); // EG1 no smearing

  // configs without NonLin settings (for validation of EMC CF)
  } else if (trainConfig == 2089){ // EMCAL+DCAL clusters standard cuts, INT7, NL , std TM,
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411790009fe32220000","0163103100b00010"); // INT7
  } else if (trainConfig == 2090){ // EMCAL+DCAL clusters standard cuts, INT7, NL , std TM,
    cuts.AddCutPCMCalo("0008e113","0dm00009f9730000dge0404000","411790009fe32220000","0163103100b00010"); // EG2
  } else if (trainConfig == 2091){ // EMCAL+DCAL clusters standard cuts, INT7, NL , std TM,
    cuts.AddCutPCMCalo("0008d113","0dm00009f9730000dge0404000","411790009fe32220000","0163103100b00010"); // EG1

  // configs with and without smearing NCell >= 2
  } else if (trainConfig == 2092){ // EMCAL+DCAL clusters standard cuts, INT7, NL , std TM,
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe32220000","0163103100b00010"); // INT7 smearing
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe32220000","0163103100000010"); // INT7 no smearing
  } else if (trainConfig == 2093){ // EMCAL+DCAL clusters standard cuts, INT7, NL , std TM,
    cuts.AddCutPCMCalo("0008e113","0dm00009f9730000dge0404000","411790109fe32220000","0163103100b00010"); // EG2 smearing
    cuts.AddCutPCMCalo("0008e113","0dm00009f9730000dge0404000","411790109fe32220000","0163103100000010"); // EG2 no smearing
  } else if (trainConfig == 2094){ // EMCAL+DCAL clusters standard cuts, INT7, NL , std TM,
    cuts.AddCutPCMCalo("0008d113","0dm00009f9730000dge0404000","411790109fe32220000","0163103100b00010"); // EG1 smearing
    cuts.AddCutPCMCalo("0008d113","0dm00009f9730000dge0404000","411790109fe32220000","0163103100000010"); // EG1 no smearing

  // configs with MC finetuning
  } else if (trainConfig == 2100){  // EMCal+DCAL clusters standard cuts, triggers, NL kSDM, tight timing, E/p TM
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0400000","4117918077032230000","0163103100000010"); // INT7
  } else if (trainConfig == 2101){  // EMCal+DCAL clusters standard cuts, triggers, NL DExp, tight timing, E/p TM
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0400000","4117928077032230000","0163103100000010"); // INT7
  } else if (trainConfig == 2102){  // EMCal+DCAL clusters standard cuts, triggers, NL kSDM sep for EMC/DMC, tight timing, E/p TM
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0400000","4117911077032230000","0163103100000010"); // INT7
  } else if (trainConfig == 2103){  // EMCal+DCAL clusters standard cuts, triggers, NL kSDM, tight timing, E/p TM
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0400000","4117917077032230000","0163103100000010"); // INT7
  } else if (trainConfig == 2104){  // EMCal+DCAL clusters standard cuts, triggers, NL kSDM, tight timing, E/p TM
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0400000","4117901077032230000","0163103100000010"); // INT7
  } else if (trainConfig == 2105){  // EMCal+DCAL clusters standard cuts, triggers, NL kSDM, tight timing, E/p TM E>0.6 for NL
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0400000","4117901077022230000","0163103100000010"); // INT7

  } else if (trainConfig == 2110){  // EMCal+DCAL clusters standard cuts, triggers, NL kSDM, tight timing, E/p TM  // with eta<0.8
    cuts.AddCutPCMCalo("00010113","0d200009327000008250400000","4117918077032230000","0163103100000010"); // INT7
  } else if (trainConfig == 2111){  // EMCal+DCAL clusters standard cuts, triggers, NL kSDM, tight timing, E/p TM  // with eta<0.8
    cuts.AddCutPCMCalo("00010113","0da00009327000008250400000","4117918077032230000","0163103100000010"); // INT7
    cuts.AddCutPCMCalo("00010113","0db00009327000008250400000","4117918077032230000","0163103100000010"); // INT7
    cuts.AddCutPCMCalo("00010113","0dc00009327000008250400000","4117918077032230000","0163103100000010"); // INT7
  } else if (trainConfig == 2112){  // EMCal+DCAL clusters standard cuts, triggers, NL kSDM, tight timing, E/p TM  // with eta<0.8
    cuts.AddCutPCMCalo("00010113","0dh00009327000008250400000","4117918077032230000","0163103100000010"); // INT7
    cuts.AddCutPCMCalo("00010113","0di00009327000008250400000","4117918077032230000","0163103100000010"); // INT7

  } else if (trainConfig == 2120){  // EMCal+DCAL clusters standard cuts, triggers, NL kSDM, tight timing, E/p TM  // with eta<0.8 // to be used with MBW
    cuts.AddCutPCMCalo("00010113","0d200009327000008250400000","4117918077032230000","0163103100000010"); // INT7
  } else if (trainConfig == 2121){  // EMCal+DCAL clusters standard cuts, triggers, NL kSDM, tight timing, E/p TM  // with eta<0.8
    cuts.AddCutPCMCalo("00010113","0da00009327000008250400000","4117918077032230000","0163103100000010"); // INT7
    cuts.AddCutPCMCalo("00010113","0db00009327000008250400000","4117918077032230000","0163103100000010"); // INT7
    cuts.AddCutPCMCalo("00010113","0dc00009327000008250400000","4117918077032230000","0163103100000010"); // INT7
  } else if (trainConfig == 2122){  // EMCal+DCAL clusters standard cuts, triggers, NL kSDM, tight timing, E/p TM  // with eta<0.8
    cuts.AddCutPCMCalo("00010113","0dh00009327000008250400000","4117918077032230000","0163103100000010"); // INT7
    cuts.AddCutPCMCalo("00010113","0di00009327000008250400000","4117918077032230000","0163103100000010"); // INT7

  // NCell efficiency pp 13TeV
  } else if (trainConfig == 2130){ // EMCAL+DCAL clusters standard cuts, INT7, NL , std TM,
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe3h220000","0163103100000010"); // NCell effi on MC
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe3i220000","0163103100000010"); // NCell effi on MC
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe3l220000","0163103100000010"); // NCell effi on MC
  } else if (trainConfig == 2131){ // EMCAL+DCAL clusters standard cuts, INT7, NL , std TM,
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe30220000","0163103100000010"); // No NCell cut (NCell >= 1)
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe32220000","0163103100000010"); // NCell cut >= 2
  } else if (trainConfig == 2132){ // EMCAL+DCAL clusters standard cuts, INT7, NL , std TM, new NCell calculation
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe3m220000","0163103100000010"); // NCell effi on MC
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe3n220000","0163103100000010"); // NCell effi on MC
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe3o220000","0163103100000010"); // NCell effi on MC
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe3p220000","0163103100000010"); // NCell effi on MC
  } else if (trainConfig == 2133){ // EMCAL+DCAL clusters standard cuts, INT7, NL , std TM, new NCell calculation TB variations
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe3q220000","0163103100000010"); // NCell effi on MC
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe3r220000","0163103100000010"); // NCell effi on MC
  } else if (trainConfig == 2134){ // 1 cell clusters only  (switch off exotics + M02 cut)
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411790109f03s000000","0163103100000010"); // NL 01
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411793709f03s000000","0163103100000010"); // TBNL
  } else if (trainConfig == 2135){ // 2 cell clusters only  (switch off exotics + M02 cut)
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411790109f03t000000","0163103100000010"); // NL 01
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411793709f03t000000","0163103100000010"); // TBNL
  } else if (trainConfig == 2136){ // 2 cell clusters only  (with exotics + M02 cut)
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe3t220000","0163103100000010"); // NL 01
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411793709fe3t220000","0163103100000010"); // TBNL
  } else if (trainConfig == 2137){ // >=2 cell clusters  (switch off exotics + M02 cut)
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411790109f032000000","0163103100000010"); // NL 01
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411793709f032000000","0163103100000010"); // TBNL
  } else if (trainConfig == 2138){ // >=2 cell clusters (with exotics + M02 cut)
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe32220000","0163103100000010"); // NL 01
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411793709fe32220000","0163103100000010"); // TBNL

    // different clusterization settings
  } else if (trainConfig == 2139){ // std. cuts, Seed 500, agg 85
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe32220000","0163103100000010"); // NL 01
  } else if (trainConfig == 2140){ // std. cuts, Seed 500, agg 90
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe32220000","0163103100000010"); // NL 01
  } else if (trainConfig == 2141){ // std. cuts, Seed 475, Agg 95
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe32220000","0163103100000010"); // NL 01

  } else if (trainConfig == 2142){ // Effi on isolated clusters only
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe3u220000","0163103100000010"); // TB correction
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe3v220000","0163103100000010"); // P2 correction
  } else if (trainConfig == 2143){ // Effi on isolated clusters only
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fh3m220000","0163103100000010"); //  Exotics > 2 GeV
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fi3m220000","0163103100000010"); //  Exotics > 3 GeV
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fj3m220000","0163103100000010"); //  Exotics > 4 GeV
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fk3m220000","0163103100000010"); //  Exotics > 5 GeV
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fl3m220000","0163103100000010"); //  Exotics > 6 GeV

  // PCM-EDC systematics
  } else if (trainConfig == 2200){ // PCM based systematics
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","4117918077032230000","0163103100000010"); // std
    cuts.AddCutPCMCalo("00010013","00200009327000008250400000","4117918077032230000","0163103100000010"); // std - no pileup cut
    cuts.AddCutPCMCalo("00010013","00202209327000008250400000","4117918077032230000","0163103100000010"); // restrict acceptance to EMCAL loose
    cuts.AddCutPCMCalo("00010013","00204409327000008250400000","4117918077032230000","0163103100000010"); // restrict acceptance to EMCAL tight
  } else if (trainConfig == 2201){ // PCM based systematics
    cuts.AddCutPCMCalo("00010113","00200009327400008250400000","4117918077032230000","0163103100000010"); // dEdx pi: 1: 0.4-3, -10: 3. ->
    cuts.AddCutPCMCalo("00010113","00200009367400008250400000","4117918077032230000","0163103100000010"); // dEdx pi: 2: 0.4-3, 0.5: 3. ->
    cuts.AddCutPCMCalo("00010113","00200009227000008250400000","4117918077032230000","0163103100000010"); // dEdx e -3, 5
    cuts.AddCutPCMCalo("00010113","00200009127000008250400000","4117918077032230000","0163103100000010"); // dEdx e -5, 5
  } else if (trainConfig == 2202){ // PCM based systematics
    cuts.AddCutPCMCalo("00010113","00200049327000008250400000","4117918077032230000","0163103100000010"); // single pt > 0.075
    cuts.AddCutPCMCalo("00010113","00200019327000008250400000","4117918077032230000","0163103100000010"); // single pt > 0.1
    cuts.AddCutPCMCalo("00010113","00200009327000009250400000","4117918077032230000","0163103100000010"); // qt 2D 0.03
    cuts.AddCutPCMCalo("00010113","00200009327000003250400000","4117918077032230000","0163103100000010"); // qt 1D 0.05
    cuts.AddCutPCMCalo("00010113","00200009327000002250400000","4117918077032230000","0163103100000010"); // qt 1D 0.07
  } else if (trainConfig == 2203){ // PCM based systematics
    cuts.AddCutPCMCalo("00010113","00200009327000008850400000","4117918077032230000","0163103100000010"); // 2D psi pair chi2 var
    cuts.AddCutPCMCalo("00010113","00200009327000008260400000","4117918077032230000","0163103100000010"); // 2D psi pair chi2 var
    cuts.AddCutPCMCalo("00010113","00200009327000008860400000","4117918077032230000","0163103100000010"); // 2D psi pair chi2 var
    cuts.AddCutPCMCalo("00010113","00200009327000008280400000","4117918077032230000","0163103100000010"); // 2D psi pair chi2 var
    cuts.AddCutPCMCalo("00010113","00200009327000008880400000","4117918077032230000","0163103100000010"); // 2D psi pair chi2 var
  } else if (trainConfig == 2204){ // PCM based systematics
    cuts.AddCutPCMCalo("00010113","00200006327000008250400000","4117918077032230000","0163103100000010"); // min TPC cl > 0.7
    cuts.AddCutPCMCalo("00010113","00200008327000008250400000","4117918077032230000","0163103100000010"); // min TPC cl > 0.35
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","4117918077032230000","0163106100000010"); // alpha < 0.8
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","4117918077032230000","0163105100000010"); // alpha < 0.75

  } else if (trainConfig == 2210){ // EMC min E variations
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","4117918077022230000","0163103100000010"); //0.6 GeV/c
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","4117918077042230000","0163103100000010"); //0.8 GeV/c
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","4117918077052230000","0163103100000010"); //0.9 GeV/c
  } else if (trainConfig == 2211){ //EMCAL minNCells variation
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","4117918077031230000","0163103100000010"); //n cells >= 1
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","4117918077033230000","0163103100000010"); //n cells >= 3
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","4117918077032200000","0163103100000010"); //no max M02 cut
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","4117918077032250000","0163103100000010"); //M02 < 0.3
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","41179180770322k0000","0163103100000010"); //M02, pT-dep
  } else if (trainConfig == 2212){ // EMCAL track matching variations
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","4117918076032230000","0163103100000010"); // track matching variations
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","4117918078032230000","0163103100000010"); //
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","411791807f032230000","0163103100000010"); //
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","4117918073032230000","0163103100000010"); //
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","4117918070032230000","0163103100000010"); //
  } else if (trainConfig == 2213){ // EMCAL clusters, timing variation
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","4117918057032230000","0163103100000010"); // time 50ns
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","4117918047032230000","0163103100000010"); // time 100ns
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","4117918037032230000","0163103100000010"); // time 200ns
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","4117918027032230000","0163103100000010"); // time 500ns
  } else if (trainConfig == 2214){ // EMCAL clusters, timing variation
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","41179180a7032230000","0163103100000010"); // time
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","4117918087032230000","0163103100000010"); // time
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","4117918097032230000","0163103100000010"); // time

  } else if (trainConfig == 2220){  // EMCal+DCAL clusters standard cuts, Sphericity
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0400000","4117918077032230000","0h63103100000010"); // INT7
    cuts.AddCutPCMCalo("h0510113","00200009f9730000dge0400000","4117918077032230000","0h63103100000010"); // INT7
    cuts.AddCutPCMCalo("h5a10113","00200009f9730000dge0400000","4117918077032230000","0h63103100000010"); // INT7
    cuts.AddCutPCMCalo("h0a10113","00200009f9730000dge0400000","4117918077032230000","0h63103100000010"); // INT7
    cuts.AddCutPCMCalo("h0310113","00200009f9730000dge0400000","4117918077032230000","0h63103100000010"); // INT7
    cuts.AddCutPCMCalo("h7a10113","00200009f9730000dge0400000","4117918077032230000","0h63103100000010"); // INT7
  } else if (trainConfig == 2221){  // EMCal+DCAL clusters standard cuts, V0M mult selections
    cuts.AddCutPCMCalo("m0110113","00200009f9730000dge0400000","4117918077032230000","0h63103100000010"); // INT7 0-1%
    cuts.AddCutPCMCalo("m1510113","00200009f9730000dge0400000","4117918077032230000","0h63103100000010"); // INT7 1-5%
    cuts.AddCutPCMCalo("m5k10113","00200009f9730000dge0400000","4117918077032230000","0h63103100000010"); // INT7 5-20%
    cuts.AddCutPCMCalo("n2410113","00200009f9730000dge0400000","4117918077032230000","0h63103100000010"); // INT7 20-40%
    cuts.AddCutPCMCalo("n4710113","00200009f9730000dge0400000","4117918077032230000","0h63103100000010"); // INT7 40-70%
    cuts.AddCutPCMCalo("n7a10113","00200009f9730000dge0400000","4117918077032230000","0h63103100000010"); // INT7 70-100%
  } else if (trainConfig == 2222){  // EMCal+DCAL clusters standard cuts, SPD mult selections
    cuts.AddCutPCMCalo("o0110113","00200009f9730000dge0400000","4117918077032230000","0h63103100000010"); // INT7 0-1%
    cuts.AddCutPCMCalo("o0210113","00200009f9730000dge0400000","4117918077032230000","0h63103100000010"); // INT7 0-2%
    cuts.AddCutPCMCalo("o0510113","00200009f9730000dge0400000","4117918077032230000","0h63103100000010"); // INT7 0-5%
    cuts.AddCutPCMCalo("o5k10113","00200009f9730000dge0400000","4117918077032230000","0h63103100000010"); // INT7 5-20%
    cuts.AddCutPCMCalo("p2610113","00200009f9730000dge0400000","4117918077032230000","0h63103100000010"); // INT7 20-60%
    cuts.AddCutPCMCalo("p6a10113","00200009f9730000dge0400000","4117918077032230000","0h63103100000010"); // INT7 60-100%

 // special V0AND configurations for PCM cut QA
  } else if (trainConfig == 2250){  // std cut
    cuts.AddCutPCMCalo("00010113","00200009327300008250400000","4117918077032230000","0163103100000010"); // INT7
  } else if (trainConfig == 2251){  //
    cuts.AddCutPCMCalo("00010113","0020000932730000b250400000","4117918077032230000","0163103100000010"); // qT<0.125pT (1D)
    cuts.AddCutPCMCalo("00010113","0020000932730000a250400000","4117918077032230000","0163103100000010"); // qT<0.125pT (2D) alpha<0.95
    cuts.AddCutPCMCalo("00010113","0020000932730000a259400000","4117918077032230000","0163103100000010"); // qT<0.125pT (2D) alpha<0.99
    cuts.AddCutPCMCalo("00010113","0020000932730000a25a400000","4117918077032230000","0163103100000010"); // qT<0.125pT (2D) alpha<1
  } else if (trainConfig == 2252){  //
    cuts.AddCutPCMCalo("00010113","0020000932730000acc0400000","4117918077032230000","0163103100000010"); // PsiPair<0.15 + chi2<40 + qT<0.125pT (2D) alpha<0.95
    cuts.AddCutPCMCalo("00010113","0020000932730000a180400000","4117918077032230000","0163103100000010"); // PsiPair<0.20 + chi2<50 + qT<0.125pT (2D) alpha<0.95
    cuts.AddCutPCMCalo("00010113","0020000932730000acca400000","4117918077032230000","0163103100000010"); // PsiPair<0.15 + chi2<40 + qT<0.125pT (2D) alpha<1
    cuts.AddCutPCMCalo("00010113","0020000932730000a18a400000","4117918077032230000","0163103100000010"); // PsiPair<0.20 + chi2<50 + qT<0.125pT (2D) alpha<1
  } else if (trainConfig == 2253){  //
    cuts.AddCutPCMCalo("00010113","00200009227300008fd0400000","4117918077032230000","0163103100000010"); // PsiPair<0.15exp(-0.065chi2)
    cuts.AddCutPCMCalo("00010113","00200009227300008ge0400000","4117918077032230000","0163103100000010"); // PsiPair<0.18exp(-0.055chi2)
    cuts.AddCutPCMCalo("00010113","00200009227300008hf0400000","4117918077032230000","0163103100000010"); // PsiPair<0.20exp(-0.050chi2)
  } else if (trainConfig == 2254){  //
    cuts.AddCutPCMCalo("00010113","0020000922730000afd0400000","4117918077032230000","0163103100000010"); // PsiPair<0.15exp(-0.065chi2) + qT<0.125pT (2D) alpha<0.95
    cuts.AddCutPCMCalo("00010113","0020000922730000age0400000","4117918077032230000","0163103100000010"); // PsiPair<0.18exp(-0.055chi2) + qT<0.125pT (2D) alpha<0.95
    cuts.AddCutPCMCalo("00010113","0020000922730000ahf0400000","4117918077032230000","0163103100000010"); // PsiPair<0.20exp(-0.050chi2) + qT<0.125pT (2D) alpha<0.95
  } else if (trainConfig == 2255){  //
    cuts.AddCutPCMCalo("00010113","0020000922730000gfd0400000","4117918077032230000","0163103100000010"); // PsiPair<0.15exp(-0.065chi2) + qT<0.140pT (2D) alpha<0.95
    cuts.AddCutPCMCalo("00010113","0020000922730000gge0400000","4117918077032230000","0163103100000010"); // PsiPair<0.18exp(-0.055chi2) + qT<0.140pT (2D) alpha<0.95
    cuts.AddCutPCMCalo("00010113","0020000922730000ghf0400000","4117918077032230000","0163103100000010"); // PsiPair<0.20exp(-0.050chi2) + qT<0.140pT (2D) alpha<0.95
  } else if (trainConfig == 2256){  //
    cuts.AddCutPCMCalo("00010113","0020000922730000gfda400000","4117918077032230000","0163103100000010"); // PsiPair<0.15exp(-0.065chi2) + qT<0.140pT (2D) alpha<1
    cuts.AddCutPCMCalo("00010113","0020000922730000ggea400000","4117918077032230000","0163103100000010"); // PsiPair<0.18exp(-0.055chi2) + qT<0.140pT (2D) alpha<1
    cuts.AddCutPCMCalo("00010113","0020000922730000ghfa400000","4117918077032230000","0163103100000010"); // PsiPair<0.20exp(-0.050chi2) + qT<0.140pT (2D) alpha<1
  // EMC7
  } else if (trainConfig == 2260){  // std cut
    cuts.AddCutPCMCalo("00052113","00200009327300008250400000","4117918077032230000","0163103100000010"); // EMC7
  } else if (trainConfig == 2261){  // std cut
    cuts.AddCutPCMCalo("00052113","0020000932730000b250400000","4117918077032230000","0163103100000010"); // qT<0.125pT (1D)
    cuts.AddCutPCMCalo("00052113","0020000932730000a250400000","4117918077032230000","0163103100000010"); // qT<0.125pT (2D) alpha<0.95
    cuts.AddCutPCMCalo("00052113","0020000932730000a259400000","4117918077032230000","0163103100000010"); // qT<0.125pT (2D) alpha<0.99
    cuts.AddCutPCMCalo("00052113","0020000932730000a25a400000","4117918077032230000","0163103100000010"); // qT<0.125pT (2D) alpha<1
  } else if (trainConfig == 2262){  // std cut
    cuts.AddCutPCMCalo("00052113","0020000932730000acc0400000","4117918077032230000","0163103100000010"); // PsiPair<0.15 + chi2<40 + qT<0.125pT (2D) alpha<0.95
    cuts.AddCutPCMCalo("00052113","0020000932730000a180400000","4117918077032230000","0163103100000010"); // PsiPair<0.20 + chi2<50 + qT<0.125pT (2D) alpha<0.95
    cuts.AddCutPCMCalo("00052113","0020000932730000acca400000","4117918077032230000","0163103100000010"); // PsiPair<0.15 + chi2<40 + qT<0.125pT (2D) alpha<1
    cuts.AddCutPCMCalo("00052113","0020000932730000a18a400000","4117918077032230000","0163103100000010"); // PsiPair<0.20 + chi2<50 + qT<0.125pT (2D) alpha<1
  // EGA
  } else if (trainConfig == 2270){  // std cut
    cuts.AddCutPCMCalo("00081113","00200009327300008250400000","4117918077032230000","0163103100000010"); // EGA
  } else if (trainConfig == 2271){  // std cut
    cuts.AddCutPCMCalo("00081113","0020000932730000b250400000","4117918077032230000","0163103100000010"); // qT<0.125pT (1D)
    cuts.AddCutPCMCalo("00081113","0020000932730000a250400000","4117918077032230000","0163103100000010"); // qT<0.125pT (2D) alpha<0.95
    cuts.AddCutPCMCalo("00081113","0020000932730000a259400000","4117918077032230000","0163103100000010"); // qT<0.125pT (2D) alpha<0.99
    cuts.AddCutPCMCalo("00081113","0020000932730000a25a400000","4117918077032230000","0163103100000010"); // qT<0.125pT (2D) alpha<1
  } else if (trainConfig == 2272){  // std cut
    cuts.AddCutPCMCalo("00081113","0020000932730000acc0400000","4117918077032230000","0163103100000010"); // PsiPair<0.15 + chi2<40 + qT<0.125pT (2D) alpha<0.95
    cuts.AddCutPCMCalo("00081113","0020000932730000a180400000","4117918077032230000","0163103100000010"); // PsiPair<0.20 + chi2<50 + qT<0.125pT (2D) alpha<0.95
    cuts.AddCutPCMCalo("00081113","0020000932730000acca400000","4117918077032230000","0163103100000010"); // PsiPair<0.15 + chi2<40 + qT<0.125pT (2D) alpha<1
    cuts.AddCutPCMCalo("00081113","0020000932730000a18a400000","4117918077032230000","0163103100000010"); // PsiPair<0.20 + chi2<50 + qT<0.125pT (2D) alpha<1
  // EG1+DG1
  } else if (trainConfig == 2270){  // std cut
    cuts.AddCutPCMCalo("0008d113","00200009327300008250400000","4117918077032230000","0163103100000010"); // EGA
  } else if (trainConfig == 2271){  // std cut
    cuts.AddCutPCMCalo("0008d113","0020000932730000b250400000","4117918077032230000","0163103100000010"); // qT<0.125pT (1D)
    cuts.AddCutPCMCalo("0008d113","0020000932730000a250400000","4117918077032230000","0163103100000010"); // qT<0.125pT (2D) alpha<0.95
    cuts.AddCutPCMCalo("0008d113","0020000932730000a259400000","4117918077032230000","0163103100000010"); // qT<0.125pT (2D) alpha<0.99
    cuts.AddCutPCMCalo("0008d113","0020000932730000a25a400000","4117918077032230000","0163103100000010"); // qT<0.125pT (2D) alpha<1
  } else if (trainConfig == 2272){  // std cut
    cuts.AddCutPCMCalo("0008d113","0020000932730000acc0400000","4117918077032230000","0163103100000010"); // PsiPair<0.15 + chi2<40 + qT<0.125pT (2D) alpha<0.95
    cuts.AddCutPCMCalo("0008d113","0020000932730000a180400000","4117918077032230000","0163103100000010"); // PsiPair<0.20 + chi2<50 + qT<0.125pT (2D) alpha<0.95
    cuts.AddCutPCMCalo("0008d113","0020000932730000acca400000","4117918077032230000","0163103100000010"); // PsiPair<0.15 + chi2<40 + qT<0.125pT (2D) alpha<1
    cuts.AddCutPCMCalo("0008d113","0020000932730000a18a400000","4117918077032230000","0163103100000010"); // PsiPair<0.20 + chi2<50 + qT<0.125pT (2D) alpha<1
  // *********************************************************************************************************
  // 5 TeV 2017 pp - Jet configurations
  // *********************************************************************************************************

  } else if ( trainConfig == 2300){ // Jet analysis pp 5 TeV 2017 EMCAL+DCal
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","411790007l032230000","2l63103100000010"); // INT7 - NO NL
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","411790607l032230000","2l63103100000010"); // INT7 - TB NL
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","411791107l032230000","2l63103100000010"); // Standard EDC
  } else if ( trainConfig == 2301){ // Jet QA
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","411791107l032230000","3l63103100000010"); //

  // TB NL testconfigs diff aggregation thresholds
  } else if ( trainConfig == 2400){ // LHC12 100 MeV aggregation TB NL
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0400000","111110106f032230000","0163103100000010"); // std
  } else if ( trainConfig == 2401){ // LHC12 100 MeV aggregation TB NL
    cuts.AddCutPCMCalo("00052113","00200009f9730000dge0400000","111110106f032230000","0163103100000010"); // std
    cuts.AddCutPCMCalo("00081113","00200009f9730000dge0400000","111110106f032230000","0163103100000010"); // std
  } else if ( trainConfig == 2402){ // LHC12 50 MeV aggregation TB NL
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0400000","111110206f032230000","0163103100000010"); // std
  } else if ( trainConfig == 2403){ // LHC12 50 MeV aggregation TB NL
    cuts.AddCutPCMCalo("00052113","00200009f9730000dge0400000","111110206f032230000","0163103100000010"); // std
    cuts.AddCutPCMCalo("00081113","00200009f9730000dge0400000","111110206f032230000","0163103100000010"); // std
  } else if ( trainConfig == 2404){ // LHC12 150 MeV aggregation TB NL
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0400000","111110306f032230000","0163103100000010"); // std
  } else if ( trainConfig == 2405){ // LHC12 150 MeV aggregation TB NL
    cuts.AddCutPCMCalo("00052113","00200009f9730000dge0400000","111110306f032230000","0163103100000010"); // std
    cuts.AddCutPCMCalo("00081113","00200009f9730000dge0400000","111110306f032230000","0163103100000010"); // std
  } else if ( trainConfig == 2406){ // LHC12 300 MeV aggregation TB NL
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0400000","111110406f032230000","0163103100000010"); // std
  } else if ( trainConfig == 2407){ // LHC12 300 MeV aggregation TB NL
    cuts.AddCutPCMCalo("00052113","00200009f9730000dge0400000","111110406f032230000","0163103100000010"); // std
    cuts.AddCutPCMCalo("00081113","00200009f9730000dge0400000","111110406f032230000","0163103100000010"); // std

  } else if (trainConfig == 2450){  // TB parametrization from Nico on Martin 100MeV points + FineTune
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0400000","411790207fg32230000","0163103100000010"); // INT7
  } else if (trainConfig == 2451){  // TB parametrization from Nico on Martin 100MeV points + FineTune
    cuts.AddCutPCMCalo("0008e113","00200009f9730000dge0400000","411790207fg32230000","0163103100000010"); // EG2
    cuts.AddCutPCMCalo("0008d113","00200009f9730000dge0400000","411790207fg32230000","0163103100000010"); // EG1
  } else if (trainConfig == 2452){  // no NL
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0400000","411790407fg32230000","0163103100000010"); // INT7
  } else if (trainConfig == 2453){  // no NL
    cuts.AddCutPCMCalo("0008e113","00200009f9730000dge0400000","411790407fg32230000","0163103100000010"); // EG2
    cuts.AddCutPCMCalo("0008d113","00200009f9730000dge0400000","411790407fg32230000","0163103100000010"); // EG1
  } else if (trainConfig == 2454){  // 150 MeV aggregation TB NL tests
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0400000","411790507fg32230000","0163103100000010"); // INT7
  } else if (trainConfig == 2455){  // 150 MeV aggregation TB NL tests
    cuts.AddCutPCMCalo("0008e113","00200009f9730000dge0400000","411790507fg32230000","0163103100000010"); // EG2
    cuts.AddCutPCMCalo("0008d113","00200009f9730000dge0400000","411790507fg32230000","0163103100000010"); // EG1
  } else if (trainConfig == 2456){  // kPi0MCv3 for MC and kTestBeamv4 for data
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0400000","411790607fg32230000","0163103100000010"); // INT7
  } else if (trainConfig == 2457){  // kPi0MCv3 for MC and kTestBeamv4 for data
    cuts.AddCutPCMCalo("0008e113","00200009f9730000dge0400000","411790607fg32230000","0163103100000010"); // EG2
    cuts.AddCutPCMCalo("0008d113","00200009f9730000dge0400000","411790607fg32230000","0163103100000010"); // EG1
  } else if (trainConfig == 2458){  // no NL
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0400000","411790007fg32230000","0163103100000010"); // INT7
  } else if (trainConfig == 2459){  // no NL
    cuts.AddCutPCMCalo("0008e113","00200009f9730000dge0400000","411790007fg32230000","0163103100000010"); // EG2
    cuts.AddCutPCMCalo("0008d113","00200009f9730000dge0400000","411790007fg32230000","0163103100000010"); // EG1
  } else if (trainConfig == 2460){  // no E-Scale NL
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0400000","411790407fg32230000","0163103100000010"); // INT7
  } else if (trainConfig == 2461){  // no E-Scale NL
    cuts.AddCutPCMCalo("0008e113","00200009f9730000dge0400000","411790407fg32230000","0163103100000010"); // EG2
    cuts.AddCutPCMCalo("0008d113","00200009f9730000dge0400000","411790407fg32230000","0163103100000010"); // EG1
  } else if (trainConfig == 2462){  // v1 Evi (used in E-calib)
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0400000","411790907fg32230000","0163103100000010"); // INT7
  } else if (trainConfig == 2463){  // v1 Evi (used in E-calib)
    cuts.AddCutPCMCalo("0008e113","00200009f9730000dge0400000","411790907fg32230000","0163103100000010"); // EG2
    cuts.AddCutPCMCalo("0008d113","00200009f9730000dge0400000","411790907fg32230000","0163103100000010"); // EG1


  //************************************************ PCM- EDC analysis 5 TeV pp INT7 sys *********************************
  } else if (trainConfig == 2510) { // PCM variations
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0400000","411793107f032230000","0163103100000010"); //New standard cut for 8TeV analysis for RpA
    cuts.AddCutPCMCalo("00010013","00200009f9730000dge0400000","411793107f032230000","0163103100000010"); // no SPD pileup cut
    cuts.AddCutPCMCalo("00010113","00100009f9730000dge0400000","411793107f032230000","0163103100000010"); // R cut 2.8 -180 cm
    cuts.AddCutPCMCalo("00010113","00500009f9730000dge0400000","411793107f032230000","0163103100000010"); // R cut 10. -180 cm
  } else if (trainConfig == 2511) {
    cuts.AddCutPCMCalo("00010113","00200069f9730000dge0400000","411793107f032230000","0163103100000010"); // min pT 40 MeV
    cuts.AddCutPCMCalo("00010113","00200049f9730000dge0400000","411793107f032230000","0163103100000010"); // min pT 75 MeV
    cuts.AddCutPCMCalo("00010113","00200019f9730000dge0400000","411793107f032230000","0163103100000010"); // min pT 100MeV
  } else if (trainConfig == 2512) {
    cuts.AddCutPCMCalo("00010113","00200068f9730000dge0400000","411793107f032230000","0163103100000010"); // TPC cluster 35%
    cuts.AddCutPCMCalo("00010113","00200066f9730000dge0400000","411793107f032230000","0163103100000010"); // TPC cluster 70%
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0600000","411793107f032230000","0163103100000010"); // cosPA 0.9
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0300000","411793107f032230000","0163103100000010"); // cosPA 0.75
  } else if (trainConfig == 2513) {
    cuts.AddCutPCMCalo("00010113","0020000939730000dge0400000","411793107f032230000","0163103100000010"); // nsig electron   -4,5
    cuts.AddCutPCMCalo("00010113","0020000969730000dge0400000","411793107f032230000","0163103100000010"); // nsig electron -2.5,4
    cuts.AddCutPCMCalo("00010113","00200009f5730000dge0400000","411793107f032230000","0163103100000010"); // nsig pion 2,-10
    cuts.AddCutPCMCalo("00010113","00200009f1730000dge0400000","411793107f032230000","0163103100000010"); // nsig pion 0,-10
  } else if (trainConfig == 2514) {
    cuts.AddCutPCMCalo("00010113","00200009f9030000dge0400000","411793107f032230000","0163103100000010"); // pion nsig min mom 0.50 GeV/c
    cuts.AddCutPCMCalo("00010113","00200009f9630000dge0400000","411793107f032230000","0163103100000010"); // pion nsig min mom 0.25 GeV/c
    cuts.AddCutPCMCalo("00010113","00200009f9760000dge0400000","411793107f032230000","0163103100000010"); // pion nsig max mom 2.00 GeV/c
    cuts.AddCutPCMCalo("00010113","00200009f9710000dge0400000","411793107f032230000","0163103100000010"); // pion nsig max mom 5.00 GeV/c
  } else if (trainConfig == 2515) {
    cuts.AddCutPCMCalo("00010113","00200009f9730000age0400000","411793107f032230000","0163103100000010"); // qT max 0.040, qt pt max 0.11
    cuts.AddCutPCMCalo("00010113","00200009f9730000ege0400000","411793107f032230000","0163103100000010"); // qT max 0.060, qt pt max 0.14
    cuts.AddCutPCMCalo("00010113","00200009f9730000fge0400000","411793107f032230000","0163103100000010"); // qT max 0.070, qt pt max 0.16
  } else if (trainConfig == 2516) {
    cuts.AddCutPCMCalo("00010113","00200009f9730000d1e0400000","411793107f032230000","0163103100000010"); // chi2 50 no chi2 dep.
    cuts.AddCutPCMCalo("00010113","00200009f9730000dfe0400000","411793107f032230000","0163103100000010"); // chi2 50 chi2 dep -0.065
    cuts.AddCutPCMCalo("00010113","00200009f9730000dhe0400000","411793107f032230000","0163103100000010"); // chi2 50 chi2 dep -0.050
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0404000","411793107f032230000","0163103100000010"); // reject close v0
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0406000","411793107f032230000","0163103100000010"); // double count with open angle 0.04
  } else if (trainConfig == 2517) {
    cuts.AddCutPCMCalo("00010113","00200009f9730000dgd0400000","411793107f032230000","0163103100000010"); // Psi pair 0.15 dep
    cuts.AddCutPCMCalo("00010113","00200009f9730000dgf0400000","411793107f032230000","0163103100000010"); // Psi pair 0.20 dep
    cuts.AddCutPCMCalo("00010113","00200009f9730000dgg0400000","411793107f032230000","0163103100000010"); // Psi pair 0.30 dep
  } else if (trainConfig == 2518) {
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0400000","411793107f032230000","0263103100000010"); // variation BG scheme track mult
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0400000","411793107f032230000","0163107100000010"); // alpha meson 0.85
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0400000","411793107f032230000","0163105100000010"); // alpha meson 0.75
    cuts.AddCutPCMCalo("00010113","00200009227300008250404000","411793107f032230000","0163103100000010"); // old cuts (run1)
  } else if (trainConfig == 2519) { // CALO variations
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0400000","411793107f022230000","0163103100000010"); // 0.6 GeV/c
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0400000","411793107f042230000","0163103100000010"); // 0.8 GeV/c
  } else if (trainConfig == 2520) {
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0400000","411793107f031230000","0163103100000010"); // n cells >= 1
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0400000","411793107f033230000","0163103100000010"); // n cells >= 3
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0400000","411793107f032200000","0163103100000010"); // no max M02 cut
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0400000","411793107f032250000","0163103100000010"); // M02 < 0.3
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0400000","411793107f0322k0000","0163103100000010"); // M02, pT-dep
  } else if (trainConfig == 2521) {
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0400000","411793107e032230000","0163103100000010"); // TM var e/p 2.0
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0400000","411793107g032230000","0163103100000010"); // TM var e/p 1.5
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0400000","411793107h032230000","0163103100000010"); // TM var e/p 1.25
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0400000","4117931077032230000","0163103100000010"); // TM var no veto
  } else if (trainConfig == 2522) {
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0400000","411793207f032230000","0163103100000010"); // NL 32
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0400000","411793307f032230000","0163103100000010"); // NL 33
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0400000","411793407f032230000","0163103100000010"); // NL 34
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0400000","411793807f032230000","0163103100000010"); // NL 38
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0400000","411793907f032230000","0163103100000010"); // NL 39
  } else if (trainConfig == 2523) {
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0400000","411793109f032230000","0163103100000010"); // 20/25ns
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0400000","411793104f032230000","0163103100000010"); // 100ns
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0400000","411793105f032230000","0163103100000010"); // 50ns
  } else if (trainConfig == 2524) { // NCell effi
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0400000","411793207f03h230000","0163103100000010"); //
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0400000","411793207f03i230000","0163103100000010"); //
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0400000","411793207f03j230000","0163103100000010"); //
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0400000","411793207f03k230000","0163103100000010"); //
  } else if (trainConfig == 2525) { // BG variations (rotation/swapping)
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0400000","411793207f032230000","0163103100000010"); //
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0400000","411793207f032230000","0r63103100000010"); //
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0400000","411793207f032230000","0s63103100000010"); //
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0400000","411793207f032230000","0t63103100000010"); //
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0400000","411793207f032230000","0u63103100000010"); //
  } else if (trainConfig == 2526) { // BG variations (rotation/swapping)
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0400000","411793207f032230000","0163103100b00010"); //
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0400000","411793207f032230000","0r63103100b00010"); //
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0400000","411793207f032230000","0s63103100b00010"); //
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0400000","411793207f032230000","0t63103100b00010"); //
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0400000","411793207f032230000","0u63103100b00010"); //
  } else if (trainConfig == 2527) { // BG variations (rotation/swapping)
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0400000","388553207f032230000","0163103100000010"); //
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0400000","388553207f032230000","0r63103100000010"); //
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0400000","388553207f032230000","0s63103100000010"); //
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0400000","388553207f032230000","0t63103100000010"); //
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0400000","388553207f032230000","0u63103100000010"); //
  } else if (trainConfig == 2528) { // BG variations (rotation/swapping)
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0400000","388553207f032230000","0163103100b00010"); //
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0400000","388553207f032230000","0r63103100b00010"); //
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0400000","388553207f032230000","0s63103100b00010"); //
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0400000","388553207f032230000","0t63103100b00010"); //
    cuts.AddCutPCMCalo("00010113","00200009f9730000dge0400000","388553207f032230000","0u63103100b00010"); //

    //*************************************************************************************************
    // 13 TeV PCM-PHOS - Systematics
    //*************************************************************************************************
    // Variations for systematics
  } else if (trainConfig == 2600){ // PHOS clusters standard cuts, INT7, NL, std 19
    //                                                               ||
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","24466000sa01cc00000","0163103100b00010"); // INT7 NoNL
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","24466110sa01cc00000","0163103100b00010"); // INT7 NL11
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","24466120sa01cc00000","0163103100b00010"); // INT7 NL12
    //cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163103100b00010"); // INT7 NL19
  } else if (trainConfig == 2601){ // PHOS clusters standard cuts, PHI7, NL, std 19
    //                     ||                                        ||
    cuts.AddCutPCMCalo("00062113","0dm00009f9730000dge0404000","24466000sa01cc00000","0163103100000010"); // PHI7  No NL
    cuts.AddCutPCMCalo("00062113","0dm00009f9730000dge0404000","24466110sa01cc00000","0163103100000010"); // PHI7  NL11
    cuts.AddCutPCMCalo("00062113","0dm00009f9730000dge0404000","24466120sa01cc00000","0163103100000010"); // PHI7  NL12
    cuts.AddCutPCMCalo("00062113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163103100000010"); // PHI7  NL19
    //----------------------------------------------------------------------------------------------------------------------------------------
    // Variations of PHOS Part
    //Standard: "24466190sa01cc00000"
  } else if (trainConfig == 2620){ // timing Cut variation  std -30+30ns
    //                                                                  |
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","244661901a01cc00000","0163103100b00010"); //1:     -1000  +1000 ns
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","244661905a01cc00000","0163103100b00010"); //5:     -50    +50   ns
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","244661907a01cc00000","0163103100b00010"); //7:     -30    +30   ns
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","244661909a01cc00000","0163103100b00010"); //9:     -20    +25   ns
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","24466190aa01cc00000","0163103100b00010"); //a:     -12.5  +13   ns
  } else if (trainConfig == 2621){ // Timing Efficiency Variations std s == 30ns
    //                                                                  |
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","24466190ra01cc00000","0163103100b00010"); //r:     LowPt from MB, 30ns
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","24466190ta01cc00000","0163103100b00010"); //t:     25ns
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","24466190ua01cc00000","0163103100b00010"); //u:     50ns
  } else if (trainConfig == 2622){ // track matching variation std a == pt dependent
    //                                                                   |
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","24466190s001cc00000","0163103100b00010"); //0
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","24466190s101cc00000","0163103100b00010"); //1
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","24466190s401cc00000","0163103100b00010"); //4
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","24466190s501cc00000","0163103100b00010"); //5
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","24466190s601cc00000","0163103100b00010"); //6
  } else if (trainConfig == 2623){ // min energy variation std 0.3 GeV/c
    //                                                                     |
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa00cc00000","0163103100b00010"); //0:     off
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa09cc00000","0163103100b00010"); //9:     0.1 GeV/c
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa02cc00000","0163103100b00010"); //2:     0.5 GeV/c
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa03cc00000","0163103100b00010"); //3:     0.6 GeV/c
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa04cc00000","0163103100b00010"); //4:     0.7 GeV/c
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa05cc00000","0163103100b00010"); //5:     0.8 GeV/c
  } else if (trainConfig == 2624){ // min nCells & M02 variation, std cc
    // std: min nCells = 2 >1GeV; M02 max=100, min=0.1
    //                                                                      ||
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa011000000","0163103100b00010"); //100:   min nCells = 1, minM02 off
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa012200000","0163103100b00010"); //220:   min nCells = 2, all E
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa013200000","0163103100b00010"); //320:   min nCells = 3, all E
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa01dc00000","0163103100b00010"); //dc0:   min nCells = 3, E>1GeV; minM02==0.1 off for E<1GeV
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa01cd00000","0163103100b00010"); //cd0:   min nCells = 2, E>1GeV; minM02==0.2 off for E<1GeV
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa01cc70000","0163103100b00010"); //cc7:   maxM02 == 1.3
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa01cc80000","0163103100b00010"); //cc8:   maxM02 == 2.5
  } else if (trainConfig == 2625){ // reconstructed conversion, std 0==off
    //                                                                          |
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa01cc00100","0163103100b00010"); //1:   rec. conv. 0.02
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa01cc00200","0163103100b00010"); //2:   rec. conv. 0.025
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa01cc00300","0163103100b00010"); //3:   rec. conv. 0.03
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa01cc00400","0163103100b00010"); //4:   rec. conv. 0.035
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa01cc00500","0163103100b00010"); //5:   rec. conv. 0.04
  //----------------------------------------------------------------------------------------------------------------------------------------
  // Meson variations
  } else if ( trainConfig == 2640){//rapidity, std 1 == <0.8
    //                                                                                    |
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163403100b00010"); //4: rapidity variation  y<0.5
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163803100b00010"); //8: rapidity variation  y<0.25
  } else if ( trainConfig == 2641){//alpha, std 3 == <=1.0
    //                                                                                      |
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163106100b00010"); //6: alpha meson variation 1 0<alpha<0.8
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163105100b00010"); //5: alpha meson variation 2 0<alpha<0.75
  } else if ( trainConfig == 2642){ //opening angle, std 1 == >0.5, 1 cell diagonal
    //                                                                                              |
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163103100b00000"); //0: min opening angle 0    -> open
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163103100b00030"); //3: min opening angle 0.01 -> 2 cell diag
  } else if ( trainConfig == 2643){  //Smearing Check
    //                                                                                          |
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163103100400010"); //4
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163103100500010"); //5
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163103100a00010"); //a
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163103100b00010"); //b
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163103100c00010"); //c
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163103100d00010"); //d
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163103100e00010"); //e
  //----------------------------------------------------------------------------------------------------------------------------------------
  //Event Cut Variations, std 00010113
  } else if ( trainConfig == 2650){  //PHOS Sphericiry Check
    //                  |||
    cuts.AddCutPCMCalo("h0510113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163103100b00010"); //h05:  0.    - 0.5
    cuts.AddCutPCMCalo("h5a10113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163103100b00010"); //h5a:  0.5    - 1.
    cuts.AddCutPCMCalo("h0a10113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163103100b00010"); //h0a:  0.    - 1.
  } else if ( trainConfig == 2651){  //PHOS Mult Check
    //                  |||
    cuts.AddCutPCMCalo("n0110113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163103100b00010"); //n01 INT7 0-10%
    cuts.AddCutPCMCalo("n1210113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163103100b00010"); //n12 INT7 10-20%
    cuts.AddCutPCMCalo("n2510113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163103100b00010"); //n25 INT7 20-50%
    cuts.AddCutPCMCalo("n5a10113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163103100b00010"); //n5a INT7 50-100%
    cuts.AddCutPCMCalo("m0110113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163103100b00010"); //m01 INT7
    cuts.AddCutPCMCalo("m1510113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163103100b00010"); //m15 INT7
    cuts.AddCutPCMCalo("m5a10113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163103100b00010"); //m5a INT7
  } else if ( trainConfig == 2652){ // SPD mult selection PCMPHOS
    //                  |||
    cuts.AddCutPCMCalo("o0110113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163103100b00010"); //o01: 0-1%
    cuts.AddCutPCMCalo("o0210113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163103100b00010"); //o02: 0-2%
    cuts.AddCutPCMCalo("o0510113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163103100b00010"); //o05: 0-5%
    cuts.AddCutPCMCalo("o5k10113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163103100b00010"); //o5k: 5-20%
    cuts.AddCutPCMCalo("p2610113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163103100b00010"); //p26: 20-60%
    cuts.AddCutPCMCalo("p6a10113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163103100b00010"); //p6a: 60-100%
  //----------------------------------------------------------------------------------------------------------------------------------------
  // Variations of PCM Part  // std 0dm00009f9730000dge0404000
  } else if (trainConfig == 2660) { //single pt, std 0 == 0.05 GeV/c
    //                                    |
    cuts.AddCutPCMCalo("00010113", "0dm00069f9730000dge0404000","24466190sa01cc00000", "0163103100b00010"); //6: min pT 40 MeV
    cuts.AddCutPCMCalo("00010113", "0dm00049f9730000dge0404000","24466190sa01cc00000", "0163103100b00010"); //4: min pT 75 MeV
    cuts.AddCutPCMCalo("00010113", "0dm00019f9730000dge0404000","24466190sa01cc00000", "0163103100b00010"); //1: min pT 100MeV
  } else if (trainConfig == 2661) { // min TPC clusters / findable clusters, std 9 == 60%
    //                                     |
    cuts.AddCutPCMCalo("00010113", "0dm00008f9730000dge0404000","24466190sa01cc00000", "0163103100b00010"); //8: TPC cluster 35%
    cuts.AddCutPCMCalo("00010113", "0dm00006f9730000dge0404000","24466190sa01cc00000", "0163103100b00010"); //6: TPC cluster 70%
  } else if (trainConfig == 2662) { //Cosine Pointing Angle, std 4 == 0.85
    //                                                  |
    cuts.AddCutPCMCalo("00010113", "0dm00009f9730000dge0604000","24466190sa01cc00000", "0163103100b00010"); //6: cosPA 0.9
    cuts.AddCutPCMCalo("00010113", "0dm00009f9730000dge0304000","24466190sa01cc00000", "0163103100b00010"); //3: cosPA 0.75
  } else if (trainConfig == 2663) { // nsigma electron; std f == -3, 4
    //                                      |
    cuts.AddCutPCMCalo("00010113", "0dm0000939730000dge0404000","24466190sa01cc00000", "0163103100b00010"); //3: nsig electron   -4, 5
    cuts.AddCutPCMCalo("00010113", "0dm0000969730000dge0404000","24466190sa01cc00000", "0163103100b00010"); //6: nsig electron -2.5, 4
  } else if (trainConfig == 2664) { // nsigma pion; std 9 == low pt 1, highpt 0.5
    //                                       |
    cuts.AddCutPCMCalo("00010113", "0dm00009f5730000dge0404000","24466190sa01cc00000", "0163103100b00010"); //5: nsig pion 2,-10
    cuts.AddCutPCMCalo("00010113", "0dm00009f2730000dge0404000","24466190sa01cc00000", "0163103100b00010"); //2: nsig pion 1,-10
    cuts.AddCutPCMCalo("00010113", "0dm00009f1730000dge0404000","24466190sa01cc00000", "0163103100b00010"); //1: nsig pion 0,-10
  } else if (trainConfig == 2665) {// pion nsig min mom, std 7 == 0.4 GeV/c
    //                                        |
    cuts.AddCutPCMCalo("00010113", "0dm00009f9030000dge0404000","24466190sa01cc00000", "0163103100b00010"); //0: pion nsig min mom 0.50 GeV/c
    cuts.AddCutPCMCalo("00010113", "0dm00009f9630000dge0404000","24466190sa01cc00000", "0163103100b00010"); //6: pion nsig min mom 0.25 GeV/c
  } else if (trainConfig == 2666) {// pion nsig max mom, std 3 == 20 GeV/c
    //                                         |
    cuts.AddCutPCMCalo("00010113", "0dm00009f9760000dge0404000","24466190sa01cc00000", "0163103100b00010"); //6: pion nsig max mom 2.00 GeV/c
    cuts.AddCutPCMCalo("00010113", "0dm00009f9710000dge0404000","24466190sa01cc00000", "0163103100b00010"); //1: pion nsig max mom 5.00 GeV/c
    cuts.AddCutPCMCalo("00010113", "0dm00009f9700000dge0404000","24466190sa01cc00000", "0163103100b00010"); //0: pion nsig max mom 100.00 GeV/c
  } else if (trainConfig == 2667){ //qT max, std d == qTPtMax=0.125, qTMax=0.050, fDo2DQt=kTRUE, fDoQtGammaSelection=2
    //                                              |
    cuts.AddCutPCMCalo("00010113", "0dm00009f97300008ge0404000","24466190sa01cc00000", "0163103100b00010"); //8: qT max 0.05 1D, fDo2DQt=kTRUE
    cuts.AddCutPCMCalo("00010113", "0dm00009f97300003ge0404000","24466190sa01cc00000", "0163103100b00010"); //3: qT max 0.05 1D, fDo2DQt=kFALSE
    cuts.AddCutPCMCalo("00010113", "0dm00009f97300002ge0404000","24466190sa01cc00000", "0163103100b00010"); //2: qT max 0.06 2D, fDo2DQt=kTRUE
    cuts.AddCutPCMCalo("00010113", "0dm00009f97300009ge0404000","24466190sa01cc00000", "0163103100b00010"); //9: qT max 0.03 2D, fDo2DQt=kTRUE
  } else if (trainConfig == 2668) { // Psi pair 1D, std e == 0.18, fDo2DPsiPairChi2 = 2
    //                                                |
    cuts.AddCutPCMCalo("00010113", "0dm00009f9730000dg50404000","24466190sa01cc00000", "0163103100b00010"); //5: Psi pair 0.1  1D
    cuts.AddCutPCMCalo("00010113", "0dm00009f9730000dg10404000","24466190sa01cc00000", "0163103100b00010"); //1: Psi pair 0.1  1D
    cuts.AddCutPCMCalo("00010113", "0dm00009f9730000dg60404000","24466190sa01cc00000", "0163103100b00010"); //6: Psi pair 0.05 2D
    cuts.AddCutPCMCalo("00010113", "0dm00009f9730000dg80404000","24466190sa01cc00000", "0163103100b00010"); //8: Psi pair 0.2  2D
  } else if (trainConfig == 2669){ // qT 2D pT dep, std qt d == qTPtMax=0.125, qTMax=0.050, fDo2DQt=kTRUE, fDoQtGammaSelection=2, std 0 == off
    //                                              |  |
    cuts.AddCutPCMCalo("00010113", "0dm00009f9730000dge0404000","24466190sa01cc00000", "0163103100b00010"); //d,0: qT<0.110pT (2D) alpha<0.99
    cuts.AddCutPCMCalo("00010113", "0dm00009f9730000c259404000","24466190sa01cc00000", "0163103100b00010"); //c,9: qT<0.110pT (2D) alpha<0.99
    cuts.AddCutPCMCalo("00010113", "0dm00009f9730000a259404000","24466190sa01cc00000", "0163103100b00010"); //a,9: qT<0.125pT (2D) alpha<0.99
    cuts.AddCutPCMCalo("00010113", "0dm00009f9730000e259404000","24466190sa01cc00000", "0163103100b00010"); //e,9: qT<0.130pT (2D) alpha<0.99
  } else if (trainConfig == 2670){ // chi2, PsiPair cont., std 2 ==
    //                                               |
    cuts.AddCutPCMCalo("00010113", "0dm00009f97300008fd0404000","24466190sa01cc00000", "0163103100b00010"); //f: PsiPair<0.15exp(-0.065chi2)
    cuts.AddCutPCMCalo("00010113", "0dm00009f97300008ge0404000","24466190sa01cc00000", "0163103100b00010"); //g: PsiPair<0.18exp(-0.055chi2)
    cuts.AddCutPCMCalo("00010113", "0dm00009f97300008hf0404000","24466190sa01cc00000", "0163103100b00010"); //h: PsiPair<0.20exp(-0.050chi2)
  //----------------------------------------------------------------------------------------------------------------------------------------
    //*************************************************************************************************
    // 13 TeV PCM-PHOS - Systematics - Trigger
    //*************************************************************************************************
    // Variations of PHOS Part
    //Standard: "24466190sa01cc00000"
  } else if (trainConfig == 2720){ // timing Cut variation  std -30+30ns
    //                                                                  |
    cuts.AddCutPCMCalo("00062113","0dm00009f9730000dge0404000","244661901a01cc00000","0163103100b00010"); //1:     -1000  +1000 ns
    cuts.AddCutPCMCalo("00062113","0dm00009f9730000dge0404000","244661905a01cc00000","0163103100b00010"); //5:     -50    +50   ns
    cuts.AddCutPCMCalo("00062113","0dm00009f9730000dge0404000","244661907a01cc00000","0163103100b00010"); //7:     -30    +30   ns
    cuts.AddCutPCMCalo("00062113","0dm00009f9730000dge0404000","244661909a01cc00000","0163103100b00010"); //9:     -20    +25   ns
    cuts.AddCutPCMCalo("00062113","0dm00009f9730000dge0404000","24466190aa01cc00000","0163103100b00010"); //a:     -12.5  +13   ns
  } else if (trainConfig == 2721){ // Timing Efficiency Variations std s == 30ns
    //                                                                  |
    cuts.AddCutPCMCalo("00062113","0dm00009f9730000dge0404000","24466190ra01cc00000","0163103100b00010"); //r:     LowPt from MB, 30ns
    cuts.AddCutPCMCalo("00062113","0dm00009f9730000dge0404000","24466190ta01cc00000","0163103100b00010"); //t:     25ns
    cuts.AddCutPCMCalo("00062113","0dm00009f9730000dge0404000","24466190ua01cc00000","0163103100b00010"); //u:     50ns
  } else if (trainConfig == 2722){ // track matching variation std a == pt dependent
    //                                                                   |
    cuts.AddCutPCMCalo("00062113","0dm00009f9730000dge0404000","24466190s001cc00000","0163103100b00010"); //0
    cuts.AddCutPCMCalo("00062113","0dm00009f9730000dge0404000","24466190s101cc00000","0163103100b00010"); //1
    cuts.AddCutPCMCalo("00062113","0dm00009f9730000dge0404000","24466190s401cc00000","0163103100b00010"); //4
    cuts.AddCutPCMCalo("00062113","0dm00009f9730000dge0404000","24466190s501cc00000","0163103100b00010"); //5
    cuts.AddCutPCMCalo("00062113","0dm00009f9730000dge0404000","24466190s601cc00000","0163103100b00010"); //6
  } else if (trainConfig == 2723){ // min energy variation std 0.3 GeV/c
    //                                                                     |
    cuts.AddCutPCMCalo("00062113","0dm00009f9730000dge0404000","24466190sa00cc00000","0163103100b00010"); //0:     off
    cuts.AddCutPCMCalo("00062113","0dm00009f9730000dge0404000","24466190sa09cc00000","0163103100b00010"); //9:     0.1 GeV/c
    cuts.AddCutPCMCalo("00062113","0dm00009f9730000dge0404000","24466190sa02cc00000","0163103100b00010"); //2:     0.5 GeV/c
    cuts.AddCutPCMCalo("00062113","0dm00009f9730000dge0404000","24466190sa03cc00000","0163103100b00010"); //3:     0.6 GeV/c
    cuts.AddCutPCMCalo("00062113","0dm00009f9730000dge0404000","24466190sa04cc00000","0163103100b00010"); //4:     0.7 GeV/c
    cuts.AddCutPCMCalo("00062113","0dm00009f9730000dge0404000","24466190sa05cc00000","0163103100b00010"); //5:     0.8 GeV/c
  } else if (trainConfig == 2724){ // min nCells & M02 variation, std cc
    // std: min nCells = 2 >1GeV; M02 max=100, min=0.1
    //                                                                      ||
    cuts.AddCutPCMCalo("00062113","0dm00009f9730000dge0404000","24466190sa011000000","0163103100b00010"); //100:   min nCells = 1, minM02 off
    cuts.AddCutPCMCalo("00062113","0dm00009f9730000dge0404000","24466190sa012200000","0163103100b00010"); //220:   min nCells = 2, all E
    cuts.AddCutPCMCalo("00062113","0dm00009f9730000dge0404000","24466190sa013200000","0163103100b00010"); //320:   min nCells = 3, all E
    cuts.AddCutPCMCalo("00062113","0dm00009f9730000dge0404000","24466190sa01dc00000","0163103100b00010"); //dc0:   min nCells = 3, E>1GeV; minM02==0.1 off for E<1GeV
    cuts.AddCutPCMCalo("00062113","0dm00009f9730000dge0404000","24466190sa01cd00000","0163103100b00010"); //cd0:   min nCells = 2, E>1GeV; minM02==0.2 off for E<1GeV
    cuts.AddCutPCMCalo("00062113","0dm00009f9730000dge0404000","24466190sa01cc70000","0163103100b00010"); //cc7:   maxM02 == 1.3
    cuts.AddCutPCMCalo("00062113","0dm00009f9730000dge0404000","24466190sa01cc80000","0163103100b00010"); //cc8:   maxM02 == 2.5
  } else if (trainConfig == 2725){ // reconstructed conversion, std 0==off
    //                                                                          |
    cuts.AddCutPCMCalo("00062113","0dm00009f9730000dge0404000","24466190sa01cc00100","0163103100b00010"); //1:   rec. conv. 0.02
    cuts.AddCutPCMCalo("00062113","0dm00009f9730000dge0404000","24466190sa01cc00200","0163103100b00010"); //2:   rec. conv. 0.025
    cuts.AddCutPCMCalo("00062113","0dm00009f9730000dge0404000","24466190sa01cc00300","0163103100b00010"); //3:   rec. conv. 0.03
    cuts.AddCutPCMCalo("00062113","0dm00009f9730000dge0404000","24466190sa01cc00400","0163103100b00010"); //4:   rec. conv. 0.035
    cuts.AddCutPCMCalo("00062113","0dm00009f9730000dge0404000","24466190sa01cc00500","0163103100b00010"); //5:   rec. conv. 0.04
  //----------------------------------------------------------------------------------------------------------------------------------------
  // Meson variations
  } else if ( trainConfig == 2740){//rapidity, std 1 == <0.8
    //                                                                                    |
    cuts.AddCutPCMCalo("00062113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163403100b00010"); //4: rapidity variation  y<0.5
    cuts.AddCutPCMCalo("00062113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163803100b00010"); //8: rapidity variation  y<0.25
  } else if ( trainConfig == 2741){//alpha, std 3 == <=1.0
    //                                                                                      |
    cuts.AddCutPCMCalo("00062113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163106100b00010"); //6: alpha meson variation 1 0<alpha<0.8
    cuts.AddCutPCMCalo("00062113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163105100b00010"); //5: alpha meson variation 2 0<alpha<0.75
  } else if ( trainConfig == 2742){ //opening angle, std 1 == >0.5, 1 cell diagonal
    //                                                                                              |
    cuts.AddCutPCMCalo("00062113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163103100b00000"); //0: min opening angle 0    -> open
    cuts.AddCutPCMCalo("00062113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163103100b00030"); //3: min opening angle 0.01 -> 2 cell diag
  } else if ( trainConfig == 2743){  //Smearing Check
    //                                                                                          |
    cuts.AddCutPCMCalo("00062113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163103100400010"); //4
    cuts.AddCutPCMCalo("00062113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163103100500010"); //5
    cuts.AddCutPCMCalo("00062113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163103100a00010"); //a
    cuts.AddCutPCMCalo("00062113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163103100b00010"); //b
    cuts.AddCutPCMCalo("00062113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163103100c00010"); //c
    cuts.AddCutPCMCalo("00062113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163103100d00010"); //d
    cuts.AddCutPCMCalo("00062113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163103100e00010"); //e
  //----------------------------------------------------------------------------------------------------------------------------------------
  // Variations of PCM Part  // std 0dm00009f9730000dge0404000
  } else if (trainConfig == 2760) { //single pt, std 0 == 0.05 GeV/c
    //                                    |
    cuts.AddCutPCMCalo("00062113", "0dm00069f9730000dge0404000","24466190sa01cc00000", "0163103100b00010"); //6: min pT 40 MeV
    cuts.AddCutPCMCalo("00062113", "0dm00049f9730000dge0404000","24466190sa01cc00000", "0163103100b00010"); //4: min pT 75 MeV
    cuts.AddCutPCMCalo("00062113", "0dm00019f9730000dge0404000","24466190sa01cc00000", "0163103100b00010"); //1: min pT 100MeV
  } else if (trainConfig == 2761) { // min TPC clusters / findable clusters, std 9 == 60%
    //                                     |
    cuts.AddCutPCMCalo("00062113", "0dm00008f9730000dge0404000","24466190sa01cc00000", "0163103100b00010"); //8: TPC cluster 35%
    cuts.AddCutPCMCalo("00062113", "0dm00006f9730000dge0404000","24466190sa01cc00000", "0163103100b00010"); //6: TPC cluster 70%
  } else if (trainConfig == 2762) { //Cosine Pointing Angle, std 4 == 0.85
    //                                                  |
    cuts.AddCutPCMCalo("00062113", "0dm00009f9730000dge0604000","24466190sa01cc00000", "0163103100b00010"); //6: cosPA 0.9
    cuts.AddCutPCMCalo("00062113", "0dm00009f9730000dge0304000","24466190sa01cc00000", "0163103100b00010"); //3: cosPA 0.75
  } else if (trainConfig == 2763) { // nsigma electron; std f == -3, 4
    //                                      |
    cuts.AddCutPCMCalo("00062113", "0dm0000939730000dge0404000","24466190sa01cc00000", "0163103100b00010"); //3: nsig electron   -4, 5
    cuts.AddCutPCMCalo("00062113", "0dm0000969730000dge0404000","24466190sa01cc00000", "0163103100b00010"); //6: nsig electron -2.5, 4
  } else if (trainConfig == 2764) { // nsigma pion; std 9 == low pt 1, highpt 0.5
    //                                       |
    cuts.AddCutPCMCalo("00062113", "0dm00009f5730000dge0404000","24466190sa01cc00000", "0163103100b00010"); //5: nsig pion 2,-10
    cuts.AddCutPCMCalo("00062113", "0dm00009f2730000dge0404000","24466190sa01cc00000", "0163103100b00010"); //2: nsig pion 1,-10
    cuts.AddCutPCMCalo("00062113", "0dm00009f1730000dge0404000","24466190sa01cc00000", "0163103100b00010"); //1: nsig pion 0,-10
  } else if (trainConfig == 2765) {// pion nsig min mom, std 7 == 0.4 GeV/c
    //                                        |
    cuts.AddCutPCMCalo("00062113", "0dm00009f9030000dge0404000","24466190sa01cc00000", "0163103100b00010"); //0: pion nsig min mom 0.50 GeV/c
    cuts.AddCutPCMCalo("00062113", "0dm00009f9630000dge0404000","24466190sa01cc00000", "0163103100b00010"); //6: pion nsig min mom 0.25 GeV/c
  } else if (trainConfig == 2766) {// pion nsig max mom, std 3 == 20 GeV/c
    //                                         |
    cuts.AddCutPCMCalo("00062113", "0dm00009f9760000dge0404000","24466190sa01cc00000", "0163103100b00010"); //6: pion nsig max mom 2.00 GeV/c
    cuts.AddCutPCMCalo("00062113", "0dm00009f9710000dge0404000","24466190sa01cc00000", "0163103100b00010"); //1: pion nsig max mom 5.00 GeV/c
    cuts.AddCutPCMCalo("00062113", "0dm00009f9700000dge0404000","24466190sa01cc00000", "0163103100b00010"); //0: pion nsig max mom 100.00 GeV/c
  } else if (trainConfig == 2767){ //qT max, std d == qTPtMax=0.125, qTMax=0.050, fDo2DQt=kTRUE, fDoQtGammaSelection=2
    //                                              |
    cuts.AddCutPCMCalo("00062113", "0dm00009f97300008ge0404000","24466190sa01cc00000", "0163103100b00010"); //8: qT max 0.05 1D, fDo2DQt=kTRUE
    cuts.AddCutPCMCalo("00062113", "0dm00009f97300003ge0404000","24466190sa01cc00000", "0163103100b00010"); //3: qT max 0.05 1D, fDo2DQt=kFALSE
    cuts.AddCutPCMCalo("00062113", "0dm00009f97300002ge0404000","24466190sa01cc00000", "0163103100b00010"); //2: qT max 0.06 2D, fDo2DQt=kTRUE
    cuts.AddCutPCMCalo("00062113", "0dm00009f97300009ge0404000","24466190sa01cc00000", "0163103100b00010"); //9: qT max 0.03 2D, fDo2DQt=kTRUE
  } else if (trainConfig == 2768) { // Psi pair 1D, std e == 0.18, fDo2DPsiPairChi2 = 2
    //                                                |
    cuts.AddCutPCMCalo("00062113", "0dm00009f9730000dg50404000","24466190sa01cc00000", "0163103100b00010"); //5: Psi pair 0.1  1D
    cuts.AddCutPCMCalo("00062113", "0dm00009f9730000dg10404000","24466190sa01cc00000", "0163103100b00010"); //1: Psi pair 0.1  1D
    cuts.AddCutPCMCalo("00062113", "0dm00009f9730000dg60404000","24466190sa01cc00000", "0163103100b00010"); //6: Psi pair 0.05 2D
    cuts.AddCutPCMCalo("00062113", "0dm00009f9730000dg80404000","24466190sa01cc00000", "0163103100b00010"); //8: Psi pair 0.2  2D
  } else if (trainConfig == 2769){ // qT 2D pT dep, std qt d == qTPtMax=0.125, qTMax=0.050, fDo2DQt=kTRUE, fDoQtGammaSelection=2, std 0 == off
    //                                              |  |
    cuts.AddCutPCMCalo("00062113", "0dm00009f9730000dge0404000","24466190sa01cc00000", "0163103100b00010"); //d,0: qT<0.110pT (2D) alpha<0.99
    cuts.AddCutPCMCalo("00062113", "0dm00009f9730000c259404000","24466190sa01cc00000", "0163103100b00010"); //c,9: qT<0.110pT (2D) alpha<0.99
    cuts.AddCutPCMCalo("00062113", "0dm00009f9730000a259404000","24466190sa01cc00000", "0163103100b00010"); //a,9: qT<0.125pT (2D) alpha<0.99
    cuts.AddCutPCMCalo("00062113", "0dm00009f9730000e259404000","24466190sa01cc00000", "0163103100b00010"); //e,9: qT<0.130pT (2D) alpha<0.99
  } else if (trainConfig == 2770){ // chi2, PsiPair cont., std 2 ==
    //                                               |
    cuts.AddCutPCMCalo("00062113", "0dm00009f97300008fd0404000","24466190sa01cc00000", "0163103100b00010"); //f: PsiPair<0.15exp(-0.065chi2)
    cuts.AddCutPCMCalo("00062113", "0dm00009f97300008ge0404000","24466190sa01cc00000", "0163103100b00010"); //g: PsiPair<0.18exp(-0.055chi2)
    cuts.AddCutPCMCalo("00062113", "0dm00009f97300008hf0404000","24466190sa01cc00000", "0163103100b00010"); //h: PsiPair<0.20exp(-0.050chi2)
  //----------------------------------------------------------------------------------------------------------------------------------------


  } else if (trainConfig == 3010){ // settings Dmitri
    cuts.AddCutPCMCalo("00010103", "0dm00009f9730000dge0404000","24466640ya09dc00000", "0163103100b00010"); // INT7
  } else if (trainConfig == 3011){ // settings Dmitri
    cuts.AddCutPCMCalo("00062103", "0dm00009f9730000dge0404000","24466640ya09dc00000", "0163103100b00010"); // PHI7


    // Variations for systematics   For Trigger
    //EG2
  } else if ( trainConfig ==3100){ // EMCAL clusters
    cuts.AddCutPCMCalo("0008e113","0dm00009f9730000dge0404000","4117901067032230000","0163103100b00010"); // Mixed bck
    cuts.AddCutPCMCalo("0008e113","0dm00009f9730000dge0404000","4117901067032230000","0r63103100b00010"); // 90 degree rotation wo evt. weight
    cuts.AddCutPCMCalo("0008e113","0dm00009f9730000dge0404000","4117901067032230000","0s63103100b00010"); // 90 degree rotation with evt. weight
    cuts.AddCutPCMCalo("0008e113","0dm00009f9730000dge0404000","4117901067032230000","0t63103100b00010"); // random angle with evt. weight
    cuts.AddCutPCMCalo("0008e113","0dm00009f9730000dge0404000","4117901067032230000","0u63103100b00010"); // random angle with multiple decays with evt. weight
  } else if (trainConfig == 3101) {  //   (fPSigSmearing, fPSigSmearingCte)
    cuts.AddCutPCMCalo("0008e113","0dm00009f9730000dge0404000","4117901067032230000","0163103100a00010"); // smearing (0.0275, 0.025)
    cuts.AddCutPCMCalo("0008e113","0dm00009f9730000dge0404000","4117901067032230000","0163103100b00010"); // smearing (0.025,  0.030)
    cuts.AddCutPCMCalo("0008e113","0dm00009f9730000dge0404000","4117901067032230000","0163103100c00010"); // smearing (0.0275, 0.020)
    cuts.AddCutPCMCalo("0008e113","0dm00009f9730000dge0404000","4117901067032230000","0163103100f00010"); // smearing (0.0275, 0.015)
    cuts.AddCutPCMCalo("0008e113","0dm00009f9730000dge0404000","4117901067032230000","0163103100g00010"); // smearing (0.025,  0.020)
    // Variations of EDC Part
  } else if (trainConfig == 3102){ // timing Cut variation  std -50+30ns
    cuts.AddCutPCMCalo("0008e113","0dm00009f9730000dge0404000","4117901017032230000","0163103100b00010"); //     -1000  +1000 ns
    cuts.AddCutPCMCalo("0008e113","0dm00009f9730000dge0404000","4117901077032230000","0163103100b00010"); //     -30    +30   ns
    cuts.AddCutPCMCalo("0008e113","0dm00009f9730000dge0404000","4117901097032230000","0163103100b00010"); //     -20    +25   ns
    cuts.AddCutPCMCalo("0008e113","0dm00009f9730000dge0404000","41179010a7032230000","0163103100b00010"); //     -12.5  +13   ns
  } else if (trainConfig == 3103){ // track matching variation
    cuts.AddCutPCMCalo("0008e113","0dm00009f9730000dge0404000","4117901060032230000","0163103100b00010"); //
    cuts.AddCutPCMCalo("0008e113","0dm00009f9730000dge0404000","4117901061032230000","0163103100b00010"); //
    cuts.AddCutPCMCalo("0008e113","0dm00009f9730000dge0404000","4117901066032230000","0163103100b00010"); //
  } else if (trainConfig == 3104){ // min nCells & M02 variation // std: min nCells = 1
    cuts.AddCutPCMCalo("0008e113","0dm00009f9730000dge0404000","4117901067031230000","0163103100b00010"); //   min nCells = 1
    cuts.AddCutPCMCalo("0008e113","0dm00009f9730000dge0404000","4117901067033230000","0163103100b00010"); //   min nCells = 3
  } else if (trainConfig == 3105){ // min energy variation std 0.7 GeV/c
    cuts.AddCutPCMCalo("0008e113","0dm00009f9730000dge0404000","4117901067002230000","0163103100b00010"); //     0.1 GeV/c
    cuts.AddCutPCMCalo("0008e113","0dm00009f9730000dge0404000","4117901067012230000","0163103100b00010"); //     0.5 GeV/c
    cuts.AddCutPCMCalo("0008e113","0dm00009f9730000dge0404000","4117901067022230000","0163103100b00010"); //     0.6 GeV/c
    cuts.AddCutPCMCalo("0008e113","0dm00009f9730000dge0404000","4117901067042230000","0163103100b00010"); //     0.8 GeV/c
    cuts.AddCutPCMCalo("0008e113","0dm00009f9730000dge0404000","4117901067052230000","0163103100b00010"); //     0.9 GeV/c
  } else if (trainConfig == 3106){ // min nCells & M02 variation  // std: M02 max=0.7, min=0.1
    cuts.AddCutPCMCalo("0008e113","0dm00009f9730000dge0404000","4117901067032210000","0163103100b00010"); //   max M02    = 1
    cuts.AddCutPCMCalo("0008e113","0dm00009f9730000dge0404000","4117901067032240000","0163103100b00010"); //   max M02    = 0.4
  } else if (trainConfig == 3107){ // exotic cluster
    cuts.AddCutPCMCalo("0008e113","0dm00009f9730000dge0404000","4117901067232230000","0163103100b00010"); // ExC 99
    cuts.AddCutPCMCalo("0008e113","0dm00009f9730000dge0404000","4117901067532230000","0163103100b00010"); // ExC 97
    cuts.AddCutPCMCalo("0008e113","0dm00009f9730000dge0404000","4117901067832230000","0163103100b00010"); // ExC 95
    cuts.AddCutPCMCalo("0008e113","0dm00009f9730000dge0404000","4117901067b32230000","0163103100b00010"); // ExC 95 + TCard > 50
    cuts.AddCutPCMCalo("0008e113","0dm00009f9730000dge0404000","4117901067e32230000","0163103100b00010"); // ExC 97 + TCard > 50

    // Variations of PCM Part  // std 0dm00009f9730000dge0404000
  } else if (trainConfig == 3108) {
    cuts.AddCutPCMCalo("0008e113","0dm00069f9730000dge0404000","4117901067032230000","0163103100b00010"); // min pT 40 MeV
    cuts.AddCutPCMCalo("0008e113","0dm00049f9730000dge0404000","4117901067032230000","0163103100b00010"); // min pT 75 MeV
    cuts.AddCutPCMCalo("0008e113","0dm00019f9730000dge0404000","4117901067032230000","0163103100b00010"); // min pT 100MeV
  } else if (trainConfig == 3109) {
    cuts.AddCutPCMCalo("0008e113","0dm00008f9730000dge0404000","4117901067032230000","0163103100b00010"); // TPC cluster 35%
    cuts.AddCutPCMCalo("0008e113","0dm00006f9730000dge0404000","4117901067032230000","0163103100b00010"); // TPC cluster 70%
    cuts.AddCutPCMCalo("0008e113","0dm00009f9730000dge0604000","4117901067032230000","0163103100b00010"); // cosPA 0.9
    cuts.AddCutPCMCalo("0008e113","0dm00009f9730000dge0304000","4117901067032230000","0163103100b00010"); // cosPA 0.75
  } else if (trainConfig == 3110) {
    cuts.AddCutPCMCalo("0008e113","0dm0000939730000dge0404000","4117901067032230000","0163103100b00010"); // nsig electron   -4,5
    cuts.AddCutPCMCalo("0008e113","0dm0000969730000dge0404000","4117901067032230000","0163103100b00010"); // nsig electron -2.5,4
    cuts.AddCutPCMCalo("0008e113","0dm00009f5730000dge0404000","4117901067032230000","0163103100b00010"); // nsig pion 2,-10
    cuts.AddCutPCMCalo("0008e113","0dm00009f1730000dge0404000","4117901067032230000","0163103100b00010"); // nsig pion 0,-10
  } else if (trainConfig == 3111) {
    cuts.AddCutPCMCalo("0008e113","0dm00009f9030000dge0404000","4117901067032230000","0163103100b00010"); // pion nsig min mom 0.50 GeV/c
    cuts.AddCutPCMCalo("0008e113","0dm00009f9630000dge0404000","4117901067032230000","0163103100b00010"); // pion nsig min mom 0.25 GeV/c
    cuts.AddCutPCMCalo("0008e113","0dm00009f9760000dge0404000","4117901067032230000","0163103100b00010"); // pion nsig max mom 2.00 GeV/c
    cuts.AddCutPCMCalo("0008e113","0dm00009f9710000dge0404000","4117901067032230000","0163103100b00010"); // pion nsig max mom 5.00 GeV/c
  } else if (trainConfig ==3112){
    cuts.AddCutPCMCalo("0008e113","0dm00009f97300008ge0404000","4117901067032230000","0163103100b00010"); // qT max 0.05 1D
    cuts.AddCutPCMCalo("0008e113","0dm00009f97300003ge0404000","4117901067032230000","0163103100b00010"); // qT max 0.05 1D
    cuts.AddCutPCMCalo("0008e113","0dm00009f97300002ge0404000","4117901067032230000","0163103100b00010"); // qT max 0.06 2D
    cuts.AddCutPCMCalo("0008e113","0dm00009f97300009ge0404000","4117901067032230000","0163103100b00010"); // qT max 0.03 2D
  } else if (trainConfig == 3113) {
    cuts.AddCutPCMCalo("0008e113","0dm00009f9730000dg50404000","4117901067032230000","0163103100b00010"); // Psi pair 0.1  1D
    cuts.AddCutPCMCalo("0008e113","0dm00009f9730000dg10404000","4117901067032230000","0163103100b00010"); // Psi pair 0.1  1D
    cuts.AddCutPCMCalo("0008e113","0dm00009f9730000dg60404000","4117901067032230000","0163103100b00010"); // Psi pair 0.05 2D
    cuts.AddCutPCMCalo("0008e113","0dm00009f9730000dg80404000","4117901067032230000","0163103100b00010"); // Psi pair 0.2  2D
  } else if ( trainConfig == 3114){ // qT 2D pT dep
    cuts.AddCutPCMCalo("0008e113","0dm00009f9730000dge0404000","4117901067032230000","0163103100b00010"); // qT<0.110pT (2D) alpha<0.99
    cuts.AddCutPCMCalo("0008e113","0dm00009f9730000c259404000","4117901067032230000","0163103100b00010"); // qT<0.110pT (2D) alpha<0.99
    cuts.AddCutPCMCalo("0008e113","0dm00009f9730000a259404000","4117901067032230000","0163103100b00010"); // qT<0.125pT (2D) alpha<0.99
    cuts.AddCutPCMCalo("0008e113","0dm00009f9730000e259404000","4117901067032230000","0163103100b00010"); // qT<0.130pT (2D) alpha<0.99
  } else if ( trainConfig == 3115){ // chi2, PsiPair cont.
    cuts.AddCutPCMCalo("0008e113","0dm00009f97300008fd0404000","4117901067032230000","0163103100b00010"); // PsiPair<0.15exp(-0.065chi2)
    cuts.AddCutPCMCalo("0008e113","0dm00009f97300008ge0404000","4117901067032230000","0163103100b00010"); // PsiPair<0.18exp(-0.055chi2)
    cuts.AddCutPCMCalo("0008e113","0dm00009f97300008hf0404000","4117901067032230000","0163103100b00010"); // PsiPair<0.20exp(-0.050chi2)
  } else if (trainConfig == 3116){ // exotic cluster in triggered data
    cuts.AddCutPCMCalo("0008e113","0dm00009f9730000dge0404000","4117901067232230000","0163103100b00010"); // ExC 99
    cuts.AddCutPCMCalo("0008e113","0dm00009f9730000dge0404000","4117901067532230000","0163103100b00010"); // ExC 97
    cuts.AddCutPCMCalo("0008e113","0dm00009f9730000dge0404000","4117901067832230000","0163103100b00010"); // ExC 95
    cuts.AddCutPCMCalo("0008e113","0dm00009f9730000dge0404000","4117901067b32230000","0163103100b00010"); // ExC 95 + TCard > 50
    cuts.AddCutPCMCalo("0008e113","0dm00009f9730000dge0404000","4117901067e32230000","0163103100b00010"); // ExC 97 + TCard > 50
  } else if ( trainConfig == 3117){ // EDC  ///   R Bins //  1
    cuts.AddCutPCMCalo("0008e113","0d200009f9730000dge0404000","4117901067032230000","0163103100b00010"); // RBins    min = 5,      max = 180
    cuts.AddCutPCMCalo("0008e113","0dm00009f9730000dge0404000","4117901067032230000","0163103100b00010"); // RBins    min = 5,      max = 180 without 55 -72
    cuts.AddCutPCMCalo("0008e113","0dh00009f9730000dge0404000","4117901067032230000","0163103100b00010"); // RBins    min = 5,      max = 13
    cuts.AddCutPCMCalo("0008e113","0di00009f9730000dge0404000","4117901067032230000","0163103100b00010"); // RBins    min = 13,     max = 33.5
  } else if ( trainConfig == 3118){ // EDC  ///   R Bins // 2
    cuts.AddCutPCMCalo("0008e113","0dj00009f9730000dge0404000","4117901067032230000","0163103100b00010"); // RBins    min = 33.5,   max = 55
    cuts.AddCutPCMCalo("0008e113","0dk00009f9730000dge0404000","4117901067032230000","0163103100b00010"); // RBins    min = 55,     max = 72
    cuts.AddCutPCMCalo("0008e113","0dl00009f9730000dge0404000","4117901067032230000","0163103100b00010"); // RBins    min = 72,     max = 95
    cuts.AddCutPCMCalo("0008e113","0dg00009f9730000dge0404000","4117901067032230000","0163103100b00010"); // RBins    min = 95,     max = 180
  } else if ( trainConfig == 3119){ // EMCAL clusters
    cuts.AddCutPCMCalo("0008e113","0dm00009f9730000dge0404000","4117901067032230000","0r63103100b00010"); // 90 degree rotation wo evt. weight
    cuts.AddCutPCMCalo("0008e113","0dm00009f9730000dge0404000","4117901067032230000","0t63103100b00010"); // random angle with evt. weight
  } else if ( trainConfig == 3120){//alpha, std 3 == <=1.0
    cuts.AddCutPCMCalo("0008e113","0dm00009f9730000dge0404000","4117901067032230000","0163106100b00010"); //6: alpha meson variation 1 0<alpha<0.8
    cuts.AddCutPCMCalo("0008e113","0dm00009f9730000dge0404000","4117901067032230000","0163105100b00010"); //5: alpha meson variation 2 0<alpha<0.75
  } else if ( trainConfig == 3121){ //opening angle, std 1 == >0.5, 1 cell diagonal
    cuts.AddCutPCMCalo("0008e113","0dm00009f9730000dge0404000","4117901067032230000","0163103100b00000"); //0: min opening angle 0    -> open
    cuts.AddCutPCMCalo("0008e113","0dm00009f9730000dge0404000","4117901067032230000","0163103100b00030"); //3: min opening angle 0.01 -> 2 cell diag



    // EG1
  } else if ( trainConfig ==3200){ // EMCAL clusters
    cuts.AddCutPCMCalo("0008d113","0dm00009f9730000dge0404000","4117901067032230000","0163103100b00010"); // Mixed bck
    cuts.AddCutPCMCalo("0008d113","0dm00009f9730000dge0404000","4117901067032230000","0r63103100b00010"); // 90 degree rotation wo evt. weight
    cuts.AddCutPCMCalo("0008d113","0dm00009f9730000dge0404000","4117901067032230000","0s63103100b00010"); // 90 degree rotation with evt. weight
    cuts.AddCutPCMCalo("0008d113","0dm00009f9730000dge0404000","4117901067032230000","0t63103100b00010"); // random angle with evt. weight
    cuts.AddCutPCMCalo("0008d113","0dm00009f9730000dge0404000","4117901067032230000","0u63103100b00010"); // random angle with multiple decays with evt. weight
  } else if (trainConfig == 3201) {  //   (fPSigSmearing, fPSigSmearingCte)
    cuts.AddCutPCMCalo("0008d113","0dm00009f9730000dge0404000","4117901067032230000","0163103100a00010"); // smearing (0.0275, 0.025)
    cuts.AddCutPCMCalo("0008d113","0dm00009f9730000dge0404000","4117901067032230000","0163103100b00010"); // smearing (0.025,  0.030)
    cuts.AddCutPCMCalo("0008d113","0dm00009f9730000dge0404000","4117901067032230000","0163103100c00010"); // smearing (0.0275, 0.020)
    cuts.AddCutPCMCalo("0008d113","0dm00009f9730000dge0404000","4117901067032230000","0163103100f00010"); // smearing (0.0275, 0.015)
    cuts.AddCutPCMCalo("0008d113","0dm00009f9730000dge0404000","4117901067032230000","0163103100g00010"); // smearing (0.025,  0.020)
    // Variations of EDC Part
  } else if (trainConfig == 3202){ // timing Cut variation  std -50+30ns
    cuts.AddCutPCMCalo("0008d113","0dm00009f9730000dge0404000","4117901017032230000","0163103100b00010"); //     -1000  +1000 ns
    cuts.AddCutPCMCalo("0008d113","0dm00009f9730000dge0404000","4117901077032230000","0163103100b00010"); //     -30    +30   ns
    cuts.AddCutPCMCalo("0008d113","0dm00009f9730000dge0404000","4117901097032230000","0163103100b00010"); //     -20    +25   ns
    cuts.AddCutPCMCalo("0008d113","0dm00009f9730000dge0404000","41179010a7032230000","0163103100b00010"); //     -12.5  +13   ns
  } else if (trainConfig == 3203){ // track matching variation
    cuts.AddCutPCMCalo("0008d113","0dm00009f9730000dge0404000","4117901060032230000","0163103100b00010"); //
    cuts.AddCutPCMCalo("0008d113","0dm00009f9730000dge0404000","4117901061032230000","0163103100b00010"); //
    cuts.AddCutPCMCalo("0008d113","0dm00009f9730000dge0404000","4117901066032230000","0163103100b00010"); //
  } else if (trainConfig == 3204){ // min nCells & M02 variation // std: min nCells = 1
    cuts.AddCutPCMCalo("0008d113","0dm00009f9730000dge0404000","4117901067031230000","0163103100b00010"); //   min nCells = 1
    cuts.AddCutPCMCalo("0008d113","0dm00009f9730000dge0404000","4117901067033230000","0163103100b00010"); //   min nCells = 3
  } else if (trainConfig == 3205){ // min energy variation std 0.7 GeV/c
    cuts.AddCutPCMCalo("0008d113","0dm00009f9730000dge0404000","4117901067002230000","0163103100b00010"); //     0.1 GeV/c
    cuts.AddCutPCMCalo("0008d113","0dm00009f9730000dge0404000","4117901067012230000","0163103100b00010"); //     0.5 GeV/c
    cuts.AddCutPCMCalo("0008d113","0dm00009f9730000dge0404000","4117901067022230000","0163103100b00010"); //     0.6 GeV/c
    cuts.AddCutPCMCalo("0008d113","0dm00009f9730000dge0404000","4117901067042230000","0163103100b00010"); //     0.8 GeV/c
    cuts.AddCutPCMCalo("0008d113","0dm00009f9730000dge0404000","4117901067052230000","0163103100b00010"); //     0.9 GeV/c
  } else if (trainConfig == 3206){ // min nCells & M02 variation  // std: M02 max=0.7, min=0.1
    cuts.AddCutPCMCalo("0008d113","0dm00009f9730000dge0404000","4117901067032210000","0163103100b00010"); //   max M02    = 1
    cuts.AddCutPCMCalo("0008d113","0dm00009f9730000dge0404000","4117901067032240000","0163103100b00010"); //   max M02    = 0.4
  } else if (trainConfig == 3207){ // exotic cluster
    cuts.AddCutPCMCalo("0008d113","0dm00009f9730000dge0404000","4117901067232230000","0163103100b00010"); // ExC 99
    cuts.AddCutPCMCalo("0008d113","0dm00009f9730000dge0404000","4117901067532230000","0163103100b00010"); // ExC 97
    cuts.AddCutPCMCalo("0008d113","0dm00009f9730000dge0404000","4117901067832230000","0163103100b00010"); // ExC 95
    cuts.AddCutPCMCalo("0008d113","0dm00009f9730000dge0404000","4117901067b32230000","0163103100b00010"); // ExC 95 + TCard > 50
    cuts.AddCutPCMCalo("0008d113","0dm00009f9730000dge0404000","4117901067e32230000","0163103100b00010"); // ExC 97 + TCard > 50

    // Variations of PCM Part  // std 0dm00009f9730000dge0404000
  } else if (trainConfig == 3208) {
    cuts.AddCutPCMCalo("0008d113","0dm00069f9730000dge0404000","4117901067032230000","0163103100b00010"); // min pT 40 MeV
    cuts.AddCutPCMCalo("0008d113","0dm00049f9730000dge0404000","4117901067032230000","0163103100b00010"); // min pT 75 MeV
    cuts.AddCutPCMCalo("0008d113","0dm00019f9730000dge0404000","4117901067032230000","0163103100b00010"); // min pT 100MeV
  } else if (trainConfig == 3209) {
    cuts.AddCutPCMCalo("0008d113","0dm00008f9730000dge0404000","4117901067032230000","0163103100b00010"); // TPC cluster 35%
    cuts.AddCutPCMCalo("0008d113","0dm00006f9730000dge0404000","4117901067032230000","0163103100b00010"); // TPC cluster 70%
    cuts.AddCutPCMCalo("0008d113","0dm00009f9730000dge0604000","4117901067032230000","0163103100b00010"); // cosPA 0.9
    cuts.AddCutPCMCalo("0008d113","0dm00009f9730000dge0304000","4117901067032230000","0163103100b00010"); // cosPA 0.75
  } else if (trainConfig == 3210) {
    cuts.AddCutPCMCalo("0008d113","0dm0000939730000dge0404000","4117901067032230000","0163103100b00010"); // nsig electron   -4,5
    cuts.AddCutPCMCalo("0008d113","0dm0000969730000dge0404000","4117901067032230000","0163103100b00010"); // nsig electron -2.5,4
    cuts.AddCutPCMCalo("0008d113","0dm00009f5730000dge0404000","4117901067032230000","0163103100b00010"); // nsig pion 2,-10
    cuts.AddCutPCMCalo("0008d113","0dm00009f1730000dge0404000","4117901067032230000","0163103100b00010"); // nsig pion 0,-10
  } else if (trainConfig == 3211) {
    cuts.AddCutPCMCalo("0008d113","0dm00009f9030000dge0404000","4117901067032230000","0163103100b00010"); // pion nsig min mom 0.50 GeV/c
    cuts.AddCutPCMCalo("0008d113","0dm00009f9630000dge0404000","4117901067032230000","0163103100b00010"); // pion nsig min mom 0.25 GeV/c
    cuts.AddCutPCMCalo("0008d113","0dm00009f9760000dge0404000","4117901067032230000","0163103100b00010"); // pion nsig max mom 2.00 GeV/c
    cuts.AddCutPCMCalo("0008d113","0dm00009f9710000dge0404000","4117901067032230000","0163103100b00010"); // pion nsig max mom 5.00 GeV/c
  } else if (trainConfig ==3212){
    cuts.AddCutPCMCalo("0008d113","0dm00009f97300008ge0404000","4117901067032230000","0163103100b00010"); // qT max 0.05 1D
    cuts.AddCutPCMCalo("0008d113","0dm00009f97300003ge0404000","4117901067032230000","0163103100b00010"); // qT max 0.05 1D
    cuts.AddCutPCMCalo("0008d113","0dm00009f97300002ge0404000","4117901067032230000","0163103100b00010"); // qT max 0.06 2D
    cuts.AddCutPCMCalo("0008d113","0dm00009f97300009ge0404000","4117901067032230000","0163103100b00010"); // qT max 0.03 2D
  } else if (trainConfig == 3213) {
    cuts.AddCutPCMCalo("0008d113","0dm00009f9730000dg50404000","4117901067032230000","0163103100b00010"); // Psi pair 0.1  1D
    cuts.AddCutPCMCalo("0008d113","0dm00009f9730000dg10404000","4117901067032230000","0163103100b00010"); // Psi pair 0.1  1D
    cuts.AddCutPCMCalo("0008d113","0dm00009f9730000dg60404000","4117901067032230000","0163103100b00010"); // Psi pair 0.05 2D
    cuts.AddCutPCMCalo("0008d113","0dm00009f9730000dg80404000","4117901067032230000","0163103100b00010"); // Psi pair 0.2  2D
  } else if ( trainConfig == 3214){ // qT 2D pT dep
    cuts.AddCutPCMCalo("0008d113","0dm00009f9730000dge0404000","4117901067032230000","0163103100b00010"); // qT<0.110pT (2D) alpha<0.99
    cuts.AddCutPCMCalo("0008d113","0dm00009f9730000c259404000","4117901067032230000","0163103100b00010"); // qT<0.110pT (2D) alpha<0.99
    cuts.AddCutPCMCalo("0008d113","0dm00009f9730000a259404000","4117901067032230000","0163103100b00010"); // qT<0.125pT (2D) alpha<0.99
    cuts.AddCutPCMCalo("0008d113","0dm00009f9730000e259404000","4117901067032230000","0163103100b00010"); // qT<0.130pT (2D) alpha<0.99
  } else if ( trainConfig == 3215){ // chi2, PsiPair cont.
    cuts.AddCutPCMCalo("0008d113","0dm00009f97300008fd0404000","4117901067032230000","0163103100b00010"); // PsiPair<0.15exp(-0.065chi2)
    cuts.AddCutPCMCalo("0008d113","0dm00009f97300008ge0404000","4117901067032230000","0163103100b00010"); // PsiPair<0.18exp(-0.055chi2)
    cuts.AddCutPCMCalo("0008d113","0dm00009f97300008hf0404000","4117901067032230000","0163103100b00010"); // PsiPair<0.20exp(-0.050chi2)
  } else if (trainConfig == 3216){ // exotic cluster in triggered data
    cuts.AddCutPCMCalo("0008d113","0dm00009f9730000dge0404000","4117901067232230000","0163103100b00010"); // ExC 99
    cuts.AddCutPCMCalo("0008d113","0dm00009f9730000dge0404000","4117901067532230000","0163103100b00010"); // ExC 97
    cuts.AddCutPCMCalo("0008d113","0dm00009f9730000dge0404000","4117901067832230000","0163103100b00010"); // ExC 95
    cuts.AddCutPCMCalo("0008d113","0dm00009f9730000dge0404000","4117901067b32230000","0163103100b00010"); // ExC 95 + TCard > 50
    cuts.AddCutPCMCalo("0008d113","0dm00009f9730000dge0404000","4117901067e32230000","0163103100b00010"); // ExC 97 + TCard > 50
  } else if ( trainConfig == 3217){ // EDC  ///   R Bins // 1
    cuts.AddCutPCMCalo("0008d113","0d200009f9730000dge0404000","4117901067032230000","0163103100b00010"); // RBins    min = 5,      max = 180
    cuts.AddCutPCMCalo("0008d113","0dm00009f9730000dge0404000","4117901067032230000","0163103100b00010"); // RBins    min = 5,      max = 180 without 55 -72
    cuts.AddCutPCMCalo("0008d113","0dh00009f9730000dge0404000","4117901067032230000","0163103100b00010"); // RBins    min = 5,      max = 13
    cuts.AddCutPCMCalo("0008d113","0di00009f9730000dge0404000","4117901067032230000","0163103100b00010"); // RBins    min = 13,     max = 33.5
  } else if ( trainConfig == 3218){ // EDC  ///   R Bins // 2
    cuts.AddCutPCMCalo("0008d113","0dj00009f9730000dge0404000","4117901067032230000","0163103100b00010"); // RBins    min = 33.5,   max = 55
    cuts.AddCutPCMCalo("0008d113","0dk00009f9730000dge0404000","4117901067032230000","0163103100b00010"); // RBins    min = 55,     max = 72
    cuts.AddCutPCMCalo("0008d113","0dl00009f9730000dge0404000","4117901067032230000","0163103100b00010"); // RBins    min = 72,     max = 95
    cuts.AddCutPCMCalo("0008d113","0dg00009f9730000dge0404000","4117901067032230000","0163103100b00010"); // RBins    min = 95,     max = 180
  } else if ( trainConfig == 3219){ // EMCAL clusters
    cuts.AddCutPCMCalo("0008d113","0dm00009f9730000dge0404000","4117901067032230000","0r63103100b00010"); // 90 degree rotation wo evt. weight
    cuts.AddCutPCMCalo("0008d113","0dm00009f9730000dge0404000","4117901067032230000","0t63103100b00010"); // random angle with evt. weight
  } else if ( trainConfig == 3220){//alpha, std 3 == <=1.0
    cuts.AddCutPCMCalo("0008d113","0dm00009f9730000dge0404000","4117901067032230000","0163106100b00010"); //6: alpha meson variation 1 0<alpha<0.8
    cuts.AddCutPCMCalo("0008d113","0dm00009f9730000dge0404000","4117901067032230000","0163105100b00010"); //5: alpha meson variation 2 0<alpha<0.75
  } else if ( trainConfig == 3221){ //opening angle, std 1 == >0.5, 1 cell diagonal
    cuts.AddCutPCMCalo("0008d113","0dm00009f9730000dge0404000","4117901067032230000","0163103100b00000"); //0: min opening angle 0    -> open
    cuts.AddCutPCMCalo("0008d113","0dm00009f9730000dge0404000","4117901067032230000","0163103100b00030"); //3: min opening angle 0.01 -> 2 cell diag




    // systematics for PCM-EDC  13 TeV- LowB

  } else if (trainConfig == 3300) {   // std Cut
    cuts.AddCutPCMCalo("00010113", "0dm00089f9730000iih0404000","4117901067032230000","0163103100000010");
    cuts.AddCutPCMCalo("00010113", "0dm00089f9730000iih0404000","4117901067032230000","0163103100b00010");
    cuts.AddCutPCMCalo("00010113", "0dm00089f9730000iih0404000","4117900067032230000","0163103100000010");
    cuts.AddCutPCMCalo("00010113", "0dm00089f9730000iih0404000","4117937067032230000","0163103100000010");

  } else if (trainConfig == 3301) {   // min pT variations
    cuts.AddCutPCMCalo("00010113", "0dm00079f9730000iih0404000","4117901067032230000","0163103100b00010");  // eta < 0.8  // remove  55-72 bin, min pT 0  (40) MeV
    cuts.AddCutPCMCalo("00010113", "0dm000p9f9730000iih0404000","4117901067032230000","0163103100b00010");  // eta < 0.8  // remove  55-72 bin, min pT 30 (75) MeV
    cuts.AddCutPCMCalo("00010113", "0dm00069f9730000iih0404000","4117901067032230000","0163103100b00010");  // eta < 0.8  // remove  55-72 bin, min pT    100MeV

  } else if (trainConfig == 3302) {   // TPC clusters, cosPA
    cuts.AddCutPCMCalo("00010113", "0dm00088f9730000iih0404000","4117901067032230000","0163103100b00010");  // TPC cluster 35%
    cuts.AddCutPCMCalo("00010113", "0dm00086f9730000iih0404000","4117901067032230000","0163103100b00010");  // TPC cluster 70%
    cuts.AddCutPCMCalo("00010113", "0dm00089f9730000iih0604000","4117901067032230000","0163103100b00010");  // cosPA 0.9
    cuts.AddCutPCMCalo("00010113", "0dm00089f9730000iih0304000","4117901067032230000","0163103100b00010");  // cosPA 0.75

  } else if (trainConfig == 3303) {   // TPC clusters, cosPA
    cuts.AddCutPCMCalo("00010113", "0dm0008939730000iih0404000","4117901067032230000","0163103100b00010");  // nsig electron   -4,5
    cuts.AddCutPCMCalo("00010113", "0dm0008969730000iih0404000","4117901067032230000","0163103100b00010");  // nsig electron -2.5,4
    cuts.AddCutPCMCalo("00010113", "0dm00089f5730000iih0404000","4117901067032230000","0163103100b00010");  // nsig pion 2,-10
    cuts.AddCutPCMCalo("00010113", "0dm00089f1730000iih0404000","4117901067032230000","0163103100b00010");  // nsig pion 0,-10

  } else if (trainConfig == 3304) {
    cuts.AddCutPCMCalo("00010113", "0dm00089f9030000iih0404000","4117901067032230000","0163103100b00010");  // pion nsig min mom 0.50 GeV/c
    cuts.AddCutPCMCalo("00010113", "0dm00089f9630000iih0404000","4117901067032230000","0163103100b00010");  // pion nsig min mom 0.25 GeV/c
    cuts.AddCutPCMCalo("00010113", "0dm00089f9760000iih0404000","4117901067032230000","0163103100b00010");  // pion nsig max mom 2.00 GeV/c
    cuts.AddCutPCMCalo("00010113", "0dm00089f9710000iih0404000","4117901067032230000","0163103100b00010");  // pion nsig max mom 5.00 GeV/c


  } else if (trainConfig == 3305) {   // chi2 variations
    cuts.AddCutPCMCalo("00010113", "0dm00089f9730000i1h0404000","4117901067032230000","0163103100b00010");  // chi2 50
    cuts.AddCutPCMCalo("00010113", "0dm00089f9730000ijh0404000","4117901067032230000","0163103100b00010");  // chi2 50 chi2 dep -0.085
    cuts.AddCutPCMCalo("00010113", "0dm00089f9730000ifh0404000","4117901067032230000","0163103100b00010");  // chi2 50 chi2 dep -0.065
    cuts.AddCutPCMCalo("00010113", "0dm00089f9730000iih0400000","4117901067032230000","0163103100b00010");  // remove reject close v0

  } else if (trainConfig == 3306) {   // Psi pair variations
    cuts.AddCutPCMCalo("00010113", "0dm00089f9730000iif0404000","4117901067032230000","0163103100b00010");  // Psi pair 0.2 dep
    cuts.AddCutPCMCalo("00010113", "0dm00089f9730000iii0404000","4117901067032230000","0163103100b00010");  // Psi pair 0.40 dep
    cuts.AddCutPCMCalo("00010113", "0dm00089f9730000iig0404000","4117901067032230000","0163103100b00010");  // Psi pair 0.30 dep
    cuts.AddCutPCMCalo("00010113", "0dm00089227300008250404000","4117901067032230000","0163103100b00010");  // old cuts (run1)

  } else if (trainConfig == 3307) {
    cuts.AddCutPCMCalo("00010113", "0dm00089f9730000iih0404000","4117901067032230000","0263103100b00010"); // variation BG scheme track mult
    cuts.AddCutPCMCalo("00010113", "0dm00089f9730000iih0404000","4117901067032230000","0r63103100b00010"); // variation BG scheme track mult
    cuts.AddCutPCMCalo("00010113", "0dm00089f9730000iih0404000","4117901067032230000","0163107100b00010"); // alpha meson 0.85
    cuts.AddCutPCMCalo("00010113", "0dm00089f9730000iih0404000","4117901067032230000","0163105100b00010"); // alpha meson 0.75

  } else if (trainConfig == 3308) {   // qT variations
    cuts.AddCutPCMCalo("00010113", "0dm00089f9730000jih0404000","4117901067032230000","0163103100b00010");  // qT max 0.040, qtptmax 0.25
    cuts.AddCutPCMCalo("00010113", "0dm00089f9730000kih0404000","4117901067032230000","0163103100b00010");  // qT max 0.045, qtptmax 0.3
    //    cuts.AddCutPCMCalo("00010113", "0dm00089f9730000fih0404000","4117901067032230000","0163103100b00010");  // qT max 0.070, qtptmax 0.16
    // cuts.AddCutPCMCalo("00010113", "0dm00089f9730000dih0404000", "0152101500000000"); // alpha meson pT dependent

    } else if ( trainConfig ==3309){ // EMCAL clusters
    cuts.AddCutPCMCalo("00010113","0dm00089f9730000iih0404000","4117901067032230000","0163103100b00010"); // Mixed bck
    cuts.AddCutPCMCalo("00010113","0dm00089f9730000iih0404000","4117901067032230000","0r63103100b00010"); // 90 degree rotation wo evt. weight
    cuts.AddCutPCMCalo("00010113","0dm00089f9730000iih0404000","4117901067032230000","0s63103100b00010"); // 90 degree rotation with evt. weight
    cuts.AddCutPCMCalo("00010113","0dm00089f9730000iih0404000","4117901067032230000","0t63103100b00010"); // random angle with evt. weight
    cuts.AddCutPCMCalo("00010113","0dm00089f9730000iih0404000","4117901067032230000","0u63103100b00010"); // random angle with multiple decays with evt. weight
  } else if (trainConfig == 3310) {  //   (fPSigSmearing, fPSigSmearingCte)
    cuts.AddCutPCMCalo("00010113","0dm00089f9730000iih0404000","4117901067032230000","0163103100a00010"); // smearing (0.0275, 0.025)
    cuts.AddCutPCMCalo("00010113","0dm00089f9730000iih0404000","4117901067032230000","0163103100b00010"); // smearing (0.025,  0.030)
    cuts.AddCutPCMCalo("00010113","0dm00089f9730000iih0404000","4117901067032230000","0163103100c00010"); // smearing (0.0275, 0.020)
    cuts.AddCutPCMCalo("00010113","0dm00089f9730000iih0404000","4117901067032230000","0163103100f00010"); // smearing (0.0275, 0.015)
    cuts.AddCutPCMCalo("00010113","0dm00089f9730000iih0404000","4117901067032230000","0163103100g00010"); // smearing (0.025,  0.020)
    // Variations of EDC Part
  } else if (trainConfig == 3311){ // timing Cut variation  std -50+30ns
    cuts.AddCutPCMCalo("00010113","0dm00089f9730000iih0404000","4117901017032230000","0163103100b00010"); //     -1000  +1000 ns
    cuts.AddCutPCMCalo("00010113","0dm00089f9730000iih0404000","4117901077032230000","0163103100b00010"); //     -30    +30   ns
    cuts.AddCutPCMCalo("00010113","0dm00089f9730000iih0404000","4117901097032230000","0163103100b00010"); //     -20    +25   ns
    cuts.AddCutPCMCalo("00010113","0dm00089f9730000iih0404000","41179010a7032230000","0163103100b00010"); //     -12.5  +13   ns
  } else if (trainConfig == 3312){ // track matching variation
    cuts.AddCutPCMCalo("00010113","0dm00089f9730000iih0404000","4117901060032230000","0163103100b00010"); //
    cuts.AddCutPCMCalo("00010113","0dm00089f9730000iih0404000","4117901061032230000","0163103100b00010"); //
    cuts.AddCutPCMCalo("00010113","0dm00089f9730000iih0404000","4117901066032230000","0163103100b00010"); //
  } else if (trainConfig == 3313){ // min nCells & M02 variation // std: min nCells = 1
    cuts.AddCutPCMCalo("00010113","0dm00089f9730000iih0404000","4117901067031230000","0163103100b00010"); //   min nCells = 1
    cuts.AddCutPCMCalo("00010113","0dm00089f9730000iih0404000","4117901067033230000","0163103100b00010"); //   min nCells = 3
  } else if (trainConfig == 3314){ // min energy variation std 0.7 GeV/c
    cuts.AddCutPCMCalo("00010113","0dm00089f9730000iih0404000","4117901067002230000","0163103100b00010"); //     0.1 GeV/c
    cuts.AddCutPCMCalo("00010113","0dm00089f9730000iih0404000","4117901067012230000","0163103100b00010"); //     0.5 GeV/c
    cuts.AddCutPCMCalo("00010113","0dm00089f9730000iih0404000","4117901067022230000","0163103100b00010"); //     0.6 GeV/c
    cuts.AddCutPCMCalo("00010113","0dm00089f9730000iih0404000","4117901067042230000","0163103100b00010"); //     0.8 GeV/c
    cuts.AddCutPCMCalo("00010113","0dm00089f9730000iih0404000","4117901067052230000","0163103100b00010"); //     0.9 GeV/c
  } else if (trainConfig == 3315){ // min nCells & M02 variation  // std: M02 max=0.7, min=0.1
    cuts.AddCutPCMCalo("00010113","0dm00089f9730000iih0404000","4117901067032210000","0163103100b00010"); //   max M02    = 1
    cuts.AddCutPCMCalo("00010113","0dm00089f9730000iih0404000","4117901067032240000","0163103100b00010"); //   max M02    = 0.4
  } else if (trainConfig == 3316){ // exotic cluster
    cuts.AddCutPCMCalo("00010113","0dm00089f9730000iih0404000","4117901067232230000","0163103100b00010"); // ExC 99
    cuts.AddCutPCMCalo("00010113","0dm00089f9730000iih0404000","4117901067532230000","0163103100b00010"); // ExC 97
    cuts.AddCutPCMCalo("00010113","0dm00089f9730000iih0404000","4117901067832230000","0163103100b00010"); // ExC 95
    cuts.AddCutPCMCalo("00010113","0dm00089f9730000iih0404000","4117901067b32230000","0163103100b00010"); // ExC 95 + TCard > 50
    cuts.AddCutPCMCalo("00010113","0dm00089f9730000iih0404000","4117901067e32230000","0163103100b00010"); // ExC 97 + TCard > 50
  } else if ( trainConfig == 3317){ //opening angle, std 1 == >0.5, 1 cell diagonal
    cuts.AddCutPCMCalo("00010113","0dm00089f9730000iih0404000","4117901067032230000","0163103100b00000"); //0: min opening angle 0    -> open
    cuts.AddCutPCMCalo("00010113","0dm00089f9730000iih0404000","4117901067032230000","0163103100b00030"); //3: min opening angle 0.01 -> 2 cell diag
  } else if (trainConfig == 3318){ // EMCAL+DCAL clusters standard cuts, INT7, NL , std TM, tight timing
    cuts.AddCutPCMCalo("00010113","0dm00089f9730000iih0404000","4117911067e32230000","0163103100b00010"); // INT7 NL11
    cuts.AddCutPCMCalo("00010113","0dm00089f9730000iih0404000","4117912067e32230000","0163103100b00010"); // INT7 NL12
    cuts.AddCutPCMCalo("00010113","0dm00089f9730000iih0404000","4117921067e32230000","0163103100b00010"); // INT7 NL21
    cuts.AddCutPCMCalo("00010113","0dm00089f9730000iih0404000","4117922067e32230000","0163103100b00010"); // INT7 NL22






    // systematics for PCM-PHOS  13 TeV- LowB

  } else if (trainConfig == 3400) {   // std Cut
    cuts.AddCutPCMCalo("00010113", "0dm00089f9730000iih0404000","24466190sa01cc00000","0163103100000010");

  } else if (trainConfig == 3401) {   // min pT variations
    cuts.AddCutPCMCalo("00010113", "0dm00079f9730000iih0404000","24466190sa01cc00000","0163103100000010");  // eta < 0.8  // remove  55-72 bin, min pT 0  (40) MeV
    cuts.AddCutPCMCalo("00010113", "0dm000p9f9730000iih0404000","24466190sa01cc00000","0163103100000010");  // eta < 0.8  // remove  55-72 bin, min pT 30 (75) MeV
    cuts.AddCutPCMCalo("00010113", "0dm00069f9730000iih0404000","24466190sa01cc00000","0163103100000010");  // eta < 0.8  // remove  55-72 bin, min pT    100MeV
    
  } else if (trainConfig == 3402) {   // TPC clusters, cosPA
    cuts.AddCutPCMCalo("00010113", "0dm00088f9730000iih0404000","24466190sa01cc00000","0163103100000010");  // TPC cluster 35%
    cuts.AddCutPCMCalo("00010113", "0dm00086f9730000iih0404000","24466190sa01cc00000","0163103100000010");  // TPC cluster 70%
    cuts.AddCutPCMCalo("00010113", "0dm00089f9730000iih0604000","24466190sa01cc00000","0163103100000010");  // cosPA 0.9
    cuts.AddCutPCMCalo("00010113", "0dm00089f9730000iih0304000","24466190sa01cc00000","0163103100000010");  // cosPA 0.75

  } else if (trainConfig == 3403) {   // TPC clusters, cosPA
    cuts.AddCutPCMCalo("00010113", "0dm0008939730000iih0404000","24466190sa01cc00000","0163103100000010");  // nsig electron   -4,5
    cuts.AddCutPCMCalo("00010113", "0dm0008969730000iih0404000","24466190sa01cc00000","0163103100000010");  // nsig electron -2.5,4
    cuts.AddCutPCMCalo("00010113", "0dm00089f5730000iih0404000","24466190sa01cc00000","0163103100000010");  // nsig pion 2,-10
    cuts.AddCutPCMCalo("00010113", "0dm00089f1730000iih0404000","24466190sa01cc00000","0163103100000010");  // nsig pion 0,-10

  } else if (trainConfig == 3404) {
    cuts.AddCutPCMCalo("00010113", "0dm00089f9030000iih0404000","24466190sa01cc00000","0163103100000010");  // pion nsig min mom 0.50 GeV/c
    cuts.AddCutPCMCalo("00010113", "0dm00089f9630000iih0404000","24466190sa01cc00000","0163103100000010");  // pion nsig min mom 0.25 GeV/c
    cuts.AddCutPCMCalo("00010113", "0dm00089f9760000iih0404000","24466190sa01cc00000","0163103100000010");  // pion nsig max mom 2.00 GeV/c
    cuts.AddCutPCMCalo("00010113", "0dm00089f9710000iih0404000","24466190sa01cc00000","0163103100000010");  // pion nsig max mom 5.00 GeV/c


  } else if (trainConfig == 3405) {   // chi2 variations
    cuts.AddCutPCMCalo("00010113", "0dm00089f9730000i1h0404000","24466190sa01cc00000","0163103100000010");  // chi2 50
    cuts.AddCutPCMCalo("00010113", "0dm00089f9730000ijh0404000","24466190sa01cc00000","0163103100000010");  // chi2 50 chi2 dep -0.085
    cuts.AddCutPCMCalo("00010113", "0dm00089f9730000ifh0404000","24466190sa01cc00000","0163103100000010");  // chi2 50 chi2 dep -0.065
    cuts.AddCutPCMCalo("00010113", "0dm00089f9730000iih0400000","24466190sa01cc00000","0163103100000010");  // remove reject close v0

  } else if (trainConfig == 3406) {   // Psi pair variations
    cuts.AddCutPCMCalo("00010113", "0dm00089f9730000iif0404000","24466190sa01cc00000","0163103100000010");  // Psi pair 0.2 dep
    cuts.AddCutPCMCalo("00010113", "0dm00089f9730000iii0404000","24466190sa01cc00000","0163103100000010");  // Psi pair 0.40 dep
    cuts.AddCutPCMCalo("00010113", "0dm00089f9730000iig0404000","24466190sa01cc00000","0163103100000010");  // Psi pair 0.30 dep
    cuts.AddCutPCMCalo("00010113", "0dm00089227300008250404000","24466190sa01cc00000","0163103100000010");  // old cuts (run1)

  } else if (trainConfig == 3407) {
    cuts.AddCutPCMCalo("00010113", "0dm00089f9730000iih0404000","24466190sa01cc00000","0263103100000010"); // variation BG scheme track mult
    cuts.AddCutPCMCalo("00010113", "0dm00089f9730000iih0404000","24466190sa01cc00000","0r63103100000010"); // variation BG scheme track mult
    cuts.AddCutPCMCalo("00010113", "0dm00089f9730000iih0404000","24466190sa01cc00000","0163107100000010"); // alpha meson 0.85
    cuts.AddCutPCMCalo("00010113", "0dm00089f9730000iih0404000","24466190sa01cc00000","0163105100000010"); // alpha meson 0.75

  } else if (trainConfig == 3408) {   // qT variations
    cuts.AddCutPCMCalo("00010113", "0dm00089f9730000jih0404000","24466190sa01cc00000","0163103100000010");  // qT max 0.040, qtptmax 0.25
    cuts.AddCutPCMCalo("00010113", "0dm00089f9730000kih0404000","24466190sa01cc00000","0163103100000010");  // qT max 0.045, qtptmax 0.3
    //    cuts.AddCutPCMCalo("00010113", "0dm00089f9730000fih0404000","24466190sa01cc00000","0163103100000010");  // qT max 0.070, qtptmax 0.16
    // cuts.AddCutPCMCalo("00010113", "0dm00089f9730000dih0404000", "0152101500000000"); // alpha meson pT dependent

    // Variations of PHOS Part
    //Standard: "24466190sa01cc00000"
  } else if (trainConfig == 3409){ // timing Cut variation  std -30+30ns
    //                                                                  |
    cuts.AddCutPCMCalo("00010113","0dm00089f9730000iih0404000","244661901a01cc00000","0163103100000010"); //1:     -1000  +1000 ns
    cuts.AddCutPCMCalo("00010113","0dm00089f9730000iih0404000","244661905a01cc00000","0163103100000010"); //5:     -50    +50   ns
    cuts.AddCutPCMCalo("00010113","0dm00089f9730000iih0404000","244661907a01cc00000","0163103100000010"); //7:     -30    +30   ns
    cuts.AddCutPCMCalo("00010113","0dm00089f9730000iih0404000","244661909a01cc00000","0163103100000010"); //9:     -20    +25   ns
    cuts.AddCutPCMCalo("00010113","0dm00089f9730000iih0404000","24466190aa01cc00000","0163103100000010"); //a:     -12.5  +13   ns
  } else if (trainConfig == 3410){ // Timing Efficiency Variations std s == 30ns
    //                                                                  |
    cuts.AddCutPCMCalo("00010113","0dm00089f9730000iih0404000","24466190ra01cc00000","0163103100000010"); //r:     LowPt from MB, 30ns
    cuts.AddCutPCMCalo("00010113","0dm00089f9730000iih0404000","24466190ta01cc00000","0163103100000010"); //t:     25ns
    cuts.AddCutPCMCalo("00010113","0dm00089f9730000iih0404000","24466190ua01cc00000","0163103100000010"); //u:     50ns
  } else if (trainConfig == 3411){ // track matching variation std a == pt dependent
    //                                                                   |
    cuts.AddCutPCMCalo("00010113","0dm00089f9730000iih0404000","24466190s001cc00000","0163103100000010"); //0
    cuts.AddCutPCMCalo("00010113","0dm00089f9730000iih0404000","24466190s101cc00000","0163103100000010"); //1
    cuts.AddCutPCMCalo("00010113","0dm00089f9730000iih0404000","24466190s401cc00000","0163103100000010"); //4
    cuts.AddCutPCMCalo("00010113","0dm00089f9730000iih0404000","24466190s501cc00000","0163103100000010"); //5
    cuts.AddCutPCMCalo("00010113","0dm00089f9730000iih0404000","24466190s601cc00000","0163103100000010"); //6
  } else if (trainConfig == 3412){ // min energy variation std 0.3 GeV/c
    //                                                                     |
    cuts.AddCutPCMCalo("00010113","0dm00089f9730000iih0404000","24466190sa00cc00000","0163103100000010"); //0:     off
    cuts.AddCutPCMCalo("00010113","0dm00089f9730000iih0404000","24466190sa09cc00000","0163103100000010"); //9:     0.1 GeV/c
    cuts.AddCutPCMCalo("00010113","0dm00089f9730000iih0404000","24466190sa02cc00000","0163103100000010"); //2:     0.5 GeV/c
    cuts.AddCutPCMCalo("00010113","0dm00089f9730000iih0404000","24466190sa03cc00000","0163103100000010"); //3:     0.6 GeV/c
    cuts.AddCutPCMCalo("00010113","0dm00089f9730000iih0404000","24466190sa04cc00000","0163103100000010"); //4:     0.7 GeV/c
    cuts.AddCutPCMCalo("00010113","0dm00089f9730000iih0404000","24466190sa05cc00000","0163103100000010"); //5:     0.8 GeV/c
  } else if (trainConfig == 3413){ // min nCells & M02 variation, std cc
    // std: min nCells = 2 >1GeV; M02 max=100, min=0.1
    //                                                                      ||
    cuts.AddCutPCMCalo("00010113","0dm00089f9730000iih0404000","24466190sa011000000","0163103100000010"); //100:   min nCells = 1, minM02 off
    cuts.AddCutPCMCalo("00010113","0dm00089f9730000iih0404000","24466190sa012200000","0163103100000010"); //220:   min nCells = 2, all E
    cuts.AddCutPCMCalo("00010113","0dm00089f9730000iih0404000","24466190sa013200000","0163103100000010"); //320:   min nCells = 3, all E
    cuts.AddCutPCMCalo("00010113","0dm00089f9730000iih0404000","24466190sa01dc00000","0163103100000010"); //dc0:   min nCells = 3, E>1GeV; minM02==0.1 off for E<1GeV
    cuts.AddCutPCMCalo("00010113","0dm00089f9730000iih0404000","24466190sa01cd00000","0163103100000010"); //cd0:   min nCells = 2, E>1GeV; minM02==0.2 off for E<1GeV
    cuts.AddCutPCMCalo("00010113","0dm00089f9730000iih0404000","24466190sa01cc70000","0163103100000010"); //cc7:   maxM02 == 1.3
    cuts.AddCutPCMCalo("00010113","0dm00089f9730000iih0404000","24466190sa01cc80000","0163103100000010"); //cc8:   maxM02 == 2.5
  } else if (trainConfig == 3414){ // reconstructed conversion, std 0==off       |
    cuts.AddCutPCMCalo("00010113","0dm00089f9730000iih0404000","24466190sa01cc00100","0163103100000010"); //1:   rec. conv. 0.02
    cuts.AddCutPCMCalo("00010113","0dm00089f9730000iih0404000","24466190sa01cc00200","0163103100000010"); //2:   rec. conv. 0.025
    cuts.AddCutPCMCalo("00010113","0dm00089f9730000iih0404000","24466190sa01cc00300","0163103100000010"); //3:   rec. conv. 0.03
    cuts.AddCutPCMCalo("00010113","0dm00089f9730000iih0404000","24466190sa01cc00400","0163103100000010"); //4:   rec. conv. 0.035
    cuts.AddCutPCMCalo("00010113","0dm00089f9730000iih0404000","24466190sa01cc00500","0163103100000010"); //5:   rec. conv. 0.04
  } else if ( trainConfig == 3415){ // pcmphos  ///   R Bins // 1
    cuts.AddCutPCMCalo("00010113","0d200009f9730000dge0404000","24466190sa01cc00000","0163103100000010"); // RBins    min = 5,      max = 180
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163103100000010"); // RBins    min = 5,      max = 180 without 55 -72
    cuts.AddCutPCMCalo("00010113","0dh00009f9730000dge0404000","24466190sa01cc00000","0163103100000010"); // RBins    min = 5,      max = 13
    cuts.AddCutPCMCalo("00010113","0di00009f9730000dge0404000","24466190sa01cc00000","0163103100000010"); // RBins    min = 13,     max = 33.5
  } else if ( trainConfig == 3416){ // pcmphos  ///   R Bins // 2
    cuts.AddCutPCMCalo("00062113","0dj00009f9730000dge0404000","24466190sa01cc00000","0163103100000010"); // RBins    min = 33.5,   max = 55
    cuts.AddCutPCMCalo("00062113","0dk00009f9730000dge0404000","24466190sa01cc00000","0163103100000010"); // RBins    min = 55,     max = 72
    cuts.AddCutPCMCalo("00062113","0dl00009f9730000dge0404000","24466190sa01cc00000","0163103100000010"); // RBins    min = 72,     max = 95
    cuts.AddCutPCMCalo("00062113","0dg00009f9730000dge0404000","24466190sa01cc00000","0163103100000010"); // RBins    min = 95,     max = 180





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

  Int_t numberOfCuts = cuts.GetNCuts();

  TList *EventCutList = new TList();
  TList *ConvCutList = new TList();
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
  ConvCutList->SetOwner(kTRUE);
  AliConversionPhotonCuts **analysisCuts = new AliConversionPhotonCuts*[numberOfCuts];
  ClusterCutList->SetOwner(kTRUE);
  AliCaloPhotonCuts **analysisClusterCuts = new AliCaloPhotonCuts*[numberOfCuts];
  MesonCutList->SetOwner(kTRUE);
  AliConversionMesonCuts **analysisMesonCuts = new AliConversionMesonCuts*[numberOfCuts];

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
      if (enableLightOutput == 4) fTrackMatcher->SetLightOutput(kTRUE);
      mgr->AddTask(fTrackMatcher);
      mgr->ConnectInput(fTrackMatcher,0,cinput);
    }
    analysisEventCuts[i] = new AliConvEventCuts();

    TString EventCutPos = cuts.GetEventCut(i);
    EventCutPos = EventCutPos(3,2);
    TString ClusterCutPos = cuts.GetClusterCut(i);
    ClusterCutPos = ClusterCutPos(0,1);
    TString TriggerHelperName = Form("CaloTriggerHelper_%s_%i_%i", cuts.GetEventCut(i).Data(),enableTriggerMimicking,TriggerMimickingDDLEffiFlag);
    if( (!(AliCaloTriggerMimicHelper*)mgr->GetTask(TriggerHelperName.Data())) && (!ClusterCutPos.CompareTo("2")) && ( enableTriggerMimicking==3 || enableTriggerMimicking==4 ) ){
      AliCaloTriggerMimicHelper* fMimickHelper = new AliCaloTriggerMimicHelper(TriggerHelperName.Data(), caloCutPos.Atoi(), isMC);
      task->SetCaloTriggerHelperName(TriggerHelperName.Data());
      analysisEventCuts[i]->SetCaloTriggerHelperName(TriggerHelperName.Data());
      if (enableTriggerMimicking==3 || enableTriggerMimicking==4){
          fMimickHelper->SetPHOSTrigger(AliCaloTriggerMimicHelper::kPHOSAny) ;
      } else {
          fMimickHelper->SetPHOSTrigger(AliCaloTriggerMimicHelper::kPHOSL0) ;
      }
      fMimickHelper->SetEfficiencyChoiceOption_TriggerHelper(TriggerMimickingDDLEffiFlag);
      mgr->AddTask(fMimickHelper);
      mgr->ConnectInput(fMimickHelper,0,cinput);
      if (enableLightOutput>=1){
          fMimickHelper->SetLightOutput(enableLightOutput);
      } else if (enableExtMatchAndQA>3){
          fMimickHelper->SetLightOutput(0);
      } else {
          fMimickHelper->SetLightOutput(1);
      }
      if (!EventCutPos.CompareTo("10")){
          fMimickHelper->SetTriggerHelperRunMode(3);
      }
      if (!EventCutPos.CompareTo("14")){
          fMimickHelper->SetTriggerHelperRunMode(2);
      }
    }


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

    if (doWeightingPart) analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile( kTRUE, kTRUE, kFALSE, fileNamePtWeights,
                                                                                           mcInputNamePi0, mcInputNameEta, "",fitNamePi0,fitNameEta);

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

    if (enableMultiplicityWeighting){
      cout << "enableling mult weighting" << endl;
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
    if (enableLightOutput == 1 || enableLightOutput == 2 ) analysisEventCuts[i]->SetLightOutput(1);
    if (enableLightOutput == 4) analysisEventCuts[i]->SetLightOutput(2);
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
    if (enableLightOutput == 1 || enableLightOutput == 2 || enableLightOutput ==5 ) analysisCuts[i]->SetLightOutput(1);
    if (enableLightOutput == 4) analysisCuts[i]->SetLightOutput(2);
    if (enableLightOutput == 0) analysisCuts[i]->SetPlotTrackPID(kTRUE);
    analysisCuts[i]->InitializeCutsFromCutString((cuts.GetPhotonCut(i)).Data());
    analysisCuts[i]->SetIsHeavyIon(isHeavyIon);
    ConvCutList->Add(analysisCuts[i]);
    analysisCuts[i]->SetFillCutHistograms("",kFALSE);

    analysisClusterCuts[i] = new AliCaloPhotonCuts(isMC);
    analysisClusterCuts[i]->SetHistoToModifyAcceptance(histoAcc);
    analysisClusterCuts[i]->SetV0ReaderName(V0ReaderName);
    analysisClusterCuts[i]->SetCorrectionTaskSetting(corrTaskSetting);
    analysisClusterCuts[i]->SetCaloTrackMatcherName(TrackMatcherName);
    if (enableLightOutput == 1 || enableLightOutput == 2 || enableLightOutput ==5 ) analysisClusterCuts[i]->SetLightOutput(1);
    if (enableLightOutput == 4) analysisClusterCuts[i]->SetLightOutput(2);
    analysisClusterCuts[i]->InitializeCutsFromCutString((cuts.GetClusterCut(i)).Data());
    ClusterCutList->Add(analysisClusterCuts[i]);
    analysisClusterCuts[i]->SetExtendedMatchAndQA(enableExtMatchAndQA);

    analysisMesonCuts[i] = new AliConversionMesonCuts();
    if (enableLightOutput > 0 && enableLightOutput != 4) analysisMesonCuts[i]->SetLightOutput(1);
    if (enableLightOutput == 4) analysisMesonCuts[i]->SetLightOutput(2);
    analysisMesonCuts[i]->SetRunningMode(2);
    analysisMesonCuts[i]->InitializeCutsFromCutString((cuts.GetMesonCut(i)).Data());
    MesonCutList->Add(analysisMesonCuts[i]);
    analysisMesonCuts[i]->SetFillCutHistograms("");
    analysisEventCuts[i]->SetAcceptedHeader(HeaderList);
    if(doSmear) analysisMesonCuts[i]->SetDefaultSmearing(bremSmear,smearPar,smearParConst);

    if(analysisMesonCuts[i]->DoGammaSwappForBg()) analysisClusterCuts[i]->SetUseEtaPhiMapForBackCand(kTRUE);
    analysisClusterCuts[i]->SetFillCutHistograms("");

  }

  task->SetEventCutList(numberOfCuts,EventCutList);
  task->SetConversionCutList(numberOfCuts,ConvCutList);
  task->SetCaloCutList(numberOfCuts,ClusterCutList);
  task->SetMesonCutList(numberOfCuts,MesonCutList);
  task->SetMoveParticleAccordingToVertex(kTRUE);
  task->SetDoMesonAnalysis(kTRUE);
  task->SetCorrectionTaskSetting(corrTaskSetting);
  task->SetDoMesonQA(enableQAMesonTask); //Attention new switch for Pi0 QA
  task->SetDoPhotonQA(enableQAPhotonTask);  //Attention new switch small for Photon QA
  if (enableLightOutput != 4) task->SetDoClusterQA(1);  //Attention new switch small for Cluster QA
  task->SetUseTHnSparse(enableTHnSparse);
  task->SetDoTreeInvMassShowerShape(doTreeClusterShowerShape);
  task->SetEnableSortingOfMCClusLabels(enableSortingMCLabels);
  if(enableTreeConvGammaShape) task->SetDoTreeConvGammaShowerShape(kTRUE);
  if(enableExtMatchAndQA > 1){ task->SetPlotHistsExtQA(kTRUE);}

  if (initializedMatBudWeigths_existing) {
      task->SetDoMaterialBudgetWeightingOfGammasForTrueMesons(kTRUE);
  }

  //connect containers
  AliAnalysisDataContainer *coutput =
    mgr->CreateContainer(!(corrTaskSetting.CompareTo("")) ? Form("GammaConvCalo_%i",trainConfig) : Form("GammaConvCalo_%i_%s",trainConfig,corrTaskSetting.Data()), TList::Class(),
              AliAnalysisManager::kOutputContainer, Form("GammaConvCalo_%i.root",trainConfig) );

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
