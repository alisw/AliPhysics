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
//($ALIPHYSICS/PWGGA/GammaConv/AliAnalysisTaskGammaCaloMix.cxx) for
//pp together with all supporting classes
//***************************************************************************************

//***************************************************************************************
//main function
//***************************************************************************************
void AddTask_GammaCaloMix_pp(
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
  TString   settingMaxFacPtHard           = "3.",       // maximum factor between hardest jet and ptHard generated
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


  TString addTaskName                 = "AddTask_GammaCaloMix_pp";
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
  if(rmaxFacPtHardSetting->GetEntries()<1){cout << "ERROR: AddTask_GammaCaloMix_pp during parsing of settingMaxFacPtHard String '" << settingMaxFacPtHard.Data() << "'" << endl; return;}
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
  AliAnalysisTaskGammaCaloMix *task=NULL;
  task= new AliAnalysisTaskGammaCaloMix(Form("GammaCaloMix_%i",trainConfig));
  task->SetIsHeavyIon(isHeavyIon);
  task->SetIsMC(isMC);
  task->SetV0ReaderName(V0ReaderName);
  task->SetLightOutput(enableLightOutput);
  task->SetCorrectionTaskSetting(corrTaskSetting);
  task->SetTrackMatcherRunningMode(trackMatcherRunningMode);

  // here is the order of the cluster cut string
  // usually for EMCal we start with 11111: default values for            "ClusterType", "EtaMin", "EtaMax", "PhiMin", "PhiMax"
  // then two numbers for nonlinearity, e.g. 21: this is                  "NonLinearity1", "NonLinearity2"
  // Then some cuts on the clusters, e.g. 06003222: this is               "DistanceToBadChannel", "Timing", "TrackMatching", "ExoticCell", "MinEnergy", "MinNCells", "MinM02", "MaxM02"
  // finally some for now unused cuts, usually 0000: this is              "MinMaxM20", "RecConv", "MaximumDispersion", "NLM"


  // *****************************************************************************************************
  // ******************** pp 13 TeV cuts  *****************************************************
  // *****************************************************************************************************
  if (trainConfig == 1){
    cuts.AddCutCaloCalo("00010113","24466190sa01cc00000","3885501067032220000","0163103100000050");
  } else if (trainConfig == 2){  //    Timing Variations
    cuts.AddCutCaloCalo("00010113","244661907a01cc00000","3885501067032220000","0163103100000050");
    cuts.AddCutCaloCalo("00010113","244661904a01cc00000","3885501067032220000","0163103100000050");
} else if (trainConfig == 3){  //    Clus Energy Phos  Variations
    cuts.AddCutCaloCalo("00010113","24466190sa09cc00000","3885501067032220000","0163103100000050");  // PHOS energy cut 0.5
    cuts.AddCutCaloCalo("00010113","24466190sa07cc00000","3885501067032220000","0163103100000050");  // PHOS energy cut 0.2
    cuts.AddCutCaloCalo("00010113","24466190sa01cc00000","3885501067002220000","0163103100000050");    // 0.1 EMCal energy cut
    cuts.AddCutCaloCalo("00010113","24466190sa01cc00000","3885501067012220000","0163103100000050");    // 0.5 EMCal energy cut
  } else if (trainConfig == 4){
    cuts.AddCutCaloCalo("00010113","24466190sa01cc00000","1111101067032230000","0163103100000050");    // PHOS EMCAL; PHOS energy cut 0.3
    cuts.AddCutCaloCalo("00010113","24466190sa01cc00000","4117901067032230000","0163103100000050");    // PHOS EDC; PHOS energy cut 0.3
    cuts.AddCutCaloCalo("00010113","1111112067032230000","3885501067012220000","0163103100000050");    // EMCAL DCAL

} else if (trainConfig == 5){  // EMCal+DCAL clusters standard cuts, Sphericity
    cuts.AddCutCaloCalo("h0310113","24466190sa01cc00000","3885512067032220000","0163103100b00010"); //  0.    - 0.3
    cuts.AddCutCaloCalo("h0510113","24466190sa01cc00000","3885512067032220000","0163103100b00010"); //  0.    - 0.5
    cuts.AddCutCaloCalo("h3710113","24466190sa01cc00000","3885512067032220000","0163103100b00010"); //  0.3    - 0.7
    cuts.AddCutCaloCalo("h5a10113","24466190sa01cc00000","3885512067032220000","0163103100b00010"); //  0.5   - 1.
    cuts.AddCutCaloCalo("h7a10113","24466190sa01cc00000","3885512067032220000","0163103100b00010"); //  0.7   - 1.
    cuts.AddCutCaloCalo("h0a10113","24466190sa01cc00000","3885512067032220000","0163103100b00010"); //  0.    - 1.
} else if (trainConfig == 6){  // EMCal+DCAL clusters standard cuts, V0M mult selections
    cuts.AddCutCaloCalo("n0a10113","24466190sa01cc00000","3885512067032220000","0163103100b00010"); // INT7 0-100%
    cuts.AddCutCaloCalo("m0110113","24466190sa01cc00000","3885512067032220000","0163103100b00010"); // INT7 0-1%
    cuts.AddCutCaloCalo("m1510113","24466190sa01cc00000","3885512067032220000","0163103100b00010"); // INT7 1-5%
    cuts.AddCutCaloCalo("m5a10113","24466190sa01cc00000","3885512067032220000","0163103100b00010"); // INT7 5-10%
    cuts.AddCutCaloCalo("maf10113","24466190sa01cc00000","3885512067032220000","0163103100b00010"); // INT7 10-15%
    cuts.AddCutCaloCalo("mfk10113","24466190sa01cc00000","3885512067032220000","0163103100b00010"); // INT7 15-20%
    cuts.AddCutCaloCalo("n2310113","24466190sa01cc00000","3885512067032220000","0163103100b00010"); // INT7 20-30%
    cuts.AddCutCaloCalo("n3410113","24466190sa01cc00000","3885512067032220000","0163103100b00010"); // INT7 30-40%
    cuts.AddCutCaloCalo("n4510113","24466190sa01cc00000","3885512067032220000","0163103100b00010"); // INT7 40-50%
    cuts.AddCutCaloCalo("n5a10113","24466190sa01cc00000","3885512067032220000","0163103100b00010"); // INT7 50-100%


    // ******************** pp 13 TeV Triggered Data        *****************************************
  } else if (trainConfig == 10){
    cuts.AddCutCaloCalo("0008e113","24466190pa01cc00000","3885512067032220000","0163103100000050");// EG2+DG2
    cuts.AddCutCaloCalo("0008d113","24466190pa01cc00000","3885512067032220000","0163103100000050");// EG1+DG1
    cuts.AddCutCaloCalo("0009b113","24466190pa01cc00000","3885512067032220000","0163103100000050");// EGJ+DJ1


  // Variations for systematics
} else if ( trainConfig == 20){ // EMCAL clusters
    cuts.AddCutCaloCalo("00010113","24466190sa01cc00000","3885501067032230000","0163103100b00010"); // Mixed bck
    cuts.AddCutCaloCalo("00010113","24466190sa01cc00000","3885501067032230000","0r63103100b00010"); // 90 degree rotation wo evt. weight
    cuts.AddCutCaloCalo("00010113","24466190sa01cc00000","3885501067032230000","0s63103100b00010"); // 90 degree rotation with evt. weight
    cuts.AddCutCaloCalo("00010113","24466190sa01cc00000","3885501067032230000","0t63103100b00010"); // random angle with evt. weight
    cuts.AddCutCaloCalo("00010113","24466190sa01cc00000","3885501067032230000","0u63103100b00010"); // random angle with multiple decays with evt. weight
  // Variations of EDC Part
} else if (trainConfig == 21){ // EMCAL+DCAL clusters standard cuts, INT7, NL , std TM, tight timing
    cuts.AddCutCaloCalo("00010113","24466190sa01cc00000","3885500067032230000","0163103100b00010"); // INT7 NoNL
    cuts.AddCutCaloCalo("00010113","24466190sa01cc00000","3885501067032230000","0163103100b00010"); // INT7 TBNL
    cuts.AddCutCaloCalo("00010113","24466190sa01cc00000","3885511067032230000","0163103100b00010"); // INT7 NL11
    cuts.AddCutCaloCalo("00010113","24466190sa01cc00000","3885512067032230000","0163103100b00010"); // INT7 NL12
    cuts.AddCutCaloCalo("00010113","24466190sa01cc00000","3885521067032230000","0163103100b00010"); // INT7 NL21
    cuts.AddCutCaloCalo("00010113","24466190sa01cc00000","3885522067032230000","0163103100b00010"); // INT7 NL22
} else if (trainConfig == 22){ // timing Cut variation  std -50+30ns
  cuts.AddCutCaloCalo("00010113","24466190sa01cc00000","3885501017032230000","0163103100b00010"); //     -1000  +1000 ns
  cuts.AddCutCaloCalo("00010113","24466190sa01cc00000","3885501077032230000","0163103100b00010"); //     -30    +30   ns
  cuts.AddCutCaloCalo("00010113","24466190sa01cc00000","3885501097032230000","0163103100b00010"); //     -20    +25   ns
  cuts.AddCutCaloCalo("00010113","24466190sa01cc00000","38855010a7032230000","0163103100b00010"); //     -12.5  +13   ns
} else if (trainConfig == 23){ // track matching variation
  cuts.AddCutCaloCalo("00010113","24466190sa01cc00000","3885501060032230000","0163103100b00010"); //
  cuts.AddCutCaloCalo("00010113","24466190sa01cc00000","3885501061032230000","0163103100b00010"); //
  cuts.AddCutCaloCalo("00010113","24466190sa01cc00000","3885501066032230000","0163103100b00010"); //
} else if (trainConfig == 24){ // min nCells & M02 variation // std: min nCells = 1
  cuts.AddCutCaloCalo("00010113","24466190sa01cc00000","3885501067031230000","0163103100b00010"); //   min nCells = 1
  cuts.AddCutCaloCalo("00010113","24466190sa01cc00000","3885501067033230000","0163103100b00010"); //   min nCells = 3
} else if (trainConfig == 25){ // min energy variation std 0.7 GeV/c
  cuts.AddCutCaloCalo("00010113","24466190sa01cc00000","3885501067002230000","0163103100b00010"); //     0.1 GeV/c
  cuts.AddCutCaloCalo("00010113","24466190sa01cc00000","3885501067012230000","0163103100b00010"); //     0.5 GeV/c
  cuts.AddCutCaloCalo("00010113","24466190sa01cc00000","3885501067022230000","0163103100b00010"); //     0.6 GeV/c
  cuts.AddCutCaloCalo("00010113","24466190sa01cc00000","3885501067042230000","0163103100b00010"); //     0.8 GeV/c
  cuts.AddCutCaloCalo("00010113","24466190sa01cc00000","3885501067052230000","0163103100b00010"); //     0.9 GeV/c
} else if (trainConfig == 26){ // min nCells & M02 variation  // std: M02 max=0.7, min=0.1
  cuts.AddCutCaloCalo("00010113","24466190sa01cc00000","3885501067032210000","0163103100b00010"); //   max M02    = 1
  cuts.AddCutCaloCalo("00010113","24466190sa01cc00000","3885501067032240000","0163103100b00010"); //   max M02    = 0.4
} else if (trainConfig == 27){ // exotic cluster
  cuts.AddCutCaloCalo("00010113","24466190sa01cc00000","3885501067232230000","0163103100b00010"); // ExC 99
  cuts.AddCutCaloCalo("00010113","24466190sa01cc00000","3885501067532230000","0163103100b00010"); // ExC 97
  cuts.AddCutCaloCalo("00010113","24466190sa01cc00000","3885501067832230000","0163103100b00010"); // ExC 95
  cuts.AddCutCaloCalo("00010113","24466190sa01cc00000","3885501067b32230000","0163103100b00010"); // ExC 95 + TCard > 50
  cuts.AddCutCaloCalo("00010113","24466190sa01cc00000","3885501067e32230000","0163103100b00010"); // ExC 97 + TCard > 50

  // Variations of PHOS Part
  //Standard: "24466190sa01cc00000"
} else if (trainConfig == 30){ // timing Cut variation  std -30+30ns
  //                                  |
  cuts.AddCutCaloCalo("00010113","244661901a01cc00000","3885512067032220000","0163103100000010"); //1:     -1000  +1000 ns
  cuts.AddCutCaloCalo("00010113","244661905a01cc00000","3885512067032220000","0163103100000010"); //5:     -50    +50   ns
  cuts.AddCutCaloCalo("00010113","244661907a01cc00000","3885512067032220000","0163103100000010"); //7:     -30    +30   ns
  cuts.AddCutCaloCalo("00010113","244661909a01cc00000","3885512067032220000","0163103100000010"); //9:     -20    +25   ns
  cuts.AddCutCaloCalo("00010113","24466190aa01cc00000","3885512067032220000","0163103100000010"); //a:     -12.5  +13   ns
} else if (trainConfig == 31){ // Timing Efficiency Variations std s == 30ns
  //                                  |
  cuts.AddCutCaloCalo("00010113","24466190ra01cc00000","3885512067032220000","0163103100000010"); //r:     LowPt from MB, 30ns
  cuts.AddCutCaloCalo("00010113","24466190ta01cc00000","3885512067032220000","0163103100000010"); //t:     25ns
  cuts.AddCutCaloCalo("00010113","24466190ua01cc00000","3885512067032220000","0163103100000010"); //u:     50ns
} else if (trainConfig == 32){ // track matching variation std a == pt dependent
  //                                   |
  cuts.AddCutCaloCalo("00010113","24466190s001cc00000","3885512067032220000","0163103100000010"); //0
  cuts.AddCutCaloCalo("00010113","24466190s101cc00000","3885512067032220000","0163103100000010"); //1
  cuts.AddCutCaloCalo("00010113","24466190s401cc00000","3885512067032220000","0163103100000010"); //4
  cuts.AddCutCaloCalo("00010113","24466190s501cc00000","3885512067032220000","0163103100000010"); //5
  cuts.AddCutCaloCalo("00010113","24466190s601cc00000","3885512067032220000","0163103100000010"); //6
} else if (trainConfig == 33){ // min energy variation std 0.3 GeV/c
  //                                     |
  cuts.AddCutCaloCalo("00010113","24466190sa00cc00000","3885512067032220000","0163103100000010"); //0:     off
  cuts.AddCutCaloCalo("00010113","24466190sa09cc00000","3885512067032220000","0163103100000010"); //9:     0.1 GeV/c
  cuts.AddCutCaloCalo("00010113","24466190sa02cc00000","3885512067032220000","0163103100000010"); //2:     0.5 GeV/c
  cuts.AddCutCaloCalo("00010113","24466190sa03cc00000","3885512067032220000","0163103100000010"); //3:     0.6 GeV/c
  cuts.AddCutCaloCalo("00010113","24466190sa04cc00000","3885512067032220000","0163103100000010"); //4:     0.7 GeV/c
  cuts.AddCutCaloCalo("00010113","24466190sa05cc00000","3885512067032220000","0163103100000010"); //5:     0.8 GeV/c
} else if (trainConfig == 34){ // min nCells & M02 variation, std cc
  // std: min nCells = 2 >1GeV; M02 max=100, min=0.1
  //                                     |||
  cuts.AddCutCaloCalo("00010113","24466190sa011000000","3885512067032220000","0163103100000010"); //100:   min nCells = 1, minM02 off
  cuts.AddCutCaloCalo("00010113","24466190sa012200000","3885512067032220000","0163103100000010"); //220:   min nCells = 2, all E
  cuts.AddCutCaloCalo("00010113","24466190sa013200000","3885512067032220000","0163103100000010"); //320:   min nCells = 3, all E
  cuts.AddCutCaloCalo("00010113","24466190sa01dc00000","3885512067032220000","0163103100000010"); //dc0:   min nCells = 3, E>1GeV; minM02==0.1 off for E<1GeV
  cuts.AddCutCaloCalo("00010113","24466190sa01cd00000","3885512067032220000","0163103100000010"); //cd0:   min nCells = 2, E>1GeV; minM02==0.2 off for E<1GeV
  cuts.AddCutCaloCalo("00010113","24466190sa01cc70000","3885512067032220000","0163103100000010"); //cc7:   maxM02 == 1.3
  cuts.AddCutCaloCalo("00010113","24466190sa01cc80000","3885512067032220000","0163103100000010"); //cc8:   maxM02 == 2.5
} else {
    Error(Form("GammaCaloMix_%i",trainConfig), "wrong trainConfig variable no cuts have been specified for the configuration");
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
  TList *ClusterCutList2 = new TList();
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
  ClusterCutList2->SetOwner(kTRUE);
  AliCaloPhotonCuts **analysisClusterCuts2 = new AliCaloPhotonCuts*[numberOfCuts];
  MesonCutList->SetOwner(kTRUE);
  AliConversionMesonCuts **analysisMesonCuts = new AliConversionMesonCuts*[numberOfCuts];

  for(Int_t i = 0; i<numberOfCuts; i++){
    //create AliCaloTrackMatcher instance, if there is none present
    TString caloCutPos = cuts.GetClusterCut(i);
    TString caloCutPos2 = cuts.GetClusterCut2(i);
    caloCutPos.Resize(1);
    caloCutPos2.Resize(1);
    TString TrackMatcherName = Form("CaloTrackMatcher_%s_%s_%i",caloCutPos.Data(),caloCutPos2.Data(),trackMatcherRunningMode);
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
    mcInputNamePi0 = Form("Pi0_%s%s_%s", mcName.Data(), mcNameAdd.Data(), energy.Data() );
    mcInputNameEta = Form("Eta_%s%s_%s", mcName.Data(), mcNameAdd.Data(), energy.Data() );

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
    analysisEventCuts[i]->SetLightOutput(enableLightOutput);
    // analysisEventCuts[i]->SetUseSphericityTrue(kTRUE);
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

    analysisClusterCuts2[i] = new AliCaloPhotonCuts(isMC);
    analysisClusterCuts2[i]->SetHistoToModifyAcceptance(histoAcc);
    analysisClusterCuts2[i]->SetV0ReaderName(V0ReaderName);
    analysisClusterCuts2[i]->SetCorrectionTaskSetting(corrTaskSetting);
    analysisClusterCuts2[i]->SetCaloTrackMatcherName(TrackMatcherName);
    analysisClusterCuts2[i]->SetLightOutput(enableLightOutput);
    analysisClusterCuts2[i]->InitializeCutsFromCutString((cuts.GetClusterCut2(i)).Data());
    ClusterCutList2->Add(analysisClusterCuts2[i]);
    analysisClusterCuts2[i]->SetExtendedMatchAndQA(enableExtMatchAndQA);
    analysisClusterCuts2[i]->SetFillCutHistograms("");

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
  task->SetCaloCaloCutLists(numberOfCuts,ClusterCutList,ClusterCutList2);
  task->SetMesonCutList(numberOfCuts,MesonCutList);
  task->SetDoMesonAnalysis(kTRUE);
  task->SetCorrectionTaskSetting(corrTaskSetting);
  task->SetDoMesonQA(enableQAMesonTask); //Attention new switch for Pi0 QA
  task->SetDoClusterQA(enableQAClusterTask);  //Attention new switch small for Cluster QA
  task->SetDoTHnSparse(enableTHnSparse);
  task->SetProduceTreeEOverP(doTreeEOverP);
  task->SetEnableSortingOfMCClusLabels(enableSortingMCLabels);
  if(enableExtMatchAndQA > 1){ task->SetPlotHistsExtQA(kTRUE);}
  task->SetLocalDebugFlag(localDebugFlag);

  //connect containers
  AliAnalysisDataContainer *coutput =
    mgr->CreateContainer(!(corrTaskSetting.CompareTo("")) ? Form("GammaCaloMix_%i",trainConfig) : Form("GammaCaloMix_%i_%s",trainConfig,corrTaskSetting.Data()), TList::Class(),
              AliAnalysisManager::kOutputContainer,Form("GammaCaloMix_%i.root",trainConfig));

  mgr->AddTask(task);
  mgr->ConnectInput(task,0,cinput);
  mgr->ConnectOutput(task,1,coutput);

  Int_t nContainer = 2;
  for(Int_t i = 0; i<numberOfCuts; i++){
      if(enableQAMesonTask==5){
          mgr->ConnectOutput(task,nContainer,mgr->CreateContainer(Form("%s_%s_%s ClusterTimingEff",(cuts.GetEventCut(i)).Data(),(cuts.GetClusterCut(i)).Data(),(cuts.GetMesonCut(i)).Data()), TTree::Class(), AliAnalysisManager::kOutputContainer, Form("GammaCaloMix_%i.root",trainConfig)) );
          nContainer++;
      }
  }


  return;

}
