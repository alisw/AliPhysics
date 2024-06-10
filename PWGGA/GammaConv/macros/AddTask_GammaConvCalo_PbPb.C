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
//PbPb together with all supporting classes
//***************************************************************************************

//***************************************************************************************
//main function
//***************************************************************************************
void AddTask_GammaConvCalo_PbPb(
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
  TString   settingMaxFacPtHard           = "3.",     // maximum factor between hardest jet and ptHard generated
  Int_t     debugLevel                    = 0,        // introducing debug levels for grid running
  // settings for weights
  // FPTW:fileNamePtWeights, FMUW:fileNameMultWeights, FCEF:fileNameCentFlattening,   FMAW:fileNameMatBudWeights, FEPC:fileNamedEdxPostCalib, separate with ;
  // Material Budget Weights file for Run 2
  // FMAW:alien:///alice/cern.ch/user/a/amarin//MBW/MCInputFileMaterialBudgetWeightsLHC16_Pythia_00010103_0d000009266300008850404000_date181214.root
  TString   fileNameExternalInputs        = "",
  Bool_t    doWeightingPart               = kFALSE,   // enable Weighting
  TString   generatorName                 = "DPMJET", // generator Name
  Bool_t    enableMultiplicityWeighting   = kFALSE,   //
  Int_t     enableMatBudWeightsPi0        = 0,        // 1 = three radial bins, 2 = 10 radial bins (2 is the default when using weights)
  Int_t     enableElecDeDxPostCalibration = 0,        // 0 = off, 1 = as function of TPC clusters, 2 = as function of convR. Option 2 is available only if FEPC is given.
  TString   periodNameAnchor              = "",       //
  Bool_t    enableFlattening              = kFALSE,   // switch on centrality flattening for LHC11h
  // special settings
  Bool_t    enableSortingMCLabels         = kTRUE,    // enable sorting for MC cluster labels
  Bool_t    enableTreeConvGammaShape      = kFALSE,   // enable additional tree for conversion properties for clusters
  Bool_t    doPrimaryTrackMatching        = kTRUE,    // enable basic track matching for all primary tracks to cluster
  Int_t     headerSelectionInt            = 0,        // 1 pi0 header, 2 eta header, 3 both (only for "named" boxes)
  Bool_t    enableHeaderOverlap           = kTRUE,    // allow overlapping header for the clusters
  // subwagon config
  TString   additionalTrainConfig         = "0"       // additional counter for trainconfig
) {

  AliCutHandlerPCM cuts;

  TString fileNamePtWeights     = cuts.GetSpecialFileNameFromString (fileNameExternalInputs, "FPTW:");
  TString fileNameMultWeights   = cuts.GetSpecialFileNameFromString (fileNameExternalInputs, "FMUW:");
  TString fileNameCentFlattening= cuts.GetSpecialFileNameFromString (fileNameExternalInputs, "FCEF:");
  TString fileNameMatBudWeights = cuts.GetSpecialFileNameFromString (fileNameExternalInputs, "FMAW:");
  TString fileNamedEdxPostCalib = cuts.GetSpecialFileNameFromString (fileNameExternalInputs, "FEPC:");


  TString addTaskName                 = "AddTask_GammaConvCalo_PbPb";
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
  TString strdoTreeClusterShowerShape = cuts.GetSpecialSettingFromAddConfig(additionalTrainConfig, "INVMASSCLUSTree", "", addTaskName);
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
  if(rmaxFacPtHardSetting->GetEntries()<1){cout << "ERROR: AddTask_GammaConvCalo_PbPb during parsing of settingMaxFacPtHard String '" << settingMaxFacPtHard.Data() << "'" << endl; return;}
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

  Int_t isHeavyIon = 1;

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
  TString cutnumberEvent = "10000003";
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
  if(enableLightOutput >= 20){  // eta
    enableLightOutput -= 20;
    task->SetPi0EtaSwitch(2);
  } else if(enableLightOutput >= 10){  // pi0
    enableLightOutput -= 10;
    task->SetPi0EtaSwitch(1);
  } else {
    task->SetPi0EtaSwitch(0);
  }
  if (enableLightOutput >= 1 && enableLightOutput != 4) task->SetLightOutput(1);
  else if (enableLightOutput == 4) task->SetLightOutput(2);
  if (enableLightOutput == 4 || enableLightOutput == 5) task->SetECalibOutput(1);
  if (enableLightOutput == 6) task->SetECalibOutput(2);
  task->SetDoPrimaryTrackMatching(doPrimaryTrackMatching);
  task->SetTrackMatcherRunningMode(trackMatcherRunningMode);
  if(trainConfig >= 950 && trainConfig <= 1000) task->SetDoHBTHistoOutput(kTRUE);

  // cluster cuts
  // 0 "ClusterType",  1 "EtaMin", 2 "EtaMax", 3 "PhiMin", 4 "PhiMax", 5 "DistanceToBadChannel", 6 "Timing", 7 "TrackMatching", 8 "ExoticCell",
  // 9 "MinEnergy", 10 "MinNCells", 11 "MinM02", 12 "MaxM02", 13 "MinMaxM20", 14 "RecConv", 15 "MaximumDispersion", 16 "NLM"

  //****************************************************************************************************
  // EMCal 2.76 TeV Pb-Pb  LHC11h & LHC10h MB
  //****************************************************************************************************
  if (trainConfig == 1){ // EMCAL clusters
    cuts.AddCutPCMCalo("60100013","00200009297002008250400000","1111100053032230000","0163103100000010"); // 0-5%
    cuts.AddCutPCMCalo("61200013","00200009297002008250400000","1111100053032230000","0163103100000010"); // 5-10%
    cuts.AddCutPCMCalo("50100013","00200009297002008250400000","1111100053032230000","0163103100000010"); // 0-10%
    cuts.AddCutPCMCalo("52400013","00200009297002008250400000","1111100053032230000","0163103100000010"); // 20-40%
    cuts.AddCutPCMCalo("52500013","00200009297002008250400000","1111100053032230000","0163103100000010"); // 20-50%
  } else if (trainConfig == 2){ // EMCAL clusters no timing
    cuts.AddCutPCMCalo("60100013","00200009297002008250400000","1111100003032230000","0163103100000010"); // 0-5%
    cuts.AddCutPCMCalo("61200013","00200009297002008250400000","1111100003032230000","0163103100000010"); // 5-10%
    cuts.AddCutPCMCalo("50100013","00200009297002008250400000","1111100003032230000","0163103100000010"); // 0-10%
    cuts.AddCutPCMCalo("52400013","00200009297002008250400000","1111100003032230000","0163103100000010"); // 20-40%
    cuts.AddCutPCMCalo("52500013","00200009297002008250400000","1111100003032230000","0163103100000010"); // 20-50%
  } else if (trainConfig == 3){ // EMCAL clusters
    cuts.AddCutPCMCalo("60100013","00200009297002008250400000","1111100053032230000","0163103100000010"); // 0-5%
    cuts.AddCutPCMCalo("61200013","00200009297002008250400000","1111100053032230000","0163103100000010"); // 5-10%
    cuts.AddCutPCMCalo("50100013","00200009297002008250400000","1111100053032230000","0163103100000010"); // 0-10%
    cuts.AddCutPCMCalo("51200013","00200009297002008250400000","1111100053032230000","0163103100000010"); // 10-20%
    cuts.AddCutPCMCalo("52400013","00200009297002008250400000","1111100053032230000","0163103100000010"); // 20-40%
  } else if (trainConfig == 4){ // EMCAL clusters no timing
    cuts.AddCutPCMCalo("60100013","00200009297002008250400000","1111100003032230000","0163103100000010"); // 0-5%
    cuts.AddCutPCMCalo("61200013","00200009297002008250400000","1111100003032230000","0163103100000010"); // 5-10%
    cuts.AddCutPCMCalo("50100013","00200009297002008250400000","1111100003032230000","0163103100000010"); // 0-10%
    cuts.AddCutPCMCalo("51200013","00200009297002008250400000","1111100003032230000","0163103100000010"); // 10-20%
    cuts.AddCutPCMCalo("52400013","00200009297002008250400000","1111100003032230000","0163103100000010"); // 20-40%
  } else if (trainConfig == 5){ // EMCAL clusters
    cuts.AddCutPCMCalo("54600013","00200009297002008250400000","1111100053032230000","0163103100000010"); // 40-60%
    cuts.AddCutPCMCalo("56800013","00200009297002008250400000","1111100053032230000","0163103100000010"); // 60-80%
    cuts.AddCutPCMCalo("52600013","00200009297002008250400000","1111100053032230000","0163103100000010"); // 20-60%
    cuts.AddCutPCMCalo("54800013","00200009297002008250400000","1111100053032230000","0163103100000010"); // 40-80%
    cuts.AddCutPCMCalo("52500013","00200009297002008250400000","1111100053032230000","0163103100000010"); // 20-50%
  } else if (trainConfig == 6){ // EMCAL clusters no timing
    cuts.AddCutPCMCalo("54600013","00200009297002008250400000","1111100003032230000","0163103100000010"); // 40-60%
    cuts.AddCutPCMCalo("56800013","00200009297002008250400000","1111100003032230000","0163103100000010"); // 60-80%
    cuts.AddCutPCMCalo("52600013","00200009297002008250400000","1111100003032230000","0163103100000010"); // 20-60%
    cuts.AddCutPCMCalo("54800013","00200009297002008250400000","1111100003032230000","0163103100000010"); // 40-80%
    cuts.AddCutPCMCalo("52500013","00200009297002008250400000","1111100003032230000","0163103100000010"); // 20-50%

  } else if (trainConfig == 7){ // EMCAL clusters
    cuts.AddCutPCMCalo("60100013","00200009297002008250400000","1111100053032230000","0163103100000010"); // 0-5%
    cuts.AddCutPCMCalo("61200013","00200009297002008250400000","1111100053032230000","0163103100000010"); // 5-10%
    cuts.AddCutPCMCalo("50100013","00200009297002008250400000","1111100053032230000","0163103100000010"); // 0-10%
  } else if (trainConfig == 8){ // EMCAL clusters no timing
    cuts.AddCutPCMCalo("60100013","00200009297002008250400000","1111100003032230000","0163103100000010"); // 0-5%
    cuts.AddCutPCMCalo("61200013","00200009297002008250400000","1111100003032230000","0163103100000010"); // 5-10%
    cuts.AddCutPCMCalo("50100013","00200009297002008250400000","1111100003032230000","0163103100000010"); // 0-10%
  } else if (trainConfig == 9){ // EMCAL clusters
    cuts.AddCutPCMCalo("51200013","00200009297002008250400000","1111100053032230000","0163103100000010"); // 10-20%
    cuts.AddCutPCMCalo("52400013","00200009297002008250400000","1111100053032230000","0163103100000010"); // 20-40%
    cuts.AddCutPCMCalo("52500013","00200009297002008250400000","1111100053032230000","0163103100000010"); // 20-50%
  } else if (trainConfig == 10){ // EMCAL clusters no timing
    cuts.AddCutPCMCalo("51200013","00200009297002008250400000","1111100003032230000","0163103100000010"); // 10-20%
    cuts.AddCutPCMCalo("52400013","00200009297002008250400000","1111100003032230000","0163103100000010"); // 20-40%
    cuts.AddCutPCMCalo("52500013","00200009297002008250400000","1111100003032230000","0163103100000010"); // 20-50%
  } else if (trainConfig == 11){ // EMCAL clusters
    cuts.AddCutPCMCalo("54600013","00200009297002008250400000","1111100053032230000","0163103100000010"); // 40-60%
    cuts.AddCutPCMCalo("56800013","00200009297002008250400000","1111100053032230000","0163103100000010"); // 60-80%
  } else if (trainConfig == 12){ // EMCAL clusters no timing
    cuts.AddCutPCMCalo("54600013","00200009297002008250400000","1111100003032230000","0163103100000010"); // 40-60%
    cuts.AddCutPCMCalo("56800013","00200009297002008250400000","1111100003032230000","0163103100000010"); // 60-80%


  } else if (trainConfig == 13){ // EMCAL clusters
    cuts.AddCutPCMCalo("50900013","00200009297002008250400000","1111100053032230000","0163103100000010"); // 0-90%
  } else if (trainConfig == 14){ // EMCAL clusters
    cuts.AddCutPCMCalo("50900013","00200009297002008250400000","1111100003032230000","0163103100000010"); // 0-90%

    // EMC triggers LHC15o
  } else if (trainConfig == 20){ //MB
    cuts.AddCutPCMCalo("10110013","00200009297002008250400000","1111100051032230000","0163103100000050"); //
    cuts.AddCutPCMCalo("11310013","00200009297002008250400000","1111100051032230000","0163103100000050"); //
    cuts.AddCutPCMCalo("13510013","00200009297002008250400000","1111100051032230000","0163103100000050"); //
    cuts.AddCutPCMCalo("15010013","00200009297002008250400000","1111100051032230000","0163103100000050"); //
  } else if (trainConfig == 21){ //EG2
    cuts.AddCutPCMCalo("10185013","00200009297002008250400000","1111100051032230000","0163103100000050"); //
    cuts.AddCutPCMCalo("11385013","00200009297002008250400000","1111100051032230000","0163103100000050"); //
    cuts.AddCutPCMCalo("13585013","00200009297002008250400000","1111100051032230000","0163103100000050"); //
    cuts.AddCutPCMCalo("15085013","00200009297002008250400000","1111100051032230000","0163103100000050"); //
  } else if (trainConfig == 22){ //EG1
    cuts.AddCutPCMCalo("10183013","00200009297002008250400000","1111100051032230000","0163103100000050"); //
    cuts.AddCutPCMCalo("11383013","00200009297002008250400000","1111100051032230000","0163103100000050"); //
    cuts.AddCutPCMCalo("13583013","00200009297002008250400000","1111100051032230000","0163103100000050"); //
    cuts.AddCutPCMCalo("15083013","00200009297002008250400000","1111100051032230000","0163103100000050"); //
  } else if (trainConfig == 23){ //EJ2
    cuts.AddCutPCMCalo("10195013","00200009297002008250400000","1111100051032230000","0163103100000050"); //
    cuts.AddCutPCMCalo("11395013","00200009297002008250400000","1111100051032230000","0163103100000050"); //
    cuts.AddCutPCMCalo("13595013","00200009297002008250400000","1111100051032230000","0163103100000050"); //
    cuts.AddCutPCMCalo("15095013","00200009297002008250400000","1111100051032230000","0163103100000050"); //
  } else if (trainConfig == 24){ //EJ1
    cuts.AddCutPCMCalo("10193013","00200009297002008250400000","1111100051032230000","0163103100000050"); //
    cuts.AddCutPCMCalo("11393013","00200009297002008250400000","1111100051032230000","0163103100000050"); //
    cuts.AddCutPCMCalo("13593013","00200009297002008250400000","1111100051032230000","0163103100000050"); //
    cuts.AddCutPCMCalo("15093013","00200009297002008250400000","1111100051032230000","0163103100000050"); //
    // DCAL triggers LHC15o
  } else if (trainConfig == 25){ //MB
    cuts.AddCutPCMCalo("10110013","00200009297002008250400000","3885500017032220000","0163103100000050"); //
    cuts.AddCutPCMCalo("11310013","00200009297002008250400000","3885500017032220000","0163103100000050"); //
    cuts.AddCutPCMCalo("13510013","00200009297002008250400000","3885500017032220000","0163103100000050"); //
    cuts.AddCutPCMCalo("15010013","00200009297002008250400000","3885500017032220000","0163103100000050"); //
  } else if (trainConfig == 26){ //EG2
    cuts.AddCutPCMCalo("10185013","00200009297002008250400000","3885500017032220000","0163103100000050"); //
    cuts.AddCutPCMCalo("11385013","00200009297002008250400000","3885500017032220000","0163103100000050"); //
    cuts.AddCutPCMCalo("13585013","00200009297002008250400000","3885500017032220000","0163103100000050"); //
    cuts.AddCutPCMCalo("15085013","00200009297002008250400000","3885500017032220000","0163103100000050"); //
  } else if (trainConfig == 27){ //EG1
    cuts.AddCutPCMCalo("10183013","00200009297002008250400000","3885500017032220000","0163103100000050"); //
    cuts.AddCutPCMCalo("11383013","00200009297002008250400000","3885500017032220000","0163103100000050"); //
    cuts.AddCutPCMCalo("13583013","00200009297002008250400000","3885500017032220000","0163103100000050"); //
    cuts.AddCutPCMCalo("15083013","00200009297002008250400000","3885500017032220000","0163103100000050"); //
  } else if (trainConfig == 28){ //EJ2
    cuts.AddCutPCMCalo("10195013","00200009297002008250400000","3885500017032220000","0163103100000050"); //
    cuts.AddCutPCMCalo("11395013","00200009297002008250400000","3885500017032220000","0163103100000050"); //
    cuts.AddCutPCMCalo("13595013","00200009297002008250400000","3885500017032220000","0163103100000050"); //
    cuts.AddCutPCMCalo("15095013","00200009297002008250400000","3885500017032220000","0163103100000050"); //
  } else if (trainConfig == 29){ //EJ1
    cuts.AddCutPCMCalo("10193013","00200009297002008250400000","3885500017032220000","0163103100000050"); //
    cuts.AddCutPCMCalo("11393013","00200009297002008250400000","3885500017032220000","0163103100000050"); //
    cuts.AddCutPCMCalo("13593013","00200009297002008250400000","3885500017032220000","0163103100000050"); //
    cuts.AddCutPCMCalo("15093013","00200009297002008250400000","3885500017032220000","0163103100000050"); //

  //****************************************************************************************************
  // EMCal 2.76TeV Pb-Pb for LHC11h EMC trigger
  //****************************************************************************************************
  } else if (trainConfig == 30){ // EMCAL clusters
    cuts.AddCutPCMCalo("50980013","00200009297002008250400000","1111100053032230000","0163103100000010"); // 0-90%
    cuts.AddCutPCMCalo("50180013","00200009297002008250400000","1111100053032230000","0163103100000010"); // 0-10%
    cuts.AddCutPCMCalo("51280013","00200009297002008250400000","1111100053032230000","0163103100000010"); // 10-20%
    cuts.AddCutPCMCalo("52580013","00200009297002008250400000","1111100053032230000","0163103100000010"); // 20-50%
    cuts.AddCutPCMCalo("55880013","00200009297002008250400000","1111100053032230000","0163103100000010"); // 50-80%
  } else if (trainConfig == 31){ // EMCAL clusters no timing cut
    cuts.AddCutPCMCalo("50980013","00200009297002008250400000","1111100003032230000","0163103100000010"); // 0-90%
    cuts.AddCutPCMCalo("50180013","00200009297002008250400000","1111100003032230000","0163103100000010"); // 0-10%
    cuts.AddCutPCMCalo("51280013","00200009297002008250400000","1111100003032230000","0163103100000010"); // 10-20%
    cuts.AddCutPCMCalo("52580013","00200009297002008250400000","1111100003032230000","0163103100000010"); // 20-50%
    cuts.AddCutPCMCalo("55880013","00200009297002008250400000","1111100003032230000","0163103100000010"); // 50-80%
  } else if (trainConfig == 32){ // EMCAL clusters
    cuts.AddCutPCMCalo("50980013","00200009297002008250400000","1111100053032230000","0163103100000010"); // 0-90%
  } else if (trainConfig == 33){ // EMCAL clusters no timing cut
    cuts.AddCutPCMCalo("50980013","00200009297002008250400000","1111100003032230000","0163103100000010"); // 0-90%
  } else if (trainConfig == 34){ // EMCAL clusters
    cuts.AddCutPCMCalo("50100013","00200009297002008250400000","1111102053032230000","0163103100000010"); // 0-10
    cuts.AddCutPCMCalo("50100013","00200009297002008250400000","1111171053032230000","0163103100000010"); // 0-10
  } else if (trainConfig == 35){ // EMCAL clusters
    cuts.AddCutPCMCalo("50100013","00200009297002008250400000","11111020530a2230000","0163103100000010"); // 0-10 // using cluster cuts of Astrid
    cuts.AddCutPCMCalo("50100013","00200009297002008250400000","11111710530a2230000","0163103100000010"); // 0-10
    cuts.AddCutPCMCalo("50100013","00200009297002008250400000","1111171053032230000","0163103100000010");
  } else if (trainConfig == 36){ // EMCAL clusters
    cuts.AddCutPCMCalo("52500013","00200009297002008250400000","1111102053032230000","0163103100000010"); // 20-50
    cuts.AddCutPCMCalo("52500013","00200009297002008250400000","1111172053032230000","0163103100000010"); // 20-50
  } else if (trainConfig == 37){ // EMCAL clusters
    cuts.AddCutPCMCalo("52500013","00200009297002008250400000","11111020530a2230000","0163103100000010"); // 20-50
    cuts.AddCutPCMCalo("52500013","00200009297002008250400000","11111720530a2230000","0163103100000010"); // 20-50
  } else if (trainConfig == 38){ // EMCAL clusters - added signals
    cuts.AddCutPCMCalo("50100023","00200009297002008250400000","1111171053032230000","0163103100000010"); // 0-10
    cuts.AddCutPCMCalo("50100023","00200009297002008250400000","11111710530a2230000","0163103100000010"); // 0-10
    cuts.AddCutPCMCalo("50100023","00200009297002008250400000","11111710530b2230000","0163103100000010"); // 0-10
  } else if (trainConfig == 39){ // EMCAL clusters - added signals
    cuts.AddCutPCMCalo("50100023","00200009297002008250400000","1111171053032230000","0163103100000010"); // 0-10
    cuts.AddCutPCMCalo("50100023","00200009297002008250400000","11111710530a2230000","0163103100000010"); // 0-10
    cuts.AddCutPCMCalo("50100023","00200009297002008250400000","11111710530b2230000","0163103100000010"); // 0-10
  } else if (trainConfig == 40){ // EMCAL clusters - added signals
    cuts.AddCutPCMCalo("52500023","00200009297002008250400000","1111172053032230000","0163103100000010"); // 20-50
    cuts.AddCutPCMCalo("52500023","00200009297002008250400000","11111720530a2230000","0163103100000010"); // 20-50
    cuts.AddCutPCMCalo("52500023","00200009297002008250400000","11111720530b2230000","0163103100000010"); // 20-50
  } else if (trainConfig == 41){ // EMCAL clusters - added signals
    cuts.AddCutPCMCalo("52500023","00200009297002008250400000","1111172053032230000","0163103100000010"); // 20-50
    cuts.AddCutPCMCalo("52500023","00200009297002008250400000","11111720530a2230000","0163103100000010"); // 20-50
    cuts.AddCutPCMCalo("52500023","00200009297002008250400000","11111720530b2230000","0163103100000010"); // 20-50
  } else if (trainConfig == 42){ // EMCAL clusters all headers
    cuts.AddCutPCMCalo("50100003","00200009297002008250400000","11111020530a2230000","0163103100000010"); // 0-10
    cuts.AddCutPCMCalo("50100003","00200009297002008250400000","11111710530a2230000","0163103100000010"); // 0-10
    cuts.AddCutPCMCalo("50100003","00200009297002008250400000","1111171053032230000","0163103100000010");
  } else if (trainConfig == 43){ // EMCAL clusters added signals pi0 header forseen
    cuts.AddCutPCMCalo("50100023","00200009297002008250400000","11111020530a2230000","0163103100000010"); // 0-10 // using cluster cuts of Astrid
    cuts.AddCutPCMCalo("50100023","00200009297002008250400000","11111710530a2230000","0163103100000010"); // 0-10
    cuts.AddCutPCMCalo("50100023","00200009297002008250400000","1111171053032230000","0163103100000010");
  } else if (trainConfig == 44){ // EMCAL clusters added signals eta header forseen
    cuts.AddCutPCMCalo("50100023","00200009297002008250400000","11111020530a2230000","0163103100000010"); // 0-10 // using cluster cuts of Astrid
    cuts.AddCutPCMCalo("50100023","00200009297002008250400000","11111710530a2230000","0163103100000010"); // 0-10
    cuts.AddCutPCMCalo("50100023","00200009297002008250400000","1111171053032230000","0163103100000010");
  } else if (trainConfig == 45){ // EMCAL clusters added signals other header forseen
    cuts.AddCutPCMCalo("50100023","00200009297002008250400000","11111020530a2230000","0163103100000010"); // 0-10 // using cluster cuts of Astrid
    cuts.AddCutPCMCalo("50100023","00200009297002008250400000","11111710530a2230000","0163103100000010"); // 0-10
    cuts.AddCutPCMCalo("50100023","00200009297002008250400000","1111171053032230000","0163103100000010");
  } else if (trainConfig == 46){ // EMCAL clusters added signals other header forseen
    cuts.AddCutPCMCalo("50100023","00200009297002008250400000","11111020530a2230000","0163103100000010"); // 0-10 // using cluster cuts of Astrid
    cuts.AddCutPCMCalo("50100023","00200009297002008250400000","11111710530a2230000","0163103100000010"); // 0-10
    cuts.AddCutPCMCalo("50100023","00200009297002008250400000","1111171053032230000","0163103100000010");
  } else if (trainConfig == 47){ // EMCAL clusters V0 Cent
    cuts.AddCutPCMCalo("10100013","00200009297002008250400000","11111020530a2230000","0163103100000010"); // 0-10 // using cluster cuts of Astrid
    cuts.AddCutPCMCalo("10100013","00200009297002008250400000","11111710530a2230000","0163103100000010"); // 0-10
    cuts.AddCutPCMCalo("10100013","00200009297002008250400000","1111171053032230000","0163103100000010");

  //****************************************************************************************************
  // PHOS 2.76TeV Pb-Pb LHC10h & LHC11h
  //****************************************************************************************************
  } else if (trainConfig == 101){ // PHOS clusters
    cuts.AddCutPCMCalo("60100013","00200009297002008250400000","2444400042033200000","0163103100000010"); // 0-5%
    cuts.AddCutPCMCalo("61200013","00200009297002008250400000","2444400042033200000","0163103100000010"); // 5-10%
    cuts.AddCutPCMCalo("50100013","00200009297002008250400000","2444400042033200000","0163103100000010"); // 0-10%
    cuts.AddCutPCMCalo("52400013","00200009297002008250400000","2444400042033200000","0163103100000010"); // 20-40%
    cuts.AddCutPCMCalo("52500013","00200009297002008250400000","2444400042033200000","0163103100000010"); // 20-50%
  } else if (trainConfig == 102){ // PHOS clusters
    cuts.AddCutPCMCalo("60100013","00200009297002008250400000","2444400042033200000","0163103100000010"); // 0-5%
    cuts.AddCutPCMCalo("61200013","00200009297002008250400000","2444400042033200000","0163103100000010"); // 5-10%
    cuts.AddCutPCMCalo("50100013","00200009297002008250400000","2444400042033200000","0163103100000010"); // 0-10%
    cuts.AddCutPCMCalo("51200013","00200009297002008250400000","2444400042033200000","0163103100000010"); // 10-20%
    cuts.AddCutPCMCalo("52400013","00200009297002008250400000","2444400042033200000","0163103100000010"); // 20-40%
  } else if (trainConfig == 103){ // PHOS clusters
    cuts.AddCutPCMCalo("54600013","00200009297002008250400000","2444400042033200000","0163103100000010"); // 40-60%
    cuts.AddCutPCMCalo("56800013","00200009297002008250400000","2444400042033200000","0163103100000010"); // 60-80%
    cuts.AddCutPCMCalo("52600013","00200009297002008250400000","2444400042033200000","0163103100000010"); // 20-60%
    cuts.AddCutPCMCalo("54800013","00200009297002008250400000","2444400042033200000","0163103100000010"); // 40-80%
    cuts.AddCutPCMCalo("52500013","00200009297002008250400000","2444400042033200000","0163103100000010"); // 20-50%
  } else if (trainConfig == 104){ // PHOS clusters no timing
    cuts.AddCutPCMCalo("60100013","00200009297002008250400000","2444400002033200000","0163103100000010"); // 0-5%
    cuts.AddCutPCMCalo("61200013","00200009297002008250400000","2444400002033200000","0163103100000010"); // 5-10%
    cuts.AddCutPCMCalo("50100013","00200009297002008250400000","2444400002033200000","0163103100000010"); // 0-10%
    cuts.AddCutPCMCalo("52400013","00200009297002008250400000","2444400002033200000","0163103100000010"); // 20-40%
    cuts.AddCutPCMCalo("52500013","00200009297002008250400000","2444400002033200000","0163103100000010"); // 20-50%
  } else if (trainConfig == 105){ // PHOS clusters no timing
    cuts.AddCutPCMCalo("60100013","00200009297002008250400000","2444400002033200000","0163103100000010"); // 0-5%
    cuts.AddCutPCMCalo("61200013","00200009297002008250400000","2444400002033200000","0163103100000010"); // 5-10%
    cuts.AddCutPCMCalo("50100013","00200009297002008250400000","2444400002033200000","0163103100000010"); // 0-10%
    cuts.AddCutPCMCalo("51200013","00200009297002008250400000","2444400002033200000","0163103100000010"); // 10-20%
    cuts.AddCutPCMCalo("52400013","00200009297002008250400000","2444400002033200000","0163103100000010"); // 20-40%
  } else if (trainConfig == 106){ // PHOS clusters no timing
    cuts.AddCutPCMCalo("54600013","00200009297002008250400000","2444400002033200000","0163103100000010"); // 40-60%
    cuts.AddCutPCMCalo("56800013","00200009297002008250400000","2444400002033200000","0163103100000010"); // 60-80%
    cuts.AddCutPCMCalo("52600013","00200009297002008250400000","2444400002033200000","0163103100000010"); // 20-60%
    cuts.AddCutPCMCalo("54800013","00200009297002008250400000","2444400002033200000","0163103100000010"); // 40-80%
    cuts.AddCutPCMCalo("52500013","00200009297002008250400000","2444400002033200000","0163103100000010"); // 20-50%
  } else if (trainConfig == 107){ // PHOS clusters no timing
    cuts.AddCutPCMCalo("50900013","00200009297002008250400000","2444400002033200000","0163103100000010"); // 0-90%

  } else if (trainConfig == 107){ // EMCAL clusters
    cuts.AddCutPCMCalo("50910013","00200009f9730000dge0404000","11111010500a2220000","01331031000000d0"); // 0-90%
    cuts.AddCutPCMCalo("50110013","00200009f9730000dge0404000","11111010500a2220000","01331061000000d0"); //
    cuts.AddCutPCMCalo("51310013","00200009f9730000dge0404000","11111010500b2220000","01331061000000d0"); //
    cuts.AddCutPCMCalo("53510013","00200009f9730000dge0404000","1111101050032220000","01331031000000d0"); //
    cuts.AddCutPCMCalo("55910013","00200009f9730000dge0404000","1111101050032220000","01331031000000d0"); //
  } else if (trainConfig == 108){ // EMCAL clusters no timing cut
    cuts.AddCutPCMCalo("50980013","00200009f9730000dge0404000","11111010500a2220000","01331031000000d0"); // 0-90%
    cuts.AddCutPCMCalo("50180013","00200009f9730000dge0404000","11111010500a2220000","01331061000000d0"); //
    cuts.AddCutPCMCalo("51380013","00200009f9730000dge0404000","11111010500b2220000","01331061000000d0"); //
    cuts.AddCutPCMCalo("53580013","00200009f9730000dge0404000","1111101050032220000","01331031000000d0"); //
    cuts.AddCutPCMCalo("55980013","00200009f9730000dge0404000","1111101050032220000","01331031000000d0"); //
  } else if (trainConfig == 109){ // EMCAL clusters
    cuts.AddCutPCMCalo("50900013","00200009f9730000dge0404000","11111010500a2220000","01331031000000d0"); // 0-90%
    cuts.AddCutPCMCalo("50100013","00200009f9730000dge0404000","11111010500a2220000","01331061000000d0"); //
    cuts.AddCutPCMCalo("51300013","00200009f9730000dge0404000","11111010500b2220000","01331061000000d0"); //
    cuts.AddCutPCMCalo("53500013","00200009f9730000dge0404000","1111101050032220000","01331031000000d0"); //
    cuts.AddCutPCMCalo("55900013","00200009f9730000dge0404000","1111101050032220000","01331031000000d0"); //
  } else if (trainConfig == 110){ // EMCAL clusters
    cuts.AddCutPCMCalo("50900013","00200009f9730000dge0404000","11111010500a2220000","01331031000000d0"); // 0-90%
    cuts.AddCutPCMCalo("50100013","00200009f9730000dge0404000","11111010500a2220000","01331061000000d0"); //
    cuts.AddCutPCMCalo("51300013","00200009f9730000dge0404000","11111010500b2220000","01331061000000d0"); //
    cuts.AddCutPCMCalo("53500013","00200009f9730000dge0404000","1111101050032220000","01331031000000d0"); //
    cuts.AddCutPCMCalo("55900013","00200009f9730000dge0404000","1111101050032220000","01331031000000d0"); //
  } else if (trainConfig == 111){ // EMCAL clusters
    cuts.AddCutPCMCalo("50900013","00200009f9730000dge0404000","11111010500a2220000","01331031000000d0"); // 0-90%
    cuts.AddCutPCMCalo("50100013","00200009f9730000dge0404000","11111010500a2220000","01331061000000d0"); //
    cuts.AddCutPCMCalo("51300013","00200009f9730000dge0404000","11111010500b2220000","01331061000000d0"); //
    cuts.AddCutPCMCalo("53500013","00200009f9730000dge0404000","1111101050032220000","01331031000000d0"); //
    cuts.AddCutPCMCalo("55900013","00200009f9730000dge0404000","1111101050032220000","01331031000000d0"); //
    
  //****************************************************************************************************
  // EMCal 5TeV Pb-Pb LHC15o
  //****************************************************************************************************
  } else if (trainConfig == 201){ // EMCAL clusters central
    cuts.AddCutPCMCalo("10110013","00200009327000008250400000","1111100057032230000","0163103100000010"); // 0-10
    cuts.AddCutPCMCalo("30110013","00200009327000008250400000","1111100057032230000","0163103100000010"); // 0-5
    cuts.AddCutPCMCalo("31210013","00200009327000008250400000","1111100057032230000","0163103100000010"); // 5-10
  } else if (trainConfig == 202){ // EMCAL clusters semi-central
    cuts.AddCutPCMCalo("11210013","00200009327000008250400000","1111100057032230000","0163103100000010"); // 10-20
    cuts.AddCutPCMCalo("12310013","00200009327000008250400000","1111100057032230000","0163103100000010"); // 20-30
    cuts.AddCutPCMCalo("13410013","00200009327000008250400000","1111100057032230000","0163103100000010"); // 30-40
    cuts.AddCutPCMCalo("12410013","00200009327000008250400000","1111100057032230000","0163103100000010"); // 20-40
  } else if (trainConfig == 203){ // EMCAL clusters semi-central
    cuts.AddCutPCMCalo("14510013","00200009327000008250400000","1111100057032230000","0163103100000010"); // 40-50
    cuts.AddCutPCMCalo("14610013","00200009327000008250400000","1111100057032230000","0163103100000010"); // 40-60
    cuts.AddCutPCMCalo("15610013","00200009327000008250400000","1111100057032230000","0163103100000010"); // 40-60
  } else if (trainConfig == 204){ // EMCAL clusters peripheral
    cuts.AddCutPCMCalo("16710013","00200009327000008250400000","1111100057032230000","0163103100000010"); // 60-70
    cuts.AddCutPCMCalo("17810013","00200009327000008250400000","1111100057032230000","0163103100000010"); // 70-80
    cuts.AddCutPCMCalo("18910013","00200009327000008250400000","1111100057032230000","0163103100000010"); // 80-90
    cuts.AddCutPCMCalo("16810013","00200009327000008250400000","1111100057032230000","0163103100000010"); // 60-80

  } else if (trainConfig == 205){ // EMCAL clusters central add sig
    cuts.AddCutPCMCalo("10110023","00200009327000008250400000","1111100057032230000","0163103100000010"); // 0-10
    cuts.AddCutPCMCalo("30110023","00200009327000008250400000","1111100057032230000","0163103100000010"); // 0-5
    cuts.AddCutPCMCalo("31210023","00200009327000008250400000","1111100057032230000","0163103100000010"); // 5-10
  } else if (trainConfig == 206){ // EMCAL clusters semi-central add sig
    cuts.AddCutPCMCalo("11210023","00200009327000008250400000","1111100057032230000","0163103100000010"); // 10-20
    cuts.AddCutPCMCalo("12310023","00200009327000008250400000","1111100057032230000","0163103100000010"); // 20-30
    cuts.AddCutPCMCalo("13410023","00200009327000008250400000","1111100057032230000","0163103100000010"); // 30-40
    cuts.AddCutPCMCalo("12410023","00200009327000008250400000","1111100057032230000","0163103100000010"); // 20-40
  } else if (trainConfig == 207){ // EMCAL clusters semi-central add sig
    cuts.AddCutPCMCalo("14510023","00200009327000008250400000","1111100057032230000","0163103100000010"); // 40-50
    cuts.AddCutPCMCalo("14610023","00200009327000008250400000","1111100057032230000","0163103100000010"); // 40-60
    cuts.AddCutPCMCalo("15610023","00200009327000008250400000","1111100057032230000","0163103100000010"); // 40-60
  } else if (trainConfig == 208){ // EMCAL clusters peripheral add sig
    cuts.AddCutPCMCalo("16710023","00200009327000008250400000","1111100057032230000","0163103100000010"); // 60-70
    cuts.AddCutPCMCalo("17810023","00200009327000008250400000","1111100057032230000","0163103100000010"); // 70-80
    cuts.AddCutPCMCalo("18910023","00200009327000008250400000","1111100057032230000","0163103100000010"); // 80-90
    cuts.AddCutPCMCalo("16810023","00200009327000008250400000","1111100057032230000","0163103100000010"); // 60-80

  } else if (trainConfig == 209){ // EMCAL clusters central add sig
    cuts.AddCutPCMCalo("10110023","00200009327000008250400000","1111100057032230000","0163103100000010"); // 0-10
    cuts.AddCutPCMCalo("30110023","00200009327000008250400000","1111100057032230000","0163103100000010"); // 0-5
    cuts.AddCutPCMCalo("31210023","00200009327000008250400000","1111100057032230000","0163103100000010"); // 5-10
  } else if (trainConfig == 210){ // EMCAL clusters semi-central add sig
    cuts.AddCutPCMCalo("11210023","00200009327000008250400000","1111100057032230000","0163103100000010"); // 10-20
    cuts.AddCutPCMCalo("12310023","00200009327000008250400000","1111100057032230000","0163103100000010"); // 20-30
    cuts.AddCutPCMCalo("13410023","00200009327000008250400000","1111100057032230000","0163103100000010"); // 30-40
    cuts.AddCutPCMCalo("12410023","00200009327000008250400000","1111100057032230000","0163103100000010"); // 20-40
  } else if (trainConfig == 211){ // EMCAL clusters semi-central add sig
    cuts.AddCutPCMCalo("14510023","00200009327000008250400000","1111100057032230000","0163103100000010"); // 40-50
    cuts.AddCutPCMCalo("14610023","00200009327000008250400000","1111100057032230000","0163103100000010"); // 40-60
    cuts.AddCutPCMCalo("15610023","00200009327000008250400000","1111100057032230000","0163103100000010"); // 40-60
  } else if (trainConfig == 212){ // EMCAL clusters peripheral add sig
    cuts.AddCutPCMCalo("16710023","00200009327000008250400000","1111100057032230000","0163103100000010"); // 60-70
    cuts.AddCutPCMCalo("17810023","00200009327000008250400000","1111100057032230000","0163103100000010"); // 70-80
    cuts.AddCutPCMCalo("18910023","00200009327000008250400000","1111100057032230000","0163103100000010"); // 80-90
    cuts.AddCutPCMCalo("16810023","00200009327000008250400000","1111100057032230000","0163103100000010"); // 60-80

  } else if (trainConfig == 213){ // EMCAL clusters central TB
    cuts.AddCutPCMCalo("10110013","00200009327000008250400000","1111102057032230000","0163103100000010"); // 0-10
    cuts.AddCutPCMCalo("30110013","00200009327000008250400000","1111102057032230000","0163103100000010"); // 0-5
    cuts.AddCutPCMCalo("31210013","00200009327000008250400000","1111102057032230000","0163103100000010"); // 5-10
  } else if (trainConfig == 214){ // EMCAL clusters semi-central TB
    cuts.AddCutPCMCalo("11210013","00200009327000008250400000","1111102057032230000","0163103100000010"); // 10-20
    cuts.AddCutPCMCalo("12310013","00200009327000008250400000","1111102057032230000","0163103100000010"); // 20-30
    cuts.AddCutPCMCalo("13410013","00200009327000008250400000","1111102057032230000","0163103100000010"); // 30-40
    cuts.AddCutPCMCalo("12410013","00200009327000008250400000","1111102057032230000","0163103100000010"); // 20-40
  } else if (trainConfig == 215){ // EMCAL clusters semi-central TB
    cuts.AddCutPCMCalo("14510013","00200009327000008250400000","1111102057032230000","0163103100000010"); // 40-50
    cuts.AddCutPCMCalo("14610013","00200009327000008250400000","1111102057032230000","0163103100000010"); // 40-60
    cuts.AddCutPCMCalo("15610013","00200009327000008250400000","1111102057032230000","0163103100000010"); // 40-60
  } else if (trainConfig == 216){ // EMCAL clusters peripheral TB
    cuts.AddCutPCMCalo("16710013","00200009327000008250400000","1111102057032230000","0163103100000010"); // 60-70
    cuts.AddCutPCMCalo("17810013","00200009327000008250400000","1111102057032230000","0163103100000010"); // 70-80
    cuts.AddCutPCMCalo("18910013","00200009327000008250400000","1111102057032230000","0163103100000010"); // 80-90
    cuts.AddCutPCMCalo("16810013","00200009327000008250400000","1111102057032230000","0163103100000010"); // 60-80

  //****************************************************************************************************
  // EMCal 5TeV Pb-Pb LHC15o EMC triggers
  //****************************************************************************************************
  } else if (trainConfig == 221){ // EMCAL clusters central
    cuts.AddCutPCMCalo("10183013","00200009327000008250400000","1111100057032230000","0163103100000010"); // 0-10
    cuts.AddCutPCMCalo("30183013","00200009327000008250400000","1111100057032230000","0163103100000010"); // 0-5
    cuts.AddCutPCMCalo("31283013","00200009327000008250400000","1111100057032230000","0163103100000010"); // 5-10
  } else if (trainConfig == 222){ // EMCAL clusters semi-central
    cuts.AddCutPCMCalo("11283013","00200009327000008250400000","1111100057032230000","0163103100000010"); // 10-20
    cuts.AddCutPCMCalo("12383013","00200009327000008250400000","1111100057032230000","0163103100000010"); // 20-30
    cuts.AddCutPCMCalo("13483013","00200009327000008250400000","1111100057032230000","0163103100000010"); // 30-40
    cuts.AddCutPCMCalo("12483013","00200009327000008250400000","1111100057032230000","0163103100000010"); // 20-40
  } else if (trainConfig == 223){ // EMCAL clusters semi-central
    cuts.AddCutPCMCalo("14583013","00200009327000008250400000","1111100057032230000","0163103100000010"); // 40-50
    cuts.AddCutPCMCalo("14683013","00200009327000008250400000","1111100057032230000","0163103100000010"); // 40-60
    cuts.AddCutPCMCalo("15683013","00200009327000008250400000","1111100057032230000","0163103100000010"); // 40-60
  } else if (trainConfig == 224){ // EMCAL clusters peripheral
    cuts.AddCutPCMCalo("16783013","00200009327000008250400000","1111100057032230000","0163103100000010"); // 60-70
    cuts.AddCutPCMCalo("17883013","00200009327000008250400000","1111100057032230000","0163103100000010"); // 70-80
    cuts.AddCutPCMCalo("18983013","00200009327000008250400000","1111100057032230000","0163103100000010"); // 80-90
    cuts.AddCutPCMCalo("16883013","00200009327000008250400000","1111100057032230000","0163103100000010"); // 60-80

  } else if (trainConfig == 226){ // EMCAL clusters - peripheral centrality selection for PbPb EMCal
    cuts.AddCutPCMCalo("15910013","00200009327000008250400000","11111020530a2230000","0163103100000010"); //
    cuts.AddCutPCMCalo("15910013","00200009327000008250400000","11111870530a2230000","0163103100000010"); //
    cuts.AddCutPCMCalo("15910013","00200009327000008250400000","11111020530b2230000","0163103100000010"); //
    cuts.AddCutPCMCalo("15910013","00200009327000008250400000","11111870530b2230000","0163103100000010"); //

  } else if (trainConfig == 230){ // EMCAL clusters - correction convcalo f1
    cuts.AddCutPCMCalo("10110013","00200009327000008250400000","1111181003032230000","0163103100000010"); // 0-10
    cuts.AddCutPCMCalo("11210013","00200009327000008250400000","1111181003032230000","0163103100000010"); // 10-20
    cuts.AddCutPCMCalo("12510013","00200009327000008250400000","1111181003032230000","0163103100000010"); // 20-50
    cuts.AddCutPCMCalo("15910013","00200009327000008250400000","1111181003032230000","0163103100000010"); // 50-90
  } else if (trainConfig == 231){ // EMCAL clusters - correction calocalo f2
    cuts.AddCutPCMCalo("10110013","00200009327000008250400000","1111192003032230000","0163103100000010"); // 0-10
    cuts.AddCutPCMCalo("11210013","00200009327000008250400000","1111192003032230000","0163103100000010"); // 10-20
    cuts.AddCutPCMCalo("12510013","00200009327000008250400000","1111192003032230000","0163103100000010"); // 20-50
    cuts.AddCutPCMCalo("15910013","00200009327000008250400000","1111192003032230000","0163103100000010"); // 50-90
  } else if (trainConfig == 232){ // EMCAL clusters - 0-90% centrality for PbPb EMCal cluster QA
    cuts.AddCutPCMCalo("10910013","00200009327000008250400000","1111100003032230000","0163103100000010"); // 0-90
  } else if (trainConfig == 233){ // EMCAL clusters - centrality selection for PbPb EMCal
    cuts.AddCutPCMCalo("10110013","00200009327000008250400000","1111184053032230000","0163103100000010"); //  0-10 calo correction cent dep
    cuts.AddCutPCMCalo("11210013","00200009327000008250400000","1111185053032230000","0163103100000010"); // 10-20 calo correction cent dep
    cuts.AddCutPCMCalo("12510013","00200009327000008250400000","1111186053032230000","0163103100000010"); // 20-50 calo correction cent dep
    cuts.AddCutPCMCalo("15910013","00200009327000008250400000","1111187053032230000","0163103100000010"); // 50-90 calo correction cent dep

  } else if (trainConfig == 240){ // EMCAL clusters - 0-90% centrality for PbPb EMCal cluster QA
    cuts.AddCutPCMCalo("50910013","00200009327000008250400000","1111100053032230000","0163103100000010"); //  0-90 calo correction cent dep
    cuts.AddCutPCMCalo("50910013","00200009327000008250400000","1111102053032230000","0163103100000010"); //  0-90 calo correction cent dep
    cuts.AddCutPCMCalo("50910013","00200009327000008250400000","1111183053032230000","0163103100000010"); //  0-90 calo correction cent dep
  } else if (trainConfig == 241){ // EMCAL clusters - 0-90% centrality for PbPb EMCal cluster QA
    cuts.AddCutPCMCalo("50910613","00200009327000008250400000","1111100053032230000","0163103100000010"); //  0-90 calo correction cent dep
    cuts.AddCutPCMCalo("50910613","00200009327000008250400000","1111102053032230000","0163103100000010"); //  0-90 calo correction cent dep
    cuts.AddCutPCMCalo("50910613","00200009327000008250400000","1111183053032230000","0163103100000010"); //  0-90 calo correction cent dep
  } else if (trainConfig == 242){ // EMCAL clusters - centrality selection for PbPb EMCal
    cuts.AddCutPCMCalo("50110013","00200009327000008250400000","1111100053032230000","0163103100000010"); //  0-10 calo correction cent dep
    cuts.AddCutPCMCalo("51210013","00200009327000008250400000","1111100053032230000","0163103100000010"); // 10-20 calo correction cent dep
    cuts.AddCutPCMCalo("52510013","00200009327000008250400000","1111100053032230000","0163103100000010"); // 20-50 calo correction cent dep
    cuts.AddCutPCMCalo("55910013","00200009327000008250400000","1111100053032230000","0163103100000010"); // 50-90 calo correction cent dep
  } else if (trainConfig == 243){ // EMCAL clusters - centrality selection for PbPb EMCal
    cuts.AddCutPCMCalo("50110613","00200009327000008250400000","1111100053032230000","0163103100000010"); //  0-10 calo correction cent dep
    cuts.AddCutPCMCalo("51210613","00200009327000008250400000","1111100053032230000","0163103100000010"); // 10-20 calo correction cent dep
    cuts.AddCutPCMCalo("52510613","00200009327000008250400000","1111100053032230000","0163103100000010"); // 20-50 calo correction cent dep
    cuts.AddCutPCMCalo("55910613","00200009327000008250400000","1111100053032230000","0163103100000010"); // 50-90 calo correction cent dep
  } else if (trainConfig == 244){ // EMCAL clusters - centrality selection for PbPb EMCal
    cuts.AddCutPCMCalo("50110613","00200009327000008250400000","1111102053032230000","0163103100000010"); //  0-10 calo correction cent dep
    cuts.AddCutPCMCalo("51210613","00200009327000008250400000","1111102053032230000","0163103100000010"); // 10-20 calo correction cent dep
    cuts.AddCutPCMCalo("52510613","00200009327000008250400000","1111102053032230000","0163103100000010"); // 20-50 calo correction cent dep
    cuts.AddCutPCMCalo("55910613","00200009327000008250400000","1111102053032230000","0163103100000010"); // 50-90 calo correction cent dep
  } else if (trainConfig == 245){ // EMCAL clusters - centrality selection for PbPb EMCal
    cuts.AddCutPCMCalo("50110613","00200009327000008250400000","1111184053032230000","0163103100000010"); //  0-10 calo correction cent dep
    cuts.AddCutPCMCalo("51210613","00200009327000008250400000","1111185053032230000","0163103100000010"); // 10-20 calo correction cent dep
    cuts.AddCutPCMCalo("52510613","00200009327000008250400000","1111186053032230000","0163103100000010"); // 20-50 calo correction cent dep
    cuts.AddCutPCMCalo("55910613","00200009327000008250400000","1111187053032230000","0163103100000010"); // 50-90 calo correction cent dep
  } else if (trainConfig == 246){ // EMCAL clusters - 0-90% centrality for PbPb EMCal cluster QA
    cuts.AddCutPCMCalo("50910613","00200009327000008250400000","1111183050032230000","0163103100000010"); //  0-90 calo correction cent dep
    cuts.AddCutPCMCalo("50910613","00200009327000008250400000","1111183051032230000","0163103100000010"); //  0-90 calo correction cent dep
    cuts.AddCutPCMCalo("50910613","00200009327000008250400000","1111183057032230000","0163103100000010"); //  0-90 calo correction cent dep
  } else if (trainConfig == 247){ // EMCAL clusters - centrality selection for PbPb EMCal
    cuts.AddCutPCMCalo("50110613","00200009327000008250400000","1111184050032230000","0163103100000010"); //  0-10 calo correction cent dep
    cuts.AddCutPCMCalo("51210613","00200009327000008250400000","1111185050032230000","0163103100000010"); // 10-20 calo correction cent dep
    cuts.AddCutPCMCalo("52510613","00200009327000008250400000","1111186050032230000","0163103100000010"); // 20-50 calo correction cent dep
    cuts.AddCutPCMCalo("55910613","00200009327000008250400000","1111187050032230000","0163103100000010"); // 50-90 calo correction cent dep
  } else if (trainConfig == 248){ // EMCAL clusters - centrality selection for PbPb EMCal
    cuts.AddCutPCMCalo("50110613","00200009327000008250400000","1111184051032230000","0163103100000010"); //  0-10 calo correction cent dep
    cuts.AddCutPCMCalo("51210613","00200009327000008250400000","1111185051032230000","0163103100000010"); // 10-20 calo correction cent dep
    cuts.AddCutPCMCalo("52510613","00200009327000008250400000","1111186051032230000","0163103100000010"); // 20-50 calo correction cent dep
    cuts.AddCutPCMCalo("55910613","00200009327000008250400000","1111187051032230000","0163103100000010"); // 50-90 calo correction cent dep
  } else if (trainConfig == 249){ // EMCAL clusters - centrality selection for PbPb EMCal
    cuts.AddCutPCMCalo("50110613","00200009327000008250400000","1111184057032230000","0163103100000010"); //  0-10 calo correction cent dep
    cuts.AddCutPCMCalo("51210613","00200009327000008250400000","1111185057032230000","0163103100000010"); // 10-20 calo correction cent dep
    cuts.AddCutPCMCalo("52510613","00200009327000008250400000","1111186057032230000","0163103100000010"); // 20-50 calo correction cent dep
    cuts.AddCutPCMCalo("55910613","00200009327000008250400000","1111187057032230000","0163103100000010"); // 50-90 calo correction cent dep
  } else if (trainConfig == 250){ // EMCAL clusters - 0-90% centrality for PbPb EMCal cluster QA
    cuts.AddCutPCMCalo("10910013","00200009327000008250400000","1111183051032230000","0163103100000010"); //  0-90 calo correction cent dep
    cuts.AddCutPCMCalo("10910a13","00200009327000008250400000","1111183051032230000","0163103100000010"); //  0-90 calo correction cent dep
  } else if (trainConfig == 251){ // EMCAL clusters - 0-90% centrality for PbPb EMCal cluster QA
    cuts.AddCutPCMCalo("10110013","00200009327000008250400000","1111183051032230000","0163103100000010"); //  0-90 calo correction cent dep
    cuts.AddCutPCMCalo("11210013","00200009327000008250400000","1111183051032230000","0163103100000010"); //  0-90 calo correction cent dep
    cuts.AddCutPCMCalo("12510013","00200009327000008250400000","1111183051032230000","0163103100000010"); //  0-90 calo correction cent dep
    cuts.AddCutPCMCalo("15910013","00200009327000008250400000","1111183051032230000","0163103100000010"); //  0-90 calo correction cent dep
  } else if (trainConfig == 252){ // EMCAL clusters - 0-90% centrality for PbPb EMCal cluster QA
    cuts.AddCutPCMCalo("10110a13","00200009327000008250400000","1111183051032230000","0163103100000010"); //  0-90 calo correction cent dep
    cuts.AddCutPCMCalo("11210a13","00200009327000008250400000","1111183051032230000","0163103100000010"); //  0-90 calo correction cent dep
    cuts.AddCutPCMCalo("12510a13","00200009327000008250400000","1111183051032230000","0163103100000010"); //  0-90 calo correction cent dep
    cuts.AddCutPCMCalo("15910a13","00200009327000008250400000","1111183051032230000","0163103100000010"); //  0-90 calo correction cent dep
  } else if (trainConfig == 253){ // EMCAL clusters - 20180718 - default without corrections
    cuts.AddCutPCMCalo("10110a13","00200009327000008250400000","1111100051032230000","0163103100000010"); //  0-90 calo correction cent dep
    cuts.AddCutPCMCalo("11210a13","00200009327000008250400000","1111100051032230000","0163103100000010"); //  0-90 calo correction cent dep
    cuts.AddCutPCMCalo("12410a13","00200009327000008250400000","1111100051032230000","0163103100000010"); //  0-90 calo correction cent dep
    cuts.AddCutPCMCalo("14610a13","00200009327000008250400000","1111100051032230000","0163103100000010"); //  0-90 calo correction cent dep
    cuts.AddCutPCMCalo("16810a13","00200009327000008250400000","1111100051032230000","0163103100000010"); //  0-90 calo correction cent dep
  } else if (trainConfig == 254){ // EMCAL clusters - 20180718 - default with peripheral corrections
    cuts.AddCutPCMCalo("10110a13","00200009327000008250400000","1111187051032230000","0163103100000010"); //
    cuts.AddCutPCMCalo("11210a13","00200009327000008250400000","1111187051032230000","0163103100000010"); //
    cuts.AddCutPCMCalo("12410a13","00200009327000008250400000","1111187051032230000","0163103100000010"); //
    cuts.AddCutPCMCalo("14610a13","00200009327000008250400000","1111187051032230000","0163103100000010"); //
    cuts.AddCutPCMCalo("16810a13","00200009327000008250400000","1111187051032230000","0163103100000010"); //
  } else if (trainConfig == 255){ // EMCAL clusters - 20180718 - default with peripheral corrections
    cuts.AddCutPCMCalo("12410a13","00200009327000008250400000","1111186051032230000","0163103100000010"); //
    cuts.AddCutPCMCalo("12510a13","00200009327000008250400000","1111186051032230000","0163103100000010"); //
    cuts.AddCutPCMCalo("14610a13","00200009327000008250400000","1111187051032230000","0163103100000010"); //
    cuts.AddCutPCMCalo("16810a13","00200009327000008250400000","1111187051032230000","0163103100000010"); //
  } else if (trainConfig == 256){ // EMCAL clusters - 20180718 - default with corrections
    cuts.AddCutPCMCalo("30110a13","00200009327000008250400000","1111184051032230000","0163103100000010"); //
    cuts.AddCutPCMCalo("31210a13","00200009327000008250400000","1111184051032230000","0163103100000010"); //
    cuts.AddCutPCMCalo("10110a13","00200009327000008250400000","1111184051032230000","0163103100000010"); //
    cuts.AddCutPCMCalo("11210a13","00200009327000008250400000","1111185051032230000","0163103100000010"); //
  } else if (trainConfig == 257){ // EMCAL clusters - 20180718 - default with peripheral corrections
    cuts.AddCutPCMCalo("14610a13","00200009327000008250400000","1111187050032230000","0163103100000010"); //
    cuts.AddCutPCMCalo("14610a13","00200009327000008250400000","1111187051032230000","0163103100000010"); //
    cuts.AddCutPCMCalo("14610a13","00200009327000008250400000","111118705i032230000","0163103100000010"); //
    cuts.AddCutPCMCalo("14610a13","00200009327000008250400000","111118705j032230000","0163103100000010"); //
    cuts.AddCutPCMCalo("14610a13","00200009327000008250400000","111118705k032230000","0163103100000010"); //
  } else if (trainConfig == 258){ // EMCAL clusters - 20180718 - default with corrections
    cuts.AddCutPCMCalo("10110a13","00200009327000008250400000","1111184050032230000","0163103100000010"); //
    cuts.AddCutPCMCalo("10110a13","00200009327000008250400000","1111184051032230000","0163103100000010"); //
    cuts.AddCutPCMCalo("10110a13","00200009327000008250400000","111118405i032230000","0163103100000010"); //
    cuts.AddCutPCMCalo("10110a13","00200009327000008250400000","111118405j032230000","0163103100000010"); //
    cuts.AddCutPCMCalo("10110a13","00200009327000008250400000","111118405k032230000","0163103100000010"); //
  } else if (trainConfig == 260){ // EMCAL clusters - 20181110 - No NL
    cuts.AddCutPCMCalo("12410a13","00200009327000008250400000","111110005k032230000","0h63103100000010"); //
    cuts.AddCutPCMCalo("12510a13","00200009327000008250400000","111110005k032230000","0h63103100000010"); //
    cuts.AddCutPCMCalo("14610a13","00200009327000008250400000","111110005k032230000","0h63103100000010"); //
    cuts.AddCutPCMCalo("16810a13","00200009327000008250400000","111110005k032230000","0h63103100000010"); //
  } else if (trainConfig == 261){ // EMCAL clusters - 20181110 - No NL
    cuts.AddCutPCMCalo("30110a13","00200009327000008250400000","111110005k0a2230000","0h63103100000010"); //
    cuts.AddCutPCMCalo("31210a13","00200009327000008250400000","111110005k0a2230000","0h63103100000010"); //
    cuts.AddCutPCMCalo("10110a13","00200009327000008250400000","111110005k0a2230000","0h63103100000010"); //
    cuts.AddCutPCMCalo("11210a13","00200009327000008250400000","111110005k0b2230000","0h63103100000010"); //
  } else if (trainConfig == 262){ // EMCAL clusters - 20181207 - Peri NL
    cuts.AddCutPCMCalo("12410a13","00200009327000008250400000","111118105k032230000","0h63103100000010"); //
    cuts.AddCutPCMCalo("12510a13","00200009327000008250400000","111118105k032230000","0h63103100000010"); //
    cuts.AddCutPCMCalo("14610a13","00200009327000008250400000","111118105k032230000","0h63103100000010"); //
    cuts.AddCutPCMCalo("16810a13","00200009327000008250400000","111118105k032230000","0h63103100000010"); //
  } else if (trainConfig == 263){ // EMCAL clusters - 20181207 - Peri NL
    cuts.AddCutPCMCalo("30110a13","00200009327000008250400000","111118105k0a2230000","0h63103100000010"); //
    cuts.AddCutPCMCalo("31210a13","00200009327000008250400000","111118105k0a2230000","0h63103100000010"); //
    cuts.AddCutPCMCalo("10110a13","00200009327000008250400000","111118105k0a2230000","0h63103100000010"); //
    cuts.AddCutPCMCalo("11210a13","00200009327000008250400000","111118105k0b2230000","0h63103100000010"); //
  } else if (trainConfig == 265){ // EMCAL clusters - 20181207 - Peri NL
    cuts.AddCutPCMCalo("12410a13","00200009327000008250400000","111118305k032220000","0h63103100000010"); //
    cuts.AddCutPCMCalo("12510a13","00200009327000008250400000","111118305k032220000","0h63103100000010"); //
    cuts.AddCutPCMCalo("14610a13","00200009327000008250400000","111118305k032220000","0h63103100000010"); //
    cuts.AddCutPCMCalo("16810a13","00200009327000008250400000","111118305k032220000","0h63103100000010"); //
  } else if (trainConfig == 266){ // EMCAL clusters - 20181207 - Peri NL
    cuts.AddCutPCMCalo("30110a13","00200009327000008250400000","111118305k0a2220000","0h63103100000010"); //
    cuts.AddCutPCMCalo("31210a13","00200009327000008250400000","111118305k0a2220000","0h63103100000010"); //
    cuts.AddCutPCMCalo("10110a13","00200009327000008250400000","111118305k0a2220000","0h63103100000010"); //
    cuts.AddCutPCMCalo("11210a13","00200009327000008250400000","111118305k0b2220000","0h63103100000010"); //
  } else if (trainConfig == 267){ // EMCAL + DCal clusters - 20190301
    cuts.AddCutPCMCalo("12410a13","00200009327000008250400000","411798305k032220000","0h63103100000010"); //
    cuts.AddCutPCMCalo("12510a13","00200009327000008250400000","411798305k032220000","0h63103100000010"); //
    cuts.AddCutPCMCalo("14610a13","00200009327000008250400000","411798305k032220000","0h63103100000010"); //
    cuts.AddCutPCMCalo("14810a13","00200009327000008250400000","411798305k032220000","0h63103100000010"); //
    cuts.AddCutPCMCalo("16810a13","00200009327000008250400000","411798305k032220000","0h63103100000010"); //
  } else if (trainConfig == 268){ // EMCAL + DCal clusters - 20190301
    cuts.AddCutPCMCalo("30110a13","00200009327000008250400000","411798305k0a2220000","0h63103100000010"); //
    cuts.AddCutPCMCalo("31210a13","00200009327000008250400000","411798305k0a2220000","0h63103100000010"); //
    cuts.AddCutPCMCalo("10110a13","00200009327000008250400000","411798305k0a2220000","0h63103100000010"); //
    cuts.AddCutPCMCalo("11210a13","00200009327000008250400000","411798305k0b2220000","0h63103100000010"); //
    cuts.AddCutPCMCalo("10210a13","00200009327000008250400000","411798305k0a2220000","0h63103100000010"); //


  //systematics for LHC15o 20-40%
  } else if (trainConfig == 280){ // PCM-EMCAL
    cuts.AddCutPCMCalo("12410a13","00200009327000008250400000","111118105k032230000","0h63103100000010"); //
  } else if (trainConfig == 281){ // PCM-EMCAL syst 1/7
    cuts.AddCutPCMCalo("12410a13","00200009327000008250400000","111118105k042230000","0h63103100000010"); // min energy cluster variation 2  800 MeV
    cuts.AddCutPCMCalo("12410a13","00200009327000008250400000","111118105k032230000","0h63106100000010"); // alpha meson variation 1 0<alpha<0.8
    cuts.AddCutPCMCalo("12410a13","00200009327000008250400000","111118105k032230000","0h63103100000010"); // min/max M02  0.1<M<0.5
    cuts.AddCutPCMCalo("12410a13","00200009327000008250400000","111118105k032250000","0h63103100000010"); // min/max M02  0.1<M<0.3
  } else if (trainConfig == 282){ // PCM-EMCAL syst 3/7
    cuts.AddCutPCMCalo("12410a13","00200009327000008250400000","111118105k032230000","0h63303100000010"); // rapidity variation  y<0.6
    cuts.AddCutPCMCalo("12410a13","00200009327000008250400000","111118105k032230000","0h63403100000010"); // rapidity variation  y<0.5
    cuts.AddCutPCMCalo("12410a13","00200009327000008250400000","111118105k032230000","0h63103100000060"); // min opening angle 0.017
    cuts.AddCutPCMCalo("12410a13","00200009327000008250400000","111118105k032230000","0h63103100000070"); // min opening angle 0.016
  } else if (trainConfig == 283){ // PCM-EMCAL syst 5/7
    cuts.AddCutPCMCalo("12410a13","00200009327000008250400000","1111181053032230000","0h63103100000010"); // fixed window
    cuts.AddCutPCMCalo("12410a13","00200009327000008250400000","1111181056032230000","0h63103100000010"); // tm pt dependent var 1
    cuts.AddCutPCMCalo("12410a13","00200009327000008250400000","111110105k032230000","0h63103100000010"); // NL variation, TB
    cuts.AddCutPCMCalo("12410a13","00200009327000008250400000","111118205k032230000","0h63103100000010"); // NL variation, Calo
  } else if (trainConfig == 284){ // PCM-EMCAL variations 1/6
    cuts.AddCutPCMCalo("12410a13","00200009227000008250400000","111118105k032230000","0h63103100000010"); // dEdx e -3, 5
    cuts.AddCutPCMCalo("12410a13","00200009127000008250400000","111118105k032230000","0h63103100000010"); // dEdx e -5, 5
    cuts.AddCutPCMCalo("12410a13","00200009357000008250400000","111118105k032230000","0h63103100000010"); // dEdx pi 2
    cuts.AddCutPCMCalo("12410a13","00200009317000008250400000","111118105k032230000","0h63103100000010"); // dEdx pi 0
    cuts.AddCutPCMCalo("12410a13","00200009387300008250400000","111118105k032230000","0h63103100000010"); // dEdx pi 2 high 1 (> 3.5 GeV)
  } else if (trainConfig == 285){ // PCM-EMCAL variations 2/6
    cuts.AddCutPCMCalo("12410a13","00200009327000009250400000","111118105k032230000","0h63103100000010"); // qt 2D 0.03
    cuts.AddCutPCMCalo("12410a13","00200009327000003250400000","111118105k032230000","0h63103100000010"); // qt 1D 0.05
    cuts.AddCutPCMCalo("12410a13","00200009327000002250400000","111118105k032230000","0h63103100000010"); // qt 1D 0.07
    cuts.AddCutPCMCalo("12410a13","00200049327000008250400000","111118105k032230000","0h63103100000010"); // single pt > 0.075
    cuts.AddCutPCMCalo("12410a13","00200019327000008250400000","111118105k032230000","0h63103100000010"); // single pt > 0.1
  } else if (trainConfig == 286){ // PCM-EMCAL variations 3/6
    cuts.AddCutPCMCalo("12410a13","00200009327000008850400000","111118105k032230000","0h63103100000010"); // 2D psi pair chi2 var
    cuts.AddCutPCMCalo("12410a13","00200009327000008260400000","111118105k032230000","0h63103100000010"); // 2D psi pair chi2 var
    cuts.AddCutPCMCalo("12410a13","00200009327000008860400000","111118105k032230000","0h63103100000010"); // 2D psi pair chi2 var
    cuts.AddCutPCMCalo("12410a13","00200009327000008280400000","111118105k032230000","0h63103100000010"); // 2D psi pair chi2 var
    cuts.AddCutPCMCalo("12410a13","00200009327000008880400000","111118105k032230000","0h63103100000010"); // 2D psi pair chi2 var
  } else if (trainConfig == 287){ // PCM-EMCAL variations pi dEdx 5/6
    cuts.AddCutPCMCalo("12410a13","00200006327000008250400000","111118105k032230000","0h63103100000010"); // min TPC cl > 0.7
    cuts.AddCutPCMCalo("12410a13","00200008327000008250400000","111118105k032230000","0h63103100000010"); // min TPC cl > 0.35
    cuts.AddCutPCMCalo("12410a13","00200009317300008250400000","111118105k032230000","0h63103100000010"); // dEdx pi: 0: 0.4-3.5, -10: 3.5 ->
    cuts.AddCutPCMCalo("12410a13","00200009327300008250400000","111118105k032230000","0h63103100000010"); // dEdx pi: 1: 0.4-3.5, -10: 3.5 ->
    cuts.AddCutPCMCalo("12410a13","00200009325000008250400000","111118105k032230000","0h63103100000010"); // dEdx pi: 1: 0.3 ->
  } else if (trainConfig == 288){ // PCM-EMCAL variations pi dEdx 6/6
    cuts.AddCutPCMCalo("12410a13","00200009327600008250400000","111118105k032230000","0h63103100000010"); // dEdx pi: 1: 0.4-2, -10: 2. ->
    cuts.AddCutPCMCalo("12410a13","00200009327400008250400000","111118105k032230000","0h63103100000010"); // dEdx pi: 1: 0.4-3, -10: 3. ->
    cuts.AddCutPCMCalo("12410a13","00200009315600008250400000","111118105k032230000","0h63103100000010"); // dEdx pi: 0: 0.3-2, -10: 2. ->
    cuts.AddCutPCMCalo("12410a13","00200009367400008250400000","111118105k032230000","0h63103100000010"); // dEdx pi: 2: 0.4-3, 0.5: 3. ->

  //systematics for LHC15o 60-80%
  } else if (trainConfig == 290){ // PCM-EMCAL
    cuts.AddCutPCMCalo("16810a13","00200009327000008250400000","111118105k032230000","0h63103100000010"); //
  } else if (trainConfig == 291){ // PCM-EMCAL syst 1/7
    cuts.AddCutPCMCalo("16810a13","00200009327000008250400000","111118105k042230000","0h63103100000010"); // min energy cluster variation 2  800 MeV
    cuts.AddCutPCMCalo("16810a13","00200009327000008250400000","111118105k032230000","0h63106100000010"); // alpha meson variation 1 0<alpha<0.8
    cuts.AddCutPCMCalo("16810a13","00200009327000008250400000","111118105k032230000","0h63103100000010"); // min/max M02  0.1<M<0.5
    cuts.AddCutPCMCalo("16810a13","00200009327000008250400000","111118105k032250000","0h63103100000010"); // min/max M02  0.1<M<0.3
  } else if (trainConfig == 292){ // PCM-EMCAL syst 3/7
    cuts.AddCutPCMCalo("16810a13","00200009327000008250400000","111118105k032230000","0h63303100000010"); // rapidity variation  y<0.6
    cuts.AddCutPCMCalo("16810a13","00200009327000008250400000","111118105k032230000","0h63403100000010"); // rapidity variation  y<0.5
    cuts.AddCutPCMCalo("16810a13","00200009327000008250400000","111118105k032230000","0h63103100000060"); // min opening angle 0.017
    cuts.AddCutPCMCalo("16810a13","00200009327000008250400000","111118105k032230000","0h63103100000070"); // min opening angle 0.016
  } else if (trainConfig == 293){ // PCM-EMCAL syst 5/7
    cuts.AddCutPCMCalo("16810a13","00200009327000008250400000","1111181053032230000","0h63103100000010"); // fixed window
    cuts.AddCutPCMCalo("16810a13","00200009327000008250400000","1111181056032230000","0h63103100000010"); // tm pt dependent var 1
    cuts.AddCutPCMCalo("16810a13","00200009327000008250400000","111110105k032230000","0h63103100000010"); // NL variation, TB
    cuts.AddCutPCMCalo("16810a13","00200009327000008250400000","111118205k032230000","0h63103100000010"); // NL variation, Calo
  } else if (trainConfig == 294){ // PCM-EMCAL variations 1/6
    cuts.AddCutPCMCalo("16810a13","00200009227000008250400000","111118105k032230000","0h63103100000010"); // dEdx e -3, 5
    cuts.AddCutPCMCalo("16810a13","00200009127000008250400000","111118105k032230000","0h63103100000010"); // dEdx e -5, 5
    cuts.AddCutPCMCalo("16810a13","00200009357000008250400000","111118105k032230000","0h63103100000010"); // dEdx pi 2
    cuts.AddCutPCMCalo("16810a13","00200009317000008250400000","111118105k032230000","0h63103100000010"); // dEdx pi 0
    cuts.AddCutPCMCalo("16810a13","00200009387300008250400000","111118105k032230000","0h63103100000010"); // dEdx pi 2 high 1 (> 3.5 GeV)
  } else if (trainConfig == 295){ // PCM-EMCAL variations 2/6
    cuts.AddCutPCMCalo("16810a13","00200009327000009250400000","111118105k032230000","0h63103100000010"); // qt 2D 0.03
    cuts.AddCutPCMCalo("16810a13","00200009327000003250400000","111118105k032230000","0h63103100000010"); // qt 1D 0.05
    cuts.AddCutPCMCalo("16810a13","00200009327000002250400000","111118105k032230000","0h63103100000010"); // qt 1D 0.07
    cuts.AddCutPCMCalo("16810a13","00200049327000008250400000","111118105k032230000","0h63103100000010"); // single pt > 0.075
    cuts.AddCutPCMCalo("16810a13","00200019327000008250400000","111118105k032230000","0h63103100000010"); // single pt > 0.1
  } else if (trainConfig == 296){ // PCM-EMCAL variations 3/6
    cuts.AddCutPCMCalo("16810a13","00200009327000008850400000","111118105k032230000","0h63103100000010"); // 2D psi pair chi2 var
    cuts.AddCutPCMCalo("16810a13","00200009327000008260400000","111118105k032230000","0h63103100000010"); // 2D psi pair chi2 var
    cuts.AddCutPCMCalo("16810a13","00200009327000008860400000","111118105k032230000","0h63103100000010"); // 2D psi pair chi2 var
    cuts.AddCutPCMCalo("16810a13","00200009327000008280400000","111118105k032230000","0h63103100000010"); // 2D psi pair chi2 var
    cuts.AddCutPCMCalo("16810a13","00200009327000008880400000","111118105k032230000","0h63103100000010"); // 2D psi pair chi2 var
  } else if (trainConfig == 297){ // PCM-EMCAL variations pi dEdx 5/6
    cuts.AddCutPCMCalo("16810a13","00200006327000008250400000","111118105k032230000","0h63103100000010"); // min TPC cl > 0.7
    cuts.AddCutPCMCalo("16810a13","00200008327000008250400000","111118105k032230000","0h63103100000010"); // min TPC cl > 0.35
    cuts.AddCutPCMCalo("16810a13","00200009317300008250400000","111118105k032230000","0h63103100000010"); // dEdx pi: 0: 0.4-3.5, -10: 3.5 ->
    cuts.AddCutPCMCalo("16810a13","00200009327300008250400000","111118105k032230000","0h63103100000010"); // dEdx pi: 1: 0.4-3.5, -10: 3.5 ->
    cuts.AddCutPCMCalo("16810a13","00200009325000008250400000","111118105k032230000","0h63103100000010"); // dEdx pi: 1: 0.3 ->
  } else if (trainConfig == 298){ // PCM-EMCAL variations pi dEdx 6/6
    cuts.AddCutPCMCalo("16810a13","00200009327600008250400000","111118105k032230000","0h63103100000010"); // dEdx pi: 1: 0.4-2, -10: 2. ->
    cuts.AddCutPCMCalo("16810a13","00200009327400008250400000","111118105k032230000","0h63103100000010"); // dEdx pi: 1: 0.4-3, -10: 3. ->
    cuts.AddCutPCMCalo("16810a13","00200009315600008250400000","111118105k032230000","0h63103100000010"); // dEdx pi: 0: 0.3-2, -10: 2. ->
    cuts.AddCutPCMCalo("16810a13","00200009367400008250400000","111118105k032230000","0h63103100000010"); // dEdx pi: 2: 0.4-3, 0.5: 3. ->

  //****************************************************************************************************
  // EMCal 5TeV Xe-Xe LHC17n
  //****************************************************************************************************
  } else if (trainConfig == 300){ // EMCAL clusters - 0-80% centrality - optimized PCM cuts
    cuts.AddCutPCMCalo("10810013","00200089f9730000iih0400000","4117901017032230000","0163103100000010"); // 0-80
  } else if (trainConfig == 301){ // EMCAL clusters  - optimized PCM cuts
    cuts.AddCutPCMCalo("10210013","00200089f9730000iih0400000","4117901017032230000","0163103100000010"); // 0-20
    cuts.AddCutPCMCalo("12410013","00200089f9730000iih0400000","4117901017032230000","0163103100000010"); // 20-40
    cuts.AddCutPCMCalo("10410013","00200089f9730000iih0400000","4117901017032230000","0163103100000010"); // 0-40
    cuts.AddCutPCMCalo("14810013","00200089f9730000iih0400000","4117901017032230000","0163103100000010"); // 40-80
  } else if (trainConfig == 302){ // EMCAL clusters -  EMCal cluster QA TB calib - optimized PCM cuts
    cuts.AddCutPCMCalo("10810013","00200089f9730000iih0400000","4117931017032230000","0163103100000010"); // 0-80
    cuts.AddCutPCMCalo("10210013","00200089f9730000iih0400000","4117931017032230000","0163103100000010"); // 0-20
    cuts.AddCutPCMCalo("12410013","00200089f9730000iih0400000","4117931017032230000","0163103100000010"); // 20-40
    cuts.AddCutPCMCalo("10410013","00200089f9730000iih0400000","4117931017032230000","0163103100000010"); // 0-40
    cuts.AddCutPCMCalo("14810013","00200089f9730000iih0400000","4117931017032230000","0163103100000010"); // 40-80

  //****************************************************************************************************
  // PHOS 5TeV Xe-Xe LHC17n
  //****************************************************************************************************
  } else if (trainConfig == 400){ // PHOS clusters - 0-90% centrality
    cuts.AddCutPCMCalo("10810013","00200089f9730000iih0400000","2446600000012200000","0163103100000010"); // 0-80%
    cuts.AddCutPCMCalo("10810013","00200089f9730000iih0400000","2446600004012200000","0163103100000010"); // 0-80%
  } else if (trainConfig == 401) {  // PHOS Cent dep Xe-Xe
    cuts.AddCutPCMCalo("10210013","00200089f9730000iih0400000","2446600011012200000","0163103100000010"); // 0-20%
    cuts.AddCutPCMCalo("12410013","00200089f9730000iih0400000","2446600011012200000","0163103100000010"); // 20-40%
    cuts.AddCutPCMCalo("10410013","00200089f9730000iih0400000","2446600011012200000","0163103100000010"); // 0-40%
    cuts.AddCutPCMCalo("14810013","00200089f9730000iih0400000","2446600011012200000","0163103100000010"); // 40-80%
  } else if (trainConfig == 402) {  // PHOS Cent dep Xe-Xe - optimized PCM cuts
    cuts.AddCutPCMCalo("10210013","00200089f9730000iih0400000","2446600000012200000","0163103100000010"); // 0-20%
    cuts.AddCutPCMCalo("12410013","00200089f9730000iih0400000","2446600000012200000","0163103100000010"); // 20-40%
    cuts.AddCutPCMCalo("10410013","00200089f9730000iih0400000","2446600000012200000","0163103100000010"); // 0-40%
    cuts.AddCutPCMCalo("14810013","00200089f9730000iih0400000","2446600000012200000","0163103100000010"); // 40-80%
  } else if (trainConfig == 403) {  // PHOS Cent dep Xe-Xe - optimized PCM cuts
    cuts.AddCutPCMCalo("10210013","00200089f9730000iih0400000","2446600007012200000","0163103100000010"); // 0-20%
    cuts.AddCutPCMCalo("12410013","00200089f9730000iih0400000","2446600007012200000","0163103100000010"); // 20-40%
    cuts.AddCutPCMCalo("10410013","00200089f9730000iih0400000","2446600007012200000","0163103100000010"); // 0-40%
    cuts.AddCutPCMCalo("14810013","00200089f9730000iih0400000","2446600007012200000","0163103100000010"); // 40-80%
  } else if (trainConfig == 404) {  // PHOS Cent dep Xe-Xe - optimized PCM cuts
    cuts.AddCutPCMCalo("10210013","00200089f9730000iih0400000","2446600008012200000","0163103100000010"); // 0-20%
    cuts.AddCutPCMCalo("12410013","00200089f9730000iih0400000","2446600008012200000","0163103100000010"); // 20-40%
    cuts.AddCutPCMCalo("10410013","00200089f9730000iih0400000","2446600008012200000","0163103100000010"); // 0-40%
    cuts.AddCutPCMCalo("14810013","00200089f9730000iih0400000","2446600008012200000","0163103100000010"); // 40-80%
  } else if (trainConfig == 405) {  // PHOS Cent dep Xe-Xe - optimized PCM cuts
    cuts.AddCutPCMCalo("10210013","00200089f9730000iih0400000","2446600004012200000","0163103100000010"); // 0-20%
    cuts.AddCutPCMCalo("12410013","00200089f9730000iih0400000","2446600004012200000","0163103100000010"); // 20-40%
    cuts.AddCutPCMCalo("10410013","00200089f9730000iih0400000","2446600004012200000","0163103100000010"); // 0-40%
    cuts.AddCutPCMCalo("14810013","00200089f9730000iih0400000","2446600004012200000","0163103100000010"); // 40-80%

  //****************************************************************************************************
  // PHOS 5TeV Pb-Pb LHC15o
  //****************************************************************************************************
  } else if (trainConfig == 601){ // EMCAL clusters central
    cuts.AddCutPCMCalo("10110013","00200009327000008250400000","2446600051012200000","0163103100000010"); // 0-10
    cuts.AddCutPCMCalo("30110013","00200009327000008250400000","2446600051012200000","0163103100000010"); // 0-5
    cuts.AddCutPCMCalo("31210013","00200009327000008250400000","2446600051012200000","0163103100000010"); // 5-10
  } else if (trainConfig == 602){ // EMCAL clusters semi-central
    cuts.AddCutPCMCalo("11210013","00200009327000008250400000","2446600051012200000","0163103100000010"); // 10-20
    cuts.AddCutPCMCalo("12310013","00200009327000008250400000","2446600051012200000","0163103100000010"); // 20-30
    cuts.AddCutPCMCalo("13410013","00200009327000008250400000","2446600051012200000","0163103100000010"); // 30-40
    cuts.AddCutPCMCalo("12410013","00200009327000008250400000","2446600051012200000","0163103100000010"); // 20-40
  } else if (trainConfig == 603){ // EMCAL clusters semi-central
    cuts.AddCutPCMCalo("14510013","00200009327000008250400000","2446600051012200000","0163103100000010"); // 40-50
    cuts.AddCutPCMCalo("14610013","00200009327000008250400000","2446600051012200000","0163103100000010"); // 40-60
    cuts.AddCutPCMCalo("15610013","00200009327000008250400000","2446600051012200000","0163103100000010"); // 40-60
  } else if (trainConfig == 604){ // EMCAL clusters peripheral
    cuts.AddCutPCMCalo("16710013","00200009327000008250400000","2446600051012200000","0163103100000010"); // 60-70
    cuts.AddCutPCMCalo("17810013","00200009327000008250400000","2446600051012200000","0163103100000010"); // 70-80
    cuts.AddCutPCMCalo("18910013","00200009327000008250400000","2446600051012200000","0163103100000010"); // 80-90
    cuts.AddCutPCMCalo("16810013","00200009327000008250400000","2446600051012200000","0163103100000010"); // 60-80
  } else if (trainConfig == 605){ // EMCAL clusters central
    cuts.AddCutPCMCalo("10110013","00200009327000008250400000","2446601051012200000","0163103100000010"); // 0-10
    cuts.AddCutPCMCalo("30110013","00200009327000008250400000","2446601051012200000","0163103100000010"); // 0-5
    cuts.AddCutPCMCalo("31210013","00200009327000008250400000","2446601051012200000","0163103100000010"); // 5-10
  } else if (trainConfig == 606){ // EMCAL clusters semi-central
    cuts.AddCutPCMCalo("11210013","00200009327000008250400000","2446601051012200000","0163103100000010"); // 10-20
    cuts.AddCutPCMCalo("12310013","00200009327000008250400000","2446601051012200000","0163103100000010"); // 20-30
    cuts.AddCutPCMCalo("13410013","00200009327000008250400000","2446601051012200000","0163103100000010"); // 30-40
    cuts.AddCutPCMCalo("12410013","00200009327000008250400000","2446601051012200000","0163103100000010"); // 20-40
  } else if (trainConfig == 607){ // EMCAL clusters semi-central
    cuts.AddCutPCMCalo("14510013","00200009327000008250400000","2446601051012200000","0163103100000010"); // 40-50
    cuts.AddCutPCMCalo("14610013","00200009327000008250400000","2446601051012200000","0163103100000010"); // 40-60
    cuts.AddCutPCMCalo("15610013","00200009327000008250400000","2446601051012200000","0163103100000010"); // 40-60
  } else if (trainConfig == 608){ // EMCAL clusters peripheral
    cuts.AddCutPCMCalo("16710013","00200009327000008250400000","2446601051012200000","0163103100000010"); // 60-70
    cuts.AddCutPCMCalo("17810013","00200009327000008250400000","2446601051012200000","0163103100000010"); // 70-80
    cuts.AddCutPCMCalo("18910013","00200009327000008250400000","2446601051012200000","0163103100000010"); // 80-90
    cuts.AddCutPCMCalo("16810013","00200009327000008250400000","2446601051012200000","0163103100000010"); // 60-80

  } else if (trainConfig == 610){ // PHOS clusters - 20181018
    cuts.AddCutPCMCalo("12410a13","00200009327000008250400000","2446600054012200000","0h63103100000010"); //
    cuts.AddCutPCMCalo("12510a13","00200009327000008250400000","2446600054012200000","0h63103100000010"); //
    cuts.AddCutPCMCalo("14610a13","00200009327000008250400000","2446600054012200000","0h63103100000010"); //
    cuts.AddCutPCMCalo("14810a13","00200009327000008250400000","2446600054012200000","0h63103100000010"); //
    cuts.AddCutPCMCalo("16810a13","00200009327000008250400000","2446600054012200000","0h63103100000010"); //
  } else if (trainConfig == 611){ // PHOS clusters - 20181018
    cuts.AddCutPCMCalo("30110a13","00200009327000008250400000","2446600054012200000","0h63103100000010"); //
    cuts.AddCutPCMCalo("31210a13","00200009327000008250400000","2446600054012200000","0h63103100000010"); //
    cuts.AddCutPCMCalo("10110a13","00200009327000008250400000","2446600054012200000","0h63103100000010"); //
    cuts.AddCutPCMCalo("11210a13","00200009327000008250400000","2446600054012200000","0h63103100000010"); //
    cuts.AddCutPCMCalo("10210a13","00200009327000008250400000","2446600054012200000","0h63103100000010"); //
  } else if (trainConfig == 612){ // PHOS clusters - 20190118
    cuts.AddCutPCMCalo("12410a13","00200009327000008250400000","244668105b012200000","0h63103100000010"); //
    cuts.AddCutPCMCalo("12510a13","00200009327000008250400000","244668105b012200000","0h63103100000010"); //
    cuts.AddCutPCMCalo("14610a13","00200009327000008250400000","244668105b012200000","0h63103100000010"); //
    cuts.AddCutPCMCalo("14810a13","00200009327000008250400000","244668105b012200000","0h63103100000010"); //
    cuts.AddCutPCMCalo("16810a13","00200009327000008250400000","244668105b012200000","0h63103100000010"); //
  } else if (trainConfig == 613){ // PHOS clusters - 20190118
    cuts.AddCutPCMCalo("30110a13","00200009327000008250400000","244668105b012200000","0h63103100000010"); //
    cuts.AddCutPCMCalo("31210a13","00200009327000008250400000","244668105b012200000","0h63103100000010"); //
    cuts.AddCutPCMCalo("10110a13","00200009327000008250400000","244668105b012200000","0h63103100000010"); //
    cuts.AddCutPCMCalo("11210a13","00200009327000008250400000","244668105b012200000","0h63103100000010"); //
    cuts.AddCutPCMCalo("10210a13","00200009327000008250400000","244668105b012200000","0h63103100000010"); //

  } else if (trainConfig == 620){ // PHOS clusters - 20190301 - HBT cuts
    cuts.AddCutPCMCalo("12410a13","00200009327000008250400000","244668105b012200000","0h63103100000010"); // std meson cut
    cuts.AddCutPCMCalo("12410a13","00200009327000008250400000","244668105b022200000","0h63103100000010"); // Eclustmin
    cuts.AddCutPCMCalo("12410a13","00600009327000008250400000","244668105b012200000","0h63103100000010"); // R
    cuts.AddCutPCMCalo("12410a13","00200009327000006250400000","244668105b012200000","0h63103100000010"); // Qt
    cuts.AddCutPCMCalo("12410a13","00200009327000008250e00000","244668105b012200000","0h63103100000010"); // cosP
    cuts.AddCutPCMCalo("12410a13","00600009327000006250e00000","244668105b022200000","0h63103100000010"); // R, Qt, cosP, Eclustmin
  } else if (trainConfig == 621){ // PHOS clusters - 20190301 - HBT cuts
    cuts.AddCutPCMCalo("10110a13","00200009327000008250400000","244668105b012200000","0h63103100000010"); // std meson cut
    cuts.AddCutPCMCalo("10110a13","00200009327000008250400000","244668105b022200000","0h63103100000010"); // Eclustmin
    cuts.AddCutPCMCalo("10110a13","00600009327000008250400000","244668105b012200000","0h63103100000010"); // R
    cuts.AddCutPCMCalo("10110a13","00200009327000006250400000","244668105b012200000","0h63103100000010"); // Qt
    cuts.AddCutPCMCalo("10110a13","00200009327000008250e00000","244668105b012200000","0h63103100000010"); // cosP
    cuts.AddCutPCMCalo("10110a13","00600009327000006250e00000","244668105b022200000","0h63103100000010"); // R, Qt, cosP, Eclustmin

  //PCM-PHOS PbPb5TeV 2015 20-40% systematics
  } else if(trainConfig == 630){//std cut
    cuts.AddCutPCMCalo("12410a13","00200009327000008250400000","244668105b012200000","0h63103100000010"); //
  } else if(trainConfig == 631){//dEdx e-line variation and dE/dx pi-line variation
    cuts.AddCutPCMCalo("12410a13","00200009127000008250400000","244668105b012200000","0h63103100000010"); //-5 < sigma < 5
    cuts.AddCutPCMCalo("12410a13","00200009227000008250400000","244668105b012200000","0h63103100000010"); //-3 < sigma < 5
    cuts.AddCutPCMCalo("12410a13","00200009327400008250400000","244668105b012200000","0h63103100000010"); //1, -10, 0.4, 3
    cuts.AddCutPCMCalo("12410a13","00200009367400008250400000","244668105b012200000","0h63103100000010"); //2, 0.5, 0.4, 3
  } else if(trainConfig == 632){//dE/dx pi-line variation and single pT variation and 2D triangular chi2 and psi pair
    cuts.AddCutPCMCalo("12410a13","00200009317400008250400000","244668105b012200000","0h63103100000010"); //0, -10, 0.4, 3
    cuts.AddCutPCMCalo("12410a13","00200049327000008250400000","244668105b012200000","0h63103100000010"); //single pT 0.075 GeV/c
    cuts.AddCutPCMCalo("12410a13","00200019327000008250400000","244668105b012200000","0h63103100000010"); //single pT 0.1   GeV/c
    cuts.AddCutPCMCalo("12410a13","00200009327000008850400000","244668105b012200000","0h63103100000010"); //20 & 0.1
  } else if(trainConfig == 633){//2D triangular chi2 and psi pair
    cuts.AddCutPCMCalo("12410a13","00200009327000008260400000","244668105b012200000","0h63103100000010"); //30 & 0.05
    cuts.AddCutPCMCalo("12410a13","00200009327000008860400000","244668105b012200000","0h63103100000010"); //20 & 0.05
    cuts.AddCutPCMCalo("12410a13","00200009327000008280400000","244668105b012200000","0h63103100000010"); //30 & 0.2
    cuts.AddCutPCMCalo("12410a13","00200009327000008880400000","244668105b012200000","0h63103100000010"); //20 & 0.2
  } else if(trainConfig == 634){//min TPC clusters and //elliptic cut Qt and alpha
    cuts.AddCutPCMCalo("12410a13","00200000327000008250400000","244668105b012200000","0h63103100000010"); //0
    cuts.AddCutPCMCalo("12410a13","00200008327000008250400000","244668105b012200000","0h63103100000010"); //0.35
    cuts.AddCutPCMCalo("12410a13","00200009327000009250400000","244668105b012200000","0h63103100000010"); // qT 0.03 no quadratic
    cuts.AddCutPCMCalo("12410a13","00200009327000003250400000","244668105b012200000","0h63103100000010"); // qT 0.05 y  quadratic
  } else if(trainConfig == 635){// and min phi max phi and NL variations
    cuts.AddCutPCMCalo("12410a13","00209909327000008250400000","244668105b012200000","0h63103100000010"); //4.54 > phi > 5.59
    cuts.AddCutPCMCalo("12410a13","00200009327000008250400000","244660105b012200000","0h63103100000010"); // PHOS people NL
  } else if(trainConfig == 636){ // first set of variations CLUSTER
    cuts.AddCutPCMCalo("12410a13","00200009327000008250400000","244668105b072200000","0h63103100000010"); // min energy 0.2 GeV
    cuts.AddCutPCMCalo("12410a13","00200009327000008250400000","244668105b082200000","0h63103100000010"); // min energy 0.4 GeV
    cuts.AddCutPCMCalo("12410a13","00200009327000008250400000","244668105b022200000","0h63103100000010"); // min energy 0.5 GeV
  } else if(trainConfig == 637){ // second set of variations CLUSTER
    cuts.AddCutPCMCalo("12410a13","00200009327000008250400000","244668105b012270000","0h63103100000010"); // min/max M02
    cuts.AddCutPCMCalo("12410a13","00200009327000008250400000","244668105b012280000","0h63103100000010"); // min/max M02
    cuts.AddCutPCMCalo("12410a13","00200009327000008250400000","244668105b013200000","0h63103100000010"); // min number 3 cells
    cuts.AddCutPCMCalo("12410a13","00200009327000008250400000","244668105b014200000","0h63103100000010"); // min number 4 cells
  } else if(trainConfig == 638){ // MESON
    cuts.AddCutPCMCalo("12410a13","00200009327000008250400000","244668105b012200000","0h63403100000010"); // rapidity variation  y<0.5
    cuts.AddCutPCMCalo("12410a13","00200009327000008250400000","244668105b012200000","0h63803100000010"); // rapidity variation  y<0.25
    cuts.AddCutPCMCalo("12410a13","00200009327000008250400000","244668105b012200000","0h63106100000010"); // alpha meson variation 1 0<alpha<0.8
  } else if(trainConfig == 639){ // fourth set of variations
    cuts.AddCutPCMCalo("12410a13","00200009327000008250400000","2446681051012200000","0h63103100000010"); // tm variation
    cuts.AddCutPCMCalo("12410a13","00200009327000008250400000","2446681054012200000","0h63103100000010"); // tm variation
    cuts.AddCutPCMCalo("12410a13","00200009327000008250400000","244668105b012200000","0h63103100000000"); // min opening angle 0    -> open
    cuts.AddCutPCMCalo("12410a13","00200009327000008250400000","244668105b012200000","0h63103100000030"); // min opening angle 0.01 -> 2 cell diag

  // **********************************************************************************************************
  // ***************************** PCM-EMC configurations PbPb run 2 2018 *************************************
  // **********************************************************************************************************
  } else if (trainConfig == 750){ // EMCAL clusters - 0-100%
    cuts.AddCutPCMCalo("10910013","0dm00009ab770c00amd0404000","411790105ke30220000","0133103100000010"); //  0-100%
  } else if (trainConfig == 751){ // EMCAL clusters - 4 cent classes with 0-10 and 30-50 triggered with neutral overlap correction via N matched tracks per cluster
    cuts.AddCutPCMCalo("10130e03","0dm00009ab770c00amd0404000","411790105te30220000","0133103100000010"); //  0-10%
    cuts.AddCutPCMCalo("11310e03","0dm00009ab770c00amd0404000","411790105te30220000","0133103100000010"); // 10-30%
    cuts.AddCutPCMCalo("13530e03","0dm00009ab770c00amd0404000","411790105te30220000","0133103100000010"); // 30-50%
    cuts.AddCutPCMCalo("15910e03","0dm00009ab770c00amd0404000","411790105te30220000","0133103100000010"); // 50-90%
  } else if (trainConfig == 752){ // EMCAL clusters - 4 cent classes with 0-10 and 30-50 triggered -standard PbPb TM for secondaries
    cuts.AddCutPCMCalo("10130e03","0dm00009ab770c00amd0404000","411790105ke30220000","0133103100000010"); //  0-10%
    cuts.AddCutPCMCalo("11310e03","0dm00009ab770c00amd0404000","411790105ke30220000","0133103100000010"); // 10-30%
    cuts.AddCutPCMCalo("13530e03","0dm00009ab770c00amd0404000","411790105ke30220000","0133103100000010"); // 30-50%
    cuts.AddCutPCMCalo("15910e03","0dm00009ab770c00amd0404000","411790105ke30220000","0133103100000010"); // 50-90%
  } else if (trainConfig == 753){ // EMCAL clusters - cen4 cent classes with 0-10 and 30-50 triggered without OOB pile up correction for MC with neutral overlap correction via N matched tracks per cluster
    cuts.AddCutPCMCalo("10130053","0dm00009ab770c00amd0404000","411790105te30220000","0133103100000010"); //  0-10%
    cuts.AddCutPCMCalo("11310053","0dm00009ab770c00amd0404000","411790105te30220000","0133103100000010"); // 10-30%
    cuts.AddCutPCMCalo("13530053","0dm00009ab770c00amd0404000","411790105te30220000","0133103100000010"); // 30-50%
    cuts.AddCutPCMCalo("15910053","0dm00009ab770c00amd0404000","411790105te30220000","0133103100000010"); // 50-90%
  } else if (trainConfig == 754){ // EMCAL clusters - cen4 cent classes with 0-10 and 30-50 triggered without OOB pile up correction for MC -standard PbPb TM for secondaries
    cuts.AddCutPCMCalo("10130053","0dm00009ab770c00amd0404000","411790105ke30220000","0133103100000010"); //  0-10%
    cuts.AddCutPCMCalo("11310053","0dm00009ab770c00amd0404000","411790105ke30220000","0133103100000010"); // 10-30%
    cuts.AddCutPCMCalo("13530053","0dm00009ab770c00amd0404000","411790105ke30220000","0133103100000010"); // 30-50%
    cuts.AddCutPCMCalo("15910053","0dm00009ab770c00amd0404000","411790105ke30220000","0133103100000010"); // 50-90%
  } else if (trainConfig == 755){ // EMCAL clusters - cen4 cent classes with 0-10 and 30-50 triggered without OOB pile up correction for MC with added Signal! with neutral overlap correction via N matched tracks per cluster
    cuts.AddCutPCMCalo("10130023","0dm00009ab770c00amd0404000","411790105te30220000","0133103100000010"); //  0-10%
    cuts.AddCutPCMCalo("11310023","0dm00009ab770c00amd0404000","411790105te30220000","0133103100000010"); // 10-30%
    cuts.AddCutPCMCalo("13530023","0dm00009ab770c00amd0404000","411790105te30220000","0133103100000010"); // 30-50%
    cuts.AddCutPCMCalo("15910023","0dm00009ab770c00amd0404000","411790105te30220000","0133103100000010"); // 50-90%
  } else if (trainConfig == 756){ // EMCAL clusters - cen4 cent classes with 0-10 and 30-50 triggered without OOB pile up correction for MC with added Signal! -standard PbPb TM for secondaries
    cuts.AddCutPCMCalo("10130023","0dm00009ab770c00amd0404000","411790105ke30220000","0133103100000010"); //  0-10%
    cuts.AddCutPCMCalo("11310023","0dm00009ab770c00amd0404000","411790105ke30220000","0133103100000010"); // 10-30%
    cuts.AddCutPCMCalo("13530023","0dm00009ab770c00amd0404000","411790105ke30220000","0133103100000010"); // 30-50%
    cuts.AddCutPCMCalo("15910023","0dm00009ab770c00amd0404000","411790105ke30220000","0133103100000010"); // 50-90%
  } else if (trainConfig == 757){ // EMCAL clusters - cen4 cent classes with 0-10 and 30-50 triggered with neutral overlap correction NonLin like
    cuts.AddCutPCMCalo("10130e03","0dm00009ab770c00amd0404000","411790105we30220000","0133103100000010"); //  0-10%
    cuts.AddCutPCMCalo("11310e03","0dm00009ab770c00amd0404000","411790105we30220000","0133103100000010"); // 10-30%
    cuts.AddCutPCMCalo("13530e03","0dm00009ab770c00amd0404000","411790105we30220000","0133103100000010"); // 30-50%
    cuts.AddCutPCMCalo("15910e03","0dm00009ab770c00amd0404000","411790105we30220000","0133103100000010"); // 50-90%
  } else if (trainConfig == 758){ // EMCAL clusters - cen4 cent classes with 0-10 and 30-50 triggered without OOB pile up correction for MC with added Signal! with neutral overlap correction NonLin like
    cuts.AddCutPCMCalo("10130023","0dm00009ab770c00amd0404000","411790105we30220000","0133103100000010"); //  0-10%
    cuts.AddCutPCMCalo("11310023","0dm00009ab770c00amd0404000","411790105we30220000","0133103100000010"); // 10-30%
    cuts.AddCutPCMCalo("13530023","0dm00009ab770c00amd0404000","411790105we30220000","0133103100000010"); // 30-50%
    cuts.AddCutPCMCalo("15910023","0dm00009ab770c00amd0404000","411790105we30220000","0133103100000010"); // 50-90%
  } else if (trainConfig == 759){ // EMCAL clusters - cen4 cent classes with 0-10 and 30-50 triggered without OOB pile up correction for MC - with neutral overlap correction NonLin like
    cuts.AddCutPCMCalo("10130053","0dm00009ab770c00amd0404000","411790105we30220000","0133103100000010"); //  0-10%
    cuts.AddCutPCMCalo("11310053","0dm00009ab770c00amd0404000","411790105we30220000","0133103100000010"); // 10-30%
    cuts.AddCutPCMCalo("13530053","0dm00009ab770c00amd0404000","411790105we30220000","0133103100000010"); // 30-50%
    cuts.AddCutPCMCalo("15910053","0dm00009ab770c00amd0404000","411790105we30220000","0133103100000010"); // 50-90%
  } else if (trainConfig == 760){ // EMCAL clusters - cen4 cent classes with 0-10 and 30-50 triggered with neutral overlap correction charged particle density
    cuts.AddCutPCMCalo("10130e03","0dm00009ab770c00amd0404000","411790105ve30220000","0133103100000010"); //  0-10%
    cuts.AddCutPCMCalo("11310e03","0dm00009ab770c00amd0404000","411790105ve30220000","0133103100000010"); // 10-30%
    cuts.AddCutPCMCalo("13530e03","0dm00009ab770c00amd0404000","411790105ve30220000","0133103100000010"); // 30-50%
    cuts.AddCutPCMCalo("15910e03","0dm00009ab770c00amd0404000","411790105ve30220000","0133103100000010"); // 50-90%
  } else if (trainConfig == 761){ // EMCAL clusters - cen4 cent classes with 0-10 and 30-50 triggered without OOB pile up correction for MC with added Signal! with neutral overlap correction charged particle density
    cuts.AddCutPCMCalo("10130023","0dm00009ab770c00amd0404000","411790105ve30220000","0133103100000010"); //  0-10%
    cuts.AddCutPCMCalo("11310023","0dm00009ab770c00amd0404000","411790105ve30220000","0133103100000010"); // 10-30%
    cuts.AddCutPCMCalo("13530023","0dm00009ab770c00amd0404000","411790105ve30220000","0133103100000010"); // 30-50%
    cuts.AddCutPCMCalo("15910023","0dm00009ab770c00amd0404000","411790105ve30220000","0133103100000010"); // 50-90%
  } else if (trainConfig == 762){ // EMCAL clusters - cen4 cent classes with 0-10 and 30-50 triggered without OOB pile up correction for MC - with neutral overlap correction charged particle density
    cuts.AddCutPCMCalo("10130053","0dm00009ab770c00amd0404000","411790105ve30220000","0133103100000010"); //  0-10%
    cuts.AddCutPCMCalo("11310053","0dm00009ab770c00amd0404000","411790105ve30220000","0133103100000010"); // 10-30%
    cuts.AddCutPCMCalo("13530053","0dm00009ab770c00amd0404000","411790105ve30220000","0133103100000010"); // 30-50%
    cuts.AddCutPCMCalo("15910053","0dm00009ab770c00amd0404000","411790105ve30220000","0133103100000010"); // 50-90%
  } else if (trainConfig == 764){ // EMCAL+DCal clusters - 13 TeV PCM-EMC cuts
    cuts.AddCutPCMCalo("10130e03","0dm00009ab770c00amd0404000","411790105fe30220000","0r63103100000010"); //  0-10%
    cuts.AddCutPCMCalo("11310e03","0dm00009ab770c00amd0404000","411790105fe30220000","0r63103100000010"); // 10-30%
    cuts.AddCutPCMCalo("13530e03","0dm00009ab770c00amd0404000","411790105fe30220000","0r63103100000010"); // 30-50%
    cuts.AddCutPCMCalo("15910e03","0dm00009ab770c00amd0404000","411790105fe30220000","0r63103100000010"); // 50-90%
  } else if (trainConfig == 771){ // EMCAL+DCal clusters - cent
    cuts.AddCutPCMCalo("10110013","0dm00009ab770c00amd0404000","411790105fe30220000","0r63103100000010"); // 00-10%
    cuts.AddCutPCMCalo("30110013","0dm00009ab770c00amd0404000","411790105fe30220000","0r63103100000010"); // 00-05%
    cuts.AddCutPCMCalo("31210013","0dm00009ab770c00amd0404000","411790105fe30220000","0r63103100000010"); // 05-10%
  } else if (trainConfig == 772){ // EMCAL+DCal clusters - semi-central
    cuts.AddCutPCMCalo("11210013","0dm00009ab770c00amd0404000","411790105fe30220000","0r63103100000010"); // 10-20%
    cuts.AddCutPCMCalo("12310013","0dm00009ab770c00amd0404000","411790105fe30220000","0r63103100000010"); // 20-30%
    cuts.AddCutPCMCalo("13410013","0dm00009ab770c00amd0404000","411790105fe30220000","0r63103100000010"); // 30-40%
    cuts.AddCutPCMCalo("12410013","0dm00009ab770c00amd0404000","411790105fe30220000","0r63103100000010"); // 20-40%
  } else if (trainConfig == 773){ // EMCAL+DCal clusters - semi peripheral
    cuts.AddCutPCMCalo("14510013","0dm00009ab770c00amd0404000","411790105fe30220000","0r63103100000010"); // 40-50%
    cuts.AddCutPCMCalo("14610013","0dm00009ab770c00amd0404000","411790105fe30220000","0r63103100000010"); // 40-60%
    cuts.AddCutPCMCalo("15610013","0dm00009ab770c00amd0404000","411790105fe30220000","0r63103100000010"); // 50-60%
  } else if (trainConfig == 774){ // EMCAL+DCal clusters - peripheral
    cuts.AddCutPCMCalo("16710013","0dm00009ab770c00amd0404000","411790105fe30220000","0r63103100000010"); // 60-70%
    cuts.AddCutPCMCalo("17810013","0dm00009ab770c00amd0404000","411790105fe30220000","0r63103100000010"); // 70-80%
    cuts.AddCutPCMCalo("18910013","0dm00009ab770c00amd0404000","411790105fe30220000","0r63103100000010"); // 80-90%
    cuts.AddCutPCMCalo("16810013","0dm00009ab770c00amd0404000","411790105fe30220000","0r63103100000010"); // 60-80%
  // Same as above but with neutral energy overlap correction for EMCal!
  } else if (trainConfig == 775){ // EMCAL+DCal clusters - cent
    cuts.AddCutPCMCalo("10110013","0dm00009ab770c00amd0404000","411790105te30220000","0r63103100000010"); // 00-10%
    cuts.AddCutPCMCalo("30110013","0dm00009ab770c00amd0404000","411790105te30220000","0r63103100000010"); // 00-05%
    cuts.AddCutPCMCalo("31210013","0dm00009ab770c00amd0404000","411790105te30220000","0r63103100000010"); // 05-10%
  } else if (trainConfig == 776){ // EMCAL+DCal clusters - semi-central
    cuts.AddCutPCMCalo("11210013","0dm00009ab770c00amd0404000","411790105te30220000","0r63103100000010"); // 10-20%
    cuts.AddCutPCMCalo("12310013","0dm00009ab770c00amd0404000","411790105te30220000","0r63103100000010"); // 20-30%
    cuts.AddCutPCMCalo("13410013","0dm00009ab770c00amd0404000","411790105te30220000","0r63103100000010"); // 30-40%
    cuts.AddCutPCMCalo("12410013","0dm00009ab770c00amd0404000","411790105te30220000","0r63103100000010"); // 20-40%
  } else if (trainConfig == 777){ // EMCAL+DCal clusters - semi peripheral
    cuts.AddCutPCMCalo("14510013","0dm00009ab770c00amd0404000","411790105te30220000","0r63103100000010"); // 40-50%
    cuts.AddCutPCMCalo("14610013","0dm00009ab770c00amd0404000","411790105te30220000","0r63103100000010"); // 40-60%
    cuts.AddCutPCMCalo("15610013","0dm00009ab770c00amd0404000","411790105te30220000","0r63103100000010"); // 50-60%
  } else if (trainConfig == 778){ // EMCAL+DCal clusters - peripheral
    cuts.AddCutPCMCalo("16710013","0dm00009ab770c00amd0404000","411790105te30220000","0r63103100000010"); // 60-70%
    cuts.AddCutPCMCalo("17810013","0dm00009ab770c00amd0404000","411790105te30220000","0r63103100000010"); // 70-80%
    cuts.AddCutPCMCalo("18910013","0dm00009ab770c00amd0404000","411790105te30220000","0r63103100000010"); // 80-90%
    cuts.AddCutPCMCalo("16810013","0dm00009ab770c00amd0404000","411790105te30220000","0r63103100000010"); // 60-80%
  } else if (trainConfig == 779){ // EMCAL+DCal clusters - with cent and semi cent trigger with OOB pileup cut
    cuts.AddCutPCMCalo("10130e13","0dm00009ab770c00amd0404000","411790105te30220000","0r63103100000010"); // 00-10%
    cuts.AddCutPCMCalo("11310e13","0dm00009ab770c00amd0404000","411790105te30220000","0r63103100000010"); // 10-30%
    cuts.AddCutPCMCalo("13530e13","0dm00009ab770c00amd0404000","411790105te30220000","0r63103100000010"); // 30-50%
    cuts.AddCutPCMCalo("15910e13","0dm00009ab770c00amd0404000","411790105te30220000","0r63103100000010"); // 50-90%
  } else if (trainConfig == 780){ // EMCAL+DCal clusters - with cent and semi cent trigger without OOB pileup cut
    cuts.AddCutPCMCalo("10130013","0dm00009ab770c00amd0404000","411790105te30220000","0r63103100000010"); // 00-10%
    cuts.AddCutPCMCalo("11310013","0dm00009ab770c00amd0404000","411790105te30220000","0r63103100000010"); // 10-30%
    cuts.AddCutPCMCalo("13530013","0dm00009ab770c00amd0404000","411790105te30220000","0r63103100000010"); // 30-50%
    cuts.AddCutPCMCalo("15910013","0dm00009ab770c00amd0404000","411790105te30220000","0r63103100000010"); // 50-90%
  } else if (trainConfig == 781){ // EMCAL+DCal clusters - with cent and semi cent trigger with OOB pileup cut w/o TM
    cuts.AddCutPCMCalo("10130e13","0dm00009ab770c00amd0404000","411790105ke30220000","0r63103100000010"); // 00-10%
    cuts.AddCutPCMCalo("11310e13","0dm00009ab770c00amd0404000","411790105ke30220000","0r63103100000010"); // 10-30%
    cuts.AddCutPCMCalo("13530e13","0dm00009ab770c00amd0404000","411790105ke30220000","0r63103100000010"); // 30-50%
    cuts.AddCutPCMCalo("15910e13","0dm00009ab770c00amd0404000","411790105ke30220000","0r63103100000010"); // 50-90%
  } else if (trainConfig == 782){ // EMCAL+DCal clusters - with cent and semi cent trigger without OOB pileup cut w/o TM
    cuts.AddCutPCMCalo("10130013","0dm00009ab770c00amd0404000","411790105ke30220000","0r63103100000010"); // 00-10%
    cuts.AddCutPCMCalo("11310013","0dm00009ab770c00amd0404000","411790105ke30220000","0r63103100000010"); // 10-30%
    cuts.AddCutPCMCalo("13530013","0dm00009ab770c00amd0404000","411790105ke30220000","0r63103100000010"); // 30-50%
    cuts.AddCutPCMCalo("15910013","0dm00009ab770c00amd0404000","411790105ke30220000","0r63103100000010"); // 50-90%
  } else if (trainConfig == 783){ // EMCAL+DCal clusters - cent without TM
    cuts.AddCutPCMCalo("30110013","0dm00009ab770c00amd0404000","411790105ke30220000","0r63103100000010"); // 00-05%
    cuts.AddCutPCMCalo("31210013","0dm00009ab770c00amd0404000","411790105ke30220000","0r63103100000010"); // 05-10%
  } else if (trainConfig == 784){ // EMCAL+DCal clusters - semi-central without TM
    cuts.AddCutPCMCalo("11210013","0dm00009ab770c00amd0404000","411790105ke30220000","0r63103100000010"); // 10-20%
    cuts.AddCutPCMCalo("12310013","0dm00009ab770c00amd0404000","411790105ke30220000","0r63103100000010"); // 20-30%
    cuts.AddCutPCMCalo("13410013","0dm00009ab770c00amd0404000","411790105ke30220000","0r63103100000010"); // 30-40%
  } else if (trainConfig == 785){ // EMCAL+DCal clusters - semi peripheral without TM
    cuts.AddCutPCMCalo("14510013","0dm00009ab770c00amd0404000","411790105ke30220000","0r63103100000010"); // 40-50%
    cuts.AddCutPCMCalo("15610013","0dm00009ab770c00amd0404000","411790105ke30220000","0r63103100000010"); // 50-60%
  } else if (trainConfig == 786){ // EMCAL+DCal clusters - peripheral without TM
    cuts.AddCutPCMCalo("16710013","0dm00009ab770c00amd0404000","411790105ke30220000","0r63103100000010"); // 60-70%
    cuts.AddCutPCMCalo("17810013","0dm00009ab770c00amd0404000","411790105ke30220000","0r63103100000010"); // 70-80%
    cuts.AddCutPCMCalo("18910013","0dm00009ab770c00amd0404000","411790105ke30220000","0r63103100000010"); // 80-90%
  } else if (trainConfig == 787){ // EMCAL clusters - cen4 cent classes with 0-10 and 30-50 triggered with neutral overlap correction maxwell boltzmann
    cuts.AddCutPCMCalo("10130e03","0dm00009ab770c00amd0404000","411790105xe30220000","0133103100000010"); //  0-10%
    cuts.AddCutPCMCalo("11310e03","0dm00009ab770c00amd0404000","411790105xe30220000","0133103100000010"); // 10-30%
    cuts.AddCutPCMCalo("13530e03","0dm00009ab770c00amd0404000","411790105xe30220000","0133103100000010"); // 30-50%
    cuts.AddCutPCMCalo("15910e03","0dm00009ab770c00amd0404000","411790105xe30220000","0133103100000010"); // 50-90%
  } else if (trainConfig == 788){ // EMCAL clusters - cen4 cent classes with 0-10 and 30-50 triggered without OOB pile up correction for MC with added Signal! with neutral overlap correction maxwell boltzmann
    cuts.AddCutPCMCalo("10130023","0dm00009ab770c00amd0404000","411790105xe30220000","0133103100000010"); //  0-10%
    cuts.AddCutPCMCalo("11310023","0dm00009ab770c00amd0404000","411790105xe30220000","0133103100000010"); // 10-30%
    cuts.AddCutPCMCalo("13530023","0dm00009ab770c00amd0404000","411790105xe30220000","0133103100000010"); // 30-50%
    cuts.AddCutPCMCalo("15910023","0dm00009ab770c00amd0404000","411790105xe30220000","0133103100000010"); // 50-90%
  } else if (trainConfig == 789){ // EMCAL clusters - cen4 cent classes with 0-10 and 30-50 triggered without OOB pile up correction for MC - with neutral overlap correction maxwell boltzmann
    cuts.AddCutPCMCalo("10130053","0dm00009ab770c00amd0404000","411790105xe30220000","0133103100000010"); //  0-10%
    cuts.AddCutPCMCalo("11310053","0dm00009ab770c00amd0404000","411790105xe30220000","0133103100000010"); // 10-30%
    cuts.AddCutPCMCalo("13530053","0dm00009ab770c00amd0404000","411790105xe30220000","0133103100000010"); // 30-50%
    cuts.AddCutPCMCalo("15910053","0dm00009ab770c00amd0404000","411790105xe30220000","0133103100000010"); // 50-90%

  // **********************************************************************************************************
  // ***************************** PCM-PHOS       QA configurations PbPb run 2 2018 ***************************
  // **********************************************************************************************************
  } else if (trainConfig == 850){ // PHOS clusters - centrality selection for PbPb
    cuts.AddCutPCMCalo("10930013","00200009327000008250400000","24466000ha082200000","0133103100000050"); //  0-100%
  } else if (trainConfig == 851){ // PHOS clusters - centrality selection for PbPb
    cuts.AddCutPCMCalo("10130013","00200009327000008250400000","24466000ha082200000","0133103100000050"); //  0-10%
    cuts.AddCutPCMCalo("11330013","00200009327000008250400000","24466000ha082200000","0133103100000050"); // 10-30%
    cuts.AddCutPCMCalo("13530013","00200009327000008250400000","24466000ha082200000","0133103100000050"); // 30-50%
    cuts.AddCutPCMCalo("15930013","00200009327000008250400000","24466000ha082200000","0133103100000050"); // 50-90%
  } else if (trainConfig == 852){ // PHOS clusters - centrality selection for PbPb
    cuts.AddCutPCMCalo("10130a13","00200009327000008250400000","24466000ha082200000","0h33103100000010"); //
    cuts.AddCutPCMCalo("11230a13","00200009327000008250400000","24466000ha082200000","0h33103100000010"); //
    cuts.AddCutPCMCalo("12430a13","00200009327000008250400000","24466000ha082200000","0h33103100000010"); //
    cuts.AddCutPCMCalo("14630a13","00200009327000008250400000","24466000ha082200000","0h33103100000010"); //
    cuts.AddCutPCMCalo("16830a13","00200009327000008250400000","24466000ha082200000","0h33103100000010"); //
  } else if (trainConfig == 853){ // PHOS clusters - centrality selection for PbPb
    cuts.AddCutPCMCalo("10130a13","00200009327000008250400000","24466000ha082200000","0h33103100000010"); //
    cuts.AddCutPCMCalo("11310a13","00200009327000008250400000","24466000ha082200000","0h33103100000010"); //
    cuts.AddCutPCMCalo("13530a13","00200009327000008250400000","24466000ha082200000","0h33103100000010"); //
    cuts.AddCutPCMCalo("15910a13","00200009327000008250400000","24466000ha082200000","0h33103100000010"); //
  } else if (trainConfig == 854){ // PHOS clusters - centrality selection for PbPb
    cuts.AddCutPCMCalo("10130a13","00200009f9730000dge0400000","24466810ha082200000","0h33103100000010"); //
    cuts.AddCutPCMCalo("11310a13","00200009f9730000dge0400000","24466810ha082200000","0h33103100000010"); //
    cuts.AddCutPCMCalo("13530a13","00200009f9730000dge0400000","24466810ha082200000","0h33103100000010"); //
    cuts.AddCutPCMCalo("15910a13","00200009f9730000dge0400000","24466810ha082200000","0h33103100000010"); //
  } else if (trainConfig == 855){ // PHOS clusters - centrality selection for PbPb
    cuts.AddCutPCMCalo("30130a13","00200009f9730000dge0400000","24466810ha082200000","0h43103100000010"); //
    cuts.AddCutPCMCalo("31230a13","00200009f9730000dge0400000","24466810ha082200000","0h43103100000010"); //
    cuts.AddCutPCMCalo("11210a13","00200009f9730000dge0400000","24466810ha082200000","0h43103100000010"); //
    cuts.AddCutPCMCalo("12310a13","00200009f9730000dge0400000","24466810ha082200000","0h43103100000010"); //
  } else if (trainConfig == 856){ // PHOS clusters - centrality selection for PbPb
    cuts.AddCutPCMCalo("13430a13","00200009f9730000dge0400000","24466810ha082200000","0h43103100000010"); //
    cuts.AddCutPCMCalo("14530a13","00200009f9730000dge0400000","24466810ha082200000","0h43103100000010"); //
    cuts.AddCutPCMCalo("15610a13","00200009f9730000dge0400000","24466810ha082200000","0h43103100000010"); //
    cuts.AddCutPCMCalo("16710a13","00200009f9730000dge0400000","24466810ha082200000","0h43103100000010"); //
    cuts.AddCutPCMCalo("17810a13","00200009f9730000dge0400000","24466810ha082200000","0h43103100000010"); //
    cuts.AddCutPCMCalo("18910a13","00200009f9730000dge0400000","24466810ha082200000","0h43103100000010"); //
  } else if (trainConfig == 857){ // PHOS clusters - centrality selection for PbPb
    cuts.AddCutPCMCalo("30130a13","0dm00009f9730000dge0404000","24466810ha082200000","0h43103100000010"); //
    cuts.AddCutPCMCalo("31230a13","0dm00009f9730000dge0404000","24466810ha082200000","0h43103100000010"); //
    cuts.AddCutPCMCalo("11210a13","0dm00009f9730000dge0404000","24466810ha082200000","0h43103100000010"); //
    cuts.AddCutPCMCalo("12310a13","0dm00009f9730000dge0404000","24466810ha082200000","0h43103100000010"); //
  } else if (trainConfig == 858){ // PHOS clusters - centrality selection for PbPb
    cuts.AddCutPCMCalo("13430a13","0dm00009f9730000dge0404000","24466810ha082200000","0h43103100000010"); //
    cuts.AddCutPCMCalo("14530a13","0dm00009f9730000dge0404000","24466810ha082200000","0h43103100000010"); //
    cuts.AddCutPCMCalo("15610a13","0dm00009f9730000dge0404000","24466810ha082200000","0h43103100000010"); //
    cuts.AddCutPCMCalo("16710a13","0dm00009f9730000dge0404000","24466810ha082200000","0h43103100000010"); //
    cuts.AddCutPCMCalo("17810a13","0dm00009f9730000dge0404000","24466810ha082200000","0h43103100000010"); //
    cuts.AddCutPCMCalo("18910a13","0dm00009f9730000dge0404000","24466810ha082200000","0h43103100000010"); //
  } else if (trainConfig == 860){ // PHOS clusters - trigger 20181018
    cuts.AddCutPCMCalo("30162a13","00200009f9730000dge0400000","24466810ha082200000","0h43103100000010"); //
    cuts.AddCutPCMCalo("31262a13","00200009f9730000dge0400000","24466810ha082200000","0h43103100000010"); //
    cuts.AddCutPCMCalo("11262a13","00200009f9730000dge0400000","24466810ha082200000","0h43103100000010"); //
    cuts.AddCutPCMCalo("12362a13","00200009f9730000dge0400000","24466810ha082200000","0h43103100000010"); //
    
  } else if (trainConfig == 880){ // PHOS clusters - similar to EMCal 780
    cuts.AddCutPCMCalo("10130e03","0dm00009f9730000dge0404000","24466810ha082200000","0133103100000010"); // 00-10%
    cuts.AddCutPCMCalo("11310e03","0dm00009f9730000dge0404000","24466810ha082200000","0133103100000010"); // 10-30%
    cuts.AddCutPCMCalo("13530e03","0dm00009f9730000dge0404000","24466810ha082200000","0133103100000010"); // 30-50%
    cuts.AddCutPCMCalo("15910e03","0dm00009f9730000dge0404000","24466810ha082200000","0133103100000010"); // 50-90%
  } else if (trainConfig == 881){ // PHOS clusters - similar to EMCal 780 only 0-10 triggered
    cuts.AddCutPCMCalo("10130e03","0dm00009f9730000dge0404000","24466810ha082200000","0133103100000010"); // 00-10%
  } else if (trainConfig == 882){ // PHOS clusters - similar to EMCal 780 only 30-50 triggered
    cuts.AddCutPCMCalo("13530e03","0dm00009f9730000dge0404000","24466810ha082200000","0133103100000010"); // 30-50%
  // **********************************************************************************************************
  // ***************************** PCM-PHOS  HBT configurations PbPb run 2 2018 *******************************
  // **********************************************************************************************************
  // if(trainConfig >= 950 && trainConfig <= 1000) <--- RESERVED FOR HBT STUDY ENABLING
  } else if (trainConfig == 950){ // PHOS clusters - centrality selection for PbPb
    cuts.AddCutPCMCalo("10910a13","00200009327000008250400000","24466000ha082200000","0h33103100000010"); //
    cuts.AddCutPCMCalo("10130a13","00200009327000008250400000","24466000ha082200000","0h33103100000010"); //
    cuts.AddCutPCMCalo("13530a13","00200009327000008250400000","24466000ha082200000","0h33103100000010"); //
  } else if (trainConfig == 951){ // PHOS clusters - centrality selection for PbPb
    cuts.AddCutPCMCalo("10130a13","00200009f9730000dge0404000","24466810ha082200000","0h33103100000010"); //
    cuts.AddCutPCMCalo("13530a13","00200009f9730000dge0404000","24466810ha082200000","0h33103100000010"); //
  } else if (trainConfig == 952){ // PHOS clusters - centrality selection for PbPb
    cuts.AddCutPCMCalo("10910a13","00600009a27000006250800000","24466810ha082200000","0h33103100000010"); //
  } else if (trainConfig == 953){ // PHOS clusters - centrality selection for PbPb - 400 MeV cluster thresholds
    cuts.AddCutPCMCalo("10130e03","0dm00009f9730000dge0404000","24466810ha082200000","0h33103100000010"); //
    cuts.AddCutPCMCalo("11310e03","0dm00009f9730000dge0404000","24466810ha082200000","0h33103100000010"); //
    cuts.AddCutPCMCalo("13530e03","0dm00009f9730000dge0404000","24466810ha082200000","0h33103100000010"); //
    cuts.AddCutPCMCalo("15910e03","0dm00009f9730000dge0404000","24466810ha082200000","0h33103100000010"); //
  } else if (trainConfig == 954){ // PHOS clusters - centrality selection for PbPb - 300 MeV cluster thresholds
    cuts.AddCutPCMCalo("10130e03","0dm00009f9730000dge0404000","24466810ha01ee00000","0h33103100000010"); //
    cuts.AddCutPCMCalo("13530e03","0dm00009f9730000dge0404000","24466810ha01ee00000","0h33103100000010"); //
  } else if (trainConfig == 955){ // PHOS clusters - centrality selection for PbPb - 100 MeV cluster thresholds
    cuts.AddCutPCMCalo("10130e03","0dm00009f9730000dge0404000","24466810ha09ee00000","0h33103100000010"); //
    cuts.AddCutPCMCalo("11310e03","0dm00009f9730000dge0404000","24466810ha09ee00000","0h33103100000010"); //
    cuts.AddCutPCMCalo("13530e03","0dm00009f9730000dge0404000","24466810ha09ee00000","0h33103100000010"); //
    cuts.AddCutPCMCalo("15910e03","0dm00009f9730000dge0404000","24466810ha09ee00000","0h33103100000010"); //
  } else if (trainConfig == 956){ // PHOS clusters - centrality selection for PbPb - 300 MeV cluster thresholds - eventplane cuts check 1
    cuts.AddCutPCMCalo("10130e03","0dm00009f9730000dge0404001","24466810ha01ee00000","0h33103100000010"); //
    cuts.AddCutPCMCalo("11310e03","0dm00009f9730000dge0404001","24466810ha01ee00000","0h33103100000010"); //
    cuts.AddCutPCMCalo("13530e03","0dm00009f9730000dge0404001","24466810ha01ee00000","0h33103100000010"); //
    cuts.AddCutPCMCalo("15910e03","0dm00009f9730000dge0404001","24466810ha01ee00000","0h33103100000010"); //
  } else if (trainConfig == 957){ // PHOS clusters - centrality selection for PbPb - 300 MeV cluster thresholds - eventplane cuts check 2
    cuts.AddCutPCMCalo("10130e03","0dm00009f9730000dge0404002","24466810ha01ee00000","0h33103100000010"); //
    cuts.AddCutPCMCalo("11310e03","0dm00009f9730000dge0404002","24466810ha01ee00000","0h33103100000010"); //
    cuts.AddCutPCMCalo("13530e03","0dm00009f9730000dge0404002","24466810ha01ee00000","0h33103100000010"); //
    cuts.AddCutPCMCalo("15910e03","0dm00009f9730000dge0404002","24466810ha01ee00000","0h33103100000010"); //
  //PCM-PHOS PbPb HBT systematics studies 30-50
  } else if ( trainConfig == 960){ // Systematics (1/4) TM variations & opening angle
    cuts.AddCutPCMCalo("13530e03","0dm00009f9730000dge0404000","24466810h101ee00000","0h33103100000010"); // TM stricter eta and phi
    cuts.AddCutPCMCalo("13530e03","0dm00009f9730000dge0404000","24466810hb01ee00000","0h33103100000010"); // TM stricter phi
    cuts.AddCutPCMCalo("13530e03","0dm00009f9730000dge0404000","24466810ha01ee00000","0h33103100000000"); // min opening angle 0
    cuts.AddCutPCMCalo("13530e03","0dm00009f9730000dge0404000","24466810ha01ee00000","0h33103100000030"); // min opening angle 0.01 (2cell diag)
  } else if ( trainConfig == 961){ // Systematics (1/4) rapidity variations & alpha meson
    cuts.AddCutPCMCalo("13530e03","0dm00009f9730000dge0404000","24466810ha01ee00000","0h33403100000010"); // rapidity variation  y<0.5
    cuts.AddCutPCMCalo("13530e03","0dm00009f9730000dge0404000","24466810ha01ee00000","0h33803100000010"); // rapidity variation  y<0.25
    cuts.AddCutPCMCalo("13530e03","0dm00009f9730000dge0404000","24466810ha01ee00000","0h33106100000010"); // alpha meson
    cuts.AddCutPCMCalo("13530e03","0dm00009f9730000dge0404000","24466810ha01ee00000","0h33105100000010"); // alpha meson
  } else if ( trainConfig == 962){ // Systematics (1/4) dE/dx variations
    cuts.AddCutPCMCalo("13530e03","0dm0000919730000dge0404000","24466810ha01ee00000","0h33103100000010"); // -5 < sigma < 5
    cuts.AddCutPCMCalo("13530e03","0dm0000929730000dge0404000","24466810ha01ee00000","0h33103100000010"); // -3 < sigma < 5
    cuts.AddCutPCMCalo("13530e03","0dm00049f9730000dge0404000","24466810ha01ee00000","0h33103100000010"); // single pT 0.075 GeV/c
    cuts.AddCutPCMCalo("13530e03","0dm00019f9730000dge0404000","24466810ha01ee00000","0h33103100000010"); // single pT 0.1   GeV/c
  } else if ( trainConfig == 963){ // Systematics (1/4) 2D triangular chi2 and psi pair, min TPC clusters
    cuts.AddCutPCMCalo("13530e03","0dm00009f9730000d250404000","24466810ha01ee00000","0h33103100000010"); // chi2 and psi pair
    cuts.AddCutPCMCalo("13530e03","0dm00009f9730000d880404000","24466810ha01ee00000","0h33103100000010"); // chi2 and psi pair
    cuts.AddCutPCMCalo("13530e03","0dm00000f9730000dge0404000","24466810ha01ee00000","0h33103100000010"); // min TPC clusters
    cuts.AddCutPCMCalo("13530e03","0dm00008f9730000dge0404000","24466810ha01ee00000","0h33103100000010"); // min TPC clusters
    //PCM-PHOS PbPb HBT systematics studies 0-10
  } else if ( trainConfig == 964){ // Systematics (1/4) TM variations & opening angle
    cuts.AddCutPCMCalo("10130e03","0dm00009f9730000dge0404000","24466810h101ee00000","0h33103100000010"); // TM stricter eta and phi
    cuts.AddCutPCMCalo("10130e03","0dm00009f9730000dge0404000","24466810hb01ee00000","0h33103100000010"); // TM stricter phi
    cuts.AddCutPCMCalo("10130e03","0dm00009f9730000dge0404000","24466810ha01ee00000","0h33103100000000"); // min opening angle 0
    cuts.AddCutPCMCalo("10130e03","0dm00009f9730000dge0404000","24466810ha01ee00000","0h33103100000030"); // min opening angle 0.01 (2cell diag)
  } else if ( trainConfig == 965){ // Systematics (1/4) rapidity variations & alpha meson
    cuts.AddCutPCMCalo("10130e03","0dm00009f9730000dge0404000","24466810ha01ee00000","0h33403100000010"); // rapidity variation  y<0.5
    cuts.AddCutPCMCalo("10130e03","0dm00009f9730000dge0404000","24466810ha01ee00000","0h33803100000010"); // rapidity variation  y<0.25
    cuts.AddCutPCMCalo("10130e03","0dm00009f9730000dge0404000","24466810ha01ee00000","0h33106100000010"); // alpha meson
    cuts.AddCutPCMCalo("10130e03","0dm00009f9730000dge0404000","24466810ha01ee00000","0h33105100000010"); // alpha meson
  } else if ( trainConfig == 966){ // Systematics (1/4) dE/dx variations
    cuts.AddCutPCMCalo("10130e03","0dm0000919730000dge0404000","24466810ha01ee00000","0h33103100000010"); // -5 < sigma < 5
    cuts.AddCutPCMCalo("10130e03","0dm0000929730000dge0404000","24466810ha01ee00000","0h33103100000010"); // -3 < sigma < 5
    cuts.AddCutPCMCalo("10130e03","0dm00049f9730000dge0404000","24466810ha01ee00000","0h33103100000010"); // single pT 0.075 GeV/c
    cuts.AddCutPCMCalo("10130e03","0dm00019f9730000dge0404000","24466810ha01ee00000","0h33103100000010"); // single pT 0.1   GeV/c
  } else if ( trainConfig == 967){ // Systematics (1/4) 2D triangular chi2 and psi pair, min TPC clusters
    cuts.AddCutPCMCalo("10130e03","0dm00009f9730000d250404000","24466810ha01ee00000","0h33103100000010"); // chi2 and psi pair
    cuts.AddCutPCMCalo("10130e03","0dm00009f9730000d880404000","24466810ha01ee00000","0h33103100000010"); // chi2 and psi pair
    cuts.AddCutPCMCalo("10130e03","0dm00000f9730000dge0404000","24466810ha01ee00000","0h33103100000010"); // min TPC clusters
    cuts.AddCutPCMCalo("10130e03","0dm00008f9730000dge0404000","24466810ha01ee00000","0h33103100000010"); // min TPC clusters

  // **********************************************************************************************************
  // ***************************** PCM-EMC configurations PbPb run 2 2018 pass 3 *******************************
  // **********************************************************************************************************
  // Standard: EMCal + DCal, only 2ndary cluster track matching, mixed event
  } else if (trainConfig == 1000){ // 4 cents
    cuts.AddCutPCMCalo("10130e03","0dm00009ab770c00amd0404000","411790105ke30220000","0133103100000010"); // 00-10%
    cuts.AddCutPCMCalo("11310e03","0dm00009ab770c00amd0404000","411790105ke30220000","0133103100000010"); // 10-30%
    cuts.AddCutPCMCalo("13530e03","0dm00009ab770c00amd0404000","411790105ke30220000","0133103100000010"); // 30-50%
    cuts.AddCutPCMCalo("15910e03","0dm00009ab770c00amd0404000","411790105ke30220000","0133103100000010"); // 50-90%
  } else if (trainConfig == 1001){ // 4 cents, no OOB Pileup correction for MC
    cuts.AddCutPCMCalo("10130053","0dm00009ab770c00amd0404000","411790105ke30220000","0133103100000010"); // 00-10%
    cuts.AddCutPCMCalo("11310053","0dm00009ab770c00amd0404000","411790105ke30220000","0133103100000010"); // 10-30%
    cuts.AddCutPCMCalo("13530053","0dm00009ab770c00amd0404000","411790105ke30220000","0133103100000010"); // 30-50%
    cuts.AddCutPCMCalo("15910053","0dm00009ab770c00amd0404000","411790105ke30220000","0133103100000010"); // 50-90%
  } else if (trainConfig == 1002){ // 4 cents, no OOB Pileup correction for MC with added signal
    cuts.AddCutPCMCalo("10130023","0dm00009ab770c00amd0404000","411790105ke30220000","0133103100000010"); // 00-10%
    cuts.AddCutPCMCalo("11310023","0dm00009ab770c00amd0404000","411790105ke30220000","0133103100000010"); // 10-30%
    cuts.AddCutPCMCalo("13530023","0dm00009ab770c00amd0404000","411790105ke30220000","0133103100000010"); // 30-50%
    cuts.AddCutPCMCalo("15910023","0dm00009ab770c00amd0404000","411790105ke30220000","0133103100000010"); // 50-90%
  } else if (trainConfig == 1003){ // 4 cents, no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutPCMCalo("10130023","0dm00009ab770c00amd0404000","411790105ke30220000","0133103100000010"); // 00-10%
    cuts.AddCutPCMCalo("11310023","0dm00009ab770c00amd0404000","411790105ke30220000","0133103100000010"); // 10-30%
    cuts.AddCutPCMCalo("13530023","0dm00009ab770c00amd0404000","411790105ke30220000","0133103100000010"); // 30-50%
    cuts.AddCutPCMCalo("15910023","0dm00009ab770c00amd0404000","411790105ke30220000","0133103100000010"); // 50-90%
  } else if (trainConfig == 1004){ // 4 cents, no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutPCMCalo("10130023","0dm00009ab770c00amd0404000","411790105ke30220000","0133103100000010"); // 00-10%
    cuts.AddCutPCMCalo("11310023","0dm00009ab770c00amd0404000","411790105ke30220000","0133103100000010"); // 10-30%
    cuts.AddCutPCMCalo("13530023","0dm00009ab770c00amd0404000","411790105ke30220000","0133103100000010"); // 30-50%
    cuts.AddCutPCMCalo("15910023","0dm00009ab770c00amd0404000","411790105ke30220000","0133103100000010"); // 50-90%
  } else if (trainConfig == 1005){ // 4 cents with NOC NonLin like data
    cuts.AddCutPCMCalo("10130e03","0dm00009ab770c00amd0404000","411790105we30220000","0133103100000010"); // 00-10%
    cuts.AddCutPCMCalo("11310e03","0dm00009ab770c00amd0404000","411790105we30220000","0133103100000010"); // 10-30%
    cuts.AddCutPCMCalo("13530e03","0dm00009ab770c00amd0404000","411790105we30220000","0133103100000010"); // 30-50%
    cuts.AddCutPCMCalo("15910e03","0dm00009ab770c00amd0404000","411790105we30220000","0133103100000010"); // 50-90%
  } else if (trainConfig == 1006){ // 4 cents, with NOC NonLin like MC no OOB Pileup correction for MC
    cuts.AddCutPCMCalo("10130053","0dm00009ab770c00amd0404000","411790105xe30220000","0133103100000010"); // 00-10%
    cuts.AddCutPCMCalo("11310053","0dm00009ab770c00amd0404000","411790105xe30220000","0133103100000010"); // 10-30%
    cuts.AddCutPCMCalo("13530053","0dm00009ab770c00amd0404000","411790105xe30220000","0133103100000010"); // 30-50%
    cuts.AddCutPCMCalo("15910053","0dm00009ab770c00amd0404000","411790105xe30220000","0133103100000010"); // 50-90%
  } else if (trainConfig == 1007){ // 4 cents, with NOC NonLin like MC no OOB Pileup correction for MC with added signal
    cuts.AddCutPCMCalo("10130023","0dm00009ab770c00amd0404000","411790105xe30220000","0133103100000010"); // 00-10%
    cuts.AddCutPCMCalo("11310023","0dm00009ab770c00amd0404000","411790105xe30220000","0133103100000010"); // 10-30%
    cuts.AddCutPCMCalo("13530023","0dm00009ab770c00amd0404000","411790105xe30220000","0133103100000010"); // 30-50%
    cuts.AddCutPCMCalo("15910023","0dm00009ab770c00amd0404000","411790105xe30220000","0133103100000010"); // 50-90%
  } else if (trainConfig == 1008){ // 4 cents, with NOC NonLin like MC no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutPCMCalo("10130023","0dm00009ab770c00amd0404000","411790105xe30220000","0133103100000010"); // 00-10%
    cuts.AddCutPCMCalo("11310023","0dm00009ab770c00amd0404000","411790105xe30220000","0133103100000010"); // 10-30%
    cuts.AddCutPCMCalo("13530023","0dm00009ab770c00amd0404000","411790105xe30220000","0133103100000010"); // 30-50%
    cuts.AddCutPCMCalo("15910023","0dm00009ab770c00amd0404000","411790105xe30220000","0133103100000010"); // 50-90%
  } else if (trainConfig == 1009){ // 4 cents, with NOC NonLin like MC no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutPCMCalo("10130023","0dm00009ab770c00amd0404000","411790105xe30220000","0133103100000010"); // 00-10%
    cuts.AddCutPCMCalo("11310023","0dm00009ab770c00amd0404000","411790105xe30220000","0133103100000010"); // 10-30%
    cuts.AddCutPCMCalo("13530023","0dm00009ab770c00amd0404000","411790105xe30220000","0133103100000010"); // 30-50%
    cuts.AddCutPCMCalo("15910023","0dm00009ab770c00amd0404000","411790105xe30220000","0133103100000010"); // 50-90%
  } else if (trainConfig == 1010){ // 4 cents with mean N matched tracks per cluster data
    cuts.AddCutPCMCalo("10130e03","0dm00009ab770c00amd0404000","411790105te30220000","0133103100000010"); // 00-10%
    cuts.AddCutPCMCalo("11310e03","0dm00009ab770c00amd0404000","411790105te30220000","0133103100000010"); // 10-30%
    cuts.AddCutPCMCalo("13530e03","0dm00009ab770c00amd0404000","411790105te30220000","0133103100000010"); // 30-50%
    cuts.AddCutPCMCalo("15910e03","0dm00009ab770c00amd0404000","411790105te30220000","0133103100000010"); // 50-90%
  } else if (trainConfig == 1011){ // 4 cents, with mean N matched tracks per cluster MC no OOB Pileup correction for MC
    cuts.AddCutPCMCalo("10130053","0dm00009ab770c00amd0404000","411790105te30220000","0133103100000010"); // 00-10%
    cuts.AddCutPCMCalo("11310053","0dm00009ab770c00amd0404000","411790105te30220000","0133103100000010"); // 10-30%
    cuts.AddCutPCMCalo("13530053","0dm00009ab770c00amd0404000","411790105te30220000","0133103100000010"); // 30-50%
    cuts.AddCutPCMCalo("15910053","0dm00009ab770c00amd0404000","411790105te30220000","0133103100000010"); // 50-90%
  } else if (trainConfig == 1012){ // 4 cents, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal
    cuts.AddCutPCMCalo("10130023","0dm00009ab770c00amd0404000","411790105te30220000","0133103100000010"); // 00-10%
    cuts.AddCutPCMCalo("11310023","0dm00009ab770c00amd0404000","411790105te30220000","0133103100000010"); // 10-30%
    cuts.AddCutPCMCalo("13530023","0dm00009ab770c00amd0404000","411790105te30220000","0133103100000010"); // 30-50%
    cuts.AddCutPCMCalo("15910023","0dm00009ab770c00amd0404000","411790105te30220000","0133103100000010"); // 50-90%
  } else if (trainConfig == 1013){ // 4 cents, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutPCMCalo("10130023","0dm00009ab770c00amd0404000","411790105te30220000","0133103100000010"); // 00-10%
    cuts.AddCutPCMCalo("11310023","0dm00009ab770c00amd0404000","411790105te30220000","0133103100000010"); // 10-30%
    cuts.AddCutPCMCalo("13530023","0dm00009ab770c00amd0404000","411790105te30220000","0133103100000010"); // 30-50%
    cuts.AddCutPCMCalo("15910023","0dm00009ab770c00amd0404000","411790105te30220000","0133103100000010"); // 50-90%
  } else if (trainConfig == 1014){ // 4 cents, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutPCMCalo("10130023","0dm00009ab770c00amd0404000","411790105te30220000","0133103100000010"); // 00-10%
    cuts.AddCutPCMCalo("11310023","0dm00009ab770c00amd0404000","411790105te30220000","0133103100000010"); // 10-30%
    cuts.AddCutPCMCalo("13530023","0dm00009ab770c00amd0404000","411790105te30220000","0133103100000010"); // 30-50%
    cuts.AddCutPCMCalo("15910023","0dm00009ab770c00amd0404000","411790105te30220000","0133103100000010"); // 50-90%
  } else if (trainConfig == 1015){ // 4 cents with mean N matched tracks per cluster data
    cuts.AddCutPCMCalo("10130e03","0dm00009ab770c00amd0404000","411790105ye30220000","0133103100000010"); // 00-10%
    cuts.AddCutPCMCalo("11310e03","0dm00009ab770c00amd0404000","411790105ye30220000","0133103100000010"); // 10-30%
    cuts.AddCutPCMCalo("13530e03","0dm00009ab770c00amd0404000","411790105ye30220000","0133103100000010"); // 30-50%
    cuts.AddCutPCMCalo("15910e03","0dm00009ab770c00amd0404000","411790105ye30220000","0133103100000010"); // 50-90%
  } else if (trainConfig == 1016){ // 4 cents, with mean N matched tracks per cluster MC no OOB Pileup correction for MC
    cuts.AddCutPCMCalo("10130053","0dm00009ab770c00amd0404000","411790105ye30220000","0133103100000010"); // 00-10%
    cuts.AddCutPCMCalo("11310053","0dm00009ab770c00amd0404000","411790105ye30220000","0133103100000010"); // 10-30%
    cuts.AddCutPCMCalo("13530053","0dm00009ab770c00amd0404000","411790105ye30220000","0133103100000010"); // 30-50%
    cuts.AddCutPCMCalo("15910053","0dm00009ab770c00amd0404000","411790105ye30220000","0133103100000010"); // 50-90%
  // with rotation background
  } else if (trainConfig == 1050){ // 4 cents
    cuts.AddCutPCMCalo("10130e03","0dm00009ab770c00amd0404000","411790105ke30220000","0s331031000000d0"); // 00-10%
    cuts.AddCutPCMCalo("11310e03","0dm00009ab770c00amd0404000","411790105ke30220000","0s331031000000d0"); // 10-30%
    cuts.AddCutPCMCalo("13530e03","0dm00009ab770c00amd0404000","411790105ke30220000","0s331031000000d0"); // 30-50%
    cuts.AddCutPCMCalo("15910e03","0dm00009ab770c00amd0404000","411790105ke30220000","0s331031000000d0"); // 50-90%
  } else if (trainConfig == 1051){ // 4 cents, no OOB Pileup correction for MC
    cuts.AddCutPCMCalo("10130053","0dm00009ab770c00amd0404000","411790105ke30220000","0s331031000000d0"); // 00-10%
    cuts.AddCutPCMCalo("11310053","0dm00009ab770c00amd0404000","411790105ke30220000","0s331031000000d0"); // 10-30%
    cuts.AddCutPCMCalo("13530053","0dm00009ab770c00amd0404000","411790105ke30220000","0s331031000000d0"); // 30-50%
    cuts.AddCutPCMCalo("15910053","0dm00009ab770c00amd0404000","411790105ke30220000","0s331031000000d0"); // 50-90%
  } else if (trainConfig == 1052){ // 4 cents, no OOB Pileup correction for MC with added signal
    cuts.AddCutPCMCalo("10130023","0dm00009ab770c00amd0404000","411790105ke30220000","0s331031000000d0"); // 00-10%
    cuts.AddCutPCMCalo("11310023","0dm00009ab770c00amd0404000","411790105ke30220000","0s331031000000d0"); // 10-30%
    cuts.AddCutPCMCalo("13530023","0dm00009ab770c00amd0404000","411790105ke30220000","0s331031000000d0"); // 30-50%
    cuts.AddCutPCMCalo("15910023","0dm00009ab770c00amd0404000","411790105ke30220000","0s331031000000d0"); // 50-90%
  } else if (trainConfig == 1053){ // 4 cents, no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutPCMCalo("10130023","0dm00009ab770c00amd0404000","411790105ke30220000","0s331031000000d0"); // 00-10%
    cuts.AddCutPCMCalo("11310023","0dm00009ab770c00amd0404000","411790105ke30220000","0s331031000000d0"); // 10-30%
    cuts.AddCutPCMCalo("13530023","0dm00009ab770c00amd0404000","411790105ke30220000","0s331031000000d0"); // 30-50%
    cuts.AddCutPCMCalo("15910023","0dm00009ab770c00amd0404000","411790105ke30220000","0s331031000000d0"); // 50-90%
  } else if (trainConfig == 1054){ // 4 cents, no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutPCMCalo("10130023","0dm00009ab770c00amd0404000","411790105ke30220000","0s331031000000d0"); // 00-10%
    cuts.AddCutPCMCalo("11310023","0dm00009ab770c00amd0404000","411790105ke30220000","0s331031000000d0"); // 10-30%
    cuts.AddCutPCMCalo("13530023","0dm00009ab770c00amd0404000","411790105ke30220000","0s331031000000d0"); // 30-50%
    cuts.AddCutPCMCalo("15910023","0dm00009ab770c00amd0404000","411790105ke30220000","0s331031000000d0"); // 50-90%
  } else if (trainConfig == 1055){ // 4 cents with NOC NonLin like data
    cuts.AddCutPCMCalo("10130e03","0dm00009ab770c00amd0404000","411790105we30220000","0s331031000000d0"); // 00-10%
    cuts.AddCutPCMCalo("11310e03","0dm00009ab770c00amd0404000","411790105we30220000","0s331031000000d0"); // 10-30%
    cuts.AddCutPCMCalo("13530e03","0dm00009ab770c00amd0404000","411790105we30220000","0s331031000000d0"); // 30-50%
    cuts.AddCutPCMCalo("15910e03","0dm00009ab770c00amd0404000","411790105we30220000","0s331031000000d0"); // 50-90%
  } else if (trainConfig == 1056){ // 4 cents, with NOC NonLin like MC no OOB Pileup correction for MC
    cuts.AddCutPCMCalo("10130053","0dm00009ab770c00amd0404000","411790105xe30220000","0s331031000000d0"); // 00-10%
    cuts.AddCutPCMCalo("11310053","0dm00009ab770c00amd0404000","411790105xe30220000","0s331031000000d0"); // 10-30%
    cuts.AddCutPCMCalo("13530053","0dm00009ab770c00amd0404000","411790105xe30220000","0s331031000000d0"); // 30-50%
    cuts.AddCutPCMCalo("15910053","0dm00009ab770c00amd0404000","411790105xe30220000","0s331031000000d0"); // 50-90%
  } else if (trainConfig == 1057){ // 4 cents, with NOC NonLin like MC no OOB Pileup correction for MC with added signal
    cuts.AddCutPCMCalo("10130023","0dm00009ab770c00amd0404000","411790105xe30220000","0s331031000000d0"); // 00-10%
    cuts.AddCutPCMCalo("11310023","0dm00009ab770c00amd0404000","411790105xe30220000","0s331031000000d0"); // 10-30%
    cuts.AddCutPCMCalo("13530023","0dm00009ab770c00amd0404000","411790105xe30220000","0s331031000000d0"); // 30-50%
    cuts.AddCutPCMCalo("15910023","0dm00009ab770c00amd0404000","411790105xe30220000","0s331031000000d0"); // 50-90%
  } else if (trainConfig == 1058){ // 4 cents, with NOC NonLin like MC no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutPCMCalo("10130023","0dm00009ab770c00amd0404000","411790105xe30220000","0s331031000000d0"); // 00-10%
    cuts.AddCutPCMCalo("11310023","0dm00009ab770c00amd0404000","411790105xe30220000","0s331031000000d0"); // 10-30%
    cuts.AddCutPCMCalo("13530023","0dm00009ab770c00amd0404000","411790105xe30220000","0s331031000000d0"); // 30-50%
    cuts.AddCutPCMCalo("15910023","0dm00009ab770c00amd0404000","411790105xe30220000","0s331031000000d0"); // 50-90%
  } else if (trainConfig == 1059){ // 4 cents, with NOC NonLin like MC no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutPCMCalo("10130023","0dm00009ab770c00amd0404000","411790105xe30220000","0s331031000000d0"); // 00-10%
    cuts.AddCutPCMCalo("11310023","0dm00009ab770c00amd0404000","411790105xe30220000","0s331031000000d0"); // 10-30%
    cuts.AddCutPCMCalo("13530023","0dm00009ab770c00amd0404000","411790105xe30220000","0s331031000000d0"); // 30-50%
    cuts.AddCutPCMCalo("15910023","0dm00009ab770c00amd0404000","411790105xe30220000","0s331031000000d0"); // 50-90%
  } else if (trainConfig == 1060){ // 4 cents with mean N matched tracks per cluster data
    cuts.AddCutPCMCalo("10130e03","0dm00009ab770c00amd0404000","411790105te30220000","0s331031000000d0"); // 00-10%
    cuts.AddCutPCMCalo("11310e03","0dm00009ab770c00amd0404000","411790105te30220000","0s331031000000d0"); // 10-30%
    cuts.AddCutPCMCalo("13530e03","0dm00009ab770c00amd0404000","411790105te30220000","0s331031000000d0"); // 30-50%
    cuts.AddCutPCMCalo("15910e03","0dm00009ab770c00amd0404000","411790105te30220000","0s331031000000d0"); // 50-90%
  } else if (trainConfig == 1061){ // 4 cents, with mean N matched tracks per cluster MC no OOB Pileup correction for MC
    cuts.AddCutPCMCalo("10130053","0dm00009ab770c00amd0404000","411790105te30220000","0s331031000000d0"); // 00-10%
    cuts.AddCutPCMCalo("11310053","0dm00009ab770c00amd0404000","411790105te30220000","0s331031000000d0"); // 10-30%
    cuts.AddCutPCMCalo("13530053","0dm00009ab770c00amd0404000","411790105te30220000","0s331031000000d0"); // 30-50%
    cuts.AddCutPCMCalo("15910053","0dm00009ab770c00amd0404000","411790105te30220000","0s331031000000d0"); // 50-90%
  } else if (trainConfig == 1062){ // 4 cents, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal
    cuts.AddCutPCMCalo("10130023","0dm00009ab770c00amd0404000","411790105te30220000","0s331031000000d0"); // 00-10%
    cuts.AddCutPCMCalo("11310023","0dm00009ab770c00amd0404000","411790105te30220000","0s331031000000d0"); // 10-30%
    cuts.AddCutPCMCalo("13530023","0dm00009ab770c00amd0404000","411790105te30220000","0s331031000000d0"); // 30-50%
    cuts.AddCutPCMCalo("15910023","0dm00009ab770c00amd0404000","411790105te30220000","0s331031000000d0"); // 50-90%
  } else if (trainConfig == 1063){ // 4 cents, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutPCMCalo("10130023","0dm00009ab770c00amd0404000","411790105te30220000","0s331031000000d0"); // 00-10%
    cuts.AddCutPCMCalo("11310023","0dm00009ab770c00amd0404000","411790105te30220000","0s331031000000d0"); // 10-30%
    cuts.AddCutPCMCalo("13530023","0dm00009ab770c00amd0404000","411790105te30220000","0s331031000000d0"); // 30-50%
    cuts.AddCutPCMCalo("15910023","0dm00009ab770c00amd0404000","411790105te30220000","0s331031000000d0"); // 50-90%
  } else if (trainConfig == 1064){ // 4 cents, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutPCMCalo("10130023","0dm00009ab770c00amd0404000","411790105te30220000","0s331031000000d0"); // 00-10%
    cuts.AddCutPCMCalo("11310023","0dm00009ab770c00amd0404000","411790105te30220000","0s331031000000d0"); // 10-30%
    cuts.AddCutPCMCalo("13530023","0dm00009ab770c00amd0404000","411790105te30220000","0s331031000000d0"); // 30-50%
    cuts.AddCutPCMCalo("15910023","0dm00009ab770c00amd0404000","411790105te30220000","0s331031000000d0"); // 50-90%
  } else if (trainConfig == 1065){ // 4 cents with mean N matched tracks per cluster data
    cuts.AddCutPCMCalo("10130e03","0dm00009ab770c00amd0404000","411790105ye30220000","0s33103100000010"); // 00-10%
    cuts.AddCutPCMCalo("11310e03","0dm00009ab770c00amd0404000","411790105ye30220000","0s33103100000010"); // 10-30%
    cuts.AddCutPCMCalo("13530e03","0dm00009ab770c00amd0404000","411790105ye30220000","0s33103100000010"); // 30-50%
    cuts.AddCutPCMCalo("15910e03","0dm00009ab770c00amd0404000","411790105ye30220000","0s33103100000010"); // 50-90%
  } else if (trainConfig == 1066){ // 4 cents, with mean N matched tracks per cluster MC no OOB Pileup correction for MC
    cuts.AddCutPCMCalo("10130053","0dm00009ab770c00amd0404000","411790105ye30220000","0s33103100000010"); // 00-10%
    cuts.AddCutPCMCalo("11310053","0dm00009ab770c00amd0404000","411790105ye30220000","0s33103100000010"); // 10-30%
    cuts.AddCutPCMCalo("13530053","0dm00009ab770c00amd0404000","411790105ye30220000","0s33103100000010"); // 30-50%
    cuts.AddCutPCMCalo("15910053","0dm00009ab770c00amd0404000","411790105ye30220000","0s33103100000010"); // 50-90%
  // **********************************************************************************************************
  // **************************** PCM-EMC configurations PbPb run 2 2018 pass 3 cent ******************************
  // **********************************************************************************************************
  // Standard: EMCal + DCal, only 2ndary cluster track matching, mixed event
  } else if (trainConfig == 1100){ // central
    cuts.AddCutPCMCalo("10130e03","0dm00009ab770c00amd0404000","411790105ke30220000","0133103100000010"); // 00-10%
  } else if (trainConfig == 1101){ // central, no OOB Pileup correction for MC
    cuts.AddCutPCMCalo("10130053","0dm00009ab770c00amd0404000","411790105ke30220000","0133103100000010"); // 00-10%
  } else if (trainConfig == 1102){ // central, no OOB Pileup correction for MC with added signal
    cuts.AddCutPCMCalo("10130023","0dm00009ab770c00amd0404000","411790105ke30220000","0133103100000010"); // 00-10%
  } else if (trainConfig == 1103){ // central, no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutPCMCalo("10130023","0dm00009ab770c00amd0404000","411790105ke30220000","0133103100000010"); // 00-10%
  } else if (trainConfig == 1104){ // central, no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutPCMCalo("10130023","0dm00009ab770c00amd0404000","411790105ke30220000","0133103100000010"); // 00-10%
  } else if (trainConfig == 1105){ // central with NOC NonLin like data
    cuts.AddCutPCMCalo("10130e03","0dm00009ab770c00amd0404000","411790105we30220000","0133103100000010"); // 00-10%
  } else if (trainConfig == 1106){ // central, with NOC NonLin like MC no OOB Pileup correction for MC
    cuts.AddCutPCMCalo("10130053","0dm00009ab770c00amd0404000","411790105xe30220000","0133103100000010"); // 00-10%
  } else if (trainConfig == 1107){ // central, with NOC NonLin like MC no OOB Pileup correction for MC with added signal
    cuts.AddCutPCMCalo("10130023","0dm00009ab770c00amd0404000","411790105xe30220000","0133103100000010"); // 00-10%
  } else if (trainConfig == 1108){ // central, with NOC NonLin like MC no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutPCMCalo("10130023","0dm00009ab770c00amd0404000","411790105xe30220000","0133103100000010"); // 00-10%
  } else if (trainConfig == 1109){ // central, with NOC NonLin like MC no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutPCMCalo("10130023","0dm00009ab770c00amd0404000","411790105xe30220000","0133103100000010"); // 00-10%
  } else if (trainConfig == 1110){ // central with mean N matched tracks per cluster data
    cuts.AddCutPCMCalo("10130e03","0dm00009ab770c00amd0404000","411790105te30220000","0133103100000010"); // 00-10%
  } else if (trainConfig == 1111){ // central, with mean N matched tracks per cluster MC no OOB Pileup correction for MC
    cuts.AddCutPCMCalo("10130053","0dm00009ab770c00amd0404000","411790105te30220000","0133103100000010"); // 00-10%
  } else if (trainConfig == 1112){ // central, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal
    cuts.AddCutPCMCalo("10130023","0dm00009ab770c00amd0404000","411790105te30220000","0133103100000010"); // 00-10%
  } else if (trainConfig == 1113){ // central, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutPCMCalo("10130023","0dm00009ab770c00amd0404000","411790105te30220000","0133103100000010"); // 00-10%
  } else if (trainConfig == 1114){ // central, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutPCMCalo("10130023","0dm00009ab770c00amd0404000","411790105te30220000","0133103100000010"); // 00-10%
  } else if (trainConfig == 1115){ // central with mean N matched tracks per cluster data
    cuts.AddCutPCMCalo("10130e03","0dm00009ab770c00amd0404000","411790105ye30220000","0133103100000010"); // 00-10%
  } else if (trainConfig == 1116){ // central, with mean N matched tracks per cluster MC no OOB Pileup correction for MC
    cuts.AddCutPCMCalo("10130053","0dm00009ab770c00amd0404000","411790105ye30220000","0133103100000010"); // 00-10%
  // with rotation background
  } else if (trainConfig == 1150){ // central
    cuts.AddCutPCMCalo("10130e03","0dm00009ab770c00amd0404000","411790105ke30220000","0s331031000000d0"); // 00-10%
  } else if (trainConfig == 1151){ // central, no OOB Pileup correction for MC
    cuts.AddCutPCMCalo("10130053","0dm00009ab770c00amd0404000","411790105ke30220000","0s331031000000d0"); // 00-10%
  } else if (trainConfig == 1152){ // central, no OOB Pileup correction for MC with added signal
    cuts.AddCutPCMCalo("10130023","0dm00009ab770c00amd0404000","411790105ke30220000","0s331031000000d0"); // 00-10%
  } else if (trainConfig == 1153){ // central, no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutPCMCalo("10130023","0dm00009ab770c00amd0404000","411790105ke30220000","0s331031000000d0"); // 00-10%
  } else if (trainConfig == 1154){ // central, no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutPCMCalo("10130023","0dm00009ab770c00amd0404000","411790105ke30220000","0s331031000000d0"); // 00-10%
  } else if (trainConfig == 1155){ // central with NOC NonLin like data
    cuts.AddCutPCMCalo("10130e03","0dm00009ab770c00amd0404000","411790105we30220000","0s331031000000d0"); // 00-10%
  } else if (trainConfig == 1156){ // central, with NOC NonLin like MC no OOB Pileup correction for MC
    cuts.AddCutPCMCalo("10130053","0dm00009ab770c00amd0404000","411790105xe30220000","0s331031000000d0"); // 00-10%
  } else if (trainConfig == 1157){ // central, with NOC NonLin like MC no OOB Pileup correction for MC with added signal
    cuts.AddCutPCMCalo("10130023","0dm00009ab770c00amd0404000","411790105xe30220000","0s331031000000d0"); // 00-10%
  } else if (trainConfig == 1158){ // central, with NOC NonLin like MC no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutPCMCalo("10130023","0dm00009ab770c00amd0404000","411790105xe30220000","0s331031000000d0"); // 00-10%
  } else if (trainConfig == 1159){ // central, with NOC NonLin like MC no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutPCMCalo("10130023","0dm00009ab770c00amd0404000","411790105xe30220000","0s331031000000d0"); // 00-10%
  } else if (trainConfig == 1160){ // central with mean N matched tracks per cluster data
    cuts.AddCutPCMCalo("10130e03","0dm00009ab770c00amd0404000","411790105te30220000","0s331031000000d0"); // 00-10%
  } else if (trainConfig == 1161){ // central, with mean N matched tracks per cluster MC no OOB Pileup correction for MC
    cuts.AddCutPCMCalo("10130053","0dm00009ab770c00amd0404000","411790105te30220000","0s331031000000d0"); // 00-10%
  } else if (trainConfig == 1162){ // central, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal
    cuts.AddCutPCMCalo("10130023","0dm00009ab770c00amd0404000","411790105te30220000","0s331031000000d0"); // 00-10%
  } else if (trainConfig == 1163){ // central, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutPCMCalo("10130023","0dm00009ab770c00amd0404000","411790105te30220000","0s331031000000d0"); // 00-10%
  } else if (trainConfig == 1164){ // central, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutPCMCalo("10130023","0dm00009ab770c00amd0404000","411790105te30220000","0s331031000000d0"); // 00-10%
  } else if (trainConfig == 1165){ // central with mean N matched tracks per cluster data
    cuts.AddCutPCMCalo("10130e03","0dm00009ab770c00amd0404000","411790105ye30220000","0s33103100000010"); // 00-10%
  } else if (trainConfig == 1166){ // central, with mean N matched tracks per cluster MC no OOB Pileup correction for MC
    cuts.AddCutPCMCalo("10130053","0dm00009ab770c00amd0404000","411790105ye30220000","0s33103100000010"); // 00-10%
  // **********************************************************************************************************
  // ************************** PCM-EMC configurations PbPb run 2 2018 pass 3 semicent ****************************
  // **********************************************************************************************************
  // Standard: EMCal + DCal, only 2ndary cluster track matching, mixed event
  } else if (trainConfig == 1300){ // semicent
    cuts.AddCutPCMCalo("11310e03","0dm00009ab770c00amd0404000","411790105ke30220000","0133103100000010"); // 10-30%
  } else if (trainConfig == 1301){ // semicent, no OOB Pileup correction for MC
    cuts.AddCutPCMCalo("11310053","0dm00009ab770c00amd0404000","411790105ke30220000","0133103100000010"); // 10-30%
  } else if (trainConfig == 1302){ // semicent, no OOB Pileup correction for MC with added signal
    cuts.AddCutPCMCalo("11310023","0dm00009ab770c00amd0404000","411790105ke30220000","0133103100000010"); // 10-30%
  } else if (trainConfig == 1303){ // semicent, no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutPCMCalo("11310023","0dm00009ab770c00amd0404000","411790105ke30220000","0133103100000010"); // 10-30%
  } else if (trainConfig == 1304){ // semicent, no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutPCMCalo("11310023","0dm00009ab770c00amd0404000","411790105ke30220000","0133103100000010"); // 10-30%
  } else if (trainConfig == 1305){ // semicent with NOC NonLin like data
    cuts.AddCutPCMCalo("11310e03","0dm00009ab770c00amd0404000","411790105we30220000","0133103100000010"); // 10-30%
  } else if (trainConfig == 1306){ // semicent, with NOC NonLin like MC no OOB Pileup correction for MC
    cuts.AddCutPCMCalo("11310053","0dm00009ab770c00amd0404000","411790105xe30220000","0133103100000010"); // 10-30%
  } else if (trainConfig == 1307){ // semicent, with NOC NonLin like MC no OOB Pileup correction for MC with added signal
    cuts.AddCutPCMCalo("11310023","0dm00009ab770c00amd0404000","411790105xe30220000","0133103100000010"); // 10-30%
  } else if (trainConfig == 1308){ // semicent, with NOC NonLin like MC no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutPCMCalo("11310023","0dm00009ab770c00amd0404000","411790105xe30220000","0133103100000010"); // 10-30%
  } else if (trainConfig == 1309){ // semicent, with NOC NonLin like MC no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutPCMCalo("11310023","0dm00009ab770c00amd0404000","411790105xe30220000","0133103100000010"); // 10-30%
  } else if (trainConfig == 1310){ // semicent with mean N matched tracks per cluster data
    cuts.AddCutPCMCalo("11310e03","0dm00009ab770c00amd0404000","411790105te30220000","0133103100000010"); // 10-30%
  } else if (trainConfig == 1311){ // semicent, with mean N matched tracks per cluster MC no OOB Pileup correction for MC
    cuts.AddCutPCMCalo("11310053","0dm00009ab770c00amd0404000","411790105te30220000","0133103100000010"); // 10-30%
  } else if (trainConfig == 1312){ // semicent, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal
    cuts.AddCutPCMCalo("11310023","0dm00009ab770c00amd0404000","411790105te30220000","0133103100000010"); // 10-30%
  } else if (trainConfig == 1313){ // semicent, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutPCMCalo("11310023","0dm00009ab770c00amd0404000","411790105te30220000","0133103100000010"); // 10-30%
  } else if (trainConfig == 1314){ // semicent, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutPCMCalo("11310023","0dm00009ab770c00amd0404000","411790105te30220000","0133103100000010"); // 10-30%
  } else if (trainConfig == 1315){ // semicent with mean N matched tracks per cluster data
    cuts.AddCutPCMCalo("11310e03","0dm00009ab770c00amd0404000","411790105ye30220000","0133103100000010"); // 10-30%
  } else if (trainConfig == 1316){ // semicent, with mean N matched tracks per cluster MC no OOB Pileup correction for MC
    cuts.AddCutPCMCalo("11310053","0dm00009ab770c00amd0404000","411790105ye30220000","0133103100000010"); // 10-30%
  // with rotation background
  } else if (trainConfig == 1350){ // semicent
    cuts.AddCutPCMCalo("11310e03","0dm00009ab770c00amd0404000","411790105ke30220000","0s33103100000010"); // 10-30%
  } else if (trainConfig == 1351){ // semicent, no OOB Pileup correction for MC
    cuts.AddCutPCMCalo("11310053","0dm00009ab770c00amd0404000","411790105ke30220000","0s33103100000010"); // 10-30%
  } else if (trainConfig == 1352){ // semicent, no OOB Pileup correction for MC with added signal
    cuts.AddCutPCMCalo("11310023","0dm00009ab770c00amd0404000","411790105ke30220000","0s33103100000010"); // 10-30%
  } else if (trainConfig == 1353){ // semicent, no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutPCMCalo("11310023","0dm00009ab770c00amd0404000","411790105ke30220000","0s33103100000010"); // 10-30%
  } else if (trainConfig == 1354){ // semicent, no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutPCMCalo("11310023","0dm00009ab770c00amd0404000","411790105ke30220000","0s33103100000010"); // 10-30%
  } else if (trainConfig == 1355){ // semicent with NOC NonLin like data
    cuts.AddCutPCMCalo("11310e03","0dm00009ab770c00amd0404000","411790105we30220000","0s33103100000010"); // 10-30%
  } else if (trainConfig == 1356){ // semicent, with NOC NonLin like MC no OOB Pileup correction for MC
    cuts.AddCutPCMCalo("11310053","0dm00009ab770c00amd0404000","411790105xe30220000","0s33103100000010"); // 10-30%
  } else if (trainConfig == 1357){ // semicent, with NOC NonLin like MC no OOB Pileup correction for MC with added signal
    cuts.AddCutPCMCalo("11310023","0dm00009ab770c00amd0404000","411790105xe30220000","0s33103100000010"); // 10-30%
  } else if (trainConfig == 1358){ // semicent, with NOC NonLin like MC no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutPCMCalo("11310023","0dm00009ab770c00amd0404000","411790105xe30220000","0s33103100000010"); // 10-30%
  } else if (trainConfig == 1359){ // semicent, with NOC NonLin like MC no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutPCMCalo("11310023","0dm00009ab770c00amd0404000","411790105xe30220000","0s33103100000010"); // 10-30%
  } else if (trainConfig == 1360){ // semicent with mean N matched tracks per cluster data
    cuts.AddCutPCMCalo("11310e03","0dm00009ab770c00amd0404000","411790105te30220000","0s33103100000010"); // 10-30%
  } else if (trainConfig == 1361){ // semicent, with mean N matched tracks per cluster MC no OOB Pileup correction for MC
    cuts.AddCutPCMCalo("11310053","0dm00009ab770c00amd0404000","411790105te30220000","0s33103100000010"); // 10-30%
  } else if (trainConfig == 1362){ // semicent, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal
    cuts.AddCutPCMCalo("11310023","0dm00009ab770c00amd0404000","411790105te30220000","0s33103100000010"); // 10-30%
  } else if (trainConfig == 1363){ // semicent, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutPCMCalo("11310023","0dm00009ab770c00amd0404000","411790105te30220000","0s33103100000010"); // 10-30%
  } else if (trainConfig == 1364){ // semicent, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutPCMCalo("11310023","0dm00009ab770c00amd0404000","411790105te30220000","0s33103100000010"); // 10-30%
  } else if (trainConfig == 1365){ // semicent with mean N matched tracks per cluster data
    cuts.AddCutPCMCalo("11310e03","0dm00009ab770c00amd0404000","411790105ye30220000","0s33103100000010"); // 10-30%
  } else if (trainConfig == 1366){ // semicent, with mean N matched tracks per cluster MC no OOB Pileup correction for MC
    cuts.AddCutPCMCalo("11310053","0dm00009ab770c00amd0404000","411790105ye30220000","0s33103100000010"); // 10-30%
  // **********************************************************************************************************
  // *********************** PCM-EMC configurations PbPb run 2 2018 pass 3 semiperipheral *************************
  // **********************************************************************************************************
  // Standard: EMCal + DCal, only 2ndary cluster track matching, mixed event
  } else if (trainConfig == 1500){ // semiperipheral
    cuts.AddCutPCMCalo("13530e03","0dm00009ab770c00amd0404000","411790105ke30220000","0133103100000010"); // 30-50%
  } else if (trainConfig == 1501){ // semiperipheral, no OOB Pileup correction for MC
    cuts.AddCutPCMCalo("13530053","0dm00009ab770c00amd0404000","411790105ke30220000","0133103100000010"); // 30-50%
  } else if (trainConfig == 1502){ // semiperipheral, no OOB Pileup correction for MC with added signal
    cuts.AddCutPCMCalo("13530023","0dm00009ab770c00amd0404000","411790105ke30220000","0133103100000010"); // 30-50%
  } else if (trainConfig == 1503){ // semiperipheral, no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutPCMCalo("13530023","0dm00009ab770c00amd0404000","411790105ke30220000","0133103100000010"); // 30-50%
  } else if (trainConfig == 1504){ // semiperipheral, no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutPCMCalo("13530023","0dm00009ab770c00amd0404000","411790105ke30220000","0133103100000010"); // 30-50%
  } else if (trainConfig == 1505){ // semiperipheral with NOC NonLin like data
    cuts.AddCutPCMCalo("13530e03","0dm00009ab770c00amd0404000","411790105we30220000","0133103100000010"); // 30-50%
  } else if (trainConfig == 1506){ // semiperipheral, with NOC NonLin like MC no OOB Pileup correction for MC
    cuts.AddCutPCMCalo("13530053","0dm00009ab770c00amd0404000","411790105xe30220000","0133103100000010"); // 30-50%
  } else if (trainConfig == 1507){ // semiperipheral, with NOC NonLin like MC no OOB Pileup correction for MC with added signal
    cuts.AddCutPCMCalo("13530023","0dm00009ab770c00amd0404000","411790105xe30220000","0133103100000010"); // 30-50%
  } else if (trainConfig == 1508){ // semiperipheral, with NOC NonLin like MC no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutPCMCalo("13530023","0dm00009ab770c00amd0404000","411790105xe30220000","0133103100000010"); // 30-50%
  } else if (trainConfig == 1509){ // semiperipheral, with NOC NonLin like MC no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutPCMCalo("13530023","0dm00009ab770c00amd0404000","411790105xe30220000","0133103100000010"); // 30-50%
  } else if (trainConfig == 1510){ // semiperipheral with mean N matched tracks per cluster data
    cuts.AddCutPCMCalo("13530e03","0dm00009ab770c00amd0404000","411790105te30220000","0133103100000010"); // 30-50%
  } else if (trainConfig == 1511){ // semiperipheral, with mean N matched tracks per cluster MC no OOB Pileup correction for MC
    cuts.AddCutPCMCalo("13530053","0dm00009ab770c00amd0404000","411790105te30220000","0133103100000010"); // 30-50%
  } else if (trainConfig == 1512){ // semiperipheral, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal
    cuts.AddCutPCMCalo("13530023","0dm00009ab770c00amd0404000","411790105te30220000","0133103100000010"); // 30-50%
  } else if (trainConfig == 1513){ // semiperipheral, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutPCMCalo("13530023","0dm00009ab770c00amd0404000","411790105te30220000","0133103100000010"); // 30-50%
  } else if (trainConfig == 1514){ // semiperipheral, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutPCMCalo("13530023","0dm00009ab770c00amd0404000","411790105te30220000","0133103100000010"); // 30-50%
  } else if (trainConfig == 1515){ // semiperipheral with mean N matched tracks per cluster data
    cuts.AddCutPCMCalo("13530e03","0dm00009ab770c00amd0404000","411790105ye30220000","0133103100000010"); // 30-50%
  } else if (trainConfig == 1516){ // semiperipheral, with mean N matched tracks per cluster MC no OOB Pileup correction for MC
    cuts.AddCutPCMCalo("13530053","0dm00009ab770c00amd0404000","411790105ye30220000","0133103100000010"); // 30-50%
  // with rotation background
  } else if (trainConfig == 1550){ // semiperipheral
    cuts.AddCutPCMCalo("13530e03","0dm00009ab770c00amd0404000","411790105ke30220000","0s33103100000010"); // 30-50%
  } else if (trainConfig == 1551){ // semiperipheral, no OOB Pileup correction for MC
    cuts.AddCutPCMCalo("13530053","0dm00009ab770c00amd0404000","411790105ke30220000","0s33103100000010"); // 30-50%
  } else if (trainConfig == 1552){ // semiperipheral, no OOB Pileup correction for MC with added signal
    cuts.AddCutPCMCalo("13530023","0dm00009ab770c00amd0404000","411790105ke30220000","0s33103100000010"); // 30-50%
  } else if (trainConfig == 1553){ // semiperipheral, no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutPCMCalo("13530023","0dm00009ab770c00amd0404000","411790105ke30220000","0s33103100000010"); // 30-50%
  } else if (trainConfig == 1554){ // semiperipheral, no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutPCMCalo("13530023","0dm00009ab770c00amd0404000","411790105ke30220000","0s33103100000010"); // 30-50%
  } else if (trainConfig == 1555){ // semiperipheral with NOC NonLin like data
    cuts.AddCutPCMCalo("13530e03","0dm00009ab770c00amd0404000","411790105we30220000","0s33103100000010"); // 30-50%
  } else if (trainConfig == 1556){ // semiperipheral, with NOC NonLin like MC no OOB Pileup correction for MC
    cuts.AddCutPCMCalo("13530053","0dm00009ab770c00amd0404000","411790105xe30220000","0s33103100000010"); // 30-50%
  } else if (trainConfig == 1557){ // semiperipheral, with NOC NonLin like MC no OOB Pileup correction for MC with added signal
    cuts.AddCutPCMCalo("13530023","0dm00009ab770c00amd0404000","411790105xe30220000","0s33103100000010"); // 30-50%
  } else if (trainConfig == 1558){ // semiperipheral, with NOC NonLin like MC no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutPCMCalo("13530023","0dm00009ab770c00amd0404000","411790105xe30220000","0s33103100000010"); // 30-50%
  } else if (trainConfig == 1559){ // semiperipheral, with NOC NonLin like MC no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutPCMCalo("13530023","0dm00009ab770c00amd0404000","411790105xe30220000","0s33103100000010"); // 30-50%
  } else if (trainConfig == 1560){ // semiperipheral with mean N matched tracks per cluster data
    cuts.AddCutPCMCalo("13530e03","0dm00009ab770c00amd0404000","411790105te30220000","0s33103100000010"); // 30-50%
  } else if (trainConfig == 1561){ // semiperipheral, with mean N matched tracks per cluster MC no OOB Pileup correction for MC
    cuts.AddCutPCMCalo("13530053","0dm00009ab770c00amd0404000","411790105te30220000","0s33103100000010"); // 30-50%
  } else if (trainConfig == 1562){ // semiperipheral, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal
    cuts.AddCutPCMCalo("13530023","0dm00009ab770c00amd0404000","411790105te30220000","0s33103100000010"); // 30-50%
  } else if (trainConfig == 1563){ // semiperipheral, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutPCMCalo("13530023","0dm00009ab770c00amd0404000","411790105te30220000","0s33103100000010"); // 30-50%
  } else if (trainConfig == 1564){ // semiperipheral, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutPCMCalo("13530023","0dm00009ab770c00amd0404000","411790105te30220000","0s33103100000010"); // 30-50%
  } else if (trainConfig == 1565){ // semiperipheral with mean N matched tracks per cluster data
    cuts.AddCutPCMCalo("13530e03","0dm00009ab770c00amd0404000","411790105ye30220000","0s33103100000010"); // 30-50%
  } else if (trainConfig == 1566){ // semiperipheral, with mean N matched tracks per cluster MC no OOB Pileup correction for MC
    cuts.AddCutPCMCalo("13530053","0dm00009ab770c00amd0404000","411790105ye30220000","0s33103100000010"); // 30-50%
  // **********************************************************************************************************
  // ************************* PCM-EMC configurations PbPb run 2 2018 pass 3 peripheral ***************************
  // **********************************************************************************************************
  // Standard: EMCal + DCal, only 2ndary cluster track matching, mixed event
  } else if (trainConfig == 1900){ // peripheral
    cuts.AddCutPCMCalo("15910e03","0dm00009ab770c00amd0404000","411790105ke30220000","0133103100000010"); // 50-90%
  } else if (trainConfig == 1901){ // peripheral, no OOB Pileup correction for MC
    cuts.AddCutPCMCalo("15910053","0dm00009ab770c00amd0404000","411790105ke30220000","0133103100000010"); // 50-90%
  } else if (trainConfig == 1902){ // peripheral, no OOB Pileup correction for MC with added signal
    cuts.AddCutPCMCalo("15910023","0dm00009ab770c00amd0404000","411790105ke30220000","0133103100000010"); // 50-90%
  } else if (trainConfig == 1903){ // peripheral, no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutPCMCalo("15910023","0dm00009ab770c00amd0404000","411790105ke30220000","0133103100000010"); // 50-90%
  } else if (trainConfig == 1904){ // peripheral, no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutPCMCalo("15910023","0dm00009ab770c00amd0404000","411790105ke30220000","0133103100000010"); // 50-90%
  } else if (trainConfig == 1905){ // peripheral with NOC NonLin like data
    cuts.AddCutPCMCalo("15910e03","0dm00009ab770c00amd0404000","411790105we30220000","0133103100000010"); // 50-90%
  } else if (trainConfig == 1906){ // peripheral, with NOC NonLin like MC no OOB Pileup correction for MC
    cuts.AddCutPCMCalo("15910053","0dm00009ab770c00amd0404000","411790105xe30220000","0133103100000010"); // 50-90%
  } else if (trainConfig == 1907){ // peripheral, with NOC NonLin like MC no OOB Pileup correction for MC with added signal
    cuts.AddCutPCMCalo("15910023","0dm00009ab770c00amd0404000","411790105xe30220000","0133103100000010"); // 50-90%
  } else if (trainConfig == 1908){ // peripheral, with NOC NonLin like MC no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutPCMCalo("15910023","0dm00009ab770c00amd0404000","411790105xe30220000","0133103100000010"); // 50-90%
  } else if (trainConfig == 1909){ // peripheral, with NOC NonLin like MC no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutPCMCalo("15910023","0dm00009ab770c00amd0404000","411790105xe30220000","0133103100000010"); // 50-90%
  } else if (trainConfig == 1910){ // peripheral with mean N matched tracks per cluster data
    cuts.AddCutPCMCalo("15910e03","0dm00009ab770c00amd0404000","411790105te30220000","0133103100000010"); // 50-90%
  } else if (trainConfig == 1911){ // peripheral, with mean N matched tracks per cluster MC no OOB Pileup correction for MC
    cuts.AddCutPCMCalo("15910053","0dm00009ab770c00amd0404000","411790105te30220000","0133103100000010"); // 50-90%
  } else if (trainConfig == 1912){ // peripheral, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal
    cuts.AddCutPCMCalo("15910023","0dm00009ab770c00amd0404000","411790105te30220000","0133103100000010"); // 50-90%
  } else if (trainConfig == 1913){ // peripheral, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutPCMCalo("15910023","0dm00009ab770c00amd0404000","411790105te30220000","0133103100000010"); // 50-90%
  } else if (trainConfig == 1914){ // peripheral, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutPCMCalo("15910023","0dm00009ab770c00amd0404000","411790105te30220000","0133103100000010"); // 50-90%
  } else if (trainConfig == 1915){ // peripheral with mean N matched tracks per cluster data
    cuts.AddCutPCMCalo("15910e03","0dm00009ab770c00amd0404000","411790105ye30220000","0133103100000010"); // 50-90%
  } else if (trainConfig == 1916){ // peripheral, with mean N matched tracks per cluster MC no OOB Pileup correction for MC
    cuts.AddCutPCMCalo("15910053","0dm00009ab770c00amd0404000","411790105ye30220000","0133103100000010"); // 50-90%
  // with rotation background
  } else if (trainConfig == 1950){ // peripheral
    cuts.AddCutPCMCalo("15910e03","0dm00009ab770c00amd0404000","411790105ke30220000","0s33103100000010"); // 50-90%
  } else if (trainConfig == 1951){ // peripheral, no OOB Pileup correction for MC
    cuts.AddCutPCMCalo("15910053","0dm00009ab770c00amd0404000","411790105ke30220000","0s33103100000010"); // 50-90%
  } else if (trainConfig == 1952){ // peripheral, no OOB Pileup correction for MC with added signal
    cuts.AddCutPCMCalo("15910023","0dm00009ab770c00amd0404000","411790105ke30220000","0s33103100000010"); // 50-90%
  } else if (trainConfig == 1953){ // peripheral, no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutPCMCalo("15910023","0dm00009ab770c00amd0404000","411790105ke30220000","0s33103100000010"); // 50-90%
  } else if (trainConfig == 1954){ // peripheral, no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutPCMCalo("15910023","0dm00009ab770c00amd0404000","411790105ke30220000","0s33103100000010"); // 50-90%
  } else if (trainConfig == 1955){ // peripheral with NOC NonLin like data
    cuts.AddCutPCMCalo("15910e03","0dm00009ab770c00amd0404000","411790105we30220000","0s33103100000010"); // 50-90%
  } else if (trainConfig == 1956){ // peripheral, with NOC NonLin like MC no OOB Pileup correction for MC
    cuts.AddCutPCMCalo("15910053","0dm00009ab770c00amd0404000","411790105xe30220000","0s33103100000010"); // 50-90%
  } else if (trainConfig == 1957){ // peripheral, with NOC NonLin like MC no OOB Pileup correction for MC with added signal
    cuts.AddCutPCMCalo("15910023","0dm00009ab770c00amd0404000","411790105xe30220000","0s33103100000010"); // 50-90%
  } else if (trainConfig == 1958){ // peripheral, with NOC NonLin like MC no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutPCMCalo("15910023","0dm00009ab770c00amd0404000","411790105xe30220000","0s33103100000010"); // 50-90%
  } else if (trainConfig == 1959){ // peripheral, with NOC NonLin like MC no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutPCMCalo("15910023","0dm00009ab770c00amd0404000","411790105xe30220000","0s33103100000010"); // 50-90%
  } else if (trainConfig == 1960){ // peripheral with mean N matched tracks per cluster data
    cuts.AddCutPCMCalo("15910e03","0dm00009ab770c00amd0404000","411790105te30220000","0s33103100000010"); // 50-90%
  } else if (trainConfig == 1961){ // peripheral, with mean N matched tracks per cluster MC no OOB Pileup correction for MC
    cuts.AddCutPCMCalo("15910053","0dm00009ab770c00amd0404000","411790105te30220000","0s33103100000010"); // 50-90%
  } else if (trainConfig == 1962){ // peripheral, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal
    cuts.AddCutPCMCalo("15910023","0dm00009ab770c00amd0404000","411790105te30220000","0s33103100000010"); // 50-90%
  } else if (trainConfig == 1963){ // peripheral, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutPCMCalo("15910023","0dm00009ab770c00amd0404000","411790105te30220000","0s33103100000010"); // 50-90%
  } else if (trainConfig == 1964){ // peripheral, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutPCMCalo("15910023","0dm00009ab770c00amd0404000","411790105te30220000","0s33103100000010"); // 50-90%
  } else if (trainConfig == 1965){ // peripheral with mean N matched tracks per cluster data
    cuts.AddCutPCMCalo("15910e03","0dm00009ab770c00amd0404000","411790105ye30220000","0s33103100000010"); // 50-90%
  } else if (trainConfig == 1966){ // peripheral, with mean N matched tracks per cluster MC no OOB Pileup correction for MC
    cuts.AddCutPCMCalo("15910053","0dm00009ab770c00amd0404000","411790105ye30220000","0s33103100000010"); // 50-90%

  // **********************************************************************************************************
  // ***************************** PCM-EMC configurations PbPb run 2 2015 pass 3 *******************************
  // **********************************************************************************************************
  // Standard: EMCal + DCal, only 2ndary cluster track matching, mixed event
  } else if (trainConfig == 2000){ // 4 cents
    cuts.AddCutPCMCalo("10110e03","0dm00009ab770c00amd0404000","411790105ke30220000","0133103100000010"); // 00-10%
    cuts.AddCutPCMCalo("11310e03","0dm00009ab770c00amd0404000","411790105ke30220000","0133103100000010"); // 10-30%
    cuts.AddCutPCMCalo("13510e03","0dm00009ab770c00amd0404000","411790105ke30220000","0133103100000010"); // 30-50%
    cuts.AddCutPCMCalo("15910e03","0dm00009ab770c00amd0404000","411790105ke30220000","0133103100000010"); // 50-90%
  } else if (trainConfig == 2001){ // 4 cents, no OOB Pileup correction for MC
    cuts.AddCutPCMCalo("10110053","0dm00009ab770c00amd0404000","411790105ke30220000","0133103100000010"); // 00-10%
    cuts.AddCutPCMCalo("11310053","0dm00009ab770c00amd0404000","411790105ke30220000","0133103100000010"); // 10-30%
    cuts.AddCutPCMCalo("13510053","0dm00009ab770c00amd0404000","411790105ke30220000","0133103100000010"); // 30-50%
    cuts.AddCutPCMCalo("15910053","0dm00009ab770c00amd0404000","411790105ke30220000","0133103100000010"); // 50-90%
  } else if (trainConfig == 2002){ // 4 cents, no OOB Pileup correction for MC with added signal
    cuts.AddCutPCMCalo("10110023","0dm00009ab770c00amd0404000","411790105ke30220000","0133103100000010"); // 00-10%
    cuts.AddCutPCMCalo("11310023","0dm00009ab770c00amd0404000","411790105ke30220000","0133103100000010"); // 10-30%
    cuts.AddCutPCMCalo("13510023","0dm00009ab770c00amd0404000","411790105ke30220000","0133103100000010"); // 30-50%
    cuts.AddCutPCMCalo("15910023","0dm00009ab770c00amd0404000","411790105ke30220000","0133103100000010"); // 50-90%
  } else if (trainConfig == 2003){ // 4 cents, no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutPCMCalo("10110023","0dm00009ab770c00amd0404000","411790105ke30220000","0133103100000010"); // 00-10%
    cuts.AddCutPCMCalo("11310023","0dm00009ab770c00amd0404000","411790105ke30220000","0133103100000010"); // 10-30%
    cuts.AddCutPCMCalo("13510023","0dm00009ab770c00amd0404000","411790105ke30220000","0133103100000010"); // 30-50%
    cuts.AddCutPCMCalo("15910023","0dm00009ab770c00amd0404000","411790105ke30220000","0133103100000010"); // 50-90%
  } else if (trainConfig == 2004){ // 4 cents, no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutPCMCalo("10110023","0dm00009ab770c00amd0404000","411790105ke30220000","0133103100000010"); // 00-10%
    cuts.AddCutPCMCalo("11310023","0dm00009ab770c00amd0404000","411790105ke30220000","0133103100000010"); // 10-30%
    cuts.AddCutPCMCalo("13510023","0dm00009ab770c00amd0404000","411790105ke30220000","0133103100000010"); // 30-50%
    cuts.AddCutPCMCalo("15910023","0dm00009ab770c00amd0404000","411790105ke30220000","0133103100000010"); // 50-90%
  } else if (trainConfig == 2005){ // 4 cents with NOC NonLin like data
    cuts.AddCutPCMCalo("10110e03","0dm00009ab770c00amd0404000","411790105we30220000","0133103100000010"); // 00-10%
    cuts.AddCutPCMCalo("11310e03","0dm00009ab770c00amd0404000","411790105we30220000","0133103100000010"); // 10-30%
    cuts.AddCutPCMCalo("13510e03","0dm00009ab770c00amd0404000","411790105we30220000","0133103100000010"); // 30-50%
    cuts.AddCutPCMCalo("15910e03","0dm00009ab770c00amd0404000","411790105we30220000","0133103100000010"); // 50-90%
  } else if (trainConfig == 2006){ // 4 cents, with NOC NonLin like MC no OOB Pileup correction for MC
    cuts.AddCutPCMCalo("10110053","0dm00009ab770c00amd0404000","411790105xe30220000","0133103100000010"); // 00-10%
    cuts.AddCutPCMCalo("11310053","0dm00009ab770c00amd0404000","411790105xe30220000","0133103100000010"); // 10-30%
    cuts.AddCutPCMCalo("13510053","0dm00009ab770c00amd0404000","411790105xe30220000","0133103100000010"); // 30-50%
    cuts.AddCutPCMCalo("15910053","0dm00009ab770c00amd0404000","411790105xe30220000","0133103100000010"); // 50-90%
  } else if (trainConfig == 2007){ // 4 cents, with NOC NonLin like MC no OOB Pileup correction for MC with added signal
    cuts.AddCutPCMCalo("10110023","0dm00009ab770c00amd0404000","411790105xe30220000","0133103100000010"); // 00-10%
    cuts.AddCutPCMCalo("11310023","0dm00009ab770c00amd0404000","411790105xe30220000","0133103100000010"); // 10-30%
    cuts.AddCutPCMCalo("13510023","0dm00009ab770c00amd0404000","411790105xe30220000","0133103100000010"); // 30-50%
    cuts.AddCutPCMCalo("15910023","0dm00009ab770c00amd0404000","411790105xe30220000","0133103100000010"); // 50-90%
  } else if (trainConfig == 2008){ // 4 cents, with NOC NonLin like MC no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutPCMCalo("10110023","0dm00009ab770c00amd0404000","411790105xe30220000","0133103100000010"); // 00-10%
    cuts.AddCutPCMCalo("11310023","0dm00009ab770c00amd0404000","411790105xe30220000","0133103100000010"); // 10-30%
    cuts.AddCutPCMCalo("13510023","0dm00009ab770c00amd0404000","411790105xe30220000","0133103100000010"); // 30-50%
    cuts.AddCutPCMCalo("15910023","0dm00009ab770c00amd0404000","411790105xe30220000","0133103100000010"); // 50-90%
  } else if (trainConfig == 2009){ // 4 cents, with NOC NonLin like MC no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutPCMCalo("10110023","0dm00009ab770c00amd0404000","411790105xe30220000","0133103100000010"); // 00-10%
    cuts.AddCutPCMCalo("11310023","0dm00009ab770c00amd0404000","411790105xe30220000","0133103100000010"); // 10-30%
    cuts.AddCutPCMCalo("13510023","0dm00009ab770c00amd0404000","411790105xe30220000","0133103100000010"); // 30-50%
    cuts.AddCutPCMCalo("15910023","0dm00009ab770c00amd0404000","411790105xe30220000","0133103100000010"); // 50-90%
  } else if (trainConfig == 2010){ // 4 cents with mean N matched tracks per cluster data
    cuts.AddCutPCMCalo("10110e03","0dm00009ab770c00amd0404000","411790105te30220000","0133103100000010"); // 00-10%
    cuts.AddCutPCMCalo("11310e03","0dm00009ab770c00amd0404000","411790105te30220000","0133103100000010"); // 10-30%
    cuts.AddCutPCMCalo("13510e03","0dm00009ab770c00amd0404000","411790105te30220000","0133103100000010"); // 30-50%
    cuts.AddCutPCMCalo("15910e03","0dm00009ab770c00amd0404000","411790105te30220000","0133103100000010"); // 50-90%
  } else if (trainConfig == 2011){ // 4 cents, with mean N matched tracks per cluster MC no OOB Pileup correction for MC
    cuts.AddCutPCMCalo("10110053","0dm00009ab770c00amd0404000","411790105te30220000","0133103100000010"); // 00-10%
    cuts.AddCutPCMCalo("11310053","0dm00009ab770c00amd0404000","411790105te30220000","0133103100000010"); // 10-30%
    cuts.AddCutPCMCalo("13510053","0dm00009ab770c00amd0404000","411790105te30220000","0133103100000010"); // 30-50%
    cuts.AddCutPCMCalo("15910053","0dm00009ab770c00amd0404000","411790105te30220000","0133103100000010"); // 50-90%
  } else if (trainConfig == 2012){ // 4 cents, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal
    cuts.AddCutPCMCalo("10110023","0dm00009ab770c00amd0404000","411790105te30220000","0133103100000010"); // 00-10%
    cuts.AddCutPCMCalo("11310023","0dm00009ab770c00amd0404000","411790105te30220000","0133103100000010"); // 10-30%
    cuts.AddCutPCMCalo("13510023","0dm00009ab770c00amd0404000","411790105te30220000","0133103100000010"); // 30-50%
    cuts.AddCutPCMCalo("15910023","0dm00009ab770c00amd0404000","411790105te30220000","0133103100000010"); // 50-90%
  } else if (trainConfig == 2013){ // 4 cents, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutPCMCalo("10110023","0dm00009ab770c00amd0404000","411790105te30220000","0133103100000010"); // 00-10%
    cuts.AddCutPCMCalo("11310023","0dm00009ab770c00amd0404000","411790105te30220000","0133103100000010"); // 10-30%
    cuts.AddCutPCMCalo("13510023","0dm00009ab770c00amd0404000","411790105te30220000","0133103100000010"); // 30-50%
    cuts.AddCutPCMCalo("15910023","0dm00009ab770c00amd0404000","411790105te30220000","0133103100000010"); // 50-90%
  } else if (trainConfig == 2014){ // 4 cents, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutPCMCalo("10110023","0dm00009ab770c00amd0404000","411790105te30220000","0133103100000010"); // 00-10%
    cuts.AddCutPCMCalo("11310023","0dm00009ab770c00amd0404000","411790105te30220000","0133103100000010"); // 10-30%
    cuts.AddCutPCMCalo("13510023","0dm00009ab770c00amd0404000","411790105te30220000","0133103100000010"); // 30-50%
    cuts.AddCutPCMCalo("15910023","0dm00009ab770c00amd0404000","411790105te30220000","0133103100000010"); // 50-90%
  } else if (trainConfig == 2015){ // 4 cents with mean N matched tracks per cluster data
    cuts.AddCutPCMCalo("10110e03","0dm00009ab770c00amd0404000","411790105ye30220000","0133103100000010"); // 00-10%
    cuts.AddCutPCMCalo("11310e03","0dm00009ab770c00amd0404000","411790105ye30220000","0133103100000010"); // 10-30%
    cuts.AddCutPCMCalo("13510e03","0dm00009ab770c00amd0404000","411790105ye30220000","0133103100000010"); // 30-50%
    cuts.AddCutPCMCalo("15910e03","0dm00009ab770c00amd0404000","411790105ye30220000","0133103100000010"); // 50-90%
  } else if (trainConfig == 2016){ // 4 cents, with mean N matched tracks per cluster MC no OOB Pileup correction for MC
    cuts.AddCutPCMCalo("10110053","0dm00009ab770c00amd0404000","411790105ye30220000","0133103100000010"); // 00-10%
    cuts.AddCutPCMCalo("11310053","0dm00009ab770c00amd0404000","411790105ye30220000","0133103100000010"); // 10-30%
    cuts.AddCutPCMCalo("13510053","0dm00009ab770c00amd0404000","411790105ye30220000","0133103100000010"); // 30-50%
    cuts.AddCutPCMCalo("15910053","0dm00009ab770c00amd0404000","411790105ye30220000","0133103100000010"); // 50-90%
  // with rotation background
  } else if (trainConfig == 2050){ // 4 cents
    cuts.AddCutPCMCalo("10110e03","0dm00009ab770c00amd0404000","411790105ke30220000","0s331031000000d0"); // 00-10%
    cuts.AddCutPCMCalo("11310e03","0dm00009ab770c00amd0404000","411790105ke30220000","0s331031000000d0"); // 10-30%
    cuts.AddCutPCMCalo("13510e03","0dm00009ab770c00amd0404000","411790105ke30220000","0s331031000000d0"); // 30-50%
    cuts.AddCutPCMCalo("15910e03","0dm00009ab770c00amd0404000","411790105ke30220000","0s331031000000d0"); // 50-90%
  } else if (trainConfig == 2051){ // 4 cents, no OOB Pileup correction for MC
    cuts.AddCutPCMCalo("10110053","0dm00009ab770c00amd0404000","411790105ke30220000","0s331031000000d0"); // 00-10%
    cuts.AddCutPCMCalo("11310053","0dm00009ab770c00amd0404000","411790105ke30220000","0s331031000000d0"); // 10-30%
    cuts.AddCutPCMCalo("13510053","0dm00009ab770c00amd0404000","411790105ke30220000","0s331031000000d0"); // 30-50%
    cuts.AddCutPCMCalo("15910053","0dm00009ab770c00amd0404000","411790105ke30220000","0s331031000000d0"); // 50-90%
  } else if (trainConfig == 2052){ // 4 cents, no OOB Pileup correction for MC with added signal
    cuts.AddCutPCMCalo("10110023","0dm00009ab770c00amd0404000","411790105ke30220000","0s331031000000d0"); // 00-10%
    cuts.AddCutPCMCalo("11310023","0dm00009ab770c00amd0404000","411790105ke30220000","0s331031000000d0"); // 10-30%
    cuts.AddCutPCMCalo("13510023","0dm00009ab770c00amd0404000","411790105ke30220000","0s331031000000d0"); // 30-50%
    cuts.AddCutPCMCalo("15910023","0dm00009ab770c00amd0404000","411790105ke30220000","0s331031000000d0"); // 50-90%
  } else if (trainConfig == 2053){ // 4 cents, no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutPCMCalo("10110023","0dm00009ab770c00amd0404000","411790105ke30220000","0s331031000000d0"); // 00-10%
    cuts.AddCutPCMCalo("11310023","0dm00009ab770c00amd0404000","411790105ke30220000","0s331031000000d0"); // 10-30%
    cuts.AddCutPCMCalo("13510023","0dm00009ab770c00amd0404000","411790105ke30220000","0s331031000000d0"); // 30-50%
    cuts.AddCutPCMCalo("15910023","0dm00009ab770c00amd0404000","411790105ke30220000","0s331031000000d0"); // 50-90%
  } else if (trainConfig == 2054){ // 4 cents, no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutPCMCalo("10110023","0dm00009ab770c00amd0404000","411790105ke30220000","0s331031000000d0"); // 00-10%
    cuts.AddCutPCMCalo("11310023","0dm00009ab770c00amd0404000","411790105ke30220000","0s331031000000d0"); // 10-30%
    cuts.AddCutPCMCalo("13510023","0dm00009ab770c00amd0404000","411790105ke30220000","0s331031000000d0"); // 30-50%
    cuts.AddCutPCMCalo("15910023","0dm00009ab770c00amd0404000","411790105ke30220000","0s331031000000d0"); // 50-90%
  } else if (trainConfig == 2055){ // 4 cents with NOC NonLin like data
    cuts.AddCutPCMCalo("10110e03","0dm00009ab770c00amd0404000","411790105we30220000","0s331031000000d0"); // 00-10%
    cuts.AddCutPCMCalo("11310e03","0dm00009ab770c00amd0404000","411790105we30220000","0s331031000000d0"); // 10-30%
    cuts.AddCutPCMCalo("13510e03","0dm00009ab770c00amd0404000","411790105we30220000","0s331031000000d0"); // 30-50%
    cuts.AddCutPCMCalo("15910e03","0dm00009ab770c00amd0404000","411790105we30220000","0s331031000000d0"); // 50-90%
  } else if (trainConfig == 2056){ // 4 cents, with NOC NonLin like MC no OOB Pileup correction for MC
    cuts.AddCutPCMCalo("10110053","0dm00009ab770c00amd0404000","411790105xe30220000","0s331031000000d0"); // 00-10%
    cuts.AddCutPCMCalo("11310053","0dm00009ab770c00amd0404000","411790105xe30220000","0s331031000000d0"); // 10-30%
    cuts.AddCutPCMCalo("13510053","0dm00009ab770c00amd0404000","411790105xe30220000","0s331031000000d0"); // 30-50%
    cuts.AddCutPCMCalo("15910053","0dm00009ab770c00amd0404000","411790105xe30220000","0s331031000000d0"); // 50-90%
  } else if (trainConfig == 2057){ // 4 cents, with NOC NonLin like MC no OOB Pileup correction for MC with added signal
    cuts.AddCutPCMCalo("10110023","0dm00009ab770c00amd0404000","411790105xe30220000","0s331031000000d0"); // 00-10%
    cuts.AddCutPCMCalo("11310023","0dm00009ab770c00amd0404000","411790105xe30220000","0s331031000000d0"); // 10-30%
    cuts.AddCutPCMCalo("13510023","0dm00009ab770c00amd0404000","411790105xe30220000","0s331031000000d0"); // 30-50%
    cuts.AddCutPCMCalo("15910023","0dm00009ab770c00amd0404000","411790105xe30220000","0s331031000000d0"); // 50-90%
  } else if (trainConfig == 2058){ // 4 cents, with NOC NonLin like MC no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutPCMCalo("10110023","0dm00009ab770c00amd0404000","411790105xe30220000","0s331031000000d0"); // 00-10%
    cuts.AddCutPCMCalo("11310023","0dm00009ab770c00amd0404000","411790105xe30220000","0s331031000000d0"); // 10-30%
    cuts.AddCutPCMCalo("13510023","0dm00009ab770c00amd0404000","411790105xe30220000","0s331031000000d0"); // 30-50%
    cuts.AddCutPCMCalo("15910023","0dm00009ab770c00amd0404000","411790105xe30220000","0s331031000000d0"); // 50-90%
  } else if (trainConfig == 2059){ // 4 cents, with NOC NonLin like MC no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutPCMCalo("10110023","0dm00009ab770c00amd0404000","411790105xe30220000","0s331031000000d0"); // 00-10%
    cuts.AddCutPCMCalo("11310023","0dm00009ab770c00amd0404000","411790105xe30220000","0s331031000000d0"); // 10-30%
    cuts.AddCutPCMCalo("13510023","0dm00009ab770c00amd0404000","411790105xe30220000","0s331031000000d0"); // 30-50%
    cuts.AddCutPCMCalo("15910023","0dm00009ab770c00amd0404000","411790105xe30220000","0s331031000000d0"); // 50-90%
  } else if (trainConfig == 2060){ // 4 cents with mean N matched tracks per cluster data
    cuts.AddCutPCMCalo("10110e03","0dm00009ab770c00amd0404000","411790105te30220000","0s331031000000d0"); // 00-10%
    cuts.AddCutPCMCalo("11310e03","0dm00009ab770c00amd0404000","411790105te30220000","0s331031000000d0"); // 10-30%
    cuts.AddCutPCMCalo("13510e03","0dm00009ab770c00amd0404000","411790105te30220000","0s331031000000d0"); // 30-50%
    cuts.AddCutPCMCalo("15910e03","0dm00009ab770c00amd0404000","411790105te30220000","0s331031000000d0"); // 50-90%
  } else if (trainConfig == 2061){ // 4 cents, with mean N matched tracks per cluster MC no OOB Pileup correction for MC
    cuts.AddCutPCMCalo("10110053","0dm00009ab770c00amd0404000","411790105te30220000","0s331031000000d0"); // 00-10%
    cuts.AddCutPCMCalo("11310053","0dm00009ab770c00amd0404000","411790105te30220000","0s331031000000d0"); // 10-30%
    cuts.AddCutPCMCalo("13510053","0dm00009ab770c00amd0404000","411790105te30220000","0s331031000000d0"); // 30-50%
    cuts.AddCutPCMCalo("15910053","0dm00009ab770c00amd0404000","411790105te30220000","0s331031000000d0"); // 50-90%
  } else if (trainConfig == 2062){ // 4 cents, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal
    cuts.AddCutPCMCalo("10110023","0dm00009ab770c00amd0404000","411790105te30220000","0s331031000000d0"); // 00-10%
    cuts.AddCutPCMCalo("11310023","0dm00009ab770c00amd0404000","411790105te30220000","0s331031000000d0"); // 10-30%
    cuts.AddCutPCMCalo("13510023","0dm00009ab770c00amd0404000","411790105te30220000","0s331031000000d0"); // 30-50%
    cuts.AddCutPCMCalo("15910023","0dm00009ab770c00amd0404000","411790105te30220000","0s331031000000d0"); // 50-90%
  } else if (trainConfig == 2063){ // 4 cents, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutPCMCalo("10110023","0dm00009ab770c00amd0404000","411790105te30220000","0s331031000000d0"); // 00-10%
    cuts.AddCutPCMCalo("11310023","0dm00009ab770c00amd0404000","411790105te30220000","0s331031000000d0"); // 10-30%
    cuts.AddCutPCMCalo("13510023","0dm00009ab770c00amd0404000","411790105te30220000","0s331031000000d0"); // 30-50%
    cuts.AddCutPCMCalo("15910023","0dm00009ab770c00amd0404000","411790105te30220000","0s331031000000d0"); // 50-90%
  } else if (trainConfig == 2064){ // 4 cents, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutPCMCalo("10110023","0dm00009ab770c00amd0404000","411790105te30220000","0s331031000000d0"); // 00-10%
    cuts.AddCutPCMCalo("11310023","0dm00009ab770c00amd0404000","411790105te30220000","0s331031000000d0"); // 10-30%
    cuts.AddCutPCMCalo("13510023","0dm00009ab770c00amd0404000","411790105te30220000","0s331031000000d0"); // 30-50%
    cuts.AddCutPCMCalo("15910023","0dm00009ab770c00amd0404000","411790105te30220000","0s331031000000d0"); // 50-90%
  } else if (trainConfig == 2065){ // 4 cents with mean N matched tracks per cluster data
    cuts.AddCutPCMCalo("10110e03","0dm00009ab770c00amd0404000","411790105ye30220000","0s33103100000010"); // 00-10%
    cuts.AddCutPCMCalo("11310e03","0dm00009ab770c00amd0404000","411790105ye30220000","0s33103100000010"); // 10-30%
    cuts.AddCutPCMCalo("13510e03","0dm00009ab770c00amd0404000","411790105ye30220000","0s33103100000010"); // 30-50%
    cuts.AddCutPCMCalo("15910e03","0dm00009ab770c00amd0404000","411790105ye30220000","0s33103100000010"); // 50-90%
  } else if (trainConfig == 2066){ // 4 cents, with mean N matched tracks per cluster MC no OOB Pileup correction for MC
    cuts.AddCutPCMCalo("10110053","0dm00009ab770c00amd0404000","411790105ye30220000","0s33103100000010"); // 00-10%
    cuts.AddCutPCMCalo("11310053","0dm00009ab770c00amd0404000","411790105ye30220000","0s33103100000010"); // 10-30%
    cuts.AddCutPCMCalo("13510053","0dm00009ab770c00amd0404000","411790105ye30220000","0s33103100000010"); // 30-50%
    cuts.AddCutPCMCalo("15910053","0dm00009ab770c00amd0404000","411790105ye30220000","0s33103100000010"); // 50-90%
  // **********************************************************************************************************
  // **************************** PCM-EMC configurations PbPb run 2 2015 pass 3 cent ******************************
  // **********************************************************************************************************
  // Standard: EMCal + DCal, only 2ndary cluster track matching, mixed event
  } else if (trainConfig == 2100){ // central
    cuts.AddCutPCMCalo("10110e03","0dm00009ab770c00amd0404000","411790105ke30220000","0133103100000010"); // 00-10%
  } else if (trainConfig == 2101){ // central, no OOB Pileup correction for MC
    cuts.AddCutPCMCalo("10110053","0dm00009ab770c00amd0404000","411790105ke30220000","0133103100000010"); // 00-10%
  } else if (trainConfig == 2102){ // central, no OOB Pileup correction for MC with added signal
    cuts.AddCutPCMCalo("10110023","0dm00009ab770c00amd0404000","411790105ke30220000","0133103100000010"); // 00-10%
  } else if (trainConfig == 2103){ // central, no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutPCMCalo("10110023","0dm00009ab770c00amd0404000","411790105ke30220000","0133103100000010"); // 00-10%
  } else if (trainConfig == 2104){ // central, no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutPCMCalo("10110023","0dm00009ab770c00amd0404000","411790105ke30220000","0133103100000010"); // 00-10%
  } else if (trainConfig == 2105){ // central with NOC NonLin like data
    cuts.AddCutPCMCalo("10110e03","0dm00009ab770c00amd0404000","411790105we30220000","0133103100000010"); // 00-10%
  } else if (trainConfig == 2106){ // central, with NOC NonLin like MC no OOB Pileup correction for MC
    cuts.AddCutPCMCalo("10110053","0dm00009ab770c00amd0404000","411790105xe30220000","0133103100000010"); // 00-10%
  } else if (trainConfig == 2107){ // central, with NOC NonLin like MC no OOB Pileup correction for MC with added signal
    cuts.AddCutPCMCalo("10110023","0dm00009ab770c00amd0404000","411790105xe30220000","0133103100000010"); // 00-10%
  } else if (trainConfig == 2108){ // central, with NOC NonLin like MC no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutPCMCalo("10110023","0dm00009ab770c00amd0404000","411790105xe30220000","0133103100000010"); // 00-10%
  } else if (trainConfig == 2109){ // central, with NOC NonLin like MC no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutPCMCalo("10110023","0dm00009ab770c00amd0404000","411790105xe30220000","0133103100000010"); // 00-10%
  } else if (trainConfig == 2110){ // central with mean N matched tracks per cluster data
    cuts.AddCutPCMCalo("10110e03","0dm00009ab770c00amd0404000","411790105te30220000","0133103100000010"); // 00-10%
  } else if (trainConfig == 2111){ // central, with mean N matched tracks per cluster MC no OOB Pileup correction for MC
    cuts.AddCutPCMCalo("10110053","0dm00009ab770c00amd0404000","411790105te30220000","0133103100000010"); // 00-10%
  } else if (trainConfig == 2112){ // central, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal
    cuts.AddCutPCMCalo("10110023","0dm00009ab770c00amd0404000","411790105te30220000","0133103100000010"); // 00-10%
  } else if (trainConfig == 2113){ // central, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutPCMCalo("10110023","0dm00009ab770c00amd0404000","411790105te30220000","0133103100000010"); // 00-10%
  } else if (trainConfig == 2114){ // central, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutPCMCalo("10110023","0dm00009ab770c00amd0404000","411790105te30220000","0133103100000010"); // 00-10%
  } else if (trainConfig == 2115){ // central with mean N matched tracks per cluster data
    cuts.AddCutPCMCalo("10110e03","0dm00009ab770c00amd0404000","411790105ye30220000","0133103100000010"); // 00-10%
  } else if (trainConfig == 2116){ // central, with mean N matched tracks per cluster MC no OOB Pileup correction for MC
    cuts.AddCutPCMCalo("10110053","0dm00009ab770c00amd0404000","411790105ye30220000","0133103100000010"); // 00-10%
  // with rotation background
  } else if (trainConfig == 2150){ // central
    cuts.AddCutPCMCalo("10110e03","0dm00009ab770c00amd0404000","411790105ke30220000","0s331031000000d0"); // 00-10%
  } else if (trainConfig == 2151){ // central, no OOB Pileup correction for MC
    cuts.AddCutPCMCalo("10110053","0dm00009ab770c00amd0404000","411790105ke30220000","0s331031000000d0"); // 00-10%
  } else if (trainConfig == 2152){ // central, no OOB Pileup correction for MC with added signal
    cuts.AddCutPCMCalo("10110023","0dm00009ab770c00amd0404000","411790105ke30220000","0s331031000000d0"); // 00-10%
  } else if (trainConfig == 2153){ // central, no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutPCMCalo("10110023","0dm00009ab770c00amd0404000","411790105ke30220000","0s331031000000d0"); // 00-10%
  } else if (trainConfig == 2154){ // central, no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutPCMCalo("10110023","0dm00009ab770c00amd0404000","411790105ke30220000","0s331031000000d0"); // 00-10%
  } else if (trainConfig == 2155){ // central with NOC NonLin like data
    cuts.AddCutPCMCalo("10110e03","0dm00009ab770c00amd0404000","411790105we30220000","0s331031000000d0"); // 00-10%
  } else if (trainConfig == 2156){ // central, with NOC NonLin like MC no OOB Pileup correction for MC
    cuts.AddCutPCMCalo("10110053","0dm00009ab770c00amd0404000","411790105xe30220000","0s331031000000d0"); // 00-10%
  } else if (trainConfig == 2157){ // central, with NOC NonLin like MC no OOB Pileup correction for MC with added signal
    cuts.AddCutPCMCalo("10110023","0dm00009ab770c00amd0404000","411790105xe30220000","0s331031000000d0"); // 00-10%
  } else if (trainConfig == 2158){ // central, with NOC NonLin like MC no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutPCMCalo("10110023","0dm00009ab770c00amd0404000","411790105xe30220000","0s331031000000d0"); // 00-10%
  } else if (trainConfig == 2159){ // central, with NOC NonLin like MC no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutPCMCalo("10110023","0dm00009ab770c00amd0404000","411790105xe30220000","0s331031000000d0"); // 00-10%
  } else if (trainConfig == 2160){ // central with mean N matched tracks per cluster data
    cuts.AddCutPCMCalo("10110e03","0dm00009ab770c00amd0404000","411790105te30220000","0s331031000000d0"); // 00-10%
  } else if (trainConfig == 2161){ // central, with mean N matched tracks per cluster MC no OOB Pileup correction for MC
    cuts.AddCutPCMCalo("10110053","0dm00009ab770c00amd0404000","411790105te30220000","0s331031000000d0"); // 00-10%
  } else if (trainConfig == 2162){ // central, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal
    cuts.AddCutPCMCalo("10110023","0dm00009ab770c00amd0404000","411790105te30220000","0s331031000000d0"); // 00-10%
  } else if (trainConfig == 2163){ // central, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutPCMCalo("10110023","0dm00009ab770c00amd0404000","411790105te30220000","0s331031000000d0"); // 00-10%
  } else if (trainConfig == 2164){ // central, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutPCMCalo("10110023","0dm00009ab770c00amd0404000","411790105te30220000","0s331031000000d0"); // 00-10%
  } else if (trainConfig == 2165){ // central with mean N matched tracks per cluster data
    cuts.AddCutPCMCalo("10110e03","0dm00009ab770c00amd0404000","411790105ye30220000","0s33103100000010"); // 00-10%
  } else if (trainConfig == 2166){ // central, with mean N matched tracks per cluster MC no OOB Pileup correction for MC
    cuts.AddCutPCMCalo("10110053","0dm00009ab770c00amd0404000","411790105ye30220000","0s33103100000010"); // 00-10%
  // **********************************************************************************************************
  // ************************** PCM-EMC configurations PbPb run 2 2015 pass 3 semicent ****************************
  // **********************************************************************************************************
  // Standard: EMCal + DCal, only 2ndary cluster track matching, mixed event
  } else if (trainConfig == 2300){ // semicent
    cuts.AddCutPCMCalo("11310e03","0dm00009ab770c00amd0404000","411790105ke30220000","0133103100000010"); // 10-30%
  } else if (trainConfig == 2301){ // semicent, no OOB Pileup correction for MC
    cuts.AddCutPCMCalo("11310053","0dm00009ab770c00amd0404000","411790105ke30220000","0133103100000010"); // 10-30%
  } else if (trainConfig == 2302){ // semicent, no OOB Pileup correction for MC with added signal
    cuts.AddCutPCMCalo("11310023","0dm00009ab770c00amd0404000","411790105ke30220000","0133103100000010"); // 10-30%
  } else if (trainConfig == 2303){ // semicent, no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutPCMCalo("11310023","0dm00009ab770c00amd0404000","411790105ke30220000","0133103100000010"); // 10-30%
  } else if (trainConfig == 2304){ // semicent, no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutPCMCalo("11310023","0dm00009ab770c00amd0404000","411790105ke30220000","0133103100000010"); // 10-30%
  } else if (trainConfig == 2305){ // semicent with NOC NonLin like data
    cuts.AddCutPCMCalo("11310e03","0dm00009ab770c00amd0404000","411790105we30220000","0133103100000010"); // 10-30%
  } else if (trainConfig == 2306){ // semicent, with NOC NonLin like MC no OOB Pileup correction for MC
    cuts.AddCutPCMCalo("11310053","0dm00009ab770c00amd0404000","411790105xe30220000","0133103100000010"); // 10-30%
  } else if (trainConfig == 2307){ // semicent, with NOC NonLin like MC no OOB Pileup correction for MC with added signal
    cuts.AddCutPCMCalo("11310023","0dm00009ab770c00amd0404000","411790105xe30220000","0133103100000010"); // 10-30%
  } else if (trainConfig == 2308){ // semicent, with NOC NonLin like MC no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutPCMCalo("11310023","0dm00009ab770c00amd0404000","411790105xe30220000","0133103100000010"); // 10-30%
  } else if (trainConfig == 2309){ // semicent, with NOC NonLin like MC no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutPCMCalo("11310023","0dm00009ab770c00amd0404000","411790105xe30220000","0133103100000010"); // 10-30%
  } else if (trainConfig == 2310){ // semicent with mean N matched tracks per cluster data
    cuts.AddCutPCMCalo("11310e03","0dm00009ab770c00amd0404000","411790105te30220000","0133103100000010"); // 10-30%
  } else if (trainConfig == 2311){ // semicent, with mean N matched tracks per cluster MC no OOB Pileup correction for MC
    cuts.AddCutPCMCalo("11310053","0dm00009ab770c00amd0404000","411790105te30220000","0133103100000010"); // 10-30%
  } else if (trainConfig == 2312){ // semicent, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal
    cuts.AddCutPCMCalo("11310023","0dm00009ab770c00amd0404000","411790105te30220000","0133103100000010"); // 10-30%
  } else if (trainConfig == 2313){ // semicent, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutPCMCalo("11310023","0dm00009ab770c00amd0404000","411790105te30220000","0133103100000010"); // 10-30%
  } else if (trainConfig == 2314){ // semicent, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutPCMCalo("11310023","0dm00009ab770c00amd0404000","411790105te30220000","0133103100000010"); // 10-30%
  } else if (trainConfig == 2315){ // semicent with mean N matched tracks per cluster data
    cuts.AddCutPCMCalo("11310e03","0dm00009ab770c00amd0404000","411790105ye30220000","0133103100000010"); // 10-30%
  } else if (trainConfig == 2316){ // semicent, with mean N matched tracks per cluster MC no OOB Pileup correction for MC
    cuts.AddCutPCMCalo("11310053","0dm00009ab770c00amd0404000","411790105ye30220000","0133103100000010"); // 10-30%
  // with rotation background
  } else if (trainConfig == 2350){ // semicent
    cuts.AddCutPCMCalo("11310e03","0dm00009ab770c00amd0404000","411790105ke30220000","0s33103100000010"); // 10-30%
  } else if (trainConfig == 2351){ // semicent, no OOB Pileup correction for MC
    cuts.AddCutPCMCalo("11310053","0dm00009ab770c00amd0404000","411790105ke30220000","0s33103100000010"); // 10-30%
  } else if (trainConfig == 2352){ // semicent, no OOB Pileup correction for MC with added signal
    cuts.AddCutPCMCalo("11310023","0dm00009ab770c00amd0404000","411790105ke30220000","0s33103100000010"); // 10-30%
  } else if (trainConfig == 2353){ // semicent, no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutPCMCalo("11310023","0dm00009ab770c00amd0404000","411790105ke30220000","0s33103100000010"); // 10-30%
  } else if (trainConfig == 2354){ // semicent, no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutPCMCalo("11310023","0dm00009ab770c00amd0404000","411790105ke30220000","0s33103100000010"); // 10-30%
  } else if (trainConfig == 2355){ // semicent with NOC NonLin like data
    cuts.AddCutPCMCalo("11310e03","0dm00009ab770c00amd0404000","411790105we30220000","0s33103100000010"); // 10-30%
  } else if (trainConfig == 2356){ // semicent, with NOC NonLin like MC no OOB Pileup correction for MC
    cuts.AddCutPCMCalo("11310053","0dm00009ab770c00amd0404000","411790105xe30220000","0s33103100000010"); // 10-30%
  } else if (trainConfig == 2357){ // semicent, with NOC NonLin like MC no OOB Pileup correction for MC with added signal
    cuts.AddCutPCMCalo("11310023","0dm00009ab770c00amd0404000","411790105xe30220000","0s33103100000010"); // 10-30%
  } else if (trainConfig == 2358){ // semicent, with NOC NonLin like MC no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutPCMCalo("11310023","0dm00009ab770c00amd0404000","411790105xe30220000","0s33103100000010"); // 10-30%
  } else if (trainConfig == 2359){ // semicent, with NOC NonLin like MC no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutPCMCalo("11310023","0dm00009ab770c00amd0404000","411790105xe30220000","0s33103100000010"); // 10-30%
  } else if (trainConfig == 2360){ // semicent with mean N matched tracks per cluster data
    cuts.AddCutPCMCalo("11310e03","0dm00009ab770c00amd0404000","411790105te30220000","0s33103100000010"); // 10-30%
  } else if (trainConfig == 2361){ // semicent, with mean N matched tracks per cluster MC no OOB Pileup correction for MC
    cuts.AddCutPCMCalo("11310053","0dm00009ab770c00amd0404000","411790105te30220000","0s33103100000010"); // 10-30%
  } else if (trainConfig == 2362){ // semicent, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal
    cuts.AddCutPCMCalo("11310023","0dm00009ab770c00amd0404000","411790105te30220000","0s33103100000010"); // 10-30%
  } else if (trainConfig == 2363){ // semicent, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutPCMCalo("11310023","0dm00009ab770c00amd0404000","411790105te30220000","0s33103100000010"); // 10-30%
  } else if (trainConfig == 2364){ // semicent, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutPCMCalo("11310023","0dm00009ab770c00amd0404000","411790105te30220000","0s33103100000010"); // 10-30%
  } else if (trainConfig == 2365){ // semicent with mean N matched tracks per cluster data
    cuts.AddCutPCMCalo("11310e03","0dm00009ab770c00amd0404000","411790105ye30220000","0s33103100000010"); // 10-30%
  } else if (trainConfig == 2366){ // semicent, with mean N matched tracks per cluster MC no OOB Pileup correction for MC
    cuts.AddCutPCMCalo("11310053","0dm00009ab770c00amd0404000","411790105ye30220000","0s33103100000010"); // 10-30%
  // **********************************************************************************************************
  // *********************** PCM-EMC configurations PbPb run 2 2015 pass 3 semiperipheral *************************
  // **********************************************************************************************************
  // Standard: EMCal + DCal, only 2ndary cluster track matching, mixed event
  } else if (trainConfig == 2500){ // semiperipheral
    cuts.AddCutPCMCalo("13510e03","0dm00009ab770c00amd0404000","411790105ke30220000","0133103100000010"); // 30-50%
  } else if (trainConfig == 2501){ // semiperipheral, no OOB Pileup correction for MC
    cuts.AddCutPCMCalo("13510053","0dm00009ab770c00amd0404000","411790105ke30220000","0133103100000010"); // 30-50%
  } else if (trainConfig == 2502){ // semiperipheral, no OOB Pileup correction for MC with added signal
    cuts.AddCutPCMCalo("13510023","0dm00009ab770c00amd0404000","411790105ke30220000","0133103100000010"); // 30-50%
  } else if (trainConfig == 2503){ // semiperipheral, no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutPCMCalo("13510023","0dm00009ab770c00amd0404000","411790105ke30220000","0133103100000010"); // 30-50%
  } else if (trainConfig == 2504){ // semiperipheral, no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutPCMCalo("13510023","0dm00009ab770c00amd0404000","411790105ke30220000","0133103100000010"); // 30-50%
  } else if (trainConfig == 2505){ // semiperipheral with NOC NonLin like data
    cuts.AddCutPCMCalo("13510e03","0dm00009ab770c00amd0404000","411790105we30220000","0133103100000010"); // 30-50%
  } else if (trainConfig == 2506){ // semiperipheral, with NOC NonLin like MC no OOB Pileup correction for MC
    cuts.AddCutPCMCalo("13510053","0dm00009ab770c00amd0404000","411790105xe30220000","0133103100000010"); // 30-50%
  } else if (trainConfig == 2507){ // semiperipheral, with NOC NonLin like MC no OOB Pileup correction for MC with added signal
    cuts.AddCutPCMCalo("13510023","0dm00009ab770c00amd0404000","411790105xe30220000","0133103100000010"); // 30-50%
  } else if (trainConfig == 2508){ // semiperipheral, with NOC NonLin like MC no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutPCMCalo("13510023","0dm00009ab770c00amd0404000","411790105xe30220000","0133103100000010"); // 30-50%
  } else if (trainConfig == 2509){ // semiperipheral, with NOC NonLin like MC no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutPCMCalo("13510023","0dm00009ab770c00amd0404000","411790105xe30220000","0133103100000010"); // 30-50%
  } else if (trainConfig == 2510){ // semiperipheral with mean N matched tracks per cluster data
    cuts.AddCutPCMCalo("13510e03","0dm00009ab770c00amd0404000","411790105te30220000","0133103100000010"); // 30-50%
  } else if (trainConfig == 2511){ // semiperipheral, with mean N matched tracks per cluster MC no OOB Pileup correction for MC
    cuts.AddCutPCMCalo("13510053","0dm00009ab770c00amd0404000","411790105te30220000","0133103100000010"); // 30-50%
  } else if (trainConfig == 2512){ // semiperipheral, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal
    cuts.AddCutPCMCalo("13510023","0dm00009ab770c00amd0404000","411790105te30220000","0133103100000010"); // 30-50%
  } else if (trainConfig == 2513){ // semiperipheral, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutPCMCalo("13510023","0dm00009ab770c00amd0404000","411790105te30220000","0133103100000010"); // 30-50%
  } else if (trainConfig == 2514){ // semiperipheral, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutPCMCalo("13510023","0dm00009ab770c00amd0404000","411790105te30220000","0133103100000010"); // 30-50%
  } else if (trainConfig == 2515){ // semiperipheral with mean N matched tracks per cluster data
    cuts.AddCutPCMCalo("13510e03","0dm00009ab770c00amd0404000","411790105ye30220000","0133103100000010"); // 30-50%
  } else if (trainConfig == 2515){ // semiperipheral, with mean N matched tracks per cluster MC no OOB Pileup correction for MC
    cuts.AddCutPCMCalo("13510053","0dm00009ab770c00amd0404000","411790105ye30220000","0133103100000010"); // 30-50%
  // with rotation background
  } else if (trainConfig == 2550){ // semiperipheral
    cuts.AddCutPCMCalo("13510e03","0dm00009ab770c00amd0404000","411790105ke30220000","0s33103100000010"); // 30-50%
  } else if (trainConfig == 2551){ // semiperipheral, no OOB Pileup correction for MC
    cuts.AddCutPCMCalo("13510053","0dm00009ab770c00amd0404000","411790105ke30220000","0s33103100000010"); // 30-50%
  } else if (trainConfig == 2552){ // semiperipheral, no OOB Pileup correction for MC with added signal
    cuts.AddCutPCMCalo("13510023","0dm00009ab770c00amd0404000","411790105ke30220000","0s33103100000010"); // 30-50%
  } else if (trainConfig == 2553){ // semiperipheral, no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutPCMCalo("13510023","0dm00009ab770c00amd0404000","411790105ke30220000","0s33103100000010"); // 30-50%
  } else if (trainConfig == 2554){ // semiperipheral, no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutPCMCalo("13510023","0dm00009ab770c00amd0404000","411790105ke30220000","0s33103100000010"); // 30-50%
  } else if (trainConfig == 2555){ // semiperipheral with NOC NonLin like data
    cuts.AddCutPCMCalo("13510e03","0dm00009ab770c00amd0404000","411790105we30220000","0s33103100000010"); // 30-50%
  } else if (trainConfig == 2556){ // semiperipheral, with NOC NonLin like MC no OOB Pileup correction for MC
    cuts.AddCutPCMCalo("13510053","0dm00009ab770c00amd0404000","411790105xe30220000","0s33103100000010"); // 30-50%
  } else if (trainConfig == 2557){ // semiperipheral, with NOC NonLin like MC no OOB Pileup correction for MC with added signal
    cuts.AddCutPCMCalo("13510023","0dm00009ab770c00amd0404000","411790105xe30220000","0s33103100000010"); // 30-50%
  } else if (trainConfig == 2558){ // semiperipheral, with NOC NonLin like MC no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutPCMCalo("13510023","0dm00009ab770c00amd0404000","411790105xe30220000","0s33103100000010"); // 30-50%
  } else if (trainConfig == 2559){ // semiperipheral, with NOC NonLin like MC no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutPCMCalo("13510023","0dm00009ab770c00amd0404000","411790105xe30220000","0s33103100000010"); // 30-50%
  } else if (trainConfig == 2560){ // semiperipheral with mean N matched tracks per cluster data
    cuts.AddCutPCMCalo("13510e03","0dm00009ab770c00amd0404000","411790105te30220000","0s33103100000010"); // 30-50%
  } else if (trainConfig == 2561){ // semiperipheral, with mean N matched tracks per cluster MC no OOB Pileup correction for MC
    cuts.AddCutPCMCalo("13510053","0dm00009ab770c00amd0404000","411790105te30220000","0s33103100000010"); // 30-50%
  } else if (trainConfig == 2562){ // semiperipheral, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal
    cuts.AddCutPCMCalo("13510023","0dm00009ab770c00amd0404000","411790105te30220000","0s33103100000010"); // 30-50%
  } else if (trainConfig == 2563){ // semiperipheral, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutPCMCalo("13510023","0dm00009ab770c00amd0404000","411790105te30220000","0s33103100000010"); // 30-50%
  } else if (trainConfig == 2564){ // semiperipheral, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutPCMCalo("13510023","0dm00009ab770c00amd0404000","411790105te30220000","0s33103100000010"); // 30-50%
  } else if (trainConfig == 2565){ // semiperipheral with mean N matched tracks per cluster data
    cuts.AddCutPCMCalo("13510e03","0dm00009ab770c00amd0404000","411790105ye30220000","0s33103100000010"); // 30-50%
  } else if (trainConfig == 2566){ // semiperipheral, with mean N matched tracks per cluster MC no OOB Pileup correction for MC
    cuts.AddCutPCMCalo("13510053","0dm00009ab770c00amd0404000","411790105ye30220000","0s33103100000010"); // 30-50%
  // **********************************************************************************************************
  // ************************* PCM-EMC configurations PbPb run 2 2015 pass 3 peripheral ***************************
  // **********************************************************************************************************
  // Standard: EMCal + DCal, only 2ndary cluster track matching, mixed event
  } else if (trainConfig == 2900){ // peripheral
    cuts.AddCutPCMCalo("15910e03","0dm00009ab770c00amd0404000","411790105ke30220000","0133103100000010"); // 50-90%
  } else if (trainConfig == 2901){ // peripheral, no OOB Pileup correction for MC
    cuts.AddCutPCMCalo("15910053","0dm00009ab770c00amd0404000","411790105ke30220000","0133103100000010"); // 50-90%
  } else if (trainConfig == 2902){ // peripheral, no OOB Pileup correction for MC with added signal
    cuts.AddCutPCMCalo("15910023","0dm00009ab770c00amd0404000","411790105ke30220000","0133103100000010"); // 50-90%
  } else if (trainConfig == 2903){ // peripheral, no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutPCMCalo("15910023","0dm00009ab770c00amd0404000","411790105ke30220000","0133103100000010"); // 50-90%
  } else if (trainConfig == 2904){ // peripheral, no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutPCMCalo("15910023","0dm00009ab770c00amd0404000","411790105ke30220000","0133103100000010"); // 50-90%
  } else if (trainConfig == 2905){ // peripheral with NOC NonLin like data
    cuts.AddCutPCMCalo("15910e03","0dm00009ab770c00amd0404000","411790105we30220000","0133103100000010"); // 50-90%
  } else if (trainConfig == 2906){ // peripheral, with NOC NonLin like MC no OOB Pileup correction for MC
    cuts.AddCutPCMCalo("15910053","0dm00009ab770c00amd0404000","411790105xe30220000","0133103100000010"); // 50-90%
  } else if (trainConfig == 2907){ // peripheral, with NOC NonLin like MC no OOB Pileup correction for MC with added signal
    cuts.AddCutPCMCalo("15910023","0dm00009ab770c00amd0404000","411790105xe30220000","0133103100000010"); // 50-90%
  } else if (trainConfig == 2908){ // peripheral, with NOC NonLin like MC no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutPCMCalo("15910023","0dm00009ab770c00amd0404000","411790105xe30220000","0133103100000010"); // 50-90%
  } else if (trainConfig == 2909){ // peripheral, with NOC NonLin like MC no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutPCMCalo("15910023","0dm00009ab770c00amd0404000","411790105xe30220000","0133103100000010"); // 50-90%
  } else if (trainConfig == 2910){ // peripheral with mean N matched tracks per cluster data
    cuts.AddCutPCMCalo("15910e03","0dm00009ab770c00amd0404000","411790105te30220000","0133103100000010"); // 50-90%
  } else if (trainConfig == 2911){ // peripheral, with mean N matched tracks per cluster MC no OOB Pileup correction for MC
    cuts.AddCutPCMCalo("15910053","0dm00009ab770c00amd0404000","411790105te30220000","0133103100000010"); // 50-90%
  } else if (trainConfig == 2912){ // peripheral, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal
    cuts.AddCutPCMCalo("15910023","0dm00009ab770c00amd0404000","411790105te30220000","0133103100000010"); // 50-90%
  } else if (trainConfig == 2913){ // peripheral, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutPCMCalo("15910023","0dm00009ab770c00amd0404000","411790105te30220000","0133103100000010"); // 50-90%
  } else if (trainConfig == 2914){ // peripheral, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutPCMCalo("15910023","0dm00009ab770c00amd0404000","411790105te30220000","0133103100000010"); // 50-90%
  } else if (trainConfig == 2915){ // peripheral with mean N matched tracks per cluster data
    cuts.AddCutPCMCalo("15910e03","0dm00009ab770c00amd0404000","411790105ye30220000","0133103100000010"); // 50-90%
  } else if (trainConfig == 2916){ // peripheral, with mean N matched tracks per cluster MC no OOB Pileup correction for MC
    cuts.AddCutPCMCalo("15910053","0dm00009ab770c00amd0404000","411790105ye30220000","0133103100000010"); // 50-90%
  // with rotation background
  } else if (trainConfig == 2950){ // peripheral
    cuts.AddCutPCMCalo("15910e03","0dm00009ab770c00amd0404000","411790105ke30220000","0s33103100000010"); // 50-90%
  } else if (trainConfig == 2951){ // peripheral, no OOB Pileup correction for MC
    cuts.AddCutPCMCalo("15910053","0dm00009ab770c00amd0404000","411790105ke30220000","0s33103100000010"); // 50-90%
  } else if (trainConfig == 2952){ // peripheral, no OOB Pileup correction for MC with added signal
    cuts.AddCutPCMCalo("15910023","0dm00009ab770c00amd0404000","411790105ke30220000","0s33103100000010"); // 50-90%
  } else if (trainConfig == 2953){ // peripheral, no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutPCMCalo("15910023","0dm00009ab770c00amd0404000","411790105ke30220000","0s33103100000010"); // 50-90%
  } else if (trainConfig == 2954){ // peripheral, no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutPCMCalo("15910023","0dm00009ab770c00amd0404000","411790105ke30220000","0s33103100000010"); // 50-90%
  } else if (trainConfig == 2955){ // peripheral with NOC NonLin like data
    cuts.AddCutPCMCalo("15910e03","0dm00009ab770c00amd0404000","411790105we30220000","0s33103100000010"); // 50-90%
  } else if (trainConfig == 2956){ // peripheral, with NOC NonLin like MC no OOB Pileup correction for MC
    cuts.AddCutPCMCalo("15910053","0dm00009ab770c00amd0404000","411790105xe30220000","0s33103100000010"); // 50-90%
  } else if (trainConfig == 2957){ // peripheral, with NOC NonLin like MC no OOB Pileup correction for MC with added signal
    cuts.AddCutPCMCalo("15910023","0dm00009ab770c00amd0404000","411790105xe30220000","0s33103100000010"); // 50-90%
  } else if (trainConfig == 2958){ // peripheral, with NOC NonLin like MC no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutPCMCalo("15910023","0dm00009ab770c00amd0404000","411790105xe30220000","0s33103100000010"); // 50-90%
  } else if (trainConfig == 2959){ // peripheral, with NOC NonLin like MC no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutPCMCalo("15910023","0dm00009ab770c00amd0404000","411790105xe30220000","0s33103100000010"); // 50-90%
  } else if (trainConfig == 2960){ // peripheral with mean N matched tracks per cluster data
    cuts.AddCutPCMCalo("15910e03","0dm00009ab770c00amd0404000","411790105te30220000","0s33103100000010"); // 50-90%
  } else if (trainConfig == 2961){ // peripheral, with mean N matched tracks per cluster MC no OOB Pileup correction for MC
    cuts.AddCutPCMCalo("15910053","0dm00009ab770c00amd0404000","411790105te30220000","0s33103100000010"); // 50-90%
  } else if (trainConfig == 2962){ // peripheral, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal
    cuts.AddCutPCMCalo("15910023","0dm00009ab770c00amd0404000","411790105te30220000","0s33103100000010"); // 50-90%
  } else if (trainConfig == 2963){ // peripheral, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutPCMCalo("15910023","0dm00009ab770c00amd0404000","411790105te30220000","0s33103100000010"); // 50-90%
  } else if (trainConfig == 2964){ // peripheral, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutPCMCalo("15910023","0dm00009ab770c00amd0404000","411790105te30220000","0s33103100000010"); // 50-90%
  } else if (trainConfig == 2965){ // peripheral with mean N matched tracks per cluster data
    cuts.AddCutPCMCalo("15910e03","0dm00009ab770c00amd0404000","411790105ye30220000","0s33103100000010"); // 50-90%
  } else if (trainConfig == 2966){ // peripheral, with mean N matched tracks per cluster MC no OOB Pileup correction for MC
    cuts.AddCutPCMCalo("15910053","0dm00009ab770c00amd0404000","411790105ye30220000","0s33103100000010"); // 50-90%
  
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

  TList *EventCutList     = new TList();
  TList *ConvCutList      = new TList();
  TList *ClusterCutList   = new TList();
  TList *MesonCutList     = new TList();

  TList *HeaderList       = new TList();
  if (generatorName.CompareTo("LHC13d2")==0){
    TObjString *Header1 = new TObjString("pi0_1");
    HeaderList->Add(Header1);
  //    TObjString *Header3 = new TObjString("eta_2");
  //    HeaderList->Add(Header3);

  } else if (generatorName.CompareTo("LHC12a17x_fix")==0){
    TObjString *Header1 = new TObjString("PARAM");
    HeaderList->Add(Header1);
  } else if (generatorName.CompareTo("LHC14a1a")==0){
    if (headerSelectionInt == 1){
      TObjString *Header1 = new TObjString("pi0_1");
      HeaderList->Add(Header1);
    } else if (headerSelectionInt == 2){
      TObjString *Header1 = new TObjString("eta_2");
      HeaderList->Add(Header1);
    } else if (headerSelectionInt == 3){
      TString nameHeaders[2]    = { "pi0_1", "eta_2" };
      for (Int_t iHead = 0; iHead < 2; iHead++ ){
        TObjString *Header = new TObjString(nameHeaders[iHead]);
        HeaderList->Add(Header);
      }
    } else if (headerSelectionInt == 4){
      TObjString *Header1 = new TObjString("pi0EMC_3");
      HeaderList->Add(Header1);
    } else if (headerSelectionInt == 5){
      TObjString *Header1 = new TObjString("etaEMC_5");
      HeaderList->Add(Header1);
    } else if (headerSelectionInt == 6){
      TString nameHeaders[2]    = { "pi0EMC_3", "etaEMC_5" };
      for (Int_t iHead = 0; iHead < 2; iHead++ ){
        TObjString *Header = new TObjString(nameHeaders[iHead]);
        HeaderList->Add(Header);
      }
    } else if (headerSelectionInt == 7){
      TString nameHeaders[4]    = { "pi0_1", "eta_2", "pi0EMC_3", "etaEMC_5" };
      for (Int_t iHead = 0; iHead < 4; iHead++ ){
        TObjString *Header = new TObjString(nameHeaders[iHead]);
        HeaderList->Add(Header);
      }
    } else if (headerSelectionInt == 8){
      TObjString *Header1 = new TObjString("gEMCPhoton_7");
      HeaderList->Add(Header1);
    } else if (headerSelectionInt == 9){
      TString nameHeaders[10]   = { "Pythia_Jets_PtHard_1_10", "Pythia_Jets_PtHard_2_10", "Pythia_Jets_PtHard_3_10", "Pythia_Jets_PtHard_4_10", "Pythia_Jets_PtHard_5_10",
        "Pythia_Jets_PtHard_6_10", "Pythia_Jets_PtHard_7_10", "Pythia_Jets_PtHard_8_10", "Pythia_Jets_PtHard_9_10", "Pythia_Jets_PtHard_10_10"
      };
      for (Int_t iHead = 0; iHead < 10; iHead++ ){
        TObjString *Header = new TObjString(nameHeaders[iHead]);
        HeaderList->Add(Header);
      }
    } else if (headerSelectionInt == 10){
      TObjString *Header1 = new TObjString("pythia_bele_10_10");
      HeaderList->Add(Header1);
    } else if (headerSelectionInt == 11){
      TString nameHeaders[34]   = { "Pythia_Jets_PtHard_1_10", "Pythia_Jets_PtHard_2_10", "Pythia_Jets_PtHard_3_10", "Pythia_Jets_PtHard_4_10", "Pythia_Jets_PtHard_5_10",
        "Pythia_Jets_PtHard_6_10", "Pythia_Jets_PtHard_7_10", "Pythia_Jets_PtHard_8_10", "Pythia_Jets_PtHard_9_10", "Pythia_Jets_PtHard_10_10",
        "gEMCPhoton_7", "flat pt kstar_8", "flat pt kstarbar_9", "pythia_cele_10_10", "pythia_cele_18_10",
        "pythia_cele_30_10", "pythia_cele_50_10", "pythia_ccbar_10_10", "pythia_ccbar_18_10", "pythia_ccbar_30_10",
        "pythia_ccbar_50_10", "pythia_bele_10_10", "pythia_bele_18_10", "pythia_bele_30_10", "pythia_bele_50_10",
        "pythia_bbbar_10_10", "pythia_bbbar_18_10", "pythia_bbbar_30_10", "pythia_bbbar_50_10", "pi0_1",
        "eta_2", "pi0EMC_3", "etaEMC_5", "hijing_0"
      };
      for (Int_t iHead = 0; iHead < 34; iHead++ ){
        TObjString *Header = new TObjString(nameHeaders[iHead]);
        HeaderList->Add(Header);
      }
    } else if (headerSelectionInt == 12){
      TObjString *Header1 = new TObjString("pi0PHS_4");
      HeaderList->Add(Header1);
    } else if (headerSelectionInt == 13){
      TObjString *Header1 = new TObjString("etaPHS_6");
      HeaderList->Add(Header1);
    } else if (headerSelectionInt == 14){
      TString nameHeaders[2]    = { "pi0PHS_4", "etaPHS_6" };
      for (Int_t iHead = 0; iHead < 2; iHead++ ){
        TObjString *Header = new TObjString(nameHeaders[iHead]);
        HeaderList->Add(Header);
      }
    } else {
      TString nameHeaders[2]    = { "pi0_1", "eta_2" };
      for (Int_t iHead = 0; iHead < 2; iHead++ ){
        TObjString *Header = new TObjString(nameHeaders[iHead]);
        HeaderList->Add(Header);
      }
    }
  } else if (generatorName.CompareTo("LHC14a1b")==0 || generatorName.CompareTo("LHC14a1c")==0){
    if (headerSelectionInt == 1 || headerSelectionInt == 2 || headerSelectionInt == 3 ){
      TObjString *Header1 = new TObjString("BOX");
      HeaderList->Add(Header1);
    } if (headerSelectionInt == 4 || headerSelectionInt == 5 || headerSelectionInt == 6 ){
      TObjString *Header1 = new TObjString("PARAM_EMC");
      HeaderList->Add(Header1);
    } if (headerSelectionInt == 12 || headerSelectionInt == 13 || headerSelectionInt == 14 ){
      TObjString *Header1 = new TObjString("PARAM_PHOS");
      HeaderList->Add(Header1);
    }
  } else if (generatorName.CompareTo("LHC16h4")==0 || generatorName.CompareTo("LHC19h3")==0){
    if (headerSelectionInt == 1){
      TObjString *Header1 = new TObjString("Injector (pi0)_1");
      HeaderList->Add(Header1);
    } else if (headerSelectionInt == 2){
      TObjString *Header1 = new TObjString("Injector (eta)_2");
      HeaderList->Add(Header1);
    }
  } else if ( (generatorName.CompareTo("LHC20g10")==0) || generatorName.BeginsWith("LHC24a1") ){

    auto fillSingle = [&HeaderList](Size_t theSwitch){
      if (theSwitch == 1 ) { HeaderList->Add(new TObjString("Injector (pi0)"));  } // PCM Pi0
      if (theSwitch == 2 ) { HeaderList->Add(new TObjString("Injector (pi0a)")); } // PHOS Pi0
      if (theSwitch == 3 ) { HeaderList->Add(new TObjString("Injector (pi0b)")); } // PHOS Pi0
      if (theSwitch == 4 ) { HeaderList->Add(new TObjString("Injector (pi0c)")); } // PHOS Pi0
      if (theSwitch == 5 ) { HeaderList->Add(new TObjString("Injector (pi0d)")); } // PHOS Pi0
      if (theSwitch == 6 ) { HeaderList->Add(new TObjString("Injector (eta)"));  } // PCM Eta
      if (theSwitch == 7 ) { HeaderList->Add(new TObjString("Injector (etaa)")); } // PHOS Eta
      if (theSwitch == 8 ) { HeaderList->Add(new TObjString("Hijing")); }          // HIJING
      if (theSwitch == 9 ) { HeaderList->Add(new TObjString("Pileup")); }          // PileUp
      if (theSwitch == 10 ) { HeaderList->Add(new TObjString("Injector (pi0e)")); } // EMCAL 1 Pi0
      if (theSwitch == 11 ) { HeaderList->Add(new TObjString("Injector (pi0f)")); } // EMCAL 2 Pi0
      if (theSwitch == 12 ) { HeaderList->Add(new TObjString("Injector (pi0g)")); } // EMCAL 3 Pi0
      if (theSwitch == 13 ) { HeaderList->Add(new TObjString("Injector (pi0h)")); } // DCAL 1 Pi0
      if (theSwitch == 14 ) { HeaderList->Add(new TObjString("Injector (pi0i)")); } // DCal 2 Pi0
      if (theSwitch == 15 ) { HeaderList->Add(new TObjString("Injector (etab)")); } // EMCAL 1 Eta
      if (theSwitch == 16 ) { HeaderList->Add(new TObjString("Injector (etac)")); } // EMCAL 2 Eta
      if (theSwitch == 17 ) { HeaderList->Add(new TObjString("Injector (etad)")); } // EMCAL 3 Eta
      if (theSwitch == 18 ) { HeaderList->Add(new TObjString("Injector (etae)")); } // DCAL 1 Eta
      if (theSwitch == 19 ) { HeaderList->Add(new TObjString("Injector (etaf)")); } // DCal 2 Eta
    };

    auto fillFromTo = [&fillSingle](Size_t theStart, Size_t theInclEnd){
      for (Size_t i=theStart; i<=theInclEnd; ++i) { fillSingle(i); }
    };

    if      (headerSelectionInt == 0)  { fillFromTo(1, 9); }            // all PCM + PHOS + HIJING + PileUp
    else if (headerSelectionInt <= 19)  { fillSingle(headerSelectionInt); } // a single one
    else if (headerSelectionInt == 20) { fillFromTo(1, 5); }            // pi0) + pi0x)
    else if (headerSelectionInt == 21) { fillFromTo(2, 5); }            // all pi0x)
    else if (headerSelectionInt == 22) { fillFromTo(6, 7); }            // eta) + etaa)
    else if (headerSelectionInt == 23) {                                // pi0) + eta)
      fillSingle(1);
      fillSingle(6);
    } 
    else if (headerSelectionInt == 24) { fillFromTo(10, 19); }            // all EMCal
    else if (headerSelectionInt == 25) { fillFromTo(1, 19); }             // all
    else if (headerSelectionInt == 26) {                                  // PCM + all EMCal Pi0 + Eta
      fillSingle(1);
      fillSingle(6);
      fillFromTo(10, 19);
    }
    else if (headerSelectionInt == 27) {                                  // PCM + all EMCal + HIJING Pi0 + Eta
      fillSingle(1);
      fillSingle(6);
      fillSingle(8);
      fillFromTo(10, 19);
    }
    else if (headerSelectionInt == 28) {                                  // PCM + PileUp Pi0 + Eta
      fillSingle(1);
      fillSingle(6);
      fillSingle(8);
      fillSingle(9);
    }
    else if (headerSelectionInt == 29) {                                  // PCM + all EMCal + HIJING + PileUp Pi0 + Eta
      fillSingle(1);
      fillSingle(6);
      fillFromTo(8, 19);
    }
    else if (headerSelectionInt == 30) {                                  // all EMCal Pi0
      fillFromTo(10, 14);
    }
    else if (headerSelectionInt == 31) {                                  // all EMCal + HIJING Pi0
      fillSingle(8);
      fillFromTo(10, 14);
    }
    else if (headerSelectionInt == 32) {                                  // all EMCal + HIJING + PileUp Pi0
      fillFromTo(8, 14);
    }
    else if (headerSelectionInt == 33) {                                  // PCM + all EMCal Pi0
      fillSingle(1);
      fillFromTo(10, 14);
    }
    else if (headerSelectionInt == 34) {                                  // PCM + all EMCal + HIJING Pi0
      fillSingle(1);
      fillSingle(8);
      fillFromTo(10, 14);
    }
    else if (headerSelectionInt == 35) {                                  // PCM + all EMCal + HIJING + PileUp Pi0
      fillSingle(1);
      fillFromTo(8, 14);
    }
    else if (headerSelectionInt == 36) {                                  // PCM + HIJING + PileUp Pi0
      fillSingle(1);
      fillSingle(8);
      fillSingle(9);
    }
    else if (headerSelectionInt == 37) {                                  // PCM + HIJING Pi0
      fillSingle(1);
      fillSingle(8);
      fillSingle(9);
    }
    else if (headerSelectionInt == 38) {                                  // all EMCal Eta
      fillFromTo(15, 19);
    }
    else if (headerSelectionInt == 39) {                                  // all EMCal + HIJING Eta
      fillSingle(8);
      fillFromTo(15, 19);
    }
    else if (headerSelectionInt == 40) {                                  // all EMCal + HIJING + PileUp Eta
      fillSingle(8);
      fillSingle(9);
      fillFromTo(15, 14);
    }
    else if (headerSelectionInt == 41) {                                  // PCM + all EMCal Eta
      fillSingle(6);
      fillFromTo(15, 19);
    }
    else if (headerSelectionInt == 42) {                                  // PCM + all EMCal + HIJING Eta
      fillSingle(6);
      fillSingle(8);
      fillFromTo(15, 19);
    }
    else if (headerSelectionInt == 43) {                                  // PCM + all EMCal + HIJING + PileUp Eta
      fillSingle(6);
      fillSingle(8);
      fillSingle(9);
      fillFromTo(15, 19);
    }
    else if (headerSelectionInt == 44) {                                  // PCM + HIJING + PileUp Eta
      fillSingle(6);
      fillSingle(8);
      fillSingle(9);
    }
    else if (headerSelectionInt == 45) {                                  // PCM + HIJING Eta
      fillSingle(6);
      fillSingle(8);
    }
    else {
      Error(Form("%s_%i", addTaskName.Data(),  trainConfig), "No valid header selection");
      return;
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
      if (enableLightOutput == 4) fTrackMatcher->SetLightOutput(kTRUE);
      mgr->AddTask(fTrackMatcher);
      mgr->ConnectInput(fTrackMatcher,0,cinput);
    }

    analysisEventCuts[i] = new AliConvEventCuts();
    // if ( trainConfig == 1){
    //   if (generatorName.CompareTo("LHC14a1a") ==0 || generatorName.CompareTo("LHC14a1b") ==0 || generatorName.CompareTo("LHC14a1c") ==0 ){
    //     if ( i == 0 && doWeightingPart)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNamePtWeights, Form("Pi0_Hijing_%s_PbPb_2760GeV_0005TPC",generatorName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_0005TPC",generatorName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0005V0M","Eta_Fit_Data_PbPb_2760GeV_0005V0M");
    //     if ( i == 1 && doWeightingPart)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNamePtWeights, Form("Pi0_Hijing_%s_PbPb_2760GeV_0510TPC",generatorName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_0510TPC",generatorName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0510V0M","Eta_Fit_Data_PbPb_2760GeV_0510V0M");
    //     if ( i == 2 && doWeightingPart)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNamePtWeights, Form("Pi0_Hijing_%s_PbPb_2760GeV_0010TPC",generatorName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_0010TPC",generatorName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0010V0M","Eta_Fit_Data_PbPb_2760GeV_0010V0M");
    //     if ( i == 3 && doWeightingPart)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNamePtWeights, Form("Pi0_Hijing_%s_PbPb_2760GeV_2040TPC",generatorName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_2040TPC",generatorName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_2040V0M","Eta_Fit_Data_PbPb_2760GeV_2040V0M");
    //     if ( i == 4 && doWeightingPart)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNamePtWeights, Form("Pi0_Hijing_%s_PbPb_2760GeV_2050TPC",generatorName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_2050TPC",generatorName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_2050V0M","Eta_Fit_Data_PbPb_2760GeV_2050V0M");
    //   }
    // }
    analysisEventCuts[i]->SetV0ReaderName(V0ReaderName);
    if (periodNameV0Reader.CompareTo("") != 0) analysisEventCuts[i]->SetPeriodEnum(periodNameV0Reader);

    if (trainConfig == 109 ){
      analysisEventCuts[i]->SelectSpecialTrigger(AliVEvent::kMB, "AliVEvent::kMB" );
    }
    if (trainConfig == 110 ){
      analysisEventCuts[i]->SelectSpecialTrigger(AliVEvent::kCentral,"AliVEvent::kCentral" );
    }
    if (trainConfig == 111 ){
      analysisEventCuts[i]->SelectSpecialTrigger(AliVEvent::kSemiCentral,"AliVEvent::kSemiCentral" );
    }
    
    
    if(periodNameAnchor.Contains("LHC11h") && enableFlattening){
      cout << "entering the cent. flattening loop -> searching for file: " << fileNameCentFlattening.Data() << endl;

      if( fileNameCentFlattening.Contains("FlatFile") ){
        analysisEventCuts[i]->SetUseWeightFlatCentralityFromFile(enableFlattening, fileNameCentFlattening, "Cent");
      } else if( fileNameCentFlattening.Contains("Good") ){
        analysisEventCuts[i]->SetUseWeightFlatCentralityFromFile(enableFlattening, fileNameCentFlattening, "CentGoodRuns");
      }else if( fileNameCentFlattening.Contains("SemiGood") ){
        analysisEventCuts[i]->SetUseWeightFlatCentralityFromFile(enableFlattening, fileNameCentFlattening, "CentSemiGoodRuns");
      }else {
        analysisEventCuts[i]->SetUseWeightFlatCentralityFromFile(enableFlattening, fileNameCentFlattening, "CentTotalRuns");
      }
    }

    TString histoNameMCPi0PT = "";  TString histoNameMCEtaPT = "";  TString histoNameMCK0sPT = "";    // pT spectra to be weighted
    TString fitNamePi0PT = "";      TString fitNameEtaPT = "";      TString fitNameK0sPT = "";        // fit to correct shape of pT spectra
    Bool_t weightPi0 = kFALSE;      Bool_t weightEta = kFALSE;      Bool_t weightK0s = kFALSE;
    if(doWeightingPart){
      cout << "INFO enabeling pT weighting" << endl;
      if(periodNameAnchor.CompareTo("LHC18q")==0){
        TString eventCutString  = cuts.GetEventCut(i);
        TString eventCutShort   = eventCutString(0,6);   // first six digits
        weightPi0         = kTRUE;
        histoNameMCPi0PT  = Form("Pi0_%s_5TeV_%s",   periodNameV0Reader.Data(), eventCutString.Data());  // MC
        fitNamePi0PT      = Form("Pi0_Data_5TeV_%s", eventCutShort.Data());                              // fit to data
        weightEta         = kTRUE;
        histoNameMCEtaPT  = Form("Eta_%s_5TeV_%s",   periodNameV0Reader.Data(), eventCutString.Data());
        fitNameEtaPT      = Form("Eta_Data_5TeV_%s", eventCutShort.Data());
      }
      analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(weightPi0, weightEta, weightK0s, fileNamePtWeights, histoNameMCPi0PT, histoNameMCEtaPT, histoNameMCK0sPT, fitNamePi0PT, fitNameEtaPT, fitNameK0sPT);
    }

    TString dataInputMultHisto  = "";
    TString mcInputMultHisto    = "";
    if (enableMultiplicityWeighting){
      cout << "INFO enableling mult weighting" << endl;
      if(periodNameAnchor.CompareTo("LHC15o")==0){
        TString cutNumber = cuts.GetEventCut(i);
        TString centCut = cutNumber(0,3);  // first three digits of event cut
        dataInputMultHisto = Form("%s_%s", periodNameAnchor.Data(), centCut.Data());
        mcInputMultHisto   = Form("%s_%s", generatorName.Data(), centCut.Data());
        cout << "INFO read " << dataInputMultHisto.Data() << " and " <<  mcInputMultHisto.Data() << " from " << fileNameMultWeights.Data() << endl;
      }
      analysisEventCuts[i]->SetUseWeightMultiplicityFromFile(kTRUE, fileNameMultWeights, dataInputMultHisto, mcInputMultHisto );
    }

    if (trainConfig == 34 || trainConfig == 35){
      if (generatorName.CompareTo("LHC14a1a") ==0 || generatorName.CompareTo("LHC14a1b") ==0 || generatorName.CompareTo("LHC14a1c") ==0 ){
        if ( doWeightingPart)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNamePtWeights, Form("Pi0_Hijing_%s_PbPb_2760GeV_0010TPC",generatorName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_0010TPC",generatorName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0010V0M","Eta_Fit_Data_PbPb_2760GeV_0010V0M");
      }
    }
    if (trainConfig == 36 || trainConfig == 37){
      if (generatorName.CompareTo("LHC14a1a") ==0 || generatorName.CompareTo("LHC14a1b") ==0 || generatorName.CompareTo("LHC14a1c") ==0 ){
        if ( doWeightingPart)   analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNamePtWeights, Form("Pi0_Hijing_%s_PbPb_2760GeV_2050TPC",generatorName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_2050TPC",generatorName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_2050V0M","Eta_Fit_Data_PbPb_2760GeV_2050V0M");
      }
    }

    if (trainConfig == 38 || trainConfig == 39 || trainConfig == 43 || trainConfig == 44){
      if (generatorName.CompareTo("LHC14a1a") ==0 || generatorName.CompareTo("LHC14a1b") ==0 || generatorName.CompareTo("LHC14a1c") ==0 ){
        if ( doWeightingPart)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNamePtWeights, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_0010TPC",generatorName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_0010TPC",generatorName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0010V0M","Eta_Fit_Data_PbPb_2760GeV_0010V0M");
      }
    }

    if (trainConfig == 40 || trainConfig == 41){
      if (generatorName.CompareTo("LHC14a1a") ==0 || generatorName.CompareTo("LHC14a1b") ==0 || generatorName.CompareTo("LHC14a1c") ==0 ){
        if ( doWeightingPart)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNamePtWeights, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_2050TPC",generatorName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_2050TPC",generatorName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_2050V0M","Eta_Fit_Data_PbPb_2760GeV_2050V0M");
      }
    }

    if (enableLightOutput == 1 || enableLightOutput == 2 ) analysisEventCuts[i]->SetLightOutput(1);
    if (enableLightOutput == 4) analysisEventCuts[i]->SetLightOutput(2);
    analysisEventCuts[i]->InitializeCutsFromCutString((cuts.GetEventCut(i)).Data());
    if (generatorName.CompareTo("LHC14a1b") ==0 || generatorName.CompareTo("LHC14a1c") ==0 ){
      if (headerSelectionInt == 1 || headerSelectionInt == 4 || headerSelectionInt == 12 ) analysisEventCuts[i]->SetAddedSignalPDGCode(111);
      if (headerSelectionInt == 2 || headerSelectionInt == 5 || headerSelectionInt == 13 ) analysisEventCuts[i]->SetAddedSignalPDGCode(221);
    }
    analysisEventCuts[i]->SetCorrectionTaskSetting(corrTaskSetting);
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
    EventCutList->Add(analysisEventCuts[i]);
    analysisEventCuts[i]->SetFillCutHistograms("",kFALSE);

    analysisCuts[i] = new AliConversionPhotonCuts();

    if (enableMatBudWeightsPi0 > 0){
      if (isMC > 0){
	      if (analysisCuts[i]->InitializeMaterialBudgetWeights(enableMatBudWeightsPi0,fileNameMatBudWeights)){
	        initializedMatBudWeigths_existing = kTRUE;
          } else {
            cout << "ERROR The initialization of the materialBudgetWeights did not work out." << endl;
          }
      } else {
        cout << "ERROR 'enableMatBudWeightsPi0'-flag was set > 0 even though this is not a MC task. It was automatically reset to 0." << endl;
        }
    }

    analysisCuts[i]->SetV0ReaderName(V0ReaderName);

    // extract single filenames from fileNamedEdxPostCalib
    TObjArray *lArrFnamesdEdxPostCalib = nullptr;
    if (enableElecDeDxPostCalibration == 2){
      if (isMC){
        cout << "ERROR enableElecDeDxPostCalibration set to 2 even if MC file. Automatically reset to 0"<< endl;
        enableElecDeDxPostCalibration = 0;
      } else {
        lArrFnamesdEdxPostCalib = fileNamedEdxPostCalib.Tokenize("+");
        if (!lArrFnamesdEdxPostCalib){
          cout << "ERROR fileNamedEdxPostCalib.Tokenize() returned nullptr\n";
          return;
        }
        int lNFnames = lArrFnamesdEdxPostCalib->GetEntriesFast();
        if (!(lNFnames == 1) && !(lNFnames == numberOfCuts)){
          cout << "ERROR: FEPC either has to be one filename or several separated by a '+' where the number of filenames has to match the number of cuts in the trainconfig.\nOr no filename at all if the file from the OADB is to be taken\n";
          return;
        }
      }
    }
    if (enableElecDeDxPostCalibration == 1){
      if (isMC == 0){
        if(fileNamedEdxPostCalib.CompareTo("") != 0){
          analysisCuts[i]->SetElecDeDxPostCalibrationCustomFile(fileNamedEdxPostCalib);
          cout << "Setting custom dEdx recalibration file: " << fileNamedEdxPostCalib.Data() << endl;
        }
        analysisCuts[i]->SetDoElecDeDxPostCalibration(kTRUE);
        cout << "Enabled TPC dEdx recalibration." << endl;
      } else{
        cout << "ERROR enableElecDeDxPostCalibration set to True even if MC file. Automatically reset to 0"<< endl;
        enableElecDeDxPostCalibration=kFALSE;
        analysisCuts[i]->SetDoElecDeDxPostCalibration(kFALSE);
      }
    } else if(enableElecDeDxPostCalibration == 2){
      int lNFnames = lArrFnamesdEdxPostCalib->GetEntriesFast();
      if (lNFnames){
        if (enableElecDeDxPostCalibration==2){
          analysisCuts[i]->ForceTPCRecalibrationAsFunctionOfConvR();
        }
        const TString &lFname = (static_cast<TObjString*>(lArrFnamesdEdxPostCalib->At(lNFnames > 1 ? i : 0)))->GetString();
        printf("Cut config %s_%s_%s:\nSetting custom dEdx recalibration file: %s\n", cuts.GetEventCut(i).Data(), cuts.GetPhotonCut(i).Data(), cuts.GetMesonCut(i).Data(), lFname.Data());
        Bool_t lSuccess = analysisCuts[i]->InitializeElecDeDxPostCalibration(lFname);
        if (!lSuccess) {
          return;
        }
      }
      analysisCuts[i]->SetDoElecDeDxPostCalibration(kTRUE);
      cout << "Enabled TPC dEdx recalibration." << endl;
    }
    if (enableLightOutput == 1 || enableLightOutput == 2 || enableLightOutput == 5 ) analysisCuts[i]->SetLightOutput(1);
    if (enableLightOutput == 4) analysisCuts[i]->SetLightOutput(2);
    if (enableLightOutput == 0) analysisCuts[i]->SetPlotTrackPID(kTRUE);
    analysisCuts[i]->InitializeCutsFromCutString((cuts.GetPhotonCut(i)).Data());
    ConvCutList->Add(analysisCuts[i]);
    analysisCuts[i]->SetFillCutHistograms("",kFALSE);

    analysisClusterCuts[i] = new AliCaloPhotonCuts(isMC);
    analysisClusterCuts[i]->SetHistoToModifyAcceptance(histoAcc);
    analysisClusterCuts[i]->SetV0ReaderName(V0ReaderName);
    analysisClusterCuts[i]->SetCorrectionTaskSetting(corrTaskSetting);
    analysisClusterCuts[i]->SetCaloTrackMatcherName(TrackMatcherName);
    if (enableLightOutput == 1 || enableLightOutput == 2 || enableLightOutput == 5 ) analysisClusterCuts[i]->SetLightOutput(1);
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

    if(analysisMesonCuts[i]->DoGammaSwappForBg()) analysisClusterCuts[i]->SetUseEtaPhiMapForBackCand(kTRUE);
    analysisClusterCuts[i]->SetFillCutHistograms("");

  }

  task->SetAllowOverlapHeaders(enableHeaderOverlap);
  task->SetEventCutList(numberOfCuts,EventCutList);
  task->SetConversionCutList(numberOfCuts,ConvCutList);
  task->SetCaloCutList(numberOfCuts,ClusterCutList);
  task->SetMesonCutList(numberOfCuts,MesonCutList);
  task->SetCorrectionTaskSetting(corrTaskSetting);
  task->SetMoveParticleAccordingToVertex(kTRUE);
  task->SetDoMesonAnalysis(kTRUE);
  task->SetDoMesonQA(enableQAMesonTask); //Attention new switch for Pi0 QA
  task->SetDoPhotonQA(enableQAPhotonTask);  //Attention new switch small for Photon QA
  if (enableLightOutput != 4) task->SetDoClusterQA(1);  //Attention new switch small for Cluster QA
  task->SetEnableSortingOfMCClusLabels(enableSortingMCLabels);
  task->SetDoTreeInvMassShowerShape(doTreeClusterShowerShape);
  task->SetUseTHnSparse(enableTHnSparse);
  if(enableExtMatchAndQA > 1){ task->SetPlotHistsExtQA(kTRUE);}

  if (initializedMatBudWeigths_existing) {
      task->SetDoMaterialBudgetWeightingOfGammasForTrueMesons(kTRUE);
  }


  //connect containers
  AliAnalysisDataContainer *coutput =
    mgr->CreateContainer(!(corrTaskSetting.CompareTo("")) ? Form("GammaConvCalo_%i",trainConfig) : Form("GammaConvCalo_%i_%s",trainConfig,corrTaskSetting.Data()), TList::Class(),
              AliAnalysisManager::kOutputContainer,Form("GCoCa_%i.root",trainConfig));

  mgr->AddTask(task);
  mgr->ConnectInput(task,0,cinput);
  mgr->ConnectOutput(task,1,coutput);
  Int_t nContainer = 2;
  for(Int_t i = 0; i<numberOfCuts; i++){
    if(enableQAPhotonTask>1){
      mgr->ConnectOutput(task,nContainer,mgr->CreateContainer(Form("%s_%s_%s_%s Photon DCA tree",(cuts.GetEventCut(i)).Data(),(cuts.GetPhotonCut(i)).Data(),(cuts.GetClusterCut(i)).Data(),(cuts.GetMesonCut(i)).Data()), TTree::Class(), AliAnalysisManager::kOutputContainer, Form("GCoCa_%i.root",trainConfig)) );
      nContainer++;
    }
    if(enableQAMesonTask>1){
	    mgr->ConnectOutput(task,nContainer,mgr->CreateContainer(Form("%s_%s_%s_%s Meson DCA tree",(cuts.GetEventCut(i)).Data(),(cuts.GetPhotonCut(i)).Data(),(cuts.GetClusterCut(i)).Data(),(cuts.GetMesonCut(i)).Data()), TTree::Class(), AliAnalysisManager::kOutputContainer, Form("GCoCa_%i.root",trainConfig)) );
      nContainer++;
    }
  }
  return;

}
