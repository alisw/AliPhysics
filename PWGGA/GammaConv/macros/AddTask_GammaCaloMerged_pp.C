/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: Friederike Bock, Daniel Mühlheim                               *
 * Version 1.0                                                            *
 *                                                                        *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

//***************************************************************************************
//This AddTask is supposed to set up the main task
//($ALIPHYSICS/PWGGA/GammaConv/AliAnalysisTaskGammaCaloMerged.cxx) for
//pp together with all supporting classes
//***************************************************************************************

//***************************************************************************************
//main function
//***************************************************************************************
void AddTask_GammaCaloMerged_pp(
  Int_t     trainConfig                   = 1,        // change different set of cuts
  Int_t     isMC                          = 0,        // run MC
  TString   photonCutNumberV0Reader       = "",       // 00000008400000000100000000 nom. B, 00000088400000000100000000 low B
  TString   periodNameV0Reader            = "",
  // general setting for task
  Int_t     enableQAMesonTask             = 0,        // enable QA in AliAnalysisTaskGammaConvV1
  Int_t     enableQAClusterTask           = 0,        // enable additional QA task
  Int_t     enableExtMatchAndQA           = 0,        // disabled (0), extMatch (1), extQA_noCellQA (2), extMatch+extQA_noCellQA (3), extQA+cellQA (4), extMatch+extQA+cellQA (5)
  Int_t     enableLightOutput             = 0,        // switch to run light output (only essential histograms for afterburner)
  Int_t     enableTriggerMimicking        = 0,        // enable trigger mimicking
  Bool_t    enableTriggerOverlapRej       = kFALSE,   // enable trigger overlap rejection
  TString   settingMaxFacPtHard           = "3.",     // maximum factor between hardest jet and ptHard generated
  // settings for weights
  // FPTW:fileNamePtWeights, separate with ;
  TString   fileNameExternalInputs        = "",
  Bool_t    doWeightingPart               = kFALSE,   // enable Weighting
  TString   generatorName                 = "",       // generator name
  // special settings
  Int_t     selectedMeson                 = 1,        // put flag for selected meson
  Bool_t    enableSortingMCLabels         = kTRUE,    // enable sorting for MC cluster labels
  Bool_t    enableDetailedPrintout        = kFALSE,   // enable detailed printout
  Double_t  minEnergyForExoticsCut        = 1.0,      // minimum energy to be used for exotics CutHandler
  Bool_t    enableExoticsQA               = kFALSE,   // switch to run QA for exotic clusters
  Bool_t    runDetailedM02                = kFALSE,   // switch on very detailed M02 distribution
  Int_t     minAllowedPi0Overlaps         = -1,   // set maximum number of Pi0 overlaps in MC
  Int_t     maxAllowedPi0Overlaps         = -1,   // set maximum number of Pi0 overlaps in MC
  Bool_t    doOverlapsFromCluster         = kFALSE, // overlaps as event criteria (false) or from single clusters (true)
  // subwagon config
  TString   additionalTrainConfig         = "0"       // additional counter for trainconfig
  ) {

  AliCutHandlerPCM cuts;

  TString fileNamePtWeights     = cuts.GetSpecialFileNameFromString (fileNameExternalInputs, "FPTW:");
  TString fileNameCustomTriggerMimicOADB   = cuts.GetSpecialFileNameFromString (fileNameExternalInputs, "FTRM:");

  TString addTaskName                 = "AddTask_GammaMerged_pp";
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

  TObjArray *rmaxFacPtHardSetting = settingMaxFacPtHard.Tokenize("_");
  if(rmaxFacPtHardSetting->GetEntries()<1){cout << "ERROR: AddTask_GammaCaloMerged_pp during parsing of settingMaxFacPtHard String '" << settingMaxFacPtHard.Data() << "'" << endl; return;}
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
  AliAnalysisManager *mgr           = AliAnalysisManager::GetAnalysisManager();
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
  AliAnalysisTaskGammaCaloMerged *task  = NULL;
  task                                  = new AliAnalysisTaskGammaCaloMerged(Form("GammaCaloMerged_%i",trainConfig));
  task->SetIsHeavyIon(isHeavyIon);
  task->SetIsMC(isMC);
  task->SetV0ReaderName(V0ReaderName);
  task->SetLightOutput(enableLightOutput);
  task->SetCorrectionTaskSetting(corrTaskSetting);
  task->SetTrackMatcherRunningMode(trackMatcherRunningMode);

  // cluster cuts
  // 0 "ClusterType",  1 "EtaMin", 2 "EtaMax", 3 "PhiMin", 4 "PhiMax", 5 "DistanceToBadChannel", 6 "Timing", 7 "TrackMatching", 8 "ExoticCell",
  // 9 "MinEnergy", 10 "MinNCells", 11 "MinM02", 12 "MaxM02", 13 "MinM20", 14 "MaxM20", 15 "MaximumDispersion", 16 "NLM"

  // ************************************* EMCAL cuts ****************************************************
  // LHC11a
  if (trainConfig == 1){ // all defaults for LHC11a
    cuts.AddCutMergedCalo("00003113","1111121053032200000","1111121053022210001","0163301100000000"); // INT1
    cuts.AddCutMergedCalo("00051013","1111121053032200000","1111121053022210001","0163301100000000"); // EMC1
  } else if (trainConfig == 2){ // no TM in basis cut
    cuts.AddCutMergedCalo("00003113","1111121050032200000","1111121053022210001","0163301100000000"); // INT1
    cuts.AddCutMergedCalo("00051013","1111121050032200000","1111121053022210001","0163301100000000"); // EMC1
  } else if (trainConfig == 3){ // open M02, open Mass, open Alpha
    cuts.AddCutMergedCalo("00003113","1111121053032200000","1111121053022000001","0163300000000000"); // INT1
    cuts.AddCutMergedCalo("00051013","1111121053032200000","1111121053022000001","0163300000000000"); // EMC1
  } else if (trainConfig == 4){  //NLM exotics default frac = 0.97
    cuts.AddCutMergedCalo("00003113","1111121053532200000","1111121053522110001","0163301100000000"); // INT1
    cuts.AddCutMergedCalo("00051013","1111121053532200000","1111121053522110001","0163301100000000"); // EMC1
  } else if (trainConfig == 5){  // new default
    cuts.AddCutMergedCalo("00003113","1111121053032200000","1111121053022700001","0163300700000000"); // Mass only band at 0, no explicit exotics cut, M02 cut at 0.27
    cuts.AddCutMergedCalo("00051013","1111121053032200000","1111121053022700001","0163300700000000"); // Mass only band at 0, no explicit exotics cut, M02 cut at 0.27
    cuts.AddCutMergedCalo("00003113","1111121053032200000","1111121053022700001","0163300000000000"); // INT1
    cuts.AddCutMergedCalo("00051013","1111121053032200000","1111121053022700001","0163300000000000"); // EMC1
  } else if (trainConfig == 6){  // new default
    cuts.AddCutMergedCalo("00003113","1111121050032200000","1111121053022000001","0163300000000000"); // no M02, no exotics, TM only for second cut
    cuts.AddCutMergedCalo("00051013","1111121050032200000","1111121053022000001","0163300000000000"); // no M02, no exotics, TM only for second cut
    cuts.AddCutMergedCalo("00003113","1111121050032200000","1111121053022700001","0163300000000000"); // M02 > 0.27, no exotics, TM only for second cut
    cuts.AddCutMergedCalo("00051013","1111121050032200000","1111121053022700001","0163300000000000"); // M02 > 0.27, no exotics, TM only for second cut
  } else if (trainConfig == 7){  // new default , pt dep TM
    cuts.AddCutMergedCalo("00003113","1111121057032200000","1111121057022700001","0163300000000000"); // M02 > 0.27, no exotics
    cuts.AddCutMergedCalo("00051013","1111121057032200000","1111121057022700001","0163300000000000"); // M02 > 0.27, no exotics
    cuts.AddCutMergedCalo("00003113","1111121050032200000","1111121057022700001","0163300000000000"); // M02 > 0.27, no exotics, TM only for second
    cuts.AddCutMergedCalo("00051013","1111121050032200000","1111121057022700001","0163300000000000"); // M02 > 0.27, no exotics, TM only for second
  } else if (trainConfig == 8){  // new default, with eta < 0.7, y < 0.7
    cuts.AddCutMergedCalo("00003113","1551121053032200000","1551121053022700001","0163200000000000"); // Mass only band at 0, no explicit exotics cut, M02 cut at 0.27
    cuts.AddCutMergedCalo("00051013","1551121053032200000","1551121053022700001","0163200000000000"); // Mass only band at 0, no explicit exotics cut, M02 cut at 0.27
    cuts.AddCutMergedCalo("00003113","1551121057032200000","1551121057022700001","0163200000000000"); // Mass only band at 0, no explicit exotics cut, M02 cut at 0.27, pt dep TM
    cuts.AddCutMergedCalo("00051013","1551121057032200000","1551121057022700001","0163200000000000"); // Mass only band at 0, no explicit exotics cut, M02 cut at 0.27, pt dep TM
  } else if (trainConfig == 9){  // new default , pt dep TM, no M02
    cuts.AddCutMergedCalo("00003113","1111121057032200000","1111121057022000001","0163300000000000"); // no M02, no exotics
    cuts.AddCutMergedCalo("00051013","1111121057032200000","1111121057022000001","0163300000000000"); // no M02, no exotics

  // INT1 variations
  } else if (trainConfig == 10){ // M02 var
    cuts.AddCutMergedCalo("00003113","1111121053032200000","1111121053022110001","0163301100000000"); // min 0.3 function default
    cuts.AddCutMergedCalo("00003113","1111121053032200000","1111121053022310001","0163301100000000"); // min 0.25 function default
    cuts.AddCutMergedCalo("00003113","1111121053032200000","1111121053022410001","0163301100000000"); // min 0.27, tighter lower func
    cuts.AddCutMergedCalo("00003113","1111121053032200000","1111121053022520001","0163301100000000"); // min 0.27, looser func up and low
    cuts.AddCutMergedCalo("00003113","1111121053032200000","1111121053022430001","0163301100000000"); // min 0.27, tighter func up and low
  } else if (trainConfig == 11){  // variation track matching to cluster & mass variations
    cuts.AddCutMergedCalo("00003113","1111121050032200000","1111121050022210001","0163301100000000"); // no TM
    cuts.AddCutMergedCalo("00003113","1111121051032200000","1111121051022210001","0163301100000000"); // looser TM
    cuts.AddCutMergedCalo("00003113","1111121056032200000","1111121056022210001","0163301100000000"); // tighter TM
    cuts.AddCutMergedCalo("00003113","1111121053032200000","1111121053022210001","0163301300000000"); // tighter mass cut
    cuts.AddCutMergedCalo("00003113","1111121053032200000","1111121053022210001","0163301500000000"); // looser mass cut
  } else if (trainConfig == 12){ // NL var
    cuts.AddCutMergedCalo("00003113","1111122053032200000","1111122053022210001","0163301100000000"); // SDM loose time
    cuts.AddCutMergedCalo("00003113","1111111053032200000","1111111053022210001","0163301100000000"); // conv calo tight time
    cuts.AddCutMergedCalo("00003113","1111101053032200000","1111101053022210001","0163301100000000"); // SDM Jason
    cuts.AddCutMergedCalo("00003113","1111100053032200000","1111100053022210001","0163301100000000"); // none
  } else if (trainConfig == 13){ // Alpha cut variations & TRD material
    cuts.AddCutMergedCalo("00003113","1111121053032200000","1111121053022210001","0163303100000000"); // NLM 1 looser
    cuts.AddCutMergedCalo("00003113","1111121053032200000","1111121053022210001","0163305100000000"); // NLM 1 tighter
    cuts.AddCutMergedCalo("00003113","1113121053032200000","1113121053022210001","0163301100000000");// TRD infront
    cuts.AddCutMergedCalo("00003113","1111221053032200000","1111221053022210001","0163301100000000");// no TRD infront
  } else if (trainConfig == 14){ // varied exoctics
    cuts.AddCutMergedCalo("00003113","1111121053032200000","1111121053022700001","0163300000000000"); // no exotics
    cuts.AddCutMergedCalo("00003113","1111121053232200000","1111121053222700001","0163300000000000"); // frac = 0.99
    cuts.AddCutMergedCalo("00003113","1111121053332200000","1111121053322700001","0163300000000000"); // frac = 0.98
    cuts.AddCutMergedCalo("00003113","1111121053532200000","1111121053522700001","0163300000000000"); // frac = 0.97
    cuts.AddCutMergedCalo("00003113","1111121053732200000","1111121053722700001","0163300000000000"); // frac = 0.96
    cuts.AddCutMergedCalo("00003113","1111121053932200000","1111121053922700001","0163300000000000"); // frac = 0.94
  } else if (trainConfig == 15){  // variation of mass and alpha cut
    cuts.AddCutMergedCalo("00003113","1111121053032200000","1111121053022700001","0163301700000000"); // only band at 0
    cuts.AddCutMergedCalo("00003113","1111121053032200000","1111121053022700001","0163300700000000"); // only band at 0, no alpha
    cuts.AddCutMergedCalo("00003113","1111121053032200000","1111121053022700001","0163301000000000"); // no mass cut
    cuts.AddCutMergedCalo("00003113","1111121053032200000","1111121053022700001","0163300000000000"); // no mass cut, no alpha
    cuts.AddCutMergedCalo("00003113","1111121053032200000","1111121053022700001","0163300100000000"); // no alpha
  } else if (trainConfig == 16){ // varied M02
    cuts.AddCutMergedCalo("00003113","1111121053032200000","1111121053022600001","0163300700000000"); // min M02= 0.3
    cuts.AddCutMergedCalo("00003113","1111121053032200000","1111121053022800001","0163300700000000"); // min M02= 0.25
    cuts.AddCutMergedCalo("00003113","1111121053032200000","1111121053022600001","0163300000000000"); // no mass, min M02 = 0.3
    cuts.AddCutMergedCalo("00003113","1111121053032200000","1111121053022700001","0163300000000000"); // no mass, min M02 = 0.27
    cuts.AddCutMergedCalo("00003113","1111121053032200000","1111121053022800001","0163300000000000"); // no mass, min M02 = 0.25
  } else if (trainConfig == 17){  // variation track matching to cluster & mass variations new defaults
    cuts.AddCutMergedCalo("00003113","1111121050032200000","1111121050022700001","0163300000000000"); // no TM
    cuts.AddCutMergedCalo("00003113","1111121051032200000","1111121051022700001","0163300000000000"); // looser TM
    cuts.AddCutMergedCalo("00003113","1111121056032200000","1111121056022700001","0163300000000000"); // tighter TM
    cuts.AddCutMergedCalo("00003113","1113121053032200000","1113121053022700001","0163300000000000");// TRD infront
    cuts.AddCutMergedCalo("00003113","1111221053032200000","1111221053022700001","0163300000000000");// no TRD infront
  } else if (trainConfig == 18){ // NL var new defaults
    cuts.AddCutMergedCalo("00003113","1111122053032200000","1111122053022700001","0163300000000000"); // SDM loose time
    cuts.AddCutMergedCalo("00003113","1111111053032200000","1111111053022700001","0163300000000000"); // conv calo tight time
    cuts.AddCutMergedCalo("00003113","1111101053032200000","1111101053022700001","0163300000000000"); // SDM Jason
    cuts.AddCutMergedCalo("00003113","1111100053032200000","1111100053022700001","0163300000000000"); // none

  // EMC1 variations
  } else if (trainConfig == 20){ // M02 var
    cuts.AddCutMergedCalo("00051013","1111121053032200000","1111121053022110001","0163301100000000"); // min 0.3 function default
    cuts.AddCutMergedCalo("00051013","1111121053032200000","1111121053022310001","0163301100000000"); // min 0.25 function default
    cuts.AddCutMergedCalo("00051013","1111121053032200000","1111121053022410001","0163301100000000"); // min 0.27, tighter lower func
    cuts.AddCutMergedCalo("00051013","1111121053032200000","1111121053022520001","0163301100000000"); // min 0.27, looser func up and low
    cuts.AddCutMergedCalo("00051013","1111121053032200000","1111121053022430001","0163301100000000"); // min 0.27, tighter func up and low
  } else if (trainConfig == 21){ // variation track matching to cluster & mass variations
    cuts.AddCutMergedCalo("00051013","1111121050032200000","1111121050022210001","0163301100000000"); // no TM
    cuts.AddCutMergedCalo("00051013","1111121051032200000","1111121051022210001","0163301100000000"); // looser TM
    cuts.AddCutMergedCalo("00051013","1111121056032200000","1111121056022210001","0163301100000000"); // tighter TM
    cuts.AddCutMergedCalo("00051013","1111121053032200000","1111121053022210001","0163301300000000"); // tighter mass cut
    cuts.AddCutMergedCalo("00051013","1111121053032200000","1111121053022210001","0163301500000000"); // looser mass cut
  } else if (trainConfig == 22){ // NL var
    cuts.AddCutMergedCalo("00051013","1111122053032200000","1111122053022210001","0163301100000000"); // SDM loose time
    cuts.AddCutMergedCalo("00051013","1111111053032200000","1111111053022210001","0163301100000000"); // conv calo tight time
    cuts.AddCutMergedCalo("00051013","1111101053032200000","1111101053022210001","0163301100000000"); // SDM Jason
    cuts.AddCutMergedCalo("00051013","1111100053032200000","1111100053022210001","0163301100000000"); // none
  } else if (trainConfig == 23){ // Alpha cut variations & TRD material
    cuts.AddCutMergedCalo("00051013","1111121053032200000","1111121053022210001","0163303100000000"); // NLM 1 looser
    cuts.AddCutMergedCalo("00051013","1111121053032200000","1111121053022210001","0163305100000000"); // NLM 1 tighter
    cuts.AddCutMergedCalo("00051013","1113121053032200000","1113121053022210001","0163301100000000");// TRD infront
    cuts.AddCutMergedCalo("00051013","1111221053032200000","1111221053022210001","0163301100000000");// no TRD infront
  } else if (trainConfig == 24){ // varied exoctics
    cuts.AddCutMergedCalo("00051013","1111121053032200000","1111121053022700001","0163300000000000"); // no exotics
    cuts.AddCutMergedCalo("00051013","1111121053232200000","1111121053222700001","0163300000000000"); // frac = 0.99
    cuts.AddCutMergedCalo("00051013","1111121053332200000","1111121053322700001","0163300000000000"); // frac = 0.98
    cuts.AddCutMergedCalo("00051013","1111121053532200000","1111121053522700001","0163300000000000"); // frac = 0.97
    cuts.AddCutMergedCalo("00051013","1111121053732200000","1111121053722700001","0163300000000000"); // frac = 0.96
    cuts.AddCutMergedCalo("00051013","1111121053932200000","1111121053922700001","0163300000000000"); // frac = 0.94
  } else if (trainConfig == 25){  // variation of mass and alpha cut
    cuts.AddCutMergedCalo("00051013","1111121053032200000","1111121053022700001","0163301700000000"); // only band at 0
    cuts.AddCutMergedCalo("00051013","1111121053032200000","1111121053022700001","0163300700000000"); // only band at 0, no alpha
    cuts.AddCutMergedCalo("00051013","1111121053032200000","1111121053022700001","0163301000000000"); // no mass cut
    cuts.AddCutMergedCalo("00051013","1111121053032200000","1111121053022700001","0163300000000000"); // no mass cut, no alpha
    cuts.AddCutMergedCalo("00051013","1111121053032200000","1111121053022700001","0163300100000000"); // no alpha
  } else if (trainConfig == 26){ // varied M02
    cuts.AddCutMergedCalo("00051013","1111121053032200000","1111121053022600001","0163300700000000"); // min M02= 0.3
    cuts.AddCutMergedCalo("00051013","1111121053032200000","1111121053022800001","0163300700000000"); // min M02= 0.25
    cuts.AddCutMergedCalo("00051013","1111121053032200000","1111121053022600001","0163300000000000"); // no mass, min M02 = 0.3
    cuts.AddCutMergedCalo("00051013","1111121053032200000","1111121053022700001","0163300000000000"); // no mass, min M02 = 0.27
    cuts.AddCutMergedCalo("00051013","1111121053032200000","1111121053022800001","0163300000000000"); // no mass, min M02 = 0.25
  } else if (trainConfig == 27){  // variation track matching to cluster & mass variations new defaults
    cuts.AddCutMergedCalo("00051013","1111121050032200000","1111121050022700001","0163300000000000"); // no TM
    cuts.AddCutMergedCalo("00051013","1111121051032200000","1111121051022700001","0163300000000000"); // looser TM
    cuts.AddCutMergedCalo("00051013","1111121056032200000","1111121056022700001","0163300000000000"); // tighter TM
    cuts.AddCutMergedCalo("00051013","1113121053032200000","1113121053022700001","0163300000000000");// TRD infront
    cuts.AddCutMergedCalo("00051013","1111221053032200000","1111221053022700001","0163300000000000");// no TRD infront
  } else if (trainConfig == 28){ // NL var new defaults
    cuts.AddCutMergedCalo("00051013","1111122053032200000","1111122053022700001","0163300000000000"); // SDM loose time
    cuts.AddCutMergedCalo("00051013","1111111053032200000","1111111053022700001","0163300000000000"); // conv calo tight time
    cuts.AddCutMergedCalo("00051013","1111101053032200000","1111101053022700001","0163300000000000"); // SDM Jason
    cuts.AddCutMergedCalo("00051013","1111100053032200000","1111100053022700001","0163300000000000"); // none

  } else if (trainConfig == 38){     // V1 clusterizer no NLM
    cuts.AddCutMergedCalo("00003113","1111121053032200000","1111121053022000000","0163300000000000"); // INT1
    cuts.AddCutMergedCalo("00003113","1111121053032200000","1111121053022700000","0163300000000000"); // INT1
    cuts.AddCutMergedCalo("00051013","1111121053032200000","1111121053022000000","0163300000000000"); // EMC1
    cuts.AddCutMergedCalo("00051013","1111121053032200000","1111121053022700000","0163300000000000"); // EMC1
  } else if (trainConfig == 39){     // V1 clusterizer NLM2
    cuts.AddCutMergedCalo("00003113","1111121053032200000","1111121053022000002","0163300000000000"); // INT1
    cuts.AddCutMergedCalo("00003113","1111121053032200000","1111121053022700002","0163300700000000"); // INT1
    cuts.AddCutMergedCalo("00003113","1111121053032200000","1111121053022210002","0163302200000000"); // INT1
    cuts.AddCutMergedCalo("00051013","1111121053032200000","1111121053022000002","0163300000000000"); // EMC1
    cuts.AddCutMergedCalo("00051013","1111121053032200000","1111121053022700002","0163300700000000"); // EMC1
    cuts.AddCutMergedCalo("00051013","1111121053032200000","1111121053022210002","0163302200000000"); // EMC1

  // LHC13g
  } else if (trainConfig == 40){  // NLM1 with mass cuts and alpha, TM fixed vs pt
    cuts.AddCutMergedCalo("00010113","1111121063032200000","1111121063022210001","0163301100000000"); // INT7
    cuts.AddCutMergedCalo("00052013","1111121063032200000","1111121063022210001","0163301100000000"); // EMC7
    cuts.AddCutMergedCalo("00085013","1111121063032200000","1111121063022210001","0163301100000000"); // EG2
    cuts.AddCutMergedCalo("00083013","1111121063032200000","1111121063022210001","0163301100000000"); // EG1
  } else if (trainConfig == 41){   // NLM1 with mass cuts and alpha, TM fixed vs pt, no TM in basis cut
    cuts.AddCutMergedCalo("00010113","1111121060032200000","1111121063022210001","0163301100000000"); // INT7
    cuts.AddCutMergedCalo("00052013","1111121060032200000","1111121063022210001","0163301100000000"); // EMC7
    cuts.AddCutMergedCalo("00085013","1111121060032200000","1111121063022210001","0163301100000000"); // EG2
    cuts.AddCutMergedCalo("00083013","1111121060032200000","1111121063022210001","0163301100000000"); // EG1
  } else if (trainConfig == 42){  // NLM1 no M02, no mass, no alpha, TM fixed vs pt
    cuts.AddCutMergedCalo("00010113","1111121063032200000","1111121063022000001","0163300000000000"); // INT7
    cuts.AddCutMergedCalo("00052013","1111121063032200000","1111121063022000001","0163300000000000"); // EMC7
    cuts.AddCutMergedCalo("00085013","1111121063032200000","1111121063022000001","0163300000000000"); // EG2
    cuts.AddCutMergedCalo("00083013","1111121063032200000","1111121063022000001","0163300000000000"); // EG1
  } else if (trainConfig == 43){  // NLM1 Mass only band at 0, no explicit exotics cut, M02 cut at 0.27, TM fixed vs pt
    cuts.AddCutMergedCalo("00010113","1111121063032200000","1111121063022700001","0163300700000000"); // INT7
    cuts.AddCutMergedCalo("00052013","1111121063032200000","1111121063022700001","0163300700000000"); // EMC7
    cuts.AddCutMergedCalo("00085013","1111121063032200000","1111121063022700001","0163300700000000"); // EG2
    cuts.AddCutMergedCalo("00083013","1111121063032200000","1111121063022700001","0163300700000000"); // EG1
  } else if (trainConfig == 44){  // NLM1 no explicit exotics cut, M02 cut at 0.27, TM fixed vs pt
    cuts.AddCutMergedCalo("00010113","1111121063032200000","1111121063022700001","0163300000000000"); // INT7
    cuts.AddCutMergedCalo("00052013","1111121063032200000","1111121063022700001","0163300000000000"); // EMC7
    cuts.AddCutMergedCalo("00085013","1111121063032200000","1111121063022700001","0163300000000000"); // EG2
    cuts.AddCutMergedCalo("00083013","1111121063032200000","1111121063022700001","0163300000000000"); // EG1
  } else if (trainConfig == 45){  // NLM1 exotics default frac = 0.97
    cuts.AddCutMergedCalo("00010113","1111121063532200000","1111121063522110001","0163301100000000"); // INT7
    cuts.AddCutMergedCalo("00052013","1111121063532200000","1111121063522110001","0163301100000000"); // EMC7
    cuts.AddCutMergedCalo("00085013","1111121063532200000","1111121063522110001","0163301100000000"); // EG2
    cuts.AddCutMergedCalo("00083013","1111121063532200000","1111121063522110001","0163301100000000"); // EG1
  } else if (trainConfig == 46){  // NLM1 no explicit exotics cut, M02 cut at 0.27, TM fixed vs pt, TM only in merged
    cuts.AddCutMergedCalo("00010113","1111121060032200000","1111121063022700001","0163300000000000"); // INT7
    cuts.AddCutMergedCalo("00052013","1111121060032200000","1111121063022700001","0163300000000000"); // EMC7
    cuts.AddCutMergedCalo("00085013","1111121060032200000","1111121063022700001","0163300000000000"); // EG2
    cuts.AddCutMergedCalo("00083013","1111121060032200000","1111121063022700001","0163300000000000"); // EG1
  } else if (trainConfig == 47){  // NLM1 no explicit exotics cut, no M02 cut, TM fixed vs pt
    cuts.AddCutMergedCalo("00010113","1111121060032200000","1111121063022000001","0163300000000000"); // INT7
    cuts.AddCutMergedCalo("00052013","1111121060032200000","1111121063022000001","0163300000000000"); // EMC7
    cuts.AddCutMergedCalo("00085013","1111121060032200000","1111121063022000001","0163300000000000"); // EG2
    cuts.AddCutMergedCalo("00083013","1111121060032200000","1111121063022000001","0163300000000000"); // EG1
  } else if (trainConfig == 48){  // pp 2.76TeV paper cuts
    cuts.AddCutMergedCalo("00010113","1111121067032200000","1111121067022700001","0163300000000000"); // INT7
    cuts.AddCutMergedCalo("00052013","1111121067032200000","1111121067022700001","0163300000000000"); // EMC7
    cuts.AddCutMergedCalo("00085013","1111121067032200000","1111121067022700001","0163300000000000"); // EG2
    cuts.AddCutMergedCalo("00083013","1111121067032200000","1111121067022700001","0163300000000000"); // EG1
  } else if (trainConfig == 49){  // pp 2.76TeV paper cuts  w/o mass
    cuts.AddCutMergedCalo("00010113","1111121067032200000","1111121067022000001","0163300000000000");
    cuts.AddCutMergedCalo("00052013","1111121067032200000","1111121067022000001","0163300000000000"); // EMC7
    cuts.AddCutMergedCalo("00085013","1111121067032200000","1111121067022000001","0163300000000000"); // EG2
    cuts.AddCutMergedCalo("00083013","1111121067032200000","1111121067022000001","0163300000000000"); // EG1

  // INT7 variations
  } else if (trainConfig == 50){ // NLM 1 M02 var
    cuts.AddCutMergedCalo("00010113","1111121063032200000","1111121063022110001","0163301100000000"); // min 0.3 function default
    cuts.AddCutMergedCalo("00010113","1111121063032200000","1111121063022310001","0163301100000000"); // min 0.25 function default
    cuts.AddCutMergedCalo("00010113","1111121063032200000","1111121063022410001","0163301100000000"); // min 0.27, tighter lower func
    cuts.AddCutMergedCalo("00010113","1111121063032200000","1111121063022520001","0163301100000000"); // min 0.27, looser func up and low
    cuts.AddCutMergedCalo("00010113","1111121063032200000","1111121063022430001","0163301100000000"); // min 0.27, tighter func up and low
  } else if (trainConfig == 51){  // EMCAL clusters, variation track matching to cluster & Mass INT7 NLM 1
    cuts.AddCutMergedCalo("00010113","1111121060032200000","1111121060022210001","0163301100000000"); // no TM
    cuts.AddCutMergedCalo("00010113","1111121061032200000","1111121061022210001","0163301100000000"); // looser TM
    cuts.AddCutMergedCalo("00010113","1111121066032200000","1111121066022210001","0163301100000000"); // tighter TM
    cuts.AddCutMergedCalo("00010113","1111121063032200000","1111121063022210001","0163301300000000"); // tighter mass cut
    cuts.AddCutMergedCalo("00010113","1111121063032200000","1111121063022210001","0163301500000000"); // looser mass cut
  } else if (trainConfig == 52){  // NL variations INT7
    cuts.AddCutMergedCalo("00010113","1111122063032200000","1111122063022210001","0163301100000000"); // SDM loose time
    cuts.AddCutMergedCalo("00010113","1111111063032200000","1111111063022210001","0163301100000000"); // conv calo tight time
    cuts.AddCutMergedCalo("00010113","1111101063032200000","1111101063022210001","0163301100000000"); // SDM Jason
    cuts.AddCutMergedCalo("00010113","1111100063032200000","1111100063022210001","0163301100000000"); // none
  } else if (trainConfig == 53){  // Alpha cut variations & TRD material INT7 NLM 1
    cuts.AddCutMergedCalo("00010113","1111121063032200000","1111121063022210001","0163303100000000"); // NLM 1 looser
    cuts.AddCutMergedCalo("00010113","1111121063032200000","1111121063022210001","0163305100000000"); // NLM 1 tighter
    cuts.AddCutMergedCalo("00010113","1112121063032200000","1112121063022210001","0163301100000000"); // TRD infront
    cuts.AddCutMergedCalo("00010113","1111321063032200000","1111321063022210001","0163301100000000"); // no TRD infront
  } else if (trainConfig == 54){ // varied exoctics
    cuts.AddCutMergedCalo("00010113","1111121063032200000","1111121063022700001","0163300000000000"); // no exotics
    cuts.AddCutMergedCalo("00010113","1111121063232200000","1111121063222700001","0163300000000000"); // frac = 0.99
    cuts.AddCutMergedCalo("00010113","1111121063332200000","1111121063322700001","0163300000000000"); // frac = 0.98
    cuts.AddCutMergedCalo("00010113","1111121063532200000","1111121063522700001","0163300000000000"); // frac = 0.97
    cuts.AddCutMergedCalo("00010113","1111121063732200000","1111121063722700001","0163300000000000"); // frac = 0.96
    cuts.AddCutMergedCalo("00010113","1111121063932200000","1111121063922700001","0163300000000000"); // frac = 0.94
  } else if (trainConfig == 55){  // variation of mass and alpha cut
    cuts.AddCutMergedCalo("00010113","1111121063032200000","1111121063022700001","0163301700000000"); // only band at 0
    cuts.AddCutMergedCalo("00010113","1111121063032200000","1111121063022700001","0163300700000000"); // only band at 0, no alpha
    cuts.AddCutMergedCalo("00010113","1111121063032200000","1111121063022700001","0163301000000000"); // no mass cut
    cuts.AddCutMergedCalo("00010113","1111121063032200000","1111121063022700001","0163300000000000"); // no mass cut, no alpha
    cuts.AddCutMergedCalo("00010113","1111121063032200000","1111121063022700001","0163300100000000"); // no alpha
  } else if (trainConfig == 56){ // varied M02
    cuts.AddCutMergedCalo("00010113","1111121063032200000","1111121063022600001","0163300700000000"); // min M02= 0.3
    cuts.AddCutMergedCalo("00010113","1111121063032200000","1111121063022800001","0163300700000000"); // min M02= 0.25
    cuts.AddCutMergedCalo("00010113","1111121063032200000","1111121063022600001","0163300000000000"); // no mass, min M02 = 0.3
    cuts.AddCutMergedCalo("00010113","1111121063032200000","1111121063022700001","0163300000000000"); // no mass, min M02 = 0.27
    cuts.AddCutMergedCalo("00010113","1111121063032200000","1111121063022800001","0163300000000000"); // no mass, min M02 = 0.25
  } else if (trainConfig == 57){  // variation track matching to cluster & mass variations new defaults
    cuts.AddCutMergedCalo("00010113","1111121060032200000","1111121060022700001","0163300000000000"); // no TM
    cuts.AddCutMergedCalo("00010113","1111121061032200000","1111121061022700001","0163300000000000"); // looser TM
    cuts.AddCutMergedCalo("00010113","1111121066032200000","1111121066022700001","0163300000000000"); // tighter TM
    cuts.AddCutMergedCalo("00010113","1112121063032200000","1112121063022700001","0163300000000000"); // TRD infront
    cuts.AddCutMergedCalo("00010113","1111321063032200000","1111321063022700001","0163300000000000"); // no TRD infront
  } else if (trainConfig == 58){ // NL var new defaults
    cuts.AddCutMergedCalo("00010113","1111122063032200000","1111122063022700001","0163300000000000"); // SDM loose time
    cuts.AddCutMergedCalo("00010113","1111111063032200000","1111111063022700001","0163300000000000"); // conv calo tight time
    cuts.AddCutMergedCalo("00010113","1111101063032200000","1111101063022700001","0163300000000000"); // SDM Jason
    cuts.AddCutMergedCalo("00010113","1111100063032200000","1111100063022700001","0163300000000000"); // none

  // EMC7 variations
  } else if (trainConfig == 60){ // NLM 1 M02 var
    cuts.AddCutMergedCalo("00052013","1111121063032200000","1111121063022110001","0163301100000000"); // min 0.3 function default
    cuts.AddCutMergedCalo("00052013","1111121063032200000","1111121063022310001","0163301100000000"); // min 0.25 function default
    cuts.AddCutMergedCalo("00052013","1111121063032200000","1111121063022410001","0163301100000000"); // min 0.27, tighter lower func
    cuts.AddCutMergedCalo("00052013","1111121063032200000","1111121063022520001","0163301100000000"); // min 0.27, looser func up and low
    cuts.AddCutMergedCalo("00052013","1111121063032200000","1111121063022430001","0163301100000000"); // min 0.27, tighter func up and low
  } else if (trainConfig == 61){  // EMCAL clusters, variation track matching to cluster & Mass EMC7 NLM 1
    cuts.AddCutMergedCalo("00052013","1111121060032200000","1111121060022210001","0163301100000000"); // no TM
    cuts.AddCutMergedCalo("00052013","1111121061032200000","1111121061022210001","0163301100000000"); // looser TM
    cuts.AddCutMergedCalo("00052013","1111121066032200000","1111121066022210001","0163301100000000"); // tighter TM
    cuts.AddCutMergedCalo("00052013","1111121063032200000","1111121063022210001","0163301300000000"); // tighter mass cut
    cuts.AddCutMergedCalo("00052013","1111121063032200000","1111121063022210001","0163301500000000"); // looser mass cut
  } else if (trainConfig == 62){  // NL variations EMC7
    cuts.AddCutMergedCalo("00052013","1111122063032200000","1111122063022210001","0163301100000000"); // SDM loose time
    cuts.AddCutMergedCalo("00052013","1111111063032200000","1111111063022210001","0163301100000000"); // conv calo tight time
    cuts.AddCutMergedCalo("00052013","1111101063032200000","1111101063022210001","0163301100000000"); // SDM Jason
    cuts.AddCutMergedCalo("00052013","1111100063032200000","1111100063022210001","0163301100000000"); // none
  } else if (trainConfig == 63){  // Alpha cut variations & TRD material EMC7 NLM 1
    cuts.AddCutMergedCalo("00052013","1111121063032200000","1111121063022210001","0163303100000000"); // NLM 1 looser
    cuts.AddCutMergedCalo("00052013","1111121063032200000","1111121063022210001","0163305100000000"); // NLM 1 tighter
    cuts.AddCutMergedCalo("00052013","1112121063032200000","1112121063022210001","0163301100000000"); // TRD infront
    cuts.AddCutMergedCalo("00052013","1111321063032200000","1111321063022210001","0163301100000000"); // no TRD infront
  } else if (trainConfig == 64){ // varied exoctics
    cuts.AddCutMergedCalo("00052013","1111121063032200000","1111121063022700001","0163300000000000"); // no exotics
    cuts.AddCutMergedCalo("00052013","1111121063232200000","1111121063222700001","0163300000000000"); // frac = 0.99
    cuts.AddCutMergedCalo("00052013","1111121063332200000","1111121063322700001","0163300000000000"); // frac = 0.98
    cuts.AddCutMergedCalo("00052013","1111121063532200000","1111121063522700001","0163300000000000"); // frac = 0.97
    cuts.AddCutMergedCalo("00052013","1111121063732200000","1111121063722700001","0163300000000000"); // frac = 0.96
    cuts.AddCutMergedCalo("00052013","1111121063932200000","1111121063922700001","0163300000000000"); // frac = 0.94
  } else if (trainConfig == 65){  // variation of mass and alpha cut
    cuts.AddCutMergedCalo("00052013","1111121063032200000","1111121063022700001","0163301700000000"); // only band at 0
    cuts.AddCutMergedCalo("00052013","1111121063032200000","1111121063022700001","0163300700000000"); // only band at 0, no alpha
    cuts.AddCutMergedCalo("00052013","1111121063032200000","1111121063022700001","0163301000000000"); // no mass cut
    cuts.AddCutMergedCalo("00052013","1111121063032200000","1111121063022700001","0163300000000000"); // no mass cut, no alpha
    cuts.AddCutMergedCalo("00052013","1111121063032200000","1111121063022700001","0163300100000000"); // no alpha
  } else if (trainConfig == 66){ // varied M02
    cuts.AddCutMergedCalo("00052013","1111121063032200000","1111121063022600001","0163300700000000"); // min M02= 0.3
    cuts.AddCutMergedCalo("00052013","1111121063032200000","1111121063022800001","0163300700000000"); // min M02= 0.25
    cuts.AddCutMergedCalo("00052013","1111121063032200000","1111121063022600001","0163300000000000"); // no mass, min M02 = 0.3
    cuts.AddCutMergedCalo("00052013","1111121063032200000","1111121063022700001","0163300000000000"); // no mass, min M02 = 0.27
    cuts.AddCutMergedCalo("00052013","1111121063032200000","1111121063022800001","0163300000000000"); // no mass, min M02 = 0.25
  } else if (trainConfig == 67){  // variation track matching to cluster & mass variations new defaults
    cuts.AddCutMergedCalo("00052013","1111121060032200000","1111121060022700001","0163300000000000"); // no TM
    cuts.AddCutMergedCalo("00052013","1111121061032200000","1111121061022700001","0163300000000000"); // looser TM
    cuts.AddCutMergedCalo("00052013","1111121066032200000","1111121066022700001","0163300000000000"); // tighter TM
    cuts.AddCutMergedCalo("00052013","1112121063032200000","1112121063022700001","0163300000000000"); // TRD infront
    cuts.AddCutMergedCalo("00052013","1111321063032200000","1111321063022700001","0163300000000000"); // no TRD infront
  } else if (trainConfig == 68){ // NL var new defaults
    cuts.AddCutMergedCalo("00052013","1111122063032200000","1111122063022700001","0163300000000000"); // SDM loose time
    cuts.AddCutMergedCalo("00052013","1111111063032200000","1111111063022700001","0163300000000000"); // conv calo tight time
    cuts.AddCutMergedCalo("00052013","1111101063032200000","1111101063022700001","0163300000000000"); // SDM Jason
    cuts.AddCutMergedCalo("00052013","1111100063032200000","1111100063022700001","0163300000000000"); // none

  // EG2 variations
  } else if (trainConfig == 70){ // NLM 1 M02 var
    cuts.AddCutMergedCalo("00085013","1111121063032200000","1111121063022110001","0163301100000000"); // min 0.3 function default
    cuts.AddCutMergedCalo("00085013","1111121063032200000","1111121063022310001","0163301100000000"); // min 0.25 function default
    cuts.AddCutMergedCalo("00085013","1111121063032200000","1111121063022410001","0163301100000000"); // min 0.27, tighter lower func
    cuts.AddCutMergedCalo("00085013","1111121063032200000","1111121063022520001","0163301100000000"); // min 0.27, looser func up and low
    cuts.AddCutMergedCalo("00085013","1111121063032200000","1111121063022430001","0163301100000000"); // min 0.27, tighter func up and low
  } else if (trainConfig == 71){  // EMCAL clusters, variation track matching to cluster & Mass EG2 NLM 1
    cuts.AddCutMergedCalo("00085013","1111121060032200000","1111121060022210001","0163301100000000"); // no TM
    cuts.AddCutMergedCalo("00085013","1111121061032200000","1111121061022210001","0163301100000000"); // looser TM
    cuts.AddCutMergedCalo("00085013","1111121066032200000","1111121066022210001","0163301100000000"); // tighter TM
    cuts.AddCutMergedCalo("00085013","1111121063032200000","1111121063022210001","0163301300000000"); // tighter mass cut
    cuts.AddCutMergedCalo("00085013","1111121063032200000","1111121063022210001","0163301500000000"); // looser mass cut
  } else if (trainConfig == 72){  // NL variations EG2
    cuts.AddCutMergedCalo("00085013","1111122063032200000","1111122063022210001","0163301100000000"); // SDM loose time
    cuts.AddCutMergedCalo("00085013","1111111063032200000","1111111063022210001","0163301100000000"); // conv calo tight time
    cuts.AddCutMergedCalo("00085013","1111101063032200000","1111101063022210001","0163301100000000"); // SDM Jason
    cuts.AddCutMergedCalo("00085013","1111100063032200000","1111100063022210001","0163301100000000"); // none
  } else if (trainConfig == 73){  // Alpha cut variations & TRD material EG2 NLM 1
    cuts.AddCutMergedCalo("00085013","1111121063032200000","1111121063022210001","0163303100000000"); // NLM 1 looser
    cuts.AddCutMergedCalo("00085013","1111121063032200000","1111121063022210001","0163305100000000"); // NLM 1 tighter
    cuts.AddCutMergedCalo("00085013","1112121063032200000","1112121063022210001","0163301100000000"); // TRD infront
    cuts.AddCutMergedCalo("00085013","1111321063032200000","1111321063022210001","0163301100000000"); // no TRD infront
  } else if (trainConfig == 74){ // varied exoctics
    cuts.AddCutMergedCalo("00085013","1111121063032200000","1111121063022700001","0163300000000000"); // no exotics
    cuts.AddCutMergedCalo("00085013","1111121063232200000","1111121063222700001","0163300000000000"); // frac = 0.99
    cuts.AddCutMergedCalo("00085013","1111121063332200000","1111121063322700001","0163300000000000"); // frac = 0.98
    cuts.AddCutMergedCalo("00085013","1111121063532200000","1111121063522700001","0163300000000000"); // frac = 0.97
    cuts.AddCutMergedCalo("00085013","1111121063732200000","1111121063722700001","0163300000000000"); // frac = 0.96
    cuts.AddCutMergedCalo("00085013","1111121063932200000","1111121063922700001","0163300000000000"); // frac = 0.94
  } else if (trainConfig == 75){  // variation of mass and alpha cut
    cuts.AddCutMergedCalo("00085013","1111121063032200000","1111121063022700001","0163301700000000"); // only band at 0
    cuts.AddCutMergedCalo("00085013","1111121063032200000","1111121063022700001","0163300700000000"); // only band at 0, no alpha
    cuts.AddCutMergedCalo("00085013","1111121063032200000","1111121063022700001","0163301000000000"); // no mass cut
    cuts.AddCutMergedCalo("00085013","1111121063032200000","1111121063022700001","0163300000000000"); // no mass cut, no alpha
    cuts.AddCutMergedCalo("00085013","1111121063032200000","1111121063022700001","0163300100000000"); // no alpha
  } else if (trainConfig == 76){ // varied M02
    cuts.AddCutMergedCalo("00085013","1111121063032200000","1111121063022600001","0163300700000000"); // min M02= 0.3
    cuts.AddCutMergedCalo("00085013","1111121063032200000","1111121063022800001","0163300700000000"); // min M02= 0.25
    cuts.AddCutMergedCalo("00085013","1111121063032200000","1111121063022600001","0163300000000000"); // no mass, min M02 = 0.3
    cuts.AddCutMergedCalo("00085013","1111121063032200000","1111121063022700001","0163300000000000"); // no mass, min M02 = 0.27
    cuts.AddCutMergedCalo("00085013","1111121063032200000","1111121063022800001","0163300000000000"); // no mass, min M02 = 0.25
  } else if (trainConfig == 77){  // variation track matching to cluster & mass variations new defaults
    cuts.AddCutMergedCalo("00085013","1111121060032200000","1111121060022700001","0163300000000000"); // no TM
    cuts.AddCutMergedCalo("00085013","1111121061032200000","1111121061022700001","0163300000000000"); // looser TM
    cuts.AddCutMergedCalo("00085013","1111121066032200000","1111121066022700001","0163300000000000"); // tighter TM
    cuts.AddCutMergedCalo("00085013","1112121063032200000","1112121063022700001","0163300000000000"); // TRD infront
    cuts.AddCutMergedCalo("00085013","1111321063032200000","1111321063022700001","0163300000000000"); // no TRD infront
  } else if (trainConfig == 78){ // NL var new defaults
    cuts.AddCutMergedCalo("00085013","1111122063032200000","1111122063022700001","0163300000000000"); // SDM loose time
    cuts.AddCutMergedCalo("00085013","1111111063032200000","1111111063022700001","0163300000000000"); // conv calo tight time
    cuts.AddCutMergedCalo("00085013","1111101063032200000","1111101063022700001","0163300000000000"); // SDM Jason
    cuts.AddCutMergedCalo("00085013","1111100063032200000","1111100063022700001","0163300000000000"); // none

  // EG1 variations
  } else if (trainConfig == 80){ // NLM 1 M02 var
    cuts.AddCutMergedCalo("00083013","1111121063032200000","1111121063022110001","0163301100000000"); // min 0.3 function default
    cuts.AddCutMergedCalo("00083013","1111121063032200000","1111121063022310001","0163301100000000"); // min 0.25 function default
    cuts.AddCutMergedCalo("00083013","1111121063032200000","1111121063022410001","0163301100000000"); // min 0.27, tighter lower func
    cuts.AddCutMergedCalo("00083013","1111121063032200000","1111121063022520001","0163301100000000"); // min 0.27, looser func up and low
    cuts.AddCutMergedCalo("00083013","1111121063032200000","1111121063022430001","0163301100000000"); // min 0.27, tighter func up and low
  } else if (trainConfig == 81){  // EMCAL clusters, variation track matching to cluster & Mass EG1 NLM 1
    cuts.AddCutMergedCalo("00083013","1111121060032200000","1111121060022210001","0163301100000000"); // no TM
    cuts.AddCutMergedCalo("00083013","1111121061032200000","1111121061022210001","0163301100000000"); // looser TM
    cuts.AddCutMergedCalo("00083013","1111121066032200000","1111121066022210001","0163301100000000"); // tighter TM
    cuts.AddCutMergedCalo("00083013","1111121063032200000","1111121063022210001","0163301300000000"); // tighter mass cut
    cuts.AddCutMergedCalo("00083013","1111121063032200000","1111121063022210001","0163301500000000"); // looser mass cut
  } else if (trainConfig == 82){  // NL variations EG1
    cuts.AddCutMergedCalo("00083013","1111122063032200000","1111122063022210001","0163301100000000"); // SDM loose time
    cuts.AddCutMergedCalo("00083013","1111111063032200000","1111111063022210001","0163301100000000"); // conv calo tight time
    cuts.AddCutMergedCalo("00083013","1111101063032200000","1111101063022210001","0163301100000000"); // SDM Jason
    cuts.AddCutMergedCalo("00083013","1111100063032200000","1111100063022210001","0163301100000000"); // none
  } else if (trainConfig == 83){  // Alpha cut variations & TRD material EG1 NLM 1
    cuts.AddCutMergedCalo("00083013","1111121063032200000","1111121063022210001","0163303100000000"); // NLM 1 looser
    cuts.AddCutMergedCalo("00083013","1111121063032200000","1111121063022210001","0163305100000000"); // NLM 1 tighter
    cuts.AddCutMergedCalo("00083013","1112121063032200000","1112121063022210001","0163301100000000"); // TRD infront
    cuts.AddCutMergedCalo("00083013","1111321063032200000","1111321063022210001","0163301100000000"); // no TRD infront
  } else if (trainConfig == 84){ // varied exoctics
    cuts.AddCutMergedCalo("00083013","1111121063032200000","1111121063022700001","0163300000000000"); // no exotics
    cuts.AddCutMergedCalo("00083013","1111121063232200000","1111121063222700001","0163300000000000"); // frac = 0.99
    cuts.AddCutMergedCalo("00083013","1111121063332200000","1111121063322700001","0163300000000000"); // frac = 0.98
    cuts.AddCutMergedCalo("00083013","1111121063532200000","1111121063522700001","0163300000000000"); // frac = 0.97
    cuts.AddCutMergedCalo("00083013","1111121063732200000","1111121063722700001","0163300000000000"); // frac = 0.96
    cuts.AddCutMergedCalo("00083013","1111121063932200000","1111121063922700001","0163300000000000"); // frac = 0.94
  } else if (trainConfig == 85){  // variation of mass and alpha cut
    cuts.AddCutMergedCalo("00083013","1111121063032200000","1111121063022700001","0163301700000000"); // only band at 0
    cuts.AddCutMergedCalo("00083013","1111121063032200000","1111121063022700001","0163300700000000"); // only band at 0, no alpha
    cuts.AddCutMergedCalo("00083013","1111121063032200000","1111121063022700001","0163301000000000"); // no mass cut
    cuts.AddCutMergedCalo("00083013","1111121063032200000","1111121063022700001","0163300000000000"); // no mass cut, no alpha
    cuts.AddCutMergedCalo("00083013","1111121063032200000","1111121063022700001","0163300100000000"); // no alpha
  } else if (trainConfig == 86){ // varied M02
    cuts.AddCutMergedCalo("00083013","1111121063032200000","1111121063022600001","0163300700000000"); // min M02= 0.3
    cuts.AddCutMergedCalo("00083013","1111121063032200000","1111121063022800001","0163300700000000"); // min M02= 0.25
    cuts.AddCutMergedCalo("00083013","1111121063032200000","1111121063022600001","0163300000000000"); // no mass, min M02 = 0.3
    cuts.AddCutMergedCalo("00083013","1111121063032200000","1111121063022700001","0163300000000000"); // no mass, min M02 = 0.27
    cuts.AddCutMergedCalo("00083013","1111121063032200000","1111121063022800001","0163300000000000"); // no mass, min M02 = 0.25
  } else if (trainConfig == 87){  // variation track matching to cluster & mass variations new defaults
    cuts.AddCutMergedCalo("00083013","1111121060032200000","1111121060022700001","0163300000000000"); // no TM
    cuts.AddCutMergedCalo("00083013","1111121061032200000","1111121061022700001","0163300000000000"); // looser TM
    cuts.AddCutMergedCalo("00083013","1111121066032200000","1111121066022700001","0163300000000000"); // tighter TM
    cuts.AddCutMergedCalo("00083013","1112121063032200000","1112121063022700001","0163300000000000"); // TRD infront
    cuts.AddCutMergedCalo("00083013","1111321063032200000","1111321063022700001","0163300000000000"); // no TRD infront
  } else if (trainConfig == 88){ // NL var new defaults
    cuts.AddCutMergedCalo("00083013","1111122063032200000","1111122063022700001","0163300000000000"); // SDM loose time
    cuts.AddCutMergedCalo("00083013","1111111063032200000","1111111063022700001","0163300000000000"); // conv calo tight time
    cuts.AddCutMergedCalo("00083013","1111101063032200000","1111101063022700001","0163300000000000"); // SDM Jason
    cuts.AddCutMergedCalo("00083013","1111100063032200000","1111100063022700001","0163300000000000"); // none

  } else if (trainConfig == 90){  // new default without mass, eta < 0.7, y < 0.7
    cuts.AddCutMergedCalo("00010113","1551121063032200000","1551121063022700001","0163200000000000");
    cuts.AddCutMergedCalo("00052013","1551121063032200000","1551121063022700001","0163200000000000");
    cuts.AddCutMergedCalo("00085013","1551121063032200000","1551121063022700001","0163200000000000");
    cuts.AddCutMergedCalo("00083013","1551121063032200000","1551121063022700001","0163200000000000");
  } else if (trainConfig == 91){  // new default without mass, TM pt dep, eta < 0.7, y < 0.7
    cuts.AddCutMergedCalo("00010113","1551121067032200000","1551121067022700001","0163200000000000");
    cuts.AddCutMergedCalo("00052013","1551121067032200000","1551121067022700001","0163200000000000");
    cuts.AddCutMergedCalo("00085013","1551121067032200000","1551121067022700001","0163200000000000");
    cuts.AddCutMergedCalo("00083013","1551121067032200000","1551121067022700001","0163200000000000");


  } else if (trainConfig == 95){  // new defaults LHC13g no NLM no mass, no alpha, no M02
    cuts.AddCutMergedCalo("00010113","1111121063032200000","1111121063022000000","0163300000000000"); // INT7
    cuts.AddCutMergedCalo("00052013","1111121063032200000","1111121063022000000","0163300000000000"); // EMC7
    cuts.AddCutMergedCalo("00085013","1111121063032200000","1111121063022000000","0163300000000000"); // EG2
    cuts.AddCutMergedCalo("00083013","1111121063032200000","1111121063022000000","0163300000000000"); // EG1
  } else if (trainConfig == 96){  // new defaults LHC13g no NLM, no mass
    cuts.AddCutMergedCalo("00010113","1111121063032200000","1111121063022700000","0163300000000000"); // INT7
    cuts.AddCutMergedCalo("00052013","1111121063032200000","1111121063022700000","0163300000000000"); // EMC7
    cuts.AddCutMergedCalo("00085013","1111121063032200000","1111121063022700000","0163300000000000"); // EG2
    cuts.AddCutMergedCalo("00083013","1111121063032200000","1111121063022700000","0163300000000000"); // EG1
  } else if (trainConfig == 97){  // new defaults LHC13g NLM2
    cuts.AddCutMergedCalo("00010113","1111121063032200000","1111121063022700002","0163300700000000"); // INT7
    cuts.AddCutMergedCalo("00052013","1111121063032200000","1111121063022700002","0163300700000000"); // EMC7
    cuts.AddCutMergedCalo("00085013","1111121063032200000","1111121063022700002","0163300700000000"); // EG2
    cuts.AddCutMergedCalo("00083013","1111121063032200000","1111121063022700002","0163300700000000"); // EG1
  } else if (trainConfig == 98){  // new defaults LHC13g NLM2
    cuts.AddCutMergedCalo("00010113","1111121063032200000","1111121063022210002","0163302200000000"); // INT7
    cuts.AddCutMergedCalo("00052013","1111121063032200000","1111121063022210002","0163302200000000"); // EMC7
    cuts.AddCutMergedCalo("00085013","1111121063032200000","1111121063022210002","0163302200000000"); // EG2
    cuts.AddCutMergedCalo("00083013","1111121063032200000","1111121063022210002","0163302200000000"); // EG1
  } else if (trainConfig == 99){  // NLM2 no M02, no mass, no alpah
    cuts.AddCutMergedCalo("00010113","1111121063032200000","1111121063022000002","0163300000000000"); // INT7
    cuts.AddCutMergedCalo("00052013","1111121063032200000","1111121063022000002","0163300000000000"); // EMC7
    cuts.AddCutMergedCalo("00085013","1111121063032200000","1111121063022000002","0163300000000000"); // EG2
    cuts.AddCutMergedCalo("00083013","1111121063032200000","1111121063022000002","0163300000000000"); // EG1

  // LHC12
  // default
  } else if (trainConfig == 101){ // pPb cuts : open timing, TB nonlin, with TM
    cuts.AddCutMergedCalo("00010113","1111101017032200000","1111101017022700001","0163300000000000"); // INT7
    cuts.AddCutMergedCalo("00052113","1111101017032200000","1111101017022700001","0163300000000000"); // EMC7
    cuts.AddCutMergedCalo("00081113","1111101017032200000","1111101017022700001","0163300000000000"); // EGA
  } else if (trainConfig == 102){ //  pPb cuts : open timing, TB nonlin, without TM
    cuts.AddCutMergedCalo("00010113","1111101010032200000","1111101010022700001","0163300000000000"); // INT7
    cuts.AddCutMergedCalo("00052113","1111101010032200000","1111101010022700001","0163300000000000"); // EMC7
    cuts.AddCutMergedCalo("00081113","1111101010032200000","1111101010022700001","0163300000000000"); // EGA
  } else if (trainConfig == 103){ //  pPb cuts : open timing, TB nonlin, without TM
    cuts.AddCutMergedCalo("00010113","1111101010032200000","1111101010022700001","0163300000000000"); // INT7
  } else if (trainConfig == 104){ //  pPb cuts : open timing, TB nonlin, without TM
    cuts.AddCutMergedCalo("00052113","1111101010032200000","1111101010022700001","0163300000000000"); // EMC7
  } else if (trainConfig == 105){ //  pPb cuts : open timing, TB nonlin, without TM
    cuts.AddCutMergedCalo("00081113","1111101010032200000","1111101010022700001","0163300000000000"); // EGA
  } else if (trainConfig == 107){  // no M02, pt dep TM
    cuts.AddCutMergedCalo("00010113","1111111067032200000","1111111067022000001","0163300000000000"); // INT7
    cuts.AddCutMergedCalo("00052113","1111111067032200000","1111111067022000001","0163300000000000"); // EMC7
    cuts.AddCutMergedCalo("00081113","1111111067032200000","1111111067022000001","0163300000000000"); // EGA
  } else if (trainConfig == 108){  // no M02, TM only in merged, pt dep TM
    cuts.AddCutMergedCalo("00010113","1111111060032200000","1111111067022000001","0163300000000000"); // INT7
    cuts.AddCutMergedCalo("00052113","1111111060032200000","1111111067022000001","0163300000000000"); // EMC7
    cuts.AddCutMergedCalo("00081113","1111111060032200000","1111111067022000001","0163300000000000"); // EGA
  } else if (trainConfig == 109){  // M02 cut at 0.27, pt dep TM
    cuts.AddCutMergedCalo("00010113","1111111067032200000","1111111067022700001","0163300000000000"); // INT7
    cuts.AddCutMergedCalo("00052113","1111111067032200000","1111111067022700001","0163300000000000"); // EMC7
    cuts.AddCutMergedCalo("00081113","1111111067032200000","1111111067022700001","0163300000000000"); // EGA
  } else if (trainConfig == 110){  // M02 cut at 0.27, TM only in merged, pt dep TM
    cuts.AddCutMergedCalo("00010113","1111111060032200000","1111111067022700001","0163300000000000"); // INT7
    cuts.AddCutMergedCalo("00052113","1111111060032200000","1111111067022700001","0163300000000000"); // EMC7
    cuts.AddCutMergedCalo("00081113","1111111060032200000","1111111067022700001","0163300000000000"); // EGA

  } else if (trainConfig == 111){  // EMCAL clusters, different triggers no NonLinearity
    cuts.AddCutMergedCalo("00010113","1111100060032200000","1111100060022210001","0163301100000000"); // INT7
    cuts.AddCutMergedCalo("00052113","1111100060032200000","1111100060022210001","0163301100000000"); // EMC7
    cuts.AddCutMergedCalo("00081113","1111100060032200000","1111100060022110001","0163301100000000"); // EMCEGA,

  } else if (trainConfig == 112){  // EMCAL clusters, different triggers with NonLinearity
    cuts.AddCutMergedCalo("00010113","1111111060032200000","1111111060022210001","0163301100000000"); // INT7
    cuts.AddCutMergedCalo("00052113","1111111060032200000","1111111060022210001","0163301100000000"); // EMC7
    cuts.AddCutMergedCalo("00081113","1111111060032200000","1111111060022210001","0163301100000000"); // EMCEGA,
  } else if (trainConfig == 113){  // EMCAL clusters, different triggers with NonLinearity track matching to cluster
    cuts.AddCutMergedCalo("00010113","1111111067032200000","1111111067022210001","0163301100000000"); // INT7
    cuts.AddCutMergedCalo("00052113","1111111067032200000","1111111067022210001","0163301100000000"); // EMC7
    cuts.AddCutMergedCalo("00081113","1111111067032200000","1111111067022210001","0163301100000000"); // EMCEGA,

  // new default
  } else if (trainConfig == 114){  // NLM1 no M02, no mass, no alpah
    cuts.AddCutMergedCalo("00010113","1111111067032200000","1111111067022000001","0163300000000000"); // INT7
    cuts.AddCutMergedCalo("00052113","1111111067032200000","1111111067022000001","0163300000000000"); // EMC7
    cuts.AddCutMergedCalo("00081113","1111111067032200000","1111111067022000001","0163300000000000"); // EGA
    cuts.AddCutMergedCalo("00091113","1111111067032200000","1111111067022000001","0163300000000000"); // EJE
  } else if (trainConfig == 115){  // Mass only band at 0, M02 cut at 0.27
    cuts.AddCutMergedCalo("00010113","1111111067032200000","1111111067022700001","0163300700000000"); // INT7
    cuts.AddCutMergedCalo("00052113","1111111067032200000","1111111067022700001","0163300700000000"); // EMC7
    cuts.AddCutMergedCalo("00081113","1111111067032200000","1111111067022700001","0163300700000000"); // EGA
  } else if (trainConfig == 116){  // M02 cut at 0.27
    cuts.AddCutMergedCalo("00010113","1111111067032200000","1111111067022700001","0163300000000000"); // INT7
    cuts.AddCutMergedCalo("00052113","1111111067032200000","1111111067022700001","0163300000000000"); // EMC7
    cuts.AddCutMergedCalo("00081113","1111111067032200000","1111111067022700001","0163300000000000"); // EGA
    cuts.AddCutMergedCalo("00091113","1111111067032200000","1111111067022700001","0163300000000000"); // EJE
  } else if (trainConfig == 117){  // M02 cut at 0.27, TM only in merged
    cuts.AddCutMergedCalo("00010113","1111111060032200000","1111111067022700001","0163300000000000"); // INT7
    cuts.AddCutMergedCalo("00052113","1111111060032200000","1111111067022700001","0163300000000000"); // EMC7
    cuts.AddCutMergedCalo("00081113","1111111060032200000","1111111067022700001","0163300000000000"); // EGA
  } else if (trainConfig == 118){  // NLM1 no M02, TM only in merged
    cuts.AddCutMergedCalo("00010113","1111111060032200000","1111111067022000001","0163300000000000"); // INT7
    cuts.AddCutMergedCalo("00052113","1111111060032200000","1111111067022000001","0163300000000000"); // EMC7
    cuts.AddCutMergedCalo("00081113","1111111060032200000","1111111067022000001","0163300000000000"); // EGA
  } else if (trainConfig == 119){  // new default, with eta < 0.7, y < 0.7
    cuts.AddCutMergedCalo("00010113","1551111067032200000","1551111067022700001","0163200000000000"); // INT7
    cuts.AddCutMergedCalo("00052113","1551111067032200000","1551111067022700001","0163200000000000"); // EMC7
    cuts.AddCutMergedCalo("00081113","1551111067032200000","1551111067022700001","0163200000000000"); // EGA

    //kINT7 - NLM1
  } else if (trainConfig == 120){  // TRD material INT7 NLM 1
    cuts.AddCutMergedCalo("00010113","1112111067032200000","1112111067022700001","0163300000000000"); // TRD infront
    cuts.AddCutMergedCalo("00010113","1111311067032200000","1111311067022700001","0163300000000000"); // no TRD infront
  } else if (trainConfig == 121){ // varied exoctics
    cuts.AddCutMergedCalo("00010113","1111111067032200000","1111111067022700001","0163300000000000"); // no exotics
    cuts.AddCutMergedCalo("00010113","1111111067232200000","1111111067222700001","0163300000000000"); // frac = 0.99
    cuts.AddCutMergedCalo("00010113","1111111067332200000","1111111067322700001","0163300000000000"); // frac = 0.98
    cuts.AddCutMergedCalo("00010113","1111111067532200000","1111111067522700001","0163300000000000"); // frac = 0.97
    cuts.AddCutMergedCalo("00010113","1111111067732200000","1111111067722700001","0163300000000000"); // frac = 0.96
    cuts.AddCutMergedCalo("00010113","1111111067932200000","1111111067922700001","0163300000000000"); // frac = 0.94
  } else if (trainConfig == 122){ // varied M02 part 1
    cuts.AddCutMergedCalo("00010113","1111111067032200000","1111111067022600001","0163300000000000"); // min M02 = 0.3
    cuts.AddCutMergedCalo("00010113","1111111067032200000","1111111067022c00001","0163300000000000"); // min M02 = 0.29
    cuts.AddCutMergedCalo("00010113","1111111067032200000","1111111067022b00001","0163300000000000"); // min M02 = 0.28
  } else if (trainConfig == 123){ // varied M02 part 2
    cuts.AddCutMergedCalo("00010113","1111111067032200000","1111111067022700001","0163300000000000"); // min M02 = 0.27
    cuts.AddCutMergedCalo("00010113","1111111067032200000","1111111067022a00001","0163300000000000"); // min M02 = 0.26
    cuts.AddCutMergedCalo("00010113","1111111067032200000","1111111067022800001","0163300000000000"); // min M02 = 0.25
  } else if (trainConfig == 124){  // variation track matching to cluster part 1
    cuts.AddCutMergedCalo("00010113","1111111060032200000","1111111060022700001","0163300000000000"); // no TM
    cuts.AddCutMergedCalo("00010113","1111111063032200000","1111111063022700001","0163300000000000"); // old matching std
  } else if (trainConfig == 125){  // variation track matching to cluster part 2
    cuts.AddCutMergedCalo("00010113","111111106b032200000","111111106b022700001","0163300000000000"); // looser TM
    cuts.AddCutMergedCalo("00010113","111111106a032200000","111111106a022700001","0163300000000000"); // tighter TM
    cuts.AddCutMergedCalo("00010113","1111111068032200000","1111111068022700001","0163300000000000"); // even looser TM
    cuts.AddCutMergedCalo("00010113","1111111066032200000","1111111066022700001","0163300000000000"); // even tighter TM
  } else if (trainConfig == 126){ // NL var new defaults
    cuts.AddCutMergedCalo("00010113","1111112067032200000","1111112067022700001","0163300000000000"); // SDM loose time
    cuts.AddCutMergedCalo("00010113","1111101067032200000","1111101067022700001","0163300000000000"); // SDM Jason
    cuts.AddCutMergedCalo("00010113","1111100067032200000","1111100067022700001","0163300000000000"); // none
    cuts.AddCutMergedCalo("00010113","1111121067032200000","1111121067022700001","0163300000000000"); // DPOW NL ConvCalo
    cuts.AddCutMergedCalo("00010113","1111122067032200000","1111122067022700001","0163300000000000"); // DPOW NL Calo
  } else if (trainConfig == 127){  // std pT dep, |eta| < 0.3, y < 0.3
    cuts.AddCutMergedCalo("00010113","1661111067032200000","1661111067022700001","0163700000000000"); // diff eta/rap cuts
  } else if (trainConfig == 128){ // NL var 1
    cuts.AddCutMergedCalo("00010113","1111100067032200000","1111100067022700001","0163300000000000"); // NL off
    cuts.AddCutMergedCalo("00010113","1111102067032200000","1111102067022700001","0163300000000000"); // Testbeam (v3) for data, Pi0MCv3 for MC
    cuts.AddCutMergedCalo("00010113","1111106067032200000","1111106067022700001","0163300000000000"); // Testbeam (v4) for data, Pi0MCv3 for MC
    cuts.AddCutMergedCalo("00010113","1111111067032200000","1111111067022700001","0163300000000000"); // default (CCRF)
  } else if (trainConfig == 129){ // NL var 2
    cuts.AddCutMergedCalo("00010113","1111112067032200000","1111112067022700001","0163300000000000"); // CRF
    cuts.AddCutMergedCalo("00010113","1111121067032200000","1111121067022700001","0163300000000000"); // CCMF
    cuts.AddCutMergedCalo("00010113","1111122067032200000","1111122067022700001","0163300000000000"); // CMF

    //kEMC7 - NLM1
  } else if (trainConfig == 130){  // TRD material INT7 NLM 1
    cuts.AddCutMergedCalo("00052113","1112111067032200000","1112111067022700001","0163300000000000"); // TRD infront
    cuts.AddCutMergedCalo("00052113","1111311067032200000","1111311067022700001","0163300000000000"); // no TRD infront
  } else if (trainConfig == 131){ // varied exoctics
    cuts.AddCutMergedCalo("00052113","1111111067032200000","1111111067022700001","0163300000000000"); // no exotics
    cuts.AddCutMergedCalo("00052113","1111111067232200000","1111111067222700001","0163300000000000"); // frac = 0.99
    cuts.AddCutMergedCalo("00052113","1111111067332200000","1111111067322700001","0163300000000000"); // frac = 0.98
    cuts.AddCutMergedCalo("00052113","1111111067532200000","1111111067522700001","0163300000000000"); // frac = 0.97
    cuts.AddCutMergedCalo("00052113","1111111067732200000","1111111067722700001","0163300000000000"); // frac = 0.96
    cuts.AddCutMergedCalo("00052113","1111111067932200000","1111111067922700001","0163300000000000"); // frac = 0.94
  } else if (trainConfig == 132){ // varied M02 part 1
    cuts.AddCutMergedCalo("00052113","1111111067032200000","1111111067022600001","0163300000000000"); // min M02 = 0.3
    cuts.AddCutMergedCalo("00052113","1111111067032200000","1111111067022c00001","0163300000000000"); // min M02 = 0.29
    cuts.AddCutMergedCalo("00052113","1111111067032200000","1111111067022b00001","0163300000000000"); // min M02 = 0.28
  } else if (trainConfig == 133){ // varied M02 part 2
    cuts.AddCutMergedCalo("00052113","1111111067032200000","1111111067022700001","0163300000000000"); // min M02 = 0.27
    cuts.AddCutMergedCalo("00052113","1111111067032200000","1111111067022a00001","0163300000000000"); // min M02 = 0.26
    cuts.AddCutMergedCalo("00052113","1111111067032200000","1111111067022800001","0163300000000000"); // min M02 = 0.25
  } else if (trainConfig == 134){  // variation track matching to cluster part 1
    cuts.AddCutMergedCalo("00052113","1111111060032200000","1111111060022700001","0163300000000000"); // no TM
    cuts.AddCutMergedCalo("00052113","1111111063032200000","1111111063022700001","0163300000000000"); // old matching std
  } else if (trainConfig == 135){  // variation track matching to cluster part 2
    cuts.AddCutMergedCalo("00052113","111111106b032200000","111111106b022700001","0163300000000000"); // looser TM
    cuts.AddCutMergedCalo("00052113","111111106a032200000","111111106a022700001","0163300000000000"); // tighter TM
    cuts.AddCutMergedCalo("00052113","1111111068032200000","1111111068022700001","0163300000000000"); // even looser TM
    cuts.AddCutMergedCalo("00052113","1111111066032200000","1111111066022700001","0163300000000000"); // even tighter TM
  } else if (trainConfig == 136){ // NL var new defaults
    cuts.AddCutMergedCalo("00052113","1111112067032200000","1111112067022700001","0163300000000000"); // SDM loose time
    cuts.AddCutMergedCalo("00052113","1111101067032200000","1111101067022700001","0163300000000000"); // SDM Jason
    cuts.AddCutMergedCalo("00052113","1111100067032200000","1111100067022700001","0163300000000000"); // none
    cuts.AddCutMergedCalo("00052113","1111121067032200000","1111121067022700001","0163300000000000"); // DPOW NL ConvCalo
    cuts.AddCutMergedCalo("00052113","1111122067032200000","1111122067022700001","0163300000000000"); // DPOW NL Calo
  } else if (trainConfig == 137){  // std pT dep, |eta| < 0.3, y < 0.3
    cuts.AddCutMergedCalo("00052113","1661111067032200000","1661111067022700001","0163700000000000"); // diff eta/rap cuts
  } else if (trainConfig == 138){ // NL var 1
    cuts.AddCutMergedCalo("00052113","1111100067032200000","1111100067022700001","0163300000000000"); // NL off
    cuts.AddCutMergedCalo("00052113","1111102067032200000","1111102067022700001","0163300000000000"); // Testbeam (v3) for data, Pi0MCv3 for MC
    cuts.AddCutMergedCalo("00052113","1111106067032200000","1111106067022700001","0163300000000000"); // Testbeam (v4) for data, Pi0MCv3 for MC
    cuts.AddCutMergedCalo("00052113","1111111067032200000","1111111067022700001","0163300000000000"); // default (CCRF)
  } else if (trainConfig == 139){ // NL var 2
    cuts.AddCutMergedCalo("00052113","1111112067032200000","1111112067022700001","0163300000000000"); // CRF
    cuts.AddCutMergedCalo("00052113","1111121067032200000","1111121067022700001","0163300000000000"); // CCMF
    cuts.AddCutMergedCalo("00052113","1111122067032200000","1111122067022700001","0163300000000000"); // CMF

  //kEMCEGA - NLM1
  } else if (trainConfig == 140){  // TRD material INT7 NLM 1
    cuts.AddCutMergedCalo("00081113","1112111067032200000","1112111067022700001","0163300000000000"); // TRD infront
    cuts.AddCutMergedCalo("00081113","1111311067032200000","1111311067022700001","0163300000000000"); // no TRD infront
  } else if (trainConfig == 141){ // varied exoctics
    cuts.AddCutMergedCalo("00081113","1111111067032200000","1111111067022700001","0163300000000000"); // no exotics
    cuts.AddCutMergedCalo("00081113","1111111067232200000","1111111067222700001","0163300000000000"); // frac = 0.99
    cuts.AddCutMergedCalo("00081113","1111111067332200000","1111111067322700001","0163300000000000"); // frac = 0.98
    cuts.AddCutMergedCalo("00081113","1111111067532200000","1111111067522700001","0163300000000000"); // frac = 0.97
    cuts.AddCutMergedCalo("00081113","1111111067732200000","1111111067722700001","0163300000000000"); // frac = 0.96
    cuts.AddCutMergedCalo("00081113","1111111067932200000","1111111067922700001","0163300000000000"); // frac = 0.94
  } else if (trainConfig == 142){ // varied M02 part 1
    cuts.AddCutMergedCalo("00081113","1111111067032200000","1111111067022600001","0163300000000000"); // min M02 = 0.3
    cuts.AddCutMergedCalo("00081113","1111111067032200000","1111111067022c00001","0163300000000000"); // min M02 = 0.29
    cuts.AddCutMergedCalo("00081113","1111111067032200000","1111111067022b00001","0163300000000000"); // min M02 = 0.28
  } else if (trainConfig == 143){ // varied M02 part 2
    cuts.AddCutMergedCalo("00081113","1111111067032200000","1111111067022700001","0163300000000000"); // min M02 = 0.27
    cuts.AddCutMergedCalo("00081113","1111111067032200000","1111111067022a00001","0163300000000000"); // min M02 = 0.26
    cuts.AddCutMergedCalo("00081113","1111111067032200000","1111111067022800001","0163300000000000"); // min M02 = 0.25
  } else if (trainConfig == 144){  // variation track matching to cluster part 1
    cuts.AddCutMergedCalo("00081113","1111111060032200000","1111111060022700001","0163300000000000"); // no TM
    cuts.AddCutMergedCalo("00081113","1111111063032200000","1111111063022700001","0163300000000000"); // old matching std
  } else if (trainConfig == 145){  // variation track matching to cluster part 2
    cuts.AddCutMergedCalo("00081113","111111106b032200000","111111106b022700001","0163300000000000"); // looser TM
    cuts.AddCutMergedCalo("00081113","111111106a032200000","111111106a022700001","0163300000000000"); // tighter TM
    cuts.AddCutMergedCalo("00081113","1111111068032200000","1111111068022700001","0163300000000000"); // even looser TM
    cuts.AddCutMergedCalo("00081113","1111111066032200000","1111111066022700001","0163300000000000"); // even tighter TM
  } else if (trainConfig == 146){ // NL var new defaults
    cuts.AddCutMergedCalo("00081113","1111112067032200000","1111112067022700001","0163300000000000"); // SDM loose time
    cuts.AddCutMergedCalo("00081113","1111101067032200000","1111101067022700001","0163300000000000"); // SDM Jason
    cuts.AddCutMergedCalo("00081113","1111100067032200000","1111100067022700001","0163300000000000"); // none
    cuts.AddCutMergedCalo("00081113","1111121067032200000","1111121067022700001","0163300000000000"); // DPOW NL ConvCalo
    cuts.AddCutMergedCalo("00081113","1111122067032200000","1111122067022700001","0163300000000000"); // DPOW NL Calo
  } else if (trainConfig == 147){  // std pT dep, |eta| < 0.3, y < 0.3
    cuts.AddCutMergedCalo("00081113","1661111067032200000","1661111067022700001","0163700000000000"); // diff eta/rap cuts
  } else if (trainConfig == 148){ // NL var 1
    cuts.AddCutMergedCalo("00081113","1111100067032200000","1111100067022700001","0163300000000000"); // NL off
    cuts.AddCutMergedCalo("00081113","1111102067032200000","1111102067022700001","0163300000000000"); // Testbeam (v3) for data, Pi0MCv3 for MC
    cuts.AddCutMergedCalo("00081113","1111106067032200000","1111106067022700001","0163300000000000"); // Testbeam (v4) for data, Pi0MCv3 for MC
    cuts.AddCutMergedCalo("00081113","1111111067032200000","1111111067022700001","0163300000000000"); // default (CCRF)
  } else if (trainConfig == 149){ // NL var 2
    cuts.AddCutMergedCalo("00081113","1111112067032200000","1111112067022700001","0163300000000000"); // CRF
    cuts.AddCutMergedCalo("00081113","1111121067032200000","1111121067022700001","0163300000000000"); // CCMF
    cuts.AddCutMergedCalo("00081113","1111122067032200000","1111122067022700001","0163300000000000"); // CMF
  // NLM2 cuts
  } else if (trainConfig == 160){  //no NLM no explicit exotics cut, no M02 cut, no mass
    cuts.AddCutMergedCalo("00010113","1111111067032200000","1111111067022000000","0163300000000000"); // INT7
    cuts.AddCutMergedCalo("00052113","1111111067032200000","1111111067022000000","0163300000000000"); // EMC7
    cuts.AddCutMergedCalo("00081113","1111111067032200000","1111111067022000000","0163300000000000"); // EGA
  } else if (trainConfig == 161){  // no NLM no mass, no explicit exotics cut, M02 cut at 0.27
    cuts.AddCutMergedCalo("00010113","1111111067032200000","1111111067022700000","0163300000000000"); // INT7
    cuts.AddCutMergedCalo("00052113","1111111067032200000","1111111067022700000","0163300000000000"); // EMC7
    cuts.AddCutMergedCalo("00081113","1111111067032200000","1111111067022700000","0163300000000000"); // EGA

  } else if (trainConfig == 162){  // NLM2 no M02, no mass, no alpha
    cuts.AddCutMergedCalo("00010113","1111111067032200000","1111111067022000002","0163300000000000"); // INT7
    cuts.AddCutMergedCalo("00052113","1111111067032200000","1111111067022000002","0163300000000000"); // EMC7
    cuts.AddCutMergedCalo("00081113","1111111067032200000","1111111067022000002","0163300000000000"); // EGA
  } else if (trainConfig == 163){  // NLM2 Mass only band at 0, no explicit exotics cut, M02 cut at 0.27
    cuts.AddCutMergedCalo("00010113","1111111067032200000","1111111067022700002","0163300700000000"); // INT7
    cuts.AddCutMergedCalo("00052113","1111111067032200000","1111111067022700002","0163300700000000"); // EMC7
    cuts.AddCutMergedCalo("00081113","1111111067032200000","1111111067022700002","0163300700000000"); // EGA
  } else if (trainConfig == 164){  // EMCAL clusters, different triggers with NonLinearity track matching to cluster
    cuts.AddCutMergedCalo("00010113","1111111067032200000","1111111067022210002","0163302200000000"); // INT7
    cuts.AddCutMergedCalo("00052113","1111111067032200000","1111111067022210002","0163302200000000"); // EMC7
    cuts.AddCutMergedCalo("00081113","1111111067032200000","1111111067022210002","0163302200000000"); // EMCEGA,

    // use EOverP veto in track matching
    // INT7
  } else if (trainConfig == 171){  // default TM with different E/p cuts
    cuts.AddCutMergedCalo("00010113","111111106c032200000","111111106c022700001","0163300000000000"); // no E/p cut, but reference histograms for variations
    cuts.AddCutMergedCalo("00010113","111111106d032200000","111111106d022700001","0163300000000000"); // loosest cut: EOverPMax= 3.0
    cuts.AddCutMergedCalo("00010113","111111106e032200000","111111106e022700001","0163300000000000"); // EOverPMax= 2.0
  } else if (trainConfig == 172){  // default TM with different E/p cuts
    cuts.AddCutMergedCalo("00010113","111111106f032200000","111111106f022700001","0163300000000000"); // EOverPMax= 1.75
    cuts.AddCutMergedCalo("00010113","111111106g032200000","111111106g022700001","0163300000000000"); // EOverPMax= 1.5
    cuts.AddCutMergedCalo("00010113","111111106h032200000","111111106h022700001","0163300000000000"); // hardest cut: EOverPMax= 1.25
    // EMC7
  } else if (trainConfig == 174){  // default TM with different E/p cuts
    cuts.AddCutMergedCalo("00052113","111111106c032200000","111111106c022700001","0163300000000000"); // no E/p cut, but reference histograms for variations
    cuts.AddCutMergedCalo("00052113","111111106d032200000","111111106d022700001","0163300000000000"); // loosest cut: EOverPMax= 3.0
    cuts.AddCutMergedCalo("00052113","111111106e032200000","111111106e022700001","0163300000000000"); // EOverPMax= 2.0
  } else if (trainConfig == 175){  // default TM with different E/p cuts
    cuts.AddCutMergedCalo("00052113","111111106f032200000","111111106f022700001","0163300000000000"); // EOverPMax= 1.75
    cuts.AddCutMergedCalo("00052113","111111106g032200000","111111106g022700001","0163300000000000"); // EOverPMax= 1.5
    cuts.AddCutMergedCalo("00052113","111111106h032200000","111111106h022700001","0163300000000000"); // hardest cut: EOverPMax= 1.25
    // EGA
  } else if (trainConfig == 177){  // default TM with different E/p cuts
    cuts.AddCutMergedCalo("00081113","111111106c032200000","111111106c022700001","0163300000000000"); // no E/p cut, but reference histograms for variations
    cuts.AddCutMergedCalo("00081113","111111106d032200000","111111106d022700001","0163300000000000"); // loosest cut: EOverPMax= 3.0
    cuts.AddCutMergedCalo("00081113","111111106e032200000","111111106e022700001","0163300000000000"); // EOverPMax= 2.0
  } else if (trainConfig == 178){  // default TM with different E/p cuts
    cuts.AddCutMergedCalo("00081113","111111106f032200000","111111106f022700001","0163300000000000"); // EOverPMax= 1.75
    cuts.AddCutMergedCalo("00081113","111111106g032200000","111111106g022700001","0163300000000000"); // EOverPMax= 1.5
    cuts.AddCutMergedCalo("00081113","111111106h032200000","111111106h022700001","0163300000000000"); // hardest cut: EOverPMax= 1.25

  // T0AND cuts
  } else if (trainConfig == 181){  // EMCAL clusters, different triggers with NonLinearity track matching to cluster
    cuts.AddCutMergedCalo("00011113","1111111067032200000","1111111067022110001","0163301100000000"); // INT8
    cuts.AddCutMergedCalo("00053113","1111111067032200000","1111111067022110001","0163301100000000"); // EMC8
    cuts.AddCutMergedCalo("00082113","1111111067032200000","1111111067022110001","0163301100000000"); // EMC8EGA,
  } else if (trainConfig == 182){  // EMCAL clusters, different triggers with NonLinearity track matching to cluster
    cuts.AddCutMergedCalo("00011113","1111111067032200000","1111111067022110002","0163302200000000"); // INT8
    cuts.AddCutMergedCalo("00053113","1111111067032200000","1111111067022110002","0163302200000000"); // EMC8
    cuts.AddCutMergedCalo("00082113","1111111067032200000","1111111067022110002","0163302200000000"); // EMC8EGA,

  // standard cuts (with EOverP reference plots) to test effect of new track matching modes (cf. AliCaloTrackMatcher.h)
  } else if (trainConfig == 183){ // fRunningMode = 5
    cuts.AddCutMergedCalo("00010113","111111106c032200000","111111106c022700001","0163300000000000"); // INT7
    cuts.AddCutMergedCalo("00052113","111111106c032200000","111111106c022700001","0163300000000000"); // EMC7
    cuts.AddCutMergedCalo("00081113","111111106c032200000","111111106c022700001","0163300000000000"); // EGA
  } else if (trainConfig == 184){ // fRunningMode = 6
    cuts.AddCutMergedCalo("00010113","111111106c032200000","111111106c022700001","0163300000000000"); // INT7
    cuts.AddCutMergedCalo("00052113","111111106c032200000","111111106c022700001","0163300000000000"); // EMC7
    cuts.AddCutMergedCalo("00081113","111111106c032200000","111111106c022700001","0163300000000000"); // EGA


    // shape definition study (variation of w0): using std cuts + TM off + strong M02 variation
  } else if (trainConfig == 185){
    cuts.AddCutMergedCalo("00081113","1111111060032200000","1111111060022900001","0163300000000000"); // min M02 = 0.10 + min M02 = 0.10
    cuts.AddCutMergedCalo("00081113","1111111060032200000","1111111060022700001","0163300000000000"); // min M02 = 0.10 + min M02 = 0.27
    cuts.AddCutMergedCalo("00081113","1111111060032200000","1111111060022600001","0163300000000000"); // min M02 = 0.10 + min M02 = 0.30
  } else if (trainConfig == 186){
    cuts.AddCutMergedCalo("00081113","1111111060032200000","1111111060022d00001","0163300000000000"); // min M02 = 0.10 + min M02 = 0.33
    cuts.AddCutMergedCalo("00081113","1111111060032200000","1111111060022e00001","0163300000000000"); // min M02 = 0.10 + min M02 = 0.36
    cuts.AddCutMergedCalo("00081113","1111111060032200000","1111111060022f00001","0163300000000000"); // min M02 = 0.10 + min M02 = 0.39


    // standard cuts for cell time variation
  } else if (trainConfig == 187){ // 200 ns
    cuts.AddCutMergedCalo("00010113","1111111067032200000","1111111067022700001","0163300000000000"); // INT7
    cuts.AddCutMergedCalo("00052113","1111111067032200000","1111111067022700001","0163300000000000"); // EMC7
    cuts.AddCutMergedCalo("00081113","1111111067032200000","1111111067022700001","0163300000000000"); // EGA
  } else if (trainConfig == 188){ // 100 ns
    cuts.AddCutMergedCalo("00010113","1111111067032200000","1111111067022700001","0163300000000000"); // INT7
    cuts.AddCutMergedCalo("00052113","1111111067032200000","1111111067022700001","0163300000000000"); // EMC7
    cuts.AddCutMergedCalo("00081113","1111111067032200000","1111111067022700001","0163300000000000"); // EGA
    // minimum cell aggregation energy variation
  } else if (trainConfig == 189){ // cell min energy 75 MeV
    cuts.AddCutMergedCalo("00010113","1111111067032200000","1111111067022700001","0163300000000000"); // INT7
    cuts.AddCutMergedCalo("00052113","1111111067032200000","1111111067022700001","0163300000000000"); // EMC7
    cuts.AddCutMergedCalo("00081113","1111111067032200000","1111111067022700001","0163300000000000"); // EGA
  } else if (trainConfig == 190){ // cell min energy 125 MeV
    cuts.AddCutMergedCalo("00010113","1111111067032200000","1111111067022700001","0163300000000000"); // INT7
    cuts.AddCutMergedCalo("00052113","1111111067032200000","1111111067022700001","0163300000000000"); // EMC7
    cuts.AddCutMergedCalo("00081113","1111111067032200000","1111111067022700001","0163300000000000"); // EGA

    // multiple standard cuts for supermodule-wise analysis
  } else if (trainConfig == 191){  // M02 cut at 0.27
    cuts.AddCutMergedCalo("00010113","1111111067032200000","1111111067022700001","0163300000000000"); // INT7
    cuts.AddCutMergedCalo("00052113","1111111067032200000","1111111067022700001","0163300000000000"); // EMC7
    cuts.AddCutMergedCalo("00081113","1111111067032200000","1111111067022700001","0163300000000000"); // EGA
  } else if (trainConfig == 192){  // M02 cut at 0.27
    cuts.AddCutMergedCalo("00010113","1111111067032200000","1111111067022700001","0163300000000000"); // INT7
    cuts.AddCutMergedCalo("00052113","1111111067032200000","1111111067022700001","0163300000000000"); // EMC7
    cuts.AddCutMergedCalo("00081113","1111111067032200000","1111111067022700001","0163300000000000"); // EGA
  } else if (trainConfig == 193){  // M02 cut at 0.27
    cuts.AddCutMergedCalo("00010113","1111111067032200000","1111111067022700001","0163300000000000"); // INT7
    cuts.AddCutMergedCalo("00052113","1111111067032200000","1111111067022700001","0163300000000000"); // EMC7
    cuts.AddCutMergedCalo("00081113","1111111067032200000","1111111067022700001","0163300000000000"); // EGA
  } else if (trainConfig == 194){  // M02 cut at 0.27
    cuts.AddCutMergedCalo("00010113","1111111067032200000","1111111067022700001","0163300000000000"); // INT7
    cuts.AddCutMergedCalo("00052113","1111111067032200000","1111111067022700001","0163300000000000"); // EMC7
    cuts.AddCutMergedCalo("00081113","1111111067032200000","1111111067022700001","0163300000000000"); // EGA
  } else if (trainConfig == 195){  // M02 cut at 0.27
    cuts.AddCutMergedCalo("00010113","1111111067032200000","1111111067022700001","0163300000000000"); // INT7
    cuts.AddCutMergedCalo("00052113","1111111067032200000","1111111067022700001","0163300000000000"); // EMC7
    cuts.AddCutMergedCalo("00081113","1111111067032200000","1111111067022700001","0163300000000000"); // EGA
  } else if (trainConfig == 196){  // M02 cut at 0.27
    cuts.AddCutMergedCalo("00010113","1111111067032200000","1111111067022700001","0163300000000000"); // INT7
    cuts.AddCutMergedCalo("00052113","1111111067032200000","1111111067022700001","0163300000000000"); // EMC7
    cuts.AddCutMergedCalo("00081113","1111111067032200000","1111111067022700001","0163300000000000"); // EGA
  } else if (trainConfig == 197){  // M02 cut at 0.27
    cuts.AddCutMergedCalo("00010113","1111111067032200000","1111111067022700001","0163300000000000"); // INT7
    cuts.AddCutMergedCalo("00052113","1111111067032200000","1111111067022700001","0163300000000000"); // EMC7
    cuts.AddCutMergedCalo("00081113","1111111067032200000","1111111067022700001","0163300000000000"); // EGA
  } else if (trainConfig == 198){  // M02 cut at 0.27
    cuts.AddCutMergedCalo("00010113","1111111067032200000","1111111067022700001","0163300000000000"); // INT7
    cuts.AddCutMergedCalo("00052113","1111111067032200000","1111111067022700001","0163300000000000"); // EMC7
    cuts.AddCutMergedCalo("00081113","1111111067032200000","1111111067022700001","0163300000000000"); // EGA
  } else if (trainConfig == 199){  // M02 cut at 0.27
    cuts.AddCutMergedCalo("00010113","1111111067032200000","1111111067022700001","0163300000000000"); // INT7
    cuts.AddCutMergedCalo("00052113","1111111067032200000","1111111067022700001","0163300000000000"); // EMC7
    cuts.AddCutMergedCalo("00081113","1111111067032200000","1111111067022700001","0163300000000000"); // EGA
  } else if (trainConfig == 200){  // M02 cut at 0.27
    cuts.AddCutMergedCalo("00010113","1111111067032200000","1111111067022700001","0163300000000000"); // INT7
    cuts.AddCutMergedCalo("00052113","1111111067032200000","1111111067022700001","0163300000000000"); // EMC7
    cuts.AddCutMergedCalo("00081113","1111111067032200000","1111111067022700001","0163300000000000"); // EGA


    // EOverP cuts using the "old" TrackMatching now called fRunningMode = 5 (matching with all tracks) (copied from configs 177,178)
  } else if (trainConfig == 201){  // fRunningMode = 5 with different E/p cuts
    cuts.AddCutMergedCalo("00081113","111111106c032200000","111111106c022700001","0163300000000000"); // no E/p cut, but reference histograms for variations
    cuts.AddCutMergedCalo("00081113","111111106d032200000","111111106d022700001","0163300000000000"); // loosest cut: EOverPMax= 3.0
    cuts.AddCutMergedCalo("00081113","111111106e032200000","111111106e022700001","0163300000000000"); // EOverPMax= 2.0
  } else if (trainConfig == 202){  // fRunningMode = 5 with different E/p cuts
    cuts.AddCutMergedCalo("00081113","111111106f032200000","111111106f022700001","0163300000000000"); // EOverPMax= 1.75
    cuts.AddCutMergedCalo("00081113","111111106g032200000","111111106g022700001","0163300000000000"); // EOverPMax= 1.5
    cuts.AddCutMergedCalo("00081113","111111106h032200000","111111106h022700001","0163300000000000"); // hardest cut: EOverPMax= 1.25

  } else if (trainConfig == 210){  // new standard
    cuts.AddCutMergedCalo("00010113","111111106f032200000","111111106f022700001","0163300000000000"); // INT7
    cuts.AddCutMergedCalo("00052113","111111106f032200000","111111106f022700001","0163300000000000"); // EMC7
    cuts.AddCutMergedCalo("00081113","111111106f032200000","111111106f022700001","0163300000000000"); // EGA
  } else if (trainConfig == 211){  // new standard with TB NL
    cuts.AddCutMergedCalo("00010113","1111101060032200000","1111101060022700001","0163300000000000"); // INT7
    cuts.AddCutMergedCalo("00052113","1111101060032200000","1111101060022700001","0163300000000000"); // EMC7
    cuts.AddCutMergedCalo("00081113","1111101060032200000","1111101060022700001","0163300000000000"); // EGA
  } else if (trainConfig == 212){  // new standard with exotic cut 0.97
    cuts.AddCutMergedCalo("00010113","111111106f532200000","111111106f522700001","0163300000000000"); // INT7
    cuts.AddCutMergedCalo("00052113","111111106f532200000","111111106f522700001","0163300000000000"); // EMC7
    cuts.AddCutMergedCalo("00081113","111111106f532200000","111111106f522700001","0163300000000000"); // EGA
  } else if (trainConfig == 213){  // new standard with exotic cut 0.97 and open M02>0.1
    cuts.AddCutMergedCalo("00010113","111111106f532200000","111111106f522900001","0163300000000000"); // INT7
    cuts.AddCutMergedCalo("00052113","111111106f532200000","111111106f522900001","0163300000000000"); // EMC7
    cuts.AddCutMergedCalo("00081113","111111106f532200000","111111106f522900001","0163300000000000"); // EGA
  } else if (trainConfig == 214){  // new standard with V1
    cuts.AddCutMergedCalo("00010113","111111106f032200000","111111106f022700002","0163300000000000"); // INT7
    cuts.AddCutMergedCalo("00052113","111111106f032200000","111111106f022700002","0163300000000000"); // EMC7
    cuts.AddCutMergedCalo("00081113","111111106f032200000","111111106f022700002","0163300000000000"); // EGA

  } else if (trainConfig == 215){  // new TB
    cuts.AddCutMergedCalo("00010113","111110106f032200000","111110106f022700001","0163300000000000"); // INT7
    cuts.AddCutMergedCalo("00052113","111110106f032200000","111110106f022700001","0163300000000000"); // EMC7
    cuts.AddCutMergedCalo("00081113","111110106f032200000","111110106f022700001","0163300000000000"); // EGA

  } else if (trainConfig == 216){  // new TB+finetuning CCRF
    cuts.AddCutMergedCalo("00010113","111113106f032200000","111113106f022700001","0163300000000000"); // INT7
    cuts.AddCutMergedCalo("00052113","111113106f032200000","111113106f022700001","0163300000000000"); // EMC7
    cuts.AddCutMergedCalo("00081113","111113106f032200000","111113106f022700001","0163300000000000"); // EGA
  } else if (trainConfig == 217){  // new TB+finetuning CRF
    cuts.AddCutMergedCalo("00010113","111113206f032200000","111113206f022700001","0163300000000000"); // INT7
    cuts.AddCutMergedCalo("00052113","111113206f032200000","111113206f022700001","0163300000000000"); // EMC7
    cuts.AddCutMergedCalo("00081113","111113206f032200000","111113206f022700001","0163300000000000"); // EGA
  } else if (trainConfig == 218){  // new TB+finetuning CCMF
    cuts.AddCutMergedCalo("00010113","111113306f032200000","111113306f022700001","0163300000000000"); // INT7
    cuts.AddCutMergedCalo("00052113","111113306f032200000","111113306f022700001","0163300000000000"); // EMC7
    cuts.AddCutMergedCalo("00081113","111113306f032200000","111113306f022700001","0163300000000000"); // EGA
  } else if (trainConfig == 219){  // new TB+finetuning CMF
    cuts.AddCutMergedCalo("00010113","111113406f032200000","111113406f022700001","0163300000000000"); // INT7
    cuts.AddCutMergedCalo("00052113","111113406f032200000","111113406f022700001","0163300000000000"); // EMC7
    cuts.AddCutMergedCalo("00081113","111113406f032200000","111113406f022700001","0163300000000000"); // EGA
  // systematics pp 8 TeV
  } else if (trainConfig == 220){ // M02 var 1
    cuts.AddCutMergedCalo("00081113","111111106f032200000","111111106f022600001","0163300000000000"); // min 0.3
    cuts.AddCutMergedCalo("00081113","111111106f032200000","111111106f022800001","0163300000000000"); // min 0.25
  } else if (trainConfig == 221){ // M02 var 2
    cuts.AddCutMergedCalo("00081113","111111106f032200000","111111106f022900001","0163300000000000"); // min 0.1
    cuts.AddCutMergedCalo("00081113","111111106f032200000","111111106f022a00001","0163300000000000"); // min 0.26
  } else if (trainConfig == 222){ // M02 var 3
    cuts.AddCutMergedCalo("00081113","111111106f032200000","111111106f022d00001","0163300000000000"); // min 0.33
    cuts.AddCutMergedCalo("00081113","111111106f032200000","111111106f022e00001","0163300000000000"); // min 0.36
  } else if (trainConfig == 223){  // NL variations 1
    cuts.AddCutMergedCalo("00081113","111112206f032200000","111112206f022700001","0163300000000000"); // SDM loose time
    cuts.AddCutMergedCalo("00081113","111111206f032200000","111111206f022700001","0163300000000000"); // conv calo tight time
  } else if (trainConfig == 224){ // NL variations 2
    cuts.AddCutMergedCalo("00081113","111116506f032200000","111116506f022700001","0163300000000000"); // tb
    cuts.AddCutMergedCalo("00081113","111110006f032200000","111110006f022700001","0163300000000000"); // none
  } else if (trainConfig == 225){  // TRD material
    cuts.AddCutMergedCalo("00081113","111211106f032200000","111211106f022700001","0163300000000000"); // TRD infront
    cuts.AddCutMergedCalo("00081113","111131106f032200000","111131106f022700001","0163300000000000"); // no TRD infront
  } else if (trainConfig == 226){ // varied exoctics 1
    cuts.AddCutMergedCalo("00081113","111111106f232200000","111111106f222700001","0163300000000000"); // frac = 0.99
    cuts.AddCutMergedCalo("00081113","111111106f332200000","111111106f322700001","0163300000000000"); // frac = 0.98
  } else if (trainConfig == 227){ // varied exoctics 2
    cuts.AddCutMergedCalo("00081113","111111106f532200000","111111106f522700001","0163300000000000"); // frac = 0.97
    cuts.AddCutMergedCalo("00081113","111111106f732200000","111111106f722700001","0163300000000000"); // frac = 0.96
  } else if (trainConfig == 228){  // variation of mass and alpha cut
    cuts.AddCutMergedCalo("00081113","111111106f032200000","111111106f022700001","0163301700000000"); // only band at 0
    cuts.AddCutMergedCalo("00081113","111111106f032200000","111111106f022700001","0163300700000000"); // only band at 0, no alpha
  } else if (trainConfig == 229){ // mass cut
    cuts.AddCutMergedCalo("00081113","111111106f032200000","111111106f022700001","0163301000000000"); // no mass cut
    cuts.AddCutMergedCalo("00081113","111111106f032200000","111111106f022700001","0163300100000000"); // no alphat
  } else if (trainConfig == 230){  // variation track matching to cluster & mass variations new defaults
    cuts.AddCutMergedCalo("00081113","1111111060032200000","1111111060022700001","0163300000000000"); // no TM
    cuts.AddCutMergedCalo("00081113","1111111061032200000","1111111061022700001","0163300000000000"); // looser TM
  } else if (trainConfig == 231){ // variation track matching to cluster & mass variations new defaults
    cuts.AddCutMergedCalo("00081113","1111111066032200000","1111111066022700001","0163300000000000"); // tighter TM
    cuts.AddCutMergedCalo("00081113","1111111067032200000","1111111067022700001","0163300000000000"); // EMC TM
  } else if (trainConfig == 232){  // dist to bad channel
    cuts.AddCutMergedCalo("00081113","111111106f032200000","111111116f022700001","0163300000000000"); // 1 <= coll+row
    cuts.AddCutMergedCalo("00081113","111111106f032200000","111111126f022700001","0163300000000000"); // 2 <= coll+row
  } else if (trainConfig == 233){ // dist to bad channel
    cuts.AddCutMergedCalo("00081113","111111106f032200000","111111156f022700001","0163300000000000"); // 1 <= row || coll || 0.5*(coll+row)
    cuts.AddCutMergedCalo("00081113","111111106f032200000","111111166f022700001","0163300000000000"); // 2 <= row || coll || 0.5*(coll+row)
  } else if (trainConfig == 234){  // timing cuts
    cuts.AddCutMergedCalo("00081113","111111109f032200000","111111109f022700001","0163300000000000"); // -20 to 25
    cuts.AddCutMergedCalo("00081113","111111105f032200000","111111105f022700001","0163300000000000"); // -50 to 50
  } else if (trainConfig == 235){  // timing cuts 2
    cuts.AddCutMergedCalo("00081113","111111103f032200000","111111103f022700001","0163300000000000"); // -200 to 200
    cuts.AddCutMergedCalo("00081113","111111102f032200000","111111102f022700001","0163300000000000"); // -500 to 500
  } else if (trainConfig == 236){  // min E cuts
    cuts.AddCutMergedCalo("00081113","111111106f032200000","111111106f042700001","0163300000000000"); // E>7GeV
    cuts.AddCutMergedCalo("00081113","111111106f032200000","111111106f092700001","0163300000000000"); // E>9.5GeV
  } else if (trainConfig == 237){  // TRD material
    cuts.AddCutMergedCalo("00081113","111211106f032200000","111211106f022700001","0163300000000000"); // TRD infront
    cuts.AddCutMergedCalo("00081113","111131106f032200000","111131106f022700001","0163300000000000"); // no TRD infront
  } else if (trainConfig == 238){ // variation track matching to cluster & mass variations new defaults
    cuts.AddCutMergedCalo("00081113","111111106d032200000","111111106d022700001","0163300000000000"); // fEOverPMax = 3.0
    cuts.AddCutMergedCalo("00081113","111111106e032200000","111111106e022700001","0163300000000000"); // fEOverPMax = 2.0
  } else if (trainConfig == 239){ // variation track matching to cluster & mass variations new defaults
    cuts.AddCutMergedCalo("00081113","111111106g032200000","111111106g022700001","0163300000000000"); // fEOverPMax = 1.5
    cuts.AddCutMergedCalo("00081113","111111106h032200000","111111106h022700001","0163300000000000"); // fEOverPMax = 1.25

  } else if (trainConfig == 250){  // EMCAL clusters 7 TeV LHC11 TM on
    cuts.AddCutMergedCalo("00010113","1111100067032200000","1111100067022700001","0163300000000000"); // INT7
    cuts.AddCutMergedCalo("00052113","1111100067032200000","1111100067022700001","0163300000000000"); // EMC7
  } else if (trainConfig == 251){  // EMCAL clusters 7 TeV LHC11 TM off
    cuts.AddCutMergedCalo("00010113","1111100067032200000","1111100060022700001","0163300000000000"); // INT7
    cuts.AddCutMergedCalo("00052113","1111100067032200000","1111100060022700001","0163300000000000"); // EMC7
  } else if (trainConfig == 252){  // EMCAL clusters 7 TeV LHC11 TM on + TB NonLin
    cuts.AddCutMergedCalo("00010113","1111102067032200000","1111102067022700001","0163300000000000"); // INT7
    cuts.AddCutMergedCalo("00052113","1111102067032200000","1111102067022700001","0163300000000000"); // EMC7
  } else if (trainConfig == 253){  // EMCAL clusters 7 TeV LHC11 TM off + TB NonLin
    cuts.AddCutMergedCalo("00010113","1111102067032200000","1111102060022700001","0163300000000000"); // INT7
    cuts.AddCutMergedCalo("00052113","1111102067032200000","1111102060022700001","0163300000000000"); // EMC7

  } else if (trainConfig == 260){  // EMCAL clusters 7 TeV LHC11 TM on + PCMEDMC NonLin
    cuts.AddCutMergedCalo("00010113","1111111067032200000","1111111067022700001","0163300000000000"); // INT7
    cuts.AddCutMergedCalo("00010113","1111121067032200000","1111121067022700001","0163300000000000"); // INT7
    cuts.AddCutMergedCalo("00052113","1111111067032200000","111111106f022700001","0163300000000000"); // EMC7
    cuts.AddCutMergedCalo("00052113","1111121067032200000","111112106f022700001","0163300000000000"); // EMC7
  } else if (trainConfig == 261){  // EMCAL clusters 7 TeV LHC10 TM on + PCMEDMC NonLin
    cuts.AddCutMergedCalo("00000113","1111111067032200000","1111111067022700001","0163300000000000"); // INT7
    cuts.AddCutMergedCalo("00000113","1111121067032200000","1111121067022700001","0163300000000000"); // INT7

  } else if (trainConfig == 262){  // EMCAL clusters 7 TeV LHC11 TM on + TB NonLin
    cuts.AddCutMergedCalo("00010113","111110106f032200000","111110106f022700001","0163300000000000"); // INT7
    cuts.AddCutMergedCalo("00052113","111110106f032200000","111110106f022700001","0163300000000000"); // EMC7
  } else if (trainConfig == 263){  // EMCAL clusters 7 TeV LHC11 TM on + TB NonLin
    cuts.AddCutMergedCalo("00010113","111113106f032200000","111113106f022700001","0163300000000000"); // INT7
    cuts.AddCutMergedCalo("00052113","111113106f032200000","111113106f022700001","0163300000000000"); // EMC7
  } else if (trainConfig == 264){  // EMCAL clusters 7 TeV LHC11 TM on + TB NonLin
    cuts.AddCutMergedCalo("00010113","111113206f032200000","111113206f022700001","0163300000000000"); // INT7
    cuts.AddCutMergedCalo("00052113","111113206f032200000","111113206f022700001","0163300000000000"); // EMC7
  } else if (trainConfig == 265){  // EMCAL clusters 7 TeV LHC11 TM on + TB NonLin
    cuts.AddCutMergedCalo("00010113","111113306f032200000","111113306f022700001","0163300000000000"); // INT7
    cuts.AddCutMergedCalo("00052113","111113306f032200000","111113306f022700001","0163300000000000"); // EMC7
  } else if (trainConfig == 266){  // EMCAL clusters 7 TeV LHC11 TM on + TB NonLin
    cuts.AddCutMergedCalo("00010113","111113406f032200000","111113406f022700001","0163300000000000"); // INT7
    cuts.AddCutMergedCalo("00052113","111113406f032200000","111113406f022700001","0163300000000000"); // EMC7
    // 13 TeV & 5 TeV
  } else if (trainConfig == 401){ // pp 2.76TeV paper cuts : open timing, TB nonlin
    cuts.AddCutMergedCalo("00010113","1111101017032200000","1111101017022700001","0163300000000000"); // INT7
    cuts.AddCutMergedCalo("00052113","1111101017032200000","1111101017022700001","0163300000000000"); // EMC7
    cuts.AddCutMergedCalo("00085113","1111101017032200000","1111101017022700001","0163300000000000"); // EG2
    cuts.AddCutMergedCalo("00083113","1111101017032200000","1111101017022700001","0163300000000000"); // EG1
  } else if (trainConfig == 402){  // pp 2.76TeV paper cuts:  w/o mass, open timing, TB nonlin
    cuts.AddCutMergedCalo("00010113","1111101017032200000","1111101017022000001","0163300000000000");
    cuts.AddCutMergedCalo("00052113","1111101017032200000","1111101017022000001","0163300000000000"); // EMC7
    cuts.AddCutMergedCalo("00085113","1111101017032200000","1111101017022000001","0163300000000000"); // EG2
    cuts.AddCutMergedCalo("00083113","1111101017032200000","1111101017022000001","0163300000000000"); // EG1
  } else if (trainConfig == 403){ // pp 2.76TeV paper cuts : open timing, TB nonlin
    cuts.AddCutMergedCalo("00010113","1111100017032200000","1111100017022700001","0163300000000000"); // INT7
    cuts.AddCutMergedCalo("00052113","1111100017032200000","1111100017022700001","0163300000000000"); // EMC7
    cuts.AddCutMergedCalo("00085113","1111100017032200000","1111100017022700001","0163300000000000"); // EG2
    cuts.AddCutMergedCalo("00083113","1111100017032200000","1111100017022700001","0163300000000000"); // EG1
  } else if (trainConfig == 404){  // pp 2.76TeV paper cuts:  w/o mass, open timing, TB nonlin
    cuts.AddCutMergedCalo("00010113","1111100017032200000","1111100017022000001","0163300000000000");
    cuts.AddCutMergedCalo("00052113","1111100017032200000","1111100017022000001","0163300000000000"); // EMC7
    cuts.AddCutMergedCalo("00085113","1111100017032200000","1111100017022000001","0163300000000000"); // EG2
    cuts.AddCutMergedCalo("00083113","1111100017032200000","1111100017022000001","0163300000000000"); // EG1


  } else if (trainConfig == 410){ // no nonlin
    cuts.AddCutMergedCalo("00010113","1111100067032200000","1111100067022700001","0163300000000000"); // INT7
    cuts.AddCutMergedCalo("00052113","1111100067032200000","1111100067022700001","0163300000000000"); // EMC7
    cuts.AddCutMergedCalo("00085113","1111100067032200000","1111100067022700001","0163300000000000"); // EG2
    cuts.AddCutMergedCalo("00083113","1111100067032200000","1111100067022700001","0163300000000000"); // EG1
  } else if (trainConfig == 411){ // standard kPi0MCv5 for MC and kSDMv5 for data from Jason
    cuts.AddCutMergedCalo("00010113","1111101067032200000","1111101067022700001","0163300000000000"); // INT7
    cuts.AddCutMergedCalo("00052113","1111101067032200000","1111101067022700001","0163300000000000"); // EMC7
    cuts.AddCutMergedCalo("00085113","1111101067032200000","1111101067022700001","0163300000000000"); // EG2
    cuts.AddCutMergedCalo("00083113","1111101067032200000","1111101067022700001","0163300000000000"); // EG1
  } else if (trainConfig == 412){ // NonLinearity pp Calo - only shifting MC
    cuts.AddCutMergedCalo("00010113","1111112067032200000","1111112067022700001","0163300000000000"); // INT7
    cuts.AddCutMergedCalo("00052113","1111112067032200000","1111112067022700001","0163300000000000"); // EMC7
    cuts.AddCutMergedCalo("00085113","1111112067032200000","1111112067022700001","0163300000000000"); // EG2
    cuts.AddCutMergedCalo("00083113","1111112067032200000","1111112067022700001","0163300000000000"); // EG1
  } else if (trainConfig == 413){ // nonlin vars MB
    cuts.AddCutMergedCalo("00010113","1111111067032200000","1111111067022700001","0163300000000000"); // INT7
    cuts.AddCutMergedCalo("00010113","1111112067032200000","1111112067022700001","0163300000000000"); // INT7
    cuts.AddCutMergedCalo("00010113","1111121067032200000","1111121067022700001","0163300000000000"); // INT7
    cuts.AddCutMergedCalo("00010113","1111122067032200000","1111122067022700001","0163300000000000"); // INT7
  } else if (trainConfig == 414){ // nonlin vars EMC7
    cuts.AddCutMergedCalo("00052113","1111111067032200000","1111111067022700001","0163300000000000"); // EMC7
    cuts.AddCutMergedCalo("00052113","1111112067032200000","1111112067022700001","0163300000000000"); // EMC7
    cuts.AddCutMergedCalo("00052113","1111121067032200000","1111121067022700001","0163300000000000"); // EMC7
    cuts.AddCutMergedCalo("00052113","1111122067032200000","1111122067022700001","0163300000000000"); // EMC7
  } else if (trainConfig == 415){ // nonlin vars EG2
    cuts.AddCutMergedCalo("00085113","1111111067032200000","1111111067022700001","0163300000000000"); // EG2
    cuts.AddCutMergedCalo("00085113","1111112067032200000","1111112067022700001","0163300000000000"); // EG2
    cuts.AddCutMergedCalo("00085113","1111121067032200000","1111121067022700001","0163300000000000"); // EG2
    cuts.AddCutMergedCalo("00085113","1111122067032200000","1111122067022700001","0163300000000000"); // EG2

  } else if (trainConfig == 460){ // INT7 EMCAL standard cut but with E/p TM veto
    cuts.AddCutMergedCalo("00010113","1111111067032200000","1111111067022700001","0163300000000000"); // std INT7
    cuts.AddCutMergedCalo("00010113","111111106c032200000","111111106c022700001","0163300000000000"); // fEOverPMax = 9e9
    cuts.AddCutMergedCalo("00010113","111111106d032200000","111111106d022700001","0163300000000000"); // fEOverPMax = 3.0
    cuts.AddCutMergedCalo("00010113","111111106e032200000","111111106e022700001","0163300000000000"); // fEOverPMax = 2.0
  } else if (trainConfig == 461){ // INT7 EMCAL standard cut but with E/p TM veto
    cuts.AddCutMergedCalo("00010113","111111106f032200000","111111106f022700001","0163300000000000"); // fEOverPMax = 1.75
    cuts.AddCutMergedCalo("00010113","111111106g032200000","111111106g022700001","0163300000000000"); // fEOverPMax = 1.5
    cuts.AddCutMergedCalo("00010113","111111106h032200000","111111106h022700001","0163300000000000"); // fEOverPMax = 1.25
  } else if (trainConfig == 462){ // EMC7 EMCAL standard cut but with E/p TM veto
    cuts.AddCutMergedCalo("00052113","1111111067032200000","1111111067022700001","0163300000000000"); // std EMC7
    cuts.AddCutMergedCalo("00052113","111111106c032200000","111111106c022700001","0163300000000000"); // fEOverPMax = 9e9
    cuts.AddCutMergedCalo("00052113","111111106d032200000","111111106d022700001","0163300000000000"); // fEOverPMax = 3.0
    cuts.AddCutMergedCalo("00052113","111111106e032200000","111111106e022700001","0163300000000000"); // fEOverPMax = 2.0
  } else if (trainConfig == 463){ // EMC7 EMCAL standard cut but with E/p TM veto
    cuts.AddCutMergedCalo("00052113","111111106f032200000","111111106f022700001","0163300000000000"); // fEOverPMax = 1.75
    cuts.AddCutMergedCalo("00052113","111111106g032200000","111111106g022700001","0163300000000000"); // fEOverPMax = 1.5
    cuts.AddCutMergedCalo("00052113","111111106h032200000","111111106h022700001","0163300000000000"); // fEOverPMax = 1.25
  } else if (trainConfig == 464){ // EG1 EMCAL standard cut but with E/p TM veto
    cuts.AddCutMergedCalo("00083113","1111111067032200000","1111111067022700001","0163300000000000"); // std EG1
    cuts.AddCutMergedCalo("00083113","111111106c032200000","111111106c022700001","0163300000000000"); // fEOverPMax = 9e9
    cuts.AddCutMergedCalo("00083113","111111106d032200000","111111106d022700001","0163300000000000"); // fEOverPMax = 3.0
    cuts.AddCutMergedCalo("00083113","111111106e032200000","111111106e022700001","0163300000000000"); // fEOverPMax = 2.0
  } else if (trainConfig == 465){ // EG1 EMCAL standard cut but with E/p TM veto
    cuts.AddCutMergedCalo("00083113","111111106f032200000","111111106f022700001","0163300000000000"); // fEOverPMax = 1.75
    cuts.AddCutMergedCalo("00083113","111111106g032200000","111111106g022700001","0163300000000000"); // fEOverPMax = 1.5
    cuts.AddCutMergedCalo("00083113","111111106h032200000","111111106h022700001","0163300000000000"); // fEOverPMax = 1.25
  } else if (trainConfig == 466){ // EG2 EMCAL standard cut but with E/p TM veto
    cuts.AddCutMergedCalo("00085113","1111111067032200000","1111111067022700001","0163300000000000"); // std EG2
    cuts.AddCutMergedCalo("00085113","111111106c032200000","111111106c022700001","0163300000000000"); // fEOverPMax = 9e9
    cuts.AddCutMergedCalo("00085113","111111106d032200000","111111106d022700001","0163300000000000"); // fEOverPMax = 3.0
    cuts.AddCutMergedCalo("00085113","111111106e032200000","111111106e022700001","0163300000000000"); // fEOverPMax = 2.0
  } else if (trainConfig == 467){ // EG2 EMCAL standard cut but with E/p TM veto
    cuts.AddCutMergedCalo("00085113","111111106f032200000","111111106f022700001","0163300000000000"); // fEOverPMax = 1.75
    cuts.AddCutMergedCalo("00085113","111111106g032200000","111111106g022700001","0163300000000000"); // fEOverPMax = 1.5
    cuts.AddCutMergedCalo("00085113","111111106h032200000","111111106h022700001","0163300000000000"); // fEOverPMax = 1.25
  } else if (trainConfig == 480){ // no TM - CALO+CALOFAST triggers
    cuts.AddCutMergedCalo("000a0113","1111111060032200000","1111111060022700001","0163300000000000"); // INT7
    cuts.AddCutMergedCalo("000a1113","1111111060032200000","1111111060022700001","0163300000000000"); // EMC7
    cuts.AddCutMergedCalo("000a2113","1111111060032200000","1111111060022700001","0163300000000000"); // EG2
    cuts.AddCutMergedCalo("000a3113","1111111060032200000","1111111060022700001","0163300000000000"); // EG1
  } else if (trainConfig == 481){ // no TM
    cuts.AddCutMergedCalo("00010113","1111111060032200000","1111111060022700001","0163300000000000"); // INT7
    cuts.AddCutMergedCalo("00052113","1111111060032200000","1111111060022700001","0163300000000000"); // EMC7
    cuts.AddCutMergedCalo("00085113","1111111060032200000","1111111060022700001","0163300000000000"); // EG2
    cuts.AddCutMergedCalo("00083113","1111111060032200000","1111111060022700001","0163300000000000"); // EG1
  } else if (trainConfig == 482){ // no TM - CALO+CALOFAST triggers
    cuts.AddCutMergedCalo("000a0113","4117931060032200000","4117931060022700001","0163300000000000"); // INT7
    cuts.AddCutMergedCalo("000a1113","4117931060032200000","4117931060022700001","0163300000000000"); // EMC7
    cuts.AddCutMergedCalo("000a2113","4117931060032200000","4117931060022700001","0163300000000000"); // EG2
    cuts.AddCutMergedCalo("000a3113","4117931060032200000","4117931060022700001","0163300000000000"); // EG1
  } else if (trainConfig == 483){ // no TM - CALO+CALOFAST triggers
    cuts.AddCutMergedCalo("00010113","4117931060032200000","4117931060022700001","0163300000000000"); // INT7
    cuts.AddCutMergedCalo("000a1113","4117931060032200000","4117931060022700001","0163300000000000"); // EMC7
    cuts.AddCutMergedCalo("000a2113","4117931060032200000","4117931060022700001","0163300000000000"); // EG2
  } else if (trainConfig == 484){ // no TM - CALO+CALOFAST triggers with exotics cut 0.97
    cuts.AddCutMergedCalo("00010113","4117931060532200000","4117931060522700001","0163300000000000"); // INT7
    cuts.AddCutMergedCalo("000a1113","4117931060532200000","4117931060522700001","0163300000000000"); // EMC7
    cuts.AddCutMergedCalo("000a2113","4117931060532200000","4117931060522700001","0163300000000000"); // EG2
  } else if (trainConfig == 485){ // no TM - CALO+CALOFAST triggers with exotics cut 0.97 and open M02 (>0.1)
    cuts.AddCutMergedCalo("00010113","4117931060532200000","4117931060522900001","0163300000000000"); // INT7
    cuts.AddCutMergedCalo("000a1113","4117931060532200000","4117931060522900001","0163300000000000"); // EMC7
    cuts.AddCutMergedCalo("000a2113","4117931060532200000","4117931060522900001","0163300000000000"); // EG2
  } else if (trainConfig == 486){ // no TM - CALO+CALOFAST triggers with exotics cut 0.97 and open M02 (>0.1)
    cuts.AddCutMergedCalo("00010113","4117901060532200000","4117901060522900001","0163300000000000"); // INT7
    cuts.AddCutMergedCalo("000a1113","4117901060532200000","4117901060522900001","0163300000000000"); // EMC7
    cuts.AddCutMergedCalo("000a2113","4117901060532200000","4117901060522900001","0163300000000000"); // EG2
  } else if (trainConfig == 487){ // no TM - CALO+CALOFAST triggers - TB NL
    cuts.AddCutMergedCalo("00010113","4117931060032200000","4117931060022700001","0163300000000000"); // INT7
    cuts.AddCutMergedCalo("000a1113","4117931060032200000","4117931060022700001","0163300000000000"); // EMC7
    cuts.AddCutMergedCalo("000a2113","4117931060032200000","4117931060022700001","0163300000000000"); // EG2
  } else if (trainConfig == 488){ // no TM - CALO+CALOFAST triggers - TB NL
    cuts.AddCutMergedCalo("00010113","4117932060032200000","4117932060022700001","0163300000000000"); // INT7
    cuts.AddCutMergedCalo("000a1113","4117932060032200000","4117932060022700001","0163300000000000"); // EMC7
    cuts.AddCutMergedCalo("000a2113","4117932060032200000","4117932060022700001","0163300000000000"); // EG2
  } else if (trainConfig == 489){ // no TM - CALO+CALOFAST triggers - TB NL
    cuts.AddCutMergedCalo("00010113","4117931060032200000","4117931060022700001","0163300000000000"); // INT7
    cuts.AddCutMergedCalo("000a1113","4117931060032200000","4117931060022700001","0163300000000000"); // EMC7
    cuts.AddCutMergedCalo("000a2113","4117931060032200000","4117931060022700001","0163300000000000"); // EG2
  } else if (trainConfig == 490){ // no TM - CALO+CALOFAST triggers - TB NL
    cuts.AddCutMergedCalo("00010113","4117931060032200000","4117931060022700002","0163300000000000"); // INT7
    cuts.AddCutMergedCalo("000a1113","4117931060032200000","4117931060022700002","0163300000000000"); // EMC7
    cuts.AddCutMergedCalo("000a2113","4117931060032200000","4117931060022700002","0163300000000000"); // EG2

  } else if (trainConfig == 500){ // EMCAL clusters pp 5 TeV V0M mult cuts
    cuts.AddCutMergedCalo("m0110113","4117931060032200000","4117931060022700001","0163300000000000");; // std 0-1%
    cuts.AddCutMergedCalo("m1510113","4117931060032200000","4117931060022700001","0163300000000000");; // std 1-5%
    cuts.AddCutMergedCalo("m5k10113","4117931060032200000","4117931060022700001","0163300000000000");; // std 5-20%
  } else if (trainConfig == 501){ // EMCAL clusters pp 5 TeV V0M mult cuts
    cuts.AddCutMergedCalo("n2410113","4117931060032200000","4117931060022700001","0163300000000000");; // std 20-40%
    cuts.AddCutMergedCalo("n4710113","4117931060032200000","4117931060022700001","0163300000000000");; // std 40-70%
    cuts.AddCutMergedCalo("n7a10113","4117931060032200000","4117931060022700001","0163300000000000");; // std 70-100%
  } else if (trainConfig == 505){ // EMCAL clusters pp 5 TeV SPD mult cuts
    cuts.AddCutMergedCalo("o0110113","4117931060032200000","4117931060022700001","0163300000000000");; // std 0-1%
    cuts.AddCutMergedCalo("o0210113","4117931060032200000","4117931060022700001","0163300000000000");; // std 0-2%
    cuts.AddCutMergedCalo("o0510113","4117931060032200000","4117931060022700001","0163300000000000");; // std 2-5%
  } else if (trainConfig == 506){ // EMCAL clusters pp 5 TeV SPD mult cuts
    cuts.AddCutMergedCalo("o5k10113","4117931060032200000","4117931060022700001","0163300000000000");; // std 5-20%
    cuts.AddCutMergedCalo("p2610113","4117931060032200000","4117931060022700001","0163300000000000");; // std 20-60%
    cuts.AddCutMergedCalo("p6a10113","4117931060032200000","4117931060022700001","0163300000000000");; // std 60-100%

  } else if (trainConfig == 510){ // EMCAL clusters pp 5 TeV V0M mult cuts
    cuts.AddCutMergedCalo("m01a1113","4117931060032200000","4117931060022700001","0163300000000000");; // std 0-1%
    cuts.AddCutMergedCalo("m15a1113","4117931060032200000","4117931060022700001","0163300000000000");; // std 1-5%
    cuts.AddCutMergedCalo("m5ka1113","4117931060032200000","4117931060022700001","0163300000000000");; // std 5-20%
  } else if (trainConfig == 511){ // EMCAL clusters pp 5 TeV V0M mult cuts
    cuts.AddCutMergedCalo("n24a1113","4117931060032200000","4117931060022700001","0163300000000000");; // std 20-40%
    cuts.AddCutMergedCalo("n47a1113","4117931060032200000","4117931060022700001","0163300000000000");; // std 40-70%
    cuts.AddCutMergedCalo("n7aa1113","4117931060032200000","4117931060022700001","0163300000000000");; // std 70-100%
  } else if (trainConfig == 515){ // EMCAL clusters pp 5 TeV SPD mult cuts
    cuts.AddCutMergedCalo("o01a1113","4117931060032200000","4117931060022700001","0163300000000000");; // std 0-1%
    cuts.AddCutMergedCalo("o02a1113","4117931060032200000","4117931060022700001","0163300000000000");; // std 0-2%
    cuts.AddCutMergedCalo("o05a1113","4117931060032200000","4117931060022700001","0163300000000000");; // std 2-5%
  } else if (trainConfig == 516){ // EMCAL clusters pp 5 TeV SPD mult cuts
    cuts.AddCutMergedCalo("o5ka1113","4117931060032200000","4117931060022700001","0163300000000000");; // std 5-20%
    cuts.AddCutMergedCalo("p26a1113","4117931060032200000","4117931060022700001","0163300000000000");; // std 20-60%
    cuts.AddCutMergedCalo("p6aa1113","4117931060032200000","4117931060022700001","0163300000000000");; // std 60-100%

  } else if (trainConfig == 520){ // EMCAL clusters pp 5 TeV V0M mult cuts
    cuts.AddCutMergedCalo("m01a2113","4117931060032200000","4117931060022700001","0163300000000000");; // std 0-1%
    cuts.AddCutMergedCalo("m15a2113","4117931060032200000","4117931060022700001","0163300000000000");; // std 1-5%
    cuts.AddCutMergedCalo("m5ka2113","4117931060032200000","4117931060022700001","0163300000000000");; // std 5-20%
  } else if (trainConfig == 521){ // EMCAL clusters pp 5 TeV V0M mult cuts
    cuts.AddCutMergedCalo("n24a2113","4117931060032200000","4117931060022700001","0163300000000000");; // std 20-40%
    cuts.AddCutMergedCalo("n47a2113","4117931060032200000","4117931060022700001","0163300000000000");; // std 40-70%
    cuts.AddCutMergedCalo("n7aa2113","4117931060032200000","4117931060022700001","0163300000000000");; // std 70-100%
  } else if (trainConfig == 525){ // EMCAL clusters pp 5 TeV SPD mult cuts
    cuts.AddCutMergedCalo("o01a2113","4117931060032200000","4117931060022700001","0163300000000000");; // std 0-1%
    cuts.AddCutMergedCalo("o02a2113","4117931060032200000","4117931060022700001","0163300000000000");; // std 0-2%
    cuts.AddCutMergedCalo("o05a2113","4117931060032200000","4117931060022700001","0163300000000000");; // std 2-5%
  } else if (trainConfig == 526){ // EMCAL clusters pp 5 TeV SPD mult cuts
    cuts.AddCutMergedCalo("o5ka2113","4117931060032200000","4117931060022700001","0163300000000000");; // std 5-20%
    cuts.AddCutMergedCalo("p26a2113","4117931060032200000","4117931060022700001","0163300000000000");; // std 20-60%
    cuts.AddCutMergedCalo("p6aa2113","4117931060032200000","4117931060022700001","0163300000000000");; // std 60-100%

  } else if (trainConfig == 530){ // EMCAL clusters pp 5 TeV V0M mult cuts
    cuts.AddCutMergedCalo("m0181113","4117931060032200000","4117931060022700001","0163300000000000");; // std 0-1%
    cuts.AddCutMergedCalo("m0281113","4117931060032200000","4117931060022700001","0163300000000000");; // std 0-2%
    cuts.AddCutMergedCalo("m0581113","4117931060032200000","4117931060022700001","0163300000000000");; // std 2-5%
  } else if (trainConfig == 531){ // EMCAL clusters pp 5 TeV V0M mult cuts
    cuts.AddCutMergedCalo("m5k81113","4117931060032200000","4117931060022700001","0163300000000000");; // std 5-20%
    cuts.AddCutMergedCalo("n2681113","4117931060032200000","4117931060022700001","0163300000000000");; // std 20-60%
    cuts.AddCutMergedCalo("n6a81113","4117931060032200000","4117931060022700001","0163300000000000");; // std 60-100%
  } else if (trainConfig == 535){ // EMCAL clusters pp 5 TeV SPD mult cuts
    cuts.AddCutMergedCalo("o0181113","4117931060032200000","4117931060022700001","0163300000000000");; // std 0-1%
    cuts.AddCutMergedCalo("o0281113","4117931060032200000","4117931060022700001","0163300000000000");; // std 0-2%
    cuts.AddCutMergedCalo("o0581113","4117931060032200000","4117931060022700001","0163300000000000");; // std 2-5%
  } else if (trainConfig == 536){ // EMCAL clusters pp 5 TeV SPD mult cuts
    cuts.AddCutMergedCalo("o5k81113","4117931060032200000","4117931060022700001","0163300000000000");; // std 5-20%
    cuts.AddCutMergedCalo("p2681113","4117931060032200000","4117931060022700001","0163300000000000");; // std 20-60%
    cuts.AddCutMergedCalo("p6a81113","4117931060032200000","4117931060022700001","0163300000000000");; // std 60-100%

  //Analysis cuts pp 8 TeV with TPC, no TM, +-50ns, Nico TB NL
  } else if (trainConfig == 1400){
    cuts.AddCutMergedCalo("00010113","1111165060032200000","1111165060022700001","0163300000000000"); // INT7
  } else if (trainConfig == 1401){
    cuts.AddCutMergedCalo("00052113","1111165060032200000","1111165060022700001","0163300000000000"); // EMC7
    cuts.AddCutMergedCalo("00081113","1111165060032200000","1111165060022700001","0163300000000000"); // EGA
  } else if (trainConfig == 1402){
    cuts.AddCutMergedCalo("00091113","1111165060032200000","1111165060022700001","0163300000000000"); // EJE
  //Analysis cuts pp 8 TeV with TPC, std TM, +-50ns, Nico TB NL
  } else if (trainConfig == 1410){
    cuts.AddCutMergedCalo("00010113","1111165067032200000","1111165067022700001","0163300000000000"); // INT7
  } else if (trainConfig == 1411){
    cuts.AddCutMergedCalo("00052113","1111165067032200000","1111165067022700001","0163300000000000"); // EMC7
    cuts.AddCutMergedCalo("00081113","1111165067032200000","1111165067022700001","0163300000000000"); // EGA
  } else if (trainConfig == 1412){
    cuts.AddCutMergedCalo("00091113","1111165067032200000","1111165067022700001","0163300000000000"); // EJE
  //Analysis cuts pp 8 TeV with TPC, fEOverPMax = 1.75, +-50ns, Nico TB NL
  } else if (trainConfig == 1420){
    cuts.AddCutMergedCalo("00010113","111116506f032200000","111116506f022700001","0163300000000000"); // INT7
  } else if (trainConfig == 1421){
    cuts.AddCutMergedCalo("00052113","111116506f032200000","111116506f022700001","0163300000000000"); // EMC7
    cuts.AddCutMergedCalo("00081113","111116506f032200000","111116506f022700001","0163300000000000"); // EGA
  } else if (trainConfig == 1422){
    cuts.AddCutMergedCalo("00091113","111116506f032200000","111116506f022700001","0163300000000000"); // EJE

  //pp 13 TeV with triggers
  } else if (trainConfig == 1500){ //no NL, tight timing, E/p TM
    cuts.AddCutMergedCalo("00010113","41179650af032200000","41179650af022700001","0163300000000000"); // INT7
  } else if (trainConfig == 1501){
    cuts.AddCutMergedCalo("0008e113","41179650af032200000","41179650af022700001","0163300000000000"); // EG2+DG2
    cuts.AddCutMergedCalo("0008d113","41179650af032200000","41179650af022700001","0163300000000000"); // EG1+DG1
  } else if (trainConfig == 1502){
    cuts.AddCutMergedCalo("0009c113","41179650af032200000","41179650af022700001","0163300000000000"); // EJ2+DJ2
    cuts.AddCutMergedCalo("0009b113","41179650af032200000","41179650af022700001","0163300000000000"); // EJ1+DJ1
  } else if (trainConfig == 1503){ //no NL, -50ns, 30ns timing cut, E/p TM
    cuts.AddCutMergedCalo("00010113","411796505f032200000","411796505f022700001","0163300000000000"); // INT7
  } else if (trainConfig == 1504){
    cuts.AddCutMergedCalo("0008e113","411796505f032200000","411796505f022700001","0163300000000000"); // EG2+DG2
    cuts.AddCutMergedCalo("0008d113","411796505f032200000","411796505f022700001","0163300000000000"); // EG1+DG1
  } else if (trainConfig == 1505){
    cuts.AddCutMergedCalo("0009c113","411796505f032200000","411796505f022700001","0163300000000000"); // EJ2+DJ2
    cuts.AddCutMergedCalo("0009b113","411796505f032200000","411796505f022700001","0163300000000000"); // EJ1+DJ1
  } else if (trainConfig == 1506){ //no NL, open timing, E/p TM
    cuts.AddCutMergedCalo("00010113","411793100f032200000","411793100f022700001","0163300000000000"); // INT7
  } else if (trainConfig == 1507){
    cuts.AddCutMergedCalo("0008e113","411793100f032200000","411793100f022700001","0163300000000000"); // EG2+DG2
    cuts.AddCutMergedCalo("0008d113","411793100f032200000","411793100f022700001","0163300000000000"); // EG1+DG1
  } else if (trainConfig == 1508){
    cuts.AddCutMergedCalo("0009c113","411793100f032200000","411793100f022700001","0163300000000000"); // EJ2+DJ2
    cuts.AddCutMergedCalo("0009b113","411793100f032200000","411793100f022700001","0163300000000000"); // EJ1+DJ1

  } else if (trainConfig == 1511){ // -50ns, 30ns timing cut, E/p TM, exotics variation
    cuts.AddCutMergedCalo("0008d113","411796505f232200000","411796505f222700001","0163300000000000"); // EG1+DG1
    cuts.AddCutMergedCalo("0008d113","411796505f532200000","411796505f522700001","0163300000000000"); // EG1+DG1
    cuts.AddCutMergedCalo("0008d113","411796505f832200000","411796505f822700001","0163300000000000"); // EG1+DG1
  } else if (trainConfig == 1512){ // -50ns, 30ns timing cut, E/p TM, exotics variation with T-Card
    cuts.AddCutMergedCalo("0008d113","411796505fa32200000","411796505fa22700001","0163300000000000"); // EG1+DG1
    cuts.AddCutMergedCalo("0008d113","411796505fb32200000","411796505fb22700001","0163300000000000"); // EG1+DG1
    cuts.AddCutMergedCalo("0008d113","411796505fc32200000","411796505fc22700001","0163300000000000"); // EG1+DG1
    cuts.AddCutMergedCalo("0008d113","411796505fd32200000","411796505fd22700001","0163300000000000"); // EG1+DG1
    cuts.AddCutMergedCalo("0008d113","411796505fe32200000","411796505fe22700001","0163300000000000"); // EG1+DG1
    cuts.AddCutMergedCalo("0008d113","411796505ff32200000","411796505ff22700001","0163300000000000"); // EG1+DG1

  } else if (trainConfig == 1513){ //TB NL, -30ns, 35ns timing cut, E/p TM, with exotic cut (F+=0.95, TCard requirement > 50GeV)
    cuts.AddCutMergedCalo("00010113","411790106fe32200000","411790106fe22700001","0163300000000000"); // INT7
  } else if (trainConfig == 1514){
    cuts.AddCutMergedCalo("0008e113","411790106fe32200000","411790106fe22700001","0163300000000000"); // EG2+DG2
    cuts.AddCutMergedCalo("0008d113","411790106fe32200000","411790106fe22700001","0163300000000000"); // EG1+DG1
  } else if (trainConfig == 1515){ //TB NL, -30ns, 35ns timing cut, No TM, with exotic cut (F+=0.95, TCard requirement > 50GeV)
    cuts.AddCutMergedCalo("00010113","4117901060e32200000","4117901060e22700001","0163300000000000"); // INT7 No TM
  } else if (trainConfig == 1516){
    cuts.AddCutMergedCalo("0008e113","4117901060e32200000","4117901060e22700001","0163300000000000"); // EG2+DG2 No TM
    cuts.AddCutMergedCalo("0008d113","4117901060e32200000","4117901060e22700001","0163300000000000"); // EG1+DG1 No TM

    // NLM variations for V1 clusterizer studies
  } else if (trainConfig == 1517){ //TB NL, -30ns, 35ns timing cut, E/p TM, with exotic cut (F+=0.95, TCard requirement > 50GeV)
    cuts.AddCutMergedCalo("00010113","411790106fe32200000","411790106fe22700001","0163300000000000"); // INT7 NLM = 1
    cuts.AddCutMergedCalo("00010113","411790106fe32200000","411790106fe22700002","0163300000000000"); // INT7 NLM = 2
    cuts.AddCutMergedCalo("00010113","411790106fe32200000","411790106fe22700004","0163300000000000"); // INT7 NLM > 2
  } else if (trainConfig == 1518){
    cuts.AddCutMergedCalo("0008e113","411790106fe32200000","411790106fe22700001","0163300000000000"); // EG2+DG2 INT7 NLM = 1
    cuts.AddCutMergedCalo("0008e113","411790106fe32200000","411790106fe22700002","0163300000000000"); // EG2+DG2 INT7 NLM = 2
    cuts.AddCutMergedCalo("0008e113","411790106fe32200000","411790106fe22700004","0163300000000000"); // EG2+DG2 NLM > 2
  } else if (trainConfig == 1519){
    cuts.AddCutMergedCalo("0008d113","411790106fe32200000","411790106fe22700001","0163300000000000"); // EG1+DG1 INT7 NLM = 1
    cuts.AddCutMergedCalo("0008d113","411790106fe32200000","411790106fe22700002","0163300000000000"); // EG1+DG1 INT7 NLM = 2
    cuts.AddCutMergedCalo("0008d113","411790106fe32200000","411790106fe22700004","0163300000000000"); // EG1+DG1 NLM > 2
    // NLM variations for V1 clusterizer studies w/o track matching
  } else if (trainConfig == 1520){ //TB NL, -30ns, 35ns timing cut, no TM, with exotic cut (F+=0.95, TCard requirement > 50GeV)
    cuts.AddCutMergedCalo("00010113","4117901060e32200000","4117901060e22700001","0163300000000000"); // INT7 NLM = 1
    cuts.AddCutMergedCalo("00010113","4117901060e32200000","4117901060e22700002","0163300000000000"); // INT7 NLM = 2
    cuts.AddCutMergedCalo("00010113","4117901060e32200000","4117901060e22700003","0163300000000000"); // INT7 NLM <= 2
    cuts.AddCutMergedCalo("00010113","4117901060e32200000","4117901060e22700004","0163300000000000"); // INT7 NLM > 2
  } else if (trainConfig == 1521){
    cuts.AddCutMergedCalo("0008e113","4117901060e32200000","4117901060e22700001","0163300000000000"); // EG2+DG2 INT7 NLM = 1
    cuts.AddCutMergedCalo("0008e113","4117901060e32200000","4117901060e22700002","0163300000000000"); // EG2+DG2 INT7 NLM = 2
    cuts.AddCutMergedCalo("0008e113","4117901060e32200000","4117901060e22700003","0163300000000000"); // EG2+DG2 INT7 NLM <= 2
    cuts.AddCutMergedCalo("0008e113","4117901060e32200000","4117901060e22700004","0163300000000000"); // EG2+DG2 NLM > 2
  } else if (trainConfig == 1522){
    cuts.AddCutMergedCalo("0008d113","4117901060e32200000","4117901060e22700001","0163300000000000"); // EG1+DG1 INT7 NLM = 1
    cuts.AddCutMergedCalo("0008d113","4117901060e32200000","4117901060e22700002","0163300000000000"); // EG1+DG1 INT7 NLM = 2
    cuts.AddCutMergedCalo("0008d113","4117901060e32200000","4117901060e22700003","0163300000000000"); // EG1+DG1 INT7 NLM <= 2
    cuts.AddCutMergedCalo("0008d113","4117901060e32200000","4117901060e22700004","0163300000000000"); // EG1+DG1 NLM > 2


  // No M02 and M20 cut for QA
  } else if (trainConfig == 1523){ //TB NL, -30ns, 35ns timing cut, E/p TM, with exotic cut (F+=0.95, TCard requirement > 50GeV)
    cuts.AddCutMergedCalo("00010113","411790106fe32200000","411790106fe02000001","0163300000000000"); // INT7
  } else if (trainConfig == 1524){
    cuts.AddCutMergedCalo("0008e113","411790106fe32000000","411790106fe02000001","0163300000000000"); // EG2+DG2
    cuts.AddCutMergedCalo("0008d113","411790106fe32000000","411790106fe02000001","0163300000000000"); // EG1+DG1

  // No M02 and M20 cut for QA (for V1 Clusterizer without NLM cut)
  } else if (trainConfig == 1525){ //TB NL, -30ns, 35ns timing cut, E/p TM, with exotic cut (F+=0.95, TCard requirement > 50GeV)
    cuts.AddCutMergedCalo("00010113","411790106fe32200000","411790106fe02000000","0163300000000000"); // INT7
  } else if (trainConfig == 1526){
    cuts.AddCutMergedCalo("0008e113","411790106fe32000000","411790106fe02000000","0163300000000000"); // EG2+DG2
    cuts.AddCutMergedCalo("0008d113","411790106fe32000000","411790106fe02000000","0163300000000000"); // EG1+DG1

  // No M02 and M20 cut for QA (for V1 Clusterizer NLM cut 1 + 2 NLM)
  } else if (trainConfig == 1527){ //TB NL, -30ns, 35ns timing cut, E/p TM, with exotic cut (F+=0.95, TCard requirement > 50GeV)
    cuts.AddCutMergedCalo("00010113","411790106fe32200000","411790106fe02000003","0163300000000000"); // INT7
  } else if (trainConfig == 1528){
    cuts.AddCutMergedCalo("0008e113","411790106fe32000000","411790106fe02000003","0163300000000000"); // EG2+DG2
    cuts.AddCutMergedCalo("0008d113","411790106fe32000000","411790106fe02000003","0163300000000000"); // EG1+DG1


  } else if (trainConfig == 1530){ //TB NL, -30ns, 35ns timing cut, TM without E/p, with exotic cut (F+=0.95, TCard requirement > 50GeV)
    cuts.AddCutMergedCalo("00010113","411790106ce32200000","411790106ce22700001","0163300000000000"); // INT7
  } else if (trainConfig == 1531){
    cuts.AddCutMergedCalo("0008e113","411790106ce32200000","411790106ce22700001","0163300000000000"); // EG2+DG2
    cuts.AddCutMergedCalo("0008d113","411790106ce32200000","411790106ce22700001","0163300000000000"); // EG1+DG1


  } else if (trainConfig == 1532){ //TB NL, -30ns, 35ns timing cut, E/p TM variation, with exotic cut (F+=0.95, TCard requirement > 50GeV)
    cuts.AddCutMergedCalo("00010113","411790106de32200000","411790106de22700001","0163300000000000"); // INT7
    cuts.AddCutMergedCalo("00010113","411790106de32200000","411790106he22700001","0163300000000000"); // INT7
  } else if (trainConfig == 1533){
    cuts.AddCutMergedCalo("0008e113","411790106de32200000","411790106de22700001","0163300000000000"); // EG2+DG2
    cuts.AddCutMergedCalo("0008e113","411790106he32200000","411790106he22700001","0163300000000000"); // EG2+DG2
    cuts.AddCutMergedCalo("0008d113","411790106de32200000","411790106de22700001","0163300000000000"); // EG1+DG1
    cuts.AddCutMergedCalo("0008d113","411790106he32200000","411790106he22700001","0163300000000000"); // EG1+DG1

    // cuts for overlapp study
  } else if (trainConfig == 1534){ //TB NL, -30ns, 35ns timing cut, E/p TM, with exotic cut (F+=0.95, TCard requirement > 50GeV)
    cuts.AddCutMergedCalo("00010113","411790106fe32200000","411790106fe22700001","0163300000000000"); // INT7 0 overlaps
  } else if (trainConfig == 1535){ //TB NL, -30ns, 35ns timing cut, E/p TM, with exotic cut (F+=0.95, TCard requirement > 50GeV)
    cuts.AddCutMergedCalo("00010113","411790106fe32200000","411790106fe22700001","0163300000000000"); // INT7  1 overlaps
  } else if (trainConfig == 1536){ //TB NL, -30ns, 35ns timing cut, E/p TM, with exotic cut (F+=0.95, TCard requirement > 50GeV)
    cuts.AddCutMergedCalo("00010113","411790106fe32200000","411790106fe22700001","0163300000000000"); // INT7  2 overlaps
  } else if (trainConfig == 1537){ //TB NL, -30ns, 35ns timing cut, E/p TM, with exotic cut (F+=0.95, TCard requirement > 50GeV)
    cuts.AddCutMergedCalo("00010113","411790106fe32200000","411790106fe22700001","0163300000000000"); // INT7  > 1 overlaps
  } else if (trainConfig == 1538){ //TB NL, -30ns, 35ns timing cut, E/p TM, with exotic cut (F+=0.95, TCard requirement > 50GeV)
    cuts.AddCutMergedCalo("00010113","411790106fe32200000","411790106fe22700001","0163300000000000"); // INT7  > 2 overlaps
  } else if (trainConfig == 1539){ //TB NL, -30ns, 35ns timing cut, E/p TM, with exotic cut (F+=0.95, TCard requirement > 50GeV)
    cuts.AddCutMergedCalo("00010113","411790106fe32200000","411790106fe22700001","0163300000000000"); // INT7  > 3 overlaps


  // systematics pp 13 TeV INT7
} else if (trainConfig == 1550){ //TB NL variations
  cuts.AddCutMergedCalo("00010113","411793306fe32200000","411793306fe22700001","0163300000000000"); // INT7 NL 33
  cuts.AddCutMergedCalo("00010113","411793406fe32200000","411793406fe22700001","0163300000000000"); // INT7 NL 34
  cuts.AddCutMergedCalo("00010113","411793806fe32200000","411793806fe22700001","0163300000000000"); // INT7 NL 38
  cuts.AddCutMergedCalo("00010113","411793906fe32200000","411793906fe22700001","0163300000000000"); // INT7 NL 39
} else if (trainConfig == 1551){ //M02 variation
  cuts.AddCutMergedCalo("00010113","411790106fe32200000","411790106fe22d00001","0163300000000000"); // INT7 M02 < 0.33
  cuts.AddCutMergedCalo("00010113","411790106fe32200000","411790106fe22600001","0163300000000000"); // INT7 M02 < 0.3
  cuts.AddCutMergedCalo("00010113","411790106fe32200000","411790106fe22800001","0163300000000000"); // INT7 M02 < 0.25
  cuts.AddCutMergedCalo("00010113","411790106fe32200000","411790106fe22k00001","0163300000000000"); // INT7 M02 < 0.2
} else if (trainConfig == 1552){  // dist to bad channel
  cuts.AddCutMergedCalo("00010113","411790106fe32200000","411790116fe22700001","0163300000000000"); // INT7 dist to bad channel 1
  cuts.AddCutMergedCalo("00010113","411790106fe32200000","411790126fe22700001","0163300000000000"); // INT7 dist to bad channel 2
} else if (trainConfig == 1553){  // exotics
  cuts.AddCutMergedCalo("00010113","411790106fb32200000","411790106fb22700001","0163300000000000"); // INT7 F+ < 0.95
  cuts.AddCutMergedCalo("00010113","411790106f532200000","411790106f522700001","0163300000000000"); // INT7 F+ < 0.97 no TCard
  cuts.AddCutMergedCalo("00010113","411790106fg32200000","411790106fg22700001","0163300000000000"); // INT7 exotics in corr framework
} else if (trainConfig == 1554){  // timing
  cuts.AddCutMergedCalo("00010113","411790105fe32200000","411790105fe22700001","0163300000000000"); // INT7 timing: -50 to 50
  cuts.AddCutMergedCalo("00010113","411790109fe32200000","411790109fe22700001","0163300000000000"); // INT7 timing: -20 to 25
  cuts.AddCutMergedCalo("00010113","41179010afe32200000","41179010afe22700001","0163300000000000"); // INT7 timing: -12.5 to 13
} else if (trainConfig == 1555){  // track matching
  cuts.AddCutMergedCalo("00010113","4117901060e32200000","4117901060e22700001","0163300000000000"); // INT7 no TM
  cuts.AddCutMergedCalo("00010113","4117901067e32200000","4117901067e22700001","0163300000000000"); // INT7 no E/p TM
  cuts.AddCutMergedCalo("00010113","411790106he32200000","411790106he22700001","0163300000000000"); // INT7 E/p = 1.25 (strict)
  cuts.AddCutMergedCalo("00010113","411790106de32200000","411790106de22700001","0163300000000000"); // INT7 E/p = 3.0 (loose)
} else if (trainConfig == 1556){  // rapidity
  cuts.AddCutMergedCalo("00010113","411790106fe32200000","411790106fe22700001","0163100000000000"); // INT7 rap = 0.8
  cuts.AddCutMergedCalo("00010113","411790106fe32200000","411790106fe22700001","0163200000000000"); // INT7 rap = 0.7
} else if (trainConfig == 1557){  // EMCal/DCal
  cuts.AddCutMergedCalo("00010113","111110106fe32200000","411790106fe22700001","0163300000000000"); // INT7 EMCal
  cuts.AddCutMergedCalo("00010113","388550106fe32200000","411790106fe22700001","0163300000000000"); // INT7 DCal
} else if (trainConfig == 1558){  // clusterization settings
  cuts.AddCutMergedCalo("00010113","411790106fe32200000","411790106fe22700001","0163300000000000"); // INT7 Eagg = 150MeV
} else if (trainConfig == 1559){  // clusterization settings
  cuts.AddCutMergedCalo("00010113","411790106fe32200000","411790106fe22700001","0163300000000000"); // INT7 Eagg = 75MeV


} else if (trainConfig == 1570){ //TB NL variations
  cuts.AddCutMergedCalo("0008e113","411793306fe32200000","411793306fe22700001","0163300000000000"); // EG2 NL 33
  cuts.AddCutMergedCalo("0008e113","411793406fe32200000","411793406fe22700001","0163300000000000"); // EG2 NL 34
  cuts.AddCutMergedCalo("0008e113","411793806fe32200000","411793806fe22700001","0163300000000000"); // EG2 NL 38
  cuts.AddCutMergedCalo("0008e113","411793906fe32200000","411793906fe22700001","0163300000000000"); // EG2 NL 39
} else if (trainConfig == 1571){ //M02 variation
  cuts.AddCutMergedCalo("0008e113","411790106fe32200000","411790106fe22d00001","0163300000000000"); // EG2 M02 < 0.33
  cuts.AddCutMergedCalo("0008e113","411790106fe32200000","411790106fe22600001","0163300000000000"); // EG2 M02 < 0.3
  cuts.AddCutMergedCalo("0008e113","411790106fe32200000","411790106fe22800001","0163300000000000"); // EG2 M02 < 0.25
  cuts.AddCutMergedCalo("0008e113","411790106fe32200000","411790106fe22k00001","0163300000000000"); // EG2 M02 < 0.2
} else if (trainConfig == 1572){  // dist to bad channel
  cuts.AddCutMergedCalo("0008e113","411790106fe32200000","411790116fe22700001","0163300000000000"); // EG2 dist to bad channel 1
  cuts.AddCutMergedCalo("0008e113","411790106fe32200000","411790126fe22700001","0163300000000000"); // EG2 dist to bad channel 2
} else if (trainConfig == 1573){  // exotics
  cuts.AddCutMergedCalo("0008e113","411790106fb32200000","411790106fb22700001","0163300000000000"); // EG2 F+ < 0.95
  cuts.AddCutMergedCalo("0008e113","411790106f532200000","411790106f522700001","0163300000000000"); // EG2 F+ < 0.97 no TCard
  cuts.AddCutMergedCalo("0008e113","411790106fg32200000","411790106fg22700001","0163300000000000"); // EG2 exotics in corr framework
} else if (trainConfig == 1574){  // timing
  cuts.AddCutMergedCalo("0008e113","411790105fe32200000","411790105fe22700001","0163300000000000"); // EG2 timing: -50 to 50
  cuts.AddCutMergedCalo("0008e113","411790109fe32200000","411790109fe22700001","0163300000000000"); // EG2 timing: -20 to 25
  cuts.AddCutMergedCalo("0008e113","41179010afe32200000","41179010afe22700001","0163300000000000"); // EG2 timing: -12.5 to 13
} else if (trainConfig == 1575){  // track matching
  cuts.AddCutMergedCalo("0008e113","4117901060e32200000","4117901060e22700001","0163300000000000"); // EG2 no TM
  cuts.AddCutMergedCalo("0008e113","4117901067e32200000","4117901067e22700001","0163300000000000"); // EG2 no E/p TM
  cuts.AddCutMergedCalo("0008e113","411790106he32200000","411790106he22700001","0163300000000000"); // EG2 E/p = 1.25 (strict)
  cuts.AddCutMergedCalo("0008e113","411790106de32200000","411790106de22700001","0163300000000000"); // EG2 E/p = 3.0 (loose)
} else if (trainConfig == 1576){  // rapidity
  cuts.AddCutMergedCalo("0008e113","411790106fe32200000","411790106fe22700001","0163100000000000"); // EG2 rap = 0.8
  cuts.AddCutMergedCalo("0008e113","411790106fe32200000","411790106fe22700001","0163200000000000"); // EG2 rap = 0.7
} else if (trainConfig == 1577){  // EMCal/DCal
  cuts.AddCutMergedCalo("00085113","111110106fe32200000","411790106fe22700001","0163300000000000"); // EG2 EMCal
  cuts.AddCutMergedCalo("00089113","388550106fe32200000","411790106fe22700001","0163300000000000"); // EG2 DCal
} else if (trainConfig == 1578){  // clusterization settings
  cuts.AddCutMergedCalo("0008e113","411790106fe32200000","411790106fe22700001","0163300000000000"); // EG2 Eagg = 150MeV
} else if (trainConfig == 1579){  // clusterization settings
  cuts.AddCutMergedCalo("0008e113","411790106fe32200000","411790106fe22700001","0163300000000000"); // EG2 Eagg = 75MeV


} else if (trainConfig == 1590){ //TB NL variations
  cuts.AddCutMergedCalo("0008d113","411793306fe32200000","411793306fe22700001","0163300000000000"); // EG1 NL 33
  cuts.AddCutMergedCalo("0008d113","411793406fe32200000","411793406fe22700001","0163300000000000"); // EG1 NL 34
  cuts.AddCutMergedCalo("0008d113","411793806fe32200000","411793806fe22700001","0163300000000000"); // EG1 NL 38
  cuts.AddCutMergedCalo("0008d113","411793906fe32200000","411793906fe22700001","0163300000000000"); // EG1 NL 39
} else if (trainConfig == 1591){ //M02 variation
  cuts.AddCutMergedCalo("0008d113","411790106fe32200000","411790106fe22d00001","0163300000000000"); // EG1 M02 < 0.33
  cuts.AddCutMergedCalo("0008d113","411790106fe32200000","411790106fe22600001","0163300000000000"); // EG1 M02 < 0.3
  cuts.AddCutMergedCalo("0008d113","411790106fe32200000","411790106fe22800001","0163300000000000"); // EG1 M02 < 0.25
  cuts.AddCutMergedCalo("0008d113","411790106fe32200000","411790106fe22k00001","0163300000000000"); // EG1 M02 < 0.2
} else if (trainConfig == 1592){  // dist to bad channel
  cuts.AddCutMergedCalo("0008d113","411790106fe32200000","411790116fe22700001","0163300000000000"); // EG1 dist to bad channel 1
  cuts.AddCutMergedCalo("0008d113","411790106fe32200000","411790126fe22700001","0163300000000000"); // EG1 dist to bad channel 2
} else if (trainConfig == 1593){  // exotics
  cuts.AddCutMergedCalo("0008d113","411790106fb32200000","411790106fb22700001","0163300000000000"); // EG1 F+ < 0.95
  cuts.AddCutMergedCalo("0008d113","411790106f532200000","411790106f522700001","0163300000000000"); // EG1 F+ < 0.97 no TCard
  cuts.AddCutMergedCalo("0008d113","411790106fg32200000","411790106fg22700001","0163300000000000"); // EG1 exotics in corr framework
} else if (trainConfig == 1594){  // timing
  cuts.AddCutMergedCalo("0008d113","411790105fe32200000","411790105fe22700001","0163300000000000"); // EG1 timing: -50 to 50
  cuts.AddCutMergedCalo("0008d113","411790109fe32200000","411790109fe22700001","0163300000000000"); // EG1 timing: -20 to 25
  cuts.AddCutMergedCalo("0008d113","41179010afe32200000","41179010afe22700001","0163300000000000"); // EG1 timing: -12.5 to 13
} else if (trainConfig == 1595){  // track matching
  cuts.AddCutMergedCalo("0008d113","4117901060e32200000","4117901060e22700001","0163300000000000"); // EG1 no TM
  cuts.AddCutMergedCalo("0008d113","4117901067e32200000","4117901067e22700001","0163300000000000"); // EG1 no E/p TM
  cuts.AddCutMergedCalo("0008d113","411790106he32200000","411790106he22700001","0163300000000000"); // EG1 E/p = 1.25 (strict)
  cuts.AddCutMergedCalo("0008d113","411790106de32200000","411790106de22700001","0163300000000000"); // EG1 E/p = 3.0 (loose)
} else if (trainConfig == 1596){  // rapidity
  cuts.AddCutMergedCalo("0008d113","411790106fe32200000","411790106fe22700001","0163100000000000"); // EG1 rap = 0.8
  cuts.AddCutMergedCalo("0008d113","411790106fe32200000","411790106fe22700001","0163200000000000"); // EG1 rap = 0.7
} else if (trainConfig == 1597){  // EMCal/DCal
  cuts.AddCutMergedCalo("00083113","111110106fe32200000","411790106fe22700001","0163300000000000"); // EG1 EMCal
  cuts.AddCutMergedCalo("0008b113","388550106fe32200000","411790106fe22700001","0163300000000000"); // EG1 DCal
} else if (trainConfig == 1598){  // clusterization settings
  cuts.AddCutMergedCalo("0008d113","411790106fe32200000","411790106fe22700001","0163300000000000"); // EG1 Eagg = 150MeV
} else if (trainConfig == 1599){  // clusterization settings
  cuts.AddCutMergedCalo("0008d113","411790106fe32200000","411790106fe22700001","0163300000000000"); // EG1 Eagg = 75MeV

// In Jet
} else if (trainConfig == 1700){ //TB NL, -30ns, 35ns timing cut, E/p TM, with exotic cut (F+=0.95, TCard requirement > 50GeV)
  cuts.AddCutMergedCalo("00010113","411790106fe32200000","411790106fe22700001","2163300000000000"); // INT7
} else if (trainConfig == 1701){
  cuts.AddCutMergedCalo("0008e113","411790106fe32200000","411790106fe22700001","2163300000000000"); // EG2+DG2
  cuts.AddCutMergedCalo("0008d113","411790106fe32200000","411790106fe22700001","2163300000000000"); // EG1+DG1
  // Out of jet
} else if (trainConfig == 1702){ //TB NL, -30ns, 35ns timing cut, E/p TM, with exotic cut (F+=0.95, TCard requirement > 50GeV)
  cuts.AddCutMergedCalo("00010113","411790106fe32200000","411790106fe22700001","6163300000000000"); // INT7, out of Jet
} else if (trainConfig == 1703){
  cuts.AddCutMergedCalo("0008e113","411790106fe32200000","411790106fe22700001","6163300000000000"); // EG2+DG2, out of Jet
  cuts.AddCutMergedCalo("0008d113","411790106fe32200000","411790106fe22700001","6163300000000000"); // EG1+DG1, out of Jet
} else if (trainConfig == 1704){ //TB NL, -30ns, 35ns timing cut, E/p TM, with exotic cut (F+=0.95, TCard requirement > 50GeV)
  cuts.AddCutMergedCalo("00010113","411790106fe32200000","411790106fe22700001","7163300000000000"); // INT7, out of Jet on away side
} else if (trainConfig == 1705){
  cuts.AddCutMergedCalo("0008e113","411790106fe32200000","411790106fe22700001","7163300000000000"); // EG2+DG2, out of Jet on away side
  cuts.AddCutMergedCalo("0008d113","411790106fe32200000","411790106fe22700001","7163300000000000"); // EG1+DG1, out of Jet on away side
} else if (trainConfig == 1706){ //TB NL, -30ns, 35ns timing cut, E/p TM, with exotic cut (F+=0.95, TCard requirement > 50GeV)
  cuts.AddCutMergedCalo("00010113","411790106fe32200000","411790106fe22700001","8163300000000000"); // INT7, out of Jet with donut shape around jet axis
} else if (trainConfig == 1707){
  cuts.AddCutMergedCalo("0008e113","411790106fe32200000","411790106fe22700001","8163300000000000"); // EG2+DG2, out of Jet with donut shape around jet axis
  cuts.AddCutMergedCalo("0008d113","411790106fe32200000","411790106fe22700001","8163300000000000"); // EG1+DG1, out of Jet with donut shape around jet axis

  // systematics pp 8 TeV no TM
  // MB configs
  } else if (trainConfig == 3900){ // std
    cuts.AddCutMergedCalo("00010113","1111132060032200000","1111132060g22700001","0163300000000000"); // std
  } else if (trainConfig == 3901){ // M02 var 1
    cuts.AddCutMergedCalo("00010113","1111132060g32200000","1111132060g22700001","0163300000000000"); // exotics from CF
    cuts.AddCutMergedCalo("00010113","1111132060e32200000","1111132060e22700001","0163300000000000"); // exotics F097 TCard 50
    cuts.AddCutMergedCalo("00010113","1111132060g32200000","1111132060g22900001","0163300000000000"); // exotics from CF M02 0.10
  } else if (trainConfig == 3902){ // M02 var 2
    cuts.AddCutMergedCalo("00010113","1111132060032200000","1111132060022d00001","0163300000000000"); // min 0.33
    cuts.AddCutMergedCalo("00010113","1111132060032200000","1111132060022e00001","0163300000000000"); // min 0.36
  } else if (trainConfig == 3903){ // M02 var 3
    cuts.AddCutMergedCalo("00010113","1111132060032200000","1111132060022800001","0163300000000000"); // min 0.25
    cuts.AddCutMergedCalo("00010113","1111132060032200000","1111132060022a00001","0163300000000000"); // min 0.26
  } else if (trainConfig == 3904){  // NL variations 1
    cuts.AddCutMergedCalo("00010113","1111133060032200000","1111133060022700001","0163300000000000"); // 32
    cuts.AddCutMergedCalo("00010113","1111134060032200000","1111134060022700001","0163300000000000"); // 33
  } else if (trainConfig == 3905){ // NL variations 2
    cuts.AddCutMergedCalo("00010113","1111138060032200000","1111138060022700001","0163300000000000"); // 38
    cuts.AddCutMergedCalo("00010113","1111139060032200000","1111139060022700001","0163300000000000"); // 39
  } else if (trainConfig == 3906){  // dist to bad channel
    cuts.AddCutMergedCalo("00010113","1111132060032200000","1111132160022700001","0163300000000000"); // 1 <= coll+row
    cuts.AddCutMergedCalo("00010113","1111132060032200000","1111132260022700001","0163300000000000"); // 2 <= coll+row
  } else if (trainConfig == 3907){ // M02 var 1
    cuts.AddCutMergedCalo("00010113","1111132060e32200000","1111132060e22k00001","0163300000000000"); // exotics F097 TCard 50 M02 0.20
    cuts.AddCutMergedCalo("00010113","1111132060d32200000","1111132060d22900001","0163300000000000"); // exotics F097 TCard 40 M02 0.10
    cuts.AddCutMergedCalo("00010113","1111132060a32200000","1111132060a22900001","0163300000000000"); // exotics F095 TCard 40 M02 0.10
  } else if (trainConfig == 3908){  // timing cuts 3
    cuts.AddCutMergedCalo("00010113","1111132040032200000","1111132040022700001","0163300000000000"); // 10ns
    cuts.AddCutMergedCalo("00010113","1111132070032200000","1111132070022700001","0163300000000000"); // 30ns
  } else if (trainConfig == 3909){  // timing cuts
    cuts.AddCutMergedCalo("00010113","1111132090032200000","1111132090022700001","0163300000000000"); // -20 to 25
    cuts.AddCutMergedCalo("00010113","1111132050032200000","1111132050022700001","0163300000000000"); // -50 to 50
  } else if (trainConfig == 3910){  // timing cuts 2
    cuts.AddCutMergedCalo("00010113","1111132030032200000","1111132030022700001","0163300000000000"); // -200 to 200
    cuts.AddCutMergedCalo("00010113","1111132020032200000","1111132020022700001","0163300000000000"); // -500 to 500
  } else if (trainConfig == 3911){
    cuts.AddCutMergedCalo("00010113","1111132060032200000","1111132060022700001","0163200000000000"); // rapidity < 0.7
    cuts.AddCutMergedCalo("00010113","1111132060032200000","1111132060022700001","0163100000000000"); // rapidity < 0.8
  } else if (trainConfig == 3912){ // M02 var 3
    cuts.AddCutMergedCalo("00010113","1111132060032200000","1111132060022g00001","0163300000000000"); // min 0.24
    cuts.AddCutMergedCalo("00010113","1111132060032200000","1111132060022h00001","0163300000000000"); // min 0.23
    cuts.AddCutMergedCalo("00010113","1111132060032200000","1111132060022j00001","0163300000000000"); // min 0.21
  } else if (trainConfig == 3913){ // M02 var 4
    cuts.AddCutMergedCalo("00010113","1111132060032200000","1111132060022k00001","0163300000000000"); // min 0.20
    cuts.AddCutMergedCalo("00010113","1111132060032200000","1111132060022900001","0163300000000000"); // min 0.10
    cuts.AddCutMergedCalo("00010113","1111132060032200000","1111132060022m00001","0163300000000000"); // min 0.18
  } else if (trainConfig == 3914){ // M02 var 4
    cuts.AddCutMergedCalo("00010113","1111132060032200000","1111132060022700001","0163400000000000"); // std
    cuts.AddCutMergedCalo("00010113","1111132060032200000","1111132060022700001","0163c00000000000"); // std
    cuts.AddCutMergedCalo("00010113","1111132060032200000","1111132060022700001","0163d00000000000"); // std

  // EMC7 configs
  } else if (trainConfig == 3920){ // std
    cuts.AddCutMergedCalo("00052113","1111132060032200000","1111132060g22700001","0163300000000000"); // std
  } else if (trainConfig == 3921){ // M02 var 1
    cuts.AddCutMergedCalo("00052113","1111132060g32200000","1111132060g22700001","0163300000000000"); // exotics from CF
    cuts.AddCutMergedCalo("00052113","1111132060e32200000","1111132060e22700001","0163300000000000"); // exotics F097 TCard 50
    cuts.AddCutMergedCalo("00052113","1111132060g32200000","1111132060g22900001","0163300000000000"); // exotics from CF M02 0.10
  } else if (trainConfig == 3922){ // M02 var 2
    cuts.AddCutMergedCalo("00052113","1111132060032200000","1111132060022d00001","0163300000000000"); // min 0.33
    cuts.AddCutMergedCalo("00052113","1111132060032200000","1111132060022e00001","0163300000000000"); // min 0.36
  } else if (trainConfig == 3923){ // M02 var 3
    cuts.AddCutMergedCalo("00052113","1111132060032200000","1111132060022800001","0163300000000000"); // min 0.25
    cuts.AddCutMergedCalo("00052113","1111132060032200000","1111132060022a00001","0163300000000000"); // min 0.26
  } else if (trainConfig == 3924){  // NL variations 1
    cuts.AddCutMergedCalo("00052113","1111133060032200000","1111133060022700001","0163300000000000"); // 32
    cuts.AddCutMergedCalo("00052113","1111134060032200000","1111134060022700001","0163300000000000"); // 33
  } else if (trainConfig == 3925){ // NL variations 2
    cuts.AddCutMergedCalo("00052113","1111138060032200000","1111138060022700001","0163300000000000"); // 38
    cuts.AddCutMergedCalo("00052113","1111139060032200000","1111139060022700001","0163300000000000"); // 39
  } else if (trainConfig == 3926){  // dist to bad channel
    cuts.AddCutMergedCalo("00052113","1111132060032200000","1111132160022700001","0163300000000000"); // 1 <= coll+row
    cuts.AddCutMergedCalo("00052113","1111132060032200000","1111132260022700001","0163300000000000"); // 2 <= coll+row
  } else if (trainConfig == 3927){ // M02 var 1
    cuts.AddCutMergedCalo("00052113","1111132060e32200000","1111132060e22k00001","0163300000000000"); // exotics F097 TCard 50 M02 0.20
    cuts.AddCutMergedCalo("00052113","1111132060d32200000","1111132060d22900001","0163300000000000"); // exotics F097 TCard 40 M02 0.10
    cuts.AddCutMergedCalo("00052113","1111132060a32200000","1111132060a22900001","0163300000000000"); // exotics F095 TCard 40 M02 0.10
  } else if (trainConfig == 3928){  // timing cuts 3
    cuts.AddCutMergedCalo("00052113","1111132040032200000","1111132040022700001","0163300000000000"); // 10ns
    cuts.AddCutMergedCalo("00052113","1111132070032200000","1111132070022700001","0163300000000000"); // 30ns
  } else if (trainConfig == 3929){  // timing cuts
    cuts.AddCutMergedCalo("00052113","1111132090032200000","1111132090022700001","0163300000000000"); // -20 to 25
    cuts.AddCutMergedCalo("00052113","1111132050032200000","1111132050022700001","0163300000000000"); // -50 to 50
  } else if (trainConfig == 3930){  // timing cuts 2
    cuts.AddCutMergedCalo("00052113","1111132030032200000","1111132030022700001","0163300000000000"); // -200 to 200
    cuts.AddCutMergedCalo("00052113","1111132020032200000","1111132020022700001","0163300000000000"); // -500 to 500
  } else if (trainConfig == 3931){
    cuts.AddCutMergedCalo("00052113","1111132060032200000","1111132060022700001","0163200000000000"); // rapidity < 0.7
    cuts.AddCutMergedCalo("00052113","1111132060032200000","1111132060022700001","0163100000000000"); // rapidity < 0.8
  } else if (trainConfig == 3932){ // M02 var 3
    cuts.AddCutMergedCalo("00052113","1111132060032200000","1111132060022g00001","0163300000000000"); // min 0.24
    cuts.AddCutMergedCalo("00052113","1111132060032200000","1111132060022h00001","0163300000000000"); // min 0.23
    cuts.AddCutMergedCalo("00052113","1111132060032200000","1111132060022j00001","0163300000000000"); // min 0.21
  } else if (trainConfig == 3933){ // M02 var 4
    cuts.AddCutMergedCalo("00052113","1111132060032200000","1111132060022k00001","0163300000000000"); // min 0.20
    cuts.AddCutMergedCalo("00052113","1111132060032200000","1111132060022900001","0163300000000000"); // min 0.10
    cuts.AddCutMergedCalo("00052113","1111132060032200000","1111132060022m00001","0163300000000000"); // min 0.18
  } else if (trainConfig == 3934){ // M02 var 4
    cuts.AddCutMergedCalo("00052113","1111132060032200000","1111132060022700001","0163400000000000"); // std
    cuts.AddCutMergedCalo("00052113","1111132060032200000","1111132060022700001","0163c00000000000"); // std
    cuts.AddCutMergedCalo("00052113","1111132060032200000","1111132060022700001","0163d00000000000"); // std
  // EGA configs
  } else if (trainConfig == 3940){ // std
    cuts.AddCutMergedCalo("00081113","1111132060032200000","1111132060g22700001","0163300000000000"); // std
  } else if (trainConfig == 3941){ // M02 var 1
    cuts.AddCutMergedCalo("00081113","1111132060032200000","1111132060g22700001","0163300000000000"); // exotics from CF
    cuts.AddCutMergedCalo("00081113","1111132060032200000","1111132060e22700001","0163300000000000"); // exotics F097 TCard 50
    cuts.AddCutMergedCalo("00081113","1111132060g32200000","1111132060g22900001","0163300000000000"); // exotics from CF M02 0.10
  } else if (trainConfig == 3942){ // M02 var 2
    cuts.AddCutMergedCalo("00081113","1111132060032200000","1111132060022d00001","0163300000000000"); // min 0.33
    cuts.AddCutMergedCalo("00081113","1111132060032200000","1111132060022e00001","0163300000000000"); // min 0.36
  } else if (trainConfig == 3943){ // M02 var 3
    cuts.AddCutMergedCalo("00081113","1111132060032200000","1111132060022800001","0163300000000000"); // min 0.25
    cuts.AddCutMergedCalo("00081113","1111132060032200000","1111132060022a00001","0163300000000000"); // min 0.26
  } else if (trainConfig == 3944){  // NL variations 1
    cuts.AddCutMergedCalo("00081113","1111133060032200000","1111133060022700001","0163300000000000"); // 32
    cuts.AddCutMergedCalo("00081113","1111134060032200000","1111134060022700001","0163300000000000"); // 33
  } else if (trainConfig == 3945){ // NL variations 2
    cuts.AddCutMergedCalo("00081113","1111138060032200000","1111138060022700001","0163300000000000"); // 38
    cuts.AddCutMergedCalo("00081113","1111139060032200000","1111139060022700001","0163300000000000"); // 39
  } else if (trainConfig == 3946){  // dist to bad channel
    cuts.AddCutMergedCalo("00081113","1111132060032200000","1111132160022700001","0163300000000000"); // 1 <= coll+row
    cuts.AddCutMergedCalo("00081113","1111132060032200000","1111132260022700001","0163300000000000"); // 2 <= coll+row
  } else if (trainConfig == 3947){ // M02 var 1
    cuts.AddCutMergedCalo("00081113","1111132060e32200000","1111132060e22k00001","0163300000000000"); // exotics F097 TCard 50 M02 0.20
    cuts.AddCutMergedCalo("00081113","1111132060d32200000","1111132060d22900001","0163300000000000"); // exotics F097 TCard 40 M02 0.10
    cuts.AddCutMergedCalo("00081113","1111132060a32200000","1111132060a22900001","0163300000000000"); // exotics F095 TCard 40 M02 0.10
  } else if (trainConfig == 3948){  // timing cuts 3
    cuts.AddCutMergedCalo("00081113","1111132040032200000","1111132040022700001","0163300000000000"); // 10ns
    cuts.AddCutMergedCalo("00081113","1111132070032200000","1111132070022700001","0163300000000000"); // 30ns
  } else if (trainConfig == 3949){  // timing cuts
    cuts.AddCutMergedCalo("00081113","1111132090032200000","1111132090022700001","0163300000000000"); // -20 to 25
    cuts.AddCutMergedCalo("00081113","1111132050032200000","1111132050022700001","0163300000000000"); // -50 to 50
  } else if (trainConfig == 3950){  // timing cuts 2
    cuts.AddCutMergedCalo("00081113","1111132030032200000","1111132030022700001","0163300000000000"); // -200 to 200
    cuts.AddCutMergedCalo("00081113","1111132020032200000","1111132020022700001","0163300000000000"); // -500 to 500
  } else if (trainConfig == 3951){
    cuts.AddCutMergedCalo("00081113","1111132060032200000","1111132060022700001","0163200000000000"); // rapidity < 0.7
    cuts.AddCutMergedCalo("00081113","1111132060032200000","1111132060022700001","0163100000000000"); // rapidity < 0.8
  } else if (trainConfig == 3952){ // M02 var 3
    cuts.AddCutMergedCalo("00081113","1111132060032200000","1111132060022g00001","0163300000000000"); // min 0.24
    cuts.AddCutMergedCalo("00081113","1111132060032200000","1111132060022h00001","0163300000000000"); // min 0.23
    cuts.AddCutMergedCalo("00081113","1111132060032200000","1111132060022j00001","0163300000000000"); // min 0.21
  } else if (trainConfig == 3953){ // M02 var 4
    cuts.AddCutMergedCalo("00081113","1111132060032200000","1111132060022k00001","0163300000000000"); // min 0.20
    cuts.AddCutMergedCalo("00081113","1111132060032200000","1111132060022900001","0163300000000000"); // min 0.10
    cuts.AddCutMergedCalo("00081113","1111132060032200000","1111132060022m00001","0163300000000000"); // min 0.18
  } else if (trainConfig == 3954){ // M02 var 4
    cuts.AddCutMergedCalo("00081113","1111132060032200000","1111132060022700001","0163400000000000"); // std
    cuts.AddCutMergedCalo("00081113","1111132060032200000","1111132060022700001","0163c00000000000"); // std
    cuts.AddCutMergedCalo("00081113","1111132060032200000","1111132060022700001","0163d00000000000"); // std

  } else if (trainConfig == 4000){ // INT7-based configs
    cuts.AddCutMergedCalo("00010113","1111132060g32200000","1111132060g22700001","0163300000000000"); // std
    cuts.AddCutMergedCalo("00052113","1111132060g32200000","1111132060g22700001","0163300000000000"); // std
    cuts.AddCutMergedCalo("00081113","1111132060g32200000","1111132060g22700001","0163300000000000"); // std
  } else if (trainConfig == 4001){ // INT8-based configs
    cuts.AddCutMergedCalo("00011113","1111132060g32200000","1111132060g22700001","0163300000000000"); // std
    cuts.AddCutMergedCalo("00053113","1111132060g32200000","1111132060g22700001","0163300000000000"); // std
    cuts.AddCutMergedCalo("00082113","1111132060g32200000","1111132060g22700001","0163300000000000"); // std
  } else if (trainConfig == 4010){ // INT7-based configs with TM
    cuts.AddCutMergedCalo("00010113","111113206fg32200000","111113206fg22700001","0163300000000000"); // std
    cuts.AddCutMergedCalo("00052113","111113206fg32200000","111113206fg22700001","0163300000000000"); // std
    cuts.AddCutMergedCalo("00081113","111113206fg32200000","111113206fg22700001","0163300000000000"); // std
  } else if (trainConfig == 4011){ // INT8-based configs with TM
    cuts.AddCutMergedCalo("00011113","111113206fg32200000","111113206fg22700001","0163300000000000"); // std
    cuts.AddCutMergedCalo("00053113","111113206fg32200000","111113206fg22700001","0163300000000000"); // std
    cuts.AddCutMergedCalo("00082113","111113206fg32200000","111113206fg22700001","0163300000000000"); // std
  } else if (trainConfig == 4012){ // INT8-based configs with TM
    cuts.AddCutMergedCalo("00081113","1111132060g32200000","1111132060g22700001","0163300000000000"); // std
  } else if (trainConfig == 4013){ // INT8-based configs with TM
    cuts.AddCutMergedCalo("00081113","1111132060g32200000","1111132060g22700001","0163300000000000"); // std
  } else if (trainConfig == 4014){ // INT8-based configs with TM
    cuts.AddCutMergedCalo("00081113","1111132060g32200000","1111132060g22700001","0163300000000000"); // std
  } else if (trainConfig == 4015){ // INT8-based configs with TM
    cuts.AddCutMergedCalo("00081113","1111132060g32200000","1111132060g22700001","0163300000000000"); // std
  } else if (trainConfig == 4016){ // INT8-based configs with TM
    cuts.AddCutMergedCalo("00081113","1111132060g32200000","1111132060g22700001","0163300000000000"); // std
  } else if (trainConfig == 4017){ // INT8-based configs with TM
    cuts.AddCutMergedCalo("00081113","1111132060g32200000","1111132060g22700001","0163300000000000"); // std
  } else if (trainConfig == 4022){ // INT8-based configs with TM
    cuts.AddCutMergedCalo("00081113","1111132060g32200000","1111132060g22700001","0163300000000000"); // std
  } else if (trainConfig == 4023){ // INT8-based configs with TM
    cuts.AddCutMergedCalo("00081113","1111132060g32200000","1111132060g22700001","0163300000000000"); // std
  } else if (trainConfig == 4024){ // INT8-based configs with TM
    cuts.AddCutMergedCalo("00081113","1111132060g32200000","1111132060g22700001","0163300000000000"); // std
  } else if (trainConfig == 4025){ // INT8-based configs with TM
    cuts.AddCutMergedCalo("00081113","1111132060g32200000","1111132060g22700001","0163300000000000"); // std
  } else if (trainConfig == 4026){ // INT8-based configs with TM
    cuts.AddCutMergedCalo("00081113","1111132060g32200000","1111132060g22700001","0163300000000000"); // std
  } else if (trainConfig == 4027){ // INT8-based configs with TM
    cuts.AddCutMergedCalo("00081113","1111132060g32200000","1111132060g22700001","0163300000000000"); // std

  // systematics pp 5 TeV no TM
  // MB configs
  } else if (trainConfig == 4100){ // M02 var 1
    cuts.AddCutMergedCalo("000a0113","4117901060g32200000","4117901060g22700001","0163300000000000"); // std
  } else if (trainConfig == 4101){ // M02 var 1
    cuts.AddCutMergedCalo("000a0113","4117901060g32200000","4117901060g22600001","0163300000000000"); // min 0.3
    cuts.AddCutMergedCalo("000a0113","4117901060g32200000","4117901060g22900001","0163300000000000"); // min 0.1
  } else if (trainConfig == 4102){ // M02 var 2
    cuts.AddCutMergedCalo("000a0113","4117901060g32200000","4117901060g22d00001","0163300000000000"); // min 0.33
    cuts.AddCutMergedCalo("000a0113","4117901060g32200000","4117901060g22e00001","0163300000000000"); // min 0.36
  } else if (trainConfig == 4103){ // M02 var 3
    cuts.AddCutMergedCalo("000a0113","4117901060g32200000","4117901060g22j00001","0163300000000000"); // min 0.21
    cuts.AddCutMergedCalo("000a0113","4117901060g32200000","4117901060g22k00001","0163300000000000"); // min 0.20
  } else if (trainConfig == 4104){  // NL variations 1
    cuts.AddCutMergedCalo("000a0113","4117902060g32200000","4117902060g22700001","0163300000000000"); // 02 Run2 FineTune
    // cuts.AddCutMergedCalo("000a0113","4117934060g32200000","4117934060g22700001","0163300000000000"); // 33
  } else if (trainConfig == 4105){ // NL variations 2
    cuts.AddCutMergedCalo("000a0113","4117938060g32200000","4117938060g22700001","0163300000000000"); // 38
    cuts.AddCutMergedCalo("000a0113","4117939060g32200000","4117939060g22700001","0163300000000000"); // 39
  } else if (trainConfig == 4106){  // dist to bad channel
    cuts.AddCutMergedCalo("000a0113","4117901060g32200000","4117901160g22700001","0163300000000000"); // 1 <= coll+row
    cuts.AddCutMergedCalo("000a0113","4117901060g32200000","4117901260g22700001","0163300000000000"); // 2 <= coll+row
  } else if (trainConfig == 4107){ // dist to bad channel
    cuts.AddCutMergedCalo("000a0113","4117901060g32200000","4117901560g22700001","0163300000000000"); // 1 <= row || coll || 0.5*(coll+row)
    cuts.AddCutMergedCalo("000a0113","4117901060g32200000","4117901660g22700001","0163300000000000"); // 2 <= row || coll || 0.5*(coll+row)
  } else if (trainConfig == 4108){  // timing cuts 3
    cuts.AddCutMergedCalo("000a0113","4117901040g32200000","4117901040g22700001","0163300000000000"); // 10ns
    cuts.AddCutMergedCalo("000a0113","4117901070g32200000","4117901070g22700001","0163300000000000"); // 30ns
  } else if (trainConfig == 4109){  // timing cuts
    cuts.AddCutMergedCalo("000a0113","4117901090g32200000","4117901090g22700001","0163300000000000"); // -20 to 25
    cuts.AddCutMergedCalo("000a0113","4117901050g32200000","4117901050g22700001","0163300000000000"); // -50 to 50
  } else if (trainConfig == 4110){  // timing cuts 2
    cuts.AddCutMergedCalo("000a0113","4117901030g32200000","4117901030g22700001","0163300000000000"); // -200 to 200
    cuts.AddCutMergedCalo("000a0113","4117901020g32200000","4117901020g22700001","0163300000000000"); // -500 to 500
  } else if (trainConfig == 4111){  // rapidity
    cuts.AddCutMergedCalo("000a0113","4117901060g32200000","4117901060g22700001","0163200000000000"); // rapidity < 0.7
    cuts.AddCutMergedCalo("000a0113","4117901060g32200000","4117901060g22700001","0163100000000000"); // rapidity < 0.8
  } else if (trainConfig == 4112){ // EMCal-only
    cuts.AddCutMergedCalo("000a0113","1111101060g32200000","1111101060g22700001","0163300000000000"); // std
  } else if (trainConfig == 4113){ // mult
    cuts.AddCutMergedCalo("m01a0113","4117901060g32200000","4117901060g22700001","0163300000000000"); // std
    cuts.AddCutMergedCalo("m15a0113","4117901060g32200000","4117901060g22700001","0163300000000000"); // std
    cuts.AddCutMergedCalo("m5ka0113","4117901060g32200000","4117901060g22700001","0163300000000000"); // std
  } else if (trainConfig == 4114){ // mult
    cuts.AddCutMergedCalo("n24a0113","4117901060g32200000","4117901060g22700001","0163300000000000"); // std
    cuts.AddCutMergedCalo("n47a0113","4117901060g32200000","4117901060g22700001","0163300000000000"); // std
    cuts.AddCutMergedCalo("n7aa0113","4117901060g32200000","4117901060g22700001","0163300000000000"); // std


  // EMC7 configs
  } else if (trainConfig == 4120){ // M02 var 1
    cuts.AddCutMergedCalo("000aq113","4117901060g32200000","4117901060g22700001","0163300000000000"); // std
  } else if (trainConfig == 4121){ // M02 var 1
    cuts.AddCutMergedCalo("000aq113","4117901060g32200000","4117901060g22600001","0163300000000000"); // min 0.3
    cuts.AddCutMergedCalo("000aq113","4117901060g32200000","4117901060g22900001","0163300000000000"); // min 0.1
  } else if (trainConfig == 4122){ // M02 var 2
    cuts.AddCutMergedCalo("000aq113","4117901060g32200000","4117901060g22d00001","0163300000000000"); // min 0.33
    cuts.AddCutMergedCalo("000aq113","4117901060g32200000","4117901060g22e00001","0163300000000000"); // min 0.36
  } else if (trainConfig == 4123){ // M02 var 3
    cuts.AddCutMergedCalo("000aq113","4117901060g32200000","4117901060g22j00001","0163300000000000"); // min 0.21
    cuts.AddCutMergedCalo("000aq113","4117901060g32200000","4117901060g22k00001","0163300000000000"); // min 0.20
  } else if (trainConfig == 4124){  // NL variations 1
    cuts.AddCutMergedCalo("000aq113","4117902060g32200000","4117902060g22700001","0163300000000000"); // 02 Run2 FineTune
    // cuts.AddCutMergedCalo("000aq113","4117934060g32200000","4117934060g22700001","0163300000000000"); // 33
  } else if (trainConfig == 4125){ // NL variations 2
    cuts.AddCutMergedCalo("000aq113","4117938060g32200000","4117938060g22700001","0163300000000000"); // 38
    cuts.AddCutMergedCalo("000aq113","4117939060g32200000","4117939060g22700001","0163300000000000"); // 39
  } else if (trainConfig == 4126){  // dist to bad channel
    cuts.AddCutMergedCalo("000aq113","4117901060g32200000","4117901160g22700001","0163300000000000"); // 1 <= coll+row
    cuts.AddCutMergedCalo("000aq113","4117901060g32200000","4117901260g22700001","0163300000000000"); // 2 <= coll+row
  } else if (trainConfig == 4127){ // dist to bad channel
    cuts.AddCutMergedCalo("000aq113","4117901060g32200000","4117901560g22700001","0163300000000000"); // 1 <= row || coll || 0.5*(coll+row)
    cuts.AddCutMergedCalo("000aq113","4117901060g32200000","4117901660g22700001","0163300000000000"); // 2 <= row || coll || 0.5*(coll+row)
  } else if (trainConfig == 4128){  // timing cuts 3
    cuts.AddCutMergedCalo("000aq113","4117901040g32200000","4117901040g22700001","0163300000000000"); // 10ns
    cuts.AddCutMergedCalo("000aq113","4117901070g32200000","4117901070g22700001","0163300000000000"); // 30ns
  } else if (trainConfig == 4129){  // timing cuts
    cuts.AddCutMergedCalo("000aq113","4117901090g32200000","4117901090g22700001","0163300000000000"); // -20 to 25
    cuts.AddCutMergedCalo("000aq113","4117901050g32200000","4117901050g22700001","0163300000000000"); // -50 to 50
  } else if (trainConfig == 4130){  // timing cuts 2
    cuts.AddCutMergedCalo("000aq113","4117901030g32200000","4117901030g22700001","0163300000000000"); // -200 to 200
    cuts.AddCutMergedCalo("000aq113","4117901020g32200000","4117901020g22700001","0163300000000000"); // -500 to 500
  } else if (trainConfig == 4131){  // rapidity
    cuts.AddCutMergedCalo("000aq113","4117901060g32200000","4117901060g22700001","0163200000000000"); // rapidity < 0.7
    cuts.AddCutMergedCalo("000aq113","4117901060g32200000","4117901060g22700001","0163100000000000"); // rapidity < 0.8
  } else if (trainConfig == 4132){ // EMCal-only
    cuts.AddCutMergedCalo("000aq113","1111101060g32200000","1111101060g22700001","0163300000000000"); // std
  } else if (trainConfig == 4133){ // mult
    cuts.AddCutMergedCalo("m01aq113","4117901060g32200000","4117901060g22700001","0163300000000000"); // std
    cuts.AddCutMergedCalo("m15aq113","4117901060g32200000","4117901060g22700001","0163300000000000"); // std
    cuts.AddCutMergedCalo("m5kaq113","4117901060g32200000","4117901060g22700001","0163300000000000"); // std
  } else if (trainConfig == 4134){ // mult
    cuts.AddCutMergedCalo("n24aq113","4117901060g32200000","4117901060g22700001","0163300000000000"); // std
    cuts.AddCutMergedCalo("n47aq113","4117901060g32200000","4117901060g22700001","0163300000000000"); // std
    cuts.AddCutMergedCalo("n7aaq113","4117901060g32200000","4117901060g22700001","0163300000000000"); // std

  // EG1 configs
  } else if (trainConfig == 4140){ // M02 var 1
    cuts.AddCutMergedCalo("000am113","4117901060g32200000","4117901060g22700001","0163300000000000"); // std 0.27
  } else if (trainConfig == 4141){ // M02 var 1
    cuts.AddCutMergedCalo("000am113","4117901060g32200000","4117901060g22600001","0163300000000000"); // min 0.3
    cuts.AddCutMergedCalo("000am113","4117901060g32200000","4117901060g22900001","0163300000000000"); // min 0.1
  } else if (trainConfig == 4142){ // M02 var 2
    cuts.AddCutMergedCalo("000am113","4117901060g32200000","4117901060g22d00001","0163300000000000"); // min 0.33
    cuts.AddCutMergedCalo("000am113","4117901060g32200000","4117901060g22e00001","0163300000000000"); // min 0.36
  } else if (trainConfig == 4143){ // M02 var 3
    cuts.AddCutMergedCalo("000am113","4117901060g32200000","4117901060g22j00001","0163300000000000"); // min 0.21
    cuts.AddCutMergedCalo("000am113","4117901060g32200000","4117901060g22k00001","0163300000000000"); // min 0.20
  } else if (trainConfig == 4144){  // NL variations 1
    cuts.AddCutMergedCalo("000am113","4117902060g32200000","4117902060g22700001","0163300000000000"); // 02 Run2 FineTune
    // cuts.AddCutMergedCalo("000am113","4117934060g32200000","4117934060g22700001","0163300000000000"); // 33
  } else if (trainConfig == 4145){ // NL variations 2
    cuts.AddCutMergedCalo("000am113","4117938060g32200000","4117938060g22700001","0163300000000000"); // 38
    cuts.AddCutMergedCalo("000am113","4117939060g32200000","4117939060g22700001","0163300000000000"); // 39
  } else if (trainConfig == 4146){  // dist to bad channel
    cuts.AddCutMergedCalo("000am113","4117901060g32200000","4117901160g22700001","0163300000000000"); // 1 <= coll+row
    cuts.AddCutMergedCalo("000am113","4117901060g32200000","4117901260g22700001","0163300000000000"); // 2 <= coll+row
  } else if (trainConfig == 4147){ // dist to bad channel
    cuts.AddCutMergedCalo("000am113","4117901060g32200000","4117901560g22700001","0163300000000000"); // 1 <= row || coll || 0.5*(coll+row)
    cuts.AddCutMergedCalo("000am113","4117901060g32200000","4117901660g22700001","0163300000000000"); // 2 <= row || coll || 0.5*(coll+row)
  } else if (trainConfig == 4148){  // timing cuts 3
    cuts.AddCutMergedCalo("000am113","4117901040g32200000","4117901040g22700001","0163300000000000"); // 10ns
    cuts.AddCutMergedCalo("000am113","4117901070g32200000","4117901070g22700001","0163300000000000"); // 30ns
  } else if (trainConfig == 4149){  // timing cuts
    cuts.AddCutMergedCalo("000am113","4117901090g32200000","4117901090g22700001","0163300000000000"); // -20 to 25
    cuts.AddCutMergedCalo("000am113","4117901050g32200000","4117901050g22700001","0163300000000000"); // -50 to 50
  } else if (trainConfig == 4150){  // timing cuts 2
    cuts.AddCutMergedCalo("000am113","4117901030g32200000","4117901030g22700001","0163300000000000"); // -200 to 200
    cuts.AddCutMergedCalo("000am113","4117901020g32200000","4117901020g22700001","0163300000000000"); // -500 to 500
  } else if (trainConfig == 4151){  // rapidity
    cuts.AddCutMergedCalo("000am113","4117901060g32200000","4117901060g22700001","0163200000000000"); // rapidity < 0.7
    cuts.AddCutMergedCalo("000am113","4117901060g32200000","4117901060g22700001","0163100000000000"); // rapidity < 0.8
  } else if (trainConfig == 4152){ // EMCal-only
    cuts.AddCutMergedCalo("000am113","1111101060g32200000","1111101060g22700001","0163300000000000"); // std
  } else if (trainConfig == 4153){ // mult
    cuts.AddCutMergedCalo("m01am113","4117901060g32200000","4117901060g22700001","0163300000000000"); // std
    cuts.AddCutMergedCalo("m15am113","4117901060g32200000","4117901060g22700001","0163300000000000"); // std
    cuts.AddCutMergedCalo("m5kam113","4117901060g32200000","4117901060g22700001","0163300000000000"); // std
  } else if (trainConfig == 4154){ // mult
    cuts.AddCutMergedCalo("n24am113","4117901060g32200000","4117901060g22700001","0163300000000000"); // std
    cuts.AddCutMergedCalo("n47am113","4117901060g32200000","4117901060g22700001","0163300000000000"); // std
    cuts.AddCutMergedCalo("n7aam113","4117901060g32200000","4117901060g22700001","0163300000000000"); // std


  } else {
    Error(Form("GammaCaloMerged_%i",trainConfig), "wrong trainConfig variable no cuts have been specified for the configuration");
    return;
  }

  if(!cuts.AreValid()){
    cout << "\n\n****************************************************" << endl;
    cout << "ERROR: No valid cuts stored in CutHandlerCaloMerged! Returning..." << endl;
    cout << "****************************************************\n\n" << endl;
    return;
  }

  Int_t numberOfCuts          = cuts.GetNCuts();

  TList *EventCutList         = new TList();
  TList *ClusterCutList       = new TList();
  TList *ClusterMergedCutList = new TList();
  TList *MesonCutList         = new TList();

  TList *HeaderList           = new TList();
  if (generatorName.Contains("LHC12i3")){
    TObjString *Header2       = new TObjString("BOX");
    HeaderList->Add(Header2);
  } else if (generatorName.CompareTo("LHC14e2b")==0){
    TObjString *Header2       = new TObjString("pi0_1");
    HeaderList->Add(Header2);
    TObjString *Header3       = new TObjString("eta_2");
    HeaderList->Add(Header3);
  }

  TString energy      = "";
  TString mcName      = "";
  TString mcNameAdd   = "";
  if (generatorName.Contains("WOSDD")){
    mcNameAdd         = "_WOSDD";
  } else if (generatorName.Contains("WSDD")){
    mcNameAdd         = "_WSDD";
  }
  if (generatorName.Contains("LHC12i3")){
    energy            = "2760GeV";
    mcName            = "Pythia8_LHC12i3";
  } else if (generatorName.Contains("LHC12f1a")){
    energy            = "2760GeV";
    mcName            = "Pythia8_LHC12f1a";
  } else if (generatorName.Contains("LHC12f1b")){
    energy            = "2760GeV";
    mcName            = "Phojet_LHC12f1b";
  } else if (generatorName.Contains("LHC14e2a")){
    energy            = "8TeV";
    mcName            = "Pythia8_LHC14e2a";
  } else if (generatorName.Contains("LHC14e2b")){
    energy            = "8TeV";
    mcName            = "Pythia8_LHC14e2b";
  } else if (generatorName.Contains("LHC14e2c")){
    energy            = "8TeV";
    mcName            = "Phojet_LHC14e2c";
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
  AliConvEventCuts **analysisEventCuts          = new AliConvEventCuts*[numberOfCuts];
  ClusterCutList->SetOwner(kTRUE);
  AliCaloPhotonCuts **analysisClusterCuts       = new AliCaloPhotonCuts*[numberOfCuts];
  ClusterMergedCutList->SetOwner(kTRUE);
  AliCaloPhotonCuts **analysisClusterMergedCuts = new AliCaloPhotonCuts*[numberOfCuts];
  MesonCutList->SetOwner(kTRUE);
  AliConversionMesonCuts **analysisMesonCuts    = new AliConversionMesonCuts*[numberOfCuts];

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

    analysisEventCuts[i]    = new AliConvEventCuts();

    // definition of weighting input
    TString fitNamePi0      = Form("Pi0_Fit_Data_%s",energy.Data());
    TString fitNameEta      = Form("Eta_Fit_Data_%s",energy.Data());

    TString mcInputNamePi0  = "";
    TString mcInputNameEta  = "";
    mcInputNamePi0          = Form("Pi0_%s%s_%s", mcName.Data(), mcNameAdd.Data(), energy.Data() );
    mcInputNameEta          = Form("Eta_%s%s_%s", mcName.Data(), mcNameAdd.Data(), energy.Data() );

    if (doWeightingPart) analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNamePtWeights, mcInputNamePi0, mcInputNameEta, "",fitNamePi0,fitNameEta);

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
    if(periodNameV0Reader.CompareTo("") != 0) analysisEventCuts[i]->SetPeriodEnum(periodNameV0Reader);
    analysisEventCuts[i]->SetLightOutput(enableLightOutput);
    analysisEventCuts[i]->InitializeCutsFromCutString((cuts.GetEventCut(i)).Data());
    EventCutList->Add(analysisEventCuts[i]);
    analysisEventCuts[i]->SetFillCutHistograms("",kFALSE);

    analysisClusterCuts[i]        = new AliCaloPhotonCuts(isMC);
    analysisClusterCuts[i]->SetIsPureCaloCut(2);
    analysisClusterCuts[i]->SetHistoToModifyAcceptance(histoAcc);
    analysisClusterCuts[i]->SetV0ReaderName(V0ReaderName);
    analysisClusterCuts[i]->SetCorrectionTaskSetting(corrTaskSetting);
    analysisClusterCuts[i]->SetCaloTrackMatcherName(TrackMatcherName);
    analysisClusterCuts[i]->SetLightOutput(enableLightOutput);
    analysisClusterCuts[i]->SetExoticsMinCellEnergyCut(minEnergyForExoticsCut);
    analysisClusterCuts[i]->SetExoticsQA(enableExoticsQA);
    analysisClusterCuts[i]->InitializeCutsFromCutString((cuts.GetClusterCut(i)).Data());
    ClusterCutList->Add(analysisClusterCuts[i]);
    analysisClusterCuts[i]->SetExtendedMatchAndQA(enableExtMatchAndQA);
    analysisClusterCuts[i]->SetFillCutHistograms("");

    analysisClusterMergedCuts[i]  = new AliCaloPhotonCuts(isMC);
    analysisClusterMergedCuts[i]->SetIsPureCaloCut(1);
    analysisClusterMergedCuts[i]->SetHistoToModifyAcceptance(histoAcc);
    analysisClusterMergedCuts[i]->SetV0ReaderName(V0ReaderName);
    analysisClusterMergedCuts[i]->SetCorrectionTaskSetting(corrTaskSetting);
    analysisClusterMergedCuts[i]->SetCaloTrackMatcherName(TrackMatcherName);
    analysisClusterMergedCuts[i]->SetLightOutput(enableLightOutput);
    analysisClusterMergedCuts[i]->SetExoticsMinCellEnergyCut(minEnergyForExoticsCut);
//     analysisClusterMergedCuts[i]->SetExoticsQA(enableExoticsQA);
    analysisClusterMergedCuts[i]->InitializeCutsFromCutString((cuts.GetClusterMergedCut(i)).Data());
    ClusterMergedCutList->Add(analysisClusterMergedCuts[i]);
    analysisClusterMergedCuts[i]->SetExtendedMatchAndQA(enableExtMatchAndQA);
    analysisClusterMergedCuts[i]->SetFillCutHistograms("");

    analysisMesonCuts[i]          = new AliConversionMesonCuts();
    analysisMesonCuts[i]->SetEnableOpeningAngleCut(kFALSE);
    analysisMesonCuts[i]->SetIsMergedClusterCut(1);
    analysisMesonCuts[i]->SetLightOutput(enableLightOutput);
    analysisMesonCuts[i]->InitializeCutsFromCutString((cuts.GetMesonCut(i)).Data());
    MesonCutList->Add(analysisMesonCuts[i]);
    analysisMesonCuts[i]->SetFillCutHistograms("");
    analysisEventCuts[i]->SetAcceptedHeader(HeaderList);
  }
  if(minAllowedPi0Overlaps>-1){ task->SetMinNeutralPionOverlapsMC(minAllowedPi0Overlaps);}
  if(maxAllowedPi0Overlaps>-1){ task->SetMaxNeutralPionOverlapsMC(maxAllowedPi0Overlaps);}
  task->SetEnableDetailedM02Distribtuon(runDetailedM02);
  task->SetSelectedMesonID(selectedMeson);
  task->SetCorrectionTaskSetting(corrTaskSetting);
  task->SetEventCutList(numberOfCuts,EventCutList);
  task->SetCaloCutList(numberOfCuts,ClusterCutList);
  task->SetCaloMergedCutList(numberOfCuts,ClusterMergedCutList);
  task->SetMesonCutList(numberOfCuts,MesonCutList);
  task->SetDoMesonQA(enableQAMesonTask); //Attention new switch for Pi0 QA
  task->SetDoClusterQA(enableQAClusterTask);//Attention new switch small for Cluster QA
  task->SetEnableSortingOfMCClusLabels(enableSortingMCLabels);
  task->SetOverlapFromCluster(doOverlapsFromCluster);
  if(enableExtMatchAndQA > 1){ task->SetPlotHistsExtQA(kTRUE);}
  if (enableDetailedPrintout) task->SetEnableDetailedPrintout(enableDetailedPrintout);//Attention new switch small for Cluster QA
  //connect containers
  AliAnalysisDataContainer *coutput =
    mgr->CreateContainer(!(corrTaskSetting.CompareTo("")) ? Form("GammaCaloMerged_%i",trainConfig) : Form("GammaCaloMerged_%i_%s",trainConfig,corrTaskSetting.Data()), TList::Class(),
              AliAnalysisManager::kOutputContainer, Form("GammaCaloMerged_%i.root",trainConfig) );

  mgr->AddTask(task);
  mgr->ConnectInput(task, 0, cinput);
  mgr->ConnectOutput(task, 1, coutput);

  return;
}
