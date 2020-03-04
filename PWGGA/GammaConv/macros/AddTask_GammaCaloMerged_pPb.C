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
//pPb together with all supporting classes
//***************************************************************************************

//***************************************************************************************
//main function
//***************************************************************************************
void AddTask_GammaCaloMerged_pPb(
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
  TString   settingMaxFacPtHard           = "3.",       // maximum factor between hardest jet and ptHard generated
  // settings for weights
  // FPTW:fileNamePtWeights, separate with ;
  TString   fileNameExternalInputs        = "",
  Int_t     doWeightingPart               = 0,        // enable Weighting
  TString   generatorName                 = "",       // period name
  // special settings
  Int_t     selectedMeson                 = 1,        // put flag for selected meson
  Bool_t    enableSortingMCLabels         = kTRUE,    // enable sorting for MC cluster labels
  Bool_t    enableDetailedPrintout        = kFALSE,   // enable detailed printout
  Double_t  minEnergyForExoticsCut        = 1.0,      // minimum energy to be used for exotics CutHandler
  Bool_t    enableExoticsQA               = kFALSE,   // switch to run QA for exotic clusters
  Bool_t    runDetailedM02                = kFALSE,   // switch on very detailed M02 distribution
  // subwagon config
  Int_t     maxAllowedPi0Overlaps         = -1,   // set maximum number of Pi0 overlaps in MC
  TString   additionalTrainConfig         = "0"       // additional counter for trainconfig
  ) {

  AliCutHandlerPCM cuts;


  TString addTaskName                 = "AddTask_GammaMerged_pPb";
  TString sAdditionalTrainConfig      = cuts.GetSpecialSettingFromAddConfig(additionalTrainConfig, "", "", addTaskName);
  if (sAdditionalTrainConfig.Atoi() > 0){
    trainConfig = trainConfig + sAdditionalTrainConfig.Atoi();
    cout << "INFO: " << addTaskName.Data() << " running additionalTrainConfig '" << sAdditionalTrainConfig.Atoi() << "', train config: '" << trainConfig << "'" << endl;
  }

  TString fileNamePtWeights           = cuts.GetSpecialFileNameFromString (fileNameExternalInputs, "FPTW:");
  TString fileNameCustomTriggerMimicOADB   = cuts.GetSpecialFileNameFromString (fileNameExternalInputs, "FTRM:");

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
  if(rmaxFacPtHardSetting->GetEntries()<1){cout << "ERROR: AddTask_GammaCaloMerged_pPb during parsing of settingMaxFacPtHard String '" << settingMaxFacPtHard.Data() << "'" << endl; return;}
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
  AliAnalysisManager *mgr           = AliAnalysisManager::GetAnalysisManager();
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
  AliAnalysisTaskGammaCaloMerged *task  = NULL;
  task                                  = new AliAnalysisTaskGammaCaloMerged(Form("GammaCaloMerged_%i",trainConfig));
  task->SetIsHeavyIon(isHeavyIon);
  task->SetIsMC(isMC);
  task->SetV0ReaderName(V0ReaderName);
  task->SetCorrectionTaskSetting(corrTaskSetting);
  task->SetLightOutput(enableLightOutput);
  task->SetTrackMatcherRunningMode(trackMatcherRunningMode);

  // cluster cuts
  // 0 "ClusterType",  1 "EtaMin", 2 "EtaMax", 3 "PhiMin", 4 "PhiMax", 5 "DistanceToBadChannel", 6 "Timing", 7 "TrackMatching", 8 "ExoticCell",
  // 9 "MinEnergy", 10 "MinNCells", 11 "MinM02", 12 "MaxM02", 13 "MinM20", 14 "MaxM20", 15 "MaximumDispersion", 16 "NLM"

  // ************************************* EMCAL cuts ****************************************************
  // LHC13b-d
  if (trainConfig == 1){ // NLM 1 no non linearity
    cuts.AddCutMergedCalo("80010013","1111100050032200000","1111100050022110001","0163300000000000"); //
    cuts.AddCutMergedCalo("80052013","1111100050032200000","1111100050022110001","0163300000000000"); //
    cuts.AddCutMergedCalo("80085013","1111100050032200000","1111100050022110001","0163300000000000"); //
    cuts.AddCutMergedCalo("80083013","1111100050032200000","1111100050022110001","0163300000000000"); //
  } else if (trainConfig == 2){ // NLM 2 no non linearity
    cuts.AddCutMergedCalo("80010013","1111100050032200000","1111100050022110002","0163302200000000"); //
    cuts.AddCutMergedCalo("80052013","1111100050032200000","1111100050022110002","0163302200000000"); //
    cuts.AddCutMergedCalo("80085013","1111100050032200000","1111100050022110002","0163302200000000"); //
    cuts.AddCutMergedCalo("80083013","1111100050032200000","1111100050022110002","0163302200000000"); //
  } else if (trainConfig == 3){ // NLM 1 conv non linearity
    cuts.AddCutMergedCalo("80010013","1111141053032200000","1111141053022110001","0163301100000000"); //
    cuts.AddCutMergedCalo("80052013","1111141053032200000","1111141053022110001","0163301100000000"); //
    cuts.AddCutMergedCalo("80085013","1111141053032200000","1111141053022110001","0163301100000000"); //
    cuts.AddCutMergedCalo("80083013","1111141053032200000","1111141053022110001","0163301100000000"); //
  } else if (trainConfig == 4){ // NLM 2 conv non linearity
    cuts.AddCutMergedCalo("80010013","1111141053032200000","1111141053022110002","0163302200000000"); //
    cuts.AddCutMergedCalo("80052013","1111141053032200000","1111141053022110002","0163302200000000"); //
    cuts.AddCutMergedCalo("80085013","1111141053032200000","1111141053022110002","0163302200000000"); //
    cuts.AddCutMergedCalo("80083013","1111141053032200000","1111141053022110002","0163302200000000"); //
  } else if (trainConfig == 5){ // NLM 1 conv non linearity
    cuts.AddCutMergedCalo("80010013","1111141053032200000","1111141053022210001","0163301100000000"); //
    cuts.AddCutMergedCalo("80052013","1111141053032200000","1111141053022210001","0163301100000000"); //
    cuts.AddCutMergedCalo("80085013","1111141053032200000","1111141053022210001","0163301100000000"); //
    cuts.AddCutMergedCalo("80083013","1111141053032200000","1111141053022210001","0163301100000000"); //
  } else if (trainConfig == 6){ // NLM 2 conv non linearity
    cuts.AddCutMergedCalo("80010013","1111141053032200000","1111141053022210002","0163302200000000"); //
    cuts.AddCutMergedCalo("80052013","1111141053032200000","1111141053022210002","0163302200000000"); //
    cuts.AddCutMergedCalo("80085013","1111141053032200000","1111141053022210002","0163302200000000"); //
    cuts.AddCutMergedCalo("80083013","1111141053032200000","1111141053022210002","0163302200000000"); //

  } else if (trainConfig == 10){ // pp 2.76TeV paper cuts
    cuts.AddCutMergedCalo("00010113","1111001057032200000","1111100057022700001","0163300000000000"); // w/o NL,
    cuts.AddCutMergedCalo("00010113","1111001057032200000","1111100057022000001","0163300000000000"); // w/o NL, w/o M02
    cuts.AddCutMergedCalo("00010113","1111411057032200000","1111141057022700001","0163300000000000"); //
    cuts.AddCutMergedCalo("00010113","1111411057032200000","1111141057022000001","0163300000000000"); // w/o M02
  } else if (trainConfig == 11){ // pp 2.76TeV paper cuts
    cuts.AddCutMergedCalo("00210113","1111411057032200000","1111141057022700001","0163300000000000"); //
    cuts.AddCutMergedCalo("02410113","1111411057032200000","1111141057022700001","0163300000000000"); //
    cuts.AddCutMergedCalo("04610113","1111411057032200000","1111141057022700001","0163300000000000"); //
    cuts.AddCutMergedCalo("06010113","1111411057032200000","1111141057022700001","0163300000000000"); //

  } else if (trainConfig == 30){ // pp 2.76TeV paper cuts
    cuts.AddCutMergedCalo("80052013","1111001057032200000","1111100057022700001","0163300000000000"); // w/o NL,
    cuts.AddCutMergedCalo("80052013","1111001057032200000","1111100057022000001","0163300000000000"); // w/o NL, w/o M02
    cuts.AddCutMergedCalo("80052013","1111411057032200000","1111141057022700001","0163300000000000"); //
    cuts.AddCutMergedCalo("80052013","1111411057032200000","1111141057022000001","0163300000000000"); // w/o M02
  } else if (trainConfig == 31){ // pp 2.76TeV paper cuts
    cuts.AddCutMergedCalo("80252013","1111001057032200000","1111100057022700001","0163300000000000"); // w/o NL,
    cuts.AddCutMergedCalo("82452013","1111001057032200000","1111100057022000001","0163300000000000"); // w/o NL, w/o M02
    cuts.AddCutMergedCalo("84652013","1111411057032200000","1111141057022700001","0163300000000000"); //
    cuts.AddCutMergedCalo("86052013","1111411057032200000","1111141057022000001","0163300000000000"); // w/o M02

  } else if (trainConfig == 50){ // pp 2.76TeV paper cuts
    cuts.AddCutMergedCalo("80085013","1111001057032200000","1111100057022700001","0163300000000000"); // w/o NL,
    cuts.AddCutMergedCalo("80085013","1111001057032200000","1111100057022000001","0163300000000000"); // w/o NL, w/o M02
    cuts.AddCutMergedCalo("80085013","1111411057032200000","1111141057022700001","0163300000000000"); //
    cuts.AddCutMergedCalo("80085013","1111411057032200000","1111141057022000001","0163300000000000"); // w/o M02
  } else if (trainConfig == 51){ // pp 2.76TeV paper cuts
    cuts.AddCutMergedCalo("80285013","1111001057032200000","1111100057022700001","0163300000000000"); // w/o NL,
    cuts.AddCutMergedCalo("82485013","1111001057032200000","1111100057022000001","0163300000000000"); // w/o NL, w/o M02
    cuts.AddCutMergedCalo("84685013","1111411057032200000","1111141057022700001","0163300000000000"); //
    cuts.AddCutMergedCalo("86085013","1111411057032200000","1111141057022000001","0163300000000000"); // w/o M02

  } else if (trainConfig == 70){ // pp 2.76TeV paper cuts
    cuts.AddCutMergedCalo("80083013","1111001057032200000","1111100057022700001","0163300000000000"); // w/o NL,
    cuts.AddCutMergedCalo("80083013","1111001057032200000","1111100057022000001","0163300000000000"); // w/o NL, w/o M02
    cuts.AddCutMergedCalo("80083013","1111411057032200000","1111141057022700001","0163300000000000"); //
    cuts.AddCutMergedCalo("80083013","1111411057032200000","1111141057022000001","0163300000000000"); // w/o M02
  } else if (trainConfig == 71){ // pp 2.76TeV paper cuts
    cuts.AddCutMergedCalo("80283013","1111001057032200000","1111100057022700001","0163300000000000"); // w/o NL,
    cuts.AddCutMergedCalo("82483013","1111001057032200000","1111100057022000001","0163300000000000"); // w/o NL, w/o M02
    cuts.AddCutMergedCalo("84683013","1111411057032200000","1111141057022700001","0163300000000000"); //
    cuts.AddCutMergedCalo("86083013","1111411057032200000","1111141057022000001","0163300000000000"); // w/o M02



  // run 2 data
  } else if (trainConfig == 101){ // pp 2.76TeV paper cuts : open timing, TB nonlin
    cuts.AddCutMergedCalo("80010113","1111101017032200000","1111101017022700001","0163300000000000"); // INT7
    cuts.AddCutMergedCalo("80085113","1111101017032200000","1111101017022700001","0163300000000000"); // EG2
    cuts.AddCutMergedCalo("80083113","1111101017032200000","1111101017022700001","0163300000000000"); // EG1
  } else if (trainConfig == 102){  // pp 2.76TeV paper cuts:  w/o mass, open timing, TB nonlin
    cuts.AddCutMergedCalo("80010113","1111101017032200000","1111101017022000001","0163300000000000");
    cuts.AddCutMergedCalo("80052113","1111101017032200000","1111101017022000001","0163300000000000"); // EMC7
    cuts.AddCutMergedCalo("80085113","1111101017032200000","1111101017022000001","0163300000000000"); // EG2
    cuts.AddCutMergedCalo("80083113","1111101017032200000","1111101017022000001","0163300000000000"); // EG1
  } else if (trainConfig == 103){ // pp 2.76TeV paper cuts : open timing, TB nonlin
    cuts.AddCutMergedCalo("80010113","1111100017032200000","1111100017022700001","0163300000000000"); // INT7
    cuts.AddCutMergedCalo("80052113","1111100017032200000","1111100017022700001","0163300000000000"); // EMC7
    cuts.AddCutMergedCalo("80085113","1111100017032200000","1111100017022700001","0163300000000000"); // EG2
    cuts.AddCutMergedCalo("80083113","1111100017032200000","1111100017022700001","0163300000000000"); // EG1
  } else if (trainConfig == 104){  // pp 2.76TeV paper cuts:  w/o mass, open timing, TB nonlin
    cuts.AddCutMergedCalo("80010113","1111100017032200000","1111100017022000001","0163300000000000");
    cuts.AddCutMergedCalo("80052113","1111100017032200000","1111100017022000001","0163300000000000"); // EMC7
    cuts.AddCutMergedCalo("80085113","1111100017032200000","1111100017022000001","0163300000000000"); // EG2
    cuts.AddCutMergedCalo("80083113","1111100017032200000","1111100017022000001","0163300000000000"); // EG1
  } else if (trainConfig == 105){ // pp 2.76TeV paper cuts : open timing, TB nonlin, no TM
    cuts.AddCutMergedCalo("80010113","1111101010032200000","1111101010022700001","0163300000000000"); // INT7
    cuts.AddCutMergedCalo("80085113","1111101010032200000","1111101010022700001","0163300000000000"); // EG2
    cuts.AddCutMergedCalo("80083113","1111101010032200000","1111101010022700001","0163300000000000"); // EG1
  } else if (trainConfig == 106){ // pp 2.76TeV paper cuts : open timing, TB nonlin, no TM
    cuts.AddCutMergedCalo("80010103","1111101010032200000","1111101010022700001","0163300000000000"); // INT7
    cuts.AddCutMergedCalo("80085103","1111101010032200000","1111101010022700001","0163300000000000"); // EG2
    cuts.AddCutMergedCalo("80083103","1111101010032200000","1111101010022700001","0163300000000000"); // EG1
  } else if (trainConfig == 107){ // pp 2.76TeV paper cuts : open timing, TB nonlin, no TM
    cuts.AddCutMergedCalo("80010123","1111101010032200000","1111101010022700001","0163300000000000"); // INT7
    cuts.AddCutMergedCalo("80085123","1111101010032200000","1111101010022700001","0163300000000000"); // EG2
    cuts.AddCutMergedCalo("80083123","1111101010032200000","1111101010022700001","0163300000000000"); // EG1
  } else if (trainConfig == 108){ // pp 2.76TeV paper cuts : open timing, TB nonlin, no TM
    cuts.AddCutMergedCalo("80010103","1111101010032200000","1111101010022700001","0163300000000000"); // INT7
    cuts.AddCutMergedCalo("80010113","1111101010032200000","1111101010022700001","0163300000000000"); // INT7
    cuts.AddCutMergedCalo("80010123","1111101010032200000","1111101010022700001","0163300000000000"); // INT7
  } else if (trainConfig == 109){ // pp 2.76TeV paper cuts : open timing, TB nonlin, no TM
    cuts.AddCutMergedCalo("80010143","1111101010032200000","1111101010022700001","0163300000000000"); // INT7
    cuts.AddCutMergedCalo("80085143","1111101010032200000","1111101010022700001","0163300000000000"); // EG2
    cuts.AddCutMergedCalo("80083143","1111101010032200000","1111101010022700001","0163300000000000"); // EG1
    // pPb 8 TeV LHC16 periods
  } else if (trainConfig == 200){ // open timing, TB nonlin, no TM
    cuts.AddCutMergedCalo("80010113","1111101010032200000","1111101010022700001","0163300000000000"); // INT7
  } else if (trainConfig == 201){ // open timing, TB nonlin, no TM
    cuts.AddCutMergedCalo("80052113","1111101010032200000","1111101010022700001","0163300000000000"); // EMC7
  } else if (trainConfig == 202){ // open timing, TB nonlin, no TM
    cuts.AddCutMergedCalo("80085113","1111101010032200000","1111101010022700001","0163300000000000"); // EG2
  } else if (trainConfig == 203){ // open timing, TB nonlin, no TM
    cuts.AddCutMergedCalo("80083113","1111101010032200000","1111101010022700001","0163300000000000"); // EG1
  } else if (trainConfig == 204){ // open timing, TB nonlin
    cuts.AddCutMergedCalo("80010113","1111101017032200000","1111101017022700001","0163300000000000"); // INT7
  } else if (trainConfig == 205){ // open timing, TB nonlin
    cuts.AddCutMergedCalo("80052113","1111101017032200000","1111101017022700001","0163300000000000"); // EMC7
  } else if (trainConfig == 206){ // open timing, TB nonlin
    cuts.AddCutMergedCalo("80085113","1111101017032200000","1111101017022700001","0163300000000000"); // EG2
  } else if (trainConfig == 207){ // open timing, TB nonlin
    cuts.AddCutMergedCalo("80083113","1111101017032200000","1111101017022700001","0163300000000000"); // EG1
  } else if (trainConfig == 410){ // NL kSDM
    cuts.AddCutMergedCalo("80010113","1111142010032200000","1111142010022700001","0163300000000000"); // INT7
    cuts.AddCutMergedCalo("80085113","1111142010032200000","1111142010022700001","0163300000000000"); // EG2
    cuts.AddCutMergedCalo("80083113","1111142010032200000","1111142010022700001","0163300000000000"); // EG1
  } else if (trainConfig == 411){ // NL DPow
    cuts.AddCutMergedCalo("80010113","1111152010032200000","1111152010022700001","0163300000000000"); // INT7
    cuts.AddCutMergedCalo("80085113","1111152010032200000","1111152010022700001","0163300000000000"); // EG2
    cuts.AddCutMergedCalo("80083113","1111152010032200000","1111152010022700001","0163300000000000"); // EG1
  } else if (trainConfig == 413){ // NL TB
    cuts.AddCutMergedCalo("80010113","1111106010032200000","1111106010022700001","0163300000000000"); // INT7
    cuts.AddCutMergedCalo("80085113","1111106010032200000","1111106010022700001","0163300000000000"); // EG2
    cuts.AddCutMergedCalo("80083113","1111106010032200000","1111106010022700001","0163300000000000"); // EG1

  } else if (trainConfig == 421){ // varied exoctics
    cuts.AddCutMergedCalo("80010113","1111142067032200000","1111142067022700001","0163300000000000"); // no exotics
    cuts.AddCutMergedCalo("80010113","1111142067232200000","1111142067222700001","0163300000000000"); // frac = 0.99
    cuts.AddCutMergedCalo("80010113","1111142067332200000","1111142067322700001","0163300000000000"); // frac = 0.98
    cuts.AddCutMergedCalo("80010113","1111142067532200000","1111142067522700001","0163300000000000"); // frac = 0.97
    cuts.AddCutMergedCalo("80010113","1111142067732200000","1111142067722700001","0163300000000000"); // frac = 0.96
    cuts.AddCutMergedCalo("80010113","1111142067932200000","1111142067922700001","0163300000000000"); // frac = 0.94
  } else if (trainConfig == 422){ // varied M02 part 1
    cuts.AddCutMergedCalo("80010113","1111142067032200000","1111142067022600001","0163300000000000"); // min M02 = 0.3
    cuts.AddCutMergedCalo("80010113","1111142067032200000","1111142067022c00001","0163300000000000"); // min M02 = 0.29
    cuts.AddCutMergedCalo("80010113","1111142067032200000","1111142067022b00001","0163300000000000"); // min M02 = 0.28
  } else if (trainConfig == 423){ // varied M02 part 2
    cuts.AddCutMergedCalo("80010113","1111142067032200000","1111142067022700001","0163300000000000"); // min M02 = 0.27
    cuts.AddCutMergedCalo("80010113","1111142067032200000","1111142067022a00001","0163300000000000"); // min M02 = 0.26
    cuts.AddCutMergedCalo("80010113","1111142067032200000","1111142067022800001","0163300000000000"); // min M02 = 0.25
  } else if (trainConfig == 424){  // variation track matching to cluster part 1
    cuts.AddCutMergedCalo("80010113","1111142060032200000","1111142060022700001","0163300000000000"); // no TM
    cuts.AddCutMergedCalo("80010113","1111142063032200000","1111142063022700001","0163300000000000"); // old matching std
  } else if (trainConfig == 425){  // variation track matching to cluster part 2
    cuts.AddCutMergedCalo("80010113","111114206b032200000","111114206b022700001","0163300000000000"); // looser TM
    cuts.AddCutMergedCalo("80010113","111114206a032200000","111114206a022700001","0163300000000000"); // tighter TM
    cuts.AddCutMergedCalo("80010113","1111142068032200000","1111142068022700001","0163300000000000"); // even looser TM
    cuts.AddCutMergedCalo("80010113","1111142066032200000","1111142066022700001","0163300000000000"); // even tighter TM
  } else if (trainConfig == 427){  // std pT dep, |eta| < 0.3, y < 0.3
    cuts.AddCutMergedCalo("80010113","1661142067032200000","1661142067022700001","0163700000000000"); // diff eta/rap cuts
  } else if (trainConfig == 428){ // NL var 1
    cuts.AddCutMergedCalo("80010113","1111100067032200000","1111100067022700001","0163300000000000"); // NL off
    cuts.AddCutMergedCalo("80010113","1111102067032200000","1111102067022700001","0163300000000000"); // Testbeam (v3) for data, Pi0MCv3 for MC
    cuts.AddCutMergedCalo("80010113","1111106067032200000","1111106067022700001","0163300000000000"); // Testbeam (v4) for data, Pi0MCv3 for MC
    cuts.AddCutMergedCalo("80010113","1111142067032200000","1111142067022700001","0163300000000000"); // default (CCRF)
  } else if (trainConfig == 429){ // NL var 2
    cuts.AddCutMergedCalo("80010113","1111141067032200000","1111141067022700001","0163300000000000"); // CRF
    cuts.AddCutMergedCalo("80010113","1111151067032200000","1111151067022700001","0163300000000000"); // CCMF
    cuts.AddCutMergedCalo("80010113","1111152067032200000","1111152067022700001","0163300000000000"); // CMF

  } else if (trainConfig == 480){ // no TM - CALO+CALOFAST triggers
    cuts.AddCutMergedCalo("800a0113","1111101010032200000","1111101010022700001","0163300000000000"); // INT7
    cuts.AddCutMergedCalo("800a1113","1111101010032200000","1111101010022700001","0163300000000000"); // EMC7
    cuts.AddCutMergedCalo("800a2113","1111101010032200000","1111101010022700001","0163300000000000"); // EG2
    cuts.AddCutMergedCalo("800a3113","1111101010032200000","1111101010022700001","0163300000000000"); // EG1
  } else if (trainConfig == 481){ // no TM
    cuts.AddCutMergedCalo("80010113","1111101010032200000","1111101010022700001","0163300000000000"); // INT7
    cuts.AddCutMergedCalo("80052113","1111101010032200000","1111101010022700001","0163300000000000"); // EMC7
    cuts.AddCutMergedCalo("80085113","1111101010032200000","1111101010022700001","0163300000000000"); // EG2
    cuts.AddCutMergedCalo("80083113","1111101010032200000","1111101010022700001","0163300000000000"); // EG1

  } else if (trainConfig == 490){ // no TM - CALO+CALOFAST triggers + testbeam
    cuts.AddCutMergedCalo("800a0113","1111106060032200000","1111106060022700001","0163300000000000"); // INT7
    cuts.AddCutMergedCalo("800a2113","1111106060032200000","1111106060022700001","0163300000000000"); // EG2
    cuts.AddCutMergedCalo("800a3113","1111106060032200000","1111106060022700001","0163300000000000"); // EG1
  } else if (trainConfig == 491){ // no TM - CALO+CALOFAST triggers + testbeam
    cuts.AddCutMergedCalo("800a4113","1111106060032200000","1111106060022700001","0163300000000000"); // EJ2
    cuts.AddCutMergedCalo("800a5113","1111106060032200000","1111106060022700001","0163300000000000"); // EJ1

  // TM OFF + NL kSDM
  } else if (trainConfig == 800){ // MB header only
    cuts.AddCutMergedCalo("80010113","1111142060032200000","1111142060022700001","0163300000000000"); //
  } else if (trainConfig == 801){ // Special header (use with doWeightingPart=4 for both headers)
    cuts.AddCutMergedCalo("80010123","1111142060032200000","1111142060022700001","0163300000000000"); //
  } else if (trainConfig == 802){ // Special header (use with doWeightingPart=0 for JJ header only)
    cuts.AddCutMergedCalo("80010123","1111142060032200000","1111142060022700001","0163300000000000"); //
  } else if (trainConfig == 803){ // Special header (use with doWeightingPart=4 for both headers)
    cuts.AddCutMergedCalo("80010123","1111142060032200000","1111142060022700002","0163300000000000"); // NLM 2 for V1
  } else if (trainConfig == 804){ // Special header (use with doWeightingPart=0 for JJ header only)
    cuts.AddCutMergedCalo("80010123","1111142060032200000","1111142060022700002","0163300000000000"); // NLM 2 for V1

  // TM ON + NL kSDM
  // standard TM (7)
  } else if (trainConfig == 900){ // MB header only
    cuts.AddCutMergedCalo("80010113","1111142067032200000","1111142067022700001","0163300000000000"); //
  } else if (trainConfig == 901){ // Special header (use with doWeightingPart=4 for both headers)
    cuts.AddCutMergedCalo("80010123","1111142067032200000","1111142067022700001","0163300000000000"); //
  } else if (trainConfig == 902){ // Special header (use with doWeightingPart=0 for JJ header only)
    cuts.AddCutMergedCalo("80010123","1111142067032200000","1111142067022700001","0163300000000000"); //
  } else if (trainConfig == 903){ // Special header (use with doWeightingPart=4 for both headers)
    cuts.AddCutMergedCalo("80010123","1111142067032200000","1111142067022700002","0163300000000000"); // NLM 2 for V1
  } else if (trainConfig == 904){ // Special header (use with doWeightingPart=0 for JJ header only)
    cuts.AddCutMergedCalo("80010123","1111142067032200000","1111142067022700002","0163300000000000"); // NLM 2 for V1

  // TM E/p vars
  } else if (trainConfig == 910){ // Special header (use with doWeightingPart=4 for both headers)
    cuts.AddCutMergedCalo("80010123","111114206f032200000","111114206f022700001","0163300000000000"); // fEOverPMax = 1.75
    cuts.AddCutMergedCalo("80010123","111114206g032200000","111114206g022700001","0163300000000000"); // fEOverPMax = 1.5
    cuts.AddCutMergedCalo("80010123","111114206h032200000","111114206h022700001","0163300000000000"); // fEOverPMax = 1.25

  // M02 vars
  } else if (trainConfig == 920){ // Special header (use with doWeightingPart=4 for both headers)
    cuts.AddCutMergedCalo("80010123","1111142067032200000","1111142067022600001","0163300000000000"); // min M02 = 0.3
    cuts.AddCutMergedCalo("80010123","1111142067032200000","1111142067022c00001","0163300000000000"); // min M02 = 0.3
    cuts.AddCutMergedCalo("80010123","1111142067032200000","1111142067022b00001","0163300000000000"); // min M02 = 0.28

  // EG1 triggered configs (same as 800 and 900 but for triggers)
  // TM OFF + NL kSDM
  } else if (trainConfig == 1000){ // MB header only
    cuts.AddCutMergedCalo("80083113","1111142060032200000","1111142060022700001","0163300000000000"); //
  } else if (trainConfig == 1001){ // Special header (use with doWeightingPart=4 for both headers)
    cuts.AddCutMergedCalo("80083123","1111142060032200000","1111142060022700001","0163300000000000"); //
  } else if (trainConfig == 1002){ // Special header (use with doWeightingPart=0 for JJ header only)
    cuts.AddCutMergedCalo("80083123","1111142060032200000","1111142060022700001","0163300000000000"); //
  } else if (trainConfig == 1003){ // Special header (use with doWeightingPart=4 for both headers)
    cuts.AddCutMergedCalo("80083123","1111142060032200000","1111142060022700002","0163300000000000"); // NLM 2 for V1
  } else if (trainConfig == 1004){ // Special header (use with doWeightingPart=0 for JJ header only)
    cuts.AddCutMergedCalo("80083123","1111142060032200000","1111142060022700002","0163300000000000"); // NLM 2 for V1

  // TM ON + NL kSDM
  // standard TM (7)
  } else if (trainConfig == 1100){ // MB header only
    cuts.AddCutMergedCalo("80083113","1111142067032200000","1111142067022700001","0163300000000000"); //
  } else if (trainConfig == 1101){ // Special header (use with doWeightingPart=4 for both headers)
    cuts.AddCutMergedCalo("80083123","1111142067032200000","1111142067022700001","0163300000000000"); //
  } else if (trainConfig == 1102){ // Special header (use with doWeightingPart=0 for JJ header only)
    cuts.AddCutMergedCalo("80083123","1111142067032200000","1111142067022700001","0163300000000000"); //
  } else if (trainConfig == 1103){ // Special header (use with doWeightingPart=4 for both headers)
    cuts.AddCutMergedCalo("80083123","1111142067032200000","1111142067022700002","0163300000000000"); // NLM 2 for V1
  } else if (trainConfig == 1104){ // Special header (use with doWeightingPart=0 for JJ header only)
    cuts.AddCutMergedCalo("80083123","1111142067032200000","1111142067022700002","0163300000000000"); // NLM 2 for V1

  // TM E/p vars
  } else if (trainConfig == 1110){ // Special header (use with doWeightingPart=4 for both headers)
    cuts.AddCutMergedCalo("80083123","111116505f032200000","111116505f022700001","0163300000000000"); // fEOverPMax = 1.75
    cuts.AddCutMergedCalo("80083123","111116505g032200000","111116505g022700001","0163300000000000"); // fEOverPMax = 1.5
    cuts.AddCutMergedCalo("80083123","111116505h032200000","111116505h022700001","0163300000000000"); // fEOverPMax = 1.25
  } else if (trainConfig == 1111){ // Special header (use with doWeightingPart=4 for both headers)
    cuts.AddCutMergedCalo("80085123","111116505f032200000","111116505f022700001","0163300000000000"); // fEOverPMax = 1.75
    cuts.AddCutMergedCalo("80085123","111116505g032200000","111116505g022700001","0163300000000000"); // fEOverPMax = 1.5
    cuts.AddCutMergedCalo("80085123","111116505h032200000","111116505h022700001","0163300000000000"); // fEOverPMax = 1.25

  // 31 NL with TM
  } else if (trainConfig == 1410){
    cuts.AddCutMergedCalo("80010103","111113105f032200000","111113105f022700001","0163300000000000"); // INT7
  } else if (trainConfig == 1411){
    cuts.AddCutMergedCalo("80052103","111113105f032200000","111113105f022700001","0163300000000000"); // INT7
  } else if (trainConfig == 1412){
    cuts.AddCutMergedCalo("80083103","111113105f032200000","111113105f022700001","0163300000000000"); // EG1
    cuts.AddCutMergedCalo("80085103","111113105f032200000","111113105f022700001","0163300000000000"); // EG2
  } else if (trainConfig == 1413){
    cuts.AddCutMergedCalo("80010103","111113205f032200000","111113205f022700001","0163300000000000"); // INT7
  } else if (trainConfig == 1414){
    cuts.AddCutMergedCalo("80083103","111113205f032200000","111113205f022700001","0163300000000000"); // EG1
    cuts.AddCutMergedCalo("80085103","111113205f032200000","111113205f022700001","0163300000000000"); // EG2
  // 31 NL w/o TM
  } else if (trainConfig == 1420){
    cuts.AddCutMergedCalo("80010103","1111131050032200000","1111131050022700001","0163300000000000"); // INT7
  } else if (trainConfig == 1421){
    cuts.AddCutMergedCalo("80052103","1111131050032200000","1111131050022700001","0163300000000000"); // INT7
  } else if (trainConfig == 1422){
    cuts.AddCutMergedCalo("80083103","1111131050032200000","1111131050022700001","0163300000000000"); // EG1
    cuts.AddCutMergedCalo("80085103","1111131050032200000","1111131050022700001","0163300000000000"); // EG2
  } else if (trainConfig == 1423){
    cuts.AddCutMergedCalo("80010103","1111132050032200000","1111132050022700001","0163300000000000"); // INT7
  } else if (trainConfig == 1424){
    cuts.AddCutMergedCalo("80083103","1111132050032200000","1111132050022700001","0163300000000000"); // EG1
    cuts.AddCutMergedCalo("80085103","1111132050032200000","1111132050022700001","0163300000000000"); // EG2

  // 31 NL with TM for special MC header
  } else if (trainConfig == 1430){
    cuts.AddCutMergedCalo("80010123","111113105f032200000","111113105f022700001","0163300000000000"); // INT7
  } else if (trainConfig == 1431){
    cuts.AddCutMergedCalo("80052123","111113105f032200000","111113105f022700001","0163300000000000"); // INT7
  } else if (trainConfig == 1432){
    cuts.AddCutMergedCalo("80083123","111113105f032200000","111113105f022700001","0163300000000000"); // EG1
    cuts.AddCutMergedCalo("80085123","111113105f032200000","111113105f022700001","0163300000000000"); // EG2
  // 31 NL w/o TM for special MC header
  } else if (trainConfig == 1440){
    cuts.AddCutMergedCalo("80010123","1111131050032200000","1111131050022700001","0163300000000000"); // INT7
  } else if (trainConfig == 1441){
    cuts.AddCutMergedCalo("80052123","1111131050032200000","1111131050022700001","0163300000000000"); // INT7
  } else if (trainConfig == 1442){
    cuts.AddCutMergedCalo("80083123","1111131050032200000","1111131050022700001","0163300000000000"); // EG1
    cuts.AddCutMergedCalo("80085123","1111131050032200000","1111131050022700001","0163300000000000"); // EG2


  //EMCAL+DCAL
  // std cut with TB (Nico) NL
  } else if (trainConfig == 3400){
    cuts.AddCutMergedCalo("80010123","411796505f032200000","411796505f022700001","0163300000000000"); // INT7
  } else if (trainConfig == 3401){
    cuts.AddCutMergedCalo("8008e123","411796505f032200000","411796505f022700001","0163300000000000"); // EG2+DG2
    cuts.AddCutMergedCalo("8008d123","411796505f032200000","411796505f022700001","0163300000000000"); // EG1+DG1
  } else if (trainConfig == 3402){
    cuts.AddCutMergedCalo("8009c123","411796505f032200000","411796505f022700001","0163300000000000"); // EJ2+DJ2
    cuts.AddCutMergedCalo("8009b123","411796505f032200000","411796505f022700001","0163300000000000"); // EJ1+DJ1
  } else if (trainConfig == 3403){
    cuts.AddCutMergedCalo("80010123","411796505f032200000","411796505f022700001","0163300000000000"); // INT7
  } else if (trainConfig == 3404){
    cuts.AddCutMergedCalo("8008e123","411796505f032200000","411796505f022700001","0163300000000000"); // EG2+DG2
    cuts.AddCutMergedCalo("8008d123","411796505f032200000","411796505f022700001","0163300000000000"); // EG1+DG1
  } else if (trainConfig == 3405){
    cuts.AddCutMergedCalo("8009c123","411796505f032200000","411796505f022700001","0163300000000000"); // EJ2+DJ2
    cuts.AddCutMergedCalo("8009b123","411796505f032200000","411796505f022700001","0163300000000000"); // EJ1+DJ1

  // std cut with TB+FineTuning NL
  } else if (trainConfig == 3410){
    cuts.AddCutMergedCalo("80010123","411794705f032200000","411794705f022700001","0163300000000000"); // INT7
  } else if (trainConfig == 3411){
    cuts.AddCutMergedCalo("8008e123","411794705f032200000","411794705f022700001","0163300000000000"); // EG2+DG2
    cuts.AddCutMergedCalo("8008d123","411794705f032200000","411794705f022700001","0163300000000000"); // EG1+DG1
  } else if (trainConfig == 3412){
    cuts.AddCutMergedCalo("8009c123","411794705f032200000","411794705f022700001","0163300000000000"); // EJ2+DJ2
    cuts.AddCutMergedCalo("8009b123","411794705f032200000","411794705f022700001","0163300000000000"); // EJ1+DJ1
  } else if (trainConfig == 3413){
    cuts.AddCutMergedCalo("80010123","411794705f032200000","411794705f022700001","0163300000000000"); // INT7
  } else if (trainConfig == 3414){
    cuts.AddCutMergedCalo("8008e123","411794705f032200000","411794705f022700001","0163300000000000"); // EG2+DG2
    cuts.AddCutMergedCalo("8008d123","411794705f032200000","411794705f022700001","0163300000000000"); // EG1+DG1
  } else if (trainConfig == 3415){
    cuts.AddCutMergedCalo("8009c123","411794705f032200000","411794705f022700001","0163300000000000"); // EJ2+DJ2
    cuts.AddCutMergedCalo("8009b123","411794705f032200000","411794705f022700001","0163300000000000"); // EJ1+DJ1

  } else if (trainConfig == 3420){
    cuts.AddCutMergedCalo("80010123","411794805f032200000","411794805f022700001","0163300000000000"); // INT7
  } else if (trainConfig == 3421){
    cuts.AddCutMergedCalo("8008e123","411794805f032200000","411794805f022700001","0163300000000000"); // EG2+DG2
    cuts.AddCutMergedCalo("8008d123","411794805f032200000","411794805f022700001","0163300000000000"); // EG1+DG1
  } else if (trainConfig == 3422){
    cuts.AddCutMergedCalo("8009c123","411794805f032200000","411794805f022700001","0163300000000000"); // EJ2+DJ2
    cuts.AddCutMergedCalo("8009b123","411794805f032200000","411794805f022700001","0163300000000000"); // EJ1+DJ1
  } else if (trainConfig == 3423){ // standard cut
    cuts.AddCutMergedCalo("80010123","411794805f032200000","411794805f022700001","0163300000000000"); // INT7
  } else if (trainConfig == 3424){ // standard cut
    cuts.AddCutMergedCalo("8008e123","411794805f032200000","411794805f022700001","0163300000000000"); // EG2+DG2
    cuts.AddCutMergedCalo("8008d123","411794805f032200000","411794805f022700001","0163300000000000"); // EG1+DG1
  } else if (trainConfig == 3425){ // standard cut
    cuts.AddCutMergedCalo("8009c123","411794805f032200000","411794805f022700001","0163300000000000"); // EJ2+DJ2
    cuts.AddCutMergedCalo("8009b123","411794805f032200000","411794805f022700001","0163300000000000"); // EJ1+DJ1
  } else if (trainConfig == 3426){ // open M02 for QA
    cuts.AddCutMergedCalo("80010123","411794805f032000000","411794805f022000001","0163300000000000"); // INT7
  } else if (trainConfig == 3427){ // open M02 for QA
    cuts.AddCutMergedCalo("8008e123","411794805f032000000","411794805f022000001","0163300000000000"); // EG2+DG2
    cuts.AddCutMergedCalo("8008d123","411794805f032000000","411794805f022000001","0163300000000000"); // EG1+DG1
  } else if (trainConfig == 3428){ // open M02 for QA
    cuts.AddCutMergedCalo("8009c123","411794805f032000000","411794805f022000001","0163300000000000"); // EJ2+DJ2
    cuts.AddCutMergedCalo("8009b123","411794805f032000000","411794805f022000001","0163300000000000"); // EJ1+DJ1

  // configs for pass1 without track matching. three configs always have same NL
  } else if (trainConfig == 3430){
    cuts.AddCutMergedCalo("80010123","4117947050032200000","4117947050022700001","0163300000000000"); // INT7
  } else if (trainConfig == 3431){
    cuts.AddCutMergedCalo("8008e123","4117947050032200000","4117947050022700001","0163300000000000"); // EG2+DG2
    cuts.AddCutMergedCalo("8008d123","4117947050032200000","4117947050022700001","0163300000000000"); // EG1+DG1
  } else if (trainConfig == 3432){
    cuts.AddCutMergedCalo("8009c123","4117947050032200000","4117947050022700001","0163300000000000"); // EJ2+DJ2
    cuts.AddCutMergedCalo("8009b123","4117947050032200000","4117947050022700001","0163300000000000"); // EJ1+DJ1
  } else if (trainConfig == 3433){
    cuts.AddCutMergedCalo("80010123","4117948050032200000","4117948050022700001","0163300000000000"); // INT7
  } else if (trainConfig == 3434){
    cuts.AddCutMergedCalo("8008e123","4117948050032200000","4117948050022700001","0163300000000000"); // EG2+DG2
    cuts.AddCutMergedCalo("8008d123","4117948050032200000","4117948050022700001","0163300000000000"); // EG1+DG1
  } else if (trainConfig == 3435){
    cuts.AddCutMergedCalo("8009c123","4117948050032200000","4117948050022700001","0163300000000000"); // EJ2+DJ2
    cuts.AddCutMergedCalo("8009b123","4117948050032200000","4117948050022700001","0163300000000000"); // EJ1+DJ1

  } else if (trainConfig == 3440){
    cuts.AddCutMergedCalo("80010123","4117957050032200000","4117957050022700001","0163300000000000"); // INT7
  } else if (trainConfig == 3441){
    cuts.AddCutMergedCalo("8008e123","4117957050032200000","4117957050022700001","0163300000000000"); // EG2+DG2
    cuts.AddCutMergedCalo("8008d123","4117957050032200000","4117957050022700001","0163300000000000"); // EG1+DG1
  } else if (trainConfig == 3442){
    cuts.AddCutMergedCalo("8009c123","4117957050032200000","4117957050022700001","0163300000000000"); // EJ2+DJ2
    cuts.AddCutMergedCalo("8009b123","4117957050032200000","4117957050022700001","0163300000000000"); // EJ1+DJ1
  } else if (trainConfig == 3443){
    cuts.AddCutMergedCalo("80010123","4117958050032200000","4117958050022700001","0163300000000000"); // INT7
  } else if (trainConfig == 3444){
    cuts.AddCutMergedCalo("8008e123","4117958050032200000","4117958050022700001","0163300000000000"); // EG2+DG2
    cuts.AddCutMergedCalo("8008d123","4117958050032200000","4117958050022700001","0163300000000000"); // EG1+DG1
  } else if (trainConfig == 3445){
    cuts.AddCutMergedCalo("8009c123","4117958050032200000","4117958050022700001","0163300000000000"); // EJ2+DJ2
    cuts.AddCutMergedCalo("8009b123","4117958050032200000","4117958050022700001","0163300000000000"); // EJ1+DJ1



  // no TM configs for mcp2
  // std cut with TB (Nico) NL
  } else if (trainConfig == 3501){
    cuts.AddCutMergedCalo("80010123","4117900050032200000","4117900050022700001","0163300000000000"); // INT7
  } else if (trainConfig == 3502){
    cuts.AddCutMergedCalo("8008e123","4117900050032200000","4117900050022700001","0163300000000000"); // EG2+DG2
    cuts.AddCutMergedCalo("8008d123","4117900050032200000","4117900050022700001","0163300000000000"); // EG1+DG1
  } else if (trainConfig == 3504){
    cuts.AddCutMergedCalo("80010123","4117965050032200000","4117965050022700001","0163300000000000"); // INT7
  } else if (trainConfig == 3505){
    cuts.AddCutMergedCalo("8008e123","4117965050032200000","4117965050022700001","0163300000000000"); // EG2+DG2
    cuts.AddCutMergedCalo("8008d123","4117965050032200000","4117965050022700001","0163300000000000"); // EG1+DG1
  } else if (trainConfig == 3506){
    cuts.AddCutMergedCalo("80010123","4117901050032200000","4117901050022700001","0163300000000000"); // INT7
  } else if (trainConfig == 3507){
    cuts.AddCutMergedCalo("8008e123","4117901050032200000","4117901050022700001","0163300000000000"); // EG2+DG2
    cuts.AddCutMergedCalo("8008d123","4117901050032200000","4117901050022700001","0163300000000000"); // EG1+DG1

  // std cut with TB+FineTuning NL
  } else if (trainConfig == 3510){
    cuts.AddCutMergedCalo("80010123","4117931050032200000","4117931050022700001","0163300000000000"); // INT7
  } else if (trainConfig == 3511){
    cuts.AddCutMergedCalo("8008e123","4117931050032200000","4117931050022700001","0163300000000000"); // EG2+DG2
    cuts.AddCutMergedCalo("8008d123","4117931050032200000","4117931050022700001","0163300000000000"); // EG1+DG1
  } else if (trainConfig == 3513){
    cuts.AddCutMergedCalo("80010123","4117932050032200000","4117932050022700001","0163300000000000"); // INT7
  } else if (trainConfig == 3514){
    cuts.AddCutMergedCalo("8008e123","4117932050032200000","4117932050022700001","0163300000000000"); // EG2+DG2
    cuts.AddCutMergedCalo("8008d123","4117932050032200000","4117932050022700001","0163300000000000"); // EG1+DG1

  } else if (trainConfig == 3520){
    cuts.AddCutMergedCalo("80010123","4117933050032200000","4117933050022700001","0163300000000000"); // INT7
  } else if (trainConfig == 3521){
    cuts.AddCutMergedCalo("8008e123","4117933050032200000","4117933050022700001","0163300000000000"); // EG2+DG2
    cuts.AddCutMergedCalo("8008d123","4117933050032200000","4117933050022700001","0163300000000000"); // EG1+DG1
  } else if (trainConfig == 3523){
    cuts.AddCutMergedCalo("80010123","4117934050032200000","4117934050022700001","0163300000000000"); // INT7
  } else if (trainConfig == 3524){
    cuts.AddCutMergedCalo("8008e123","4117934050032200000","4117934050022700001","0163300000000000"); // EG2+DG2
    cuts.AddCutMergedCalo("8008d123","4117934050032200000","4117934050022700001","0163300000000000"); // EG1+DG1
  } else if (trainConfig == 3526){ // config for V1 clusterizer
    cuts.AddCutMergedCalo("80010123","4117957050032200000","4117957050022700002","0163300000000000"); // INT7
  } else if (trainConfig == 3527){ // config for V1 clusterizer
    cuts.AddCutMergedCalo("8008e123","4117957050032200000","4117957050022700002","0163300000000000"); // EG2+DG2
    cuts.AddCutMergedCalo("8008d123","4117957050032200000","4117957050022700002","0163300000000000"); // EG1+DG1

  // systematics MB
  } else if (trainConfig == 3600){
    cuts.AddCutMergedCalo("80010123","4117947050032200000","4117947050022700001","0163300000000000"); // default
  } else if (trainConfig == 3602){ // varied M02 part 1
    cuts.AddCutMergedCalo("80010123","4117947050032200000","4117947050022600001","0163300000000000"); // min M02 = 0.3
    cuts.AddCutMergedCalo("80010123","4117947050032200000","4117947050022c00001","0163300000000000"); // min M02 = 0.29
    cuts.AddCutMergedCalo("80010123","4117947050032200000","4117947050022b00001","0163300000000000"); // min M02 = 0.28
    cuts.AddCutMergedCalo("80010123","4117947050032200000","4117947050022800001","0163300000000000"); // min M02 = 0.25
    cuts.AddCutMergedCalo("80010123","4117947050032000000","4117947050022000001","0163300000000000"); // open M02
  } else if (trainConfig == 3603){ // |eta| < 0.3, y < 0.3
    cuts.AddCutMergedCalo("80010123","4667947050032200000","4667947050022700001","0163300000000000"); // diff eta/rap cuts
  } else if (trainConfig == 3604){ // NL vars
    cuts.AddCutMergedCalo("80010123","4117900050032200000","4117900050022700001","0163300000000000"); // default
    cuts.AddCutMergedCalo("80010123","4117947050032200000","4117947050022700001","0163300000000000"); // 47
    cuts.AddCutMergedCalo("80010123","4117948050032200000","4117948050022700001","0163300000000000"); // 48
    cuts.AddCutMergedCalo("80010123","4117957050032200000","4117957050022700001","0163300000000000"); // 57
    cuts.AddCutMergedCalo("80010123","4117958050032200000","4117958050022700001","0163300000000000"); // 58
  } else if (trainConfig == 3605){ // timing
    cuts.AddCutMergedCalo("80010123","4117947040032200000","4117947040022700001","0163300000000000"); //
    cuts.AddCutMergedCalo("80010123","4117947070032200000","4117947070022700001","0163300000000000"); //
    cuts.AddCutMergedCalo("80010123","41179470a0032200000","41179470a0022700001","0163300000000000"); //
    cuts.AddCutMergedCalo("80010123","4117947000032200000","4117947000022700001","0163300000000000"); //

  // systematics EG2
  } else if (trainConfig == 3700){
    cuts.AddCutMergedCalo("8008e123","4117947050032200000","4117947050022700001","0163300000000000"); // default
  } else if (trainConfig == 3702){ // varied M02 part 1
    cuts.AddCutMergedCalo("8008e123","4117947050032200000","4117947050022600001","0163300000000000"); // min M02 = 0.3
    cuts.AddCutMergedCalo("8008e123","4117947050032200000","4117947050022c00001","0163300000000000"); // min M02 = 0.29
    cuts.AddCutMergedCalo("8008e123","4117947050032200000","4117947050022b00001","0163300000000000"); // min M02 = 0.28
    cuts.AddCutMergedCalo("8008e123","4117947050032200000","4117947050022800001","0163300000000000"); // min M02 = 0.25
    cuts.AddCutMergedCalo("8008e123","4117947050032000000","4117947050022000001","0163300000000000"); // open M02
  } else if (trainConfig == 3703){ // |eta| < 0.3, y < 0.3
    cuts.AddCutMergedCalo("8008e123","4667947050032200000","4667947050022700001","0163300000000000"); // diff eta/rap cuts
  } else if (trainConfig == 3704){ // NL vars
    cuts.AddCutMergedCalo("8008e123","4117900050032200000","4117900050022700001","0163300000000000"); // default
    cuts.AddCutMergedCalo("8008e123","4117947050032200000","4117947050022700001","0163300000000000"); // 47
    cuts.AddCutMergedCalo("8008e123","4117948050032200000","4117948050022700001","0163300000000000"); // 48
    cuts.AddCutMergedCalo("8008e123","4117957050032200000","4117957050022700001","0163300000000000"); // 57
    cuts.AddCutMergedCalo("8008e123","4117958050032200000","4117958050022700001","0163300000000000"); // 58
  } else if (trainConfig == 3705){ // timing
    cuts.AddCutMergedCalo("8008e123","4117947040032200000","4117947040022700001","0163300000000000"); //
    cuts.AddCutMergedCalo("8008e123","4117947070032200000","4117947070022700001","0163300000000000"); //
    cuts.AddCutMergedCalo("8008e123","41179470a0032200000","41179470a0022700001","0163300000000000"); //
    cuts.AddCutMergedCalo("8008e123","4117947000032200000","4117947000022700001","0163300000000000"); //

  // systematics pPb8TeV w/o TM for muon calo pass
  // MB configs
  } else if (trainConfig == 3800){
    cuts.AddCutMergedCalo("80010123","4117932050032200000","4117932050022700001","0163300000000000"); // std
    cuts.AddCutMergedCalo("80010123","1111132050032200000","1111132050022700001","0163300000000000"); // std
  } else if (trainConfig == 3801){ // varied M02 part 1
    cuts.AddCutMergedCalo("80010123","4117932050032200000","4117932050022600001","0163300000000000"); // min 0.3
    cuts.AddCutMergedCalo("80010123","4117932050032200000","4117932050022900001","0163300000000000"); // min 0.1
  } else if (trainConfig == 3802){ // varied M02 part 1
    cuts.AddCutMergedCalo("80010123","4117932050032200000","4117932050022d00001","0163300000000000"); // min 0.33
    cuts.AddCutMergedCalo("80010123","4117932050032200000","4117932050022e00001","0163300000000000"); // min 0.36
  } else if (trainConfig == 3803){ // varied M02 part 2
    cuts.AddCutMergedCalo("80010123","4117932050032200000","4117932050022800001","0163300000000000"); // min 0.25
    cuts.AddCutMergedCalo("80010123","4117932050032200000","4117932050022a00001","0163300000000000"); // min 0.26
  } else if (trainConfig == 3804){ // NL vars
    cuts.AddCutMergedCalo("80010123","4117933050032200000","4117933050022700001","0163300000000000"); // 32
    cuts.AddCutMergedCalo("80010123","4117934050032200000","4117934050022700001","0163300000000000"); // 33
  } else if (trainConfig == 3805){ // varied M02 part 1
    cuts.AddCutMergedCalo("80010123","4117938050032200000","4117938050022700001","0163300000000000"); // 38
    cuts.AddCutMergedCalo("80010123","4117939050032200000","4117939050022700001","0163300000000000"); // 39
  } else if (trainConfig == 3806){ // distance to bad channel
    cuts.AddCutMergedCalo("80010123","4117932050032200000","4117932150022700001","0163300000000000"); // 1 <= coll+row
    cuts.AddCutMergedCalo("80010123","4117932050032200000","4117932250022700001","0163300000000000"); // 2 <= coll+row
  } else if (trainConfig == 3807){ // varied M02 part 1
    cuts.AddCutMergedCalo("80010123","4117932050032200000","4117932550022700001","0163300000000000"); // 1 <= row || coll || 0.5*(coll+row)
    cuts.AddCutMergedCalo("80010123","4117932050032200000","4117932650022700001","0163300000000000"); // 2 <= row || coll || 0.5*(coll+row)
  } else if (trainConfig == 3808){ // timing
    cuts.AddCutMergedCalo("80010123","4117932040032200000","4117932040022700001","0163300000000000"); // 10ns
    cuts.AddCutMergedCalo("80010123","4117932070032200000","4117932070022700001","0163300000000000"); // 30ns
  } else if (trainConfig == 3809){ // timing
    cuts.AddCutMergedCalo("80010123","4117932090032200000","4117932090022700001","0163300000000000"); // -20 to 25
    cuts.AddCutMergedCalo("80010123","4117932050032200000","4117932050022700001","0163300000000000"); // -50 to 50
  } else if (trainConfig == 3810){ // timing
    cuts.AddCutMergedCalo("80010123","4117932030032200000","4117932030022700001","0163300000000000"); // -200 to 200
    cuts.AddCutMergedCalo("80010123","4117932020032200000","4117932020022700001","0163300000000000"); // -500 to 500
  } else if (trainConfig == 3811){
    cuts.AddCutMergedCalo("80010123","4117932050032200000","4117932050022700001","0163200000000000"); // rapidity < 0.7
    cuts.AddCutMergedCalo("80010123","4117932050032200000","4117932050022700001","0163100000000000"); // rapidity < 0.8
  } else if (trainConfig == 3812){ // varied M02 part 1
    cuts.AddCutMergedCalo("80010123","4117932050032200000","4117932050022g00001","0163300000000000"); // min 0.24
    cuts.AddCutMergedCalo("80010123","4117932050032200000","4117932050022h00001","0163300000000000"); // min 0.23
    cuts.AddCutMergedCalo("80010123","4117932050032200000","4117932050022j00001","0163300000000000"); // min 0.21
  } else if (trainConfig == 3813){ // varied M02 part 2
    cuts.AddCutMergedCalo("80010123","4117932050032200000","4117932050022k00001","0163300000000000"); // min 0.20
    cuts.AddCutMergedCalo("80010123","4117932050032200000","4117932050022l00001","0163300000000000"); // min 0.19
    cuts.AddCutMergedCalo("80010123","4117932050032200000","4117932050022m00001","0163300000000000"); // min 0.18

   // EG2 configs
  } else if (trainConfig == 3820){
    cuts.AddCutMergedCalo("8008e123","4117932050032200000","4117932050022700001","0163300000000000"); // std
    cuts.AddCutMergedCalo("80085123","1111132050032200000","1111132050022700001","0163300000000000"); // std
  } else if (trainConfig == 3821){ // varied M02 part 1
    cuts.AddCutMergedCalo("8008e123","4117932050032200000","4117932050022600001","0163300000000000"); // min 0.3
    cuts.AddCutMergedCalo("8008e123","4117932050032200000","4117932050022900001","0163300000000000"); // min 0.1
  } else if (trainConfig == 3822){ // varied M02 part 1
    cuts.AddCutMergedCalo("8008e123","4117932050032200000","4117932050022d00001","0163300000000000"); // min 0.33
    cuts.AddCutMergedCalo("8008e123","4117932050032200000","4117932050022e00001","0163300000000000"); // min 0.36
  } else if (trainConfig == 3823){ // varied M02 part 2
    cuts.AddCutMergedCalo("8008e123","4117932050032200000","4117932050022800001","0163300000000000"); // min 0.25
    cuts.AddCutMergedCalo("8008e123","4117932050032200000","4117932050022a00001","0163300000000000"); // min 0.26
  } else if (trainConfig == 3824){ // NL vars
    cuts.AddCutMergedCalo("8008e123","4117933050032200000","4117933050022700001","0163300000000000"); // 32
    cuts.AddCutMergedCalo("8008e123","4117934050032200000","4117934050022700001","0163300000000000"); // 33
  } else if (trainConfig == 3825){ // varied M02 part 1
    cuts.AddCutMergedCalo("8008e123","4117938050032200000","4117938050022700001","0163300000000000"); // 38
    cuts.AddCutMergedCalo("8008e123","4117939050032200000","4117939050022700001","0163300000000000"); // 39
  } else if (trainConfig == 3826){ // distance to bad channel
    cuts.AddCutMergedCalo("8008e123","4117932050032200000","4117932150022700001","0163300000000000"); // 1 <= coll+row
    cuts.AddCutMergedCalo("8008e123","4117932050032200000","4117932250022700001","0163300000000000"); // 2 <= coll+row
  } else if (trainConfig == 3827){ // varied M02 part 1
    cuts.AddCutMergedCalo("8008e123","4117932050032200000","4117932550022700001","0163300000000000"); // 1 <= row || coll || 0.5*(coll+row)
    cuts.AddCutMergedCalo("8008e123","4117932050032200000","4117932650022700001","0163300000000000"); // 2 <= row || coll || 0.5*(coll+row)
  } else if (trainConfig == 3828){ // timing
    cuts.AddCutMergedCalo("8008e123","4117932040032200000","4117932040022700001","0163300000000000"); // 10ns
    cuts.AddCutMergedCalo("8008e123","4117932070032200000","4117932070022700001","0163300000000000"); // 30ns
  } else if (trainConfig == 3829){ // timing
    cuts.AddCutMergedCalo("8008e123","4117932090032200000","4117932090022700001","0163300000000000"); // -20 to 25
    cuts.AddCutMergedCalo("8008e123","4117932050032200000","4117932050022700001","0163300000000000"); // -50 to 50
  } else if (trainConfig == 3830){ // timing
    cuts.AddCutMergedCalo("8008e123","4117932030032200000","4117932030022700001","0163300000000000"); // -200 to 200
    cuts.AddCutMergedCalo("8008e123","4117932020032200000","4117932020022700001","0163300000000000"); // -500 to 500
  } else if (trainConfig == 3831){
    cuts.AddCutMergedCalo("8008e123","4117932050032200000","4117932050022700001","0163200000000000"); // rapidity < 0.7
    cuts.AddCutMergedCalo("8008e123","4117932050032200000","4117932050022700001","0163100000000000"); // rapidity < 0.8
  } else if (trainConfig == 3832){ // varied M02 part 1
    cuts.AddCutMergedCalo("8008e123","4117932050032200000","4117932050022g00001","0163300000000000"); // min 0.24
    cuts.AddCutMergedCalo("8008e123","4117932050032200000","4117932050022h00001","0163300000000000"); // min 0.23
    cuts.AddCutMergedCalo("8008e123","4117932050032200000","4117932050022j00001","0163300000000000"); // min 0.21
  } else if (trainConfig == 3833){ // varied M02 part 2
    cuts.AddCutMergedCalo("8008e123","4117932050032200000","4117932050022k00001","0163300000000000"); // min 0.20
    cuts.AddCutMergedCalo("8008e123","4117932050032200000","4117932050022l00001","0163300000000000"); // min 0.19
    cuts.AddCutMergedCalo("8008e123","4117932050032200000","4117932050022m00001","0163300000000000"); // min 0.18

  // EG1 configs
  } else if (trainConfig == 3840){
    cuts.AddCutMergedCalo("8008d123","4117932050032200000","4117932050022700001","0163300000000000"); // std
    cuts.AddCutMergedCalo("80083123","1111132050032200000","1111132050022700001","0163300000000000"); // std
  } else if (trainConfig == 3841){ // varied M02 part 1
    cuts.AddCutMergedCalo("8008d123","4117932050032200000","4117932050022600001","0163300000000000"); // min 0.3
    cuts.AddCutMergedCalo("8008d123","4117932050032200000","4117932050022900001","0163300000000000"); // min 0.1
  } else if (trainConfig == 3842){ // varied M02 part 1
    cuts.AddCutMergedCalo("8008d123","4117932050032200000","4117932050022d00001","0163300000000000"); // min 0.33
    cuts.AddCutMergedCalo("8008d123","4117932050032200000","4117932050022e00001","0163300000000000"); // min 0.36
  } else if (trainConfig == 3843){ // varied M02 part 2
    cuts.AddCutMergedCalo("8008d123","4117932050032200000","4117932050022800001","0163300000000000"); // min 0.25
    cuts.AddCutMergedCalo("8008d123","4117932050032200000","4117932050022a00001","0163300000000000"); // min 0.26
  } else if (trainConfig == 3844){ // NL vars
    cuts.AddCutMergedCalo("8008d123","4117933050032200000","4117933050022700001","0163300000000000"); // 32
    cuts.AddCutMergedCalo("8008d123","4117934050032200000","4117934050022700001","0163300000000000"); // 33
  } else if (trainConfig == 3845){ // varied M02 part 1
    cuts.AddCutMergedCalo("8008d123","4117938050032200000","4117938050022700001","0163300000000000"); // 38
    cuts.AddCutMergedCalo("8008d123","4117939050032200000","4117939050022700001","0163300000000000"); // 39
  } else if (trainConfig == 3846){ // distance to bad channel
    cuts.AddCutMergedCalo("8008d123","4117932050032200000","4117932150022700001","0163300000000000"); // 1 <= coll+row
    cuts.AddCutMergedCalo("8008d123","4117932050032200000","4117932250022700001","0163300000000000"); // 2 <= coll+row
  } else if (trainConfig == 3847){ // varied M02 part 1
    cuts.AddCutMergedCalo("8008d123","4117932050032200000","4117932550022700001","0163300000000000"); // 1 <= row || coll || 0.5*(coll+row)
    cuts.AddCutMergedCalo("8008d123","4117932050032200000","4117932650022700001","0163300000000000"); // 2 <= row || coll || 0.5*(coll+row)
  } else if (trainConfig == 3848){ // timing
    cuts.AddCutMergedCalo("8008d123","4117932040032200000","4117932040022700001","0163300000000000"); // 10ns
    cuts.AddCutMergedCalo("8008d123","4117932070032200000","4117932070022700001","0163300000000000"); // 30ns
  } else if (trainConfig == 3849){ // timing
    cuts.AddCutMergedCalo("8008d123","4117932090032200000","4117932090022700001","0163300000000000"); // -20 to 25
    cuts.AddCutMergedCalo("8008d123","4117932050032200000","4117932050022700001","0163300000000000"); // -50 to 50
  } else if (trainConfig == 3850){ // timing
    cuts.AddCutMergedCalo("8008d123","4117932030032200000","4117932030022700001","0163300000000000"); // -200 to 200
    cuts.AddCutMergedCalo("8008d123","4117932020032200000","4117932020022700001","0163300000000000"); // -500 to 500
  } else if (trainConfig == 3851){
    cuts.AddCutMergedCalo("8008d123","4117932050032200000","4117932050022700001","0163200000000000"); // rapidity < 0.7
    cuts.AddCutMergedCalo("8008d123","4117932050032200000","4117932050022700001","0163100000000000"); // rapidity < 0.8
  } else if (trainConfig == 3852){ // varied M02 part 1
    cuts.AddCutMergedCalo("8008d123","4117932050032200000","4117932050022g00001","0163300000000000"); // min 0.24
    cuts.AddCutMergedCalo("8008d123","4117932050032200000","4117932050022h00001","0163300000000000"); // min 0.23
    cuts.AddCutMergedCalo("8008d123","4117932050032200000","4117932050022j00001","0163300000000000"); // min 0.21
  } else if (trainConfig == 3853){ // varied M02 part 2
    cuts.AddCutMergedCalo("8008d123","4117932050032200000","4117932050022k00001","0163300000000000"); // min 0.20
    cuts.AddCutMergedCalo("8008d123","4117932050032200000","4117932050022l00001","0163300000000000"); // min 0.19
    cuts.AddCutMergedCalo("8008d123","4117932050032200000","4117932050022m00001","0163300000000000"); // min 0.18


  } else if (trainConfig == 4000){ // for dec. gamma MC
    cuts.AddCutMergedCalo("80010103","1111131050032200000","1111131050022700001","0163300000000000"); // std
    cuts.AddCutMergedCalo("80085103","1111131050032200000","1111131050022700001","0163300000000000"); // std
    cuts.AddCutMergedCalo("80083103","1111131050032200000","1111131050022700001","0163300000000000"); // std
  } else if (trainConfig == 4001){ // for JJ MC
    cuts.AddCutMergedCalo("80010123","1111131050032200000","1111131050022700001","0163300000000000"); // std
    cuts.AddCutMergedCalo("80085123","1111131050032200000","1111131050022700001","0163300000000000"); // std
    cuts.AddCutMergedCalo("80083123","1111131050032200000","1111131050022700001","0163300000000000"); // std
  } else if (trainConfig == 4002){ // cent dep.
    cuts.AddCutMergedCalo("a0110103","1111131050032200000","1111131050022700001","0163300000000000"); // std
    cuts.AddCutMergedCalo("a1210103","1111131050032200000","1111131050022700001","0163300000000000"); // std
    cuts.AddCutMergedCalo("81210103","1111131050032200000","1111131050022700001","0163300000000000"); // std
  } else if (trainConfig == 4003){ // cent dep.
    cuts.AddCutMergedCalo("82410103","1111131050032200000","1111131050022700001","0163300000000000"); // std
    cuts.AddCutMergedCalo("84610103","1111131050032200000","1111131050022700001","0163300000000000"); // std
    cuts.AddCutMergedCalo("86810103","1111131050032200000","1111131050022700001","0163300000000000"); // std
    cuts.AddCutMergedCalo("88010103","1111131050032200000","1111131050022700001","0163300000000000"); // std
  } else if (trainConfig == 4004){ // cent dep.
    cuts.AddCutMergedCalo("a0185103","1111131050032200000","1111131050022700001","0163300000000000"); // std
    cuts.AddCutMergedCalo("a1285103","1111131050032200000","1111131050022700001","0163300000000000"); // std
    cuts.AddCutMergedCalo("81285103","1111131050032200000","1111131050022700001","0163300000000000"); // std
  } else if (trainConfig == 4005){ // cent dep.
    cuts.AddCutMergedCalo("82485103","1111131050032200000","1111131050022700001","0163300000000000"); // std
    cuts.AddCutMergedCalo("84685103","1111131050032200000","1111131050022700001","0163300000000000"); // std
    cuts.AddCutMergedCalo("86885103","1111131050032200000","1111131050022700001","0163300000000000"); // std
    cuts.AddCutMergedCalo("88085103","1111131050032200000","1111131050022700001","0163300000000000"); // std
  } else if (trainConfig == 4006){ // cent dep.
    cuts.AddCutMergedCalo("a0183103","1111131050032200000","1111131050022700001","0163300000000000"); // std
    cuts.AddCutMergedCalo("a1283103","1111131050032200000","1111131050022700001","0163300000000000"); // std
    cuts.AddCutMergedCalo("81283103","1111131050032200000","1111131050022700001","0163300000000000"); // std
  } else if (trainConfig == 4007){ // cent dep.
    cuts.AddCutMergedCalo("82483103","1111131050032200000","1111131050022700001","0163300000000000"); // std
    cuts.AddCutMergedCalo("84683103","1111131050032200000","1111131050022700001","0163300000000000"); // std
    cuts.AddCutMergedCalo("86883103","1111131050032200000","1111131050022700001","0163300000000000"); // std
    cuts.AddCutMergedCalo("88083103","1111131050032200000","1111131050022700001","0163300000000000"); // std
  } else if (trainConfig == 4010){ // for dec. gamma MC
    cuts.AddCutMergedCalo("80010103","3885531050032200000","3885531050022700001","0163300000000000"); // std
    cuts.AddCutMergedCalo("8008b103","3885531050032200000","3885531050022700001","0163300000000000"); // std
    cuts.AddCutMergedCalo("80089103","3885531050032200000","3885531050022700001","0163300000000000"); // std
  } else if (trainConfig == 4011){ // for JJ MC
    cuts.AddCutMergedCalo("80010123","3885531050032200000","3885531050022700001","0163300000000000"); // std
    cuts.AddCutMergedCalo("8008b123","3885531050032200000","3885531050022700001","0163300000000000"); // std
    cuts.AddCutMergedCalo("80089123","3885531050032200000","3885531050022700001","0163300000000000"); // std
  } else if (trainConfig == 4020){ // for dec. gamma MC
    cuts.AddCutMergedCalo("80010103","4117931050032200000","4117931050022700001","0163300000000000"); // std
    cuts.AddCutMergedCalo("8008e103","4117931050032200000","4117931050022700001","0163300000000000"); // std
    cuts.AddCutMergedCalo("8008d103","4117931050032200000","4117931050022700001","0163300000000000"); // std
  } else if (trainConfig == 4021){ // for JJ MC
    cuts.AddCutMergedCalo("80010123","4117931050032200000","4117931050022700001","0163300000000000"); // std
    cuts.AddCutMergedCalo("8008e123","4117931050032200000","4117931050022700001","0163300000000000"); // std
    cuts.AddCutMergedCalo("8008d123","4117931050032200000","4117931050022700001","0163300000000000"); // std

  } else if (trainConfig == 4022){ // for JJ MC, T0-based
    cuts.AddCutMergedCalo("80010123","4117939050032200000","4117939050022700001","0163300000000000"); // std
    cuts.AddCutMergedCalo("8008e123","4117939050032200000","4117939050022700001","0163300000000000"); // std
    cuts.AddCutMergedCalo("8008d123","4117939050032200000","4117939050022700001","0163300000000000"); // std
  } else if (trainConfig == 4023){ // for JJ MC, T0-based
    cuts.AddCutMergedCalo("80011123","4117939050032200000","4117939050022700001","0163300000000000"); // std
    cuts.AddCutMergedCalo("8008g123","4117939050032200000","4117939050022700001","0163300000000000"); // std
    cuts.AddCutMergedCalo("8008f123","4117939050032200000","4117939050022700001","0163300000000000"); // std
  } else if (trainConfig == 4024){ // for JJ MC, T0-based
    cuts.AddCutMergedCalo("80011123","4117931050032200000","4117931050022700001","0163300000000000"); // std
    cuts.AddCutMergedCalo("8008g123","4117931050032200000","4117931050022700001","0163300000000000"); // std
    cuts.AddCutMergedCalo("8008f123","4117931050032200000","4117931050022700001","0163300000000000"); // std


  // systematics pPb5TeV w/o TM
  // MB configs
  } else if (trainConfig == 4100){
    cuts.AddCutMergedCalo("80010103","1111132050032200000","1111132050022700001","0163300000000000"); // std
  } else if (trainConfig == 4101){ // varied M02 part 1
    cuts.AddCutMergedCalo("80010103","1111132050032200000","1111132050022600001","0163300000000000"); // min 0.3
    cuts.AddCutMergedCalo("80010103","1111132050032200000","1111132050022900001","0163300000000000"); // min 0.1
  } else if (trainConfig == 4102){ // varied M02 part 1
    cuts.AddCutMergedCalo("80010103","1111132050032200000","1111132050022d00001","0163300000000000"); // min 0.33
    cuts.AddCutMergedCalo("80010103","1111132050032200000","1111132050022e00001","0163300000000000"); // min 0.36
  } else if (trainConfig == 4103){ // varied M02 part 2
    cuts.AddCutMergedCalo("80010103","1111132050032200000","1111132050022800001","0163300000000000"); // min 0.25
    cuts.AddCutMergedCalo("80010103","1111132050032200000","1111132050022a00001","0163300000000000"); // min 0.26
  } else if (trainConfig == 4104){ // NL vars
    cuts.AddCutMergedCalo("80010103","1111133050032200000","1111133050022700001","0163300000000000"); // 32
    cuts.AddCutMergedCalo("80010103","1111134050032200000","1111134050022700001","0163300000000000"); // 33
  } else if (trainConfig == 4105){ // varied M02 part 1
    cuts.AddCutMergedCalo("80010103","1111138050032200000","1111138050022700001","0163300000000000"); // 38
    cuts.AddCutMergedCalo("80010103","1111139050032200000","1111139050022700001","0163300000000000"); // 39
  } else if (trainConfig == 4106){ // distance to bad channel
    cuts.AddCutMergedCalo("80010103","1111132050032200000","1111132150022700001","0163300000000000"); // 1 <= coll+row
    cuts.AddCutMergedCalo("80010103","1111132050032200000","1111132250022700001","0163300000000000"); // 2 <= coll+row
  } else if (trainConfig == 4107){ // varied M02 part 1
    cuts.AddCutMergedCalo("80010103","1111132050032200000","1111132550022700001","0163300000000000"); // 1 <= row || coll || 0.5*(coll+row)
    cuts.AddCutMergedCalo("80010103","1111132050032200000","1111132650022700001","0163300000000000"); // 2 <= row || coll || 0.5*(coll+row)
  } else if (trainConfig == 4108){ // timing
    cuts.AddCutMergedCalo("80010103","1111132040032200000","1111132040022700001","0163300000000000"); // 10ns
    cuts.AddCutMergedCalo("80010103","1111132070032200000","1111132070022700001","0163300000000000"); // 30ns
  } else if (trainConfig == 4109){ // timing
    cuts.AddCutMergedCalo("80010103","1111132090032200000","1111132090022700001","0163300000000000"); // -20 to 25
    cuts.AddCutMergedCalo("80010103","1111132050032200000","1111132050022700001","0163300000000000"); // -50 to 50
  } else if (trainConfig == 4110){ // timing
    cuts.AddCutMergedCalo("80010103","1111132030032200000","1111132030022700001","0163300000000000"); // -200 to 200
    cuts.AddCutMergedCalo("80010103","1111132020032200000","1111132020022700001","0163300000000000"); // -500 to 500
  } else if (trainConfig == 4111){
    cuts.AddCutMergedCalo("80010103","1111132050032200000","1111132050022700001","0163200000000000"); // rapidity < 0.7
    cuts.AddCutMergedCalo("80010103","1111132050032200000","1111132050022700001","0163100000000000"); // rapidity < 0.8

   // EG2 configs
  } else if (trainConfig == 4120){
    cuts.AddCutMergedCalo("80085103","1111132050032200000","1111132050022700001","0163300000000000"); // std
  } else if (trainConfig == 4121){ // varied M02 part 1
    cuts.AddCutMergedCalo("80085103","1111132050032200000","1111132050022600001","0163300000000000"); // min 0.3
    cuts.AddCutMergedCalo("80085103","1111132050032200000","1111132050022900001","0163300000000000"); // min 0.1
  } else if (trainConfig == 4122){ // varied M02 part 1
    cuts.AddCutMergedCalo("80085103","1111132050032200000","1111132050022d00001","0163300000000000"); // min 0.33
    cuts.AddCutMergedCalo("80085103","1111132050032200000","1111132050022e00001","0163300000000000"); // min 0.36
  } else if (trainConfig == 4123){ // varied M02 part 2
    cuts.AddCutMergedCalo("80085103","1111132050032200000","1111132050022800001","0163300000000000"); // min 0.25
    cuts.AddCutMergedCalo("80085103","1111132050032200000","1111132050022a00001","0163300000000000"); // min 0.26
  } else if (trainConfig == 4124){ // NL vars
    cuts.AddCutMergedCalo("80085103","1111133050032200000","1111133050022700001","0163300000000000"); // 32
    cuts.AddCutMergedCalo("80085103","1111134050032200000","1111134050022700001","0163300000000000"); // 33
  } else if (trainConfig == 4125){ // varied M02 part 1
    cuts.AddCutMergedCalo("80085103","1111138050032200000","1111138050022700001","0163300000000000"); // 38
    cuts.AddCutMergedCalo("80085103","1111139050032200000","1111139050022700001","0163300000000000"); // 39
  } else if (trainConfig == 4126){ // distance to bad channel
    cuts.AddCutMergedCalo("80085103","1111132050032200000","1111132150022700001","0163300000000000"); // 1 <= coll+row
    cuts.AddCutMergedCalo("80085103","1111132050032200000","1111132250022700001","0163300000000000"); // 2 <= coll+row
  } else if (trainConfig == 4127){ // varied M02 part 1
    cuts.AddCutMergedCalo("80085103","1111132050032200000","1111132550022700001","0163300000000000"); // 1 <= row || coll || 0.5*(coll+row)
    cuts.AddCutMergedCalo("80085103","1111132050032200000","1111132650022700001","0163300000000000"); // 2 <= row || coll || 0.5*(coll+row)
  } else if (trainConfig == 4128){ // timing
    cuts.AddCutMergedCalo("80085103","1111132040032200000","1111132040022700001","0163300000000000"); // 10ns
    cuts.AddCutMergedCalo("80085103","1111132070032200000","1111132070022700001","0163300000000000"); // 30ns
  } else if (trainConfig == 4129){ // timing
    cuts.AddCutMergedCalo("80085103","1111132090032200000","1111132090022700001","0163300000000000"); // -20 to 25
    cuts.AddCutMergedCalo("80085103","1111132050032200000","1111132050022700001","0163300000000000"); // -50 to 50
  } else if (trainConfig == 4130){ // timing
    cuts.AddCutMergedCalo("80085103","1111132030032200000","1111132030022700001","0163300000000000"); // -200 to 200
    cuts.AddCutMergedCalo("80085103","1111132020032200000","1111132020022700001","0163300000000000"); // -500 to 500
  } else if (trainConfig == 4131){
    cuts.AddCutMergedCalo("80085103","1111132050032200000","1111132050022700001","0163200000000000"); // rapidity < 0.7
    cuts.AddCutMergedCalo("80085103","1111132050032200000","1111132050022700001","0163100000000000"); // rapidity < 0.8

  // EG1 configs
  } else if (trainConfig == 4140){
    cuts.AddCutMergedCalo("80083103","1111132050032200000","1111132050022700001","0163300000000000"); // std
  } else if (trainConfig == 4141){ // varied M02 part 1
    cuts.AddCutMergedCalo("80083103","1111132050032200000","1111132050022600001","0163300000000000"); // min 0.3
    cuts.AddCutMergedCalo("80083103","1111132050032200000","1111132050022900001","0163300000000000"); // min 0.1
  } else if (trainConfig == 4142){ // varied M02 part 1
    cuts.AddCutMergedCalo("80083103","1111132050032200000","1111132050022d00001","0163300000000000"); // min 0.33
    cuts.AddCutMergedCalo("80083103","1111132050032200000","1111132050022e00001","0163300000000000"); // min 0.36
  } else if (trainConfig == 4143){ // varied M02 part 2
    cuts.AddCutMergedCalo("80083103","1111132050032200000","1111132050022800001","0163300000000000"); // min 0.25
    cuts.AddCutMergedCalo("80083103","1111132050032200000","1111132050022a00001","0163300000000000"); // min 0.26
  } else if (trainConfig == 4144){ // NL vars
    cuts.AddCutMergedCalo("80083103","1111133050032200000","1111133050022700001","0163300000000000"); // 32
    cuts.AddCutMergedCalo("80083103","1111134050032200000","1111134050022700001","0163300000000000"); // 33
  } else if (trainConfig == 4145){ // varied M02 part 1
    cuts.AddCutMergedCalo("80083103","1111138050032200000","1111138050022700001","0163300000000000"); // 38
    cuts.AddCutMergedCalo("80083103","1111139050032200000","1111139050022700001","0163300000000000"); // 39
  } else if (trainConfig == 4146){ // distance to bad channel
    cuts.AddCutMergedCalo("80083103","1111132050032200000","1111132150022700001","0163300000000000"); // 1 <= coll+row
    cuts.AddCutMergedCalo("80083103","1111132050032200000","1111132250022700001","0163300000000000"); // 2 <= coll+row
  } else if (trainConfig == 4147){ // varied M02 part 1
    cuts.AddCutMergedCalo("80083103","1111132050032200000","1111132550022700001","0163300000000000"); // 1 <= row || coll || 0.5*(coll+row)
    cuts.AddCutMergedCalo("80083103","1111132050032200000","1111132650022700001","0163300000000000"); // 2 <= row || coll || 0.5*(coll+row)
  } else if (trainConfig == 4148){ // timing
    cuts.AddCutMergedCalo("80083103","1111132040032200000","1111132040022700001","0163300000000000"); // 10ns
    cuts.AddCutMergedCalo("80083103","1111132070032200000","1111132070022700001","0163300000000000"); // 30ns
  } else if (trainConfig == 4149){ // timing
    cuts.AddCutMergedCalo("80083103","1111132090032200000","1111132090022700001","0163300000000000"); // -20 to 25
    cuts.AddCutMergedCalo("80083103","1111132050032200000","1111132050022700001","0163300000000000"); // -50 to 50
  } else if (trainConfig == 4150){ // timing
    cuts.AddCutMergedCalo("80083103","1111132030032200000","1111132030022700001","0163300000000000"); // -200 to 200
    cuts.AddCutMergedCalo("80083103","1111132020032200000","1111132020022700001","0163300000000000"); // -500 to 500
  } else if (trainConfig == 4151){
    cuts.AddCutMergedCalo("80083103","1111132050032200000","1111132050022700001","0163200000000000"); // rapidity < 0.7
    cuts.AddCutMergedCalo("80083103","1111132050032200000","1111132050022700001","0163100000000000"); // rapidity < 0.8


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

  Int_t numberOfCuts            = cuts.GetNCuts();
  TList *EventCutList           = new TList();
  TList *ClusterCutList         = new TList();
  TList *ClusterMergedCutList   = new TList();
  TList *MesonCutList           = new TList();

  TList *HeaderList             = new TList();
  if (doWeightingPart==1) {
    TObjString *Header1         = new TObjString("pi0_1");
    HeaderList->Add(Header1);
  }
  if (doWeightingPart==2){
    TObjString *Header3         = new TObjString("eta_2");
    HeaderList->Add(Header3);
  }
  if (doWeightingPart==3) {
    TObjString *Header1         = new TObjString("pi0_1");
    HeaderList->Add(Header1);
    TObjString *Header3         = new TObjString("eta_2");
    HeaderList->Add(Header3);
  }

  if (periodNameV0Reader.Contains("LHC18b9")||periodNameV0Reader.Contains("LHC17g8")){
    TObjString *HeaderPMB = new TObjString("EPOSLHC_0");
    TObjString *HeaderP8J = new TObjString("Pythia8Jets_1");
    if (doWeightingPart==4 || doWeightingPart==6) { // all headers
      HeaderList->Add(HeaderPMB);
      HeaderList->Add(HeaderP8J);
    } else if (doWeightingPart==5 || doWeightingPart==7) { // only MB header
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

  TString energy      = "";
  TString mcName      = "";
  TString mcNameAdd   = "";
  if (generatorName.Contains("LHC18b9")){
    energy            = "8.16TeV";
    mcName            = "LHC18b9";
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

    analysisEventCuts[i]          = new AliConvEventCuts();

    // definition of weighting input
    TString fitNamePi0      = Form("Pi0_Fit_Data_%s",energy.Data());
    TString fitNameEta      = Form("Eta_Fit_Data_%s",energy.Data());

    TString mcInputNamePi0  = "";
    TString mcInputNameEta  = "";
    mcInputNamePi0          = Form("Pi0_%s%s_%s", mcName.Data(), mcNameAdd.Data(), energy.Data() );
    mcInputNameEta          = Form("Eta_%s%s_%s", mcName.Data(), mcNameAdd.Data(), energy.Data() );

    if (doWeightingPart > 5) analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNamePtWeights, mcInputNamePi0, mcInputNameEta, "",fitNamePi0,fitNameEta);

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
  if(maxAllowedPi0Overlaps>-1){ task->SetMaxNeutralPionOverlapsMC(maxAllowedPi0Overlaps);}
  task->SetEnableDetailedM02Distribtuon(runDetailedM02);
  task->SetEventCutList(numberOfCuts,EventCutList);
  task->SetCaloCutList(numberOfCuts,ClusterCutList);
  task->SetCorrectionTaskSetting(corrTaskSetting);
  task->SetCaloMergedCutList(numberOfCuts,ClusterMergedCutList);
  task->SetMesonCutList(numberOfCuts,MesonCutList);
  task->SetDoMesonQA(enableQAMesonTask); //Attention new switch for Pi0 QA
  task->SetDoClusterQA(enableQAClusterTask);  //Attention new switch small for Cluster QA
  task->SetEnableSortingOfMCClusLabels(enableSortingMCLabels);
  if(enableExtMatchAndQA > 1){ task->SetPlotHistsExtQA(kTRUE);}

  //connect containers
  AliAnalysisDataContainer *coutput =
    mgr->CreateContainer(!(corrTaskSetting.CompareTo("")) ? Form("GammaCaloMerged_%i",trainConfig) : Form("GammaCaloMerged_%i_%s",trainConfig,corrTaskSetting.Data()), TList::Class(),
              AliAnalysisManager::kOutputContainer,Form("GammaCaloMerged_%i.root",trainConfig) );

  mgr->AddTask(task);
  mgr->ConnectInput(task,0,cinput);
  mgr->ConnectOutput(task,1,coutput);

  return;

}
