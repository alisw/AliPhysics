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
//PbPb together with all supporting classes
//***************************************************************************************

//***************************************************************************************
//main function
//***************************************************************************************
void AddTask_GammaCaloMerged_PbPb(
  Int_t     trainConfig                   = 1,        // change different set of cuts
  Int_t     isMC                          = 0,        // run MC
  TString   photonCutNumberV0Reader       = "",       // 00000008400000000100000000 nom. B, 00000088400000000100000000 low B
  TString   periodNameV0Reader            = "",
  // general setting for task
  Int_t     enableQAMesonTask             = 0,        // enable QA in AliAnalysisTaskGammaConvV1
  Int_t     enableQAClusterTask           = 0,        // enable additional QA task
  Int_t     enableExtMatchAndQA           = 0,        // disabled (0), extMatch (1), extQA_noCellQA (2), extMatch+extQA_noCellQA (3), extQA+cellQA (4), extMatch+extQA+cellQA (5)
  Int_t     enableLightOutput             = 0,        // switch to run light output (only essential histograms for afterburner)
  Int_t     enableTriggerMimicking        = 0,   // enable trigger mimicking
  Bool_t    enableTriggerOverlapRej       = kFALSE,   // enable trigger overlap rejection
  TString   settingMaxFacPtHard           = "3.",       // maximum factor between hardest jet and ptHard generated
  // settings for weights
  // FPTW:fileNamePtWeights, FCEF:fileNameCentFlattening, separate with ;
  TString   fileNameExternalInputs        = "",
  Bool_t    doWeightingPart               = kFALSE,   // enable Weighting
  Int_t     enableCentFlattening          = 0,        // enable centrality flattening
  Int_t     headerSelectionInt            = 0,        // set header selection
  TString   generatorName                 = "",
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

  TString fileNamePtWeights     = cuts.GetSpecialFileNameFromString (fileNameExternalInputs, "FPTW:");
  TString fileNameCentFlattening= cuts.GetSpecialFileNameFromString (fileNameExternalInputs, "FCEF:");

  TString addTaskName                 = "AddTask_GammaMerged_PbPb";
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
  if(rmaxFacPtHardSetting->GetEntries()<1){cout << "ERROR: AddTask_GammaCaloMerged_PbPb during parsing of settingMaxFacPtHard String '" << settingMaxFacPtHard.Data() << "'" << endl; return;}
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

  Bool_t runJetJetAndQAwithCaloPhotonCuts = kFALSE;
  if(isMC == 2 && enableExtMatchAndQA > 1) runJetJetAndQAwithCaloPhotonCuts = kTRUE;

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
  AliAnalysisTaskGammaCaloMerged *task  = NULL;
  task                                  = new AliAnalysisTaskGammaCaloMerged(Form("GammaCaloMerged_%i",trainConfig));
  task->SetIsHeavyIon(isHeavyIon);
  task->SetIsMC(isMC);
  task->SetV0ReaderName(V0ReaderName);
  task->SetLightOutput(enableLightOutput);
  task->SetTrackMatcherRunningMode(trackMatcherRunningMode);

  // cluster cuts
  // 0 "ClusterType",  1 "EtaMin", 2 "EtaMax", 3 "PhiMin", 4 "PhiMax", 5 "DistanceToBadChannel", 6 "Timing", 7 "TrackMatching", 8 "ExoticCell",
  // 9 "MinEnergy", 10 "MinNCells", 11 "MinM02", 12 "MaxM02", 13 "MinM20", 14 "MaxM20", 15 "MaximumDispersion", 16 "NLM"

  // ************************************* EMCAL cuts ****************************************************
  // LHC13b-d
  if (trainConfig == 1){ // NLM 1 no non linearity, no mass/alpha
    cuts.AddCutMergedCalo("60100013","1111100050032200000","1111100050022110001","0163300000000000"); //
    cuts.AddCutMergedCalo("61200013","1111100050032200000","1111100050022110001","0163300000000000"); //
    cuts.AddCutMergedCalo("50100013","1111100050032200000","1111100050022110001","0163300000000000"); //
    cuts.AddCutMergedCalo("51200013","1111100050032200000","1111100050022110001","0163300000000000"); //
  } else if (trainConfig == 2){ // NLM 1 open cuts
    cuts.AddCutMergedCalo("60100013","1111100050032200000","1111100050022000001","0163300000000000"); //
    cuts.AddCutMergedCalo("61200013","1111100050032200000","1111100050022000001","0163300000000000"); //
    cuts.AddCutMergedCalo("50100013","1111100050032200000","1111100050022000001","0163300000000000"); //
    cuts.AddCutMergedCalo("51200013","1111100050032200000","1111100050022000001","0163300000000000"); //
  } else if (trainConfig == 3){ // NLM 1 no non linearity, no mass/alpha
    cuts.AddCutMergedCalo("52400013","1111100050032200000","1111100050022110001","0163300000000000"); //
    cuts.AddCutMergedCalo("54600013","1111100050032200000","1111100050022110001","0163300000000000"); //
    cuts.AddCutMergedCalo("52500013","1111100050032200000","1111100050022110001","0163300000000000"); //
    cuts.AddCutMergedCalo("56800013","1111100050032200000","1111100050022110001","0163300000000000"); //
  } else if (trainConfig == 4){ // NLM 1 open cuts
    cuts.AddCutMergedCalo("52400013","1111100050032200000","1111100050022000001","0163300000000000"); //
    cuts.AddCutMergedCalo("54600013","1111100050032200000","1111100050022000001","0163300000000000"); //
    cuts.AddCutMergedCalo("52500013","1111100050032200000","1111100050022000001","0163300000000000"); //
    cuts.AddCutMergedCalo("56800013","1111100050032200000","1111100050022000001","0163300000000000"); //


  // run 2 configs
  // first look
    // CL cents
  } else if (trainConfig == 201){ // NLM 1 no non linearity, no mass/alpha
    cuts.AddCutMergedCalo("20110113","1111100050032200000","1111100050022110001","0163300000000000"); //
    cuts.AddCutMergedCalo("21210113","1111100050032200000","1111100050022110001","0163300000000000"); //
    cuts.AddCutMergedCalo("22510113","1111100050032200000","1111100050022110001","0163300000000000"); //
    cuts.AddCutMergedCalo("25910113","1111100050032200000","1111100050022110001","0163300000000000"); //
    cuts.AddCutMergedCalo("20010113","1111100050032200000","1111100050022110001","0163300000000000"); //
  } else if (trainConfig == 202){ // EMCAL clusters
    cuts.AddCutMergedCalo("20110113","1111100050032200000","1111100050022000001","0163300000000000"); //
    cuts.AddCutMergedCalo("21210113","1111100050032200000","1111100050022000001","0163300000000000"); //
    cuts.AddCutMergedCalo("22510113","1111100050032200000","1111100050022000001","0163300000000000"); //
    cuts.AddCutMergedCalo("25910113","1111100050032200000","1111100050022000001","0163300000000000"); //
    cuts.AddCutMergedCalo("20010113","1111100050032200000","1111100050022000001","0163300000000000"); //

    // V0 cents
  } else if (trainConfig == 203){ // NLM 1 no non linearity, no mass/alpha
    cuts.AddCutMergedCalo("50110113","1111100050032200000","1111100050022110001","0163300000000000"); //
    cuts.AddCutMergedCalo("51210113","1111100050032200000","1111100050022110001","0163300000000000"); //
    cuts.AddCutMergedCalo("52510113","1111100050032200000","1111100050022110001","0163300000000000"); //
    cuts.AddCutMergedCalo("55910113","1111100050032200000","1111100050022110001","0163300000000000"); //
    cuts.AddCutMergedCalo("50010113","1111100050032200000","1111100050022110001","0163300000000000"); //
  } else if (trainConfig == 204){ // NLM 1 no non linearity, no mass/alpha
    cuts.AddCutMergedCalo("50110113","1111100050032200000","1111100050022000001","0163300000000000"); //
    cuts.AddCutMergedCalo("51210113","1111100050032200000","1111100050022000001","0163300000000000"); //
    cuts.AddCutMergedCalo("52510113","1111100050032200000","1111100050022000001","0163300000000000"); //
    cuts.AddCutMergedCalo("55910113","1111100050032200000","1111100050022000001","0163300000000000"); //
    cuts.AddCutMergedCalo("50010113","1111100050032200000","1111100050022000001","0163300000000000"); //

  } else if (trainConfig == 300){ // mEDC configs 0-90%
    cuts.AddCutMergedCalo("10910a13","4117931050032200000","4117931050022700001","0163300000000000"); // INT7
  } else if (trainConfig == 301){
    cuts.AddCutMergedCalo("1098ea13","4117931050032200000","4117931050022700001","0163300000000000"); // EG2+DG2
    cuts.AddCutMergedCalo("1098da13","4117931050032200000","4117931050022700001","0163300000000000"); // EG1+DG1

  } else if (trainConfig == 320){ // mEDC configs 0-10%
    cuts.AddCutMergedCalo("10110a13","4117931050032200000","4117931050022700001","0163300000000000"); // INT7
  } else if (trainConfig == 321){
    cuts.AddCutMergedCalo("1018ea13","4117931050032200000","4117931050022700001","0163300000000000"); // EG2+DG2
    cuts.AddCutMergedCalo("1018da13","4117931050032200000","4117931050022700001","0163300000000000"); // EG1+DG1

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
    }else {
      TObjString *Header1 = new TObjString("pi0_1");
      HeaderList->Add(Header1);
      TObjString *Header2 = new TObjString("eta_2");
      HeaderList->Add(Header2);
    }
  } else if (generatorName.CompareTo("LHC14a1b")==0 || generatorName.CompareTo("LHC14a1c")==0){
    TObjString *Header1 = new TObjString("BOX");
    HeaderList->Add(Header1);
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
      mgr->AddTask(fTrackMatcher);
      mgr->ConnectInput(fTrackMatcher,0,cinput);
    }

    analysisEventCuts[i]          = new AliConvEventCuts();

    // switch on centrality flattening
    if(generatorName.CompareTo("LHC11h") && (enableCentFlattening > 0)){
      cout << "entering the flattening loop -> searching for file: " << fileNameCentFlattening.Data() << endl;
      if( fileNameCentFlattening.Contains("Low") ){
        analysisEventCuts[i]->SetUseWeightFlatCentralityFromFile(enableCentFlattening, fileNameCentFlattening, "CentLowRange");
      }else if( fileNameCentFlattening.Contains("Middle") ){
        analysisEventCuts[i]->SetUseWeightFlatCentralityFromFile(enableCentFlattening, fileNameCentFlattening, "CentMiddleRange");
      }else if( fileNameCentFlattening.Contains("High") ){
        analysisEventCuts[i]->SetUseWeightFlatCentralityFromFile(enableCentFlattening, fileNameCentFlattening, "CentHighRange");
      }else {
        analysisEventCuts[i]->SetUseWeightFlatCentralityFromFile(enableCentFlattening, fileNameCentFlattening, "Cent");
      }
    }

    analysisEventCuts[i]->SetTriggerMimicking(enableTriggerMimicking);
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
    if(periodNameV0Reader.CompareTo("") != 0) analysisEventCuts[i]->SetPeriodEnum(periodNameV0Reader);
    analysisEventCuts[i]->SetLightOutput(enableLightOutput);
    analysisEventCuts[i]->InitializeCutsFromCutString((cuts.GetEventCut(i)).Data());
    EventCutList->Add(analysisEventCuts[i]);
    analysisEventCuts[i]->SetFillCutHistograms("",kFALSE);

    analysisClusterCuts[i]        = new AliCaloPhotonCuts(runJetJetAndQAwithCaloPhotonCuts);
    analysisClusterCuts[i]->SetIsPureCaloCut(2);
    analysisClusterCuts[i]->SetHistoToModifyAcceptance(histoAcc);
    analysisClusterCuts[i]->SetV0ReaderName(V0ReaderName);
    analysisClusterCuts[i]->SetCaloTrackMatcherName(TrackMatcherName);
    analysisClusterCuts[i]->SetLightOutput(enableLightOutput);
    analysisClusterCuts[i]->InitializeCutsFromCutString((cuts.GetClusterCut(i)).Data());
    ClusterCutList->Add(analysisClusterCuts[i]);
    analysisClusterCuts[i]->SetExtendedMatchAndQA(enableExtMatchAndQA);
    analysisClusterCuts[i]->SetFillCutHistograms("");

    analysisClusterMergedCuts[i]  = new AliCaloPhotonCuts(runJetJetAndQAwithCaloPhotonCuts);
    analysisClusterMergedCuts[i]->SetIsPureCaloCut(1);
    analysisClusterMergedCuts[i]->SetHistoToModifyAcceptance(histoAcc);
    analysisClusterMergedCuts[i]->SetV0ReaderName(V0ReaderName);
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
  task->SetSelectedMesonID(selectedMeson);
  task->SetEventCutList(numberOfCuts,EventCutList);
  task->SetCaloCutList(numberOfCuts,ClusterCutList);
  task->SetCaloMergedCutList(numberOfCuts,ClusterMergedCutList);
  task->SetMesonCutList(numberOfCuts,MesonCutList);
  task->SetDoMesonQA(enableQAMesonTask); //Attention new switch for Pi0 QA
  task->SetDoClusterQA(enableQAClusterTask);  //Attention new switch small for Cluster QA
  task->SetEnableSortingOfMCClusLabels(enableSortingMCLabels);
  if(enableExtMatchAndQA > 1){ task->SetPlotHistsExtQA(kTRUE);}
  if (enableDetailedPrintout) task->SetEnableDetailedPrintout(enableDetailedPrintout);//Attention new switch small for Cluster QA

  //connect containers
  AliAnalysisDataContainer *coutput =
    mgr->CreateContainer(Form("GammaCaloMerged_%i",trainConfig), TList::Class(),
              AliAnalysisManager::kOutputContainer,Form("GammaCaloMerged_%i.root",trainConfig));

  mgr->AddTask(task);
  mgr->ConnectInput(task,0,cinput);
  mgr->ConnectOutput(task,1,coutput);

  return;

}
