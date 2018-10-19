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
  TString   periodname                    = "",        // period name
  // general setting for task
  Int_t     enableQAMesonTask             = 0,        // enable QA in AliAnalysisTaskGammaConvV1
  Int_t     enableQAClusterTask           = 0,        // enable additional QA task
  Int_t     enableExtMatchAndQA           = 0,                            // disabled (0), extMatch (1), extQA_noCellQA (2), extMatch+extQA_noCellQA (3), extQA+cellQA (4), extMatch+extQA+cellQA (5)
  Int_t     enableLightOutput             = 0,   // switch to run light output (only essential histograms for afterburner)
  Bool_t    enableTriggerMimicking        = kFALSE,   // enable trigger mimicking
  Bool_t    enableTriggerOverlapRej       = kFALSE,   // enable trigger overlap rejection
  Float_t   maxFacPtHard                  = 3.,       // maximum factor between hardest jet and ptHard generated
  // settings for weights
  // FPTW:fileNamePtWeights, separate with ;
  TString   fileNameExternalInputs        = "",
  Int_t    doWeightingPart                = 0,        // enable Weighting
  // special settings
  Int_t     selectedMeson                 = 1,                  // put flag for selected meson
  Bool_t    enableSortingMCLabels         = kTRUE,    // enable sorting for MC cluster labels
  Bool_t    enableDetailedPrintout        = kFALSE,             // enable detailed printout
  Double_t  minEnergyForExoticsCut        = 1.0,                // minimum energy to be used for exotics CutHandler
  Bool_t    enableExoticsQA               = kFALSE,             // switch to run QA for exotic clusters
  Bool_t    runDetailedM02                = kFALSE,             // switch on very detailed M02 distribution
  // subwagon config
  TString   additionalTrainConfig         = "0"       // additional counter for trainconfig
  ) {

  AliCutHandlerPCM cuts;

  TString fileNamePtWeights     = cuts.GetSpecialFileNameFromString (fileNameExternalInputs, "FPTW:");

  TString sAdditionalTrainConfig      = cuts.GetSpecialSettingFromAddConfig(additionalTrainConfig, "", "");
  if (sAdditionalTrainConfig.Atoi() > 0){
    trainConfig = trainConfig + sAdditionalTrainConfig.Atoi();
    cout << "INFO: AddTask_GammaConvV1_pPb running additionalTrainConfig '" << sAdditionalTrainConfig.Atoi() << "', train config: '" << trainConfig << "'" << endl;
  }
  TString corrTaskSetting             = cuts.GetSpecialSettingFromAddConfig(additionalTrainConfig, "CF", "");
  if(corrTaskSetting.CompareTo(""))
    cout << "corrTaskSetting: " << corrTaskSetting.Data() << endl;

  Int_t trackMatcherRunningMode = 0; // CaloTrackMatcher running mode
  TString strTrackMatcherRunningMode  = cuts.GetSpecialSettingFromAddConfig(additionalTrainConfig, "TM", "");
  if(additionalTrainConfig.Contains("TM"))
    trackMatcherRunningMode = strTrackMatcherRunningMode.Atoi();

  TH1S* histoAcc = 0x0;         // histo for modified acceptance
  TString strModifiedAcc              = cuts.GetSpecialSettingFromAddConfig(additionalTrainConfig, "MODIFYACC", "");
  if(strModifiedAcc.Contains("MODIFYACC")){
    cout << "INFO: AddTask_GammaCalo_pp activating 'MODIFYacc'" << endl;
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

  Int_t isHeavyIon = 2;

  // ================== GetAnalysisManager ===============================
  AliAnalysisManager *mgr           = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error(Form("AddTask_GammaCaloMerged_pPb_%i",trainConfig), "No analysis manager found.");
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
    if (doWeightingPart==4) {
      HeaderList->Add(HeaderPMB);
      HeaderList->Add(HeaderP8J);
    } else {
      HeaderList->Add(HeaderP8J);
    }
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
    if( !(AliCaloTrackMatcher*)mgr->GetTask(TrackMatcherName.Data()) ){
      AliCaloTrackMatcher* fTrackMatcher = new AliCaloTrackMatcher(TrackMatcherName.Data(),caloCutPos.Atoi(),trackMatcherRunningMode);
      fTrackMatcher->SetV0ReaderName(V0ReaderName);
      fTrackMatcher->SetCorrectionTaskSetting(corrTaskSetting);
      mgr->AddTask(fTrackMatcher);
      mgr->ConnectInput(fTrackMatcher,0,cinput);
    }

    analysisEventCuts[i]          = new AliConvEventCuts();

    analysisEventCuts[i]->SetTriggerMimicking(enableTriggerMimicking);
    analysisEventCuts[i]->SetTriggerOverlapRejecion(enableTriggerOverlapRej);
    analysisEventCuts[i]->SetMaxFacPtHard(maxFacPtHard);
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
