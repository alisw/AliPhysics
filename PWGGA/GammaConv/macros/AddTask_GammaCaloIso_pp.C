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
//($ALIPHYSICS/PWGGA/GammaConv/AliAnalysisTaskGammaCaloIso.cxx) for
//pp together with all supporting classes
//***************************************************************************************

//***************************************************************************************
//main function
//***************************************************************************************
void AddTask_GammaCaloIso_pp(
  Int_t     trainConfig                   = 1,        // change different set of cuts
  Int_t     isMC                          = 0,        // run MC
  TString   photonCutNumberV0Reader       = "",       // 00000008400000000100000000 nom. B, 00000088400000000100000000 low B
  TString   periodNameV0Reader            = "",
  // general setting for task
  Int_t     enableQAMesonTask             = 0,        // enable QA in AliAnalysisTaskGammaConvV1
  Int_t     enableQAClusterTask           = 0,        // enable additional QA task
  Int_t     enableExtMatchAndQA           = 0,                            // disabled (0), extMatch (1), extQA_noCellQA (2), extMatch+extQA_noCellQA (3), extQA+cellQA (4), extMatch+extQA+cellQA (5)
  Bool_t    enableLightOutput             = kFALSE,   // switch to run light output (only essential histograms for afterburner)
  Bool_t    enableTHnSparse               = kFALSE,   // switch on THNsparse
  Int_t     enableTriggerMimicking        = 0,        // enable trigger mimicking
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


  TString addTaskName                 = "AddTask_GammaCaloIso_pp";
  TString sAdditionalTrainConfig      = cuts.GetSpecialSettingFromAddConfig(additionalTrainConfig, "", "", addTaskName);
  if (sAdditionalTrainConfig.Atoi() > 0){
    trainConfig = trainConfig + sAdditionalTrainConfig.Atoi();
    cout << "INFO: " << addTaskName.Data() << " running additionalTrainConfig '" << sAdditionalTrainConfig.Atoi() << "', train config: '" << trainConfig << "'" << endl;
  }

  TString fileNamePtWeights           = cuts.GetSpecialFileNameFromString (fileNameExternalInputs, "FPTW:");
  TString fileNameMultWeights         = cuts.GetSpecialFileNameFromString (fileNameExternalInputs, "FMUW:");

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
  if(rmaxFacPtHardSetting->GetEntries()<1){cout << "ERROR: AddTask_GammaCaloIso_pp during parsing of settingMaxFacPtHard String '" << settingMaxFacPtHard.Data() << "'" << endl; return;}
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
  AliAnalysisTaskGammaCaloIso *task=NULL;
  task= new AliAnalysisTaskGammaCaloIso(Form("GammaCaloIso_%i",trainConfig));
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
  // *************************** pp 13 TeV cuts paper EMC ************************************************
  // *****************************************************************************************************
  //*************************** STANDARD EMC EMC MESON CUTS **********************************************
  //******************************************************************************************************

  //*********************************NO NONLINEARITY CUTS*************************************************

  //EG2Trigger

//adopted to Isolated Photon Analysis 13TeV
  if (trainConfig == 22){ //Isolation Cut variations, Trackmatching, maxM02
  cuts.AddCutCalo("00085113","e2211000600322l0000","0163103100000050"); //Standard Meson EmCal Cuts
  cuts.AddCutCalo("00085113","e221100360512260000","0163103100000050"); //Isolation Analysis Cuts
  cuts.AddCutCalo("00085113","e221100060512260000","0163103100000050"); //test distance to bad channel
  cuts.AddCutCalo("00085113","e221100360012260000","0163103100000050"); //test exotic cluster rejection
  cuts.AddCutCalo("00085113","e2211003605122l0000","0163103100000050"); //test shape cut

  //standart Meson EmCal Cuts + Isolation
  } else if (trainConfig == 23){ // with Time Cut
  cuts.AddCutCalo("00085113","12211000600322l0000","0163103100000050");
  cuts.AddCutCalo("00085113","b2211000600322l0000","0163103100000050"); // R=0.1
  cuts.AddCutCalo("00085113","c2211000600322l0000","0163103100000050"); // R=0.2
  cuts.AddCutCalo("00085113","d2211000600322l0000","0163103100000050"); // R=0.3
  cuts.AddCutCalo("00085113","e2211000600322l0000","0163103100000050"); // R=0.4

  // *****************************************************************************************************
  // ******************************* pp 13 TeV DATA 16k,l ************************************************

  //standart Meson EmCal Cuts + Isolation
  } else if (trainConfig == 25){ // without Timing Cut for V2Clusterizer
  cuts.AddCutCalo("00085113","12211000000322l0000","0163103100000050");
  cuts.AddCutCalo("00085113","b2211000000322l0000","0163103100000050"); // R=0.1
  cuts.AddCutCalo("00085113","c2211000000322l0000","0163103100000050"); // R=0.2
  cuts.AddCutCalo("00085113","d2211000000322l0000","0163103100000050"); // R=0.3
  cuts.AddCutCalo("00085113","e2211000000322l0000","0163103100000050"); // R=0.4

  //adopted to Isolated Photon Analysis 13TeV
  } else if (trainConfig == 26){ // without TrackMatching
  cuts.AddCutCalo("00085113","e2211000000322l0000","0163103100000050"); //Standard Meson EmCal Cuts
  cuts.AddCutCalo("00085113","e221100300512260000","0163103100000050"); //Isolation Analysis Cuts
  cuts.AddCutCalo("00085113","e221100000512260000","0163103100000050"); //test distance to bad channel
  cuts.AddCutCalo("00085113","e221100300012260000","0163103100000050"); //test exotic cluster rejection
  cuts.AddCutCalo("00085113","e2211003005122l0000","0163103100000050"); //test shape cut

  //EG1 Trigger

  //adopted to Isolated Photon Analysis 13TeV
   } else if (trainConfig == 32){ //Isolation Cut variations, Trackmatching, maxM02
    cuts.AddCutCalo("00083113","e2211000600322l0000","0163103100000050"); //Standard Meson EmCal Cuts
    cuts.AddCutCalo("00083113","e221100360512260000","0163103100000050"); //Isolation Analysis Cuts
    cuts.AddCutCalo("00083113","e221100060512260000","0163103100000050"); //test distance to bad channel
    cuts.AddCutCalo("00083113","e221100360012260000","0163103100000050"); //test exotic cluster rejection
    cuts.AddCutCalo("00083113","e2211003605122l0000","0163103100000050"); //test shape cut

    //standart Meson EmCal Cuts + Isolation
    } else if (trainConfig == 33){ // with Time Cut
    cuts.AddCutCalo("00083113","12211000600322l0000","0163103100000050");
    cuts.AddCutCalo("00083113","b2211000600322l0000","0163103100000050"); // R=0.1
    cuts.AddCutCalo("00083113","c2211000600322l0000","0163103100000050"); // R=0.2
    cuts.AddCutCalo("00083113","d2211000600322l0000","0163103100000050"); // R=0.3
    cuts.AddCutCalo("00083113","e2211000600322l0000","0163103100000050"); // R=0.4

    // *****************************************************************************************************
    // ******************************* pp 13 TeV DATA 16k,l ************************************************

    //standart Meson EmCal Cuts + Isolation
    } else if (trainConfig == 35){ // without Timing Cut for V2Clusterizer
    cuts.AddCutCalo("00083113","12211000a00322l0000","0163103100000050");
    cuts.AddCutCalo("00083113","b2211000a00322l0000","0163103100000050"); // R=0.1
    cuts.AddCutCalo("00083113","c2211000a00322l0000","0163103100000050"); // R=0.2
    cuts.AddCutCalo("00083113","d2211000a00322l0000","0163103100000050"); // R=0.3
    cuts.AddCutCalo("00083113","e2211000a00322l0000","0163103100000050"); // R=0.4

    //adopted to Isolated Photon Analysis 13TeV
    } else if (trainConfig == 36){ // without Timing Cut
    cuts.AddCutCalo("00083113","e2211000100322l0000","0163103100000050"); //Standard Meson EmCal Cuts
    cuts.AddCutCalo("00083113","e221100310512260000","0163103100000050"); //Isolation Analysis Cuts
    cuts.AddCutCalo("00083113","e221100010512260000","0163103100000050"); //test distance to bad channel
    cuts.AddCutCalo("00083113","e221100310012260000","0163103100000050"); //test exotic cluster rejection
    cuts.AddCutCalo("00083113","e2211003105122l0000","0163103100000050"); //test shape cut

  } else {
    Error(Form("GammaCaloIso_%i",trainConfig), "wrong trainConfig variable no cuts have been specified for the configuration");
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
  if (generatorName.Contains("LHC17i3a1")){
    energy = "13TeV";
    mcName = "Pythia8_LHC17i3a1";
  } else if (generatorName.Contains("LHC17i3c1")){
      energy = "13TeV";
      mcName = "Pythia8_LHC17i3c1";
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
    if( !(AliCaloTrackMatcher*)mgr->GetTask(TrackMatcherName.Data()) ){
      AliCaloTrackMatcher* fTrackMatcher = new AliCaloTrackMatcher(TrackMatcherName.Data(),caloCutPos.Atoi(),trackMatcherRunningMode);
      fTrackMatcher->SetV0ReaderName(V0ReaderName);
      fTrackMatcher->SetCorrectionTaskSetting(corrTaskSetting);
      mgr->AddTask(fTrackMatcher);
      mgr->ConnectInput(fTrackMatcher,0,cinput);
    }

    //create AliPhotonIsolation instance, if there is none present
    TString PhotonIsolationName = "PhotonIsolation";
    if( !(AliPhotonIsolation*)mgr->GetTask(PhotonIsolationName.Data()) ){
      AliPhotonIsolation* fPhotonIsolation = new AliPhotonIsolation(PhotonIsolationName.Data());
     // fPhotonIsolation->SetV0ReaderName(V0ReaderName);
     // fPhotonIsolation->SetCorrectionTaskSetting(corrTaskSetting);
      mgr->AddTask(fPhotonIsolation);
      mgr->ConnectInput(fPhotonIsolation,0,cinput);
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
    analysisEventCuts[i]->SetTriggerOverlapRejecion(enableTriggerOverlapRej);
    if(fMinPtHardSet)
      analysisEventCuts[i]->SetMinFacPtHard(minFacPtHard);
    if(fMaxPtHardSet)
      analysisEventCuts[i]->SetMaxFacPtHard(maxFacPtHard);
    if(fSingleMaxPtHardSet)
      analysisEventCuts[i]->SetMaxFacPtHardSingleParticle(maxFacPtHardSingle);    analysisEventCuts[i]->SetV0ReaderName(V0ReaderName);
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
  task->SetDoMesonAnalysis(kTRUE);
  task->SetCorrectionTaskSetting(corrTaskSetting);
  task->SetDoMesonQA(enableQAMesonTask); //Attention new switch for Pi0 QA
  task->SetDoClusterQA(enableQAClusterTask);  //Attention new switch small for Cluster QA
  task->SetDoTHnSparse(enableTHnSparse);
  task->SetProduceTreeEOverP(doTreeEOverP);
  task->SetEnableSortingOfMCClusLabels(enableSortingMCLabels);
  if(trainConfig == 446) task->SetSoftAnalysis(kTRUE);
  if(enableExtMatchAndQA > 1){ task->SetPlotHistsExtQA(kTRUE);}
  if(trainConfig == 106 || trainConfig == 125 || trainConfig == 145){
    task->SetInOutTimingCluster(-30e-9,35e-9);
  }
  task->SetLocalDebugFlag(localDebugFlag);

  //connect containers
  AliAnalysisDataContainer *coutput =
    mgr->CreateContainer(!(corrTaskSetting.CompareTo("")) ? Form("GammaCaloIso_%i",trainConfig) : Form("GammaCaloIso_%i_%s",trainConfig,corrTaskSetting.Data()), TList::Class(),
              AliAnalysisManager::kOutputContainer,Form("GammaCaloIso_%i.root",trainConfig));

  mgr->AddTask(task);
  mgr->ConnectInput(task,0,cinput);
  mgr->ConnectOutput(task,1,coutput);

  return;

}
