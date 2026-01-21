/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: Joshua Koenig                                                  *
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
// This AddTask is supposed to set up the main task
//($ALIPHYSICS/PWGGA/GammaConv/AliAnalysisTaskMesonJetCorrelation.cxx) for
// pp together with all supporting classes
//***************************************************************************************

//***************************************************************************************
// main function
//***************************************************************************************

void AddTask_MesonJetCorr_K0Lambda(
  Int_t trainConfig = 1,                // change different set of cuts
  Int_t isMC = 0,                       // run MC
  int meson = 0,                        // meson: 2 = K0, 3 = Lambda, 4 = anti-Lambda
  TString photonCutNumberV0Reader = "", // 00000008400000000100000000 nom. B, 00000088400000000100000000 low B
  TString periodNameV0Reader = "",
  // general setting for task
  Int_t enableQAK0Lambda = 0,             // enable QA in AliAnalysisTaskGammaConvV1
  Int_t enableQAJets = 0,                  // enable additional QA for jets
  int enableLightOutput = kFALSE,          // switch to run light output (only essential histograms for afterburner)
  Bool_t enableTHnSparse = kFALSE,         // switch on THNsparse
  Int_t enableTriggerMimicking = 0,        // enable trigger mimicking
  Bool_t enableTriggerOverlapRej = kFALSE, // enable trigger overlap rejection
  TString settingMaxFacPtHard = "3.",      // maximum factor between hardest jet and ptHard generated
  Int_t debugLevel = 0,                    // introducing debug levels for grid running
  // settings for weights
  // FPTW:fileNamePtWeights, FMUW:fileNameMultWeights, FMAW:fileNameMatBudWeights, FEPC:fileNamedEdxPostCalib, separate with ;
  TString fileNameExternalInputs = "",
  Int_t doWeightingPart = 0,                   // enable Weighting
  TString generatorName = "DPMJET",            // generator Name
  Bool_t enableMultiplicityWeighting = kFALSE, //
  TString periodNameAnchor = "",               //
  Int_t enableMatBudWeightsPi0 = 0,            // 1 = three radial bins, 2 = 10 radial bins
  Bool_t enableElecDeDxPostCalibration = kFALSE,
  // special settings
  Bool_t enableChargedPrimary = kFALSE,
  bool doFillMesonDCATree = false, // switch to enable filling the meson DCA tree for pile-up estimation
  bool useCentralEvtSelection = true,
  bool setPi0Unstable = false,
  bool enableAddBackground = false,
  bool enableRadiusDep = false,
  int runOnlyZPt = 0,           // if 0, bot pt and z histograms will be filled, if 1, pt histograms will be filled, if 2, only z histograms will be filled
  bool doTrackingStudies = false,
  bool cutOnJetEnergyAsymm = false,
  // subwagon config
  TString additionalTrainConfig = "0" // additional counter for trainconfig + special settings
)
{

  AliCutHandlerPCM cuts(13);

  TString addTaskName = "AddTask_MesonJetCorr_ConvCalo";
  TString sAdditionalTrainConfig = cuts.GetSpecialSettingFromAddConfig(additionalTrainConfig, "", "", addTaskName);
  if (sAdditionalTrainConfig.Atoi() > 0) {
    trainConfig = trainConfig + sAdditionalTrainConfig.Atoi();
    cout << "INFO: " << addTaskName.Data() << " running additionalTrainConfig '" << sAdditionalTrainConfig.Atoi() << "', train config: '" << trainConfig << "'" << endl;
  }

  TString fileNamePtWeights = cuts.GetSpecialFileNameFromString(fileNameExternalInputs, "FPTW:");
  TString fileNameMultWeights = cuts.GetSpecialFileNameFromString(fileNameExternalInputs, "FMUW:");
  TString fileNamedEdxPostCalib = cuts.GetSpecialFileNameFromString(fileNameExternalInputs, "FEPC:");
  TString fileNameCustomTriggerMimicOADB = cuts.GetSpecialFileNameFromString(fileNameExternalInputs, "FTRM:");
  TString fileNameMatBudWeights = cuts.GetSpecialFileNameFromString(fileNameExternalInputs, "FMAW:");
  TString fileNameArmPodExpectation = cuts.GetSpecialFileNameFromString(fileNameExternalInputs, "FAPE:");

  Int_t trackMatcherRunningMode = 0; // CaloTrackMatcher running mode
  TString strTrackMatcherRunningMode = cuts.GetSpecialSettingFromAddConfig(additionalTrainConfig, "TM", "", addTaskName);
  if (additionalTrainConfig.Contains("TM"))
    trackMatcherRunningMode = strTrackMatcherRunningMode.Atoi();

  TString nameJetFinder = (additionalTrainConfig.Contains("JET") == true) ? cuts.GetSpecialSettingFromAddConfig(additionalTrainConfig, "JET", "", addTaskName) : "";
  printf("nameJetFinder: %s\n", nameJetFinder.Data());
  
  TObjArray* rmaxFacPtHardSetting = settingMaxFacPtHard.Tokenize("_");
  if (rmaxFacPtHardSetting->GetEntries() < 1) {
    cout << "ERROR: AddTask_MesonJetCorr_pp during parsing of settingMaxFacPtHard String '" << settingMaxFacPtHard.Data() << "'" << endl;
    return;
  }
  Bool_t fMinPtHardSet = kFALSE;
  Double_t minFacPtHard = -1;
  Bool_t fMaxPtHardSet = kFALSE;
  Double_t maxFacPtHard = 100;
  Bool_t fSingleMaxPtHardSet = kFALSE;
  Double_t maxFacPtHardSingle = 100;
  Bool_t fJetFinderUsage = kFALSE;
  Bool_t fUsePtHardFromFile = kFALSE;
  Bool_t fUseAddOutlierRej = kFALSE;
  for (Int_t i = 0; i < rmaxFacPtHardSetting->GetEntries(); i++) {
    TObjString* tempObjStrPtHardSetting = (TObjString*)rmaxFacPtHardSetting->At(i);
    TString strTempSetting = tempObjStrPtHardSetting->GetString();
    if (strTempSetting.BeginsWith("MINPTHFAC:")) {
      strTempSetting.Replace(0, 10, "");
      minFacPtHard = strTempSetting.Atof();
      cout << "running with min pT hard jet fraction of: " << minFacPtHard << endl;
      fMinPtHardSet = kTRUE;
    } else if (strTempSetting.BeginsWith("MAXPTHFAC:")) {
      strTempSetting.Replace(0, 10, "");
      maxFacPtHard = strTempSetting.Atof();
      cout << "running with max pT hard jet fraction of: " << maxFacPtHard << endl;
      fMaxPtHardSet = kTRUE;
    } else if (strTempSetting.BeginsWith("MAXPTHFACSINGLE:")) {
      strTempSetting.Replace(0, 16, "");
      maxFacPtHardSingle = strTempSetting.Atof();
      cout << "running with max single particle pT hard fraction of: " << maxFacPtHardSingle << endl;
      fSingleMaxPtHardSet = kTRUE;
    } else if (strTempSetting.BeginsWith("USEJETFINDER:")) {
      strTempSetting.Replace(0, 13, "");
      if (strTempSetting.Atoi() == 1) {
        cout << "using MC jet finder for outlier removal" << endl;
        fJetFinderUsage = kTRUE;
      }
    } else if (strTempSetting.BeginsWith("PTHFROMFILE:")) {
      strTempSetting.Replace(0, 12, "");
      if (strTempSetting.Atoi() == 1) {
        cout << "using MC jet finder for outlier removal" << endl;
        fUsePtHardFromFile = kTRUE;
      }
    } else if (strTempSetting.BeginsWith("ADDOUTLIERREJ:")) {
      strTempSetting.Replace(0, 14, "");
      if (strTempSetting.Atoi() == 1) {
        cout << "using path based outlier removal" << endl;
        fUseAddOutlierRej = kTRUE;
      }
    } else if (rmaxFacPtHardSetting->GetEntries() == 1 && strTempSetting.Atof() > 0) {
      maxFacPtHard = strTempSetting.Atof();
      cout << "running with max pT hard jet fraction of: " << maxFacPtHard << endl;
      fMaxPtHardSet = kTRUE;
    }
  }

  Int_t isHeavyIon = 0;

  // ================== GetAnalysisManager ===============================
  AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error(Form("%s_%i", addTaskName.Data(), trainConfig), "No analysis manager found.");
    return;
  }

  // ================== GetInputEventHandler =============================
  AliVEventHandler* inputHandler = mgr->GetInputEventHandler();

  //=========  Set Cutnumber for V0Reader ================================
  TString cutnumberPhoton = photonCutNumberV0Reader.Data();
  TString cutnumberEvent = "00000003";
  AliAnalysisDataContainer* cinput = mgr->GetCommonInputContainer();

  //========= Add V0 Reader to  ANALYSIS manager if not yet existent =====
  TString V0ReaderName = Form("V0ReaderV1_%s_%s", cutnumberEvent.Data(), cutnumberPhoton.Data());
  AliV0ReaderV1* fV0ReaderV1 = NULL;
  if (!(AliV0ReaderV1*)mgr->GetTask(V0ReaderName.Data())) {
    cout << "V0Reader: " << V0ReaderName.Data() << " not found!!" << endl;
    return;
  } else {
    cout << "V0Reader: " << V0ReaderName.Data() << " found!!" << endl;
  }

  //================================================
  //========= Add task to the ANALYSIS manager =====
  //================================================
  AliAnalysisTaskMesonJetCorrelation* task = NULL;
  task = new AliAnalysisTaskMesonJetCorrelation(Form("MesonJetCorrelation_%i_%i", meson, trainConfig));
  task->SetIsHeavyIon(isHeavyIon);
  task->SetIsMC(isMC);
  task->SetV0ReaderName(V0ReaderName);
  task->SetTrackMatcherRunningMode(trackMatcherRunningMode); // have to do this!


  if (trainConfig == 1) { // default cuts but no armenteros
    cuts.AddCutPCM("00010103", "10200011500000000000000000", "2r52103l00000000");
  } else if (trainConfig == 2) { // default cuts but no armenteros
    cuts.AddCutPCM("000fc103", "10200011500000000000000000", "2r52103l00000000");
  } else if (trainConfig == 3) { // default cuts but no armenteros
    cuts.AddCutPCM("000fb103", "10200011500000000000000000", "2r52103l00000000"); 

  } else if (trainConfig == 4) { // default cuts, loose armenteros
    cuts.AddCutPCM("00010103", "10200011500000001000000000", "2r52103l00000000");
  } else if (trainConfig == 5) { // default cuts, loose armenteros
    cuts.AddCutPCM("000fc103", "10200011500000001000000000", "2r52103l00000000");
  } else if (trainConfig == 6) { // default cuts, loose armenteros
    cuts.AddCutPCM("000fb103", "10200011500000001000000000", "2r52103l00000000"); 

  } else if (trainConfig == 7) { // default cuts but no armenteros, offline V0 reader
    cuts.AddCutPCM("00010103", "00200011500000000000000000", "2r52103l00000000");

  } else if (trainConfig == 10) { // default cuts with cut around armenteros of expected position within 0.02
    cuts.AddCutPCM("00010103", "1020001150000000b000000000", "2r52103l00000000");
  } else if (trainConfig == 11) { // default cuts with cut around armenteros of expected position within 0.02
    cuts.AddCutPCM("000fc103", "1020001150000000b000000000", "2r52103l00000000");
  } else if (trainConfig == 12) { // default cuts with cut around armenteros of expected position within 0.02
    cuts.AddCutPCM("000fb103", "1020001150000000b000000000", "2r52103l00000000"); 

  //---------------------------------------
  // configs for K0s meson pp 13 TeV
  //---------------------------------------
  } else if (trainConfig == 101) { // K0s default cuts, min bias
    cuts.AddCutPCM("00010103", "10200011500000001000000000", "2r52103l00000000"); 
  } else if (trainConfig == 102) { // K0s default cuts, EJ2
    cuts.AddCutPCM("000fc103", "10200011500000001000000000", "2r52103l00000000"); 
  } else if (trainConfig == 103) { // K0s default cuts, EJ1
    cuts.AddCutPCM("000fb103", "10200011500000001000000000", "2r52103l00000000"); 

  //---------------------------------------
  // configs for Lambda pp 13 TeV
  //---------------------------------------
  } else if (trainConfig == 201) { // K0s default cuts, min bias
    cuts.AddCutPCM("00010103", "10200011500000002000000000", "2r52103l00000000"); 
  } else if (trainConfig == 202) { // K0s default cuts, EJ2
    cuts.AddCutPCM("000fc103", "10200011500000002000000000", "2r52103l00000000"); 
  } else if (trainConfig == 203) { // K0s default cuts, EJ1
    cuts.AddCutPCM("000fb103", "10200011500000002000000000", "2r52103l00000000"); 

  //---------------------------------------
  // configs for Anti-Lambda pp 13 TeV
  //---------------------------------------
  } else if (trainConfig == 301) { // K0s default cuts, min bias
    cuts.AddCutPCM("00010103", "10200011500000003000000000", "2r52103l00000000"); 
  } else if (trainConfig == 302) { // K0s default cuts, EJ2
    cuts.AddCutPCM("000fc103", "10200011500000003000000000", "2r52103l00000000"); 
  } else if (trainConfig == 303) { // K0s default cuts, EJ1
    cuts.AddCutPCM("000fb103", "10200011500000003000000000", "2r52103l00000000"); 

  } else {
    Error(Form("MesonJetCorrelation_%i", trainConfig), "wrong trainConfig variable no cuts have been specified for the configuration");
    return;
  }

  if (!cuts.AreValid()) {
    cout << "\n\n****************************************************" << endl;
    cout << "ERROR: No valid cuts stored in CutHandlerCalo! Returning..." << endl;
    cout << "****************************************************\n\n"
         << endl;
    return;
  }

  Int_t numberOfCuts = cuts.GetNCuts();

  TList* EventCutList = new TList();
  TList* K0LambdaCutList = new TList();
  TList* MesonCutList = new TList();

  EventCutList->SetOwner(kTRUE);
  AliConvEventCuts** analysisEventCuts = new AliConvEventCuts*[numberOfCuts];
  K0LambdaCutList->SetOwner(kTRUE);
  AliConvK0LambdaCuts** analysisK0LambdaCuts = new AliConvK0LambdaCuts*[numberOfCuts];
  MesonCutList->SetOwner(kTRUE);
  AliConversionMesonCuts** analysisMesonCuts = new AliConversionMesonCuts*[numberOfCuts];

  for (Int_t i = 0; i < numberOfCuts; i++) {

    //---------------------------------------------------------//
    //------------------------ Event Cuts ---------------------//
    //---------------------------------------------------------//
    analysisEventCuts[i] = new AliConvEventCuts();
    // analysisEventCuts[i]->SetCaloTriggerHelperName(TriggerHelperName.Data());
    analysisEventCuts[i]->SetTriggerMimicking(enableTriggerMimicking);
    analysisEventCuts[i]->SetTriggerOverlapRejecion(enableTriggerOverlapRej);
    analysisEventCuts[i]->SetV0ReaderName(V0ReaderName);
    if (enableLightOutput > 0)
      analysisEventCuts[i]->SetLightOutput(kTRUE);
    if (fMinPtHardSet)
      analysisEventCuts[i]->SetMinFacPtHard(minFacPtHard);
    if (fMaxPtHardSet)
      analysisEventCuts[i]->SetMaxFacPtHard(maxFacPtHard);
    if (fSingleMaxPtHardSet)
      analysisEventCuts[i]->SetMaxFacPtHardSingleParticle(maxFacPtHardSingle);
    if (fJetFinderUsage)
      analysisEventCuts[i]->SetUseJetFinderForOutliers(kTRUE);
    if (fUsePtHardFromFile)
      analysisEventCuts[i]->SetUsePtHardBinFromFile(kTRUE);
    if (fUseAddOutlierRej)
      analysisEventCuts[i]->SetUseAdditionalOutlierRejection(kTRUE);
    if (periodNameV0Reader.CompareTo("") != 0)
      analysisEventCuts[i]->SetPeriodEnum(periodNameV0Reader);
    analysisEventCuts[i]->InitializeCutsFromCutString((cuts.GetEventCut(i)).Data());
    analysisEventCuts[i]->SetFillCutHistograms("", kFALSE);
    EventCutList->Add(analysisEventCuts[i]);

    //---------------------------------------------------------//
    //------------------- K0+Lambdan Cuts --------------------//
    //---------------------------------------------------------//
    analysisK0LambdaCuts[i] = new AliConvK0LambdaCuts();
    if (enableLightOutput > 0)
      analysisK0LambdaCuts[i]->SetLightOutput(1);
    analysisK0LambdaCuts[i]->SetDoQA(enableQAK0Lambda); 
    analysisK0LambdaCuts[i]->SetIsMC(isMC); 
    if(!fileNameArmPodExpectation.EqualTo("")){
      analysisK0LambdaCuts[i]->InitArmPodRefHistos(fileNameArmPodExpectation);
    }
    analysisK0LambdaCuts[i]->InitializeCutsFromCutString((cuts.GetPhotonCut(i)).Data());
    analysisK0LambdaCuts[i]->SetFillCutHistograms("K0LambdaCuts", true);
    K0LambdaCutList->Add(analysisK0LambdaCuts[i]);

    //---------------------------------------------------------//
    //------------------------ Meson Cuts ---------------------//
    //---------------------------------------------------------//
    analysisMesonCuts[i] = new AliConversionMesonCuts();
    analysisMesonCuts[i]->InitializeCutsFromCutString((cuts.GetMesonCut(i)).Data());
    analysisMesonCuts[i]->SetFillCutHistograms("");
    // analysisEventCuts[i]->SetAcceptedHeader(HeaderList);
    analysisMesonCuts[i]->SetLightOutput(kTRUE);
    MesonCutList->Add(analysisMesonCuts[i]);

  }

  task->SetMesonKind(meson);
  task->SetMesonZPt(runOnlyZPt);
  if(meson < 2) task->SetIsConv(true);
  if (enableQAK0Lambda) task->SetDoMesonQA(true);
  task->SetJetContainerAddName(nameJetFinder);
  task->SetEventCutList(numberOfCuts, EventCutList);
  task->SetK0LambdaCutList(numberOfCuts, K0LambdaCutList);
  task->SetMesonCutList(numberOfCuts, MesonCutList);
  task->SetDoJetQA(enableQAJets);
  task->SetUseTHnSparseForResponse(enableTHnSparse);
  if(enableMatBudWeightsPi0) task->SetDoMaterialBudgetWeightingOfGammasForTrueMesons(true);
  if(doFillMesonDCATree) task->SetFillMesonDCATree(true);
  task->SetDoUseCentralEvtSelection(useCentralEvtSelection);
  task->SetForcePi0Unstable(setPi0Unstable);
  task->SetUseMixedBackAdd(enableAddBackground);
  task->SetDoRadiusDependence(enableRadiusDep);
  task->SetDoTrackingEff(doTrackingStudies);
  task->SetCutJetEnergyAsymm(cutOnJetEnergyAsymm);

  //connect containers
  TString nameContainer = Form("MesonJetCorrelation_Conv_%i_%i%s", meson, trainConfig, nameJetFinder.EqualTo("") == true ? "" : Form("_%s", nameJetFinder.Data()) );
  AliAnalysisDataContainer* coutput = mgr->CreateContainer(nameContainer, TList::Class(), AliAnalysisManager::kOutputContainer, Form("MJC_K0Lambda_%i_%i.root", meson, trainConfig));
  
  mgr->AddTask(task);
  mgr->ConnectInput(task, 0, cinput);
  mgr->ConnectOutput(task, 1, coutput);

  if(doFillMesonDCATree){
    for(int i = 0; i<numberOfCuts; i++){
      mgr->ConnectOutput(task,2+i,mgr->CreateContainer(Form("%s_%s_%s %s Meson DCA tree",(cuts.GetEventCut(i)).Data(),(cuts.GetPhotonCut(i)).Data(),(cuts.GetMesonCut(i)).Data(), nameJetFinder.Data()), TTree::Class(), AliAnalysisManager::kOutputContainer, Form("MJC_Co_%i_%i.root", meson, trainConfig)) );
    }
  }

  return;
}

