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

void AddTask_MesonJetCorr_Conv(
  Int_t trainConfig = 1,                // change different set of cuts
  Int_t isMC = 0,                       // run MC
  int meson = 0,                        // meson: 0=pi0, 1 = eta
  TString photonCutNumberV0Reader = "", // 00000008400000000100000000 nom. B, 00000088400000000100000000 low B
  TString periodNameV0Reader = "",
  // general setting for task
  Int_t enableQAMesonTask = 0,             // enable QA in AliAnalysisTaskGammaConvV1
  Int_t enableQAPhotonTask = 0,            // enable additional QA task
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

  //---------------------------------------
  // configs for pi0 meson pp 13 TeV
  //---------------------------------------
  if (trainConfig == 1) {
    cuts.AddCutPCM("00010103", "0dm00009f9730000dge0404000", "0152103500000000"); // config without jet requirement
  } else if (trainConfig == 2) {
    cuts.AddCutPCM("00010103", "0dm00009f9730000dge0404000", "2s52103500000000"); // in-Jet mass cut around pi0: 0.1-0.15, rotation back
  } else if (trainConfig == 3) {
    cuts.AddCutPCM("00010103", "0dm00009f9730000dge0404000", "2152103500000000"); // in-Jet mass cut around pi0: 0.1-0.15, mixed jet back
  } else if (trainConfig == 6) {
    cuts.AddCutPCM("000fc103", "0dm00009f9730000dge0404000", "2s52103500000000"); // Jet-low trigg in-Jet mass cut around pi0: 0.1-0.15, rotation back
  } else if (trainConfig == 7) {
    cuts.AddCutPCM("000fb103", "0dm00009f9730000dge0404000", "2s52103500000000"); // Jet-high trigg in-Jet mass cut around pi0: 0.1-0.15, rotation back
  } else if (trainConfig == 16) { // same as 6 but with mixed jet back
    cuts.AddCutPCM("000fc103", "0dm00009f9730000dge0404000", "2152103500000000"); // Jet-low trigg in-Jet mass cut around pi0: 0.1-0.15, mixed jet back
  } else if (trainConfig == 17) { // same as 7 but with mixed jet back
    cuts.AddCutPCM("000fb103", "0dm00009f9730000dge0404000", "2152103500000000"); // Jet-high trigg in-Jet mass cut around pi0: 0.1-0.15, mixed jet back
  
  } else if (trainConfig == 20) {
    cuts.AddCutPCM("00010103", "0dm00009f9730000dge0404000", "es52103500000000"); // decay daughters inside jet
  } else if (trainConfig == 21) {
    cuts.AddCutPCM("000fc103", "0dm00009f9730000dge0404000", "es52103500000000"); // decay daughters inside jet
    cuts.AddCutPCM("000fb103", "0dm00009f9730000dge0404000", "es52103500000000"); // decay daughters inside jet
    
  // configs with TRD/ITS conversion requirement
  } else if (trainConfig == 22) {
    cuts.AddCutPCM("00010103", "0dm00009f9730000dge0474000", "2s52103500000000"); // in-Jet mass cut around pi0: 0.1-0.15, rotation back
  } else if (trainConfig == 23) {
    cuts.AddCutPCM("000fc103", "0dm00009f9730000dge0474000", "2s52103500000000"); // in-Jet mass cut around pi0: 0.1-0.15, rotation back
  } else if (trainConfig == 24) {
    cuts.AddCutPCM("000fb103", "0dm00009f9730000dge0474000", "2s52103500000000"); // in-Jet mass cut around pi0: 0.1-0.15, rotation back
  
  // qt cut variations
  } else if (trainConfig == 25) {
    cuts.AddCutPCM("00010103", "0dm00009f97300003ge0404000", "2s52103500000000"); // qT max 0.05 1D
    cuts.AddCutPCM("00010103", "0dm00009f97300002ge0404000", "2s52103500000000"); // qT max 0.06 2D
    cuts.AddCutPCM("00010103", "0dm00009f97300009ge0404000", "2s52103500000000"); // qT max 0.03 2D


  // configs with eta < 0.5
  } else if (trainConfig == 30) {
    cuts.AddCutPCM("00010103", "0dm00009f9730000dge0404000", "2s52403500000000"); // in-Jet mass cut around pi0: 0.1-0.15, rotation back

    // configs with eta < 1.35
  } else if (trainConfig == 40) {
    cuts.AddCutPCM("00010103", "0dm00009f9730000dge0404000", "2s52003500000000"); // in-Jet mass cut around pi0: 0.1-0.15, rotation back
  

  //--- Systamtic variations for INT7 trigger
  } else if (trainConfig == 100) { // standard cut
    cuts.AddCutPCM("00010103", "0dm00009f9730000dge0404000", "2s52103500000000"); 
  } else if (trainConfig == 101) { // background variation
    cuts.AddCutPCM("00010103", "0dm00009f9730000dge0404000", "2s52103500000000");  // mixed back
  } else if (trainConfig == 102) { // alpha cut variation
    cuts.AddCutPCM("00010103", "0dm00009f9730000dge0404000", "2s52105500000000"); // alpha cut 0-0.75
    cuts.AddCutPCM("00010103", "0dm00009f9730000dge0404000", "2s52108500000000"); // alpha cut 0-0.65
  } else if (trainConfig == 103) { // opening angle var.
    cuts.AddCutPCM("00010103", "0dm00009f9730000dge0404000", "2s63103400000000"); // no opening angle cut

  //-- PCM variations
  } else if (trainConfig == 111) { // min pt electron variation
    cuts.AddCutPCM("00010103", "0dm00069f9730000dge0404000", "2s52103500000000"); // min pT 40 MeV
    cuts.AddCutPCM("00010103", "0dm00049f9730000dge0404000", "2s52103500000000"); // min pT 50 MeV
    cuts.AddCutPCM("00010103", "0dm00019f9730000dge0404000", "2s52103500000000"); // min pT 100 MeV
  } else if (trainConfig == 112) { // min pt electron variation
    cuts.AddCutPCM("00010103", "0dm00008f9730000dge0404000", "2s52103500000000"); // TPC cluster 35%
    cuts.AddCutPCM("00010103", "0dm00006f9730000dge0404000", "2s52103500000000"); // TPC cluster 70%
    cuts.AddCutPCM("00010103", "0dm00009f9730000dge0604000", "2s52103500000000"); // cosPA 0.9
    cuts.AddCutPCM("00010103", "0dm00009f9730000dge0304000", "2s52103500000000"); // cosPA 0.75
  } else if (trainConfig == 113) {
    cuts.AddCutPCM("00010103", "0dm0000939730000dge0404000", "2s52103500000000"); // nsig electron   -4,5
    cuts.AddCutPCM("00010103", "0dm0000969730000dge0404000", "2s52103500000000"); // nsig electron -2.5,4
    cuts.AddCutPCM("00010103", "0dm00009f5730000dge0404000", "2s52103500000000"); // nsig pion 2,-10
    cuts.AddCutPCM("00010103", "0dm00009f1730000dge0404000", "2s52103500000000"); // nsig pion 0,-10
  } else if (trainConfig == 114) {
    cuts.AddCutPCM("00010103", "0dm00009f9030000dge0404000", "2s52103500000000"); // pion nsig min mom 0.50 GeV/c
    cuts.AddCutPCM("00010103", "0dm00009f9630000dge0404000", "2s52103500000000"); // pion nsig min mom 0.25 GeV/c
    cuts.AddCutPCM("00010103", "0dm00009f9760000dge0404000", "2s52103500000000"); // pion nsig max mom 2.00 GeV/c
    cuts.AddCutPCM("00010103", "0dm00009f9710000dge0404000", "2s52103500000000"); // pion nsig max mom 5.00 GeV/c
  } else if (trainConfig == 115) {
    cuts.AddCutPCM("00010103", "0dm00009f97300008ge0404000", "2s52103500000000"); // qT max 0.05 1D
    cuts.AddCutPCM("00010103", "0dm00009f97300003ge0404000", "2s52103500000000"); // qT max 0.05 1D
    cuts.AddCutPCM("00010103", "0dm00009f97300002ge0404000", "2s52103500000000"); // qT max 0.06 1D
    cuts.AddCutPCM("00010103", "0dm00009f97300009ge0404000", "2s52103500000000"); // qT max 0.03 1D
  } else if (trainConfig == 116) {
    cuts.AddCutPCM("00010103", "0dm00009f9730000dg50404000", "2s52103500000000"); // Psi pair 0.1  1D
    cuts.AddCutPCM("00010103", "0dm00009f9730000dg10404000", "2s52103500000000"); // Psi pair 0.1  1D
    cuts.AddCutPCM("00010103", "0dm00009f9730000dg60404000", "2s52103500000000"); // Psi pair 0.05  1D
    cuts.AddCutPCM("00010103", "0dm00009f9730000dg80404000", "2s52103500000000"); // Psi pair 0.2  1D
  } else if (trainConfig == 117) {
    cuts.AddCutPCM("00010103", "0dm00009f9730000c259404000", "2s52103500000000"); // qT<0.110pT (2D) alpha<0.99
    cuts.AddCutPCM("00010103", "0dm00009f9730000a259404000", "2s52103500000000"); // qT<0.125pT (2D) alpha<0.99
    cuts.AddCutPCM("00010103", "0dm00009f9730000e259404000", "2s52103500000000"); // qT<0.130pT (2D) alpha<0.99
  } else if (trainConfig == 118) {
    cuts.AddCutPCM("00010103", "0dm00009f9730000dge0404000", "2s52103500000000"); // PsiPair<0.15exp(-0.065chi2)
    cuts.AddCutPCM("00010103", "0dm00009f9730000dge0404000", "2s52103500000000"); // PsiPair<0.18exp(-0.055chi2)
    cuts.AddCutPCM("00010103", "0dm00009f9730000dge0404000", "2s52103500000000"); // PsiPair<0.20exp(-0.050chi2)
  } else if (trainConfig == 119) {
    cuts.AddCutPCM("00010103", "0dm00009f9730000dge0474000", "2s52103500000000"); // config with ITS+TRDrequirement for electron tracks

    //---------------------------------------
    // configs for eta meson pp 13 TeV
    //---------------------------------------

  } else if (trainConfig == 1002) {
    cuts.AddCutPCM("00010103", "0dm00009f9730000dge0474000", "2r52103l00000000"); // in-Jet mass cut around eta: 0.5-0.6, rotation back
  } else if (trainConfig == 1003) {
    cuts.AddCutPCM("00010103", "0dm00009f9730000dge0474000", "2152103l00000000"); // in-Jet mass cut around pi0: 0.5-0.6, mixed jet back

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
  TList* ConvCutList = new TList();
  TList* MesonCutList = new TList();

  EventCutList->SetOwner(kTRUE);
  AliConvEventCuts** analysisEventCuts = new AliConvEventCuts*[numberOfCuts];
  ConvCutList->SetOwner(kTRUE);
  AliConversionPhotonCuts** analysisConvCuts = new AliConversionPhotonCuts*[numberOfCuts];
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
    //------------------- Conv Photon Cuts --------------------//
    //---------------------------------------------------------//
    analysisConvCuts[i] = new AliConversionPhotonCuts();

    if (enableMatBudWeightsPi0 > 0){
      if (isMC > 0){
        if (!analysisConvCuts[i]->InitializeMaterialBudgetWeights(enableMatBudWeightsPi0,fileNameMatBudWeights)){
          cout << "ERROR The initialization of the materialBudgetWeights did not work out." << endl;
          enableMatBudWeightsPi0 = false;
        }
      } else {
        cout << "ERROR 'enableMatBudWeightsPi0'-flag was set > 0 even though this is not a MC task. It was automatically reset to 0." << endl;
      }
    }

    analysisConvCuts[i]->SetV0ReaderName(V0ReaderName);
    if (enableElecDeDxPostCalibration) {
      if (isMC == 0) {
        if (fileNamedEdxPostCalib.CompareTo("") != 0) {
          analysisConvCuts[i]->SetElecDeDxPostCalibrationCustomFile(fileNamedEdxPostCalib);
          cout << "Setting custom dEdx recalibration file: " << fileNamedEdxPostCalib.Data() << endl;
        }
        analysisConvCuts[i]->SetDoElecDeDxPostCalibration(enableElecDeDxPostCalibration);
        cout << "Enabled TPC dEdx recalibration." << endl;
      } else {
        cout << "ERROR enableElecDeDxPostCalibration set to True even if MC file. Automatically reset to 0" << endl;
        enableElecDeDxPostCalibration = kFALSE;
        analysisConvCuts[i]->SetDoElecDeDxPostCalibration(kFALSE);
      }
    }
    if (enableLightOutput == 1 || enableLightOutput == 2 || enableLightOutput == 5)
      analysisConvCuts[i]->SetLightOutput(1);
    if (enableLightOutput == 4)
      analysisConvCuts[i]->SetLightOutput(2);
    if (enableLightOutput == 0)
      analysisConvCuts[i]->SetPlotTrackPID(kTRUE);
    analysisConvCuts[i]->InitializeCutsFromCutString((cuts.GetPhotonCut(i)).Data());
    analysisConvCuts[i]->SetIsHeavyIon(isHeavyIon);
    analysisConvCuts[i]->SetFillCutHistograms("", kFALSE);
    ConvCutList->Add(analysisConvCuts[i]);

    //---------------------------------------------------------//
    //------------------------ Meson Cuts ---------------------//
    //---------------------------------------------------------//
    analysisMesonCuts[i] = new AliConversionMesonCuts();
    analysisMesonCuts[i]->InitializeCutsFromCutString((cuts.GetMesonCut(i)).Data());
    analysisMesonCuts[i]->SetFillCutHistograms("");
    // analysisEventCuts[i]->SetAcceptedHeader(HeaderList);
    if (enableLightOutput > 0)
      analysisMesonCuts[i]->SetLightOutput(kTRUE);
    MesonCutList->Add(analysisMesonCuts[i]);
  }

  task->SetMesonKind(meson);
  task->SetMesonZPt(runOnlyZPt);
  task->SetIsConv(true);
  task->SetJetContainerAddName(nameJetFinder);
  task->SetEventCutList(numberOfCuts, EventCutList);
  task->SetMesonCutList(numberOfCuts, MesonCutList);
  task->SetConversionCutList(numberOfCuts, ConvCutList);
  task->SetDoMesonQA(enableQAMesonTask); 
  task->SetUseTHnSparseForResponse(enableTHnSparse);
  if(enableMatBudWeightsPi0) task->SetDoMaterialBudgetWeightingOfGammasForTrueMesons(true);
  if(doFillMesonDCATree) task->SetFillMesonDCATree(true);
  task->SetDoUseCentralEvtSelection(useCentralEvtSelection);
  task->SetForcePi0Unstable(setPi0Unstable);
  task->SetUseMixedBackAdd(enableAddBackground);
  task->SetDoRadiusDependence(enableRadiusDep);

  //connect containers
  TString nameContainer = Form("MesonJetCorrelation_Conv_%i_%i%s", meson, trainConfig, nameJetFinder.EqualTo("") == true ? "" : Form("_%s", nameJetFinder.Data()) );
  AliAnalysisDataContainer* coutput = mgr->CreateContainer(nameContainer, TList::Class(), AliAnalysisManager::kOutputContainer, Form("MesonJetCorrelation_Conv_%i_%i.root", meson, trainConfig));
  
  mgr->AddTask(task);
  mgr->ConnectInput(task, 0, cinput);
  mgr->ConnectOutput(task, 1, coutput);

  if(doFillMesonDCATree){
    for(int i = 0; i<numberOfCuts; i++){
      mgr->ConnectOutput(task,2+i,mgr->CreateContainer(Form("%s_%s_%s %s Meson DCA tree",(cuts.GetEventCut(i)).Data(),(cuts.GetPhotonCut(i)).Data(),(cuts.GetMesonCut(i)).Data(), nameJetFinder.Data()), TTree::Class(), AliAnalysisManager::kOutputContainer, Form("MesonJetCorrelation_Conv_%i_%i.root", meson, trainConfig)) );
    }
  }

  return;
}

