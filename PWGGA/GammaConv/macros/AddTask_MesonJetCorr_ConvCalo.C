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

void AddTask_MesonJetCorr_ConvCalo(
  Int_t trainConfig = 1,                // change different set of cuts
  Int_t isMC = 0,                       // run MC
  int meson = 0,                        // meson: 0=pi0, 1 = eta
  TString photonCutNumberV0Reader = "", // 00000008400000000100000000 nom. B, 00000088400000000100000000 low B
  TString periodNameV0Reader = "",
  // general setting for task
  Int_t enableQAMesonTask = 0,             // enable QA in AliAnalysisTaskGammaConvV1
  Int_t enableQAPhotonTask = 0,            // enable additional QA task
  Int_t enableExtMatchAndQA = 0,           // disabled (0), extMatch (1), extQA_noCellQA (2), extMatch+extQA_noCellQA (3), extQA+cellQA (4), extMatch+extQA+cellQA (5)
  Int_t enableLightOutput = 0,             // switch to run light output (only essential histograms for afterburner)
  Bool_t enableTHnSparse = kFALSE,         // switch on THNsparse
  Int_t enableTriggerMimicking = 0,        // enable trigger mimicking
  Bool_t enableTriggerOverlapRej = kFALSE, // enable trigger overlap rejection
  TString settingMaxFacPtHard = "3.",      // maximum factor between hardest jet and ptHard generated
  Int_t debugLevel = 0,                    // introducing debug levels for grid running
  // settings for weights
  // FPTW:fileNamePtWeights, FMUW:fileNameMultWeights,  FMAW:fileNameMatBudWeights,  separate with ;
  // Material Budget Weights file for Run 2
  // FMAW:alien:///alice/cern.ch/user/a/amarin//MBW/MCInputFileMaterialBudgetWeightsLHC16_Pythia_00010103_0d000009266300008850404000_date181214.root
  TString fileNameExternalInputs = "",
  Int_t doWeightingPart = 0,                   // enable Weighting
  TString generatorName = "DPMJET",            // generator Name
  Bool_t enableMultiplicityWeighting = kFALSE, //
  Int_t enableMatBudWeightsPi0 = 0,            // 1 = three radial bins, 2 = 10 radial bins (2 is the default when using weights)
  Bool_t enableElecDeDxPostCalibration = kFALSE,
  TString periodNameAnchor = "", //
  // special settings
  Bool_t enableSortingMCLabels = kTRUE,     // enable sorting for MC cluster labels
  Bool_t enableTreeConvGammaShape = kFALSE, // enable additional tree for conversion properties for clusters
  Bool_t doSmear = kFALSE,                  // switches to run user defined smearing
  Double_t bremSmear = 1.,
  Double_t smearPar = 0.,                // conv photon smearing params
  Double_t smearParConst = 0.,           // conv photon smearing params
  Bool_t doPrimaryTrackMatching = kTRUE, // enable basic track matching for all primary tracks to cluster
  bool useCentralEvtSelection = true,
  bool setPi0Unstable = false,
  bool enableAddBackground = false,
  bool enableRadiusDep = false,
  int runOnlyZPt = 0,           // if 0, bot pt and z histograms will be filled, if 1, pt histograms will be filled, if 2, only z histograms will be filled
  // subwagon config
  TString additionalTrainConfig = "0" // additional counter for trainconfig
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

  TString corrTaskSetting = cuts.GetSpecialSettingFromAddConfig(additionalTrainConfig, "CF", "", addTaskName);
  if (corrTaskSetting.CompareTo(""))
    cout << "corrTaskSetting: " << corrTaskSetting.Data() << endl;

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
  task = new AliAnalysisTaskMesonJetCorrelation(Form("MesonJetCorrelation_ConvCalo_%i_%i", meson, trainConfig));
  task->SetIsHeavyIon(isHeavyIon);
  task->SetIsMC(isMC);
  task->SetV0ReaderName(V0ReaderName);
  task->SetCorrectionTaskSetting(corrTaskSetting);
  task->SetTrackMatcherRunningMode(trackMatcherRunningMode); // have to do this!

  //---------------------------------------
  // configs for pi0 meson pp 13 TeV
  //---------------------------------------
  // configs without NonLinearity as the NonLin is already applied in the CF
  if (trainConfig == 1) {
    cuts.AddCutPCMCalo("00010103", "0dm00009f9730000dge0404000", "411790009fe30230000", "0r63103100000010"); // test config
  } else if (trainConfig == 2) {
    cuts.AddCutPCMCalo("00010103", "0dm00009f9730000dge0404000", "411790009fe30230000", "2s63103400000010"); // in-Jet, mass cut pi0: 0.1-0.15, rotation back
  } else if (trainConfig == 3) {
    cuts.AddCutPCMCalo("00010103", "0dm00009f9730000dge0404000", "411790009fe30230000", "2163103400000010"); // in-Jet, mass cut pi0: 0.1-0.15, mixed jet back
  } else if (trainConfig == 4) {
    cuts.AddCutPCMCalo("0008e103", "0dm00009f9730000dge0404000", "411790009fe30230000", "2s63103400000010"); // EG2 in-Jet, mass cut pi0: 0.1-0.15, rotation back
  } else if (trainConfig == 5) {
    cuts.AddCutPCMCalo("0008d103", "0dm00009f9730000dge0404000", "411790009fe30230000", "2s63103400000010"); // EG1 in-Jet, mass cut pi0: 0.1-0.15, rotation back
  } else if (trainConfig == 6) {
    cuts.AddCutPCMCalo("000fc103", "0dm00009f9730000dge0404000", "411790009fe30230000", "2s63103400000010"); // Jet low trigg in-Jet, mass cut pi0: 0.1-0.15, rotation back
  } else if (trainConfig == 7) {
    cuts.AddCutPCMCalo("000fb103", "0dm00009f9730000dge0404000", "411790009fe30230000", "2s63103400000010"); // Jet high trigg in-Jet, mass cut pi0: 0.1-0.15, rotation back

  // EJ1 and EJ2 only EMCal triggers
  } else if (trainConfig == 8) {
    cuts.AddCutPCMCalo("000f5103", "0dm00009f9730000dge0404000", "411790009fe30230000", "2s631034000000d0"); // Jet-low trigg in-jet, pi0 mass: 0.1-0.15, rotation back
  } else if (trainConfig == 9) {
    cuts.AddCutPCMCalo("000f3103", "0dm00009f9730000dge0404000", "411790009fe30230000", "2s631034000000d0"); // Jet-high trigg in-jet, pi0 mass: 0.1-0.15, rotation back

  // configs with Mesons only in EMCal, EMCal triggers only
  } else if (trainConfig == 10) {
    cuts.AddCutPCMCalo("00010103", "0dm00009f9730000dge0404000", "111110009fe30230000", "2s631034000000d0"); // in-jet, pi0 mass: 0.1-0.15, rotation back
  } else if (trainConfig == 11) {
    cuts.AddCutPCMCalo("000f5103", "0dm00009f9730000dge0404000", "111110009fe30230000", "2s631034000000d0"); // Jet-low trigg in-jet, pi0 mass: 0.1-0.15, rotation back
  } else if (trainConfig == 12) {
    cuts.AddCutPCMCalo("000f3103", "0dm00009f9730000dge0404000", "111110009fe30230000", "2s631034000000d0"); // Jet-high trigg in-jet, pi0 mass: 0.1-0.15, rotation back

  } else if (trainConfig == 14) { // same as 4 but with jet mixing back
    cuts.AddCutPCMCalo("0008e103", "0dm00009f9730000dge0404000", "411790009fe30230000", "2163103400000010"); // EG2 in-Jet, mass cut pi0: 0.1-0.15, mixed jet back
  } else if (trainConfig == 15) { // same as 5 but with jet mixing back
    cuts.AddCutPCMCalo("0008d103", "0dm00009f9730000dge0404000", "411790009fe30230000", "2163103400000010"); // EG1 in-Jet, mass cut pi0: 0.1-0.15, mixed jet back
  } else if (trainConfig == 16) { // same as 6 but with jet mixing back
    cuts.AddCutPCMCalo("000fc103", "0dm00009f9730000dge0404000", "411790009fe30230000", "2163103400000010"); // Jet low trigg in-Jet, mass cut pi0: 0.1-0.15, mixed jet back
  } else if (trainConfig == 17) { // same as 7 but with jet mixing back
    cuts.AddCutPCMCalo("000fb103", "0dm00009f9730000dge0404000", "411790009fe30230000", "2163103400000010"); // Jet high trigg in-Jet, mass cut pi0: 0.1-0.15, mixed jet back
  

  } else if (trainConfig == 20) {
    cuts.AddCutPCMCalo("00010103", "0dm00009f9730000dge0404000", "411790009fe30230000", "es63103400000010"); // decay daughters also inside jet
  } else if (trainConfig == 21) {
    cuts.AddCutPCMCalo("000fc103", "0dm00009f9730000dge0404000", "411790009fe30230000", "es63103400000010"); // decay daughters also inside jet
    cuts.AddCutPCMCalo("000fb103", "0dm00009f9730000dge0404000", "411790009fe30230000", "es63103400000010"); // decay daughters also inside jet


  // configs with NonLinearity 
  } else if (trainConfig == 22) {
    cuts.AddCutPCMCalo("00010103", "0dm00009f9730000dge0404000", "411790109fe30230000", "2s63103400000010"); // in-Jet, mass cut pi0: 0.1-0.15, rotation back
  } else if (trainConfig == 23) {
    cuts.AddCutPCMCalo("00010103", "0dm00009f9730000dge0404000", "411790109fe30230000", "2163103400000010"); // in-Jet, mass cut pi0: 0.1-0.15, mixed jet back
  } else if (trainConfig == 24) {
    cuts.AddCutPCMCalo("0008e103", "0dm00009f9730000dge0404000", "411790109fe30230000", "2s63103400000010"); // EG2 in-Jet, mass cut pi0: 0.1-0.15, rotation back
  } else if (trainConfig == 25) {
    cuts.AddCutPCMCalo("0008d103", "0dm00009f9730000dge0404000", "411790109fe30230000", "2s63103400000010"); // EG1 in-Jet, mass cut pi0: 0.1-0.15, rotation back
  } else if (trainConfig == 26) {
    cuts.AddCutPCMCalo("000fc103", "0dm00009f9730000dge0404000", "411790109fe30230000", "2s63103400000010"); // Jet low trigg in-Jet, mass cut pi0: 0.1-0.15, rotation back
  } else if (trainConfig == 27) {
    cuts.AddCutPCMCalo("000fb103", "0dm00009f9730000dge0404000", "411790109fe30230000", "2s63103400000010"); // Jet high trigg in-Jet, mass cut pi0: 0.1-0.15, rotation back


   // configs with eta < 0.5
  } else if (trainConfig == 30) {
    cuts.AddCutPCMCalo("00010103", "0dm00009f9730000dge0404000", "411790009fe30230000", "2s634034000000d0"); // in-jet, pi0 mass: 0.1-0.15, rotation back
  } else if (trainConfig == 31) {
    cuts.AddCutPCMCalo("000fc103", "0dm00009f9730000dge0404000", "411790009fe30230000", "2s634034000000d0"); // Jet-low trigg in-jet, pi0 mass: 0.1-0.15, rotation back
  } else if (trainConfig == 32) {
    cuts.AddCutPCMCalo("000fb103", "0dm00009f9730000dge0404000", "411790009fe30230000", "2s634034000000d0"); // Jet-high trigg in-jet, pi0 mass: 0.1-0.15, rotation back

  // configs with eta < 0.5, only EMCal triggers (not DCal)
  } else if (trainConfig == 33) {
    cuts.AddCutPCMCalo("000f5103", "0dm00009f9730000dge0404000", "411790009fe30230000", "2s634034000000d0"); // Jet-low trigg in-jet, pi0 mass: 0.1-0.15, rotation back
  } else if (trainConfig == 34) {
    cuts.AddCutPCMCalo("000f3103", "0dm00009f9730000dge0404000", "411790009fe30230000", "2s634034000000d0"); // Jet-high trigg in-jet, pi0 mass: 0.1-0.15, rotation back

  // configs with eta < 0.5, only EMCal triggers (not DCal), Mesons only with EMCal
  } else if (trainConfig == 35) {
    cuts.AddCutPCMCalo("000f5103", "0dm00009f9730000dge0404000", "111110009fe30230000", "2s634034000000d0"); // Jet-low trigg in-jet, pi0 mass: 0.1-0.15, rotation back
  } else if (trainConfig == 36) {
    cuts.AddCutPCMCalo("000f3103", "0dm00009f9730000dge0404000", "111110009fe30230000", "2s634034000000d0"); // Jet-high trigg in-jet, pi0 mass: 0.1-0.15, rotation back

   // configs with eta < 1.35
  } else if (trainConfig == 40) {
    cuts.AddCutPCMCalo("00010103", "0dm00009f9730000dge0404000", "411790009fe30230000", "2s630034000000d0"); // in-jet, pi0 mass: 0.1-0.15, rotation back
  } else if (trainConfig == 41) {
    cuts.AddCutPCMCalo("000fc103", "0dm00009f9730000dge0404000", "411790009fe30230000", "2s630034000000d0"); // Jet-low trigg in-jet, pi0 mass: 0.1-0.15, rotation back
  } else if (trainConfig == 42) {
    cuts.AddCutPCMCalo("000fb103", "0dm00009f9730000dge0404000", "411790009fe30230000", "2s630034000000d0"); // Jet-high trigg in-jet, pi0 mass: 0.1-0.15, rotation back



  // ------ cut variations
  } else if (trainConfig == 100) { // 
    cuts.AddCutPCMCalo("00010103", "0dm00009f9730000dge0404000", "411790009fe30230000", "2s63103400000010"); // NL var. etc which is handled in correction framework
  } else if (trainConfig == 101) { // background variation
    cuts.AddCutPCMCalo("00010103", "0dm00009f9730000dge0404000", "411790009fe30230000", "2163103400000010"); // jet mixing back
  } else if (trainConfig == 102) { // alpha cut variation
    cuts.AddCutPCMCalo("00010103", "0dm00009f9730000dge0404000", "411790009fe30230000", "2s63105400000010"); // alpha cut 0-0.75
    cuts.AddCutPCMCalo("00010103", "0dm00009f9730000dge0404000", "411790009fe30230000", "2s63108400000010"); // alpha cut 0-0.65
  } else if (trainConfig == 103) { // opening angle var.
    cuts.AddCutPCMCalo("00010103", "0dm00009f9730000dge0404000", "411790009fe30230000", "2s63103400000000"); // no opening angle cut


  } else if (trainConfig == 105) { // M02 variation (std M02 = 0.5)
    cuts.AddCutPCMCalo("00010103", "0dm00009f9730000dge0404000", "411790009fe30240000", "2s63103400000010"); // M02 = 0.4
    cuts.AddCutPCMCalo("00010103", "0dm00009f9730000dge0404000", "411790009fe30220000", "2s63103400000010"); // M02 = 0.7
    cuts.AddCutPCMCalo("00010103", "0dm00009f9730000dge0404000", "411790009fe30210000", "2s63103400000010"); // M02 = 1.0
    cuts.AddCutPCMCalo("00010103", "0dm00009f9730000dge0404000", "411790009fe302v0000", "2s63103400000010"); // 0.5 < M02 < 0.7
  } else if (trainConfig == 106) { // TM variations for mesons
    cuts.AddCutPCMCalo("00010103", "0dm00009f9730000dge0404000", "4117900090e30230000", "2s63103400000010"); // no TM
    cuts.AddCutPCMCalo("00010103", "0dm00009f9730000dge0404000", "411790009ee30230000", "2s63103400000010"); // TM var EoverP 2.00
    cuts.AddCutPCMCalo("00010103", "0dm00009f9730000dge0404000", "411790009ge30230000", "2s63103400000010"); // TM var EoverP 1.5
    cuts.AddCutPCMCalo("00010103", "0dm00009f9730000dge0404000", "4117900097e30230000", "2s63103400000010"); // No E/p
    cuts.AddCutPCMCalo("00010103", "0dm00009f9730000dge0404000", "411790009le30230000", "2s63103400000010"); // std + sec TM
    cuts.AddCutPCMCalo("00010103", "0dm00009f9730000dge0404000", "411790009ne30230000", "2s63103400000010"); // TM var, Eta (0.035, 0.010, 2.5); Phi (0.085, 0.015, 2.)
  } else if (trainConfig == 107) { // cluster time
    cuts.AddCutPCMCalo("00010103", "0dm00009f9730000dge0404000", "411790005fe30230000", "2s63103400000010"); // -50 - 50 ns
    cuts.AddCutPCMCalo("00010103", "0dm00009f9730000dge0404000", "411790006fe30230000", "2s63103400000010"); // -30 - 35 ns
    cuts.AddCutPCMCalo("00010103", "0dm00009f9730000dge0404000", "41179000afe30230000", "2s63103400000010"); // -12.5 - 13 ns
  } else if (trainConfig == 108) { // NCell
    cuts.AddCutPCMCalo("00010103", "0dm00009f9730000dge0404000", "411790009fe3r230000", "2s63103400000010"); // EDC pi0 tagging for gamma clus, Gaussian Fit
    cuts.AddCutPCMCalo("00010103", "0dm00009f9730000dge0404000", "411790009fe3n230000", "2s63103400000010"); // PCMEDC pi0 tagging for gamma clus, Gaussian Fit
    cuts.AddCutPCMCalo("00010103", "0dm00009f9730000dge0404000", "411790009fe3m230000", "2s63103400000010"); // PCMEDC pi0 tagging for all clus, Gaussian Fit
    cuts.AddCutPCMCalo("00010103", "0dm00009f9730000dge0404000", "411790009fe3l230000", "2s63103400000010"); // PCMEDC pi0 tagging for gamma clus, pol2 fit
  } else if (trainConfig == 109) { // min energy
    cuts.AddCutPCMCalo("00010103", "0dm00009f9730000dge0404000", "411790009fe10230000", "2s63103400000010"); // INT7, minE = 0.5
    cuts.AddCutPCMCalo("00010103", "0dm00009f9730000dge0404000", "411790009fe20230000", "2s63103400000010"); // INT7, minE = 0.6
    cuts.AddCutPCMCalo("00010103", "0dm00009f9730000dge0404000", "411790009fe40230000", "2s63103400000010"); // INT7, minE = 0.8
  } else if (trainConfig == 110) { // Exotics
    cuts.AddCutPCMCalo("00010103", "0dm00009f9730000dge0404000", "411790009f030230000", "2s63103400000010"); // no exotics
    cuts.AddCutPCMCalo("00010103", "0dm00009f9730000dge0404000", "411790009fb30230000", "2s63103400000010"); // F+ < 0.95

  } else if (trainConfig == 111) { // min pt electron variation
    cuts.AddCutPCMCalo("00010103", "0dm00069f9730000dge0404000", "411790009fe30230000", "2s63103400000010"); // min pT 40 MeV
    cuts.AddCutPCMCalo("00010103", "0dm00049f9730000dge0404000", "411790009fe30230000", "2s63103400000010"); // min pT 50 MeV
    cuts.AddCutPCMCalo("00010103", "0dm00019f9730000dge0404000", "411790009fe30230000", "2s63103400000010"); // min pT 100 MeV
  } else if (trainConfig == 112) { // min pt electron variation
    cuts.AddCutPCMCalo("00010103", "0dm00008f9730000dge0404000", "411790009fe30230000", "2s63103400000010"); // TPC cluster 35%
    cuts.AddCutPCMCalo("00010103", "0dm00006f9730000dge0404000", "411790009fe30230000", "2s63103400000010"); // TPC cluster 70%
    cuts.AddCutPCMCalo("00010103", "0dm00009f9730000dge0604000", "411790009fe30230000", "2s63103400000010"); // cosPA 0.9
    cuts.AddCutPCMCalo("00010103", "0dm00009f9730000dge0304000", "411790009fe30230000", "2s63103400000010"); // cosPA 0.75
  } else if (trainConfig == 113) {
    cuts.AddCutPCMCalo("00010103", "0dm0000939730000dge0404000", "411790009fe30230000", "2s63103400000010"); // nsig electron   -4,5
    cuts.AddCutPCMCalo("00010103", "0dm0000969730000dge0404000", "411790009fe30230000", "2s63103400000010"); // nsig electron -2.5,4
    cuts.AddCutPCMCalo("00010103", "0dm00009f5730000dge0404000", "411790009fe30230000", "2s63103400000010"); // nsig pion 2,-10
    cuts.AddCutPCMCalo("00010103", "0dm00009f1730000dge0404000", "411790009fe30230000", "2s63103400000010"); // nsig pion 0,-10
  } else if (trainConfig == 114) {
    cuts.AddCutPCMCalo("00010103", "0dm00009f9030000dge0404000", "411790009fe30230000", "2s63103400000010"); // pion nsig min mom 0.50 GeV/c
    cuts.AddCutPCMCalo("00010103", "0dm00009f9630000dge0404000", "411790009fe30230000", "2s63103400000010"); // pion nsig min mom 0.25 GeV/c
    cuts.AddCutPCMCalo("00010103", "0dm00009f9760000dge0404000", "411790009fe30230000", "2s63103400000010"); // pion nsig max mom 2.00 GeV/c
    cuts.AddCutPCMCalo("00010103", "0dm00009f9710000dge0404000", "411790009fe30230000", "2s63103400000010"); // pion nsig max mom 5.00 GeV/c
  } else if (trainConfig == 115) {
    cuts.AddCutPCMCalo("00010103", "0dm00009f97300008ge0404000", "411790009fe30230000", "2s63103400000010"); // qT max 0.05 1D
    cuts.AddCutPCMCalo("00010103", "0dm00009f97300003ge0404000", "411790009fe30230000", "2s63103400000010"); // qT max 0.05 1D
    cuts.AddCutPCMCalo("00010103", "0dm00009f97300002ge0404000", "411790009fe30230000", "2s63103400000010"); // qT max 0.06 1D
    cuts.AddCutPCMCalo("00010103", "0dm00009f97300009ge0404000", "411790009fe30230000", "2s63103400000010"); // qT max 0.03 1D
  } else if (trainConfig == 116) {
    cuts.AddCutPCMCalo("00010103", "0dm00009f9730000dg50404000", "411790009fe30230000", "2s63103400000010"); // Psi pair 0.1  1D
    cuts.AddCutPCMCalo("00010103", "0dm00009f9730000dg10404000", "411790009fe30230000", "2s63103400000010"); // Psi pair 0.1  1D
    cuts.AddCutPCMCalo("00010103", "0dm00009f9730000dg60404000", "411790009fe30230000", "2s63103400000010"); // Psi pair 0.05  1D
    cuts.AddCutPCMCalo("00010103", "0dm00009f9730000dg80404000", "411790009fe30230000", "2s63103400000010"); // Psi pair 0.2  1D
  } else if (trainConfig == 117) {
    cuts.AddCutPCMCalo("00010103", "0dm00009f9730000c259404000", "411790009fe30230000", "2s63103400000010"); // qT<0.110pT (2D) alpha<0.99
    cuts.AddCutPCMCalo("00010103", "0dm00009f9730000a259404000", "411790009fe30230000", "2s63103400000010"); // qT<0.125pT (2D) alpha<0.99
    cuts.AddCutPCMCalo("00010103", "0dm00009f9730000e259404000", "411790009fe30230000", "2s63103400000010"); // qT<0.130pT (2D) alpha<0.99
  } else if (trainConfig == 118) {
    cuts.AddCutPCMCalo("00010103", "0dm00009f9730000dge0404000", "411790009fe30230000", "2s63103400000010"); // PsiPair<0.15exp(-0.065chi2)
    cuts.AddCutPCMCalo("00010103", "0dm00009f9730000dge0404000", "411790009fe30230000", "2s63103400000010"); // PsiPair<0.18exp(-0.055chi2)
    cuts.AddCutPCMCalo("00010103", "0dm00009f9730000dge0404000", "411790009fe30230000", "2s63103400000010"); // PsiPair<0.20exp(-0.050chi2)


  // ------ EJ2 cut variations
  } else if (trainConfig == 130) { // 
    cuts.AddCutPCMCalo("000fc103", "0dm00009f9730000dge0404000", "411790009fe30230000", "2s63103400000010"); // NL var. etc which is handled in correction framework
  } else if (trainConfig == 131) { // background variation
    cuts.AddCutPCMCalo("000fc103", "0dm00009f9730000dge0404000", "411790009fe30230000", "2163103400000010"); // jet mixing back
  } else if (trainConfig == 132) { // alpha cut variation
    cuts.AddCutPCMCalo("000fc103", "0dm00009f9730000dge0404000", "411790009fe30230000", "2s63105400000010"); // alpha cut 0-0.75
    cuts.AddCutPCMCalo("000fc103", "0dm00009f9730000dge0404000", "411790009fe30230000", "2s63108400000010"); // alpha cut 0-0.65
  } else if (trainConfig == 133) { // opening angle var.
    cuts.AddCutPCMCalo("000fc103", "0dm00009f9730000dge0404000", "411790009fe30230000", "2s63103400000000"); // no opening angle cut


  } else if (trainConfig == 135) { // M02 variation (std M02 = 0.5)
    cuts.AddCutPCMCalo("000fc103", "0dm00009f9730000dge0404000", "411790009fe30240000", "2s63103400000010"); // M02 = 0.4
    cuts.AddCutPCMCalo("000fc103", "0dm00009f9730000dge0404000", "411790009fe30220000", "2s63103400000010"); // M02 = 0.7
    cuts.AddCutPCMCalo("000fc103", "0dm00009f9730000dge0404000", "411790009fe30210000", "2s63103400000010"); // M02 = 1.0
    cuts.AddCutPCMCalo("000fc103", "0dm00009f9730000dge0404000", "411790009fe302v0000", "2s63103400000010"); // 0.5 < M02 < 0.7
  } else if (trainConfig == 136) { // TM variations for mesons
    cuts.AddCutPCMCalo("000fc103", "0dm00009f9730000dge0404000", "4117900090e30230000", "2s63103400000010"); // no TM
    cuts.AddCutPCMCalo("000fc103", "0dm00009f9730000dge0404000", "411790009ee30230000", "2s63103400000010"); // TM var EoverP 2.00
    cuts.AddCutPCMCalo("000fc103", "0dm00009f9730000dge0404000", "411790009ge30230000", "2s63103400000010"); // TM var EoverP 1.5
    cuts.AddCutPCMCalo("000fc103", "0dm00009f9730000dge0404000", "4117900097e30230000", "2s63103400000010"); // No E/p
    cuts.AddCutPCMCalo("000fc103", "0dm00009f9730000dge0404000", "411790009le30230000", "2s63103400000010"); // std + sec TM
    cuts.AddCutPCMCalo("000fc103", "0dm00009f9730000dge0404000", "411790009ne30230000", "2s63103400000010"); // TM var, Eta (0.035, 0.010, 2.5); Phi (0.085, 0.015, 2.)
  } else if (trainConfig == 137) { // cluster time
    cuts.AddCutPCMCalo("000fc103", "0dm00009f9730000dge0404000", "411790005fe30230000", "2s63103400000010"); // -50 - 50 ns
    cuts.AddCutPCMCalo("000fc103", "0dm00009f9730000dge0404000", "411790006fe30230000", "2s63103400000010"); // -30 - 35 ns
    cuts.AddCutPCMCalo("000fc103", "0dm00009f9730000dge0404000", "41179000afe30230000", "2s63103400000010"); // -12.5 - 13 ns
  } else if (trainConfig == 138) { // NCell
    cuts.AddCutPCMCalo("000fc103", "0dm00009f9730000dge0404000", "411790009fe3r230000", "2s63103400000010"); // EDC pi0 tagging for gamma clus, Gaussian Fit
    cuts.AddCutPCMCalo("000fc103", "0dm00009f9730000dge0404000", "411790009fe3n230000", "2s63103400000010"); // PCMEDC pi0 tagging for gamma clus, Gaussian Fit
    cuts.AddCutPCMCalo("000fc103", "0dm00009f9730000dge0404000", "411790009fe3m230000", "2s63103400000010"); // PCMEDC pi0 tagging for all clus, Gaussian Fit
    cuts.AddCutPCMCalo("000fc103", "0dm00009f9730000dge0404000", "411790009fe3l230000", "2s63103400000010"); // PCMEDC pi0 tagging for gamma clus, pol2 fit
  } else if (trainConfig == 139) { // min energy
    cuts.AddCutPCMCalo("000fc103", "0dm00009f9730000dge0404000", "411790009fe10230000", "2s63103400000010"); // INT7, minE = 0.5
    cuts.AddCutPCMCalo("000fc103", "0dm00009f9730000dge0404000", "411790009fe20230000", "2s63103400000010"); // INT7, minE = 0.6
    cuts.AddCutPCMCalo("000fc103", "0dm00009f9730000dge0404000", "411790009fe40230000", "2s63103400000010"); // INT7, minE = 0.8
  } else if (trainConfig == 140) { // Exotics
    cuts.AddCutPCMCalo("000fc103", "0dm00009f9730000dge0404000", "411790009f030230000", "2s63103400000010"); // no exotics
    cuts.AddCutPCMCalo("000fc103", "0dm00009f9730000dge0404000", "411790009fb30230000", "2s63103400000010"); // F+ < 0.95

  } else if (trainConfig == 141) { // min pt electron variation
    cuts.AddCutPCMCalo("000fc103", "0dm00069f9730000dge0404000", "411790009fe30230000", "2s63103400000010"); // min pT 40 MeV
    cuts.AddCutPCMCalo("000fc103", "0dm00049f9730000dge0404000", "411790009fe30230000", "2s63103400000010"); // min pT 50 MeV
    cuts.AddCutPCMCalo("000fc103", "0dm00019f9730000dge0404000", "411790009fe30230000", "2s63103400000010"); // min pT 100 MeV
  } else if (trainConfig == 142) { // min pt electron variation
    cuts.AddCutPCMCalo("000fc103", "0dm00008f9730000dge0404000", "411790009fe30230000", "2s63103400000010"); // TPC cluster 35%
    cuts.AddCutPCMCalo("000fc103", "0dm00006f9730000dge0404000", "411790009fe30230000", "2s63103400000010"); // TPC cluster 70%
    cuts.AddCutPCMCalo("000fc103", "0dm00009f9730000dge0604000", "411790009fe30230000", "2s63103400000010"); // cosPA 0.9
    cuts.AddCutPCMCalo("000fc103", "0dm00009f9730000dge0304000", "411790009fe30230000", "2s63103400000010"); // cosPA 0.75
  } else if (trainConfig == 143) {
    cuts.AddCutPCMCalo("000fc103", "0dm0000939730000dge0404000", "411790009fe30230000", "2s63103400000010"); // nsig electron   -4,5
    cuts.AddCutPCMCalo("000fc103", "0dm0000969730000dge0404000", "411790009fe30230000", "2s63103400000010"); // nsig electron -2.5,4
    cuts.AddCutPCMCalo("000fc103", "0dm00009f5730000dge0404000", "411790009fe30230000", "2s63103400000010"); // nsig pion 2,-10
    cuts.AddCutPCMCalo("000fc103", "0dm00009f1730000dge0404000", "411790009fe30230000", "2s63103400000010"); // nsig pion 0,-10
  } else if (trainConfig == 144) {
    cuts.AddCutPCMCalo("000fc103", "0dm00009f9030000dge0404000", "411790009fe30230000", "2s63103400000010"); // pion nsig min mom 0.50 GeV/c
    cuts.AddCutPCMCalo("000fc103", "0dm00009f9630000dge0404000", "411790009fe30230000", "2s63103400000010"); // pion nsig min mom 0.25 GeV/c
    cuts.AddCutPCMCalo("000fc103", "0dm00009f9760000dge0404000", "411790009fe30230000", "2s63103400000010"); // pion nsig max mom 2.00 GeV/c
    cuts.AddCutPCMCalo("000fc103", "0dm00009f9710000dge0404000", "411790009fe30230000", "2s63103400000010"); // pion nsig max mom 5.00 GeV/c
  } else if (trainConfig == 145) {
    cuts.AddCutPCMCalo("000fc103", "0dm00009f97300008ge0404000", "411790009fe30230000", "2s63103400000010"); // qT max 0.05 1D
    cuts.AddCutPCMCalo("000fc103", "0dm00009f97300003ge0404000", "411790009fe30230000", "2s63103400000010"); // qT max 0.05 1D
    cuts.AddCutPCMCalo("000fc103", "0dm00009f97300002ge0404000", "411790009fe30230000", "2s63103400000010"); // qT max 0.06 1D
    cuts.AddCutPCMCalo("000fc103", "0dm00009f97300009ge0404000", "411790009fe30230000", "2s63103400000010"); // qT max 0.03 1D
  } else if (trainConfig == 146) {
    cuts.AddCutPCMCalo("000fc103", "0dm00009f9730000dg50404000", "411790009fe30230000", "2s63103400000010"); // Psi pair 0.1  1D
    cuts.AddCutPCMCalo("000fc103", "0dm00009f9730000dg10404000", "411790009fe30230000", "2s63103400000010"); // Psi pair 0.1  1D
    cuts.AddCutPCMCalo("000fc103", "0dm00009f9730000dg60404000", "411790009fe30230000", "2s63103400000010"); // Psi pair 0.05  1D
    cuts.AddCutPCMCalo("000fc103", "0dm00009f9730000dg80404000", "411790009fe30230000", "2s63103400000010"); // Psi pair 0.2  1D
  } else if (trainConfig == 147) {
    cuts.AddCutPCMCalo("000fc103", "0dm00009f9730000c259404000", "411790009fe30230000", "2s63103400000010"); // qT<0.110pT (2D) alpha<0.99
    cuts.AddCutPCMCalo("000fc103", "0dm00009f9730000a259404000", "411790009fe30230000", "2s63103400000010"); // qT<0.125pT (2D) alpha<0.99
    cuts.AddCutPCMCalo("000fc103", "0dm00009f9730000e259404000", "411790009fe30230000", "2s63103400000010"); // qT<0.130pT (2D) alpha<0.99
  } else if (trainConfig == 148) {
    cuts.AddCutPCMCalo("000fc103", "0dm00009f9730000dge0404000", "411790009fe30230000", "2s63103400000010"); // PsiPair<0.15exp(-0.065chi2)
    cuts.AddCutPCMCalo("000fc103", "0dm00009f9730000dge0404000", "411790009fe30230000", "2s63103400000010"); // PsiPair<0.18exp(-0.055chi2)
    cuts.AddCutPCMCalo("000fc103", "0dm00009f9730000dge0404000", "411790009fe30230000", "2s63103400000010"); // PsiPair<0.20exp(-0.050chi2)




  // ------ EJ1 cut variations
  } else if (trainConfig == 160) { // 
    cuts.AddCutPCMCalo("000fb103", "0dm00009f9730000dge0404000", "411790009fe30230000", "2s63103400000010"); // NL var. etc which is handled in correction framework
  } else if (trainConfig == 161) { // background variation
    cuts.AddCutPCMCalo("000fb103", "0dm00009f9730000dge0404000", "411790009fe30230000", "2163103400000010"); // jet mixing back
  } else if (trainConfig == 162) { // alpha cut variation
    cuts.AddCutPCMCalo("000fb103", "0dm00009f9730000dge0404000", "411790009fe30230000", "2s63105400000010"); // alpha cut 0-0.75
    cuts.AddCutPCMCalo("000fb103", "0dm00009f9730000dge0404000", "411790009fe30230000", "2s63108400000010"); // alpha cut 0-0.65
  } else if (trainConfig == 163) { // opening angle var.
    cuts.AddCutPCMCalo("000fb103", "0dm00009f9730000dge0404000", "411790009fe30230000", "2s63103400000000"); // no opening angle cut


  } else if (trainConfig == 165) { // M02 variation (std M02 = 0.5)
    cuts.AddCutPCMCalo("000fb103", "0dm00009f9730000dge0404000", "411790009fe30240000", "2s63103400000010"); // M02 = 0.4
    cuts.AddCutPCMCalo("000fb103", "0dm00009f9730000dge0404000", "411790009fe30220000", "2s63103400000010"); // M02 = 0.7
    cuts.AddCutPCMCalo("000fb103", "0dm00009f9730000dge0404000", "411790009fe30210000", "2s63103400000010"); // M02 = 1.0
    cuts.AddCutPCMCalo("000fb103", "0dm00009f9730000dge0404000", "411790009fe302v0000", "2s63103400000010"); // 0.5 < M02 < 0.7
  } else if (trainConfig == 166) { // TM variations for mesons
    cuts.AddCutPCMCalo("000fb103", "0dm00009f9730000dge0404000", "4117900090e30230000", "2s63103400000010"); // no TM
    cuts.AddCutPCMCalo("000fb103", "0dm00009f9730000dge0404000", "411790009ee30230000", "2s63103400000010"); // TM var EoverP 2.00
    cuts.AddCutPCMCalo("000fb103", "0dm00009f9730000dge0404000", "411790009ge30230000", "2s63103400000010"); // TM var EoverP 1.5
    cuts.AddCutPCMCalo("000fb103", "0dm00009f9730000dge0404000", "4117900097e30230000", "2s63103400000010"); // No E/p
    cuts.AddCutPCMCalo("000fb103", "0dm00009f9730000dge0404000", "411790009le30230000", "2s63103400000010"); // std + sec TM
    cuts.AddCutPCMCalo("000fb103", "0dm00009f9730000dge0404000", "411790009ne30230000", "2s63103400000010"); // TM var, Eta (0.035, 0.010, 2.5); Phi (0.085, 0.015, 2.)
  } else if (trainConfig == 167) { // cluster time
    cuts.AddCutPCMCalo("000fb103", "0dm00009f9730000dge0404000", "411790005fe30230000", "2s63103400000010"); // -50 - 50 ns
    cuts.AddCutPCMCalo("000fb103", "0dm00009f9730000dge0404000", "411790006fe30230000", "2s63103400000010"); // -30 - 35 ns
    cuts.AddCutPCMCalo("000fb103", "0dm00009f9730000dge0404000", "41179000afe30230000", "2s63103400000010"); // -12.5 - 13 ns
  } else if (trainConfig == 168) { // NCell
    cuts.AddCutPCMCalo("000fb103", "0dm00009f9730000dge0404000", "411790009fe3r230000", "2s63103400000010"); // EDC pi0 tagging for gamma clus, Gaussian Fit
    cuts.AddCutPCMCalo("000fb103", "0dm00009f9730000dge0404000", "411790009fe3n230000", "2s63103400000010"); // PCMEDC pi0 tagging for gamma clus, Gaussian Fit
    cuts.AddCutPCMCalo("000fb103", "0dm00009f9730000dge0404000", "411790009fe3m230000", "2s63103400000010"); // PCMEDC pi0 tagging for all clus, Gaussian Fit
    cuts.AddCutPCMCalo("000fb103", "0dm00009f9730000dge0404000", "411790009fe3l230000", "2s63103400000010"); // PCMEDC pi0 tagging for gamma clus, pol2 fit
  } else if (trainConfig == 169) { // min energy
    cuts.AddCutPCMCalo("000fb103", "0dm00009f9730000dge0404000", "411790009fe10230000", "2s63103400000010"); // INT7, minE = 0.5
    cuts.AddCutPCMCalo("000fb103", "0dm00009f9730000dge0404000", "411790009fe20230000", "2s63103400000010"); // INT7, minE = 0.6
    cuts.AddCutPCMCalo("000fb103", "0dm00009f9730000dge0404000", "411790009fe40230000", "2s63103400000010"); // INT7, minE = 0.8
  } else if (trainConfig == 170) { // Exotics
    cuts.AddCutPCMCalo("000fb103", "0dm00009f9730000dge0404000", "411790009f030230000", "2s63103400000010"); // no exotics
    cuts.AddCutPCMCalo("000fb103", "0dm00009f9730000dge0404000", "411790009fb30230000", "2s63103400000010"); // F+ < 0.95

  } else if (trainConfig == 171) { // min pt electron variation
    cuts.AddCutPCMCalo("000fb103", "0dm00069f9730000dge0404000", "411790009fe30230000", "2s63103400000010"); // min pT 40 MeV
    cuts.AddCutPCMCalo("000fb103", "0dm00049f9730000dge0404000", "411790009fe30230000", "2s63103400000010"); // min pT 50 MeV
    cuts.AddCutPCMCalo("000fb103", "0dm00019f9730000dge0404000", "411790009fe30230000", "2s63103400000010"); // min pT 100 MeV
  } else if (trainConfig == 172) { // min pt electron variation
    cuts.AddCutPCMCalo("000fb103", "0dm00008f9730000dge0404000", "411790009fe30230000", "2s63103400000010"); // TPC cluster 35%
    cuts.AddCutPCMCalo("000fb103", "0dm00006f9730000dge0404000", "411790009fe30230000", "2s63103400000010"); // TPC cluster 70%
    cuts.AddCutPCMCalo("000fb103", "0dm00009f9730000dge0604000", "411790009fe30230000", "2s63103400000010"); // cosPA 0.9
    cuts.AddCutPCMCalo("000fb103", "0dm00009f9730000dge0304000", "411790009fe30230000", "2s63103400000010"); // cosPA 0.75
  } else if (trainConfig == 173) {
    cuts.AddCutPCMCalo("000fb103", "0dm0000939730000dge0404000", "411790009fe30230000", "2s63103400000010"); // nsig electron   -4,5
    cuts.AddCutPCMCalo("000fb103", "0dm0000969730000dge0404000", "411790009fe30230000", "2s63103400000010"); // nsig electron -2.5,4
    cuts.AddCutPCMCalo("000fb103", "0dm00009f5730000dge0404000", "411790009fe30230000", "2s63103400000010"); // nsig pion 2,-10
    cuts.AddCutPCMCalo("000fb103", "0dm00009f1730000dge0404000", "411790009fe30230000", "2s63103400000010"); // nsig pion 0,-10
  } else if (trainConfig == 174) {
    cuts.AddCutPCMCalo("000fb103", "0dm00009f9030000dge0404000", "411790009fe30230000", "2s63103400000010"); // pion nsig min mom 0.50 GeV/c
    cuts.AddCutPCMCalo("000fb103", "0dm00009f9630000dge0404000", "411790009fe30230000", "2s63103400000010"); // pion nsig min mom 0.25 GeV/c
    cuts.AddCutPCMCalo("000fb103", "0dm00009f9760000dge0404000", "411790009fe30230000", "2s63103400000010"); // pion nsig max mom 2.00 GeV/c
    cuts.AddCutPCMCalo("000fb103", "0dm00009f9710000dge0404000", "411790009fe30230000", "2s63103400000010"); // pion nsig max mom 5.00 GeV/c
  } else if (trainConfig == 175) {
    cuts.AddCutPCMCalo("000fb103", "0dm00009f97300008ge0404000", "411790009fe30230000", "2s63103400000010"); // qT max 0.05 1D
    cuts.AddCutPCMCalo("000fb103", "0dm00009f97300003ge0404000", "411790009fe30230000", "2s63103400000010"); // qT max 0.05 1D
    cuts.AddCutPCMCalo("000fb103", "0dm00009f97300002ge0404000", "411790009fe30230000", "2s63103400000010"); // qT max 0.06 1D
    cuts.AddCutPCMCalo("000fb103", "0dm00009f97300009ge0404000", "411790009fe30230000", "2s63103400000010"); // qT max 0.03 1D
  } else if (trainConfig == 176) {
    cuts.AddCutPCMCalo("000fb103", "0dm00009f9730000dg50404000", "411790009fe30230000", "2s63103400000010"); // Psi pair 0.1  1D
    cuts.AddCutPCMCalo("000fb103", "0dm00009f9730000dg10404000", "411790009fe30230000", "2s63103400000010"); // Psi pair 0.1  1D
    cuts.AddCutPCMCalo("000fb103", "0dm00009f9730000dg60404000", "411790009fe30230000", "2s63103400000010"); // Psi pair 0.05  1D
    cuts.AddCutPCMCalo("000fb103", "0dm00009f9730000dg80404000", "411790009fe30230000", "2s63103400000010"); // Psi pair 0.2  1D
  } else if (trainConfig == 177) {
    cuts.AddCutPCMCalo("000fb103", "0dm00009f9730000c259404000", "411790009fe30230000", "2s63103400000010"); // qT<0.110pT (2D) alpha<0.99
    cuts.AddCutPCMCalo("000fb103", "0dm00009f9730000a259404000", "411790009fe30230000", "2s63103400000010"); // qT<0.125pT (2D) alpha<0.99
    cuts.AddCutPCMCalo("000fb103", "0dm00009f9730000e259404000", "411790009fe30230000", "2s63103400000010"); // qT<0.130pT (2D) alpha<0.99
  } else if (trainConfig == 178) {
    cuts.AddCutPCMCalo("000fb103", "0dm00009f9730000dge0404000", "411790009fe30230000", "2s63103400000010"); // PsiPair<0.15exp(-0.065chi2)
    cuts.AddCutPCMCalo("000fb103", "0dm00009f9730000dge0404000", "411790009fe30230000", "2s63103400000010"); // PsiPair<0.18exp(-0.055chi2)
    cuts.AddCutPCMCalo("000fb103", "0dm00009f9730000dge0404000", "411790009fe30230000", "2s63103400000010"); // PsiPair<0.20exp(-0.050chi2)


    //---------------------------------------
    // configs for eta meson pp 13 TeV
    //---------------------------------------
  } else if (trainConfig == 1002) {
    cuts.AddCutPCMCalo("00010103", "0dm00009f9730000dge0404000", "411790009fe30230000", "2r63103l00000010"); // in-Jet, mass cut eta: 0.5-0.6, rotation back
  } else if (trainConfig == 1003) {
    cuts.AddCutPCMCalo("00010103", "0dm00009f9730000dge0404000", "411790009fe30230000", "2163103l00000010"); // in-Jet, mass cut eta: 0.5-0.6, mixed jet back
  } else {
    Error(Form("MesonJetCorrelation_ConvCalo_%i", trainConfig), "wrong trainConfig variable no cuts have been specified for the configuration");
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
  TList* ClusterCutList = new TList();
  TList* MesonCutList = new TList();

  EventCutList->SetOwner(kTRUE);
  AliConvEventCuts** analysisEventCuts = new AliConvEventCuts*[numberOfCuts];
  ConvCutList->SetOwner(kTRUE);
  AliConversionPhotonCuts** analysisConvCuts = new AliConversionPhotonCuts*[numberOfCuts];
  ClusterCutList->SetOwner(kTRUE);
  AliCaloPhotonCuts** analysisClusterCuts = new AliCaloPhotonCuts*[numberOfCuts];
  MesonCutList->SetOwner(kTRUE);
  AliConversionMesonCuts** analysisMesonCuts = new AliConversionMesonCuts*[numberOfCuts];

  for (Int_t i = 0; i < numberOfCuts; i++) {
    //create AliCaloTrackMatcher instance, if there is none present
    TString caloCutPos = cuts.GetClusterCut(i);
    caloCutPos.Resize(1);
    TString TrackMatcherName = Form("CaloTrackMatcher_%s_%i", caloCutPos.Data(), trackMatcherRunningMode);
    if (corrTaskSetting.CompareTo("")) {
      TrackMatcherName = TrackMatcherName + "_" + corrTaskSetting.Data();
      cout << "Using separate track matcher for correction framework setting: " << TrackMatcherName.Data() << endl;
    }
    if (!(AliCaloTrackMatcher*)mgr->GetTask(TrackMatcherName.Data())) {
      AliCaloTrackMatcher* fTrackMatcher = new AliCaloTrackMatcher(TrackMatcherName.Data(), caloCutPos.Atoi(), trackMatcherRunningMode);
      fTrackMatcher->SetV0ReaderName(V0ReaderName);
      fTrackMatcher->SetCorrectionTaskSetting(corrTaskSetting);
      mgr->AddTask(fTrackMatcher);
      mgr->ConnectInput(fTrackMatcher, 0, cinput);
    }

    //---------------------------------------------------------//
    //------------------------ Event Cuts ---------------------//
    //---------------------------------------------------------//
    analysisEventCuts[i] = new AliConvEventCuts();
    // analysisEventCuts[i]->SetCaloTriggerHelperName(TriggerHelperName.Data());
    analysisEventCuts[i]->SetTriggerMimicking(enableTriggerMimicking);
    analysisEventCuts[i]->SetTriggerOverlapRejecion(enableTriggerOverlapRej);
    analysisEventCuts[i]->SetV0ReaderName(V0ReaderName);
    analysisEventCuts[i]->SetCorrectionTaskSetting(corrTaskSetting);
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
    //------------------- Calo Photon Cuts --------------------//
    //---------------------------------------------------------//
    analysisClusterCuts[i] = new AliCaloPhotonCuts(isMC);
    // analysisClusterCuts[i]->SetHistoToModifyAcceptance(histoAcc);
    analysisClusterCuts[i]->SetV0ReaderName(V0ReaderName);
    analysisClusterCuts[i]->SetCorrectionTaskSetting(corrTaskSetting);
    analysisClusterCuts[i]->SetCaloTrackMatcherName(TrackMatcherName);
    if (enableLightOutput > 0)
      analysisClusterCuts[i]->SetLightOutput(kTRUE);
    analysisClusterCuts[i]->InitializeCutsFromCutString((cuts.GetClusterCut(i)).Data());
    analysisClusterCuts[i]->SetExtendedMatchAndQA(enableExtMatchAndQA);
    ClusterCutList->Add(analysisClusterCuts[i]);

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
    analysisMesonCuts[i]->SetIsMergedClusterCut(2);
    analysisMesonCuts[i]->SetCaloMesonCutsObject(analysisClusterCuts[i]);
    analysisMesonCuts[i]->SetFillCutHistograms("");
    // analysisEventCuts[i]->SetAcceptedHeader(HeaderList);
    analysisClusterCuts[i]->SetFillCutHistograms("");
    if (enableLightOutput > 0)
      analysisMesonCuts[i]->SetLightOutput(kTRUE);
    if (analysisMesonCuts[i]->DoGammaSwappForBg())
      analysisClusterCuts[i]->SetUseEtaPhiMapForBackCand(kTRUE); // needed in case of rotation background QA histos
    MesonCutList->Add(analysisMesonCuts[i]);
  }

  task->SetMesonKind(meson);
  task->SetMesonZPt(runOnlyZPt);
  task->SetIsConvCalo(true);
  task->SetJetContainerAddName(nameJetFinder);
  task->SetEventCutList(numberOfCuts, EventCutList);
  task->SetCaloCutList(numberOfCuts, ClusterCutList);
  task->SetMesonCutList(numberOfCuts, MesonCutList);
  task->SetConversionCutList(numberOfCuts, ConvCutList);
  if(enableMatBudWeightsPi0) task->SetDoMaterialBudgetWeightingOfGammasForTrueMesons(true);
  task->SetCorrectionTaskSetting(corrTaskSetting);
  task->SetDoMesonQA(enableQAMesonTask); //Attention new switch for Pi0 QA
                                         //   task->SetDoClusterQA(enableQAClusterTask);  //Attention new switch small for Cluster QA
  task->SetUseTHnSparseForResponse(enableTHnSparse);
  task->SetDoUseCentralEvtSelection(useCentralEvtSelection);
  task->SetForcePi0Unstable(setPi0Unstable);
  task->SetUseMixedBackAdd(enableAddBackground);
  task->SetDoRadiusDependence(enableRadiusDep);

  //connect containers
  TString nameContainer = Form("MesonJetCorrelation_ConvCalo_%i_%i%s%s", meson, trainConfig, corrTaskSetting.EqualTo("") == true ? "" : Form("_%s", corrTaskSetting.Data()), nameJetFinder.EqualTo("") == true ? "" : Form("_%s", nameJetFinder.Data()) );
  AliAnalysisDataContainer* coutput = mgr->CreateContainer(nameContainer, TList::Class(), AliAnalysisManager::kOutputContainer, Form("MesonJetCorrelation_ConvCalo_%i_%i.root", meson, trainConfig));

  mgr->AddTask(task);
  mgr->ConnectInput(task, 0, cinput);
  mgr->ConnectOutput(task, 1, coutput);

  return;
}

