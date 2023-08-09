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

void AddTask_MesonJetCorr_Calo(
  Int_t trainConfig = 1,                // change different set of cuts
  Int_t isMC = 0,                       // run MC
  int meson = 0,                        // meson: 0=pi0, 1 = eta
  TString photonCutNumberV0Reader = "", // 00000008400000000100000000 nom. B, 00000088400000000100000000 low B
  TString periodNameV0Reader = "",
  // general setting for task
  Int_t enableQAMesonTask = 0,             // enable QA in AliAnalysisTaskGammaConvV1
  Int_t enableQAClusterTask = 0,           // enable additional QA task
  Int_t enableExtMatchAndQA = 0,           // disabled (0), extMatch (1), extQA_noCellQA (2), extMatch+extQA_noCellQA (3), extQA+cellQA (4), extMatch+extQA+cellQA (5)
  Int_t enableLightOutput = kFALSE,        // switch to run light output (only essential histograms for afterburner) (pi0 only mode: lighOutput + 10)
  Bool_t enableTHnSparse = kFALSE,         // switch on THNsparse
  Int_t enableTriggerMimicking = 0,        // enable trigger mimicking
  Bool_t enableTriggerOverlapRej = kFALSE, // enable trigger overlap rejection
  TString settingMaxFacPtHard = "3.",      // maximum factor between hardest jet and ptHard generated
  Int_t debugLevel = 0,                    // introducing debug levels for grid running
  // settings for weights
  // FPTW:fileNamePtWeights, FMUW:fileNameMultWeights, separate with ;
  TString fileNameExternalInputs = "",
  Int_t doWeightingPart = 0,                   // enable Weighting
  TString generatorName = "DPMJET",            // generator Name
  Bool_t enableMultiplicityWeighting = kFALSE, //
  TString periodNameAnchor = "",               //
  // special settings
  Bool_t enableSortingMCLabels = kTRUE, // enable sorting for MC cluster labels
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

  TString addTaskName = "AddTask_MesonJetCorr_Calo";
  TString sAdditionalTrainConfig = cuts.GetSpecialSettingFromAddConfig(additionalTrainConfig, "", "", addTaskName);
  if (sAdditionalTrainConfig.Atoi() > 0) {
    trainConfig = trainConfig + sAdditionalTrainConfig.Atoi();
    cout << "INFO: " << addTaskName.Data() << " running additionalTrainConfig '" << sAdditionalTrainConfig.Atoi() << "', train config: '" << trainConfig << "'" << endl;
  }

  TString fileNamePtWeights = cuts.GetSpecialFileNameFromString(fileNameExternalInputs, "FPTW:");
  TString fileNameMultWeights = cuts.GetSpecialFileNameFromString(fileNameExternalInputs, "FMUW:");
  TString fileNameCustomTriggerMimicOADB = cuts.GetSpecialFileNameFromString(fileNameExternalInputs, "FTRM:");

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
  task = new AliAnalysisTaskMesonJetCorrelation(Form("MesonJetCorrelation_Calo_%i_%i", meson, trainConfig));
  task->SetIsHeavyIon(isHeavyIon);
  task->SetIsMC(isMC);
  task->SetV0ReaderName(V0ReaderName);
  task->SetCorrectionTaskSetting(corrTaskSetting);
  task->SetTrackMatcherRunningMode(trackMatcherRunningMode); // have to do this!

  //---------------------------------------
  // configs for pi0 meson pp 13 TeV
  //---------------------------------------
  // configs without NonLinearity as NonLin is applied in CF
  if (trainConfig == 1) {
    cuts.AddCutCalo("00010103", "411790009fe30230000", "0s631031000000d0"); // test config without in-jet selection
  } else if (trainConfig == 2) {
    cuts.AddCutCalo("00010103", "411790009fe30230000", "2s631034000000d0"); // in-jet, pi0 mass: 0.1-0.15, rotation back
  } else if (trainConfig == 3) {
    cuts.AddCutCalo("00010103", "411790009fe30230000", "21631034000000d0"); // in-jet, pi0 mass: 0.1-0.15, mixed jet back
  } else if (trainConfig == 4) {
    cuts.AddCutCalo("0008e103", "411790009fe30230000", "2s631034000000d0"); // EG2 in-jet, pi0 mass: 0.1-0.15, rotation back
  } else if (trainConfig == 5) {
    cuts.AddCutCalo("0008d103", "411790009fe30230000", "2s631034000000d0"); // EG1 in-jet, pi0 mass: 0.1-0.15, rotation back
  } else if (trainConfig == 6) {
    cuts.AddCutCalo("000fc103", "411790009fe30230000", "2s631034000000d0"); // Jet-low trigg in-jet, pi0 mass: 0.1-0.15, rotation back
  } else if (trainConfig == 7) {
    cuts.AddCutCalo("000fb103", "411790009fe30230000", "2s631034000000d0"); // Jet-high trigg in-jet, pi0 mass: 0.1-0.15, rotation back

  // EJ1 and EJ2 only EMCal triggers
  } else if (trainConfig == 8) {
    cuts.AddCutCalo("000f5103", "411790009fe30230000", "2s631034000000d0"); // Jet-low trigg in-jet, pi0 mass: 0.1-0.15, rotation back
  } else if (trainConfig == 9) {
    cuts.AddCutCalo("000f3103", "411790009fe30230000", "2s631034000000d0"); // Jet-high trigg in-jet, pi0 mass: 0.1-0.15, rotation back

  // configs with Mesons only in EMCal, EMCal triggers only
  } else if (trainConfig == 10) {
    cuts.AddCutCalo("00010103", "111110009fe30230000", "2s631034000000d0"); // in-jet, pi0 mass: 0.1-0.15, rotation back
  } else if (trainConfig == 11) {
    cuts.AddCutCalo("000f5103", "111110009fe30230000", "2s631034000000d0"); // Jet-low trigg in-jet, pi0 mass: 0.1-0.15, rotation back
  } else if (trainConfig == 12) {
    cuts.AddCutCalo("000f3103", "111110009fe30230000", "2s631034000000d0"); // Jet-high trigg in-jet, pi0 mass: 0.1-0.15, rotation back


  } else if (trainConfig == 14) { // same as 4 but with mixed jet back
    cuts.AddCutCalo("0008e103", "411790009fe30230000", "21631034000000d0"); // EG2 in-jet, pi0 mass: 0.1-0.15, mixed jet back
  } else if (trainConfig == 15) {// same as 5 but with mixed jet back
    cuts.AddCutCalo("0008d103", "411790009fe30230000", "21631034000000d0"); // EG1 in-jet, pi0 mass: 0.1-0.15, mixed jet back
  } else if (trainConfig == 16) {// same as 6 but with mixed jet back
    cuts.AddCutCalo("000fc103", "411790009fe30230000", "21631034000000d0"); // Jet-low trigg in-jet, pi0 mass: 0.1-0.15, mixed jet back
  } else if (trainConfig == 17) {// same as 7 but with mixed jet back
    cuts.AddCutCalo("000fb103", "411790009fe30230000", "21631034000000d0"); // Jet-high trigg in-jet, pi0 mass: 0.1-0.15, mixed jet back

  } else if (trainConfig == 20) {
    cuts.AddCutCalo("00010103", "411790009fe30230000", "es631034000000d0"); // decay daughters also inside jet
  } else if (trainConfig == 21) {
    cuts.AddCutCalo("000fc103", "411790009fe30230000", "es631034000000d0"); // decay daughters also inside jet
    cuts.AddCutCalo("000fb103", "411790009fe30230000", "es631034000000d0"); // decay daughters also inside jet

  // configs with NonLinearity 
  } else if (trainConfig == 22) {
    cuts.AddCutCalo("00010103", "411790109fe30230000", "2s631034000000d0"); // in-jet, pi0 mass: 0.1-0.15, rotation back
  } else if (trainConfig == 23) {
    cuts.AddCutCalo("00010103", "411790109fe30230000", "21631034000000d0"); // in-jet, pi0 mass: 0.1-0.15, mixed jet back
  } else if (trainConfig == 24) {
    cuts.AddCutCalo("0008e103", "411790109fe30230000", "2s631034000000d0"); // EG2 in-jet, pi0 mass: 0.1-0.15, rotation back
  } else if (trainConfig == 25) {
    cuts.AddCutCalo("0008d103", "411790109fe30230000", "2s631034000000d0"); // EG1 in-jet, pi0 mass: 0.1-0.15, rotation back
  } else if (trainConfig == 26) {
    cuts.AddCutCalo("000fc103", "411790109fe30230000", "2s631034000000d0"); // Jet-low trigg in-jet, pi0 mass: 0.1-0.15, rotation back
  } else if (trainConfig == 27) {
    cuts.AddCutCalo("000fb103", "411790109fe30230000", "2s631034000000d0"); // Jet-high trigg in-jet, pi0 mass: 0.1-0.15, rotation back


  // configs with eta < 0.5
  } else if (trainConfig == 30) {
    cuts.AddCutCalo("00010103", "411790009fe30230000", "2s634034000000d0"); // in-jet, pi0 mass: 0.1-0.15, rotation back
  } else if (trainConfig == 31) {
    cuts.AddCutCalo("000fc103", "411790009fe30230000", "2s634034000000d0"); // Jet-low trigg in-jet, pi0 mass: 0.1-0.15, rotation back
  } else if (trainConfig == 32) {
    cuts.AddCutCalo("000fb103", "411790009fe30230000", "2s634034000000d0"); // Jet-high trigg in-jet, pi0 mass: 0.1-0.15, rotation back

  // configs with eta < 0.5, only EMCal triggers (not DCal)
  } else if (trainConfig == 33) {
    cuts.AddCutCalo("000f5103", "411790009fe30230000", "2s634034000000d0"); // Jet-low trigg in-jet, pi0 mass: 0.1-0.15, rotation back
  } else if (trainConfig == 34) {
    cuts.AddCutCalo("000f3103", "411790009fe30230000", "2s634034000000d0"); // Jet-high trigg in-jet, pi0 mass: 0.1-0.15, rotation back

  // configs with eta < 0.5, only EMCal triggers (not DCal), Mesons only with EMCal
  } else if (trainConfig == 35) {
    cuts.AddCutCalo("000f5103", "111110009fe30230000", "2s634034000000d0"); // Jet-low trigg in-jet, pi0 mass: 0.1-0.15, rotation back
  } else if (trainConfig == 36) {
    cuts.AddCutCalo("000f3103", "111110009fe30230000", "2s634034000000d0"); // Jet-high trigg in-jet, pi0 mass: 0.1-0.15, rotation back

  // configs with eta < 1.35
  } else if (trainConfig == 40) {
    cuts.AddCutCalo("00010103", "411790009fe30230000", "2s630034000000d0"); // in-jet, pi0 mass: 0.1-0.15, rotation back
  } else if (trainConfig == 41) {
    cuts.AddCutCalo("000fc103", "411790009fe30230000", "2s630034000000d0"); // Jet-low trigg in-jet, pi0 mass: 0.1-0.15, rotation back
  } else if (trainConfig == 42) {
    cuts.AddCutCalo("000fb103", "411790009fe30230000", "2s630034000000d0"); // Jet-high trigg in-jet, pi0 mass: 0.1-0.15, rotation back

  //---------------------------------------
  // Cut variations for standard cut (2)
  //---------------------------------------

  //--------  INT7 Meson cut variations
  } else if (trainConfig == 100) { // 
    cuts.AddCutCalo("00010103", "411790009fe30230000", "2s631034000000d0"); // NL var. etc which is handled in correction framework
  } else if (trainConfig == 101) { // background variation
    cuts.AddCutCalo("00010103", "411790009fe30230000", "21631034000000d0"); // jet mixing back
  } else if (trainConfig == 102) { // alpha cut variation
    cuts.AddCutCalo("00010103", "411790009fe30230000", "2s631054000000d0"); // alpha cut 0-0.75
    cuts.AddCutCalo("00010103", "411790009fe30230000", "2s631084000000d0"); // alpha cut 0-0.65
  } else if (trainConfig == 103) { // opening angle var.
    cuts.AddCutCalo("00010103", "411790009fe30230000", "2s631034000000b0"); // INT7 Op. Ang. var 1 cell dist + 0.0152
    cuts.AddCutCalo("00010103", "411790009fe30230000", "2s631034000000g0"); // INT7 Op. Ang. var 1 cell dist + 0.0202
    cuts.AddCutCalo("00010103", "411790009fe30230000", "2s631034000000a0"); // INT7 Op. Ang. var 1 cell dist + 0


  } else if (trainConfig == 105) { // M02 variation (std M02 = 0.5)
    cuts.AddCutCalo("00010103", "411790009fe30240000", "2s631034000000d0"); // M02 = 0.4
    cuts.AddCutCalo("00010103", "411790009fe30220000", "2s631034000000d0"); // M02 = 0.7
    cuts.AddCutCalo("00010103", "411790009fe30210000", "2s631034000000d0"); // M02 = 1.0
    cuts.AddCutCalo("00010103", "411790009fe302v0000", "2s631034000000d0"); // 0.5 < M02 < 0.7
  } else if (trainConfig == 106) { // TM variations for mesons
    cuts.AddCutCalo("00010103", "4117900090e30230000", "2s631034000000d0"); // INT7 no TM
    cuts.AddCutCalo("00010103", "411790009ee30230000", "2s631034000000d0"); // TM var EoverP 2.00
    cuts.AddCutCalo("00010103", "411790009ge30230000", "2s631034000000d0"); // TM var EoverP 1.5
    cuts.AddCutCalo("00010103", "4117900097e30230000", "2s631034000000d0"); // No E/p
    cuts.AddCutCalo("00010103", "411790009le30230000", "2s631034000000d0"); // std + sec TM
    cuts.AddCutCalo("00010103", "411790009ne30230000", "2s631034000000d0"); // INT7 TM var, Eta (0.035, 0.010, 2.5); Phi (0.085, 0.015, 2.)
  } else if (trainConfig == 107) { // cluster time
    cuts.AddCutCalo("00010103", "411790005fe30230000", "2s631034000000d0"); // -50 - 50 ns
    cuts.AddCutCalo("00010103", "411790006fe30230000", "2s631034000000d0"); // -30 - 35 ns
    cuts.AddCutCalo("00010103", "41179000afe30230000", "2s631034000000d0"); // -12.5 - 13 ns
  } else if (trainConfig == 108) { // NCell
    cuts.AddCutCalo("00010103", "411790009fe3r230000", "2s631034000000d0"); // INT7 EDC pi0 tagging for gamma clus, Gaussian Fit
    cuts.AddCutCalo("00010103", "411790009fe3n230000", "2s631034000000d0"); // INT7 PCMEDC pi0 tagging for gamma clus, Gaussian Fit
    cuts.AddCutCalo("00010103", "411790009fe3m230000", "2s631034000000d0"); // INT7 PCMEDC pi0 tagging for all clus, Gaussian Fit
    cuts.AddCutCalo("00010103", "411790009fe3l230000", "2s631034000000d0"); // INT7 PCMEDC pi0 tagging for gamma clus, pol2 fit
  } else if (trainConfig == 109) { // min energy
    cuts.AddCutCalo("00010103", "411790009fe10230000", "2s631034000000d0"); // INT7, minE = 0.5
    cuts.AddCutCalo("00010103", "411790009fe20230000", "2s631034000000d0"); // INT7, minE = 0.6
    cuts.AddCutCalo("00010103", "411790009fe40230000", "2s631034000000d0"); // INT7, minE = 0.8
  } else if (trainConfig == 110) { // Exotics
    cuts.AddCutCalo("00010103", "411790009f030230000", "2s631034000000d0"); // no exotics
    cuts.AddCutCalo("00010103", "411790009fb30230000", "2s631034000000d0"); // F+ < 0.95



    //--------  EJ2 Meson cut variations
  } else if (trainConfig == 130) { // 
    cuts.AddCutCalo("000fc103", "411790009fe30230000", "2s631034000000d0"); // NL var. etc which is handled in correction framework
  } else if (trainConfig == 131) { // background variation
    cuts.AddCutCalo("000fc103", "411790009fe30230000", "21631034000000d0"); // jet mixing back
  } else if (trainConfig == 132) { // alpha cut variation
    cuts.AddCutCalo("000fc103", "411790009fe30230000", "2s631054000000d0"); // alpha cut 0-0.75
    cuts.AddCutCalo("000fc103", "411790009fe30230000", "2s631084000000d0"); // alpha cut 0-0.65
  } else if (trainConfig == 133) { // opening angle var.
    cuts.AddCutCalo("000fc103", "411790009fe30230000", "2s631034000000b0"); // INT7 Op. Ang. var 1 cell dist + 0.0152
    cuts.AddCutCalo("000fc103", "411790009fe30230000", "2s631034000000g0"); // INT7 Op. Ang. var 1 cell dist + 0.0202
    cuts.AddCutCalo("000fc103", "411790009fe30230000", "2s631034000000a0"); // INT7 Op. Ang. var 1 cell dist + 0


  } else if (trainConfig == 135) { // M02 variation (std M02 = 0.5)
    cuts.AddCutCalo("000fc103", "411790009fe30240000", "2s631034000000d0"); // M02 = 0.4
    cuts.AddCutCalo("000fc103", "411790009fe30220000", "2s631034000000d0"); // M02 = 0.7
    cuts.AddCutCalo("000fc103", "411790009fe30210000", "2s631034000000d0"); // M02 = 1.0
    cuts.AddCutCalo("000fc103", "411790009fe302v0000", "2s631034000000d0"); // 0.5 < M02 < 0.7
  } else if (trainConfig == 136) { // TM variations for mesons
    cuts.AddCutCalo("000fc103", "4117900090e30230000", "2s631034000000d0"); // INT7 no TM
    cuts.AddCutCalo("000fc103", "411790009ee30230000", "2s631034000000d0"); // TM var EoverP 2.00
    cuts.AddCutCalo("000fc103", "411790009ge30230000", "2s631034000000d0"); // TM var EoverP 1.5
    cuts.AddCutCalo("000fc103", "4117900097e30230000", "2s631034000000d0"); // No E/p
    cuts.AddCutCalo("000fc103", "411790009le30230000", "2s631034000000d0"); // std + sec TM
    cuts.AddCutCalo("000fc103", "411790009ne30230000", "2s631034000000d0"); // INT7 TM var, Eta (0.035, 0.010, 2.5); Phi (0.085, 0.015, 2.)
  } else if (trainConfig == 137) { // cluster time
    cuts.AddCutCalo("000fc103", "411790005fe30230000", "2s631034000000d0"); // -50 - 50 ns
    cuts.AddCutCalo("000fc103", "411790006fe30230000", "2s631034000000d0"); // -30 - 35 ns
    cuts.AddCutCalo("000fc103", "41179000afe30230000", "2s631034000000d0"); // -12.5 - 13 ns
  } else if (trainConfig == 138) { // NCell
    cuts.AddCutCalo("000fc103", "411790009fe3r230000", "2s631034000000d0"); // INT7 EDC pi0 tagging for gamma clus, Gaussian Fit
    cuts.AddCutCalo("000fc103", "411790009fe3n230000", "2s631034000000d0"); // INT7 PCMEDC pi0 tagging for gamma clus, Gaussian Fit
    cuts.AddCutCalo("000fc103", "411790009fe3m230000", "2s631034000000d0"); // INT7 PCMEDC pi0 tagging for all clus, Gaussian Fit
    cuts.AddCutCalo("000fc103", "411790009fe3l230000", "2s631034000000d0"); // INT7 PCMEDC pi0 tagging for gamma clus, pol2 fit
  } else if (trainConfig == 139) { // min energy
    cuts.AddCutCalo("000fc103", "411790009fe10230000", "2s631034000000d0"); // INT7, minE = 0.5
    cuts.AddCutCalo("000fc103", "411790009fe20230000", "2s631034000000d0"); // INT7, minE = 0.6
    cuts.AddCutCalo("000fc103", "411790009fe40230000", "2s631034000000d0"); // INT7, minE = 0.8
  } else if (trainConfig == 140) { // Exotics
    cuts.AddCutCalo("000fc103", "411790009f030230000", "2s631034000000d0"); // no exotics
    cuts.AddCutCalo("000fc103", "411790009fb30230000", "2s631034000000d0"); // F+ < 0.95

  //--------  EJ1 Meson cut variations
  } else if (trainConfig == 160) { // 
    cuts.AddCutCalo("000fb103", "411790009fe30230000", "2s631034000000d0"); // NL var. etc which is handled in correction framework
  } else if (trainConfig == 161) { // background variation
    cuts.AddCutCalo("000fb103", "411790009fe30230000", "21631034000000d0"); // jet mixing back
  } else if (trainConfig == 162) { // alpha cut variation
    cuts.AddCutCalo("000fb103", "411790009fe30230000", "2s631054000000d0"); // alpha cut 0-0.75
    cuts.AddCutCalo("000fb103", "411790009fe30230000", "2s631084000000d0"); // alpha cut 0-0.65
  } else if (trainConfig == 163) { // opening angle var.
    cuts.AddCutCalo("000fb103", "411790009fe30230000", "2s631034000000b0"); // INT7 Op. Ang. var 1 cell dist + 0.0152
    cuts.AddCutCalo("000fb103", "411790009fe30230000", "2s631034000000g0"); // INT7 Op. Ang. var 1 cell dist + 0.0202
    cuts.AddCutCalo("000fb103", "411790009fe30230000", "2s631034000000a0"); // INT7 Op. Ang. var 1 cell dist + 0


  } else if (trainConfig == 165) { // M02 variation (std M02 = 0.5)
    cuts.AddCutCalo("000fb103", "411790009fe30240000", "2s631034000000d0"); // M02 = 0.4
    cuts.AddCutCalo("000fb103", "411790009fe30220000", "2s631034000000d0"); // M02 = 0.7
    cuts.AddCutCalo("000fb103", "411790009fe30210000", "2s631034000000d0"); // M02 = 1.0
    cuts.AddCutCalo("000fb103", "411790009fe302v0000", "2s631034000000d0"); // 0.5 < M02 < 0.7
  } else if (trainConfig == 166) { // TM variations for mesons
    cuts.AddCutCalo("000fb103", "4117900090e30230000", "2s631034000000d0"); // INT7 no TM
    cuts.AddCutCalo("000fb103", "411790009ee30230000", "2s631034000000d0"); // TM var EoverP 2.00
    cuts.AddCutCalo("000fb103", "411790009ge30230000", "2s631034000000d0"); // TM var EoverP 1.5
    cuts.AddCutCalo("000fb103", "4117900097e30230000", "2s631034000000d0"); // No E/p
    cuts.AddCutCalo("000fb103", "411790009le30230000", "2s631034000000d0"); // std + sec TM
    cuts.AddCutCalo("000fb103", "411790009ne30230000", "2s631034000000d0"); // INT7 TM var, Eta (0.035, 0.010, 2.5); Phi (0.085, 0.015, 2.)
  } else if (trainConfig == 167) { // cluster time
    cuts.AddCutCalo("000fb103", "411790005fe30230000", "2s631034000000d0"); // -50 - 50 ns
    cuts.AddCutCalo("000fb103", "411790006fe30230000", "2s631034000000d0"); // -30 - 35 ns
    cuts.AddCutCalo("000fb103", "41179000afe30230000", "2s631034000000d0"); // -12.5 - 13 ns
  } else if (trainConfig == 168) { // NCell
    cuts.AddCutCalo("000fb103", "411790009fe3r230000", "2s631034000000d0"); // INT7 EDC pi0 tagging for gamma clus, Gaussian Fit
    cuts.AddCutCalo("000fb103", "411790009fe3n230000", "2s631034000000d0"); // INT7 PCMEDC pi0 tagging for gamma clus, Gaussian Fit
    cuts.AddCutCalo("000fb103", "411790009fe3m230000", "2s631034000000d0"); // INT7 PCMEDC pi0 tagging for all clus, Gaussian Fit
    cuts.AddCutCalo("000fb103", "411790009fe3l230000", "2s631034000000d0"); // INT7 PCMEDC pi0 tagging for gamma clus, pol2 fit
  } else if (trainConfig == 169) { // min energy
    cuts.AddCutCalo("000fb103", "411790009fe10230000", "2s631034000000d0"); // INT7, minE = 0.5
    cuts.AddCutCalo("000fb103", "411790009fe20230000", "2s631034000000d0"); // INT7, minE = 0.6
    cuts.AddCutCalo("000fb103", "411790009fe40230000", "2s631034000000d0"); // INT7, minE = 0.8
  } else if (trainConfig == 170) { // Exotics
    cuts.AddCutCalo("000fb103", "411790009f030230000", "2s631034000000d0"); // no exotics
    cuts.AddCutCalo("000fb103", "411790009fb30230000", "2s631034000000d0"); // F+ < 0.95

    //---------------------------------------
    // configs for eta meson pp 13 TeV
    //---------------------------------------
  } else if (trainConfig == 1002) {
    cuts.AddCutCalo("00010103", "411790109fe30230000", "2s631034000000d0"); // in-jet, pi0 mass: 0.1-0.15, rotation back
  } else if (trainConfig == 1003) {
    cuts.AddCutCalo("00010103", "411790109fe30230000", "21631034000000d0"); // in-jet, pi0 mass: 0.1-0.15, mixed jet back
  } else {
    Error(Form("MesonJetCorrelation_Calo_%i", trainConfig), "wrong trainConfig variable no cuts have been specified for the configuration");
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
  TList* ClusterCutList = new TList();
  TList* MesonCutList = new TList();

  EventCutList->SetOwner(kTRUE);
  AliConvEventCuts** analysisEventCuts = new AliConvEventCuts*[numberOfCuts];
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
  task->SetIsCalo(true);
  if(additionalTrainConfig.Contains("JET")){task->SetJetContainerAddName(nameJetFinder);}
  task->SetEventCutList(numberOfCuts, EventCutList);
  task->SetCaloCutList(numberOfCuts, ClusterCutList);
  task->SetMesonCutList(numberOfCuts, MesonCutList);
  //   task->SetDoMesonAnalysis(kTRUE); // I think we dont need that!
  task->SetCorrectionTaskSetting(corrTaskSetting);
  task->SetDoMesonQA(enableQAMesonTask); //Attention new switch for Pi0 QA
                                         //   task->SetDoClusterQA(enableQAClusterTask);  //Attention new switch small for Cluster QA
  task->SetUseTHnSparseForResponse(enableTHnSparse);
  task->SetDoUseCentralEvtSelection(useCentralEvtSelection);
  task->SetForcePi0Unstable(setPi0Unstable);
  task->SetUseMixedBackAdd(enableAddBackground);
  task->SetDoRadiusDependence(enableRadiusDep);

  //connect containers
  TString nameContainer = Form("MesonJetCorrelation_Calo_%i_%i%s%s", meson, trainConfig, corrTaskSetting.EqualTo("") == true ? "" : Form("_%s", corrTaskSetting.Data()), nameJetFinder.EqualTo("") == true ? "" : Form("_%s", nameJetFinder.Data()) );
  AliAnalysisDataContainer* coutput = mgr->CreateContainer(nameContainer, TList::Class(), AliAnalysisManager::kOutputContainer, Form("MesonJetCorrelation_Calo_%i_%i.root", meson, trainConfig));

  mgr->AddTask(task);
  mgr->ConnectInput(task, 0, cinput);
  mgr->ConnectOutput(task, 1, coutput);

  return;
}

