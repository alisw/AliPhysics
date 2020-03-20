/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: Friederike Bock, Annika Passfeld                               *
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
//($ALIPHYSICS/PWGGA/GammaConv/AliAnalysisTaskGammaConvV1.cxx) for
//pp together with all supporting classes
//***************************************************************************************

void AddTask_GammaConvV1_pPb(
    Int_t     trainConfig                   = 1,        // change different set of cuts
    Int_t     isMC                          = 0,        // run MC
    TString   photonCutNumberV0Reader       = "",
    TString   periodNameV0Reader            = "",
    // general setting for task
    Int_t     enableQAMesonTask             = 0,        // enable QA in AliAnalysisTaskGammaConvV1
    Int_t     enableQAPhotonTask            = 0,        // enable additional QA task
    Bool_t    enableLightOutput             = kFALSE,   // switch to run light output (only essential histograms for afterburner)
    Bool_t    enableTHnSparse               = kFALSE,   // switch on THNsparse
    Int_t     enableTriggerMimicking        = 0,        // enable trigger mimicking
    Bool_t    enableTriggerOverlapRej       = kFALSE,   // enable trigger overlap rejection
    TString   settingMaxFacPtHard           = "3.",       // maximum factor between hardest jet and ptHard generated
    Int_t     debugLevel                    = 0,        // introducing debug levels for grid running
    // settings for weights
    // FPTW:fileNamePtWeights, FMUW:fileNameMultWeights, FMAW:fileNameMatBudWeights, FEPC:fileNamedEdxPostCalib, separate with ;
    TString   fileNameExternalInputs        = "",
    Int_t     doWeightingPart               = 0,        // enable Weighting
    TString   generatorName                 = "DPMJET", // generator Name
    Bool_t    enableMultiplicityWeighting   = kFALSE,   //
    TString   periodNameAnchor              = "",       //
    Int_t     enableMatBudWeightsPi0        = 0,        // 1 = three radial bins, 2 = 10 radial bins
    Bool_t    enableElecDeDxPostCalibration = kFALSE,
    // special settings
    Bool_t    enablePlotVsCentrality        = kFALSE,
    // subwagon config
    TString   additionalTrainConfig         = "0"       // additional counter for trainconfig + special settings
    ) {

  AliCutHandlerPCM cuts;

  Int_t trackMatcherRunningMode       = 0; // CaloTrackMatcher running mode
  TString fileNamePtWeights           = cuts.GetSpecialFileNameFromString (fileNameExternalInputs, "FPTW:");
  TString fileNameMultWeights         = cuts.GetSpecialFileNameFromString (fileNameExternalInputs, "FMUW:");
  TString fileNameMatBudWeights       = cuts.GetSpecialFileNameFromString (fileNameExternalInputs, "FMAW:");
  TString fileNamedEdxPostCalib       = cuts.GetSpecialFileNameFromString (fileNameExternalInputs, "FEPC:");

  TString addTaskName                 = "AddTask_GammaConvV1_pPb";
  TString sAdditionalTrainConfig      = cuts.GetSpecialSettingFromAddConfig(additionalTrainConfig, "", "", addTaskName);
  if (sAdditionalTrainConfig.Atoi() > 0){
    trainConfig = trainConfig + sAdditionalTrainConfig.Atoi();
    cout << "INFO: " << addTaskName.Data() << " running additionalTrainConfig '" << sAdditionalTrainConfig.Atoi() << "', train config: '" << trainConfig << "'" << endl;
  }
  TString corrTaskSetting         = cuts.GetSpecialSettingFromAddConfig(additionalTrainConfig, "CF","", addTaskName);
  if(corrTaskSetting.CompareTo(""))
    cout << "corrTaskSetting: " << corrTaskSetting.Data() << endl;
  if(additionalTrainConfig.Contains("MaterialBudgetWeights"))
    fileNameMatBudWeights         = cuts.GetSpecialSettingFromAddConfig(additionalTrainConfig, "MaterialBudgetWeights",fileNameMatBudWeights, addTaskName);

  TObjArray *rmaxFacPtHardSetting = settingMaxFacPtHard.Tokenize("_");
  if(rmaxFacPtHardSetting->GetEntries()<1){cout << "ERROR: AddTask_GammaConvV1_pPb during parsing of settingMaxFacPtHard String '" << settingMaxFacPtHard.Data() << "'" << endl; return;}
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

  if (debugLevel > 0){
    cout << "enabled debugging for trainconfig: " << debugLevel << endl;
  }
  // ================== GetAnalysisManager ===============================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error(Form("%s_%i", addTaskName.Data(),  trainConfig), "No analysis manager found.");
    return ;
  }

  // ================== GetInputEventHandler =============================
  AliVEventHandler *inputHandler=mgr->GetInputEventHandler();

  //=========  Set Cutnumber for V0Reader ================================
  TString cutnumberPhoton     = photonCutNumberV0Reader.Data();
  TString cutnumberEvent      = "80000003";

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
  //            find input container
  AliAnalysisTaskGammaConvV1 *task=NULL;
  task= new AliAnalysisTaskGammaConvV1(Form("GammaConvV1_%i",trainConfig));
  task->SetIsHeavyIon(isHeavyIon);
  task->SetIsMC(isMC);
  task->SetV0ReaderName(V0ReaderName);
  task->SetLightOutput(enableLightOutput);

  Bool_t doEtaShiftIndCuts = kFALSE;
  TString stringShift = "";

  // new standard configurations  MB
  if(trainConfig == 1){
    cuts.AddCutPCM("80010113", "00200009f9730000dge0400000", "0162103500900000"); // new default
    cuts.AddCutPCM("80010113", "00200009327000008250400000", "0162103500900000"); // new default, no to close
  } else if (trainConfig == 2){
    cuts.AddCutPCM("80010113", "00200009f9730000dge0400000", "0162103500000000"); // new default
    cuts.AddCutPCM("80010113", "00200009327000008250400000", "0162103500000000"); // new default, no to close
  } else if (trainConfig == 3){
    cuts.AddCutPCM("80010123", "00200009f9730000dge0400000", "0162103500900000"); // new default
    cuts.AddCutPCM("80010123", "00200009327000008250400000", "0162103500900000"); // new default, no to close
  } else if (trainConfig == 4){
    cuts.AddCutPCM("80010123", "00200009f9730000dge0400000", "0162103500000000"); // new default
    cuts.AddCutPCM("80010123", "00200009327000008250400000", "0162103500000000"); // new default, no to close
  } else if (trainConfig == 5){ // past-future protection
    cuts.AddCutPCM("80010213", "00200009f9730000dge0400000", "0162103500900000"); // new default, +-2.225\mus no other interaction
    cuts.AddCutPCM("80010213", "00200009327000008250400000", "0162103500900000"); // new default, no to close, +-2.225\mus no other interaction
    cuts.AddCutPCM("80010513", "00200009f9730000dge0400000", "0162103500900000"); // new default, +-1.075\mus no other interaction
    cuts.AddCutPCM("80010513", "00200009327000008250400000", "0162103500900000"); // new default, no to close, +-1.075\mus no other interaction
  } else if (trainConfig == 6){ // past-future protection
    cuts.AddCutPCM("80010213", "00200009f9730000dge0400000", "0162103500000000"); // new default,  +-2.225\mus no other interaction
    cuts.AddCutPCM("80010213", "00200009327000008250400000", "0162103500000000"); // new default, no to close, +-2.225\mus no other interaction
    cuts.AddCutPCM("80010513", "00200009f9730000dge0400000", "0162103500000000"); // new default,  +-1.075\mus no other interaction
    cuts.AddCutPCM("80010513", "00200009327000008250400000", "0162103500000000"); // new default, no to close,  +-1.075\mus no other interaction
  } else if (trainConfig == 7){
    cuts.AddCutPCM("80010123", "00200009f9730000dge0400000", "0162103500900000"); // new default
    cuts.AddCutPCM("80010123", "00200009327000008250400000", "0162103500900000"); // new default, no to close
  } else if (trainConfig == 8){
    cuts.AddCutPCM("80010123", "00200009f9730000dge0400000", "0162103500000000"); // new default
    cuts.AddCutPCM("80010123", "00200009327000008250400000", "0162103500000000"); // new default, no to close

  // default cut all cents without smearing and to close V0
  } else if (trainConfig == 10) {
    cuts.AddCutPCM("80210113", "00200009327000008250400000", "0162103500000000"); // 0-20%
    cuts.AddCutPCM("82410113", "00200009327000008250400000", "0162103500000000"); // 20-40%
    cuts.AddCutPCM("84610113", "00200009327000008250400000", "0162103500000000"); // 40-60%
    cuts.AddCutPCM("86010113", "00200009327000008250400000", "0162103500000000"); // 60-100%
  } else if (trainConfig == 11) { // past future protection: +-2.225\mus no other interaction
    cuts.AddCutPCM("80210213", "00200009327000008250400000", "0162103500000000"); // 0-20%
    cuts.AddCutPCM("82410213", "00200009327000008250400000", "0162103500000000"); // 20-40%
    cuts.AddCutPCM("84610213", "00200009327000008250400000", "0162103500000000"); // 40-60%
    cuts.AddCutPCM("86010213", "00200009327000008250400000", "0162103500000000"); // 60-100%
  } else if (trainConfig == 12) { // past future protection: +-1.075\mus no other interaction
    cuts.AddCutPCM("80210513", "00200009327000008250400000", "0162103500000000"); // 0-20%
    cuts.AddCutPCM("82410513", "00200009327000008250400000", "0162103500000000"); // 20-40%
    cuts.AddCutPCM("84610513", "00200009327000008250400000", "0162103500000000"); // 40-60%
    cuts.AddCutPCM("86010513", "00200009327000008250400000", "0162103500000000"); // 60-100%
  } else if (trainConfig == 13) {
    cuts.AddCutPCM("90010113", "00200009327000008250400000", "0162103500900000"); // new default, no to close
    cuts.AddCutPCM("90210113", "00200009327000008250400000", "0162103500900000"); // new default, no to close
    cuts.AddCutPCM("92410113", "00200009327000008250400000", "0162103500900000"); // new default, no to close
    cuts.AddCutPCM("94610113", "00200009327000008250400000", "0162103500900000"); // new default, no to close
    cuts.AddCutPCM("96010113", "00200009327000008250400000", "0162103500900000"); // new default, no to close
  } else if (trainConfig == 14) {
    cuts.AddCutPCM("e0010113", "00200009327000008250400000", "0162103500900000"); // new default, no to close
    cuts.AddCutPCM("e0210113", "00200009327000008250400000", "0162103500900000"); // new default, no to close
    cuts.AddCutPCM("e2410113", "00200009327000008250400000", "0162103500900000"); // new default, no to close
    cuts.AddCutPCM("e4610113", "00200009327000008250400000", "0162103500900000"); // new default, no to close
    cuts.AddCutPCM("e6010113", "00200009327000008250400000", "0162103500900000"); // new default, no to close

  // new standard configurations 0-20
  } else if (trainConfig == 20){
    cuts.AddCutPCM("80210113", "00200009f9730000dge0400000", "0162103500900000"); //new default
    cuts.AddCutPCM("80210113", "00200009327000008250400000", "0162103500900000"); //new default, no to close
  } else if (trainConfig == 21){
    cuts.AddCutPCM("80210123", "00200009f9730000dge0400000", "0162103500900000"); //new default
    cuts.AddCutPCM("80210123", "00200009327000008250400000", "0162103500900000"); //new default, no to close
  } else if (trainConfig == 22){
    cuts.AddCutPCM("80210123", "00200009f9730000dge0400000", "0162103500000000"); //new default
    cuts.AddCutPCM("80210123", "00200009327000008250400000", "0162103500000000"); //new default, no to close

  // new standard configurations 20-40
  } else if (trainConfig == 30){
    cuts.AddCutPCM("82410113", "00200009f9730000dge0400000", "0162103500900000"); //new default
    cuts.AddCutPCM("82410113", "00200009327000008250400000", "0162103500900000"); //new default, no to close
  } else if (trainConfig == 31){
    cuts.AddCutPCM("82410123", "00200009f9730000dge0400000", "0162103500900000"); //new default
    cuts.AddCutPCM("82410123", "00200009327000008250400000", "0162103500900000"); //new default, no to close
  } else if (trainConfig == 32){
    cuts.AddCutPCM("82410123", "00200009f9730000dge0400000", "0162103500000000"); //new default
    cuts.AddCutPCM("82410123", "00200009327000008250400000", "0162103500000000"); //new default, no to close

  // new standard configurations 40-60
  } else if (trainConfig == 40){
    cuts.AddCutPCM("84610113", "00200009f9730000dge0400000", "0162103500900000"); //new default
    cuts.AddCutPCM("84610113", "00200009327000008250400000", "0162103500900000"); //new default, no to close
  } else if (trainConfig == 41){
    cuts.AddCutPCM("84610123", "00200009f9730000dge0400000", "0162103500900000"); //new default
    cuts.AddCutPCM("84610123", "00200009327000008250400000", "0162103500900000"); //new default, no to close
  } else if (trainConfig == 42){
    cuts.AddCutPCM("84610123", "00200009f9730000dge0400000", "0162103500000000"); //new default
    cuts.AddCutPCM("84610123", "00200009327000008250400000", "0162103500000000"); //new default, no to close

  // new standard configurations 60-80
  } else if (trainConfig == 50){
    cuts.AddCutPCM("86010113", "00200009f9730000dge0400000", "0162103500900000"); //new default
    cuts.AddCutPCM("86010113", "00200009327000008250400000", "0162103500900000"); //new default, no to close
  } else if (trainConfig == 51){
    cuts.AddCutPCM("86010123", "00200009f9730000dge0400000", "0162103500900000"); //new default
    cuts.AddCutPCM("86010123", "00200009327000008250400000", "0162103500900000"); //new default, no to close
  } else if (trainConfig == 52){
    cuts.AddCutPCM("86010123", "00200009f9730000dge0400000", "0162103500000000"); //new default
    cuts.AddCutPCM("86010123", "00200009327000008250400000", "0162103500000000"); //new default, no to close

  //--------------------------------------------------------------------------
  // Systematics variations for standard ana w/o to close V0, wo smearing
  //--------------------------------------------------------------------------
  } else if (trainConfig == 100) { //default + dEdx Variation
    cuts.AddCutPCM("80010113", "00200009327000008250400000", "0162103500000000"); // new default
    cuts.AddCutPCM("80010113", "00200009217000008260400000", "0162103500000000"); // old standard cut
    cuts.AddCutPCM("80010113", "00200009f9730000dge0400000", "0162103500000000"); // new default w/ to clos v0
    cuts.AddCutPCM("80010113", "00200009f9730000d250400000", "0162103500000000"); //
    cuts.AddCutPCM("80010113", "00200009f97300008ge0400000", "0162103500000000"); //
  } else if (trainConfig == 101) { // gamma eta & meson rap var
    cuts.AddCutPCM("80010113", "03200009327000008250400000", "0162303500000000"); // |eta| < 0.65, |y| < 0.6
    cuts.AddCutPCM("80010113", "04200009327000008250400000", "0162203500000000"); // |eta| < 0.75, |y| < 0.7
    cuts.AddCutPCM("80010113", "01200009327000008250400000", "0162403500000000"); // |eta| < 0.6, |y| < 0.5
  } else if (trainConfig == 102) { // minR and single pt var
    cuts.AddCutPCM("80010113", "00100009327000008250400000", "0162103500000000"); // minR 2.8
    cuts.AddCutPCM("80010113", "00900009327000008250400000", "0162103500000000"); // minR 7.5
    cuts.AddCutPCM("80010113", "00200079327000008250400000", "0162103500000000"); // single pT 0. GeV/c
    cuts.AddCutPCM("80010113", "00200019327000008250400000", "0162103500000000"); // single pT 0.1 GeV/c
  } else if (trainConfig == 103) { // TPC cluster & edEdx var
    cuts.AddCutPCM("80010113", "00200008327000008250400000", "0162103500000000"); // TPC Cluster 0.35
    cuts.AddCutPCM("80010113", "00200006327000008250400000", "0162103500000000"); // TPC Cluster 0.7
    cuts.AddCutPCM("80010113", "00200009227000008250400000", "0162103500000000"); // edEdx -4,5
    cuts.AddCutPCM("80010113", "00200009627000008250400000", "0162103500000000"); // edEdx -2.5,4
    cuts.AddCutPCM("80010113", "00200009127000008250400000", "0162103500000000"); // edEdx 5,5
  } else if (trainConfig == 104) { //PidEdx Variation
    cuts.AddCutPCM("80010113", "00200009357000008250400000", "0162103500000000"); // PidEdx(2, >0.4GeV)
    cuts.AddCutPCM("80010113", "00200009387300008250400000", "0162103500000000"); // PidEdx(2, >0.4GeV; 1>3.5GeV)
    cuts.AddCutPCM("80010113", "00200009320000008250400000", "0162103500000000"); // PidEdx(1, >0.5GeV)
    cuts.AddCutPCM("80010113", "00200009325000008250400000", "0162103500000000"); // PidEdx(1, >0.3GeV)
  } else if (trainConfig == 105) { //PidEdx Variation 2
    cuts.AddCutPCM("80010113", "00200009327300008250400000", "0162103500000000"); // PidEdx(1, >0.4GeV; -10>3.5GeV)
    cuts.AddCutPCM("80010113", "00200009326000008250400000", "0162103500000000"); // PidEdx(1, >0.25GeV)
    cuts.AddCutPCM("80010113", "00200009326200008250400000", "0162103500000000"); // PidEdx(1, >0.25GeV; -10>4GeV)
    cuts.AddCutPCM("80010113", "00200009327200008250400000", "0162103500000000"); // PidEdx(1, >0.4GeV; -10>4GeV)
  } else if (trainConfig == 106) { // qt & psipair variation
    cuts.AddCutPCM("80010113", "00200009327000003250400000", "0162103500000000"); // qT 0.05 1D
    cuts.AddCutPCM("80010113", "00200009327000009250400000", "0162103500000000"); // qT 0.03 2D
    cuts.AddCutPCM("80010113", "00200009327000002250400000", "0162103500000000"); // qT 0.07 1D
    cuts.AddCutPCM("80010113", "00200009327000008240400000", "0162103500000000"); // Psi Pair: 1D 0.2
  } else if (trainConfig == 107) { // PsiPair Variation
    cuts.AddCutPCM("80010113", "00200009327000008210400000", "0162103500000000"); // Psi Pair: 1D 0.1
    cuts.AddCutPCM("80010113", "00200009327000008150400000", "0162103500000000"); // chi2 50  2D
    cuts.AddCutPCM("80010113", "00200009327000008850400000", "0162103500000000"); // chi2 20  2D
    cuts.AddCutPCM("80010113", "00200009327000008260400000", "0162103500000000"); // psi pair 0.05
  } else if (trainConfig == 108) { // cos point & meson alpha variation
    cuts.AddCutPCM("80010113", "00200009327000008250300000", "0162103500000000"); // cos pointing angle 0.75
    cuts.AddCutPCM("80010113", "00200009327000008250600000", "0162103500000000"); // cos pointing angle 0.9
    cuts.AddCutPCM("80010113", "00200009327000008250400000", "0162106500000000"); // alpha meson cut 0.8
  } else if (trainConfig == 109) { // BG mixing & smear variations
    cuts.AddCutPCM("80010113", "00200009327000008250400000", "0262103500000000"); // BG track multiplicity
    cuts.AddCutPCM("80010113", "00200009327000008250400000", "0162103500800000"); // fPSigSmearingCte=0.014;
    cuts.AddCutPCM("80010113", "00200009327000008250400000", "0162103500900000"); // fPSigSmearingCte=0.014;

  //--------------------------------------------------------------------------
  // Systematics variations for standard ana w/o to close V0, wo smearing added signals
  //--------------------------------------------------------------------------
  } else if (trainConfig == 120) { //default + dEdx Variation added sig
    cuts.AddCutPCM("80010123", "00200009327000008250400000", "0162103500000000"); // new default
    cuts.AddCutPCM("80010123", "00200009217000008260400000", "0162103500000000"); // old standard cut
    cuts.AddCutPCM("80010123", "00200009f9730000dge0400000", "0162103500000000"); // new default w/ to close V0
  } else if (trainConfig == 121) {  // gamma eta & meson rap var added sig
    cuts.AddCutPCM("80010123", "03200009327000008250400000", "0162303500000000"); // |eta| < 0.65, |y| < 0.6
    cuts.AddCutPCM("80010123", "04200009327000008250400000", "0162203500000000"); // |eta| < 0.75, |y| < 0.7
    cuts.AddCutPCM("80010123", "01200009327000008250400000", "0162403500000000"); // |eta| < 0.6, |y| < 0.5
  } else if (trainConfig == 122) { // minR and single pt var added sig
    cuts.AddCutPCM("80010123", "00100009327000008250400000", "0162103500000000"); // minR 2.8
    cuts.AddCutPCM("80010123", "00900009327000008250400000", "0162103500000000"); // minR 7.5
    cuts.AddCutPCM("80010123", "00200079327000008250400000", "0162103500000000"); // single pT 0. GeV/c
    cuts.AddCutPCM("80010123", "00200019327000008250400000", "0162103500000000"); // single pT 0.1 GeV/c
  } else if (trainConfig == 123) { // TPC cluster & edEdx var add sig
    cuts.AddCutPCM("80010123", "00200008327000008250400000", "0162103500000000"); // TPC Cluster 0.35
    cuts.AddCutPCM("80010123", "00200006327000008250400000", "0162103500000000"); // TPC Cluster 0.7
    cuts.AddCutPCM("80010123", "00200009227000008250400000", "0162103500000000"); // edEdx -3,5
    cuts.AddCutPCM("80010123", "00200009627000008250400000", "0162103500000000"); // edEdx -2.5,4
    cuts.AddCutPCM("80010123", "00200009127000008250400000", "0162103500000000"); // edEdx -5,5
  } else if (trainConfig == 124) { //PidEdx Variation
    cuts.AddCutPCM("80010123", "00200009357000008250400000", "0162103500000000"); // PidEdx(2, >0.4GeV)
    cuts.AddCutPCM("80010123", "00200009387300008250400000", "0162103500000000"); // PidEdx(2, >0.4GeV; 1>3.5GeV)
    cuts.AddCutPCM("80010123", "00200009320000008250400000", "0162103500000000"); // PidEdx(1, >0.5GeV)
    cuts.AddCutPCM("80010123", "00200009325000008250400000", "0162103500000000"); // PidEdx(1, >0.3GeV)
  } else if (trainConfig == 125) { //PidEdx Variation 2
    cuts.AddCutPCM("80010123", "00200009327300008250400000", "0162103500000000"); // PidEdx(1, >0.4GeV; -10>3.5GeV)
    cuts.AddCutPCM("80010123", "00200009326000008250400000", "0162103500000000"); // PidEdx(1, >0.25GeV)
    cuts.AddCutPCM("80010123", "00200009326200008250400000", "0162103500000000"); // PidEdx(1, >0.25GeV; -10>4GeV)
    cuts.AddCutPCM("80010123", "00200009327200008250400000", "0162103500000000"); // PidEdx(1, >0.4GeV; -10>4GeV)
  } else if (trainConfig == 126) { // qt & psipair variation
    cuts.AddCutPCM("80010123", "00200009327000003250400000", "0162103500000000"); // qT 0.05 1D
    cuts.AddCutPCM("80010123", "00200009327000009250400000", "0162103500000000"); // qT 0.03 2D
    cuts.AddCutPCM("80010123", "00200009327000002250400000", "0162103500000000"); // qT 0.07 1D
    cuts.AddCutPCM("80010123", "00200009327000008240400000", "0162103500000000"); // Psi Pair: 1D 0.2
  } else if (trainConfig == 127) { // PsiPair Variation
    cuts.AddCutPCM("80010123", "00200009327000008210400000", "0162103500000000"); // Psi Pair: 1D 0.1
    cuts.AddCutPCM("80010123", "00200009327000008150400000", "0162103500000000"); // chi2 50  2D
    cuts.AddCutPCM("80010123", "00200009327000008850400000", "0162103500000000"); // chi2 20  2D
    cuts.AddCutPCM("80010123", "00200009327000008260400000", "0162103500000000"); // psi pair 0.05
  } else if (trainConfig == 128) { // cos point & meson alpha variation
    cuts.AddCutPCM("80010123", "00200009327000008250300000", "0162103500000000"); // cos pointing angle 0.75
    cuts.AddCutPCM("80010123", "00200009327000008250600000", "0162103500000000"); // cos pointing angle 0.9
    cuts.AddCutPCM("80010123", "00200009327000008250400000", "0162106500000000"); // alpha meson cut 0.8
  } else if (trainConfig == 129) { // BG mixing & smear variations
    cuts.AddCutPCM("80010123", "00200009327000008250400000", "0262103500000000"); // BG track multiplicity
    cuts.AddCutPCM("80010123", "00200009327000008250400000", "0162103500800000"); // fPSigSmearingCte=0.014;
    cuts.AddCutPCM("80010123", "00200009327000008250400000", "0162103500900000"); // fPSigSmearingCte=0.014;
    // pseudorapidity studies
    //-----------------------------------------------------------------------
  } else if (trainConfig == 130) {
    cuts.AddCutPCM("80010113", "0a200009327000008250400000", "0162103500000000"); //Eta cut -0.9 - -0.2 and 0.2 - 0.9
    cuts.AddCutPCM("80010113", "0b200009327000008250400000", "0162103500000000"); //Eta cut -0.9 - -0.2 and 0.2 - 0.9 with LineCut
  } else if (trainConfig == 131) {
    cuts.AddCutPCM("80010123", "0a200009327000008250400000", "0162103500000000"); //Eta cut -0.9 - -0.2 and 0.2 - 0.9 Added Signals
    cuts.AddCutPCM("80010123", "0b200009327000008250400000", "0162103500000000"); //Eta cut -0.9 - -0.2 and 0.2 - 0.9 Added Signals

  //--------------------------------------------------------------------------
  // purity studies (kappa cut)
  //--------------------------------------------------------------------------
  } else if (trainConfig == 150) {
    cuts.AddCutPCM("80010113", "00200009300000008250400000", "0162103500000000"); // -3 < kappa <  5
    cuts.AddCutPCM("80010113", "00200009500000008250400000", "0162103500000000"); // -5 < kappa < 10
    cuts.AddCutPCM("80010113", "00200009600000008250400000", "0162103500000000"); // -3 < kappa < 10
    cuts.AddCutPCM("80010113", "00200009700000008250400000", "0162103500000000"); //  0 < kappa < 10
  } else if (trainConfig == 152) {
    cuts.AddCutPCM("80210113", "00200009300000008250400000", "0162103500000000"); // -3 < kappa <  5
    cuts.AddCutPCM("80210113", "00200009500000008250400000", "0162103500000000"); // -5 < kappa < 10
    cuts.AddCutPCM("80210113", "00200009600000008250400000", "0162103500000000"); // -3 < kappa < 10
    cuts.AddCutPCM("80210113", "00200009700000008250400000", "0162103500000000"); //  0 < kappa < 10
  } else if (trainConfig == 153) {
    cuts.AddCutPCM("82410113", "00200009300000008250400000", "0162103500000000"); // -3 < kappa <  5
    cuts.AddCutPCM("82410113", "00200009500000008250400000", "0162103500000000"); // -5 < kappa < 10
    cuts.AddCutPCM("82410113", "00200009600000008250400000", "0162103500000000"); // -3 < kappa < 10
    cuts.AddCutPCM("82410113", "00200009700000008250400000", "0162103500000000"); //  0 < kappa < 10
  } else if (trainConfig == 154) {
    cuts.AddCutPCM("84610113", "00200009300000008250400000", "0162103500000000"); // -3 < kappa <  5
    cuts.AddCutPCM("84610113", "00200009500000008250400000", "0162103500000000"); // -5 < kappa < 10
    cuts.AddCutPCM("84610113", "00200009600000008250400000", "0162103500000000"); // -3 < kappa < 10
    cuts.AddCutPCM("84610113", "00200009700000008250400000", "0162103500000000"); //  0 < kappa < 10
  } else if (trainConfig == 155) {
    cuts.AddCutPCM("86010113", "00200009300000008250400000", "0162103500000000"); // -3 < kappa <  5
    cuts.AddCutPCM("86010113", "00200009500000008250400000", "0162103500000000"); // -5 < kappa < 10
    cuts.AddCutPCM("86010113", "00200009600000008250400000", "0162103500000000"); // -3 < kappa < 10
    cuts.AddCutPCM("86010113", "00200009700000008250400000", "0162103500000000"); //  0 < kappa < 10

  //--------------------------------------------------------------------------
  // Material weight studies
  //--------------------------------------------------------------------------
  // standard cuts
  } else if (trainConfig == 200) {
    cuts.AddCutPCM("80010113", "00200009327000008250404000", "0162103500000000");
  } else if (trainConfig == 201) {
    cuts.AddCutPCM("80010113", "00200009327000008250404000", "0162103500000000");
  } else if (trainConfig == 202) {
    cuts.AddCutPCM("80010113", "00200009327000008250404000", "0162103500000000");

  // standard cuts added signals
  } else if (trainConfig == 210) {
    cuts.AddCutPCM("80010123", "00200009327000008250404000", "0162103500000000");
  } else if (trainConfig == 211) {
    cuts.AddCutPCM("80010123", "00200009327000008250404000", "0162103500000000");
  } else if (trainConfig == 212) {
    cuts.AddCutPCM("80010123", "00200009327000008250404000", "0162103500000000");

  //--------------------------------------------------------------------------
  // 2016 pPb w/o past future protection
  //--------------------------------------------------------------------------
  // Min Bias
  } else if (trainConfig == 300) {
    cuts.AddCutPCM("80010113", "00200009397302001280004000", "0162103500000000"); // Min Bias
    cuts.AddCutPCM("80010113", "00200009397302001280004000", "0162101500000000"); // Min Bias , alpha pT dependent
  // TRD trigger
  } else if (trainConfig == 301) {
    cuts.AddCutPCM("80047113", "00200009397302001280004000", "0162103500000000"); // TRD trigger HQU for 8TeV
    cuts.AddCutPCM("80043113", "00200009397302001280004000", "0162103500000000"); // TRD trigger HSE for 8TeV
  // Calo triggers
  } else if (trainConfig == 302) {  // EMC triggers
    cuts.AddCutPCM("80010113", "00200009f9730000dge0400000", "0162103500000000", "1111100007032230000"); // Min Bias
    cuts.AddCutPCM("80052113", "00200009f9730000dge0400000", "0162103500000000", "1111100007032230000"); // EMC7
    cuts.AddCutPCM("80085113", "00200009f9730000dge0400000", "0162103500000000", "1111100007032230000"); // EG2
    cuts.AddCutPCM("80083113", "00200009f9730000dge0400000", "0162103500000000", "1111100007032230000"); // EG1
  } else if (trainConfig == 303) {  // PHOS triggers
    cuts.AddCutPCM("80010113", "00200009f9730000dge0400000", "0162103500000000", "2444400041013200000"); // MinBias
    cuts.AddCutPCM("80052113", "00200009f9730000dge0400000", "0162103500000000", "2444400041013200000"); // PHI7
  //
  } else if (trainConfig == 304) {
    cuts.AddCutPCM("80110113", "00200009f9730000dge0400000", "0162103500000000"); // 0-10
    cuts.AddCutPCM("81210113", "00200009f9730000dge0400000", "0162103500000000"); // 0-20
    cuts.AddCutPCM("82410113", "00200009f9730000dge0400000", "0162103500000000"); // 20-40
    cuts.AddCutPCM("84610113", "00200009f9730000dge0400000", "0162103500000000"); // 40-60
    cuts.AddCutPCM("86810113", "00200009f9730000dge0400000", "0162103500000000"); // 60-80
    cuts.AddCutPCM("88010113", "00200009f9730000dge0400000", "0162103500000000"); // 80-100
  } else if (trainConfig == 305) {
    cuts.AddCutPCM("80010113", "00200009f9730000dge0400000", "0162103500000000"); // new default for 5TeV
    cuts.AddCutPCM("80210113", "00200009f9730000dge0400000", "0162103500000000"); // 0-20
    cuts.AddCutPCM("86010113", "00200009f9730000dge0400000", "0162103500000000"); // 60-100
    cuts.AddCutPCM("a0110113", "00200009f9730000dge0400000", "0162103500000000"); // 0-5
    cuts.AddCutPCM("a1210113", "00200009f9730000dge0400000", "0162103500000000"); // 5-10

  } else if (trainConfig == 306){ // pPb 2013 TeV defaults
    cuts.AddCutPCM("80110113", "00200009327000008250400000", "0162103500000000"); // 0-10
    cuts.AddCutPCM("81210113", "00200009327000008250400000", "0162103500000000"); // 10-20
    cuts.AddCutPCM("82410113", "00200009327000008250400000", "0162103500000000"); // 20-40
    cuts.AddCutPCM("84610113", "00200009327000008250400000", "0162103500000000"); // 40-60
    cuts.AddCutPCM("86810113", "00200009327000008250400000", "0162103500000000"); // 60-80
    cuts.AddCutPCM("88010113", "00200009327000008250400000", "0162103500000000"); // 80-100
  } else if (trainConfig == 307){ // pPb 2013 TeV defaults
    cuts.AddCutPCM("80210113", "00200009327000008250400000", "0162103500000000"); // 0-20
    cuts.AddCutPCM("86010113", "00200009327000008250400000", "0162103500000000"); // 60-100
    cuts.AddCutPCM("a0110113", "00200009327000008250400000", "0162103500000000"); // 0-5
    cuts.AddCutPCM("a1210113", "00200009327000008250400000", "0162103500000000"); // 5-10

  } else if (trainConfig == 308) { // pPb 5 TeV Run2 narrow cent
    cuts.AddCutPCM("c0110113", "00200009f9730000dge0400000", "0162103500000000"); // 0-1
    cuts.AddCutPCM("c0210113", "00200009f9730000dge0400000", "0162103500000000"); // 0-2
  //--------------------------------------------------------------------------
  // 2016 pPb w/ past future protection 2.24 \mus protected
  //--------------------------------------------------------------------------
  // Min Bias
  } else if (trainConfig == 310) {
    cuts.AddCutPCM("80010213", "00200009397302001280004000", "0162103500000000"); // Min Bias
    cuts.AddCutPCM("80010213", "00200009397302001280004000", "0162101500000000"); // Min Bias , alpha pT dependent
  // TRD trigger
  } else if (trainConfig == 311) {
    cuts.AddCutPCM("80047213", "00200009397302001280004000", "0162103500000000"); // TRD trigger HQU for 8TeV
    cuts.AddCutPCM("80043213", "00200009397302001280004000", "0162103500000000"); // TRD trigger HSE for 8TeV
  // Calo triggers
  } else if (trainConfig == 312) {  // EMC triggers
    cuts.AddCutPCM("80010213", "00200009f9730000dge0400000", "0162103500000000", "1111100007032230000"); // Min Bias
    cuts.AddCutPCM("80052213", "00200009f9730000dge0400000", "0162103500000000", "1111100007032230000"); // EMC7
    cuts.AddCutPCM("80085213", "00200009f9730000dge0400000", "0162103500000000", "1111100007032230000"); // EG2
    cuts.AddCutPCM("80083213", "00200009f9730000dge0400000", "0162103500000000", "1111100007032230000"); // EG1
  } else if (trainConfig == 313) {  // PHOS triggers
    cuts.AddCutPCM("80010213", "00200009f9730000dge0400000", "0162103500000000", "2444400041013200000"); // MinBias
    cuts.AddCutPCM("80052213", "00200009f9730000dge0400000", "0162103500000000", "2444400041013200000"); // PHI7
  } else if (trainConfig == 314) {
    cuts.AddCutPCM("80010213", "00200009f9730000dge0400000", "0162103500000000"); // new default for 5TeV
    cuts.AddCutPCM("80010213", "00200009f9730000dge0400000", "0162101500000000"); // new default, alpha pT dependent for 5TeV


  } else if (trainConfig == 360) {
    cuts.AddCutPCM("80010113", "00200009f9730000dge0400000", "0152103500000000"); //New standard cut for 8TeV analysis V0AND with double counting cut, TOF removed
    cuts.AddCutPCM("00010013", "00200009f9730000dge0400000", "0152103500000000"); // no SPD pileup cut
    cuts.AddCutPCM("80010113", "00100009f9730000dge0400000", "0152103500000000"); // R cut 2.8 -180 cm
    cuts.AddCutPCM("80010113", "00500009f9730000dge0400000", "0152103500000000"); // R cut 10. -180 cm
  } else if (trainConfig == 361) {
    cuts.AddCutPCM("80010113", "00200069f9730000dge0400000", "0152103500000000"); // min pT 40 MeV
    cuts.AddCutPCM("80010113", "00200049f9730000dge0400000", "0152103500000000"); // min pT 75 MeV
    cuts.AddCutPCM("80010113", "00200019f9730000dge0400000", "0152103500000000"); // min pT 100MeV
  } else if (trainConfig == 362) {
    cuts.AddCutPCM("80010113", "00200068f9730000dge0400000", "0152103500000000"); // TPC cluster 35%
    cuts.AddCutPCM("80010113", "00200066f9730000dge0400000", "0152103500000000"); // TPC cluster 70%
    cuts.AddCutPCM("80010113", "00200009f9730000dge0600000", "0152103500000000"); // cosPA 0.9
    cuts.AddCutPCM("80010113", "00200009f9730000dge0300000", "0152103500000000"); // cosPA 0.75
  } else if (trainConfig == 363) {
    cuts.AddCutPCM("80010113", "0020000939730000dge0400000", "0152103500000000"); // nsig electron   -4,5
    cuts.AddCutPCM("80010113", "0020000969730000dge0400000", "0152103500000000"); // nsig electron -2.5,4
    cuts.AddCutPCM("80010113", "00200009f5730000dge0400000", "0152103500000000"); // nsig pion 2,-10
    cuts.AddCutPCM("80010113", "00200009f1730000dge0400000", "0152103500000000"); // nsig pion 0,-10
  } else if (trainConfig == 364) {
    cuts.AddCutPCM("80010113", "00200009f9030000dge0400000", "0152103500000000"); // pion nsig min mom 0.50 GeV/c
    cuts.AddCutPCM("80010113", "00200009f9630000dge0400000", "0152103500000000"); // pion nsig min mom 0.25 GeV/c
    cuts.AddCutPCM("80010113", "00200009f9760000dge0400000", "0152103500000000"); // pion nsig max mom 2.00 GeV/c
    cuts.AddCutPCM("80010113", "00200009f9710000dge0400000", "0152103500000000"); // pion nsig max mom 5.00 GeV/c
  } else if (trainConfig == 365) {
    cuts.AddCutPCM("80010113", "00200009f9730000age0400000", "0152103500000000"); // qT max 0.040, qtptmax 0.11
    cuts.AddCutPCM("80010113", "00200009f9730000ege0400000", "0152103500000000"); // qT max 0.060, qtptmax 0.14
    cuts.AddCutPCM("80010113", "00200009f9730000fge0400000", "0152103500000000"); // qT max 0.070, qtptmax 0.16
  } else if (trainConfig == 366) {
    cuts.AddCutPCM("80010113", "00200009f9730000d1e0400000", "0152103500000000"); // chi2 50
    cuts.AddCutPCM("80010113", "00200009f9730000dfe0400000", "0152103500000000"); // chi2 50 chi2 dep -0.065
    cuts.AddCutPCM("80010113", "00200009f9730000dhe0400000", "0152103500000000"); // chi2 50 chi2 dep -0.050
    cuts.AddCutPCM("80010113", "00200009f9730000dge0404000", "0152103500000000"); // reject close v0
    cuts.AddCutPCM("80010113", "00200009f9730000dge0406000", "0152103500000000"); // double count with open angle 0.04
  } else if (trainConfig == 367) {
    cuts.AddCutPCM("80010113", "00200009f9730000dgd0400000", "0152103500000000"); // Psi pair 0.15 dep
    cuts.AddCutPCM("80010113", "00200009f9730000dgf0400000", "0152103500000000"); // Psi pair 0.20 dep
    cuts.AddCutPCM("80010113", "00200009f9730000dgg0400000", "0152103500000000"); // Psi pair 0.30 dep
  } else if (trainConfig == 368) {
    cuts.AddCutPCM("80010113", "00200009f9730000dge0400000", "0252103500000000"); // variation BG scheme track mult
    cuts.AddCutPCM("80010113", "00200009f9730000dge0400000", "0152107500000000"); // alpha meson 0.85
    cuts.AddCutPCM("80010113", "00200009f9730000dge0400000", "0152105500000000"); // alpha meson 0.75
    cuts.AddCutPCM("80010113", "00200009227300008250404000", "0152103500000000"); // old cuts (run1)
  } else if (trainConfig == 369) {
    cuts.AddCutPCM("80010213", "00200009f9730000dge0400000", "0152103500000000"); //same as std + maximum past future rejection
    cuts.AddCutPCM("80010513", "00200009f9730000dge0400000", "0152103500000000"); //same as std + medium past future rejection
  } else if (trainConfig == 370) {
    cuts.AddCutPCM("80010113", "00200009f9730000dge0400000", "0152103500000000"); 
    cuts.AddCutPCM("80010113", "00200009f9730000dge0400000", "0r52103500000000"); 
    cuts.AddCutPCM("80010113", "00200009f9730000dge0400000", "0s52103500000000"); 
    cuts.AddCutPCM("80010113", "00200009f9730000dge0400000", "0t52103500000000"); 
    cuts.AddCutPCM("80010113", "00200009f9730000dge0400000", "0u52103500000000"); 
  } else if (trainConfig == 371) {
    cuts.AddCutPCM("80010113", "00200009f9730000dge0400000", "0152103500900000"); 
    cuts.AddCutPCM("80010113", "00200009f9730000dge0400000", "0r52103500900000"); 
    cuts.AddCutPCM("80010113", "00200009f9730000dge0400000", "0s52103500900000"); 
    cuts.AddCutPCM("80010113", "00200009f9730000dge0400000", "0t52103500900000"); 
    cuts.AddCutPCM("80010113", "00200009f9730000dge0400000", "0u52103500900000"); 
  } else if (trainConfig == 372) {
    cuts.AddCutPCM("80010113", "00200009f9730000dge0400000", "0152103500800000"); 
    cuts.AddCutPCM("80010113", "00200009f9730000dge0400000", "0r52103500800000"); 
    cuts.AddCutPCM("80010113", "00200009f9730000dge0400000", "0s52103500800000"); 
    cuts.AddCutPCM("80010113", "00200009f9730000dge0400000", "0t52103500800000"); 
    cuts.AddCutPCM("80010113", "00200009f9730000dge0400000", "0u52103500800000"); 

  //--------------------------------------------------------------------------
  // Systematics variations for direct photons pPb 5 TeV Run1 and Run2 - 0-100%
  //--------------------------------------------------------------------------
  } else if (trainConfig == 400) { //default + dEdx Variation
    cuts.AddCutPCM("80010113", "00200009327000008250400000", "0162103500000000"); // new default
    cuts.AddCutPCM("80010113", "00200009217000008260400000", "0162103500000000"); // old standard cut
    cuts.AddCutPCM("80010113", "00200009f9730000dge0400000", "0162103500000000"); // new default w/ to clos v0
    cuts.AddCutPCM("80010113", "00200009f9730000d250400000", "0162103500000000"); //
    cuts.AddCutPCM("80010113", "00200009f97300008ge0400000", "0162103500000000"); //
  } else if (trainConfig == 401) { // gamma eta & meson rap var
    cuts.AddCutPCM("80010113", "03200009327000008250400000", "0162303500000000"); // |eta| < 0.65, |y| < 0.6
    cuts.AddCutPCM("80010113", "04200009327000008250400000", "0162203500000000"); // |eta| < 0.75, |y| < 0.7
    cuts.AddCutPCM("80010113", "01200009327000008250400000", "0162403500000000"); // |eta| < 0.6, |y| < 0.5
  } else if (trainConfig == 402) { // minR and single pt var
    cuts.AddCutPCM("80010113", "00100009327000008250400000", "0162103500000000"); // minR 2.8
    cuts.AddCutPCM("80010113", "00900009327000008250400000", "0162103500000000"); // minR 7.5
    cuts.AddCutPCM("80010113", "00200079327000008250400000", "0162103500000000"); // single pT 0. GeV/c
    cuts.AddCutPCM("80010113", "00200019327000008250400000", "0162103500000000"); // single pT 0.1 GeV/c
  } else if (trainConfig == 403) { // TPC cluster & edEdx var
    cuts.AddCutPCM("80010113", "00200008327000008250400000", "0162103500000000"); // TPC Cluster 0.35
    cuts.AddCutPCM("80010113", "00200006327000008250400000", "0162103500000000"); // TPC Cluster 0.7
    cuts.AddCutPCM("80010113", "00200009227000008250400000", "0162103500000000"); // edEdx -4,5
    cuts.AddCutPCM("80010113", "00200009627000008250400000", "0162103500000000"); // edEdx -2.5,4
    cuts.AddCutPCM("80010113", "00200009127000008250400000", "0162103500000000"); // edEdx 5,5
  } else if (trainConfig == 404) { //PidEdx Variation
    cuts.AddCutPCM("80010113", "00200009357000008250400000", "0162103500000000"); // PidEdx(2, >0.4GeV)
    cuts.AddCutPCM("80010113", "00200009387300008250400000", "0162103500000000"); // PidEdx(2, >0.4GeV; 1>3.5GeV)
    cuts.AddCutPCM("80010113", "00200009320000008250400000", "0162103500000000"); // PidEdx(1, >0.5GeV)
    cuts.AddCutPCM("80010113", "00200009325000008250400000", "0162103500000000"); // PidEdx(1, >0.3GeV)
  } else if (trainConfig == 405) { //PidEdx Variation 2
    cuts.AddCutPCM("80010113", "00200009327300008250400000", "0162103500000000"); // PidEdx(1, >0.4GeV; -10>3.5GeV)
    cuts.AddCutPCM("80010113", "00200009326000008250400000", "0162103500000000"); // PidEdx(1, >0.25GeV)
    cuts.AddCutPCM("80010113", "00200009326200008250400000", "0162103500000000"); // PidEdx(1, >0.25GeV; -10>4GeV)
    cuts.AddCutPCM("80010113", "00200009327200008250400000", "0162103500000000"); // PidEdx(1, >0.4GeV; -10>4GeV)
  } else if (trainConfig == 406) { // qt & psipair variation
    cuts.AddCutPCM("80010113", "00200009327000003250400000", "0162103500000000"); // qT 0.05 1D
    cuts.AddCutPCM("80010113", "00200009327000009250400000", "0162103500000000"); // qT 0.03 2D
    cuts.AddCutPCM("80010113", "00200009327000002250400000", "0162103500000000"); // qT 0.07 1D
    cuts.AddCutPCM("80010113", "00200009327000008240400000", "0162103500000000"); // Psi Pair: 1D 0.2
  } else if (trainConfig == 407) { // PsiPair Variation
    cuts.AddCutPCM("80010113", "00200009327000008210400000", "0162103500000000"); // Psi Pair: 1D 0.1
    cuts.AddCutPCM("80010113", "00200009327000008150400000", "0162103500000000"); // chi2 50  2D
    cuts.AddCutPCM("80010113", "00200009327000008850400000", "0162103500000000"); // chi2 20  2D
    cuts.AddCutPCM("80010113", "00200009327000008260400000", "0162103500000000"); // psi pair 0.05
  } else if (trainConfig == 408) { // cos point & meson alpha variation
    cuts.AddCutPCM("80010113", "00200009327000008250300000", "0162103500000000"); // cos pointing angle 0.75
    cuts.AddCutPCM("80010113", "00200009327000008250600000", "0162103500000000"); // cos pointing angle 0.9
    cuts.AddCutPCM("80010113", "00200009327000008250400000", "0162106500000000"); // alpha meson cut 0.8
  } else if (trainConfig == 409) { // BG mixing & smear variations
    cuts.AddCutPCM("80010113", "00200009327000008250400000", "0262103500000000"); // BG track multiplicity
    cuts.AddCutPCM("80010113", "00200009327000008250400000", "0162103500800000"); // fPSigSmearingCte=0.014;
    cuts.AddCutPCM("80010113", "00200009327000008250400000", "0162103500900000"); // fPSigSmearingCte=0.014;

  //--------------------------------------------------------------------------
  // Systematics variations for direct photons pPb 5 TeV Run1 and Run2 - 0-20%
  //--------------------------------------------------------------------------
  } else if (trainConfig == 410) { //default + dEdx Variation
    cuts.AddCutPCM("80210113", "00200009327000008250400000", "0162103500000000"); // new default
    cuts.AddCutPCM("80210113", "00200009217000008260400000", "0162103500000000"); // old standard cut
    cuts.AddCutPCM("80210113", "00200009f9730000dge0400000", "0162103500000000"); // new default w/ to clos v0
  } else if (trainConfig == 411) { // gamma eta & meson rap var
    cuts.AddCutPCM("80210113", "03200009327000008250400000", "0162303500000000"); // |eta| < 0.65, |y| < 0.6
    cuts.AddCutPCM("80210113", "04200009327000008250400000", "0162203500000000"); // |eta| < 0.75, |y| < 0.7
    cuts.AddCutPCM("80210113", "01200009327000008250400000", "0162403500000000"); // |eta| < 0.6, |y| < 0.5
  } else if (trainConfig == 412) { // minR and single pt var
    cuts.AddCutPCM("80210113", "00100009327000008250400000", "0162103500000000"); // minR 2.8
    cuts.AddCutPCM("80210113", "00900009327000008250400000", "0162103500000000"); // minR 7.5
    cuts.AddCutPCM("80210113", "00200079327000008250400000", "0162103500000000"); // single pT 0. GeV/c
    cuts.AddCutPCM("80210113", "00200019327000008250400000", "0162103500000000"); // single pT 0.1 GeV/c
  } else if (trainConfig == 413) { // TPC cluster & edEdx var
    cuts.AddCutPCM("80210113", "00200008327000008250400000", "0162103500000000"); // TPC Cluster 0.35
    cuts.AddCutPCM("80210113", "00200006327000008250400000", "0162103500000000"); // TPC Cluster 0.7
    cuts.AddCutPCM("80210113", "00200009227000008250400000", "0162103500000000"); // edEdx -4,5
    cuts.AddCutPCM("80210113", "00200009627000008250400000", "0162103500000000"); // edEdx -2.5,4
    cuts.AddCutPCM("80210113", "00200009127000008250400000", "0162103500000000"); // edEdx 5,5
  } else if (trainConfig == 414) { //PidEdx Variation
    cuts.AddCutPCM("80210113", "00200009357000008250400000", "0162103500000000"); // PidEdx(2, >0.4GeV)
    cuts.AddCutPCM("80210113", "00200009387300008250400000", "0162103500000000"); // PidEdx(2, >0.4GeV; 1>3.5GeV)
    cuts.AddCutPCM("80210113", "00200009320000008250400000", "0162103500000000"); // PidEdx(1, >0.5GeV)
    cuts.AddCutPCM("80210113", "00200009325000008250400000", "0162103500000000"); // PidEdx(1, >0.3GeV)
  } else if (trainConfig == 415) { //PidEdx Variation 2
    cuts.AddCutPCM("80210113", "00200009327300008250400000", "0162103500000000"); // PidEdx(1, >0.4GeV; -10>3.5GeV)
    cuts.AddCutPCM("80210113", "00200009326000008250400000", "0162103500000000"); // PidEdx(1, >0.25GeV)
    cuts.AddCutPCM("80210113", "00200009326200008250400000", "0162103500000000"); // PidEdx(1, >0.25GeV; -10>4GeV)
    cuts.AddCutPCM("80210113", "00200009327200008250400000", "0162103500000000"); // PidEdx(1, >0.4GeV; -10>4GeV)
  } else if (trainConfig == 416) { // qt & psipair variation
    cuts.AddCutPCM("80210113", "00200009327000003250400000", "0162103500000000"); // qT 0.05 1D
    cuts.AddCutPCM("80210113", "00200009327000009250400000", "0162103500000000"); // qT 0.03 2D
    cuts.AddCutPCM("80210113", "00200009327000002250400000", "0162103500000000"); // qT 0.07 1D
    cuts.AddCutPCM("80210113", "00200009327000008240400000", "0162103500000000"); // Psi Pair: 1D 0.2
  } else if (trainConfig == 417) { // PsiPair Variation
    cuts.AddCutPCM("80210113", "00200009327000008210400000", "0162103500000000"); // Psi Pair: 1D 0.1
    cuts.AddCutPCM("80210113", "00200009327000008150400000", "0162103500000000"); // chi2 50  2D
    cuts.AddCutPCM("80210113", "00200009327000008850400000", "0162103500000000"); // chi2 20  2D
    cuts.AddCutPCM("80210113", "00200009327000008260400000", "0162103500000000"); // psi pair 0.05
  } else if (trainConfig == 418) { // cos point & meson alpha variation
    cuts.AddCutPCM("80210113", "00200009327000008250300000", "0162103500000000"); // cos pointing angle 0.75
    cuts.AddCutPCM("80210113", "00200009327000008250600000", "0162103500000000"); // cos pointing angle 0.9
    cuts.AddCutPCM("80210113", "00200009327000008250400000", "0162106500000000"); // alpha meson cut 0.8
  } else if (trainConfig == 419) { // BG mixing & smear variations
    cuts.AddCutPCM("80210113", "00200009327000008250400000", "0262103500000000"); // BG track multiplicity
    cuts.AddCutPCM("80210113", "00200009327000008250400000", "0162103500800000"); // fPSigSmearingCte=0.014;
    cuts.AddCutPCM("80210113", "00200009327000008250400000", "0162103500900000"); // fPSigSmearingCte=0.014;

  //--------------------------------------------------------------------------
  // Systematics variations for direct photons pPb 5 TeV Run1 and Run2 - 20-40%
  //--------------------------------------------------------------------------
  } else if (trainConfig == 420) { //default + dEdx Variation
    cuts.AddCutPCM("82410113", "00200009327000008250400000", "0162103500000000"); // new default
    cuts.AddCutPCM("82410113", "00200009217000008260400000", "0162103500000000"); // old standard cut
    cuts.AddCutPCM("82410113", "00200009f9730000dge0400000", "0162103500000000"); // new default w/ to clos v0
  } else if (trainConfig == 421) { // gamma eta & meson rap var
    cuts.AddCutPCM("82410113", "03200009327000008250400000", "0162303500000000"); // |eta| < 0.65, |y| < 0.6
    cuts.AddCutPCM("82410113", "04200009327000008250400000", "0162203500000000"); // |eta| < 0.75, |y| < 0.7
    cuts.AddCutPCM("82410113", "01200009327000008250400000", "0162403500000000"); // |eta| < 0.6, |y| < 0.5
  } else if (trainConfig == 422) { // minR and single pt var
    cuts.AddCutPCM("82410113", "00100009327000008250400000", "0162103500000000"); // minR 2.8
    cuts.AddCutPCM("82410113", "00900009327000008250400000", "0162103500000000"); // minR 7.5
    cuts.AddCutPCM("82410113", "00200079327000008250400000", "0162103500000000"); // single pT 0. GeV/c
    cuts.AddCutPCM("82410113", "00200019327000008250400000", "0162103500000000"); // single pT 0.1 GeV/c
  } else if (trainConfig == 423) { // TPC cluster & edEdx var
    cuts.AddCutPCM("82410113", "00200008327000008250400000", "0162103500000000"); // TPC Cluster 0.35
    cuts.AddCutPCM("82410113", "00200006327000008250400000", "0162103500000000"); // TPC Cluster 0.7
    cuts.AddCutPCM("82410113", "00200009227000008250400000", "0162103500000000"); // edEdx -4,5
    cuts.AddCutPCM("82410113", "00200009627000008250400000", "0162103500000000"); // edEdx -2.5,4
    cuts.AddCutPCM("82410113", "00200009127000008250400000", "0162103500000000"); // edEdx 5,5
  } else if (trainConfig == 424) { //PidEdx Variation
    cuts.AddCutPCM("82410113", "00200009357000008250400000", "0162103500000000"); // PidEdx(2, >0.4GeV)
    cuts.AddCutPCM("82410113", "00200009387300008250400000", "0162103500000000"); // PidEdx(2, >0.4GeV; 1>3.5GeV)
    cuts.AddCutPCM("82410113", "00200009320000008250400000", "0162103500000000"); // PidEdx(1, >0.5GeV)
    cuts.AddCutPCM("82410113", "00200009325000008250400000", "0162103500000000"); // PidEdx(1, >0.3GeV)
  } else if (trainConfig == 425) { //PidEdx Variation 2
    cuts.AddCutPCM("82410113", "00200009327300008250400000", "0162103500000000"); // PidEdx(1, >0.4GeV; -10>3.5GeV)
    cuts.AddCutPCM("82410113", "00200009326000008250400000", "0162103500000000"); // PidEdx(1, >0.25GeV)
    cuts.AddCutPCM("82410113", "00200009326200008250400000", "0162103500000000"); // PidEdx(1, >0.25GeV; -10>4GeV)
    cuts.AddCutPCM("82410113", "00200009327200008250400000", "0162103500000000"); // PidEdx(1, >0.4GeV; -10>4GeV)
  } else if (trainConfig == 426) { // qt & psipair variation
    cuts.AddCutPCM("82410113", "00200009327000003250400000", "0162103500000000"); // qT 0.05 1D
    cuts.AddCutPCM("82410113", "00200009327000009250400000", "0162103500000000"); // qT 0.03 2D
    cuts.AddCutPCM("82410113", "00200009327000002250400000", "0162103500000000"); // qT 0.07 1D
    cuts.AddCutPCM("82410113", "00200009327000008240400000", "0162103500000000"); // Psi Pair: 1D 0.2
  } else if (trainConfig == 427) { // PsiPair Variation
    cuts.AddCutPCM("82410113", "00200009327000008210400000", "0162103500000000"); // Psi Pair: 1D 0.1
    cuts.AddCutPCM("82410113", "00200009327000008150400000", "0162103500000000"); // chi2 50  2D
    cuts.AddCutPCM("82410113", "00200009327000008850400000", "0162103500000000"); // chi2 20  2D
    cuts.AddCutPCM("82410113", "00200009327000008260400000", "0162103500000000"); // psi pair 0.05
  } else if (trainConfig == 428) { // cos point & meson alpha variation
    cuts.AddCutPCM("82410113", "00200009327000008250300000", "0162103500000000"); // cos pointing angle 0.75
    cuts.AddCutPCM("82410113", "00200009327000008250600000", "0162103500000000"); // cos pointing angle 0.9
    cuts.AddCutPCM("82410113", "00200009327000008250400000", "0162106500000000"); // alpha meson cut 0.8
  } else if (trainConfig == 429) { // BG mixing & smear variations
    cuts.AddCutPCM("82410113", "00200009327000008250400000", "0262103500000000"); // BG track multiplicity
    cuts.AddCutPCM("82410113", "00200009327000008250400000", "0162103500800000"); // fPSigSmearingCte=0.014;
    cuts.AddCutPCM("82410113", "00200009327000008250400000", "0162103500900000"); // fPSigSmearingCte=0.014;

  //--------------------------------------------------------------------------
  // Systematics variations for direct photons pPb 5 TeV Run1 and Run2 - 40-60%
  //--------------------------------------------------------------------------
  } else if (trainConfig == 430) { //default + dEdx Variation
    cuts.AddCutPCM("84610113", "00200009327000008250400000", "0162103500000000"); // new default
    cuts.AddCutPCM("84610113", "00200009217000008260400000", "0162103500000000"); // old standard cut
    cuts.AddCutPCM("84610113", "00200009f9730000dge0400000", "0162103500000000"); // new default w/ to clos v0
  } else if (trainConfig == 431) { // gamma eta & meson rap var
    cuts.AddCutPCM("84610113", "03200009327000008250400000", "0162303500000000"); // |eta| < 0.65, |y| < 0.6
    cuts.AddCutPCM("84610113", "04200009327000008250400000", "0162203500000000"); // |eta| < 0.75, |y| < 0.7
    cuts.AddCutPCM("84610113", "01200009327000008250400000", "0162403500000000"); // |eta| < 0.6, |y| < 0.5
  } else if (trainConfig == 432) { // minR and single pt var
    cuts.AddCutPCM("84610113", "00100009327000008250400000", "0162103500000000"); // minR 2.8
    cuts.AddCutPCM("84610113", "00900009327000008250400000", "0162103500000000"); // minR 7.5
    cuts.AddCutPCM("84610113", "00200079327000008250400000", "0162103500000000"); // single pT 0. GeV/c
    cuts.AddCutPCM("84610113", "00200019327000008250400000", "0162103500000000"); // single pT 0.1 GeV/c
  } else if (trainConfig == 433) { // TPC cluster & edEdx var
    cuts.AddCutPCM("84610113", "00200008327000008250400000", "0162103500000000"); // TPC Cluster 0.35
    cuts.AddCutPCM("84610113", "00200006327000008250400000", "0162103500000000"); // TPC Cluster 0.7
    cuts.AddCutPCM("84610113", "00200009227000008250400000", "0162103500000000"); // edEdx -4,5
    cuts.AddCutPCM("84610113", "00200009627000008250400000", "0162103500000000"); // edEdx -2.5,4
    cuts.AddCutPCM("84610113", "00200009127000008250400000", "0162103500000000"); // edEdx 5,5
  } else if (trainConfig == 434) { //PidEdx Variation
    cuts.AddCutPCM("84610113", "00200009357000008250400000", "0162103500000000"); // PidEdx(2, >0.4GeV)
    cuts.AddCutPCM("84610113", "00200009387300008250400000", "0162103500000000"); // PidEdx(2, >0.4GeV; 1>3.5GeV)
    cuts.AddCutPCM("84610113", "00200009320000008250400000", "0162103500000000"); // PidEdx(1, >0.5GeV)
    cuts.AddCutPCM("84610113", "00200009325000008250400000", "0162103500000000"); // PidEdx(1, >0.3GeV)
  } else if (trainConfig == 435) { //PidEdx Variation 2
    cuts.AddCutPCM("84610113", "00200009327300008250400000", "0162103500000000"); // PidEdx(1, >0.4GeV; -10>3.5GeV)
    cuts.AddCutPCM("84610113", "00200009326000008250400000", "0162103500000000"); // PidEdx(1, >0.25GeV)
    cuts.AddCutPCM("84610113", "00200009326200008250400000", "0162103500000000"); // PidEdx(1, >0.25GeV; -10>4GeV)
    cuts.AddCutPCM("84610113", "00200009327200008250400000", "0162103500000000"); // PidEdx(1, >0.4GeV; -10>4GeV)
  } else if (trainConfig == 436) { // qt & psipair variation
    cuts.AddCutPCM("84610113", "00200009327000003250400000", "0162103500000000"); // qT 0.05 1D
    cuts.AddCutPCM("84610113", "00200009327000009250400000", "0162103500000000"); // qT 0.03 2D
    cuts.AddCutPCM("84610113", "00200009327000002250400000", "0162103500000000"); // qT 0.07 1D
    cuts.AddCutPCM("84610113", "00200009327000008240400000", "0162103500000000"); // Psi Pair: 1D 0.2
  } else if (trainConfig == 437) { // PsiPair Variation
    cuts.AddCutPCM("84610113", "00200009327000008210400000", "0162103500000000"); // Psi Pair: 1D 0.1
    cuts.AddCutPCM("84610113", "00200009327000008150400000", "0162103500000000"); // chi2 50  2D
    cuts.AddCutPCM("84610113", "00200009327000008850400000", "0162103500000000"); // chi2 20  2D
    cuts.AddCutPCM("84610113", "00200009327000008260400000", "0162103500000000"); // psi pair 0.05
  } else if (trainConfig == 438) { // cos point & meson alpha variation
    cuts.AddCutPCM("84610113", "00200009327000008250300000", "0162103500000000"); // cos pointing angle 0.75
    cuts.AddCutPCM("84610113", "00200009327000008250600000", "0162103500000000"); // cos pointing angle 0.9
    cuts.AddCutPCM("84610113", "00200009327000008250400000", "0162106500000000"); // alpha meson cut 0.8
  } else if (trainConfig == 439) { // BG mixing & smear variations
    cuts.AddCutPCM("84610113", "00200009327000008250400000", "0262103500000000"); // BG track multiplicity
    cuts.AddCutPCM("84610113", "00200009327000008250400000", "0162103500800000"); // fPSigSmearingCte=0.014;
    cuts.AddCutPCM("84610113", "00200009327000008250400000", "0162103500900000"); // fPSigSmearingCte=0.014;

  //--------------------------------------------------------------------------
  // Systematics variations for direct photons pPb 5 TeV Run1 and Run2 - 60-100%
  //--------------------------------------------------------------------------
  } else if (trainConfig == 440) { //default + dEdx Variation
    cuts.AddCutPCM("86010113", "00200009327000008250400000", "0162103500000000"); // new default
    cuts.AddCutPCM("86010113", "00200009217000008260400000", "0162103500000000"); // old standard cut
    cuts.AddCutPCM("86010113", "00200009f9730000dge0400000", "0162103500000000"); // new default w/ to clos v0
  } else if (trainConfig == 441) { // gamma eta & meson rap var
    cuts.AddCutPCM("86010113", "03200009327000008250400000", "0162303500000000"); // |eta| < 0.65, |y| < 0.6
    cuts.AddCutPCM("86010113", "04200009327000008250400000", "0162203500000000"); // |eta| < 0.75, |y| < 0.7
    cuts.AddCutPCM("86010113", "01200009327000008250400000", "0162403500000000"); // |eta| < 0.6, |y| < 0.5
  } else if (trainConfig == 442) { // minR and single pt var
    cuts.AddCutPCM("86010113", "00100009327000008250400000", "0162103500000000"); // minR 2.8
    cuts.AddCutPCM("86010113", "00900009327000008250400000", "0162103500000000"); // minR 7.5
    cuts.AddCutPCM("86010113", "00200079327000008250400000", "0162103500000000"); // single pT 0. GeV/c
    cuts.AddCutPCM("86010113", "00200019327000008250400000", "0162103500000000"); // single pT 0.1 GeV/c
  } else if (trainConfig == 443) { // TPC cluster & edEdx var
    cuts.AddCutPCM("86010113", "00200008327000008250400000", "0162103500000000"); // TPC Cluster 0.35
    cuts.AddCutPCM("86010113", "00200006327000008250400000", "0162103500000000"); // TPC Cluster 0.7
    cuts.AddCutPCM("86010113", "00200009227000008250400000", "0162103500000000"); // edEdx -4,5
    cuts.AddCutPCM("86010113", "00200009627000008250400000", "0162103500000000"); // edEdx -2.5,4
    cuts.AddCutPCM("86010113", "00200009127000008250400000", "0162103500000000"); // edEdx 5,5
  } else if (trainConfig == 444) { //PidEdx Variation
    cuts.AddCutPCM("86010113", "00200009357000008250400000", "0162103500000000"); // PidEdx(2, >0.4GeV)
    cuts.AddCutPCM("86010113", "00200009387300008250400000", "0162103500000000"); // PidEdx(2, >0.4GeV; 1>3.5GeV)
    cuts.AddCutPCM("86010113", "00200009320000008250400000", "0162103500000000"); // PidEdx(1, >0.5GeV)
    cuts.AddCutPCM("86010113", "00200009325000008250400000", "0162103500000000"); // PidEdx(1, >0.3GeV)
  } else if (trainConfig == 445) { //PidEdx Variation 2
    cuts.AddCutPCM("86010113", "00200009327300008250400000", "0162103500000000"); // PidEdx(1, >0.4GeV; -10>3.5GeV)
    cuts.AddCutPCM("86010113", "00200009326000008250400000", "0162103500000000"); // PidEdx(1, >0.25GeV)
    cuts.AddCutPCM("86010113", "00200009326200008250400000", "0162103500000000"); // PidEdx(1, >0.25GeV; -10>4GeV)
    cuts.AddCutPCM("86010113", "00200009327200008250400000", "0162103500000000"); // PidEdx(1, >0.4GeV; -10>4GeV)
  } else if (trainConfig == 446) { // qt & psipair variation
    cuts.AddCutPCM("86010113", "00200009327000003250400000", "0162103500000000"); // qT 0.05 1D
    cuts.AddCutPCM("86010113", "00200009327000009250400000", "0162103500000000"); // qT 0.03 2D
    cuts.AddCutPCM("86010113", "00200009327000002250400000", "0162103500000000"); // qT 0.07 1D
    cuts.AddCutPCM("86010113", "00200009327000008240400000", "0162103500000000"); // Psi Pair: 1D 0.2
  } else if (trainConfig == 447) { // PsiPair Variation
    cuts.AddCutPCM("86010113", "00200009327000008210400000", "0162103500000000"); // Psi Pair: 1D 0.1
    cuts.AddCutPCM("86010113", "00200009327000008150400000", "0162103500000000"); // chi2 50  2D
    cuts.AddCutPCM("86010113", "00200009327000008850400000", "0162103500000000"); // chi2 20  2D
    cuts.AddCutPCM("86010113", "00200009327000008260400000", "0162103500000000"); // psi pair 0.05
  } else if (trainConfig == 448) { // cos point & meson alpha variation
    cuts.AddCutPCM("86010113", "00200009327000008250300000", "0162103500000000"); // cos pointing angle 0.75
    cuts.AddCutPCM("86010113", "00200009327000008250600000", "0162103500000000"); // cos pointing angle 0.9
    cuts.AddCutPCM("86010113", "00200009327000008250400000", "0162106500000000"); // alpha meson cut 0.8
  } else if (trainConfig == 449) { // BG mixing & smear variations
    cuts.AddCutPCM("86010113", "00200009327000008250400000", "0262103500000000"); // BG track multiplicity
    cuts.AddCutPCM("86010113", "00200009327000008250400000", "0162103500800000"); // fPSigSmearingCte=0.014;
    cuts.AddCutPCM("86010113", "00200009327000008250400000", "0162103500900000"); // fPSigSmearingCte=0.014;


  //Systematics HM
  } else if (trainConfig == 800) {//MB
    cuts.AddCutPCM("80010113", "00200009a27000008250a04120", "0162103500000000"); // default
    cuts.AddCutPCM("80010113", "00200079a27000008250a04120", "0162103500000000"); // min pT no cut
    cuts.AddCutPCM("80010113", "00200069a27000008250a04120", "0162103500000000"); // min pT 40 cut
    cuts.AddCutPCM("80010113", "00200049a27000008250a04120", "0162103500000000"); // min pT 75 cut
    cuts.AddCutPCM("80010113", "00200008a27000008250a04120", "0162103500000000"); // TPC cluster 35%
    cuts.AddCutPCM("80010113", "00200007a27000008250a04120", "0162103500000000"); // TPC cluster 70%
    cuts.AddCutPCM("80010113", "00200009b27000008250a04120", "0162103500000000"); // edEdx -3.2,3.2
    cuts.AddCutPCM("80010113", "00200009c27000008250a04120", "0162103500000000"); // edEdx -2.8,2.8
 } else if (trainConfig == 801) {//0-10
    cuts.AddCutPCM("80110113", "00200009a27000008250a04120", "0162103500000000"); // default
    cuts.AddCutPCM("80110113", "00200079a27000008250a04120", "0162103500000000"); // min pT no cut
    cuts.AddCutPCM("80110113", "00200069a27000008250a04120", "0162103500000000"); // min pT 40 cut
    cuts.AddCutPCM("80110113", "00200049a27000008250a04120", "0162103500000000"); // min pT 75 cut
    cuts.AddCutPCM("80110113", "00200008a27000008250a04120", "0162103500000000"); // TPC cluster 35%
    cuts.AddCutPCM("80110113", "00200006a27000008250a04120", "0162103500000000"); // TPC cluster 70%
    cuts.AddCutPCM("80110113", "00200009b27000008250a04120", "0162103500000000"); // edEdx -3,2,3.2
    cuts.AddCutPCM("80110113", "00200009c27000008250a04120", "0162103500000000"); // edEdx -2.8,2.8
  } else if (trainConfig == 802) {//0-20
    cuts.AddCutPCM("80210113", "00200009a27000008250a04120", "0162103500000000"); // default
    cuts.AddCutPCM("80210113", "00200079a27000008250a04120", "0162103500000000"); // min pT no cut
    cuts.AddCutPCM("80210113", "00200069a27000008250a04120", "0162103500000000"); // min pT 40 cut
    cuts.AddCutPCM("80210113", "00200049a27000008250a04120", "0162103500000000"); // min pT 75 cut
    cuts.AddCutPCM("80210113", "00200008a27000008250a04120", "0162103500000000"); // TPC cluster 35%
    cuts.AddCutPCM("80210113", "00200006a27000008250a04120", "0162103500000000"); // TPC cluster 70%
    cuts.AddCutPCM("80210113", "00200009b27000008250a04120", "0162103500000000"); // edEdx -3,2,3.2
    cuts.AddCutPCM("80210113", "00200009c27000008250a04120", "0162103500000000"); // edEdx -2.8,2.8
  } else if (trainConfig == 803) {//20-40
    cuts.AddCutPCM("82410113", "00200009a27000008250a04120", "0162103500000000"); // default
    cuts.AddCutPCM("82410113", "00200079a27000008250a04120", "0162103500000000"); // min pT no cut
    cuts.AddCutPCM("82410113", "00200069a27000008250a04120", "0162103500000000"); // min pT 40 cut
    cuts.AddCutPCM("82410113", "00200049a27000008250a04120", "0162103500000000"); // min pT 75 cut
    cuts.AddCutPCM("82410113", "00200008a27000008250a04120", "0162103500000000"); // TPC cluster 35%
    cuts.AddCutPCM("82410113", "00200006a27000008250a04120", "0162103500000000"); // TPC cluster 70%
    cuts.AddCutPCM("82410113", "00200009b27000008250a04120", "0162103500000000"); // edEdx -3.2,3.2
    cuts.AddCutPCM("82410113", "00200009c27000008250a04120", "0162103500000000"); // edEdx -2.8,2.8
  } else if (trainConfig == 804) {//40-60
    cuts.AddCutPCM("84610113", "00200009a27000008250a04120", "0162103500000000"); // default
    cuts.AddCutPCM("84610113", "00200079a27000008250a04120", "0162103500000000"); // min pT no cut
    cuts.AddCutPCM("84610113", "00200069a27000008250a04120", "0162103500000000"); // min pT 40 cut
    cuts.AddCutPCM("84610113", "00200049a27000008250a04120", "0162103500000000"); // min pT 75 cut
    cuts.AddCutPCM("84610113", "00200008a27000008250a04120", "0162103500000000"); // TPC cluster 35%
    cuts.AddCutPCM("84610113", "00200006a27000008250a04120", "0162103500000000"); // TPC cluster 70%
    cuts.AddCutPCM("84610113", "00200009b27000008250a04120", "0162103500000000"); // edEdx -3.2,3.2
    cuts.AddCutPCM("84610113", "00200009c27000008250a04120", "0162103500000000"); // edEdx -2.8,2.8
  } else if (trainConfig == 805) {//60-80
    cuts.AddCutPCM("86810113", "00200009a27000008250a04120", "0162103500000000"); // default
    cuts.AddCutPCM("86810113", "00200079a27000008250a04120", "0162103500000000"); // min pT no cut
    cuts.AddCutPCM("86810113", "00200069a27000008250a04120", "0162103500000000"); // min pT 40 cut
    cuts.AddCutPCM("86810113", "00200049a27000008250a04120", "0162103500000000"); // min pT 75 cut
    cuts.AddCutPCM("86810113", "00200008a27000008250a04120", "0162103500000000"); // TPC cluster 35%
    cuts.AddCutPCM("86810113", "00200006a27000008250a04120", "0162103500000000"); // TPC cluster 70%
    cuts.AddCutPCM("86810113", "00200009b27000008250a04120", "0162103500000000"); // edEdx -3.2,3.2
    cuts.AddCutPCM("86810113", "00200009c27000008250a04120", "0162103500000000"); // edEdx -2.8,2.8
  } else if (trainConfig == 806) {//80-100
    cuts.AddCutPCM("88010113", "00200009a27000008250a04120", "0162103500000000"); // default
    cuts.AddCutPCM("88010113", "00200079a27000008250a04120", "0162103500000000"); // min pT no cut
    cuts.AddCutPCM("88010113", "00200069a27000008250a04120", "0162103500000000"); // min pT 40 cut
    cuts.AddCutPCM("88010113", "00200049a27000008250a04120", "0162103500000000"); // min pT 75 cut
    cuts.AddCutPCM("88010113", "00200008a27000008250a04120", "0162103500000000"); // TPC cluster 35%
    cuts.AddCutPCM("88010113", "00200006a27000008250a04120", "0162103500000000"); // TPC cluster 70%
    cuts.AddCutPCM("88010113", "00200009b27000008250a04120", "0162103500000000"); // edEdx -3.2,3.2
    cuts.AddCutPCM("88010113", "00200009c27000008250a04120", "0162103500000000"); // edEdx -2.8,2.8
  } else if (trainConfig == 807) {//MB
    cuts.AddCutPCM("80010113", "00200009a57000008250a04120", "0162103500000000"); // pidEdx 2,-10
    cuts.AddCutPCM("80010113", "00200009a17000008250a04120", "0162103500000000"); // pidEdx 0,-10
    cuts.AddCutPCM("80010113", "00200009a27000003250a04120", "0162103500000000"); // qT max 0.05 1D
    cuts.AddCutPCM("80010113", "00200009a27000002250a04120", "0162103500000000"); // qT max 0.06 2D
    cuts.AddCutPCM("80010113", "00200009a27000008250904120", "0162103500000000"); // CosPA 0.99
    cuts.AddCutPCM("80010113", "00200009a27000008250b04120", "0162103500000000"); // CosPA 0.985
  } else if (trainConfig == 808) {//0-10
    cuts.AddCutPCM("80110113", "00200009a57000008250a04120", "0162103500000000"); // pidEdx 2,-10
    cuts.AddCutPCM("80110113", "00200009a17000008250a04120", "0162103500000000"); // pidEdx 0,-10
    cuts.AddCutPCM("80110113", "00200009a27000003250a04120", "0162103500000000"); // qT max 0.05 1D
    cuts.AddCutPCM("80110113", "00200009a27000002250a04120", "0162103500000000"); // qT max 0.06 2D
    cuts.AddCutPCM("80110113", "00200009a27000008250904120", "0162103500000000"); // CosPA 0.99
    cuts.AddCutPCM("80110113", "00200009a27000008250b04120", "0162103500000000"); // CosPA 0.985
  } else if (trainConfig == 809) {//0-20
    cuts.AddCutPCM("80210113", "00200009a57000008250a04120", "0162103500000000"); // pidEdx 2,-10
    cuts.AddCutPCM("80210113", "00200009a17000008250a04120", "0162103500000000"); // pidEdx 0,-10
    cuts.AddCutPCM("80210113", "00200009a27000003250a04120", "0162103500000000"); // qT max 0.05 1D
    cuts.AddCutPCM("80210113", "00200009a27000002250a04120", "0162103500000000"); // qT max 0.06 2D
    cuts.AddCutPCM("80210113", "00200009a27000008250904120", "0162103500000000"); // CosPA 0.99
    cuts.AddCutPCM("80210113", "00200009a27000008250b04120", "0162103500000000"); // CosPA 0.985
  } else if (trainConfig == 810) {//20-40
    cuts.AddCutPCM("82410113", "00200009a57000008250a04120", "0162103500000000"); // pidEdx 2,-10
    cuts.AddCutPCM("82410113", "00200009a17000008250a04120", "0162103500000000"); // pidEdx 0,-10
    cuts.AddCutPCM("82410113", "00200009a27000003250a04120", "0162103500000000"); // qT max 0.05 1D
    cuts.AddCutPCM("82410113", "00200009a27000002250a04120", "0162103500000000"); // qT max 0.06 2D
    cuts.AddCutPCM("82410113", "00200009a27000008250904120", "0162103500000000"); // CosPA 0.99
    cuts.AddCutPCM("82410113", "00200009a27000008250b04120", "0162103500000000"); // CosPA 0.985
  } else if (trainConfig == 811) {//40-60
    cuts.AddCutPCM("84610113", "00200009a57000008250a04120", "0162103500000000"); // pidEdx 2,-10
    cuts.AddCutPCM("84610113", "00200009a17000008250a04120", "0162103500000000"); // pidEdx 0,-10
    cuts.AddCutPCM("84610113", "00200009a27000003250a04120", "0162103500000000"); // qT max 0.05 1D
    cuts.AddCutPCM("84610113", "00200009a27000002250a04120", "0162103500000000"); // qT max 0.06 2D
    cuts.AddCutPCM("84610113", "00200009a27000008250904120", "0162103500000000"); // CosPA 0.99
    cuts.AddCutPCM("84610113", "00200009a27000008250b04120", "0162103500000000"); // CosPA 0.985
  } else if (trainConfig == 812) {//60-80
    cuts.AddCutPCM("86810113", "00200009a27000008250a04120", "0162103500000000"); // pidEdx 2,-10
    cuts.AddCutPCM("86810113", "00200009a27000008250a04120", "0162103500000000"); // pidEdx 0,-10
    cuts.AddCutPCM("86810113", "00200009a27000003250a04120", "0162103500000000"); // qT max 0.05 1D
    cuts.AddCutPCM("86810113", "00200009a27000002250a04120", "0162103500000000"); // qT max 0.06 2D
    cuts.AddCutPCM("86810113", "00200009a27000008250904120", "0162103500000000"); // CosPA 0.99
    cuts.AddCutPCM("86810113", "00200009a27000008250b04120", "0162103500000000"); // CosPA 0.985
  } else if (trainConfig == 813) {//80-100
    cuts.AddCutPCM("88010113", "00200009a27000008250a04120", "0162103500000000"); // pidEdx 2,-10
    cuts.AddCutPCM("88010113", "00200009a27000008250a04120", "0162103500000000"); // pidEdx 0,-10
    cuts.AddCutPCM("88010113", "00200009a27000003250a04120", "0162103500000000"); // qT max 0.05 1D
    cuts.AddCutPCM("88010113", "00200009a27000002250a04120", "0162103500000000"); // qT max 0.06 2D
    cuts.AddCutPCM("88010113", "00200009a27000008250904120", "0162103500000000"); // CosPA 0.99
    cuts.AddCutPCM("88010113", "00200009a27000008250b04120", "0162103500000000"); // CosPA 0.985
  } else if (trainConfig == 814) {//MB
    cuts.AddCutPCM("80010113", "00200009a27000008250a04120", "0162103500000000"); // default
    cuts.AddCutPCM("80010113", "00200009a27000008250c04120", "0162103500000000"); // cosPA 0.996
    cuts.AddCutPCM("80010113", "00200009a27000008250d04120", "0162103500000000"); // cosPA 0.997
    cuts.AddCutPCM("80010113", "00200009a27000008250e04120", "0162103500000000"); // cosPA 0.998
    cuts.AddCutPCM("80010113", "00200009a27000008250f04120", "0162103500000000"); // cosPA 0.999
  } else if (trainConfig == 815) {//0-10
    cuts.AddCutPCM("80110113", "00200009a27000008250a04120", "0162103500000000"); // default
    cuts.AddCutPCM("80110113", "00200009a27000008250c04120", "0162103500000000"); // cosPA 0.996
    cuts.AddCutPCM("80110113", "00200009a27000008250d04120", "0162103500000000"); // cosPA 0.997
    cuts.AddCutPCM("80110113", "00200009a27000008250e04120", "0162103500000000"); // cosPA 0.998
    cuts.AddCutPCM("80110113", "00200009a27000008250f04120", "0162103500000000"); // cosPA 0.999
  } else if (trainConfig == 816) {//0-20
    cuts.AddCutPCM("80210113", "00200009a27000008250a04120", "0162103500000000"); // default
    cuts.AddCutPCM("80210113", "00200009a27000008250c04120", "0162103500000000"); // cosPA 0.996
    cuts.AddCutPCM("80210113", "00200009a27000008250d04120", "0162103500000000"); // cosPA 0.997
    cuts.AddCutPCM("80210113", "00200009a27000008250e04120", "0162103500000000"); // cosPA 0.998
    cuts.AddCutPCM("80210113", "00200009a27000008250f04120", "0162103500000000"); // cosPA 0.999
  } else if (trainConfig == 817) {//20-40
    cuts.AddCutPCM("82410113", "00200009a27000008250a04120", "0162103500000000"); // default
    cuts.AddCutPCM("82410113", "00200009a27000008250c04120", "0162103500000000"); // cosPA 0.996
    cuts.AddCutPCM("82410113", "00200009a27000008250d04120", "0162103500000000"); // cosPA 0.997
    cuts.AddCutPCM("82410113", "00200009a27000008250e04120", "0162103500000000"); // cosPA 0.998
    cuts.AddCutPCM("82410113", "00200009a27000008250f04120", "0162103500000000"); // cosPA 0.999
  } else if (trainConfig == 818) {//40-60
    cuts.AddCutPCM("84610113", "00200009a27000008250a04120", "0162103500000000"); // default
    cuts.AddCutPCM("84610113", "00200009a27000008250c04120", "0162103500000000"); // cosPA 0.996
    cuts.AddCutPCM("84610113", "00200009a27000008250d04120", "0162103500000000"); // cosPA 0.997
    cuts.AddCutPCM("84610113", "00200009a27000008250e04120", "0162103500000000"); // cosPA 0.998
    cuts.AddCutPCM("84610113", "00200009a27000008250f04120", "0162103500000000"); // cosPA 0.999
  } else if (trainConfig == 819) {//60-80
    cuts.AddCutPCM("86810113", "00200009a27000008250a04120", "0162103500000000"); // default
    cuts.AddCutPCM("86810113", "00200009a27000008250c04120", "0162103500000000"); // cosPA 0.996
    cuts.AddCutPCM("86810113", "00200009a27000008250d04120", "0162103500000000"); // cosPA 0.997
    cuts.AddCutPCM("86810113", "00200009a27000008250e04120", "0162103500000000"); // cosPA 0.998
    cuts.AddCutPCM("86810113", "00200009a27000008250f04120", "0162103500000000"); // cosPA 0.999
  } else if (trainConfig == 820) {//80-100
    cuts.AddCutPCM("88010113", "00200009a27000008250a04120", "0162103500000000"); // default
    cuts.AddCutPCM("88010113", "00200009a27000008250c04120", "0162103500000000"); // cosPA 0.996
    cuts.AddCutPCM("88010113", "00200009a27000008250d04120", "0162103500000000"); // cosPA 0.997
    cuts.AddCutPCM("88010113", "00200009a27000008250e04120", "0162103500000000"); // cosPA 0.998
    cuts.AddCutPCM("88010113", "00200009a27000008250f04120", "0162103500000000"); // cosPA 0.999


  // OLD CONFIGURATIONS WITH DCA TREE
  } else if (trainConfig == 1001){
    cuts.AddCutPCM("80010113", "00200009f9730000dge0400000", "0162103500900000"); // new standard pPb MB
  } else if (trainConfig == 1002) {
    cuts.AddCutPCM("80010113", "00200009327000008250400000", "0162103500900000"); // new standard pPb MB
  } else if (trainConfig == 1003) {
    cuts.AddCutPCM("80210113", "00200009f9730000dge0400000", "0162103500900000"); // new standard pPb 0-20
  } else if (trainConfig == 1004) {
    cuts.AddCutPCM("80210113", "00200009327000008250400000", "0162103500900000"); // new standard pPb 0-20
  } else if (trainConfig == 1005) {
    cuts.AddCutPCM("82410113", "00200009f9730000dge0400000", "0162103500900000"); // new standard pPb 20-40
  } else if (trainConfig == 1006) {
    cuts.AddCutPCM("82410113", "00200009327000008250400000", "0162103500900000"); // new standard pPb 20-40
  } else if (trainConfig == 1007) {
    cuts.AddCutPCM("84610113", "00200009f9730000dge0400000", "0162103500900000"); // new standard pPb 40-60
  } else if (trainConfig == 1008) {
    cuts.AddCutPCM("84610113", "00200009327000008250400000", "0162103500900000"); // new standard pPb 40-60
  } else if (trainConfig == 1009) {
    cuts.AddCutPCM("86010113", "00200009f9730000dge0400000", "0162103500900000"); // new standard pPb 60-80
  } else if (trainConfig == 1010) {
    cuts.AddCutPCM("86010113", "00200009327000008250400000", "0162103500900000"); // new standard pPb 60-80

  // configurations with past future protection (2.25 \mus protected)
  } else if (trainConfig == 1011){
    cuts.AddCutPCM("80010213", "00200009f9730000dge0400000", "0162103500900000"); // new standard pPb MB
  } else if (trainConfig == 1012) {
    cuts.AddCutPCM("80010213", "00200009327000008250400000", "0162103500900000"); // new standard pPb MB
  } else if (trainConfig == 1013) {
    cuts.AddCutPCM("80210213", "00200009f9730000dge0400000", "0162103500900000"); // new standard pPb 0-20
  } else if (trainConfig == 1014) {
    cuts.AddCutPCM("80210213", "00200009327000008250400000", "0162103500900000"); // new standard pPb 0-20
  } else if (trainConfig == 1015) {
    cuts.AddCutPCM("82410213", "00200009f9730000dge0400000", "0162103500900000"); // new standard pPb 20-40
  } else if (trainConfig == 1016) {
    cuts.AddCutPCM("82410213", "00200009327000008250400000", "0162103500900000"); // new standard pPb 20-40
  } else if (trainConfig == 1017) {
    cuts.AddCutPCM("84610213", "00200009f9730000dge0400000", "0162103500900000"); // new standard pPb 40-60
  } else if (trainConfig == 1018) {
    cuts.AddCutPCM("84610213", "00200009327000008250400000", "0162103500900000"); // new standard pPb 40-60
  } else if (trainConfig == 1019) {
    cuts.AddCutPCM("86010213", "00200009f9730000dge0400000", "0162103500900000"); // new standard pPb 60-80
  } else if (trainConfig == 1020) {
    cuts.AddCutPCM("86010213", "00200009327000008250400000", "0162103500900000"); // new standard pPb 60-80

  // configurations for eta cuts
  } else if (trainConfig == 1021) {
    cuts.AddCutPCM("80010113", "0a200009327000008250400000", "0162103500000000"); //Eta cut -0.9 - -0.2 and 0.2 - 0.9
  } else if (trainConfig == 1022) {
    cuts.AddCutPCM("80010113", "0b200009327000008250400000", "0162103500000000"); //Eta cut -0.9 - -0.2 and 0.2 - 0.9 with LineCut

  // ZNA and CL1 configs
  } else if (trainConfig == 1030) {
    cuts.AddCutPCM("90010113", "00200009327000008250400000", "0162103500900000");
  } else if (trainConfig == 1031) {
    cuts.AddCutPCM("90210113", "00200009327000008250400000", "0162103500900000");
  } else if (trainConfig == 1032) {
    cuts.AddCutPCM("92410113", "00200009327000008250400000", "0162103500900000");
  } else if (trainConfig == 1033) {
    cuts.AddCutPCM("94610113", "00200009327000008250400000", "0162103500900000");
  } else if (trainConfig == 1034) {
    cuts.AddCutPCM("96010113", "00200009327000008250400000", "0162103500900000");

  } else if (trainConfig == 1035) {
    cuts.AddCutPCM("e0010113", "00200009327000008250400000", "0162103500900000");
  } else if (trainConfig == 1036) {
    cuts.AddCutPCM("e0210113", "00200009327000008250400000", "0162103500900000");
  } else if (trainConfig == 1037) {
    cuts.AddCutPCM("e2410113", "00200009327000008250400000", "0162103500900000");
  } else if (trainConfig == 1038) {
    cuts.AddCutPCM("e4610113", "00200009327000008250400000", "0162103500900000");
  } else if (trainConfig == 1039) {
    cuts.AddCutPCM("e6010113", "00200009327000008250400000", "0162103500900000");

  //Run 2 pPb
  } else if (trainConfig == 1100) {
    cuts.AddCutPCM("80010113", "00200009f9730000dge0400000", "0162103500000000"); // new default for 5TeV
  } else if (trainConfig == 1101) {
    cuts.AddCutPCM("80110113", "00200009f9730000dge0400000", "0162103500000000"); // 0-10
  } else if (trainConfig == 1102) {
    cuts.AddCutPCM("81210113", "00200009f9730000dge0400000", "0162103500000000"); // 0-20
  } else if (trainConfig == 1103) {
    cuts.AddCutPCM("82410113", "00200009f9730000dge0400000", "0162103500000000"); // 20-40
  } else if (trainConfig == 1104) {
    cuts.AddCutPCM("84610113", "00200009f9730000dge0400000", "0162103500000000"); // 40-60
  } else if (trainConfig == 1105) {
    cuts.AddCutPCM("86810113", "00200009f9730000dge0400000", "0162103500000000"); // 60-80
  } else if (trainConfig == 1106) {
    cuts.AddCutPCM("88010113", "00200009f9730000dge0400000", "0162103500000000"); // 80-100
  } else if (trainConfig == 1107) {
    cuts.AddCutPCM("80210113", "00200009f9730000dge0400000", "0162103500000000"); // 0-20
  } else if (trainConfig == 1108) {
    cuts.AddCutPCM("86010113", "00200009f9730000dge0400000", "0162103500000000"); // 60-100
  } else if (trainConfig == 1109) {
    cuts.AddCutPCM("a0110113", "00200009f9730000dge0400000", "0162103500000000"); // 0-5
  } else if (trainConfig == 1110) {
    cuts.AddCutPCM("a1210113", "00200009f9730000dge0400000", "0162103500000000"); // 5-10
  } else if (trainConfig == 1111) {
    cuts.AddCutPCM("c0110113", "00200009f9730000dge0400000", "0162103500000000"); // 0-1
  } else if (trainConfig == 1112) {
    cuts.AddCutPCM("c0210113", "00200009f9730000dge0400000", "0162103500000000"); // 0-2

  } else if (trainConfig == 1120) {
    cuts.AddCutPCM("80010113", "0dm00009f9730000dge0404000", "0152103500000000"); // new default (R region rej. + eta<0.8 + DC)
  } else if (trainConfig == 1121) { // TOF single leg cut
    cuts.AddCutPCM("80010113", "0dm00009f9730600dge0404000", "0152103500000000"); // new default (R region rej. + eta<0.8 + DC)
  } else if (trainConfig == 1122) { // TOF both leg cut
    cuts.AddCutPCM("80010113", "0dm00009f9730700dge0404000", "0152103500000000"); // new default (R region rej. + eta<0.8 + DC)
  } else if (trainConfig == 1123) { // TOF single leg cut
    cuts.AddCutPCM("80010113", "0dm00009f9730800dge0404000", "0152103500000000"); // new default (R region rej. + eta<0.8 + DC)
  } else if (trainConfig == 1124) { // TOF both leg cut
    cuts.AddCutPCM("80010113", "0dm00009f9730900dge0404000", "0152103500000000"); // new default (R region rej. + eta<0.8 + DC)

  } else if (trainConfig == 1125) { // T0-based cuts
    cuts.AddCutPCM("80011103", "0dm00009f9730000dge0404000", "0152103500000000"); // new default (R region rej. + eta<0.8 + DC)

  } else if (trainConfig == 1130) {
    cuts.AddCutPCM("80210103", "0dm00009f9730000dge0404000", "0152103500000000"); // 00-20
    cuts.AddCutPCM("82410103", "0dm00009f9730000dge0404000", "0152103500000000"); // 20-40
    cuts.AddCutPCM("84610103", "0dm00009f9730000dge0404000", "0152103500000000"); // 40-60
  } else if (trainConfig == 1131) {
    cuts.AddCutPCM("86810103", "0dm00009f9730000dge0404000", "0152103500000000"); // 60-80
    cuts.AddCutPCM("88010103", "0dm00009f9730000dge0404000", "0152103500000000"); // 80-100
    cuts.AddCutPCM("89010103", "0dm00009f9730000dge0404000", "0152103500000000"); // 90-100

  } else if (trainConfig == 1150) {
    cuts.AddCutPCM("80010123", "00200009f9730000dge0400000", "0162103500000000", "4117901050032230000"); // new default for 8TeV+triggers
  } else if (trainConfig == 1151) {
    cuts.AddCutPCM("8008e123", "00200009f9730000dge0400000", "0162103500000000", "4117901050032230000"); // new default for 8TeV+triggers
    cuts.AddCutPCM("8008d123", "00200009f9730000dge0400000", "0162103500000000", "4117901050032230000"); // new default for 8TeV+triggers

  // triggers EMC7
  } else if (trainConfig == 1200) {
    cuts.AddCutPCM("80052113", "00200009327000008250400000", "0162103500900000", "1111100007032230000"); // new standard pPb MB
  } else if (trainConfig == 1201) {
    cuts.AddCutPCM("80252113", "00200009327000008250400000", "0162103500900000", "1111100007032230000"); // new standard pPb 0-20
  } else if (trainConfig == 1202) {
    cuts.AddCutPCM("82452113", "00200009327000008250400000", "0162103500900000", "1111100007032230000"); // new standard pPb 20-40
  } else if (trainConfig == 1203) {
    cuts.AddCutPCM("84652113", "00200009327000008250400000", "0162103500900000", "1111100007032230000"); // new standard pPb 40-60
  } else if (trainConfig == 1204) {
    cuts.AddCutPCM("86052113", "00200009327000008250400000", "0162103500900000", "1111100007032230000"); // new standard pPb 60-100
  // triggers EG2
  } else if (trainConfig == 1210) {
    cuts.AddCutPCM("80085113", "00200009327000008250400000", "0162103500900000", "1111100007032230000"); // new standard pPb MB
  } else if (trainConfig == 1211) {
    cuts.AddCutPCM("80285113", "00200009327000008250400000", "0162103500900000", "1111100007032230000"); // new standard pPb 0-20
  } else if (trainConfig == 1212) {
    cuts.AddCutPCM("82485113", "00200009327000008250400000", "0162103500900000", "1111100007032230000"); // new standard pPb 20-40
  } else if (trainConfig == 1213) {
    cuts.AddCutPCM("84685113", "00200009327000008250400000", "0162103500900000", "1111100007032230000"); // new standard pPb 40-60
  } else if (trainConfig == 1214) {
    cuts.AddCutPCM("86085113", "00200009327000008250400000", "0162103500900000", "1111100007032230000"); // new standard pPb 60-100
  // triggers EG1
  } else if (trainConfig == 1220) {
    cuts.AddCutPCM("80083113", "00200009327000008250400000", "0162103500900000", "1111100007032230000"); // new standard pPb MB
  } else if (trainConfig == 1221) {
    cuts.AddCutPCM("80283113", "00200009327000008250400000", "0162103500900000", "1111100007032230000"); // new standard pPb 0-20
  } else if (trainConfig == 1222) {
    cuts.AddCutPCM("82483113", "00200009327000008250400000", "0162103500900000", "1111100007032230000"); // new standard pPb 20-40
  } else if (trainConfig == 1223) {
    cuts.AddCutPCM("84683113", "00200009327000008250400000", "0162103500900000", "1111100007032230000"); // new standard pPb 40-60
  } else if (trainConfig == 1224) {
    cuts.AddCutPCM("86083113", "00200009327000008250400000", "0162103500900000", "1111100007032230000"); // new standard pPb 60-100

  //Run 2 pPb HM
  } else if (trainConfig == 1500) {
    cuts.AddCutPCM("80010113", "00200009a27000008250a04120", "0162103500000000"); // test 1 for 5TeV
  } else if (trainConfig == 1501) {
    cuts.AddCutPCM("80110113", "00200009a27000008250a04120", "0162103500000000"); // 0-10
  } else if (trainConfig == 1502) {
    cuts.AddCutPCM("80210113", "00200009a27000008250a04120", "0162103500000000"); // 0-20
  } else if (trainConfig == 1503) {
    cuts.AddCutPCM("82410113", "00200009a27000008250a04120", "0162103500000000"); // 20-40
  } else if (trainConfig == 1504) {
    cuts.AddCutPCM("84610113", "00200009a27000008250a04120", "0162103500000000"); // 40-60
  } else if (trainConfig == 1505) {
    cuts.AddCutPCM("86810113", "00200009a27000008250a04120", "0162103500000000"); // 60-80
  } else if (trainConfig == 1506) {
    cuts.AddCutPCM("88010113", "00200009a27000008250a04120", "0162103500000000"); // 80-100
  } else if (trainConfig == 1507) {
    cuts.AddCutPCM("80010113", "0d200009a27000008250a04120", "0162103500000000"); // test 2 for 5TeV
  } else if (trainConfig == 1508) {
    cuts.AddCutPCM("80110113", "0d200009a27000008250a04120", "0162103500000000"); // 0-10
  } else if (trainConfig == 1509) {
    cuts.AddCutPCM("80210113", "0d200009a27000008250a04120", "0162103500000000"); // 0-20
  } else if (trainConfig == 1510) {
    cuts.AddCutPCM("82410113", "0d200009a27000008250a04120", "0162103500000000"); // 20-40
  } else if (trainConfig == 1511) {
    cuts.AddCutPCM("84610113", "0d200009a27000008250a04120", "0162103500000000"); // 40-60
  } else if (trainConfig == 1512) {
    cuts.AddCutPCM("86810113", "0d200009a27000008250a04120", "0162103500000000"); // 60-80
  } else if (trainConfig == 1513) {
    cuts.AddCutPCM("88010113", "0d200009a27000008250a04120", "0162103500000000"); // 80-100
  } else if (trainConfig == 1514) {
    cuts.AddCutPCM("90010113", "0d200009a27000008250a04120", "0162103500000000"); // CL1
  } else if (trainConfig == 1515) {
    cuts.AddCutPCM("90110113", "0d200009a27000008250a04120", "0162103500000000"); // 0-10
  } else if (trainConfig == 1516) {
    cuts.AddCutPCM("90210113", "0d200009a27000008250a04120", "0162103500000000"); // 0-20
  } else if (trainConfig == 1517) {
    cuts.AddCutPCM("92410113", "0d200009a27000008250a04120", "0162103500000000"); // 20-40
  } else if (trainConfig == 1518) {
    cuts.AddCutPCM("94610113", "0d200009a27000008250a04120", "0162103500000000"); // 40-60
  } else if (trainConfig == 1519) {
    cuts.AddCutPCM("96810113", "0d200009a27000008250a04120", "0162103500000000"); // 60-80
  } else if (trainConfig == 1520) {
    cuts.AddCutPCM("98010113", "0d200009a27000008250a04120", "0162103500000000"); // 80-100
  } else if (trainConfig == 1521) {
    cuts.AddCutPCM("e0010113", "0d200009a27000008250a04120", "0162103500000000"); // ZNA
  } else if (trainConfig == 1522) {
    cuts.AddCutPCM("e0110113", "0d200009a27000008250a04120", "0162103500000000"); // 0-10
  } else if (trainConfig == 1523) {
    cuts.AddCutPCM("e0210113", "0d200009a27000008250a04120", "0162103500000000"); // 0-20
  } else if (trainConfig == 1524) {
    cuts.AddCutPCM("e2410113", "0d200009a27000008250a04120", "0162103500000000"); // 20-40
  } else if (trainConfig == 1525) {
    cuts.AddCutPCM("e4610113", "0d200009a27000008250a04120", "0162103500000000"); // 40-60
  } else if (trainConfig == 1526) {
    cuts.AddCutPCM("e6810113", "0d200009a27000008250a04120", "0162103500000000"); // 60-80
  } else if (trainConfig == 1527) {
    cuts.AddCutPCM("e8010113", "0d200009a27000008250a04120", "0162103500000000"); // 80-100

  //--------------------------------------------------------------------------
  // Configurations for Jet analysis for pPb 5.02 TeV
  //--------------------------------------------------------------------------
  } else if (trainConfig == 1600) {
    cuts.AddCutPCM("80010113", "00200009f9730000dge0400000", "0162103500000000"); // new default for 5TeV
    cuts.AddCutPCM("80010113", "00200009f9730000dge0400000", "2162103500000000"); // Injet 
    cuts.AddCutPCM("80010113", "00200009f9730000dge0400000", "3162103500000000"); // Injet with JetQA


  } else {
    Error(Form("GammaConvV1_%i",trainConfig), "wrong trainConfig variable no cuts have been specified for the configuration");
    return;
  }

  if(!cuts.AreValid()){
    cout << "\n\n****************************************************" << endl;
    cout << "ERROR: No valid cuts stored in CutHandlerConv! Returning..." << endl;
    cout << "****************************************************\n\n" << endl;
    return;
  }

  Int_t numberOfCuts = cuts.GetNCuts();

  TList *EventCutList = new TList();
  TList *ConvCutList = new TList();
  TList *MesonCutList = new TList();
  TList *ClusterCutList = new TList();

  TList *HeaderList = new TList();
  if (doWeightingPart==1) {
    TObjString *Header1 = new TObjString("pi0_1");
    HeaderList->Add(Header1);
  }
  if (doWeightingPart==2){
    TObjString *Header3 = new TObjString("eta_2");
    HeaderList->Add(Header3);
  }
  if (doWeightingPart==3) {
    TObjString *Header1 = new TObjString("pi0_1");
    HeaderList->Add(Header1);
    TObjString *Header3 = new TObjString("eta_2");
    HeaderList->Add(Header3);
  }
  if (periodNameV0Reader.Contains("LHC18b9")||periodNameV0Reader.Contains("LHC17g8")){
    TObjString *HeaderPMB = new TObjString("EPOSLHC_0");
    TObjString *HeaderP8J = new TObjString("Pythia8Jets_1");
    if (doWeightingPart==4) { // all headers
      HeaderList->Add(HeaderPMB);
      HeaderList->Add(HeaderP8J);
    } else if (doWeightingPart==5) { // only MB header
      HeaderList->Add(HeaderPMB);
    } else { // only JJ header
      HeaderList->Add(HeaderP8J);
    }
  }

  Bool_t doWeighting = kFALSE;
  if (doWeightingPart == 1 || doWeightingPart == 2 || doWeightingPart == 3) doWeighting = kTRUE;

  EventCutList->SetOwner(kTRUE);
  AliConvEventCuts **analysisEventCuts = new AliConvEventCuts*[numberOfCuts];
  ConvCutList->SetOwner(kTRUE);
  AliConversionPhotonCuts **analysisCuts = new AliConversionPhotonCuts*[numberOfCuts];
  MesonCutList->SetOwner(kTRUE);
  AliConversionMesonCuts **analysisMesonCuts = new AliConversionMesonCuts*[numberOfCuts];
  ClusterCutList->SetOwner(kTRUE);
  AliCaloPhotonCuts **analysisClusterCuts     = new AliCaloPhotonCuts*[numberOfCuts];
  Bool_t enableClustersForTrigger             = kFALSE;
  Bool_t initializedMatBudWeigths_existing    = kFALSE;

  if (doWeighting) Printf("weighting has been switched on");

  for(Int_t i = 0; i<numberOfCuts; i++){
    cout << "initialization of cutnumber: " << i << endl;

    analysisEventCuts[i] = new AliConvEventCuts();
    if ( trainConfig == 1 || trainConfig == 2 ||  trainConfig == 5 || trainConfig == 6 ||
        ( trainConfig > 99 && trainConfig < 120 ) ||
        ( trainConfig > 199 && trainConfig < 210 ) || trainConfig == 1001 ||  trainConfig == 1002 ||  trainConfig == 1011 || trainConfig == 1012  ||  trainConfig == 1021 || trainConfig == 1022 ||   trainConfig == 1100){
      if (doWeighting){
        if (generatorName.CompareTo("DPMJET")==0){
          analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE, fileNamePtWeights, "Pi0_DPMJET_LHC13b2_efix_pPb_5023GeV_MBV0A",
                                        "Eta_DPMJET_LHC13b2_efix_pPb_5023GeV_MBV0A", "","Pi0_Fit_Data_pPb_5023GeV_MBV0A","Eta_Fit_Data_pPb_5023GeV_MBV0A");
        } else if (generatorName.CompareTo("HIJING")==0){
          analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE, fileNamePtWeights, "Pi0_Hijing_LHC13e7_pPb_5023GeV_MBV0A",
                                        "Eta_Hijing_LHC13e7_pPb_5023GeV_MBV0A", "","Pi0_Fit_Data_pPb_5023GeV_MBV0A","Eta_Fit_Data_pPb_5023GeV_MBV0A");
        }
      }
    }
    if ( trainConfig == 3 || trainConfig == 4 || trainConfig == 7 || trainConfig == 8 ||
        ( trainConfig > 119 && trainConfig < 140 ) ||
        ( trainConfig > 209 && trainConfig < 220 ) ){
      if (doWeighting){
        analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE, fileNamePtWeights, "Pi0_Hijing_LHC13e7_addSig_pPb_5023GeV_MBV0A",
                                      "Eta_Hijing_LHC13e7_addSig_pPb_5023GeV_MBV0A", "","Pi0_Fit_Data_pPb_5023GeV_MBV0A","Eta_Fit_Data_pPb_5023GeV_MBV0A");
      }
    }

    TString dataInputMultHisto  = "";
    TString mcInputMultHisto    = "";
    TString triggerString       = (cuts.GetEventCut(i)).Data();
    triggerString               = triggerString(3,2);

    dataInputMultHisto          = Form("%s_%s", periodNameAnchor.Data(), triggerString.Data());
    mcInputMultHisto            = Form("%s_%s", periodNameV0Reader.Data(), triggerString.Data());

    if (enableMultiplicityWeighting){
      cout << "enabling mult weighting" << endl;
      analysisEventCuts[i]->SetUseWeightMultiplicityFromFile( kTRUE, fileNameMultWeights, dataInputMultHisto, mcInputMultHisto );
    }

    if (debugLevel > 0) analysisEventCuts[i]->SetDebugLevel(debugLevel);
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
    analysisEventCuts[i]->SetCorrectionTaskSetting(corrTaskSetting);
    analysisEventCuts[i]->SetV0ReaderName(V0ReaderName);
    if (periodNameV0Reader.CompareTo("") != 0) analysisEventCuts[i]->SetPeriodEnum(periodNameV0Reader);
    analysisEventCuts[i]->SetLightOutput(enableLightOutput);
    analysisEventCuts[i]->InitializeCutsFromCutString((cuts.GetEventCut(i)).Data());
    if (doEtaShiftIndCuts) {
      analysisEventCuts[i]->DoEtaShift(doEtaShiftIndCuts);
      analysisEventCuts[i]->SetEtaShift(stringShift);
    }

    EventCutList->Add(analysisEventCuts[i]);
    analysisEventCuts[i]->SetFillCutHistograms("",kFALSE);
    cout << "initialized event cut: " << (cuts.GetEventCut(i)).Data() << endl;

    if(!cuts.GetClusterCut(i).CompareTo("")){
      cout << "\nNo cluster cut set, not filling cluster histograms for triggers ...\n" << endl;
    } else {
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

        enableClustersForTrigger  = kTRUE;
        analysisClusterCuts[i]    = new AliCaloPhotonCuts();
        analysisClusterCuts[i]->SetV0ReaderName(V0ReaderName);
        analysisClusterCuts[i]->SetCorrectionTaskSetting(corrTaskSetting);
        analysisClusterCuts[i]->SetLightOutput(enableLightOutput);
        analysisClusterCuts[i]->SetCaloTrackMatcherName(TrackMatcherName);
        analysisClusterCuts[i]->InitializeCutsFromCutString((cuts.GetClusterCut(i)).Data());
        ClusterCutList->Add(analysisClusterCuts[i]);
        analysisClusterCuts[i]->SetFillCutHistograms("");
    }


    analysisCuts[i] = new AliConversionPhotonCuts();
    if ( trainConfig > 149 && trainConfig < 156 ){
      analysisCuts[i]->SetSwitchToKappaInsteadOfNSigdEdxTPC(kTRUE);
    }

    if (enableMatBudWeightsPi0 > 0){
      if (isMC > 0){
        if (analysisCuts[i]->InitializeMaterialBudgetWeights(enableMatBudWeightsPi0,fileNameMatBudWeights)){
          initializedMatBudWeigths_existing = kTRUE;}
        else {cout << "ERROR The initialization of the materialBudgetWeights did not work out." << endl;}
      } else {cout << "ERROR 'enableMatBudWeightsPi0'-flag was set > 0 even though this is not a MC task. It was automatically reset to 0." << endl;}
    }
    if (enableElecDeDxPostCalibration){
      if (isMC == 0){
        if(fileNamedEdxPostCalib.CompareTo("") != 0){
          analysisCuts[i]->SetElecDeDxPostCalibrationCustomFile(fileNamedEdxPostCalib);
          cout << "Setting custom dEdx recalibration file: " << fileNamedEdxPostCalib.Data() << endl;
        }
        analysisCuts[i]->SetDoElecDeDxPostCalibration(enableElecDeDxPostCalibration);
        cout << "Enabled TPC dEdx recalibration." << endl;
      } else{
        cout << "ERROR enableElecDeDxPostCalibration set to True even if MC file. Automatically reset to 0"<< endl;
        enableElecDeDxPostCalibration=kFALSE;
        analysisCuts[i]->SetDoElecDeDxPostCalibration(kFALSE);
      }
    }
    analysisCuts[i]->SetV0ReaderName(V0ReaderName);
    analysisCuts[i]->SetLightOutput(enableLightOutput);
    analysisCuts[i]->InitializeCutsFromCutString((cuts.GetPhotonCut(i)).Data());

    analysisCuts[i]->SetIsHeavyIon(isHeavyIon);
    if (trainConfig == 1015 || trainConfig==1016 || trainConfig==1017  || trainConfig==1018) {
            analysisCuts[i]->SetDodEdxSigmaCut(kFALSE);
    }
    ConvCutList->Add(analysisCuts[i]);
    analysisCuts[i]->SetFillCutHistograms("",kFALSE);
    cout << "initialized photon cut: " << (cuts.GetPhotonCut(i)).Data() << endl;

    analysisMesonCuts[i] = new AliConversionMesonCuts();
    analysisMesonCuts[i]->SetLightOutput(enableLightOutput);
    if (trainConfig ==1013 || trainConfig ==1014){
      analysisMesonCuts[i]->SetOpeningAngleCut(0.000);
    }
    analysisMesonCuts[i]->InitializeCutsFromCutString((cuts.GetMesonCut(i)).Data());
    MesonCutList->Add(analysisMesonCuts[i]);
    analysisMesonCuts[i]->SetFillCutHistograms("");
    analysisEventCuts[i]->SetAcceptedHeader(HeaderList);
    cout << "initialized meson cut: " << (cuts.GetMesonCut(i)).Data() << endl;
  }

  task->SetDoTHnSparse(enableTHnSparse);
  task->SetEventCutList(numberOfCuts,EventCutList);
  task->SetConversionCutList(numberOfCuts,ConvCutList);
  task->SetMesonCutList(numberOfCuts,MesonCutList);
  task->SetMoveParticleAccordingToVertex(kTRUE);
  task->SetDoMesonAnalysis(kTRUE);
  task->SetDoMesonQA(enableQAMesonTask); //Attention new switch for Pi0 QA
  task->SetDoPhotonQA(enableQAPhotonTask); //Attention new switch small for Photon QA
  task->SetDoPlotVsCentrality(enablePlotVsCentrality);
  if (trainConfig ==1013 || trainConfig ==1014){
          task->SetDoTHnSparse(0);
  }
  if (enableClustersForTrigger){
    task->SetDoClusterSelectionForTriggerNorm(enableClustersForTrigger);
    task->SetClusterCutList(numberOfCuts,ClusterCutList);
  }

  if (initializedMatBudWeigths_existing) {
      task->SetDoMaterialBudgetWeightingOfGammasForTrueMesons(kTRUE);
  }


  //connect containers
  AliAnalysisDataContainer *coutput =
      mgr->CreateContainer(!(corrTaskSetting.CompareTo("")) ? Form("GammaConvV1_%i",trainConfig) : Form("GammaConvV1_%i_%s",trainConfig,corrTaskSetting.Data()), TList::Class(),
                AliAnalysisManager::kOutputContainer, Form("GammaConvV1_%i.root",trainConfig) );

  mgr->AddTask(task);
  mgr->ConnectInput(task,0,cinput);
  mgr->ConnectOutput(task,1,coutput);
  Int_t nContainer = 2;
  for(Int_t i = 0; i<numberOfCuts; i++){
    if(enableQAPhotonTask>1){
      if (initializedMatBudWeigths_existing) {
	mgr->ConnectOutput(task,nContainer,mgr->CreateContainer(Form("%s_%s_%s MBW Photon DCA tree",(cuts.GetEventCut(i)).Data(),(cuts.GetPhotonCut(i)).Data(),(cuts.GetMesonCut(i)).Data()), TTree::Class(), AliAnalysisManager::kOutputContainer, Form("GammaConvV1_%i.root",trainConfig)) );
      }else{
	mgr->ConnectOutput(task,nContainer,mgr->CreateContainer(Form("%s_%s_%s Photon DCA tree",(cuts.GetEventCut(i)).Data(),(cuts.GetPhotonCut(i)).Data(),(cuts.GetMesonCut(i)).Data()), TTree::Class(), AliAnalysisManager::kOutputContainer, Form("GammaConvV1_%i.root",trainConfig)) );
      }
      nContainer++;
    }
    if(enableQAMesonTask>1){
      if (initializedMatBudWeigths_existing) {
	mgr->ConnectOutput(task,nContainer,mgr->CreateContainer(Form("%s_%s_%s MBW Meson DCA tree",(cuts.GetEventCut(i)).Data(),(cuts.GetPhotonCut(i)).Data(),(cuts.GetMesonCut(i)).Data()), TTree::Class(), AliAnalysisManager::kOutputContainer, Form("GammaConvV1_%i.root",trainConfig)) );
      }else{
	mgr->ConnectOutput(task,nContainer,mgr->CreateContainer(Form("%s_%s_%s Meson DCA tree",(cuts.GetEventCut(i)).Data(),(cuts.GetPhotonCut(i)).Data(),(cuts.GetMesonCut(i)).Data()), TTree::Class(), AliAnalysisManager::kOutputContainer, Form("GammaConvV1_%i.root",trainConfig)) );
      }
      nContainer++;
    }
  }

  return;

}
