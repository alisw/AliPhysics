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
//($ALIPHYSICS/PWGGA/GammaConv/AliAnalysisTaskNeutralMesonToPiPlPiMiNeutralMeson.cxx) for
//pp together with all supporting classes
//***************************************************************************************

//***************************************************************************************
//main function
//***************************************************************************************
void AddTask_GammaConvNeutralMesonPiPlPiMiNeutralMeson_CaloMode_pp(
    Int_t     trainConfig                 = 1,
    Int_t     isMC                        = 0,                        //run MC
    TString   photonCutNumberV0Reader     = "",                       // 00000008400000000100000000 nom. B, 00000088400000000100000000 low B
    Int_t     selectHeavyNeutralMeson     = 0,                        //run eta prime instead of omega
    Int_t     enableQAMesonTask           = 1,                        //enable QA in AliAnalysisTaskNeutralMesonToPiPlPiMiNeutralMeson
    Int_t     enableExtMatchAndQA         = 0,                        // disabled (0), extMatch (1), extQA_noCellQA (2), extMatch+extQA_noCellQA (3), extQA+cellQA (4), extMatch+extQA+cellQA (5)
    Int_t     enableTriggerMimicking      = 0,                        // enable trigger mimicking
    Bool_t    enableTriggerOverlapRej     = kFALSE,                   // enable trigger overlap rejection
    TString   settingMaxFacPtHard         = "3.",                     // maximum factor between hardest jet and ptHard generated
    TString   fileNameExternalInputs      = "",                       // path to file for weigting input
    Bool_t    doWeighting                 = kFALSE,                   //enable Weighting
    TString   generatorName               = "HIJING",
    Double_t  tolerance                   = -1,
    TString   periodNameV0Reader          = "",                       // period Name for V0Reader
    Int_t     runLightOutput              = 0,                        // run light output option 0: no light output 1: most cut histos stiched off 2: unecessary omega hists turned off as well
    Int_t     prefilterRunFlag            = 1500,                     // flag to change the prefiltering of ESD tracks. See SetHybridTrackCutsAODFiltering() in AliPrimaryPionCuts
    Bool_t    usePtDepSelectionWindowCut  = kFALSE,                   // use pt dependent meson selection window cut
    TString   additionalTrainConfig       = "0"                       // additional counter for trainconfig, this has to be always the last parameter
  ) {

  AliCutHandlerPCM cuts(13);

  //parse additionalTrainConfig flag
  Int_t trackMatcherRunningMode = 0; // CaloTrackMatcher running mode

  TObjArray *rAddConfigArr = additionalTrainConfig.Tokenize("_");
  if(rAddConfigArr->GetEntries()<1){std::cout << "ERROR during parsing of additionalTrainConfig String '" << additionalTrainConfig.Data() << "'" << std::endl; return;}
  TObjString* rAdditionalTrainConfig;
  for(Int_t i = 0; i<rAddConfigArr->GetEntries() ; i++){
    if(i==0) rAdditionalTrainConfig = (TObjString*)rAddConfigArr->At(i);
    else{
      TObjString* temp = (TObjString*) rAddConfigArr->At(i);
      TString tempStr = temp->GetString();
      if(tempStr.BeginsWith("TM")){
        TString tempType = tempStr;
        tempType.Replace(0,2,"");
        trackMatcherRunningMode = tempType.Atoi();
        cout << Form("INFO: AddTask_GammaConvNeutralMesonPiPlPiMiPiZero_CaloMode_pp will use running mode '%i' for the TrackMatcher!",trackMatcherRunningMode) << endl;
      }
    }
  }
  TString sAdditionalTrainConfig = rAdditionalTrainConfig->GetString();
  if (sAdditionalTrainConfig.Atoi() > 0){
    trainConfig = trainConfig + sAdditionalTrainConfig.Atoi();
    std::cout << "INFO: AddTask_GammaConvNeutralMesonPiPlPiMiNeutralMeson_CaloMode_pp running additionalTrainConfig '" << sAdditionalTrainConfig.Atoi() << "', train config: '" << trainConfig << "'" << std::endl;
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

  Int_t isHeavyIon = 0;
  Int_t neutralPionMode = 2;

  // Get additional inputs

  TString fileNameCustomTriggerMimicOADB   = cuts.GetSpecialFileNameFromString (fileNameExternalInputs, "FTRM:");

  // ================== GetAnalysisManager ===============================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error(Form("AddTask_GammaConvNeutralMesonPiPlPiMiNeutralMeson_CaloMode_pp_%i",trainConfig), "No analysis manager found.");
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
    std::cout << "V0Reader: " << V0ReaderName.Data() << " not found!!"<< std::endl;
    return;
  } else {
    std::cout << "V0Reader: " << V0ReaderName.Data() << " found!!"<< std::endl;
  }

  TString PionCuts      = "000000200";
  if(periodNameV0Reader.Contains("LHC11") ||
     periodNameV0Reader.Contains("LHC14k1") ||
     periodNameV0Reader.Contains("LHC14b7") ||
     periodNameV0Reader.Contains("LHC16c2")){
       PionCuts      = "000000000";
  }
  //================================================
  //========= Add Pion Selector ====================
  TString PionSelectorName  =  Form("PionSelector_%s", PionCuts.Data());
  if( !(AliPrimaryPionSelector*)mgr->GetTask(PionSelectorName.Data()) ){
    AliPrimaryPionSelector *fPionSelector = new AliPrimaryPionSelector(PionSelectorName.Data());
    AliPrimaryPionCuts *fPionCuts=0;
    if( PionCuts!=""){
      fPionCuts= new AliPrimaryPionCuts(PionCuts.Data(),PionCuts.Data());
      fPionCuts->SetPrefilterRunFlag(prefilterRunFlag);
      if(runLightOutput>0) fPionCuts->SetLightOutput(kTRUE);
      fPionCuts->SetPeriodName(periodNameV0Reader);
      if(fPionCuts->InitializeCutsFromCutString(PionCuts.Data())){
        fPionSelector->SetPrimaryPionCuts(fPionCuts);
        fPionCuts->SetFillCutHistograms("",kTRUE);
      }
    }

    fPionSelector->Init();
    mgr->AddTask(fPionSelector);

    AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer();
    mgr->ConnectInput (fPionSelector,0,cinput1);
  }

  AliAnalysisTaskNeutralMesonToPiPlPiMiNeutralMeson *task=NULL;
  task= new AliAnalysisTaskNeutralMesonToPiPlPiMiNeutralMeson(Form("GammaConvNeutralMesonPiPlPiMiNeutralMeson_%i_%i",neutralPionMode, trainConfig));
  task->SetPionSelectorName(PionSelectorName.Data());
  task->SetIsHeavyIon(isHeavyIon);
  task->SetIsMC(isMC);
  task->SetV0ReaderName(V0ReaderName);
  if(runLightOutput>1) task->SetLightOutput(kTRUE);
  task->SetTolerance(tolerance);
  task->SetTrackMatcherRunningMode(trackMatcherRunningMode);




  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //                                          ETA MESON
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  if( trainConfig == 1 ) {
    // everything open, min pt charged pi = 100 MeV
    cuts.AddCutHeavyMesonCalo("00000113","1111113047032230000","000010400","0103503a00000000","0103503000000000");
    cuts.AddCutHeavyMesonCalo("00000113","1111113047032230000","000010400","0103503a00000000","0103503000000000");
  } else if( trainConfig == 2 ) {
    // closing charged pion cuts, minimum TPC cluster = 80, TPC dEdx pi = \pm 3 sigma, min pt charged pi = 100 MeV
    cuts.AddCutHeavyMesonCalo("00000113","1111100047032230000","002010700","0103503a00000000","0103503000000000");
  } else if( trainConfig == 3) {
    // eta < 0.9
    // closing charged pion cuts, minimum TPC cluster = 80, TPC dEdx pi = \pm 3 sigma, pi+pi- mass cut of 0.65, min pt charged pi = 100 MeV
    // closing neural pion cuts, 0.1 < M_gamma,gamma < 0.15
    // maxChi2 per cluster TPC <4, require TPC refit, DCA XY pT dependend 0.0182+0.0350/pt^1.01, DCA_Z = 3.0
    // timing cluster cut open
    cuts.AddCutHeavyMesonCalo("00010113","1111100047032230000","30a330708","0103503400000000","0153503000000000"); // all of the above
    cuts.AddCutHeavyMesonCalo("00052113","1111100047032230000","30a330708","0103503400000000","0153503000000000"); // all of the above
    cuts.AddCutHeavyMesonCalo("00062113","1111100047032230000","30a330708","0103503400000000","0153503000000000"); // all of the above
    cuts.AddCutHeavyMesonCalo("00083113","1111100047032230000","30a330708","0103503400000000","0153503000000000"); // all of the above
    cuts.AddCutHeavyMesonCalo("00085113","1111100047032230000","30a330708","0103503400000000","0153503000000000"); // all of the above
  } else if( trainConfig == 4) {
    // same as 3 but only MB
    cuts.AddCutHeavyMesonCalo("00010113","1111100047032230000","30a330708","0103503400000000","0153503000000000"); // all of the above
  } else if( trainConfig == 10)  { // Standard EMCal (7 TeV)
    cuts.AddCutHeavyMesonCalo("00000113","1111111047032230000","302010708","0103603900000000","0153503000000000");
    cuts.AddCutHeavyMesonCalo("00000113","1111111047032230000","322010708","0103603900000000","0153503000000000"); // with ITS requirement
  } else if( trainConfig == 11)  { // Test for EMCal (5TeV)
    cuts.AddCutHeavyMesonCalo("00010113","1111111047032230000","302010708","0103603900000000","0153503000000000");
    cuts.AddCutHeavyMesonCalo("00010113","1111111047032230000","322010708","0103603900000000","0153503000000000"); // with ITS requirement
  } else if( trainConfig == 50)  { // Standard PHOS
    cuts.AddCutHeavyMesonCalo("00000113","2444411043012300000","302010708","0103603200000000","0153503000000000"); // PCM-PHOS nonLin
    cuts.AddCutHeavyMesonCalo("00000113","2444411043012300000","322010708","0103603200000000","0153503000000000"); // with ITS requirement
  } else if( trainConfig == 51)  { // Test for PHOS (5 TeV)
    cuts.AddCutHeavyMesonCalo("00010113","2444411043012300000","302010708","0103603200000000","0153503000000000"); // PCM-PHOS nonLin
    cuts.AddCutHeavyMesonCalo("00010113","2444411043012300000","322010708","0103603200000000","0153503000000000"); // with ITS requirement

  } else if( trainConfig == 60) { //pp 13 TeV EDC
    cuts.AddCutHeavyMesonCalo("00010113","411791206f032230000","30a210708","01631031000000d0","0153503000000000"); // INT7 NL 12
    cuts.AddCutHeavyMesonCalo("00010113","411791206f032230000","32a210708","01631031000000d0","0153503000000000"); // INT7 NL 12
  } else if( trainConfig == 61) { //pp 13 TeV EDC Trigger
    cuts.AddCutHeavyMesonCalo("0008e113","411791206f032230000","30a210708","01631031000000d0","0153503000000000"); // EG2 NL 12
    cuts.AddCutHeavyMesonCalo("0008d113","411791206f032230000","30a210708","01631031000000d0","0153503000000000"); // EG1 NL 12
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //                                          OMEGA MESON
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  } else if( trainConfig == 100 ) {
    // everything open, min pt charged pi = 100 MeV
    cuts.AddCutHeavyMesonCalo("00000113","1111100047032230000","000010400","0103503a00000000","0103503000000000");
  } else if( trainConfig == 101 ) {
    // closing charged pion cuts, minimum TPC cluster = 80, TPC dEdx pi = \pm 3 sigma, min pt charged pi = 100 MeV
    cuts.AddCutHeavyMesonCalo("00000113","1111100047032230000","002010700","0103503a00000000","0103503000000000");
  } else if( trainConfig == 102) {
    // eta < 0.9
    // closing charged pion cuts, minimum TPC cluster = 80, TPC dEdx pi = \pm 3 sigma, pi+pi- mass cut of 0.65, min pt charged pi = 100 MeV
    // closing neural pion cuts, 0.1 < M_gamma,gamma < 0.15
    // maxChi2 per cluster TPC <4, require TPC refit, DCA XY pT dependend 0.0182+0.0350/pt^1.01, DCA_Z = 3.0
    // timing cluster cut open
    cuts.AddCutHeavyMesonCalo("00010113","1111100047032230000","30a330708","0103503400000000","0153503000000000"); // all of the above
    cuts.AddCutHeavyMesonCalo("00052113","1111100047032230000","30a330708","0103503400000000","0153503000000000"); // all of the above
    cuts.AddCutHeavyMesonCalo("00062113","1111100047032230000","30a330708","0103503400000000","0153503000000000"); // all of the above
    cuts.AddCutHeavyMesonCalo("00083113","1111100047032230000","30a330708","0103503400000000","0153503000000000"); // all of the above
    cuts.AddCutHeavyMesonCalo("00085113","1111100047032230000","30a330708","0103503400000000","0153503000000000"); // all of the above

  } else if( trainConfig == 103) {
    // same as 102 but only MB
    cuts.AddCutHeavyMesonCalo("00010113","1111100047032230000","30a330708","0103503400000000","0153503000000000"); // all of the above

  } else if( trainConfig == 110)  { // Standard EMCal
    cuts.AddCutHeavyMesonCalo("00000113","1111111047032230000","32c010708","0103603700000000","0153503000000000");

    // EMCal pp 5 TeV
  } else if( trainConfig == 111)  { // Test for EMCal (5 TeV) with TPC refit + ITS requirement
    cuts.AddCutHeavyMesonCalo("00010113","1111111047032230000","32c010708","0103603700000000","0153503000000000"); // INT7
    cuts.AddCutHeavyMesonCalo("00052113","1111111047032230000","32c010708","0103603700000000","0153503000000000"); // EMC7
    cuts.AddCutHeavyMesonCalo("00085113","1111111047032230000","32c010708","0103603700000000","0153503000000000"); // EG2
    cuts.AddCutHeavyMesonCalo("00083113","1111111047032230000","32c010708","0103603700000000","0153503000000000"); // EG1
    // EMCal pp 13 TeV
  } else if( trainConfig == 112)  { // Test for EMCal (13 TeV) with TPC refit + ITS requirement
    cuts.AddCutHeavyMesonCalo("00010113","1111111047032230000","32c51070a","0103603700000000","0153503000000000"); // INT7
    cuts.AddCutHeavyMesonCalo("00085113","1111111047032230000","32c51070a","0103603700000000","0153503000000000"); // EG2
    cuts.AddCutHeavyMesonCalo("00083113","1111111047032230000","32c51070a","0103603700000000","0153503000000000"); // EG1
  } else if( trainConfig == 113)  { // AOD and ESD Comparison
    cuts.AddCutHeavyMesonCalo("00010113","1111111047032230000","32c510008","0103603700000000","0153503000000000"); // INT7

    // EMCal LHC11 pp 7TeV
  } else if( trainConfig == 115){ // EMCal LHC11 std with background use with pt dep mass cut
    cuts.AddCutHeavyMesonCalo("00010113","111113106f032230000","32c51070a","0103603600000000","0153503000000000"); // INT7 (LHC11 acc) 3 sigma
    cuts.AddCutHeavyMesonCalo("00052113","111113106f032230000","32c51070a","0103603600000000","0153503000000000"); // EMC7 (LHC11 acc)
  } else if( trainConfig == 116){ // with new SPD pileup cut
    cuts.AddCutHeavyMesonCalo("00010c13","111113106f032230000","32c51070a","0103603600000000","0153503000000000"); // INT7 (LHC11 acc) 3 sigma
    cuts.AddCutHeavyMesonCalo("00052c13","111113106f032230000","32c51070a","0103683600000000","0153503000000000"); // EMC7 min Pt cut 5 GeV
    cuts.AddCutHeavyMesonCalo("00052c13","111113106f032230000","32c51070a","01036d3600000000","0153503000000000"); // 10.0
  // inv mass cut variation
  } else if( trainConfig == 117){ // pt cut variation old SPD pileup cut
    cuts.AddCutHeavyMesonCalo("00052113","111113106f032230000","32c51070a","01036c3600000000","0153503000000000"); // 8.0
    cuts.AddCutHeavyMesonCalo("00052113","111113106f032230000","32c51070a","01036e3600000000","0153503000000000"); // new Energy cut 5 GeV
    cuts.AddCutHeavyMesonCalo("00052113","111113106f032230000","32c51070a","01036c3600000000","0153503000000000"); // new Energy cut 7.5 GeV
    cuts.AddCutHeavyMesonCalo("00052113","111113106f032230000","32c51070a","01036g3600000000","0153503000000000"); // new Energy cut 6.0 GeV
  } else if( trainConfig == 118){ // pt cut variation new SPD pileup cut
    cuts.AddCutHeavyMesonCalo("00052c13","111113106f032230000","32c51070a","01036c3600000000","0153503000000000"); // 8.0
    cuts.AddCutHeavyMesonCalo("00052c13","111113106f032230000","32c51070a","01036e3600000000","0153503000000000"); // new Energy cut 5 GeV
    cuts.AddCutHeavyMesonCalo("00052c13","111113106f032230000","32c51070a","01036f3600000000","0153503000000000"); // new Energy cut 7.5 GeV
  } else if( trainConfig == 119){ // no PID cut
    cuts.AddCutHeavyMesonCalo("00010113","111113106f032230000","32c51070a","0103603600000000","0153503000000000"); // INT7 (LHC11 acc) 3 sigma
    cuts.AddCutHeavyMesonCalo("00052113","111113106f032230000","32c51070a","01036c3600000000","0153503000000000"); // 8.0
    cuts.AddCutHeavyMesonCalo("00052113","111113106f032230000","32c51070a","01036e3600000000","0153503000000000"); // new Energy cut 5 GeV
    cuts.AddCutHeavyMesonCalo("00052113","111113106f032230000","32c51070a","01036c3600000000","0153503000000000"); // new Energy cut 7.5 GeV
    cuts.AddCutHeavyMesonCalo("00052113","111113106f032230000","32c51070a","01036g3600000000","0153503000000000"); // new Energy cut 6.0 GeV
    // ---------------------------------
    // systematic studies 7 TeV (EMCal)
    // ---------------------------------

    // charged pion cuts
  } else if( trainConfig == 120)   {
    cuts.AddCutHeavyMesonCalo("00000113","1111a3104f032230000","32c51070a","0103603600000000","0153503000000000"); // with TPC refit + ITS requirement
    cuts.AddCutHeavyMesonCalo("00000113","1111a3104f032230000","32c51070a","0103603c00000000","0153503000000000"); // with TPC refit + ITS requirement
  } else if (trainConfig == 121){ // remove pileup
    cuts.AddCutHeavyMesonCalo("00000013","1111a3104f032230000","32c51070a","0103603c00000000","0153503000000000"); // rmeove pileup
  } else if (trainConfig == 122){ // nonlin
    cuts.AddCutHeavyMesonCalo("00000113","1111a12047032230000","32c51070a","0103603c00000000","0153503000000000"); // NonLinearity ConvCalo - kTestBeamv3 + shifting MC
    cuts.AddCutHeavyMesonCalo("00000113","1111a13047032230000","32c51070a","0103603c00000000","0153503000000000"); // NonLinearity pp Calo - only shifting MC - no timing cut
    cuts.AddCutHeavyMesonCalo("00000113","1111a21047032230000","32c51070a","0103603c00000000","0153503000000000"); // NonLinearity pp ConvCalo - only shifting MC - no timing cut (Fits)
    cuts.AddCutHeavyMesonCalo("00000113","1111a22047032230000","32c51070a","0103603c00000000","0153503000000000"); // NonLinearity pp ConvCalo - only shifting MC - no timing cut (Fits)
    cuts.AddCutHeavyMesonCalo("00000113","1111a23047032230000","32c51070a","0103603c00000000","0153503000000000");  // NonLinearity ConvCalo - kTestBeamv3 + shifting MC
  } else if (trainConfig == 123){ // timing diff
    cuts.AddCutHeavyMesonCalo("00000113","1111a11037032230000","32c51070a","0103603c00000000","0153503000000000"); // timing diff
    cuts.AddCutHeavyMesonCalo("00000113","1111a11057032230000","32c51070a","0103603c00000000","0153503000000000"); // timing diff
    cuts.AddCutHeavyMesonCalo("00000113","1111a11077032230000","32c51070a","0103603c00000000","0153503000000000"); // timing diff
    cuts.AddCutHeavyMesonCalo("00000113","1111a11097032230000","32c51070a","0103603c00000000","0153503000000000"); // timing diff
  } else if (trainConfig == 124){ // TrackMatching
    cuts.AddCutHeavyMesonCalo("00000113","1111a11046032230000","32c51070a","0103603c00000000","0153503000000000");
    cuts.AddCutHeavyMesonCalo("00000113","1111a11048032230000","32c51070a","0103603c00000000","0153503000000000");
    cuts.AddCutHeavyMesonCalo("00000113","1111a11049032230000","32c51070a","0103603c00000000","0153503000000000");
    cuts.AddCutHeavyMesonCalo("00000113","1111a1104a032230000","32c51070a","0103603c00000000","0153503000000000");
    cuts.AddCutHeavyMesonCalo("00000113","1111a1104b032230000","32c51070a","0103603c00000000","0153503000000000");
    cuts.AddCutHeavyMesonCalo("00000113","1111a11043032230000","32c51070a","0103603c00000000","0153503000000000");
    cuts.AddCutHeavyMesonCalo("00000113","1111a1104f032230000","32c51070a","0103603c00000000","0153503000000000"); // E/p
  } else if (trainConfig == 125){ // MinEnergy (of cluster)
    cuts.AddCutHeavyMesonCalo("00000113","1111a11047022230000","32c51070a","0103603c00000000","0153503000000000"); // 0.6
    cuts.AddCutHeavyMesonCalo("00000113","1111a11047042230000","32c51070a","0103603c00000000","0153503000000000"); // 0.8
    cuts.AddCutHeavyMesonCalo("00000113","1111a11047052230000","32c51070a","0103603c00000000","0153503000000000"); // 0.9
  } else if (trainConfig == 126){ // MinNCells
    cuts.AddCutHeavyMesonCalo("00000113","1111a11047031230000","32c51070a","0103603c00000000","0153503000000000"); // 1
    cuts.AddCutHeavyMesonCalo("00000113","1111a11047033230000","32c51070a","0103603c00000000","0153503000000000"); // 3
  } else if (trainConfig == 127){ // MinMaxM02
    cuts.AddCutHeavyMesonCalo("00000113","1111a11047032330000","32c51070a","0103603c00000000","0153503000000000"); // 0.2 - 0.5
    cuts.AddCutHeavyMesonCalo("00000113","1111a11047032130000","32c51070a","0103603c00000000","0153503000000000"); // 0.002 - 0.5
    cuts.AddCutHeavyMesonCalo("00000113","1111a11047032240000","32c51070a","0103603c00000000","0153503000000000"); // 0.1 - 0.4
    cuts.AddCutHeavyMesonCalo("00000113","1111a11047032220000","32c51070a","0103603c00000000","0153503000000000"); // 0.1 - 0.7

  // other variations
  } else if (trainConfig == 140) { // pT Cut
    cuts.AddCutHeavyMesonCalo("00000113","1111a11047032230000","32c50070a","0103603700000000","0153503000000000"); // pt>0.075
    cuts.AddCutHeavyMesonCalo("00000113","1111a11047032230000","32c52070a","0103603700000000","0153503000000000"); // pt>0.125
    cuts.AddCutHeavyMesonCalo("00000113","1111a11047032230000","32c53070a","0103603700000000","0153503000000000"); // pt>0.15
    cuts.AddCutHeavyMesonCalo("00000113","1111a11047032230000","32c54070a","0103603700000000","0153503000000000"); // pt>0.4
  } else if (trainConfig == 141) { // TPCdEdxCutPion
    cuts.AddCutHeavyMesonCalo("00000113","1111a11047032230000","32c51050a","0103603700000000","0153503000000000"); // -4,4
    cuts.AddCutHeavyMesonCalo("00000113","1111a11047032230000","32c51080a","0103603700000000","0153503000000000"); // -2,3
    cuts.AddCutHeavyMesonCalo("00000113","1111a11047032230000","32c51020a","0103603700000000","0153503000000000"); // -6,7
    cuts.AddCutHeavyMesonCalo("00000113","1111a11047032230000","32c51030a","0103603700000000","0153503000000000"); // -5,5
    cuts.AddCutHeavyMesonCalo("00000113","1111a11047032230000","32c51040a","0103603700000000","0153503000000000"); // -4,5
  // neutral pion cuts
  } else if (trainConfig == 142) { // Mass window
    cuts.AddCutHeavyMesonCalo("00000113","1111a11047032230000","32c51070a","0103603100000000","0153503000000000"); // 100-145
    cuts.AddCutHeavyMesonCalo("00000113","1111a11047032230000","32c51070a","0103603200000000","0153503000000000"); // 110-145
    cuts.AddCutHeavyMesonCalo("00000113","1111a11047032230000","32c51070a","0103603300000000","0153503000000000"); // 120-145
    cuts.AddCutHeavyMesonCalo("00000113","1111a11047032230000","32c51070a","0103603400000000","0153503000000000"); // 100-150
    cuts.AddCutHeavyMesonCalo("00000113","1111a11047032230000","32c51070a","0103603500000000","0153503000000000"); // 110-150
    cuts.AddCutHeavyMesonCalo("00000113","1111a11047032230000","32c51070a","0103603600000000","0153503000000000"); // 120-150
    cuts.AddCutHeavyMesonCalo("00000113","1111a11047032230000","32c51070a","0103603a00000000","0153503000000000"); // 80-145

  } else if( trainConfig == 145)   { // no background calculation
    cuts.AddCutHeavyMesonCalo("00000113","1111a11047032230000","32c51070a","0103603700000000","0453503000000000");

    // ---------------------------------
    // systematic studies 7 TeV (PHOS)
    // ---------------------------------

  } else if(trainConfig == 150)  { // Standard PHOS
    cuts.AddCutHeavyMesonCalo("00000113","2444411044013300000","32c51070a","0103603n00000000","0153503000000000"); //  with TPC refit + ITS requirement

    // *************Variations in AliConvEventCuts**************************
  } else if(trainConfig == 151)  { // removePileUp
    cuts.AddCutHeavyMesonCalo("00000013","2444411044013300000","32c01070a","0103603n00000000","0153503000000000"); // PC

    // *************Variations in AliCaloPhotonsCut**************************
  } else if(trainConfig == 152)  { // Timing diff(std is -100ns to 100ns)
    cuts.AddCutHeavyMesonCalo("00000113","24444110b4013300000","32c51070a","0103603n00000000","0153503000000000"); // 130ns
    cuts.AddCutHeavyMesonCalo("00000113","24444110c4013300000","32c51070a","0103603n00000000","0153503000000000"); // 110ns
    cuts.AddCutHeavyMesonCalo("00000113","24444110d4013300000","32c51070a","0103603n00000000","0153503000000000"); // 120ns
    cuts.AddCutHeavyMesonCalo("00000113","24444110e4013300000","32c51070a","0103603n00000000","0153503000000000"); // 90ns
    cuts.AddCutHeavyMesonCalo("00000113","24444110f4013300000","32c51070a","0103603n00000000","0153503000000000"); // 80ns
  } else if(trainConfig == 153)  { // Track matching
    cuts.AddCutHeavyMesonCalo("00000113","2444411041013300000","32c51070a","0103603n00000000","0153503000000000"); //
    cuts.AddCutHeavyMesonCalo("00000113","2444411043013300000","32c51070a","0103603n00000000","0153503000000000"); //
    cuts.AddCutHeavyMesonCalo("00000113","2444411045013300000","32c51070a","0103603n00000000","0153503000000000"); //
    cuts.AddCutHeavyMesonCalo("00000113","2444411046013300000","32c51070a","0103603n00000000","0153503000000000"); //
    cuts.AddCutHeavyMesonCalo("00000113","2444411047013300000","32c51070a","0103603n00000000","0153503000000000"); //
    cuts.AddCutHeavyMesonCalo("00000113","2444411048013300000","32c51070a","0103603n00000000","0153503000000000"); //
    cuts.AddCutHeavyMesonCalo("00000113","2444411049013300000","32c51070a","0103603n00000000","0153503000000000"); //
  } else if(trainConfig == 154)  { // MinEnergy (of cluster) (std is 0.5 GeV)
    cuts.AddCutHeavyMesonCalo("00000113","2444411044023300000","32c51070a","0103603n00000000","0153503000000000"); // 0.3
    cuts.AddCutHeavyMesonCalo("00000113","2444411044033300000","32c51070a","0103603n00000000","0153503000000000"); // 0.7
    cuts.AddCutHeavyMesonCalo("00000113","2444411044043300000","32c51070a","0103603n00000000","0153503000000000"); // 0.8
    cuts.AddCutHeavyMesonCalo("00000113","2444411044083300000","32c51070a","0103603n00000000","0153503000000000"); // 0.4
  } else if(trainConfig == 155)  { // Min N of cells (std is 2)
    cuts.AddCutHeavyMesonCalo("00000113","2444411044011300000","32c51070a","0103603n00000000","0153503000000000"); // 1
    cuts.AddCutHeavyMesonCalo("00000113","2444411044013300000","32c51070a","0103603n00000000","0153503000000000"); // 3
    cuts.AddCutHeavyMesonCalo("00000113","2444411044014300000","32c51070a","0103603n00000000","0153503000000000"); // 4
  } else if(trainConfig == 156)  { // MinMaxM02 (std is >0.2)
    cuts.AddCutHeavyMesonCalo("00000113","2444411044013200000","32c51070a","0103603n00000000","0153503000000000"); // >0.1
    cuts.AddCutHeavyMesonCalo("00000113","2444411044013100000","32c51070a","0103603n00000000","0153503000000000"); // >0.002

    // *************Variations in AliPrimaryPionCuts******************
  } else if( trainConfig == 157)  { // ClsTPCCut
    cuts.AddCutHeavyMesonCalo("00000113","2444411044013300000","32151070a","0103603n00000000","0153503000000000"); // 70
    cuts.AddCutHeavyMesonCalo("00000113","2444411044013300000","32351070a","0103603n00000000","0153503000000000"); // 100
    cuts.AddCutHeavyMesonCalo("00000113","2444411044013300000","32551070a","0103603n00000000","0153503000000000"); // 35% of findable clusters
    cuts.AddCutHeavyMesonCalo("00000113","2444411044013300000","32651070a","0103603n00000000","0153503000000000"); // 60% of findable clusters
    cuts.AddCutHeavyMesonCalo("00000113","2444411044013300000","32b51070a","0103603n00000000","0153503000000000"); // PHOS public note

  } else if ( trainConfig == 158) { // DCACut
    cuts.AddCutHeavyMesonCalo("00000113","2444411044013300000","32c11070a","0103603n00000000","0153503000000000"); // XYPtDep("0.0182+0.0350/pt^1.01");
    cuts.AddCutHeavyMesonCalo("00000113","2444411044013300000","32c21070a","0103603n00000000","0153503000000000"); // z=2cm xy=1cm
    cuts.AddCutHeavyMesonCalo("00000113","2444411044013300000","32c31070a","0103603n00000000","0153503000000000"); // z=3cm XYPtDep("0.0182+0.0350/pt^1.01");
    cuts.AddCutHeavyMesonCalo("00000113","2444411044013300000","32c41070a","0103603n00000000","0153503000000000"); // z=3cm xy=0.5
  } else if ( trainConfig == 159) { // pT cut
    cuts.AddCutHeavyMesonCalo("00000113","2444411044013300000","32c50070a","0103603n00000000","0153503000000000"); // pt>0.075
    cuts.AddCutHeavyMesonCalo("00000113","2444411044013300000","32c52070a","0103603n00000000","0153503000000000"); // pt>0.125
    cuts.AddCutHeavyMesonCalo("00000113","2444411044013300000","32c53070a","0103603n00000000","0153503000000000"); // pt>0.15
    cuts.AddCutHeavyMesonCalo("00000113","2444411044013300000","32c54070a","0103603n00000000","0153503000000000"); // pt>0.4
  } else if ( trainConfig == 160) { // TPDdEdxCutPion
    cuts.AddCutHeavyMesonCalo("00000113","2444411044013300000","32c51050a","0103603n00000000","0153503000000000"); // -4,4
    cuts.AddCutHeavyMesonCalo("00000113","2444411044013300000","32c51080a","0103603n00000000","0153503000000000"); // -2,3
    cuts.AddCutHeavyMesonCalo("00000113","2444411044013300000","32c51020a","0103603n00000000","0153503000000000"); // -6,7
    cuts.AddCutHeavyMesonCalo("00000113","2444411044013300000","32c51030a","0103603n00000000","0153503000000000"); // -5,5
    cuts.AddCutHeavyMesonCalo("00000113","2444411044013300000","32c51040a","0103603n00000000","0153503000000000"); // -4,5
  } else if ( trainConfig == 161) { // Mass cut
    cuts.AddCutHeavyMesonCalo("00000113","2444411044013300000","32c510707","0103603n00000000","0153503000000000"); // 0.7 GeV
    cuts.AddCutHeavyMesonCalo("00000113","2444411044013300000","32c510706","0103603n00000000","0153503000000000"); // 0.65 GeV
    cuts.AddCutHeavyMesonCalo("00000113","2444411044013300000","32c510701","0103603n00000000","0153503000000000"); // 1 GeV
    cuts.AddCutHeavyMesonCalo("00000113","2444411044013300000","32c510702","0103603n00000000","0153503000000000"); // 0.75 GeV
    cuts.AddCutHeavyMesonCalo("00000113","2444411044013300000","32c510704","0103603n00000000","0153503000000000"); // 0.54 eta mass
    // *************Variations in AliConversionMesonCuts (NeutralPion) ******************
  } else if ( trainConfig == 162) { // RapidityMesonCut
    cuts.AddCutHeavyMesonCalo("00000113","2444411044013300000","32c51070a","0103503n00000000","0153503000000000"); // 0.85
    cuts.AddCutHeavyMesonCalo("00000113","2444411044013300000","32c51070a","0103303n00000000","0153503000000000"); // 0.6
    cuts.AddCutHeavyMesonCalo("00000113","2444411044013300000","32c51070a","0103203n00000000","0153503000000000"); // 0.7
  } else if ( trainConfig == 163) { // pT
    cuts.AddCutHeavyMesonCalo("00000113","2444411044013300000","32c51070a","0103613n00000000","0153503000000000"); // 0.4
    cuts.AddCutHeavyMesonCalo("00000113","2444411044013300000","32c51070a","0103623n00000000","0153503000000000"); // 0.7
    cuts.AddCutHeavyMesonCalo("00000113","2444411044013300000","32c51070a","0103673n00000000","0153503000000000"); // 0.5
  } else if ( trainConfig == 164) { // max opening angle cut
    cuts.AddCutHeavyMesonCalo("00000113","2444411044013300000","32c51070a","0103603n00000001","0153503000000000"); //
    cuts.AddCutHeavyMesonCalo("00000113","2444411044013300000","32c51070a","0103603n00000002","0153503000000000"); //
    cuts.AddCutHeavyMesonCalo("00000113","2444411044013300000","32c51070a","0103603n00000003","0153503000000000"); //
  } else if ( trainConfig == 165) { // selectionWindow (std is 120-160)
    cuts.AddCutHeavyMesonCalo("00000113","2444411044013300000","32c51070a","0103603100000000","0153503000000000"); // 0.1-0.145
    cuts.AddCutHeavyMesonCalo("00000113","2444411044013300000","32c51070a","0103603200000000","0153503000000000"); // 0.11-0.145
    cuts.AddCutHeavyMesonCalo("00000113","2444411044013300000","32c51070a","0103603300000000","0153503000000000"); // 0.12-0.145
    cuts.AddCutHeavyMesonCalo("00000113","2444411044013300000","32c51070a","0103603600000000","0153503000000000"); // 0.12-0.5
    cuts.AddCutHeavyMesonCalo("00000113","2444411044013300000","32c51070a","0103603700000000","0153503000000000"); // 0.1 -0.155
    cuts.AddCutHeavyMesonCalo("00000113","2444411044013300000","32c51070a","0103603900000000","0153503000000000"); // 0.11 -0.155
    cuts.AddCutHeavyMesonCalo("00000113","2444411044013300000","32c51070a","0103603a00000000","0153503000000000"); // 0.08 -0.145
    // *************Variations in AliConversionMesonCuts (omega) ******************
  } else if ( trainConfig == 166) { // selectionWindow (std is 120-160)
    cuts.AddCutHeavyMesonCalo("00000113","2444411044013300000","32c51070a","0103603n00000000","0a53503000000000"); // likesign
    cuts.AddCutHeavyMesonCalo("00000113","2444411044013300000","32c51070a","0103603n00000000","0b53503000000000"); // sideband right
    cuts.AddCutHeavyMesonCalo("00000113","2444411044013300000","32c51070a","0103603n00000000","0c53503000000000"); // sideband left
    cuts.AddCutHeavyMesonCalo("00000113","2444411044013300000","32c51070a","0103603n00000000","0d53503000000000"); // sideband both sides
  } else if ( trainConfig == 167) { // Number of BckEvents
    cuts.AddCutHeavyMesonCalo("00000113","2444411044013300000","32c51070a","0103603n00000000","0133503000000000"); // 20
    cuts.AddCutHeavyMesonCalo("00000113","2444411044013300000","32c51070a","0103603n00000000","0163503000000000"); // 80
  } else if ( trainConfig == 168) { // rapidity cut
    cuts.AddCutHeavyMesonCalo("00000113","2444411044013300000","32c51070a","0103603n00000000","0153203000000000"); // 0.7
    cuts.AddCutHeavyMesonCalo("00000113","2444411044013300000","32c51070a","0103603n00000000","0153003000000000"); // 1.35
    cuts.AddCutHeavyMesonCalo("00000113","2444411044013300000","32c51070a","0103603n00000000","0153303000000000"); // 0.6
  } else if ( trainConfig == 169) { // alpha max cut
    cuts.AddCutHeavyMesonCalo("00000113","2444411044013300000","32c51070a","0103603n00000000","0153503000000000"); // 0-0.85
    cuts.AddCutHeavyMesonCalo("00000113","2444411044013300000","32c51070a","0103603n00000000","0153503000000000"); // 0-0.75
    cuts.AddCutHeavyMesonCalo("00000113","2444411044013300000","32c51070a","0103603n00000000","0153503000000000"); // 0-0.7
  } else if(trainConfig == 170)  { // no back calculation
    cuts.AddCutHeavyMesonCalo("00000113","2444411044013300000","32c51070a","0103603n00000000","0453503000000000"); // no bck calculation

  } else if(trainConfig == 180)  { // Calo non lin
    cuts.AddCutHeavyMesonCalo("00000113","2444412044013300000","32c51070a","0103603n00000000","0153503000000000");
    cuts.AddCutHeavyMesonCalo("00000113","2444421044013300000","32c51070a","0103603n00000000","0153503000000000");

    // pp 7 TeV LHC11 sys
 } else if( trainConfig == 190){ // pileup variation
    cuts.AddCutHeavyMesonCalo("00052113","111113106f032230000","32c51070a","01036c3c00000000","0153503000000000"); // standard
    cuts.AddCutHeavyMesonCalo("00052c13","111113106f032230000","32c51070a","01036c3c00000000","0153503000000000"); // loose
 } else if( trainConfig == 191){ // mass cut variation
    cuts.AddCutHeavyMesonCalo("00052113","111113106f032230000","32c51070a","01036c3600000000","0153503000000000"); // 2 sigma
    cuts.AddCutHeavyMesonCalo("00052113","111113106f032230000","32c51070a","01036c3b00000000","0153503000000000"); // 1 sigma
    cuts.AddCutHeavyMesonCalo("00052113","111113106f032230000","32c51070a","01036c3d00000000","0153503000000000"); // 4 sigma
 } else if( trainConfig == 192){ // pT cut variation
    cuts.AddCutHeavyMesonCalo("00052113","111113106f032230000","32c51070a","01036d3c00000000","0153503000000000"); // 10 GeV
    cuts.AddCutHeavyMesonCalo("00052113","111113106f032230000","32c51070a","01036b3c00000000","0153503000000000"); // 6 GeV
    cuts.AddCutHeavyMesonCalo("00052113","111113106f032230000","32c51070a","01036a3c00000000","0153503000000000"); // 4 GeV
 } else if( trainConfig == 193){ // cluster timing
    cuts.AddCutHeavyMesonCalo("00052113","111113107f032230000","32c51070a","01036c3c00000000","0153503000000000"); // -30 / 30ns
    cuts.AddCutHeavyMesonCalo("00052113","111113108f032230000","32c51070a","01036c3c00000000","0153503000000000"); // -20 / 30ns
    cuts.AddCutHeavyMesonCalo("00052113","111113105f032230000","32c51070a","01036c3c00000000","0153503000000000"); // -50 / 50ns
 } else if( trainConfig == 194){ // PID
    cuts.AddCutHeavyMesonCalo("00052113","111113106f032230000","32c51050a","01036c3c00000000","0153503000000000"); // standard
    cuts.AddCutHeavyMesonCalo("00052113","111113106f032230000","32c51080a","01036c3c00000000","0153503000000000"); // standard
    cuts.AddCutHeavyMesonCalo("00052113","111113106f032230000","32c51020a","01036c3c00000000","0153503000000000"); // standard
    cuts.AddCutHeavyMesonCalo("00052113","111113106f032230000","32c51030a","01036c3c00000000","0153503000000000"); // standard
 } else if( trainConfig == 195){ // background description
    cuts.AddCutHeavyMesonCalo("00052113","111113106f032230000","32c51070a","01036c3c00000000","0a53503000000000"); // standard
    cuts.AddCutHeavyMesonCalo("00052113","111113106f032230000","32c51070a","01036c3c00000000","0b53503000000000"); // standard
    cuts.AddCutHeavyMesonCalo("00052113","111113106f032230000","32c51070a","01036c3c00000000","0c53503000000000"); // standard
    cuts.AddCutHeavyMesonCalo("00052113","111113106f032230000","32c51070a","01036c3c00000000","0d53503000000000"); // standard
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //                                          ETA PRIME MESON
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  } else if( trainConfig == 200 ) {
    // everything open, min pt charged pi = 100 MeV
    cuts.AddCutHeavyMesonCalo("00000113","1111100047032230000","000010400","0103503m00000000","0103503000000000");
  } else if( trainConfig == 201 ) {
    // closing charged pion cuts, minimum TPC cluster = 80, TPC dEdx pi = \pm 3 sigma, min pt charged pi = 100 MeV
    cuts.AddCutHeavyMesonCalo("00000113","1111100047032230000","002010700","0103503m00000000","0103503000000000");
  } else if( trainConfig == 202) {
    // eta < 0.9
    // closing charged pion cuts, minimum TPC cluster = 80, TPC dEdx pi = \pm 3 sigma, pi+pi- mass cut of 1.5, min pt charged pi = 100 MeV
    // closing neural pion cuts, 0.5 < M_gamma,gamma < 0.6
    // maxChi2 per cluster TPC <4, require TPC refit, DCA XY pT dependend 0.0182+0.0350/pt^1.01, DCA_Z = 3.0
    // timing cluster cut open
    cuts.AddCutHeavyMesonCalo("00010113","1111100047032230000","30a330709","0103503l00000000","0153503000000000"); // all of the above
    cuts.AddCutHeavyMesonCalo("00052113","1111100047032230000","30a330709","0103503l00000000","0153503000000000"); // all of the above
    cuts.AddCutHeavyMesonCalo("00062113","1111100047032230000","30a330709","0103503l00000000","0153503000000000"); // all of the above
    cuts.AddCutHeavyMesonCalo("00083113","1111100047032230000","30a330709","0103503l00000000","0153503000000000"); // all of the above
    cuts.AddCutHeavyMesonCalo("00085113","1111100047032230000","30a330709","0103503l00000000","0153503000000000"); // all of the above
  } else if( trainConfig == 203) {
    // same as 202 but only MB
    cuts.AddCutHeavyMesonCalo("00010113","1111100047032230000","30a330709","0103503l00000000","0153503000000000"); // 0.5-0.6 eta mass cut
    cuts.AddCutHeavyMesonCalo("00010113","1111100047032230000","30a330709","0103503m00000000","0153503000000000"); // 0.4-0.7 eta mass cut
  } else if( trainConfig == 204) {
    // same as 202 but with mass cut variations
    cuts.AddCutHeavyMesonCalo("00010113","1111100047032230000","30a330700","0103503l00000000","0153503000000000"); // pi+pi- mass cut of 10
    cuts.AddCutHeavyMesonCalo("00010113","1111100047032230000","30a330701","0103503l00000000","0153503000000000"); // pi+pi- mass cut of 1
    cuts.AddCutHeavyMesonCalo("00010113","1111100047032230000","30a330708","0103503l00000000","0153503000000000"); // pi+pi- mass cut of 0.85

  } else if ( trainConfig == 205 ) {
    cuts.AddCutHeavyMesonCalo("00010113","411791106f032220000","32c510700","0103603l000000d0","01535030000000d0"); // INT7
    cuts.AddCutHeavyMesonCalo("0008e113","411791106f032220000","32c510700","0103603l000000d0","01535030000000d0"); // EMC7
    cuts.AddCutHeavyMesonCalo("0008d113","411791106f032220000","32c510700","0103603l000000d0","01535030000000d0"); // EMC7
    cuts.AddCutHeavyMesonCalo("0009b113","411791106f032220000","32c510700","0103603l000000d0","01535030000000d0"); // EMC7
  } else if ( trainConfig == 206 ) { // no event mixing
    cuts.AddCutHeavyMesonCalo("00010113","411790106f032220000","32c510700","0103603l00000000","0453503000000000"); // INT7
    cuts.AddCutHeavyMesonCalo("0008e113","411790106f032220000","32c510700","0103603l00000000","0453503000000000"); // EMC7
    cuts.AddCutHeavyMesonCalo("0008d113","411790106f032220000","32c510700","0103603l00000000","0453503000000000"); // EMC7
    cuts.AddCutHeavyMesonCalo("0009b113","411790106f032220000","32c510700","0103603l00000000","0453503000000000"); // EMC7
  } else if ( trainConfig == 207 ) { // no event mixing only tiggers
    cuts.AddCutHeavyMesonCalo("0008e113","411790106f032220000","32c510700","0103603l00000000","0453503000000000"); // EMC7
    cuts.AddCutHeavyMesonCalo("0008d113","411790106f032220000","32c510700","0103603l00000000","0453503000000000"); // EMC7

    // PCM-PHOS
  } else if ( trainConfig == 250 ) { // INT7 + PHI7
    cuts.AddCutHeavyMesonCalo("00010113","24466190wa01cc00000","32c510700","0103603l00000010","0453503000000010"); // INT7
    cuts.AddCutHeavyMesonCalo("00062113","24466190wa01cc00000","32c510700","0103603l00000010","0453503000000010"); // PHI7
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //                                          D0 MESON
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  } else if( trainConfig == 300 ) {
    // everything open, min pt charged pi = 100 MeV
    cuts.AddCutHeavyMesonCalo("00010113","1111100047032230000","000010400","0103503a00000000","0103503000000000");
  } else if( trainConfig == 301 ) {
    // closing charged pion cuts, minimum TPC cluster = 80, TPC dEdx pi = \pm 3 sigma, min pt charged pi = 100 MeV
    cuts.AddCutHeavyMesonCalo("00010113","1111100047032230000","002010700","0103503a00000000","0103503000000000");
  } else if( trainConfig == 302) {
    // eta < 0.9
    // closing charged pion cuts, minimum TPC cluster = 80, TPC dEdx pi = \pm 3 sigma, pi+pi- mass cut of 0.65, min pt charged pi = 100 MeV
    // closing neural pion cuts, 0.1 < M_gamma,gamma < 0.15
    // maxChi2 per cluster TPC <4, require TPC refit, DCA XY pT dependend 0.0182+0.0350/pt^1.01, DCA_Z = 3.0
    // timing cluster cut open
    cuts.AddCutHeavyMesonCalo("00010113","1111100047032230000","30a330700","0103503400000000","0153503000000000"); // all of the above
    cuts.AddCutHeavyMesonCalo("00052113","1111100047032230000","30a330700","0103503400000000","0153503000000000"); // all of the above
    cuts.AddCutHeavyMesonCalo("00062113","1111100047032230000","30a330700","0103503400000000","0153503000000000"); // all of the above
    cuts.AddCutHeavyMesonCalo("00083113","1111100047032230000","30a330700","0103503400000000","0153503000000000"); // all of the above
    cuts.AddCutHeavyMesonCalo("00085113","1111100047032230000","30a330700","0103503400000000","0153503000000000"); // all of the above
  } else if( trainConfig == 303) {
    // same as 102 but only MB
    cuts.AddCutHeavyMesonCalo("00010113","1111100047032230000","30a330708","0103503400000000","0153503000000000"); // all of the above

    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //                                          OMEGA MESON pp 13 TeV
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // PHOS pp 13 TeV
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    //Standard Cuts
  } else if(trainConfig == 400)  { //Standard PHOS 13TeV MB
    cuts.AddCutHeavyMesonCalo("00010113","24466190sa01cc00000","32c51070a","0103103300000000","0453503000000000"); // INT7
  } else if(trainConfig == 401)  { //Standard PHOS 13TeV PHI7
    cuts.AddCutHeavyMesonCalo("00062113","24466190sa01cc00000","32c51070a","0103103300000000","0453503000000000"); // PHI7

    //Gamma Energy Cuts
  } else if(trainConfig == 403)  { //Standard PHOS 13TeV Trigger PHI7, GammaCut > 4GeV
    cuts.AddCutHeavyMesonCalo("00062113","24466190sa01cc00000","32c51070a","01031k3300000000","0453503000000000"); // PHI7, new Gamma Energy cut 4 GeV
  } else if(trainConfig == 404)  { //Standard PHOS 13TeV Trigger PHI7, GammaCut > 5GeV
    cuts.AddCutHeavyMesonCalo("00062113","24466190sa01cc00000","32c51070a","01031e3300000000","0453503000000000"); // PHI7, new Gamma Energy cut 5. GeV
  } else if(trainConfig == 405)  { //Standard PHOS 13TeV Trigger PHI7, GammaCut > 6GeV
    cuts.AddCutHeavyMesonCalo("00062113","24466190sa01cc00000","32c51070a","01031g3300000000","0453503000000000"); // PHI7, new Gamma Energy cut 6 GeV
  } else if(trainConfig == 406)  { //Standard PHOS 13TeV Trigger PHI7, GammaCut > 7.5GeV
    cuts.AddCutHeavyMesonCalo("00062113","24466190sa01cc00000","32c51070a","01031f3300000000","0453503000000000"); // PHI7, new Gamma Energy cut 7.5 GeV

    //QA Plots
  } else if(trainConfig == 415)  { //Standard PHOS 13TeV MB, no shared TPC clusters, Shared cluster Fraction =0
    cuts.AddCutHeavyMesonCalo("00010113","24466190sa01cc00000","32e51070a","0103103300000000","0453503000000000"); // INT7
  } else if(trainConfig == 416)  { //Standard PHOS 13TeV MB, no shared TPC clusters, Shared cluster Fraction <=0.4
    cuts.AddCutHeavyMesonCalo("00010113","24466190sa01cc00000","32f51070a","0103103300000000","0453503000000000"); // INT7


    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // EMC pp 13 TeV
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    //Standard Cuts
  } else if(trainConfig == 430)  { // EDC 13 TeV
    cuts.AddCutHeavyMesonCalo("00010113","411792109fe32220000","32c51070a","0103103100000000","0453503000000000"); // INT7
  } else if(trainConfig == 431)  { // EDC 13 TeV EG2
    cuts.AddCutHeavyMesonCalo("0008e113","411792109fe32220000","32c51070a","0103103100000000","0453503000000000"); // EG2
  } else if(trainConfig == 432)  { // EDC 13 TeV EG1
    cuts.AddCutHeavyMesonCalo("0008d113","411792109fe32220000","32c51070a","0103103100000000","0453503000000000"); // EG1

    //Gamma Energy Cuts EG2
  } else if(trainConfig == 433)  { //Standard EMCal 13TeV Trigger EG2, GammaCut > 4GeV
    cuts.AddCutHeavyMesonCalo("0008e113","411792109fe32220000","32c51070a","01031k3100000000","0453503000000000"); // EG2, new Gamma Energy cut 4 GeV
  } else if(trainConfig == 434)  { //Standard EMCal 13TeV Trigger EG2, GammaCut > 5GeV
    cuts.AddCutHeavyMesonCalo("0008e113","411792109fe32220000","32c51070a","01031e3100000000","0453503000000000"); // EG2, new Gamma Energy cut 5 GeV
  } else if(trainConfig == 435)  { //Standard EMCal 13TeV Trigger EG2, GammaCut > 6GeV
    cuts.AddCutHeavyMesonCalo("0008e113","411792109fe32220000","32c51070a","01031g3100000000","0453503000000000"); // EG2, new Gamma Energy cut 6 GeV
  } else if(trainConfig == 436)  { //Standard EMCal 13TeV Trigger EG2, GammaCut > 7.5GeV
    cuts.AddCutHeavyMesonCalo("0008e113","411792109fe32220000","32c51070a","01031f3100000000","0453503000000000"); // EG2, new Gamma Energy cut 7.5 GeV

    //Gamma Energy Cuts EG1
  } else if(trainConfig == 437)  { //EMCal Trig Pt Cut Variations EG1, GammaCut > 8GeV
    cuts.AddCutHeavyMesonCalo("0008d113","411792109fe32220000","32c51070a","01031l3100000000","0453503000000000"); // EG1, new Gamma Energy cut 8 GeV
  } else if(trainConfig == 438)  { //EMCal Trig Pt Cut Variations EG1, GammaCut > 9GeV
    cuts.AddCutHeavyMesonCalo("0008d113","411792109fe32220000","32c51070a","01031m3100000000","0453503000000000"); // EG1, new Gamma Energy cut 9 GeV
  } else if(trainConfig == 439)  { //EMCal Trig Pt Cut Variations EG1, GammaCut > 10GeV
    cuts.AddCutHeavyMesonCalo("0008d113","411792109fe32220000","32c51070a","01031h3100000000","0453503000000000"); // EG1, new Gamma Energy cut 10 GeV

    //no shared TPC clusters, Shared cluster Fraction =0
  } else if(trainConfig == 440)  { // EDC 13 TeV, Shared cluster Fraction =0
    cuts.AddCutHeavyMesonCalo("00010113","411792109fe32220000","32e51070a","0103103100000000","0453503000000000"); // INT7
  } else if(trainConfig == 441)  { //Standard EMCal 13TeV Trigger EG2, GammaCut > 6GeV, Shared cluster Fraction =0
    cuts.AddCutHeavyMesonCalo("0008e113","411792109fe32220000","32e51070a","01031g3100000000","0453503000000000"); // EG2, new Gamma Energy cut 6 GeV
  } else if(trainConfig == 442)  { //Standard EMCal 13TeV Trigger EG1, GammaCut > 10GeV, Shared cluster Fraction =0
    cuts.AddCutHeavyMesonCalo("0008d113","411792109fe32220000","32e51070a","01031h3100000000","0453503000000000"); // EG1, new Gamma Energy cut 10 GeV

    //no shared TPC clusters, Shared cluster Fraction <=0.4
  } else if(trainConfig == 443)  { // EDC 13 TeV, Shared cluster Fraction <=0.4
    cuts.AddCutHeavyMesonCalo("00010113","411792109fe32220000","32f51070a","0103103100000000","0453503000000000"); // INT7
  } else if(trainConfig == 444)  { //Standard EMCal 13TeV Trigger EG2, GammaCut > 6GeV, Shared cluster Fraction <=0.4
    cuts.AddCutHeavyMesonCalo("0008e113","411792109fe32220000","32f51070a","01031g3100000000","0453503000000000"); // EG2, new Gamma Energy cut 6 GeV
  } else if(trainConfig == 445)  { //Standard EMCal 13TeV Trigger EG1, GammaCut > 10GeV, Shared cluster Fraction <=0.4
    cuts.AddCutHeavyMesonCalo("0008d113","411792109fe32220000","32f51070a","01031h3100000000","0453503000000000"); // EG1, new Gamma Energy cut 10 GeV

    // Variations on 13 TeV for 7 TeV systematics
    // EMC (without nonlin)
  } else if(trainConfig == 450)  { // std
    cuts.AddCutHeavyMesonCalo("00010113","411791106f032220000","32c51070a","0103603600000000","0153503000000000"); // 2 sigma mass cut
  } else if(trainConfig == 451)  { // mass window cut
    cuts.AddCutHeavyMesonCalo("00010113","411791106f032220000","32c51070a","0103603b00000000","0153503000000000"); // 1 sigma mass cut
    cuts.AddCutHeavyMesonCalo("00010113","411791106f032220000","32c51070a","0103603c00000000","0153503000000000"); // 3 sigma mass cut
    cuts.AddCutHeavyMesonCalo("00010113","411791106f032220000","32c51070a","0103603d00000000","0153503000000000"); // 4 sigma mass cut
  } else if(trainConfig == 452)  { // background description
    cuts.AddCutHeavyMesonCalo("00010113","411791106f032220000","32c51070a","0103603600000000","0a53503000000000"); // likesign
    cuts.AddCutHeavyMesonCalo("00010113","411791106f032220000","32c51070a","0103603600000000","0d53503000000000"); // sideband mixing
  } else if(trainConfig == 453)  { // track matching
    cuts.AddCutHeavyMesonCalo("00010113","4117911061032220000","32c51070a","0103603600000000","0153503000000000");
    cuts.AddCutHeavyMesonCalo("00010113","4117911063032220000","32c51070a","0103603600000000","0153503000000000");
    cuts.AddCutHeavyMesonCalo("00010113","4117911065032220000","32c51070a","0103603600000000","0153503000000000");
    cuts.AddCutHeavyMesonCalo("00010113","4117911066032220000","32c51070a","0103603600000000","0153503000000000");

    //Standard Cuts of Pi0 Analysis: ("00010113","411792106fe32220000","0r631031000000d0")
    //MesonCut r63==Background->ignored, d==OpeningAngle for Background->ignored
  } else if(trainConfig == 454)  { // EDC 13 TeV
    cuts.AddCutHeavyMesonCalo("00010113","411792106fe32220000","32c51070a","0103103100000000","0453503000000000"); // INT7
  } else if(trainConfig == 455)  { // EDC 13 TeV EG2
    cuts.AddCutHeavyMesonCalo("0008e113","411792106fe32220000","32c51070a","0103103100000000","0453503000000000"); // EG2
  } else if(trainConfig == 456)  { // EDC 13 TeV EG1
    cuts.AddCutHeavyMesonCalo("0008d113","411792106fe32220000","32c51070a","0103103100000000","0453503000000000"); // EG1

    //Charged Pion Mass Cut <600 MeV
  } else if(trainConfig == 460)  { // EDC 13 TeV
    cuts.AddCutHeavyMesonCalo("00010113","411792109fe32220000","32c51070b","0103103100000000","0453503000000000"); // INT7
  } else if(trainConfig == 461)  { // EDC 13 TeV EG2
    cuts.AddCutHeavyMesonCalo("0008e113","411792109fe32220000","32c51070b","0103103100000000","0453503000000000"); // EG2
  } else if(trainConfig == 462)  { // EDC 13 TeV EG1
    cuts.AddCutHeavyMesonCalo("0008d113","411792109fe32220000","32c51070b","0103103100000000","0453503000000000"); // EG1
    //Charged Pion Mass Cut <850 MeV, Neutral Charged Pion Cut 1000< MeV
  } else if(trainConfig == 463)  { // EDC 13 TeV
    cuts.AddCutHeavyMesonCalo("00010113","411792109fe32220000","32c51070c","0103103100000000","0453503000000000"); // INT7
  } else if(trainConfig == 464)  { // EDC 13 TeV EG2
    cuts.AddCutHeavyMesonCalo("0008e113","411792109fe32220000","32c51070c","0103103100000000","0453503000000000"); // EG2
  } else if(trainConfig == 465)  { // EDC 13 TeV EG1
    cuts.AddCutHeavyMesonCalo("0008d113","411792109fe32220000","32c51070c","0103103100000000","0453503000000000"); // EG1
    //Charged Pion Mass Cut <850 MeV, Neutral Charged Pion Cut <850 MeV
  } else if(trainConfig == 466)  { // EDC 13 TeV
    cuts.AddCutHeavyMesonCalo("00010113","411792109fe32220000","32c51070d","0103103100000000","0453503000000000"); // INT7
  } else if(trainConfig == 467)  { // EDC 13 TeV EG2
    cuts.AddCutHeavyMesonCalo("0008e113","411792109fe32220000","32c51070d","0103103100000000","0453503000000000"); // EG2
  } else if(trainConfig == 468)  { // EDC 13 TeV EG1
    cuts.AddCutHeavyMesonCalo("0008d113","411792109fe32220000","32c51070d","0103103100000000","0453503000000000"); // EG1
    //Charged Pion Mass Cut <600 MeV, Neutral Charged Pion Cut <600 MeV
  } else if(trainConfig == 469)  { // EDC 13 TeV
    cuts.AddCutHeavyMesonCalo("00010113","411792109fe32220000","32c51070e","0103103100000000","0453503000000000"); // INT7
  } else if(trainConfig == 470)  { // EDC 13 TeV EG2
    cuts.AddCutHeavyMesonCalo("0008e113","411792109fe32220000","32c51070e","0103103100000000","0453503000000000"); // EG2
  } else if(trainConfig == 471)  { // EDC 13 TeV EG1
    cuts.AddCutHeavyMesonCalo("0008d113","411792109fe32220000","32c51070e","0103103100000000","0453503000000000"); // EG1

    // Variations on 13 TeV for 7 TeV systematics
    // PHOS (without nonlin)
  } else if(trainConfig == 500)  { //Standard PHOS 13TeV Standard
    cuts.AddCutHeavyMesonCalo("00010113","24466510ga012200000","32c51070a","0103603800000000","0153503000000000"); // INT7
  } else if(trainConfig == 501)  { // mass window cut
    cuts.AddCutHeavyMesonCalo("00010113","24466510ga012200000","32c51070a","0103603h00000000","0153503000000000"); // 1 sigma
    cuts.AddCutHeavyMesonCalo("00010113","24466510ga012200000","32c51070a","0103603i00000000","0153503000000000"); // 3 sigma
    cuts.AddCutHeavyMesonCalo("00010113","24466510ga012200000","32c51070a","0103603j00000000","0153503000000000"); // 4 sigma
  } else if(trainConfig == 502)  { // background description
    cuts.AddCutHeavyMesonCalo("00010113","24466510ga012200000","32c51070a","0103603800000000","0a53503000000000"); // likesign
    cuts.AddCutHeavyMesonCalo("00010113","24466510ga012200000","32c51070a","0103603800000000","0d53503000000000"); // sideband mixing
  } else if(trainConfig == 503)  { // track matching
    cuts.AddCutHeavyMesonCalo("00010113","24466510g1012200000","32c51070a","0103603800000000","0153503000000000");
    cuts.AddCutHeavyMesonCalo("00010113","24466510g3012200000","32c51070a","0103603800000000","0153503000000000");
    cuts.AddCutHeavyMesonCalo("00010113","24466510g5012200000","32c51070a","0103603800000000","0153503000000000");
    cuts.AddCutHeavyMesonCalo("00010113","24466510g6012200000","32c51070a","0103603800000000","0153503000000000");
    cuts.AddCutHeavyMesonCalo("00010113","24466510g7012200000","32c51070a","0103603800000000","0153503000000000");
    cuts.AddCutHeavyMesonCalo("00010113","24466510g9012200000","32c51070a","0103603800000000","0153503000000000");


  } else if(trainConfig == 510)  { // EDC 13 TeV, nSigma Variation: -2.5,2.5
    cuts.AddCutHeavyMesonCalo("00010113","411792109fe32220000","32c51090a","0103103100000000","0453503000000000"); // INT7
  } else if(trainConfig == 511)  { // EDC 13 TeV, nSigma Variation: -3.5,3.5
    cuts.AddCutHeavyMesonCalo("00010113","411792109fe32220000","32c510a0a","0103103100000000","0453503000000000"); // INT7
  } else if(trainConfig == 512)  { // EDC 13 TeV, nSigma Variation: -2.0,2.0
    cuts.AddCutHeavyMesonCalo("00010113","411792109fe32220000","32c510b0a","0103103100000000","0453503000000000"); // INT7
  } else if(trainConfig == 513)  { // EDC 13 TeV, nSigma Variation: -4.0,4.0
    cuts.AddCutHeavyMesonCalo("00010113","411792109fe32220000","32c51050a","0103103100000000","0453503000000000"); // INT7
  } else if(trainConfig == 514)  { // EDC 13 TeV, nSigma Variation: -5.0,5.0
    cuts.AddCutHeavyMesonCalo("00010113","411792109fe32220000","32c51030a","0103103100000000","0453503000000000"); // INT7

  } else if(trainConfig == 520)  { // EDC 13 TeV, min opening angle 5mrad
    cuts.AddCutHeavyMesonCalo("00010113","411792109fe32220000","32c510b0a","0103103100000010","0453503000000000"); // INT7

  } else if(trainConfig == 525)  { // EDC 13 TeV, change of mass cut u // EMC-EMC - 1 sigma - gamma selection
    cuts.AddCutHeavyMesonCalo("00010113","411792109fe32220000","32c51070a","0103103u00000000","0453503000000000"); // INT7
  } else if(trainConfig == 526)  { // EDC 13 TeV, change of mass cut v // EMC-EMC - 3 sigma - gamma selection
    cuts.AddCutHeavyMesonCalo("00010113","411792109fe32220000","32c51070a","0103103v00000000","0453503000000000"); // INT7
  } else if(trainConfig == 527)  { // EDC 13 TeV, change of mass cut w // EMC-EMC - 4 sigma - gamma selection
    cuts.AddCutHeavyMesonCalo("00010113","411792109fe32220000","32c51070a","0103103w00000000","0453503000000000"); // INT7
  } else if(trainConfig == 528)  { // EDC 13 TeV, change of mass cut x // EMC-EMC - 2 sigma - gamma selection
    cuts.AddCutHeavyMesonCalo("00010113","411792109fe32220000","32c51070a","0103103x00000000","0453503000000000"); // INT7



  } else if(trainConfig == 600)  { //Standard PHOS 13TeV MB
    cuts.AddCutHeavyMesonCalo("00010113","24466190sa01cc00000","32c51070a","0103103300000000","0153503000000000"); // INT7
  } else if(trainConfig == 601)  { //Standard PHOS 13TeV PHI7
    cuts.AddCutHeavyMesonCalo("00062113","24466190sa01cc00000","32c51070a","0103103300000000","0153503000000000"); // PHI7

    //Gamma Energy Cuts
  } else if(trainConfig == 603)  { //Standard PHOS 13TeV Trigger PHI7, GammaCut > 4GeV
    cuts.AddCutHeavyMesonCalo("00062113","24466190sa01cc00000","32c51070a","01031k3300000000","0153503000000000"); // PHI7, new Gamma Energy cut 4 GeV
  } else if(trainConfig == 604)  { //Standard PHOS 13TeV Trigger PHI7, GammaCut > 5GeV
    cuts.AddCutHeavyMesonCalo("00062113","24466190sa01cc00000","32c51070a","01031e3300000000","0153503000000000"); // PHI7, new Gamma Energy cut 5. GeV
  } else if(trainConfig == 605)  { //Standard PHOS 13TeV Trigger PHI7, GammaCut > 6GeV
    cuts.AddCutHeavyMesonCalo("00062113","24466190sa01cc00000","32c51070a","01031g3300000000","0153503000000000"); // PHI7, new Gamma Energy cut 6 GeV
  } else if(trainConfig == 606)  { //Standard PHOS 13TeV Trigger PHI7, GammaCut > 7.5GeV
    cuts.AddCutHeavyMesonCalo("00062113","24466190sa01cc00000","32c51070a","01031f3300000000","0153503000000000"); // PHI7, new Gamma Energy cut 7.5 GeV

    //QA Plots
  } else if(trainConfig == 615)  { //Standard PHOS 13TeV MB, no shared TPC clusters, Shared cluster Fraction =0
    cuts.AddCutHeavyMesonCalo("00010113","24466190sa01cc00000","32e51070a","0103103300000000","0153503000000000"); // INT7
  } else if(trainConfig == 616)  { //Standard PHOS 13TeV MB, no shared TPC clusters, Shared cluster Fraction <=0.4
    cuts.AddCutHeavyMesonCalo("00010113","24466190sa01cc00000","32f51070a","0103103300000000","0153503000000000"); // INT7


    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // EMC pp 13 TeV
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    //Standard Cuts
  } else if(trainConfig == 630)  { // EDC 13 TeV
    cuts.AddCutHeavyMesonCalo("00010113","411792109fe32220000","32c51070a","0103103100000000","0153503000000000"); // INT7
  } else if(trainConfig == 631)  { // EDC 13 TeV EG2
    cuts.AddCutHeavyMesonCalo("0008e113","411792109fe32220000","32c51070a","0103103100000000","0153503000000000"); // EG2
  } else if(trainConfig == 632)  { // EDC 13 TeV EG1
    cuts.AddCutHeavyMesonCalo("0008d113","411792109fe32220000","32c51070a","0103103100000000","0153503000000000"); // EG1

    //Gamma Energy Cuts EG2
  } else if(trainConfig == 633)  { //Standard EMCal 13TeV Trigger EG2, GammaCut > 4GeV
    cuts.AddCutHeavyMesonCalo("0008e113","411792109fe32220000","32c51070a","01031k3100000000","0153503000000000"); // EG2, new Gamma Energy cut 4 GeV
  } else if(trainConfig == 634)  { //Standard EMCal 13TeV Trigger EG2, GammaCut > 5GeV
    cuts.AddCutHeavyMesonCalo("0008e113","411792109fe32220000","32c51070a","01031e3100000000","0153503000000000"); // EG2, new Gamma Energy cut 5 GeV
  } else if(trainConfig == 635)  { //Standard EMCal 13TeV Trigger EG2, GammaCut > 6GeV
    cuts.AddCutHeavyMesonCalo("0008e113","411792109fe32220000","32c51070a","01031g3100000000","0153503000000000"); // EG2, new Gamma Energy cut 6 GeV
  } else if(trainConfig == 636)  { //Standard EMCal 13TeV Trigger EG2, GammaCut > 7.5GeV
    cuts.AddCutHeavyMesonCalo("0008e113","411792109fe32220000","32c51070a","01031f3100000000","0153503000000000"); // EG2, new Gamma Energy cut 7.5 GeV

    //Gamma Energy Cuts EG1
  } else if(trainConfig == 637)  { //EMCal Trig Pt Cut Variations EG1, GammaCut > 8GeV
    cuts.AddCutHeavyMesonCalo("0008d113","411792109fe32220000","32c51070a","01031l3100000000","0153503000000000"); // EG1, new Gamma Energy cut 8 GeV
  } else if(trainConfig == 638)  { //EMCal Trig Pt Cut Variations EG1, GammaCut > 9GeV
    cuts.AddCutHeavyMesonCalo("0008d113","411792109fe32220000","32c51070a","01031m3100000000","0153503000000000"); // EG1, new Gamma Energy cut 9 GeV
  } else if(trainConfig == 639)  { //EMCal Trig Pt Cut Variations EG1, GammaCut > 10GeV
    cuts.AddCutHeavyMesonCalo("0008d113","411792109fe32220000","32c51070a","01031h3100000000","0153503000000000"); // EG1, new Gamma Energy cut 10 GeV

    //no shared TPC clusters, Shared cluster Fraction =0
  } else if(trainConfig == 640)  { // EDC 13 TeV, Shared cluster Fraction =0
    cuts.AddCutHeavyMesonCalo("00010113","411792109fe32220000","32e51070a","0103103100000000","0153503000000000"); // INT7
  } else if(trainConfig == 641)  { //Standard EMCal 13TeV Trigger EG2, GammaCut > 6GeV, Shared cluster Fraction =0
    cuts.AddCutHeavyMesonCalo("0008e113","411792109fe32220000","32e51070a","01031g3100000000","0153503000000000"); // EG2, new Gamma Energy cut 6 GeV
  } else if(trainConfig == 642)  { //Standard EMCal 13TeV Trigger EG1, GammaCut > 10GeV, Shared cluster Fraction =0
    cuts.AddCutHeavyMesonCalo("0008d113","411792109fe32220000","32e51070a","01031h3100000000","0153503000000000"); // EG1, new Gamma Energy cut 10 GeV

    //no shared TPC clusters, Shared cluster Fraction <=0.4
  } else if(trainConfig == 643)  { // EDC 13 TeV, Shared cluster Fraction <=0.4
    cuts.AddCutHeavyMesonCalo("00010113","411792109fe32220000","32f51070a","0103103100000000","0153503000000000"); // INT7
  } else if(trainConfig == 644)  { //Standard EMCal 13TeV Trigger EG2, GammaCut > 6GeV, Shared cluster Fraction <=0.4
    cuts.AddCutHeavyMesonCalo("0008e113","411792109fe32220000","32f51070a","01031g3100000000","0153503000000000"); // EG2, new Gamma Energy cut 6 GeV
  } else if(trainConfig == 645)  { //Standard EMCal 13TeV Trigger EG1, GammaCut > 10GeV, Shared cluster Fraction <=0.4
    cuts.AddCutHeavyMesonCalo("0008d113","411792109fe32220000","32f51070a","01031h3100000000","0153503000000000"); // EG1, new Gamma Energy cut 10 GeV

    //Charged Pion Mass Cut <600 MeV
  } else if(trainConfig == 646)  { // EDC 13 TeV
    cuts.AddCutHeavyMesonCalo("00010113","411792109fe32220000","32c51070b","0103103100000000","0153503000000000"); // INT7
  } else if(trainConfig == 647)  { // EDC 13 TeV EG2
    cuts.AddCutHeavyMesonCalo("0008e113","411792109fe32220000","32c51070b","0103103100000000","0153503000000000"); // EG2
  } else if(trainConfig == 648)  { // EDC 13 TeV EG1
    cuts.AddCutHeavyMesonCalo("0008d113","411792109fe32220000","32c51070b","0103103100000000","0153503000000000"); // EG1

    //Standard Cuts of Pi0 Analysis: ("00010113","411792106fe32220000","0r631031000000d0")
    //MesonCut r63==Background->ignored, d==OpeningAngle for Background->ignored
  } else if(trainConfig == 654)  { // EDC 13 TeV
    cuts.AddCutHeavyMesonCalo("00010113","411792106fe32220000","32c51070a","0103103100000000","0153503000000000"); // INT7
  } else if(trainConfig == 655)  { // EDC 13 TeV EG2
    cuts.AddCutHeavyMesonCalo("0008e113","411792106fe32220000","32c51070a","0103103100000000","0153503000000000"); // EG2
  } else if(trainConfig == 656)  { // EDC 13 TeV EG1
    cuts.AddCutHeavyMesonCalo("0008d113","411792106fe32220000","32c51070a","0103103100000000","0153503000000000"); // EG1

    //Charged Pion Mass Cut <600 MeV
  } else if(trainConfig == 660)  { // EDC 13 TeV
    cuts.AddCutHeavyMesonCalo("00010113","411792109fe32220000","32c51070b","0103103100000000","0153503000000000"); // INT7
  } else if(trainConfig == 661)  { // EDC 13 TeV EG2
    cuts.AddCutHeavyMesonCalo("0008e113","411792109fe32220000","32c51070b","0103103100000000","0153503000000000"); // EG2
  } else if(trainConfig == 662)  { // EDC 13 TeV EG1
    cuts.AddCutHeavyMesonCalo("0008d113","411792109fe32220000","32c51070b","0103103100000000","0153503000000000"); // EG1
    //Charged Pion Mass Cut <850 MeV, Neutral Charged Pion Cut 1000< MeV
  } else if(trainConfig == 663)  { // EDC 13 TeV
    cuts.AddCutHeavyMesonCalo("00010113","411792109fe32220000","32c51070c","0103103100000000","0153503000000000"); // INT7
  } else if(trainConfig == 664)  { // EDC 13 TeV EG2
    cuts.AddCutHeavyMesonCalo("0008e113","411792109fe32220000","32c51070c","0103103100000000","0153503000000000"); // EG2
  } else if(trainConfig == 665)  { // EDC 13 TeV EG1
    cuts.AddCutHeavyMesonCalo("0008d113","411792109fe32220000","32c51070c","0103103100000000","0153503000000000"); // EG1
    //Charged Pion Mass Cut <850 MeV, Neutral Charged Pion Cut <850 MeV
  } else if(trainConfig == 666)  { // EDC 13 TeV
    cuts.AddCutHeavyMesonCalo("00010113","411792109fe32220000","32c51070d","0103103100000000","0153503000000000"); // INT7
  } else if(trainConfig == 667)  { // EDC 13 TeV EG2
    cuts.AddCutHeavyMesonCalo("0008e113","411792109fe32220000","32c51070d","0103103100000000","0153503000000000"); // EG2
  } else if(trainConfig == 668)  { // EDC 13 TeV EG1
    cuts.AddCutHeavyMesonCalo("0008d113","411792109fe32220000","32c51070d","0103103100000000","0153503000000000"); // EG1
    //Charged Pion Mass Cut <600 MeV, Neutral Charged Pion Cut <600 MeV
  } else if(trainConfig == 669)  { // EDC 13 TeV
    cuts.AddCutHeavyMesonCalo("00010113","411792109fe32220000","32c51070e","0103103100000000","0153503000000000"); // INT7
  } else if(trainConfig == 670)  { // EDC 13 TeV EG2
    cuts.AddCutHeavyMesonCalo("0008e113","411792109fe32220000","32c51070e","0103103100000000","0153503000000000"); // EG2
  } else if(trainConfig == 671)  { // EDC 13 TeV EG1
    cuts.AddCutHeavyMesonCalo("0008d113","411792109fe32220000","32c51070e","0103103100000000","0153503000000000"); // EG1

  } else if(trainConfig == 700)  { //Standard PHOS 13TeV MB
    cuts.AddCutHeavyMesonCalo("00010113","24466190sa01cc00000","32c51070a","0103103300000000","0a53503000000000"); // INT7
  } else if(trainConfig == 701)  { //Standard PHOS 13TeV PHI7
    cuts.AddCutHeavyMesonCalo("00062113","24466190sa01cc00000","32c51070a","0103103300000000","0a53503000000000"); // PHI7

    //Gamma Energy Cuts
  } else if(trainConfig == 703)  { //Standard PHOS 13TeV Trigger PHI7, GammaCut > 4GeV
    cuts.AddCutHeavyMesonCalo("00062113","24466190sa01cc00000","32c51070a","01031k3300000000","0a53503000000000"); // PHI7, new Gamma Energy cut 4 GeV
  } else if(trainConfig == 704)  { //Standard PHOS 13TeV Trigger PHI7, GammaCut > 5GeV
    cuts.AddCutHeavyMesonCalo("00062113","24466190sa01cc00000","32c51070a","01031e3300000000","0a53503000000000"); // PHI7, new Gamma Energy cut 5. GeV
  } else if(trainConfig == 705)  { //Standard PHOS 13TeV Trigger PHI7, GammaCut > 6GeV
    cuts.AddCutHeavyMesonCalo("00062113","24466190sa01cc00000","32c51070a","01031g3300000000","0a53503000000000"); // PHI7, new Gamma Energy cut 6 GeV
  } else if(trainConfig == 706)  { //Standard PHOS 13TeV Trigger PHI7, GammaCut > 7.5GeV
    cuts.AddCutHeavyMesonCalo("00062113","24466190sa01cc00000","32c51070a","01031f3300000000","0a53503000000000"); // PHI7, new Gamma Energy cut 7.5 GeV

    //QA Plots
  } else if(trainConfig == 715)  { //Standard PHOS 13TeV MB, no shared TPC clusters, Shared cluster Fraction =0
    cuts.AddCutHeavyMesonCalo("00010113","24466190sa01cc00000","32e51070a","0103103300000000","0a53503000000000"); // INT7
  } else if(trainConfig == 716)  { //Standard PHOS 13TeV MB, no shared TPC clusters, Shared cluster Fraction <=0.4
    cuts.AddCutHeavyMesonCalo("00010113","24466190sa01cc00000","32f51070a","0103103300000000","0a53503000000000"); // INT7


    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // EMC pp 13 TeV
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    //Standard Cuts
  } else if(trainConfig == 730)  { // EDC 13 TeV
    cuts.AddCutHeavyMesonCalo("00010113","411792109fe32220000","32c51070a","0103103100000000","0a53503000000000"); // INT7
  } else if(trainConfig == 731)  { // EDC 13 TeV EG2
    cuts.AddCutHeavyMesonCalo("0008e113","411792109fe32220000","32c51070a","0103103100000000","0a53503000000000"); // EG2
  } else if(trainConfig == 732)  { // EDC 13 TeV EG1
    cuts.AddCutHeavyMesonCalo("0008d113","411792109fe32220000","32c51070a","0103103100000000","0a53503000000000"); // EG1

    //Gamma Energy Cuts EG2
  } else if(trainConfig == 733)  { //Standard EMCal 13TeV Trigger EG2, GammaCut > 4GeV
    cuts.AddCutHeavyMesonCalo("0008e113","411792109fe32220000","32c51070a","01031k3100000000","0a53503000000000"); // EG2, new Gamma Energy cut 4 GeV
  } else if(trainConfig == 734)  { //Standard EMCal 13TeV Trigger EG2, GammaCut > 5GeV
    cuts.AddCutHeavyMesonCalo("0008e113","411792109fe32220000","32c51070a","01031e3100000000","0a53503000000000"); // EG2, new Gamma Energy cut 5 GeV
  } else if(trainConfig == 735)  { //Standard EMCal 13TeV Trigger EG2, GammaCut > 6GeV
    cuts.AddCutHeavyMesonCalo("0008e113","411792109fe32220000","32c51070a","01031g3100000000","0a53503000000000"); // EG2, new Gamma Energy cut 6 GeV
  } else if(trainConfig == 736)  { //Standard EMCal 13TeV Trigger EG2, GammaCut > 7.5GeV
    cuts.AddCutHeavyMesonCalo("0008e113","411792109fe32220000","32c51070a","01031f3100000000","0a53503000000000"); // EG2, new Gamma Energy cut 7.5 GeV

    //Gamma Energy Cuts EG1
  } else if(trainConfig == 737)  { //EMCal Trig Pt Cut Variations EG1, GammaCut > 8GeV
    cuts.AddCutHeavyMesonCalo("0008d113","411792109fe32220000","32c51070a","01031l3100000000","0a53503000000000"); // EG1, new Gamma Energy cut 8 GeV
  } else if(trainConfig == 738)  { //EMCal Trig Pt Cut Variations EG1, GammaCut > 9GeV
    cuts.AddCutHeavyMesonCalo("0008d113","411792109fe32220000","32c51070a","01031m3100000000","0a53503000000000"); // EG1, new Gamma Energy cut 9 GeV
  } else if(trainConfig == 739)  { //EMCal Trig Pt Cut Variations EG1, GammaCut > 10GeV
    cuts.AddCutHeavyMesonCalo("0008d113","411792109fe32220000","32c51070a","01031h3100000000","0a53503000000000"); // EG1, new Gamma Energy cut 10 GeV

    //no shared TPC clusters, Shared cluster Fraction =0
  } else if(trainConfig == 740)  { // EDC 13 TeV, Shared cluster Fraction =0
    cuts.AddCutHeavyMesonCalo("00010113","411792109fe32220000","32e51070a","0103103100000000","0a53503000000000"); // INT7
  } else if(trainConfig == 741)  { //Standard EMCal 13TeV Trigger EG2, GammaCut > 6GeV, Shared cluster Fraction =0
    cuts.AddCutHeavyMesonCalo("0008e113","411792109fe32220000","32e51070a","01031g3100000000","0a53503000000000"); // EG2, new Gamma Energy cut 6 GeV
  } else if(trainConfig == 742)  { //Standard EMCal 13TeV Trigger EG1, GammaCut > 10GeV, Shared cluster Fraction =0
    cuts.AddCutHeavyMesonCalo("0008d113","411792109fe32220000","32e51070a","01031h3100000000","0a53503000000000"); // EG1, new Gamma Energy cut 10 GeV

    //no shared TPC clusters, Shared cluster Fraction <=0.4
  } else if(trainConfig == 743)  { // EDC 13 TeV, Shared cluster Fraction <=0.4
    cuts.AddCutHeavyMesonCalo("00010113","411792109fe32220000","32f51070a","0103103100000000","0a53503000000000"); // INT7
  } else if(trainConfig == 744)  { //Standard EMCal 13TeV Trigger EG2, GammaCut > 6GeV, Shared cluster Fraction <=0.4
    cuts.AddCutHeavyMesonCalo("0008e113","411792109fe32220000","32f51070a","01031g3100000000","0a53503000000000"); // EG2, new Gamma Energy cut 6 GeV
  } else if(trainConfig == 745)  { //Standard EMCal 13TeV Trigger EG1, GammaCut > 10GeV, Shared cluster Fraction <=0.4
    cuts.AddCutHeavyMesonCalo("0008d113","411792109fe32220000","32f51070a","01031h3100000000","0a53503000000000"); // EG1, new Gamma Energy cut 10 GeV

    //Charged Pion Mass Cut <600 MeV
  } else if(trainConfig == 746)  { // EDC 13 TeV
    cuts.AddCutHeavyMesonCalo("00010113","411792109fe32220000","32c51070b","0103103100000000","0a53503000000000"); // INT7
  } else if(trainConfig == 747)  { // EDC 13 TeV EG2
    cuts.AddCutHeavyMesonCalo("0008e113","411792109fe32220000","32c51070b","0103103100000000","0a53503000000000"); // EG2
  } else if(trainConfig == 748)  { // EDC 13 TeV EG1
    cuts.AddCutHeavyMesonCalo("0008d113","411792109fe32220000","32c51070b","0103103100000000","0a53503000000000"); // EG1

    //Standard Cuts of Pi0 Analysis: ("00010113","411792106fe32220000","0r631031000000d0")
    //MesonCut r63==Background->ignored, d==OpeningAngle for Background->ignored
  } else if(trainConfig == 754)  { // EDC 13 TeV
    cuts.AddCutHeavyMesonCalo("00010113","411792106fe32220000","32c51070a","0103103100000000","0a53503000000000"); // INT7
  } else if(trainConfig == 755)  { // EDC 13 TeV EG2
    cuts.AddCutHeavyMesonCalo("0008e113","411792106fe32220000","32c51070a","0103103100000000","0a53503000000000"); // EG2
  } else if(trainConfig == 756)  { // EDC 13 TeV EG1
    cuts.AddCutHeavyMesonCalo("0008d113","411792106fe32220000","32c51070a","0103103100000000","0a53503000000000"); // EG1



    //Charged Pion Mass Cut <600 MeV
  } else if(trainConfig == 760)  { // EDC 13 TeV
    cuts.AddCutHeavyMesonCalo("00010113","411792109fe32220000","32c51070b","0103103100000000","0a53503000000000"); // INT7
  } else if(trainConfig == 761)  { // EDC 13 TeV EG2
    cuts.AddCutHeavyMesonCalo("0008e113","411792109fe32220000","32c51070b","0103103100000000","0a53503000000000"); // EG2
  } else if(trainConfig == 762)  { // EDC 13 TeV EG1
    cuts.AddCutHeavyMesonCalo("0008d113","411792109fe32220000","32c51070b","0103103100000000","0a53503000000000"); // EG1
    //Charged Pion Mass Cut <850 MeV, Neutral Charged Pion Cut 1000< MeV
  } else if(trainConfig == 763)  { // EDC 13 TeV
    cuts.AddCutHeavyMesonCalo("00010113","411792109fe32220000","32c51070c","0103103100000000","0a53503000000000"); // INT7
  } else if(trainConfig == 764)  { // EDC 13 TeV EG2
    cuts.AddCutHeavyMesonCalo("0008e113","411792109fe32220000","32c51070c","0103103100000000","0a53503000000000"); // EG2
  } else if(trainConfig == 765)  { // EDC 13 TeV EG1
    cuts.AddCutHeavyMesonCalo("0008d113","411792109fe32220000","32c51070c","0103103100000000","0a53503000000000"); // EG1
    //Charged Pion Mass Cut <850 MeV, Neutral Charged Pion Cut <850 MeV
  } else if(trainConfig == 766)  { // EDC 13 TeV
    cuts.AddCutHeavyMesonCalo("00010113","411792109fe32220000","32c51070d","0103103100000000","0a53503000000000"); // INT7
  } else if(trainConfig == 767)  { // EDC 13 TeV EG2
    cuts.AddCutHeavyMesonCalo("0008e113","411792109fe32220000","32c51070d","0103103100000000","0a53503000000000"); // EG2
  } else if(trainConfig == 768)  { // EDC 13 TeV EG1
    cuts.AddCutHeavyMesonCalo("0008d113","411792109fe32220000","32c51070d","0103103100000000","0a53503000000000"); // EG1
    //Charged Pion Mass Cut <600 MeV, Neutral Charged Pion Cut <600 MeV
  } else if(trainConfig == 769)  { // EDC 13 TeV
    cuts.AddCutHeavyMesonCalo("00010113","411792109fe32220000","32c51070e","0103103100000000","0a53503000000000"); // INT7
  } else if(trainConfig == 770)  { // EDC 13 TeV EG2
    cuts.AddCutHeavyMesonCalo("0008e113","411792109fe32220000","32c51070e","0103103100000000","0a53503000000000"); // EG2
  } else if(trainConfig == 771)  { // EDC 13 TeV EG1
    cuts.AddCutHeavyMesonCalo("0008d113","411792109fe32220000","32c51070e","0103103100000000","0a53503000000000"); // EG1

    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // EMC pp 13 TeV Sideband Mixing
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //Standard Cuts
  } else if(trainConfig == 830)  { // EDC 13 TeV
    cuts.AddCutHeavyMesonCalo("00010113","411792109fe32220000","32c51070a","0103103100000000","0d53503000000000"); // INT7, Sideband, both sides
  } else if(trainConfig == 831)  { // EDC 13 TeV EG2
    cuts.AddCutHeavyMesonCalo("0008e113","411792109fe32220000","32c51070a","0103103100000000","0d53503000000000"); // EG2, Sideband, both sides
  } else if(trainConfig == 832)  { // EDC 13 TeV EG1
    cuts.AddCutHeavyMesonCalo("0008d113","411792109fe32220000","32c51070a","0103103100000000","0d53503000000000"); // EG1, Sideband, both sides

  } else if(trainConfig == 833)  { // EDC 13 TeV
    cuts.AddCutHeavyMesonCalo("00010113","411792109fe32220000","32c51070a","0103103100000000","0b53503000000000"); // INT7, Sideband, right side
  } else if(trainConfig == 834)  { // EDC 13 TeV
    cuts.AddCutHeavyMesonCalo("00010113","411792109fe32220000","32c51070a","0103103100000000","0c53503000000000"); // INT7, Sideband, left side



    //Charged Pion Mass Cut <600 MeV
  } else if(trainConfig == 860)  { // EDC 13 TeV
    cuts.AddCutHeavyMesonCalo("00010113","411792109fe32220000","32c51070b","0103103100000000","0b53503000000000"); // INT7
  } else if(trainConfig == 861)  { // EDC 13 TeV EG2
    cuts.AddCutHeavyMesonCalo("0008e113","411792109fe32220000","32c51070b","0103103100000000","0b53503000000000"); // EG2
  } else if(trainConfig == 862)  { // EDC 13 TeV EG1
    cuts.AddCutHeavyMesonCalo("0008d113","411792109fe32220000","32c51070b","0103103100000000","0b53503000000000"); // EG1
    //Charged Pion Mass Cut <850 MeV, Neutral Charged Pion Cut 1000< MeV
  } else if(trainConfig == 863)  { // EDC 13 TeV
    cuts.AddCutHeavyMesonCalo("00010113","411792109fe32220000","32c51070c","0103103100000000","0b53503000000000"); // INT7
  } else if(trainConfig == 864)  { // EDC 13 TeV EG2
    cuts.AddCutHeavyMesonCalo("0008e113","411792109fe32220000","32c51070c","0103103100000000","0b53503000000000"); // EG2
  } else if(trainConfig == 865)  { // EDC 13 TeV EG1
    cuts.AddCutHeavyMesonCalo("0008d113","411792109fe32220000","32c51070c","0103103100000000","0b53503000000000"); // EG1
    //Charged Pion Mass Cut <850 MeV, Neutral Charged Pion Cut <850 MeV
  } else if(trainConfig == 866)  { // EDC 13 TeV
    cuts.AddCutHeavyMesonCalo("00010113","411792109fe32220000","32c51070d","0103103100000000","0b53503000000000"); // INT7
  } else if(trainConfig == 867)  { // EDC 13 TeV EG2
    cuts.AddCutHeavyMesonCalo("0008e113","411792109fe32220000","32c51070d","0103103100000000","0b53503000000000"); // EG2
  } else if(trainConfig == 868)  { // EDC 13 TeV EG1
    cuts.AddCutHeavyMesonCalo("0008d113","411792109fe32220000","32c51070d","0103103100000000","0b53503000000000"); // EG1
    //Charged Pion Mass Cut <600 MeV, Neutral Charged Pion Cut <600 MeV
  } else if(trainConfig == 869)  { // EDC 13 TeV
    cuts.AddCutHeavyMesonCalo("00010113","411792109fe32220000","32c51070e","0103103100000000","0b53503000000000"); // INT7
  } else if(trainConfig == 870)  { // EDC 13 TeV EG2
    cuts.AddCutHeavyMesonCalo("0008e113","411792109fe32220000","32c51070e","0103103100000000","0b53503000000000"); // EG2
  } else if(trainConfig == 871)  { // EDC 13 TeV EG1
    cuts.AddCutHeavyMesonCalo("0008d113","411792109fe32220000","32c51070e","0103103100000000","0b53503000000000"); // EG1

  } else {
    Error(Form("GammaConvNeutralMeson_CaloMode_%i",trainConfig), "wrong trainConfig variable no cuts have been specified for the configuration");
    return;
  }


  if(!cuts.AreValid()){
    std::cout << "\n\n****************************************************" << std::endl;
    std::cout << "ERROR: No valid cuts stored in CutHandlerNeutralCalo! Returning..." << std::endl;
    std::cout << "****************************************************\n\n" << std::endl;
    return;
  }

  Int_t numberOfCuts = cuts.GetNCuts();

  TList *EventCutList = new TList();
  TList *ClusterCutList  = new TList();
  TList *NeutralPionCutList = new TList();
  TList *MesonCutList = new TList();
  TList *PionCutList  = new TList();

  TList *HeaderList = new TList();
  TObjString *Header1 = new TObjString("pi0_1");
  HeaderList->Add(Header1);
  TObjString *Header3 = new TObjString("eta_2");
  HeaderList->Add(Header3);

  EventCutList->SetOwner(kTRUE);
  AliConvEventCuts **analysisEventCuts = new AliConvEventCuts*[numberOfCuts];
  ClusterCutList->SetOwner(kTRUE);
  AliCaloPhotonCuts **analysisClusterCuts = new AliCaloPhotonCuts*[numberOfCuts];
  NeutralPionCutList->SetOwner(kTRUE);
  AliConversionMesonCuts **analysisNeutralPionCuts   = new AliConversionMesonCuts*[numberOfCuts];
  MesonCutList->SetOwner(kTRUE);
  AliConversionMesonCuts **analysisMesonCuts   = new AliConversionMesonCuts*[numberOfCuts];
  PionCutList->SetOwner(kTRUE);
  AliPrimaryPionCuts **analysisPionCuts     = new AliPrimaryPionCuts*[numberOfCuts];

  for(Int_t i = 0; i<numberOfCuts; i++){
    //create AliCaloTrackMatcher instance, if there is none present
    TString caloCutPos = cuts.GetClusterCut(i);
    caloCutPos.Resize(1);
    TString TrackMatcherName = Form("CaloTrackMatcher_%s_%i",caloCutPos.Data(),trackMatcherRunningMode);
    if( !(AliCaloTrackMatcher*)mgr->GetTask(TrackMatcherName.Data()) ){
      AliCaloTrackMatcher* fTrackMatcher = new AliCaloTrackMatcher(TrackMatcherName.Data(),caloCutPos.Atoi(),trackMatcherRunningMode);
      fTrackMatcher->SetV0ReaderName(V0ReaderName);
      mgr->AddTask(fTrackMatcher);
      mgr->ConnectInput(fTrackMatcher,0,cinput);
    }

    analysisEventCuts[i] = new AliConvEventCuts();
    analysisEventCuts[i]->SetV0ReaderName(V0ReaderName);
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
    if(runLightOutput>0) analysisEventCuts[i]->SetLightOutput(kTRUE);
    analysisEventCuts[i]->InitializeCutsFromCutString((cuts.GetEventCut(i)).Data());
    if (periodNameV0Reader.CompareTo("") != 0) analysisEventCuts[i]->SetPeriodEnum(periodNameV0Reader);
    EventCutList->Add(analysisEventCuts[i]);
    analysisEventCuts[i]->SetFillCutHistograms("",kFALSE);


    analysisClusterCuts[i] = new AliCaloPhotonCuts();
    analysisClusterCuts[i]->SetV0ReaderName(V0ReaderName);
    if(runLightOutput>0) analysisClusterCuts[i]->SetLightOutput(kTRUE);

    analysisClusterCuts[i]->SetExtendedMatchAndQA(enableExtMatchAndQA);
    analysisClusterCuts[i]->SetCaloTrackMatcherName(TrackMatcherName);
    if( ! analysisClusterCuts[i]->InitializeCutsFromCutString((cuts.GetClusterCut(i)).Data()) ) {
      std::cout<<"ERROR: analysisClusterCuts [" <<i<<"]"<<std::endl;
      return 0;
    } else {
      analysisClusterCuts[i]->InitializeCutsFromCutString((cuts.GetClusterCut(i)).Data());
      ClusterCutList->Add(analysisClusterCuts[i]);
      analysisClusterCuts[i]->SetFillCutHistograms("");
    }

    analysisNeutralPionCuts[i] = new AliConversionMesonCuts();
    analysisNeutralPionCuts[i]->SetUsePtDepSelectionWindow(usePtDepSelectionWindowCut);
    if(runLightOutput>0) analysisNeutralPionCuts[i]->SetLightOutput(kTRUE);
    if( ! analysisNeutralPionCuts[i]->InitializeCutsFromCutString((cuts.GetNDMCut(i)).Data()) ) {
      std::cout<<"ERROR: analysisMesonCuts [ " <<i<<" ] "<<std::endl;
      return 0;
    } else {
      NeutralPionCutList->Add(analysisNeutralPionCuts[i]);
      analysisNeutralPionCuts[i]->SetFillCutHistograms("");
    }

    analysisMesonCuts[i] = new AliConversionMesonCuts();
    if(runLightOutput>0) analysisMesonCuts[i]->SetLightOutput(kTRUE);
    if( ! analysisMesonCuts[i]->InitializeCutsFromCutString((cuts.GetMesonCut(i)).Data()) ) {
      std::cout<<"ERROR: analysisMesonCuts [ " <<i<<" ] "<<std::endl;
      return 0;
    } else {
      MesonCutList->Add(analysisMesonCuts[i]);
      analysisMesonCuts[i]->SetFillCutHistograms("");
    }
    analysisEventCuts[i]->SetAcceptedHeader(HeaderList);

    TString cutName( Form("%s_%s_%s_%s_%s",(cuts.GetEventCut(i)).Data(), (cuts.GetClusterCut(i)).Data(),(cuts.GetPionCut(i)).Data(),(cuts.GetNDMCut(i)).Data(), (cuts.GetMesonCut(i)).Data() ) );
    analysisPionCuts[i] = new AliPrimaryPionCuts();
    analysisPionCuts[i]->SetPrefilterRunFlag(prefilterRunFlag);
    analysisPionCuts[i]->SetPeriodName(periodNameV0Reader);
    if(runLightOutput>0) analysisPionCuts[i]->SetLightOutput(kTRUE);

        if( !analysisPionCuts[i]->InitializeCutsFromCutString((cuts.GetPionCut(i)).Data())) {
      std::cout<< "ERROR:  analysisPionCuts [ " <<i<<" ] "<<std::endl;
      return 0;
    } else {
      PionCutList->Add(analysisPionCuts[i]);
      analysisPionCuts[i]->SetFillCutHistograms("",kFALSE,cutName);
    }
  }

  task->SetNDMRecoMode(neutralPionMode);
  task->SetEventCutList(numberOfCuts,EventCutList);
  task->SetClusterCutList(ClusterCutList);
  task->SetNeutralPionCutList(NeutralPionCutList);
  task->SetMesonCutList(MesonCutList);
  task->SetPionCutList(PionCutList);

  task->SetMoveParticleAccordingToVertex(kFALSE);
  task->SetSelectedHeavyNeutralMeson(selectHeavyNeutralMeson);

  task->SetDoMesonQA(enableQAMesonTask );

  //connect containers
  AliAnalysisDataContainer *coutput =
  mgr->CreateContainer(Form("GammaConvNeutralMesonPiPlPiMiNeutralMeson_%i_%i_%i.root",selectHeavyNeutralMeson,neutralPionMode, trainConfig), TList::Class(),
              AliAnalysisManager::kOutputContainer,Form("GammaConvNeutralMesonPiPlPiMiNeutralMeson_%i_%i_%i.root",selectHeavyNeutralMeson,neutralPionMode, trainConfig));

  mgr->AddTask(task);
  mgr->ConnectInput(task,0,cinput);
  mgr->ConnectOutput(task,1,coutput);

  return;

}
