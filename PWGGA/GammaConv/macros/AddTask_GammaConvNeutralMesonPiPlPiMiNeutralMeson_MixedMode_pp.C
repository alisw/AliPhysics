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
void AddTask_GammaConvNeutralMesonPiPlPiMiNeutralMeson_MixedMode_pp(
    Int_t     trainConfig                   = 1,
    Int_t     isMC                          = 0,                        //run MC
    TString   photonCutNumberV0Reader       = "",                       // 00000008400000000100000000 nom. B, 00000088400000000100000000 low B
    Int_t     selectHeavyNeutralMeson       = 0,                        //run eta prime instead of omega
    Int_t     enableQAMesonTask             = 1,                        //enable QA in AliAnalysisTaskNeutralMesonToPiPlPiMiNeutralMeson
    Int_t     enableExtMatchAndQA           = 0,                        // disabled (0), extMatch (1), extQA_noCellQA (2), extMatch+extQA_noCellQA (3), extQA+cellQA (4), extMatch+extQA+cellQA (5)
    Int_t     enableTriggerMimicking        = 0,                        // enable trigger mimicking
    Bool_t    enableTriggerOverlapRej       = kFALSE,                   // enable trigger overlap rejection
    TString   settingMaxFacPtHard           = "3.",                     // maximum factor between hardest jet and ptHard generated
    TString   fileNameExternalInputs        = "MCSpectraInput.root",    // path to file for weigting input
    Bool_t    doWeighting                   = kFALSE,                   //enable Weighting
    Bool_t    enableElecDeDxPostCalibration = kFALSE,                   // enable post calibration of elec pos dEdX
    TString   generatorName               = "HIJING",
    Double_t  tolerance                   = -1,
    TString   periodNameV0Reader          = "",                       // period Name for V0Reader
    Int_t     runLightOutput              = 0,                        // run light output option 0: no light output 1: most cut histos stiched off 2: unecessary omega hists turned off as well
    Int_t     prefilterRunFlag            = 1500,                     // flag to change the prefiltering of ESD tracks. See SetHybridTrackCutsAODFiltering() in AliPrimaryPionCuts
    Bool_t    usePtDepSelectionWindowCut  = kFALSE,                   // use pt dependent meson selection window cut
    TString   additionalTrainConfig       = "0"                       // additional counter for trainconfig, this has to be always the last parameter
  ) {

  AliCutHandlerPCM cuts(13);
  TString fileNamedEdxPostCalib = cuts.GetSpecialFileNameFromString (fileNameExternalInputs, "FEPC:");
  TString fileNameCustomTriggerMimicOADB   = cuts.GetSpecialFileNameFromString (fileNameExternalInputs, "FTRM:");

  Int_t trackMatcherRunningMode = 0; // CaloTrackMatcher running mode
  //parse additionalTrainConfig flag
  TObjArray *rAddConfigArr = additionalTrainConfig.Tokenize("_");
  if(rAddConfigArr->GetEntries()<1){cout << "ERROR during parsing of additionalTrainConfig String '" << additionalTrainConfig.Data() << "'" << endl; return;}
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
        cout << Form("INFO: AddTask_GammaConvNeutralMesonPiPlPiMiPiZero_MixedMode_pp will use running mode '%i' for the TrackMatcher!",trackMatcherRunningMode) << endl;
      }
    }
  }
  TString sAdditionalTrainConfig = rAdditionalTrainConfig->GetString();
  if (sAdditionalTrainConfig.Atoi() > 0){
    trainConfig = trainConfig + sAdditionalTrainConfig.Atoi();
    cout << "INFO: AddTask_GammaConvNeutralMesonPiPlPiMiNeutralMeson_MixedMode_pp running additionalTrainConfig '" << sAdditionalTrainConfig.Atoi() << "', train config: '" << trainConfig << "'" << endl;
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
  Int_t neutralPionMode = 1;

  // ================== GetAnalysisManager ===============================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error(Form("AddTask_GammaConvNeutralMesonPiPlPiMiNeutralMeson_MixedMode_pp_%i",trainConfig), "No analysis manager found.");
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
  //========= Add Pion Selector ====================
  TString PionCuts          = "000000200";
  if( !(AliPrimaryPionSelector*)mgr->GetTask("PionSelector") ){
    AliPrimaryPionSelector *fPionSelector = new AliPrimaryPionSelector("PionSelector");
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
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009327000008250400000","1111100047032230000","000010400","0103503a00000000","0103503000000000");
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009327000008250400000","1111100047032230000","000010400","0103503a00000000","0103503000000000");
  } else if( trainConfig == 2 ) {
    // closing charged pion cuts, minimum TPC cluster = 80, TPC dEdx pi = \pm 3 sigma, min pt charged pi = 100 MeV
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009327000008250400000","1111100047032230000","002010700","0103503a00000000","0103503000000000");
  } else if( trainConfig == 3) {
    // eta < 0.9
    // closing charged pion cuts, minimum TPC cluster = 80, TPC dEdx pi = \pm 3 sigma, pi+pi- mass cut of 0.65, min pt charged pi = 100 MeV
    // closing neural pion cuts, 0.1 < M_gamma,gamma < 0.15
    // maxChi2 per cluster TPC <4, require TPC refit, DCA XY pT dependend 0.0182+0.0350/pt^1.01, DCA_Z = 3.0
    // timing cluster cut open
    cuts.AddCutHeavyMesonPCMCalo("00010113","00200009327000008250400000","1111100047032230000","30a330708","0103503400000000","0153503000000000"); // all of the above
    cuts.AddCutHeavyMesonPCMCalo("00052113","00200009327000008250400000","1111100047032230000","30a330708","0103503400000000","0153503000000000"); // all of the above
    cuts.AddCutHeavyMesonPCMCalo("00062113","00200009327000008250400000","1111100047032230000","30a330708","0103503400000000","0153503000000000"); // all of the above
    cuts.AddCutHeavyMesonPCMCalo("00083113","00200009327000008250400000","1111100047032230000","30a330708","0103503400000000","0153503000000000"); // all of the above
    cuts.AddCutHeavyMesonPCMCalo("00085113","00200009327000008250400000","1111100047032230000","30a330708","0103503400000000","0153503000000000"); // all of the above
  } else if( trainConfig == 4) {
    // same as 3 but only MB
    cuts.AddCutHeavyMesonPCMCalo("00010113","00200009327000008250400000","1111100047032230000","30a330708","0103503400000000","0153503000000000"); // all of the above
  } else if( trainConfig == 10)  { // EMCal 7 TeV
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","1111111047032230000","302010708","0103603800000000","0153503000000000");
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","1111111047032230000","322010708","0103603800000000","0153503000000000"); // with ITS requirement
  } else if( trainConfig == 50 ) { // PHOS 7 TeV
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","2444411043012300000","302010708","0103603600000000","0153503000000000"); //PCM-PHOS non lin
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","2444411043012300000","322010708","0103603600000000","0153503000000000"); //with ITS requirement

    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //                                          OMEGA MESON
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  } else if( trainConfig == 100 ) {
    // everything open, min pt charged pi = 100 MeV
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009327000008250400000","1111100047032230000","000010400","0103503a00000000","0103503000000000");
  } else if( trainConfig == 101 ) {
    // closing charged pion cuts, minimum TPC cluster = 80, TPC dEdx pi = \pm 3 sigma, min pt charged pi = 100 MeV
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009327000008250400000","1111100047032230000","002010700","0103503a00000000","0103503000000000");
  } else if( trainConfig == 102) {
    // eta < 0.9
    // closing charged pion cuts, minimum TPC cluster = 80, TPC dEdx pi = \pm 3 sigma, pi+pi- mass cut of 0.65, min pt charged pi = 100 MeV
    // closing neural pion cuts, 0.1 < M_gamma,gamma < 0.15
    // maxChi2 per cluster TPC <4, require TPC refit, DCA XY pT dependend 0.0182+0.0350/pt^1.01, DCA_Z = 3.0
    // timing cluster cut open
    cuts.AddCutHeavyMesonPCMCalo("00010113","00200009327000008250400000","1111100047032230000","30a330708","0103503400000000","0153503000000000"); // all of the above
    cuts.AddCutHeavyMesonPCMCalo("00052113","00200009327000008250400000","1111100047032230000","30a330708","0103503400000000","0153503000000000"); // all of the above
    cuts.AddCutHeavyMesonPCMCalo("00062113","00200009327000008250400000","1111100047032230000","30a330708","0103503400000000","0153503000000000"); // all of the above
    cuts.AddCutHeavyMesonPCMCalo("00083113","00200009327000008250400000","1111100047032230000","30a330708","0103503400000000","0153503000000000"); // all of the above
    cuts.AddCutHeavyMesonPCMCalo("00085113","00200009327000008250400000","1111100047032230000","30a330708","0103503400000000","0153503000000000"); // all of the above

  } else if( trainConfig == 103) {
    // same as 102 but only MB
    cuts.AddCutHeavyMesonPCMCalo("00010113","00200009327000008250400000","1111100047032230000","30a330708","0103503400000000","0153503000000000"); // all of the above
  } else if( trainConfig == 110)  { // EMCal 7 TeV
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","1111111047032230000","302010708","0103603700000000","0153503000000000");

    // pp 5 TeV test
  } else if ( trainConfig == 111) { // with TPC refit + ITS requirement
    cuts.AddCutHeavyMesonPCMCalo("00010113","00200009227000008250400000","1111111047032230000","32c010708","0103603700000000","0153503000000000"); // INT7
    cuts.AddCutHeavyMesonPCMCalo("00052113","00200009227000008250400000","1111111047032230000","32c010708","0103603700000000","0153503000000000"); // EMC7
    cuts.AddCutHeavyMesonPCMCalo("00085113","00200009227000008250400000","1111111047032230000","32c010708","0103603700000000","0153503000000000"); // EG2
    cuts.AddCutHeavyMesonPCMCalo("00083113","00200009227000008250400000","1111111047032230000","32c010708","0103603700000000","0153503000000000"); // EG1

    // pp 13 TeV test
  } else if ( trainConfig == 112) { // with TPC refit + ITS requirement
    cuts.AddCutHeavyMesonPCMCalo("00010113","00200009227000008250400000","1111111047032230000","32c51070a","0103603700000000","0153503000000000"); // INT7
    cuts.AddCutHeavyMesonPCMCalo("00085113","00200009227000008250400000","1111111047032230000","32c51070a","0103603700000000","0153503000000000"); // EG2
    cuts.AddCutHeavyMesonPCMCalo("00083113","00200009227000008250400000","1111111047032230000","32c51070a","0103603700000000","0153503000000000"); // EG1
  } else if ( trainConfig == 113) { // pp13 TeV AOD and ESD Comparison
    cuts.AddCutHeavyMesonPCMCalo("00010113","00200009227000008250400000","1111111047032230000","32c510708","0103603700000000","0153503000000000"); // INT7
    // EMCal LHC11pp 7 TeV
  } else if ( trainConfig == 115) { // EMCal LHC11 std (no background)
    cuts.AddCutHeavyMesonPCMCalo("00010113","00200009227000008250400000","111111105f032230000","32c51070a","0103603700000000","0453503000000000"); // INT7 (LHC10 acc)
    cuts.AddCutHeavyMesonPCMCalo("00052113","00200009227000008250400000","111111105f032230000","32c51070a","0103603700000000","0453503000000000"); // EMC7 (LHC10 acc)
    cuts.AddCutHeavyMesonPCMCalo("00052113","00200009227000008250400000","111111105f032230000","32c51070a","0103653700000000","0453503000000000"); // 5 geV min pT pi0
  } else if ( trainConfig == 116) { // EMCal LHC11 with LHC10 acc cuts
    cuts.AddCutHeavyMesonPCMCalo("00010113","00200009227000008250400000","1111a1105f032230000","32c51070a","0103603700000000","0453503000000000"); // INT7 (LHC10 acc)
    cuts.AddCutHeavyMesonPCMCalo("00052113","00200009227000008250400000","1111a1105f032230000","32c51070a","0103603700000000","0453503000000000"); // EMC7 (LHC10 acc)
   // Test for EMCal (13 TeV) without background calculation
  } else if ( trainConfig == 117) { // with TPC refit + ITS requirement
    cuts.AddCutHeavyMesonPCMCalo("00010113","00200009227000008250400000","111111104f032230000","32c51070a","0103603700000000","0453503000000000"); // INT7
    cuts.AddCutHeavyMesonPCMCalo("00085113","00200009227000008250400000","111111104f032230000","32c51070a","0103603700000000","0453503000000000"); // EG2
    cuts.AddCutHeavyMesonPCMCalo("00083113","00200009227000008250400000","111111104f032230000","32c51070a","0103603700000000","0453503000000000"); // EG1
    // ---------------------------------
    // systematic studies 7 TeV (EMCal)
    // ---------------------------------

  } else if( trainConfig == 120)  { //std
    cuts.AddCutHeavyMesonPCMCalo("00000113","0dm0000922700000dge0404000","1111a3104f032230000","32c51070a","0103603700000000","0153503000000000"); // with TPC refit (new standard)
    cuts.AddCutHeavyMesonPCMCalo("00000113","0dm0000922700000dge0404000","1111a3104f032230000","32c51070a","0103603f00000000","0153503000000000"); // with TPC refit (new standard)
  } else if( trainConfig == 121)  { //pileup
    cuts.AddCutHeavyMesonPCMCalo("00000113","0dm0000922700000dge0404000","1111a3104f032230000","32c51070a","0103603700000000","0153503000000000"); // pileup
  } else if ( trainConfig == 122){ // singlept
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200019227000008250400000","1111a3104f032230000","32c51070a","0103603700000000","0153503000000000"); // 0.100 GeV
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200049227000008250400000","1111a3104f032230000","32c51070a","0103603700000000","0153503000000000"); // 0.075 GeV
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200069227000008250400000","1111a3104f032230000","32c51070a","0103603700000000","0153503000000000"); // 0.04 GeV
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200059227000008250400000","1111a3104f032230000","32c51070a","0103603700000000","0153503000000000"); // 0.125 GeV
  } else if ( trainConfig == 123){ // clstpc
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200008227000008250400000","1111a3104f032230000","32c51070a","0103603700000000","0153503000000000"); // fMinClsTPCToF= 0.35;fUseCorrectedTPCClsInfo=1;
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200006227000008250400000","1111a3104f032230000","32c51070a","0103603700000000","0153503000000000"); // fMinClsTPCToF= 0.70;fUseCorrectedTPCClsInfo=1;
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200001227000008250400000","1111a3104f032230000","32c51070a","0103603700000000","0153503000000000"); // fMinClsTPCToF= 60;fUseCorrectedTPCClsInfo=0;
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200002227000008250400000","1111a3104f032230000","32c51070a","0103603700000000","0153503000000000"); // fMinClsTPCToF= 80;fUseCorrectedTPCClsInfo=0;
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200003227000008250400000","1111a3104f032230000","32c51070a","0103603700000000","0153503000000000"); // fMinClsTPCToF= 100;fUseCorrectedTPCClsInfo=0;
  } else if ( trainConfig == 124){ // TPCdEdxCutElectron
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009327000008250400000","1111a3104f032230000","32c51070a","0103603700000000","0153503000000000"); // -4,5
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009627000008250400000","1111a3104f032230000","32c51070a","0103603700000000","0153503000000000"); // -2.5,4
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009427000008250400000","1111a3104f032230000","32c51070a","0103603700000000","0153503000000000"); // -6,7
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009527000008250400000","1111a3104f032230000","32c51070a","0103603700000000","0153503000000000"); // -4,4
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009627000008250400000","1111a3104f032230000","32c51070a","0103603700000000","0153503000000000"); // -2.5,4
  } else if ( trainConfig == 125){ // TPCdEdxCutPion
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009217000008250400000","1111a3104f032230000","32c51070a","0103603700000000","0153503000000000"); // fPIDnSigmaAbovePionLine=0; fPIDnSigmaAbovePionLineHighPt=-10;
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009237000008250400000","1111a3104f032230000","32c51070a","0103603700000000","0153503000000000"); // fPIDnSigmaAbovePionLine=2.5; fPIDnSigmaAbovePionLineHighPt=-10;
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009247000008250400000","1111a3104f032230000","32c51070a","0103603700000000","0153503000000000"); // fPIDnSigmaAbovePionLine=0.5; fPIDnSigmaAbovePionLineHighPt=-10;
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009257000008250400000","1111a3104f032230000","32c51070a","0103603700000000","0153503000000000"); // fPIDnSigmaAbovePionLine=2; fPIDnSigmaAbovePionLineHighPt=-10;
  } else if ( trainConfig == 126){ // QtMaxCut
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000003250400000","1111a3104f032230000","32c51070a","0103603700000000","0153503000000000");
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000009250400000","1111a3104f032230000","32c51070a","0103603700000000","0153503000000000");
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000002250400000","1111a3104f032230000","32c51070a","0103603700000000","0153503000000000");
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000006250400000","1111a3104f032230000","32c51070a","0103603700000000","0153503000000000");
  } else if ( trainConfig == 127){ // Chi2GammaCut
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008150400000","1111a3104f032230000","32c51070a","0103603700000000","0153503000000000");
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008850400000","1111a3104f032230000","32c51070a","0103603700000000","0153503000000000");
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008a50400000","1111a3104f032230000","32c51070a","0103603700000000","0153503000000000");
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008950400000","1111a3104f032230000","32c51070a","0103603700000000","0153503000000000");
  } else if ( trainConfig == 128){ // PsiPair
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008260400000","1111a3104f032230000","32c51070a","0103603700000000","0153503000000000");
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008280400000","1111a3104f032230000","32c51070a","0103603700000000","0153503000000000");
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008210400000","1111a3104f032230000","32c51070a","0103603700000000","0153503000000000");
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008220400000","1111a3104f032230000","32c51070a","0103603700000000","0153503000000000");

  } else if ( trainConfig == 129){ // NonLin
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","1111a12047032230000","32c51070a","0103603700000000","0153503000000000");
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","1111a13047032230000","32c51070a","0103603700000000","0153503000000000");
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","1111a21047032230000","32c51070a","0103603700000000","0153503000000000");
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","1111a22047032230000","32c51070a","0103603700000000","0153503000000000");
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","1111a23047032230000","32c51070a","0103603700000000","0153503000000000");
  } else if ( trainConfig == 130){ // Timing diff
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","1111a110b7032230000","32c51070a","0103603700000000","0153503000000000");
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","1111a110c7032230000","32c51070a","0103603700000000","0153503000000000");
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","1111a110d7032230000","32c51070a","0103603700000000","0153503000000000");
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","1111a110e7032230000","32c51070a","0103603700000000","0153503000000000");
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","1111a110f7032230000","32c51070a","0103603700000000","0153503000000000");
  } else if ( trainConfig == 131){ // Track Matching
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","1111a11046032230000","32c51070a","0103603700000000","0153503000000000");
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","1111a11048032230000","32c51070a","0103603700000000","0153503000000000");
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","1111a11049032230000","32c51070a","0103603700000000","0153503000000000");
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","1111a1104a032230000","32c51070a","0103603700000000","0153503000000000");
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","1111a1104b032230000","32c51070a","0103603700000000","0153503000000000");
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","1111a11043032230000","32c51070a","0103603700000000","0153503000000000");
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","1111a1104f032230000","32c51070a","0103603700000000","0153503000000000"); // E/p
  } else if ( trainConfig == 132){ // MinEnergy (of cluster)
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","1111a11047022230000","32c51070a","0103603700000000","0153503000000000");
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","1111a11047042230000","32c51070a","0103603700000000","0153503000000000");
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","1111a11047052230000","32c51070a","0103603700000000","0153503000000000");
  } else if ( trainConfig == 133){ // MinNCells
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","1111a11047032130000","32c51070a","0103603700000000","0153503000000000");
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","1111a11047032330000","32c51070a","0103603700000000","0153503000000000");
  } else if ( trainConfig == 134){ // MinMaxM02
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","1111a11047032330000","32c51070a","0103603700000000","0153503000000000");
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","1111a11047032130000","32c51070a","0103603700000000","0153503000000000");
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","1111a11047032240000","32c51070a","0103603700000000","0153503000000000");
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","1111a11047032220000","32c51070a","0103603700000000","0153503000000000");

// other variations
  } else if( trainConfig == 140)  {  // pT Cut
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","1111a11047032230000","32c50070a","0103603700000000","0153503000000000"); // pt>0.075
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","1111a11047032230000","32c52070a","0103603700000000","0153503000000000"); // pt>0.125
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","1111a11047032230000","32c53070a","0103603700000000","0153503000000000"); // pt>0.15
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","1111a11047032230000","32c54070a","0103603700000000","0153503000000000"); // pt>0.4
  } else if (trainConfig == 141) { // TPCdEdxCutPion
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","1111a11047032230000","32c51050a","0103603700000000","0153503000000000"); // -4,4
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","1111a11047032230000","32c51080a","0103603700000000","0153503000000000"); // -2,3
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","1111a11047032230000","32c51020a","0103603700000000","0153503000000000"); // -6,7
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","1111a11047032230000","32c51030a","0103603700000000","0153503000000000"); // -5,5
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","1111a11047032230000","32c51040a","0103603700000000","0153503000000000"); // -4,5
  } else if (trainConfig == 142) { // Mass window
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","1111a11047032230000","32c51070a","0103603100000000","0153503000000000"); // 100-145
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","1111a11047032230000","32c51070a","0103603200000000","0153503000000000"); // 110-145
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","1111a11047032230000","32c51070a","0103603300000000","0153503000000000"); // 120-145
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","1111a11047032230000","32c51070a","0103603400000000","0153503000000000"); // 100-150
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","1111a11047032230000","32c51070a","0103603500000000","0153503000000000"); // 110-150
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","1111a11047032230000","32c51070a","0103603600000000","0153503000000000"); // 120-150
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","1111a11047032230000","32c51070a","0103603a00000000","0153503000000000"); //  80-145
  } else if( trainConfig == 145)  { // no event mixing (used for trigger LHC11)
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","1111111047032230000","32b11070a","0103603700000000","0453503000000000"); // without background
    // ---------------------------------
    // systematic studies 7 TeV (PCM-PHOS)
    // ---------------------------------

  } else if( trainConfig == 150 ) { // PHOS 7 TeV
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","2444411044013300000","32c51070a","0103603n00000000","0153503000000000"); // with TPC refit + ITS requirement

    // *************Variations in AliConvEventCuts**************************
  } else if (trainConfig == 151) { // remove pileup
    cuts.AddCutHeavyMesonPCMCalo("00000013","00200009227000008250400000","2444411044013300000","32c51070a","0103603n00000000","0153503000000000"); // off


    // *************Variations in AliCaloPhotonsCut**************************
  } else if(trainConfig == 152)  { // Timing diff(std is -100ns to 100ns)
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","24444110b4013300000","32c51070a","0103603n00000000","0153503000000000"); // 130ns
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","24444110c4013300000","32c51070a","0103603n00000000","0153503000000000"); // 110ns
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","24444110d4013300000","32c51070a","0103603n00000000","0153503000000000"); // 120ns
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","24444110e4013300000","32c51070a","0103603n00000000","0153503000000000"); // 90ns
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","24444110f4013300000","32c51070a","0103603n00000000","0153503000000000"); // 80ns
  } else if(trainConfig == 153)  { // Track matching
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","2444411041013300000","32c51070a","0103603n00000000","0153503000000000"); //
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","2444411043013300000","32c51070a","0103603n00000000","0153503000000000"); //
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","2444411045013300000","32c51070a","0103603n00000000","0153503000000000"); //
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","2444411046013300000","32c51070a","0103603n00000000","0153503000000000"); //
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","2444411047013300000","32c51070a","0103603n00000000","0153503000000000"); //
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","2444411048013300000","32c51070a","0103603n00000000","0153503000000000"); //
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","2444411049013300000","32c51070a","0103603n00000000","0153503000000000"); //
  } else if(trainConfig == 154)  { // MinEnergy (of cluster) (std is 0.5 GeV)
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","2444411044023300000","32c51070a","0103603n00000000","0153503000000000"); // 0.6
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","2444411044033300000","32c51070a","0103603n00000000","0153503000000000"); // 0.1
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","2444411044043300000","32c51070a","0103603n00000000","0153503000000000"); // 0.7
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","2444411044083300000","32c51070a","0103603n00000000","0153503000000000"); // 0.8
  } else if(trainConfig == 155)  { // Min N of cells (std is 2)
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","2444411044011300000","32c51070a","0103603n00000000","0153503000000000"); // 1
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","2444411044012300000","32c51070a","0103603n00000000","0153503000000000"); // 3
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","2444411044014300000","32c51070a","0103603n00000000","0153503000000000"); // 4
  } else if(trainConfig == 156)  { // MinMaxM02 (std is >0.2)
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","2444411044013200000","32c51070a","0103603n00000000","0153503000000000"); // >0.1
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","2444411044013100000","32c51070a","0103603n00000000","0153503000000000"); // >0.002

    // *************Variations in AliPrimaryPionCuts******************
  } else if( trainConfig == 157)  { // ClsTPCCut
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","2444411044013300000","32151070a","0103603n00000000","0153503000000000"); // 70
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","2444411044013300000","32351070a","0103603n00000000","0153503000000000"); // 100
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","2444411044013300000","32551070a","0103603n00000000","0153503000000000"); // 35% of findable clusters
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","2444411044013300000","32651070a","0103603n00000000","0153503000000000"); // 60% of findable clusters
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","2444411044013300000","32b51070a","0103603n00000000","0153503000000000"); // PHOS public note

  } else if ( trainConfig == 158) { // DCACut
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","2444411044013300000","32c51070a","0103603n00000000","0153503000000000"); // XYPtDep("0.0182+0.0350/pt^1.01");
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","2444411044013300000","32c51070a","0103603n00000000","0153503000000000"); // z=2cm xy=1cm
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","2444411044013300000","32c51070a","0103603n00000000","0153503000000000"); // z=3cm XYPtDep("0.0182+0.0350/pt^1.01");
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","2444411044013300000","32c51070a","0103603n00000000","0153503000000000"); // z=3cm xy=0.5
  } else if ( trainConfig == 159) { // pT cut
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","2444411044013300000","32c50070a","0103603n00000000","0153503000000000"); // pt>0.075
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","2444411044013300000","32c52070a","0103603n00000000","0153503000000000"); // pt>0.125
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","2444411044013300000","32c53070a","0103603n00000000","0153503000000000"); // pt>0.15
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","2444411044013300000","32c54070a","0103603n00000000","0153503000000000"); // pt>0.4
  } else if ( trainConfig == 160) { // TPDdEdxCutPion
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","2444411044013300000","32c51050a","0103603n00000000","0153503000000000"); // -4,4
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","2444411044013300000","32c51080a","0103603n00000000","0153503000000000"); // -2,3
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","2444411044013300000","32c51020a","0103603n00000000","0153503000000000"); // -6,7
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","2444411044013300000","32c51030a","0103603n00000000","0153503000000000"); // -5,5
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","2444411044013300000","32c51040a","0103603n00000000","0153503000000000"); // -4,5
  } else if ( trainConfig == 161) { // Mass cut
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","2444411044013300000","32c510707","0103603n00000000","0153503000000000"); // 0.7 GeV
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","2444411044013300000","32c510706","0103603n00000000","0153503000000000"); // 0.65 GeV
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","2444411044013300000","32c510701","0103603n00000000","0153503000000000"); // 1 GeV
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","2444411044013300000","32c510702","0103603n00000000","0153503000000000"); // 0.75 GeV
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","2444411044013300000","32c510704","0103603n00000000","0153503000000000"); // 0.54 eta mass
    // *************Variations in AliConversionMesonCuts (NeutralPion) ******************
  } else if ( trainConfig == 162) { // RapidityMesonCut
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","2444411044013300000","32c51070a","0103503n00000000","0153503000000000"); // 0.85
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","2444411044013300000","32c51070a","0103303n00000000","0153503000000000"); // 0.6
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","2444411044013300000","32c51070a","0103203n00000000","0153503000000000"); // 0.7
  } else if ( trainConfig == 163) { // pT
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","2444411044013300000","32c51070a","0103613n00000000","0153503000000000"); // 0.4
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","2444411044013300000","32c51070a","0103623n00000000","0153503000000000"); // 0.7
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","2444411044013300000","32c51070a","0103673n00000000","0153503000000000"); // 0.5
  } else if ( trainConfig == 164) { // max opening angle
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","2444411044013300000","32c51070a","0103607n00000001","0153503000000000"); //
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","2444411044013300000","32c51070a","0103605n00000002","0153503000000000"); //
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","2444411044013300000","32c51070a","0103600n00000003","0153503000000000"); //
  } else if ( trainConfig == 165) { // selectionWindow (std is 120-160)
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","2444411044013300000","32c51070a","0103603100000000","0153503000000000"); // 0.1-0.145
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","2444411044013300000","32c51070a","0103603200000000","0153503000000000"); // 0.11-0.145
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","2444411044013300000","32c51070a","0103603300000000","0153503000000000"); // 0.12-0.145
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","2444411044013300000","32c51070a","0103603600000000","0153503000000000"); // 0.12-0.5
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","2444411044013300000","32c51070a","0103603700000000","0153503000000000"); // 0.1 -0.155
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","2444411044013300000","32c51070a","0103603900000000","0153503000000000"); // 0.11 -0.155
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","2444411044013300000","32c51070a","0103603a00000000","0153503000000000"); // 0.08 -0.145
    // *************Variations in AliConversionMesonCuts (omega) ******************
  } else if ( trainConfig == 166) { // selectionWindow (std is 120-160)
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","2444411044013300000","32c51070a","0103603n00000000","0a53503000000000"); // likesign
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","2444411044013300000","32c51070a","0103603n00000000","0b53503000000000"); // sideband right
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","2444411044013300000","32c51070a","0103603n00000000","0c53503000000000"); // sideband left
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","2444411044013300000","32c51070a","0103603n00000000","0d53503000000000"); // sideband both sides
  } else if ( trainConfig == 167) { // Number of BckEvents
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","2444411044013300000","32c51070a","0103603n00000000","0133503000000000"); // 20
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","2444411044013300000","32c51070a","0103603n00000000","0163503000000000"); // 80
  } else if ( trainConfig == 168) { // rapidity cut
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","2444411044013300000","32c51070a","0103603n00000000","0153203000000000"); // 0.7
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","2444411044013300000","32c51070a","0103603n00000000","0153003000000000"); // 1.35
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","2444411044013300000","32c51070a","0103603n00000000","0153303000000000"); // 0.6
  } else if ( trainConfig == 169) { // alphacut
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","2444411044013300000","32c51070a","0103603n00000000","0153507000000000"); // 0-0.85
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","2444411044013300000","32c51070a","0103603n00000000","0153505000000000"); // 0-0.75
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","2444411044013300000","32c51070a","0103603n00000000","0153500000000000"); // 0-0.7

    // *************Variations in AliConversionPhotonCuts*******************
  } else if ( trainConfig == 170) { // EtaCut
    cuts.AddCutHeavyMesonPCMCalo("00000113","01200009227000008250400000","2444411044013300000","32c51070a","0103603n00000000","0153503000000000"); // eta 0.6
    cuts.AddCutHeavyMesonPCMCalo("00000113","04200009227000008250400000","2444411044013300000","32c51070a","0103603n00000000","0153503000000000"); // eta 0.75
    cuts.AddCutHeavyMesonPCMCalo("00000113","0d200009227000008250400000","2444411044013300000","32c51070a","0103603n00000000","0153503000000000"); // eta 0.8
  } else if ( trainConfig == 171) { // SinglePtCut
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200019227000008250400000","2444411044013300000","32c51070a","0103603n00000000","0153503000000000"); // 0.100 GeV
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200049227000008250400000","2444411044013300000","32c51070a","0103603n00000000","0153503000000000"); // 0.075 GeV
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200069227000008250400000","2444411044013300000","32c51070a","0103603n00000000","0153503000000000"); // 0.04 GeV
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200059227000008250400000","2444411044013300000","32c51070a","0103603n00000000","0153503000000000"); // 0.125 GeV
  } else if ( trainConfig == 172) { // ClsTPCCut
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200008227000008250400000","2444411044013300000","32c51070a","0103603n00000000","0153503000000000"); // fMinClsTPCToF= 0.35;fUseCorrectedTPCClsInfo=1;
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200006227000008250400000","2444411044013300000","32c51070a","0103603n00000000","0153503000000000"); // fMinClsTPCToF= 0.70;fUseCorrectedTPCClsInfo=1;
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200001227000008250400000","2444411044013300000","32c51070a","0103603n00000000","0153503000000000"); // fMinClsTPCToF= 60;fUseCorrectedTPCClsInfo=0;
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200002227000008250400000","2444411044013300000","32c51070a","0103603n00000000","0153503000000000"); // fMinClsTPCToF= 80;fUseCorrectedTPCClsInfo=0;
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200003227000008250400000","2444411044013300000","32c51070a","0103603n00000000","0153503000000000"); // fMinClsTPCToF= 100;fUseCorrectedTPCClsInfo=0;
  } else if ( trainConfig == 173) { // TPCdEdxCutElectron
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009327000008250400000","2444411044013300000","32c51070a","0103603n00000000","0153503000000000"); // -4,5
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009627000008250400000","2444411044013300000","32c51070a","0103603n00000000","0153503000000000"); // -2.5,4
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009427000008250400000","2444411044013300000","32c51070a","0103603n00000000","0153503000000000"); // -6,7
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009527000008250400000","2444411044013300000","32c51070a","0103603n00000000","0153503000000000"); // -4,4
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009627000008250400000","2444411044013300000","32c51070a","0103603n00000000","0153503000000000"); // -2.5,4
  } else if ( trainConfig == 174) { // TPCdEdxCutPion
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009217000008250400000","2444411044013300000","32c51070a","0103603n00000000","0153503000000000"); // fPIDnSigmaAbovePionLine=0; fPIDnSigmaAbovePionLineHighPt=-10;
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009237000008250400000","2444411044013300000","32c51070a","0103603n00000000","0153503000000000"); // fPIDnSigmaAbovePionLine=2.5; fPIDnSigmaAbovePionLineHighPt=-10;
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009247000008250400000","2444411044013300000","32c51070a","0103603n00000000","0153503000000000"); // fPIDnSigmaAbovePionLine=0.5; fPIDnSigmaAbovePionLineHighPt=-10;
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009257000008250400000","2444411044013300000","32c51070a","0103603n00000000","0153503000000000"); // fPIDnSigmaAbovePionLine=2; fPIDnSigmaAbovePionLineHighPt=-10;
  } else if ( trainConfig == 175) { // piMomdedxSigmaCut
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009225000008250400000","2444411044013300000","32c51070a","0103603n00000000","0153503000000000"); // 0.3 GeV
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009220000008250400000","2444411044013300000","32c51070a","0103603n00000000","0153503000000000"); // 0.5 GeV
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009228000008250400000","2444411044013300000","32c51070a","0103603n00000000","0153503000000000"); // 0.2 GeV
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009226000008250400000","2444411044013300000","32c51070a","0103603n00000000","0153503000000000"); // 0.25 GeV
  } else if ( trainConfig == 176) { // QtMaxCut
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000003250400000","2444411044013300000","32c51070a","0103603n00000000","0153503000000000"); // fQtMax=0.05; fDo2DQt=kFALSE;
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000009250400000","2444411044013300000","32c51070a","0103603n00000000","0153503000000000"); // fQtMax=0.03; fDo2DQt=kTRUE;
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000002250400000","2444411044013300000","32c51070a","0103603n00000000","0153503000000000"); // fQtMax=0.06; fDo2DQt=kTRUE;
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000006250400000","2444411044013300000","32c51070a","0103603n00000000","0153503000000000"); // fQtMax=0.02; fDo2DQt=kTRUE;
  } else if ( trainConfig == 177) { // Chi2GammaCut
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008150400000","2444411044013300000","32c51070a","0103603n00000000","0153503000000000"); // fChi2CutConversion = 50.;
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008850400000","2444411044013300000","32c51070a","0103603n00000000","0153503000000000"); // fChi2CutConversion = 20.;
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008a50400000","2444411044013300000","32c51070a","0103603n00000000","0153503000000000"); // fChi2CutConversion = 25.;
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008950400000","2444411044013300000","32c51070a","0103603n00000000","0153503000000000"); // fChi2CutConversion = 15.;
  } else if ( trainConfig == 178) { // PsiPair
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008260400000","2444411044013300000","32c51070a","0103603n00000000","0153503000000000"); // 0.05 and 2D
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008280400000","2444411044013300000","32c51070a","0103603n00000000","0153503000000000"); // 0.2 and 2D
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008210400000","2444411044013300000","32c51070a","0103603n00000000","0153503000000000"); // 0.1 and 1D
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008220400000","2444411044013300000","32c51070a","0103603n00000000","0153503000000000"); // 0.05 and 1D
  } else if ( trainConfig == 179) { // CosinePointingAngle
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250700000","2444411044013300000","32c51070a","0103603n00000000","0153503000000000"); // fCosPAngleCut = 0.95;
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250300000","2444411044013300000","32c51070a","0103603n00000000","0153503000000000"); // fCosPAngleCut = 0.75;
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250600000","2444411044013300000","32c51070a","0103603n00000000","0153503000000000"); // fCosPAngleCut = 0.9;
  } else if( trainConfig == 180 ) { // PHOS 7 TeV
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","2444411044013300000","32c51070a","0103603n00000000","0453503000000000"); // without background
  } else if( trainConfig == 181 ) { // nonlin variation
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","2444421044012300000","32c51070a","0103603n00000000","0153503000000000");
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009227000008250400000","2444412044012300000","32c51070a","0103603n00000000","0153503000000000");
    // PHOS pp 5 TeV
  } else if(trainConfig == 190)  { // Standard PHOS  with TPC refit + ITS requirement
    cuts.AddCutHeavyMesonPCMCalo("00010113","00200009227000008250400000","2444411044012300000","32c01070a","0103603n00000000","0153503000000000"); // INT7
    cuts.AddCutHeavyMesonPCMCalo("00062113","00200009227000008250400000","2444411044012300000","32c01070a","0103603n00000000","0153503000000000"); // PHI7

    // PHOS LHC11 pp 7 TeV
  } else if(trainConfig == 195)  { // with TPC refit + ITS requirement
    cuts.AddCutHeavyMesonPCMCalo("00010113","00200009227000008250400000","2444400053012300000","32c010708","0103603n00000000","0153503000000000"); // INT7
    cuts.AddCutHeavyMesonPCMCalo("00062113","00200009227000008250400000","2444400053012300000","32c010708","0103603n00000000","0153503000000000"); // PHI7


  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //                                          ETA PRIME MESON
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  } else if( trainConfig == 200 ) {
    // everything open, min pt charged pi = 100 MeV
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009327000008250400000","1111100047032230000","000010400","0103503m00000000","0103503000000000");
  } else if( trainConfig == 201 ) {
    // closing charged pion cuts, minimum TPC cluster = 80, TPC dEdx pi = \pm 3 sigma, min pt charged pi = 100 MeV
    cuts.AddCutHeavyMesonPCMCalo("00000113","00200009327000008250400000","1111100047032230000","002010700","0103503m00000000","0103503000000000");
  } else if( trainConfig == 202) {
    // eta < 0.9
    // closing charged pion cuts, minimum TPC cluster = 80, TPC dEdx pi = \pm 3 sigma, pi+pi- mass cut of 1.5, min pt charged pi = 100 MeV
    // closing neural pion cuts, 0.5 < M_gamma,gamma < 0.6
    // maxChi2 per cluster TPC <4, require TPC refit, DCA XY pT dependend 0.0182+0.0350/pt^1.01, DCA_Z = 3.0
    // timing cluster cut open
    cuts.AddCutHeavyMesonPCMCalo("00010113","00200009327000008250400000","1111100047032230000","30a330709","0103503l00000000","0153503000000000"); // all of the above
    cuts.AddCutHeavyMesonPCMCalo("00052113","00200009327000008250400000","1111100047032230000","30a330709","0103503l00000000","0153503000000000"); // all of the above
    cuts.AddCutHeavyMesonPCMCalo("00062113","00200009327000008250400000","1111100047032230000","30a330709","0103503l00000000","0153503000000000"); // all of the above
    cuts.AddCutHeavyMesonPCMCalo("00083113","00200009327000008250400000","1111100047032230000","30a330709","0103503l00000000","0153503000000000"); // all of the above
    cuts.AddCutHeavyMesonPCMCalo("00085113","00200009327000008250400000","1111100047032230000","30a330709","0103503l00000000","0153503000000000"); // all of the above
  } else if( trainConfig == 203) {
    // same as 202 but only MB
    cuts.AddCutHeavyMesonPCMCalo("00010113","00200009327000008250400000","1111100047032230000","30a330709","0103503l00000000","0153503000000000"); // 0.5-0.6 eta mass cut
    cuts.AddCutHeavyMesonPCMCalo("00010113","00200009327000008250400000","1111100047032230000","30a330709","0103503m00000000","0153503000000000"); // 0.4-0.7 eta mass cut
  } else if( trainConfig == 204) {
    // same as 202 but with mass cut variations
    cuts.AddCutHeavyMesonPCMCalo("00010113","00200009327000008250400000","1111100047032230000","30a330700","0103503l00000000","0153503000000000"); // pi+pi- mass cut of 10
    cuts.AddCutHeavyMesonPCMCalo("00010113","00200009327000008250400000","1111100047032230000","30a330701","0103503l00000000","0153503000000000"); // pi+pi- mass cut of 1
    cuts.AddCutHeavyMesonPCMCalo("00010113","00200009327000008250400000","1111100047032230000","30a330708","0103503l00000000","0153503000000000"); // pi+pi- mass cut of 0.85
  } else if ( trainConfig == 205 ) {
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","411791106f032220000","32c510700","0103603l00000010","0153503000000010"); // INT7
    cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411791106f032220000","32c510700","0103603l00000010","0153503000000010"); // EMC7
    cuts.AddCutHeavyMesonPCMCalo("0008d113","0dm00009f9730000dge0404000","411791106f032220000","32c510700","0103603l00000010","0153503000000010"); // EMC7
    cuts.AddCutHeavyMesonPCMCalo("0009b113","0dm00009f9730000dge0404000","411791106f032220000","32c510700","0103603l00000010","0153503000000010"); // EMC7
  } else if ( trainConfig == 206 ) { // no event mixing
    cuts.AddCutHeavyMesonPCMCalo("00010113","00200009227000008250400000","411790106f032220000","32c510700","0103603l00000000","0453503000000000"); // INT7
    cuts.AddCutHeavyMesonPCMCalo("0008e113","00200009227000008250400000","411790106f032220000","32c510700","0103603l00000000","0453503000000000"); // EMC7
    cuts.AddCutHeavyMesonPCMCalo("0008d113","00200009227000008250400000","411790106f032220000","32c510700","0103603l00000000","0453503000000000"); // EMC7
    cuts.AddCutHeavyMesonPCMCalo("0009b113","00200009227000008250400000","411790106f032220000","32c510700","0103603l00000000","0453503000000000"); // EMC7
  } else if ( trainConfig == 207 ) { // no event mixing only tiggers
    cuts.AddCutHeavyMesonPCMCalo("0008e113","00200009227000008250400000","411790106f032220000","32c510700","0103603l00000000","0453503000000000"); // EMC7
    cuts.AddCutHeavyMesonPCMCalo("0008d113","00200009227000008250400000","411790106f032220000","32c510700","0103603l00000000","0453503000000000"); // EMC7

  // PCM-PHOS
  } else if ( trainConfig == 250 ) { // INT7 + PHI7
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","24466190wa01cc00000","32c510700","0103603l00000000","0453503000000000"); // INT7
    cuts.AddCutHeavyMesonPCMCalo("00062113","0dm00009f9730000dge0404000","24466190wa01cc00000","32c510700","0103603l00000000","0453503000000000"); // PHI7

  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //                                          OMEGA MESON pp 13 TeV
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // PHOS pp 13 TeV
  } else if(trainConfig == 400)  { // pp13 TeV AOD and ESD Comparison
    cuts.AddCutHeavyMesonPCMCalo("00010113","00200009227000008250400000","2444411044012300000","32c51070a","0103603400000000","0153503000000000"); // INT7
  } else if(trainConfig == 401)  { // Standard PHOS MB
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","24466190wa01cc00000","32c51070a","0103603400000000","0153503000000000"); // INT7
  } else if(trainConfig == 402)  { //Standard PHOS Trigger PHI7
    cuts.AddCutHeavyMesonPCMCalo("00062113","0dm00009f9730000dge0404000","24466190wa01cc00000","32c51070a","0103603400000000","0153503000000000"); // PHI7
  } else if(trainConfig == 405)  { //Standard EMCal 13TeV
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","411791106f032220000","32c51070a","0103603200000000","0153503000000000"); // INT7
  } else if(trainConfig == 406)  { //Standard EMCal 13TeV + Triggers
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","411791106f032220000","32c51070a","0103603200000000","0153503000000000"); // INT7
    cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411791106f032220000","32c51070a","0103603200000000","0153503000000000"); // PHI7
    cuts.AddCutHeavyMesonPCMCalo("0008d113","0dm00009f9730000dge0404000","411791106f032220000","32c51070a","0103603200000000","0153503000000000"); // INT7
    cuts.AddCutHeavyMesonPCMCalo("0009b113","0dm00009f9730000dge0404000","411791106f032220000","32c51070a","0103603200000000","0153503000000000"); // PHI7
  } else if(trainConfig == 407)  { //Standard EMCal 13TeV MB, testbeam nl
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","411790106f032220000","32c51070a","0103603200000000","0153503000000000"); // INT7
  } else if(trainConfig == 408)  { //Standard EMCal 13TeV Triggers, testbeam nl
    cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411790106f032220000","32c51070a","0103603200000000","0453503000000000"); // EG2
    cuts.AddCutHeavyMesonPCMCalo("0008d113","0dm00009f9730000dge0404000","411790106f032220000","32c51070a","0103603200000000","0453503000000000"); // EG1
  } else if(trainConfig == 409)  { //Standard EMCal 13TeV, no testbeam nl
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","411790006f032220000","32c51070a","0103603200000000","0153503000000000"); // INT7
    cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411790006f032220000","32c51070a","0103603200000000","0453503000000000"); // EG2
    cuts.AddCutHeavyMesonPCMCalo("0008d113","0dm00009f9730000dge0404000","411790006f032220000","32c51070a","0103603200000000","0453503000000000"); // EG1
  } else if(trainConfig == 410)  { //PHOS Trig Pt Cut Variations
      cuts.AddCutHeavyMesonPCMCalo("00062113","0dm00009f9730000dge0404000","24466190wa01cc00000","32c51070a","01036c3400000000","0153503000000000"); // PHI7, Pion 8 GeV
      cuts.AddCutHeavyMesonPCMCalo("00062113","0dm00009f9730000dge0404000","24466190wa01cc00000","32c51070a","01036g3400000000","0153503000000000"); // PHI7, new Gamma Energy cut 6 GeV
  } else if(trainConfig == 411)  { //EMCal Trig Pt Cut Variations EG2
      cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411790106f032220000","32c51070a","01036c3200000000","0453503000000000"); // EG2, Pion 8 GeV
      cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411790106f032220000","32c51070a","01036g3200000000","0453503000000000"); // EG2, new Gamma Energy cut 6 GeV
  } else if(trainConfig == 412)  { //EMCal Trig Pt Cut Variations EG1
      cuts.AddCutHeavyMesonPCMCalo("0008d113","0dm00009f9730000dge0404000","411790106f032220000","32c51070a","01036q3200000000","0453503000000000"); // EG1, Pion 12 GeV
      cuts.AddCutHeavyMesonPCMCalo("0008d113","0dm00009f9730000dge0404000","411790106f032220000","32c51070a","01036h3200000000","0453503000000000"); // EG1, new Gamma Energy cut 10 GeV
  } else if(trainConfig == 413)  { //PHOS Trig Pt Cut Variations
      cuts.AddCutHeavyMesonPCMCalo("00062113","0dm00009f9730000dge0404000","24466190wa01cc00000","32c51070a","01036e3400000000","0153503000000000"); // PHI7, new Gamma Energy cut 5. GeV
      cuts.AddCutHeavyMesonPCMCalo("00062113","0dm00009f9730000dge0404000","24466190wa01cc00000","32c51070a","01036f3400000000","0153503000000000"); // PHI7, new Gamma Energy cut 7.5 GeV
  } else if(trainConfig == 414)  { //EMCal Trig Pt Cut Variations EG2
      cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411790106f032220000","32c51070a","01036e3200000000","0453503000000000"); // EG2, new Gamma Energy cut 5 GeV
      cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411790106f032220000","32c51070a","01036f3200000000","0453503000000000"); // EG2, new Gamma Energy cut 7.5 GeV
  } else if(trainConfig == 415)  { //EMCal Trig Pt Cut Variations EG1
      cuts.AddCutHeavyMesonPCMCalo("0008d113","0dm00009f9730000dge0404000","411790106f032220000","32c51070a","01036f3200000000","0453503000000000"); // EG1, new Gamma Energy cut 7.5 GeV
      cuts.AddCutHeavyMesonPCMCalo("0008d113","0dm00009f9730000dge0404000","411790106f032220000","32c51070a","01036i3200000000","0453503000000000"); // EG1, new Gamma Energy cut 12 GeV
  // Variations on 5 TeV for 7 TeV systematics
  // PCM-EMC (without nonlin)
  } else if(trainConfig == 450)  { //Standard PCM-EMCal 13TeV, no testbeam nl
    cuts.AddCutHeavyMesonPCMCalo("00010113","00200009f9730000dge0400000","411791106f032220000","32c51070a","0103603700000000","0153503000000000"); // INT7
  } else if(trainConfig == 451)  { // mass selection window
    cuts.AddCutHeavyMesonPCMCalo("00010113","00200009f9730000dge0400000","411791106f032220000","32c51070a","0103603e00000000","0153503000000000"); // 1 sigma
    cuts.AddCutHeavyMesonPCMCalo("00010113","00200009f9730000dge0400000","411791106f032220000","32c51070a","0103603f00000000","0153503000000000"); // 3 sigma
    cuts.AddCutHeavyMesonPCMCalo("00010113","00200009f9730000dge0400000","411791106f032220000","32c51070a","0103603g00000000","0153503000000000"); // 4 sigma
  } else if(trainConfig == 452)  { // background description
    cuts.AddCutHeavyMesonPCMCalo("00010113","00200009f9730000dge0400000","411791106f032220000","32c51070a","0103603700000000","0a53503000000000"); // likesign mixing
    cuts.AddCutHeavyMesonPCMCalo("00010113","00200009f9730000dge0400000","411791106f032220000","32c51070a","0103603700000000","0d53503000000000"); // sideband migxin
  } else if(trainConfig == 453)  { // track matching
    cuts.AddCutHeavyMesonPCMCalo("00010113","00200009f9730000dge0400000","411791106f032220000","32c51070a","0103603700000000","0153503000000000"); // INT7
    cuts.AddCutHeavyMesonPCMCalo("00010113","00200009f9730000dge0400000","4117911061032220000","32c51070a","0103603700000000","0153503000000000"); // INT7
    cuts.AddCutHeavyMesonPCMCalo("00010113","00200009f9730000dge0400000","4117911063032220000","32c51070a","0103603700000000","0153503000000000"); // INT7
    cuts.AddCutHeavyMesonPCMCalo("00010113","00200009f9730000dge0400000","4117911065032220000","32c51070a","0103603700000000","0153503000000000"); // INT7
    cuts.AddCutHeavyMesonPCMCalo("00010113","00200009f9730000dge0400000","4117911066032220000","32c51070a","0103603700000000","0153503000000000"); // INT7
  } else if(trainConfig == 454)  { // TPC dEdxCutPion
    cuts.AddCutHeavyMesonPCMCalo("00010113","00200009f9730000dge0400000","411791106f032220000","32c51080a","0103603700000000","0153503000000000");
    cuts.AddCutHeavyMesonPCMCalo("00010113","00200009f9730000dge0400000","411791106f032220000","32c51020a","0103603700000000","0153503000000000");
    cuts.AddCutHeavyMesonPCMCalo("00010113","00200009f9730000dge0400000","411791106f032220000","32c51030a","0103603700000000","0153503000000000");
    cuts.AddCutHeavyMesonPCMCalo("00010113","00200009f9730000dge0400000","411791106f032220000","32c51040a","0103603700000000","0153503000000000");
  // Variations on 13 TeV for 7 TeV systematics
  // PCM-PHOS (without nonlin)
  } else if(trainConfig == 500)  { // Standard PHOS
    cuts.AddCutHeavyMesonPCMCalo("00010113","00200009f9730000dge0400000","24466510ga012200000","32c51070a","0103603900000000","0153503000000000"); // INT7
  } else if(trainConfig == 501)  { // Selection window cut
    cuts.AddCutHeavyMesonPCMCalo("00010113","00200009f9730000dge0400000","24466510ga012200000","32c51070a","0103603h00000000","0153503000000000"); // 1 sigma
    cuts.AddCutHeavyMesonPCMCalo("00010113","00200009f9730000dge0400000","24466510ga012200000","32c51070a","0103603i00000000","0153503000000000"); // 3 sigma
    cuts.AddCutHeavyMesonPCMCalo("00010113","00200009f9730000dge0400000","24466510ga012200000","32c51070a","0103603j00000000","0153503000000000"); // 4 sigma

  // *************Variations in AliCaloPhotonsCut**************************
  } else if(trainConfig == 502)  { // track matching
    cuts.AddCutHeavyMesonPCMCalo("00010113","00200009f9730000dge0400000","24466510g1012200000","32c51070a","0103603900000000","0153503000000000"); // INT7
    cuts.AddCutHeavyMesonPCMCalo("00010113","00200009f9730000dge0400000","24466510g3012200000","32c51070a","0103603900000000","0153503000000000"); // INT7
    cuts.AddCutHeavyMesonPCMCalo("00010113","00200009f9730000dge0400000","24466510g5012200000","32c51070a","0103603900000000","0153503000000000"); // INT7
    cuts.AddCutHeavyMesonPCMCalo("00010113","00200009f9730000dge0400000","24466510g6012200000","32c51070a","0103603900000000","0153503000000000"); // INT7
    cuts.AddCutHeavyMesonPCMCalo("00010113","00200009f9730000dge0400000","24466510g7012200000","32c51070a","0103603900000000","0153503000000000"); // INT7
    cuts.AddCutHeavyMesonPCMCalo("00010113","00200009f9730000dge0400000","24466510g9012200000","32c51070a","0103603900000000","0153503000000000"); // INT7
  } else if(trainConfig == 503)  { // min energy clusters
    cuts.AddCutHeavyMesonPCMCalo("00010113","00200009f9730000dge0400000","24466510ga022200000","32c51070a","0103603900000000","0153503000000000"); // INT7
    cuts.AddCutHeavyMesonPCMCalo("00010113","00200009f9730000dge0400000","24466510ga032200000","32c51070a","0103603900000000","0153503000000000"); // INT7
    cuts.AddCutHeavyMesonPCMCalo("00010113","00200009f9730000dge0400000","24466510ga042200000","32c51070a","0103603900000000","0153503000000000"); // INT7
    cuts.AddCutHeavyMesonPCMCalo("00010113","00200009f9730000dge0400000","24466510ga082200000","32c51070a","0103603900000000","0153503000000000"); // INT7
  } else if(trainConfig == 504)  { // nmb of cells
    cuts.AddCutHeavyMesonPCMCalo("00010113","00200009f9730000dge0400000","24466510ga011200000","32c51070a","0103603900000000","0153503000000000"); // INT7
    cuts.AddCutHeavyMesonPCMCalo("00010113","00200009f9730000dge0400000","24466510ga013200000","32c51070a","0103603900000000","0153503000000000"); // INT7
    cuts.AddCutHeavyMesonPCMCalo("00010113","00200009f9730000dge0400000","24466510ga014200000","32c51070a","0103603900000000","0153503000000000"); // INT7
  } else if(trainConfig == 505)  { // minmax M02
    cuts.AddCutHeavyMesonPCMCalo("00010113","00200009f9730000dge0400000","24466510ga012100000","32c51070a","0103603900000000","0153503000000000"); // INT7
    cuts.AddCutHeavyMesonPCMCalo("00010113","00200009f9730000dge0400000","24466510ga012300000","32c51070a","0103603900000000","0153503000000000"); // INT7
  } else {
    Error(Form("GammaConvNeutralMeson_MixedMode_%i",trainConfig), "wrong trainConfig variable no cuts have been specified for the configuration");
    return;
  }

  if(!cuts.AreValid()){
    cout << "\n\n****************************************************" << endl;
    cout << "ERROR: No valid cuts stored in CutHandlerNeutralMixed! Returning..." << endl;
    cout << "****************************************************\n\n" << endl;
    return;
  }

  Int_t numberOfCuts = cuts.GetNCuts();

  TList *EventCutList = new TList();
  TList *ConvCutList  = new TList();
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
  ConvCutList->SetOwner(kTRUE);
  AliConversionPhotonCuts **analysisCuts = new AliConversionPhotonCuts*[numberOfCuts];
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

    analysisCuts[i] = new AliConversionPhotonCuts();
    analysisCuts[i]->SetV0ReaderName(V0ReaderName);

    // post calibration of dEdx energy loss
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

    if(runLightOutput>0) analysisCuts[i]->SetLightOutput(kTRUE);
    if( ! analysisCuts[i]->InitializeCutsFromCutString((cuts.GetPhotonCut(i)).Data()) ) {
      cout<<"ERROR: analysisCuts [" <<i<<"]"<<endl;
      return 0;
    } else {
      ConvCutList->Add(analysisCuts[i]);
      analysisCuts[i]->SetFillCutHistograms("",kFALSE);
    }

    analysisClusterCuts[i] = new AliCaloPhotonCuts();
    analysisClusterCuts[i]->SetV0ReaderName(V0ReaderName);
    if(runLightOutput>0) analysisClusterCuts[i]->SetLightOutput(kTRUE);
    analysisClusterCuts[i]->SetCaloTrackMatcherName(TrackMatcherName);
    analysisClusterCuts[i]->SetExtendedMatchAndQA(enableExtMatchAndQA);
    if( ! analysisClusterCuts[i]->InitializeCutsFromCutString((cuts.GetClusterCut(i)).Data()) ) {
      cout<<"ERROR: analysisClusterCuts [" <<i<<"]"<<endl;
      return 0;
    } else {
      ClusterCutList->Add(analysisClusterCuts[i]);
      analysisClusterCuts[i]->SetFillCutHistograms("");
    }

    analysisNeutralPionCuts[i] = new AliConversionMesonCuts();
    analysisNeutralPionCuts[i]->SetUsePtDepSelectionWindow(usePtDepSelectionWindowCut);
    if(runLightOutput>0) analysisNeutralPionCuts[i]->SetLightOutput(kTRUE);
    if( ! analysisNeutralPionCuts[i]->InitializeCutsFromCutString((cuts.GetNDMCut(i)).Data()) ) {
      cout<<"ERROR: analysisMesonCuts [ " <<i<<" ] "<<endl;
      return 0;
    } else {
      NeutralPionCutList->Add(analysisNeutralPionCuts[i]);
      analysisNeutralPionCuts[i]->SetFillCutHistograms("");
    }

    analysisMesonCuts[i] = new AliConversionMesonCuts();
    if(runLightOutput>0) analysisMesonCuts[i]->SetLightOutput(kTRUE);
    if( ! analysisMesonCuts[i]->InitializeCutsFromCutString((cuts.GetMesonCut(i)).Data()) ) {
      cout<<"ERROR: analysisMesonCuts [ " <<i<<" ] "<<endl;
      return 0;
    } else {
      MesonCutList->Add(analysisMesonCuts[i]);
      analysisMesonCuts[i]->SetFillCutHistograms("");
    }
    analysisEventCuts[i]->SetAcceptedHeader(HeaderList);

    TString cutName( Form("%s_%s_%s_%s_%s_%s",(cuts.GetEventCut(i)).Data(), (cuts.GetPhotonCut(i)).Data(), (cuts.GetClusterCut(i)).Data(),(cuts.GetPionCut(i)).Data(),(cuts.GetNDMCut(i)).Data(), (cuts.GetMesonCut(i)).Data() ) );
    analysisPionCuts[i] = new AliPrimaryPionCuts();
    analysisPionCuts[i]->SetPrefilterRunFlag(prefilterRunFlag);
    analysisPionCuts[i]->SetPeriodName(periodNameV0Reader);
    if(runLightOutput>0) analysisPionCuts[i]->SetLightOutput(kTRUE);

    if( !analysisPionCuts[i]->InitializeCutsFromCutString((cuts.GetPionCut(i)).Data())) {
      cout<< "ERROR:  analysisPionCuts [ " <<i<<" ] "<<endl;
      return 0;
    } else {
      PionCutList->Add(analysisPionCuts[i]);
      analysisPionCuts[i]->SetFillCutHistograms("",kFALSE,cutName);
    }
  }

  task->SetNDMRecoMode(neutralPionMode);
  task->SetEventCutList(numberOfCuts,EventCutList);
  task->SetConversionCutList(ConvCutList);
  task->SetClusterCutList(ClusterCutList);
  task->SetNeutralPionCutList(NeutralPionCutList);
  task->SetMesonCutList(MesonCutList);
  task->SetPionCutList(PionCutList);

  task->SetMoveParticleAccordingToVertex(kTRUE);
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
