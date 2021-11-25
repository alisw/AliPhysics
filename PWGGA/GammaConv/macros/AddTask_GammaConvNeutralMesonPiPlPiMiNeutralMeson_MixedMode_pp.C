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
    // settings for weights
    // FPTW:fileNamePtWeights, FMUW:fileNameMultWeights,  FMAW:fileNameMatBudWeights,  separate with ;
    // Material Budget Weights file for Run 2
    // FMAW:alien:///alice/cern.ch/user/a/amarin//MBW/MCInputFileMaterialBudgetWeightsLHC16_Pythia_00010103_0d000009266300008850404000_date181214.root
    TString   fileNameExternalInputs        = "MCSpectraInput.root",    // path to file for weigting input
    Bool_t    doWeighting                   = kFALSE,                   //enable Weighting
    Bool_t    enableElecDeDxPostCalibration = kFALSE,                   // enable post calibration of elec pos dEdX
    TString   generatorName               = "HIJING",
    Double_t  tolerance                   = -1,
    TString   periodNameV0Reader          = "",                       // period Name for V0Reader
    Int_t     runLightOutput              = 0,                        // run light output option 0: no light output 1: most cut histos stiched off 2: unecessary omega hists turned off as well
    Int_t     prefilterRunFlag            = 1500,                     // flag to change the prefiltering of ESD tracks. See SetHybridTrackCutsAODFiltering() in AliPrimaryPionCuts
    Bool_t    usePtDepSelectionWindowCut  = kFALSE,                   // use pt dependent meson selection window cut
    Bool_t    enableSortingMCLabels       = kTRUE,                    // enable sorting for MC cluster labels
    Int_t     enableMatBudWeightsPi0      = 0,                        // 1 = three radial bins, 2 = 10 radial bins (2 is the default when using weights)
    TString   additionalTrainConfig       = "0"                       // additional counter for trainconfig, this has to be always the last parameter
  ) {

  AliCutHandlerPCM cuts(13);
  TString addTaskName                       = "AddTask_GammaConvNeutralMesonPiPlPiMiNeutralMeson_MixedMode_pp";
  TString fileNamePtWeights                 = cuts.GetSpecialFileNameFromString (fileNameExternalInputs, "FPTW:");
  TString fileNameMultWeights               = cuts.GetSpecialFileNameFromString (fileNameExternalInputs, "FMUW:");
  TString fileNameMatBudWeights             = cuts.GetSpecialFileNameFromString (fileNameExternalInputs, "FMAW:");
  TString fileNamedEdxPostCalib             = cuts.GetSpecialFileNameFromString (fileNameExternalInputs, "FEPC:");
  TString fileNameCustomTriggerMimicOADB    = cuts.GetSpecialFileNameFromString (fileNameExternalInputs, "FTRM:");

  if(additionalTrainConfig.Contains("MaterialBudgetWeights"))
    fileNameMatBudWeights         = cuts.GetSpecialSettingFromAddConfig(additionalTrainConfig, "MaterialBudgetWeights",fileNameMatBudWeights, addTaskName);

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


    //PCM-EMC configs with newest NL (apply only with cell scale!)
  } else if(trainConfig == 210)  { //EDC 13TeV MB, NCell >=2 + NCell efficiency
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe3n220000","32c51070m","0103603l00000000","0453503000000000"); // INT7
  } else if(trainConfig == 211)  { //EDC 13TeV MB, NCell >=2 + NCell efficiency
    cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411790109fe3n220000","32c51070m","0103603l00000000","0453503000000000"); // INT7
  } else if(trainConfig == 212)  { //EDC 13TeV MB, NCell >=2 + NCell efficiency
    cuts.AddCutHeavyMesonPCMCalo("0008d113","0dm00009f9730000dge0404000","411790109fe3n220000","32c51070m","0103603l00000000","0453503000000000"); // INT7
    //PCM-EMC configs with newest NL (apply only with cell scale!)
  } else if(trainConfig == 213)  { //EDC 13TeV MB, NCell >=1
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe30220000","32c51070m","0103603l00000000","0453503000000000"); // INT7
  } else if(trainConfig == 214)  { //EDC 13TeV MB, NCell >=1
    cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411790109fe30220000","32c51070m","0103603l00000000","0453503000000000"); // INT7
  } else if(trainConfig == 215)  { //EDC 13TeV MB, NCell >=1
    cuts.AddCutHeavyMesonPCMCalo("0008d113","0dm00009f9730000dge0404000","411790109fe30220000","32c51070m","0103603l00000000","0453503000000000"); // INT7

  // PCM-PHOS
  } else if ( trainConfig == 250 ) { // INT7 + PHI7
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","24466190wa01cc00000","32c510700","0103603l00000000","0453503000000000"); // INT7
    cuts.AddCutHeavyMesonPCMCalo("00062113","0dm00009f9730000dge0404000","24466190wa01cc00000","32c510700","0103603l00000000","0453503000000000"); // PHI7


  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //                                          OMEGA MESON pp 13 TeV
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // PHOS pp 13 TeV
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  //Standard Cuts
  } else if(trainConfig == 400)  { // Standard PHOS MB
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa01cc00000","32c51070a","0103103400000000","0453503000000000"); // INT7
  } else if(trainConfig == 401)  { //Standard PHOS Trigger PHI7
    cuts.AddCutHeavyMesonPCMCalo("00062113","0dm00009f9730000dge0404000","24466190sa01cc00000","32c51070a","0103103400000000","0453503000000000"); // PHI7

  //Gamma Energy Cuts
  } else if(trainConfig == 403)  { //Standard PHOS 13TeV Trigger PHI7, GammaCut > 4GeV
    cuts.AddCutHeavyMesonPCMCalo("00062113","0dm00009f9730000dge0404000","24466190sa01cc00000","32c51070a","01031k3400000000","0453503000000000"); // PHI7, new Gamma Energy cut 4 GeV
  } else if(trainConfig == 404)  { //Standard PHOS 13TeV Trigger PHI7, GammaCut > 5GeV
    cuts.AddCutHeavyMesonPCMCalo("00062113","0dm00009f9730000dge0404000","24466190sa01cc00000","32c51070a","01031e3400000000","0453503000000000"); // PHI7, new Gamma Energy cut 5 GeV
  } else if(trainConfig == 405)  { //Standard PHOS 13TeV Trigger PHI7, GammaCut > 6GeV
    cuts.AddCutHeavyMesonPCMCalo("00062113","0dm00009f9730000dge0404000","24466190sa01cc00000","32c51070a","01031g3400000000","0453503000000000"); // PHI7, new Gamma Energy cut 6 GeV
  } else if(trainConfig == 406)  { //Standard PHOS 13TeV Trigger PHI7, GammaCut > 7.5GeV
    cuts.AddCutHeavyMesonPCMCalo("00062113","0dm00009f9730000dge0404000","24466190sa01cc00000","32c51070a","01031f3400000000","0453503000000000"); // PHI7, new Gamma Energy cut 7.5 GeV

  //QA Plots
  } else if(trainConfig == 415)  { // Standard PHOS MB, no shared TPC clusters, Shared cluster Fraction =0
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa01cc00000","32e51070a","0103103400000000","0453503000000000"); // INT7
  } else if(trainConfig == 416)  { // Standard PHOS MB, no shared TPC clusters, Shared cluster Fraction <=0.4
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa01cc00000","32f51070a","0103103400000000","0453503000000000"); // INT7


    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // EMC pp 13 TeV
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    //Standard Cuts
  } else if(trainConfig == 430)  { //Standard EMCal 13TeV
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","411792109fe32220000","32c51070a","0103103200000000","0453503000000000"); // INT7
  } else if(trainConfig == 431)  { // EDC 13 TeV EG2
    cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411792109fe32220000","32c51070a","0103103200000000","0453503000000000"); // EG2
  } else if(trainConfig == 432)  { // EDC 13 TeV EG1
    cuts.AddCutHeavyMesonPCMCalo("0008d113","0dm00009f9730000dge0404000","411792109fe32220000","32c51070a","0103103200000000","0453503000000000"); // EG1

    //Gamma Energy Cuts EG2
  } else if(trainConfig == 433)  { //Standard EMCal 13TeV Trigger EG2, GammaCut > 4GeV
    cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411792109fe32220000","32c51070a","01031k3200000000","0453503000000000"); // EG2, new Gamma Energy cut 4 GeV
  } else if(trainConfig == 434)  { //Standard EMCal 13TeV Trigger EG2, GammaCut > 5GeV
    cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411792109fe32220000","32c51070a","01031e3200000000","0453503000000000"); // EG2, new Gamma Energy cut 5 GeV
  } else if(trainConfig == 435)  { //Standard EMCal 13TeV Trigger EG2, GammaCut > 6GeV
    cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411792109fe32220000","32c51070a","01031g3200000000","0453503000000000"); // EG2, new Gamma Energy cut 6 GeV
  } else if(trainConfig == 436)  { //Standard EMCal 13TeV Trigger EG2, GammaCut > 7.5GeV
    cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411792109fe32220000","32c51070a","01031f3200000000","0453503000000000"); // EG2, new Gamma Energy cut 7.5 GeV

    //Gamma Energy Cuts EG1
  } else if(trainConfig == 437)  { //Standard EMCal 13TeV Trigger EG1, GammaCut > 8GeV
    cuts.AddCutHeavyMesonPCMCalo("0008d113","0dm00009f9730000dge0404000","411792109fe32220000","32c51070a","01031l3200000000","0453503000000000"); // EG1, new Gamma Energy cut 8 GeV
  } else if(trainConfig == 438)  { //Standard EMCal 13TeV Trigger EG1, GammaCut > 9GeV
    cuts.AddCutHeavyMesonPCMCalo("0008d113","0dm00009f9730000dge0404000","411792109fe32220000","32c51070a","01031m3200000000","0453503000000000"); // EG1, new Gamma Energy cut 9 GeV
  } else if(trainConfig == 439)  { //Standard EMCal 13TeV Trigger EG1, GammaCut > 10GeV
    cuts.AddCutHeavyMesonPCMCalo("0008d113","0dm00009f9730000dge0404000","411792109fe32220000","32c51070a","01031h3200000000","0453503000000000"); // EG1, new Gamma Energy cut 10 GeV

    //no shared TPC clusters, Shared cluster Fraction =0
  } else if(trainConfig == 440)  { //Standard EMCal 13TeV, Shared cluster Fraction =0
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","411792106fe32220000","32e51070a","0103103200000000","0453503000000000"); // INT7
  } else if(trainConfig == 441)  { //Standard EMCal 13TeV Trigger EG2, GammaCut > 6GeV, Shared cluster Fraction =0
    cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411792106fe32220000","32e51070a","01031g3200000000","0453503000000000"); // EG2, new Gamma Energy cut 6 GeV
  } else if(trainConfig == 442)  { //Standard EMCal 13TeV Trigger EG1, GammaCut > 10GeV, Shared cluster Fraction =0
    cuts.AddCutHeavyMesonPCMCalo("0008d113","0dm00009f9730000dge0404000","411792106fe32220000","32e51070a","01031h3200000000","0453503000000000"); // EG1, new Gamma Energy cut 10 GeV

    //no shared TPC clusters, Shared cluster Fraction <=0.4
  } else if(trainConfig == 443)  { //Standard EMCal 13TeV, Shared cluster Fraction <=0.4
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","411792106fe32220000","32f51070a","0103103200000000","0453503000000000"); // INT7
  } else if(trainConfig == 444)  { //Standard EMCal 13TeV Trigger EG2, GammaCut > 6GeV, Shared cluster Fraction <=0.4
    cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411792106fe32220000","32f51070a","01031g3200000000","0453503000000000"); // EG2, new Gamma Energy cut 6 GeV
  } else if(trainConfig == 445)  { //Standard EMCal 13TeV Trigger EG1, GammaCut > 10GeV, Shared cluster Fraction <=0.4
    cuts.AddCutHeavyMesonPCMCalo("0008d113","0dm00009f9730000dge0404000","411792106fe32220000","32f51070a","01031h3200000000","0453503000000000"); // EG1, new Gamma Energy cut 10 GeV

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
  } else if(trainConfig == 600)  { // Standard PHOS MB
      cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa01cc00000","32c51070a","0103103400000000","0153503000000000"); // INT7
  } else if(trainConfig == 601)  { //Standard PHOS Trigger PHI7
      cuts.AddCutHeavyMesonPCMCalo("00062113","0dm00009f9730000dge0404000","24466190sa01cc00000","32c51070a","0103103400000000","0153503000000000"); // PHI7

      //Gamma Energy Cuts
  } else if(trainConfig == 603)  { //Standard PHOS 13TeV Trigger PHI7, GammaCut > 4GeV
      cuts.AddCutHeavyMesonPCMCalo("00062113","0dm00009f9730000dge0404000","24466190sa01cc00000","32c51070a","01031k3400000000","0153503000000000"); // PHI7, new Gamma Energy cut 4 GeV
  } else if(trainConfig == 604)  { //Standard PHOS 13TeV Trigger PHI7, GammaCut > 5GeV
      cuts.AddCutHeavyMesonPCMCalo("00062113","0dm00009f9730000dge0404000","24466190sa01cc00000","32c51070a","01031e3400000000","0153503000000000"); // PHI7, new Gamma Energy cut 5 GeV
  } else if(trainConfig == 605)  { //Standard PHOS 13TeV Trigger PHI7, GammaCut > 6GeV
      cuts.AddCutHeavyMesonPCMCalo("00062113","0dm00009f9730000dge0404000","24466190sa01cc00000","32c51070a","01031g3400000000","0153503000000000"); // PHI7, new Gamma Energy cut 6 GeV
  } else if(trainConfig == 606)  { //Standard PHOS 13TeV Trigger PHI7, GammaCut > 7.5GeV
      cuts.AddCutHeavyMesonPCMCalo("00062113","0dm00009f9730000dge0404000","24466190sa01cc00000","32c51070a","01031f3400000000","0153503000000000"); // PHI7, new Gamma Energy cut 7.5 GeV

      //QA Plots
  } else if(trainConfig == 615)  { // Standard PHOS MB, no shared TPC clusters, Shared cluster Fraction =0
      cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa01cc00000","32e51070a","0103103400000000","0153503000000000"); // INT7
  } else if(trainConfig == 616)  { // Standard PHOS MB, no shared TPC clusters, Shared cluster Fraction <=0.4
      cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa01cc00000","32f51070a","0103103400000000","0153503000000000"); // INT7


      // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      // EMC pp 13 TeV
      // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      //Standard Cuts
  } else if(trainConfig == 630)  { //Standard EMCal 13TeV
      cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","411792109fe32220000","32c51070a","0103103200000000","0153503000000000"); // INT7
  } else if(trainConfig == 631)  { // EDC 13 TeV EG2
      cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411792109fe32220000","32c51070a","0103103200000000","0153503000000000"); // EG2
  } else if(trainConfig == 632)  { // EDC 13 TeV EG1
      cuts.AddCutHeavyMesonPCMCalo("0008d113","0dm00009f9730000dge0404000","411792109fe32220000","32c51070a","0103103200000000","0153503000000000"); // EG1

      //Gamma Energy Cuts EG2
  } else if(trainConfig == 633)  { //Standard EMCal 13TeV Trigger EG2, GammaCut > 4GeV
      cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411792109fe32220000","32c51070a","01031k3200000000","0153503000000000"); // EG2, new Gamma Energy cut 4 GeV
  } else if(trainConfig == 634)  { //Standard EMCal 13TeV Trigger EG2, GammaCut > 5GeV
      cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411792109fe32220000","32c51070a","01031e3200000000","0153503000000000"); // EG2, new Gamma Energy cut 5 GeV
  } else if(trainConfig == 635)  { //Standard EMCal 13TeV Trigger EG2, GammaCut > 6GeV
      cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411792109fe32220000","32c51070a","01031g3200000000","0153503000000000"); // EG2, new Gamma Energy cut 6 GeV
  } else if(trainConfig == 636)  { //Standard EMCal 13TeV Trigger EG2, GammaCut > 7.5GeV
      cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411792109fe32220000","32c51070a","01031f3200000000","0153503000000000"); // EG2, new Gamma Energy cut 7.5 GeV

      //Gamma Energy Cuts EG1
  } else if(trainConfig == 637)  { //Standard EMCal 13TeV Trigger EG1, GammaCut > 8GeV
      cuts.AddCutHeavyMesonPCMCalo("0008d113","0dm00009f9730000dge0404000","411792109fe32220000","32c51070a","01031l3200000000","0153503000000000"); // EG1, new Gamma Energy cut 8 GeV
  } else if(trainConfig == 638)  { //Standard EMCal 13TeV Trigger EG1, GammaCut > 9GeV
      cuts.AddCutHeavyMesonPCMCalo("0008d113","0dm00009f9730000dge0404000","411792109fe32220000","32c51070a","01031m3200000000","0153503000000000"); // EG1, new Gamma Energy cut 9 GeV
  } else if(trainConfig == 639)  { //Standard EMCal 13TeV Trigger EG1, GammaCut > 10GeV
      cuts.AddCutHeavyMesonPCMCalo("0008d113","0dm00009f9730000dge0404000","411792109fe32220000","32c51070a","01031h3200000000","0153503000000000"); // EG1, new Gamma Energy cut 10 GeV

      //no shared TPC clusters, Shared cluster Fraction =0
  } else if(trainConfig == 640)  { //Standard EMCal 13TeV, Shared cluster Fraction =0
      cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","411792106fe32220000","32e51070a","0103103200000000","0153503000000000"); // INT7
  } else if(trainConfig == 641)  { //Standard EMCal 13TeV Trigger EG2, GammaCut > 6GeV, Shared cluster Fraction =0
      cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411792106fe32220000","32e51070a","01031g3200000000","0153503000000000"); // EG2, new Gamma Energy cut 6 GeV
  } else if(trainConfig == 642)  { //Standard EMCal 13TeV Trigger EG1, GammaCut > 10GeV, Shared cluster Fraction =0
      cuts.AddCutHeavyMesonPCMCalo("0008d113","0dm00009f9730000dge0404000","411792106fe32220000","32e51070a","01031h3200000000","0153503000000000"); // EG1, new Gamma Energy cut 10 GeV

      //no shared TPC clusters, Shared cluster Fraction <=0.4
  } else if(trainConfig == 643)  { //Standard EMCal 13TeV, Shared cluster Fraction <=0.4
      cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","411792106fe32220000","32f51070a","0103103200000000","0153503000000000"); // INT7
  } else if(trainConfig == 644)  { //Standard EMCal 13TeV Trigger EG2, GammaCut > 6GeV, Shared cluster Fraction <=0.4
      cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411792106fe32220000","32f51070a","01031g3200000000","0153503000000000"); // EG2, new Gamma Energy cut 6 GeV
  } else if(trainConfig == 645)  { //Standard EMCal 13TeV Trigger EG1, GammaCut > 10GeV, Shared cluster Fraction <=0.4
      cuts.AddCutHeavyMesonPCMCalo("0008d113","0dm00009f9730000dge0404000","411792106fe32220000","32f51070a","01031h3200000000","0153503000000000"); // EG1, new Gamma Energy cut 10 GeV
  } else if(trainConfig == 700)  { // Standard PHOS MB
      cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa01cc00000","32c51070a","0103103400000000","0a53503000000000"); // INT7
  } else if(trainConfig == 701)  { //Standard PHOS Trigger PHI7
      cuts.AddCutHeavyMesonPCMCalo("00062113","0dm00009f9730000dge0404000","24466190sa01cc00000","32c51070a","0103103400000000","0a53503000000000"); // PHI7

      //Gamma Energy Cuts
  } else if(trainConfig == 703)  { //Standard PHOS 13TeV Trigger PHI7, GammaCut > 4GeV
      cuts.AddCutHeavyMesonPCMCalo("00062113","0dm00009f9730000dge0404000","24466190sa01cc00000","32c51070a","01031k3400000000","0a53503000000000"); // PHI7, new Gamma Energy cut 4 GeV
  } else if(trainConfig == 704)  { //Standard PHOS 13TeV Trigger PHI7, GammaCut > 5GeV
      cuts.AddCutHeavyMesonPCMCalo("00062113","0dm00009f9730000dge0404000","24466190sa01cc00000","32c51070a","01031e3400000000","0a53503000000000"); // PHI7, new Gamma Energy cut 5 GeV
  } else if(trainConfig == 705)  { //Standard PHOS 13TeV Trigger PHI7, GammaCut > 6GeV
      cuts.AddCutHeavyMesonPCMCalo("00062113","0dm00009f9730000dge0404000","24466190sa01cc00000","32c51070a","01031g3400000000","0a53503000000000"); // PHI7, new Gamma Energy cut 6 GeV
  } else if(trainConfig == 706)  { //Standard PHOS 13TeV Trigger PHI7, GammaCut > 7.5GeV
      cuts.AddCutHeavyMesonPCMCalo("00062113","0dm00009f9730000dge0404000","24466190sa01cc00000","32c51070a","01031f3400000000","0a53503000000000"); // PHI7, new Gamma Energy cut 7.5 GeV

      //QA Plots
  } else if(trainConfig == 715)  { // Standard PHOS MB, no shared TPC clusters, Shared cluster Fraction =0
      cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa01cc00000","32e51070a","0103103400000000","0a53503000000000"); // INT7
  } else if(trainConfig == 716)  { // Standard PHOS MB, no shared TPC clusters, Shared cluster Fraction <=0.4
      cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa01cc00000","32f51070a","0103103400000000","0a53503000000000"); // INT7


      // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      // EMC pp 13 TeV
      // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      //Standard Cuts
  } else if(trainConfig == 730)  { //Standard EMCal 13TeV
      cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","411792109fe32220000","32c51070a","0103103200000000","0a53503000000000"); // INT7
  } else if(trainConfig == 731)  { // EDC 13 TeV EG2
      cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411792109fe32220000","32c51070a","0103103200000000","0a53503000000000"); // EG2
  } else if(trainConfig == 732)  { // EDC 13 TeV EG1
      cuts.AddCutHeavyMesonPCMCalo("0008d113","0dm00009f9730000dge0404000","411792109fe32220000","32c51070a","0103103200000000","0a53503000000000"); // EG1

      //Gamma Energy Cuts EG2
  } else if(trainConfig == 733)  { //Standard EMCal 13TeV Trigger EG2, GammaCut > 4GeV
      cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411792109fe32220000","32c51070a","01031k3200000000","0a53503000000000"); // EG2, new Gamma Energy cut 4 GeV
  } else if(trainConfig == 734)  { //Standard EMCal 13TeV Trigger EG2, GammaCut > 5GeV
      cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411792109fe32220000","32c51070a","01031e3200000000","0a53503000000000"); // EG2, new Gamma Energy cut 5 GeV
  } else if(trainConfig == 735)  { //Standard EMCal 13TeV Trigger EG2, GammaCut > 6GeV
      cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411792109fe32220000","32c51070a","01031g3200000000","0a53503000000000"); // EG2, new Gamma Energy cut 6 GeV
  } else if(trainConfig == 736)  { //Standard EMCal 13TeV Trigger EG2, GammaCut > 7.5GeV
      cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411792109fe32220000","32c51070a","01031f3200000000","0a53503000000000"); // EG2, new Gamma Energy cut 7.5 GeV

      //Gamma Energy Cuts EG1
  } else if(trainConfig == 737)  { //Standard EMCal 13TeV Trigger EG1, GammaCut > 8GeV
      cuts.AddCutHeavyMesonPCMCalo("0008d113","0dm00009f9730000dge0404000","411792109fe32220000","32c51070a","01031l3200000000","0a53503000000000"); // EG1, new Gamma Energy cut 8 GeV
  } else if(trainConfig == 738)  { //Standard EMCal 13TeV Trigger EG1, GammaCut > 9GeV
      cuts.AddCutHeavyMesonPCMCalo("0008d113","0dm00009f9730000dge0404000","411792109fe32220000","32c51070a","01031m3200000000","0a53503000000000"); // EG1, new Gamma Energy cut 9 GeV
  } else if(trainConfig == 739)  { //Standard EMCal 13TeV Trigger EG1, GammaCut > 10GeV
      cuts.AddCutHeavyMesonPCMCalo("0008d113","0dm00009f9730000dge0404000","411792109fe32220000","32c51070a","01031h3200000000","0a53503000000000"); // EG1, new Gamma Energy cut 10 GeV

      //no shared TPC clusters, Shared cluster Fraction =0
  } else if(trainConfig == 740)  { //Standard EMCal 13TeV, Shared cluster Fraction =0
      cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","411792106fe32220000","32e51070a","0103103200000000","0a53503000000000"); // INT7
  } else if(trainConfig == 741)  { //Standard EMCal 13TeV Trigger EG2, GammaCut > 6GeV, Shared cluster Fraction =0
      cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411792106fe32220000","32e51070a","01031g3200000000","0a53503000000000"); // EG2, new Gamma Energy cut 6 GeV
  } else if(trainConfig == 742)  { //Standard EMCal 13TeV Trigger EG1, GammaCut > 10GeV, Shared cluster Fraction =0
      cuts.AddCutHeavyMesonPCMCalo("0008d113","0dm00009f9730000dge0404000","411792106fe32220000","32e51070a","01031h3200000000","0a53503000000000"); // EG1, new Gamma Energy cut 10 GeV

      //no shared TPC clusters, Shared cluster Fraction <=0.4
  } else if(trainConfig == 743)  { //Standard EMCal 13TeV, Shared cluster Fraction <=0.4
      cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","411792106fe32220000","32f51070a","0103103200000000","0a53503000000000"); // INT7
  } else if(trainConfig == 744)  { //Standard EMCal 13TeV Trigger EG2, GammaCut > 6GeV, Shared cluster Fraction <=0.4
      cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411792106fe32220000","32f51070a","01031g3200000000","0a53503000000000"); // EG2, new Gamma Energy cut 6 GeV
  } else if(trainConfig == 745)  { //Standard EMCal 13TeV Trigger EG1, GammaCut > 10GeV, Shared cluster Fraction <=0.4
      cuts.AddCutHeavyMesonPCMCalo("0008d113","0dm00009f9730000dge0404000","411792106fe32220000","32f51070a","01031h3200000000","0a53503000000000"); // EG1, new Gamma Energy cut 10 GeV

      // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      // EMC pp 13 TeV Fitting, Cut Variations very low EMCAL Pt
      // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  } else if(trainConfig == 900)  { //Standard EMCal 13TeV
      cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","411792109fe32220000","32c51070a","0103103200000000","0453503000000000"); // INT7
      //-----
      //Calo Variations
      //-----
  } else if(trainConfig == 904)  { //EDC 13TeV MB, NCell: 0 (no cut up to 4 GeV ) -> 411792109fe30220000
      cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","411792109fe30220000","32c51070a","0103103200000000","0453503000000000"); // INT7
  } else if(trainConfig == 907)  { //EDC 13TeV MB, NCell: v (NCell Cut 2, but with probability in MC to let clusters through ) -> 411792109fe3n220000
      cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","411792109fe3n220000","32c51070a","0103103200000000","0453503000000000"); // INT7

      //NL Changes
      //NCell 2, NL 21
  } else if(trainConfig == 960)  { //Standard EMCal 13TeV, NCell 2, NL 21
      cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","411792109fe32220000","32c51070a","0103103200000000","0453503000000000"); // INT7
      //NCell 2, NL 01
  } else if(trainConfig == 963)  { //Standard EMCal 13TeV, NCell 2, NL 01
      cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe32220000","32c51070a","0103103200000000","0453503000000000"); // INT7
      //NCell 0, NL 21
  } else if(trainConfig == 970)  { //Standard EMCal 13TeV, NCell 0, NL 21
      cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","411792109fe30220000","32c51070a","0103103200000000","0453503000000000"); // INT7
      //NCell 0, NL 01
  } else if(trainConfig == 973)  { //Standard EMCal 13TeV, NCell 0, NL 01
      cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe30220000","32c51070a","0103103200000000","0453503000000000"); // INT7
      //NCell n (2cell + NCell effi), NL 01
  } else if(trainConfig == 983)  { //Standard EMCal 13TeV, NCell n (NCell >=2 + NCell effi), NL 01
      cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe3n220000","32c51070a","0103103200000000","0453503000000000"); // INT7
      //NCell n (2cell + NCell effi), NL 01
  } else if(trainConfig == 986)  { //Standard EMCal 13TeV, NCell n (NCell >=2 + NCell effi), NL 01 + PCM photon smearing
      cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe3n220000","32c51070a","0103103200b00000","0453503000000000"); // INT7
      //-----
      //Trigger EG2: Calo Variations
      //-----
  } else if(trainConfig == 1004)  { //EDC 13TeV MB, NCell: 0 (no cut up to 4 GeV ) -> 411792109fe30220000
      cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411792109fe30220000","32c51070a","0103103200000000","0453503000000000"); // INT7
  } else if(trainConfig == 1007)  { //EDC 13TeV MB, NCell: v (NCell Cut 2, but with probability in MC to let clusters through ) -> 411792109fe3n220000
      cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411792109fe3n220000","32c51070a","0103103200000000","0453503000000000"); // INT7

      //NL Changes
      //NCell 2, NL 21
  } else if(trainConfig == 1060)  { //EDC 13TeV MB, NCell: 2
    cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411792109fe32220000","32c51070a","0103103200000000","0453503000000000"); // INT7
  } else if(trainConfig == 1061)  { //EDC 13TeV MB, min/max pt cut s (no min pt, max pt  20gev)
    cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411792109fe32220000","32c51070a","01031s3200000000","0453503000000000"); // INT7
  } else if(trainConfig == 1062)  { //EDC 13TeV MB, min/max pt cut v (no min pt, max pt  25gev)
    cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411792109fe32220000","32c51070a","01031v3200000000","0453503000000000"); // INT7
    //NCell 0, NL 01 (new std)
  } else if(trainConfig == 1063)  { //EDC 13TeV MB, NCell: 0 (no cut up to 4 GeV ) -> 411792109fe32220000
    cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411790109fe32220000","32c51070a","0103103200000000","0453503000000000"); // INT7
  } else if(trainConfig == 1064)  { //EDC 13TeV MB, min/max pt cut s (no min pt, max pt  20gev)
    cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411790109fe32220000","32c51070a","01031s3200000000","0453503000000000"); // INT7
  } else if(trainConfig == 1065)  { //EDC 13TeV MB, min/max pt cut v (no min pt, max pt  25gev)
    cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411790109fe32220000","32c51070a","01031v3200000000","0453503000000000"); // INT7

    //NCell 0, NL 21
  } else if(trainConfig == 1070)  { //EDC 13TeV MB, NCell: 0 (no cut up to 4 GeV ) -> 411792109fe30220000
    cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411792109fe30220000","32c51070a","0103103200000000","0453503000000000"); // INT7
  } else if(trainConfig == 1071)  { //EDC 13TeV MB, min/max pt cut s (no min pt, max pt  20gev)
    cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411792109fe30220000","32c51070a","01031s3200000000","0453503000000000"); // INT7
  } else if(trainConfig == 1072)  { //EDC 13TeV MB, min/max pt cut v (no min pt, max pt  25gev)
    cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411792109fe30220000","32c51070a","01031v3200000000","0453503000000000"); // INT7
  //NCell 0, NL 01 (new std)
  } else if(trainConfig == 1073)  { //EDC 13TeV MB, NCell: 2
    cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411790109fe30220000","32c51070a","0103103200000000","0453503000000000"); // INT7
  } else if(trainConfig == 1074)  { //EDC 13TeV MB, min/max pt cut s (no min pt, max pt  20gev)
    cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411790109fe30220000","32c51070a","01031s3200000000","0453503000000000"); // INT7
  } else if(trainConfig == 1075)  { //EDC 13TeV MB, min/max pt cut v (no min pt, max pt  25gev)
    cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411790109fe30220000","32c51070a","01031v3200000000","0453503000000000"); // INT7

    //NCell n (2cell + NCell effi), NL 01
  } else if(trainConfig == 1083)  { //Standard EMCal 13TeV, NCell n (NCell >=2 + NCell effi), NL 01
    cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411790109fe3n220000","32c51070a","0103103200000000","0453503000000000"); // EG2
    //NCell n (2cell + NCell effi), NL 01
  } else if(trainConfig == 1086)  { //Standard EMCal 13TeV, NCell n (NCell >=2 + NCell effi), NL 01 + PCM photon smearing
    cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411790109fe3n220000","32c51070a","0103103200b00000","0453503000000000"); // EG2
      //-----
      //Trigger EG1: Calo Variations
      //-----
  } else if(trainConfig == 1104)  { //EDC 13TeV MB, NCell: 0 (no cut up to 4 GeV ) -> 411792109fe30220000
    cuts.AddCutHeavyMesonPCMCalo("0008d113","0dm00009f9730000dge0404000","411792109fe30220000","32c51070a","0103103200000000","0453503000000000"); // INT7
  } else if(trainConfig == 1107)  { //EDC 13TeV MB, NCell: v (NCell Cut 2, but with probability in MC to let clusters through ) -> 411792109fe3n220000
    cuts.AddCutHeavyMesonPCMCalo("0008d113","0dm00009f9730000dge0404000","411792109fe3n220000","32c51070a","0103103200000000","0453503000000000"); // INT7

    //NL Changes
    //NCell 2, NL 21
  } else if(trainConfig == 1160)  { //EDC 13TeV MB, NCell: 2
    cuts.AddCutHeavyMesonPCMCalo("0008d113","0dm00009f9730000dge0404000","411792109fe32220000","32c51070a","0103103200000000","0453503000000000"); // INT7
  } else if(trainConfig == 1161)  { //EDC 13TeV MB, min/max pt cut s (no min pt, max pt  20gev)
    cuts.AddCutHeavyMesonPCMCalo("0008d113","0dm00009f9730000dge0404000","411792109fe32220000","32c51070a","01031s3200000000","0453503000000000"); // INT7
  } else if(trainConfig == 1162)  { //EDC 13TeV MB, min/max pt cut v (no min pt, max pt  25gev)
    cuts.AddCutHeavyMesonPCMCalo("0008d113","0dm00009f9730000dge0404000","411792109fe32220000","32c51070a","01031v3200000000","0453503000000000"); // INT7
    //NCell 2, NL 01 (new std)
  } else if(trainConfig == 1163)  { //EDC 13TeV MB, NCell: 2
    cuts.AddCutHeavyMesonPCMCalo("0008d113","0dm00009f9730000dge0404000","411790109fe32220000","32c51070a","0103103200000000","0453503000000000"); // INT7
  } else if(trainConfig == 1164)  { //EDC 13TeV MB, min/max pt cut s (no min pt, max pt  20gev)
    cuts.AddCutHeavyMesonPCMCalo("0008d113","0dm00009f9730000dge0404000","411790109fe32220000","32c51070a","01031s3200000000","0453503000000000"); // INT7
  } else if(trainConfig == 1165)  { //EDC 13TeV MB, min/max pt cut v (no min pt, max pt  25gev)
    cuts.AddCutHeavyMesonPCMCalo("0008d113","0dm00009f9730000dge0404000","411790109fe32220000","32c51070a","01031v3200000000","0453503000000000"); // INT7

    //NCell 0, NL 21
  } else if(trainConfig == 1170)  { //EDC 13TeV MB, NCell: 0 (no cut up to 4 GeV ) -> 411792109fe30220000
    cuts.AddCutHeavyMesonPCMCalo("0008d113","0dm00009f9730000dge0404000","411792109fe30220000","32c51070a","0103103200000000","0453503000000000"); // INT7
  } else if(trainConfig == 1171)  { //EDC 13TeV MB, min/max pt cut s (no min pt, max pt  20gev)
    cuts.AddCutHeavyMesonPCMCalo("0008d113","0dm00009f9730000dge0404000","411792109fe30220000","32c51070a","01031s3200000000","0453503000000000"); // INT7
  } else if(trainConfig == 1172)  { //EDC 13TeV MB, min/max pt cut v (no min pt, max pt  25gev)
    cuts.AddCutHeavyMesonPCMCalo("0008d113","0dm00009f9730000dge0404000","411792109fe30220000","32c51070a","01031v3200000000","0453503000000000"); // INT7
    //NCell 0, NL 01
  } else if(trainConfig == 1173)  { //EDC 13TeV MB, NCell: 0 (no cut up to 4 GeV ) -> 411792109fe30220000
    cuts.AddCutHeavyMesonPCMCalo("0008d113","0dm00009f9730000dge0404000","411790109fe30220000","32c51070a","0103103200000000","0453503000000000"); // INT7
  } else if(trainConfig == 1174)  { //EDC 13TeV MB, min/max pt cut s (no min pt, max pt  20gev)
    cuts.AddCutHeavyMesonPCMCalo("0008d113","0dm00009f9730000dge0404000","411790109fe30220000","32c51070a","01031s3200000000","0453503000000000"); // INT7
  } else if(trainConfig == 1175)  { //EDC 13TeV MB, min/max pt cut v (no min pt, max pt  25gev)
    cuts.AddCutHeavyMesonPCMCalo("0008d113","0dm00009f9730000dge0404000","411790109fe30220000","32c51070a","01031v3200000000","0453503000000000"); // INT7

    //NCell n (2cell + NCell effi), NL 01
  } else if(trainConfig == 1183)  { //Standard EMCal 13TeV, NCell n (NCell >=2 + NCell effi), NL 01
    cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411790109fe3n220000","32c51070a","0103103200000000","0453503000000000"); // EG1
    //NCell n (2cell + NCell effi), NL 01
  } else if(trainConfig == 1186)  { //Standard EMCal 13TeV, NCell n (NCell >=2 + NCell effi), NL 01 + PCM photon smearing
    cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411790109fe3n220000","32c51070a","0103103200b00000","0453503000000000"); // EG1

    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // EMC pp 13 TeV Fitting, Systematics
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //To do: Adjust MesonCut string for pi0 for smearing and variation of parametrization
    //Standard Cuts of Pi0 Analysis: ("00010113","0dm00009f9730000dge0404000","411790109fe30220000","0r631031000000d0")
    //MesonCut r63==Background->ignored, d==OpeningAngle for Background->ignored =>0453503000000000
  } else if(trainConfig == 2000)  { //EMCal + DCal INT7 Standard
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe30220000","32c51070a","0103103200000000","0453503000000000"); // INT7 Standard
    //-----
    //INT7: Event Variations
    //-----
    //Std: 00010113
  } else if(trainConfig == 2001)  { //EMCal + DCal INT7, Event cut var. Remove Pileup, Std 1-> True
    //                            00010113
    //                                 |
    cuts.AddCutHeavyMesonPCMCalo("00010013","0dm00009f9730000dge0404000","411790109fe30220000","32c51070a","0103103200000000","0453503000000000"); // INT7 Pileup not removed
    //-----
    //INT7: Calo Variations
    //-----
    //Std: 411790109fe30220000
  } else if(trainConfig == 2101)  { //EMCal + DCal INT7, Calo cut var. NonLins, Std 01
    //                                                                    411790109fe30220000
    //                                                                         ||
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","411799609fe30220000","32c51070a","0103103200000000","0453503000000000"); // INT7 no FT applied
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","411799709fe30220000","32c51070a","0103103200000000","0453503000000000"); // INT7 EMC fine tuning applied
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","411799809fe30220000","32c51070a","0103103200000000","0453503000000000"); // INT7 PCM-EMC fine tuning applied
  } else if(trainConfig == 2102)  { //EMCal + DCal INT7, Calo cut var. time, Std 9 -> -20+25
    //                                                                    411790109fe30220000
    //                                                                            |
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","411790105fe30220000","32c51070a","0103103200000000","0453503000000000"); // INT7 time -50+50
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","411790106fe30220000","32c51070a","0103103200000000","0453503000000000"); // INT7 time -30+35
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","411790108fe30220000","32c51070a","0103103200000000","0453503000000000"); // INT7 time -20+30
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","41179010afe30220000","32c51070a","0103103200000000","0453503000000000"); // INT7 time -12.5+13
  } else if(trainConfig == 2103)  { //EMCal + DCal INT7, Calo cut var. energy, Std 3 -> 0.7 GeV
    //                                                                    411790109fe30220000
    //                                                                               |
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe20220000","32c51070a","0103103200000000","0453503000000000"); // INT7 energy 0.6 GeV
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe40220000","32c51070a","0103103200000000","0453503000000000"); // INT7 energy 0.8 GeV
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe50220000","32c51070a","0103103200000000","0453503000000000"); // INT7 energy 0.9 GeV // only meaningfull at higher pTs
  } else if(trainConfig == 2104)  { //EMCal + DCal INT7, Calo cut var. NCell, Std 0 -> Turned Off until 4GeV; then min 2 Cells
    //                                                                    411790109fe30220000
    //                                                                                |
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe32220000","32c51070a","0103103200000000","0453503000000000"); // INT7 NCells 2
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe3n220000","32c51070a","0103103200000000","0453503000000000"); // INT7 NCells 2 var (PCM-EMCal tagging corr)
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe3r220000","32c51070a","0103103200000000","0453503000000000"); // INT7 NCells 2 var (EMCal tagging corr)
  } else if(trainConfig == 2105)  { //EMCal + DCal INT7, Calo cut var. max M02, 2 -> INT7 M02 0.7
    //                                                                    411790109fe30220000
    //                                                                                  |
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe30230000","32c51070a","0103103200000000","0453503000000000"); // INT7 M02 0.5
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe30210000","32c51070a","0103103200000000","0453503000000000"); // INT7 M02 1.0
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe302k0000","32c51070a","0103103200000000","0453503000000000"); // INT7 M02 E dep
  } else if(trainConfig == 2106)  { //EMCal + DCal INT7, Calo cut var. TM, Std f
    //                                                                    411790109fe30220000
    //                                                                             |
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","411790109ee30220000","32c51070a","0103103200000000","0453503000000000"); // INT7 TM var
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","411790109ge30220000","32c51070a","0103103200000000","0453503000000000"); // INT7 TM var
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","411790109he30220000","32c51070a","0103103200000000","0453503000000000"); // INT7 TM var
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","4117901097e30220000","32c51070a","0103103200000000","0453503000000000"); // INT7 TM var
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","4117901090e30220000","32c51070a","0103103200000000","0453503000000000"); // INT7 TM var, no TM
  } else if(trainConfig == 2107)  { //EMCal + DCal INT7, Calo cut var. Exotics, Std e, active F+ < 0.97
    //                                                                    411790109fe30220000
    //                                                                              |
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","411790109f030220000","32c51070a","0103103200000000","0453503000000000"); // INT7 no exotics cut
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fb30220000","32c51070a","0103103200000000","0453503000000000"); // INT7 F+ < 0.95
    //-----
    //INT7: Primary Pion / Charged Pion (Pi+ Pi-) Variations
    //-----
    //Std: 32c51070a
  } else if(trainConfig == 2201)  { //EMCal + DCal INT7, Ch.Pi cut var. Ch.Pi ITS Requirement, Std 2 -> first or second SPD cluster required
    //                                                                                          32c51070a
    //                                                                                           |
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe32220000","30c51070a","0103103200000000","0453503000000000"); // INT7, Ch.Pi ITS, no SPD cluster requirement
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe32220000","31c51070a","0103103200000000","0453503000000000"); // INT7, Ch.Pi ITS, first SPD cluster required
  } else if(trainConfig == 2202)  { //EMCal + DCal INT7, Ch.Pi cut var. Ch.Pi Cls TPC, Std c -> MinClsTPC 80. + Refit
    //                                                                                          32c51070a
    //                                                                                            |
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe32220000","32251070a","0103103200000000","0453503000000000"); // INT7, Ch.Pi, MinClsTPC 80. + no Refit
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe32220000","32g51070a","0103103200000000","0453503000000000"); // INT7, Ch.Pi, MinClsTPC 70. + Refit
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe32220000","32h51070a","0103103200000000","0453503000000000"); // INT7, Ch.Pi, MinClsTPC 100. + Refit
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe32220000","32i51070a","0103103200000000","0453503000000000"); // INT7, Ch.Pi, MinClsTPC 80. + Refit, 35% find. cls
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe32220000","32j51070a","0103103200000000","0453503000000000"); // INT7, Ch.Pi, MinClsTPC 80. + Refit, 60% find. cls
  } else if(trainConfig == 2203)  { //EMCal + DCal INT7, Ch.Pi cut var. Ch.Pi pT, Std 1 -> pt>0.1
    //                                                                                          32c51070a
    //                                                                                              |
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe32220000","32c50070a","0103103200000000","0453503000000000"); // INT7, Ch.Pi pt>0.075
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe32220000","32c52070a","0103103200000000","0453503000000000"); // INT7, Ch.Pi pt>0.125
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe32220000","32c53070a","0103103200000000","0453503000000000"); // INT7, Ch.Pi pt>0.15
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe32220000","32c54070a","0103103200000000","0453503000000000"); // INT7, Ch.Pi pt>0.4
  } else if(trainConfig == 2204)  { //EMCal + DCal INT7, Ch.Pi cut var. Ch.Pi TPC dEdx, Std 7 -> -3,3
    //                                                                                          32c51070a
    //                                                                                                |
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe32220000","32c51050a","0103103200000000","0453503000000000"); // INT7, Ch.Pi -4,4
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe32220000","32c51080a","0103103200000000","0453503000000000"); // INT7, Ch.Pi -2,3
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe32220000","32c51020a","0103103200000000","0453503000000000"); // INT7, Ch.Pi -6,7
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe32220000","32c51030a","0103103200000000","0453503000000000"); // INT7, Ch.Pi -5,5
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe32220000","32c51040a","0103103200000000","0453503000000000"); // INT7, Ch.Pi -4,5
  } else if(trainConfig == 2205)  { //EMCal + DCal INT7, Ch.Pi cut var. Ch.Pi Mass, Std a -> Ch.Pi<850MeV
    //                                                                                          32c51070a
    //                                                                                                  |
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe32220000","32c51070f","0103103200000000","0453503000000000"); // INT7, Ch.Pi<650MeV
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe32220000","32c51070c","0103103200000000","0453503000000000"); // INT7, Ch.Pi<850MeV, Neu.Pi<1000MeV
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe32220000","32c51070h","0103103200000000","0453503000000000"); // INT7, Ch.Pi<650MeV, Neu.Pi<1000MeV
    //-----
    //INT7: Neutral Meson (Pi0) Cut Variations
    //-----
    //Std: 0103103200000000
  } else if(trainConfig == 2302)  { //EMCal + DCal INT7, N.Pi cut var. rapidity, Std 1 -> -0.8, 0.8
    //                                                                                                      0103103200000000
    //                                                                                                          |
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe30220000","32c51070a","0103503200000000","0453503000000000"); // INT7, N.Pi rap. -0.85, 0.85
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe30220000","32c51070a","0103603200000000","0453503000000000"); // INT7, N.Pi rap. -0.75, 0.75
  } else if(trainConfig == 2304)  { //EMCal + DCal INT7, N.Pi cut var. alpha, Std 3 -> 0.0-1.0
    //                                                                                                      0103103200000000
    //                                                                                                            |
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe30220000","32c51070a","0103105100000000","0453503000000000"); // INT7 alpha 0-0.75
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe30220000","32c51070a","0103108100000000","0453503000000000"); // INT7 alpha 0-0.6
  } else if(trainConfig == 2305)  { //EMCal + DCal INT7, N.Pi cut var. Selection Window, Std 1 -> 2 sigma
    //                                                                                                      0103103200000000
    //                                                                                                             |
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe30220000","32c51070a","0103103u00000000","0453503000000000"); // INT7, 1 sigma, NEEDS CHANGE AS THIS IS EMCEMC
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe30220000","32c51070a","0103103v00000000","0453503000000000"); // INT7, 3 sigma, NEEDS CHANGE AS THIS IS EMCEMC
  } else if(trainConfig == 2306)  { //EMCal + DCal INT7, N.Pi cut var. open. angle, Std 0 -> off
    //                                                                                                      0103103200000000
    //                                                                                                                    |
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe30220000","32c51070a","01031031000000d0","0453503000000000"); // INT7 Op. Ang. var 1 cell dist + 0.017
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe30220000","32c51070a","01031031000000b0","0453503000000000"); // INT7 Op. Ang. var 1 cell dist + 0.0152
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe30220000","32c51070a","01031031000000g0","0453503000000000"); // INT7 Op. Ang. var 1 cell dist + 0.0202
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe30220000","32c51070a","01031031000000a0","0453503000000000"); // INT7 Op. Ang. var 1 cell dist + 0

    //-----
    //INT7: Omega Meson Cut Variations
    //-----
    //Std: 0453503000000000
  } else if(trainConfig == 2401)  { //EMCal + DCal INT7, Omega cut var. Background Scheme, Std 4 -> off
    //                                                                                                                         0453503000000000
    //                                                                                                                          |
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe30220000","32c51070a","0103103200000000","0153103000000000"); // INT7, Om Event Mixing
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe30220000","32c51070a","0103103200000000","0a53503000000000"); // INT7, Om LikeSignMixing
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe30220000","32c51070a","0103103200000000","0d53503000000000"); // INT7, Om SideBandMixing
  } else if(trainConfig == 2402)  { //EMCal + DCal INT7, Omega cut var. rapidity, Std 5 -> -0.85, 0.85
    //                                                                                                                         0453503000000000
    //                                                                                                                             |
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe30220000","32c51070a","0103103200000000","0453103000000000"); // INT7, Om rap. -0.8, 0.8
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe30220000","32c51070a","0103103200000000","0453603000000000"); // INT7, Om rap. -0.75, 0.75
  } else if(trainConfig == 2404)  { //EMCal + DCal INT7, Omega cut var. alpha, Std 3 -> 0.0-1.0
    //                                                                                                                         0453503000000000
    //                                                                                                                               |
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe30220000","32c51070a","0103103200000000","0453505000000000"); // INT7 alpha 0-0.75
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe30220000","32c51070a","0103103200000000","0453508000000000"); // INT7 alpha 0-0.6
  } else if(trainConfig == 2410)  { //EMCal + DCal INT7, Omega cut var. Background Scheme single cfg, Std 4 -> off, EventMixing
    //                                                                                                                         0453503000000000
    //                                                                                                                          |
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe30220000","32c51070a","0103103200000000","0153503000000000"); // INT7, Om Event Mixing
  } else if(trainConfig == 2411)  { //EMCal + DCal INT7, Omega cut var. Background Scheme single cfg, Std 4 -> off, LikeSignMixing
    //                                                                                                                         0453503000000000
    //                                                                                                                          |
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe30220000","32c51070a","0103103200000000","0a53503000000000"); // INT7, Om LikeSignMixing
  } else if(trainConfig == 2412)  { //EMCal + DCal INT7, Omega cut var. Background Scheme single cfg, Std 4 -> off, SideBandMixing
    //                                                                                                                         0453503000000000
    //                                                                                                                          |
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe30220000","32c51070a","0103103200000000","0d53503000000000"); // INT7, Om SideBandMixing

    //-----
    //INT7: PCM Conversion Cut
    //-----
    //Std: 0dm00009f9730000dge0404000
  } else if (trainConfig == 2501) {   // min pT variations
    //                                        0dm00009f9730000dge0404000
    //                                              |
    cuts.AddCutHeavyMesonPCMCalo("00010113", "0dm00069f9730000dge0404000", "411790109fe30220000", "32c51070a", "0103103200000000", "0453503000000000"); // eta < 0.8  // remove  55-72 bin, min pT 40 MeV
    cuts.AddCutHeavyMesonPCMCalo("00010113", "0dm00049f9730000dge0404000", "411790109fe30220000", "32c51070a", "0103103200000000", "0453503000000000"); // eta < 0.8  // remove  55-72 bin, min pT 75 MeV
    cuts.AddCutHeavyMesonPCMCalo("00010113", "0dm00019f9730000dge0404000", "411790109fe30220000", "32c51070a", "0103103200000000", "0453503000000000"); // eta < 0.8  // remove  55-72 bin, min pT 100MeV

  } else if (trainConfig == 2502) {   // TPC clusters, cosPA
    //                                        0dm00009f9730000dge0404000
    //                                               |
    cuts.AddCutHeavyMesonPCMCalo("00010113", "0dm00008f9730000dge0404000", "411790109fe30220000", "32c51070a", "0103103200000000", "0453503000000000"); // TPC cluster 35%
    cuts.AddCutHeavyMesonPCMCalo("00010113", "0dm00006f9730000dge0404000", "411790109fe30220000", "32c51070a", "0103103200000000", "0453503000000000"); // TPC cluster 70%
    cuts.AddCutHeavyMesonPCMCalo("00010113", "0dm00009f9730000dge0604000", "411790109fe30220000", "32c51070a", "0103103200000000", "0453503000000000"); // cosPA 0.9
    cuts.AddCutHeavyMesonPCMCalo("00010113", "0dm00009f9730000dge0304000", "411790109fe30220000", "32c51070a", "0103103200000000", "0453503000000000"); // cosPA 0.75

  } else if (trainConfig == 2503) {   // TPC clusters, cosPA
    //                                        0dm00009f9730000dge0404000
    //                                                |
    cuts.AddCutHeavyMesonPCMCalo("00010113", "0dm0000939730000dge0404000", "411790109fe30220000", "32c51070a", "0103103200000000", "0453503000000000"); // nsig electron   -4,5
    cuts.AddCutHeavyMesonPCMCalo("00010113", "0dm0000969730000dge0404000", "411790109fe30220000", "32c51070a", "0103103200000000", "0453503000000000"); // nsig electron -2.5,4
    cuts.AddCutHeavyMesonPCMCalo("00010113", "0dm00009f5730000dge0404000", "411790109fe30220000", "32c51070a", "0103103200000000", "0453503000000000"); // nsig pion 2,-10
    cuts.AddCutHeavyMesonPCMCalo("00010113", "0dm00009f1730000dge0404000", "411790109fe30220000", "32c51070a", "0103103200000000", "0453503000000000"); // nsig pion 0,-10

  } else if (trainConfig == 2504) {
    //                                        0dm00009f9730000dge0404000
    //                                                  |
    cuts.AddCutHeavyMesonPCMCalo("00010113", "0dm00009f9030000dge0404000", "411790109fe30220000", "32c51070a", "0103103200000000", "0453503000000000"); // pion nsig min mom 0.50 GeV/c
    cuts.AddCutHeavyMesonPCMCalo("00010113", "0dm00009f9630000dge0404000", "411790109fe30220000", "32c51070a", "0103103200000000", "0453503000000000"); // pion nsig min mom 0.25 GeV/c
    cuts.AddCutHeavyMesonPCMCalo("00010113", "0dm00009f9760000dge0404000", "411790109fe30220000", "32c51070a", "0103103200000000", "0453503000000000"); // pion nsig max mom 2.00 GeV/c
    cuts.AddCutHeavyMesonPCMCalo("00010113", "0dm00009f9710000dge0404000", "411790109fe30220000", "32c51070a", "0103103200000000", "0453503000000000"); // pion nsig max mom 5.00 GeV/c


  } else if (trainConfig == 2505) {   // chi2 variations
    //                                        0dm00009f9730000dge0404000
    //                                                         |
    cuts.AddCutHeavyMesonPCMCalo("00010113", "0dm00009f9730000d1e0404000", "411790109fe30220000", "32c51070a", "0103103200000000", "0453503000000000"); // chi2 50
    cuts.AddCutHeavyMesonPCMCalo("00010113", "0dm00009f9730000dfe0404000", "411790109fe30220000", "32c51070a", "0103103200000000", "0453503000000000"); // chi2 50 chi2 dep -0.065
    cuts.AddCutHeavyMesonPCMCalo("00010113", "0dm00009f9730000dhe0404000", "411790109fe30220000", "32c51070a", "0103103200000000", "0453503000000000"); // chi2 50 chi2 dep -0.050
    cuts.AddCutHeavyMesonPCMCalo("00010113", "0dm00009f9730000dge0400000", "411790109fe30220000", "32c51070a", "0103103200000000", "0453503000000000"); // remove reject close v0

  } else if (trainConfig == 2506) {   // Psi pair variations
    //                                        0dm00009f9730000dge0404000
    //                                                         ||
    cuts.AddCutHeavyMesonPCMCalo("00010113", "0dm00009f9730000dgd0404000", "411790109fe30220000", "32c51070a", "0103103200000000", "0453503000000000"); // Psi pair 0.15 dep
    cuts.AddCutHeavyMesonPCMCalo("00010113", "0dm00009f9730000dgf0404000", "411790109fe30220000", "32c51070a", "0103103200000000", "0453503000000000"); // Psi pair 0.20 dep
    cuts.AddCutHeavyMesonPCMCalo("00010113", "0dm00009f9730000dgg0404000", "411790109fe30220000", "32c51070a", "0103103200000000", "0453503000000000"); // Psi pair 0.30 dep
    cuts.AddCutHeavyMesonPCMCalo("00010113", "0dm00009227300008250404000", "411790109fe30220000", "32c51070a", "0103103200000000", "0453503000000000"); // old cuts (run1)

  } else if(trainConfig == 2511)  { //EMCal + DCal INT7 Standard 13TeV, Material Budget Studies
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe30220000","32c51070a","0103103200000000","0453503000000000"); // INT7 Standard
  } else if(trainConfig == 2512)  { //EMCal + DCal INT7 Standard 13TeV, Material Budget Studies
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe30220000","32c51070a","0103103200000000","0453503000000000"); // INT7 Standard
  } else if(trainConfig == 2513)  { //EMCal + DCal INT7 Standard 13TeV, Material Budget Studies
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe30220000","32c51070a","0103103200000000","0453503000000000"); // INT7 Standard
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //Standard Cuts of Pi0 Analysis: ("0008e113","0dm00009f9730000dge0404000","411790109fe30220000","0r631031000000d0")
    //MesonCut r63==Background->ignored, d==OpeningAngle for Background->ignored =>0453503000000000
  } else if(trainConfig == 3000)  { //EMCal + DCal EG2 Standard
    cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411790109fe30220000","32c51070a","01031v3200000000","0453503000000000"); // EG2 Standard
    //-----
    //EG2: Event Variations
    //-----
    //Std: 0008e113
  } else if(trainConfig == 3001)  { //EMCal + DCal EG2, Event cut var. Remove Pileup, Std 1-> True
    //                            0008e113
    //                                 |
    cuts.AddCutHeavyMesonPCMCalo("0008e013","0dm00009f9730000dge0404000","411790109fe30220000","32c51070a","01031v3200000000","0453503000000000"); // INT7 Pileup not removed
    //-----
    //EG2: Calo Variations
    //-----
    //Std: 411790109fe30220000
  } else if(trainConfig == 3101)  { //EMCal + DCal EG2, Calo cut var. NonLins, Std 01
    //                                                                    411790109fe30220000
    //                                                                         ||
    cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411799609fe30220000","32c51070a","01031v3200000000","0453503000000000"); // EG2 no FT applied
    cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411799709fe30220000","32c51070a","01031v3200000000","0453503000000000"); // EG2 EMC fine tuning applied
    cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411799809fe30220000","32c51070a","01031v3200000000","0453503000000000"); // EG2 PCM-EMC fine tuning applied
  } else if(trainConfig == 3102)  { //EMCal + DCal EG2, Calo cut var. time, Std 9 -> -20+25
    //                                                                    411790109fe30220000
    //                                                                            |
    cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411790105fe30220000","32c51070a","01031v3200000000","0453503000000000"); // EG2 time -50+50
    cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411790106fe30220000","32c51070a","01031v3200000000","0453503000000000"); // EG2 time -30+35
    cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411790108fe30220000","32c51070a","01031v3200000000","0453503000000000"); // EG2 time -20+30
    cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","41179010afe30220000","32c51070a","01031v3200000000","0453503000000000"); // EG2 time -12.5+13
  } else if(trainConfig == 3103)  { //EMCal + DCal EG2, Calo cut var. energy, Std 3 -> 0.7 GeV
    //                                                                    411790109fe30220000
    //                                                                               |
    cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411790109fe20220000","32c51070a","01031v3200000000","0453503000000000"); // EG2 energy 0.6 GeV
    cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411790109fe40220000","32c51070a","01031v3200000000","0453503000000000"); // EG2 energy 0.8 GeV
    cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411790109fe50220000","32c51070a","01031v3200000000","0453503000000000"); // EG2 energy 0.9 GeV // only meaningfull at higher pTs
  } else if(trainConfig == 3104)  { //EMCal + DCal EG2, Calo cut var. NCell, Std 0 -> Turned Off until 4GeV; then min 2 Cells
    //                                                                    411790109fe30220000
    //                                                                                |
    cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411790109fe32220000","32c51070a","01031v3200000000","0453503000000000"); // EG2 NCells 2
    cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411790109fe3n220000","32c51070a","01031v3200000000","0453503000000000"); // EG2 NCells 2 var (PCM-EMCal tagging corr)
    cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411790109fe3r220000","32c51070a","01031v3200000000","0453503000000000"); // EG2 NCells 2 var (EMCal tagging corr)
  } else if(trainConfig == 3105)  { //EMCal + DCal EG2, Calo cut var. max M02, 2 -> INT7 M02 0.7
    //                                                                    411790109fe30220000
    //                                                                                  |
    cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411790109fe30230000","32c51070a","01031v3200000000","0453503000000000"); // EG2 M02 0.5
    cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411790109fe30210000","32c51070a","01031v3200000000","0453503000000000"); // EG2 M02 1.0
    cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411790109fe302k0000","32c51070a","01031v3200000000","0453503000000000"); // EG2 M02 E dep
  } else if(trainConfig == 3106)  { //EMCal + DCal EG2, Calo cut var. TM, Std f
    //                                                                    411790109fe30220000
    //                                                                             |
    cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411790109ee30220000","32c51070a","01031v3200000000","0453503000000000"); // EG2 TM var
    cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411790109ge30220000","32c51070a","01031v3200000000","0453503000000000"); // EG2 TM var
    cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411790109he30220000","32c51070a","01031v3200000000","0453503000000000"); // EG2 TM var
    cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","4117901097e30220000","32c51070a","01031v3200000000","0453503000000000"); // EG2 TM var
    cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","4117901090e30220000","32c51070a","01031v3200000000","0453503000000000"); // EG2 TM var, no TM
  } else if(trainConfig == 3107)  { //EMCal + DCal EG2, Calo cut var. Exotics, Std e, active F+ < 0.97
    //                                                                    411790109fe30220000
    //                                                                              |
    cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411790109f030220000","32c51070a","01031v3200000000","0453503000000000"); // EG2 no exotics cut
    cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411790109fb30220000","32c51070a","01031v3200000000","0453503000000000"); // EG2 F+ < 0.95
    //-----
    //EG2: Primary Pion / Charged Pion (Pi+ Pi-) Variations
    //-----
    //Std: 32c51070a
  } else if(trainConfig == 3201)  { //EMCal + DCal EG2, Ch.Pi cut var. Ch.Pi ITS Requirement, Std 2 -> first or second SPD cluster required
    //                                                                                          32c51070a
    //                                                                                           |
    cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411790109fe32220000","30c51070a","01031v3200000000","0453503000000000"); // EG2, Ch.Pi ITS, no SPD cluster requirement
    cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411790109fe32220000","31c51070a","01031v3200000000","0453503000000000"); // EG2, Ch.Pi ITS, first SPD cluster required
  } else if(trainConfig == 3202)  { //EMCal + DCal EG2, Ch.Pi cut var. Ch.Pi Cls TPC, Std c -> MinClsTPC 80. + Refit
    //                                                                                          32c51070a
    //                                                                                            |
    cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411790109fe32220000","32251070a","01031v3200000000","0453503000000000"); // EG2, Ch.Pi, MinClsTPC 80. + no Refit
    cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411790109fe32220000","32g51070a","01031v3200000000","0453503000000000"); // EG2, Ch.Pi, MinClsTPC 70. + Refit
    cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411790109fe32220000","32h51070a","01031v3200000000","0453503000000000"); // EG2, Ch.Pi, MinClsTPC 100. + Refit
    cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411790109fe32220000","32i51070a","01031v3200000000","0453503000000000"); // EG2, Ch.Pi, MinClsTPC 80. + Refit, 35% find. cls
    cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411790109fe32220000","32j51070a","01031v3200000000","0453503000000000"); // EG2, Ch.Pi, MinClsTPC 80. + Refit, 60% find. cls
  } else if(trainConfig == 3203)  { //EMCal + DCal EG2, Ch.Pi cut var. Ch.Pi pT, Std 1 -> pt>0.1
    //                                                                                          32c51070a
    //                                                                                              |
    cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411790109fe32220000","32c50070a","01031v3200000000","0453503000000000"); // EG2, Ch.Pi pt>0.075
    cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411790109fe32220000","32c52070a","01031v3200000000","0453503000000000"); // EG2, Ch.Pi pt>0.125
    cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411790109fe32220000","32c53070a","01031v3200000000","0453503000000000"); // EG2, Ch.Pi pt>0.15
    cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411790109fe32220000","32c54070a","01031v3200000000","0453503000000000"); // EG2, Ch.Pi pt>0.4
  } else if(trainConfig == 3204)  { //EMCal + DCal EG2, Ch.Pi cut var. Ch.Pi TPC dEdx, Std 7 -> -3,3
    //                                                                                          32c51070a
    //                                                                                                |
    cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411790109fe32220000","32c51050a","01031v3200000000","0453503000000000"); // EG2, Ch.Pi -4,4
    cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411790109fe32220000","32c51080a","01031v3200000000","0453503000000000"); // EG2, Ch.Pi -2,3
    cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411790109fe32220000","32c51020a","01031v3200000000","0453503000000000"); // EG2, Ch.Pi -6,7
    cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411790109fe32220000","32c51030a","01031v3200000000","0453503000000000"); // EG2, Ch.Pi -5,5
    cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411790109fe32220000","32c51040a","01031v3200000000","0453503000000000"); // EG2, Ch.Pi -4,5
  } else if(trainConfig == 3205)  { //EMCal + DCal EG2, Ch.Pi cut var. Ch.Pi Mass, Std a -> Ch.Pi<850MeV
    //                                                                                          32c51070a
    //                                                                                                  |
    cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411790109fe32220000","32c51070f","01031v3200000000","0453503000000000"); // EG2, Ch.Pi<650MeV
    cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411790109fe32220000","32c51070c","01031v3200000000","0453503000000000"); // EG2, Ch.Pi<850MeV, Neu.Pi<1000MeV
    cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411790109fe32220000","32c51070h","01031v3200000000","0453503000000000"); // EG2, Ch.Pi<650MeV, Neu.Pi<1000MeV

    //-----
    //EG2: Neutral Meson (Pi0) Cut Variations
    //-----
    //Std: 01031v3200000000
  } else if(trainConfig == 3302)  { //EMCal + DCal EG2, N.Pi cut var. rapidity, Std 1 -> -0.8, 0.8
    //                                                                                                      01031v3200000000
    //                                                                                                          |
    cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411790109fe30220000","32c51070a","01035v3200000000","0453503000000000"); // EG2, N.Pi rap. -0.85, 0.85
    cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411790109fe30220000","32c51070a","01036v3200000000","0453503000000000"); // EG2, N.Pi rap. -0.75, 0.75
  } else if(trainConfig == 3303)  { //EMCal + DCal EG2, N.Pi cut var. maxMass, Std v -> 25GeV
    //                                                                                                      01031v3200000000
    //                                                                                                           |
    cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411790109fe30220000","32c51070a","0103103200000000","0453503000000000"); // EG2, N.Pi maxMass off
    cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411790109fe30220000","32c51070a","01031s3200000000","0453503000000000"); // EG2, N.Pi maxMass 20GeV
    cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411790109fe30220000","32c51070a","01031t3200000000","0453503000000000"); // EG2, N.Pi maxMass 12GeV, maxMass 20GeV
    cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411790109fe30220000","32c51070a","01031w3200000000","0453503000000000"); // EG2, N.Pi minMass 12GeV, maxMass 25GeV
  } else if(trainConfig == 3304)  { //EMCal + DCal EG2, N.Pi cut var. alpha, Std 3 -> 0.0-1.0
    //                                                                                                      01031v3200000000
    //                                                                                                            |
    cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411790109fe30220000","32c51070a","01031v5100000000","0453503000000000"); // EG2 alpha 0-0.75
    cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411790109fe30220000","32c51070a","01031v8100000000","0453503000000000"); // EG2 alpha 0-0.6
  } else if(trainConfig == 3305)  { //EMCal + DCal EG2, N.Pi cut var. Selection Window, Std 1 -> 2 sigma
    //                                                                                                      01031v3200000000
    //                                                                                                             |
    cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411790109fe30220000","32c51070a","01031v3u00000000","0453503000000000"); // EG2, 1 sigma, NEEDS CHANGE AS THIS IS EMCEMC
    cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411790109fe30220000","32c51070a","01031v3v00000000","0453503000000000"); // EG2, 3 sigma, NEEDS CHANGE AS THIS IS EMCEMC
  } else if(trainConfig == 3306)  { //EMCal + DCal EG2, N.Pi cut var. open. angle, Std 0 -> off
    //                                                                                                      01031v3200000000
    //                                                                                                                    |
    cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411790109fe30220000","32c51070a","01031v31000000d0","0453503000000000"); // EG2 Op. Ang. var 1 cell dist + 0.017
    cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411790109fe30220000","32c51070a","01031v31000000b0","0453503000000000"); // EG2 Op. Ang. var 1 cell dist + 0.0152
    cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411790109fe30220000","32c51070a","01031v31000000g0","0453503000000000"); // EG2 Op. Ang. var 1 cell dist + 0.0202
    cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411790109fe30220000","32c51070a","01031v31000000a0","0453503000000000"); // EG2 Op. Ang. var 1 cell dist + 0

    //-----
    //EG2: Omega Meson Cut Variations
    //-----
    //Std: 0453503000000000
  } else if(trainConfig == 3401)  { //EMCal + DCal EG2, Omega cut var. Background Scheme, Std 4 -> off
    //                                                                                                                         0453503000000000
    //                                                                                                                          |
    cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411790109fe30220000","32c51070a","01031v3200000000","0153103000000000"); // EG2, Om Event Mixing
    cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411790109fe30220000","32c51070a","01031v3200000000","0a53603000000000"); // EG2, Om LikeSignMixing
    cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411790109fe30220000","32c51070a","01031v3200000000","0d53603000000000"); // EG2, Om SideBandMixing
  } else if(trainConfig == 3402)  { //EMCal + DCal EG2, Omega cut var. rapidity, Std 5 -> -0.85, 0.85
    //                                                                                                                         0453503000000000
    //                                                                                                                             |
    cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411790109fe30220000","32c51070a","01031v3200000000","0453103000000000"); // EG2, Om rap. -0.8, 0.8
    cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411790109fe30220000","32c51070a","01031v3200000000","0453603000000000"); // EG2, Om rap. -0.75, 0.75
  } else if(trainConfig == 3404)  { //EMCal + DCal EG2, Omega cut var. alpha, Std 3 -> 0.0-1.0
    //                                                                                                                         0453503000000000
    //                                                                                                                               |
    cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411790109fe30220000","32c51070a","01031v3200000000","0453505000000000"); // EG2 alpha 0-0.75
    cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411790109fe30220000","32c51070a","01031v3200000000","0453508000000000"); // EG2 alpha 0-0.6
  } else if(trainConfig == 3410)  { //EMCal + DCal EG2, Omega cut var. Background Scheme single cfg, Std 4 -> off, EventMixing
    //                                                                                                                         0453503000000000
    //                                                                                                                          |
    cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411790109fe30220000","32c51070a","01031v3200000000","0153103000000000"); // EG2, Om Event Mixing
  } else if(trainConfig == 3411)  { //EMCal + DCal EG2, Omega cut var. Background Scheme single cfg, Std 4 -> off, LikeSignMixing
    //                                                                                                                         0453503000000000
    //                                                                                                                          |
    cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411790109fe30220000","32c51070a","01031v3200000000","0a53603000000000"); // EG2, Om LikeSignMixing
  } else if(trainConfig == 3412)  { //EMCal + DCal EG2, Omega cut var. Background Scheme single cfg, Std 4 -> off, SideBandMixing
    //                                                                                                                         0453503000000000
    //                                                                                                                          |
    cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411790109fe30220000","32c51070a","01031v3200000000","0d53603000000000"); // EG2, Om SideBandMixing

    //-----
    //INT7: PCM Conversion Cut
    //-----
    //Std: 0dm00009f9730000dge0404000
  } else if (trainConfig == 3501) {   // min pT variations
    //                                        0dm00009f9730000dge0404000
    //                                              |
    cuts.AddCutHeavyMesonPCMCalo("00010113", "0dm00069f9730000dge0404000", "411790109fe30220000", "32c51070a", "01031v3200000000", "0453503000000000"); // eta < 0.8  // remove  55-72 bin, min pT 40 MeV
    cuts.AddCutHeavyMesonPCMCalo("00010113", "0dm00049f9730000dge0404000", "411790109fe30220000", "32c51070a", "01031v3200000000", "0453503000000000"); // eta < 0.8  // remove  55-72 bin, min pT 75 MeV
    cuts.AddCutHeavyMesonPCMCalo("00010113", "0dm00019f9730000dge0404000", "411790109fe30220000", "32c51070a", "01031v3200000000", "0453503000000000"); // eta < 0.8  // remove  55-72 bin, min pT 100MeV

  } else if (trainConfig == 3502) {   // TPC clusters, cosPA
    //                                        0dm00009f9730000dge0404000
    //                                               |
    cuts.AddCutHeavyMesonPCMCalo("00010113", "0dm00008f9730000dge0404000", "411790109fe30220000", "32c51070a", "01031v3200000000", "0453503000000000"); // TPC cluster 35%
    cuts.AddCutHeavyMesonPCMCalo("00010113", "0dm00006f9730000dge0404000", "411790109fe30220000", "32c51070a", "01031v3200000000", "0453503000000000"); // TPC cluster 70%
    cuts.AddCutHeavyMesonPCMCalo("00010113", "0dm00009f9730000dge0604000", "411790109fe30220000", "32c51070a", "01031v3200000000", "0453503000000000"); // cosPA 0.9
    cuts.AddCutHeavyMesonPCMCalo("00010113", "0dm00009f9730000dge0304000", "411790109fe30220000", "32c51070a", "01031v3200000000", "0453503000000000"); // cosPA 0.75

  } else if (trainConfig == 3503) {   // TPC clusters, cosPA
    //                                        0dm00009f9730000dge0404000
    //                                                |
    cuts.AddCutHeavyMesonPCMCalo("00010113", "0dm0000939730000dge0404000", "411790109fe30220000", "32c51070a", "01031v3200000000", "0453503000000000"); // nsig electron   -4,5
    cuts.AddCutHeavyMesonPCMCalo("00010113", "0dm0000969730000dge0404000", "411790109fe30220000", "32c51070a", "01031v3200000000", "0453503000000000"); // nsig electron -2.5,4
    cuts.AddCutHeavyMesonPCMCalo("00010113", "0dm00009f5730000dge0404000", "411790109fe30220000", "32c51070a", "01031v3200000000", "0453503000000000"); // nsig pion 2,-10
    cuts.AddCutHeavyMesonPCMCalo("00010113", "0dm00009f1730000dge0404000", "411790109fe30220000", "32c51070a", "01031v3200000000", "0453503000000000"); // nsig pion 0,-10

  } else if (trainConfig == 3504) {
    //                                        0dm00009f9730000dge0404000
    //                                                  |
    cuts.AddCutHeavyMesonPCMCalo("00010113", "0dm00009f9030000dge0404000", "411790109fe30220000", "32c51070a", "01031v3200000000", "0453503000000000"); // pion nsig min mom 0.50 GeV/c
    cuts.AddCutHeavyMesonPCMCalo("00010113", "0dm00009f9630000dge0404000", "411790109fe30220000", "32c51070a", "01031v3200000000", "0453503000000000"); // pion nsig min mom 0.25 GeV/c
    cuts.AddCutHeavyMesonPCMCalo("00010113", "0dm00009f9760000dge0404000", "411790109fe30220000", "32c51070a", "01031v3200000000", "0453503000000000"); // pion nsig max mom 2.00 GeV/c
    cuts.AddCutHeavyMesonPCMCalo("00010113", "0dm00009f9710000dge0404000", "411790109fe30220000", "32c51070a", "01031v3200000000", "0453503000000000"); // pion nsig max mom 5.00 GeV/c


  } else if (trainConfig == 3505) {   // chi2 variations
    //                                        0dm00009f9730000dge0404000
    //                                                         |
    cuts.AddCutHeavyMesonPCMCalo("00010113", "0dm00009f9730000d1e0404000", "411790109fe30220000", "32c51070a", "01031v3200000000", "0453503000000000"); // chi2 50
    cuts.AddCutHeavyMesonPCMCalo("00010113", "0dm00009f9730000dfe0404000", "411790109fe30220000", "32c51070a", "01031v3200000000", "0453503000000000"); // chi2 50 chi2 dep -0.065
    cuts.AddCutHeavyMesonPCMCalo("00010113", "0dm00009f9730000dhe0404000", "411790109fe30220000", "32c51070a", "01031v3200000000", "0453503000000000"); // chi2 50 chi2 dep -0.050
    cuts.AddCutHeavyMesonPCMCalo("00010113", "0dm00009f9730000dge0400000", "411790109fe30220000", "32c51070a", "01031v3200000000", "0453503000000000"); // remove reject close v0

  } else if (trainConfig == 3506) {   // Psi pair variations
    //                                        0dm00009f9730000dge0404000
    //                                                         ||
    cuts.AddCutHeavyMesonPCMCalo("00010113", "0dm00009f9730000dgd0404000", "411790109fe30220000", "32c51070a", "01031v3200000000", "0453503000000000"); // Psi pair 0.15 dep
    cuts.AddCutHeavyMesonPCMCalo("00010113", "0dm00009f9730000dgf0404000", "411790109fe30220000", "32c51070a", "01031v3200000000", "0453503000000000"); // Psi pair 0.20 dep
    cuts.AddCutHeavyMesonPCMCalo("00010113", "0dm00009f9730000dgg0404000", "411790109fe30220000", "32c51070a", "01031v3200000000", "0453503000000000"); // Psi pair 0.30 dep
    cuts.AddCutHeavyMesonPCMCalo("00010113", "0dm00009227300008250404000", "411790109fe30220000", "32c51070a", "01031v3200000000", "0453503000000000"); // old cuts (run1)
  } else if(trainConfig == 3511)  { //EMCal + DCal EG2  Standard 13TeV, Material Budget Studies
    cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411790109fe30220000","32c51070a","01031v3200000000","0453503000000000"); // EG2 Standard
  } else if(trainConfig == 3512)  { //EMCal + DCal EG2  Standard 13TeV, Material Budget Studies
    cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411790109fe30220000","32c51070a","01031v3200000000","0453503000000000"); // EG2 Standard
  } else if(trainConfig == 3513)  { //EMCal + DCal EG2  Standard 13TeV, Material Budget Studies
    cuts.AddCutHeavyMesonPCMCalo("0008e113","0dm00009f9730000dge0404000","411790109fe30220000","32c51070a","01031v3200000000","0453503000000000"); // EG2 Standard
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //Standard Cuts of Pi0 Analysis: ("0008d113","0dm00009f9730000dge0404000","411790109fe30220000","0r631031000000d0")
    //MesonCut r63==Background->ignored, d==OpeningAngle for Background->ignored =>0453503000000000
  } else if(trainConfig == 4000)  { //EMCal + DCal EG1 Standard
    cuts.AddCutHeavyMesonPCMCalo("0008d113","0dm00009f9730000dge0404000","411790109fe30220000","32c51070a","01031v3200000000","0453503000000000"); // EG1 Standard
    //-----
    //EG1: Event Variations
    //-----
    //Std: 0008d113
  } else if(trainConfig == 4001)  { //EMCal + DCal EG1, Event cut var. Remove Pileup, Std 1-> True
    //                            0008d113
    //                                 |
    cuts.AddCutHeavyMesonPCMCalo("0008d013","0dm00009f9730000dge0404000","411790109fe30220000","32c51070a","01031v3200000000","0453503000000000"); // INT7 Pileup not removed
    //-----
    //EG1: Calo Variations
    //-----
    //Std: 411790109fe30220000
  } else if(trainConfig == 4101)  { //EMCal + DCal EG1, Calo cut var. NonLins, Std 01
    //                                                                    411790109fe30220000
    //                                                                         ||
    cuts.AddCutHeavyMesonPCMCalo("0008d113","0dm00009f9730000dge0404000","411799609fe30220000","32c51070a","01031v3200000000","0453503000000000"); // EG1 no FT applied
    cuts.AddCutHeavyMesonPCMCalo("0008d113","0dm00009f9730000dge0404000","411799709fe30220000","32c51070a","01031v3200000000","0453503000000000"); // EG1 EMC fine tuning applied
    cuts.AddCutHeavyMesonPCMCalo("0008d113","0dm00009f9730000dge0404000","411799809fe30220000","32c51070a","01031v3200000000","0453503000000000"); // EG1 PCM-EMC fine tuning applied
  } else if(trainConfig == 4102)  { //EMCal + DCal EG1, Calo cut var. time, Std 9 -> -20+25
    //                                                                    411790109fe30220000
    //                                                                            |
    cuts.AddCutHeavyMesonPCMCalo("0008d113","0dm00009f9730000dge0404000","411790105fe30220000","32c51070a","01031v3200000000","0453503000000000"); // EG1 time -50+50
    cuts.AddCutHeavyMesonPCMCalo("0008d113","0dm00009f9730000dge0404000","411790106fe30220000","32c51070a","01031v3200000000","0453503000000000"); // EG1 time -30+35
    cuts.AddCutHeavyMesonPCMCalo("0008d113","0dm00009f9730000dge0404000","411790108fe30220000","32c51070a","01031v3200000000","0453503000000000"); // EG1 time -20+30
    cuts.AddCutHeavyMesonPCMCalo("0008d113","0dm00009f9730000dge0404000","41179010afe30220000","32c51070a","01031v3200000000","0453503000000000"); // EG1 time -12.5+13
  } else if(trainConfig == 4103)  { //EMCal + DCal EG1, Calo cut var. energy, Std 3 -> 0.7 GeV
    //                                                                    411790109fe30220000
    //                                                                          ||
    cuts.AddCutHeavyMesonPCMCalo("0008d113","0dm00009f9730000dge0404000","411790109fe20220000","32c51070a","01031v3200000000","0453503000000000"); // EG1 energy 0.6 GeV
    cuts.AddCutHeavyMesonPCMCalo("0008d113","0dm00009f9730000dge0404000","411790109fe40220000","32c51070a","01031v3200000000","0453503000000000"); // EG1 energy 0.8 GeV
    cuts.AddCutHeavyMesonPCMCalo("0008d113","0dm00009f9730000dge0404000","411790109fe50220000","32c51070a","01031v3200000000","0453503000000000"); // EG1 energy 0.9 GeV // only meaningfull at higher pTs
  } else if(trainConfig == 4104)  { //EMCal + DCal EG1, Calo cut var. NCell, Std 0 -> Turned Off until 4GeV; then min 2 Cells
    //                                                                    411790109fe30220000
    //                                                                                |
    cuts.AddCutHeavyMesonPCMCalo("0008d113","0dm00009f9730000dge0404000","411790109fe32220000","32c51070a","01031v3200000000","0453503000000000"); // EG1 NCells 2
    cuts.AddCutHeavyMesonPCMCalo("0008d113","0dm00009f9730000dge0404000","411790109fe3n220000","32c51070a","01031v3200000000","0453503000000000"); // EG1 NCells 2 var (PCM-EMCal tagging corr)
    cuts.AddCutHeavyMesonPCMCalo("0008d113","0dm00009f9730000dge0404000","411790109fe3r220000","32c51070a","01031v3200000000","0453503000000000"); // EG1 NCells 2 var (EMCal tagging corr)
  } else if(trainConfig == 4105)  { //EMCal + DCal EG1, Calo cut var. max M02, 2 -> INT7 M02 0.7
    //                                                                    411790109fe30220000
    //                                                                                  |
    cuts.AddCutHeavyMesonPCMCalo("0008d113","0dm00009f9730000dge0404000","411790109fe30230000","32c51070a","01031v3200000000","0453503000000000"); // EG1 M02 0.5
    cuts.AddCutHeavyMesonPCMCalo("0008d113","0dm00009f9730000dge0404000","411790109fe30210000","32c51070a","01031v3200000000","0453503000000000"); // EG1 M02 1.0
    cuts.AddCutHeavyMesonPCMCalo("0008d113","0dm00009f9730000dge0404000","411790109fe302k0000","32c51070a","01031v3200000000","0453503000000000"); // EG1 M02 E dep
  } else if(trainConfig == 4106)  { //EMCal + DCal EG1, Calo cut var. TM, Std f
    //                                                                    411790109fe30220000
    //                                                                             |
    cuts.AddCutHeavyMesonPCMCalo("0008d113","0dm00009f9730000dge0404000","411790109ee30220000","32c51070a","01031v3200000000","0453503000000000"); // EG1 TM var
    cuts.AddCutHeavyMesonPCMCalo("0008d113","0dm00009f9730000dge0404000","411790109ge30220000","32c51070a","01031v3200000000","0453503000000000"); // EG1 TM var
    cuts.AddCutHeavyMesonPCMCalo("0008d113","0dm00009f9730000dge0404000","411790109he30220000","32c51070a","01031v3200000000","0453503000000000"); // EG1 TM var
    cuts.AddCutHeavyMesonPCMCalo("0008d113","0dm00009f9730000dge0404000","4117901097e30220000","32c51070a","01031v3200000000","0453503000000000"); // EG1 TM var
    cuts.AddCutHeavyMesonPCMCalo("0008d113","0dm00009f9730000dge0404000","4117901090e30220000","32c51070a","01031v3200000000","0453503000000000"); // EG1 TM var, no TM
  } else if(trainConfig == 4107)  { //EMCal + DCal EG1, Calo cut var. Exotics, Std e, active F+ < 0.97
    //                                                                    411790109fe30220000
    //                                                                              |
    cuts.AddCutHeavyMesonPCMCalo("0008d113","0dm00009f9730000dge0404000","411790109f030220000","32c51070a","01031v3200000000","0453503000000000"); // EG1 no exotics cut
    cuts.AddCutHeavyMesonPCMCalo("0008d113","0dm00009f9730000dge0404000","411790109fb30220000","32c51070a","01031v3200000000","0453503000000000"); // EG1 F+ < 0.95
    //-----
    //EG1: Primary Pion / Charged Pion (Pi+ Pi-) Variations
    //-----
    //Std: 32c51070a
  } else if(trainConfig == 4201)  { //EMCal + DCal EG1, Ch.Pi cut var. Ch.Pi ITS Requirement, Std 2 -> first or second SPD cluster required
    //                                                                                          32c51070a
    //                                                                                           |
    cuts.AddCutHeavyMesonPCMCalo("0008d113","0dm00009f9730000dge0404000","411790109fe32220000","30c51070a","01031v3200000000","0453503000000000"); // EG1, Ch.Pi ITS, no SPD cluster requirement
    cuts.AddCutHeavyMesonPCMCalo("0008d113","0dm00009f9730000dge0404000","411790109fe32220000","31c51070a","01031v3200000000","0453503000000000"); // EG1, Ch.Pi ITS, first SPD cluster required
  } else if(trainConfig == 4202)  { //EMCal + DCal EG1, Ch.Pi cut var. Ch.Pi Cls TPC, Std c -> MinClsTPC 80. + Refit
    //                                                                                          32c51070a
    //                                                                                            |
    cuts.AddCutHeavyMesonPCMCalo("0008d113","0dm00009f9730000dge0404000","411790109fe32220000","32251070a","01031v3200000000","0453503000000000"); // EG1, Ch.Pi, MinClsTPC 80. + no Refit
    cuts.AddCutHeavyMesonPCMCalo("0008d113","0dm00009f9730000dge0404000","411790109fe32220000","32g51070a","01031v3200000000","0453503000000000"); // EG1, Ch.Pi, MinClsTPC 70. + Refit
    cuts.AddCutHeavyMesonPCMCalo("0008d113","0dm00009f9730000dge0404000","411790109fe32220000","32h51070a","01031v3200000000","0453503000000000"); // EG1, Ch.Pi, MinClsTPC 100. + Refit
    cuts.AddCutHeavyMesonPCMCalo("0008d113","0dm00009f9730000dge0404000","411790109fe32220000","32i51070a","01031v3200000000","0453503000000000"); // EG1, Ch.Pi, MinClsTPC 80. + Refit, 35% find. cls
    cuts.AddCutHeavyMesonPCMCalo("0008d113","0dm00009f9730000dge0404000","411790109fe32220000","32j51070a","01031v3200000000","0453503000000000"); // EG1, Ch.Pi, MinClsTPC 80. + Refit, 60% find. cls
  } else if(trainConfig == 4203)  { //EMCal + DCal EG1, Ch.Pi cut var. Ch.Pi pT, Std 1 -> pt>0.1
    //                                                                                          32c51070a
    //                                                                                              |
    cuts.AddCutHeavyMesonPCMCalo("0008d113","0dm00009f9730000dge0404000","411790109fe32220000","32c50070a","01031v3200000000","0453503000000000"); // EG1, Ch.Pi pt>0.075
    cuts.AddCutHeavyMesonPCMCalo("0008d113","0dm00009f9730000dge0404000","411790109fe32220000","32c52070a","01031v3200000000","0453503000000000"); // EG1, Ch.Pi pt>0.125
    cuts.AddCutHeavyMesonPCMCalo("0008d113","0dm00009f9730000dge0404000","411790109fe32220000","32c53070a","01031v3200000000","0453503000000000"); // EG1, Ch.Pi pt>0.15
    cuts.AddCutHeavyMesonPCMCalo("0008d113","0dm00009f9730000dge0404000","411790109fe32220000","32c54070a","01031v3200000000","0453503000000000"); // EG1, Ch.Pi pt>0.4
  } else if(trainConfig == 4204)  { //EMCal + DCal EG1, Ch.Pi cut var. Ch.Pi TPC dEdx, Std 7 -> -3,3
    //                                                                                          32c51070a
    //                                                                                                |
    cuts.AddCutHeavyMesonPCMCalo("0008d113","0dm00009f9730000dge0404000","411790109fe32220000","32c51050a","01031v3200000000","0453503000000000"); // EG1, Ch.Pi -4,4
    cuts.AddCutHeavyMesonPCMCalo("0008d113","0dm00009f9730000dge0404000","411790109fe32220000","32c51080a","01031v3200000000","0453503000000000"); // EG1, Ch.Pi -2,3
    cuts.AddCutHeavyMesonPCMCalo("0008d113","0dm00009f9730000dge0404000","411790109fe32220000","32c51020a","01031v3200000000","0453503000000000"); // EG1, Ch.Pi -6,7
    cuts.AddCutHeavyMesonPCMCalo("0008d113","0dm00009f9730000dge0404000","411790109fe32220000","32c51030a","01031v3200000000","0453503000000000"); // EG1, Ch.Pi -5,5
    cuts.AddCutHeavyMesonPCMCalo("0008d113","0dm00009f9730000dge0404000","411790109fe32220000","32c51040a","01031v3200000000","0453503000000000"); // EG1, Ch.Pi -4,5
  } else if(trainConfig == 4205)  { //EMCal + DCal EG1, Ch.Pi cut var. Ch.Pi Mass, Std a -> Ch.Pi<850MeV
    //                                                                                          32c51070a
    //                                                                                                  |
    cuts.AddCutHeavyMesonPCMCalo("0008d113","0dm00009f9730000dge0404000","411790109fe32220000","32c51070f","01031v3200000000","0453503000000000"); // EG1, Ch.Pi<650MeV
    cuts.AddCutHeavyMesonPCMCalo("0008d113","0dm00009f9730000dge0404000","411790109fe32220000","32c51070c","01031v3200000000","0453503000000000"); // EG1, Ch.Pi<850MeV, Neu.Pi<1000MeV
    cuts.AddCutHeavyMesonPCMCalo("0008d113","0dm00009f9730000dge0404000","411790109fe32220000","32c51070h","01031v3200000000","0453503000000000"); // EG1, Ch.Pi<650MeV, Neu.Pi<1000MeV

    //-----
    //EG1: Neutral Meson (Pi0) Cut Variations
    //-----
    //Std: 01031v3200000000
  } else if(trainConfig == 4302)  { //EMCal + DCal EG1, N.Pi cut var. rapidity, Std 1 -> -0.8, 0.8
    //                                                                                                      01031v3200000000
    //                                                                                                          |
    cuts.AddCutHeavyMesonPCMCalo("0008d113","0dm00009f9730000dge0404000","411790109fe30220000","32c51070a","01035v3200000000","0453503000000000"); // EG1, N.Pi rap. -0.85, 0.85
    cuts.AddCutHeavyMesonPCMCalo("0008d113","0dm00009f9730000dge0404000","411790109fe30220000","32c51070a","01036v3200000000","0453503000000000"); // EG1, N.Pi rap. -0.75, 0.75
  } else if(trainConfig == 4303)  { //EMCal + DCal EG1, N.Pi cut var. maxMass, Std v -> 25GeV
    //                                                                                                      01031v3200000000
    //                                                                                                           |
    cuts.AddCutHeavyMesonPCMCalo("0008d113","0dm00009f9730000dge0404000","411790109fe30220000","32c51070a","0103103200000000","0453503000000000"); // EG1, N.Pi maxMass off
    cuts.AddCutHeavyMesonPCMCalo("0008d113","0dm00009f9730000dge0404000","411790109fe30220000","32c51070a","01031s3200000000","0453503000000000"); // EG1, N.Pi maxMass 20GeV
    cuts.AddCutHeavyMesonPCMCalo("0008d113","0dm00009f9730000dge0404000","411790109fe30220000","32c51070a","01031t3200000000","0453503000000000"); // EG1, N.Pi maxMass 12GeV, maxMass 20GeV
    cuts.AddCutHeavyMesonPCMCalo("0008d113","0dm00009f9730000dge0404000","411790109fe30220000","32c51070a","01031w3200000000","0453503000000000"); // EG1, N.Pi minMass 12GeV, maxMass 25GeV
  } else if(trainConfig == 4304)  { //EMCal + DCal EG1, N.Pi cut var. alpha, Std 3 -> 0.0-1.0
    //                                                                                                      01031v3200000000
    //                                                                                                            |
    cuts.AddCutHeavyMesonPCMCalo("0008d113","0dm00009f9730000dge0404000","411790109fe30220000","32c51070a","01031v5100000000","0453503000000000"); // EG1 alpha 0-0.75
    cuts.AddCutHeavyMesonPCMCalo("0008d113","0dm00009f9730000dge0404000","411790109fe30220000","32c51070a","01031v8100000000","0453503000000000"); // EG1 alpha 0-0.6
  } else if(trainConfig == 4305)  { //EMCal + DCal EG1, N.Pi cut var. Selection Window, Std 1 -> 2 sigma
    //                                                                                                      01031v3200000000
    //                                                                                                             |
    cuts.AddCutHeavyMesonPCMCalo("0008d113","0dm00009f9730000dge0404000","411790109fe30220000","32c51070a","01031v3u00000000","0453503000000000"); // EG1, 1 sigma, NEEDS CHANGE AS THIS IS EMCEMC
    cuts.AddCutHeavyMesonPCMCalo("0008d113","0dm00009f9730000dge0404000","411790109fe30220000","32c51070a","01031v3v00000000","0453503000000000"); // EG1, 3 sigma, NEEDS CHANGE AS THIS IS EMCEMC
  } else if(trainConfig == 4306)  { //EMCal + DCal EG1, N.Pi cut var. open. angle, Std 0 -> off
    //                                                                                                      01031v3200000000
    //                                                                                                                    |
    cuts.AddCutHeavyMesonPCMCalo("0008d113","0dm00009f9730000dge0404000","411790109fe30220000","32c51070a","01031v31000000d0","0453503000000000"); // EG1 Op. Ang. var 1 cell dist + 0.017
    cuts.AddCutHeavyMesonPCMCalo("0008d113","0dm00009f9730000dge0404000","411790109fe30220000","32c51070a","01031v31000000b0","0453503000000000"); // EG1 Op. Ang. var 1 cell dist + 0.0152
    cuts.AddCutHeavyMesonPCMCalo("0008d113","0dm00009f9730000dge0404000","411790109fe30220000","32c51070a","01031v31000000g0","0453503000000000"); // EG1 Op. Ang. var 1 cell dist + 0.0202
    cuts.AddCutHeavyMesonPCMCalo("0008d113","0dm00009f9730000dge0404000","411790109fe30220000","32c51070a","01031v31000000a0","0453503000000000"); // EG1 Op. Ang. var 1 cell dist + 0

    //-----
    //EG1: Omega Meson Cut Variations
    //-----
    //Std: 0453503000000000
  } else if(trainConfig == 4401)  { //EMCal + DCal EG1, Omega cut var. Background Scheme, Std 4 -> off
    //                                                                                                                         0453503000000000
    //                                                                                                                          |
    cuts.AddCutHeavyMesonPCMCalo("0008d113","0dm00009f9730000dge0404000","411790109fe30220000","32c51070a","01031v3200000000","0153503000000000"); // EG1, Om Event Mixing
    cuts.AddCutHeavyMesonPCMCalo("0008d113","0dm00009f9730000dge0404000","411790109fe30220000","32c51070a","01031v3200000000","0a53503000000000"); // EG1, Om LikeSignMixing
    cuts.AddCutHeavyMesonPCMCalo("0008d113","0dm00009f9730000dge0404000","411790109fe30220000","32c51070a","01031v3200000000","0d53503000000000"); // EG1, Om SideBandMixing
  } else if(trainConfig == 4402)  { //EMCal + DCal EG1, Omega cut var. rapidity, Std 5 -> -0.85, 0.85
    //                                                                                                                         0453503000000000
    //                                                                                                                             |
    cuts.AddCutHeavyMesonPCMCalo("0008d113","0dm00009f9730000dge0404000","411790109fe30220000","32c51070a","01031v3200000000","0453103000000000"); // EG1, Om rap. -0.8, 0.8
    cuts.AddCutHeavyMesonPCMCalo("0008d113","0dm00009f9730000dge0404000","411790109fe30220000","32c51070a","01031v3200000000","0453603000000000"); // EG1, Om rap. -0.75, 0.75
  } else if(trainConfig == 4404)  { //EMCal + DCal EG1, Omega cut var. alpha, Std 3 -> 0.0-1.0
    //                                                                                                                         0453503000000000
    //                                                                                                                               |
    cuts.AddCutHeavyMesonPCMCalo("0008d113","0dm00009f9730000dge0404000","411790109fe30220000","32c51070a","01031v3200000000","0453505000000000"); // EG1 alpha 0-0.75
    cuts.AddCutHeavyMesonPCMCalo("0008d113","0dm00009f9730000dge0404000","411790109fe30220000","32c51070a","01031v3200000000","0453508000000000"); // EG1 alpha 0-0.6
  } else if(trainConfig == 4410)  { //EMCal + DCal EG1, Omega cut var. Background Scheme single cfg, Std 4 -> off, EventMixing
    //                                                                                                                         0453503000000000
    //                                                                                                                          |
    cuts.AddCutHeavyMesonPCMCalo("0008d113","0dm00009f9730000dge0404000","411790109fe30220000","32c51070a","01031v3200000000","0153503000000000"); // EG1, Om Event Mixing
  } else if(trainConfig == 4411)  { //EMCal + DCal EG1, Omega cut var. Background Scheme single cfg, Std 4 -> off, LikeSignMixing
    //                                                                                                                         0453503000000000
    //                                                                                                                          |
    cuts.AddCutHeavyMesonPCMCalo("0008d113","0dm00009f9730000dge0404000","411790109fe30220000","32c51070a","01031v3200000000","0a53503000000000"); // EG1, Om LikeSignMixing
  } else if(trainConfig == 4412)  { //EMCal + DCal EG1, Omega cut var. Background Scheme single cfg, Std 4 -> off, SideBandMixing
    //                                                                                                                         0453503000000000
    //                                                                                                                          |
    cuts.AddCutHeavyMesonPCMCalo("0008d113","0dm00009f9730000dge0404000","411790109fe30220000","32c51070a","01031v3200000000","0d53503000000000"); // EG1, Om SideBandMixing
    //-----
    //INT7: PCM Conversion Cut
    //-----
    //Std: 0dm00009f9730000dge0404000
  } else if (trainConfig == 4501) {   // min pT variations
    //                                        0dm00009f9730000dge0404000
    //                                              |
    cuts.AddCutHeavyMesonPCMCalo("00010113", "0dm00069f9730000dge0404000", "411790109fe30220000", "32c51070a", "01031v3200000000", "0453503000000000"); // eta < 0.8  // remove  55-72 bin, min pT 40 MeV
    cuts.AddCutHeavyMesonPCMCalo("00010113", "0dm00049f9730000dge0404000", "411790109fe30220000", "32c51070a", "01031v3200000000", "0453503000000000"); // eta < 0.8  // remove  55-72 bin, min pT 75 MeV
    cuts.AddCutHeavyMesonPCMCalo("00010113", "0dm00019f9730000dge0404000", "411790109fe30220000", "32c51070a", "01031v3200000000", "0453503000000000"); // eta < 0.8  // remove  55-72 bin, min pT 100MeV

  } else if (trainConfig == 4502) {   // TPC clusters, cosPA
    //                                        0dm00009f9730000dge0404000
    //                                               |
    cuts.AddCutHeavyMesonPCMCalo("00010113", "0dm00008f9730000dge0404000", "411790109fe30220000", "32c51070a", "01031v3200000000", "0453503000000000"); // TPC cluster 35%
    cuts.AddCutHeavyMesonPCMCalo("00010113", "0dm00006f9730000dge0404000", "411790109fe30220000", "32c51070a", "01031v3200000000", "0453503000000000"); // TPC cluster 70%
    cuts.AddCutHeavyMesonPCMCalo("00010113", "0dm00009f9730000dge0604000", "411790109fe30220000", "32c51070a", "01031v3200000000", "0453503000000000"); // cosPA 0.9
    cuts.AddCutHeavyMesonPCMCalo("00010113", "0dm00009f9730000dge0304000", "411790109fe30220000", "32c51070a", "01031v3200000000", "0453503000000000"); // cosPA 0.75

  } else if (trainConfig == 4503) {   // TPC clusters, cosPA
    //                                        0dm00009f9730000dge0404000
    //                                                |
    cuts.AddCutHeavyMesonPCMCalo("00010113", "0dm0000939730000dge0404000", "411790109fe30220000", "32c51070a", "01031v3200000000", "0453503000000000"); // nsig electron   -4,5
    cuts.AddCutHeavyMesonPCMCalo("00010113", "0dm0000969730000dge0404000", "411790109fe30220000", "32c51070a", "01031v3200000000", "0453503000000000"); // nsig electron -2.5,4
    cuts.AddCutHeavyMesonPCMCalo("00010113", "0dm00009f5730000dge0404000", "411790109fe30220000", "32c51070a", "01031v3200000000", "0453503000000000"); // nsig pion 2,-10
    cuts.AddCutHeavyMesonPCMCalo("00010113", "0dm00009f1730000dge0404000", "411790109fe30220000", "32c51070a", "01031v3200000000", "0453503000000000"); // nsig pion 0,-10

  } else if (trainConfig == 4504) {
    //                                        0dm00009f9730000dge0404000
    //                                                  |
    cuts.AddCutHeavyMesonPCMCalo("00010113", "0dm00009f9030000dge0404000", "411790109fe30220000", "32c51070a", "01031v3200000000", "0453503000000000"); // pion nsig min mom 0.50 GeV/c
    cuts.AddCutHeavyMesonPCMCalo("00010113", "0dm00009f9630000dge0404000", "411790109fe30220000", "32c51070a", "01031v3200000000", "0453503000000000"); // pion nsig min mom 0.25 GeV/c
    cuts.AddCutHeavyMesonPCMCalo("00010113", "0dm00009f9760000dge0404000", "411790109fe30220000", "32c51070a", "01031v3200000000", "0453503000000000"); // pion nsig max mom 2.00 GeV/c
    cuts.AddCutHeavyMesonPCMCalo("00010113", "0dm00009f9710000dge0404000", "411790109fe30220000", "32c51070a", "01031v3200000000", "0453503000000000"); // pion nsig max mom 5.00 GeV/c


  } else if (trainConfig == 4505) {   // chi2 variations
    //                                        0dm00009f9730000dge0404000
    //                                                         |
    cuts.AddCutHeavyMesonPCMCalo("00010113", "0dm00009f9730000d1e0404000", "411790109fe30220000", "32c51070a", "01031v3200000000", "0453503000000000"); // chi2 50
    cuts.AddCutHeavyMesonPCMCalo("00010113", "0dm00009f9730000dfe0404000", "411790109fe30220000", "32c51070a", "01031v3200000000", "0453503000000000"); // chi2 50 chi2 dep -0.065
    cuts.AddCutHeavyMesonPCMCalo("00010113", "0dm00009f9730000dhe0404000", "411790109fe30220000", "32c51070a", "01031v3200000000", "0453503000000000"); // chi2 50 chi2 dep -0.050
    cuts.AddCutHeavyMesonPCMCalo("00010113", "0dm00009f9730000dge0400000", "411790109fe30220000", "32c51070a", "01031v3200000000", "0453503000000000"); // remove reject close v0

  } else if (trainConfig == 4506) {   // Psi pair variations
    //                                        0dm00009f9730000dge0404000
    //                                                         ||
    cuts.AddCutHeavyMesonPCMCalo("00010113", "0dm00009f9730000dgd0404000", "411790109fe30220000", "32c51070a", "01031v3200000000", "0453503000000000"); // Psi pair 0.15 dep
    cuts.AddCutHeavyMesonPCMCalo("00010113", "0dm00009f9730000dgf0404000", "411790109fe30220000", "32c51070a", "01031v3200000000", "0453503000000000"); // Psi pair 0.20 dep
    cuts.AddCutHeavyMesonPCMCalo("00010113", "0dm00009f9730000dgg0404000", "411790109fe30220000", "32c51070a", "01031v3200000000", "0453503000000000"); // Psi pair 0.30 dep
    cuts.AddCutHeavyMesonPCMCalo("00010113", "0dm00009227300008250404000", "411790109fe30220000", "32c51070a", "01031v3200000000", "0453503000000000"); // old cuts (run1)
  } else if(trainConfig == 4511)  { //EMCal + DCal EG1  Standard 13TeV, Material Budget Studies
    cuts.AddCutHeavyMesonPCMCalo("0008d113","0dm00009f9730000dge0404000","411790109fe30220000","32c51070a","01031v3200000000","0453503000000000"); // EG1 Standard
  } else if(trainConfig == 4512)  { //EMCal + DCal EG1  Standard 13TeV, Material Budget Studies
    cuts.AddCutHeavyMesonPCMCalo("0008d113","0dm00009f9730000dge0404000","411790109fe30220000","32c51070a","01031v3200000000","0453503000000000"); // EG1 Standard
  } else if(trainConfig == 4513)  { //EMCal + DCal EG1  Standard 13TeV, Material Budget Studies
    cuts.AddCutHeavyMesonPCMCalo("0008d113","0dm00009f9730000dge0404000","411790109fe30220000","32c51070a","01031v3200000000","0453503000000000"); // EG1 Standard
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // PHOS pp 13 TeV Fitting, Systematics
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //To do: Adjust MesonCut string for pi0 for smearing and variation of parametrization
    //Standard Cuts of Pi0 Analysis: ("00010113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163103100000010")
    //                                                                                                 |    |      |
    //                                                                                              "0103103400000000"
    //MesonCut 063==Background->ignored, 4==Pi0 Mass Window->ignored, 1==OpeningAngle for Background->ignored =>0103103200000000
  } else if(trainConfig == 6000)  { //PHOS INT7 Standard
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa01cc00000","32c51070a","0103103400000000","0453503000000000"); // INT7 Standard
    //-----
    //INT7: Event Variations
    //-----
    //Std: 00010113
  } else if(trainConfig == 6001)  { //PHOS INT7, Event cut var. Remove Pileup, Std 1-> True
    //                            00010113
    //                                 |
    cuts.AddCutHeavyMesonPCMCalo("00010013","0dm00009f9730000dge0404000","24466190sa01cc00000","32c51070a","0103103400000000","0453503000000000"); // INT7 Pileup not removed
    //-----
    //INT7: Calo Variations
    //-----
    //Std: 24466190sa01cc00000
  } else if(trainConfig == 6101)  { //PHOS INT7, Calo cut var. NonLins, Std 01
    //                                                                    24466190sa01cc00000
    //                                                                         ||
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","24466000sa01cc00000","32c51070a","0103103400000000","0453503000000000"); // INT7 NoNL
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","24466110sa01cc00000","32c51070a","0103103400000000","0453503000000000"); // INT7 NL11
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","24466120sa01cc00000","32c51070a","0103103400000000","0453503000000000"); // INT7 NL12
  } else if(trainConfig == 6102)  { //PHOS INT7, Calo cut var. time, Std 9 -> -20+25
    //                                                                    24466190sa01cc00000
    //                                                                            |
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","24466190ra01cc00000","32c51070a","0103103400000000","0453503000000000"); // INT7 r: LowPt from MB, 30ns
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","24466190ta01cc00000","32c51070a","0103103400000000","0453503000000000"); // INT7 t: 25ns
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","24466190ua01cc00000","32c51070a","0103103400000000","0453503000000000"); // INT7 u: 50ns
  } else if(trainConfig == 6103)  { //PHOS INT7, Calo cut var. energy, Std 1 -> 0.6 GeV
    //                                                                    24466190sa01cc00000
    //                                                                               |
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa09cc00000","32c51070a","0103103400000000","0453503000000000"); // INT7 9: 0.1 GeV/c
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa02cc00000","32c51070a","0103103400000000","0453503000000000"); // INT7 2: 0.5 GeV/c
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa03cc00000","32c51070a","0103103400000000","0453503000000000"); // INT7 3: 0.6 GeV/c
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa04cc00000","32c51070a","0103103400000000","0453503000000000"); // INT7 4: 0.7 GeV/c
  } else if(trainConfig == 6104)  { //PHOS INT7, Calo cut var. NCell & M02, Std cc -> min nCells = 2 >1GeV; M02 max=100, min=0.1, part 2
      //                                                                  24466190sa01cc00000
      //                                                                              |||
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa012200000","32c51070a","0103103400000000","0453503000000000"); // INT7 220: min nCells = 2, all E
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa013200000","32c51070a","0103103400000000","0453503000000000"); // INT7 320: min nCells = 3, all E
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa01dc00000","32c51070a","0103103400000000","0453503000000000"); // INT7 dc0: min nCells = 3, E>1GeV; minM02==0.1 off for E<1GeV
  } else if(trainConfig == 6105)  { //PHOS INT7, Calo cut var. NCell & M02, Std cc -> min nCells = 2 >1GeV; M02 max=100, min=0.1, part 2
      //                                                                  24466190sa01cc00000
      //                                                                              |||
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa011000000","32c51070a","0103103400000000","0453503000000000"); // INT7 100: min nCells = 1, minM02 off
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa01cc70000","32c51070a","0103103400000000","0453503000000000"); // INT7 cc7: maxM02 == 1.3
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa01cc80000","32c51070a","0103103400000000","0453503000000000"); // INT7 cc8: maxM02 == 2.5
  } else if(trainConfig == 6106)  { //PHOS INT7, Calo cut var. TM, Std a
    //                                                                    24466190sa01cc00000
    //                                                                             |
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","24466190s001cc00000","32c51070a","0103103400000000","0453503000000000"); // INT7 TM var 0
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","24466190s101cc00000","32c51070a","0103103400000000","0453503000000000"); // INT7 TM var 1
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","24466190s401cc00000","32c51070a","0103103400000000","0453503000000000"); // INT7 TM var 4
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","24466190s501cc00000","32c51070a","0103103400000000","0453503000000000"); // INT7 TM var 5
  } else if(trainConfig == 6000)  { //PHOS INT7 Calo cut var. reconstructed conversion
      //                                                                  24466190sa01cc00000
      //                                                                                  |
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa01cc00100","32c51070a","0103103400000000","0453503000000000"); // INT7 0: 0.02
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa01cc00200","32c51070a","0103103400000000","0453503000000000"); // INT7 1: 0.025
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa01cc00300","32c51070a","0103103400000000","0453503000000000"); // INT7 3: 0.03
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa01cc00400","32c51070a","0103103400000000","0453503000000000"); // INT7 4: 0.035
    //-----
    //INT7: Primary Pion / Charged Pion (Pi+ Pi-) Variations
    //-----
    //Std: 32c51070a
  } else if(trainConfig == 6201)  { //PHOS INT7, Ch.Pi cut var. Ch.Pi ITS Requirement, Std 2 -> first or second SPD cluster required
    //                                                                                          32c51070a
    //                                                                                           |
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa01cc00000","30c51070a","0103103400000000","0453503000000000"); // INT7, Ch.Pi ITS, no SPD cluster requirement
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa01cc00000","31c51070a","0103103400000000","0453503000000000"); // INT7, Ch.Pi ITS, first SPD cluster required
  } else if(trainConfig == 6202)  { //PHOS INT7, Ch.Pi cut var. Ch.Pi Cls TPC, Std c -> MinClsTPC 80. + Refit
    //                                                                                          32c51070a
    //                                                                                            |
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa01cc00000","32251070a","0103103400000000","0453503000000000"); // INT7, Ch.Pi, MinClsTPC 80. + no Refit
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa01cc00000","32g51070a","0103103400000000","0453503000000000"); // INT7, Ch.Pi, MinClsTPC 70. + Refit
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa01cc00000","32h51070a","0103103400000000","0453503000000000"); // INT7, Ch.Pi, MinClsTPC 100. + Refit
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa01cc00000","32i51070a","0103103400000000","0453503000000000"); // INT7, Ch.Pi, MinClsTPC 80. + Refit, 35% find. cls
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa01cc00000","32j51070a","0103103400000000","0453503000000000"); // INT7, Ch.Pi, MinClsTPC 80. + Refit, 60% find. cls
  } else if(trainConfig == 6203)  { //PHOS INT7, Ch.Pi cut var. Ch.Pi pT, Std 1 -> pt>0.1
    //                                                                                          32c51070a
    //                                                                                              |
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa01cc00000","32c50070a","0103103400000000","0453503000000000"); // INT7, Ch.Pi pt>0.075
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa01cc00000","32c52070a","0103103400000000","0453503000000000"); // INT7, Ch.Pi pt>0.125
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa01cc00000","32c53070a","0103103400000000","0453503000000000"); // INT7, Ch.Pi pt>0.15
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa01cc00000","32c54070a","0103103400000000","0453503000000000"); // INT7, Ch.Pi pt>0.4
  } else if(trainConfig == 6204)  { //PHOS INT7, Ch.Pi cut var. Ch.Pi TPC dEdx, Std 7 -> -3,3
    //                                                                                          32c51070a
    //                                                                                                |
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa01cc00000","32c51050a","0103103400000000","0453503000000000"); // INT7, Ch.Pi -4,4
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa01cc00000","32c51080a","0103103400000000","0453503000000000"); // INT7, Ch.Pi -2,3
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa01cc00000","32c51020a","0103103400000000","0453503000000000"); // INT7, Ch.Pi -6,7
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa01cc00000","32c51030a","0103103400000000","0453503000000000"); // INT7, Ch.Pi -5,5
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa01cc00000","32c51040a","0103103400000000","0453503000000000"); // INT7, Ch.Pi -4,5
  } else if(trainConfig == 6205)  { //PHOS INT7, Ch.Pi cut var. Ch.Pi Mass, Std a -> Ch.Pi<850MeV
    //                                                                                          32c51070a
    //                                                                                                  |
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa01cc00000","32c51070f","0103103400000000","0453503000000000"); // INT7, Ch.Pi<650MeV
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa01cc00000","32c51070c","0103103400000000","0453503000000000"); // INT7, Ch.Pi<850MeV, Neu.Pi<1000MeV
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa01cc00000","32c51070h","0103103400000000","0453503000000000"); // INT7, Ch.Pi<650MeV, Neu.Pi<1000MeV
    //-----
    //INT7: Neutral Meson (Pi0) Cut Variations
    //-----
    //Std: 0103103400000000
  } else if(trainConfig == 6302)  { //PHOS INT7, N.Pi cut var. rapidity, Std 1 -> -0.8, 0.8
    //                                                                                                      0103103400000000
    //                                                                                                          |
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa01cc00000","32c51070a","0103503400000000","0453503000000000"); // INT7, N.Pi rap. 5: -0.85, 0.85
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa01cc00000","32c51070a","0103603400000000","0453503000000000"); // INT7, N.Pi rap. 6: -0.75, 0.75
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa01cc00000","32c51070a","0103403400000000","0453503000000000"); // INT7, N.Pi rap. 4: -0.5, 0.5
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa01cc00000","32c51070a","0103803400000000","0453503000000000"); // INT7, N.Pi rap. 8: -0.25, 0.25
  } else if(trainConfig == 6304)  { //PHOS INT7, N.Pi cut var. alpha, Std 3 -> 0.0-1.0
    //                                                                                                      0103103400000000
    //                                                                                                            |
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa01cc00000","32c51070a","0103105400000000","0453503000000000"); // INT7 alpha 0-0.75
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa01cc00000","32c51070a","0103108400000000","0453503000000000"); // INT7 alpha 0-0.6
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa01cc00000","32c51070a","0103107400000000","0453503000000000"); // INT7 alpha 0-0.85
  } else if(trainConfig == 6305)  { //PHOS INT7, N.Pi cut var. Selection Window, Std 1 -> 2 sigma
    //                                                                                                      0103103400000000
    //                                                                                                             |
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa01cc00000","32c51070a","0103103k00000000","0453503000000000"); // INT7, 1.5 sigma
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa01cc00000","32c51070a","0103103l00000000","0453503000000000"); // INT7, 2.5 sigma
  } else if(trainConfig == 6306)  { //PHOS INT7, N.Pi cut var. open. angle, Std 0 -> off
    //                                                                                                      0103103400000000
    //                                                                                                                    |
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa01cc00000","32c51070a","0103103400000010","0453503000000000"); // INT7 Op. Ang. var 1: min opening angle 0.005
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa01cc00000","32c51070a","0103103400000020","0453503000000000"); // INT7 Op. Ang. var 3: min opening angle 0.01 -> 2 cell


    //-----
    //INT7: Omega Meson Cut Variations
    //-----
    //Std: 0453503000000000
  } else if(trainConfig == 6401)  { //PHOS INT7, Omega cut var. Background Scheme, Std 4 -> off
    //                                                                                                                         0453503000000000
    //                                                                                                                          |
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa01cc00000","32c51070a","0103103400000000","0153503000000000"); // INT7, Om Event Mixing
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa01cc00000","32c51070a","0103103400000000","0a53503000000000"); // INT7, Om LikeSignMixing
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa01cc00000","32c51070a","0103103400000000","0d53503000000000"); // INT7, Om SideBandMixing
  } else if(trainConfig == 6402)  { //PHOS INT7, Omega cut var. rapidity, Std 5 -> -0.85, 0.85
    //                                                                                                                         0453503000000000
    //                                                                                                                             |
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa01cc00000","32c51070a","0103103400000000","0453103000000000"); // INT7, Om rap. -0.8, 0.8
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa01cc00000","32c51070a","0103103400000000","0453603000000000"); // INT7, Om rap. -0.75, 0.75
  } else if(trainConfig == 6404)  { //PHOS INT7, Omega cut var. alpha, Std 3 -> 0.0-1.0
    //                                                                                                                         0453503000000000
    //                                                                                                                               |
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa01cc00000","32c51070a","0103103400000000","0453505000000000"); // INT7 alpha 0-0.75
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa01cc00000","32c51070a","0103103400000000","0453508000000000"); // INT7 alpha 0-0.6
  } else if(trainConfig == 6410)  { //PHOS INT7, Omega cut var. Background Scheme single cfg, Std 4 -> off, EventMixing
    //                                                                                                                         0453503000000000
    //                                                                                                                          |
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa01cc00000","32c51070a","0103103400000000","0153503000000000"); // INT7, Om Event Mixing
  } else if(trainConfig == 6411)  { //PHOS INT7, Omega cut var. Background Scheme single cfg, Std 4 -> off, LikeSignMixing
    //                                                                                                                         0453503000000000
    //                                                                                                                          |
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa01cc00000","32c51070a","0103103400000000","0a53503000000000"); // INT7, Om LikeSignMixing
  } else if(trainConfig == 6412)  { //PHOS INT7, Omega cut var. Background Scheme single cfg, Std 4 -> off, SideBandMixing
    //                                                                                                                         0453503000000000
    //                                                                                                                          |
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa01cc00000","32c51070a","0103103400000000","0d53503000000000"); // INT7, Om SideBandMixing

    //-----
    //INT7: PCM Conversion Cut
    //-----
    //Std: 0dm00009f9730000dge0404000
  } else if (trainConfig == 6501) {   // min pT variations
    //                                        0dm00009f9730000dge0404000
    //                                              |
    cuts.AddCutHeavyMesonPCMCalo("00010113", "0dm00069f9730000dge0404000", "24466190sa01cc00000", "32c51070a", "0103103400000000", "0453503000000000"); // eta < 0.8  // remove  55-72 bin, min pT 40 MeV
    cuts.AddCutHeavyMesonPCMCalo("00010113", "0dm00049f9730000dge0404000", "24466190sa01cc00000", "32c51070a", "0103103400000000", "0453503000000000"); // eta < 0.8  // remove  55-72 bin, min pT 75 MeV
    cuts.AddCutHeavyMesonPCMCalo("00010113", "0dm00019f9730000dge0404000", "24466190sa01cc00000", "32c51070a", "0103103400000000", "0453503000000000"); // eta < 0.8  // remove  55-72 bin, min pT 100MeV

  } else if (trainConfig == 6502) {   // TPC clusters, cosPA
    //                                        0dm00009f9730000dge0404000
    //                                               |
    cuts.AddCutHeavyMesonPCMCalo("00010113", "0dm00008f9730000dge0404000", "24466190sa01cc00000", "32c51070a", "0103103400000000", "0453503000000000"); // TPC cluster 35%
    cuts.AddCutHeavyMesonPCMCalo("00010113", "0dm00006f9730000dge0404000", "24466190sa01cc00000", "32c51070a", "0103103400000000", "0453503000000000"); // TPC cluster 70%
    cuts.AddCutHeavyMesonPCMCalo("00010113", "0dm00009f9730000dge0604000", "24466190sa01cc00000", "32c51070a", "0103103400000000", "0453503000000000"); // cosPA 0.9
    cuts.AddCutHeavyMesonPCMCalo("00010113", "0dm00009f9730000dge0304000", "24466190sa01cc00000", "32c51070a", "0103103400000000", "0453503000000000"); // cosPA 0.75

  } else if (trainConfig == 6503) {   // TPC clusters, cosPA
    //                                        0dm00009f9730000dge0404000
    //                                                |
    cuts.AddCutHeavyMesonPCMCalo("00010113", "0dm0000939730000dge0404000", "24466190sa01cc00000", "32c51070a", "0103103400000000", "0453503000000000"); // nsig electron   -4,5
    cuts.AddCutHeavyMesonPCMCalo("00010113", "0dm0000969730000dge0404000", "24466190sa01cc00000", "32c51070a", "0103103400000000", "0453503000000000"); // nsig electron -2.5,4
    cuts.AddCutHeavyMesonPCMCalo("00010113", "0dm00009f5730000dge0404000", "24466190sa01cc00000", "32c51070a", "0103103400000000", "0453503000000000"); // nsig pion 2,-10
    cuts.AddCutHeavyMesonPCMCalo("00010113", "0dm00009f1730000dge0404000", "24466190sa01cc00000", "32c51070a", "0103103400000000", "0453503000000000"); // nsig pion 0,-10

  } else if (trainConfig == 6504) {
    //                                        0dm00009f9730000dge0404000
    //                                                  |
    cuts.AddCutHeavyMesonPCMCalo("00010113", "0dm00009f9030000dge0404000", "24466190sa01cc00000", "32c51070a", "0103103400000000", "0453503000000000"); // pion nsig min mom 0.50 GeV/c
    cuts.AddCutHeavyMesonPCMCalo("00010113", "0dm00009f9630000dge0404000", "24466190sa01cc00000", "32c51070a", "0103103400000000", "0453503000000000"); // pion nsig min mom 0.25 GeV/c
    cuts.AddCutHeavyMesonPCMCalo("00010113", "0dm00009f9760000dge0404000", "24466190sa01cc00000", "32c51070a", "0103103400000000", "0453503000000000"); // pion nsig max mom 2.00 GeV/c
    cuts.AddCutHeavyMesonPCMCalo("00010113", "0dm00009f9710000dge0404000", "24466190sa01cc00000", "32c51070a", "0103103400000000", "0453503000000000"); // pion nsig max mom 5.00 GeV/c


  } else if (trainConfig == 6505) {   // chi2 variations
    //                                        0dm00009f9730000dge0404000
    //                                                         |
    cuts.AddCutHeavyMesonPCMCalo("00010113", "0dm00009f9730000d1e0404000", "24466190sa01cc00000", "32c51070a", "0103103400000000", "0453503000000000"); // chi2 50
    cuts.AddCutHeavyMesonPCMCalo("00010113", "0dm00009f9730000dfe0404000", "24466190sa01cc00000", "32c51070a", "0103103400000000", "0453503000000000"); // chi2 50 chi2 dep -0.065
    cuts.AddCutHeavyMesonPCMCalo("00010113", "0dm00009f9730000dhe0404000", "24466190sa01cc00000", "32c51070a", "0103103400000000", "0453503000000000"); // chi2 50 chi2 dep -0.050
    cuts.AddCutHeavyMesonPCMCalo("00010113", "0dm00009f9730000dge0400000", "24466190sa01cc00000", "32c51070a", "0103103400000000", "0453503000000000"); // remove reject close v0

  } else if (trainConfig == 6506) {   // Psi pair variations
    //                                        0dm00009f9730000dge0404000
    //                                                         ||
    cuts.AddCutHeavyMesonPCMCalo("00010113", "0dm00009f9730000dgd0404000", "24466190sa01cc00000", "32c51070a", "0103103400000000", "0453503000000000"); // Psi pair 0.15 dep
    cuts.AddCutHeavyMesonPCMCalo("00010113", "0dm00009f9730000dgf0404000", "24466190sa01cc00000", "32c51070a", "0103103400000000", "0453503000000000"); // Psi pair 0.20 dep
    cuts.AddCutHeavyMesonPCMCalo("00010113", "0dm00009f9730000dgg0404000", "24466190sa01cc00000", "32c51070a", "0103103400000000", "0453503000000000"); // Psi pair 0.30 dep
    cuts.AddCutHeavyMesonPCMCalo("00010113", "0dm00009227300008250404000", "24466190sa01cc00000", "32c51070a", "0103103400000000", "0453503000000000"); // old cuts (run1)

  } else if(trainConfig == 6511)  { //PHOS INT7 Standard 13TeV, Material Budget Studies
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa01cc00000","32c51070a","0103103400000000","0453503000000000"); // INT7 Standard
  } else if(trainConfig == 6512)  { //PHOS INT7 Standard 13TeV, Material Budget Studies
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa01cc00000","32c51070a","0103103400000000","0453503000000000"); // INT7 Standard
  } else if(trainConfig == 6513)  { //PHOS INT7 Standard 13TeV, Material Budget Studies
    cuts.AddCutHeavyMesonPCMCalo("00010113","0dm00009f9730000dge0404000","24466190sa01cc00000","32c51070a","0103103400000000","0453503000000000"); // INT7 Standard

    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

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

  Bool_t initializedMatBudWeigths_existing    = kFALSE;

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

    if (enableMatBudWeightsPi0 > 0){
      if (isMC > 0){
        Int_t FlagMatBudWeightsPi0=enableMatBudWeightsPi0;
        if (enableMatBudWeightsPi0>=10){
          FlagMatBudWeightsPi0-=10;
        }
        if (analysisCuts[i]->InitializeMaterialBudgetWeights(FlagMatBudWeightsPi0,fileNameMatBudWeights)){
          initializedMatBudWeigths_existing = kTRUE;
        }
        else {cout << "ERROR The initialization of the materialBudgetWeights did not work out." << endl;}
      }
      else {cout << "ERROR 'enableMatBudWeightsPi0'-flag was set > 0 even though this is not a MC task. It was automatically reset to 0." << endl;}
    }

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

  task->SetEnableSortingOfMCClusLabels(enableSortingMCLabels);

  if (initializedMatBudWeigths_existing) {
      task->SetDoMaterialBudgetWeightingOfGammasForTrueMesons(kTRUE);
      if (enableMatBudWeightsPi0>=10){
          task->SetDoMaterialBudgetWeightingOfGammasForInvMassHistogram(kTRUE);
      }
  }

  //connect containers
  AliAnalysisDataContainer *coutput =
  mgr->CreateContainer(Form("GammaConvNeutralMesonPiPlPiMiNeutralMeson_%i_%i_%i.root",selectHeavyNeutralMeson,neutralPionMode, trainConfig), TList::Class(),
              AliAnalysisManager::kOutputContainer,Form("GammaConvNeutralMesonPiPlPiMiNeutralMeson_%i_%i_%i.root",selectHeavyNeutralMeson,neutralPionMode, trainConfig));

  mgr->AddTask(task);
  mgr->ConnectInput(task,0,cinput);
  mgr->ConnectOutput(task,1,coutput);

  return;

}
