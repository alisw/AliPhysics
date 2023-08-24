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
    Int_t     enableQAMesonTask           = 1,                        // enable QA in AliAnalysisTaskNeutralMesonToPiPlPiMiNeutralMeson; 0: no QA, 1: general meson QA, 2: background QA, 3: 3D histogram, 4: Dalitz plots, 5: trees, 23: enable background calculations; 
                                                                      //    combinations: 6: 1+2, 7: 1+2+3, 8: 1+2+3+5, 9: 2+3, 10: 2+3+5, 11: 1+2+3+4+5
                                                                      //    QA can't be run with light output! 
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
    Bool_t    enableSortingMCLabels       = kTRUE,                    // enable sorting for MC cluster labels
    Bool_t    enableMLBckRedStudy         = kFALSE,                   // enable saving the output as tree for ML reduction study
    Int_t     MLBckRedStudyCutOff         = 10,                       // every which case that is not true meson should be saved
    TString   additionalTrainConfig       = "0"                       // additional counter for trainconfig, this has to be always the last parameter
  ) {

  AliCutHandlerPCM cuts(13);
  Bool_t usePionPreselection = kTRUE;
  //parse additionalTrainConfig flag
  Int_t trackMatcherRunningMode = 0; // CaloTrackMatcher running mode
  TString unsmearingoutputs = "0123"; // 0: No correction, 1: One pi0 mass errer subtracted, 2: pz of pi0 corrected to fix its mass, 3: Lambda(alpha)*DeltaPi0 subtracted

  TString addTaskName                       = "AddTask_GammaConvNeutralMesonPiPlPiMiNeutralMeson_CaloMode_pp";

  TString corrTaskSetting             = cuts.GetSpecialSettingFromAddConfig(additionalTrainConfig, "CF", "", addTaskName);
  if(corrTaskSetting.CompareTo(""))
    cout << "corrTaskSetting: " << corrTaskSetting.Data() << endl;

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
      } else if(tempStr.Contains("NoPionPreselection")){
        usePionPreselection = kFALSE;
      }
      if(tempStr.BeginsWith("UNSMEARING")){ // 0: No correction, 1: One pi0 mass errer subtracted, 2: pz of pi0 corrected to fix its mass, 3: Lambda(alpha)*DeltaPi0 subtracted
        TString tempType = tempStr;
        tempType.Replace(0,9,"");
        unsmearingoutputs = tempType;
        cout << "INFO: AddTask_GammaConvNeutralMesonPiPlPiMiPiZero_CaloMode_pp will output the following minv_pT histograms:" << endl;
        if(unsmearingoutputs.Contains("0")) cout << "- Uncorrected" << endl;
        if(unsmearingoutputs.Contains("1")) cout << "- SubNDM" << endl;
        if(unsmearingoutputs.Contains("2")) cout << "- Fixpz" << endl;
        if(unsmearingoutputs.Contains("3")) cout << "- SubLambda" << endl;
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

  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();

  //========= Add V0 Reader to  ANALYSIS manager if not yet existent =====
  TString cutnumberEvent = "00000003";
  TString V0ReaderName        = Form("V0ReaderV1_%s_%s",cutnumberEvent.Data(),photonCutNumberV0Reader.Data());
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
     periodNameV0Reader.Contains("LHC16c2") ||
     !usePionPreselection){
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
      if(runLightOutput>=1){
          fPionCuts->SetLightOutput(kTRUE);
      }
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
  if(runLightOutput>=2) {
      task->SetLightOutput(2);
  } else if(runLightOutput>=1) {
      task->SetLightOutput(1);
  }
  task->SetTolerance(tolerance);
  task->SetTrackMatcherRunningMode(trackMatcherRunningMode);
  task->SetCorrectionTaskSetting(corrTaskSetting);

  if(enableMLBckRedStudy && !isMC){
    cout << "Error: Trees for ML studies implemented only for MC. Returning..." << endl;
    return;
  }
  if(enableMLBckRedStudy){
    task->SetBckgReductionTree(MLBckRedStudyCutOff);
  }



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


  } else if(trainConfig == 210)  { //EDC 13TeV MB, NCell: v (NCell Cut 2, but with probability in MC to let clusters through )
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe3n230000","32c51070m","0103603l00000000","0453503000000000"); // INT7
  } else if(trainConfig == 211)  { //EDC 13TeV MB, NCell: v (NCell Cut 2, but with probability in MC to let clusters through )
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe3n230000","32c51070m","0103603l00000000","0453503000000000"); // EG2
  } else if(trainConfig == 212)  { //EDC 13TeV MB, NCell: v (NCell Cut 2, but with probability in MC to let clusters through )
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe3n230000","32c51070m","0103603l00000000","0453503000000000"); // EG1
  } else if(trainConfig == 213)  { //EDC 13TeV MB, NCell: 0
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe30220000","32c51070m","0103603l00000000","0453503000000000"); // INT7
  } else if(trainConfig == 214)  { //EDC 13TeV MB, NCell: 0
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe30220000","32c51070m","0103603l00000000","0453503000000000"); // EG2
  } else if(trainConfig == 215)  { //EDC 13TeV MB, NCell: 0
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe30220000","32c51070m","0103603l00000000","0453503000000000"); // EG1

    // PCM-PHOS
  } else if ( trainConfig == 250 ) { // INT7 + PHI7
    cuts.AddCutHeavyMesonCalo("00010113","24466190wa01cc00000","32c510700","0103603l00000010","0453503000000010"); // INT7
    cuts.AddCutHeavyMesonCalo("00062113","24466190wa01cc00000","32c510700","0103603l00000010","0453503000000010"); // PHI7

  // pp 5 TeV
  } else if( trainConfig == 260 ) {
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe30220000","32c51070a","0103603l00000000","0453503000000000");
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe30220000","32c51070a","0103603l00000000","0453503000000000"); // EG2
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe30220000","32c51070a","0103603l00000000","0453503000000000"); // EG1

  } else if( trainConfig == 261 ) { // Pion Mass cut as for Omega
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe30220000","32l51070a","0103603l00000000","0453503000000000"); // trainConfig == 2000
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe30220000","32l51070a","0103603l00000000","0453503000000000"); // trainConfig == 3000
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe30220000","32l51070a","0103603l00000000","0453503000000000"); // trainConfig == 4000

  } else if( trainConfig == 262 ) {
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe30220000","32c51070a","0103603l00000000","0453503000000000");
  } else if( trainConfig == 263 ) {
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe30220000","32c51070a","0103603l00000000","0453503000000000"); // EG2
  } else if( trainConfig == 264 ) {
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe30220000","32c51070a","0103603l00000000","0453503000000000"); // EG1


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

  } else if( trainConfig == 310) { // No Pion Mass cut
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe30220000","32l510700","0103103x00000000","0453503000000000"); // trainConfig == 2000
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe30220000","32l510700","01031v3x00000000","0453503000000000"); // trainConfig == 3000
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe30220000","32l510700","01031v3x00000000","0453503000000000"); // trainConfig == 4000

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
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe32220000","32c51070a","0103103100000000","0153503000000000"); // INT7
  } else if(trainConfig == 631)  { // EDC 13 TeV EG2
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe32220000","32c51070a","0103103100000000","0153503000000000"); // EG2
  } else if(trainConfig == 632)  { // EDC 13 TeV EG1
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe32220000","32c51070a","0103103100000000","0153503000000000"); // EG1

    //Gamma Energy Cuts EG2
  } else if(trainConfig == 633)  { //Standard EMCal 13TeV Trigger EG2, GammaCut > 4GeV
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe32220000","32c51070a","01031k3100000000","0153503000000000"); // EG2, new Gamma Energy cut 4 GeV
  } else if(trainConfig == 634)  { //Standard EMCal 13TeV Trigger EG2, GammaCut > 5GeV
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe32220000","32c51070a","01031e3100000000","0153503000000000"); // EG2, new Gamma Energy cut 5 GeV
  } else if(trainConfig == 635)  { //Standard EMCal 13TeV Trigger EG2, GammaCut > 6GeV
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe32220000","32c51070a","01031g3100000000","0153503000000000"); // EG2, new Gamma Energy cut 6 GeV
  } else if(trainConfig == 636)  { //Standard EMCal 13TeV Trigger EG2, GammaCut > 7.5GeV
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe32220000","32c51070a","01031f3100000000","0153503000000000"); // EG2, new Gamma Energy cut 7.5 GeV

    //Gamma Energy Cuts EG1
  } else if(trainConfig == 637)  { //EMCal Trig Pt Cut Variations EG1, GammaCut > 8GeV
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe32220000","32c51070a","01031l3100000000","0153503000000000"); // EG1, new Gamma Energy cut 8 GeV
  } else if(trainConfig == 638)  { //EMCal Trig Pt Cut Variations EG1, GammaCut > 9GeV
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe32220000","32c51070a","01031m3100000000","0153503000000000"); // EG1, new Gamma Energy cut 9 GeV
  } else if(trainConfig == 639)  { //EMCal Trig Pt Cut Variations EG1, GammaCut > 10GeV
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe32220000","32c51070a","01031h3100000000","0153503000000000"); // EG1, new Gamma Energy cut 10 GeV

    //no shared TPC clusters, Shared cluster Fraction =0
  } else if(trainConfig == 640)  { // EDC 13 TeV, Shared cluster Fraction =0
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe32220000","32e51070a","0103103100000000","0153503000000000"); // INT7
  } else if(trainConfig == 641)  { //Standard EMCal 13TeV Trigger EG2, GammaCut > 6GeV, Shared cluster Fraction =0
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe32220000","32e51070a","01031g3100000000","0153503000000000"); // EG2, new Gamma Energy cut 6 GeV
  } else if(trainConfig == 642)  { //Standard EMCal 13TeV Trigger EG1, GammaCut > 10GeV, Shared cluster Fraction =0
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe32220000","32e51070a","01031h3100000000","0153503000000000"); // EG1, new Gamma Energy cut 10 GeV

    //no shared TPC clusters, Shared cluster Fraction <=0.4
  } else if(trainConfig == 643)  { // EDC 13 TeV, Shared cluster Fraction <=0.4
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe32220000","32f51070a","0103103100000000","0153503000000000"); // INT7
  } else if(trainConfig == 644)  { //Standard EMCal 13TeV Trigger EG2, GammaCut > 6GeV, Shared cluster Fraction <=0.4
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe32220000","32f51070a","01031g3100000000","0153503000000000"); // EG2, new Gamma Energy cut 6 GeV
  } else if(trainConfig == 645)  { //Standard EMCal 13TeV Trigger EG1, GammaCut > 10GeV, Shared cluster Fraction <=0.4
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe32220000","32f51070a","01031h3100000000","0153503000000000"); // EG1, new Gamma Energy cut 10 GeV

    //Charged Pion Mass Cut <600 MeV
  } else if(trainConfig == 646)  { // EDC 13 TeV
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe32220000","32c51070b","0103103100000000","0153503000000000"); // INT7
  } else if(trainConfig == 647)  { // EDC 13 TeV EG2
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe32220000","32c51070b","0103103100000000","0153503000000000"); // EG2
  } else if(trainConfig == 648)  { // EDC 13 TeV EG1
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe32220000","32c51070b","0103103100000000","0153503000000000"); // EG1

    //Standard Cuts of Pi0 Analysis: ("00010113","411790106fe32220000","0r631031000000d0")
    //MesonCut r63==Background->ignored, d==OpeningAngle for Background->ignored
  } else if(trainConfig == 654)  { // EDC 13 TeV
    cuts.AddCutHeavyMesonCalo("00010113","411790106fe32220000","32c51070a","0103103100000000","0153503000000000"); // INT7
  } else if(trainConfig == 655)  { // EDC 13 TeV EG2
    cuts.AddCutHeavyMesonCalo("0008e113","411790106fe32220000","32c51070a","0103103100000000","0153503000000000"); // EG2
  } else if(trainConfig == 656)  { // EDC 13 TeV EG1
    cuts.AddCutHeavyMesonCalo("0008d113","411790106fe32220000","32c51070a","0103103100000000","0153503000000000"); // EG1

    //Charged Pion Mass Cut <600 MeV
  } else if(trainConfig == 660)  { // EDC 13 TeV
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe32220000","32c51070b","0103103100000000","0153503000000000"); // INT7
  } else if(trainConfig == 661)  { // EDC 13 TeV EG2
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe32220000","32c51070b","0103103100000000","0153503000000000"); // EG2
  } else if(trainConfig == 662)  { // EDC 13 TeV EG1
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe32220000","32c51070b","0103103100000000","0153503000000000"); // EG1
    //Charged Pion Mass Cut <850 MeV, Neutral Charged Pion Cut 1000< MeV
  } else if(trainConfig == 663)  { // EDC 13 TeV
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe32220000","32c51070c","0103103100000000","0153503000000000"); // INT7
  } else if(trainConfig == 664)  { // EDC 13 TeV EG2
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe32220000","32c51070c","0103103100000000","0153503000000000"); // EG2
  } else if(trainConfig == 665)  { // EDC 13 TeV EG1
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe32220000","32c51070c","0103103100000000","0153503000000000"); // EG1
    //Charged Pion Mass Cut <850 MeV, Neutral Charged Pion Cut <850 MeV
  } else if(trainConfig == 666)  { // EDC 13 TeV
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe32220000","32c51070d","0103103100000000","0153503000000000"); // INT7
  } else if(trainConfig == 667)  { // EDC 13 TeV EG2
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe32220000","32c51070d","0103103100000000","0153503000000000"); // EG2
  } else if(trainConfig == 668)  { // EDC 13 TeV EG1
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe32220000","32c51070d","0103103100000000","0153503000000000"); // EG1
    //Charged Pion Mass Cut <600 MeV, Neutral Charged Pion Cut <600 MeV
  } else if(trainConfig == 669)  { // EDC 13 TeV
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe32220000","32c51070e","0103103100000000","0153503000000000"); // INT7
  } else if(trainConfig == 670)  { // EDC 13 TeV EG2
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe32220000","32c51070e","0103103100000000","0153503000000000"); // EG2
  } else if(trainConfig == 671)  { // EDC 13 TeV EG1
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe32220000","32c51070e","0103103100000000","0153503000000000"); // EG1

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
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe32220000","32c51070a","0103103100000000","0a53503000000000"); // INT7
  } else if(trainConfig == 731)  { // EDC 13 TeV EG2
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe32220000","32c51070a","0103103100000000","0a53503000000000"); // EG2
  } else if(trainConfig == 732)  { // EDC 13 TeV EG1
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe32220000","32c51070a","0103103100000000","0a53503000000000"); // EG1

    //Gamma Energy Cuts EG2
  } else if(trainConfig == 733)  { //Standard EMCal 13TeV Trigger EG2, GammaCut > 4GeV
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe32220000","32c51070a","01031k3100000000","0a53503000000000"); // EG2, new Gamma Energy cut 4 GeV
  } else if(trainConfig == 734)  { //Standard EMCal 13TeV Trigger EG2, GammaCut > 5GeV
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe32220000","32c51070a","01031e3100000000","0a53503000000000"); // EG2, new Gamma Energy cut 5 GeV
  } else if(trainConfig == 735)  { //Standard EMCal 13TeV Trigger EG2, GammaCut > 6GeV
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe32220000","32c51070a","01031g3100000000","0a53503000000000"); // EG2, new Gamma Energy cut 6 GeV
  } else if(trainConfig == 736)  { //Standard EMCal 13TeV Trigger EG2, GammaCut > 7.5GeV
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe32220000","32c51070a","01031f3100000000","0a53503000000000"); // EG2, new Gamma Energy cut 7.5 GeV

    //Gamma Energy Cuts EG1
  } else if(trainConfig == 737)  { //EMCal Trig Pt Cut Variations EG1, GammaCut > 8GeV
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe32220000","32c51070a","01031l3100000000","0a53503000000000"); // EG1, new Gamma Energy cut 8 GeV
  } else if(trainConfig == 738)  { //EMCal Trig Pt Cut Variations EG1, GammaCut > 9GeV
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe32220000","32c51070a","01031m3100000000","0a53503000000000"); // EG1, new Gamma Energy cut 9 GeV
  } else if(trainConfig == 739)  { //EMCal Trig Pt Cut Variations EG1, GammaCut > 10GeV
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe32220000","32c51070a","01031h3100000000","0a53503000000000"); // EG1, new Gamma Energy cut 10 GeV

    //no shared TPC clusters, Shared cluster Fraction =0
  } else if(trainConfig == 740)  { // EDC 13 TeV, Shared cluster Fraction =0
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe32220000","32e51070a","0103103100000000","0a53503000000000"); // INT7
  } else if(trainConfig == 741)  { //Standard EMCal 13TeV Trigger EG2, GammaCut > 6GeV, Shared cluster Fraction =0
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe32220000","32e51070a","01031g3100000000","0a53503000000000"); // EG2, new Gamma Energy cut 6 GeV
  } else if(trainConfig == 742)  { //Standard EMCal 13TeV Trigger EG1, GammaCut > 10GeV, Shared cluster Fraction =0
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe32220000","32e51070a","01031h3100000000","0a53503000000000"); // EG1, new Gamma Energy cut 10 GeV

    //no shared TPC clusters, Shared cluster Fraction <=0.4
  } else if(trainConfig == 743)  { // EDC 13 TeV, Shared cluster Fraction <=0.4
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe32220000","32f51070a","0103103100000000","0a53503000000000"); // INT7
  } else if(trainConfig == 744)  { //Standard EMCal 13TeV Trigger EG2, GammaCut > 6GeV, Shared cluster Fraction <=0.4
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe32220000","32f51070a","01031g3100000000","0a53503000000000"); // EG2, new Gamma Energy cut 6 GeV
  } else if(trainConfig == 745)  { //Standard EMCal 13TeV Trigger EG1, GammaCut > 10GeV, Shared cluster Fraction <=0.4
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe32220000","32f51070a","01031h3100000000","0a53503000000000"); // EG1, new Gamma Energy cut 10 GeV

    //Charged Pion Mass Cut <600 MeV
  } else if(trainConfig == 746)  { // EDC 13 TeV
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe32220000","32c51070b","0103103100000000","0a53503000000000"); // INT7
  } else if(trainConfig == 747)  { // EDC 13 TeV EG2
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe32220000","32c51070b","0103103100000000","0a53503000000000"); // EG2
  } else if(trainConfig == 748)  { // EDC 13 TeV EG1
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe32220000","32c51070b","0103103100000000","0a53503000000000"); // EG1

    //Standard Cuts of Pi0 Analysis: ("00010113","411790106fe32220000","0r631031000000d0")
    //MesonCut r63==Background->ignored, d==OpeningAngle for Background->ignored
  } else if(trainConfig == 754)  { // EDC 13 TeV
    cuts.AddCutHeavyMesonCalo("00010113","411790106fe32220000","32c51070a","0103103100000000","0a53503000000000"); // INT7
  } else if(trainConfig == 755)  { // EDC 13 TeV EG2
    cuts.AddCutHeavyMesonCalo("0008e113","411790106fe32220000","32c51070a","0103103100000000","0a53503000000000"); // EG2
  } else if(trainConfig == 756)  { // EDC 13 TeV EG1
    cuts.AddCutHeavyMesonCalo("0008d113","411790106fe32220000","32c51070a","0103103100000000","0a53503000000000"); // EG1



    //Charged Pion Mass Cut <600 MeV
  } else if(trainConfig == 760)  { // EDC 13 TeV
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe32220000","32c51070b","0103103100000000","0a53503000000000"); // INT7
  } else if(trainConfig == 761)  { // EDC 13 TeV EG2
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe32220000","32c51070b","0103103100000000","0a53503000000000"); // EG2
  } else if(trainConfig == 762)  { // EDC 13 TeV EG1
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe32220000","32c51070b","0103103100000000","0a53503000000000"); // EG1
    //Charged Pion Mass Cut <850 MeV, Neutral Charged Pion Cut 1000< MeV
  } else if(trainConfig == 763)  { // EDC 13 TeV
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe32220000","32c51070c","0103103100000000","0a53503000000000"); // INT7
  } else if(trainConfig == 764)  { // EDC 13 TeV EG2
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe32220000","32c51070c","0103103100000000","0a53503000000000"); // EG2
  } else if(trainConfig == 765)  { // EDC 13 TeV EG1
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe32220000","32c51070c","0103103100000000","0a53503000000000"); // EG1
    //Charged Pion Mass Cut <850 MeV, Neutral Charged Pion Cut <850 MeV
  } else if(trainConfig == 766)  { // EDC 13 TeV
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe32220000","32c51070d","0103103100000000","0a53503000000000"); // INT7
  } else if(trainConfig == 767)  { // EDC 13 TeV EG2
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe32220000","32c51070d","0103103100000000","0a53503000000000"); // EG2
  } else if(trainConfig == 768)  { // EDC 13 TeV EG1
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe32220000","32c51070d","0103103100000000","0a53503000000000"); // EG1
    //Charged Pion Mass Cut <600 MeV, Neutral Charged Pion Cut <600 MeV
  } else if(trainConfig == 769)  { // EDC 13 TeV
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe32220000","32c51070e","0103103100000000","0a53503000000000"); // INT7
  } else if(trainConfig == 770)  { // EDC 13 TeV EG2
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe32220000","32c51070e","0103103100000000","0a53503000000000"); // EG2
  } else if(trainConfig == 771)  { // EDC 13 TeV EG1
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe32220000","32c51070e","0103103100000000","0a53503000000000"); // EG1

    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // EMC pp 13 TeV Sideband Mixing
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //Standard Cuts
  } else if(trainConfig == 830)  { // EDC 13 TeV
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe32220000","32c51070a","0103103100000000","0d53503000000000"); // INT7, Sideband, both sides
  } else if(trainConfig == 831)  { // EDC 13 TeV EG2
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe32220000","32c51070a","0103103100000000","0d53503000000000"); // EG2, Sideband, both sides
  } else if(trainConfig == 832)  { // EDC 13 TeV EG1
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe32220000","32c51070a","0103103100000000","0d53503000000000"); // EG1, Sideband, both sides

  } else if(trainConfig == 833)  { // EDC 13 TeV
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe32220000","32c51070a","0103103100000000","0b53503000000000"); // INT7, Sideband, right side
  } else if(trainConfig == 834)  { // EDC 13 TeV
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe32220000","32c51070a","0103103100000000","0c53503000000000"); // INT7, Sideband, left side



    //Charged Pion Mass Cut <600 MeV
  } else if(trainConfig == 860)  { // EDC 13 TeV
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe32220000","32c51070b","0103103100000000","0b53503000000000"); // INT7
  } else if(trainConfig == 861)  { // EDC 13 TeV EG2
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe32220000","32c51070b","0103103100000000","0b53503000000000"); // EG2
  } else if(trainConfig == 862)  { // EDC 13 TeV EG1
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe32220000","32c51070b","0103103100000000","0b53503000000000"); // EG1
    //Charged Pion Mass Cut <850 MeV, Neutral Charged Pion Cut 1000< MeV
  } else if(trainConfig == 863)  { // EDC 13 TeV
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe32220000","32c51070c","0103103100000000","0b53503000000000"); // INT7
  } else if(trainConfig == 864)  { // EDC 13 TeV EG2
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe32220000","32c51070c","0103103100000000","0b53503000000000"); // EG2
  } else if(trainConfig == 865)  { // EDC 13 TeV EG1
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe32220000","32c51070c","0103103100000000","0b53503000000000"); // EG1
    //Charged Pion Mass Cut <850 MeV, Neutral Charged Pion Cut <850 MeV
  } else if(trainConfig == 866)  { // EDC 13 TeV
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe32220000","32c51070d","0103103100000000","0b53503000000000"); // INT7
  } else if(trainConfig == 867)  { // EDC 13 TeV EG2
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe32220000","32c51070d","0103103100000000","0b53503000000000"); // EG2
  } else if(trainConfig == 868)  { // EDC 13 TeV EG1
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe32220000","32c51070d","0103103100000000","0b53503000000000"); // EG1
    //Charged Pion Mass Cut <600 MeV, Neutral Charged Pion Cut <600 MeV
  } else if(trainConfig == 869)  { // EDC 13 TeV
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe32220000","32c51070e","0103103100000000","0b53503000000000"); // INT7
  } else if(trainConfig == 870)  { // EDC 13 TeV EG2
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe32220000","32c51070e","0103103100000000","0b53503000000000"); // EG2
  } else if(trainConfig == 871)  { // EDC 13 TeV EG1
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe32220000","32c51070e","0103103100000000","0b53503000000000"); // EG1


    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // EMC pp 13 TeV Fitting, Cut Variations very low EMCAL Pt
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  } else if(trainConfig == 900)  { //Standard EDC 13TeV MB, Std Cut
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe32220000","32c51070a","0103103100000000","0453503000000000"); // INT7
    //-----
    //Calo Variations
    //-----
  } else if(trainConfig == 901)  { //EDC 13TeV MB, timing: 6 (-30+35) -> 411790106fe32220000
    cuts.AddCutHeavyMesonCalo("00010113","411790106fe32220000","32c51070a","0103103100000000","0453503000000000"); // INT7
  } else if(trainConfig == 902)  { //EDC 13TeV MB, MinE: c (0.65 GeV) -> 411790109fec2220000
    cuts.AddCutHeavyMesonCalo("00010113","411790109fec2220000","32c51070a","0103103100000000","0453503000000000"); // INT7
  } else if(trainConfig == 903)  { //EDC 13TeV MB, MinE: c (0.8 GeV) -> 411790109fe42220000
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe42220000","32c51070a","0103103100000000","0453503000000000"); // INT7
  } else if(trainConfig == 904)  { //EDC 13TeV MB, NCell: 0 (no cut up to 4 GeV ) -> 411790109fe30220000
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe30220000","32c51070a","0103103100000000","0453503000000000"); // INT7
  } else if(trainConfig == 905)  { //EDC 13TeV MB, M02: 3 (0.5)  -> 411790109fe32230000
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe32230000","32c51070a","0103103100000000","0453503000000000"); // INT7
  } else if(trainConfig == 906)  { //EDC 13TeV MB, M02: 1 (1.0) -> 411790109fe32210000
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe32210000","32c51070a","0103103100000000","0453503000000000"); // INT7
  } else if(trainConfig == 907)  { //EDC 13TeV MB, NCell: v (NCell Cut 2, but with probability in MC to let clusters through ) -> 411790109fe3n220000
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe3n220000","32c51070a","0103103100000000","0453503000000000"); // INT7
    //-----
    //Meson Cut Variations
    //-----
  } else if(trainConfig == 910)  { //EDC 13TeV MB, rapidity 5 (0.85)
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe32220000","32c51070a","0103503100000000","0453503000000000"); // INT7
  } else if(trainConfig == 911)  { //EDC 13TeV MB, rapidity 6 (0.75)
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe32220000","32c51070a","0103603100000000","0453503000000000"); // INT7
  } else if(trainConfig == 912)  { //EDC 13TeV MB, non pt dependent pi0 mass 1 (0.10-0.145)
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe32220000","32c51070a","0103103100000000","0453503000000000"); // INT7
  } else if(trainConfig == 913)  { //EDC 13TeV MB, non pt dependent pi0 mass 2 (0.11-0.145)
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe32220000","32c51070a","0103103200000000","0453503000000000"); // INT7
  } else if(trainConfig == 914)  { //EDC 13TeV MB, non pt dependent pi0 mass 3 (0.12-0.145)
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe32220000","32c51070a","0103103300000000","0453503000000000"); // INT7
  } else if(trainConfig == 915)  { //EDC 13TeV MB, non pt dependent pi0 mass 4 (0.10-0.150)
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe32220000","32c51070a","0103103400000000","0453503000000000"); // INT7

    //-----
    //Primary Pion Variations
    //-----
  } else if(trainConfig == 920)  { //EDC 13TeV MB, Ch.Pi 650MeV
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe32220000","32c51070f","0103103100000000","0453503000000000"); // INT7
  } else if(trainConfig == 921)  { //EDC 13TeV MB, Ch.Pi 700MeV
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe32220000","32c51070g","0103103100000000","0453503000000000"); // INT7
  } else if(trainConfig == 922)  { //EDC 13TeV MB, Ch.Pi 750MeV
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe32220000","32c51070a","0103103100000000","0453503000000000"); // INT7
  } else if(trainConfig == 923)  { //EDC 13TeV MB, Ch.Pi 850MeV, Ch.Neu.Pi 1000MeV
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe32220000","32c51070c","0103103100000000","0453503000000000"); // INT7
  } else if(trainConfig == 924)  { //EDC 13TeV MB, Ch.Pi 650MeV, Ch.Neu.Pi 1000MeV
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe32220000","32c51070h","0103103100000000","0453503000000000"); // INT7
  } else if(trainConfig == 925)  { //EDC 13TeV MB, Ch.Pi 650MeV, Ch.Neu.Pi 850MeV
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe32220000","32c51070i","0103103100000000","0453503000000000"); // INT7
  } else if(trainConfig == 926)  { //EDC 13TeV MB, Ch.Pi 460MeV
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe32220000","32c51070j","0103103100000000","0453503000000000"); // INT7
  } else if(trainConfig == 927)  { //EDC 13TeV MB, Ch.Pi 480MeV
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe32220000","32c51070k","0103103100000000","0453503000000000"); // INT7
  } else if(trainConfig == 928)  { //EDC 13TeV MB, Ch.Pi 520MeV
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe32220000","32c51070l","0103103100000000","0453503000000000"); // INT7
  } else if(trainConfig == 930)  { //EDC 13TeV MB, Ch.Pi 650MeV, TOF Cut 3
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe32220000","32c51073f","0103103100000000","0453503000000000"); // INT7
  } else if(trainConfig == 931)  { //EDC 13TeV MB, Ch.Pi 460MeV, TOF Cut 3
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe32220000","32c51073j","0103103100000000","0453503000000000"); // INT7
  } else if(trainConfig == 932)  { //EDC 13TeV MB, Ch.Pi 480MeV, TOF Cut 3
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe32220000","32c51073k","0103103100000000","0453503000000000"); // INT7
  } else if(trainConfig == 933)  { //EDC 13TeV MB, Ch.Pi 520MeV, TOF Cut 3
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe32220000","32c51073l","0103103100000000","0453503000000000"); // INT7

    //-----
    //Omega Meson Cut Variations
    //-----
  } else if(trainConfig == 950)  { //EDC 13TeV MB, rapidity 5 (0.85)
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe32220000","32c51070a","0103103100000000","045350j000000000"); // INT7
  } else if(trainConfig == 951)  { //EDC 13TeV MB, rapidity 5 (0.85)
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe32220000","32c51070a","0103103100000000","045350k000000000"); // INT7

    //NL Changes
    //NCell 2, NL 21
  } else if(trainConfig == 960)  { //Standard EDC 13TeV MB, Std Cut
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe32220000","32c51070a","0103103100000000","0453503000000000"); // INT7
    //NCell 2, NL 01
  } else if(trainConfig == 963)  { //Standard EDC 13TeV MB, Std Cut
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe32220000","32c51070a","0103103100000000","0453503000000000"); // INT7
    //NCell 0, NL 21
  } else if(trainConfig == 970)  { //Standard EDC 13TeV MB, Std Cut
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe30220000","32c51070a","0103103100000000","0453503000000000"); // INT7
    //NCell 0, NL 01
  } else if(trainConfig == 973)  { //Standard EDC 13TeV MB, Std Cut
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe30220000","32c51070a","0103103100000000","0453503000000000"); // INT7
    //NCell 2 + NCell eficiency, NL 01
  } else if(trainConfig == 983)  { //Standard EDC 13TeV MB, Std Cut
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe3n220000","32c51070a","0103103100000000","0453503000000000"); // INT7
    //-----
    //Trigger EG2: Calo Variations
    //-----
  } else if(trainConfig == 1001)  { //EDC 13TeV MB, timing: 6 (-30+35) -> 411790106fe32220000
    cuts.AddCutHeavyMesonCalo("0008e113","411790106fe32220000","32c51070a","0103103100000000","0453503000000000"); // INT7
  } else if(trainConfig == 1002)  { //EDC 13TeV MB, MinE: c (0.65 GeV) -> 411790109fec2220000
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fec2220000","32c51070a","0103103100000000","0453503000000000"); // INT7
  } else if(trainConfig == 1003)  { //EDC 13TeV MB, MinE: c (0.8 GeV) -> 411790109fe42220000
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe42220000","32c51070a","0103103100000000","0453503000000000"); // INT7
  } else if(trainConfig == 1004)  { //EDC 13TeV MB, NCell: 0 (no cut up to 4 GeV ) -> 411790109fe30220000
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe30220000","32c51070a","0103103100000000","0453503000000000"); // INT7
  } else if(trainConfig == 1005)  { //EDC 13TeV MB, M02: 3 (0.5)  -> 411790109fe32230000
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe32230000","32c51070a","0103103100000000","0453503000000000"); // INT7
  } else if(trainConfig == 1006)  { //EDC 13TeV MB, M02: 1 (1.0) -> 411790109fe32210000
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe32210000","32c51070a","0103103100000000","0453503000000000"); // INT7
  } else if(trainConfig == 1007)  { //EDC 13TeV MB, NCell: v (NCell Cut 2, but with probability in MC to let clusters through ) -> 411790109fe3n220000
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe3n220000","32c51070a","0103103100000000","0453503000000000"); // INT7
    //-----
    //Trigger EG2: Meson Cut Variations
    //-----
  } else if(trainConfig == 1010)  { //EDC 13TeV MB, rapidity 5 (0.85)
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe32220000","32c51070a","0103503100000000","0453503000000000"); // INT7
  } else if(trainConfig == 1011)  { //EDC 13TeV MB, rapidity 6 (0.75)
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe32220000","32c51070a","0103603100000000","0453503000000000"); // INT7
  } else if(trainConfig == 1012)  { //EDC 13TeV MB, non pt dependent pi0 mass 1 (0.10-0.145)
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe32220000","32c51070a","0103103100000000","0453503000000000"); // INT7
  } else if(trainConfig == 1013)  { //EDC 13TeV MB, non pt dependent pi0 mass 2 (0.11-0.145)
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe32220000","32c51070a","0103103200000000","0453503000000000"); // INT7
  } else if(trainConfig == 1014)  { //EDC 13TeV MB, non pt dependent pi0 mass 3 (0.12-0.145)
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe32220000","32c51070a","0103103300000000","0453503000000000"); // INT7
  } else if(trainConfig == 1015)  { //EDC 13TeV MB, non pt dependent pi0 mass 4 (0.10-0.150)
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe32220000","32c51070a","0103103400000000","0453503000000000"); // INT7
  } else if(trainConfig == 1016)  { //EDC 13TeV MB, min/max pt cut n (min gamma 5gev, max pt 20gev)
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe32220000","32c51070a","01031n3100000000","0453503000000000"); // INT7
  } else if(trainConfig == 1017)  { //EDC 13TeV MB, min/max pt cut s (no min pt, max pt 20gev)
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe32220000","32c51070a","01031s3100000000","0453503000000000"); // INT7
  } else if(trainConfig == 1018)  { //EDC 13TeV MB, min/max pt cut v (no min pt, max pt 20gev)
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe32220000","32c51070a","01031v3100000000","0453503000000000"); // INT7

    //-----
    //Trigger EG2: Primary Pion Variations
    //-----
  } else if(trainConfig == 1020)  { //EDC 13TeV MB, Ch.Pi 650MeV
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe32220000","32c51070f","0103103100000000","0453503000000000"); // INT7
  } else if(trainConfig == 1021)  { //EDC 13TeV MB, Ch.Pi 700MeV
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe32220000","32c51070g","0103103100000000","0453503000000000"); // INT7
  } else if(trainConfig == 1022)  { //EDC 13TeV MB, Ch.Pi 750MeV
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe32220000","32c51070a","0103103100000000","0453503000000000"); // INT7
  } else if(trainConfig == 1023)  { //EDC 13TeV MB, Ch.Pi 850MeV, Ch.Neu.Pi 1000MeV
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe32220000","32c51070c","0103103100000000","0453503000000000"); // INT7
  } else if(trainConfig == 1024)  { //EDC 13TeV MB, Ch.Pi 650MeV, Ch.Neu.Pi 1000MeV
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe32220000","32c51070h","0103103100000000","0453503000000000"); // INT7
  } else if(trainConfig == 1025)  { //EDC 13TeV MB, Ch.Pi 650MeV, Ch.Neu.Pi 850MeV
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe32220000","32c51070i","0103103100000000","0453503000000000"); // INT7

    //NL Changes
    //NCell 2, NL 21
  } else if(trainConfig == 1060)  { //EDC 13TeV MB, NCell: 2
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe32220000","32c51070a","0103103100000000","0453503000000000"); // INT7
  } else if(trainConfig == 1061)  { //EDC 13TeV MB, min/max pt cut s (no min pt, max pt 20gev), NCell Cut 2
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe32220000","32c51070a","01031s3100000000","0453503000000000"); // INT7
  } else if(trainConfig == 1062)  { //EDC 13TeV MB, min/max pt cut v (no min pt, max pt  25gev), NCell Cut 2
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe32220000","32c51070a","01031v3100000000","0453503000000000"); // INT7
    //NCell 2, NL 01 (new std)
  } else if(trainConfig == 1063)  { //EDC 13TeV MB, NCell: 2
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe32220000","32c51070a","0103103100000000","0453503000000000"); // INT7
  } else if(trainConfig == 1064)  { //EDC 13TeV MB, min/max pt cut s (no min pt, max pt 20gev), NCell Cut 2
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe32220000","32c51070a","01031s3100000000","0453503000000000"); // INT7
  } else if(trainConfig == 1065)  { //EDC 13TeV MB, min/max pt cut v (no min pt, max pt  25gev), NCell Cut 2
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe32220000","32c51070a","01031v3100000000","0453503000000000"); // INT7

    //NCell 0, NL 21
  } else if(trainConfig == 1070)  { //EDC 13TeV MB, NCell: 0 (no cut up to 4 GeV ) -> 411790109fe30220000
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe30220000","32c51070a","0103103100000000","0453503000000000"); // INT7
  } else if(trainConfig == 1071)  { //EDC 13TeV MB, min/max pt cut s (no min pt, max pt 20gev), NCell Cut 0
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe30220000","32c51070a","01031s3100000000","0453503000000000"); // INT7
  } else if(trainConfig == 1072)  { //EDC 13TeV MB, min/max pt cut v (no min pt, max pt  25gev), NCell Cut 0
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe30220000","32c51070a","01031v3100000000","0453503000000000"); // INT7
    //NCell 0, NL 01 (new std)
  } else if(trainConfig == 1073)  { //EDC 13TeV MB, NCell: 0 (no cut up to 4 GeV ) -> 411790109fe30220000
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe30220000","32c51070a","0103103100000000","0453503000000000"); // INT7
  } else if(trainConfig == 1074)  { //EDC 13TeV MB, min/max pt cut s (no min pt, max pt 20gev), NCell Cut 0
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe30220000","32c51070a","01031s3100000000","0453503000000000"); // INT7
  } else if(trainConfig == 1075)  { //EDC 13TeV MB, min/max pt cut v (no min pt, max pt  25gev), NCell Cut 0
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe30220000","32c51070a","01031v3100000000","0453503000000000"); // INT7

    //NCell 2 + NCell effi, NL 01
  } else if(trainConfig == 1083)  { //EDC 13TeV EG2, NCell: n NCell2 + NCell effi in MC
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe3n220000","32c51070a","0103103100000000","0453503000000000"); // EG2
  } else if(trainConfig == 1084)  { //EDC 13TeV EG2, min/max pt cut s (no min pt, max pt 20gev), NCell: n NCell2 + NCell effi in MC
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe3n220000","32c51070a","01031s3100000000","0453503000000000"); // EG2
  } else if(trainConfig == 1085)  { //EDC 13TeV EG2, min/max pt cut v (no min pt, max pt  25gev), NCell: n NCell2 + NCell effi in MC
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe3n220000","32c51070a","01031v3100000000","0453503000000000"); // EG2


    //-----
    //Trigger EG1: Calo Variations
    //-----
  } else if(trainConfig == 1101)  { //EDC 13TeV MB, timing: 6 (-30+35) -> 411790106fe32220000
    cuts.AddCutHeavyMesonCalo("0008d113","411790106fe32220000","32c51070a","0103103100000000","0453503000000000"); // INT7
  } else if(trainConfig == 1102)  { //EDC 13TeV MB, MinE: c (0.65 GeV) -> 411790109fec2220000
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fec2220000","32c51070a","0103103100000000","0453503000000000"); // INT7
  } else if(trainConfig == 1103)  { //EDC 13TeV MB, MinE: c (0.8 GeV) -> 411790109fe42220000
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe42220000","32c51070a","0103103100000000","0453503000000000"); // INT7
  } else if(trainConfig == 1104)  { //EDC 13TeV MB, NCell: 0 (no cut up to 4 GeV ) -> 411790109fe30220000
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe30220000","32c51070a","0103103100000000","0453503000000000"); // INT7
  } else if(trainConfig == 1105)  { //EDC 13TeV MB, M02: 3 (0.5)  -> 411790109fe32230000
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe32230000","32c51070a","0103103100000000","0453503000000000"); // INT7
  } else if(trainConfig == 1106)  { //EDC 13TeV MB, M02: 1 (1.0) -> 411790109fe32210000
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe32210000","32c51070a","0103103100000000","0453503000000000"); // INT7
  } else if(trainConfig == 1107)  { //EDC 13TeV MB, NCell: v (NCell Cut 2, but with probability in MC to let clusters through ) -> 411790109fe3n220000
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe3n220000","32c51070a","0103103100000000","0453503000000000"); // INT7
    //-----
    //Trigger EG1: Meson Cut Variations
    //-----
  } else if(trainConfig == 1110)  { //EDC 13TeV MB, rapidity 5 (0.85)
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe32220000","32c51070a","0103503100000000","0453503000000000"); // INT7
  } else if(trainConfig == 1111)  { //EDC 13TeV MB, rapidity 6 (0.75)
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe32220000","32c51070a","0103603100000000","0453503000000000"); // INT7
  } else if(trainConfig == 1112)  { //EDC 13TeV MB, non pt dependent pi0 mass 1 (0.10-0.145)
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe32220000","32c51070a","0103103100000000","0453503000000000"); // INT7
  } else if(trainConfig == 1113)  { //EDC 13TeV MB, non pt dependent pi0 mass 2 (0.11-0.145)
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe32220000","32c51070a","0103103200000000","0453503000000000"); // INT7
  } else if(trainConfig == 1114)  { //EDC 13TeV MB, non pt dependent pi0 mass 3 (0.12-0.145)
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe32220000","32c51070a","0103103300000000","0453503000000000"); // INT7
  } else if(trainConfig == 1115)  { //EDC 13TeV MB, non pt dependent pi0 mass 4 (0.10-0.150)
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe32220000","32c51070a","0103103400000000","0453503000000000"); // INT7
  } else if(trainConfig == 1116)  { //EDC 13TeV MB, min/max pt cut o (min gamma 10gev, max pt 20gev)
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe32220000","32c51070a","01031o3100000000","0453503000000000"); // INT7
  } else if(trainConfig == 1117)  { //EDC 13TeV MB, min/max pt cut p (min gamma 10gev, max pt 25gev)
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe32220000","32c51070a","01031p3100000000","0453503000000000"); // INT7
  } else if(trainConfig == 1118)  { //EDC 13TeV MB, min/max pt cut s (no min pt, max pt 20gev)
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe32220000","32c51070a","01031s3100000000","0453503000000000"); // INT7
  } else if(trainConfig == 1119)  { //EDC 13TeV MB, min/max pt cut v (no min pt, max pt  25gev)
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe32220000","32c51070a","01031v3100000000","0453503000000000"); // INT7

    //-----
    //Trigger EG1: Primary Pion Variations
    //-----
  } else if(trainConfig == 1120)  { //EDC 13TeV MB, Ch.Pi 650MeV
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe32220000","32c51070f","0103103100000000","0453503000000000"); // INT7
  } else if(trainConfig == 1121)  { //EDC 13TeV MB, Ch.Pi 700MeV
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe32220000","32c51070g","0103103100000000","0453503000000000"); // INT7
  } else if(trainConfig == 1122)  { //EDC 13TeV MB, Ch.Pi 750MeV
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe32220000","32c51070a","0103103100000000","0453503000000000"); // INT7
  } else if(trainConfig == 1123)  { //EDC 13TeV MB, Ch.Pi 850MeV, Ch.Neu.Pi 1000MeV
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe32220000","32c51070c","0103103100000000","0453503000000000"); // INT7
  } else if(trainConfig == 1124)  { //EDC 13TeV MB, Ch.Pi 650MeV, Ch.Neu.Pi 1000MeV
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe32220000","32c51070h","0103103100000000","0453503000000000"); // INT7
  } else if(trainConfig == 1125)  { //EDC 13TeV MB, Ch.Pi 650MeV, Ch.Neu.Pi 850MeV
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe32220000","32c51070i","0103103100000000","0453503000000000"); // INT7

    //NL Changes
    //NCell 2, NL 21
  } else if(trainConfig == 1160)  { //EDC 13TeV MB, NCell: 2
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe32220000","32c51070a","0103103100000000","0453503000000000"); // INT7
  } else if(trainConfig == 1161)  { //EDC 13TeV MB, min/max pt cut s (no min pt, max pt 20gev), NCell Cut 2
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe32220000","32c51070a","01031s3100000000","0453503000000000"); // INT7
  } else if(trainConfig == 1162)  { //EDC 13TeV MB, min/max pt cut v (no min pt, max pt  25gev), NCell Cut 2
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe32220000","32c51070a","01031v3100000000","0453503000000000"); // INT7
    //NCell 2, NL 01 (new std)
  } else if(trainConfig == 1163)  { //EDC 13TeV MB, NCell: 2
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe32220000","32c51070a","0103103100000000","0453503000000000"); // INT7
  } else if(trainConfig == 1164)  { //EDC 13TeV MB, min/max pt cut s (no min pt, max pt 20gev), NCell Cut 2
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe32220000","32c51070a","01031s3100000000","0453503000000000"); // INT7
  } else if(trainConfig == 1165)  { //EDC 13TeV MB, min/max pt cut v (no min pt, max pt  25gev), NCell Cut 2
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe32220000","32c51070a","01031v3100000000","0453503000000000"); // INT7

    //NCell 0, NL 21
  } else if(trainConfig == 1170)  { //EDC 13TeV MB, NCell: 0 (no cut up to 4 GeV ) -> 411790109fe30220000
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe30220000","32c51070a","0103103100000000","0453503000000000"); // INT7
  } else if(trainConfig == 1171)  { //EDC 13TeV MB, min/max pt cut s (no min pt, max pt 20gev), NCell Cut 0
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe30220000","32c51070a","01031s3100000000","0453503000000000"); // INT7
  } else if(trainConfig == 1172)  { //EDC 13TeV MB, min/max pt cut v (no min pt, max pt  25gev), NCell Cut 0
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe30220000","32c51070a","01031v3100000000","0453503000000000"); // INT7
    //NCell 0, NL 01 (new std)
  } else if(trainConfig == 1173)  { //EDC 13TeV MB, NCell: 0 (no cut up to 4 GeV ) -> 411790109fe30220000
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe30220000","32c51070a","0103103100000000","0453503000000000"); // INT7
  } else if(trainConfig == 1174)  { //EDC 13TeV MB, min/max pt cut s (no min pt, max pt 20gev), NCell Cut 0
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe30220000","32c51070a","01031s3100000000","0453503000000000"); // INT7
  } else if(trainConfig == 1175)  { //EDC 13TeV MB, min/max pt cut v (no min pt, max pt  25gev), NCell Cut 0
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe30220000","32c51070a","01031v3100000000","0453503000000000"); // INT7

    //NCell 2 + NCell effi, NL 01
  } else if(trainConfig == 1183)  { //EDC 13TeV EG1, NCell: n NCell2 + NCell effi in MC
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe3n220000","32c51070a","0103103100000000","0453503000000000"); // EG1
  } else if(trainConfig == 1184)  { //EDC 13TeV EG1, min/max pt cut s (no min pt, max pt 20gev), NCell: n NCell2 + NCell effi in MC
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe3n220000","32c51070a","01031s3100000000","0453503000000000"); // EG1
  } else if(trainConfig == 1185)  { //EDC 13TeV EG1, min/max pt cut v (no min pt, max pt  25gev), NCell: n NCell2 + NCell effi in MC
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe3n220000","32c51070a","01031v3100000000","0453503000000000"); // EG1
  } else if(trainConfig == 1200)  { //EMCal + DCal INT7 Ch Pi MassCut function o -> 0<Mass<470 & 510<Mass<650
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe30220000","32c51070o","0103103100000000","0453503000000000"); //INT7 MassCut o
  } else if(trainConfig == 1201)  { //EMCal + DCal INT7 Ch Pi MassCut function p -> 0<Mass<470 & 510<Mass<700
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe30220000","32c51070p","0103103100000000","0453503000000000"); //INT7 MassCut p
  } else if(trainConfig == 1202)  { //EMCal + DCal INT7 Ch Pi MassCut function q -> 0<Mass<470 & 510<Mass<850
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe30220000","32c51070q","0103103100000000","0453503000000000"); //INT7 MassCut q
  } else if(trainConfig == 1203)  { //EMCal + DCal INT7 Ch Pi MassCut function o -> 0<Mass<470 & 510<Mass<650, ToF 3 -> -3, 5 Sigma
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe30220000","32c51073o","0103103100000000","0453503000000000"); //INT7 MassCut o, ToF 3
  } else if(trainConfig == 1204)  { //EMCal + DCal INT7 Ch Pi MassCut function p -> 0<Mass<470 & 510<Mass<700, ToF 3 -> -3, 5 Sigma
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe30220000","32c51073p","0103103100000000","0453503000000000"); //INT7 MassCut p, ToF 3
  } else if(trainConfig == 1205)  { //EMCal + DCal INT7 Ch Pi MassCut function q -> 0<Mass<470 & 510<Mass<850, ToF 3 -> -3, 5 Sigma
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe30220000","32c51073q","0103103100000000","0453503000000000"); //INT7 MassCut q, ToF 3
  } else if(trainConfig == 1210)  { //EMCal + DCal EG2 Ch Pi MassCut function o -> 0<Mass<470 & 510<Mass<650
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe30220000","32c51070o","01031v3100000000","0453503000000000"); //EG2 MassCut o
  } else if(trainConfig == 1211)  { //EMCal + DCal EG2 Ch Pi MassCut function p -> 0<Mass<470 & 510<Mass<700
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe30220000","32c51070p","01031v3100000000","0453503000000000"); //EG2 MassCut p
  } else if(trainConfig == 1212)  { //EMCal + DCal EG2 Ch Pi MassCut function q -> 0<Mass<470 & 510<Mass<850
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe30220000","32c51070q","01031v3100000000","0453503000000000"); //EG2 MassCut q
  } else if(trainConfig == 1213)  { //EMCal + DCal EG2 Ch Pi MassCut function o -> 0<Mass<470 & 510<Mass<650, ToF 3 -> -3, 5 Sigma
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe30220000","32c51073o","01031v3100000000","0453503000000000"); //EG2 MassCut o, ToF 3
  } else if(trainConfig == 1214)  { //EMCal + DCal EG2 Ch Pi MassCut function p -> 0<Mass<470 & 510<Mass<700, ToF 3 -> -3, 5 Sigma
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe30220000","32c51073p","01031v3100000000","0453503000000000"); //EG2 MassCut p, ToF 3
  } else if(trainConfig == 1215)  { //EMCal + DCal EG2 Ch Pi MassCut function q -> 0<Mass<470 & 510<Mass<850, ToF 3 -> -3, 5 Sigma
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe30220000","32c51073q","01031v3100000000","0453503000000000"); //EG2 MassCut q, ToF 3
  } else if(trainConfig == 1220)  { //EMCal + DCal EG1 Ch Pi MassCut function o -> 0<Mass<470 & 510<Mass<650
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe30220000","32c51070o","01031v3100000000","0453503000000000"); //EG1 MassCut o
  } else if(trainConfig == 1221)  { //EMCal + DCal EG1 Ch Pi MassCut function p -> 0<Mass<470 & 510<Mass<700
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe30220000","32c51070p","01031v3100000000","0453503000000000"); //EG1 MassCut p
  } else if(trainConfig == 1222)  { //EMCal + DCal EG1 Ch Pi MassCut function q -> 0<Mass<470 & 510<Mass<850
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe30220000","32c51070q","01031v3100000000","0453503000000000"); //EG1 MassCut q
  } else if(trainConfig == 1223)  { //EMCal + DCal EG1 Ch Pi MassCut function o -> 0<Mass<470 & 510<Mass<650, ToF 3 -> -3, 5 Sigma
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe30220000","32c51073o","01031v3100000000","0453503000000000"); //EG1 MassCut o, ToF 3
  } else if(trainConfig == 1224)  { //EMCal + DCal EG1 Ch Pi MassCut function p -> 0<Mass<470 & 510<Mass<700, ToF 3 -> -3, 5 Sigma
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe30220000","32c51073p","01031v3100000000","0453503000000000"); //EG1 MassCut p, ToF 3
  } else if(trainConfig == 1225)  { //EMCal + DCal EG1 Ch Pi MassCut function q -> 0<Mass<470 & 510<Mass<850, ToF 3 -> -3, 5 Sigma
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe30220000","32c51073q","01031v3100000000","0453503000000000"); //EG1 MassCut q, ToF 3

  } else if(trainConfig == 1230)  { //EMCal + DCal INT7, N.Pi cut var. Selection Window, Std 1 -> 2 sigma
    //                         00010113   411790109fe30220000   32c51070a   0103103100000000   0453503000000000
    //                                                                             |
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe30220000","32c51070a","0103103100000000","0453503000000000"); // INT7, 2.0 sigma
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe30220000","32c51070a","0103103x00000000","0453503000000000"); // INT7, 3.0 sigma
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe30220000","32c51070a","0103103w00000000","0453503000000000"); // INT7, 4.0 sigma
  } else if(trainConfig == 1231)  { //EMCal + DCal EG2, N.Pi cut var. Selection Window, Std 1 -> 2 sigma
    //                         0008e013   411790109fe30220000   32c51070a   01031v3100000000   0453503000000000
    //                                                                             |
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe30220000","32c51070a","01031v3100000000","0453503000000000"); // EG2, 2.0 sigma
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe30220000","32c51070a","01031v3x00000000","0453503000000000"); // EG2, 3.0 sigma
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe30220000","32c51070a","01031v3w00000000","0453503000000000"); // EG2, 4.0 sigma
  } else if(trainConfig == 1232)  { //EMCal + DCal EG1, N.Pi cut var. Selection Window, Std 1 -> 2 sigma
    //                         0008d113   411790109fe30220000   32c51070a   01031v3100000000   0453503000000000
    //                                                                             |
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe30220000","32c51070a","01031v3100000000","0453503000000000"); // EG1, 2.0 sigma
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe30220000","32c51070a","01031v3x00000000","0453503000000000"); // EG1, 3.0 sigma
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe30220000","32c51070a","01031v3w00000000","0453503000000000"); // EG1, 4.0 sigma
  } else if(trainConfig == 1233)  { //PHOS INT7, N.Pi cut var. Selection Window, Std 3 -> 2 sigma
    //                                                                      0103103300000000
    //                                                                             |
    cuts.AddCutHeavyMesonCalo("00010113","24466190sa01cc00000","32c51070a","0103103300000000","0453503000000000"); // INT7, 2.0 sigma
    cuts.AddCutHeavyMesonCalo("00010113","24466190sa01cc00000","32c51070a","0103103s00000000","0453503000000000"); // INT7, 3.0 sigma
    cuts.AddCutHeavyMesonCalo("00010113","24466190sa01cc00000","32c51070a","0103103q00000000","0453503000000000"); // INT7, 4.0 sigma

  } else if(trainConfig == 1234)  { //EMCal + DCal EG2, N.Pi cut var. Selection Window, Std 1 -> 2 sigma, max Pi0Pt 20GeV
    //                         0008e013   411790109fe30220000   32c51070a   01031v3100000000   0453503000000000
    //                                                                           | |
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe30220000","32c51070a","01031s3100000000","0453503000000000"); // EG2, 2.0 sigma
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe30220000","32c51070a","01031s3x00000000","0453503000000000"); // EG2, 3.0 sigma
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe30220000","32c51070a","01031s3w00000000","0453503000000000"); // EG2, 4.0 sigma
  } else if(trainConfig == 1235)  { //EMCal + DCal EG1, N.Pi cut var. Selection Window, Std 1 -> 2 sigma, max Pi0Pt 20GeV
    //                         0008d113   411790109fe30220000   32c51070a   01031v3100000000   0453503000000000
    //                                                                           | |
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe30220000","32c51070a","01031s3100000000","0453503000000000"); // EG1, 2.0 sigma
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe30220000","32c51070a","01031s3x00000000","0453503000000000"); // EG1, 3.0 sigma
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe30220000","32c51070a","01031s3w00000000","0453503000000000"); // EG1, 4.0 sigma


  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // Cuts for the omega pp 5TeV analysis
  // ++++++++++++++++++++++++++++++++++++++++  Minimum Bias  ++++++++++++++++++++++++++++++++++++++++++++++
  }else if (trainConfig == 1500){ // Standard 5 TeV omega INT7 cutstring
    cuts.AddCutHeavyMesonCalo("00010113","411790105fe30220000","32c51079a","0000003100000000","0400503000000000");
  }else if (trainConfig == 1501){ // Non linearity variations
    cuts.AddCutHeavyMesonCalo("00010113","411799705fe30220000","32c51079a","0000003100000000","0400503000000000"); // 97: CRF
    cuts.AddCutHeavyMesonCalo("00010113","411799805fe30220000","32c51079a","0000003100000000","0400503000000000"); // 98: CCRF
  }else if (trainConfig == 1502){ // Cluster timing variations
    cuts.AddCutHeavyMesonCalo("00010113","411790106fe30220000","32c51079a","0000003100000000","0400503000000000"); // 6: -30 - 35
    cuts.AddCutHeavyMesonCalo("00010113","411790107fe30220000","32c51079a","0000003100000000","0400503000000000"); // 7: -30 - 30
    cuts.AddCutHeavyMesonCalo("00010113","411790108fe30220000","32c51079a","0000003100000000","0400503000000000"); // 8: -20 - 30
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe30220000","32c51079a","0000003100000000","0400503000000000"); // 9: -20 - 25
    cuts.AddCutHeavyMesonCalo("00010113","41179010afe30220000","32c51079a","0000003100000000","0400503000000000"); // a: -12.5 - 13
  }else if (trainConfig == 1503){ // Track matching variations
    cuts.AddCutHeavyMesonCalo("00010113","411790105ce30220000","32c51079a","0000003100000000","0400503000000000"); // c: No E/p cut
    cuts.AddCutHeavyMesonCalo("00010113","411790105ee30220000","32c51079a","0000003100000000","0400503000000000"); // e: E/p < 2
    cuts.AddCutHeavyMesonCalo("00010113","411790105ge30220000","32c51079a","0000003100000000","0400503000000000"); // g: E/p < 1.5
    cuts.AddCutHeavyMesonCalo("00010113","411790105ne30220000","32c51079a","0000003100000000","0400503000000000"); // n: E/p < 1.75 (standard), but eta and phi varied
    cuts.AddCutHeavyMesonCalo("00010113","411790105oe30220000","32c51079a","0000003100000000","0400503000000000"); // o: E/p < 1.75 (standard), but eta and phi varied
  }else if (trainConfig == 1504){ // Exotic cluster variations
    cuts.AddCutHeavyMesonCalo("00010113","411790105f030220000","32c51079a","0000003100000000","0400503000000000"); // 0: No exotics cut
    cuts.AddCutHeavyMesonCalo("00010113","411790105fb30220000","32c51079a","0000003100000000","0400503000000000"); // b: fExoticEnergyFracCluster = 0.95
    cuts.AddCutHeavyMesonCalo("00010113","411790105fi30220000","32c51079a","0000003100000000","0400503000000000"); // i: fExoticMinEnergyCell = 3
  }else if (trainConfig == 1505){ // Min cluster energy variations
    cuts.AddCutHeavyMesonCalo("00010113","411790105fe20220000","32c51079a","0000003100000000","0400503000000000"); // 2: E > 0.6
    cuts.AddCutHeavyMesonCalo("00010113","411790105fe40220000","32c51079a","0000003100000000","0400503000000000"); // 4: E > 0.8
  }else if (trainConfig == 1506){ // NCell variation
    cuts.AddCutHeavyMesonCalo("00010113","411790105fe3n220000","32c51079a","0000003100000000","0400503000000000");
  }else if (trainConfig == 1507){ // M02 variations
    cuts.AddCutHeavyMesonCalo("00010113","411790105fe30210000","32c51079a","0000003100000000","0400503000000000"); // 1: M02 < 1
    cuts.AddCutHeavyMesonCalo("00010113","411790105fe30230000","32c51079a","0000003100000000","0400503000000000"); // 3: M02 < 0.5
  }else if (trainConfig == 1508){ // pi0 mass selection window variations
    cuts.AddCutHeavyMesonCalo("00010113","411790105fe30220000","32c51079a","0000003u00000000","0400503000000000"); // u: 1.5 sigma
    cuts.AddCutHeavyMesonCalo("00010113","411790105fe30220000","32c51079a","0000003v00000000","0400503000000000"); // v: 2.5 sigma
    cuts.AddCutHeavyMesonCalo("00010113","411790105fe30220000","32c51079a","0000003x00000000","0400503000000000"); // x: 3 sigma
    cuts.AddCutHeavyMesonCalo("00010113","411790105fe30220000","32c51079a","0000003r00000000","0400503000000000"); // r: 3.5 sigma, no gamma selection
    cuts.AddCutHeavyMesonCalo("00010113","411790105fe30220000","32c51079a","0000003w00000000","0400503000000000"); // w: 4 sigma
  }else if (trainConfig == 1509){ // pi0 asymmetry variations
    cuts.AddCutHeavyMesonCalo("00010113","411790105fe30220000","32c51079a","0000005100000000","0400503000000000"); // 5: alpha < 0.75
    cuts.AddCutHeavyMesonCalo("00010113","411790105fe30220000","32c51079a","0000006100000000","0400503000000000"); // 6: alpha < 0.8
    cuts.AddCutHeavyMesonCalo("00010113","411790105fe30220000","32c51079a","0000007100000000","0400503000000000"); // 7: alpha < 0.85
  }else if (trainConfig == 1510){ // ITS cluster requirement variation
    cuts.AddCutHeavyMesonCalo("00010113","411790105fe30220000","34c51079a","0000003100000000","0400503000000000"); // 4: min 3 ITS cluster
  }else if (trainConfig == 1511){ // TPC cluster requirement variation
    cuts.AddCutHeavyMesonCalo("00010113","411790105fe30220000","32e51079a","0000003100000000","0400503000000000"); // e: No shared clusters
  }else if (trainConfig == 1512){ // Charged pion DCA
    cuts.AddCutHeavyMesonCalo("00010113","411790105fe30220000","32c01079a","0000003100000000","0400503000000000"); // 0: No cut
    cuts.AddCutHeavyMesonCalo("00010113","411790105fe30220000","32c31079a","0000003100000000","0400503000000000"); // 5: xy<2.4, z<3.2
    cuts.AddCutHeavyMesonCalo("00010113","411790105fe30220000","32c61079a","0000003100000000","0400503000000000"); // 6: z,xy < 0.5 (very tight)
  }else if (trainConfig == 1513){ // TOF requirement
    cuts.AddCutHeavyMesonCalo("00010113","411790105fe30220000","32c51070a","0000003100000000","0400503000000000"); // 0: No TOF
    cuts.AddCutHeavyMesonCalo("00010113","411790105fe30220000","32c51076a","0000003100000000","0400503000000000"); // 6: stricter Kp rejection
    cuts.AddCutHeavyMesonCalo("00010113","411790105fe30220000","32c51072a","0000003100000000","0400503000000000"); // 2: Pion selection
  }else if (trainConfig == 1514){ // min pT
    cuts.AddCutHeavyMesonCalo("00010113","411790105fe30220000","32c50079a","0000003100000000","0400503000000000"); // 0: pT > 0.075 GeV
    cuts.AddCutHeavyMesonCalo("00010113","411790105fe30220000","32c52079a","0000003100000000","0400503000000000"); // 2: pT > 0.125 GeV
    cuts.AddCutHeavyMesonCalo("00010113","411790105fe30220000","32c53079a","0000003100000000","0400503000000000"); // 3: pT > 0.15 GeV
  }else if (trainConfig == 1515){ // TPC dEdx sigma
    cuts.AddCutHeavyMesonCalo("00010113","411790105fe30220000","32c510b9a","0000003100000000","0400503000000000"); // b: -2,2
    cuts.AddCutHeavyMesonCalo("00010113","411790105fe30220000","32c51099a","0000003100000000","0400503000000000"); // 9: -2.5,2.5
    cuts.AddCutHeavyMesonCalo("00010113","411790105fe30220000","32c510a9a","0000003100000000","0400503000000000"); // a: -3.5,3.5
    cuts.AddCutHeavyMesonCalo("00010113","411790105fe30220000","32c51059a","0000003100000000","0400503000000000"); // 5: -4,4
    cuts.AddCutHeavyMesonCalo("00010113","411790105fe30220000","32c51039a","0000003100000000","0400503000000000"); // 3: -5,5
  }else if (trainConfig == 1516){ // PiPlPiMi Mass
    cuts.AddCutHeavyMesonCalo("00010113","411790105fe30220000","32c51079r","0000003100000000","0400503000000000"); // r: 0.8 GeV/c^2
    cuts.AddCutHeavyMesonCalo("00010113","411790105fe30220000","32c51079s","0000003100000000","0400503000000000"); // s: 0.825 GeV/c^2
    cuts.AddCutHeavyMesonCalo("00010113","411790105fe30220000","32c51079t","0000003100000000","0400503000000000"); // t: 0.875 GeV/c^2
    cuts.AddCutHeavyMesonCalo("00010113","411790105fe30220000","32c51079u","0000003100000000","0400503000000000"); // u: 0.9 GeV/c^2






    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // EMC pp 13 TeV Fitting, Systematics
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //Standard Cuts of Pi0 Analysis: ("00010113","411790109fe30220000","0r631031000000d0")
    //                                                                   ||           |
    //                                                                 "0103103x00000000"
    //MesonCut r63==Background->ignored, d==OpeningAngle for Background->ignored =>0103103x00000000
  } else if(trainConfig == 2000)  { //EMCal + DCal INT7 Standard
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe30220000","32l51070a","0103103x00000000","0453503000000000"); // INT7 Standard
    //-----
    //INT7: Event Variations
    //-----
    //Std: 00010113
  } else if(trainConfig == 2001)  { //EMCal + DCal INT7, Event cut var. Remove Pileup, Std 1-> True
    //                         00010113   411790109fe30220000   32l51070a   0103103x00000000   0453503000000000
    //                              |
    cuts.AddCutHeavyMesonCalo("00010013","411790109fe30220000","32l51070a","0103103x00000000","0453503000000000"); // INT7 Pileup not removed
    //-----
    //INT7: Calo Variations
    //-----
    //Std: 411790109fe30220000
  } else if(trainConfig == 2101)  { //EMCal + DCal INT7, Calo cut var. NonLins, Std 01
    //                         00010113   411790109fe30220000   32l51070a   0103103x00000000   0453503000000000
    //                                         ||
    cuts.AddCutHeavyMesonCalo("00010113","411799609fe30220000","32l51070a","0103103x00000000","0453503000000000"); // INT7 no FT applied
    cuts.AddCutHeavyMesonCalo("00010113","411799709fe30220000","32l51070a","0103103x00000000","0453503000000000"); // INT7 EMC fine tuning applied
    cuts.AddCutHeavyMesonCalo("00010113","411799809fe30220000","32l51070a","0103103x00000000","0453503000000000"); // INT7 PCM-EMC fine tuning applied
  } else if(trainConfig == 2102)  { //EMCal + DCal INT7, Calo cut var. time, Std 9 -> -20+25
    //                         00010113   411790109fe30220000   32l51070a   0103103x00000000   0453503000000000
    //                                            |
    cuts.AddCutHeavyMesonCalo("00010113","411790105fe30220000","32l51070a","0103103x00000000","0453503000000000"); // INT7 time -50+50
    cuts.AddCutHeavyMesonCalo("00010113","411790106fe30220000","32l51070a","0103103x00000000","0453503000000000"); // INT7 time -30+35
    cuts.AddCutHeavyMesonCalo("00010113","411790108fe30220000","32l51070a","0103103x00000000","0453503000000000"); // INT7 time -20+30
    cuts.AddCutHeavyMesonCalo("00010113","41179010afe30220000","32l51070a","0103103x00000000","0453503000000000"); // INT7 time -12.5+13
  } else if(trainConfig == 2103)  { //EMCal + DCal INT7, Calo cut var. energy, Std 3 -> 0.7 GeV
    //                         00010113   411790109fe30220000   32l51070a   0103103x00000000   0453503000000000
    //                                               |
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe20220000","32l51070a","0103103x00000000","0453503000000000"); // INT7 energy 0.6 GeV
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe40220000","32l51070a","0103103x00000000","0453503000000000"); // INT7 energy 0.8 GeV
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe50220000","32l51070a","0103103x00000000","0453503000000000"); // INT7 energy 0.9 GeV // only meaningfull at higher pTs
  } else if(trainConfig == 2104)  { //EMCal + DCal INT7, Calo cut var. NCell, Std 0 -> Turned Off until 4GeV; then min 2 Cells
    //                         00010113   411790109fe30220000   32l51070a   0103103x00000000   0453503000000000
    //                                                |
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe32220000","32l51070a","0103103x00000000","0453503000000000"); // INT7 NCells 2
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe3n220000","32l51070a","0103103x00000000","0453503000000000"); // INT7 NCells 2 var (PCM-EMCal tagging corr)
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe3r220000","32l51070a","0103103x00000000","0453503000000000"); // INT7 NCells 2 var (EMCal tagging corr)
  } else if(trainConfig == 2105)  { //EMCal + DCal INT7, Calo cut var. max M02, 2 -> INT7 M02 0.7
    //                         00010113   411790109fe30220000   32l51070a   0103103x00000000   0453503000000000
    //                                                  |
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe30230000","32l51070a","0103103x00000000","0453503000000000"); // INT7 M02 0.5
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe30210000","32l51070a","0103103x00000000","0453503000000000"); // INT7 M02 1.0
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe302k0000","32l51070a","0103103x00000000","0453503000000000"); // INT7 M02 E dep
  } else if(trainConfig == 2106)  { //EMCal + DCal INT7, Calo cut var. TM, Std f,
    //  [1] + 1 / pow(x + pow(1 / ([0] - [1]), 1 / [2]), [2]; Eta (0.04, 0.010, 2.5); Phi (0.09, 0.015, 2.); EoverP 1.75
    //                         00010113   411790109fe30220000   32l51070a   0103103x00000000   0453503000000000
    //                                             |
    cuts.AddCutHeavyMesonCalo("00010113","411790109ee30220000","32l51070a","0103103x00000000","0453503000000000"); // INT7 TM var e: EoverP 2.00
    cuts.AddCutHeavyMesonCalo("00010113","411790109ge30220000","32l51070a","0103103x00000000","0453503000000000"); // INT7 TM var g: EoverP 1.5
    cuts.AddCutHeavyMesonCalo("00010113","4117901097e30220000","32l51070a","0103103x00000000","0453503000000000"); // INT7 TM var 7: no EoverP
    cuts.AddCutHeavyMesonCalo("00010113","411790109ne30220000","32l51070a","0103103x00000000","0453503000000000"); // INT7 TM var n: Eta (0.035, 0.010, 2.5); Phi (0.085, 0.015, 2.)
    cuts.AddCutHeavyMesonCalo("00010113","411790109oe30220000","32l51070a","0103103x00000000","0453503000000000"); // INT7 TM var, Eta (0.045, 0.010, 2.5); Phi (0.095, 0.015, 2.)
  } else if(trainConfig == 2107)  { //EMCal + DCal INT7, Calo cut var. Exotics, Std e, active F+ < 0.97
    //                         00010113   411790109fe30220000   32l51070a   0103103x00000000   0453503000000000
    //                                              |
    cuts.AddCutHeavyMesonCalo("00010113","411790109f030220000","32l51070a","0103103x00000000","0453503000000000"); // INT7 no exotics cut
    cuts.AddCutHeavyMesonCalo("00010113","411790109fb30220000","32l51070a","0103103x00000000","0453503000000000"); // INT7 F+ < 0.95
    //-----
    //INT7: Primary Pion / Charged Pion (Pi+ Pi-) Variations
    //-----
    //Std: 32l51070a
  } else if(trainConfig == 2201)  { //EMCal + DCal INT7, Ch.Pi cut var. Ch.Pi ITS Requirement, Std 2 -> first or second SPD cluster required
    //                         00010113   411790109fe30220000   32l51070a   0103103x00000000   0453503000000000
    //                                                           |
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe30220000","34l51070a","0103103x00000000","0453503000000000"); // INT7, Ch.Pi ITS, first or second SPD cluster required, min number of ITS clusters = 3
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe30220000","35l51070a","0103103x00000000","0453503000000000"); // INT7, Ch.Pi ITS, first or second SPD cluster required, min number of ITS clusters = 4
  } else if(trainConfig == 2202)  { //EMCal + DCal INT7, Ch.Pi cut var. Ch.Pi Cls TPC, Std l -> max shared clusters 10 + min TPC PID clusters 50
    //                         00010113   411790109fe30220000   32l51070a   0103103x00000000   0453503000000000
    //                                                            |
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe30220000","32e51070a","0103103x00000000","0453503000000000"); // INT7, Ch.Pi, MinClsTPC 80. + Refit MaxSharedClsTPCFrac=0.
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe30220000","32k51070a","0103103x00000000","0453503000000000"); // INT7 max shared clusters 10
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe30220000","32c51070a","0103103x00000000","0453503000000000"); // INT7 c -> MinClsTPC 80. + Refit
  } else if(trainConfig == 2203)  { //EMCal + DCal INT7, Ch.Pi cut var. Ch.Pi pT, Std 1 -> pt>0.1
    //                         00010113   411790109fe30220000   32l51070a   0103103x00000000   0453503000000000
    //                                                              |
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe30220000","32l50070a","0103103x00000000","0453503000000000"); // INT7, Ch.Pi pt>0.075
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe30220000","32l52070a","0103103x00000000","0453503000000000"); // INT7, Ch.Pi pt>0.125
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe30220000","32l53070a","0103103x00000000","0453503000000000"); // INT7, Ch.Pi pt>0.15
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe30220000","32l54070a","0103103x00000000","0453503000000000"); // INT7, Ch.Pi pt>0.4
  } else if(trainConfig == 2204)  { //EMCal + DCal INT7, Ch.Pi cut var. Ch.Pi TPC dEdx Low, Std 7 -> -3,3
    //                         00010113   411790109fe30220000   32l51070a   0103103x00000000   0453503000000000
    //                                                                |
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe30220000","32l510c0a","0103103x00000000","0453503000000000"); // INT7, Ch.Pi -3.5,3.0
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe30220000","32l510d0a","0103103x00000000","0453503000000000"); // INT7, Ch.Pi -3.25,3.0
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe30220000","32l510e0a","0103103x00000000","0453503000000000"); // INT7, Ch.Pi -2.75,3.0
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe30220000","32l510f0a","0103103x00000000","0453503000000000"); // INT7, Ch.Pi -2.5,3.0
  } else if(trainConfig == 2205)  { //EMCal + DCal INT7, Ch.Pi cut var. Ch.Pi TPC dEdx High, Std 7 -> -3,3
    //                         00010113   411790109fe30220000   32l51070a   0103103x00000000   0453503000000000
    //                                                                |
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe30220000","32l510g0a","0103103x00000000","0453503000000000"); // INT7, Ch.Pi -3.0,2.5
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe30220000","32l510h0a","0103103x00000000","0453503000000000"); // INT7, Ch.Pi -3.0,2.75
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe30220000","32l510i0a","0103103x00000000","0453503000000000"); // INT7, Ch.Pi -3.0,3.25
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe30220000","32l510j0a","0103103x00000000","0453503000000000"); // INT7, Ch.Pi -3.0,3.5
  } else if(trainConfig == 2206)  { //EMCal + DCal INT7, Ch.Pi cut var. Ch.Pi Mass, Std a -> Ch.Pi<850MeV
    //                         00010113   411790109fe30220000   32l51070a   0103103x00000000   0453503000000000
    //                                                                  |
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe30220000","32l51070r","0103103x00000000","0453503000000000"); // INT7, Ch.Pi<800MeV
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe30220000","32l51070s","0103103x00000000","0453503000000000"); // INT7, Ch.Pi<825MeV
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe30220000","32l51070t","0103103x00000000","0453503000000000"); // INT7, Ch.Pi<875MeV
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe30220000","32l51070u","0103103x00000000","0453503000000000"); // INT7, Ch.Pi<900MeV
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe30220000","32l51070n","0103103x00000000","0453503000000000"); // INT7, Ch.Pi<1000MeV
  } else if(trainConfig == 2207)  { //EMCal + DCal INT7 shared cluster + TPC pid clusters
    // WARNING: this variation has NO cut on nSigma and is used for general dEdx plot with only track cuts
    // WARNING: this variation has NO cut on pion mass
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe30220000","32c510000","0103103x00000000","0453503000000000"); // INT7 default cut
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe30220000","32k510000","0103103x00000000","0453503000000000"); // INT7 max shared clusters 10
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe30220000","32l510000","0103103x00000000","0453503000000000"); // INT7 max shared clusters 10 + min TPC PID clusters 50
  } else if(trainConfig == 2208)  { //EMCal + DCal INT7 shared cluster + TPC pid clusters
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe30220000","32k51070a","0103103x00000000","0453503000000000"); // INT7 max shared clusters 10
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe30220000","32l51070a","0103103x00000000","0453503000000000"); // INT7 max shared clusters 10 + min TPC PID clusters 50
    //-----
    //INT7: Neutral Meson (Pi0) Cut Variations
    //-----
    //Std: 0103103x00000000
  } else if(trainConfig == 2302)  { //EMCal + DCal INT7, N.Pi cut var. rapidity, Std 1 -> -0.8, 0.8
    //                         00010113   411790109fe30220000   32l51070a   0103103x00000000   0453503000000000
    //                                                                          |
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe30220000","32l51070a","0103503x00000000","0453503000000000"); // INT7, N.Pi rap. -0.85, 0.85
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe30220000","32l51070a","0103603x00000000","0453503000000000"); // INT7, N.Pi rap. -0.75, 0.75
  } else if(trainConfig == 2304)  { //EMCal + DCal INT7, N.Pi cut var. alpha, Std 3 -> 0.0-1.0
    //                         00010113   411790109fe30220000   32l51070a   0103103x00000000   0453503000000000
    //                                                                            |
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe30220000","32l51070a","0103105x00000000","0453503000000000"); // INT7 alpha 0-0.75
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe30220000","32l51070a","0103106x00000000","0453503000000000"); // INT7 alpha 0-0.8
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe30220000","32l51070a","0103107x00000000","0453503000000000"); // INT7 alpha 0-0.85
  } else if(trainConfig == 2305)  { //EMCal + DCal INT7, N.Pi cut var. Selection Window, Std x -> 3 sigma
    //                         00010113   411790109fe30220000   32l51070a   0103103x00000000   0453503000000000
    //                                                                             |
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe30220000","32l51070a","0103103100000000","0453503000000000"); // INT7, 2.0 sigma
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe30220000","32l51070a","0103103w00000000","0453503000000000"); // INT7, 4.0 sigma
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe30220000","32l51070a","0103103v00000000","0453503000000000"); // INT7, 2.5 sigma
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe30220000","32l51070a","0103103r00000000","0453503000000000"); // INT7, 3.5 sigma
  } else if(trainConfig == 2306)  { //EMCal + DCal INT7, N.Pi cut var. open. angle, Std 0 -> off
    //                         00010113   411790109fe30220000   32l51070a   0103103x00000000   0453503000000000
    //                                                                                    |
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe30220000","32l51070a","0103103x000000d0","0453503000000000"); // INT7 Op. Ang. var 1 cell dist + 0.017
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe30220000","32l51070a","0103103x000000b0","0453503000000000"); // INT7 Op. Ang. var 1 cell dist + 0.0152
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe30220000","32l51070a","0103103x000000g0","0453503000000000"); // INT7 Op. Ang. var 1 cell dist + 0.0202
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe30220000","32l51070a","0103103x000000a0","0453503000000000"); // INT7 Op. Ang. var 1 cell dist + 0

    //-----
    //INT7: Omega Meson Cut Variations
    //-----
    //Std: 0453503000000000
  } else if(trainConfig == 2401)  { //EMCal + DCal INT7, Omega cut var. Background Scheme, Std 4 -> off
    //                         00010113   411790109fe30220000   32l51070a   0103103x00000000   0453503000000000
    //                                                                                          |
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe30220000","32l51070a","0103103x00000000","0153503000000000"); // INT7, Om Event Mixing
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe30220000","32l51070a","0103103x00000000","0a53503000000000"); // INT7, Om LikeSignMixing
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe30220000","32l51070a","0103103x00000000","0d53503000000000"); // INT7, Om SideBandMixing
  } else if(trainConfig == 2402)  { //EMCal + DCal INT7, Omega cut var. rapidity, Std 5 -> -0.85, 0.85
    //                         00010113   411790109fe30220000   32l51070a   0103103x00000000   0453503000000000
    //                                                                                             |
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe30220000","32l51070a","0103103x00000000","0453103000000000"); // INT7, Om rap. -0.8, 0.8
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe30220000","32l51070a","0103103x00000000","0453603000000000"); // INT7, Om rap. -0.75, 0.75
  } else if(trainConfig == 2404)  { //EMCal + DCal INT7, Omega cut var. alpha, Std 3 -> 0.0-1.0
    //                         00010113   411790109fe30220000   32l51070a   0103103x00000000   0453503000000000
    //                                                                                               |
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe30220000","32l51070a","0103103x00000000","0453505000000000"); // INT7 alpha 0-0.75
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe30220000","32l51070a","0103103x00000000","0453506000000000"); // INT7 alpha 0-0.8
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe30220000","32l51070a","0103103x00000000","0453507000000000"); // INT7 alpha 0-0.85
  } else if(trainConfig == 2410)  { //EMCal + DCal INT7, Omega cut var. Background Scheme single cfg, Std 4 -> off, EventMixing
    //                         00010113   411790109fe30220000   32l51070a   0103103x00000000   0453503000000000
    //                                                                                          |
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe30220000","32l51070a","0103103x00000000","0153503000000000"); // INT7, Om Event Mixing
  } else if(trainConfig == 2411)  { //EMCal + DCal INT7, Omega cut var. Background Scheme single cfg, Std 4 -> off, LikeSignMixing
    //                         00010113   411790109fe30220000   32l51070a   0103103x00000000   0453503000000000
    //                                                                                          |
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe30220000","32l51070a","0103103x00000000","0a53503000000000"); // INT7, Om LikeSignMixing
  } else if(trainConfig == 2412)  { //EMCal + DCal INT7, Omega cut var. Background Scheme single cfg, Std 4 -> off, SideBandMixing
    //                         00010113   411790109fe30220000   32l51070a   0103103x00000000   0453503000000000
    //                                                                                          |
    cuts.AddCutHeavyMesonCalo("00010113","411790109fe30220000","32l51070a","0103103x00000000","0d53503000000000"); // INT7, Om SideBandMixing
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //Standard Cuts of Pi0 Analysis: ("0008e113","411790109fe30220000","0r631031000000d0")
    //                                                                   ||  |        |
    //                                                                 "01031v3x00000000"
    //MesonCut r63==Background->ignored, d==OpeningAngle for Background->ignored, Pi0 MaxPt Cut v->ignored =>0453503000000000
  } else if(trainConfig == 3000)  { //EMCal + DCal EG2 Standard
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe30220000","32l51070a","01031v3x00000000","0453503000000000"); // EG2 Standard
    //-----
    //EG2: Event Variations
    //-----
    //Std: 0008e113
  } else if(trainConfig == 3001)  { //EMCal + DCal EG2, Event cut var. Remove Pileup, Std 1-> True
    //                         0008e013   411790109fe30220000   32l51070a   01031v3x00000000   0453503000000000
    //                              |
    cuts.AddCutHeavyMesonCalo("0008e013","411790109fe30220000","32l51070a","01031v3x00000000","0453503000000000"); // INT7 Pileup not removed
    //-----
    //EG2: Calo Variations
    //-----
    //Std: 411790109fe30220000
  } else if(trainConfig == 3101)  { //EMCal + DCal EG2, Calo cut var. NonLins, Std 01
    //                         0008e013   411790109fe30220000   32l51070a   01031v3x00000000   0453503000000000
    //                                         ||
    cuts.AddCutHeavyMesonCalo("0008e113","411799609fe30220000","32l51070a","01031v3x00000000","0453503000000000"); // EG2 no FT applied
    cuts.AddCutHeavyMesonCalo("0008e113","411799709fe30220000","32l51070a","01031v3x00000000","0453503000000000"); // EG2 EMC fine tuning applied
    cuts.AddCutHeavyMesonCalo("0008e113","411799809fe30220000","32l51070a","01031v3x00000000","0453503000000000"); // EG2 PCM-EMC fine tuning applied
  } else if(trainConfig == 3102)  { //EMCal + DCal EG2, Calo cut var. time, Std 9 -> -20+25
    //                         0008e013   411790109fe30220000   32l51070a   01031v3x00000000   0453503000000000
    //                                            |
    cuts.AddCutHeavyMesonCalo("0008e113","411790105fe30220000","32l51070a","01031v3x00000000","0453503000000000"); // EG2 time -50+50
    cuts.AddCutHeavyMesonCalo("0008e113","411790106fe30220000","32l51070a","01031v3x00000000","0453503000000000"); // EG2 time -30+35
    cuts.AddCutHeavyMesonCalo("0008e113","411790108fe30220000","32l51070a","01031v3x00000000","0453503000000000"); // EG2 time -20+30
    cuts.AddCutHeavyMesonCalo("0008e113","41179010afe30220000","32l51070a","01031v3x00000000","0453503000000000"); // EG2 time -12.5+13
  } else if(trainConfig == 3103)  { //EMCal + DCal EG2, Calo cut var. energy, Std 3 -> 0.7 GeV
    //                         0008e013   411790109fe30220000   32l51070a   01031v3x00000000   0453503000000000
    //                                               |
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe20220000","32l51070a","01031v3x00000000","0453503000000000"); // EG2 energy 0.6 GeV
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe40220000","32l51070a","01031v3x00000000","0453503000000000"); // EG2 energy 0.8 GeV
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe50220000","32l51070a","01031v3x00000000","0453503000000000"); // EG2 energy 0.9 GeV // only meaningfull at higher pTs
  } else if(trainConfig == 3104)  { //EMCal + DCal EG2, Calo cut var. NCell, Std 0 -> Turned Off until 4GeV; then min 2 Cells
    //                         0008e013   411790109fe30220000   32l51070a   01031v3x00000000   0453503000000000
    //                                                |
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe32220000","32l51070a","01031v3x00000000","0453503000000000"); // EG2 NCells 2
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe3n220000","32l51070a","01031v3x00000000","0453503000000000"); // EG2 NCells 2 var (PCM-EMCal tagging corr)
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe3r220000","32l51070a","01031v3x00000000","0453503000000000"); // EG2 NCells 2 var (EMCal tagging corr)
  } else if(trainConfig == 3105)  { //EMCal + DCal EG2, Calo cut var. max M02, 2 -> INT7 M02 0.7
    //                         0008e013   411790109fe30220000   32l51070a   01031v3x00000000   0453503000000000
    //                                                  |
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe30230000","32l51070a","01031v3x00000000","0453503000000000"); // EG2 M02 0.5
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe30210000","32l51070a","01031v3x00000000","0453503000000000"); // EG2 M02 1.0
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe302k0000","32l51070a","01031v3x00000000","0453503000000000"); // EG2 M02 E dep
  } else if(trainConfig == 3106)  { //EMCal + DCal EG2, Calo cut var. TM, Std f,
    //  [1] + 1 / pow(x + pow(1 / ([0] - [1]), 1 / [2]), [2]; Eta (0.04, 0.010, 2.5); Phi (0.09, 0.015, 2.); EoverP 1.75
    //                         0008e013   411790109fe30220000   32l51070a   01031v3x00000000   0453503000000000
    //                                             |
    cuts.AddCutHeavyMesonCalo("0008e113","411790109ee30220000","32l51070a","01031v3x00000000","0453503000000000"); // EG2 TM var EoverP 2.00
    cuts.AddCutHeavyMesonCalo("0008e113","411790109ge30220000","32l51070a","01031v3x00000000","0453503000000000"); // EG2 TM var EoverP 1.5
    cuts.AddCutHeavyMesonCalo("0008e113","4117901097e30220000","32l51070a","01031v3x00000000","0453503000000000"); // EG2 TM var no EoverP
    cuts.AddCutHeavyMesonCalo("0008e113","411790109ne30220000","32l51070a","01031v3x00000000","0453503000000000"); // EG2 TM var, Eta (0.035, 0.010, 2.5); Phi (0.085, 0.015, 2.)
    cuts.AddCutHeavyMesonCalo("0008e113","411790109oe30220000","32l51070a","01031v3x00000000","0453503000000000"); // EG2 TM var, Eta (0.045, 0.010, 2.5); Phi (0.095, 0.015, 2.)
  } else if(trainConfig == 3107)  { //EMCal + DCal EG2, Calo cut var. Exotics, Std e, active F+ < 0.97
    //                         0008e013   411790109fe30220000   32l51070a   01031v3x00000000   0453503000000000
    //                                              |
    cuts.AddCutHeavyMesonCalo("0008e113","411790109f030220000","32l51070a","01031v3x00000000","0453503000000000"); // EG2 no exotics cut
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fb30220000","32l51070a","01031v3x00000000","0453503000000000"); // EG2 F+ < 0.95
    //-----
    //EG2: Primary Pion / Charged Pion (Pi+ Pi-) Variations
    //-----
    //Std: 32l51070a
  } else if(trainConfig == 3201)  { //EMCal + DCal EG2, Ch.Pi cut var. Ch.Pi ITS Requirement, Std 2 -> first or second SPD cluster required
    //                         0008e013   411790109fe30220000   32l51070a   01031v3x00000000   0453503000000000
    //                                                           |
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe30220000","34l51070a","01031v3x00000000","0453503000000000"); // EG2, Ch.Pi ITS, first or second SPD cluster required, min number of ITS clusters = 3
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe30220000","35l51070a","01031v3x00000000","0453503000000000"); // EG2, Ch.Pi ITS, first or second SPD cluster required, min number of ITS clusters = 4
  } else if(trainConfig == 3202)  { //EMCal + DCal EG2, Ch.Pi cut var. Ch.Pi Cls TPC, Std l -> max shared clusters 10 + min TPC PID clusters 50
    //                         0008e013   411790109fe30220000   32l51070a   01031v3x00000000   0453503000000000
    //                                                            |
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe30220000","32e51070a","01031v3x00000000","0453503000000000"); // EG2, Ch.Pi, MinClsTPC 80. + Refit MaxSharedClsTPCFrac=0.
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe30220000","32k51070a","01031v3x00000000","0453503000000000"); // EG2, Ch.Pi, max shared clusters 10
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe30220000","32c51070a","01031v3x00000000","0453503000000000"); // EG2, Ch.Pi, c -> MinClsTPC 80. + Refit
  } else if(trainConfig == 3203)  { //EMCal + DCal EG2, Ch.Pi cut var. Ch.Pi pT, Std 1 -> pt>0.1
    //                         0008e013   411790109fe30220000   32l51070a   01031v3x00000000   0453503000000000
    //                                                              |
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe30220000","32l50070a","01031v3x00000000","0453503000000000"); // EG2, Ch.Pi pt>0.075
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe30220000","32l52070a","01031v3x00000000","0453503000000000"); // EG2, Ch.Pi pt>0.125
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe30220000","32l53070a","01031v3x00000000","0453503000000000"); // EG2, Ch.Pi pt>0.15
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe30220000","32l54070a","01031v3x00000000","0453503000000000"); // EG2, Ch.Pi pt>0.4
  } else if(trainConfig == 3204)  { //EMCal + DCal EG2, Ch.Pi cut var. Ch.Pi TPC dEdx Low, Std 7 -> -3,3
    //                         0008e013   411790109fe30220000   32l51070a   01031v3x00000000   0453503000000000
    //                                                                |
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe30220000","32l510c0a","01031v3x00000000","0453503000000000"); // EG2, Ch.Pi -3.5,3.0
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe30220000","32l510d0a","01031v3x00000000","0453503000000000"); // EG2, Ch.Pi -3.25,3.0
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe30220000","32l510e0a","01031v3x00000000","0453503000000000"); // EG2, Ch.Pi -2.75,3.0
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe30220000","32l510f0a","01031v3x00000000","0453503000000000"); // EG2, Ch.Pi -2.5,3.0
  } else if(trainConfig == 3205)  { //EMCal + DCal EG2, Ch.Pi cut var. Ch.Pi TPC dEdx High, Std 7 -> -3,3
    //                         0008e013   411790109fe30220000   32l51070a   01031v3x00000000   0453503000000000
    //                                                                |
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe30220000","32l510g0a","01031v3x00000000","0453503000000000"); // EG2, Ch.Pi -3.0,2.5
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe30220000","32l510h0a","01031v3x00000000","0453503000000000"); // EG2, Ch.Pi -3.0,2.25
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe30220000","32l510i0a","01031v3x00000000","0453503000000000"); // EG2, Ch.Pi -3.0,3.25
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe30220000","32l510j0a","01031v3x00000000","0453503000000000"); // EG2, Ch.Pi -3.0,3.5
  } else if(trainConfig == 3206)  { //EMCal + DCal EG2, Ch.Pi cut var. Ch.Pi Mass, Std a -> Ch.Pi<850MeV
    //                         0008e013   411790109fe30220000   32l51070a   01031v3x00000000   0453503000000000
    //                                                                  |
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe30220000","32l51070r","01031v3x00000000","0453503000000000"); // EG2, Ch.Pi<800MeV
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe30220000","32l51070s","01031v3x00000000","0453503000000000"); // EG2, Ch.Pi<825MeV
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe30220000","32l51070t","01031v3x00000000","0453503000000000"); // EG2, Ch.Pi<875MeV
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe30220000","32l51070u","01031v3x00000000","0453503000000000"); // EG2, Ch.Pi<900MeV
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe30220000","32l51070n","01031v3x00000000","0453503000000000"); // EG2, Ch.Pi<1000MeV
  } else if(trainConfig == 3208)  { //EMCal + DCal shared cluster + TPC pid clusters
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe30220000","32l51070a","01031v3x00000000","0453503000000000"); // EG2 Standard
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe30220000","32l51070a","01031v3x00000000","0453503000000000"); // EG2 Standard
    //-----
    //EG2: Neutral Meson (Pi0) Cut Variations
    //-----
    //Std: 01031v3x00000000
  } else if(trainConfig == 3302)  { //EMCal + DCal EG2, N.Pi cut var. rapidity, Std 1 -> -0.8, 0.8
    //                         0008e013   411790109fe30220000   32l51070a   01031v3x00000000   0453503000000000
    //                                                                          |
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe30220000","32l51070a","01035v3x00000000","0453503000000000"); // EG2, N.Pi rap. -0.85, 0.85
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe30220000","32l51070a","01036v3x00000000","0453503000000000"); // EG2, N.Pi rap. -0.75, 0.75
  } else if(trainConfig == 3303)  { //EMCal + DCal EG2, N.Pi cut var. maxMass, Std v -> 25GeV
    //                         0008e013   411790109fe30220000   32l51070a   01031v3x00000000   0453503000000000
    //                                                                           |
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe30220000","32l51070a","0103103x00000000","0453503000000000"); // EG2, N.Pi maxMass off
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe30220000","32l51070a","01031s3x00000000","0453503000000000"); // EG2, N.Pi maxMass 20GeV
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe30220000","32l51070a","01031y3x00000000","0453503000000000"); // EG2, N.Pi maxMass 22GeV
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe30220000","32l51070a","01031o3x00000000","0453503000000000"); // EG2, N.Pi maxMass 25GeV, 1 Gamma >5.GeV
  } else if(trainConfig == 3304)  { //EMCal + DCal EG2, N.Pi cut var. alpha, Std 3 -> 0.0-1.0
    //                         0008e013   411790109fe30220000   32l51070a   01031v3x00000000   0453503000000000
    //                                                                            |
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe30220000","32l51070a","01031v5x00000000","0453503000000000"); // EG2 alpha 0-0.75
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe30220000","32l51070a","01031v6x00000000","0453503000000000"); // EG2 alpha 0-0.8
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe30220000","32l51070a","01031v7x00000000","0453503000000000"); // EG2 alpha 0-0.85
  } else if(trainConfig == 3305)  { //EMCal + DCal EG2, N.Pi cut var. Selection Window, Std x -> 3 sigma
    //                         0008e013   411790109fe30220000   32l51070a   01031v3x00000000   0453503000000000
    //                                                                             |
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe30220000","32l51070a","01031v3100000000","0453503000000000"); // EG2, 2.0 sigma
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe30220000","32l51070a","01031v3w00000000","0453503000000000"); // EG2, 4.0 sigma
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe30220000","32l51070a","01031v3v00000000","0453503000000000"); // EG2, 2.5 sigma
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe30220000","32l51070a","01031v3r00000000","0453503000000000"); // EG2, 3.5 sigma
  } else if(trainConfig == 3306)  { //EMCal + DCal EG2, N.Pi cut var. open. angle, Std 0 -> off
    //                         0008e013   411790109fe30220000   32l51070a   01031v3x00000000   0453503000000000
    //                                                                                    |
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe30220000","32l51070a","01031v3x000000d0","0453503000000000"); // EG2 Op. Ang. var 1 cell dist + 0.017
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe30220000","32l51070a","01031v3x000000b0","0453503000000000"); // EG2 Op. Ang. var 1 cell dist + 0.0152
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe30220000","32l51070a","01031v3x000000g0","0453503000000000"); // EG2 Op. Ang. var 1 cell dist + 0.0202
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe30220000","32l51070a","01031v3x000000a0","0453503000000000"); // EG2 Op. Ang. var 1 cell dist + 0

    //-----
    //EG2: Omega Meson Cut Variations
    //-----
    //Std: 0453503000000000
  } else if(trainConfig == 3401)  { //EMCal + DCal EG2, Omega cut var. Background Scheme, Std 4 -> off
    //                         0008e013   411790109fe30220000   32l51070a   01031v3x00000000   0453503000000000
    //                                                                                          |
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe30220000","32l51070a","01031v3x00000000","0153503000000000"); // EG2, Om Event Mixing
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe30220000","32l51070a","01031v3x00000000","0a53503000000000"); // EG2, Om LikeSignMixing
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe30220000","32l51070a","01031v3x00000000","0d53503000000000"); // EG2, Om SideBandMixing
  } else if(trainConfig == 3402)  { //EMCal + DCal EG2, Omega cut var. rapidity, Std 5 -> -0.85, 0.85
    //                         0008e013   411790109fe30220000   32l51070a   01031v3x00000000   0453503000000000
    //                                                                                             |
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe30220000","32l51070a","01031v3x00000000","0453103000000000"); // EG2, Om rap. -0.8, 0.8
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe30220000","32l51070a","01031v3x00000000","0453603000000000"); // EG2, Om rap. -0.75, 0.75
  } else if(trainConfig == 3404)  { //EMCal + DCal EG2, Omega cut var. alpha, Std 3 -> 0.0-1.0
    //                         0008e013   411790109fe30220000   32l51070a   01031v3x00000000   0453503000000000
    //                                                                                               |
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe30220000","32l51070a","01031v3x00000000","0453505000000000"); // EG2 alpha 0-0.75
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe30220000","32l51070a","01031v3x00000000","0453506000000000"); // EG2 alpha 0-0.8
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe30220000","32l51070a","01031v3x00000000","0453507000000000"); // EG2 alpha 0-0.85
  } else if(trainConfig == 3410)  { //EMCal + DCal EG2, Omega cut var. Background Scheme single cfg, Std 4 -> off, EventMixing
    //                         0008e013   411790109fe30220000   32l51070a   01031v3x00000000   0453503000000000
    //                                                                                          |
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe30220000","32l51070a","01031v3x00000000","0153503000000000"); // EG2, Om Event Mixing
  } else if(trainConfig == 3411)  { //EMCal + DCal EG2, Omega cut var. Background Scheme single cfg, Std 4 -> off, LikeSignMixing
    //                         0008e013   411790109fe30220000   32l51070a   01031v3x00000000   0453503000000000
    //                                                                                          |
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe30220000","32l51070a","01031v3x00000000","0a53503000000000"); // EG2, Om LikeSignMixing
  } else if(trainConfig == 3412)  { //EMCal + DCal EG2, Omega cut var. Background Scheme single cfg, Std 4 -> off, SideBandMixing
    //                         0008e013   411790109fe30220000   32l51070a   01031v3x00000000   0453503000000000
    //                                                                                          |
    cuts.AddCutHeavyMesonCalo("0008e113","411790109fe30220000","32l51070a","01031v3x00000000","0d53503000000000"); // EG2, Om SideBandMixing
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //Standard Cuts of Pi0 Analysis: ("0008d113","411790109fe30220000","0r631031000000d0")
    //                                                                   ||  |        |
    //                                                                 "01031v3x00000000"
    //MesonCut r63==Background->ignored, d==OpeningAngle for Background->ignored, Pi0 MaxPt Cut v->ignored =>0453503000000000
  } else if(trainConfig == 4000)  { //EMCal + DCal EG1 Standard
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe30220000","32l51070a","01031v3x00000000","0453503000000000"); // EG1 Standard
    //-----
    //EG1: Event Variations
    //-----
    //Std: 0008d113
  } else if(trainConfig == 4001)  { //EMCal + DCal EG1, Event cut var. Remove Pileup, Std 1-> True
    //                         0008d113   411790109fe30220000   32l51070a   01031v3x00000000   0453503000000000
    //                              |
    cuts.AddCutHeavyMesonCalo("0008d013","411790109fe30220000","32l51070a","01031v3x00000000","0453503000000000"); // INT7 Pileup not removed
    //-----
    //EG1: Calo Variations
    //-----
    //Std: 411790109fe30220000
  } else if(trainConfig == 4101)  { //EMCal + DCal EG1, Calo cut var. NonLins, Std 01
    //                         0008d113   411790109fe30220000   32l51070a   01031v3x00000000   0453503000000000
    //                                         ||
    cuts.AddCutHeavyMesonCalo("0008d113","411799609fe30220000","32l51070a","01031v3x00000000","0453503000000000"); // EG1 no FT applied
    cuts.AddCutHeavyMesonCalo("0008d113","411799709fe30220000","32l51070a","01031v3x00000000","0453503000000000"); // EG1 EMC fine tuning applied
    cuts.AddCutHeavyMesonCalo("0008d113","411799809fe30220000","32l51070a","01031v3x00000000","0453503000000000"); // EG1 PCM-EMC fine tuning applied
  } else if(trainConfig == 4102)  { //EMCal + DCal EG1, Calo cut var. time, Std 9 -> -20+25
    //                         0008d113   411790109fe30220000   32l51070a   01031v3x00000000   0453503000000000
    //                                            |
    cuts.AddCutHeavyMesonCalo("0008d113","411790105fe30220000","32l51070a","01031v3x00000000","0453503000000000"); // EG1 time -50+50
    cuts.AddCutHeavyMesonCalo("0008d113","411790106fe30220000","32l51070a","01031v3x00000000","0453503000000000"); // EG1 time -30+35
    cuts.AddCutHeavyMesonCalo("0008d113","411790108fe30220000","32l51070a","01031v3x00000000","0453503000000000"); // EG1 time -20+30
    cuts.AddCutHeavyMesonCalo("0008d113","41179010afe30220000","32l51070a","01031v3x00000000","0453503000000000"); // EG1 time -12.5+13
  } else if(trainConfig == 4103)  { //EMCal + DCal EG1, Calo cut var. energy, Std 3 -> 0.7 GeV
    //                         0008d113   411790109fe30220000   32l51070a   01031v3x00000000   0453503000000000
    //                                          ||
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe20220000","32l51070a","01031v3x00000000","0453503000000000"); // EG1 energy 0.6 GeV
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe40220000","32l51070a","01031v3x00000000","0453503000000000"); // EG1 energy 0.8 GeV
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe50220000","32l51070a","01031v3x00000000","0453503000000000"); // EG1 energy 0.9 GeV // only meaningfull at higher pTs
  } else if(trainConfig == 4104)  { //EMCal + DCal EG1, Calo cut var. NCell, Std 0 -> Turned Off until 4GeV; then min 2 Cells
    //                         0008d113   411790109fe30220000   32l51070a   01031v3x00000000   0453503000000000
    //                                                |
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe32220000","32l51070a","01031v3x00000000","0453503000000000"); // EG1 NCells 2
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe3n220000","32l51070a","01031v3x00000000","0453503000000000"); // EG1 NCells 2 var (PCM-EMCal tagging corr)
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe3r220000","32l51070a","01031v3x00000000","0453503000000000"); // EG1 NCells 2 var (EMCal tagging corr)
  } else if(trainConfig == 4105)  { //EMCal + DCal EG1, Calo cut var. max M02, 2 -> INT7 M02 0.7
    //                         0008d113   411790109fe30220000   32l51070a   01031v3x00000000   0453503000000000
    //                                                  |
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe30230000","32l51070a","01031v3x00000000","0453503000000000"); // EG1 M02 0.5
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe30210000","32l51070a","01031v3x00000000","0453503000000000"); // EG1 M02 1.0
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe302k0000","32l51070a","01031v3x00000000","0453503000000000"); // EG1 M02 E dep
  } else if(trainConfig == 4106)  { //EMCal + DCal EG1, Calo cut var. TM, Std f,
    //  [1] + 1 / pow(x + pow(1 / ([0] - [1]), 1 / [2]), [2]; Eta (0.04, 0.010, 2.5); Phi (0.09, 0.015, 2.); EoverP 1.75
    //                         0008d113   411790109fe30220000   32l51070a   01031v3x00000000   0453503000000000
    //                                             |
    cuts.AddCutHeavyMesonCalo("0008d113","411790109ee30220000","32l51070a","01031v3x00000000","0453503000000000"); // EG1 TM var EoverP 2.00
    cuts.AddCutHeavyMesonCalo("0008d113","411790109ge30220000","32l51070a","01031v3x00000000","0453503000000000"); // EG1 TM var EoverP 1.50
    cuts.AddCutHeavyMesonCalo("0008d113","4117901097e30220000","32l51070a","01031v3x00000000","0453503000000000"); // EG1 TM var no EoverP
    cuts.AddCutHeavyMesonCalo("0008d113","411790109ne30220000","32l51070a","01031v3x00000000","0453503000000000"); // EG1 TM var, Eta (0.035, 0.010, 2.5); Phi (0.085, 0.015, 2.) EoverP 1.75
    cuts.AddCutHeavyMesonCalo("0008d113","411790109oe30220000","32l51070a","01031v3x00000000","0453503000000000"); // EG1 TM var, Eta (0.045, 0.010, 2.5); Phi (0.095, 0.015, 2.) EoverP 1.75
  } else if(trainConfig == 4107)  { //EMCal + DCal EG1, Calo cut var. Exotics, Std e, active F+ < 0.97
    //                         0008d113   411790109fe30220000   32l51070a   01031v3x00000000   0453503000000000
    //                                              |
    cuts.AddCutHeavyMesonCalo("0008d113","411790109f030220000","32l51070a","01031v3x00000000","0453503000000000"); // EG1 no exotics cut
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fb30220000","32l51070a","01031v3x00000000","0453503000000000"); // EG1 F+ < 0.95
    //-----
    //EG1: Primary Pion / Charged Pion (Pi+ Pi-) Variations
    //-----
    //Std: 32l51070a
  } else if(trainConfig == 4201)  { //EMCal + DCal EG1, Ch.Pi cut var. Ch.Pi ITS Requirement, Std 2 -> first or second SPD cluster required
    //                         0008d113   411790109fe30220000   32l51070a   01031v3x00000000   0453503000000000
    //                                                           |
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe30220000","34l51070a","01031v3x00000000","0453503000000000"); // EG1, Ch.Pi ITS, first or second SPD cluster required, min number of ITS clusters = 3
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe30220000","35l51070a","01031v3x00000000","0453503000000000"); // EG1, Ch.Pi ITS, ffirst or second SPD cluster required, min number of ITS clusters = 4
  } else if(trainConfig == 4202)  { //EMCal + DCal EG1, Ch.Pi cut var. Ch.Pi Cls TPC, Std l -> max shared clusters 10 + min TPC PID clusters 50
    //                         0008d113   411790109fe30220000   32l51070a   01031v3x00000000   0453503000000000
    //                                                            |
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe30220000","32e51070a","01031v3x00000000","0453503000000000"); // EG1, Ch.Pi, MinClsTPC 80. + Refit MaxSharedClsTPCFrac=0.
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe30220000","32k51070a","01031v3x00000000","0453503000000000"); // EG1, Ch.Pi, max shared clusters 10
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe30220000","32c51070a","01031v3x00000000","0453503000000000"); // EG1, Ch.Pi, MinClsTPC 80. + Refit
  } else if(trainConfig == 4203)  { //EMCal + DCal EG1, Ch.Pi cut var. Ch.Pi pT, Std 1 -> pt>0.1
    //                         0008d113   411790109fe30220000   32l51070a   01031v3x00000000   0453503000000000
    //                                                              |
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe30220000","32l50070a","01031v3x00000000","0453503000000000"); // EG1, Ch.Pi pt>0.075
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe30220000","32l52070a","01031v3x00000000","0453503000000000"); // EG1, Ch.Pi pt>0.125
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe30220000","32l53070a","01031v3x00000000","0453503000000000"); // EG1, Ch.Pi pt>0.15
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe30220000","32l54070a","01031v3x00000000","0453503000000000"); // EG1, Ch.Pi pt>0.4
  } else if(trainConfig == 4204)  { //EMCal + DCal EG1, Ch.Pi cut var. Ch.Pi TPC dEdx Low, Std 7 -> -3,3
    //                         0008d113   411790109fe30220000   32l51070a   01031v3x00000000   0453503000000000
    //                                                                |
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe30220000","32l510c0a","01031v3x00000000","0453503000000000"); // EG1, Ch.Pi -3.5,3.0
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe30220000","32l510d0a","01031v3x00000000","0453503000000000"); // EG1, Ch.Pi -3.25,3.0
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe30220000","32l510e0a","01031v3x00000000","0453503000000000"); // EG1, Ch.Pi -2.75,3.0
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe30220000","32l510f0a","01031v3x00000000","0453503000000000"); // EG1, Ch.Pi -2.5,3.0
  } else if(trainConfig == 4205)  { //EMCal + DCal EG1, Ch.Pi cut var. Ch.Pi TPC dEdx High, Std 7 -> -3,3
    //                         0008d113   411790109fe30220000   32l51070a   01031v3x00000000   0453503000000000
    //                                                                |
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe30220000","32l510g0a","01031v3x00000000","0453503000000000"); // EG1, Ch.Pi -3.0,2.5
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe30220000","32l510h0a","01031v3x00000000","0453503000000000"); // EG1, Ch.Pi -3.0,2.25
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe30220000","32l510i0a","01031v3x00000000","0453503000000000"); // EG1, Ch.Pi -3.0,3.25
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe30220000","32l510j0a","01031v3x00000000","0453503000000000"); // EG1, Ch.Pi -3.0,3.5
  } else if(trainConfig == 4206)  { //EMCal + DCal EG1, Ch.Pi cut var. Ch.Pi Mass, Std a -> Ch.Pi<850MeV
    //                         0008d113   411790109fe30220000   32l51070a   01031v3x00000000   0453503000000000
    //                                                                  |
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe30220000","32l51070r","01031v3x00000000","0453503000000000"); // EG1, Ch.Pi<800MeV
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe30220000","32l51070s","01031v3x00000000","0453503000000000"); // EG1, Ch.Pi<825MeV
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe30220000","32l51070t","01031v3x00000000","0453503000000000"); // EG1, Ch.Pi<875MeV
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe30220000","32l51070u","01031v3x00000000","0453503000000000"); // EG1, Ch.Pi<900MeV
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe30220000","32l51070n","01031v3x00000000","0453503000000000"); // EG1, Ch.Pi<1000MeV
    //-----
    //EG1: Neutral Meson (Pi0) Cut Variations
    //-----
    //Std: 01031v3x00000000
  } else if(trainConfig == 4302)  { //EMCal + DCal EG1, N.Pi cut var. rapidity, Std 1 -> -0.8, 0.8
    //                         0008d113   411790109fe30220000   32l51070a   01031v3x00000000   0453503000000000
    //                                                                          |
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe30220000","32l51070a","01035v3x00000000","0453503000000000"); // EG1, N.Pi rap. -0.85, 0.85
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe30220000","32l51070a","01036v3x00000000","0453503000000000"); // EG1, N.Pi rap. -0.75, 0.75
  } else if(trainConfig == 4303)  { //EMCal + DCal EG1, N.Pi cut var. maxMass, Std v -> 25GeV
    //                         0008d113   411790109fe30220000   32l51070a   01031v3x00000000   0453503000000000
    //                                                                           |
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe30220000","32l51070a","0103103x00000000","0453503000000000"); // EG1, N.Pi maxMass off
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe30220000","32l51070a","01031s3x00000000","0453503000000000"); // EG1, N.Pi maxMass 20GeV
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe30220000","32l51070a","01031y3x00000000","0453503000000000"); // EG1, N.Pi maxMass 22GeV
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe30220000","32l51070a","01031p3x00000000","0453503000000000"); // EG1, N.Pi maxMass 25GeV, 1 Gamma >11.GeV
  } else if(trainConfig == 4304)  { //EMCal + DCal EG1, N.Pi cut var. alpha, Std 3 -> 0.0-1.0
    //                         0008d113   411790109fe30220000   32l51070a   01031v3x00000000   0453503000000000
    //                                                                            |
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe30220000","32l51070a","01031v5x00000000","0453503000000000"); // EG1 alpha 0-0.75
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe30220000","32l51070a","01031v6x00000000","0453503000000000"); // EG1 alpha 0-0.8
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe30220000","32l51070a","01031v7x00000000","0453503000000000"); // EG1 alpha 0-0.85
  } else if(trainConfig == 4305)  { //EMCal + DCal EG1, N.Pi cut var. Selection Window, Std x -> 3 sigma
    //                         0008d113   411790109fe30220000   32l51070a   01031v3x00000000   0453503000000000
    //                                                                             |
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe30220000","32l51070a","01031v3100000000","0453503000000000"); // EG1, 2.0 sigma
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe30220000","32l51070a","01031v3w00000000","0453503000000000"); // EG1, 4.0 sigma
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe30220000","32l51070a","01031v3v00000000","0453503000000000"); // EG1, 2.5 sigma
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe30220000","32l51070a","01031v3r00000000","0453503000000000"); // EG1, 3.5 sigma
  } else if(trainConfig == 4306)  { //EMCal + DCal EG1, N.Pi cut var. open. angle, Std 0 -> off
    //                         0008d113   411790109fe30220000   32l51070a   01031v3x00000000   0453503000000000
    //                                                                                    |
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe30220000","32l51070a","01031v3x000000d0","0453503000000000"); // EG1 Op. Ang. var 1 cell dist + 0.017
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe30220000","32l51070a","01031v3x000000b0","0453503000000000"); // EG1 Op. Ang. var 1 cell dist + 0.0152
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe30220000","32l51070a","01031v3x000000g0","0453503000000000"); // EG1 Op. Ang. var 1 cell dist + 0.0202
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe30220000","32l51070a","01031v3x000000a0","0453503000000000"); // EG1 Op. Ang. var 1 cell dist + 0

    //-----
    //EG1: Omega Meson Cut Variations
    //-----
    //Std: 0453503000000000
  } else if(trainConfig == 4401)  { //EMCal + DCal EG1, Omega cut var. Background Scheme, Std 4 -> off
    //                         0008d113   411790109fe30220000   32l51070a   01031v3x00000000   0453503000000000
    //                                                                                          |
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe30220000","32l51070a","01031v3x00000000","0153503000000000"); // EG1, Om Event Mixing
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe30220000","32l51070a","01031v3x00000000","0a53503000000000"); // EG1, Om LikeSignMixing
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe30220000","32l51070a","01031v3x00000000","0d53503000000000"); // EG1, Om SideBandMixing
  } else if(trainConfig == 4402)  { //EMCal + DCal EG1, Omega cut var. rapidity, Std 5 -> -0.85, 0.85
    //                         0008d113   411790109fe30220000   32l51070a   01031v3x00000000   0453503000000000
    //                                                                                             |
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe30220000","32l51070a","01031v3x00000000","0453103000000000"); // EG1, Om rap. -0.8, 0.8
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe30220000","32l51070a","01031v3x00000000","0453603000000000"); // EG1, Om rap. -0.75, 0.75
  } else if(trainConfig == 4404)  { //EMCal + DCal EG1, Omega cut var. alpha, Std 3 -> 0.0-1.0
    //                         0008d113   411790109fe30220000   32l51070a   01031v3x00000000   0453503000000000
    //                                                                                               |
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe30220000","32l51070a","01031v3x00000000","0453505000000000"); // EG1 alpha 0-0.75
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe30220000","32l51070a","01031v3x00000000","0453506000000000"); // EG1 alpha 0-0.8
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe30220000","32l51070a","01031v3x00000000","0453507000000000"); // EG1 alpha 0-0.85
  } else if(trainConfig == 4410)  { //EMCal + DCal EG1, Omega cut var. Background Scheme single cfg, Std 4 -> off, EventMixing
    //                         0008d113   411790109fe30220000   32l51070a   01031v3x00000000   0453503000000000
    //                                                                                          |
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe30220000","32l51070a","01031v3x00000000","0153503000000000"); // EG1, Om Event Mixing
  } else if(trainConfig == 4411)  { //EMCal + DCal EG1, Omega cut var. Background Scheme single cfg, Std 4 -> off, LikeSignMixing
    //                         0008d113   411790109fe30220000   32l51070a   01031v3x00000000   0453503000000000
    //                                                                                          |
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe30220000","32l51070a","01031v3x00000000","0a53503000000000"); // EG1, Om LikeSignMixing
  } else if(trainConfig == 4412)  { //EMCal + DCal EG1, Omega cut var. Background Scheme single cfg, Std 4 -> off, SideBandMixing
    //                         0008d113   411790109fe30220000   32l51070a   01031v3x00000000   0453503000000000
    //                                                                                          |
    cuts.AddCutHeavyMesonCalo("0008d113","411790109fe30220000","32l51070a","01031v3x00000000","0d53503000000000"); // EG1, Om SideBandMixing
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // EMC pp 13 TeV Fitting, Systematics
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //Standard Cuts of Pi0 Analysis: ("00010113","24466190sa01cc00000","0163103100000010");
    //                                                                    |    |      |
    //                                                                 "0103103s00000000"
    //MesonCut 163==Background->ignored, s==Pi0 Mass Window->ignored, 1==OpeningAngle for Background->ignored =>0103103s00000000
  } else if(trainConfig == 6000)  { //PHOS INT7 Standard
    cuts.AddCutHeavyMesonCalo("00010113","24466190sa01cc00000","32l51070a","0103103s00000000","0453503000000000"); // INT7 Standard
    //-----
    //INT7: Event Variations
    //-----
    //Std: 00010113
  } else if(trainConfig == 6001)  { //PHOS INT7, Event cut var. Remove Pileup, Std 1-> True
    //                         00010113
    //                              |
    cuts.AddCutHeavyMesonCalo("00010013","24466190sa01cc00000","32l51070a","0103103s00000000","0453503000000000"); // INT7 Pileup not removed
    //-----
    //INT7: Calo Variations
    //-----
    //Std: 24466190sa01cc00000
  } else if(trainConfig == 6101)  { //PHOS INT7, Calo cut var. NonLins, Std 19
    //                                    24466190sa01cc00000
    //                                         ||
    cuts.AddCutHeavyMesonCalo("00010113","24466000sa01cc00000","32l51070a","0103103s00000000","0453503000000000"); // INT7 NoNL
    cuts.AddCutHeavyMesonCalo("00010113","24466110sa01cc00000","32l51070a","0103103s00000000","0453503000000000"); // INT7 NL11
    cuts.AddCutHeavyMesonCalo("00010113","24466120sa01cc00000","32l51070a","0103103s00000000","0453503000000000"); // INT7 NL12
  } else if(trainConfig == 6102)  { //PHOS INT7, Calo cut var. time, Std s -> -30+30 with timing efficiency
    //                                    24466190sa01cc00000
    //                                            |
    cuts.AddCutHeavyMesonCalo("00010113","24466190ra01cc00000","32l51070a","0103103s00000000","0453503000000000"); // INT7 r: LowPt from MB, 30ns
    cuts.AddCutHeavyMesonCalo("00010113","24466190ta01cc00000","32l51070a","0103103s00000000","0453503000000000"); // INT7 t: 25ns
    cuts.AddCutHeavyMesonCalo("00010113","24466190ua01cc00000","32l51070a","0103103s00000000","0453503000000000"); // INT7 u: 50ns
  } else if(trainConfig == 6103)  { //PHOS INT7, Calo cut var. energy, Std 1 -> 0.3 GeV
    //                                    24466190sa01cc00000
    //                                               |
    cuts.AddCutHeavyMesonCalo("00010113","24466190sa09cc00000","32l51070a","0103103s00000000","0453503000000000"); // INT7 9: energy 0.1 GeV
    cuts.AddCutHeavyMesonCalo("00010113","24466190sa02cc00000","32l51070a","0103103s00000000","0453503000000000"); // INT7 2: energy 0.5 GeV
    cuts.AddCutHeavyMesonCalo("00010113","24466190sa07cc00000","32l51070a","0103103s00000000","0453503000000000"); // INT7 7: energy 0.2 GeV
    cuts.AddCutHeavyMesonCalo("00010113","24466190sa08cc00000","32l51070a","0103103s00000000","0453503000000000"); // INT7 8: energy 0.4 GeV
  } else if(trainConfig == 6104)  { //PHOS INT7, Calo cut var. NCell & M02, Std cc -> min nCells = 2 >1GeV; M02 max=100, min=0.1, part 1
    //                                    24466190sa01cc00000
    //                                                |||
    cuts.AddCutHeavyMesonCalo("00010113","24466190sa012200000","32l51070a","0103103s00000000","0453503000000000"); // INT7 220: min nCells = 2, all E
    cuts.AddCutHeavyMesonCalo("00010113","24466190sa013200000","32l51070a","0103103s00000000","0453503000000000"); // INT7 320: min nCells = 3, all E
    cuts.AddCutHeavyMesonCalo("00010113","24466190sa01dc00000","32l51070a","0103103s00000000","0453503000000000"); // INT7 dc0: min nCells = 3, E>1GeV; minM02==0.1 off for E<1GeV
    cuts.AddCutHeavyMesonCalo("00010113","24466190sa01cd00000","32l51070a","0103103s00000000","0453503000000000"); // INT7 cd0: min nCells = 2, E>1GeV; minM02==0.2 off for E<1GeV
  } else if(trainConfig == 6105)  { //PHOS INT7, Calo cut var. NCell & M02, Std cc -> min nCells = 2 >1GeV; M02 max=100, min=0.1, part 2
    //                                    24466190sa01cc00000
    //                                                |||
    cuts.AddCutHeavyMesonCalo("00010113","24466190sa011000000","32l51070a","0103103s00000000","0453503000000000");// INT7 100: min nCells = 1, minM02 off
    cuts.AddCutHeavyMesonCalo("00010113","24466190sa01cc70000","32l51070a","0103103s00000000","0453503000000000");// INT7 cc7: maxM02 == 1.3
    cuts.AddCutHeavyMesonCalo("00010113","24466190sa01cc80000","32l51070a","0103103s00000000","0453503000000000");// INT7 cc8: maxM02 == 2.5
  } else if(trainConfig == 6106)  { //PHOS INT7, Calo cut var. TM, Std a
    //                                    24466190sa01cc00000
    //                                             |
    cuts.AddCutHeavyMesonCalo("00010113","24466190s001cc00000","32l51070a","0103103s00000000","0453503000000000"); // 0
    cuts.AddCutHeavyMesonCalo("00010113","24466190s101cc00000","32l51070a","0103103s00000000","0453503000000000"); // 1
    cuts.AddCutHeavyMesonCalo("00010113","24466190s401cc00000","32l51070a","0103103s00000000","0453503000000000"); // 4
    cuts.AddCutHeavyMesonCalo("00010113","24466190s501cc00000","32l51070a","0103103s00000000","0453503000000000"); // 5
  } else if(trainConfig == 6108)  { //PHOS INT7, Calo cut var. reconstructed conversion
    //                                    24466190sa01cc00000
    //                                                    |
    cuts.AddCutHeavyMesonCalo("00010113","24466190sa01cc00100","32l51070a","0103103s00000000","0453503000000000"); // 0: 0.02
    cuts.AddCutHeavyMesonCalo("00010113","24466190sa01cc00200","32l51070a","0103103s00000000","0453503000000000"); // 1: 0.025
    cuts.AddCutHeavyMesonCalo("00010113","24466190sa01cc00300","32l51070a","0103103s00000000","0453503000000000"); // 3: 0.03
    cuts.AddCutHeavyMesonCalo("00010113","24466190sa01cc00400","32l51070a","0103103s00000000","0453503000000000"); // 4: 0.035
    //-----
    //INT7: Primary Pion / Charged Pion (Pi+ Pi-) Variations
    //-----
    //Std: 32l51070a
  } else if(trainConfig == 6201)  { //PHOS INT7, Ch.Pi cut var. Ch.Pi ITS Requirement, Std 2 -> first or second SPD cluster required
    //                                                          32l51070a
    //                                                           |
    cuts.AddCutHeavyMesonCalo("00010113","24466190sa01cc00000","34l51070a","0103103s00000000","0453503000000000"); // INT7, Ch.Pi ITS, first or second SPD cluster required, min number of ITS clusters = 3
    cuts.AddCutHeavyMesonCalo("00010113","24466190sa01cc00000","35l51070a","0103103s00000000","0453503000000000"); // INT7, Ch.Pi ITS, first or second SPD cluster required, min number of ITS clusters = 4
  } else if(trainConfig == 6202)  { //PHOS INT7, Ch.Pi cut var. Ch.Pi Cls TPC, Std l -> max shared clusters 10 + min TPC PID clusters 50
    //                                                          32l51070a
    //                                                            |
    cuts.AddCutHeavyMesonCalo("00010113","24466190sa01cc00000","32e51070a","0103103s00000000","0453503000000000"); // INT7, Ch.Pi, MinClsTPC 80. + Refit MaxSharedClsTPCFrac=0.
    cuts.AddCutHeavyMesonCalo("00010113","24466190sa01cc00000","32k51070a","0103103s00000000","0453503000000000"); // INT7, Ch.Pi, max shared clusters 10
    cuts.AddCutHeavyMesonCalo("00010113","24466190sa01cc00000","32c51070a","0103103s00000000","0453503000000000"); // INT7, Ch.Pi, c -> MinClsTPC 80. + Refit
  } else if(trainConfig == 6203)  { //PHOS INT7, Ch.Pi cut var. Ch.Pi pT, Std 1 -> pt>0.1
    //                                                          32l51070a
    //                                                              |
    cuts.AddCutHeavyMesonCalo("00010113","24466190sa01cc00000","32l50070a","0103103s00000000","0453503000000000"); // INT7, Ch.Pi pt>0.075
    cuts.AddCutHeavyMesonCalo("00010113","24466190sa01cc00000","32l52070a","0103103s00000000","0453503000000000"); // INT7, Ch.Pi pt>0.125
    cuts.AddCutHeavyMesonCalo("00010113","24466190sa01cc00000","32l53070a","0103103s00000000","0453503000000000"); // INT7, Ch.Pi pt>0.15
    cuts.AddCutHeavyMesonCalo("00010113","24466190sa01cc00000","32l54070a","0103103s00000000","0453503000000000"); // INT7, Ch.Pi pt>0.4
  } else if(trainConfig == 6204)  { //PHOS INT7, Ch.Pi cut var. Ch.Pi TPC dEdx Low, Std 7 -> -3,3
    //                                                          32l51070a
    //                                                                |
    cuts.AddCutHeavyMesonCalo("00010113","24466190sa01cc00000","32l510c0a","0103103s00000000","0453503000000000"); // INT7, Ch.Pi -3.5,3.0
    cuts.AddCutHeavyMesonCalo("00010113","24466190sa01cc00000","32l510d0a","0103103s00000000","0453503000000000"); // INT7, Ch.Pi -3.25,3.0
    cuts.AddCutHeavyMesonCalo("00010113","24466190sa01cc00000","32l510e0a","0103103s00000000","0453503000000000"); // INT7, Ch.Pi -2.75,3.0
    cuts.AddCutHeavyMesonCalo("00010113","24466190sa01cc00000","32l510f0a","0103103s00000000","0453503000000000"); // INT7, Ch.Pi -2.5,3.0
  } else if(trainConfig == 6205)  { //PHOS INT7, Ch.Pi cut var. Ch.Pi TPC dEdx High, Std 7 -> -3,3
    //                                                          32l51070a
    //                                                                |
    cuts.AddCutHeavyMesonCalo("00010113","24466190sa01cc00000","32l510g0a","0103103s00000000","0453503000000000"); // INT7, Ch.Pi -3.0,2.5
    cuts.AddCutHeavyMesonCalo("00010113","24466190sa01cc00000","32l510h0a","0103103s00000000","0453503000000000"); // INT7, Ch.Pi -3.0,2.75
    cuts.AddCutHeavyMesonCalo("00010113","24466190sa01cc00000","32l510i0a","0103103s00000000","0453503000000000"); // INT7, Ch.Pi -3.0,3.25
    cuts.AddCutHeavyMesonCalo("00010113","24466190sa01cc00000","32l510j0a","0103103s00000000","0453503000000000"); // INT7, Ch.Pi -3.0,3.5
  } else if(trainConfig == 6206)  { //PHOS INT7, Ch.Pi cut var. Ch.Pi Mass, Std a -> Ch.Pi<850MeV
    //                                                          32l51070a
    //                                                                  |
    cuts.AddCutHeavyMesonCalo("00010113","24466190sa01cc00000","32l51070r","0103103s00000000","0453503000000000"); // INT7, Ch.Pi<800MeV
    cuts.AddCutHeavyMesonCalo("00010113","24466190sa01cc00000","32l51070s","0103103s00000000","0453503000000000"); // INT7, Ch.Pi<825MeV
    cuts.AddCutHeavyMesonCalo("00010113","24466190sa01cc00000","32l51070t","0103103s00000000","0453503000000000"); // INT7, Ch.Pi<875MeV
    cuts.AddCutHeavyMesonCalo("00010113","24466190sa01cc00000","32l51070u","0103103s00000000","0453503000000000"); // INT7, Ch.Pi<900MeV
    cuts.AddCutHeavyMesonCalo("00010113","24466190sa01cc00000","32l51070n","0103103s00000000","0453503000000000"); // INT7, Ch.Pi<1000MeV
    //-----
    //INT7: Neutral Meson (Pi0) Cut Variations
    //-----
    //Std: 0103103s00000000
  } else if(trainConfig == 6302)  { //PHOS INT7, N.Pi cut var. rapidity, Std 1 -> -0.8, 0.8
    //                                                                      0103103s00000000
    //                                                                          |
    cuts.AddCutHeavyMesonCalo("00010113","24466190sa01cc00000","32l51070a","0103503s00000000","0453503000000000"); // INT7, N.Pi rap. 5: -0.85, 0.85
    cuts.AddCutHeavyMesonCalo("00010113","24466190sa01cc00000","32l51070a","0103603s00000000","0453503000000000"); // INT7, N.Pi rap. 6: -0.75, 0.75
    cuts.AddCutHeavyMesonCalo("00010113","24466190sa01cc00000","32l51070a","0103403s00000000","0453503000000000"); // INT7, N.Pi rap. 4: -0.5, 0.5
    cuts.AddCutHeavyMesonCalo("00010113","24466190sa01cc00000","32l51070a","0103803s00000000","0453503000000000"); // INT7, N.Pi rap. 8: -0.25, 0.25
  } else if(trainConfig == 6304)  { //PHOS INT7, N.Pi cut var. alpha, Std 3 -> 0.0-1.0
    //                                                                      0103103s00000000
    //                                                                            |
    cuts.AddCutHeavyMesonCalo("00010113","24466190sa01cc00000","32l51070a","0103105s00000000","0453503000000000"); // INT7 alpha 0-0.75
    cuts.AddCutHeavyMesonCalo("00010113","24466190sa01cc00000","32l51070a","0103106s00000000","0453503000000000"); // INT7 alpha 0-0.8
    cuts.AddCutHeavyMesonCalo("00010113","24466190sa01cc00000","32l51070a","0103107s00000000","0453503000000000"); // INT7 alpha 0-0.85
  } else if(trainConfig == 6305)  { //PHOS INT7, N.Pi cut var. Selection Window, Std s -> 3 sigma
    //                                                                      0103103s00000000
    //                                                                             |
    cuts.AddCutHeavyMesonCalo("00010113","24466190sa01cc00000","32l51070a","0103103300000000","0453503000000000"); // INT7, 2.0 sigma
    cuts.AddCutHeavyMesonCalo("00010113","24466190sa01cc00000","32l51070a","0103103q00000000","0453503000000000"); // INT7, 4.0 sigma
    cuts.AddCutHeavyMesonCalo("00010113","24466190sa01cc00000","32l51070a","0103103l00000000","0453503000000000"); // INT7, 2.5 sigma
    cuts.AddCutHeavyMesonCalo("00010113","24466190sa01cc00000","32l51070a","0103103k00000000","0453503000000000"); // INT7, 3.5 sigma
  } else if(trainConfig == 6306)  { //PHOS INT7, N.Pi cut var. open. angle, Std 0 -> off
    //                                                                      0103103s00000000
    //                                                                                    |
    cuts.AddCutHeavyMesonCalo("00010113","24466190sa01cc00000","32l51070a","0103103s00000010","0453503000000000"); // INT7 Op. Ang. var 1: min opening angle 0.005
    cuts.AddCutHeavyMesonCalo("00010113","24466190sa01cc00000","32l51070a","0103103s00000030","0453503000000000"); // INT7 Op. Ang. var 3: min opening angle 0.01 -> 2 cell

    //-----
    //INT7: Omega Meson Cut Variations
    //-----
    //Std: 0453503000000000
  } else if(trainConfig == 6401)  { //PHOS INT7, Omega cut var. Background Scheme, Std 4 -> off
    //                                                                                         0453503000000000
    //                                                                                          |
    cuts.AddCutHeavyMesonCalo("00010113","24466190sa01cc00000","32l51070a","0103103s00000000","0153503000000000"); // INT7, Om Event Mixing
    cuts.AddCutHeavyMesonCalo("00010113","24466190sa01cc00000","32l51070a","0103103s00000000","0a53503000000000"); // INT7, Om LikeSignMixing
    cuts.AddCutHeavyMesonCalo("00010113","24466190sa01cc00000","32l51070a","0103103s00000000","0d53503000000000"); // INT7, Om SideBandMixing
  } else if(trainConfig == 6402)  { //PHOS INT7, Omega cut var. rapidity, Std 5 -> -0.85, 0.85
    //                                                                                         0453503000000000
    //                                                                                             |
    cuts.AddCutHeavyMesonCalo("00010113","24466190sa01cc00000","32l51070a","0103103s00000000","0453103000000000"); // INT7, Om rap. -0.8, 0.8
    cuts.AddCutHeavyMesonCalo("00010113","24466190sa01cc00000","32l51070a","0103103s00000000","0453603000000000"); // INT7, Om rap. -0.75, 0.75
  } else if(trainConfig == 6404)  { //PHOS INT7, Omega cut var. alpha, Std 3 -> 0.0-1.0
    //                                                                                         0453503000000000
    //                                                                                               |
    cuts.AddCutHeavyMesonCalo("00010113","24466190sa01cc00000","32l51070a","0103103s00000000","0453505000000000"); // INT7 alpha 0-0.75
    cuts.AddCutHeavyMesonCalo("00010113","24466190sa01cc00000","32l51070a","0103103s00000000","0453506000000000"); // INT7 alpha 0-0.8
    cuts.AddCutHeavyMesonCalo("00010113","24466190sa01cc00000","32l51070a","0103103s00000000","0453507000000000"); // INT7 alpha 0-0.85
  } else if(trainConfig == 6410)  { //PHOS INT7, Omega cut var. Background Scheme single cfg, Std 4 -> off, EventMixing
    //                                                                                         0453503000000000
    //                                                                                          |
    cuts.AddCutHeavyMesonCalo("00010113","24466190sa01cc00000","32l51070a","0103103s00000000","0153503000000000"); // INT7, Om Event Mixing
  } else if(trainConfig == 6411)  { //PHOS INT7, Omega cut var. Background Scheme single cfg, Std 4 -> off, LikeSignMixing
    //                                                                                         0453503000000000
    //                                                                                          |
    cuts.AddCutHeavyMesonCalo("00010113","24466190sa01cc00000","32l51070a","0103103s00000000","0a53503000000000"); // INT7, Om LikeSignMixing
  } else if(trainConfig == 6412)  { //PHOS INT7, Omega cut var. Background Scheme single cfg, Std 4 -> off, SideBandMixing
    //                                                                                         0453503000000000
    //                                                                                          |
    cuts.AddCutHeavyMesonCalo("00010113","24466190sa01cc00000","32l51070a","0103103s00000000","0d53503000000000"); // INT7, Om SideBandMixing

    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



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
  if( numberOfCuts > 1 && enableMLBckRedStudy) {
    cout << "\n\n****************************************************" << endl;
    cout << "ERROR: Trees for ML studies implemented only for one cut at the time. Returning..." << endl;
    cout << "****************************************************\n\n" << endl;
    return;
  }

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
    if(corrTaskSetting.CompareTo("")){
      TrackMatcherName = TrackMatcherName+"_"+corrTaskSetting.Data();
      cout << "Using separate track matcher for correction framework setting: " << TrackMatcherName.Data() << endl;
    }
    if( !(AliCaloTrackMatcher*)mgr->GetTask(TrackMatcherName.Data()) ){
      AliCaloTrackMatcher* fTrackMatcher = new AliCaloTrackMatcher(TrackMatcherName.Data(),caloCutPos.Atoi(),trackMatcherRunningMode);
      fTrackMatcher->SetV0ReaderName(V0ReaderName);
      mgr->AddTask(fTrackMatcher);
      mgr->ConnectInput(fTrackMatcher,0,cinput);
    }

    analysisEventCuts[i] = new AliConvEventCuts();
    analysisEventCuts[i]->SetV0ReaderName(V0ReaderName);
    analysisEventCuts[i]->SetCorrectionTaskSetting(corrTaskSetting);
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
    if(runLightOutput>=1) analysisEventCuts[i]->SetLightOutput(kTRUE);
    analysisEventCuts[i]->InitializeCutsFromCutString((cuts.GetEventCut(i)).Data());
    if (periodNameV0Reader.CompareTo("") != 0) analysisEventCuts[i]->SetPeriodEnum(periodNameV0Reader);
    EventCutList->Add(analysisEventCuts[i]);
    analysisEventCuts[i]->SetFillCutHistograms("",kFALSE);


    analysisClusterCuts[i] = new AliCaloPhotonCuts();
    analysisClusterCuts[i]->SetV0ReaderName(V0ReaderName);
    analysisClusterCuts[i]->SetCorrectionTaskSetting(corrTaskSetting);

    if(runLightOutput>=4) {
        analysisClusterCuts[i]->SetLightOutput(2);
    } else if(runLightOutput>=1) {
        analysisClusterCuts[i]->SetLightOutput(1);
    }

    analysisClusterCuts[i]->SetExtendedMatchAndQA(enableExtMatchAndQA);
    analysisClusterCuts[i]->SetCaloTrackMatcherName(TrackMatcherName);
    if( ! analysisClusterCuts[i]->InitializeCutsFromCutString((cuts.GetClusterCut(i)).Data()) ) {
      std::cout<<"ERROR: analysisClusterCuts [" <<i<<"]"<<std::endl;
      return ;
    } else {
      analysisClusterCuts[i]->InitializeCutsFromCutString((cuts.GetClusterCut(i)).Data());
      ClusterCutList->Add(analysisClusterCuts[i]);
      analysisClusterCuts[i]->SetFillCutHistograms("");
    }

    analysisNeutralPionCuts[i] = new AliConversionMesonCuts();
    analysisNeutralPionCuts[i]->SetUsePtDepSelectionWindow(usePtDepSelectionWindowCut);
    if(runLightOutput>=4) {
        analysisNeutralPionCuts[i]->SetLightOutput(2);
    } else if(runLightOutput>=1) {
        analysisNeutralPionCuts[i]->SetLightOutput(1);
    }
    if( ! analysisNeutralPionCuts[i]->InitializeCutsFromCutString((cuts.GetNDMCut(i)).Data()) ) {
      std::cout<<"ERROR: analysisMesonCuts [ " <<i<<" ] "<<std::endl;
      return ;
    } else {
      NeutralPionCutList->Add(analysisNeutralPionCuts[i]);
      analysisNeutralPionCuts[i]->SetFillCutHistograms("");
    }

    analysisMesonCuts[i] = new AliConversionMesonCuts();
    if(runLightOutput>=4) {
        analysisMesonCuts[i]->SetLightOutput(2);
    } else if(runLightOutput>=1) {
        analysisMesonCuts[i]->SetLightOutput(1);
    }
    if( ! analysisMesonCuts[i]->InitializeCutsFromCutString((cuts.GetMesonCut(i)).Data()) ) {
      std::cout<<"ERROR: analysisMesonCuts [ " <<i<<" ] "<<std::endl;
      return ;
    } else {
      MesonCutList->Add(analysisMesonCuts[i]);
      analysisMesonCuts[i]->SetFillCutHistograms("");
    }
    analysisEventCuts[i]->SetAcceptedHeader(HeaderList);

    TString cutName( Form("%s_%s_%s_%s_%s",(cuts.GetEventCut(i)).Data(), (cuts.GetClusterCut(i)).Data(),(cuts.GetPionCut(i)).Data(),(cuts.GetNDMCut(i)).Data(), (cuts.GetMesonCut(i)).Data() ) );
    analysisPionCuts[i] = new AliPrimaryPionCuts();
    analysisPionCuts[i]->SetPrefilterRunFlag(prefilterRunFlag);
    analysisPionCuts[i]->SetPeriodName(periodNameV0Reader);
    if(runLightOutput>=1) analysisPionCuts[i]->SetLightOutput(kTRUE);

        if( !analysisPionCuts[i]->InitializeCutsFromCutString((cuts.GetPionCut(i)).Data())) {
      std::cout<< "ERROR:  analysisPionCuts [ " <<i<<" ] "<<std::endl;
      return ;
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
  task->SetCorrectionTaskSetting(corrTaskSetting);

  task->SetMoveParticleAccordingToVertex(kFALSE);
  task->SetSelectedHeavyNeutralMeson(selectHeavyNeutralMeson);

  task->SetDoMesonQA(enableQAMesonTask );

  task->SetUnsmearedOutputs(unsmearingoutputs);

  task->SetEnableSortingOfMCClusLabels(enableSortingMCLabels);

  //connect containers
  AliAnalysisDataContainer *couttree  = 0;
  AliAnalysisDataContainer *coutput =
  mgr->CreateContainer(!(corrTaskSetting.CompareTo("")) ? Form("GammaConvNeutralMesonPiPlPiMiNeutralMeson_%i_%i_%i.root",selectHeavyNeutralMeson,neutralPionMode, trainConfig) : Form("GammaConvNeutralMesonPiPlPiMiNeutralMeson_%i_%i_%i_%s.root",selectHeavyNeutralMeson,neutralPionMode, trainConfig, corrTaskSetting.Data()), TList::Class(),
              AliAnalysisManager::kOutputContainer,Form("GammaConvNeutralMesonPiPlPiMiNeutralMeson_%i_%i_%i.root",selectHeavyNeutralMeson,neutralPionMode, trainConfig));

  if(enableMLBckRedStudy){
    couttree = mgr->CreateContainer( Form("GammaConvNeutralMesonPiPlPiMiNeutralMesonTree_%i_%i_%i.root",selectHeavyNeutralMeson,neutralPionMode, trainConfig),
                                    TTree::Class(),
                                    AliAnalysisManager::kOutputContainer,
                                    Form("GammaConvNeutralMesonPiPlPiMiNeutralMesonTree_%i_%i_%i.root",selectHeavyNeutralMeson,neutralPionMode, trainConfig) );
  }
  
  mgr->AddTask(task);
  mgr->ConnectInput(task,0,cinput);
  mgr->ConnectOutput(task,1,coutput);
  if(enableMLBckRedStudy) mgr->ConnectOutput(task,2,couttree);

  return;

}
