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
void AddTask_GammaConvNeutralMesonPiPlPiMiNeutralMeson_ConvMode_pp(
    Int_t     trainConfig                   = 1,
    Int_t     isMC                          = 0,                        //run MC
    TString   photonCutNumberV0Reader       = "",                       // 00000008400000000100000000 nom. B, 00000088400000000100000000 low B
    Int_t     selectHeavyNeutralMeson       = 0,                        //run eta prime instead of omega
    Int_t     enableQAMesonTask             = 1,                        //enable QA in AliAnalysisTaskNeutralMesonToPiPlPiMiNeutralMeson
    Int_t     enableTriggerMimicking        = 0,                        // enable trigger mimicking
    Bool_t    enableTriggerOverlapRej       = kFALSE,                   // enable trigger overlap rejection    
    // settings for weights
    // FPTW:fileNamePtWeights, FMUW:fileNameMultWeights,  FMAW:fileNameMatBudWeights,  separate with ;
    // Material Budget Weights file for Run 2
    // FMAW:alien:///alice/cern.ch/user/a/amarin//MBW/MCInputFileMaterialBudgetWeightsLHC16_Pythia_00010103_0d000009266300008850404000_date181214.root
    TString   fileNameExternalInputs        = "MCSpectraInput.root",    //
    Bool_t    doWeighting                   = kFALSE,                   //enable Weighting
    Bool_t    enableElecDeDxPostCalibration = kFALSE,                   // enable post calibration of elec pos dEdX
    TString   generatorName                 = "HIJING",
    Double_t  tolerance                     = -1,
    TString   periodNameV0Reader            = "",                       // period Name for V0Reader
    Int_t     runLightOutput                = 0,                        // run light output option 0: no light output 1: most cut histos stiched off 2: unecessary omega hists turned off as well
    Int_t     prefilterRunFlag              = 1500,                     // flag to change the prefiltering of ESD tracks. See SetHybridTrackCutsAODFiltering() in AliPrimaryPionCuts
    Bool_t    usePtDepSelectionWindowCut    = kFALSE,                   // use pt dependent meson selection window cut
    Bool_t    usePreSelection               = kTRUE,
    Int_t     enableMatBudWeightsPi0        = 0,                        // 1 = three radial bins, 2 = 10 radial bins (2 is the default when using weights)
    TString   additionalTrainConfig         = "0"                       // additional counter for trainconfig, this has to be always the last parameter
  ) {
  AliCutHandlerPCM cuts(13);

  TString addTaskName                       = "AddTask_GammaConvNeutralMesonPiPlPiMiNeutralMeson_ConvMode_pp";
  TString fileNamePtWeights                 = cuts.GetSpecialFileNameFromString (fileNameExternalInputs, "FPTW:");
  TString fileNameMultWeights               = cuts.GetSpecialFileNameFromString (fileNameExternalInputs, "FMUW:");
  TString fileNameMatBudWeights             = cuts.GetSpecialFileNameFromString (fileNameExternalInputs, "FMAW:");
  TString fileNamedEdxPostCalib             = cuts.GetSpecialFileNameFromString (fileNameExternalInputs, "FEPC:");
  TString fileNameCustomTriggerMimicOADB    = cuts.GetSpecialFileNameFromString (fileNameExternalInputs, "FTRM:");

  if(additionalTrainConfig.Contains("MaterialBudgetWeights"))
    fileNameMatBudWeights         = cuts.GetSpecialSettingFromAddConfig(additionalTrainConfig, "MaterialBudgetWeights",fileNameMatBudWeights, addTaskName);
  //parse additionalTrainConfig flag
  TObjArray *rAddConfigArr = additionalTrainConfig.Tokenize("_");
  if(rAddConfigArr->GetEntries()<1){cout << "ERROR during parsing of additionalTrainConfig String '" << additionalTrainConfig.Data() << "'" << endl; return;}
  TObjString* rAdditionalTrainConfig;
  for(Int_t i = 0; i<rAddConfigArr->GetEntries() ; i++){
    if(i==0) rAdditionalTrainConfig = (TObjString*)rAddConfigArr->At(i);
    else{
      TObjString* temp = (TObjString*) rAddConfigArr->At(i);
      TString tempStr = temp->GetString();
      cout << "INFO: nothing to do, no definition available!" << endl;
    }
  }
  TString sAdditionalTrainConfig = rAdditionalTrainConfig->GetString();
  if (sAdditionalTrainConfig.Atoi() > 0){
    trainConfig = trainConfig + sAdditionalTrainConfig.Atoi();
    cout << "INFO: AddTask_GammaConvNeutralMesonPiPlPiMiNeutralMeson_ConvMode_pp running additionalTrainConfig '" << sAdditionalTrainConfig.Atoi() << "', train config: '" << trainConfig << "'" << endl;
  }

  Int_t isHeavyIon = 0;
  Int_t neutralPionMode = 0;

  // ================== GetAnalysisManager ===============================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error(Form("AddTask_GammaConvNeutralMesonPiPlPiMiNeutralMeson_ConvMode_pp_%i",trainConfig), "No analysis manager found.");
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
  TString PionCuts          = "000000200";            //Electron Cuts
  if (!usePreSelection){                              //no PreSelection dEdx applied
      PionCuts= "000000000";
      cout<<"Preselection Disabled, PionCuts set to "<<PionCuts<<endl;
  } else {                                            //PreSelection dEdx applied
      cout<<"Preselection Enabled, PionCuts set to "<<PionCuts<<endl;
  }

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
    //connect input V0Reader
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

  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //                                          ETA MESON
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  if( trainConfig == 1 ) {
    // everything open, min pt charged pi = 100 MeV
    cuts.AddCutHeavyMesonPCM("00000113","00200009327000008250400000","000010400","0103503a00000000","0103503000000000");
  } else if( trainConfig == 2 ) {
    // closing charged pion cuts, minimum TPC cluster = 80, TPC dEdx pi = \pm 3 sigma, min pt charged pi = 100 MeV
    cuts.AddCutHeavyMesonPCM("00000113","00200009327000008250400000","002010700","0103503a00000000","0103503000000000");
  } else if( trainConfig == 3) {
    // eta < 0.9
    // closing charged pion cuts, minimum TPC cluster = 80, TPC dEdx pi = \pm 3 sigma, pi+pi- mass cut of 0.65, min pt charged pi = 100 MeV
    // closing neural pion cuts, 0.1 < M_gamma,gamma < 0.15
    // maxChi2 per cluster TPC <4, require TPC refit, DCA XY pT dependend 0.0182+0.0350/pt^1.01, DCA_Z = 3.0
    // timing cluster cut open
    cuts.AddCutHeavyMesonPCM("00010113","00200009327000008250400000","30a330708","0103503400000000","0153503000000000"); // all of the above
    cuts.AddCutHeavyMesonPCM("00052113","00200009327000008250400000","30a330708","0103503400000000","0153503000000000"); // all of the above
    cuts.AddCutHeavyMesonPCM("00062113","00200009327000008250400000","30a330708","0103503400000000","0153503000000000"); // all of the above
    cuts.AddCutHeavyMesonPCM("00083113","00200009327000008250400000","30a330708","0103503400000000","0153503000000000"); // all of the above
    cuts.AddCutHeavyMesonPCM("00085113","00200009327000008250400000","30a330708","0103503400000000","0153503000000000"); // all of the above
  } else if( trainConfig == 4) {
    // same as 3 but only MB
    cuts.AddCutHeavyMesonPCM("00010113","00200009327000008250400000","30a330708","0103503400000000","0153503000000000"); // all of the above

  } else if ( trainConfig == 10) { // Standard cut (for now)
    cuts.AddCutHeavyMesonPCM("00000113","00200009227000008250400000","302010708","0103603800000000","0153503000000000");
    cuts.AddCutHeavyMesonPCM("00000113","00200009227000008250400000","322010708","0103603800000000","0153503000000000"); // with ITS requirement
  } else if ( trainConfig == 11) { // test 5 TeV
    cuts.AddCutHeavyMesonPCM("00010113","00200009227000008250400000","302010708","0103603800000000","0153503000000000");
    cuts.AddCutHeavyMesonPCM("00010113","00200009227000008250400000","322010708","0103603800000000","0153503000000000"); // with ITS requirement

    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //                                          OMEGA MESON
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  } else if( trainConfig == 100 ) {
    // everything open, min pt charged pi = 100 MeV
    cuts.AddCutHeavyMesonPCM("00000113","00200009327000008250400000","000010400","0103503a00000000","0103503000000000");
  } else if( trainConfig == 101 ) {
    // closing charged pion cuts, minimum TPC cluster = 80, TPC dEdx pi = \pm 3 sigma, min pt charged pi = 100 MeV
    cuts.AddCutHeavyMesonPCM("00000113","00200009327000008250400000","002010700","0103503a00000000","0103503000000000");
  } else if( trainConfig == 102) {
    // eta < 0.9
    // closing charged pion cuts, minimum TPC cluster = 80, TPC dEdx pi = \pm 3 sigma, pi+pi- mass cut of 0.65, min pt charged pi = 100 MeV
    // closing neural pion cuts, 0.1 < M_gamma,gamma < 0.15
    // maxChi2 per cluster TPC <4, require TPC refit, DCA XY pT dependend 0.0182+0.0350/pt^1.01, DCA_Z = 3.0
    // timing cluster cut open
    cuts.AddCutHeavyMesonPCM("00010113","00200009327000008250400000","30a330708","0103503400000000","0153503000000000"); // all of the above
    cuts.AddCutHeavyMesonPCM("00052113","00200009327000008250400000","30a330708","0103503400000000","0153503000000000"); // all of the above
    cuts.AddCutHeavyMesonPCM("00062113","00200009327000008250400000","30a330708","0103503400000000","0153503000000000"); // all of the above
    cuts.AddCutHeavyMesonPCM("00083113","00200009327000008250400000","30a330708","0103503400000000","0153503000000000"); // all of the above
    cuts.AddCutHeavyMesonPCM("00085113","00200009327000008250400000","30a330708","0103503400000000","0153503000000000"); // all of the above

  } else if( trainConfig == 103) {
    // same as 102 but only MB
    cuts.AddCutHeavyMesonPCM("00010113","00200009327000008250400000","30a330708","0103503400000000","0153503000000000"); // all of the above

    // pp 7 TeV
  } else if ( trainConfig == 110) { // Standard cut (for now)
    cuts.AddCutHeavyMesonPCM("00000113","00200009227000008250400000","32c010708","0103603500000000","0153503000000000");
    // pp 5 TeV test
  } else if ( trainConfig == 111) { // with TPC refit + ITS requirement
    cuts.AddCutHeavyMesonPCM("00010113","00200009227000008250400000","32c010708","0103603500000000","0153503000000000"); // V0AND
    cuts.AddCutHeavyMesonPCM("00052113","00200009227000008250400000","32c010708","0103603500000000","0153503000000000"); // EMC7
    cuts.AddCutHeavyMesonPCM("00083113","00200009227000008250400000","32c010708","0103603500000000","0153503000000000"); // EG1
    cuts.AddCutHeavyMesonPCM("00085113","00200009227000008250400000","32c010708","0103603500000000","0153503000000000"); // EG2
    cuts.AddCutHeavyMesonPCM("00062113","00200009227000008250400000","32c010708","0103603500000000","0153503000000000"); // PHI7
    // pp 13 TeV test
  } else if ( trainConfig == 112) { // with TPC refit + ITS requirement
    cuts.AddCutHeavyMesonPCM("00010113","00200009227000008250400000","32c51070a","0103603500000000","0153503000000000"); // V0AND
    cuts.AddCutHeavyMesonPCM("00052113","00200009227000008250400000","32c51070a","0103603500000000","0153503000000000"); // EMC7
    cuts.AddCutHeavyMesonPCM("00083113","00200009227000008250400000","32c51070a","0103603500000000","0153503000000000"); // EG1
    cuts.AddCutHeavyMesonPCM("00085113","00200009227000008250400000","32c51070a","0103603500000000","0153503000000000"); // EG2
    cuts.AddCutHeavyMesonPCM("00062113","00200009227000008250400000","32c51070a","0103603500000000","0153503000000000"); // PHI7
 } else if ( trainConfig == 113) { // pp13 TeV AOD and ESD comparison
    cuts.AddCutHeavyMesonPCM("00010113","00200009227000008250400000","32c510708","0103603500000000","0153503000000000"); // V0AND
    // LHC11 7 TeV triggered test
  } else if ( trainConfig == 115) { // with TPC refit + ITS requirement
    cuts.AddCutHeavyMesonPCM("00010113","00200009227000008250400000","32c010708","0103603500000000","0153503000000000"); // V0AND
    cuts.AddCutHeavyMesonPCM("00052113","00200009227000008250400000","32c010708","0103603500000000","0153503000000000"); // EMC7
    cuts.AddCutHeavyMesonPCM("00062113","00200009227000008250400000","32c010708","0103603500000000","0153503000000000"); // PHI7

    // --------------------------
    // systematic studies 7 TeV
    // --------------------------

    // charged pion cuts
  } else if ( trainConfig == 120) {
    cuts.AddCutHeavyMesonPCM("00000113","0dm0000922700000dge0404000","32c51070a","0103603a00000000","0153503000000000"); // use 4 vec
    cuts.AddCutHeavyMesonPCM("00000113","0dm0000922700000dge0404000","32c51070a","0103603o00000000","0153503000000000"); // use 4 vec
  } else if ( trainConfig == 121) { // pTCut
    cuts.AddCutHeavyMesonPCM("00000113","0dm0000922700000dge0404000","32c50070a","0103603500000000","0153503000000000"); // pt>0.075
    cuts.AddCutHeavyMesonPCM("00000113","0dm0000922700000dge0404000","32c52070a","0103603500000000","0153503000000000"); // pt>0.125
    cuts.AddCutHeavyMesonPCM("00000113","0dm0000922700000dge0404000","32c53070a","0103603500000000","0153503000000000"); // pt>0.15
    cuts.AddCutHeavyMesonPCM("00000113","0dm0000922700000dge0404000","32c54070a","0103603500000000","0153503000000000"); // pt>0.4
  } else if ( trainConfig == 122) { // TPCdEdxCutPion
    cuts.AddCutHeavyMesonPCM("00000113","0dm0000922700000dge0404000","32c51050a","0103603500000000","0153503000000000"); // -4,5
    cuts.AddCutHeavyMesonPCM("00000113","0dm0000922700000dge0404000","32c51080a","0103603500000000","0153503000000000"); // -2.5,4
    cuts.AddCutHeavyMesonPCM("00000113","0dm0000922700000dge0404000","32c51020a","0103603500000000","0153503000000000"); // -6,7
    cuts.AddCutHeavyMesonPCM("00000113","0dm0000922700000dge0404000","32c51030a","0103603500000000","0153503000000000"); // -4,4
    cuts.AddCutHeavyMesonPCM("00000113","0dm0000922700000dge0404000","32c51040a","0103603500000000","0153503000000000"); // -2.5,4
    // Neutral Pion Cuts
  } else if ( trainConfig == 123) { // invariant mass cut (std 110-150)
    cuts.AddCutHeavyMesonPCM("00000113","0dm0000922700000dge0404000","32c51070a","0103603100000000","0153503000000000"); // 100-145
    cuts.AddCutHeavyMesonPCM("00000113","0dm0000922700000dge0404000","32c51070a","0103603200000000","0153503000000000"); // 110-145
    cuts.AddCutHeavyMesonPCM("00000113","0dm0000922700000dge0404000","32c51070a","0103603300000000","0153503000000000"); // 120-145
    cuts.AddCutHeavyMesonPCM("00000113","0dm0000922700000dge0404000","32c51070a","0103603800000000","0153503000000000"); // 125-145
    cuts.AddCutHeavyMesonPCM("00000113","0dm0000922700000dge0404000","32c51070a","0103603a00000000","0153503000000000"); //  80-145
  } else if ( trainConfig == 130) {
    cuts.AddCutHeavyMesonPCM("00000113","0dm0000922700000dge0404000","32b51070a","0103603500000000","0453503000000000"); // without background

  // Conversion Systematics
  } else if ( trainConfig == 140){ // pileup
    cuts.AddCutHeavyMesonPCM("00000113","00200009227000008250400000","32c51070a","0103603500000000","0153503000000000");
  } else if ( trainConfig == 141){ // singlept
    cuts.AddCutHeavyMesonPCM("00000113","00200019227000008250400000","32c51070a","0103603500000000","0153503000000000"); // 0.100 GeV
    cuts.AddCutHeavyMesonPCM("00000113","00200049227000008250400000","32c51070a","0103603500000000","0153503000000000"); // 0.075 GeV
    cuts.AddCutHeavyMesonPCM("00000113","00200069227000008250400000","32c51070a","0103603500000000","0153503000000000"); // 0.04 GeV
    cuts.AddCutHeavyMesonPCM("00000113","00200059227000008250400000","32c51070a","0103603500000000","0153503000000000"); // 0.125 GeV
  } else if ( trainConfig == 142){ // clstpc
    cuts.AddCutHeavyMesonPCM("00000113","00200008227000008250400000","32c51070a","0103603500000000","0153503000000000"); // fMinClsTPCToF= 0.35;fUseCorrectedTPCClsInfo=1;
    cuts.AddCutHeavyMesonPCM("00000113","00200006227000008250400000","32c51070a","0103603500000000","0153503000000000"); // fMinClsTPCToF= 0.70;fUseCorrectedTPCClsInfo=1;
    cuts.AddCutHeavyMesonPCM("00000113","00200001227000008250400000","32c51070a","0103603500000000","0153503000000000"); // fMinClsTPCToF= 60;fUseCorrectedTPCClsInfo=0;
    cuts.AddCutHeavyMesonPCM("00000113","00200002227000008250400000","32c51070a","0103603500000000","0153503000000000"); // fMinClsTPCToF= 80;fUseCorrectedTPCClsInfo=0;
    cuts.AddCutHeavyMesonPCM("00000113","00200003227000008250400000","32c51070a","0103603500000000","0153503000000000"); // fMinClsTPCToF= 100;fUseCorrectedTPCClsInfo=0;
  } else if ( trainConfig == 143){ // TPCdEdxCutElectron
    cuts.AddCutHeavyMesonPCM("00000113","00200009327000008250400000","32c51070a","0103603500000000","0153503000000000"); // -4,5
    cuts.AddCutHeavyMesonPCM("00000113","00200009627000008250400000","32c51070a","0103603500000000","0153503000000000"); // -2.5,4
    cuts.AddCutHeavyMesonPCM("00000113","00200009427000008250400000","32c51070a","0103603500000000","0153503000000000"); // -6,7
    cuts.AddCutHeavyMesonPCM("00000113","00200009527000008250400000","32c51070a","0103603500000000","0153503000000000"); // -4,4
    cuts.AddCutHeavyMesonPCM("00000113","00200009627000008250400000","32c51070a","0103603500000000","0153503000000000"); // -2.5,4
  } else if ( trainConfig == 144){ // TPCdEdxCutPion
    cuts.AddCutHeavyMesonPCM("00000113","00200009217000008250400000","32c51070a","0103603500000000","0153503000000000"); // fPIDnSigmaAbovePionLine=0; fPIDnSigmaAbovePionLineHighPt=-10;
    cuts.AddCutHeavyMesonPCM("00000113","00200009237000008250400000","32c51070a","0103603500000000","0153503000000000"); // fPIDnSigmaAbovePionLine=2.5; fPIDnSigmaAbovePionLineHighPt=-10;
    cuts.AddCutHeavyMesonPCM("00000113","00200009247000008250400000","32c51070a","0103603500000000","0153503000000000"); // fPIDnSigmaAbovePionLine=0.5; fPIDnSigmaAbovePionLineHighPt=-10;
    cuts.AddCutHeavyMesonPCM("00000113","00200009257000008250400000","32c51070a","0103603500000000","0153503000000000"); // fPIDnSigmaAbovePionLine=2; fPIDnSigmaAbovePionLineHighPt=-10;
  } else if ( trainConfig == 145){ // QtMaxCut
    cuts.AddCutHeavyMesonPCM("00000113","00200009227000003250400000","32c51070a","0103603500000000","0153503000000000");
    cuts.AddCutHeavyMesonPCM("00000113","00200009227000009250400000","32c51070a","0103603500000000","0153503000000000");
    cuts.AddCutHeavyMesonPCM("00000113","00200009227000002250400000","32c51070a","0103603500000000","0153503000000000");
    cuts.AddCutHeavyMesonPCM("00000113","00200009227000006250400000","32c51070a","0103603500000000","0153503000000000");
  } else if ( trainConfig == 146){ // Chi2GammaCut
    cuts.AddCutHeavyMesonPCM("00000113","00200009227000008150400000","32c51070a","0103603500000000","0153503000000000");
    cuts.AddCutHeavyMesonPCM("00000113","00200009227000008850400000","32c51070a","0103603500000000","0153503000000000");
    cuts.AddCutHeavyMesonPCM("00000113","00200009227000008a50400000","32c51070a","0103603500000000","0153503000000000");
    cuts.AddCutHeavyMesonPCM("00000113","00200009227000008950400000","32c51070a","0103603500000000","0153503000000000");
  } else if ( trainConfig == 147){ // PsiPair
    cuts.AddCutHeavyMesonPCM("00000113","00200009227000008260400000","32c51070a","0103603500000000","0153503000000000");
    cuts.AddCutHeavyMesonPCM("00000113","00200009227000008280400000","32c51070a","0103603500000000","0153503000000000");
    cuts.AddCutHeavyMesonPCM("00000113","00200009227000008210400000","32c51070a","0103603500000000","0153503000000000");
    cuts.AddCutHeavyMesonPCM("00000113","00200009227000008220400000","32c51070a","0103603500000000","0153503000000000");
  } else if ( trainConfig == 150){ // super loose track cuts
    cuts.AddCutHeavyMesonPCM("00000113","00200009227000008250400000","302010708","0103603500000000","0153503000000000"); // no ITS cut
  } else if ( trainConfig == 151){ // ITS clusterCuts track
    cuts.AddCutHeavyMesonPCM("00000113","00200009227000008250400000","302010708","0103603500000000","0153503000000000"); // no ITS cut
    cuts.AddCutHeavyMesonPCM("00000113","00200009227000008250400000","312010708","0103603500000000","0153503000000000"); // hot first layer of SPD
    cuts.AddCutHeavyMesonPCM("00000113","00200009227000008250400000","322010708","0103603500000000","0153503000000000"); // hit in any layer SPD and ITS refit
    cuts.AddCutHeavyMesonPCM("00000113","00200009227000008250400000","372010708","0103603500000000","0153503000000000"); // hit in any layer SPD min 1 cls in ITS total
  } else if ( trainConfig == 152){ // TPC cluster cut
    cuts.AddCutHeavyMesonPCM("00000113","00200009227000008250400000","302010708","0103603500000000","0153503000000000"); // min 80 clusters
    cuts.AddCutHeavyMesonPCM("00000113","00200009227000008250400000","30c010708","0103603500000000","0153503000000000"); // min 80 clusters + TPC refit
    cuts.AddCutHeavyMesonPCM("00000113","00200009227000008250400000","30a010708","0103603500000000","0153503000000000"); // min 80 clusters + TPC refit + chi2 etc
  } else if ( trainConfig == 153){ // DCA cut
    cuts.AddCutHeavyMesonPCM("00000113","00200009227000008250400000","30c510708","0103603500000000","0153503000000000"); // no DCA cut
  } else if ( trainConfig == 154){ // pt cut
    cuts.AddCutHeavyMesonPCM("00000113","00200009227000008250400000","30c030708","0103603500000000","0153503000000000"); // stricter pt cut of 150MeV like in AOD preselection
  } else if ( trainConfig == 155){ // mass cut use 4 vec
    cuts.AddCutHeavyMesonPCM("00000113","00200009227000008250400000","30201070a","0103603500000000","0153503000000000"); // no ITS cut
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //                                          ETA PRIME MESON
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  } else if( trainConfig == 200 ) {
    // everything open, min pt charged pi = 100 MeV
    cuts.AddCutHeavyMesonPCM("00000113","00200009327000008250400000","000010400","0103503m00000000","0103503000000000");
  } else if( trainConfig == 201 ) {
    // closing charged pion cuts, minimum TPC cluster = 80, TPC dEdx pi = \pm 3 sigma, min pt charged pi = 100 MeV
    cuts.AddCutHeavyMesonPCM("00000113","00200009327000008250400000","002010700","0103503m00000000","0103503000000000");
  } else if( trainConfig == 202) {
    // eta < 0.9
    // closing charged pion cuts, minimum TPC cluster = 80, TPC dEdx pi = \pm 3 sigma, pi+pi- mass cut of 1.5, min pt charged pi = 100 MeV
    // closing neural pion cuts, 0.5 < M_gamma,gamma < 0.6
    // maxChi2 per cluster TPC <4, require TPC refit, DCA XY pT dependend 0.0182+0.0350/pt^1.01, DCA_Z = 3.0
    // timing cluster cut open
    cuts.AddCutHeavyMesonPCM("00010113","00200009327000008250400000","30a330709","0103503l00000000","0153503000000000"); // all of the above
    cuts.AddCutHeavyMesonPCM("00052113","00200009327000008250400000","30a330709","0103503l00000000","0153503000000000"); // all of the above
    cuts.AddCutHeavyMesonPCM("00062113","00200009327000008250400000","30a330709","0103503l00000000","0153503000000000"); // all of the above
    cuts.AddCutHeavyMesonPCM("00083113","00200009327000008250400000","30a330709","0103503l00000000","0153503000000000"); // all of the above
    cuts.AddCutHeavyMesonPCM("00085113","00200009327000008250400000","30a330709","0103503l00000000","0153503000000000"); // all of the above
  } else if( trainConfig == 203) {
    // same as 202 but only MB
    cuts.AddCutHeavyMesonPCM("00010113","00200009327000008250400000","30a330709","0103503l00000000","0153503000000000"); // 0.5-0.6 eta mass cut
    cuts.AddCutHeavyMesonPCM("00010113","00200009327000008250400000","30a330709","0103503m00000000","0153503000000000"); // 0.4-0.7 eta mass cut
  } else if( trainConfig == 204) {
    // same as 202 but with mass cut variations
    cuts.AddCutHeavyMesonPCM("00010113","00200009327000008250400000","30a330700","0103503l00000000","0153503000000000"); // pi+pi- mass cut of 10
    cuts.AddCutHeavyMesonPCM("00010113","00200009327000008250400000","30a330701","0103503l00000000","0153503000000000"); // pi+pi- mass cut of 1
    cuts.AddCutHeavyMesonPCM("00010113","00200009327000008250400000","30a330708","0103503l00000000","0153503000000000"); // pi+pi- mass cut of 0.85
  } else if ( trainConfig == 205 ) { // min bias only
    cuts.AddCutHeavyMesonPCM("00010113","0dm00009f9730000dge0404000","32c51070m","0103603l00000000","0453503000000000"); // INT7
  } else if ( trainConfig == 206 ) { // no event mixing
    cuts.AddCutHeavyMesonPCM("00010113","00200009227000008250400000","32c510700","0103603l00000000","0453503000000000"); // INT7
    cuts.AddCutHeavyMesonPCM("0008e113","00200009227000008250400000","32c510700","0103603l00000000","0453503000000000"); // EG1
    cuts.AddCutHeavyMesonPCM("0008d113","00200009227000008250400000","32c510700","0103603l00000000","0453503000000000"); // EG2

    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //                                          OMEGA MESON (pp @ 13TeV)
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  } else if ( trainConfig == 401) { //Standard 13TeV
    cuts.AddCutHeavyMesonPCM("00010113","0dm00009f9730000dge0404000","32c51070a","0103103500000000","0153503000000000"); // V0AND
  } else if ( trainConfig == 402) { //Standard 13TeV + PHI7
    cuts.AddCutHeavyMesonPCM("00010113","0dm00009f9730000dge0404000","32c51070a","0103103500000000","0153503000000000"); // V0AND
    cuts.AddCutHeavyMesonPCM("00062113","0dm00009f9730000dge0404000","32c51070a","0103103500000000","0153503000000000"); // PHI7
  } else if ( trainConfig == 403) { //Standard 13TeV + PHOS Triggers
    cuts.AddCutHeavyMesonPCM("00010113","0dm00009f9730000dge0404000","32c51070a","0103103500000000","0153503000000000"); // V0AND
    cuts.AddCutHeavyMesonPCM("00061113","0dm00009f9730000dge0404000","32c51070a","0103103500000000","0153503000000000"); // PHI1
    cuts.AddCutHeavyMesonPCM("00062113","0dm00009f9730000dge0404000","32c51070a","0103103500000000","0153503000000000"); // PHI7
    cuts.AddCutHeavyMesonPCM("00063113","0dm00009f9730000dge0404000","32c51070a","0103103500000000","0153503000000000"); // PHI8
  } else if ( trainConfig == 404) { //Standard 13TeV, and some Triggers
    cuts.AddCutHeavyMesonPCM("00010113","0dm00009f9730000dge0404000","32c51070a","0103103500000000","0153503000000000"); // V0AND
    cuts.AddCutHeavyMesonPCM("00052113","0dm00009f9730000dge0404000","32c51070a","0103103500000000","0153503000000000"); // EMC7
    cuts.AddCutHeavyMesonPCM("00083113","0dm00009f9730000dge0404000","32c51070a","0103103500000000","0153503000000000"); // EG1
    cuts.AddCutHeavyMesonPCM("00085113","0dm00009f9730000dge0404000","32c51070a","0103103500000000","0153503000000000"); // EG2
    cuts.AddCutHeavyMesonPCM("00062113","0dm00009f9730000dge0404000","32c51070a","0103103500000000","0153503000000000"); // PHI7
  } else if ( trainConfig == 406) { //Standard 13TeV, and some more Triggers
    cuts.AddCutHeavyMesonPCM("00010113","0dm00009f9730000dge0404000","32c51070a","0103103500000000","0153503000000000"); // V0AND
    cuts.AddCutHeavyMesonPCM("00061113","0dm00009f9730000dge0404000","32c51070a","0103103500000000","0153503000000000"); // PHI1
    cuts.AddCutHeavyMesonPCM("00062113","0dm00009f9730000dge0404000","32c51070a","0103103500000000","0153503000000000"); // PHI7
    cuts.AddCutHeavyMesonPCM("00063113","0dm00009f9730000dge0404000","32c51070a","0103103500000000","0153503000000000"); // PHI8
    cuts.AddCutHeavyMesonPCM("0008e113","0dm00009f9730000dge0404000","32c51070a","0103103500000000","0153503000000000"); // PHI1
    cuts.AddCutHeavyMesonPCM("0008d113","0dm00009f9730000dge0404000","32c51070a","0103103500000000","0153503000000000"); // PHI7
    cuts.AddCutHeavyMesonPCM("0009b113","0dm00009f9730000dge0404000","32c51070a","0103103500000000","0153503000000000"); // PHI8
  } else if ( trainConfig == 416) { //Standard 13TeV, no shared TPC clusters
    cuts.AddCutHeavyMesonPCM("00010113","0dm00009f9730000dge0404000","32e51070a","0103103500000000","0153503000000000"); // V0AND
  } else if ( trainConfig == 420) { //Standard 13TeV, no shared TPC clusters
    cuts.AddCutHeavyMesonPCM("00010113","0dm00009f9730000dge0404000","32f51070a","0103103500000000","0153503000000000"); // V0AND

    // Variations on 5 TeV for 7 TeV systematics
  } else if ( trainConfig == 450) { //Standard 13TeV
    cuts.AddCutHeavyMesonPCM("00010113","00200009f9730000dge0400000","32c51070a","0103603a00000000","0153503000000000"); // V0AND
  } else if ( trainConfig == 451) { // mass window cut
    cuts.AddCutHeavyMesonPCM("00010113","00200009f9730000dge0400000","32c51070a","0103603n00000000","0153503000000000"); // 1 sigma
    cuts.AddCutHeavyMesonPCM("00010113","00200009f9730000dge0400000","32c51070a","0103603o00000000","0153503000000000"); // 3 sigma
    cuts.AddCutHeavyMesonPCM("00010113","00200009f9730000dge0400000","32c51070a","0103603p00000000","0153503000000000"); // 4 sigma
  } else if ( trainConfig == 452) { // background description
    cuts.AddCutHeavyMesonPCM("00010113","00200009f9730000dge0400000","32c51070a","0103603a00000000","0a53503000000000"); // likesign
    cuts.AddCutHeavyMesonPCM("00010113","00200009f9730000dge0400000","32c51070a","0103603a00000000","0d53503000000000"); // sideband

  } else if ( trainConfig == 701) { //Standard 13TeV
    cuts.AddCutHeavyMesonPCM("00010113","0dm00009f9730000dge0404000","32c51070a","0103103500000000","0a53503000000000"); // V0AND
  } else if ( trainConfig == 702) { //Standard 13TeV + PHI7
    cuts.AddCutHeavyMesonPCM("00010113","0dm00009f9730000dge0404000","32c51070a","0103103500000000","0a53503000000000"); // V0AND
    cuts.AddCutHeavyMesonPCM("00062113","0dm00009f9730000dge0404000","32c51070a","0103103500000000","0a53503000000000"); // PHI7
  } else if ( trainConfig == 703) { //Standard 13TeV + PHOS Triggers
    cuts.AddCutHeavyMesonPCM("00010113","0dm00009f9730000dge0404000","32c51070a","0103103500000000","0a53503000000000"); // V0AND
    cuts.AddCutHeavyMesonPCM("00061113","0dm00009f9730000dge0404000","32c51070a","0103103500000000","0a53503000000000"); // PHI1
    cuts.AddCutHeavyMesonPCM("00062113","0dm00009f9730000dge0404000","32c51070a","0103103500000000","0a53503000000000"); // PHI7
    cuts.AddCutHeavyMesonPCM("00063113","0dm00009f9730000dge0404000","32c51070a","0103103500000000","0a53503000000000"); // PHI8
  } else if ( trainConfig == 704) { //Standard 13TeV, and some Triggers
    cuts.AddCutHeavyMesonPCM("00010113","0dm00009f9730000dge0404000","32c51070a","0103103500000000","0a53503000000000"); // V0AND
    cuts.AddCutHeavyMesonPCM("00052113","0dm00009f9730000dge0404000","32c51070a","0103103500000000","0a53503000000000"); // EMC7
    cuts.AddCutHeavyMesonPCM("00083113","0dm00009f9730000dge0404000","32c51070a","0103103500000000","0a53503000000000"); // EG1
    cuts.AddCutHeavyMesonPCM("00085113","0dm00009f9730000dge0404000","32c51070a","0103103500000000","0a53503000000000"); // EG2
    cuts.AddCutHeavyMesonPCM("00062113","0dm00009f9730000dge0404000","32c51070a","0103103500000000","0a53503000000000"); // PHI7
  } else if ( trainConfig == 706) { //Standard 13TeV, and some more Triggers
    cuts.AddCutHeavyMesonPCM("00010113","0dm00009f9730000dge0404000","32c51070a","0103103500000000","0a53503000000000"); // V0AND
    cuts.AddCutHeavyMesonPCM("00061113","0dm00009f9730000dge0404000","32c51070a","0103103500000000","0a53503000000000"); // PHI1
    cuts.AddCutHeavyMesonPCM("00062113","0dm00009f9730000dge0404000","32c51070a","0103103500000000","0a53503000000000"); // PHI7
    cuts.AddCutHeavyMesonPCM("00063113","0dm00009f9730000dge0404000","32c51070a","0103103500000000","0a53503000000000"); // PHI8
    cuts.AddCutHeavyMesonPCM("0008e113","0dm00009f9730000dge0404000","32c51070a","0103103500000000","0a53503000000000"); // PHI1
    cuts.AddCutHeavyMesonPCM("0008d113","0dm00009f9730000dge0404000","32c51070a","0103103500000000","0a53503000000000"); // PHI7
    cuts.AddCutHeavyMesonPCM("0009b113","0dm00009f9730000dge0404000","32c51070a","0103103500000000","0a53503000000000"); // PHI8
  } else if ( trainConfig == 716) { //Standard 13TeV, no shared TPC clusters
    cuts.AddCutHeavyMesonPCM("00010113","0dm00009f9730000dge0404000","32e51070a","0103103500000000","0a53503000000000"); // V0AND
  } else if ( trainConfig == 720) { //Standard 13TeV, no shared TPC clusters
    cuts.AddCutHeavyMesonPCM("00010113","0dm00009f9730000dge0404000","32f51070a","0103103500000000","0a53503000000000"); // V0AND

    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // EMC pp 13 TeV Fitting, Systematics
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //Standard Cuts of Pi0 Analysis: ("00010113","0dm00009f9730000dge0404000","0r631031000000d0")
    //MesonCut r63==Background->ignored, d==OpeningAngle for Background->ignored =>0453503000000000
  } else if(trainConfig == 2000)  { //PCM INT7 Standard
    cuts.AddCutHeavyMesonPCM("00010113","0dm00009f9730000dge0404000","32c51070a","0103103500000000","0453503000000000"); // INT7 Standard
    //-----
    //INT7: Event Variations
    //-----
    //Std: 00010113
  } else if(trainConfig == 2001)  { //PCM INT7, Event cut var. Remove Pileup, Std 1-> True
    //                         00010113
    //                              |
    cuts.AddCutHeavyMesonPCM("00010013","0dm00009f9730000dge0404000","32c51070a","0103103500000000","0453503000000000"); // INT7 Pileup not removed
    //-----
    //INT7: Primary Pion / Charged Pion (Pi+ Pi-) Variations
    //-----
    //Std: 32c51070a
  } else if(trainConfig == 2201)  { //PCM INT7, Ch.Pi cut var. Ch.Pi ITS Requirement, Std 2 -> first or second SPD cluster required
    //                                                                32c51070a
    //                                                                 |
    cuts.AddCutHeavyMesonPCM("00010113","0dm00009f9730000dge0404000","30c51070a","0103103500000000","0453503000000000"); // INT7, Ch.Pi ITS, no SPD cluster requirement
    cuts.AddCutHeavyMesonPCM("00010113","0dm00009f9730000dge0404000","31c51070a","0103103500000000","0453503000000000"); // INT7, Ch.Pi ITS, first SPD cluster required
  } else if(trainConfig == 2202)  { //PCM INT7, Ch.Pi cut var. Ch.Pi Cls TPC, Std c -> MinClsTPC 80. + Refit
    //                                                                32c51070a
    //                                                                  |
    cuts.AddCutHeavyMesonPCM("00010113","0dm00009f9730000dge0404000","32251070a","0103103500000000","0453503000000000"); // INT7, Ch.Pi, MinClsTPC 80. + no Refit
    cuts.AddCutHeavyMesonPCM("00010113","0dm00009f9730000dge0404000","32g51070a","0103103500000000","0453503000000000"); // INT7, Ch.Pi, MinClsTPC 70. + Refit
    cuts.AddCutHeavyMesonPCM("00010113","0dm00009f9730000dge0404000","32h51070a","0103103500000000","0453503000000000"); // INT7, Ch.Pi, MinClsTPC 100. + Refit
    cuts.AddCutHeavyMesonPCM("00010113","0dm00009f9730000dge0404000","32i51070a","0103103500000000","0453503000000000"); // INT7, Ch.Pi, MinClsTPC 80. + Refit, 35% find. cls
    cuts.AddCutHeavyMesonPCM("00010113","0dm00009f9730000dge0404000","32j51070a","0103103500000000","0453503000000000"); // INT7, Ch.Pi, MinClsTPC 80. + Refit, 60% find. cls
  } else if(trainConfig == 2203)  { //PCM INT7, Ch.Pi cut var. Ch.Pi pT, Std 1 -> pt>0.1
    //                                                                32c51070a
    //                                                                    |
    cuts.AddCutHeavyMesonPCM("00010113","0dm00009f9730000dge0404000","32c50070a","0103103500000000","0453503000000000"); // INT7, Ch.Pi pt>0.075
    cuts.AddCutHeavyMesonPCM("00010113","0dm00009f9730000dge0404000","32c52070a","0103103500000000","0453503000000000"); // INT7, Ch.Pi pt>0.125
    cuts.AddCutHeavyMesonPCM("00010113","0dm00009f9730000dge0404000","32c53070a","0103103500000000","0453503000000000"); // INT7, Ch.Pi pt>0.15
    cuts.AddCutHeavyMesonPCM("00010113","0dm00009f9730000dge0404000","32c54070a","0103103500000000","0453503000000000"); // INT7, Ch.Pi pt>0.4
  } else if(trainConfig == 2204)  { //PCM INT7, Ch.Pi cut var. Ch.Pi TPC dEdx, Std 7 -> -3,3
    //                                                                32c51070a
    //                                                                      |
    cuts.AddCutHeavyMesonPCM("00010113","0dm00009f9730000dge0404000","32c51050a","0103103500000000","0453503000000000"); // INT7, Ch.Pi -4,4
    cuts.AddCutHeavyMesonPCM("00010113","0dm00009f9730000dge0404000","32c51080a","0103103500000000","0453503000000000"); // INT7, Ch.Pi -2,3
    cuts.AddCutHeavyMesonPCM("00010113","0dm00009f9730000dge0404000","32c51020a","0103103500000000","0453503000000000"); // INT7, Ch.Pi -6,7
    cuts.AddCutHeavyMesonPCM("00010113","0dm00009f9730000dge0404000","32c51030a","0103103500000000","0453503000000000"); // INT7, Ch.Pi -5,5
    cuts.AddCutHeavyMesonPCM("00010113","0dm00009f9730000dge0404000","32c51040a","0103103500000000","0453503000000000"); // INT7, Ch.Pi -4,5
  } else if(trainConfig == 2205)  { //PCM INT7, Ch.Pi cut var. Ch.Pi Mass, Std a -> Ch.Pi<850MeV
    //                                                                32c51070a
    //                                                                        |
    cuts.AddCutHeavyMesonPCM("00010113","0dm00009f9730000dge0404000","32c51070f","0103103500000000","0453503000000000"); // INT7, Ch.Pi<650MeV
    cuts.AddCutHeavyMesonPCM("00010113","0dm00009f9730000dge0404000","32c51070c","0103103500000000","0453503000000000"); // INT7, Ch.Pi<850MeV, Neu.Pi<1000MeV
    cuts.AddCutHeavyMesonPCM("00010113","0dm00009f9730000dge0404000","32c51070h","0103103500000000","0453503000000000"); // INT7, Ch.Pi<650MeV, Neu.Pi<1000MeV
    //-----
    //INT7: Neutral Meson (Pi0) Cut Variations
    //-----
    //Std: 0103103500000000
  } else if(trainConfig == 2302)  { //PCM INT7, N.Pi cut var. rapidity, Std 1 -> -0.8, 0.8
    //                                                                            0103103500000000
    //                                                                                |
    cuts.AddCutHeavyMesonPCM("00010113","0dm00009f9730000dge0404000","32c51070a","0103503500000000","0453503000000000"); // INT7, N.Pi rap. -0.85, 0.85
    cuts.AddCutHeavyMesonPCM("00010113","0dm00009f9730000dge0404000","32c51070a","0103603500000000","0453503000000000"); // INT7, N.Pi rap. -0.75, 0.75
  } else if(trainConfig == 2304)  { //PCM INT7, N.Pi cut var. alpha, Std 3 -> 0.0-1.0
    //                                                                            0103103500000000
    //                                                                                  |
    cuts.AddCutHeavyMesonPCM("00010113","0dm00009f9730000dge0404000","32c51070a","0103105100000000","0453503000000000"); // INT7 alpha 0-0.75
    cuts.AddCutHeavyMesonPCM("00010113","0dm00009f9730000dge0404000","32c51070a","0103108100000000","0453503000000000"); // INT7 alpha 0-0.6
    cuts.AddCutHeavyMesonPCM("00010113","0dm00009f9730000dge0404000","32c51070a","0103101100000000","0453503000000000"); // alpha meson pT dependent
  } else if(trainConfig == 2305)  { //PCM INT7, N.Pi cut var. Selection Window, Std 1 -> 2 sigma
    //                                                                            0103103500000000
    //                                                                                   |
    cuts.AddCutHeavyMesonPCM("00010113","0dm00009f9730000dge0404000","32c51070a","0103103u00000000","0453503000000000"); // INT7, 1 sigma, NEEDS CHANGE AS THIS IS EMCEMC
    cuts.AddCutHeavyMesonPCM("00010113","0dm00009f9730000dge0404000","32c51070a","0103103v00000000","0453503000000000"); // INT7, 3 sigma, NEEDS CHANGE AS THIS IS EMCEMC
  } else if(trainConfig == 2306)  { //PCM INT7, N.Pi cut var. open. angle, Std 0 -> off
    //                                                                            0103103500000000
    //                                                                                          |
    cuts.AddCutHeavyMesonPCM("00010113","0dm00009f9730000dge0404000","32c51070a","01031031000000d0","0453503000000000"); // INT7 Op. Ang. var 1 cell dist + 0.017
    cuts.AddCutHeavyMesonPCM("00010113","0dm00009f9730000dge0404000","32c51070a","01031031000000b0","0453503000000000"); // INT7 Op. Ang. var 1 cell dist + 0.0152
    cuts.AddCutHeavyMesonPCM("00010113","0dm00009f9730000dge0404000","32c51070a","01031031000000g0","0453503000000000"); // INT7 Op. Ang. var 1 cell dist + 0.0202
    cuts.AddCutHeavyMesonPCM("00010113","0dm00009f9730000dge0404000","32c51070a","01031031000000a0","0453503000000000"); // INT7 Op. Ang. var 1 cell dist + 0

    //-----
    //INT7: Omega Meson Cut Variations
    //-----
    //Std: 0453503000000000
  } else if(trainConfig == 2401)  { //PCM INT7, Omega cut var. Background Scheme, Std 4 -> off
    //                                                                                               0453503000000000
    //                                                                                                |
    cuts.AddCutHeavyMesonPCM("00010113","0dm00009f9730000dge0404000","32c51070a","0103103500000000","0153103000000000"); // INT7, Om Event Mixing
    cuts.AddCutHeavyMesonPCM("00010113","0dm00009f9730000dge0404000","32c51070a","0103103500000000","0a53603000000000"); // INT7, Om LikeSignMixing
    cuts.AddCutHeavyMesonPCM("00010113","0dm00009f9730000dge0404000","32c51070a","0103103500000000","0d53603000000000"); // INT7, Om SideBandMixing
  } else if(trainConfig == 2402)  { //PCM INT7, Omega cut var. rapidity, Std 5 -> -0.85, 0.85
    //                                                                                               0453503000000000
    //                                                                                                   |
    cuts.AddCutHeavyMesonPCM("00010113","0dm00009f9730000dge0404000","32c51070a","0103103500000000","0453103000000000"); // INT7, Om rap. -0.8, 0.8
    cuts.AddCutHeavyMesonPCM("00010113","0dm00009f9730000dge0404000","32c51070a","0103103500000000","0453603000000000"); // INT7, Om rap. -0.75, 0.75
  } else if(trainConfig == 2404)  { //PCM INT7, Omega cut var. alpha, Std 3 -> 0.0-1.0
    //                                                                                               0453503000000000
    //                                                                                                     |
    cuts.AddCutHeavyMesonPCM("00010113","0dm00009f9730000dge0404000","32c51070a","0103103500000000","0453505000000000"); // INT7 alpha 0-0.75
    cuts.AddCutHeavyMesonPCM("00010113","0dm00009f9730000dge0404000","32c51070a","0103103500000000","0453508000000000"); // INT7 alpha 0-0.6
  } else if(trainConfig == 2410)  { //PCM INT7, Omega cut var. Background Scheme single cfg, Std 4 -> off, EventMixing
    //                                                                                               0453503000000000
    //                                                                                                |
    cuts.AddCutHeavyMesonPCM("00010113","0dm00009f9730000dge0404000","32c51070a","0103103500000000","0153103000000000"); // INT7, Om Event Mixing
  } else if(trainConfig == 2411)  { //PCM INT7, Omega cut var. Background Scheme single cfg, Std 4 -> off, LikeSignMixing
    //                                                                                               0453503000000000
    //                                                                                                |
    cuts.AddCutHeavyMesonPCM("00010113","0dm00009f9730000dge0404000","32c51070a","0103103500000000","0a53603000000000"); // INT7, Om LikeSignMixing
  } else if(trainConfig == 2412)  { //PCM INT7, Omega cut var. Background Scheme single cfg, Std 4 -> off, SideBandMixing
    //                                                                                               0453503000000000
    //                                                                                                |
    cuts.AddCutHeavyMesonPCM("00010113","0dm00009f9730000dge0404000","32c51070a","0103103500000000","0d53603000000000"); // INT7, Om SideBandMixing

    //-----
    //INT7: PCM Conversion Cut
    //-----
    //Std: 0dm00009f9730000dge0404000
  } else if (trainConfig == 2501) {   // min pT variations
    //                                    0dm00009f9730000dge0404000
    //                                          |
    cuts.AddCutHeavyMesonPCM("00010113", "0dm00069f9730000dge0404000", "32c51070a", "0103103500000000", "0453503000000000"); // eta < 0.8  // remove  55-72 bin, min pT 40 MeV
    cuts.AddCutHeavyMesonPCM("00010113", "0dm00049f9730000dge0404000", "32c51070a", "0103103500000000", "0453503000000000"); // eta < 0.8  // remove  55-72 bin, min pT 75 MeV
    cuts.AddCutHeavyMesonPCM("00010113", "0dm00019f9730000dge0404000", "32c51070a", "0103103500000000", "0453503000000000"); // eta < 0.8  // remove  55-72 bin, min pT 100MeV

  } else if (trainConfig == 2502) {   // TPC clusters, cosPA
    //                                    0dm00009f9730000dge0404000
    //                                           |
    cuts.AddCutHeavyMesonPCM("00010113", "0dm00008f9730000dge0404000", "32c51070a", "0103103500000000", "0453503000000000"); // TPC cluster 35%
    cuts.AddCutHeavyMesonPCM("00010113", "0dm00006f9730000dge0404000", "32c51070a", "0103103500000000", "0453503000000000"); // TPC cluster 70%
    cuts.AddCutHeavyMesonPCM("00010113", "0dm00009f9730000dge0604000", "32c51070a", "0103103500000000", "0453503000000000"); // cosPA 0.9
    cuts.AddCutHeavyMesonPCM("00010113", "0dm00009f9730000dge0304000", "32c51070a", "0103103500000000", "0453503000000000"); // cosPA 0.75

  } else if (trainConfig == 2503) {   // TPC clusters, cosPA
    //                                    0dm00009f9730000dge0404000
    //                                            |
    cuts.AddCutHeavyMesonPCM("00010113", "0dm0000939730000dge0404000", "32c51070a", "0103103500000000", "0453503000000000"); // nsig electron   -4,5
    cuts.AddCutHeavyMesonPCM("00010113", "0dm0000969730000dge0404000", "32c51070a", "0103103500000000", "0453503000000000"); // nsig electron -2.5,4
    cuts.AddCutHeavyMesonPCM("00010113", "0dm00009f5730000dge0404000", "32c51070a", "0103103500000000", "0453503000000000"); // nsig pion 2,-10
    cuts.AddCutHeavyMesonPCM("00010113", "0dm00009f1730000dge0404000", "32c51070a", "0103103500000000", "0453503000000000"); // nsig pion 0,-10

  } else if (trainConfig == 2504) {
    //                                    0dm00009f9730000dge0404000
    //                                              |
    cuts.AddCutHeavyMesonPCM("00010113", "0dm00009f9030000dge0404000", "32c51070a", "0103103500000000", "0453503000000000"); // pion nsig min mom 0.50 GeV/c
    cuts.AddCutHeavyMesonPCM("00010113", "0dm00009f9630000dge0404000", "32c51070a", "0103103500000000", "0453503000000000"); // pion nsig min mom 0.25 GeV/c
    cuts.AddCutHeavyMesonPCM("00010113", "0dm00009f9760000dge0404000", "32c51070a", "0103103500000000", "0453503000000000"); // pion nsig max mom 2.00 GeV/c
    cuts.AddCutHeavyMesonPCM("00010113", "0dm00009f9710000dge0404000", "32c51070a", "0103103500000000", "0453503000000000"); // pion nsig max mom 5.00 GeV/c


  } else if (trainConfig == 2505) {   // chi2 variations
    //                                    0dm00009f9730000dge0404000
    //                                                     ||
    cuts.AddCutHeavyMesonPCM("00010113", "0dm00009f9730000d1e0404000", "32c51070a", "0103103500000000", "0453503000000000"); // chi2 50
    cuts.AddCutHeavyMesonPCM("00010113", "0dm00009f9730000dfe0404000", "32c51070a", "0103103500000000", "0453503000000000"); // chi2 50 chi2 dep -0.065
    cuts.AddCutHeavyMesonPCM("00010113", "0dm00009f9730000dhe0404000", "32c51070a", "0103103500000000", "0453503000000000"); // chi2 50 chi2 dep -0.050
    cuts.AddCutHeavyMesonPCM("00010113", "0dm00009f9730000dge0400000", "32c51070a", "0103103500000000", "0453503000000000"); // remove reject close v0

  } else if (trainConfig == 2506) {   // Psi pair variations
    //                                    0dm00009f9730000dge0404000
    //                                                     ||
    cuts.AddCutHeavyMesonPCM("00010113", "0dm00009f9730000dgd0404000", "32c51070a", "0103103500000000", "0453503000000000"); // Psi pair 0.15 dep
    cuts.AddCutHeavyMesonPCM("00010113", "0dm00009f9730000dgf0404000", "32c51070a", "0103103500000000", "0453503000000000"); // Psi pair 0.20 dep
    cuts.AddCutHeavyMesonPCM("00010113", "0dm00009f9730000dgg0404000", "32c51070a", "0103103500000000", "0453503000000000"); // Psi pair 0.30 dep
    cuts.AddCutHeavyMesonPCM("00010113", "0dm00009227300008250404000", "32c51070a", "0103103500000000", "0453503000000000"); // old cuts (run1)


  } else if ( trainConfig == 2511) { //Standard 13TeV, Material Budget Studies
    cuts.AddCutHeavyMesonPCM("00010113", "0dm00009f9730000dge0404000", "32c51070a", "0103103500000000", "0453503000000000");

  } else if ( trainConfig == 2512) { //Standard 13TeV, Material Budget Studies
    cuts.AddCutHeavyMesonPCM("00010113", "0dm00009f9730000dge0404000", "32c51070a", "0103103500000000", "0453503000000000");

  } else if ( trainConfig == 2513) { //Standard 13TeV, Material Budget Studies
    cuts.AddCutHeavyMesonPCM("00010113", "0dm00009f9730000dge0404000", "32c51070a", "0103103500000000", "0453503000000000");

    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  } else {
    Error(Form("GammaConvNeutralMeson_ConvMode_%i",trainConfig), "wrong trainConfig variable no cuts have been specified for the configuration");
    return;
  }

  if(!cuts.AreValid()){
    cout << "\n\n****************************************************" << endl;
    cout << "ERROR: No valid cuts stored in CutHandlerNeutralConv! Returning..." << endl;
    cout << "****************************************************\n\n" << endl;
    return;
  }

  Int_t numberOfCuts = cuts.GetNCuts();

  TList *EventCutList = new TList();
  TList *ConvCutList  = new TList();
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
  NeutralPionCutList->SetOwner(kTRUE);
  AliConversionMesonCuts **analysisNeutralPionCuts   = new AliConversionMesonCuts*[numberOfCuts];
  MesonCutList->SetOwner(kTRUE);
  AliConversionMesonCuts **analysisMesonCuts   = new AliConversionMesonCuts*[numberOfCuts];
  PionCutList->SetOwner(kTRUE);
  AliPrimaryPionCuts **analysisPionCuts     = new AliPrimaryPionCuts*[numberOfCuts];
  Bool_t initializedMatBudWeigths_existing    = kFALSE;

  for(Int_t i = 0; i<numberOfCuts; i++){
    analysisEventCuts[i] = new AliConvEventCuts();
    analysisEventCuts[i]->SetV0ReaderName(V0ReaderName);
    analysisEventCuts[i]->SetTriggerMimicking(enableTriggerMimicking);
    if(fileNameCustomTriggerMimicOADB.CompareTo("") != 0)
      analysisEventCuts[i]->SetCustomTriggerMimicOADBFile(fileNameCustomTriggerMimicOADB);
    analysisEventCuts[i]->SetTriggerOverlapRejecion(enableTriggerOverlapRej);
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
          cout << "MBW properly initialized" << endl;
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

    TString cutName( Form("%s_%s_%s_%s_%s",(cuts.GetEventCut(i)).Data(), (cuts.GetPhotonCut(i)).Data(),(cuts.GetPionCut(i)).Data(),(cuts.GetNDMCut(i)).Data(), (cuts.GetMesonCut(i)).Data() ) );
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
  task->SetNeutralPionCutList(NeutralPionCutList);
  task->SetMesonCutList(MesonCutList);
  task->SetPionCutList(PionCutList);

  task->SetMoveParticleAccordingToVertex(kTRUE);

  task->SetSelectedHeavyNeutralMeson(selectHeavyNeutralMeson);

  task->SetDoMesonQA(enableQAMesonTask);
  if (initializedMatBudWeigths_existing) {
      task->SetDoMaterialBudgetWeightingOfGammasForTrueMesons(kTRUE);
      if (enableMatBudWeightsPi0>=10){
          task->SetDoMaterialBudgetWeightingOfGammasForInvMassHistogram(kTRUE);
      }
  }

  //connect containers
  AliAnalysisDataContainer *coutput =
  mgr->CreateContainer( (usePreSelection) ? Form("GammaConvNeutralMesonPiPlPiMiNeutralMeson_%i_%i_%i.root",selectHeavyNeutralMeson,neutralPionMode, trainConfig) : Form("GammaConvNeutralMesonPiPlPiMiNeutralMeson_%i_%i_%i_PreSel%s.root",selectHeavyNeutralMeson,neutralPionMode, trainConfig, PionCuts.Data()), TList::Class(),
              AliAnalysisManager::kOutputContainer,Form("GammaConvNeutralMesonPiPlPiMiNeutralMeson_%i_%i_%i.root",selectHeavyNeutralMeson,neutralPionMode, trainConfig));

  mgr->AddTask(task);
  mgr->ConnectInput(task,0,cinput);
  mgr->ConnectOutput(task,1,coutput);

  return;

}
