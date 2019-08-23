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
    Int_t     trainConfig                 = 1,
    Int_t     isMC                        = 0,                        //run MC
    TString   photonCutNumberV0Reader     = "",                       // 00000008400000000100000000 nom. B, 00000088400000000100000000 low B
    Int_t     selectHeavyNeutralMeson     = 0,                        //run eta prime instead of omega
    Int_t     enableQAMesonTask           = 1,                        //enable QA in AliAnalysisTaskNeutralMesonToPiPlPiMiNeutralMeson
    Bool_t    enableTriggerMimicking      = kFALSE,                   // enable trigger mimicking
    Bool_t    enableTriggerOverlapRej     = kFALSE,                   // enable trigger overlap rejection
    TString   fileNameInputForWeighting   = "MCSpectraInput.root",    // path to file for weigting input
    Bool_t    doWeighting                 = kFALSE,                   //enable Weighting
    TString   generatorName               = "HIJING",
    Double_t  tolerance                   = -1,
    TString   periodNameV0Reader          = "",                       // period Name for V0Reader
    Int_t     runLightOutput              = 0,                        // run light output option 0: no light output 1: most cut histos stiched off 2: unecessary omega hists turned off as well
    TString   additionalTrainConfig       = "0"                       // additional counter for trainconfig, this has to be always the last parameter
  ) {

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

  if( !(AliPrimaryPionSelector*)mgr->GetTask("PionSelector") ){
    AliPrimaryPionSelector *fPionSelector = new AliPrimaryPionSelector("PionSelector");
    AliPrimaryPionCuts *fPionCuts=0;
    if( PionCuts!=""){
      fPionCuts= new AliPrimaryPionCuts(PionCuts.Data(),PionCuts.Data());
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
  task->SetIsHeavyIon(isHeavyIon);
  task->SetIsMC(isMC);
  task->SetV0ReaderName(V0ReaderName);
  if(runLightOutput>1) task->SetLightOutput(kTRUE);
  task->SetTolerance(tolerance);
  AliCutHandlerPCM cuts(13);

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
    cuts.AddCutHeavyMesonPCM("00000113","00200009227000008250400000","32c510708","0103603500000000","0153503000000000"); // std
    cuts.AddCutHeavyMesonPCM("00000113","00200009227000008250400000","32c51070a","0103603500000000","0153503000000000"); // use 4 vec
  } else if ( trainConfig == 121) { // pTCut
    cuts.AddCutHeavyMesonPCM("00000113","00200009227000008250400000","32c50070a","0103603500000000","0153503000000000"); // pt>0.075
    cuts.AddCutHeavyMesonPCM("00000113","00200009227000008250400000","32c52070a","0103603500000000","0153503000000000"); // pt>0.125
    cuts.AddCutHeavyMesonPCM("00000113","00200009227000008250400000","32c53070a","0103603500000000","0153503000000000"); // pt>0.15
    cuts.AddCutHeavyMesonPCM("00000113","00200009227000008250400000","32c54070a","0103603500000000","0153503000000000"); // pt>0.4
  } else if ( trainConfig == 122) { // TPCdEdxCutPion
    cuts.AddCutHeavyMesonPCM("00000113","00200009227000008250400000","32c51050a","0103603500000000","0153503000000000"); // -4,5
    cuts.AddCutHeavyMesonPCM("00000113","00200009227000008250400000","32c51080a","0103603500000000","0153503000000000"); // -2.5,4
    cuts.AddCutHeavyMesonPCM("00000113","00200009227000008250400000","32c51020a","0103603500000000","0153503000000000"); // -6,7
    cuts.AddCutHeavyMesonPCM("00000113","00200009227000008250400000","32c51030a","0103603500000000","0153503000000000"); // -4,4
    cuts.AddCutHeavyMesonPCM("00000113","00200009227000008250400000","32c51040a","0103603500000000","0153503000000000"); // -2.5,4
    // Neutral Pion Cuts
  } else if ( trainConfig == 123) { // invariant mass cut (std 110-150)
    cuts.AddCutHeavyMesonPCM("00000113","00200009227000008250400000","32c51070a","0103603100000000","0153503000000000"); // 100-145
    cuts.AddCutHeavyMesonPCM("00000113","00200009227000008250400000","32c51070a","0103603200000000","0153503000000000"); // 110-145
    cuts.AddCutHeavyMesonPCM("00000113","00200009227000008250400000","32c51070a","0103603300000000","0153503000000000"); // 120-145
    cuts.AddCutHeavyMesonPCM("00000113","00200009227000008250400000","32c51070a","0103603800000000","0153503000000000"); // 125-145
    cuts.AddCutHeavyMesonPCM("00000113","00200009227000008250400000","32c51070a","0103603a00000000","0153503000000000"); //  80-145
  } else if ( trainConfig == 130) {
    cuts.AddCutHeavyMesonPCM("00000113","00200009227000008250400000","32b51070a","0103603500000000","0453503000000000"); // without background

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
    
    // pp 13 TeV
  } else if ( trainConfig == 401) { //Standard 13TeV
    cuts.AddCutHeavyMesonPCM("00010113","0d200009327000008250404000","32c51070a","0103603s00000000","0153503000000000"); // V0AND
  } else if ( trainConfig == 402) { //Standard 13TeV + PHI7
    cuts.AddCutHeavyMesonPCM("00010113","0d200009327000008250404000","32c51070a","0103603s00000000","0153503000000000"); // V0AND
    cuts.AddCutHeavyMesonPCM("00062113","0d200009327000008250404000","32c51070a","0103603s00000000","0153503000000000"); // PHI7
  } else if ( trainConfig == 403) { //Standard 13TeV + PHOS Triggers
    cuts.AddCutHeavyMesonPCM("00010113","0d200009327000008250404000","32c51070a","0103603s00000000","0153503000000000"); // V0AND
    cuts.AddCutHeavyMesonPCM("00061113","0d200009327000008250404000","32c51070a","0103603s00000000","0153503000000000"); // PHI1
    cuts.AddCutHeavyMesonPCM("00062113","0d200009327000008250404000","32c51070a","0103603s00000000","0153503000000000"); // PHI7
    cuts.AddCutHeavyMesonPCM("00063113","0d200009327000008250404000","32c51070a","0103603s00000000","0153503000000000"); // PHI8
  } else if ( trainConfig == 404) { //Standard 13TeV, and some Triggers
    cuts.AddCutHeavyMesonPCM("00010113","0d200009327000008250404000","32c51070a","0103603500000000","0153503000000000"); // V0AND
    cuts.AddCutHeavyMesonPCM("00052113","0d200009327000008250404000","32c51070a","0103603500000000","0153503000000000"); // EMC7
    cuts.AddCutHeavyMesonPCM("00083113","0d200009327000008250404000","32c51070a","0103603500000000","0153503000000000"); // EG1
    cuts.AddCutHeavyMesonPCM("00085113","0d200009327000008250404000","32c51070a","0103603500000000","0153503000000000"); // EG2
    cuts.AddCutHeavyMesonPCM("00062113","0d200009327000008250404000","32c51070a","0103603500000000","0153503000000000"); // PHI7
  } else if ( trainConfig == 406) { //Standard 13TeV, and some more Triggers
    cuts.AddCutHeavyMesonPCM("00010113","0d200009327000008250404000","32c51070a","0103603s00000000","0153503000000000"); // V0AND
    cuts.AddCutHeavyMesonPCM("00061113","0d200009327000008250404000","32c51070a","0103603s00000000","0153503000000000"); // PHI1
    cuts.AddCutHeavyMesonPCM("00062113","0d200009327000008250404000","32c51070a","0103603s00000000","0153503000000000"); // PHI7
    cuts.AddCutHeavyMesonPCM("00063113","0d200009327000008250404000","32c51070a","0103603s00000000","0153503000000000"); // PHI8
    cuts.AddCutHeavyMesonPCM("0008e113","0d200009327000008250404000","32c51070a","0103603s00000000","0153503000000000"); // PHI1
    cuts.AddCutHeavyMesonPCM("0008d113","0d200009327000008250404000","32c51070a","0103603s00000000","0153503000000000"); // PHI7
    cuts.AddCutHeavyMesonPCM("0009b113","0d200009327000008250404000","32c51070a","0103603s00000000","0153503000000000"); // PHI8
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

  for(Int_t i = 0; i<numberOfCuts; i++){
    analysisEventCuts[i] = new AliConvEventCuts();
    analysisEventCuts[i]->SetV0ReaderName(V0ReaderName);
    analysisEventCuts[i]->SetTriggerMimicking(enableTriggerMimicking);
    analysisEventCuts[i]->SetTriggerOverlapRejecion(enableTriggerOverlapRej);
    if(runLightOutput>0) analysisEventCuts[i]->SetLightOutput(kTRUE);
    analysisEventCuts[i]->InitializeCutsFromCutString((cuts.GetEventCut(i)).Data());
    if (periodNameV0Reader.CompareTo("") != 0) analysisEventCuts[i]->SetPeriodEnum(periodNameV0Reader);
    EventCutList->Add(analysisEventCuts[i]);
    analysisEventCuts[i]->SetFillCutHistograms("",kFALSE);

    analysisCuts[i] = new AliConversionPhotonCuts();
    if(runLightOutput>0) analysisCuts[i]->SetLightOutput(kTRUE);
    analysisCuts[i]->SetV0ReaderName(V0ReaderName);
    if( ! analysisCuts[i]->InitializeCutsFromCutString((cuts.GetPhotonCut(i)).Data()) ) {
      cout<<"ERROR: analysisCuts [" <<i<<"]"<<endl;
      return 0;
    } else {
      ConvCutList->Add(analysisCuts[i]);
      analysisCuts[i]->SetFillCutHistograms("",kFALSE);

    }

    analysisNeutralPionCuts[i] = new AliConversionMesonCuts();
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

  //connect containers
  AliAnalysisDataContainer *coutput =
  mgr->CreateContainer(Form("GammaConvNeutralMesonPiPlPiMiNeutralMeson_%i_%i_%i.root",selectHeavyNeutralMeson,neutralPionMode, trainConfig), TList::Class(),
              AliAnalysisManager::kOutputContainer,Form("GammaConvNeutralMesonPiPlPiMiNeutralMeson_%i_%i_%i.root",selectHeavyNeutralMeson,neutralPionMode, trainConfig));

  mgr->AddTask(task);
  mgr->ConnectInput(task,0,cinput);
  mgr->ConnectOutput(task,1,coutput);

  return;

}
