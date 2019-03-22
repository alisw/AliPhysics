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
    Int_t trainConfig                 = 1,
    Int_t isMC                        = 0,                                  //run MC
    TString   photonCutNumberV0Reader       = "",       // 00000008400000000100000000 nom. B, 00000088400000000100000000 low B
    Int_t selectHeavyNeutralMeson     = 0,                                  //run eta prime instead of omega
    Int_t enableQAMesonTask           = 1,                                  //enable QA in AliAnalysisTaskNeutralMesonToPiPlPiMiNeutralMeson
    TString fileNameInputForWeighting = "MCSpectraInput.root",              // path to file for weigting input
    Bool_t doWeighting                = kFALSE,                             //enable Weighting
    TString generatorName             = "HIJING",
    Double_t tolerance                = -1,
    TString periodNameV0Reader        = "",                                 // period Name for V0Reader
    Int_t runLightOutput              = 0,                                  // run light output option 0: no light output 1: most cut histos stiched off 2: unecessary omega hists turned off as well
    TString additionalTrainConfig     = "0"                                 // additional counter for trainconfig, this has to be always the last parameter
  ) {

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

  Int_t isHeavyIon = 0;
  Int_t neutralPionMode = 2;

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
  //================================================
  //========= Add Pion Selector ====================
  if( !(AliPrimaryPionSelector*)mgr->GetTask("PionSelector") ){
    AliPrimaryPionSelector *fPionSelector = new AliPrimaryPionSelector("PionSelector");
    AliPrimaryPionCuts *fPionCuts=0;
    if( PionCuts!=""){
      fPionCuts= new AliPrimaryPionCuts(PionCuts.Data(),PionCuts.Data());
      if(runLightOutput>0) fPionCuts->SetLightOutput(kTRUE);
      //if(runLightOutput>0) fPionCuts->SetLightOutput(kTRUE);
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


  AliCutHandlerPCM cuts(13);


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
    cuts.AddCutHeavyMesonCalo("00000113","1111111047032230000","302010708","0103603700000000","0153503000000000");
    cuts.AddCutHeavyMesonCalo("00000113","1111111047032230000","322010708","0103603700000000","0153503000000000"); // with ITS requirement

    // EMCal pp 5 TeV
  } else if( trainConfig == 111)  { // Test for EMCal (5 TeV) with TPC refit + ITS requirement
    cuts.AddCutHeavyMesonCalo("00010113","1111111047032230000","32c010708","0103603700000000","0153503000000000"); // INT7
    cuts.AddCutHeavyMesonCalo("00052113","1111111047032230000","32c010708","0103603700000000","0153503000000000"); // EMC7
    cuts.AddCutHeavyMesonCalo("00085113","1111111047032230000","32c010708","0103603700000000","0153503000000000"); // EG2
    cuts.AddCutHeavyMesonCalo("00083113","1111111047032230000","32c010708","0103603700000000","0153503000000000"); // EG1
    // EMCal pp 13 TeV
  } else if( trainConfig == 112)  { // Test for EMCal (13 TeV) with TPC refit + ITS requirement
    cuts.AddCutHeavyMesonCalo("00010113","1111111047032230000","32c010708","0103603700000000","0153503000000000"); // INT7
    cuts.AddCutHeavyMesonCalo("00052113","1111111047032230000","32c010708","0103603700000000","0153503000000000"); // EMC7
    cuts.AddCutHeavyMesonCalo("00085113","1111111047032230000","32c010708","0103603700000000","0153503000000000"); // EG2
    cuts.AddCutHeavyMesonCalo("00083113","1111111047032230000","32c010708","0103603700000000","0153503000000000"); // EG1
    // EMCal LHC11 pp 7TeV
  } else if( trainConfig == 115){ // with TPC refit + ITS requirement
    cuts.AddCutHeavyMesonCalo("00010113","1111111057032230000","32c010708","0103603700000000","0153503000000000"); // INT7
    cuts.AddCutHeavyMesonCalo("00052113","1111111057032230000","32c010708","0103603700000000","0153503000000000"); // EMC7
    // ---------------------------------
    // systematic studies 7 TeV (EMCal)
    // ---------------------------------

    // charged pion cuts
  } else if( trainConfig == 120)   {
    cuts.AddCutHeavyMesonCalo("00000113","1111111047032230000","30c010708","0103603700000000","0153503000000000"); // with TPC refit (new standard)
    cuts.AddCutHeavyMesonCalo("00000113","1111111047032230000","322010708","0103603700000000","0153503000000000"); // with ITS requirement
    cuts.AddCutHeavyMesonCalo("00000113","1111111047032230000","32c010708","0103603700000000","0153503000000000"); // with TPC refit + ITS requirement
  } else if (trainConfig == 121) { // pT Cut
    cuts.AddCutHeavyMesonCalo("00000113","1111111047032230000","30c000708","0103603700000000","0153503000000000"); // pt>0.075
    cuts.AddCutHeavyMesonCalo("00000113","1111111047032230000","30c020708","0103603700000000","0153503000000000"); // pt>0.125
    cuts.AddCutHeavyMesonCalo("00000113","1111111047032230000","30c030708","0103603700000000","0153503000000000"); // pt>0.15
    cuts.AddCutHeavyMesonCalo("00000113","1111111047032230000","30c040708","0103603700000000","0153503000000000"); // pt>0.4
  } else if (trainConfig == 122) { // TPCdEdxCutPion
    cuts.AddCutHeavyMesonCalo("00000113","1111111047032230000","30c010508","0103603700000000","0153503000000000"); // -4,4
    cuts.AddCutHeavyMesonCalo("00000113","1111111047032230000","30c010808","0103603700000000","0153503000000000"); // -2,3
    cuts.AddCutHeavyMesonCalo("00000113","1111111047032230000","30c010208","0103603700000000","0153503000000000"); // -6,7
    cuts.AddCutHeavyMesonCalo("00000113","1111111047032230000","30c010308","0103603700000000","0153503000000000"); // -5,5
    cuts.AddCutHeavyMesonCalo("00000113","1111111047032230000","30c010408","0103603700000000","0153503000000000"); // -4,5
  // neutral pion cuts
  } else if (trainConfig == 123) { // Mass window
    cuts.AddCutHeavyMesonCalo("00000113","1111111047032230000","302010708","0103603100000000","0153503000000000"); // 100-145
    cuts.AddCutHeavyMesonCalo("00000113","1111111047032230000","302010708","0103603200000000","0153503000000000"); // 110-145
    cuts.AddCutHeavyMesonCalo("00000113","1111111047032230000","302010708","0103603300000000","0153503000000000"); // 120-145
    cuts.AddCutHeavyMesonCalo("00000113","1111111047032230000","302010708","0103603400000000","0153503000000000"); // 100-150
    cuts.AddCutHeavyMesonCalo("00000113","1111111047032230000","302010708","0103603500000000","0153503000000000"); // 110-150
    cuts.AddCutHeavyMesonCalo("00000113","1111111047032230000","302010708","0103603600000000","0153503000000000"); // 120-150
    cuts.AddCutHeavyMesonCalo("00000113","1111111047032230000","302010708","0103603a00000000","0153503000000000"); // 80-145

    // ---------------------------------
    // systematic studies 7 TeV (PHOS)
    // ---------------------------------

  } else if(trainConfig == 150)  { // Standard PHOS
    cuts.AddCutHeavyMesonCalo("00000113","2444411043012300000","302010708","0103603n00000000","0153503000000000"); // PCM-PHOS nonLin (standard)
    cuts.AddCutHeavyMesonCalo("00000113","2444411043012300000","30c010708","0103603n00000000","0153503000000000"); //  with TPC refit
    cuts.AddCutHeavyMesonCalo("00000113","2444411043012300000","322010708","0103603n00000000","0153503000000000"); //  with ITS requirement
    cuts.AddCutHeavyMesonCalo("00000113","2444411043012300000","32c010708","0103603n00000000","0153503000000000"); //  with TPC refit + ITS requirement

    // *************Variations in AliConvEventCuts**************************
  } else if(trainConfig == 151)  { // removePileUp
    cuts.AddCutHeavyMesonCalo("00000013","2444411043012300000","302010708","0103603n00000000","0153503000000000"); // PC

    // *************Variations in AliCaloPhotonsCut**************************
  } else if(trainConfig == 152)  { // Timing diff(std is -100ns to 100ns)
    cuts.AddCutHeavyMesonCalo("00000113","2444411033012300000","302010708","0103603n00000000","0153503000000000"); // 200ns
    cuts.AddCutHeavyMesonCalo("00000113","2444411053012300000","302010708","0103603n00000000","0153503000000000"); // 50ns
    cuts.AddCutHeavyMesonCalo("00000113","2444411073012300000","302010708","0103603n00000000","0153503000000000"); // 30ns
    cuts.AddCutHeavyMesonCalo("00000113","2444411093012300000","302010708","0103603n00000000","0153503000000000"); // 20-25ns
  } else if(trainConfig == 153)  { // Track matching
    cuts.AddCutHeavyMesonCalo("00000113","2444411041012300000","302010708","0103603n00000000","0153503000000000"); //
    cuts.AddCutHeavyMesonCalo("00000113","2444411042012300000","302010708","0103603n00000000","0153503000000000"); //
    cuts.AddCutHeavyMesonCalo("00000113","2444411044012300000","302010708","0103603n00000000","0153503000000000"); //
    cuts.AddCutHeavyMesonCalo("00000113","2444411046012300000","302010708","0103603n00000000","0153503000000000"); //
    cuts.AddCutHeavyMesonCalo("00000113","2444411047012300000","302010708","0103603n00000000","0153503000000000"); //
    cuts.AddCutHeavyMesonCalo("00000113","2444411048012300000","302010708","0103603n00000000","0153503000000000"); //
  } else if(trainConfig == 154)  { // MinEnergy (of cluster) (std is 0.5 GeV)
    cuts.AddCutHeavyMesonCalo("00000113","2444411043022300000","302010708","0103603n00000000","0153503000000000"); // 0.6
    cuts.AddCutHeavyMesonCalo("00000113","2444411043002300000","302010708","0103603n00000000","0153503000000000"); // 0.1
    cuts.AddCutHeavyMesonCalo("00000113","2444411043032300000","302010708","0103603n00000000","0153503000000000"); // 0.7
    cuts.AddCutHeavyMesonCalo("00000113","2444411043042300000","302010708","0103603n00000000","0153503000000000"); // 0.8
  } else if(trainConfig == 155)  { // Min N of cells (std is 2)
    cuts.AddCutHeavyMesonCalo("00000113","2444411043011300000","302010708","0103603n00000000","0153503000000000"); // 1
    cuts.AddCutHeavyMesonCalo("00000113","2444411043031300000","302010708","0103603n00000000","0153503000000000"); // 3
    cuts.AddCutHeavyMesonCalo("00000113","2444411043041300000","302010708","0103603n00000000","0153503000000000"); // 4
  } else if(trainConfig == 156)  { // MinMaxM02 (std is >0.2)
    cuts.AddCutHeavyMesonCalo("00000113","2444411043012200000","302010708","0103603n00000000","0153503000000000"); // >0.1
    cuts.AddCutHeavyMesonCalo("00000113","2444411043012100000","302010708","0103603n00000000","0153503000000000"); // >0.002
    cuts.AddCutHeavyMesonCalo("00000113","2444411043012330000","302010708","0103603n00000000","0153503000000000"); // 0.2-0.5

    // *************Variations in AliPrimaryPionCuts******************
  } else if( trainConfig == 157)  { // ClsTPCCut
    cuts.AddCutHeavyMesonCalo("00000113","2444411043012300000","301010708","0103603n00000000","0153503000000000"); // 70
    cuts.AddCutHeavyMesonCalo("00000113","2444411043012300000","303010708","0103603n00000000","0153503000000000"); // 100
    cuts.AddCutHeavyMesonCalo("00000113","2444411043012300000","305010708","0103603n00000000","0153503000000000"); // 35% of findable clusters
    cuts.AddCutHeavyMesonCalo("00000113","2444411043012300000","306010708","0103603n00000000","0153503000000000"); // 60% of findable clusters
    cuts.AddCutHeavyMesonCalo("00000113","2444411043012300000","30b010708","0103603n00000000","0153503000000000"); // PHOS public note

  } else if ( trainConfig == 158) { // DCACut
    cuts.AddCutHeavyMesonCalo("00000113","2444411043012300000","302110708","0103603n00000000","0153503000000000"); // XYPtDep("0.0182+0.0350/pt^1.01");
    cuts.AddCutHeavyMesonCalo("00000113","2444411043012300000","302210708","0103603n00000000","0153503000000000"); // z=2cm xy=1cm
    cuts.AddCutHeavyMesonCalo("00000113","2444411043012300000","302310708","0103603n00000000","0153503000000000"); // z=3cm XYPtDep("0.0182+0.0350/pt^1.01");
    cuts.AddCutHeavyMesonCalo("00000113","2444411043012300000","302410708","0103603n00000000","0153503000000000"); // z=3cm xy=0.5
  } else if ( trainConfig == 159) { // pT cut
    cuts.AddCutHeavyMesonCalo("00000113","2444411043012300000","302000708","0103603n00000000","0153503000000000"); // pt>0.075
    cuts.AddCutHeavyMesonCalo("00000113","2444411043012300000","302020708","0103603n00000000","0153503000000000"); // pt>0.125
    cuts.AddCutHeavyMesonCalo("00000113","2444411043012300000","302030708","0103603n00000000","0153503000000000"); // pt>0.15
    cuts.AddCutHeavyMesonCalo("00000113","2444411043012300000","302040708","0103603n00000000","0153503000000000"); // pt>0.4
  } else if ( trainConfig == 160) { // TPDdEdxCutPion
    cuts.AddCutHeavyMesonCalo("00000113","2444411043012300000","302010508","0103603n00000000","0153503000000000"); // -4,4
    cuts.AddCutHeavyMesonCalo("00000113","2444411043012300000","302010808","0103603n00000000","0153503000000000"); // -2,3
    cuts.AddCutHeavyMesonCalo("00000113","2444411043012300000","302010208","0103603n00000000","0153503000000000"); // -6,7
    cuts.AddCutHeavyMesonCalo("00000113","2444411043012300000","302010308","0103603n00000000","0153503000000000"); // -5,5
    cuts.AddCutHeavyMesonCalo("00000113","2444411043012300000","302010408","0103603n00000000","0153503000000000"); // -4,5
  } else if ( trainConfig == 161) { // Mass cut
    cuts.AddCutHeavyMesonCalo("00000113","2444411043012300000","302010707","0103603n00000000","0153503000000000"); // 0.7 GeV
    cuts.AddCutHeavyMesonCalo("00000113","2444411043012300000","302010706","0103603n00000000","0153503000000000"); // 0.65 GeV
    cuts.AddCutHeavyMesonCalo("00000113","2444411043012300000","302010701","0103603n00000000","0153503000000000"); // 1 GeV
    cuts.AddCutHeavyMesonCalo("00000113","2444411043012300000","302010702","0103603n00000000","0153503000000000"); // 0.75 GeV
    cuts.AddCutHeavyMesonCalo("00000113","2444411043012300000","302010704","0103603n00000000","0153503000000000"); // 0.54 eta mass
    // *************Variations in AliConversionMesonCuts (NeutralPion) ******************
  } else if ( trainConfig == 162) { // RapidityMesonCut
    cuts.AddCutHeavyMesonCalo("00000113","2444411043012300000","302010708","0103503n00000000","0153503000000000"); // 0.85
    cuts.AddCutHeavyMesonCalo("00000113","2444411043012300000","302010708","0103303n00000000","0153503000000000"); // 0.6
    cuts.AddCutHeavyMesonCalo("00000113","2444411043012300000","302010708","0103203n00000000","0153503000000000"); // 0.7
  } else if ( trainConfig == 163) { // pT
    cuts.AddCutHeavyMesonCalo("00000113","2444411043012300000","302010708","0103613n00000000","0153503000000000"); // 0.4
    cuts.AddCutHeavyMesonCalo("00000113","2444411043012300000","302010708","0103623n00000000","0153503000000000"); // 0.7
    cuts.AddCutHeavyMesonCalo("00000113","2444411043012300000","302010708","0103673n00000000","0153503000000000"); // 0.5
  } else if ( trainConfig == 164) { // alphaMesonCut
    cuts.AddCutHeavyMesonCalo("00000113","2444411043012300000","302010708","0103607n00000000","0153503000000000"); // 0-0.85
    cuts.AddCutHeavyMesonCalo("00000113","2444411043012300000","302010708","0103605n00000000","0153503000000000"); // 0-0.75
    cuts.AddCutHeavyMesonCalo("00000113","2444411043012300000","302010708","0103600n00000000","0153503000000000"); // 0-0.7
    cuts.AddCutHeavyMesonCalo("00000113","2444411043012300000","302010708","0103604n00000000","0153503000000000"); // 0-0.65
  } else if ( trainConfig == 165) { // selectionWindow (std is 120-160)
    cuts.AddCutHeavyMesonCalo("00000113","2444411043012300000","302010708","0103603100000000","0153503000000000"); // 0.1-0.145
    cuts.AddCutHeavyMesonCalo("00000113","2444411043012300000","302010708","0103603200000000","0153503000000000"); // 0.11-0.145
    cuts.AddCutHeavyMesonCalo("00000113","2444411043012300000","302010708","0103603300000000","0153503000000000"); // 0.12-0.145
    cuts.AddCutHeavyMesonCalo("00000113","2444411043012300000","302010708","0103603600000000","0153503000000000"); // 0.12-0.5
    cuts.AddCutHeavyMesonCalo("00000113","2444411043012300000","302010708","0103603700000000","0153503000000000"); // 0.1 -0.155
    cuts.AddCutHeavyMesonCalo("00000113","2444411043012300000","302010708","0103603900000000","0153503000000000"); // 0.11 -0.155
    cuts.AddCutHeavyMesonCalo("00000113","2444411043012300000","302010708","0103603a00000000","0153503000000000"); // 0.08 -0.145
    // *************Variations in AliConversionMesonCuts (omega) ******************
  } else if ( trainConfig == 166) { // selectionWindow (std is 120-160)
    cuts.AddCutHeavyMesonCalo("00000113","2444411043012300000","302010708","0103603n00000000","0a53503000000000"); // likesign
    cuts.AddCutHeavyMesonCalo("00000113","2444411043012300000","302010708","0103603n00000000","0b53503000000000"); // sideband right
    cuts.AddCutHeavyMesonCalo("00000113","2444411043012300000","302010708","0103603n00000000","0c53503000000000"); // sideband left
    cuts.AddCutHeavyMesonCalo("00000113","2444411043012300000","302010708","0103603n00000000","0d53503000000000"); // sideband both sides
  } else if ( trainConfig == 167) { // Number of BckEvents
    cuts.AddCutHeavyMesonCalo("00000113","2444411043012300000","302010708","0103603n00000000","0133503000000000"); // 20
    cuts.AddCutHeavyMesonCalo("00000113","2444411043012300000","302010708","0103603n00000000","0163503000000000"); // 80
  } else if ( trainConfig == 168) { // rapidity cut
    cuts.AddCutHeavyMesonCalo("00000113","2444411043012300000","302010708","0103603n00000000","0153203000000000"); // 0.7
    cuts.AddCutHeavyMesonCalo("00000113","2444411043012300000","302010708","0103603n00000000","0153003000000000"); // 1.35
    cuts.AddCutHeavyMesonCalo("00000113","2444411043012300000","302010708","0103603n00000000","0153303000000000"); // 0.6
  } else if ( trainConfig == 169) { // alphacut
    cuts.AddCutHeavyMesonCalo("00000113","2444411043012300000","302010708","0103603n00000000","0153507000000000"); // 0-0.85
    cuts.AddCutHeavyMesonCalo("00000113","2444411043012300000","302010708","0103603n00000000","0153505000000000"); // 0-0.75
    cuts.AddCutHeavyMesonCalo("00000113","2444411043012300000","302010708","0103603n00000000","0153500000000000"); // 0-0.7
    // PHOS pp 5 TeV
  } else if(trainConfig == 190)  { // Standard PHOS  with TPC refit + ITS requirement
    cuts.AddCutHeavyMesonCalo("00010113","2444411043012300000","32c010708","0103603n00000000","0153503000000000"); // INT7
    cuts.AddCutHeavyMesonCalo("00062113","2444411043012300000","32c010708","0103603n00000000","0153503000000000"); // PHI7
    // PHOS pp 13 TeV
  } else if(trainConfig == 191)  { // Standard PHOS  with TPC refit + ITS requirement
    cuts.AddCutHeavyMesonCalo("00010113","2444411043012300000","32c010708","0103603n00000000","0153503000000000"); // INT7
    cuts.AddCutHeavyMesonCalo("00062113","2444411043012300000","32c010708","0103603n00000000","0153503000000000"); // PHI7
    // PHOS LHC11 pp 7 TeV
  } else if(trainConfig == 195)  { // with TPC refit + ITS requirement
    cuts.AddCutHeavyMesonCalo("00010113","2444400053012300000","32c010708","0103603n00000000","0153503000000000"); // INT7
    cuts.AddCutHeavyMesonCalo("00062113","2444400053012300000","32c010708","0103603n00000000","0153503000000000"); // PHI7
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
    if(runLightOutput>0) analysisEventCuts[i]->SetLightOutput(kTRUE);
    analysisEventCuts[i]->InitializeCutsFromCutString((cuts.GetEventCut(i)).Data());
    if (periodNameV0Reader.CompareTo("") != 0) analysisEventCuts[i]->SetPeriodEnum(periodNameV0Reader);
    EventCutList->Add(analysisEventCuts[i]);
    analysisEventCuts[i]->SetFillCutHistograms("",kFALSE);

    analysisClusterCuts[i] = new AliCaloPhotonCuts();
    analysisClusterCuts[i]->SetV0ReaderName(V0ReaderName);
    if(runLightOutput>0) analysisClusterCuts[i]->SetLightOutput(kTRUE);
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
