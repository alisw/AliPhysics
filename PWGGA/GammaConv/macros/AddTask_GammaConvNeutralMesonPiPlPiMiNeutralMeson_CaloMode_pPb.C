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
//pPb together with all supporting classes
//***************************************************************************************

//***************************************************************************************
//main function
//***************************************************************************************
void AddTask_GammaConvNeutralMesonPiPlPiMiNeutralMeson_CaloMode_pPb(
    Int_t     trainConfig                 = 1,
    Int_t     isMC                        = 0,                        //run MC
    TString   photonCutNumberV0Reader     = "",                       // 00000008400000000100000000 nom. B, 00000088400000000100000000 low B
    Int_t     selectHeavyNeutralMeson     = 0,                        //run eta prime instead of omega
    Int_t     enableQAMesonTask           = 1,                        //enable QA in AliAnalysisTaskNeutralMesonToPiPlPiMiNeutralMeson
    Int_t     enableExtMatchAndQA         = 0,                        // disabled (0), extMatch (1), extQA_noCellQA (2), extMatch+extQA_noCellQA (3), extQA+cellQA (4), extMatch+extQA+cellQA (5)
    Int_t     enableTriggerMimicking      = 0,                        // enable trigger mimicking
    Bool_t    enableTriggerOverlapRej     = kFALSE,                   // enable trigger overlap rejection
    TString   fileNameInputForWeighting   = "MCSpectraInput.root",    // path to file for weigting input
    Int_t     doWeightingPart             = 0,                        //enable Weighting
    TString   generatorName               = "HIJING",
    Double_t  tolerance                   = -1,
    TString   periodNameV0Reader          = "",                       // period Name for V0Reader
    Int_t     runLightOutput              = 0,                        // run light output option 0: no light output 1: most cut histos stiched off 2: unecessary omega hists turned off as well
    Int_t     prefilterRunFlag            = 1500,                     // flag to change the prefiltering of ESD tracks. See SetHybridTrackCutsAODFiltering() in AliPrimaryPionCuts
    TString   additionalTrainConfig       = "0"                       // additional counter for trainconfig, this has to be always the last parameter
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
        cout << Form("INFO: AddTask_GammaConvNeutralMesonPiPlPiMiPiZero_CaloMode_pPb will use running mode '%i' for the TrackMatcher!",trackMatcherRunningMode) << endl;
      }
    }
  }
  TString sAdditionalTrainConfig = rAdditionalTrainConfig->GetString();
  if (sAdditionalTrainConfig.Atoi() > 0){
    trainConfig = trainConfig + sAdditionalTrainConfig.Atoi();
    std::cout << "INFO: AddTask_GammaConvNeutralMesonPiPlPiMiNeutralMeson_CaloMode_pPb running additionalTrainConfig '" << sAdditionalTrainConfig.Atoi() << "', train config: '" << trainConfig << "'" << std::endl;
  }

  Int_t isHeavyIon = 2;
  Int_t neutralPionMode = 2;

  // ================== GetAnalysisManager ===============================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error(Form("AddTask_GammaConvNeutralMesonPiPlPiMiNeutralMeson_CaloMode_pPb_%i",trainConfig), "No analysis manager found.");
    return ;
  }

  // ================== GetInputEventHandler =============================
  AliVEventHandler *inputHandler=mgr->GetInputEventHandler();


  //=========  Set Cutnumber for V0Reader ================================
  TString cutnumberPhoton = photonCutNumberV0Reader.Data();
  TString cutnumberEvent = "80000003";
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


  AliCutHandlerPCM cuts(13);


  // ******************************************************************************************************
  // ++++++++++++++++++++++++++++++   N A M I N G  C O N V E N T I O N   ++++++++++++++++++++++++++++++++++
  //
  //  config = meson * 1000  + method * 100 + variation
  //
  //  meson:  eta = 0; omega = 1; etaPrime = 2
  //  method: EMC = 0; PHOS = 5;
  //
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // ******************************************************************************************************


  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //                                          ETA MESON
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   //************************************************ PCM- EDC analysis 5 TeV pPb *********************************************
  // no event mixing background
  if (trainConfig == 1){ // EMC  INT7 run1 & run2
    cuts.AddCutHeavyMesonCalo("80010113","411790105f032230000","32c51070a","0103603o00000000","0453503000000000"); // 0-100% without NL
    cuts.AddCutHeavyMesonCalo("80010113","111110105f032230000","32c51070a","0103603o00000000","0453503000000000"); // 0-100% without NL just EMC
  } else if (trainConfig == 2){ // EMC  INT7 run1 & run2
    cuts.AddCutHeavyMesonCalo("80010113","411793105f032230000","32c51070a","0103603o00000000","0453503000000000"); // 0-100% PCM NL
    cuts.AddCutHeavyMesonCalo("80010113","111113105f032230000","32c51070a","0103603o00000000","0453503000000000"); // 0-100% PCM NL just EMC
  } else if (trainConfig == 3){ // EMC EMC triggers
    cuts.AddCutHeavyMesonCalo("80052113","111110105f032230000","32c51070a","0103603o00000000","0453503000000000"); // 0-100% without NL just EMC, EMC7
    cuts.AddCutHeavyMesonCalo("80052113","111113105f032230000","32c51070a","0103603o00000000","0453503000000000"); // 0-100% PCM NL just EMC, EMC7
  } else if (trainConfig == 4){ // EMC EMC triggers
    cuts.AddCutHeavyMesonCalo("80083113","111110105f032230000","32c51070a","0103603o00000000","0453503000000000"); // 0-100% without NL just EMC, EG1
    cuts.AddCutHeavyMesonCalo("80083113","111113105f032230000","32c51070a","0103603o00000000","0453503000000000"); // 0-100% PCM NL just EMC, EG1
  } else if (trainConfig == 5){ // EMC EMC triggers
    cuts.AddCutHeavyMesonCalo("80085113","111110105f032230000","32c51070a","0103603o00000000","0453503000000000"); // 0-100% without NL just EMC, EG2
    cuts.AddCutHeavyMesonCalo("80085113","111113105f032230000","32c51070a","0103603o00000000","0453503000000000"); // 0-100% PCM NL just EMC, EG2

    // EMCal + EDC triggers
  } else if (trainConfig == 6){ // EMC EMC triggers
    cuts.AddCutHeavyMesonCalo("8008d113","411790105f032230000","32c51070a","0103603o00000000","0453503000000000"); // 0-100% without NL just EMC, EG1
    cuts.AddCutHeavyMesonCalo("8008d113","411793105f032230000","32c51070a","0103603o00000000","0453503000000000"); // 0-100% PCM NL just EMC, EG1
  } else if (trainConfig == 7){ // EMC EMC triggers
    cuts.AddCutHeavyMesonCalo("8008e113","411790105f032230000","32c51070a","0103603o00000000","0453503000000000"); // 0-100% without NL just EMC, EG2
    cuts.AddCutHeavyMesonCalo("8008e113","411793105f032230000","32c51070a","0103603o00000000","0453503000000000"); // 0-100% PCM NL just EMC, EG2

 //************************************************ PCM- PHOS analysis 5 TeV pPb ********************************************
  } else if (trainConfig == 501){ // PHOS  INT7 run1
    cuts.AddCutHeavyMesonCalo("80010113","244440004a013200000","32c51070a","0103603q00000000","0453503000000000");  // 0-100% without NL
  } else if (trainConfig == 502){ // PHOS  PHI7 run1
    cuts.AddCutHeavyMesonCalo("80062113","244440004a013200000","32c51070a","0103603q00000000","0453503000000000");  // 0-100% without NL
  } else if (trainConfig == 503) {  // PHOS  INT7 run2
    cuts.AddCutHeavyMesonCalo("80010113","24466000ha012200000","32c51070a","0103603q00000000","0453503000000000");  // 0-100% without NL

  } else if (trainConfig == 504){ // PHOS  INT7 run1
    cuts.AddCutHeavyMesonCalo("80010113","244445104a013200000","32c51070a","0103603q00000000","0453503000000000");  // 0-100% PCM NL
  } else if (trainConfig == 505){ // PHOS  PHI7 run1
    cuts.AddCutHeavyMesonCalo("80062113","244445104a013200000","32c51070a","0103603q00000000","0453503000000000");  // 0-100% PCM NL
  } else if (trainConfig == 506) {  // PHOS  INT7 run2
    cuts.AddCutHeavyMesonCalo("80010113","24466410ha012200000","32c51070a","0103603q00000000","0453503000000000");  // 0-100% PCM NL

    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //                                          OMEGA MESON
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  //************************************************ PCM- EDC analysis 5 TeV pPb *********************************************

  // no event mixing background
  }else if (trainConfig == 1001){ // EMC  INT7 run1 & run2
    cuts.AddCutHeavyMesonCalo("80010113","411790105f032230000","32c51070a","0103603o00000000","0453503000000000"); // 0-100% without NL
    cuts.AddCutHeavyMesonCalo("80010113","111110105f032230000","32c51070a","0103603o00000000","0453503000000000"); // 0-100% without NL just EMC
  } else if (trainConfig == 1002){ // EMC  INT7 run1 & run2
    cuts.AddCutHeavyMesonCalo("80010113","411793105f032230000","32c51070a","0103603o00000000","0453503000000000"); // 0-100% PCM NL
    cuts.AddCutHeavyMesonCalo("80010113","111113105f032230000","32c51070a","0103603o00000000","0453503000000000"); // 0-100% PCM NL just EMC
  } else if (trainConfig == 1003){ // EMC EMC triggers
    cuts.AddCutHeavyMesonCalo("80052113","111110105f032230000","32c51070a","0103603o00000000","0453503000000000"); // 0-100% without NL just EMC, EMC7
    cuts.AddCutHeavyMesonCalo("80052113","111113105f032230000","32c51070a","0103603o00000000","0453503000000000"); // 0-100% PCM NL just EMC, EMC7
  } else if (trainConfig == 1004){ // EMC EMC triggers
    cuts.AddCutHeavyMesonCalo("80083113","111110105f032230000","32c51070a","0103603o00000000","0453503000000000"); // 0-100% without NL just EMC, EG1
    cuts.AddCutHeavyMesonCalo("80083113","111113105f032230000","32c51070a","0103603o00000000","0453503000000000"); // 0-100% PCM NL just EMC, EG1
  } else if (trainConfig == 1005){ // EMC EMC triggers
    cuts.AddCutHeavyMesonCalo("80085113","111110105f032230000","32c51070a","0103603o00000000","0453503000000000"); // 0-100% without NL just EMC, EG2
    cuts.AddCutHeavyMesonCalo("80085113","111113105f032230000","32c51070a","0103603o00000000","0453503000000000"); // 0-100% PCM NL just EMC, EG2

    // EMCal + EDC triggers
  } else if (trainConfig == 1006){ // EMC EMC triggers
    cuts.AddCutHeavyMesonCalo("8008d113","411790105f032230000","32c51070a","0103603o00000000","0453503000000000"); // 0-100% without NL just EMC, EG1
    cuts.AddCutHeavyMesonCalo("8008d113","411793105f032230000","32c51070a","0103603o00000000","0453503000000000"); // 0-100% PCM NL just EMC, EG1
  } else if (trainConfig == 1007){ // EMC EMC triggers
    cuts.AddCutHeavyMesonCalo("8008e113","411790105f032230000","32c51070a","0103603o00000000","0453503000000000"); // 0-100% without NL just EMC, EG2
    cuts.AddCutHeavyMesonCalo("8008e113","411793105f032230000","32c51070a","0103603o00000000","0453503000000000"); // 0-100% PCM NL just EMC, EG2

  }else if (trainConfig == 1008){ // EMC  INT7 mixed background study pi+ pi- from same
    cuts.AddCutHeavyMesonCalo("80010113","411790105f032230000","32c51070a","0103603o00000000","0o53503000000000"); // Mixed event pi+ pi- from same

  }else if (trainConfig == 1009){ // EMC  INT7 mixed background study pi+ pi0 from same
    cuts.AddCutHeavyMesonCalo("80010113","411790105f032230000","32c51070a","0103603o00000000","0p53503000000000"); // Mixed event pi+ pi0 from same

  }else if (trainConfig == 1010){ // EMC  INT7 background study
    cuts.AddCutHeavyMesonCalo("80010113","411790105f032230000","32c51070a","0103603o00000000","0r53503000000000"); // Rotation around pi0
    cuts.AddCutHeavyMesonCalo("80010113","411790105f032230000","32c51070a","0103603o00000000","0a53503000000000"); // Likesign method

 //************************************************ PCM- PHOS analysis 5 TeV pPb ********************************************
  } else if (trainConfig == 1501){ // PHOS  INT7 run1
    cuts.AddCutHeavyMesonCalo("80010113","244440004a013200000","32c51070a","0103603q00000000","0453503000000000");  // 0-100% without NL
  } else if (trainConfig == 1502){ // PHOS  PHI7 run1
    cuts.AddCutHeavyMesonCalo("80062113","244440004a013200000","32c51070a","0103603q00000000","0453503000000000");  // 0-100% without NL
  } else if (trainConfig == 1503) {  // PHOS  INT7 run2
    cuts.AddCutHeavyMesonCalo("80010113","24466000ha012200000","32c51070a","0103603q00000000","0453503000000000");  // 0-100% without NL

  } else if (trainConfig == 1504){ // PHOS  INT7 run1
    cuts.AddCutHeavyMesonCalo("80010113","244445104a013200000","32c51070a","0103603q00000000","0453503000000000");  // 0-100% PCM NL
  } else if (trainConfig == 1505){ // PHOS  PHI7 run1
    cuts.AddCutHeavyMesonCalo("80062113","244445104a013200000","32c51070a","0103603q00000000","0453503000000000");  // 0-100% PCM NL
  } else if (trainConfig == 1506) {  // PHOS  INT7 run2
    cuts.AddCutHeavyMesonCalo("80010113","24466410ha012200000","32c51070a","0103603q00000000","0453503000000000");  // 0-100% PCM NL

  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //                                          ETA PRIME MESON
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  //************************************************ PCM- EDC analysis 5 TeV pPb *********************************************

  // no event mixing background
  } else if (trainConfig == 2001){ // EMC  INT7 run1 & run2
    cuts.AddCutHeavyMesonCalo("80010113","411790105f032230000","32c510700","0103603l00000000","0453503000000000"); // 0-100% without NL
    cuts.AddCutHeavyMesonCalo("80010113","111110105f032230000","32c510700","0103603l00000000","0453503000000000"); // 0-100% without NL just EMC
  } else if (trainConfig == 2002){ // EMC  INT7 run1 & run2
    cuts.AddCutHeavyMesonCalo("80010113","411793105f032230000","32c510700","0103603l00000000","0453503000000000"); // 0-100% PCM NL
    cuts.AddCutHeavyMesonCalo("80010113","111113105f032230000","32c510700","0103603l00000000","0453503000000000"); // 0-100% PCM NL just EMC
  } else if (trainConfig == 2003){ // EMC EMC triggers
    cuts.AddCutHeavyMesonCalo("80052113","111110105f032230000","32c510700","0103603l00000000","0453503000000000"); // 0-100% without NL just EMC, EMC7
    cuts.AddCutHeavyMesonCalo("80052113","111113105f032230000","32c510700","0103603l00000000","0453503000000000"); // 0-100% PCM NL just EMC, EMC7
  } else if (trainConfig == 2004){ // EMC EMC triggers
    cuts.AddCutHeavyMesonCalo("80083113","111110105f032230000","32c510700","0103603l00000000","0453503000000000"); // 0-100% without NL just EMC, EG1
    cuts.AddCutHeavyMesonCalo("80083113","111113105f032230000","32c510700","0103603l00000000","0453503000000000"); // 0-100% PCM NL just EMC, EG1
  } else if (trainConfig == 2005){ // EMC EMC triggers
    cuts.AddCutHeavyMesonCalo("80085113","111110105f032230000","32c510700","0103603l00000000","0453503000000000"); // 0-100% without NL just EMC, EG2
    cuts.AddCutHeavyMesonCalo("80085113","111113105f032230000","32c510700","0103603l00000000","0453503000000000"); // 0-100% PCM NL just EMC, EG2

    // EMCal + EDC triggers
  } else if (trainConfig == 2006){ // EMC EMC triggers
    cuts.AddCutHeavyMesonCalo("8008d113","411790105f032230000","32c510700","0103603l00000000","0453503000000000"); // 0-100% without NL just EMC, EG1
    cuts.AddCutHeavyMesonCalo("8008d113","411793105f032230000","32c510700","0103603l00000000","0453503000000000"); // 0-100% PCM NL just EMC, EG1
  } else if (trainConfig == 2007){ // EMC EMC triggers
    cuts.AddCutHeavyMesonCalo("8008e113","411790105f032230000","32c510700","0103603l00000000","0453503000000000"); // 0-100% without NL just EMC, EG2
    cuts.AddCutHeavyMesonCalo("8008e113","411793105f032230000","32c510700","0103603l00000000","0453503000000000"); // 0-100% PCM NL just EMC, EG2

  //************************************************ PCM- PHOS analysis 5 TeV pPb ********************************************
  } else if (trainConfig == 2501){ // PHOS  INT7 run1
    cuts.AddCutHeavyMesonCalo("80010113","244440004a013200000","32c510700","0103603l00000000","0453503000000000"); // 0-100% without NL
  } else if (trainConfig == 2502){ // PHOS  PHI7 run1
    cuts.AddCutHeavyMesonCalo("80062113","244440004a013200000","32c510700","0103603l00000000","0453503000000000"); // 0-100% without NL
  } else if (trainConfig == 2503) {  // PHOS  INT7 run2
    cuts.AddCutHeavyMesonCalo("80010113","24466000ha012200000","32c510700","0103603l00000000","0453503000000000"); // 0-100% without NL

  } else if (trainConfig == 2504){ // PHOS  INT7 run1
    cuts.AddCutHeavyMesonCalo("80010113","244445104a013200000","32c510700","0103603l00000000","0453503000000000"); // 0-100% PCM NL
  } else if (trainConfig == 2505){ // PHOS  PHI7 run1
    cuts.AddCutHeavyMesonCalo("80062113","244445104a013200000","32c510700","0103603l00000000","0453503000000000"); // 0-100% PCM NL
  } else if (trainConfig == 2506) {  // PHOS  INT7 run2
    cuts.AddCutHeavyMesonCalo("80010113","24466410ha012200000","32c510700","0103603l00000000","0453503000000000"); // 0-100% PCM NL
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

  if (periodNameV0Reader.Contains("LHC17g6a2") || periodNameV0Reader.Contains("LHC17g6a3") ){
    TObjString *HeaderPMB = new TObjString("Dpmjet_0");
    TObjString *HeaderP8J = new TObjString("Pythia8JetsGammaTrg_1");
    if (doWeightingPart==4) { // all headers
      HeaderList->Add(HeaderPMB);
      HeaderList->Add(HeaderP8J);
    } else if (doWeightingPart==5) { // only MB header
      HeaderList->Add(HeaderPMB);
    } else { // only JJ header
      HeaderList->Add(HeaderP8J);
    }
  }

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
    analysisEventCuts[i]->SetTriggerOverlapRejecion(enableTriggerOverlapRej);
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
