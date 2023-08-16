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
    Bool_t    usePtDepSelectionWindowCut  = kFALSE,                   // use pt dependent meson selection window cut
    Bool_t    enableSortingMCLabels       = kTRUE,                    // enable sorting for MC cluster labels
    TString   additionalTrainConfig       = "0"                       // additional counter for trainconfig, this has to be always the last parameter
  ) {

  AliCutHandlerPCM cuts(13);
  TString addTaskName                       = "AddTask_GammaConvNeutralMesonPiPlPiMiNeutralMeson_CaloMode_pPb";
  
  //parse additionalTrainConfig flag
  Int_t trackMatcherRunningMode = 0; // CaloTrackMatcher running mode
  TString unsmearingoutputs = "012"; // 0: No correction, 1: One pi0 mass errer subtracted, 2: pz of pi0 corrected to fix its mass, 3: Lambda(alpha)*DeltaPi0 subtracted

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
        cout << Form("INFO: AddTask_GammaConvNeutralMesonPiPlPiMiPiZero_CaloMode_pPb will use running mode '%i' for the TrackMatcher!",trackMatcherRunningMode) << endl;
      }
      if(tempStr.BeginsWith("UNSMEARING")){
        TString tempType = tempStr;
        tempType.Replace(0,9,"");
        unsmearingoutputs = tempType;
        cout << "INFO: AddTask_GammaConvNeutralMesonPiPlPiMiPiZero_CaloMode_pPb will output the following minv_pT histograms:" << endl;
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
  task->SetCorrectionTaskSetting(corrTaskSetting);


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

  }else if (trainConfig == 1011){ // EMC  INT7 standard cut study (1)
    cuts.AddCutHeavyMesonCalo("80010113","411790105fe32220000","32c51070a","0103603o00000000","0453503000000000"); // M02 max 0.7 (new standard)
  }else if (trainConfig == 1012){ // EMC  INT7 cut study 2
    cuts.AddCutHeavyMesonCalo("80010113","411790104fe32220000","32c51070a","0103603o00000000","0453503000000000"); // timing +-100
    cuts.AddCutHeavyMesonCalo("80010113","411790106fe32220000","32c51070a","0103603o00000000","0453503000000000"); // timing +-35
  }else if (trainConfig == 1013){ // EMC  INT7 cut study 3
    cuts.AddCutHeavyMesonCalo("80010113","4117901051e32220000","32c51070a","0103603o00000000","0453503000000000"); // fixed detadphi cut (not pT dep.)
    cuts.AddCutHeavyMesonCalo("80010113","4117901057e32220000","32c51070a","0103603o00000000","0453503000000000"); // pt dep. detadphi cut, no E/p cut
    cuts.AddCutHeavyMesonCalo("80010113","411790105fe32220000","32c51070a","0103603700000000","0453503000000000"); // large pi0 selection window
  }else if (trainConfig == 1014){ // EMC  INT7 cut study 4
    cuts.AddCutHeavyMesonCalo("80010113","411793305fe32220000","32c51070a","0103603o00000000","0453503000000000"); // NonLinearity 33
    cuts.AddCutHeavyMesonCalo("80010113","411793405fe32220000","32c51070a","0103603o00000000","0453503000000000"); // NonLinearity 34
  }else if (trainConfig == 1015){ // EMC  INT7 cut study 5
    cuts.AddCutHeavyMesonCalo("80010113","411790105fe32220000","32c51070a","0103603o00000000","0453603000000000"); // y < 0.75
    cuts.AddCutHeavyMesonCalo("80010113","411790105f032220000","32c51070a","0103603o00000000","0453503000000000"); // no exotic cut
  }else if (trainConfig == 1016){ // EMC  INT7 cut study 6
    cuts.AddCutHeavyMesonCalo("80010113","411790105fe32240000","32c51070a","0103603o00000000","0453503000000000"); // M02 max 0.4
    cuts.AddCutHeavyMesonCalo("80010113","411790105fe32240000","32c51070a","0103603o00000000","0453503000000000"); // M02 max 1
  }else if (trainConfig == 1017){ // EMC  INT7 cut study 7
    cuts.AddCutHeavyMesonCalo("80010113","411790105fe32220000","32c51050a","0103603o00000000","0453503000000000"); // dEdx +-4 sigma
    cuts.AddCutHeavyMesonCalo("80010113","411790105fe32220000","32c51080a","0103603o00000000","0453503000000000"); // dEdx -2+3 sigma

    //***********EMCal Trigger****************
  }else if (trainConfig == 1018){ // EMC MinBias reference for run1 triggered
    cuts.AddCutHeavyMesonCalo("80010113","111110105f030230000","32c51070a","0103603o00000000","0453503000000000"); // Min Bias run1 reference
  }else if (trainConfig == 1019){ // EMC run1 EG1 trigger
    cuts.AddCutHeavyMesonCalo("80083113","111110105f030230000","32c51070a","0103603o00000000","0453503000000000"); // 0-100% without NL just EMC, EG1
  }else if (trainConfig == 1020){ // EMC run1 EG2 trigger
    cuts.AddCutHeavyMesonCalo("80085113","111110105f030230000","32c51070a","0103603o00000000","0453503000000000"); // 0-100% without NL just EMC, EG2

    //***********New cutstrings without NCell cut*******************************
  }else if (trainConfig == 1021){ // EMC  INT7 new standard
    cuts.AddCutHeavyMesonCalo("80010113","411790105fe30220000","32c51070a","0103603o00000000","0453503000000000");
  }else if (trainConfig == 1022){ // Cutvariation for NCell cut
    cuts.AddCutHeavyMesonCalo("80010113","411790105fe3n220000","32c51070a","0103603o00000000","0453503000000000");
  }else if (trainConfig == 1023){ // Cutvariation for pi0 pT
    cuts.AddCutHeavyMesonCalo("80010113","411790105fe30220000","32c51070a","0103663o00000000","0453503000000000"); // pi0 pT > 1.5
    cuts.AddCutHeavyMesonCalo("80010113","411790105fe30220000","32c51070a","01036z3o00000000","0453503000000000"); // pi0 pT > 2
  }else if (trainConfig == 1024){ // Asymmetry cut on omega
    cuts.AddCutHeavyMesonCalo("80010113","411790105fe30220000","32c51070a","0103603o00000000","045350l000000000"); 


  // ********************************************************************
  // ************** Cuts for the omega p-Pb 5TeV analysis ***************
  // ****************** Minimum Bias - Run 2 (LHC16qt) ******************
  }else if (trainConfig == 1100){ // Standard 5 TeV omega INT7 cutstring
    cuts.AddCutHeavyMesonCalo("80010113","411790105fe30220000","32c51079a","0000003100000000","0400503000000000");
  }else if (trainConfig == 1101){ // Non linearity variations
    cuts.AddCutHeavyMesonCalo("80010113","411799705fe30220000","32c51079a","0000003100000000","0400503000000000"); // 97: CRF
    cuts.AddCutHeavyMesonCalo("80010113","411799805fe30220000","32c51079a","0000003100000000","0400503000000000"); // 98: CCRF
  }else if (trainConfig == 1102){ // Cluster timing variations
    cuts.AddCutHeavyMesonCalo("80010113","411790106fe30220000","32c51079a","0000003100000000","0400503000000000"); // 6: -30 - 35
    cuts.AddCutHeavyMesonCalo("80010113","411790107fe30220000","32c51079a","0000003100000000","0400503000000000"); // 7: -30 - 30
    cuts.AddCutHeavyMesonCalo("80010113","411790108fe30220000","32c51079a","0000003100000000","0400503000000000"); // 8: -20 - 30
    cuts.AddCutHeavyMesonCalo("80010113","411790109fe30220000","32c51079a","0000003100000000","0400503000000000"); // 9: -20 - 25
    cuts.AddCutHeavyMesonCalo("80010113","41179010afe30220000","32c51079a","0000003100000000","0400503000000000"); // a: -12.5 - 13
  }else if (trainConfig == 1103){ // Track matching variations
    cuts.AddCutHeavyMesonCalo("80010113","411790105ce30220000","32c51079a","0000003100000000","0400503000000000"); // c: No E/p cut
    cuts.AddCutHeavyMesonCalo("80010113","411790105ee30220000","32c51079a","0000003100000000","0400503000000000"); // e: E/p < 2
    cuts.AddCutHeavyMesonCalo("80010113","411790105ge30220000","32c51079a","0000003100000000","0400503000000000"); // g: E/p < 1.5
    cuts.AddCutHeavyMesonCalo("80010113","411790105ne30220000","32c51079a","0000003100000000","0400503000000000"); // n: E/p < 1.75 (standard), but eta and phi varied
    cuts.AddCutHeavyMesonCalo("80010113","411790105oe30220000","32c51079a","0000003100000000","0400503000000000"); // o: E/p < 1.75 (standard), but eta and phi varied
  }else if (trainConfig == 1104){ // Exotic cluster variations
    cuts.AddCutHeavyMesonCalo("80010113","411790105f030220000","32c51079a","0000003100000000","0400503000000000"); // 0: No exotics cut
    cuts.AddCutHeavyMesonCalo("80010113","411790105fb30220000","32c51079a","0000003100000000","0400503000000000"); // b: fExoticEnergyFracCluster = 0.95
    cuts.AddCutHeavyMesonCalo("80010113","411790105fi30220000","32c51079a","0000003100000000","0400503000000000"); // i: fExoticMinEnergyCell = 3
  }else if (trainConfig == 1105){ // Min cluster energy variations
    cuts.AddCutHeavyMesonCalo("80010113","411790105fe20220000","32c51079a","0000003100000000","0400503000000000"); // 2: E > 0.6
    cuts.AddCutHeavyMesonCalo("80010113","411790105fe40220000","32c51079a","0000003100000000","0400503000000000"); // 4: E > 0.8
  }else if (trainConfig == 1106){ // NCell variation
    cuts.AddCutHeavyMesonCalo("80010113","411790105fe3n220000","32c51079a","0000003100000000","0400503000000000");
  }else if (trainConfig == 1107){ // M02 variations
    cuts.AddCutHeavyMesonCalo("80010113","411790105fe30210000","32c51079a","0000003100000000","0400503000000000"); // 1: M02 < 1
    cuts.AddCutHeavyMesonCalo("80010113","411790105fe30230000","32c51079a","0000003100000000","0400503000000000"); // 3: M02 < 0.5
  }else if (trainConfig == 1108){ // pi0 mass selection window variations
    cuts.AddCutHeavyMesonCalo("80010113","411790105fe30220000","32c51079a","0000003u00000000","0400503000000000"); // u: 1.5 sigma
    cuts.AddCutHeavyMesonCalo("80010113","411790105fe30220000","32c51079a","0000003v00000000","0400503000000000"); // v: 2.5 sigma
    cuts.AddCutHeavyMesonCalo("80010113","411790105fe30220000","32c51079a","0000003x00000000","0400503000000000"); // x: 3 sigma
    cuts.AddCutHeavyMesonCalo("80010113","411790105fe30220000","32c51079a","0000003r00000000","0400503000000000"); // r: 3.5 sigma, no gamma selection
    cuts.AddCutHeavyMesonCalo("80010113","411790105fe30220000","32c51079a","0000003w00000000","0400503000000000"); // w: 4 sigma
  }else if (trainConfig == 1109){ // pi0 asymmetry variations
    cuts.AddCutHeavyMesonCalo("80010113","411790105fe30220000","32c51079a","0000005100000000","0400503000000000"); // 5: alpha < 0.75
    cuts.AddCutHeavyMesonCalo("80010113","411790105fe30220000","32c51079a","0000006100000000","0400503000000000"); // 6: alpha < 0.8
    cuts.AddCutHeavyMesonCalo("80010113","411790105fe30220000","32c51079a","0000007100000000","0400503000000000"); // 7: alpha < 0.85
  }else if (trainConfig == 1110){ // ITS cluster requirement variation
    cuts.AddCutHeavyMesonCalo("80010113","411790105fe30220000","34c51079a","0000003100000000","0400503000000000"); // 4: min 3 ITS cluster
  }else if (trainConfig == 1111){ // TPC cluster requirement variation
    cuts.AddCutHeavyMesonCalo("80010113","411790105fe30220000","32e51079a","0000003100000000","0400503000000000"); // e: No shared clusters
  }else if (trainConfig == 1112){ // Charged pion DCA
    cuts.AddCutHeavyMesonCalo("80010113","411790105fe30220000","32c01079a","0000003100000000","0400503000000000"); // 0: No cut
    cuts.AddCutHeavyMesonCalo("80010113","411790105fe30220000","32c31079a","0000003100000000","0400503000000000"); // 5: Strict pT dep. cut
    cuts.AddCutHeavyMesonCalo("80010113","411790105fe30220000","32c61079a","0000003100000000","0400503000000000"); // 6: z,xy < 0.5 (very tight)
  }else if (trainConfig == 1113){ // TOF requirement
    cuts.AddCutHeavyMesonCalo("80010113","411790105fe30220000","32c51070a","0000003100000000","0400503000000000"); // 0: No TOF
    cuts.AddCutHeavyMesonCalo("80010113","411790105fe30220000","32c51076a","0000003100000000","0400503000000000"); // 6: stricter Kp rejection
    cuts.AddCutHeavyMesonCalo("80010113","411790105fe30220000","32c51072a","0000003100000000","0400503000000000"); // 2: Pion selection
  }else if (trainConfig == 1114){ // min pT
    cuts.AddCutHeavyMesonCalo("80010113","411790105fe30220000","32c50079a","0000003100000000","0400503000000000"); // 0: pT > 0.075 GeV
    cuts.AddCutHeavyMesonCalo("80010113","411790105fe30220000","32c52079a","0000003100000000","0400503000000000"); // 2: pT > 0.125 GeV
    cuts.AddCutHeavyMesonCalo("80010113","411790105fe30220000","32c53079a","0000003100000000","0400503000000000"); // 3: pT > 0.15 GeV
  }else if (trainConfig == 1115){ // TPC dEdx sigma
    cuts.AddCutHeavyMesonCalo("80010113","411790105fe30220000","32c510b9a","0000003100000000","0400503000000000"); // b: -2,2
    cuts.AddCutHeavyMesonCalo("80010113","411790105fe30220000","32c51099a","0000003100000000","0400503000000000"); // 9: -2.5,2.5
    cuts.AddCutHeavyMesonCalo("80010113","411790105fe30220000","32c510a9a","0000003100000000","0400503000000000"); // a: -3.5,3.5
    cuts.AddCutHeavyMesonCalo("80010113","411790105fe30220000","32c51059a","0000003100000000","0400503000000000"); // 5: -4,4
    cuts.AddCutHeavyMesonCalo("80010113","411790105fe30220000","32c51039a","0000003100000000","0400503000000000"); // 3: -5,5
  }else if (trainConfig == 1116){ // PiPlPiMi Mass
    cuts.AddCutHeavyMesonCalo("80010113","411790105fe30220000","32c51079r","0000003100000000","0400503000000000"); // r: 0.8 GeV/c^2
    cuts.AddCutHeavyMesonCalo("80010113","411790105fe30220000","32c51079s","0000003100000000","0400503000000000"); // s: 0.825 GeV/c^2
    cuts.AddCutHeavyMesonCalo("80010113","411790105fe30220000","32c51079t","0000003100000000","0400503000000000"); // t: 0.875 GeV/c^2
    cuts.AddCutHeavyMesonCalo("80010113","411790105fe30220000","32c51079u","0000003100000000","0400503000000000"); // u: 0.9 GeV/c^2
 
  //****************** Minimum Bias - Run 1 (LHC13bcdef) ******************
  }else if (trainConfig == 1200){ // Standard 5 TeV omega EG1 cutstring
    cuts.AddCutHeavyMesonCalo("80010113","111110105fe30220000","32c51072a","0000003100000000","0400503000000000");
  }else if (trainConfig == 1201){ // Non linearity variations
    cuts.AddCutHeavyMesonCalo("80010113","111119705fe30220000","32c51072a","0000003100000000","0400503000000000"); // CRF
    cuts.AddCutHeavyMesonCalo("80010113","111119805fe30220000","32c51072a","0000003100000000","0400503000000000"); // CCRF
  }else if (trainConfig == 1202){ // Cluster timing variations
    cuts.AddCutHeavyMesonCalo("80010113","111110106fe30220000","32c51072a","0000003100000000","0400503000000000"); // 6: -30 - 35
    cuts.AddCutHeavyMesonCalo("80010113","111110107fe30220000","32c51072a","0000003100000000","0400503000000000"); // 7: -30 - 30
    cuts.AddCutHeavyMesonCalo("80010113","111110108fe30220000","32c51072a","0000003100000000","0400503000000000"); // 8: -20 - 30
    cuts.AddCutHeavyMesonCalo("80010113","111110109fe30220000","32c51072a","0000003100000000","0400503000000000"); // 9: -20 - 25
    cuts.AddCutHeavyMesonCalo("80010113","11111010afe30220000","32c51072a","0000003100000000","0400503000000000"); // a: -12.5 - 13
  }else if (trainConfig == 1203){ // Track matching variations
    cuts.AddCutHeavyMesonCalo("80010113","111110105ce30220000","32c51072a","0000003100000000","0400503000000000"); // c: No E/p cut
    cuts.AddCutHeavyMesonCalo("80010113","111110105de30220000","32c51072a","0000003100000000","0400503000000000"); // d: E/p < 3
    cuts.AddCutHeavyMesonCalo("80010113","111110105ee30220000","32c51072a","0000003100000000","0400503000000000"); // e: E/p < 2
    cuts.AddCutHeavyMesonCalo("80010113","111110105ge30220000","32c51072a","0000003100000000","0400503000000000"); // g: E/p < 1.5
    cuts.AddCutHeavyMesonCalo("80010113","111110105he30220000","32c51072a","0000003100000000","0400503000000000"); // h: E/p < 1.25
    cuts.AddCutHeavyMesonCalo("80010113","111110105ne30220000","32c51072a","0000003100000000","0400503000000000"); // o: E/p < 1.75 (standard), but eta and phi varied
    cuts.AddCutHeavyMesonCalo("80010113","111110105oe30220000","32c51072a","0000003100000000","0400503000000000"); // n: E/p < 1.75 (standard), but eta and phi varied
  }else if (trainConfig == 1204){ // Exotic cluster variations
    cuts.AddCutHeavyMesonCalo("80010113","111110105f030220000","32c51072a","0000003100000000","0400503000000000"); // 0: No exotics cut
    cuts.AddCutHeavyMesonCalo("80010113","111110105fb30220000","32c51072a","0000003100000000","0400503000000000"); // b: fExoticEnergyFracCluster = 0.95
    cuts.AddCutHeavyMesonCalo("80010113","111110105fi30220000","32c51072a","0000003100000000","0400503000000000"); // i: fExoticMinEnergyCell = 3
  }else if (trainConfig == 1205){ // Min cluster energy variations
    cuts.AddCutHeavyMesonCalo("80010113","111110105fe20220000","32c51072a","0000003100000000","0400503000000000"); // 2: E > 0.6
    cuts.AddCutHeavyMesonCalo("80010113","111110105fe40220000","32c51072a","0000003100000000","0400503000000000"); // 3: E > 0.8
  }else if (trainConfig == 1206){ // NCell variation
    cuts.AddCutHeavyMesonCalo("80010113","111110105fe3n220000","32c51072a","0000003100000000","0400503000000000");
  }else if (trainConfig == 1207){ // M02 variations
    cuts.AddCutHeavyMesonCalo("80010113","111110105fe30210000","32c51072a","0000003100000000","0400503000000000"); // 1: M02 < 1
    cuts.AddCutHeavyMesonCalo("80010113","111110105fe30230000","32c51072a","0000003100000000","0400503000000000"); // 3: M02 < 0.5
  }else if (trainConfig == 1208){ // pi0 mass selection window variations
    cuts.AddCutHeavyMesonCalo("80010113","111110105fe30220000","32c51072a","0000003r00000000","0400503000000000"); // r: 3.5 sigma, no gamma selection
    cuts.AddCutHeavyMesonCalo("80010113","111110105fe30220000","32c51072a","0000003u00000000","0400503000000000"); // u: 1.5 sigma
    cuts.AddCutHeavyMesonCalo("80010113","111110105fe30220000","32c51072a","0000003v00000000","0400503000000000"); // v: 2.5 sigma
    cuts.AddCutHeavyMesonCalo("80010113","111110105fe30220000","32c51072a","0000003w00000000","0400503000000000"); // w: 4 sigma
    cuts.AddCutHeavyMesonCalo("80010113","111110105fe30220000","32c51072a","0000003x00000000","0400503000000000"); // x: 3 sigma
  }else if (trainConfig == 1209){ // pi0 asymmetry variations
    cuts.AddCutHeavyMesonCalo("80010113","111110105fe30220000","32c51072a","0000005100000000","0400503000000000"); // alpha < 0.75
    cuts.AddCutHeavyMesonCalo("80010113","111110105fe30220000","32c51072a","0000006100000000","0400503000000000"); // alpha < 0.8
    cuts.AddCutHeavyMesonCalo("80010113","111110105fe30220000","32c51072a","0000007100000000","0400503000000000"); // alpha < 0.85
  }else if (trainConfig == 1210){ // ITS cluster requirement variation
    cuts.AddCutHeavyMesonCalo("80010113","111110105fe30220000","34c51072a","0000003100000000","0400503000000000"); // 4: min 3 ITS cluster
  }else if (trainConfig == 1211){ // TPC cluster requirement variation
    cuts.AddCutHeavyMesonCalo("80010113","111110105fe30220000","32e51072a","0000003100000000","0400503000000000"); // e: No shared clusters
  }else if (trainConfig == 1212){ // Charged pion DCA
    cuts.AddCutHeavyMesonCalo("80010113","111110105fe30220000","32c61072a","0000003100000000","0400503000000000"); // 6: z,xy < 0.5 (very tight)
    cuts.AddCutHeavyMesonCalo("80010113","111110105fe30220000","32c01072a","0000003100000000","0400503000000000"); // 0: No cut
  }else if (trainConfig == 1213){ // TOF requirement
    cuts.AddCutHeavyMesonCalo("80010113","111110105fe30220000","32c51070a","0000003100000000","0400503000000000"); // 0: No TOF
    cuts.AddCutHeavyMesonCalo("80010113","111110105fe30220000","32c51073a","0000003100000000","0400503000000000"); // 3: -3<sigma<5
  }else if (trainConfig == 1214){ // min pT
    cuts.AddCutHeavyMesonCalo("80010113","111110105fe30220000","32c50072a","0000003100000000","0400503000000000"); // 0: pT > 0.075 GeV
    cuts.AddCutHeavyMesonCalo("80010113","111110105fe30220000","32c52072a","0000003100000000","0400503000000000"); // 2: pT > 0.125 GeV
    cuts.AddCutHeavyMesonCalo("80010113","111110105fe30220000","32c53072a","0000003100000000","0400503000000000"); // 3: pT > 0.15 GeV
  }else if (trainConfig == 1215){ // TPC dEdx sigma
    cuts.AddCutHeavyMesonCalo("80010113","111110105fe30220000","32c51032a","0000003100000000","0400503000000000"); // 3: -5,5
    cuts.AddCutHeavyMesonCalo("80010113","111110105fe30220000","32c51042a","0000003100000000","0400503000000000"); // 4: -4,5
    cuts.AddCutHeavyMesonCalo("80010113","111110105fe30220000","32c51052a","0000003100000000","0400503000000000"); // 5: -4,4
    cuts.AddCutHeavyMesonCalo("80010113","111110105fe30220000","32c51062a","0000003100000000","0400503000000000"); // 6: -3,4
    cuts.AddCutHeavyMesonCalo("80010113","111110105fe30220000","32c510a2a","0000003100000000","0400503000000000"); // a: -3.5,3.5
    cuts.AddCutHeavyMesonCalo("80010113","111110105fe30220000","32c51092a","0000003100000000","0400503000000000"); // 9: -2.5,2.5
    cuts.AddCutHeavyMesonCalo("80010113","111110105fe30220000","32c510b2a","0000003100000000","0400503000000000"); // b: -2,2
  }else if (trainConfig == 1216){ // PiPlPiMi Mass
    cuts.AddCutHeavyMesonCalo("80010113","111110105fe30220000","32c51072r","0000003100000000","0400503000000000"); // r: 0.8 GeV/c^2
    cuts.AddCutHeavyMesonCalo("80010113","111110105fe30220000","32c51072s","0000003100000000","0400503000000000"); // s: 0.825 GeV/c^2
    cuts.AddCutHeavyMesonCalo("80010113","111110105fe30220000","32c51072t","0000003100000000","0400503000000000"); // t: 0.875 GeV/c^2
    cuts.AddCutHeavyMesonCalo("80010113","111110105fe30220000","32c51072u","0000003100000000","0400503000000000"); // u: 0.9 GeV/c^2

  //****************** EG1 - Run 1 (LHC13bcdef) ******************
  }else if (trainConfig == 1300){ // Standard 5 TeV omega EG1 cutstring
    cuts.AddCutHeavyMesonCalo("80083113","111110105fe30220000","32c51072a","0000003100000000","0400503000000000");
  }else if (trainConfig == 1301){ // Non linearity variations
    cuts.AddCutHeavyMesonCalo("80083113","111119705fe30220000","32c51072a","0000003100000000","0400503000000000"); // CRF
    cuts.AddCutHeavyMesonCalo("80083113","111119805fe30220000","32c51072a","0000003100000000","0400503000000000"); // CCRF
  }else if (trainConfig == 1302){ // Cluster timing variations
    cuts.AddCutHeavyMesonCalo("80083113","111110106fe30220000","32c51072a","0000003100000000","0400503000000000"); // 6: -30 - 35
    cuts.AddCutHeavyMesonCalo("80083113","111110107fe30220000","32c51072a","0000003100000000","0400503000000000"); // 7: -30 - 30
    cuts.AddCutHeavyMesonCalo("80083113","111110108fe30220000","32c51072a","0000003100000000","0400503000000000"); // 8: -20 - 30
    cuts.AddCutHeavyMesonCalo("80083113","111110109fe30220000","32c51072a","0000003100000000","0400503000000000"); // 9: -20 - 25
    cuts.AddCutHeavyMesonCalo("80083113","11111010afe30220000","32c51072a","0000003100000000","0400503000000000"); // a: -12.5 - 13
  }else if (trainConfig == 1303){ // Track matching variations
    cuts.AddCutHeavyMesonCalo("80083113","111110105ce30220000","32c51072a","0000003100000000","0400503000000000"); // c: No E/p cut
    cuts.AddCutHeavyMesonCalo("80083113","111110105de30220000","32c51072a","0000003100000000","0400503000000000"); // d: E/p < 3
    cuts.AddCutHeavyMesonCalo("80083113","111110105ee30220000","32c51072a","0000003100000000","0400503000000000"); // e: E/p < 2
    cuts.AddCutHeavyMesonCalo("80083113","111110105ge30220000","32c51072a","0000003100000000","0400503000000000"); // g: E/p < 1.5
    cuts.AddCutHeavyMesonCalo("80083113","111110105he30220000","32c51072a","0000003100000000","0400503000000000"); // h: E/p < 1.25
    cuts.AddCutHeavyMesonCalo("80083113","111110105ne30220000","32c51072a","0000003100000000","0400503000000000"); // o: E/p < 1.75 (standard), but eta and phi varied
    cuts.AddCutHeavyMesonCalo("80083113","111110105oe30220000","32c51072a","0000003100000000","0400503000000000"); // n: E/p < 1.75 (standard), but eta and phi varied
  }else if (trainConfig == 1304){ // Exotic cluster variations
    cuts.AddCutHeavyMesonCalo("80083113","111110105f030220000","32c51072a","0000003100000000","0400503000000000"); // 0: No exotics cut
    cuts.AddCutHeavyMesonCalo("80083113","111110105fb30220000","32c51072a","0000003100000000","0400503000000000"); // b: fExoticEnergyFracCluster = 0.95
    cuts.AddCutHeavyMesonCalo("80083113","111110105fi30220000","32c51072a","0000003100000000","0400503000000000"); // i: fExoticMinEnergyCell = 3
  }else if (trainConfig == 1305){ // Min cluster energy variations
    cuts.AddCutHeavyMesonCalo("80083113","111110105fe20220000","32c51072a","0000003100000000","0400503000000000"); // 2: E > 0.6
    cuts.AddCutHeavyMesonCalo("80083113","111110105fe40220000","32c51072a","0000003100000000","0400503000000000"); // 3: E > 0.8
  }else if (trainConfig == 1306){ // NCell variation
    cuts.AddCutHeavyMesonCalo("80083113","111110105fe3n220000","32c51072a","0000003100000000","0400503000000000");
  }else if (trainConfig == 1307){ // M02 variations
    cuts.AddCutHeavyMesonCalo("80083113","111110105fe30210000","32c51072a","0000003100000000","0400503000000000"); // 1: M02 < 1
    cuts.AddCutHeavyMesonCalo("80083113","111110105fe30230000","32c51072a","0000003100000000","0400503000000000"); // 3: M02 < 0.5
  }else if (trainConfig == 1308){ // pi0 mass selection window variations
    cuts.AddCutHeavyMesonCalo("80083113","111110105fe30220000","32c51072a","0000003r00000000","0400503000000000"); // r: 3.5 sigma, no gamma selection
    cuts.AddCutHeavyMesonCalo("80083113","111110105fe30220000","32c51072a","0000003u00000000","0400503000000000"); // u: 1.5 sigma
    cuts.AddCutHeavyMesonCalo("80083113","111110105fe30220000","32c51072a","0000003v00000000","0400503000000000"); // v: 2.5 sigma
    cuts.AddCutHeavyMesonCalo("80083113","111110105fe30220000","32c51072a","0000003w00000000","0400503000000000"); // w: 4 sigma
    cuts.AddCutHeavyMesonCalo("80083113","111110105fe30220000","32c51072a","0000003x00000000","0400503000000000"); // x: 3 sigma
  }else if (trainConfig == 1309){ // pi0 asymmetry variations
    cuts.AddCutHeavyMesonCalo("80083113","111110105fe30220000","32c51072a","0000005100000000","0400503000000000"); // alpha < 0.75
    cuts.AddCutHeavyMesonCalo("80083113","111110105fe30220000","32c51072a","0000006100000000","0400503000000000"); // alpha < 0.8
    cuts.AddCutHeavyMesonCalo("80083113","111110105fe30220000","32c51072a","0000007100000000","0400503000000000"); // alpha < 0.85
  }else if (trainConfig == 1310){ // ITS cluster requirement variation
    cuts.AddCutHeavyMesonCalo("80083113","111110105fe30220000","34c51072a","0000003100000000","0400503000000000"); // 4: min 3 ITS cluster
  }else if (trainConfig == 1311){ // TPC cluster requirement variation
    cuts.AddCutHeavyMesonCalo("80083113","111110105fe30220000","32e51072a","0000003100000000","0400503000000000"); // e: No shared clusters
  }else if (trainConfig == 1312){ // Charged pion DCA
    cuts.AddCutHeavyMesonCalo("80083113","111110105fe30220000","32c61072a","0000003100000000","0400503000000000"); // 6: z,xy < 0.5 (very tight)
    cuts.AddCutHeavyMesonCalo("80083113","111110105fe30220000","32c01072a","0000003100000000","0400503000000000"); // 0: No cut
  }else if (trainConfig == 1313){ // TOF requirement
    cuts.AddCutHeavyMesonCalo("80083113","111110105fe30220000","32c51070a","0000003100000000","0400503000000000"); // 0: No TOF
    cuts.AddCutHeavyMesonCalo("80083113","111110105fe30220000","32c51073a","0000003100000000","0400503000000000"); // 3: -3<sigma<5
  }else if (trainConfig == 1314){ // min pT
    cuts.AddCutHeavyMesonCalo("80083113","111110105fe30220000","32c50072a","0000003100000000","0400503000000000"); // 0: pT > 0.075 GeV
    cuts.AddCutHeavyMesonCalo("80083113","111110105fe30220000","32c52072a","0000003100000000","0400503000000000"); // 2: pT > 0.125 GeV
    cuts.AddCutHeavyMesonCalo("80083113","111110105fe30220000","32c53072a","0000003100000000","0400503000000000"); // 3: pT > 0.15 GeV
  }else if (trainConfig == 1315){ // TPC dEdx sigma
    cuts.AddCutHeavyMesonCalo("80083113","111110105fe30220000","32c51032a","0000003100000000","0400503000000000"); // 3: -5,5
    cuts.AddCutHeavyMesonCalo("80083113","111110105fe30220000","32c51042a","0000003100000000","0400503000000000"); // 4: -4,5
    cuts.AddCutHeavyMesonCalo("80083113","111110105fe30220000","32c51052a","0000003100000000","0400503000000000"); // 5: -4,4
    cuts.AddCutHeavyMesonCalo("80083113","111110105fe30220000","32c51062a","0000003100000000","0400503000000000"); // 6: -3,4
    cuts.AddCutHeavyMesonCalo("80083113","111110105fe30220000","32c510a2a","0000003100000000","0400503000000000"); // a: -3.5,3.5
    cuts.AddCutHeavyMesonCalo("80083113","111110105fe30220000","32c51092a","0000003100000000","0400503000000000"); // 9: -2.5,2.5
    cuts.AddCutHeavyMesonCalo("80083113","111110105fe30220000","32c510b2a","0000003100000000","0400503000000000"); // b: -2,2
  }else if (trainConfig == 1316){ // PiPlPiMi Mass
    cuts.AddCutHeavyMesonCalo("80083113","111110105fe30220000","32c51072r","0000003100000000","0400503000000000"); // r: 0.8 GeV/c^2
    cuts.AddCutHeavyMesonCalo("80083113","111110105fe30220000","32c51072s","0000003100000000","0400503000000000"); // s: 0.825 GeV/c^2
    cuts.AddCutHeavyMesonCalo("80083113","111110105fe30220000","32c51072t","0000003100000000","0400503000000000"); // t: 0.875 GeV/c^2
    cuts.AddCutHeavyMesonCalo("80083113","111110105fe30220000","32c51072u","0000003100000000","0400503000000000"); // u: 0.9 GeV/c^2

  //****************** EG2 - Run 1 (LHC13bcdef) ******************
  }else if (trainConfig == 1400){ // Standard 5 TeV omega EG2 cutstring
    cuts.AddCutHeavyMesonCalo("80085113","111110105fe30220000","32c51072a","0000003100000000","0400503000000000");
  }else if (trainConfig == 1401){ // Non linearity variations
    cuts.AddCutHeavyMesonCalo("80085113","111119705fe30220000","32c51072a","0000003100000000","0400503000000000"); // CRF
    cuts.AddCutHeavyMesonCalo("80085113","111119805fe30220000","32c51072a","0000003100000000","0400503000000000"); // CCRF
  }else if (trainConfig == 1402){ // Cluster timing variations
    cuts.AddCutHeavyMesonCalo("80085113","111110106fe30220000","32c51072a","0000003100000000","0400503000000000"); // 6: -30 - 35
    cuts.AddCutHeavyMesonCalo("80085113","111110107fe30220000","32c51072a","0000003100000000","0400503000000000"); // 7: -30 - 30
    cuts.AddCutHeavyMesonCalo("80085113","111110108fe30220000","32c51072a","0000003100000000","0400503000000000"); // 8: -20 - 30
    cuts.AddCutHeavyMesonCalo("80085113","111110109fe30220000","32c51072a","0000003100000000","0400503000000000"); // 9: -20 - 25
    cuts.AddCutHeavyMesonCalo("80085113","11111010afe30220000","32c51072a","0000003100000000","0400503000000000"); // a: -12.5 - 13
  }else if (trainConfig == 1403){ // Track matching variations
    cuts.AddCutHeavyMesonCalo("80085113","111110105ce30220000","32c51072a","0000003100000000","0400503000000000"); // c: No E/p cut
    cuts.AddCutHeavyMesonCalo("80085113","111110105de30220000","32c51072a","0000003100000000","0400503000000000"); // d: E/p < 3
    cuts.AddCutHeavyMesonCalo("80085113","111110105ee30220000","32c51072a","0000003100000000","0400503000000000"); // e: E/p < 2
    cuts.AddCutHeavyMesonCalo("80085113","111110105ge30220000","32c51072a","0000003100000000","0400503000000000"); // g: E/p < 1.5
    cuts.AddCutHeavyMesonCalo("80085113","111110105he30220000","32c51072a","0000003100000000","0400503000000000"); // h: E/p < 1.25
    cuts.AddCutHeavyMesonCalo("80085113","111110105ne30220000","32c51072a","0000003100000000","0400503000000000"); // o: E/p < 1.75 (standard), but eta and phi varied
    cuts.AddCutHeavyMesonCalo("80085113","111110105oe30220000","32c51072a","0000003100000000","0400503000000000"); // n: E/p < 1.75 (standard), but eta and phi varied
  }else if (trainConfig == 1404){ // Exotic cluster variations
    cuts.AddCutHeavyMesonCalo("80085113","111110105f030220000","32c51072a","0000003100000000","0400503000000000"); // 0: No exotics cut
    cuts.AddCutHeavyMesonCalo("80085113","111110105fb30220000","32c51072a","0000003100000000","0400503000000000"); // b: fExoticEnergyFracCluster = 0.95
    cuts.AddCutHeavyMesonCalo("80085113","111110105fi30220000","32c51072a","0000003100000000","0400503000000000"); // i: fExoticMinEnergyCell = 3
  }else if (trainConfig == 1405){ // Min cluster energy variations
    cuts.AddCutHeavyMesonCalo("80085113","111110105fe20220000","32c51072a","0000003100000000","0400503000000000"); // 2: E > 0.6
    cuts.AddCutHeavyMesonCalo("80085113","111110105fe40220000","32c51072a","0000003100000000","0400503000000000"); // 3: E > 0.8
  }else if (trainConfig == 1406){ // NCell variation
    cuts.AddCutHeavyMesonCalo("80085113","111110105fe3n220000","32c51072a","0000003100000000","0400503000000000");
  }else if (trainConfig == 1407){ // M02 variations
    cuts.AddCutHeavyMesonCalo("80085113","111110105fe30210000","32c51072a","0000003100000000","0400503000000000"); // 1: M02 < 1
    cuts.AddCutHeavyMesonCalo("80085113","111110105fe30230000","32c51072a","0000003100000000","0400503000000000"); // 3: M02 < 0.5
  }else if (trainConfig == 1408){ // pi0 mass selection window variations
    cuts.AddCutHeavyMesonCalo("80085113","111110105fe30220000","32c51072a","0000003r00000000","0400503000000000"); // r: 3.5 sigma, no gamma selection
    cuts.AddCutHeavyMesonCalo("80085113","111110105fe30220000","32c51072a","0000003u00000000","0400503000000000"); // u: 1.5 sigma
    cuts.AddCutHeavyMesonCalo("80085113","111110105fe30220000","32c51072a","0000003v00000000","0400503000000000"); // v: 2.5 sigma
    cuts.AddCutHeavyMesonCalo("80085113","111110105fe30220000","32c51072a","0000003w00000000","0400503000000000"); // w: 4 sigma
    cuts.AddCutHeavyMesonCalo("80085113","111110105fe30220000","32c51072a","0000003x00000000","0400503000000000"); // x: 3 sigma
  }else if (trainConfig == 1409){ // pi0 asymmetry variations
    cuts.AddCutHeavyMesonCalo("80085113","111110105fe30220000","32c51072a","0000005100000000","0400503000000000"); // alpha < 0.75
    cuts.AddCutHeavyMesonCalo("80085113","111110105fe30220000","32c51072a","0000006100000000","0400503000000000"); // alpha < 0.8
    cuts.AddCutHeavyMesonCalo("80085113","111110105fe30220000","32c51072a","0000007100000000","0400503000000000"); // alpha < 0.85
  }else if (trainConfig == 1410){ // ITS cluster requirement variation
    cuts.AddCutHeavyMesonCalo("80085113","111110105fe30220000","34c51072a","0000003100000000","0400503000000000"); // 4: min 3 ITS cluster
  }else if (trainConfig == 1411){ // TPC cluster requirement variation
    cuts.AddCutHeavyMesonCalo("80085113","111110105fe30220000","32e51072a","0000003100000000","0400503000000000"); // e: No shared clusters
  }else if (trainConfig == 1412){ // Charged pion DCA
    cuts.AddCutHeavyMesonCalo("80085113","111110105fe30220000","32c61072a","0000003100000000","0400503000000000"); // 6: z,xy < 0.5 (very tight)
    cuts.AddCutHeavyMesonCalo("80085113","111110105fe30220000","32c01072a","0000003100000000","0400503000000000"); // 0: No cut
  }else if (trainConfig == 1413){ // TOF requirement
    cuts.AddCutHeavyMesonCalo("80085113","111110105fe30220000","32c51070a","0000003100000000","0400503000000000"); // 0: No TOF
    cuts.AddCutHeavyMesonCalo("80085113","111110105fe30220000","32c51073a","0000003100000000","0400503000000000"); // 3: -3<sigma<5
  }else if (trainConfig == 1414){ // min pT
    cuts.AddCutHeavyMesonCalo("80085113","111110105fe30220000","32c50072a","0000003100000000","0400503000000000"); // 0: pT > 0.075 GeV
    cuts.AddCutHeavyMesonCalo("80085113","111110105fe30220000","32c52072a","0000003100000000","0400503000000000"); // 2: pT > 0.125 GeV
    cuts.AddCutHeavyMesonCalo("80085113","111110105fe30220000","32c53072a","0000003100000000","0400503000000000"); // 3: pT > 0.15 GeV
  }else if (trainConfig == 1415){ // TPC dEdx sigma
    cuts.AddCutHeavyMesonCalo("80085113","111110105fe30220000","32c51032a","0000003100000000","0400503000000000"); // 3: -5,5
    cuts.AddCutHeavyMesonCalo("80085113","111110105fe30220000","32c51042a","0000003100000000","0400503000000000"); // 4: -4,5
    cuts.AddCutHeavyMesonCalo("80085113","111110105fe30220000","32c51052a","0000003100000000","0400503000000000"); // 5: -4,4
    cuts.AddCutHeavyMesonCalo("80085113","111110105fe30220000","32c51062a","0000003100000000","0400503000000000"); // 6: -3,4
    cuts.AddCutHeavyMesonCalo("80085113","111110105fe30220000","32c510a2a","0000003100000000","0400503000000000"); // a: -3.5,3.5
    cuts.AddCutHeavyMesonCalo("80085113","111110105fe30220000","32c51092a","0000003100000000","0400503000000000"); // 9: -2.5,2.5
    cuts.AddCutHeavyMesonCalo("80085113","111110105fe30220000","32c510b2a","0000003100000000","0400503000000000"); // b: -2,2
  }else if (trainConfig == 1416){ // PiPlPiMi Mass
    cuts.AddCutHeavyMesonCalo("80085113","111110105fe30220000","32c51072r","0000003100000000","0400503000000000"); // r: 0.8 GeV/c^2
    cuts.AddCutHeavyMesonCalo("80085113","111110105fe30220000","32c51072s","0000003100000000","0400503000000000"); // s: 0.825 GeV/c^2
    cuts.AddCutHeavyMesonCalo("80085113","111110105fe30220000","32c51072t","0000003100000000","0400503000000000"); // t: 0.875 GeV/c^2
    cuts.AddCutHeavyMesonCalo("80085113","111110105fe30220000","32c51072u","0000003100000000","0400503000000000"); // u: 0.9 GeV/c^2



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

  } else if (trainConfig == 1507) {  // PHOS  INT7 run2
    cuts.AddCutHeavyMesonCalo("80010113","24466190sa01cc00000","32c51070a","0103603q00000000","0453503000000000");  // 0-100% NL19

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


  //************************************************ PCM- EDC analysis 13 TeV pPb *********************************************
  } else if (trainConfig == 2010 ){ //EDC 13TeV MB, NCell: v (NCell Cut 2, but with probability in MC to let clusters through )
    cuts.AddCutHeavyMesonCalo("80010113","411790109fe3n230000","32c51070m","0103603l00000000","0453503000000000"); // INT7
  } else if (trainConfig == 2011){ //EDC 13TeV MB, NCell: v (NCell Cut 2, but with probability in MC to let clusters through )
    cuts.AddCutHeavyMesonCalo("8008e113","411790109fe3n230000","32c51070m","0103603l00000000","0453503000000000"); // EG2
  } else if (trainConfig == 2012){ //EDC 13TeV MB, NCell: v (NCell Cut 2, but with probability in MC to let clusters through )
    cuts.AddCutHeavyMesonCalo("8008d113","411790109fe3n230000","32c51070m","0103603l00000000","0453503000000000"); // EG1
  } else if (trainConfig == 2013 ){ //EDC 13TeV MB, NCell: v (NCell Cut 2, but with probability in MC to let clusters through ), no NL
    cuts.AddCutHeavyMesonCalo("80010113","411790009fe3n230000","32c51070m","0103603l00000000","0453503000000000"); // INT7
  } else if (trainConfig == 2014){ //EDC 13TeV MB, NCell: v (NCell Cut 2, but with probability in MC to let clusters through ), no NL
    cuts.AddCutHeavyMesonCalo("8008e113","411790009fe3n230000","32c51070m","0103603l00000000","0453503000000000"); // EG2
  } else if (trainConfig == 2015){ //EDC 13TeV MB, NCell: v (NCell Cut 2, but with probability in MC to let clusters through ), no NL
    cuts.AddCutHeavyMesonCalo("8008d113","411790009fe3n230000","32c51070m","0103603l00000000","0453503000000000"); // EG1
  } else {
    Error(Form("GammaConvNeutralMeson_CaloMode_%i",trainConfig), "wrong trainConfig variable no cuts have been specified for the configuration");
    return;
  }

  if(!cuts.AreValid()){
    std::cout << "\n\n****************************************************" << std::endl;
    std::cout << "ERROR: No valid cuts stored in CutHandlerNeutralCalo! Returning..." << std::endl;
    std::cout << "****************************************************\n\n" << std::endl;
    return ;
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

    analysisEventCuts[i] = new AliConvEventCuts();
    analysisEventCuts[i]->SetV0ReaderName(V0ReaderName);
    analysisEventCuts[i]->SetCorrectionTaskSetting(corrTaskSetting);
    analysisEventCuts[i]->SetTriggerMimicking(enableTriggerMimicking);
    analysisEventCuts[i]->SetTriggerOverlapRejecion(enableTriggerOverlapRej);
    if(runLightOutput>0) analysisEventCuts[i]->SetLightOutput(kTRUE);
    analysisEventCuts[i]->InitializeCutsFromCutString((cuts.GetEventCut(i)).Data());
    if (periodNameV0Reader.CompareTo("") != 0) analysisEventCuts[i]->SetPeriodEnum(periodNameV0Reader);
    EventCutList->Add(analysisEventCuts[i]);
    analysisEventCuts[i]->SetFillCutHistograms("",kFALSE);


    analysisClusterCuts[i] = new AliCaloPhotonCuts();
    analysisClusterCuts[i]->SetV0ReaderName(V0ReaderName);
    analysisClusterCuts[i]->SetCorrectionTaskSetting(corrTaskSetting);
    if(runLightOutput>0) analysisClusterCuts[i]->SetLightOutput(kTRUE);

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
    if(runLightOutput>0) analysisNeutralPionCuts[i]->SetLightOutput(kTRUE);
    if( ! analysisNeutralPionCuts[i]->InitializeCutsFromCutString((cuts.GetNDMCut(i)).Data()) ) {
      std::cout<<"ERROR: analysisMesonCuts [ " <<i<<" ] "<<std::endl;
      return ;
    } else {
      NeutralPionCutList->Add(analysisNeutralPionCuts[i]);
      analysisNeutralPionCuts[i]->SetFillCutHistograms("");
    }

    analysisMesonCuts[i] = new AliConversionMesonCuts();
    if(runLightOutput>0) analysisMesonCuts[i]->SetLightOutput(kTRUE);
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
    if(runLightOutput>0) analysisPionCuts[i]->SetLightOutput(kTRUE);

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
  AliAnalysisDataContainer *coutput =
  mgr->CreateContainer(!(corrTaskSetting.CompareTo("")) ? Form("GammaConvNeutralMesonPiPlPiMiNeutralMeson_%i_%i_%i.root",selectHeavyNeutralMeson,neutralPionMode, trainConfig) : Form("GammaConvNeutralMesonPiPlPiMiNeutralMeson_%i_%i_%i_%s.root",selectHeavyNeutralMeson,neutralPionMode, trainConfig, corrTaskSetting.Data()), TList::Class(),
              AliAnalysisManager::kOutputContainer,Form("GammaConvNeutralMesonPiPlPiMiNeutralMeson_%i_%i_%i.root",selectHeavyNeutralMeson,neutralPionMode, trainConfig));

  mgr->AddTask(task);
  mgr->ConnectInput(task,0,cinput);
  mgr->ConnectOutput(task,1,coutput);

  return;

}
