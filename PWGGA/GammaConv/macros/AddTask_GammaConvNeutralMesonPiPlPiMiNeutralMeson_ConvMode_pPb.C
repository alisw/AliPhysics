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
void AddTask_GammaConvNeutralMesonPiPlPiMiNeutralMeson_ConvMode_pPb(
    Int_t     trainConfig                   = 1,
    Int_t     isMC                          = 0,                        //run MC
    TString   photonCutNumberV0Reader       = "",                       // 00000008400000000100000000 nom. B, 00000088400000000100000000 low B
    Int_t     selectHeavyNeutralMeson       = 0,                        //run eta prime instead of omega
    Int_t     enableQAMesonTask             = 1,                        //enable QA in AliAnalysisTaskNeutralMesonToPiPlPiMiNeutralMeson
    Int_t     enableTriggerMimicking        = 0,                        // enable trigger mimicking
    Bool_t    enableTriggerOverlapRej       = kFALSE,                   // enable trigger overlap rejection
    TString   fileNameExternalInputs        = "MCSpectraInput.root",    // path to file for weigting input
    Int_t     doWeighting                   = kFALSE,                       //enable Weighting
    Bool_t    enableElecDeDxPostCalibration = kFALSE,                 // enable post calibration of elec pos dEdX
    TString   generatorName                 = "HIJING",
    Double_t  tolerance                     = -1,
    TString   periodNameV0Reader            = "",                       // period Name for V0Reader
    Int_t     runLightOutput                = 0,                        // run light output option 0: no light output 1: most cut histos stiched off 2: unecessary omega hists turned off as well
    Int_t     prefilterRunFlag              = 1500,                     // flag to change the prefiltering of ESD tracks. See SetHybridTrackCutsAODFiltering() in AliPrimaryPionCuts
    Bool_t    usePtDepSelectionWindowCut    = kFALSE,                   // use pt dependent meson selection window cut
    Int_t     enableMatBudWeightsPi0        = 0,                        // 1 = three radial bins, 2 = 10 radial bins (2 is the default when using weights)
    TString   additionalTrainConfig         = "0"                       // additional counter for trainconfig, this has to be always the last parameter
  ){

  AliCutHandlerPCM cuts(13);
  TString addTaskName                       = "AddTask_GammaConvNeutralMesonPiPlPiMiNeutralMeson_ConvMode_pPb";
  TString fileNamePtWeights                 = cuts.GetSpecialFileNameFromString (fileNameExternalInputs, "FPTW:");
  TString fileNameMultWeights               = cuts.GetSpecialFileNameFromString (fileNameExternalInputs, "FMUW:");
  TString fileNameMatBudWeights             = cuts.GetSpecialFileNameFromString (fileNameExternalInputs, "FMAW:");
  TString fileNamedEdxPostCalib             = cuts.GetSpecialFileNameFromString (fileNameExternalInputs, "FEPC:");
  TString fileNameCustomTriggerMimicOADB    = cuts.GetSpecialFileNameFromString (fileNameExternalInputs, "FTRM:");

  if(additionalTrainConfig.Contains("MaterialBudgetWeights"))
    fileNameMatBudWeights         = cuts.GetSpecialSettingFromAddConfig(additionalTrainConfig, "MaterialBudgetWeights",fileNameMatBudWeights, addTaskName);
  
  //parse additionalTrainConfig flag
  TString unsmearingoutputs = "012"; // 0: No correction, 1: One pi0 mass errer subtracted, 2: pz of pi0 corrected to fix its mass, 3: Lambda(alpha)*DeltaPi0 subtracted

  TObjArray *rAddConfigArr = additionalTrainConfig.Tokenize("_");
  if(rAddConfigArr->GetEntries()<1){cout << "ERROR during parsing of additionalTrainConfig String '" << additionalTrainConfig.Data() << "'" << endl; return;}
  TObjString* rAdditionalTrainConfig;
  for(Int_t i = 0; i<rAddConfigArr->GetEntries() ; i++){
    if(i==0) rAdditionalTrainConfig = (TObjString*)rAddConfigArr->At(i);
    else{
      TObjString* temp = (TObjString*) rAddConfigArr->At(i);
      TString tempStr = temp->GetString();
      if(tempStr.BeginsWith("UNSMEARING")){ // 0: No correction, 1: One pi0 mass errer subtracted, 2: pz of pi0 corrected to fix its mass, 3: Lambda(alpha)*DeltaPi0 subtracted
        TString tempType = tempStr;
        tempType.Replace(0,9,"");
        unsmearingoutputs = tempType;
        cout << "INFO: AddTask_GammaConvNeutralMesonPiPlPiMiPiZero_ConvMode_pPb will output the following minv_pT histograms:" << endl;
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
    cout << "INFO: AddTask_GammaConvNeutralMesonPiPlPiMiNeutralMeson_ConvMode_pPb running additionalTrainConfig '" << sAdditionalTrainConfig.Atoi() << "', train config: '" << trainConfig << "'" << endl;
  }

  Int_t isHeavyIon = 2;
  Int_t neutralPionMode = 0;

  // ================== GetAnalysisManager ===============================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error(Form("AddTask_GammaConvNeutralMesonPiPlPiMiNeutralMeson_ConvMode_pPb_%i",trainConfig), "No analysis manager found.");
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
  if (trainConfig == 1){ // INT7 run1 & run2
    cuts.AddCutHeavyMesonPCM("80010113","00200009f9730000dge0400000","32c51070a","0103603s00000000","0453503000000000"); // 0-100%, INT7
  } else if (trainConfig == 2){ // EMC EMC triggers
    cuts.AddCutHeavyMesonPCM("80052113","00200009f9730000dge0400000","32c51070a","0103603s00000000","0453503000000000"); // 0-100%, EMC7
  } else if (trainConfig == 3){ // EMC EMC triggers
    cuts.AddCutHeavyMesonPCM("80083113","00200009f9730000dge0400000","32c51070a","0103603s00000000","0453503000000000"); // 0-100%, EG1
  } else if (trainConfig == 4){ // EMC EMC triggers
    cuts.AddCutHeavyMesonPCM("80085113","00200009f9730000dge0400000","32c51070a","0103603s00000000","0453503000000000"); // 0-100%, EG2

    // EMCal + EDC triggers
  } else if (trainConfig == 5){ // EMC EMC triggers
    cuts.AddCutHeavyMesonPCM("8008d113","00200009f9730000dge0400000","32c51070a","0103603s00000000","0453503000000000"); // 0-100%, EG1
  } else if (trainConfig == 6){ // EMC EMC triggers
    cuts.AddCutHeavyMesonPCM("8008e113","00200009f9730000dge0400000","32c51070a","0103603s00000000","0453503000000000"); // 0-100%, EG2

 //************************************************ PCM- PHOS analysis 5 TeV pPb ********************************************
  } else if (trainConfig == 501){ // PHOS  PHI7 run1
    cuts.AddCutHeavyMesonPCM("80062113","00200009f9730000dge0400000","32c51070a","0103603s00000000","0453503000000000");  // 0-100%

    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //                                          OMEGA MESON
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  //************************************************ PCM- EDC analysis 5 TeV pPb *********************************************

  // no event mixing background
  }else if (trainConfig == 1001){ // EMC  INT7 run1 & run2
    cuts.AddCutHeavyMesonPCM("80010113","00200009f9730000dge0400000","32c51070a","0103603s00000000","0453503000000000"); // 0-100%
  } else if (trainConfig == 1002){ // EMC EMC triggers
    cuts.AddCutHeavyMesonPCM("80052113","00200009f9730000dge0400000","32c51070a","0103603s00000000","0453503000000000"); // 0-100%, EMC7
  } else if (trainConfig == 1003){ // EMC EMC triggers
    cuts.AddCutHeavyMesonPCM("80083113","00200009f9730000dge0400000","32c51070a","0103603s00000000","0453503000000000"); // 0-100%, EG1
  } else if (trainConfig == 1004){ // EMC EMC triggers
    cuts.AddCutHeavyMesonPCM("80085113","00200009f9730000dge0400000","32c51070a","0103603s00000000","0453503000000000"); // 0-100%, EG2

    // EMCal + EDC triggers
  } else if (trainConfig == 1005){ // EMC EMC triggers
    cuts.AddCutHeavyMesonPCM("8008d113","00200009f9730000dge0400000","32c51070a","0103603s00000000","0453503000000000"); // 0-100%, EG1
  } else if (trainConfig == 1006){ // EMC EMC triggers
    cuts.AddCutHeavyMesonPCM("8008e113","00200009f9730000dge0400000","32c51070a","0103603s00000000","0453503000000000"); // 0-100%, EG2

  }else if (trainConfig == 1007){ // PCM  INT7 standard cut study guesstimate
    cuts.AddCutHeavyMesonPCM("80010113","00200009227000008250400000","32c51070a","0103603o00000000","0453503000000000"); //  First converstion cut guesstimate
  }else if (trainConfig == 1008){ // PCM  INT7 standard cut study guesstimate
    cuts.AddCutHeavyMesonPCM("80010113","0dm00009f9730000dge0404000","32c51070a","0103103500000000","0153503000000000"); // Second converstion cut guesstimate
  }else if (trainConfig == 1009){ // PCM  INT7 with restricted DCA
    cuts.AddCutHeavyMesonPCM("80010113","0dm00009f9730000dge0404000","32c31070a","0103103500000000","0153503000000000"); // DCA resticted
  }else if (trainConfig == 1010){ // PCM  INT7 with different pi0 selection windows
    cuts.AddCutHeavyMesonPCM("80010113","0dm00009f9730000dge0404000","32c51070a","0103103s00000000","0153503000000000"); // PCM 2 sigma
    cuts.AddCutHeavyMesonPCM("80010113","0dm00009f9730000dge0404000","32c51070a","0103113s00000000","0153503000000000"); // PCM 2 sigma, pT > 0.4
    cuts.AddCutHeavyMesonPCM("80010113","0dm00009f9730000dge0404000","32c51070a","0103123s00000000","0153503000000000"); // PCM 2 sigma, pT > 0.7
    cuts.AddCutHeavyMesonPCM("80010113","0dm00009f9730000dge0404000","32c51070a","0103133s00000000","0153503000000000"); // PCM 2 sigma, pT > 0.9
  }else if (trainConfig == 1011){ // TOF variations to increase charged pion purity
    cuts.AddCutHeavyMesonPCM("80010113","0dm00009f9730000dge0404000","32c51076a","0000003z00000000","0400503000000000"); // TOF Kaon Proton rejection
    cuts.AddCutHeavyMesonPCM("80010113","0dm00009f9730000dge0404000","32c51070a","0000003z00000000","0400503000000000"); // No cut
    cuts.AddCutHeavyMesonPCM("80010113","0dm00009f9730000dge0404000","32c51072a","0000003z00000000","0400503000000000"); // TOF Pion selection
    cuts.AddCutHeavyMesonPCM("80010113","0dm00009f9730000dge0404000","32c51077a","0000003z00000000","0400503000000000"); // TOF Kaon Proton rejection
    cuts.AddCutHeavyMesonPCM("80010113","0dm00009f9730000dge0404000","32c51078a","0000003z00000000","0400503000000000"); // TOF Kaon Proton rejection
    cuts.AddCutHeavyMesonPCM("80010113","0dm00009f9730000dge0404000","32c51079a","0000003z00000000","0400503000000000"); // TOF Kaon Proton rejection
  }else if (trainConfig == 1012){ // DCA variations to increase charged pion purity
    cuts.AddCutHeavyMesonPCM("80010113","0dm00009f9730000dge0404000","32c41076a","0000003z00000000","0400503000000000"); // TOF Kaon Proton rejection + DCA 0.5/3
    cuts.AddCutHeavyMesonPCM("80010113","0dm00009f9730000dge0404000","32c31076a","0000003z00000000","0400503000000000"); // TOF Kaon Proton rejection + DCA pTdep/3
    cuts.AddCutHeavyMesonPCM("80010113","0dm00009f9730000dge0404000","32c71076a","0000003z00000000","0400503000000000"); // TOF Kaon Proton rejection + DCA pTdep/pTdep


    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //                                        OMEGA MESON (p-Pb @ 5 TeV)
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  }else if (trainConfig == 1100){ // Standard 5 TeV omega INT7 cutstring
    cuts.AddCutHeavyMesonPCM("80010113","0dm00009f9730000dge0404000","32c51079a","0000003z00000000","0400503000000000");
  }else if (trainConfig == 1101){ // Min track pT variations
    cuts.AddCutHeavyMesonPCM("80010113","0dm00069f9730000dge0404000","32c51079a","0000003z00000000","0400503000000000"); // 6: e_pT > 40 MeV
    cuts.AddCutHeavyMesonPCM("80010113","0dm00049f9730000dge0404000","32c51079a","0000003z00000000","0400503000000000"); // 4: e_pT > 75 MeV
    cuts.AddCutHeavyMesonPCM("80010113","0dm00019f9730000dge0404000","32c51079a","0000003z00000000","0400503000000000"); // 1: e_pT > 100 MeV
  }else if (trainConfig == 1102){ // min TPC cluster variations
    cuts.AddCutHeavyMesonPCM("80010113","0dm00008f8730000dge0404000","32c51079a","0000003z00000000","0400503000000000"); // 8: > 35% of findable
    cuts.AddCutHeavyMesonPCM("80010113","0dm00006f6730000dge0404000","32c51079a","0000003z00000000","0400503000000000"); // 6: > 70% of findable
  }else if (trainConfig == 1103){ // electron PID variations
    cuts.AddCutHeavyMesonPCM("80010113","0dm0000939730000dge0404000","32c51079a","0000003z00000000","0400503000000000"); // 3: -4 - 5
    cuts.AddCutHeavyMesonPCM("80010113","0dm0000969730000dge0404000","32c51079a","0000003z00000000","0400503000000000"); // 6: -2.5 - 4
  }else if (trainConfig == 1104){ // pion rejection variations
    cuts.AddCutHeavyMesonPCM("80010113","0dm00009f1730000dge0404000","32c51079a","0000003z00000000","0400503000000000"); // 1: nsigma > 0
    cuts.AddCutHeavyMesonPCM("80010113","0dm00009f5730000dge0404000","32c51079a","0000003z00000000","0400503000000000"); // 5: nsigma > 2
  }else if (trainConfig == 1105){ // qT variations
    cuts.AddCutHeavyMesonPCM("80010113","0dm00009f9730000lge0404000","32c51079a","0000003z00000000","0400503000000000"); // l: const = 0.03, pT = 0.11
    cuts.AddCutHeavyMesonPCM("80010113","0dm00009f9730000jge0404000","32c51079a","0000003z00000000","0400503000000000"); // j: const = 0.04, pT = 0.25
    cuts.AddCutHeavyMesonPCM("80010113","0dm00009f9730000kge0404000","32c51079a","0000003z00000000","0400503000000000"); // k: const = 0.045, pT = 0.3
    cuts.AddCutHeavyMesonPCM("80010113","0dm00009f9730000ege0404000","32c51079a","0000003z00000000","0400503000000000"); // e: const = 0.06, pT = 0.14
    cuts.AddCutHeavyMesonPCM("80010113","0dm00009f9730000fge0404000","32c51079a","0000003z00000000","0400503000000000"); // f: const = 0.07, pT = 0.16
  }else if (trainConfig == 1106){ // PsiPair variations
    cuts.AddCutHeavyMesonPCM("80010113","0dm00009f9730000de20404000","32c51079a","0000003z00000000","0400503000000000"); // e2: const PsiPair
    cuts.AddCutHeavyMesonPCM("80010113","0dm00009f9730000d290404000","32c51079a","0000003z00000000","0400503000000000"); // 29 : chi2 = 30 PsiPair = 0.1 triangle (Nico)
    cuts.AddCutHeavyMesonPCM("80010113","0dm00009f9730000d860404000","32c51079a","0000003z00000000","0400503000000000"); // 86 : chi2 = 20 PsiPair = 0.05 triangle (Nico var)
    cuts.AddCutHeavyMesonPCM("80010113","0dm00009f9730000dfd0404000","32c51079a","0000003z00000000","0400503000000000"); // fd: 0.065 exp, chi = 0.18
    cuts.AddCutHeavyMesonPCM("80010113","0dm00009f9730000dif0404000","32c51079a","0000003z00000000","0400503000000000"); // if: 0.075 exp, chi = 0.2
  }else if (trainConfig == 1107){ // Cos point angle
    cuts.AddCutHeavyMesonPCM("80010113","0dm00009f9730000dge0304000","32c51079a","0000003z00000000","0400503000000000"); // 3: 0.75
    cuts.AddCutHeavyMesonPCM("80010113","0dm00009f9730000dge0504000","32c51079a","0000003z00000000","0400503000000000"); // 5: 0.88
  }else if (trainConfig == 1108){ // pi0 mass selection window variations
    cuts.AddCutHeavyMesonPCM("80010113","0dm00009f9730000dge0404000","32c51079a","0000003500000000","0400503000000000"); // 5: 2.0 sigma, no gamma selection
    cuts.AddCutHeavyMesonPCM("80010113","0dm00009f9730000dge0404000","32c51079a","0000003o00000000","0400503000000000"); // o: 2.5 sigma, no gamma selection
    cuts.AddCutHeavyMesonPCM("80010113","0dm00009f9730000dge0404000","32c51079a","0000003n00000000","0400503000000000"); // n: 3.5 sigma, no gamma selection
    cuts.AddCutHeavyMesonPCM("80010113","0dm00009f9730000dge0404000","32c51079a","0000003p00000000","0400503000000000"); // p: 4 sigma, no gamma selection
  }else if (trainConfig == 1109){ // pi0 asymmetry variations
    cuts.AddCutHeavyMesonPCM("80010113","0dm00009f9730000dge0404000","32c51079a","0000005z00000000","0400503000000000"); // 5: alpha < 0.75
    cuts.AddCutHeavyMesonPCM("80010113","0dm00009f9730000dge0404000","32c51079a","0000006z00000000","0400503000000000"); // 6: alpha < 0.8
    cuts.AddCutHeavyMesonPCM("80010113","0dm00009f9730000dge0404000","32c51079a","0000007z00000000","0400503000000000"); // 7: alpha < 0.85
  }else if (trainConfig == 1110){ // ITS cluster requirement variation
    cuts.AddCutHeavyMesonPCM("80010113","0dm00009f9730000dge0404000","34c51079a","0000003z00000000","0400503000000000"); // 4: min 3 ITS cluster
  }else if (trainConfig == 1111){ // TPC cluster requirement variation
    cuts.AddCutHeavyMesonPCM("80010113","0dm00009f9730000dge0404000","32e51079a","0000003z00000000","0400503000000000"); // e: No shared clusters
  }else if (trainConfig == 1112){ // Charged pion DCA
    cuts.AddCutHeavyMesonPCM("80010113","0dm00009f9730000dge0404000","32c01079a","0000003z00000000","0400503000000000"); // 0: No cut
    cuts.AddCutHeavyMesonPCM("80010113","0dm00009f9730000dge0404000","32c31079a","0000003z00000000","0400503000000000"); // 5: strict pT dep
    cuts.AddCutHeavyMesonPCM("80010113","0dm00009f9730000dge0404000","32c61079a","0000003z00000000","0400503000000000"); // 6: z,xy < 0.5 (very tight)
  }else if (trainConfig == 1113){ // TOF requirement
    cuts.AddCutHeavyMesonPCM("80010113","0dm00009f9730000dge0404000","32c51070a","0000003z00000000","0400503000000000"); // 0: No TOF
    cuts.AddCutHeavyMesonPCM("80010113","0dm00009f9730000dge0404000","32c51076a","0000003z00000000","0400503000000000"); // 6: stricter Kp rejection
    cuts.AddCutHeavyMesonPCM("80010113","0dm00009f9730000dge0404000","32c51072a","0000003z00000000","0400503000000000"); // 2: Pion selection
  }else if (trainConfig == 1114){ // min pT
    cuts.AddCutHeavyMesonPCM("80010113","0dm00009f9730000dge0404000","32c50079a","0000003z00000000","0400503000000000"); // 0: pT > 0.075 GeV
    cuts.AddCutHeavyMesonPCM("80010113","0dm00009f9730000dge0404000","32c52079a","0000003z00000000","0400503000000000"); // 2: pT > 0.125 GeV
    cuts.AddCutHeavyMesonPCM("80010113","0dm00009f9730000dge0404000","32c53079a","0000003z00000000","0400503000000000"); // 3: pT > 0.15 GeV
  }else if (trainConfig == 1115){ // TPC dEdx sigma
    cuts.AddCutHeavyMesonPCM("80010113","0dm00009f9730000dge0404000","32c510b9a","0000003z00000000","0400503000000000"); // b: -2,2
    cuts.AddCutHeavyMesonPCM("80010113","0dm00009f9730000dge0404000","32c51099a","0000003z00000000","0400503000000000"); // 9: -2.5,2.5
    cuts.AddCutHeavyMesonPCM("80010113","0dm00009f9730000dge0404000","32c510a9a","0000003z00000000","0400503000000000"); // a: -3.5,3.5
    cuts.AddCutHeavyMesonPCM("80010113","0dm00009f9730000dge0404000","32c51059a","0000003z00000000","0400503000000000"); // 5: -4,4
    cuts.AddCutHeavyMesonPCM("80010113","0dm00009f9730000dge0404000","32c51039a","0000003z00000000","0400503000000000"); // 3: -5,5
  }else if (trainConfig == 1116){ // PiPlPiMi Mass
    cuts.AddCutHeavyMesonPCM("80010113","0dm00009f9730000dge0404000","32c51079r","0000003z00000000","0400503000000000"); // r: 0.8 GeV/c^2
    cuts.AddCutHeavyMesonPCM("80010113","0dm00009f9730000dge0404000","32c51079s","0000003z00000000","0400503000000000"); // s: 0.825 GeV/c^2
    cuts.AddCutHeavyMesonPCM("80010113","0dm00009f9730000dge0404000","32c51079t","0000003z00000000","0400503000000000"); // t: 0.875 GeV/c^2
    cuts.AddCutHeavyMesonPCM("80010113","0dm00009f9730000dge0404000","32c51079u","0000003z00000000","0400503000000000"); // u: 0.9 GeV/c^2


 //************************************************ PCM- PHOS analysis 5 TeV pPb ********************************************
  } else if (trainConfig == 1501){ // PHOS  PHI7 run1
    cuts.AddCutHeavyMesonPCM("80062113","00200009f9730000dge0400000","32c51070a","0103603s00000000","0453503000000000");  // 0-100%


  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //                                          ETA PRIME MESON
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  //************************************************ PCM- EDC analysis 5 TeV pPb *********************************************

  // no event mixing background
  }else if (trainConfig == 2001){ // EMC  INT7 run1 & run2
    cuts.AddCutHeavyMesonPCM("80010113","00200009f9730000dge0400000","32c510700","0103603l00000000","0453503000000000"); // 0-100% without NL
  } else if (trainConfig == 2002){ // EMC EMC triggers
    cuts.AddCutHeavyMesonPCM("80052113","00200009f9730000dge0400000","32c510700","0103603l00000000","0453503000000000"); // 0-100% without NL just EMC, EMC7
  } else if (trainConfig == 2003){ // EMC EMC triggers
    cuts.AddCutHeavyMesonPCM("80083113","00200009f9730000dge0400000","32c510700","0103603l00000000","0453503000000000"); // 0-100% without NL just EMC, EG1
  } else if (trainConfig == 2004){ // EMC EMC triggers
    cuts.AddCutHeavyMesonPCM("80085113","00200009f9730000dge0400000","32c510700","0103603l00000000","0453503000000000"); // 0-100% without NL just EMC, EG2

    // EMCal + EDC triggers
  } else if (trainConfig == 2005){ // EMC EMC triggers
    cuts.AddCutHeavyMesonPCM("8008d113","00200009f9730000dge0400000","32c510700","0103603l00000000","0453503000000000"); // 0-100% without NL just EMC, EG1
  } else if (trainConfig == 2006){ // EMC EMC triggers
    cuts.AddCutHeavyMesonPCM("8008e113","00200009f9730000dge0400000","32c510700","0103603l00000000","0453503000000000"); // 0-100% without NL just EMC, EG2

  //************************************************ PCM- PHOS analysis 5 TeV pPb ********************************************
  } else if (trainConfig == 2501){ // PHOS  PHI7 run1
    cuts.AddCutHeavyMesonPCM("80062113","00200009f9730000dge0400000","32c510700","0103603l00000000","0453503000000000"); // 0-100% without NL
 
  //************************************************ PCM- EDC analysis 13 TeV pPb *********************************************
  } else if (trainConfig == 2010) { // min bias only
    cuts.AddCutHeavyMesonPCM("80010113","0dm00009f9730000dge0404000","32c51070m","0103603l00000000","0453503000000000"); // INT7
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

  if (periodNameV0Reader.Contains("LHC17g6a2") || periodNameV0Reader.Contains("LHC17g6a3") ){
    TObjString *HeaderPMB = new TObjString("Dpmjet_0");
    TObjString *HeaderP8J = new TObjString("Pythia8JetsGammaTrg_1");
    if (doWeighting==4) { // all headers
      HeaderList->Add(HeaderPMB);
      HeaderList->Add(HeaderP8J);
    } else if (doWeighting==5) { // only MB header
      HeaderList->Add(HeaderPMB);
    } else { // only JJ header
      HeaderList->Add(HeaderP8J);
    }
  }

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
    analysisEventCuts[i]->SetTriggerOverlapRejecion(enableTriggerOverlapRej);
    if(runLightOutput>0) analysisEventCuts[i]->SetLightOutput(kTRUE);
    analysisEventCuts[i]->InitializeCutsFromCutString((cuts.GetEventCut(i)).Data());
    if (periodNameV0Reader.CompareTo("") != 0) analysisEventCuts[i]->SetPeriodEnum(periodNameV0Reader);
    EventCutList->Add(analysisEventCuts[i]);
    analysisEventCuts[i]->SetFillCutHistograms("",kFALSE);

    analysisCuts[i] = new AliConversionPhotonCuts();
    if(runLightOutput>0) analysisCuts[i]->SetLightOutput(kTRUE);
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

  task->SetUnsmearedOutputs(unsmearingoutputs);

  //connect containers
  AliAnalysisDataContainer *coutput =
  mgr->CreateContainer(Form("GammaConvNeutralMesonPiPlPiMiNeutralMeson_%i_%i_%i.root",selectHeavyNeutralMeson,neutralPionMode, trainConfig), TList::Class(),
              AliAnalysisManager::kOutputContainer,Form("GammaConvNeutralMesonPiPlPiMiNeutralMeson_%i_%i_%i.root",selectHeavyNeutralMeson,neutralPionMode, trainConfig));

  mgr->AddTask(task);
  mgr->ConnectInput(task,0,cinput);
  mgr->ConnectOutput(task,1,coutput);

  return;

}
