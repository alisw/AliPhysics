/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: Friederike Bock, Lucia Leardini,                               *
 * Version 1.0                                                            *
 *                                                                        *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
AliAnalysisTaskGammaConvV1* AddTask_hikari(
   Int_t   trainConfig            = 1,
   Int_t   isMC                   = 0,    
   TString photonCutNumberV0Reader= "00000008400100001500000000", // 00000008400000000100000000 nom. B, 00000088400000000100000000 low B
   TString periodNameV0Reader     = "",
   Int_t   enableQAMesonTask      = 0,    
   Int_t   enableQAPhotonTask     = 0,
   TString fileNameExternalInputs = "",//FPTW:fileNamePtWeights, FMUW:fileNameMultWeights, FMAW:fileNameMatBudWeights, FEPC:fileNamedEdxPostCalib, separate with ;
   Int_t   enableMatBudWeightsPi0 = 0,        // 1 = three radial bins, 2 = 10 radial bins
   TString periodNameAnchor       = "",//Nch weighting & Ptweighting
   Bool_t  doPartWeight           = kFALSE,//partcle weighting
   Bool_t  doMultWeight           = kFALSE,//Nch weighting
   Bool_t  doPostCalibration      = kFALSE,//post calib
   Int_t   isHeavyIon             = 0,
   TString additionalTrainConfig  = "0"// additional counter for trainconfig + special settings
   ) {


  
  Bool_t    enableLightOutput = kFALSE;
  Bool_t    enableTHnSparse   = kFALSE;
  // TString   generatorName                 = "DPMJET", // generator Name
  // Bool_t    enableTriggerMimicking        = kFALSE,   // enable trigger mimicking
  // Bool_t    enableTriggerOverlapRej       = kFALSE,   // enable trigger overlap rejection
  // Int_t     debugLevel                    = 0,        // introducing debug levels for grid running
  // Float_t   maxFacPtHard                  = 3.,       // maximum factor between hardest jet and ptHard generated
  // // special settings
  // Bool_t    enableChargedPrimary          = kFALSE,
  // Bool_t    doSmear                       = kFALSE,   // switches to run user defined smearing
  // Double_t  bremSmear                     = 1.,
  // Double_t  smearPar                      = 0.,       // conv photon smearing params
  // Double_t  smearParConst                 = 0.,       // conv photon smearing params

  AliCutHandlerPCM cuts;

  TString fileNamePtWeights       = cuts.GetSpecialFileNameFromString (fileNameExternalInputs, "FPTW:");
  TString fileNameMultWeights     = cuts.GetSpecialFileNameFromString (fileNameExternalInputs, "FMUW:");
  TString fileNameMatBudWeights   = cuts.GetSpecialFileNameFromString (fileNameExternalInputs, "FMAW:");
  TString fileNamedEdxPostCalib   = cuts.GetSpecialFileNameFromString (fileNameExternalInputs, "FEPC:");

  TString addTaskName             = "AddTask_hikari";
  TString sAdditionalTrainConfig  = cuts.GetSpecialSettingFromAddConfig(additionalTrainConfig, "", "", addTaskName);
  if (sAdditionalTrainConfig.Atoi() > 0){
    trainConfig = trainConfig + sAdditionalTrainConfig.Atoi();
    cout << "INFO: " << addTaskName.Data() << " running additionalTrainConfig '" << sAdditionalTrainConfig.Atoi() << "', train config: '" << trainConfig << "'" << endl;
  }
  if(additionalTrainConfig.Contains("MaterialBudgetWeights"))
    fileNameMatBudWeights         = cuts.GetSpecialSettingFromAddConfig(additionalTrainConfig, "MaterialBudgetWeights",fileNameMatBudWeights, addTaskName);


  // ================== GetAnalysisManager ===============================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error(Form("%s_%i", addTaskName.Data(),  trainConfig), "No analysis manager found.");
    return NULL;
  }

  // ================== GetInputEventHandler =============================
  AliVEventHandler *inputHandler=mgr->GetInputEventHandler();

  //========= Check whether PID Reponse is there ====
  if(!(AliPIDResponse*)mgr->GetTask("PIDResponseTask")){
    Error(Form("AddTask_GammaConvV1_%i",trainConfig), "No PID response has been initialized aborting.");
    return NULL;
  }

  //=========  Set Cutnumber for V0Reader ================================
  TString cutnumberPhoton     = photonCutNumberV0Reader.Data();
  TString cutnumberEvent      = "00000003";
  if(isHeavyIon==2){
    cutnumberEvent      = "80000003";
  }
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  
  //========= Check V0 Reader in  ANALYSIS manager  =====
  TString V0ReaderName        = Form("V0ReaderV1_%s_%s",cutnumberEvent.Data(),cutnumberPhoton.Data());
  AliV0ReaderV1 *fV0ReaderV1  =  NULL;
  if( !(AliV0ReaderV1*)mgr->GetTask(V0ReaderName.Data()) ){
    cout << "V0Reader: " << V0ReaderName.Data() << " not found!!"<< endl;
    return NULL;
  } else {
    cout << "V0Reader: " << V0ReaderName.Data() << " found!!"<< endl;
  }

  //========= Add task to the ANALYSIS manager find input container
  AliAnalysisTaskGammaConvV1 *task=NULL;
  task= new AliAnalysisTaskGammaConvV1(Form("GammaConvV1_%i",trainConfig));
  task->SetIsHeavyIon(isHeavyIon);
  task->SetIsMC(isMC);
  task->SetV0ReaderName(V0ReaderName);
  task->SetLightOutput(enableLightOutput);
  if(isHeavyIon==0){
    if (trainConfig == 1){
      cuts.AddCutPCM("00010113", "00200009227300008250404000", "0152101500000000"); // alpha cut on 
      cuts.AddCutPCM("00010113", "00200009227300008250404000", "0152103500000000"); // alpha cut off   //A
      cuts.AddCutPCM("00010113", "00200009a27300008250404000", "0152103500000000"); // dEdx            //B
      cuts.AddCutPCM("00010113", "00200009227300008250a04000", "0152103500000000"); // pa 
      cuts.AddCutPCM("00010113", "00200009227300008250404120", "0152103500000000"); // dca
    } else if (trainConfig == 2){
    //cuts.AddCutPCM("00010113", "00200009227300008250404000", "0152103500000000"); // eta 0.9         //A
      cuts.AddCutPCM("00010113", "0c200009227300008250404000", "0152103500000000"); // eta 0.85
      cuts.AddCutPCM("00010113", "0d200009227300008250404000", "0152103500000000"); // eta 0.8         //D
      cuts.AddCutPCM("00010113", "04200009227300008250404000", "0152103500000000"); // eta 0.75    
    } else if (trainConfig == 3){
    //cuts.AddCutPCM("00010113", "00200009a27300008250404000", "0152103500000000"); // eta 0.9         //B
      cuts.AddCutPCM("00010113", "0c200009a27300008250404000", "0152103500000000"); // eta 0.85
      cuts.AddCutPCM("00010113", "0d200009a27300008250404000", "0152103500000000"); // eta 0.8         //C
      cuts.AddCutPCM("00010113", "04200009a27300008250404000", "0152103500000000"); // eta 0.75    
    } else if (trainConfig == 4){
    //cuts.AddCutPCM("00010113", "0d200009227300008250404000", "0152103500000000"); // pa 0.85         //D
      cuts.AddCutPCM("00010113", "0d200009227300008250b04000", "0152103500000000"); // pa 0.985
      cuts.AddCutPCM("00010113", "0d200009227300008250904000", "0152103500000000"); // pa 0.99
      cuts.AddCutPCM("00010113", "0d200009227300008250a04000", "0152103500000000"); // pa 0.995
    } else if (trainConfig == 5){
    //cuts.AddCutPCM("00010113", "0d200009a27300008250404000", "0152103500000000"); // pa 0.85         //C
      cuts.AddCutPCM("00010113", "0d200009a27300008250b04000", "0152103500000000"); // pa 0.985
      cuts.AddCutPCM("00010113", "0d200009a27300008250904000", "0152103500000000"); // pa 0.99
      cuts.AddCutPCM("00010113", "0d200009a27300008250a04000", "0152103500000000"); // pa 0.995        //E
    } else if (trainConfig == 6){
      cuts.AddCutPCM("00010113", "0d200009a27300008250a04000", "0152103500000000"); // AllCat          //E
      cuts.AddCutPCM("00010113", "0d200009a27300008250a24000", "0152103500000000"); // Cat1   TPC
      cuts.AddCutPCM("00010113", "0d200009a27300008250a54000", "0152103500000000"); // Cat2,3 ITS
    } else if (trainConfig == 7){
      cuts.AddCutPCM("00010113", "0d200009a27300008250a04000", "0152103500000000"); // AllCat          //E
      cuts.AddCutPCM("00010113", "0d200009a27300008250a24000", "0152103500000000"); // Cat1   TPC
      cuts.AddCutPCM("00010113", "0d200009a27300008250a54000", "0152103500000000"); // Cat2,3 ITS
    } else if (trainConfig == 8){
    //cuts.AddCutPCM("00010113", "0d200009227300008250a04000", "0152103500000000"); // standard cut for pp 5 TeV analysis VAND //E
      cuts.AddCutPCM("00010113", "0d200008227300008250a04000", "0152103500000000"); // TPC cluster 35%
      cuts.AddCutPCM("00010113", "0d200006227300008250a04000", "0152103500000000"); // TPC cluster 70%
    } else if (trainConfig == 9){
      cuts.AddCutPCM("00010113", "0d200079227300008250a04000", "0152103500000000"); // min pT no cut
      cuts.AddCutPCM("00010113", "0d200069227300008250a04000", "0152103500000000"); // min pT 40 MeV
      cuts.AddCutPCM("00010113", "0d200049227300008250a04000", "0152103500000000"); // min pT 75 MeV
    } else if (trainConfig == 10){
      cuts.AddCutPCM("00010113", "0d200009a27300008250a04000", "0152103500000000"); // standard cut for pp 5 TeV analysis VAND
      cuts.AddCutPCM("00010113", "0d200008a27300008250a04000", "0152103500000000"); // TPC cluster 35%
      cuts.AddCutPCM("00010113", "0d200006a27300008250a04000", "0152103500000000"); // TPC cluster 70%
    } else if (trainConfig == 11){
      cuts.AddCutPCM("00010113", "0d200079a27300008250a04000", "0152103500000000"); // min pT no cut
      cuts.AddCutPCM("00010113", "0d200069a27300008250a04000", "0152103500000000"); // min pT 40 MeV
      cuts.AddCutPCM("00010113", "0d200049a27300008250a04000", "0152103500000000"); // min pT 75 MeV    
    } else if (trainConfig == 12){
      cuts.AddCutPCM("00010113", "0d200009a27300008250a04000", "0152103500000000"); // edEdx -3,3
      cuts.AddCutPCM("00010113", "0d200009b27300008250a04000", "0152103500000000"); // edEdx -3,2,3.2
      cuts.AddCutPCM("00010113", "0d200009c27300008250a04000", "0152103500000000"); // edEdx -2.8,2.8
      cuts.AddCutPCM("00010113", "0d200009a57300008250a04000", "0152103500000000"); // pidEdx 2,-10
      cuts.AddCutPCM("00010113", "0d200009a17300008250a04000", "0152103500000000"); // pidEdx 0,-10
    } else if (trainConfig == 13){
      cuts.AddCutPCM("00010113", "0d200009a27300008250a04000", "0152103500000000"); // edEdx -3,3
      cuts.AddCutPCM("00010113", "0d200009b27300008250a04000", "0152103500000000"); // edEdx -3,2,3.2
      cuts.AddCutPCM("00010113", "0d200009c27300008250a04000", "0152103500000000"); // edEdx -2.8,2.8
      cuts.AddCutPCM("00010113", "0d200009a57300008250a04000", "0152103500000000"); // pidEdx 2,-10
      cuts.AddCutPCM("00010113", "0d200009a17300008250a04000", "0152103500000000"); // pidEdx 0,-10
    } else if (trainConfig == 14){
      cuts.AddCutPCM("00010113", "0d200009220300008250a04000", "0152103500000000"); // pion nsig min mom 0.50 GeV/c
      cuts.AddCutPCM("00010113", "0d200009226300008250a04000", "0152103500000000"); // pion nsig min mom 0.25 GeV/c
      cuts.AddCutPCM("00010113", "0d200009227600008250a04000", "0152103500000000"); // pion nsig max mom 2.00 GeV/c
      cuts.AddCutPCM("00010113", "0d200009227100008250a04000", "0152103500000000"); // pion nsig max mom 5.00 GeV/c
    } else if (trainConfig == 15){
      cuts.AddCutPCM("00010113", "0d200009a20300008250a04000", "0152103500000000"); // pion nsig min mom 0.50 GeV/c
      cuts.AddCutPCM("00010113", "0d200009a26300008250a04000", "0152103500000000"); // pion nsig min mom 0.25 GeV/c
      cuts.AddCutPCM("00010113", "0d200009a27600008250a04000", "0152103500000000"); // pion nsig max mom 2.00 GeV/c
      cuts.AddCutPCM("00010113", "0d200009a27100008250a04000", "0152103500000000"); // pion nsig max mom 5.00 GeV/c
    } else if (trainConfig == 16){
      cuts.AddCutPCM("00010113", "0d200009227300003250a04000", "0152103500000000"); // qT max 0.05 1D
      cuts.AddCutPCM("00010113", "0d200009227300002250a04000", "0152103500000000"); // qT max 0.06 2D
      cuts.AddCutPCM("00010113", "0d200009227300008210a04000", "0152103500000000"); // Psi pair 0.1  1D
      cuts.AddCutPCM("00010113", "0d200009227300008260a04000", "0152103500000000"); // Psi pair 0.05 2D
      cuts.AddCutPCM("00010113", "0d020009227300008280a04000", "0152103500000000"); // Psi pair 0.2  2D
    } else if (trainConfig == 17){
      cuts.AddCutPCM("00010113", "0d200009a27300003250a04000", "0152103500000000"); // qT max 0.05 1D
      cuts.AddCutPCM("00010113", "0d200009a27300002250a04000", "0152103500000000"); // qT max 0.06 2D
      cuts.AddCutPCM("00010113", "0d200009a27300008210a04000", "0152103500000000"); // Psi pair 0.1  1D
      cuts.AddCutPCM("00010113", "0d200009a27300008260a04000", "0152103500000000"); // Psi pair 0.05 2D
      cuts.AddCutPCM("00010113", "0d020009a27300008280a04000", "0152103500000000"); // Psi pair 0.2  2D
    } else if (trainConfig == 18){
      cuts.AddCutPCM("00010113", "0d200009227300008210a04000", "0152103500000000"); // variation chi2 30 psi pair 0.1 1D
      cuts.AddCutPCM("00010113", "0d200009227300008860a04000", "0152103500000000"); // variation chi2 20 psi pair 0.05 2D
      cuts.AddCutPCM("00010113", "0d200009227300008180a04000", "0152103500000000"); // variation chi2 50 psi pair 0.2 2D
      cuts.AddCutPCM("00010113", "0d200009227300008254a04000", "0152103500000000"); // Photon Asymmetry Cut
    } else if (trainConfig == 19){
      cuts.AddCutPCM("00010113", "0d200009a27300008210a04000", "0152103500000000"); // variation chi2 30 psi pair 0.1 1D
      cuts.AddCutPCM("00010113", "0d200009a27300008860a04000", "0152103500000000"); // variation chi2 20 psi pair 0.05 2D
      cuts.AddCutPCM("00010113", "0d200009a27300008180a04000", "0152103500000000"); // variation chi2 50 psi pair 0.2 2D
      cuts.AddCutPCM("00010113", "0d200009a27300008254a04000", "0152103500000000"); // Photon Asymmetry Cut
    } else if (trainConfig == 20){
      cuts.AddCutPCM("00010113", "0d200009a27300008250a00000", "0152103500000000"); // no double counting
      cuts.AddCutPCM("00010113", "0d200009a27300008250a04000", "0152105500000000"); // meson alpha < 0.75
      cuts.AddCutPCM("00010113", "0d200009a27300008250a04000", "0152107500000000"); // meson alpha < 0.85
    } else if (trainConfig == 21){
      cuts.AddCutPCM("00010113", "0d200009227300008250a04030", "0152103500000000"); // dcaz 4cm
      cuts.AddCutPCM("00010113", "0d200009227300008250a04040", "0152103500000000"); // dcaz 3cm
      cuts.AddCutPCM("00010113", "0d200009227300008250a04200", "0152103500000000"); // dcar 5cm
      cuts.AddCutPCM("00010113", "0d200009227300008250a04300", "0152105500000000"); // dcar 4cm
    } else if (trainConfig == 22){
      cuts.AddCutPCM("00010113", "0d200009a27300008250a04030", "0152103500000000"); // dcaz 4cm
      cuts.AddCutPCM("00010113", "0d200009a27300008250a04040", "0152103500000000"); // dcaz 3cm
      cuts.AddCutPCM("00010113", "0d200009a27300008250a04200", "0152103500000000"); // dcar 5cm
      cuts.AddCutPCM("00010113", "0d200009a27300008250a04300", "0152105500000000"); // dcar 4cm

    } else if (trainConfig == 100){//cut selection for pp 5 TeV 2017 ------------------------------------
      cuts.AddCutPCM("00010113", "0d200009227300008250404000", "0152103500000000"); // Standard cut for pp 5 TeV analysis VAND
      cuts.AddCutPCM("00010113", "0d200008227300008250404000", "0152103500000000"); // TPC cluster 35%
      cuts.AddCutPCM("00010113", "0d200006227300008250404000", "0152103500000000"); // TPC cluster 70%
    } else if (trainConfig == 101){
      cuts.AddCutPCM("00010113", "0d200079227300008250404000", "0152103500000000"); // min pT no cut
      cuts.AddCutPCM("00010113", "0d200049227300008250404000", "0152103500000000"); // min pT 75 cut
      cuts.AddCutPCM("00010113", "0d200019227300008250404000", "0152103500000000"); // min pT 100 cut
      cuts.AddCutPCM("00010113", "0d200029227300008250404000", "0152103500000000"); // min pT 150 cut
    } else if (trainConfig == 102){
      cuts.AddCutPCM("00010113", "0d200009327300008250404000", "0152103500000000"); // edEdx -4,5
      cuts.AddCutPCM("00010113", "0d200009627300008250404000", "0152103500000000"); // edEdx -2.5,4
      cuts.AddCutPCM("00010113", "0d200009257300008250404000", "0152103500000000"); // pidEdx 2,-10
      cuts.AddCutPCM("00010113", "0d200009217300008250404000", "0152103500000000"); // pidEdx 0,-10
    } else if (trainConfig == 103){
      cuts.AddCutPCM("00010113", "0d200009247300008250404000", "0152103500000000"); // pidEdx 3, 1
      cuts.AddCutPCM("00010113", "0d200009226300008250404000", "0152103500000000"); // pion nsig min mom 0.25 GeV/c
      cuts.AddCutPCM("00010113", "0d200009227600008250404000", "0152103500000000"); // pion nsig max mom 2.00 GeV/c
      cuts.AddCutPCM("00010113", "0d200009227100008250404000", "0152103500000000"); // pion nsig max mom 5.00 GeV/c
    } else if (trainConfig == 104){
      cuts.AddCutPCM("00010113", "0d200009227300002250404000", "0152103500000000"); // qT max 0.06 2D
      cuts.AddCutPCM("00010113", "0d200009227300009250404000", "0152103500000000"); // qT max 0.03 2D
      cuts.AddCutPCM("00010113", "0d200009227300008260404000", "0152103500000000"); // Psi pair 0.05 2D, chi2 30.
      cuts.AddCutPCM("00010113", "0d200009227300008280404000", "0152103500000000"); // Psi pair 0.2  2D, chi2 30.
    } else if (trainConfig == 105){
      cuts.AddCutPCM("00010113", "0d200009227300008850404000", "0152103500000000"); // chi2 20. psi pair 0.1 2D
      cuts.AddCutPCM("00010113", "0d200009227300008150404000", "0152103500000000"); // chi2 50. psi pair 0.1 2D
      cuts.AddCutPCM("00010113", "0d200009227300008254404000", "0152103500000000"); // Photon Asymmetry Cut
      cuts.AddCutPCM("00010113", "0d200009227300008250604000", "0152103500000000"); // CosPA 0.9
    } else if (trainConfig == 106){
      cuts.AddCutPCM("00010113", "0d200009227300008860404000", "0152103500000000"); // variation chi2 20 psi pair 0.2 2D
      cuts.AddCutPCM("00010113", "0d200009227300002252404000", "0152103500000000"); // variation qT max 0.06 2D, asym var 1
      cuts.AddCutPCM("00010113", "0d200009227300002254404000", "0152103500000000"); // variation qT max 0.06 2D, asym vat 2 pt dep
      cuts.AddCutPCM("00010113", "0d200009227300002256404000", "0152103500000000"); // variation qT max 0.06 2D, asym var 3
    } else if (trainConfig == 107){
      cuts.AddCutPCM("00010113", "0d200009227300008250004000", "0152103500000000"); // no CosPA
      cuts.AddCutPCM("00010113", "0d200009227300008250400000", "0152103500000000"); // no double counting
      cuts.AddCutPCM("00010113", "0d200009227300008250404000", "0152101500000000"); // meson alpha pt dep
      cuts.AddCutPCM("00010113", "0d200009227300008250404000", "0152107500000000"); // meson alpha < 0.85

    } else if (trainConfig == 400){//cut selection for pp 5 TeV 2017 ------------------------------------
      cuts.AddCutPCM("00010113", "0d200009227300008250404000", "0152103500000000"); // Standard cut for pp 5 TeV analysis VAND
    } else if (trainConfig == 401){
      cuts.AddCutPCM("00010113", "0d200009227300008250424000", "0152103500000000"); // Standard cut for pp 5 TeV analysis VAND Cat1
    } else if (trainConfig == 402){
      cuts.AddCutPCM("00010113", "0d200009227300008250454000", "0152103500000000"); // Standard cut for pp 5 TeV analysis VAND Cat23
    } else if (trainConfig == 403){
      cuts.AddCutPCM("00010113", "00200009227300008250404000", "0152103500000000"); // eta < 0.9
      cuts.AddCutPCM("00010113", "0c200009227300008250404000", "0152103500000000"); // eta < 0.85
      cuts.AddCutPCM("00010113", "0d200009247300008250404000", "0152103500000000"); // pidEdx 3, 1
      cuts.AddCutPCM("00010113", "0d200009227100008250404000", "0152103500000000"); // pion nsig max mom 5.00 GeV/c
    } else if (trainConfig == 404){
      cuts.AddCutPCM("00010113", "0d200009227300008860404000", "0152103500000000"); // variation chi2 20 psi pair 0.2 2D
      cuts.AddCutPCM("00010113", "0d200009227300002252404000", "0152103500000000"); // variation qT max 0.06 2D, asym var 1
      cuts.AddCutPCM("00010113", "0d200009227300002254404000", "0152103500000000"); // variation qT max 0.06 2D, asym vat 2 pt dep
      cuts.AddCutPCM("00010113", "0d200009227300002256404000", "0152103500000000"); // variation qT max 0.06 2D, asym var 3
    } else if (trainConfig == 405){
      cuts.AddCutPCM("00010113", "0d200009a27300008250904120", "0152103500000000"); //cosPA, 0.99 eta 0.9
      cuts.AddCutPCM("00010113", "0d200079227300008250404000", "0152103500000000"); // min pT no cut
      cuts.AddCutPCM("00010113", "0d200049227300008250404000", "0152103500000000"); // min pT 75 MeV
      cuts.AddCutPCM("00010113", "0d200019227300008250404000", "0152103500000000"); // min pT 100 MeV
      cuts.AddCutPCM("00010113", "0d200008227300008250404000", "0152103500000000"); // TPC cluster 35%
      cuts.AddCutPCM("00010113", "0d200006227300008250404000", "0152103500000000"); // TPC cluster 70%
    } else if (trainConfig == 406){
      cuts.AddCutPCM("00010113", "0d200009327300008250404000", "0152103500000000"); // edEdx -4,5
      cuts.AddCutPCM("00010113", "0d200009627300008250404000", "0152103500000000"); // edEdx -2.5,4
      cuts.AddCutPCM("00010113", "0d200009257300008250404000", "0152103500000000"); // pidEdx 2,-10
      cuts.AddCutPCM("00010113", "0d200009217300008250404000", "0152103500000000"); // pidEdx 0,-10
      cuts.AddCutPCM("00010113", "0d200009226300008250404000", "0152103500000000"); // pion nsig min mom 0.25 GeV/c
      cuts.AddCutPCM("00010113", "0d200009227600008250404000", "0152103500000000"); // pion nsig max mom 2.00 GeV/c
    } else if (trainConfig == 407){
      cuts.AddCutPCM("00010113", "0d200009227300002250404000", "0152103500000000"); // qT max 0.06 2D
      cuts.AddCutPCM("00010113", "0d200009227300009250404000", "0152103500000000"); // qT max 0.03 2D
      cuts.AddCutPCM("00010113", "0d200009227300008260404000", "0152103500000000"); // Psi pair 0.05 2D, chi2 30.
      cuts.AddCutPCM("00010113", "0d200009227300008280404000", "0152103500000000"); // Psi pair 0.2  2D, chi2 30.
      cuts.AddCutPCM("00010113", "0d200009227300008850404000", "0152103500000000"); // chi2 20. psi pair 0.1 2D
      cuts.AddCutPCM("00010113", "0d200009227300008150404000", "0152103500000000"); // chi2 50. psi pair 0.1 2D
    } else if (trainConfig == 408){
      cuts.AddCutPCM("00010113", "0d200009227300008254404000", "0152103500000000"); // Photon Asymmetry Cut
      cuts.AddCutPCM("00010113", "0d200009227300008250604000", "0152103500000000"); // CosPA 0.9
      cuts.AddCutPCM("00010113", "0d200009227300008250004000", "0152103500000000"); // no CosPA
      cuts.AddCutPCM("00010113", "0d200009227300008250400000", "0152103500000000"); // no double counting
      cuts.AddCutPCM("00010113", "0d200009227300008250404000", "0152101500000000"); // meson alpha pt dep
      cuts.AddCutPCM("00010113", "0d200009227300008250404000", "0152107500000000"); // meson alpha < 0.85

    } else if (trainConfig == 409){
      cuts.AddCutPCM("00010113", "0d200009247000008250404000", "0152103500000000"); //
      cuts.AddCutPCM("00010113", "0d200009287000008250404000", "0152103500000000"); //
      cuts.AddCutPCM("00010113", "0d200009297000008250404000", "0152103500000000"); //
    } else if (trainConfig == 410){
      cuts.AddCutPCM("00010113", "0d200009247000008250404000", "0152103500000000"); // to be used for MBW
      cuts.AddCutPCM("00010113", "0d200009287000008250404000", "0152103500000000"); // to be used for MBW
      cuts.AddCutPCM("00010113", "0d200009297000008250404000", "0152103500000000"); // to be used for MBW

    } else if (trainConfig == 440){ // as 400 to be used MBW eta
      cuts.AddCutPCM("00010113", "00200009227300008250404000", "0152103500000000"); // 
      cuts.AddCutPCM("00010113", "0c200009227300008250404000", "0152103500000000"); // 
      cuts.AddCutPCM("00010113", "0d200009227300008250404000", "0152103500000000"); // 
    } else if (trainConfig == 441){// as 440 to be used MBW
      cuts.AddCutPCM("00010113", "0da00009227300008250404000", "0152103500000000"); // R 5-33.5 cm
      cuts.AddCutPCM("00010113", "0db00009227300008250404000", "0152103500000000"); // R 33.5-72 cm
      cuts.AddCutPCM("00010113", "0dc00009227300008250404000", "0152103500000000"); // R 72-180 cm
    } else if (trainConfig == 442){// small R
      cuts.AddCutPCM("00010113", "0dh00009227300008250404000", "0152103500000000"); // R 5-13
      cuts.AddCutPCM("00010113", "0di00009227300008250404000", "0152103500000000"); // R 13-33.5.
      cuts.AddCutPCM("00010113", "0dj00009227300008250404000", "0152103500000000"); // R 33-55
    } else if (trainConfig == 443){// large R
      cuts.AddCutPCM("00010113", "0dk00009227300008250404000", "0152103500000000"); // R 55-72
      cuts.AddCutPCM("00010113", "0dl00009227300008250404000", "0152103500000000"); // R 72-95
      cuts.AddCutPCM("00010113", "0dg00009227300008250404000", "0152103500000000"); // R 95-180
    } else if (trainConfig == 444){// Ana's request
      cuts.AddCutPCM("00010113", "0dm00009227300008250404000", "0152103500000000"); // R 5-180, exclude 55-72
    } else if (trainConfig == 445){// default + category
      cuts.AddCutPCM("00010113", "0d200009227300008250404000", "0152103500000000"); // AllCat          //E
      cuts.AddCutPCM("00010113", "0d200009227300008250424000", "0152103500000000"); // Cat1   TPC
      cuts.AddCutPCM("00010113", "0d200009227300008250454000", "0152103500000000"); // Cat2,3 ITS
    } else if (trainConfig == 446){//
      cuts.AddCutPCM("00010113", "0d200079227300008250404000", "0152103500000000"); // min pT no cut
      cuts.AddCutPCM("00010113", "0d200049227300008250404000", "0152103500000000"); // min pT 75 MeV
      cuts.AddCutPCM("00010113", "0d200019227300008250404000", "0152103500000000"); // min pT 100 MeV
      cuts.AddCutPCM("00010113", "0d200029227300008250404000", "0152103500000000"); // min pT 150 MeV
    } else if (trainConfig == 447){//
      cuts.AddCutPCM("00010113", "0d200008227300008250404000", "0152103500000000"); // TPC cluster 35%
      cuts.AddCutPCM("00010113", "0d200006227300008250404000", "0152103500000000"); // TPC cluster 70%
    } else if (trainConfig == 448){//
      cuts.AddCutPCM("00010113", "0d200009a27300008250404000", "0152103500000000"); // edEdx -3,3
      cuts.AddCutPCM("00010113", "0d200009b27300008250404000", "0152103500000000"); // edEdx -3,2,3.2
      cuts.AddCutPCM("00010113", "0d200009c27300008250404000", "0152103500000000"); // edEdx -2.8,2.8
    } else if (trainConfig == 449){//R 5-33.5 cm 
      cuts.AddCutPCM("00010113", "0da00009227300008250404000", "0152103500000000"); // vertex +-10
      cuts.AddCutPCM("00010114", "0da00009227300008250404000", "0152103500000000"); // vertex +-7.5
      cuts.AddCutPCM("00010115", "0da00009227300008250404000", "0152103500000000"); // vertex +-5.
      cuts.AddCutPCM("00010116", "0da00009227300008250404000", "0152103500000000"); // vertex +-2.5
    } else if (trainConfig == 450){//cat 1 Meson selection
      cuts.AddCutPCM("00010113", "0d200009227300008250404000", "0152103510000000"); // shared electron
      cuts.AddCutPCM("00010113", "0d200009227300008250404000", "0152103520000000"); // veto cat 1 Meson
    } else if (trainConfig == 451){//R max 180
      cuts.AddCutPCM("00010113", "0d200009227300008250404000", "0152103500000000"); // 5  < R < 180 (default)
      cuts.AddCutPCM("00010113", "0dn00009227300008250404000", "0152103500000000"); // 10 < R < 180
      cuts.AddCutPCM("00010113", "0do00009227300008250404000", "0152103500000000"); // 15 < R < 180
      cuts.AddCutPCM("00010113", "0dp00009227300008250404000", "0152103500000000"); // 20 < R < 180
    } else if (trainConfig == 452){//R max 95
      cuts.AddCutPCM("00010113", "0dq00009227300008250404000", "0152103500000000"); // 5  < R < 95
      cuts.AddCutPCM("00010113", "0dr00009227300008250404000", "0152103500000000"); // 10 < R < 95
      cuts.AddCutPCM("00010113", "0ds00009227300008250404000", "0152103500000000"); // 15 < R < 95
      cuts.AddCutPCM("00010113", "0dt00009227300008250404000", "0152103500000000"); // 20 < R < 95
    } else {
      Error(Form("GammaConvV1_%i",trainConfig), "wrong trainConfig variable no cuts have been specified for the configuration");
      return NULL;
    }

  }else if(isHeavyIon==2){
    if (trainConfig == 1){
      cuts.AddCutPCM("80010113", "0d200009a27300008250404000", "0162103500000000"); // default 0100 +
    } else if (trainConfig == 2){
      cuts.AddCutPCM("a0110113", "0d200009a27300008250404000", "0162103500000000"); // default 0005
    } else if (trainConfig == 3){
      cuts.AddCutPCM("80210113", "0d200009a27300008250404000", "0162103500000000"); // default 0020
    } else if (trainConfig == 4){
      cuts.AddCutPCM("80010113", "00200009a27300008250404000", "0162103500000000"); // eta 0.9
      cuts.AddCutPCM("80010113", "0c200009a27300008250404000", "0162103500000000"); // eta 0.85
    //cuts.AddCutPCM("80010113", "0d200009a27300008250404000", "0162103500000000"); // eta 0.8      +
      cuts.AddCutPCM("80010113", "04200009a27300008250404000", "0162103500000000"); // eta 0.75
    } else if (trainConfig == 5){
      cuts.AddCutPCM("80010113", "0d200009327300008250404000", "0162103500000000"); // dEdx -3 5
      cuts.AddCutPCM("80010113", "0d200009227300008250404000", "0162103500000000"); // dEdx -3 4
    //cuts.AddCutPCM("80010113", "0d200009a27300008250404000", "0162103500000000"); // dEdx -3 3
    } else if (trainConfig == 6){
    //cuts.AddCutPCM("80010113", "0d200009a27300008250404000", "0162103500000000"); // pa 0.85
      cuts.AddCutPCM("80010113", "0d200009a27300008250b04000", "0162103500000000"); // pa 0.985
      cuts.AddCutPCM("80010113", "0d200009a27300008250904000", "0162103500000000"); // pa 0.99
    } else if (trainConfig == 7){
      cuts.AddCutPCM("80010113", "0d200009a27300008250404000", "0162103500000000"); // All          +
      cuts.AddCutPCM("80010113", "0d200009a27300008250424000", "0162103500000000"); // Cat1
      cuts.AddCutPCM("80010113", "0d200009a27300008250454000", "0162103500000000"); // Cat23
    } else if (trainConfig == 8){
      cuts.AddCutPCM("80010113", "0d200079a27000008250404000", "0162103500000000"); // min pT no cut
      cuts.AddCutPCM("80010113", "0d200069a27000008250404000", "0162103500000000"); // min pT 40 cut
      cuts.AddCutPCM("80010113", "0d200049a27000008250404000", "0162103500000000"); // min pT 75 cut
    } else if (trainConfig == 9){
      cuts.AddCutPCM("80010113", "0d200008a27000008250404000", "0162103500000000"); // TPC cluster 35%
      cuts.AddCutPCM("80010113", "0d200007a27000008250404000", "0162103500000000"); // TPC cluster 70%
      cuts.AddCutPCM("80010113", "0d200009b27000008250404000", "0162103500000000"); // edEdx -3.2,3.2
      cuts.AddCutPCM("80010113", "0d200009c27000008250404000", "0162103500000000"); // edEdx -2.8,2.8
    } else if (trainConfig == 10){
      cuts.AddCutPCM("80010113", "0d200009a57000008250404000", "0162103500000000"); // pidEdx 2,-10
      cuts.AddCutPCM("80010113", "0d200009a17000008250404000", "0162103500000000"); // pidEdx 0,-10
      cuts.AddCutPCM("80010113", "0d200009a27000003250404000", "0162103500000000"); // qT max 0.05 1D
      cuts.AddCutPCM("80010113", "0d200009a27000002250404000", "0162103500000000"); // qT max 0.06 2D
    } else if (trainConfig == 11){
      cuts.AddCutPCM("a0110113", "0d200079a27000008250404000", "0162103500000000"); // min pT no cut
      cuts.AddCutPCM("a0110113", "0d200069a27000008250404000", "0162103500000000"); // min pT 40 cut
      cuts.AddCutPCM("a0110113", "0d200049a27000008250404000", "0162103500000000"); // min pT 75 cut
    } else if (trainConfig == 12){
      cuts.AddCutPCM("a0110113", "0d200008a27000008250404000", "0162103500000000"); // TPC cluster 35%
      cuts.AddCutPCM("a0110113", "0d200007a27000008250404000", "0162103500000000"); // TPC cluster 70%
      cuts.AddCutPCM("a0110113", "0d200009b27000008250404000", "0162103500000000"); // edEdx -3.2,3.2
      cuts.AddCutPCM("a0110113", "0d200009c27000008250404000", "0162103500000000"); // edEdx -2.8,2.8
    } else if (trainConfig == 13){
      cuts.AddCutPCM("a0110113", "0d200009a57000008250404000", "0162103500000000"); // pidEdx 2,-10
      cuts.AddCutPCM("a0110113", "0d200009a17000008250404000", "0162103500000000"); // pidEdx 0,-10
      cuts.AddCutPCM("a0110113", "0d200009a27000003250404000", "0162103500000000"); // qT max 0.05 1D
      cuts.AddCutPCM("a0110113", "0d200009a27000002250404000", "0162103500000000"); // qT max 0.06 2D
    } else if (trainConfig == 14){
      cuts.AddCutPCM("80210113", "0d200079a27000008250404000", "0162103500000000"); // min pT no cut
      cuts.AddCutPCM("80210113", "0d200069a27000008250404000", "0162103500000000"); // min pT 40 cut
      cuts.AddCutPCM("80210113", "0d200049a27000008250404000", "0162103500000000"); // min pT 75 cut
    } else if (trainConfig == 15){
      cuts.AddCutPCM("80210113", "0d200008a27000008250404000", "0162103500000000"); // TPC cluster 35%
      cuts.AddCutPCM("80210113", "0d200007a27000008250404000", "0162103500000000"); // TPC cluster 70%
      cuts.AddCutPCM("80210113", "0d200009b27000008250404000", "0162103500000000"); // edEdx -3.2,3.2
      cuts.AddCutPCM("80210113", "0d200009c27000008250404000", "0162103500000000"); // edEdx -2.8,2.8
    } else if (trainConfig == 16){
      cuts.AddCutPCM("80210113", "0d200009a57000008250404000", "0162103500000000"); // pidEdx 2,-10
      cuts.AddCutPCM("80210113", "0d200009a17000008250404000", "0162103500000000"); // pidEdx 0,-10
      cuts.AddCutPCM("80210113", "0d200009a27000003250404000", "0162103500000000"); // qT max 0.05 1D
      cuts.AddCutPCM("80210113", "0d200009a27000002250404000", "0162103500000000"); // qT max 0.06 2D
    } else if (trainConfig == 17){
      cuts.AddCutPCM("82410113", "0d200009a27300008250404000", "0162103500000000"); // default 2040 
    } else if (trainConfig == 18){
      cuts.AddCutPCM("84610113", "0d200009a27300008250404000", "0162103500000000"); // default 4060
    } else if (trainConfig == 19){
      cuts.AddCutPCM("86810113", "0d200009a27300008250404000", "0162103500000000"); // default 6080
    } else if (trainConfig == 20){
      cuts.AddCutPCM("88010113", "0d200009a27300008250404000", "0162103500000000"); // default 80100


      //MB
    } else if (trainConfig == 100){
      cuts.AddCutPCM("80010113", "0d200009a27300008250404000", "0162103500000000"); // default 0100 +
      cuts.AddCutPCM("80010113", "0d200008a27300008250404000", "0162103500000000"); // TPC cluster 35%
      cuts.AddCutPCM("80010113", "0d200006a27300008250404000", "0162103500000000"); // TPC cluster 70%
    } else if (trainConfig == 101){
      cuts.AddCutPCM("80010113", "0d200079a27300008250404000", "0162103500000000"); // min pT no cut
      cuts.AddCutPCM("80010113", "0d200069a27300008250404000", "0162103500000000"); // min pT 40 cut
      cuts.AddCutPCM("80010113", "0d200049a27300008250404000", "0162103500000000"); // min pT 75 cut
    } else if (trainConfig == 102){
      cuts.AddCutPCM("80010113", "0d200019a27300008250404000", "0162103500000000"); // min pT 100 cut
      cuts.AddCutPCM("80010113", "0d200029a27300008250404000", "0162103500000000"); // min pT 150 cut
    } else if (trainConfig == 103){
      cuts.AddCutPCM("80010113", "0d200009327300008250404000", "0162103500000000"); // edEdx -4,5
      cuts.AddCutPCM("80010113", "0d200009627300008250404000", "0162103500000000"); // edEdx -2.5,4
      cuts.AddCutPCM("80010113", "0d200009227300008250404000", "0162103500000000"); // edEdx -3,5
    } else if (trainConfig == 104){
      cuts.AddCutPCM("80010113", "0d200009a57300008250404000", "0162103500000000"); // pidEdx 2,-10
      cuts.AddCutPCM("80010113", "0d200009a17300008250404000", "0162103500000000"); // pidEdx 0,-10
      cuts.AddCutPCM("80010113", "0d200009a47300008250404000", "0162103500000000"); // pidEdx 3, 1
    } else if (trainConfig == 105){
      cuts.AddCutPCM("80010113", "0d200009b27300008250404000", "0162103500000000"); // edEdx -3,2,3.2
      cuts.AddCutPCM("80010113", "0d200009c27300008250404000", "0162103500000000"); // edEdx -2.8,2.8
      cuts.AddCutPCM("80010113", "0d200009a27100008250404000", "0162103500000000"); // pion nsig max mom 5.00 GeV/c
    } else if (trainConfig == 106){
      cuts.AddCutPCM("80010113", "0d200009a20300008250404000", "0162103500000000"); // pion nsig min mom 0.50 GeV/c
      cuts.AddCutPCM("80010113", "0d200009a26300008250404000", "0162103500000000"); // pion nsig min mom 0.25 GeV/c
      cuts.AddCutPCM("80010113", "0d200009a27600008250404000", "0162103500000000"); // pion nsig max mom 2.00 GeV/c
    } else if (trainConfig == 107){
      cuts.AddCutPCM("80010113", "0d200009a27300002250404000", "0162103500000000"); // qT max 0.06 2D
      cuts.AddCutPCM("80010113", "0d200009a27300009250404000", "0162103500000000"); // qT max 0.03 2D
      cuts.AddCutPCM("80010113", "0d200009a27000003250404000", "0162103500000000"); // qT max 0.05 1D
    } else if (trainConfig == 108){
      cuts.AddCutPCM("80010113", "0d200009a27300008260404000", "0162103500000000"); // Psi pair 0.05 2D, chi2 30.
      cuts.AddCutPCM("80010113", "0d200009a27300008280404000", "0162103500000000"); // Psi pair 0.2  2D, chi2 30.
      cuts.AddCutPCM("80010113", "0d200009a27300008254404000", "0162103500000000"); // Photon Asymmetry Cut
    } else if (trainConfig == 109){
      cuts.AddCutPCM("80010113", "0d200009a27300008850404000", "0162103500000000"); // chi2 20. psi pair 0.1 2D
      cuts.AddCutPCM("80010113", "0d200009a27300008860404000", "0162103500000000"); // variation chi2 20 psi pair 0.2 2D
      cuts.AddCutPCM("80010113", "0d200009a27300008150404000", "0162103500000000"); // chi2 50. psi pair 0.1 2D
    } else if (trainConfig == 110){
      cuts.AddCutPCM("80010113", "0d200009a27300008250400000", "0162103500000000"); // no double counting
      cuts.AddCutPCM("80010113", "0d200009a27300008250004000", "0162103500000000"); // no CosPA
      cuts.AddCutPCM("80010113", "0d200009a27300008250404000", "0162101500000000"); // meson alpha pt dep
      cuts.AddCutPCM("80010113", "0d200009a27300008250404000", "0162107500000000"); // meson alpha < 0.85
    } else if (trainConfig == 111){
      cuts.AddCutPCM("80010113", "0d200009a27300002252404000", "0162103500000000"); // variation qT max 0.06 2D, asym var 1
      cuts.AddCutPCM("80010113", "0d200009a27300002254404000", "0162103500000000"); // variation qT max 0.06 2D, asym vat 2 pt dep
      cuts.AddCutPCM("80010113", "0d200009a27300002256404000", "0162103500000000"); // variation qT max 0.06 2D, asym var 3
      //0-5
    } else if (trainConfig == 200){
      cuts.AddCutPCM("a0110113", "0d200009a27300008250404000", "0162103500000000"); // default 0005
      //0-20
    } else if (trainConfig == 300){
      cuts.AddCutPCM("80210113", "0d200009a27300008250404000", "0162103500000000"); // default 0020
      //20-40
    } else if (trainConfig == 400){
      cuts.AddCutPCM("82410113", "0d200009a27300008250404000", "0162103500000000"); // default 2040
      //40-60
    } else if (trainConfig == 500){
      cuts.AddCutPCM("84610113", "0d200009a27300008250404000", "0162103500000000"); // default 4060
      //60-80
    } else if (trainConfig == 600){
      cuts.AddCutPCM("86810113", "0d200009a27300008250404000", "0162103500000000"); // default 6080
      //80-100
    } else if (trainConfig == 700){
      cuts.AddCutPCM("88010113", "0d200009a27300008250404000", "0162103500000000"); // default 80100



    } else if (trainConfig == 440){// as 400 to be used MBW eta
      cuts.AddCutPCM("80010113", "00200009a27300008250404000", "0162103500000000"); // 
      cuts.AddCutPCM("80010113", "0c200009a27300008250404000", "0162103500000000"); // 
      cuts.AddCutPCM("80010113", "0d200009a27300008250404000", "0162103500000000"); // 
    } else if (trainConfig == 441){// as 440 to be used MBW
      cuts.AddCutPCM("80010113", "0da00009a27300008250404000", "0162103500000000"); // R 5-33.5   cm
      cuts.AddCutPCM("80010113", "0db00009a27300008250404000", "0162103500000000"); // R 33.5-72  cm
      cuts.AddCutPCM("80010113", "0dc00009a27300008250404000", "0162103500000000"); // R 72-180   cm
    } else if (trainConfig == 442){// further study material 
      cuts.AddCutPCM("80010113", "0dh00009a27300008250404000", "0162103500000000"); // R 95.0-180 cm Gas Volume
      cuts.AddCutPCM("80010113", "0di00009a27300008250404000", "0162103500000000"); // R 5-13.0   cm SPD
    } else if (trainConfig == 443){// default + category
      cuts.AddCutPCM("80010113", "0d200009a27300008250424000", "0162103500000000"); // default 0100 + Cat1
    } else if (trainConfig == 444){// default + category
      cuts.AddCutPCM("80010113", "0d200009a27300008250454000", "0162103500000000"); // default 0100 + Cat23
    } else if (trainConfig == 447){// small R
      cuts.AddCutPCM("80010113", "0dh00009a27300008250404000", "0162103500000000"); // R 5-13
      cuts.AddCutPCM("80010113", "0di00009a27300008250404000", "0162103500000000"); // R 13-33.5.
      cuts.AddCutPCM("80010113", "0dj00009a27300008250404000", "0162103500000000"); // R 33-55
    } else if (trainConfig == 448){// large R
      cuts.AddCutPCM("80010113", "0dk00009a27300008250404000", "0162103500000000"); // R 55-72
      cuts.AddCutPCM("80010113", "0dl00009a27300008250404000", "0162103500000000"); // R 72-95
      cuts.AddCutPCM("80010113", "0dg00009a27300008250404000", "0162103500000000"); // R 95-180
    } else if (trainConfig == 449){// Ana's request
      cuts.AddCutPCM("80010113", "0dm00009a27300008250404000", "0162103500000000"); // R 5-180, exclude 55-72
    } else if (trainConfig == 450){//
      cuts.AddCutPCM("80010113", "0d200079a27300008250404000", "0162103500000000"); // min pT no cut
      cuts.AddCutPCM("80010113", "0d200049a27300008250404000", "0162103500000000"); // min pT 75 cut
    } else if (trainConfig == 451){//
      cuts.AddCutPCM("80010113", "0d200019a27300008250404000", "0162103500000000"); // min pT 100 cut
      cuts.AddCutPCM("80010113", "0d200029a27300008250404000", "0162103500000000"); // min pT 150 cut
    } else if (trainConfig == 452){//
      cuts.AddCutPCM("80010113", "0d200001a27300008250404000", "0162103500000000"); // TPC cluster 35%
      cuts.AddCutPCM("80010113", "0d200002a27300008250404000", "0162103500000000"); // TPC cluster 70%
    } else if (trainConfig == 453){//
      cuts.AddCutPCM("80010113", "0d200009227300008250404000", "0162103500000000"); // edEdx -3,4
      cuts.AddCutPCM("80010113", "0d200009b27300008250404000", "0162103500000000"); // edEdx -3,2,3.2
      cuts.AddCutPCM("80010113", "0d200009c27300008250404000", "0162103500000000"); // edEdx -2.8,2.8
    } else if (trainConfig == 454){//R 5-33.5 cm
      cuts.AddCutPCM("80010113", "0da00009a27300008250404000", "0162103500000000"); // +-10
      cuts.AddCutPCM("80010114", "0da00009a27300008250404000", "0162103500000000"); // +-7.5
    } else if (trainConfig == 455){//R 5-33.5 cm
      cuts.AddCutPCM("80010115", "0da00009a27300008250404000", "0162103500000000"); // +-5.
      cuts.AddCutPCM("80010116", "0da00009a27300008250404000", "0162103500000000"); // +-2.5     
    } else if (trainConfig == 456){//cat 1 Meson selection
      cuts.AddCutPCM("80010113", "0d200009a27300008250404000", "0162103500000000"); // 5 < R < 180 (default)
      cuts.AddCutPCM("80010113", "0d200009a27300008250404000", "0162103510000000"); // shared electron
      cuts.AddCutPCM("80010113", "0d200009a27300008250404000", "0162103520000000"); // veto cat 1 Meson
    } else if (trainConfig == 457){//R max 180
      cuts.AddCutPCM("80010113", "0dn00009a27300008250404000", "0162103500000000"); // 10 < R < 180
      cuts.AddCutPCM("80010113", "0do00009a27300008250404000", "0162103500000000"); // 15 < R < 180
      cuts.AddCutPCM("80010113", "0dp00009a27300008250404000", "0162103500000000"); // 20 < R < 180
    } else if (trainConfig == 458){//R max 95
      cuts.AddCutPCM("80010113", "0dq00009a27300008250404000", "0162103500000000"); // 5  < R < 95
      cuts.AddCutPCM("80010113", "0dr00009a27300008250404000", "0162103500000000"); // 10 < R < 95
    } else if (trainConfig == 459){//R max 95
      cuts.AddCutPCM("80010113", "0ds00009a27300008250404000", "0162103500000000"); // 15 < R < 95
      cuts.AddCutPCM("80010113", "0dt00009a27300008250404000", "0162103500000000"); // 20 < R < 95
    } else {
      Error(Form("GammaConvV1_%i",trainConfig), "wrong trainConfig variable no cuts have been specified for the configuration");
      return NULL; 
    }
  }

  
  if(!cuts.AreValid()){
    cout << "\n\n****************************************************" << endl;
    cout << "ERROR: No valid cuts stored in CutHandlerPCM! Returning..." << endl;
    cout << "****************************************************\n\n" << endl;
    return NULL;
  }

  Int_t numberOfCuts = cuts.GetNCuts();

  TList *EventCutList = new TList();
  TList *ConvCutList = new TList();
  TList *MesonCutList = new TList();

  EventCutList->SetOwner(kTRUE);
  AliConvEventCuts **analysisEventCuts        = new AliConvEventCuts*[numberOfCuts];
  ConvCutList->SetOwner(kTRUE);
  AliConversionPhotonCuts **analysisCuts      = new AliConversionPhotonCuts*[numberOfCuts];
  MesonCutList->SetOwner(kTRUE);
  AliConversionMesonCuts **analysisMesonCuts  = new AliConversionMesonCuts*[numberOfCuts];
  Bool_t initializedMatBudWeigths_existing    = kFALSE;

  for(Int_t i = 0; i<numberOfCuts; i++){
    analysisEventCuts[i]          = new AliConvEventCuts();
    TString fitNamePi0            = "Pi0_Fit_Data_5TeV2017";
    TString fitNameEta            = "Eta_Fit_Data_5TeV2017";
    TString mcInputNamePi0        = "Pi0_Pythia8_LHC17pq_fastwoSDD_5TeV2017";
    TString mcInputNameEta        = "Eta_Pythia8_LHC17pq_fastwoSDD_5TeV2017";

    if (doPartWeight){
      analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE, fileNamePtWeights,
								   mcInputNamePi0, mcInputNameEta, "",
								   fitNamePi0,fitNameEta,"");
    }

    TString dataInputMultHisto    = "";
    TString mcInputMultHisto      = "";
    if (doMultWeight){
      cout << "INFO enableling mult weighting" << endl;
      if( periodNameAnchor.CompareTo("LHC15n")==0  ||
	  periodNameAnchor.CompareTo("LHC16d")==0  ||
	  periodNameAnchor.CompareTo("LHC17p")==0  ||
	  periodNameAnchor.CompareTo("LHC17q")==0  ||
	  periodNameAnchor.CompareTo("LHC16qt")==0  ){
	TString cutNumber = cuts.GetEventCut(i);
	TString centCut = cutNumber(0,3);  // first three digits of event cut
	dataInputMultHisto = Form("%s_%s", periodNameAnchor.Data(), centCut.Data());
	mcInputMultHisto   = Form("%s_%s", periodNameV0Reader.Data(), centCut.Data());
	cout<< "Histogram names data/MC:: "<< dataInputMultHisto.Data()<< " " << mcInputMultHisto.Data()<< endl;
	analysisEventCuts[i]->SetUseWeightMultiplicityFromFile( kTRUE, fileNameMultWeights, dataInputMultHisto, mcInputMultHisto );
      }
    }
    analysisEventCuts[i]->SetV0ReaderName(V0ReaderName);
    if (periodNameV0Reader.CompareTo("") != 0) analysisEventCuts[i]->SetPeriodEnum(periodNameV0Reader);
    analysisEventCuts[i]->SetLightOutput(enableLightOutput);
    analysisEventCuts[i]->InitializeCutsFromCutString((cuts.GetEventCut(i)).Data());
    EventCutList->Add(analysisEventCuts[i]);
    analysisEventCuts[i]->SetFillCutHistograms("",kFALSE);
  
    analysisCuts[i] = new AliConversionPhotonCuts();
    if (enableMatBudWeightsPi0 > 0){
      if (isMC > 0){
	if (analysisCuts[i]->InitializeMaterialBudgetWeights(enableMatBudWeightsPi0,fileNameMatBudWeights)){
	  initializedMatBudWeigths_existing = kTRUE;}
	else {cout << "ERROR The initialization of the materialBudgetWeights did not work out." << endl;}
      }
      else {cout << "ERROR 'enableMatBudWeightsPi0'-flag was set > 0 even though this is not a MC task. It was automatically reset to 0." << endl;}
    }

    analysisCuts[i]->ForceTPCRecalibrationAsFunctionOfConvR();

    if (doPostCalibration){
      if (isMC == 0){
	if(analysisCuts[i]->InitializeElecDeDxPostCalibration(fileNamedEdxPostCalib)){
	  analysisCuts[i]->SetDoElecDeDxPostCalibration(doPostCalibration);
	  cout << "Setting TPC dEdx post calibration file: " << fileNamedEdxPostCalib.Data() << endl;
        } else{
          doPostCalibration=kFALSE;
	  analysisCuts[i]->SetDoElecDeDxPostCalibration(doPostCalibration);
        }	
      }else{
        cout << "ERROR enableElecDeDxPostCalibration set to True even if MC file. Automatically reset to 0"<< endl;
        doPostCalibration=kFALSE;
        analysisCuts[i]->SetDoElecDeDxPostCalibration(kFALSE);
      }
    }

    analysisCuts[i]->SetV0ReaderName(V0ReaderName);
    analysisCuts[i]->SetLightOutput(enableLightOutput);
    analysisCuts[i]->InitializeCutsFromCutString((cuts.GetPhotonCut(i)).Data());

    ConvCutList->Add(analysisCuts[i]);
    analysisCuts[i]->SetFillCutHistograms("",kFALSE);
    analysisMesonCuts[i] = new AliConversionMesonCuts();
    analysisMesonCuts[i]->SetLightOutput(enableLightOutput);
    analysisMesonCuts[i]->InitializeCutsFromCutString((cuts.GetMesonCut(i)).Data());
    MesonCutList->Add(analysisMesonCuts[i]);
    analysisMesonCuts[i]->SetFillCutHistograms("");
  }

  task->SetDoTHnSparse(enableTHnSparse);
  task->SetEventCutList(numberOfCuts,EventCutList);
  task->SetConversionCutList(numberOfCuts,ConvCutList);
  task->SetMesonCutList(numberOfCuts,MesonCutList);
  task->SetMoveParticleAccordingToVertex(kTRUE);
  task->SetDoMesonAnalysis(kTRUE);
  task->SetDoMesonQA(enableQAMesonTask); //Attention new switch for Pi0 QA
  task->SetDoPhotonQA(enableQAPhotonTask);  //Attention new switch small for Photon QA
  if (initializedMatBudWeigths_existing) {
    task->SetDoMaterialBudgetWeightingOfGammasForTrueMesons(kTRUE);
  }

  //connect containers
  AliAnalysisDataContainer *coutput =
    mgr->CreateContainer(Form("GammaConvV1_%i",trainConfig), TList::Class(),
			 AliAnalysisManager::kOutputContainer, Form("GammaConvV1_%i.root",trainConfig) );
  mgr->AddTask(task);
  mgr->ConnectInput(task,0,cinput);
  mgr->ConnectOutput(task,1,coutput);
  Int_t nContainer = 2;
  for(Int_t i = 0; i<numberOfCuts; i++){
    if(enableQAPhotonTask>1){
      mgr->ConnectOutput(task,nContainer,mgr->CreateContainer(Form("%s_%s_%s Photon DCA tree",(cuts.GetEventCut(i)).Data(),(cuts.GetPhotonCut(i)).Data(),(cuts.GetMesonCut(i)).Data()), TTree::Class(), AliAnalysisManager::kOutputContainer, Form("GammaConvV1_%i.root",trainConfig)) );
      nContainer++;
    }
    if(enableQAMesonTask>1){
      mgr->ConnectOutput(task,nContainer,mgr->CreateContainer(Form("%s_%s_%s Meson DCA tree",(cuts.GetEventCut(i)).Data(),(cuts.GetPhotonCut(i)).Data(),(cuts.GetMesonCut(i)).Data()), TTree::Class(), AliAnalysisManager::kOutputContainer, Form("GammaConvV1_%i.root",trainConfig)) );
      nContainer++;
    }
  }
  return task;
}
