/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: Friederike Bock, Sandro Wenzel                                 *
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

//***************************************************************************************
//This AddTask is supposed to set up the main task
//($ALIPHYSICS/PWGGA/GammaConv/AliAnalysisTaskHeavyNeutralMesonToGG.cxx) for
//pPb together with all supporting classes
//***************************************************************************************

//***************************************************************************************
//main function
//***************************************************************************************
void AddTask_GammaHeavyMeson_CaloMode_pp(
  Int_t     selectedMeson                 = 0,        // select the corresponding meson: 0 pi0, 1 eta, 2 eta'
  Int_t     trainConfig                   = 1,        // change different set of cuts
  Int_t     isMC                          = 0,        // run MC
  TString   photonCutNumberV0Reader       = "",       // 00000008400000000100000000 nom. B, 00000088400000000100000000 low B
  TString   periodNameV0Reader            = "",
  // general setting for task
  Int_t     enableQAMesonTask             = 0,        // enable QA in AliAnalysisTaskGammaConvV1
  Int_t     enableQAClusterTask           = 0,        // enable additional QA task
  Int_t     enableExtMatchAndQA           = 0,        // disabled (0), extMatch (1), extQA_noCellQA (2), extMatch+extQA_noCellQA (3), extQA+cellQA (4), extMatch+extQA+cellQA (5)
  Int_t     enableLightOutput             = 0,        // switch to run light output (only essential histograms for afterburner)
  Bool_t    enableTHnSparse               = kFALSE,   // switch on THNsparse
  Int_t     enableTriggerMimicking        = 0,        // enable trigger mimicking
  Bool_t    enableTriggerOverlapRej       = kFALSE,   // enable trigger overlap rejection
  TString   settingMaxFacPtHard           = "3.",       // maximum factor between hardest jet and ptHard generated
  Int_t     debugLevel                    = 0,        // introducing debug levels for grid running
  // settings for weights
  // FPTW:fileNamePtWeights, FMUW:fileNameMultWeights, separate with ;
  TString   fileNameExternalInputs        = "",
  Int_t     doWeightingPart               = 0,        // enable Weighting
  TString   generatorName                 = "DPMJET", // generator Name
  Bool_t    enableMultiplicityWeighting   = kFALSE,   //
  TString   periodNameAnchor              = "",       //
  // special settings
  Bool_t    enableSortingMCLabels         = kTRUE,    // enable sorting for MC cluster labels
  // subwagon config
  TString   additionalTrainConfig         = "0"       // additional counter for trainconfig
) {

  AliCutHandlerPCM cuts(13);

  TString addTaskName                 = "AddTask_GammaHeavyMeson_CaloMode_pp";
  TString sAdditionalTrainConfig      = cuts.GetSpecialSettingFromAddConfig(additionalTrainConfig, "", "", addTaskName);
  if (sAdditionalTrainConfig.Atoi() > 0){
    trainConfig = trainConfig + sAdditionalTrainConfig.Atoi();
    cout << "INFO: " << addTaskName.Data() << " running additionalTrainConfig '" << sAdditionalTrainConfig.Atoi() << "', train config: '" << trainConfig << "'" << endl;
  }

  TString fileNamePtWeights           = cuts.GetSpecialFileNameFromString (fileNameExternalInputs, "FPTW:");
  TString fileNameMultWeights         = cuts.GetSpecialFileNameFromString (fileNameExternalInputs, "FMUW:");

  TString corrTaskSetting             = cuts.GetSpecialSettingFromAddConfig(additionalTrainConfig, "CF", "", addTaskName);
  if(corrTaskSetting.CompareTo(""))
    cout << "corrTaskSetting: " << corrTaskSetting.Data() << endl;

  Int_t trackMatcherRunningMode = 0; // CaloTrackMatcher running mode
  TString strTrackMatcherRunningMode  = cuts.GetSpecialSettingFromAddConfig(additionalTrainConfig, "TM", "", addTaskName);
  if(additionalTrainConfig.Contains("TM"))
    trackMatcherRunningMode = strTrackMatcherRunningMode.Atoi();

  Bool_t doTreeEOverP = kFALSE; // switch to produce EOverP tree
  TString strdoTreeEOverP             = cuts.GetSpecialSettingFromAddConfig(additionalTrainConfig, "EPCLUSTree", "", addTaskName);
  if(strdoTreeEOverP.Atoi()==1)
    doTreeEOverP = kTRUE;

  TH1S* histoAcc = 0x0;         // histo for modified acceptance
  TString strModifiedAcc              = cuts.GetSpecialSettingFromAddConfig(additionalTrainConfig, "MODIFYACC", "", addTaskName);
  if(strModifiedAcc.Contains("MODIFYACC")){
    TString tempType = strModifiedAcc;
    tempType.Replace(0,9,"");
    cout << "INFO: connecting to alien..." << endl;
    TGrid::Connect("alien://");
    cout << "done!" << endl;
    TFile *w = TFile::Open(fileNamePtWeights.Data());
    if(!w){cout << "ERROR: Could not open file: " << fileNamePtWeights.Data() << endl;return;}
    histoAcc = (TH1S*) w->Get(tempType.Data());
    if(!histoAcc) {cout << "ERROR: Could not find histo: " << tempType.Data() << endl;return;}
    cout << "found: " << histoAcc << endl;
  }

  Int_t localDebugFlag = 0;
  TString strLocalDebugFlag              = cuts.GetSpecialSettingFromAddConfig(additionalTrainConfig, "LOCALDEBUGFLAG", "", addTaskName);
  if(strLocalDebugFlag.Atoi()>0)
    localDebugFlag = strLocalDebugFlag.Atoi();

  TObjArray *rmaxFacPtHardSetting = settingMaxFacPtHard.Tokenize("_");
  if(rmaxFacPtHardSetting->GetEntries()<1){cout << "ERROR: AddTask_GammaHeavyMeson_CaloMode_pp during parsing of settingMaxFacPtHard String '" << settingMaxFacPtHard.Data() << "'" << endl; return;}
  Bool_t fMinPtHardSet        = kFALSE;
  Double_t minFacPtHard       = -1;
  Bool_t fMaxPtHardSet        = kFALSE;
  Double_t maxFacPtHard       = 100;
  Bool_t fSingleMaxPtHardSet  = kFALSE;
  Double_t maxFacPtHardSingle = 100;
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
    } else if(rmaxFacPtHardSetting->GetEntries()==1 && strTempSetting.Atof()>0){
      maxFacPtHard               = strTempSetting.Atof();
      cout << "running with max pT hard jet fraction of: " << maxFacPtHard << endl;
      fMaxPtHardSet        = kTRUE;
    }
  }

  Int_t isHeavyIon = 0;
  // meson reco mode: 0 - PCM-PCM, 1 - PCM-Calo, 2 - Calo-Calo
  Int_t mesonRecoMode = 2;

  // ================== GetAnalysisManager ===============================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error(Form("%s_%i", addTaskName.Data(),  trainConfig), "No analysis manager found.");
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
  //========= Add task to the ANALYSIS manager =====
  //================================================
  AliAnalysisTaskHeavyNeutralMesonToGG *task=NULL;
  task= new AliAnalysisTaskHeavyNeutralMesonToGG(Form("HeavyNeutralMesonToGG_%i_%i_%i", mesonRecoMode, selectedMeson, trainConfig));
  task->SetIsHeavyIon(isHeavyIon);
  task->SetIsMC(isMC);
  task->SetV0ReaderName(V0ReaderName);
  task->SetCorrectionTaskSetting(corrTaskSetting);

  task->SetMesonRecoMode(mesonRecoMode); // meson reco mode: 0 - PCM-PCM, 1 - PCM-Calo, 2 - Calo-Calo
  if (enableLightOutput > 1) task->SetLightOutput(kTRUE);
  task->SetDoPrimaryTrackMatching(kTRUE);

  // *********************************************************************************************************
  // 2.76 TeV  pp Run1 - EMC configurations
  // *********************************************************************************************************
  if (trainConfig == 1){ // pp 2.76 TeV LHC11a paper cuts
    cuts.AddCutCalo("00003013","1111121057032220000","0163103000000050"); // MB w/o pileup
    cuts.AddCutCalo("00003113","1111121057032220000","0163103000000050"); // MB
    cuts.AddCutCalo("00051013","1111121057032220000","0163103000000050"); // EMC1
  } else if (trainConfig == 2){  // pp 2.76 TeV LHC13g paper cuts
    cuts.AddCutCalo("00010113","1111121067032220000","0163103000000050");
    cuts.AddCutCalo("00010013","1111121067032220000","0163103000000050"); // without pile-up correction
    cuts.AddCutCalo("00052013","1111121067032220000","0163103000000050"); // EMC7
    cuts.AddCutCalo("00083013","1111121067032220000","0163103000000050"); // EMCEG1,
    cuts.AddCutCalo("00085013","1111121067032220000","0163103000000050"); // EMCEG2,
  } else if (trainConfig == 3){ // pp 2.76 TeV LHC11a paper cuts
    cuts.AddCutCalo("00003013","1111121057032220000","0163103b00000050"); // MB w/o pileup
    cuts.AddCutCalo("00003113","1111121057032220000","0163103b00000050"); // MB
    cuts.AddCutCalo("00051013","1111121057032220000","0163103b00000050"); // EMC1
  } else if (trainConfig == 4){  // pp 2.76 TeV LHC13g paper cuts
    cuts.AddCutCalo("00010113","1111121067032220000","0163103b00000050");
    cuts.AddCutCalo("00010013","1111121067032220000","0163103b00000050"); // without pile-up correction
    cuts.AddCutCalo("00052013","1111121067032220000","0163103b00000050"); // EMC7
    cuts.AddCutCalo("00083013","1111121067032220000","0163103b00000050"); // EMCEG1,
    cuts.AddCutCalo("00085013","1111121067032220000","0163103b00000050"); // EMCEG2,
  // *********************************************************************************************************
  // pp 8 TeV pp Run1 - EMC configurations
  // *********************************************************************************************************
  } else if (trainConfig == 100){ // EMCAL clusters pp 8 TeV paper cuts
    cuts.AddCutCalo("00010113","1111111067032220000","01631030000000d0"); // std
    cuts.AddCutCalo("00052113","1111111067032220000","01631030000000d0"); // std
    cuts.AddCutCalo("00081113","1111111067032220000","01631030000000d0"); // std
  } else if (trainConfig == 101){ // EMCAL clusters pp 8 TeV paper cuts
    cuts.AddCutCalo("00010113","1111111067032220000","01631030000000d0"); // std
    cuts.AddCutCalo("00052113","1111111067032220000","01631030000000d0"); // std
    cuts.AddCutCalo("00081113","1111111067032220000","01631030000000d0"); // std
  } else if (trainConfig == 102){ // EMCAL clusters pp 8 TeV paper cuts
    cuts.AddCutCalo("00010113","1111111067032220000","0163103b000000d0"); // std
    cuts.AddCutCalo("00052113","1111111067032220000","0163103b000000d0"); // std
    cuts.AddCutCalo("00081113","1111111067032220000","0163103b000000d0"); // std
  } else if (trainConfig == 103){ // EMCAL clusters pp 8 TeV paper cuts
    cuts.AddCutCalo("00010113","1111111067032220000","0163103b000000d0"); // std
    cuts.AddCutCalo("00052113","1111111067032220000","0163103b000000d0"); // std
    cuts.AddCutCalo("00081113","1111111067032220000","0163103b000000d0"); // std

  // *********************************************************************************************************
  // 7 TeV  pp Run1 LHC10x - EMC configurations
  // *********************************************************************************************************
  } else if (trainConfig == 200){ // EMCAL clusters pp 7 TeV, pT dep matching
    cuts.AddCutCalo("00000113","11111110b7032220000","01631030000000d0"); // std
    cuts.AddCutCalo("00000113","1111111007032220000","01631030000000d0"); // std
  } else if (trainConfig == 201){ // EMCAL clusters pp 7 TeV, pT dep matching
    cuts.AddCutCalo("00000113","11111110b7032220000","0163103b000000d0"); // std
    cuts.AddCutCalo("00000113","1111111007032220000","0163103b000000d0"); // std

  // *********************************************************************************************************
  // 5 TeV  pp Run2 - EMC configurations
  // *********************************************************************************************************
  } else if (trainConfig == 300){ // EMCAL clusters // -50ns, 30ns timing cut // no NL
    cuts.AddCutCalo("00010113","1111100067032220000","01631030000000d0");
    cuts.AddCutCalo("00052013","1111100067032220000","01631030000000d0");
    cuts.AddCutCalo("00085013","1111100067032220000","01631030000000d0");
    cuts.AddCutCalo("00083013","1111100067032220000","01631030000000d0");
  } else if (trainConfig == 301){ // EMCAL clusters // -50ns, 30ns timing cut // PCM-EMC NL
    cuts.AddCutCalo("00010113","1111111067032220000","01631030000000d0");
    cuts.AddCutCalo("00052013","1111111067032220000","01631030000000d0");
    cuts.AddCutCalo("00085013","1111111067032220000","01631030000000d0");
    cuts.AddCutCalo("00083013","1111111067032220000","01631030000000d0");
  } else if (trainConfig == 302){ // EMCAL clusters // -50ns, 30ns timing cut // no NL
    cuts.AddCutCalo("00010113","1111100067032220000","0163103b000000d0");
    cuts.AddCutCalo("00052013","1111100067032220000","0163103b000000d0");
    cuts.AddCutCalo("00085013","1111100067032220000","0163103b000000d0");
    cuts.AddCutCalo("00083013","1111100067032220000","0163103b000000d0");
  } else if (trainConfig == 303){ // EMCAL clusters // -50ns, 30ns timing cut // PCM-EMC NL
    cuts.AddCutCalo("00010113","1111111067032220000","0163103b000000d0");
    cuts.AddCutCalo("00052013","1111111067032220000","0163103b000000d0");
    cuts.AddCutCalo("00085013","1111111067032220000","0163103b000000d0");
    cuts.AddCutCalo("00083013","1111111067032220000","0163103b000000d0");

    // *********************************************************************************************************
    // 5 TeV  pp Run2 - DMC configurations
    // *********************************************************************************************************
  } else if (trainConfig == 350){ // DCAL clusters // -50ns, 30ns timing cut // no NL
    cuts.AddCutCalo("00010113","3885500067032220000","01631030000000d0");
    cuts.AddCutCalo("00055113","3885500067032220000","01631030000000d0");
    cuts.AddCutCalo("00089113","3885500067032220000","01631030000000d0");
    cuts.AddCutCalo("0008b113","3885500067032220000","01631030000000d0");
  } else if (trainConfig == 351){ // DCAL clusters // -50ns, 30ns timing cut // no NL
    cuts.AddCutCalo("00010113","3885500067032220000","0163103b000000d0");
    cuts.AddCutCalo("00055113","3885500067032220000","0163103b000000d0");
    cuts.AddCutCalo("00089113","3885500067032220000","0163103b000000d0");
    cuts.AddCutCalo("0008b113","3885500067032220000","0163103b000000d0");

  // *********************************************************************************************************
  // 13 TeV  pp Run2 - EMC configurations
  // *********************************************************************************************************
  } else if (trainConfig == 400){ // EMCAL clusters
    cuts.AddCutCalo("00010113","1111100017032220000","01631030000000d0"); // 1000ns timing cut, no NL INT7
    cuts.AddCutCalo("00010113","1111100067032220000","01631030000000d0"); // -30ns, 35ns timing cut, no NL INT7
  } else if (trainConfig == 401){ // EMCAL clusters
    cuts.AddCutCalo("00052113","1111100017032220000","01631030000000d0"); // 1000ns timing cut, no NL EMC7
    cuts.AddCutCalo("00085113","1111100017032220000","01631030000000d0"); // 1000ns timing cut, no NL EG2
    cuts.AddCutCalo("00083113","1111100017032220000","01631030000000d0"); // 1000ns timing cut, no NL EG1
  } else if (trainConfig == 402){ // EMCAL clusters
    cuts.AddCutCalo("00052113","1111100067032220000","01631030000000d0"); // -50ns, 30ns timing cut, no NL EMC7
    cuts.AddCutCalo("00085113","1111100067032220000","01631030000000d0"); // -50ns, 30ns timing cut, no NL EG2
    cuts.AddCutCalo("00083113","1111100067032220000","01631030000000d0"); // -50ns, 30ns timing cut, no NL EG1
  } else if (trainConfig == 403){ // EMCAL clusters, PCM-EMC NL // -50ns, 30ns timing cut
    cuts.AddCutCalo("00010113","1111111067032220000","01631030000000d0");
    cuts.AddCutCalo("00085113","1111111067032220000","01631030000000d0");
    cuts.AddCutCalo("00083113","1111111067032220000","01631030000000d0");
  } else if (trainConfig == 404){ // EMCAL clusters
    cuts.AddCutCalo("00010113","1111100017032220000","0163103b000000d0"); // 1000ns timing cut, no NL INT7
    cuts.AddCutCalo("00010113","1111100067032220000","0163103b000000d0"); // -30ns, 35ns timing cut, no NL INT7
  } else if (trainConfig == 405){ // EMCAL clusters
    cuts.AddCutCalo("00052113","1111100017032220000","0163103b000000d0"); // 1000ns timing cut, no NL EMC7
    cuts.AddCutCalo("00085113","1111100017032220000","0163103b000000d0"); // 1000ns timing cut, no NL EG2
    cuts.AddCutCalo("00083113","1111100017032220000","0163103b000000d0"); // 1000ns timing cut, no NL EG1
  } else if (trainConfig == 406){ // EMCAL clusters
    cuts.AddCutCalo("00052113","1111100067032220000","0163103b000000d0"); // -50ns, 30ns timing cut, no NL EMC7
    cuts.AddCutCalo("00085113","1111100067032220000","0163103b000000d0"); // -50ns, 30ns timing cut, no NL EG2
    cuts.AddCutCalo("00083113","1111100067032220000","0163103b000000d0"); // -50ns, 30ns timing cut, no NL EG1
  } else if (trainConfig == 407){ // EMCAL clusters, PCM-EMC NL // -50ns, 30ns timing cut
    cuts.AddCutCalo("00010113","1111111067032220000","0163103b000000d0");
    cuts.AddCutCalo("00085113","1111111067032220000","0163103b000000d0");
    cuts.AddCutCalo("00083113","1111111067032220000","0163103b000000d0");

  // *********************************************************************************************************
  // 13 TeV  DMC configurations
  // *********************************************************************************************************
  } else if (trainConfig == 450){ //DCAL
    cuts.AddCutCalo("00010113","3885500017032220000","01631030000000d0"); // 1000ns timing cut, no NL INT7
    cuts.AddCutCalo("00010113","3885500067032220000","01631030000000d0"); // -30ns, 35ns timing cut, no NL INT7
  } else if (trainConfig == 451){ // DCAL clusters
    cuts.AddCutCalo("00055113","3885500017032220000","01631030000000d0"); // 1000ns timing cut, no NL DMC7
    cuts.AddCutCalo("00089113","3885500017032220000","01631030000000d0"); // 1000ns timing cut, no NL DG2
    cuts.AddCutCalo("0008b113","3885500017032220000","01631030000000d0"); // 1000ns timing cut, no NL DG1
  } else if (trainConfig == 452){ // DCAL clusters
    cuts.AddCutCalo("00055113","3885500067032220000","01631030000000d0"); // -50ns, 30ns timing cut, no NL DMC7
    cuts.AddCutCalo("00089113","3885500067032220000","01631030000000d0"); // -50ns, 30ns timing cut, no NL DG2
    cuts.AddCutCalo("0008b113","3885500067032220000","01631030000000d0"); // -50ns, 30ns timing cut, no NL DG1
  } else if (trainConfig == 453){ //DCAL
    cuts.AddCutCalo("00010113","3885500017032220000","0163103b000000d0"); // 1000ns timing cut, no NL INT7
    cuts.AddCutCalo("00010113","3885500067032220000","0163103b000000d0"); // -30ns, 35ns timing cut, no NL INT7
  } else if (trainConfig == 454){ // DCAL clusters
    cuts.AddCutCalo("00055113","3885500017032220000","0163103b000000d0"); // 1000ns timing cut, no NL DMC7
    cuts.AddCutCalo("00089113","3885500017032220000","0163103b000000d0"); // 1000ns timing cut, no NL DG2
    cuts.AddCutCalo("0008b113","3885500017032220000","0163103b000000d0"); // 1000ns timing cut, no NL DG1
  } else if (trainConfig == 455){ // DCAL clusters
    cuts.AddCutCalo("00055113","3885500067032220000","0163103b000000d0"); // -50ns, 30ns timing cut, no NL DMC7
    cuts.AddCutCalo("00089113","3885500067032220000","0163103b000000d0"); // -50ns, 30ns timing cut, no NL DG2
    cuts.AddCutCalo("0008b113","3885500067032220000","0163103b000000d0"); // -50ns, 30ns timing cut, no NL DG1

  // *****************************************************************************************************
  // 2.76 TeV pp Run1 - PHOS configurations
  // *****************************************************************************************************
  } else if (trainConfig == 500) { //PHOS clusters
    cuts.AddCutCalo("00003113","2444400040033200000","0163803000000010"); //pp LHC11a with SDD, PHOS
    cuts.AddCutCalo("00010113","2444400040033200000","0163803000000010"); //pp LHC13g default MB
    cuts.AddCutCalo("00061113","2444400040033200000","0163803000000010"); //pp LHC11a PHI1
    cuts.AddCutCalo("00062113","2444400040033200000","0163803000000010"); //pp LHC11a PHI7
  } else if (trainConfig == 501) { //PHOS clusters
    cuts.AddCutCalo("00003113","2444400040033200000","0163803b00000010"); //pp LHC11a with SDD, PHOS
    cuts.AddCutCalo("00010113","2444400040033200000","0163803b00000010"); //pp LHC13g default MB
    cuts.AddCutCalo("00061113","2444400040033200000","0163803b00000010"); //pp LHC11a PHI1
    cuts.AddCutCalo("00062113","2444400040033200000","0163803b00000010"); //pp LHC11a PHI7

  // *****************************************************************************************************
  // 8 TeV pp Run1 - PHOS configurations
  // *****************************************************************************************************
  } else if (trainConfig == 600){ // PHOS clusters
    cuts.AddCutCalo("00010113","2444400040013300000","0163103000000010");
    cuts.AddCutCalo("00062113","2444400040013300000","0163103000000010");
  } else if (trainConfig == 601){ // PHOS clusters
    cuts.AddCutCalo("00010113","2444400040013300000","0163103b00000010");
    cuts.AddCutCalo("00062113","2444400040013300000","0163103b00000010");

  // *****************************************************************************************************
  // 7 TeV pp Run1 - PHOS configurations
  // *****************************************************************************************************
  } else if (trainConfig == 700){
    cuts.AddCutCalo("00000113","2444400000013300000","0163803000000010"); // QA
    cuts.AddCutCalo("00000113","2444400040013300000","0163803000000010"); // 100ns timing cut, no track matching
    cuts.AddCutCalo("00000113","2444400043013300000","0163803000000010"); // 100ns timing cut
  } else if (trainConfig == 701){
    cuts.AddCutCalo("00062113","2444400000013300000","0163803000000010"); // QA
    cuts.AddCutCalo("00062113","2444400040013300000","0163803000000010"); // 100ns timing cut, no track matching
    cuts.AddCutCalo("00062113","2444400043013300000","0163803000000010"); // 100ns timing cut
  } else if (trainConfig == 702){ // train config for bad channels and NonLin Variation
    cuts.AddCutCalo("00000113","2444400000013300000","0163803000000010"); // no NonLin
    cuts.AddCutCalo("00000113","2444401000013300000","0163803000000010"); // extern PHOS NonLin
    cuts.AddCutCalo("00000113","2444412000013300000","0163803000000010"); // constant non Lin first iteration
  } else if (trainConfig == 703){
    cuts.AddCutCalo("00000113","2444400000013300000","0163803b00000010"); // QA
    cuts.AddCutCalo("00000113","2444400040013300000","0163803b00000010"); // 100ns timing cut, no track matching
    cuts.AddCutCalo("00000113","2444400043013300000","0163803b00000010"); // 100ns timing cut
  } else if (trainConfig == 704){
    cuts.AddCutCalo("00062113","2444400000013300000","0163803b00000010"); // QA
    cuts.AddCutCalo("00062113","2444400040013300000","0163803b00000010"); // 100ns timing cut, no track matching
    cuts.AddCutCalo("00062113","2444400043013300000","0163803b00000010"); // 100ns timing cut
  } else if (trainConfig == 705){ // train config for bad channels and NonLin Variation
    cuts.AddCutCalo("00000113","2444400000013300000","0163803b00000010"); // no NonLin
    cuts.AddCutCalo("00000113","2444401000013300000","0163803b00000010"); // extern PHOS NonLin
    cuts.AddCutCalo("00000113","2444412000013300000","0163803b00000010"); // constant non Lin first iteration

  // *****************************************************************************************************
    // 5 TeV pp Run2 - PHOS configurations
  // *****************************************************************************************************
  } else if (trainConfig == 800){ // PHOS clusters with larger acceptance
    cuts.AddCutCalo("00010113","2446600040013300000","0163103000000010"); // INT7
    cuts.AddCutCalo("00062113","2446600040013300000","0163103000000010"); // PHI7
  } else if (trainConfig == 801){ // PHOS clusters with larger acceptance
    cuts.AddCutCalo("00010113","2446600040013300000","0163103b00000010"); // INT7
    cuts.AddCutCalo("00062113","2446600040013300000","0163103b00000010"); // PHI7

  // *****************************************************************************************************
  // 13 TeV pp Run2 - PHOS configurations
  // *****************************************************************************************************
  } else if (trainConfig == 900){ // PHOS clusters with larger acceptance
    cuts.AddCutCalo("00010113","2446600040013300000","0163103000000010"); // INT7
    cuts.AddCutCalo("00062113","2446600040013300000","0163103000000010"); // PHI7
  } else if (trainConfig == 901){ // PHOS clusters with larger acceptance
    cuts.AddCutCalo("00010113","2446600040013300000","0163103b00000010"); // INT7
    cuts.AddCutCalo("00062113","2446600040013300000","0163103b00000010"); // PHI7
  } else if (trainConfig == 902){ // PHOS for eta prime
    cuts.AddCutCalo("00010113","2446600043013300000","01631030000000d0"); // INT7
    cuts.AddCutCalo("00062113","2446600043013300000","01631030000000d0"); // PHI7

    // *****************************************************************************************************
    // 13 TeV pp Run2 - EDC (EMCal + DCal) configurations
    // *****************************************************************************************************

  } else if (trainConfig == 2020){ // EMCAL+DCAL clusters standard cuts, INT7, NL , std TM
    cuts.AddCutCalo("00010113","411792106f032230000","0r631031000000d0"); // INT7 NL 12 + TB
  } else if (trainConfig == 2021){ // EMCAL+DCAL clusters standard cuts, EG2, NL , std TM
    cuts.AddCutCalo("0008e113","411792106f032230000","0r631031000000d0"); // EG2  NL 12 + TB
  } else if (trainConfig == 2022){ // EMCAL+DCAL clusters standard cuts, EG1, NL , std TM
    cuts.AddCutCalo("0008d113","411792106f032230000","0r631031000000d0"); // EG1  NL 12 + TB
  } else if (trainConfig == 2023){ //EMCal + DCal INT7 cut var. NonLins
    cuts.AddCutCalo("00010113","411791106f032230000","01631031000000d0"); // INT7 NL11
    cuts.AddCutCalo("00010113","411792106f032230000","01631031000000d0"); // INT7 NL12
    cuts.AddCutCalo("00010113","411792206f032230000","01631031000000d0"); // INT7 NL22
    cuts.AddCutCalo("00010113","411793306f032230000","01631031000000d0"); // INT7 NL33
    cuts.AddCutCalo("00010113","411793406f032230000","01631031000000d0"); // INT7 NL34
  } else if (trainConfig == 2024){ //EMCal + DCal INT7 cut var. time
    cuts.AddCutCalo("00010113","411792105f032230000","01631031000000d0"); // INT7 time -50+50
    cuts.AddCutCalo("00010113","411792109f032230000","01631031000000d0"); // INT7 time -20+25
    cuts.AddCutCalo("00010113","411792108f032230000","01631031000000d0"); // INT7 time -20+30
  } else if (trainConfig == 2025){ //EMCal + DCal INT7 cut var. energy and NCell
    cuts.AddCutCalo("00010113","411792106f022230000","01631031000000d0"); // INT7 energy 0.6 GeV
    cuts.AddCutCalo("00010113","411792106f042230000","01631031000000d0"); // INT7 energy 0.8 GeV
    cuts.AddCutCalo("00010113","411792106f052230000","01631031000000d0"); // INT7 energy 0.9 GeV
    cuts.AddCutCalo("00010113","411792106f031230000","01631031000000d0"); // INT7 NCells 1
    cuts.AddCutCalo("00010113","411792106f033230000","01631031000000d0"); // INT7 NCells 3
  } else if (trainConfig == 2026){ //EMCal + DCal INT7 cut var. M02 and TM
    cuts.AddCutCalo("00010113","411792106f032220000","01631031000000d0"); // INT7 M02 0.7
    cuts.AddCutCalo("00010113","411792106f032250000","01631031000000d0"); // INT7 M02 0.3
    cuts.AddCutCalo("00010113","411792106f0322k0000","01631031000000d0"); // INT7 M02 E dep
    cuts.AddCutCalo("00010113","411792106e032230000","01631031000000d0"); // INT7 TM var
    cuts.AddCutCalo("00010113","411792106g032230000","01631031000000d0"); // INT7 TM var
    cuts.AddCutCalo("00010113","411792106h032230000","01631031000000d0"); // INT7 TM var
    cuts.AddCutCalo("00010113","4117921067032230000","01631031000000d0"); // INT7 TM var
  } else if (trainConfig == 2027){ //EMCal + DCal INT7 cut var. open. angle and alpha
    cuts.AddCutCalo("00010113","411792106f032230000","01631031000000b0"); // INT7 Op. Ang. var 1 cell dist + 0.0152
    cuts.AddCutCalo("00010113","411792106f032230000","01631031000000g0"); // INT7 Op. Ang. var 1 cell dist + 0.0202
    cuts.AddCutCalo("00010113","411792106f032230000","01631031000000a0"); // INT7 Op. Ang. var 1 cell dist + 0
    cuts.AddCutCalo("00010113","411792106f032230000","01631031000000h0"); // INT7 Op. Ang. var pT dep.
    cuts.AddCutCalo("00010113","411792106f032230000","01631051000000d0"); // INT7 alpha cut 0-0.75
    cuts.AddCutCalo("00010113","411792106f032230000","01631081000000d0"); // INT7 alpha cut 0-0.6

  } else if (trainConfig == 2028){ //EMCal + DCal EG2 cut var. NonLins
    cuts.AddCutCalo("0008e113","411791106f032230000","01631031000000d0"); // EG2 NL11
    cuts.AddCutCalo("0008e113","411791206f032230000","01631031000000d0"); // EG2 NL12
    cuts.AddCutCalo("0008e113","411792206f032230000","01631031000000d0"); // EG2 NL22
    cuts.AddCutCalo("0008e113","411793306f032230000","01631031000000d0"); // EG2 NL33
    cuts.AddCutCalo("0008e113","411793406f032230000","01631031000000d0"); // EG2 NL34
  } else if (trainConfig == 2029){ //EMCal + DCal EG2 cut var. time
    cuts.AddCutCalo("0008e113","411792105f032230000","01631031000000d0"); // EG2 time -50+50
    cuts.AddCutCalo("0008e113","411792109f032230000","01631031000000d0"); // EG2 time -20+25
    cuts.AddCutCalo("0008e113","411792108f032230000","01631031000000d0"); // EG2 time -20+30
  } else if (trainConfig == 2030){ //EMCal + DCal EG2 cut var. energy and NCell
    cuts.AddCutCalo("0008e113","411792106f022230000","01631031000000d0"); // EG2 energy 0.6 GeV
    cuts.AddCutCalo("0008e113","411792106f042230000","01631031000000d0"); // EG2 energy 0.8 GeV
    cuts.AddCutCalo("0008e113","411792106f052230000","01631031000000d0"); // EG2 energy 0.9 GeV
    cuts.AddCutCalo("0008e113","411792106f031230000","01631031000000d0"); // EG2 NCells 1
    cuts.AddCutCalo("0008e113","411792106f033230000","01631031000000d0"); // EG2 NCells 3
  } else if (trainConfig == 2031){ //EMCal + DCal EG2 cut var. M02 and TM
    cuts.AddCutCalo("0008e113","411792106f032220000","01631031000000d0"); // EG2 M02 0.7
    cuts.AddCutCalo("0008e113","411792106f032250000","01631031000000d0"); // EG2 M02 0.3
    cuts.AddCutCalo("0008e113","411792106f0322k0000","01631031000000d0"); // EG2 M02 E dep
    cuts.AddCutCalo("0008e113","411792106e032230000","01631031000000d0"); // EG2 TM var
    cuts.AddCutCalo("0008e113","411792106g032230000","01631031000000d0"); // EG2 TM var
    cuts.AddCutCalo("0008e113","411792106h032230000","01631031000000d0"); // EG2 TM var
    cuts.AddCutCalo("0008e113","4117921067032230000","01631031000000d0"); // EG2 TM var
  } else if (trainConfig == 2032){ //EMCal + DCal EG2 cut var. open. angle and alpha
    cuts.AddCutCalo("0008e113","411792106f032230000","01631031000000b0"); // EG2 Op. Ang. var 1 cell dist + 0.0152
    cuts.AddCutCalo("0008e113","411792106f032230000","01631031000000g0"); // EG2 Op. Ang. var 1 cell dist + 0.0202
    cuts.AddCutCalo("0008e113","411792106f032230000","01631031000000a0"); // EG2 Op. Ang. var 1 cell dist + 0
    cuts.AddCutCalo("0008e113","411792106f032230000","01631031000000h0"); // EG2 Op. Ang. var pT dep.
    cuts.AddCutCalo("0008e113","411792106f032230000","01631051000000d0"); // EG2 alpha cut 0-0.75
    cuts.AddCutCalo("0008e113","411792106f032230000","01631081000000d0"); // EG2 alpha cut 0-0.6

  } else if (trainConfig == 2033){ //EMCal + DCal EG1 cut var. NonLins
    cuts.AddCutCalo("0008d113","411791106f032230000","01631031000000d0"); // EG1 NL11
    cuts.AddCutCalo("0008d113","411791206f032230000","01631031000000d0"); // EG1 NL12
    cuts.AddCutCalo("0008d113","411792206f032230000","01631031000000d0"); // EG1 NL22
    cuts.AddCutCalo("0008d113","411793306f032230000","01631031000000d0"); // EG1 NL33
    cuts.AddCutCalo("0008d113","411793406f032230000","01631031000000d0"); // EG1 NL34
  } else if (trainConfig == 2034){ //EMCal + DCal EG1 cut var. time
    cuts.AddCutCalo("0008d113","411792105f032230000","01631031000000d0"); // EG1 time -50+50
    cuts.AddCutCalo("0008d113","411792109f032230000","01631031000000d0"); // EG1 time -20+25
    cuts.AddCutCalo("0008d113","411792108f032230000","01631031000000d0"); // EG1 time -20+30
  } else if (trainConfig == 2035){ //EMCal + DCal EG1 cut var. energy and NCell
    cuts.AddCutCalo("0008d113","411792106f022230000","01631031000000d0"); // EG1 energy 0.6 GeV
    cuts.AddCutCalo("0008d113","411792106f042230000","01631031000000d0"); // EG1 energy 0.8 GeV
    cuts.AddCutCalo("0008d113","411792106f052230000","01631031000000d0"); // EG1 energy 0.9 GeV
    cuts.AddCutCalo("0008d113","411792106f031230000","01631031000000d0"); // EG1 NCells 1
    cuts.AddCutCalo("0008d113","411792106f033230000","01631031000000d0"); // EG1 NCells 3
  } else if (trainConfig == 2036){ //EMCal + DCal EG1 cut var. M02 and TM
    cuts.AddCutCalo("0008d113","411792106f032220000","01631031000000d0"); // EG1 M02 0.7
    cuts.AddCutCalo("0008d113","411792106f032250000","01631031000000d0"); // EG1 M02 0.3
    cuts.AddCutCalo("0008d113","411792106f0322k0000","01631031000000d0"); // EG1 M02 E dep
    cuts.AddCutCalo("0008d113","411792106e032230000","01631031000000d0"); // EG1 TM var
    cuts.AddCutCalo("0008d113","411792106g032230000","01631031000000d0"); // EG1 TM var
    cuts.AddCutCalo("0008d113","411792106h032230000","01631031000000d0"); // EG1 TM var
    cuts.AddCutCalo("0008d113","4117921067032230000","01631031000000d0"); // EG1 TM var
  } else if (trainConfig == 2037){ //EMCal + DCal EG1 cut var. open. angle and alpha
    cuts.AddCutCalo("0008d113","411792106f032230000","01631031000000b0"); // EG1 Op. Ang. var 1 cell dist + 0.0152
    cuts.AddCutCalo("0008d113","411792106f032230000","01631031000000g0"); // EG1 Op. Ang. var 1 cell dist + 0.0202
    cuts.AddCutCalo("0008d113","411792106f032230000","01631031000000a0"); // EG1 Op. Ang. var 1 cell dist + 0
    cuts.AddCutCalo("0008d113","411792106f032230000","01631031000000h0"); // EG1 Op. Ang. var pT dep.
    cuts.AddCutCalo("0008d113","411792106f032230000","01631051000000d0"); // EG1 alpha cut 0-0.75
    cuts.AddCutCalo("0008d113","411792106f032230000","01631081000000d0"); // EG1 alpha cut 0-0.6
  } else if (trainConfig == 2038){ // //EMCal + DCal EG1 cut var. trigger mimick
    cuts.AddCutCalo("0008e113","411792106f032230000","01631031000000d0"); // EG2  NL 12 + TB
    cuts.AddCutCalo("0008d113","411792106f032230000","01631031000000d0"); // EG1  NL 12 + TB


  } else if (trainConfig == 2040){ //EMCal + DCal EG1 mult. diff.
    cuts.AddCutCalo("m0110113","411792106f032230000","01631031000000d0"); // INT7, NL12, mult. dep 0.0% - 1.0%
    cuts.AddCutCalo("m1210113","411792106f032230000","01631031000000d0"); // INT7, NL12, mult. dep 1.0% - 2.0%
    cuts.AddCutCalo("m2310113","411792106f032230000","01631031000000d0"); // INT7, NL12, mult. dep 2.0% - 3.0%
    cuts.AddCutCalo("m3410113","411792106f032230000","01631031000000d0"); // INT7, NL12, mult. dep 3.0% - 4.0%
    cuts.AddCutCalo("m4510113","411792106f032230000","01631031000000d0"); // INT7, NL12, mult. dep 4.0% - 5.0%
    cuts.AddCutCalo("m5710113","411792106f032230000","01631031000000d0"); // INT7, NL12, mult. dep 5.0% - 7.0%
    cuts.AddCutCalo("m7a10113","411792106f032230000","01631031000000d0"); // INT7, NL12, mult. dep 7.0% - 10.0%
  } else if (trainConfig == 2041){
    cuts.AddCutCalo("n1210113","411792106f032230000","01631031000000d0"); // INT7, NL12, mult. dep 10.0% - 20.0%
    cuts.AddCutCalo("n2510113","411792106f032230000","01631031000000d0"); // INT7, NL12, mult. dep 20.0% - 50.0%
    cuts.AddCutCalo("n5a10113","411792106f032230000","01631031000000d0"); // INT7, NL12, mult. dep 50.0% - 100.0%

  } else if (trainConfig == 2042){  // high mult trigger
    cuts.AddCutCalo("00076113","411792106f032230000","01631031000000d0"); // INT7, NL12, high mult V0M,
  } else if (trainConfig == 2043){ //EMCal + DCal EG1 mult. diff.
    cuts.AddCutCalo("q0176113","411792106f032230000","01631031000000d0"); // INT7, NL12, high mult V0M, mult. dep 0% - 0.1%
    cuts.AddCutCalo("q1276113","411792106f032230000","01631031000000d0"); // INT7, NL12, high mult V0M, mult. dep 0.1% - 0.2%
    cuts.AddCutCalo("q2376113","411792106f032230000","01631031000000d0"); // INT7, NL12, high mult V0M, mult. dep 0.2% - 0.3%
    cuts.AddCutCalo("q3476113","411792106f032230000","01631031000000d0"); // INT7, NL12, high mult V0M, mult. dep 0.3% - 0.4%
    cuts.AddCutCalo("q4576113","411792106f032230000","01631031000000d0"); // INT7, NL12, high mult V0M, mult. dep 0.4% - 0.5%
    cuts.AddCutCalo("q5776113","411792106f032230000","01631031000000d0"); // INT7, NL12, high mult V0M, mult. dep 0.5% - 0.7%
    cuts.AddCutCalo("q7a76113","411792106f032230000","01631031000000d0"); // INT7, NL12, high mult V0M, mult. dep 0.7% - 1.0%
  } else if (trainConfig == 2044){
    cuts.AddCutCalo("m0176113","411792106f032230000","01631031000000d0"); // INT7, NL12, high mult V0M, mult. dep 0.0% - 1.0%
    cuts.AddCutCalo("m1276113","411792106f032230000","01631031000000d0"); // INT7, NL12, high mult V0M, mult. dep 1.0% - 2.0%
    cuts.AddCutCalo("m2376113","411792106f032230000","01631031000000d0"); // INT7, NL12, high mult V0M, mult. dep 2.0% - 3.0%
    cuts.AddCutCalo("m3476113","411792106f032230000","01631031000000d0"); // INT7, NL12, high mult V0M, mult. dep 3.0% - 4.0%
    cuts.AddCutCalo("m4576113","411792106f032230000","01631031000000d0"); // INT7, NL12, high mult V0M, mult. dep 4.0% - 5.0%
    cuts.AddCutCalo("m5776113","411792106f032230000","01631031000000d0"); // INT7, NL12, high mult V0M, mult. dep 5.0% - 7.0%
    cuts.AddCutCalo("m7a76113","411792106f032230000","01631031000000d0"); // INT7, NL12, high mult V0M, mult. dep 7.0% - 10.0%

  } else if (trainConfig == 2045){ //EMCal + DCal EG1 sphericity.
    cuts.AddCutCalo("h0510113","411791106f032230000","01631031000000d0"); // INT7, NL12, sphericity
    cuts.AddCutCalo("h5a10113","411791106f032230000","01631031000000d0"); // INT7, NL12, sphericity
  } else if (trainConfig == 2046){ //EMCal + DCal EG1 sphericity.
    cuts.AddCutCalo("h058e113","411791106f032230000","01631031000000d0"); // EG2, NL12, sphericity
    cuts.AddCutCalo("h5a8e113","411791106f032230000","01631031000000d0"); // EG2, NL12, sphericity
  } else if (trainConfig == 2047){ //EMCal + DCal EG1 sphericity.
    cuts.AddCutCalo("h058d113","411791106f032230000","01631031000000d0"); // EG1, NL12, sphericity
    cuts.AddCutCalo("h5a8d113","411791106f032230000","01631031000000d0"); // EG1, NL12, sphericity



  } else if (trainConfig == 2100){ // EMCAL+DCAL clusters, for eta prime
    cuts.AddCutCalo("00010113","4117911067032230000","01631030000000d0"); // no NL INT7
    cuts.AddCutCalo("0008e113","4117911067032230000","01631030000000d0"); // no NL EG2
    cuts.AddCutCalo("0008d113","4117911067032230000","01631030000000d0"); // no NL EG1
  } else {
    Error(Form("HeavyNeutralMesonToGG_%i_%i", mesonRecoMode, trainConfig), "wrong trainConfig variable no cuts have been specified for the configuration");
    return;
  }

  if(!cuts.AreValid()){
    cout << "\n\n****************************************************" << endl;
    cout << "ERROR: No valid cuts stored in CutHandlerHeavyMesonCalo! Returning..." << endl;
    cout << "****************************************************\n\n" << endl;
    return;
  }

  TList *EventCutList   = new TList();
  TList *ConvCutList    = new TList();
  TList *ClusterCutList = new TList();
  TList *MesonCutList   = new TList();

  Int_t numberOfCuts = cuts.GetNCuts();

  TList *HeaderList = new TList();
  if (generatorName.Contains("LHC12i3")){
    TObjString *Header2 = new TObjString("BOX");
    HeaderList->Add(Header2);
  } else if (generatorName.CompareTo("LHC14e2b")==0){
    TObjString *Header2 = new TObjString("pi0_1");
    HeaderList->Add(Header2);
    TObjString *Header3 = new TObjString("eta_2");
    HeaderList->Add(Header3);
  }

  TString energy = "";
  TString mcName = "";
  TString mcNameAdd = "";
  if (generatorName.Contains("WOSDD")){
    mcNameAdd = "_WOSDD";
  } else if (generatorName.Contains("WSDD")){
    mcNameAdd = "_WSDD";
  }
  if (generatorName.Contains("LHC12i3")){
    energy = "2760GeV";
    mcName = "Pythia8_LHC12i3";
  } else if (generatorName.Contains("LHC12f1a")){
    energy = "2760GeV";
    mcName = "Pythia8_LHC12f1a";
  } else if (generatorName.Contains("LHC12f1b")){
    energy = "2760GeV";
    mcName = "Phojet_LHC12f1b";
  } else if (generatorName.Contains("LHC14e2a")){
    energy = "8TeV";
    mcName = "Pythia8_LHC14e2a";
  } else if (generatorName.Contains("LHC14e2b")){
    energy = "8TeV";
    mcName = "Pythia8_LHC14e2b";
  } else if (generatorName.Contains("LHC14e2c")){
    energy = "8TeV";
    mcName = "Phojet_LHC14e2c";
  }

  EventCutList->SetOwner(kTRUE);
  AliConvEventCuts **analysisEventCuts        = new AliConvEventCuts*[numberOfCuts];
  ConvCutList->SetOwner(kTRUE);
  AliCaloPhotonCuts **analysisClusterCuts     = new AliCaloPhotonCuts*[numberOfCuts];
  MesonCutList->SetOwner(kTRUE);
  AliConversionMesonCuts **analysisMesonCuts  = new AliConversionMesonCuts*[numberOfCuts];

  for(Int_t i = 0; i<numberOfCuts; i++){
    //create AliCaloTrackMatcher instance, if there is none present
    TString caloCutPos = cuts.GetClusterCut(i);
    caloCutPos.Resize(1);
    TString TrackMatcherName = Form("CaloTrackMatcher_%s",caloCutPos.Data());
    if( !(AliCaloTrackMatcher*)mgr->GetTask(TrackMatcherName.Data()) ){
      AliCaloTrackMatcher* fTrackMatcher = new AliCaloTrackMatcher(TrackMatcherName.Data(),caloCutPos.Atoi());
      fTrackMatcher->SetV0ReaderName(V0ReaderName);
      fTrackMatcher->SetCorrectionTaskSetting(corrTaskSetting);
      mgr->AddTask(fTrackMatcher);
      mgr->ConnectInput(fTrackMatcher,0,cinput);
    }

    analysisEventCuts[i] = new AliConvEventCuts();

    TString dataInputMultHisto  = "";
    TString mcInputMultHisto    = "";
    TString triggerString       = cuts.GetEventCut(i);
    triggerString               = triggerString(3,2);

    dataInputMultHisto          = Form("%s_%s", periodNameAnchor.Data(), triggerString.Data());
    mcInputMultHisto            = Form("%s_%s", periodNameV0Reader.Data(), triggerString.Data());

    if (enableMultiplicityWeighting){
      cout << "enabling mult weighting" << endl;
      analysisEventCuts[i]->SetUseWeightMultiplicityFromFile( kTRUE, fileNameMultWeights, dataInputMultHisto, mcInputMultHisto );
    }

    // definition of weighting input
    TString fitNamePi0 = Form("Pi0_Fit_Data_%s",energy.Data());
    TString fitNameEta = Form("Eta_Fit_Data_%s",energy.Data());
    TString fAddedSignalString = cuts.GetEventCut(i);
    fAddedSignalString = fAddedSignalString(6,1);
    Bool_t fAddedSignal = kFALSE;
    if (fAddedSignalString.CompareTo("2") == 0) fAddedSignal = kTRUE;
    TString mcInputNamePi0 = "";
    TString mcInputNameEta = "";
    if (fAddedSignal && (generatorName.Contains("LHC12i3") || generatorName.CompareTo("LHC14e2b")==0)){
      mcInputNamePi0 = Form("Pi0_%s%s_addSig_%s", mcName.Data(), mcNameAdd.Data(), energy.Data() );
      mcInputNameEta = Form("Eta_%s%s_addSig_%s", mcName.Data(), mcNameAdd.Data(), energy.Data() );
    } else {
      mcInputNamePi0 = Form("Pi0_%s%s_%s", mcName.Data(), mcNameAdd.Data(), energy.Data() );
      mcInputNameEta = Form("Eta_%s%s_%s", mcName.Data(), mcNameAdd.Data(), energy.Data() );
    }

    if (doWeightingPart) analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile( kTRUE, kTRUE, kFALSE, fileNamePtWeights,
      mcInputNamePi0, mcInputNameEta, "",fitNamePi0,fitNameEta);


    analysisEventCuts[i]->SetTriggerMimicking(enableTriggerMimicking);
    analysisEventCuts[i]->SetTriggerOverlapRejecion(enableTriggerOverlapRej);
    if(fMinPtHardSet)
      analysisEventCuts[i]->SetMinFacPtHard(minFacPtHard);
    if(fMaxPtHardSet)
      analysisEventCuts[i]->SetMaxFacPtHard(maxFacPtHard);
    if(fSingleMaxPtHardSet)
      analysisEventCuts[i]->SetMaxFacPtHardSingleParticle(maxFacPtHardSingle);
    analysisEventCuts[i]->SetV0ReaderName(V0ReaderName);
    analysisEventCuts[i]->SetCorrectionTaskSetting(corrTaskSetting);
    if (periodNameV0Reader.CompareTo("") != 0) analysisEventCuts[i]->SetPeriodEnum(periodNameV0Reader);
    if (enableLightOutput > 0) analysisEventCuts[i]->SetLightOutput(kTRUE);
    analysisEventCuts[i]->InitializeCutsFromCutString((cuts.GetEventCut(i)).Data());
    EventCutList->Add(analysisEventCuts[i]);
    analysisEventCuts[i]->SetFillCutHistograms("",kFALSE);

    analysisClusterCuts[i] = new AliCaloPhotonCuts(isMC);
    analysisClusterCuts[i]->SetHistoToModifyAcceptance(histoAcc);
    analysisClusterCuts[i]->SetV0ReaderName(V0ReaderName);
    analysisClusterCuts[i]->SetCorrectionTaskSetting(corrTaskSetting);
    analysisClusterCuts[i]->SetCaloTrackMatcherName(TrackMatcherName);
    if (enableLightOutput > 0) analysisClusterCuts[i]->SetLightOutput(kTRUE);
    analysisClusterCuts[i]->InitializeCutsFromCutString((cuts.GetClusterCut(i)).Data());
    ClusterCutList->Add(analysisClusterCuts[i]);
    analysisClusterCuts[i]->SetExtendedMatchAndQA(enableExtMatchAndQA);
    analysisClusterCuts[i]->SetFillCutHistograms("");

    analysisMesonCuts[i] = new AliConversionMesonCuts();
    analysisMesonCuts[i]->SetLightOutput(enableLightOutput);
    analysisMesonCuts[i]->InitializeCutsFromCutString((cuts.GetMesonCut(i)).Data());
    analysisMesonCuts[i]->SetIsMergedClusterCut(2);
    analysisMesonCuts[i]->SetCaloMesonCutsObject(analysisClusterCuts[i]);
    MesonCutList->Add(analysisMesonCuts[i]);
    analysisMesonCuts[i]->SetFillCutHistograms("");
    analysisEventCuts[i]->SetAcceptedHeader(HeaderList);
  }

  task->SetEventCutList(numberOfCuts,EventCutList);
  task->SetCaloCutList(numberOfCuts,ClusterCutList);
  task->SetMesonCutList(numberOfCuts,MesonCutList);
  task->SetMoveParticleAccordingToVertex(kTRUE);
  task->SetCorrectionTaskSetting(corrTaskSetting);
  task->SetMesonType(selectedMeson);
  task->SetDoMesonQA(enableQAMesonTask); //Attention new switch for Pi0 QA
  task->SetDoPhotonQA(enableQAClusterTask);  //Attention new switch small for Photon QA
  task->SetDoClusterQA(1);  //Attention new switch small for Cluster QA
  task->SetUseTHnSparse(enableTHnSparse);
  task->SetEnableSortingOfMCClusLabels(enableSortingMCLabels);
  if(enableExtMatchAndQA > 1){ task->SetPlotHistsExtQA(kTRUE);}

  //connect containers
  AliAnalysisDataContainer *coutput =
  mgr->CreateContainer(!(corrTaskSetting.CompareTo("")) ? Form("HeavyNeutralMesonToGG_%i_%i_%i",mesonRecoMode, selectedMeson, trainConfig)
                                                        : Form("HeavyNeutralMesonToGG_%i_%i_%i_%s", mesonRecoMode, selectedMeson, trainConfig, corrTaskSetting.Data()), TList::Class(),
                        AliAnalysisManager::kOutputContainer,Form("HeavyNeutralMesonToGG_%i_%i.root",mesonRecoMode,trainConfig) );

  mgr->AddTask(task);
  mgr->ConnectInput(task,0,cinput);
  mgr->ConnectOutput(task,1,coutput);

  return;

}
