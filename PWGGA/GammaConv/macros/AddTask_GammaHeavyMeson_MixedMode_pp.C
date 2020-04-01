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
void AddTask_GammaHeavyMeson_MixedMode_pp(
  Int_t     selectedMeson                 = 0,        // select the corresponding meson: 0 pi0, 1 eta, 2 eta'
  Int_t     trainConfig                   = 1,        // change different set of cuts
  Int_t     isMC                          = 0,        // run MC
  TString   photonCutNumberV0Reader       = "",       // 00000008400000000100000000 nom. B, 00000088400000000100000000 low B
  TString   periodNameV0Reader            = "",
  // general setting for task
  Int_t     enableQAMesonTask             = 0,        // enable QA in AliAnalysisTaskGammaConvV1
  Int_t     enableQAPhotonTask            = 0,        // enable additional QA task
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
  Bool_t    enableElecDeDxPostCalibration = kFALSE,
  TString   periodNameAnchor              = "",       //
  // special settings
  Bool_t    enableSortingMCLabels         = kTRUE,    // enable sorting for MC cluster labels
  Bool_t    enableTreeConvGammaShape      = kFALSE,   // enable additional tree for conversion properties for clusters
  Bool_t    doSmear                       = kFALSE,   // switches to run user defined smearing
  Double_t  bremSmear                     = 1.,
  Double_t  smearPar                      = 0.,       // conv photon smearing params
  Double_t  smearParConst                 = 0.,       // conv photon smearing params
  Bool_t    doPrimaryTrackMatching        = kTRUE,    // enable basic track matching for all primary tracks to cluster
  // subwagon config
  TString   additionalTrainConfig         = "0"       // additional counter for trainconfig
) {

  AliCutHandlerPCM cuts;

  TString fileNamePtWeights           = cuts.GetSpecialFileNameFromString (fileNameExternalInputs, "FPTW:");
  TString fileNameMultWeights         = cuts.GetSpecialFileNameFromString (fileNameExternalInputs, "FMUW:");
  TString fileNamedEdxPostCalib       = cuts.GetSpecialFileNameFromString (fileNameExternalInputs, "FEPC:");

  TString addTaskName                 = "AddTask_GammaHeavyMeson_MixedMode_pp";
  TString sAdditionalTrainConfig      = cuts.GetSpecialSettingFromAddConfig(additionalTrainConfig, "", "", addTaskName);
  if (sAdditionalTrainConfig.Atoi() > 0){
    trainConfig = trainConfig + sAdditionalTrainConfig.Atoi();
    cout << "INFO: " << addTaskName.Data() << " running additionalTrainConfig '" << sAdditionalTrainConfig.Atoi() << "', train config: '" << trainConfig << "'" << endl;
  }
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

  Bool_t doTreeClusterShowerShape = kFALSE; // switch to produce EOverP tree
  TString strdoTreeClusterShowerShape             = cuts.GetSpecialSettingFromAddConfig(additionalTrainConfig, "INVMASSCLUSTree", "", addTaskName);
  if(strdoTreeClusterShowerShape.Atoi()==1)
    doTreeClusterShowerShape = kTRUE;

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
  if(rmaxFacPtHardSetting->GetEntries()<1){cout << "ERROR: AddTask_GammaHeavyMeson_MixedMode_pp during parsing of settingMaxFacPtHard String '" << settingMaxFacPtHard.Data() << "'" << endl; return;}
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
  Int_t mesonRecoMode = 1;

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
  task->SetDoPrimaryTrackMatching(doPrimaryTrackMatching);

  // *****************************************************************************************************
  // pp 2.76 TeV EMC configurations, pi0/eta paper cuts -
  // *****************************************************************************************************
  if (trainConfig == 1){ // EMCAL clusters 2.76 TeV LHC11a, with SDD (0), kEMC1 (1) final analysis cuts
    cuts.AddCutPCMCalo("00003113","00200009327000008250400000","1111121057032230000","0163103000000010");
    cuts.AddCutPCMCalo("00051013","00200009327000008250400000","1111121057032230000","0163103000000010");
  } else if (trainConfig == 2){  // LHC13g final analysis cuts
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111121067032230000","0163103000000010"); // INT7
    cuts.AddCutPCMCalo("00052013","00200009327000008250400000","1111121067032230000","0163103000000010"); // EMC7
    cuts.AddCutPCMCalo("00083013","00200009327000008250400000","1111121067032230000","0163103000000010"); // EMCEG1,
    cuts.AddCutPCMCalo("00085013","00200009327000008250400000","1111121067032230000","0163103000000010"); // EMCEG2,
  } else if (trainConfig == 3){ // EMCAL clusters 2.76 TeV LHC11a, with SDD (0), kEMC1 (1) final analysis cuts
    cuts.AddCutPCMCalo("00003113","00200009327000008250400000","1111121057032230000","0163103b00000010");
    cuts.AddCutPCMCalo("00051013","00200009327000008250400000","1111121057032230000","0163103b00000010");
  } else if (trainConfig == 4){  // LHC13g final analysis cuts
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111121067032230000","0163103b00000010"); // INT7
    cuts.AddCutPCMCalo("00052013","00200009327000008250400000","1111121067032230000","0163103b00000010"); // EMC7
    cuts.AddCutPCMCalo("00083013","00200009327000008250400000","1111121067032230000","0163103b00000010"); // EMCEG1,
    cuts.AddCutPCMCalo("00085013","00200009327000008250400000","1111121067032230000","0163103b00000010"); // EMCEG2,

  // *****************************************************************************************************
  // pp 8 TeV EMC configurations, pi0/eta paper cuts
  // *****************************************************************************************************
  } else if (trainConfig == 100){ // EMCAL clusters 8 TeV LHC12
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111111067032230000","0163103000000010"); // std
    cuts.AddCutPCMCalo("00052113","00200009327000008250400000","1111111067032230000","0163103000000010"); // std
    cuts.AddCutPCMCalo("00081113","00200009327000008250400000","1111111067032230000","0163103000000010"); // std
  } else if (trainConfig == 101){ // EMCAL clusters 8 TeV LHC12
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111111067032230000","0163103b00000010"); // std
    cuts.AddCutPCMCalo("00052113","00200009327000008250400000","1111111067032230000","0163103b00000010"); // std
    cuts.AddCutPCMCalo("00081113","00200009327000008250400000","1111111067032230000","0163103b00000010"); // std

  // *****************************************************************************************************
  // pp 7 TeV LHC10x EMC configurations
  // *****************************************************************************************************
  } else if (trainConfig == 200){
    cuts.AddCutPCMCalo("00000113","00200009327000008250400000","11111110b7032230000","0163103000000010"); // std
    cuts.AddCutPCMCalo("00000113","00200009327000008250400000","1111111007032230000","0163103000000010"); // std
  } else if (trainConfig == 201){
    cuts.AddCutPCMCalo("00000113","00200009327000008250400000","11111110b7032230000","0163103b00000010"); // std
    cuts.AddCutPCMCalo("00000113","00200009327000008250400000","1111111007032230000","0163103b00000010"); // std

  // *********************************************************************************************************
  // 5 TeV  pp Run2 - EMC configurations
  // *********************************************************************************************************
  } else if (trainConfig == 300){ // EMCAL clusters // -50ns, 30ns timing cut
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111100017032230000","0163103000000010");
    cuts.AddCutPCMCalo("00052013","00200009327000008250400000","1111100017032230000","0163103000000010");
    cuts.AddCutPCMCalo("00085013","00200009327000008250400000","1111100017032230000","0163103000000010");
    cuts.AddCutPCMCalo("00083013","00200009327000008250400000","1111100017032230000","0163103000000010");
  } else if (trainConfig == 301){ // EMCAL clusters // -50ns, 30ns timing cut
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111100067032230000","0163103000000010");
    cuts.AddCutPCMCalo("00052013","00200009327000008250400000","1111100067032230000","0163103000000010");
    cuts.AddCutPCMCalo("00085013","00200009327000008250400000","1111100067032230000","0163103000000010");
    cuts.AddCutPCMCalo("00083013","00200009327000008250400000","1111100067032230000","0163103000000010");
  } else if (trainConfig == 302){ // EMCAL clusters // -50ns, 30ns timing cut
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111100017032230000","0163103b00000010");
    cuts.AddCutPCMCalo("00052013","00200009327000008250400000","1111100017032230000","0163103b00000010");
    cuts.AddCutPCMCalo("00085013","00200009327000008250400000","1111100017032230000","0163103b00000010");
    cuts.AddCutPCMCalo("00083013","00200009327000008250400000","1111100017032230000","0163103b00000010");
  } else if (trainConfig == 303){ // EMCAL clusters // -50ns, 30ns timing cut
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111100067032230000","0163103b00000010");
    cuts.AddCutPCMCalo("00052013","00200009327000008250400000","1111100067032230000","0163103b00000010");
    cuts.AddCutPCMCalo("00085013","00200009327000008250400000","1111100067032230000","0163103b00000010");
    cuts.AddCutPCMCalo("00083013","00200009327000008250400000","1111100067032230000","0163103b00000010");

  // *********************************************************************************************************
  // 5 TeV  pp Run2 - DMC configurations
  // *********************************************************************************************************
  } else if (trainConfig == 350){ // DCAL clusters
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","3885500081041220000","0163103000000010"); //
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","3885500011041220000","0163103000000010"); //
  } else if (trainConfig == 351){ // DCAL clusters // -1000ns, 1000ns timing cut // no NL
    cuts.AddCutPCMCalo("00055113","00200009327000008250400000","3885500011041220000","0163103000000010");
    cuts.AddCutPCMCalo("00089113","00200009327000008250400000","3885500011041220000","0163103000000010");
    cuts.AddCutPCMCalo("0008b113","00200009327000008250400000","3885500011041220000","0163103000000010");
  } else if (trainConfig == 352){ // DCAL clusters // -50ns, 30ns timing cut // no NL
    cuts.AddCutPCMCalo("00055113","00200009327000008250400000","3885500081041220000","0163103000000010");
    cuts.AddCutPCMCalo("00089113","00200009327000008250400000","3885500011041220000","0163103000000010");
    cuts.AddCutPCMCalo("0008b113","00200009327000008250400000","3885500081041220000","0163103000000010");
  } else if (trainConfig == 353){ // DCAL clusters
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","3885500081041220000","0163103b00000010"); //
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","3885500011041220000","0163103b00000010"); //
  } else if (trainConfig == 354){ // DCAL clusters // -1000ns, 1000ns timing cut // no NL
    cuts.AddCutPCMCalo("00055113","00200009327000008250400000","3885500011041220000","0163103b00000010");
    cuts.AddCutPCMCalo("00089113","00200009327000008250400000","3885500011041220000","0163103b00000010");
    cuts.AddCutPCMCalo("0008b113","00200009327000008250400000","3885500011041220000","0163103b00000010");
  } else if (trainConfig == 355){ // DCAL clusters // -50ns, 30ns timing cut // no NL
    cuts.AddCutPCMCalo("00055113","00200009327000008250400000","3885500081041220000","0163103b00000010");
    cuts.AddCutPCMCalo("00089113","00200009327000008250400000","3885500011041220000","0163103b00000010");
    cuts.AddCutPCMCalo("0008b113","00200009327000008250400000","3885500081041220000","0163103b00000010");

  // *********************************************************************************************************
  // 13 TeV  pp Run2 - EMC configurations
  // *********************************************************************************************************
  } else if (trainConfig == 400){ // EMCAL clusters
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111100017032230000","0163103000000010");
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111100067032230000","0163103000000010");
  } else if (trainConfig == 401){ // EMCAL clusters
    cuts.AddCutPCMCalo("00052013","00200009327000008250400000","1111100017032230000","0163103000000010");
    cuts.AddCutPCMCalo("00085013","00200009327000008250400000","1111100017032230000","0163103000000010");
    cuts.AddCutPCMCalo("00083013","00200009327000008250400000","1111100017032230000","0163103000000010");
  } else if (trainConfig == 402){ // EMCAL clusters
    cuts.AddCutPCMCalo("00052013","00200009327000008250400000","1111100067032230000","0163103000000010");
    cuts.AddCutPCMCalo("00085013","00200009327000008250400000","1111100067032230000","0163103000000010");
    cuts.AddCutPCMCalo("00083013","00200009327000008250400000","1111100067032230000","0163103000000010");
  } else if (trainConfig == 403){ // EMCAL clusters, PCM-EMC NL // -50ns, 30ns timing cut
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111111067032230000","0163103000000010");
    cuts.AddCutPCMCalo("00052013","00200009327000008250400000","1111111067032230000","0163103000000010");
    cuts.AddCutPCMCalo("00085013","00200009327000008250400000","1111111067032230000","0163103000000010");
    cuts.AddCutPCMCalo("00083013","00200009327000008250400000","1111111067032230000","0163103000000010");
  } else if (trainConfig == 404){ // EMCAL clusters
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111100017032230000","0163103b00000010");
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","1111100067032230000","0163103b00000010");
  } else if (trainConfig == 405){ // EMCAL clusters
    cuts.AddCutPCMCalo("00052013","00200009327000008250400000","1111100017032230000","0163103b00000010");
    cuts.AddCutPCMCalo("00085013","00200009327000008250400000","1111100017032230000","0163103b00000010");
    cuts.AddCutPCMCalo("00083013","00200009327000008250400000","1111100017032230000","0163103b00000010");
  } else if (trainConfig == 406){ // EMCAL clusters
    cuts.AddCutPCMCalo("00052013","00200009327000008250400000","1111100067032230000","0163103b00000010");
    cuts.AddCutPCMCalo("00085013","00200009327000008250400000","1111100067032230000","0163103b00000010");
    cuts.AddCutPCMCalo("00083013","00200009327000008250400000","1111100067032230000","0163103b00000010");

  // *********************************************************************************************************
  // 13 TeV  DMC configurations
  // *********************************************************************************************************
  } else if (trainConfig == 450){ // DCAL clusters
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","3885500081041220000","0163103000000010"); //
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","3885500011041220000","0163103000000010"); //
  } else if (trainConfig == 451){ // DCAL clusters // -1000ns, 1000ns timing cut // no NL
    cuts.AddCutPCMCalo("00055113","00200009327000008250400000","3885500011041220000","0163103000000010");
    cuts.AddCutPCMCalo("00089113","00200009327000008250400000","3885500011041220000","0163103000000010");
    cuts.AddCutPCMCalo("0008b113","00200009327000008250400000","3885500011041220000","0163103000000010");
  } else if (trainConfig == 452){ // DCAL clusters // -50ns, 30ns timing cut // no NL
    cuts.AddCutPCMCalo("00055113","00200009327000008250400000","3885500081041220000","0163103000000010");
    cuts.AddCutPCMCalo("00089113","00200009327000008250400000","3885500011041220000","0163103000000010");
    cuts.AddCutPCMCalo("0008b113","00200009327000008250400000","3885500081041220000","0163103000000010");
  } else if (trainConfig == 453){ // DCAL clusters
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","3885500081041220000","0163103b00000010"); //
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","3885500011041220000","0163103b00000010"); //
  } else if (trainConfig == 454){ // DCAL clusters // -1000ns, 1000ns timing cut // no NL
    cuts.AddCutPCMCalo("00055113","00200009327000008250400000","3885500011041220000","0163103b00000010");
    cuts.AddCutPCMCalo("00089113","00200009327000008250400000","3885500011041220000","0163103b00000010");
    cuts.AddCutPCMCalo("0008b113","00200009327000008250400000","3885500011041220000","0163103b00000010");
  } else if (trainConfig == 455){ // DCAL clusters // -50ns, 30ns timing cut // no NL
    cuts.AddCutPCMCalo("00055113","00200009327000008250400000","3885500081041220000","0163103b00000010");
    cuts.AddCutPCMCalo("00089113","00200009327000008250400000","3885500011041220000","0163103b00000010");
    cuts.AddCutPCMCalo("0008b113","00200009327000008250400000","3885500081041220000","0163103b00000010");

  // *****************************************************************************************************
  // 2.76 TeV pp Run1 - PHOS configurations
  // *****************************************************************************************************
  } else if (trainConfig == 500) { //PHOS clusters
    cuts.AddCutPCMCalo("00003113","00200009327000008250400000","2444400041033200000","0163103000000010");
    cuts.AddCutPCMCalo("00061113","00200009327000008250400000","2444400041033200000","0163103000000010");
    // LHC13g & LHC12x
  } else if (trainConfig == 501) { //PHOS clusters
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","2444400041033200000","0163103000000010");
    cuts.AddCutPCMCalo("00062113","00200009327000008250400000","2444400041033200000","0163103000000010");
  } else if (trainConfig == 502) { //PHOS clusters
    cuts.AddCutPCMCalo("00003113","00200009327000008250400000","2444400041033200000","0163103b00000010");
    cuts.AddCutPCMCalo("00061113","00200009327000008250400000","2444400041033200000","0163103b00000010");
    // LHC13g & LHC12x
  } else if (trainConfig == 503) { //PHOS clusters
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","2444400041033200000","0163103b00000010");
    cuts.AddCutPCMCalo("00062113","00200009327000008250400000","2444400041033200000","0163103b00000010");

  // *****************************************************************************************************
  // 8 TeV pp Run1 - PHOS configurations
  // *****************************************************************************************************
  } else if (trainConfig == 600){ // PHOS clusters
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","2444400072023200000","0163103000000010"); // 600 MeV cluster min energy, t<|30|ns
    cuts.AddCutPCMCalo("00062113","00200009327000008250400000","2444400072023200000","0163103000000010"); // 600 MeV cluster min energy, t<|30|ns
  } else if (trainConfig == 601){ // PHOS clusters
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","2444400072023200000","0163103b00000010"); // 600 MeV cluster min energy, t<|30|ns
    cuts.AddCutPCMCalo("00062113","00200009327000008250400000","2444400072023200000","0163103b00000010"); // 600 MeV cluster min energy, t<|30|ns

  // *****************************************************************************************************
  // 7 TeV pp Run1 - PHOS configurations
  // *****************************************************************************************************
  } else if (trainConfig == 700){
    cuts.AddCutPCMCalo("00000113","00200009327000008250400000","2444400000013300000","0163103000000010"); // QA
  } else if (trainConfig == 701){
    cuts.AddCutPCMCalo("00000113","00200009327000008250400000","2444400040013300000","0163103000000010"); // 100ns timing cut, no track matching
    cuts.AddCutPCMCalo("00000113","00200009327000008250400000","2444400043013300000","0163103000000010"); // 100ns timing cut
  } else if (trainConfig == 702){ // train config for bad channels and NonLin Variation
    cuts.AddCutPCMCalo("00000113","00200009327000008250400000","2444411000013300000","0163103000000010"); // own constant ConvCalo NonLin
    cuts.AddCutPCMCalo("00000113","00200009327000008250400000","2444421000013300000","0163103000000010"); // ex
  } else if (trainConfig == 703){
    cuts.AddCutPCMCalo("00000113","00200009327000008250400000","2444400000013300000","0163103b00000010"); // QA
  } else if (trainConfig == 704){
    cuts.AddCutPCMCalo("00000113","00200009327000008250400000","2444400040013300000","0163103b00000010"); // 100ns timing cut, no track matching
    cuts.AddCutPCMCalo("00000113","00200009327000008250400000","2444400043013300000","0163103b00000010"); // 100ns timing cut
  } else if (trainConfig == 705){ // train config for bad channels and NonLin Variation
    cuts.AddCutPCMCalo("00000113","00200009327000008250400000","2444411000013300000","0163103b00000010"); // own constant ConvCalo NonLin
    cuts.AddCutPCMCalo("00000113","00200009327000008250400000","2444421000013300000","0163103b00000010"); // ex

  // *****************************************************************************************************
  // 5 TeV pp Run2 - PHOS configurations
  // *****************************************************************************************************
  } else if (trainConfig == 800){ // PHOS clusters with larger acceptance
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","2446600040013300000","0163103000000010"); // INT7
    cuts.AddCutPCMCalo("00062113","00200009327000008250400000","2446600040013300000","0163103000000010"); // PHI7
  } else if (trainConfig == 801){ // PHOS clusters with larger acceptance
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","2446600040013300000","0163103b00000010"); // INT7
    cuts.AddCutPCMCalo("00062113","00200009327000008250400000","2446600040013300000","0163103b00000010"); // PHI7

  // *****************************************************************************************************
  // 13 TeV pp Run2 - PHOS configurations
  // *****************************************************************************************************
  } else if (trainConfig == 900){ // PHOS clusters with larger acceptance
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","2446600000013300000","0163103000000010"); // QA
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","2446600040013300000","0163103000000010"); // QA, 100ns timing
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","2446600043013300000","0163103000000010"); // QA, 100ns timing, TM on with default EMC params
  } else if (trainConfig == 901){ // PHI7
    cuts.AddCutPCMCalo("00062113","00200009327000008250400000","2446600000013300000","0163103000000010"); // QA
    cuts.AddCutPCMCalo("00062113","00200009327000008250400000","2446600040013300000","0163103000000010"); // QA, 100ns timing
    cuts.AddCutPCMCalo("00062113","00200009327000008250400000","2446600043013300000","0163103000000010"); // QA, 100ns timing, TM on with default EMC params
  } else if (trainConfig == 902){ // PHOS clusters with larger acceptance
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","2446600000013300000","0163103b00000010"); // QA
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","2446600040013300000","0163103b00000010"); // QA, 100ns timing
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","2446600043013300000","0163103b00000010"); // QA, 100ns timing, TM on with default EMC params
  } else if (trainConfig == 903){ // PHI7
    cuts.AddCutPCMCalo("00062113","00200009327000008250400000","2446600000013300000","0163103b00000010"); // QA
    cuts.AddCutPCMCalo("00062113","00200009327000008250400000","2446600040013300000","0163103b00000010"); // QA, 100ns timing
    cuts.AddCutPCMCalo("00062113","00200009327000008250400000","2446600043013300000","0163103b00000010"); // QA, 100ns timing, TM on with default EMC params
  } else if (trainConfig == 904){ // for eta prime
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","2446600043013300000","01631030000000d0"); // INT7
    cuts.AddCutPCMCalo("00062113","00200009327000008250400000","2446600043013300000","01631030000000d0"); // PHI7


  // *********************************************************************************************************
  // 13 TeV  pp Run2 - EDC configurations
  // *********************************************************************************************************
  } else if (trainConfig == 1000){ // EMCAL + DCal clusters
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","4117911067032230000","01631030000000d0"); //INT7
    cuts.AddCutPCMCalo("0008e113","00200009327000008250400000","4117911067032230000","01631030000000d0"); //EG2
    cuts.AddCutPCMCalo("0008d113","00200009327000008250400000","4117911067032230000","01631030000000d0"); //EG1

  } else {
    Error(Form("HeavyNeutralMesonToGG_%i",trainConfig), "wrong trainConfig variable no cuts have been specified for the configuration");
    return;
  }

  if(!cuts.AreValid()){
    cout << "\n\n****************************************************" << endl;
    cout << "ERROR: No valid cuts stored in CutHandlerHeavyMesonConvCalo! Returning..." << endl;
    cout << "****************************************************\n\n" << endl;
    return;
  }

  TList *EventCutList   = new TList();
  TList *ConvCutList    = new TList();
  TList *ClusterCutList = new TList();
  TList *MesonCutList   = new TList();

  Int_t numberOfCuts = cuts.GetNCuts();

  TList *HeaderList     = new TList();
  if (doWeightingPart==1) {
    TObjString *Header1 = new TObjString("pi0_1");
    HeaderList->Add(Header1);
  }
  if (doWeightingPart==2){
    TObjString *Header3 = new TObjString("eta_2");
    HeaderList->Add(Header3);
  }
  if (doWeightingPart==3) {
    TObjString *Header1 = new TObjString("pi0_1");
    HeaderList->Add(Header1);
    TObjString *Header3 = new TObjString("eta_2");
    HeaderList->Add(Header3);
  }

  EventCutList->SetOwner(kTRUE);
  AliConvEventCuts **analysisEventCuts        = new AliConvEventCuts*[numberOfCuts];
  ConvCutList->SetOwner(kTRUE);
  AliConversionPhotonCuts **analysisCuts      = new AliConversionPhotonCuts*[numberOfCuts];
  ClusterCutList->SetOwner(kTRUE);
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

    analysisCuts[i] = new AliConversionPhotonCuts();
    analysisCuts[i]->SetV0ReaderName(V0ReaderName);
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
    if (enableLightOutput > 0) analysisCuts[i]->SetLightOutput(kTRUE);
    analysisCuts[i]->InitializeCutsFromCutString((cuts.GetPhotonCut(i)).Data());
    analysisCuts[i]->SetIsHeavyIon(isHeavyIon);
    ConvCutList->Add(analysisCuts[i]);
    analysisCuts[i]->SetFillCutHistograms("",kFALSE);

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
    if (enableLightOutput > 0) analysisMesonCuts[i]->SetLightOutput(kTRUE);
    analysisMesonCuts[i]->SetRunningMode(2);
    analysisMesonCuts[i]->InitializeCutsFromCutString((cuts.GetMesonCut(i)).Data());
    MesonCutList->Add(analysisMesonCuts[i]);
    analysisMesonCuts[i]->SetFillCutHistograms("");
    analysisEventCuts[i]->SetAcceptedHeader(HeaderList);
  }

  task->SetEventCutList(numberOfCuts,EventCutList);
  task->SetConversionCutList(numberOfCuts,ConvCutList);
  task->SetCaloCutList(numberOfCuts,ClusterCutList);
  task->SetMesonCutList(numberOfCuts,MesonCutList);
  task->SetMoveParticleAccordingToVertex(kTRUE);
  task->SetCorrectionTaskSetting(corrTaskSetting);
  task->SetMesonType(selectedMeson);
  task->SetDoMesonQA(enableQAMesonTask); //Attention new switch for Pi0 QA
  task->SetDoPhotonQA(enableQAPhotonTask);  //Attention new switch small for Photon QA
  task->SetDoClusterQA(1);  //Attention new switch small for Cluster QA
  task->SetUseTHnSparse(enableTHnSparse);
  task->SetEnableSortingOfMCClusLabels(enableSortingMCLabels);
  if(enableExtMatchAndQA > 1){ task->SetPlotHistsExtQA(kTRUE);}

  //connect containers
  AliAnalysisDataContainer *coutput =
  mgr->CreateContainer(!(corrTaskSetting.CompareTo("")) ? Form("HeavyNeutralMesonToGG_%i_%i_%i", mesonRecoMode, selectedMeson, trainConfig)
                                                        :  Form("HeavyNeutralMesonToGG_%i_%i_%i_%s", mesonRecoMode, selectedMeson, trainConfig, corrTaskSetting.Data()), TList::Class(),
                        AliAnalysisManager::kOutputContainer,Form("HeavyNeutralMesonToGG_%i_%i.root",mesonRecoMode,trainConfig) );

  mgr->AddTask(task);
  mgr->ConnectInput(task,0,cinput);
  mgr->ConnectOutput(task,1,coutput);

  return;

}
