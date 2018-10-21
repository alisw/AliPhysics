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
//($ALIPHYSICS/PWGGA/GammaConv/AliAnalysisTaskGammaConvCalo.cxx) for
//PbPb together with all supporting classes
//***************************************************************************************

//***************************************************************************************
//main function
//***************************************************************************************
void AddTask_GammaConvCalo_PbPb(
  Int_t     trainConfig                   = 1,        // change different set of cuts
  Int_t     isMC                          = 0,        // run MC
  TString   photonCutNumberV0Reader       = "",       // 00000008400000000100000000 nom. B, 00000088400000000100000000 low B
  TString   periodNameV0Reader            = "",
  TString   periodName                    = "",        // period name
  // general setting for task
  Int_t     enableQAMesonTask             = 0,        // enable QA in AliAnalysisTaskGammaConvV1
  Int_t     enableQAPhotonTask            = 0,        // enable additional QA task
  Int_t     enableExtMatchAndQA           = 0,                            // disabled (0), extMatch (1), extQA_noCellQA (2), extMatch+extQA_noCellQA (3), extQA+cellQA (4), extMatch+extQA+cellQA (5)
  Int_t     enableLightOutput             = 0,   // switch to run light output (only essential histograms for afterburner)
  Bool_t    enableTHnSparse               = kFALSE,   // switch on THNsparse
  Bool_t    enableTriggerMimicking        = kFALSE,   // enable trigger mimicking
  Bool_t    enableTriggerOverlapRej       = kFALSE,   // enable trigger overlap rejection
  Float_t   maxFacPtHard                  = 3.,       // maximum factor between hardest jet and ptHard generated
  Int_t     debugLevel                    = 0,        // introducing debug levels for grid running
  // settings for weights
  // FPTW:fileNamePtWeights, FMUW:fileNameMultWeights, FCEF:fileNameCentFlattening, separate with ;
  TString   fileNameExternalInputs        = "",
  Bool_t    doWeightingPart               = kFALSE,        // enable Weighting
  TString   generatorName                 = "DPMJET", // generator Name
  Bool_t    enableMultiplicityWeighting   = kFALSE,   //
  TString   periodNameAnchor              = "",       //
  Bool_t    enableFlattening              = kFALSE,                 // switch on centrality flattening for LHC11h
  // special settings
  Bool_t    enableSortingMCLabels         = kTRUE,    // enable sorting for MC cluster labels
  Bool_t    enableTreeConvGammaShape      = kFALSE,                 // enable additional tree for conversion properties for clusters
  Bool_t    doPrimaryTrackMatching        = kTRUE,                  // enable basic track matching for all primary tracks to cluster
  Int_t     headerSelectionInt            = 0,                      // 1 pi0 header, 2 eta header, 3 both (only for "named" boxes)
  Bool_t    enableHeaderOverlap           = kTRUE,                // allow overlapping header for the clusters
  // subwagon config
  TString   additionalTrainConfig         = "0"       // additional counter for trainconfig
  ) {

  AliCutHandlerPCM cuts;

  TString fileNamePtWeights     = cuts.GetSpecialFileNameFromString (fileNameExternalInputs, "FPTW:");
  TString fileNameMultWeights   = cuts.GetSpecialFileNameFromString (fileNameExternalInputs, "FMUW:");
  TString fileNameCentFlattening= cuts.GetSpecialFileNameFromString (fileNameExternalInputs, "FCEF:");

  TString sAdditionalTrainConfig      = cuts.GetSpecialSettingFromAddConfig(additionalTrainConfig, "", "");
  if (sAdditionalTrainConfig.Atoi() > 0){
    trainConfig = trainConfig + sAdditionalTrainConfig.Atoi();
    cout << "INFO: AddTask_GammaConvV1_pPb running additionalTrainConfig '" << sAdditionalTrainConfig.Atoi() << "', train config: '" << trainConfig << "'" << endl;
  }
  TString corrTaskSetting             = cuts.GetSpecialSettingFromAddConfig(additionalTrainConfig, "CF", "");
  if(corrTaskSetting.CompareTo(""))
    cout << "corrTaskSetting: " << corrTaskSetting.Data() << endl;

  Int_t trackMatcherRunningMode = 0; // CaloTrackMatcher running mode
  TString strTrackMatcherRunningMode  = cuts.GetSpecialSettingFromAddConfig(additionalTrainConfig, "TM", "");
  if(additionalTrainConfig.Contains("TM"))
    trackMatcherRunningMode = strTrackMatcherRunningMode.Atoi();

  Bool_t doTreeEOverP = kFALSE; // switch to produce EOverP tree
  TString strdoTreeEOverP             = cuts.GetSpecialSettingFromAddConfig(additionalTrainConfig, "EPCLUSTree", "");
  if(strdoTreeEOverP.Atoi()==1)
    doTreeEOverP = kTRUE;

  Bool_t doTreeClusterShowerShape = kFALSE; // switch to produce EOverP tree
  TString strdoTreeClusterShowerShape             = cuts.GetSpecialSettingFromAddConfig(additionalTrainConfig, "INVMASSCLUSTree", "");
  if(strdoTreeClusterShowerShape.Atoi()==1)
    doTreeClusterShowerShape = kTRUE;

  TH1S* histoAcc = 0x0;         // histo for modified acceptance
  TString strModifiedAcc              = cuts.GetSpecialSettingFromAddConfig(additionalTrainConfig, "MODIFYACC", "");
  if(strModifiedAcc.Contains("MODIFYACC")){
    cout << "INFO: AddTask_GammaCalo_pp activating 'MODIFYacc'" << endl;
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
  TString strLocalDebugFlag              = cuts.GetSpecialSettingFromAddConfig(additionalTrainConfig, "LOCALDEBUGFLAG", "");
  if(strLocalDebugFlag.Atoi()>0)
    localDebugFlag = strLocalDebugFlag.Atoi();

  Int_t isHeavyIon = 1;

  // ================== GetAnalysisManager ===============================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error(Form("AddTask_GammaConvCalo_PbPb_%i",trainConfig), "No analysis manager found.");
    return ;
  }

  // ================== GetInputEventHandler =============================
  AliVEventHandler *inputHandler=mgr->GetInputEventHandler();

  //=========  Set Cutnumber for V0Reader ================================
  TString cutnumberPhoton = photonCutNumberV0Reader.Data();
  TString cutnumberEvent = "10000003";
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
  AliAnalysisTaskGammaConvCalo *task=NULL;
  task= new AliAnalysisTaskGammaConvCalo(Form("GammaConvCalo_%i",trainConfig));
  task->SetIsHeavyIon(isHeavyIon);
  task->SetIsMC(isMC);
  task->SetV0ReaderName(V0ReaderName);
  task->SetCorrectionTaskSetting(corrTaskSetting);
  if (enableLightOutput > 1) task->SetLightOutput(kTRUE);
  task->SetDoPrimaryTrackMatching(doPrimaryTrackMatching);
  task->SetTrackMatcherRunningMode(trackMatcherRunningMode);

  // cluster cuts
  // 0 "ClusterType",  1 "EtaMin", 2 "EtaMax", 3 "PhiMin", 4 "PhiMax", 5 "DistanceToBadChannel", 6 "Timing", 7 "TrackMatching", 8 "ExoticCell",
  // 9 "MinEnergy", 10 "MinNCells", 11 "MinM02", 12 "MaxM02", 13 "MinMaxM20", 14 "RecConv", 15 "MaximumDispersion", 16 "NLM"

  //****************************************************************************************************
  // EMCal 2.76 TeV Pb-Pb  LHC11h & LHC10h MB
  //****************************************************************************************************
  if (trainConfig == 1){ // EMCAL clusters
    cuts.AddCutPCMCalo("60100013","00200009297002008250400000","1111100053032230000","0163103100000010"); // 0-5%
    cuts.AddCutPCMCalo("61200013","00200009297002008250400000","1111100053032230000","0163103100000010"); // 5-10%
    cuts.AddCutPCMCalo("50100013","00200009297002008250400000","1111100053032230000","0163103100000010"); // 0-10%
    cuts.AddCutPCMCalo("52400013","00200009297002008250400000","1111100053032230000","0163103100000010"); // 20-40%
    cuts.AddCutPCMCalo("52500013","00200009297002008250400000","1111100053032230000","0163103100000010"); // 20-50%
  } else if (trainConfig == 2){ // EMCAL clusters no timing
    cuts.AddCutPCMCalo("60100013","00200009297002008250400000","1111100003032230000","0163103100000010"); // 0-5%
    cuts.AddCutPCMCalo("61200013","00200009297002008250400000","1111100003032230000","0163103100000010"); // 5-10%
    cuts.AddCutPCMCalo("50100013","00200009297002008250400000","1111100003032230000","0163103100000010"); // 0-10%
    cuts.AddCutPCMCalo("52400013","00200009297002008250400000","1111100003032230000","0163103100000010"); // 20-40%
    cuts.AddCutPCMCalo("52500013","00200009297002008250400000","1111100003032230000","0163103100000010"); // 20-50%
  } else if (trainConfig == 3){ // EMCAL clusters
    cuts.AddCutPCMCalo("60100013","00200009297002008250400000","1111100053032230000","0163103100000010"); // 0-5%
    cuts.AddCutPCMCalo("61200013","00200009297002008250400000","1111100053032230000","0163103100000010"); // 5-10%
    cuts.AddCutPCMCalo("50100013","00200009297002008250400000","1111100053032230000","0163103100000010"); // 0-10%
    cuts.AddCutPCMCalo("51200013","00200009297002008250400000","1111100053032230000","0163103100000010"); // 10-20%
    cuts.AddCutPCMCalo("52400013","00200009297002008250400000","1111100053032230000","0163103100000010"); // 20-40%
  } else if (trainConfig == 4){ // EMCAL clusters no timing
    cuts.AddCutPCMCalo("60100013","00200009297002008250400000","1111100003032230000","0163103100000010"); // 0-5%
    cuts.AddCutPCMCalo("61200013","00200009297002008250400000","1111100003032230000","0163103100000010"); // 5-10%
    cuts.AddCutPCMCalo("50100013","00200009297002008250400000","1111100003032230000","0163103100000010"); // 0-10%
    cuts.AddCutPCMCalo("51200013","00200009297002008250400000","1111100003032230000","0163103100000010"); // 10-20%
    cuts.AddCutPCMCalo("52400013","00200009297002008250400000","1111100003032230000","0163103100000010"); // 20-40%
  } else if (trainConfig == 5){ // EMCAL clusters
    cuts.AddCutPCMCalo("54600013","00200009297002008250400000","1111100053032230000","0163103100000010"); // 40-60%
    cuts.AddCutPCMCalo("56800013","00200009297002008250400000","1111100053032230000","0163103100000010"); // 60-80%
    cuts.AddCutPCMCalo("52600013","00200009297002008250400000","1111100053032230000","0163103100000010"); // 20-60%
    cuts.AddCutPCMCalo("54800013","00200009297002008250400000","1111100053032230000","0163103100000010"); // 40-80%
    cuts.AddCutPCMCalo("52500013","00200009297002008250400000","1111100053032230000","0163103100000010"); // 20-50%
  } else if (trainConfig == 6){ // EMCAL clusters no timing
    cuts.AddCutPCMCalo("54600013","00200009297002008250400000","1111100003032230000","0163103100000010"); // 40-60%
    cuts.AddCutPCMCalo("56800013","00200009297002008250400000","1111100003032230000","0163103100000010"); // 60-80%
    cuts.AddCutPCMCalo("52600013","00200009297002008250400000","1111100003032230000","0163103100000010"); // 20-60%
    cuts.AddCutPCMCalo("54800013","00200009297002008250400000","1111100003032230000","0163103100000010"); // 40-80%
    cuts.AddCutPCMCalo("52500013","00200009297002008250400000","1111100003032230000","0163103100000010"); // 20-50%

  } else if (trainConfig == 7){ // EMCAL clusters
    cuts.AddCutPCMCalo("60100013","00200009297002008250400000","1111100053032230000","0163103100000010"); // 0-5%
    cuts.AddCutPCMCalo("61200013","00200009297002008250400000","1111100053032230000","0163103100000010"); // 5-10%
    cuts.AddCutPCMCalo("50100013","00200009297002008250400000","1111100053032230000","0163103100000010"); // 0-10%
  } else if (trainConfig == 8){ // EMCAL clusters no timing
    cuts.AddCutPCMCalo("60100013","00200009297002008250400000","1111100003032230000","0163103100000010"); // 0-5%
    cuts.AddCutPCMCalo("61200013","00200009297002008250400000","1111100003032230000","0163103100000010"); // 5-10%
    cuts.AddCutPCMCalo("50100013","00200009297002008250400000","1111100003032230000","0163103100000010"); // 0-10%
  } else if (trainConfig == 9){ // EMCAL clusters
    cuts.AddCutPCMCalo("51200013","00200009297002008250400000","1111100053032230000","0163103100000010"); // 10-20%
    cuts.AddCutPCMCalo("52400013","00200009297002008250400000","1111100053032230000","0163103100000010"); // 20-40%
    cuts.AddCutPCMCalo("52500013","00200009297002008250400000","1111100053032230000","0163103100000010"); // 20-50%
  } else if (trainConfig == 10){ // EMCAL clusters no timing
    cuts.AddCutPCMCalo("51200013","00200009297002008250400000","1111100003032230000","0163103100000010"); // 10-20%
    cuts.AddCutPCMCalo("52400013","00200009297002008250400000","1111100003032230000","0163103100000010"); // 20-40%
    cuts.AddCutPCMCalo("52500013","00200009297002008250400000","1111100003032230000","0163103100000010"); // 20-50%
  } else if (trainConfig == 11){ // EMCAL clusters
    cuts.AddCutPCMCalo("54600013","00200009297002008250400000","1111100053032230000","0163103100000010"); // 40-60%
    cuts.AddCutPCMCalo("56800013","00200009297002008250400000","1111100053032230000","0163103100000010"); // 60-80%
  } else if (trainConfig == 12){ // EMCAL clusters no timing
    cuts.AddCutPCMCalo("54600013","00200009297002008250400000","1111100003032230000","0163103100000010"); // 40-60%
    cuts.AddCutPCMCalo("56800013","00200009297002008250400000","1111100003032230000","0163103100000010"); // 60-80%


  } else if (trainConfig == 13){ // EMCAL clusters
    cuts.AddCutPCMCalo("50900013","00200009297002008250400000","1111100053032230000","0163103100000010"); // 0-90%
  } else if (trainConfig == 14){ // EMCAL clusters
    cuts.AddCutPCMCalo("50900013","00200009297002008250400000","1111100003032230000","0163103100000010"); // 0-90%

    // EMC triggers LHC15o
  } else if (trainConfig == 20){ //MB
    cuts.AddCutPCMCalo("10110013","00200009297002008250400000","1111100051032230000","0163103100000050"); //
    cuts.AddCutPCMCalo("11310013","00200009297002008250400000","1111100051032230000","0163103100000050"); //
    cuts.AddCutPCMCalo("13510013","00200009297002008250400000","1111100051032230000","0163103100000050"); //
    cuts.AddCutPCMCalo("15010013","00200009297002008250400000","1111100051032230000","0163103100000050"); //
  } else if (trainConfig == 21){ //EG2
    cuts.AddCutPCMCalo("10185013","00200009297002008250400000","1111100051032230000","0163103100000050"); //
    cuts.AddCutPCMCalo("11385013","00200009297002008250400000","1111100051032230000","0163103100000050"); //
    cuts.AddCutPCMCalo("13585013","00200009297002008250400000","1111100051032230000","0163103100000050"); //
    cuts.AddCutPCMCalo("15085013","00200009297002008250400000","1111100051032230000","0163103100000050"); //
  } else if (trainConfig == 22){ //EG1
    cuts.AddCutPCMCalo("10183013","00200009297002008250400000","1111100051032230000","0163103100000050"); //
    cuts.AddCutPCMCalo("11383013","00200009297002008250400000","1111100051032230000","0163103100000050"); //
    cuts.AddCutPCMCalo("13583013","00200009297002008250400000","1111100051032230000","0163103100000050"); //
    cuts.AddCutPCMCalo("15083013","00200009297002008250400000","1111100051032230000","0163103100000050"); //
  } else if (trainConfig == 23){ //EJ2
    cuts.AddCutPCMCalo("10195013","00200009297002008250400000","1111100051032230000","0163103100000050"); //
    cuts.AddCutPCMCalo("11395013","00200009297002008250400000","1111100051032230000","0163103100000050"); //
    cuts.AddCutPCMCalo("13595013","00200009297002008250400000","1111100051032230000","0163103100000050"); //
    cuts.AddCutPCMCalo("15095013","00200009297002008250400000","1111100051032230000","0163103100000050"); //
  } else if (trainConfig == 24){ //EJ1
    cuts.AddCutPCMCalo("10193013","00200009297002008250400000","1111100051032230000","0163103100000050"); //
    cuts.AddCutPCMCalo("11393013","00200009297002008250400000","1111100051032230000","0163103100000050"); //
    cuts.AddCutPCMCalo("13593013","00200009297002008250400000","1111100051032230000","0163103100000050"); //
    cuts.AddCutPCMCalo("15093013","00200009297002008250400000","1111100051032230000","0163103100000050"); //
    // DCAL triggers LHC15o
  } else if (trainConfig == 25){ //MB
    cuts.AddCutPCMCalo("10110013","00200009297002008250400000","3885500017032220000","0163103100000050"); //
    cuts.AddCutPCMCalo("11310013","00200009297002008250400000","3885500017032220000","0163103100000050"); //
    cuts.AddCutPCMCalo("13510013","00200009297002008250400000","3885500017032220000","0163103100000050"); //
    cuts.AddCutPCMCalo("15010013","00200009297002008250400000","3885500017032220000","0163103100000050"); //
  } else if (trainConfig == 26){ //EG2
    cuts.AddCutPCMCalo("10185013","00200009297002008250400000","3885500017032220000","0163103100000050"); //
    cuts.AddCutPCMCalo("11385013","00200009297002008250400000","3885500017032220000","0163103100000050"); //
    cuts.AddCutPCMCalo("13585013","00200009297002008250400000","3885500017032220000","0163103100000050"); //
    cuts.AddCutPCMCalo("15085013","00200009297002008250400000","3885500017032220000","0163103100000050"); //
  } else if (trainConfig == 27){ //EG1
    cuts.AddCutPCMCalo("10183013","00200009297002008250400000","3885500017032220000","0163103100000050"); //
    cuts.AddCutPCMCalo("11383013","00200009297002008250400000","3885500017032220000","0163103100000050"); //
    cuts.AddCutPCMCalo("13583013","00200009297002008250400000","3885500017032220000","0163103100000050"); //
    cuts.AddCutPCMCalo("15083013","00200009297002008250400000","3885500017032220000","0163103100000050"); //
  } else if (trainConfig == 28){ //EJ2
    cuts.AddCutPCMCalo("10195013","00200009297002008250400000","3885500017032220000","0163103100000050"); //
    cuts.AddCutPCMCalo("11395013","00200009297002008250400000","3885500017032220000","0163103100000050"); //
    cuts.AddCutPCMCalo("13595013","00200009297002008250400000","3885500017032220000","0163103100000050"); //
    cuts.AddCutPCMCalo("15095013","00200009297002008250400000","3885500017032220000","0163103100000050"); //
  } else if (trainConfig == 29){ //EJ1
    cuts.AddCutPCMCalo("10193013","00200009297002008250400000","3885500017032220000","0163103100000050"); //
    cuts.AddCutPCMCalo("11393013","00200009297002008250400000","3885500017032220000","0163103100000050"); //
    cuts.AddCutPCMCalo("13593013","00200009297002008250400000","3885500017032220000","0163103100000050"); //
    cuts.AddCutPCMCalo("15093013","00200009297002008250400000","3885500017032220000","0163103100000050"); //

  //****************************************************************************************************
  // EMCal 2.76TeV Pb-Pb for LHC11h EMC trigger
  //****************************************************************************************************
  } else if (trainConfig == 30){ // EMCAL clusters
    cuts.AddCutPCMCalo("50980013","00200009297002008250400000","1111100053032230000","0163103100000010"); // 0-90%
    cuts.AddCutPCMCalo("50180013","00200009297002008250400000","1111100053032230000","0163103100000010"); // 0-10%
    cuts.AddCutPCMCalo("51280013","00200009297002008250400000","1111100053032230000","0163103100000010"); // 10-20%
    cuts.AddCutPCMCalo("52580013","00200009297002008250400000","1111100053032230000","0163103100000010"); // 20-50%
    cuts.AddCutPCMCalo("55880013","00200009297002008250400000","1111100053032230000","0163103100000010"); // 50-80%
  } else if (trainConfig == 31){ // EMCAL clusters no timing cut
    cuts.AddCutPCMCalo("50980013","00200009297002008250400000","1111100003032230000","0163103100000010"); // 0-90%
    cuts.AddCutPCMCalo("50180013","00200009297002008250400000","1111100003032230000","0163103100000010"); // 0-10%
    cuts.AddCutPCMCalo("51280013","00200009297002008250400000","1111100003032230000","0163103100000010"); // 10-20%
    cuts.AddCutPCMCalo("52580013","00200009297002008250400000","1111100003032230000","0163103100000010"); // 20-50%
    cuts.AddCutPCMCalo("55880013","00200009297002008250400000","1111100003032230000","0163103100000010"); // 50-80%
  } else if (trainConfig == 32){ // EMCAL clusters
    cuts.AddCutPCMCalo("50980013","00200009297002008250400000","1111100053032230000","0163103100000010"); // 0-90%
  } else if (trainConfig == 33){ // EMCAL clusters no timing cut
    cuts.AddCutPCMCalo("50980013","00200009297002008250400000","1111100003032230000","0163103100000010"); // 0-90%
  } else if (trainConfig == 34){ // EMCAL clusters
    cuts.AddCutPCMCalo("50100013","00200009297002008250400000","1111102053032230000","0163103100000010"); // 0-10
    cuts.AddCutPCMCalo("50100013","00200009297002008250400000","1111171053032230000","0163103100000010"); // 0-10
  } else if (trainConfig == 35){ // EMCAL clusters
    cuts.AddCutPCMCalo("50100013","00200009297002008250400000","11111020530a2230000","0163103100000010"); // 0-10 // using cluster cuts of Astrid
    cuts.AddCutPCMCalo("50100013","00200009297002008250400000","11111710530a2230000","0163103100000010"); // 0-10
    cuts.AddCutPCMCalo("50100013","00200009297002008250400000","1111171053032230000","0163103100000010");
  } else if (trainConfig == 36){ // EMCAL clusters
    cuts.AddCutPCMCalo("52500013","00200009297002008250400000","1111102053032230000","0163103100000010"); // 20-50
    cuts.AddCutPCMCalo("52500013","00200009297002008250400000","1111172053032230000","0163103100000010"); // 20-50
  } else if (trainConfig == 37){ // EMCAL clusters
    cuts.AddCutPCMCalo("52500013","00200009297002008250400000","11111020530a2230000","0163103100000010"); // 20-50
    cuts.AddCutPCMCalo("52500013","00200009297002008250400000","11111720530a2230000","0163103100000010"); // 20-50
  } else if (trainConfig == 38){ // EMCAL clusters - added signals
    cuts.AddCutPCMCalo("50100023","00200009297002008250400000","1111171053032230000","0163103100000010"); // 0-10
    cuts.AddCutPCMCalo("50100023","00200009297002008250400000","11111710530a2230000","0163103100000010"); // 0-10
    cuts.AddCutPCMCalo("50100023","00200009297002008250400000","11111710530b2230000","0163103100000010"); // 0-10
  } else if (trainConfig == 39){ // EMCAL clusters - added signals
    cuts.AddCutPCMCalo("50100023","00200009297002008250400000","1111171053032230000","0163103100000010"); // 0-10
    cuts.AddCutPCMCalo("50100023","00200009297002008250400000","11111710530a2230000","0163103100000010"); // 0-10
    cuts.AddCutPCMCalo("50100023","00200009297002008250400000","11111710530b2230000","0163103100000010"); // 0-10
  } else if (trainConfig == 40){ // EMCAL clusters - added signals
    cuts.AddCutPCMCalo("52500023","00200009297002008250400000","1111172053032230000","0163103100000010"); // 20-50
    cuts.AddCutPCMCalo("52500023","00200009297002008250400000","11111720530a2230000","0163103100000010"); // 20-50
    cuts.AddCutPCMCalo("52500023","00200009297002008250400000","11111720530b2230000","0163103100000010"); // 20-50
  } else if (trainConfig == 41){ // EMCAL clusters - added signals
    cuts.AddCutPCMCalo("52500023","00200009297002008250400000","1111172053032230000","0163103100000010"); // 20-50
    cuts.AddCutPCMCalo("52500023","00200009297002008250400000","11111720530a2230000","0163103100000010"); // 20-50
    cuts.AddCutPCMCalo("52500023","00200009297002008250400000","11111720530b2230000","0163103100000010"); // 20-50
  } else if (trainConfig == 42){ // EMCAL clusters all headers
    cuts.AddCutPCMCalo("50100003","00200009297002008250400000","11111020530a2230000","0163103100000010"); // 0-10
    cuts.AddCutPCMCalo("50100003","00200009297002008250400000","11111710530a2230000","0163103100000010"); // 0-10
    cuts.AddCutPCMCalo("50100003","00200009297002008250400000","1111171053032230000","0163103100000010");
  } else if (trainConfig == 43){ // EMCAL clusters added signals pi0 header forseen
    cuts.AddCutPCMCalo("50100023","00200009297002008250400000","11111020530a2230000","0163103100000010"); // 0-10 // using cluster cuts of Astrid
    cuts.AddCutPCMCalo("50100023","00200009297002008250400000","11111710530a2230000","0163103100000010"); // 0-10
    cuts.AddCutPCMCalo("50100023","00200009297002008250400000","1111171053032230000","0163103100000010");
  } else if (trainConfig == 44){ // EMCAL clusters added signals eta header forseen
    cuts.AddCutPCMCalo("50100023","00200009297002008250400000","11111020530a2230000","0163103100000010"); // 0-10 // using cluster cuts of Astrid
    cuts.AddCutPCMCalo("50100023","00200009297002008250400000","11111710530a2230000","0163103100000010"); // 0-10
    cuts.AddCutPCMCalo("50100023","00200009297002008250400000","1111171053032230000","0163103100000010");
  } else if (trainConfig == 45){ // EMCAL clusters added signals other header forseen
    cuts.AddCutPCMCalo("50100023","00200009297002008250400000","11111020530a2230000","0163103100000010"); // 0-10 // using cluster cuts of Astrid
    cuts.AddCutPCMCalo("50100023","00200009297002008250400000","11111710530a2230000","0163103100000010"); // 0-10
    cuts.AddCutPCMCalo("50100023","00200009297002008250400000","1111171053032230000","0163103100000010");
  } else if (trainConfig == 46){ // EMCAL clusters added signals other header forseen
    cuts.AddCutPCMCalo("50100023","00200009297002008250400000","11111020530a2230000","0163103100000010"); // 0-10 // using cluster cuts of Astrid
    cuts.AddCutPCMCalo("50100023","00200009297002008250400000","11111710530a2230000","0163103100000010"); // 0-10
    cuts.AddCutPCMCalo("50100023","00200009297002008250400000","1111171053032230000","0163103100000010");
  } else if (trainConfig == 47){ // EMCAL clusters V0 Cent
    cuts.AddCutPCMCalo("10100013","00200009297002008250400000","11111020530a2230000","0163103100000010"); // 0-10 // using cluster cuts of Astrid
    cuts.AddCutPCMCalo("10100013","00200009297002008250400000","11111710530a2230000","0163103100000010"); // 0-10
    cuts.AddCutPCMCalo("10100013","00200009297002008250400000","1111171053032230000","0163103100000010");

  //****************************************************************************************************
  // PHOS 2.76TeV Pb-Pb LHC10h & LHC11h
  //****************************************************************************************************
  } else if (trainConfig == 101){ // PHOS clusters
    cuts.AddCutPCMCalo("60100013","00200009297002008250400000","2444400042033200000","0163103100000010"); // 0-5%
    cuts.AddCutPCMCalo("61200013","00200009297002008250400000","2444400042033200000","0163103100000010"); // 5-10%
    cuts.AddCutPCMCalo("50100013","00200009297002008250400000","2444400042033200000","0163103100000010"); // 0-10%
    cuts.AddCutPCMCalo("52400013","00200009297002008250400000","2444400042033200000","0163103100000010"); // 20-40%
    cuts.AddCutPCMCalo("52500013","00200009297002008250400000","2444400042033200000","0163103100000010"); // 20-50%
  } else if (trainConfig == 102){ // PHOS clusters
    cuts.AddCutPCMCalo("60100013","00200009297002008250400000","2444400042033200000","0163103100000010"); // 0-5%
    cuts.AddCutPCMCalo("61200013","00200009297002008250400000","2444400042033200000","0163103100000010"); // 5-10%
    cuts.AddCutPCMCalo("50100013","00200009297002008250400000","2444400042033200000","0163103100000010"); // 0-10%
    cuts.AddCutPCMCalo("51200013","00200009297002008250400000","2444400042033200000","0163103100000010"); // 10-20%
    cuts.AddCutPCMCalo("52400013","00200009297002008250400000","2444400042033200000","0163103100000010"); // 20-40%
  } else if (trainConfig == 103){ // PHOS clusters
    cuts.AddCutPCMCalo("54600013","00200009297002008250400000","2444400042033200000","0163103100000010"); // 40-60%
    cuts.AddCutPCMCalo("56800013","00200009297002008250400000","2444400042033200000","0163103100000010"); // 60-80%
    cuts.AddCutPCMCalo("52600013","00200009297002008250400000","2444400042033200000","0163103100000010"); // 20-60%
    cuts.AddCutPCMCalo("54800013","00200009297002008250400000","2444400042033200000","0163103100000010"); // 40-80%
    cuts.AddCutPCMCalo("52500013","00200009297002008250400000","2444400042033200000","0163103100000010"); // 20-50%
  } else if (trainConfig == 104){ // PHOS clusters no timing
    cuts.AddCutPCMCalo("60100013","00200009297002008250400000","2444400002033200000","0163103100000010"); // 0-5%
    cuts.AddCutPCMCalo("61200013","00200009297002008250400000","2444400002033200000","0163103100000010"); // 5-10%
    cuts.AddCutPCMCalo("50100013","00200009297002008250400000","2444400002033200000","0163103100000010"); // 0-10%
    cuts.AddCutPCMCalo("52400013","00200009297002008250400000","2444400002033200000","0163103100000010"); // 20-40%
    cuts.AddCutPCMCalo("52500013","00200009297002008250400000","2444400002033200000","0163103100000010"); // 20-50%
  } else if (trainConfig == 105){ // PHOS clusters no timing
    cuts.AddCutPCMCalo("60100013","00200009297002008250400000","2444400002033200000","0163103100000010"); // 0-5%
    cuts.AddCutPCMCalo("61200013","00200009297002008250400000","2444400002033200000","0163103100000010"); // 5-10%
    cuts.AddCutPCMCalo("50100013","00200009297002008250400000","2444400002033200000","0163103100000010"); // 0-10%
    cuts.AddCutPCMCalo("51200013","00200009297002008250400000","2444400002033200000","0163103100000010"); // 10-20%
    cuts.AddCutPCMCalo("52400013","00200009297002008250400000","2444400002033200000","0163103100000010"); // 20-40%
  } else if (trainConfig == 106){ // PHOS clusters no timing
    cuts.AddCutPCMCalo("54600013","00200009297002008250400000","2444400002033200000","0163103100000010"); // 40-60%
    cuts.AddCutPCMCalo("56800013","00200009297002008250400000","2444400002033200000","0163103100000010"); // 60-80%
    cuts.AddCutPCMCalo("52600013","00200009297002008250400000","2444400002033200000","0163103100000010"); // 20-60%
    cuts.AddCutPCMCalo("54800013","00200009297002008250400000","2444400002033200000","0163103100000010"); // 40-80%
    cuts.AddCutPCMCalo("52500013","00200009297002008250400000","2444400002033200000","0163103100000010"); // 20-50%
  } else if (trainConfig == 107){ // PHOS clusters no timing
    cuts.AddCutPCMCalo("50900013","00200009297002008250400000","2444400002033200000","0163103100000010"); // 0-90%

  //****************************************************************************************************
  // EMCal 5TeV Pb-Pb LHC15o
  //****************************************************************************************************
  } else if (trainConfig == 201){ // EMCAL clusters central
    cuts.AddCutPCMCalo("10110013","00200009327000008250400000","1111100057032230000","0163103100000010"); // 0-10
    cuts.AddCutPCMCalo("30110013","00200009327000008250400000","1111100057032230000","0163103100000010"); // 0-5
    cuts.AddCutPCMCalo("31210013","00200009327000008250400000","1111100057032230000","0163103100000010"); // 5-10
  } else if (trainConfig == 202){ // EMCAL clusters semi-central
    cuts.AddCutPCMCalo("11210013","00200009327000008250400000","1111100057032230000","0163103100000010"); // 10-20
    cuts.AddCutPCMCalo("12310013","00200009327000008250400000","1111100057032230000","0163103100000010"); // 20-30
    cuts.AddCutPCMCalo("13410013","00200009327000008250400000","1111100057032230000","0163103100000010"); // 30-40
    cuts.AddCutPCMCalo("12410013","00200009327000008250400000","1111100057032230000","0163103100000010"); // 20-40
  } else if (trainConfig == 203){ // EMCAL clusters semi-central
    cuts.AddCutPCMCalo("14510013","00200009327000008250400000","1111100057032230000","0163103100000010"); // 40-50
    cuts.AddCutPCMCalo("14610013","00200009327000008250400000","1111100057032230000","0163103100000010"); // 40-60
    cuts.AddCutPCMCalo("15610013","00200009327000008250400000","1111100057032230000","0163103100000010"); // 40-60
  } else if (trainConfig == 204){ // EMCAL clusters peripheral
    cuts.AddCutPCMCalo("16710013","00200009327000008250400000","1111100057032230000","0163103100000010"); // 60-70
    cuts.AddCutPCMCalo("17810013","00200009327000008250400000","1111100057032230000","0163103100000010"); // 70-80
    cuts.AddCutPCMCalo("18910013","00200009327000008250400000","1111100057032230000","0163103100000010"); // 80-90
    cuts.AddCutPCMCalo("16810013","00200009327000008250400000","1111100057032230000","0163103100000010"); // 60-80

  } else if (trainConfig == 205){ // EMCAL clusters central add sig
    cuts.AddCutPCMCalo("10110023","00200009327000008250400000","1111100057032230000","0163103100000010"); // 0-10
    cuts.AddCutPCMCalo("30110023","00200009327000008250400000","1111100057032230000","0163103100000010"); // 0-5
    cuts.AddCutPCMCalo("31210023","00200009327000008250400000","1111100057032230000","0163103100000010"); // 5-10
  } else if (trainConfig == 206){ // EMCAL clusters semi-central add sig
    cuts.AddCutPCMCalo("11210023","00200009327000008250400000","1111100057032230000","0163103100000010"); // 10-20
    cuts.AddCutPCMCalo("12310023","00200009327000008250400000","1111100057032230000","0163103100000010"); // 20-30
    cuts.AddCutPCMCalo("13410023","00200009327000008250400000","1111100057032230000","0163103100000010"); // 30-40
    cuts.AddCutPCMCalo("12410023","00200009327000008250400000","1111100057032230000","0163103100000010"); // 20-40
  } else if (trainConfig == 207){ // EMCAL clusters semi-central add sig
    cuts.AddCutPCMCalo("14510023","00200009327000008250400000","1111100057032230000","0163103100000010"); // 40-50
    cuts.AddCutPCMCalo("14610023","00200009327000008250400000","1111100057032230000","0163103100000010"); // 40-60
    cuts.AddCutPCMCalo("15610023","00200009327000008250400000","1111100057032230000","0163103100000010"); // 40-60
  } else if (trainConfig == 208){ // EMCAL clusters peripheral add sig
    cuts.AddCutPCMCalo("16710023","00200009327000008250400000","1111100057032230000","0163103100000010"); // 60-70
    cuts.AddCutPCMCalo("17810023","00200009327000008250400000","1111100057032230000","0163103100000010"); // 70-80
    cuts.AddCutPCMCalo("18910023","00200009327000008250400000","1111100057032230000","0163103100000010"); // 80-90
    cuts.AddCutPCMCalo("16810023","00200009327000008250400000","1111100057032230000","0163103100000010"); // 60-80

  } else if (trainConfig == 209){ // EMCAL clusters central add sig
    cuts.AddCutPCMCalo("10110023","00200009327000008250400000","1111100057032230000","0163103100000010"); // 0-10
    cuts.AddCutPCMCalo("30110023","00200009327000008250400000","1111100057032230000","0163103100000010"); // 0-5
    cuts.AddCutPCMCalo("31210023","00200009327000008250400000","1111100057032230000","0163103100000010"); // 5-10
  } else if (trainConfig == 210){ // EMCAL clusters semi-central add sig
    cuts.AddCutPCMCalo("11210023","00200009327000008250400000","1111100057032230000","0163103100000010"); // 10-20
    cuts.AddCutPCMCalo("12310023","00200009327000008250400000","1111100057032230000","0163103100000010"); // 20-30
    cuts.AddCutPCMCalo("13410023","00200009327000008250400000","1111100057032230000","0163103100000010"); // 30-40
    cuts.AddCutPCMCalo("12410023","00200009327000008250400000","1111100057032230000","0163103100000010"); // 20-40
  } else if (trainConfig == 211){ // EMCAL clusters semi-central add sig
    cuts.AddCutPCMCalo("14510023","00200009327000008250400000","1111100057032230000","0163103100000010"); // 40-50
    cuts.AddCutPCMCalo("14610023","00200009327000008250400000","1111100057032230000","0163103100000010"); // 40-60
    cuts.AddCutPCMCalo("15610023","00200009327000008250400000","1111100057032230000","0163103100000010"); // 40-60
  } else if (trainConfig == 212){ // EMCAL clusters peripheral add sig
    cuts.AddCutPCMCalo("16710023","00200009327000008250400000","1111100057032230000","0163103100000010"); // 60-70
    cuts.AddCutPCMCalo("17810023","00200009327000008250400000","1111100057032230000","0163103100000010"); // 70-80
    cuts.AddCutPCMCalo("18910023","00200009327000008250400000","1111100057032230000","0163103100000010"); // 80-90
    cuts.AddCutPCMCalo("16810023","00200009327000008250400000","1111100057032230000","0163103100000010"); // 60-80

  } else if (trainConfig == 213){ // EMCAL clusters central TB
    cuts.AddCutPCMCalo("10110013","00200009327000008250400000","1111102057032230000","0163103100000010"); // 0-10
    cuts.AddCutPCMCalo("30110013","00200009327000008250400000","1111102057032230000","0163103100000010"); // 0-5
    cuts.AddCutPCMCalo("31210013","00200009327000008250400000","1111102057032230000","0163103100000010"); // 5-10
  } else if (trainConfig == 214){ // EMCAL clusters semi-central TB
    cuts.AddCutPCMCalo("11210013","00200009327000008250400000","1111102057032230000","0163103100000010"); // 10-20
    cuts.AddCutPCMCalo("12310013","00200009327000008250400000","1111102057032230000","0163103100000010"); // 20-30
    cuts.AddCutPCMCalo("13410013","00200009327000008250400000","1111102057032230000","0163103100000010"); // 30-40
    cuts.AddCutPCMCalo("12410013","00200009327000008250400000","1111102057032230000","0163103100000010"); // 20-40
  } else if (trainConfig == 215){ // EMCAL clusters semi-central TB
    cuts.AddCutPCMCalo("14510013","00200009327000008250400000","1111102057032230000","0163103100000010"); // 40-50
    cuts.AddCutPCMCalo("14610013","00200009327000008250400000","1111102057032230000","0163103100000010"); // 40-60
    cuts.AddCutPCMCalo("15610013","00200009327000008250400000","1111102057032230000","0163103100000010"); // 40-60
  } else if (trainConfig == 216){ // EMCAL clusters peripheral TB
    cuts.AddCutPCMCalo("16710013","00200009327000008250400000","1111102057032230000","0163103100000010"); // 60-70
    cuts.AddCutPCMCalo("17810013","00200009327000008250400000","1111102057032230000","0163103100000010"); // 70-80
    cuts.AddCutPCMCalo("18910013","00200009327000008250400000","1111102057032230000","0163103100000010"); // 80-90
    cuts.AddCutPCMCalo("16810013","00200009327000008250400000","1111102057032230000","0163103100000010"); // 60-80

  } else if (trainConfig == 226){ // EMCAL clusters - peripheral centrality selection for PbPb EMCal
    cuts.AddCutPCMCalo("15910013","00200009327000008250400000","11111020530a2230000","0163103100000010"); //
    cuts.AddCutPCMCalo("15910013","00200009327000008250400000","11111870530a2230000","0163103100000010"); //
    cuts.AddCutPCMCalo("15910013","00200009327000008250400000","11111020530b2230000","0163103100000010"); //
    cuts.AddCutPCMCalo("15910013","00200009327000008250400000","11111870530b2230000","0163103100000010"); //

  } else if (trainConfig == 240){ // EMCAL clusters - 0-90% centrality for PbPb EMCal cluster QA
    cuts.AddCutPCMCalo("50910013","00200009327000008250400000","1111100053032230000","0163103100000010"); //  0-90 calo correction cent dep
    cuts.AddCutPCMCalo("50910013","00200009327000008250400000","1111102053032230000","0163103100000010"); //  0-90 calo correction cent dep
    cuts.AddCutPCMCalo("50910013","00200009327000008250400000","1111183053032230000","0163103100000010"); //  0-90 calo correction cent dep
  } else if (trainConfig == 241){ // EMCAL clusters - 0-90% centrality for PbPb EMCal cluster QA
    cuts.AddCutPCMCalo("50910613","00200009327000008250400000","1111100053032230000","0163103100000010"); //  0-90 calo correction cent dep
    cuts.AddCutPCMCalo("50910613","00200009327000008250400000","1111102053032230000","0163103100000010"); //  0-90 calo correction cent dep
    cuts.AddCutPCMCalo("50910613","00200009327000008250400000","1111183053032230000","0163103100000010"); //  0-90 calo correction cent dep
  } else if (trainConfig == 242){ // EMCAL clusters - centrality selection for PbPb EMCal
    cuts.AddCutPCMCalo("50110013","00200009327000008250400000","1111100053032230000","0163103100000010"); //  0-10 calo correction cent dep
    cuts.AddCutPCMCalo("51210013","00200009327000008250400000","1111100053032230000","0163103100000010"); // 10-20 calo correction cent dep
    cuts.AddCutPCMCalo("52510013","00200009327000008250400000","1111100053032230000","0163103100000010"); // 20-50 calo correction cent dep
    cuts.AddCutPCMCalo("55910013","00200009327000008250400000","1111100053032230000","0163103100000010"); // 50-90 calo correction cent dep
  } else if (trainConfig == 243){ // EMCAL clusters - centrality selection for PbPb EMCal
    cuts.AddCutPCMCalo("50110613","00200009327000008250400000","1111100053032230000","0163103100000010"); //  0-10 calo correction cent dep
    cuts.AddCutPCMCalo("51210613","00200009327000008250400000","1111100053032230000","0163103100000010"); // 10-20 calo correction cent dep
    cuts.AddCutPCMCalo("52510613","00200009327000008250400000","1111100053032230000","0163103100000010"); // 20-50 calo correction cent dep
    cuts.AddCutPCMCalo("55910613","00200009327000008250400000","1111100053032230000","0163103100000010"); // 50-90 calo correction cent dep
  } else if (trainConfig == 244){ // EMCAL clusters - centrality selection for PbPb EMCal
    cuts.AddCutPCMCalo("50110613","00200009327000008250400000","1111102053032230000","0163103100000010"); //  0-10 calo correction cent dep
    cuts.AddCutPCMCalo("51210613","00200009327000008250400000","1111102053032230000","0163103100000010"); // 10-20 calo correction cent dep
    cuts.AddCutPCMCalo("52510613","00200009327000008250400000","1111102053032230000","0163103100000010"); // 20-50 calo correction cent dep
    cuts.AddCutPCMCalo("55910613","00200009327000008250400000","1111102053032230000","0163103100000010"); // 50-90 calo correction cent dep
  } else if (trainConfig == 245){ // EMCAL clusters - centrality selection for PbPb EMCal
    cuts.AddCutPCMCalo("50110613","00200009327000008250400000","1111184053032230000","0163103100000010"); //  0-10 calo correction cent dep
    cuts.AddCutPCMCalo("51210613","00200009327000008250400000","1111185053032230000","0163103100000010"); // 10-20 calo correction cent dep
    cuts.AddCutPCMCalo("52510613","00200009327000008250400000","1111186053032230000","0163103100000010"); // 20-50 calo correction cent dep
    cuts.AddCutPCMCalo("55910613","00200009327000008250400000","1111187053032230000","0163103100000010"); // 50-90 calo correction cent dep
  } else if (trainConfig == 246){ // EMCAL clusters - 0-90% centrality for PbPb EMCal cluster QA
    cuts.AddCutPCMCalo("50910613","00200009327000008250400000","1111183050032230000","0163103100000010"); //  0-90 calo correction cent dep
    cuts.AddCutPCMCalo("50910613","00200009327000008250400000","1111183051032230000","0163103100000010"); //  0-90 calo correction cent dep
    cuts.AddCutPCMCalo("50910613","00200009327000008250400000","1111183057032230000","0163103100000010"); //  0-90 calo correction cent dep
  } else if (trainConfig == 247){ // EMCAL clusters - centrality selection for PbPb EMCal
    cuts.AddCutPCMCalo("50110613","00200009327000008250400000","1111184050032230000","0163103100000010"); //  0-10 calo correction cent dep
    cuts.AddCutPCMCalo("51210613","00200009327000008250400000","1111185050032230000","0163103100000010"); // 10-20 calo correction cent dep
    cuts.AddCutPCMCalo("52510613","00200009327000008250400000","1111186050032230000","0163103100000010"); // 20-50 calo correction cent dep
    cuts.AddCutPCMCalo("55910613","00200009327000008250400000","1111187050032230000","0163103100000010"); // 50-90 calo correction cent dep
  } else if (trainConfig == 248){ // EMCAL clusters - centrality selection for PbPb EMCal
    cuts.AddCutPCMCalo("50110613","00200009327000008250400000","1111184051032230000","0163103100000010"); //  0-10 calo correction cent dep
    cuts.AddCutPCMCalo("51210613","00200009327000008250400000","1111185051032230000","0163103100000010"); // 10-20 calo correction cent dep
    cuts.AddCutPCMCalo("52510613","00200009327000008250400000","1111186051032230000","0163103100000010"); // 20-50 calo correction cent dep
    cuts.AddCutPCMCalo("55910613","00200009327000008250400000","1111187051032230000","0163103100000010"); // 50-90 calo correction cent dep
  } else if (trainConfig == 249){ // EMCAL clusters - centrality selection for PbPb EMCal
    cuts.AddCutPCMCalo("50110613","00200009327000008250400000","1111184057032230000","0163103100000010"); //  0-10 calo correction cent dep
    cuts.AddCutPCMCalo("51210613","00200009327000008250400000","1111185057032230000","0163103100000010"); // 10-20 calo correction cent dep
    cuts.AddCutPCMCalo("52510613","00200009327000008250400000","1111186057032230000","0163103100000010"); // 20-50 calo correction cent dep
    cuts.AddCutPCMCalo("55910613","00200009327000008250400000","1111187057032230000","0163103100000010"); // 50-90 calo correction cent dep
  } else if (trainConfig == 250){ // EMCAL clusters - 0-90% centrality for PbPb EMCal cluster QA
    cuts.AddCutPCMCalo("10910013","00200009327000008250400000","1111183051032230000","0163103100000010"); //  0-90 calo correction cent dep
    cuts.AddCutPCMCalo("10910a13","00200009327000008250400000","1111183051032230000","0163103100000010"); //  0-90 calo correction cent dep
  } else if (trainConfig == 251){ // EMCAL clusters - 0-90% centrality for PbPb EMCal cluster QA
    cuts.AddCutPCMCalo("10110013","00200009327000008250400000","1111183051032230000","0163103100000010"); //  0-90 calo correction cent dep
    cuts.AddCutPCMCalo("11210013","00200009327000008250400000","1111183051032230000","0163103100000010"); //  0-90 calo correction cent dep
    cuts.AddCutPCMCalo("12510013","00200009327000008250400000","1111183051032230000","0163103100000010"); //  0-90 calo correction cent dep
    cuts.AddCutPCMCalo("15910013","00200009327000008250400000","1111183051032230000","0163103100000010"); //  0-90 calo correction cent dep
  } else if (trainConfig == 252){ // EMCAL clusters - 0-90% centrality for PbPb EMCal cluster QA
    cuts.AddCutPCMCalo("10110a13","00200009327000008250400000","1111183051032230000","0163103100000010"); //  0-90 calo correction cent dep
    cuts.AddCutPCMCalo("11210a13","00200009327000008250400000","1111183051032230000","0163103100000010"); //  0-90 calo correction cent dep
    cuts.AddCutPCMCalo("12510a13","00200009327000008250400000","1111183051032230000","0163103100000010"); //  0-90 calo correction cent dep
    cuts.AddCutPCMCalo("15910a13","00200009327000008250400000","1111183051032230000","0163103100000010"); //  0-90 calo correction cent dep
  } else if (trainConfig == 253){ // EMCAL clusters - 20180718 - default without corrections
    cuts.AddCutPCMCalo("10110a13","00200009327000008250400000","1111100051032230000","0163103100000010"); //  0-90 calo correction cent dep
    cuts.AddCutPCMCalo("11210a13","00200009327000008250400000","1111100051032230000","0163103100000010"); //  0-90 calo correction cent dep
    cuts.AddCutPCMCalo("12410a13","00200009327000008250400000","1111100051032230000","0163103100000010"); //  0-90 calo correction cent dep
    cuts.AddCutPCMCalo("14610a13","00200009327000008250400000","1111100051032230000","0163103100000010"); //  0-90 calo correction cent dep
    cuts.AddCutPCMCalo("16810a13","00200009327000008250400000","1111100051032230000","0163103100000010"); //  0-90 calo correction cent dep
  } else if (trainConfig == 254){ // EMCAL clusters - 20180718 - default with peripheral corrections
    cuts.AddCutPCMCalo("10110a13","00200009327000008250400000","1111187051032230000","0163103100000010"); //
    cuts.AddCutPCMCalo("11210a13","00200009327000008250400000","1111187051032230000","0163103100000010"); //
    cuts.AddCutPCMCalo("12410a13","00200009327000008250400000","1111187051032230000","0163103100000010"); //
    cuts.AddCutPCMCalo("14610a13","00200009327000008250400000","1111187051032230000","0163103100000010"); //
    cuts.AddCutPCMCalo("16810a13","00200009327000008250400000","1111187051032230000","0163103100000010"); //
  } else if (trainConfig == 255){ // EMCAL clusters - 20180718 - default with peripheral corrections
    cuts.AddCutPCMCalo("12410a13","00200009327000008250400000","1111186051032230000","0163103100000010"); //
    cuts.AddCutPCMCalo("12510a13","00200009327000008250400000","1111186051032230000","0163103100000010"); //
    cuts.AddCutPCMCalo("14610a13","00200009327000008250400000","1111187051032230000","0163103100000010"); //
    cuts.AddCutPCMCalo("16810a13","00200009327000008250400000","1111187051032230000","0163103100000010"); //
  } else if (trainConfig == 256){ // EMCAL clusters - 20180718 - default with corrections
    cuts.AddCutPCMCalo("30110a13","00200009327000008250400000","1111184051032230000","0163103100000010"); //
    cuts.AddCutPCMCalo("31210a13","00200009327000008250400000","1111184051032230000","0163103100000010"); //
    cuts.AddCutPCMCalo("10110a13","00200009327000008250400000","1111184051032230000","0163103100000010"); //
    cuts.AddCutPCMCalo("11210a13","00200009327000008250400000","1111185051032230000","0163103100000010"); //


  } else if (trainConfig == 290){ // EMCAL clusters - correction convcalo f1
    cuts.AddCutPCMCalo("10110013","00200009327000008250400000","1111181003032230000","0163103100000010"); // 0-10
    cuts.AddCutPCMCalo("11210013","00200009327000008250400000","1111181003032230000","0163103100000010"); // 10-20
    cuts.AddCutPCMCalo("12510013","00200009327000008250400000","1111181003032230000","0163103100000010"); // 20-50
    cuts.AddCutPCMCalo("15910013","00200009327000008250400000","1111181003032230000","0163103100000010"); // 50-90
  } else if (trainConfig == 291){ // EMCAL clusters - correction calocalo f2
    cuts.AddCutPCMCalo("10110013","00200009327000008250400000","1111192003032230000","0163103100000010"); // 0-10
    cuts.AddCutPCMCalo("11210013","00200009327000008250400000","1111192003032230000","0163103100000010"); // 10-20
    cuts.AddCutPCMCalo("12510013","00200009327000008250400000","1111192003032230000","0163103100000010"); // 20-50
    cuts.AddCutPCMCalo("15910013","00200009327000008250400000","1111192003032230000","0163103100000010"); // 50-90
  } else if (trainConfig == 292){ // EMCAL clusters - 0-90% centrality for PbPb EMCal cluster QA
    cuts.AddCutPCMCalo("10910013","00200009327000008250400000","1111100003032230000","0163103100000010"); // 0-90
  } else if (trainConfig == 293){ // EMCAL clusters - centrality selection for PbPb EMCal
    cuts.AddCutPCMCalo("10110013","00200009327000008250400000","1111184053032230000","0163103100000010"); //  0-10 calo correction cent dep
    cuts.AddCutPCMCalo("11210013","00200009327000008250400000","1111185053032230000","0163103100000010"); // 10-20 calo correction cent dep
    cuts.AddCutPCMCalo("12510013","00200009327000008250400000","1111186053032230000","0163103100000010"); // 20-50 calo correction cent dep
    cuts.AddCutPCMCalo("15910013","00200009327000008250400000","1111187053032230000","0163103100000010"); // 50-90 calo correction cent dep


  //****************************************************************************************************
  // EMCal 5TeV Xe-Xe LHC17n
  //****************************************************************************************************
  } else if (trainConfig == 300){ // EMCAL clusters - 0-80% centrality
    cuts.AddCutPCMCalo("10810013","00200009327000008250400000","1111100017032230000","0163103100000010"); // 0-80
  } else if (trainConfig == 301){ // EMCAL clusters
    cuts.AddCutPCMCalo("10210013","00200009327000008250400000","1111100017032230000","0163103100000010"); // 0-20
    cuts.AddCutPCMCalo("12410013","00200009327000008250400000","1111100017032230000","0163103100000010"); // 20-40
    cuts.AddCutPCMCalo("10410013","00200009327000008250400000","1111100017032230000","0163103100000010"); // 0-40
    cuts.AddCutPCMCalo("14810013","00200009327000008250400000","1111100017032230000","0163103100000010"); // 40-80
  } else if (trainConfig == 302){ // EMCAL clusters - 0-80% centrality for EMCal cluster QA TB calib
    cuts.AddCutPCMCalo("10810013","00200009327000008250400000","1111102017032230000","0163103100000010"); // 0-80
  } else if (trainConfig == 303){ // EMCAL clusters TB calib
    cuts.AddCutPCMCalo("10210013","00200009327000008250400000","1111102017032230000","0163103100000010"); // 0-20
    cuts.AddCutPCMCalo("12410013","00200009327000008250400000","1111102017032230000","0163103100000010"); // 20-40
    cuts.AddCutPCMCalo("10410013","00200009327000008250400000","1111102017032230000","0163103100000010"); // 0-40
    cuts.AddCutPCMCalo("14810013","00200009327000008250400000","1111102017032230000","0163103100000010"); // 40-80
  } else if (trainConfig == 304){ // EMCAL clusters -  EMCal cluster QA TB calib, PCM min pt 0.02
    cuts.AddCutPCMCalo("10810013","00200089327000008250400000","1111102017032230000","0163103100000010"); // 0-80
    cuts.AddCutPCMCalo("10210013","00200089327000008250400000","1111102017032230000","0163103100000010"); // 0-20
    cuts.AddCutPCMCalo("12410013","00200089327000008250400000","1111102017032230000","0163103100000010"); // 20-40
    cuts.AddCutPCMCalo("10410013","00200089327000008250400000","1111102017032230000","0163103100000010"); // 0-40
    cuts.AddCutPCMCalo("14810013","00200089327000008250400000","1111102017032230000","0163103100000010"); // 40-80

  } else if (trainConfig == 310){ // EMCAL clusters  non lin vars 0-20
    cuts.AddCutPCMCalo("10210013","00200009327000008250400000","1111171017032230000","0163103100000010");
    cuts.AddCutPCMCalo("10210013","00200009327000008250400000","1111172017032230000","0163103100000010");
    cuts.AddCutPCMCalo("10210013","00200009327000008250400000","1111181017032230000","0163103100000010");
    cuts.AddCutPCMCalo("10210013","00200009327000008250400000","1111182017032230000","0163103100000010");
  } else if (trainConfig == 320){ // EMCAL clusters  non lin vars 20-40
    cuts.AddCutPCMCalo("12410013","00200009327000008250400000","1111171007032230000","0163103100000010");
    cuts.AddCutPCMCalo("12410013","00200009327000008250400000","1111172007032230000","0163103100000010");
    cuts.AddCutPCMCalo("12410013","00200009327000008250400000","1111181007032230000","0163103100000010");
    cuts.AddCutPCMCalo("12410013","00200009327000008250400000","1111182007032230000","0163103100000010");
  } else if (trainConfig == 330){ // EMCAL clusters  non lin vars 40-80
    cuts.AddCutPCMCalo("14810013","00200009327000008250400000","1111171007032230000","0163103100000010");
    cuts.AddCutPCMCalo("14810013","00200009327000008250400000","1111172007032230000","0163103100000010");
    cuts.AddCutPCMCalo("14810013","00200009327000008250400000","1111181007032230000","0163103100000010");
    cuts.AddCutPCMCalo("14810013","00200009327000008250400000","1111182007032230000","0163103100000010");
  } else if (trainConfig == 340){ // EMCAL clusters  non lin vars 0-40
    cuts.AddCutPCMCalo("10410013","00200009327000008250400000","1111171007032230000","0163103100000010");
    cuts.AddCutPCMCalo("10410013","00200009327000008250400000","1111172007032230000","0163103100000010");
    cuts.AddCutPCMCalo("10410013","00200009327000008250400000","1111181007032230000","0163103100000010");
    cuts.AddCutPCMCalo("10410013","00200009327000008250400000","1111182007032230000","0163103100000010");
  } else if (trainConfig == 350){ // EMCAL clusters  non lin vars 0-80
    cuts.AddCutPCMCalo("10810013","00200009327000008250400000","1111171007032230000","0163103100000010");
    cuts.AddCutPCMCalo("10810013","00200009327000008250400000","1111172007032230000","0163103100000010");
    cuts.AddCutPCMCalo("10810013","00200009327000008250400000","1111181007032230000","0163103100000010");
    cuts.AddCutPCMCalo("10810013","00200009327000008250400000","1111182007032230000","0163103100000010");

  //****************************************************************************************************
  // PHOS 5TeV Xe-Xe LHC17n
  //****************************************************************************************************
  } else if (trainConfig == 400){ // PHOS clusters - 0-90% centrality
    cuts.AddCutPCMCalo("10810013","00200009327000008250400000","2446600011013200000","0163103100000010"); // 0-80%
    cuts.AddCutPCMCalo("10810013","00200009327000008250400000","2446601011013200000","0163103100000010"); // 0-80%
  } else if (trainConfig == 401) {  // PHOS Cent dep Xe-Xe
    cuts.AddCutPCMCalo("10210013","00200009327000008250400000","2446600011013200000","0163103100000010"); // 0-20%
    cuts.AddCutPCMCalo("12410013","00200009327000008250400000","2446600011013200000","0163103100000010"); // 20-40%
    cuts.AddCutPCMCalo("10410013","00200009327000008250400000","2446600011013200000","0163103100000010"); // 0-40%
    cuts.AddCutPCMCalo("14810013","00200009327000008250400000","2446600011013200000","0163103100000010"); // 40-80%
  } else if (trainConfig == 402) {  // PHOS Cent dep diff timing
    cuts.AddCutPCMCalo("10210013","00200009327000008250400000","2446601011013200000","0163103100000010"); // 0-20%
    cuts.AddCutPCMCalo("12410013","00200009327000008250400000","2446601011013200000","0163103100000010"); // 20-40%
    cuts.AddCutPCMCalo("10410013","00200009327000008250400000","2446601011013200000","0163103100000010"); // 0-40%
    cuts.AddCutPCMCalo("14810013","00200009327000008250400000","2446601011013200000","0163103100000010"); // 40-80%
  } else if (trainConfig == 403) {  // PHOS Cent dep diff timing, min pT PCM 0.02
    cuts.AddCutPCMCalo("10210013","00200089327000008250400000","2446601011013200000","0163103100000010"); // 0-20%
    cuts.AddCutPCMCalo("12410013","00200089327000008250400000","2446601011013200000","0163103100000010"); // 20-40%
    cuts.AddCutPCMCalo("10410013","00200089327000008250400000","2446601011013200000","0163103100000010"); // 0-40%
    cuts.AddCutPCMCalo("14810013","00200089327000008250400000","2446601011013200000","0163103100000010"); // 40-80%
  } else if (trainConfig == 410) {  // PHOS Cluster non lin vars 0-20
    cuts.AddCutPCMCalo("10210013","00200009327000008250400000","2446601011013200000","0163103100000010");
    cuts.AddCutPCMCalo("10210013","00200009327000008250400000","2446671011013200000","0163103100000010");
    cuts.AddCutPCMCalo("10210013","00200009327000008250400000","2446672011013200000","0163103100000010");
    cuts.AddCutPCMCalo("10210013","00200009327000008250400000","2446681011013200000","0163103100000010");
    cuts.AddCutPCMCalo("10210013","00200009327000008250400000","2446682011013200000","0163103100000010");
  } else if (trainConfig == 420) {  // PHOS Cluster non lin vars 20-40
    cuts.AddCutPCMCalo("12410013","00200009327000008250400000","2446601011013200000","0163103100000010");
    cuts.AddCutPCMCalo("12410013","00200009327000008250400000","2446671011013200000","0163103100000010");
    cuts.AddCutPCMCalo("12410013","00200009327000008250400000","2446672011013200000","0163103100000010");
    cuts.AddCutPCMCalo("12410013","00200009327000008250400000","2446681011013200000","0163103100000010");
    cuts.AddCutPCMCalo("12410013","00200009327000008250400000","2446682011013200000","0163103100000010");
  } else if (trainConfig == 430) {  // PHOS Cluster non lin vars 40-80
    cuts.AddCutPCMCalo("14810013","00200009327000008250400000","2446601011013200000","0163103100000010");
    cuts.AddCutPCMCalo("14810013","00200009327000008250400000","2446671011013200000","0163103100000010");
    cuts.AddCutPCMCalo("14810013","00200009327000008250400000","2446672011013200000","0163103100000010");
    cuts.AddCutPCMCalo("14810013","00200009327000008250400000","2446681011013200000","0163103100000010");
    cuts.AddCutPCMCalo("14810013","00200009327000008250400000","2446682011013200000","0163103100000010");
  } else if (trainConfig == 440) {  // PHOS Cluster non lin vars 0-40
    cuts.AddCutPCMCalo("10410013","00200009327000008250400000","2446601011013200000","0163103100000010");
    cuts.AddCutPCMCalo("10410013","00200009327000008250400000","2446671011013200000","0163103100000010");
    cuts.AddCutPCMCalo("10410013","00200009327000008250400000","2446672011013200000","0163103100000010");
    cuts.AddCutPCMCalo("10410013","00200009327000008250400000","2446681011013200000","0163103100000010");
    cuts.AddCutPCMCalo("10410013","00200009327000008250400000","2446682011013200000","0163103100000010");
  } else if (trainConfig == 450) {  // PHOS Cluster non lin vars 0-80
    cuts.AddCutPCMCalo("10810013","00200009327000008250400000","2446601011013200000","0163103100000010");
    cuts.AddCutPCMCalo("10810013","00200009327000008250400000","2446671011013200000","0163103100000010");
    cuts.AddCutPCMCalo("10810013","00200009327000008250400000","2446672011013200000","0163103100000010");
    cuts.AddCutPCMCalo("10810013","00200009327000008250400000","2446681011013200000","0163103100000010");
    cuts.AddCutPCMCalo("10810013","00200009327000008250400000","2446682011013200000","0163103100000010");

  //****************************************************************************************************
  // PHOS 5TeV Pb-Pb LHC15o
  //****************************************************************************************************
  } else if (trainConfig == 601){ // EMCAL clusters central
    cuts.AddCutPCMCalo("10110013","00200009327000008250400000","2446600051013200000","0163103100000010"); // 0-10
    cuts.AddCutPCMCalo("30110013","00200009327000008250400000","2446600051013200000","0163103100000010"); // 0-5
    cuts.AddCutPCMCalo("31210013","00200009327000008250400000","2446600051013200000","0163103100000010"); // 5-10
  } else if (trainConfig == 602){ // EMCAL clusters semi-central
    cuts.AddCutPCMCalo("11210013","00200009327000008250400000","2446600051013200000","0163103100000010"); // 10-20
    cuts.AddCutPCMCalo("12310013","00200009327000008250400000","2446600051013200000","0163103100000010"); // 20-30
    cuts.AddCutPCMCalo("13410013","00200009327000008250400000","2446600051013200000","0163103100000010"); // 30-40
    cuts.AddCutPCMCalo("12410013","00200009327000008250400000","2446600051013200000","0163103100000010"); // 20-40
  } else if (trainConfig == 603){ // EMCAL clusters semi-central
    cuts.AddCutPCMCalo("14510013","00200009327000008250400000","2446600051013200000","0163103100000010"); // 40-50
    cuts.AddCutPCMCalo("14610013","00200009327000008250400000","2446600051013200000","0163103100000010"); // 40-60
    cuts.AddCutPCMCalo("15610013","00200009327000008250400000","2446600051013200000","0163103100000010"); // 40-60
  } else if (trainConfig == 604){ // EMCAL clusters peripheral
    cuts.AddCutPCMCalo("16710013","00200009327000008250400000","2446600051013200000","0163103100000010"); // 60-70
    cuts.AddCutPCMCalo("17810013","00200009327000008250400000","2446600051013200000","0163103100000010"); // 70-80
    cuts.AddCutPCMCalo("18910013","00200009327000008250400000","2446600051013200000","0163103100000010"); // 80-90
    cuts.AddCutPCMCalo("16810013","00200009327000008250400000","2446600051013200000","0163103100000010"); // 60-80
  } else if (trainConfig == 605){ // EMCAL clusters central
    cuts.AddCutPCMCalo("10110013","00200009327000008250400000","2446601051013200000","0163103100000010"); // 0-10
    cuts.AddCutPCMCalo("30110013","00200009327000008250400000","2446601051013200000","0163103100000010"); // 0-5
    cuts.AddCutPCMCalo("31210013","00200009327000008250400000","2446601051013200000","0163103100000010"); // 5-10
  } else if (trainConfig == 606){ // EMCAL clusters semi-central
    cuts.AddCutPCMCalo("11210013","00200009327000008250400000","2446601051013200000","0163103100000010"); // 10-20
    cuts.AddCutPCMCalo("12310013","00200009327000008250400000","2446601051013200000","0163103100000010"); // 20-30
    cuts.AddCutPCMCalo("13410013","00200009327000008250400000","2446601051013200000","0163103100000010"); // 30-40
    cuts.AddCutPCMCalo("12410013","00200009327000008250400000","2446601051013200000","0163103100000010"); // 20-40
  } else if (trainConfig == 607){ // EMCAL clusters semi-central
    cuts.AddCutPCMCalo("14510013","00200009327000008250400000","2446601051013200000","0163103100000010"); // 40-50
    cuts.AddCutPCMCalo("14610013","00200009327000008250400000","2446601051013200000","0163103100000010"); // 40-60
    cuts.AddCutPCMCalo("15610013","00200009327000008250400000","2446601051013200000","0163103100000010"); // 40-60
  } else if (trainConfig == 608){ // EMCAL clusters peripheral
    cuts.AddCutPCMCalo("16710013","00200009327000008250400000","2446601051013200000","0163103100000010"); // 60-70
    cuts.AddCutPCMCalo("17810013","00200009327000008250400000","2446601051013200000","0163103100000010"); // 70-80
    cuts.AddCutPCMCalo("18910013","00200009327000008250400000","2446601051013200000","0163103100000010"); // 80-90
    cuts.AddCutPCMCalo("16810013","00200009327000008250400000","2446601051013200000","0163103100000010"); // 60-80

  } else if (trainConfig == 610){ // PHOS clusters - 20181018
    cuts.AddCutPCMCalo("12410a13","00200009327000008250400000","2446600051013200000","0163103100000010"); //
    cuts.AddCutPCMCalo("12510a13","00200009327000008250400000","2446600051013200000","0163103100000010"); //
    cuts.AddCutPCMCalo("14610a13","00200009327000008250400000","2446600051013200000","0163103100000010"); //
    cuts.AddCutPCMCalo("16810a13","00200009327000008250400000","2446600051013200000","0163103100000010"); //
  } else if (trainConfig == 611){ // PHOS clusters - 20181018
    cuts.AddCutPCMCalo("30110a13","00200009327000008250400000","2446600051013200000","0163103100000010"); //
    cuts.AddCutPCMCalo("31210a13","00200009327000008250400000","2446600051013200000","0163103100000010"); //
    cuts.AddCutPCMCalo("10110a13","00200009327000008250400000","2446600051013200000","0163103100000010"); //
    cuts.AddCutPCMCalo("11210a13","00200009327000008250400000","2446600051013200000","0163103100000010"); //


  //****************************************************************************************************
  // EMCal 5TeV Pb-Pb LHC15o EMC triggers
  //****************************************************************************************************
  } else if (trainConfig == 201){ // EMCAL clusters central
    cuts.AddCutPCMCalo("10183013","00200009327000008250400000","1111100057032230000","0163103100000010"); // 0-10
    cuts.AddCutPCMCalo("30183013","00200009327000008250400000","1111100057032230000","0163103100000010"); // 0-5
    cuts.AddCutPCMCalo("31283013","00200009327000008250400000","1111100057032230000","0163103100000010"); // 5-10
  } else if (trainConfig == 202){ // EMCAL clusters semi-central
    cuts.AddCutPCMCalo("11283013","00200009327000008250400000","1111100057032230000","0163103100000010"); // 10-20
    cuts.AddCutPCMCalo("12383013","00200009327000008250400000","1111100057032230000","0163103100000010"); // 20-30
    cuts.AddCutPCMCalo("13483013","00200009327000008250400000","1111100057032230000","0163103100000010"); // 30-40
    cuts.AddCutPCMCalo("12483013","00200009327000008250400000","1111100057032230000","0163103100000010"); // 20-40
  } else if (trainConfig == 203){ // EMCAL clusters semi-central
    cuts.AddCutPCMCalo("14583013","00200009327000008250400000","1111100057032230000","0163103100000010"); // 40-50
    cuts.AddCutPCMCalo("14683013","00200009327000008250400000","1111100057032230000","0163103100000010"); // 40-60
    cuts.AddCutPCMCalo("15683013","00200009327000008250400000","1111100057032230000","0163103100000010"); // 40-60
  } else if (trainConfig == 204){ // EMCAL clusters peripheral
    cuts.AddCutPCMCalo("16783013","00200009327000008250400000","1111100057032230000","0163103100000010"); // 60-70
    cuts.AddCutPCMCalo("17883013","00200009327000008250400000","1111100057032230000","0163103100000010"); // 70-80
    cuts.AddCutPCMCalo("18983013","00200009327000008250400000","1111100057032230000","0163103100000010"); // 80-90
    cuts.AddCutPCMCalo("16883013","00200009327000008250400000","1111100057032230000","0163103100000010"); // 60-80


  } else {
    Error(Form("GammaConvCalo_%i",trainConfig), "wrong trainConfig variable no cuts have been specified for the configuration");
    return;
  }

  if(!cuts.AreValid()){
    cout << "\n\n****************************************************" << endl;
    cout << "ERROR: No valid cuts stored in CutHandlerConvCalo! Returning..." << endl;
    cout << "****************************************************\n\n" << endl;
    return;
  }

  Int_t numberOfCuts = cuts.GetNCuts();

  TList *EventCutList     = new TList();
  TList *ConvCutList      = new TList();
  TList *ClusterCutList   = new TList();
  TList *MesonCutList     = new TList();

  TList *HeaderList       = new TList();
  if (periodName.CompareTo("LHC13d2")==0){
    TObjString *Header1 = new TObjString("pi0_1");
    HeaderList->Add(Header1);
  //    TObjString *Header3 = new TObjString("eta_2");
  //    HeaderList->Add(Header3);

  } else if (periodName.CompareTo("LHC12a17x_fix")==0){
    TObjString *Header1 = new TObjString("PARAM");
    HeaderList->Add(Header1);
  } else if (periodName.CompareTo("LHC14a1a")==0){
    if (headerSelectionInt == 1){
      TObjString *Header1 = new TObjString("pi0_1");
      HeaderList->Add(Header1);
    } else if (headerSelectionInt == 2){
      TObjString *Header1 = new TObjString("eta_2");
      HeaderList->Add(Header1);
    } else if (headerSelectionInt == 3){
      TString nameHeaders[2]    = { "pi0_1", "eta_2" };
      for (Int_t iHead = 0; iHead < 2; iHead++ ){
        TObjString *Header = new TObjString(nameHeaders[iHead]);
        HeaderList->Add(Header);
      }
    } else if (headerSelectionInt == 4){
      TObjString *Header1 = new TObjString("pi0EMC_3");
      HeaderList->Add(Header1);
    } else if (headerSelectionInt == 5){
      TObjString *Header1 = new TObjString("etaEMC_5");
      HeaderList->Add(Header1);
    } else if (headerSelectionInt == 6){
      TString nameHeaders[2]    = { "pi0EMC_3", "etaEMC_5" };
      for (Int_t iHead = 0; iHead < 2; iHead++ ){
        TObjString *Header = new TObjString(nameHeaders[iHead]);
        HeaderList->Add(Header);
      }
    } else if (headerSelectionInt == 7){
      TString nameHeaders[4]    = { "pi0_1", "eta_2", "pi0EMC_3", "etaEMC_5" };
      for (Int_t iHead = 0; iHead < 4; iHead++ ){
        TObjString *Header = new TObjString(nameHeaders[iHead]);
        HeaderList->Add(Header);
      }
    } else if (headerSelectionInt == 8){
      TObjString *Header1 = new TObjString("gEMCPhoton_7");
      HeaderList->Add(Header1);
    } else if (headerSelectionInt == 9){
      TString nameHeaders[10]   = { "Pythia_Jets_PtHard_1_10", "Pythia_Jets_PtHard_2_10", "Pythia_Jets_PtHard_3_10", "Pythia_Jets_PtHard_4_10", "Pythia_Jets_PtHard_5_10",
        "Pythia_Jets_PtHard_6_10", "Pythia_Jets_PtHard_7_10", "Pythia_Jets_PtHard_8_10", "Pythia_Jets_PtHard_9_10", "Pythia_Jets_PtHard_10_10"
      };
      for (Int_t iHead = 0; iHead < 10; iHead++ ){
        TObjString *Header = new TObjString(nameHeaders[iHead]);
        HeaderList->Add(Header);
      }
    } else if (headerSelectionInt == 10){
      TObjString *Header1 = new TObjString("pythia_bele_10_10");
      HeaderList->Add(Header1);
    } else if (headerSelectionInt == 11){
      TString nameHeaders[34]   = { "Pythia_Jets_PtHard_1_10", "Pythia_Jets_PtHard_2_10", "Pythia_Jets_PtHard_3_10", "Pythia_Jets_PtHard_4_10", "Pythia_Jets_PtHard_5_10",
        "Pythia_Jets_PtHard_6_10", "Pythia_Jets_PtHard_7_10", "Pythia_Jets_PtHard_8_10", "Pythia_Jets_PtHard_9_10", "Pythia_Jets_PtHard_10_10",
        "gEMCPhoton_7", "flat pt kstar_8", "flat pt kstarbar_9", "pythia_cele_10_10", "pythia_cele_18_10",
        "pythia_cele_30_10", "pythia_cele_50_10", "pythia_ccbar_10_10", "pythia_ccbar_18_10", "pythia_ccbar_30_10",
        "pythia_ccbar_50_10", "pythia_bele_10_10", "pythia_bele_18_10", "pythia_bele_30_10", "pythia_bele_50_10",
        "pythia_bbbar_10_10", "pythia_bbbar_18_10", "pythia_bbbar_30_10", "pythia_bbbar_50_10", "pi0_1",
        "eta_2", "pi0EMC_3", "etaEMC_5", "hijing_0"
      };
      for (Int_t iHead = 0; iHead < 34; iHead++ ){
        TObjString *Header = new TObjString(nameHeaders[iHead]);
        HeaderList->Add(Header);
      }
    } else if (headerSelectionInt == 12){
      TObjString *Header1 = new TObjString("pi0PHS_4");
      HeaderList->Add(Header1);
    } else if (headerSelectionInt == 13){
      TObjString *Header1 = new TObjString("etaPHS_6");
      HeaderList->Add(Header1);
    } else if (headerSelectionInt == 14){
      TString nameHeaders[2]    = { "pi0PHS_4", "etaPHS_6" };
      for (Int_t iHead = 0; iHead < 2; iHead++ ){
        TObjString *Header = new TObjString(nameHeaders[iHead]);
        HeaderList->Add(Header);
      }
    } else {
      TString nameHeaders[2]    = { "pi0_1", "eta_2" };
      for (Int_t iHead = 0; iHead < 2; iHead++ ){
        TObjString *Header = new TObjString(nameHeaders[iHead]);
        HeaderList->Add(Header);
      }
    }
  } else if (periodName.CompareTo("LHC14a1b")==0 || periodName.CompareTo("LHC14a1c")==0){
    if (headerSelectionInt == 1 || headerSelectionInt == 2 || headerSelectionInt == 3 ){
      TObjString *Header1 = new TObjString("BOX");
      HeaderList->Add(Header1);
    } if (headerSelectionInt == 4 || headerSelectionInt == 5 || headerSelectionInt == 6 ){
      TObjString *Header1 = new TObjString("PARAM_EMC");
      HeaderList->Add(Header1);
    } if (headerSelectionInt == 12 || headerSelectionInt == 13 || headerSelectionInt == 14 ){
      TObjString *Header1 = new TObjString("PARAM_PHOS");
      HeaderList->Add(Header1);
    }
  } else if (periodName.CompareTo("LHC16h4b")==0 || periodName.CompareTo("LHC16h4b2")==0){
    if (headerSelectionInt == 1){
      TObjString *Header1 = new TObjString("Injector (pi0)_1");
      HeaderList->Add(Header1);
    } else if (headerSelectionInt == 2){
      TObjString *Header1 = new TObjString("Injector (eta)_2");
      HeaderList->Add(Header1);
    } else {
      TObjString *Header1 = new TObjString("Injector (pi0)_1");
      HeaderList->Add(Header1);
      TObjString *Header2 = new TObjString("Injector (eta)_2");
      HeaderList->Add(Header2);
    }
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
    TString TrackMatcherName = Form("CaloTrackMatcher_%s_%i",caloCutPos.Data(),trackMatcherRunningMode);
    if( !(AliCaloTrackMatcher*)mgr->GetTask(TrackMatcherName.Data()) ){
      AliCaloTrackMatcher* fTrackMatcher = new AliCaloTrackMatcher(TrackMatcherName.Data(),caloCutPos.Atoi(),trackMatcherRunningMode);
      fTrackMatcher->SetV0ReaderName(V0ReaderName);
      fTrackMatcher->SetCorrectionTaskSetting(corrTaskSetting);
      mgr->AddTask(fTrackMatcher);
      mgr->ConnectInput(fTrackMatcher,0,cinput);
    }

    analysisEventCuts[i] = new AliConvEventCuts();
//     if ( trainConfig == 1){
//       if (periodName.CompareTo("LHC14a1a") ==0 || periodName.CompareTo("LHC14a1b") ==0 || periodName.CompareTo("LHC14a1c") ==0 ){
//         if ( i == 0 && doWeightingPart)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNamePtWeights, Form("Pi0_Hijing_%s_PbPb_2760GeV_0005TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_0005TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0005V0M","Eta_Fit_Data_PbPb_2760GeV_0005V0M");
//         if ( i == 1 && doWeightingPart)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNamePtWeights, Form("Pi0_Hijing_%s_PbPb_2760GeV_0510TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_0510TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0510V0M","Eta_Fit_Data_PbPb_2760GeV_0510V0M");
//         if ( i == 2 && doWeightingPart)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNamePtWeights, Form("Pi0_Hijing_%s_PbPb_2760GeV_0010TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_0010TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0010V0M","Eta_Fit_Data_PbPb_2760GeV_0010V0M");
//         if ( i == 3 && doWeightingPart)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNamePtWeights, Form("Pi0_Hijing_%s_PbPb_2760GeV_2040TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_2040TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_2040V0M","Eta_Fit_Data_PbPb_2760GeV_2040V0M");
//         if ( i == 4 && doWeightingPart)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNamePtWeights, Form("Pi0_Hijing_%s_PbPb_2760GeV_2050TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_2050TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_2050V0M","Eta_Fit_Data_PbPb_2760GeV_2050V0M");
//       }
//     }
    analysisEventCuts[i]->SetV0ReaderName(V0ReaderName);
    if (periodNameV0Reader.CompareTo("") != 0) analysisEventCuts[i]->SetPeriodEnum(periodNameV0Reader);

    if(periodName.Contains("LHC11h") && enableFlattening){
      cout << "entering the cent. flattening loop -> searching for file: " << fileNameCentFlattening.Data() << endl;

      if( fileNameCentFlattening.Contains("FlatFile") ){
        analysisEventCuts[i]->SetUseWeightFlatCentralityFromFile(enableFlattening, fileNameCentFlattening, "Cent");
      } else if( fileNameCentFlattening.Contains("Good") ){
        analysisEventCuts[i]->SetUseWeightFlatCentralityFromFile(enableFlattening, fileNameCentFlattening, "CentGoodRuns");
      }else if( fileNameCentFlattening.Contains("SemiGood") ){
        analysisEventCuts[i]->SetUseWeightFlatCentralityFromFile(enableFlattening, fileNameCentFlattening, "CentSemiGoodRuns");
      }else {
        analysisEventCuts[i]->SetUseWeightFlatCentralityFromFile(enableFlattening, fileNameCentFlattening, "CentTotalRuns");
      }
    }

    TString dataInputMultHisto  = "";
    TString mcInputMultHisto    = "";
    if (enableMultiplicityWeighting){
      cout << "INFO enableling mult weighting" << endl;
      if(periodNameAnchor.CompareTo("LHC15o")==0){
        TString cutNumber = cuts.GetEventCut(i);
        TString centCut = cutNumber(0,3);  // first three digits of event cut
        dataInputMultHisto = Form("%s_%s", periodNameAnchor.Data(), centCut.Data());
        mcInputMultHisto   = Form("%s_%s", periodName.Data(), centCut.Data());
        cout << "INFO read " << dataInputMultHisto.Data() << " and " <<  mcInputMultHisto.Data() << " from " << fileNameMultWeights.Data() << endl;
      }
      analysisEventCuts[i]->SetUseWeightMultiplicityFromFile(kTRUE, fileNameMultWeights, dataInputMultHisto, mcInputMultHisto );
    }

    if (trainConfig == 34 || trainConfig == 35){
      if (periodName.CompareTo("LHC14a1a") ==0 || periodName.CompareTo("LHC14a1b") ==0 || periodName.CompareTo("LHC14a1c") ==0 ){
        if ( doWeightingPart)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNamePtWeights, Form("Pi0_Hijing_%s_PbPb_2760GeV_0010TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_0010TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0010V0M","Eta_Fit_Data_PbPb_2760GeV_0010V0M");
      }
    }
    if (trainConfig == 36 || trainConfig == 37){
      if (periodName.CompareTo("LHC14a1a") ==0 || periodName.CompareTo("LHC14a1b") ==0 || periodName.CompareTo("LHC14a1c") ==0 ){
        if ( doWeightingPart)   analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNamePtWeights, Form("Pi0_Hijing_%s_PbPb_2760GeV_2050TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_2050TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_2050V0M","Eta_Fit_Data_PbPb_2760GeV_2050V0M");
      }
    }

    if (trainConfig == 38 || trainConfig == 39 || trainConfig == 43 || trainConfig == 44){
      if (periodName.CompareTo("LHC14a1a") ==0 || periodName.CompareTo("LHC14a1b") ==0 || periodName.CompareTo("LHC14a1c") ==0 ){
        if ( doWeightingPart)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNamePtWeights, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_0010TPC",periodName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_0010TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0010V0M","Eta_Fit_Data_PbPb_2760GeV_0010V0M");
      }
    }

    if (trainConfig == 40 || trainConfig == 41){
      if (periodName.CompareTo("LHC14a1a") ==0 || periodName.CompareTo("LHC14a1b") ==0 || periodName.CompareTo("LHC14a1c") ==0 ){
        if ( doWeightingPart)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNamePtWeights, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_2050TPC",periodName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_2050TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_2050V0M","Eta_Fit_Data_PbPb_2760GeV_2050V0M");
      }
    }

    if (enableLightOutput > 0) analysisEventCuts[i]->SetLightOutput(kTRUE);
    analysisEventCuts[i]->InitializeCutsFromCutString((cuts.GetEventCut(i)).Data());
    if (periodName.CompareTo("LHC14a1b") ==0 || periodName.CompareTo("LHC14a1c") ==0 ){
      if (headerSelectionInt == 1 || headerSelectionInt == 4 || headerSelectionInt == 12 ) analysisEventCuts[i]->SetAddedSignalPDGCode(111);
      if (headerSelectionInt == 2 || headerSelectionInt == 5 || headerSelectionInt == 13 ) analysisEventCuts[i]->SetAddedSignalPDGCode(221);
    }
    analysisEventCuts[i]->SetCorrectionTaskSetting(corrTaskSetting);
    EventCutList->Add(analysisEventCuts[i]);
    analysisEventCuts[i]->SetFillCutHistograms("",kFALSE);

    analysisCuts[i] = new AliConversionPhotonCuts();
    analysisCuts[i]->SetV0ReaderName(V0ReaderName);
    if (enableLightOutput > 0) analysisCuts[i]->SetLightOutput(kTRUE);
    analysisCuts[i]->InitializeCutsFromCutString((cuts.GetPhotonCut(i)).Data());
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

  task->SetAllowOverlapHeaders(enableHeaderOverlap);
  task->SetEventCutList(numberOfCuts,EventCutList);
  task->SetConversionCutList(numberOfCuts,ConvCutList);
  task->SetCaloCutList(numberOfCuts,ClusterCutList);
  task->SetMesonCutList(numberOfCuts,MesonCutList);
  task->SetCorrectionTaskSetting(corrTaskSetting);
  task->SetMoveParticleAccordingToVertex(kTRUE);
  task->SetDoMesonAnalysis(kTRUE);
  task->SetDoMesonQA(enableQAMesonTask); //Attention new switch for Pi0 QA
  task->SetDoPhotonQA(enableQAPhotonTask);  //Attention new switch small for Photon QA
  task->SetDoClusterQA(1);  //Attention new switch small for Cluster QA
  task->SetEnableSortingOfMCClusLabels(enableSortingMCLabels);
  task->SetDoTreeInvMassShowerShape(doTreeClusterShowerShape);
  task->SetUseTHnSparse(enableTHnSparse);
  if(enableExtMatchAndQA > 1){ task->SetPlotHistsExtQA(kTRUE);}

  //connect containers
  AliAnalysisDataContainer *coutput =
    mgr->CreateContainer(!(corrTaskSetting.CompareTo("")) ? Form("GammaConvCalo_%i",trainConfig) : Form("GammaConvCalo_%i_%s",trainConfig,corrTaskSetting.Data()), TList::Class(),
              AliAnalysisManager::kOutputContainer,Form("GammaConvCalo_%i.root",trainConfig));

  mgr->AddTask(task);
  mgr->ConnectInput(task,0,cinput);
  mgr->ConnectOutput(task,1,coutput);
  if(enableQAPhotonTask>1){
    for(Int_t i = 0; i<numberOfCuts; i++){
      mgr->ConnectOutput(task,2+i,mgr->CreateContainer(Form("%s_%s_%s_%s Photon DCA tree",(cuts.GetEventCut(i)).Data(),(cuts.GetPhotonCut(i)).Data(),cuts.GetClusterCut(i).Data(),(cuts.GetMesonCut(i)).Data()), TTree::Class(), AliAnalysisManager::kOutputContainer, Form("GammaConvCalo_%i.root",trainConfig)) );
    }
  }

  return;

}
