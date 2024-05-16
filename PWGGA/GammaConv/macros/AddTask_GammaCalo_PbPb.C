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
//($ALIPHYSICS/PWGGA/GammaConv/AliAnalysisTaskGammaCalo.cxx) for
//PbPb together with all supporting classes
//***************************************************************************************

//***************************************************************************************
//main function
//***************************************************************************************
void AddTask_GammaCalo_PbPb(
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
  Bool_t    enableHeaderOverlap           = kTRUE,    // enable trigger overlap rejection
  TString   settingMaxFacPtHard           = "3.",     // maximum factor between hardest jet and ptHard generated
  Int_t     debugLevel                    = 0,        // introducing debug levels for grid running
  // settings for weights
  // FPTW:fileNamePtWeights, FMUW:fileNameMultWeights, FCEF:fileNameCentFlattening, separate with ;
  TString   fileNameExternalInputs        = "",
  Bool_t    doWeighting                   = kFALSE,               // enable Weighting
  Int_t     headerSelectionInt            = 0,        // enable Weighting
  TString   generatorName                 = "DPMJET", // generator Name
  Bool_t    enableMultiplicityWeighting   = kFALSE,   //
  TString   periodNameAnchor              = "",       //
  // special settings
  Bool_t    enableSortingMCLabels         = kTRUE,    // enable sorting for MC cluster labels
  Bool_t    doFlattening                  = kFALSE,   // switch on centrality flattening for LHC11h
  Bool_t    doPrimaryTrackMatching        = kTRUE,    // enable basic track matching for all primary tracks to cluster    
  // subwagon config
  TString   additionalTrainConfig         = "0"       // additional counter for trainconfig
  ) {

  AliCutHandlerPCM cuts;


  TString fileNamePtWeights     = cuts.GetSpecialFileNameFromString (fileNameExternalInputs, "FPTW:");
  TString fileNameMultWeights   = cuts.GetSpecialFileNameFromString (fileNameExternalInputs, "FMUW:");
  TString fileNameCentFlattening= cuts.GetSpecialFileNameFromString (fileNameExternalInputs, "FCEF:");

  TString addTaskName                 = "AddTask_GammaCalo_PbPb";
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
  if(rmaxFacPtHardSetting->GetEntries()<1){cout << "ERROR: AddTask_GammaCalo_PbPb during parsing of settingMaxFacPtHard String '" << settingMaxFacPtHard.Data() << "'" << endl; return;}
  Bool_t fMinPtHardSet        = kFALSE;
  Double_t minFacPtHard       = -1;
  Bool_t fMaxPtHardSet        = kFALSE;
  Double_t maxFacPtHard       = 100;
  Bool_t fSingleMaxPtHardSet  = kFALSE;
  Double_t maxFacPtHardSingle = 100;
  Bool_t fJetFinderUsage      = kFALSE;
  Bool_t fUsePtHardFromFile      = kFALSE;
  Bool_t fUseAddOutlierRej      = kFALSE;
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
    } else if(strTempSetting.BeginsWith("USEJETFINDER:")){
      strTempSetting.Replace(0,13,"");
      if(strTempSetting.Atoi()==1){
        cout << "using MC jet finder for outlier removal" << endl;
        fJetFinderUsage        = kTRUE;
      }
    } else if(strTempSetting.BeginsWith("PTHFROMFILE:")){
      strTempSetting.Replace(0,12,"");
      if(strTempSetting.Atoi()==1){
        cout << "using MC jet finder for outlier removal" << endl;
        fUsePtHardFromFile        = kTRUE;
      }
    } else if(strTempSetting.BeginsWith("ADDOUTLIERREJ:")){
      strTempSetting.Replace(0,14,"");
      if(strTempSetting.Atoi()==1){
        cout << "using path based outlier removal" << endl;
        fUseAddOutlierRej        = kTRUE;
      }
    } else if(rmaxFacPtHardSetting->GetEntries()==1 && strTempSetting.Atof()>0){
      maxFacPtHard               = strTempSetting.Atof();
      cout << "running with max pT hard jet fraction of: " << maxFacPtHard << endl;
      fMaxPtHardSet        = kTRUE;
    }
  }


  Int_t isHeavyIon = 1;

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
  AliAnalysisTaskGammaCalo *task=NULL;
  task= new AliAnalysisTaskGammaCalo(Form("GammaCalo_%i",trainConfig));
  task->SetIsHeavyIon(isHeavyIon);
  task->SetIsMC(isMC);
  task->SetV0ReaderName(V0ReaderName);
  if (enableLightOutput >= 4) task->SetECalibOutput(kTRUE);
  task->SetCorrectionTaskSetting(corrTaskSetting);
  task->SetLightOutput(enableLightOutput);
  task->SetDoPrimaryTrackMatching(doPrimaryTrackMatching);
  task->SetTrackMatcherRunningMode(trackMatcherRunningMode);

  // meson cuts
  // meson type (Dalitz or not), BG scheme, pool depth, rotation degrees, rapidity cut, radius cut, alpha, chi2, shared electrons, reject to close v0, MC smearing, dca, dca, dca

  if (trainConfig == 1){ // EMCAL clusters
    cuts.AddCutCalo("60100013","1111100053032230000","0163103100000050"); // 0-5%
    cuts.AddCutCalo("61200013","1111100053032230000","0163103100000050"); // 5-10%
    cuts.AddCutCalo("50100013","1111100053032230000","0163103100000050"); // 0-10%
    cuts.AddCutCalo("52400013","1111100053032230000","0163103100000050"); // 20-40%
    cuts.AddCutCalo("52500013","1111100053032230000","0163103100000050"); // 20-50%
  } else if (trainConfig == 2){ // EMCAL clusters no timing cut
    cuts.AddCutCalo("60100013","1111100003032230000","0163103100000050"); // 0-5%
    cuts.AddCutCalo("61200013","1111100003032230000","0163103100000050"); // 5-10%
    cuts.AddCutCalo("50100013","1111100003032230000","0163103100000050"); // 0-10%
    cuts.AddCutCalo("52400013","1111100003032230000","0163103100000050"); // 20-40%
    cuts.AddCutCalo("52500013","1111100003032230000","0163103100000050"); // 20-50%
  } else if (trainConfig == 3){ // EMCAL clusters
    cuts.AddCutCalo("60100013","1111100053032230000","0163103100000050"); // 0-5%
    cuts.AddCutCalo("61200013","1111100053032230000","0163103100000050"); // 5-10%
    cuts.AddCutCalo("50100013","1111100053032230000","0163103100000050"); // 0-10%
    cuts.AddCutCalo("51200013","1111100053032230000","0163103100000050"); // 10-20%
    cuts.AddCutCalo("52400013","1111100053032230000","0163103100000050"); // 20-40%
  } else if (trainConfig == 4){ // EMCAL clusters no timing cut
    cuts.AddCutCalo("60100013","1111100003032230000","0163103100000050"); // 0-5%
    cuts.AddCutCalo("61200013","1111100003032230000","0163103100000050"); // 5-10%
    cuts.AddCutCalo("50100013","1111100003032230000","0163103100000050"); // 0-10%
    cuts.AddCutCalo("51200013","1111100003032230000","0163103100000050"); // 10-20%
    cuts.AddCutCalo("52400013","1111100003032230000","0163103100000050"); // 20-40%
  } else if (trainConfig == 5){ // EMCAL clusters
    cuts.AddCutCalo("54600013","1111100053032230000","0163103100000050"); // 40-60%
    cuts.AddCutCalo("56800013","1111100053032230000","0163103100000050"); // 60-80%
    cuts.AddCutCalo("52600013","1111100053032230000","0163103100000050"); // 20-60%
    cuts.AddCutCalo("54800013","1111100053032230000","0163103100000050"); // 40-80%
    cuts.AddCutCalo("52500013","1111100053032230000","0163103100000050"); // 20-50%
  } else if (trainConfig == 6){ // EMCAL clusters  no timing cut
    cuts.AddCutCalo("54600013","1111100003032230000","0163103100000050"); // 40-60%
    cuts.AddCutCalo("56800013","1111100003032230000","0163103100000050"); // 60-80%
    cuts.AddCutCalo("52600013","1111100003032230000","0163103100000050"); // 20-60%
    cuts.AddCutCalo("54800013","1111100003032230000","0163103100000050"); // 40-80%
    cuts.AddCutCalo("52500013","1111100003032230000","0163103100000050"); // 20-50%
  } else if (trainConfig == 7){ // EMCAL clusters
    cuts.AddCutCalo("60100013","1111100053032230000","0163103100000050"); // 0-5%
    cuts.AddCutCalo("61200013","1111100053032230000","0163103100000050"); // 5-10%
    cuts.AddCutCalo("50100013","1111100053032230000","0163103100000050"); // 0-10%
  } else if (trainConfig == 8){ // EMCAL clusters no timing cut
    cuts.AddCutCalo("60100013","1111100003032230000","0163103100000050"); // 0-5%
    cuts.AddCutCalo("61200013","1111100003032230000","0163103100000050"); // 5-10%
    cuts.AddCutCalo("50100013","1111100003032230000","0163103100000050"); // 0-10%
  } else if (trainConfig == 9){ // EMCAL clusters
    cuts.AddCutCalo("51200013","1111100053032230000","0163103100000050"); // 10-20%
    cuts.AddCutCalo("52400013","1111100053032230000","0163103100000050"); // 20-40%
    cuts.AddCutCalo("52500013","1111100053032230000","0163103100000050"); // 20-50%
  } else if (trainConfig == 10){ // EMCAL clusters no timing cut
    cuts.AddCutCalo("51200013","1111100003032230000","0163103100000050"); // 10-20%
    cuts.AddCutCalo("52400013","1111100003032230000","0163103100000050"); // 20-40%
    cuts.AddCutCalo("52500013","1111100003032230000","0163103100000050"); // 20-50%
  } else if (trainConfig == 11){ // EMCAL clusters
    cuts.AddCutCalo("54600013","1111100053032230000","0163103100000050"); // 40-60%
    cuts.AddCutCalo("56800013","1111100053032230000","0163103100000050"); // 60-80%
  } else if (trainConfig == 12){ // EMCAL clusters no timing cut
    cuts.AddCutCalo("54600013","1111100003032230000","0163103100000050"); // 40-60%
    cuts.AddCutCalo("56800013","1111100003032230000","0163103100000050"); // 60-80%

  } else if (trainConfig == 13){ // EMCAL clusters
    cuts.AddCutCalo("50900013","1111100053032230000","0163103100000050"); // 0-90%
  } else if (trainConfig == 14){ // EMCAL clusters no timing cut
    cuts.AddCutCalo("50900013","1111100003032230000","0163103100000050"); // 0-90%
  } else if (trainConfig == 15){ // EMCAL clusters
    cuts.AddCutCalo("50100013","1111100053032230000","0163103100000050"); // 0-90%
    cuts.AddCutCalo("51300013","1111100053032230000","0163103100000050"); // 0-90%
    cuts.AddCutCalo("53500013","1111100053032230000","0163103100000050"); // 0-90%
    cuts.AddCutCalo("55900013","1111100053032230000","0163103100000050"); // 0-90%
  } else if (trainConfig == 16){ //EG2
    cuts.AddCutCalo("50185013","1111100053032230000","0163103100000050"); //
    cuts.AddCutCalo("51385013","1111100053032230000","0163103100000050"); //
    cuts.AddCutCalo("53585013","1111100053032230000","0163103100000050"); //
    cuts.AddCutCalo("55985013","1111100053032230000","0163103100000050"); //
  } else if (trainConfig == 17){ //EG1
    cuts.AddCutCalo("50183013","1111100053032230000","0163103100000050"); //
    cuts.AddCutCalo("51383013","1111100053032230000","0163103100000050"); //
    cuts.AddCutCalo("53583013","1111100053032230000","0163103100000050"); //
    cuts.AddCutCalo("55983013","1111100053032230000","0163103100000050"); //
  } else if (trainConfig == 18){ //EJ2
    cuts.AddCutCalo("50195013","1111100053032230000","0163103100000050"); //
    cuts.AddCutCalo("51395013","1111100053032230000","0163103100000050"); //
    cuts.AddCutCalo("53595013","1111100053032230000","0163103100000050"); //
    cuts.AddCutCalo("55995013","1111100053032230000","0163103100000050"); //
  } else if (trainConfig == 19){ //EJ1
    cuts.AddCutCalo("50193013","1111100053032230000","0163103100000050"); //
    cuts.AddCutCalo("51393013","1111100053032230000","0163103100000050"); //
    cuts.AddCutCalo("53593013","1111100053032230000","0163103100000050"); //
    cuts.AddCutCalo("55993013","1111100053032230000","0163103100000050"); //

    // EMC triggers LHC15o
  } else if (trainConfig == 20){ //MB
    cuts.AddCutCalo("10110013","1111100051032230000","0163103100000050"); //
    cuts.AddCutCalo("11310013","1111100051032230000","0163103100000050"); //
    cuts.AddCutCalo("13510013","1111100051032230000","0163103100000050"); //
    cuts.AddCutCalo("15010013","1111100051032230000","0163103100000050"); //
  } else if (trainConfig == 21){ //EG2
    cuts.AddCutCalo("10185013","1111100051032230000","0163103100000050"); //
    cuts.AddCutCalo("11385013","1111100051032230000","0163103100000050"); //
    cuts.AddCutCalo("13585013","1111100051032230000","0163103100000050"); //
    cuts.AddCutCalo("15085013","1111100051032230000","0163103100000050"); //
  } else if (trainConfig == 22){ //EG1
    cuts.AddCutCalo("10183013","1111100051032230000","0163103100000050"); //
    cuts.AddCutCalo("11383013","1111100051032230000","0163103100000050"); //
    cuts.AddCutCalo("13583013","1111100051032230000","0163103100000050"); //
    cuts.AddCutCalo("15083013","1111100051032230000","0163103100000050"); //
  } else if (trainConfig == 23){ //EJ2
    cuts.AddCutCalo("10195013","1111100051032230000","0163103100000050"); //
    cuts.AddCutCalo("11395013","1111100051032230000","0163103100000050"); //
    cuts.AddCutCalo("13595013","1111100051032230000","0163103100000050"); //
    cuts.AddCutCalo("15095013","1111100051032230000","0163103100000050"); //
  } else if (trainConfig == 24){ //EJ1
    cuts.AddCutCalo("10193013","1111100051032230000","0163103100000050"); //
    cuts.AddCutCalo("11393013","1111100051032230000","0163103100000050"); //
    cuts.AddCutCalo("13593013","1111100051032230000","0163103100000050"); //
    cuts.AddCutCalo("15093013","1111100051032230000","0163103100000050"); //
    // DCAL triggers LHC15o
  } else if (trainConfig == 25){ //MB
    cuts.AddCutCalo("10110013","3885500017032220000","0163103100000050"); //
    cuts.AddCutCalo("11310013","3885500017032220000","0163103100000050"); //
    cuts.AddCutCalo("13510013","3885500017032220000","0163103100000050"); //
    cuts.AddCutCalo("15010013","3885500017032220000","0163103100000050"); //
  } else if (trainConfig == 26){ //EG2
    cuts.AddCutCalo("10185013","3885500017032220000","0163103100000050"); //
    cuts.AddCutCalo("11385013","3885500017032220000","0163103100000050"); //
    cuts.AddCutCalo("13585013","3885500017032220000","0163103100000050"); //
    cuts.AddCutCalo("15085013","3885500017032220000","0163103100000050"); //
  } else if (trainConfig == 27){ //EG1
    cuts.AddCutCalo("10183013","3885500017032220000","0163103100000050"); //
    cuts.AddCutCalo("11383013","3885500017032220000","0163103100000050"); //
    cuts.AddCutCalo("13583013","3885500017032220000","0163103100000050"); //
    cuts.AddCutCalo("15083013","3885500017032220000","0163103100000050"); //
  } else if (trainConfig == 28){ //EJ2
    cuts.AddCutCalo("10195013","3885500017032220000","0163103100000050"); //
    cuts.AddCutCalo("11395013","3885500017032220000","0163103100000050"); //
    cuts.AddCutCalo("13595013","3885500017032220000","0163103100000050"); //
    cuts.AddCutCalo("15095013","3885500017032220000","0163103100000050"); //
  } else if (trainConfig == 29){ //EJ1
    cuts.AddCutCalo("10193013","3885500017032220000","0163103100000050"); //
    cuts.AddCutCalo("11393013","3885500017032220000","0163103100000050"); //
    cuts.AddCutCalo("13593013","3885500017032220000","0163103100000050"); //
    cuts.AddCutCalo("15093013","3885500017032220000","0163103100000050"); //

  // EMCal trigger for LHC11h
  } else if (trainConfig == 30){ // EMCAL clusters
    cuts.AddCutCalo("50900013","11111010500a2220000","01331031000000d0"); // 0-90%
    cuts.AddCutCalo("50100013","11111010500a2220000","01331061000000d0"); //
    cuts.AddCutCalo("51300013","11111010500b2220000","01331061000000d0"); //
    cuts.AddCutCalo("53500013","1111101050032220000","01331031000000d0"); //
    cuts.AddCutCalo("55900013","1111101050032220000","01331031000000d0"); //
  } else if (trainConfig == 31){ // EMCAL clusters no timing cut
    cuts.AddCutCalo("50980013","11111010500a2220000","01331031000000d0"); // 0-90%
    cuts.AddCutCalo("50180013","11111010500a2220000","01331061000000d0"); //
    cuts.AddCutCalo("51380013","11111010500b2220000","01331061000000d0"); //
    cuts.AddCutCalo("53580013","1111101050032220000","01331031000000d0"); //
    cuts.AddCutCalo("55980013","1111101050032220000","01331031000000d0"); //
   } else if (trainConfig == 32){ // EMCAL clusters
    cuts.AddCutCalo("50980013","1111100053032230000","0163103100000050"); // 0-90%
  } else if (trainConfig == 33){ // EMCAL clusters no timing cut
    cuts.AddCutCalo("50980013","1111100003032230000","0163103100000050"); // 0-90%
  } else if (trainConfig == 34){ // EMCAL clusters
    cuts.AddCutCalo("50100013","1111102053032230000","01631031000000d0"); // 0-10
    cuts.AddCutCalo("50100013","1111171053032230000","01631031000000d0"); // 0-10
  } else if (trainConfig == 35){ // EMCAL clusters
    cuts.AddCutCalo("50100013","11111020530a2230000","01631031000000d0"); // 0-10 // reproduce Astrids cuts with opening angle cut
    cuts.AddCutCalo("50100013","11111020530a2230000","0163103100000000"); // 0-10 // reproduce Astrids cuts without opening angle cut
    cuts.AddCutCalo("50100013","11111710530a2230000","01631031000000d0"); // 0-10
  } else if (trainConfig == 36){ // EMCAL clusters
    cuts.AddCutCalo("52500013","1111102053032230000","01631031000000d0"); // 20-50
    cuts.AddCutCalo("52500013","1111172053032230000","01631031000000d0"); // 20-50
  } else if (trainConfig == 37){ // EMCAL clusters
    cuts.AddCutCalo("52500013","11111020530a2230000","01631031000000d0"); // 20-50 // reproduce Astrids cuts with opening angle cut
    cuts.AddCutCalo("52500013","11111020530a2230000","0163103100000000"); // 20-50 // reproduce Astrids cuts without opening angle cut
    cuts.AddCutCalo("52500013","11111720530a2230000","01631031000000d0"); // 20-50
  } else if (trainConfig == 38){ // EMCAL clusters - added signals
    cuts.AddCutCalo("50100023","1111171053032230000","01631031000000d0"); // 0-10
    cuts.AddCutCalo("50100023","11111710530a2230000","01631031000000d0"); // 0-10 // reproduce Astrids cuts with opening angle cut
    cuts.AddCutCalo("50100023","11111710530b2230000","01631031000000d0"); // 0-10
  } else if (trainConfig == 39){ // EMCAL clusters - added signals
    cuts.AddCutCalo("50100023","1111171053032230000","01631031000000d0"); // 0-10
    cuts.AddCutCalo("50100023","11111710530a2230000","01631031000000d0"); // 0-10
    cuts.AddCutCalo("50100023","11111710530b2230000","01631031000000d0"); // 0-10
  } else if (trainConfig == 40){ // EMCAL clusters - added signals
    cuts.AddCutCalo("52500023","1111172053032230000","01631031000000d0"); // 20-50
    cuts.AddCutCalo("52500023","11111720530a2230000","01631031000000d0"); // 20-50
    cuts.AddCutCalo("52500023","11111720530b2230000","01631031000000d0"); // 20-50
  } else if (trainConfig == 41){ // EMCAL clusters - added signals
    cuts.AddCutCalo("52500023","1111172053032230000","01631031000000d0"); // 20-50
    cuts.AddCutCalo("52500023","11111720530a2230000","01631031000000d0"); // 20-50
    cuts.AddCutCalo("52500023","11111720530b2230000","01631031000000d0"); // 20-50
  } else if (trainConfig == 42){ // EMCAL clusters - all headers
    cuts.AddCutCalo("50100003","11111020530a2230000","01631031000000d0"); // 0-10   // reproduce Astrids cuts with opening angle cut
    cuts.AddCutCalo("50100003","11111020530a2230000","0163103100000000"); // 0-10   // reproduce Astrids cuts with opening angle cut
    cuts.AddCutCalo("50100003","11111710530a2230000","01631031000000d0"); // 0-10   // new calib
  } else if (trainConfig == 43){ // EMCAL clusters - added signals pi0 forseen
    cuts.AddCutCalo("50100023","11111020530a2230000","01631031000000d0"); // 0-10  // reproduce Astrids cuts with opening angle cut
    cuts.AddCutCalo("50100023","11111020530a2230000","0163103100000000"); // 0-10  // reproduce Astrids cuts without opening angle cut
    cuts.AddCutCalo("50100023","11111710530a2230000","01631031000000d0");
  } else if (trainConfig == 44){ // EMCAL clusters - added signals eta forseen
    cuts.AddCutCalo("50100023","11111020530a2230000","01631031000000d0"); // 0-10  // reproduce Astrids cuts with opening angle cut
    cuts.AddCutCalo("50100023","11111020530a2230000","0163103100000000"); // 0-10  // reproduce Astrids cuts without opening angle cut
    cuts.AddCutCalo("50100023","11111710530a2230000","01631031000000d0");
  } else if (trainConfig == 45){ // EMCAL clusters - added signals other
    cuts.AddCutCalo("50100023","11111020530a2230000","01631031000000d0"); // 0-10  // reproduce Astrids cuts with opening angle cut
    cuts.AddCutCalo("50100023","11111020530a2230000","0163103100000000"); // 0-10  // reproduce Astrids cuts without opening angle cut
    cuts.AddCutCalo("50100023","11111710530a2230000","01631031000000d0");
  } else if (trainConfig == 46){ // EMCAL clusters - added signals other
    cuts.AddCutCalo("50100023","11111020530a2230000","01631031000000d0"); // 0-10  // reproduce Astrids cuts with opening angle cut
    cuts.AddCutCalo("50100023","11111020530a2230000","0163103100000000"); // 0-10  // reproduce Astrids cuts without opening angle cut
    cuts.AddCutCalo("50100023","11111710530a2230000","01631031000000d0");
  } else if (trainConfig == 47){ // EMCAL clusters V0 Cent
    cuts.AddCutCalo("10100013","11111020530a2230000","01631031000000d0"); // 0-10 // reproduce Astrids cuts with opening angle cut
    cuts.AddCutCalo("10100013","11111020530a2230000","0163103100000000"); // 0-10 // reproduce Astrids cuts without opening angle cut
    cuts.AddCutCalo("10100013","11111710530a2230000","01631031000000d0"); // 0-10

  // trainconfig for PbPb studies in 2.76 TeV with TB  nonlin
  } else if (trainConfig == 48){
    cuts.AddCutCalo("10100013","11111020530a2230000","01631031000000d0"); // 0-10 V0M  // reproduce Astrids cuts with opening angle cut
    cuts.AddCutCalo("11200013","11111020530a2230000","01631031000000d0"); // 10-20 V0M // reproduce Astrids cuts with opening angle cut
  } else if (trainConfig == 49){
    cuts.AddCutCalo("30100013","11111020530a2230000","01631031000000d0"); // 0-5 V0M   // reproduce Astrids cuts with opening angle cut
    cuts.AddCutCalo("31200013","11111020530a2230000","01631031000000d0"); // 5-10 V0M  // reproduce Astrids cuts with opening angle cut
  } else if (trainConfig == 50){
    cuts.AddCutCalo("12300013","11111020530a2230000","01631031000000d0"); // 20-30 V0M // reproduce Astrids cuts with opening angle cut
    cuts.AddCutCalo("13400013","11111020530a2230000","01631031000000d0"); // 30-40 V0M // reproduce Astrids cuts with opening angle cut
    cuts.AddCutCalo("14500013","11111020530a2230000","01631031000000d0"); // 40-50 V0M // reproduce Astrids cuts with opening angle cut
    cuts.AddCutCalo("15600013","11111020530a2230000","01631031000000d0"); // 50-60 V0M // reproduce Astrids cuts with opening angle cut
  } else if (trainConfig == 51){
    cuts.AddCutCalo("16700013","11111020530a2230000","01631031000000d0"); // 60-70 V0M // reproduce Astrids cuts with opening angle cut
    cuts.AddCutCalo("17800013","11111020530a2230000","01631031000000d0"); // 70-80 V0M // reproduce Astrids cuts with opening angle cut
    cuts.AddCutCalo("18900013","11111020530a2230000","01631031000000d0"); // 80-90 V0M // reproduce Astrids cuts with opening angle cut
    cuts.AddCutCalo("14600013","11111020530a2230000","01631031000000d0"); // 40-60 V0M // reproduce Astrids cuts with opening angle cut
    cuts.AddCutCalo("16800013","11111020530a2230000","01631031000000d0"); // 60-80 V0M // reproduce Astrids cuts with opening angle cut
    // trainconfig for PbPb studies in 2.76 TeV with no  nonlin
  } else if (trainConfig == 52){
    cuts.AddCutCalo("10100013","11111000530a2230000","01631031000000d0"); // 0-10 V0M  // reproduce Astrids cuts with opening angle cut
    cuts.AddCutCalo("11200013","11111000530a2230000","01631031000000d0"); // 10-20 V0M // reproduce Astrids cuts with opening angle cut
  } else if (trainConfig == 53){
    cuts.AddCutCalo("30100013","11111000530a2230000","01631031000000d0"); // 0-5 V0M   // reproduce Astrids cuts with opening angle cut
    cuts.AddCutCalo("31200013","11111000530a2230000","01631031000000d0"); // 5-10 V0M  // reproduce Astrids cuts with opening angle cut
  } else if (trainConfig == 54){
    cuts.AddCutCalo("12300013","11111000530a2230000","01631031000000d0"); // 20-30 V0M // reproduce Astrids cuts with opening angle cut
    cuts.AddCutCalo("13400013","11111000530a2230000","01631031000000d0"); // 30-40 V0M // reproduce Astrids cuts with opening angle cut
    cuts.AddCutCalo("14500013","11111000530a2230000","01631031000000d0"); // 40-50 V0M // reproduce Astrids cuts with opening angle cut
    cuts.AddCutCalo("15600013","11111000530a2230000","01631031000000d0"); // 50-60 V0M // reproduce Astrids cuts with opening angle cut
  } else if (trainConfig == 55){
    cuts.AddCutCalo("16700013","11111000530a2230000","01631031000000d0"); // 60-70 V0M // reproduce Astrids cuts with opening angle cut
    cuts.AddCutCalo("17800013","11111000530a2230000","01631031000000d0"); // 70-80 V0M // reproduce Astrids cuts with opening angle cut
    cuts.AddCutCalo("18900013","11111000530a2230000","01631031000000d0"); // 80-90 V0M // reproduce Astrids cuts with opening angle cut
    cuts.AddCutCalo("14600013","11111000530a2230000","01631031000000d0"); // 40-60 V0M // reproduce Astrids cuts with opening angle cut
    cuts.AddCutCalo("16800013","11111000530a2230000","01631031000000d0"); // 60-80 V0M // reproduce Astrids cuts with opening angle cut
    // trainconfig for PbPb studies in 2.76 TeV with TB nonlin
  } else if (trainConfig == 56){
    cuts.AddCutCalo("50100013","11111020530a2230000","01631031000000d0"); // 0-10 TM  // reproduce Astrids cuts with opening angle cut
    cuts.AddCutCalo("51200013","11111020530a2230000","01631031000000d0"); // 10-20 TM // reproduce Astrids cuts with opening angle cut
  } else if (trainConfig == 57){
    cuts.AddCutCalo("60100013","11111020530a2230000","01631031000000d0"); // 0-5 TM   // reproduce Astrids cuts with opening angle cut
    cuts.AddCutCalo("61100013","11111020530a2230000","01631031000000d0"); // 5-10 TM  // reproduce Astrids cuts with opening angle cut
  } else if (trainConfig == 58){
    cuts.AddCutCalo("52300013","11111020530a2230000","01631031000000d0"); // 20-30 TM // reproduce Astrids cuts with opening angle cut
    cuts.AddCutCalo("53400013","11111020530a2230000","01631031000000d0"); // 30-40 TM // reproduce Astrids cuts with opening angle cut
    cuts.AddCutCalo("54500013","11111020530a2230000","01631031000000d0"); // 40-50 TM // reproduce Astrids cuts with opening angle cut
    cuts.AddCutCalo("55600013","11111020530a2230000","01631031000000d0"); // 50-60 TM // reproduce Astrids cuts with opening angle cut
  } else if (trainConfig == 59){
    cuts.AddCutCalo("56700013","11111020530a2230000","01631031000000d0"); // 60-70 TM // reproduce Astrids cuts with opening angle cut
    cuts.AddCutCalo("57800013","11111020530a2230000","01631031000000d0"); // 70-80 TM // reproduce Astrids cuts with opening angle cut
    cuts.AddCutCalo("58900013","11111020530a2230000","01631031000000d0"); // 80-90 TM // reproduce Astrids cuts with opening angle cut
    cuts.AddCutCalo("54600013","11111020530a2230000","01631031000000d0"); // 40-60 TM // reproduce Astrids cuts with opening angle cut
    cuts.AddCutCalo("56800013","11111020530a2230000","01631031000000d0"); // 60-80 TM // reproduce Astrids cuts with opening angle cut

  // trainconfig for PbPb studies in 2.76 TeV with TB  nonlin
  } else if (trainConfig == 60){
    cuts.AddCutCalo("10100023","11111020530a2230000","01631031000000d0"); // 0-10 V0M  // reproduce Astrids cuts with opening angle cut
    cuts.AddCutCalo("11200023","11111020530a2230000","01631031000000d0"); // 10-20 V0M // reproduce Astrids cuts with opening angle cut
  } else if (trainConfig == 61){
    cuts.AddCutCalo("30100023","11111020530a2230000","01631031000000d0"); // 0-5 V0M   // reproduce Astrids cuts with opening angle cut
    cuts.AddCutCalo("31200023","11111020530a2230000","01631031000000d0"); // 5-10 V0M  // reproduce Astrids cuts with opening angle cut
  } else if (trainConfig == 62){
    cuts.AddCutCalo("10100023","11111020530a2230000","01631031000000d0"); // 0-10 V0M  // reproduce Astrids cuts with opening angle cut
    cuts.AddCutCalo("11200023","11111020530a2230000","01631031000000d0"); // 10-20 V0M // reproduce Astrids cuts with opening angle cut
    cuts.AddCutCalo("30100023","11111020530a2230000","01631031000000d0"); // 0-5 V0M   // reproduce Astrids cuts with opening angle cut
    cuts.AddCutCalo("31200023","11111020530a2230000","01631031000000d0"); // 5-10 V0M  // reproduce Astrids cuts with opening angle cut
  } else if (trainConfig == 63){
    cuts.AddCutCalo("12300023","11111020530a2230000","01631031000000d0"); // 20-30 V0M // reproduce Astrids cuts with opening angle cut
    cuts.AddCutCalo("13400023","11111020530a2230000","01631031000000d0"); // 30-40 V0M // reproduce Astrids cuts with opening angle cut
    cuts.AddCutCalo("14500023","11111020530a2230000","01631031000000d0"); // 40-50 V0M // reproduce Astrids cuts with opening angle cut
    cuts.AddCutCalo("15600023","11111020530a2230000","01631031000000d0"); // 50-60 V0M // reproduce Astrids cuts with opening angle cut
  } else if (trainConfig == 64){
    cuts.AddCutCalo("12300023","11111020530a2230000","01631031000000d0"); // 20-30 V0M // reproduce Astrids cuts with opening angle cut
    cuts.AddCutCalo("13400023","11111020530a2230000","01631031000000d0"); // 30-40 V0M // reproduce Astrids cuts with opening angle cut
    cuts.AddCutCalo("14500023","11111020530a2230000","01631031000000d0"); // 40-50 V0M // reproduce Astrids cuts with opening angle cut
    cuts.AddCutCalo("15600023","11111020530a2230000","01631031000000d0"); // 50-60 V0M // reproduce Astrids cuts with opening angle cut
  } else if (trainConfig == 65){
    cuts.AddCutCalo("16700023","11111020530a2230000","01631031000000d0"); // 60-70 V0M // reproduce Astrids cuts with opening angle cut
    cuts.AddCutCalo("17800023","11111020530a2230000","01631031000000d0"); // 70-80 V0M // reproduce Astrids cuts with opening angle cut
    cuts.AddCutCalo("18900023","11111020530a2230000","01631031000000d0"); // 80-90 V0M // reproduce Astrids cuts with opening angle cut
    cuts.AddCutCalo("14600023","11111020530a2230000","01631031000000d0"); // 40-60 V0M // reproduce Astrids cuts with opening angle cut
    cuts.AddCutCalo("16800023","11111020530a2230000","01631031000000d0"); // 60-80 V0M // reproduce Astrids cuts with opening angle cut
  } else if (trainConfig == 66){
    cuts.AddCutCalo("16700023","11111020530a2230000","01631031000000d0"); // 60-70 V0M // reproduce Astrids cuts with opening angle cut
    cuts.AddCutCalo("17800023","11111020530a2230000","01631031000000d0"); // 70-80 V0M // reproduce Astrids cuts with opening angle cut
    cuts.AddCutCalo("18900023","11111020530a2230000","01631031000000d0"); // 80-90 V0M // reproduce Astrids cuts with opening angle cut
    cuts.AddCutCalo("14600023","11111020530a2230000","01631031000000d0"); // 40-60 V0M // reproduce Astrids cuts with opening angle cut
    cuts.AddCutCalo("16800023","11111020530a2230000","01631031000000d0"); // 60-80 V0M // reproduce Astrids cuts with opening angle cut

  } else if (trainConfig == 67){
    cuts.AddCutCalo("40100013","11111020530a2230000","01631031000000d0"); // 0-10 new track array cuts in MC  // reproduce Astrids cuts with opening angle cut
    cuts.AddCutCalo("41200013","11111020530a2230000","01631031000000d0"); // 10-20 new track array cuts in MC // reproduce Astrids cuts with opening angle cut
    cuts.AddCutCalo("42300013","11111020530a2230000","01631031000000d0"); // 20-30 new track array cuts in MC // reproduce Astrids cuts with opening angle cut
  } else if (trainConfig == 68){
    cuts.AddCutCalo("70100013","11111020530a2230000","01631031000000d0"); // 0-5 new track array cuts in MC   // reproduce Astrids cuts with opening angle cut
    cuts.AddCutCalo("71200013","11111020530a2230000","01631031000000d0"); // 5-10 new track array cuts in MC  // reproduce Astrids cuts with opening angle cut
    cuts.AddCutCalo("72300013","11111020530a2230000","01631031000000d0"); // 10-15 new track array cuts in MC  // reproduce Astrids cuts with opening angle cut
    cuts.AddCutCalo("73400013","11111020530a2230000","01631031000000d0"); // 15-20 new track array cuts in MC  // reproduce Astrids cuts with opening angle cut
  } else if (trainConfig == 69){
    cuts.AddCutCalo("43400013","11111020530a2230000","01631031000000d0"); // 30-40 new track array cuts in MC  // reproduce Astrids cuts with opening angle cut
    cuts.AddCutCalo("44500013","11111020530a2230000","01631031000000d0"); // 40-50 new track array cuts in MC // reproduce Astrids cuts with opening angle cut
    cuts.AddCutCalo("45600013","11111020530a2230000","01631031000000d0"); // 50-60 new track array cuts in MC // reproduce Astrids cuts with opening angle cut
  } else if (trainConfig == 70){
    cuts.AddCutCalo("46700013","11111020530a2230000","01631031000000d0"); // 60-70 new track array cuts in MC  // reproduce Astrids cuts with opening angle cut
    cuts.AddCutCalo("47800013","11111020530a2230000","01631031000000d0"); // 70-80 new track array cuts in MC // reproduce Astrids cuts with opening angle cut
    cuts.AddCutCalo("48900013","11111020530a2230000","01631031000000d0"); // 80-90 new track array cuts in MC // reproduce Astrids cuts with opening angle cut

  } else if (trainConfig == 71){ // added signals
    cuts.AddCutCalo("40100023","11111020530a2230000","01631031000000d0"); // 0-10 new track array cuts in MC  // reproduce Astrids cuts with opening angle cut
    cuts.AddCutCalo("41200023","11111020530a2230000","01631031000000d0"); // 10-20 new track array cuts in MC // reproduce Astrids cuts with opening angle cut
    cuts.AddCutCalo("42300023","11111020530a2230000","01631031000000d0"); // 20-30 new track array cuts in MC // reproduce Astrids cuts with opening angle cut
  } else if (trainConfig == 72){ // added signals
    cuts.AddCutCalo("70100023","11111020530a2230000","01631031000000d0"); // 0-5 new track array cuts in MC   // reproduce Astrids cuts with opening angle cut
    cuts.AddCutCalo("71200023","11111020530a2230000","01631031000000d0"); // 5-10 new track array cuts in MC  // reproduce Astrids cuts with opening angle cut
    cuts.AddCutCalo("72300023","11111020530a2230000","01631031000000d0"); // 10-15 new track array cuts in MC  // reproduce Astrids cuts with opening angle cut
    cuts.AddCutCalo("73400023","11111020530a2230000","01631031000000d0"); // 15-20 new track array cuts in MC  // reproduce Astrids cuts with opening angle cut
  } else if (trainConfig == 73){ // added signals
    cuts.AddCutCalo("43400023","11111020530a2230000","01631031000000d0"); // 30-40 new track array cuts in MC  // reproduce Astrids cuts with opening angle cut
    cuts.AddCutCalo("44500023","11111020530a2230000","01631031000000d0"); // 40-50 new track array cuts in MC // reproduce Astrids cuts with opening angle cut
    cuts.AddCutCalo("45600023","11111020530a2230000","01631031000000d0"); // 60-70 new track array cuts in MC // reproduce Astrids cuts with opening angle cut
  } else if (trainConfig == 74){ // added signals
    cuts.AddCutCalo("46700023","11111020530a2230000","01631031000000d0"); // 60-70 new track array cuts in MC  // reproduce Astrids cuts with opening angle cut
    cuts.AddCutCalo("47800023","11111020530a2230000","01631031000000d0"); // 70-80 new track array cuts in MC // reproduce Astrids cuts with opening angle cut
    cuts.AddCutCalo("48900023","11111020530a2230000","01631031000000d0"); // 80-90 new track array cuts in MC // reproduce Astrids cuts with opening angle cut
  } else if (trainConfig == 75){ // added signals duplicate
    cuts.AddCutCalo("40100023","11111020530a2230000","01631031000000d0"); // 0-10 new track array cuts in MC  // reproduce Astrids cuts with opening angle cut
    cuts.AddCutCalo("41200023","11111020530a2230000","01631031000000d0"); // 10-20 new track array cuts in MC // reproduce Astrids cuts with opening angle cut
    cuts.AddCutCalo("42300023","11111020530a2230000","01631031000000d0"); // 20-30 new track array cuts in MC // reproduce Astrids cuts with opening angle cut
  } else if (trainConfig == 76){ // added signals duplicate
    cuts.AddCutCalo("70100023","11111020530a2230000","01631031000000d0"); // 0-5 new track array cuts in MC   // reproduce Astrids cuts with opening angle cut
    cuts.AddCutCalo("71200023","11111020530a2230000","01631031000000d0"); // 5-10 new track array cuts in MC  // reproduce Astrids cuts with opening angle cut
    cuts.AddCutCalo("72300023","11111020530a2230000","01631031000000d0"); // 10-15 new track array cuts in MC  // reproduce Astrids cuts with opening angle cut
    cuts.AddCutCalo("73400023","11111020530a2230000","01631031000000d0"); // 15-20 new track array cuts in MC  // reproduce Astrids cuts with opening angle cut
  } else if (trainConfig == 78){ // added signals duplicate
    cuts.AddCutCalo("43400023","11111020530a2230000","01631031000000d0"); // 30-40 new track array cuts in MC  // reproduce Astrids cuts with opening angle cut
    cuts.AddCutCalo("44500023","11111020530a2230000","01631031000000d0"); // 40-50 new track array cuts in MC // reproduce Astrids cuts with opening angle cut
    cuts.AddCutCalo("45600023","11111020530a2230000","01631031000000d0"); // 60-70 new track array cuts in MC // reproduce Astrids cuts with opening angle cut
  } else if (trainConfig == 79){ // added signals duplicate
    cuts.AddCutCalo("46700023","11111020530a2230000","01631031000000d0"); // 60-70 new track array cuts in MC  // reproduce Astrids cuts with opening angle cut
    cuts.AddCutCalo("47800023","11111020530a2230000","01631031000000d0"); // 70-80 new track array cuts in MC // reproduce Astrids cuts with opening angle cut
    cuts.AddCutCalo("48900023","11111020530a2230000","01631031000000d0"); // 80-90 new track array cuts in MC // reproduce Astrids cuts with opening angle cut


  // **********************************************************************************************************
  // ***************************** PHOS configurations run 1 **************************************************
  // **********************************************************************************************************
  } else if (trainConfig == 101){ // PHOS clusters
    cuts.AddCutCalo("60100013","2444400040033200000","0163103100000030"); // 0-5%
    cuts.AddCutCalo("61200013","2444400040033200000","0163103100000030"); // 5-10%
    cuts.AddCutCalo("50100013","2444400040033200000","0163103100000030"); // 0-10%
    cuts.AddCutCalo("52400013","2444400040033200000","0163103100000030"); // 20-40%
    cuts.AddCutCalo("52500013","2444400040033200000","0163103100000030"); // 20-50%
  } else if (trainConfig == 102){ // PHOS clusters
    cuts.AddCutCalo("60100013","2444400040033200000","0163103100000030"); // 0-5%
    cuts.AddCutCalo("61200013","2444400040033200000","0163103100000030"); // 5-10%
    cuts.AddCutCalo("50100013","2444400040033200000","0163103100000030"); // 0-10%
    cuts.AddCutCalo("51200013","2444400040033200000","0163103100000030"); // 10-20%
    cuts.AddCutCalo("52400013","2444400040033200000","0163103100000030"); // 20-40%
  } else if (trainConfig == 103){ // PHOS clusters
    cuts.AddCutCalo("54600013","2444400040033200000","0163103100000030"); // 40-60%
    cuts.AddCutCalo("56800013","2444400040033200000","0163103100000030"); // 60-80%
    cuts.AddCutCalo("52600013","2444400040033200000","0163103100000030"); // 20-60%
    cuts.AddCutCalo("54800013","2444400040033200000","0163103100000030"); // 40-80%
    cuts.AddCutCalo("52500013","2444400040033200000","0163103100000030"); // 20-50%
  } else if (trainConfig == 104){ // PHOS clusters with TM
    cuts.AddCutCalo("60100013","2444400042033200000","0163103100000030"); // 0-5%
    cuts.AddCutCalo("61200013","2444400042033200000","0163103100000030"); // 5-10%
    cuts.AddCutCalo("50100013","2444400042033200000","0163103100000030"); // 0-10%
    cuts.AddCutCalo("52400013","2444400042033200000","0163103100000030"); // 20-40%
    cuts.AddCutCalo("52500013","2444400042033200000","0163103100000030"); // 20-50%
  } else if (trainConfig == 105){ // PHOS clusters with TM
    cuts.AddCutCalo("60100013","2444400042033200000","0163103100000030"); // 0-5%
    cuts.AddCutCalo("61200013","2444400042033200000","0163103100000030"); // 5-10%
    cuts.AddCutCalo("50100013","2444400042033200000","0163103100000030"); // 0-10%
    cuts.AddCutCalo("51200013","2444400042033200000","0163103100000030"); // 10-20%
    cuts.AddCutCalo("52400013","2444400042033200000","0163103100000030"); // 20-40%
  } else if (trainConfig == 106){ // PHOS clusters with TM
    cuts.AddCutCalo("54600013","2444400042033200000","0163103100000030"); // 40-60%
    cuts.AddCutCalo("56800013","2444400042033200000","0163103100000030"); // 60-80%
    cuts.AddCutCalo("52600013","2444400042033200000","0163103100000030"); // 20-60%
    cuts.AddCutCalo("54800013","2444400042033200000","0163103100000030"); // 40-80%
    cuts.AddCutCalo("52500013","2444400042033200000","0163103100000030"); // 20-50%
  } else if (trainConfig == 107){ // PHOS clusters with TM
    cuts.AddCutCalo("50900013","2444400002033200000","0163103100000030"); // 0-90%


  // **********************************************************************************************************
  // ***************************** EMC configurations PbPb run 2 **********************************************
  // **********************************************************************************************************
  } else if (trainConfig == 201){ // EMCAL clusters central
    cuts.AddCutCalo("10110013","1111100057032230000","01631031000000d0"); // 0-10
    cuts.AddCutCalo("30110013","1111100057032230000","01631031000000d0"); // 0-5
    cuts.AddCutCalo("31210013","1111100057032230000","01631031000000d0"); // 5-10
  } else if (trainConfig == 202){ // EMCAL clusters semi-central
    cuts.AddCutCalo("11210013","1111100057032230000","01631031000000d0"); // 10-20
    cuts.AddCutCalo("12310013","1111100057032230000","01631031000000d0"); // 20-30
    cuts.AddCutCalo("13410013","1111100057032230000","01631031000000d0"); // 30-40
    cuts.AddCutCalo("12410013","1111100057032230000","01631031000000d0"); // 20-40
  } else if (trainConfig == 203){ // EMCAL clusters semi-central
    cuts.AddCutCalo("14510013","1111100057032230000","01631031000000d0"); // 40-50
    cuts.AddCutCalo("14610013","1111100057032230000","01631031000000d0"); // 40-60
    cuts.AddCutCalo("15610013","1111100057032230000","01631031000000d0"); // 50-60
  } else if (trainConfig == 204){ // EMCAL clusters peripheral
    cuts.AddCutCalo("16710013","1111100057032230000","01631031000000d0"); // 60-70
    cuts.AddCutCalo("17810013","1111100057032230000","01631031000000d0"); // 70-80
    cuts.AddCutCalo("18910013","1111100057032230000","01631031000000d0"); // 80-90
    cuts.AddCutCalo("16810013","1111100057032230000","01631031000000d0"); // 60-80

  } else if (trainConfig == 205){ // EMCAL clusters central add sig
    cuts.AddCutCalo("10110023","1111100057032230000","01631031000000d0"); // 0-10
    cuts.AddCutCalo("30110023","1111100057032230000","01631031000000d0"); // 0-5
    cuts.AddCutCalo("31210023","1111100057032230000","01631031000000d0"); // 5-10
  } else if (trainConfig == 206){ // EMCAL clusters semi-central add sig
    cuts.AddCutCalo("11210023","1111100057032230000","01631031000000d0"); // 10-20
    cuts.AddCutCalo("12310023","1111100057032230000","01631031000000d0"); // 20-30
    cuts.AddCutCalo("13410023","1111100057032230000","01631031000000d0"); // 30-40
    cuts.AddCutCalo("12410023","1111100057032230000","01631031000000d0"); // 20-40
  } else if (trainConfig == 207){ // EMCAL clusters semi-central add sig
    cuts.AddCutCalo("14510023","1111100057032230000","01631031000000d0"); // 40-50
    cuts.AddCutCalo("14610023","1111100057032230000","01631031000000d0"); // 40-60
    cuts.AddCutCalo("15610023","1111100057032230000","01631031000000d0"); // 50-60
  } else if (trainConfig == 208){ // EMCAL clusters peripheral add sig
    cuts.AddCutCalo("16710023","1111100057032230000","01631031000000d0"); // 60-70
    cuts.AddCutCalo("17810023","1111100057032230000","01631031000000d0"); // 70-80
    cuts.AddCutCalo("18910023","1111100057032230000","01631031000000d0"); // 80-90
    cuts.AddCutCalo("16810023","1111100057032230000","01631031000000d0"); // 60-80

  } else if (trainConfig == 209){ // EMCAL clusters - correction convcalo f1
    cuts.AddCutCalo("10110013","1111181053032230000","0163103100000050"); // 0-10
    cuts.AddCutCalo("11210013","1111181053032230000","0163103100000050"); // 10-20
    cuts.AddCutCalo("12510013","1111181053032230000","0163103100000050"); // 20-50
    cuts.AddCutCalo("15910013","1111181053032230000","0163103100000050"); // 50-90
  } else if (trainConfig == 210){ // EMCAL clusters - correction calocalo f2
    cuts.AddCutCalo("10110013","1111192053032230000","0163103100000050"); // 0-10
    cuts.AddCutCalo("11210013","1111192053032230000","0163103100000050"); // 10-20
    cuts.AddCutCalo("12510013","1111192053032230000","0163103100000050"); // 20-50
    cuts.AddCutCalo("15910013","1111192053032230000","0163103100000050"); // 50-90
  } else if (trainConfig == 211){ // EMCAL clusters - 0-90% centrality for PbPb EMCal cluster QA
    cuts.AddCutCalo("10910013","1111100003032230000","0163103100000050"); // 0-90

  } else if (trainConfig == 212){ // EMCAL clusters - 0-90% centrality for PbPb EMCal
    cuts.AddCutCalo("50910113","1111183053032230000","0163103100000050"); //  0-90 calo correction cent dep
  } else if (trainConfig == 213){ // EMCAL clusters - centrality selection for PbPb EMCal
    cuts.AddCutCalo("50110113","1111184053032230000","0163103100000050"); //  0-10 calo correction cent dep
    cuts.AddCutCalo("51210113","1111185053032230000","0163103100000050"); // 10-20 calo correction cent dep
    cuts.AddCutCalo("52510113","1111186053032230000","0163103100000050"); // 20-50 calo correction cent dep
    cuts.AddCutCalo("55910113","1111187053032230000","0163103100000050"); // 50-90 calo correction cent dep

  // trainconfig for PbPb studies in 5 TeV with TB nonlin
  } else if (trainConfig == 214){
    cuts.AddCutCalo("10110113","11111020530a2230000","01631031000000d0"); // 0-10 V0M  // reproduce Astrids cuts with opening angle cut
    cuts.AddCutCalo("11210113","11111020530a2230000","01631031000000d0"); // 10-20 V0M // reproduce Astrids cuts with opening angle cut
  } else if (trainConfig == 215){
    cuts.AddCutCalo("30110113","11111020530a2230000","01631031000000d0"); // 0-5 V0M   // reproduce Astrids cuts with opening angle cut
    cuts.AddCutCalo("31210113","11111020530a2230000","01631031000000d0"); // 5-10 V0M  // reproduce Astrids cuts with opening angle cut
  } else if (trainConfig == 216){
    cuts.AddCutCalo("12310113","11111020530a2230000","01631031000000d0"); // 20-30 V0M // reproduce Astrids cuts with opening angle cut
    cuts.AddCutCalo("13410113","11111020530a2230000","01631031000000d0"); // 30-40 V0M // reproduce Astrids cuts with opening angle cut
    cuts.AddCutCalo("14510113","11111020530a2230000","01631031000000d0"); // 40-50 V0M // reproduce Astrids cuts with opening angle cut
    cuts.AddCutCalo("15610113","11111020530a2230000","01631031000000d0"); // 50-60 V0M // reproduce Astrids cuts with opening angle cut
  } else if (trainConfig == 217){
    cuts.AddCutCalo("16710113","11111020530a2230000","01631031000000d0"); // 60-70 V0M // reproduce Astrids cuts with opening angle cut
    cuts.AddCutCalo("17810113","11111020530a2230000","01631031000000d0"); // 70-80 V0M // reproduce Astrids cuts with opening angle cut
    cuts.AddCutCalo("18910113","11111020530a2230000","01631031000000d0"); // 80-90 V0M // reproduce Astrids cuts with opening angle cut
    cuts.AddCutCalo("14610113","11111020530a2230000","01631031000000d0"); // 40-60 V0M // reproduce Astrids cuts with opening angle cut
    cuts.AddCutCalo("16810113","11111020530a2230000","01631031000000d0"); // 60-80 V0M // reproduce Astrids cuts with opening angle cut
  // trainconfig for PbPb studies in 5 TeV with no  nonlin
  } else if (trainConfig == 218){
    cuts.AddCutCalo("10110113","11111000530a2230000","01631031000000d0"); // 0-10 V0M  // reproduce Astrids cuts with opening angle cut
    cuts.AddCutCalo("11210113","11111000530a2230000","01631031000000d0"); // 10-20 V0M // reproduce Astrids cuts with opening angle cut
  } else if (trainConfig == 219){
    cuts.AddCutCalo("30110113","11111000530a2230000","01631031000000d0"); // 0-5 V0M   // reproduce Astrids cuts with opening angle cut
    cuts.AddCutCalo("31210113","11111000530a2230000","01631031000000d0"); // 5-10 V0M  // reproduce Astrids cuts with opening angle cut
  } else if (trainConfig == 220){
    cuts.AddCutCalo("12310113","11111000530a2230000","01631031000000d0"); // 20-30 V0M // reproduce Astrids cuts with opening angle cut
    cuts.AddCutCalo("13410113","11111000530a2230000","01631031000000d0"); // 30-40 V0M // reproduce Astrids cuts with opening angle cut
    cuts.AddCutCalo("14510113","11111000530a2230000","01631031000000d0"); // 40-50 V0M // reproduce Astrids cuts with opening angle cut
    cuts.AddCutCalo("15610113","11111000530a2230000","01631031000000d0"); // 50-60 V0M // reproduce Astrids cuts with opening angle cut
  } else if (trainConfig == 221){
    cuts.AddCutCalo("16710113","11111000530a2230000","01631031000000d0"); // 60-70 V0M // reproduce Astrids cuts with opening angle cut
    cuts.AddCutCalo("17810113","11111000530a2230000","01631031000000d0"); // 70-80 V0M // reproduce Astrids cuts with opening angle cut
    cuts.AddCutCalo("18910113","11111000530a2230000","01631031000000d0"); // 80-90 V0M // reproduce Astrids cuts with opening angle cut
    cuts.AddCutCalo("14610113","11111000530a2230000","01631031000000d0"); // 40-60 V0M // reproduce Astrids cuts with opening angle cut
    cuts.AddCutCalo("16810113","11111000530a2230000","01631031000000d0"); // 60-80 V0M // reproduce Astrids cuts with opening angle cut
    // trainconfig for PbPb studies in 5 TeV with TB nonlin
  } else if (trainConfig == 222){
    cuts.AddCutCalo("50110113","11111020530a2230000","01631031000000d0"); // 0-10 TM  // reproduce Astrids cuts with opening angle cut
    cuts.AddCutCalo("51210113","11111020530a2230000","01631031000000d0"); // 10-20 TM // reproduce Astrids cuts with opening angle cut
  } else if (trainConfig == 223){
    cuts.AddCutCalo("60110113","11111020530a2230000","01631031000000d0"); // 0-5 TM   // reproduce Astrids cuts with opening angle cut
    cuts.AddCutCalo("61110113","11111020530a2230000","01631031000000d0"); // 5-10 TM  // reproduce Astrids cuts with opening angle cut
  } else if (trainConfig == 224){
    cuts.AddCutCalo("52310113","11111020530a2230000","01631031000000d0"); // 20-30 TM // reproduce Astrids cuts with opening angle cut
    cuts.AddCutCalo("53410113","11111020530a2230000","01631031000000d0"); // 30-40 TM // reproduce Astrids cuts with opening angle cut
    cuts.AddCutCalo("54510113","11111020530a2230000","01631031000000d0"); // 40-50 TM // reproduce Astrids cuts with opening angle cut
    cuts.AddCutCalo("55610113","11111020530a2230000","01631031000000d0"); // 50-60 TM // reproduce Astrids cuts with opening angle cut
  } else if (trainConfig == 225){
    cuts.AddCutCalo("56710113","11111020530a2230000","01631031000000d0"); // 60-70 TM // reproduce Astrids cuts with opening angle cut
    cuts.AddCutCalo("57810113","11111020530a2230000","01631031000000d0"); // 70-80 TM // reproduce Astrids cuts with opening angle cut
    cuts.AddCutCalo("58910113","11111020530a2230000","01631031000000d0"); // 80-90 TM // reproduce Astrids cuts with opening angle cut
    cuts.AddCutCalo("54610113","11111020530a2230000","01631031000000d0"); // 40-60 TM // reproduce Astrids cuts with opening angle cut
    cuts.AddCutCalo("56810113","11111020530a2230000","01631031000000d0"); // 60-80 TM // reproduce Astrids cuts with opening angle cut

  } else if (trainConfig == 226){ // EMCAL clusters - peripheral centrality selection for PbPb EMCal
    cuts.AddCutCalo("15910113","11111020530a2230000","01631031000000d0"); //
    cuts.AddCutCalo("15910113","11111870530a2230000","01631031000000d0"); //
    cuts.AddCutCalo("15910113","11111020530b2230000","01631031000000d0"); //
    cuts.AddCutCalo("15910113","11111870530b2230000","01631031000000d0"); //

  } else if (trainConfig == 230){ // EMCAL clusters central add sig
    cuts.AddCutCalo("10110023","1111100057032230000","01631031000000d0"); // 0-10
    cuts.AddCutCalo("30110023","1111100057032230000","01631031000000d0"); // 0-5
    cuts.AddCutCalo("31210023","1111100057032230000","01631031000000d0"); // 5-10
  } else if (trainConfig == 231){ // EMCAL clusters semi-central add sig
    cuts.AddCutCalo("11210023","1111100057032230000","01631031000000d0"); // 10-20
    cuts.AddCutCalo("12310023","1111100057032230000","01631031000000d0"); // 20-30
    cuts.AddCutCalo("13410023","1111100057032230000","01631031000000d0"); // 30-40
    cuts.AddCutCalo("12410023","1111100057032230000","01631031000000d0"); // 20-40
  } else if (trainConfig == 232){ // EMCAL clusters semi-central add sig
    cuts.AddCutCalo("14510023","1111100057032230000","01631031000000d0"); // 40-50
    cuts.AddCutCalo("14610023","1111100057032230000","01631031000000d0"); // 40-60
    cuts.AddCutCalo("15610023","1111100057032230000","01631031000000d0"); // 50-60
  } else if (trainConfig == 233){ // EMCAL clusters peripheral add sig
    cuts.AddCutCalo("16710023","1111100057032230000","01631031000000d0"); // 60-70
    cuts.AddCutCalo("17810023","1111100057032230000","01631031000000d0"); // 70-80
    cuts.AddCutCalo("18910023","1111100057032230000","01631031000000d0"); // 80-90
    cuts.AddCutCalo("16810023","1111100057032230000","01631031000000d0"); // 60-80


  } else if (trainConfig == 234){ // EMCAL clusters central, TB
    cuts.AddCutCalo("10110113","1111102057032230000","01631031000000d0"); // 0-10
    cuts.AddCutCalo("30110113","1111102057032230000","01631031000000d0"); // 0-5
    cuts.AddCutCalo("31210113","1111102057032230000","01631031000000d0"); // 5-10
  } else if (trainConfig == 235){ // EMCAL clusters semi-central, TB
    cuts.AddCutCalo("11210113","1111102057032230000","01631031000000d0"); // 10-20
    cuts.AddCutCalo("12310113","1111102057032230000","01631031000000d0"); // 20-30
    cuts.AddCutCalo("13410113","1111102057032230000","01631031000000d0"); // 30-40
    cuts.AddCutCalo("12410113","1111102057032230000","01631031000000d0"); // 20-40
  } else if (trainConfig == 236){ // EMCAL clusters semi-central, TB
    cuts.AddCutCalo("14510113","1111102057032230000","01631031000000d0"); // 40-50
    cuts.AddCutCalo("14610113","1111102057032230000","01631031000000d0"); // 40-60
    cuts.AddCutCalo("15610113","1111102057032230000","01631031000000d0"); // 50-60
  } else if (trainConfig == 237){ // EMCAL clusters peripheral, TB
    cuts.AddCutCalo("16710113","1111102057032230000","01631031000000d0"); // 60-70
    cuts.AddCutCalo("17810113","1111102057032230000","01631031000000d0"); // 70-80
    cuts.AddCutCalo("18910113","1111102057032230000","01631031000000d0"); // 80-90
    cuts.AddCutCalo("16810113","1111102057032230000","01631031000000d0"); // 60-80


  } else if (trainConfig == 240){ // EMCAL clusters - 0-90% centrality for PbPb EMCal
    cuts.AddCutCalo("50910013","1111100053032230000","0163103100000050"); //  0-90 calo correction cent dep
    cuts.AddCutCalo("50910013","1111102053032230000","0163103100000050"); //  0-90 calo correction cent dep
    cuts.AddCutCalo("50910013","1111183053032230000","0163103100000050"); //  0-90 calo correction cent dep
  } else if (trainConfig == 241){ // EMCAL clusters - 0-90% centrality for PbPb EMCal
    cuts.AddCutCalo("50910613","1111100053032230000","0163103100000050"); //  0-90 calo correction cent dep
    cuts.AddCutCalo("50910613","1111102053032230000","0163103100000050"); //  0-90 calo correction cent dep
    cuts.AddCutCalo("50910613","1111183053032230000","0163103100000050"); //  0-90 calo correction cent dep
  } else if (trainConfig == 242){ // EMCAL clusters - centrality selection for PbPb EMCal
    cuts.AddCutCalo("50110013","1111100053032230000","0163103100000050"); //  0-10 calo correction cent dep
    cuts.AddCutCalo("51210013","1111100053032230000","0163103100000050"); // 10-20 calo correction cent dep
    cuts.AddCutCalo("52510013","1111100053032230000","0163103100000050"); // 20-50 calo correction cent dep
    cuts.AddCutCalo("55910013","1111100053032230000","0163103100000050"); // 50-90 calo correction cent dep
  } else if (trainConfig == 243){ // EMCAL clusters - centrality selection for PbPb EMCal
    cuts.AddCutCalo("50110613","1111100053032230000","0163103100000050"); //  0-10 calo correction cent dep
    cuts.AddCutCalo("51210613","1111100053032230000","0163103100000050"); // 10-20 calo correction cent dep
    cuts.AddCutCalo("52510613","1111100053032230000","0163103100000050"); // 20-50 calo correction cent dep
    cuts.AddCutCalo("55910613","1111100053032230000","0163103100000050"); // 50-90 calo correction cent dep
  } else if (trainConfig == 244){ // EMCAL clusters - centrality selection for PbPb EMCal
    cuts.AddCutCalo("50110613","1111102053032230000","0163103100000050"); //  0-10 calo correction cent dep
    cuts.AddCutCalo("51210613","1111102053032230000","0163103100000050"); // 10-20 calo correction cent dep
    cuts.AddCutCalo("52510613","1111102053032230000","0163103100000050"); // 20-50 calo correction cent dep
    cuts.AddCutCalo("55910613","1111102053032230000","0163103100000050"); // 50-90 calo correction cent dep
  } else if (trainConfig == 245){ // EMCAL clusters - centrality selection for PbPb EMCal
    cuts.AddCutCalo("50110613","1111184053032230000","0163103100000050"); //  0-10 calo correction cent dep
    cuts.AddCutCalo("51210613","1111185053032230000","0163103100000050"); // 10-20 calo correction cent dep
    cuts.AddCutCalo("52510613","1111186053032230000","0163103100000050"); // 20-50 calo correction cent dep
    cuts.AddCutCalo("55910613","1111187053032230000","0163103100000050"); // 50-90 calo correction cent dep
  } else if (trainConfig == 246){ // EMCAL clusters - 0-90% centrality for PbPb EMCal
    cuts.AddCutCalo("50910613","1111183050032230000","0163103100000050"); //  0-90 calo correction cent dep
    cuts.AddCutCalo("50910613","1111183051032230000","0163103100000050"); //  0-90 calo correction cent dep
    cuts.AddCutCalo("50910613","1111183057032230000","0163103100000050"); //  0-90 calo correction cent dep
  } else if (trainConfig == 247){ // EMCAL clusters - centrality selection for PbPb EMCal
    cuts.AddCutCalo("50110613","1111184050032230000","0163103100000050"); //  0-10 calo correction cent dep
    cuts.AddCutCalo("51210613","1111185050032230000","0163103100000050"); // 10-20 calo correction cent dep
    cuts.AddCutCalo("52510613","1111186050032230000","0163103100000050"); // 20-50 calo correction cent dep
    cuts.AddCutCalo("55910613","1111187050032230000","0163103100000050"); // 50-90 calo correction cent dep
  } else if (trainConfig == 248){ // EMCAL clusters - centrality selection for PbPb EMCal
    cuts.AddCutCalo("50110613","1111184051032230000","0163103100000050"); //  0-10 calo correction cent dep
    cuts.AddCutCalo("51210613","1111185051032230000","0163103100000050"); // 10-20 calo correction cent dep
    cuts.AddCutCalo("52510613","1111186051032230000","0163103100000050"); // 20-50 calo correction cent dep
    cuts.AddCutCalo("55910613","1111187051032230000","0163103100000050"); // 50-90 calo correction cent dep
  } else if (trainConfig == 249){ // EMCAL clusters - centrality selection for PbPb EMCal
    cuts.AddCutCalo("50110613","1111184057032230000","0163103100000050"); //  0-10 calo correction cent dep
    cuts.AddCutCalo("51210613","1111185057032230000","0163103100000050"); // 10-20 calo correction cent dep
    cuts.AddCutCalo("52510613","1111186057032230000","0163103100000050"); // 20-50 calo correction cent dep
    cuts.AddCutCalo("55910613","1111187057032230000","0163103100000050"); // 50-90 calo correction cent dep
  } else if (trainConfig == 250){ // EMCAL clusters - centrality selection for PbPb EMCal
    cuts.AddCutCalo("10910013","1111183051032230000","0163103100000050"); //  0-90 calo correction
    cuts.AddCutCalo("10910a13","1111183051032230000","0163103100000050"); //  0-90 calo correction
  } else if (trainConfig == 251){ // EMCAL clusters - centrality selection for PbPb EMCal
    cuts.AddCutCalo("10110013","1111183051032230000","0163103100000050"); //  0-90 calo correction
    cuts.AddCutCalo("11210013","1111183051032230000","0163103100000050"); //  0-90 calo correction
    cuts.AddCutCalo("12510013","1111183051032230000","0163103100000050"); //  0-90 calo correction
    cuts.AddCutCalo("15910013","1111183051032230000","0163103100000050"); //  0-90 calo correction
  } else if (trainConfig == 252){ // EMCAL clusters - centrality selection for PbPb EMCal
    cuts.AddCutCalo("10110a13","1111183051032230000","0163103100000050"); //  0-90 calo correction
    cuts.AddCutCalo("11210a13","1111183051032230000","0163103100000050"); //  0-90 calo correction
    cuts.AddCutCalo("12510a13","1111183051032230000","0163103100000050"); //  0-90 calo correction
    cuts.AddCutCalo("15910a13","1111183051032230000","0163103100000050"); //  0-90 calo correction
  } else if (trainConfig == 253){ // EMCAL clusters - 20180718 - default without corrections
    cuts.AddCutCalo("10110a13","1111100051032230000","0163103100000050"); //  0-90 calo correction
    cuts.AddCutCalo("11210a13","1111100051032230000","0163103100000050"); //  0-90 calo correction
    cuts.AddCutCalo("12410a13","1111100051032230000","0163103100000050"); //  0-90 calo correction
    cuts.AddCutCalo("14610a13","1111100051032230000","0163103100000050"); //  0-90 calo correction
    cuts.AddCutCalo("16810a13","1111100051032230000","0163103100000050"); //  0-90 calo correction
  } else if (trainConfig == 254){ // EMCAL clusters - 20180718 - default with peripheral corrections
    cuts.AddCutCalo("10110a13","1111187051032230000","0163103100000050"); //
    cuts.AddCutCalo("11210a13","1111187051032230000","0163103100000050"); //
    cuts.AddCutCalo("12410a13","1111187051032230000","0163103100000050"); //
    cuts.AddCutCalo("14610a13","1111187051032230000","0163103100000050"); //
    cuts.AddCutCalo("16810a13","1111187051032230000","0163103100000050"); //
  } else if (trainConfig == 255){ // EMCAL clusters - 20180718 - default with peripheral corrections
    cuts.AddCutCalo("12410a13","1111186051032230000","0163103100000050"); //
    cuts.AddCutCalo("12510a13","1111186051032230000","0163103100000050"); //
    cuts.AddCutCalo("14610a13","1111187051032230000","0163103100000050"); //
    cuts.AddCutCalo("16810a13","1111187051032230000","0163103100000050"); //
  } else if (trainConfig == 256){ // EMCAL clusters - 20180718 - default with corrections
    cuts.AddCutCalo("30110a13","1111184051032230000","0163103100000050"); //
    cuts.AddCutCalo("31210a13","1111184051032230000","0163103100000050"); //
    cuts.AddCutCalo("10110a13","1111184051032230000","0163103100000050"); //
    cuts.AddCutCalo("11210a13","1111185051032230000","0163103100000050"); //
  } else if (trainConfig == 257){ // EMCAL clusters - 20180718 - default with peripheral corrections
    cuts.AddCutCalo("14610a13","1111187050032230000","0163103100000050"); //
    cuts.AddCutCalo("14610a13","1111187051032230000","0163103100000050"); //
    cuts.AddCutCalo("14610a13","111118705i032230000","0163103100000050"); //
    cuts.AddCutCalo("14610a13","111118705j032230000","0163103100000050"); //
    cuts.AddCutCalo("14610a13","111118705k032230000","0163103100000050"); //
  } else if (trainConfig == 258){ // EMCAL clusters - TM studies with MIP subtraction
    cuts.AddCutCalo("10110a13","1111184050032230000","0163103100000050"); //
    cuts.AddCutCalo("10110a13","1111184051032230000","0163103100000050"); //
    cuts.AddCutCalo("10110a13","111118405i032230000","0163103100000050"); //
    cuts.AddCutCalo("10110a13","111118405j032230000","0163103100000050"); //
    cuts.AddCutCalo("10110a13","111118405k032230000","0163103100000050"); //
  } else if (trainConfig == 259){ // EMCAL clusters - studies for flat energy subtraction
    cuts.AddCutCalo("14610a13","1111187051032230000","0163103100000050"); //
    cuts.AddCutCalo("14610a13","11111870500f2230000","0163103100000050"); //
    cuts.AddCutCalo("10110a13","1111184051032230000","0163103100000050"); //
    cuts.AddCutCalo("10110a13","11111840500f2230000","0163103100000050"); //
  } else if (trainConfig == 260){ // EMCAL clusters - 20181110 - No NL
    cuts.AddCutCalo("12410a13","111110005k032230000","0163103100000050"); //
    cuts.AddCutCalo("12510a13","111110005k032230000","0163103100000050"); //
    cuts.AddCutCalo("14610a13","111110005k032230000","0163103100000050"); //
    cuts.AddCutCalo("16810a13","111110005k032230000","0163103100000050"); //
  } else if (trainConfig == 261){ // EMCAL clusters - 20181110 - No NL
    cuts.AddCutCalo("30110a13","111110005k0a2230000","0163103100000050"); //
    cuts.AddCutCalo("31210a13","111110005k0a2230000","0163103100000050"); //
    cuts.AddCutCalo("10110a13","111110005k0a2230000","0163103100000050"); //
    cuts.AddCutCalo("11210a13","111110005k0b2230000","0163103100000050"); //
  } else if (trainConfig == 262){ // EMCAL clusters - 20181207 - Peri NL
    cuts.AddCutCalo("12410a13","111118105k032230000","0163103100000050"); //
    cuts.AddCutCalo("12510a13","111118105k032230000","0163103100000050"); //
    cuts.AddCutCalo("14610a13","111118105k032230000","0163103100000050"); //
    cuts.AddCutCalo("16810a13","111118105k032230000","0163103100000050"); //
  } else if (trainConfig == 263){ // EMCAL clusters - 20181207 - Peri NL, asym cut on 0.1
    cuts.AddCutCalo("30110a13","111118105k032230000","016310d100000050"); //
    cuts.AddCutCalo("31210a13","111118105k032230000","016310d100000050"); //
    cuts.AddCutCalo("10110a13","111118105k032230000","016310d100000050"); //
    cuts.AddCutCalo("11210a13","111118105k032230000","016310d100000050"); //
  } else if (trainConfig == 264){ // EMCal + DCal clusters - 0-90% centrality
    cuts.AddCutCalo("10910a13","4117900050032230000","0163103100000050"); //
  } else if (trainConfig == 265){ // EMCAL clusters - 20190215
    cuts.AddCutCalo("12410a13","111118305k032220000","01631031000000d0"); //
    cuts.AddCutCalo("12510a13","111118305k032220000","01631031000000d0"); //
    cuts.AddCutCalo("14610a13","111118305k032220000","01631031000000d0"); //
    cuts.AddCutCalo("16810a13","111118305k032220000","01631031000000d0"); //
  } else if (trainConfig == 266){ // EMCAL clusters - 20190215, asym cut on 0.8
    cuts.AddCutCalo("30110a13","111118305k0a2220000","01631061000000d0"); //
    cuts.AddCutCalo("31210a13","111118305k0a2220000","01631061000000d0"); //
    cuts.AddCutCalo("10110a13","111118305k0a2220000","01631061000000d0"); //
    cuts.AddCutCalo("11210a13","111118305k0b2220000","01631061000000d0"); //
  } else if (trainConfig == 267){ // EMCAL + DCal clusters - 20190301
    cuts.AddCutCalo("12410a13","411798305k032220000","01631031000000d0"); //
    cuts.AddCutCalo("12510a13","411798305k032220000","01631031000000d0"); //
    cuts.AddCutCalo("14610a13","411798305k032220000","01631031000000d0"); //
    cuts.AddCutCalo("14810a13","411798305k032220000","01631031000000d0"); //
    cuts.AddCutCalo("16810a13","411798305k032220000","01631031000000d0"); //
  } else if (trainConfig == 268){ // EMCAL + DCal clusters - 20190301
    cuts.AddCutCalo("30110a13","411798305k0a2220000","01631061000000d0"); //
    cuts.AddCutCalo("31210a13","411798305k0a2220000","01631061000000d0"); //
    cuts.AddCutCalo("10110a13","411798305k0a2220000","01631061000000d0"); //
    cuts.AddCutCalo("11210a13","411798305k0b2220000","01631061000000d0"); //
    cuts.AddCutCalo("10210a13","411798305k0a2220000","01631061000000d0"); //

  //systematics for LHC15o 0-10%
  } else if (trainConfig == 270){ // EMCAL
    cuts.AddCutCalo("10110a13","111118105k032230000","016310d1000000d0"); //
  } else if (trainConfig == 271){ // EMCAL syst 1/7
    cuts.AddCutCalo("10110a13","111118105k022230000","016310d1000000d0"); // min energy cluster variation 1  600 MeV
    cuts.AddCutCalo("10110a13","111118105k042230000","016310d1000000d0"); // min energy cluster variation 2  800 MeV
    cuts.AddCutCalo("10110a13","111118105k052230000","016310d1000000d0"); // min energy cluster variation 3  900 MeV
    cuts.AddCutCalo("10110a13","111118105k032230000","01631031000000d0"); // alpha meson variation 1 0<alpha<1.0
    cuts.AddCutCalo("10110a13","111118105k032230000","01631061000000d0"); // alpha meson variation 1 0<alpha<0.8
  } else if (trainConfig == 272){ // EMCAL syst 2/7
    cuts.AddCutCalo("10110a13","111118105k032230000","016310d1000000d0"); // min/max M02  0.1<M<0.5
    cuts.AddCutCalo("10110a13","111118105k032200000","016310d1000000d0"); // min/max M02  0.1<M<100
    cuts.AddCutCalo("10110a13","111118105k032250000","016310d1000000d0"); // min/max M02  0.1<M<0.3
    cuts.AddCutCalo("10110a13","111118105k032260000","016310d1000000d0"); // min/max M02  0.1<M<0.27
  } else if (trainConfig == 273){ // EMCAL syst 3/7
    cuts.AddCutCalo("10110a13","111218105k032230000","016310d1000000d0"); // only modules with TRD infront
    cuts.AddCutCalo("10110a13","111138105k032230000","016310d1000000d0"); // no modules with TRD infront
    cuts.AddCutCalo("10110a13","111118105k032230000","016330d1000000d0"); // rapidity variation  y<0.6
    cuts.AddCutCalo("10110a13","111118105k032230000","016340d1000000d0"); // rapidity variation  y<0.5
  } else if (trainConfig == 274){ // EMCAL syst 4/7
    cuts.AddCutCalo("10110a13","111118105k032230000","016310d100000040"); // min opening angle 0.0152
    cuts.AddCutCalo("10110a13","111118105k032230000","016310d100000060"); // min opening angle 0.017
    cuts.AddCutCalo("10110a13","111118105k032230000","016310d100000070"); // min opening angle 0.016
    cuts.AddCutCalo("10110a13","111118105k032230000","016310d100000080"); // min opening angle 0.018
    cuts.AddCutCalo("10110a13","111118105k032230000","016310d100000090"); // min opening angle 0.019
  } else if (trainConfig == 275){ // EMCAL syst 5/7
    cuts.AddCutCalo("10110a13","1111181053032230000","016310d1000000d0"); // fixed window
    cuts.AddCutCalo("10110a13","1111181056032230000","016310d1000000d0"); // tm pt dependent var 1
    cuts.AddCutCalo("10110a13","1111181058032230000","016310d1000000d0"); // tm pt dependent var 2
    cuts.AddCutCalo("10110a13","1111181059032230000","016310d1000000d0"); // tm pt dependent var 3
  } else if (trainConfig == 276){ // EMCAL syst 6/7
    cuts.AddCutCalo("10110a13","111110105k032230000","016310d1000000d0"); // NL variation, TB
    cuts.AddCutCalo("10110a13","111118105k032230000","016310d1000000d0"); // NL variation, ConvCalo
    cuts.AddCutCalo("10110a13","111118205k032230000","016310d1000000d0"); // NL variation, Calo
  } else if (trainConfig == 277){ // EMCAL syst 7/7
    cuts.AddCutCalo("10110a13","111118104k032230000","016310d1000000d0"); // cluster timing cut
    cuts.AddCutCalo("10110a13","111118107k032230000","016310d1000000d0"); // cluster timing cut
    cuts.AddCutCalo("10110a13","11111810ak032230000","016310d1000000d0"); // cluster timing cut

  //systematics for LHC15o 20-40%
  } else if (trainConfig == 280){ // EMCAL
    cuts.AddCutCalo("12410a13","111118105k032230000","01631031000000d0"); //
  } else if (trainConfig == 281){ // EMCAL syst 1/7
    cuts.AddCutCalo("12410a13","111118105k022230000","01631031000000d0"); // min energy cluster variation 1  600 MeV
    cuts.AddCutCalo("12410a13","111118105k042230000","01631031000000d0"); // min energy cluster variation 2  800 MeV
    cuts.AddCutCalo("12410a13","111118105k052230000","01631031000000d0"); // min energy cluster variation 3  900 MeV
    cuts.AddCutCalo("12410a13","111118105k032230000","01631061000000d0"); // alpha meson variation 1 0<alpha<0.8
    cuts.AddCutCalo("12410a13","111118105k032230000","01631051000000d0"); // alpha meson variation 2 0<alpha<0.75
  } else if (trainConfig == 282){ // EMCAL syst 2/7
    cuts.AddCutCalo("12410a13","111118105k032230000","01631031000000d0"); // min/max M02  0.1<M<0.5
    cuts.AddCutCalo("12410a13","111118105k032200000","01631031000000d0"); // min/max M02  0.1<M<100
    cuts.AddCutCalo("12410a13","111118105k032250000","01631031000000d0"); // min/max M02  0.1<M<0.3
    cuts.AddCutCalo("12410a13","111118105k032260000","01631031000000d0"); // min/max M02  0.1<M<0.27
  } else if (trainConfig == 283){ // EMCAL syst 3/7
    cuts.AddCutCalo("12410a13","111218105k032230000","01631031000000d0"); // only modules with TRD infront
    cuts.AddCutCalo("12410a13","111138105k032230000","01631031000000d0"); // no modules with TRD infront
    cuts.AddCutCalo("12410a13","111118105k032230000","01633031000000d0"); // rapidity variation  y<0.6
    cuts.AddCutCalo("12410a13","111118105k032230000","01634031000000d0"); // rapidity variation  y<0.5
  } else if (trainConfig == 284){ // EMCAL syst 4/7
    cuts.AddCutCalo("12410a13","111118105k032230000","0163103100000040"); // min opening angle 0.0152
    cuts.AddCutCalo("12410a13","111118105k032230000","0163103100000060"); // min opening angle 0.017
    cuts.AddCutCalo("12410a13","111118105k032230000","0163103100000070"); // min opening angle 0.016
    cuts.AddCutCalo("12410a13","111118105k032230000","0163103100000080"); // min opening angle 0.018
    cuts.AddCutCalo("12410a13","111118105k032230000","0163103100000090"); // min opening angle 0.019
  } else if (trainConfig == 285){ // EMCAL syst 5/7
    cuts.AddCutCalo("12410a13","1111181053032230000","01631031000000d0"); // fixed window
    cuts.AddCutCalo("12410a13","1111181056032230000","01631031000000d0"); // tm pt dependent var 1
    cuts.AddCutCalo("12410a13","1111181058032230000","01631031000000d0"); // tm pt dependent var 2
    cuts.AddCutCalo("12410a13","1111181059032230000","01631031000000d0"); // tm pt dependent var 3
  } else if (trainConfig == 286){ // EMCAL syst 6/7
    cuts.AddCutCalo("12410a13","111110105k032230000","01631031000000d0"); // NL variation, TB
    cuts.AddCutCalo("12410a13","111118105k032230000","01631031000000d0"); // NL variation, ConvCalo
    cuts.AddCutCalo("12410a13","111118205k032230000","01631031000000d0"); // NL variation, Calo
  } else if (trainConfig == 287){ // EMCAL syst 7/7
    cuts.AddCutCalo("12410a13","111118104k032230000","01631031000000d0"); // cluster timing cut
    cuts.AddCutCalo("12410a13","111118107k032230000","01631031000000d0"); // cluster timing cut
    cuts.AddCutCalo("12410a13","11111810ak032230000","01631031000000d0"); // cluster timing cut
  } else if (trainConfig == 288){ // EMCAL M02 studies
    cuts.AddCutCalo("12410a13","111118105k032230000","01631031000000d0"); // min/max M02  0.1<M<0.5
    cuts.AddCutCalo("12410a13","111118105k032200000","01631031000000d0"); // min/max M02  0.1<M<1000
    cuts.AddCutCalo("12410a13","111118105k032240000","01631031000000d0"); // min/max M02  0.1<M<0.4
    cuts.AddCutCalo("12410a13","111118105k032220000","01631031000000d0"); // min/max M02  0.1<M<0.7

  //systematics for LHC15o 60-80%
  } else if (trainConfig == 290){ // EMCAL
    cuts.AddCutCalo("16810a13","111118105k032230000","01631031000000d0"); //
  } else if (trainConfig == 291){ // EMCAL syst 1/7
    cuts.AddCutCalo("16810a13","111118105k022230000","01631031000000d0"); // min energy cluster variation 1  600 MeV
    cuts.AddCutCalo("16810a13","111118105k042230000","01631031000000d0"); // min energy cluster variation 2  800 MeV
    cuts.AddCutCalo("16810a13","111118105k052230000","01631031000000d0"); // min energy cluster variation 3  900 MeV
    cuts.AddCutCalo("16810a13","111118105k032230000","01631061000000d0"); // alpha meson variation 1 0<alpha<0.8
    cuts.AddCutCalo("16810a13","111118105k032230000","01631051000000d0"); // alpha meson variation 2 0<alpha<0.75
  } else if (trainConfig == 292){ // EMCAL syst 2/7
    cuts.AddCutCalo("16810a13","111118105k032230000","01631031000000d0"); // min/max M02  0.1<M<0.5
    cuts.AddCutCalo("16810a13","111118105k032200000","01631031000000d0"); // min/max M02  0.1<M<100
    cuts.AddCutCalo("16810a13","111118105k032250000","01631031000000d0"); // min/max M02  0.1<M<0.3
    cuts.AddCutCalo("16810a13","111118105k032260000","01631031000000d0"); // min/max M02  0.1<M<0.27
  } else if (trainConfig == 293){ // EMCAL syst 3/7
    cuts.AddCutCalo("16810a13","111218105k032230000","01631031000000d0"); // only modules with TRD infront
    cuts.AddCutCalo("16810a13","111138105k032230000","01631031000000d0"); // no modules with TRD infront
    cuts.AddCutCalo("16810a13","111118105k032230000","01633031000000d0"); // rapidity variation  y<0.6
    cuts.AddCutCalo("16810a13","111118105k032230000","01634031000000d0"); // rapidity variation  y<0.5
  } else if (trainConfig == 294){ // EMCAL syst 4/7
    cuts.AddCutCalo("16810a13","111118105k032230000","0163103100000040"); // min opening angle 0.0152
    cuts.AddCutCalo("16810a13","111118105k032230000","0163103100000060"); // min opening angle 0.017
    cuts.AddCutCalo("16810a13","111118105k032230000","0163103100000070"); // min opening angle 0.016
    cuts.AddCutCalo("16810a13","111118105k032230000","0163103100000080"); // min opening angle 0.018
    cuts.AddCutCalo("16810a13","111118105k032230000","0163103100000090"); // min opening angle 0.019
  } else if (trainConfig == 295){ // EMCAL syst 5/7
    cuts.AddCutCalo("16810a13","1111181053032230000","01631031000000d0"); // fixed window
    cuts.AddCutCalo("16810a13","1111181056032230000","01631031000000d0"); // tm pt dependent var 1
    cuts.AddCutCalo("16810a13","1111181058032230000","01631031000000d0"); // tm pt dependent var 2
    cuts.AddCutCalo("16810a13","1111181059032230000","01631031000000d0"); // tm pt dependent var 3
  } else if (trainConfig == 296){ // EMCAL syst 6/7
    cuts.AddCutCalo("16810a13","111110105k032230000","01631031000000d0"); // NL variation, TB
    cuts.AddCutCalo("16810a13","111118105k032230000","01631031000000d0"); // NL variation, ConvCalo
    cuts.AddCutCalo("16810a13","111118205k032230000","01631031000000d0"); // NL variation, Calo
  } else if (trainConfig == 297){ // EMCAL syst 7/7
    cuts.AddCutCalo("16810a13","111118104k032230000","01631031000000d0"); // cluster timing cut
    cuts.AddCutCalo("16810a13","111118107k032230000","01631031000000d0"); // cluster timing cut
    cuts.AddCutCalo("16810a13","11111810ak032230000","01631031000000d0"); // cluster timing cut

  // **********************************************************************************************************
  // ***************************** EMC configurations XeXe run 2 **********************************************
  // **********************************************************************************************************
  } else if (trainConfig == 300){ // EMCAL clusters - 0-80% centrality for XeXe EMCal cluster QA
    cuts.AddCutCalo("10810013","411790000k032220000","01631031000000d0"); // 0-80
  } else if (trainConfig == 301){ // EMCAL clusters - 0-80% centrality for XeXe EMCal cluster QA
    cuts.AddCutCalo("10210013","411790000k032220000","01631031000000d0"); // 0-20
    cuts.AddCutCalo("12410013","411790000k032220000","01631031000000d0"); // 20-40
    cuts.AddCutCalo("10410013","411790000k032220000","01631031000000d0"); // 0-40
    cuts.AddCutCalo("14810013","411790000k032220000","01631031000000d0"); // 40-80
  } else if (trainConfig == 302){ // EMCAL clusters - 0-80% centrality for XeXe EMCal cluster QA - TB nl
    cuts.AddCutCalo("10810013","411790105k032220000","01631031000000d0"); // 0-80
  } else if (trainConfig == 303){ // EMCAL clusters - diff centralities for XeXe EMCal cluster QA - TB nl
    cuts.AddCutCalo("10210013","411790105k032220000","01631031000000d0"); // 0-20
    cuts.AddCutCalo("12410013","411790105k032220000","01631031000000d0"); // 20-40
    cuts.AddCutCalo("10410013","411790105k032220000","01631031000000d0"); // 0-40
    cuts.AddCutCalo("14810013","411790105k032220000","01631031000000d0"); // 40-80
  } else if (trainConfig == 304){ // EMCAL clusters - 0-80% centrality for XeXe EMCal cluster QA - TB nl + pp PCM-EMC tuning
    cuts.AddCutCalo("10810013","411793105k032220000","01631031000000d0"); // 0-80
  } else if (trainConfig == 305){ // EMCAL clusters - diff centralities for XeXe EMCal cluster QA - TB nl + pp PCM-EMC tuning
    cuts.AddCutCalo("10210013","411793105k032220000","01631031000000d0"); // 0-20
    cuts.AddCutCalo("12410013","411793105k032220000","01631031000000d0"); // 20-40
    cuts.AddCutCalo("10410013","411793105k032220000","01631031000000d0"); // 0-40
    cuts.AddCutCalo("14810013","411793105k032220000","01631031000000d0"); // 40-80

  // **********************************************************************************************************
  // ***************************** PHOS configurations XeXe run 2 *********************************************
  // **********************************************************************************************************
  } else if (trainConfig == 400){ // PHOS clusters - 0-80% centrality for XeXe PHOS cluster QA
    cuts.AddCutCalo("10810013","2446600000012200000","0163103100000010"); // 0-80
  } else if (trainConfig == 401){ // PHOS clusters - centrality for XeXe PHOS cluster QA
    cuts.AddCutCalo("10210013","2446600000012200000","0163103100000010"); // 0-20
    cuts.AddCutCalo("12410013","2446600000012200000","0163103100000010"); // 20-40
    cuts.AddCutCalo("10410013","2446600000012200000","0163103100000010"); // 0-40
    cuts.AddCutCalo("14810013","2446600000012200000","0163103100000010"); // 40-80
  } else if (trainConfig == 403){ // PHOS clusters - 0-80% centrality for XeXe PHOS cluster QA, PHOS TM sigma 2
    cuts.AddCutCalo("10810013","2446600007012200000","0163103100000010"); // 0-80
  } else if (trainConfig == 404){ // PHOS clusters - centrality for XeXe PHOS cluster QA
    cuts.AddCutCalo("10210013","2446600007012200000","0163103100000010"); // 0-20
    cuts.AddCutCalo("12410013","2446600007012200000","0163103100000010"); // 20-40
    cuts.AddCutCalo("10410013","2446600007012200000","0163103100000010"); // 0-40
    cuts.AddCutCalo("14810013","2446600007012200000","0163103100000010"); // 40-80
  } else if (trainConfig == 405){ // PHOS clusters - 0-80% centrality for XeXe PHOS cluster QA, PHOS TM sigma 2.5
    cuts.AddCutCalo("10810013","2446600008012200000","0163103100000010"); // 0-80
  } else if (trainConfig == 406){ // PHOS clusters - centrality for XeXe PHOS cluster QA
    cuts.AddCutCalo("10210013","2446600008012200000","0163103100000010"); // 0-20
    cuts.AddCutCalo("12410013","2446600008012200000","0163103100000010"); // 20-40
    cuts.AddCutCalo("10410013","2446600008012200000","0163103100000010"); // 0-40
    cuts.AddCutCalo("14810013","2446600008012200000","0163103100000010"); // 40-80
  } else if (trainConfig == 407){ // PHOS clusters - centrality for XeXe PHOS cluster QA
    cuts.AddCutCalo("10210013","2446600004012200000","0163103100000010"); // 0-20
    cuts.AddCutCalo("12410013","2446600004012200000","0163103100000010"); // 20-40
    cuts.AddCutCalo("10410013","2446600004012200000","0163103100000010"); // 0-40
    cuts.AddCutCalo("14810013","2446600004012200000","0163103100000010"); // 40-80


  // **********************************************************************************************************
  // ***************************** PHOS configurations PbPb run 2**********************************************
  // **********************************************************************************************************
  } else if (trainConfig == 610){ // PHOS clusters - 20181018
    cuts.AddCutCalo("12410a13","2446600051012200000","0163103100000010"); //
    cuts.AddCutCalo("12510a13","2446600051012200000","0163103100000010"); //
    cuts.AddCutCalo("14610a13","2446600051012200000","0163103100000010"); //
    cuts.AddCutCalo("14810a13","2446600051012200000","0163103100000010"); //
    cuts.AddCutCalo("16810a13","2446600051012200000","0163103100000010"); //
  } else if (trainConfig == 611){ // PHOS clusters - 20181018
    cuts.AddCutCalo("30110a13","2446600051012200000","0163103100000010"); //
    cuts.AddCutCalo("31210a13","2446600051012200000","0163103100000010"); //
    cuts.AddCutCalo("10110a13","2446600051012200000","0163103100000010"); //
    cuts.AddCutCalo("11210a13","2446600051012200000","0163103100000010"); //
    cuts.AddCutCalo("10210a13","2446600051012200000","0163103100000010"); //
  } else if (trainConfig == 612){ // PHOS clusters - 20181018
    cuts.AddCutCalo("12410a13","2446600054012200000","0163103100000010"); //
    cuts.AddCutCalo("12510a13","2446600054012200000","0163103100000010"); //
    cuts.AddCutCalo("14610a13","2446600054012200000","0163103100000010"); //
    cuts.AddCutCalo("14810a13","2446600054012200000","0163103100000010"); //
    cuts.AddCutCalo("16810a13","2446600054012200000","0163103100000010"); //
  } else if (trainConfig == 613){ // PHOS clusters - 20181018
    cuts.AddCutCalo("30110a13","2446600054012200000","0163103100000010"); //
    cuts.AddCutCalo("31210a13","2446600054012200000","0163103100000010"); //
    cuts.AddCutCalo("10110a13","2446600054012200000","0163103100000010"); //
    cuts.AddCutCalo("11210a13","2446600054012200000","0163103100000010"); //
    cuts.AddCutCalo("10210a13","2446600054012200000","0163103100000010"); //
  // **********************************************************************************************************
  // ***************************** EMC configurations PbPb run 2 EMC trigger **********************************************
  // **********************************************************************************************************
  } else if (trainConfig == 701){ // EMCAL clusters central
    cuts.AddCutCalo("10183013","1111100057032230000","01631031000000d0"); // 0-10
    cuts.AddCutCalo("30183013","1111100057032230000","01631031000000d0"); // 0-5
    cuts.AddCutCalo("31283013","1111100057032230000","01631031000000d0"); // 5-10
  } else if (trainConfig == 702){ // EMCAL clusters semi-central
    cuts.AddCutCalo("11283013","1111100057032230000","01631031000000d0"); // 10-20
    cuts.AddCutCalo("12383013","1111100057032230000","01631031000000d0"); // 20-30
    cuts.AddCutCalo("13483013","1111100057032230000","01631031000000d0"); // 30-40
    cuts.AddCutCalo("12483013","1111100057032230000","01631031000000d0"); // 20-40
  } else if (trainConfig == 703){ // EMCAL clusters semi-central
    cuts.AddCutCalo("14583013","1111100057032230000","01631031000000d0"); // 40-50
    cuts.AddCutCalo("14683013","1111100057032230000","01631031000000d0"); // 40-60
    cuts.AddCutCalo("15683013","1111100057032230000","01631031000000d0"); // 50-60
  } else if (trainConfig == 704){ // EMCAL clusters peripheral
    cuts.AddCutCalo("16783013","1111100057032230000","01631031000000d0"); // 60-70
    cuts.AddCutCalo("17883013","1111100057032230000","01631031000000d0"); // 70-80
    cuts.AddCutCalo("18983013","1111100057032230000","01631031000000d0"); // 80-90
    cuts.AddCutCalo("16883013","1111100057032230000","01631031000000d0"); // 60-80

  //EG2 + DG2
  } else if (trainConfig == 710){ // EMCAL + DCal clusters - 20190301
    cuts.AddCutCalo("1248ea13","411798305k032220000","01631031000000d0"); //
    cuts.AddCutCalo("1258ea13","411798305k032220000","01631031000000d0"); //
    cuts.AddCutCalo("1468ea13","411798305k032220000","01631031000000d0"); //
    cuts.AddCutCalo("1688ea13","411798305k032220000","01631031000000d0"); //
  } else if (trainConfig == 711){ // EMCAL + DCal clusters - 20190301
    cuts.AddCutCalo("1018ea13","411798305k0a2220000","01631061000000d0"); //
    cuts.AddCutCalo("1128ea13","411798305k0b2220000","01631061000000d0"); //
  //EG1 + DG1
  } else if (trainConfig == 712){ // EMCAL + DCal clusters - 20190301
    cuts.AddCutCalo("1248da13","411798305k032220000","01631031000000d0"); //
    cuts.AddCutCalo("1258da13","411798305k032220000","01631031000000d0"); //
    cuts.AddCutCalo("1468da13","411798305k032220000","01631031000000d0"); //
    cuts.AddCutCalo("1688da13","411798305k032220000","01631031000000d0"); //
  } else if (trainConfig == 713){ // EMCAL + DCal clusters - 20190301
    cuts.AddCutCalo("1018da13","411798305k0a2220000","01631061000000d0"); //
    cuts.AddCutCalo("1128da13","411798305k0b2220000","01631061000000d0"); //
  //EJ1 + DJ1
  } else if (trainConfig == 714){ // EMCAL + DCal clusters - 20190301
    cuts.AddCutCalo("1249ba13","411798305k032220000","01631031000000d0"); //
    cuts.AddCutCalo("1259ba13","411798305k032220000","01631031000000d0"); //
    cuts.AddCutCalo("1469ba13","411798305k032220000","01631031000000d0"); //
    cuts.AddCutCalo("1689ba13","411798305k032220000","01631031000000d0"); //
  } else if (trainConfig == 715){ // EMCAL + DCal clusters - 20190301
    cuts.AddCutCalo("1019ba13","411798305k0a2220000","01631061000000d0"); //
    cuts.AddCutCalo("1129ba13","411798305k0b2220000","01631061000000d0"); //
  //EJ2 + DJ2
  } else if (trainConfig == 716){ // EMCAL + DCal clusters - 20190301
    cuts.AddCutCalo("1249ca13","411798305k032220000","01631031000000d0"); //
    cuts.AddCutCalo("1259ca13","411798305k032220000","01631031000000d0"); //
    cuts.AddCutCalo("1469ca13","411798305k032220000","01631031000000d0"); //
    cuts.AddCutCalo("1689ca13","411798305k032220000","01631031000000d0"); //
  } else if (trainConfig == 717){ // EMCAL + DCal clusters - 20190301
    cuts.AddCutCalo("1019ca13","411798305k0a2220000","01631061000000d0"); //
    cuts.AddCutCalo("1129ca13","411798305k0b2220000","01631061000000d0"); //
  //EMC7 + DMC7
  } else if (trainConfig == 718){ // EMCAL + DCal clusters - 20190301
    cuts.AddCutCalo("12457a13","411798305k032220000","01631031000000d0"); //
    cuts.AddCutCalo("12557a13","411798305k032220000","01631031000000d0"); //
    cuts.AddCutCalo("14657a13","411798305k032220000","01631031000000d0"); //
    cuts.AddCutCalo("16857a13","411798305k032220000","01631031000000d0"); //
  } else if (trainConfig == 719){ // EMCAL + DCal clusters - 20190301
    cuts.AddCutCalo("10157a13","411798305k0a2220000","01631061000000d0"); //
    cuts.AddCutCalo("11257a13","411798305k0b2220000","01631061000000d0"); //

  // **********************************************************************************************************
  // ***************************** EMCal+DCal configurations PbPb run 2 2018 *******************************
  // **********************************************************************************************************
  } else if (trainConfig == 750){ // EMCAL+DCal clusters
    cuts.AddCutCalo("10910013","411790105ke30220000","01331031000000d0"); //  0-90%
  } else if (trainConfig == 751){ // EMCAL+DCal clusters - cent
    cuts.AddCutCalo("10110013","4117901050e30220000","0s631031000000d0"); // 00-10%
    cuts.AddCutCalo("30110013","4117901050e30220000","0s631031000000d0"); // 00-05%
    cuts.AddCutCalo("31210013","4117901050e30220000","0s631031000000d0"); // 05-10%
  } else if (trainConfig == 752){ // EMCAL+DCal clusters - semi-central
    cuts.AddCutCalo("11210013","4117901050e30220000","0s631031000000d0"); // 10-20%
    cuts.AddCutCalo("12310013","4117901050e30220000","0s631031000000d0"); // 20-30%
    cuts.AddCutCalo("13410013","4117901050e30220000","0s631031000000d0"); // 30-40%
    cuts.AddCutCalo("12410013","4117901050e30220000","0s631031000000d0"); // 20-40%
  } else if (trainConfig == 753){ // EMCAL+DCal clusters - semi peripheral
    cuts.AddCutCalo("14510013","4117901050e30220000","0s631031000000d0"); // 40-50%
    cuts.AddCutCalo("14610013","4117901050e30220000","0s631031000000d0"); // 40-60%
    cuts.AddCutCalo("15610013","4117901050e30220000","0s631031000000d0"); // 50-60%
  } else if (trainConfig == 754){ // EMCAL+DCal clusters - peripheral
    cuts.AddCutCalo("16710013","4117901050e30220000","0s631031000000d0"); // 60-70%
    cuts.AddCutCalo("17810013","4117901050e30220000","0s631031000000d0"); // 70-80%
    cuts.AddCutCalo("18910013","4117901050e30220000","0s631031000000d0"); // 80-90%
    cuts.AddCutCalo("16810013","4117901050e30220000","0s631031000000d0"); // 60-80%
  } else if (trainConfig == 755){ // EMCAL+DCal clusters - cent
    cuts.AddCutCalo("10110013","411790105fe30220000","0s631031000000d0"); // 00-10%
    cuts.AddCutCalo("30110013","411790105fe30220000","0s631031000000d0"); // 00-05%
    cuts.AddCutCalo("31210013","411790105fe30220000","0s631031000000d0"); // 05-10%
  } else if (trainConfig == 756){ // EMCAL+DCal clusters - semi-central
    cuts.AddCutCalo("11210013","411790105fe30220000","0s631031000000d0"); // 10-20%
    cuts.AddCutCalo("12310013","411790105fe30220000","0s631031000000d0"); // 20-30%
    cuts.AddCutCalo("13410013","411790105fe30220000","0s631031000000d0"); // 30-40%
    cuts.AddCutCalo("12410013","411790105fe30220000","0s631031000000d0"); // 20-40%
  } else if (trainConfig == 757){ // EMCAL+DCal clusters - semi peripheral
    cuts.AddCutCalo("14510013","411790105fe30220000","0s631031000000d0"); // 40-50%
    cuts.AddCutCalo("14610013","411790105fe30220000","0s631031000000d0"); // 40-60%
    cuts.AddCutCalo("15610013","411790105fe30220000","0s631031000000d0"); // 50-60%
  } else if (trainConfig == 758){ // EMCAL+DCal clusters - peripheral
    cuts.AddCutCalo("16710013","411790105fe30220000","0s631031000000d0"); // 60-70%
    cuts.AddCutCalo("17810013","411790105fe30220000","0s631031000000d0"); // 70-80%
    cuts.AddCutCalo("18910013","411790105fe30220000","0s631031000000d0"); // 80-90%
    cuts.AddCutCalo("16810013","411790105fe30220000","0s631031000000d0"); // 60-80%
  } else if (trainConfig == 759){ // EMCAL+DCal clusters - MIP sub
    cuts.AddCutCalo("10130a03","411790105k0a2220000","01331061000000d0"); //
    cuts.AddCutCalo("11310a03","411790105k0b2220000","01331061000000d0"); //
    cuts.AddCutCalo("13530a03","411790105k032220000","01331031000000d0"); //
    cuts.AddCutCalo("15910a03","411790105k032220000","01331031000000d0"); //

  } else if (trainConfig == 760){ // EMCAL clusters - centrality selection for PbPb EMCal
    cuts.AddCutCalo("10130a13","4117900050032220000","01331031000000d0"); //
    cuts.AddCutCalo("13530a13","4117900050032220000","01331031000000d0"); //
  } else if (trainConfig == 761){ // EMCAL clusters - centrality selection for PbPb EMCal
    cuts.AddCutCalo("10110a13","4117900050032220000","01331031000000d0"); //
    cuts.AddCutCalo("13510a13","4117900050032220000","01331031000000d0"); //
  } else if (trainConfig == 762){ // EG1-CENT- EMCAL clusters - centrality selection for PbPb EMCal
    cuts.AddCutCalo("1018da13","411790105fga2220000","01431031000000d0"); //
    cuts.AddCutCalo("1138da13","411790105fgb2220000","01431031000000d0"); //
    cuts.AddCutCalo("1358da13","411790105fg32220000","01431031000000d0"); //
    cuts.AddCutCalo("1598da13","411790105fg32220000","01431031000000d0"); //
  } else if (trainConfig == 763){ // EG2-CENT- EMCAL clusters - centrality selection for PbPb EMCal
    cuts.AddCutCalo("1018ea13","411790105fga2220000","01431031000000d0"); //
    cuts.AddCutCalo("1138ea13","411790105fgb2220000","01431031000000d0"); //
    cuts.AddCutCalo("1358ea13","411790105fg32220000","01431031000000d0"); //
    cuts.AddCutCalo("1598ea13","411790105fg32220000","01431031000000d0"); //
  } else if (trainConfig == 764){ // EG1-CALO- EMCAL clusters - centrality selection for PbPb EMCal
    cuts.AddCutCalo("101ala13","411790105fga2220000","01431031000000d0"); //
    cuts.AddCutCalo("113ala13","411790105fgb2220000","01431031000000d0"); //
    cuts.AddCutCalo("135ala13","411790105fg32220000","01431031000000d0"); //
    cuts.AddCutCalo("159ala13","411790105fg32220000","01431031000000d0"); //
  } else if (trainConfig == 765){ // EG2-CALO- EMCAL clusters - centrality selection for PbPb EMCal
    cuts.AddCutCalo("101ama13","411790105fga2220000","01431031000000d0"); //
    cuts.AddCutCalo("113ama13","411790105fgb2220000","01431031000000d0"); //
    cuts.AddCutCalo("135ama13","411790105fg32220000","01431031000000d0"); //
    cuts.AddCutCalo("159ama13","411790105fg32220000","01431031000000d0"); //

  } else if (trainConfig == 770){ // EMCAL+DCal clusters - QA
    cuts.AddCutCalo("10930013","4117900000032220000","01331031000000d0"); //  0-90% QA (open timing)

  } else if (trainConfig == 780){ // EMCAL+DCal standard with 4 cents
    cuts.AddCutCalo("10130e03","411790105ke30220000","01331031000000d0"); // 00-10%
    cuts.AddCutCalo("11310e03","411790105ke30220000","01331031000000d0"); // 10-30%
    cuts.AddCutCalo("13530e03","411790105ke30220000","01331031000000d0"); // 30-50%
    cuts.AddCutCalo("15910e03","411790105ke30220000","01331031000000d0"); // 50-90%
  } else if (trainConfig == 781){ // EMCAL+DCal standard with rotation method
    cuts.AddCutCalo("10130e03","411790105ke30220000","0s631031000000d0"); // 00-10%
    cuts.AddCutCalo("11310e03","411790105ke30220000","0s631031000000d0"); // 10-30%
    cuts.AddCutCalo("13530e03","411790105ke30220000","0s631031000000d0"); // 30-50%
    cuts.AddCutCalo("15910e03","411790105ke30220000","0s631031000000d0"); // 50-90%
  } else if (trainConfig == 782){ // EMCAL+DCal standard with 4 cents without OOB Pileup correction
    cuts.AddCutCalo("10130053","411790105ke30220000","01331031000000d0"); // 00-10%
    cuts.AddCutCalo("11310053","411790105ke30220000","01331031000000d0"); // 10-30%
    cuts.AddCutCalo("13530053","411790105ke30220000","01331031000000d0"); // 30-50%
    cuts.AddCutCalo("15910053","411790105ke30220000","01331031000000d0"); // 50-90%
  } else if (trainConfig == 783){ // EMCAL+DCal standard with rotation method without OOB Pileup correction
    cuts.AddCutCalo("10130053","411790105ke30220000","0s631031000000d0"); // 00-10%
    cuts.AddCutCalo("11310053","411790105ke30220000","0s631031000000d0"); // 10-30%
    cuts.AddCutCalo("13530053","411790105ke30220000","0s631031000000d0"); // 30-50%
    cuts.AddCutCalo("15910053","411790105ke30220000","0s631031000000d0"); // 50-90%
  } else if (trainConfig == 784){ // EMCAL+DCal clusters - Neutral Overlap with mean N matched tracks per cluster
    cuts.AddCutCalo("10130e03","411790105te30220000","01331031000000d0"); // 00-10%
    cuts.AddCutCalo("11310e03","411790105te30220000","01331031000000d0"); // 10-30%
    cuts.AddCutCalo("13530e03","411790105te30220000","01331031000000d0"); // 30-50%
    cuts.AddCutCalo("15910e03","411790105te30220000","01331031000000d0"); // 50-90%
  } else if (trainConfig == 785){ // EMCAL+DCal clusters - Neutral Overlap with mean N matched tracks per cluster without OOB Pileup correction
    cuts.AddCutCalo("10130053","411790105te30220000","01331031000000d0"); // 00-10%
    cuts.AddCutCalo("11310053","411790105te30220000","01331031000000d0"); // 10-30%
    cuts.AddCutCalo("13530053","411790105te30220000","01331031000000d0"); // 30-50%
    cuts.AddCutCalo("15910053","411790105te30220000","01331031000000d0"); // 50-90%
  } else if (trainConfig == 786){ // EMCAL+DCal clusters - Neutral Overlap NonLin like
    cuts.AddCutCalo("10130e03","411790105we30220000","01331031000000d0"); // 00-10%
    cuts.AddCutCalo("11310e03","411790105we30220000","01331031000000d0"); // 10-30%
    cuts.AddCutCalo("13530e03","411790105we30220000","01331031000000d0"); // 30-50%
    cuts.AddCutCalo("15910e03","411790105we30220000","01331031000000d0"); // 50-90%
  } else if (trainConfig == 787){ // EMCAL+DCal clusters -  Neutral Overlap NonLin like without OOB Pileup correction
    cuts.AddCutCalo("10130053","411790105we30220000","01331031000000d0"); // 00-10%
    cuts.AddCutCalo("11310053","411790105we30220000","01331031000000d0"); // 10-30%
    cuts.AddCutCalo("13530053","411790105we30220000","01331031000000d0"); // 30-50%
    cuts.AddCutCalo("15910053","411790105we30220000","01331031000000d0"); // 50-90%
  } else if (trainConfig == 788){ // EMCAL+DCal clusters - Neutral Overlap with mean charged tracks per cell
    cuts.AddCutCalo("10130e03","411790105ve30220000","01331031000000d0"); // 00-10%
    cuts.AddCutCalo("11310e03","411790105ve30220000","01331031000000d0"); // 10-30%
    cuts.AddCutCalo("13530e03","411790105ve30220000","01331031000000d0"); // 30-50%
    cuts.AddCutCalo("15910e03","411790105ve30220000","01331031000000d0"); // 50-90%
  } else if (trainConfig == 789){ // EMCAL+DCal clusters -  Neutral Overlap with mean charged tracks per cell without OOB Pileup correction
    cuts.AddCutCalo("10130053","411790105ve30220000","01331031000000d0"); // 00-10%
    cuts.AddCutCalo("11310053","411790105ve30220000","01331031000000d0"); // 10-30%
    cuts.AddCutCalo("13530053","411790105ve30220000","01331031000000d0"); // 30-50%
    cuts.AddCutCalo("15910053","411790105ve30220000","01331031000000d0"); // 50-90%
  // Random neutral energy overlap correction!
  } else if (trainConfig == 790){ // EMCAL+DCal clusters - 13 TeV emcal cuts w/o TM
    cuts.AddCutCalo("10130e03","411790105ue30220000","0s631031000000d0"); // 00-10%
    cuts.AddCutCalo("11310e03","411790105ue30220000","0s631031000000d0"); // 10-30%
    cuts.AddCutCalo("13530e03","411790105ue30220000","0s631031000000d0"); // 30-50%
    cuts.AddCutCalo("15910e03","411790105ue30220000","0s631031000000d0"); // 50-90%
  } else if (trainConfig == 791){ // EMCAL+DCal clusters - 13 TeV emcal cuts w/o TM with Mean neutral energy overlap correction, cent and semi cent trigger with OOB pileup cut
    cuts.AddCutCalo("10130e13","411790105te30220000","0s631031000000d0"); // 00-10%
    cuts.AddCutCalo("11310e13","411790105te30220000","0s631031000000d0"); // 10-30%
    cuts.AddCutCalo("13530e13","411790105te30220000","0s631031000000d0"); // 30-50%
    cuts.AddCutCalo("15910e13","411790105te30220000","0s631031000000d0"); // 50-90%
  // Mean neutral energy overlap correction!
  } else if (trainConfig == 792){ // EMCAL+DCal clusters - 13 TeV emcal cuts w/o TM with Mean neutral energy overlap correction, cent and semi cent trigger without OOB pileup cut
    cuts.AddCutCalo("10130013","411790105te30220000","0s631031000000d0"); // 00-10%
    cuts.AddCutCalo("11310013","411790105te30220000","0s631031000000d0"); // 10-30%
    cuts.AddCutCalo("13530013","411790105te30220000","0s631031000000d0"); // 30-50%
    cuts.AddCutCalo("15910013","411790105te30220000","0s631031000000d0"); // 50-90%
  } else if (trainConfig == 793){ // EMCAL+DCal clusters - 13 TeV emcal cuts w/o TM with Mean neutral energy overlap correction using Ncharge per cell, cent and semi cent trigger with OOB pileup cut
    cuts.AddCutCalo("10130e13","411790105ve30220000","0s631031000000d0"); // 00-10%
    cuts.AddCutCalo("11310e13","411790105ve30220000","0s631031000000d0"); // 10-30%
    cuts.AddCutCalo("13530e13","411790105ve30220000","0s631031000000d0"); // 30-50%
    cuts.AddCutCalo("15910e13","411790105ve30220000","0s631031000000d0"); // 50-90%
  } else if (trainConfig == 794){ // EMCAL+DCal standard with 4 cents without OOB Pileup correction, Added Signal
    cuts.AddCutCalo("10130023","411790105ke30220000","01331031000000d0"); // 00-10%
    cuts.AddCutCalo("11310023","411790105ke30220000","01331031000000d0"); // 10-30%
    cuts.AddCutCalo("13530023","411790105ke30220000","01331031000000d0"); // 30-50%
    cuts.AddCutCalo("15910023","411790105ke30220000","01331031000000d0"); // 50-90%
  } else if (trainConfig == 795){ // EMCAL+DCal clusters -  Neutral Overlap NonLin like without OOB Pileup correction, Added Signal
    cuts.AddCutCalo("10130023","411790105we30220000","01331031000000d0"); // 00-10%
    cuts.AddCutCalo("11310023","411790105we30220000","01331031000000d0"); // 10-30%
    cuts.AddCutCalo("13530023","411790105we30220000","01331031000000d0"); // 30-50%
    cuts.AddCutCalo("15910023","411790105we30220000","01331031000000d0"); // 50-90%
  } else if (trainConfig == 796){ // EMCAL+DCal clusters -  Neutral Overlap with mean charged tracks per cell without OOB Pileup correction, Added Signal
    cuts.AddCutCalo("10130023","411790105ve30220000","01331031000000d0"); // 00-10%
    cuts.AddCutCalo("11310023","411790105ve30220000","01331031000000d0"); // 10-30%
    cuts.AddCutCalo("13530023","411790105ve30220000","01331031000000d0"); // 30-50%
    cuts.AddCutCalo("15910023","411790105ve30220000","01331031000000d0"); // 50-90%
  } else if (trainConfig == 797){ // EMCAL+DCal clusters - Neutral Overlap NonLin like - OAC 0.0152
    cuts.AddCutCalo("10130e03","411790105we30220000","01331031000000b0"); // 00-10%
    cuts.AddCutCalo("11310e03","411790105we30220000","01331031000000b0"); // 10-30%
    cuts.AddCutCalo("13530e03","411790105we30220000","01331031000000b0"); // 30-50%
    cuts.AddCutCalo("15910e03","411790105we30220000","01331031000000b0"); // 50-90%
  } else if (trainConfig == 798){ // EMCAL+DCal clusters -  Neutral Overlap NonLin like without OOB Pileup correction - OAC 0.0152
    cuts.AddCutCalo("10130053","411790105we30220000","01331031000000b0"); // 00-10%
    cuts.AddCutCalo("11310053","411790105we30220000","01331031000000b0"); // 10-30%
    cuts.AddCutCalo("13530053","411790105we30220000","01331031000000b0"); // 30-50%
    cuts.AddCutCalo("15910053","411790105we30220000","01331031000000b0"); // 50-90%
  } else if (trainConfig == 799){ // EMCAL+DCal clusters -  Neutral Overlap NonLin like without OOB Pileup correction, Added Signal - OAC 0.0152
    cuts.AddCutCalo("10130023","411790105we30220000","01331031000000b0"); // 00-10%
    cuts.AddCutCalo("11310023","411790105we30220000","01331031000000b0"); // 10-30%
    cuts.AddCutCalo("13530023","411790105we30220000","01331031000000b0"); // 30-50%
    cuts.AddCutCalo("15910023","411790105we30220000","01331031000000b0"); // 50-90%
  } else if (trainConfig == 800){ // EMCAL+DCal clusters - Neutral Overlap NonLin like - OAC 0.0202
    cuts.AddCutCalo("10130e03","411790105we30220000","01331031000000g0"); // 00-10%
    cuts.AddCutCalo("11310e03","411790105we30220000","01331031000000g0"); // 10-30%
    cuts.AddCutCalo("13530e03","411790105we30220000","01331031000000g0"); // 30-50%
    cuts.AddCutCalo("15910e03","411790105we30220000","01331031000000g0"); // 50-90%
  } else if (trainConfig == 801){ // EMCAL+DCal clusters -  Neutral Overlap NonLin like without OOB Pileup correction - OAC 0.0202
    cuts.AddCutCalo("10130053","411790105we30220000","01331031000000g0"); // 00-10%
    cuts.AddCutCalo("11310053","411790105we30220000","01331031000000g0"); // 10-30%
    cuts.AddCutCalo("13530053","411790105we30220000","01331031000000g0"); // 30-50%
    cuts.AddCutCalo("15910053","411790105we30220000","01331031000000g0"); // 50-90%
  } else if (trainConfig == 802){ // EMCAL+DCal clusters -  Neutral Overlap NonLin like without OOB Pileup correction, Added Signal - OAC 0.0202
    cuts.AddCutCalo("10130023","411790105we30220000","01331031000000g0"); // 00-10%
    cuts.AddCutCalo("11310023","411790105we30220000","01331031000000g0"); // 10-30%
    cuts.AddCutCalo("13530023","411790105we30220000","01331031000000g0"); // 30-50%
    cuts.AddCutCalo("15910023","411790105we30220000","01331031000000g0"); // 50-90%
  } else if (trainConfig == 803){ // EMCAL+DCal clusters - Neutral Overlap NonLin like - rotation background method
    cuts.AddCutCalo("10130e03","411790105we30220000","0s331031000000d0"); // 00-10%
    cuts.AddCutCalo("11310e03","411790105we30220000","0s331031000000d0"); // 10-30%
    cuts.AddCutCalo("13530e03","411790105we30220000","0s331031000000d0"); // 30-50%
    cuts.AddCutCalo("15910e03","411790105we30220000","0s331031000000d0"); // 50-90%
  } else if (trainConfig == 804){ // EMCAL+DCal clusters -  Neutral Overlap NonLin like without OOB Pileup correction - rotation background method
    cuts.AddCutCalo("10130053","411790105we30220000","0s331031000000d0"); // 00-10%
    cuts.AddCutCalo("11310053","411790105we30220000","0s331031000000d0"); // 10-30%
    cuts.AddCutCalo("13530053","411790105we30220000","0s331031000000d0"); // 30-50%
    cuts.AddCutCalo("15910053","411790105we30220000","0s331031000000d0"); // 50-90%
  } else if (trainConfig == 805){ // EMCAL+DCal clusters -  Neutral Overlap NonLin like without OOB Pileup correction, Added Signal - rotation background method
    cuts.AddCutCalo("10130023","411790105we30220000","0s331031000000d0"); // 00-10%
    cuts.AddCutCalo("11310023","411790105we30220000","0s331031000000d0"); // 10-30%
    cuts.AddCutCalo("13530023","411790105we30220000","0s331031000000d0"); // 30-50%
    cuts.AddCutCalo("15910023","411790105we30220000","0s331031000000d0"); // 50-90%
  } else if (trainConfig == 806){ // EMCAL+DCal clusters - Neutral Overlap Correction maxwell boltzmann
    cuts.AddCutCalo("10130e03","411790105xe30220000","01331031000000d0"); // 00-10%
    cuts.AddCutCalo("11310e03","411790105xe30220000","01331031000000d0"); // 10-30%
    cuts.AddCutCalo("13530e03","411790105xe30220000","01331031000000d0"); // 30-50%
    cuts.AddCutCalo("15910e03","411790105xe30220000","01331031000000d0"); // 50-90%
  } else if (trainConfig == 807){ // EMCAL+DCal clusters -  Neutral Overlap Correction maxwell boltzmann without OOB Pileup correction
    cuts.AddCutCalo("10130053","411790105xe30220000","01331031000000d0"); // 00-10%
    cuts.AddCutCalo("11310053","411790105xe30220000","01331031000000d0"); // 10-30%
    cuts.AddCutCalo("13530053","411790105xe30220000","01331031000000d0"); // 30-50%
    cuts.AddCutCalo("15910053","411790105xe30220000","01331031000000d0"); // 50-90%
  } else if (trainConfig == 808){ // EMCAL+DCal clusters -  Neutral Overlap Correction maxwell boltzmann without OOB Pileup correction, Added Signal
    cuts.AddCutCalo("10130023","411790105xe30220000","01331031000000d0"); // 00-10%
    cuts.AddCutCalo("11310023","411790105xe30220000","01331031000000d0"); // 10-30%
    cuts.AddCutCalo("13530023","411790105xe30220000","01331031000000d0"); // 30-50%
    cuts.AddCutCalo("15910023","411790105xe30220000","01331031000000d0"); // 50-90%
  // **********************************************************************************************************
  // ***************************** PHOS       QA configurations PbPb run 2 2018 *******************************
  // **********************************************************************************************************
  } else if (trainConfig == 850){ // PHOS clusters - 20181018
    cuts.AddCutCalo("10930013","24466000h0082300000","0133103100000010"); // 0-90%
  } else if (trainConfig == 851){ // PHOS clusters - 20181018
    cuts.AddCutCalo("10130013","24466000h0082300000","0133103100000010"); //
    cuts.AddCutCalo("11530013","24466000h0082300000","0133103100000010"); //
    cuts.AddCutCalo("15930013","24466000h0082300000","0133103100000010"); //
  } else if (trainConfig == 852){ // PHOS clusters - 20181018
    cuts.AddCutCalo("10130a13","24466000ha082200000","0133103100000010"); //
    cuts.AddCutCalo("11230a13","24466000ha082200000","0133103100000010"); //
    cuts.AddCutCalo("12430a13","24466000ha082200000","0133103100000010"); //
    cuts.AddCutCalo("14630a13","24466000ha082200000","0133103100000010"); //
    cuts.AddCutCalo("16830a13","24466000ha082200000","0133103100000010"); //
  } else if (trainConfig == 853){ // PHOS clusters - 20181018
    cuts.AddCutCalo("10130a13","24466000ha082200000","0133103100000010"); //
    cuts.AddCutCalo("11310a13","24466000ha082200000","0133103100000010"); //
    cuts.AddCutCalo("13530a13","24466000ha082200000","0133103100000010"); //
    cuts.AddCutCalo("15910a13","24466000ha082200000","0133103100000010"); //
  } else if (trainConfig == 854){ // PHOS clusters - 20181018
    cuts.AddCutCalo("10130a13","24466810ha082200000","0133103100000010"); //
    cuts.AddCutCalo("11310a13","24466810ha082200000","0133103100000010"); //
    cuts.AddCutCalo("13530a13","24466810ha082200000","0133103100000010"); //
    cuts.AddCutCalo("15910a13","24466810ha082200000","0133103100000010"); //
  } else if (trainConfig == 855){ // PHOS clusters - 20181018
    cuts.AddCutCalo("30130a13","24466810ha082200000","0143103100000010"); //
    cuts.AddCutCalo("31230a13","24466810ha082200000","0143103100000010"); //
    cuts.AddCutCalo("11210a13","24466810ha082200000","0143103100000010"); //
    cuts.AddCutCalo("12310a13","24466810ha082200000","0143103100000010"); //
  } else if (trainConfig == 856){ // PHOS clusters - 20181018
    cuts.AddCutCalo("13430a13","24466810ha082200000","0143103100000010"); //
    cuts.AddCutCalo("14530a13","24466810ha082200000","0143103100000010"); //
    cuts.AddCutCalo("15610a13","24466810ha082200000","0143103100000010"); //
    cuts.AddCutCalo("16710a13","24466810ha082200000","0143103100000010"); //
    cuts.AddCutCalo("17810a13","24466810ha082200000","0143103100000010"); //
    cuts.AddCutCalo("18910a13","24466810ha082200000","0143103100000010"); //
  } else if (trainConfig == 860){ // PHOS clusters - trigger 20181018
    cuts.AddCutCalo("10162a13","24466810ha082200000","0143103100000010"); //
    cuts.AddCutCalo("11362a13","24466810ha082200000","0143103100000010"); //
    cuts.AddCutCalo("13562a13","24466810ha082200000","0143103100000010"); //
    cuts.AddCutCalo("15962a13","24466810ha082200000","0143103100000010"); //
  } else if (trainConfig == 880){ // PHOS clusters - similar to EMCal 780
    cuts.AddCutCalo("10130e03","24466810ha082200000","0133103100000010"); // 00-10%
    cuts.AddCutCalo("11310e03","24466810ha082200000","0133103100000010"); // 10-30%
    cuts.AddCutCalo("13530e03","24466810ha082200000","0133103100000010"); // 30-50%
    cuts.AddCutCalo("15910e03","24466810ha082200000","0133103100000010"); // 50-90%
  } else if (trainConfig == 881){ // PHOS clusters - similar to EMCal 780 only 0-10 trigger
    cuts.AddCutCalo("10130e03","24466810ha082200000","0133103100000010"); // 00-10%
  } else if (trainConfig == 882){ // PHOS clusters - similar to EMCal 780 only 30-50 trigger
    cuts.AddCutCalo("13530e03","24466810ha082200000","0133103100000010"); // 30-50%

  // **********************************************************************************************************
  // ****************************** EMCal+DCal configurations PbPb run 2 2018 *********************************
  // ********************************** only use clusters with matched track **********************************
  // **********************************************************************************************************
  } else if (trainConfig == 955){ // EMCAL+DCal clusters - cent
    cuts.AddCutCalo("10110013","411790105qe30220000","0s631031000000d0"); // 00-10%
    cuts.AddCutCalo("30110013","411790105qe30220000","0s631031000000d0"); // 00-05%
    cuts.AddCutCalo("31210013","411790105qe30220000","0s631031000000d0"); // 05-10%
  } else if (trainConfig == 956){ // EMCAL+DCal clusters - semi-central
    cuts.AddCutCalo("11210013","411790105qe30220000","0s631031000000d0"); // 10-20%
    cuts.AddCutCalo("12310013","411790105qe30220000","0s631031000000d0"); // 20-30%
    cuts.AddCutCalo("13410013","411790105qe30220000","0s631031000000d0"); // 30-40%
    cuts.AddCutCalo("12410013","411790105qe30220000","0s631031000000d0"); // 20-40%
  } else if (trainConfig == 957){ // EMCAL+DCal clusters - semi peripheral
    cuts.AddCutCalo("14510013","411790105qe30220000","0s631031000000d0"); // 40-50%
    cuts.AddCutCalo("14610013","411790105qe30220000","0s631031000000d0"); // 40-60%
    cuts.AddCutCalo("15610013","411790105qe30220000","0s631031000000d0"); // 50-60%
  } else if (trainConfig == 958){ // EMCAL+DCal clusters - peripheral
    cuts.AddCutCalo("16710013","411790105qe30220000","0s631031000000d0"); // 60-70%
    cuts.AddCutCalo("17810013","411790105qe30220000","0s631031000000d0"); // 70-80%
    cuts.AddCutCalo("18910013","411790105qe30220000","0s631031000000d0"); // 80-90%
    cuts.AddCutCalo("16810013","411790105qe30220000","0s631031000000d0"); // 60-80%

  // **********************************************************************************************************
  // ***************************** EMC configurations PbPb run 2 2018 pass 3 *******************************
  // **********************************************************************************************************
  // Standard: EMCal + DCal, only 2ndary cluster track matching, mixed event
  } else if (trainConfig == 1000){ // 4 cents
    cuts.AddCutCalo("10130e03","411790105ke30220000","01331031000000d0"); // 00-10%
    cuts.AddCutCalo("11310e03","411790105ke30220000","01331031000000d0"); // 10-30%
    cuts.AddCutCalo("13530e03","411790105ke30220000","01331031000000d0"); // 30-50%
    cuts.AddCutCalo("15910e03","411790105ke30220000","01331031000000d0"); // 50-90%
  } else if (trainConfig == 1001){ // 4 cents, no OOB Pileup correction for MC
    cuts.AddCutCalo("10130053","411790105ke30220000","01331031000000d0"); // 00-10%
    cuts.AddCutCalo("11310053","411790105ke30220000","01331031000000d0"); // 10-30%
    cuts.AddCutCalo("13530053","411790105ke30220000","01331031000000d0"); // 30-50%
    cuts.AddCutCalo("15910053","411790105ke30220000","01331031000000d0"); // 50-90%
  } else if (trainConfig == 1002){ // 4 cents, no OOB Pileup correction for MC with added signal
    cuts.AddCutCalo("10130023","411790105ke30220000","01331031000000d0"); // 00-10%
    cuts.AddCutCalo("11310023","411790105ke30220000","01331031000000d0"); // 10-30%
    cuts.AddCutCalo("13530023","411790105ke30220000","01331031000000d0"); // 30-50%
    cuts.AddCutCalo("15910023","411790105ke30220000","01331031000000d0"); // 50-90%
  } else if (trainConfig == 1003){ // 4 cents, no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutCalo("10130023","411790105ke30220000","01331031000000d0"); // 00-10%
    cuts.AddCutCalo("11310023","411790105ke30220000","01331031000000d0"); // 10-30%
    cuts.AddCutCalo("13530023","411790105ke30220000","01331031000000d0"); // 30-50%
    cuts.AddCutCalo("15910023","411790105ke30220000","01331031000000d0"); // 50-90%
  } else if (trainConfig == 1004){ // 4 cents, no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutCalo("10130023","411790105ke30220000","01331031000000d0"); // 00-10%
    cuts.AddCutCalo("11310023","411790105ke30220000","01331031000000d0"); // 10-30%
    cuts.AddCutCalo("13530023","411790105ke30220000","01331031000000d0"); // 30-50%
    cuts.AddCutCalo("15910023","411790105ke30220000","01331031000000d0"); // 50-90%
  } else if (trainConfig == 1005){ // 4 cents with NOC NonLin like data
    cuts.AddCutCalo("10130e03","411790105we30220000","01331031000000d0"); // 00-10%
    cuts.AddCutCalo("11310e03","411790105we30220000","01331031000000d0"); // 10-30%
    cuts.AddCutCalo("13530e03","411790105we30220000","01331031000000d0"); // 30-50%
    cuts.AddCutCalo("15910e03","411790105we30220000","01331031000000d0"); // 50-90%
  } else if (trainConfig == 1006){ // 4 cents, with NOC NonLin like MC no OOB Pileup correction for MC
    cuts.AddCutCalo("10130053","411790105xe30220000","01331031000000d0"); // 00-10%
    cuts.AddCutCalo("11310053","411790105xe30220000","01331031000000d0"); // 10-30%
    cuts.AddCutCalo("13530053","411790105xe30220000","01331031000000d0"); // 30-50%
    cuts.AddCutCalo("15910053","411790105xe30220000","01331031000000d0"); // 50-90%
  } else if (trainConfig == 1007){ // 4 cents, with NOC NonLin like MC no OOB Pileup correction for MC with added signal
    cuts.AddCutCalo("10130023","411790105xe30220000","01331031000000d0"); // 00-10%
    cuts.AddCutCalo("11310023","411790105xe30220000","01331031000000d0"); // 10-30%
    cuts.AddCutCalo("13530023","411790105xe30220000","01331031000000d0"); // 30-50%
    cuts.AddCutCalo("15910023","411790105xe30220000","01331031000000d0"); // 50-90%
  } else if (trainConfig == 1008){ // 4 cents, with NOC NonLin like MC no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutCalo("10130023","411790105xe30220000","01331031000000d0"); // 00-10%
    cuts.AddCutCalo("11310023","411790105xe30220000","01331031000000d0"); // 10-30%
    cuts.AddCutCalo("13530023","411790105xe30220000","01331031000000d0"); // 30-50%
    cuts.AddCutCalo("15910023","411790105xe30220000","01331031000000d0"); // 50-90%
  } else if (trainConfig == 1009){ // 4 cents, with NOC NonLin like MC no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutCalo("10130023","411790105xe30220000","01331031000000d0"); // 00-10%
    cuts.AddCutCalo("11310023","411790105xe30220000","01331031000000d0"); // 10-30%
    cuts.AddCutCalo("13530023","411790105xe30220000","01331031000000d0"); // 30-50%
    cuts.AddCutCalo("15910023","411790105xe30220000","01331031000000d0"); // 50-90%
  } else if (trainConfig == 1010){ // 4 cents with mean N matched tracks per cluster data
    cuts.AddCutCalo("10130e03","411790105te30220000","01331031000000d0"); // 00-10%
    cuts.AddCutCalo("11310e03","411790105te30220000","01331031000000d0"); // 10-30%
    cuts.AddCutCalo("13530e03","411790105te30220000","01331031000000d0"); // 30-50%
    cuts.AddCutCalo("15910e03","411790105te30220000","01331031000000d0"); // 50-90%
  } else if (trainConfig == 1011){ // 4 cents, with mean N matched tracks per cluster MC no OOB Pileup correction for MC
    cuts.AddCutCalo("10130053","411790105te30220000","01331031000000d0"); // 00-10%
    cuts.AddCutCalo("11310053","411790105te30220000","01331031000000d0"); // 10-30%
    cuts.AddCutCalo("13530053","411790105te30220000","01331031000000d0"); // 30-50%
    cuts.AddCutCalo("15910053","411790105te30220000","01331031000000d0"); // 50-90%
  } else if (trainConfig == 1012){ // 4 cents, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal
    cuts.AddCutCalo("10130023","411790105te30220000","01331031000000d0"); // 00-10%
    cuts.AddCutCalo("11310023","411790105te30220000","01331031000000d0"); // 10-30%
    cuts.AddCutCalo("13530023","411790105te30220000","01331031000000d0"); // 30-50%
    cuts.AddCutCalo("15910023","411790105te30220000","01331031000000d0"); // 50-90%
  } else if (trainConfig == 1013){ // 4 cents, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutCalo("10130023","411790105te30220000","01331031000000d0"); // 00-10%
    cuts.AddCutCalo("11310023","411790105te30220000","01331031000000d0"); // 10-30%
    cuts.AddCutCalo("13530023","411790105te30220000","01331031000000d0"); // 30-50%
    cuts.AddCutCalo("15910023","411790105te30220000","01331031000000d0"); // 50-90%
  } else if (trainConfig == 1014){ // 4 cents, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutCalo("10130023","411790105te30220000","01331031000000d0"); // 00-10%
    cuts.AddCutCalo("11310023","411790105te30220000","01331031000000d0"); // 10-30%
    cuts.AddCutCalo("13530023","411790105te30220000","01331031000000d0"); // 30-50%
    cuts.AddCutCalo("15910023","411790105te30220000","01331031000000d0"); // 50-90%
  } else if (trainConfig == 1015){ // 4 cents with mean N matched tracks per cluster data
    cuts.AddCutCalo("10130e03","411790105ye30220000","01331031000000d0"); // 00-10%
    cuts.AddCutCalo("11310e03","411790105ye30220000","01331031000000d0"); // 10-30%
    cuts.AddCutCalo("13530e03","411790105ye30220000","01331031000000d0"); // 30-50%
    cuts.AddCutCalo("15910e03","411790105ye30220000","01331031000000d0"); // 50-90%
  } else if (trainConfig == 1016){ // 4 cents, with mean N matched tracks per cluster MC no OOB Pileup correction for MC
    cuts.AddCutCalo("10130053","411790105ye30220000","01331031000000d0"); // 00-10%
    cuts.AddCutCalo("11310053","411790105ye30220000","01331031000000d0"); // 10-30%
    cuts.AddCutCalo("13530053","411790105ye30220000","01331031000000d0"); // 30-50%
    cuts.AddCutCalo("15910053","411790105ye30220000","01331031000000d0"); // 50-90%
  // with rotation background
  } else if (trainConfig == 1050){ // 4 cents
    cuts.AddCutCalo("10130e03","411790105ke30220000","0s331031000000d0"); // 00-10%
    cuts.AddCutCalo("11310e03","411790105ke30220000","0s331031000000d0"); // 10-30%
    cuts.AddCutCalo("13530e03","411790105ke30220000","0s331031000000d0"); // 30-50%
    cuts.AddCutCalo("15910e03","411790105ke30220000","0s331031000000d0"); // 50-90%
  } else if (trainConfig == 1051){ // 4 cents, no OOB Pileup correction for MC
    cuts.AddCutCalo("10130053","411790105ke30220000","0s331031000000d0"); // 00-10%
    cuts.AddCutCalo("11310053","411790105ke30220000","0s331031000000d0"); // 10-30%
    cuts.AddCutCalo("13530053","411790105ke30220000","0s331031000000d0"); // 30-50%
    cuts.AddCutCalo("15910053","411790105ke30220000","0s331031000000d0"); // 50-90%
  } else if (trainConfig == 1052){ // 4 cents, no OOB Pileup correction for MC with added signal
    cuts.AddCutCalo("10130023","411790105ke30220000","0s331031000000d0"); // 00-10%
    cuts.AddCutCalo("11310023","411790105ke30220000","0s331031000000d0"); // 10-30%
    cuts.AddCutCalo("13530023","411790105ke30220000","0s331031000000d0"); // 30-50%
    cuts.AddCutCalo("15910023","411790105ke30220000","0s331031000000d0"); // 50-90%
  } else if (trainConfig == 1053){ // 4 cents, no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutCalo("10130023","411790105ke30220000","0s331031000000d0"); // 00-10%
    cuts.AddCutCalo("11310023","411790105ke30220000","0s331031000000d0"); // 10-30%
    cuts.AddCutCalo("13530023","411790105ke30220000","0s331031000000d0"); // 30-50%
    cuts.AddCutCalo("15910023","411790105ke30220000","0s331031000000d0"); // 50-90%
  } else if (trainConfig == 1054){ // 4 cents, no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutCalo("10130023","411790105ke30220000","0s331031000000d0"); // 00-10%
    cuts.AddCutCalo("11310023","411790105ke30220000","0s331031000000d0"); // 10-30%
    cuts.AddCutCalo("13530023","411790105ke30220000","0s331031000000d0"); // 30-50%
    cuts.AddCutCalo("15910023","411790105ke30220000","0s331031000000d0"); // 50-90%
  } else if (trainConfig == 1055){ // 4 cents with NOC NonLin like data
    cuts.AddCutCalo("10130e03","411790105we30220000","0s331031000000d0"); // 00-10%
    cuts.AddCutCalo("11310e03","411790105we30220000","0s331031000000d0"); // 10-30%
    cuts.AddCutCalo("13530e03","411790105we30220000","0s331031000000d0"); // 30-50%
    cuts.AddCutCalo("15910e03","411790105we30220000","0s331031000000d0"); // 50-90%
  } else if (trainConfig == 1056){ // 4 cents, with NOC NonLin like MC no OOB Pileup correction for MC
    cuts.AddCutCalo("10130053","411790105xe30220000","0s331031000000d0"); // 00-10%
    cuts.AddCutCalo("11310053","411790105xe30220000","0s331031000000d0"); // 10-30%
    cuts.AddCutCalo("13530053","411790105xe30220000","0s331031000000d0"); // 30-50%
    cuts.AddCutCalo("15910053","411790105xe30220000","0s331031000000d0"); // 50-90%
  } else if (trainConfig == 1057){ // 4 cents, with NOC NonLin like MC no OOB Pileup correction for MC with added signal
    cuts.AddCutCalo("10130023","411790105xe30220000","0s331031000000d0"); // 00-10%
    cuts.AddCutCalo("11310023","411790105xe30220000","0s331031000000d0"); // 10-30%
    cuts.AddCutCalo("13530023","411790105xe30220000","0s331031000000d0"); // 30-50%
    cuts.AddCutCalo("15910023","411790105xe30220000","0s331031000000d0"); // 50-90%
  } else if (trainConfig == 1058){ // 4 cents, with NOC NonLin like MC no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutCalo("10130023","411790105xe30220000","0s331031000000d0"); // 00-10%
    cuts.AddCutCalo("11310023","411790105xe30220000","0s331031000000d0"); // 10-30%
    cuts.AddCutCalo("13530023","411790105xe30220000","0s331031000000d0"); // 30-50%
    cuts.AddCutCalo("15910023","411790105xe30220000","0s331031000000d0"); // 50-90%
  } else if (trainConfig == 1059){ // 4 cents, with NOC NonLin like MC no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutCalo("10130023","411790105xe30220000","0s331031000000d0"); // 00-10%
    cuts.AddCutCalo("11310023","411790105xe30220000","0s331031000000d0"); // 10-30%
    cuts.AddCutCalo("13530023","411790105xe30220000","0s331031000000d0"); // 30-50%
    cuts.AddCutCalo("15910023","411790105xe30220000","0s331031000000d0"); // 50-90%
  } else if (trainConfig == 1060){ // 4 cents with mean N matched tracks per cluster data
    cuts.AddCutCalo("10130e03","411790105te30220000","0s331031000000d0"); // 00-10%
    cuts.AddCutCalo("11310e03","411790105te30220000","0s331031000000d0"); // 10-30%
    cuts.AddCutCalo("13530e03","411790105te30220000","0s331031000000d0"); // 30-50%
    cuts.AddCutCalo("15910e03","411790105te30220000","0s331031000000d0"); // 50-90%
  } else if (trainConfig == 1061){ // 4 cents, with mean N matched tracks per cluster MC no OOB Pileup correction for MC
    cuts.AddCutCalo("10130053","411790105te30220000","0s331031000000d0"); // 00-10%
    cuts.AddCutCalo("11310053","411790105te30220000","0s331031000000d0"); // 10-30%
    cuts.AddCutCalo("13530053","411790105te30220000","0s331031000000d0"); // 30-50%
    cuts.AddCutCalo("15910053","411790105te30220000","0s331031000000d0"); // 50-90%
  } else if (trainConfig == 1062){ // 4 cents, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal
    cuts.AddCutCalo("10130023","411790105te30220000","0s331031000000d0"); // 00-10%
    cuts.AddCutCalo("11310023","411790105te30220000","0s331031000000d0"); // 10-30%
    cuts.AddCutCalo("13530023","411790105te30220000","0s331031000000d0"); // 30-50%
    cuts.AddCutCalo("15910023","411790105te30220000","0s331031000000d0"); // 50-90%
  } else if (trainConfig == 1063){ // 4 cents, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutCalo("10130023","411790105te30220000","0s331031000000d0"); // 00-10%
    cuts.AddCutCalo("11310023","411790105te30220000","0s331031000000d0"); // 10-30%
    cuts.AddCutCalo("13530023","411790105te30220000","0s331031000000d0"); // 30-50%
    cuts.AddCutCalo("15910023","411790105te30220000","0s331031000000d0"); // 50-90%
  } else if (trainConfig == 1064){ // 4 cents, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutCalo("10130023","411790105te30220000","0s331031000000d0"); // 00-10%
    cuts.AddCutCalo("11310023","411790105te30220000","0s331031000000d0"); // 10-30%
    cuts.AddCutCalo("13530023","411790105te30220000","0s331031000000d0"); // 30-50%
    cuts.AddCutCalo("15910023","411790105te30220000","0s331031000000d0"); // 50-90%
  } else if (trainConfig == 1065){ // 4 cents with mean N matched tracks per cluster data
    cuts.AddCutCalo("10130e03","411790105ye30220000","0s331031000000d0"); // 00-10%
    cuts.AddCutCalo("11310e03","411790105ye30220000","0s331031000000d0"); // 10-30%
    cuts.AddCutCalo("13530e03","411790105ye30220000","0s331031000000d0"); // 30-50%
    cuts.AddCutCalo("15910e03","411790105ye30220000","0s331031000000d0"); // 50-90%
  } else if (trainConfig == 1066){ // 4 cents, with mean N matched tracks per cluster MC no OOB Pileup correction for MC
    cuts.AddCutCalo("10130053","411790105ye30220000","0s331031000000d0"); // 00-10%
    cuts.AddCutCalo("11310053","411790105ye30220000","0s331031000000d0"); // 10-30%
    cuts.AddCutCalo("13530053","411790105ye30220000","0s331031000000d0"); // 30-50%
    cuts.AddCutCalo("15910053","411790105ye30220000","0s331031000000d0"); // 50-90%
  // **********************************************************************************************************
  // **************************** EMC configurations PbPb run 2 2018 pass 3 cent ******************************
  // **********************************************************************************************************
  // Standard: EMCal + DCal, only 2ndary cluster track matching, mixed event
  } else if (trainConfig == 1100){ // central
    cuts.AddCutCalo("10130e03","411790105ke30220000","01331031000000d0"); // 00-10%
  } else if (trainConfig == 1101){ // central, no OOB Pileup correction for MC
    cuts.AddCutCalo("10130053","411790105ke30220000","01331031000000d0"); // 00-10%
  } else if (trainConfig == 1102){ // central, no OOB Pileup correction for MC with added signal
    cuts.AddCutCalo("10130023","411790105ke30220000","01331031000000d0"); // 00-10%
  } else if (trainConfig == 1103){ // central, no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutCalo("10130023","411790105ke30220000","01331031000000d0"); // 00-10%
  } else if (trainConfig == 1104){ // central, no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutCalo("10130023","411790105ke30220000","01331031000000d0"); // 00-10%
  } else if (trainConfig == 1105){ // central with NOC NonLin like data
    cuts.AddCutCalo("10130e03","411790105we30220000","01331031000000d0"); // 00-10%
  } else if (trainConfig == 1106){ // central, with NOC NonLin like MC no OOB Pileup correction for MC
    cuts.AddCutCalo("10130053","411790105xe30220000","01331031000000d0"); // 00-10%
  } else if (trainConfig == 1107){ // central, with NOC NonLin like MC no OOB Pileup correction for MC with added signal
    cuts.AddCutCalo("10130023","411790105xe30220000","01331031000000d0"); // 00-10%
  } else if (trainConfig == 1108){ // central, with NOC NonLin like MC no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutCalo("10130023","411790105xe30220000","01331031000000d0"); // 00-10%
  } else if (trainConfig == 1109){ // central, with NOC NonLin like MC no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutCalo("10130023","411790105xe30220000","01331031000000d0"); // 00-10%
  } else if (trainConfig == 1110){ // central with mean N matched tracks per cluster data
    cuts.AddCutCalo("10130e03","411790105te30220000","01331031000000d0"); // 00-10%
  } else if (trainConfig == 1111){ // central, with mean N matched tracks per cluster MC no OOB Pileup correction for MC
    cuts.AddCutCalo("10130053","411790105te30220000","01331031000000d0"); // 00-10%
  } else if (trainConfig == 1112){ // central, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal
    cuts.AddCutCalo("10130023","411790105te30220000","01331031000000d0"); // 00-10%
  } else if (trainConfig == 1113){ // central, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutCalo("10130023","411790105te30220000","01331031000000d0"); // 00-10%
  } else if (trainConfig == 1114){ // central, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutCalo("10130023","411790105te30220000","01331031000000d0"); // 00-10%
  } else if (trainConfig == 1115){ // central with mean N matched tracks per cluster data
    cuts.AddCutCalo("10130e03","411790105ye30220000","01331031000000d0"); // 00-10%
  } else if (trainConfig == 1116){ // central, with mean N matched tracks per cluster MC no OOB Pileup correction for MC
    cuts.AddCutCalo("10130053","411790105ye30220000","01331031000000d0"); // 00-10%
  // with rotation background
  } else if (trainConfig == 1150){ // central
    cuts.AddCutCalo("10130e03","411790105ke30220000","0s331031000000d0"); // 00-10%
  } else if (trainConfig == 1151){ // central, no OOB Pileup correction for MC
    cuts.AddCutCalo("10130053","411790105ke30220000","0s331031000000d0"); // 00-10%
  } else if (trainConfig == 1152){ // central, no OOB Pileup correction for MC with added signal
    cuts.AddCutCalo("10130023","411790105ke30220000","0s331031000000d0"); // 00-10%
  } else if (trainConfig == 1153){ // central, no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutCalo("10130023","411790105ke30220000","0s331031000000d0"); // 00-10%
  } else if (trainConfig == 1154){ // central, no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutCalo("10130023","411790105ke30220000","0s331031000000d0"); // 00-10%
  } else if (trainConfig == 1155){ // central with NOC NonLin like data
    cuts.AddCutCalo("10130e03","411790105we30220000","0s331031000000d0"); // 00-10%
  } else if (trainConfig == 1156){ // central, with NOC NonLin like MC no OOB Pileup correction for MC
    cuts.AddCutCalo("10130053","411790105xe30220000","0s331031000000d0"); // 00-10%
  } else if (trainConfig == 1157){ // central, with NOC NonLin like MC no OOB Pileup correction for MC with added signal
    cuts.AddCutCalo("10130023","411790105xe30220000","0s331031000000d0"); // 00-10%
  } else if (trainConfig == 1158){ // central, with NOC NonLin like MC no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutCalo("10130023","411790105xe30220000","0s331031000000d0"); // 00-10%
  } else if (trainConfig == 1159){ // central, with NOC NonLin like MC no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutCalo("10130023","411790105xe30220000","0s331031000000d0"); // 00-10%
  } else if (trainConfig == 1160){ // central with mean N matched tracks per cluster data
    cuts.AddCutCalo("10130e03","411790105te30220000","0s331031000000d0"); // 00-10%
  } else if (trainConfig == 1161){ // central, with mean N matched tracks per cluster MC no OOB Pileup correction for MC
    cuts.AddCutCalo("10130053","411790105te30220000","0s331031000000d0"); // 00-10%
  } else if (trainConfig == 1162){ // central, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal
    cuts.AddCutCalo("10130023","411790105te30220000","0s331031000000d0"); // 00-10%
  } else if (trainConfig == 1163){ // central, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutCalo("10130023","411790105te30220000","0s331031000000d0"); // 00-10%
  } else if (trainConfig == 1164){ // central, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutCalo("10130023","411790105te30220000","0s331031000000d0"); // 00-10%
  } else if (trainConfig == 1165){ // central with mean N matched tracks per cluster data
    cuts.AddCutCalo("10130e03","411790105ye30220000","0s331031000000d0"); // 00-10%
  } else if (trainConfig == 1166){ // central, with mean N matched tracks per cluster MC no OOB Pileup correction for MC
    cuts.AddCutCalo("10130053","411790105ye30220000","0s331031000000d0"); // 00-10%
  // **********************************************************************************************************
  // ************************** EMC configurations PbPb run 2 2018 pass 3 semicent ****************************
  // **********************************************************************************************************
  // Standard: EMCal + DCal, only 2ndary cluster track matching, mixed event
  } else if (trainConfig == 1300){ // semicent
    cuts.AddCutCalo("11310e03","411790105ke30220000","01331031000000d0"); // 10-30%
  } else if (trainConfig == 1301){ // semicent, no OOB Pileup correction for MC
    cuts.AddCutCalo("11310053","411790105ke30220000","01331031000000d0"); // 10-30%
  } else if (trainConfig == 1302){ // semicent, no OOB Pileup correction for MC with added signal
    cuts.AddCutCalo("11310023","411790105ke30220000","01331031000000d0"); // 10-30%
  } else if (trainConfig == 1303){ // semicent, no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutCalo("11310023","411790105ke30220000","01331031000000d0"); // 10-30%
  } else if (trainConfig == 1304){ // semicent, no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutCalo("11310023","411790105ke30220000","01331031000000d0"); // 10-30%
  } else if (trainConfig == 1305){ // semicent with NOC NonLin like data
    cuts.AddCutCalo("11310e03","411790105we30220000","01331031000000d0"); // 10-30%
  } else if (trainConfig == 1306){ // semicent, with NOC NonLin like MC no OOB Pileup correction for MC
    cuts.AddCutCalo("11310053","411790105xe30220000","01331031000000d0"); // 10-30%
  } else if (trainConfig == 1307){ // semicent, with NOC NonLin like MC no OOB Pileup correction for MC with added signal
    cuts.AddCutCalo("11310023","411790105xe30220000","01331031000000d0"); // 10-30%
  } else if (trainConfig == 1308){ // semicent, with NOC NonLin like MC no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutCalo("11310023","411790105xe30220000","01331031000000d0"); // 10-30%
  } else if (trainConfig == 1309){ // semicent, with NOC NonLin like MC no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutCalo("11310023","411790105xe30220000","01331031000000d0"); // 10-30%
  } else if (trainConfig == 1310){ // semicent with mean N matched tracks per cluster data
    cuts.AddCutCalo("11310e03","411790105te30220000","01331031000000d0"); // 10-30%
  } else if (trainConfig == 1311){ // semicent, with mean N matched tracks per cluster MC no OOB Pileup correction for MC
    cuts.AddCutCalo("11310053","411790105te30220000","01331031000000d0"); // 10-30%
  } else if (trainConfig == 1312){ // semicent, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal
    cuts.AddCutCalo("11310023","411790105te30220000","01331031000000d0"); // 10-30%
  } else if (trainConfig == 1313){ // semicent, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutCalo("11310023","411790105te30220000","01331031000000d0"); // 10-30%
  } else if (trainConfig == 1314){ // semicent, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutCalo("11310023","411790105te30220000","01331031000000d0"); // 10-30%
  } else if (trainConfig == 1315){ // semicent with mean N matched tracks per cluster data
    cuts.AddCutCalo("11310e03","411790105ye30220000","01331031000000d0"); // 10-30%
  } else if (trainConfig == 1316){ // semicent, with mean N matched tracks per cluster MC no OOB Pileup correction for MC
    cuts.AddCutCalo("11310053","411790105ye30220000","01331031000000d0"); // 10-30%
  // with rotation background
  } else if (trainConfig == 1350){ // semicent
    cuts.AddCutCalo("11310e03","411790105ke30220000","0s331031000000d0"); // 10-30%
  } else if (trainConfig == 1351){ // semicent, no OOB Pileup correction for MC
    cuts.AddCutCalo("11310053","411790105ke30220000","0s331031000000d0"); // 10-30%
  } else if (trainConfig == 1352){ // semicent, no OOB Pileup correction for MC with added signal
    cuts.AddCutCalo("11310023","411790105ke30220000","0s331031000000d0"); // 10-30%
  } else if (trainConfig == 1353){ // semicent, no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutCalo("11310023","411790105ke30220000","0s331031000000d0"); // 10-30%
  } else if (trainConfig == 1354){ // semicent, no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutCalo("11310023","411790105ke30220000","0s331031000000d0"); // 10-30%
  } else if (trainConfig == 1355){ // semicent with NOC NonLin like data
    cuts.AddCutCalo("11310e03","411790105we30220000","0s331031000000d0"); // 10-30%
  } else if (trainConfig == 1356){ // semicent, with NOC NonLin like MC no OOB Pileup correction for MC
    cuts.AddCutCalo("11310053","411790105xe30220000","0s331031000000d0"); // 10-30%
  } else if (trainConfig == 1357){ // semicent, with NOC NonLin like MC no OOB Pileup correction for MC with added signal
    cuts.AddCutCalo("11310023","411790105xe30220000","0s331031000000d0"); // 10-30%
  } else if (trainConfig == 1358){ // semicent, with NOC NonLin like MC no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutCalo("11310023","411790105xe30220000","0s331031000000d0"); // 10-30%
  } else if (trainConfig == 1359){ // semicent, with NOC NonLin like MC no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutCalo("11310023","411790105xe30220000","0s331031000000d0"); // 10-30%
  } else if (trainConfig == 1360){ // semicent with mean N matched tracks per cluster data
    cuts.AddCutCalo("11310e03","411790105te30220000","0s331031000000d0"); // 10-30%
  } else if (trainConfig == 1361){ // semicent, with mean N matched tracks per cluster MC no OOB Pileup correction for MC
    cuts.AddCutCalo("11310053","411790105te30220000","0s331031000000d0"); // 10-30%
  } else if (trainConfig == 1362){ // semicent, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal
    cuts.AddCutCalo("11310023","411790105te30220000","0s331031000000d0"); // 10-30%
  } else if (trainConfig == 1363){ // semicent, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutCalo("11310023","411790105te30220000","0s331031000000d0"); // 10-30%
  } else if (trainConfig == 1364){ // semicent, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutCalo("11310023","411790105te30220000","0s331031000000d0"); // 10-30%
  } else if (trainConfig == 1365){ // semicent with mean N matched tracks per cluster data
    cuts.AddCutCalo("11310e03","411790105ye30220000","0s331031000000d0"); // 10-30%
  } else if (trainConfig == 1366){ // semicent, with mean N matched tracks per cluster MC no OOB Pileup correction for MC
    cuts.AddCutCalo("11310053","411790105ye30220000","0s331031000000d0"); // 10-30%
  // **********************************************************************************************************
  // *********************** EMC configurations PbPb run 2 2018 pass 3 semiperipheral *************************
  // **********************************************************************************************************
  // Standard: EMCal + DCal, only 2ndary cluster track matching, mixed event
  } else if (trainConfig == 1500){ // semiperipheral
    cuts.AddCutCalo("13530e03","411790105ke30220000","01331031000000d0"); // 30-50%
  } else if (trainConfig == 1501){ // semiperipheral, no OOB Pileup correction for MC
    cuts.AddCutCalo("13530053","411790105ke30220000","01331031000000d0"); // 30-50%
  } else if (trainConfig == 1502){ // semiperipheral, no OOB Pileup correction for MC with added signal
    cuts.AddCutCalo("13530023","411790105ke30220000","01331031000000d0"); // 30-50%
  } else if (trainConfig == 1503){ // semiperipheral, no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutCalo("13530023","411790105ke30220000","01331031000000d0"); // 30-50%
  } else if (trainConfig == 1504){ // semiperipheral, no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutCalo("13530023","411790105ke30220000","01331031000000d0"); // 30-50%
  } else if (trainConfig == 1505){ // semiperipheral with NOC NonLin like data
    cuts.AddCutCalo("13530e03","411790105we30220000","01331031000000d0"); // 30-50%
  } else if (trainConfig == 1506){ // semiperipheral, with NOC NonLin like MC no OOB Pileup correction for MC
    cuts.AddCutCalo("13530053","411790105xe30220000","01331031000000d0"); // 30-50%
  } else if (trainConfig == 1507){ // semiperipheral, with NOC NonLin like MC no OOB Pileup correction for MC with added signal
    cuts.AddCutCalo("13530023","411790105xe30220000","01331031000000d0"); // 30-50%
  } else if (trainConfig == 1508){ // semiperipheral, with NOC NonLin like MC no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutCalo("13530023","411790105xe30220000","01331031000000d0"); // 30-50%
  } else if (trainConfig == 1509){ // semiperipheral, with NOC NonLin like MC no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutCalo("13530023","411790105xe30220000","01331031000000d0"); // 30-50%
  } else if (trainConfig == 1510){ // semiperipheral with mean N matched tracks per cluster data
    cuts.AddCutCalo("13530e03","411790105te30220000","01331031000000d0"); // 30-50%
  } else if (trainConfig == 1511){ // semiperipheral, with mean N matched tracks per cluster MC no OOB Pileup correction for MC
    cuts.AddCutCalo("13530053","411790105te30220000","01331031000000d0"); // 30-50%
  } else if (trainConfig == 1512){ // semiperipheral, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal
    cuts.AddCutCalo("13530023","411790105te30220000","01331031000000d0"); // 30-50%
  } else if (trainConfig == 1513){ // semiperipheral, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutCalo("13530023","411790105te30220000","01331031000000d0"); // 30-50%
  } else if (trainConfig == 1514){ // semiperipheral, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutCalo("13530023","411790105te30220000","01331031000000d0"); // 30-50%
  } else if (trainConfig == 1515){ // semiperipheral with mean N matched tracks per cluster data
    cuts.AddCutCalo("13530e03","411790105ye30220000","01331031000000d0"); // 30-50%
  } else if (trainConfig == 1516){ // semiperipheral, with mean N matched tracks per cluster MC no OOB Pileup correction for MC
    cuts.AddCutCalo("13530053","411790105ye30220000","01331031000000d0"); // 30-50%
  // with rotation background
  } else if (trainConfig == 1550){ // semiperipheral
    cuts.AddCutCalo("13530e03","411790105ke30220000","0s331031000000d0"); // 30-50%
  } else if (trainConfig == 1551){ // semiperipheral, no OOB Pileup correction for MC
    cuts.AddCutCalo("13530053","411790105ke30220000","0s331031000000d0"); // 30-50%
  } else if (trainConfig == 1552){ // semiperipheral, no OOB Pileup correction for MC with added signal
    cuts.AddCutCalo("13530023","411790105ke30220000","0s331031000000d0"); // 30-50%
  } else if (trainConfig == 1553){ // semiperipheral, no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutCalo("13530023","411790105ke30220000","0s331031000000d0"); // 30-50%
  } else if (trainConfig == 1554){ // semiperipheral, no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutCalo("13530023","411790105ke30220000","0s331031000000d0"); // 30-50%
  } else if (trainConfig == 1555){ // semiperipheral with NOC NonLin like data
    cuts.AddCutCalo("13530e03","411790105we30220000","0s331031000000d0"); // 30-50%
  } else if (trainConfig == 1556){ // semiperipheral, with NOC NonLin like MC no OOB Pileup correction for MC
    cuts.AddCutCalo("13530053","411790105xe30220000","0s331031000000d0"); // 30-50%
  } else if (trainConfig == 1557){ // semiperipheral, with NOC NonLin like MC no OOB Pileup correction for MC with added signal
    cuts.AddCutCalo("13530023","411790105xe30220000","0s331031000000d0"); // 30-50%
  } else if (trainConfig == 1558){ // semiperipheral, with NOC NonLin like MC no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutCalo("13530023","411790105xe30220000","0s331031000000d0"); // 30-50%
  } else if (trainConfig == 1559){ // semiperipheral, with NOC NonLin like MC no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutCalo("13530023","411790105xe30220000","0s331031000000d0"); // 30-50%
  } else if (trainConfig == 1560){ // semiperipheral with mean N matched tracks per cluster data
    cuts.AddCutCalo("13530e03","411790105te30220000","0s331031000000d0"); // 30-50%
  } else if (trainConfig == 1561){ // semiperipheral, with mean N matched tracks per cluster MC no OOB Pileup correction for MC
    cuts.AddCutCalo("13530053","411790105te30220000","0s331031000000d0"); // 30-50%
  } else if (trainConfig == 1562){ // semiperipheral, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal
    cuts.AddCutCalo("13530023","411790105te30220000","0s331031000000d0"); // 30-50%
  } else if (trainConfig == 1563){ // semiperipheral, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutCalo("13530023","411790105te30220000","0s331031000000d0"); // 30-50%
  } else if (trainConfig == 1564){ // semiperipheral, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutCalo("13530023","411790105te30220000","0s331031000000d0"); // 30-50%
  } else if (trainConfig == 1565){ // semiperipheral with mean N matched tracks per cluster data
    cuts.AddCutCalo("13530e03","411790105ye30220000","0s331031000000d0"); // 30-50%
  } else if (trainConfig == 1566){ // semiperipheral, with mean N matched tracks per cluster MC no OOB Pileup correction for MC
    cuts.AddCutCalo("13530053","411790105ye30220000","0s331031000000d0"); // 30-50%
  // **********************************************************************************************************
  // ************************* EMC configurations PbPb run 2 2018 pass 3 peripheral ***************************
  // **********************************************************************************************************
  // Standard: EMCal + DCal, only 2ndary cluster track matching, mixed event
  } else if (trainConfig == 1900){ // peripheral
    cuts.AddCutCalo("15910e03","411790105ke30220000","01331031000000d0"); // 50-90%
  } else if (trainConfig == 1901){ // peripheral, no OOB Pileup correction for MC
    cuts.AddCutCalo("15910053","411790105ke30220000","01331031000000d0"); // 50-90%
  } else if (trainConfig == 1902){ // peripheral, no OOB Pileup correction for MC with added signal
    cuts.AddCutCalo("15910023","411790105ke30220000","01331031000000d0"); // 50-90%
  } else if (trainConfig == 1903){ // peripheral, no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutCalo("15910023","411790105ke30220000","01331031000000d0"); // 50-90%
  } else if (trainConfig == 1904){ // peripheral, no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutCalo("15910023","411790105ke30220000","01331031000000d0"); // 50-90%
  } else if (trainConfig == 1905){ // peripheral with NOC NonLin like data
    cuts.AddCutCalo("15910e03","411790105we30220000","01331031000000d0"); // 50-90%
  } else if (trainConfig == 1906){ // peripheral, with NOC NonLin like MC no OOB Pileup correction for MC
    cuts.AddCutCalo("15910053","411790105xe30220000","01331031000000d0"); // 50-90%
  } else if (trainConfig == 1907){ // peripheral, with NOC NonLin like MC no OOB Pileup correction for MC with added signal
    cuts.AddCutCalo("15910023","411790105xe30220000","01331031000000d0"); // 50-90%
  } else if (trainConfig == 1908){ // peripheral, with NOC NonLin like MC no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutCalo("15910023","411790105xe30220000","01331031000000d0"); // 50-90%
  } else if (trainConfig == 1909){ // peripheral, with NOC NonLin like MC no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutCalo("15910023","411790105xe30220000","01331031000000d0"); // 50-90%
  } else if (trainConfig == 1910){ // peripheral with mean N matched tracks per cluster data
    cuts.AddCutCalo("15910e03","411790105te30220000","01331031000000d0"); // 50-90%
  } else if (trainConfig == 1911){ // peripheral, with mean N matched tracks per cluster MC no OOB Pileup correction for MC
    cuts.AddCutCalo("15910053","411790105te30220000","01331031000000d0"); // 50-90%
  } else if (trainConfig == 1912){ // peripheral, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal
    cuts.AddCutCalo("15910023","411790105te30220000","01331031000000d0"); // 50-90%
  } else if (trainConfig == 1913){ // peripheral, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutCalo("15910023","411790105te30220000","01331031000000d0"); // 50-90%
  } else if (trainConfig == 1914){ // peripheral, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutCalo("15910023","411790105te30220000","01331031000000d0"); // 50-90%
  } else if (trainConfig == 1915){ // peripheral with mean N matched tracks per cluster data
    cuts.AddCutCalo("15910e03","411790105ye30220000","01331031000000d0"); // 50-90%
  } else if (trainConfig == 1916){ // peripheral, with mean N matched tracks per cluster MC no OOB Pileup correction for MC
    cuts.AddCutCalo("15910053","411790105ye30220000","01331031000000d0"); // 50-90%
  // with rotation background
  } else if (trainConfig == 1950){ // peripheral
    cuts.AddCutCalo("15910e03","411790105ke30220000","0s331031000000d0"); // 50-90%
  } else if (trainConfig == 1951){ // peripheral, no OOB Pileup correction for MC
    cuts.AddCutCalo("15910053","411790105ke30220000","0s331031000000d0"); // 50-90%
  } else if (trainConfig == 1952){ // peripheral, no OOB Pileup correction for MC with added signal
    cuts.AddCutCalo("15910023","411790105ke30220000","0s331031000000d0"); // 50-90%
  } else if (trainConfig == 1953){ // peripheral, no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutCalo("15910023","411790105ke30220000","0s331031000000d0"); // 50-90%
  } else if (trainConfig == 1954){ // peripheral, no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutCalo("15910023","411790105ke30220000","0s331031000000d0"); // 50-90%
  } else if (trainConfig == 1955){ // peripheral with NOC NonLin like data
    cuts.AddCutCalo("15910e03","411790105we30220000","0s331031000000d0"); // 50-90%
  } else if (trainConfig == 1956){ // peripheral, with NOC NonLin like MC no OOB Pileup correction for MC
    cuts.AddCutCalo("15910053","411790105xe30220000","0s331031000000d0"); // 50-90%
  } else if (trainConfig == 1957){ // peripheral, with NOC NonLin like MC no OOB Pileup correction for MC with added signal
    cuts.AddCutCalo("15910023","411790105xe30220000","0s331031000000d0"); // 50-90%
  } else if (trainConfig == 1958){ // peripheral, with NOC NonLin like MC no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutCalo("15910023","411790105xe30220000","0s331031000000d0"); // 50-90%
  } else if (trainConfig == 1959){ // peripheral, with NOC NonLin like MC no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutCalo("15910023","411790105xe30220000","0s331031000000d0"); // 50-90%
  } else if (trainConfig == 1960){ // peripheral with mean N matched tracks per cluster data
    cuts.AddCutCalo("15910e03","411790105te30220000","0s331031000000d0"); // 50-90%
  } else if (trainConfig == 1961){ // peripheral, with mean N matched tracks per cluster MC no OOB Pileup correction for MC
    cuts.AddCutCalo("15910053","411790105te30220000","0s331031000000d0"); // 50-90%
  } else if (trainConfig == 1962){ // peripheral, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal
    cuts.AddCutCalo("15910023","411790105te30220000","0s331031000000d0"); // 50-90%
  } else if (trainConfig == 1963){ // peripheral, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutCalo("15910023","411790105te30220000","0s331031000000d0"); // 50-90%
  } else if (trainConfig == 1964){ // peripheral, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutCalo("15910023","411790105te30220000","0s331031000000d0"); // 50-90%
  } else if (trainConfig == 1965){ // peripheral with mean N matched tracks per cluster data
    cuts.AddCutCalo("15910e03","411790105ye30220000","0s331031000000d0"); // 50-90%
  } else if (trainConfig == 1966){ // peripheral, with mean N matched tracks per cluster MC no OOB Pileup correction for MC
    cuts.AddCutCalo("15910053","411790105ye30220000","0s331031000000d0"); // 50-90%
    // **********************************************************************************************************
  // ***************************** EMC configurations PbPb run 2 2015 pass 3 *******************************
  // **********************************************************************************************************
  // Standard: EMCal + DCal, only 2ndary cluster track matching, mixed event
  } else if (trainConfig == 2000){ // 4 cents
    cuts.AddCutCalo("10110e03","411790105ke30220000","01331031000000d0"); // 00-10%
    cuts.AddCutCalo("11310e03","411790105ke30220000","01331031000000d0"); // 10-30%
    cuts.AddCutCalo("13510e03","411790105ke30220000","01331031000000d0"); // 30-50%
    cuts.AddCutCalo("15910e03","411790105ke30220000","01331031000000d0"); // 50-90%
  } else if (trainConfig == 2001){ // 4 cents, no OOB Pileup correction for MC
    cuts.AddCutCalo("10110053","411790105ke30220000","01331031000000d0"); // 00-10%
    cuts.AddCutCalo("11310053","411790105ke30220000","01331031000000d0"); // 10-30%
    cuts.AddCutCalo("13510053","411790105ke30220000","01331031000000d0"); // 30-50%
    cuts.AddCutCalo("15910053","411790105ke30220000","01331031000000d0"); // 50-90%
  } else if (trainConfig == 2002){ // 4 cents, no OOB Pileup correction for MC with added signal
    cuts.AddCutCalo("10110023","411790105ke30220000","01331031000000d0"); // 00-10%
    cuts.AddCutCalo("11310023","411790105ke30220000","01331031000000d0"); // 10-30%
    cuts.AddCutCalo("13510023","411790105ke30220000","01331031000000d0"); // 30-50%
    cuts.AddCutCalo("15910023","411790105ke30220000","01331031000000d0"); // 50-90%
  } else if (trainConfig == 2003){ // 4 cents, no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutCalo("10110023","411790105ke30220000","01331031000000d0"); // 00-10%
    cuts.AddCutCalo("11310023","411790105ke30220000","01331031000000d0"); // 10-30%
    cuts.AddCutCalo("13510023","411790105ke30220000","01331031000000d0"); // 30-50%
    cuts.AddCutCalo("15910023","411790105ke30220000","01331031000000d0"); // 50-90%
  } else if (trainConfig == 2004){ // 4 cents, no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutCalo("10110023","411790105ke30220000","01331031000000d0"); // 00-10%
    cuts.AddCutCalo("11310023","411790105ke30220000","01331031000000d0"); // 10-30%
    cuts.AddCutCalo("13510023","411790105ke30220000","01331031000000d0"); // 30-50%
    cuts.AddCutCalo("15910023","411790105ke30220000","01331031000000d0"); // 50-90%
  } else if (trainConfig == 2005){ // 4 cents with NOC NonLin like data
    cuts.AddCutCalo("10110e03","411790105we30220000","01331031000000d0"); // 00-10%
    cuts.AddCutCalo("11310e03","411790105we30220000","01331031000000d0"); // 10-30%
    cuts.AddCutCalo("13510e03","411790105we30220000","01331031000000d0"); // 30-50%
    cuts.AddCutCalo("15910e03","411790105we30220000","01331031000000d0"); // 50-90%
  } else if (trainConfig == 2006){ // 4 cents, with NOC NonLin like MC no OOB Pileup correction for MC
    cuts.AddCutCalo("10110053","411790105xe30220000","01331031000000d0"); // 00-10%
    cuts.AddCutCalo("11310053","411790105xe30220000","01331031000000d0"); // 10-30%
    cuts.AddCutCalo("13510053","411790105xe30220000","01331031000000d0"); // 30-50%
    cuts.AddCutCalo("15910053","411790105xe30220000","01331031000000d0"); // 50-90%
  } else if (trainConfig == 2007){ // 4 cents, with NOC NonLin like MC no OOB Pileup correction for MC with added signal
    cuts.AddCutCalo("10110023","411790105xe30220000","01331031000000d0"); // 00-10%
    cuts.AddCutCalo("11310023","411790105xe30220000","01331031000000d0"); // 10-30%
    cuts.AddCutCalo("13510023","411790105xe30220000","01331031000000d0"); // 30-50%
    cuts.AddCutCalo("15910023","411790105xe30220000","01331031000000d0"); // 50-90%
  } else if (trainConfig == 2008){ // 4 cents, with NOC NonLin like MC no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutCalo("10110023","411790105xe30220000","01331031000000d0"); // 00-10%
    cuts.AddCutCalo("11310023","411790105xe30220000","01331031000000d0"); // 10-30%
    cuts.AddCutCalo("13510023","411790105xe30220000","01331031000000d0"); // 30-50%
    cuts.AddCutCalo("15910023","411790105xe30220000","01331031000000d0"); // 50-90%
  } else if (trainConfig == 2009){ // 4 cents, with NOC NonLin like MC no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutCalo("10110023","411790105xe30220000","01331031000000d0"); // 00-10%
    cuts.AddCutCalo("11310023","411790105xe30220000","01331031000000d0"); // 10-30%
    cuts.AddCutCalo("13510023","411790105xe30220000","01331031000000d0"); // 30-50%
    cuts.AddCutCalo("15910023","411790105xe30220000","01331031000000d0"); // 50-90%
  } else if (trainConfig == 2010){ // 4 cents with mean N matched tracks per cluster data
    cuts.AddCutCalo("10110e03","411790105te30220000","01331031000000d0"); // 00-10%
    cuts.AddCutCalo("11310e03","411790105te30220000","01331031000000d0"); // 10-30%
    cuts.AddCutCalo("13510e03","411790105te30220000","01331031000000d0"); // 30-50%
    cuts.AddCutCalo("15910e03","411790105te30220000","01331031000000d0"); // 50-90%
  } else if (trainConfig == 2011){ // 4 cents, with mean N matched tracks per cluster MC no OOB Pileup correction for MC
    cuts.AddCutCalo("10110053","411790105te30220000","01331031000000d0"); // 00-10%
    cuts.AddCutCalo("11310053","411790105te30220000","01331031000000d0"); // 10-30%
    cuts.AddCutCalo("13510053","411790105te30220000","01331031000000d0"); // 30-50%
    cuts.AddCutCalo("15910053","411790105te30220000","01331031000000d0"); // 50-90%
  } else if (trainConfig == 2012){ // 4 cents, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal
    cuts.AddCutCalo("10110023","411790105te30220000","01331031000000d0"); // 00-10%
    cuts.AddCutCalo("11310023","411790105te30220000","01331031000000d0"); // 10-30%
    cuts.AddCutCalo("13510023","411790105te30220000","01331031000000d0"); // 30-50%
    cuts.AddCutCalo("15910023","411790105te30220000","01331031000000d0"); // 50-90%
  } else if (trainConfig == 2013){ // 4 cents, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutCalo("10110023","411790105te30220000","01331031000000d0"); // 00-10%
    cuts.AddCutCalo("11310023","411790105te30220000","01331031000000d0"); // 10-30%
    cuts.AddCutCalo("13510023","411790105te30220000","01331031000000d0"); // 30-50%
    cuts.AddCutCalo("15910023","411790105te30220000","01331031000000d0"); // 50-90%
  } else if (trainConfig == 2014){ // 4 cents, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutCalo("10110023","411790105te30220000","01331031000000d0"); // 00-10%
    cuts.AddCutCalo("11310023","411790105te30220000","01331031000000d0"); // 10-30%
    cuts.AddCutCalo("13510023","411790105te30220000","01331031000000d0"); // 30-50%
    cuts.AddCutCalo("15910023","411790105te30220000","01331031000000d0"); // 50-90%
  } else if (trainConfig == 2015){ // 4 cents with NOC via <E_clus-E>/E
    cuts.AddCutCalo("10110e03","411790105ye30220000","01331031000000d0"); // 00-10%
    cuts.AddCutCalo("11310e03","411790105ye30220000","01331031000000d0"); // 10-30%
    cuts.AddCutCalo("13510e03","411790105ye30220000","01331031000000d0"); // 30-50%
    cuts.AddCutCalo("15910e03","411790105ye30220000","01331031000000d0"); // 50-90%
  } else if (trainConfig == 2016){ // 4 cents, with NOC via <E_clus-E>/E no OOB Pileup correction for MC
    cuts.AddCutCalo("10110053","411790105ye30220000","01331031000000d0"); // 00-10%
    cuts.AddCutCalo("11310053","411790105ye30220000","01331031000000d0"); // 10-30%
    cuts.AddCutCalo("13510053","411790105ye30220000","01331031000000d0"); // 30-50%
    cuts.AddCutCalo("15910053","411790105ye30220000","01331031000000d0"); // 50-90%
  // with rotation background
  } else if (trainConfig == 2050){ // 4 cents
    cuts.AddCutCalo("10110e03","411790105ke30220000","0s331031000000d0"); // 00-10%
    cuts.AddCutCalo("11310e03","411790105ke30220000","0s331031000000d0"); // 10-30%
    cuts.AddCutCalo("13510e03","411790105ke30220000","0s331031000000d0"); // 30-50%
    cuts.AddCutCalo("15910e03","411790105ke30220000","0s331031000000d0"); // 50-90%
  } else if (trainConfig == 2051){ // 4 cents, no OOB Pileup correction for MC
    cuts.AddCutCalo("10110053","411790105ke30220000","0s331031000000d0"); // 00-10%
    cuts.AddCutCalo("11310053","411790105ke30220000","0s331031000000d0"); // 10-30%
    cuts.AddCutCalo("13510053","411790105ke30220000","0s331031000000d0"); // 30-50%
    cuts.AddCutCalo("15910053","411790105ke30220000","0s331031000000d0"); // 50-90%
  } else if (trainConfig == 2052){ // 4 cents, no OOB Pileup correction for MC with added signal
    cuts.AddCutCalo("10110023","411790105ke30220000","0s331031000000d0"); // 00-10%
    cuts.AddCutCalo("11310023","411790105ke30220000","0s331031000000d0"); // 10-30%
    cuts.AddCutCalo("13510023","411790105ke30220000","0s331031000000d0"); // 30-50%
    cuts.AddCutCalo("15910023","411790105ke30220000","0s331031000000d0"); // 50-90%
  } else if (trainConfig == 2053){ // 4 cents, no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutCalo("10110023","411790105ke30220000","0s331031000000d0"); // 00-10%
    cuts.AddCutCalo("11310023","411790105ke30220000","0s331031000000d0"); // 10-30%
    cuts.AddCutCalo("13510023","411790105ke30220000","0s331031000000d0"); // 30-50%
    cuts.AddCutCalo("15910023","411790105ke30220000","0s331031000000d0"); // 50-90%
  } else if (trainConfig == 2054){ // 4 cents, no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutCalo("10110023","411790105ke30220000","0s331031000000d0"); // 00-10%
    cuts.AddCutCalo("11310023","411790105ke30220000","0s331031000000d0"); // 10-30%
    cuts.AddCutCalo("13510023","411790105ke30220000","0s331031000000d0"); // 30-50%
    cuts.AddCutCalo("15910023","411790105ke30220000","0s331031000000d0"); // 50-90%
  } else if (trainConfig == 2055){ // 4 cents with NOC NonLin like data
    cuts.AddCutCalo("10110e03","411790105we30220000","0s331031000000d0"); // 00-10%
    cuts.AddCutCalo("11310e03","411790105we30220000","0s331031000000d0"); // 10-30%
    cuts.AddCutCalo("13510e03","411790105we30220000","0s331031000000d0"); // 30-50%
    cuts.AddCutCalo("15910e03","411790105we30220000","0s331031000000d0"); // 50-90%
  } else if (trainConfig == 2056){ // 4 cents, with NOC NonLin like MC no OOB Pileup correction for MC
    cuts.AddCutCalo("10110053","411790105xe30220000","0s331031000000d0"); // 00-10%
    cuts.AddCutCalo("11310053","411790105xe30220000","0s331031000000d0"); // 10-30%
    cuts.AddCutCalo("13510053","411790105xe30220000","0s331031000000d0"); // 30-50%
    cuts.AddCutCalo("15910053","411790105xe30220000","0s331031000000d0"); // 50-90%
  } else if (trainConfig == 2057){ // 4 cents, with NOC NonLin like MC no OOB Pileup correction for MC with added signal
    cuts.AddCutCalo("10110023","411790105xe30220000","0s331031000000d0"); // 00-10%
    cuts.AddCutCalo("11310023","411790105xe30220000","0s331031000000d0"); // 10-30%
    cuts.AddCutCalo("13510023","411790105xe30220000","0s331031000000d0"); // 30-50%
    cuts.AddCutCalo("15910023","411790105xe30220000","0s331031000000d0"); // 50-90%
  } else if (trainConfig == 2058){ // 4 cents, with NOC NonLin like MC no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutCalo("10110023","411790105xe30220000","0s331031000000d0"); // 00-10%
    cuts.AddCutCalo("11310023","411790105xe30220000","0s331031000000d0"); // 10-30%
    cuts.AddCutCalo("13510023","411790105xe30220000","0s331031000000d0"); // 30-50%
    cuts.AddCutCalo("15910023","411790105xe30220000","0s331031000000d0"); // 50-90%
  } else if (trainConfig == 2059){ // 4 cents, with NOC NonLin like MC no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutCalo("10110023","411790105xe30220000","0s331031000000d0"); // 00-10%
    cuts.AddCutCalo("11310023","411790105xe30220000","0s331031000000d0"); // 10-30%
    cuts.AddCutCalo("13510023","411790105xe30220000","0s331031000000d0"); // 30-50%
    cuts.AddCutCalo("15910023","411790105xe30220000","0s331031000000d0"); // 50-90%
  } else if (trainConfig == 2060){ // 4 cents with mean N matched tracks per cluster data
    cuts.AddCutCalo("10110e03","411790105te30220000","0s331031000000d0"); // 00-10%
    cuts.AddCutCalo("11310e03","411790105te30220000","0s331031000000d0"); // 10-30%
    cuts.AddCutCalo("13510e03","411790105te30220000","0s331031000000d0"); // 30-50%
    cuts.AddCutCalo("15910e03","411790105te30220000","0s331031000000d0"); // 50-90%
  } else if (trainConfig == 2061){ // 4 cents, with mean N matched tracks per cluster MC no OOB Pileup correction for MC
    cuts.AddCutCalo("10110053","411790105te30220000","0s331031000000d0"); // 00-10%
    cuts.AddCutCalo("11310053","411790105te30220000","0s331031000000d0"); // 10-30%
    cuts.AddCutCalo("13510053","411790105te30220000","0s331031000000d0"); // 30-50%
    cuts.AddCutCalo("15910053","411790105te30220000","0s331031000000d0"); // 50-90%
  } else if (trainConfig == 2062){ // 4 cents, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal
    cuts.AddCutCalo("10110023","411790105te30220000","0s331031000000d0"); // 00-10%
    cuts.AddCutCalo("11310023","411790105te30220000","0s331031000000d0"); // 10-30%
    cuts.AddCutCalo("13510023","411790105te30220000","0s331031000000d0"); // 30-50%
    cuts.AddCutCalo("15910023","411790105te30220000","0s331031000000d0"); // 50-90%
  } else if (trainConfig == 2063){ // 4 cents, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutCalo("10110023","411790105te30220000","0s331031000000d0"); // 00-10%
    cuts.AddCutCalo("11310023","411790105te30220000","0s331031000000d0"); // 10-30%
    cuts.AddCutCalo("13510023","411790105te30220000","0s331031000000d0"); // 30-50%
    cuts.AddCutCalo("15910023","411790105te30220000","0s331031000000d0"); // 50-90%
  } else if (trainConfig == 2064){ // 4 cents, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutCalo("10110023","411790105te30220000","0s331031000000d0"); // 00-10%
    cuts.AddCutCalo("11310023","411790105te30220000","0s331031000000d0"); // 10-30%
    cuts.AddCutCalo("13510023","411790105te30220000","0s331031000000d0"); // 30-50%
    cuts.AddCutCalo("15910023","411790105te30220000","0s331031000000d0"); // 50-90%
  } else if (trainConfig == 2065){ // 4 cents with NOC via <E_clus-E>/E
    cuts.AddCutCalo("10110e03","411790105ye30220000","0s331031000000d0"); // 00-10%
    cuts.AddCutCalo("11310e03","411790105ye30220000","0s331031000000d0"); // 10-30%
    cuts.AddCutCalo("13510e03","411790105ye30220000","0s331031000000d0"); // 30-50%
    cuts.AddCutCalo("15910e03","411790105ye30220000","0s331031000000d0"); // 50-90%
  } else if (trainConfig == 2066){ // 4 cents, with NOC via <E_clus-E>/E no OOB Pileup correction for MC
    cuts.AddCutCalo("10110053","411790105ye30220000","0s331031000000d0"); // 00-10%
    cuts.AddCutCalo("11310053","411790105ye30220000","0s331031000000d0"); // 10-30%
    cuts.AddCutCalo("13510053","411790105ye30220000","0s331031000000d0"); // 30-50%
    cuts.AddCutCalo("15910053","411790105ye30220000","0s331031000000d0"); // 50-90%
  // **********************************************************************************************************
  // **************************** EMC configurations PbPb run 2 2015 pass 3 cent ******************************
  // **********************************************************************************************************
  // Standard: EMCal + DCal, only 2ndary cluster track matching, mixed event
  } else if (trainConfig == 2100){ // central
    cuts.AddCutCalo("10110e03","411790105ke30220000","01331031000000d0"); // 00-10%
  } else if (trainConfig == 2101){ // central, no OOB Pileup correction for MC
    cuts.AddCutCalo("10110053","411790105ke30220000","01331031000000d0"); // 00-10%
  } else if (trainConfig == 2102){ // central, no OOB Pileup correction for MC with added signal
    cuts.AddCutCalo("10110023","411790105ke30220000","01331031000000d0"); // 00-10%
  } else if (trainConfig == 2103){ // central, no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutCalo("10110023","411790105ke30220000","01331031000000d0"); // 00-10%
  } else if (trainConfig == 2104){ // central, no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutCalo("10110023","411790105ke30220000","01331031000000d0"); // 00-10%
  } else if (trainConfig == 2105){ // central with NOC NonLin like data
    cuts.AddCutCalo("10110e03","411790105we30220000","01331031000000d0"); // 00-10%
  } else if (trainConfig == 2106){ // central, with NOC NonLin like MC no OOB Pileup correction for MC
    cuts.AddCutCalo("10110053","411790105xe30220000","01331031000000d0"); // 00-10%
  } else if (trainConfig == 2107){ // central, with NOC NonLin like MC no OOB Pileup correction for MC with added signal
    cuts.AddCutCalo("10110023","411790105xe30220000","01331031000000d0"); // 00-10%
  } else if (trainConfig == 2108){ // central, with NOC NonLin like MC no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutCalo("10110023","411790105xe30220000","01331031000000d0"); // 00-10%
  } else if (trainConfig == 2109){ // central, with NOC NonLin like MC no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutCalo("10110023","411790105xe30220000","01331031000000d0"); // 00-10%
  } else if (trainConfig == 2110){ // central with mean N matched tracks per cluster data
    cuts.AddCutCalo("10110e03","411790105te30220000","01331031000000d0"); // 00-10%
  } else if (trainConfig == 2111){ // central, with mean N matched tracks per cluster MC no OOB Pileup correction for MC
    cuts.AddCutCalo("10110053","411790105te30220000","01331031000000d0"); // 00-10%
  } else if (trainConfig == 2112){ // central, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal
    cuts.AddCutCalo("10110023","411790105te30220000","01331031000000d0"); // 00-10%
  } else if (trainConfig == 2113){ // central, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutCalo("10110023","411790105te30220000","01331031000000d0"); // 00-10%
  } else if (trainConfig == 2114){ // central, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutCalo("10110023","411790105te30220000","01331031000000d0"); // 00-10%
  } else if (trainConfig == 2115){ // central NOC via <E_clus-E>/E data
    cuts.AddCutCalo("10110e03","411790105ye30220000","01331031000000d0"); // 00-10%
  } else if (trainConfig == 2116){ // central, NOC via <E_clus-E>/E MC no OOB Pileup correction for MC
    cuts.AddCutCalo("10110053","411790105ye30220000","01331031000000d0"); // 00-10%
  // with rotation background
  } else if (trainConfig == 2150){ // central
    cuts.AddCutCalo("10110e03","411790105ke30220000","0s331031000000d0"); // 00-10%
  } else if (trainConfig == 2151){ // central, no OOB Pileup correction for MC
    cuts.AddCutCalo("10110053","411790105ke30220000","0s331031000000d0"); // 00-10%
  } else if (trainConfig == 2152){ // central, no OOB Pileup correction for MC with added signal
    cuts.AddCutCalo("10110023","411790105ke30220000","0s331031000000d0"); // 00-10%
  } else if (trainConfig == 2153){ // central, no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutCalo("10110023","411790105ke30220000","0s331031000000d0"); // 00-10%
  } else if (trainConfig == 2154){ // central, no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutCalo("10110023","411790105ke30220000","0s331031000000d0"); // 00-10%
  } else if (trainConfig == 2155){ // central with NOC NonLin like data
    cuts.AddCutCalo("10110e03","411790105we30220000","0s331031000000d0"); // 00-10%
  } else if (trainConfig == 2156){ // central, with NOC NonLin like MC no OOB Pileup correction for MC
    cuts.AddCutCalo("10110053","411790105xe30220000","0s331031000000d0"); // 00-10%
  } else if (trainConfig == 2157){ // central, with NOC NonLin like MC no OOB Pileup correction for MC with added signal
    cuts.AddCutCalo("10110023","411790105xe30220000","0s331031000000d0"); // 00-10%
  } else if (trainConfig == 2158){ // central, with NOC NonLin like MC no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutCalo("10110023","411790105xe30220000","0s331031000000d0"); // 00-10%
  } else if (trainConfig == 2159){ // central, with NOC NonLin like MC no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutCalo("10110023","411790105xe30220000","0s331031000000d0"); // 00-10%
  } else if (trainConfig == 2160){ // central with mean N matched tracks per cluster data
    cuts.AddCutCalo("10110e03","411790105te30220000","0s331031000000d0"); // 00-10%
  } else if (trainConfig == 2161){ // central, with mean N matched tracks per cluster MC no OOB Pileup correction for MC
    cuts.AddCutCalo("10110053","411790105te30220000","0s331031000000d0"); // 00-10%
  } else if (trainConfig == 2162){ // central, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal
    cuts.AddCutCalo("10110023","411790105te30220000","0s331031000000d0"); // 00-10%
  } else if (trainConfig == 2163){ // central, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutCalo("10110023","411790105te30220000","0s331031000000d0"); // 00-10%
  } else if (trainConfig == 2164){ // central, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutCalo("10110023","411790105te30220000","0s331031000000d0"); // 00-10%
  } else if (trainConfig == 2165){ // central NOC via <E_clus-E>/E data
    cuts.AddCutCalo("10110e03","411790105ye30220000","0s331031000000d0"); // 00-10%
  } else if (trainConfig == 2166){ // central, NOC via <E_clus-E>/E MC no OOB Pileup correction for MC
    cuts.AddCutCalo("10110053","411790105ye30220000","0s331031000000d0"); // 00-10%
  // **********************************************************************************************************
  // ************************** EMC configurations PbPb run 2 2015 pass 3 semicent ****************************
  // **********************************************************************************************************
  // Standard: EMCal + DCal, only 2ndary cluster track matching, mixed event
  } else if (trainConfig == 2300){ // semicent
    cuts.AddCutCalo("11310e03","411790105ke30220000","01331031000000d0"); // 10-30%
  } else if (trainConfig == 2301){ // semicent, no OOB Pileup correction for MC
    cuts.AddCutCalo("11310053","411790105ke30220000","01331031000000d0"); // 10-30%
  } else if (trainConfig == 2302){ // semicent, no OOB Pileup correction for MC with added signal
    cuts.AddCutCalo("11310023","411790105ke30220000","01331031000000d0"); // 10-30%
  } else if (trainConfig == 2303){ // semicent, no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutCalo("11310023","411790105ke30220000","01331031000000d0"); // 10-30%
  } else if (trainConfig == 2304){ // semicent, no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutCalo("11310023","411790105ke30220000","01331031000000d0"); // 10-30%
  } else if (trainConfig == 2305){ // semicent with NOC NonLin like data
    cuts.AddCutCalo("11310e03","411790105we30220000","01331031000000d0"); // 10-30%
  } else if (trainConfig == 2306){ // semicent, with NOC NonLin like MC no OOB Pileup correction for MC
    cuts.AddCutCalo("11310053","411790105xe30220000","01331031000000d0"); // 10-30%
  } else if (trainConfig == 2307){ // semicent, with NOC NonLin like MC no OOB Pileup correction for MC with added signal
    cuts.AddCutCalo("11310023","411790105xe30220000","01331031000000d0"); // 10-30%
  } else if (trainConfig == 2308){ // semicent, with NOC NonLin like MC no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutCalo("11310023","411790105xe30220000","01331031000000d0"); // 10-30%
  } else if (trainConfig == 2309){ // semicent, with NOC NonLin like MC no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutCalo("11310023","411790105xe30220000","01331031000000d0"); // 10-30%
  } else if (trainConfig == 2310){ // semicent with mean N matched tracks per cluster data
    cuts.AddCutCalo("11310e03","411790105te30220000","01331031000000d0"); // 10-30%
  } else if (trainConfig == 2311){ // semicent, with mean N matched tracks per cluster MC no OOB Pileup correction for MC
    cuts.AddCutCalo("11310053","411790105te30220000","01331031000000d0"); // 10-30%
  } else if (trainConfig == 2312){ // semicent, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal
    cuts.AddCutCalo("11310023","411790105te30220000","01331031000000d0"); // 10-30%
  } else if (trainConfig == 2313){ // semicent, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutCalo("11310023","411790105te30220000","01331031000000d0"); // 10-30%
  } else if (trainConfig == 2314){ // semicent, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutCalo("11310023","411790105te30220000","01331031000000d0"); // 10-30%
  } else if (trainConfig == 2315){ // semicent NOC via <E_clus-E>/E data
    cuts.AddCutCalo("11310e03","411790105ye30220000","01331031000000d0"); // 10-30%
  } else if (trainConfig == 2316){ // semicent, NOC via <E_clus-E>/E MC no OOB Pileup correction for MC
    cuts.AddCutCalo("11310053","411790105ye30220000","01331031000000d0"); // 10-30%
  // with rotation background
  } else if (trainConfig == 2350){ // semicent
    cuts.AddCutCalo("11310e03","411790105ke30220000","0s331031000000d0"); // 10-30%
  } else if (trainConfig == 2351){ // semicent, no OOB Pileup correction for MC
    cuts.AddCutCalo("11310053","411790105ke30220000","0s331031000000d0"); // 10-30%
  } else if (trainConfig == 2352){ // semicent, no OOB Pileup correction for MC with added signal
    cuts.AddCutCalo("11310023","411790105ke30220000","0s331031000000d0"); // 10-30%
  } else if (trainConfig == 2353){ // semicent, no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutCalo("11310023","411790105ke30220000","0s331031000000d0"); // 10-30%
  } else if (trainConfig == 2354){ // semicent, no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutCalo("11310023","411790105ke30220000","0s331031000000d0"); // 10-30%
  } else if (trainConfig == 2355){ // semicent with NOC NonLin like data
    cuts.AddCutCalo("11310e03","411790105we30220000","0s331031000000d0"); // 10-30%
  } else if (trainConfig == 2356){ // semicent, with NOC NonLin like MC no OOB Pileup correction for MC
    cuts.AddCutCalo("11310053","411790105xe30220000","0s331031000000d0"); // 10-30%
  } else if (trainConfig == 2357){ // semicent, with NOC NonLin like MC no OOB Pileup correction for MC with added signal
    cuts.AddCutCalo("11310023","411790105xe30220000","0s331031000000d0"); // 10-30%
  } else if (trainConfig == 2358){ // semicent, with NOC NonLin like MC no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutCalo("11310023","411790105xe30220000","0s331031000000d0"); // 10-30%
  } else if (trainConfig == 2359){ // semicent, with NOC NonLin like MC no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutCalo("11310023","411790105xe30220000","0s331031000000d0"); // 10-30%
  } else if (trainConfig == 2360){ // semicent with mean N matched tracks per cluster data
    cuts.AddCutCalo("11310e03","411790105te30220000","0s331031000000d0"); // 10-30%
  } else if (trainConfig == 2361){ // semicent, with mean N matched tracks per cluster MC no OOB Pileup correction for MC
    cuts.AddCutCalo("11310053","411790105te30220000","0s331031000000d0"); // 10-30%
  } else if (trainConfig == 2362){ // semicent, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal
    cuts.AddCutCalo("11310023","411790105te30220000","0s331031000000d0"); // 10-30%
  } else if (trainConfig == 2363){ // semicent, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutCalo("11310023","411790105te30220000","0s331031000000d0"); // 10-30%
  } else if (trainConfig == 2364){ // semicent, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutCalo("11310023","411790105te30220000","0s331031000000d0"); // 10-30%
  } else if (trainConfig == 2365){ // semicent NOC via <E_clus-E>/E data
    cuts.AddCutCalo("11310e03","411790105ye30220000","0s331031000000d0"); // 10-30%
  } else if (trainConfig == 2366){ // semicent, NOC via <E_clus-E>/E MC no OOB Pileup correction for MC
    cuts.AddCutCalo("11310053","411790105ye30220000","0s331031000000d0"); // 10-30%
  // **********************************************************************************************************
  // *********************** EMC configurations PbPb run 2 2015 pass 3 semiperipheral *************************
  // **********************************************************************************************************
  // Standard: EMCal + DCal, only 2ndary cluster track matching, mixed event
  } else if (trainConfig == 2500){ // semiperipheral
    cuts.AddCutCalo("13510e03","411790105ke30220000","01331031000000d0"); // 30-50%
  } else if (trainConfig == 2501){ // semiperipheral, no OOB Pileup correction for MC
    cuts.AddCutCalo("13510053","411790105ke30220000","01331031000000d0"); // 30-50%
  } else if (trainConfig == 2502){ // semiperipheral, no OOB Pileup correction for MC with added signal
    cuts.AddCutCalo("13510023","411790105ke30220000","01331031000000d0"); // 30-50%
  } else if (trainConfig == 2503){ // semiperipheral, no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutCalo("13510023","411790105ke30220000","01331031000000d0"); // 30-50%
  } else if (trainConfig == 2504){ // semiperipheral, no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutCalo("13510023","411790105ke30220000","01331031000000d0"); // 30-50%
  } else if (trainConfig == 2505){ // semiperipheral with NOC NonLin like data
    cuts.AddCutCalo("13510e03","411790105we30220000","01331031000000d0"); // 30-50%
  } else if (trainConfig == 2506){ // semiperipheral, with NOC NonLin like MC no OOB Pileup correction for MC
    cuts.AddCutCalo("13510053","411790105xe30220000","01331031000000d0"); // 30-50%
  } else if (trainConfig == 2507){ // semiperipheral, with NOC NonLin like MC no OOB Pileup correction for MC with added signal
    cuts.AddCutCalo("13510023","411790105xe30220000","01331031000000d0"); // 30-50%
  } else if (trainConfig == 2508){ // semiperipheral, with NOC NonLin like MC no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutCalo("13510023","411790105xe30220000","01331031000000d0"); // 30-50%
  } else if (trainConfig == 2509){ // semiperipheral, with NOC NonLin like MC no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutCalo("13510023","411790105xe30220000","01331031000000d0"); // 30-50%
  } else if (trainConfig == 2510){ // semiperipheral with mean N matched tracks per cluster data
    cuts.AddCutCalo("13510e03","411790105te30220000","01331031000000d0"); // 30-50%
  } else if (trainConfig == 2511){ // semiperipheral, with mean N matched tracks per cluster MC no OOB Pileup correction for MC
    cuts.AddCutCalo("13510053","411790105te30220000","01331031000000d0"); // 30-50%
  } else if (trainConfig == 2512){ // semiperipheral, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal
    cuts.AddCutCalo("13510023","411790105te30220000","01331031000000d0"); // 30-50%
  } else if (trainConfig == 2513){ // semiperipheral, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutCalo("13510023","411790105te30220000","01331031000000d0"); // 30-50%
  } else if (trainConfig == 2514){ // semiperipheral, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutCalo("13510023","411790105te30220000","01331031000000d0"); // 30-50%
  } else if (trainConfig == 2515){ // semiperipheral NOC via <E_clus-E>/E data
    cuts.AddCutCalo("13510e03","411790105ye30220000","01331031000000d0"); // 30-50%
  } else if (trainConfig == 2516){ // semiperipheral, NOC via <E_clus-E>/E MC no OOB Pileup correction for MC
    cuts.AddCutCalo("13510053","411790105ye30220000","01331031000000d0"); // 30-50%
  // with rotation background
  } else if (trainConfig == 2550){ // semiperipheral
    cuts.AddCutCalo("13510e03","411790105ke30220000","0s331031000000d0"); // 30-50%
  } else if (trainConfig == 2551){ // semiperipheral, no OOB Pileup correction for MC
    cuts.AddCutCalo("13510053","411790105ke30220000","0s331031000000d0"); // 30-50%
  } else if (trainConfig == 2552){ // semiperipheral, no OOB Pileup correction for MC with added signal
    cuts.AddCutCalo("13510023","411790105ke30220000","0s331031000000d0"); // 30-50%
  } else if (trainConfig == 2553){ // semiperipheral, no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutCalo("13510023","411790105ke30220000","0s331031000000d0"); // 30-50%
  } else if (trainConfig == 2554){ // semiperipheral, no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutCalo("13510023","411790105ke30220000","0s331031000000d0"); // 30-50%
  } else if (trainConfig == 2555){ // semiperipheral with NOC NonLin like data
    cuts.AddCutCalo("13510e03","411790105we30220000","0s331031000000d0"); // 30-50%
  } else if (trainConfig == 2556){ // semiperipheral, with NOC NonLin like MC no OOB Pileup correction for MC
    cuts.AddCutCalo("13510053","411790105xe30220000","0s331031000000d0"); // 30-50%
  } else if (trainConfig == 2557){ // semiperipheral, with NOC NonLin like MC no OOB Pileup correction for MC with added signal
    cuts.AddCutCalo("13510023","411790105xe30220000","0s331031000000d0"); // 30-50%
  } else if (trainConfig == 2558){ // semiperipheral, with NOC NonLin like MC no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutCalo("13510023","411790105xe30220000","0s331031000000d0"); // 30-50%
  } else if (trainConfig == 2559){ // semiperipheral, with NOC NonLin like MC no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutCalo("13510023","411790105xe30220000","0s331031000000d0"); // 30-50%
  } else if (trainConfig == 2560){ // semiperipheral with mean N matched tracks per cluster data
    cuts.AddCutCalo("13510e03","411790105te30220000","0s331031000000d0"); // 30-50%
  } else if (trainConfig == 2561){ // semiperipheral, with mean N matched tracks per cluster MC no OOB Pileup correction for MC
    cuts.AddCutCalo("13510053","411790105te30220000","0s331031000000d0"); // 30-50%
  } else if (trainConfig == 2562){ // semiperipheral, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal
    cuts.AddCutCalo("13510023","411790105te30220000","0s331031000000d0"); // 30-50%
  } else if (trainConfig == 2563){ // semiperipheral, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutCalo("13510023","411790105te30220000","0s331031000000d0"); // 30-50%
  } else if (trainConfig == 2564){ // semiperipheral, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutCalo("13510023","411790105te30220000","0s331031000000d0"); // 30-50%
  } else if (trainConfig == 2565){ // semiperipheral NOC via <E_clus-E>/E data
    cuts.AddCutCalo("13510e03","411790105ye30220000","0s331031000000d0"); // 30-50%
  } else if (trainConfig == 2566){ // semiperipheral, NOC via <E_clus-E>/E MC no OOB Pileup correction for MC
    cuts.AddCutCalo("13510053","411790105ye30220000","0s331031000000d0"); // 30-50%
  // **********************************************************************************************************
  // ************************* EMC configurations PbPb run 2 2015 pass 3 peripheral ***************************
  // **********************************************************************************************************
  // Standard: EMCal + DCal, only 2ndary cluster track matching, mixed event
  } else if (trainConfig == 2900){ // peripheral
    cuts.AddCutCalo("15910e03","411790105ke30220000","01331031000000d0"); // 50-90%
  } else if (trainConfig == 2901){ // peripheral, no OOB Pileup correction for MC
    cuts.AddCutCalo("15910053","411790105ke30220000","01331031000000d0"); // 50-90%
  } else if (trainConfig == 2902){ // peripheral, no OOB Pileup correction for MC with added signal
    cuts.AddCutCalo("15910023","411790105ke30220000","01331031000000d0"); // 50-90%
  } else if (trainConfig == 2903){ // peripheral, no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutCalo("15910023","411790105ke30220000","01331031000000d0"); // 50-90%
  } else if (trainConfig == 2904){ // peripheral, no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutCalo("15910023","411790105ke30220000","01331031000000d0"); // 50-90%
  } else if (trainConfig == 2905){ // peripheral with NOC NonLin like data
    cuts.AddCutCalo("15910e03","411790105we30220000","01331031000000d0"); // 50-90%
  } else if (trainConfig == 2906){ // peripheral, with NOC NonLin like MC no OOB Pileup correction for MC
    cuts.AddCutCalo("15910053","411790105xe30220000","01331031000000d0"); // 50-90%
  } else if (trainConfig == 2907){ // peripheral, with NOC NonLin like MC no OOB Pileup correction for MC with added signal
    cuts.AddCutCalo("15910023","411790105xe30220000","01331031000000d0"); // 50-90%
  } else if (trainConfig == 2908){ // peripheral, with NOC NonLin like MC no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutCalo("15910023","411790105xe30220000","01331031000000d0"); // 50-90%
  } else if (trainConfig == 2909){ // peripheral, with NOC NonLin like MC no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutCalo("15910023","411790105xe30220000","01331031000000d0"); // 50-90%
  } else if (trainConfig == 2910){ // peripheral with mean N matched tracks per cluster data
    cuts.AddCutCalo("15910e03","411790105te30220000","01331031000000d0"); // 50-90%
  } else if (trainConfig == 2911){ // peripheral, with mean N matched tracks per cluster MC no OOB Pileup correction for MC
    cuts.AddCutCalo("15910053","411790105te30220000","01331031000000d0"); // 50-90%
  } else if (trainConfig == 2912){ // peripheral, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal
    cuts.AddCutCalo("15910023","411790105te30220000","01331031000000d0"); // 50-90%
  } else if (trainConfig == 2913){ // peripheral, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutCalo("15910023","411790105te30220000","01331031000000d0"); // 50-90%
  } else if (trainConfig == 2914){ // peripheral, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutCalo("15910023","411790105te30220000","01331031000000d0"); // 50-90%
  } else if (trainConfig == 2915){ // peripheral NOC via <E_clus-E>/E data
    cuts.AddCutCalo("15910e03","411790105ye30220000","01331031000000d0"); // 50-90%
  } else if (trainConfig == 2916){ // peripheral, NOC via <E_clus-E>/E MC no OOB Pileup correction for MC
    cuts.AddCutCalo("15910053","411790105ye30220000","01331031000000d0"); // 50-90%
  // with rotation background
  } else if (trainConfig == 2950){ // peripheral
    cuts.AddCutCalo("15910e03","411790105ke30220000","0s331031000000d0"); // 50-90%
  } else if (trainConfig == 2951){ // peripheral, no OOB Pileup correction for MC
    cuts.AddCutCalo("15910053","411790105ke30220000","0s331031000000d0"); // 50-90%
  } else if (trainConfig == 2952){ // peripheral, no OOB Pileup correction for MC with added signal
    cuts.AddCutCalo("15910023","411790105ke30220000","0s331031000000d0"); // 50-90%
  } else if (trainConfig == 2953){ // peripheral, no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutCalo("15910023","411790105ke30220000","0s331031000000d0"); // 50-90%
  } else if (trainConfig == 2954){ // peripheral, no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutCalo("15910023","411790105ke30220000","0s331031000000d0"); // 50-90%
  } else if (trainConfig == 2955){ // peripheral with NOC NonLin like data
    cuts.AddCutCalo("15910e03","411790105we30220000","0s331031000000d0"); // 50-90%
  } else if (trainConfig == 2956){ // peripheral, with NOC NonLin like MC no OOB Pileup correction for MC
    cuts.AddCutCalo("15910053","411790105xe30220000","0s331031000000d0"); // 50-90%
  } else if (trainConfig == 2957){ // peripheral, with NOC NonLin like MC no OOB Pileup correction for MC with added signal
    cuts.AddCutCalo("15910023","411790105xe30220000","0s331031000000d0"); // 50-90%
  } else if (trainConfig == 2958){ // peripheral, with NOC NonLin like MC no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutCalo("15910023","411790105xe30220000","0s331031000000d0"); // 50-90%
  } else if (trainConfig == 2959){ // peripheral, with NOC NonLin like MC no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutCalo("15910023","411790105xe30220000","0s331031000000d0"); // 50-90%
  } else if (trainConfig == 2960){ // peripheral with mean N matched tracks per cluster data
    cuts.AddCutCalo("15910e03","411790105te30220000","0s331031000000d0"); // 50-90%
  } else if (trainConfig == 2961){ // peripheral, with mean N matched tracks per cluster MC no OOB Pileup correction for MC
    cuts.AddCutCalo("15910053","411790105te30220000","0s331031000000d0"); // 50-90%
  } else if (trainConfig == 2962){ // peripheral, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal
    cuts.AddCutCalo("15910023","411790105te30220000","0s331031000000d0"); // 50-90%
  } else if (trainConfig == 2963){ // peripheral, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal (copy for eta)
    cuts.AddCutCalo("15910023","411790105te30220000","0s331031000000d0"); // 50-90%
  } else if (trainConfig == 2964){ // peripheral, with mean N matched tracks per cluster MC no OOB Pileup correction for MC with added signal (copy for pi0 + eta)
    cuts.AddCutCalo("15910023","411790105te30220000","0s331031000000d0"); // 50-90%
  } else if (trainConfig == 2965){ // peripheral NOC via <E_clus-E>/E data
    cuts.AddCutCalo("15910e03","411790105ye30220000","0s331031000000d0"); // 50-90%
  } else if (trainConfig == 2966){ // peripheral, NOC via <E_clus-E>/E MC no OOB Pileup correction for MC
    cuts.AddCutCalo("15910053","411790105ye30220000","0s331031000000d0"); // 50-90%

  } else {
    Error(Form("GammaConvCalo_%i",trainConfig), "wrong trainConfig variable no cuts have been specified for the configuration");
    return;
  }

  if(!cuts.AreValid()){
    cout << "\n\n****************************************************" << endl;
    cout << "ERROR: No valid cuts stored in CutHandlerCalo! Returning..." << endl;
    cout << "****************************************************\n\n" << endl;
    return;
  }

  Int_t numberOfCuts = cuts.GetNCuts();

  TList *HeaderList = new TList();
  if (generatorName.CompareTo("LHC13d2")==0){
    TObjString *Header1 = new TObjString("pi0_1");
    HeaderList->Add(Header1);
  //    TObjString *Header3 = new TObjString("eta_2");
  //    HeaderList->Add(Header3);

  } else if (generatorName.CompareTo("LHC12a17x_fix")==0){
    TObjString *Header1 = new TObjString("PARAM");
    HeaderList->Add(Header1);
  } else if (generatorName.CompareTo("LHC14a1a")==0){
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
  } else if (generatorName.CompareTo("LHC14a1b")==0 || generatorName.CompareTo("LHC14a1c")==0){
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
  } else if (generatorName.CompareTo("LHC16h4")==0 || generatorName.CompareTo("LHC19h3")==0){
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
  } else if ( (generatorName.CompareTo("LHC20g10")==0) || generatorName.BeginsWith("LHC24a1") ) {

    auto fillSingle = [&HeaderList](Size_t theSwitch){
      if (theSwitch == 1 ) { HeaderList->Add(new TObjString("Injector (pi0)"));  } // PCM Pi0
      if (theSwitch == 2 ) { HeaderList->Add(new TObjString("Injector (pi0a)")); } // PHOS Pi0
      if (theSwitch == 3 ) { HeaderList->Add(new TObjString("Injector (pi0b)")); } // PHOS Pi0
      if (theSwitch == 4 ) { HeaderList->Add(new TObjString("Injector (pi0c)")); } // PHOS Pi0
      if (theSwitch == 5 ) { HeaderList->Add(new TObjString("Injector (pi0d)")); } // PHOS Pi0
      if (theSwitch == 6 ) { HeaderList->Add(new TObjString("Injector (eta)"));  } // PCM Eta
      if (theSwitch == 7 ) { HeaderList->Add(new TObjString("Injector (etaa)")); } // PHOS Eta
      if (theSwitch == 8 ) { HeaderList->Add(new TObjString("Hijing")); }          // HIJING
      if (theSwitch == 9 ) { HeaderList->Add(new TObjString("Pileup")); }          // PileUp
      if (theSwitch == 10 ) { HeaderList->Add(new TObjString("Injector (pi0e)")); } // EMCAL 1 Pi0
      if (theSwitch == 11 ) { HeaderList->Add(new TObjString("Injector (pi0f)")); } // EMCAL 2 Pi0
      if (theSwitch == 12 ) { HeaderList->Add(new TObjString("Injector (pi0g)")); } // EMCAL 3 Pi0
      if (theSwitch == 13 ) { HeaderList->Add(new TObjString("Injector (pi0h)")); } // DCAL 1 Pi0
      if (theSwitch == 14 ) { HeaderList->Add(new TObjString("Injector (pi0i)")); } // DCal 2 Pi0
      if (theSwitch == 15 ) { HeaderList->Add(new TObjString("Injector (etab)")); } // EMCAL 1 Eta
      if (theSwitch == 16 ) { HeaderList->Add(new TObjString("Injector (etac)")); } // EMCAL 2 Eta
      if (theSwitch == 17 ) { HeaderList->Add(new TObjString("Injector (etad)")); } // EMCAL 3 Eta
      if (theSwitch == 18 ) { HeaderList->Add(new TObjString("Injector (etae)")); } // DCAL 1 Eta
      if (theSwitch == 19 ) { HeaderList->Add(new TObjString("Injector (etaf)")); } // DCal 2 Eta
    };

    auto fillFromTo = [&fillSingle](Size_t theStart, Size_t theInclEnd){
      for (Size_t i=theStart; i<=theInclEnd; ++i) { fillSingle(i); }
    };

    if      (headerSelectionInt == 0)  { fillFromTo(1, 9); }            // all PCM + PHOS + HIJING + PileUp
    else if (headerSelectionInt <= 19)  { fillSingle(headerSelectionInt); } // a single one
    else if (headerSelectionInt == 20) { fillFromTo(1, 5); }            // pi0) + pi0x)
    else if (headerSelectionInt == 21) { fillFromTo(2, 5); }            // all pi0x)
    else if (headerSelectionInt == 22) { fillFromTo(6, 7); }            // eta) + etaa)
    else if (headerSelectionInt == 23) {                                // pi0) + eta)
      fillSingle(1);
      fillSingle(6);
    } 
    else if (headerSelectionInt == 24) { fillFromTo(10, 19); }            // all EMCal
    else if (headerSelectionInt == 25) { fillFromTo(1, 19); }             // all
    else if (headerSelectionInt == 26) {                                  // PCM + all EMCal Pi0 + Eta
      fillSingle(1);
      fillSingle(6);
      fillFromTo(10, 19);
    }
    else if (headerSelectionInt == 27) {                                  // PCM + all EMCal + HIJING Pi0 + Eta
      fillSingle(1);
      fillSingle(6);
      fillSingle(8);
      fillFromTo(10, 19);
    }
    else if (headerSelectionInt == 28) {                                  // PCM + PileUp Pi0 + Eta
      fillSingle(1);
      fillSingle(6);
      fillSingle(8);
      fillSingle(9);
    }
    else if (headerSelectionInt == 29) {                                  // PCM + all EMCal + HIJING + PileUp Pi0 + Eta
      fillSingle(1);
      fillSingle(6);
      fillFromTo(8, 19);
    }
    else if (headerSelectionInt == 30) {                                  // all EMCal Pi0
      fillFromTo(10, 14);
    }
    else if (headerSelectionInt == 31) {                                  // all EMCal + HIJING Pi0
      fillSingle(8);
      fillFromTo(10, 14);
    }
    else if (headerSelectionInt == 32) {                                  // all EMCal + HIJING + PileUp Pi0
      fillFromTo(8, 14);
    }
    else if (headerSelectionInt == 33) {                                  // PCM + all EMCal Pi0
      fillSingle(1);
      fillFromTo(10, 14);
    }
    else if (headerSelectionInt == 34) {                                  // PCM + all EMCal + HIJING Pi0
      fillSingle(1);
      fillSingle(8);
      fillFromTo(10, 14);
    }
    else if (headerSelectionInt == 35) {                                  // PCM + all EMCal + HIJING + PileUp Pi0
      fillSingle(1);
      fillFromTo(8, 14);
    }
    else if (headerSelectionInt == 36) {                                  // PCM + HIJING + PileUp Pi0
      fillSingle(1);
      fillSingle(8);
      fillSingle(9);
    }
    else if (headerSelectionInt == 37) {                                  // PCM + HIJING Pi0
      fillSingle(1);
      fillSingle(8);
      fillSingle(9);
    }
    else if (headerSelectionInt == 38) {                                  // all EMCal Eta
      fillFromTo(15, 19);
    }
    else if (headerSelectionInt == 39) {                                  // all EMCal + HIJING Eta
      fillSingle(8);
      fillFromTo(15, 19);
    }
    else if (headerSelectionInt == 40) {                                  // all EMCal + HIJING + PileUp Eta
      fillSingle(8);
      fillSingle(9);
      fillFromTo(15, 14);
    }
    else if (headerSelectionInt == 41) {                                  // PCM + all EMCal Eta
      fillSingle(6);
      fillFromTo(15, 19);
    }
    else if (headerSelectionInt == 42) {                                  // PCM + all EMCal + HIJING Eta
      fillSingle(6);
      fillSingle(8);
      fillFromTo(15, 19);
    }
    else if (headerSelectionInt == 43) {                                  // PCM + all EMCal + HIJING + PileUp Eta
      fillSingle(6);
      fillSingle(8);
      fillSingle(9);
      fillFromTo(15, 19);
    }
    else if (headerSelectionInt == 44) {                                  // PCM + HIJING + PileUp Eta
      fillSingle(6);
      fillSingle(8);
      fillSingle(9);
    }
    else if (headerSelectionInt == 45) {                                  // PCM + HIJING Eta
      fillSingle(6);
      fillSingle(8);
    }
    else {
      Error(Form("%s_%i", addTaskName.Data(),  trainConfig), "No valid header selection");
      return;
    }
 }

  TList *EventCutList   = new TList();
  TList *ClusterCutList = new TList();
  TList *MesonCutList   = new TList();


  EventCutList->SetOwner(kTRUE);
  AliConvEventCuts **analysisEventCuts        = new AliConvEventCuts*[numberOfCuts];
  ClusterCutList->SetOwner(kTRUE);
  AliCaloPhotonCuts **analysisClusterCuts     = new AliCaloPhotonCuts*[numberOfCuts];
  MesonCutList->SetOwner(kTRUE);
  AliConversionMesonCuts **analysisMesonCuts  = new AliConversionMesonCuts*[numberOfCuts];

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
      if (enableLightOutput == 4) fTrackMatcher->SetLightOutput(kTRUE);
      mgr->AddTask(fTrackMatcher);
      mgr->ConnectInput(fTrackMatcher,0,cinput);
    }

    analysisEventCuts[i]    = new AliConvEventCuts();
    analysisEventCuts[i]->SetV0ReaderName(V0ReaderName);
    if (periodNameV0Reader.CompareTo("") != 0) analysisEventCuts[i]->SetPeriodEnum(periodNameV0Reader);

    if(generatorName.Contains("LHC11h") && doFlattening){
      cout << "entering the cent. flattening loop -> searching for file: " << fileNameCentFlattening.Data() << endl;

      if( fileNameCentFlattening.Contains("FlatFile") ){
        analysisEventCuts[i]->SetUseWeightFlatCentralityFromFile(doFlattening, fileNameCentFlattening, "Cent");
      } else if( fileNameCentFlattening.Contains("Good") ){
        analysisEventCuts[i]->SetUseWeightFlatCentralityFromFile(doFlattening, fileNameCentFlattening, "CentGoodRuns");
      }else if( fileNameCentFlattening.Contains("SemiGood") ){
        analysisEventCuts[i]->SetUseWeightFlatCentralityFromFile(doFlattening, fileNameCentFlattening, "CentSemiGoodRuns");
      }else {
        analysisEventCuts[i]->SetUseWeightFlatCentralityFromFile(doFlattening, fileNameCentFlattening, "CentTotalRuns");
      }
    }

    TString histoNameMCPi0PT = "";  TString histoNameMCEtaPT = "";  TString histoNameMCK0sPT = "";    // pT spectra to be weighted
    TString fitNamePi0PT = "";      TString fitNameEtaPT = "";      TString fitNameK0sPT = "";        // fit to correct shape of pT spectra
    Bool_t weightPi0 = kFALSE;      Bool_t weightEta = kFALSE;      Bool_t weightK0s = kFALSE;
    if(doWeighting){
      cout << "INFO enabeling pT weighting" << endl;
      if(periodNameAnchor.CompareTo("LHC18q")==0){
        TString eventCutString  = cuts.GetEventCut(i);
        TString eventCutShort   = eventCutString(0,6);   // first six digits
        weightPi0         = kTRUE;
        histoNameMCPi0PT  = Form("Pi0_%s_5TeV_%s",   periodNameV0Reader.Data(), eventCutString.Data());  // MC
        fitNamePi0PT      = Form("Pi0_Data_5TeV_%s", eventCutShort.Data());                              // fit to data
        weightEta         = kTRUE;
        histoNameMCEtaPT  = Form("Eta_%s_5TeV_%s",   periodNameV0Reader.Data(), eventCutString.Data());
        fitNameEtaPT      = Form("Eta_Data_5TeV_%s", eventCutShort.Data());
      }
      analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(weightPi0, weightEta, weightK0s, fileNamePtWeights, histoNameMCPi0PT, histoNameMCEtaPT, histoNameMCK0sPT, fitNamePi0PT, fitNameEtaPT, fitNameK0sPT);
    }

    TString dataInputMultHisto  = "";
    TString mcInputMultHisto    = "";
    if (enableMultiplicityWeighting){
      cout << "INFO enableling mult weighting" << endl;
      if(periodNameAnchor.CompareTo("LHC15o")==0){
        TString cutNumber = cuts.GetEventCut(i);
        TString centCut = cutNumber(0,3);  // first three digits of event cut
        dataInputMultHisto = Form("%s_%s", periodNameAnchor.Data(), centCut.Data());
        mcInputMultHisto   = Form("%s_%s", generatorName.Data(), centCut.Data());
        cout << "INFO read " << dataInputMultHisto.Data() << " and " <<  mcInputMultHisto.Data() << " from " << fileNameMultWeights.Data() << endl;
      }
      analysisEventCuts[i]->SetUseWeightMultiplicityFromFile(kTRUE, fileNameMultWeights, dataInputMultHisto, mcInputMultHisto );
    }

    if (trainConfig == 34 || trainConfig == 35 ){
      if (generatorName.CompareTo("LHC14a1a") ==0 || generatorName.CompareTo("LHC14a1b") ==0 || generatorName.CompareTo("LHC14a1c") ==0 ){
        if ( i == 0 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNamePtWeights, Form("Pi0_Hijing_%s_PbPb_2760GeV_0010TPC",generatorName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_0010TPC",generatorName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0010V0M","Eta_Fit_Data_PbPb_2760GeV_0010V0M");
        if ( i == 1 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNamePtWeights, Form("Pi0_Hijing_%s_PbPb_2760GeV_0010TPC",generatorName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_0010TPC",generatorName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0010V0M","Eta_Fit_Data_PbPb_2760GeV_0010V0M");
      }
    }
    if (trainConfig == 36 || trainConfig == 37){
      if (generatorName.CompareTo("LHC14a1a") ==0 || generatorName.CompareTo("LHC14a1b") ==0 || generatorName.CompareTo("LHC14a1c") ==0 ){
        if ( i == 0 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNamePtWeights, Form("Pi0_Hijing_%s_PbPb_2760GeV_2050TPC",generatorName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_2050TPC",generatorName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_2050V0M","Eta_Fit_Data_PbPb_2760GeV_2050V0M");
        if ( i == 1 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNamePtWeights, Form("Pi0_Hijing_%s_PbPb_2760GeV_2050TPC",generatorName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_2050TPC",generatorName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_2050V0M","Eta_Fit_Data_PbPb_2760GeV_2050V0M");
      }
    }

    if (trainConfig == 38 || trainConfig == 39 || trainConfig == 43 || trainConfig == 44 ){
      if (generatorName.CompareTo("LHC14a1a") ==0 || generatorName.CompareTo("LHC14a1b") ==0 || generatorName.CompareTo("LHC14a1c") ==0 ){
        if ( i == 0 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNamePtWeights, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_0010TPC",generatorName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_0010TPC",generatorName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0010V0M","Eta_Fit_Data_PbPb_2760GeV_0010V0M");
        if ( i == 1 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNamePtWeights, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_0010TPC",generatorName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_0010TPC",generatorName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0010V0M","Eta_Fit_Data_PbPb_2760GeV_0010V0M");
        if ( i == 2 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNamePtWeights, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_0010TPC",generatorName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_0010TPC",generatorName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0010V0M","Eta_Fit_Data_PbPb_2760GeV_0010V0M");
      }
    }

    if (trainConfig == 40 || trainConfig == 41){
      if (generatorName.CompareTo("LHC14a1a") ==0 || generatorName.CompareTo("LHC14a1b") ==0 || generatorName.CompareTo("LHC14a1c") ==0 ){
        if ( i == 0 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNamePtWeights, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_2050TPC",generatorName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_2050TPC",generatorName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_2050V0M","Eta_Fit_Data_PbPb_2760GeV_2050V0M");
        if ( i == 1 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNamePtWeights, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_2050TPC",generatorName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_2050TPC",generatorName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_2050V0M","Eta_Fit_Data_PbPb_2760GeV_2050V0M");
        if ( i == 2 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNamePtWeights, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_2050TPC",generatorName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_2050TPC",generatorName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_2050V0M","Eta_Fit_Data_PbPb_2760GeV_2050V0M");
      }
    }

    analysisEventCuts[i]->SetDebugLevel(0);
    if (enableLightOutput > 0) analysisEventCuts[i]->SetLightOutput(1);
    if(fMinPtHardSet)
      analysisEventCuts[i]->SetMinFacPtHard(minFacPtHard);
    if(fMaxPtHardSet)
      analysisEventCuts[i]->SetMaxFacPtHard(maxFacPtHard);
    if(fSingleMaxPtHardSet)
      analysisEventCuts[i]->SetMaxFacPtHardSingleParticle(maxFacPtHardSingle);
    if(fJetFinderUsage)
      analysisEventCuts[i]->SetUseJetFinderForOutliers(kTRUE);
    if(fUsePtHardFromFile)
      analysisEventCuts[i]->SetUsePtHardBinFromFile(kTRUE);
    if(fUseAddOutlierRej)
      analysisEventCuts[i]->SetUseAdditionalOutlierRejection(kTRUE);
    analysisEventCuts[i]->InitializeCutsFromCutString((cuts.GetEventCut(i)).Data());
    EventCutList->Add(analysisEventCuts[i]);
    analysisEventCuts[i]->SetCorrectionTaskSetting(corrTaskSetting);
    analysisEventCuts[i]->SetFillCutHistograms("",kFALSE);
    if (generatorName.CompareTo("LHC14a1b") ==0 || generatorName.CompareTo("LHC14a1c") ==0 ){
      if (headerSelectionInt == 1 || headerSelectionInt == 4 || headerSelectionInt == 12 ) analysisEventCuts[i]->SetAddedSignalPDGCode(111);
      if (headerSelectionInt == 2 || headerSelectionInt == 5 || headerSelectionInt == 13 ) analysisEventCuts[i]->SetAddedSignalPDGCode(221);
    }
    if (enableLightOutput == 1 || enableLightOutput == 2 ) analysisEventCuts[i]->SetLightOutput(1);
    if (enableLightOutput == 4) analysisEventCuts[i]->SetLightOutput(2);

    analysisClusterCuts[i]  = new AliCaloPhotonCuts(isMC);
    analysisClusterCuts[i]->SetHistoToModifyAcceptance(histoAcc);
    analysisClusterCuts[i]->SetV0ReaderName(V0ReaderName);
    analysisClusterCuts[i]->SetCorrectionTaskSetting(corrTaskSetting);
    analysisClusterCuts[i]->SetCaloTrackMatcherName(TrackMatcherName);
    if (enableLightOutput == 1 || enableLightOutput == 2 || enableLightOutput ==5 ) analysisClusterCuts[i]->SetLightOutput(1);
    if (enableLightOutput == 4) analysisClusterCuts[i]->SetLightOutput(2);
    analysisClusterCuts[i]->InitializeCutsFromCutString((cuts.GetClusterCut(i)).Data());
    ClusterCutList->Add(analysisClusterCuts[i]);
    analysisClusterCuts[i]->SetExtendedMatchAndQA(enableExtMatchAndQA);

    analysisMesonCuts[i]    = new AliConversionMesonCuts();
    if (enableLightOutput > 0 && enableLightOutput != 4) analysisMesonCuts[i]->SetLightOutput(1);
    if (enableLightOutput == 4) analysisMesonCuts[i]->SetLightOutput(2);
    analysisMesonCuts[i]->InitializeCutsFromCutString((cuts.GetMesonCut(i)).Data());
    analysisMesonCuts[i]->SetIsMergedClusterCut(2);
    analysisMesonCuts[i]->SetCaloMesonCutsObject(analysisClusterCuts[i]);
    MesonCutList->Add(analysisMesonCuts[i]);
    analysisMesonCuts[i]->SetFillCutHistograms("");
    analysisEventCuts[i]->SetAcceptedHeader(HeaderList);

    if(analysisMesonCuts[i]->DoGammaSwappForBg()) analysisClusterCuts[i]->SetUseEtaPhiMapForBackCand(kTRUE);
    analysisClusterCuts[i]->SetFillCutHistograms("");
  }
  task->SetAllowOverlapHeaders(enableHeaderOverlap);
  task->SetEventCutList(numberOfCuts,EventCutList);
  task->SetCaloCutList(numberOfCuts,ClusterCutList);
  task->SetMesonCutList(numberOfCuts,MesonCutList);
  task->SetDoMesonAnalysis(kTRUE);
  task->SetCorrectionTaskSetting(corrTaskSetting);
  task->SetDoMesonQA(enableQAMesonTask);        //Attention new switch for Pi0 QA
  task->SetDoClusterQA(enableQAClusterTask);    //Attention new switch small for Cluster QA
  task->SetDoTHnSparse(enableTHnSparse);
  task->SetProduceTreeEOverP(doTreeEOverP);
  task->SetEnableSortingOfMCClusLabels(enableSortingMCLabels);
  if(enableExtMatchAndQA > 1){ task->SetPlotHistsExtQA(kTRUE);}

  //connect containers
  AliAnalysisDataContainer *coutput =
    mgr->CreateContainer(!(corrTaskSetting.CompareTo("")) ? Form("GammaCalo_%i",trainConfig) : Form("GammaCalo_%i_%s",trainConfig,corrTaskSetting.Data()), TList::Class(),
              AliAnalysisManager::kOutputContainer,Form("GCa_%i.root",trainConfig) );

  mgr->AddTask(task);
  mgr->ConnectInput(task,0,cinput);
  mgr->ConnectOutput(task,1,coutput);

  Int_t nContainer = 2;
  for(Int_t i = 0; i<numberOfCuts; i++){
      if(enableQAMesonTask==5){
          mgr->ConnectOutput(task,nContainer,mgr->CreateContainer(Form("%s_%s_%s ClusterTimingEff",(cuts.GetEventCut(i)).Data(),(cuts.GetClusterCut(i)).Data(),(cuts.GetMesonCut(i)).Data()), TTree::Class(), AliAnalysisManager::kOutputContainer, Form("GCa_%i.root",trainConfig)) );
          nContainer++;
      }
  }


  return;
}
