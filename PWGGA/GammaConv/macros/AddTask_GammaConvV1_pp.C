/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: Friederike Bock, Lucia Leardini                                *
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
//($ALIPHYSICS/PWGGA/GammaConv/AliAnalysisTaskGammaConvV1.cxx) for
//pp together with all supporting classes
//***************************************************************************************

void AddTask_GammaConvV1_pp(
    Int_t     trainConfig                   = 1,        // change different set of cuts
    Int_t     isMC                          = 0,        // run MC
    TString   photonCutNumberV0Reader       = "",       // 00000008400000000100000000 nom. B, 00000088400000000100000000 low B
    TString   periodNameV0Reader            = "",
    // general setting for task
    Int_t     enableQAMesonTask             = 0,        // enable QA in AliAnalysisTaskGammaConvV1
    Int_t     enableQAPhotonTask            = 0,        // enable additional QA task
    Bool_t    enableLightOutput             = kFALSE,   // switch to run light output (only essential histograms for afterburner)
    Bool_t    enableTHnSparse               = kFALSE,   // switch on THNsparse
    Int_t     enableTriggerMimicking        = 0,        // enable trigger mimicking
    Bool_t    enableTriggerOverlapRej       = kFALSE,   // enable trigger overlap rejection
    TString   settingMaxFacPtHard           = "3.",       // maximum factor between hardest jet and ptHard generated
    Int_t     debugLevel                    = 0,        // introducing debug levels for grid running
    // settings for weights
    // FPTW:fileNamePtWeights, FMUW:fileNameMultWeights, FMAW:fileNameMatBudWeights, FEPC:fileNamedEdxPostCalib, separate with ;
    TString   fileNameExternalInputs        = "",
    Int_t     doWeightingPart               = 0,        // enable Weighting
    TString   generatorName                 = "DPMJET", // generator Name
    Bool_t    enableMultiplicityWeighting   = kFALSE,   //
    TString   periodNameAnchor              = "",       //
    Int_t     enableMatBudWeightsPi0        = 0,        // 1 = three radial bins, 2 = 10 radial bins
    Bool_t    enableElecDeDxPostCalibration = kFALSE,
    // special settings
    Bool_t    enableChargedPrimary          = kFALSE,
    Bool_t    doSmear                       = kFALSE,   // switches to run user defined smearing
    Double_t  bremSmear                     = 1.,
    Double_t  smearPar                      = 0.,       // conv photon smearing params
    Double_t  smearParConst                 = 0.,       // conv photon smearing params
    // subwagon config
    TString   additionalTrainConfig         = "0"       // additional counter for trainconfig + special settings
                            ) {

  Int_t trackMatcherRunningMode = 0; // CaloTrackMatcher running mode
  Int_t isHeavyIon = 0;

  AliCutHandlerPCM cuts;

  TString fileNamePtWeights           = cuts.GetSpecialFileNameFromString (fileNameExternalInputs, "FPTW:");
  TString fileNameMultWeights         = cuts.GetSpecialFileNameFromString (fileNameExternalInputs, "FMUW:");
  TString fileNameMatBudWeights       = cuts.GetSpecialFileNameFromString (fileNameExternalInputs, "FMAW:");
  TString fileNamedEdxPostCalib       = cuts.GetSpecialFileNameFromString (fileNameExternalInputs, "FEPC:");

  TString addTaskName                 = "AddTask_GammaConvV1_pp";
  TString sAdditionalTrainConfig      = cuts.GetSpecialSettingFromAddConfig(additionalTrainConfig, "", "", addTaskName);
  if (sAdditionalTrainConfig.Atoi() > 0){
    trainConfig = trainConfig + sAdditionalTrainConfig.Atoi();
    cout << "INFO: " << addTaskName.Data() << " running additionalTrainConfig '" << sAdditionalTrainConfig.Atoi() << "', train config: '" << trainConfig << "'" << endl;
  }
  TString corrTaskSetting         = cuts.GetSpecialSettingFromAddConfig(additionalTrainConfig, "CF","", addTaskName);
  if(corrTaskSetting.CompareTo(""))
    cout << "corrTaskSetting: " << corrTaskSetting.Data() << endl;
  if(additionalTrainConfig.Contains("MaterialBudgetWeights"))
    fileNameMatBudWeights         = cuts.GetSpecialSettingFromAddConfig(additionalTrainConfig, "MaterialBudgetWeights",fileNameMatBudWeights, addTaskName);
  TString strTrackMatcherRunningMode         = cuts.GetSpecialSettingFromAddConfig(additionalTrainConfig, "TM","", addTaskName);
  if(additionalTrainConfig.Contains("TM"))
    trackMatcherRunningMode = strTrackMatcherRunningMode.Atoi();


 TObjArray *rmaxFacPtHardSetting = settingMaxFacPtHard.Tokenize("_");
  if(rmaxFacPtHardSetting->GetEntries()<1){cout << "ERROR: AddTask_GammaConvV1_pp during parsing of settingMaxFacPtHard String '" << settingMaxFacPtHard.Data() << "'" << endl; return;}
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

  // ================== GetAnalysisManager ===============================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error(Form("%s_%i", addTaskName.Data(),  trainConfig), "No analysis manager found.");
    return ;
  }

  // ================== GetInputEventHandler =============================
  AliVEventHandler *inputHandler=mgr->GetInputEventHandler();

  //========= Check whether PID Reponse is there ====
  if(!(AliPIDResponse*)mgr->GetTask("PIDResponseTask")){
    Error(Form("AddTask_GammaConvV1_%i",trainConfig), "No PID response has been initialized aborting.");
    return;
  }

  //=========  Set Cutnumber for V0Reader ================================
  TString cutnumberPhoton     = photonCutNumberV0Reader.Data();
  TString cutnumberEvent      = "00000003";
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();

  //========= Check V0 Reader in  ANALYSIS manager  =====
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
  //            find input container
  AliAnalysisTaskGammaConvV1 *task=NULL;
  task= new AliAnalysisTaskGammaConvV1(Form("GammaConvV1_%i",trainConfig));
  task->SetIsHeavyIon(isHeavyIon);
  task->SetIsMC(isMC);
  task->SetV0ReaderName(V0ReaderName);
  task->SetLightOutput(enableLightOutput);
  // Cut Numbers to use in Analysis


  //---------  standard configurations for 2.76TeV V00R without SDD --------------------------
  if(trainConfig == 1){    // various standard cuts
    cuts.AddCutPCM("00000113", "00200009397300008250400000", "0163103100900000"); // new pi0/eta cut 2.76TeV
    cuts.AddCutPCM("00000113", "00200009397300008250400000", "0163103100000000"); // new pi0/eta cut 2.76TeV without MC smearing
    cuts.AddCutPCM("00000113", "00200009366300003800000000", "0163103100900000"); // standard cut Pi0 pp 2.76TeV PbPb paper 2012
  } else if (trainConfig == 2) { // various standard cuts added signals
    cuts.AddCutPCM("00000123", "00200009397300008250400000", "0163103100900000"); // new pi0/eta cut 2.76TeV
    cuts.AddCutPCM("00000123", "00200009397300008250400000", "0163103100000000"); // new pi0/eta cut 2.76TeV without MC smearing
    cuts.AddCutPCM("00000123", "00200009366300003800000000", "0163103100900000"); // standard cut Pi0 pp 2.76TeV PbPb paper 2012
  } else if (trainConfig == 3) { // additional standards
    cuts.AddCutPCM("00000113", "00200009297002008250400000", "0163103100900000"); // standard cut LHC11h pp 2.76TeV
    cuts.AddCutPCM("00000113", "00200009227302008250404000", "0163101500000000"); // Ana eta analysis prefered
    cuts.AddCutPCM("00000113", "00200008366300000200000000", "0163103100900000"); // standard cut Pi0 pp 7TeV, all photon qualities
  } else if (trainConfig == 4) { // additional standards added signals
    cuts.AddCutPCM("00000123", "00200009297002008250400000", "0163103100900000"); // standard cut LHC11h pp 2.76TeV
    cuts.AddCutPCM("00000123", "00200009227302008250404000", "0163101500000000"); // Ana eta analysis prefered
    cuts.AddCutPCM("00000123", "00200008366300000200000000", "0163103100900000"); // standard cut Pi0 pp 7TeV, all photon qualities
  } else if (trainConfig == 5) {
    // variations to different standards
    cuts.AddCutPCM("00000113", "00200009327000008250400000", "0163103100900000"); // go with 1sigma pi rejec to infty
    cuts.AddCutPCM("00000113", "00200009317000008250400000", "0163103100900000"); // go with 0sigma pi rejec to infty
    cuts.AddCutPCM("00000113", "00200009357000008250400000", "0163103100900000"); // go with 2sigma pi reject to infy
    // eta cut variation
    cuts.AddCutPCM("00000113", "03200009397300008250400000", "0163103100900000"); // eta 0.65
    cuts.AddCutPCM("00000113", "04200009397300008250400000", "0163103100900000"); // eta 0.75
  } else if (trainConfig == 6) {
    // variations to different standards added signals
    cuts.AddCutPCM("00000123", "00200009327000008250400000", "0163103100900000"); // go with 1sigma pi rejec to infty
    cuts.AddCutPCM("00000123", "00200009317000008250400000", "0163103100900000"); // go with 0sigma pi rejec to infty
    cuts.AddCutPCM("00000123", "00200009357000008250400000", "0163103100900000"); // go with 2sigma pi reject to infy
    // eta cut variation added signals
    cuts.AddCutPCM("00000123", "03200009397300008250400000", "0163103100900000"); // eta 0.65
    cuts.AddCutPCM("00000123", "04200009397300008250400000", "0163103100900000"); // eta 0.75
  } else if (trainConfig == 7) {  //pion rejection variations
    cuts.AddCutPCM("00000113", "00200009395300008250400000", "0163103100900000"); // min for pi dEdx 0.3GeV
    cuts.AddCutPCM("00000113", "00200009390300008250400000", "0163103100900000"); // min for pi dEdx  0.5GeV
    cuts.AddCutPCM("00000113", "00200009397400008250400000", "0163103100900000"); // max for pi dEdx  3 GeV
    cuts.AddCutPCM("00000113", "00200009397200008250400000", "0163103100900000"); // new pi0/eta cut 4 GeV
  } else if (trainConfig == 8) { //pion rejection variations added signals
    cuts.AddCutPCM("00000123", "00200009395300008250400000", "0163103100900000"); // min for pi dEdx 0.3GeV
    cuts.AddCutPCM("00000123", "00200009390300008250400000", "0163103100900000"); // min for pi dEdx  0.5GeV
    cuts.AddCutPCM("00000123", "00200009397400008250400000", "0163103100900000"); // max for pi dEdx  3 GeV
    cuts.AddCutPCM("00000123", "00200009397200008250400000", "0163103100900000"); // new pi0/eta cut 4 GeV
  } else if (trainConfig == 9) {  //electron rejection variations
    cuts.AddCutPCM("00000113", "00200009197300008250400000", "0163103100900000"); // dEdx e +-5 sigma
    cuts.AddCutPCM("00000113", "00200009497300008250400000", "0163103100900000"); // dEdx e -6,+7 sigma
    cuts.AddCutPCM("00000113", "00200009297300008250400000", "0163103100900000"); // dEdx e -3,+5 sigma
    cuts.AddCutPCM("00000113", "00200009597300008250400000", "0163103100900000"); // dEdx e +-4 sigma
  } else if (trainConfig == 10) { //electron rejection variations added signals
    cuts.AddCutPCM("00000123", "00200009197300008250400000", "0163103100900000"); // dEdx e +-5 sigma
    cuts.AddCutPCM("00000123", "00200009497300008250400000", "0163103100900000"); // dEdx e -6,+7 sigma
    cuts.AddCutPCM("00000123", "00200009297300008250400000", "0163103100900000"); // dEdx e -3,+5 sigma
    cuts.AddCutPCM("00000123", "00200009597300008250400000", "0163103100900000"); // dEdx e +-4 sigma
  } else if (trainConfig == 11) { // single leg cuts
    cuts.AddCutPCM("00000113", "00200049397300008250400000", "0163103100900000"); // variation pt 0.075
    cuts.AddCutPCM("00000113", "00200019397300008250400000", "0163103100900000"); // variation pt 0.1
    cuts.AddCutPCM("00000113", "00200006397300008250400000", "0163103100900000"); // variation TPC cls 0.7
    cuts.AddCutPCM("00000113", "00200008397300008250400000", "0163103100900000"); // variation TPC cls 0.35
  } else if (trainConfig == 12) { // single leg cuts added signals
    cuts.AddCutPCM("00000123", "00200049397300008250400000", "0163103100900000"); // variation pt 0.075
    cuts.AddCutPCM("00000123", "00200019397300008250400000", "0163103100900000"); // variation pt 0.1
    cuts.AddCutPCM("00000123", "00200006397300008250400000", "0163103100900000"); // variation TPC cls 0.7
    cuts.AddCutPCM("00000123", "00200008397300008250400000", "0163103100900000"); // variation TPC cls 0.35
  } else if (trainConfig == 13) { // Qt variations
    cuts.AddCutPCM("00000113", "00200009397300009250400000", "0163103100900000"); // variation qt 0.03
    cuts.AddCutPCM("00000113", "00200009397300002250400000", "0163103100900000"); // variation qt 0.06
    cuts.AddCutPCM("00000113", "00200009397300003250400000", "0163103100900000"); // variation qt 0.05 no 2D
    cuts.AddCutPCM("00000113", "00200009397300008250400000", "0163105100900000"); // tighter alpha meson 0.75
  } else if (trainConfig == 14) { // Qt variations added signals
    cuts.AddCutPCM("00000123", "00200009397300009250400000", "0163103100900000"); // variation qt 0.03
    cuts.AddCutPCM("00000123", "00200009397300002250400000", "0163103100900000"); // variation qt 0.06
    cuts.AddCutPCM("00000123", "00200009397300003250400000", "0163103100900000"); // variation qt 0.05 no 2D
    cuts.AddCutPCM("00000123", "00200009397300008250400000", "0163105100900000"); // tighter alpha meson 0.75
  } else if (trainConfig == 15) { // chi2 - Psi pair variations
    cuts.AddCutPCM("00000113", "00200009397300008150400000", "0163103100900000"); // chi2 50 with psi pair 0.1
    cuts.AddCutPCM("00000113", "00200009397300008850400000", "0163103100900000"); // chi2 20 with psi pair 0.1
    cuts.AddCutPCM("00000113", "00200009397300008280400000", "0163103100900000"); // chi2 30 with psi pair 0.2
    cuts.AddCutPCM("00000113", "00200009397300008260400000", "0163103100900000"); // chi2 30 with psi pair 0.05
  } else if (trainConfig == 16) { // chi2 - Psi pair variations
    cuts.AddCutPCM("00000123", "00200009397300008150400000", "0163103100900000"); // chi2 50 with psi pair 0.1
    cuts.AddCutPCM("00000123", "00200009397300008850400000", "0163103100900000"); // chi2 20 with psi pair 0.1
    cuts.AddCutPCM("00000123", "00200009397300008280400000", "0163103100900000"); // chi2 30 with psi pair 0.2
    cuts.AddCutPCM("00000123", "00200009397300008260400000", "0163103100900000"); // chi2 30 with psi pair 0.05
  } else if (trainConfig == 17) { // chi2 - Psi pair variations
    cuts.AddCutPCM("00000113", "00200009397300008180400000", "0163103100900000"); // chi2 50 with psi pair 0.2
    cuts.AddCutPCM("00000113", "00200009397300008860400000", "0163103100900000"); // chi2 20 with psi pair 0.05
    cuts.AddCutPCM("00000113", "00200009397300008000400000", "0163103100900000"); // chi2 100, no psi pair
  } else if (trainConfig == 18) { // chi2 - Psi pair variations
    cuts.AddCutPCM("00000123", "00200009397300008180400000", "0163103100900000"); // chi2 50 with psi pair 0.2
    cuts.AddCutPCM("00000123", "00200009397300008860400000", "0163103100900000"); // chi2 20 with psi pair 0.05
    cuts.AddCutPCM("00000123", "00200009397300008000400000", "0163103100900000"); // chi2 100, no psi pair
  } else if (trainConfig == 19) { // photon quality
    cuts.AddCutPCM("00000113", "00200009397300008250420000", "0163103100900000"); // photon quality 1
    cuts.AddCutPCM("00000113", "00200009397300008250430000", "0163103100900000"); // photon quality 2
    cuts.AddCutPCM("00000113", "00200009397300008250440000", "0163103100900000"); // photon quality 3
    cuts.AddCutPCM("00000113", "00200009397300008250400000", "0163106100900000"); // tighter alpha meson 0.85
  } else if (trainConfig == 20) { // photon quality added signals
    cuts.AddCutPCM("00000123", "00200009397300008250420000", "0163103100900000"); // photon quality 1
    cuts.AddCutPCM("00000123", "00200009397300008250430000", "0163103100900000"); // photon quality 2
    cuts.AddCutPCM("00000123", "00200009397300008250440000", "0163103100900000"); // photon quality 3
    cuts.AddCutPCM("00000123", "00200009397300008250400000", "0163106100900000"); // tighter alpha meson 0.85
  } else if (trainConfig == 21) { // much stricter min R cut
    cuts.AddCutPCM("00000113", "00700009397300008250400000", "0163103100900000"); // min R = 35 cm
    cuts.AddCutPCM("00000113", "00700009397300008250420000", "0163103100900000"); // photon quality 1, min R = 35 cm
    cuts.AddCutPCM("00000113", "00700009397300008250430000", "0163103100900000"); // photon quality 2, min R = 35 cm
    cuts.AddCutPCM("00000113", "00700009397300008250440000", "0163103100900000"); // photon quality 3, min R = 35 cm
  } else if (trainConfig == 22) { // much stricter min R cut added signals
    cuts.AddCutPCM("00000123", "00700009397300008250400000", "0163103100900000"); // min R = 35 cm
    cuts.AddCutPCM("00000123", "00700009397300008250420000", "0163103100900000"); // photon quality 1, min R = 35 cm
    cuts.AddCutPCM("00000123", "00700009397300008250430000", "0163103100900000"); // photon quality 2, min R = 35 cm
    cuts.AddCutPCM("00000123", "00700009397300008250440000", "0163103100900000"); // photon quality 3, min R = 35 cm
  } else if (trainConfig == 23) { // meson cut variations
    cuts.AddCutPCM("00000113", "00200009397300008250400000", "0163203100900000"); // y meson < 0.7
    cuts.AddCutPCM("00000113", "00200009397300008250400000", "0163503100900000"); // y meson < 0.85
    cuts.AddCutPCM("00000113", "00200009397300008250400000", "0263203100900000"); // mixed event with track mult
    cuts.AddCutPCM("00000113", "00200009397300008250400000", "0063503100900000"); // BG with rotation
  } else if (trainConfig == 24) { // meson cut variations added signals
    cuts.AddCutPCM("00000123", "00200009397300008250400000", "0163203100900000"); // y meson < 0.7
    cuts.AddCutPCM("00000123", "00200009397300008250400000", "0163503100900000"); // y meson < 0.85
    cuts.AddCutPCM("00000123", "00200009397300008250400000", "0263203100900000"); // mixed event with track mult
    cuts.AddCutPCM("00000123", "00200009397300008250400000", "0063503100900000"); // BG with rotation

  //--------- testing triggers for EMC in 2.76TeV ----------------------------------------------
  } else if (trainConfig == 25) { //LHC11a
    cuts.AddCutPCM("00003113", "00200009397300008250400000", "0163103100000000","1111121057032220000"); //WSDD test with triggers INT1
    cuts.AddCutPCM("00003113", "00202209397300008250400000", "0163103100000000","1111121057032220000"); //WSDD test with triggers INT1 restricted wide EMC range
    cuts.AddCutPCM("00003113", "00204409397300008250400000", "0163103100000000","1111121057032220000"); //WSDD test with triggers INT1 restricted tight EMC range
    cuts.AddCutPCM("00051113", "00200009397300008250400000", "0163103100000000","1111121057032220000"); //WSDD test with triggers EMC1
    cuts.AddCutPCM("00051113", "00202209397300008250400000", "0163103100000000","1111121057032220000"); //WSDD test with triggers EMC1 restricted wide EMC range
    cuts.AddCutPCM("00051113", "00204409397300008250400000", "0163103100000000","1111121057032220000"); //WSDD test with triggers EMC1 restricted tight EMC range
  } else if (trainConfig == 26) { //LHC13g full acceptance
    cuts.AddCutPCM("00010113", "00200009397300008250400000", "0163103100000000","1111121067032220000"); //INT7
    cuts.AddCutPCM("00052113", "00200009397300008250400000", "0163103100000000","1111121067032220000"); //EMC7
    cuts.AddCutPCM("00085113", "00200009397300008250400000", "0163103100000000","1111121067032220000"); //EG2
    cuts.AddCutPCM("00083113", "00200009397300008250400000", "0163103100000000","1111121067032220000"); //EG1
  } else if (trainConfig == 27) { //LHC13g loose EMC acceptance for PCM photons
    cuts.AddCutPCM("00010113", "00202209397300008250400000", "0163103100000000","1111121067032220000"); //INT7
    cuts.AddCutPCM("00052113", "00202209397300008250400000", "0163103100000000","1111121067032220000"); //EMC7
    cuts.AddCutPCM("00085113", "00202209397300008250400000", "0163103100000000","1111121067032220000"); //EG2
    cuts.AddCutPCM("00083113", "00202209397300008250400000", "0163103100000000","1111121067032220000"); //EG1
  } else if (trainConfig == 28) { //LHC13g tight EMC acceptance for PCM photons
    cuts.AddCutPCM("00010113", "00204409397300008250400000", "0163103100000000","1111121067032220000"); //INT7
    cuts.AddCutPCM("00052113", "00204409397300008250400000", "0163103100000000","1111121067032220000"); //EMC7
    cuts.AddCutPCM("00085113", "00204409397300008250400000", "0163103100000000","1111121067032220000"); //EG2
    cuts.AddCutPCM("00083113", "00204409397300008250400000", "0163103100000000","1111121067032220000"); //EG1

  //---------  standard configurations for 2.76TeV V00R with SDD -----------------------------
  } else if(trainConfig == 30){    // various standard cuts
    cuts.AddCutPCM("00003113", "00200009397300008250400000", "0163103100900000"); // new pi0/eta cut 2.76TeV
    cuts.AddCutPCM("00003113", "00200009397300008250400000", "0163103100000000"); // new pi0/eta cut 2.76TeV without MC smearing
    cuts.AddCutPCM("00003113", "00200009366300003800000000", "0163103100900000"); // standard cut Pi0 pp 2.76TeV PbPb paper 2012
  } else if (trainConfig == 31) {
     cuts.AddCutPCM("00003113", "00200009366300003800000000", "0163103100900000"); //standard cut Pi0 pp 2.76TeV with SDD , only Minbias MC
  } else if (trainConfig == 32) {
     cuts.AddCutPCM("00003123", "00200009366300003800000000", "0163103100900000"); //standard cut Pi0 pp 2.76TeV with SDD , only Boxes MC


  //--------- Ana marin: variations for eta reanlysis 2.76TeV 2015 ----------------------------
  } else if (trainConfig == 40) {
    cuts.AddCutPCM("00000113", "00200009217302008250404000", "0152101500000000"); //New standard cut for eta: alpha pT dependent
    cuts.AddCutPCM("00000113", "00200009217302008250404000", "0152106500000000"); // alpha variation 0.8
    cuts.AddCutPCM("00000113", "00200009217302008250404000", "0152103500000000"); // alpha variation  1
    cuts.AddCutPCM("00000113", "00200049217302008250404000", "0152101500000000"); // single pT cut 0.075
  } else if (trainConfig == 41) {
    cuts.AddCutPCM("00000123", "00200009217302008250404000", "0152101500000000"); //New standard cut for eta: alpha pT dependent
    cuts.AddCutPCM("00000123", "00200009217302008250404000", "0152106500000000"); // alpha variation  0.8
    cuts.AddCutPCM("00000123", "00200009217302008250404000", "0152103500000000"); // alpha variation  1
    cuts.AddCutPCM("00000123", "00200049217302008250404000", "0152101500000000"); // single pT cut 0.075
  } else if (trainConfig == 42) {
    cuts.AddCutPCM("00000113", "00200019217302008250404000", "0152101500000000"); // single pT cut 0.1
    cuts.AddCutPCM("00000113", "00200008217302008250404000", "0152101500000000"); // TPC cls 0.35
    cuts.AddCutPCM("00000113", "00200006217302008250404000", "0152101500000000"); // TPC cls 0.7
    cuts.AddCutPCM("00000113", "00200009217302009250404000", "0152101500000000"); // qT cut 0.03 2D
  } else if (trainConfig == 43) {
    cuts.AddCutPCM("00000123", "00200019217302008250404000", "0152101500000000"); // single pT cut 0.1
    cuts.AddCutPCM("00000123", "00200008217302008250404000", "0152101500000000"); // TPC cls 0.35
    cuts.AddCutPCM("00000123", "00200006217302008250404000", "0152101500000000"); // TPC cls 0.7
    cuts.AddCutPCM("00000123", "00200009217302009250404000", "0152101500000000"); // qT cut 0.03 2D
  } else if (trainConfig == 44) {
    cuts.AddCutPCM("00000113", "00200009217302008210404000", "0152101500000000"); // variation chi2 30 psi pair 0.1 1D
    cuts.AddCutPCM("00000113", "00200009217302008180404000", "0152101500000000"); // variation chi2 50 psi pair 0.2 2D
    cuts.AddCutPCM("00000113", "00200009217302008860404000", "0152101500000000"); // variation chi2 20 psi pair 0.05 2D
    cuts.AddCutPCM("00000113", "00200009217302009250404000", "0252101500000000"); // variation BG scheme track mult
  } else if (trainConfig == 45) {
    cuts.AddCutPCM("00000123", "00200009217302008210404000", "0152101500000000"); // variation chi2 30 psi pair 0.1 1D
    cuts.AddCutPCM("00000123", "00200009217302008180404000", "0152101500000000"); // variation chi2 50 psi pair 0.2 2D
    cuts.AddCutPCM("00000123", "00200009217302008860404000", "0152101500000000"); // variation chi2 20 psi pair 0.05 2D
    cuts.AddCutPCM("00000123", "00200009217302009250404000", "0252101500000000"); // variation BG scheme track mult
  } else if (trainConfig == 46) {
    cuts.AddCutPCM("00000113", "00200009227302008250404000", "0152101500000000"); //New standard cut for eta: alpha pT dependent
    cuts.AddCutPCM("00000113", "00200009217302002250404000", "0152101500000000"); // qT
    cuts.AddCutPCM("00000113", "00200009217302009250404000", "0152101500000000"); // qT
    cuts.AddCutPCM("00000113", "00200009217302008250004000", "0152101500000000"); // cosPA
  } else if (trainConfig == 47) {
    cuts.AddCutPCM("00000123", "00200009227302008250404000", "0152101500000000"); //New standard cut for eta: alpha pT dependent
    cuts.AddCutPCM("00000123", "00200009217302002250404000", "0152101500000000"); // qT
    cuts.AddCutPCM("00000123", "00200009217302009250404000", "0152101500000000"); // qT
    cuts.AddCutPCM("00000123", "00200009217302008250004000", "0152101500000000"); // cosPA
  } else if (trainConfig == 48) {
    cuts.AddCutPCM("00000113", "00200009317302008250404000", "0152101500000000"); //dEdx variation
    cuts.AddCutPCM("00000113", "00200009617302008250404000", "0152101500000000"); //dEdx variation
    cuts.AddCutPCM("00000113", "00200009215302008250404000", "0152101500000000"); //dEdx variation
    cuts.AddCutPCM("00000113", "00200009210302008250404000", "0152101500000000"); //dEdx variation
  } else if (trainConfig == 49) {
    cuts.AddCutPCM("00000123", "00200009317302008250404000", "0152101500000000"); //dEdx variation
    cuts.AddCutPCM("00000123", "00200009617302008250404000", "0152101500000000"); //dEdx variation
    cuts.AddCutPCM("00000123", "00200009215302008250404000", "0152101500000000"); //dEdx variation
    cuts.AddCutPCM("00000123", "00200009210302008250404000", "0152101500000000"); //dEdx variation


    //---------configs for V0AND 8TeV --------------------------//
  } else if (trainConfig == 60) {
    cuts.AddCutPCM("00010113", "00200009227300008250404000", "0152103500000000"); //New standard cut for 8TeV analysis V0AND with double counting cut, TOF removed
    cuts.AddCutPCM("00010013", "00200009227300008250404000", "0152103500000000"); // no SPD pileup cut
    cuts.AddCutPCM("00010113", "00100009227300008250404000", "0152103500000000"); // R cut 2.8 -180 cm
    cuts.AddCutPCM("00010113", "00500009227300008250404000", "0152103500000000"); // R cut 10. -180 cm
  } else if (trainConfig == 61) {
    cuts.AddCutPCM("00010113", "00200069227300008250404000", "0152103500000000"); // min pT 40 MeV
    cuts.AddCutPCM("00010113", "00200049227300008250404000", "0152103500000000"); // min pT 75 MeV
    cuts.AddCutPCM("00010113", "00200019227300008250404000", "0152103500000000"); // min pT 100MeV
  } else if (trainConfig == 62) {
    cuts.AddCutPCM("00010113", "00200008227300008250404000", "0152103500000000"); // TPC cluster 35%
    cuts.AddCutPCM("00010113", "00200006227300008250404000", "0152103500000000"); // TPC cluster 70%
    cuts.AddCutPCM("00010113", "00200009227300008250604000", "0152103500000000"); // cosPA 0.9
    cuts.AddCutPCM("00010113", "00200009227300008250304000", "0152103500000000"); // cosPA 0.75
  } else if (trainConfig == 63) {
    cuts.AddCutPCM("00010113", "00200009327300008250404000", "0152103500000000"); // nsig electron   -4,5
    cuts.AddCutPCM("00010113", "00200009627300008250404000", "0152103500000000"); // nsig electron -2.5,4
    cuts.AddCutPCM("00010113", "00200009257300008250404000", "0152103500000000"); // nsig pion 2,-10
    cuts.AddCutPCM("00010113", "00200009217300008250404000", "0152103500000000"); // nsig pion 0,-10
  } else if (trainConfig == 64) {
    cuts.AddCutPCM("00010113", "00200009220300008250404000", "0152103500000000"); // pion nsig min mom 0.50 GeV/c
    cuts.AddCutPCM("00010113", "00200009226300008250404000", "0152103500000000"); // pion nsig min mom 0.25 GeV/c
    cuts.AddCutPCM("00010113", "00200009227600008250404000", "0152103500000000"); // pion nsig max mom 2.00 GeV/c
    cuts.AddCutPCM("00010113", "00200009227100008250404000", "0152103500000000"); // pion nsig max mom 5.00 GeV/c
  } else if (trainConfig == 65) {
    cuts.AddCutPCM("00010113", "00200009227300003250404000", "0152103500000000"); // qT max 0.05 1D
    cuts.AddCutPCM("00010113", "00200009227300002250404000", "0152103500000000"); // qT max 0.06 2D
    cuts.AddCutPCM("00010113", "00200009227300009250404000", "0152103500000000"); // qT max 0.03 2D
  } else if (trainConfig == 66) {
    cuts.AddCutPCM("00010113", "00200009227300008150404000", "0152103500000000"); // chi2 50
    cuts.AddCutPCM("00010113", "00200009227300008850404000", "0152103500000000"); // chi2 20
    cuts.AddCutPCM("00010113", "00200009227300008250400000", "0152103500000000"); // no double counting
    cuts.AddCutPCM("00010113", "00200009227300008250406000", "0152103500000000"); // double count with open angle 0.04
  } else if (trainConfig == 67) {
    cuts.AddCutPCM("00010113", "00200009227300008210404000", "0152103500000000"); // Psi pair 0.1  1D
    cuts.AddCutPCM("00010113", "00200009227300008260404000", "0152103500000000"); // Psi pair 0.05 2D
    cuts.AddCutPCM("00010113", "00200009227300008280404000", "0152103500000000"); // Psi pair 0.2  2D
  } else if (trainConfig == 68) {
    cuts.AddCutPCM("00010113", "00200009227300008250404000", "0252103500000000"); // variation BG scheme track mult
    cuts.AddCutPCM("00010113", "00200009227300008250404000", "0152107500000000"); // alpha meson 0.85
    cuts.AddCutPCM("00010113", "00200009227300008250404000", "0152105500000000"); // alpha meson 0.75
  } else if (trainConfig == 69) {
    cuts.AddCutPCM("00010113", "00200009227300008250404000", "0152103500000000"); //New standard cut for 8TeV analysis V0AND with double counting cut, TOF removed
    cuts.AddCutPCM("00010213", "00200009227300008250404000", "0152103500000000"); //same as above + maximum past future rejection
    cuts.AddCutPCM("00010513", "00200009227300008250404000", "0152103500000000"); //same as above + medium past future rejection

  //---------configs for 8TeV triggers --------------------------//
  } else if (trainConfig == 70) { //pp 8TeV cuts with EMC triggers
    cuts.AddCutPCM("00010113", "00200009227302008250400000", "0152103500000000","1111111067032220000"); // standard cut 8tev
    cuts.AddCutPCM("00052113", "00200009227302008250400000", "0152103500000000","1111111067032220000"); // trigger kEMC7
    cuts.AddCutPCM("00081113", "00200009227302008250400000", "0152103500000000","1111111067032220000"); // trigger kEGA7
  } else if (trainConfig == 71) { //pp 8TeV cuts with EMC triggers, restricted phi regio EMC tight
    cuts.AddCutPCM("00010113", "00204409227302008250400000", "0152103500000000","1111111067032220000"); // standard cut 8tev
    cuts.AddCutPCM("00052113", "00204409227302008250400000", "0152103500000000","1111111067032220000"); // trigger kEMC7
    cuts.AddCutPCM("00081113", "00204409227302008250400000", "0152103500000000","1111111067032220000"); // trigger kEGA7
  } else if (trainConfig == 72) { //pp 8TeV cuts with EMC triggers, restricted phi regio EMC wide
    cuts.AddCutPCM("00010113", "00202209227302008250400000", "0152103500000000","1111111067032220000"); // standard cut 8tev
    cuts.AddCutPCM("00052113", "00202209227302008250400000", "0152103500000000","1111111067032220000"); // trigger kEMC7
    cuts.AddCutPCM("00081113", "00202209227302008250400000", "0152103500000000","1111111067032220000"); // trigger kEGA7
  } else if (trainConfig == 74) { //pp 8TeV std cuts with EMC triggers
    cuts.AddCutPCM("00010113", "00200009227300008250404000", "0152103500000000","1111111067032220000"); // new standard cut 8tev
    cuts.AddCutPCM("00052113", "00200009227300008250404000", "0152103500000000","1111111067032220000"); // trigger kEMC7
    cuts.AddCutPCM("00081113", "00200009227300008250404000", "0152103500000000","1111111067032220000"); // trigger kEGA7
  } else if (trainConfig == 75) { //pp 8TeV cuts with smearing
    cuts.AddCutPCM("00010113", "00200009227300008250404000", "0152103500000000"); //New standard cut for 8TeV analysis V0AND with double counting cut, TOF removed
    cuts.AddCutPCM("00010113", "00200009227300008250404000", "0152103500800000"); //smearing 8
    cuts.AddCutPCM("00010113", "00200009227300008250404000", "0152103500900000"); //smearing 9


    //---------systematic studies mesons and direct photon 2016 pp 7TeV------------------//
  } else if (trainConfig == 80) {
    cuts.AddCutPCM("00000113", "00200009227300008250404000", "0152103500000000"); //New standard cut for 7TeV analysis V0OR with double counting cut, TOF removed
    cuts.AddCutPCM("00000013", "00200009227300008250404000", "0152103500000000"); // no SPD pileup cut
    cuts.AddCutPCM("00000113", "00100009227300008250404000", "0152103500000000"); // R cut 2.8 -180 cm
    cuts.AddCutPCM("00000113", "00500009227300008250404000", "0152103500000000"); // R cut 10. -180 cm
  } else if (trainConfig == 81) {
    cuts.AddCutPCM("00000113", "00200069227300008250404000", "0152103500000000"); // min pT 40 MeV
    cuts.AddCutPCM("00000113", "00200049227300008250404000", "0152103500000000"); // min pT 75 MeV
    cuts.AddCutPCM("00000113", "00200019227300008250404000", "0152103500000000"); // min pT 100MeV
  } else if (trainConfig == 82) {
    cuts.AddCutPCM("00000113", "00200008227300008250404000", "0152103500000000"); // TPC cluster 35%
    cuts.AddCutPCM("00000113", "00200006227300008250404000", "0152103500000000"); // TPC cluster 70%
    cuts.AddCutPCM("00000113", "00200009227300008250604000", "0152103500000000"); // cosPA 0.9
    cuts.AddCutPCM("00000113", "00200009227300008250304000", "0152103500000000"); // cosPA 0.75
  } else if (trainConfig == 83) {
    cuts.AddCutPCM("00000113", "00200009327300008250404000", "0152103500000000"); // nsig electron   -4,5
    cuts.AddCutPCM("00000113", "00200009627300008250404000", "0152103500000000"); // nsig electron -2.5,4
    cuts.AddCutPCM("00000113", "00200009257300008250404000", "0152103500000000"); // nsig pion 2,-10
    cuts.AddCutPCM("00000113", "00200009217300008250404000", "0152103500000000"); // nsig pion 0,-10
  } else if (trainConfig == 84) {
    cuts.AddCutPCM("00000113", "00200009220300008250404000", "0152103500000000"); // pion nsig min mom 0.50 GeV/c
    cuts.AddCutPCM("00000113", "00200009226300008250404000", "0152103500000000"); // pion nsig min mom 0.25 GeV/c
    cuts.AddCutPCM("00000113", "00200009227600008250404000", "0152103500000000"); // pion nsig max mom 2.00 GeV/c
    cuts.AddCutPCM("00000113", "00200009227100008250404000", "0152103500000000"); // pion nsig max mom 5.00 GeV/c
  } else if (trainConfig == 85) {
    cuts.AddCutPCM("00000113", "00200009227300003250404000", "0152103500000000"); // qT max 0.05 1D
    cuts.AddCutPCM("00000113", "00200009227300002250404000", "0152103500000000"); // qT max 0.06 2D
    cuts.AddCutPCM("00000113", "00200009227300009250404000", "0152103500000000"); // qT max 0.03 2D
  } else if (trainConfig == 86) {
    cuts.AddCutPCM("00000113", "00200009227300008150404000", "0152103500000000"); // chi2 50
    cuts.AddCutPCM("00000113", "00200009227300008850404000", "0152103500000000"); // chi2 20
    cuts.AddCutPCM("00000113", "00200009227300008250400000", "0152103500000000"); // no double counting
    cuts.AddCutPCM("00000113", "00200009227300008250406000", "0152103500000000"); // double count with open angle 0.04
  } else if (trainConfig == 87) {
    cuts.AddCutPCM("00000113", "00200009227300008210404000", "0152103500000000"); // Psi pair 0.1  1D
    cuts.AddCutPCM("00000113", "00200009227300008260404000", "0152103500000000"); // Psi pair 0.05 2D
    cuts.AddCutPCM("00000113", "00200009227300008280404000", "0152103500000000"); // Psi pair 0.2  2D
  } else if (trainConfig == 88) {
    cuts.AddCutPCM("00000113", "00200009227300008250404000", "0252103500000000"); // variation BG scheme track mult
    cuts.AddCutPCM("00000113", "00200009227300008250404000", "0152107500000000"); // alpha meson 0.85
    cuts.AddCutPCM("00000113", "00200009227300008250404000", "0152105500000000"); // alpha meson 0.75

   //--------- pp7TeV purity studies (kappa cut)    -------------------------//
  } else if (trainConfig == 89) {
    cuts.AddCutPCM("00000113", "00200009300000008250404000", "0152103500000000"); // -3 < kappa <  5
    cuts.AddCutPCM("00000113", "00200009500000008250404000", "0152103500000000"); // -5 < kappa < 10
    cuts.AddCutPCM("00000113", "00200009600000008250404000", "0152103500000000"); // -3 < kappa < 10
    cuts.AddCutPCM("00000113", "00200009700000008250404000", "0152103500000000"); //  0 < kappa < 10

  // --------- testing validity of old 7 TeV result -------------------------//
  } else if (trainConfig == 90) {
    cuts.AddCutPCM("00000113", "00200009227302008250400000", "0152103500000000"); //New standard cut for 7TeV analysis V0OR
    cuts.AddCutPCM("00000113", "00200008366300000200000000", "0163103100900000"); //Old standard cut for 7TeV analysis V0OR
    cuts.AddCutPCM("00000113", "00200009360300007800004000", "0263103100900000"); //dalitz: New Standard Only MB, standard pp7Tev cut dalitz
  } else if (trainConfig == 91) {
    cuts.AddCutPCM("00000113", "00200009227302008250400000", "0152103500000000"); //New standard cut for 7TeV analysis V0OR
    cuts.AddCutPCM("00000113", "00200009227302008254400000", "0152103500000000"); //asym pt dep
    cuts.AddCutPCM("00000113", "00200009227302008255400000", "0152103500000000"); //asym tight pt dep
  } else if (trainConfig == 92) {
    cuts.AddCutPCM("00000113", "00200009227302008250404000", "0152103500000000"); //New standard cut for 7TeV analysis V0OR with double counting cut
    cuts.AddCutPCM("00000113", "00200009227302008250404000", "0152503500000000"); //y < 0.85
    cuts.AddCutPCM("00000113", "00200009227302008250404000", "0152303500000000"); //y < 0.60
  } else if (trainConfig == 93) {
    cuts.AddCutPCM("00000113", "00200009227300008250404000", "0152103500000000"); //New standard cut for 7TeV analysis V0OR with double counting cut, TOF removed

  //----------------QA omega analysis---------------------------------------//
  } else if (trainConfig == 94) {
    cuts.AddCutPCM("00000113", "00200009227000008250400000", "0152103500000000"); // std conv cut used in omega analysis

    //  Test RBins for 7TeV and new cuts
  } else if (trainConfig == 95) {
    cuts.AddCutPCM("00000113", "00200009227300008250404000", "0152103500000000"); //New standard cut for 7TeV analysis V0OR with double counting cut, TOF removed
    cuts.AddCutPCM("00000113", "00200009f9730000dge0404000", "0152103500000000"); //New standard cut for 7TeV analysis V0OR with double counting cut, TOF removed
    cuts.AddCutPCM("00000113", "00m00009f9730000dge0404000", "0152103500000000"); //New standard cut for 7TeV analysis V0OR with double counting cut, TOF removed
  } else if (trainConfig == 96) {
    cuts.AddCutPCM("00000113", "00a00009f9730000dge0404000", "0152103500000000"); //New standard cut for 7TeV analysis V0OR with double counting cut, TOF removed
    cuts.AddCutPCM("00000113", "00b00009f9730000dge0404000", "0152103500000000"); //New standard cut for 7TeV analysis V0OR with double counting cut, TOF removed
    cuts.AddCutPCM("00000113", "00c00009f9730000dge0404000", "0152103500000000"); //New standard cut for 7TeV analysis V0OR with double counting cut, TOF removed
  } else if (trainConfig == 97) {  // To be used with MBW
    cuts.AddCutPCM("00000113", "00a00009f9730000dge0404000", "0152103500000000"); //New standard cut for 7TeV analysis V0OR with double counting cut, TOF removed
    cuts.AddCutPCM("00000113", "00b00009f9730000dge0404000", "0152103500000000"); //New standard cut for 7TeV analysis V0OR with double counting cut, TOF removed
    cuts.AddCutPCM("00000113", "00c00009f9730000dge0404000", "0152103500000000"); //New standard cut for 7TeV analysis V0OR with double counting cut, TOF removed



  // ------------------------- run 2 High mult triggers --------------------------------------
  } else if (trainConfig == 100) {
    cuts.AddCutPCM("00074113", "00200009227302008250404000", "0152103500000000"); // for V0 High-Mult trigger
    cuts.AddCutPCM("00075113", "00200009227302008250404000", "0152103500000000"); // for SPD High-Mult trigger
    cuts.AddCutPCM("00074013", "00200009227302008250400000", "0152103500000000"); // check # of entries w/ pileup rejection cut for V0HM
    cuts.AddCutPCM("00075013", "00200009227302008250400000", "0152103500000000"); // check # of entries w/o pileup rejection cut for SPHM


  // ------------------------- run 2 configurations ------------------------------------------
  } else if (trainConfig == 111) {  // standard config for run 2
    cuts.AddCutPCM("00010113", "00200009227302008250400000", "0152103500000000"); //New standard cut for eta analysis
    cuts.AddCutPCM("00010113", "00200009227302008250400000", "0152101500000000"); //variation alpha pT dependent
    cuts.AddCutPCM("00010113", "00200009227302008250400000", "0152109500000000"); //variation alpha
    cuts.AddCutPCM("00010113", "00200009227302008250400000", "0152101500000002"); //variation alpha opan max
  } else if (trainConfig == 112) {  // standard config for run 2 with
    cuts.AddCutPCM("00010113", "00200009227302008250400000", "0152103500900000"); //New standard cut for eta analysis
    cuts.AddCutPCM("00010113", "00200009227302008250400000", "0152101500900000"); //variation alpha pT dependent
    cuts.AddCutPCM("00010113", "00200009227302008250400000", "0152109500900000"); //variation alpha
    cuts.AddCutPCM("00010113", "00200009227302008250400000", "0152101500900002"); //variation alpha opan max
  } else if (trainConfig == 113) {
    cuts.AddCutPCM("00010113", "00200009227302008250400000", "0152101500000000");
    cuts.AddCutPCM("00010113", "00200009227302008250404000", "0152101500000000"); // double counting cut
    cuts.AddCutPCM("00010113", "00200009327302008250400000", "0152101500000000"); // dEdx 4 sigma below e
    cuts.AddCutPCM("00010113", "00200009327302008250404000", "0152101500000000");
  } else if (trainConfig == 114) {
    cuts.AddCutPCM("00010113", "00200009227300008254404000", "0152101500000000"); // 13TeV with asymmetry and pT dep alpha cut
    cuts.AddCutPCM("00010113", "00200009227300008254404000", "0152103500000000"); // 13TeV with asymmetry
    cuts.AddCutPCM("00010113", "00200009227300008250404000", "0152101500000000"); // 13TeV pT dep alpha cut
    cuts.AddCutPCM("00010113", "00200009227300008250404000", "0152103500000000"); // 13TeV
  } else if (trainConfig == 115) {
    cuts.AddCutPCM("00010113", "00200009227302008250404000", "0152101500000000");
    cuts.AddCutPCM("00010113", "00200009327302008250404000", "0152101500000000"); // dEdx 4 sigma below e
    cuts.AddCutPCM("00010113", "00200079227302008250404000", "0152101500000000"); // pT cut at 0
    cuts.AddCutPCM("00010113", "00200079327302008250404000", "0152101500000000");
  } else if (trainConfig == 116) {
    cuts.AddCutPCM("00010113", "00200009227302008254404000", "0152101500000000"); // standard cut Gamma pp 13TeV
    cuts.AddCutPCM("30110113", "00200009227302008254404000", "0152101500000000"); // mult.: 0-5%
    cuts.AddCutPCM("31310113", "00200009227302008254404000", "0152101500000000"); // mult.: 5-15%
    cuts.AddCutPCM("33610113", "00200009227302008254404000", "0152101500000000"); // mult.: 15-30%
  } else if (trainConfig == 117) {
    cuts.AddCutPCM("13510113", "00200009227302008254404000", "0152101500000000"); // mult.: 30-50%
    cuts.AddCutPCM("15010113", "00200009227302008254404000", "0152101500000000"); // mult.: 50-100%
    cuts.AddCutPCM("10110113", "00200009227302008254404000", "0152101500000000"); // mult.: 0-10%
    cuts.AddCutPCM("11010113", "00200009227302008254404000", "0152101500000000"); // mult.: 10-100%
  } else if (trainConfig == 118){
    cuts.AddCutPCM("00010113", "00200009a27300008250a04120", "0152103500000000"); // New standard cut for pp 5 TeV analysis VAND
    cuts.AddCutPCM("00010113", "00200079a27300008250a04120", "0152103500000000"); // min pT no cut
    cuts.AddCutPCM("00010113", "00200069a27300008250a04120", "0152103500000000"); // min pT 40 MeV
    cuts.AddCutPCM("00010113", "00200049a27300008250a04120", "0152103500000000"); // min pT 75 MeV
    cuts.AddCutPCM("00010113", "00200008a27300008250a04120", "0152103500000000"); // TPC cluster 35%
    cuts.AddCutPCM("00010113", "00200006a27300008250a04120", "0152103500000000"); // TPC cluster 70%
  } else if (trainConfig == 119){
    cuts.AddCutPCM("00010113", "00200009a27300008250a04120", "0152103500000000"); // edEdx -3,3
    cuts.AddCutPCM("00010113", "00200009b27300008250a04120", "0152103500000000"); // edEdx -3,2,3.2
    cuts.AddCutPCM("00010113", "00200009c27300008250a04120", "0152103500000000"); // edEdx -2.8,2.8
    cuts.AddCutPCM("00010113", "00200009a57300008250a04120", "0152103500000000"); // pidEdx 2,-10
    cuts.AddCutPCM("00010113", "00200009a17300008250a04120", "0152103500000000"); // pidEdx 0,-10
  } else if (trainConfig == 120){
    cuts.AddCutPCM("00010113", "00200009a20300008250a04120", "0152103500000000"); // pion nsig min mom 0.50 GeV/c
    cuts.AddCutPCM("00010113", "00200009a26300008250a04120", "0152103500000000"); // pion nsig min mom 0.25 GeV/c
    cuts.AddCutPCM("00010113", "00200009a27600008250a04120", "0152103500000000"); // pion nsig max mom 2.00 GeV/c
    cuts.AddCutPCM("00010113", "00200009a27100008250a04120", "0152103500000000"); // pion nsig max mom 5.00 GeV/c
  } else if (trainConfig == 121){
    cuts.AddCutPCM("00010113", "00200009a27300003250a04120", "0152103500000000"); // qT max 0.05 1D
    cuts.AddCutPCM("00010113", "00200009a27300002250a04120", "0152103500000000"); // qT max 0.06 2D
    cuts.AddCutPCM("00010113", "00200009a27300008210a04120", "0152103500000000"); // Psi pair 0.1  1D
    cuts.AddCutPCM("00010113", "00200009a27300008260a04120", "0152103500000000"); // Psi pair 0.05 2D
    cuts.AddCutPCM("00010113", "00020009a27300008280a04120", "0152103500000000"); // Psi pair 0.2  2D
  } else if (trainConfig == 122){
    cuts.AddCutPCM("00010113", "00200009a27300008210a04120", "0152103500000000"); // variation chi2 30 psi pair 0.1 1D
    cuts.AddCutPCM("00010113", "00200009a27300008860a04120", "0152103500000000"); // variation chi2 20 psi pair 0.05 2D
    cuts.AddCutPCM("00010113", "00200009a27300008180a04120", "0152103500000000"); // variation chi2 50 psi pair 0.2 2D
    cuts.AddCutPCM("00010113", "00200009a27300008254a04120", "0152103500000000"); // Photon Asymmetry Cut
  } else if (trainConfig == 123){
    cuts.AddCutPCM("00010113", "00200009a27300008250904120", "0152103500000000"); // CosPA 0.99
    cuts.AddCutPCM("00010113", "00200009a27300008250b04120", "0152103500000000"); // CosPA 0.985
    cuts.AddCutPCM("00010113", "00200009a27300008250a00120", "0152103500000000"); // no double counting
    cuts.AddCutPCM("00010113", "00200009a27300008250a04120", "0152105500000000"); // meson alpha < 0.75
    cuts.AddCutPCM("00010113", "00200009a27300008250a04120", "0152107500000000"); // meson alpha < 0.85
  } else if (trainConfig == 124){
    cuts.AddCutPCM("00010113", "00200009a27300008250a04130", "0152103500000000"); // dcaz 4cm
    cuts.AddCutPCM("00010113", "00200009a27300008250a04140", "0152103500000000"); // dcaz 3cm
    cuts.AddCutPCM("00010113", "00200009a27300008250a04220", "0152103500000000"); // dcar 5cm
    cuts.AddCutPCM("00010113", "00200009a27300008250a04320", "0152105500000000"); // dcar 4cm
  } else if (trainConfig == 125){// w/o map
    cuts.AddCutPCM("00010113", "00200009a27300008250404120", "0152103500000000"); // cosPA 0.85
    cuts.AddCutPCM("00010113", "00200009a27300008250a04120", "0152103500000000"); // cosPA 0.995
    cuts.AddCutPCM("00010113", "0d200009a27300008250404120", "0152103500000000"); // cosPA 0.85  |eta| < 0.8
    cuts.AddCutPCM("00010113", "0d200009a27300008250a04120", "0152103500000000"); // cosPA 0.995 |eta| < 0.8
  } else if (trainConfig == 126){// w/ map
    cuts.AddCutPCM("00010113", "00200009a27300008250404120", "0152103500000000"); // cosPA 0.85
    cuts.AddCutPCM("00010113", "00200009a27300008250a04120", "0152103500000000"); // cosPA 0.995
    cuts.AddCutPCM("00010113", "0d200009a27300008250404120", "0152103500000000"); // cosPA 0.85  |eta| < 0.8
    cuts.AddCutPCM("00010113", "0d200009a27300008250a04120", "0152103500000000"); // cosPA 0.995 |eta| < 0.8
  } else if (trainConfig == 127){// w/ map w/ MBW
    cuts.AddCutPCM("00010113", "00200009a27300008250404120", "0152103500000000"); // cosPA 0.85
    cuts.AddCutPCM("00010113", "00200009a27300008250a04120", "0152103500000000"); // cosPA 0.995
    cuts.AddCutPCM("00010113", "0d200009a27300008250404120", "0152103500000000"); // cosPA 0.85  |eta| < 0.8
    cuts.AddCutPCM("00010113", "0d200009a27300008250a04120", "0152103500000000"); // cosPA 0.995 |eta| < 0.8
    // eta cut variation to exclude central cathode
  } else if (trainConfig == 130) {
    cuts.AddCutPCM("00010113", "00200009227300008250404000", "0152103500000000"); //standard cut from 8TeV ana
    cuts.AddCutPCM("00010113", "0a200009227300008250404000", "0152103500000000"); //eta cut 0.2 < |eta| < 0.9
  } else if (trainConfig == 131) {
    cuts.AddCutPCM("00000113", "00200009227300008250404000", "0152103500000000"); //standard cut from 7TeV ana
    cuts.AddCutPCM("00000113", "0a200009227300008250404000", "0152103500000000"); //eta cut 0.2 < |eta| < 0.9
    // minimum pT studies on electrons and photons
  } else if (trainConfig == 132) {
    cuts.AddCutPCM("00010113", "00200009227300008250404000", "0152103500000000"); //standard cut from 8TeV ana
    cuts.AddCutPCM("00010113", "00200099227300008250404000", "0152103500000000"); // gamma pT > 0.1 GeV/c
    cuts.AddCutPCM("00010113", "002000a9227300008250404000", "0152103500000000"); // gamma pT > 0.15 GeV/c
    cuts.AddCutPCM("00010113", "002000b9227300008250404000", "0152103500000000"); // gamma pT > 0.2 GeV/c
  } else if (trainConfig == 133) {
    cuts.AddCutPCM("00010113", "002000c9227300008250404000", "0152103500000000"); // 8 TeV std, but min electron pT > 0.6 for all configs
    cuts.AddCutPCM("00010113", "002000d9227300008250404000", "0152103500000000"); // gamma pT > 0.1 GeV/c
    cuts.AddCutPCM("00010113", "002000e9227300008250404000", "0152103500000000"); // gamma pT > 0.15 GeV/c
    cuts.AddCutPCM("00010113", "002000f9227300008250404000", "0152103500000000"); // gamma pT > 0.2 GeV/c
  } else if (trainConfig == 134) {
    cuts.AddCutPCM("00000113", "00200009227300008250404000", "0152103500000000"); //standard cut from 7TeV ana
    cuts.AddCutPCM("00000113", "00200099227300008250404000", "0152103500000000"); // gamma pT > 0.1 GeV/c
    cuts.AddCutPCM("00000113", "002000a9227300008250404000", "0152103500000000"); // gamma pT > 0.15 GeV/c
    cuts.AddCutPCM("00000113", "002000b9227300008250404000", "0152103500000000"); // gamma pT > 0.2 GeV/c
  } else if (trainConfig == 135) {
    cuts.AddCutPCM("00000113", "002000c9227300008250404000", "0152103500000000"); // 7 TeV std, but min electron pT > 0.6 for all configs
    cuts.AddCutPCM("00000113", "002000d9227300008250404000", "0152103500000000"); // gamma pT > 0.1 GeV/c
    cuts.AddCutPCM("00000113", "002000e9227300008250404000", "0152103500000000"); // gamma pT > 0.15 GeV/c
    cuts.AddCutPCM("00000113", "002000f9227300008250404000", "0152103500000000"); // gamma pT > 0.2 GeV/c
  } else if (trainConfig == 136) { // 8 TeV
    cuts.AddCutPCM("00010113", "002000g9227300008250404000", "0152103500000000"); // 0.075 GeV + min gamma pT cut of 150 MeV
    cuts.AddCutPCM("00010113", "002000h9227300008250404000", "0152103500000000"); // 0.100 GeV + min gamma pT cut of 200 MeV
    cuts.AddCutPCM("00010113", "002000i9227300008250404000", "0152103500000000"); // 0.150 GeV + min gamma pT cut of 300 MeV
    cuts.AddCutPCM("00010113", "002000j9227300008250404000", "0152103500000000"); // 0.150 GeV
  } else if (trainConfig == 137) { // 7 TeV
    cuts.AddCutPCM("00000113", "002000g9227300008250404000", "0152103500000000"); // 0.075 GeV + min gamma pT cut of 150 MeV
    cuts.AddCutPCM("00000113", "002000h9227300008250404000", "0152103500000000"); // 0.100 GeV + min gamma pT cut of 200 MeV
    cuts.AddCutPCM("00000113", "002000i9227300008250404000", "0152103500000000"); // 0.150 GeV + min gamma pT cut of 300 MeV
    cuts.AddCutPCM("00000113", "002000j9227300008250404000", "0152103500000000"); // 0.150 GeV

    // -------- Material weights -several configs with same sets for running with various mat weights ----------
  } else if (trainConfig == 160) { // like last last two in 70 and dalitz standard 7TeV
    cuts.AddCutPCM("00000113", "00200009227302008250400000", "0152103500000000"); //New standard cut for 7TeV analysis V0OR
    cuts.AddCutPCM("00000113", "00200008366300000200000000", "0163103100900000"); //Old standard cut for 7TeV analysis V0OR
    cuts.AddCutPCM("00000113", "00200009360300007800004000", "0263103100900000"); //dalitz: New Standard Only MB, standard pp7Tev cut dalitz
  } else if (trainConfig == 161) { // like last last two in 70 and dalitz standard 7TeV
    cuts.AddCutPCM("00000113", "00200009227302008250400000", "0152103500000000"); //New standard cut for 7TeV analysis V0OR
    cuts.AddCutPCM("00000113", "00200008366300000200000000", "0163103100900000"); //Old standard cut for 7TeV analysis V0OR
    cuts.AddCutPCM("00000113", "00200009360300007800004000", "0263103100900000"); //dalitz: New Standard Only MB, standard pp7Tev cut dalitz
  } else if (trainConfig == 162) { // like last last two in 70 and dalitz standard 7TeV
    cuts.AddCutPCM("00000113", "00200009227302008250400000", "0152103500000000"); //New standard cut for 7TeV analysis V0OR
    cuts.AddCutPCM("00000113", "00200008366300000200000000", "0163103100900000"); //Old standard cut for 7TeV analysis V0OR
    cuts.AddCutPCM("00000113", "00200009360300007800004000", "0263103100900000"); //dalitz: New Standard Only MB, standard pp7Tev cut dalitz
  } else if (trainConfig == 163) { // like last last two in 70 and dalitz standard 7TeV
    cuts.AddCutPCM("00000113", "00200009227302008250400000", "0152103500000000"); //New standard cut for 7TeV analysis V0OR
    cuts.AddCutPCM("00000113", "00200008366300000200000000", "0163103100900000"); //Old standard cut for 7TeV analysis V0OR
    cuts.AddCutPCM("00000113", "00200009360300007800004000", "0263103100900000"); //dalitz: New Standard Only MB, standard pp7Tev cut dalitz
  } else if (trainConfig == 164) { // like last last two in 70 and dalitz standard 7TeV
    cuts.AddCutPCM("00000113", "00200009227302008250400000", "0152103500000000"); //New standard cut for 7TeV analysis V0OR
    cuts.AddCutPCM("00000113", "00200008366300000200000000", "0163103100900000"); //Old standard cut for 7TeV analysis V0OR
    cuts.AddCutPCM("00000113", "00200009360300007800004000", "0263103100900000"); //dalitz: New Standard Only MB, standard pp7Tev cut dalitz
  } else if (trainConfig == 165) { // like last last two in 70 and dalitz standard 7TeV
    cuts.AddCutPCM("00000113", "00200009227302008250400000", "0152103500000000"); //New standard cut for 7TeV analysis V0OR
    cuts.AddCutPCM("00000113", "00200008366300000200000000", "0163103100900000"); //Old standard cut for 7TeV analysis V0OR
    cuts.AddCutPCM("00000113", "00200009360300007800004000", "0263103100900000"); //dalitz: New Standard Only MB, standard pp7Tev cut dalitz
  } else if (trainConfig == 166) { // like last last two in 70 and dalitz standard 7TeV
    cuts.AddCutPCM("00000113", "00200009227302008250400000", "0152103500000000"); //New standard cut for 7TeV analysis V0OR
    cuts.AddCutPCM("00000113", "00200008366300000200000000", "0163103100900000"); //Old standard cut for 7TeV analysis V0OR
    cuts.AddCutPCM("00000113", "00200009360300007800004000", "0263103100900000"); //dalitz: New Standard Only MB, standard pp7Tev cut dalitz
  } else if (trainConfig == 167) { // like last last two in 70 and dalitz standard 7TeV
    cuts.AddCutPCM("00000113", "00200009227302008250400000", "0152103500000000"); //New standard cut for 7TeV analysis V0OR
    cuts.AddCutPCM("00000113", "00200008366300000200000000", "0163103100900000"); //Old standard cut for 7TeV analysis V0OR
    cuts.AddCutPCM("00000113", "00200009360300007800004000", "0263103100900000"); //dalitz: New Standard Only MB, standard pp7Tev cut dalitz
  } else if (trainConfig == 168) { // like last last two in 70 and dalitz standard 7TeV
    cuts.AddCutPCM("00000113", "00200009227302008250400000", "0152103500000000"); //New standard cut for 7TeV analysis V0OR
    cuts.AddCutPCM("00000113", "00200008366300000200000000", "0163103100900000"); //Old standard cut for 7TeV analysis V0OR
    cuts.AddCutPCM("00000113", "00200009360300007800004000", "0263103100900000"); //dalitz: New Standard Only MB, standard pp7Tev cut dalitz
  } else if (trainConfig == 169) { // like last last two in 70 and dalitz standard 7TeV
    cuts.AddCutPCM("00000113", "00200009227302008250400000", "0152103500000000"); //New standard cut for 7TeV analysis V0OR
    cuts.AddCutPCM("00000113", "00200008366300000200000000", "0163103100900000"); //Old standard cut for 7TeV analysis V0OR
    cuts.AddCutPCM("00000113", "00200009360300007800004000", "0263103100900000"); //dalitz: New Standard Only MB, standard pp7Tev cut dalitz


  // -------------------------- mult cut studies --------------------------------------------
  } else if (trainConfig == 200) { // kMB
    cuts.AddCutPCM("00100113", "00200009227302008250400000", "0152103500000000"); // 0 -2
    cuts.AddCutPCM("01200113", "00200009227302008250400000", "0152103500000000"); // 2 -5
    cuts.AddCutPCM("02300113", "00200009227302008250400000", "0152103500000000"); // 5 -10
    cuts.AddCutPCM("03400113", "00200009227302008250400000", "0152103500000000"); // 10-30
    cuts.AddCutPCM("04500113", "00200009227302008250400000", "0152103500000000"); // 30-100


  // **************************************************************************************************
  // ************************************ 13 TeV configuarations **************************************
  // **************************************************************************************************
    // Min Bias
  } else if (trainConfig == 300) {
    cuts.AddCutPCM("00010113", "00200009297000001280000000", "0152103500000000"); // Min Bias more open cuts
    cuts.AddCutPCM("00010113", "00200009227300008250400000", "0152101500000000"); // Min Bias default cuts 2.76 TeV
  } else if (trainConfig == 301) { // low B
    cuts.AddCutPCM("00010113", "00200089297000001280000000", "0152103500000000"); // Min Bias more open cuts
    cuts.AddCutPCM("00010113", "00200089227300008250400000", "0152101500000000"); // Min Bias default cuts 2.76 TeV
  } else if (trainConfig == 302) { //MB JL.
    cuts.AddCutPCM("00010113", "0d200009327000008250404000", "0163103100000010"); // Min Bias Same Cuts as PCMCalo for PHOS
  } else if (trainConfig == 303) { // low B optimized chi2,PsiPair,qt
    cuts.AddCutPCM("00010113", "00200089297000001280000000", "0152103500000000"); // pt dep Qt, chi2-psipair exp
    cuts.AddCutPCM("00010113", "0020008929700000iih0400000", "0152103500000000"); // pt dep Qt, chi2-psipair exp
    cuts.AddCutPCM("00010113", "0020008929700000i280400000", "0152103500000000"); // pt dep Qt, open chi2-psipair
    cuts.AddCutPCM("00010113", "00200089297000001ih0400000", "0152103500000000"); // chi2-psipair, open Qt

    // High mult triggers
  } else if (trainConfig == 310) {
    cuts.AddCutPCM("00074113", "00200009227300008250404000", "0152103500000000"); // for V0 High-Mult trigger
    cuts.AddCutPCM("00074013", "00200009227300008250404000", "0152103500000000"); // check # of entries w/ pileup rejection cut for V0HM
    cuts.AddCutPCM("00076113", "00200009227300008250404000", "0152103500000000"); // for V0 High-Mult trigger
    cuts.AddCutPCM("00076013", "00200009227300008250404000", "0152103500000000"); // check # of entries w/ pileup rejection cut for V0HM
  } else if (trainConfig == 311) { // low B
    cuts.AddCutPCM("00074113", "00200089227300008250404000", "0152103500000000"); // for V0 High-Mult trigger
    cuts.AddCutPCM("00074013", "00200089227300008250404000", "0152103500000000"); // check # of entries w/ pileup rejection cut for V0HM
    cuts.AddCutPCM("00076113", "00200089227300008250404000", "0152103500000000"); // for V0 High-Mult trigger
    cuts.AddCutPCM("00076013", "00200089227300008250404000", "0152103500000000"); // check # of entries w/ pileup rejection cut for V0HM

  } else if (trainConfig == 312) { // dEdx recalib MB
    cuts.AddCutPCM("00010113", "00200009f9730000dge0400000", "0163103100000010"); // MB
  } else if (trainConfig == 313) { // dEdx recalib TRD Triggers
    cuts.AddCutPCM("00049113", "00200009f9730000dge0400000", "0163103100000010"); // HSE
    cuts.AddCutPCM("0004a113", "00200009f9730000dge0400000", "0163103100000010"); // HQU
    cuts.AddCutPCM("0004b113", "00200009f9730000dge0400000", "0163103100000010"); // HJT
  } else if (trainConfig == 314) { // dEdx recalib PHOS Trigger
    cuts.AddCutPCM("00062113", "00200009f9730000dge0400000", "0163103100000010"); // PHI7
  } else if (trainConfig == 315) { // dEdx recalib PHOS Trigger, /w calo cut number
    cuts.AddCutPCM("00062113", "00200009f9730000dge0400000", "0163103100000010","2446600040012300000"); // PHI7
  } else if (trainConfig == 316) { // dEdx recalib V0M
    cuts.AddCutPCM("00074113", "00200009f9730000dge0400000", "0163103100000010"); // V0M
    cuts.AddCutPCM("00076113", "00200009f9730000dge0400000", "0163103100000010"); // V0M /w pileup condition
  } else if (trainConfig == 317) { // dEdx recalib EMC Triggers
    cuts.AddCutPCM("0008d113", "00200009f9730000dge0400000", "0163103100000010"); // EG1
    cuts.AddCutPCM("0008e113", "00200009f9730000dge0400000", "0163103100000010"); // EG2
  } else if (trainConfig == 318) { // dEdx recalib EMC Triggers, /w calo cut number
    cuts.AddCutPCM("0008d113", "00200009f9730000dge0400000", "0163103100000010", "411791106f032230000"); // EG1
    cuts.AddCutPCM("0008e113", "00200009f9730000dge0400000", "0163103100000010", "411791106f032230000"); // EG2

  // EMCal triggered sets
  } else if (trainConfig == 320) { // EMC triggers +-1000 ns
    cuts.AddCutPCM("00010113", "00200009227300008250404000", "0163103100000000","1111100010032220000"); //INT7
    cuts.AddCutPCM("00052113", "00200009227300008250404000", "0163103100000000","1111100010032220000"); //EMC7
    cuts.AddCutPCM("00085113", "00200009227300008250404000", "0163103100000000","1111100010032220000"); //EG2
    cuts.AddCutPCM("00083113", "00200009227300008250404000", "0163103100000000","1111100010032220000"); //EG1
  } else if (trainConfig == 321) { // EMC triggers -50, +30 ns
    cuts.AddCutPCM("00010113", "00200009227300008250404000", "0163103100000000","1111100060032220000"); //INT7
    cuts.AddCutPCM("00052113", "00200009227300008250404000", "0163103100000000","1111100060032220000"); //EMC7
    cuts.AddCutPCM("00085113", "00200009227300008250404000", "0163103100000000","1111100060032220000"); //EG2
    cuts.AddCutPCM("00083113", "00200009227300008250404000", "0163103100000000","1111100060032220000"); //EG1
  } else if (trainConfig == 322) { // EMC triggers -50, +30 ns
    cuts.AddCutPCM("00010113", "00200009227300008250404000", "0163103100000000","1111100060032220000"); //INT7
    cuts.AddCutPCM("00085113", "00200009227300008250404000", "0163103100000000","1111100060032220000"); //EG2
    cuts.AddCutPCM("00083113", "00200009227300008250404000", "0163103100000000","1111100060032220000"); //EG1
  } else if (trainConfig == 325) { // EMC triggers +-1000 ns low B
    cuts.AddCutPCM("00010113", "00200089227300008250404000", "0163103100000000","1111100010032220000"); //INT7
    cuts.AddCutPCM("00052113", "00200089227300008250404000", "0163103100000000","1111100010032220000"); //EMC7
    cuts.AddCutPCM("00085113", "00200089227300008250404000", "0163103100000000","1111100010032220000"); //EG2
    cuts.AddCutPCM("00083113", "00200089227300008250404000", "0163103100000000","1111100010032220000"); //EG1
  } else if (trainConfig == 326) { // EMC triggers -50, +30 ns low B
    cuts.AddCutPCM("00010113", "00200089227300008250404000", "0163103100000000","1111100060032220000"); //INT7
    cuts.AddCutPCM("00052113", "00200089227300008250404000", "0163103100000000","1111100060032220000"); //EMC7
    cuts.AddCutPCM("00085113", "00200089227300008250404000", "0163103100000000","1111100060032220000"); //EG2
    cuts.AddCutPCM("00083113", "00200089227300008250404000", "0163103100000000","1111100060032220000"); //EG1
  } else if (trainConfig == 327) { // EMC triggers -50, +30 ns low B
    cuts.AddCutPCM("00010113", "00200089227300008250404000", "0163103100000000","1111100060032220000"); //INT7
    cuts.AddCutPCM("00085113", "00200089227300008250404000", "0163103100000000","1111100060032220000"); //EG2
    cuts.AddCutPCM("00083113", "00200089227300008250404000", "0163103100000000","1111100060032220000"); //EG1

  // DCal triggered sets
  } else if (trainConfig == 330){ //DCAL triggers +- 1000ns
    cuts.AddCutPCM("00010113", "00200009227300008250404000", "0163103100000000","3885500010032220000"); //INT7
    cuts.AddCutPCM("00055113", "00200009227300008250404000", "0163103100000000","3885500010032220000"); //DMC7
    cuts.AddCutPCM("00089113", "00200009227300008250404000", "0163103100000000","3885500010032220000"); //DG2
    cuts.AddCutPCM("0008b113", "00200009227300008250404000", "0163103100000000","3885500010032220000"); //DG1
  } else if (trainConfig == 331){ //DCAL triggers -50, +30 ns
    cuts.AddCutPCM("00010113", "00200009227300008250404000", "0163103100000000","3885500060032220000"); //INT7
    cuts.AddCutPCM("00055113", "00200009227300008250404000", "0163103100000000","3885500060032220000"); //DMC7
    cuts.AddCutPCM("00089113", "00200009227300008250404000", "0163103100000000","3885500060032220000"); //DG2
    cuts.AddCutPCM("0008b113", "00200009227300008250404000", "0163103100000000","3885500060032220000"); //DG1
  } else if (trainConfig == 332){ //DCAL triggers -50, +30 ns
    cuts.AddCutPCM("00010113", "00200009227300008250404000", "0163103100000000","3885500060032220000"); //INT7
    cuts.AddCutPCM("00089113", "00200009227300008250404000", "0163103100000000","3885500060032220000"); //DG2
    cuts.AddCutPCM("0008b113", "00200009227300008250404000", "0163103100000000","3885500060032220000"); //DG1
  } else if (trainConfig == 335){ //DCAL triggers +- 1000ns low B
    cuts.AddCutPCM("00010113", "00200089227300008250404000", "0163103100000000","3885500010032220000"); //INT7
    cuts.AddCutPCM("00055113", "00200089227300008250404000", "0163103100000000","3885500010032220000"); //DMC7
    cuts.AddCutPCM("00089113", "00200089227300008250404000", "0163103100000000","3885500010032220000"); //DG2
    cuts.AddCutPCM("0008b113", "00200089227300008250404000", "0163103100000000","3885500010032220000"); //DG1
  } else if (trainConfig == 336){ //DCAL triggers -50, +30 ns low B
    cuts.AddCutPCM("00010113", "00200089227300008250404000", "0163103100000000","3885500060032220000"); //INT7
    cuts.AddCutPCM("00055113", "00200089227300008250404000", "0163103100000000","3885500060032220000"); //DMC7
    cuts.AddCutPCM("00089113", "00200089227300008250404000", "0163103100000000","3885500060032220000"); //DG2
    cuts.AddCutPCM("0008b113", "00200089227300008250404000", "0163103100000000","3885500060032220000"); //DG1

  // PHOS trigered sets
  } else if (trainConfig == 340) { // PHOS triggers
    cuts.AddCutPCM("00010113", "00200009227300008250404000", "0163103100000000","2444400000013300000"); //INT7
    cuts.AddCutPCM("00062113", "00200009227300008250404000", "0163103100000000","2444400000013300000"); //PHI7
  } else if (trainConfig == 341) { // PHOS triggers low B
    cuts.AddCutPCM("00010113", "00200089227300008250404000", "0163103100000000","2444400000013300000"); //INT7
    cuts.AddCutPCM("00062113", "00200089227300008250404000", "0163103100000000","2444400000013300000"); //PHI7

  } else if (trainConfig == 350) {
    cuts.AddCutPCM("00010113", "00200009227300008250400000", "0152101500000000"); // Min Bias default cuts 2.76 TeV
    cuts.AddCutPCM("00049113", "00200009227300008250400000", "0152101500000000"); // TRD HSE trigger with INT7
    cuts.AddCutPCM("0004a113", "00200009227300008250400000", "0152101500000000"); // TRD HQU trigger with INT7
    cuts.AddCutPCM("0004b113", "00200009227300008250400000", "0152101500000000"); // TRD HJT trigger with INT7
    cuts.AddCutPCM("0004c113", "00200009227300008250400000", "0152101500000000"); // TRD HNU trigger with INT7
  } else if (trainConfig == 352) {//Testing clean0 and clean1
    cuts.AddCutPCM("00010113", "00200009ea7300008250404000", "0163103100000000","1111100060032220000"); //INT7
    //Low B
    cuts.AddCutPCM("00010113", "00200089ea7300008250404000", "0163103100000000","1111100060032220000"); //INT7
  } else if (trainConfig == 353) {//Testing clean0 and clean1; same as above to be able to run with different QA settings
    cuts.AddCutPCM("00010113", "00200009ea7300008250404000", "0163103100000000","1111100060032220000"); //INT7
    //Low B
    cuts.AddCutPCM("00010113", "00200089ea7300008250404000", "0163103100000000","1111100060032220000"); //INT7

  //---------configs for V0AND 13TeVLowB --------------------------//
  } else if (trainConfig == 360) {
    cuts.AddCutPCM("00010113", "00200089227300008280404000", "0152103500000000"); //New standard cut for 13TeVLowB analysis V0AND with double counting cut, TOF removed
    cuts.AddCutPCM("00010013", "00200089227300008280404000", "0152103500000000"); // no SPD pileup cut
    cuts.AddCutPCM("00010113", "00100089227300008280404000", "0152103500000000"); // R cut 2.8 -180 cm
    cuts.AddCutPCM("00010113", "00500089227300008280404000", "0152103500000000"); // R cut 10. -180 cm
  } else if (trainConfig == 361) {
    cuts.AddCutPCM("00010113", "00200069227300008280404000", "0152103500000000"); // min pT 40 MeV
    cuts.AddCutPCM("00010113", "00200049227300008280404000", "0152103500000000"); // min pT 75 MeV
    cuts.AddCutPCM("00010113", "00200019227300008280404000", "0152103500000000"); // min pT 100MeV
    cuts.AddCutPCM("00010113", "002000892273000082c0404000", "0152103500000000"); // Psi pair 0.15  2D
  } else if (trainConfig == 362) {
    cuts.AddCutPCM("00010113", "00200088227300008280404000", "0152103500000000"); // TPC cluster 35%
    cuts.AddCutPCM("00010113", "00200086227300008280404000", "0152103500000000"); // TPC cluster 70%
    cuts.AddCutPCM("00010113", "00200089227300008280604000", "0152103500000000"); // cosPA 0.9
    cuts.AddCutPCM("00010113", "00200089227300008280304000", "0152103500000000"); // cosPA 0.75
  } else if (trainConfig == 363) {
    cuts.AddCutPCM("00010113", "00200089327300008280404000", "0152103500000000"); // nsig electron   -4,5
    cuts.AddCutPCM("00010113", "00200089627300008280404000", "0152103500000000"); // nsig electron -2.5,4
    cuts.AddCutPCM("00010113", "00200089257300008280404000", "0152103500000000"); // nsig pion 2,-10
    cuts.AddCutPCM("00010113", "00200089217300008280404000", "0152103500000000"); // nsig pion 0,-10
  } else if (trainConfig == 364) {
    cuts.AddCutPCM("00010113", "00200089220300008280404000", "0152103500000000"); // pion nsig min mom 0.50 GeV/c
    cuts.AddCutPCM("00010113", "00200089226300008280404000", "0152103500000000"); // pion nsig min mom 0.25 GeV/c
    cuts.AddCutPCM("00010113", "00200089227600008280404000", "0152103500000000"); // pion nsig max mom 2.00 GeV/c
    cuts.AddCutPCM("00010113", "00200089227100008280404000", "0152103500000000"); // pion nsig max mom 5.00 GeV/c
  } else if (trainConfig == 365) {
    cuts.AddCutPCM("00010113", "00200089227300003280404000", "0152103500000000"); // qT max 0.05 1D
    cuts.AddCutPCM("00010113", "00200089227300002280404000", "0152103500000000"); // qT max 0.06 2D
    cuts.AddCutPCM("00010113", "00200089227300009280404000", "0152103500000000"); // qT max 0.03 2D
  } else if (trainConfig == 366) {
    cuts.AddCutPCM("00010113", "00200089227300008180404000", "0152103500000000"); // chi2 50
    cuts.AddCutPCM("00010113", "00200089227300008880404000", "0152103500000000"); // chi2 20
    cuts.AddCutPCM("00010113", "00200089227300008280404000", "0252103500000000"); // variation BG scheme track mult
    cuts.AddCutPCM("00010113", "00200089227300008280406000", "0152103500000000"); // double count with open angle 0.04
  } else if (trainConfig == 367) {
    cuts.AddCutPCM("00010113", "00200089227300008210404000", "0152103500000000"); // Psi pair 0.1  1D
    cuts.AddCutPCM("00010113", "00200089227300008280404000", "0152103500000000"); // Psi pair 0.2  2D
    cuts.AddCutPCM("00010113", "002000892273000082a0404000", "0152103500000000"); // Psi pair 0.25  2D
    cuts.AddCutPCM("00010113", "002000892273000082b0404000", "0152103500000000"); // Psi pair 0.3  2D
  } else if (trainConfig == 368) {
    cuts.AddCutPCM("00010113", "00200089227300008280404000", "0152103500000000"); //New standard cut for 13TeVLowB analysis V0AND with double counting cut, TOF removed
    cuts.AddCutPCM("00010213", "00200089227300008280404000", "0152103500000000"); //same as above + maximum past future rejection
    cuts.AddCutPCM("00010513", "00200089227300008280404000", "0152103500000000"); //same as above + medium past future rejection
    cuts.AddCutPCM("00010113", "0d200089227300008280404000", "0152103500000000"); //eta < 0.8

  // ---------------------------------- cut selection for pp 5 TeV 2015+2017 ------------------------------------
  } else if (trainConfig == 400){
    cuts.AddCutPCM("00010113", "00200009f9730000dge0400000", "0152103500000000"); // new standard dEdx, pt dep Qt, chi2-psipair exp
  } else if (trainConfig == 401) {
    cuts.AddCutPCM("00010013", "00200009f9730000dge0400000", "0152103500000000"); // no SPD pileup cut
  // systematic variations for 5 TeV 2015+2017
  } else if (trainConfig == 402) {
    cuts.AddCutPCM("00010113", "0d200009f9730000dge0400000", "0152103500000000"); // eta 0.8
    cuts.AddCutPCM("00010113", "0c200009f9730000dge0400000", "0152103500000000"); // eta 0.85
    cuts.AddCutPCM("00010113", "00100009f9730000dge0400000", "0152103500000000"); // R cut 2.8 -180 cm
    cuts.AddCutPCM("00010113", "00500009f9730000dge0400000", "0152103500000000"); // R cut 10. -180 cm
  } else if (trainConfig == 403) {
    cuts.AddCutPCM("00010113", "00200069f9730000dge0400000", "0152103500000000"); // min pT 40 MeV
    cuts.AddCutPCM("00010113", "00200049f9730000dge0400000", "0152103500000000"); // min pT 75 MeV
    cuts.AddCutPCM("00010113", "00200019f9730000dge0400000", "0152103500000000"); // min pT 100MeV
  } else if (trainConfig == 404) {
    cuts.AddCutPCM("00010113", "00200068f9730000dge0400000", "0152103500000000"); // TPC cluster 35%
    cuts.AddCutPCM("00010113", "00200066f9730000dge0400000", "0152103500000000"); // TPC cluster 70%
    cuts.AddCutPCM("00010113", "00200009f9730000dge0600000", "0152103500000000"); // cosPA 0.9
    cuts.AddCutPCM("00010113", "00200009f9730000dge0300000", "0152103500000000"); // cosPA 0.75
  } else if (trainConfig == 405) {
    cuts.AddCutPCM("00010113", "0020000939730000dge0400000", "0152103500000000"); // nsig electron   -4,5
    cuts.AddCutPCM("00010113", "0020000969730000dge0400000", "0152103500000000"); // nsig electron -2.5,4
    cuts.AddCutPCM("00010113", "00200009f5730000dge0400000", "0152103500000000"); // nsig pion 2,-10
    cuts.AddCutPCM("00010113", "00200009f1730000dge0400000", "0152103500000000"); // nsig pion 0,-10
  } else if (trainConfig == 406) {
    cuts.AddCutPCM("00010113", "00200009f9030000dge0400000", "0152103500000000"); // pion nsig min mom 0.50 GeV/c
    cuts.AddCutPCM("00010113", "00200009f9630000dge0400000", "0152103500000000"); // pion nsig min mom 0.25 GeV/c
    cuts.AddCutPCM("00010113", "00200009f9760000dge0400000", "0152103500000000"); // pion nsig max mom 2.00 GeV/c
    cuts.AddCutPCM("00010113", "00200009f9710000dge0400000", "0152103500000000"); // pion nsig max mom 5.00 GeV/c
  } else if (trainConfig == 407) {
    cuts.AddCutPCM("00010113", "00200009f9730000age0400000", "0152103500000000"); // qT max 0.040, qtptmax 0.11
    cuts.AddCutPCM("00010113", "00200009f9730000ege0400000", "0152103500000000"); // qT max 0.060, qtptmax 0.14
    cuts.AddCutPCM("00010113", "00200009f9730000fge0400000", "0152103500000000"); // qT max 0.070, qtptmax 0.16
  } else if (trainConfig == 408) {
    cuts.AddCutPCM("00010113", "00200009f9730000d1e0400000", "0152103500000000"); // chi2 50
    cuts.AddCutPCM("00010113", "00200009f9730000dfe0400000", "0152103500000000"); // chi2 50 chi2 dep -0.065
    cuts.AddCutPCM("00010113", "00200009f9730000dhe0400000", "0152103500000000"); // chi2 50 chi2 dep -0.050
    cuts.AddCutPCM("00010113", "00200009f9730000dge0404000", "0152103500000000"); // reject close v0
    cuts.AddCutPCM("00010113", "00200009f9730000dge0406000", "0152103500000000"); // double count with open angle 0.04
  } else if (trainConfig == 409) {
    cuts.AddCutPCM("00010113", "00200009f9730000dgd0400000", "0152103500000000"); // Psi pair 0.15 dep
    cuts.AddCutPCM("00010113", "00200009f9730000dgf0400000", "0152103500000000"); // Psi pair 0.20 dep
    cuts.AddCutPCM("00010113", "00200009f9730000dgg0400000", "0152103500000000"); // Psi pair 0.30 dep
  } else if (trainConfig == 410) {
    cuts.AddCutPCM("00010113", "00200009f9730000dge0400000", "0252103500000000"); // variation BG scheme track mult
    cuts.AddCutPCM("00010113", "00200009f9730000dge0400000", "0152107500000000"); // alpha meson 0.85
    cuts.AddCutPCM("00010113", "00200009f9730000dge0400000", "0152105500000000"); // alpha meson 0.75
    cuts.AddCutPCM("00010113", "00200009227300008250404000", "0152103500000000"); // old cuts (run1)
  } else if (trainConfig == 411){
    cuts.AddCutPCM("00010113", "0dm00009f9730000dge0404000", "0152103500000000"); // excluding 55-72 cm, eta 0.8
  } else if (trainConfig == 412){
    cuts.AddCutPCM("00010113", "00200009f9730000dge0400000", "0162103500000000");  
    cuts.AddCutPCM("00010113", "00200009f9730000dge0400000", "0r62103500000000");  
    cuts.AddCutPCM("00010113", "00200009f9730000dge0400000", "0s62103500000000");    
    cuts.AddCutPCM("00010113", "00200009f9730000dge0400000", "0t62103500000000");    
    cuts.AddCutPCM("00010113", "00200009f9730000dge0400000", "0u62103500000000");    
  } else if (trainConfig == 413){
    cuts.AddCutPCM("00010113", "00200009f9730000dge0400000", "0162103500800000");  
    cuts.AddCutPCM("00010113", "00200009f9730000dge0400000", "0r62103500800000");  
    cuts.AddCutPCM("00010113", "00200009f9730000dge0400000", "0s62103500800000");    
    cuts.AddCutPCM("00010113", "00200009f9730000dge0400000", "0t62103500800000");    
    cuts.AddCutPCM("00010113", "00200009f9730000dge0400000", "0u62103500800000");    
  } else if (trainConfig == 414){
    cuts.AddCutPCM("00010113", "00200009f9730000dge0400000", "0162103500900000");  
    cuts.AddCutPCM("00010113", "00200009f9730000dge0400000", "0r62103500900000");  
    cuts.AddCutPCM("00010113", "00200009f9730000dge0400000", "0s62103500900000");    
    cuts.AddCutPCM("00010113", "00200009f9730000dge0400000", "0t62103500900000");    
    cuts.AddCutPCM("00010113", "00200009f9730000dge0400000", "0u62103500900000");    

  } else if (trainConfig == 422){   // AM changed from 441 to 422  // updated 190923 with better 2D cuts
    cuts.AddCutPCM("00010113", "0da00009f9730000dge0404000", "0152103500000000"); // Standard cut for pp 5 TeV analysis VAND, R 5-33.5
    cuts.AddCutPCM("00010113", "0db00009f9730000dge0404000", "0152103500000000"); // Standard cut for pp 5 TeV analysis VAND  R 33.5-72
    cuts.AddCutPCM("00010113", "0dc00009f9730000dge0404000", "0152103500000000"); // Standard cut for pp 5 TeV analysis VAND  R 72-180
  } else if (trainConfig == 423){   // AM split the RBin "a" in two
    cuts.AddCutPCM("00010113", "0dh00009f9730000dge0404000", "0152103500000000"); // Standard cut for pp 5 TeV analysis VAND, R 5-13
    cuts.AddCutPCM("00010113", "0di00009f9730000dge0404000", "0152103500000000"); // Standard cut for pp 5 TeV analysis VAND  R 13-33.5
    cuts.AddCutPCM("00010113", "0dl00009f9730000dge0404000", "0152103500000000"); // Standard cut for pp 5 TeV analysis VAND  R 72-95
  } else if (trainConfig == 424){   // AM split the RBin "a" in two
    cuts.AddCutPCM("00010113", "0dj00009f9730000dge0404000", "0152103500000000"); // Standard cut for pp 5 TeV analysis VAND, R 33-55
    cuts.AddCutPCM("00010113", "0dk00009f9730000dge0404000", "0152103500000000"); // Standard cut for pp 5 TeV analysis VAND  R 55-72
    cuts.AddCutPCM("00010113", "0dg00009f9730000dge0404000", "0152103500000000"); // Standard cut for pp 5 TeV analysis VAND  R 95-180
  } else if (trainConfig == 425){   // AM exclude 55-72 region
    cuts.AddCutPCM("00010113", "0dm00009f9730000dge0404000", "0152103500000000"); // Standard cut for pp 5 TeV analysis VAND, R  5-180, 55-72 excluded
    cuts.AddCutPCM("00010113", "0dd00009f9730000dge0404000", "0152103500000000"); // Standard cut for pp 5 TeV analysis VAND, R 5-55
    cuts.AddCutPCM("00010113", "0d200009f9730000dge0404000", "0152103500000000"); // Standard cut for pp 5 TeV analysis VAND, R 5-180

  } else if (trainConfig == 440){ // as 400 to be used MBW   // updated 190923 with better 2D cuts
    cuts.AddCutPCM("00010113", "00200009f9730000dge0404000", "0152103500000000"); // Standard cut for pp 5 TeV analysis VAND
    cuts.AddCutPCM("00010113", "0c200009f9730000dge0404000", "0152103500000000"); // Standard cut for pp 5 TeV analysis VAND
    cuts.AddCutPCM("00010113", "0d200009f9730000dge0404000", "0152103500000000"); // Standard cut for pp 5 TeV analysis VAND
  } else if (trainConfig == 442){ // as 422 (before 441) to be used MBW
    cuts.AddCutPCM("00010113", "0da00009f9730000dge0404000", "0152103500000000"); // Standard cut for pp 5 TeV analysis VAND  R 5-33.5
    cuts.AddCutPCM("00010113", "0db00009f9730000dge0404000", "0152103500000000"); // Standard cut for pp 5 TeV analysis VAND  R 33.5-72.
    cuts.AddCutPCM("00010113", "0dc00009f9730000dge0404000", "0152103500000000"); // Standard cut for pp 5 TeV analysis VAND  R 72-180
  } else if (trainConfig == 443){ // as 423 to be used MBW
    cuts.AddCutPCM("00010113", "0dh00009f9730000dge0404000", "0152103500000000"); // Standard cut for pp 5 TeV analysis VAND  R 5-13
    cuts.AddCutPCM("00010113", "0di00009f9730000dge0404000", "0152103500000000"); // Standard cut for pp 5 TeV analysis VAND  R 13-33.5.
    cuts.AddCutPCM("00010113", "0dl00009f9730000dge0404000", "0152103500000000"); // Standard cut for pp 5 TeV analysis VAND  R 72-95
  } else if (trainConfig == 444){   // as 424 to be used with MBW
    cuts.AddCutPCM("00010113", "0dj00009f9730000dge0404000", "0152103500000000"); // Standard cut for pp 5 TeV analysis VAND, R 33-55
    cuts.AddCutPCM("00010113", "0dk00009f9730000dge0404000", "0152103500000000"); // Standard cut for pp 5 TeV analysis VAND  R 55-72
    cuts.AddCutPCM("00010113", "0dg00009f9730000dge0404000", "0152103500000000"); // Standard cut for pp 5 TeV analysis VAND  R 95-180
  } else if (trainConfig == 445){   // AM exclude 55-72 region
    cuts.AddCutPCM("00010113", "0dm00009f9730000dge0404000", "0152103500000000"); // Standard cut for pp 5 TeV analysis VAND, R 5-180 55-72 excluded
    cuts.AddCutPCM("00010113", "0dd00009f9730000dge0404000", "0152103500000000"); // Standard cut for pp 5 TeV analysis VAND, R 5-55
    cuts.AddCutPCM("00010113", "0d200009f9730000dge0404000", "0152103500000000"); // Standard cut for pp 5 TeV analysis VAND, R 5-180

  } else if (trainConfig == 449){ // PCM standard
    cuts.AddCutPCM("00010113", "00200009227300008250404000", "0152103500000000"); // old standard
    cuts.AddCutPCM("00010113", "00200009f97300008250400000", "0152103500000000"); // new standard dEdx
    cuts.AddCutPCM("00010113", "00200009f9730000dge0400000", "0152103500000000"); // new standard dEdx, pt dep Qt, chi2-psipair exp
    cuts.AddCutPCM("00010113", "00200009f9730000d250400000", "0152103500000000"); // new standard dEdx, pt dep Qt, open chi2-psipair
    cuts.AddCutPCM("00010113", "00200009f97300008ge0400000", "0152103500000000"); // new standard dEdx, chi2-psipair, open Qt

  } else if (trainConfig == 450){ // PCM sphericity
    cuts.AddCutPCM("00010113", "00200009f9730000dge0400000", "0152103500000000");
    cuts.AddCutPCM("h0510113", "00200009f9730000dge0400000", "0152103500000000");
    cuts.AddCutPCM("h5a10113", "00200009f9730000dge0400000", "0152103500000000");
    cuts.AddCutPCM("h0a10113", "00200009f9730000dge0400000", "0152103500000000");
    cuts.AddCutPCM("h0310113", "00200009f9730000dge0400000", "0152103500000000");
    cuts.AddCutPCM("h7a10113", "00200009f9730000dge0400000", "0152103500000000");
  } else if (trainConfig == 451){ // PCM V0M multiplicity
    cuts.AddCutPCM("m0110113", "00200009f9730000dge0400000", "0152103500000000"); // 0-1%
    cuts.AddCutPCM("m1510113", "00200009f9730000dge0400000", "0152103500000000"); // 1-5%
    cuts.AddCutPCM("m5k10113", "00200009f9730000dge0400000", "0152103500000000"); // 5-20%
    cuts.AddCutPCM("n2410113", "00200009f9730000dge0400000", "0152103500000000"); // 20-40%
    cuts.AddCutPCM("n4710113", "00200009f9730000dge0400000", "0152103500000000"); // 40-70%
    cuts.AddCutPCM("n7a10113", "00200009f9730000dge0400000", "0152103500000000"); // 70-100%
  } else if (trainConfig == 452){ // PCM SPD multiplicity
    cuts.AddCutPCM("o0110113", "00200009f9730000dge0400000", "0152103500000000"); // 0-1%
    cuts.AddCutPCM("o0210113", "00200009f9730000dge0400000", "0152103500000000"); // 0-2%
    cuts.AddCutPCM("o0510113", "00200009f9730000dge0400000", "0152103500000000"); // 0-5%
    cuts.AddCutPCM("o5k10113", "00200009f9730000dge0400000", "0152103500000000"); // 5-20%
    cuts.AddCutPCM("p2610113", "00200009f9730000dge0400000", "0152103500000000"); // 20-60%
    cuts.AddCutPCM("p6a10113", "00200009f9730000dge0400000", "0152103500000000"); // 60-100%
  } else if (trainConfig == 455){ // PCM Isolated Pi0 analysis
    cuts.AddCutPCM("00010113", "08200009f9730000dge0400000", "4152103500000000");
  } else if (trainConfig == 456){ // PCM HighPtHadron analysis
    cuts.AddCutPCM("00010113", "00200009f9730000dge0400000", "5152103500000000");

  //---------configs for V0AND 8TeV --------------------------//
  } else if (trainConfig == 460) {
    cuts.AddCutPCM("00010113", "00200009f9730000dge0400000", "0152103500000000"); //New standard cut for 8TeV analysis V0AND with double counting cut, TOF removed
    cuts.AddCutPCM("00010013", "00200009f9730000dge0400000", "0152103500000000"); // no SPD pileup cut
    cuts.AddCutPCM("00010113", "00100009f9730000dge0400000", "0152103500000000"); // R cut 2.8 -180 cm
    cuts.AddCutPCM("00010113", "00500009f9730000dge0400000", "0152103500000000"); // R cut 10. -180 cm
  } else if (trainConfig == 461) {
    cuts.AddCutPCM("00010113", "00200069f9730000dge0400000", "0152103500000000"); // min pT 40 MeV
    cuts.AddCutPCM("00010113", "00200049f9730000dge0400000", "0152103500000000"); // min pT 75 MeV
    cuts.AddCutPCM("00010113", "00200019f9730000dge0400000", "0152103500000000"); // min pT 100MeV
  } else if (trainConfig == 462) {
    cuts.AddCutPCM("00010113", "00200068f9730000dge0400000", "0152103500000000"); // TPC cluster 35%
    cuts.AddCutPCM("00010113", "00200066f9730000dge0400000", "0152103500000000"); // TPC cluster 70%
    cuts.AddCutPCM("00010113", "00200009f9730000dge0600000", "0152103500000000"); // cosPA 0.9
    cuts.AddCutPCM("00010113", "00200009f9730000dge0300000", "0152103500000000"); // cosPA 0.75
  } else if (trainConfig == 463) {
    cuts.AddCutPCM("00010113", "0020000939730000dge0400000", "0152103500000000"); // nsig electron   -4,5
    cuts.AddCutPCM("00010113", "0020000969730000dge0400000", "0152103500000000"); // nsig electron -2.5,4
    cuts.AddCutPCM("00010113", "00200009f5730000dge0400000", "0152103500000000"); // nsig pion 2,-10
    cuts.AddCutPCM("00010113", "00200009f1730000dge0400000", "0152103500000000"); // nsig pion 0,-10
  } else if (trainConfig == 464) {
    cuts.AddCutPCM("00010113", "00200009f9030000dge0400000", "0152103500000000"); // pion nsig min mom 0.50 GeV/c
    cuts.AddCutPCM("00010113", "00200009f9630000dge0400000", "0152103500000000"); // pion nsig min mom 0.25 GeV/c
    cuts.AddCutPCM("00010113", "00200009f9760000dge0400000", "0152103500000000"); // pion nsig max mom 2.00 GeV/c
    cuts.AddCutPCM("00010113", "00200009f9710000dge0400000", "0152103500000000"); // pion nsig max mom 5.00 GeV/c
  } else if (trainConfig == 465) {
    cuts.AddCutPCM("00010113", "00200009f9730000age0400000", "0152103500000000"); // qT max 0.040, qt pt max 0.11
    cuts.AddCutPCM("00010113", "00200009f9730000ege0400000", "0152103500000000"); // qT max 0.060, qt pt max 0.14
    cuts.AddCutPCM("00010113", "00200009f9730000fge0400000", "0152103500000000"); // qT max 0.070, qt pt max 0.16
  } else if (trainConfig == 466) {
    cuts.AddCutPCM("00010113", "00200009f9730000d1e0400000", "0152103500000000"); // chi2 50 no chi2 dep.
    cuts.AddCutPCM("00010113", "00200009f9730000dfe0400000", "0152103500000000"); // chi2 50 chi2 dep -0.065
    cuts.AddCutPCM("00010113", "00200009f9730000dhe0400000", "0152103500000000"); // chi2 50 chi2 dep -0.050
    cuts.AddCutPCM("00010113", "00200009f9730000dge0404000", "0152103500000000"); // reject close v0
    cuts.AddCutPCM("00010113", "00200009f9730000dge0406000", "0152103500000000"); // double count with open angle 0.04
  } else if (trainConfig == 467) {
    cuts.AddCutPCM("00010113", "00200009f9730000dgd0400000", "0152103500000000"); // Psi pair 0.15 dep
    cuts.AddCutPCM("00010113", "00200009f9730000dgf0400000", "0152103500000000"); // Psi pair 0.20 dep
    cuts.AddCutPCM("00010113", "00200009f9730000dgg0400000", "0152103500000000"); // Psi pair 0.30 dep
  } else if (trainConfig == 468) {
    cuts.AddCutPCM("00010113", "00200009f9730000dge0400000", "0252103500000000"); // variation BG scheme track mult
    cuts.AddCutPCM("00010113", "00200009f9730000dge0400000", "0152107500000000"); // alpha meson 0.85
    cuts.AddCutPCM("00010113", "00200009f9730000dge0400000", "0152105500000000"); // alpha meson 0.75
    cuts.AddCutPCM("00010113", "00200009227300008250404000", "0152103500000000"); // old cuts (run1)
  } else if (trainConfig == 469) {
    cuts.AddCutPCM("00010213", "00200009f9730000dge0400000", "0152103500000000"); //same as std + maximum past future rejection
    cuts.AddCutPCM("00010513", "00200009f9730000dge0400000", "0152103500000000"); //same as std + medium past future rejection


  //----------------------------- configuration for Jet analysis ----------------------------------------------------
  } else if ( trainConfig == 500){ // Jet analysis pp 5 TeV 2017
    cuts.AddCutPCM("00010113","00200009327000008250400000","2152103500000000"); //
  } else if ( trainConfig == 501){
    cuts.AddCutPCM("00010113","00200009327000008250400000","3152103500000000"); // Jet QA

   //----------------------Cuts by A. Marin for 13 TeV-----------------

 // Low B Field
  } else if (trainConfig == 600) {
    cuts.AddCutPCM("00010113", "00200089227302001280004000", "0152103500000000"); // Min Bias
    cuts.AddCutPCM("00010113", "00200089327302001280004000", "0152103500000000"); // Open dEdx
  } else if (trainConfig == 601) {
    cuts.AddCutPCM("00010113", "00200089327302001280004000", "0152103500000000"); // Min Bias
    cuts.AddCutPCM("00010113", "00200089127302001280004000", "0152103500000000"); // Open dEdx
  } else if (trainConfig == 602) {
    cuts.AddCutPCM("00010113", "00a00089267300008250404000", "0152103500000000"); //
    cuts.AddCutPCM("00010113", "00b00089267300008250404000", "0152103500000000"); //
    cuts.AddCutPCM("00010113", "00c00089267300008250404000", "0152103500000000"); //
  } else if (trainConfig == 603) {    // asymetry cut removed from configs 602-603-604  and 653-654-655 on 14.06.2018 (cuts low pT)
    cuts.AddCutPCM("00010113", "00200089267300008250404000", "0152103500000000"); // Min Bias with photon asym and dedx at high pT, 0.02 lowB
  } else if (trainConfig == 604) {
    cuts.AddCutPCM("00010113", "002000p9267300008250404000", "0152103500000000"); // Min Bias with photon asym and dedx at high pT, pt 0.03
    cuts.AddCutPCM("00010113", "002000s9267300008250404000", "0152103500000000"); // Min Bias with photon asym and dedx at high pT, pT 0.04
    cuts.AddCutPCM("00010113", "00200009267300008250404000", "0152103500000000"); // Min Bias with photon asym and dedx at high pT, pT 0.04
  } else if (trainConfig == 605) {     // pT scan
    cuts.AddCutPCM("00010113", "002000o9267300008250404000", "0152103500000000"); // Min Bias with photon asym and dedx at high pT, pT 0.024
    cuts.AddCutPCM("00010113", "002000q9267300008250404000", "0152103500000000"); // Min Bias with photon asym and dedx at high pT, pT 0.032
    cuts.AddCutPCM("00010113", "002000r9267300008250404000", "0152103500000000"); // Min Bias with photon asym and dedx at high pT, pT 0.036

  } else if (trainConfig == 606) { // asymetry cut removed from configs 702-703-704  and 753-754-755 on 14.06.2018 (cuts low pT)
    cuts.AddCutPCM("00010113", "0d200089267300008250404000", "0152103500000000"); // eta < 0.8
  } else if (trainConfig == 607) { // asymetry cut removed from configs 702-703-704  and 753-754-755 on 14.06.2018 (cuts low pT)
    cuts.AddCutPCM("00010113", "0da00089267300008250404000", "0152103500000000"); // eta < 0.8
    cuts.AddCutPCM("00010113", "0db00089267300008250404000", "0152103500000000"); // eta < 0.8
    cuts.AddCutPCM("00010113", "0dc00089267300008250404000", "0152103500000000"); // eta < 0.8
    // Configs for systematics
  } else if (trainConfig == 608) {    // Systematics of pT
    cuts.AddCutPCM("00010113", "0d200089267300008250404000", "0152103500000000"); // eta < 0.8 standard for 13 TeV
    cuts.AddCutPCM("00010113", "0d200079267300008250404000", "0152103500000000"); // min pT no cut
    cuts.AddCutPCM("00010113", "0d2000p9267300008250404000", "0152103500000000"); // min pT 30 MeV
    cuts.AddCutPCM("00010113", "0d2000s9267300008250404000", "0152103500000000"); // min pT 40 MeV
    cuts.AddCutPCM("00010113", "0d200009267300008250404000", "0152103500000000"); // min pT 50 MeV
  } else if (trainConfig == 609) {
    cuts.AddCutPCM("00010113", "0d200088267300008250404000", "0152103500000000"); // TPC cluster 35%
    cuts.AddCutPCM("00010113", "0d200086267300008250404000", "0152103500000000"); // TPC cluster 70%
  } else if (trainConfig == 610) {
    cuts.AddCutPCM("00010113", "0d200089367300008250404000", "0152103500000000"); // edEdx -4,5
    cuts.AddCutPCM("00010113", "0d200089667300008250404000", "0152103500000000"); // edEdx -2.5,4
    cuts.AddCutPCM("00010113", "0d200089257300008250404000", "0152103500000000"); // pidEdx 2,-10
    cuts.AddCutPCM("00010113", "0d200089217300008250404000", "0152103500000000"); // pidEdx 0,-10
  } else if (trainConfig == 611) {
    cuts.AddCutPCM("00010113", "0d200089260300008250404000", "0152103500000000"); // pion nsig min mom 0.50 GeV/c
    cuts.AddCutPCM("00010113", "0d200089266300008250404000", "0152103500000000"); // pion nsig min mom 0.25 GeV/c
    cuts.AddCutPCM("00010113", "0d200089267600008250404000", "0152103500000000"); // pion nsig max mom 2.00 GeV/c
    cuts.AddCutPCM("00010113", "0d200089267100008250404000", "0152103500000000"); // pion nsig max mom 5.00 GeV/c
  } else if (trainConfig == 612){
    cuts.AddCutPCM("00010113", "0d200089267300003250404000", "0152103500000000"); // qT max 0.05 1D
    cuts.AddCutPCM("00010113", "0d200089267300002250404000", "0152103500000000"); // qT max 0.06 2D
    cuts.AddCutPCM("00010113", "0d200089267300009250404000", "0152103500000000"); // qT max 0.03 2D
    cuts.AddCutPCM("00010113", "0d200089267300008210404000", "0152103500000000"); // Psi pair 0.1  1D
    cuts.AddCutPCM("00010113", "0d200089267300008260404000", "0152103500000000"); // Psi pair 0.05 2D
    cuts.AddCutPCM("00010113", "0d200089267300008280404000", "0152103500000000"); // Psi pair 0.2  2D
  } else if (trainConfig == 614) { // ITS area with large weights split in two
    cuts.AddCutPCM("00010113", "0dh00089267300008250404000", "0152103500000000"); // eta < 0.8   5-13
    cuts.AddCutPCM("00010113", "0di00089267300008250404000", "0152103500000000"); // eta < 0.8  13-33.5
    cuts.AddCutPCM("00010113", "0dl00089267300008250404000", "0152103500000000"); // eta < 0.8  72-95
  } else if (trainConfig == 615) { // scan of meson alpha cut
    cuts.AddCutPCM("00010113", "0d200089267300008250404000", "0152100500000000"); // alpha <0.7
    cuts.AddCutPCM("00010113", "0d200089267300008250404000", "0152108500000000"); // alpha <0.6
    cuts.AddCutPCM("00010113", "0d200089267300008250404000", "015210g500000000"); // alpha <0.5
    cuts.AddCutPCM("00010113", "0d200089267300008250404000", "015210f500000000"); // alpha <0.4
  } else if (trainConfig == 616) {
    cuts.AddCutPCM("00010113", "0d200089267300008250404000", "015210e500000000"); // alpha <0.3
    cuts.AddCutPCM("00010113", "0d200089267300008250404000", "015210a500000000"); // alpha <0.2
    cuts.AddCutPCM("00010113", "0d200089267300008250404000", "015210d500000000"); // alpha <0.1
  } else if (trainConfig == 617) { // RBins studies,
    cuts.AddCutPCM("00010113", "0dj00089267300008250404000", "0152103500000000"); // eta < 0.8  33.5-55 cm
    cuts.AddCutPCM("00010113", "0dk00089267300008250404000", "0152103500000000"); // eta < 0.8  55-72 cm
    cuts.AddCutPCM("00010113", "0dg00089267300008250404000", "0152103500000000"); // eta < 0.8  72-95 cm
  } else if (trainConfig == 618) { // R 5-180 and remove r bin 55-72
    cuts.AddCutPCM("00010113", "0dm00089267300008250404000", "0152103500000000"); // eta < 0.8
    cuts.AddCutPCM("00010113", "0d200089327000008250404000", "0152103500000000"); // eta < 0.8

  } else if (trainConfig == 619) { // R 5-180 and remove r bin 55-72
    cuts.AddCutPCM("00010113", "0d200089f9730000iih0404000", "0152101500000000"); // eta < 0.8  // Test alpha meson pT dependent
    cuts.AddCutPCM("00010113", "0dm00089f9730000iih0404000", "0152103500000000"); // eta < 0.8  // remove  55-72 bin
    cuts.AddCutPCM("00010113", "0dd00089f9730000iih0404000", "0152103500000000"); // eta < 0.8  // use 5-55 bin only
  } else if (trainConfig == 620) { // R 5-180
    cuts.AddCutPCM("00010113", "0d200089f9730000iih0404000", "0152103500000000"); // eta < 0.8  // Test improved cuts
  } else if (trainConfig == 621) { // R 5-180
    cuts.AddCutPCM("00010113", "0da00089f9730000iih0404000", "0152103500000000"); // eta < 0.8  // Test improved cuts
    cuts.AddCutPCM("00010113", "0db00089f9730000iih0404000", "0152103500000000"); // eta < 0.8  // Test improved cuts
    cuts.AddCutPCM("00010113", "0dc00089f9730000iih0404000", "0152103500000000"); // eta < 0.8  // Test improved cuts
  } else if (trainConfig == 622) { // R 5-180
    cuts.AddCutPCM("00010113", "0dh00089f9730000iih0404000", "0152103500000000"); // eta < 0.8  // Test improved cuts
    cuts.AddCutPCM("00010113", "0di00089f9730000iih0404000", "0152103500000000"); // eta < 0.8  // Test improved cuts
    cuts.AddCutPCM("00010113", "0dj00089f9730000iih0404000", "0152103500000000"); // eta < 0.8  // Test improved cuts
  } else if (trainConfig == 623) { // R 5-180
    cuts.AddCutPCM("00010113", "0dk00089f9730000iih0404000", "0152103500000000"); // eta < 0.8  // Test improved cuts
    cuts.AddCutPCM("00010113", "0dl00089f9730000iih0404000", "0152103500000000"); // eta < 0.8  // Test improved cuts
    cuts.AddCutPCM("00010113", "0dg00089f9730000iih0404000", "0152103500000000"); // eta < 0.8  // Test improved cuts
  } else if (trainConfig == 628) { // R 5-180  // Cat 1, cat 2+3
    cuts.AddCutPCM("00010113", "0d200089f9730000iih0424000", "0152103500000000"); // eta < 0.8  // Test improved cuts
    cuts.AddCutPCM("00010113", "0d200089f9730000iih0454000", "0152103500000000"); // eta < 0.8  // Test improved cuts
    cuts.AddCutPCM("00010113", "0dm00089f9730000iih0424000", "0152103500000000"); // eta < 0.8  // Test improved cuts
    cuts.AddCutPCM("00010113", "0dm00089f9730000iih0454000", "0152103500000000"); // eta < 0.8  // Test improved cuts
    cuts.AddCutPCM("00010113", "0dm00089f9730000iih0404000", "0152103520000000"); // eta < 0.8  // Test improved cuts



 // Low B Field to be used with MBW
  } else if (trainConfig == 652) {
    cuts.AddCutPCM("00010113", "00a00089267300008250404000", "0152103500000000"); //
    cuts.AddCutPCM("00010113", "00b00089267300008250404000", "0152103500000000"); //
    cuts.AddCutPCM("00010113", "00c00089267300008250404000", "0152103500000000"); //
  } else if (trainConfig == 653) {
    cuts.AddCutPCM("00010113", "00200089267300008250404000", "0152103500000000"); // Min Bias with photon asym and dedx at high pT, 0.02 lowB
  } else if (trainConfig == 654) {
    cuts.AddCutPCM("00010113", "002000p9267300008250404000", "0152103500000000"); // Min Bias with photon asym and dedx at high pT, pt 0.03
    cuts.AddCutPCM("00010113", "002000s9267300008250404000", "0152103500000000"); // Min Bias with photon asym and dedx at high pT, pT 0.04
    cuts.AddCutPCM("00010113", "00200009267300008250404000", "0152103500000000"); // Min Bias with photon asym and dedx at high pT, pT 0.04
  } else if (trainConfig == 655) {
    cuts.AddCutPCM("00010113", "002000o9267300008250404000", "0152103500000000"); // Min Bias with photon asym and dedx at high pT, pT 0.024
    cuts.AddCutPCM("00010113", "002000q9267300008250404000", "0152103500000000"); // Min Bias with photon asym and dedx at high pT, pT 0.032
    cuts.AddCutPCM("00010113", "002000r9267300008250404000", "0152103500000000"); // Min Bias with photon asym and dedx at high pT, pT 0.036
  } else if (trainConfig == 656) { // asymetry cut removed from configs 702-703-704  and 753-754-755 on 14.06.2018 (cuts low pT)
    cuts.AddCutPCM("00010113", "0d200089267300008250404000", "0152103500000000"); // eta < 0.8
  } else if (trainConfig == 657) { // asymetry cut removed from configs 702-703-704  and 753-754-755 on 14.06.2018 (cuts low pT)
    cuts.AddCutPCM("00010113", "0da00089267300008250404000", "0152103500000000"); // eta < 0.8
    cuts.AddCutPCM("00010113", "0db00089267300008250404000", "0152103500000000"); // eta < 0.8
    cuts.AddCutPCM("00010113", "0dc00089267300008250404000", "0152103500000000"); // eta < 0.8
    // Configs for systematics
  } else if (trainConfig == 658) {    // Systematics of pT
    cuts.AddCutPCM("00010113", "0d200089267300008250404000", "0152103500000000"); // eta < 0.8 standard for 13 TeV
    cuts.AddCutPCM("00010113", "0d200079267300008250404000", "0152103500000000"); // min pT no cut
    cuts.AddCutPCM("00010113", "0d2000p9267300008250404000", "0152103500000000"); // min pT 30 MeV
    cuts.AddCutPCM("00010113", "0d2000s9267300008250404000", "0152103500000000"); // min pT 40 MeV
    cuts.AddCutPCM("00010113", "0d200009267300008250404000", "0152103500000000"); // min pT 50 MeV
  } else if (trainConfig == 659) {
    cuts.AddCutPCM("00010113", "0d200088267300008250404000", "0152103500000000"); // TPC cluster 35%
    cuts.AddCutPCM("00010113", "0d200086267300008250404000", "0152103500000000"); // TPC cluster 70%
  } else if (trainConfig == 660) {
    cuts.AddCutPCM("00010113", "0d200089367300008250404000", "0152103500000000"); // edEdx -4,5
    cuts.AddCutPCM("00010113", "0d200089667300008250404000", "0152103500000000"); // edEdx -2.5,4
    cuts.AddCutPCM("00010113", "0d200089257300008250404000", "0152103500000000"); // pidEdx 2,-10
    cuts.AddCutPCM("00010113", "0d200089217300008250404000", "0152103500000000"); // pidEdx 0,-10
  } else if (trainConfig == 661) {
    cuts.AddCutPCM("00010113", "0d200089260300008250404000", "0152103500000000"); // pion nsig min mom 0.50 GeV/c
    cuts.AddCutPCM("00010113", "0d200089266300008250404000", "0152103500000000"); // pion nsig min mom 0.25 GeV/c
    cuts.AddCutPCM("00010113", "0d200089267600008250404000", "0152103500000000"); // pion nsig max mom 2.00 GeV/c
    cuts.AddCutPCM("00010113", "0d200089267100008250404000", "0152103500000000"); // pion nsig max mom 5.00 GeV/c
  } else if (trainConfig == 662){
    cuts.AddCutPCM("00010113", "0d200089267300003250404000", "0152103500000000"); // qT max 0.05 1D
    cuts.AddCutPCM("00010113", "0d200089267300002250404000", "0152103500000000"); // qT max 0.06 2D
    cuts.AddCutPCM("00010113", "0d200089267300009250404000", "0152103500000000"); // qT max 0.03 2D
    cuts.AddCutPCM("00010113", "0d200089267300008210404000", "0152103500000000"); // Psi pair 0.1  1D
    cuts.AddCutPCM("00010113", "0d200089267300008260404000", "0152103500000000"); // Psi pair 0.05 2D
    cuts.AddCutPCM("00010113", "0d200089267300008280404000", "0152103500000000"); // Psi pair 0.2  2D
  } else if (trainConfig == 664) { // asymetry cut removed from configs 702-703-704  and 753-754-755 on 14.06.2018 (cuts low pT)
    cuts.AddCutPCM("00010113", "0dh00089267300008250404000", "0152103500000000"); // eta < 0.8
    cuts.AddCutPCM("00010113", "0di00089267300008250404000", "0152103500000000"); // eta < 0.8
    cuts.AddCutPCM("00010113", "0dl00089267300008250404000", "0152103500000000"); // eta < 0.8
  } else if (trainConfig == 665) { // scan of meson alpha cut
    cuts.AddCutPCM("00010113", "0d200089267300008250404000", "0152100500000000"); // alpha <0.7
    cuts.AddCutPCM("00010113", "0d200089267300008250404000", "0152108500000000"); // alpha <0.6
    cuts.AddCutPCM("00010113", "0d200089267300008250404000", "015210g500000000"); // alpha <0.5
    cuts.AddCutPCM("00010113", "0d200089267300008250404000", "015210f500000000"); // alpha <0.4
  } else if (trainConfig == 666) {
    cuts.AddCutPCM("00010113", "0d200089267300008250404000", "015210e500000000"); // alpha <0.3
    cuts.AddCutPCM("00010113", "0d200089267300008250404000", "015210a500000000"); // alpha <0.2
    cuts.AddCutPCM("00010113", "0d200089267300008250404000", "015210d500000000"); // alpha <0.1
  } else if (trainConfig == 667) { // RBins studies,
    cuts.AddCutPCM("00010113", "0dj00089267300008250404000", "0152103500000000"); // eta < 0.8  33.5-55 cm
    cuts.AddCutPCM("00010113", "0dk00089267300008250404000", "0152103500000000"); // eta < 0.8  55-72 cm
    cuts.AddCutPCM("00010113", "0dg00089267300008250404000", "0152103500000000"); // eta < 0.8  72-95 cm
  } else if (trainConfig == 668) { // R 5-180 and remove r bin 55-72
    cuts.AddCutPCM("00010113", "0dm00089267300008250404000", "0152103500000000"); // eta < 0.8
    cuts.AddCutPCM("00010113", "0d200089327000008250404000", "0152103500000000"); // eta < 0.8
  } else if (trainConfig == 669) { // R 5-180 and remove r bin 55-72
    cuts.AddCutPCM("00010113", "0d200089f9730000iih0404000", "0152101500000000"); // eta < 0.8  // Test alpha meson pT dependent
    cuts.AddCutPCM("00010113", "0dm00089f9730000iih0404000", "0152103500000000"); // eta < 0.8  // remove  55-72 bin
    cuts.AddCutPCM("00010113", "0dd00089f9730000iih0404000", "0152103500000000"); // eta < 0.8  // use 5-55 bin only
  } else if (trainConfig == 670) { // R 5-180
    cuts.AddCutPCM("00010113", "0d200089f9730000iih0404000", "0152103500000000"); // eta < 0.8  // Test improved cuts
  } else if (trainConfig == 671) { // R 5-180
    cuts.AddCutPCM("00010113", "0da00089f9730000iih0404000", "0152103500000000"); // eta < 0.8  // Test improved cuts
    cuts.AddCutPCM("00010113", "0db00089f9730000iih0404000", "0152103500000000"); // eta < 0.8  // Test improved cuts
    cuts.AddCutPCM("00010113", "0dc00089f9730000iih0404000", "0152103500000000"); // eta < 0.8  // Test improved cuts
  } else if (trainConfig == 672) { // R 5-180
    cuts.AddCutPCM("00010113", "0dh00089f9730000iih0404000", "0152103500000000"); // eta < 0.8  // Test improved cuts
    cuts.AddCutPCM("00010113", "0di00089f9730000iih0404000", "0152103500000000"); // eta < 0.8  // Test improved cuts
    cuts.AddCutPCM("00010113", "0dj00089f9730000iih0404000", "0152103500000000"); // eta < 0.8  // Test improved cuts
  } else if (trainConfig == 673) { // R 5-180
    cuts.AddCutPCM("00010113", "0dk00089f9730000iih0404000", "0152103500000000"); // eta < 0.8  // Test improved cuts
    cuts.AddCutPCM("00010113", "0dl00089f9730000iih0404000", "0152103500000000"); // eta < 0.8  // Test improved cuts
    cuts.AddCutPCM("00010113", "0dg00089f9730000iih0404000", "0152103500000000"); // eta < 0.8  // Test improved cuts
  } else if (trainConfig == 678) { // R 5-180  // Cat 1, cat 2+3
    cuts.AddCutPCM("00010113", "0d200089f9730000iih0424000", "0152103500000000"); // eta < 0.8  // Test improved cuts
    cuts.AddCutPCM("00010113", "0d200089f9730000iih0454000", "0152103500000000"); // eta < 0.8  // Test improved cuts
    cuts.AddCutPCM("00010113", "0dm00089f9730000iih0424000", "0152103500000000"); // eta < 0.8  // Test improved cuts
    cuts.AddCutPCM("00010113", "0dm00089f9730000iih0454000", "0152103500000000"); // eta < 0.8  // Test improved cuts
    cuts.AddCutPCM("00010113", "0dm00089f9730000iih0404000", "0152103520000000"); // eta < 0.8  // Test improved cuts





  // Material studies Ana-----nomB
  } else if (trainConfig == 700) {
    cuts.AddCutPCM("00010113", "00200009266300008854404000", "0152103500000000"); // Min Bias
    cuts.AddCutPCM("00010113", "00200009266300008854404000", "0152101500000000"); // alpha pT dependent and gamma asym cut
    cuts.AddCutPCM("00010113", "00200009266300008284404000", "0152101500000000"); // alpha pT dependent and gamma asym cut , chi2 30
  } else if (trainConfig == 701) { // to be used with MBW
    cuts.AddCutPCM("00010113", "00200009266300008854404000", "0152103500000000"); // Min Bias
    cuts.AddCutPCM("00010113", "00200009266300008854404000", "0152101500000000"); // alpha pT dependent and gamma asym cut
    cuts.AddCutPCM("00010113", "00200009266300008284404000", "0152101500000000"); // alpha pT dependent and gamma asym cut , chi2 30
  } else if (trainConfig == 702) {
    //    cuts.AddCutPCM("00010113", "00200009267300008254404000", "0152103500000000"); // Min Bias with photon asym and dedx at high pT
    //    cuts.AddCutPCM("00010113", "00200009267300008250404000", "0152103500000000"); // Min Bias with  dedx at high pT
    cuts.AddCutPCM("00010113", "00a00009267300008250404000", "0152103500000000"); //
    cuts.AddCutPCM("00010113", "00b00009267300008250404000", "0152103500000000"); //
    cuts.AddCutPCM("00010113", "00c00009267300008250404000", "0152103500000000"); //
  } else if (trainConfig == 703) { // asymetry cut removed from configs 702-703-704  and 753-754-755 on 14.06.2018 (cuts low pT)
    cuts.AddCutPCM("00010113", "00200009267300008250404000", "0152103500000000"); // Min Bias with photon asym and dedx at high pT
  } else if (trainConfig == 704) {
    cuts.AddCutPCM("00010113", "00200049267300008250404000", "0152103500000000"); // Min Bias with photon asym and dedx at high pT, pt 0.075
    cuts.AddCutPCM("00010113", "00200019267300008250404000", "0152103500000000"); // Min Bias with photon asym and dedx at high pT, pT 0.1
    cuts.AddCutPCM("00010113", "00200059267300008250404000", "0152103500000000"); // Min Bias with photon asym and dedx at high pT, pT 0.125
  } else if (trainConfig == 705) {
    cuts.AddCutPCM("00010113", "002000d9267300008250404000", "0152103500000000"); // Min Bias with photon asym and dedx at high pT, pT 0.06
    cuts.AddCutPCM("00010113", "002000m9267300008250404000", "0152103500000000"); // Min Bias with photon asym and dedx at high pT, pT 0.08
    cuts.AddCutPCM("00010113", "002000n9267300008250404000", "0152103500000000"); // Min Bias with photon asym and dedx at high pT, pT 0.09

  } else if (trainConfig == 706) { // asymetry cut removed from configs 702-703-704  and 753-754-755 on 14.06.2018 (cuts low pT)
    cuts.AddCutPCM("00010113", "0d200009267300008250404000", "0152103500000000"); // eta < 0.8
  } else if (trainConfig == 707) { // asymetry cut removed from configs 702-703-704  and 753-754-755 on 14.06.2018 (cuts low pT)
    cuts.AddCutPCM("00010113", "0da00009267300008250404000", "0152103500000000"); // eta < 0.8
    cuts.AddCutPCM("00010113", "0db00009267300008250404000", "0152103500000000"); // eta < 0.8
    cuts.AddCutPCM("00010113", "0dc00009267300008250404000", "0152103500000000"); // eta < 0.8
    // Configs for systematics
  } else if (trainConfig == 708) {    // Systematics of pT
    cuts.AddCutPCM("00010113", "0d200009267300008250404000", "0152103500000000"); // eta < 0.8 standard for 13 TeV
    cuts.AddCutPCM("00010113", "0d200079267300008250404000", "0152103500000000"); // min pT no cut
    cuts.AddCutPCM("00010113", "0d200069267300008250404000", "0152103500000000"); // min pT 40 MeV
    cuts.AddCutPCM("00010113", "0d200049267300008250404000", "0152103500000000"); // min pT 75 MeV
    cuts.AddCutPCM("00010113", "0d200019267300008250404000", "0152103500000000"); // min pT 100 MeV
  } else if (trainConfig == 709) {
    cuts.AddCutPCM("00010113", "0d200008267300008250404000", "0152103500000000"); // TPC cluster 35%
    cuts.AddCutPCM("00010113", "0d200006267300008250404000", "0152103500000000"); // TPC cluster 70%
  } else if (trainConfig == 710) {
    cuts.AddCutPCM("00010113", "0d200009367300008250404000", "0152103500000000"); // edEdx -4,5
    cuts.AddCutPCM("00010113", "0d200009667300008250404000", "0152103500000000"); // edEdx -2.5,4
    cuts.AddCutPCM("00010113", "0d200009257300008250404000", "0152103500000000"); // pidEdx 2,-10
    cuts.AddCutPCM("00010113", "0d200009217300008250404000", "0152103500000000"); // pidEdx 0,-10
  } else if (trainConfig == 711) {
    cuts.AddCutPCM("00010113", "0d200009260300008250404000", "0152103500000000"); // pion nsig min mom 0.50 GeV/c
    cuts.AddCutPCM("00010113", "0d200009266300008250404000", "0152103500000000"); // pion nsig min mom 0.25 GeV/c
    cuts.AddCutPCM("00010113", "0d200009267600008250404000", "0152103500000000"); // pion nsig max mom 2.00 GeV/c
    cuts.AddCutPCM("00010113", "0d200009267100008250404000", "0152103500000000"); // pion nsig max mom 5.00 GeV/c
  } else if (trainConfig == 712){
    cuts.AddCutPCM("00010113", "0d200009267300003250404000", "0152103500000000"); // qT max 0.05 1D
    cuts.AddCutPCM("00010113", "0d200009267300002250404000", "0152103500000000"); // qT max 0.06 2D
    cuts.AddCutPCM("00010113", "0d200009267300009250404000", "0152103500000000"); // qT max 0.03 2D
    cuts.AddCutPCM("00010113", "0d200009267300008210404000", "0152103500000000"); // Psi pair 0.1  1D
    cuts.AddCutPCM("00010113", "0d200009267300008260404000", "0152103500000000"); // Psi pair 0.05 2D
    cuts.AddCutPCM("00010113", "0d200009267300008280404000", "0152103500000000"); // Psi pair 0.2  2D
  } else if (trainConfig == 714) { // RBins studies, ITS with large weights spplited
    cuts.AddCutPCM("00010113", "0dh00009267300008250404000", "0152103500000000"); // eta < 0.8   5-13 cm
    cuts.AddCutPCM("00010113", "0di00009267300008250404000", "0152103500000000"); // eta < 0.8  13-33.5 cm
    cuts.AddCutPCM("00010113", "0dl00009267300008250404000", "0152103500000000"); // eta < 0.8  72-95 cm
  } else if (trainConfig == 715) { // scan of meson alpha cut
    cuts.AddCutPCM("00010113", "0d200009267300008250404000", "0152100500000000"); // alpha <0.7
    cuts.AddCutPCM("00010113", "0d200009267300008250404000", "0152108500000000"); // alpha <0.6
    cuts.AddCutPCM("00010113", "0d200009267300008250404000", "015210g500000000"); // alpha <0.5
    cuts.AddCutPCM("00010113", "0d200009267300008250404000", "015210f500000000"); // alpha <0.4
  } else if (trainConfig == 716) {
    cuts.AddCutPCM("00010113", "0d200009267300008250404000", "015210e500000000"); // alpha <0.3
    cuts.AddCutPCM("00010113", "0d200009267300008250404000", "015210a500000000"); // alpha <0.2
    cuts.AddCutPCM("00010113", "0d200009267300008250404000", "015210d500000000"); // alpha <0.1
  } else if (trainConfig == 717) { // RBins studies,
    cuts.AddCutPCM("00010113", "0dj00009267300008250404000", "0152103500000000"); // eta < 0.8  33.5-55 cm
    cuts.AddCutPCM("00010113", "0dk00009267300008250404000", "0152103500000000"); // eta < 0.8  55-72 cm
    cuts.AddCutPCM("00010113", "0dg00009267300008250404000", "0152103500000000"); // eta < 0.8  95-180 cm
  } else if (trainConfig == 718) { // R 5-180 and remove r bin 55-72
    cuts.AddCutPCM("00010113", "0dm00009267300008250404000", "0152103500000000"); // eta < 0.8
    cuts.AddCutPCM("00010113", "0d200009327000008250404000", "0152103500000000"); // eta < 0.8
  } else if (trainConfig == 719) { // R 5-180 and remove r bin 55-72
    cuts.AddCutPCM("00010113", "0d200009f9730000dge0404000", "0152101500000000"); // eta < 0.8  // Test alpha meson pT dependent
    cuts.AddCutPCM("00010113", "0dm00009f9730000dge0404000", "0152103500000000"); // eta < 0.8  // remove  55-72 bin
    cuts.AddCutPCM("00010113", "0dd00009f9730000dge0404000", "0152103500000000"); // eta < 0.8  // use 5-55 bin only
  } else if (trainConfig == 720) { // R 5-180
    cuts.AddCutPCM("00010113", "0d200009f9730000dge0404000", "0152103500000000"); // eta < 0.8  // Test improved cuts
  } else if (trainConfig == 721) { // R 5-180
    cuts.AddCutPCM("00010113", "0da00009f9730000dge0404000", "0152103500000000"); // eta < 0.8  // Test improved cuts
    cuts.AddCutPCM("00010113", "0db00009f9730000dge0404000", "0152103500000000"); // eta < 0.8  // Test improved cuts
    cuts.AddCutPCM("00010113", "0dc00009f9730000dge0404000", "0152103500000000"); // eta < 0.8  // Test improved cuts
  } else if (trainConfig == 722) { // R 5-180
    cuts.AddCutPCM("00010113", "0dh00009f9730000dge0404000", "0152103500000000"); // eta < 0.8  // Test improved cuts
    cuts.AddCutPCM("00010113", "0di00009f9730000dge0404000", "0152103500000000"); // eta < 0.8  // Test improved cuts
    cuts.AddCutPCM("00010113", "0dj00009f9730000dge0404000", "0152103500000000"); // eta < 0.8  // Test improved cuts
  } else if (trainConfig == 723) { // R 5-180
    cuts.AddCutPCM("00010113", "0dk00009f9730000dge0404000", "0152103500000000"); // eta < 0.8  // Test improved cuts
    cuts.AddCutPCM("00010113", "0dl00009f9730000dge0404000", "0152103500000000"); // eta < 0.8  // Test improved cuts
    cuts.AddCutPCM("00010113", "0dg00009f9730000dge0404000", "0152103500000000"); // eta < 0.8  // Test improved cuts
  } else if (trainConfig == 724) { // R 5-180  // Test smearing parameters. Changed smaring cut 08.11.2019
    cuts.AddCutPCM("00010113", "0d200009f9730000dge0404000", "0152103500900000"); // eta < 0.8  // Test improved cuts
    cuts.AddCutPCM("00010113", "0d200009f9730000dge0404000", "0152103500m00000"); // eta < 0.8  // Test improved cuts
    cuts.AddCutPCM("00010113", "0d200009f9730000dge0404000", "0152103500n00000"); // eta < 0.8  // Test improved cuts
  } else if (trainConfig == 725) { // R 5-180  // Test smearing parameters -1
    cuts.AddCutPCM("00010113", "0d200009f9730000dge0404000", "0152103500d00000"); // eta < 0.8  // Test improved cuts
    cuts.AddCutPCM("00010113", "0d200009f9730000dge0404000", "0152103500e00000"); // eta < 0.8  // Test improved cuts
    cuts.AddCutPCM("00010113", "0d200009f9730000dge0404000", "0152103500f00000"); // eta < 0.8  // Test improved cuts
    cuts.AddCutPCM("00010113", "0d200009f9730000dge0404000", "0152103500g00000"); // eta < 0.8  // Test improved cuts
  } else if (trainConfig == 726) { // R 5-180  // Test smearing parameters -1
    cuts.AddCutPCM("00010113", "0d200009f9730000dge0404000", "0152103500h00000"); // eta < 0.8  // Test improved cuts
    cuts.AddCutPCM("00010113", "0d200009f9730000dge0404000", "0152103500i00000"); // eta < 0.8  // Test improved cuts
    cuts.AddCutPCM("00010113", "0d200009f9730000dge0404000", "0152103500j00000"); // eta < 0.8  // Test improved cuts
  } else if (trainConfig == 727) { // R 5-180  // Test smearing parameters -1
    cuts.AddCutPCM("00010113", "0d200009f9730000dge0404000", "0152103500k00000"); // eta < 0.8  // Test improved cuts
    cuts.AddCutPCM("00010113", "0d200009f9730000dge0404000", "0152103500l00000"); // eta < 0.8  // Test improved cuts
    cuts.AddCutPCM("00010113", "0d200009f9730000dge0404000", "0152103500d00000"); // eta < 0.8  // Test improved cuts
  } else if (trainConfig == 728) { // R 5-180  // Cat 1, cat 2+3
    cuts.AddCutPCM("00010113", "0d200009f9730000dge0424000", "0152103500000000"); // eta < 0.8  // Test improved cuts
    cuts.AddCutPCM("00010113", "0d200009f9730000dge0454000", "0152103500000000"); // eta < 0.8  // Test improved cuts
    cuts.AddCutPCM("00010113", "0dm00009f9730000dge0424000", "0152103500000000"); // eta < 0.8  // Test improved cuts
    cuts.AddCutPCM("00010113", "0dm00009f9730000dge0454000", "0152103500000000"); // eta < 0.8  // Test improved cuts
    cuts.AddCutPCM("00010113", "0dm00009f9730000dge0404000", "0152103520000000"); // eta < 0.8  // Test improved cuts
  } else if (trainConfig == 729) { // R 5-180  // Multiplicity selections
    cuts.AddCutPCM("n0110113", "0dm00009f9730000dge0404000", "0152103500000000"); // eta < 0.8  // Test improved cuts
    cuts.AddCutPCM("n1210113", "0dm00009f9730000dge0404000", "0152103500000000"); // eta < 0.8  // Test improved cuts
    cuts.AddCutPCM("n2510113", "0dm00009f9730000dge0404000", "0152103500000000"); // eta < 0.8  // Test improved cuts
    cuts.AddCutPCM("n5a10113", "0dm00009f9730000dge0404000", "0152103500000000"); // eta < 0.8  // Test improved cuts
  } else if (trainConfig == 730) { // R 5-180  // Multiplicity selections
    cuts.AddCutPCM("m0110113", "0dm00009f9730000dge0404000", "0152103500000000"); // eta < 0.8  // Test improved cuts
    cuts.AddCutPCM("m1510113", "0dm00009f9730000dge0404000", "0152103500000000"); // eta < 0.8  // Test improved cuts
    cuts.AddCutPCM("m5a10113", "0dm00009f9730000dge0404000", "0152103500000000"); // eta < 0.8  // Test improved cuts
  } else if (trainConfig == 731) { // R 5-180  // Sphericity selections
    cuts.AddCutPCM("h0510113", "0dm00009f9730000dge0404000", "0152103500000000"); // eta < 0.8  // Test improved cuts
    cuts.AddCutPCM("h5a10113", "0dm00009f9730000dge0404000", "0152103500000000"); // eta < 0.8  // Test improved cuts
    cuts.AddCutPCM("h0a10113", "0dm00009f9730000dge0404000", "0152103500000000"); // eta < 0.8  // Test improved cuts

    // config like 70X but to be used with weights +50
  } else if (trainConfig == 752) { // as iConfig 702 to be used with MBW
    cuts.AddCutPCM("00010113", "00a00009267300008250404000", "0152103500000000"); //
    cuts.AddCutPCM("00010113", "00b00009267300008250404000", "0152103500000000"); //
    cuts.AddCutPCM("00010113", "00c00009267300008250404000", "0152103500000000"); //
  } else if (trainConfig == 753) {
    cuts.AddCutPCM("00010113", "00200009267300008250404000", "0152103500000000"); // Min Bias with photon asym and dedx at high pT
  } else if (trainConfig == 754) {
    cuts.AddCutPCM("00010113", "00200049267300008250404000", "0152103500000000"); // Min Bias with photon asym and dedx at high pT, pt 0.075
    cuts.AddCutPCM("00010113", "00200019267300008250404000", "0152103500000000"); // Min Bias with photon asym and dedx at high pT, pT 0.1
    cuts.AddCutPCM("00010113", "00200059267300008250404000", "0152103500000000"); // Min Bias with photon asym and dedx at high pT, pT 0.125
  } else if (trainConfig == 755) {
    cuts.AddCutPCM("00010113", "002000d9267300008250404000", "0152103500000000"); // Min Bias with photon asym and dedx at high pT, pT 0.06
    cuts.AddCutPCM("00010113", "002000m9267300008250404000", "0152103500000000"); // Min Bias with photon asym and dedx at high pT, pT 0.08
    cuts.AddCutPCM("00010113", "002000n9267300008250404000", "0152103500000000"); // Min Bias with photon asym and dedx at high pT, pT 0.09

  } else if (trainConfig == 756) { // asymetry cut removed from configs 702-703-704  and 753-754-755 on 14.06.2018 (cuts low pT)
    cuts.AddCutPCM("00010113", "0d200009267300008250404000", "0152103500000000"); // eta < 0.8
  } else if (trainConfig == 757) { // asymetry cut removed from configs 702-703-704  and 753-754-755 on 14.06.2018 (cuts low pT)
    cuts.AddCutPCM("00010113", "0da00009267300008250404000", "0152103500000000"); // eta < 0.8
    cuts.AddCutPCM("00010113", "0db00009267300008250404000", "0152103500000000"); // eta < 0.8
    cuts.AddCutPCM("00010113", "0dc00009267300008250404000", "0152103500000000"); // eta < 0.8
    // Configs for systematics
  } else if (trainConfig == 758) {    // Systematics of pT
    cuts.AddCutPCM("00010113", "0d200009267300008250404000", "0152103500000000"); // eta < 0.8 standard for 13 TeV
    cuts.AddCutPCM("00010113", "0d200079267300008250404000", "0152103500000000"); // min pT no cut
    cuts.AddCutPCM("00010113", "0d200069267300008250404000", "0152103500000000"); // min pT 40 MeV
    cuts.AddCutPCM("00010113", "0d200049267300008250404000", "0152103500000000"); // min pT 75 MeV
    cuts.AddCutPCM("00010113", "0d200019267300008250404000", "0152103500000000"); // min pT 100 MeV
  } else if (trainConfig == 759) {
    cuts.AddCutPCM("00010113", "0d200008267300008250404000", "0152103500000000"); // TPC cluster 35%
    cuts.AddCutPCM("00010113", "0d200006267300008250404000", "0152103500000000"); // TPC cluster 70%
  } else if (trainConfig == 760) {
    cuts.AddCutPCM("00010113", "0d200009367300008250404000", "0152103500000000"); // edEdx -4,5
    cuts.AddCutPCM("00010113", "0d200009667300008250404000", "0152103500000000"); // edEdx -2.5,4
    cuts.AddCutPCM("00010113", "0d200009257300008250404000", "0152103500000000"); // pidEdx 2,-10
    cuts.AddCutPCM("00010113", "0d200009217300008250404000", "0152103500000000"); // pidEdx 0,-10
  } else if (trainConfig == 761) {
    cuts.AddCutPCM("00010113", "0d200009260300008250404000", "0152103500000000"); // pion nsig min mom 0.50 GeV/c
    cuts.AddCutPCM("00010113", "0d200009266300008250404000", "0152103500000000"); // pion nsig min mom 0.25 GeV/c
    cuts.AddCutPCM("00010113", "0d200009267600008250404000", "0152103500000000"); // pion nsig max mom 2.00 GeV/c
    cuts.AddCutPCM("00010113", "0d200009267100008250404000", "0152103500000000"); // pion nsig max mom 5.00 GeV/c
  } else if (trainConfig == 762){
    cuts.AddCutPCM("00010113", "0d200009267300003250404000", "0152103500000000"); // qT max 0.05 1D
    cuts.AddCutPCM("00010113", "0d200009267300002250404000", "0152103500000000"); // qT max 0.06 2D
    cuts.AddCutPCM("00010113", "0d200009267300009250404000", "0152103500000000"); // qT max 0.03 2D
    cuts.AddCutPCM("00010113", "0d200009267300008210404000", "0152103500000000"); // Psi pair 0.1  1D
    cuts.AddCutPCM("00010113", "0d200009267300008260404000", "0152103500000000"); // Psi pair 0.05 2D
    cuts.AddCutPCM("00010113", "0d200009267300008280404000", "0152103500000000"); // Psi pair 0.2  2D
  } else if (trainConfig == 764) { // RBins studies, ITS with large weights spplited
    cuts.AddCutPCM("00010113", "0dh00009267300008250404000", "0152103500000000"); // eta < 0.8   5-13 cm
    cuts.AddCutPCM("00010113", "0di00009267300008250404000", "0152103500000000"); // eta < 0.8  13-33.5
    cuts.AddCutPCM("00010113", "0dl00009267300008250404000", "0152103500000000"); // eta < 0.8  72-95
  } else if (trainConfig == 765) { // scan of meson alpha cut
    cuts.AddCutPCM("00010113", "0d200009267300008250404000", "0152100500000000"); // alpha <0.7
    cuts.AddCutPCM("00010113", "0d200009267300008250404000", "0152108500000000"); // alpha <0.6
    cuts.AddCutPCM("00010113", "0d200009267300008250404000", "015210g500000000"); // alpha <0.5
    cuts.AddCutPCM("00010113", "0d200009267300008250404000", "015210f500000000"); // alpha <0.4
  } else if (trainConfig == 766) {
    cuts.AddCutPCM("00010113", "0d200009267300008250404000", "015210e500000000"); // alpha <0.3
    cuts.AddCutPCM("00010113", "0d200009267300008250404000", "015210a500000000"); // alpha <0.2
    cuts.AddCutPCM("00010113", "0d200009267300008250404000", "015210d500000000"); // alpha <0.1
  } else if (trainConfig == 767) { // RBins studies,
    cuts.AddCutPCM("00010113", "0dj00009267300008250404000", "0152103500000000"); // eta < 0.8  33.5-55 cm
    cuts.AddCutPCM("00010113", "0dk00009267300008250404000", "0152103500000000"); // eta < 0.8  55-72 cm
    cuts.AddCutPCM("00010113", "0dg00009267300008250404000", "0152103500000000"); // eta < 0.8  95-180 cm
  } else if (trainConfig == 768) { // R 5-180 and remove r bin 55-72
    cuts.AddCutPCM("00010113", "0dm00009267300008250404000", "0152103500000000"); // eta < 0.8
    cuts.AddCutPCM("00010113", "0d200009327000008250404000", "0152103500000000"); // eta < 0.8
  } else if (trainConfig == 769) { // R 5-180 and remove r bin 55-72
    cuts.AddCutPCM("00010113", "0d200009f9730000dge0404000", "0152101500000000"); // eta < 0.8  // Test alpha meson pT dependent
    cuts.AddCutPCM("00010113", "0dm00009f9730000dge0404000", "0152103500000000"); // eta < 0.8  // remove  55-72 bin
    cuts.AddCutPCM("00010113", "0dd00009f9730000dge0404000", "0152103500000000"); // eta < 0.8  // use 5-55 bin only
  } else if (trainConfig == 770) { // R 5-180
    cuts.AddCutPCM("00010113", "0d200009f9730000dge0404000", "0152103500000000"); // eta < 0.8  // Test improved cuts
  } else if (trainConfig == 771) { // R 5-180
    cuts.AddCutPCM("00010113", "0da00009f9730000dge0404000", "0152103500000000"); // eta < 0.8  // Test improved cuts
    cuts.AddCutPCM("00010113", "0db00009f9730000dge0404000", "0152103500000000"); // eta < 0.8  // Test improved cuts
    cuts.AddCutPCM("00010113", "0dc00009f9730000dge0404000", "0152103500000000"); // eta < 0.8  // Test improved cuts

  } else if (trainConfig == 772) { // R 5-180
    cuts.AddCutPCM("00010113", "0dh00009f9730000dge0404000", "0152103500000000"); // eta < 0.8  // Test improved cuts
    cuts.AddCutPCM("00010113", "0di00009f9730000dge0404000", "0152103500000000"); // eta < 0.8  // Test improved cuts
    cuts.AddCutPCM("00010113", "0dj00009f9730000dge0404000", "0152103500000000"); // eta < 0.8  // Test improved cuts

  } else if (trainConfig == 773) { // R 5-180
    cuts.AddCutPCM("00010113", "0dk00009f9730000dge0404000", "0152103500000000"); // eta < 0.8  // Test improved cuts
    cuts.AddCutPCM("00010113", "0dl00009f9730000dge0404000", "0152103500000000"); // eta < 0.8  // Test improved cuts
    cuts.AddCutPCM("00010113", "0dg00009f9730000dge0404000", "0152103500000000"); // eta < 0.8  // Test improved cuts


  } else if (trainConfig == 778) { // R 5-180  // Cat 1, cat 2+3
    cuts.AddCutPCM("00010113", "0d200009f9730000dge0424000", "0152103500000000"); // eta < 0.8  // Test improved cuts
    cuts.AddCutPCM("00010113", "0d200009f9730000dge0454000", "0152103500000000"); // eta < 0.8  // Test improved cuts
    cuts.AddCutPCM("00010113", "0dm00009f9730000dge0424000", "0152103500000000"); // eta < 0.8  // Test improved cuts
    cuts.AddCutPCM("00010113", "0dm00009f9730000dge0454000", "0152103500000000"); // eta < 0.8  // Test improved cuts
    cuts.AddCutPCM("00010113", "0dm00009f9730000dge0404000", "0152103520000000"); // eta < 0.8  // Test improved cuts
  } else if (trainConfig == 779) { // R 5-180  // Multiplicity selections
    cuts.AddCutPCM("n0110113", "0dm00009f9730000dge0404000", "0152103500000000"); // eta < 0.8  // Test improved cuts
    cuts.AddCutPCM("n1210113", "0dm00009f9730000dge0404000", "0152103500000000"); // eta < 0.8  // Test improved cuts
    cuts.AddCutPCM("n2510113", "0dm00009f9730000dge0404000", "0152103500000000"); // eta < 0.8  // Test improved cuts
    cuts.AddCutPCM("n5a10113", "0dm00009f9730000dge0404000", "0152103500000000"); // eta < 0.8  // Test improved cuts
  } else if (trainConfig == 780) { // R 5-180  // Sphericity selections
    cuts.AddCutPCM("m0110113", "0dm00009f9730000dge0404000", "0152103500000000"); // eta < 0.8  // Test improved cuts
    cuts.AddCutPCM("m1510113", "0dm00009f9730000dge0404000", "0152103500000000"); // eta < 0.8  // Test improved cuts
    cuts.AddCutPCM("m5a10113", "0dm00009f9730000dge0404000", "0152103500000000"); // eta < 0.8  // Test improved cuts
  } else if (trainConfig == 781) { // R 5-180  // Sphericity selections
    cuts.AddCutPCM("h0510113", "0dm00009f9730000dge0404000", "0152103500000000"); // eta < 0.8  // Test improved cuts
    cuts.AddCutPCM("h5a10113", "0dm00009f9730000dge0404000", "0152103500000000"); // eta < 0.8  // Test improved cuts
    cuts.AddCutPCM("h0a10113", "0dm00009f9730000dge0404000", "0152103500000000"); // eta < 0.8  // Test improved cuts


    //-----------------same as 7XX to be used with MBW extracted from 5TeV  


  } else if (trainConfig == 819) { // R 5-180 and remove r bin 55-72
    cuts.AddCutPCM("00010113", "0d200009f9730000dge0404000", "0152101500000000"); // eta < 0.8  // Test alpha meson pT dependent
    cuts.AddCutPCM("00010113", "0dm00009f9730000dge0404000", "0152103500000000"); // eta < 0.8  // remove  55-72 bin
    cuts.AddCutPCM("00010113", "0dd00009f9730000dge0404000", "0152103500000000"); // eta < 0.8  // use 5-55 bin only
  } else if (trainConfig == 820) { // R 5-180
    cuts.AddCutPCM("00010113", "0d200009f9730000dge0404000", "0152103500000000"); // eta < 0.8  // Test improved cuts
  } else if (trainConfig == 821) { // R 5-180
    cuts.AddCutPCM("00010113", "0da00009f9730000dge0404000", "0152103500000000"); // eta < 0.8  // Test improved cuts
    cuts.AddCutPCM("00010113", "0db00009f9730000dge0404000", "0152103500000000"); // eta < 0.8  // Test improved cuts
    cuts.AddCutPCM("00010113", "0dc00009f9730000dge0404000", "0152103500000000"); // eta < 0.8  // Test improved cuts
  } else if (trainConfig == 822) { // R 5-180
    cuts.AddCutPCM("00010113", "0dh00009f9730000dge0404000", "0152103500000000"); // eta < 0.8  // Test improved cuts
    cuts.AddCutPCM("00010113", "0di00009f9730000dge0404000", "0152103500000000"); // eta < 0.8  // Test improved cuts
    cuts.AddCutPCM("00010113", "0dj00009f9730000dge0404000", "0152103500000000"); // eta < 0.8  // Test improved cuts
  } else if (trainConfig == 823) { // R 5-180
    cuts.AddCutPCM("00010113", "0dk00009f9730000dge0404000", "0152103500000000"); // eta < 0.8  // Test improved cuts
    cuts.AddCutPCM("00010113", "0dl00009f9730000dge0404000", "0152103500000000"); // eta < 0.8  // Test improved cuts
    cuts.AddCutPCM("00010113", "0dg00009f9730000dge0404000", "0152103500000000"); // eta < 0.8  // Test improved cuts
  } else if (trainConfig == 828) { // R 5-180  // Cat 1, cat 2+3   Meson Cat >=2
    cuts.AddCutPCM("00010113", "0d200009f9730000dge0424000", "0152103500000000"); // eta < 0.8  // Test improved cuts
    cuts.AddCutPCM("00010113", "0d200009f9730000dge0454000", "0152103500000000"); // eta < 0.8  // Test improved cuts
    cuts.AddCutPCM("00010113", "0dm00009f9730000dge0424000", "0152103500000000"); // eta < 0.8  // Test improved cuts
    cuts.AddCutPCM("00010113", "0dm00009f9730000dge0454000", "0152103500000000"); // eta < 0.8  // Test improved cuts
    cuts.AddCutPCM("00010113", "0dm00009f9730000dge0404000", "0152103520000000"); // eta < 0.8  // Test improved cuts

    //-----------------same as 7XX to be used with MBW extracted from 5TeV Nch  


  } else if (trainConfig == 869) { // R 5-180 and remove r bin 55-72
    cuts.AddCutPCM("00010113", "0d200009f9730000dge0404000", "0152101500000000"); // eta < 0.8  // Test alpha meson pT dependent
    cuts.AddCutPCM("00010113", "0dm00009f9730000dge0404000", "0152103500000000"); // eta < 0.8  // remove  55-72 bin
    cuts.AddCutPCM("00010113", "0dd00009f9730000dge0404000", "0152103500000000"); // eta < 0.8  // use 5-55 bin only
  } else if (trainConfig == 870) { // R 5-180
    cuts.AddCutPCM("00010113", "0d200009f9730000dge0404000", "0152103500000000"); // eta < 0.8  // Test improved cuts
  } else if (trainConfig == 871) { // R 5-180
    cuts.AddCutPCM("00010113", "0da00009f9730000dge0404000", "0152103500000000"); // eta < 0.8  // Test improved cuts
    cuts.AddCutPCM("00010113", "0db00009f9730000dge0404000", "0152103500000000"); // eta < 0.8  // Test improved cuts
    cuts.AddCutPCM("00010113", "0dc00009f9730000dge0404000", "0152103500000000"); // eta < 0.8  // Test improved cuts
  } else if (trainConfig == 872) { // R 5-180
    cuts.AddCutPCM("00010113", "0dh00009f9730000dge0404000", "0152103500000000"); // eta < 0.8  // Test improved cuts
    cuts.AddCutPCM("00010113", "0di00009f9730000dge0404000", "0152103500000000"); // eta < 0.8  // Test improved cuts
    cuts.AddCutPCM("00010113", "0dj00009f9730000dge0404000", "0152103500000000"); // eta < 0.8  // Test improved cuts
  } else if (trainConfig == 873) { // R 5-180
    cuts.AddCutPCM("00010113", "0dk00009f9730000dge0404000", "0152103500000000"); // eta < 0.8  // Test improved cuts
    cuts.AddCutPCM("00010113", "0dl00009f9730000dge0404000", "0152103500000000"); // eta < 0.8  // Test improved cuts
    cuts.AddCutPCM("00010113", "0dg00009f9730000dge0404000", "0152103500000000"); // eta < 0.8  // Test improved cuts
  } else if (trainConfig == 878) { // R 5-180  // Cat 1, cat 2+3, Meson Cat >= 2
    cuts.AddCutPCM("00010113", "0d200009f9730000dge0424000", "0152103500000000"); // eta < 0.8  // Test improved cuts
    cuts.AddCutPCM("00010113", "0d200009f9730000dge0454000", "0152103500000000"); // eta < 0.8  // Test improved cuts
    cuts.AddCutPCM("00010113", "0dm00009f9730000dge0424000", "0152103500000000"); // eta < 0.8  // Test improved cuts
    cuts.AddCutPCM("00010113", "0dm00009f9730000dge0454000", "0152103500000000"); // eta < 0.8  // Test improved cuts
    cuts.AddCutPCM("00010113", "0dm00009f9730000dge0404000", "0152103520000000"); // eta < 0.8  // Test improved cuts
  } else if (trainConfig == 879) { // R 5-180  // Multiplicity selections
    cuts.AddCutPCM("n0110113", "0dm00009f9730000dge0404000", "0152103500000000"); // eta < 0.8  // Test improved cuts
    cuts.AddCutPCM("n1210113", "0dm00009f9730000dge0404000", "0152103500000000"); // eta < 0.8  // Test improved cuts
    cuts.AddCutPCM("n2510113", "0dm00009f9730000dge0404000", "0152103500000000"); // eta < 0.8  // Test improved cuts
    cuts.AddCutPCM("n5a10113", "0dm00009f9730000dge0404000", "0152103500000000"); // eta < 0.8  // Test improved cuts
  } else if (trainConfig == 880) { // R 5-180  // Multiplicity selections
    cuts.AddCutPCM("m0110113", "0dm00009f9730000dge0404000", "0152103500000000"); // eta < 0.8  // Test improved cuts
    cuts.AddCutPCM("m1510113", "0dm00009f9730000dge0404000", "0152103500000000"); // eta < 0.8  // Test improved cuts
    cuts.AddCutPCM("m5a10113", "0dm00009f9730000dge0404000", "0152103500000000"); // eta < 0.8  // Test improved cuts


    // systematics for 13 TeV
    // defauul cut  AM 
    //  .AddCutPCM("00010113", "0dm00009f9730000dge0404000", "0152103500000000"); // eta < 0.8  // Test improved cuts

  } else if (trainConfig == 881) {   // min pT variations
    cuts.AddCutPCM("00010113", "0dm00069f9730000dge0404000", "0152103500000000"); // eta < 0.8  // remove  55-72 bin, min pT 40 MeV
    cuts.AddCutPCM("00010113", "0dm00049f9730000dge0404000", "0152103500000000"); // eta < 0.8  // remove  55-72 bin, min pT 75 MeV
    cuts.AddCutPCM("00010113", "0dm00019f9730000dge0404000", "0152103500000000"); // eta < 0.8  // remove  55-72 bin, min pT 100MeV

  } else if (trainConfig == 882) {   // TPC clusters, cosPA
    cuts.AddCutPCM("00010113", "0dm00008f9730000dge0404000", "0152103500000000"); // TPC cluster 35%
    cuts.AddCutPCM("00010113", "0dm00006f9730000dge0404000", "0152103500000000"); // TPC cluster 70%
    cuts.AddCutPCM("00010113", "0dm00009f9730000dge0604000", "0152103500000000"); // cosPA 0.9
    cuts.AddCutPCM("00010113", "0dm00009f9730000dge0304000", "0152103500000000"); // cosPA 0.75

  } else if (trainConfig == 883) {   // TPC clusters, cosPA
    cuts.AddCutPCM("00010113", "0dm0000939730000dge0404000", "0152103500000000"); // nsig electron   -4,5
    cuts.AddCutPCM("00010113", "0dm0000969730000dge0404000", "0152103500000000"); // nsig electron -2.5,4
    cuts.AddCutPCM("00010113", "0dm00009f5730000dge0404000", "0152103500000000"); // nsig pion 2,-10
    cuts.AddCutPCM("00010113", "0dm00009f1730000dge0404000", "0152103500000000"); // nsig pion 0,-10

  } else if (trainConfig == 884) {
    cuts.AddCutPCM("00010113", "0dm00009f9030000dge0404000", "0152103500000000"); // pion nsig min mom 0.50 GeV/c
    cuts.AddCutPCM("00010113", "0dm00009f9630000dge0404000", "0152103500000000"); // pion nsig min mom 0.25 GeV/c
    cuts.AddCutPCM("00010113", "0dm00009f9760000dge0404000", "0152103500000000"); // pion nsig max mom 2.00 GeV/c
    cuts.AddCutPCM("00010113", "0dm00009f9710000dge0404000", "0152103500000000"); // pion nsig max mom 5.00 GeV/c


  } else if (trainConfig == 885) {   // chi2 variations
    cuts.AddCutPCM("00010113", "0dm00009f9730000d1e0404000", "0152103500000000"); // chi2 50
    cuts.AddCutPCM("00010113", "0dm00009f9730000dfe0404000", "0152103500000000"); // chi2 50 chi2 dep -0.065
    cuts.AddCutPCM("00010113", "0dm00009f9730000dhe0404000", "0152103500000000"); // chi2 50 chi2 dep -0.050
    cuts.AddCutPCM("00010113", "0dm00009f9730000dge0400000", "0152103500000000"); // remove reject close v0

  } else if (trainConfig == 886) {   // Psi pair variations
    cuts.AddCutPCM("00010113", "0dm00009f9730000dgd0404000", "0152103500000000"); // Psi pair 0.15 dep
    cuts.AddCutPCM("00010113", "0dm00009f9730000dgf0404000", "0152103500000000"); // Psi pair 0.20 dep
    cuts.AddCutPCM("00010113", "0dm00009f9730000dgg0404000", "0152103500000000"); // Psi pair 0.30 dep
    cuts.AddCutPCM("00010113", "0dm00009227300008250404000", "0152103500000000"); // old cuts (run1)

  } else if (trainConfig == 887) {   
    cuts.AddCutPCM("00010113", "0dm00009f9730000dge0404000", "0252103500000000"); // variation BG scheme track mult
    cuts.AddCutPCM("00010113", "0dm00009f9730000dge0404000", "0r52103500000000"); // variation BG scheme track mult
    cuts.AddCutPCM("00010113", "0dm00009f9730000dge0404000", "0152107500000000"); // alpha meson 0.85
    cuts.AddCutPCM("00010113", "0dm00009f9730000dge0404000", "0152105500000000"); // alpha meson 0.75

  } else if (trainConfig == 888) {   // qT variations
    cuts.AddCutPCM("00010113", "0dm00009f9730000age0404000", "0152103500000000"); // qT max 0.040, qtptmax 0.11
    cuts.AddCutPCM("00010113", "0dm00009f9730000ege0404000", "0152103500000000"); // qT max 0.060, qtptmax 0.14
    cuts.AddCutPCM("00010113", "0dm00009f9730000fge0404000", "0152103500000000"); // qT max 0.070, qtptmax 0.16
    cuts.AddCutPCM("00010113", "0dm00009f9730000dge0404000", "0152101500000000"); // alpha meson pT dependent



    // ---------same as 6XX  low B with MBW extracted from 5TeV

  } else if (trainConfig == 919) { // R 5-180 and remove r bin 55-72
    cuts.AddCutPCM("00010113", "0d200089f9730000iih0404000", "0152101500000000"); // eta < 0.8  // Test alpha meson pT dependent
    cuts.AddCutPCM("00010113", "0dm00089f9730000iih0404000", "0152103500000000"); // eta < 0.8  // remove  55-72 bin
    cuts.AddCutPCM("00010113", "0dd00089f9730000iih0404000", "0152103500000000"); // eta < 0.8  // use 5-55 bin only
  } else if (trainConfig == 920) { // R 5-180
    cuts.AddCutPCM("00010113", "0d200089f9730000iih0404000", "0152103500000000"); // eta < 0.8  // Test improved cuts
  } else if (trainConfig == 921) { // R 5-180
    cuts.AddCutPCM("00010113", "0da00089f9730000iih0404000", "0152103500000000"); // eta < 0.8  // Test improved cuts
    cuts.AddCutPCM("00010113", "0db00089f9730000iih0404000", "0152103500000000"); // eta < 0.8  // Test improved cuts
    cuts.AddCutPCM("00010113", "0dc00089f9730000iih0404000", "0152103500000000"); // eta < 0.8  // Test improved cuts
  } else if (trainConfig == 922) { // R 5-180
    cuts.AddCutPCM("00010113", "0dh00089f9730000iih0404000", "0152103500000000"); // eta < 0.8  // Test improved cuts
    cuts.AddCutPCM("00010113", "0di00089f9730000iih0404000", "0152103500000000"); // eta < 0.8  // Test improved cuts
    cuts.AddCutPCM("00010113", "0dj00089f9730000iih0404000", "0152103500000000"); // eta < 0.8  // Test improved cuts
  } else if (trainConfig == 923) { // R 5-180
    cuts.AddCutPCM("00010113", "0dk00089f9730000iih0404000", "0152103500000000"); // eta < 0.8  // Test improved cuts
    cuts.AddCutPCM("00010113", "0dl00089f9730000iih0404000", "0152103500000000"); // eta < 0.8  // Test improved cuts
    cuts.AddCutPCM("00010113", "0dg00089f9730000iih0404000", "0152103500000000"); // eta < 0.8  // Test improved cuts
  } else if (trainConfig == 928) { // R 5-180  // Cat 1, cat 2+3   Meson Cat >=2
    cuts.AddCutPCM("00010113", "0d200089f9730000iih0424000", "0152103500000000"); // eta < 0.8  // Test improved cuts
    cuts.AddCutPCM("00010113", "0d200089f9730000iih0454000", "0152103500000000"); // eta < 0.8  // Test improved cuts
    cuts.AddCutPCM("00010113", "0dm00089f9730000iih0424000", "0152103500000000"); // eta < 0.8  // Test improved cuts
    cuts.AddCutPCM("00010113", "0dm00089f9730000iih0454000", "0152103500000000"); // eta < 0.8  // Test improved cuts
    cuts.AddCutPCM("00010113", "0dm00089f9730000iih0404000", "0152103520000000"); // eta < 0.8  // Test improved cuts

    //-----------------same as 6XX to be used with MBW extracted from 5TeV Nch  


  } else if (trainConfig == 969) { // R 5-180 and remove r bin 55-72
    cuts.AddCutPCM("00010113", "0d200089f9730000iih0404000", "0152101500000000"); // eta < 0.8  // Test alpha meson pT dependent
    cuts.AddCutPCM("00010113", "0dm00089f9730000iih0404000", "0152103500000000"); // eta < 0.8  // remove  55-72 bin
    cuts.AddCutPCM("00010113", "0dd00089f9730000iih0404000", "0152103500000000"); // eta < 0.8  // use 5-55 bin only
  } else if (trainConfig == 970) { // R 5-180
    cuts.AddCutPCM("00010113", "0d200089f9730000iih0404000", "0152103500000000"); // eta < 0.8  // Test improved cuts
  } else if (trainConfig == 971) { // R 5-180
    cuts.AddCutPCM("00010113", "0da00089f9730000iih0404000", "0152103500000000"); // eta < 0.8  // Test improved cuts
    cuts.AddCutPCM("00010113", "0db00089f9730000iih0404000", "0152103500000000"); // eta < 0.8  // Test improved cuts
    cuts.AddCutPCM("00010113", "0dc00089f9730000iih0404000", "0152103500000000"); // eta < 0.8  // Test improved cuts
  } else if (trainConfig == 972) { // R 5-180
    cuts.AddCutPCM("00010113", "0dh00089f9730000iih0404000", "0152103500000000"); // eta < 0.8  // Test improved cuts
    cuts.AddCutPCM("00010113", "0di00089f9730000iih0404000", "0152103500000000"); // eta < 0.8  // Test improved cuts
    cuts.AddCutPCM("00010113", "0dj00089f9730000iih0404000", "0152103500000000"); // eta < 0.8  // Test improved cuts
  } else if (trainConfig == 973) { // R 5-180
    cuts.AddCutPCM("00010113", "0dk00089f9730000iih0404000", "0152103500000000"); // eta < 0.8  // Test improved cuts
    cuts.AddCutPCM("00010113", "0dl00089f9730000iih0404000", "0152103500000000"); // eta < 0.8  // Test improved cuts
    cuts.AddCutPCM("00010113", "0dg00089f9730000iih0404000", "0152103500000000"); // eta < 0.8  // Test improved cuts
  } else if (trainConfig == 978) { // R 5-180  // Cat 1, cat 2+3, Meson Cat >= 2
    cuts.AddCutPCM("00010113", "0d200089f9730000iih0424000", "0152103500000000"); // eta < 0.8  // Test improved cuts
    cuts.AddCutPCM("00010113", "0d200089f9730000iih0454000", "0152103500000000"); // eta < 0.8  // Test improved cuts
    cuts.AddCutPCM("00010113", "0dm00089f9730000iih0424000", "0152103500000000"); // eta < 0.8  // Test improved cuts
    cuts.AddCutPCM("00010113", "0dm00089f9730000iih0454000", "0152103500000000"); // eta < 0.8  // Test improved cuts
    cuts.AddCutPCM("00010113", "0dm00089f9730000iih0404000", "0152103520000000"); // eta < 0.8  // Test improved cuts



  //----------------------------- configuration for 2.76TeV standard cuts ----------------------------------------------------
  } else if (trainConfig == 1001){
    cuts.AddCutPCM("00000113", "00200009397300008250400000", "0163103100900000"); // new pi0/eta cut 2.76TeV
  } else if (trainConfig == 1002) {
    cuts.AddCutPCM("00000113", "00200009397300008250400000", "0163103100000000"); // new pi0/eta cut 2.76TeV without MC smearing
  } else if (trainConfig == 1003) {
    cuts.AddCutPCM("00000113", "00200009366300003800000000", "0163103100900000"); // standard cut Pi0 pp 2.76TeV PbPb paper 2012
  } else if (trainConfig == 1004) {
    cuts.AddCutPCM("00000113", "00200009297002008250400000", "0163103100900000"); // standard cut LHC11h pp 2.76TeV
  } else if (trainConfig == 1005) {
    cuts.AddCutPCM("00000113", "00200009227302008250404000", "0163101500000000"); // Ana eta analysis prefered 2.76TeV
  } else if (trainConfig == 1006) {
    cuts.AddCutPCM("00000113", "00200009327000008250400000", "0163103100900000"); // go with 1sigma pi rejec to infty 2.76TeV
  } else if (trainConfig == 1007) {
    cuts.AddCutPCM("00000113", "00200009317000008250400000", "0163103100900000"); // go with 0sigma pi rejec to infty 2.76TeV
  } else if (trainConfig == 1008) {
    cuts.AddCutPCM("00000113", "00200009357000008250400000", "0163103100900000"); // go with 2sigma pi reject to infy 2.76TeV
  } else if (trainConfig == 1009) {
    cuts.AddCutPCM("00003113", "00200009397300008250400000", "0163103100900000"); // new pi0/eta cut 2.76TeV wSDD
  } else if (trainConfig == 1010) {
    cuts.AddCutPCM("00051013", "00200009397300008250400000", "0163103100900000"); // new pi0/eta cut 2.76TeV wSDD & EMC1

  //----------------------------- configuration for  8 TeV standard  --------------------------------------------------------
  } else if (trainConfig == 1020) {
    cuts.AddCutPCM("00010113", "00200009227300008250404000", "0152103500000000"); //standard cut pp 8 TeV
  } else if (trainConfig == 1021) {
    cuts.AddCutPCM("00052113", "00200009227302008250400000", "0152103500000000"); //standard cut pp 8 TeV EMC7
  } else if (trainConfig == 1022) {
    cuts.AddCutPCM("00081113", "00200009227302008250400000", "0152103500000000"); //standard cut pp 8 TeV EGA
  } else if (trainConfig == 1023) {
    cuts.AddCutPCM("00010213", "00200009227300008250404000", "0152103500000000"); //standard cut pp 8 TeV + past future max rejection
  } else if (trainConfig == 1024) {
    cuts.AddCutPCM("00010513", "00200009227300008250404000", "0152103500000000"); //standard cut pp 8 TeV + past future medium rejection
  } else if (trainConfig == 1025) {
    cuts.AddCutPCM("00010113", "0a200009227300008250404000", "0152103500000000"); //eta cut 0.2 < |eta| < 0.9
  } else if (trainConfig == 1026) {
    cuts.AddCutPCM("00010613", "00200009227300008250404000", "0152103500000000"); //V0M vs TPCout cut6
  } else if (trainConfig == 1027) {
    cuts.AddCutPCM("00010113", "002000j9227300008250404000", "0152103500000000"); //asym pT cut: 0.100 GeV and 0.075 GeV
  } else if (trainConfig == 1028) {
    cuts.AddCutPCM("00010113", "002000l9227300008250404000", "0152103500000000"); //asym pT cut: 0.200 GeV and 0.075 GeV
  } else if (trainConfig == 1029) {
    cuts.AddCutPCM("00010113", "0dm00009f9730000dge0404000", "0152103500000000"); // new default (R region rej. + eta<0.8 + DC)
  //----------------------------- configuration for  7 TeV standard cuts -----------------------------------------------------
  } else if (trainConfig == 1030) {
    cuts.AddCutPCM("00000113", "00200009227300008250404000", "0152103500000000"); //New standard cut pp 7 TeV direct photon
  } else if (trainConfig == 1031) {
    cuts.AddCutPCM("00000113", "00200009227302008250400000", "0152103500000000"); //standard cut pp 7 TeV
  } else if (trainConfig == 1032) {
    cuts.AddCutPCM("00000113", "00200008366300000200000000", "0163103100900000"); //old standard cut pp 7 TeV
  } else if (trainConfig == 1033) {
    cuts.AddCutPCM("00000113", "002000c9227300008250404000", "0152103500000000"); // 7 TeV std, but min electron pT > 0.6 for all configs
  } else if (trainConfig == 1034) {
    cuts.AddCutPCM("00000113", "002000b9227300008250404000", "0152103500000000"); // gamma pT > 0.1 GeV/c
  } else if (trainConfig == 1035) {
    cuts.AddCutPCM("00000113", "002000e9227300008250404000", "0152103500000000"); // gamma pT > 0.15 GeV/c
  } else if (trainConfig == 1036) {
    cuts.AddCutPCM("00000113", "002000f9227300008250404000", "0152103500000000"); // gamma pT > 0.2 GeV/c
  } else if (trainConfig == 1037) {
    cuts.AddCutPCM("00000113", "0dm0000922700000dge0404000", "0152103500000000"); // new cuts consistent with omega analysis

  } else if (trainConfig == 1039) { // T0-based MB trigger
    cuts.AddCutPCM("00011113", "0dm00009f9730000dge0404000", "0152103500000000"); // new default (R region rej. + eta<0.8 + DC)

  //----------------------------- configuration for run 2 analysis 13 TeV ----------------------------------------------------
  } else if (trainConfig == 1040){
    cuts.AddCutPCM("00010113", "00200009227302008254404000", "0152101500000000"); //standard cut Gamma pp 13TeV, V0AND
  } else if (trainConfig == 1041){
    cuts.AddCutPCM("00074113", "00200009227302008254404000", "0152101500000000"); //standard cut Gamma pp 13TeV, V0 HM
  } else if (trainConfig == 1042){
    cuts.AddCutPCM("00075113", "00200009227302008254404000", "0152101500000000"); //standard cut Gamma pp 13TeV, SPD HM
  } else if (trainConfig == 1043){
    cuts.AddCutPCM("00010113", "00200009227300008250404000", "0152103500000000"); //New standard cut Gamma Pi0 Eta pp 13TeV, V0AND
  } else if (trainConfig == 1044){
    cuts.AddCutPCM("00010113", "00200009266300008854404000", "0152101500000000"); // A. Marin alpha pT dependent and gamma asym cut
  } else if (trainConfig == 1045){
    cuts.AddCutPCM("00010113", "00200009267300008254404000", "0152103500000000"); // A. Marin alpha pT dependent and gamma asym cut
  } else if (trainConfig == 1046){
    cuts.AddCutPCM("00010113", "00a00009267300008254404000", "0152103500000000"); // A. Marin alpha pT dependent and gamma asym cut
  } else if (trainConfig == 1047){
    cuts.AddCutPCM("00010113", "00b00009267300008254404000", "0152103500000000"); // A. Marin alpha pT dependent and gamma asym cut
  } else if (trainConfig == 1048){
    cuts.AddCutPCM("00010113", "00c00009267300008254404000", "0152103500000000"); // A. Marin alpha pT dependent and gamma asym cut
  } else if (trainConfig == 1049){
    cuts.AddCutPCM("00010113", "00200009227300008250404000", "0163103100000000"); // J. Luehder AOD Compare

  //----------------------------- configuration for run 2 analysis 5 TeV ----------------------------------------------------
  } else if (trainConfig == 1050){
    cuts.AddCutPCM("00010113", "00200009227300008250404000", "0152101500000000"); //old standard cut pp 5 TeV VAND
  } else if (trainConfig == 1051){
    cuts.AddCutPCM("00010113", "00200009227300008250404000", "0152103500000000"); //new standard cut pp 5 TeV VAND
  } else if (trainConfig == 1052){
    cuts.AddCutPCM("00010113", "00a00009227300008250404000", "0152103500000000"); //new standard cut pp 5 TeV VAND
  } else if (trainConfig == 1053){
    cuts.AddCutPCM("00010113", "00b00009227300008250404000", "0152103500000000"); //new standard cut pp 5 TeV VAND
  } else if (trainConfig == 1054){
    cuts.AddCutPCM("00010113", "00c00009227300008250404000", "0152103500000000"); //new standard cut pp 5 TeV VAND
  } else if (trainConfig == 1055){
    cuts.AddCutPCM("00010113", "00200009227300008250704000", "0152103500000000"); //test cosPA scan
  } else if (trainConfig == 1056){
    cuts.AddCutPCM("00010113", "00200009227300008250804000", "0152103500000000"); //test cosPA scan
  } else if (trainConfig == 1057){
    cuts.AddCutPCM("00010113", "00200009227300008250904000", "0152103500000000"); //test cosPA scan
  } else if (trainConfig == 1058){
    cuts.AddCutPCM("00010113", "00200009227300008250a04000", "0152103500000000"); //test cosPA scan
  } else if (trainConfig == 1059){
    cuts.AddCutPCM("00010113", "00200009a27300008250904120", "0152103500000000"); //cosPA, 0.99 eta 0.9
  } else if (trainConfig == 1060){
    cuts.AddCutPCM("00010113", "0d200009a27300008250904120", "0152103500000000"); //cosPA, 0.99 eta 0.8
  } else if (trainConfig == 1061){
    cuts.AddCutPCM("00010113", "00200009a27300008250a04120", "0152103500000000"); //cosPA, 0.995 eta 0.9
  } else if (trainConfig == 1062){
    cuts.AddCutPCM("00010113", "0d200009a27300008250a04120", "0152103500000000"); //cosPA, 0.995 eta 0.8
  //----------------------------- configuration for run 2 analysis 13 TeV Triggers --------------------------------------------
  } else if (trainConfig == 1070) { // EMC triggers -50, +30 ns
    cuts.AddCutPCM("00010113", "00200009227300008250404000", "0163103100000000","1111100060032220000"); //INT7
  } else if (trainConfig == 1071) { // EMC triggers -50, +30 ns
    cuts.AddCutPCM("00085113", "00200009227300008250404000", "0163103100000000","1111100060032220000"); //EG2
  } else if (trainConfig == 1072) { // EMC triggers -50, +30 ns
    cuts.AddCutPCM("00083113", "00200009227300008250404000", "0163103100000000","1111100060032220000"); //EG1
  } else if (trainConfig == 1073) { // DCAL triggers -50, +30 ns
    cuts.AddCutPCM("00089113", "00200009227300008250404000", "0163103100000000","3885500060032220000"); //DG2
  } else if (trainConfig == 1074) { // DCAL triggers -50, +30 ns
    cuts.AddCutPCM("0008b113", "00200009227300008250404000", "0163103100000000","3885500060032220000"); //DG1

   //----------------------------- configuration for run 2 analysis 13 TeVLowB --------------------------------------------
  } else if (trainConfig == 1080) { //
    cuts.AddCutPCM("00010113", "00200089227300008280404000", "0152103500000000"); //standard cut

  } else if (trainConfig == 1090) { //Standard cut for pp 5 TeV analysis VAND
    cuts.AddCutPCM("00010113", "0d200009227300008250404000", "0152103500000000"); //

  } else if (trainConfig == 1100) { //new cut for pp 8 TeV analysis VAND triggers and MB for RpA calcluation
    cuts.AddCutPCM("00010113", "00200009f9730000dge0400000", "0152103500000000","1111101060032230000"); // new standard dEdx, pt dep Qt, chi2-psipair exp
  } else if (trainConfig == 1101) { //new cut for pp 8 TeV analysis VAND triggers and MB for RpA calcluation
    cuts.AddCutPCM("00052113", "00200009f9730000dge0400000", "0152103500000000","1111101060032230000"); // new standard dEdx, pt dep Qt, chi2-psipair exp
    cuts.AddCutPCM("00081113", "00200009f9730000dge0400000", "0152103500000000","1111101060032230000"); // new standard dEdx, pt dep Qt, chi2-psipair exp


  } else if (trainConfig == 1500) { // TOF single leg cut
    cuts.AddCutPCM("00010113", "0dm00009f9730600dge0404000", "0152103500000000"); // new default (R region rej. + eta<0.8 + DC)
  } else if (trainConfig == 1501) { // TOF both leg cut
    cuts.AddCutPCM("00010113", "0dm00009f9730700dge0404000", "0152103500000000"); // new default (R region rej. + eta<0.8 + DC)
  } else if (trainConfig == 1502) { // TOF single leg cut
    cuts.AddCutPCM("00010113", "0dm00009f9730800dge0404000", "0152103500000000"); // new default (R region rej. + eta<0.8 + DC)
  } else if (trainConfig == 1503) { // TOF both leg cut
    cuts.AddCutPCM("00010113", "0dm00009f9730900dge0404000", "0152103500000000"); // new default (R region rej. + eta<0.8 + DC)


  } else if (trainConfig == 2001) { // Double Gap event selection  special trigger 11
    cuts.AddCutPCM("000b0113", "0d200009267300008250404000", "0152103500000000"); // eta < 0.8

 // special V0AND configurations for cut QA
  } else if (trainConfig == 2400){
    // std: |eta|<0.8, 5<R<180cm, TPCcls>0.6, -3<nsigE<5, nsigPi<1 (0.4-3.5GeV), nsigPi<-10 (3.5GeV+)
    // qT<0.05, chi2<30, PsiPair<0.1 (2D), photon asym < 0.95, cosPA>0.85, toocloseV0+doublecounting
    cuts.AddCutPCM("00010113", "0d200009227300008250404000", "0152103500000000"); // Standard
  } else if (trainConfig == 2401){ // chi2, PsiPair
    cuts.AddCutPCM("00010113", "0d200009227300008c50404000", "0152103500000000"); // chi2<40
    cuts.AddCutPCM("00010113", "0d200009227300008150404000", "0152103500000000"); // chi2<50
    cuts.AddCutPCM("00010113", "0d2000092273000082c0404000", "0152103500000000"); // PsiPair<0.15
    cuts.AddCutPCM("00010113", "0d200009227300008280404000", "0152103500000000"); // PsiPair<0.20
    cuts.AddCutPCM("00010113", "0d200009227300008cc0404000", "0152103500000000"); // PsiPair<0.15 + chi2<40
    cuts.AddCutPCM("00010113", "0d200009227300008180404000", "0152103500000000"); // PsiPair<0.20 + chi2<50
  } else if (trainConfig == 2402){ // qT 1D pT dep
    cuts.AddCutPCM("00010113", "0d200009227300001250404000", "0152103500000000"); // qT<0.1 (1D)
    cuts.AddCutPCM("00010113", "0d20000922730000d250404000", "0152103500000000"); // qT<0.110pT (1D)
    cuts.AddCutPCM("00010113", "0d20000922730000b250404000", "0152103500000000"); // qT<0.125pT (1D)
    cuts.AddCutPCM("00010113", "0d20000922730000f250404000", "0152103500000000"); // qT<0.130pT (1D)
  } else if (trainConfig == 2403){ // qT 2D pT dep
    cuts.AddCutPCM("00010113", "0d20000922730000c250404000", "0152103500000000"); // qT<0.110pT (2D) alpha<0.95
    cuts.AddCutPCM("00010113", "0d20000922730000a250404000", "0152103500000000"); // qT<0.125pT (2D) alpha<0.95
    cuts.AddCutPCM("00010113", "0d20000922730000e250404000", "0152103500000000"); // qT<0.130pT (2D) alpha<0.95
  } else if (trainConfig == 2404){ // qT 2D pT dep
    cuts.AddCutPCM("00010113", "0d20000922730000c259404000", "0152103500000000"); // qT<0.110pT (2D) alpha<0.99
    cuts.AddCutPCM("00010113", "0d20000922730000a259404000", "0152103500000000"); // qT<0.125pT (2D) alpha<0.99
    cuts.AddCutPCM("00010113", "0d20000922730000e259404000", "0152103500000000"); // qT<0.130pT (2D) alpha<0.99
  } else if (trainConfig == 2405){ // qT 2D pT dep
    cuts.AddCutPCM("00010113", "0d20000922730000c25a404000", "0152103500000000"); // qT<0.110pT (2D) alpha<1
    cuts.AddCutPCM("00010113", "0d20000922730000a25a404000", "0152103500000000"); // qT<0.125pT (2D) alpha<1
    cuts.AddCutPCM("00010113", "0d20000922730000e25a404000", "0152103500000000"); // qT<0.130pT (2D) alpha<1
  } else if (trainConfig == 2406){ // combined configs
    cuts.AddCutPCM("00010113", "0d20000922730000acc0404000", "0152103500000000"); // PsiPair<0.15 + chi2<40 + qT<0.125pT (2D) alpha<0.95
    cuts.AddCutPCM("00010113", "0d20000922730000a180404000", "0152103500000000"); // PsiPair<0.20 + chi2<50 + qT<0.125pT (2D) alpha<0.95
    cuts.AddCutPCM("00010113", "0d20000922730000acca404000", "0152103500000000"); // PsiPair<0.15 + chi2<40 + qT<0.125pT (2D) alpha<1
    cuts.AddCutPCM("00010113", "0d20000922730000a18a404000", "0152103500000000"); // PsiPair<0.20 + chi2<50 + qT<0.125pT (2D) alpha<1
  } else if (trainConfig == 2407){ // chi2, PsiPair cont.
    cuts.AddCutPCM("00010113", "0d200009227300008fd0404000", "0152103500000000"); // PsiPair<0.15exp(-0.065chi2)
    cuts.AddCutPCM("00010113", "0d200009227300008ge0404000", "0152103500000000"); // PsiPair<0.18exp(-0.055chi2)
    cuts.AddCutPCM("00010113", "0d200009227300008hf0404000", "0152103500000000"); // PsiPair<0.20exp(-0.050chi2)
  } else if (trainConfig == 2408){ // combined configs cont.
    cuts.AddCutPCM("00010113", "0d20000922730000afd0404000", "0152103500000000"); // PsiPair<0.15exp(-0.065chi2) + qT<0.125pT (2D) alpha<0.95
    cuts.AddCutPCM("00010113", "0d20000922730000age0404000", "0152103500000000"); // PsiPair<0.18exp(-0.055chi2) + qT<0.125pT (2D) alpha<0.95
    cuts.AddCutPCM("00010113", "0d20000922730000ahf0404000", "0152103500000000"); // PsiPair<0.20exp(-0.050chi2) + qT<0.125pT (2D) alpha<0.95
  } else if (trainConfig == 2409){ // combined configs cont.
    cuts.AddCutPCM("00010113", "0d20000922730000gfda404000", "0152103500000000"); // PsiPair<0.15exp(-0.065chi2) + qT<0.140pT (2D) alpha<1
    cuts.AddCutPCM("00010113", "0d20000922730000ggea404000", "0152103500000000"); // PsiPair<0.18exp(-0.055chi2) + qT<0.140pT (2D) alpha<1
    cuts.AddCutPCM("00010113", "0d20000922730000ghfa404000", "0152103500000000"); // PsiPair<0.20exp(-0.050chi2) + qT<0.140pT (2D) alpha<1
  // EGA
  } else if (trainConfig == 2410){
    // std: |eta|<0.8, 5<R<180cm, TPCcls>0.6, -3<nsigE<5, nsigPi<1 (0.4-3.5GeV), nsigPi<-10 (3.5GeV+)
    // qT<0.05, chi2<30, PsiPair<0.1 (2D), photon asym < 0.95, cosPA>0.85, toocloseV0+doublecounting
    cuts.AddCutPCM("00081113", "0d200009227300008250404000", "0152103500000000"); // Standard
  } else if (trainConfig == 2411){ // chi2, PsiPair
    cuts.AddCutPCM("00081113", "0d200009227300008c50404000", "0152103500000000"); // chi2<40
    cuts.AddCutPCM("00081113", "0d200009227300008150404000", "0152103500000000"); // chi2<50
    cuts.AddCutPCM("00081113", "0d2000092273000082c0404000", "0152103500000000"); // PsiPair<0.15
    cuts.AddCutPCM("00081113", "0d200009227300008280404000", "0152103500000000"); // PsiPair<0.20
    cuts.AddCutPCM("00081113", "0d200009227300008cc0404000", "0152103500000000"); // PsiPair<0.15 + chi2<40
    cuts.AddCutPCM("00081113", "0d200009227300008180404000", "0152103500000000"); // PsiPair<0.20 + chi2<50
  } else if (trainConfig == 2412){ // qT 1D pT dep
    cuts.AddCutPCM("00081113", "0d200009227300001250404000", "0152103500000000"); // qT<0.1 (1D)
    cuts.AddCutPCM("00081113", "0d20000922730000d250404000", "0152103500000000"); // qT<0.110pT (1D)
    cuts.AddCutPCM("00081113", "0d20000922730000b250404000", "0152103500000000"); // qT<0.125pT (1D)
    cuts.AddCutPCM("00081113", "0d20000922730000f250404000", "0152103500000000"); // qT<0.130pT (1D)
  } else if (trainConfig == 2413){ // qT 2D pT dep
    cuts.AddCutPCM("00081113", "0d20000922730000c250404000", "0152103500000000"); // qT<0.110pT (2D) alpha<0.95
    cuts.AddCutPCM("00081113", "0d20000922730000a250404000", "0152103500000000"); // qT<0.125pT (2D) alpha<0.95
    cuts.AddCutPCM("00081113", "0d20000922730000e250404000", "0152103500000000"); // qT<0.130pT (2D) alpha<0.95
  } else if (trainConfig == 2414){ // qT 2D pT dep
    cuts.AddCutPCM("00081113", "0d20000922730000c259404000", "0152103500000000"); // qT<0.110pT (2D) alpha<0.99
    cuts.AddCutPCM("00081113", "0d20000922730000a259404000", "0152103500000000"); // qT<0.125pT (2D) alpha<0.99
    cuts.AddCutPCM("00081113", "0d20000922730000e259404000", "0152103500000000"); // qT<0.130pT (2D) alpha<0.99
  } else if (trainConfig == 2415){ // qT 2D pT dep
    cuts.AddCutPCM("00081113", "0d20000922730000c25a404000", "0152103500000000"); // qT<0.110pT (2D) alpha<1
    cuts.AddCutPCM("00081113", "0d20000922730000a25a404000", "0152103500000000"); // qT<0.125pT (2D) alpha<1
    cuts.AddCutPCM("00081113", "0d20000922730000e25a404000", "0152103500000000"); // qT<0.130pT (2D) alpha<1
  } else if (trainConfig == 2416){ // combined configs
    cuts.AddCutPCM("00081113", "0d20000922730000acc0404000", "0152103500000000"); // PsiPair<0.15 + chi2<40 + qT<0.125pT (2D) alpha<0.95
    cuts.AddCutPCM("00081113", "0d20000922730000a180404000", "0152103500000000"); // PsiPair<0.20 + chi2<50 + qT<0.125pT (2D) alpha<0.95
    cuts.AddCutPCM("00081113", "0d20000922730000acca404000", "0152103500000000"); // PsiPair<0.15 + chi2<40 + qT<0.125pT (2D) alpha<1
    cuts.AddCutPCM("00081113", "0d20000922730000a18a404000", "0152103500000000"); // PsiPair<0.20 + chi2<50 + qT<0.125pT (2D) alpha<1
  } else if (trainConfig == 2417){ // chi2, PsiPair cont.
    cuts.AddCutPCM("00081113", "0d200009227300008fd0404000", "0152103500000000"); // PsiPair<0.15exp(-0.065chi2)
    cuts.AddCutPCM("00081113", "0d200009227300008ge0404000", "0152103500000000"); // PsiPair<0.18exp(-0.055chi2)
    cuts.AddCutPCM("00081113", "0d200009227300008hf0404000", "0152103500000000"); // PsiPair<0.20exp(-0.050chi2)
  } else if (trainConfig == 2418){ // combined configs cont.
    cuts.AddCutPCM("00081113", "0d20000922730000afd0404000", "0152103500000000"); // PsiPair<0.15exp(-0.065chi2) + qT<0.125pT (2D) alpha<0.95
    cuts.AddCutPCM("00081113", "0d20000922730000age0404000", "0152103500000000"); // PsiPair<0.18exp(-0.055chi2) + qT<0.125pT (2D) alpha<0.95
    cuts.AddCutPCM("00081113", "0d20000922730000ahf0404000", "0152103500000000"); // PsiPair<0.20exp(-0.050chi2) + qT<0.125pT (2D) alpha<0.95
  } else if (trainConfig == 2419){ // combined configs cont.
    cuts.AddCutPCM("00081113", "0d20000922730000gfda404000", "0152103500000000"); // PsiPair<0.15exp(-0.065chi2) + qT<0.140pT (2D) alpha<1
    cuts.AddCutPCM("00081113", "0d20000922730000ggea404000", "0152103500000000"); // PsiPair<0.18exp(-0.055chi2) + qT<0.140pT (2D) alpha<1
    cuts.AddCutPCM("00081113", "0d20000922730000ghfa404000", "0152103500000000"); // PsiPair<0.20exp(-0.050chi2) + qT<0.140pT (2D) alpha<1
  // EMC7
  } else if (trainConfig == 2426){ // combined configs
    cuts.AddCutPCM("00052113", "0d20000922730000acc0404000", "0152103500000000"); // PsiPair<0.15 + chi2<40 + qT<0.125pT (2D) alpha<0.95
    cuts.AddCutPCM("00052113", "0d20000922730000a180404000", "0152103500000000"); // PsiPair<0.20 + chi2<50 + qT<0.125pT (2D) alpha<0.95
    cuts.AddCutPCM("00052113", "0d20000922730000acca404000", "0152103500000000"); // PsiPair<0.15 + chi2<40 + qT<0.125pT (2D) alpha<1
    cuts.AddCutPCM("00052113", "0d20000922730000a18a404000", "0152103500000000"); // PsiPair<0.20 + chi2<50 + qT<0.125pT (2D) alpha<1
  // EG1+DG1
  } else if (trainConfig == 2436){ // combined configs
    cuts.AddCutPCM("0008d113", "0d20000922730000acc0404000", "0152103500000000"); // PsiPair<0.15 + chi2<40 + qT<0.125pT (2D) alpha<0.95
    cuts.AddCutPCM("0008d113", "0d20000922730000a180404000", "0152103500000000"); // PsiPair<0.20 + chi2<50 + qT<0.125pT (2D) alpha<0.95
    cuts.AddCutPCM("0008d113", "0d20000922730000acca404000", "0152103500000000"); // PsiPair<0.15 + chi2<40 + qT<0.125pT (2D) alpha<1
    cuts.AddCutPCM("0008d113", "0d20000922730000a18a404000", "0152103500000000"); // PsiPair<0.20 + chi2<50 + qT<0.125pT (2D) alpha<1

  //*************************************************************************************************
  // 7 TeV PCM sys for omega analysis variations
  //*************************************************************************************************
  } else if ( trainConfig == 2600){ // std
    cuts.AddCutPCM("00000113","00200009227000008250400000","0163103100000010");
  } else if ( trainConfig == 2601){ // pileup
    cuts.AddCutPCM("00000113","00200009227000008250400000","0163103100000010");
  } else if ( trainConfig == 2602){ // singlept
    cuts.AddCutPCM("00000113","00200019227000008250400000","0163103100000010"); // 0.100 GeV
    cuts.AddCutPCM("00000113","00200049227000008250400000","0163103100000010"); // 0.075 GeV
    cuts.AddCutPCM("00000113","00200069227000008250400000","0163103100000010"); // 0.04 GeV
    cuts.AddCutPCM("00000113","00200059227000008250400000","0163103100000010"); // 0.125 GeV
  } else if ( trainConfig == 2603){ // clstpc
    cuts.AddCutPCM("00000113","00200008227000008250400000","0163103100000010"); // fMinClsTPCToF= 0.35;fUseCorrectedTPCClsInfo=1;
    cuts.AddCutPCM("00000113","00200006227000008250400000","0163103100000010"); // fMinClsTPCToF= 0.70;fUseCorrectedTPCClsInfo=1;
    cuts.AddCutPCM("00000113","00200001227000008250400000","0163103100000010"); // fMinClsTPCToF= 60;fUseCorrectedTPCClsInfo=0;
    cuts.AddCutPCM("00000113","00200002227000008250400000","0163103100000010"); // fMinClsTPCToF= 80;fUseCorrectedTPCClsInfo=0;
    cuts.AddCutPCM("00000113","00200003227000008250400000","0163103100000010"); // fMinClsTPCToF= 100;fUseCorrectedTPCClsInfo=0;
  } else if ( trainConfig == 2604){ // TPCdEdxCutElectron
    cuts.AddCutPCM("00000113","00200009327000008250400000","0163103100000010"); // -4,5
    cuts.AddCutPCM("00000113","00200009627000008250400000","0163103100000010"); // -2.5,4
    cuts.AddCutPCM("00000113","00200009427000008250400000","0163103100000010"); // -6,7
    cuts.AddCutPCM("00000113","00200009527000008250400000","0163103100000010"); // -4,4
    cuts.AddCutPCM("00000113","00200009627000008250400000","0163103100000010"); // -2.5,4
  } else if ( trainConfig == 2605){ // TPCdEdxCutPion
    cuts.AddCutPCM("00000113","00200009217000008250400000","0163103100000010"); // fPIDnSigmaAbovePionLine=0; fPIDnSigmaAbovePionLineHighPt=-10;
    cuts.AddCutPCM("00000113","00200009237000008250400000","0163103100000010"); // fPIDnSigmaAbovePionLine=2.5; fPIDnSigmaAbovePionLineHighPt=-10;
    cuts.AddCutPCM("00000113","00200009247000008250400000","0163103100000010"); // fPIDnSigmaAbovePionLine=0.5; fPIDnSigmaAbovePionLineHighPt=-10;
    cuts.AddCutPCM("00000113","00200009257000008250400000","0163103100000010"); // fPIDnSigmaAbovePionLine=2; fPIDnSigmaAbovePionLineHighPt=-10;
  } else if ( trainConfig == 2606){ // QtMaxCut
    cuts.AddCutPCM("00000113","00200009227000003250400000","0163103100000010");
    cuts.AddCutPCM("00000113","00200009227000009250400000","0163103100000010");
    cuts.AddCutPCM("00000113","00200009227000002250400000","0163103100000010");
    cuts.AddCutPCM("00000113","00200009227000006250400000","0163103100000010");
  } else if ( trainConfig == 2607){ // Chi2GammaCut
    cuts.AddCutPCM("00000113","00200009227000008150400000","0163103100000010");
    cuts.AddCutPCM("00000113","00200009227000008850400000","0163103100000010");
    cuts.AddCutPCM("00000113","00200009227000008a50400000","0163103100000010");
    cuts.AddCutPCM("00000113","00200009227000008950400000","0163103100000010");
  } else if ( trainConfig == 2608){ // PsiPair
    cuts.AddCutPCM("00000113","00200009227000008260400000","0163103100000010");
    cuts.AddCutPCM("00000113","00200009227000008280400000","0163103100000010");
    cuts.AddCutPCM("00000113","00200009227000008210400000","0163103100000010");
    cuts.AddCutPCM("00000113","00200009227000008220400000","0163103100000010");

  } else {
    Error(Form("GammaConvV1_%i",trainConfig), "wrong trainConfig variable no cuts have been specified for the configuration");
    return;
  }

if(!cuts.AreValid()){
    cout << "\n\n****************************************************" << endl;
    cout << "ERROR: No valid cuts stored in CutHandlerPCM! Returning..." << endl;
    cout << "****************************************************\n\n" << endl;
    return;
  }

  Int_t numberOfCuts = cuts.GetNCuts();

  TList *EventCutList = new TList();
  TList *ConvCutList = new TList();
  TList *MesonCutList = new TList();
  TList *ClusterCutList = new TList();

  TList *HeaderList = new TList();
  if (periodNameAnchor.Contains("LHC12i3")){
    TObjString *Header2 = new TObjString("BOX");
    HeaderList->Add(Header2);
  } else if (periodNameAnchor.CompareTo("LHC14e2b")==0){
    TObjString *Header2 = new TObjString("pi0_1");
    HeaderList->Add(Header2);
    TObjString *Header3 = new TObjString("eta_2");
    HeaderList->Add(Header3);
  }

  TString energy = "";
  TString mcName = "";
  TString mcNameAdd = "";
  if (periodNameAnchor.Contains("WOSDD")){
    mcNameAdd = "_WOSDD";
  } else if (periodNameAnchor.Contains("WSDD")){
    mcNameAdd = "_WSDD";
  }
  if (periodNameAnchor.Contains("LHC12i3")){
    energy = "2760GeV";
    mcName = "Pythia8_LHC12i3";
  } else if (periodNameAnchor.Contains("LHC12f1a")){
    energy = "2760GeV";
    mcName = "Pythia8_LHC12f1a";
  } else if (periodNameAnchor.Contains("LHC12f1b")){
    energy = "2760GeV";
    mcName = "Phojet_LHC12f1b";
  } else if (periodNameAnchor.Contains("LHC14e2a")){
    energy = "8TeV";
    mcName = "Pythia8_LHC14e2a";
  } else if (periodNameAnchor.Contains("LHC14e2b")){
    energy = "8TeV";
    mcName = "Pythia8_LHC14e2b";
  } else if (periodNameAnchor.Contains("LHC14e2c")){
    energy = "8TeV";
    mcName = "Phojet_LHC14e2c";
  }

  EventCutList->SetOwner(kTRUE);
  AliConvEventCuts **analysisEventCuts        = new AliConvEventCuts*[numberOfCuts];
  ConvCutList->SetOwner(kTRUE);
  AliConversionPhotonCuts **analysisCuts      = new AliConversionPhotonCuts*[numberOfCuts];
  MesonCutList->SetOwner(kTRUE);
  AliConversionMesonCuts **analysisMesonCuts  = new AliConversionMesonCuts*[numberOfCuts];
  ClusterCutList->SetOwner(kTRUE);
  AliCaloPhotonCuts **analysisClusterCuts     = new AliCaloPhotonCuts*[numberOfCuts];
  Bool_t enableClustersForTrigger             = kFALSE;
  Bool_t initializedMatBudWeigths_existing    = kFALSE;

  for(Int_t i = 0; i<numberOfCuts; i++){
    analysisEventCuts[i]          = new AliConvEventCuts();
    TString fitNamePi0            = Form("Pi0_Fit_Data_%s",energy.Data());
    TString fitNameEta            = Form("Eta_Fit_Data_%s",energy.Data());
    TString fAddedSignalString    = (cuts.GetEventCut(i)).Data();
    fAddedSignalString            = fAddedSignalString(6,1);
    Bool_t fAddedSignal           = kFALSE;
    if (fAddedSignalString.CompareTo("2") == 0)
      fAddedSignal                = kTRUE;

    TString mcInputNamePi0        = "";
    TString mcInputNameEta        = "";
    if (fAddedSignal && (periodNameAnchor.Contains("LHC12i3") || periodNameAnchor.CompareTo("LHC14e2b")==0)){
      mcInputNamePi0              = Form("Pi0_%s%s_addSig_%s", mcName.Data(), mcNameAdd.Data(), energy.Data() );
      mcInputNameEta              = Form("Eta_%s%s_addSig_%s", mcName.Data(), mcNameAdd.Data(), energy.Data() );
    } else {
      mcInputNamePi0              = Form("Pi0_%s%s_%s", mcName.Data(), mcNameAdd.Data(), energy.Data() );
      mcInputNameEta              = Form("Eta_%s%s_%s", mcName.Data(), mcNameAdd.Data(), energy.Data() );
    }
    //if (doWeightingPart) analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kFALSE, kFALSE, kFALSE, fileNamePtWeights, mcInputNamePi0, mcInputNameEta, "",fitNamePi0,fitNameEta);
    if (doWeightingPart) analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE, fileNamePtWeights, mcInputNamePi0, mcInputNameEta, "",fitNamePi0,fitNameEta);

    TString dataInputMultHisto    = "";
    TString mcInputMultHisto      = "";
    TString triggerString         = (cuts.GetEventCut(i)).Data();
    triggerString                 = triggerString(3,2);

    if (triggerString.CompareTo("03")==0)
      triggerString               = "00";
    if (periodNameAnchor.CompareTo("LHC13g") == 0 && triggerString.CompareTo("10")== 0 )
      triggerString               = "00";
    if ( ((periodNameAnchor.CompareTo("LHC16d") == 0) || (periodNameAnchor.CompareTo("LHC17c") == 0) || (periodNameAnchor.CompareTo("LHC18b") == 0)) && triggerString.CompareTo("10")== 0 )
      triggerString               = "000";


    dataInputMultHisto            = Form("%s_%s", periodNameAnchor.Data(), triggerString.Data());
    mcInputMultHisto              = Form("%s_%s", periodNameV0Reader.Data(), triggerString.Data());

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
    if(fJetFinderUsage)
      analysisEventCuts[i]->SetUseJetFinderForOutliers(kTRUE);
    if(fUsePtHardFromFile)
      analysisEventCuts[i]->SetUsePtHardBinFromFile(kTRUE);
    if(fUseAddOutlierRej)
      analysisEventCuts[i]->SetUseAdditionalOutlierRejection(kTRUE);
    analysisEventCuts[i]->SetV0ReaderName(V0ReaderName);
    analysisEventCuts[i]->SetCorrectionTaskSetting(corrTaskSetting);
    if (periodNameV0Reader.CompareTo("") != 0) analysisEventCuts[i]->SetPeriodEnum(periodNameV0Reader);
    analysisEventCuts[i]->SetLightOutput(enableLightOutput);
    analysisEventCuts[i]->InitializeCutsFromCutString((cuts.GetEventCut(i)).Data());

    EventCutList->Add(analysisEventCuts[i]);
    analysisEventCuts[i]->SetFillCutHistograms("",kFALSE);
    if(!cuts.GetClusterCut(i).CompareTo("")){
      cout << "\nNo cluster cut set, not filling cluster histograms for triggers ...\n" << endl;
    } else {
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

        enableClustersForTrigger  = kTRUE;
        analysisClusterCuts[i]    = new AliCaloPhotonCuts();
        analysisClusterCuts[i]->SetV0ReaderName(V0ReaderName);
        analysisClusterCuts[i]->SetCorrectionTaskSetting(corrTaskSetting);
        analysisClusterCuts[i]->SetCaloTrackMatcherName(TrackMatcherName);
        analysisClusterCuts[i]->SetLightOutput(enableLightOutput);
        analysisClusterCuts[i]->InitializeCutsFromCutString((cuts.GetClusterCut(i)).Data());
        ClusterCutList->Add(analysisClusterCuts[i]);
        analysisClusterCuts[i]->SetFillCutHistograms("");
    }

    analysisCuts[i]               = new AliConversionPhotonCuts();
    if ( trainConfig == 89 ){
      analysisCuts[i]->SetSwitchToKappaInsteadOfNSigdEdxTPC(kTRUE);
    }
    if (enableMatBudWeightsPi0 > 0){
        if (isMC > 0){
            if (analysisCuts[i]->InitializeMaterialBudgetWeights(enableMatBudWeightsPi0,fileNameMatBudWeights)){
                initializedMatBudWeigths_existing = kTRUE;
		cout << "MBW properly initialized" << endl;
	    }
            else {cout << "ERROR The initialization of the materialBudgetWeights did not work out." << endl;}
        }
        else {cout << "ERROR 'enableMatBudWeightsPi0'-flag was set > 0 even though this is not a MC task. It was automatically reset to 0." << endl;}
    }
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

    analysisCuts[i]->SetV0ReaderName(V0ReaderName);
    analysisCuts[i]->SetLightOutput(enableLightOutput);
    analysisCuts[i]->InitializeCutsFromCutString((cuts.GetPhotonCut(i)).Data());
    if (trainConfig == 1021){
      analysisCuts[i]->SetDodEdxSigmaCut(kFALSE);
    }

    ConvCutList->Add(analysisCuts[i]);
    analysisCuts[i]->SetFillCutHistograms("",kFALSE);

    analysisMesonCuts[i]          = new AliConversionMesonCuts();
    analysisMesonCuts[i]->SetLightOutput(enableLightOutput);
    analysisMesonCuts[i]->InitializeCutsFromCutString((cuts.GetMesonCut(i)).Data());
    MesonCutList->Add(analysisMesonCuts[i]);
    analysisMesonCuts[i]->SetFillCutHistograms("");
    analysisEventCuts[i]->SetAcceptedHeader(HeaderList);
    if(doSmear) analysisMesonCuts[i]->SetDefaultSmearing(bremSmear,smearPar,smearParConst);

  }

  task->SetEventCutList(numberOfCuts,EventCutList);
  task->SetConversionCutList(numberOfCuts,ConvCutList);
  task->SetMesonCutList(numberOfCuts,MesonCutList);
  task->SetMoveParticleAccordingToVertex(kTRUE);
  task->SetDoMesonAnalysis(kTRUE);
  task->SetDoMesonQA(enableQAMesonTask); //Attention new switch for Pi0 QA
  task->SetDoPhotonQA(enableQAPhotonTask);  //Attention new switch small for Photon QA
  task->SetDoChargedPrimary(enableChargedPrimary);
  if (enableClustersForTrigger){
    task->SetDoClusterSelectionForTriggerNorm(enableClustersForTrigger);
    task->SetClusterCutList(numberOfCuts,ClusterCutList);
  }
  if (initializedMatBudWeigths_existing) {
      task->SetDoMaterialBudgetWeightingOfGammasForTrueMesons(kTRUE);
  }

  //connect containers
  AliAnalysisDataContainer *coutput =
    mgr->CreateContainer(!(corrTaskSetting.CompareTo("")) ? Form("GammaConvV1_%i",trainConfig) : Form("GammaConvV1_%i_%s",trainConfig,corrTaskSetting.Data()), TList::Class(),
              AliAnalysisManager::kOutputContainer, Form("GammaConvV1_%i.root",trainConfig) );

  mgr->AddTask(task);
  mgr->ConnectInput(task,0,cinput);
  mgr->ConnectOutput(task,1,coutput);
  Int_t nContainer = 2;
  for(Int_t i = 0; i<numberOfCuts; i++){
    if(enableQAPhotonTask>1){
      if (initializedMatBudWeigths_existing) {
	mgr->ConnectOutput(task,nContainer,mgr->CreateContainer(Form("%s_%s_%s MBW Photon DCA tree",(cuts.GetEventCut(i)).Data(),(cuts.GetPhotonCut(i)).Data(),(cuts.GetMesonCut(i)).Data()), TTree::Class(), AliAnalysisManager::kOutputContainer, Form("GammaConvV1_%i.root",trainConfig)) );
      }else{
	mgr->ConnectOutput(task,nContainer,mgr->CreateContainer(Form("%s_%s_%s Photon DCA tree",(cuts.GetEventCut(i)).Data(),(cuts.GetPhotonCut(i)).Data(),(cuts.GetMesonCut(i)).Data()), TTree::Class(), AliAnalysisManager::kOutputContainer, Form("GammaConvV1_%i.root",trainConfig)) );
      }
      nContainer++;
    }
    if(enableQAMesonTask>1){
      if (initializedMatBudWeigths_existing) {
	mgr->ConnectOutput(task,nContainer,mgr->CreateContainer(Form("%s_%s_%s MBW Meson DCA tree",(cuts.GetEventCut(i)).Data(),(cuts.GetPhotonCut(i)).Data(),(cuts.GetMesonCut(i)).Data()), TTree::Class(), AliAnalysisManager::kOutputContainer, Form("GammaConvV1_%i.root",trainConfig)) );
      }else{
	mgr->ConnectOutput(task,nContainer,mgr->CreateContainer(Form("%s_%s_%s Meson DCA tree",(cuts.GetEventCut(i)).Data(),(cuts.GetPhotonCut(i)).Data(),(cuts.GetMesonCut(i)).Data()), TTree::Class(), AliAnalysisManager::kOutputContainer, Form("GammaConvV1_%i.root",trainConfig)) );
      }
      nContainer++;
    }
  }

  return;

}
