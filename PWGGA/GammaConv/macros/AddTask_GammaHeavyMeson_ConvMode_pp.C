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
void AddTask_GammaHeavyMeson_ConvMode_pp(
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
  // meson reco mode: 0 - PCM-PCM, 1 - PCM-Calo, 2 - Calo-Calo
  Int_t mesonRecoMode = 0;

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
  if(rmaxFacPtHardSetting->GetEntries()<1){cout << "ERROR: AddTask_GammaHeavyMeson_ConvMode_pp during parsing of settingMaxFacPtHard String '" << settingMaxFacPtHard.Data() << "'" << endl; return;}
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

  // ================== GetAnalysisManager ===============================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error(Form("%s_%i", addTaskName.Data(),  trainConfig), "No analysis manager found.");
    return ;
  }

  // ================== GetInputEventHandler =============================
  AliVEventHandler *inputHandler=mgr->GetInputEventHandler();

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
  //================================================
  AliAnalysisTaskHeavyNeutralMesonToGG *task=NULL;
  task= new AliAnalysisTaskHeavyNeutralMesonToGG(Form("HeavyNeutralMesonToGG_%i_%i_%i", mesonRecoMode, selectedMeson, trainConfig));
  task->SetIsHeavyIon(isHeavyIon);
  task->SetIsMC(isMC);
  task->SetV0ReaderName(V0ReaderName);
  task->SetCorrectionTaskSetting(corrTaskSetting);

  task->SetMesonRecoMode(mesonRecoMode); // meson reco mode: 0 - PCM-PCM, 1 - PCM-Calo, 2 - Calo-Calo
  if (enableLightOutput > 1) task->SetLightOutput(kTRUE);

  // *********************************************************************************************************
  // 2.76 TeV  pp Run1
  // *********************************************************************************************************
  if(trainConfig == 1){    // various standard cuts
    cuts.AddCutPCM("00000113", "00200009397300008250400000", "0163103000900000"); // new pi0/eta cut 2.76TeV
    cuts.AddCutPCM("00000113", "00200009397300008250400000", "0163103000000000"); // new pi0/eta cut 2.76TeV without MC smearing
    cuts.AddCutPCM("00000113", "00200009366300003800000000", "0163103000900000"); // standard cut Pi0 pp 2.76TeV PbPb paper 2012
  } else if (trainConfig == 2) { // various standard cuts added signals
    cuts.AddCutPCM("00000123", "00200009397300008250400000", "0163103000900000"); // new pi0/eta cut 2.76TeV
    cuts.AddCutPCM("00000123", "00200009397300008250400000", "0163103000000000"); // new pi0/eta cut 2.76TeV without MC smearing
    cuts.AddCutPCM("00000123", "00200009366300003800000000", "0163103000900000"); // standard cut Pi0 pp 2.76TeV PbPb paper 2012
  } else if (trainConfig == 3) { // additional standards
    cuts.AddCutPCM("00000113", "00200009297002008250400000", "0163103000900000"); // standard cut LHC11h pp 2.76TeV
    cuts.AddCutPCM("00000113", "00200009227302008250404000", "0163101000000000"); // Ana eta analysis prefered
    cuts.AddCutPCM("00000113", "00200008366300000200000000", "0163103000900000"); // standard cut Pi0 pp 7TeV, all photon qualities
  } else if (trainConfig == 4) { // additional standards added signals
    cuts.AddCutPCM("00000123", "00200009297002008250400000", "0163103000900000"); // standard cut LHC11h pp 2.76TeV
    cuts.AddCutPCM("00000123", "00200009227302008250404000", "0163101000000000"); // Ana eta analysis prefered
    cuts.AddCutPCM("00000123", "00200008366300000200000000", "0163103000900000"); // standard cut Pi0 pp 7TeV, all photon qualities
  } else if(trainConfig == 5){    // various standard cuts with rej pi0 window
    cuts.AddCutPCM("00000113", "00200009397300008250400000", "0163103b00900000"); // new pi0/eta cut 2.76TeV
    cuts.AddCutPCM("00000113", "00200009397300008250400000", "0163103b00000000"); // new pi0/eta cut 2.76TeV without MC smearing
    cuts.AddCutPCM("00000113", "00200009366300003800000000", "0163103b00900000"); // standard cut Pi0 pp 2.76TeV PbPb paper 2012

  // *********************************************************************************************************
  // 8 TeV  pp Run1
  // *********************************************************************************************************
  } else if (trainConfig == 100) {
    cuts.AddCutPCM("00010113", "00200009227300008250404000", "0152103000000000"); // New standard cut for 8TeV analysis V0AND with double counting cut, TOF removed
    cuts.AddCutPCM("00010013", "00200009227300008250404000", "0152103000000000"); // no SPD pileup cut
  } else if (trainConfig == 101) { //with rej pi0 window
    cuts.AddCutPCM("00010113", "00200009227300008250404000", "0152103b00000000"); // New standard cut for 8TeV analysis V0AND with double counting cut, TOF removed
    cuts.AddCutPCM("00010013", "00200009227300008250404000", "0152103b00000000"); // no SPD pileup cut

  // *********************************************************************************************************
  // 7 TeV  pp Run1
  // *********************************************************************************************************
  } else if (trainConfig == 200) {
    cuts.AddCutPCM("00000113", "00200009227300008250404000", "0152103000000000"); //New standard cut for 7TeV analysis V0OR with double counting cut, TOF removed
    cuts.AddCutPCM("00000013", "00200009227300008250404000", "0152103000000000"); // no SPD pileup cut
    cuts.AddCutPCM("00000113", "00200009227302008250400000", "0152103000000000"); //New standard cut for 7TeV analysis V0OR
    cuts.AddCutPCM("00000113", "00200008366300000200000000", "0163103000900000"); //Old standard cut for 7TeV analysis V0OR
    cuts.AddCutPCM("00000113", "00200009360300007800004000", "0263103000900000"); //dalitz: New Standard Only MB, standard pp7Tev cut dalitz
  } else if (trainConfig == 201) { //with rej pi0 window
    cuts.AddCutPCM("00000113", "00200009227300008250404000", "0152103b00000000"); //New standard cut for 7TeV analysis V0OR with double counting cut, TOF removed
    cuts.AddCutPCM("00000113", "00200009227302008250400000", "0152103b00000000"); //New standard cut for 7TeV analysis V0OR

  // *********************************************************************************************************
  // 5 TeV  pp Run2
  // *********************************************************************************************************
  } else if (trainConfig == 300) {
    cuts.AddCutPCM("00010113", "00200009227302008250400000", "0163103000000000"); //
  } else if (trainConfig == 301) {
    cuts.AddCutPCM("00010113", "00200009227300008250404000", "0163103000000000","1111100017032220000"); //INT7
    cuts.AddCutPCM("00052113", "00200009227300008250404000", "0163103000000000","1111100017032220000"); //EMC7
    cuts.AddCutPCM("00085113", "00200009227300008250404000", "0163103000000000","1111100017032220000"); //EG2
    cuts.AddCutPCM("00083113", "00200009227300008250404000", "0163103000000000","1111100017032220000"); //EG1
  } else if (trainConfig == 302) { //with rej pi0 window
    cuts.AddCutPCM("00010113", "00200009227300008250404000", "0163103b00000000","1111100017032220000"); //INT7
    cuts.AddCutPCM("00052113", "00200009227300008250404000", "0163103b00000000","1111100017032220000"); //EMC7
    cuts.AddCutPCM("00085113", "00200009227300008250404000", "0163103b00000000","1111100017032220000"); //EG2
    cuts.AddCutPCM("00083113", "00200009227300008250404000", "0163103b00000000","1111100017032220000"); //EG1

  // *********************************************************************************************************
  // 13 TeV  pp Run2
  // *********************************************************************************************************
  } else if (trainConfig == 400) {
    cuts.AddCutPCM("00010113", "00200009227300008250404000", "0163103000000000"); //INT7
  } else if (trainConfig == 401) {
    cuts.AddCutPCM("00010113", "00200009227300008250404000", "0163103000000000","1111100017032220000"); //INT7
    cuts.AddCutPCM("00052113", "00200009227300008250404000", "0163103000000000","1111100017032220000"); //EMC7
    cuts.AddCutPCM("00085113", "00200009227300008250404000", "0163103000000000","1111100017032220000"); //EG2
    cuts.AddCutPCM("00083113", "00200009227300008250404000", "0163103000000000","1111100017032220000"); //EG1
  } else if (trainConfig == 402) {
    cuts.AddCutPCM("00010113", "00200009227300008250404000", "0163103000000000","1111100067032220000"); //INT7
    cuts.AddCutPCM("00052113", "00200009227300008250404000", "0163103000000000","1111100067032220000"); //EMC7
    cuts.AddCutPCM("00085113", "00200009227300008250404000", "0163103000000000","1111100067032220000"); //EG2
    cuts.AddCutPCM("00083113", "00200009227300008250404000", "0163103000000000","1111100067032220000"); //EG1
  } else if (trainConfig == 403) { //with rej pi0 window
    cuts.AddCutPCM("00010113", "00200009227300008250404000", "0163103b00000000","1111100017032220000"); //INT7
    cuts.AddCutPCM("00052113", "00200009227300008250404000", "0163103b00000000","1111100017032220000"); //EMC7
    cuts.AddCutPCM("00085113", "00200009227300008250404000", "0163103b00000000","1111100017032220000"); //EG2
    cuts.AddCutPCM("00083113", "00200009227300008250404000", "0163103b00000000","1111100017032220000"); //EG1
  } else if (trainConfig == 404) { // with rej pi0 window
    cuts.AddCutPCM("00010113", "00200009227300008250404000", "0163103b00000000","1111100067032220000"); //INT7
    cuts.AddCutPCM("00052113", "00200009227300008250404000", "0163103b00000000","1111100067032220000"); //EMC7
    cuts.AddCutPCM("00085113", "00200009227300008250404000", "0163103b00000000","1111100067032220000"); //EG2
    cuts.AddCutPCM("00083113", "00200009227300008250404000", "0163103b00000000","1111100067032220000"); //EG1
  } else if (trainConfig == 405) {
    cuts.AddCutPCM("00010113", "00200009227300008250404000", "0163103000000000","1111100017032220000"); //INT7
    cuts.AddCutPCM("00010113", "00200009227300008250404000", "0163103b00000000","1111100017032220000"); //INT7
  } else if (trainConfig == 406) {
    cuts.AddCutPCM("00010113", "00200009227300008250404000", "0163103000000000","1111100067032220000"); //INT7
    cuts.AddCutPCM("00010113", "00200009227300008250404000", "0163103b00000000","1111100067032220000"); //INT7
  }  else if (trainConfig == 407) { //for eta prime
    cuts.AddCutPCM("00010113", "00200009227300008250404000", "01631030000000d0"); //INT7

  } else {
    Error(Form("HeavyNeutralMesonToGG_%i",trainConfig), "wrong trainConfig variable no cuts have been specified for the configuration");
    return;
  }

  if(!cuts.AreValid()){
    cout << "\n\n****************************************************" << endl;
    cout << "ERROR: No valid cuts stored in CutHandlerHeavyMesonConv! Returning..." << endl;
    cout << "****************************************************\n\n" << endl;
    return;
  }

  TList *EventCutList   = new TList();
  TList *ConvCutList    = new TList();
  TList *MesonCutList   = new TList();
  TList *ClusterCutList = new TList();

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
  AliConversionPhotonCuts **analysisCuts      = new AliConversionPhotonCuts*[numberOfCuts];
  MesonCutList->SetOwner(kTRUE);
  AliConversionMesonCuts **analysisMesonCuts  = new AliConversionMesonCuts*[numberOfCuts];
  ClusterCutList->SetOwner(kTRUE);
  AliCaloPhotonCuts **analysisClusterCuts     = new AliCaloPhotonCuts*[numberOfCuts];
  Bool_t enableClustersForTrigger             = kFALSE;

  for(Int_t i = 0; i<numberOfCuts; i++){
    analysisEventCuts[i] = new AliConvEventCuts();
    TString fitNamePi0            = Form("Pi0_Fit_Data_%s",energy.Data());
    TString fitNameEta            = Form("Eta_Fit_Data_%s",energy.Data());
    TString fAddedSignalString    = (cuts.GetEventCut(i)).Data();
    fAddedSignalString            = fAddedSignalString(6,1);
    Bool_t fAddedSignal           = kFALSE;
    if (fAddedSignalString.CompareTo("2") == 0)
      fAddedSignal                = kTRUE;

    TString mcInputNamePi0        = "";
    TString mcInputNameEta        = "";
    if (fAddedSignal && (generatorName.Contains("LHC12i3") || generatorName.CompareTo("LHC14e2b")==0)){
      mcInputNamePi0              = Form("Pi0_%s%s_addSig_%s", mcName.Data(), mcNameAdd.Data(), energy.Data() );
      mcInputNameEta              = Form("Eta_%s%s_addSig_%s", mcName.Data(), mcNameAdd.Data(), energy.Data() );
    } else {
      mcInputNamePi0              = Form("Pi0_%s%s_%s", mcName.Data(), mcNameAdd.Data(), energy.Data() );
      mcInputNameEta              = Form("Eta_%s%s_%s", mcName.Data(), mcNameAdd.Data(), energy.Data() );
    }
    if (doWeightingPart) analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE, fileNamePtWeights, mcInputNamePi0, mcInputNameEta, "",fitNamePi0,fitNameEta);


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

    if ( trainConfig == 301 || trainConfig == 401 || trainConfig == 402   ){
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

      enableClustersForTrigger  = kTRUE;
      analysisClusterCuts[i]    = new AliCaloPhotonCuts();
      analysisClusterCuts[i]->SetV0ReaderName(V0ReaderName);
      analysisClusterCuts[i]->SetCorrectionTaskSetting(corrTaskSetting);
      analysisClusterCuts[i]->SetLightOutput(enableLightOutput);
      analysisClusterCuts[i]->InitializeCutsFromCutString((cuts.GetClusterCut(i)).Data());
      ClusterCutList->Add(analysisClusterCuts[i]);
      analysisClusterCuts[i]->SetFillCutHistograms("");
    }


    analysisCuts[i] = new AliConversionPhotonCuts();
    analysisCuts[i]->SetV0ReaderName(V0ReaderName);
    if (enableLightOutput > 0) analysisCuts[i]->SetLightOutput(kTRUE);
    analysisCuts[i]->InitializeCutsFromCutString((cuts.GetPhotonCut(i)).Data());
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
    analysisCuts[i]->SetIsHeavyIon(isHeavyIon);
    ConvCutList->Add(analysisCuts[i]);
    analysisCuts[i]->SetFillCutHistograms("",kFALSE);

    analysisMesonCuts[i] = new AliConversionMesonCuts();
    if (enableLightOutput > 0) analysisMesonCuts[i]->SetLightOutput(kTRUE);
    analysisMesonCuts[i]->SetRunningMode(0);
    analysisMesonCuts[i]->InitializeCutsFromCutString((cuts.GetMesonCut(i)).Data());
    MesonCutList->Add(analysisMesonCuts[i]);
    analysisMesonCuts[i]->SetFillCutHistograms("");
    analysisEventCuts[i]->SetAcceptedHeader(HeaderList);
  }

  task->SetEventCutList(numberOfCuts,EventCutList);
  task->SetConversionCutList(numberOfCuts,ConvCutList);
  task->SetMesonCutList(numberOfCuts,MesonCutList);
  task->SetMoveParticleAccordingToVertex(kTRUE);
  task->SetCorrectionTaskSetting(corrTaskSetting);
  if (enableClustersForTrigger){
    task->SetDoClusterSelectionForTriggerNorm(enableClustersForTrigger);
    task->SetCaloCutList(numberOfCuts,ClusterCutList);
    task->SetDoClusterQA(1);  //Attention new switch small for Cluster QA
    if(enableExtMatchAndQA > 1){ task->SetPlotHistsExtQA(kTRUE);}
  }
  task->SetMesonType(selectedMeson);
  task->SetDoMesonQA(enableQAMesonTask); //Attention new switch for Pi0 QA
  task->SetDoPhotonQA(enableQAPhotonTask);  //Attention new switch small for Photon QA
  task->SetUseTHnSparse(enableTHnSparse);

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
