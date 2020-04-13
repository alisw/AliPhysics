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
//($ALIPHYSICS/PWGGA/GammaConv/AliAnalysisTaskGammaConvCaloIso.cxx) for
//pp together with all supporting classes
//***************************************************************************************

//***************************************************************************************
//main function
//***************************************************************************************
void AddTask_GammaConvCaloIso_pp(
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
  Float_t   maxFacPtHard                  = 3.,       // maximum factor between hardest jet and ptHard generated
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

  TString fileNamePtWeights     = cuts.GetSpecialFileNameFromString (fileNameExternalInputs, "FPTW:");
  TString fileNameMultWeights   = cuts.GetSpecialFileNameFromString (fileNameExternalInputs, "FMUW:");

  TString addTaskName                 = "AddTask_GammaConvCaloIso_pp";
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

  Int_t isHeavyIon = 0;

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
  AliAnalysisTaskGammaConvCaloIso *task=NULL;
  task= new AliAnalysisTaskGammaConvCaloIso(Form("GammaConvCaloIso_%i",trainConfig));
  task->SetIsHeavyIon(isHeavyIon);
  task->SetIsMC(isMC);
  task->SetV0ReaderName(V0ReaderName);
  if (enableLightOutput > 1) task->SetLightOutput(kTRUE);
  task->SetDoPrimaryTrackMatching(doPrimaryTrackMatching);
  task->SetTrackMatcherRunningMode(trackMatcherRunningMode);

  // cluster cuts
  // 0 "ClusterType",  1 "EtaMin", 2 "EtaMax", 3 "PhiMin", 4 "PhiMax", 5 "DistanceToBadChannel", 6 "Timing", 7 "TrackMatching", 8 "ExoticCell",
  // 9 "MinEnergy", 10 "MinNCells", 11 "MinM02", 12 "MaxM02", 13 "MinMaxM20", 14 "RecConv", 15 "MaximumDispersion", 16 "NLM"

  // *****************************************************************************************************
  // pp 13 TeV setup, Isolation
  // *****************************************************************************************************
  //PhotonIsolation Variations, LHC17i3 MC without non linearity
  if (trainConfig == 33){ // without TrackMatching
    cuts.AddCutPCMCalo("00083113","00200009327000008250400000","1111100067032230000","0163103100000010"); //without TrackMatching
   /* cuts.AddCutPCMCalo("00083113","00200009327000008250400000","b221100060032230000","0163103100000010"); //without TrackMatching
    cuts.AddCutPCMCalo("00083113","00200009327000008250400000","c221100060032230000","0163103100000010"); //without TrackMatching
    cuts.AddCutPCMCalo("00083113","00200009327000008250400000","d221100060032230000","0163103100000010"); //without TrackMatching
    cuts.AddCutPCMCalo("00083113","00200009327000008250400000","e221100060032230000","0163103100000010"); //without TrackMatching
*/
  } else {
    Error(Form("GammaConvCaloIso_%i",trainConfig), "wrong trainConfig variable no cuts have been specified for the configuration");
    return;
  }

  if(!cuts.AreValid()){
    cout << "\n\n****************************************************" << endl;
    cout << "ERROR: No valid cuts stored in CutHandlerConvCalo! Returning..." << endl;
    cout << "****************************************************\n\n" << endl;
    return;
  }

  Int_t numberOfCuts = cuts.GetNCuts();

  TList *EventCutList = new TList();
  TList *ConvCutList = new TList();
  TList *ClusterCutList = new TList();
  TList *MesonCutList = new TList();

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
  if (generatorName.Contains("LHC17i3a1")){
    energy = "13TeV";
    mcName = "Pythia8_LHC17i3a1";
  } else if (generatorName.Contains("LHC17i3c1")){
      energy = "13TeV";
      mcName = "Pythia8_LHC17i3c1";
  }


  EventCutList->SetOwner(kTRUE);
  AliConvEventCuts **analysisEventCuts = new AliConvEventCuts*[numberOfCuts];
  ConvCutList->SetOwner(kTRUE);
  AliConversionPhotonCuts **analysisCuts = new AliConversionPhotonCuts*[numberOfCuts];
  ClusterCutList->SetOwner(kTRUE);
  AliCaloPhotonCuts **analysisClusterCuts = new AliCaloPhotonCuts*[numberOfCuts];
  MesonCutList->SetOwner(kTRUE);
  AliConversionMesonCuts **analysisMesonCuts = new AliConversionMesonCuts*[numberOfCuts];

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

    //create AliPhotonIsolation instance, if there is none present
    TString PhotonIsolationName = "PhotonIsolation";
    if( !(AliPhotonIsolation*)mgr->GetTask(PhotonIsolationName.Data()) ){
      AliPhotonIsolation* fPhotonIsolation = new AliPhotonIsolation(PhotonIsolationName.Data());
     // fPhotonIsolation->SetV0ReaderName(V0ReaderName);
     // fPhotonIsolation->SetCorrectionTaskSetting(corrTaskSetting);
      mgr->AddTask(fPhotonIsolation);
      mgr->ConnectInput(fPhotonIsolation,0,cinput);
    }

    analysisEventCuts[i] = new AliConvEventCuts();

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

    TString dataInputMultHisto  = "";
    TString mcInputMultHisto    = "";
    TString triggerString   = cuts.GetEventCut(i);
    triggerString           = triggerString(3,2);
    if (triggerString.CompareTo("03")==0)
      triggerString         = "00";
    if (periodNameAnchor.CompareTo("LHC13g") == 0 && triggerString.CompareTo("10")== 0 )
      triggerString         = "00";

    dataInputMultHisto      = Form("%s_%s", periodNameAnchor.Data(), triggerString.Data());
    mcInputMultHisto        = Form("%s_%s", periodNameV0Reader.Data(), triggerString.Data());

    if (enableMultiplicityWeighting){
      cout << "enableling mult weighting" << endl;
      analysisEventCuts[i]->SetUseWeightMultiplicityFromFile( kTRUE, fileNameMultWeights, dataInputMultHisto, mcInputMultHisto );
    }

    analysisEventCuts[i]->SetTriggerMimicking(enableTriggerMimicking);
    analysisEventCuts[i]->SetTriggerOverlapRejecion(enableTriggerOverlapRej);
    analysisEventCuts[i]->SetMaxFacPtHard(maxFacPtHard);
    analysisEventCuts[i]->SetV0ReaderName(V0ReaderName);
    analysisEventCuts[i]->SetCorrectionTaskSetting(corrTaskSetting);
    if (periodNameV0Reader.CompareTo("") != 0) analysisEventCuts[i]->SetPeriodEnum(periodNameV0Reader);
    if (enableLightOutput > 0 && enableLightOutput < 3) analysisEventCuts[i]->SetLightOutput(kTRUE);
    analysisEventCuts[i]->InitializeCutsFromCutString((cuts.GetEventCut(i)).Data());
    EventCutList->Add(analysisEventCuts[i]);
    analysisEventCuts[i]->SetFillCutHistograms("",kFALSE);

    analysisCuts[i] = new AliConversionPhotonCuts();
    analysisCuts[i]->SetV0ReaderName(V0ReaderName);
    if (enableLightOutput > 0 && enableLightOutput < 3) analysisCuts[i]->SetLightOutput(kTRUE);
    analysisCuts[i]->InitializeCutsFromCutString((cuts.GetPhotonCut(i)).Data());
    analysisCuts[i]->SetIsHeavyIon(isHeavyIon);
    ConvCutList->Add(analysisCuts[i]);
    analysisCuts[i]->SetFillCutHistograms("",kFALSE);

    analysisClusterCuts[i] = new AliCaloPhotonCuts(isMC);
    analysisClusterCuts[i]->SetHistoToModifyAcceptance(histoAcc);
    analysisClusterCuts[i]->SetV0ReaderName(V0ReaderName);
    analysisClusterCuts[i]->SetCorrectionTaskSetting(corrTaskSetting);
    analysisClusterCuts[i]->SetCaloTrackMatcherName(TrackMatcherName);
    if (enableLightOutput > 0 && enableLightOutput < 3) analysisClusterCuts[i]->SetLightOutput(kTRUE);
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
    if(doSmear) analysisMesonCuts[i]->SetDefaultSmearing(bremSmear,smearPar,smearParConst);
  }

  task->SetEventCutList(numberOfCuts,EventCutList);
  task->SetConversionCutList(numberOfCuts,ConvCutList);
  task->SetCaloCutList(numberOfCuts,ClusterCutList);
  task->SetMesonCutList(numberOfCuts,MesonCutList);
  task->SetMoveParticleAccordingToVertex(kTRUE);
  task->SetDoMesonAnalysis(kTRUE);
  task->SetCorrectionTaskSetting(corrTaskSetting);
  task->SetDoMesonQA(enableQAMesonTask); //Attention new switch for Pi0 QA
  task->SetDoPhotonQA(enableQAPhotonTask);  //Attention new switch small for Photon QA
  task->SetDoClusterQA(1);  //Attention new switch small for Cluster QA
  task->SetUseTHnSparse(enableTHnSparse);
  task->SetDoTreeInvMassShowerShape(doTreeClusterShowerShape);
  task->SetEnableSortingOfMCClusLabels(enableSortingMCLabels);
  if(enableTreeConvGammaShape) task->SetDoTreeConvGammaShowerShape(kTRUE);
  if(enableExtMatchAndQA > 1){ task->SetPlotHistsExtQA(kTRUE);}

  //connect containers
  AliAnalysisDataContainer *coutput =
    mgr->CreateContainer(!(corrTaskSetting.CompareTo("")) ? Form("GammaConvCaloIso_%i",trainConfig) : Form("GammaConvCaloIso_%i_%s",trainConfig,corrTaskSetting.Data()), TList::Class(),
              AliAnalysisManager::kOutputContainer, Form("GammaConvCaloIso_%i.root",trainConfig) );

  mgr->AddTask(task);
  mgr->ConnectInput(task,0,cinput);
  mgr->ConnectOutput(task,1,coutput);
  if(enableQAPhotonTask>1){
    for(Int_t i = 0; i<numberOfCuts; i++){
      mgr->ConnectOutput(task,2+i,mgr->CreateContainer(Form("%s_%s_%s_%s Photon DCA tree",(cuts.GetEventCut(i)).Data(),(cuts.GetPhotonCut(i)).Data(),cuts.GetClusterCut(i).Data(),(cuts.GetMesonCut(i)).Data()), TTree::Class(), AliAnalysisManager::kOutputContainer, Form("GammaConvCaloIso_%i.root",trainConfig)) );
    }
  }

  return;

}
