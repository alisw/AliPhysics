/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                       *
 * Author: Steven Merkel, Adrian Mechler                     *
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
void AddTask_SigmaPlToProtonPiZeroAOD(
  Int_t     trainConfig                   = 1,        // change different set of cuts
  Int_t     isMC                          = 0,        // run MC
  TString   photonCutNumberV0Reader       = "00000008400100001500000000",       // 00000008400000000100000000 nom. B, 00000088400000000100000000 low B
  TString   periodNameV0Reader            = "LHC18b",
  // general setting for task
  Int_t     enableExtMatchAndQA           = 0,        // disabled (0), extMatch (1), extQA_noCellQA (2), extMatch+extQA_noCellQA (3), extQA+cellQA (4), extMatch+extQA+cellQA (5)
  Int_t     enableLightOutput             = kTRUE,   // switch to run light output (only essential histograms for afterburner)
  Int_t     enableTriggerMimicking        = 0,        // enable trigger mimicking
  Bool_t    enableTriggerOverlapRej       = kFALSE,   // enable trigger overlap rejection
  // settings for weights
  // FPTW:fileNamePtWeights, FMUW:fileNameMultWeights, separate with ;
  TString   fileNameExternalInputs        = "FMUW:fileNameMultWeights;FEPC:fileNamedEdxPostCalib",
  Int_t     doWeightingPart               = 0,        // enable Weighting
  TString   generatorName                 = "LHC18b", //"DPMJET", // generator Name
  Bool_t    enableMultiplicityWeighting   = kFALSE,   //
  TString   periodNameAnchor              = "LHC18b",       //
  // subwagon config
  TString   additionalTrainConfig         = "0"       // additional counter for trainconfig)
  )
{
   AliCutHandlerPCM cuts(13);


  TString addTaskName                 = "AddTask_SigmaPlToProtonPiZeroAOD";
    TString sAdditionalTrainConfig      = cuts.GetSpecialSettingFromAddConfig(additionalTrainConfig, "", "", addTaskName);
  if (sAdditionalTrainConfig.Atoi() > 0){
    trainConfig = trainConfig + sAdditionalTrainConfig.Atoi();
    cout << "INFO: " << addTaskName.Data() << " running additionalTrainConfig '" << sAdditionalTrainConfig.Atoi() << "', train config: '" << trainConfig << "'" << endl;
  }

  TString fileNamePtWeights           = cuts.GetSpecialFileNameFromString (fileNameExternalInputs, "FPTW:");
  TString fileNameMultWeights         = cuts.GetSpecialFileNameFromString (fileNameExternalInputs, "FMUW:");
  TString fileNameCustomTriggerMimicOADB   = cuts.GetSpecialFileNameFromString (fileNameExternalInputs, "FTRM:");

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

  Int_t isHeavyIon = 0;

  // ================== GetAnalysisManager ===============================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    // Error(Form("%s_%i", addTaskName.Data(),  trainConfig), "No analysis manager found.");
    return;
  }

  // // // // ================== GetInputEventHandler =============================
  AliVEventHandler *inputHandler=mgr->GetInputEventHandler();

  // //=========  Set Cutnumber for V0Reader ================================
  TString cutnumberPhoton = photonCutNumberV0Reader.Data();
  TString cutnumberEvent = "00000003";
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();

  //   //========= Add V0 Reader to  ANALYSIS manager if not yet existent =====
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
  AliAnalysisTaskSigmaPlToProtonPiZeroAOD* task = NULL;
  task = new AliAnalysisTaskSigmaPlToProtonPiZeroAOD(addTaskName.Data());
  task->SelectCollisionCandidates(AliVEvent::kAnyINT);
  task->SetIsHeavyIon(isHeavyIon);
  task->SetIsMC(isMC);
  task->SetV0ReaderName(V0ReaderName);
  task->SetCorrectionTaskSetting(corrTaskSetting);
  task->SetTrackMatcherRunningMode(trackMatcherRunningMode);

  // here is the order of the cluster cut string
  // usually for EMCal we start with 11111: default values for            "ClusterType", "EtaMin", "EtaMax", "PhiMin", "PhiMax"
  // then two numbers for nonlinearity, e.g. 21: this is                  "NonLinearity1", "NonLinearity2"
  // Then some cuts on the clusters, e.g. 06003222: this is               "DistanceToBadChannel", "Timing", "TrackMatching", "ExoticCell", "MinEnergy", "MinNCells", "MinM02", "MaxM02"
  // finally some for now unused cuts, usually 0000: this is              "MinMaxM20", "RecConv", "MaximumDispersion", "NLM"


  // *****************************************************************************************************
  // ******************** pp 13 TeV cuts             *****************************************************
  // *****************************************************************************************************
  if (trainConfig == 1){   // PHOS Standard
    cuts.AddCutHeavyMesonSigma("00010113","24466190ra09cc00000","0163103100000010","20120003211311011"); 
  } else if (trainConfig == 2){   // PHOS PID
    cuts.AddCutHeavyMesonSigma("00010113","24466190ra09cc00000","0163103100000010","20220003211311011"); // NCluster TPC 50
    cuts.AddCutHeavyMesonSigma("00010113","24466190ra09cc00000","0163103100000010","20320003211311011"); // NCluster TPC 65
    cuts.AddCutHeavyMesonSigma("00010113","24466190ra09cc00000","0163103100000010","20110003211311011"); // Nchi2 TPC 4
  } else if (trainConfig == 3){   // PHOS
    cuts.AddCutHeavyMesonSigma("00010113","24466190sa09cc00000","0163103100000010","20120003211311011");
    cuts.AddCutHeavyMesonSigma("00010113","244661904a09cc00000","0163103100000010","20120003211311011");
    cuts.AddCutHeavyMesonSigma("00010113","24466190ra01cc00000","0163103100000010","20120003211311011");
    cuts.AddCutHeavyMesonSigma("00010113","24466190ra07cc00000","0163103100000010","20120003211311011");
  } else if (trainConfig == 4){   // PHOS DCA-Variation
    cuts.AddCutHeavyMesonSigma("00010113","24466190ra09cc00000","0163103100000010","20120002211311011"); // DCAXY 0.02
    cuts.AddCutHeavyMesonSigma("00010113","24466190ra09cc00000","0163103100000010","20120004211311011"); // DCAXY 0.015
  } else if (trainConfig == 5){   // PHOS Pion mass-Variation
    cuts.AddCutHeavyMesonSigma("00010113","24466190ra09cc00000","0163103100000010","20120003211322011"); // tighter mass cut
    cuts.AddCutHeavyMesonSigma("00010113","24466190ra09cc00000","0163103100000010","20120003211342011"); // on left side tighter
    cuts.AddCutHeavyMesonSigma("00010113","24466190ra09cc00000","0163103100000010","20120003211333011"); // looser mass cut
  } else if (trainConfig == 6){   // PHOS Opening angle-Variation
    cuts.AddCutHeavyMesonSigma("00010113","24466190ra09cc00000","0163103100000010","20120003211311021"); // wider cut
    cuts.AddCutHeavyMesonSigma("00010113","24466190ra09cc00000","0163103100000010","20120003211311031"); // tighter cut
    cuts.AddCutHeavyMesonSigma("00010113","24466190ra09cc00000","0163103100000010","20120003211311041"); // tightest cut
  } else if (trainConfig == 7){   // Background describtion
    cuts.AddCutHeavyMesonSigma("00010113","24466190ra09cc00000","0163103100000010","20120003211311010"); // changed to pion rotation  
  } else if (trainConfig == 10){    // PHI 7
    cuts.AddCutHeavyMesonSigma("00062113","24466190sa09cc00000","0163103100000010","20120003211311011");



  } else if (trainConfig == 100){    // EDC
    cuts.AddCutHeavyMesonSigma("00010113","4117912067032230000","0163103100000050","20120003211311011");
  } else if (trainConfig == 101){    // PHI7
    cuts.AddCutHeavyMesonSigma("0008e113","4117912067032230000","0163103100000010","20120003211311011");
    cuts.AddCutHeavyMesonSigma("0008d113","4117912067032230000","0163103100000010","20120003211311011");
  }

  if(!cuts.AreValid()){
    cout << "\n\n****************************************************" << endl;
    cout << "ERROR: No valid cuts stored in CutHandlerCalo! Returning..." << endl;
    cout << "****************************************************\n\n" << endl;
    return;
  }

  Int_t numberOfCuts = cuts.GetNCuts();

  TList *EventCutList = new TList();
  TList *ClusterCutList = new TList();
  TList *MesonCutList = new TList();
  TList *SigmaCutList = new TList();

  TList *HeaderList = new TList();

  TString energy = "";
  TString mcName = "";
  TString mcNameAdd = "";
  if (generatorName.Contains("WOSDD")){
    mcNameAdd = "_WOSDD";
  } else if (generatorName.Contains("WSDD")){
    mcNameAdd = "_WSDD";
  }


  EventCutList->SetOwner(kTRUE);
  AliConvEventCuts **analysisEventCuts = new AliConvEventCuts*[numberOfCuts];
  ClusterCutList->SetOwner(kTRUE);
  AliCaloPhotonCuts **analysisClusterCuts = new AliCaloPhotonCuts*[numberOfCuts];
  MesonCutList->SetOwner(kTRUE);
  AliConversionMesonCuts **analysisMesonCuts = new AliConversionMesonCuts*[numberOfCuts];
  SigmaCutList->SetOwner(kTRUE);
  AliCaloSigmaCuts **analysisSigmaCuts = new AliCaloSigmaCuts*[numberOfCuts];

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

    // definition of weighting input
    // TString fitNamePi0 = Form("Pi0_Fit_Data_%s",energy.Data());
    // TString fitNameEta = Form("Eta_Fit_Data_%s",energy.Data());
    // TString fAddedSignalString = cuts.GetEventCut(i);
    // fAddedSignalString = fAddedSignalString(6,1);
    // Bool_t fAddedSignal = kFALSE;
    // if (fAddedSignalString.CompareTo("2") == 0) fAddedSignal = kTRUE;

    // TString mcInputNamePi0 = "";
    // TString mcInputNameEta = "";
    // if (fAddedSignal && (generatorName.Contains("LHC12i3") || generatorName.CompareTo("LHC14e2b")==0)){
    //   mcInputNamePi0 = Form("Pi0_%s%s_addSig_%s", mcName.Data(), mcNameAdd.Data(), energy.Data() );
    //   mcInputNameEta = Form("Eta_%s%s_addSig_%s", mcName.Data(), mcNameAdd.Data(), energy.Data() );
    // } else {
    //   mcInputNamePi0 = Form("Pi0_%s%s_%s", mcName.Data(), mcNameAdd.Data(), energy.Data() );
    //   mcInputNameEta = Form("Eta_%s%s_%s", mcName.Data(), mcNameAdd.Data(), energy.Data() );
    // }

    // if (doWeightingPart) analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE, fileNamePtWeights, mcInputNamePi0, mcInputNameEta, "",fitNamePi0,fitNameEta);

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

    // if (enableMultiplicityWeighting) analysisEventCuts[i]->SetUseWeightMultiplicityFromFile( kTRUE, fileNameMultWeights, dataInputMultHisto, mcInputMultHisto );

    analysisEventCuts[i]->SetTriggerMimicking(enableTriggerMimicking);
    if(fileNameCustomTriggerMimicOADB.CompareTo("") != 0)
      analysisEventCuts[i]->SetCustomTriggerMimicOADBFile(fileNameCustomTriggerMimicOADB);
    analysisEventCuts[i]->SetTriggerOverlapRejecion(enableTriggerOverlapRej);
    analysisEventCuts[i]->SetV0ReaderName(V0ReaderName);
    analysisEventCuts[i]->SetCorrectionTaskSetting(corrTaskSetting);
    if (enableLightOutput > 0) analysisEventCuts[i]->SetLightOutput(kTRUE);
    if(trainConfig == 447) analysisEventCuts[i]->SetUseSphericityTrue(kTRUE);
    if (periodNameV0Reader.CompareTo("") != 0) analysisEventCuts[i]->SetPeriodEnum(periodNameV0Reader);
    analysisEventCuts[i]->InitializeCutsFromCutString((cuts.GetEventCut(i)).Data());
    EventCutList->Add(analysisEventCuts[i]);
    if (enableLightOutput > 0) analysisEventCuts[i]->SetFillCutHistograms("",kFALSE);
    else analysisEventCuts[i]->SetFillCutHistograms("",kTRUE);

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
    // analysisMesonCuts[i]->SetEnableSigmaOpenCut(kTRUE);                       // enables the overloaded DCAGammaGamma to work as an AP like cut
    analysisMesonCuts[i]->InitializeCutsFromCutString((cuts.GetMesonCut(i)).Data());
    analysisMesonCuts[i]->SetIsMergedClusterCut(2);
    analysisMesonCuts[i]->SetCaloMesonCutsObject(analysisClusterCuts[i]);
    MesonCutList->Add(analysisMesonCuts[i]);
    analysisMesonCuts[i]->SetFillCutHistograms("");


    analysisSigmaCuts[i] = new AliCaloSigmaCuts();
    analysisSigmaCuts[i]->InitializeCutsFromCutString((cuts.GetSigmaCut(i)).Data());
    SigmaCutList->Add(analysisSigmaCuts[i]);
    analysisSigmaCuts[i]->SetFillCutHistograms("");
    analysisEventCuts[i]->SetAcceptedHeader(HeaderList);

   
  }
  task->SetEventCutList(numberOfCuts,EventCutList);
  task->SetCaloCutList(numberOfCuts,ClusterCutList);
  task->SetMesonCutList(numberOfCuts,MesonCutList);
  task->SetSigmaCutList(numberOfCuts,SigmaCutList);
  // task->SetCorrectionTaskSetting(corrTaskSetting);


  //connect containers
    AliAnalysisDataContainer *coutput =
    mgr->CreateContainer(!(corrTaskSetting.CompareTo("")) ? Form("SigmaPlus%i",trainConfig) : Form("SigmaPlus%i_%s",trainConfig,corrTaskSetting.Data()), TList::Class(),
              AliAnalysisManager::kOutputContainer,Form("SigmaPlus%i.root",trainConfig));
  mgr->AddTask(task);
  mgr->ConnectInput(task,0,cinput);
  mgr->ConnectOutput(task,1,coutput);

  return;

}
