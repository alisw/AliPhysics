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
//($ALIPHYSICS/PWGGA/GammaConv/AliAnalysisTaskPi0v2Calo.cxx) for
//PbPb together with all supporting classes
//***************************************************************************************

//***************************************************************************************
//main function
//***************************************************************************************
void AddTask_Pi0v2Calo_PbPb(
  int     trainConfig                   = 1,        // change different set of cuts
  int     isMC                          = 0,        // run MC
  TString   photonCutNumberV0Reader     = "",       // 00000008400000000100000000 nom. B, 00000088400000000100000000 low B
  TString   periodNameV0Reader          = "",
  // general setting for task
  int     enableLightOutput             = 0,        // switch to run light output (only essential histograms for afterburner)
  int     enableTriggerMimicking        = 0,        // enable trigger mimicking
  bool    enableHeaderOverlap           = true,     // enable trigger overlap rejection
  TString   settingMaxFacPtHard         = "3.",     // maximum factor between hardest jet and ptHard generated
  int     debugLevel                    = 0,        // introducing debug levels for grid running
  // settings for weights
  // FPTW:fileNamePtWeights, FMUW:fileNameMultWeights, FCEF:fileNameCentFlattening, separate with ;
  TString   fileNameExternalInputs      = "",
  TString   periodName                  = "",       // name of the period, e.g. LHC15o
  TString   generatorName               = "DPMJET", // generator Name
  bool    enableMultiplicityWeighting   = false,    //
  TString   periodNameAnchor            = "",       //
  // special settings
  bool    enableSortingMCLabels         = true,     // enable sorting for MC cluster labels
  bool    doFlattening                  = false,    // switch on centrality flattening for LHC11h
  bool    doPrimaryTrackMatching        = true,     // enable basic track matching for all primary tracks to cluster
  bool    doInAndOutPlane               = false,    // enable in and out of plane method instead of dPhi(Pi0-Plane)
  bool    doFlowQA                      = false,    // enables flow QA histos
  // subwagon config
  TString   additionalTrainConfig       = "0"       // additional counter for trainconfig
  ) {

  AliCutHandlerPCM cuts;


  TString fileNamePtWeights     = cuts.GetSpecialFileNameFromString (fileNameExternalInputs, "FPTW:");
  TString fileNameMultWeights   = cuts.GetSpecialFileNameFromString (fileNameExternalInputs, "FMUW:");
  TString fileNameCentFlattening= cuts.GetSpecialFileNameFromString (fileNameExternalInputs, "FCEF:");

  TString addTaskName                 = "AddTask_Pi0v2Calo_PbPb";
  TString sAdditionalTrainConfig      = cuts.GetSpecialSettingFromAddConfig(additionalTrainConfig, "", "", addTaskName);
  if (sAdditionalTrainConfig.Atoi() > 0){
    trainConfig = trainConfig + sAdditionalTrainConfig.Atoi();
    cout << "INFO: " << addTaskName.Data() << " running additionalTrainConfig '" << sAdditionalTrainConfig.Atoi() << "', train config: '" << trainConfig << "'" << endl;
  }
  TString corrTaskSetting             = cuts.GetSpecialSettingFromAddConfig(additionalTrainConfig, "CF", "", addTaskName);
  if(corrTaskSetting.CompareTo(""))
    cout << "corrTaskSetting: " << corrTaskSetting.Data() << endl;

  int trackMatcherRunningMode = 0; // CaloTrackMatcher running mode
  TString strTrackMatcherRunningMode  = cuts.GetSpecialSettingFromAddConfig(additionalTrainConfig, "TM", "", addTaskName);
  if(additionalTrainConfig.Contains("TM"))
    trackMatcherRunningMode = strTrackMatcherRunningMode.Atoi();

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

  int localDebugFlag = 0;
  TString strLocalDebugFlag              = cuts.GetSpecialSettingFromAddConfig(additionalTrainConfig, "LOCALDEBUGFLAG", "", addTaskName);
  if(strLocalDebugFlag.Atoi()>0)
    localDebugFlag = strLocalDebugFlag.Atoi();

  TObjArray *rmaxFacPtHardSetting = settingMaxFacPtHard.Tokenize("_");
  if(rmaxFacPtHardSetting->GetEntries()<1){cout << "ERROR: AddTask_Pi0v2Calo_PbPb during parsing of settingMaxFacPtHard String '" << settingMaxFacPtHard.Data() << "'" << endl; return;}
  bool fMinPtHardSet        = false;
  Double_t minFacPtHard       = -1;
  bool fMaxPtHardSet        = false;
  Double_t maxFacPtHard       = 100;
  bool fSingleMaxPtHardSet  = false;
  Double_t maxFacPtHardSingle = 100;
  bool fJetFinderUsage      = false;
  bool fUsePtHardFromFile      = false;
  bool fUseAddOutlierRej      = false;
  for(int i = 0; i<rmaxFacPtHardSetting->GetEntries() ; i++){
    TObjString* tempObjStrPtHardSetting     = (TObjString*) rmaxFacPtHardSetting->At(i);
    TString strTempSetting                  = tempObjStrPtHardSetting->GetString();
    if(strTempSetting.BeginsWith("MINPTHFAC:")){
      strTempSetting.Replace(0,10,"");
      minFacPtHard               = strTempSetting.Atof();
      cout << "running with min pT hard jet fraction of: " << minFacPtHard << endl;
      fMinPtHardSet        = true;
    } else if(strTempSetting.BeginsWith("MAXPTHFAC:")){
      strTempSetting.Replace(0,10,"");
      maxFacPtHard               = strTempSetting.Atof();
      cout << "running with max pT hard jet fraction of: " << maxFacPtHard << endl;
      fMaxPtHardSet        = true;
    } else if(strTempSetting.BeginsWith("MAXPTHFACSINGLE:")){
      strTempSetting.Replace(0,16,"");
      maxFacPtHardSingle         = strTempSetting.Atof();
      cout << "running with max single particle pT hard fraction of: " << maxFacPtHardSingle << endl;
      fSingleMaxPtHardSet        = true;
    } else if(strTempSetting.BeginsWith("USEJETFINDER:")){
      strTempSetting.Replace(0,13,"");
      if(strTempSetting.Atoi()==1){
        cout << "using MC jet finder for outlier removal" << endl;
        fJetFinderUsage        = true;
      }
    } else if(strTempSetting.BeginsWith("PTHFROMFILE:")){
      strTempSetting.Replace(0,12,"");
      if(strTempSetting.Atoi()==1){
        cout << "using MC jet finder for outlier removal" << endl;
        fUsePtHardFromFile        = true;
      }
    } else if(strTempSetting.BeginsWith("ADDOUTLIERREJ:")){
      strTempSetting.Replace(0,14,"");
      if(strTempSetting.Atoi()==1){
        cout << "using path based outlier removal" << endl;
        fUseAddOutlierRej        = true;
      }
    } else if(rmaxFacPtHardSetting->GetEntries()==1 && strTempSetting.Atof()>0){
      maxFacPtHard               = strTempSetting.Atof();
      cout << "running with max pT hard jet fraction of: " << maxFacPtHard << endl;
      fMaxPtHardSet        = true;
    }
  }


  int isHeavyIon = 1;

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
  AliAnalysisTaskPi0v2Calo *task=NULL;
  task= new AliAnalysisTaskPi0v2Calo(Form("Pi0v2Calo_%i",trainConfig));
  task->SetIsHeavyIon(isHeavyIon);
  task->SetIsMC(isMC);
  task->SetV0ReaderName(V0ReaderName);
  task->SetCorrectionTaskSetting(corrTaskSetting);
  task->SetLightOutput(enableLightOutput);
  task->SetDoPrimaryTrackMatching(doPrimaryTrackMatching);
  task->SetTrackMatcherRunningMode(trackMatcherRunningMode);
  task->SetPeriodName(periodName);
  task->SetInAndOutPlane(doInAndOutPlane);

  // meson cuts
  // meson type (Dalitz or not), BG scheme, pool depth, rotation degrees, rapidity cut, radius cut, alpha, chi2, shared electrons, reject to close v0, MC smearing, dca, dca, dca

  // **********************************************************************************************************
  // ****************************** EMC configurations PbPb run 2 2018 pass 3 *********************************
  // **********************************************************************************************************
  // LHC18qr with 0-10 and 30-50 cent triggers
  if (trainConfig == 1){ // EMCAL clusters - 4 cent classes with 0-10 and 30-50 triggered -standard PbPb TM for secondaries
    cuts.AddCutCalo("10130e13","411790105ke30220000","01331031000000d0"); //  0-10%
    cuts.AddCutCalo("11310e13","411790105ke30220000","01331031000000d0"); // 10-30%
    cuts.AddCutCalo("13530e13","411790105ke30220000","01331031000000d0"); // 30-50%
    cuts.AddCutCalo("15910e13","411790105ke30220000","01331031000000d0"); // 50-90%
  } else if (trainConfig == 2){ // EMCAL clusters - cent -standard PbPb TM for secondaries
    cuts.AddCutCalo("10130e13","411790105ke30220000","01331031000000d0"); //  0-10%
  } else if (trainConfig == 3){ // EMCAL clusters - semi cent -standard PbPb TM for secondaries
    cuts.AddCutCalo("11310e13","411790105ke30220000","01331031000000d0"); // 10-30%
  } else if (trainConfig == 4){ // EMCAL clusters - semi peripheral -standard PbPb TM for secondaries
    cuts.AddCutCalo("13530e13","411790105ke30220000","01331031000000d0"); // 30-50%
  } else if (trainConfig == 5){ // EMCAL clusters - peripheral -standard PbPb TM for secondaries
    cuts.AddCutCalo("15910e13","411790105ke30220000","01331031000000d0"); // 50-90%

  // **********************************************************************************************************
  // ****************************** EMC configurations PbPb run 2 2015 pass 2 *********************************
  // **********************************************************************************************************
  } else if (trainConfig == 10){ // EMCAL clusters - 4 cent classes with 0-10 and 30-50 triggered -standard PbPb TM for secondaries
    cuts.AddCutCalo("10110e13","411790105ke30220000","01331031000000d0"); //  0-10%
    cuts.AddCutCalo("11310e13","411790105ke30220000","01331031000000d0"); // 10-30%
    cuts.AddCutCalo("13510e13","411790105ke30220000","01331031000000d0"); // 30-50%
    cuts.AddCutCalo("15910e13","411790105ke30220000","01331031000000d0"); // 50-90%
  } else if (trainConfig == 11){ // EMCAL clusters - cent -standard PbPb TM for secondaries
    cuts.AddCutCalo("10110e13","411790105ke30220000","01331031000000d0"); //  0-10%
  } else if (trainConfig == 12){ // EMCAL clusters - semi cent -standard PbPb TM for secondaries
    cuts.AddCutCalo("11310e13","411790105ke30220000","01331031000000d0"); // 10-30%
  } else if (trainConfig == 13){ // EMCAL clusters - semi peripheral -standard PbPb TM for secondaries
    cuts.AddCutCalo("13510e13","411790105ke30220000","01331031000000d0"); // 30-50%
  } else if (trainConfig == 14){ // EMCAL clusters - peripheral -standard PbPb TM for secondaries
    cuts.AddCutCalo("15910e13","411790105ke30220000","01331031000000d0"); // 50-90%

  // **********************************************************************************************************
  // *************************** EMC configurations PbPb run 2 primary track mult *****************************
  // **********************************************************************************************************
  } else if (trainConfig == 100){ // EMCAL+DCal clusters triggerd cents
    cuts.AddCutCalo("50130e13","411790105ke30220000","01331031000000d0"); //  00-10%
    cuts.AddCutCalo("53530e13","411790105ke30220000","01331031000000d0"); //  30-50%
  } else if (trainConfig == 101){ // EMCAL+DCal clusters
    cuts.AddCutCalo("50110e13","411790105ke30220000","01331031000000d0"); //  00-10%
    cuts.AddCutCalo("53510e13","411790105ke30220000","01331031000000d0"); //  30-50%
  } else if (trainConfig == 102){ // EMCAL+DCal clusters for MC testing
    cuts.AddCutCalo("60310e13","411790105ke30220000","01331031000000d0"); //  00-15%
  } else {
    Error(Form("Pi0v2Ca_%i",trainConfig), "wrong trainConfig variable no cuts have been specified for the configuration");
    return;
  }

  if(!cuts.AreValid()){
    cout << "\n\n****************************************************" << endl;
    cout << "ERROR: No valid cuts stored in CutHandlerCalo! Returning..." << endl;
    cout << "****************************************************\n\n" << endl;
    return;
  }
  int numberOfCuts = cuts.GetNCuts();

  TList *EventCutList   = new TList();
  TList *ClusterCutList = new TList();
  TList *MesonCutList   = new TList();


  EventCutList->SetOwner(true);
  AliConvEventCuts **analysisEventCuts        = new AliConvEventCuts*[numberOfCuts];
  ClusterCutList->SetOwner(true);
  AliCaloPhotonCuts **analysisClusterCuts     = new AliCaloPhotonCuts*[numberOfCuts];
  MesonCutList->SetOwner(true);
  AliConversionMesonCuts **analysisMesonCuts  = new AliConversionMesonCuts*[numberOfCuts];

  for(int i = 0; i<numberOfCuts; i++){
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
      if (enableLightOutput == 4) fTrackMatcher->SetLightOutput(true);
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
      analysisEventCuts[i]->SetUseWeightMultiplicityFromFile(true, fileNameMultWeights, dataInputMultHisto, mcInputMultHisto );
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
      analysisEventCuts[i]->SetUseJetFinderForOutliers(true);
    if(fUsePtHardFromFile)
      analysisEventCuts[i]->SetUsePtHardBinFromFile(true);
    if(fUseAddOutlierRej)
      analysisEventCuts[i]->SetUseAdditionalOutlierRejection(true);
    analysisEventCuts[i]->InitializeCutsFromCutString((cuts.GetEventCut(i)).Data());
    EventCutList->Add(analysisEventCuts[i]);
    analysisEventCuts[i]->SetCorrectionTaskSetting(corrTaskSetting);
    analysisEventCuts[i]->SetFillCutHistograms("",false);
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

    analysisMesonCuts[i]    = new AliConversionMesonCuts();
    if (enableLightOutput > 0 && enableLightOutput != 4) analysisMesonCuts[i]->SetLightOutput(1);
    if (enableLightOutput == 4) analysisMesonCuts[i]->SetLightOutput(2);
    analysisMesonCuts[i]->InitializeCutsFromCutString((cuts.GetMesonCut(i)).Data());
    analysisMesonCuts[i]->SetIsMergedClusterCut(2);
    analysisMesonCuts[i]->SetCaloMesonCutsObject(analysisClusterCuts[i]);
    MesonCutList->Add(analysisMesonCuts[i]);
    analysisMesonCuts[i]->SetFillCutHistograms("");

    if(analysisMesonCuts[i]->DoGammaSwappForBg()) analysisClusterCuts[i]->SetUseEtaPhiMapForBackCand(true);
    analysisClusterCuts[i]->SetFillCutHistograms("");
  }
  task->SetAllowOverlapHeaders(enableHeaderOverlap);
  task->SetEventCutList(numberOfCuts,EventCutList);
  task->SetCaloCutList(numberOfCuts,ClusterCutList);
  task->SetMesonCutList(numberOfCuts,MesonCutList);
  task->SetDoMesonAnalysis(true);
  task->SetCorrectionTaskSetting(corrTaskSetting);
  task->SetEnableSortingOfMCClusLabels(enableSortingMCLabels);
  task->SetFlowQa(doFlowQA);

  //connect containers
  AliAnalysisDataContainer *coutput =
    mgr->CreateContainer(!(corrTaskSetting.CompareTo("")) ? Form("Pi0v2Calo_%i",trainConfig) : Form("Pi0v2Calo_%i_%s",trainConfig,corrTaskSetting.Data()), TList::Class(),
              AliAnalysisManager::kOutputContainer,Form("Pi0v2Ca_%i.root",trainConfig) );

  mgr->AddTask(task);
  mgr->ConnectInput(task,0,cinput);
  mgr->ConnectOutput(task,1,coutput);

  return;
}
