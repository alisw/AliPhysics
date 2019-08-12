/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: Friederike Bock                                                *
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
//($ALIPHYSICS/PWGGA/GammaConv/AliAnalysisTaskGammaTriggerQA.cxx) for
//PbPb together with all supporting classes
//***************************************************************************************

//***************************************************************************************
//CutHandler contains all cuts for a certain analysis and trainconfig,
//it automatically checks length of cutStrings and takes care of the number of added cuts,
//no specification of the variable 'numberOfCuts' needed anymore.
//***************************************************************************************
class CutHandlerTriggerQA{
  public:
    CutHandlerTriggerQA(Int_t nMax=10){
      nCuts=0; nMaxCuts=nMax; validCuts = true;
      eventCutArray = new TString[nMaxCuts]; clusterCutArray = new TString[nMaxCuts];
      for(Int_t i=0; i<nMaxCuts; i++) {eventCutArray[i] = ""; clusterCutArray[i] = "";}
    }

    void AddCut(TString eventCut, TString clusterCut){
      if(nCuts>=nMaxCuts) {cout << "ERROR in CutHandlerTriggerQA: Exceeded maximum number of cuts!" << endl; validCuts = false; return;}
      if( eventCut.Length()!=8 || clusterCut.Length()!=19  ) {cout << "ERROR in CutHandlerTriggerQA: Incorrect length of cut string!" << endl; validCuts = false; return;}
      eventCutArray[nCuts]=eventCut; clusterCutArray[nCuts]=clusterCut;
      nCuts++;
      return;
    }
    Bool_t AreValid(){return validCuts;}
    Int_t GetNCuts(){if(validCuts) return nCuts; else return 0;}
    TString GetEventCut(Int_t i){if(validCuts&&i<nMaxCuts&&i>=0) return eventCutArray[i]; else{cout << "ERROR in CutHandlerTriggerQA: GetEventCut wrong index i" << endl;return "";}}
    TString GetClusterCut(Int_t i){if(validCuts&&i<nMaxCuts&&i>=0) return clusterCutArray[i]; else {cout << "ERROR in CutHandlerTriggerQA: GetClusterCut wrong index i" << endl;return "";}}
  private:
    Bool_t validCuts;
    Int_t nCuts; Int_t nMaxCuts;
    TString* eventCutArray;
    TString* clusterCutArray;
};

//***************************************************************************************
//main function
//***************************************************************************************
void AddTask_GammaTriggerQA_PbPb( Int_t     trainConfig                     = 1,                    // change different set of cuts
                                  Int_t     isMC                            = 0,                    // run MC
                                  TString   cutnumberAODBranch              = "111110006008400000001500000",
                                  TString   periodName                      = "LHC13d2",            // name of the period for added signals and weighting
                                  TString   periodNameV0Reader              = "",                   // name of period for V0Reader
                                  Int_t     enableExtMatchAndQA             = 0,                    // disabled (0), extMatch (1), extQA_noCellQA (2), extMatch+extQA_noCellQA (3), extQA+cellQA (4), extMatch+extQA+cellQA (5)
                                  Bool_t    runLightOutput                  = kFALSE,               // switch to run light output (only essential histograms for afterburner)
                                  Bool_t    doFlattening                    = kFALSE,               // switch on centrality flattening for LHC11h
                                  Int_t     QAflag                          = 0,                    // switch for QA
                                  TString   fileNameInputForCentFlattening  = "",                   // file name for centrality flattening
                                  TString   additionalTrainConfig           = "0"                   // additional counter for trainconfig
                ) {

  TString corrTaskSetting = ""; // select which correction task setting to use
  //parse additionalTrainConfig flag
  TObjArray *rAddConfigArr = additionalTrainConfig.Tokenize("_");
  if(rAddConfigArr->GetEntries()<1){cout << "ERROR: AddTask_GammaTriggerQA_PbPb during parsing of additionalTrainConfig String '" << additionalTrainConfig.Data() << "'" << endl; return;}
  TObjString* rAdditionalTrainConfig;
  for(Int_t i = 0; i<rAddConfigArr->GetEntries() ; i++){
    if(i==0) rAdditionalTrainConfig = (TObjString*)rAddConfigArr->At(i);
    else{
      TObjString* temp = (TObjString*) rAddConfigArr->At(i);
      TString tempStr = temp->GetString();
    }
  }
  TString sAdditionalTrainConfig = rAdditionalTrainConfig->GetString();
  if (sAdditionalTrainConfig.Atoi() > 0){
    trainConfig = trainConfig + sAdditionalTrainConfig.Atoi();
    cout << "INFO: AddTask_GammaTriggerQA_PbPb running additionalTrainConfig '" << sAdditionalTrainConfig.Atoi() << "', train config: '" << trainConfig << "'" << endl;
  }
  cout << "corrTaskSetting: " << corrTaskSetting.Data() << endl;

  Int_t isHeavyIon = 1;

  // ================== GetAnalysisManager ===============================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error(Form("AddTask_GammaTriggerQA_PbPb_%i",trainConfig), "No analysis manager found.");
    return ;
  }

  // ================== GetInputEventHandler =============================
  AliVEventHandler *inputHandler=mgr->GetInputEventHandler();

  Bool_t isMCForOtherSettings = 0;
  if (isMC > 0) isMCForOtherSettings = 1;
  //========= Add PID Reponse to ANALYSIS manager ====
  if(!(AliPIDResponse*)mgr->GetTask("PIDResponseTask")){
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
    AddTaskPIDResponse(isMCForOtherSettings);
  }

  //=========  Set Cutnumber for V0Reader ================================
  TString cutnumberPhoton   = "00000008400100001500000000";
  if (periodNameV0Reader.CompareTo("LHC17n") == 0 || periodNameV0Reader.Contains("LHC17j7"))
    cutnumberPhoton         = "00000088400100001500000000";
  TString cutnumberEvent    = "10000003";

  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();

  //========= Add V0 Reader to  ANALYSIS manager if not yet existent =====
  TString V0ReaderName = Form("V0ReaderV1_%s_%s",cutnumberEvent.Data(),cutnumberPhoton.Data());
  if( !(AliV0ReaderV1*)mgr->GetTask(V0ReaderName.Data()) ){
    AliV0ReaderV1 *fV0ReaderV1 = new AliV0ReaderV1(V0ReaderName.Data());
    if (periodNameV0Reader.CompareTo("") != 0) fV0ReaderV1->SetPeriodName(periodNameV0Reader);
    fV0ReaderV1->SetUseOwnXYZCalculation(kTRUE);
    fV0ReaderV1->SetCreateAODs(kFALSE);// AOD Output
    fV0ReaderV1->SetUseAODConversionPhoton(kTRUE);

    if (!mgr) {
      Error("AddTask_V0ReaderV1", "No analysis manager found.");
      return;
    }

    AliConvEventCuts *fEventCuts=NULL;
    if(cutnumberEvent!=""){
      fEventCuts= new AliConvEventCuts(cutnumberEvent.Data(),cutnumberEvent.Data());
      fEventCuts->SetPreSelectionCutFlag(kTRUE);
      fEventCuts->SetV0ReaderName(V0ReaderName);
      if (periodNameV0Reader.CompareTo("") != 0) fEventCuts->SetPeriodEnum(periodNameV0Reader);
      fEventCuts->SetLightOutput(runLightOutput);
      if(fEventCuts->InitializeCutsFromCutString(cutnumberEvent.Data())){
        fV0ReaderV1->SetEventCuts(fEventCuts);
        fEventCuts->SetFillCutHistograms("",kTRUE);
      }
    }

    // Set AnalysisCut Number
    AliConversionPhotonCuts *fCuts=NULL;
    if(cutnumberPhoton!=""){
      fCuts= new AliConversionPhotonCuts(cutnumberPhoton.Data(),cutnumberPhoton.Data());
      fCuts->SetPreSelectionCutFlag(kTRUE);
      fCuts->SetIsHeavyIon(isHeavyIon);
      fCuts->SetV0ReaderName(V0ReaderName);
      fCuts->SetLightOutput(runLightOutput);
      if(fCuts->InitializeCutsFromCutString(cutnumberPhoton.Data())){
        fV0ReaderV1->SetConversionCuts(fCuts);
        fCuts->SetFillCutHistograms("",kTRUE);
      }
    }

    if(inputHandler->IsA()==AliAODInputHandler::Class()){
    // AOD mode
      fV0ReaderV1->AliV0ReaderV1::SetDeltaAODBranchName(Form("GammaConv_%s_gamma",cutnumberAODBranch.Data()));
    }
    fV0ReaderV1->Init();

    AliLog::SetGlobalLogLevel(AliLog::kFatal);

    //connect input V0Reader
    mgr->AddTask(fV0ReaderV1);
    mgr->ConnectInput(fV0ReaderV1,0,cinput);

  }

  //================================================
  //========= Add task to the ANALYSIS manager =====
  //================================================
  AliAnalysisTaskGammaTriggerQA *task=NULL;
  task= new AliAnalysisTaskGammaTriggerQA(Form("GammaTriggerQA_%i",trainConfig));
  task->SetIsHeavyIon(isHeavyIon);
  task->SetIsMC(isMC);
  task->SetV0ReaderName(V0ReaderName);
  task->SetCorrectionTaskSetting(corrTaskSetting);
  task->SetLightOutput(runLightOutput);

  //create cut handler
  CutHandlerTriggerQA cuts;

  // **********************************************************************************************************
  // ***************************** EMCAL configurations *******************************************************
  // **********************************************************************************************************
  if (trainConfig == 1){ // EMCAL clusters
    cuts.AddCut("30100013","1111100053032230000"); // 0-5%
    cuts.AddCut("150100013","1111100053032230000"); // 0-10%
    cuts.AddCut("12500013","1111100053032230000"); // 20-50%
  } else if (trainConfig == 2){ // EMCAL clusters
    cuts.AddCut("10910013","1111100053032230000"); // 0-90%
    cuts.AddCut("1093a013","1111100053032230000"); // 0-90%
    cuts.AddCut("1093b013","1111100053032230000"); // 0-90%
    cuts.AddCut("1093c013","1111100053032230000"); // 0-90%
    cuts.AddCut("1093d013","1111100053032230000"); // 0-90%
    cuts.AddCut("1093e013","1111100053032230000"); // 0-90%
  } else if (trainConfig == 3){ // EMCAL clusters
    cuts.AddCut("12510013","1111100053032230000"); // 20-50%
    cuts.AddCut("1253a013","1111100053032230000"); // 20-50%
    cuts.AddCut("1253b013","1111100053032230000"); // 20-50%
    cuts.AddCut("1253c013","1111100053032230000"); // 20-50%
    cuts.AddCut("1253d013","1111100053032230000"); // 20-50%
    cuts.AddCut("1253e013","1111100053032230000"); // 20-50%
  } else if (trainConfig == 4){ // EMCAL clusters
    cuts.AddCut("10110013","1111100053032230000"); // 0-10%
    cuts.AddCut("1013a013","1111100053032230000"); // 0-10%
    cuts.AddCut("1013b013","1111100053032230000"); // 0-10%
    cuts.AddCut("1013c013","1111100053032230000"); // 0-10%
    cuts.AddCut("1013d013","1111100053032230000"); // 0-10%
    cuts.AddCut("1013e013","1111100053032230000"); // 0-10%
  } else if (trainConfig == 5){ // EMCAL clusters
    cuts.AddCut("13510013","1111100053032230000"); // 30-50%
    cuts.AddCut("1353a013","1111100053032230000"); // 30-50%
    cuts.AddCut("1353b013","1111100053032230000"); // 30-50%
    cuts.AddCut("1353c013","1111100053032230000"); // 30-50%
    cuts.AddCut("1353d013","1111100053032230000"); // 30-50%
    cuts.AddCut("1353e013","1111100053032230000"); // 30-50%
  } else if (trainConfig == 6){ // EMCAL clusters
    cuts.AddCut("10910013","1111100053032230000"); // 0-90%
    // **********************************************************************************************************
  // ***************************** PHOS configurations ********************************************************
  // **********************************************************************************************************
  } else if (trainConfig == 101){ // PHOS clusters
    cuts.AddCut("30100013","2444400040033200000"); // 0-5%
    cuts.AddCut("10100013","2444400040033200000"); // 0-10%
    cuts.AddCut("12500013","2444400040033200000"); // 20-50%
  } else if (trainConfig == 107){ // PHOS clusters with TM
    cuts.AddCut("10900013","2444400002033200000"); // 0-90%


  } else {
    Error(Form("GammaConvCalo_%i",trainConfig), "wrong trainConfig variable no cuts have been specified for the configuration");
    return;
  }

  if(!cuts.AreValid()){
    cout << "\n\n****************************************************" << endl;
    cout << "ERROR: No valid cuts stored in CutHandlerTriggerQA! Returning..." << endl;
    cout << "****************************************************\n\n" << endl;
    return;
  }

  Int_t numberOfCuts = cuts.GetNCuts();

  TList *EventCutList   = new TList();
  TList *ClusterCutList = new TList();

  EventCutList->SetOwner(kTRUE);
  AliConvEventCuts **analysisEventCuts        = new AliConvEventCuts*[numberOfCuts];
  ClusterCutList->SetOwner(kTRUE);
  AliCaloPhotonCuts **analysisClusterCuts     = new AliCaloPhotonCuts*[numberOfCuts];

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

    analysisEventCuts[i]    = new AliConvEventCuts();
    analysisEventCuts[i]->SetV0ReaderName(V0ReaderName);
    if (periodNameV0Reader.CompareTo("") != 0) analysisEventCuts[i]->SetPeriodEnum(periodNameV0Reader);

    if(periodName.Contains("LHC11h") && doFlattening){
      cout << "entering the cent. flattening loop -> searching for file: " << fileNameInputForCentFlattening.Data() << endl;

      if( fileNameInputForCentFlattening.Contains("FlatFile") ){
        analysisEventCuts[i]->SetUseWeightFlatCentralityFromFile(doFlattening, fileNameInputForCentFlattening, "Cent");
      } else if( fileNameInputForCentFlattening.Contains("Good") ){
        analysisEventCuts[i]->SetUseWeightFlatCentralityFromFile(doFlattening, fileNameInputForCentFlattening, "CentGoodRuns");
      }else if( fileNameInputForCentFlattening.Contains("SemiGood") ){
        analysisEventCuts[i]->SetUseWeightFlatCentralityFromFile(doFlattening, fileNameInputForCentFlattening, "CentSemiGoodRuns");
      }else {
        analysisEventCuts[i]->SetUseWeightFlatCentralityFromFile(doFlattening, fileNameInputForCentFlattening, "CentTotalRuns");
      }
    }

    analysisEventCuts[i]->SetDebugLevel(0);
    analysisEventCuts[i]->SetLightOutput(runLightOutput);
    analysisEventCuts[i]->InitializeCutsFromCutString((cuts.GetEventCut(i)).Data());
    EventCutList->Add(analysisEventCuts[i]);
    analysisEventCuts[i]->SetCorrectionTaskSetting(corrTaskSetting);
    analysisEventCuts[i]->SetFillCutHistograms("",kTRUE);

    analysisClusterCuts[i]  = new AliCaloPhotonCuts(isMC);
    analysisClusterCuts[i]->SetV0ReaderName(V0ReaderName);
    analysisClusterCuts[i]->SetCorrectionTaskSetting(corrTaskSetting);
    analysisClusterCuts[i]->SetCaloTrackMatcherName(TrackMatcherName);
    analysisClusterCuts[i]->SetLightOutput(runLightOutput);
    analysisClusterCuts[i]->InitializeCutsFromCutString((cuts.GetClusterCut(i)).Data());
    ClusterCutList->Add(analysisClusterCuts[i]);
    analysisClusterCuts[i]->SetExtendedMatchAndQA(enableExtMatchAndQA);
    analysisClusterCuts[i]->SetFillCutHistograms("");
  }
  task->SetEventCutList(numberOfCuts,EventCutList);
  task->SetCaloCutList(numberOfCuts,ClusterCutList);
  task->SetCorrectionTaskSetting(corrTaskSetting);
  task->SetDetailedQAFlag(QAflag);

  //connect containers
  TString outputName                        = Form("GammaTriggerQA_%i",trainConfig);
  if (corrTaskSetting.CompareTo("") != 0)
    outputName                              = Form("GammaTriggerQA_%i_%s",trainConfig,corrTaskSetting.Data());

  AliAnalysisDataContainer *coutput =  mgr->CreateContainer(outputName.Data(), TList::Class(),
                                                            AliAnalysisManager::kOutputContainer,"GammaTriggerQA.root" );

  mgr->AddTask(task);
  mgr->ConnectInput(task,0,cinput);
  mgr->ConnectOutput(task,1,coutput);
  return;
}
