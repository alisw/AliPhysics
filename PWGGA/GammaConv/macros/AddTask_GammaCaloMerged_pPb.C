/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: Friederike Bock, Daniel Mühlheim                               *
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
//($ALIPHYSICS/PWGGA/GammaConv/AliAnalysisTaskGammaCaloMerged.cxx) for
//pPb together with all supporting classes
//***************************************************************************************

//***************************************************************************************
//CutHandler contains all cuts for a certain analysis and trainconfig,
//it automatically checks length of cutStrings and takes care of the number of added cuts,
//no specification of the variable 'numberOfCuts' needed anymore.
//***************************************************************************************
class CutHandlerCaloMerged{
  public:
    CutHandlerCaloMerged(Int_t nMax=10){
      nCuts=0; nMaxCuts=nMax; validCuts = true;
      eventCutArray = new TString[nMaxCuts]; clusterCutArray = new TString[nMaxCuts]; clusterMergedCutArray = new TString[nMaxCuts]; mesonCutArray = new TString[nMaxCuts];
      for(Int_t i=0; i<nMaxCuts; i++) {eventCutArray[i] = ""; clusterCutArray[i] = ""; clusterMergedCutArray[i] = ""; mesonCutArray[i] = "";}
    }

    void AddCut(TString eventCut, TString clusterCut, TString clusterMergedCut, TString mesonCut){
      if(nCuts>=nMaxCuts) {cout << "ERROR in CutHandlerCaloMerged: Exceeded maximum number of cuts!" << endl; validCuts = false; return;}
      if( eventCut.Length()!=8 || clusterCut.Length()!=19 || clusterMergedCut.Length()!=19 || mesonCut.Length()!=16 ) {cout << "ERROR in CutHandlerCaloMerged: Incorrect length of cut string!" << endl; validCuts = false; return;}
      eventCutArray[nCuts]=eventCut; clusterCutArray[nCuts]=clusterCut; clusterMergedCutArray[nCuts]=clusterMergedCut; mesonCutArray[nCuts]=mesonCut;
      nCuts++;
      return;
    }
    Bool_t AreValid(){return validCuts;}
    Int_t GetNCuts(){if(validCuts) return nCuts; else return 0;}
    TString GetEventCut(Int_t i){if(validCuts&&i<nMaxCuts&&i>=0) return eventCutArray[i]; else{cout << "ERROR in CutHandlerCaloMerged: GetEventCut wrong index i" << endl;return "";}}
    TString GetClusterCut(Int_t i){if(validCuts&&i<nMaxCuts&&i>=0) return clusterCutArray[i]; else {cout << "ERROR in CutHandlerCaloMerged: GetClusterCut wrong index i" << endl;return "";}}
    TString GetClusterMergedCut(Int_t i){if(validCuts&&i<nMaxCuts&&i>=0) return clusterMergedCutArray[i]; else {cout << "ERROR in CutHandlerCaloMerged: GetClusterMergedCut wrong index i" << endl;return "";}}
    TString GetMesonCut(Int_t i){if(validCuts&&i<nMaxCuts&&i>=0) return mesonCutArray[i]; else {cout << "ERROR in CutHandlerCaloMerged: GetMesonCut wrong index i" << endl;return "";}}
  private:
    Bool_t validCuts;
    Int_t nCuts; Int_t nMaxCuts;
    TString* eventCutArray;
    TString* clusterCutArray;
    TString* clusterMergedCutArray;
    TString* mesonCutArray;
};

//***************************************************************************************
//main function
//***************************************************************************************
void AddTask_GammaCaloMerged_pPb( Int_t     trainConfig                 = 1,                  // change different set of cuts
                                  Int_t     isMC                        = 0,                  // run MC
                                  Int_t     enableQAMesonTask           = 0,                  // enable QA in AliAnalysisTaskGammaCalo
                                  Int_t     enableQAClusterTask         = 0,                  // enable additional QA task
                                  TString   fileNameInputForWeighting   = "MCSpectraInput.root",       // path to file for weigting input / modified acceptance
                                  TString   cutnumberAODBranch          = "000000006008400001001500000",
                                  TString   periodname                  = "LHC12f1x",         // period name
                                  Int_t     doWeightingPart             = 0,                  // enables weighting
                                  Int_t     enableExtMatchAndQA         = 0,                  // disabled (0), extMatch (1), extQA_noCellQA (2), extMatch+extQA_noCellQA (3), extQA+cellQA (4), extMatch+extQA+cellQA (5)
                                  Bool_t    enableTriggerMimicking      = kFALSE,             // enable trigger mimicking
                                  Bool_t    enableTriggerOverlapRej     = kFALSE,             // enable trigger overlap rejection
                                  Float_t   maxFacPtHard                = 3.,                 // maximum factor between hardest jet and ptHard generated
                                  TString   periodNameV0Reader          = "",                 // set period name for V0 Reader
                                  Bool_t    enableSortingMCLabels       = kTRUE,              // enable sorting for MC cluster labels
                                  Bool_t    runLightOutput              = kFALSE,             // switch to run light output (only essential histograms for afterburner)
                                  Bool_t    runDetailedM02              = kFALSE,             // switch on very detailed M02 distribution
                                  TString   additionalTrainConfig       = "0"                 // additional counter for trainconfig, this has to be always the last parameter
) {

  Bool_t doTreeEOverP = kFALSE; // switch to produce EOverP tree
  TH1S* histoAcc = 0x0;         // histo for modified acceptance
  TString corrTaskSetting = ""; // sets which correction task setting to use
  //parse additionalTrainConfig flag
  TObjArray *rAddConfigArr = additionalTrainConfig.Tokenize("_");
  if(rAddConfigArr->GetEntries()<1){cout << "ERROR: AddTask_GammaCaloMerged_pPb during parsing of additionalTrainConfig String '" << additionalTrainConfig.Data() << "'" << endl; return;}
  TObjString* rAdditionalTrainConfig;
  for(Int_t i = 0; i<rAddConfigArr->GetEntries() ; i++){
    if(i==0) rAdditionalTrainConfig = (TObjString*)rAddConfigArr->At(i);
    else{
      TObjString* temp = (TObjString*) rAddConfigArr->At(i);
      TString tempStr = temp->GetString();
      if(tempStr.CompareTo("EPCLUSTree") == 0){
        cout << "INFO: AddTask_GammaCaloMerged_pPb activating 'EPCLUSTree'" << endl;
        doTreeEOverP = kTRUE;
      }else if(tempStr.BeginsWith("MODIFYACC")){
        cout << "INFO: AddTask_GammaCaloMerged_pPb activating 'MODIFYACC'" << endl;
        TString tempType = tempStr;
        tempType.Replace(0,9,"");
        cout << "INFO: connecting to alien..." << endl;
        TGrid::Connect("alien://");
        cout << "done!" << endl;
        TFile *w = TFile::Open(fileNameInputForWeighting.Data());
        if(!w){cout << "ERROR: Could not open file: " << fileNameInputForWeighting.Data() << endl;return;}
        histoAcc = (TH1S*) w->Get(tempType.Data());
        if(!histoAcc) {cout << "ERROR: Could not find histo: " << tempType.Data() << endl;return;}
        cout << "found: " << histoAcc << endl;
      }else if(tempStr.BeginsWith("CF")){
        cout << "INFO: AddTask_GammaCaloMerged_pPb will use custom branch from Correction Framework!" << endl;
        corrTaskSetting = tempStr;
        corrTaskSetting.Replace(0,2,"");
      }
    }
  }

  TString sAdditionalTrainConfig = rAdditionalTrainConfig->GetString();
  if (sAdditionalTrainConfig.Atoi() > 0){
    trainConfig = trainConfig + sAdditionalTrainConfig.Atoi();
    cout << "INFO: AddTask_GammaCaloMerged_pPb running additionalTrainConfig '" << sAdditionalTrainConfig.Atoi() << "', train config: '" << trainConfig << "'" << endl;
  }
  cout << "corrTaskSetting: " << corrTaskSetting.Data() << endl;

  Int_t isHeavyIon = 2;

  // ================== GetAnalysisManager ===============================
  AliAnalysisManager *mgr           = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error(Form("AddTask_GammaCaloMerged_pp_%i",trainConfig), "No analysis manager found.");
    return ;
  }

  // ================== GetInputEventHandler =============================
  AliVEventHandler *inputHandler    = mgr->GetInputEventHandler();
  Bool_t isMCForOtherSettings       = 0;
  if (isMC > 0)
    isMCForOtherSettings            = 1;
  //========= Add PID Reponse to ANALYSIS manager ====
  if(!(AliPIDResponse*)mgr->GetTask("PIDResponseTask")){
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
    AddTaskPIDResponse(isMCForOtherSettings);
  }

  Printf("here \n");

  //=========  Set Cutnumber for V0Reader ================================
  TString cutnumberPhoton           = "00000008400100001500000000";
  TString cutnumberEvent            = "80000003";
  Bool_t doEtaShift                 = kFALSE;
  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();

  //========= Add V0 Reader to  ANALYSIS manager if not yet existent =====
  TString V0ReaderName = Form("V0ReaderV1_%s_%s",cutnumberEvent.Data(),cutnumberPhoton.Data());
  if( !(AliV0ReaderV1*)mgr->GetTask(V0ReaderName.Data()) ){
    AliV0ReaderV1 *fV0ReaderV1      = new AliV0ReaderV1(V0ReaderName.Data());
    if (periodNameV0Reader.CompareTo("") != 0) fV0ReaderV1->SetPeriodName(periodNameV0Reader);
    fV0ReaderV1->SetUseOwnXYZCalculation(kTRUE);
    fV0ReaderV1->SetCreateAODs(kFALSE);// AOD Output
    fV0ReaderV1->SetUseAODConversionPhoton(kTRUE);

    AliConvEventCuts *fEventCuts    = NULL;
    if(cutnumberEvent!=""){
      fEventCuts                    = new AliConvEventCuts(cutnumberEvent.Data(),cutnumberEvent.Data());
      fEventCuts->SetPreSelectionCutFlag(kTRUE);
      fEventCuts->SetV0ReaderName(V0ReaderName);
      if(periodNameV0Reader.CompareTo("") != 0) fEventCuts->SetPeriodEnum(periodNameV0Reader);
      fEventCuts->SetLightOutput(runLightOutput);
      if(fEventCuts->InitializeCutsFromCutString(cutnumberEvent.Data())){
        fEventCuts->DoEtaShift(doEtaShift);
        fV0ReaderV1->SetEventCuts(fEventCuts);
        fEventCuts->SetFillCutHistograms("",kTRUE);
      }
    }
    // Set AnalysisCut Number
    AliConversionPhotonCuts *fCuts  = NULL;
    if(cutnumberPhoton!=""){
      fCuts                         = new AliConversionPhotonCuts(cutnumberPhoton.Data(),cutnumberPhoton.Data());
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
  AliAnalysisTaskGammaCaloMerged *task  = NULL;
  task                                  = new AliAnalysisTaskGammaCaloMerged(Form("GammaCaloMerged_%i",trainConfig));
  task->SetIsHeavyIon(isHeavyIon);
  task->SetIsMC(isMC);
  task->SetV0ReaderName(V0ReaderName);
  task->SetCorrectionTaskSetting(corrTaskSetting);
  task->SetLightOutput(runLightOutput);

  //create cut handler
  CutHandlerCaloMerged cuts;

  // cluster cuts
  // 0 "ClusterType",  1 "EtaMin", 2 "EtaMax", 3 "PhiMin", 4 "PhiMax", 5 "DistanceToBadChannel", 6 "Timing", 7 "TrackMatching", 8 "ExoticCell",
  // 9 "MinEnergy", 10 "MinNCells", 11 "MinM02", 12 "MaxM02", 13 "MinM20", 14 "MaxM20", 15 "MaximumDispersion", 16 "NLM"

  // ************************************* EMCAL cuts ****************************************************
  // LHC13b-d
  if (trainConfig == 1){ // NLM 1 no non linearity
    cuts.AddCut("80010013","1111100050032200000","1111100050022110001","0163300000000000"); //
    cuts.AddCut("80052013","1111100050032200000","1111100050022110001","0163300000000000"); //
    cuts.AddCut("80085013","1111100050032200000","1111100050022110001","0163300000000000"); //
    cuts.AddCut("80083013","1111100050032200000","1111100050022110001","0163300000000000"); //
  } else if (trainConfig == 2){ // NLM 2 no non linearity
    cuts.AddCut("80010013","1111100050032200000","1111100050022110002","0163302200000000"); //
    cuts.AddCut("80052013","1111100050032200000","1111100050022110002","0163302200000000"); //
    cuts.AddCut("80085013","1111100050032200000","1111100050022110002","0163302200000000"); //
    cuts.AddCut("80083013","1111100050032200000","1111100050022110002","0163302200000000"); //
  } else if (trainConfig == 3){ // NLM 1 conv non linearity
    cuts.AddCut("80010013","1111141053032200000","1111141053022110001","0163301100000000"); //
    cuts.AddCut("80052013","1111141053032200000","1111141053022110001","0163301100000000"); //
    cuts.AddCut("80085013","1111141053032200000","1111141053022110001","0163301100000000"); //
    cuts.AddCut("80083013","1111141053032200000","1111141053022110001","0163301100000000"); //
  } else if (trainConfig == 4){ // NLM 2 conv non linearity
    cuts.AddCut("80010013","1111141053032200000","1111141053022110002","0163302200000000"); //
    cuts.AddCut("80052013","1111141053032200000","1111141053022110002","0163302200000000"); //
    cuts.AddCut("80085013","1111141053032200000","1111141053022110002","0163302200000000"); //
    cuts.AddCut("80083013","1111141053032200000","1111141053022110002","0163302200000000"); //
  } else if (trainConfig == 5){ // NLM 1 conv non linearity
    cuts.AddCut("80010013","1111141053032200000","1111141053022210001","0163301100000000"); //
    cuts.AddCut("80052013","1111141053032200000","1111141053022210001","0163301100000000"); //
    cuts.AddCut("80085013","1111141053032200000","1111141053022210001","0163301100000000"); //
    cuts.AddCut("80083013","1111141053032200000","1111141053022210001","0163301100000000"); //
  } else if (trainConfig == 6){ // NLM 2 conv non linearity
    cuts.AddCut("80010013","1111141053032200000","1111141053022210002","0163302200000000"); //
    cuts.AddCut("80052013","1111141053032200000","1111141053022210002","0163302200000000"); //
    cuts.AddCut("80085013","1111141053032200000","1111141053022210002","0163302200000000"); //
    cuts.AddCut("80083013","1111141053032200000","1111141053022210002","0163302200000000"); //

  // run 2 data
  } else if (trainConfig == 101){ // pp 2.76TeV paper cuts : open timing, TB nonlin
    cuts.AddCut("80010113","1111101017032200000","1111101017022700001","0163300000000000"); // INT7
    cuts.AddCut("80085113","1111101017032200000","1111101017022700001","0163300000000000"); // EG2
    cuts.AddCut("80083113","1111101017032200000","1111101017022700001","0163300000000000"); // EG1
  } else if (trainConfig == 102){  // pp 2.76TeV paper cuts:  w/o mass, open timing, TB nonlin
    cuts.AddCut("80010113","1111101017032200000","1111101017022000001","0163300000000000");
    cuts.AddCut("80052113","1111101017032200000","1111101017022000001","0163300000000000"); // EMC7
    cuts.AddCut("80085113","1111101017032200000","1111101017022000001","0163300000000000"); // EG2
    cuts.AddCut("80083113","1111101017032200000","1111101017022000001","0163300000000000"); // EG1
  } else if (trainConfig == 103){ // pp 2.76TeV paper cuts : open timing, TB nonlin
    cuts.AddCut("80010113","1111100017032200000","1111100017022700001","0163300000000000"); // INT7
    cuts.AddCut("80052113","1111100017032200000","1111100017022700001","0163300000000000"); // EMC7
    cuts.AddCut("80085113","1111100017032200000","1111100017022700001","0163300000000000"); // EG2
    cuts.AddCut("80083113","1111100017032200000","1111100017022700001","0163300000000000"); // EG1
  } else if (trainConfig == 104){  // pp 2.76TeV paper cuts:  w/o mass, open timing, TB nonlin
    cuts.AddCut("80010113","1111100017032200000","1111100017022000001","0163300000000000");
    cuts.AddCut("80052113","1111100017032200000","1111100017022000001","0163300000000000"); // EMC7
    cuts.AddCut("80085113","1111100017032200000","1111100017022000001","0163300000000000"); // EG2
    cuts.AddCut("80083113","1111100017032200000","1111100017022000001","0163300000000000"); // EG1
  } else if (trainConfig == 105){ // pp 2.76TeV paper cuts : open timing, TB nonlin, no TM
    cuts.AddCut("80010113","1111101010032200000","1111101010022700001","0163300000000000"); // INT7
    cuts.AddCut("80085113","1111101010032200000","1111101010022700001","0163300000000000"); // EG2
    cuts.AddCut("80083113","1111101010032200000","1111101010022700001","0163300000000000"); // EG1
  } else if (trainConfig == 106){ // pp 2.76TeV paper cuts : open timing, TB nonlin, no TM
    cuts.AddCut("80010103","1111101010032200000","1111101010022700001","0163300000000000"); // INT7
    cuts.AddCut("80085103","1111101010032200000","1111101010022700001","0163300000000000"); // EG2
    cuts.AddCut("80083103","1111101010032200000","1111101010022700001","0163300000000000"); // EG1
  } else if (trainConfig == 107){ // pp 2.76TeV paper cuts : open timing, TB nonlin, no TM
    cuts.AddCut("80010123","1111101010032200000","1111101010022700001","0163300000000000"); // INT7
    cuts.AddCut("80085123","1111101010032200000","1111101010022700001","0163300000000000"); // EG2
    cuts.AddCut("80083123","1111101010032200000","1111101010022700001","0163300000000000"); // EG1
  } else if (trainConfig == 108){ // pp 2.76TeV paper cuts : open timing, TB nonlin, no TM
    cuts.AddCut("80010103","1111101010032200000","1111101010022700001","0163300000000000"); // INT7
    cuts.AddCut("80010113","1111101010032200000","1111101010022700001","0163300000000000"); // INT7
    cuts.AddCut("80010123","1111101010032200000","1111101010022700001","0163300000000000"); // INT7

    // pPb 8 TeV LHC16 periods
  } else if (trainConfig == 200){ // open timing, TB nonlin, no TM
    cuts.AddCut("80010113","1111101010032200000","1111101010022700001","0163300000000000"); // INT7
  } else if (trainConfig == 201){ // open timing, TB nonlin, no TM
    cuts.AddCut("80052113","1111101010032200000","1111101010022700001","0163300000000000"); // EMC7
  } else if (trainConfig == 202){ // open timing, TB nonlin, no TM
    cuts.AddCut("80085113","1111101010032200000","1111101010022700001","0163300000000000"); // EG2
  } else if (trainConfig == 203){ // open timing, TB nonlin, no TM
    cuts.AddCut("80083113","1111101010032200000","1111101010022700001","0163300000000000"); // EG1
  } else if (trainConfig == 204){ // open timing, TB nonlin
    cuts.AddCut("80010113","1111101017032200000","1111101017022700001","0163300000000000"); // INT7
  } else if (trainConfig == 205){ // open timing, TB nonlin
    cuts.AddCut("80052113","1111101017032200000","1111101017022700001","0163300000000000"); // EMC7
  } else if (trainConfig == 206){ // open timing, TB nonlin
    cuts.AddCut("80085113","1111101017032200000","1111101017022700001","0163300000000000"); // EG2
  } else if (trainConfig == 207){ // open timing, TB nonlin
    cuts.AddCut("80083113","1111101017032200000","1111101017022700001","0163300000000000"); // EG1
  } else {
    Error(Form("GammaCaloMerged_%i",trainConfig), "wrong trainConfig variable no cuts have been specified for the configuration");
    return;
  }

  if(!cuts.AreValid()){
    cout << "\n\n****************************************************" << endl;
    cout << "ERROR: No valid cuts stored in CutHandlerCaloMerged! Returning..." << endl;
    cout << "****************************************************\n\n" << endl;
    return;
  }

  Int_t numberOfCuts            = cuts.GetNCuts();
  TList *EventCutList           = new TList();
  TList *ClusterCutList         = new TList();
  TList *ClusterMergedCutList   = new TList();
  TList *MesonCutList           = new TList();

  TList *HeaderList             = new TList();
  if (doWeightingPart==1) {
    TObjString *Header1         = new TObjString("pi0_1");
    HeaderList->Add(Header1);
  }
  if (doWeightingPart==2){
    TObjString *Header3         = new TObjString("eta_2");
    HeaderList->Add(Header3);
  }
  if (doWeightingPart==3) {
    TObjString *Header1         = new TObjString("pi0_1");
    HeaderList->Add(Header1);
    TObjString *Header3         = new TObjString("eta_2");
    HeaderList->Add(Header3);
  }

  if (periodNameV0Reader.Contains("LHC18b9")||periodNameV0Reader.Contains("LHC17g8")){
    TObjString *HeaderPMB = new TObjString("EPOSLHC_0");
    TObjString *HeaderP8J = new TObjString("Pythia8Jets_1");
    if (doWeightingPart==4) {
      HeaderList->Add(HeaderPMB);
      HeaderList->Add(HeaderP8J);
    } else {
      HeaderList->Add(HeaderP8J);
    }
  }
  EventCutList->SetOwner(kTRUE);
  AliConvEventCuts **analysisEventCuts          = new AliConvEventCuts*[numberOfCuts];
  ClusterCutList->SetOwner(kTRUE);
  AliCaloPhotonCuts **analysisClusterCuts       = new AliCaloPhotonCuts*[numberOfCuts];
  ClusterMergedCutList->SetOwner(kTRUE);
  AliCaloPhotonCuts **analysisClusterMergedCuts = new AliCaloPhotonCuts*[numberOfCuts];
  MesonCutList->SetOwner(kTRUE);
  AliConversionMesonCuts **analysisMesonCuts    = new AliConversionMesonCuts*[numberOfCuts];

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

    analysisEventCuts[i]          = new AliConvEventCuts();

    analysisEventCuts[i]->SetTriggerMimicking(enableTriggerMimicking);
    analysisEventCuts[i]->SetTriggerOverlapRejecion(enableTriggerOverlapRej);
    analysisEventCuts[i]->SetMaxFacPtHard(maxFacPtHard);
    analysisEventCuts[i]->SetV0ReaderName(V0ReaderName);
    analysisEventCuts[i]->SetCorrectionTaskSetting(corrTaskSetting);
    if(periodNameV0Reader.CompareTo("") != 0) analysisEventCuts[i]->SetPeriodEnum(periodNameV0Reader);
    analysisEventCuts[i]->SetLightOutput(runLightOutput);
    analysisEventCuts[i]->InitializeCutsFromCutString((cuts.GetEventCut(i)).Data());
    EventCutList->Add(analysisEventCuts[i]);
    analysisEventCuts[i]->SetFillCutHistograms("",kFALSE);

    analysisClusterCuts[i]        = new AliCaloPhotonCuts(isMC);
    analysisClusterCuts[i]->SetIsPureCaloCut(2);
    analysisClusterCuts[i]->SetHistoToModifyAcceptance(histoAcc);
    analysisClusterCuts[i]->SetV0ReaderName(V0ReaderName);
    analysisClusterCuts[i]->SetCorrectionTaskSetting(corrTaskSetting);
    analysisClusterCuts[i]->SetCaloTrackMatcherName(TrackMatcherName);
    analysisClusterCuts[i]->SetLightOutput(runLightOutput);
    analysisClusterCuts[i]->InitializeCutsFromCutString((cuts.GetClusterCut(i)).Data());
    ClusterCutList->Add(analysisClusterCuts[i]);
    analysisClusterCuts[i]->SetExtendedMatchAndQA(enableExtMatchAndQA);
    analysisClusterCuts[i]->SetFillCutHistograms("");

    analysisClusterMergedCuts[i]  = new AliCaloPhotonCuts(isMC);
    analysisClusterMergedCuts[i]->SetIsPureCaloCut(1);
    analysisClusterMergedCuts[i]->SetHistoToModifyAcceptance(histoAcc);
    analysisClusterMergedCuts[i]->SetV0ReaderName(V0ReaderName);
    analysisClusterMergedCuts[i]->SetCorrectionTaskSetting(corrTaskSetting);
    analysisClusterMergedCuts[i]->SetCaloTrackMatcherName(TrackMatcherName);
    analysisClusterMergedCuts[i]->SetLightOutput(runLightOutput);
    analysisClusterMergedCuts[i]->InitializeCutsFromCutString((cuts.GetClusterMergedCut(i)).Data());
    ClusterMergedCutList->Add(analysisClusterMergedCuts[i]);
    analysisClusterMergedCuts[i]->SetExtendedMatchAndQA(enableExtMatchAndQA);
    analysisClusterMergedCuts[i]->SetFillCutHistograms("");

    analysisMesonCuts[i]          = new AliConversionMesonCuts();
    analysisMesonCuts[i]->SetEnableOpeningAngleCut(kFALSE);
    analysisMesonCuts[i]->SetIsMergedClusterCut(1);
    analysisMesonCuts[i]->SetLightOutput(runLightOutput);
    analysisMesonCuts[i]->InitializeCutsFromCutString((cuts.GetMesonCut(i)).Data());
    MesonCutList->Add(analysisMesonCuts[i]);
    analysisMesonCuts[i]->SetFillCutHistograms("");
    analysisEventCuts[i]->SetAcceptedHeader(HeaderList);
  }
  task->SetEnableDetailedM02Distribtuon(runDetailedM02);
  task->SetEventCutList(numberOfCuts,EventCutList);
  task->SetCaloCutList(numberOfCuts,ClusterCutList);
  task->SetCorrectionTaskSetting(corrTaskSetting);
  task->SetCaloMergedCutList(numberOfCuts,ClusterMergedCutList);
  task->SetMesonCutList(numberOfCuts,MesonCutList);
  task->SetDoMesonQA(enableQAMesonTask); //Attention new switch for Pi0 QA
  task->SetDoClusterQA(enableQAClusterTask);  //Attention new switch small for Cluster QA
  task->SetEnableSortingOfMCClusLabels(enableSortingMCLabels);
  if(enableExtMatchAndQA > 1){ task->SetPlotHistsExtQA(kTRUE);}

  //connect containers
  AliAnalysisDataContainer *coutput =
    mgr->CreateContainer(!(corrTaskSetting.CompareTo("")) ? Form("GammaCaloMerged_%i",trainConfig) : Form("GammaCaloMerged_%i_%s",trainConfig,corrTaskSetting.Data()), TList::Class(),
              AliAnalysisManager::kOutputContainer,Form("GammaCaloMerged_%i.root",trainConfig) );

  mgr->AddTask(task);
  mgr->ConnectInput(task,0,cinput);
  mgr->ConnectOutput(task,1,coutput);

  return;

}
