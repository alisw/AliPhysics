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
                                  TString   fileNameInputForWeighting   = "MCSpectraInput.root",       // path to file for weigting input
                                  TString   cutnumberAODBranch          = "000000006008400001001500000",
                                  TString   periodname                  = "LHC12f1x",         // period name
                                  Int_t     doWeightingPart             = 0,                  // enables weighting
                                  Int_t     enableExtMatchAndQA         = 0,                  // enable QA(3), extMatch+QA(2), extMatch(1), disabled (0)
                                  Bool_t    enableTriggerMimicking      = kFALSE,             // enable trigger mimicking
                                  Bool_t    enableTriggerOverlapRej     = kFALSE,             // enable trigger overlap rejection
                                  Float_t   maxFacPtHard                = 3.,                 // maximum factor between hardest jet and ptHard generated
                                  TString   periodNameV0Reader          = "",                 // set period name for V0 Reader
                                  Bool_t    enableSortingMCLabels       = kTRUE,              // enable sorting for MC cluster labels
                                  Bool_t    runLightOutput              = kFALSE              // switch to run light output (only essential histograms for afterburner)
) {
  
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

    if (!mgr) {
      Error("AddTask_V0ReaderV1", "No analysis manager found.");
      return;
    }
    AliConvEventCuts *fEventCuts    = NULL;
    if(cutnumberEvent!=""){
      fEventCuts                    = new AliConvEventCuts(cutnumberEvent.Data(),cutnumberEvent.Data());
      fEventCuts->SetPreSelectionCutFlag(kTRUE);
      fEventCuts->SetV0ReaderName(V0ReaderName);
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
      fV0ReaderV1->SetDeltaAODBranchName(Form("GammaConv_%s_gamma",cutnumberAODBranch.Data()));
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
  task->SetLightOutput(runLightOutput);

  //create cut handler
  CutHandlerCaloMerged cuts;

  // cluster cuts
  // 0 "ClusterType",  1 "EtaMin", 2 "EtaMax", 3 "PhiMin", 4 "PhiMax", 5 "DistanceToBadChannel", 6 "Timing", 7 "TrackMatching", 8 "ExoticCell",
  // 9 "MinEnergy", 10 "MinNCells", 11 "MinM02", 12 "MaxM02", 13 "MinM20", 14 "MaxM20", 15 "MaximumDispersion", 16 "NLM"
  
  // ************************************* EMCAL cuts ****************************************************
  // LHC13b-d
  if (trainConfig == 1){ // NLM 1 no non linearity
    cuts.AddCut("80000013","1111100050032200000","1111100050022110001","0163300000000000"); //
    cuts.AddCut("80052013","1111100050032200000","1111100050022110001","0163300000000000"); //
    cuts.AddCut("80085013","1111100050032200000","1111100050022110001","0163300000000000"); //
    cuts.AddCut("80083013","1111100050032200000","1111100050022110001","0163300000000000"); //
  } else if (trainConfig == 2){ // NLM 2 no non linearity
    cuts.AddCut("80000013","1111100050032200000","1111100050022110002","0163302200000000"); //
    cuts.AddCut("80052013","1111100050032200000","1111100050022110002","0163302200000000"); //
    cuts.AddCut("80085013","1111100050032200000","1111100050022110002","0163302200000000"); //
    cuts.AddCut("80083013","1111100050032200000","1111100050022110002","0163302200000000"); //
  } else if (trainConfig == 3){ // NLM 1 conv non linearity
    cuts.AddCut("80000013","1111141053032200000","1111141053022110001","0163301100000000"); //
    cuts.AddCut("80052013","1111141053032200000","1111141053022110001","0163301100000000"); //
    cuts.AddCut("80085013","1111141053032200000","1111141053022110001","0163301100000000"); //
    cuts.AddCut("80083013","1111141053032200000","1111141053022110001","0163301100000000"); //
  } else if (trainConfig == 4){ // NLM 2 conv non linearity
    cuts.AddCut("80000013","1111141053032200000","1111141053022110002","0163302200000000"); //
    cuts.AddCut("80052013","1111141053032200000","1111141053022110002","0163302200000000"); //
    cuts.AddCut("80085013","1111141053032200000","1111141053022110002","0163302200000000"); //
    cuts.AddCut("80083013","1111141053032200000","1111141053022110002","0163302200000000"); //
  } else if (trainConfig == 5){ // NLM 1 conv non linearity
    cuts.AddCut("80000013","1111141053032200000","1111141053022210001","0163301100000000"); //
    cuts.AddCut("80052013","1111141053032200000","1111141053022210001","0163301100000000"); //
    cuts.AddCut("80085013","1111141053032200000","1111141053022210001","0163301100000000"); //
    cuts.AddCut("80083013","1111141053032200000","1111141053022210001","0163301100000000"); //
  } else if (trainConfig == 6){ // NLM 2 conv non linearity
    cuts.AddCut("80000013","1111141053032200000","1111141053022210002","0163302200000000"); //
    cuts.AddCut("80052013","1111141053032200000","1111141053022210002","0163302200000000"); //
    cuts.AddCut("80085013","1111141053032200000","1111141053022210002","0163302200000000"); //
    cuts.AddCut("80083013","1111141053032200000","1111141053022210002","0163302200000000"); //
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
  
  EventCutList->SetOwner(kTRUE);
  AliConvEventCuts **analysisEventCuts          = new AliConvEventCuts*[numberOfCuts];
  ClusterCutList->SetOwner(kTRUE);
  AliCaloPhotonCuts **analysisClusterCuts       = new AliCaloPhotonCuts*[numberOfCuts];
  ClusterMergedCutList->SetOwner(kTRUE);
  AliCaloPhotonCuts **analysisClusterMergedCuts = new AliCaloPhotonCuts*[numberOfCuts];
  MesonCutList->SetOwner(kTRUE);
  AliConversionMesonCuts **analysisMesonCuts    = new AliConversionMesonCuts*[numberOfCuts];

  for(Int_t i = 0; i<numberOfCuts; i++){
    analysisEventCuts[i]          = new AliConvEventCuts();

    analysisEventCuts[i]->SetTriggerMimicking(enableTriggerMimicking);
    analysisEventCuts[i]->SetTriggerOverlapRejecion(enableTriggerOverlapRej);
    analysisEventCuts[i]->SetMaxFacPtHard(maxFacPtHard);
    analysisEventCuts[i]->SetV0ReaderName(V0ReaderName);
    analysisEventCuts[i]->SetLightOutput(runLightOutput);
    analysisEventCuts[i]->InitializeCutsFromCutString((cuts.GetEventCut(i)).Data());
    EventCutList->Add(analysisEventCuts[i]);
    analysisEventCuts[i]->SetFillCutHistograms("",kFALSE);
    
    analysisClusterCuts[i]        = new AliCaloPhotonCuts();
    analysisClusterCuts[i]->SetIsPureCaloCut(2);
    analysisClusterCuts[i]->SetV0ReaderName(V0ReaderName);
    analysisClusterCuts[i]->SetLightOutput(runLightOutput);
    analysisClusterCuts[i]->InitializeCutsFromCutString((cuts.GetClusterCut(i)).Data());
    ClusterCutList->Add(analysisClusterCuts[i]);
    analysisClusterCuts[i]->SetExtendedMatchAndQA(enableExtMatchAndQA);
    analysisClusterCuts[i]->SetFillCutHistograms("");
    
    analysisClusterMergedCuts[i]  = new AliCaloPhotonCuts();
    analysisClusterMergedCuts[i]->SetIsPureCaloCut(1);
    analysisClusterMergedCuts[i]->SetV0ReaderName(V0ReaderName);
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
  task->SetEventCutList(numberOfCuts,EventCutList);
  task->SetCaloCutList(numberOfCuts,ClusterCutList);
  task->SetCaloMergedCutList(numberOfCuts,ClusterMergedCutList);
  task->SetMesonCutList(numberOfCuts,MesonCutList);
  task->SetDoMesonQA(enableQAMesonTask); //Attention new switch for Pi0 QA
  task->SetDoClusterQA(enableQAClusterTask);  //Attention new switch small for Cluster QA
  task->SetEnableSortingOfMCClusLabels(enableSortingMCLabels);
  if(enableExtMatchAndQA == 2 || enableExtMatchAndQA == 3){ task->SetPlotHistsExtQA(kTRUE);}
  
  //connect containers
  AliAnalysisDataContainer *coutput =
    mgr->CreateContainer(Form("GammaCaloMerged_%i",trainConfig), TList::Class(),
              AliAnalysisManager::kOutputContainer,Form("GammaCaloMerged_%i.root",trainConfig));

  mgr->AddTask(task);
  mgr->ConnectInput(task,0,cinput);
  mgr->ConnectOutput(task,1,coutput);

  return;

}
