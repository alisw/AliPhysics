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
//pp together with all supporting classes
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
void AddTask_GammaCaloMerged_pp(  Int_t     trainConfig                 = 1,                  // change different set of cuts
                                  Int_t     isMC                        = 0,                 // run MC
                                  Int_t     enableQAMesonTask           = 0,                 // enable QA in AliAnalysisTaskGammaCalo
                                  Int_t     enableQAClusterTask         = 0,                 // enable additional QA task
                                  TString   fileNameInputForWeighting   = "MCSpectraInput.root",       // path to file for weigting input
                                  TString   cutnumberAODBranch          = "000000006008400001001500000",
                                  TString   periodname                  = "LHC12f1x",             // period name
                                  Bool_t    doWeighting                 = kFALSE,              // enables weighting
                                  Int_t     enableExtMatchAndQA         = 0,                // enable QA(3), extMatch+QA(2), extMatch(1), disabled (0)
                                  Bool_t    enableTriggerMimicking      = kFALSE,              // enable trigger mimicking
                                  Bool_t    enableTriggerOverlapRej     = kFALSE,              // enable trigger overlap rejection
                                  Float_t   maxFacPtHard                = 3.,                // maximum factor between hardest jet and ptHard generated
                                  TString   periodNameV0Reader          = "",                // period Name for respective period selected in V0Reader
                                  Int_t     selectedMeson               = 1 
) {
  
  Int_t isHeavyIon = 0;
  
  // ================== GetAnalysisManager ===============================
  AliAnalysisManager *mgr           = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error(Form("AddTask_GammaCaloMerged_pp_%i",trainConfig), "No analysis manager found.");
    return ;
  }

  // ================== GetInputEventHandler =============================
  AliVEventHandler *inputHandler    = mgr->GetInputEventHandler();

  Bool_t isMCForOtherSettings = 0;
  if (isMC > 0) isMCForOtherSettings = 1;
  //========= Add PID Reponse to ANALYSIS manager ====
  if(!(AliPIDResponse*)mgr->GetTask("PIDResponseTask")){
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
    AddTaskPIDResponse(isMCForOtherSettings);
  }
  
  Printf("here \n");
  
  //=========  Set Cutnumber for V0Reader ================================
  TString cutnumberPhoton           = "00000008400100001500000000";
  TString cutnumberEvent            = "00000003";
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

  //create cut handler
  CutHandlerCaloMerged cuts;

  // cluster cuts
  // 0 "ClusterType",  1 "EtaMin", 2 "EtaMax", 3 "PhiMin", 4 "PhiMax", 5 "DistanceToBadChannel", 6 "Timing", 7 "TrackMatching", 8 "ExoticCell",
  // 9 "MinEnergy", 10 "MinNCells", 11 "MinM02", 12 "MaxM02", 13 "MinM20", 14 "MaxM20", 15 "MaximumDispersion", 16 "NLM"
  
  // ************************************* EMCAL cuts ****************************************************
  // LHC11a
  if (trainConfig == 1){ // all defaults for LHC11a
    cuts.AddCut("00003113","1111121053032200000","1111121053022210001","0163301100000000"); // INT1
    cuts.AddCut("00003113","1111121053032200000","1111121053022210002","0163302200000000"); // INT1
    cuts.AddCut("00051013","1111121053032200000","1111121053022210001","0163301100000000"); // EMC1
    cuts.AddCut("00051013","1111121053032200000","1111121053022210002","0163302200000000"); // EMC1
  } else if (trainConfig == 2){ // NLM 1 M02 var  
    cuts.AddCut("00003113","1111121053032200000","1111121053022110001","0163301100000000"); // min 0.3 function default
    cuts.AddCut("00003113","1111121053032200000","1111121053022310001","0163301100000000"); // min 0.25 function default
    cuts.AddCut("00003113","1111121053032200000","1111121053022410001","0163301100000000"); // min 0.27, tighter lower func
    cuts.AddCut("00003113","1111121053032200000","1111121053022520001","0163301100000000"); // min 0.27, looser func up and low
    cuts.AddCut("00003113","1111121053032200000","1111121053022430001","0163301100000000"); // min 0.27, tighter func up and low
  } else if (trainConfig == 3){ // NLM 2 M02 var
    cuts.AddCut("00003113","1111121053032200000","1111121053022110002","0163302200000000"); // min 0.3 function default
    cuts.AddCut("00003113","1111121053032200000","1111121053022310002","0163302200000000"); // min 0.25 function default
    cuts.AddCut("00003113","1111121053032200000","1111121053022410002","0163302200000000"); // min 0.27, tighter lower func
    cuts.AddCut("00003113","1111121053032200000","1111121053022520002","0163302200000000"); // min 0.27, looser func up and low
    cuts.AddCut("00003113","1111121053032200000","1111121053022430002","0163302200000000"); // min 0.27, tighter func up and low
  } else if (trainConfig == 4){ // NLM 1 M02 var  
    cuts.AddCut("00051013","1111121053032200000","1111121053022110001","0163301100000000"); // min 0.3 function default
    cuts.AddCut("00051013","1111121053032200000","1111121053022310001","0163301100000000"); // min 0.25 function default
    cuts.AddCut("00051013","1111121053032200000","1111121053022410001","0163301100000000"); // min 0.27, tighter lower func
    cuts.AddCut("00051013","1111121053032200000","1111121053022520001","0163301100000000"); // min 0.27, looser func up and low
    cuts.AddCut("00051013","1111121053032200000","1111121053022430001","0163301100000000"); // min 0.27, tighter func up and low
  } else if (trainConfig == 5){ // NLM 2 M02 var
    cuts.AddCut("00051013","1111121053032200000","1111121053022110002","0163302200000000"); // min 0.3 function default
    cuts.AddCut("00051013","1111121053032200000","1111121053022310002","0163302200000000"); // min 0.25 function default
    cuts.AddCut("00051013","1111121053032200000","1111121053022410002","0163302200000000"); // min 0.27, tighter lower func
    cuts.AddCut("00051013","1111121053032200000","1111121053022520002","0163302200000000"); // min 0.27, looser func up and low
    cuts.AddCut("00051013","1111121053032200000","1111121053022430002","0163302200000000"); // min 0.27, tighter func up and low
  } else if (trainConfig == 6){  // EMCAL clusters, variation track matching to cluster & mass variations INT1 NLM1 
    cuts.AddCut("00003113","1111121051032200000","1111121051022210001","0163301100000000"); // looser TM
    cuts.AddCut("00003113","1111121056032200000","1111121056022210001","0163301100000000"); // tighter TM
    cuts.AddCut("00003113","1111121053032200000","1111121053022210001","0163301300000000"); // tighter mass cut
    cuts.AddCut("00003113","1111121053032200000","1111121053022210001","0163301500000000"); // looser mass cut
  } else if (trainConfig == 7){  // EMCAL clusters, variation track matching to cluster & mass variations INT1 NLM2 
    cuts.AddCut("00003113","1111121051032200000","1111121051022210002","0163302200000000"); // looser TM
    cuts.AddCut("00003113","1111121056032200000","1111121056022210002","0163302200000000"); // tighter TM
    cuts.AddCut("00003113","1111121053032200000","1111121053022210002","0163302400000000"); // tighter mass cut
    cuts.AddCut("00003113","1111121053032200000","1111121053022210002","0163302600000000"); // looser mass cut
  } else if (trainConfig == 8){ // EMCAL clusters, variation track matching to cluster & mass variations EMC1 NLM2 
    cuts.AddCut("00051013","1111121051032200000","1111121051022210001","0163301100000000"); // looser TM
    cuts.AddCut("00051013","1111121056032200000","1111121056022210001","0163301100000000"); // tighter TM
    cuts.AddCut("00051013","1111121053032200000","1111121053022210001","0163301300000000"); // tighter mass cut
    cuts.AddCut("00051013","1111121053032200000","1111121053022210001","0163301500000000"); // looser mass cut    
  } else if (trainConfig == 9){ // EMCAL clusters, variation track matching to cluster & mass variations EMC1 NLM2 
    cuts.AddCut("00051013","1111121051032200000","1111121051022210002","0163302200000000"); // looser TM
    cuts.AddCut("00051013","1111121056032200000","1111121056022210002","0163302200000000"); // tighter TM
    cuts.AddCut("00051013","1111121053032200000","1111121053022210002","0163302400000000"); // tighter mass cut
    cuts.AddCut("00051013","1111121053032200000","1111121053022210002","0163302600000000"); // looser mass cut
  } else if (trainConfig == 10){ // NL var
    cuts.AddCut("00003113","1111122053032200000","1111122053022210001","0163301100000000"); // SDM loose time
    cuts.AddCut("00003113","1111111053032200000","1111111053022210001","0163301100000000"); // conv calo tight time
    cuts.AddCut("00003113","1111101053032200000","1111101053022210001","0163301100000000"); // SDM Jason
    cuts.AddCut("00003113","1111100053032200000","1111100053022210001","0163301100000000"); // none
  } else if (trainConfig == 11){ // NL var
    cuts.AddCut("00003113","1111122053032200000","1111122053022210002","0163302200000000"); // SDM loose time
    cuts.AddCut("00003113","1111111053032200000","1111111053022210002","0163302200000000"); // conv calo tight time
    cuts.AddCut("00003113","1111101053032200000","1111101053022210002","0163302200000000"); // SDM Jason
    cuts.AddCut("00003113","1111100053032200000","1111100053022210002","0163302200000000"); // none
  } else if (trainConfig == 12){ // NL var
    cuts.AddCut("00051013","1111122053032200000","1111122053022210001","0163301100000000"); // SDM loose time
    cuts.AddCut("00051013","1111111053032200000","1111111053022210001","0163301100000000"); // conv calo tight time
    cuts.AddCut("00051013","1111101053032200000","1111101053022210001","0163301100000000"); // SDM Jason
    cuts.AddCut("00051013","1111100053032200000","1111100053022210001","0163301100000000"); // none
  } else if (trainConfig == 13){ // NL var
    cuts.AddCut("00051013","1111122053032200000","1111122053022210002","0163302200000000"); // SDM loose time
    cuts.AddCut("00051013","1111111053032200000","1111111053022210002","0163302200000000"); // conv calo tight time
    cuts.AddCut("00051013","1111101053032200000","1111101053022210002","0163302200000000"); // SDM Jason
    cuts.AddCut("00051013","1111100053032200000","1111100053022210002","0163302200000000"); // none
  } else if (trainConfig == 14){ // Alpha cut variations & TRD material INT1, NLM1
    cuts.AddCut("00003113","1111121053032200000","1111121053022210001","0163303100000000"); // NLM 1 looser
    cuts.AddCut("00003113","1111121053032200000","1111121053022210001","0163305100000000"); // NLM 1 tighter
    cuts.AddCut("00003113","1113121053032200000","1113121053022210001","0163301100000000");// TRD infront
    cuts.AddCut("00003113","1111221053032200000","1111221053022210001","0163301100000000");// no TRD infront
  } else if (trainConfig == 15){ // Alpha cut variations & TRD material INT1, NLM2
    cuts.AddCut("00003113","1111121053032200000","1111121053022210002","0163304200000000"); // NLM 2 looser
    cuts.AddCut("00003113","1111121053032200000","1111121053022210002","0163306200000000"); // NLM 2 tighter
    cuts.AddCut("00003113","1113121053032200000","1113121053022210002","0163302200000000");// TRD infront
    cuts.AddCut("00003113","1111221053032200000","1111221053022210002","0163302200000000");// no TRD infront
  } else if (trainConfig == 16){ // Alpha cut variations & TRD material EMC1, NLM1
    cuts.AddCut("00051013","1111121053032200000","1111121053022210001","0163303100000000"); // NLM 1 looser
    cuts.AddCut("00051013","1111121053032200000","1111121053022210001","0163305100000000"); // NLM 1 tighter
    cuts.AddCut("00051013","1113121053032200000","1113121053022210001","0163301100000000");// TRD infront
    cuts.AddCut("00051013","1111221053032200000","1111221053022210001","0163301100000000");// no TRD infront
  } else if (trainConfig == 17){ // Alpha cut variations & TRD material EMC1, NLM2
    cuts.AddCut("00051013","1113121053032200000","1113121053022210002","0163302200000000");// TRD infront
    cuts.AddCut("00051013","1111221053032200000","1111221053022210002","0163302200000000");// no TRD infront
    cuts.AddCut("00051013","1111121053032200000","1111121053022210002","0163304200000000"); // NLM 2 looser
    cuts.AddCut("00051013","1111121053032200000","1111121053022210002","0163306200000000"); // NLM 2 tighter
  } else if (trainConfig == 18){ // no TM in basis cut
    cuts.AddCut("00003113","1111121050032200000","1111121053022210001","0163301100000000"); // INT1
    cuts.AddCut("00003113","1111121050032200000","1111121053022210002","0163302200000000"); // INT1
    cuts.AddCut("00051013","1111121050032200000","1111121053022210001","0163301100000000"); // EMC1
    cuts.AddCut("00051013","1111121050032200000","1111121053022210002","0163302200000000"); // EMC1
  } else if (trainConfig == 19){ // open M02, open Mass, open Alpha
    cuts.AddCut("00003113","1111121053032200000","1111121053022000001","0163300000000000"); // INT1
    cuts.AddCut("00003113","1111121053032200000","1111121053022000002","0163300000000000"); // INT1
    cuts.AddCut("00051013","1111121053032200000","1111121053022000001","0163300000000000"); // EMC1
    cuts.AddCut("00051013","1111121053032200000","1111121053022000002","0163300000000000"); // EMC1
    
  // LHC13g
  } else if (trainConfig == 40){  // new defaults LHC13g NLM2
    cuts.AddCut("00000113","1111121063032200000","1111121063022210002","0163302200000000"); // INT7
    cuts.AddCut("00052013","1111121063032200000","1111121063022210002","0163302200000000"); // EMC7
    cuts.AddCut("00085013","1111121063032200000","1111121063022210002","0163302200000000"); // EG2
    cuts.AddCut("00083013","1111121063032200000","1111121063022210002","0163302200000000"); // EG1
  } else if (trainConfig == 41){  // new defaults LHC13g NLM1
    cuts.AddCut("00000113","1111121063032200000","1111121063022210001","0163301100000000"); // INT7
    cuts.AddCut("00052013","1111121063032200000","1111121063022210001","0163301100000000"); // EMC7
    cuts.AddCut("00085013","1111121063032200000","1111121063022210001","0163301100000000"); // EG2
    cuts.AddCut("00083013","1111121063032200000","1111121063022210001","0163301100000000"); // EG1
  } else if (trainConfig == 42){ // NLM 1 M02 var  
    cuts.AddCut("00000113","1111121053032200000","1111121063022110001","0163301100000000"); // min 0.3 function default
    cuts.AddCut("00000113","1111121053032200000","1111121063022310001","0163301100000000"); // min 0.25 function default
    cuts.AddCut("00000113","1111121053032200000","1111121063022410001","0163301100000000"); // min 0.27, tighter lower func
    cuts.AddCut("00000113","1111121053032200000","1111121063022520001","0163301100000000"); // min 0.27, looser func up and low
    cuts.AddCut("00000113","1111121053032200000","1111121063022430001","0163301100000000"); // min 0.27, tighter func up and low
  } else if (trainConfig == 43){ // NLM 2 M02 var
    cuts.AddCut("00000113","1111121053032200000","1111121063022110002","0163302200000000"); // min 0.3 function default
    cuts.AddCut("00000113","1111121053032200000","1111121063022310002","0163302200000000"); // min 0.25 function default
    cuts.AddCut("00000113","1111121053032200000","1111121063022410002","0163302200000000"); // min 0.27, tighter lower func
    cuts.AddCut("00000113","1111121053032200000","1111121063022520002","0163302200000000"); // min 0.27, looser func up and low
    cuts.AddCut("00000113","1111121053032200000","1111121063022430002","0163302200000000"); // min 0.27, tighter func up and low
  } else if (trainConfig == 44){ // NLM 1 M02 var  
    cuts.AddCut("00052013","1111121053032200000","1111121063022110001","0163301100000000"); // min 0.3 function default
    cuts.AddCut("00052013","1111121053032200000","1111121063022310001","0163301100000000"); // min 0.25 function default
    cuts.AddCut("00052013","1111121053032200000","1111121063022410001","0163301100000000"); // min 0.27, tighter lower func
    cuts.AddCut("00052013","1111121053032200000","1111121063022520001","0163301100000000"); // min 0.27, looser func up and low
    cuts.AddCut("00052013","1111121053032200000","1111121063022430001","0163301100000000"); // min 0.27, tighter func up and low
  } else if (trainConfig == 45){ // NLM 2 M02 var
    cuts.AddCut("00052013","1111121053032200000","1111121063022110002","0163302200000000"); // min 0.3 function default
    cuts.AddCut("00052013","1111121053032200000","1111121063022310002","0163302200000000"); // min 0.25 function default
    cuts.AddCut("00052013","1111121053032200000","1111121063022410002","0163302200000000"); // min 0.27, tighter lower func
    cuts.AddCut("00052013","1111121053032200000","1111121063022520002","0163302200000000"); // min 0.27, looser func up and low
    cuts.AddCut("00052013","1111121053032200000","1111121063022430002","0163302200000000"); // min 0.27, tighter func up and low
  } else if (trainConfig == 46){ // NLM 1 M02 var  
    cuts.AddCut("00085013","1111121053032200000","1111121063022110001","0163301100000000"); // min 0.3 function default
    cuts.AddCut("00085013","1111121053032200000","1111121063022310001","0163301100000000"); // min 0.25 function default
    cuts.AddCut("00085013","1111121053032200000","1111121063022410001","0163301100000000"); // min 0.27, tighter lower func
    cuts.AddCut("00085013","1111121053032200000","1111121063022520001","0163301100000000"); // min 0.27, looser func up and low
    cuts.AddCut("00085013","1111121053032200000","1111121063022430001","0163301100000000"); // min 0.27, tighter func up and low
  } else if (trainConfig == 47){ // NLM 2 M02 var
    cuts.AddCut("00085013","1111121053032200000","1111121063022110002","0163302200000000"); // min 0.3 function default
    cuts.AddCut("00085013","1111121053032200000","1111121063022310002","0163302200000000"); // min 0.25 function default
    cuts.AddCut("00085013","1111121053032200000","1111121063022410002","0163302200000000"); // min 0.27, tighter lower func
    cuts.AddCut("00085013","1111121053032200000","1111121063022520002","0163302200000000"); // min 0.27, looser func up and low
    cuts.AddCut("00085013","1111121053032200000","1111121063022430002","0163302200000000"); // min 0.27, tighter func up and low
  } else if (trainConfig == 48){ // NLM 1 M02 var  
    cuts.AddCut("00083013","1111121053032200000","1111121063022110001","0163301100000000"); // min 0.3 function default
    cuts.AddCut("00083013","1111121053032200000","1111121063022310001","0163301100000000"); // min 0.25 function default
    cuts.AddCut("00083013","1111121053032200000","1111121063022410001","0163301100000000"); // min 0.27, tighter lower func
    cuts.AddCut("00083013","1111121053032200000","1111121063022520001","0163301100000000"); // min 0.27, looser func up and low
    cuts.AddCut("00083013","1111121053032200000","1111121063022430001","0163301100000000"); // min 0.27, tighter func up and low
  } else if (trainConfig == 49){ // NLM 2 M02 var
    cuts.AddCut("00083013","1111121053032200000","1111121063022110002","0163302200000000"); // min 0.3 function default
    cuts.AddCut("00083013","1111121053032200000","1111121063022310002","0163302200000000"); // min 0.25 function default
    cuts.AddCut("00083013","1111121053032200000","1111121063022410002","0163302200000000"); // min 0.27, tighter lower func
    cuts.AddCut("00083013","1111121053032200000","1111121063022520002","0163302200000000"); // min 0.27, looser func up and low
    cuts.AddCut("00083013","1111121053032200000","1111121063022430002","0163302200000000"); // min 0.27, tighter func up and low
  } else if (trainConfig == 50){  // EMCAL clusters, variation track matching to cluster & Mass INT7 NLM 1
    cuts.AddCut("00000113","1111121061032200000","1111121061022210001","0163301100000000"); // looser TM
    cuts.AddCut("00000113","1111121066032200000","1111121066022210001","0163301100000000"); // tighter TM
    cuts.AddCut("00000113","1111121063032200000","1111121063022210001","0163301300000000"); // tighter mass cut
    cuts.AddCut("00000113","1111121063032200000","1111121063022210001","0163301500000000"); // looser mass cut
  } else if (trainConfig == 51){  // EMCAL clusters, variation track matching to cluster & Mass INT7 NLM 2
    cuts.AddCut("00000113","1111121061032200000","1111121061022210002","0163302200000000"); // looser TM
    cuts.AddCut("00000113","1111121066032200000","1111121066022210002","0163302200000000"); // tighter TM
    cuts.AddCut("00000113","1111121063032200000","1111121063022210002","0163302400000000"); // tighter mass cut
    cuts.AddCut("00000113","1111121063032200000","1111121063022210002","0163302600000000"); // looser mass cut
  } else if (trainConfig == 52){  // EMCAL clusters, variation track matching to cluster & Mass EMC7 NLM 1
    cuts.AddCut("00052013","1111121061032200000","1111121061022210001","0163301100000000"); // looser TM
    cuts.AddCut("00052013","1111121066032200000","1111121066022210001","0163301100000000"); // tighter TM
    cuts.AddCut("00052013","1111121063032200000","1111121063022210001","0163301300000000"); // tighter mass cut
    cuts.AddCut("00052013","1111121063032200000","1111121063022210001","0163301500000000"); // looser mass cut
  } else if (trainConfig == 53){  // EMCAL clusters, variation track matching to cluster & Mass EMC7 NLM 2
    cuts.AddCut("00052013","1111121061032200000","1111121061022210002","0163302200000000"); // looser TM
    cuts.AddCut("00052013","1111121066032200000","1111121066022210002","0163302200000000"); // tighter TM
    cuts.AddCut("00052013","1111121063032200000","1111121063022210002","0163302400000000"); // tighter mass cut
    cuts.AddCut("00052013","1111121063032200000","1111121063022210002","0163302600000000"); // looser mass cut
  } else if (trainConfig == 54){  // EMCAL clusters, variation track matching to cluster & Mass EG2 NLM 1
    cuts.AddCut("00085013","1111121061032200000","1111121061022210001","0163301100000000"); // looser TM
    cuts.AddCut("00085013","1111121066032200000","1111121066022210001","0163301100000000"); // tighter TM
    cuts.AddCut("00085013","1111121063032200000","1111121063022210001","0163301300000000"); // tighter mass cut
    cuts.AddCut("00085013","1111121063032200000","1111121063022210001","0163301500000000"); // looser mass cut
  } else if (trainConfig == 55){  // EMCAL clusters, variation track matching to cluster & Mass EG2 NLM 2
    cuts.AddCut("00085013","1111121061032200000","1111121061022210002","0163302200000000"); // looser TM
    cuts.AddCut("00085013","1111121066032200000","1111121066022210002","0163302200000000"); // tighter TM
    cuts.AddCut("00085013","1111121063032200000","1111121063022210002","0163302400000000"); // tighter mass cut
    cuts.AddCut("00085013","1111121063032200000","1111121063022210002","0163302600000000"); // looser mass cut
  } else if (trainConfig == 56){  // EMCAL clusters, variation track matching to cluster & Mass EG1 NLM 1
    cuts.AddCut("00083013","1111121061032200000","1111121061022210001","0163301100000000"); // looser TM
    cuts.AddCut("00083013","1111121066032200000","1111121066022210001","0163301100000000"); // tighter TM
    cuts.AddCut("00083013","1111121063032200000","1111121063022210001","0163301300000000"); // tighter mass cut
    cuts.AddCut("00083013","1111121063032200000","1111121063022210001","0163301500000000"); // looser mass cut
  } else if (trainConfig == 57){  // EMCAL clusters, variation track matching to cluster & Mass EG1 NLM 2
    cuts.AddCut("00083013","1111121061032200000","1111121061022210002","0163302200000000"); // looser TM
    cuts.AddCut("00083013","1111121066032200000","1111121066022210002","0163302200000000"); // tighter TM
    cuts.AddCut("00083013","1111121063032200000","1111121063022210002","0163302400000000"); // tighter mass cut
    cuts.AddCut("00083013","1111121063032200000","1111121063022210002","0163302600000000"); // looser mass cut
  } else if (trainConfig == 58){  // NL variations INT7
    cuts.AddCut("00000113","1111122063032200000","1111122063022210001","0163301100000000"); // SDM loose time
    cuts.AddCut("00000113","1111111063032200000","1111111063022210001","0163301100000000"); // conv calo tight time
    cuts.AddCut("00000113","1111101063032200000","1111101063022210001","0163301100000000"); // SDM Jason
    cuts.AddCut("00000113","1111100063032200000","1111100063022210001","0163301100000000"); // none
  } else if (trainConfig == 59){  // NL variations INT7
    cuts.AddCut("00000113","1111122063032200000","1111122063022210002","0163302200000000"); // SDM loose time
    cuts.AddCut("00000113","1111111063032200000","1111111063022210002","0163302200000000"); // conv calo tight time
    cuts.AddCut("00000113","1111101063032200000","1111101063022210002","0163302200000000"); // SDM Jason
    cuts.AddCut("00000113","1111100063032200000","1111100063022210002","0163302200000000"); // none
  } else if (trainConfig == 60){  // NL variations EMC7
    cuts.AddCut("00052013","1111122063032200000","1111122063022210001","0163301100000000"); // SDM loose time
    cuts.AddCut("00052013","1111111063032200000","1111111063022210001","0163301100000000"); // conv calo tight time
    cuts.AddCut("00052013","1111101063032200000","1111101063022210001","0163301100000000"); // SDM Jason
    cuts.AddCut("00052013","1111100063032200000","1111100063022210001","0163301100000000"); // none
  } else if (trainConfig == 61){  // NL variations EMC7
    cuts.AddCut("00052013","1111122063032200000","1111122063022210002","0163302200000000"); // SDM loose time
    cuts.AddCut("00052013","1111111063032200000","1111111063022210002","0163302200000000"); // conv calo tight time
    cuts.AddCut("00052013","1111101063032200000","1111101063022210002","0163302200000000"); // SDM Jason
    cuts.AddCut("00052013","1111100063032200000","1111100063022210002","0163302200000000"); // none
  } else if (trainConfig == 62){  // NL variations EG2
    cuts.AddCut("00085013","1111122063032200000","1111122063022210001","0163301100000000"); // SDM loose time
    cuts.AddCut("00085013","1111111063032200000","1111111063022210001","0163301100000000"); // conv calo tight time
    cuts.AddCut("00085013","1111101063032200000","1111101063022210001","0163301100000000"); // SDM Jason
    cuts.AddCut("00085013","1111100063032200000","1111100063022210001","0163301100000000"); // none
  } else if (trainConfig == 63){  // NL variations EG2
    cuts.AddCut("00085013","1111122063032200000","1111122063022210002","0163302200000000"); // SDM loose time
    cuts.AddCut("00085013","1111111063032200000","1111111063022210002","0163302200000000"); // conv calo tight time
    cuts.AddCut("00085013","1111101063032200000","1111101063022210002","0163302200000000"); // SDM Jason
    cuts.AddCut("00085013","1111100063032200000","1111100063022210002","0163302200000000"); // none
  } else if (trainConfig == 64){  // NL variations EG1
    cuts.AddCut("00083013","1111122063032200000","1111122063022210001","0163301100000000"); // SDM loose time
    cuts.AddCut("00083013","1111111063032200000","1111111063022210001","0163301100000000"); // conv calo tight time
    cuts.AddCut("00083013","1111101063032200000","1111101063022210001","0163301100000000"); // SDM Jason
    cuts.AddCut("00083013","1111100063032200000","1111100063022210001","0163301100000000"); // none
  } else if (trainConfig == 65){  // NL variations EG1
    cuts.AddCut("00083013","1111122063032200000","1111122063022210002","0163302200000000"); // SDM loose time
    cuts.AddCut("00083013","1111111063032200000","1111111063022210002","0163302200000000"); // conv calo tight time
    cuts.AddCut("00083013","1111101063032200000","1111101063022210002","0163302200000000"); // SDM Jason
    cuts.AddCut("00083013","1111100063032200000","1111100063022210002","0163302200000000"); // none
  } else if (trainConfig == 66){  // Alpha cut variations & TRD material INT7 NLM 1
    cuts.AddCut("00000113","1111121063032200000","1111121063022210001","0163303100000000"); // NLM 1 looser
    cuts.AddCut("00000113","1111121063032200000","1111121063022210001","0163305100000000"); // NLM 1 tighter
    cuts.AddCut("00000113","1112121063032200000","1112121063022210001","0163301100000000"); // TRD infront
    cuts.AddCut("00000113","1111321063032200000","1111321063022210001","0163301100000000"); // no TRD infront    
  } else if (trainConfig == 67){  // Alpha cut variations & TRD material INT7 NLM 2
    cuts.AddCut("00000113","1111121063032200000","1111121063022210002","0163304200000000"); // NLM 2 looser
    cuts.AddCut("00000113","1111121063032200000","1111121063022210002","0163306200000000"); // NLM 2 tighter
    cuts.AddCut("00000113","1112121063032200000","1112121063022210002","0163302200000000"); // TRD infront
    cuts.AddCut("00000113","1111321063032200000","1111321063022210002","0163302200000000"); // no TRD infront
  } else if (trainConfig == 68){  // Alpha cut variations & TRD material EMC7 NLM 1
    cuts.AddCut("00052013","1111121063032200000","1111121063022210001","0163303100000000"); // NLM 1 looser
    cuts.AddCut("00052013","1111121063032200000","1111121063022210001","0163305100000000"); // NLM 1 tighter
    cuts.AddCut("00052013","1112121063032200000","1112121063022210001","0163301100000000"); // TRD infront
    cuts.AddCut("00052013","1111321063032200000","1111321063022210001","0163301100000000"); // no TRD infront
  } else if (trainConfig == 69){  // Alpha cut variations & TRD material EMC7 NLM 2
    cuts.AddCut("00052013","1111121063032200000","1111121063022210002","0163304200000000"); // NLM 2 looser
    cuts.AddCut("00052013","1111121063032200000","1111121063022210002","0163306200000000"); // NLM 2 tighter
    cuts.AddCut("00052013","1112121063032200000","1112121063022210002","0163302200000000"); // TRD infront
    cuts.AddCut("00052013","1111321063032200000","1111321063022210002","0163302200000000"); // no TRD infront
  } else if (trainConfig == 70){  // Alpha cut variations & TRD material EG2 NLM 1
    cuts.AddCut("00085013","1111121063032200000","1111121063022210001","0163303100000000"); // NLM 1 looser
    cuts.AddCut("00085013","1111121063032200000","1111121063022210001","0163305100000000"); // NLM 1 tighter
    cuts.AddCut("00085013","1112121063032200000","1112121063022210001","0163301100000000"); // TRD infront
    cuts.AddCut("00085013","1111321063032200000","1111321063022210001","0163301100000000"); // no TRD infront
  } else if (trainConfig == 71){  // Alpha cut variations & TRD material EG2 NLM 2
    cuts.AddCut("00085013","1111121063032200000","1111121063022210002","0163304200000000"); // NLM 2 looser
    cuts.AddCut("00085013","1111121063032200000","1111121063022210002","0163306200000000"); // NLM 2 tighter
    cuts.AddCut("00085013","1112121063032200000","1112121063022210002","0163302200000000"); // TRD infront
    cuts.AddCut("00085013","1111321063032200000","1111321063022210002","0163302200000000"); // no TRD infront
  } else if (trainConfig == 72){  // Alpha cut variations & TRD material EG1 NLM 1
    cuts.AddCut("00083013","1111121063032200000","1111121063022210001","0163303100000000"); // NLM 1 looser
    cuts.AddCut("00083013","1111121063032200000","1111121063022210001","0163305100000000"); // NLM 1 tighter
    cuts.AddCut("00083013","1112121063032200000","1112121063022210001","0163301100000000"); // TRD infront
    cuts.AddCut("00083013","1111321063032200000","1111321063022210001","0163301100000000"); // no TRD infront
  } else if (trainConfig == 73){  // Alpha cut variations & TRD material EG1 NLM 2
    cuts.AddCut("00083013","1111121063032200000","1111121063022210002","0163304200000000"); // NLM 2 looser
    cuts.AddCut("00083013","1111121063032200000","1111121063022210002","0163306200000000"); // NLM 2 tighter
    cuts.AddCut("00083013","1112121063032200000","1112121063022210002","0163302200000000"); // TRD infront
    cuts.AddCut("00083013","1111321063032200000","1111321063022210002","0163302200000000"); // no TRD infront
  } else if (trainConfig == 74){  // no TM in basis cuts LHC13g NLM2
    cuts.AddCut("00000113","1111121060032200000","1111121063022210002","0163302200000000"); // INT7
    cuts.AddCut("00052013","1111121060032200000","1111121063022210002","0163302200000000"); // EMC7
    cuts.AddCut("00085013","1111121060032200000","1111121063022210002","0163302200000000"); // EG2
    cuts.AddCut("00083013","1111121060032200000","1111121063022210002","0163302200000000"); // EG1
  } else if (trainConfig == 75){  // no TM in basis cuts LHC13g NLM1
    cuts.AddCut("00000113","1111121060032200000","1111121063022210001","0163301100000000"); // INT7
    cuts.AddCut("00052013","1111121060032200000","1111121063022210001","0163301100000000"); // EMC7
    cuts.AddCut("00085013","1111121060032200000","1111121063022210001","0163301100000000"); // EG2
    cuts.AddCut("00083013","1111121060032200000","1111121063022210001","0163301100000000"); // EG1
  } else if (trainConfig == 76){  // NLM2 no M02, no mass, no alpah
    cuts.AddCut("00000113","1111121063032200000","1111121063022000002","0163300000000000"); // INT7
    cuts.AddCut("00052013","1111121063032200000","1111121063022000002","0163300000000000"); // EMC7
    cuts.AddCut("00085013","1111121063032200000","1111121063022000002","0163300000000000"); // EG2
    cuts.AddCut("00083013","1111121063032200000","1111121063022000002","0163300000000000"); // EG1
  } else if (trainConfig == 77){  // NLM1 no M02, no mass, no alpah
    cuts.AddCut("00000113","1111121063032200000","1111121063022000001","0163300000000000"); // INT7
    cuts.AddCut("00052013","1111121063032200000","1111121063022000001","0163300000000000"); // EMC7
    cuts.AddCut("00085013","1111121063032200000","1111121063022000001","0163300000000000"); // EG2
    cuts.AddCut("00083013","1111121063032200000","1111121063022000001","0163300000000000"); // EG1
    
    // LHC12
    // default with three cuts
  } else if (trainConfig == 111){  // EMCAL clusters, different triggers no NonLinearity
    cuts.AddCut("00000113","1111100060032200000","1111100060022110001","0163301100000000"); // INT7
    cuts.AddCut("00052013","1111100060032200000","1111100060022110001","0163301100000000"); // EMC7
    cuts.AddCut("00081113","1111100060032200000","1111100060022110001","0163301100000000"); // EMCEGA,
  } else if (trainConfig == 112){  // EMCAL clusters, different triggers no NonLinearity
    cuts.AddCut("00000113","1111100060032200000","1111100060022110002","0163302200000000"); // INT7
    cuts.AddCut("00052013","1111100060032200000","1111100060022110002","0163302200000000"); // EMC7
    cuts.AddCut("00081113","1111100060032200000","1111100060022110002","0163302200000000"); // EMCEGA,
  } else if (trainConfig == 113){  // EMCAL clusters, different triggers with NonLinearity
    cuts.AddCut("00000113","1111111060032200000","1111111060022110001","0163301100000000"); // INT7
    cuts.AddCut("00052013","1111111060032200000","1111111060022110001","0163301100000000"); // EMC7
    cuts.AddCut("00081113","1111111060032200000","1111111060022110001","0163301100000000"); // EMCEGA,
  } else if (trainConfig == 114){  // EMCAL clusters, different triggers with NonLinearity
    cuts.AddCut("00000113","1111111060032200000","1111111060022110002","0163302200000000"); // INT7
    cuts.AddCut("00052013","1111111060032200000","1111111060022110002","0163302200000000"); // EMC7
    cuts.AddCut("00081113","1111111060032200000","1111111060022110002","0163302200000000"); // EMCEGA,
  } else if (trainConfig == 115){  // EMCAL clusters, different triggers with NonLinearity track matching to cluster
    cuts.AddCut("00000113","1111111063032200000","1111111063022110001","0163301100000000"); // INT7
    cuts.AddCut("00052013","1111111063032200000","1111111063022110001","0163301100000000"); // EMC7
    cuts.AddCut("00081113","1111111063032200000","1111111063022110001","0163301100000000"); // EMCEGA,
  } else if (trainConfig == 116){  // EMCAL clusters, different triggers with NonLinearity track matching to cluster
    cuts.AddCut("00000113","1111111063032200000","1111111063022110002","0163302200000000"); // INT7
    cuts.AddCut("00052013","1111111063032200000","1111111063022110002","0163302200000000"); // EMC7
    cuts.AddCut("00081113","1111111063032200000","1111111063022110002","0163302200000000"); // EMCEGA,
  } else if (trainConfig == 117){  // EMCAL clusters, different triggers with NonLinearity track matching to cluster
    cuts.AddCut("00011113","1111111063032200000","1111111063022110001","0163301100000000"); // INT8
    cuts.AddCut("00053113","1111111063032200000","1111111063022110001","0163301100000000"); // EMC8
    cuts.AddCut("00082113","1111111063032200000","1111111063022110001","0163301100000000"); // EMC8EGA,
  } else if (trainConfig == 118){  // EMCAL clusters, different triggers with NonLinearity track matching to cluster
    cuts.AddCut("00011113","1111111063032200000","1111111063022110002","0163302200000000"); // INT8
    cuts.AddCut("00053113","1111111063032200000","1111111063022110002","0163302200000000"); // EMC8
    cuts.AddCut("00082113","1111111063032200000","1111111063022110002","0163302200000000"); // EMC8EGA,

 
    // 13 TeV & 5 TeV
  } else if (trainConfig == 401){ // EMCAL clusters pp 13 TeV
    cuts.AddCut("00000113","1111101013032200000","1111101013022210002","0163302200000000"); //INT7 1000ns timing cut, std NL NLM2
    cuts.AddCut("00000113","1111101013032200000","1111101013022210002","0163302200000000"); //INT7 1000ns timing cut, std NL NLM2
    cuts.AddCut("00000113","1111101013032200000","1111101013022110001","0163301100000000"); //INT7 1000ns timing cut, std NL NLM1
    cuts.AddCut("00052013","1111101013032200000","1111101013022110001","0163301100000000"); //EMC7 1000ns timing cut, std NL NLM1
 
    
    // all the cut variations
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

  Int_t numberOfCuts          = cuts.GetNCuts();

  TList *EventCutList         = new TList();
  TList *ClusterCutList       = new TList();
  TList *ClusterMergedCutList = new TList();
  TList *MesonCutList         = new TList();

  TList *HeaderList           = new TList();
  if (periodname.Contains("LHC12i3")){  
    TObjString *Header2       = new TObjString("BOX");
    HeaderList->Add(Header2);
  } else if (periodname.CompareTo("LHC14e2b")==0){
    TObjString *Header2       = new TObjString("pi0_1");
    HeaderList->Add(Header2);
    TObjString *Header3       = new TObjString("eta_2");
    HeaderList->Add(Header3);
  }  
  
  TString energy      = "";
  TString mcName      = "";
  TString mcNameAdd   = "";
  if (periodname.Contains("WOSDD")){
    mcNameAdd         = "_WOSDD";
  } else if (periodname.Contains("WSDD")){
    mcNameAdd         = "_WSDD";
  }   
  if (periodname.Contains("LHC12i3")){
    energy            = "2760GeV";
    mcName            = "Pythia8_LHC12i3";
  } else if (periodname.Contains("LHC12f1a")){  
    energy            = "2760GeV";
    mcName            = "Pythia8_LHC12f1a";
  } else if (periodname.Contains("LHC12f1b")){  
    energy            = "2760GeV";
    mcName            = "Phojet_LHC12f1b";    
  } else if (periodname.Contains("LHC14e2a")){  
    energy            = "8TeV";
    mcName            = "Pythia8_LHC14e2a";    
  } else if (periodname.Contains("LHC14e2b")){  
    energy            = "8TeV";
    mcName            = "Pythia8_LHC14e2b";      
  } else if (periodname.Contains("LHC14e2c")){    
    energy            = "8TeV";
    mcName            = "Phojet_LHC14e2c";        
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
    analysisEventCuts[i]    = new AliConvEventCuts();
    
    // definition of weighting input
    TString fitNamePi0      = Form("Pi0_Fit_Data_%s",energy.Data());
    TString fitNameEta      = Form("Eta_Fit_Data_%s",energy.Data());

    TString mcInputNamePi0  = "";
    TString mcInputNameEta  = "";
    mcInputNamePi0          = Form("Pi0_%s%s_%s", mcName.Data(), mcNameAdd.Data(), energy.Data() );
    mcInputNameEta          = Form("Eta_%s%s_%s", mcName.Data(), mcNameAdd.Data(), energy.Data() );
    
    if (doWeighting) analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE, fileNameInputForWeighting, mcInputNamePi0, mcInputNameEta, "",fitNamePi0,fitNameEta);

    analysisEventCuts[i]->SetTriggerMimicking(enableTriggerMimicking);
    analysisEventCuts[i]->SetTriggerOverlapRejecion(enableTriggerOverlapRej);
    analysisEventCuts[i]->SetMaxFacPtHard(maxFacPtHard);
    analysisEventCuts[i]->SetV0ReaderName(V0ReaderName);
    analysisEventCuts[i]->InitializeCutsFromCutString((cuts.GetEventCut(i)).Data());
    EventCutList->Add(analysisEventCuts[i]);
    analysisEventCuts[i]->SetFillCutHistograms("",kFALSE);
    
    analysisClusterCuts[i]        = new AliCaloPhotonCuts();
    analysisClusterCuts[i]->SetIsPureCaloCut(2);
    analysisClusterCuts[i]->SetV0ReaderName(V0ReaderName);
    analysisClusterCuts[i]->InitializeCutsFromCutString((cuts.GetClusterCut(i)).Data());
    ClusterCutList->Add(analysisClusterCuts[i]);
    analysisClusterCuts[i]->SetExtendedMatchAndQA(enableExtMatchAndQA);
    analysisClusterCuts[i]->SetFillCutHistograms("");

    analysisClusterMergedCuts[i]  = new AliCaloPhotonCuts();
    analysisClusterMergedCuts[i]->SetIsPureCaloCut(1);
    analysisClusterMergedCuts[i]->SetV0ReaderName(V0ReaderName);
    analysisClusterMergedCuts[i]->InitializeCutsFromCutString((cuts.GetClusterMergedCut(i)).Data());
    ClusterMergedCutList->Add(analysisClusterMergedCuts[i]);
    analysisClusterMergedCuts[i]->SetExtendedMatchAndQA(enableExtMatchAndQA);
    analysisClusterMergedCuts[i]->SetFillCutHistograms("");

    analysisMesonCuts[i]          = new AliConversionMesonCuts();
    analysisMesonCuts[i]->SetEnableOpeningAngleCut(kFALSE);
    analysisMesonCuts[i]->SetIsMergedClusterCut(1);
    analysisMesonCuts[i]->InitializeCutsFromCutString((cuts.GetMesonCut(i)).Data());
    MesonCutList->Add(analysisMesonCuts[i]);
    analysisMesonCuts[i]->SetFillCutHistograms("");
    analysisEventCuts[i]->SetAcceptedHeader(HeaderList);
  }
  task->SetSelectedMesonID(selectedMeson);
  task->SetEventCutList(numberOfCuts,EventCutList);
  task->SetCaloCutList(numberOfCuts,ClusterCutList);
  task->SetCaloMergedCutList(numberOfCuts,ClusterMergedCutList);
  task->SetMesonCutList(numberOfCuts,MesonCutList);
  task->SetDoMesonQA(enableQAMesonTask); //Attention new switch for Pi0 QA
  task->SetDoClusterQA(enableQAClusterTask);//Attention new switch small for Cluster QA
  if(enableExtMatchAndQA == 2 || enableExtMatchAndQA == 3){ task->SetPlotHistsExtQA(kTRUE);}
  
  //connect containers
  AliAnalysisDataContainer *coutput =
    mgr->CreateContainer(Form("GammaCaloMerged_%i",trainConfig), TList::Class(),
              AliAnalysisManager::kOutputContainer,Form("GammaCaloMerged_%i.root",trainConfig));

  mgr->AddTask(task);
  mgr->ConnectInput(task, 0, cinput);
  mgr->ConnectOutput(task, 1, coutput);

  return;
}
