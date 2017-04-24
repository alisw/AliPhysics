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
                                  Int_t     isMC                        = 0,                  // run MC
                                  Int_t     enableQAMesonTask           = 0,                  // enable QA in AliAnalysisTaskGammaCalo
                                  Int_t     enableQAClusterTask         = 0,                  // enable additional QA task
                                  TString   fileNameInputForWeighting   = "MCSpectraInput.root",       // path to file for weigting input / modified acceptance
                                  TString   cutnumberAODBranch          = "000000006008400001001500000",
                                  TString   periodname                  = "LHC12f1x",         // period name
                                  Bool_t    doWeighting                 = kFALSE,             // enables weighting
                                  Int_t     enableExtMatchAndQA         = 0,                  // disabled (0), extMatch (1), extQA_noCellQA (2), extMatch+extQA_noCellQA (3), extQA+cellQA (4), extMatch+extQA+cellQA (5)
                                  Bool_t    enableTriggerMimicking      = kFALSE,             // enable trigger mimicking
                                  Bool_t    enableTriggerOverlapRej     = kFALSE,             // enable trigger overlap rejection
                                  Float_t   maxFacPtHard                = 3.,                 // maximum factor between hardest jet and ptHard generated
                                  TString   periodNameV0Reader          = "",                 // period Name for respective period selected in V0Reader
                                  Int_t     selectedMeson               = 1,                  // put flag for selected meson
                                  Bool_t    enableDetailedPrintout      = kFALSE,             // enable detailed printout
                                  Bool_t    enableSortingMCLabels       = kTRUE,              // enable sorting for MC cluster labels
                                  Bool_t    runLightOutput              = kFALSE,             // switch to run light output (only essential histograms for afterburner)
                                  Double_t  minEnergyForExoticsCut      = 1.0,                // minimum energy to be used for exotics CutHandler
                                  Bool_t    runQAForExotics             = kFALSE,             // switch to run QA for exotic clusters
                                  TString   additionalTrainConfig       = "0"                 // additional counter for trainconfig, always has to be last parameter
) {
  
  TH1S* histoAcc = 0x0;         // histo for modified acceptance
  //parse additionalTrainConfig flag
  TObjArray *rAddConfigArr = additionalTrainConfig.Tokenize("_");
  if(rAddConfigArr->GetEntries()<1){cout << "ERROR: AddTask_GammaCaloMerged_pp during parsing of additionalTrainConfig String '" << additionalTrainConfig.Data() << "'" << endl; return;}
  TObjString* rAdditionalTrainConfig;
  for(Int_t i = 0; i<rAddConfigArr->GetEntries() ; i++){
    if(i==0) rAdditionalTrainConfig = (TObjString*)rAddConfigArr->At(i);
    else{
      TObjString* temp = (TObjString*) rAddConfigArr->At(i);
      TString tempStr = temp->GetString();
      if(tempStr.CompareTo("EPCLUSTree") == 0){
        cout << "INFO: AddTask_GammaCaloMerged_pp activating 'EPCLUSTree'" << endl;
        doTreeEOverP = kTRUE;
      }else if(tempStr.BeginsWith("MODIFYACC")){
        cout << "INFO: AddTask_GammaCaloMerged_pp activating 'MODIFYACC'" << endl;
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
      }
    }
  }

  if (additionalTrainConfig.Atoi() > 0){
    trainConfig = trainConfig + additionalTrainConfig.Atoi();
    cout << "INFO: AddTask_GammaCaloMerged_pp running additionalTrainConfig '" << sAdditionalTrainConfig.Atoi() << "', train config: '" << trainConfig << "'" << endl;
  }  

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
    if(trainConfig>=100 && trainConfig<200) fV0ReaderV1->SetImprovedPsiPair(0); //switch off for 8TeV as AODs are used for which improved psipair is not available

    if (!mgr) {
      Error("AddTask_V0ReaderV1", "No analysis manager found.");
      return;
    }
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
  // LHC11a
  if (trainConfig == 1){ // all defaults for LHC11a
    cuts.AddCut("00003113","1111121053032200000","1111121053022210001","0163301100000000"); // INT1
    cuts.AddCut("00051013","1111121053032200000","1111121053022210001","0163301100000000"); // EMC1
  } else if (trainConfig == 2){ // no TM in basis cut
    cuts.AddCut("00003113","1111121050032200000","1111121053022210001","0163301100000000"); // INT1
    cuts.AddCut("00051013","1111121050032200000","1111121053022210001","0163301100000000"); // EMC1
  } else if (trainConfig == 3){ // open M02, open Mass, open Alpha
    cuts.AddCut("00003113","1111121053032200000","1111121053022000001","0163300000000000"); // INT1
    cuts.AddCut("00051013","1111121053032200000","1111121053022000001","0163300000000000"); // EMC1
  } else if (trainConfig == 4){  //NLM exotics default frac = 0.97
    cuts.AddCut("00003113","1111121053532200000","1111121053522110001","0163301100000000"); // INT1
    cuts.AddCut("00051013","1111121053532200000","1111121053522110001","0163301100000000"); // EMC1
  } else if (trainConfig == 5){  // new default 
    cuts.AddCut("00003113","1111121053032200000","1111121053022700001","0163300700000000"); // Mass only band at 0, no explicit exotics cut, M02 cut at 0.27  
    cuts.AddCut("00051013","1111121053032200000","1111121053022700001","0163300700000000"); // Mass only band at 0, no explicit exotics cut, M02 cut at 0.27  
    cuts.AddCut("00003113","1111121053032200000","1111121053022700001","0163300000000000"); // INT1
    cuts.AddCut("00051013","1111121053032200000","1111121053022700001","0163300000000000"); // EMC1
  } else if (trainConfig == 6){  // new default 
    cuts.AddCut("00003113","1111121050032200000","1111121053022000001","0163300000000000"); // no M02, no exotics, TM only for second cut
    cuts.AddCut("00051013","1111121050032200000","1111121053022000001","0163300000000000"); // no M02, no exotics, TM only for second cut
    cuts.AddCut("00003113","1111121050032200000","1111121053022700001","0163300000000000"); // M02 < 0.27, no exotics, TM only for second cut
    cuts.AddCut("00051013","1111121050032200000","1111121053022700001","0163300000000000"); // M02 < 0.27, no exotics, TM only for second cut
  } else if (trainConfig == 7){  // new default , pt dep TM
    cuts.AddCut("00003113","1111121057032200000","1111121057022700001","0163300000000000"); // M02 < 0.27, no exotics
    cuts.AddCut("00051013","1111121057032200000","1111121057022700001","0163300000000000"); // M02 < 0.27, no exotics
    cuts.AddCut("00003113","1111121050032200000","1111121057022700001","0163300000000000"); // M02 < 0.27, no exotics, TM only for second
    cuts.AddCut("00051013","1111121050032200000","1111121057022700001","0163300000000000"); // M02 < 0.27, no exotics, TM only for second
  } else if (trainConfig == 8){  // new default, with eta < 0.7, y < 0.7
    cuts.AddCut("00003113","1551121053032200000","1551121053022700001","0163200000000000"); // Mass only band at 0, no explicit exotics cut, M02 cut at 0.27  
    cuts.AddCut("00051013","1551121053032200000","1551121053022700001","0163200000000000"); // Mass only band at 0, no explicit exotics cut, M02 cut at 0.27  
    cuts.AddCut("00003113","1551121057032200000","1551121057022700001","0163200000000000"); // Mass only band at 0, no explicit exotics cut, M02 cut at 0.27, pt dep TM  
    cuts.AddCut("00051013","1551121057032200000","1551121057022700001","0163200000000000"); // Mass only band at 0, no explicit exotics cut, M02 cut at 0.27, pt dep TM  

  // INT1 variations
  } else if (trainConfig == 10){ // M02 var  
    cuts.AddCut("00003113","1111121053032200000","1111121053022110001","0163301100000000"); // min 0.3 function default
    cuts.AddCut("00003113","1111121053032200000","1111121053022310001","0163301100000000"); // min 0.25 function default
    cuts.AddCut("00003113","1111121053032200000","1111121053022410001","0163301100000000"); // min 0.27, tighter lower func
    cuts.AddCut("00003113","1111121053032200000","1111121053022520001","0163301100000000"); // min 0.27, looser func up and low
    cuts.AddCut("00003113","1111121053032200000","1111121053022430001","0163301100000000"); // min 0.27, tighter func up and low
  } else if (trainConfig == 11){  // variation track matching to cluster & mass variations 
    cuts.AddCut("00003113","1111121050032200000","1111121050022210001","0163301100000000"); // no TM
    cuts.AddCut("00003113","1111121051032200000","1111121051022210001","0163301100000000"); // looser TM
    cuts.AddCut("00003113","1111121056032200000","1111121056022210001","0163301100000000"); // tighter TM
    cuts.AddCut("00003113","1111121053032200000","1111121053022210001","0163301300000000"); // tighter mass cut
    cuts.AddCut("00003113","1111121053032200000","1111121053022210001","0163301500000000"); // looser mass cut
  } else if (trainConfig == 12){ // NL var
    cuts.AddCut("00003113","1111122053032200000","1111122053022210001","0163301100000000"); // SDM loose time
    cuts.AddCut("00003113","1111111053032200000","1111111053022210001","0163301100000000"); // conv calo tight time
    cuts.AddCut("00003113","1111101053032200000","1111101053022210001","0163301100000000"); // SDM Jason
    cuts.AddCut("00003113","111110005303O2200000","1111100053022210001","0163301100000000"); // none
  } else if (trainConfig == 13){ // Alpha cut variations & TRD material 
    cuts.AddCut("00003113","1111121053032200000","1111121053022210001","0163303100000000"); // NLM 1 looser
    cuts.AddCut("00003113","1111121053032200000","1111121053022210001","0163305100000000"); // NLM 1 tighter
    cuts.AddCut("00003113","1113121053032200000","1113121053022210001","0163301100000000");// TRD infront
    cuts.AddCut("00003113","1111221053032200000","1111221053022210001","0163301100000000");// no TRD infront
  } else if (trainConfig == 14){ // varied exoctics
    cuts.AddCut("00003113","1111121053032200000","1111121053022700001","0163300000000000"); // no exotics
    cuts.AddCut("00003113","1111121053232200000","1111121053222700001","0163300000000000"); // frac = 0.99
    cuts.AddCut("00003113","1111121053332200000","1111121053322700001","0163300000000000"); // frac = 0.98
    cuts.AddCut("00003113","1111121053532200000","1111121053522700001","0163300000000000"); // frac = 0.97
    cuts.AddCut("00003113","1111121053732200000","1111121053722700001","0163300000000000"); // frac = 0.96
    cuts.AddCut("00003113","1111121053932200000","1111121053922700001","0163300000000000"); // frac = 0.94
  } else if (trainConfig == 15){  // variation of mass and alpha cut 
    cuts.AddCut("00003113","1111121053032200000","1111121053022700001","0163301700000000"); // only band at 0
    cuts.AddCut("00003113","1111121053032200000","1111121053022700001","0163300700000000"); // only band at 0, no alpha
    cuts.AddCut("00003113","1111121053032200000","1111121053022700001","0163301000000000"); // no mass cut
    cuts.AddCut("00003113","1111121053032200000","1111121053022700001","0163300000000000"); // no mass cut, no alpha
    cuts.AddCut("00003113","1111121053032200000","1111121053022700001","0163300100000000"); // no alpha
  } else if (trainConfig == 16){ // varied M02
    cuts.AddCut("00003113","1111121053032200000","1111121053022600001","0163300700000000"); // min M02= 0.3
    cuts.AddCut("00003113","1111121053032200000","1111121053022800001","0163300700000000"); // min M02= 0.25
    cuts.AddCut("00003113","1111121053032200000","1111121053022600001","0163300000000000"); // no mass, min M02 = 0.3
    cuts.AddCut("00003113","1111121053032200000","1111121053022700001","0163300000000000"); // no mass, min M02 = 0.27
    cuts.AddCut("00003113","1111121053032200000","1111121053022800001","0163300000000000"); // no mass, min M02 = 0.25
  } else if (trainConfig == 17){  // variation track matching to cluster & mass variations new defaults
    cuts.AddCut("00003113","1111121050032200000","1111121050022700001","0163300000000000"); // no TM
    cuts.AddCut("00003113","1111121051032200000","1111121051022700001","0163300000000000"); // looser TM
    cuts.AddCut("00003113","1111121056032200000","1111121056022700001","0163300000000000"); // tighter TM
    cuts.AddCut("00003113","1113121053032200000","1113121053022700001","0163300000000000");// TRD infront
    cuts.AddCut("00003113","1111221053032200000","1111221053022700001","0163300000000000");// no TRD infront
  } else if (trainConfig == 18){ // NL var new defaults
    cuts.AddCut("00003113","1111122053032200000","1111122053022700001","0163300000000000"); // SDM loose time
    cuts.AddCut("00003113","1111111053032200000","1111111053022700001","0163300000000000"); // conv calo tight time
    cuts.AddCut("00003113","1111101053032200000","1111101053022700001","0163300000000000"); // SDM Jason
    cuts.AddCut("00003113","1111100053032200000","1111100053022700001","0163300000000000"); // none
    
  // EMC1 variations  
  } else if (trainConfig == 20){ // M02 var  
    cuts.AddCut("00051013","1111121053032200000","1111121053022110001","0163301100000000"); // min 0.3 function default
    cuts.AddCut("00051013","1111121053032200000","1111121053022310001","0163301100000000"); // min 0.25 function default
    cuts.AddCut("00051013","1111121053032200000","1111121053022410001","0163301100000000"); // min 0.27, tighter lower func
    cuts.AddCut("00051013","1111121053032200000","1111121053022520001","0163301100000000"); // min 0.27, looser func up and low
    cuts.AddCut("00051013","1111121053032200000","1111121053022430001","0163301100000000"); // min 0.27, tighter func up and low
  } else if (trainConfig == 21){ // variation track matching to cluster & mass variations 
    cuts.AddCut("00051013","1111121050032200000","1111121050022210001","0163301100000000"); // no TM
    cuts.AddCut("00051013","1111121051032200000","1111121051022210001","0163301100000000"); // looser TM
    cuts.AddCut("00051013","1111121056032200000","1111121056022210001","0163301100000000"); // tighter TM
    cuts.AddCut("00051013","1111121053032200000","1111121053022210001","0163301300000000"); // tighter mass cut
    cuts.AddCut("00051013","1111121053032200000","1111121053022210001","0163301500000000"); // looser mass cut    
  } else if (trainConfig == 22){ // NL var
    cuts.AddCut("00051013","1111122053032200000","1111122053022210001","0163301100000000"); // SDM loose time
    cuts.AddCut("00051013","1111111053032200000","1111111053022210001","0163301100000000"); // conv calo tight time
    cuts.AddCut("00051013","1111101053032200000","1111101053022210001","0163301100000000"); // SDM Jason
    cuts.AddCut("00051013","1111100053032200000","1111100053022210001","0163301100000000"); // none
  } else if (trainConfig == 23){ // Alpha cut variations & TRD material
    cuts.AddCut("00051013","1111121053032200000","1111121053022210001","0163303100000000"); // NLM 1 looser
    cuts.AddCut("00051013","1111121053032200000","1111121053022210001","0163305100000000"); // NLM 1 tighter
    cuts.AddCut("00051013","1113121053032200000","1113121053022210001","0163301100000000");// TRD infront
    cuts.AddCut("00051013","1111221053032200000","1111221053022210001","0163301100000000");// no TRD infront
  } else if (trainConfig == 24){ // varied exoctics
    cuts.AddCut("00051013","1111121053032200000","1111121053022700001","0163300000000000"); // no exotics
    cuts.AddCut("00051013","1111121053232200000","1111121053222700001","0163300000000000"); // frac = 0.99
    cuts.AddCut("00051013","1111121053332200000","1111121053322700001","0163300000000000"); // frac = 0.98
    cuts.AddCut("00051013","1111121053532200000","1111121053522700001","0163300000000000"); // frac = 0.97
    cuts.AddCut("00051013","1111121053732200000","1111121053722700001","0163300000000000"); // frac = 0.96
    cuts.AddCut("00051013","1111121053932200000","1111121053922700001","0163300000000000"); // frac = 0.94
  } else if (trainConfig == 25){  // variation of mass and alpha cut 
    cuts.AddCut("00051013","1111121053032200000","1111121053022700001","0163301700000000"); // only band at 0
    cuts.AddCut("00051013","1111121053032200000","1111121053022700001","0163300700000000"); // only band at 0, no alpha
    cuts.AddCut("00051013","1111121053032200000","1111121053022700001","0163301000000000"); // no mass cut
    cuts.AddCut("00051013","1111121053032200000","1111121053022700001","0163300000000000"); // no mass cut, no alpha
    cuts.AddCut("00051013","1111121053032200000","1111121053022700001","0163300100000000"); // no alpha
  } else if (trainConfig == 26){ // varied M02
    cuts.AddCut("00051013","1111121053032200000","1111121053022600001","0163300700000000"); // min M02= 0.3
    cuts.AddCut("00051013","1111121053032200000","1111121053022800001","0163300700000000"); // min M02= 0.25
    cuts.AddCut("00051013","1111121053032200000","1111121053022600001","0163300000000000"); // no mass, min M02 = 0.3
    cuts.AddCut("00051013","1111121053032200000","1111121053022700001","0163300000000000"); // no mass, min M02 = 0.27
    cuts.AddCut("00051013","1111121053032200000","1111121053022800001","0163300000000000"); // no mass, min M02 = 0.25
  } else if (trainConfig == 27){  // variation track matching to cluster & mass variations new defaults
    cuts.AddCut("00051013","1111121050032200000","1111121050022700001","0163300000000000"); // no TM
    cuts.AddCut("00051013","1111121051032200000","1111121051022700001","0163300000000000"); // looser TM
    cuts.AddCut("00051013","1111121056032200000","1111121056022700001","0163300000000000"); // tighter TM
    cuts.AddCut("00051013","1113121053032200000","1113121053022700001","0163300000000000");// TRD infront
    cuts.AddCut("00051013","1111221053032200000","1111221053022700001","0163300000000000");// no TRD infront
  } else if (trainConfig == 28){ // NL var new defaults
    cuts.AddCut("00051013","1111122053032200000","1111122053022700001","0163300000000000"); // SDM loose time
    cuts.AddCut("00051013","1111111053032200000","1111111053022700001","0163300000000000"); // conv calo tight time
    cuts.AddCut("00051013","1111101053032200000","1111101053022700001","0163300000000000"); // SDM Jason
    cuts.AddCut("00051013","1111100053032200000","1111100053022700001","0163300000000000"); // none
  
  } else if (trainConfig == 38){     // V1 clusterizer no NLM
    cuts.AddCut("00003113","1111121053032200000","1111121053022000000","0163300000000000"); // INT1
    cuts.AddCut("00003113","1111121053032200000","1111121053022700000","0163300000000000"); // INT1
    cuts.AddCut("00051013","1111121053032200000","1111121053022000000","0163300000000000"); // EMC1
    cuts.AddCut("00051013","1111121053032200000","1111121053022700000","0163300000000000"); // EMC1
  } else if (trainConfig == 39){     // V1 clusterizer NLM2
    cuts.AddCut("00003113","1111121053032200000","1111121053022000002","0163300000000000"); // INT1
    cuts.AddCut("00003113","1111121053032200000","1111121053022700002","0163300700000000"); // INT1
    cuts.AddCut("00003113","1111121053032200000","1111121053022210002","0163302200000000"); // INT1
    cuts.AddCut("00051013","1111121053032200000","1111121053022000002","0163300000000000"); // EMC1
    cuts.AddCut("00051013","1111121053032200000","1111121053022700002","0163300700000000"); // EMC1
    cuts.AddCut("00051013","1111121053032200000","1111121053022210002","0163302200000000"); // EMC1
    
  // LHC13g
  } else if (trainConfig == 40){  // new defaults LHC13g NLM1
    cuts.AddCut("00010113","1111121063032200000","1111121063022210001","0163301100000000"); // INT7
    cuts.AddCut("00052013","1111121063032200000","1111121063022210001","0163301100000000"); // EMC7
    cuts.AddCut("00085013","1111121063032200000","1111121063022210001","0163301100000000"); // EG2
    cuts.AddCut("00083013","1111121063032200000","1111121063022210001","0163301100000000"); // EG1
  } else if (trainConfig == 41){  // no TM in basis cuts LHC13g NLM1
    cuts.AddCut("00010113","1111121060032200000","1111121063022210001","0163301100000000"); // INT7
    cuts.AddCut("00052013","1111121060032200000","1111121063022210001","0163301100000000"); // EMC7
    cuts.AddCut("00085013","1111121060032200000","1111121063022210001","0163301100000000"); // EG2
    cuts.AddCut("00083013","1111121060032200000","1111121063022210001","0163301100000000"); // EG1
  } else if (trainConfig == 42){  // NLM1 no M02, no mass, no alpah
    cuts.AddCut("00010113","1111121063032200000","1111121063022000001","0163300000000000"); // INT7
    cuts.AddCut("00052013","1111121063032200000","1111121063022000001","0163300000000000"); // EMC7
    cuts.AddCut("00085013","1111121063032200000","1111121063022000001","0163300000000000"); // EG2
    cuts.AddCut("00083013","1111121063032200000","1111121063022000001","0163300000000000"); // EG1
  } else if (trainConfig == 43){  // new default
    cuts.AddCut("00010113","1111121063032200000","1111121063022700001","0163300700000000"); // Mass only band at 0, no explicit exotics cut, M02 cut at 0.27  
    cuts.AddCut("00052013","1111121063032200000","1111121063022700001","0163300700000000"); // Mass only band at 0, no explicit exotics cut, M02 cut at 0.27  
    cuts.AddCut("00085013","1111121063032200000","1111121063022700001","0163300700000000"); // Mass only band at 0, no explicit exotics cut, M02 cut at 0.27  
    cuts.AddCut("00083013","1111121063032200000","1111121063022700001","0163300700000000"); // Mass only band at 0, no explicit exotics cut, M02 cut at 0.27  
  } else if (trainConfig == 44){  // new default without mass
    cuts.AddCut("00010113","1111121063032200000","1111121063022700001","0163300000000000"); 
    cuts.AddCut("00052013","1111121063032200000","1111121063022700001","0163300000000000"); 
    cuts.AddCut("00085013","1111121063032200000","1111121063022700001","0163300000000000"); 
    cuts.AddCut("00083013","1111121063032200000","1111121063022700001","0163300000000000"); 
  } else if (trainConfig == 45){  //NLM exotics default frac = 0.97
    cuts.AddCut("00010113","1111121063532200000","1111121063522110001","0163301100000000"); // INT7
    cuts.AddCut("00052013","1111121063532200000","1111121063522110001","0163301100000000"); // EMC7
    cuts.AddCut("00085013","1111121063532200000","1111121063522110001","0163301100000000"); // EG2
    cuts.AddCut("00083013","1111121063532200000","1111121063522110001","0163301100000000"); // EG1
  } else if (trainConfig == 46){  // new default without mass, TM only in merged
    cuts.AddCut("00010113","1111121060032200000","1111121063022700001","0163300000000000"); 
    cuts.AddCut("00052013","1111121060032200000","1111121063022700001","0163300000000000"); 
    cuts.AddCut("00085013","1111121060032200000","1111121063022700001","0163300000000000"); 
    cuts.AddCut("00083013","1111121060032200000","1111121063022700001","0163300000000000"); 
  } else if (trainConfig == 47){  // NLM1 no M02, TM only in second
    cuts.AddCut("00010113","1111121060032200000","1111121063022000001","0163300000000000"); // INT7
    cuts.AddCut("00052013","1111121060032200000","1111121063022000001","0163300000000000"); // EMC7
    cuts.AddCut("00085013","1111121060032200000","1111121063022000001","0163300000000000"); // EG2
    cuts.AddCut("00083013","1111121060032200000","1111121063022000001","0163300000000000"); // EG1
  } else if (trainConfig == 48){  // new default without mass, pt dependent TM
    cuts.AddCut("00010113","1111121067032200000","1111121067022700001","0163300000000000"); 
    cuts.AddCut("00052013","1111121067032200000","1111121067022700001","0163300000000000"); 
    cuts.AddCut("00085013","1111121067032200000","1111121067022700001","0163300000000000"); 
    cuts.AddCut("00083013","1111121067032200000","1111121067022700001","0163300000000000"); 
  } else if (trainConfig == 49){  // new default without mass, TM only in merged, pt dependent TM
    cuts.AddCut("00010113","1111121060032200000","1111121067022700001","0163300000000000"); 
    cuts.AddCut("00052013","1111121060032200000","1111121067022700001","0163300000000000"); 
    cuts.AddCut("00085013","1111121060032200000","1111121067022700001","0163300000000000"); 
    cuts.AddCut("00083013","1111121060032200000","1111121067022700001","0163300000000000"); 

  // INT7 variations  
  } else if (trainConfig == 50){ // NLM 1 M02 var  
    cuts.AddCut("00010113","1111121063032200000","1111121063022110001","0163301100000000"); // min 0.3 function default
    cuts.AddCut("00010113","1111121063032200000","1111121063022310001","0163301100000000"); // min 0.25 function default
    cuts.AddCut("00010113","1111121063032200000","1111121063022410001","0163301100000000"); // min 0.27, tighter lower func
    cuts.AddCut("00010113","1111121063032200000","1111121063022520001","0163301100000000"); // min 0.27, looser func up and low
    cuts.AddCut("00010113","1111121063032200000","1111121063022430001","0163301100000000"); // min 0.27, tighter func up and low
  } else if (trainConfig == 51){  // EMCAL clusters, variation track matching to cluster & Mass INT7 NLM 1
    cuts.AddCut("00010113","1111121060032200000","1111121060022210001","0163301100000000"); // no TM
    cuts.AddCut("00010113","1111121061032200000","1111121061022210001","0163301100000000"); // looser TM
    cuts.AddCut("00010113","1111121066032200000","1111121066022210001","0163301100000000"); // tighter TM
    cuts.AddCut("00010113","1111121063032200000","1111121063022210001","0163301300000000"); // tighter mass cut
    cuts.AddCut("00010113","1111121063032200000","1111121063022210001","0163301500000000"); // looser mass cut
  } else if (trainConfig == 52){  // NL variations INT7
    cuts.AddCut("00010113","1111122063032200000","1111122063022210001","0163301100000000"); // SDM loose time
    cuts.AddCut("00010113","1111111063032200000","1111111063022210001","0163301100000000"); // conv calo tight time
    cuts.AddCut("00010113","1111101063032200000","1111101063022210001","0163301100000000"); // SDM Jason
    cuts.AddCut("00010113","1111100063032200000","1111100063022210001","0163301100000000"); // none
  } else if (trainConfig == 53){  // Alpha cut variations & TRD material INT7 NLM 1
    cuts.AddCut("00010113","1111121063032200000","1111121063022210001","0163303100000000"); // NLM 1 looser
    cuts.AddCut("00010113","1111121063032200000","1111121063022210001","0163305100000000"); // NLM 1 tighter
    cuts.AddCut("00010113","1112121063032200000","1112121063022210001","0163301100000000"); // TRD infront
    cuts.AddCut("00010113","1111321063032200000","1111321063022210001","0163301100000000"); // no TRD infront    
  } else if (trainConfig == 54){ // varied exoctics
    cuts.AddCut("00010113","1111121063032200000","1111121063022700001","0163300000000000"); // no exotics
    cuts.AddCut("00010113","1111121063232200000","1111121063222700001","0163300000000000"); // frac = 0.99
    cuts.AddCut("00010113","1111121063332200000","1111121063322700001","0163300000000000"); // frac = 0.98
    cuts.AddCut("00010113","1111121063532200000","1111121063522700001","0163300000000000"); // frac = 0.97
    cuts.AddCut("00010113","1111121063732200000","1111121063722700001","0163300000000000"); // frac = 0.96
    cuts.AddCut("00010113","1111121063932200000","1111121063922700001","0163300000000000"); // frac = 0.94
  } else if (trainConfig == 55){  // variation of mass and alpha cut 
    cuts.AddCut("00010113","1111121063032200000","1111121063022700001","0163301700000000"); // only band at 0
    cuts.AddCut("00010113","1111121063032200000","1111121063022700001","0163300700000000"); // only band at 0, no alpha
    cuts.AddCut("00010113","1111121063032200000","1111121063022700001","0163301000000000"); // no mass cut
    cuts.AddCut("00010113","1111121063032200000","1111121063022700001","0163300000000000"); // no mass cut, no alpha
    cuts.AddCut("00010113","1111121063032200000","1111121063022700001","0163300100000000"); // no alpha
  } else if (trainConfig == 56){ // varied M02
    cuts.AddCut("00010113","1111121063032200000","1111121063022600001","0163300700000000"); // min M02= 0.3
    cuts.AddCut("00010113","1111121063032200000","1111121063022800001","0163300700000000"); // min M02= 0.25
    cuts.AddCut("00010113","1111121063032200000","1111121063022600001","0163300000000000"); // no mass, min M02 = 0.3
    cuts.AddCut("00010113","1111121063032200000","1111121063022700001","0163300000000000"); // no mass, min M02 = 0.27
    cuts.AddCut("00010113","1111121063032200000","1111121063022800001","0163300000000000"); // no mass, min M02 = 0.25
  } else if (trainConfig == 57){  // variation track matching to cluster & mass variations new defaults
    cuts.AddCut("00010113","1111121060032200000","1111121060022700001","0163300000000000"); // no TM
    cuts.AddCut("00010113","1111121061032200000","1111121061022700001","0163300000000000"); // looser TM
    cuts.AddCut("00010113","1111121066032200000","1111121066022700001","0163300000000000"); // tighter TM
    cuts.AddCut("00010113","1112121063032200000","1112121063022700001","0163300000000000"); // TRD infront
    cuts.AddCut("00010113","1111321063032200000","1111321063022700001","0163300000000000"); // no TRD infront    
  } else if (trainConfig == 58){ // NL var new defaults
    cuts.AddCut("00010113","1111122063032200000","1111122063022700001","0163300000000000"); // SDM loose time
    cuts.AddCut("00010113","1111111063032200000","1111111063022700001","0163300000000000"); // conv calo tight time
    cuts.AddCut("00010113","1111101063032200000","1111101063022700001","0163300000000000"); // SDM Jason
    cuts.AddCut("00010113","1111100063032200000","1111100063022700001","0163300000000000"); // none
  
  // EMC7 variations  
  } else if (trainConfig == 60){ // NLM 1 M02 var  
    cuts.AddCut("00052013","1111121063032200000","1111121063022110001","0163301100000000"); // min 0.3 function default
    cuts.AddCut("00052013","1111121063032200000","1111121063022310001","0163301100000000"); // min 0.25 function default
    cuts.AddCut("00052013","1111121063032200000","1111121063022410001","0163301100000000"); // min 0.27, tighter lower func
    cuts.AddCut("00052013","1111121063032200000","1111121063022520001","0163301100000000"); // min 0.27, looser func up and low
    cuts.AddCut("00052013","1111121063032200000","1111121063022430001","0163301100000000"); // min 0.27, tighter func up and low
  } else if (trainConfig == 61){  // EMCAL clusters, variation track matching to cluster & Mass EMC7 NLM 1
    cuts.AddCut("00052013","1111121060032200000","1111121060022210001","0163301100000000"); // no TM
    cuts.AddCut("00052013","1111121061032200000","1111121061022210001","0163301100000000"); // looser TM
    cuts.AddCut("00052013","1111121066032200000","1111121066022210001","0163301100000000"); // tighter TM
    cuts.AddCut("00052013","1111121063032200000","1111121063022210001","0163301300000000"); // tighter mass cut
    cuts.AddCut("00052013","1111121063032200000","1111121063022210001","0163301500000000"); // looser mass cut
  } else if (trainConfig == 62){  // NL variations EMC7
    cuts.AddCut("00052013","1111122063032200000","1111122063022210001","0163301100000000"); // SDM loose time
    cuts.AddCut("00052013","1111111063032200000","1111111063022210001","0163301100000000"); // conv calo tight time
    cuts.AddCut("00052013","1111101063032200000","1111101063022210001","0163301100000000"); // SDM Jason
    cuts.AddCut("00052013","1111100063032200000","1111100063022210001","0163301100000000"); // none
  } else if (trainConfig == 63){  // Alpha cut variations & TRD material EMC7 NLM 1
    cuts.AddCut("00052013","1111121063032200000","1111121063022210001","0163303100000000"); // NLM 1 looser
    cuts.AddCut("00052013","1111121063032200000","1111121063022210001","0163305100000000"); // NLM 1 tighter
    cuts.AddCut("00052013","1112121063032200000","1112121063022210001","0163301100000000"); // TRD infront
    cuts.AddCut("00052013","1111321063032200000","1111321063022210001","0163301100000000"); // no TRD infront
  } else if (trainConfig == 64){ // varied exoctics
    cuts.AddCut("00052013","1111121063032200000","1111121063022700001","0163300000000000"); // no exotics
    cuts.AddCut("00052013","1111121063232200000","1111121063222700001","0163300000000000"); // frac = 0.99
    cuts.AddCut("00052013","1111121063332200000","1111121063322700001","0163300000000000"); // frac = 0.98
    cuts.AddCut("00052013","1111121063532200000","1111121063522700001","0163300000000000"); // frac = 0.97
    cuts.AddCut("00052013","1111121063732200000","1111121063722700001","0163300000000000"); // frac = 0.96
    cuts.AddCut("00052013","1111121063932200000","1111121063922700001","0163300000000000"); // frac = 0.94
  } else if (trainConfig == 65){  // variation of mass and alpha cut 
    cuts.AddCut("00052013","1111121063032200000","1111121063022700001","0163301700000000"); // only band at 0
    cuts.AddCut("00052013","1111121063032200000","1111121063022700001","0163300700000000"); // only band at 0, no alpha
    cuts.AddCut("00052013","1111121063032200000","1111121063022700001","0163301000000000"); // no mass cut
    cuts.AddCut("00052013","1111121063032200000","1111121063022700001","0163300000000000"); // no mass cut, no alpha
    cuts.AddCut("00052013","1111121063032200000","1111121063022700001","0163300100000000"); // no alpha
  } else if (trainConfig == 66){ // varied M02
    cuts.AddCut("00052013","1111121063032200000","1111121063022600001","0163300700000000"); // min M02= 0.3
    cuts.AddCut("00052013","1111121063032200000","1111121063022800001","0163300700000000"); // min M02= 0.25
    cuts.AddCut("00052013","1111121063032200000","1111121063022600001","0163300000000000"); // no mass, min M02 = 0.3
    cuts.AddCut("00052013","1111121063032200000","1111121063022700001","0163300000000000"); // no mass, min M02 = 0.27
    cuts.AddCut("00052013","1111121063032200000","1111121063022800001","0163300000000000"); // no mass, min M02 = 0.25
  } else if (trainConfig == 67){  // variation track matching to cluster & mass variations new defaults
    cuts.AddCut("00052013","1111121060032200000","1111121060022700001","0163300000000000"); // no TM
    cuts.AddCut("00052013","1111121061032200000","1111121061022700001","0163300000000000"); // looser TM
    cuts.AddCut("00052013","1111121066032200000","1111121066022700001","0163300000000000"); // tighter TM
    cuts.AddCut("00052013","1112121063032200000","1112121063022700001","0163300000000000"); // TRD infront
    cuts.AddCut("00052013","1111321063032200000","1111321063022700001","0163300000000000"); // no TRD infront    
  } else if (trainConfig == 68){ // NL var new defaults
    cuts.AddCut("00052013","1111122063032200000","1111122063022700001","0163300000000000"); // SDM loose time
    cuts.AddCut("00052013","1111111063032200000","1111111063022700001","0163300000000000"); // conv calo tight time
    cuts.AddCut("00052013","1111101063032200000","1111101063022700001","0163300000000000"); // SDM Jason
    cuts.AddCut("00052013","1111100063032200000","1111100063022700001","0163300000000000"); // none
    
  // EG2 variations
  } else if (trainConfig == 70){ // NLM 1 M02 var  
    cuts.AddCut("00085013","1111121063032200000","1111121063022110001","0163301100000000"); // min 0.3 function default
    cuts.AddCut("00085013","1111121063032200000","1111121063022310001","0163301100000000"); // min 0.25 function default
    cuts.AddCut("00085013","1111121063032200000","1111121063022410001","0163301100000000"); // min 0.27, tighter lower func
    cuts.AddCut("00085013","1111121063032200000","1111121063022520001","0163301100000000"); // min 0.27, looser func up and low
    cuts.AddCut("00085013","1111121063032200000","1111121063022430001","0163301100000000"); // min 0.27, tighter func up and low
  } else if (trainConfig == 71){  // EMCAL clusters, variation track matching to cluster & Mass EG2 NLM 1
    cuts.AddCut("00085013","1111121060032200000","1111121060022210001","0163301100000000"); // no TM
    cuts.AddCut("00085013","1111121061032200000","1111121061022210001","0163301100000000"); // looser TM
    cuts.AddCut("00085013","1111121066032200000","1111121066022210001","0163301100000000"); // tighter TM
    cuts.AddCut("00085013","1111121063032200000","1111121063022210001","0163301300000000"); // tighter mass cut
    cuts.AddCut("00085013","1111121063032200000","1111121063022210001","0163301500000000"); // looser mass cut
  } else if (trainConfig == 72){  // NL variations EG2
    cuts.AddCut("00085013","1111122063032200000","1111122063022210001","0163301100000000"); // SDM loose time
    cuts.AddCut("00085013","1111111063032200000","1111111063022210001","0163301100000000"); // conv calo tight time
    cuts.AddCut("00085013","1111101063032200000","1111101063022210001","0163301100000000"); // SDM Jason
    cuts.AddCut("00085013","1111100063032200000","1111100063022210001","0163301100000000"); // none
  } else if (trainConfig == 73){  // Alpha cut variations & TRD material EG2 NLM 1
    cuts.AddCut("00085013","1111121063032200000","1111121063022210001","0163303100000000"); // NLM 1 looser
    cuts.AddCut("00085013","1111121063032200000","1111121063022210001","0163305100000000"); // NLM 1 tighter
    cuts.AddCut("00085013","1112121063032200000","1112121063022210001","0163301100000000"); // TRD infront
    cuts.AddCut("00085013","1111321063032200000","1111321063022210001","0163301100000000"); // no TRD infront
  } else if (trainConfig == 74){ // varied exoctics
    cuts.AddCut("00085013","1111121063032200000","1111121063022700001","0163300000000000"); // no exotics
    cuts.AddCut("00085013","1111121063232200000","1111121063222700001","0163300000000000"); // frac = 0.99
    cuts.AddCut("00085013","1111121063332200000","1111121063322700001","0163300000000000"); // frac = 0.98
    cuts.AddCut("00085013","1111121063532200000","1111121063522700001","0163300000000000"); // frac = 0.97
    cuts.AddCut("00085013","1111121063732200000","1111121063722700001","0163300000000000"); // frac = 0.96
    cuts.AddCut("00085013","1111121063932200000","1111121063922700001","0163300000000000"); // frac = 0.94
  } else if (trainConfig == 75){  // variation of mass and alpha cut 
    cuts.AddCut("00085013","1111121063032200000","1111121063022700001","0163301700000000"); // only band at 0
    cuts.AddCut("00085013","1111121063032200000","1111121063022700001","0163300700000000"); // only band at 0, no alpha
    cuts.AddCut("00085013","1111121063032200000","1111121063022700001","0163301000000000"); // no mass cut
    cuts.AddCut("00085013","1111121063032200000","1111121063022700001","0163300000000000"); // no mass cut, no alpha
    cuts.AddCut("00085013","1111121063032200000","1111121063022700001","0163300100000000"); // no alpha
  } else if (trainConfig == 76){ // varied M02
    cuts.AddCut("00085013","1111121063032200000","1111121063022600001","0163300700000000"); // min M02= 0.3
    cuts.AddCut("00085013","1111121063032200000","1111121063022800001","0163300700000000"); // min M02= 0.25
    cuts.AddCut("00085013","1111121063032200000","1111121063022600001","0163300000000000"); // no mass, min M02 = 0.3
    cuts.AddCut("00085013","1111121063032200000","1111121063022700001","0163300000000000"); // no mass, min M02 = 0.27
    cuts.AddCut("00085013","1111121063032200000","1111121063022800001","0163300000000000"); // no mass, min M02 = 0.25
  } else if (trainConfig == 77){  // variation track matching to cluster & mass variations new defaults
    cuts.AddCut("00085013","1111121060032200000","1111121060022700001","0163300000000000"); // no TM
    cuts.AddCut("00085013","1111121061032200000","1111121061022700001","0163300000000000"); // looser TM
    cuts.AddCut("00085013","1111121066032200000","1111121066022700001","0163300000000000"); // tighter TM
    cuts.AddCut("00085013","1112121063032200000","1112121063022700001","0163300000000000"); // TRD infront
    cuts.AddCut("00085013","1111321063032200000","1111321063022700001","0163300000000000"); // no TRD infront    
  } else if (trainConfig == 78){ // NL var new defaults
    cuts.AddCut("00085013","1111122063032200000","1111122063022700001","0163300000000000"); // SDM loose time
    cuts.AddCut("00085013","1111111063032200000","1111111063022700001","0163300000000000"); // conv calo tight time
    cuts.AddCut("00085013","1111101063032200000","1111101063022700001","0163300000000000"); // SDM Jason
    cuts.AddCut("00085013","1111100063032200000","1111100063022700001","0163300000000000"); // none
        
  // EG1 variations   
  } else if (trainConfig == 80){ // NLM 1 M02 var  
    cuts.AddCut("00083013","1111121063032200000","1111121063022110001","0163301100000000"); // min 0.3 function default
    cuts.AddCut("00083013","1111121063032200000","1111121063022310001","0163301100000000"); // min 0.25 function default
    cuts.AddCut("00083013","1111121063032200000","1111121063022410001","0163301100000000"); // min 0.27, tighter lower func
    cuts.AddCut("00083013","1111121063032200000","1111121063022520001","0163301100000000"); // min 0.27, looser func up and low
    cuts.AddCut("00083013","1111121063032200000","1111121063022430001","0163301100000000"); // min 0.27, tighter func up and low
  } else if (trainConfig == 81){  // EMCAL clusters, variation track matching to cluster & Mass EG1 NLM 1
    cuts.AddCut("00083013","1111121060032200000","1111121060022210001","0163301100000000"); // no TM
    cuts.AddCut("00083013","1111121061032200000","1111121061022210001","0163301100000000"); // looser TM
    cuts.AddCut("00083013","1111121066032200000","1111121066022210001","0163301100000000"); // tighter TM
    cuts.AddCut("00083013","1111121063032200000","1111121063022210001","0163301300000000"); // tighter mass cut
    cuts.AddCut("00083013","1111121063032200000","1111121063022210001","0163301500000000"); // looser mass cut
  } else if (trainConfig == 82){  // NL variations EG1
    cuts.AddCut("00083013","1111122063032200000","1111122063022210001","0163301100000000"); // SDM loose time
    cuts.AddCut("00083013","1111111063032200000","1111111063022210001","0163301100000000"); // conv calo tight time
    cuts.AddCut("00083013","1111101063032200000","1111101063022210001","0163301100000000"); // SDM Jason
    cuts.AddCut("00083013","1111100063032200000","1111100063022210001","0163301100000000"); // none
  } else if (trainConfig == 83){  // Alpha cut variations & TRD material EG1 NLM 1
    cuts.AddCut("00083013","1111121063032200000","1111121063022210001","0163303100000000"); // NLM 1 looser
    cuts.AddCut("00083013","1111121063032200000","1111121063022210001","0163305100000000"); // NLM 1 tighter
    cuts.AddCut("00083013","1112121063032200000","1112121063022210001","0163301100000000"); // TRD infront
    cuts.AddCut("00083013","1111321063032200000","1111321063022210001","0163301100000000"); // no TRD infront
  } else if (trainConfig == 84){ // varied exoctics
    cuts.AddCut("00083013","1111121063032200000","1111121063022700001","0163300000000000"); // no exotics
    cuts.AddCut("00083013","1111121063232200000","1111121063222700001","0163300000000000"); // frac = 0.99
    cuts.AddCut("00083013","1111121063332200000","1111121063322700001","0163300000000000"); // frac = 0.98
    cuts.AddCut("00083013","1111121063532200000","1111121063522700001","0163300000000000"); // frac = 0.97
    cuts.AddCut("00083013","1111121063732200000","1111121063722700001","0163300000000000"); // frac = 0.96
    cuts.AddCut("00083013","1111121063932200000","1111121063922700001","0163300000000000"); // frac = 0.94
  } else if (trainConfig == 85){  // variation of mass and alpha cut 
    cuts.AddCut("00083013","1111121063032200000","1111121063022700001","0163301700000000"); // only band at 0
    cuts.AddCut("00083013","1111121063032200000","1111121063022700001","0163300700000000"); // only band at 0, no alpha
    cuts.AddCut("00083013","1111121063032200000","1111121063022700001","0163301000000000"); // no mass cut
    cuts.AddCut("00083013","1111121063032200000","1111121063022700001","0163300000000000"); // no mass cut, no alpha
    cuts.AddCut("00083013","1111121063032200000","1111121063022700001","0163300100000000"); // no alpha
  } else if (trainConfig == 86){ // varied M02
    cuts.AddCut("00083013","1111121063032200000","1111121063022600001","0163300700000000"); // min M02= 0.3
    cuts.AddCut("00083013","1111121063032200000","1111121063022800001","0163300700000000"); // min M02= 0.25
    cuts.AddCut("00083013","1111121063032200000","1111121063022600001","0163300000000000"); // no mass, min M02 = 0.3
    cuts.AddCut("00083013","1111121063032200000","1111121063022700001","0163300000000000"); // no mass, min M02 = 0.27
    cuts.AddCut("00083013","1111121063032200000","1111121063022800001","0163300000000000"); // no mass, min M02 = 0.25
  } else if (trainConfig == 87){  // variation track matching to cluster & mass variations new defaults
    cuts.AddCut("00083013","1111121060032200000","1111121060022700001","0163300000000000"); // no TM
    cuts.AddCut("00083013","1111121061032200000","1111121061022700001","0163300000000000"); // looser TM
    cuts.AddCut("00083013","1111121066032200000","1111121066022700001","0163300000000000"); // tighter TM
    cuts.AddCut("00083013","1112121063032200000","1112121063022700001","0163300000000000"); // TRD infront
    cuts.AddCut("00083013","1111321063032200000","1111321063022700001","0163300000000000"); // no TRD infront    
  } else if (trainConfig == 88){ // NL var new defaults
    cuts.AddCut("00083013","1111122063032200000","1111122063022700001","0163300000000000"); // SDM loose time
    cuts.AddCut("00083013","1111111063032200000","1111111063022700001","0163300000000000"); // conv calo tight time
    cuts.AddCut("00083013","1111101063032200000","1111101063022700001","0163300000000000"); // SDM Jason
    cuts.AddCut("00083013","1111100063032200000","1111100063022700001","0163300000000000"); // none
    
  } else if (trainConfig == 90){  // new default without mass, eta < 0.7, y < 0.7
    cuts.AddCut("00010113","1551121063032200000","1551121063022700001","0163200000000000"); 
    cuts.AddCut("00052013","1551121063032200000","1551121063022700001","0163200000000000"); 
    cuts.AddCut("00085013","1551121063032200000","1551121063022700001","0163200000000000"); 
    cuts.AddCut("00083013","1551121063032200000","1551121063022700001","0163200000000000"); 
  } else if (trainConfig == 91){  // new default without mass, TM pt dep, eta < 0.7, y < 0.7
    cuts.AddCut("00010113","1551121067032200000","1551121067022700001","0163200000000000"); 
    cuts.AddCut("00052013","1551121067032200000","1551121067022700001","0163200000000000"); 
    cuts.AddCut("00085013","1551121067032200000","1551121067022700001","0163200000000000"); 
    cuts.AddCut("00083013","1551121067032200000","1551121067022700001","0163200000000000"); 
  
    
  } else if (trainConfig == 95){  // new defaults LHC13g no NLM no mass, no alpha, no M02
    cuts.AddCut("00010113","1111121063032200000","1111121063022000000","0163300000000000"); // INT7
    cuts.AddCut("00052013","1111121063032200000","1111121063022000000","0163300000000000"); // EMC7
    cuts.AddCut("00085013","1111121063032200000","1111121063022000000","0163300000000000"); // EG2
    cuts.AddCut("00083013","1111121063032200000","1111121063022000000","0163300000000000"); // EG1    
  } else if (trainConfig == 96){  // new defaults LHC13g no NLM, no mass
    cuts.AddCut("00010113","1111121063032200000","1111121063022700000","0163300000000000"); // INT7
    cuts.AddCut("00052013","1111121063032200000","1111121063022700000","0163300000000000"); // EMC7
    cuts.AddCut("00085013","1111121063032200000","1111121063022700000","0163300000000000"); // EG2
    cuts.AddCut("00083013","1111121063032200000","1111121063022700000","0163300000000000"); // EG1    
  } else if (trainConfig == 97){  // new defaults LHC13g NLM2
    cuts.AddCut("00010113","1111121063032200000","1111121063022700002","0163300700000000"); // INT7
    cuts.AddCut("00052013","1111121063032200000","1111121063022700002","0163300700000000"); // EMC7
    cuts.AddCut("00085013","1111121063032200000","1111121063022700002","0163300700000000"); // EG2
    cuts.AddCut("00083013","1111121063032200000","1111121063022700002","0163300700000000"); // EG1    
  } else if (trainConfig == 98){  // new defaults LHC13g NLM2
    cuts.AddCut("00010113","1111121063032200000","1111121063022210002","0163302200000000"); // INT7
    cuts.AddCut("00052013","1111121063032200000","1111121063022210002","0163302200000000"); // EMC7
    cuts.AddCut("00085013","1111121063032200000","1111121063022210002","0163302200000000"); // EG2
    cuts.AddCut("00083013","1111121063032200000","1111121063022210002","0163302200000000"); // EG1
  } else if (trainConfig == 99){  // NLM2 no M02, no mass, no alpah
    cuts.AddCut("00010113","1111121063032200000","1111121063022000002","0163300000000000"); // INT7
    cuts.AddCut("00052013","1111121063032200000","1111121063022000002","0163300000000000"); // EMC7
    cuts.AddCut("00085013","1111121063032200000","1111121063022000002","0163300000000000"); // EG2
    cuts.AddCut("00083013","1111121063032200000","1111121063022000002","0163300000000000"); // EG1
    
  // LHC12
  // default
  } else if (trainConfig == 107){  // no M02, pt dep TM
    cuts.AddCut("00010113","1111111067032200000","1111111067022000001","0163300000000000"); // INT7
    cuts.AddCut("00052113","1111111067032200000","1111111067022000001","0163300000000000"); // EMC7
    cuts.AddCut("00081113","1111111067032200000","1111111067022000001","0163300000000000"); // EGA
  } else if (trainConfig == 108){  // no M02, TM only in merged, pt dep TM
    cuts.AddCut("00010113","1111111060032200000","1111111067022000001","0163300000000000"); // INT7
    cuts.AddCut("00052113","1111111060032200000","1111111067022000001","0163300000000000"); // EMC7
    cuts.AddCut("00081113","1111111060032200000","1111111067022000001","0163300000000000"); // EGA
  } else if (trainConfig == 109){  // M02 cut at 0.27, pt dep TM
    cuts.AddCut("00010113","1111111067032200000","1111111067022700001","0163300000000000"); // INT7
    cuts.AddCut("00052113","1111111067032200000","1111111067022700001","0163300000000000"); // EMC7
    cuts.AddCut("00081113","1111111067032200000","1111111067022700001","0163300000000000"); // EGA
  } else if (trainConfig == 110){  // M02 cut at 0.27, TM only in merged, pt dep TM
    cuts.AddCut("00010113","1111111060032200000","1111111067022700001","0163300000000000"); // INT7
    cuts.AddCut("00052113","1111111060032200000","1111111067022700001","0163300000000000"); // EMC7
    cuts.AddCut("00081113","1111111060032200000","1111111067022700001","0163300000000000"); // EGA

  } else if (trainConfig == 111){  // EMCAL clusters, different triggers no NonLinearity
    cuts.AddCut("00010113","1111100060032200000","1111100060022210001","0163301100000000"); // INT7
    cuts.AddCut("00052113","1111100060032200000","1111100060022210001","0163301100000000"); // EMC7
    cuts.AddCut("00081113","1111100060032200000","1111100060022110001","0163301100000000"); // EMCEGA,

  } else if (trainConfig == 112){  // EMCAL clusters, different triggers with NonLinearity
    cuts.AddCut("00010113","1111111060032200000","1111111060022210001","0163301100000000"); // INT7
    cuts.AddCut("00052113","1111111060032200000","1111111060022210001","0163301100000000"); // EMC7
    cuts.AddCut("00081113","1111111060032200000","1111111060022210001","0163301100000000"); // EMCEGA,
  } else if (trainConfig == 113){  // EMCAL clusters, different triggers with NonLinearity track matching to cluster
    cuts.AddCut("00010113","1111111067032200000","1111111067022210001","0163301100000000"); // INT7
    cuts.AddCut("00052113","1111111067032200000","1111111067022210001","0163301100000000"); // EMC7
    cuts.AddCut("00081113","1111111067032200000","1111111067022210001","0163301100000000"); // EMCEGA,

  // new default
  } else if (trainConfig == 114){  // NLM1 no M02, no mass, no alpah
    cuts.AddCut("00010113","1111111067032200000","1111111067022000001","0163300000000000"); // INT7
    cuts.AddCut("00052113","1111111067032200000","1111111067022000001","0163300000000000"); // EMC7
    cuts.AddCut("00081113","1111111067032200000","1111111067022000001","0163300000000000"); // EGA
  } else if (trainConfig == 115){  // Mass only band at 0, M02 cut at 0.27  
    cuts.AddCut("00010113","1111111067032200000","1111111067022700001","0163300700000000"); // INT7
    cuts.AddCut("00052113","1111111067032200000","1111111067022700001","0163300700000000"); // EMC7
    cuts.AddCut("00081113","1111111067032200000","1111111067022700001","0163300700000000"); // EGA
  } else if (trainConfig == 116){  // M02 cut at 0.27  
    cuts.AddCut("00010113","1111111067032200000","1111111067022700001","0163300000000000"); // INT7
    cuts.AddCut("00052113","1111111067032200000","1111111067022700001","0163300000000000"); // EMC7
    cuts.AddCut("00081113","1111111067032200000","1111111067022700001","0163300000000000"); // EGA
  } else if (trainConfig == 117){  // M02 cut at 0.27, TM only in merged
    cuts.AddCut("00010113","1111111060032200000","1111111067022700001","0163300000000000"); // INT7
    cuts.AddCut("00052113","1111111060032200000","1111111067022700001","0163300000000000"); // EMC7
    cuts.AddCut("00081113","1111111060032200000","1111111067022700001","0163300000000000"); // EGA
  } else if (trainConfig == 118){  // NLM1 no M02, TM only in merged
    cuts.AddCut("00010113","1111111060032200000","1111111067022000001","0163300000000000"); // INT7
    cuts.AddCut("00052113","1111111060032200000","1111111067022000001","0163300000000000"); // EMC7
    cuts.AddCut("00081113","1111111060032200000","1111111067022000001","0163300000000000"); // EGA
  } else if (trainConfig == 119){  // new default, with eta < 0.7, y < 0.7
    cuts.AddCut("00010113","1551111067032200000","1551111067022700001","0163200000000000"); // INT7
    cuts.AddCut("00052113","1551111067032200000","1551111067022700001","0163200000000000"); // EMC7
    cuts.AddCut("00081113","1551111067032200000","1551111067022700001","0163200000000000"); // EGA
    
    //kINT7 - NLM1
  } else if (trainConfig == 120){  // TRD material INT7 NLM 1
    cuts.AddCut("00010113","1112111067032200000","1112111067022700001","0163300000000000"); // TRD infront
    cuts.AddCut("00010113","1111311067032200000","1111311067022700001","0163300000000000"); // no TRD infront
  } else if (trainConfig == 121){ // varied exoctics
    cuts.AddCut("00010113","1111111067032200000","1111111067022700001","0163300000000000"); // no exotics
    cuts.AddCut("00010113","1111111067232200000","1111111067222700001","0163300000000000"); // frac = 0.99
    cuts.AddCut("00010113","1111111067332200000","1111111067322700001","0163300000000000"); // frac = 0.98
    cuts.AddCut("00010113","1111111067532200000","1111111067522700001","0163300000000000"); // frac = 0.97
    cuts.AddCut("00010113","1111111067732200000","1111111067722700001","0163300000000000"); // frac = 0.96
    cuts.AddCut("00010113","1111111067932200000","1111111067922700001","0163300000000000"); // frac = 0.94
  } else if (trainConfig == 122){  // variation of mass and alpha cut 
    cuts.AddCut("00010113","1111111067032200000","1111111067022700001","0163301700000000"); // only band at 0
    cuts.AddCut("00010113","1111111067032200000","1111111067022700001","0163300700000000"); // only band at 0, no alpha
    cuts.AddCut("00010113","1111111067032200000","1111111067022700001","0163301000000000"); // no mass cut
    cuts.AddCut("00010113","1111111067032200000","1111111067022700001","0163300000000000"); // no mass cut, no alpha
    cuts.AddCut("00010113","1111111067032200000","1111111067022700001","0163300100000000"); // no alpha
  } else if (trainConfig == 123){ // varied M02
    cuts.AddCut("00010113","1111111067032200000","1111111067022600001","0163300700000000"); // min M02= 0.3
    cuts.AddCut("00010113","1111111067032200000","1111111067022800001","0163300700000000"); // min M02= 0.25
    cuts.AddCut("00010113","1111111067032200000","1111111067022600001","0163300000000000"); // no mass, min M02 = 0.3
    cuts.AddCut("00010113","1111111067032200000","1111111067022700001","0163300000000000"); // no mass, min M02 = 0.27
    cuts.AddCut("00010113","1111111067032200000","1111111067022800001","0163300000000000"); // no mass, min M02 = 0.25
  } else if (trainConfig == 124){  // variation track matching to cluster & mass variations new defaults
    cuts.AddCut("00010113","1111111060032200000","1111121060022700001","0163300000000000"); // no TM
    cuts.AddCut("00010113","1111111063032200000","1111121063022700001","0163300000000000"); // old matching std
    cuts.AddCut("00010113","1111111068032200000","1111121068022700001","0163300000000000"); // looser TM
    cuts.AddCut("00010113","1111111066032200000","1111121066022700001","0163300000000000"); // tighter TM
    cuts.AddCut("00010113","1111111067032200000","1111111067022700001","0163300800000000"); // looser mass cut
    cuts.AddCut("00010113","1111111067032200000","1111111067022700001","0163300900000000"); // tighter mass cut
  } else if (trainConfig == 125){ // NL var new defaults
    cuts.AddCut("00010113","1111112067032200000","1111112067022700001","0163300000000000"); // SDM loose time
    cuts.AddCut("00010113","1111101067032200000","1111101067022700001","0163300000000000"); // SDM Jason
    cuts.AddCut("00010113","1111100067032200000","1111100067022700001","0163300000000000"); // none

    //kEMC7 - NLM1
  } else if (trainConfig == 130){  // TRD material INT7 NLM 1
    cuts.AddCut("00052113","1112111067032200000","1112111067022700001","0163300000000000"); // TRD infront
    cuts.AddCut("00052113","1111311067032200000","1111311067022700001","0163300000000000"); // no TRD infront
  } else if (trainConfig == 131){ // varied exoctics
    cuts.AddCut("00052113","1111111067032200000","1111111067022700001","0163300000000000"); // no exotics
    cuts.AddCut("00052113","1111111067232200000","1111111067222700001","0163300000000000"); // frac = 0.99
    cuts.AddCut("00052113","1111111067332200000","1111111067322700001","0163300000000000"); // frac = 0.98
    cuts.AddCut("00052113","1111111067532200000","1111111067522700001","0163300000000000"); // frac = 0.97
    cuts.AddCut("00052113","1111111067732200000","1111111067722700001","0163300000000000"); // frac = 0.96
    cuts.AddCut("00052113","1111111067932200000","1111111067922700001","0163300000000000"); // frac = 0.94
  } else if (trainConfig == 132){  // variation of mass and alpha cut 
    cuts.AddCut("00052113","1111111067032200000","1111111067022700001","0163301700000000"); // only band at 0
    cuts.AddCut("00052113","1111111067032200000","1111111067022700001","0163300700000000"); // only band at 0, no alpha
    cuts.AddCut("00052113","1111111067032200000","1111111067022700001","0163301000000000"); // no mass cut
    cuts.AddCut("00052113","1111111067032200000","1111111067022700001","0163300000000000"); // no mass cut, no alpha
    cuts.AddCut("00052113","1111111067032200000","1111111067022700001","0163300100000000"); // no alpha
  } else if (trainConfig == 133){ // varied M02
    cuts.AddCut("00052113","1111111067032200000","1111111067022600001","0163300700000000"); // min M02= 0.3
    cuts.AddCut("00052113","1111111067032200000","1111111067022800001","0163300700000000"); // min M02= 0.25
    cuts.AddCut("00052113","1111111067032200000","1111111067022600001","0163300000000000"); // no mass, min M02 = 0.3
    cuts.AddCut("00052113","1111111067032200000","1111111067022700001","0163300000000000"); // no mass, min M02 = 0.27
    cuts.AddCut("00052113","1111111067032200000","1111111067022800001","0163300000000000"); // no mass, min M02 = 0.25
  } else if (trainConfig == 134){  // variation track matching to cluster & mass variations new defaults
    cuts.AddCut("00052113","1111111060032200000","1111121060022700001","0163300000000000"); // no TM
    cuts.AddCut("00052113","1111111063032200000","1111121063022700001","0163300000000000"); // old matching std
    cuts.AddCut("00052113","1111111068032200000","1111121068022700001","0163300000000000"); // looser TM
    cuts.AddCut("00052113","1111111066032200000","1111121066022700001","0163300000000000"); // tighter TM
    cuts.AddCut("00052113","1111111067032200000","1111111067022700001","0163300800000000"); // looser mass cut
    cuts.AddCut("00052113","1111111067032200000","1111111067022700001","0163300900000000"); // tighter mass cut
  } else if (trainConfig == 135){ // NL var new defaults
    cuts.AddCut("00052113","1111112067032200000","1111112067022700001","0163300000000000"); // SDM loose time
    cuts.AddCut("00052113","1111101067032200000","1111101067022700001","0163300000000000"); // SDM Jason
    cuts.AddCut("00052113","1111100067032200000","1111100067022700001","0163300000000000"); // none
  
  //kEMCEGA - NLM1
  } else if (trainConfig == 140){  // TRD material INT7 NLM 1
    cuts.AddCut("00081113","1112111067032200000","1112111067022700001","0163300000000000"); // TRD infront
    cuts.AddCut("00081113","1111311067032200000","1111311067022700001","0163300000000000"); // no TRD infront
  } else if (trainConfig == 141){ // varied exoctics
    cuts.AddCut("00081113","1111111067032200000","1111111067022700001","0163300000000000"); // no exotics
    cuts.AddCut("00081113","1111111067232200000","1111111067222700001","0163300000000000"); // frac = 0.99
    cuts.AddCut("00081113","1111111067332200000","1111111067322700001","0163300000000000"); // frac = 0.98
    cuts.AddCut("00081113","1111111067532200000","1111111067522700001","0163300000000000"); // frac = 0.97
    cuts.AddCut("00081113","1111111067732200000","1111111067722700001","0163300000000000"); // frac = 0.96
    cuts.AddCut("00081113","1111111067932200000","1111111067922700001","0163300000000000"); // frac = 0.94
  } else if (trainConfig == 142){  // variation of mass and alpha cut 
    cuts.AddCut("00081113","1111111067032200000","1111111067022700001","0163301700000000"); // only band at 0
    cuts.AddCut("00081113","1111111067032200000","1111111067022700001","0163300700000000"); // only band at 0, no alpha
    cuts.AddCut("00081113","1111111067032200000","1111111067022700001","0163301000000000"); // no mass cut
    cuts.AddCut("00081113","1111111067032200000","1111111067022700001","0163300000000000"); // no mass cut, no alpha
    cuts.AddCut("00081113","1111111067032200000","1111111067022700001","0163300100000000"); // no alpha
  } else if (trainConfig == 143){ // varied M02
    cuts.AddCut("00081113","1111111067032200000","1111111067022600001","0163300700000000"); // min M02= 0.3
    cuts.AddCut("00081113","1111111067032200000","1111111067022800001","0163300700000000"); // min M02= 0.25
    cuts.AddCut("00081113","1111111067032200000","1111111067022600001","0163300000000000"); // no mass, min M02 = 0.3
    cuts.AddCut("00081113","1111111067032200000","1111111067022700001","0163300000000000"); // no mass, min M02 = 0.27
    cuts.AddCut("00081113","1111111067032200000","1111111067022800001","0163300000000000"); // no mass, min M02 = 0.25
  } else if (trainConfig == 144){  // variation track matching to cluster & mass variations new defaults
    cuts.AddCut("00081113","1111111060032200000","1111121060022700001","0163300000000000"); // no TM
    cuts.AddCut("00081113","1111111063032200000","1111121063022700001","0163300000000000"); // old matching std
    cuts.AddCut("00081113","1111111068032200000","1111121068022700001","0163300000000000"); // looser TM
    cuts.AddCut("00081113","1111111066032200000","1111121066022700001","0163300000000000"); // tighter TM
    cuts.AddCut("00081113","1111111067032200000","1111111067022700001","0163300800000000"); // looser mass cut
    cuts.AddCut("00081113","1111111067032200000","1111111067022700001","0163300900000000"); // tighter mass cut
  } else if (trainConfig == 145){ // NL var new defaults
    cuts.AddCut("00081113","1111112067032200000","1111112067022700001","0163300000000000"); // SDM loose time
    cuts.AddCut("00081113","1111101067032200000","1111101067022700001","0163300000000000"); // SDM Jason
    cuts.AddCut("00081113","1111100067032200000","1111100067022700001","0163300000000000"); // none
    
  // NLM2 cuts  
  } else if (trainConfig == 160){  //no NLM no explicit exotics cut, no M02 cut, no mass
    cuts.AddCut("00010113","1111111067032200000","1111111067022000000","0163300000000000"); // INT7
    cuts.AddCut("00052113","1111111067032200000","1111111067022000000","0163300000000000"); // EMC7
    cuts.AddCut("00081113","1111111067032200000","1111111067022000000","0163300000000000"); // EGA
  } else if (trainConfig == 161){  // no NLM no mass, no explicit exotics cut, M02 cut at 0.27  
    cuts.AddCut("00010113","1111111067032200000","1111111067022700000","0163300000000000"); // INT7
    cuts.AddCut("00052113","1111111067032200000","1111111067022700000","0163300000000000"); // EMC7
    cuts.AddCut("00081113","1111111067032200000","1111111067022700000","0163300000000000"); // EGA

  } else if (trainConfig == 162){  // NLM2 no M02, no mass, no alpha
    cuts.AddCut("00010113","1111111067032200000","1111111067022000002","0163300000000000"); // INT7
    cuts.AddCut("00052113","1111111067032200000","1111111067022000002","0163300000000000"); // EMC7
    cuts.AddCut("00081113","1111111067032200000","1111111067022000002","0163300000000000"); // EGA
  } else if (trainConfig == 163){  // NLM2 Mass only band at 0, no explicit exotics cut, M02 cut at 0.27  
    cuts.AddCut("00010113","1111111067032200000","1111111067022700002","0163300700000000"); // INT7
    cuts.AddCut("00052113","1111111067032200000","1111111067022700002","0163300700000000"); // EMC7
    cuts.AddCut("00081113","1111111067032200000","1111111067022700002","0163300700000000"); // EGA
  } else if (trainConfig == 164){  // EMCAL clusters, different triggers with NonLinearity track matching to cluster
    cuts.AddCut("00010113","1111111067032200000","1111111067022210002","0163302200000000"); // INT7
    cuts.AddCut("00052113","1111111067032200000","1111111067022210002","0163302200000000"); // EMC7
    cuts.AddCut("00081113","1111111067032200000","1111111067022210002","0163302200000000"); // EMCEGA,

  // diff eta/rap cuts
  } else if (trainConfig == 171){  // std pT dep, |eta| < 0.3, y < 0.3
    cuts.AddCut("00010113","1661111067032200000","1661111067022700001","0163700000000000"); // INT7
    cuts.AddCut("00052113","1661111067032200000","1661111067022700001","0163700000000000"); // EMC7
    cuts.AddCut("00081113","1661111067032200000","1661111067022700001","0163700000000000"); // EGA
  } else if (trainConfig == 172){  // std, |eta| < 0.3, y < 0.3
    cuts.AddCut("00010113","1661111063032200000","1661111063022700001","0163700000000000"); // INT7
    cuts.AddCut("00052113","1661111063032200000","1661111063022700001","0163700000000000"); // EMC7
    cuts.AddCut("00081113","1661111063032200000","1661111063022700001","0163700000000000"); // EGA
  } else if (trainConfig == 173){  // no M02 pT dep, |eta| < 0.3, y < 0.3
    cuts.AddCut("00010113","1661111067032200000","1661111067022000001","0163700000000000"); // INT7
    cuts.AddCut("00052113","1661111067032200000","1661111067022000001","0163700000000000"); // EMC7
    cuts.AddCut("00081113","1661111067032200000","1661111067022000001","0163700000000000"); // EGA
  } else if (trainConfig == 174){  // no M02, |eta| < 0.3, y < 0.3
    cuts.AddCut("00010113","1661111063032200000","1661111063022000001","0163700000000000"); // INT7
    cuts.AddCut("00052113","1661111063032200000","1661111063022000001","0163700000000000"); // EMC7
    cuts.AddCut("00081113","1661111063032200000","1661111063022000001","0163700000000000"); // EGA

  // T0AND cuts
  } else if (trainConfig == 181){  // EMCAL clusters, different triggers with NonLinearity track matching to cluster
    cuts.AddCut("00011113","1111111067032200000","1111111067022110001","0163301100000000"); // INT8
    cuts.AddCut("00053113","1111111067032200000","1111111067022110001","0163301100000000"); // EMC8
    cuts.AddCut("00082113","1111111067032200000","1111111067022110001","0163301100000000"); // EMC8EGA,
  } else if (trainConfig == 182){  // EMCAL clusters, different triggers with NonLinearity track matching to cluster
    cuts.AddCut("00011113","1111111067032200000","1111111067022110002","0163302200000000"); // INT8
    cuts.AddCut("00053113","1111111067032200000","1111111067022110002","0163302200000000"); // EMC8
    cuts.AddCut("00082113","1111111067032200000","1111111067022110002","0163302200000000"); // EMC8EGA,

    // 13 TeV & 5 TeV
  } else if (trainConfig == 401){ // EMCAL clusters pp 13 TeV
    cuts.AddCut("00010113","1111101013032200000","1111101013022210002","0163302200000000"); //INT7 1000ns timing cut, std NL NLM2
    cuts.AddCut("00052013","1111101013032200000","1111101013022210002","0163302200000000"); //INT7 1000ns timing cut, std NL NLM2
    cuts.AddCut("00010113","1111101013032200000","1111101013022110001","0163301100000000"); //INT7 1000ns timing cut, std NL NLM1
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
    //create AliCaloTrackMatcher instance, if there is none present
    TString caloCutPos = cuts.GetClusterCut(i);
    caloCutPos.Resize(1);
    TString TrackMatcherName = Form("CaloTrackMatcher_%s",caloCutPos.Data());
    if( !(AliCaloTrackMatcher*)mgr->GetTask(TrackMatcherName.Data()) ){
      AliCaloTrackMatcher* fTrackMatcher = new AliCaloTrackMatcher(TrackMatcherName.Data(),caloCutPos.Atoi());
      fTrackMatcher->SetV0ReaderName(V0ReaderName);
      mgr->AddTask(fTrackMatcher);
      mgr->ConnectInput(fTrackMatcher,0,cinput);
    }

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
    if(periodNameV0Reader.CompareTo("") != 0) analysisEventCuts[i]->SetPeriodEnum(periodNameV0Reader);
    analysisEventCuts[i]->SetLightOutput(runLightOutput);
    analysisEventCuts[i]->InitializeCutsFromCutString((cuts.GetEventCut(i)).Data());
    EventCutList->Add(analysisEventCuts[i]);
    analysisEventCuts[i]->SetFillCutHistograms("",kFALSE);
    
    analysisClusterCuts[i]        = new AliCaloPhotonCuts(isMC);
    analysisClusterCuts[i]->SetIsPureCaloCut(2);
    analysisClusterCuts[i]->SetHistoToModifyAcceptance(histoAcc);
    analysisClusterCuts[i]->SetV0ReaderName(V0ReaderName);
    analysisClusterCuts[i]->SetCaloTrackMatcherName(TrackMatcherName);
    analysisClusterCuts[i]->SetLightOutput(runLightOutput);
    analysisClusterCuts[i]->SetExoticsMinCellEnergyCut(minEnergyForExoticsCut);
    analysisClusterCuts[i]->SetExoticsQA(runQAForExotics);
    analysisClusterCuts[i]->InitializeCutsFromCutString((cuts.GetClusterCut(i)).Data());
    ClusterCutList->Add(analysisClusterCuts[i]);
    analysisClusterCuts[i]->SetExtendedMatchAndQA(enableExtMatchAndQA);
    analysisClusterCuts[i]->SetFillCutHistograms("");
    
    analysisClusterMergedCuts[i]  = new AliCaloPhotonCuts(isMC);
    analysisClusterMergedCuts[i]->SetIsPureCaloCut(1);
    analysisClusterMergedCuts[i]->SetHistoToModifyAcceptance(histoAcc);
    analysisClusterMergedCuts[i]->SetV0ReaderName(V0ReaderName);
    analysisClusterMergedCuts[i]->SetCaloTrackMatcherName(TrackMatcherName);
    analysisClusterMergedCuts[i]->SetLightOutput(runLightOutput);
    analysisClusterMergedCuts[i]->SetExoticsMinCellEnergyCut(minEnergyForExoticsCut);
//     analysisClusterMergedCuts[i]->SetExoticsQA(runQAForExotics);
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
  task->SetSelectedMesonID(selectedMeson);
  task->SetEventCutList(numberOfCuts,EventCutList);
  task->SetCaloCutList(numberOfCuts,ClusterCutList);
  task->SetCaloMergedCutList(numberOfCuts,ClusterMergedCutList);
  task->SetMesonCutList(numberOfCuts,MesonCutList);
  task->SetDoMesonQA(enableQAMesonTask); //Attention new switch for Pi0 QA
  task->SetDoClusterQA(enableQAClusterTask);//Attention new switch small for Cluster QA
  task->SetEnableSortingOfMCClusLabels(enableSortingMCLabels);
  if(enableExtMatchAndQA > 1){ task->SetPlotHistsExtQA(kTRUE);}
  if (enableDetailedPrintout) task->SetEnableDetailedPrintout(enableDetailedPrintout);//Attention new switch small for Cluster QA
  //connect containers
  AliAnalysisDataContainer *coutput =
    mgr->CreateContainer(Form("GammaCaloMerged_%i",trainConfig), TList::Class(),
              AliAnalysisManager::kOutputContainer,Form("GammaCaloMerged_%i.root",trainConfig));

  mgr->AddTask(task);
  mgr->ConnectInput(task, 0, cinput);
  mgr->ConnectOutput(task, 1, coutput);

  return;
}
