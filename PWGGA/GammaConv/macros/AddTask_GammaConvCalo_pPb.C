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
//($ALIPHYSICS/PWGGA/GammaConv/AliAnalysisTaskGammaConvCalo.cxx) for
//pPb together with all supporting classes
//***************************************************************************************

//***************************************************************************************
//CutHandler contains all cuts for a certain analysis and trainconfig,
//it automatically checks length of cutStrings and takes care of the number of added cuts,
//no specification of the variable 'numberOfCuts' needed anymore.
//***************************************************************************************
class CutHandlerConvCalo{
  public:
    CutHandlerConvCalo(Int_t nMax=10){
      nCuts=0; nMaxCuts=nMax; validCuts = true;
      eventCutArray = new TString[nMaxCuts]; photonCutArray = new TString[nMaxCuts]; clusterCutArray = new TString[nMaxCuts]; mesonCutArray = new TString[nMaxCuts];
      for(Int_t i=0; i<nMaxCuts; i++) {eventCutArray[i] = ""; photonCutArray[i] = ""; clusterCutArray[i] = ""; mesonCutArray[i] = "";}
    }

    void AddCut(TString eventCut, TString photonCut, TString clusterCut, TString mesonCut){
      if(nCuts>=nMaxCuts) {cout << "ERROR in CutHandlerConvCalo: Exceeded maximum number of cuts!" << endl; validCuts = false; return;}
      if( eventCut.Length()!=8 || photonCut.Length()!=26 || clusterCut.Length()!=19 || mesonCut.Length()!=16 ) {cout << "ERROR in CutHandlerConvCalo: Incorrect length of cut string!" << endl; validCuts = false; return;}
      eventCutArray[nCuts]=eventCut; photonCutArray[nCuts]=photonCut; clusterCutArray[nCuts]=clusterCut; mesonCutArray[nCuts]=mesonCut;
      nCuts++;
      return;
    }
    Bool_t AreValid(){return validCuts;}
    Int_t GetNCuts(){if(validCuts) return nCuts; else return 0;}
    TString GetEventCut(Int_t i){if(validCuts&&i<nMaxCuts&&i>=0) return eventCutArray[i]; else{cout << "ERROR in CutHandlerConvCalo: GetEventCut wrong index i" << endl;return "";}}
    TString GetPhotonCut(Int_t i){if(validCuts&&i<nMaxCuts&&i>=0) return photonCutArray[i]; else {cout << "ERROR in CutHandlerConvCalo: GetPhotonCut wrong index i" << endl;return "";}}
    TString GetClusterCut(Int_t i){if(validCuts&&i<nMaxCuts&&i>=0) return clusterCutArray[i]; else {cout << "ERROR in CutHandlerConvCalo: GetClusterCut wrong index i" << endl;return "";}}
    TString GetMesonCut(Int_t i){if(validCuts&&i<nMaxCuts&&i>=0) return mesonCutArray[i]; else {cout << "ERROR in CutHandlerConvCalo: GetMesonCut wrong index i" << endl;return "";}}
  private:
    Bool_t validCuts;
    Int_t nCuts; Int_t nMaxCuts;
    TString* eventCutArray;
    TString* photonCutArray;
    TString* clusterCutArray;
    TString* mesonCutArray;
};

//***************************************************************************************
//main function
//***************************************************************************************
void AddTask_GammaConvCalo_pPb( Int_t     trainConfig                 = 1,                  // change different set of cuts
                                Int_t     isMC                        = 0,               // run MC
                                Int_t     enableQAMesonTask           = 0,                 // enable QA in AliAnalysisTaskGammaConvV1
                                Int_t     enableQAPhotonTask          = 0,                 // enable additional QA task
                                TString   fileNameInputForWeighting   = "MCSpectraInput.root",       // path to file for weigting input
                                Int_t     doWeightingPart             = 0,                  // enable Weighting
                                TString   generatorName               = "DPMJET",              // generator Name  
                                TString   cutnumberAODBranch          = "800000006008400000001500000",  // cutnumber for AOD branch
                                Int_t     enableExtMatchAndQA         = 0,                // enable matching histograms (1) and extended QA (2), only QA(3), all disabled (0)
                                Bool_t    isUsingTHnSparse            = kTRUE,               // enable or disable usage of THnSparses for background estimation
                                Bool_t    enableV0findingEffi         = kFALSE,              // enables V0finding efficiency histograms
                                Bool_t    enableTriggerMimicking      = kFALSE,              // enable trigger mimicking
                                Bool_t    enableTriggerOverlapRej     = kFALSE,              // enable trigger overlap rejection                  
                                Float_t   maxFacPtHard                = 3,                // maximum factor between hardest jet and ptHard generated
                                TString   periodNameV0Reader          = ""
) {

  // ================= Load Librariers =================================
  gSystem->Load("libCore");  
  gSystem->Load("libTree");
  gSystem->Load("libGeom");
  gSystem->Load("libVMC");
  gSystem->Load("libPhysics");
  gSystem->Load("libMinuit");
  gSystem->Load("libSTEERBase");
  gSystem->Load("libESD");
  gSystem->Load("libAOD");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");  
  gSystem->Load("libCDB");
  gSystem->Load("libSTEER");
  gSystem->Load("libSTEERBase");
  gSystem->Load("libTender");
  gSystem->Load("libTenderSupplies");
  gSystem->Load("libPWGflowBase");
  gSystem->Load("libPWGflowTasks");
  gSystem->Load("libPWGGAGammaConv");

  Int_t isHeavyIon = 2;
  
  // ================== GetAnalysisManager ===============================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error(Form("AddTask_GammaConvCalo_pPb_%i",trainConfig), "No analysis manager found.");
    return ;
  }

  // ================== GetInputEventHandler =============================
  AliVEventHandler *inputHandler=mgr->GetInputEventHandler();

  Bool_t isMCForOtherTasks = kFALSE;
  if (isMC > 0) isMCForOtherTasks = kTRUE;

  
  //========= Add PID Reponse to ANALYSIS manager ====
  if(!(AliPIDResponse*)mgr->GetTask("PIDResponseTask")){
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
    AddTaskPIDResponse(isMCForOtherTasks);
  }
  
  Printf("here \n");
  
  //=========  Set Cutnumber for V0Reader ================================
  TString cutnumberPhoton   = "00000008400100001500000000";
  TString cutnumberEvent    = "80000003";
  Bool_t doEtaShift         = kFALSE;
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();

  //========= Add V0 Reader to  ANALYSIS manager if not yet existent =====
  TString V0ReaderName = Form("V0ReaderV1_%s_%s",cutnumberEvent.Data(),cutnumberPhoton.Data());
  if( !(AliV0ReaderV1*)mgr->GetTask(V0ReaderName.Data()) ){
    AliV0ReaderV1 *fV0ReaderV1 = new AliV0ReaderV1(V0ReaderName.Data());
    if (periodNameV0Reader.CompareTo("") != 0) fV0ReaderV1->SetPeriodName(periodNameV0Reader);
    fV0ReaderV1->SetUseOwnXYZCalculation(kTRUE);
    fV0ReaderV1->SetCreateAODs(kFALSE);// AOD Output
    fV0ReaderV1->SetUseAODConversionPhoton(kTRUE);
    fV0ReaderV1->SetProduceV0FindingEfficiency(enableV0findingEffi);

    if (!mgr) {
      Error("AddTask_V0ReaderV1", "No analysis manager found.");
      return;
    }

    AliConvEventCuts *fEventCuts=NULL;
    if(cutnumberEvent!=""){
      fEventCuts= new AliConvEventCuts(cutnumberEvent.Data(),cutnumberEvent.Data());
      fEventCuts->SetPreSelectionCutFlag(kTRUE);
      fEventCuts->SetV0ReaderName(V0ReaderName);
      if(fEventCuts->InitializeCutsFromCutString(cutnumberEvent.Data())){
        fEventCuts->DoEtaShift(doEtaShift);
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
  AliAnalysisTaskGammaConvCalo *task=NULL;
  task= new AliAnalysisTaskGammaConvCalo(Form("GammaConvCalo_%i",trainConfig));
  task->SetIsHeavyIon(isHeavyIon);
  task->SetIsMC(isMC);
  task->SetV0ReaderName(V0ReaderName);

  //create cut handler
  CutHandlerConvCalo cuts;

  // cluster cuts
  // 0 "ClusterType",  1 "EtaMin", 2 "EtaMax", 3 "PhiMin", 4 "PhiMax", 5 "DistanceToBadChannel", 6 "Timing", 7 "TrackMatching", 8 "ExoticCell",
  // 9 "MinEnergy", 10 "MinNCells", 11 "MinM02", 12 "MaxM02", 13 "MinM20", 14 "MaxM20", 15 "MaximumDispersion", 16 "NLM"
  
  //************************************************ EMCAL clusters **********************************************************
  if (trainConfig == 1){ // min energy = 0.3 GeV/c
    cuts.AddCut("80000013","00200009327000008250400000","1111100053022230000","0163103100000010"); //standart cut, kINT7
    cuts.AddCut("80052013","00200009327000008250400000","1111100053022230000","0163103100000010"); //standard cut, kEMC7
  } else if (trainConfig == 2){  // min energy = 0.3 GeV/c
    cuts.AddCut("80083013","00200009327000008250400000","1111100053022230000","0163103100000010"); //standard cut, kEMCEG1 based on INT7
    cuts.AddCut("80085013","00200009327000008250400000","1111100053022230000","0163103100000010"); //standard cut, kEMCEG2 based on INT7
  } else if (trainConfig == 3){ // min energy = 0.4 GeV/c
    cuts.AddCut("80000013","00200009327000008250400000","1111100053032230000","0163103100000010"); //standart cut, kINT7
    cuts.AddCut("80052013","00200009327000008250400000","1111100053032230000","0163103100000010"); //standard cut, kEMC7
  } else if (trainConfig == 4){ // min energy = 0.4 GeV/
    cuts.AddCut("80083013","00200009327000008250400000","1111100053032230000","0163103100000010"); //standard cut, kEMCEG1 based on INT7
    cuts.AddCut("80085013","00200009327000008250400000","1111100053032230000","0163103100000010"); //standard cut, kEMCEG2 based on INT7
  } else if (trainConfig == 5){ //EMCAL variation of track matching
    cuts.AddCut("80000013","00200009327000008250400000","1111141051032230000","0163103100000010"); //
    cuts.AddCut("80000013","00200009327000008250400000","1111141052032230000","0163103100000010");
    cuts.AddCut("80000013","00200009327000008250400000","1111141053032230000","0163103100000010");
    cuts.AddCut("80000013","00200009327000008250400000","1111141054032230000","0163103100000010");
    cuts.AddCut("80000013","00200009327000008250400000","1111141055032230000","0163103100000010");
    cuts.AddCut("80000013","00200009327000008250400000","1111141056032230000","0163103100000010");
  } else if (trainConfig == 6){ 
    cuts.AddCut("80052013","00200009327000008250400000","1111141053062230000","0163103100000010");
    cuts.AddCut("80052013","00200009327000008250400000","1111141053072230000","0163103100000010");
    cuts.AddCut("80052013","00200009327000008250400000","1111141053082230000","0163103100000010");
    cuts.AddCut("80052013","00200009327000008250400000","1111141053092230000","0163103100000010");

  } else if (trainConfig == 20){ // EMCAL clusters standard cuts
    cuts.AddCut("80000013","00200009327000008250400000","1111141053032230000","0163103100000010"); // INT7
    cuts.AddCut("80052013","00200009327000008250400000","1111141053032230000","0163103100000010"); // EMC7
    cuts.AddCut("80085013","00200009327000008250400000","1111141053022230000","0163103100000010"); // EG2
    cuts.AddCut("80083013","00200009327000008250400000","1111141053022230000","0163103100000010"); // EG1
    
    // minimum bias variations  
  } else if (trainConfig == 21){ //EMCAL minEnergy variation
    cuts.AddCut("80000013","00200009327000008250400000","1111141053032230000","0163103100000010"); //0.6 GeV/c default
    cuts.AddCut("80000013","00200009327000008250400000","1111141053042230000","0163103100000010"); //0.7 GeV/c
    cuts.AddCut("80000013","00200009327000008250400000","1111141053052230000","0163103100000010"); //0.9 GeV/c
  } else if (trainConfig == 22){ //EMCAL minNCells variation
    cuts.AddCut("80000013","00200009327000008250400000","1111141053031230000","0163103100000010"); //n cells >= 1
    cuts.AddCut("80000013","00200009327000008250400000","1111141053033230000","0163103100000010"); //n cells >= 3
    cuts.AddCut("80000013","00200009327000008250400000","1111141053032200000","0163103100000010"); //no max M02 cut
    cuts.AddCut("80000013","00200009327000008250400000","1112141053032230000","0163103100000010"); //only modules with TRD infront
    cuts.AddCut("80000013","00200009327000008250400000","1111341053032230000","0163103100000010"); //no modules with TRD infront
  } else if (trainConfig == 23){ // EMCAL track matching variations 
    cuts.AddCut("80000013","00200009327000008250400000","1111141051032230000","0163103100000010"); //
    cuts.AddCut("80000013","00200009327000008250400000","1111141052032230000","0163103100000010"); //
    cuts.AddCut("80000013","00200009327000008250400000","1111141053032230000","0163103100000010"); //
    cuts.AddCut("80000013","00200009327000008250400000","1111141054032230000","0163103100000010"); //
    cuts.AddCut("80000013","00200009327000008250400000","1111141055032230000","0163103100000010"); //
    cuts.AddCut("80000013","00200009327000008250400000","1111141056032230000","0163103100000010"); //
  } else if (trainConfig == 24){ // PCM variations
    cuts.AddCut("80000013","00200009227000008250400000","1111141053032230000","0163103100000010"); // dEdx e -3, 5
    cuts.AddCut("80000013","00200009127000008250400000","1111141053032230000","0163103100000010"); // dEdx e -5, 5
    cuts.AddCut("80000013","00200009357000008250400000","1111141053032230000","0163103100000010"); // dEdx pi 2
    cuts.AddCut("80000013","00200009317000008250400000","1111141053032230000","0163103100000010"); // dEdx pi 0
    cuts.AddCut("80000013","00200009387300008250400000","1111141053032230000","0163103100000010"); // dEdx pi 2 high 1 (> 3.5 GeV)
  } else if (trainConfig == 25){ // PCM variations
    cuts.AddCut("80000013","00200009327000009250400000","1111141053032230000","0163103100000010"); // qt 2D 0.03
    cuts.AddCut("80000013","00200009327000003250400000","1111141053032230000","0163103100000010"); // qt 1D 0.05
    cuts.AddCut("80000013","00200009327000002250400000","1111141053032230000","0163103100000010"); // qt 1D 0.07
    cuts.AddCut("80000013","00200049327000008250400000","1111141053032230000","0163103100000010"); // single pt > 0.075
    cuts.AddCut("80000013","00200019327000008250400000","1111141053032230000","0163103100000010"); // single pt > 0.1
  } else if (trainConfig == 26){ // PCM variations
    cuts.AddCut("80000013","00200009327000008850400000","1111141053032230000","0163103100000010"); // 2D psi pair chi2 var
    cuts.AddCut("80000013","00200009327000008260400000","1111141053032230000","0163103100000010"); // 2D psi pair chi2 var
    cuts.AddCut("80000013","00200009327000008860400000","1111141053032230000","0163103100000010"); // 2D psi pair chi2 var
    cuts.AddCut("80000013","00200009327000008280400000","1111141053032230000","0163103100000010"); // 2D psi pair chi2 var
    cuts.AddCut("80000013","00200009327000008880400000","1111141053032230000","0163103100000010"); // 2D psi pair chi2 var
  } else if (trainConfig == 27){ // PCM variations
    cuts.AddCut("80000013","00200006327000008250400000","1111141053032230000","0163103100000010"); // min TPC cl > 0.7
    cuts.AddCut("80000013","00200008327000008250400000","1111141053032230000","0163103100000010"); // min TPC cl > 0.35
    cuts.AddCut("80000013","00200009327000008250400000","1111141053032230000","0163106100000010"); // alpha < 0.8
    cuts.AddCut("80000013","00200009327000008250400000","1111141053032230000","0163105100000010"); // alpha < 0.75
  } else if (trainConfig == 28){ // PCM variations
    cuts.AddCut("80000013","00202209327000008250400000","1111141053032230000","0163103100000010"); // restrict acceptance to EMCAL loose
    cuts.AddCut("80000013","00204409327000008250400000","1111141053032230000","0163103100000010"); // restrict acceptance to EMCAL tight
  } else if (trainConfig == 29){ // PCM variations pi dEdx
    cuts.AddCut("80000013","00200009317300008250400000","1111141053032230000","0163103100000010"); // dEdx pi: 0: 0.4-3.5, -10: 3.5 ->
    cuts.AddCut("80000013","00200009327300008250400000","1111141053032230000","0163103100000010"); // dEdx pi: 1: 0.4-3.5, -10: 3.5 ->
    cuts.AddCut("80000013","00200009325000008250400000","1111141053032230000","0163103100000010"); // dEdx pi: 1: 0.3 ->
    cuts.AddCut("80000013","00200009320000008250400000","1111141053032230000","0163103100000010"); // dEdx pi: 1: 0.5 ->
  } else if (trainConfig == 30){ // PCM variations pi dEdx  
    cuts.AddCut("80000013","00200009327600008250400000","1111141053032230000","0163103100000010"); // dEdx pi: 1: 0.4-2, -10: 2. ->
    cuts.AddCut("80000013","00200009327400008250400000","1111141053032230000","0163103100000010"); // dEdx pi: 1: 0.4-3, -10: 3. ->
    cuts.AddCut("80000013","00200009315600008250400000","1111141053032230000","0163103100000010"); // dEdx pi: 0: 0.3-2, -10: 2. ->
    cuts.AddCut("80000013","00200009367400008250400000","1111141053032230000","0163103100000010"); // dEdx pi: 2: 0.4-3, 0.5: 3. ->
    cuts.AddCut("80000013","00200009347400008250400000","1111141053032230000","0163103100000010"); // dEdx pi: 3: 0.4-3, 1: 3. ->
  } else if (trainConfig == 31){ // PCM variations to close V0s  
    cuts.AddCut("80000013","00200009327000008250400000","1111141053032230000","0163103100000010"); 
    cuts.AddCut("80000013","00200009327000008250401000","1111141053032230000","0163103100000010"); 
    cuts.AddCut("80000013","00200009327000008250402000","1111141053032230000","0163103100000010"); 
    cuts.AddCut("80000013","00200009327000008250403000","1111141053032230000","0163103100000010"); 
  } else if (trainConfig == 32){ // EMCal cluster, non lin variations
    cuts.AddCut("80000013","00200009327000008250400000","1111142053032230000","0163103100000010");
    cuts.AddCut("80000013","00200009327000008250400000","1111143053032230000","0163103100000010");
    cuts.AddCut("80000013","00200009327000008250400000","1111151053032230000","0163103100000010");
    cuts.AddCut("80000013","00200009327000008250400000","1111152053032230000","0163103100000010");
  } else if (trainConfig == 33){  // EMCal cluster, non lin variations
    cuts.AddCut("80000013","00200009327000008250400000","1111101053032230000","0163103100000010"); 
    cuts.AddCut("80000013","00200009327000008250400000","1111102053032230000","0163103100000010"); 
    cuts.AddCut("80000013","00200009327000008250400000","1111144053032230000","0163103100000010"); 

  } else if (trainConfig == 40){ // non linearity variations EMC7
    cuts.AddCut("80052013","00200009327000008250400000","1111100053032230000","0163103100000010"); // no non linearity
    cuts.AddCut("80052013","00200009327000008250400000","1111101053032230000","0163103100000010"); // kSDM
    cuts.AddCut("80052013","00200009327000008250400000","1111141053032230000","0163103100000010"); // conv calo
    cuts.AddCut("80052013","00200009327000008250400000","1111142053032230000","0163103100000010"); // calo
    cuts.AddCut("80052013","00200009327000008250400000","1111151053032230000","0163103100000010"); // conv calo
    cuts.AddCut("80052013","00200009327000008250400000","1111152053032230000","0163103100000010"); // calo    
  } else if (trainConfig == 41){ // non linearity variations EG2
    cuts.AddCut("80085013","00200009327000008250400000","1111100053032230000","0163103100000010"); // no non linearity
    cuts.AddCut("80085013","00200009327000008250400000","1111101053032230000","0163103100000010"); // kSDM
    cuts.AddCut("80085013","00200009327000008250400000","1111141053032230000","0163103100000010"); // conv calo
    cuts.AddCut("80085013","00200009327000008250400000","1111142053032230000","0163103100000010"); // calo
    cuts.AddCut("80085013","00200009327000008250400000","1111151053032230000","0163103100000010"); // conv calo
    cuts.AddCut("80085013","00200009327000008250400000","1111152053032230000","0163103100000010"); // calo
  } else if (trainConfig == 42){ // non linearity variations EG1
    cuts.AddCut("80083013","00200009327000008250400000","1111100053032230000","0163103100000010"); // no non linearity
    cuts.AddCut("80083013","00200009327000008250400000","1111101053032230000","0163103100000010"); // kSDM
    cuts.AddCut("80083013","00200009327000008250400000","1111141053032230000","0163103100000010"); // conv calo
    cuts.AddCut("80083013","00200009327000008250400000","1111142053032230000","0163103100000010"); // calo
    cuts.AddCut("80083013","00200009327000008250400000","1111151053032230000","0163103100000010"); // conv calo
    cuts.AddCut("80083013","00200009327000008250400000","1111152053032230000","0163103100000010"); // calo


    //************************************************ PHOS clusters **********************************************************  
  } else if (trainConfig == 301) {  // min energy = 0.3 GeV/c
    cuts.AddCut("80000013","00200009327000008250400000","2444400048033200000","0163103100000010"); //standart cut, kINT7
    cuts.AddCut("80062013","00200009327000008250400000","2444400048033200000","0163103100000010"); //standard cut, kPHI7
  } else if (trainConfig == 302) { //PHOS
    cuts.AddCut("80000013","00200009327000008250400000","2444400047033200000","0163103100000010");
    cuts.AddCut("80000013","00200009327000008250400000","2444400048033200000","0163103100000010");
    cuts.AddCut("80000013","00200009327000008250400000","2444400049033200000","0163103100000010");
  } else if (trainConfig == 303) { //PHOS
    cuts.AddCut("80000023","00200009327000008250400000","2444400047033200000","0163103100000010");
    cuts.AddCut("80000023","00200009327000008250400000","2444400048033200000","0163103100000010");
    cuts.AddCut("80000023","00200009327000008250400000","2444400049033200000","0163103100000010");
  } else {
    Error(Form("GammaConvCalo_%i",trainConfig), "wrong trainConfig variable no cuts have been specified for the configuration");
    return;
  }

  if(!cuts.AreValid()){
    cout << "\n\n****************************************************" << endl;
    cout << "ERROR: No valid cuts stored in CutHandlerConvCalo! Returning..." << endl;
    cout << "****************************************************\n\n" << endl;
    return;
  }

  TList *EventCutList   = new TList();
  TList *ConvCutList    = new TList();
  TList *ClusterCutList = new TList();
  TList *MesonCutList   = new TList();

  Int_t numberOfCuts = cuts.GetNCuts();

  TList *HeaderList     = new TList();
  if (doWeightingPart==1) {
    TObjString *Header1 = new TObjString("pi0_1");
    HeaderList->Add(Header1);
  }
  if (doWeightingPart==2){
    TObjString *Header3 = new TObjString("eta_2");
    HeaderList->Add(Header3);
  }
  if (doWeightingPart==3) {
    TObjString *Header1 = new TObjString("pi0_1");
    HeaderList->Add(Header1);
    TObjString *Header3 = new TObjString("eta_2");
    HeaderList->Add(Header3);
  }

  EventCutList->SetOwner(kTRUE);
  AliConvEventCuts **analysisEventCuts        = new AliConvEventCuts*[numberOfCuts];
  ConvCutList->SetOwner(kTRUE);
  AliConversionPhotonCuts **analysisCuts      = new AliConversionPhotonCuts*[numberOfCuts];
  ClusterCutList->SetOwner(kTRUE);
  AliCaloPhotonCuts **analysisClusterCuts     = new AliCaloPhotonCuts*[numberOfCuts];
  MesonCutList->SetOwner(kTRUE);
  AliConversionMesonCuts **analysisMesonCuts  = new AliConversionMesonCuts*[numberOfCuts];

  for(Int_t i = 0; i<numberOfCuts; i++){
    analysisEventCuts[i] = new AliConvEventCuts();   
    analysisEventCuts[i]->SetTriggerMimicking(enableTriggerMimicking);
    analysisEventCuts[i]->SetTriggerOverlapRejecion(enableTriggerOverlapRej);
    analysisEventCuts[i]->SetMaxFacPtHard(maxFacPtHard);
    analysisEventCuts[i]->SetV0ReaderName(V0ReaderName);
    analysisEventCuts[i]->InitializeCutsFromCutString((cuts.GetEventCut(i)).Data());
    EventCutList->Add(analysisEventCuts[i]);
    analysisEventCuts[i]->SetFillCutHistograms("",kFALSE);
    
    analysisCuts[i] = new AliConversionPhotonCuts();
    analysisCuts[i]->SetV0ReaderName(V0ReaderName);
    analysisCuts[i]->InitializeCutsFromCutString((cuts.GetPhotonCut(i)).Data());
    analysisCuts[i]->SetIsHeavyIon(isHeavyIon);
    ConvCutList->Add(analysisCuts[i]);
    analysisCuts[i]->SetFillCutHistograms("",kFALSE);
  
    analysisClusterCuts[i] = new AliCaloPhotonCuts((isMC==2));
    analysisClusterCuts[i]->SetV0ReaderName(V0ReaderName);
    analysisClusterCuts[i]->InitializeCutsFromCutString((cuts.GetClusterCut(i)).Data());
    ClusterCutList->Add(analysisClusterCuts[i]);
    analysisClusterCuts[i]->SetExtendedMatchAndQA(enableExtMatchAndQA);
    analysisClusterCuts[i]->SetFillCutHistograms("");
    
    analysisMesonCuts[i] = new AliConversionMesonCuts();
    analysisMesonCuts[i]->InitializeCutsFromCutString((cuts.GetMesonCut(i)).Data());
    MesonCutList->Add(analysisMesonCuts[i]);
    analysisMesonCuts[i]->SetFillCutHistograms("");
    analysisEventCuts[i]->SetAcceptedHeader(HeaderList);
  }

  task->SetEventCutList(numberOfCuts,EventCutList);
  task->SetConversionCutList(numberOfCuts,ConvCutList);
  task->SetCaloCutList(numberOfCuts,ClusterCutList);
  task->SetMesonCutList(numberOfCuts,MesonCutList);
  task->SetMoveParticleAccordingToVertex(kTRUE);
  task->SetDoMesonAnalysis(kTRUE);
  task->SetDoMesonQA(enableQAMesonTask); //Attention new switch for Pi0 QA
  task->SetDoPhotonQA(enableQAPhotonTask);  //Attention new switch small for Photon QA
  task->SetDoClusterQA(1);  //Attention new switch small for Cluster QA
  task->SetUseTHnSparse(isUsingTHnSparse);
  if(enableExtMatchAndQA == 2 || enableExtMatchAndQA == 3){ task->SetPlotHistsExtQA(kTRUE);}

  //connect containers
  AliAnalysisDataContainer *coutput =
    mgr->CreateContainer(Form("GammaConvCalo_%i",trainConfig), TList::Class(),
              AliAnalysisManager::kOutputContainer,Form("GammaConvCalo_%i.root",trainConfig));

  mgr->AddTask(task);
  mgr->ConnectInput(task,0,cinput);
  mgr->ConnectOutput(task,1,coutput);

  return;

}
