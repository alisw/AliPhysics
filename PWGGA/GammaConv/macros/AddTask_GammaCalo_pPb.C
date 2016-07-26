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
//($ALIPHYSICS/PWGGA/GammaConv/AliAnalysisTaskGammaCalo.cxx) for
//pPb together with all supporting classes
//***************************************************************************************

//***************************************************************************************
//CutHandler contains all cuts for a certain analysis and trainconfig,
//it automatically checks length of cutStrings and takes care of the number of added cuts,
//no specification of the variable 'numberOfCuts' needed anymore.
//***************************************************************************************
class CutHandlerCalo{
  public:
    CutHandlerCalo(Int_t nMax=10){
      nCuts=0; nMaxCuts=nMax; validCuts = true;
      eventCutArray = new TString[nMaxCuts]; clusterCutArray = new TString[nMaxCuts]; mesonCutArray = new TString[nMaxCuts];
      for(Int_t i=0; i<nMaxCuts; i++) {eventCutArray[i] = ""; clusterCutArray[i] = ""; mesonCutArray[i] = "";}
    }

    void AddCut(TString eventCut, TString clusterCut, TString mesonCut){
      if(nCuts>=nMaxCuts) {cout << "ERROR in CutHandlerCalo: Exceeded maximum number of cuts!" << endl; validCuts = false; return;}
      if( eventCut.Length()!=8 || clusterCut.Length()!=19 || mesonCut.Length()!=16 ) {cout << "ERROR in CutHandlerCalo: Incorrect length of cut string!" << endl; validCuts = false; return;}
      eventCutArray[nCuts]=eventCut; clusterCutArray[nCuts]=clusterCut; mesonCutArray[nCuts]=mesonCut;
      nCuts++;
      return;
    }
    Bool_t AreValid(){return validCuts;}
    Int_t GetNCuts(){if(validCuts) return nCuts; else return 0;}
    TString GetEventCut(Int_t i){if(validCuts&&i<nMaxCuts&&i>=0) return eventCutArray[i]; else{cout << "ERROR in CutHandlerCalo: GetEventCut wrong index i" << endl;return "";}}
    TString GetClusterCut(Int_t i){if(validCuts&&i<nMaxCuts&&i>=0) return clusterCutArray[i]; else {cout << "ERROR in CutHandlerCalo: GetClusterCut wrong index i" << endl;return "";}}
    TString GetMesonCut(Int_t i){if(validCuts&&i<nMaxCuts&&i>=0) return mesonCutArray[i]; else {cout << "ERROR in CutHandlerCalo: GetMesonCut wrong index i" << endl;return "";}}
  private:
    Bool_t validCuts;
    Int_t nCuts; Int_t nMaxCuts;
    TString* eventCutArray;
    TString* clusterCutArray;
    TString* mesonCutArray;
};

//***************************************************************************************
//main function
//***************************************************************************************
void AddTask_GammaCalo_pPb(  
                            Int_t      trainConfig                = 1,                // change different set of cuts
                            Int_t      isMC                       = 0,                // run MC
                            Int_t      enableQAMesonTask          = 0,                // enable QA in AliAnalysisTaskGammaConvV1
                            Int_t      enableQAClusterTask        = 0,                // enable additional QA task
                            TString    fileNameInputForWeighting  = "MCSpectraInput.root",       // path to file for weigting input
                            Int_t      doWeightingPart            = 0,                // enable Weighting
                            TString    generatorName              = "DPMJET",
                            TString    cutnumberAODBranch         = "800000006008400000001500000",   // cutnumber for AOD branch
                            Bool_t     isUsingTHnSparse           = kFALSE,           // enable or disable usage of THnSparses for background estimation
                            Int_t      enableExtMatchAndQA        = 0,                // enable QA(3), extMatch+QA(2), extMatch(1), disabled (0)
                            Bool_t     enableTriggerMimicking     = kFALSE,           // enable trigger mimicking
                            Bool_t     enableTriggerOverlapRej    = kTRUE,            // enable trigger overlap rejection
                            Float_t    maxFacPtHard               = 3,                // maximum factor between hardest jet and ptHard generated
                            TString    periodNameV0Reader         = ""
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
    Error(Form("AddTask_GammaCalo_pPb_%i",trainConfig), "No analysis manager found.");
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
  
  Printf("here \n");
  
  //=========  Set Cutnumber for V0Reader ================================
  TString cutnumberPhoton   = "00000008400100001500000000";
  TString cutnumberEvent    = "80000003";
  Bool_t doEtaShift         = kFALSE;
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();

  //========= Add V0 Reader to  ANALYSIS manager if not yet existent =====
  if( !(AliV0ReaderV1*)mgr->GetTask("V0ReaderV1") ){
    AliV0ReaderV1 *fV0ReaderV1 = new AliV0ReaderV1("V0ReaderV1");
    if (periodNameV0Reader.CompareTo("") != 0) fV0ReaderV1->SetPeriodName(periodNameV0Reader);
    fV0ReaderV1->SetCreateAODs(kFALSE);// AOD Output
    fV0ReaderV1->SetUseAODConversionPhoton(kTRUE);

    if (!mgr) {
      Error("AddTask_V0ReaderV1", "No analysis manager found.");
      return;
    }
    AliConvEventCuts *fEventCuts=NULL;
    if(cutnumberEvent!=""){
      fEventCuts = new AliConvEventCuts(cutnumberEvent.Data(),cutnumberEvent.Data());
      fEventCuts->SetPreSelectionCutFlag(kTRUE);
      if(fEventCuts->InitializeCutsFromCutString(cutnumberEvent.Data())){
        fEventCuts->DoEtaShift(doEtaShift);
        fV0ReaderV1->SetEventCuts(fEventCuts);
        fEventCuts->SetFillCutHistograms("",kTRUE);
      }
    }
    // Set AnalysisCut Number
    AliConversionPhotonCuts *fCuts=NULL;
    if(cutnumberPhoton!=""){
      fCuts = new AliConversionPhotonCuts(cutnumberPhoton.Data(),cutnumberPhoton.Data());
      fCuts->SetPreSelectionCutFlag(kTRUE);
      fCuts->SetIsHeavyIon(isHeavyIon);
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
  AliAnalysisTaskGammaCalo *task=NULL;
  task= new AliAnalysisTaskGammaCalo(Form("GammaCalo_%i",trainConfig));
  task->SetIsHeavyIon(isHeavyIon);
  task->SetIsMC(isMC);

  //create cut handler
  CutHandlerCalo cuts;

  // cluster cuts
  // 0 "ClusterType",  1 "EtaMin", 2 "EtaMax", 3 "PhiMin", 4 "PhiMax", 5 "DistanceToBadChannel", 6 "Timing", 7 "TrackMatching", 8 "ExoticCell",
  // 9 "MinEnergy", 10 "MinNCells", 11 "MinM02", 12 "MaxM02", 13 "MinM20", 14 "MaxM20", 15 "MaximumDispersion", 16 "NLM"
  
  //************************************************ EMCAL clusters *************************************************
  if (trainConfig == 1){ // min energy = 0.3 GeV/c
    cuts.AddCut("80000013","1111141050022230000","0163103100000050"); //standart cut, kINT7 // EMCAL clusters
    cuts.AddCut("80052013","1111141050022230000","0163103100000050"); //standard cut, kEMC7 // EMCAL clusters
  } else if (trainConfig == 2){  // min energy = 0.3 GeV/c
    cuts.AddCut("80083013","1111141050022230000","0163103100000050"); //standard cut, kEMCEG1 based on INT7 // EMCAL clusters
    cuts.AddCut("80085013","1111141050022230000","0163103100000050"); //standard cut, kEMCEG2 based on INT7 // EMCAL clusters
  } else if (trainConfig == 3){ // min energy = 0.4 GeV/c
    cuts.AddCut("80000013","1111141050032230000","0163103100000050"); //standart cut, kINT7 // EMCAL clusters
    cuts.AddCut("80052013","1111141050032230000","0163103100000050"); //standard cut, kEMC7 // EMCAL clusters
  } else if (trainConfig == 4){ // min energy = 0.4 GeV/c
    cuts.AddCut("80083013","1111141050032230000","0163103100000050"); //standard cut, kEMCEG1 based on INT7 // EMCAL clusters
    cuts.AddCut("80085013","1111141050032230000","0163103100000050"); //standard cut, kEMCEG2 based on INT7 // EMCAL clusters
  } else if (trainConfig == 5){ //EMCAL minEnergy variation
    cuts.AddCut("80000013","1111141050012230000","0163103100000050"); //0.5 GeV/c
    cuts.AddCut("80000013","1111141050022230000","0163103100000050"); //0.6 GeV/c
    cuts.AddCut("80000013","1111141050032230000","0163103100000050"); //0.7 GeV/c default
    cuts.AddCut("80000013","1111141050042230000","0163103100000050"); //0.8 GeV/c
    cuts.AddCut("80000013","1111141050052230000","0163103100000050"); //0.9 GeV/c
  } else if (trainConfig == 6){ //EMCAL minNCells variation
    cuts.AddCut("80000013","1111141050031230000","0163103100000050"); //n cells >= 1
    cuts.AddCut("80000013","1111141050033230000","0163103100000050"); //n cells >= 3
    cuts.AddCut("80000013","1111141050032000000","0163103100000050"); //no M02 cut
    cuts.AddCut("80000013","1112141050032230000","0163103100000050"); //only modules with TRD infront
    cuts.AddCut("80000013","1111341050032230000","0163103100000050"); //no modules with TRD infront
  } else if (trainConfig == 7){ // Validation EMCAL
    cuts.AddCut("80000013","1111141050032000000","0163103100000050");
  } else if (trainConfig == 8){ // Validation EMCAL, only added signals
    cuts.AddCut("80000023","1111141050032000000","0163103100000050");
  } else if (trainConfig == 9){ // non linearity variations INT7
    cuts.AddCut("80000013","1111100050032230000","0163103100000050"); // non nonlinearity
    cuts.AddCut("80000013","1111101050032230000","0163103100000050"); // kSDM
    cuts.AddCut("80000013","1111141050032230000","0163103100000050"); // conv calo
    cuts.AddCut("80000013","1111142050032230000","0163103100000050"); // calo
    cuts.AddCut("80000013","1111149050032230000","0163103100000050"); // calo
  } else if (trainConfig == 10){ // non linearity variations EMC7
    cuts.AddCut("80052013","1111100050032230000","0163103100000050"); // non nonlinearity
    cuts.AddCut("80052013","1111101050032230000","0163103100000050"); // kSDM
    cuts.AddCut("80052013","1111141050032230000","0163103100000050"); // conv calo
    cuts.AddCut("80052013","1111142050032230000","0163103100000050"); // calo
  } else if (trainConfig == 11){ // non linearity variations EG2
    cuts.AddCut("80085013","1111100050032230000","0163103100000050"); // non nonlinearity
    cuts.AddCut("80085013","1111101050032230000","0163103100000050"); // kSDM
    cuts.AddCut("80085013","1111141050032230000","0163103100000050"); // conv calo
    cuts.AddCut("80085013","1111142050032230000","0163103100000050"); // calo
  } else if (trainConfig == 12){ // non linearity variations EG1
    cuts.AddCut("80083013","1111100050032230000","0163103100000050"); // non nonlinearity
    cuts.AddCut("80083013","1111101050032230000","0163103100000050"); // kSDM
    cuts.AddCut("80083013","1111141050032230000","0163103100000050"); // conv calo
    cuts.AddCut("80083013","1111142050032230000","0163103100000050"); // calo
  } else if (trainConfig == 13){ // no non linearity
    cuts.AddCut("80000013","1111100050032230000","0163103100000050"); // kINT7 // EMCAL clusters
    cuts.AddCut("80052013","1111100050032230000","0163103100000050"); // kEMC7 // EMCAL clusters
    cuts.AddCut("80083013","1111100050032230000","0163103100000050"); // kEMCEG1 based on INT7 // EMCAL clusters
    cuts.AddCut("80085013","1111100050032230000","0163103100000050"); // kEMCEG2 based on INT7 // EMCAL clusters
  } else if(trainConfig == 14){ // variation opening angle
    cuts.AddCut("80000013","1111141050022230000","0163103100000050"); // standard
    cuts.AddCut("80000013","1111141050022230000","0163103100000060"); // 2 EMCal cell diagonals
    cuts.AddCut("80000013","1111141050022230000","0163103100000040"); // 0.75 EMCal cell diagonals
  } else if(trainConfig == 15){
    cuts.AddCut("80000013","1111141050032230000","0163403100000050"); // MB
    cuts.AddCut("80052013","1111141050032230000","0163403100000050"); // EMC7
    cuts.AddCut("80083013","1111141050032230000","0163403100000050"); // EG1
    cuts.AddCut("80085013","1111141050032230000","0163403100000050"); // EG2
  } else if (trainConfig == 16){ // non linearity variations INT7 with TM
    cuts.AddCut("80000013","1111100051032230000","0163103100000050"); // non nonlinearity
    cuts.AddCut("80000013","1111101051032230000","0163103100000050"); // kSDM
    cuts.AddCut("80000013","1111141051032230000","0163103100000050"); // conv calo
    cuts.AddCut("80000013","1111142051032230000","0163103100000050"); // calo
    cuts.AddCut("80000013","1111149051032230000","0163103100000050"); // calo
  } else if (trainConfig == 17){ // non linearity variations INT7 
    cuts.AddCut("80000013","1111141051032230000","0163103100000050"); // TM variations
    cuts.AddCut("80000013","1111141052032230000","0163103100000050");
    cuts.AddCut("80000013","1111141053032230000","0163103100000050");
    cuts.AddCut("80000013","1111141054032230000","0163103100000050");
    cuts.AddCut("80000013","1111141055032230000","0163103100000050");
  } else if (trainConfig == 18){ // no non linearity with TM
    cuts.AddCut("80000013","1111100051032230000","0163103100000050"); // kINT7 // EMCAL clusters
    cuts.AddCut("80052013","1111100051032230000","0163103100000050"); // kEMC7 // EMCAL clusters
    cuts.AddCut("80083013","1111100051032230000","0163103100000050"); // kEMCEG1 based on INT7 // EMCAL clusters
    cuts.AddCut("80085013","1111100051032230000","0163103100000050"); // kEMCEG2 based on INT7 // EMCAL clusters

    
  // SYSTEMATIC STUDY NEUTRAl MESON MEASUREMENTS MIKE SAS 02-02-2016
  } else if(trainConfig == 40){ // default cutstring and first set of variations nonlinearity
    cuts.AddCut("80000013","1111141051032230000","0163403100000050"); // default
    cuts.AddCut("80000013","1111142051032230000","0163403100000050"); // calo nonlinearity variation
    cuts.AddCut("80000013","1111151051032230000","0163403100000050"); // calo nonlinearity variation
    cuts.AddCut("80000013","1111152051032230000","0163403100000050"); // calo nonlinearity variation
  } else if(trainConfig == 41){ // second set of variations CLUSTER
    cuts.AddCut("80000013","1111141051022230000","0163403100000050"); // min energy cluster variation 1  600 MeV
    cuts.AddCut("80000013","1111141051042230000","0163403100000050"); // min energy cluster variation 2  800 MeV
    cuts.AddCut("80000013","1111141051052230000","0163403100000050"); // min energy cluster variation 3  900 MeV
  } else if(trainConfig == 42){ // third set of variations CLUSTER
    cuts.AddCut("80000013","1111141051032030000","0163403100000050"); // min/max M02  0<M<0.5
    cuts.AddCut("80000013","1111141051032200000","0163403100000050"); // min/max M02  0.1<M<100
    cuts.AddCut("80000013","1111141051031230000","0163403100000050"); // min number of cells variation 1  1 cell
    cuts.AddCut("80000013","1112141051032230000","0163403100000050"); // only modules with TRD infront
    cuts.AddCut("80000013","1111341051032230000","0163403100000050"); // no modules with TRD infront
  } else if(trainConfig == 43){ // third set of variations MESON
    cuts.AddCut("80000013","1111141051032230000","0163303100000050"); // rapidity variation  y<0.6
    cuts.AddCut("80000013","1111141051032230000","0163103100000050"); // rapidity variation  y<0.8
    cuts.AddCut("80000013","1111141051032230000","0163406100000050"); // alpha meson variation 1   0<alpha<0.8
    cuts.AddCut("80000013","1111141051032230000","0163405100000050"); // alpha meson variation 2  0<alpha<0.75
  } else if(trainConfig == 44){ // fourth set of variations
    cuts.AddCut("80000013","1111141050032230000","0163403100000050"); // tm variation
    cuts.AddCut("80000013","1111141055032230000","0163403100000050"); // tm variation
    cuts.AddCut("80000013","1111141051032230000","0163403100000040"); // min opening angle 0.75 cell diag
    cuts.AddCut("80000013","1111141051032230000","0163403100000060"); // min opening angle 2 cell diag
  
  //EMC7
  } else if(trainConfig == 50){ // default cutstring and first set of variations nonlinearity
    cuts.AddCut("80052013","1111141051032230000","0163403100000050"); // default
    cuts.AddCut("80052013","1111142051032230000","0163403100000050"); // calo nonlinearity variation
    cuts.AddCut("80052013","1111151051032230000","0163403100000050"); // calo nonlinearity variation
    cuts.AddCut("80052013","1111152051032230000","0163403100000050"); // calo nonlinearity variation
  } else if(trainConfig == 51){ // second set of variations CLUSTER
    cuts.AddCut("80052013","1111141051022230000","0163403100000050"); // min energy cluster variation 1  600 MeV
    cuts.AddCut("80052013","1111141051042230000","0163403100000050"); // min energy cluster variation 2  800 MeV
    cuts.AddCut("80052013","1111141051052230000","0163403100000050"); // min energy cluster variation 3  900 MeV
  } else if(trainConfig == 52){ // third set of variations CLUSTER
    cuts.AddCut("80052013","1111141051032030000","0163403100000050"); // min/max M02  0<M<0.5
    cuts.AddCut("80052013","1111141051032200000","0163403100000050"); // min/max M02  0.1<M<100
    cuts.AddCut("80052013","1111141051031230000","0163403100000050"); // min number of cells variation 1  1 cell
    cuts.AddCut("80052013","1112141051032230000","0163403100000050"); // only modules with TRD infront
    cuts.AddCut("80052013","1111341051032230000","0163403100000050"); // no modules with TRD infront
  } else if(trainConfig == 53){ // third set of variations MESON
    cuts.AddCut("80052013","1111141051032230000","0163303100000050"); // rapidity variation  y<0.6
    cuts.AddCut("80052013","1111141051032230000","0163103100000050"); // rapidity variation  y<0.8
    cuts.AddCut("80052013","1111141051032230000","0163401100000050"); // alpha meson variation 1  0.5<alpha<1
    cuts.AddCut("80052013","1111141051032230000","0163402100000050"); // alpha meson variation 2  0.6<alpha<1
  } else if(trainConfig == 54){ // default cutstring for different tender settings
    cuts.AddCut("80052013","1111141050032230000","0163403100000050"); // tm variation
    cuts.AddCut("80052013","1111141055032230000","0163403100000050"); // tm variation
    cuts.AddCut("80052013","1111141051032230000","0163403100000040"); // min opening angle 0.75 cell diag
    cuts.AddCut("80052013","1111141051032230000","0163403100000060"); // min opening angle 2 cell diag

  //EG1
  } else if(trainConfig == 60){ // default cutstring and first set of variations nonlinearity
    cuts.AddCut("80083013","1111141051032230000","0163403100000050"); // default
    cuts.AddCut("80083013","1111142051032230000","0163403100000050"); // calo nonlinearity variation
    cuts.AddCut("80083013","1111151051032230000","0163403100000050"); // calo nonlinearity variation
    cuts.AddCut("80083013","1111152051032230000","0163403100000050"); // calo nonlinearity variation
  } else if(trainConfig == 61){ // second set of variations CLUSTER
    cuts.AddCut("80083013","1111141051022230000","0163403100000050"); // min energy cluster variation 1  600 MeV
    cuts.AddCut("80083013","1111141051042230000","0163403100000050"); // min energy cluster variation 2  800 MeV
    cuts.AddCut("80083013","1111141051052230000","0163403100000050"); // min energy cluster variation 3  900 MeV
  } else if(trainConfig == 62){ // third set of variations CLUSTER
    cuts.AddCut("80083013","1111141051032030000","0163403100000050"); // min/max M02  0<M<0.5
    cuts.AddCut("80083013","1111141051032200000","0163403100000050"); // min/max M02  0.1<M<100
    cuts.AddCut("80083013","1111141051031230000","0163403100000050"); // min number of cells variation 1  1 cell
    cuts.AddCut("80083013","1112141051032230000","0163403100000050"); // only modules with TRD infront
    cuts.AddCut("80083013","1111341051032230000","0163403100000050"); // no modules with TRD infront
  } else if(trainConfig == 63){ // third set of variations MESON
    cuts.AddCut("80083013","1111141051032230000","0163303100000050"); // rapidity variation  y<0.6
    cuts.AddCut("80083013","1111141051032230000","0163103100000050"); // rapidity variation  y<0.8
    cuts.AddCut("80083013","1111141051032230000","0163401100000050"); // alpha meson variation 1  0.5<alpha<1
    cuts.AddCut("80083013","1111141051032230000","0163402100000050"); // alpha meson variation 2  0.6<alpha<1
  } else if(trainConfig == 64){ // default cutstring for different tender settings
    cuts.AddCut("80083013","1111141050032230000","0163403100000050"); // tm variation
    cuts.AddCut("80083013","1111141055032230000","0163403100000050"); // tm variation
    cuts.AddCut("80083013","1111141051032230000","0163403100000040"); // min opening angle 0.75 cell diag
    cuts.AddCut("80083013","1111141051032230000","0163403100000060"); // min opening angle 2 cell diag
   //EG2
  } else if(trainConfig == 70){ // default cutstring and first set of variations nonlinearity
    cuts.AddCut("80085013","1111141051032230000","0163403100000050"); // default
    cuts.AddCut("80085013","1111142051032230000","0163403100000050"); // calo nonlinearity variation
    cuts.AddCut("80085013","1111151051032230000","0163403100000050"); // calo nonlinearity variation
    cuts.AddCut("80085013","1111152051032230000","0163403100000050"); // calo nonlinearity variation
  } else if(trainConfig == 71){ // second set of variations CLUSTER
    cuts.AddCut("80085013","1111141051022230000","0163403100000050"); // min energy cluster variation 1  600 MeV
    cuts.AddCut("80085013","1111141051042230000","0163403100000050"); // min energy cluster variation 2  800 MeV
    cuts.AddCut("80085013","1111141051052230000","0163403100000050"); // min energy cluster variation 3  900 MeV
  } else if(trainConfig == 72){ // third set of variations CLUSTER
    cuts.AddCut("80085013","1111141051032030000","0163403100000050"); // min/max M02  0<M<0.5
    cuts.AddCut("80085013","1111141051032200000","0163403100000050"); // min/max M02  0.1<M<100
    cuts.AddCut("80085013","1111141051031230000","0163403100000050"); // min number of cells variation 1  1 cell
    cuts.AddCut("80085013","1112141051032230000","0163403100000050"); // only modules with TRD infront
    cuts.AddCut("80085013","1111341051032230000","0163403100000050"); // no modules with TRD infront
  } else if(trainConfig == 73){ // third set of variations MESON
    cuts.AddCut("80085013","1111141051032230000","0163303100000050"); // rapidity variation  y<0.6
    cuts.AddCut("80085013","1111141051032230000","0163103100000050"); // rapidity variation  y<0.8
    cuts.AddCut("80085013","1111141051032230000","0163401100000050"); // alpha meson variation 1  0.5<alpha<1
    cuts.AddCut("80085013","1111141051032230000","0163402100000050"); // alpha meson variation 2  0.6<alpha<1
  } else if(trainConfig == 74){ // default cutstring for different tender settings
    cuts.AddCut("80085013","1111141050032230000","0163403100000050"); // tm variation
    cuts.AddCut("80085013","1111141055032230000","0163403100000050"); // tm variation
    cuts.AddCut("80085013","1111141051032230000","0163403100000040"); // min opening angle 0.75 cell diag
    cuts.AddCut("80085013","1111141051032230000","0163403100000060"); // min opening angle 2 cell diag

  //all default triggers
  } else if(trainConfig == 80){
    cuts.AddCut("80000013","1111141051032230000","0163403100000050"); // default MB
    cuts.AddCut("80052013","1111141051032230000","0163403100000050"); // default EMC7
    cuts.AddCut("80083013","1111141051032230000","0163403100000050"); // default EG1
    cuts.AddCut("80085013","1111141051032230000","0163403100000050"); // default EG2
    

  //************************************************ PHOS clusters *************************************************
  } else if (trainConfig == 31) {  // min energy = 0.3 GeV/c
    cuts.AddCut("80000013","2444400040033200000","0163103100000050"); //standart cut, kINT7 // PHOS clusters
    cuts.AddCut("80062013","2444400040033200000","0163103100000050"); //standard cut, kPHI7  // PHOS clusters
  } else if (trainConfig == 32){ // Validation PHOS
    cuts.AddCut("80000013","2444400040053200000","0163103100000050");
  } else if (trainConfig == 33){ // Validation PHOS, only added signals
    cuts.AddCut("80000023","2444400040053200000","0163103100000050");

    
  } else {
    Error(Form("GammaCalo_%i",trainConfig), "wrong trainConfig variable no cuts have been specified for the configuration");
    return;
  }

  if(!cuts.AreValid()){
    cout << "\n\n****************************************************" << endl;
    cout << "ERROR: No valid cuts stored in CutHandlerCalo! Returning..." << endl;
    cout << "****************************************************\n\n" << endl;
    return;
  }

  Int_t numberOfCuts = cuts.GetNCuts();

  TList *EventCutList     = new TList();
  TList *ClusterCutList   = new TList();
  TList *MesonCutList     = new TList();

  TList *HeaderList = new TList();
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
  ClusterCutList->SetOwner(kTRUE);
  AliCaloPhotonCuts **analysisClusterCuts     = new AliCaloPhotonCuts*[numberOfCuts];
  MesonCutList->SetOwner(kTRUE);
  AliConversionMesonCuts **analysisMesonCuts  = new AliConversionMesonCuts*[numberOfCuts];

  for(Int_t i = 0; i<numberOfCuts; i++){
    analysisEventCuts[i] = new AliConvEventCuts();   
    analysisEventCuts[i]->SetTriggerMimicking(enableTriggerMimicking);
    analysisEventCuts[i]->SetTriggerOverlapRejecion(enableTriggerOverlapRej);
    analysisEventCuts[i]->SetMaxFacPtHard(maxFacPtHard);
    analysisEventCuts[i]->InitializeCutsFromCutString((cuts.GetEventCut(i)).Data());
    EventCutList->Add(analysisEventCuts[i]);
    analysisEventCuts[i]->SetFillCutHistograms("",kFALSE);
      
    analysisClusterCuts[i] = new AliCaloPhotonCuts((isMC==2));
    analysisClusterCuts[i]->SetIsPureCaloCut(2);
    analysisClusterCuts[i]->InitializeCutsFromCutString((cuts.GetClusterCut(i)).Data());
    ClusterCutList->Add(analysisClusterCuts[i]);
    analysisClusterCuts[i]->SetExtendedMatchAndQA(enableExtMatchAndQA);
    analysisClusterCuts[i]->SetFillCutHistograms("");
    
    analysisMesonCuts[i] = new AliConversionMesonCuts();
    analysisMesonCuts[i]->InitializeCutsFromCutString((cuts.GetMesonCut(i)).Data());
    analysisMesonCuts[i]->SetIsMergedClusterCut(2);
    MesonCutList->Add(analysisMesonCuts[i]);
    analysisMesonCuts[i]->SetFillCutHistograms("");
    analysisEventCuts[i]->SetAcceptedHeader(HeaderList);
  }
  task->SetEventCutList(numberOfCuts,EventCutList);
  task->SetCaloCutList(numberOfCuts,ClusterCutList);
  task->SetMesonCutList(numberOfCuts,MesonCutList);
  task->SetDoMesonAnalysis(kTRUE);
  task->SetDoMesonQA(enableQAMesonTask); //Attention new switch for Pi0 QA
  task->SetDoClusterQA(enableQAClusterTask);  //Attention new switch small for Cluster QA
  task->SetDoTHnSparse(isUsingTHnSparse);
  if(enableExtMatchAndQA == 2 || enableExtMatchAndQA == 3){ task->SetPlotHistsExtQA(kTRUE);}

  //connect containers
  AliAnalysisDataContainer *coutput =
    mgr->CreateContainer(Form("GammaCalo_%i",trainConfig), TList::Class(),
          AliAnalysisManager::kOutputContainer,Form("GammaCalo_%i.root",trainConfig));

  mgr->AddTask(task);
  mgr->ConnectInput(task,0,cinput);
  mgr->ConnectOutput(task,1,coutput);

  return;

}
