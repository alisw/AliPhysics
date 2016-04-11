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
//pp together with all supporting classes
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
void AddTask_GammaConvCalo_pp(  Int_t     trainConfig                   = 1,                      //change different set of cuts
                                Int_t     isMC                          = 0,                      //run MC
                                Int_t     enableQAMesonTask             = 1,                      //enable QA in AliAnalysisTaskGammaConvV1
                                Int_t     enableQAPhotonTask            = 1,                      // enable additional QA task
                                TString   fileNameInputForPartWeighting = "MCSpectraInput.root",  // path to file for weigting input
                                TString   cutnumberAODBranch            = "000000006008400001001500000",
                                Int_t     enableExtMatchAndQA           = 0,                      // enable matching histograms (1) and extended QA (2), only QA(3), all disabled (0)
                                TString   periodname                    = "LHC12f1x",             // period name
                                Bool_t    doParticleWeighting           = kFALSE,                 // enables weighting
                                Bool_t    enableV0findingEffi           = kFALSE,                 // enables V0finding efficiency histograms
                                Bool_t    isUsingTHnSparse              = kTRUE,                  // enable or disable usage of THnSparses for background estimation
                                Bool_t    enableTriggerMimicking        = kFALSE,                 // enable trigger mimicking
                                Bool_t    enableTriggerOverlapRej       = kFALSE,                 // enable trigger overlap rejection
                                Float_t   maxFacPtHard                  = 3.,                     // maximum factor between hardest jet and ptHard generated
                                TString   periodNameV0Reader            = "",                     // 
                                Bool_t    doTreeConvGammaShape          = kFALSE,                 //
                                Bool_t    doMultiplicityWeighting       = kFALSE,                  //
                                TString   fileNameInputForMultWeighing  = "Multiplicity.root",    //
                                TString   periodNameAnchor              = ""
              ) {
  
  Int_t isHeavyIon = 0;
  
  // ================== GetAnalysisManager ===============================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error(Form("AddTask_GammaConvCalo_pp_%i",trainConfig), "No analysis manager found.");
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
  TString cutnumberPhoton = "00000008400100001500000000";
  TString cutnumberEvent = "00000003";
  Bool_t doEtaShift = kFALSE;
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

  // ************************************* EMCAL cuts ****************************************************
  // LHC11a with new non linearities
  if (trainConfig == 1){ // EMCAL clusters 2.76 TeV LHC11a, with SDD (0), kEMC1 (1)
    cuts.AddCut("00003113","00200009327000008250400000","1111121053032230000","0163103100000010"); // 700 MeV cluster min energy
    cuts.AddCut("00003013","00200009327000008250400000","1111121053032230000","0163103100000010"); // 700 MeV cluster min energy without pileup rejection
    cuts.AddCut("00051013","00200009327000008250400000","1111121053032230000","0163103100000010"); // 700 MeV cluster min energy
    
  // minimum bias variations  
  } else if (trainConfig == 2){ //EMCAL minEnergy variation
    cuts.AddCut("00003113","00200009327000008250400000","1111121053032230000","0163103100000010"); //0.6 GeV/c default
    cuts.AddCut("00003113","00200009327000008250400000","1111121053042230000","0163103100000010"); //0.7 GeV/c
    cuts.AddCut("00003113","00200009327000008250400000","1111121053052230000","0163103100000010"); //0.9 GeV/c
  } else if (trainConfig == 3){ //EMCAL minNCells variation
    cuts.AddCut("00003113","00200009327000008250400000","1111121053031230000","0163103100000010"); //n cells >= 1
    cuts.AddCut("00003113","00200009327000008250400000","1111121053033230000","0163103100000010"); //n cells >= 3
    cuts.AddCut("00003113","00200009327000008250400000","1111121053032200000","0163103100000010"); //no max M02 cut
    cuts.AddCut("00003113","00200009327000008250400000","1113111053032230000","0163103100000010"); //only modules with TRD infront
    cuts.AddCut("00003113","00200009327000008250400000","1111211053032230000","0163103100000010"); //no modules with TRD infront
  } else if (trainConfig == 4){ // EMCAL track matching variations 
    cuts.AddCut("00003113","00200009327000008250400000","1111121051032230000","0163103100000010"); //
    cuts.AddCut("00003113","00200009327000008250400000","1111121052032230000","0163103100000010"); //
    cuts.AddCut("00003113","00200009327000008250400000","1111121053032230000","0163103100000010"); //
    cuts.AddCut("00003113","00200009327000008250400000","1111121054032230000","0163103100000010"); //
    cuts.AddCut("00003113","00200009327000008250400000","1111121055032230000","0163103100000010"); //
    cuts.AddCut("00003113","00200009327000008250400000","1111121056032230000","0163103100000010"); //
  } else if (trainConfig == 5){ // PCM variations
    cuts.AddCut("00003113","00200009227000008250400000","1111121053032230000","0163103100000010"); // dEdx e -3, 5
    cuts.AddCut("00003113","00200009127000008250400000","1111121053032230000","0163103100000010"); // dEdx e -5, 5
    cuts.AddCut("00003113","00200009357000008250400000","1111121053032230000","0163103100000010"); // dEdx pi 2
    cuts.AddCut("00003113","00200009317000008250400000","1111121053032230000","0163103100000010"); // dEdx pi 0
    cuts.AddCut("00003113","00200009387300008250400000","1111121053032230000","0163103100000010"); // dEdx pi 2 high 1 (> 3.5 GeV)
  } else if (trainConfig == 6){ // PCM variations
    cuts.AddCut("00003113","00200009327000009250400000","1111121053032230000","0163103100000010"); // qt 2D 0.03
    cuts.AddCut("00003113","00200009327000003250400000","1111121053032230000","0163103100000010"); // qt 1D 0.05
    cuts.AddCut("00003113","00200009327000002250400000","1111121053032230000","0163103100000010"); // qt 1D 0.07
    cuts.AddCut("00003113","00200049327000008250400000","1111121053032230000","0163103100000010"); // single pt > 0.075
    cuts.AddCut("00003113","00200019327000008250400000","1111121053032230000","0163103100000010"); // single pt > 0.1
  } else if (trainConfig == 7){ // PCM variations
    cuts.AddCut("00003113","00200009327000008850400000","1111121053032230000","0163103100000010"); // 2D psi pair chi2 var
    cuts.AddCut("00003113","00200009327000008260400000","1111121053032230000","0163103100000010"); // 2D psi pair chi2 var
    cuts.AddCut("00003113","00200009327000008860400000","1111121053032230000","0163103100000010"); // 2D psi pair chi2 var
    cuts.AddCut("00003113","00200009327000008280400000","1111121053032230000","0163103100000010"); // 2D psi pair chi2 var
    cuts.AddCut("00003113","00200009327000008880400000","1111121053032230000","0163103100000010"); // 2D psi pair chi2 var
  } else if (trainConfig == 8){ // PCM variations
    cuts.AddCut("00003113","00200006327000008250400000","1111121053032230000","0163103100000010"); // min TPC cl > 0.7
    cuts.AddCut("00003113","00200008327000008250400000","1111121053032230000","0163103100000010"); // min TPC cl > 0.35
    cuts.AddCut("00003113","00200009327000008250400000","1111121053032230000","0163106100000010"); // alpha < 0.8
    cuts.AddCut("00003113","00200009327000008250400000","1111121053032230000","0163105100000010"); // alpha < 0.75
  } else if (trainConfig == 9){ // PCM variations
    cuts.AddCut("00003113","00202209327000008250400000","1111121053032230000","0163103100000010"); // restrict acceptance to EMCAL loose
    cuts.AddCut("00003113","00204409327000008250400000","1111121053032230000","0163103100000010"); // restrict acceptance to EMCAL tight
  } else if (trainConfig == 10){ // PCM variations pi dEdx
    cuts.AddCut("00003113","00200009317300008250400000","1111121053032230000","0163103100000010"); // dEdx pi: 0: 0.4-3.5, -10: 3.5 ->
    cuts.AddCut("00003113","00200009327300008250400000","1111121053032230000","0163103100000010"); // dEdx pi: 1: 0.4-3.5, -10: 3.5 ->
    cuts.AddCut("00003113","00200009325000008250400000","1111121053032230000","0163103100000010"); // dEdx pi: 1: 0.3 ->
    cuts.AddCut("00003113","00200009320000008250400000","1111121053032230000","0163103100000010"); // dEdx pi: 1: 0.5 ->
  } else if (trainConfig == 11){ // PCM variations pi dEdx  
    cuts.AddCut("00003113","00200009327600008250400000","1111121053032230000","0163103100000010"); // dEdx pi: 1: 0.4-2, -10: 2. ->
    cuts.AddCut("00003113","00200009327400008250400000","1111121053032230000","0163103100000010"); // dEdx pi: 1: 0.4-3, -10: 3. ->
    cuts.AddCut("00003113","00200009315600008250400000","1111121053032230000","0163103100000010"); // dEdx pi: 0: 0.3-2, -10: 2. ->
    cuts.AddCut("00003113","00200009367400008250400000","1111121053032230000","0163103100000010"); // dEdx pi: 2: 0.4-3, 0.5: 3. ->
    cuts.AddCut("00003113","00200009347400008250400000","1111121053032230000","0163103100000010"); // dEdx pi: 3: 0.4-3, 1: 3. ->
  } else if (trainConfig == 12){ // PCM variations to close V0s  
    cuts.AddCut("00003113","00200009327000008250400000","1111121053032230000","0163103100000010"); // dEdx pi: 1: 0.4-2, -10: 2. ->
    cuts.AddCut("00003113","00200009327000008250401000","1111121053032230000","0163103100000010"); // dEdx pi: 1: 0.4-3, -10: 3. ->
    cuts.AddCut("00003113","00200009327000008250402000","1111121053032230000","0163103100000010"); // dEdx pi: 0: 0.3-2, -10: 2. ->
    cuts.AddCut("00003113","00200009327000008250403000","1111121053032230000","0163103100000010"); // dEdx pi: 2: 0.4-3, 0.5: 3. ->
  } else if (trainConfig == 13){ // EMCAL clusters 2.76 TeV LHC11a, timing variation
    cuts.AddCut("00003113","00200009327000008250400000","1111131053032230000","0163103100000010");
    cuts.AddCut("00003113","00200009327000008250400000","1111132053032230000","0163103100000010");
    cuts.AddCut("00003113","00200009327000008250400000","1111111053032230000","0163103100000010");
    cuts.AddCut("00003113","00200009327000008250400000","1111112053032230000","0163103100000010");
  }else if (trainConfig == 14){  //LHC11a NonLinearity variations
    cuts.AddCut("00003113","00200009327000008250400000","1111100053032230000","0163103100000010"); // NonLinearity none
    cuts.AddCut("00003113","00200009327000008250400000","1111101053032230000","0163103100000010"); // NonLinearity kSDMv5
    cuts.AddCut("00003113","00200009327000008250400000","1111122053032230000","0163103100000010"); // NonLinearity LHC11a Calo
    cuts.AddCut("00003113","00200009327000008250400000","1111123053032230000","0163103100000010"); // NonLinearity LHC11a Calo
    
  // EMC1 variations  
  } else if (trainConfig == 15){ //EMCAL minEnergy variation
    cuts.AddCut("00051013","00200009327000008250400000","1111121053022230000","0163103100000010"); //0.6 GeV/c
    cuts.AddCut("00051013","00200009327000008250400000","1111121053032230000","0163103100000010"); //0.7 GeV/c default
    cuts.AddCut("00051013","00200009327000008250400000","1111121053042230000","0163103100000010"); //0.8 GeV/c
    cuts.AddCut("00051013","00200009327000008250400000","1111121053052230000","0163103100000010"); //0.9 GeV/c
  } else if (trainConfig == 16){ //EMCAL minNCells variation
    cuts.AddCut("00051013","00200009327000008250400000","1111121053031230000","0163103100000010"); //n cells >= 1
    cuts.AddCut("00051013","00200009327000008250400000","1111121053033230000","0163103100000010"); //n cells >= 3
    cuts.AddCut("00051013","00200009327000008250400000","1111121053032200000","0163103100000010"); //no max M02 cut
    cuts.AddCut("00051013","00200009327000008250400000","1113111053032230000","0163103100000010"); //only modules with TRD infront
    cuts.AddCut("00051013","00200009327000008250400000","1111211053032230000","0163103100000010"); //no modules with TRD infront
  } else if (trainConfig == 17){ // EMCAL track matching variations 
    cuts.AddCut("00051013","00200009327000008250400000","1111121051032230000","0163103100000010"); //
    cuts.AddCut("00051013","00200009327000008250400000","1111121052032230000","0163103100000010"); //
    cuts.AddCut("00051013","00200009327000008250400000","1111121053032230000","0163103100000010"); //
    cuts.AddCut("00051013","00200009327000008250400000","1111121054032230000","0163103100000010"); //
    cuts.AddCut("00051013","00200009327000008250400000","1111121055032230000","0163103100000010"); //
    cuts.AddCut("00051013","00200009327000008250400000","1111121056032230000","0163103100000010"); //
  } else if (trainConfig == 18){ // PCM variations
    cuts.AddCut("00051013","00200009227000008250400000","1111121053032230000","0163103100000010"); // dEdx e -3, 5
    cuts.AddCut("00051013","00200009127000008250400000","1111121053032230000","0163103100000010"); // dEdx e -5, 5
    cuts.AddCut("00051013","00200009357000008250400000","1111121053032230000","0163103100000010"); // dEdx pi 2
    cuts.AddCut("00051013","00200009317000008250400000","1111121053032230000","0163103100000010"); // dEdx pi 0
    cuts.AddCut("00051013","00200009387300008250400000","1111121053032230000","0163103100000010"); // dEdx pi 2 high 1 (> 3.5 GeV)
  } else if (trainConfig == 19){ // PCM variations
    cuts.AddCut("00051013","00200009327000009250400000","1111121053032230000","0163103100000010"); // qt 2D 0.03
    cuts.AddCut("00051013","00200009327000003250400000","1111121053032230000","0163103100000010"); // qt 1D 0.05
    cuts.AddCut("00051013","00200009327000002250400000","1111121053032230000","0163103100000010"); // qt 1D 0.07
    cuts.AddCut("00051013","00200049327000008250400000","1111121053032230000","0163103100000010"); // single pt > 0.075
    cuts.AddCut("00051013","00200019327000008250400000","1111121053032230000","0163103100000010"); // single pt > 0.1
  } else if (trainConfig == 20){ // PCM variations
    cuts.AddCut("00051013","00200009327000008850400000","1111121053032230000","0163103100000010"); // 2D psi pair chi2 var
    cuts.AddCut("00051013","00200009327000008260400000","1111121053032230000","0163103100000010"); // 2D psi pair chi2 var
    cuts.AddCut("00051013","00200009327000008860400000","1111121053032230000","0163103100000010"); // 2D psi pair chi2 var
    cuts.AddCut("00051013","00200009327000008280400000","1111121053032230000","0163103100000010"); // 2D psi pair chi2 var
    cuts.AddCut("00051013","00200009327000008880400000","1111121053032230000","0163103100000010"); // 2D psi pair chi2 var
  } else if (trainConfig == 21){ // PCM variations
    cuts.AddCut("00051013","00200006327000008250400000","1111121053032230000","0163103100000010"); // min TPC cl > 0.7
    cuts.AddCut("00051013","00200008327000008250400000","1111121053032230000","0163103100000010"); // min TPC cl > 0.35
    cuts.AddCut("00051013","00200009327000008250400000","1111121053032230000","0163106100000010"); // alpha < 0.8
    cuts.AddCut("00051013","00200009327000008250400000","1111121053032230000","0163105100000010"); // alpha < 0.75
  } else if (trainConfig == 22){ // PCM variations
    cuts.AddCut("00051013","00202209327000008250400000","1111121053032230000","0163103100000010"); // restrict acceptance to EMCAL loose
    cuts.AddCut("00051013","00204409327000008250400000","1111121053032230000","0163103100000010"); // restrict acceptance to EMCAL tight
  } else if (trainConfig == 23){ // EMCAL clusters 2.76 TeV LHC11a, with SDD (0), kEMC1 (1)
    cuts.AddCut("00051013","00200009327000008250400000","1111121053062230000","0163103100000010"); // min Energy cluster = 4.5 GeV
    cuts.AddCut("00051013","00200009327000008250400000","1111121053072230000","0163103100000010"); // min Energy cluster = 5.0 GeV
    cuts.AddCut("00051013","00200009327000008250400000","1111121053082230000","0163103100000010"); // min Energy cluster = 5.5 GeV
    cuts.AddCut("00051013","00200009327000008250400000","1111121053092230000","0163103100000010"); // min Energy cluster = 6.0 GeV
  } else if (trainConfig == 24){ // PCM variations pi dEdx
    cuts.AddCut("00051013","00200009317300008250400000","1111121053032230000","0163103100000010"); // dEdx pi: 0: 0.4-3.5, -10: 3.5 ->
    cuts.AddCut("00051013","00200009327300008250400000","1111121053032230000","0163103100000010"); // dEdx pi: 1: 0.4-3.5, -10: 3.5 ->
    cuts.AddCut("00051013","00200009325000008250400000","1111121053032230000","0163103100000010"); // dEdx pi: 1: 0.3 ->
    cuts.AddCut("00051013","00200009320000008250400000","1111121053032230000","0163103100000010"); // dEdx pi: 1: 0.5 ->
  } else if (trainConfig == 25){ // PCM variations pi dEdx  
    cuts.AddCut("00051013","00200009327600008250400000","1111121053032230000","0163103100000010"); // dEdx pi: 1: 0.4-2, -10: 2. ->
    cuts.AddCut("00051013","00200009327400008250400000","1111121053032230000","0163103100000010"); // dEdx pi: 1: 0.4-3, -10: 3. ->
    cuts.AddCut("00051013","00200009315600008250400000","1111121053032230000","0163103100000010"); // dEdx pi: 0: 0.3-2, -10: 2. ->
    cuts.AddCut("00051013","00200009367400008250400000","1111121053032230000","0163103100000010"); // dEdx pi: 2: 0.4-3, 0.5: 3. ->
    cuts.AddCut("00051013","00200009347400008250400000","1111121053032230000","0163103100000010"); // dEdx pi: 3: 0.4-3, 1: 3. ->
  } else if (trainConfig == 26){ // PCM variations to close V0s  
    cuts.AddCut("00051013","00200009327000008250400000","1111121053032230000","0163103100000010"); // dEdx pi: 1: 0.4-2, -10: 2. ->
    cuts.AddCut("00051013","00200009327000008250401000","1111121053032230000","0163103100000010"); // dEdx pi: 1: 0.4-3, -10: 3. ->
    cuts.AddCut("00051013","00200009327000008250402000","1111121053032230000","0163103100000010"); // dEdx pi: 0: 0.3-2, -10: 2. ->
    cuts.AddCut("00051013","00200009327000008250403000","1111121053032230000","0163103100000010"); // dEdx pi: 2: 0.4-3, 0.5: 3. ->
  } else if (trainConfig == 27){ // LHC11a NonLinearity variations
    cuts.AddCut("00051013","00200009327000008250400000","1111100053032230000","0163103100000010"); // NonLinearity none
    cuts.AddCut("00051013","00200009327000008250400000","1111101053032230000","0163103100000010"); // NonLinearity kSDMv5
    cuts.AddCut("00051013","00200009327000008250400000","1111122053032230000","0163103100000010"); // NonLinearity LHC11a Calo
    cuts.AddCut("00051013","00200009327000008250400000","1111123053032230000","0163103100000010"); // NonLinearity LHC11a Calo
  } else if (trainConfig == 28){ //LHC11a NonLinearity variations
    cuts.AddCut("00051013","00200009327000008250400000","1111131053032230000","0163103100000010");
    cuts.AddCut("00051013","00200009327000008250400000","1111132053032230000","0163103100000010");
    cuts.AddCut("00051013","00200009327000008250400000","1111111053032230000","0163103100000010");
    cuts.AddCut("00051013","00200009327000008250400000","1111112053032230000","0163103100000010");
    

  }else if (trainConfig == 30){  //LHC11a additional NonLinearity variations
    cuts.AddCut("00003113","00200009327000008250400000","1111113053032230000","0163103100000010"); // NonLinearity kTestBeamv2 + ConvCalo
    cuts.AddCut("00003113","00200009327000008250400000","1111114053032230000","0163103100000010"); // NonLinearity kTestBeamv2 + Calo
    cuts.AddCut("00003113","00200009327000008250400000","1111115053032230000","0163103100000010"); // NonLinearity LHC11a ConvCalo kSDM
    cuts.AddCut("00003113","00200009327000008250400000","1111116053032230000","0163103100000010"); // NonLinearity LHC11a Calo kSDM

  } else if (trainConfig == 31){  // LHC12 without non linearity
    cuts.AddCut("00000113","00200009327000008250400000","1111100063032230000","0163103100000010"); // INT7
    cuts.AddCut("00052013","00200009327000008250400000","1111100063032230000","0163103100000010"); // EMC7
    cuts.AddCut("00081013","00200009327000008250400000","1111100063032230000","0163103100000010"); // EMCEG1,
  } else if (trainConfig == 32){  // LHC10 without non linearity
    cuts.AddCut("00000113","00200009327000008250400000","1111100013032230000","0163103100000010"); // MB
  } else if (trainConfig == 33){  // EMCal, all triggers without non linearity
    cuts.AddCut("00000113","00200009327000008250400000","1111100063032230000","0163103100000010"); // INT7
    cuts.AddCut("00052013","00200009327000008250400000","1111100063032230000","0163103100000010"); // EMC7
    cuts.AddCut("00083013","00200009327000008250400000","1111100063032230000","0163103100000010"); // EMCEG1,
    cuts.AddCut("00085013","00200009327000008250400000","1111100063032230000","0163103100000010"); // EMCEG2,
  //LHC11a EMCal no non linearity internally  
  } else if (trainConfig == 34){ 
    cuts.AddCut("00003113","00200009327000008250400000","1111100053032230000","0163103100000010"); // 700 MeV cluster min energy
    cuts.AddCut("00051013","00200009327000008250400000","1111100053032230000","0163103100000010"); // 700 MeV cluster min energy
    
  // LHC13g  
  } else if (trainConfig == 40){  // LHC13g without pileup for triggers
    cuts.AddCut("00000113","00200009327000008250400000","1111121063032230000","0163103100000010"); // INT7
    cuts.AddCut("00000013","00200009327000008250400000","1111121063032230000","0163103100000010"); // INT7
    cuts.AddCut("00052013","00200009327000008250400000","1111121063032230000","0163103100000010"); // EMC7
    cuts.AddCut("00083013","00200009327000008250400000","1111121063032230000","0163103100000010"); // EMCEG1,
    cuts.AddCut("00085013","00200009327000008250400000","1111121063032230000","0163103100000010"); // EMCEG2,
    
  // LHC13g new conv calo non lienarity with pileup
  } else if (trainConfig == 41){  // EMCal, all triggers
    cuts.AddCut("00000113","00200009327000008250400000","1111121063032230000","0163103100000010"); // INT7
    cuts.AddCut("00052113","00200009327000008250400000","1111121063032230000","0163103100000010"); // EMC7
    cuts.AddCut("00083113","00200009327000008250400000","1111121063032230000","0163103100000010"); // EMCEG1,
    cuts.AddCut("00085113","00200009327000008250400000","1111121063032230000","0163103100000010"); // EMCEG2,

  // LHC13g variations MB cut  
  } else if (trainConfig == 42){ //EMCAL minEnergy variation
    cuts.AddCut("00000113","00200009327000008250400000","1111121063022230000","0163103100000010"); //0.6 GeV/c
    cuts.AddCut("00000113","00200009327000008250400000","1111121063032230000","0163103100000010"); //0.7 GeV/c default
    cuts.AddCut("00000113","00200009327000008250400000","1111121063042230000","0163103100000010"); //0.8 GeV/c
    cuts.AddCut("00000113","00200009327000008250400000","1111121063052230000","0163103100000010"); //0.9 GeV/c
  } else if (trainConfig == 43){ //EMCAL minNCells variation
    cuts.AddCut("00000113","00200009327000008250400000","1111121063031230000","0163103100000010"); //n cells >= 1
    cuts.AddCut("00000113","00200009327000008250400000","1111121063033230000","0163103100000010"); //n cells >= 3
    cuts.AddCut("00000113","00200009327000008250400000","1111121063032200000","0163103100000010"); //no max M02 cut
    cuts.AddCut("00000113","00200009327000008250400000","1112111063032230000","0163103100000010"); //only modules with TRD infront
    cuts.AddCut("00000113","00200009327000008250400000","1111311063032230000","0163103100000010"); //no modules with TRD infront
  } else if (trainConfig == 44){ // EMCAL track matching variations 
    cuts.AddCut("00000113","00200009327000008250400000","1111121061032230000","0163103100000010"); //
    cuts.AddCut("00000113","00200009327000008250400000","1111121062032230000","0163103100000010"); //
    cuts.AddCut("00000113","00200009327000008250400000","1111121063032230000","0163103100000010"); //
    cuts.AddCut("00000113","00200009327000008250400000","1111121064032230000","0163103100000010"); //
    cuts.AddCut("00000113","00200009327000008250400000","1111121065032230000","0163103100000010"); //
    cuts.AddCut("00000113","00200009327000008250400000","1111121066032230000","0163103100000010"); //
  } else if (trainConfig == 45){ // PCM variations
    cuts.AddCut("00000113","00200009227000008250400000","1111121063032230000","0163103100000010"); // dEdx e -3, 5
    cuts.AddCut("00000113","00200009127000008250400000","1111121063032230000","0163103100000010"); // dEdx e -5, 5
    cuts.AddCut("00000113","00200009357000008250400000","1111121063032230000","0163103100000010"); // dEdx pi 2
    cuts.AddCut("00000113","00200009317000008250400000","1111121063032230000","0163103100000010"); // dEdx pi 0
    cuts.AddCut("00000113","00200009387300008250400000","1111121063032230000","0163103100000010"); // dEdx pi 2 high 1 (> 3.5 GeV)
  } else if (trainConfig == 46){ // PCM variations
    cuts.AddCut("00000113","00200009327000009250400000","1111121063032230000","0163103100000010"); // qt 2D 0.03
    cuts.AddCut("00000113","00200009327000003250400000","1111121063032230000","0163103100000010"); // qt 1D 0.05
    cuts.AddCut("00000113","00200009327000002250400000","1111121063032230000","0163103100000010"); // qt 1D 0.07
    cuts.AddCut("00000113","00200049327000008250400000","1111121063032230000","0163103100000010"); // single pt > 0.075
    cuts.AddCut("00000113","00200019327000008250400000","1111121063032230000","0163103100000010"); // single pt > 0.1
  } else if (trainConfig == 47){ // PCM variations
    cuts.AddCut("00000113","00200009327000008850400000","1111121063032230000","0163103100000010"); // 2D psi pair chi2 var
    cuts.AddCut("00000113","00200009327000008260400000","1111121063032230000","0163103100000010"); // 2D psi pair chi2 var
    cuts.AddCut("00000113","00200009327000008860400000","1111121063032230000","0163103100000010"); // 2D psi pair chi2 var
    cuts.AddCut("00000113","00200009327000008280400000","1111121063032230000","0163103100000010"); // 2D psi pair chi2 var
    cuts.AddCut("00000113","00200009327000008880400000","1111121063032230000","0163103100000010"); // 2D psi pair chi2 var
  } else if (trainConfig == 48){ // PCM variations
    cuts.AddCut("00000113","00200006327000008250400000","1111121063032230000","0163103100000010"); // min TPC cl > 0.7
    cuts.AddCut("00000113","00200008327000008250400000","1111121063032230000","0163103100000010"); // min TPC cl > 0.35
    cuts.AddCut("00000113","00200009327000008250400000","1111121063032230000","0163106100000010"); // alpha < 0.8
    cuts.AddCut("00000113","00200009327000008250400000","1111121063032230000","0163105100000010"); // alpha < 0.75
  } else if (trainConfig == 49){ // PCM variations
    cuts.AddCut("00000113","00202209327000008250400000","1111121063032230000","0163103100000010"); // restrict acceptance to EMCAL loose
    cuts.AddCut("00000113","00204409327000008250400000","1111121063032230000","0163103100000010"); // restrict acceptance to EMCAL tight
  } else if (trainConfig == 50){ // PCM variations pi dEdx
    cuts.AddCut("00000113","00200009317300008250400000","1111121063032230000","0163103100000010"); // dEdx pi: 0: 0.4-3.5, -10: 3.5 ->
    cuts.AddCut("00000113","00200009327300008250400000","1111121063032230000","0163103100000010"); // dEdx pi: 1: 0.4-3.5, -10: 3.5 ->
    cuts.AddCut("00000113","00200009325000008250400000","1111121063032230000","0163103100000010"); // dEdx pi: 1: 0.3 ->
    cuts.AddCut("00000113","00200009320000008250400000","1111121063032230000","0163103100000010"); // dEdx pi: 1: 0.5 ->
  } else if (trainConfig == 51){ // PCM variations pi dEdx  
    cuts.AddCut("00000113","00200009327600008250400000","1111121063032230000","0163103100000010"); // dEdx pi: 1: 0.4-2, -10: 2. ->
    cuts.AddCut("00000113","00200009327400008250400000","1111121063032230000","0163103100000010"); // dEdx pi: 1: 0.4-3, -10: 3. ->
    cuts.AddCut("00000113","00200009315600008250400000","1111121063032230000","0163103100000010"); // dEdx pi: 0: 0.3-2, -10: 2. ->
    cuts.AddCut("00000113","00200009367400008250400000","1111121063032230000","0163103100000010"); // dEdx pi: 2: 0.4-3, 0.5: 3. ->
    cuts.AddCut("00000113","00200009347400008250400000","1111121063032230000","0163103100000010"); // dEdx pi: 3: 0.4-3, 1: 3. ->
  } else if (trainConfig == 52){ // PCM variations to close V0s  
    cuts.AddCut("00000113","00200009327000008250400000","1111121063032230000","0163103100000010"); // dEdx pi: 1: 0.4-2, -10: 2. ->
    cuts.AddCut("00000113","00200009327000008250401000","1111121063032230000","0163103100000010"); // dEdx pi: 1: 0.4-3, -10: 3. ->
    cuts.AddCut("00000113","00200009327000008250402000","1111121063032230000","0163103100000010"); // dEdx pi: 0: 0.3-2, -10: 2. ->
    cuts.AddCut("00000113","00200009327000008250403000","1111121063032230000","0163103100000010"); // dEdx pi: 2: 0.4-3, 0.5: 3. ->
  } else if (trainConfig == 53){ // NonLinearity variations
    cuts.AddCut("00000113","00200009327000008250400000","1111131063032230000","0163103100000010"); // INT7
    cuts.AddCut("00000113","00200009327000008250400000","1111132063032230000","0163103100000010"); // INT7
    cuts.AddCut("00000113","00200009327000008250400000","1111111063032230000","0163103100000010"); // INT7
    cuts.AddCut("00000113","00200009327000008250400000","1111112063032230000","0163103100000010"); // INT7
  } else if (trainConfig == 54){  //LHC11a NonLinearity variations
    cuts.AddCut("00000113","00200009327000008250400000","1111100063032230000","0163103100000010"); // NonLinearity none
    cuts.AddCut("00000113","00200009327000008250400000","1111101063032230000","0163103100000010"); // NonLinearity kSDMv5
    cuts.AddCut("00000113","00200009327000008250400000","1111122063032230000","0163103100000010"); // NonLinearity LHC11a Calo
    cuts.AddCut("00000113","00200009327000008250400000","1111123063032230000","0163103100000010"); // NonLinearity LHC11a Calo
    
  // LHC13g variations EMC7 cut  
  } else if (trainConfig == 55){ //EMCAL minEnergy variation
    cuts.AddCut("00052013","00200009327000008250400000","1111121063022230000","0163103100000010"); //0.6 GeV/c
    cuts.AddCut("00052013","00200009327000008250400000","1111121063032230000","0163103100000010"); //0.7 GeV/c default
    cuts.AddCut("00052013","00200009327000008250400000","1111121063042230000","0163103100000010"); //0.8 GeV/c
    cuts.AddCut("00052013","00200009327000008250400000","1111121063052230000","0163103100000010"); //0.9 GeV/c
  } else if (trainConfig == 56){ //EMCAL minNCells variation
    cuts.AddCut("00052013","00200009327000008250400000","1111121063031230000","0163103100000010"); //n cells >= 1
    cuts.AddCut("00052013","00200009327000008250400000","1111121063033230000","0163103100000010"); //n cells >= 3
    cuts.AddCut("00052013","00200009327000008250400000","1111121063032200000","0163103100000010"); //no max M02 cut
    cuts.AddCut("00052013","00200009327000008250400000","1112111063032230000","0163103100000010"); //only modules with TRD infront
    cuts.AddCut("00052013","00200009327000008250400000","1111311063032230000","0163103100000010"); //no modules with TRD infront
  } else if (trainConfig == 57){ // EMCAL track matching variations 
    cuts.AddCut("00052013","00200009327000008250400000","1111121061032230000","0163103100000010"); //
    cuts.AddCut("00052013","00200009327000008250400000","1111121062032230000","0163103100000010"); //
    cuts.AddCut("00052013","00200009327000008250400000","1111121063032230000","0163103100000010"); //
    cuts.AddCut("00052013","00200009327000008250400000","1111121064032230000","0163103100000010"); //
    cuts.AddCut("00052013","00200009327000008250400000","1111121065032230000","0163103100000010"); //
    cuts.AddCut("00052013","00200009327000008250400000","1111121066032230000","0163103100000010"); //
  } else if (trainConfig == 58){ // PCM variations
    cuts.AddCut("00052013","00200009227000008250400000","1111121063032230000","0163103100000010"); // dEdx e -3, 5
    cuts.AddCut("00052013","00200009127000008250400000","1111121063032230000","0163103100000010"); // dEdx e -5, 5
    cuts.AddCut("00052013","00200009357000008250400000","1111121063032230000","0163103100000010"); // dEdx pi 2
    cuts.AddCut("00052013","00200009317000008250400000","1111121063032230000","0163103100000010"); // dEdx pi 0
    cuts.AddCut("00052013","00200009387300008250400000","1111121063032230000","0163103100000010"); // dEdx pi 2 high 1 (> 3.5 GeV)
  } else if (trainConfig == 59){ // PCM variations
    cuts.AddCut("00052013","00200009327000009250400000","1111121063032230000","0163103100000010"); // qt 2D 0.03
    cuts.AddCut("00052013","00200009327000003250400000","1111121063032230000","0163103100000010"); // qt 1D 0.05
    cuts.AddCut("00052013","00200009327000002250400000","1111121063032230000","0163103100000010"); // qt 1D 0.07
    cuts.AddCut("00052013","00200049327000008250400000","1111121063032230000","0163103100000010"); // single pt > 0.075
    cuts.AddCut("00052013","00200019327000008250400000","1111121063032230000","0163103100000010"); // single pt > 0.1
  } else if (trainConfig == 60){ // PCM variations
    cuts.AddCut("00052013","00200009327000008850400000","1111121063032230000","0163103100000010"); // 2D psi pair chi2 var
    cuts.AddCut("00052013","00200009327000008260400000","1111121063032230000","0163103100000010"); // 2D psi pair chi2 var
    cuts.AddCut("00052013","00200009327000008860400000","1111121063032230000","0163103100000010"); // 2D psi pair chi2 var
    cuts.AddCut("00052013","00200009327000008280400000","1111121063032230000","0163103100000010"); // 2D psi pair chi2 var
    cuts.AddCut("00052013","00200009327000008880400000","1111121063032230000","0163103100000010"); // 2D psi pair chi2 var
  } else if (trainConfig == 61){ // PCM variations
    cuts.AddCut("00052013","00200006327000008250400000","1111121063032230000","0163103100000010"); // min TPC cl > 0.7
    cuts.AddCut("00052013","00200008327000008250400000","1111121063032230000","0163103100000010"); // min TPC cl > 0.35
    cuts.AddCut("00052013","00200009327000008250400000","1111121063032230000","0163106100000010"); // alpha < 0.8
    cuts.AddCut("00052013","00200009327000008250400000","1111121063032230000","0163105100000010"); // alpha < 0.75
  } else if (trainConfig == 62){ // PCM variations
    cuts.AddCut("00052013","00202209327000008250400000","1111121063032230000","0163103100000010"); // restrict acceptance to EMCAL loose
    cuts.AddCut("00052013","00204409327000008250400000","1111121063032230000","0163103100000010"); // restrict acceptance to EMCAL tight
  } else if (trainConfig == 63){ // PCM variations pi dEdx
    cuts.AddCut("00052013","00200009317300008250400000","1111121063032230000","0163103100000010"); // dEdx pi: 0: 0.4-3.5, -10: 3.5 ->
    cuts.AddCut("00052013","00200009327300008250400000","1111121063032230000","0163103100000010"); // dEdx pi: 1: 0.4-3.5, -10: 3.5 ->
    cuts.AddCut("00052013","00200009325000008250400000","1111121063032230000","0163103100000010"); // dEdx pi: 1: 0.3 ->
    cuts.AddCut("00052013","00200009320000008250400000","1111121063032230000","0163103100000010"); // dEdx pi: 1: 0.5 ->
  } else if (trainConfig == 64){ // PCM variations pi dEdx  
    cuts.AddCut("00052013","00200009327600008250400000","1111121063032230000","0163103100000010"); // dEdx pi: 1: 0.4-2, -10: 2. ->
    cuts.AddCut("00052013","00200009327400008250400000","1111121063032230000","0163103100000010"); // dEdx pi: 1: 0.4-3, -10: 3. ->
    cuts.AddCut("00052013","00200009315600008250400000","1111121063032230000","0163103100000010"); // dEdx pi: 0: 0.3-2, -10: 2. ->
    cuts.AddCut("00052013","00200009367400008250400000","1111121063032230000","0163103100000010"); // dEdx pi: 2: 0.4-3, 0.5: 3. ->
    cuts.AddCut("00052013","00200009347400008250400000","1111121063032230000","0163103100000010"); // dEdx pi: 3: 0.4-3, 1: 3. ->
  } else if (trainConfig == 65){ // PCM variations to close V0s  
    cuts.AddCut("00052013","00200009327000008250400000","1111121063032230000","0163103100000010"); // dEdx pi: 1: 0.4-2, -10: 2. ->
    cuts.AddCut("00052013","00200009327000008250401000","1111121063032230000","0163103100000010"); // dEdx pi: 1: 0.4-3, -10: 3. ->
    cuts.AddCut("00052013","00200009327000008250402000","1111121063032230000","0163103100000010"); // dEdx pi: 0: 0.3-2, -10: 2. ->
    cuts.AddCut("00052013","00200009327000008250403000","1111121063032230000","0163103100000010"); // dEdx pi: 2: 0.4-3, 0.5: 3. ->
  } else if (trainConfig == 67){ // NonLinearity
    cuts.AddCut("00052013","00200009327000008250400000","1111131063032230000","0163103100000010"); // NonLinearity 
    cuts.AddCut("00052013","00200009327000008250400000","1111132063032230000","0163103100000010"); // NonLinearity 
    cuts.AddCut("00052013","00200009327000008250400000","1111111063032230000","0163103100000010"); // NonLinearity 
    cuts.AddCut("00052013","00200009327000008250400000","1111112063032230000","0163103100000010"); // NonLinearity 
  } else if (trainConfig == 68){  //LHC11a NonLinearity variations
    cuts.AddCut("00052013","00200009327000008250400000","1111100063032230000","0163103100000010"); // NonLinearity none
    cuts.AddCut("00052013","00200009327000008250400000","1111101063032230000","0163103100000010"); // NonLinearity kSDMv5
    cuts.AddCut("00052013","00200009327000008250400000","1111122063032230000","0163103100000010"); // NonLinearity LHC11a Calo
    cuts.AddCut("00052013","00200009327000008250400000","1111123063032230000","0163103100000010"); // NonLinearity LHC11a ConvCalo+TestBeamv2
    
    
  // LHC13g variations EG2 cut  
  } else if (trainConfig == 69){ //EMCAL minEnergy variation
    cuts.AddCut("00085013","00200009327000008250400000","1111121063022230000","0163103100000010"); //0.6 GeV/c
    cuts.AddCut("00085013","00200009327000008250400000","1111121063032230000","0163103100000010"); //0.7 GeV/c default
    cuts.AddCut("00085013","00200009327000008250400000","1111121063042230000","0163103100000010"); //0.8 GeV/c
    cuts.AddCut("00085013","00200009327000008250400000","1111121063052230000","0163103100000010"); //0.9 GeV/c
  } else if (trainConfig == 70){ //EMCAL minNCells variation
    cuts.AddCut("00085013","00200009327000008250400000","1111121063031230000","0163103100000010"); //n cells >= 1
    cuts.AddCut("00085013","00200009327000008250400000","1111121063033230000","0163103100000010"); //n cells >= 3
    cuts.AddCut("00085013","00200009327000008250400000","1111121063032200000","0163103100000010"); //no max M02 cut
    cuts.AddCut("00085013","00200009327000008250400000","1112111063032230000","0163103100000010"); //only modules with TRD infront
    cuts.AddCut("00085013","00200009327000008250400000","1111311063032230000","0163103100000010"); //no modules with TRD infront
  } else if (trainConfig == 71){ // EMCAL track matching variations 
    cuts.AddCut("00085013","00200009327000008250400000","1111121061032230000","0163103100000010"); //
    cuts.AddCut("00085013","00200009327000008250400000","1111121062032230000","0163103100000010"); //
    cuts.AddCut("00085013","00200009327000008250400000","1111121063032230000","0163103100000010"); //
    cuts.AddCut("00085013","00200009327000008250400000","1111121064032230000","0163103100000010"); //
    cuts.AddCut("00085013","00200009327000008250400000","1111121065032230000","0163103100000010"); //
    cuts.AddCut("00085013","00200009327000008250400000","1111121066032230000","0163103100000010"); //
  } else if (trainConfig == 72){ // PCM variations
    cuts.AddCut("00085013","00200009227000008250400000","1111121063032230000","0163103100000010"); // dEdx e -3, 5
    cuts.AddCut("00085013","00200009127000008250400000","1111121063032230000","0163103100000010"); // dEdx e -5, 5
    cuts.AddCut("00085013","00200009357000008250400000","1111121063032230000","0163103100000010"); // dEdx pi 2
    cuts.AddCut("00085013","00200009317000008250400000","1111121063032230000","0163103100000010"); // dEdx pi 0
    cuts.AddCut("00085013","00200009387300008250400000","1111121063032230000","0163103100000010"); // dEdx pi 2 high 1 (> 3.5 GeV)
  } else if (trainConfig == 73){ // PCM variations
    cuts.AddCut("00085013","00200009327000009250400000","1111121063032230000","0163103100000010"); // qt 2D 0.03
    cuts.AddCut("00085013","00200009327000003250400000","1111121063032230000","0163103100000010"); // qt 1D 0.05
    cuts.AddCut("00085013","00200009327000002250400000","1111121063032230000","0163103100000010"); // qt 1D 0.07
    cuts.AddCut("00085013","00200049327000008250400000","1111121063032230000","0163103100000010"); // single pt > 0.075
    cuts.AddCut("00085013","00200019327000008250400000","1111121063032230000","0163103100000010"); // single pt > 0.1
  } else if (trainConfig == 74){ // PCM variations
    cuts.AddCut("00085013","00200009327000008850400000","1111121063032230000","0163103100000010"); // 2D psi pair chi2 var
    cuts.AddCut("00085013","00200009327000008260400000","1111121063032230000","0163103100000010"); // 2D psi pair chi2 var
    cuts.AddCut("00085013","00200009327000008860400000","1111121063032230000","0163103100000010"); // 2D psi pair chi2 var
    cuts.AddCut("00085013","00200009327000008280400000","1111121063032230000","0163103100000010"); // 2D psi pair chi2 var
    cuts.AddCut("00085013","00200009327000008880400000","1111121063032230000","0163103100000010"); // 2D psi pair chi2 var
  } else if (trainConfig == 75){ // PCM variations
    cuts.AddCut("00085013","00200006327000008250400000","1111121063032230000","0163103100000010"); // min TPC cl > 0.7
    cuts.AddCut("00085013","00200008327000008250400000","1111121063032230000","0163103100000010"); // min TPC cl > 0.35
    cuts.AddCut("00085013","00200009327000008250400000","1111121063032230000","0163106100000010"); // alpha < 0.8
    cuts.AddCut("00085013","00200009327000008250400000","1111121063032230000","0163105100000010"); // alpha < 0.75
  } else if (trainConfig == 76){ // PCM variations
    cuts.AddCut("00085013","00202209327000008250400000","1111121063032230000","0163103100000010"); // restrict acceptance to EMCAL loose
    cuts.AddCut("00085013","00204409327000008250400000","1111121063032230000","0163103100000010"); // restrict acceptance to EMCAL tight
  } else if (trainConfig == 77){ // PCM variations pi dEdx
    cuts.AddCut("00085013","00200009317300008250400000","1111121063032230000","0163103100000010"); // dEdx pi: 0: 0.4-3.5, -10: 3.5 ->
    cuts.AddCut("00085013","00200009327300008250400000","1111121063032230000","0163103100000010"); // dEdx pi: 1: 0.4-3.5, -10: 3.5 ->
    cuts.AddCut("00085013","00200009325000008250400000","1111121063032230000","0163103100000010"); // dEdx pi: 1: 0.3 ->
    cuts.AddCut("00085013","00200009320000008250400000","1111121063032230000","0163103100000010"); // dEdx pi: 1: 0.5 ->
  } else if (trainConfig == 78){ // PCM variations pi dEdx  
    cuts.AddCut("00085013","00200009327600008250400000","1111121063032230000","0163103100000010"); // dEdx pi: 1: 0.4-2, -10: 2. ->
    cuts.AddCut("00085013","00200009327400008250400000","1111121063032230000","0163103100000010"); // dEdx pi: 1: 0.4-3, -10: 3. ->
    cuts.AddCut("00085013","00200009315600008250400000","1111121063032230000","0163103100000010"); // dEdx pi: 0: 0.3-2, -10: 2. ->
    cuts.AddCut("00085013","00200009367400008250400000","1111121063032230000","0163103100000010"); // dEdx pi: 2: 0.4-3, 0.5: 3. ->
    cuts.AddCut("00085013","00200009347400008250400000","1111121063032230000","0163103100000010"); // dEdx pi: 3: 0.4-3, 1: 3. ->
  } else if (trainConfig == 79){ // PCM variations to close V0s  
    cuts.AddCut("00085013","00200009327000008250400000","1111121063032230000","0163103100000010"); // dEdx pi: 1: 0.4-2, -10: 2. ->
    cuts.AddCut("00085013","00200009327000008250401000","1111121063032230000","0163103100000010"); // dEdx pi: 1: 0.4-3, -10: 3. ->
    cuts.AddCut("00085013","00200009327000008250402000","1111121063032230000","0163103100000010"); // dEdx pi: 0: 0.3-2, -10: 2. ->
    cuts.AddCut("00085013","00200009327000008250403000","1111121063032230000","0163103100000010"); // dEdx pi: 2: 0.4-3, 0.5: 3. ->
  } else if (trainConfig == 80){ // EMCAL clusters 2.76 TeV non lin
    cuts.AddCut("00085013","00200009327000008250400000","1111131063032230000","0163103100000010"); 
    cuts.AddCut("00085013","00200009327000008250400000","1111132063032230000","0163103100000010"); 
    cuts.AddCut("00085013","00200009327000008250400000","1111111063032230000","0163103100000010"); 
    cuts.AddCut("00085013","00200009327000008250400000","1111112063032230000","0163103100000010"); 
  } else if (trainConfig == 81){  //LHC11a NonLinearity variations
    cuts.AddCut("00085013","00200009327000008250400000","1111100063032230000","0163103100000010"); // NonLinearity none
    cuts.AddCut("00085013","00200009327000008250400000","1111101063032230000","0163103100000010"); // NonLinearity kSDMv5
    cuts.AddCut("00085013","00200009327000008250400000","1111122063032230000","0163103100000010"); // NonLinearity LHC11a Calo
    cuts.AddCut("00085013","00200009327000008250400000","1111123063032230000","0163103100000010"); // NonLinearity LHC11a ConvCalo+TestBeamv2
    

  // LHC13g variations EG1 cut
  } else if (trainConfig == 82){ //EMCAL minEnergy variation
    cuts.AddCut("00083013","00200009327000008250400000","1111121063022230000","0163103100000010"); //0.6 GeV/c
    cuts.AddCut("00083013","00200009327000008250400000","1111121063032230000","0163103100000010"); //0.7 GeV/c default
    cuts.AddCut("00083013","00200009327000008250400000","1111121063042230000","0163103100000010"); //0.8 GeV/c
    cuts.AddCut("00083013","00200009327000008250400000","1111121063052230000","0163103100000010"); //0.9 GeV/c
  } else if (trainConfig == 83){ //EMCAL minNCells variation
    cuts.AddCut("00083013","00200009327000008250400000","1111121063031230000","0163103100000010"); //n cells >= 1
    cuts.AddCut("00083013","00200009327000008250400000","1111121063033230000","0163103100000010"); //n cells >= 3
    cuts.AddCut("00083013","00200009327000008250400000","1111121063032200000","0163103100000010"); //no max M02 cut
    cuts.AddCut("00083013","00200009327000008250400000","1112111063032230000","0163103100000010"); //only modules with TRD infront
    cuts.AddCut("00083013","00200009327000008250400000","1111311063032230000","0163103100000010"); //no modules with TRD infront
  } else if (trainConfig == 84){ // EMCAL track matching variations 
    cuts.AddCut("00083013","00200009327000008250400000","1111121061032230000","0163103100000010"); //
    cuts.AddCut("00083013","00200009327000008250400000","1111121062032230000","0163103100000010"); //
    cuts.AddCut("00083013","00200009327000008250400000","1111121063032230000","0163103100000010"); //
    cuts.AddCut("00083013","00200009327000008250400000","1111121064032230000","0163103100000010"); //
    cuts.AddCut("00083013","00200009327000008250400000","1111121065032230000","0163103100000010"); //
    cuts.AddCut("00083013","00200009327000008250400000","1111121066032230000","0163103100000010"); //
  } else if (trainConfig == 85){ // PCM variations
    cuts.AddCut("00083013","00200009227000008250400000","1111121063032230000","0163103100000010"); // dEdx e -3, 5
    cuts.AddCut("00083013","00200009127000008250400000","1111121063032230000","0163103100000010"); // dEdx e -5, 5
    cuts.AddCut("00083013","00200009357000008250400000","1111121063032230000","0163103100000010"); // dEdx pi 2
    cuts.AddCut("00083013","00200009317000008250400000","1111121063032230000","0163103100000010"); // dEdx pi 0
    cuts.AddCut("00083013","00200009387300008250400000","1111121063032230000","0163103100000010"); // dEdx pi 2 high 1 (> 3.5 GeV)
  } else if (trainConfig == 86){ // PCM variations
    cuts.AddCut("00083013","00200009327000009250400000","1111121063032230000","0163103100000010"); // qt 2D 0.03
    cuts.AddCut("00083013","00200009327000003250400000","1111121063032230000","0163103100000010"); // qt 1D 0.05
    cuts.AddCut("00083013","00200009327000002250400000","1111121063032230000","0163103100000010"); // qt 1D 0.07
    cuts.AddCut("00083013","00200049327000008250400000","1111121063032230000","0163103100000010"); // single pt > 0.075
    cuts.AddCut("00083013","00200019327000008250400000","1111121063032230000","0163103100000010"); // single pt > 0.1
  } else if (trainConfig == 87){ // PCM variations
    cuts.AddCut("00083013","00200009327000008850400000","1111121063032230000","0163103100000010"); // 2D psi pair chi2 var
    cuts.AddCut("00083013","00200009327000008260400000","1111121063032230000","0163103100000010"); // 2D psi pair chi2 var
    cuts.AddCut("00083013","00200009327000008860400000","1111121063032230000","0163103100000010"); // 2D psi pair chi2 var
    cuts.AddCut("00083013","00200009327000008280400000","1111121063032230000","0163103100000010"); // 2D psi pair chi2 var
    cuts.AddCut("00083013","00200009327000008880400000","1111121063032230000","0163103100000010"); // 2D psi pair chi2 var
  } else if (trainConfig == 88){ // PCM variations
    cuts.AddCut("00083013","00200006327000008250400000","1111121063032230000","0163103100000010"); // min TPC cl > 0.7
    cuts.AddCut("00083013","00200008327000008250400000","1111121063032230000","0163103100000010"); // min TPC cl > 0.35
    cuts.AddCut("00083013","00200009327000008250400000","1111121063032230000","0163106100000010"); // alpha < 0.8
    cuts.AddCut("00083013","00200009327000008250400000","1111121063032230000","0163105100000010"); // alpha < 0.75
  } else if (trainConfig == 89){ // PCM variations
    cuts.AddCut("00083013","00202209327000008250400000","1111121063032230000","0163103100000010"); // restrict acceptance to EMCAL loose
    cuts.AddCut("00083013","00204409327000008250400000","1111121063032230000","0163103100000010"); // restrict acceptance to EMCAL tight
  } else if (trainConfig == 90){ // PCM variations pi dEdx
    cuts.AddCut("00083013","00200009317300008250400000","1111121063032230000","0163103100000010"); // dEdx pi: 0: 0.4-3.5, -10: 3.5 ->
    cuts.AddCut("00083013","00200009327300008250400000","1111121063032230000","0163103100000010"); // dEdx pi: 1: 0.4-3.5, -10: 3.5 ->
    cuts.AddCut("00083013","00200009325000008250400000","1111121063032230000","0163103100000010"); // dEdx pi: 1: 0.3 ->
    cuts.AddCut("00083013","00200009320000008250400000","1111121063032230000","0163103100000010"); // dEdx pi: 1: 0.5 ->
  } else if (trainConfig == 91){ // PCM variations pi dEdx  
    cuts.AddCut("00083013","00200009327600008250400000","1111121063032230000","0163103100000010"); // dEdx pi: 1: 0.4-2, -10: 2. ->
    cuts.AddCut("00083013","00200009327400008250400000","1111121063032230000","0163103100000010"); // dEdx pi: 1: 0.4-3, -10: 3. ->
    cuts.AddCut("00083013","00200009315600008250400000","1111121063032230000","0163103100000010"); // dEdx pi: 0: 0.3-2, -10: 2. ->
    cuts.AddCut("00083013","00200009367400008250400000","1111121063032230000","0163103100000010"); // dEdx pi: 2: 0.4-3, 0.5: 3. ->
    cuts.AddCut("00083013","00200009347400008250400000","1111121063032230000","0163103100000010"); // dEdx pi: 3: 0.4-3, 1: 3. ->
  } else if (trainConfig == 92){ // PCM variations to close V0s  
    cuts.AddCut("00083013","00200009327000008250400000","1111121063032230000","0163103100000010"); // dEdx pi: 1: 0.4-2, -10: 2. ->
    cuts.AddCut("00083013","00200009327000008250401000","1111121063032230000","0163103100000010"); // dEdx pi: 1: 0.4-3, -10: 3. ->
    cuts.AddCut("00083013","00200009327000008250402000","1111121063032230000","0163103100000010"); // dEdx pi: 0: 0.3-2, -10: 2. ->
    cuts.AddCut("00083013","00200009327000008250403000","1111121063032230000","0163103100000010"); // dEdx pi: 2: 0.4-3, 0.5: 3. ->
  } else if (trainConfig == 93){ // EMCAL NonLinearity variations
    cuts.AddCut("00083013","00200009327000008250400000","1111131063032230000","0163103100000010"); 
    cuts.AddCut("00083013","00200009327000008250400000","1111132063032230000","0163103100000010"); 
    cuts.AddCut("00083013","00200009327000008250400000","1111111063032230000","0163103100000010"); 
    cuts.AddCut("00083013","00200009327000008250400000","1111112063032230000","0163103100000010"); 
  } else if (trainConfig == 94){  //LHC11a NonLinearity variations
    cuts.AddCut("00083013","00200009327000008250400000","1111100063032230000","0163103100000010"); // NonLinearity none
    cuts.AddCut("00083013","00200009327000008250400000","1111101063032230000","0163103100000010"); // NonLinearity kSDMv5
    cuts.AddCut("00083013","00200009327000008250400000","1111122063032230000","0163103100000010"); // NonLinearity LHC11a Calo
    cuts.AddCut("00083013","00200009327000008250400000","1111123063032230000","0163103100000010"); // NonLinearity LHC11a ConvCalo+TestBeamv2

  //LHC13g
  } else if (trainConfig == 95){  // EMCAL clusters, EMCEGA triggers, track matching 0.035
    cuts.AddCut("00083013","00200009327000008250400000","1111121063032230000","0163103100000010"); // EMCEG1,
    cuts.AddCut("00085013","00200009327000008250400000","1111121063032230000","0163103100000010"); // EMCEG2,
  } else if (trainConfig == 96){  // EMCAL clusters, kEMC trigger, track matching 0.035
    cuts.AddCut("00000113","00200009327000008250400000","1111121063032230000","0163103100000010"); // INT7
    cuts.AddCut("00052013","00200009327000008250400000","1111121063032230000","0163103100000010"); // EMC7
  } else if (trainConfig == 97){  // EMCAL clusters, EMCEGA triggers, track matching 0.035
    cuts.AddCut("00000113","00200009327000008250400000","1111121063032200000","0163103100000010"); // INT7, max M02 off
    cuts.AddCut("00052013","00200009327000008250400000","1111121063032200000","0163103100000010"); // EMC7, max M02 off
    cuts.AddCut("00083013","00200009327000008250400000","1111121063032200000","0163103100000010"); // EMCEG1, max M02 off
    cuts.AddCut("00085013","00200009327000008250400000","1111121063032200000","0163103100000010"); // EMCEG2, max M02 off

  } else if (trainConfig == 98){ // MB - with multiplicity bins
    cuts.AddCut("00103113","00200009327000008250400000","1111121053032230000","0163103100000010"); // 0 -2
    cuts.AddCut("01203113","00200009327000008250400000","1111121053032230000","0163103100000010"); // 2 -5
    cuts.AddCut("02303113","00200009327000008250400000","1111121053032230000","0163103100000010"); // 5 -10
    cuts.AddCut("03403113","00200009327000008250400000","1111121053032230000","0163103100000010"); // 10 -30
    cuts.AddCut("04503113","00200009327000008250400000","1111121053032230000","0163103100000010"); // 30 -100
  } else if (trainConfig == 99){ // INT7 - with multiplicity bins
    cuts.AddCut("00100113","00200009327000008250400000","1111121063032230000","0163103100000010"); // 0 -2
    cuts.AddCut("01200113","00200009327000008250400000","1111121063032230000","0163103100000010"); // 2 -5
    cuts.AddCut("02300113","00200009327000008250400000","1111121063032230000","0163103100000010"); // 5 -10
    cuts.AddCut("03400113","00200009327000008250400000","1111121063032230000","0163103100000010"); // 10 -30
    cuts.AddCut("04500113","00200009327000008250400000","1111121063032230000","0163103100000010"); // 30 -100
    
  // ************************************* EMCAL cuts ****************************************************
  // LHC12
  } else if (trainConfig == 101){ // EMCAL clusters 8 TeV LHC12
    cuts.AddCut("00000113","00200009327000008250400000","1111111063032230000","0163103100000010"); // 700 MeV cluster min energy
  } else if (trainConfig == 102){ //EMCAL minEnergy variation
    cuts.AddCut("00000113","00200009327000008250400000","1111111063022230000","0163103100000010"); //0.6 GeV/c
    cuts.AddCut("00000113","00200009327000008250400000","1111111063032230000","0163103100000010"); //0.7 GeV/c default
    cuts.AddCut("00000113","00200009327000008250400000","1111111063042230000","0163103100000010"); //0.8 GeV/c
    cuts.AddCut("00000113","00200009327000008250400000","1111111063052230000","0163103100000010"); //0.9 GeV/c
  } else if (trainConfig == 103){ //EMCAL minNCells variation
    cuts.AddCut("00000113","00200009327000008250400000","1111111063031230000","0163103100000010"); //n cells >= 1
    cuts.AddCut("00000113","00200009327000008250400000","1111111063033230000","0163103100000010"); //n cells >= 3
    cuts.AddCut("00000113","00200009327000008250400000","1111111063032200000","0163103100000010"); //no max M02 cut
    cuts.AddCut("00000113","00200009327000008250400000","1112111063032230000","0163103100000010"); //only modules with TRD infront
    cuts.AddCut("00000113","00200009327000008250400000","1111311063032230000","0163103100000010"); //no modules with TRD infront
  } else if (trainConfig == 104){ // EMCAL track matching variations
    cuts.AddCut("00000113","00200009327000008250400000","1111111061032230000","0163103100000010"); // track matching variations
    cuts.AddCut("00000113","00200009327000008250400000","1111111062032230000","0163103100000010"); //
    cuts.AddCut("00000113","00200009327000008250400000","1111111063032230000","0163103100000010"); //
    cuts.AddCut("00000113","00200009327000008250400000","1111111064032230000","0163103100000010"); //
    cuts.AddCut("00000113","00200009327000008250400000","1111111065032230000","0163103100000010"); //
    cuts.AddCut("00000113","00200009327000008250400000","1111111066032230000","0163103100000010"); //
  } else if (trainConfig == 105){ // PCM variations
    cuts.AddCut("00000113","00200009227000008250400000","1111111063032230000","0163103100000010"); // dEdx e -3, 5
    cuts.AddCut("00000113","00200009127000008250400000","1111111063032230000","0163103100000010"); // dEdx e -5, 5
    cuts.AddCut("00000113","00200009357000008250400000","1111111063032230000","0163103100000010"); // dEdx pi 2
    cuts.AddCut("00000113","00200009317000008250400000","1111111063032230000","0163103100000010"); // dEdx pi 0
    cuts.AddCut("00000113","00200009387300008250400000","1111111063032230000","0163103100000010"); // dEdx pi 2 high 1 (> 3.5 GeV)
  } else if (trainConfig == 106){ // PCM variations
    cuts.AddCut("00000113","00200009327000009250400000","1111111063032230000","0163103100000010"); // qt 2D 0.03
    cuts.AddCut("00000113","00200009327000003250400000","1111111063032230000","0163103100000010"); // qt 1D 0.05
    cuts.AddCut("00000113","00200009327000002250400000","1111111063032230000","0163103100000010"); // qt 1D 0.07
    cuts.AddCut("00000113","00200049327000008250400000","1111111063032230000","0163103100000010"); // single pt > 0.075
    cuts.AddCut("00000113","00200019327000008250400000","1111111063032230000","0163103100000010"); // single pt > 0.1
  } else if (trainConfig == 107){ // PCM variations
    cuts.AddCut("00000113","00200009327000008850400000","1111111063032230000","0163103100000010"); // 2D psi pair chi2 var
    cuts.AddCut("00000113","00200009327000008260400000","1111111063032230000","0163103100000010"); // 2D psi pair chi2 var
    cuts.AddCut("00000113","00200009327000008860400000","1111111063032230000","0163103100000010"); // 2D psi pair chi2 var
    cuts.AddCut("00000113","00200009327000008280400000","1111111063032230000","0163103100000010"); // 2D psi pair chi2 var
    cuts.AddCut("00000113","00200009327000008880400000","1111111063032230000","0163103100000010"); // 2D psi pair chi2 var
  } else if (trainConfig == 108){ // PCM variations
    cuts.AddCut("00000113","00200006327000008250400000","1111111063032230000","0163103100000010"); // min TPC cl > 0.7
    cuts.AddCut("00000113","00200008327000008250400000","1111111063032230000","0163103100000010"); // min TPC cl > 0.35
    cuts.AddCut("00000113","00200009327000008250400000","1111111063032230000","0163106100000010"); // alpha < 0.8
    cuts.AddCut("00000113","00200009327000008250400000","1111111063032230000","0163105100000010"); // alpha < 0.75
  } else if (trainConfig == 109){ // PCM variations
    cuts.AddCut("00000113","00202209327000008250400000","1111111063032230000","0163103100000010"); // restrict acceptance to EMCAL loose
    cuts.AddCut("00000113","00204409327000008250400000","1111111063032230000","0163103100000010"); // restrict acceptance to EMCAL tight
  } else if (trainConfig == 110){ // Different NonLinearities
    cuts.AddCut("00000113","00200009327000008250400000","1111111063032230000","0163103100000010"); // NonLinearity LHC12 ConvCalo
    cuts.AddCut("00000113","00200009327000008250400000","1111112063032230000","0163103100000010"); // NonLinearity LHC12 Calo
    cuts.AddCut("00000113","00200009327000008250400000","1111121063032230000","0163103100000010"); // NonLinearity LHC12 ConvCalo MassRatioFits
    cuts.AddCut("00000113","00200009327000008250400000","1111122063032230000","0163103100000010"); // NonLinearity LHC12 Calo MassRatioFits
    cuts.AddCut("00000113","00200009327000008250400000","1111100063032230000","0163103100000010"); // NonLinearity none
  } else if (trainConfig == 111){  //Different NonLinearities part2
    cuts.AddCut("00000113","00200009327000008250400000","1111101063032230000","0163103100000010"); // NonLinearity kSDMv5
    cuts.AddCut("00000113","00200009327000008250400000","1111113063032230000","0163103100000010"); // NonLinearity LHC12 kTestBeamv2+ConvCalo
    cuts.AddCut("00000113","00200009327000008250400000","1111114063032230000","0163103100000010"); // NonLinearity LHC12 kTestBeamv2+Calo
  } else if (trainConfig == 112){ // Variations DistanceToBadChannel
    cuts.AddCut("00000113","00200009327000008250400000","1111111063032230000","0163103100000010"); //
    cuts.AddCut("00000113","00200009327000008250400000","1111111163032230000","0163103100000010"); //
    cuts.AddCut("00000113","00200009327000008250400000","1111111263032230000","0163103100000010"); //
    cuts.AddCut("00000113","00200009327000008250400000","1111111363032230000","0163103100000010"); //
    cuts.AddCut("00000113","00200009327000008250400000","1111111563032230000","0163103100000010"); //
    cuts.AddCut("00000113","00200009327000008250400000","1111111663032230000","0163103100000010"); //
  } else if (trainConfig == 113){  // EMCAL clusters, MB with minEnergy variation
    cuts.AddCut("00000113","00200009327000008250400000","1111111063032230000","0163103100000010"); // min Energy cluster = 0.7 GeV
    cuts.AddCut("00000113","00200009327000008250400000","1111111063062230000","0163103100000010"); // min Energy cluster = 4.5 GeV
    cuts.AddCut("00000113","00200009327000008250400000","1111111063092230000","0163103100000010"); // min Energy cluster = 6.0 GeV
  } else if (trainConfig == 114){ // EMCAL clusters, EMC7 with minEnergy variation
    cuts.AddCut("00052013","00200009327000008250400000","1111111063032230000","0163103100000010"); // min Energy cluster = 0.7 GeV
    cuts.AddCut("00052013","00200009327000008250400000","1111111063062230000","0163103100000010"); // min Energy cluster = 4.5 GeV
    cuts.AddCut("00052013","00200009327000008250400000","1111111063092230000","0163103100000010"); // min Energy cluster = 6.0 GeV
  } else if (trainConfig == 115){ // EMCAL clusters, EMCEGA with minEnergy variation
    cuts.AddCut("00081013","00200009327000008250400000","1111111063032230000","0163103100000010"); // min Energy cluster = 0.7 GeV
    cuts.AddCut("00081013","00200009327000008250400000","1111111063062230000","0163103100000010"); // min Energy cluster = 4.5 GeV
    cuts.AddCut("00081013","00200009327000008250400000","1111111063092230000","0163103100000010"); // min Energy cluster = 6.0 GeV
  } else if (trainConfig == 116){ // EMCAL clusters, EMCEJE with minEnergy variation
    cuts.AddCut("00091113","00200009327000008250400000","1111111063032230000","0163103100000010"); // min Energy cluster = 0.7 GeV
    cuts.AddCut("00091113","00200009327000008250400000","1111111063062230000","0163103100000010"); // min Energy cluster = 4.5 GeV
    cuts.AddCut("00091113","00200009327000008250400000","1111111063092230000","0163103100000010"); // min Energy cluster = 6.0 GeV
  } else if (trainConfig == 117){  // EMCAL clusters, EMC triggers (EMC7, EMCEGA)
    cuts.AddCut("00000113","00200009327000008250400000","1111111063032230000","0163103100000010"); // INT7
    cuts.AddCut("00052013","00200009327000008250400000","1111111063032230000","0163103100000010"); // EMC7
    cuts.AddCut("00081013","00200009327000008250400000","1111111063032230000","0163103100000010"); // EMCEGA
  } else if (trainConfig == 118){ // EMCAL clusters, timing variation
    cuts.AddCut("00000113","00200009327000008250400000","1111111053032230000","0163103100000010"); // time 50ns
    cuts.AddCut("00000113","00200009327000008250400000","1111111043032230000","0163103100000010"); // time 100ns
    cuts.AddCut("00000113","00200009327000008250400000","1111111033032230000","0163103100000010"); // time 200ns
    cuts.AddCut("00000113","00200009327000008250400000","1111111023032230000","0163103100000010"); // time 500ns
  } else if (trainConfig == 119){ // EMCAL clusters, timing variation
    cuts.AddCut("00000113","00200009327000008250400000","1111111063032230000","0163103100000010"); // time
    cuts.AddCut("00000113","00200009327000008250400000","1111111073032230000","0163103100000010"); // time
    cuts.AddCut("00000113","00200009327000008250400000","1111111083032230000","0163103100000010"); // time
    cuts.AddCut("00000113","00200009327000008250400000","1111111093032230000","0163103100000010"); // time
  } else if (trainConfig == 120){ // EMCAL clusters, MB for extendedQA
    cuts.AddCut("00000113","00200009327000008250400000","1111111063032230000","0163103100000010"); // time
  } else if (trainConfig == 121){ // EMCAL clusters kEMC for extQA
    cuts.AddCut("00052013","00200009327000008250400000","1111111063032230000","0163103100000010"); //
  } else if (trainConfig == 122){ // EMCAL clusters EMCEGA for extQA
    cuts.AddCut("00081013","00200009327000008250400000","1111111063032230000","0163103100000010"); //
  } else if (trainConfig == 123){ // EMC7 Different NonLinearities
    cuts.AddCut("00052013","00200009327000008250400000","1111100063032230000","0163103100000010"); // NonLinearity none
    cuts.AddCut("00052013","00200009327000008250400000","1111111063032230000","0163103100000010"); // NonLinearity LHC12 ConvCalo
    cuts.AddCut("00052013","00200009327000008250400000","1111112063032230000","0163103100000010"); // NonLinearity LHC12 Calo
    cuts.AddCut("00052013","00200009327000008250400000","1111121063032230000","0163103100000010"); // NonLinearity LHC12 ConvCalo MassRatioFits
    cuts.AddCut("00052013","00200009327000008250400000","1111122063032230000","0163103100000010"); // NonLinearity LHC12 Calo MassRatioFits
  } else if (trainConfig == 124){ // EGA Different NonLinearities
    cuts.AddCut("00081013","00200009327000008250400000","1111100063032230000","0163103100000010"); // NonLinearity none
    cuts.AddCut("00081013","00200009327000008250400000","1111111063032230000","0163103100000010"); // NonLinearity LHC12 ConvCalo
    cuts.AddCut("00081013","00200009327000008250400000","1111112063032230000","0163103100000010"); // NonLinearity LHC12 Calo
    cuts.AddCut("00081013","00200009327000008250400000","1111121063032230000","0163103100000010"); // NonLinearity LHC12 ConvCalo MassRatioFits
    cuts.AddCut("00081013","00200009327000008250400000","1111122063032230000","0163103100000010"); // NonLinearity LHC12 Calo MassRatioFits
  } else if (trainConfig == 126){ // PCM variations pi dEdx
    cuts.AddCut("00000113","00200009317300008250400000","1111111063032230000","0163103100000010"); // dEdx pi: 0: 0.4-3.5, -10: 3.5 ->
    cuts.AddCut("00000113","00200009327300008250400000","1111111063032230000","0163103100000010"); // dEdx pi: 1: 0.4-3.5, -10: 3.5 ->
    cuts.AddCut("00000113","00200009325000008250400000","1111111063032230000","0163103100000010"); // dEdx pi: 1: 0.3 ->
    cuts.AddCut("00000113","00200009320000008250400000","1111111063032230000","0163103100000010"); // dEdx pi: 1: 0.5 ->
  } else if (trainConfig == 127){ // PCM variations pi dEdx
    cuts.AddCut("00000113","00200009327600008250400000","1111111063032230000","0163103100000010"); // dEdx pi: 1: 0.4-2, -10: 2. ->
    cuts.AddCut("00000113","00200009327400008250400000","1111111063032230000","0163103100000010"); // dEdx pi: 1: 0.4-3, -10: 3. ->
    cuts.AddCut("00000113","00200009315600008250400000","1111111063032230000","0163103100000010"); // dEdx pi: 0: 0.3-2, -10: 2. ->
    cuts.AddCut("00000113","00200009367400008250400000","1111111063032230000","0163103100000010"); // dEdx pi: 2: 0.4-3, 0.5: 3. ->
    cuts.AddCut("00000113","00200009347400008250400000","1111111063032230000","0163103100000010"); // dEdx pi: 3: 0.4-3, 1: 3. ->
  } else if (trainConfig == 128){ // PCM variations to close V0s
    cuts.AddCut("00000113","00200009327000008250400000","1111111063032230000","0163103100000010"); //
    cuts.AddCut("00000113","00200009327000008250401000","1111111063032230000","0163103100000010"); //
    cuts.AddCut("00000113","00200009327000008250402000","1111111063032230000","0163103100000010"); //
    cuts.AddCut("00000113","00200009327000008250403000","1111111063032230000","0163103100000010"); //
  } else if (trainConfig == 129){ // EMCAL clusters No NonLinearity + Calo/ConvCalo kSDM
    cuts.AddCut("00000113","00200009327000008250400000","1111100063032230000","0163103100000010"); // NonLinearity none
    cuts.AddCut("00000113","00200009327000008250400000","1111115063032230000","0163103100000010"); // NonLinearity ConvCalo kSDM
    cuts.AddCut("00000113","00200009327000008250400000","1111116063032230000","0163103100000010"); // NonLinearity Calo kSDM
  } else if (trainConfig == 130){ // EMCAL clusters, MB INT8 for extendedQA
    cuts.AddCut("00011113","00200009327000008250400000","1111111063032230000","0163103100000010"); // INT8
  } else if (trainConfig == 131){ // EMCAL clusters kEMC8 for extQA
    cuts.AddCut("00053013","00200009327000008250400000","1111111063032230000","0163103100000010"); // EMC8
  } else if (trainConfig == 132){ // EMCAL clusters EMCEGA+INT8 for extQA
    cuts.AddCut("00082013","00200009327000008250400000","1111111063032230000","0163103100000010"); // EGA+INT8

  // ************************************* EMCAL cuts ****************************************************
  // 7 TeV
  } else if (trainConfig == 201){ // EMCAL clusters pp 7 TeV
    cuts.AddCut("00000113","00200009327000008250400000","1111111013032230000","0163103100000010"); // 1000ns timing cut, std NL
  } else if (trainConfig == 202){ // EMCAL clusters pp 7 TeV - NL variations
    cuts.AddCut("00000113","00200009327000008250400000","1111111013032230000","0163103100000010"); // NL ConvCalo
    cuts.AddCut("00000113","00200009327000008250400000","1111112013032230000","0163103100000010"); // NL Calo
    cuts.AddCut("00000113","00200009327000008250400000","1111113013032230000","0163103100000010"); // NL ConvCalo + TestBeamv3
    cuts.AddCut("00000113","00200009327000008250400000","1111114013032230000","0163103100000010"); // NL Calo + TestBeamv3
    cuts.AddCut("00000113","00200009327000008250400000","1111100013032230000","0163103100000010"); // NL off
  } else if (trainConfig == 203){ // EMCAL clusters pp 7 TeV - NL variations
    cuts.AddCut("00000113","00200009327000008250400000","1111101013032230000","0163103100000010"); // NL kSDMv5
    cuts.AddCut("00000113","00200009327000008250400000","1111102013032230000","0163103100000010"); // NL Pi0MCv3 + TestBeamv3
    cuts.AddCut("00000113","00200009327000008250400000","1111103013032230000","0163103100000010"); // NL Pi0MCv3 + TestBeamv2
    cuts.AddCut("00000113","00200009327000008250400000","1111111013032230000","0163103100000010"); // NL ConvCalo - std

  
  // ************************************* PHOS cuts ****************************************************
  // LHC11a  
  } else if (trainConfig == 301) { //PHOS clusters
    cuts.AddCut("00003113","00200009327000008250400000","2444400047033200000","0163103100000010");
    cuts.AddCut("00003113","00200009327000008250400000","2444400048033200000","0163103100000010");
    cuts.AddCut("00003113","00200009327000008250400000","2444400049033200000","0163103100000010");
    cuts.AddCut("00061113","00200009327000008250400000","2444400047033200000","0163103100000010");
    cuts.AddCut("00061113","00200009327000008250400000","2444400048033200000","0163103100000010");
    cuts.AddCut("00061113","00200009327000008250400000","2444400049033200000","0163103100000010");
  // LHC13g & LHC12x
  } else if (trainConfig == 302) { //PHOS clusters
    cuts.AddCut("00000113","00200009327000008250400000","2444400047033200000","0163103100000010");
    cuts.AddCut("00000113","00200009327000008250400000","2444400048033200000","0163103100000010");
    cuts.AddCut("00000113","00200009327000008250400000","2444400049033200000","0163103100000010");
    cuts.AddCut("00062113","00200009327000008250400000","2444400047033200000","0163103100000010");
    cuts.AddCut("00062113","00200009327000008250400000","2444400048033200000","0163103100000010");
    cuts.AddCut("00062113","00200009327000008250400000","2444400049033200000","0163103100000010");
  // LHC11a
  } else if (trainConfig == 303) { //PHOS clusters without and with added signals
    cuts.AddCut("00003113","00200009327000008250400000","2444400047033200000","0163103100000010");
    cuts.AddCut("00003123","00200009327000008250400000","2444400047033200000","0163103100000010");

  // LHC12
  } else if (trainConfig == 331){ // PHOS clusters 8 TeV LHC12
    cuts.AddCut("00000113","00200009327000008250400000","2444400078023200000","0163103100000010"); // 600 MeV cluster min energy, t<|30|ns
  } else if (trainConfig == 332){ // PHOS clusters 8 TeV LHC12
    cuts.AddCut("00062113","00200009327000008250400000","2444400078023200000","0163103100000010"); // 600 MeV cluster min energy, t<|30|ns
  } else if (trainConfig == 342){ // With/without Added Signals
    cuts.AddCut("00000113","00200009327000008250400000","2444400078023200000","0163103100000010"); //
    cuts.AddCut("00000123","00200009327000008250400000","2444400078023200000","0163103100000010"); //

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

  Int_t numberOfCuts = cuts.GetNCuts();

  TList *EventCutList = new TList();
  TList *ConvCutList = new TList();
  TList *ClusterCutList = new TList();
  TList *MesonCutList = new TList();

  TList *HeaderList = new TList();
  if (periodname.Contains("LHC12i3")){  
    TObjString *Header2 = new TObjString("BOX");
    HeaderList->Add(Header2);
  } else if (periodname.CompareTo("LHC14e2b")==0){
    TObjString *Header2 = new TObjString("pi0_1");
    HeaderList->Add(Header2);
    TObjString *Header3 = new TObjString("eta_2");
    HeaderList->Add(Header3);
  }  
  
  TString energy = "";
  TString mcName = "";
  TString mcNameAdd = "";
  if (periodname.Contains("WOSDD")){
    mcNameAdd = "_WOSDD";
  } else if (periodname.Contains("WSDD")){
    mcNameAdd = "_WSDD";
  }   
  if (periodname.Contains("LHC12i3")){
    energy = "2760GeV";
    mcName = "Pythia8_LHC12i3";
  } else if (periodname.Contains("LHC12f1a")){  
    energy = "2760GeV";
    mcName = "Pythia8_LHC12f1a";  
  } else if (periodname.Contains("LHC12f1b")){  
    energy = "2760GeV";
    mcName = "Phojet_LHC12f1b";      
  } else if (periodname.Contains("LHC14e2a")){  
    energy = "8TeV";
    mcName = "Pythia8_LHC14e2a";      
  } else if (periodname.Contains("LHC14e2b")){  
    energy = "8TeV";
    mcName = "Pythia8_LHC14e2b";        
  } else if (periodname.Contains("LHC14e2c")){    
    energy = "8TeV";
    mcName = "Phojet_LHC14e2c";          
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
    if (fAddedSignal && (periodname.Contains("LHC12i3") || periodname.CompareTo("LHC14e2b")==0)){
      mcInputNamePi0 = Form("Pi0_%s%s_addSig_%s", mcName.Data(), mcNameAdd.Data(), energy.Data() );
      mcInputNameEta = Form("Eta_%s%s_addSig_%s", mcName.Data(), mcNameAdd.Data(), energy.Data() );
    } else {
      mcInputNamePi0 = Form("Pi0_%s%s_%s", mcName.Data(), mcNameAdd.Data(), energy.Data() );
      mcInputNameEta = Form("Eta_%s%s_%s", mcName.Data(), mcNameAdd.Data(), energy.Data() );
    }  
    
    if (doParticleWeighting) analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile( kTRUE, kTRUE, kFALSE, fileNameInputForPartWeighting, 
                                                                                           mcInputNamePi0, mcInputNameEta, "",fitNamePi0,fitNameEta);

    TString dataInputMultHisto  = "";
    TString mcInputMultHisto    = "";
    TString triggerString   = cuts.GetEventCut(i);
    triggerString           = triggerString(3,2);
    if (triggerString.CompareTo("03")==0) 
      triggerString         = "00";

    dataInputMultHisto      = Form("%s_%s", periodNameAnchor.Data(), triggerString.Data());
    mcInputMultHisto        = Form("%s_%s", periodNameV0Reader.Data(), triggerString.Data());
   
    if (doMultiplicityWeighting){
      cout << "enableling mult weighting" << endl;
      analysisEventCuts[i]->SetUseWeightMultiplicityFromFile( kTRUE, fileNameInputForMultWeighing, dataInputMultHisto, mcInputMultHisto );
    }
              
    
    analysisEventCuts[i]->SetTriggerMimicking(enableTriggerMimicking);
    analysisEventCuts[i]->SetTriggerOverlapRejecion(enableTriggerOverlapRej);
    analysisEventCuts[i]->SetMaxFacPtHard(maxFacPtHard);
    analysisEventCuts[i]->InitializeCutsFromCutString((cuts.GetEventCut(i)).Data());
    analysisEventCuts[i]->SetV0ReaderName(V0ReaderName);
    EventCutList->Add(analysisEventCuts[i]);
    analysisEventCuts[i]->SetFillCutHistograms("",kFALSE);
    
    analysisCuts[i] = new AliConversionPhotonCuts();
    analysisCuts[i]->InitializeCutsFromCutString((cuts.GetPhotonCut(i)).Data());
    analysisCuts[i]->SetV0ReaderName(V0ReaderName);
    analysisCuts[i]->SetIsHeavyIon(isHeavyIon);
    ConvCutList->Add(analysisCuts[i]);
    analysisCuts[i]->SetFillCutHistograms("",kFALSE);
  
    analysisClusterCuts[i] = new AliCaloPhotonCuts((isMC==2));
    analysisClusterCuts[i]->InitializeCutsFromCutString((cuts.GetClusterCut(i)).Data());
    analysisClusterCuts[i]->SetV0ReaderName(V0ReaderName);
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
  if(doTreeConvGammaShape) task->SetDoTreeConvGammaShowerShape(kTRUE);
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
