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
//pp together with all supporting classes
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
void AddTask_GammaCalo_pp(  Int_t     trainConfig                   = 1,                            // change different set of cuts
                            Int_t     isMC                          = 0,                            // run MC
                            Int_t     enableQAMesonTask             = 0,                            // enable QA in AliAnalysisTaskGammaCalo
                            Int_t     enableQAClusterTask           = 0,                            // enable additional QA task
                            TString   fileNameInputForPartWeighting = "MCSpectraInput.root",        // path to file for weigting input
                            TString   cutnumberAODBranch            = "000000006008400001001500000",
                            TString   periodname                    = "LHC12f1x",                   // period name
                            Bool_t    doParticleWeighting           = kFALSE,                       // enables weighting
                            Bool_t    isUsingTHnSparse              = kTRUE,                        // enable or disable usage of THnSparses for background estimation
                            Int_t     enableExtMatchAndQA           = 0,                            // enable QA(3), extMatch+QA(2), extMatch(1), disabled (0)
                            Bool_t    enableTriggerMimicking        = kFALSE,                       // enable trigger mimicking
                            Bool_t    enableTriggerOverlapRej       = kFALSE,                       // enable trigger overlap rejection
                            Float_t   maxFacPtHard                  = 3.,                           // maximum factor between hardest jet and ptHard generated
                            TString   periodNameV0Reader            = "",
                            Bool_t    doMultiplicityWeighting       = kFALSE,                       //
                            TString   fileNameInputForMultWeighing  = "Multiplicity.root",          //  
                            TString   periodNameAnchor              = ""
                            
) {
  
  Int_t isHeavyIon = 0;
  
  // ================== GetAnalysisManager ===============================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error(Form("AddTask_GammaCalo_pp_%i",trainConfig), "No analysis manager found.");
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
  if( !(AliV0ReaderV1*)mgr->GetTask("V0ReaderV1") ){
    AliV0ReaderV1 *fV0ReaderV1 = new AliV0ReaderV1("V0ReaderV1");
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
  
  // ************************************* EMCAL cuts ****************************************************
  // LHC11a
  if (trainConfig == 1){ // EMCAL clusters 2.76 TeV LHC11a, with SDD (0), kEMC1 (1)
    cuts.AddCut("00003113","1111121053032220000","0163103100000050"); // 700 MeV cluster min energy
    cuts.AddCut("00051013","1111121053032220000","0163103100000050"); // 700 MeV cluster min energy
    cuts.AddCut("00003013","1111121053032220000","0163103100000050"); // 700 MeV cluster min energy
  } else if (trainConfig == 2){ //EMCAL minEnergy variation
    cuts.AddCut("00003113","1111121053022220000","0163103100000050"); //0.6 GeV/c
    cuts.AddCut("00003113","1111121053042220000","0163103100000050"); //0.8 GeV/c
    cuts.AddCut("00003113","1111121053052220000","0163103100000050"); //0.9 GeV/c
  } else if (trainConfig == 3){ //EMCAL minNCells variation
    cuts.AddCut("00003113","1111121053031220000","0163103100000050"); //n cells >= 1
    cuts.AddCut("00003113","1111121053033220000","0163103100000050"); //n cells >= 3
    cuts.AddCut("00003113","1111121053032000000","0163103100000050"); //no M02 cut
    cuts.AddCut("00003113","1113121053032220000","0163103100000050"); //only modules with TRD infront
    cuts.AddCut("00003113","1111221053032220000","0163103100000050"); //no modules with TRD infront
  } else if (trainConfig == 4 ){ // EMCAL clusters 2.76 TeV NonLinearity
    cuts.AddCut("00003113","1111122053032220000","0163103100000050"); // NonLinearity LHC11a Calo
    cuts.AddCut("00003113","1111101053032220000","0163103100000050"); // NonLinearity kSDMv5
    cuts.AddCut("00003113","1111100053032220000","0163103100000050"); // NonLinearity none
    cuts.AddCut("00003113","1111132053032220000","0163103100000050"); // 700 MeV cluster min energy
    cuts.AddCut("00003113","1111131053032220000","0163103100000050");
  } else if (trainConfig == 5){  // EMCAL clusters, MB (INT1) trigger
    cuts.AddCut("00003113","1111111053032220000","0163103100000050");
    cuts.AddCut("00003113","1111112053032220000","0163103100000050"); // MB
    cuts.AddCut("00003113","1111113053032220000","0163103100000050"); // NonLinearity kTestBeamv2 + ConvCalo
    cuts.AddCut("00003113","1111114053032220000","0163103100000050"); // NonLinearity kTestBeamv2 + Calo
    cuts.AddCut("00003113","1111115053032220000","0163103100000050"); // NonLinearity LHC11a kSDM ConvCalo
    cuts.AddCut("00003113","1111116053032220000","0163103100000050"); // NonLinearity LHC11a kSDM Calo
  } else if (trainConfig == 6 ){ // EMCAL clusters open angle variation
    cuts.AddCut("00003113","1111121053032220000","0163103100000050"); // min open angle - 0.0202
    cuts.AddCut("00003113","1111121053032220000","0163103100000030"); // min open angle - 0.01
    cuts.AddCut("00003113","1111121053032220000","0163103100000040"); // min open angle - 0.0152
    cuts.AddCut("00003113","1111121053032220000","0163103100000060"); // min open angle - 0.0404
  } else if (trainConfig == 7){  // EMCAL clusters, MB (INT1) trigger
    cuts.AddCut("00003113","1111121053031220000","0163103100000050"); // MB,                     NCells >=1
    cuts.AddCut("00003113","1111121053033220000","0163103100000050"); // MB,                     NCells >=3
    cuts.AddCut("00003113","1111121053032220000","0163103100000030"); // MB,                                                               0.01 opening
    cuts.AddCut("00003113","1111121053032220000","0163103100000040"); // MB,                                                               0.75 cell diagonal
    cuts.AddCut("00003113","1111121053032220000","0163103100000060"); // MB,                                                               2 cell diagonals
  } else if (trainConfig == 8){  // trackMatching variations
    cuts.AddCut("00003113","1111121051032220000","0163103100000050"); // MB
    cuts.AddCut("00003113","1111121052032220000","0163103100000050"); //
    cuts.AddCut("00003113","1111121053032220000","0163103100000050"); //
    cuts.AddCut("00003113","1111121054032220000","0163103100000050"); //
    cuts.AddCut("00003113","1111121055032220000","0163103100000050"); //
    cuts.AddCut("00003113","1111121056032220000","0163103100000050"); //
  } else if (trainConfig == 9){  // replay Jason
    cuts.AddCut("00003113","1111101050032220000","0163103100000050"); // my default with SDM w/o TM
    cuts.AddCut("00003113","1111101053032220000","0163103100000050"); // my default with SDM w/ TM
    cuts.AddCut("00003113","1111101020032220000","0163103100000050"); // timing to 500ns
    cuts.AddCut("00003113","1111101020032000000","0163103100000050"); // timing to 500ns, no M02 cuts == exactly Jasons cuts
    cuts.AddCut("00003113","1111101021032000000","0163103100000050"); // timing to 500ns, no M02, mild TM cut
  
    
  } else if (trainConfig == 20){  // min Energy EMC1
    cuts.AddCut("00051013","1111121053022220000","0163103100000050"); //0.6 GeV/c
    cuts.AddCut("00051013","1111121053042220000","0163103100000050"); //0.8 GeV/c
    cuts.AddCut("00051013","1111121053052220000","0163103100000050"); //0.9 GeV/c
  } else if (trainConfig == 21){ //EMCAL minNCells variation
    cuts.AddCut("00051013","1111121053031220000","0163103100000050"); //n cells >= 1
    cuts.AddCut("00051013","1111121053033220000","0163103100000050"); //n cells >= 3
    cuts.AddCut("00051013","1111121053032000000","0163103100000050"); //no M02 cut
    cuts.AddCut("00051013","1113121053032220000","0163103100000050"); //only modules with TRD infront
    cuts.AddCut("00051013","1111221053032220000","0163103100000050"); //no modules with TRD infront
  } else if (trainConfig == 22){ // EMCAL clusters 2.76 TeV NonLinearity
    cuts.AddCut("00051013","1111100053032220000","0163103100000050"); // NonLinearity none
    cuts.AddCut("00051013","1111101053032220000","0163103100000050"); // NonLinearity kSDMv5
    cuts.AddCut("00051013","1111122053032220000","0163103100000050"); // NonLinearity LHC11a Calo
    cuts.AddCut("00051013","1111131053032220000","0163103100000050");    
    cuts.AddCut("00051013","1111132053032220000","0163103100000050"); 
  } else if (trainConfig == 23){  // EMCAL clusters, Non Lin
    cuts.AddCut("00051013","1111111053032220000","0163103100000050");
    cuts.AddCut("00051013","1111112053032220000","0163103100000050"); //
    cuts.AddCut("00051013","1111113053032220000","0163103100000050"); // NonLinearity kTestBeamv2 + ConvCalo
    cuts.AddCut("00051013","1111114053032220000","0163103100000050"); // NonLinearity kTestBeamv2 + Calo
    cuts.AddCut("00051013","1111115053032220000","0163103100000050"); // NonLinearity LHC11a kSDM ConvCalo
    cuts.AddCut("00051013","1111116053032220000","0163103100000050"); // NonLinearity LHC11a kSDM Calo
  } else if (trainConfig == 24 ){ // EMCAL clusters open angle variation
    cuts.AddCut("00051013","1111121053032220000","0163103100000050"); // min open angle - 0.0202
    cuts.AddCut("00051013","1111121053032220000","0163103100000030"); // min open angle - 0.01
    cuts.AddCut("00051013","1111121053032220000","0163103100000040"); // min open angle - 0.0152
    cuts.AddCut("00051013","1111121053032220000","0163103100000060"); // min open angle - 0.0404
  } else if (trainConfig == 25){  // trackMatching variations
    cuts.AddCut("00051013","1111121051032220000","0163103100000050"); // EMC1
    cuts.AddCut("00051013","1111121052032220000","0163103100000050"); //
    cuts.AddCut("00051013","1111121053032220000","0163103100000050"); //
    cuts.AddCut("00051013","1111121054032220000","0163103100000050"); //
    cuts.AddCut("00051013","1111121055032220000","0163103100000050"); //
    cuts.AddCut("00051013","1111121056032220000","0163103100000050"); //
    
  // ************************************* Calibration configuration EMC ********************************
  } else if (trainConfig == 40){ // EMCAL clusters 2.76 TeV LHC11a, with SDD (0), kEMC1 (1) with TM
    cuts.AddCut("00003113","1111100053032220000","0163103100000050"); // MB
    cuts.AddCut("00051013","1111100053032220000","0163103100000050"); // EMC1
  } else if (trainConfig == 41){ // EMCAL clusters 2.76 TeV LHC11a, with SDD (0), kEMC1 (1) without TM
    cuts.AddCut("00003113","1111100050032220000","0163103100000050"); // 700 MeV cluster min energy
    cuts.AddCut("00051013","1111100050032220000","0163103100000050"); // 700 MeV cluster min energy
  } else if (trainConfig == 42){  // EMCAL clusters 2.76TeV LHC13g with TM
    cuts.AddCut("00000113","1111100063032220000","0163103100000050"); 
    cuts.AddCut("00052013","1111100063032220000","0163103100000050"); // EMC7
    cuts.AddCut("00085013","1111100063032220000","0163103100000050"); // EG2
    cuts.AddCut("00083013","1111100063032220000","0163103100000050"); // EG1
  } else if (trainConfig == 43){  // EMCAL clusters 2.76TeV LHC13g without TM
    cuts.AddCut("00000113","1111100060032220000","0163103100000050"); 
    cuts.AddCut("00052013","1111100060032220000","0163103100000050"); // EMC7
    cuts.AddCut("00085013","1111100060032220000","0163103100000050"); // EG2
    cuts.AddCut("00083013","1111100060032220000","0163103100000050"); // EG1
  } else if (trainConfig == 44){   // EMCAL clusters 7TeV LHC10
    cuts.AddCut("00000113","1111100010032220000","0163103100000050"); // wo TM
    cuts.AddCut("00000113","1111100013032220000","0163103100000050"); // w TM
  } else if (trainConfig == 45){  // EMCAL clusters, 8TeV LHC12 with TM
    cuts.AddCut("00000113","1111100063032230000","0163103100000050"); 
    cuts.AddCut("00052013","1111100063032230000","0163103100000050"); // EMC7
    cuts.AddCut("00081013","1111100063032230000","0163103100000050"); // EMCEGA
  } else if (trainConfig == 46){  // EMCAL clusters, 8TeV LHC12 without TM
    cuts.AddCut("00000113","1111100060032230000","0163103100000050"); 
    cuts.AddCut("00052013","1111100060032230000","0163103100000050"); // EMC7
    cuts.AddCut("00081013","1111100060032230000","0163103100000050"); // EMCEGA
    
        
  // LHC13g  
  } else if (trainConfig == 60){  // EMCAL clusters, EMCEGA triggers
    cuts.AddCut("00000113","1111121063032220000","0163103100000050"); 
    cuts.AddCut("00000013","1111121063032220000","0163103100000050"); // without pile-up correction
    cuts.AddCut("00052013","1111121063032220000","0163103100000050"); // EMC7
    cuts.AddCut("00083013","1111121063032220000","0163103100000050"); // EMCEG1,
    cuts.AddCut("00085013","1111121063032220000","0163103100000050"); // EMCEG2,
    
  // Variations INT7 trigger
  } else if (trainConfig == 61){  // min Energy 
    cuts.AddCut("00000113","1111121063022220000","0163103100000050"); //0.6 GeV/c
    cuts.AddCut("00000113","1111121063042220000","0163103100000050"); //0.8 GeV/c
    cuts.AddCut("00000113","1111121063052220000","0163103100000050"); //0.9 GeV/c
  } else if (trainConfig == 62){  // EMCAL clusters trigger
    cuts.AddCut("00000113","1111121063031220000","0163103100000050"); //                      NCells >=1
    cuts.AddCut("00000113","1111121063033220000","0163103100000050"); //                      NCells >=3
    cuts.AddCut("00000113","1111121063032000000","0163103100000050"); //                                  no M02 cut
    cuts.AddCut("00000113","1112121063032220000","0163103100000050"); //only modules with TRD infront
    cuts.AddCut("00000113","1111321063032220000","0163103100000050"); //no modules with TRD infront    
    cuts.AddCut("00000113","1111121053032220000","0163103100000050"); //  50ns timing
  } else if (trainConfig == 63){  // EMCAL clusters 2.76TeV LHC13g
    cuts.AddCut("00000113","1111100063032220000","0163103100000050"); //  NonLinearity none
    cuts.AddCut("00000113","1111101063032220000","0163103100000050"); //  standard kSDMv5
    cuts.AddCut("00000113","1111122063032220000","0163103100000050"); //  NonLinearity LHC11a Calo
    cuts.AddCut("00000113","1111131063032220000","0163103100000050"); 
    cuts.AddCut("00000113","1111132063032220000","0163103100000050"); 
  } else if (trainConfig == 64){  // EMCAL clusters 2.76TeV LHC13g
    cuts.AddCut("00000113","1111111063032220000","0163103100000050"); 
    cuts.AddCut("00000113","1111112063032220000","0163103100000050"); 
    cuts.AddCut("00000113","1111113063032220000","0163103100000050"); //  NonLinearity TB+ConvCalo
    cuts.AddCut("00000113","1111114063032220000","0163103100000050"); //               TB+Calo
    cuts.AddCut("00000113","1111115063032220000","0163103100000050"); //               kPi0MC+ConvCalo (replay of Jasons with ConvCalo)
    cuts.AddCut("00000113","1111116063032220000","0163103100000050"); //               kPi0MC+Calo (replay of Jasons with Calo)
  } else if (trainConfig == 65 ){ // EMCAL clusters open angle variation
    cuts.AddCut("00000113","1111121063032220000","0163103100000050"); // min open angle - 0.0202
    cuts.AddCut("00000113","1111121063032220000","0163103100000030"); // min open angle - 0.01
    cuts.AddCut("00000113","1111121063032220000","0163103100000040"); // min open angle - 0.0152
    cuts.AddCut("00000113","1111121063032220000","0163103100000060"); // min open angle - 0.0404
  } else if (trainConfig == 66){  // trackMatching variations
    cuts.AddCut("00000113","1111121061032220000","0163103100000050"); 
    cuts.AddCut("00000113","1111121062032220000","0163103100000050"); //
    cuts.AddCut("00000113","1111121063032220000","0163103100000050"); //
    cuts.AddCut("00000113","1111121064032220000","0163103100000050"); //
    cuts.AddCut("00000113","1111121065032220000","0163103100000050"); //
    cuts.AddCut("00000113","1111121066032220000","0163103100000050"); //

  // Variations EMC7 trigger
  } else if (trainConfig == 70){  // min Energy EMC7
    cuts.AddCut("00052013","1111121063022220000","0163103100000050"); //0.6 GeV/c
    cuts.AddCut("00052013","1111121063042220000","0163103100000050"); //0.8 GeV/c
    cuts.AddCut("00052013","1111121063052220000","0163103100000050"); //0.9 GeV/c
  } else if (trainConfig == 71){  // EMCAL clusters
    cuts.AddCut("00052013","1111121063031220000","0163103100000050"); //                      NCells >=1
    cuts.AddCut("00052013","1111121063033220000","0163103100000050"); //                      NCells >=3
    cuts.AddCut("00052013","1111121063032000000","0163103100000050"); //                                  no M02 cut
    cuts.AddCut("00052013","1112121063032220000","0163103100000050"); //only modules with TRD infront
    cuts.AddCut("00052013","1111321063032220000","0163103100000050"); //no modules with TRD infront    
    cuts.AddCut("00052013","1111121053032220000","0163103100000050"); //  50ns timing
  } else if (trainConfig == 72){  // EMCAL clusters 2.76TeV LHC13g
    cuts.AddCut("00052013","1111100063032220000","0163103100000050"); //  NonLinearity none
    cuts.AddCut("00052013","1111101063032220000","0163103100000050"); //  standard kSDMv5
    cuts.AddCut("00052013","1111122063032220000","0163103100000050"); //  NonLinearity LHC11a Calo
    cuts.AddCut("00052013","1111131063032220000","0163103100000050"); 
    cuts.AddCut("00052013","1111132063032220000","0163103100000050"); 
  } else if (trainConfig == 73){  // EMCAL clusters 2.76TeV LHC13g
    cuts.AddCut("00052013","1111111063032220000","0163103100000050"); 
    cuts.AddCut("00052013","1111112063032220000","0163103100000050"); 
    cuts.AddCut("00052013","1111113063032220000","0163103100000050"); //  NonLinearity TB+ConvCalo
    cuts.AddCut("00052013","1111114063032220000","0163103100000050"); //               TB+Calo
    cuts.AddCut("00052013","1111115063032220000","0163103100000050"); //               kPi0MC+ConvCalo (replay of Jasons with ConvCalo)
    cuts.AddCut("00052013","1111116063032220000","0163103100000050"); //               kPi0MC+Calo (replay of Jasons with Calo)
  } else if (trainConfig == 74 ){ // EMCAL clusters open angle variation
    cuts.AddCut("00052013","1111121063032220000","0163103100000050"); // min open angle - 0.0202
    cuts.AddCut("00052013","1111121063032220000","0163103100000030"); // min open angle - 0.01
    cuts.AddCut("00052013","1111121063032220000","0163103100000040"); // min open angle - 0.0152
    cuts.AddCut("00052013","1111121063032220000","0163103100000060"); // min open angle - 0.0404
  } else if (trainConfig == 75){  // trackMatching variations
    cuts.AddCut("00052013","1111121061032220000","0163103100000050"); 
    cuts.AddCut("00052013","1111121062032220000","0163103100000050"); //
    cuts.AddCut("00052013","1111121063032220000","0163103100000050"); //
    cuts.AddCut("00052013","1111121064032220000","0163103100000050"); //
    cuts.AddCut("00052013","1111121065032220000","0163103100000050"); //
    cuts.AddCut("00052013","1111121066032220000","0163103100000050"); //
    
  // Variations EG2 trigger  
  } else if (trainConfig == 80){  // min Energy 
    cuts.AddCut("00085013","1111121063022220000","0163103100000050"); //0.6 GeV/c
    cuts.AddCut("00085013","1111121063042220000","0163103100000050"); //0.8 GeV/c
    cuts.AddCut("00085013","1111121063052220000","0163103100000050"); //0.9 GeV/c
  } else if (trainConfig == 81){  // EMCAL clusters
    cuts.AddCut("00085013","1111121063031220000","0163103100000050"); //                      NCells >=1
    cuts.AddCut("00085013","1111121063033220000","0163103100000050"); //                      NCells >=3
    cuts.AddCut("00085013","1111121063032000000","0163103100000050"); //                                  no M02 cut
    cuts.AddCut("00085013","1112121063032220000","0163103100000050"); //only modules with TRD infront
    cuts.AddCut("00085013","1111321063032220000","0163103100000050"); //no modules with TRD infront    
    cuts.AddCut("00085013","1111121053032220000","0163103100000050"); //  50ns timing
  } else if (trainConfig == 82){  // EMCAL clusters 2.76TeV LHC13g
    cuts.AddCut("00085013","1111100063032220000","0163103100000050"); //  NonLinearity none
    cuts.AddCut("00085013","1111101063032220000","0163103100000050"); //  standard kSDMv5
    cuts.AddCut("00085013","1111122063032220000","0163103100000050"); //  NonLinearity LHC11a Calo
    cuts.AddCut("00085013","1111131063032220000","0163103100000050"); 
    cuts.AddCut("00085013","1111132063032220000","0163103100000050"); 
  } else if (trainConfig == 83){  // EMCAL clusters 2.76TeV LHC13g
    cuts.AddCut("00085013","1111111063032220000","0163103100000050"); 
    cuts.AddCut("00085013","1111112063032220000","0163103100000050"); 
    cuts.AddCut("00085013","1111113063032220000","0163103100000050"); //  NonLinearity TB+ConvCalo
    cuts.AddCut("00085013","1111114063032220000","0163103100000050"); //               TB+Calo
    cuts.AddCut("00085013","1111115063032220000","0163103100000050"); //               kPi0MC+ConvCalo (replay of Jasons with ConvCalo)
    cuts.AddCut("00085013","1111116063032220000","0163103100000050"); //               kPi0MC+Calo (replay of Jasons with Calo)
  } else if (trainConfig == 84 ){ // EMCAL clusters open angle variation
    cuts.AddCut("00085013","1111121063032220000","0163103100000050"); // min open angle - 0.0202
    cuts.AddCut("00085013","1111121063032220000","0163103100000030"); // min open angle - 0.01
    cuts.AddCut("00085013","1111121063032220000","0163103100000040"); // min open angle - 0.0152
    cuts.AddCut("00085013","1111121063032220000","0163103100000060"); // min open angle - 0.0404
  } else if (trainConfig == 85){  // trackMatching variations
    cuts.AddCut("00085013","1111121061032220000","0163103100000050"); 
    cuts.AddCut("00085013","1111121062032220000","0163103100000050"); //
    cuts.AddCut("00085013","1111121063032220000","0163103100000050"); //
    cuts.AddCut("00085013","1111121064032220000","0163103100000050"); //
    cuts.AddCut("00085013","1111121065032220000","0163103100000050"); //
    cuts.AddCut("00085013","1111121066032220000","0163103100000050"); //

  // Variations EG1 trigger    
  } else if (trainConfig == 91){  // min Energy EMC1
    cuts.AddCut("00083013","1111121063022220000","0163103100000050"); //0.6 GeV/c
    cuts.AddCut("00083013","1111121063042220000","0163103100000050"); //0.8 GeV/c
    cuts.AddCut("00083013","1111121063052220000","0163103100000050"); //0.9 GeV/c
  } else if (trainConfig == 92){  // EMCAL clusters, INT7 trigger
    cuts.AddCut("00083013","1111121063031220000","0163103100000050"); //                      NCells >=1
    cuts.AddCut("00083013","1111121063033220000","0163103100000050"); //                      NCells >=3
    cuts.AddCut("00083013","1111121063032000000","0163103100000050"); //                                  no M02 cut
    cuts.AddCut("00083013","1112121063032220000","0163103100000050"); //only modules with TRD infront
    cuts.AddCut("00083013","1111321063032220000","0163103100000050"); //no modules with TRD infront    
    cuts.AddCut("00083013","1111121053032220000","0163103100000050"); //  50ns timing
  } else if (trainConfig == 93){  // EMCAL clusters 2.76TeV LHC13g
    cuts.AddCut("00083013","1111100063032220000","0163103100000050"); //  NonLinearity none
    cuts.AddCut("00083013","1111101063032220000","0163103100000050"); //  standard kSDMv5
    cuts.AddCut("00083013","1111122063032220000","0163103100000050"); //  NonLinearity LHC11a Calo
    cuts.AddCut("00083013","1111131063032220000","0163103100000050"); 
    cuts.AddCut("00083013","1111132063032220000","0163103100000050"); 
  } else if (trainConfig == 94){  // EMCAL clusters 2.76TeV LHC13g
    cuts.AddCut("00083013","1111111063032220000","0163103100000050"); 
    cuts.AddCut("00083013","1111112063032220000","0163103100000050"); 
    cuts.AddCut("00083013","1111113063032220000","0163103100000050"); //  NonLinearity TB+ConvCalo
    cuts.AddCut("00083013","1111114063032220000","0163103100000050"); //               TB+Calo
    cuts.AddCut("00083013","1111115063032220000","0163103100000050"); //               kPi0MC+ConvCalo (replay of Jasons with ConvCalo)
    cuts.AddCut("00083013","1111116063032220000","0163103100000050"); //               kPi0MC+Calo (replay of Jasons with Calo)
  } else if (trainConfig == 95 ){ // EMCAL clusters open angle variation
    cuts.AddCut("00083013","1111121063032220000","0163103100000050"); // min open angle - 0.0202
    cuts.AddCut("00083013","1111121063032220000","0163103100000030"); // min open angle - 0.01
    cuts.AddCut("00083013","1111121063032220000","0163103100000040"); // min open angle - 0.0152
    cuts.AddCut("00083013","1111121063032220000","0163103100000060"); // min open angle - 0.0404
  } else if (trainConfig == 96){  // trackMatching variations
    cuts.AddCut("00083013","1111121061032220000","0163103100000050"); 
    cuts.AddCut("00083013","1111121062032220000","0163103100000050"); //
    cuts.AddCut("00083013","1111121063032220000","0163103100000050"); //
    cuts.AddCut("00083013","1111121064032220000","0163103100000050"); //
    cuts.AddCut("00083013","1111121065032220000","0163103100000050"); //
    cuts.AddCut("00083013","1111121066032220000","0163103100000050"); //
        
// 8 TeV configs

    // here is the order of the cluster cut string
    // usually for EMCal we start with 11111: default values for            "ClusterType", "EtaMin", "EtaMax", "PhiMin", "PhiMax"
    // then two numbers for nonlinearity, e.g. 21: this is                  "NonLinearity1", "NonLinearity2"
    // Then some cuts on the clusters, e.g. 06003222: this is               "DistanceToBadChannel", "Timing", "TrackMatching", "ExoticCell", "MinEnergy", "MinNCells", "MinM02", "MaxM02"
    // finally some for now unused cuts, usually 0000: this is              "MinM20", "MaxM20", "MaximumDispersion", "NLM"

    //standard cut
  } else if (trainConfig == 101){ // EMCAL clusters pp 8 TeV 
    cuts.AddCut("00000113","1111111063032230000","0163103100000050"); // 700 MeV cluster min energy
    // 8 TeV variations
  } else if (trainConfig == 102){ // EMCAL clusters pp 8 TeV, timing variation
    cuts.AddCut("00000113","1111111053032230000","0163103100000050"); // time -50ns_50ns
    cuts.AddCut("00000113","1111111063032230000","0163103100000050"); // time -30ns_35ns - standard
    cuts.AddCut("00000113","1111111073032230000","0163103100000050"); // time -30ns_30ns
    cuts.AddCut("00000113","1111111083032230000","0163103100000050"); // time -20ns_30ns
  } else if (trainConfig == 103){ //EMCAL minEnergy variation
    cuts.AddCut("00000113","1111111063012230000","0163103100000050"); //0.5 GeV/c
    cuts.AddCut("00000113","1111111063022230000","0163103100000050"); //0.6 GeV/c
    cuts.AddCut("00000113","1111111063032230000","0163103100000050"); //0.7 GeV/c default
    cuts.AddCut("00000113","1111111063042230000","0163103100000050"); //0.8 GeV/c
    cuts.AddCut("00000113","1111111063052230000","0163103100000050"); //0.9 GeV/c
  } else if (trainConfig == 104){ //EMCAL minNCells, M02, with/without TRD variation
    cuts.AddCut("00000113","1111111063031230000","0163103100000050"); //n cells >= 1
    cuts.AddCut("00000113","1111111063033230000","0163103100000050"); //n cells >= 3
    cuts.AddCut("00000113","1111111063032000000","0163103100000050"); //no M02 cut
    cuts.AddCut("00000113","1113121063032230000","0163103100000050"); //only modules with TRD infront
    cuts.AddCut("00000113","1111221063032230000","0163103100000050"); //no modules with TRD infront
  } else if (trainConfig == 105){  // trackMatching variations
    cuts.AddCut("00000113","1111111061032230000","0163103100000050"); //
    cuts.AddCut("00000113","1111111062032230000","0163103100000050"); //
    cuts.AddCut("00000113","1111111063032230000","0163103100000050"); //
    cuts.AddCut("00000113","1111111064032230000","0163103100000050"); //
    cuts.AddCut("00000113","1111111065032230000","0163103100000050"); //
    cuts.AddCut("00000113","1111111066032230000","0163103100000050"); //

  } else if (trainConfig == 108){ // EMCAL clusters pp 8 TeV, Different DistanceToBadChannels
    cuts.AddCut("00000113","1111111063032230000","0163103100000050"); //
    cuts.AddCut("00000113","1111111163032230000","0163103100000050"); //
    cuts.AddCut("00000113","1111111263032230000","0163103100000050"); //
    cuts.AddCut("00000113","1111111363032230000","0163103100000050"); //
    cuts.AddCut("00000113","1111111563032230000","0163103100000050"); //
    cuts.AddCut("00000113","1111111663032230000","0163103100000050"); //
  } else if (trainConfig == 109){ // EMCAL clusters pp 8 TeV, Different NonLinearities
    cuts.AddCut("00000113","1111101063032230000","0163103100000050"); // NonLinearity kSDMv5
    cuts.AddCut("00000113","1111113063032230000","0163103100000050"); // NonLinearity kTestBeamv2 + LHC12 ConvCalo
    cuts.AddCut("00000113","1111114063032230000","0163103100000050"); // NonLinearity kTestBeamv2 + LHC12 Calo
  } else if (trainConfig == 110){ // EMCAL clusters pp 8 TeV, Different NonLinearities
    cuts.AddCut("00000113","1111111063032230000","0163103100000050"); // NonLinearity LHC12 ConvCalo
    cuts.AddCut("00000113","1111112063032230000","0163103100000050"); // NonLinearity LHC12 Calo
    cuts.AddCut("00000113","1111121063032230000","0163103100000050"); // NonLinearity LHC12 ConvCalo MassRatioFits
    cuts.AddCut("00000113","1111122063032230000","0163103100000050"); // NonLinearity LHC12 Calo MassRatioFits
    cuts.AddCut("00000113","1111100063032230000","0163103100000050"); // NonLinearity none

   // LHC12fa-i and MC
    // default with three cuts
  } else if (trainConfig == 111){  // EMCAL clusters, different triggers no NonLinearity
    cuts.AddCut("00000113","1111100063032230000","0163103100000050"); 
    cuts.AddCut("00052013","1111100063032230000","0163103100000050"); // EMC7
    cuts.AddCut("00081013","1111100063032230000","0163103100000050"); // EMCEG1,
    
    // all the cut variations
  } else if (trainConfig == 112){  // EMCAL clusters, EMCEG1 trigger
    cuts.AddCut("00081013","1111111063032230000","0163103100000050"); // EMCEGA, 400 MeV min energy, NCells >=2, M02 default cut
    cuts.AddCut("00081013","1111111063052230000","0163103100000050"); // EMCEGA, 600 MeV min energy
  } else if (trainConfig == 113){  // EMCAL clusters, EMCEG1 trigger
    cuts.AddCut("00081013","1111111063031230000","0163103100000050"); // EMCEGA,                     NCells >=1
    cuts.AddCut("00081013","1111111063033230000","0163103100000050"); // EMCEGA,                     NCells >=3
  } else if (trainConfig == 114){  // EMCAL clusters, EMCEG1 trigger
    cuts.AddCut("00081013","1111111063032000000","0163103100000050"); // EMCEGA,                                 no M02 cut
    cuts.AddCut("00081013","1111111023032230000","0163103100000050"); // EMCEGA,                                                 500ns timing
    cuts.AddCut("00085013","1111111043032230000","0163103100000050"); // EMCEGA,                                                 100ns timing
  } else if (trainConfig == 115){  // EMCAL clusters, INT7 trigger
    cuts.AddCut("00000113","1111111063032230000","0163103100000050"); //  700 MeV min energy, NCells >=2, M02 default cut
    cuts.AddCut("00000113","1111111063052230000","0163103100000050"); //  900 MeV min energy
  } else if (trainConfig == 116){  // EMCAL clusters, INT7 trigger
    cuts.AddCut("00000113","1111111063031230000","0163103100000050"); //                        NCells >=1
    cuts.AddCut("00000113","1111111063033230000","0163103100000050"); //                        NCells >=3
  } else if (trainConfig == 117){  // EMCAL clusters, INT7 trigger
    cuts.AddCut("00000113","1111111063032000000","0163103100000050"); //                                    no M02 cut
    cuts.AddCut("00000113","1111111023032230000","0163103100000050"); //                                                    500ns timing
    cuts.AddCut("00000113","1111111043032230000","0163103100000050"); //                                                    100ns timing
  } else if (trainConfig == 118){  // EMCAL clusters, EMC7 trigger
    cuts.AddCut("00052013","1111111063032230000","0163103100000050"); // EMC7, 700 MeV min energy, NCells >=2, M02 default cut
    cuts.AddCut("00052013","1111111063052230000","0163103100000050"); // EMC7, 900 MeV min energy
  } else if (trainConfig == 119){  // EMCAL clusters, EMC7 trigger
    cuts.AddCut("00052013","1111111063031230000","0163103100000050"); // EMC7,                     NCells >=1
    cuts.AddCut("00052013","1111111063033230000","0163103100000050"); // EMC7,                     NCells >=3
  } else if (trainConfig == 120){  // EMCAL clusters, EMC7 trigger
    cuts.AddCut("00052013","1111111063032000000","0163103100000050"); // EMC7,                                 no M02 cut
    cuts.AddCut("00052013","1111111023032230000","0163103100000050"); // EMC7,                                                 500ns timing
    cuts.AddCut("00052013","1111111043032230000","0163103100000050"); // EMC7,                                                 100ns timing

  }else if (trainConfig == 121){ // EMCAL clusters, different special triggers, different conv calo non lin
    cuts.AddCut("00000113","1111111063032230000","0163103100000050"); 
    cuts.AddCut("00052013","1111111063032230000","0163103100000050"); // EMC7
    cuts.AddCut("00081013","1111111063032230000","0163103100000050"); // EMCEG1,
  }else if (trainConfig == 122){ // EMCAL clusters, different special triggers, different kSDMv5
    cuts.AddCut("00000113","1111101063032230000","0163103100000050"); 
    cuts.AddCut("00052013","1111101063032230000","0163103100000050"); // EMC7
    cuts.AddCut("00081013","1111101063032230000","0163103100000050"); // EMCEG1,
  }else if (trainConfig == 123){ // EMCAL clusters, different special triggers, different conv calo non lin
    cuts.AddCut("00000113","1111112063032230000","0163103100000050"); 
    cuts.AddCut("00052013","1111112063032230000","0163103100000050"); // EMC7
    cuts.AddCut("00081013","1111112063032230000","0163103100000050"); // EMCEG1,
  } else if (trainConfig == 124){  // EMCAL clusters, EMCEG1 trigger
    cuts.AddCut("00000113","1111111063032230000","0163103100000060"); //  700 MeV min energy, NCells >=2, M02 default cut, wider distance cut
    cuts.AddCut("00052013","1111111063032230000","0163103100000060"); // EMC7, 700 MeV min energy, NCells >=2, M02 default cut, wider distance cut
    cuts.AddCut("00081013","1111111063032230000","0163103100000060"); // EMCEGA, 700 MeV min energy, NCells >=2, M02 default cut, wider distance cut

  } else if (trainConfig == 125){ // EMCAL clusters EMC7 pp 8 TeV, Different NonLinearities
    cuts.AddCut("00052013","1111100063032230000","0163103100000050"); // NonLinearity none
    cuts.AddCut("00052013","1111111063032230000","0163103100000050"); // NonLinearity LHC12 ConvCalo
    cuts.AddCut("00052013","1111112063032230000","0163103100000050"); // NonLinearity LHC12 Calo
    cuts.AddCut("00052013","1111121063032230000","0163103100000050"); // NonLinearity LHC12 ConvCalo MassRatioFits
    cuts.AddCut("00052013","1111122063032230000","0163103100000050"); // NonLinearity LHC12 Calo MassRatioFits
  } else if (trainConfig == 126){ // EMCAL clusters EMCEGA pp 8 TeV, Different NonLinearities
    cuts.AddCut("00081013","1111100063032230000","0163103100000050"); // NonLinearity none
    cuts.AddCut("00081013","1111111063032230000","0163103100000050"); // NonLinearity LHC12 ConvCalo
    cuts.AddCut("00081013","1111112063032230000","0163103100000050"); // NonLinearity LHC12 Calo
    cuts.AddCut("00081013","1111121063032230000","0163103100000050"); // NonLinearity LHC12 ConvCalo MassRatioFits
    cuts.AddCut("00081013","1111122063032230000","0163103100000050"); // NonLinearity LHC12 Calo MassRatioFits

// PHOS @ 8 TeV
  } else if (trainConfig == 181){ // PHOS clusters
    cuts.AddCut("00000113","2444400070023200000","0163003100900050");
  } else if (trainConfig == 182){ // PHOS clusters
    cuts.AddCut("00062113","2444400070023200000","0163003100900050");



  // pp multiplicity studies
  } else if (trainConfig == 198){ // MB - with multiplicity bins
    cuts.AddCut("00103113","1111121053032220000","0163103100000050"); // 0 -2
    cuts.AddCut("01203113","1111121053032220000","0163103100000050"); // 2 -5
    cuts.AddCut("02303113","1111121053032220000","0163103100000050"); // 5 -10
    cuts.AddCut("03403113","1111121053032220000","0163103100000050"); // 10 -30
    cuts.AddCut("04503113","1111121053032220000","0163103100000050"); // 30 -100
  } else if (trainConfig == 199){ // - with multiplicity bins
    cuts.AddCut("00100113","1111121063032220000","0163103100000050"); // 0 -2
    cuts.AddCut("01200113","1111121063032220000","0163103100000050"); // 2 -5
    cuts.AddCut("02300113","1111121063032220000","0163103100000050"); // 5 -10
    cuts.AddCut("03400113","1111121063032220000","0163103100000050"); // 10 -30
    cuts.AddCut("04500113","1111121063032220000","0163103100000050"); // 30 -100
    
    

    // 7 TeV
  } else if (trainConfig == 201){ // EMCAL clusters pp 7 TeV
    cuts.AddCut("00000113","1111111013032230000","0163103100000050"); // 1000ns timing cut, std NL
  } else if (trainConfig == 202){ // EMCAL clusters pp 7 TeV - NL variations
    cuts.AddCut("00000113","1111111013032230000","0163103100000050"); // NL ConvCalo
    cuts.AddCut("00000113","1111112013032230000","0163103100000050"); // NL Calo
    cuts.AddCut("00000113","1111113013032230000","0163103100000050"); // NL ConvCalo + TestBeamv3
    cuts.AddCut("00000113","1111114013032230000","0163103100000050"); // NL Calo + TestBeamv3
    cuts.AddCut("00000113","1111100013032230000","0163103100000050"); // NL off
  } else if (trainConfig == 203){ // EMCAL clusters pp 7 TeV - NL variations
    cuts.AddCut("00000113","1111101013032230000","0163103100000050"); // NL kSDMv5
    cuts.AddCut("00000113","1111102013032230000","0163103100000050"); // NL Pi0MCv3 + TestBeamv3
    cuts.AddCut("00000113","1111103013032230000","0163103100000050"); // NL Pi0MCv3 + TestBeamv2
    cuts.AddCut("00000113","1111111013032230000","0163103100000050"); // NL ConvCalo - std

  // ************************************* PHOS cuts ****************************************************
  } else if (trainConfig == 301) { //PHOS clusters
    cuts.AddCut("00003113","2444400040033200000","0163103100000050"); //pp LHC11a with SDD, PHOS
    cuts.AddCut("00000113","2444400040033200000","0163103100000050"); //pp LHC13g default MB
    cuts.AddCut("00061113","2444400040033200000","0163103100000050"); //pp LHC11a PHI1
    cuts.AddCut("00062113","2444400040033200000","0163103100000050"); //pp LHC11a PHI7
  } else if (trainConfig == 302){ // Validation PHOS
    cuts.AddCut("00003113","2444400040033200000","0163003100900050");
  } else if (trainConfig == 303){ // PHOS clusters, without and with added signals
    cuts.AddCut("00003113","2444400040033200000","0163003100900050");
    cuts.AddCut("00003123","2444400040033200000","0163003100900050");
    
    
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

  TList *EventCutList = new TList();
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
    
    if (doParticleWeighting) analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE, fileNameInputForPartWeighting, mcInputNamePi0, mcInputNameEta, "",fitNamePi0,fitNameEta);

    TString dataInputMultHisto  = "";
    TString mcInputMultHisto    = "";
    TString triggerString   = cuts.GetEventCut(i);
    triggerString           = triggerString(3,2);
    if (triggerString.CompareTo("03")==0) 
      triggerString         = "00";

    dataInputMultHisto      = Form("%s_%s", periodNameAnchor.Data(), triggerString.Data());
    mcInputMultHisto        = Form("%s_%s", periodNameV0Reader.Data(), triggerString.Data());
   
    if (doMultiplicityWeighting) analysisEventCuts[i]->SetUseWeightMultiplicityFromFile( kTRUE, fileNameInputForMultWeighing, dataInputMultHisto, mcInputMultHisto );

    
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
