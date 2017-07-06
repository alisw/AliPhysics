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
                            TString   fileNameInputForPartWeighting = "MCSpectraInput.root",        // path to file for weigting input / modified acceptance
                            TString   cutnumberAODBranch            = "000000006008400001001500000",
                            TString   periodname                    = "LHC12f1x",                   // period name
                            Bool_t    doParticleWeighting           = kFALSE,                       // enables weighting
                            Bool_t    isUsingTHnSparse              = kTRUE,                        // enable or disable usage of THnSparses for background estimation
                            Int_t     enableExtMatchAndQA           = 0,                            // disabled (0), extMatch (1), extQA_noCellQA (2), extMatch+extQA_noCellQA (3), extQA+cellQA (4), extMatch+extQA+cellQA (5)
                            Bool_t    enableTriggerMimicking        = kFALSE,                       // enable trigger mimicking
                            Bool_t    enableTriggerOverlapRej       = kFALSE,                       // enable trigger overlap rejection
                            Float_t   maxFacPtHard                  = 3.,                           // maximum factor between hardest jet and ptHard generated
                            TString   periodNameV0Reader            = "",                           // period Name for V0 reader
                            Bool_t    doMultiplicityWeighting       = kFALSE,                       // enable multiplicity weights
                            TString   fileNameInputForMultWeighing  = "Multiplicity.root",          // file for multiplicity weights
                            TString   periodNameAnchor              = "",                           // name of anchor period for weighting
                            Bool_t    enableSortingMCLabels         = kTRUE,                        // enable sorting for MC cluster labels
                            Bool_t    runLightOutput                = kFALSE,                       // switch to run light output (only essential histograms for afterburner)
                            TString   additionalTrainConfig         = "0"                           // additional counter for trainconfig
) {
  
  Bool_t doTreeEOverP = kFALSE; // switch to produce EOverP tree
  TH1S* histoAcc = 0x0;         // histo for modified acceptance
  Int_t localDebugFlag = 0;
  //parse additionalTrainConfig flag
  TObjArray *rAddConfigArr = additionalTrainConfig.Tokenize("_");
  if(rAddConfigArr->GetEntries()<1){cout << "ERROR: AddTask_GammaCalo_pp during parsing of additionalTrainConfig String '" << additionalTrainConfig.Data() << "'" << endl; return;}
  TObjString* rAdditionalTrainConfig;
  for(Int_t i = 0; i<rAddConfigArr->GetEntries() ; i++){
    if(i==0) rAdditionalTrainConfig = (TObjString*)rAddConfigArr->At(i);
    else{
      TObjString* temp = (TObjString*) rAddConfigArr->At(i);
      TString tempStr = temp->GetString();
      if(tempStr.CompareTo("EPCLUSTree") == 0){
        cout << "INFO: AddTask_GammaCalo_pp activating 'EPCLUSTree'" << endl;
        doTreeEOverP = kTRUE;
      }else if(tempStr.BeginsWith("MODIFYACC")){
        cout << "INFO: AddTask_GammaCalo_pp activating 'MODIFYacc'" << endl;
        TString tempType = tempStr;
        tempType.Replace(0,9,"");
        cout << "INFO: connecting to alien..." << endl;
        TGrid::Connect("alien://");
        cout << "done!" << endl;
        TFile *w = TFile::Open(fileNameInputForPartWeighting.Data());
        if(!w){cout << "ERROR: Could not open file: " << fileNameInputForPartWeighting.Data() << endl;return;}
        histoAcc = (TH1S*) w->Get(tempType.Data());
        if(!histoAcc) {cout << "ERROR: Could not find histo: " << tempType.Data() << endl;return;}
        cout << "found: " << histoAcc << endl;
      }else if(tempStr.BeginsWith("LOCALDEBUGFLAG")){
        cout << "INFO: AddTask_GammaCalo_pp activating 'LOCALDEBUGFLAG'" << endl;
        TString tempType = tempStr;
        tempType.Replace(0,14,"");
        localDebugFlag = tempType.Atoi();
        cout << "INFO: debug flag set to '" << localDebugFlag << "'" << endl;
      }
    }
  }
  TString sAdditionalTrainConfig = rAdditionalTrainConfig->GetString();
  if (sAdditionalTrainConfig.Atoi() > 0){
    trainConfig = trainConfig + sAdditionalTrainConfig.Atoi();
    cout << "INFO: AddTask_GammaCalo_pp running additionalTrainConfig '" << sAdditionalTrainConfig.Atoi() << "', train config: '" << trainConfig << "'" << endl;
  }

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
  TString V0ReaderName = Form("V0ReaderV1_%s_%s",cutnumberEvent.Data(),cutnumberPhoton.Data());
  if( !(AliV0ReaderV1*)mgr->GetTask(V0ReaderName.Data()) ){
    AliV0ReaderV1 *fV0ReaderV1 = new AliV0ReaderV1(V0ReaderName.Data());
    fV0ReaderV1->SetUseOwnXYZCalculation(kTRUE);
    fV0ReaderV1->SetCreateAODs(kFALSE);// AOD Output
    if (periodNameV0Reader.CompareTo("") != 0) fV0ReaderV1->SetPeriodName(periodNameV0Reader);
    fV0ReaderV1->SetUseAODConversionPhoton(kTRUE);
    if(trainConfig>=99 && trainConfig<200) fV0ReaderV1->SetImprovedPsiPair(0); //switch off for 8TeV as AODs are used for which improved psipair is not available

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
  AliAnalysisTaskGammaCalo *task=NULL;
  task= new AliAnalysisTaskGammaCalo(Form("GammaCalo_%i",trainConfig));
  task->SetIsHeavyIon(isHeavyIon);
  task->SetIsMC(isMC);
  task->SetV0ReaderName(V0ReaderName);
  task->SetLightOutput(runLightOutput);

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
    cuts.AddCut("00003113","1111121053032200000","0163103100000050"); //no max M02 cut
    cuts.AddCut("00003113","1111121053032250000","0163103100000050"); //M02 < 0.3
    cuts.AddCut("00003113","1111121053032260000","0163103100000050"); //M02 < 0.27
    cuts.AddCut("00003113","1113121053032220000","0163103100000050"); //only modules with TRD infront
    cuts.AddCut("00003113","1111221053032220000","0163103100000050"); //no modules with TRD infront
  } else if (trainConfig == 4 ){ // EMCAL clusters 2.76 TeV NonLinearity
    cuts.AddCut("00003113","1111122053032220000","0163103100000050"); // NonLinearity LHC11a Calo
    cuts.AddCut("00003113","1111101053032220000","0163103100000050"); // NonLinearity kSDMv5
    cuts.AddCut("00003113","1111100053032220000","0163103100000050"); // NonLinearity none
    cuts.AddCut("00003113","1111111053032220000","0163103100000050"); // 700 MeV cluster min energy
    cuts.AddCut("00003113","1111112053032220000","0163103100000050");
  } else if (trainConfig == 5){  // EMCAL clusters, MB (INT1) trigger
    cuts.AddCut("00003113","1111111053032220000","0163103100000050");
    cuts.AddCut("00003113","1111112053032220000","0163103100000050"); // MB
    cuts.AddCut("00003113","1111113053032220000","0163103100000050"); // NonLinearity kTestBeamv2 + ConvCalo
    cuts.AddCut("00003113","1111114053032220000","0163103100000050"); // NonLinearity kTestBeamv2 + Calo
    cuts.AddCut("00003113","1111115053032220000","0163103100000050"); // NonLinearity LHC11a kSDM ConvCalo
    cuts.AddCut("00003113","1111116053032220000","0163103100000050"); // NonLinearity LHC11a kSDM Calo
  } else if (trainConfig == 6 ){ // EMCAL clusters open angle variation
    cuts.AddCut("00003113","1111121053032220000","0163103100000050"); // min open angle - 0.0202 - std
    cuts.AddCut("00003113","1111121053032220000","0163103100000040"); // min open angle - 0.0152
    cuts.AddCut("00003113","1111121053032220000","0163103100000060"); // min open angle - 0.017
    cuts.AddCut("00003113","1111121053032220000","0163103100000070"); // min open angle - 0.016
    cuts.AddCut("00003113","1111121053032220000","0163103100000080"); // min open angle - 0.018
  } else if (trainConfig == 7){  // EMCAL clusters, MB (INT1) trigger
    cuts.AddCut("00003113","1111121053031220000","0163103100000050"); // MB,                     NCells >=1
    cuts.AddCut("00003113","1111121053033220000","0163103100000050"); // MB,                     NCells >=3
    cuts.AddCut("00003113","1111121053032220000","0163103100000040"); // MB,                                                               0.0152 opening
    cuts.AddCut("00003113","1111121053032220000","0163103100000070"); // MB,                                                               0.018 cell diagonals
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
  } else if (trainConfig == 10){  // trackMatching variations pt dependent
    cuts.AddCut("00003113","1111121057032220000","0163103100000050"); // MB
    cuts.AddCut("00003113","1111121058032220000","0163103100000050"); //
    cuts.AddCut("00003113","1111121059032220000","0163103100000050"); //
  } else if (trainConfig == 11){  // new default
    cuts.AddCut("00003113","1111121057032220000","0163103100000050"); // MB
  
    
  } else if (trainConfig == 20){  // min Energy EMC1
    cuts.AddCut("00051013","1111121053022220000","0163103100000050"); //0.6 GeV/c
    cuts.AddCut("00051013","1111121053042220000","0163103100000050"); //0.8 GeV/c
    cuts.AddCut("00051013","1111121053052220000","0163103100000050"); //0.9 GeV/c
  } else if (trainConfig == 21){ //EMCAL minNCells variation
    cuts.AddCut("00051013","1111121053031220000","0163103100000050"); //n cells >= 1
    cuts.AddCut("00051013","1111121053033220000","0163103100000050"); //n cells >= 3
    cuts.AddCut("00051013","1111121053032200000","0163103100000050"); //no max M02 cut
    cuts.AddCut("00051013","1111121053032250000","0163103100000050"); //M02 < 0.3
    cuts.AddCut("00051013","1111121053032260000","0163103100000050"); //M02 < 0.27
    cuts.AddCut("00051013","1113121053032220000","0163103100000050"); //only modules with TRD infront
    cuts.AddCut("00051013","1111221053032220000","0163103100000050"); //no modules with TRD infront
  } else if (trainConfig == 22){ // EMCAL clusters 2.76 TeV NonLinearity
    cuts.AddCut("00051013","1111100053032220000","0163103100000050"); // NonLinearity none
    cuts.AddCut("00051013","1111101053032220000","0163103100000050"); // NonLinearity kSDMv5
    cuts.AddCut("00051013","1111122053032220000","0163103100000050"); // NonLinearity LHC11a Calo
    cuts.AddCut("00051013","1111111053032220000","0163103100000050");
    cuts.AddCut("00051013","1111112053032220000","0163103100000050");
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
    cuts.AddCut("00051013","1111121053032220000","0163103100000060"); // min open angle - 0.017
  } else if (trainConfig == 25){  // trackMatching variations
    cuts.AddCut("00051013","1111121051032220000","0163103100000050"); // EMC1
    cuts.AddCut("00051013","1111121052032220000","0163103100000050"); //
    cuts.AddCut("00051013","1111121053032220000","0163103100000050"); //
    cuts.AddCut("00051013","1111121054032220000","0163103100000050"); //
    cuts.AddCut("00051013","1111121055032220000","0163103100000050"); //
    cuts.AddCut("00051013","1111121056032220000","0163103100000050"); //
  } else if (trainConfig == 26){  // trackMatching variations pt dependent
    cuts.AddCut("00051013","1111121057032220000","0163103100000050"); // MB
    cuts.AddCut("00051013","1111121058032220000","0163103100000050"); //
    cuts.AddCut("00051013","1111121059032220000","0163103100000050"); //

  // ************************************* Calibration configuration EMC ********************************
  } else if (trainConfig == 40){ // EMCAL clusters 2.76 TeV LHC11a, with SDD (0), kEMC1 (1) with TM
    cuts.AddCut("00003113","1111100053032220000","0163103100000050"); // MB
    cuts.AddCut("00051013","1111100053032220000","0163103100000050"); // EMC1
  } else if (trainConfig == 41){ // EMCAL clusters 2.76 TeV LHC11a, with SDD (0), kEMC1 (1) without TM
    cuts.AddCut("00003113","1111100050032220000","0163103100000050"); // 700 MeV cluster min energy
    cuts.AddCut("00051013","1111100050032220000","0163103100000050"); // 700 MeV cluster min energy
  } else if (trainConfig == 42){  // EMCAL clusters 2.76TeV LHC13g with TM
    cuts.AddCut("00010113","1111100063032220000","0163103100000050"); 
    cuts.AddCut("00052013","1111100063032220000","0163103100000050"); // EMC7
    cuts.AddCut("00085013","1111100063032220000","0163103100000050"); // EG2
    cuts.AddCut("00083013","1111100063032220000","0163103100000050"); // EG1
  } else if (trainConfig == 43){  // EMCAL clusters 2.76TeV LHC13g without TM
    cuts.AddCut("00010113","1111100060032220000","0163103100000050"); 
    cuts.AddCut("00052013","1111100060032220000","0163103100000050"); // EMC7
    cuts.AddCut("00085013","1111100060032220000","0163103100000050"); // EG2
    cuts.AddCut("00083013","1111100060032220000","0163103100000050"); // EG1
  } else if (trainConfig == 44){   // EMCAL clusters 7TeV LHC10
    cuts.AddCut("00000113","1111100010032220000","0163103100000050"); // wo TM
    cuts.AddCut("00000113","1111100013032220000","0163103100000050"); // w TM
  } else if (trainConfig == 45){  // EMCAL clusters, 8TeV LHC12 with TM
    cuts.AddCut("00010113","1111100067032230000","0163103100000060");
    cuts.AddCut("00052013","1111100067032230000","0163103100000060"); // EMC7
    cuts.AddCut("00081013","1111100067032230000","0163103100000060"); // EMCEGA
  } else if (trainConfig == 46){  // EMCAL clusters, 8TeV LHC12 without TM
    cuts.AddCut("00010113","1111100060032230000","0163103100000060");
    cuts.AddCut("00052013","1111100060032230000","0163103100000060"); // EMC7
    cuts.AddCut("00081013","1111100060032230000","0163103100000060"); // EMCEGA
    
        
  // LHC13g  
  } else if (trainConfig == 60){  // EMCAL clusters, EMCEGA triggers
    cuts.AddCut("00010113","1111121063032220000","0163103100000050"); 
    cuts.AddCut("00010013","1111121063032220000","0163103100000050"); // without pile-up correction
    cuts.AddCut("00052013","1111121063032220000","0163103100000050"); // EMC7
    cuts.AddCut("00083013","1111121063032220000","0163103100000050"); // EMCEG1,
    cuts.AddCut("00085013","1111121063032220000","0163103100000050"); // EMCEG2,
    
  // Variations INT7 trigger
  } else if (trainConfig == 61){  // min Energy 
    cuts.AddCut("00010113","1111121063022220000","0163103100000050"); //0.6 GeV/c
    cuts.AddCut("00010113","1111121063042220000","0163103100000050"); //0.8 GeV/c
    cuts.AddCut("00010113","1111121063052220000","0163103100000050"); //0.9 GeV/c
  } else if (trainConfig == 62){  // EMCAL clusters trigger
    cuts.AddCut("00010113","1111121063031220000","0163103100000050"); // NCells >=1
    cuts.AddCut("00010113","1111121063033220000","0163103100000050"); // NCells >=3
    cuts.AddCut("00010113","1111121063032200000","0163103100000050"); // no max M02 cut
    cuts.AddCut("00010113","1111121063032250000","0163103100000050"); // M02 < 0.3
    cuts.AddCut("00010113","1111121063032260000","0163103100000050"); // M02 < 0.27
    cuts.AddCut("00010113","1112121063032220000","0163103100000050"); // only modules with TRD infront
    cuts.AddCut("00010113","1111321063032220000","0163103100000050"); // no modules with TRD infront    
    cuts.AddCut("00010113","1111121053032220000","0163103100000050"); //  50ns timing
  } else if (trainConfig == 63){  // EMCAL clusters 2.76TeV LHC13g
    cuts.AddCut("00010113","1111100063032220000","0163103100000050"); //  NonLinearity none
    cuts.AddCut("00010113","1111101063032220000","0163103100000050"); //  standard kSDMv5
    cuts.AddCut("00010113","1111122063032220000","0163103100000050"); //  NonLinearity LHC11a Calo
    cuts.AddCut("00010113","1111111063032220000","0163103100000050");
    cuts.AddCut("00010113","1111112063032220000","0163103100000050");
  } else if (trainConfig == 64){  // EMCAL clusters 2.76TeV LHC13g
    cuts.AddCut("00010113","1111111063032220000","0163103100000050"); 
    cuts.AddCut("00010113","1111112063032220000","0163103100000050"); 
    cuts.AddCut("00010113","1111113063032220000","0163103100000050"); //  NonLinearity TB+ConvCalo
    cuts.AddCut("00010113","1111114063032220000","0163103100000050"); //               TB+Calo
    cuts.AddCut("00010113","1111115063032220000","0163103100000050"); //               kPi0MC+ConvCalo (replay of Jasons with ConvCalo)
    cuts.AddCut("00010113","1111116063032220000","0163103100000050"); //               kPi0MC+Calo (replay of Jasons with Calo)
  } else if (trainConfig == 65 ){ // EMCAL clusters open angle variation
    cuts.AddCut("00010113","1111121063032220000","0163103100000050"); // min open angle - 0.0202
    cuts.AddCut("00010113","1111121063032220000","0163103100000030"); // min open angle - 0.01
    cuts.AddCut("00010113","1111121063032220000","0163103100000040"); // min open angle - 0.0152
    cuts.AddCut("00010113","1111121063032220000","0163103100000060"); // min open angle - 0.017
  } else if (trainConfig == 66){  // trackMatching variations
    cuts.AddCut("00010113","1111121061032220000","0163103100000050"); 
    cuts.AddCut("00010113","1111121062032220000","0163103100000050"); //
    cuts.AddCut("00010113","1111121063032220000","0163103100000050"); //
    cuts.AddCut("00010113","1111121064032220000","0163103100000050"); //
    cuts.AddCut("00010113","1111121065032220000","0163103100000050"); //
    cuts.AddCut("00010113","1111121066032220000","0163103100000050"); //
  } else if (trainConfig == 67){  // trackMatching variations pt dept
    cuts.AddCut("00010113","1111121067032220000","0163103100000050"); 
    cuts.AddCut("00010113","1111121068032220000","0163103100000050"); //
    cuts.AddCut("00010113","1111121069032220000","0163103100000050"); //

  // Variations EMC7 trigger
  } else if (trainConfig == 70){  // min Energy EMC7
    cuts.AddCut("00052013","1111121063022220000","0163103100000050"); //0.6 GeV/c
    cuts.AddCut("00052013","1111121063042220000","0163103100000050"); //0.8 GeV/c
    cuts.AddCut("00052013","1111121063052220000","0163103100000050"); //0.9 GeV/c
  } else if (trainConfig == 71){  // EMCAL clusters
    cuts.AddCut("00052013","1111121063031220000","0163103100000050"); // NCells >=1
    cuts.AddCut("00052013","1111121063033220000","0163103100000050"); // NCells >=3
    cuts.AddCut("00052013","1111121063032200000","0163103100000050"); // no max M02 cut
    cuts.AddCut("00052013","1111121063032250000","0163103100000050"); // M02 < 0.3
    cuts.AddCut("00052013","1111121063032260000","0163103100000050"); // M02 < 0.27
    cuts.AddCut("00052013","1112121063032220000","0163103100000050"); // only modules with TRD infront
    cuts.AddCut("00052013","1111321063032220000","0163103100000050"); // no modules with TRD infront    
    cuts.AddCut("00052013","1111121053032220000","0163103100000050"); //  50ns timing
  } else if (trainConfig == 72){  // EMCAL clusters 2.76TeV LHC13g
    cuts.AddCut("00052013","1111100063032220000","0163103100000050"); //  NonLinearity none
    cuts.AddCut("00052013","1111101063032220000","0163103100000050"); //  standard kSDMv5
    cuts.AddCut("00052013","1111122063032220000","0163103100000050"); //  NonLinearity LHC11a Calo
    cuts.AddCut("00052013","1111111063032220000","0163103100000050");
    cuts.AddCut("00052013","1111112063032220000","0163103100000050");
  } else if (trainConfig == 73){  // EMCAL clusters 2.76TeV LHC13g
    cuts.AddCut("00052013","1111111063032220000","0163103100000050"); 
    cuts.AddCut("00052013","1111112063032220000","0163103100000050"); 
    cuts.AddCut("00052013","1111113063032220000","0163103100000050"); //  NonLinearity TB+ConvCalo
    cuts.AddCut("00052013","1111114063032220000","0163103100000050"); //               TB+Calo
    cuts.AddCut("00052013","1111115063032220000","0163103100000050"); //               kPi0MC+ConvCalo (replay of Jasons with ConvCalo)
    cuts.AddCut("00052013","1111116063032220000","0163103100000050"); //               kPi0MC+Calo (replay of Jasons with Calo)
  } else if (trainConfig == 74 ){ // EMCAL clusters open angle variation
    cuts.AddCut("00052013","1111121063032220000","0163103100000050"); // min open angle - 0.0202 - std
    cuts.AddCut("00052013","1111121063032220000","0163103100000040"); // min open angle - 0.0152
    cuts.AddCut("00052013","1111121063032220000","0163103100000060"); // min open angle - 0.017
    cuts.AddCut("00052013","1111121063032220000","0163103100000070"); // min open angle - 0.016
    cuts.AddCut("00052013","1111121063032220000","0163103100000080"); // min open angle - 0.018
  } else if (trainConfig == 75){  // trackMatching variations
    cuts.AddCut("00052013","1111121061032220000","0163103100000050"); 
    cuts.AddCut("00052013","1111121062032220000","0163103100000050"); //
    cuts.AddCut("00052013","1111121063032220000","0163103100000050"); //
    cuts.AddCut("00052013","1111121064032220000","0163103100000050"); //
    cuts.AddCut("00052013","1111121065032220000","0163103100000050"); //
    cuts.AddCut("00052013","1111121066032220000","0163103100000050"); //
  } else if (trainConfig == 76){  // trackMatching variations pt dependent
    cuts.AddCut("00052013","1111121067032220000","0163103100000050"); 
    cuts.AddCut("00052013","1111121068032220000","0163103100000050"); //
    cuts.AddCut("00052013","1111121069032220000","0163103100000050"); //
    
  // Variations EG2 trigger  
  } else if (trainConfig == 80){  // min Energy 
    cuts.AddCut("00085013","1111121063022220000","0163103100000050"); //0.6 GeV/c
    cuts.AddCut("00085013","1111121063042220000","0163103100000050"); //0.8 GeV/c
    cuts.AddCut("00085013","1111121063052220000","0163103100000050"); //0.9 GeV/c
  } else if (trainConfig == 81){  // EMCAL clusters
    cuts.AddCut("00085013","1111121063031220000","0163103100000050"); // NCells >=1
    cuts.AddCut("00085013","1111121063033220000","0163103100000050"); // NCells >=3
    cuts.AddCut("00085013","1111121063032200000","0163103100000050"); // no max M02 cut
    cuts.AddCut("00085013","1111121063032250000","0163103100000050"); // M02 < 0.3
    cuts.AddCut("00085013","1111121063032260000","0163103100000050"); // M02 < 0.27
    cuts.AddCut("00085013","1112121063032220000","0163103100000050"); // only modules with TRD infront
    cuts.AddCut("00085013","1111321063032220000","0163103100000050"); // no modules with TRD infront    
    cuts.AddCut("00085013","1111121053032220000","0163103100000050"); //  50ns timing
  } else if (trainConfig == 82){  // EMCAL clusters 2.76TeV LHC13g
    cuts.AddCut("00085013","1111100063032220000","0163103100000050"); //  NonLinearity none
    cuts.AddCut("00085013","1111101063032220000","0163103100000050"); //  standard kSDMv5
    cuts.AddCut("00085013","1111122063032220000","0163103100000050"); //  NonLinearity LHC11a Calo
    cuts.AddCut("00085013","1111111063032220000","0163103100000050");
    cuts.AddCut("00085013","1111112063032220000","0163103100000050");
  } else if (trainConfig == 83){  // EMCAL clusters 2.76TeV LHC13g
    cuts.AddCut("00085013","1111111063032220000","0163103100000050"); 
    cuts.AddCut("00085013","1111112063032220000","0163103100000050"); 
    cuts.AddCut("00085013","1111113063032220000","0163103100000050"); //  NonLinearity TB+ConvCalo
    cuts.AddCut("00085013","1111114063032220000","0163103100000050"); //               TB+Calo
    cuts.AddCut("00085013","1111115063032220000","0163103100000050"); //               kPi0MC+ConvCalo (replay of Jasons with ConvCalo)
    cuts.AddCut("00085013","1111116063032220000","0163103100000050"); //               kPi0MC+Calo (replay of Jasons with Calo)
  } else if (trainConfig == 84 ){ // EMCAL clusters open angle variation
    cuts.AddCut("00085013","1111121063032220000","0163103100000050"); // min open angle - 0.020 - std
    cuts.AddCut("00085013","1111121063032220000","0163103100000040"); // min open angle - 0.0152
    cuts.AddCut("00085013","1111121063032220000","0163103100000060"); // min open angle - 0.017
    cuts.AddCut("00085013","1111121063032220000","0163103100000070"); // min open angle - 0.016
    cuts.AddCut("00085013","1111121063032220000","0163103100000080"); // min open angle - 0.018
  } else if (trainConfig == 85){  // trackMatching variations
    cuts.AddCut("00085013","1111121061032220000","0163103100000050"); 
    cuts.AddCut("00085013","1111121062032220000","0163103100000050"); //
    cuts.AddCut("00085013","1111121063032220000","0163103100000050"); //
    cuts.AddCut("00085013","1111121064032220000","0163103100000050"); //
    cuts.AddCut("00085013","1111121065032220000","0163103100000050"); //
    cuts.AddCut("00085013","1111121066032220000","0163103100000050"); //
  } else if (trainConfig == 86){  // trackMatching variations pt dependent
    cuts.AddCut("00085013","1111121067032220000","0163103100000050"); 
    cuts.AddCut("00085013","1111121068032220000","0163103100000050"); //
    cuts.AddCut("00085013","1111121069032220000","0163103100000050"); //

  // Variations EG1 trigger    
  } else if (trainConfig == 91){  // min Energy EMC1
    cuts.AddCut("00083013","1111121063022220000","0163103100000050"); //0.6 GeV/c
    cuts.AddCut("00083013","1111121063042220000","0163103100000050"); //0.8 GeV/c
    cuts.AddCut("00083013","1111121063052220000","0163103100000050"); //0.9 GeV/c
  } else if (trainConfig == 92){  // EMCAL clusters, INT7 trigger
    cuts.AddCut("00083013","1111121063031220000","0163103100000050"); // NCells >=1
    cuts.AddCut("00083013","1111121063033220000","0163103100000050"); // NCells >=3
    cuts.AddCut("00083013","1111121063032200000","0163103100000050"); // no max M02 cut
    cuts.AddCut("00083013","1111121063032250000","0163103100000050"); // M02 < 0.3
    cuts.AddCut("00083013","1111121063032260000","0163103100000050"); // M02 < 0.27
    cuts.AddCut("00083013","1112121063032220000","0163103100000050"); // only modules with TRD infront
    cuts.AddCut("00083013","1111321063032220000","0163103100000050"); // no modules with TRD infront    
    cuts.AddCut("00083013","1111121053032220000","0163103100000050"); // 50ns timing
  } else if (trainConfig == 93){  // EMCAL clusters 2.76TeV LHC13g
    cuts.AddCut("00083013","1111100063032220000","0163103100000050"); //  NonLinearity none
    cuts.AddCut("00083013","1111101063032220000","0163103100000050"); //  standard kSDMv5
    cuts.AddCut("00083013","1111122063032220000","0163103100000050"); //  NonLinearity LHC11a Calo
    cuts.AddCut("00083013","1111111063032220000","0163103100000050");
    cuts.AddCut("00083013","1111112063032220000","0163103100000050");
  } else if (trainConfig == 94){  // EMCAL clusters 2.76TeV LHC13g
    cuts.AddCut("00083013","1111111063032220000","0163103100000050"); 
    cuts.AddCut("00083013","1111112063032220000","0163103100000050"); 
    cuts.AddCut("00083013","1111113063032220000","0163103100000050"); //  NonLinearity TB+ConvCalo
    cuts.AddCut("00083013","1111114063032220000","0163103100000050"); //               TB+Calo
    cuts.AddCut("00083013","1111115063032220000","0163103100000050"); //               kPi0MC+ConvCalo (replay of Jasons with ConvCalo)
    cuts.AddCut("00083013","1111116063032220000","0163103100000050"); //               kPi0MC+Calo (replay of Jasons with Calo)
  } else if (trainConfig == 95 ){ // EMCAL clusters open angle variation
    cuts.AddCut("00083013","1111121063032220000","0163103100000050"); // min open angle - 0.0202 - std
    cuts.AddCut("00083013","1111121063032220000","0163103100000040"); // min open angle - 0.0152
    cuts.AddCut("00083013","1111121063032220000","0163103100000060"); // min open angle - 0.017
    cuts.AddCut("00083013","1111121063032220000","0163103100000070"); // min open angle - 0.016
    cuts.AddCut("00083013","1111121063032220000","0163103100000080"); // min open angle - 0.018
  } else if (trainConfig == 96){  // trackMatching variations
    cuts.AddCut("00083013","1111121061032220000","0163103100000050"); 
    cuts.AddCut("00083013","1111121062032220000","0163103100000050"); //
    cuts.AddCut("00083013","1111121063032220000","0163103100000050"); //
    cuts.AddCut("00083013","1111121064032220000","0163103100000050"); //
    cuts.AddCut("00083013","1111121065032220000","0163103100000050"); //
    cuts.AddCut("00083013","1111121066032220000","0163103100000050"); //
  } else if (trainConfig == 97){  // trackMatching variations pt dependent
    cuts.AddCut("00083013","1111121067032220000","0163103100000050"); 
    cuts.AddCut("00083013","1111121068032220000","0163103100000050"); //
    cuts.AddCut("00083013","1111121069032220000","0163103100000050"); //
        
// 8 TeV configs

    // here is the order of the cluster cut string
    // usually for EMCal we start with 11111: default values for            "ClusterType", "EtaMin", "EtaMax", "PhiMin", "PhiMax"
    // then two numbers for nonlinearity, e.g. 21: this is                  "NonLinearity1", "NonLinearity2"
    // Then some cuts on the clusters, e.g. 06003222: this is               "DistanceToBadChannel", "Timing", "TrackMatching", "ExoticCell", "MinEnergy", "MinNCells", "MinM02", "MaxM02"
    // finally some for now unused cuts, usually 0000: this is              "MinM20", "MaxM20", "MaximumDispersion", "NLM"

    //std, but no opening angle cut
  } else if (trainConfig == 99){ // EMCAL clusters pp 8 TeV
    cuts.AddCut("00010113","1111111067032220000","0163103100000000"); // std
    cuts.AddCut("00052113","1111111067032220000","0163103100000000"); // std
    cuts.AddCut("00081113","1111111067032220000","0163103100000000"); // std
    //standard cuts
  } else if (trainConfig == 100){ // EMCAL clusters pp 8 TeV
    cuts.AddCut("00010113","1111111067032220000","01631031000000d0"); // std
    cuts.AddCut("00052113","1111111067032220000","01631031000000d0"); // std
    cuts.AddCut("00081113","1111111067032220000","01631031000000d0"); // std
  } else if (trainConfig == 101){ // EMCAL clusters pp 8 TeV 
    cuts.AddCut("00010113","1111111067032220000","01631031000000d0"); // std
    cuts.AddCut("00052113","1111111067032220000","01631031000000d0"); // std
    cuts.AddCut("00081113","1111111067032220000","01631031000000d0"); // std

    // 8 TeV variations
  } else if (trainConfig == 102){ // EMCAL clusters pp 8 TeV, timing variation
    cuts.AddCut("00010113","1111111057032220000","01631031000000d0"); // time -50ns_50ns
    cuts.AddCut("00010113","1111111067032220000","01631031000000d0"); // time -30ns_35ns - standard
    cuts.AddCut("00010113","1111111077032220000","01631031000000d0"); // time -30ns_30ns
    cuts.AddCut("00010113","1111111087032220000","01631031000000d0"); // time -20ns_30ns
  } else if (trainConfig == 103){ //EMCAL minEnergy variation
    cuts.AddCut("00010113","1111111067012220000","01631031000000d0"); //0.5 GeV/c
    cuts.AddCut("00010113","1111111067022220000","01631031000000d0"); //0.6 GeV/c
    cuts.AddCut("00010113","1111111067032220000","01631031000000d0"); //0.7 GeV/c default
    cuts.AddCut("00010113","1111111067042220000","01631031000000d0"); //0.8 GeV/c
    cuts.AddCut("00010113","1111111067052220000","01631031000000d0"); //0.9 GeV/c
  } else if (trainConfig == 104){ //EMCAL minNCells, M02, with/without TRD variation
    cuts.AddCut("00010113","1111111067031220000","01631031000000d0"); //n cells >= 1
    cuts.AddCut("00010113","1111111067033220000","01631031000000d0"); //n cells >= 3
    cuts.AddCut("00010113","1111111067032200000","01631031000000d0"); //no max M02 cut
    cuts.AddCut("00010113","1111111067032250000","01631031000000d0"); //M02 < 0.3
    cuts.AddCut("00010113","1111111067032260000","01631031000000d0"); //M02 < 0.27
    cuts.AddCut("00010113","1112111067032220000","01631031000000d0"); //only modules with TRD infront
    cuts.AddCut("00010113","1111311067032220000","01631031000000d0"); //no modules with TRD infront
  } else if (trainConfig == 105){  // trackMatching variations
    cuts.AddCut("00010113","1111111066032220000","01631031000000d0"); //
    cuts.AddCut("00010113","1111111067032220000","01631031000000d0"); //
    cuts.AddCut("00010113","1111111068032220000","01631031000000d0"); //
    cuts.AddCut("00010113","1111111069032220000","01631031000000d0"); //
    cuts.AddCut("00010113","1111111063032220000","01631031000000d0"); //
    cuts.AddCut("00010113","1111111060032220000","01631031000000d0"); //
  } else if (trainConfig == 106){ // EMCAL clusters pp 8 TeV, combining cluster within time window and without
    cuts.AddCut("00010113","1111111007032220000","01631031000000d0"); //
  } else if (trainConfig == 107){ // EMCAL clusters open angle variation
    cuts.AddCut("00010113","1111111067032220000","0163103100000060"); // min open angle - 0.017 - std
    cuts.AddCut("00010113","1111111067032220000","0163103100000040"); // min open angle - 0.0152
    cuts.AddCut("00010113","1111111067032220000","0163103100000070"); // min open angle - 0.016
    cuts.AddCut("00010113","1111111067032220000","0163103100000080"); // min open angle - 0.018
    cuts.AddCut("00010113","1111111067032220000","0163103100000090"); // min open angle - 0.019
  } else if (trainConfig == 108){ // EMCAL clusters pp 8 TeV, Different DistanceToBadChannels
    cuts.AddCut("00010113","1111111067032220000","01631031000000d0"); //
    cuts.AddCut("00010113","1111111167032220000","01631031000000d0"); //
    cuts.AddCut("00010113","1111111267032220000","01631031000000d0"); //
    cuts.AddCut("00010113","1111111367032220000","01631031000000d0"); //
    cuts.AddCut("00010113","1111111567032220000","01631031000000d0"); //
    cuts.AddCut("00010113","1111111667032220000","01631031000000d0"); //
  } else if (trainConfig == 109){ // EMCAL clusters pp 8 TeV, Different NonLinearities
    cuts.AddCut("00010113","1111101067032220000","01631031000000d0"); // NonLinearity kSDMv5
    cuts.AddCut("00010113","1111113067032220000","01631031000000d0"); // NonLinearity kTestBeamv2 + LHC12 ConvCalo
    cuts.AddCut("00010113","1111114067032220000","01631031000000d0"); // NonLinearity kTestBeamv2 + LHC12 Calo
  } else if (trainConfig == 110){ // EMCAL clusters pp 8 TeV, Different NonLinearities
    cuts.AddCut("00010113","1111111067032220000","01631031000000d0"); // NonLinearity LHC12 ConvCalo
    cuts.AddCut("00010113","1111112067032220000","01631031000000d0"); // NonLinearity LHC12 Calo
    cuts.AddCut("00010113","1111121067032220000","01631031000000d0"); // NonLinearity LHC12 ConvCalo MassRatioFits
    cuts.AddCut("00010113","1111122067032220000","01631031000000d0"); // NonLinearity LHC12 Calo MassRatioFits
    cuts.AddCut("00010113","1111100067032220000","01631031000000d0"); // NonLinearity none
  } else if (trainConfig == 111){  // EMCAL clusters, different triggers no NonLinearity
    cuts.AddCut("00010113","1111100067032220000","01631031000000d0");
    cuts.AddCut("00052113","1111100067032220000","01631031000000d0"); // EMC7
    cuts.AddCut("00081113","1111100067032220000","01631031000000d0"); // EMCEG1,
  } else if (trainConfig == 112){ // EMCAL clusters, exotic cut var
    cuts.AddCut("00010113","1111111067032220000","01631031000000d0"); //
    cuts.AddCut("00010113","1111111067232220000","01631031000000d0"); //
    cuts.AddCut("00010113","1111111067332220000","01631031000000d0"); //
    cuts.AddCut("00010113","1111111067532220000","01631031000000d0"); //
    cuts.AddCut("00010113","1111111067732220000","01631031000000d0"); //
    cuts.AddCut("00010113","1111111067932220000","01631031000000d0"); //
  } else if (trainConfig == 113){  // trackMatching variations
    cuts.AddCut("00010113","1111111062032220000","01631031000000d0"); //
    cuts.AddCut("00010113","1111111063032220000","01631031000000d0"); //
    cuts.AddCut("00010113","1111111064032220000","01631031000000d0"); //
    cuts.AddCut("00010113","1111111065032220000","01631031000000d0"); //
    cuts.AddCut("00010113","1111111067032220000","01631031000000d0"); //
    
    // opening angle cut variations
  } else if (trainConfig == 114){ // EMCAL clusters pp 8 TeV
    cuts.AddCut("00010113","1111111067032220000","01631031000000a0"); // min open angle - 0. + 1 cell diag
  } else if (trainConfig == 115){ // EMCAL clusters pp 8 TeV
    cuts.AddCut("00010113","1111111067032220000","01631031000000b0"); // min open angle - 0.0152 + 1 cell diag
  } else if (trainConfig == 116){ // EMCAL clusters pp 8 TeV
    cuts.AddCut("00010113","1111111067032220000","01631031000000d0"); // min open angle - 0.017 + 1 cell diag
  } else if (trainConfig == 117){ // EMCAL clusters pp 8 TeV
    cuts.AddCut("00010113","1111111067032220000","01631031000000a0"); // min open angle - 0. + 1 cell diag
    cuts.AddCut("00010113","1111111067032220000","01631031000000d0"); // min open angle - 0.017 + 1 cell diag
    cuts.AddCut("00010113","1111111067032220000","01631031000000b0"); // min open angle - 0.0152 + 1 cell diag
    cuts.AddCut("00010113","1111111067032220000","01631031000000c0"); // min open angle - 0.016 + 1 cell diag
    cuts.AddCut("00010113","1111111067032220000","01631031000000e0"); // min open angle - 0.018 + 1 cell diag
    cuts.AddCut("00010113","1111111067032220000","01631031000000f0"); // min open angle - 0.019 + 1 cell diag

  } else if (trainConfig == 118){ // EMCAL clusters pp 8 TeV - no SPD PileUp
    cuts.AddCut("00010113","1111111067032220000","01631031000000d0"); // std
    cuts.AddCut("00010013","1111111067032220000","01631031000000d0"); // std - no pileup cut
    cuts.AddCut("00052113","1111111067032220000","01631031000000d0"); // std
    cuts.AddCut("00052013","1111111067032220000","01631031000000d0"); // std - no pileup cut
    cuts.AddCut("00081113","1111111067032220000","01631031000000d0"); // std
    cuts.AddCut("00081013","1111111067032220000","01631031000000d0"); // std - no pileup cut

  // only std cuts
  } else if (trainConfig == 119){ // EMCAL clusters pp 8 TeV
    cuts.AddCut("00010113","1111111067032220000","01631031000000d0"); // std
  } else if (trainConfig == 120){ // EMCAL clusters pp 8 TeV
    cuts.AddCut("00010113","1111111067032220000","01631031000000d0"); // std

    //8 TeV kEMC7 variations
  } else if (trainConfig == 121){ // EMCAL clusters pp 8 TeV, timing variation
    cuts.AddCut("00052113","1111111057032220000","01631031000000d0"); // time -50ns_50ns
    cuts.AddCut("00052113","1111111067032220000","01631031000000d0"); // time -30ns_35ns - standard
    cuts.AddCut("00052113","1111111077032220000","01631031000000d0"); // time -30ns_30ns
    cuts.AddCut("00052113","1111111087032220000","01631031000000d0"); // time -20ns_30ns
  } else if (trainConfig == 122){ //EMCAL minEnergy variation
    cuts.AddCut("00052113","1111111067012220000","01631031000000d0"); //0.5 GeV/c
    cuts.AddCut("00052113","1111111067022220000","01631031000000d0"); //0.6 GeV/c
    cuts.AddCut("00052113","1111111067032220000","01631031000000d0"); //0.7 GeV/c default
    cuts.AddCut("00052113","1111111067042220000","01631031000000d0"); //0.8 GeV/c
    cuts.AddCut("00052113","1111111067052220000","01631031000000d0"); //0.9 GeV/c
  } else if (trainConfig == 123){ //EMCAL minNCells, M02, with/without TRD variation
    cuts.AddCut("00052113","1111111067031220000","01631031000000d0"); //n cells >= 1
    cuts.AddCut("00052113","1111111067033220000","01631031000000d0"); //n cells >= 3
    cuts.AddCut("00052113","1111111067032200000","01631031000000d0"); //no max M02 cut
    cuts.AddCut("00052113","1111111067032250000","01631031000000d0"); //M02 < 0.27
    cuts.AddCut("00052113","1111111067032260000","01631031000000d0"); //M02 < 0.3
    cuts.AddCut("00052113","1112111067032220000","01631031000000d0"); //only modules with TRD infront
    cuts.AddCut("00052113","1111311067032220000","01631031000000d0"); //no modules with TRD infront
  } else if (trainConfig == 124){  // trackMatching variations
    cuts.AddCut("00052113","1111111066032220000","01631031000000d0"); //
    cuts.AddCut("00052113","1111111067032220000","01631031000000d0"); //
    cuts.AddCut("00052113","1111111068032220000","01631031000000d0"); //
    cuts.AddCut("00052113","1111111069032220000","01631031000000d0"); //
    cuts.AddCut("00052113","1111111063032220000","01631031000000d0"); //
    cuts.AddCut("00052113","1111111060032220000","01631031000000d0"); //
  } else if (trainConfig == 125){ // EMCAL clusters pp 8 TeV, combining cluster within time window and without
    cuts.AddCut("00052113","1111111007032220000","01631031000000d0"); //
  } else if (trainConfig == 126){ // EMCAL clusters open angle variation
    cuts.AddCut("00052113","1111111067032220000","0163103100000060"); // min open angle - 0.017 - std
    cuts.AddCut("00052113","1111111067032220000","0163103100000040"); // min open angle - 0.0152
    cuts.AddCut("00052113","1111111067032220000","0163103100000070"); // min open angle - 0.016
    cuts.AddCut("00052113","1111111067032220000","0163103100000080"); // min open angle - 0.018
    cuts.AddCut("00052113","1111111067032220000","0163103100000090"); // min open angle - 0.019
  } else if (trainConfig == 127){ // EMCAL clusters pp 8 TeV, Different DistanceToBadChannels
    cuts.AddCut("00052113","1111111067032220000","01631031000000d0"); //
    cuts.AddCut("00052113","1111111167032220000","01631031000000d0"); //
    cuts.AddCut("00052113","1111111267032220000","01631031000000d0"); //
    cuts.AddCut("00052113","1111111367032220000","01631031000000d0"); //
    cuts.AddCut("00052113","1111111567032220000","01631031000000d0"); //
    cuts.AddCut("00052113","1111111667032220000","01631031000000d0"); //
  } else if (trainConfig == 128){ // EMCAL clusters pp 8 TeV, Different NonLinearities
    cuts.AddCut("00052113","1111101067032220000","01631031000000d0"); // NonLinearity kSDMv5
    cuts.AddCut("00052113","1111113067032220000","01631031000000d0"); // NonLinearity kTestBeamv2 + LHC12 ConvCalo
    cuts.AddCut("00052113","1111114067032220000","01631031000000d0"); // NonLinearity kTestBeamv2 + LHC12 Calo
  } else if (trainConfig == 129){ // EMCAL clusters pp 8 TeV, Different NonLinearities
    cuts.AddCut("00052113","1111111067032220000","01631031000000d0"); // NonLinearity LHC12 ConvCalo
    cuts.AddCut("00052113","1111112067032220000","01631031000000d0"); // NonLinearity LHC12 Calo
    cuts.AddCut("00052113","1111121067032220000","01631031000000d0"); // NonLinearity LHC12 ConvCalo MassRatioFits
    cuts.AddCut("00052113","1111122067032220000","01631031000000d0"); // NonLinearity LHC12 Calo MassRatioFits
    cuts.AddCut("00052113","1111100067032220000","01631031000000d0"); // NonLinearity none
  } else if (trainConfig == 130){ // EMCAL clusters, exotic cut var
    cuts.AddCut("00052113","1111111067032220000","01631031000000d0"); //
    cuts.AddCut("00052113","1111111067232220000","01631031000000d0"); //
    cuts.AddCut("00052113","1111111067332220000","01631031000000d0"); //
    cuts.AddCut("00052113","1111111067532220000","01631031000000d0"); //
    cuts.AddCut("00052113","1111111067732220000","01631031000000d0"); //
    cuts.AddCut("00052113","1111111067932220000","01631031000000d0"); //
  } else if (trainConfig == 131){  // trackMatching variations
    cuts.AddCut("00052113","1111111067032220000","01631031000000d0"); //
    cuts.AddCut("00052113","1111111062032220000","01631031000000d0"); //
    cuts.AddCut("00052113","1111111063032220000","01631031000000d0"); //
    cuts.AddCut("00052113","1111111064032220000","01631031000000d0"); //

  } else if (trainConfig == 136){ // EMCAL clusters pp 8 TeV
    cuts.AddCut("00052113","1111111067032220000","01631031000000a0"); // min open angle - 0. + 1 cell diag
    cuts.AddCut("00052113","1111111067032220000","01631031000000d0"); // min open angle - 0.017 + 1 cell diag
    cuts.AddCut("00052113","1111111067032220000","01631031000000b0"); // min open angle - 0.0152 + 1 cell diag
    cuts.AddCut("00052113","1111111067032220000","01631031000000c0"); // min open angle - 0.016 + 1 cell diag
    cuts.AddCut("00052113","1111111067032220000","01631031000000e0"); // min open angle - 0.018 + 1 cell diag
    cuts.AddCut("00052113","1111111067032220000","01631031000000f0"); // min open angle - 0.019 + 1 cell diag

  // only std cuts
  } else if (trainConfig == 139){ // EMCAL clusters pp 8 TeV
    cuts.AddCut("00052113","1111111067032220000","01631031000000d0"); // std
  } else if (trainConfig == 140){ // EMCAL clusters pp 8 TeV
    cuts.AddCut("00052113","1111111067032220000","01631031000000d0"); // std

    //8 TeV kEMCEGA variations
  } else if (trainConfig == 141){ // EMCAL clusters pp 8 TeV, timing variation
    cuts.AddCut("00081113","1111111057032220000","01631031000000d0"); // time -50ns_50ns
    cuts.AddCut("00081113","1111111067032220000","01631031000000d0"); // time -30ns_35ns - standard
    cuts.AddCut("00081113","1111111077032220000","01631031000000d0"); // time -30ns_30ns
    cuts.AddCut("00081113","1111111087032220000","01631031000000d0"); // time -20ns_30ns
  } else if (trainConfig == 142){ //EMCAL minEnergy variation
    cuts.AddCut("00081113","1111111067012220000","01631031000000d0"); //0.5 GeV/c
    cuts.AddCut("00081113","1111111067022220000","01631031000000d0"); //0.6 GeV/c
    cuts.AddCut("00081113","1111111067032220000","01631031000000d0"); //0.7 GeV/c default
    cuts.AddCut("00081113","1111111067042220000","01631031000000d0"); //0.8 GeV/c
    cuts.AddCut("00081113","1111111067052220000","01631031000000d0"); //0.9 GeV/c
  } else if (trainConfig == 143){ //EMCAL minNCells, M02, with/without TRD variation
    cuts.AddCut("00081113","1111111067031220000","01631031000000d0"); //n cells >= 1
    cuts.AddCut("00081113","1111111067033220000","01631031000000d0"); //n cells >= 3
    cuts.AddCut("00081113","1111111067032200000","01631031000000d0"); //no max M02 cut
    cuts.AddCut("00081113","1111111067032250000","01631031000000d0"); //M02 < 0.3
    cuts.AddCut("00081113","1111111067032260000","01631031000000d0"); //M02 < 0.27
    cuts.AddCut("00081113","1112111067032220000","01631031000000d0"); //only modules with TRD infront
    cuts.AddCut("00081113","1111311067032220000","01631031000000d0"); //no modules with TRD infront
  } else if (trainConfig == 144){  // trackMatching variations
    cuts.AddCut("00081113","1111111066032220000","01631031000000d0"); //
    cuts.AddCut("00081113","1111111067032220000","01631031000000d0"); //
    cuts.AddCut("00081113","1111111068032220000","01631031000000d0"); //
    cuts.AddCut("00081113","1111111069032220000","01631031000000d0"); //
    cuts.AddCut("00081113","1111111063032220000","01631031000000d0"); //
    cuts.AddCut("00081113","1111111060032220000","01631031000000d0"); //
  } else if (trainConfig == 145){ // EMCAL clusters pp 8 TeV, combining cluster within time window and without
    cuts.AddCut("00081113","1111111007032220000","01631031000000d0"); //
  } else if (trainConfig == 146){ // EMCAL clusters open angle variation
    cuts.AddCut("00081113","1111111067032220000","0163103100000060"); // min open angle - 0.017 - std
    cuts.AddCut("00081113","1111111067032220000","0163103100000040"); // min open angle - 0.0152
    cuts.AddCut("00081113","1111111067032220000","0163103100000070"); // min open angle - 0.016
    cuts.AddCut("00081113","1111111067032220000","0163103100000080"); // min open angle - 0.018
    cuts.AddCut("00081113","1111111067032220000","0163103100000090"); // min open angle - 0.019
  } else if (trainConfig == 147){ // EMCAL clusters pp 8 TeV, Different DistanceToBadChannels
    cuts.AddCut("00081113","1111111067032220000","01631031000000d0"); //
    cuts.AddCut("00081113","1111111167032220000","01631031000000d0"); //
    cuts.AddCut("00081113","1111111267032220000","01631031000000d0"); //
    cuts.AddCut("00081113","1111111367032220000","01631031000000d0"); //
    cuts.AddCut("00081113","1111111567032220000","01631031000000d0"); //
    cuts.AddCut("00081113","1111111667032220000","01631031000000d0"); //
  } else if (trainConfig == 148){ // EMCAL clusters pp 8 TeV, Different NonLinearities
    cuts.AddCut("00081113","1111101067032220000","01631031000000d0"); // NonLinearity kSDMv5
    cuts.AddCut("00081113","1111113067032220000","01631031000000d0"); // NonLinearity kTestBeamv2 + LHC12 ConvCalo
    cuts.AddCut("00081113","1111114067032220000","01631031000000d0"); // NonLinearity kTestBeamv2 + LHC12 Calo
  } else if (trainConfig == 149){ // EMCAL clusters pp 8 TeV, Different NonLinearities
    cuts.AddCut("00081113","1111111067032220000","01631031000000d0"); // NonLinearity LHC12 ConvCalo
    cuts.AddCut("00081113","1111112067032220000","01631031000000d0"); // NonLinearity LHC12 Calo
    cuts.AddCut("00081113","1111121067032220000","01631031000000d0"); // NonLinearity LHC12 ConvCalo MassRatioFits
    cuts.AddCut("00081113","1111122067032220000","01631031000000d0"); // NonLinearity LHC12 Calo MassRatioFits
    cuts.AddCut("00081113","1111100067032220000","01631031000000d0"); // NonLinearity none
  } else if (trainConfig == 150){ // EMCAL clusters, exotic cut var
    cuts.AddCut("00081113","1111111067032220000","01631031000000d0"); //
    cuts.AddCut("00081113","1111111067232220000","01631031000000d0"); //
    cuts.AddCut("00081113","1111111067332220000","01631031000000d0"); //
    cuts.AddCut("00081113","1111111067532220000","01631031000000d0"); //
    cuts.AddCut("00081113","1111111067732220000","01631031000000d0"); //
    cuts.AddCut("00081113","1111111067932220000","01631031000000d0"); //
  } else if (trainConfig == 151){  // trackMatching variations
    cuts.AddCut("00081113","1111111067032220000","01631031000000d0"); //
    cuts.AddCut("00081113","1111111062032220000","01631031000000d0"); //
    cuts.AddCut("00081113","1111111063032220000","01631031000000d0"); //
    cuts.AddCut("00081113","1111111064032220000","01631031000000d0"); //

  } else if (trainConfig == 156){ // EMCAL clusters pp 8 TeV
    cuts.AddCut("00081113","1111111067032220000","01631031000000a0"); // min open angle - 0. + 1 cell diag
    cuts.AddCut("00081113","1111111067032220000","01631031000000d0"); // min open angle - 0.017 + 1 cell diag
    cuts.AddCut("00081113","1111111067032220000","01631031000000b0"); // min open angle - 0.0152 + 1 cell diag
    cuts.AddCut("00081113","1111111067032220000","01631031000000c0"); // min open angle - 0.016 + 1 cell diag
    cuts.AddCut("00081113","1111111067032220000","01631031000000e0"); // min open angle - 0.018 + 1 cell diag
    cuts.AddCut("00081113","1111111067032220000","01631031000000f0"); // min open angle - 0.019 + 1 cell diag

  // only std cuts
  } else if (trainConfig == 159){ // EMCAL clusters pp 8 TeV
    cuts.AddCut("00081113","1111111067032220000","01631031000000d0"); // std
  } else if (trainConfig == 160){ // EMCAL clusters pp 8 TeV
    cuts.AddCut("00081113","1111111067032220000","01631031000000d0"); // std

    //cut studies for opening angle, alpha dep. 0-0.2
  } else if (trainConfig == 161){ // EMCAL clusters open angle variation
    cuts.AddCut("00052113","1111111067032220000","016310a100000010");
    cuts.AddCut("00052113","1111111067032220000","016310a100000020");
    cuts.AddCut("00052113","1111111067032220000","016310a100000030");
    cuts.AddCut("00052113","1111111067032220000","016310a100000040");
    cuts.AddCut("00052113","1111111067032220000","016310a100000050");
  } else if (trainConfig == 162){ // EMCAL clusters open angle variation
    cuts.AddCut("00052113","1111111067032220000","016310a100000060");
    cuts.AddCut("00052113","1111111067032220000","016310a100000070");
    cuts.AddCut("00052113","1111111067032220000","016310a100000080");
    cuts.AddCut("00052113","1111111067032220000","016310a100000090");
    //cut studies for opening angle, alpha dep. 0.2-0.6
  } else if (trainConfig == 163){ // EMCAL clusters open angle variation
    cuts.AddCut("00052113","1111111067032220000","016310b100000010");
    cuts.AddCut("00052113","1111111067032220000","016310b100000020");
    cuts.AddCut("00052113","1111111067032220000","016310b100000030");
    cuts.AddCut("00052113","1111111067032220000","016310b100000040");
    cuts.AddCut("00052113","1111111067032220000","016310b100000050");
  } else if (trainConfig == 164){ // EMCAL clusters open angle variation
    cuts.AddCut("00052113","1111111067032220000","016310b100000060");
    cuts.AddCut("00052113","1111111067032220000","016310b100000070");
    cuts.AddCut("00052113","1111111067032220000","016310b100000080");
    cuts.AddCut("00052113","1111111067032220000","016310b100000090");
    //cut studies for opening angle, alpha dep. 0.6-1.0
  } else if (trainConfig == 165){ // EMCAL clusters open angle variation
    cuts.AddCut("00052113","1111111067032220000","016310c100000010");
    cuts.AddCut("00052113","1111111067032220000","016310c100000020");
    cuts.AddCut("00052113","1111111067032220000","016310c100000030");
    cuts.AddCut("00052113","1111111067032220000","016310c100000040");
    cuts.AddCut("00052113","1111111067032220000","016310c100000050");
  } else if (trainConfig == 166){ // EMCAL clusters open angle variation
    cuts.AddCut("00052113","1111111067032220000","016310c100000060");
    cuts.AddCut("00052113","1111111067032220000","016310c100000070");
    cuts.AddCut("00052113","1111111067032220000","016310c100000080");
    cuts.AddCut("00052113","1111111067032220000","016310c100000090");

    //cut studies for opening angle, alpha dep. 0-0.2
  } else if (trainConfig == 167){ // EMCAL clusters open angle variation
    cuts.AddCut("00081113","1111111067032220000","016310a100000010");
    cuts.AddCut("00081113","1111111067032220000","016310a100000020");
    cuts.AddCut("00081113","1111111067032220000","016310a100000030");
    cuts.AddCut("00081113","1111111067032220000","016310a100000040");
    cuts.AddCut("00081113","1111111067032220000","016310a100000050");
  } else if (trainConfig == 168){ // EMCAL clusters open angle variation
    cuts.AddCut("00081113","1111111067032220000","016310a100000060");
    cuts.AddCut("00081113","1111111067032220000","016310a100000070");
    cuts.AddCut("00081113","1111111067032220000","016310a100000080");
    cuts.AddCut("00081113","1111111067032220000","016310a100000090");
    //cut studies for opening angle, alpha dep. 0.2-0.6
  } else if (trainConfig == 169){ // EMCAL clusters open angle variation
    cuts.AddCut("00081113","1111111067032220000","016310b100000010");
    cuts.AddCut("00081113","1111111067032220000","016310b100000020");
    cuts.AddCut("00081113","1111111067032220000","016310b100000030");
    cuts.AddCut("00081113","1111111067032220000","016310b100000040");
    cuts.AddCut("00081113","1111111067032220000","016310b100000050");
  } else if (trainConfig == 170){ // EMCAL clusters open angle variation
    cuts.AddCut("00081113","1111111067032220000","016310b100000060");
    cuts.AddCut("00081113","1111111067032220000","016310b100000070");
    cuts.AddCut("00081113","1111111067032220000","016310b100000080");
    cuts.AddCut("00081113","1111111067032220000","016310b100000090");
    //cut studies for opening angle, alpha dep. 0.6-1.0
  } else if (trainConfig == 171){ // EMCAL clusters open angle variation
    cuts.AddCut("00081113","1111111067032220000","016310c100000010");
    cuts.AddCut("00081113","1111111067032220000","016310c100000020");
    cuts.AddCut("00081113","1111111067032220000","016310c100000030");
    cuts.AddCut("00081113","1111111067032220000","016310c100000040");
    cuts.AddCut("00081113","1111111067032220000","016310c100000050");
  } else if (trainConfig == 172){ // EMCAL clusters open angle variation
    cuts.AddCut("00081113","1111111067032220000","016310c100000060");
    cuts.AddCut("00081113","1111111067032220000","016310c100000070");
    cuts.AddCut("00081113","1111111067032220000","016310c100000080");
    cuts.AddCut("00081113","1111111067032220000","016310c100000090");

  // diff eta/rap cuts
  } else if (trainConfig == 178){ // EMCAL clusters pp 8 TeV, |eta| < 0.7, y < 0.7
    cuts.AddCut("00010113","1551111067032220000","01632031000000d0"); //
    cuts.AddCut("00052113","1551111067032220000","01632031000000d0"); //
    cuts.AddCut("00081113","1551111067032220000","01632031000000d0"); //
  } else if (trainConfig == 179){ // EMCAL clusters pp 8 TeV, |eta| < 0.3, y < 0.3
    cuts.AddCut("00010113","1661111067032220000","01637031000000d0"); //
    cuts.AddCut("00052113","1661111067032220000","01637031000000d0"); //
    cuts.AddCut("00081113","1661111067032220000","01637031000000d0"); //

  } else if (trainConfig == 180){ // EMCAL clusters pp 8 TeV, openangle 0.0152+1cell
    cuts.AddCut("00010113","1111111067032220000","01631031000000b0"); // std
    cuts.AddCut("00052113","1111111067032220000","01631031000000b0"); // std
    cuts.AddCut("00081113","1111111067032220000","01631031000000b0"); // std

  //multiple std cuts for different studies
  } else if (trainConfig == 181){ // EMCAL clusters pp 8 TeV
    cuts.AddCut("00010113","1111111067032220000","01631031000000d0"); // std
    cuts.AddCut("00052113","1111111067032220000","01631031000000d0"); // std
    cuts.AddCut("00081113","1111111067032220000","01631031000000d0"); // std
  } else if (trainConfig == 182){ // EMCAL clusters pp 8 TeV
    cuts.AddCut("00010113","1111111067032220000","01631031000000d0"); // std
    cuts.AddCut("00052113","1111111067032220000","01631031000000d0"); // std
    cuts.AddCut("00081113","1111111067032220000","01631031000000d0"); // std
  } else if (trainConfig == 183){ // EMCAL clusters pp 8 TeV
    cuts.AddCut("00010113","1111111067032220000","01631031000000d0"); // std
    cuts.AddCut("00052113","1111111067032220000","01631031000000d0"); // std
    cuts.AddCut("00081113","1111111067032220000","01631031000000d0"); // std
  } else if (trainConfig == 184){ // EMCAL clusters pp 8 TeV
    cuts.AddCut("00010113","1111111067032220000","01631031000000d0"); // std
    cuts.AddCut("00052113","1111111067032220000","01631031000000d0"); // std
    cuts.AddCut("00081113","1111111067032220000","01631031000000d0"); // std
  } else if (trainConfig == 185){ // EMCAL clusters pp 8 TeV
    cuts.AddCut("00010113","1111111067032220000","01631031000000d0"); // std
    cuts.AddCut("00052113","1111111067032220000","01631031000000d0"); // std
    cuts.AddCut("00081113","1111111067032220000","01631031000000d0"); // std
  } else if (trainConfig == 186){ // EMCAL clusters pp 8 TeV
    cuts.AddCut("00010113","1111111067032220000","01631031000000d0"); // std
    cuts.AddCut("00052113","1111111067032220000","01631031000000d0"); // std
    cuts.AddCut("00081113","1111111067032220000","01631031000000d0"); // std
  } else if (trainConfig == 187){ // EMCAL clusters pp 8 TeV
    cuts.AddCut("00010113","1111111067032220000","01631031000000d0"); // std
    cuts.AddCut("00052113","1111111067032220000","01631031000000d0"); // std
    cuts.AddCut("00081113","1111111067032220000","01631031000000d0"); // std
  } else if (trainConfig == 188){ // EMCAL clusters pp 8 TeV
    cuts.AddCut("00010113","1111111067032220000","01631031000000d0"); // std
    cuts.AddCut("00052113","1111111067032220000","01631031000000d0"); // std
    cuts.AddCut("00081113","1111111067032220000","01631031000000d0"); // std
  } else if (trainConfig == 189){ // EMCAL clusters pp 8 TeV
    cuts.AddCut("00010113","1111111067032220000","01631031000000d0"); // std
    cuts.AddCut("00052113","1111111067032220000","01631031000000d0"); // std
    cuts.AddCut("00081113","1111111067032220000","01631031000000d0"); // std
  } else if (trainConfig == 190){ // EMCAL clusters pp 8 TeV
    cuts.AddCut("00010113","1111111067032220000","01631031000000d0"); // std
    cuts.AddCut("00052113","1111111067032220000","01631031000000d0"); // std
    cuts.AddCut("00081113","1111111067032220000","01631031000000d0"); // std
  } else if (trainConfig == 191){ // EMCAL clusters pp 8 TeV
    cuts.AddCut("00010113","1111111067032220000","01631031000000d0"); // std
    cuts.AddCut("00052113","1111111067032220000","01631031000000d0"); // std
    cuts.AddCut("00081113","1111111067032220000","01631031000000d0"); // std
  } else if (trainConfig == 192){ // EMCAL clusters pp 8 TeV
    cuts.AddCut("00010113","1111111067032220000","01631031000000d0"); // std
    cuts.AddCut("00052113","1111111067032220000","01631031000000d0"); // std
    cuts.AddCut("00081113","1111111067032220000","01631031000000d0"); // std

  // pp multiplicity studies
  } else if (trainConfig == 198){ // MB - with multiplicity bins
    cuts.AddCut("00103113","1111121053032220000","0163103100000060"); // 0 -2
    cuts.AddCut("01203113","1111121053032220000","0163103100000060"); // 2 -5
    cuts.AddCut("02303113","1111121053032220000","0163103100000060"); // 5 -10
    cuts.AddCut("03403113","1111121053032220000","0163103100000060"); // 10 -30
    cuts.AddCut("04503113","1111121053032220000","0163103100000060"); // 30 -100
  } else if (trainConfig == 199){ // - with multiplicity bins
    cuts.AddCut("00100113","1111121063032220000","0163103100000060"); // 0 -2
    cuts.AddCut("01200113","1111121063032220000","0163103100000060"); // 2 -5
    cuts.AddCut("02300113","1111121063032220000","0163103100000060"); // 5 -10
    cuts.AddCut("03400113","1111121063032220000","0163103100000060"); // 10 -30
    cuts.AddCut("04500113","1111121063032220000","0163103100000060"); // 30 -100
    
    
    // 7 TeV
  } else if (trainConfig == 200){ // EMCAL clusters pp 7 TeV, pT dep matching
    cuts.AddCut("00000113","11111110b7032220000","01631031000000d0"); // std
    cuts.AddCut("00000113","1111111007032220000","01631031000000d0"); // std
  } else if (trainConfig == 201){ // EMCAL clusters pp 7 TeV, pT dep matching
    cuts.AddCut("00000113","11111110b7032220000","01631031000000d0"); // std
    cuts.AddCut("00000113","1111111007032220000","01631031000000d0"); // std
  } else if (trainConfig == 202){ // EMCAL clusters pp 7 TeV, timing variation
    cuts.AddCut("00000113","1111111037032220000","01631031000000d0"); //
    cuts.AddCut("00000113","1111111047032220000","01631031000000d0"); //
    cuts.AddCut("00000113","1111111057032220000","01631031000000d0"); //
    cuts.AddCut("00000113","1111111077032220000","01631031000000d0"); //
    cuts.AddCut("00000113","1111111087032220000","01631031000000d0"); //
  } else if (trainConfig == 203){ //EMCAL minEnergy variation
    cuts.AddCut("00000113","11111110b7012220000","01631031000000d0"); //0.5 GeV/c
    cuts.AddCut("00000113","11111110b7022220000","01631031000000d0"); //0.6 GeV/c
    cuts.AddCut("00000113","11111110b7032220000","01631031000000d0"); //0.7 GeV/c default
    cuts.AddCut("00000113","11111110b7042220000","01631031000000d0"); //0.8 GeV/c
    cuts.AddCut("00000113","11111110b7052220000","01631031000000d0"); //0.9 GeV/c
  } else if (trainConfig == 204){ //EMCAL minNCells, M02, with/without TRD variation
    cuts.AddCut("00000113","11111110b7031220000","01631031000000d0"); //n cells >= 1
    cuts.AddCut("00000113","11111110b7033220000","01631031000000d0"); //n cells >= 3
    cuts.AddCut("00000113","11111110b7032200000","01631031000000d0"); //no max M02 cut
    cuts.AddCut("00000113","11111110b7032250000","01631031000000d0"); //M02 < 0.3
    cuts.AddCut("00000113","11111110b7032260000","01631031000000d0"); //M02 < 0.27
    cuts.AddCut("00000113","11131110b7032220000","01631031000000d0"); //only modules with TRD infront
    cuts.AddCut("00000113","11112110b7032220000","01631031000000d0"); //no modules with TRD infront
  } else if (trainConfig == 205){  // trackMatching variations
    cuts.AddCut("00000113","11111110b6032220000","01631031000000d0"); //
    cuts.AddCut("00000113","11111110b7032220000","01631031000000d0"); //
    cuts.AddCut("00000113","11111110b8032220000","01631031000000d0"); //
    cuts.AddCut("00000113","11111110b9032220000","01631031000000d0"); //
    cuts.AddCut("00000113","11111110b3032220000","01631031000000d0"); //
    cuts.AddCut("00000113","11111110b0032220000","01631031000000d0"); //
  } else if (trainConfig == 206){ // EMCAL clusters pp 7 TeV, combining cluster within time window and without
    cuts.AddCut("00000113","1111111007032220000","01631031000000d0"); //
  } else if (trainConfig == 207){ // EMCAL clusters open angle variation
    cuts.AddCut("00000113","11111110b7032220000","01631031000000a0"); // min open angle - 0. + 1 cell diag
    cuts.AddCut("00000113","11111110b7032220000","01631031000000d0"); // min open angle - 0.017 + 1 cell diag
    cuts.AddCut("00000113","11111110b7032220000","01631031000000b0"); // min open angle - 0.0152 + 1 cell diag
    cuts.AddCut("00000113","11111110b7032220000","01631031000000c0"); // min open angle - 0.016 + 1 cell diag
    cuts.AddCut("00000113","11111110b7032220000","01631031000000e0"); // min open angle - 0.018 + 1 cell diag
    cuts.AddCut("00000113","11111110b7032220000","01631031000000f0"); // min open angle - 0.019 + 1 cell diag
  } else if (trainConfig == 208){ // EMCAL clusters pp 7 TeV, Different DistanceToBadChannels
    cuts.AddCut("00000113","11111110b7032220000","01631031000000d0"); //
    cuts.AddCut("00000113","11111111b7032220000","01631031000000d0"); //
    cuts.AddCut("00000113","11111112b7032220000","01631031000000d0"); //
    cuts.AddCut("00000113","11111113b7032220000","01631031000000d0"); //
    cuts.AddCut("00000113","11111115b7032220000","01631031000000d0"); //
    cuts.AddCut("00000113","11111116b7032220000","01631031000000d0"); //
  } else if (trainConfig == 209){ // EMCAL clusters pp 7 TeV, Different NonLinearities
    cuts.AddCut("00000113","11111010b7032220000","01631031000000d0"); // NonLinearity kSDMv5
    cuts.AddCut("00000113","11111130b7032220000","01631031000000d0"); // NonLinearity kTestBeamv2 + LHC12 ConvCalo
    cuts.AddCut("00000113","11111140b7032220000","01631031000000d0"); // NonLinearity kTestBeamv2 + LHC12 Calo
  } else if (trainConfig == 210){ // EMCAL clusters pp 7 TeV, Different NonLinearities
    cuts.AddCut("00000113","11111110b7032220000","01631031000000d0"); // NonLinearity LHC12 ConvCalo
    cuts.AddCut("00000113","11111120b7032220000","01631031000000d0"); // NonLinearity LHC12 Calo
    cuts.AddCut("00000113","11111210b7032220000","01631031000000d0"); // NonLinearity LHC12 ConvCalo MassRatioFits
    cuts.AddCut("00000113","11111220b7032220000","01631031000000d0"); // NonLinearity LHC12 Calo MassRatioFits
    cuts.AddCut("00000113","11111000b7032220000","01631031000000d0"); // NonLinearity none
  } else if (trainConfig == 211){  // EMCAL clusters, different triggers no NonLinearity
    cuts.AddCut("00000113","11111000b7032220000","01631031000000d0");
  } else if (trainConfig == 212){  // trackMatching variations
    cuts.AddCut("00000113","11111110b7032220000","01631031000000d0"); //
    cuts.AddCut("00000113","11111110b2032220000","01631031000000d0"); //
    cuts.AddCut("00000113","11111110b3032220000","01631031000000d0"); //
    cuts.AddCut("00000113","11111110b4032220000","01631031000000d0"); //
    cuts.AddCut("00000113","11111110b5032220000","01631031000000d0"); //

  } else if (trainConfig == 214){ // EMCAL clusters pp 8 TeV
    cuts.AddCut("00000113","11111110b7032220000","01631031000000a0"); // min open angle - 0. + 1 cell diag
  } else if (trainConfig == 215){ // EMCAL clusters pp 8 TeV
    cuts.AddCut("00000113","11111110b7032220000","01631031000000b0"); // min open angle - 0.0152 + 1 cell diag
  } else if (trainConfig == 216){ // EMCAL clusters pp 8 TeV
    cuts.AddCut("00000113","11111110b7032220000","01631031000000d0"); // min open angle - 0.017 + 1 cell diag
  } else if (trainConfig == 217){ // EMCAL clusters pp 8 TeV
    cuts.AddCut("00000113","11111110b7032220000","01631031000000a0"); // min open angle - 0. + 1 cell diag
    cuts.AddCut("00000113","11111110b7032220000","01631031000000d0"); // min open angle - 0.017 + 1 cell diag
    cuts.AddCut("00000113","11111110b7032220000","01631031000000b0"); // min open angle - 0.0152 + 1 cell diag
    cuts.AddCut("00000113","11111110b7032220000","01631031000000c0"); // min open angle - 0.016 + 1 cell diag
    cuts.AddCut("00000113","11111110b7032220000","01631031000000e0"); // min open angle - 0.018 + 1 cell diag
    cuts.AddCut("00000113","11111110b7032220000","01631031000000f0"); // min open angle - 0.019 + 1 cell diag

  } else if (trainConfig == 221){ // EMCAL clusters pp 7 TeV, std matching
    cuts.AddCut("00000113","11111110b3032220000","01631031000000d0"); // std
  } else if (trainConfig == 222){ // EMCAL clusters pp 7 TeV, no SPD pileup
    cuts.AddCut("00000113","11111110b7032220000","01631031000000d0"); // std
    cuts.AddCut("00000013","11111110b7032220000","01631031000000d0"); // std - no SPD pileup

    // pp7TeV EMCal direct photons
  } else if (trainConfig == 251) {
    cuts.AddCut("00000113","1111111063032220000","0163103100000060"); // pt const track matching, M02 < 0.7
    cuts.AddCut("00000113","1111111063032230000","0163103100000060"); // pt const track matching, M02 < 0.5
    cuts.AddCut("00000113","1111111067032220000","0163103100000060"); // pt dep track matching, M02 < 0.7
    cuts.AddCut("00000113","1111111067032230000","0163103100000060"); // pt dep track matching, M02 < 0.5
    
  } else if (trainConfig == 299){ // EMCAL clusters pp, jet triggers
    cuts.AddCut("00045113","1111111063032220000","0163103100000060"); // std
    cuts.AddCut("00046113","1111111063032220000","0163103100000060"); // std
    cuts.AddCut("00091113","1111111063032220000","0163103100000060"); // std
    cuts.AddCut("00092113","1111111063032220000","0163103100000060"); // std

  // ************************************* PHOS cuts ****************************************************
  } else if (trainConfig == 301) { //PHOS clusters
    cuts.AddCut("00003113","2444400040033200000","0163803100000010"); //pp LHC11a with SDD, PHOS
    cuts.AddCut("00010113","2444400040033200000","0163803100000010"); //pp LHC13g default MB
    cuts.AddCut("00061113","2444400040033200000","0163803100000010"); //pp LHC11a PHI1
    cuts.AddCut("00062113","2444400040033200000","0163803100000010"); //pp LHC11a PHI7
  } else if (trainConfig == 302){ // Validation PHOS
    cuts.AddCut("00003113","2444400040033200000","0163803100000010");
  } else if (trainConfig == 303){ // PHOS clusters, without and with added signals
    cuts.AddCut("00003113","2444400040033200000","0163803100000010");
    cuts.AddCut("00003123","2444400040033200000","0163803100000010");
    
  // 7 TeV direct photon PHOS
  } else if (trainConfig == 351){
    cuts.AddCut("00000113","2444400000013300000","0163803100000010"); // no nonlinearity
    cuts.AddCut("00000113","2444401000013300000","0163803100000010"); // with PHOS nonlinearity
  } else if (trainConfig == 352){
    cuts.AddCut("00000113","2444400040013300000","0163803100000010"); // 100ns timing cut, no track matching
    cuts.AddCut("00000113","2444400043013300000","0163803100000010"); // 100ns timing cut
    cuts.AddCut("00000113","2444400043013350000","0163803100000010"); // 100ns timing cut, M02<0.3
    cuts.AddCut("00000113","2444400043013330000","0163803100000010"); // 100ns timing cut, M02<0.5
    cuts.AddCut("00000113","2444400043013320000","0163803100000010"); // 100ns timing cut, M02<0.7
  } else if (trainConfig == 353){ // same as 352 but with PHOS nonlinearity
    cuts.AddCut("00000113","2444401040013300000","0163803100000010"); // 100ns timing cut, no track matching
    cuts.AddCut("00000113","2444401043013300000","0163803100000010"); // 100ns timing cut
    cuts.AddCut("00000113","2444401043013350000","0163803100000010"); // 100ns timing cut, M02<0.3
    cuts.AddCut("00000113","2444401043013330000","0163803100000010"); // 100ns timing cut, M02<0.5
    cuts.AddCut("00000113","2444401043013320000","0163803100000010"); // 100ns timing cut, M02<0.7
  
  // 7 TeV PHOS
  } else if (trainConfig == 361){
    cuts.AddCut("00000113","2444400000013300000","0163803100000010"); // QA
    cuts.AddCut("00000113","2444400040013300000","0163803100000010"); // 100ns timing cut, no track matching
    cuts.AddCut("00000113","2444400043013300000","0163803100000010"); // 100ns timing cut
  } else if (trainConfig == 362){
    cuts.AddCut("00062113","2444400000013300000","0163803100000010"); // QA
    cuts.AddCut("00062113","2444400040013300000","0163803100000010"); // 100ns timing cut, no track matching
    cuts.AddCut("00062113","2444400043013300000","0163803100000010"); // 100ns timing cut
    
  // PHOS @ 8 TeV
  } else if (trainConfig == 381){ // PHOS clusters
    cuts.AddCut("00010113","2444400040013300000","0163803100000010");
  } else if (trainConfig == 382){ // PHOS clusters
    cuts.AddCut("00062113","2444400040013300000","0163803100000010");
  // PHOS clusters RUN2 config pp 5 TeV
  } else if (trainConfig == 383){ // PHOS clusters with larger acceptance
    cuts.AddCut("00010113","2446600040013300000","0163803100000010");
    cuts.AddCut("00062113","2446600040013300000","0163803100000010");

    // 13 TeV & 5 TeV
  } else if (trainConfig == 401){ // EMCAL clusters
    cuts.AddCut("00010113","1111100017032220000","0163103100000060"); // 1000ns timing cut, no NL INT7
    cuts.AddCut("00052013","1111100017032220000","0163103100000060"); // 1000ns timing cut, no NL EMC7
    cuts.AddCut("00085013","1111100017032220000","0163103100000060"); // 1000ns timing cut, no NL EG2
    cuts.AddCut("00083013","1111100017032220000","0163103100000060"); // 1000ns timing cut, no NL EG1
  } else if (trainConfig == 402){ // EMCAL clusters
    cuts.AddCut("00010113","1111100067032220000","0163103100000060"); // -50ns, 30ns timing cut, no NL INT7
    cuts.AddCut("00052013","1111100067032220000","0163103100000060"); // -50ns, 30ns timing cut, no NL EMC7
    cuts.AddCut("00085013","1111100067032220000","0163103100000060"); // -50ns, 30ns timing cut, no NL EG2
    cuts.AddCut("00083013","1111100067032220000","0163103100000060"); // -50ns, 30ns timing cut, no NL EG1
  } else if (trainConfig == 403){ // EMCAL clusters - NonLin INT7
    cuts.AddCut("00010113","1111111017032220000","0163103100000060"); // 1000ns timing cut, NL kSDM PCMEMC
    cuts.AddCut("00010113","1111112017032220000","0163103100000060"); // 1000ns timing cut, NL kSDM EMC
    cuts.AddCut("00010113","1111121017032220000","0163103100000060"); // 1000ns timing cut, NL DExt PCMEMC
    cuts.AddCut("00010113","1111122017032220000","0163103100000060"); // 1000ns timing cut, NL DExt EMC
  } else if (trainConfig == 404){ // EMCAL clusters - NonLin INT7
    cuts.AddCut("00010113","1111111067032220000","0163103100000060"); // -50ns, 30ns timing cut, NL kSDM PCMEMC
    cuts.AddCut("00010113","1111112067032220000","0163103100000060"); // -50ns, 30ns timing cut, NL kSDM EMC
    cuts.AddCut("00010113","1111121067032220000","0163103100000060"); // -50ns, 30ns timing cut, NL DExt PCMEMC
    cuts.AddCut("00010113","1111122067032220000","0163103100000060"); // -50ns, 30ns timing cut, NL DExt EMC
  } else if (trainConfig == 405){ // EMCAL clusters - NonLin EMC7
    cuts.AddCut("00052013","1111111017032220000","0163103100000060"); // 1000ns timing cut, NL kSDM PCMEMC
    cuts.AddCut("00052013","1111112017032220000","0163103100000060"); // 1000ns timing cut, NL kSDM EMC
    cuts.AddCut("00052013","1111121017032220000","0163103100000060"); // 1000ns timing cut, NL DExt PCMEMC
    cuts.AddCut("00052013","1111122017032220000","0163103100000060"); // 1000ns timing cut, NL DExt EMC
  } else if (trainConfig == 406){ // EMCAL clusters - NonLin EMC7
    cuts.AddCut("00052013","1111111067032220000","0163103100000060"); // -50ns, 30ns timing cut, NL kSDM PCMEMC
    cuts.AddCut("00052013","1111112067032220000","0163103100000060"); // -50ns, 30ns timing cut, NL kSDM EMC
    cuts.AddCut("00052013","1111121067032220000","0163103100000060"); // -50ns, 30ns timing cut, NL DExt PCMEMC
    cuts.AddCut("00052013","1111122067032220000","0163103100000060"); // -50ns, 30ns timing cut, NL DExt EMC
  } else if (trainConfig == 407){ // EMCAL clusters - NonLin EG2
    cuts.AddCut("00085013","1111111017032220000","0163103100000060"); // 1000ns timing cut, NL kSDM PCMEMC
    cuts.AddCut("00085013","1111112017032220000","0163103100000060"); // 1000ns timing cut, NL kSDM EMC
    cuts.AddCut("00085013","1111121017032220000","0163103100000060"); // 1000ns timing cut, NL DExt PCMEMC
    cuts.AddCut("00085013","1111122017032220000","0163103100000060"); // 1000ns timing cut, NL DExt EMC
  } else if (trainConfig == 408){ // EMCAL clusters - NonLin EG2
    cuts.AddCut("00085013","1111111067032220000","0163103100000060"); // -50ns, 30ns timing cut, NL kSDM PCMEMC
    cuts.AddCut("00085013","1111112067032220000","0163103100000060"); // -50ns, 30ns timing cut, NL kSDM EMC
    cuts.AddCut("00085013","1111121067032220000","0163103100000060"); // -50ns, 30ns timing cut, NL DExt PCMEMC
    cuts.AddCut("00085013","1111122067032220000","0163103100000060"); // -50ns, 30ns timing cut, NL DExt EMC
  } else if (trainConfig == 409){ // EMCAL clusters - NonLin EG1
    cuts.AddCut("00083013","1111111017032220000","0163103100000060"); // 1000ns timing cut, NL kSDM PCMEMC
    cuts.AddCut("00083013","1111112017032220000","0163103100000060"); // 1000ns timing cut, NL kSDM EMC
    cuts.AddCut("00083013","1111121017032220000","0163103100000060"); // 1000ns timing cut, NL DExt PCMEMC
    cuts.AddCut("00083013","1111122017032220000","0163103100000060"); // 1000ns timing cut, NL DExt EMC
  } else if (trainConfig == 410){ // EMCAL clusters - NonLin EG1
    cuts.AddCut("00083013","1111111067032220000","0163103100000060"); // -50ns, 30ns timing cut, NL kSDM PCMEMC
    cuts.AddCut("00083013","1111112067032220000","0163103100000060"); // -50ns, 30ns timing cut, NL kSDM EMC
    cuts.AddCut("00083013","1111121067032220000","0163103100000060"); // -50ns, 30ns timing cut, NL DExt PCMEMC
    cuts.AddCut("00083013","1111122067032220000","0163103100000060"); // -50ns, 30ns timing cut, NL DExt EMC

  // 2.76TeV additional configurations for y range change
  } else if (trainConfig == 501){ // EMCAL clusters 2.76 TeV LHC11a, with SDD (0), kEMC1 (1)
    cuts.AddCut("00003113","1551121053032220000","0163203100000050"); // |eta| < 0.7, y < 0.7
    cuts.AddCut("00051013","1551121053032220000","0163203100000050"); // |eta| < 0.7, y < 0.7
    cuts.AddCut("00003113","1661121053032220000","0163703100000050"); // |eta| < 0.3
    cuts.AddCut("00051013","1661121053032220000","0163703100000050"); // |eta| < 0.3
  } else if (trainConfig == 502){ // EMCAL clusters 2.76 TeV LHC11a, with SDD (0), kEMC1 (1)
    cuts.AddCut("00003113","1551121057032220000","0163203100000050"); // |eta| < 0.7, y < 0.7
    cuts.AddCut("00051013","1551121057032220000","0163203100000050"); // |eta| < 0.7, y < 0.7
    cuts.AddCut("00003113","1661121057032220000","0163703100000050"); // |eta| < 0.3
    cuts.AddCut("00051013","1661121057032220000","0163703100000050"); // |eta| < 0.3
  } else if (trainConfig == 503){  // eta < 0.7
    cuts.AddCut("00010113","1551121063032220000","0163203100000050"); 
    cuts.AddCut("00052013","1551121063032220000","0163203100000050"); // EMC7
    cuts.AddCut("00083013","1551121063032220000","0163203100000050"); // EMCEG1,
    cuts.AddCut("00085013","1551121063032220000","0163203100000050"); // EMCEG2,
  } else if (trainConfig == 504){  // eta < 0.3
    cuts.AddCut("00010113","1661121063032220000","0163703100000050"); 
    cuts.AddCut("00052013","1661121063032220000","0163703100000050"); // EMC7
    cuts.AddCut("00083013","1661121063032220000","0163703100000050"); // EMCEG1,
    cuts.AddCut("00085013","1661121063032220000","0163703100000050"); // EMCEG2,
  } else if (trainConfig == 505){  // eta < 0.7, pt dependent TM
    cuts.AddCut("00010113","1551121067032220000","0163203100000050"); 
    cuts.AddCut("00052013","1551121067032220000","0163203100000050"); // EMC7
    cuts.AddCut("00083013","1551121067032220000","0163203100000050"); // EMCEG1,
    cuts.AddCut("00085013","1551121067032220000","0163203100000050"); // EMCEG2,
  } else if (trainConfig == 506){  // eta < 0.3, pt dependent TM
    cuts.AddCut("00010113","1661121067032220000","0163703100000050"); 
    cuts.AddCut("00052013","1661121067032220000","0163703100000050"); // EMC7
    cuts.AddCut("00083013","1661121067032220000","0163703100000050"); // EMCEG1,
    cuts.AddCut("00085013","1661121067032220000","0163703100000050"); // EMCEG2,
    
  } else if (trainConfig == 507){ // pt dependent TM
    cuts.AddCut("00003113","1111121057032220000","0163103100000050"); // 700 MeV cluster min energy
    cuts.AddCut("00051013","1111121057032220000","0163103100000050"); // 700 MeV cluster min energy
    cuts.AddCut("00003013","1111121057032220000","0163103100000050"); // 700 MeV cluster min energy
  } else if (trainConfig == 508){ // pt dependent TM
    cuts.AddCut("00010113","1111121067032220000","0163103100000050"); 
    cuts.AddCut("00052013","1111121067032220000","0163103100000050"); // EMC7
    cuts.AddCut("00083013","1111121067032220000","0163103100000050"); // EMCEG1,
    cuts.AddCut("00085013","1111121067032220000","0163103100000050"); // EMCEG2,
   // ************************************* DCal cuts ****************************************************
  } else if (trainConfig == 601){ // DCAL clusters pp 5.02 TeV
    cuts.AddCut("00010113","3115511081041220000","0163103100000050"); // std
    cuts.AddCut("00055113","3115511081041220000","0163103100000050"); // std
    cuts.AddCut("00089113","3115511081041220000","0163103100000050"); // std
    cuts.AddCut("0008b113","3115511081041220000","0163103100000050"); // std
  } else if (trainConfig == 602){ // timing Cut variation  std -20+50ns
    cuts.AddCut("00010113","3115511071041220000","0163103100000050"); //     +-25   ns
    cuts.AddCut("00010113","3115511091041220000","0163103100000050"); //     +-12.5 ns
  } else if (trainConfig == 603){ // track matching variation
    cuts.AddCut("00010113","3115511082041220000","0163103100000050"); // 
    cuts.AddCut("00010113","3115511083041220000","0163103100000050"); // 
    cuts.AddCut("00010113","3115511084041220000","0163103100000050"); //  
    cuts.AddCut("00010113","3115511085041220000","0163103100000050"); //  
  } else if (trainConfig == 604){ // opening angle variation
    cuts.AddCut("00010113","3115511081041220000","0163103100000000"); // 
    cuts.AddCut("00010113","3115511081041220000","0163103100000010"); // 
    cuts.AddCut("00010113","3115511081041220000","0163103100000030"); // 
    cuts.AddCut("00010113","3115511081041220000","0163103100000040"); // 
    cuts.AddCut("00010113","3115511081041220000","0163103100000080"); // 
  } else if (trainConfig == 605){ // min energy variation std 0.7 GeV/c
    cuts.AddCut("00010113","3115511081001220000","0163103100000050"); //     0.01GeV/c
    cuts.AddCut("00010113","3115511081031220000","0163103100000050"); //     0.6 GeV/c
    cuts.AddCut("00010113","3115511081061220000","0163103100000050"); //     0.8 GeV/c
    cuts.AddCut("00010113","3115511081021220000","0163103100000050"); //     0.5 GeV/c
    cuts.AddCut("00010113","3115511081081220000","0163103100000050"); //     0.9 GeV/c
  } else if (trainConfig == 606){ // min nCells & M02 variation 
    // std: min nCells = 1; M02 max=0.7, min=0.1
    cuts.AddCut("00010113","3115511081042220000","0163103100000050"); //   min nCells = 2
    cuts.AddCut("00010113","3115511081043220000","0163103100000050"); //   min nCells = 3
    cuts.AddCut("00010113","3115511081041210000","0163103100000050"); //   max M02 = 1
    cuts.AddCut("00010113","3115511081041240000","0163103100000050"); //   max M02 = 0.4
  
  // ********************************* Past future cutstudies ******************************************
  } else if (trainConfig == 700){ // EMCAL clusters pp 8 TeV MinBias
    cuts.AddCut("00010113","1111111067032220000","0163103100000060"); // std
    cuts.AddCut("00010313","1111111067032220000","0163103100000060"); // std pastfuture -100ns/175ns
    cuts.AddCut("00010413","1111111067032220000","0163103100000060"); // std pastfuture -250ns/325ns
    cuts.AddCut("00010513","1111111067032220000","0163103100000060"); // std pastfuture -1000ns/1075ns
  } else if (trainConfig == 701){ // EMCAL clusters pp 8 TeV EMC7
    cuts.AddCut("00052113","1111111067032220000","0163103100000060"); // std
    cuts.AddCut("00052313","1111111067032220000","0163103100000060"); // std pastfuture -100ns/175ns
    cuts.AddCut("00052413","1111111067032220000","0163103100000060"); // std pastfuture -250ns/325ns
    cuts.AddCut("00052513","1111111067032220000","0163103100000060"); // std pastfuture -1000ns/1075ns
  } else if (trainConfig == 702){ // EMCAL clusters pp 8 TeV EMCEGA
    cuts.AddCut("00081113","1111111067032220000","0163103100000060"); // std
    cuts.AddCut("00081313","1111111067032220000","0163103100000060"); // std pastfuture -100ns/175ns
    cuts.AddCut("00081413","1111111067032220000","0163103100000060"); // std pastfuture -250ns/325ns
    cuts.AddCut("00081513","1111111067032220000","0163103100000060"); // std pastfuture -1000ns/1075ns
  } else if (trainConfig == 703){ // EMCAL clusters pp 8 TeV MinBias, cluster time 100ns
    cuts.AddCut("00010113","1111111047032220000","0163103100000060"); // std
    cuts.AddCut("00010313","1111111047032220000","0163103100000060"); // std pastfuture -100ns/175ns
    cuts.AddCut("00010413","1111111047032220000","0163103100000060"); // std pastfuture -250ns/325ns
    cuts.AddCut("00010513","1111111047032220000","0163103100000060"); // std pastfuture -1000ns/1075ns
  } else if (trainConfig == 704){ // EMCAL clusters pp 8 TeV EMC7, cluster time 100ns
    cuts.AddCut("00052113","1111111047032220000","0163103100000060"); // std
    cuts.AddCut("00052313","1111111047032220000","0163103100000060"); // std pastfuture -100ns/175ns
    cuts.AddCut("00052413","1111111047032220000","0163103100000060"); // std pastfuture -250ns/325ns
    cuts.AddCut("00052513","1111111047032220000","0163103100000060"); // std pastfuture -1000ns/1075ns
  } else if (trainConfig == 705){ // EMCAL clusters pp 8 TeV EMCEGA, cluster time 100ns
    cuts.AddCut("00081113","1111111047032220000","0163103100000060"); // std
    cuts.AddCut("00081313","1111111047032220000","0163103100000060"); // std pastfuture -100ns/175ns
    cuts.AddCut("00081413","1111111047032220000","0163103100000060"); // std pastfuture -250ns/325ns
    cuts.AddCut("00081513","1111111047032220000","0163103100000060"); // std pastfuture -1000ns/1075ns
  } else if (trainConfig == 706){ // EMCAL clusters pp 8 TeV MinBias, cluster time 1000ns
    cuts.AddCut("00010113","1111111017032220000","0163103100000060"); // std
    cuts.AddCut("00010313","1111111017032220000","0163103100000060"); // std pastfuture -100ns/175ns
    cuts.AddCut("00010413","1111111017032220000","0163103100000060"); // std pastfuture -250ns/325ns
    cuts.AddCut("00010513","1111111017032220000","0163103100000060"); // std pastfuture -1000ns/1075ns
  } else if (trainConfig == 707){ // EMCAL clusters pp 8 TeV EMC7, cluster time 1000ns
    cuts.AddCut("00052113","1111111017032220000","0163103100000060"); // std
    cuts.AddCut("00052313","1111111017032220000","0163103100000060"); // std pastfuture -100ns/175ns
    cuts.AddCut("00052413","1111111017032220000","0163103100000060"); // std pastfuture -250ns/325ns
    cuts.AddCut("00052513","1111111017032220000","0163103100000060"); // std pastfuture -1000ns/1075ns
  } else if (trainConfig == 708){ // EMCAL clusters pp 8 TeV EMCEGA, cluster time 1000ns
    cuts.AddCut("00081113","1111111017032220000","0163103100000060"); // std
    cuts.AddCut("00081313","1111111017032220000","0163103100000060"); // std pastfuture -100ns/175ns
    cuts.AddCut("00081413","1111111017032220000","0163103100000060"); // std pastfuture -250ns/325ns
    cuts.AddCut("00081513","1111111017032220000","0163103100000060"); // std pastfuture -1000ns/1075ns

  } else if (trainConfig == 710){ // PHOS clusters pp 8 TeV MinBias
    cuts.AddCut("00010113","2444400043013300000","0163803100000010");
    cuts.AddCut("00010313","2444400043013300000","0163803100000010"); // pastfuture -100ns/175ns
    cuts.AddCut("00010413","2444400043013300000","0163803100000010"); // pastfuture -250ns/325ns
    cuts.AddCut("00010513","2444400043013300000","0163803100000010"); // pastfuture -1000ns/1075ns
  } else if (trainConfig == 711){ // PHOS clusters pp 8 TeV PHI7
    cuts.AddCut("00062113","2444400043013300000","0163803100000010");
    cuts.AddCut("00062313","2444400043013300000","0163803100000010"); // pastfuture -100ns/175ns
    cuts.AddCut("00062413","2444400043013300000","0163803100000010"); // pastfuture -250ns/325ns
    cuts.AddCut("00062513","2444400043013300000","0163803100000010"); // pastfuture -1000ns/1075ns
  } else if (trainConfig == 712){ // PHOS clusters pp 8 TeV MinBias, cluster time 1000ns
    cuts.AddCut("00010113","2444400013013300000","0163803100000010");
    cuts.AddCut("00010313","2444400013013300000","0163803100000010"); // pastfuture -100ns/175ns
    cuts.AddCut("00010413","2444400013013300000","0163803100000010"); // pastfuture -250ns/325ns
    cuts.AddCut("00010513","2444400014013300000","0163803100000010"); // pastfuture -1000ns/1075ns
  } else if (trainConfig == 713){ // PHOS clusters pp 8 TeV PHI7, cluster time 1000ns
    cuts.AddCut("00062113","2444400013013300000","0163803100000010");
    cuts.AddCut("00062313","2444400013013300000","0163803100000010"); // pastfuture -100ns/175ns
    cuts.AddCut("00062413","2444400013013300000","0163803100000010"); // pastfuture -250ns/325ns
    cuts.AddCut("00062513","2444400013013300000","0163803100000010"); // pastfuture -1000ns/1075ns
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
    if (periodNameAnchor.CompareTo("LHC13g") == 0 && triggerString.CompareTo("10")== 0 )
      triggerString         = "00";
    
    
    dataInputMultHisto      = Form("%s_%s", periodNameAnchor.Data(), triggerString.Data());
    mcInputMultHisto        = Form("%s_%s", periodNameV0Reader.Data(), triggerString.Data());
   
    if (doMultiplicityWeighting) analysisEventCuts[i]->SetUseWeightMultiplicityFromFile( kTRUE, fileNameInputForMultWeighing, dataInputMultHisto, mcInputMultHisto );

    analysisEventCuts[i]->SetTriggerMimicking(enableTriggerMimicking);
    analysisEventCuts[i]->SetTriggerOverlapRejecion(enableTriggerOverlapRej);
    analysisEventCuts[i]->SetMaxFacPtHard(maxFacPtHard);
    analysisEventCuts[i]->SetV0ReaderName(V0ReaderName);
    analysisEventCuts[i]->SetLightOutput(runLightOutput);
    if (periodNameV0Reader.CompareTo("") != 0) analysisEventCuts[i]->SetPeriodEnum(periodNameV0Reader);
    analysisEventCuts[i]->InitializeCutsFromCutString((cuts.GetEventCut(i)).Data());
    EventCutList->Add(analysisEventCuts[i]);
    analysisEventCuts[i]->SetFillCutHistograms("",kFALSE);
    
    analysisClusterCuts[i] = new AliCaloPhotonCuts(isMC);
    analysisClusterCuts[i]->SetHistoToModifyAcceptance(histoAcc);
    analysisClusterCuts[i]->SetV0ReaderName(V0ReaderName);
    analysisClusterCuts[i]->SetCaloTrackMatcherName(TrackMatcherName);
    analysisClusterCuts[i]->SetLightOutput(runLightOutput);
    analysisClusterCuts[i]->InitializeCutsFromCutString((cuts.GetClusterCut(i)).Data());
    ClusterCutList->Add(analysisClusterCuts[i]);
    analysisClusterCuts[i]->SetExtendedMatchAndQA(enableExtMatchAndQA);
    analysisClusterCuts[i]->SetFillCutHistograms("");
    
    analysisMesonCuts[i] = new AliConversionMesonCuts();
    analysisMesonCuts[i]->SetLightOutput(runLightOutput);
    analysisMesonCuts[i]->InitializeCutsFromCutString((cuts.GetMesonCut(i)).Data());
    analysisMesonCuts[i]->SetIsMergedClusterCut(2);
    analysisMesonCuts[i]->SetCaloMesonCutsObject(analysisClusterCuts[i]);
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
  task->SetProduceTreeEOverP(doTreeEOverP);
  task->SetEnableSortingOfMCClusLabels(enableSortingMCLabels);
  if(enableExtMatchAndQA > 1){ task->SetPlotHistsExtQA(kTRUE);}
  if(trainConfig == 106 || trainConfig == 125 || trainConfig == 145 || trainConfig == 206){
    task->SetInOutTimingCluster(-30e-9,35e-9);
  }
  task->SetLocalDebugFlag(localDebugFlag);
  
  //connect containers
  AliAnalysisDataContainer *coutput =
    mgr->CreateContainer(Form("GammaCalo_%i",trainConfig), TList::Class(),
              AliAnalysisManager::kOutputContainer,Form("GammaCalo_%i.root",trainConfig));

  mgr->AddTask(task);
  mgr->ConnectInput(task,0,cinput);
  mgr->ConnectOutput(task,1,coutput);

  return;

}
