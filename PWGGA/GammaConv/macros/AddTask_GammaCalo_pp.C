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
  TString corrTaskSetting = ""; // sets which correction task setting to use
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
      }else if(tempStr.BeginsWith("CF")){
        cout << "INFO: AddTask_GammaCalo_pp will use custom branch from Correction Framework!" << endl;
        corrTaskSetting = tempStr;
        corrTaskSetting.Replace(0,2,"");
      }
    }
  }
  TString sAdditionalTrainConfig = rAdditionalTrainConfig->GetString();
  if (sAdditionalTrainConfig.Atoi() > 0){
    trainConfig = trainConfig + sAdditionalTrainConfig.Atoi();
    cout << "INFO: AddTask_GammaCalo_pp running additionalTrainConfig '" << sAdditionalTrainConfig.Atoi() << "', train config: '" << trainConfig << "'" << endl;
  }
  cout << "corrTaskSetting: " << corrTaskSetting.Data() << endl;

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
  TString cutnumberPhoton = "00000008400000000100000000";
  if (  periodNameV0Reader.CompareTo("LHC16f") == 0 || periodNameV0Reader.CompareTo("LHC17g")==0 || periodNameV0Reader.CompareTo("LHC18c")==0 ||
        periodNameV0Reader.CompareTo("LHC17d1") == 0  || periodNameV0Reader.CompareTo("LHC17d12")==0 ||
        periodNameV0Reader.CompareTo("LHC17h3")==0 || periodNameV0Reader.CompareTo("LHC17k1")==0 ||
        periodNameV0Reader.CompareTo("LHC17f8b") == 0 ||
        periodNameV0Reader.CompareTo("LHC16P1JJLowB") == 0 || periodNameV0Reader.CompareTo("LHC16P1Pyt8LowB") == 0 )
    cutnumberPhoton         = "00000088400000000100000000";

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
  AliAnalysisTaskGammaCalo *task=NULL;
  task= new AliAnalysisTaskGammaCalo(Form("GammaCalo_%i",trainConfig));
  task->SetIsHeavyIon(isHeavyIon);
  task->SetIsMC(isMC);
  task->SetV0ReaderName(V0ReaderName);
  task->SetLightOutput(runLightOutput);
  task->SetCorrectionTaskSetting(corrTaskSetting);

  //create cut handler
  CutHandlerCalo cuts;

  // here is the order of the cluster cut string
  // usually for EMCal we start with 11111: default values for            "ClusterType", "EtaMin", "EtaMax", "PhiMin", "PhiMax"
  // then two numbers for nonlinearity, e.g. 21: this is                  "NonLinearity1", "NonLinearity2"
  // Then some cuts on the clusters, e.g. 06003222: this is               "DistanceToBadChannel", "Timing", "TrackMatching", "ExoticCell", "MinEnergy", "MinNCells", "MinM02", "MaxM02"
  // finally some for now unused cuts, usually 0000: this is              "MinMaxM20", "RecConv", "MaximumDispersion", "NLM"


  // *****************************************************************************************************
  // ******************** pp 2.76 TeV cuts paper EMC *****************************************************
  // *****************************************************************************************************
  if (trainConfig == 1){ // pp 2.76 TeV LHC11a paper cuts
    cuts.AddCut("00003013","1111121057032220000","0163103100000050"); // MB w/o pileup
    cuts.AddCut("00003113","1111121057032220000","0163103100000050"); // MB
    cuts.AddCut("00051013","1111121057032220000","0163103100000050"); // EMC1
  } else if (trainConfig == 2){  // pp 2.76 TeV LHC13g paper cuts
    cuts.AddCut("00010113","1111121067032220000","0163103100000050");
    cuts.AddCut("00010013","1111121067032220000","0163103100000050"); // without pile-up correction
    cuts.AddCut("00052013","1111121067032220000","0163103100000050"); // EMC7
    cuts.AddCut("00083013","1111121067032220000","0163103100000050"); // EMCEG1,
    cuts.AddCut("00085013","1111121067032220000","0163103100000050"); // EMCEG2,

  } else if (trainConfig == 3){
    cuts.AddCut("00003113","1111121057032250000","0163103100000050"); // MB
    cuts.AddCut("00003113","1111121057032260000","0163103100000050"); // MB
    cuts.AddCut("00003113","1111121057032240000","0163103100000050"); // MB
    cuts.AddCut("00003113","1111121057032290000","0163103100000050"); // MB
  } else if (trainConfig == 4){
    cuts.AddCut("00003113","11111210570322a0000","0163103100000050"); // MB
    cuts.AddCut("00003113","11111210570322b0000","0163103100000050"); // MB
    cuts.AddCut("00003113","11111210570322c0000","0163103100000050"); // MB
  } else if (trainConfig == 5){
    cuts.AddCut("00003113","1111121057032250000","01631031000000d0"); // MB
    cuts.AddCut("00003113","1111121057032260000","01631031000000d0"); // MB
    cuts.AddCut("00003113","1111121057032240000","01631031000000d0"); // MB
    cuts.AddCut("00003113","1111121057032290000","01631031000000d0"); // MB
  } else if (trainConfig == 6){
    cuts.AddCut("00003113","11111210570322a0000","01631031000000d0"); // MB
    cuts.AddCut("00003113","11111210570322b0000","01631031000000d0"); // MB
    cuts.AddCut("00003113","11111210570322c0000","01631031000000d0"); // MB

  // *****************************************************************************************************
  // ************************************* Calibration configuration EMC *********************************
  // *****************************************************************************************************
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
  // *****************************************************************************************************
  // ************************************* Direct Photon Configurations  *********************************
  // *****************************************************************************************************
  } else if (trainConfig == 60){
    cuts.AddCut("00003113","11111210570322c0000","01631031000000d0"); // MB std
    cuts.AddCut("00003113","1111121057032200000","01631031000000d0"); // no M02
    cuts.AddCut("00003113","11111210570322k0000","01631031000000d0"); // M02, pT-dep with 0.27-0.5
    cuts.AddCut("00003113","11111210570322l0000","01631031000000d0"); // M02, pT-dep with 0.32-0.5
    cuts.AddCut("00003113","11111210570322m0000","01631031000000d0"); // M02, pT-dep with 0.32-0.7
  } else if (trainConfig == 61){
    cuts.AddCut("00003113","11111210570322n0000","01631031000000d0"); // M02, pT-dep with 0.32-0.7
    cuts.AddCut("00003113","11111210570322o0000","01631031000000d0"); // M02, pT-dep with 0.27-0.7
    cuts.AddCut("00003113","11111210570322p0000","01631031000000d0"); // M02, pT-dep with 0.32-0.7
    cuts.AddCut("00003113","11111210570322q0000","01631031000000d0"); // M02, pT-dep with 0.34-0.7
  } else if (trainConfig == 62){
    cuts.AddCut("00051013","11111210570322c0000","01631031000000d0"); // MB std
    cuts.AddCut("00051013","1111121057032200000","01631031000000d0"); // no M02
    cuts.AddCut("00051013","11111210570322k0000","01631031000000d0"); // M02, pT-dep with 0.27-0.5
    cuts.AddCut("00051013","11111210570322l0000","01631031000000d0"); // M02, pT-dep with 0.32-0.5
    cuts.AddCut("00051013","11111210570322m0000","01631031000000d0"); // M02, pT-dep with 0.32-0.7
  } else if (trainConfig == 63){
    cuts.AddCut("00051013","11111210570322n0000","01631031000000d0"); // M02, pT-dep with 0.32-0.7
    cuts.AddCut("00051013","11111210570322o0000","01631031000000d0"); // M02, pT-dep with 0.27-0.7
    cuts.AddCut("00051013","11111210570322p0000","01631031000000d0"); // M02, pT-dep with 0.32-0.7
    cuts.AddCut("00051013","11111210570322q0000","01631031000000d0"); // M02, pT-dep with 0.34-0.7
  } else if (trainConfig == 64){
    cuts.AddCut("00003113","11111210570322l0000","01631031000000d0"); // M02 pt dep with new std: 0.32, 0.0072, 0.5
    cuts.AddCut("00003113","11111210570322d0000","01631031000000d0"); // M02, pt dep with  0.27, 0.0072, 0.4
    cuts.AddCut("00003113","11111210570322e0000","01631031000000d0"); // M02, pt dep with  0.31, 0.0072, 0.5
    cuts.AddCut("00003113","11111210570322f0000","01631031000000d0"); // M02, pt dep with  0.36, 0.0072, 0.7
  } else if (trainConfig == 65){
    cuts.AddCut("00003113","11111210570322m0000","01631031000000d0"); // M02, pt dep with  0.32, 0.0152, 0.5
    cuts.AddCut("00003113","11111210570322g0000","01631031000000d0"); // M02, pt dep with  0.37, 0.0072, 0.7
    cuts.AddCut("00003113","11111210570322h0000","01631031000000d0"); // M02, pT-dep with  0.30, 0.0072, 0.5
    cuts.AddCut("00003113","11111210570322i0000","01631031000000d0"); // M02, pT-dep with  0.35, 0.0072, 0.7
  } else if (trainConfig == 66){
    cuts.AddCut("00003113","11111210570322j0000","01631031000000d0"); // M02, pT-dep with  0.25, 0.0072, 0.39
    cuts.AddCut("00003113","11111210570322r0000","01631031000000d0"); // M02, pT-dep with  0.25, 0.0072, 0.5
    cuts.AddCut("00003113","11111210570322s0000","01631031000000d0"); // M02, pT-dep with  0.32, 0.0238, 0.7
    cuts.AddCut("00003113","11111210570322n0000","01631031000000d0"); // M02, pT-dep with  0.32, 0.0238, 0.7
  // *****************************************************************************************************
  // 8 TeV configs
  // *****************************************************************************************************
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

    // ATTENTION: adapted for 8 TeV dirGamma - ADJUSTED M02 -> l
    // 8 TeV variations
  } else if (trainConfig == 102){ // EMCAL clusters pp 8 TeV, timing+minEnergy variation
    cuts.AddCut("00010113","11111110570322l0000","01631031000000d0"); // time -50ns_50ns
    cuts.AddCut("00010113","11111110770322l0000","01631031000000d0"); // time -30ns_30ns
    cuts.AddCut("00010113","11111110870322l0000","01631031000000d0"); // time -20ns_30ns
    cuts.AddCut("00010113","11111110670122l0000","01631031000000d0"); //0.5 GeV/c
    cuts.AddCut("00010113","11111110670222l0000","01631031000000d0"); //0.6 GeV/c
    cuts.AddCut("00010113","11111110670422l0000","01631031000000d0"); //0.8 GeV/c
    cuts.AddCut("00010113","11111110670522l0000","01631031000000d0"); //0.9 GeV/c
  } else if (trainConfig == 103){ //EMCAL M02 variation
    cuts.AddCut("00010113","1111111067032200000","01631031000000d0"); //no max M02 cut
    cuts.AddCut("00010113","11111110670322l0000","01631031000000d0"); // M02 pt dep with new std: 0.32, 0.0072, 0.5
    cuts.AddCut("00010113","11111110670322d0000","01631031000000d0"); // M02, pt dep with  0.27, 0.0072, 0.4
    cuts.AddCut("00010113","11111110670322e0000","01631031000000d0"); // M02, pt dep with  0.31, 0.0072, 0.5
    cuts.AddCut("00010113","11111110670322f0000","01631031000000d0"); // M02, pt dep with  0.36, 0.0072, 0.7
    cuts.AddCut("00010113","11111110670322m0000","01631031000000d0"); // M02, pt dep with  0.32, 0.0152, 0.5
    cuts.AddCut("00010113","11111110670322g0000","01631031000000d0"); // M02, pt dep with  0.37, 0.0072, 0.7
    cuts.AddCut("00010113","11111110670322h0000","01631031000000d0"); // M02, pT-dep with  0.30, 0.0072, 0.5
    cuts.AddCut("00010113","11111110670322i0000","01631031000000d0"); // M02, pT-dep with  0.35, 0.0072, 0.7
    cuts.AddCut("00010113","11111110670322j0000","01631031000000d0"); // M02, pT-dep with  0.25, 0.0072, 0.39
    cuts.AddCut("00010113","11111110670322r0000","01631031000000d0"); // M02, pT-dep with  0.25, 0.0072, 0.5
    cuts.AddCut("00010113","11111110670322s0000","01631031000000d0"); // M02, pT-dep with  0.32, 0.0238, 0.7
    cuts.AddCut("00010113","11111110670322n0000","01631031000000d0"); // M02, pT-dep with  0.32, 0.0238, 0.7
  } else if (trainConfig == 104){ //EMCAL minNCells,with/without TRD variation
    cuts.AddCut("00010113","11111110670312l0000","01631031000000d0"); //n cells >= 1
    cuts.AddCut("00010113","11111110670332l0000","01631031000000d0"); //n cells >= 3
    cuts.AddCut("00010113","11121110670322l0000","01631031000000d0"); //only modules with TRD infront
    cuts.AddCut("00010113","11113110670322l0000","01631031000000d0"); //no modules with TRD infront
  } else if (trainConfig == 105){  // trackMatching variations
    cuts.AddCut("00010113","11111110660322l0000","01631031000000d0"); //
    cuts.AddCut("00010113","11111110680322l0000","01631031000000d0"); //
    cuts.AddCut("00010113","11111110690322l0000","01631031000000d0"); //
    cuts.AddCut("00010113","11111110630322l0000","01631031000000d0"); //
    cuts.AddCut("00010113","11111110600322l0000","01631031000000d0"); //
  } else if (trainConfig == 106){ // EMCAL clusters pp 8 TeV, combining cluster within time window and without
    cuts.AddCut("00010113","11111110070322l0000","01631031000000d0"); //
  } else if (trainConfig == 107){ // EMCAL clusters open angle variation
    cuts.AddCut("00010113","11111110670322l0000","01631031000000a0"); // min open angle - 0. + 1 cell diag
    cuts.AddCut("00010113","11111110670322l0000","01631031000000d0"); // min open angle - 0.017 + 1 cell diag
    cuts.AddCut("00010113","11111110670322l0000","01631031000000b0"); // min open angle - 0.0152 + 1 cell diag
    cuts.AddCut("00010113","11111110670322l0000","01631031000000c0"); // min open angle - 0.016 + 1 cell diag
    cuts.AddCut("00010113","11111110670322l0000","01631031000000e0"); // min open angle - 0.018 + 1 cell diag
    cuts.AddCut("00010113","11111110670322l0000","01631031000000f0"); // min open angle - 0.019 + 1 cell diag
  } else if (trainConfig == 108){ // EMCAL clusters pp 8 TeV, Different DistanceToBadChannels
    cuts.AddCut("00010113","11111111670322l0000","01631031000000d0"); //
    cuts.AddCut("00010113","11111112670322l0000","01631031000000d0"); //
    cuts.AddCut("00010113","11111113670322l0000","01631031000000d0"); //
    cuts.AddCut("00010113","11111115670322l0000","01631031000000d0"); //
    cuts.AddCut("00010113","11111116670322l0000","01631031000000d0"); //
  } else if (trainConfig == 109){ // EMCAL clusters pp 8 TeV, Different NonLinearities
    cuts.AddCut("00010113","11111010670322l0000","01631031000000d0"); // NonLinearity kSDMv5
    cuts.AddCut("00010113","11111130670322l0000","01631031000000d0"); // NonLinearity kTestBeamv2 + LHC12 ConvCalo
    cuts.AddCut("00010113","11111140670322l0000","01631031000000d0"); // NonLinearity kTestBeamv2 + LHC12 Calo
  } else if (trainConfig == 110){ // EMCAL clusters pp 8 TeV, Different NonLinearities
    cuts.AddCut("00010113","11111110670322l0000","01631031000000d0"); // NonLinearity LHC12 ConvCalo
    cuts.AddCut("00010113","11111120670322l0000","01631031000000d0"); // NonLinearity LHC12 Calo
    cuts.AddCut("00010113","11111210670322l0000","01631031000000d0"); // NonLinearity LHC12 ConvCalo MassRatioFits
    cuts.AddCut("00010113","11111220670322l0000","01631031000000d0"); // NonLinearity LHC12 Calo MassRatioFits
    cuts.AddCut("00010113","11111000670322l0000","01631031000000d0"); // NonLinearity none
  } else if (trainConfig == 111){  // EMCAL clusters, different triggers no NonLinearity
    cuts.AddCut("00010113","11111000670322l0000","01631031000000d0");
    cuts.AddCut("00052113","11111000670322l0000","01631031000000d0"); // EMC7
    cuts.AddCut("00081113","11111000670322l0000","01631031000000d0"); // EMCEG1,
  } else if (trainConfig == 112){ // EMCAL clusters, exotic cut var
    cuts.AddCut("00010113","11111110672322l0000","01631031000000d0"); //
    cuts.AddCut("00010113","11111110673322l0000","01631031000000d0"); //
    cuts.AddCut("00010113","11111110675322l0000","01631031000000d0"); //
    cuts.AddCut("00010113","11111110677322l0000","01631031000000d0"); //
    cuts.AddCut("00010113","11111110679322l0000","01631031000000d0"); //
  } else if (trainConfig == 113){  // trackMatching variations
    cuts.AddCut("00010113","11111110620322l0000","01631031000000d0"); //
    cuts.AddCut("00010113","11111110630322l0000","01631031000000d0"); //
    cuts.AddCut("00010113","11111110640322l0000","01631031000000d0"); //
    cuts.AddCut("00010113","11111110650322l0000","01631031000000d0"); //
  } else if (trainConfig == 114){ // EMCAL clusters pp 8 TeV
    cuts.AddCut("00010013","11111110670322l0000","01631031000000d0"); // std - no pileup cut
  } else if (trainConfig == 115){ // EMCAL clusters pp 8 TeV
    cuts.AddCut("00010113","11111110670322l0000","01631051000000d0"); // alpha
    cuts.AddCut("00010113","11111110670322l0000","01631061000000d0"); //
  } else if (trainConfig == 116){ // EMCAL clusters pp 8 TeV
    cuts.AddCut("00010113","11111110670322l0000","01633031000000d0"); // rapidity
    cuts.AddCut("00010113","11111110670322l0000","01634031000000d0"); //


  } else if (trainConfig == 118){ // EMCAL clusters pp 8 TeV - no SPD PileUp
    cuts.AddCut("00010113","11111110670322l0000","01631031000000d0"); // std
    cuts.AddCut("00010013","11111110670322l0000","01631031000000d0"); // std - no pileup cut
    cuts.AddCut("00052113","11111110670322l0000","01631031000000d0"); // std
    cuts.AddCut("00052013","11111110670322l0000","01631031000000d0"); // std - no pileup cut
    cuts.AddCut("00081113","11111110670322l0000","01631031000000d0"); // std
    cuts.AddCut("00081013","11111110670322l0000","01631031000000d0"); // std - no pileup cut

  // only std cuts
  } else if (trainConfig == 119){ // EMCAL clusters pp 8 TeV
    cuts.AddCut("00010113","1111111067032220000","01631031000000d0"); // std
  } else if (trainConfig == 120){ // EMCAL clusters pp 8 TeV
    cuts.AddCut("00010113","11111110670322l0000","01631031000000d0"); // std - dirGamma

    // ATTENTION: adapted for 8 TeV dirGamma  - ADJUSTED M02 Cut 2 -> p
    //8 TeV kEMC7 variations
  } else if (trainConfig == 121){ // EMCAL clusters pp 8 TeV, timing+minEnergy variation
    cuts.AddCut("00052113","11111110570322l0000","01631031000000d0"); // time -50ns_50ns
    cuts.AddCut("00052113","11111110770322l0000","01631031000000d0"); // time -30ns_30ns
    cuts.AddCut("00052113","11111110870322l0000","01631031000000d0"); // time -20ns_30ns
    cuts.AddCut("00052113","11111110670122l0000","01631031000000d0"); //0.5 GeV/c
    cuts.AddCut("00052113","11111110670222l0000","01631031000000d0"); //0.6 GeV/c
    cuts.AddCut("00052113","11111110670422l0000","01631031000000d0"); //0.8 GeV/c
    cuts.AddCut("00052113","11111110670522l0000","01631031000000d0"); //0.9 GeV/c
  } else if (trainConfig == 122){ //EMCAL M02 variation
    cuts.AddCut("00052113","1111111067032200000","01631031000000d0"); //no max M02 cut
    cuts.AddCut("00052113","11111110670322l0000","01631031000000d0"); // M02 pt dep with new std: 0.32, 0.0072, 0.5
    cuts.AddCut("00052113","11111110670322d0000","01631031000000d0"); // M02, pt dep with  0.27, 0.0072, 0.4
    cuts.AddCut("00052113","11111110670322e0000","01631031000000d0"); // M02, pt dep with  0.31, 0.0072, 0.5
    cuts.AddCut("00052113","11111110670322f0000","01631031000000d0"); // M02, pt dep with  0.36, 0.0072, 0.7
    cuts.AddCut("00052113","11111110670322m0000","01631031000000d0"); // M02, pt dep with  0.32, 0.0152, 0.5
    cuts.AddCut("00052113","11111110670322g0000","01631031000000d0"); // M02, pt dep with  0.37, 0.0072, 0.7
    cuts.AddCut("00052113","11111110670322h0000","01631031000000d0"); // M02, pT-dep with  0.30, 0.0072, 0.5
    cuts.AddCut("00052113","11111110670322i0000","01631031000000d0"); // M02, pT-dep with  0.35, 0.0072, 0.7
    cuts.AddCut("00052113","11111110670322j0000","01631031000000d0"); // M02, pT-dep with  0.25, 0.0072, 0.39
    cuts.AddCut("00052113","11111110670322r0000","01631031000000d0"); // M02, pT-dep with  0.25, 0.0072, 0.5
    cuts.AddCut("00052113","11111110670322s0000","01631031000000d0"); // M02, pT-dep with  0.32, 0.0238, 0.7
    cuts.AddCut("00052113","11111110670322n0000","01631031000000d0"); // M02, pT-dep with  0.32, 0.0238, 0.7
  } else if (trainConfig == 123){ //EMCAL minNCells, with/without TRD variation
    cuts.AddCut("00052113","11111110670312l0000","01631031000000d0"); //n cells >= 1
    cuts.AddCut("00052113","11111110670332l0000","01631031000000d0"); //n cells >= 3
    cuts.AddCut("00052113","11121110670322l0000","01631031000000d0"); //only modules with TRD infront
    cuts.AddCut("00052113","11113110670322l0000","01631031000000d0"); //no modules with TRD infront
  } else if (trainConfig == 124){  // trackMatching variations
    cuts.AddCut("00052113","11111110660322l0000","01631031000000d0"); //
    cuts.AddCut("00052113","11111110670322l0000","01631031000000d0"); //
    cuts.AddCut("00052113","11111110680322l0000","01631031000000d0"); //
    cuts.AddCut("00052113","11111110690322l0000","01631031000000d0"); //
    cuts.AddCut("00052113","11111110630322l0000","01631031000000d0"); //
    cuts.AddCut("00052113","11111110600322l0000","01631031000000d0"); //
  } else if (trainConfig == 125){ // EMCAL clusters pp 8 TeV, combining cluster within time window and without
    cuts.AddCut("00052113","11111110070322l0000","01631031000000d0"); //
  } else if (trainConfig == 126){ // EMCAL clusters open angle variation
    cuts.AddCut("00052113","11111110670322l0000","01631031000000a0"); // min open angle - 0. + 1 cell diag
    cuts.AddCut("00052113","11111110670322l0000","01631031000000d0"); // min open angle - 0.017 + 1 cell diag
    cuts.AddCut("00052113","11111110670322l0000","01631031000000b0"); // min open angle - 0.0152 + 1 cell diag
    cuts.AddCut("00052113","11111110670322l0000","01631031000000c0"); // min open angle - 0.016 + 1 cell diag
    cuts.AddCut("00052113","11111110670322l0000","01631031000000e0"); // min open angle - 0.018 + 1 cell diag
    cuts.AddCut("00052113","11111110670322l0000","01631031000000f0"); // min open angle - 0.019 + 1 cell diag
  } else if (trainConfig == 127){ // EMCAL clusters pp 8 TeV, Different DistanceToBadChannels
    cuts.AddCut("00052113","11111110670322l0000","01631031000000d0"); //
    cuts.AddCut("00052113","11111111670322l0000","01631031000000d0"); //
    cuts.AddCut("00052113","11111112670322l0000","01631031000000d0"); //
    cuts.AddCut("00052113","11111113670322l0000","01631031000000d0"); //
    cuts.AddCut("00052113","11111115670322l0000","01631031000000d0"); //
    cuts.AddCut("00052113","11111116670322l0000","01631031000000d0"); //
  } else if (trainConfig == 128){ // EMCAL clusters pp 8 TeV, Different NonLinearities
    cuts.AddCut("00052113","11111010670322l0000","01631031000000d0"); // NonLinearity kSDMv5
    cuts.AddCut("00052113","11111130670322l0000","01631031000000d0"); // NonLinearity kTestBeamv2 + LHC12 ConvCalo
    cuts.AddCut("00052113","11111140670322l0000","01631031000000d0"); // NonLinearity kTestBeamv2 + LHC12 Calo
  } else if (trainConfig == 129){ // EMCAL clusters pp 8 TeV, Different NonLinearities
    cuts.AddCut("00052113","11111110670322l0000","01631031000000d0"); // NonLinearity LHC12 ConvCalo
    cuts.AddCut("00052113","11111120670322l0000","01631031000000d0"); // NonLinearity LHC12 Calo
    cuts.AddCut("00052113","11111210670322l0000","01631031000000d0"); // NonLinearity LHC12 ConvCalo MassRatioFits
    cuts.AddCut("00052113","11111220670322l0000","01631031000000d0"); // NonLinearity LHC12 Calo MassRatioFits
    cuts.AddCut("00052113","11111000670322l0000","01631031000000d0"); // NonLinearity none
  } else if (trainConfig == 130){ // EMCAL clusters, exotic cut var
    cuts.AddCut("00052113","11111110670322l0000","01631031000000d0"); //
    cuts.AddCut("00052113","11111110672322l0000","01631031000000d0"); //
    cuts.AddCut("00052113","11111110673322l0000","01631031000000d0"); //
    cuts.AddCut("00052113","11111110675322l0000","01631031000000d0"); //
    cuts.AddCut("00052113","11111110677322l0000","01631031000000d0"); //
    cuts.AddCut("00052113","11111110679322l0000","01631031000000d0"); //
  } else if (trainConfig == 131){  // trackMatching variations
    cuts.AddCut("00052113","11111110670322l0000","01631031000000d0"); //
    cuts.AddCut("00052113","11111110620322l0000","01631031000000d0"); //
    cuts.AddCut("00052113","11111110630322l0000","01631031000000d0"); //
    cuts.AddCut("00052113","11111110640322l0000","01631031000000d0"); //
  } else if (trainConfig == 132){ // EMCAL clusters pp 8 TeV
    cuts.AddCut("00052013","11111110670322l0000","01631031000000d0"); // std - no pileup cut
  } else if (trainConfig == 133){ // EMCAL clusters pp 8 TeV
    cuts.AddCut("00052113","11111110670322l0000","01631051000000d0"); //
    cuts.AddCut("00052113","11111110670322l0000","01631061000000d0"); //
  } else if (trainConfig == 134){ // EMCAL clusters pp 8 TeV
    cuts.AddCut("00052113","11111110670322l0000","01633031000000d0"); //
    cuts.AddCut("00052113","11111110670322l0000","01634031000000d0"); //

  // only std cuts
  } else if (trainConfig == 139){ // EMCAL clusters pp 8 TeV
    cuts.AddCut("00052113","1111111067032220000","01631031000000d0"); // std
  } else if (trainConfig == 140){ // EMCAL clusters pp 8 TeV
    cuts.AddCut("00052113","11111110670322l0000","01631031000000d0"); // std - dirGAMMA

    // ATTENTION: adapted for 8 TeV dirGamma - ADJUSTED M02 Cut 2 -> p
    //8 TeV kEMCEGA variations
  } else if (trainConfig == 141){ // EMCAL clusters pp 8 TeV, timing+minEnergy variation
    cuts.AddCut("00081113","11111110570322l0000","01631031000000d0"); // time -50ns_50ns
    cuts.AddCut("00081113","11111110770322l0000","01631031000000d0"); // time -30ns_30ns
    cuts.AddCut("00081113","11111110870322l0000","01631031000000d0"); // time -20ns_30ns
    cuts.AddCut("00081113","11111110670122l0000","01631031000000d0"); //0.5 GeV/c
    cuts.AddCut("00081113","11111110670222l0000","01631031000000d0"); //0.6 GeV/c
    cuts.AddCut("00081113","11111110670422l0000","01631031000000d0"); //0.8 GeV/c
    cuts.AddCut("00081113","11111110670522l0000","01631031000000d0"); //0.9 GeV/c
  } else if (trainConfig == 142){ //EMCAL M02 variation
    cuts.AddCut("00081113","1111111067032200000","01631031000000d0"); //no max M02 cut
    cuts.AddCut("00081113","11111110670322l0000","01631031000000d0"); //M02, pT-dep
    cuts.AddCut("00081113","11111110670322d0000","01631031000000d0"); //M02, pT-dep
    cuts.AddCut("00081113","11111110670322e0000","01631031000000d0"); //M02, pT-dep
    cuts.AddCut("00081113","11111110670322f0000","01631031000000d0"); //M02, pT-dep
    cuts.AddCut("00081113","11111110670322m0000","01631031000000d0"); //M02, pT-dep
    cuts.AddCut("00081113","11111110670322g0000","01631031000000d0"); //M02, pT-dep
    cuts.AddCut("00081113","11111110670322h0000","01631031000000d0"); //M02, pT-dep
    cuts.AddCut("00081113","11111110670322i0000","01631031000000d0"); //M02, pT-dep
    cuts.AddCut("00081113","11111110670322j0000","01631031000000d0"); //M02, pT-dep
    cuts.AddCut("00081113","11111110670322r0000","01631031000000d0"); //M02, pT-dep
    cuts.AddCut("00081113","11111110670322s0000","01631031000000d0"); //M02, pT-dep
    cuts.AddCut("00081113","11111110670322n0000","01631031000000d0"); //M02, pT-dep
  } else if (trainConfig == 143){ //EMCAL minNCells, with/without TRD variation
    cuts.AddCut("00081113","11111110670312l0000","01631031000000d0"); //n cells >= 1
    cuts.AddCut("00081113","11111110670332l0000","01631031000000d0"); //n cells >= 3
    cuts.AddCut("00081113","11121110670322l0000","01631031000000d0"); //only modules with TRD infront
    cuts.AddCut("00081113","11113110670322l0000","01631031000000d0"); //no modules with TRD infront
  } else if (trainConfig == 144){  // trackMatching variations
    cuts.AddCut("00081113","11111110660322l0000","01631031000000d0"); //
    cuts.AddCut("00081113","11111110670322l0000","01631031000000d0"); //
    cuts.AddCut("00081113","11111110680322l0000","01631031000000d0"); //
    cuts.AddCut("00081113","11111110690322l0000","01631031000000d0"); //
    cuts.AddCut("00081113","11111110630322l0000","01631031000000d0"); //
    cuts.AddCut("00081113","11111110600322l0000","01631031000000d0"); //
  } else if (trainConfig == 145){ // EMCAL clusters pp 8 TeV, combining cluster within time window and without
    cuts.AddCut("00081113","11111110070322l0000","01631031000000d0"); //
  } else if (trainConfig == 146){ // EMCAL clusters open angle variation
    cuts.AddCut("00081113","11111110670322l0000","01631031000000a0"); // min open angle - 0. + 1 cell diag
    cuts.AddCut("00081113","11111110670322l0000","01631031000000d0"); // min open angle - 0.017 + 1 cell diag
    cuts.AddCut("00081113","11111110670322l0000","01631031000000b0"); // min open angle - 0.0152 + 1 cell diag
    cuts.AddCut("00081113","11111110670322l0000","01631031000000c0"); // min open angle - 0.016 + 1 cell diag
    cuts.AddCut("00081113","11111110670322l0000","01631031000000e0"); // min open angle - 0.018 + 1 cell diag
    cuts.AddCut("00081113","11111110670322l0000","01631031000000f0"); // min open angle - 0.019 + 1 cell diag
  } else if (trainConfig == 147){ // EMCAL clusters pp 8 TeV, Different DistanceToBadChannels
    cuts.AddCut("00081113","11111110670322l0000","01631031000000d0"); //
    cuts.AddCut("00081113","11111111670322l0000","01631031000000d0"); //
    cuts.AddCut("00081113","11111112670322l0000","01631031000000d0"); //
    cuts.AddCut("00081113","11111113670322l0000","01631031000000d0"); //
    cuts.AddCut("00081113","11111115670322l0000","01631031000000d0"); //
    cuts.AddCut("00081113","11111116670322l0000","01631031000000d0"); //
  } else if (trainConfig == 148){ // EMCAL clusters pp 8 TeV, Different NonLinearities
    cuts.AddCut("00081113","11111010670322l0000","01631031000000d0"); // NonLinearity kSDMv5
    cuts.AddCut("00081113","11111130670322l0000","01631031000000d0"); // NonLinearity kTestBeamv2 + LHC12 ConvCalo
    cuts.AddCut("00081113","11111140670322l0000","01631031000000d0"); // NonLinearity kTestBeamv2 + LHC12 Calo
  } else if (trainConfig == 149){ // EMCAL clusters pp 8 TeV, Different NonLinearities
    cuts.AddCut("00081113","11111110670322l0000","01631031000000d0"); // NonLinearity LHC12 ConvCalo
    cuts.AddCut("00081113","11111120670322l0000","01631031000000d0"); // NonLinearity LHC12 Calo
    cuts.AddCut("00081113","11111210670322l0000","01631031000000d0"); // NonLinearity LHC12 ConvCalo MassRatioFits
    cuts.AddCut("00081113","11111220670322l0000","01631031000000d0"); // NonLinearity LHC12 Calo MassRatioFits
    cuts.AddCut("00081113","11111000670322l0000","01631031000000d0"); // NonLinearity none
  } else if (trainConfig == 150){ // EMCAL clusters, exotic cut var
    cuts.AddCut("00081113","11111110670322l0000","01631031000000d0"); //
    cuts.AddCut("00081113","11111110672322l0000","01631031000000d0"); //
    cuts.AddCut("00081113","11111110673322l0000","01631031000000d0"); //
    cuts.AddCut("00081113","11111110675322l0000","01631031000000d0"); //
    cuts.AddCut("00081113","11111110677322l0000","01631031000000d0"); //
    cuts.AddCut("00081113","11111110679322l0000","01631031000000d0"); //
  } else if (trainConfig == 151){  // trackMatching variations
    cuts.AddCut("00081113","11111110670322l0000","01631031000000d0"); //
    cuts.AddCut("00081113","11111110620322l0000","01631031000000d0"); //
    cuts.AddCut("00081113","11111110630322l0000","01631031000000d0"); //
    cuts.AddCut("00081113","11111110640322l0000","01631031000000d0"); //
  } else if (trainConfig == 152){ // EMCAL clusters pp 8 TeV
    cuts.AddCut("00081013","11111110670322l0000","01631031000000d0"); // std - no pileup cut
  } else if (trainConfig == 153){ // EMCAL clusters pp 8 TeV
    cuts.AddCut("00081113","11111110670322l0000","01631051000000d0"); //
    cuts.AddCut("00081113","11111110670322l0000","01631061000000d0"); //
  } else if (trainConfig == 154){ // EMCAL clusters pp 8 TeV
    cuts.AddCut("00081113","11111110670322l0000","01633031000000d0"); //
    cuts.AddCut("00081113","11111110670322l0000","01634031000000d0"); //

  // only std cuts
  } else if (trainConfig == 159){ // EMCAL clusters pp 8 TeV
    cuts.AddCut("00081113","1111111067032220000","01631031000000d0"); // std
  } else if (trainConfig == 160){ // EMCAL clusters pp 8 TeV
    cuts.AddCut("00081113","11111110670322l0000","01631031000000d0"); // std - dirGAMMA

  // CutStudies for DirGamma
  } else if (trainConfig == 161){ //
    cuts.AddCut("00010113","11111110670322l0000","01631031000000d0"); // std
    cuts.AddCut("00052113","11111110670322l0000","01631031000000d0"); // std
    cuts.AddCut("00081113","11111110670322l0000","01631031000000d0"); // std

  } else if (trainConfig == 162){ //
    cuts.AddCut("00010113","11111110670322l0000","01631031000000d0"); // std
  } else if (trainConfig == 163){ //
    cuts.AddCut("00052113","11111110670322l0000","01631031000000d0"); // std
  } else if (trainConfig == 164){ //
    cuts.AddCut("00081113","11111110670322l0000","01631031000000d0"); // std

  } else if (trainConfig == 165){
    cuts.AddCut("00010113","11111110670322l0000","01631031000000d0"); // M02 pt dep with new std: 0.32, 0.0072, 0.5
    cuts.AddCut("00010113","11111110670322d0000","01631031000000d0"); // M02, pt dep with  0.27, 0.0072, 0.4
    cuts.AddCut("00010113","11111110670322e0000","01631031000000d0"); // M02, pt dep with  0.31, 0.0072, 0.5
    cuts.AddCut("00010113","11111110670322f0000","01631031000000d0"); // M02, pt dep with  0.36, 0.0072, 0.7
  } else if (trainConfig == 166){
    cuts.AddCut("00010113","11111110670322m0000","01631031000000d0"); // M02, pt dep with  0.32, 0.0152, 0.5
    cuts.AddCut("00010113","11111110670322g0000","01631031000000d0"); // M02, pt dep with  0.37, 0.0072, 0.7
    cuts.AddCut("00010113","11111110670322h0000","01631031000000d0"); // M02, pT-dep with  0.30, 0.0072, 0.5
    cuts.AddCut("00010113","11111110670322i0000","01631031000000d0"); // M02, pT-dep with  0.35, 0.0072, 0.7
  } else if (trainConfig == 167){
    cuts.AddCut("00010113","11111110670322j0000","01631031000000d0"); // M02, pT-dep with  0.25, 0.0072, 0.39
    cuts.AddCut("00010113","11111110670322r0000","01631031000000d0"); // M02, pT-dep with  0.25, 0.0072, 0.5
    cuts.AddCut("00010113","11111110670322s0000","01631031000000d0"); // M02, pT-dep with  0.32, 0.0238, 0.7
    cuts.AddCut("00010113","11111110670322n0000","01631031000000d0"); // M02, pT-dep with  0.32, 0.0238, 0.7

  //multiple std dirGAMMA cuts for different studies
  } else if (trainConfig == 181){ // EMCAL clusters pp 8 TeV
    cuts.AddCut("00010113","11111110670322l0000","01631031000000d0"); // std
  } else if (trainConfig == 182){ // EMCAL clusters pp 8 TeV
    cuts.AddCut("00010113","11111110670322l0000","01631031000000d0"); // std
  } else if (trainConfig == 183){ // EMCAL clusters pp 8 TeV
    cuts.AddCut("00010113","11111110670322l0000","01631031000000d0"); // std
  } else if (trainConfig == 184){ // EMCAL clusters pp 8 TeV
    cuts.AddCut("00010113","11111110670322l0000","01631031000000d0"); // std
  } else if (trainConfig == 185){ // EMCAL clusters pp 8 TeV
    cuts.AddCut("00010113","11111110670322l0000","01631031000000d0"); // std
  } else if (trainConfig == 186){ // EMCAL clusters pp 8 TeV
    cuts.AddCut("00010113","11111110670322l0000","01631031000000d0"); // std
  } else if (trainConfig == 187){ // EMCAL clusters pp 8 TeV
    cuts.AddCut("00010113","11111110670322l0000","01631031000000d0"); // std
  } else if (trainConfig == 188){ // EMCAL clusters pp 8 TeV
    cuts.AddCut("00010113","11111110670322l0000","01631031000000d0"); // std
  } else if (trainConfig == 189){ // EMCAL clusters pp 8 TeV
    cuts.AddCut("00010113","11111110670322l0000","01631031000000d0"); // std
  } else if (trainConfig == 190){ // EMCAL clusters pp 8 TeV
    cuts.AddCut("00010113","11111110670322l0000","01631031000000d0"); // std
  } else if (trainConfig == 191){ // EMCAL clusters pp 8 TeV
    cuts.AddCut("00010113","11111110670322l0000","01631031000000d0"); // std
  } else if (trainConfig == 192){ // EMCAL clusters pp 8 TeV
    cuts.AddCut("00010113","11111110670322l0000","01631031000000d0"); // std

  // 7 TeV
  } else if (trainConfig == 200){ // EMCAL clusters pp 7 TeV, pT dep matching
    cuts.AddCut("00000113","11111110b7032220000","01631031000000d0"); // std
    cuts.AddCut("00000113","1111111007032220000","01631031000000d0"); // std
  } else if (trainConfig == 201){ // EMCAL clusters pp 7 TeV, pT dep matching
    cuts.AddCut("00000113","11111110b7032220000","01631031000000d0"); // std

    // ATTENTION: adapted for dirGamma - ADJUSTED M02 -> l
  } else if (trainConfig == 202){ // EMCAL clusters pp 7 TeV, timing+minEnergy variation
    cuts.AddCut("00000113","11111110370322l0000","01631031000000d0"); // time
    cuts.AddCut("00000113","11111110470322l0000","01631031000000d0"); // time
    cuts.AddCut("00000113","11111110570322l0000","01631031000000d0"); // time
    cuts.AddCut("00000113","11111110b70122l0000","01631031000000d0"); //0.5 GeV/c
    cuts.AddCut("00000113","11111110b70222l0000","01631031000000d0"); //0.6 GeV/c
    cuts.AddCut("00000113","11111110b70422l0000","01631031000000d0"); //0.8 GeV/c
    cuts.AddCut("00000113","11111110b70522l0000","01631031000000d0"); //0.9 GeV/c
  } else if (trainConfig == 203){ //EMCAL M02 variation
    cuts.AddCut("00000113","11111110b7032200000","01631031000000d0"); //no max M02 cut
    cuts.AddCut("00000113","11111110b70322l0000","01631031000000d0"); // M02 pt dep with new std: 0.32, 0.0072, 0.5
    cuts.AddCut("00000113","11111110b70322d0000","01631031000000d0"); // M02, pt dep with  0.27, 0.0072, 0.4
    cuts.AddCut("00000113","11111110b70322e0000","01631031000000d0"); // M02, pt dep with  0.31, 0.0072, 0.5
    cuts.AddCut("00000113","11111110b70322f0000","01631031000000d0"); // M02, pt dep with  0.36, 0.0072, 0.7
    cuts.AddCut("00000113","11111110b70322m0000","01631031000000d0"); // M02, pt dep with  0.32, 0.0152, 0.5
    cuts.AddCut("00000113","11111110b70322g0000","01631031000000d0"); // M02, pt dep with  0.37, 0.0072, 0.7
    cuts.AddCut("00000113","11111110b70322h0000","01631031000000d0"); // M02, pT-dep with  0.30, 0.0072, 0.5
    cuts.AddCut("00000113","11111110b70322i0000","01631031000000d0"); // M02, pT-dep with  0.35, 0.0072, 0.7
    cuts.AddCut("00000113","11111110b70322j0000","01631031000000d0"); // M02, pT-dep with  0.25, 0.0072, 0.39
    cuts.AddCut("00000113","11111110b70322r0000","01631031000000d0"); // M02, pT-dep with  0.25, 0.0072, 0.5
    cuts.AddCut("00000113","11111110b70322s0000","01631031000000d0"); // M02, pT-dep with  0.32, 0.0238, 0.7
    cuts.AddCut("00000113","11111110b70322n0000","01631031000000d0"); // M02, pT-dep with  0.32, 0.0238, 0.7
  } else if (trainConfig == 204){ //EMCAL minNCells,with/without TRD variation
    cuts.AddCut("00000113","11111110b70312l0000","01631031000000d0"); //n cells >= 1
    cuts.AddCut("00000113","11111110b70332l0000","01631031000000d0"); //n cells >= 3
    cuts.AddCut("00000113","11121110b70322l0000","01631031000000d0"); //only modules with TRD infront
    cuts.AddCut("00000113","11113110b70322l0000","01631031000000d0"); //no modules with TRD infront
  } else if (trainConfig == 205){  // trackMatching variations
    cuts.AddCut("00000113","11111110b60322l0000","01631031000000d0"); //
    cuts.AddCut("00000113","11111110b80322l0000","01631031000000d0"); //
    cuts.AddCut("00000113","11111110b90322l0000","01631031000000d0"); //
    cuts.AddCut("00000113","11111110b30322l0000","01631031000000d0"); //
    cuts.AddCut("00000113","11111110b00322l0000","01631031000000d0"); //
  } else if (trainConfig == 207){ // EMCAL clusters open angle variation
    cuts.AddCut("00000113","11111110b70322l0000","01631031000000a0"); // min open angle - 0. + 1 cell diag
    cuts.AddCut("00000113","11111110b70322l0000","01631031000000d0"); // min open angle - 0.017 + 1 cell diag
    cuts.AddCut("00000113","11111110b70322l0000","01631031000000b0"); // min open angle - 0.0152 + 1 cell diag
    cuts.AddCut("00000113","11111110b70322l0000","01631031000000c0"); // min open angle - 0.016 + 1 cell diag
    cuts.AddCut("00000113","11111110b70322l0000","01631031000000e0"); // min open angle - 0.018 + 1 cell diag
    cuts.AddCut("00000113","11111110b70322l0000","01631031000000f0"); // min open angle - 0.019 + 1 cell diag
  } else if (trainConfig == 208){ // EMCAL clusters pp 7 TeV, Different DistanceToBadChannels
    cuts.AddCut("00000113","11111111b70322l0000","01631031000000d0"); //
    cuts.AddCut("00000113","11111112b70322l0000","01631031000000d0"); //
    cuts.AddCut("00000113","11111113b70322l0000","01631031000000d0"); //
    cuts.AddCut("00000113","11111115b70322l0000","01631031000000d0"); //
    cuts.AddCut("00000113","11111116b70322l0000","01631031000000d0"); //
  } else if (trainConfig == 209){ // EMCAL clusters pp 7 TeV, Different NonLinearities
    cuts.AddCut("00000113","11111010b70322l0000","01631031000000d0"); // NonLinearity kSDMv5
    cuts.AddCut("00000113","11111130b70322l0000","01631031000000d0"); // NonLinearity kTestBeamv2 + LHC12 ConvCalo
    cuts.AddCut("00000113","11111140b70322l0000","01631031000000d0"); // NonLinearity kTestBeamv2 + LHC12 Calo
  } else if (trainConfig == 210){ // EMCAL clusters pp 7 TeV, Different NonLinearities
    cuts.AddCut("00000113","11111110b70322l0000","01631031000000d0"); // NonLinearity LHC12 ConvCalo
    cuts.AddCut("00000113","11111120b70322l0000","01631031000000d0"); // NonLinearity LHC12 Calo
    cuts.AddCut("00000113","11111210b70322l0000","01631031000000d0"); // NonLinearity LHC12 ConvCalo MassRatioFits
    cuts.AddCut("00000113","11111220b70322l0000","01631031000000d0"); // NonLinearity LHC12 Calo MassRatioFits
    cuts.AddCut("00000113","11111000b70322l0000","01631031000000d0"); // NonLinearity none
  } else if (trainConfig == 211){  // EMCAL clusters, different triggers no NonLinearity
    cuts.AddCut("00000113","11111000b70322l0000","01631031000000d0");
    cuts.AddCut("00052113","11111000b70322l0000","01631031000000d0"); // EMC7
    cuts.AddCut("00081113","11111000b70322l0000","01631031000000d0"); // EMCEG1,
  } else if (trainConfig == 212){ // EMCAL clusters, exotic cut var
    cuts.AddCut("00000113","11111110b72322l0000","01631031000000d0"); //
    cuts.AddCut("00000113","11111110b73322l0000","01631031000000d0"); //
    cuts.AddCut("00000113","11111110b75322l0000","01631031000000d0"); //
    cuts.AddCut("00000113","11111110b77322l0000","01631031000000d0"); //
    cuts.AddCut("00000113","11111110b79322l0000","01631031000000d0"); //
  } else if (trainConfig == 213){  // trackMatching variations
    cuts.AddCut("00000113","11111110b20322l0000","01631031000000d0"); //
    cuts.AddCut("00000113","11111110b30322l0000","01631031000000d0"); //
    cuts.AddCut("00000113","11111110b40322l0000","01631031000000d0"); //
    cuts.AddCut("00000113","11111110b50322l0000","01631031000000d0"); //
  } else if (trainConfig == 214){ // EMCAL clusters pp 7 TeV
    cuts.AddCut("00010013","11111110b70322l0000","01631031000000d0"); // std - no pileup cut
  } else if (trainConfig == 215){ // EMCAL clusters pp 7 TeV
    cuts.AddCut("00000113","11111110b70322l0000","01631051000000d0"); // alpha
    cuts.AddCut("00000113","11111110b70322l0000","01631061000000d0"); //
  } else if (trainConfig == 216){ // EMCAL clusters pp 7 TeV
    cuts.AddCut("00000113","11111110b70322l0000","01633031000000d0"); // rapidity
    cuts.AddCut("00000113","11111110b70322l0000","01634031000000d0"); //


  // only std cuts
  } else if (trainConfig == 219){ // EMCAL clusters pp 7 TeV
    cuts.AddCut("00000113","11111110b7032220000","01631031000000d0"); // std
  } else if (trainConfig == 220){ // EMCAL clusters pp 7 TeV
    cuts.AddCut("00000113","11111110b70322l0000","01631031000000d0"); // std - dirGAMMA

    // std
    // ATTENTION: adapted for dirGamma - ADJUSTED M02 -> l
  } else if (trainConfig == 221){ // EMCAL clusters pp 7 TeV
    cuts.AddCut("00000113","11111110b70322l0000","01631031000000d0"); // std
  } else if (trainConfig == 222){ // EMCAL clusters pp 7 TeV
    cuts.AddCut("00000113","11111110b70322l0000","01631031000000d0"); // std
  // LHC11cd configs V0OR and V0AND
  } else if (trainConfig == 250){  // EMCAL clusters 7 TeV LHC11 TM on
    cuts.AddCut("00010113","11111000670322l0000","01631031000000d0"); // VOAND
    cuts.AddCut("00052113","11111000670322l0000","01631031000000d0"); // EMC7
    cuts.AddCut("00000113","11111000670322l0000","01631031000000d0"); // V0OR
    cuts.AddCut("00051113","11111000670322l0000","01631031000000d0"); // EMC1
  } else if (trainConfig == 251){  // EMCAL clusters 7 TeV LHC11 TM off
    cuts.AddCut("00010113","11111000600322l0000","01631031000000d0"); // VOAND
    cuts.AddCut("00052113","11111000600322l0000","01631031000000d0"); // EMC7
    cuts.AddCut("00000113","11111000600322l0000","01631031000000d0"); // V0OR
    cuts.AddCut("00051113","11111000600322l0000","01631031000000d0"); // EMC1
  } else if (trainConfig == 252){  // EMCAL clusters 7 TeV LHC11 TM on + TB NonLin
    cuts.AddCut("00010113","11111020670322l0000","01631031000000d0"); // VOAND
    cuts.AddCut("00052113","11111020670322l0000","01631031000000d0"); // EMC7
    cuts.AddCut("00000113","11111020670322l0000","01631031000000d0"); // V0OR
    cuts.AddCut("00051113","11111020670322l0000","01631031000000d0"); // EMC1
  } else if (trainConfig == 253){  // EMCAL clusters 7 TeV LHC11 TM off + TB NonLin
    cuts.AddCut("00010113","11111020600322l0000","01631031000000d0"); // VOAND
    cuts.AddCut("00052113","11111020600322l0000","01631031000000d0"); // EMC7
    cuts.AddCut("00000113","11111020600322l0000","01631031000000d0"); // V0OR
    cuts.AddCut("00051113","11111020600322l0000","01631031000000d0"); // EMC1
  } else if (trainConfig == 254){  // QA for settings of omega analysis
    cuts.AddCut("00000113","1111111047032230000","0163503900000000");

  // *****************************************************************************************************
  // ************************************* PHOS cuts ****************************************************
  // *****************************************************************************************************
  // pp 2.76 TeV
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

  // pp 7 TeV direct photon PHOS
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

  // pp 7 TeV PHOS
  } else if (trainConfig == 361){
    cuts.AddCut("00000113","2444400000013300000","0163803100000010"); // QA
    cuts.AddCut("00000113","2444400040013300000","0163803100000010"); // 100ns timing cut, no track matching
    cuts.AddCut("00000113","2444400043013300000","0163803100000010"); // 100ns timing cut
  } else if (trainConfig == 362){
    cuts.AddCut("00062113","2444400000013300000","0163803100000010"); // QA
    cuts.AddCut("00062113","2444400040013300000","0163803100000010"); // 100ns timing cut, no track matching
    cuts.AddCut("00062113","2444400043013300000","0163803100000010"); // 100ns timing cut
  } else if (trainConfig == 363){ // train config for bad channels and NonLin Variation
    cuts.AddCut("00000113","2444400000013300000","0163803100000010"); // no NonLin
    cuts.AddCut("00000113","2444401000013300000","0163803100000010"); // extern PHOS NonLin
    cuts.AddCut("00000113","2444412000013300000","0163803100000010"); // constant non Lin first iteration
  // PHOS @ 8 TeV
  } else if (trainConfig == 381){ // PHOS clusters
    cuts.AddCut("00010113","2444400040013300000","0163103100000010");
  } else if (trainConfig == 382){ // PHOS clusters
    cuts.AddCut("00062113","2444400040013300000","0163103100000010");

  // *********************************************************************************************************
  // 5 TeV  pp Run2 - EMC configurations
  // *********************************************************************************************************
  } else if (trainConfig == 401){ // EMCAL clusters
    cuts.AddCut("00010113","1111100017032220000","01631031000000d0"); // 1000ns timing cut, no NL INT7
    cuts.AddCut("00052013","1111100017032220000","01631031000000d0"); // 1000ns timing cut, no NL EMC7
    cuts.AddCut("00085013","1111100017032220000","01631031000000d0"); // 1000ns timing cut, no NL EG2
    cuts.AddCut("00083013","1111100017032220000","01631031000000d0"); // 1000ns timing cut, no NL EG1
  } else if (trainConfig == 402){ // EMCAL clusters
    cuts.AddCut("00010113","1111100067032220000","01631031000000d0"); // -50ns, 30ns timing cut, no NL INT7
    cuts.AddCut("00052013","1111100067032220000","01631031000000d0"); // -50ns, 30ns timing cut, no NL EMC7
    cuts.AddCut("00085013","1111100067032220000","01631031000000d0"); // -50ns, 30ns timing cut, no NL EG2
    cuts.AddCut("00083013","1111100067032220000","01631031000000d0"); // -50ns, 30ns timing cut, no NL EG1
  } else if (trainConfig == 403){ // EMCAL clusters - NonLin INT7
    cuts.AddCut("00010113","1111111017032220000","01631031000000d0"); // 1000ns timing cut, NL kSDM PCMEMC
    cuts.AddCut("00010113","1111112017032220000","01631031000000d0"); // 1000ns timing cut, NL kSDM EMC
    cuts.AddCut("00010113","1111121017032220000","01631031000000d0"); // 1000ns timing cut, NL DExt PCMEMC
    cuts.AddCut("00010113","1111122017032220000","01631031000000d0"); // 1000ns timing cut, NL DExt EMC
  } else if (trainConfig == 404){ // EMCAL clusters - NonLin INT7
    cuts.AddCut("00010113","1111111067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL kSDM PCMEMC
    cuts.AddCut("00010113","1111112067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL kSDM EMC
    cuts.AddCut("00010113","1111121067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL DExt PCMEMC
    cuts.AddCut("00010113","1111122067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL DExt EMC
  } else if (trainConfig == 405){ // EMCAL clusters - NonLin EMC7
    cuts.AddCut("00052013","1111111017032220000","01631031000000d0"); // 1000ns timing cut, NL kSDM PCMEMC
    cuts.AddCut("00052013","1111112017032220000","01631031000000d0"); // 1000ns timing cut, NL kSDM EMC
    cuts.AddCut("00052013","1111121017032220000","01631031000000d0"); // 1000ns timing cut, NL DExt PCMEMC
    cuts.AddCut("00052013","1111122017032220000","01631031000000d0"); // 1000ns timing cut, NL DExt EMC
  } else if (trainConfig == 406){ // EMCAL clusters - NonLin EMC7
    cuts.AddCut("00052013","1111111067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL kSDM PCMEMC
    cuts.AddCut("00052013","1111112067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL kSDM EMC
    cuts.AddCut("00052013","1111121067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL DExt PCMEMC
    cuts.AddCut("00052013","1111122067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL DExt EMC
  } else if (trainConfig == 407){ // EMCAL clusters - NonLin EG2
    cuts.AddCut("00085013","1111111017032220000","01631031000000d0"); // 1000ns timing cut, NL kSDM PCMEMC
    cuts.AddCut("00085013","1111112017032220000","01631031000000d0"); // 1000ns timing cut, NL kSDM EMC
    cuts.AddCut("00085013","1111121017032220000","01631031000000d0"); // 1000ns timing cut, NL DExt PCMEMC
    cuts.AddCut("00085013","1111122017032220000","01631031000000d0"); // 1000ns timing cut, NL DExt EMC
  } else if (trainConfig == 408){ // EMCAL clusters - NonLin EG2
    cuts.AddCut("00085013","1111111067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL kSDM PCMEMC
    cuts.AddCut("00085013","1111112067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL kSDM EMC
    cuts.AddCut("00085013","1111121067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL DExt PCMEMC
    cuts.AddCut("00085013","1111122067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL DExt EMC
  } else if (trainConfig == 409){ // EMCAL clusters - NonLin EG1
    cuts.AddCut("00083013","1111111017032220000","01631031000000d0"); // 1000ns timing cut, NL kSDM PCMEMC
    cuts.AddCut("00083013","1111112017032220000","01631031000000d0"); // 1000ns timing cut, NL kSDM EMC
    cuts.AddCut("00083013","1111121017032220000","01631031000000d0"); // 1000ns timing cut, NL DExt PCMEMC
    cuts.AddCut("00083013","1111122017032220000","01631031000000d0"); // 1000ns timing cut, NL DExt EMC
  } else if (trainConfig == 410){ // EMCAL clusters - NonLin EG1
    cuts.AddCut("00083013","1111111067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL kSDM PCMEMC
    cuts.AddCut("00083013","1111112067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL kSDM EMC
    cuts.AddCut("00083013","1111121067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL DExt PCMEMC
    cuts.AddCut("00083013","1111122067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL DExt EMC
  } else if (trainConfig == 411){ // EMCAL clusters - NonLin INT7
    cuts.AddCut("00010113","1111113017032220000","01631031000000d0"); // 1000ns timing cut, NL kSDM PCMEMC + testbeam
    cuts.AddCut("00010113","1111114017032220000","01631031000000d0"); // 1000ns timing cut, NL kSDM EMC + testbeam
    cuts.AddCut("00010113","1111123017032220000","01631031000000d0"); // 1000ns timing cut, NL DExt PCMEMC + testbeam
    cuts.AddCut("00010113","1111124017032220000","01631031000000d0"); // 1000ns timing cut, NL DExt EMC + testbeam
  } else if (trainConfig == 412){ // EMCAL clusters - NonLin INT7
    cuts.AddCut("00010113","1111113067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL kSDM PCMEMC + testbeam
    cuts.AddCut("00010113","1111114067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL kSDM EMC + testbeam
    cuts.AddCut("00010113","1111123067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL DExt PCMEMC + testbeam
    cuts.AddCut("00010113","1111124067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL DExt EMC + testbeam
  } else if (trainConfig == 413){ // EMCAL clusters - NonLin EMC7
    cuts.AddCut("00052013","1111113017032220000","01631031000000d0"); // 1000ns timing cut, NL kSDM PCMEMC
    cuts.AddCut("00052013","1111114017032220000","01631031000000d0"); // 1000ns timing cut, NL kSDM EMC
    cuts.AddCut("00052013","1111123017032220000","01631031000000d0"); // 1000ns timing cut, NL DExt PCMEMC
    cuts.AddCut("00052013","1111124017032220000","01631031000000d0"); // 1000ns timing cut, NL DExt EMC
  } else if (trainConfig == 414){ // EMCAL clusters - NonLin EMC7
    cuts.AddCut("00052013","1111113067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL kSDM PCMEMC
    cuts.AddCut("00052013","1111114067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL kSDM EMC
    cuts.AddCut("00052013","1111123067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL DExt PCMEMC
    cuts.AddCut("00052013","1111124067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL DExt EMC
  } else if (trainConfig == 415){ // EMCAL clusters - NonLin INT7
    cuts.AddCut("00010113","1111102017032220000","01631031000000d0"); // 1000ns timing cut, NL BeamTest ONLY
    cuts.AddCut("00010113","1111109017032220000","01631031000000d0"); // 1000ns timing cut, NL 8 TeV 11
    cuts.AddCut("00010113","1111117017032220000","01631031000000d0"); // 1000ns timing cut, NL kSDM PCMEMC - fixed E-dependence
    cuts.AddCut("00010113","1111127017032220000","01631031000000d0"); // 1000ns timing cut, NL DExt PCMEMC - fixed E-dependence
  } else if (trainConfig == 416){ // EMCAL clusters - NonLin INT7
    cuts.AddCut("00010113","1111102067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL BeamTest ONLY
    cuts.AddCut("00010113","1111109067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL 8 TeV 11
    cuts.AddCut("00010113","1111117067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL kSDM PCMEMC - fixed E-dependence
    cuts.AddCut("00010113","1111127067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL DExt PCMEMC - fixed E-dependence
  } else if (trainConfig == 417){ // EMCAL clusters - NonLin INT7 - with BeamTest
    cuts.AddCut("00010113","1111113017032220000","01631031000000d0"); // 1000ns timing cut, NL kSDM PCMEMC+BTv3
    cuts.AddCut("00010113","1111114017032220000","01631031000000d0"); // 1000ns timing cut, NL kSDM PCMEMC+BTv3
    cuts.AddCut("00010113","1111123017032220000","01631031000000d0"); // 1000ns timing cut, NL DExt PCMEMC+BTv3
    cuts.AddCut("00010113","1111124017032220000","01631031000000d0"); // 1000ns timing cut, NL DExt PCMEMC+BTv3
  } else if (trainConfig == 418){ // EMCAL clusters - NonLin INT7 - with BeamTest
    cuts.AddCut("00010113","1111113067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL kSDM PCMEMC+BTv3
    cuts.AddCut("00010113","1111114067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL kSDM PCMEMC+BTv3
    cuts.AddCut("00010113","1111123067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL DExt PCMEMC+BTv3
    cuts.AddCut("00010113","1111124067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL DExt PCMEMC+BTv3
    // 5 TeV JetJet configs without trackmatching
  } else if (trainConfig == 430){ // EMCAL clusters no NonLin
    cuts.AddCut("00010113","1111100060032220000","01631031000000d0"); // -50ns, 30ns timing cut, no NL INT7
    cuts.AddCut("00052013","1111100060032220000","01631031000000d0"); // -50ns, 30ns timing cut, no NL EMC7
    cuts.AddCut("00085013","1111100060032220000","01631031000000d0"); // -50ns, 30ns timing cut, no NL EG2
    cuts.AddCut("00083013","1111100060032220000","01631031000000d0"); // -50ns, 30ns timing cut, no NL EG1
  } else if (trainConfig == 431){ // EMCAL clusters - NonLin INT7 - no trackmatching for 5TeV JJ
    cuts.AddCut("00010113","1111111060032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL kSDM PCMEMC
    cuts.AddCut("00010113","1111112060032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL kSDM EMC
    cuts.AddCut("00010113","1111121060032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL DExt PCMEMC
    cuts.AddCut("00010113","1111122060032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL DExt EMC
  } else if (trainConfig == 432){ // EMCAL clusters - NonLin EMC7 - no trackmatching for 5TeV JJ
    cuts.AddCut("00052013","1111111060032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL kSDM PCMEMC
    cuts.AddCut("00052013","1111112060032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL kSDM EMC
    cuts.AddCut("00052013","1111121060032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL DExt PCMEMC
    cuts.AddCut("00052013","1111122060032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL DExt EMC
  } else if (trainConfig == 433){ // EMCAL clusters - NonLin EG2  - no trackmatching for 5TeV JJ
    cuts.AddCut("00085013","1111111060032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL kSDM PCMEMC
    cuts.AddCut("00085013","1111112060032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL kSDM EMC
    cuts.AddCut("00085013","1111121060032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL DExt PCMEMC
    cuts.AddCut("00085013","1111122060032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL DExt EMC
  } else if (trainConfig == 434){ // EMCAL clusters - NonLin EG1  - no trackmatching for 5TeV JJ
    cuts.AddCut("00083013","1111111060032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL kSDM PCMEMC
    cuts.AddCut("00083013","1111112060032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL kSDM EMC
    cuts.AddCut("00083013","1111121060032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL DExt PCMEMC
    cuts.AddCut("00083013","1111122060032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL DExt EMC

  } else if (trainConfig == 450){ // EMCAL standard cuts, different triggers
    cuts.AddCut("00010113","1111111067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL kSDM PCMEMC INT7
    cuts.AddCut("00052013","1111111067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL kSDM PCMEMC
    cuts.AddCut("00085013","1111111067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL kSDM PCMEMC EG2
    cuts.AddCut("00083013","1111111067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL kSDM PCMEMC EG1
  } else if (trainConfig == 451){ // EMCAL syst 1/7
    cuts.AddCut("00010113","1111111067022220000","01631031000000d0"); // min energy cluster variation 1  600 MeV
    cuts.AddCut("00010113","1111111067042220000","01631031000000d0"); // min energy cluster variation 2  800 MeV
    cuts.AddCut("00010113","1111111067052220000","01631031000000d0"); // min energy cluster variation 3  900 MeV
    cuts.AddCut("00010113","1111111067032220000","01631061000000d0"); // alpha meson variation 1 0<alpha<0.8
    cuts.AddCut("00010113","1111111067032220000","01631051000000d0"); // alpha meson variation 2 0<alpha<0.75
  } else if (trainConfig == 452){ // EMCAL syst 2/7
    cuts.AddCut("00010113","1111111067032230000","01631031000000d0"); // min/max M02  0.1<M<0.5
    cuts.AddCut("00010113","1111111067032200000","01631031000000d0"); // min/max M02  0.1<M<100
    cuts.AddCut("00010113","1111111067032250000","01631031000000d0"); // min/max M02  0.1<M<0.3
    cuts.AddCut("00010113","1111111067032260000","01631031000000d0"); // min/max M02  0.1<M<0.27
  } else if (trainConfig == 453){ // EMCAL syst 3/7
    cuts.AddCut("00010113","1112111067032220000","01631031000000d0"); // only modules with TRD infront
    cuts.AddCut("00010113","1111311067032220000","01631031000000d0"); // no modules with TRD infront
    cuts.AddCut("00010113","1111111067032220000","01633031000000d0"); // rapidity variation  y<0.6
    cuts.AddCut("00010113","1111111067032220000","01634031000000d0"); // rapidity variation  y<0.5
  } else if (trainConfig == 454){ // EMCAL syst 4/7
    cuts.AddCut("00010113","1111111067032220000","0163103100000040"); // min opening angle 0.0152
    cuts.AddCut("00010113","1111111067032220000","0163103100000060"); // min opening angle 0.017
    cuts.AddCut("00010113","1111111067032220000","0163103100000070"); // min opening angle 0.016
    cuts.AddCut("00010113","1111111067032220000","0163103100000080"); // min opening angle 0.018
    cuts.AddCut("00010113","1111111067032220000","0163103100000090"); // min opening angle 0.019
  } else if (trainConfig == 455){ // EMCAL syst 5/7
    cuts.AddCut("00010113","1111111063032220000","01631031000000d0"); // fixed window
    cuts.AddCut("00010113","1111111066032220000","01631031000000d0"); // tm pt dependent var 1
    cuts.AddCut("00010113","1111111068032220000","01631031000000d0"); // tm pt dependent var 2
    cuts.AddCut("00010113","1111111069032220000","01631031000000d0"); // tm pt dependent var 3
  } else if (trainConfig == 456){ // EMCAL syst 6/7
    cuts.AddCut("00010113","1111111037032230000","01631031000000d0"); // cluster timing cut
    cuts.AddCut("00010113","1111111047032230000","01631031000000d0"); // cluster timing cut
    cuts.AddCut("00010113","1111111057032230000","01631031000000d0"); // cluster timing cut
    cuts.AddCut("00010113","1111111077032230000","01631031000000d0"); // cluster timing cut
  } else if (trainConfig == 457){ // EMCAL syst 7/7
    cuts.AddCut("00010113","1111111087032230000","01631031000000d0"); // cluster timing cut
    cuts.AddCut("00010113","1111111097032230000","01631031000000d0"); // cluster timing cut
    cuts.AddCut("00010113","11111110a7032230000","01631031000000d0"); // cluster timing cut

  } else if (trainConfig == 460){ // INT7 EMCAL standard cut but with E/p TM veto
    cuts.AddCut("00010113","1111111067032220000","01631031000000d0"); // std INT7
    cuts.AddCut("00010113","111111106c032220000","01631031000000d0"); // fEOverPMax = 9e9
    cuts.AddCut("00010113","111111106d032220000","01631031000000d0"); // fEOverPMax = 3.0
    cuts.AddCut("00010113","111111106e032220000","01631031000000d0"); // fEOverPMax = 2.0
  } else if (trainConfig == 461){ // INT7 EMCAL standard cut but with E/p TM veto
    cuts.AddCut("00010113","111111106f032220000","01631031000000d0"); // fEOverPMax = 1.75
    cuts.AddCut("00010113","111111106g032220000","01631031000000d0"); // fEOverPMax = 1.5
    cuts.AddCut("00010113","111111106h032220000","01631031000000d0"); // fEOverPMax = 1.25
  } else if (trainConfig == 462){ // EMC7 EMCAL standard cut but with E/p TM veto
    cuts.AddCut("00052113","1111111067032220000","01631031000000d0"); // std EMC7
    cuts.AddCut("00052113","111111106c032220000","01631031000000d0"); // fEOverPMax = 9e9
    cuts.AddCut("00052113","111111106d032220000","01631031000000d0"); // fEOverPMax = 3.0
    cuts.AddCut("00052113","111111106e032220000","01631031000000d0"); // fEOverPMax = 2.0
  } else if (trainConfig == 463){ // EMC7 EMCAL standard cut but with E/p TM veto
    cuts.AddCut("00052113","111111106f032220000","01631031000000d0"); // fEOverPMax = 1.75
    cuts.AddCut("00052113","111111106g032220000","01631031000000d0"); // fEOverPMax = 1.5
    cuts.AddCut("00052113","111111106h032220000","01631031000000d0"); // fEOverPMax = 1.25
  } else if (trainConfig == 464){ // EG1 EMCAL standard cut but with E/p TM veto
    cuts.AddCut("00083113","1111111067032220000","01631031000000d0"); // std EG1
    cuts.AddCut("00083113","111111106c032220000","01631031000000d0"); // fEOverPMax = 9e9
    cuts.AddCut("00083113","111111106d032220000","01631031000000d0"); // fEOverPMax = 3.0
    cuts.AddCut("00083113","111111106e032220000","01631031000000d0"); // fEOverPMax = 2.0
  } else if (trainConfig == 465){ // EG1 EMCAL standard cut but with E/p TM veto
    cuts.AddCut("00083113","111111106f032220000","01631031000000d0"); // fEOverPMax = 1.75
    cuts.AddCut("00083113","111111106g032220000","01631031000000d0"); // fEOverPMax = 1.5
    cuts.AddCut("00083113","111111106h032220000","01631031000000d0"); // fEOverPMax = 1.25
  } else if (trainConfig == 466){ // EG2 EMCAL standard cut but with E/p TM veto
    cuts.AddCut("00085113","1111111067032220000","01631031000000d0"); // std EG2
    cuts.AddCut("00085113","111111106c032220000","01631031000000d0"); // fEOverPMax = 9e9
    cuts.AddCut("00085113","111111106d032220000","01631031000000d0"); // fEOverPMax = 3.0
    cuts.AddCut("00085113","111111106e032220000","01631031000000d0"); // fEOverPMax = 2.0
  } else if (trainConfig == 467){ // EG2 EMCAL standard cut but with E/p TM veto
    cuts.AddCut("00085113","111111106f032220000","01631031000000d0"); // fEOverPMax = 1.75
    cuts.AddCut("00085113","111111106g032220000","01631031000000d0"); // fEOverPMax = 1.5
    cuts.AddCut("00085113","111111106h032220000","01631031000000d0"); // fEOverPMax = 1.25
    
  } else if (trainConfig == 480){ // INT7 EMCAL standard cut with triggers - NO TM - CALO+CALOFAST readout triggers
    cuts.AddCut("000a0113","1111111060032220000","01631031000000d0"); // std INT7
    cuts.AddCut("000a1113","1111111060032220000","01631031000000d0"); // std EMC7
    cuts.AddCut("000a2113","1111111060032220000","01631031000000d0"); // std EG2
    cuts.AddCut("000a3113","1111111060032220000","01631031000000d0"); // std EG1
  } else if (trainConfig == 481){ // INT7 EMCAL standard cut with triggers - NO TM
    cuts.AddCut("00010113","1111111060032220000","01631031000000d0"); // std INT7
    cuts.AddCut("00052113","1111111060032220000","01631031000000d0"); // std EMC7
    cuts.AddCut("00085113","1111111060032220000","01631031000000d0"); // std EG2
    cuts.AddCut("00083113","1111111060032220000","01631031000000d0"); // std EG1
  // *********************************************************************************************************
  // 13 TeV  pp Run2 - EMC configurations
  // *********************************************************************************************************
  } else if (trainConfig == 500){ // EMCAL clusters
    cuts.AddCut("00010113","1111100017032220000","01631031000000d0"); // 1000ns timing cut, no NL INT7
    cuts.AddCut("00010113","1111100067032220000","01631031000000d0"); // -30ns, 35ns timing cut, no NL INT7
  } else if (trainConfig == 501){ // EMCAL clusters
    cuts.AddCut("00052113","1111100017032220000","01631031000000d0"); // 1000ns timing cut, no NL EMC7
    cuts.AddCut("00085113","1111100017032220000","01631031000000d0"); // 1000ns timing cut, no NL EG2
    cuts.AddCut("00083113","1111100017032220000","01631031000000d0"); // 1000ns timing cut, no NL EG1
  } else if (trainConfig == 502){ // EMCAL clusters
    cuts.AddCut("00052113","1111100067032220000","01631031000000d0"); // -50ns, 30ns timing cut, no NL EMC7
    cuts.AddCut("00085113","1111100067032220000","01631031000000d0"); // -50ns, 30ns timing cut, no NL EG2
    cuts.AddCut("00083113","1111100067032220000","01631031000000d0"); // -50ns, 30ns timing cut, no NL EG1
  } else if (trainConfig == 503){ // EMCAL clusters
    cuts.AddCut("00010113","1111112067032220000","01631031000000d0"); // -30ns, 35ns timing cut, NL1 INT7
    cuts.AddCut("00010113","1111122067032220000","01631031000000d0"); // -30ns, 35ns timing cut, NL2 INT7
    cuts.AddCut("00010113","1111111067032220000","01631031000000d0"); // -30ns, 35ns timing cut, NL1 INT7
    cuts.AddCut("00010113","1111121067032220000","01631031000000d0"); // -30ns, 35ns timing cut, NL2 INT7
  } else if (trainConfig == 504){ // EMCAL clusters
    cuts.AddCut("00052113","1111112067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL1 EMC7
    cuts.AddCut("00052113","1111122067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL2 EMC7
    cuts.AddCut("00052113","1111111067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL3 EMC7
    cuts.AddCut("00052113","1111121067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL4 EMC7
  } else if (trainConfig == 505){ // EMCAL clusters
    cuts.AddCut("00085113","1111112067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL1 EG2
    cuts.AddCut("00085113","1111122067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL2 EG2
    cuts.AddCut("00085113","1111111067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL3 EG2
    cuts.AddCut("00085113","1111121067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL4 EG2
  } else if (trainConfig == 506){ // EMCAL clusters
    cuts.AddCut("00083113","1111112067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL1 EG1
    cuts.AddCut("00083113","1111122067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL2 EG1
    cuts.AddCut("00083113","1111121067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL4 EG1
    cuts.AddCut("00083113","1111111067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL3 EG1
  } else if (trainConfig == 507){ // EMCAL clusters
    cuts.AddCut("00010113","1111102067032220000","01631031000000d0"); // -30ns, 35ns timing cut, NLtestbeam INT7
    cuts.AddCut("00052113","1111102067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NLtestbeam EMC7
    cuts.AddCut("00085113","1111102067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NLtestbeam EG2
    cuts.AddCut("00083113","1111102067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NLtestbeam EG1



  // *********************************************************************************************************
  // 13 TeV  DMC configurations
  // *********************************************************************************************************
  } else if (trainConfig==510){ //DCAL
    cuts.AddCut("00010113","3885500017032220000","01631031000000d0"); // 1000ns timing cut, no NL INT7
    cuts.AddCut("00010113","3885500067032220000","01631031000000d0"); // -30ns, 35ns timing cut, no NL INT7
  } else if (trainConfig == 511){ // DCAL clusters
    cuts.AddCut("00055113","3885500017032220000","01631031000000d0"); // 1000ns timing cut, no NL DMC7
    cuts.AddCut("00089113","3885500017032220000","01631031000000d0"); // 1000ns timing cut, no NL DG2
    cuts.AddCut("0008b113","3885500017032220000","01631031000000d0"); // 1000ns timing cut, no NL DG1
  } else if (trainConfig == 512){ // DCAL clusters
    cuts.AddCut("00055113","3885500067032220000","01631031000000d0"); // -50ns, 30ns timing cut, no NL DMC7
    cuts.AddCut("00089113","3885500067032220000","01631031000000d0"); // -50ns, 30ns timing cut, no NL DG2
    cuts.AddCut("0008b113","3885500067032220000","01631031000000d0"); // -50ns, 30ns timing cut, no NL DG1

  // *********************************************************************************************************
  // 13 TeV  pp Run2 - EMC configurations HM trigg
  // *********************************************************************************************************
  } else if (trainConfig == 520){ // EMCAL clusters
    cuts.AddCut("00074113","1111100017032220000","01631031000000d0"); // 1000ns timing cut, no NL VOHM
    cuts.AddCut("00076113","1111100017032220000","01631031000000d0"); // 1000ns timing cut, no NL VOHM with SPD
    cuts.AddCut("00074113","1111100067032220000","01631031000000d0"); // -50ns, 30ns timing cut, no NL VOHM
    cuts.AddCut("00076113","1111100067032220000","01631031000000d0"); // -50ns, 30ns timing cut, no NL VOHM with SPD
  } else if (trainConfig == 521){ // EMCAL clusters
    cuts.AddCut("00074113","1111112067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL1 VOHM
    cuts.AddCut("00074113","1111122067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL2 VOHM
    cuts.AddCut("00074113","1111102067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NLtestbeam VOHM
    cuts.AddCut("00074113","1111111067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL3 VOHM
    cuts.AddCut("00074113","1111121067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL4 VOHM
  } else if (trainConfig == 522){ // EMCAL clusters
    cuts.AddCut("00076113","1111121067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL4 VOHM with SPD
    cuts.AddCut("00076113","1111111067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL3 VOHM with SPD
    cuts.AddCut("00076113","1111102067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NLtestbeam VOHM with SPD
    cuts.AddCut("00076113","1111122067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL2 VOHM with SPD
    cuts.AddCut("00076113","1111112067032220000","01631031000000d0"); // -50ns, 30ns timing cut, NL1 VOHM with SPD

  // *********************************************************************************************************
  // 13 TeV  DMC configurations HM trigg
  // *********************************************************************************************************
  } else if (trainConfig==530){ //DCAL
    cuts.AddCut("00074113","3885500017032220000","01631031000000d0"); // 1000ns timing cut, no NL VOHM
    cuts.AddCut("00074113","3885500067032220000","01631031000000d0"); // -50ns, 30ns timing cut, no NL VOHM
    cuts.AddCut("00076113","3885500017032220000","01631031000000d0"); // 1000ns timing cut, no NL VOHM with SPD
    cuts.AddCut("00076113","3885500067032220000","01631031000000d0"); // -50ns, 30ns timing cut, no NL VOHM with SPD

  // *********************************************************************************************************
  // 5 TeV 2015 pp Run2 - DMC configurations
  // *********************************************************************************************************
  //                 LHC17pq
  } else if (trainConfig == 600){
    cuts.AddCut("00010113","3885500087032220000","01631031000000d0"); // std
    cuts.AddCut("00010113","3885500017032220000","01631031000000d0"); // QA
  } else if (trainConfig == 601){ // DCAL clusters pp 5.02 TeV (Triggered)
    cuts.AddCut("00055113","3885500087032220000","01631031000000d0"); // EMC7
    // Changed BGEvents to 20(80) for variations
  } else if (trainConfig == 602){ // timing Cut variation  std -20+50ns
    cuts.AddCut("00010113","3885500017032220000","01631031000000d0"); //     -1000  +1000 ns
    cuts.AddCut("00010113","3885500077032220000","01631031000000d0"); //     -30    +30   ns
    cuts.AddCut("00010113","3885500097032220000","01631031000000d0"); //     -20    +25   ns
    cuts.AddCut("00010113","38855000a7032220000","01631031000000d0"); //     -12.5  +13   ns
  } else if (trainConfig == 603){ // track matching variation
    cuts.AddCut("00010113","3885500080032220000","01631031000000d0"); //
    cuts.AddCut("00010113","3885500081032220000","01631031000000d0"); //
    cuts.AddCut("00010113","3885500086032220000","01631031000000d0"); //
  } else if (trainConfig == 604){ // opening angle variation
    cuts.AddCut("00010113","3885500087032220000","0163103100000000"); //
    cuts.AddCut("00010113","3885500087032220000","0163103100000010"); //
    cuts.AddCut("00010113","3885500087032220000","0163103100000030"); //
    cuts.AddCut("00010113","3885500087032220000","0163103100000040"); //
    cuts.AddCut("00010113","3885500087032220000","0163103100000050"); //
    cuts.AddCut("00010113","3885500087032220000","0163103100000080"); //
  } else if (trainConfig == 605){ // min energy variation std 0.7 GeV/c
    cuts.AddCut("00010113","3885500087002220000","01631031000000d0"); //     0.1GeV/c
    cuts.AddCut("00010113","3885500087012220000","01631031000000d0"); //     0.5 GeV/c
    cuts.AddCut("00010113","3885500087022220000","01631031000000d0"); //     0.6 GeV/c
    cuts.AddCut("00010113","3885500087042220000","01631031000000d0"); //     0.8 GeV/c
    cuts.AddCut("00010113","3885500087052220000","01631031000000d0"); //     0.9 GeV/c
  } else if (trainConfig == 606){ // min nCells & M02 variation
    // std: min nCells = 1; M02 max=0.7, min=0.1
    cuts.AddCut("00010113","3885500087031220000","01631031000000d0"); //   min nCells = 1
    cuts.AddCut("00010113","3885500087033220000","01631031000000d0"); //   min nCells = 3
    cuts.AddCut("00010113","3885500087032210000","01631031000000d0"); //   max M02 = 1
    cuts.AddCut("00010113","3885500087032240000","01631031000000d0"); //   max M02 = 0.4
  } else if (trainConfig == 607){ // NonLin variation
    cuts.AddCut("00010113","3885500087032220000","01631031000000d0"); // no NL
    cuts.AddCut("00010113","3885511087032220000","01631031000000d0"); // PCM-DCal kSDM
    cuts.AddCut("00010113","3885512087032220000","01631031000000d0"); // DCal kSDM
    cuts.AddCut("00010113","3885521087032220000","01631031000000d0"); // PCM-DCal DExp/DPow
    cuts.AddCut("00010113","3885522087032220000","01631031000000d0"); // DCal DExp/DPow
  } else if (trainConfig == 608){ // LHC15n
    cuts.AddCut("00010113","3885522087032220000","01631031000000d0"); // std
    cuts.AddCut("00010113","3885500017032220000","01631031000000d0"); // QA

  } else if (trainConfig == 660){ // INT7 DCAL standard cut but with E/p TM veto
    cuts.AddCut("00010113","3885511087041220000","01631031000000d0"); // std INT7
    cuts.AddCut("00010113","388551108c041220000","01631031000000d0"); // fEOverPMax = 9e9
    cuts.AddCut("00010113","388551108d041220000","01631031000000d0"); // fEOverPMax = 3.0
    cuts.AddCut("00010113","388551108e041220000","01631031000000d0"); // fEOverPMax = 2.0
  } else if (trainConfig == 661){ // INT7 DCAL standard cut but with E/p TM veto
    cuts.AddCut("00010113","388551108f041220000","01631031000000d0"); // fEOverPMax = 1.75
    cuts.AddCut("00010113","388551108g041220000","01631031000000d0"); // fEOverPMax = 1.5
    cuts.AddCut("00010113","388551108h041220000","01631031000000d0"); // fEOverPMax = 1.25
  } else if (trainConfig == 662){ // DMC7 DCAL standard cut but with E/p TM veto
    cuts.AddCut("00055113","3885511087041220000","01631031000000d0"); // std DMC7
    cuts.AddCut("00055113","388551108c041220000","01631031000000d0"); // fEOverPMax = 9e9
    cuts.AddCut("00055113","388551108d041220000","01631031000000d0"); // fEOverPMax = 3.0
    cuts.AddCut("00055113","388551108e041220000","01631031000000d0"); // fEOverPMax = 2.0
  } else if (trainConfig == 663){ // DMC7 DCAL standard cut but with E/p TM veto
    cuts.AddCut("00055113","388551108f041220000","01631031000000d0"); // fEOverPMax = 1.75
    cuts.AddCut("00055113","388551108g041220000","01631031000000d0"); // fEOverPMax = 1.5
    cuts.AddCut("00055113","388551108h041220000","01631031000000d0"); // fEOverPMax = 1.25
  } else if (trainConfig == 664){ // DG1 DCAL standard cut but with E/p TM veto
    cuts.AddCut("0008b113","3885511087041220000","01631031000000d0"); // std DG1
    cuts.AddCut("0008b113","388551108c041220000","01631031000000d0"); // fEOverPMax = 9e9
    cuts.AddCut("0008b113","388551108d041220000","01631031000000d0"); // fEOverPMax = 3.0
    cuts.AddCut("0008b113","388551108e041220000","01631031000000d0"); // fEOverPMax = 2.0
  } else if (trainConfig == 665){ // DG1 DCAL standard cut but with E/p TM veto
    cuts.AddCut("0008b113","388551108f041220000","01631031000000d0"); // fEOverPMax = 1.75
    cuts.AddCut("0008b113","388551108g041220000","01631031000000d0"); // fEOverPMax = 1.5
    cuts.AddCut("0008b113","388551108h041220000","01631031000000d0"); // fEOverPMax = 1.25
  } else if (trainConfig == 666){ // DG2 DCAL standard cut but with E/p TM veto
    cuts.AddCut("00089113","3885511087041220000","01631031000000d0"); // std DG2
    cuts.AddCut("00089113","388551108c041220000","01631031000000d0"); // fEOverPMax = 9e9
    cuts.AddCut("00089113","388551108d041220000","01631031000000d0"); // fEOverPMax = 3.0
    cuts.AddCut("00089113","388551108e041220000","01631031000000d0"); // fEOverPMax = 2.0
  } else if (trainConfig == 667){ // DG2 DCAL standard cut but with E/p TM veto
    cuts.AddCut("00089113","388551108f041220000","01631031000000d0"); // fEOverPMax = 1.75
    cuts.AddCut("00089113","388551108g041220000","01631031000000d0"); // fEOverPMax = 1.5
    cuts.AddCut("00089113","388551108h041220000","01631031000000d0"); // fEOverPMax = 1.25

  } else if (trainConfig == 680){ // DCAL standard for CALO+CALOFAST readout triggers - NO TM
    cuts.AddCut("000a0113","3885511080041220000","01631031000000d0"); // std INT7
    cuts.AddCut("000a6113","3885511080041220000","01631031000000d0"); // std DMC7
    cuts.AddCut("000a7113","3885511080041220000","01631031000000d0"); // std DG2
    cuts.AddCut("000a8113","3885511080041220000","01631031000000d0"); // std DG1
  } else if (trainConfig == 681){ // DCAL standard for standard readout triggers - NO TM
    cuts.AddCut("00010113","3885511080041220000","01631031000000d0"); // std INT7
    cuts.AddCut("00055113","3885511080041220000","01631031000000d0"); // std DMC7
    cuts.AddCut("00089113","3885511080041220000","01631031000000d0"); // std DG2
    cuts.AddCut("0008b113","3885511080041220000","01631031000000d0"); // std DG1
  // *********************************************************************************************************
  // 5 TeV 2015 pp Run2 - PHOS configurations
  // *********************************************************************************************************
  } else if (trainConfig == 700){ // PHOS clusters with larger acceptance
    cuts.AddCut("00010113","2446600040013300000","0163103100000010"); // INT7
    cuts.AddCut("00062113","2446600040013300000","0163103100000010"); // PHI7
  } else if (trainConfig == 701){ // Default cut, No TM
    cuts.AddCut("00010113","2446651040013300000","0163103100000010"); // INT7
    cuts.AddCut("00062113","2446651040013300000","0163103100000010"); // PHI7
  } else if (trainConfig == 702){ // Default cut, with TM
    cuts.AddCut("00010113","2446651044013300000","0163103100000010"); // INT7
    cuts.AddCut("00062113","2446651044013300000","0163103100000010"); // PHI7
  } else if( trainConfig == 703){ // NL variations
    cuts.AddCut("00010113","2446651044013300000","0163103100000010"); // INT7
    cuts.AddCut("00010113","2446652044013300000","0163103100000010"); // PHOS calo NL
    cuts.AddCut("00010113","2446601044013300000","0163103100000010"); // PHOS people NL

  // *********************************************************************************************************
  // 13 TeV 2015 pp Run2 - PHOS configurations
  // *********************************************************************************************************
  } else if (trainConfig == 800){ // PHOS clusters with larger acceptance NCells 3
    cuts.AddCut("00010113","2446600040013300000","0163103100000010"); // INT7
    cuts.AddCut("00062113","2446600040013300000","0163103100000010"); // PHI7
  } else if (trainConfig == 801){ // PHOS clusters with larger acceptance NCells 2
    cuts.AddCut("00010113","2446600040012300000","0163103100000010"); // INT7
    cuts.AddCut("00062113","2446600040012300000","0163103100000010"); // PHI7
  } else if (trainConfig == 802){ // PHOS clusters with larger acceptance w/ TM NCells 3
    cuts.AddCut("00010113","2446600044013300000","0163103100000010"); // INT7
    cuts.AddCut("00062113","2446600044013300000","0163103100000010"); // PHI7
  } else if (trainConfig == 803){ // PHOS clusters with larger acceptance w/ TM NCells 2
    cuts.AddCut("00010113","2446600044012300000","0163103100000010"); // INT7
    cuts.AddCut("00062113","2446600044012300000","0163103100000010"); // PHI7

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
      fTrackMatcher->SetCorrectionTaskSetting(corrTaskSetting);
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
    analysisEventCuts[i]->SetCorrectionTaskSetting(corrTaskSetting);
    analysisEventCuts[i]->SetLightOutput(runLightOutput);
    if (periodNameV0Reader.CompareTo("") != 0) analysisEventCuts[i]->SetPeriodEnum(periodNameV0Reader);
    analysisEventCuts[i]->InitializeCutsFromCutString((cuts.GetEventCut(i)).Data());
    EventCutList->Add(analysisEventCuts[i]);
    analysisEventCuts[i]->SetFillCutHistograms("",kFALSE);

    analysisClusterCuts[i] = new AliCaloPhotonCuts(isMC);
    analysisClusterCuts[i]->SetHistoToModifyAcceptance(histoAcc);
    analysisClusterCuts[i]->SetV0ReaderName(V0ReaderName);
    analysisClusterCuts[i]->SetCorrectionTaskSetting(corrTaskSetting);
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
  task->SetCorrectionTaskSetting(corrTaskSetting);
  task->SetDoMesonQA(enableQAMesonTask); //Attention new switch for Pi0 QA
  task->SetDoClusterQA(enableQAClusterTask);  //Attention new switch small for Cluster QA
  task->SetDoTHnSparse(isUsingTHnSparse);
  task->SetProduceTreeEOverP(doTreeEOverP);
  task->SetEnableSortingOfMCClusLabels(enableSortingMCLabels);
  if(enableExtMatchAndQA > 1){ task->SetPlotHistsExtQA(kTRUE);}
  if(trainConfig == 106 || trainConfig == 125 || trainConfig == 145){
    task->SetInOutTimingCluster(-30e-9,35e-9);
  }
  task->SetLocalDebugFlag(localDebugFlag);

  //connect containers
  AliAnalysisDataContainer *coutput =
    mgr->CreateContainer(!(corrTaskSetting.CompareTo("")) ? Form("GammaCalo_%i",trainConfig) : Form("GammaCalo_%i_%s",trainConfig,corrTaskSetting.Data()), TList::Class(),
              AliAnalysisManager::kOutputContainer,Form("GammaCalo_%i.root",trainConfig));

  mgr->AddTask(task);
  mgr->ConnectInput(task,0,cinput);
  mgr->ConnectOutput(task,1,coutput);

  return;

}
