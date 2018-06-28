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
                            Int_t      trainConfig                  = 1,                // change different set of cuts
                            Int_t      isMC                         = 0,                // run MC
                            Int_t      enableQAMesonTask            = 0,                // enable QA in AliAnalysisTaskGammaConvV1
                            Int_t      enableQAClusterTask          = 0,                // enable additional QA task
                            TString    fileNameInputForWeighting    = "MCSpectraInput.root",       // path to file for weigting input / modified acceptance
                            Int_t      doWeightingPart              = 0,                // enable Weighting
                            TString    generatorName                = "DPMJET",         // generator name for weighting
                            TString    cutnumberAODBranch           = "800000006008400000001500000",   // cutnumber for AOD branch
                            Bool_t     isUsingTHnSparse             = kFALSE,           // enable or disable usage of THnSparses for background estimation
                            Int_t      enableExtMatchAndQA          = 0,                // disabled (0), extMatch (1), extQA_noCellQA (2), extMatch+extQA_noCellQA (3), extQA+cellQA (4), extMatch+extQA+cellQA (5)
                            Bool_t     enableTriggerMimicking       = kFALSE,           // enable trigger mimicking
                            Bool_t     enableTriggerOverlapRej      = kTRUE,            // enable trigger overlap rejection
                            Float_t    maxFacPtHard                 = 3,                // maximum factor between hardest jet and ptHard generated
                            TString    periodNameV0Reader           = "",               // period Name for V0Reader
                            Bool_t    doMultiplicityWeighting       = kFALSE,           // enable multiplicity weights
                            TString   fileNameInputForMultWeighing  = "",               // file for multiplicity weights
                            TString   periodNameAnchor              = "",               // name of anchor period for weighting
                            Bool_t     enableSortingMCLabels        = kTRUE,            // enable sorting for MC cluster labels
                            Bool_t     runLightOutput               = kFALSE,           // switch to run light output (only essential histograms for afterburner)
                            TString    additionalTrainConfig        = "0"               // additional counter for trainconfig
                           ) {

  Bool_t doTreeEOverP = kFALSE; // switch to produce EOverP tree
  TH1S* histoAcc = 0x0;         // histo for modified acceptance
  TString corrTaskSetting = ""; // select which correction task setting to use
  //parse additionalTrainConfig flag
  TObjArray *rAddConfigArr = additionalTrainConfig.Tokenize("_");
  if(rAddConfigArr->GetEntries()<1){cout << "ERROR: AddTask_GammaCalo_pPb during parsing of additionalTrainConfig String '" << additionalTrainConfig.Data() << "'" << endl; return;}
  TObjString* rAdditionalTrainConfig;
  for(Int_t i = 0; i<rAddConfigArr->GetEntries() ; i++){
    if(i==0) rAdditionalTrainConfig = (TObjString*)rAddConfigArr->At(i);
    else{
      TObjString* temp = (TObjString*) rAddConfigArr->At(i);
      TString tempStr = temp->GetString();
      if(tempStr.CompareTo("EPCLUSTree") == 0){
        cout << "INFO: AddTask_GammaCalo_pPb activating 'EPCLUSTree'" << endl;
        doTreeEOverP = kTRUE;
      }else if(tempStr.BeginsWith("MODIFYACC")){
        cout << "INFO: AddTask_GammaCalo_pPb activating 'MODIFYACC'" << endl;
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
        cout << "INFO: AddTask_GammaCalo_pPb will use custom branch from Correction Framework!" << endl;
        corrTaskSetting = tempStr;
        corrTaskSetting.Replace(0,2,"");
      }
    }
  }
  TString sAdditionalTrainConfig = rAdditionalTrainConfig->GetString();
  if (sAdditionalTrainConfig.Atoi() > 0){
    trainConfig = trainConfig + sAdditionalTrainConfig.Atoi();
    cout << "INFO: AddTask_GammaCalo_pPb running additionalTrainConfig '" << sAdditionalTrainConfig.Atoi() << "', train config: '" << trainConfig << "'" << endl;
  }
  if(corrTaskSetting.CompareTo(""))
    cout << "corrTaskSetting: " << corrTaskSetting.Data() << endl;

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
  TString V0ReaderName = Form("V0ReaderV1_%s_%s",cutnumberEvent.Data(),cutnumberPhoton.Data());
  if( !(AliV0ReaderV1*)mgr->GetTask(V0ReaderName.Data()) ){
    AliV0ReaderV1 *fV0ReaderV1 = new AliV0ReaderV1(V0ReaderName.Data());
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
      fCuts = new AliConversionPhotonCuts(cutnumberPhoton.Data(),cutnumberPhoton.Data());
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
  task->SetCorrectionTaskSetting(corrTaskSetting);
  task->SetLightOutput(runLightOutput);

  //create cut handler
  CutHandlerCalo cuts;

  // cluster cuts
  // 0 "ClusterType",  1 "EtaMin", 2 "EtaMax", 3 "PhiMin", 4 "PhiMax", 5 "DistanceToBadChannel", 6 "Timing", 7 "TrackMatching", 8 "ExoticCell",
  // 9 "MinEnergy", 10 "MinNCells", 11 "MinM02", 12 "MaxM02", 13 "MinMaxM20", 14 "RecConv", 15 "MaximumDispersion", 16 "NLM"

  // ===============================================================================================
  // Run 1 data EMC clusters pPb 5TeV
  // ===============================================================================================
  if (trainConfig == 1){ // no non lin
    cuts.AddCut("80010113","1111100057032230000","01631031000000d0");
    cuts.AddCut("80052113","1111100057032230000","01631031000000d0");
    cuts.AddCut("80085113","1111100057032230000","01631031000000d0");
    cuts.AddCut("80083113","1111100057032230000","01631031000000d0");
  } else if (trainConfig == 2){ // no non lin
    cuts.AddCut("80010113","1111100057032230000","01631031000000d0");
  } else if (trainConfig == 3){ // no non lin
    cuts.AddCut("80210113","1111100057032230000","01631031000000d0"); // 0-20
    cuts.AddCut("82410113","1111100057032230000","01631031000000d0"); // 20-40
    cuts.AddCut("84610113","1111100057032230000","01631031000000d0"); // 40-60
    cuts.AddCut("86010113","1111100057032230000","01631031000000d0"); // 60-80
  } else if (trainConfig == 4){ // no non lin, no time cut
    cuts.AddCut("80010113","1111100007032230000","01631031000000d0");
    cuts.AddCut("80052113","1111100007032230000","01631031000000d0");
    cuts.AddCut("80085113","1111100007032230000","01631031000000d0");
    cuts.AddCut("80083113","1111100007032230000","01631031000000d0");
  } else if (trainConfig == 5){ // no non lin, no time
    cuts.AddCut("80010113","1111100007032230000","01631031000000d0");
  } else if (trainConfig == 6){ // no non lin, no time cut, cent dep
    cuts.AddCut("80210113","1111100007032230000","01631031000000d0"); // 0-20
    cuts.AddCut("82410113","1111100007032230000","01631031000000d0"); // 20-40
    cuts.AddCut("84610113","1111100007032230000","01631031000000d0"); // 40-60
    cuts.AddCut("86010113","1111100007032230000","01631031000000d0"); // 60-80

  //-----------------------------------------------------------------------------------------------
  // Standard cuts
  //-----------------------------------------------------------------------------------------------
  } else if (trainConfig == 10){ // default cutstring smaller openangle 1 cell
    cuts.AddCut("80010113","1111141057032230000","01631031000000d0"); // default with 1 cell
    cuts.AddCut("80010113","1111141057032230000","0163103100000060"); // default with 17mrad
  } else if (trainConfig == 11){ // new default cut
    cuts.AddCut("80010113","1111141057032230000","01631031000000d0"); // default with 17mrad
  } else if (trainConfig == 12){ //all default triggers
    cuts.AddCut("80010113","1111141057032230000","01631031000000d0"); // default MB
    cuts.AddCut("80052113","1111141057032230000","01631031000000d0"); // default EMC7
    cuts.AddCut("80083113","1111141057032230000","01631031000000d0"); // default EG1
    cuts.AddCut("80085113","1111141057032230000","01631031000000d0"); // default EG2
  } else if (trainConfig == 13){ // testing past future protection
    cuts.AddCut("80010113","11111410570322l0000","01631031000000d0"); // default MB
    cuts.AddCut("80052113","11111410570322l0000","01631031000000d0"); // default EMC7
    cuts.AddCut("80083113","11111410570322l0000","01631031000000d0"); // default EG1
    cuts.AddCut("80085113","11111410570322l0000","01631031000000d0"); // default EG2
  //-----------------------------------------------------------------------------------------------
  // Standard cuts cent dependent
  //-----------------------------------------------------------------------------------------------
  } else if (trainConfig == 14){ // testing past future protection, timing 1000ns
    cuts.AddCut("80210113","1111141057032230000","01631031000000d0"); // 0-20
    cuts.AddCut("80210113","11111410570322l0000","01631031000000d0"); // 0-20
    cuts.AddCut("80252113","1111141057032230000","01631031000000d0"); // 0-20 EMC7
    cuts.AddCut("80285113","1111141057032230000","01631031000000d0"); // 0-20 EG2
    cuts.AddCut("80283113","1111141057032230000","01631031000000d0"); // 0-20 EG1
  } else if (trainConfig == 15){ // testing past future protection timing 100ns
    cuts.AddCut("82410113","1111141057032230000","01631031000000d0"); // 20-40
    cuts.AddCut("82410113","11111410570322l0000","01631031000000d0"); // 20-40
    cuts.AddCut("82452113","1111141057032230000","01631031000000d0"); // 20-40 EMC7
    cuts.AddCut("82485113","1111141057032230000","01631031000000d0"); // 20-40 EG2
    cuts.AddCut("82483113","1111141057032230000","01631031000000d0"); // 20-40 EG1
  } else if (trainConfig == 16){ // testing past future protection timing 100ns
    cuts.AddCut("84610113","1111141057032230000","01631031000000d0"); // 40-60
    cuts.AddCut("84610113","11111410570322l0000","01631031000000d0"); // 40-60
    cuts.AddCut("84652113","1111141057032230000","01631031000000d0"); // 40-60 EMC7
    cuts.AddCut("84685113","1111141057032230000","01631031000000d0"); // 40-60 EG2
    cuts.AddCut("84683113","1111141057032230000","01631031000000d0"); // 40-60 EG1
  } else if (trainConfig == 17){ // testing past future protection timing 100ns
    cuts.AddCut("86010113","1111141057032230000","01631031000000d0"); // 60-80
    cuts.AddCut("86010113","11111410570322l0000","01631031000000d0"); // 60-80
    cuts.AddCut("86052113","1111141057032230000","01631031000000d0"); // 60-80 EMC7
    cuts.AddCut("86085113","1111141057032230000","01631031000000d0"); // 60-80 EG2
    cuts.AddCut("86083113","1111141057032230000","01631031000000d0"); // 60-80 EG1


  } else if (trainConfig == 21){ // default cutstring, 1cell distance lead cell
    cuts.AddCut("80010113","1111141057032230000","01631031000000a0"); // 1 cell lead cell
    cuts.AddCut("80010113","1111141057032230000","01631031000000d0"); // 1 cell lead cell, 17mrad open
    cuts.AddCut("80010113","1111141057032230000","01631031000000b0"); // 1 cell lead cell, 15mrad open
  } else if (trainConfig == 22){ // default cutstring, M02 variations
    cuts.AddCut("80010113","1111141057032250000","01631031000000d0"); // 0.3
    cuts.AddCut("80010113","1111141057032260000","01631031000000d0"); // 0.27
    cuts.AddCut("80010113","1111141057032240000","01631031000000d0"); // 0.4
  } else if (trainConfig == 23){ // default cutstring, M02 variations
    cuts.AddCut("80010113","1111141057032290000","01631031000000d0"); // 0.35
    cuts.AddCut("80010113","11111410570322a0000","01631031000000d0"); // 0.33
    cuts.AddCut("80010113","11111410570322b0000","01631031000000d0"); // 0.28
    cuts.AddCut("80010113","11111410570322c0000","01631031000000d0"); // 0.32
  } else if (trainConfig == 24){ // default cutstring, cluster energy variations, decreased tender thresholds
    cuts.AddCut("80010113","1111141057012230000","01631031000000d0"); // E cluster > 0.5
    cuts.AddCut("80010113","1111141057032230000","01631031000000d0"); // E cluster > 0.7
    cuts.AddCut("80010113","1111141057052230000","01631031000000d0"); // E cluster > 0.9
    cuts.AddCut("80010113","11111410570b2230000","01631031000000d0"); // E cluster > 1.0
    cuts.AddCut("80010113","11111410570a2230000","01631031000000d0"); // E cluster > 1.5
  } else if (trainConfig == 25){ // default cutstring, cluster energy variations, same tender thresholds
    cuts.AddCut("80010113","1111141057032230000","01631031000000d0"); // E cluster > 0.7
    cuts.AddCut("80010113","1111141057052230000","01631031000000d0"); // E cluster > 0.9
    cuts.AddCut("80010113","11111410570b2230000","01631031000000d0"); // E cluster > 1.0
    cuts.AddCut("80010113","11111410570a2230000","01631031000000d0"); // E cluster > 1.5
  } else if (trainConfig == 26){ // default cutstring, cluster energy variations, increased tender thresholds
    cuts.AddCut("80010113","1111141057052230000","01631031000000d0"); // E cluster > 0.9
    cuts.AddCut("80010113","11111410570b2230000","01631031000000d0"); // E cluster > 1.0
    cuts.AddCut("80010113","11111410570a2230000","01631031000000d0"); // E cluster > 1.5
  //-----------------------------------------------------------------------------------------------
  // Systematics variations MB
  //-----------------------------------------------------------------------------------------------
  } else if (trainConfig == 30){ // nonlinearity variations
    cuts.AddCut("80010113","1111142057032230000","01631031000000d0"); // CRF
    cuts.AddCut("80010113","1111151057032230000","01631031000000d0"); // CCMF
    cuts.AddCut("80010113","1111152057032230000","01631031000000d0"); // CMF
  } else if (trainConfig == 31){ // second set of variations CLUSTER
    cuts.AddCut("80010113","1111141057022230000","01631031000000d0"); // min energy cluster variation 1  600 MeV
    cuts.AddCut("80010113","1111141057042230000","01631031000000d0"); // min energy cluster variation 2  800 MeV
    cuts.AddCut("80010113","1111141057052230000","01631031000000d0"); // min energy cluster variation 3  900 MeV
    cuts.AddCut("80010113","1111141057032230000","01631031000000d0"); // min/max M02  0.1<M<0.5
    cuts.AddCut("80010113","1111141057032200000","01631031000000d0"); // min/max M02  0.1<M<100
  } else if (trainConfig == 32){ // third set of variations CLUSTER
    cuts.AddCut("80010113","1111141057032250000","01631031000000d0"); // min/max M02  0.1<M<0.3
    cuts.AddCut("80010113","1111141057032260000","01631031000000d0"); // min/max M02  0.1<M<0.27
    cuts.AddCut("80010113","1111141057031230000","01631031000000d0"); // min number of cells variation 1  1 cell
    cuts.AddCut("80010113","1112141057032230000","01631031000000d0"); // only modules with TRD infront
    cuts.AddCut("80010113","1111341057032230000","01631031000000d0"); // no modules with TRD infront
  } else if (trainConfig == 33){ // third set of variations MESON
    cuts.AddCut("80010113","1111141057032230000","01633031000000d0"); // rapidity variation  y<0.6
    cuts.AddCut("80010113","1111141057032230000","01634031000000d0"); // rapidity variation  y<0.5
    cuts.AddCut("80010113","1111141057032230000","01631061000000d0"); // alpha meson variation 1   0<alpha<0.8
    cuts.AddCut("80010113","1111141057032230000","01631051000000d0"); // alpha meson variation 2  0<alpha<0.75
  } else if (trainConfig == 34){ // opening angle variations
    cuts.AddCut("80010113","1111141057032230000","0163103100000040"); // min opening angle 0.0152
    cuts.AddCut("80010113","1111141057032230000","0163103100000060"); // min opening angle 0.017
    cuts.AddCut("80010113","1111141057032230000","0163103100000070"); // min opening angle 0.016
    cuts.AddCut("80010113","1111141057032230000","0163103100000080"); // min opening angle 0.018
    cuts.AddCut("80010113","1111141057032230000","0163103100000090"); // min opening angle 0.018
  } else if (trainConfig == 35){ // TM variations
    cuts.AddCut("80010113","1111141053032230000","01631031000000d0"); // fixed window
    cuts.AddCut("80010113","1111141056032230000","01631031000000d0"); // tm pt dependent var 1
    cuts.AddCut("80010113","1111141058032230000","01631031000000d0"); // tm pt dependent var 2
    cuts.AddCut("80010113","1111141059032230000","01631031000000d0"); // tm pt dependent var 3
  } else if (trainConfig == 36){ // TM variations
    cuts.AddCut("80010113","1111153057032230000","01631031000000d0"); // CCMF nonlin+testbeam
    cuts.AddCut("80010113","1111154057032230000","01631031000000d0"); // CMF nonlin+testbeam
    cuts.AddCut("80010113","1111143057032230000","01631031000000d0"); // CCRF testbeam nonlin
    cuts.AddCut("80010113","1111144057032230000","01631031000000d0"); // CRF testbeam nonlin
    cuts.AddCut("80010113","1111102057032230000","01631031000000d0"); // testbeam nonlin

  } else if (trainConfig == 37){
    cuts.AddCut("80010113","11111410570322l0000","01631031000000d0"); // M02 pt dep with new std: 0.32, 0.0072, 0.5
    cuts.AddCut("80010113","11111410570322d0000","01631031000000d0"); // M02, pt dep with  0.27, 0.0072, 0.4
    cuts.AddCut("80010113","11111410570322e0000","01631031000000d0"); // M02, pt dep with  0.31, 0.0072, 0.5
    cuts.AddCut("80010113","11111410570322f0000","01631031000000d0"); // M02, pt dep with  0.36, 0.0072, 0.7
  } else if (trainConfig == 38){
    cuts.AddCut("80010113","11111410570322m0000","01631031000000d0"); // M02, pt dep with  0.32, 0.0152, 0.5
    cuts.AddCut("80010113","11111410570322g0000","01631031000000d0"); // M02, pt dep with  0.37, 0.0072, 0.7
    cuts.AddCut("80010113","11111410570322h0000","01631031000000d0"); // M02, pT-dep with  0.30, 0.0072, 0.5
    cuts.AddCut("80010113","11111410570322i0000","01631031000000d0"); // M02, pT-dep with  0.35, 0.0072, 0.7
  } else if (trainConfig == 39){
    cuts.AddCut("80010113","11111410570322j0000","01631031000000d0"); // M02, pT-dep with  0.25, 0.0072, 0.39
    cuts.AddCut("80010113","11111410570322r0000","01631031000000d0"); // M02, pT-dep with  0.25, 0.0072, 0.5
    cuts.AddCut("80010113","11111410570322s0000","01631031000000d0"); // M02, pT-dep with  0.32, 0.0238, 0.7
    cuts.AddCut("80010113","11111410570322n0000","01631031000000d0"); // M02, pT-dep with  0.32, 0.0238, 0.7

  //-----------------------------------------------------------------------------------------------
  // Systematics variations 0-20
  //-----------------------------------------------------------------------------------------------
  } else if (trainConfig == 40){ // nonlinearity variations
    cuts.AddCut("80210113","1111142057032230000","01631031000000d0"); // CRF
    cuts.AddCut("80210113","1111151057032230000","01631031000000d0"); // CCMF
    cuts.AddCut("80210113","1111152057032230000","01631031000000d0"); // CMF
  } else if (trainConfig == 41){ // second set of variations CLUSTER
    cuts.AddCut("80210113","1111141057022230000","01631031000000d0"); // min energy cluster variation 1  600 MeV
    cuts.AddCut("80210113","1111141057042230000","01631031000000d0"); // min energy cluster variation 2  800 MeV
    cuts.AddCut("80210113","1111141057052230000","01631031000000d0"); // min energy cluster variation 3  900 MeV
    cuts.AddCut("80210113","1111141057032230000","01631031000000d0"); // min/max M02  0.1<M<0.5
    cuts.AddCut("80210113","1111141057032200000","01631031000000d0"); // min/max M02  0.1<M<100
  } else if (trainConfig == 42){ // third set of variations CLUSTER
    cuts.AddCut("80210113","1111141057032250000","01631031000000d0"); // min/max M02  0.1<M<0.3
    cuts.AddCut("80210113","1111141057032260000","01631031000000d0"); // min/max M02  0.1<M<0.27
    cuts.AddCut("80210113","1111141057031230000","01631031000000d0"); // min number of cells variation 1  1 cell
    cuts.AddCut("80210113","1112141057032230000","01631031000000d0"); // only modules with TRD infront
    cuts.AddCut("80210113","1111341057032230000","01631031000000d0"); // no modules with TRD infront
  } else if (trainConfig == 43){ // third set of variations MESON
    cuts.AddCut("80210113","1111141057032230000","01633031000000d0"); // rapidity variation  y<0.6
    cuts.AddCut("80210113","1111141057032230000","01634031000000d0"); // rapidity variation  y<0.5
    cuts.AddCut("80210113","1111141057032230000","01631061000000d0"); // alpha meson variation 1   0<alpha<0.8
    cuts.AddCut("80210113","1111141057032230000","01631051000000d0"); // alpha meson variation 2  0<alpha<0.75
  } else if (trainConfig == 44){ // opening angle variations
    cuts.AddCut("80210113","1111141057032230000","0163103100000040"); // min opening angle 0.0152
    cuts.AddCut("80210113","1111141057032230000","0163103100000060"); // min opening angle 0.017
    cuts.AddCut("80210113","1111141057032230000","0163103100000070"); // min opening angle 0.016
    cuts.AddCut("80210113","1111141057032230000","0163103100000080"); // min opening angle 0.018
    cuts.AddCut("80210113","1111141057032230000","0163103100000090"); // min opening angle 0.018
  } else if (trainConfig == 45){ // TM variations
    cuts.AddCut("80210113","1111141053032230000","01631031000000d0"); // fixed window
    cuts.AddCut("80210113","1111141056032230000","01631031000000d0"); // tm pt dependent var 1
    cuts.AddCut("80210113","1111141058032230000","01631031000000d0"); // tm pt dependent var 2
    cuts.AddCut("80210113","1111141059032230000","01631031000000d0"); // tm pt dependent var 3
  } else if (trainConfig == 46){ // TM variations
    cuts.AddCut("80210113","1111153057032230000","01631031000000d0"); // CCMF nonlin+testbeam
    cuts.AddCut("80210113","1111154057032230000","01631031000000d0"); // CMF nonlin+testbeam
    cuts.AddCut("80210113","1111143057032230000","01631031000000d0"); // CCRF testbeam nonlin
    cuts.AddCut("80210113","1111144057032230000","01631031000000d0"); // CRF testbeam nonlin
    cuts.AddCut("80210113","1111102057032230000","01631031000000d0"); // testbeam nonlin
  } else if (trainConfig == 47){ //all default triggers
    cuts.AddCut("80210113","1111141057032230000","01631031000000d0"); // default MB
    cuts.AddCut("80252113","1111141057032230000","01631031000000d0"); // default EMC7
    cuts.AddCut("80283113","1111141057032230000","01631031000000d0"); // default EG1
    cuts.AddCut("80285113","1111141057032230000","01631031000000d0"); // default EG2

  //-----------------------------------------------------------------------------------------------
  // Systematics variations 20-40
  //-----------------------------------------------------------------------------------------------
  } else if (trainConfig == 50){ // nonlinearity variations
    cuts.AddCut("82410113","1111142057032230000","01631031000000d0"); // CRF
    cuts.AddCut("82410113","1111151057032230000","01631031000000d0"); // CCMF
    cuts.AddCut("82410113","1111152057032230000","01631031000000d0"); // CMF
  } else if (trainConfig == 51){ // second set of variations CLUSTER
    cuts.AddCut("82410113","1111141057022230000","01631031000000d0"); // min energy cluster variation 1  600 MeV
    cuts.AddCut("82410113","1111141057042230000","01631031000000d0"); // min energy cluster variation 2  800 MeV
    cuts.AddCut("82410113","1111141057052230000","01631031000000d0"); // min energy cluster variation 3  900 MeV
    cuts.AddCut("82410113","1111141057032230000","01631031000000d0"); // min/max M02  0.1<M<0.5
    cuts.AddCut("82410113","1111141057032200000","01631031000000d0"); // min/max M02  0.1<M<100
  } else if (trainConfig == 52){ // third set of variations CLUSTER
    cuts.AddCut("82410113","1111141057032250000","01631031000000d0"); // min/max M02  0.1<M<0.3
    cuts.AddCut("82410113","1111141057032260000","01631031000000d0"); // min/max M02  0.1<M<0.27
    cuts.AddCut("82410113","1111141057031230000","01631031000000d0"); // min number of cells variation 1  1 cell
    cuts.AddCut("82410113","1112141057032230000","01631031000000d0"); // only modules with TRD infront
    cuts.AddCut("82410113","1111341057032230000","01631031000000d0"); // no modules with TRD infront
  } else if (trainConfig == 53){ // third set of variations MESON
    cuts.AddCut("82410113","1111141057032230000","01633031000000d0"); // rapidity variation  y<0.6
    cuts.AddCut("82410113","1111141057032230000","01634031000000d0"); // rapidity variation  y<0.5
    cuts.AddCut("82410113","1111141057032230000","01631061000000d0"); // alpha meson variation 1   0<alpha<0.8
    cuts.AddCut("82410113","1111141057032230000","01631051000000d0"); // alpha meson variation 2  0<alpha<0.75
  } else if (trainConfig == 54){ // opening angle variations
    cuts.AddCut("82410113","1111141057032230000","0163103100000040"); // min opening angle 0.0152
    cuts.AddCut("82410113","1111141057032230000","0163103100000060"); // min opening angle 0.017
    cuts.AddCut("82410113","1111141057032230000","0163103100000070"); // min opening angle 0.016
    cuts.AddCut("82410113","1111141057032230000","0163103100000080"); // min opening angle 0.018
    cuts.AddCut("82410113","1111141057032230000","0163103100000090"); // min opening angle 0.018
  } else if (trainConfig == 55){ // TM variations
    cuts.AddCut("82410113","1111141053032230000","01631031000000d0"); // fixed window
    cuts.AddCut("82410113","1111141056032230000","01631031000000d0"); // tm pt dependent var 1
    cuts.AddCut("82410113","1111141058032230000","01631031000000d0"); // tm pt dependent var 2
    cuts.AddCut("82410113","1111141059032230000","01631031000000d0"); // tm pt dependent var 3
  } else if (trainConfig == 56){ // TM variations
    cuts.AddCut("82410113","1111153057032230000","01631031000000d0"); // CCMF nonlin+testbeam
    cuts.AddCut("82410113","1111154057032230000","01631031000000d0"); // CMF nonlin+testbeam
    cuts.AddCut("82410113","1111143057032230000","01631031000000d0"); // CCRF testbeam nonlin
    cuts.AddCut("82410113","1111144057032230000","01631031000000d0"); // CRF testbeam nonlin
    cuts.AddCut("82410113","1111102057032230000","01631031000000d0"); // testbeam nonlin
  } else if (trainConfig == 57){ //all default triggers
    cuts.AddCut("82410113","1111141057032230000","01631031000000d0"); // default MB
    cuts.AddCut("82452113","1111141057032230000","01631031000000d0"); // default EMC7
    cuts.AddCut("82483113","1111141057032230000","01631031000000d0"); // default EG1
    cuts.AddCut("82485113","1111141057032230000","01631031000000d0"); // default EG2

  //-----------------------------------------------------------------------------------------------
  // Systematics variations 40-60
  //-----------------------------------------------------------------------------------------------
  } else if (trainConfig == 60){ // nonlinearity variations
    cuts.AddCut("84610113","1111142057032230000","01631031000000d0"); // CRF
    cuts.AddCut("84610113","1111151057032230000","01631031000000d0"); // CCMF
    cuts.AddCut("84610113","1111152057032230000","01631031000000d0"); // CMF
  } else if (trainConfig == 61){ // second set of variations CLUSTER
    cuts.AddCut("84610113","1111141057022230000","01631031000000d0"); // min energy cluster variation 1  600 MeV
    cuts.AddCut("84610113","1111141057042230000","01631031000000d0"); // min energy cluster variation 2  800 MeV
    cuts.AddCut("84610113","1111141057052230000","01631031000000d0"); // min energy cluster variation 3  900 MeV
    cuts.AddCut("84610113","1111141057032230000","01631031000000d0"); // min/max M02  0.1<M<0.5
    cuts.AddCut("84610113","1111141057032200000","01631031000000d0"); // min/max M02  0.1<M<100
  } else if (trainConfig == 62){ // third set of variations CLUSTER
    cuts.AddCut("84610113","1111141057032250000","01631031000000d0"); // min/max M02  0.1<M<0.3
    cuts.AddCut("84610113","1111141057032260000","01631031000000d0"); // min/max M02  0.1<M<0.27
    cuts.AddCut("84610113","1111141057031230000","01631031000000d0"); // min number of cells variation 1  1 cell
    cuts.AddCut("84610113","1112141057032230000","01631031000000d0"); // only modules with TRD infront
    cuts.AddCut("84610113","1111341057032230000","01631031000000d0"); // no modules with TRD infront
  } else if (trainConfig == 63){ // third set of variations MESON
    cuts.AddCut("84610113","1111141057032230000","01633031000000d0"); // rapidity variation  y<0.6
    cuts.AddCut("84610113","1111141057032230000","01634031000000d0"); // rapidity variation  y<0.5
    cuts.AddCut("84610113","1111141057032230000","01631061000000d0"); // alpha meson variation 1   0<alpha<0.8
    cuts.AddCut("84610113","1111141057032230000","01631051000000d0"); // alpha meson variation 2  0<alpha<0.75
  } else if (trainConfig == 64){ // opening angle variations
    cuts.AddCut("84610113","1111141057032230000","0163103100000040"); // min opening angle 0.0152
    cuts.AddCut("84610113","1111141057032230000","0163103100000060"); // min opening angle 0.017
    cuts.AddCut("84610113","1111141057032230000","0163103100000070"); // min opening angle 0.016
    cuts.AddCut("84610113","1111141057032230000","0163103100000080"); // min opening angle 0.018
    cuts.AddCut("84610113","1111141057032230000","0163103100000090"); // min opening angle 0.018
  } else if (trainConfig == 65){ // TM variations
    cuts.AddCut("84610113","1111141053032230000","01631031000000d0"); // fixed window
    cuts.AddCut("84610113","1111141056032230000","01631031000000d0"); // tm pt dependent var 1
    cuts.AddCut("84610113","1111141058032230000","01631031000000d0"); // tm pt dependent var 2
    cuts.AddCut("84610113","1111141059032230000","01631031000000d0"); // tm pt dependent var 3
  } else if (trainConfig == 66){ // TM variations
    cuts.AddCut("84610113","1111153057032230000","01631031000000d0"); // CCMF nonlin+testbeam
    cuts.AddCut("84610113","1111154057032230000","01631031000000d0"); // CMF nonlin+testbeam
    cuts.AddCut("84610113","1111143057032230000","01631031000000d0"); // CCRF testbeam nonlin
    cuts.AddCut("84610113","1111144057032230000","01631031000000d0"); // CRF testbeam nonlin
    cuts.AddCut("84610113","1111102057032230000","01631031000000d0"); // testbeam nonlin
  } else if (trainConfig == 67){ //all default triggers
    cuts.AddCut("84610113","1111141057032230000","01631031000000d0"); // default MB
    cuts.AddCut("84652113","1111141057032230000","01631031000000d0"); // default EMC7
    cuts.AddCut("84683113","1111141057032230000","01631031000000d0"); // default EG1
    cuts.AddCut("84685113","1111141057032230000","01631031000000d0"); // default EG2

  //-----------------------------------------------------------------------------------------------
  // Systematics variations 60-100
  //-----------------------------------------------------------------------------------------------
  } else if (trainConfig == 70){ // nonlinearity variations
    cuts.AddCut("86010113","1111142057032230000","01631031000000d0"); // CRF
    cuts.AddCut("86010113","1111151057032230000","01631031000000d0"); // CCMF
    cuts.AddCut("86010113","1111152057032230000","01631031000000d0"); // CMF
  } else if (trainConfig == 71){ // second set of variations CLUSTER
    cuts.AddCut("86010113","1111141057022230000","01631031000000d0"); // min energy cluster variation 1  600 MeV
    cuts.AddCut("86010113","1111141057042230000","01631031000000d0"); // min energy cluster variation 2  800 MeV
    cuts.AddCut("86010113","1111141057052230000","01631031000000d0"); // min energy cluster variation 3  900 MeV
    cuts.AddCut("86010113","1111141057032230000","01631031000000d0"); // min/max M02  0.1<M<0.5
    cuts.AddCut("86010113","1111141057032200000","01631031000000d0"); // min/max M02  0.1<M<100
  } else if (trainConfig == 72){ // third set of variations CLUSTER
    cuts.AddCut("86010113","1111141057032250000","01631031000000d0"); // min/max M02  0.1<M<0.3
    cuts.AddCut("86010113","1111141057032260000","01631031000000d0"); // min/max M02  0.1<M<0.27
    cuts.AddCut("86010113","1111141057031230000","01631031000000d0"); // min number of cells variation 1  1 cell
    cuts.AddCut("86010113","1112141057032230000","01631031000000d0"); // only modules with TRD infront
    cuts.AddCut("86010113","1111341057032230000","01631031000000d0"); // no modules with TRD infront
  } else if (trainConfig == 73){ // third set of variations MESON
    cuts.AddCut("86010113","1111141057032230000","01633031000000d0"); // rapidity variation  y<0.6
    cuts.AddCut("86010113","1111141057032230000","01634031000000d0"); // rapidity variation  y<0.5
    cuts.AddCut("86010113","1111141057032230000","01631061000000d0"); // alpha meson variation 1   0<alpha<0.8
    cuts.AddCut("86010113","1111141057032230000","01631051000000d0"); // alpha meson variation 2  0<alpha<0.75
  } else if (trainConfig == 74){ // opening angle variations
    cuts.AddCut("86010113","1111141057032230000","0163103100000040"); // min opening angle 0.0152
    cuts.AddCut("86010113","1111141057032230000","0163103100000060"); // min opening angle 0.017
    cuts.AddCut("86010113","1111141057032230000","0163103100000070"); // min opening angle 0.016
    cuts.AddCut("86010113","1111141057032230000","0163103100000080"); // min opening angle 0.018
    cuts.AddCut("86010113","1111141057032230000","0163103100000090"); // min opening angle 0.018
  } else if (trainConfig == 75){ // TM variations
    cuts.AddCut("86010113","1111141053032230000","01631031000000d0"); // fixed window
    cuts.AddCut("86010113","1111141056032230000","01631031000000d0"); // tm pt dependent var 1
    cuts.AddCut("86010113","1111141058032230000","01631031000000d0"); // tm pt dependent var 2
    cuts.AddCut("86010113","1111141059032230000","01631031000000d0"); // tm pt dependent var 3
  } else if (trainConfig == 76){ // TM variations
    cuts.AddCut("86010113","1111153057032230000","01631031000000d0"); // CCMF nonlin+testbeam
    cuts.AddCut("86010113","1111154057032230000","01631031000000d0"); // CMF nonlin+testbeam
    cuts.AddCut("86010113","1111143057032230000","01631031000000d0"); // CCRF testbeam nonlin
    cuts.AddCut("86010113","1111144057032230000","01631031000000d0"); // CRF testbeam nonlin
    cuts.AddCut("86010113","1111102057032230000","01631031000000d0"); // testbeam nonlin
  } else if (trainConfig == 77){ //all default triggers
    cuts.AddCut("86010113","1111141057032230000","01631031000000d0"); // default MB
    cuts.AddCut("86052113","1111141057032230000","01631031000000d0"); // default EMC7
    cuts.AddCut("86083113","1111141057032230000","01631031000000d0"); // default EG1
    cuts.AddCut("86085113","1111141057032230000","01631031000000d0"); // default EG2

  //-----------------------------------------------------------------------------------------------
  //testing different event mix methods
  //-----------------------------------------------------------------------------------------------
  } else if (trainConfig == 81){
    cuts.AddCut("80010113","1111141057032230000","01631031000000d0"); // default (V0 mult)
    cuts.AddCut("80010113","1111141057032230000","02631031000000d0"); // using track mult
    cuts.AddCut("80010113","1111141057032230000","09631031000000d0"); // using PtMax method

  //-----------------------------------------------------------------------------------------------
  // Systematics variations 0-20
  //-----------------------------------------------------------------------------------------------
  } else if (trainConfig == 90){ // nonlinearity variations
    cuts.AddCut("80210113","11111420570322l0000","01631031000000d0"); // CRF
    cuts.AddCut("80210113","11111510570322l0000","01631031000000d0"); // CCMF
    cuts.AddCut("80210113","11111520570322l0000","01631031000000d0"); // CMF
  } else if (trainConfig == 91){ // second set of variations CLUSTER
    cuts.AddCut("80210113","11111410570222l0000","01631031000000d0"); // min energy cluster variation 1  600 MeV
    cuts.AddCut("80210113","11111410570422l0000","01631031000000d0"); // min energy cluster variation 2  800 MeV
    cuts.AddCut("80210113","11111410570522l0000","01631031000000d0"); // min energy cluster variation 3  900 MeV
    cuts.AddCut("80210113","11111410570322l0000","01631031000000d0"); // min/max M02  0.1<M<0.5
  } else if (trainConfig == 92){ // third set of variations CLUSTER
    cuts.AddCut("80210113","1111141057032250000","01631031000000d0"); // min/max M02  0.1<M<0.3
    cuts.AddCut("80210113","1111141057032260000","01631031000000d0"); // min/max M02  0.1<M<0.27
    cuts.AddCut("80210113","11111410570312l0000","01631031000000d0"); // min number of cells variation 1  1 cell
    cuts.AddCut("80210113","11121410570322l0000","01631031000000d0"); // only modules with TRD infront
    cuts.AddCut("80210113","11113410570322l0000","01631031000000d0"); // no modules with TRD infront
  } else if (trainConfig == 93){ // third set of variations MESON
    cuts.AddCut("80210113","11111410570322l0000","01633031000000d0"); // rapidity variation  y<0.6
    cuts.AddCut("80210113","11111410570322l0000","01634031000000d0"); // rapidity variation  y<0.5
    cuts.AddCut("80210113","11111410570322l0000","01631061000000d0"); // alpha meson variation 1   0<alpha<0.8
    cuts.AddCut("80210113","11111410570322l0000","01631051000000d0"); // alpha meson variation 2  0<alpha<0.75
  } else if (trainConfig == 94){ // opening angle variations
    cuts.AddCut("80210113","11111410570322l0000","0163103100000040"); // min opening angle 0.0152
    cuts.AddCut("80210113","11111410570322l0000","0163103100000060"); // min opening angle 0.017
    cuts.AddCut("80210113","11111410570322l0000","0163103100000070"); // min opening angle 0.016
    cuts.AddCut("80210113","11111410570322l0000","0163103100000080"); // min opening angle 0.018
    cuts.AddCut("80210113","11111410570322l0000","0163103100000090"); // min opening angle 0.018
  } else if (trainConfig == 95){ // TM variations
    cuts.AddCut("80210113","11111410530322l0000","01631031000000d0"); // fixed window
    cuts.AddCut("80210113","11111410560322l0000","01631031000000d0"); // tm pt dependent var 1
    cuts.AddCut("80210113","11111410580322l0000","01631031000000d0"); // tm pt dependent var 2
    cuts.AddCut("80210113","11111410590322l0000","01631031000000d0"); // tm pt dependent var 3
  } else if (trainConfig == 96){ // TM variations
    cuts.AddCut("80210113","11111530570322l0000","01631031000000d0"); // CCMF nonlin+testbeam
    cuts.AddCut("80210113","11111540570322l0000","01631031000000d0"); // CMF nonlin+testbeam
    cuts.AddCut("80210113","11111430570322l0000","01631031000000d0"); // CCRF testbeam nonlin
    cuts.AddCut("80210113","11111440570322l0000","01631031000000d0"); // CRF testbeam nonlin
    cuts.AddCut("80210113","11111020570322l0000","01631031000000d0"); // testbeam nonlin
  } else if (trainConfig == 97){
    cuts.AddCut("80210113","1111141057032230000","01631031000000d0"); // M02 pt dep with new std: 0.32, 0.0072, 0.5
    cuts.AddCut("80210113","11111410570322d0000","01631031000000d0"); // M02, pt dep with  0.27, 0.0072, 0.4
    cuts.AddCut("80210113","11111410570322e0000","01631031000000d0"); // M02, pt dep with  0.31, 0.0072, 0.5
    cuts.AddCut("80210113","11111410570322f0000","01631031000000d0"); // M02, pt dep with  0.36, 0.0072, 0.7
  } else if (trainConfig == 98){
    cuts.AddCut("80210113","11111410570322m0000","01631031000000d0"); // M02, pt dep with  0.32, 0.0152, 0.5
    cuts.AddCut("80210113","11111410570322g0000","01631031000000d0"); // M02, pt dep with  0.37, 0.0072, 0.7
    cuts.AddCut("80210113","11111410570322h0000","01631031000000d0"); // M02, pT-dep with  0.30, 0.0072, 0.5
    cuts.AddCut("80210113","11111410570322i0000","01631031000000d0"); // M02, pT-dep with  0.35, 0.0072, 0.7
  } else if (trainConfig == 99){
    cuts.AddCut("80210113","11111410570322j0000","01631031000000d0"); // M02, pT-dep with  0.25, 0.0072, 0.39
    cuts.AddCut("80210113","11111410570322r0000","01631031000000d0"); // M02, pT-dep with  0.25, 0.0072, 0.5
    cuts.AddCut("80210113","11111410570322s0000","01631031000000d0"); // M02, pT-dep with  0.32, 0.0238, 0.7
    cuts.AddCut("80210113","11111410570322n0000","01631031000000d0"); // M02, pT-dep with  0.32, 0.0238, 0.7

  //-----------------------------------------------------------------------------------------------
  // Systematics variations 20-40
  //-----------------------------------------------------------------------------------------------
  } else if (trainConfig == 100){ // nonlinearity variations
    cuts.AddCut("82410113","11111420570322l0000","01631031000000d0"); // CRF
    cuts.AddCut("82410113","11111510570322l0000","01631031000000d0"); // CCMF
    cuts.AddCut("82410113","11111520570322l0000","01631031000000d0"); // CMF
  } else if (trainConfig == 101){ // second set of variations CLUSTER
    cuts.AddCut("82410113","11111410570222l0000","01631031000000d0"); // min energy cluster variation 1  600 MeV
    cuts.AddCut("82410113","11111410570422l0000","01631031000000d0"); // min energy cluster variation 2  800 MeV
    cuts.AddCut("82410113","11111410570522l0000","01631031000000d0"); // min energy cluster variation 3  900 MeV
    cuts.AddCut("82410113","11111410570322l0000","01631031000000d0"); // min/max M02  0.1<M<0.5
  } else if (trainConfig == 102){ // third set of variations CLUSTER
    cuts.AddCut("82410113","1111141057032250000","01631031000000d0"); // min/max M02  0.1<M<0.3
    cuts.AddCut("82410113","1111141057032260000","01631031000000d0"); // min/max M02  0.1<M<0.27
    cuts.AddCut("82410113","11111410570312l0000","01631031000000d0"); // min number of cells variation 1  1 cell
    cuts.AddCut("82410113","11121410570322l0000","01631031000000d0"); // only modules with TRD infront
    cuts.AddCut("82410113","11113410570322l0000","01631031000000d0"); // no modules with TRD infront
  } else if (trainConfig == 103){ // third set of variations MESON
    cuts.AddCut("82410113","11111410570322l0000","01633031000000d0"); // rapidity variation  y<0.6
    cuts.AddCut("82410113","11111410570322l0000","01634031000000d0"); // rapidity variation  y<0.5
    cuts.AddCut("82410113","11111410570322l0000","01631061000000d0"); // alpha meson variation 1   0<alpha<0.8
    cuts.AddCut("82410113","11111410570322l0000","01631051000000d0"); // alpha meson variation 2  0<alpha<0.75
  } else if (trainConfig == 104){ // opening angle variations
    cuts.AddCut("82410113","11111410570322l0000","0163103100000040"); // min opening angle 0.0152
    cuts.AddCut("82410113","11111410570322l0000","0163103100000060"); // min opening angle 0.017
    cuts.AddCut("82410113","11111410570322l0000","0163103100000070"); // min opening angle 0.016
    cuts.AddCut("82410113","11111410570322l0000","0163103100000080"); // min opening angle 0.018
    cuts.AddCut("82410113","11111410570322l0000","0163103100000090"); // min opening angle 0.018
  } else if (trainConfig == 105){ // TM variations
    cuts.AddCut("82410113","11111410530322l0000","01631031000000d0"); // fixed window
    cuts.AddCut("82410113","11111410560322l0000","01631031000000d0"); // tm pt dependent var 1
    cuts.AddCut("82410113","11111410580322l0000","01631031000000d0"); // tm pt dependent var 2
    cuts.AddCut("82410113","11111410590322l0000","01631031000000d0"); // tm pt dependent var 3
  } else if (trainConfig == 106){ // TM variations
    cuts.AddCut("82410113","11111530570322l0000","01631031000000d0"); // CCMF nonlin+testbeam
    cuts.AddCut("82410113","11111540570322l0000","01631031000000d0"); // CMF nonlin+testbeam
    cuts.AddCut("82410113","11111430570322l0000","01631031000000d0"); // CCRF testbeam nonlin
    cuts.AddCut("82410113","11111440570322l0000","01631031000000d0"); // CRF testbeam nonlin
    cuts.AddCut("82410113","11111020570322l0000","01631031000000d0"); // testbeam nonlin
  } else if (trainConfig == 107){
    cuts.AddCut("82410113","1111141057032230000","01631031000000d0"); // M02 pt dep with new std: 0.32, 0.0072, 0.5
    cuts.AddCut("82410113","11111410570322d0000","01631031000000d0"); // M02, pt dep with  0.27, 0.0072, 0.4
    cuts.AddCut("82410113","11111410570322e0000","01631031000000d0"); // M02, pt dep with  0.31, 0.0072, 0.5
    cuts.AddCut("82410113","11111410570322f0000","01631031000000d0"); // M02, pt dep with  0.36, 0.0072, 0.7
  } else if (trainConfig == 108){
    cuts.AddCut("82410113","11111410570322m0000","01631031000000d0"); // M02, pt dep with  0.32, 0.0152, 0.5
    cuts.AddCut("82410113","11111410570322g0000","01631031000000d0"); // M02, pt dep with  0.37, 0.0072, 0.7
    cuts.AddCut("82410113","11111410570322h0000","01631031000000d0"); // M02, pT-dep with  0.30, 0.0072, 0.5
    cuts.AddCut("82410113","11111410570322i0000","01631031000000d0"); // M02, pT-dep with  0.35, 0.0072, 0.7
  } else if (trainConfig == 109){
    cuts.AddCut("82410113","11111410570322j0000","01631031000000d0"); // M02, pT-dep with  0.25, 0.0072, 0.39
    cuts.AddCut("82410113","11111410570322r0000","01631031000000d0"); // M02, pT-dep with  0.25, 0.0072, 0.5
    cuts.AddCut("82410113","11111410570322s0000","01631031000000d0"); // M02, pT-dep with  0.32, 0.0238, 0.7
    cuts.AddCut("82410113","11111410570322n0000","01631031000000d0"); // M02, pT-dep with  0.32, 0.0238, 0.7

  //-----------------------------------------------------------------------------------------------
  // Systematics variations 40-60
  //-----------------------------------------------------------------------------------------------
  } else if (trainConfig == 110){ // nonlinearity variations
    cuts.AddCut("84610113","11111420570322l0000","01631031000000d0"); // CRF
    cuts.AddCut("84610113","11111510570322l0000","01631031000000d0"); // CCMF
    cuts.AddCut("84610113","11111520570322l0000","01631031000000d0"); // CMF
  } else if (trainConfig == 111){ // second set of variations CLUSTER
    cuts.AddCut("84610113","11111410570222l0000","01631031000000d0"); // min energy cluster variation 1  600 MeV
    cuts.AddCut("84610113","11111410570422l0000","01631031000000d0"); // min energy cluster variation 2  800 MeV
    cuts.AddCut("84610113","11111410570522l0000","01631031000000d0"); // min energy cluster variation 3  900 MeV
    cuts.AddCut("84610113","11111410570322l0000","01631031000000d0"); // min/max M02  0.1<M<0.5
  } else if (trainConfig == 112){ // third set of variations CLUSTER
    cuts.AddCut("84610113","1111141057032250000","01631031000000d0"); // min/max M02  0.1<M<0.3
    cuts.AddCut("84610113","1111141057032260000","01631031000000d0"); // min/max M02  0.1<M<0.27
    cuts.AddCut("84610113","11111410570312l0000","01631031000000d0"); // min number of cells variation 1  1 cell
    cuts.AddCut("84610113","11121410570322l0000","01631031000000d0"); // only modules with TRD infront
    cuts.AddCut("84610113","11113410570322l0000","01631031000000d0"); // no modules with TRD infront
  } else if (trainConfig == 113){ // third set of variations MESON
    cuts.AddCut("84610113","11111410570322l0000","01633031000000d0"); // rapidity variation  y<0.6
    cuts.AddCut("84610113","11111410570322l0000","01634031000000d0"); // rapidity variation  y<0.5
    cuts.AddCut("84610113","11111410570322l0000","01631061000000d0"); // alpha meson variation 1   0<alpha<0.8
    cuts.AddCut("84610113","11111410570322l0000","01631051000000d0"); // alpha meson variation 2  0<alpha<0.75
  } else if (trainConfig == 114){ // opening angle variations
    cuts.AddCut("84610113","11111410570322l0000","0163103100000040"); // min opening angle 0.0152
    cuts.AddCut("84610113","11111410570322l0000","0163103100000060"); // min opening angle 0.017
    cuts.AddCut("84610113","11111410570322l0000","0163103100000070"); // min opening angle 0.016
    cuts.AddCut("84610113","11111410570322l0000","0163103100000080"); // min opening angle 0.018
    cuts.AddCut("84610113","11111410570322l0000","0163103100000090"); // min opening angle 0.018
  } else if (trainConfig == 115){ // TM variations
    cuts.AddCut("84610113","11111410530322l0000","01631031000000d0"); // fixed window
    cuts.AddCut("84610113","11111410560322l0000","01631031000000d0"); // tm pt dependent var 1
    cuts.AddCut("84610113","11111410580322l0000","01631031000000d0"); // tm pt dependent var 2
    cuts.AddCut("84610113","11111410590322l0000","01631031000000d0"); // tm pt dependent var 3
  } else if (trainConfig == 116){ // TM variations
    cuts.AddCut("84610113","11111530570322l0000","01631031000000d0"); // CCMF nonlin+testbeam
    cuts.AddCut("84610113","11111540570322l0000","01631031000000d0"); // CMF nonlin+testbeam
    cuts.AddCut("84610113","11111430570322l0000","01631031000000d0"); // CCRF testbeam nonlin
    cuts.AddCut("84610113","11111440570322l0000","01631031000000d0"); // CRF testbeam nonlin
    cuts.AddCut("84610113","11111020570322l0000","01631031000000d0"); // testbeam nonlin
  } else if (trainConfig == 117){
    cuts.AddCut("84610113","1111141057032230000","01631031000000d0"); // M02 pt dep with new std: 0.32, 0.0072, 0.5
    cuts.AddCut("84610113","11111410570322d0000","01631031000000d0"); // M02, pt dep with  0.27, 0.0072, 0.4
    cuts.AddCut("84610113","11111410570322e0000","01631031000000d0"); // M02, pt dep with  0.31, 0.0072, 0.5
    cuts.AddCut("84610113","11111410570322f0000","01631031000000d0"); // M02, pt dep with  0.36, 0.0072, 0.7
  } else if (trainConfig == 118){
    cuts.AddCut("84610113","11111410570322m0000","01631031000000d0"); // M02, pt dep with  0.32, 0.0152, 0.5
    cuts.AddCut("84610113","11111410570322g0000","01631031000000d0"); // M02, pt dep with  0.37, 0.0072, 0.7
    cuts.AddCut("84610113","11111410570322h0000","01631031000000d0"); // M02, pT-dep with  0.30, 0.0072, 0.5
    cuts.AddCut("84610113","11111410570322i0000","01631031000000d0"); // M02, pT-dep with  0.35, 0.0072, 0.7
  } else if (trainConfig == 119){
    cuts.AddCut("84610113","11111410570322j0000","01631031000000d0"); // M02, pT-dep with  0.25, 0.0072, 0.39
    cuts.AddCut("84610113","11111410570322r0000","01631031000000d0"); // M02, pT-dep with  0.25, 0.0072, 0.5
    cuts.AddCut("84610113","11111410570322s0000","01631031000000d0"); // M02, pT-dep with  0.32, 0.0238, 0.7
    cuts.AddCut("84610113","11111410570322n0000","01631031000000d0"); // M02, pT-dep with  0.32, 0.0238, 0.7

  //-----------------------------------------------------------------------------------------------
  // Systematics variations 60-100
  //-----------------------------------------------------------------------------------------------
  } else if (trainConfig == 120){ // nonlinearity variations
    cuts.AddCut("86010113","11111420570322l0000","01631031000000d0"); // CRF
    cuts.AddCut("86010113","11111510570322l0000","01631031000000d0"); // CCMF
    cuts.AddCut("86010113","11111520570322l0000","01631031000000d0"); // CMF
  } else if (trainConfig == 121){ // second set of variations CLUSTER
    cuts.AddCut("86010113","11111410570222l0000","01631031000000d0"); // min energy cluster variation 1  600 MeV
    cuts.AddCut("86010113","11111410570422l0000","01631031000000d0"); // min energy cluster variation 2  800 MeV
    cuts.AddCut("86010113","11111410570522l0000","01631031000000d0"); // min energy cluster variation 3  900 MeV
    cuts.AddCut("86010113","11111410570322l0000","01631031000000d0"); // min/max M02  0.1<M<0.5
  } else if (trainConfig == 122){ // third set of variations CLUSTER
    cuts.AddCut("86010113","1111141057032250000","01631031000000d0"); // min/max M02  0.1<M<0.3
    cuts.AddCut("86010113","1111141057032260000","01631031000000d0"); // min/max M02  0.1<M<0.27
    cuts.AddCut("86010113","11111410570312l0000","01631031000000d0"); // min number of cells variation 1  1 cell
    cuts.AddCut("86010113","11121410570322l0000","01631031000000d0"); // only modules with TRD infront
    cuts.AddCut("86010113","11113410570322l0000","01631031000000d0"); // no modules with TRD infront
  } else if (trainConfig == 123){ // third set of variations MESON
    cuts.AddCut("86010113","11111410570322l0000","01633031000000d0"); // rapidity variation  y<0.6
    cuts.AddCut("86010113","11111410570322l0000","01634031000000d0"); // rapidity variation  y<0.5
    cuts.AddCut("86010113","11111410570322l0000","01631061000000d0"); // alpha meson variation 1   0<alpha<0.8
    cuts.AddCut("86010113","11111410570322l0000","01631051000000d0"); // alpha meson variation 2  0<alpha<0.75
  } else if (trainConfig == 124){ // opening angle variations
    cuts.AddCut("86010113","11111410570322l0000","0163103100000040"); // min opening angle 0.0152
    cuts.AddCut("86010113","11111410570322l0000","0163103100000060"); // min opening angle 0.017
    cuts.AddCut("86010113","11111410570322l0000","0163103100000070"); // min opening angle 0.016
    cuts.AddCut("86010113","11111410570322l0000","0163103100000080"); // min opening angle 0.018
    cuts.AddCut("86010113","11111410570322l0000","0163103100000090"); // min opening angle 0.018
  } else if (trainConfig == 125){ // TM variations
    cuts.AddCut("86010113","11111410530322l0000","01631031000000d0"); // fixed window
    cuts.AddCut("86010113","11111410560322l0000","01631031000000d0"); // tm pt dependent var 1
    cuts.AddCut("86010113","11111410580322l0000","01631031000000d0"); // tm pt dependent var 2
    cuts.AddCut("86010113","11111410590322l0000","01631031000000d0"); // tm pt dependent var 3
  } else if (trainConfig == 126){ // TM variations
    cuts.AddCut("86010113","11111530570322l0000","01631031000000d0"); // CCMF nonlin+testbeam
    cuts.AddCut("86010113","11111540570322l0000","01631031000000d0"); // CMF nonlin+testbeam
    cuts.AddCut("86010113","11111430570322l0000","01631031000000d0"); // CCRF testbeam nonlin
    cuts.AddCut("86010113","11111440570322l0000","01631031000000d0"); // CRF testbeam nonlin
    cuts.AddCut("86010113","11111020570322l0000","01631031000000d0"); // testbeam nonlin
  } else if (trainConfig == 127){
    cuts.AddCut("86010113","1111141057032230000","01631031000000d0"); // M02 pt dep with new std: 0.32, 0.0072, 0.5
    cuts.AddCut("86010113","11111410570322d0000","01631031000000d0"); // M02, pt dep with  0.27, 0.0072, 0.4
    cuts.AddCut("86010113","11111410570322e0000","01631031000000d0"); // M02, pt dep with  0.31, 0.0072, 0.5
    cuts.AddCut("86010113","11111410570322f0000","01631031000000d0"); // M02, pt dep with  0.36, 0.0072, 0.7
  } else if (trainConfig == 128){
    cuts.AddCut("86010113","11111410570322m0000","01631031000000d0"); // M02, pt dep with  0.32, 0.0152, 0.5
    cuts.AddCut("86010113","11111410570322g0000","01631031000000d0"); // M02, pt dep with  0.37, 0.0072, 0.7
    cuts.AddCut("86010113","11111410570322h0000","01631031000000d0"); // M02, pT-dep with  0.30, 0.0072, 0.5
    cuts.AddCut("86010113","11111410570322i0000","01631031000000d0"); // M02, pT-dep with  0.35, 0.0072, 0.7
  } else if (trainConfig == 129){
    cuts.AddCut("86010113","11111410570322j0000","01631031000000d0"); // M02, pT-dep with  0.25, 0.0072, 0.39
    cuts.AddCut("86010113","11111410570322r0000","01631031000000d0"); // M02, pT-dep with  0.25, 0.0072, 0.5
    cuts.AddCut("86010113","11111410570322s0000","01631031000000d0"); // M02, pT-dep with  0.32, 0.0238, 0.7
    cuts.AddCut("86010113","11111410570322n0000","01631031000000d0"); // M02, pT-dep with  0.32, 0.0238, 0.7

  // ===============================================================================================
  // Run 2 data EMC clusters pPb 5TeV
  // ===============================================================================================
  } else if (trainConfig == 200){ // EMCAL clusters standard cuts,
    cuts.AddCut("80010113","1111100057032230000","01631031000000d0"); // 0-100
  } else if (trainConfig == 201){ // EMCAL clusters standard cuts, no nonlin, open timing
    cuts.AddCut("80010113","1111100017032230000","01631031000000d0"); // INT7
  } else if (trainConfig == 202){ // EMCAL clusters standard cuts, no nonlin, +-50ns, trigger
    cuts.AddCut("80052113","1111100057032230000","01631031000000d0"); // EMC7
    cuts.AddCut("80085113","1111100057032230000","01631031000000d0"); // EG2
    cuts.AddCut("80083113","1111100057032230000","01631031000000d0"); // EG1
  } else if (trainConfig == 203){ // EMCAL clusters standard cuts, no nonlin, open timing, trigger
    cuts.AddCut("80052113","1111100017032230000","01631031000000d0"); // EMC7
    cuts.AddCut("80085113","1111100017032230000","01631031000000d0"); // EG2
    cuts.AddCut("80083113","1111100017032230000","01631031000000d0"); // EG1
  } else if (trainConfig == 204){ // EMCAL clusters standard cuts, triggers, no nonlin, +-50ns
    cuts.AddCut("80110113","1111100057032230000","01631031000000d0"); // 0-10
    cuts.AddCut("81210113","1111100057032230000","01631031000000d0"); // 10-20
    cuts.AddCut("82410113","1111100057032230000","01631031000000d0"); // 20-40
    cuts.AddCut("84610113","1111100057032230000","01631031000000d0"); // 40-60
    cuts.AddCut("86810113","1111100057032230000","01631031000000d0"); // 60-80
    cuts.AddCut("88010113","1111100057032230000","01631031000000d0"); // 80-100
  } else if (trainConfig == 205){ // EMCAL clusters standard cuts, triggers, no nonlin, open timing
    cuts.AddCut("80110113","1111100017032230000","01631031000000d0"); // 0-10
    cuts.AddCut("81210113","1111100017032230000","01631031000000d0"); // 10-20
    cuts.AddCut("82410113","1111100017032230000","01631031000000d0"); // 20-40
    cuts.AddCut("84610113","1111100017032230000","01631031000000d0"); // 40-60
    cuts.AddCut("86810113","1111100017032230000","01631031000000d0"); // 60-80
    cuts.AddCut("88010113","1111100017032230000","01631031000000d0"); // 80-100
  } else if (trainConfig == 206){ // EMCAL clusters standard cuts, no nonlin, +-50ns
    cuts.AddCut("80210113","1111100057032230000","01631031000000d0"); // 0-20
    cuts.AddCut("a0110113","1111100057032230000","01631031000000d0"); // 0-5
    cuts.AddCut("a1210113","1111100057032230000","01631031000000d0"); // 5-10
    cuts.AddCut("86010113","1111100057032230000","01631031000000d0"); // 60-100
  } else if (trainConfig == 207){ // EMCAL clusters standard cuts, no nonlin, open timing
    cuts.AddCut("80210113","1111100017032230000","01631031000000d0"); // 0-20
    cuts.AddCut("a0110113","1111100017032230000","01631031000000d0"); // 0-5
    cuts.AddCut("a1210113","1111100017032230000","01631031000000d0"); // 5-10
    cuts.AddCut("86010113","1111100017032230000","01631031000000d0"); // 60-100

  // non lin variations
  } else if (trainConfig == 208){ // EMCAL clusters standard cuts, non lin variations
    cuts.AddCut("80010113","1111141057032230000","01631031000000d0"); // 0-100
    cuts.AddCut("80010113","1111142057032230000","01631031000000d0"); // 0-100
    cuts.AddCut("80010113","1111151057032230000","01631031000000d0"); // 0-100
    cuts.AddCut("80010113","1111152057032230000","01631031000000d0"); // 0-100
    cuts.AddCut("80010113","1111155057032230000","01631031000000d0"); // 0-100
    cuts.AddCut("80010113","1111156057032230000","01631031000000d0"); // 0-100
  } else if (trainConfig == 209){ // EMCAL clusters standard cuts, non lin variations
    cuts.AddCut("80010113","1111102057032230000","01631031000000d0"); // 0-100
    cuts.AddCut("80010113","1111100057032230000","01631031000000d0"); // 0-100
  } else if (trainConfig == 210){ // EMCAL clusters standard cuts, non lin variations
    cuts.AddCut("80010113","1111141017032230000","01631031000000d0"); // 0-100
    cuts.AddCut("80010113","1111142017032230000","01631031000000d0"); // 0-100
    cuts.AddCut("80010113","1111151017032230000","01631031000000d0"); // 0-100
    cuts.AddCut("80010113","1111152017032230000","01631031000000d0"); // 0-100
  } else if (trainConfig == 211){ // EMCAL clusters standard cuts, non lin variations
    cuts.AddCut("80010113","1111102017032230000","01631031000000d0"); // 0-100
    cuts.AddCut("80010113","1111100017032230000","01631031000000d0"); // 0-100

  // non lin variations
  } else if (trainConfig == 212){ // EMCAL clusters standard cuts, cent vars, non lin variations
    cuts.AddCut("80110113","1111141057032230000","01631031000000d0"); // 0-10
    cuts.AddCut("81210113","1111141057032230000","01631031000000d0"); // 10-20
    cuts.AddCut("82410113","1111141057032230000","01631031000000d0"); // 20-40
    cuts.AddCut("84610113","1111141057032230000","01631031000000d0"); // 40-60
    cuts.AddCut("86810113","1111141057032230000","01631031000000d0"); // 60-80
    cuts.AddCut("88010113","1111141057032230000","01631031000000d0"); // 80-100
  } else if (trainConfig == 213){ // EMCAL clusters standard cuts, cent vars, non lin variations
    cuts.AddCut("80210113","1111141057032230000","01631031000000d0"); // 0-20
    cuts.AddCut("86010113","1111141057032230000","01631031000000d0"); // 60-100
    cuts.AddCut("a0110113","1111141057032230000","01631031000000d0"); // 0-5
    cuts.AddCut("a1210113","1111141057032230000","01631031000000d0"); // 5-10
  } else if (trainConfig == 214){ // EMCAL clusters standard cuts, cent vars, non lin variations
    cuts.AddCut("80110113","1111151057032230000","01631031000000d0"); // 0-10
    cuts.AddCut("81210113","1111151057032230000","01631031000000d0"); // 10-20
    cuts.AddCut("82410113","1111151057032230000","01631031000000d0"); // 20-40
    cuts.AddCut("84610113","1111151057032230000","01631031000000d0"); // 40-60
    cuts.AddCut("86810113","1111151057032230000","01631031000000d0"); // 60-80
    cuts.AddCut("88010113","1111151057032230000","01631031000000d0"); // 80-100
  } else if (trainConfig == 215){ // EMCAL clusters standard cuts, cent vars, non lin variations
    cuts.AddCut("80210113","1111151057032230000","01631031000000d0"); // 0-20
    cuts.AddCut("86010113","1111151057032230000","01631031000000d0"); // 60-100
    cuts.AddCut("a0110113","1111151057032230000","01631031000000d0"); // 0-5
    cuts.AddCut("a1210113","1111151057032230000","01631031000000d0"); // 5-10
  } else if (trainConfig == 216){ // EMCAL clusters standard cuts, dir gamma cut
    cuts.AddCut("80010113","11111410570322l0000","01631031000000d0"); // 0-100
    cuts.AddCut("80010113","11111420570322l0000","01631031000000d0"); // 0-100
    cuts.AddCut("80010113","11111510570322l0000","01631031000000d0"); // 0-100
    cuts.AddCut("80010113","11111520570322l0000","01631031000000d0"); // 0-100
    cuts.AddCut("80010113","11111020570322l0000","01631031000000d0"); // 0-100
  } else if (trainConfig == 217){ // EMCAL clusters standard cuts, dir gamma cut
    cuts.AddCut("80110113","11111410570322l0000","01631031000000d0");
    cuts.AddCut("81210113","11111410570322l0000","01631031000000d0");
    cuts.AddCut("82410113","11111410570322l0000","01631031000000d0");
    cuts.AddCut("84610113","11111410570322l0000","01631031000000d0");
    cuts.AddCut("86810113","11111410570322l0000","01631031000000d0");
    cuts.AddCut("88010113","11111410570322l0000","01631031000000d0");
  } else if (trainConfig == 218){ // EMCAL clusters standard cuts, dir gamma cut
    cuts.AddCut("80210113","11111410570322l0000","01631031000000d0"); // 0-20
    cuts.AddCut("86010113","11111410570322l0000","01631031000000d0"); // 60-100
    cuts.AddCut("a0110113","11111410570322l0000","01631031000000d0"); // 0-5
    cuts.AddCut("a1210113","11111410570322l0000","01631031000000d0"); // 5-10
  } else if (trainConfig == 219){ // EMCAL clusters standard cuts, dir gamma cut
    cuts.AddCut("80110113","11111510570322l0000","01631031000000d0");
    cuts.AddCut("81210113","11111510570322l0000","01631031000000d0");
    cuts.AddCut("82410113","11111510570322l0000","01631031000000d0");
    cuts.AddCut("84610113","11111510570322l0000","01631031000000d0");
    cuts.AddCut("86810113","11111510570322l0000","01631031000000d0");
    cuts.AddCut("88010113","11111510570322l0000","01631031000000d0");
  } else if (trainConfig == 220){ // EMCAL clusters standard cuts, dir gamma cut
    cuts.AddCut("80210113","11111510570322l0000","01631031000000d0"); // 0-20
    cuts.AddCut("86010113","11111510570322l0000","01631031000000d0"); // 60-100
    cuts.AddCut("a0110113","11111510570322l0000","01631031000000d0"); // 0-5
    cuts.AddCut("a1210113","11111510570322l0000","01631031000000d0"); // 5-10

  // Nonlin variations with different min E cut
  } else if (trainConfig == 221){ // EMCAL clusters, non lin variations, min E = 0.6
    cuts.AddCut("80010113","1111100017022230000","01631031000000d0"); // 0-100
    cuts.AddCut("80010113","1111141017022230000","01631031000000d0"); // 0-100
    cuts.AddCut("80010113","1111142017022230000","01631031000000d0"); // 0-100
    cuts.AddCut("80010113","1111151017022230000","01631031000000d0"); // 0-100
    cuts.AddCut("80010113","1111152017022230000","01631031000000d0"); // 0-100
  } else if (trainConfig == 222){ // EMCAL clusters, non lin variations, min E = 0.65
    cuts.AddCut("80010113","11111000170b2230000","01631031000000d0"); // 0-100
    cuts.AddCut("80010113","11111410170b2230000","01631031000000d0"); // 0-100
    cuts.AddCut("80010113","11111420170b2230000","01631031000000d0"); // 0-100
    cuts.AddCut("80010113","11111510170b2230000","01631031000000d0"); // 0-100
    cuts.AddCut("80010113","11111520170b2230000","01631031000000d0"); // 0-100
  } else if (trainConfig == 223){ // EMCAL clusters, non lin variations, min E = 0.675
    cuts.AddCut("80010113","11111000170c2230000","01631031000000d0"); // 0-100
    cuts.AddCut("80010113","11111410170c2230000","01631031000000d0"); // 0-100
    cuts.AddCut("80010113","11111420170c2230000","01631031000000d0"); // 0-100
    cuts.AddCut("80010113","11111510170c2230000","01631031000000d0"); // 0-100
    cuts.AddCut("80010113","11111520170c2230000","01631031000000d0"); // 0-100
  } else if (trainConfig == 224){ // EMCAL clusters, non lin variations, min E = 0.625
    cuts.AddCut("80010113","11111000170d2230000","01631031000000d0"); // 0-100
    cuts.AddCut("80010113","11111410170d2230000","01631031000000d0"); // 0-100
    cuts.AddCut("80010113","11111420170d2230000","01631031000000d0"); // 0-100
    cuts.AddCut("80010113","11111510170d2230000","01631031000000d0"); // 0-100
    cuts.AddCut("80010113","11111520170d2230000","01631031000000d0"); // 0-100

  // Small centrality intervals
  } else if (trainConfig == 225){
    cuts.AddCut("c0110113","1111141057032230000","01631031000000d0"); // 0-1
    cuts.AddCut("c0210113","1111141057032230000","01631031000000d0"); // 0-2
    cuts.AddCut("c0310113","1111141057032230000","01631031000000d0"); // 0-3
    cuts.AddCut("c0410113","1111141057032230000","01631031000000d0"); // 0-4
  // AOD validation
  } else if (trainConfig == 226){
    cuts.AddCut("80010113","1111151017032230000","01631031000000d0"); // 0-100
  // pPb 8 TeV EPOS+PythiaJets JJ simulation QA
  } else if (trainConfig == 227){ // same as 201 but with all headers
    cuts.AddCut("80010103","1111100017032230000","01631031000000d0"); // INT7
  } else if (trainConfig == 228){ // same as 201 but with special header(s)
    cuts.AddCut("80010123","1111100017032230000","01631031000000d0"); // INT7
  //-----------------------------------------------------------------------------------------------
  // Systematics variations MB run 2 EMC pPb std cut: cuts.AddCut("80010113","1111141057032230000","01631031000000d0"); // 0-100
  //-----------------------------------------------------------------------------------------------
  } else if (trainConfig == 230){ // set of variations CLUSTER
    cuts.AddCut("80010113","1111141057012230000","01631031000000d0"); // min energy cluster variation 1  500 MeV
    cuts.AddCut("80010113","1111141057022230000","01631031000000d0"); // min energy cluster variation 1  600 MeV
    cuts.AddCut("80010113","1111141057042230000","01631031000000d0"); // min energy cluster variation 2  800 MeV
    cuts.AddCut("80010113","1111141057052230000","01631031000000d0"); // min energy cluster variation 3  900 MeV
  } else if (trainConfig == 231){ // set of variations CLUSTER
    cuts.AddCut("80010113","1111141057032200000","01631031000000d0"); // min/max M02  0.1<M<100
    cuts.AddCut("80010113","1111141057032210000","01631031000000d0"); // min/max M02  0.1<M<1
    cuts.AddCut("80010113","1111141057032220000","01631031000000d0"); // min/max M02  0.1<M<0.7
    cuts.AddCut("80010113","1111141057032240000","01631031000000d0"); // min/max M02  0.1<M<0.4
    cuts.AddCut("80010113","1111141057032250000","01631031000000d0"); // min/max M02  0.1<M<0.3
    cuts.AddCut("80010113","1111141057032260000","01631031000000d0"); // min/max M02  0.1<M<0.27
  } else if (trainConfig == 232){ // set of variations CLUSTER
    cuts.AddCut("80010113","1111141057031230000","01631031000000d0"); // min number of cells variation 1  1 cell
    cuts.AddCut("80010113","1112141057032230000","01631031000000d0"); // only modules with TRD infront
    cuts.AddCut("80010113","1111341057032230000","01631031000000d0"); // no modules with TRD infront
  } else if (trainConfig == 233){ // set of variations MESON
    cuts.AddCut("80010113","1111141057032230000","01633031000000d0"); // rapidity variation  y<0.6
    cuts.AddCut("80010113","1111141057032230000","01634031000000d0"); // rapidity variation  y<0.5
    cuts.AddCut("80010113","1111141057032230000","01631061000000d0"); // alpha meson variation 1   0<alpha<0.8
    cuts.AddCut("80010113","1111141057032230000","01631051000000d0"); // alpha meson variation 2  0<alpha<0.75
  } else if (trainConfig == 234){ // opening angle variations
    cuts.AddCut("80010113","1111141057032230000","0163103100000040"); // min opening angle 0.0152
    cuts.AddCut("80010113","1111141057032230000","0163103100000060"); // min opening angle 0.017
    cuts.AddCut("80010113","1111141057032230000","0163103100000070"); // min opening angle 0.016
    cuts.AddCut("80010113","1111141057032230000","0163103100000080"); // min opening angle 0.018
    cuts.AddCut("80010113","1111141057032230000","0163103100000090"); // min opening angle 0.018
  } else if (trainConfig == 235){ // TM variations
    cuts.AddCut("80010113","1111141053032230000","01631031000000d0"); // fixed window
    cuts.AddCut("80010113","1111141056032230000","01631031000000d0"); // tm pt dependent var 1
    cuts.AddCut("80010113","1111141058032230000","01631031000000d0"); // tm pt dependent var 2
    cuts.AddCut("80010113","1111141059032230000","01631031000000d0"); // tm pt dependent var 3
  } else if (trainConfig == 236){ // EMCAL clusters standard cuts, non lin variations
    cuts.AddCut("80010113","1111142057032230000","01631031000000d0"); // 0-100
    cuts.AddCut("80010113","1111151057032230000","01631031000000d0"); // 0-100
    cuts.AddCut("80010113","1111152057032230000","01631031000000d0"); // 0-100
    cuts.AddCut("80010113","1111155057032230000","01631031000000d0"); // 0-100
    cuts.AddCut("80010113","1111156057032230000","01631031000000d0"); // 0-100

  //-----------------------------------------------------------------------------------------------
  // Systematics variations 0-10% run 2 EMC pPb std cut: cuts.AddCut("80010113","1111141057032230000","01631031000000d0"); // 0-10
  //-----------------------------------------------------------------------------------------------
  } else if (trainConfig == 240){ // set of variations CLUSTER
    cuts.AddCut("80110113","1111141057012230000","01631031000000d0"); // min energy cluster variation 1  500 MeV
    cuts.AddCut("80110113","1111141057022230000","01631031000000d0"); // min energy cluster variation 1  600 MeV
    cuts.AddCut("80110113","1111141057042230000","01631031000000d0"); // min energy cluster variation 2  800 MeV
    cuts.AddCut("80110113","1111141057052230000","01631031000000d0"); // min energy cluster variation 3  900 MeV
  } else if (trainConfig == 241){ // set of variations CLUSTER
    cuts.AddCut("80110113","1111141057032200000","01631031000000d0"); // min/max M02  0.1<M<100
    cuts.AddCut("80110113","1111141057032210000","01631031000000d0"); // min/max M02  0.1<M<1
    cuts.AddCut("80110113","1111141057032220000","01631031000000d0"); // min/max M02  0.1<M<0.7
    cuts.AddCut("80110113","1111141057032240000","01631031000000d0"); // min/max M02  0.1<M<0.4
    cuts.AddCut("80110113","1111141057032250000","01631031000000d0"); // min/max M02  0.1<M<0.3
    cuts.AddCut("80110113","1111141057032260000","01631031000000d0"); // min/max M02  0.1<M<0.27
  } else if (trainConfig == 242){ // set of variations CLUSTER
    cuts.AddCut("80110113","1111141057031230000","01631031000000d0"); // min number of cells variation 1  1 cell
    cuts.AddCut("80110113","1112141057032230000","01631031000000d0"); // only modules with TRD infront
    cuts.AddCut("80110113","1111341057032230000","01631031000000d0"); // no modules with TRD infront
  } else if (trainConfig == 243){ // set of variations MESON
    cuts.AddCut("80110113","1111141057032230000","01633031000000d0"); // rapidity variation  y<0.6
    cuts.AddCut("80110113","1111141057032230000","01634031000000d0"); // rapidity variation  y<0.5
    cuts.AddCut("80110113","1111141057032230000","01631061000000d0"); // alpha meson variation 1   0<alpha<0.8
    cuts.AddCut("80110113","1111141057032230000","01631051000000d0"); // alpha meson variation 2  0<alpha<0.75
  } else if (trainConfig == 244){ // opening angle variations
    cuts.AddCut("80110113","1111141057032230000","0163103100000040"); // min opening angle 0.0152
    cuts.AddCut("80110113","1111141057032230000","0163103100000060"); // min opening angle 0.017
    cuts.AddCut("80110113","1111141057032230000","0163103100000070"); // min opening angle 0.016
    cuts.AddCut("80110113","1111141057032230000","0163103100000080"); // min opening angle 0.018
    cuts.AddCut("80110113","1111141057032230000","0163103100000090"); // min opening angle 0.018
  } else if (trainConfig == 245){ // TM variations
    cuts.AddCut("80110113","1111141053032230000","01631031000000d0"); // fixed window
    cuts.AddCut("80110113","1111141056032230000","01631031000000d0"); // tm pt dependent var 1
    cuts.AddCut("80110113","1111141058032230000","01631031000000d0"); // tm pt dependent var 2
    cuts.AddCut("80110113","1111141059032230000","01631031000000d0"); // tm pt dependent var 3
  } else if (trainConfig == 246){ // EMCAL clusters standard cuts, non lin variations
    cuts.AddCut("80110113","1111142057032230000","01631031000000d0"); // 0-100
    cuts.AddCut("80110113","1111151057032230000","01631031000000d0"); // 0-100
    cuts.AddCut("80110113","1111152057032230000","01631031000000d0"); // 0-100
    cuts.AddCut("80110113","1111155057032230000","01631031000000d0"); // 0-100
    cuts.AddCut("80110113","1111156057032230000","01631031000000d0"); // 0-100

  //-----------------------------------------------------------------------------------------------
  // Systematics variations 80-100% run 2 EMC pPb std cut: cuts.AddCut("80010113","1111141057032230000","01631031000000d0"); // 80-100
  //-----------------------------------------------------------------------------------------------
  } else if (trainConfig == 250){ // set of variations CLUSTER
    cuts.AddCut("88010113","1111141057012230000","01631031000000d0"); // min energy cluster variation 1  500 MeV
    cuts.AddCut("88010113","1111141057022230000","01631031000000d0"); // min energy cluster variation 1  600 MeV
    cuts.AddCut("88010113","1111141057042230000","01631031000000d0"); // min energy cluster variation 2  800 MeV
    cuts.AddCut("88010113","1111141057052230000","01631031000000d0"); // min energy cluster variation 3  900 MeV
  } else if (trainConfig == 251){ // set of variations CLUSTER
    cuts.AddCut("88010113","1111141057032200000","01631031000000d0"); // min/max M02  0.1<M<100
    cuts.AddCut("88010113","1111141057032210000","01631031000000d0"); // min/max M02  0.1<M<1
    cuts.AddCut("88010113","1111141057032220000","01631031000000d0"); // min/max M02  0.1<M<0.7
    cuts.AddCut("88010113","1111141057032240000","01631031000000d0"); // min/max M02  0.1<M<0.4
    cuts.AddCut("88010113","1111141057032250000","01631031000000d0"); // min/max M02  0.1<M<0.3
    cuts.AddCut("88010113","1111141057032260000","01631031000000d0"); // min/max M02  0.1<M<0.27
  } else if (trainConfig == 252){ // set of variations CLUSTER
    cuts.AddCut("88010113","1111141057031230000","01631031000000d0"); // min number of cells variation 1  1 cell
    cuts.AddCut("88010113","1112141057032230000","01631031000000d0"); // only modules with TRD infront
    cuts.AddCut("88010113","1111341057032230000","01631031000000d0"); // no modules with TRD infront
  } else if (trainConfig == 253){ // set of variations MESON
    cuts.AddCut("88010113","1111141057032230000","01633031000000d0"); // rapidity variation  y<0.6
    cuts.AddCut("88010113","1111141057032230000","01634031000000d0"); // rapidity variation  y<0.5
    cuts.AddCut("88010113","1111141057032230000","01631061000000d0"); // alpha meson variation 1   0<alpha<0.8
    cuts.AddCut("88010113","1111141057032230000","01631051000000d0"); // alpha meson variation 2  0<alpha<0.75
  } else if (trainConfig == 254){ // opening angle variations
    cuts.AddCut("88010113","1111141057032230000","0163103100000040"); // min opening angle 0.0152
    cuts.AddCut("88010113","1111141057032230000","0163103100000060"); // min opening angle 0.017
    cuts.AddCut("88010113","1111141057032230000","0163103100000070"); // min opening angle 0.016
    cuts.AddCut("88010113","1111141057032230000","0163103100000080"); // min opening angle 0.018
    cuts.AddCut("88010113","1111141057032230000","0163103100000090"); // min opening angle 0.018
  } else if (trainConfig == 255){ // TM variations
    cuts.AddCut("88010113","1111141053032230000","01631031000000d0"); // fixed window
    cuts.AddCut("88010113","1111141056032230000","01631031000000d0"); // tm pt dependent var 1
    cuts.AddCut("88010113","1111141058032230000","01631031000000d0"); // tm pt dependent var 2
    cuts.AddCut("88010113","1111141059032230000","01631031000000d0"); // tm pt dependent var 3
  } else if (trainConfig == 256){ // EMCAL clusters standard cuts, non lin variations
    cuts.AddCut("88010113","1111142057032230000","01631031000000d0"); // 0-100
    cuts.AddCut("88010113","1111151057032230000","01631031000000d0"); // 0-100
    cuts.AddCut("88010113","1111152057032230000","01631031000000d0"); // 0-100
    cuts.AddCut("88010113","1111155057032230000","01631031000000d0"); // 0-100
    cuts.AddCut("88010113","1111156057032230000","01631031000000d0"); // 0-100

  // ===============================================================================================
  // Run 1 data PHOS clusters pPb 5TeV
  // ===============================================================================================
  } else if (trainConfig == 301) {  // min energy = 0.3 GeV/c
    cuts.AddCut("80010113","2444400041013200000","0163103100000010"); //standart cut, kINT7 // PHOS clusters
    cuts.AddCut("80062113","2444400041013200000","0163103100000010"); //standard cut, kPHI7  // PHOS clusters
  } else if (trainConfig == 302){ // Validation PHOS
    cuts.AddCut("80010113","2444400041013200000","0163103100000010");
  } else if (trainConfig == 303){ // Validation PHOS, only added signals
    cuts.AddCut("80010023","2444400041013200000","0163103100000010");
  } else if (trainConfig == 304){ // min energy = 0.3 GeV/c
    cuts.AddCut("80010113","2444400041013200000","0163103100000000"); // kINT7 // PHOS clusters no open angle cut
    cuts.AddCut("80062113","2444400041013200000","0163103100000000"); // kPHI7 // PHOS clusters no open angle cut
    cuts.AddCut("80010113","2444400041013200000","0163103100000030"); // kINT7 // PHOS clusters open angle cut 0.1
    cuts.AddCut("80062113","2444400041013200000","0163103100000030"); // kPHI7 // PHOS clusters  open angle cut 0.1
  } else if (trainConfig == 305){ // timing cut variations
    cuts.AddCut("80010113","2444400011013200000","0163103100000010"); // 1000ns
    cuts.AddCut("80010113","2444400031013200000","0163103100000010"); // 200ns
    cuts.AddCut("80010113","2444400051013200000","0163103100000010"); // 50ns
  } else if (trainConfig == 306) {  // PHOS non lin var INT7
    cuts.AddCut("80010113","2444401041013200000","0163103100000010"); // PHOS group standard
    cuts.AddCut("80010113","2444451041013200000","0163103100000010"); // CCMF PHOS
    cuts.AddCut("80010113","2444452041013200000","0163103100000010"); // CMF PHOS
  } else if (trainConfig == 307) {  // PHOS non lin var PHI7
    cuts.AddCut("80062113","2444401041013200000","0163103100000010"); // PHOS group standard
    cuts.AddCut("80062113","2444451041013200000","0163103100000010"); // CCMF PHOS
    cuts.AddCut("80062113","2444452041013200000","0163103100000010"); // CMF PHOS
  } else if (trainConfig == 308) {  // PHOS CCMF cent dep
    cuts.AddCut("80210113","2444451041013200000","0163103100000010"); // 0-20
    cuts.AddCut("82410113","2444451041013200000","0163103100000010"); // 20-40
    cuts.AddCut("84610113","2444451041013200000","0163103100000010"); // 40-60
    cuts.AddCut("86010113","2444451041013200000","0163103100000010"); // 60-100
  } else if (trainConfig == 309) {  // PHOS default cent dep
    cuts.AddCut("80210113","2444401041013200000","0163103100000010"); // 0-20
    cuts.AddCut("82410113","2444401041013200000","0163103100000010"); // 20-40
    cuts.AddCut("84610113","2444401041013200000","0163103100000010"); // 40-60
    cuts.AddCut("86010113","2444401041013200000","0163103100000010"); // 60-100

  } else if(trainConfig == 310){ // first set of variations CLUSTER
    cuts.AddCut("80010113","2444451041073200000","0163103100000010"); // min energy 0.2 GeV
    cuts.AddCut("80010113","2444451041083200000","0163103100000010"); // min energy 0.4 GeV
    cuts.AddCut("80010113","2444451041023200000","0163103100000010"); // min energy 0.5 GeV
  } else if(trainConfig == 311){ // second set of variations CLUSTER
    cuts.AddCut("80010113","2444451041013270000","0163103100000010"); // min/max M02  0.1<M<1
    cuts.AddCut("80010113","2444451041013280000","0163103100000010"); // min/max M02  0.1<M<1
    cuts.AddCut("80010113","2444451041012200000","0163103100000010"); // min number 2 cells
    cuts.AddCut("80010113","2444451041014200000","0163103100000010"); // min number 4 cells
  } else if(trainConfig == 312){ // MESON
    cuts.AddCut("80010113","2444451041013200000","0163403100000010"); // rapidity variation  y<0.5
    cuts.AddCut("80010113","2444451041013200000","0163803100000010"); // rapidity variation  y<0.25
    cuts.AddCut("80010113","2444451041013200000","0163106100000010"); // alpha meson variation 1 0<alpha<0.8
    cuts.AddCut("80010113","2444451041013200000","0163105100000010"); // alpha meson variation 2 0<alpha<0.75
  } else if(trainConfig == 313){ // fourth set of variations
    cuts.AddCut("80010113","2444451044013200000","0163103100000010"); // tm variation
    cuts.AddCut("80010113","2444451045013200000","0163103100000010"); // tm variation
    cuts.AddCut("80010113","2444451046013200000","0163103100000010"); // tm variation
    cuts.AddCut("80010113","2444451041013200000","0163103100000000"); // min opening angle 0    -> open
    cuts.AddCut("80010113","2444451041013200000","0163103100000030"); // min opening angle 0.01 -> 2 cell diag

  // ===============================================================================================
  // Run 2 data EMC clusters pPb 8TeV
  // ===============================================================================================
  // pPb 8TeV variations for QA
  } else if (trainConfig == 400){ // EMCAL clusters standard cuts, triggers, no nonlin, open timing (same as 201 but split)
    cuts.AddCut("80010113","1111100017032230000","01631031000000d0"); // INT7
  } else if (trainConfig == 401){ // EMCAL clusters standard cuts, triggers, no nonlin, open timing (same as 201 but split)
    cuts.AddCut("80052113","1111100017032230000","01631031000000d0"); // EMC7
  } else if (trainConfig == 402){ // EMCAL clusters standard cuts, triggers, no nonlin, open timing (same as 201 but split)
    cuts.AddCut("80085113","1111100017032230000","01631031000000d0"); // EG2
  } else if (trainConfig == 403){ // EMCAL clusters standard cuts, triggers, no nonlin, open timing (same as 201 but split)
    cuts.AddCut("80083113","1111100017032230000","01631031000000d0"); // EG1

  // ===============================================================================================
  // Run 2 data PHOS clusters pPb 5TeV
  // ===============================================================================================
  // INT7 triggers
  } else if (trainConfig == 500) {  // PHOS  INT7
    cuts.AddCut("80010113","2446600051013200000","0163103100000010"); // no non lin 0-100%
  } else if (trainConfig == 501) {  // PHOS  INT7
    cuts.AddCut("80010113","2446600041013200000","0163103100000010"); // no non lin
    cuts.AddCut("80010113","2446600011013200000","0163103100000010"); // no non lin 1000 \mus
    cuts.AddCut("80010113","2446600061013200000","0163103100000010"); // no non lin, -30, 50ns
    cuts.AddCut("80010113","24466000a1013200000","0163103100000010"); // no non lin, -12.5, 13ns
  } else if (trainConfig == 502) {  // PHOS  INT7 non lin vars
    cuts.AddCut("80010113","2446601051013200000","0163103100000010"); //
    cuts.AddCut("80010113","2446600051013200000","0163103100000010"); //
    cuts.AddCut("80010113","2446651051013200000","0163103100000010"); //
    cuts.AddCut("80010113","2446641051013200000","0163103100000010"); //
  } else if (trainConfig == 503) {  // PHOS  INT7 with cents
    cuts.AddCut("80110113","2446600051013200000","0163103100000010"); // no non lin 0-10%
    cuts.AddCut("81210113","2446600051013200000","0163103100000010"); // no non lin 10-20%
    cuts.AddCut("82410113","2446600051013200000","0163103100000010"); // no non lin 20-40%
    cuts.AddCut("84610113","2446600051013200000","0163103100000010"); // no non lin 40-60%
    cuts.AddCut("86810113","2446600051013200000","0163103100000010"); // no non lin 60-80%
    cuts.AddCut("88010113","2446600051013200000","0163103100000010"); // no non lin 80-100%
  } else if (trainConfig == 504) {  // PHOS  INT7 with cents
    cuts.AddCut("80210113","2446600051013200000","0163103100000010"); // no non lin 0-10%
    cuts.AddCut("86010113","2446600051013200000","0163103100000010"); // no non lin 60-100%
    cuts.AddCut("a0110113","2446600051013200000","0163103100000010"); // no non lin 0-5%
    cuts.AddCut("a1210113","2446600051013200000","0163103100000010"); // no non lin 5-10%
  } else if (trainConfig == 505) {  // PHOS  INT7 with cents
    cuts.AddCut("80110113","2446601051013200000","0163103100000010"); // no non lin 0-10%
    cuts.AddCut("81210113","2446601051013200000","0163103100000010"); // no non lin 10-20%
    cuts.AddCut("82410113","2446601051013200000","0163103100000010"); // no non lin 20-40%
    cuts.AddCut("84610113","2446601051013200000","0163103100000010"); // no non lin 40-60%
    cuts.AddCut("86810113","2446601051013200000","0163103100000010"); // no non lin 60-80%
    cuts.AddCut("88010113","2446601051013200000","0163103100000010"); // no non lin 80-100%
  } else if (trainConfig == 506) {  // PHOS  INT7 with cents
    cuts.AddCut("80210113","2446601051013200000","0163103100000010"); // no non lin 0-10%
    cuts.AddCut("86010113","2446601051013200000","0163103100000010"); // no non lin 60-100%
    cuts.AddCut("a0110113","2446601051013200000","0163103100000010"); // no non lin 0-5%
    cuts.AddCut("a1210113","2446601051013200000","0163103100000010"); // no non lin 5-10%
  } else if (trainConfig == 507){ // AOD validation
    cuts.AddCut("80010113","2446641051013200000","0163103100000010"); //

  // ===============================================================================================
  // Run 2 data PHOS clusters pPb 8TeV
  // ===============================================================================================
  // pPb 8TeV variations for QA
  } else if (trainConfig == 600){ // PHOS clusters standard cuts, triggers, no nonlin, +-50ns
    cuts.AddCut("80010113","2446600051013200000","0163103100000010"); // INT7
  } else if (trainConfig == 601){ // PHOS clusters standard cuts, triggers, no nonlin, +-50ns
    cuts.AddCut("80062113","2446600051013200000","0163103100000010"); // PHI7

  // ===============================================================================================
  // Run 2 data DMC clusters pPb 5TeV
  // ===============================================================================================
  // pPb 8TeV variations for QA
  } else if (trainConfig == 700){ // DCal clusters standard cuts, triggers, no nonlin, open timing
    cuts.AddCut("80010113","3885500017032230000","01631031000000d0"); // INT7
  } else if (trainConfig == 701){ // DCal clusters standard cuts, triggers, no nonlin, open timing
    cuts.AddCut("80010113","3885500057032230000","01631031000000d0"); // EMC7
  } else if (trainConfig == 702){ // DCal clusters standard cuts, triggers, no nonlin, +-50ns
    cuts.AddCut("80110113","3885500057032230000","01631031000000d0"); // 0-10
    cuts.AddCut("81210113","3885500057032230000","01631031000000d0"); // 10-20
    cuts.AddCut("82410113","3885500057032230000","01631031000000d0"); // 20-40
    cuts.AddCut("84610113","3885500057032230000","01631031000000d0"); // 40-60
    cuts.AddCut("86810113","3885500057032230000","01631031000000d0"); // 60-80
    cuts.AddCut("88010113","3885500057032230000","01631031000000d0"); // 80-100
  } else if (trainConfig == 703){ // DCal clusters standard cuts, triggers, no nonlin, open timing
    cuts.AddCut("80110113","3885500017032230000","01631031000000d0"); // 0-10
    cuts.AddCut("81210113","3885500017032230000","01631031000000d0"); // 10-20
    cuts.AddCut("82410113","3885500017032230000","01631031000000d0"); // 20-40
    cuts.AddCut("84610113","3885500017032230000","01631031000000d0"); // 40-60
    cuts.AddCut("86810113","3885500017032230000","01631031000000d0"); // 60-80
    cuts.AddCut("88010113","3885500017032230000","01631031000000d0"); // 80-100
  } else if (trainConfig == 704){ // DCal clusters standard cuts, no nonlin, +-50ns
    cuts.AddCut("80210113","3885500057032230000","01631031000000d0"); // 0-20
    cuts.AddCut("a0110113","3885500057032230000","01631031000000d0"); // 0-5
    cuts.AddCut("a1210113","3885500057032230000","01631031000000d0"); // 5-10
    cuts.AddCut("86010113","3885500057032230000","01631031000000d0"); // 60-100
  } else if (trainConfig == 705){ // DCal clusters standard cuts, no nonlin, open timing
    cuts.AddCut("80210113","3885500017032230000","01631031000000d0"); // 0-20
    cuts.AddCut("a0110113","3885500017032230000","01631031000000d0"); // 0-5
    cuts.AddCut("a1210113","3885500017032230000","01631031000000d0"); // 5-10
    cuts.AddCut("86010113","3885500017032230000","01631031000000d0"); // 60-100

  // ===============================================================================================
  // Run 2 data DMC clusters pPb 8TeV
  // ===============================================================================================
  // pPb 8TeV variations for QA
  } else if (trainConfig == 800){ // DCal clusters standard cuts, triggers, no nonlin, open timing
    cuts.AddCut("80010113","3885500017032230000","01631031000000d0"); // INT7
  } else if (trainConfig == 801){ // DCal clusters standard cuts, triggers, no nonlin, open timing
    cuts.AddCut("00055113","3885500017032230000","01631031000000d0"); // EMC7
    cuts.AddCut("00089113","3885500017032230000","01631031000000d0"); // EG2
    cuts.AddCut("0008b113","3885500017032230000","01631031000000d0"); // EG1
  } else if (trainConfig == 804){ // DCal clusters standard cuts, triggers, no nonlin, open timing
    cuts.AddCut("80010113","3885500057032230000","01631031000000d0"); // INT7
  } else if (trainConfig == 805){ // DCal clusters standard cuts, triggers, no nonlin, open timing
    cuts.AddCut("00055113","3885500057032230000","01631031000000d0"); // EMC7
    cuts.AddCut("00089113","3885500057032230000","01631031000000d0"); // EG2
    cuts.AddCut("0008b113","3885500057032230000","01631031000000d0"); // EG1


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
  if (periodNameV0Reader.Contains("LHC18b9")){
    TObjString *HeaderP8J = new TObjString("Pythia8Jets_1");
    HeaderList->Add(HeaderP8J);
  }

  EventCutList->SetOwner(kTRUE);
  AliConvEventCuts **analysisEventCuts        = new AliConvEventCuts*[numberOfCuts];
  ClusterCutList->SetOwner(kTRUE);
  AliCaloPhotonCuts **analysisClusterCuts     = new AliCaloPhotonCuts*[numberOfCuts];
  MesonCutList->SetOwner(kTRUE);
  AliConversionMesonCuts **analysisMesonCuts  = new AliConversionMesonCuts*[numberOfCuts];

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

    TString dataInputMultHisto  = "";
    TString mcInputMultHisto    = "";
    TString triggerString       = cuts.GetEventCut(i);
    triggerString               = triggerString(3,2);
    dataInputMultHisto          = Form("%s_%s", periodNameAnchor.Data(), triggerString.Data());
    mcInputMultHisto            = Form("%s_%s", periodNameV0Reader.Data(), triggerString.Data());

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
  task->SetCorrectionTaskSetting(corrTaskSetting);
  task->SetDoMesonAnalysis(kTRUE);
  task->SetDoMesonQA(enableQAMesonTask); //Attention new switch for Pi0 QA
  task->SetDoClusterQA(enableQAClusterTask);  //Attention new switch small for Cluster QA
  task->SetDoTHnSparse(isUsingTHnSparse);
  task->SetProduceTreeEOverP(doTreeEOverP);
  task->SetEnableSortingOfMCClusLabels(enableSortingMCLabels);
  if(enableExtMatchAndQA > 1){ task->SetPlotHistsExtQA(kTRUE);}

  //connect containers
  AliAnalysisDataContainer *coutput =
    mgr->CreateContainer(!(corrTaskSetting.CompareTo("")) ? Form("GammaCalo_%i",trainConfig) : Form("GammaCalo_%i_%s",trainConfig,corrTaskSetting.Data()), TList::Class(),
              AliAnalysisManager::kOutputContainer, Form("GammaCalo_%i.root",trainConfig) );

  mgr->AddTask(task);
  mgr->ConnectInput(task,0,cinput);
  mgr->ConnectOutput(task,1,coutput);

  return;

}
