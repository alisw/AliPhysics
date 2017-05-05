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
      }
    }
  }
  TString sAdditionalTrainConfig = rAdditionalTrainConfig->GetString();
  if (sAdditionalTrainConfig.Atoi() > 0){
    trainConfig = trainConfig + sAdditionalTrainConfig.Atoi();
    cout << "INFO: AddTask_GammaCalo_pPb running additionalTrainConfig '" << sAdditionalTrainConfig.Atoi() << "', train config: '" << trainConfig << "'" << endl;
  }

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
  
  //************************************************ EMCAL clusters *************************************************
  if (trainConfig == 1){ // no non lin
    cuts.AddCut("80000013","1111100057032230000","0163103100000050");
    cuts.AddCut("80052013","1111100057032230000","0163103100000050");
    cuts.AddCut("80085013","1111100057032230000","0163103100000050");
    cuts.AddCut("80083013","1111100057032230000","0163103100000050");
  } else if (trainConfig == 2){ // no non lin
    cuts.AddCut("80000013","1111100057032230000","0163103100000050");
    
  // SYSTEMATIC STUDY NEUTRAl MESON MEASUREMENTS MIKE SAS 26-04-2017
  } else if(trainConfig == 40){ // default cutstring and first set of variations nonlinearity
    cuts.AddCut("80000013","1111141057032230000","0163103100000050"); // default
    cuts.AddCut("80000013","1111142057032230000","0163103100000050"); // calo nonlinearity variation
    cuts.AddCut("80000013","1111151057032230000","0163103100000050"); // calo nonlinearity variation
    cuts.AddCut("80000013","1111152057032230000","0163103100000050"); // calo nonlinearity variation
  } else if(trainConfig == 41){ // second set of variations CLUSTER
    cuts.AddCut("80000013","1111141057022230000","0163103100000050"); // min energy cluster variation 1  600 MeV
    cuts.AddCut("80000013","1111141057042230000","0163103100000050"); // min energy cluster variation 2  800 MeV
    cuts.AddCut("80000013","1111141057052230000","0163103100000050"); // min energy cluster variation 3  900 MeV
  } else if(trainConfig == 42){ // third set of variations CLUSTER
    cuts.AddCut("80000013","1111141057032230000","0163103100000050"); // min/max M02  0.1<M<0.5
    cuts.AddCut("80000013","1111141057032200000","0163103100000050"); // min/max M02  0.1<M<100
    cuts.AddCut("80000013","1111141057032250000","0163103100000050"); // min/max M02  0.1<M<0.3
    cuts.AddCut("80000013","1111141057032260000","0163103100000050"); // min/max M02  0.1<M<0.27
    cuts.AddCut("80000013","1111141057031230000","0163103100000050"); // min number of cells variation 1  1 cell
    cuts.AddCut("80000013","1112141057032230000","0163103100000050"); // only modules with TRD infront
    cuts.AddCut("80000013","1111341057032230000","0163103100000050"); // no modules with TRD infront
  } else if(trainConfig == 43){ // third set of variations MESON
    cuts.AddCut("80000013","1111141057032230000","0163303100000050"); // rapidity variation  y<0.6
    cuts.AddCut("80000013","1111141057032230000","0163403100000050"); // rapidity variation  y<0.5
    cuts.AddCut("80000013","1111141057032230000","0163106100000050"); // alpha meson variation 1   0<alpha<0.8
    cuts.AddCut("80000013","1111141057032230000","0163105100000050"); // alpha meson variation 2  0<alpha<0.75
  } else if(trainConfig == 44){ // fourth set of variations
    cuts.AddCut("80000013","1111141053032230000","0163103100000050"); // tm variation
    cuts.AddCut("80000013","1111141055032230000","0163103100000050"); // tm variation
    cuts.AddCut("80000013","1111141053032230000","0163103100000040"); // min opening angle 0.75 cell diag
    cuts.AddCut("80000013","1111141053032230000","0163103100000050"); // min opening angle 1 cell diag
    cuts.AddCut("80000013","1111141053032230000","0163103100000060"); // min opening angle 2 cell diag
  } else if(trainConfig == 45){ // fourth set of variations
    cuts.AddCut("80000013","1111141057032230000","0163103100000050"); // default cut tm pt dependent
    cuts.AddCut("80000013","1111141056032230000","0163103100000050"); // tm pt dependent var 1
    cuts.AddCut("80000013","1111141058032230000","0163103100000050"); // tm pt dependent var 2
    cuts.AddCut("80000013","1111141059032230000","0163103100000050"); // tm pt dependent var 3
  } else if(trainConfig == 46){ // centrality dependent and NL variations
    cuts.AddCut("80200013","1111100057032230000","0163103100000050"); // 0-20% no NL
    cuts.AddCut("80200013","1111141057032230000","0163103100000050"); // 0-20% standard NL
    cuts.AddCut("82600013","1111100057032230000","0163103100000050"); // 20-60% no NL
    cuts.AddCut("82600013","1111141057032230000","0163103100000050"); // 20-60% standard NL
    cuts.AddCut("86000013","1111100057032230000","0163103100000050"); // 60-100% no NL
    cuts.AddCut("86000013","1111141057032230000","0163103100000050"); // 60-100% standard NL
  } else if(trainConfig == 47){ // new default cut
    cuts.AddCut("80000013","1111141057032230000","0163103100000050"); // default tm pt dependent
  } else if(trainConfig == 48){ // TM variations
    cuts.AddCut("80000013","1111142057032230000","0163103100000050"); // Calo nonlin tm pt dependent
    cuts.AddCut("80000013","1111143057032230000","0163103100000050"); // ConvCalo+testbeam nonlin tm pt dependent
    cuts.AddCut("80000013","1111144057032230000","0163103100000050"); // Calo+testbeam nonlin tm pt dependent
  } else if(trainConfig == 49){ // TM variations
    cuts.AddCut("80000013","1111151057032230000","0163103100000050"); // ConvCalo ratio nonlin tm pt dependent
    cuts.AddCut("80000013","1111152057032230000","0163103100000050"); // Calo ratio nonlin tm pt dependent
    cuts.AddCut("80000013","1111153057032230000","0163103100000050"); // ConvCalo ratio nonlin+testbeam tm pt dependent
    cuts.AddCut("80000013","1111154057032230000","0163103100000050"); // Calo ratio nonlin+testbeam tm pt dependent
    
  //EMC7
  } else if(trainConfig == 50){ // default cutstring and first set of variations nonlinearity
    cuts.AddCut("80052013","1111141057032230000","0163103100000050"); // default
    cuts.AddCut("80052013","1111142057032230000","0163103100000050"); // calo nonlinearity variation
    cuts.AddCut("80052013","1111151057032230000","0163103100000050"); // calo nonlinearity variation
    cuts.AddCut("80052013","1111152057032230000","0163103100000050"); // calo nonlinearity variation
  } else if(trainConfig == 51){ // second set of variations CLUSTER
    cuts.AddCut("80052013","1111141057022230000","0163103100000050"); // min energy cluster variation 1  600 MeV
    cuts.AddCut("80052013","1111141057042230000","0163103100000050"); // min energy cluster variation 2  800 MeV
    cuts.AddCut("80052013","1111141057052230000","0163103100000050"); // min energy cluster variation 3  900 MeV
  } else if(trainConfig == 52){ // third set of variations CLUSTER
    cuts.AddCut("80052013","1111141057032230000","0163103100000050"); // min/max M02  0.1<M<0.5
    cuts.AddCut("80052013","1111141057032200000","0163103100000050"); // min/max M02  0.1<M<100
    cuts.AddCut("80052013","1111141057032250000","0163103100000050"); // min/max M02  0.1<M<0.3
    cuts.AddCut("80052013","1111141057032260000","0163103100000050"); // min/max M02  0.1<M<0.27
    cuts.AddCut("80052013","1111141057031230000","0163103100000050"); // min number of cells variation 1  1 cell
    cuts.AddCut("80052013","1112141057032230000","0163103100000050"); // only modules with TRD infront
    cuts.AddCut("80052013","1111341057032230000","0163103100000050"); // no modules with TRD infront
  } else if(trainConfig == 53){ // third set of variations MESON
    cuts.AddCut("80052013","1111141057032230000","0163303100000050"); // rapidity variation  y<0.6
    cuts.AddCut("80052013","1111141057032230000","0163403100000050"); // rapidity variation  y<0.5
    cuts.AddCut("80052013","1111141057032230000","0163106100000050"); // alpha meson variation 1   0<alpha<0.8
    cuts.AddCut("80052013","1111141057032230000","0163105100000050"); // alpha meson variation 2  0<alpha<0.75
  } else if(trainConfig == 54){ // fourth set of variations
    cuts.AddCut("80052013","1111141053032230000","0163103100000050"); // tm variation
    cuts.AddCut("80052013","1111141055032230000","0163103100000050"); // tm variation
    cuts.AddCut("80052013","1111141053032230000","0163103100000040"); // min opening angle 0.75 cell diag
    cuts.AddCut("80052013","1111141053032230000","0163103100000050"); // min opening angle 1 cell diag
    cuts.AddCut("80052013","1111141053032230000","0163103100000060"); // min opening angle 2 cell diag
  } else if(trainConfig == 55){ // fourth set of variations
    cuts.AddCut("80052013","1111141057032230000","0163103100000050"); // default cut tm pt dependent
    cuts.AddCut("80052013","1111141056032230000","0163103100000050"); // tm pt dependent var 1
    cuts.AddCut("80052013","1111141058032230000","0163103100000050"); // tm pt dependent var 2
    cuts.AddCut("80052013","1111141059032230000","0163103100000050"); // tm pt dependent var 3
  } else if(trainConfig == 56){ // centrality dependent and NL variations
    cuts.AddCut("80252013","1111100057032230000","0163103100000050"); // 0-20% no NL
    cuts.AddCut("80252013","1111141057032230000","0163103100000050"); // 0-20% standard NL
    cuts.AddCut("82652013","1111100057032230000","0163103100000050"); // 20-60% no NL
    cuts.AddCut("82652013","1111141057032230000","0163103100000050"); // 20-60% standard NL
    cuts.AddCut("86052013","1111100057032230000","0163103100000050"); // 60-100% no NL
    cuts.AddCut("86052013","1111141057032230000","0163103100000050"); // 60-100% standard NL
  } else if(trainConfig == 57){ // new default cut
    cuts.AddCut("80052013","1111141057032230000","0163103100000050"); // default tm pt dependent

  //EG1
  } else if(trainConfig == 60){ // default cutstring and first set of variations nonlinearity
    cuts.AddCut("80083013","1111141057032230000","0163103100000050"); // default
    cuts.AddCut("80083013","1111142057032230000","0163103100000050"); // calo nonlinearity variation
    cuts.AddCut("80083013","1111151057032230000","0163103100000050"); // calo nonlinearity variation
    cuts.AddCut("80083013","1111152057032230000","0163103100000050"); // calo nonlinearity variation
  } else if(trainConfig == 61){ // second set of variations CLUSTER
    cuts.AddCut("80083013","1111141057022230000","0163103100000050"); // min energy cluster variation 1  600 MeV
    cuts.AddCut("80083013","1111141057042230000","0163103100000050"); // min energy cluster variation 2  800 MeV
    cuts.AddCut("80083013","1111141057052230000","0163103100000050"); // min energy cluster variation 3  900 MeV
  } else if(trainConfig == 62){ // third set of variations CLUSTER
    cuts.AddCut("80083013","1111141057032230000","0163103100000050"); // min/max M02  0.1<M<0.5
    cuts.AddCut("80083013","1111141057032200000","0163103100000050"); // min/max M02  0.1<M<100
    cuts.AddCut("80083013","1111141057032250000","0163103100000050"); // min/max M02  0.1<M<0.3
    cuts.AddCut("80083013","1111141057032260000","0163103100000050"); // min/max M02  0.1<M<0.27
    cuts.AddCut("80083013","1111141057031230000","0163103100000050"); // min number of cells variation 1  1 cell
    cuts.AddCut("80083013","1112141057032230000","0163103100000050"); // only modules with TRD infront
    cuts.AddCut("80083013","1111341057032230000","0163103100000050"); // no modules with TRD infront
  } else if(trainConfig == 63){ // third set of variations MESON
    cuts.AddCut("80083013","1111141057032230000","0163303100000050"); // rapidity variation  y<0.6
    cuts.AddCut("80083013","1111141057032230000","0163403100000050"); // rapidity variation  y<0.5
    cuts.AddCut("80083013","1111141057032230000","0163106100000050"); // alpha meson variation 1   0<alpha<0.8
    cuts.AddCut("80083013","1111141057032230000","0163105100000050"); // alpha meson variation 2  0<alpha<0.75
  } else if(trainConfig == 64){ // fourth set of variations
    cuts.AddCut("80083013","1111141053032230000","0163103100000050"); // tm variation
    cuts.AddCut("80083013","1111141055032230000","0163103100000050"); // tm variation
    cuts.AddCut("80083013","1111141053032230000","0163103100000040"); // min opening angle 0.75 cell diag
    cuts.AddCut("80083013","1111141053032230000","0163103100000050"); // min opening angle 1 cell diag
    cuts.AddCut("80083013","1111141053032230000","0163103100000060"); // min opening angle 2 cell diag
  } else if(trainConfig == 65){ // fourth set of variations
    cuts.AddCut("80083013","1111141057032230000","0163103100000050"); // default cut tm pt dependent
    cuts.AddCut("80083013","1111141056032230000","0163103100000050"); // tm pt dependent var 1
    cuts.AddCut("80083013","1111141058032230000","0163103100000050"); // tm pt dependent var 2
    cuts.AddCut("80083013","1111141059032230000","0163103100000050"); // tm pt dependent var 3
  } else if(trainConfig == 66){ // centrality dependent and NL variations
    cuts.AddCut("80283013","1111100057032230000","0163103100000050"); // 0-20% no NL
    cuts.AddCut("80283013","1111141057032230000","0163103100000050"); // 0-20% standard NL
    cuts.AddCut("82683013","1111100057032230000","0163103100000050"); // 20-60% no NL
    cuts.AddCut("82683013","1111141057032230000","0163103100000050"); // 20-60% standard NL
    cuts.AddCut("86083013","1111100057032230000","0163103100000050"); // 60-100% no NL
    cuts.AddCut("86083013","1111141057032230000","0163103100000050"); // 60-100% standard NL
  } else if(trainConfig == 67){ // new default cut
    cuts.AddCut("80083013","1111141057032230000","0163103100000050"); // default tm pt dependent
   //EG2
  } else if(trainConfig == 70){ // default cutstring and first set of variations nonlinearity
    cuts.AddCut("80085013","1111141057032230000","0163103100000050"); // default
    cuts.AddCut("80085013","1111142057032230000","0163103100000050"); // calo nonlinearity variation
    cuts.AddCut("80085013","1111151057032230000","0163103100000050"); // calo nonlinearity variation
    cuts.AddCut("80085013","1111152057032230000","0163103100000050"); // calo nonlinearity variation
  } else if(trainConfig == 71){ // second set of variations CLUSTER
    cuts.AddCut("80085013","1111141057022230000","0163103100000050"); // min energy cluster variation 1  600 MeV
    cuts.AddCut("80085013","1111141057042230000","0163103100000050"); // min energy cluster variation 2  800 MeV
    cuts.AddCut("80085013","1111141057052230000","0163103100000050"); // min energy cluster variation 3  900 MeV
  } else if(trainConfig == 72){ // third set of variations CLUSTER
    cuts.AddCut("80085013","1111141057032230000","0163103100000050"); // min/max M02  0.1<M<0.5
    cuts.AddCut("80085013","1111141057032200000","0163103100000050"); // min/max M02  0.1<M<100
    cuts.AddCut("80085013","1111141057032250000","0163103100000050"); // min/max M02  0.1<M<0.3
    cuts.AddCut("80085013","1111141057032260000","0163103100000050"); // min/max M02  0.1<M<0.27
    cuts.AddCut("80085013","1111141057031230000","0163103100000050"); // min number of cells variation 1  1 cell
    cuts.AddCut("80085013","1112141057032230000","0163103100000050"); // only modules with TRD infront
    cuts.AddCut("80085013","1111341057032230000","0163103100000050"); // no modules with TRD infront
  } else if(trainConfig == 73){ // third set of variations MESON
    cuts.AddCut("80085013","1111141057032230000","0163303100000050"); // rapidity variation  y<0.6
    cuts.AddCut("80085013","1111141057032230000","0163403100000050"); // rapidity variation  y<0.5
    cuts.AddCut("80085013","1111141057032230000","0163106100000050"); // alpha meson variation 1   0<alpha<0.8
    cuts.AddCut("80085013","1111141057032230000","0163105100000050"); // alpha meson variation 2  0<alpha<0.75
  } else if(trainConfig == 74){ // fourth set of variations
    cuts.AddCut("80085013","1111141053032230000","0163103100000050"); // tm variation
    cuts.AddCut("80085013","1111141055032230000","0163103100000050"); // tm variation
    cuts.AddCut("80085013","1111141053032230000","0163103100000040"); // min opening angle 0.75 cell diag
    cuts.AddCut("80085013","1111141053032230000","0163103100000050"); // min opening angle 1 cell diag
    cuts.AddCut("80085013","1111141053032230000","0163103100000060"); // min opening angle 2 cell diag
  } else if(trainConfig == 75){ // fourth set of variations
    cuts.AddCut("80085013","1111141057032230000","0163103100000050"); // default cut tm pt dependent
    cuts.AddCut("80085013","1111141056032230000","0163103100000050"); // tm pt dependent var 1
    cuts.AddCut("80085013","1111141058032230000","0163103100000050"); // tm pt dependent var 2
    cuts.AddCut("80085013","1111141059032230000","0163103100000050"); // tm pt dependent var 3
  } else if(trainConfig == 76){ // centrality dependent and NL variations
    cuts.AddCut("80285013","1111100057032230000","0163103100000050"); // 0-20% no NL
    cuts.AddCut("80285013","1111141057032230000","0163103100000050"); // 0-20% standard NL
    cuts.AddCut("82685013","1111100057032230000","0163103100000050"); // 20-60% no NL
    cuts.AddCut("82685013","1111141057032230000","0163103100000050"); // 20-60% standard NL
    cuts.AddCut("86085013","1111100057032230000","0163103100000050"); // 60-100% no NL
    cuts.AddCut("86085013","1111141057032230000","0163103100000050"); // 60-100% standard NL
  } else if(trainConfig == 77){ // new default cut
    cuts.AddCut("80085013","1111141057032230000","0163103100000050"); // default tm pt dependent

  //all default triggers
  } else if(trainConfig == 80){
    cuts.AddCut("80000013","1111141057032230000","0163103100000050"); // default MB
    cuts.AddCut("80052013","1111141057032230000","0163103100000050"); // default EMC7
    cuts.AddCut("80083013","1111141057032230000","0163103100000050"); // default EG1
    cuts.AddCut("80085013","1111141057032230000","0163103100000050"); // default EG2
  //testing different event mix methods
  } else if(trainConfig == 81){
    cuts.AddCut("80000013","1111141057032230000","0163103100000050"); // default (V0 mult)
    cuts.AddCut("80000013","1111141057032230000","0263103100000050"); // using track mult
    cuts.AddCut("80000013","1111141057032230000","0963103100000050"); // using PtMax method DeltaR < 0.2
    
  //************************************************ PHOS clusters *************************************************
  } else if (trainConfig == 301) {  // min energy = 0.3 GeV/c
    cuts.AddCut("80000013","2444400041013200000","0163103100000010"); //standart cut, kINT7 // PHOS clusters
    cuts.AddCut("80062013","2444400041013200000","0163103100000010"); //standard cut, kPHI7  // PHOS clusters
  } else if (trainConfig == 302){ // Validation PHOS
    cuts.AddCut("80000013","2444400041013200000","0163103100000010");
  } else if (trainConfig == 303){ // Validation PHOS, only added signals
    cuts.AddCut("80000023","2444400041013200000","0163103100000010");
  } else if (trainConfig == 304){ // min energy = 0.3 GeV/c
    cuts.AddCut("80000013","2444400041013200000","0163103100000000"); // kINT7 // PHOS clusters no open angle cut
    cuts.AddCut("80062013","2444400041013200000","0163103100000000"); // kPHI7 // PHOS clusters no open angle cut
    cuts.AddCut("80000013","2444400041013200000","0163103100000030"); // kINT7 // PHOS clusters open angle cut 0.1
    cuts.AddCut("80062013","2444400041013200000","0163103100000030"); // kPHI7 // PHOS clusters  open angle cut 0.1
  } else if (trainConfig == 305){ // timing cut variations
    cuts.AddCut("80000013","2444400011013200000","0163103100000010"); // 1000ns
    cuts.AddCut("80000013","2444400031013200000","0163103100000010"); // 200ns
    cuts.AddCut("80000013","2444400051013200000","0163103100000010"); // 50ns

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

  //connect containers
  AliAnalysisDataContainer *coutput =
    mgr->CreateContainer(Form("GammaCalo_%i",trainConfig), TList::Class(),
          AliAnalysisManager::kOutputContainer,Form("GammaCalo_%i.root",trainConfig));

  mgr->AddTask(task);
  mgr->ConnectInput(task,0,cinput);
  mgr->ConnectOutput(task,1,coutput);

  return;

}
