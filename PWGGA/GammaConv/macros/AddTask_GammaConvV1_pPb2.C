/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: Friederike Bock, Annika Passfeld                               *
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
//($ALIPHYSICS/PWGGA/GammaConv/AliAnalysisTaskGammaConvV1.cxx) for
//pp together with all supporting classes
//***************************************************************************************

//***************************************************************************************
//CutHandler contains all cuts for a certain analysis and trainconfig,
//it automatically checks length of cutStrings and takes care of the number of added cuts,
//no specification of the variable 'numberOfCuts' needed anymore.
//***************************************************************************************

class CutHandlerConv{
  public:
    CutHandlerConv(Int_t nMax=10){
      nCuts=0; nMaxCuts=nMax; validCuts = true;
      eventCutArray = new TString[nMaxCuts]; photonCutArray = new TString[nMaxCuts]; mesonCutArray = new TString[nMaxCuts];
      for(Int_t i=0; i<nMaxCuts; i++) {eventCutArray[i] = ""; photonCutArray[i] = ""; mesonCutArray[i] = "";}
    }

    void AddCut(TString eventCut, TString photonCut, TString mesonCut){
      if(nCuts>=nMaxCuts) {cout << "ERROR in CutHandlerConv: Exceeded maximum number of cuts!" << endl; validCuts = false; return;}
      if( eventCut.Length()!=8 || photonCut.Length()!=26 || mesonCut.Length()!=16 ) {cout << "ERROR in CutHandlerConv: Incorrect length of cut string!" << endl; validCuts = false; return;}
      eventCutArray[nCuts]=eventCut; photonCutArray[nCuts]=photonCut; mesonCutArray[nCuts]=mesonCut;
      nCuts++;
      return;
    }
    Bool_t AreValid(){return validCuts;}
    Int_t GetNCuts(){if(validCuts) return nCuts; else return 0;}
    TString GetEventCut(Int_t i){if(validCuts&&i<nMaxCuts&&i>=0) return eventCutArray[i]; else{cout << "ERROR in CutHandlerConv: GetEventCut wrong index i" << endl;return "";}}
    TString GetPhotonCut(Int_t i){if(validCuts&&i<nMaxCuts&&i>=0) return photonCutArray[i]; else {cout << "ERROR in CutHandlerConv: GetPhotonCut wrong index i" << endl;return "";}}
    TString GetMesonCut(Int_t i){if(validCuts&&i<nMaxCuts&&i>=0) return mesonCutArray[i]; else {cout << "ERROR in CutHandlerConv: GetMesonCut wrong index i" << endl;return "";}}
  private:
    Bool_t validCuts;
    Int_t nCuts; Int_t nMaxCuts;
    TString* eventCutArray;
    TString* photonCutArray;
    TString* mesonCutArray;
};

void AddTask_GammaConvV1_pPb2(  Int_t         trainConfig                   = 1,                                // change different set of cuts
                                Int_t         isMC                          = 0,                                // run MC
                                Int_t         enableQAMesonTask             = 0,                                // enable QA in AliAnalysisTaskGammaConvV1
                                Int_t         enableQAPhotonTask            = 0,                                // enable additional QA task
                                TString       fileNameInputForWeighting     = "MCSpectraInput.root",            // path to file for weigting input
                                Int_t         doWeightingPart               = 0,                                // enable Weighting
                                TString       generatorName                 = "DPMJET",                         // generator Name
                                TString       cutnumberAODBranch            = "800000006008400000001500000",    // cutnumber for AOD branch
                                Bool_t        enableV0findingEffi           = kFALSE,                           // enables V0finding efficiency histograms
                                Bool_t        enablePlotVsCentrality        = kFALSE,
                                Bool_t        enableTriggerMimicking        = kFALSE,                           // enable trigger mimicking
                                Bool_t        enableTriggerOverlapRej       = kFALSE,                           // enable trigger overlap rejection
                                Float_t       maxFacPtHard                  = 3.,                               // maximum factor between hardest jet and ptHard generated
                                TString       periodNameV0Reader            = "",
                                Bool_t        doMultiplicityWeighting       = kFALSE,                           //
                                TString       fileNameInputForMultWeighing  = "Multiplicity.root",              //
                                TString       periodNameAnchor              = "",
                                Bool_t        runTHnSparse                  = kTRUE,                            // switch on THNsparse
                                Bool_t        runLightOutput                = kFALSE,                           // switch to run light output (only essential histograms for afterburner)
                                TString       additionalTrainConfig         = "0"                               // additional counter for trainconfig, this has to be always the last parameter
                          ) {

  Int_t isHeavyIon = 2;
  if (additionalTrainConfig.Atoi() > 0){
    trainConfig = trainConfig + additionalTrainConfig.Atoi();
  }

  // ================== GetAnalysisManager ===============================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error(Form("AddTask_GammaConvV1_%i",trainConfig), "No analysis manager found.");
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

  //=========  Set Cutnumber for V0Reader ================================
  TString cutnumberPhoton = "00000008400100001500000000";
  TString cutnumberEvent = "80000003";
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
      if(trainConfig==15 ||trainConfig==16 ||trainConfig==17  ||trainConfig==18  ){
        fCuts->SetDodEdxSigmaCut(kFALSE);
      }
    }
    if(inputHandler->IsA()==AliAODInputHandler::Class()){
    // AOD mode
      cout << "AOD handler: adding " << cutnumberAODBranch.Data() << " as conversion branch" << endl;
      fV0ReaderV1->AliV0ReaderV1::SetDeltaAODBranchName(Form("GammaConv_%s_gamma",cutnumberAODBranch.Data()));
    }
    fV0ReaderV1->Init();

    AliLog::SetGlobalLogLevel(AliLog::kInfo);

    //connect input V0Reader
    mgr->AddTask(fV0ReaderV1);
    mgr->ConnectInput(fV0ReaderV1,0,cinput);

  }

  //================================================
  //========= Add task to the ANALYSIS manager =====
  //================================================
  //            find input container
  AliAnalysisTaskGammaConvV1 *task=NULL;
  task= new AliAnalysisTaskGammaConvV1(Form("GammaConvV1_%i",trainConfig));
  task->SetIsHeavyIon(isHeavyIon);
  task->SetIsMC(isMC);
  task->SetV0ReaderName(V0ReaderName);
  task->SetLightOutput(runLightOutput);
  // Cut Numbers to use in Analysis

  CutHandlerConv cuts;

  Bool_t doEtaShiftIndCuts = kFALSE;
  TString stringShift = "";

  // standard configurations
  if(trainConfig == 1){
    cuts.AddCut("80010113", "00200009327000008250404000", "0162103500900000"); // new standard pPb MB
  } else if (trainConfig == 2) {
    cuts.AddCut("80010113", "00200009327000008250400000", "0162103500900000"); // new standard pPb MB
  } else if (trainConfig == 3) {
    cuts.AddCut("80210113", "00200009327000008250404000", "0162103500900000"); // new standard pPb 0-20
  } else if (trainConfig == 4) {
    cuts.AddCut("80210113", "00200009327000008250400000", "0162103500900000"); // new standard pPb 0-20
  } else if (trainConfig == 5) {
    cuts.AddCut("82410113", "00200009327000008250404000", "0162103500900000"); // new standard pPb 20-40
  } else if (trainConfig == 6) {
    cuts.AddCut("82410113", "00200009327000008250400000", "0162103500900000"); // new standard pPb 20-40
  } else if (trainConfig == 7) {
    cuts.AddCut("84610113", "00200009327000008250404000", "0162103500900000"); // new standard pPb 40-60
  } else if (trainConfig == 8) {
    cuts.AddCut("84610113", "00200009327000008250400000", "0162103500900000"); // new standard pPb 40-60
  } else if (trainConfig == 9) {
    cuts.AddCut("86010113", "00200009327000008250404000", "0162103500900000"); // new standard pPb 60-80
  } else if (trainConfig == 10) {
    cuts.AddCut("86010113", "00200009327000008250400000", "0162103500900000"); // new standard pPb 60-80

  // configurations with past future protection (2.25 \mus protected)
  } else if (trainConfig == 11){
    cuts.AddCut("80010213", "00200009327000008250404000", "0162103500900000"); // new standard pPb MB
  } else if (trainConfig == 12) {
    cuts.AddCut("80010213", "00200009327000008250400000", "0162103500900000"); // new standard pPb MB
  } else if (trainConfig == 13) {
    cuts.AddCut("80210213", "00200009327000008250404000", "0162103500900000"); // new standard pPb 0-20
  } else if (trainConfig == 14) {
    cuts.AddCut("80210213", "00200009327000008250400000", "0162103500900000"); // new standard pPb 0-20
  } else if (trainConfig == 15) {
    cuts.AddCut("82410213", "00200009327000008250404000", "0162103500900000"); // new standard pPb 20-40
  } else if (trainConfig == 16) {
    cuts.AddCut("82410213", "00200009327000008250400000", "0162103500900000"); // new standard pPb 20-40
  } else if (trainConfig == 17) {
    cuts.AddCut("84610213", "00200009327000008250404000", "0162103500900000"); // new standard pPb 40-60
  } else if (trainConfig == 18) {
    cuts.AddCut("84610213", "00200009327000008250400000", "0162103500900000"); // new standard pPb 40-60
  } else if (trainConfig == 19) {
    cuts.AddCut("86010213", "00200009327000008250404000", "0162103500900000"); // new standard pPb 60-80
  } else if (trainConfig == 20) {
    cuts.AddCut("86010213", "00200009327000008250400000", "0162103500900000"); // new standard pPb 60-80

  // configurations for eta cuts
  } else if (trainConfig == 21) {
    cuts.AddCut("80010113", "0a200009327000008250400000", "0162103500000000"); //Eta cut -0.9 - -0.2 and 0.2 - 0.9
  } else if (trainConfig == 22) {
    cuts.AddCut("80010113", "0b200009327000008250400000", "0162103500000000"); //Eta cut -0.9 - -0.2 and 0.2 - 0.9 with LineCut

    // configurations for addsig (addsig part) with weighting
  } else if (trainConfig == 30) {
    cuts.AddCut("80010113", "00200009327000008250400000", "0162103500900000"); // new standard pPb MB
  } else if (trainConfig == 31) {
    cuts.AddCut("80010113", "00200009327000008250400000", "0162103500900000"); // new standard pPb MB
  } else if (trainConfig == 32) {
    cuts.AddCut("80210113", "00200009327000008250400000", "0162103500900000"); // new standard pPb 0-20
  } else if (trainConfig == 33) {
    cuts.AddCut("80210113", "00200009327000008250400000", "0162103500900000"); // new standard pPb 0-20
  } else if (trainConfig == 34) {
    cuts.AddCut("82410113", "00200009327000008250400000", "0162103500900000"); // new standard pPb 20-40
  } else if (trainConfig == 35) {
    cuts.AddCut("82410113", "00200009327000008250400000", "0162103500900000"); // new standard pPb 20-40
  } else if (trainConfig == 36) {
    cuts.AddCut("84610113", "00200009327000008250400000", "0162103500900000"); // new standard pPb 40-60
  } else if (trainConfig == 37) {
    cuts.AddCut("84610113", "00200009327000008250400000", "0162103500900000"); // new standard pPb 40-60
  } else if (trainConfig == 38) {
    cuts.AddCut("86010113", "00200009327000008250400000", "0162103500900000"); // new standard pPb 60-80
  } else if (trainConfig == 39) {
    cuts.AddCut("86010113", "00200009327000008250400000", "0162103500900000"); // new standard pPb 60-80

    // configurations for addsig (MB part) with weighting
  } else if (trainConfig == 40) {
    cuts.AddCut("80010113", "00200009327000008250400000", "0162103500900000"); // new standard pPb MB
  } else if (trainConfig == 41) {
    cuts.AddCut("80010113", "00200009327000008250400000", "0162103500900000"); // new standard pPb MB
  } else if (trainConfig == 42) {
    cuts.AddCut("80210113", "00200009327000008250400000", "0162103500900000"); // new standard pPb 0-20
  } else if (trainConfig == 43) {
    cuts.AddCut("80210113", "00200009327000008250400000", "0162103500900000"); // new standard pPb 0-20
  } else if (trainConfig == 44) {
    cuts.AddCut("82410113", "00200009327000008250400000", "0162103500900000"); // new standard pPb 20-40
  } else if (trainConfig == 45) {
    cuts.AddCut("82410113", "00200009327000008250400000", "0162103500900000"); // new standard pPb 20-40
  } else if (trainConfig == 46) {
    cuts.AddCut("84610113", "00200009327000008250400000", "0162103500900000"); // new standard pPb 40-60
  } else if (trainConfig == 47) {
    cuts.AddCut("84610113", "00200009327000008250400000", "0162103500900000"); // new standard pPb 40-60
  } else if (trainConfig == 48) {
    cuts.AddCut("86010113", "00200009327000008250400000", "0162103500900000"); // new standard pPb 60-80
  } else if (trainConfig == 49) {
    cuts.AddCut("86010113", "00200009327000008250400000", "0162103500900000"); // new standard pPb 60-80

  //Run 2 pPb
  } else if (trainConfig == 100) {
    cuts.AddCut("80010113", "00200009327000008250404000", "0162103500000000"); // new default for 5TeV
  } else if (trainConfig == 101) {
    cuts.AddCut("80110113", "00200009327000008250404000", "0162103500000000"); // 0-10
  } else if (trainConfig == 102) {
    cuts.AddCut("81210113", "00200009327000008250404000", "0162103500000000"); // 0-20
  } else if (trainConfig == 103) {
    cuts.AddCut("82410113", "00200009327000008250404000", "0162103500000000"); // 20-40
  } else if (trainConfig == 104) {
    cuts.AddCut("84610113", "00200009327000008250404000", "0162103500000000"); // 40-60
  } else if (trainConfig == 105) {
    cuts.AddCut("86810113", "00200009327000008250404000", "0162103500000000"); // 60-80
  } else if (trainConfig == 106) {
    cuts.AddCut("88010113", "00200009327000008250404000", "0162103500000000"); // 80-100
  } else if (trainConfig == 107) {
    cuts.AddCut("80210113", "00200009327000008250404000", "0162103500000000"); // 0-20
  } else if (trainConfig == 108) {
    cuts.AddCut("86010113", "00200009327000008250404000", "0162103500000000"); // 60-100
  } else if (trainConfig == 109) {
    cuts.AddCut("a0110113", "00200009327000008250404000", "0162103500000000"); // 0-5
  } else if (trainConfig == 110) {
    cuts.AddCut("a1210113", "00200009327000008250404000", "0162103500000000"); // 5-10
  } else if (trainConfig == 111) {
    cuts.AddCut("c0110113", "00200009327000008250404000", "0162103500000000"); // 0-1
  } else if (trainConfig == 112) {
    cuts.AddCut("c0210113", "00200009327000008250404000", "0162103500000000"); // 0-2
  } else {
    Error(Form("GammaConvV1_%i",trainConfig), "wrong trainConfig variable no cuts have been specified for the configuration");
    return;
  }
  if(!cuts.AreValid()){
    cout << "\n\n****************************************************" << endl;
    cout << "ERROR: No valid cuts stored in CutHandlerConv! Returning..." << endl;
    cout << "****************************************************\n\n" << endl;
    return;
  }

  Int_t numberOfCuts = cuts.GetNCuts();

  TList *EventCutList = new TList();
  TList *ConvCutList = new TList();
  TList *MesonCutList = new TList();

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
  Bool_t doWeighting = kFALSE;
  if (doWeightingPart == 1 || doWeightingPart == 2 || doWeightingPart == 3) doWeighting = kTRUE;

  EventCutList->SetOwner(kTRUE);
  AliConvEventCuts **analysisEventCuts = new AliConvEventCuts*[numberOfCuts];
  ConvCutList->SetOwner(kTRUE);
  AliConversionPhotonCuts **analysisCuts = new AliConversionPhotonCuts*[numberOfCuts];
  MesonCutList->SetOwner(kTRUE);
  AliConversionMesonCuts **analysisMesonCuts = new AliConversionMesonCuts*[numberOfCuts];

  if (doWeighting) Printf("weighting has been switched on");

  for(Int_t i = 0; i<numberOfCuts; i++){

    analysisEventCuts[i] = new AliConvEventCuts();
    if ( trainConfig == 13 || trainConfig == 15 || trainConfig == 17 || trainConfig == 19 || trainConfig == 20 || trainConfig == 22 || ( trainconfig > 39 && trainconfig < 50 )){
      if (doWeighting){
        if (generatorName.CompareTo("DPMJET")==0){
          analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE, fileNameInputForWeighting, "Pi0_DPMJET_LHC13b2_efix_pPb_5023GeV_MBV0A",
                                        "Eta_DPMJET_LHC13b2_efix_pPb_5023GeV_MBV0A", "","Pi0_Fit_Data_pPb_5023GeV_MBV0A","Eta_Fit_Data_pPb_5023GeV_MBV0A");
        } else if (generatorName.CompareTo("HIJING")==0){
          analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13e7_pPb_5023GeV_MBV0A",
                                        "Eta_Hijing_LHC13e7_pPb_5023GeV_MBV0A", "","Pi0_Fit_Data_pPb_5023GeV_MBV0A","Eta_Fit_Data_pPb_5023GeV_MBV0A");
        }
      }
    }
    if ( trainConfig == 14 || trainConfig == 16 || trainConfig == 18 || trainConfig == 21 || trainConfig == 23 || ( trainconfig > 29 && trainconfig < 40 ) ){
      if (doWeighting){
        analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13e7_addSig_pPb_5023GeV_MBV0A",
                                      "Eta_Hijing_LHC13e7_addSig_pPb_5023GeV_MBV0A", "","Pi0_Fit_Data_pPb_5023GeV_MBV0A","Eta_Fit_Data_pPb_5023GeV_MBV0A");
      }

    }


    TString dataInputMultHisto  = "";
    TString mcInputMultHisto    = "";
    TString triggerString   = (cuts.GetEventCut(i)).Data();
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
    analysisEventCuts[i]->SetV0ReaderName(V0ReaderName);
    if (periodNameV0Reader.CompareTo("") != 0) analysisEventCuts[i]->SetPeriodEnum(periodNameV0Reader);
    analysisEventCuts[i]->SetLightOutput(runLightOutput);
    analysisEventCuts[i]->InitializeCutsFromCutString((cuts.GetEventCut(i)).Data());
    if (doEtaShiftIndCuts) {
      analysisEventCuts[i]->DoEtaShift(doEtaShiftIndCuts);
      analysisEventCuts[i]->SetEtaShift(stringShift);
    }

    EventCutList->Add(analysisEventCuts[i]);
    analysisEventCuts[i]->SetFillCutHistograms("",kFALSE);

    analysisCuts[i] = new AliConversionPhotonCuts();
    analysisCuts[i]->SetV0ReaderName(V0ReaderName);
    analysisCuts[i]->SetLightOutput(runLightOutput);
    analysisCuts[i]->InitializeCutsFromCutString((cuts.GetPhotonCut(i)).Data());
    analysisCuts[i]->SetIsHeavyIon(isHeavyIon);
    if (trainConfig == 15 || trainConfig==16 || trainConfig==17  || trainConfig==18) {
            analysisCuts[i]->SetDodEdxSigmaCut(kFALSE);
    }
    ConvCutList->Add(analysisCuts[i]);
    analysisCuts[i]->SetFillCutHistograms("",kFALSE);

    analysisMesonCuts[i] = new AliConversionMesonCuts();
    if (trainConfig ==13 || trainConfig ==14){
      analysisMesonCuts[i]->SetOpeningAngleCut(0.000);
    }
    analysisMesonCuts[i]->SetLightOutput(runLightOutput);
    analysisMesonCuts[i]->InitializeCutsFromCutString((cuts.GetMesonCut(i)).Data());
    MesonCutList->Add(analysisMesonCuts[i]);
    analysisMesonCuts[i]->SetFillCutHistograms("");
    analysisEventCuts[i]->SetAcceptedHeader(HeaderList);
  }

  task->SetDoTHnSparse(runTHnSparse);
  task->SetEventCutList(numberOfCuts,EventCutList);
  task->SetConversionCutList(numberOfCuts,ConvCutList);
  task->SetMesonCutList(numberOfCuts,MesonCutList);
  task->SetMoveParticleAccordingToVertex(kTRUE);
  task->SetDoMesonAnalysis(kTRUE);
  task->SetDoMesonQA(enableQAMesonTask); //Attention new switch for Pi0 QA
  task->SetDoPhotonQA(enableQAPhotonTask);  //Attention new switch small for Photon QA
  if (trainConfig ==13 || trainConfig ==14){
          task->SetDoTHnSparse(0);
  }
  task->SetDoPlotVsCentrality(enablePlotVsCentrality);

  //connect containers
  AliAnalysisDataContainer *coutput =
    mgr->CreateContainer(Form("GammaConvV1_%i",trainConfig), TList::Class(),
              AliAnalysisManager::kOutputContainer,Form("GammaConvV1_%i.root",trainConfig));

  mgr->AddTask(task);
  mgr->ConnectInput(task,0,cinput);
  mgr->ConnectOutput(task,1,coutput);

  return;
}
