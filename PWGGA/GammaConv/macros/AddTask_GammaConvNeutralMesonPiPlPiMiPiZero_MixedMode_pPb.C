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
//($ALIPHYSICS/PWGGA/GammaConv/AliAnalysisTaskNeutralMesonToPiPlPiMiPiZero.cxx) for
//pp together with all supporting classes
//***************************************************************************************

//***************************************************************************************
//CutHandler contains all cuts for a certain analysis and trainconfig,
//it automatically checks length of cutStrings and takes care of the number of added cuts,
//no specification of the variable 'numberOfCuts' needed anymore.
//***************************************************************************************
class CutHandlerNeutralMixed{
  public:
    CutHandlerNeutralMixed(Int_t nMax=10){
      nCuts=0; nMaxCuts=nMax; validCuts = true;
      eventCutArray = new TString[nMaxCuts]; clusterCutArray = new TString[nMaxCuts]; conversionCutArray = new TString[nMaxCuts]; pionCutArray = new TString[nMaxCuts]; neutralPionCutArray = new TString[nMaxCuts]; mesonCutArray = new TString[nMaxCuts];
      for(Int_t i=0; i<nMaxCuts; i++) {eventCutArray[i] = ""; clusterCutArray[i] = ""; conversionCutArray[i] = ""; pionCutArray[i] = ""; neutralPionCutArray[i] = ""; mesonCutArray[i] = "";}
    }

    void AddCut(TString eventCut, TString clusterCut, TString conversionCut, TString pionCut, TString neutralPionCut, TString mesonCut){
      if(nCuts>=nMaxCuts) {cout << "ERROR in CutHandlerNeutralMixed: Exceeded maximum number of cuts!" << endl; validCuts = false; return;}
      if( eventCut.Length()!=8 || clusterCut.Length()!=19 || conversionCut.Length()!=26 || pionCut.Length()!=9 || neutralPionCut.Length()!=16 || mesonCut.Length()!=16 ) {cout << "ERROR in CutHandlerNeutralMixed: Incorrect length of cut string!" << endl; validCuts = false; return;}
      eventCutArray[nCuts]=eventCut; clusterCutArray[nCuts]=clusterCut; conversionCutArray[nCuts]=conversionCut; pionCutArray[nCuts]=pionCut; neutralPionCutArray[nCuts]=neutralPionCut; mesonCutArray[nCuts]=mesonCut;
      nCuts++;
      return;
    }
    Bool_t AreValid(){return validCuts;}
    Int_t GetNCuts(){if(validCuts) return nCuts; else return 0;}
    TString GetEventCut(Int_t i){if(validCuts&&i<nMaxCuts&&i>=0) return eventCutArray[i]; else{cout << "ERROR in CutHandlerNeutralMixed: GetEventCut wrong index i" << endl;return "";}}
    TString GetClusterCut(Int_t i){if(validCuts&&i<nMaxCuts&&i>=0) return clusterCutArray[i]; else {cout << "ERROR in CutHandlerNeutralMixed: GetClusterCut wrong index i" << endl;return "";}}
    TString GetConversionCut(Int_t i){if(validCuts&&i<nMaxCuts&&i>=0) return conversionCutArray[i]; else {cout << "ERROR in CutHandlerNeutralMixed: GetConversionCut wrong index i" << endl;return "";}}
    TString GetPionCut(Int_t i){if(validCuts&&i<nMaxCuts&&i>=0) return pionCutArray[i]; else {cout << "ERROR in CutHandlerNeutralMixed: GetPionCut wrong index i" << endl;return "";}}
    TString GetNeutralPionCut(Int_t i){if(validCuts&&i<nMaxCuts&&i>=0) return neutralPionCutArray[i]; else {cout << "ERROR in CutHandlerNeutralMixed: GetNeutralPionCut wrong index i" << endl;return "";}}
    TString GetMesonCut(Int_t i){if(validCuts&&i<nMaxCuts&&i>=0) return mesonCutArray[i]; else {cout << "ERROR in CutHandlerNeutralMixed: GetMesonCut wrong index i" << endl;return "";}}
  private:
    Bool_t validCuts;
    Int_t nCuts; Int_t nMaxCuts;
    TString* eventCutArray;
    TString* clusterCutArray;
    TString* conversionCutArray;
    TString* pionCutArray;
    TString* neutralPionCutArray;
    TString* mesonCutArray;
};

//***************************************************************************************
//main function
//***************************************************************************************
void AddTask_GammaConvNeutralMesonPiPlPiMiPiZero_MixedMode_pPb(    
    Int_t trainConfig                 = 1,
    Bool_t isMC                       = kFALSE,                         //run MC
    Bool_t enableQAMesonTask          = kTRUE,                          //enable QA in AliAnalysisTaskNeutralMesonToPiPlPiMiPiZero
    TString fileNameInputForWeighting = "MCSpectraInput.root",          // path to file for weigting input
    Bool_t doWeighting                = kFALSE,  //enable Weighting
    TString generatorName             = "HIJING",
    TString cutnumberAODBranch        = "000000006008400001001500000",
    Double_t tolerance                = -1,
    TString additionalTrainConfig     = "0"                              // additional counter for trainconfig, this has to be always the last parameter
  ) {

  //parse additionalTrainConfig flag
  TObjArray *rAddConfigArr = additionalTrainConfig.Tokenize("_");
  if(rAddConfigArr->GetEntries()<1){cout << "ERROR during parsing of additionalTrainConfig String '" << additionalTrainConfig.Data() << "'" << endl; return;}
  TObjString* rAdditionalTrainConfig;
  for(Int_t i = 0; i<rAddConfigArr->GetEntries() ; i++){
    if(i==0) rAdditionalTrainConfig = (TObjString*)rAddConfigArr->At(i);
    else{
      TObjString* temp = (TObjString*) rAddConfigArr->At(i);
      TString tempStr = temp->GetString();
      cout << "INFO: nothing to do, no definition available!" << endl;
    }
  }
  TString sAdditionalTrainConfig = rAdditionalTrainConfig->GetString();
  if (sAdditionalTrainConfig.Atoi() > 0){
    trainConfig = trainConfig + sAdditionalTrainConfig.Atoi();
    cout << "INFO: AddTask_GammaConvNeutralMesonPiPlPiMiPiZero_MixedMode_pPb running additionalTrainConfig '" << sAdditionalTrainConfig.Atoi() << "', train config: '" << trainConfig << "'" << endl;
  }

  Int_t isHeavyIon = 2;
  Int_t neutralPionMode = 1;
  
  // ================== GetAnalysisManager ===============================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error(Form("AddTask_GammaConvNeutralMesonPiPlPiMiPiZero_MixedMode_pPb_%i",trainConfig), "No analysis manager found.");
    return ;
  }

  // ================== GetInputEventHandler =============================
  AliVEventHandler *inputHandler=mgr->GetInputEventHandler();
  
  //========= Add PID Reponse to ANALYSIS manager ====
  if(!(AliPIDResponse*)mgr->GetTask("PIDResponseTask")){
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
    AddTaskPIDResponse(isMC);
  }
  
  //=========  Set Cutnumber for V0Reader ================================
  TString cutnumberPhoton = "06000008400100001500000000";
  TString cutnumberEvent = "80000003";
  TString PionCuts      = "000000200";            //Electron Cuts
  

  
  Bool_t doEtaShift = kFALSE;

  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  
  //========= Add V0 Reader to  ANALYSIS manager if not yet existent =====
    TString V0ReaderName = Form("V0ReaderV1_%s_%s",cutnumberEvent.Data(),cutnumberPhoton.Data());
    if( !(AliV0ReaderV1*)mgr->GetTask(V0ReaderName.Data()) ){
        AliV0ReaderV1 *fV0ReaderV1 = new AliV0ReaderV1(V0ReaderName.Data());
    
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
  //========= Add Electron Selector ================


  if( !(AliPrimaryPionSelector*)mgr->GetTask("PionSelector") ){

    AliPrimaryPionSelector *fPionSelector = new AliPrimaryPionSelector("PionSelector");
    // Set AnalysisCut Number

    AliPrimaryPionCuts *fPionCuts=0;
    if( PionCuts!=""){
      fPionCuts= new AliPrimaryPionCuts(PionCuts.Data(),PionCuts.Data());
      if(fPionCuts->InitializeCutsFromCutString(PionCuts.Data())){
        fPionSelector->SetPrimaryPionCuts(fPionCuts);
        fPionCuts->SetFillCutHistograms("",kTRUE);

      }
    }

    fPionSelector->Init();
    mgr->AddTask(fPionSelector);
    
    AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer();

    //connect input V0Reader
    mgr->ConnectInput (fPionSelector,0,cinput1);

  }

  
  
  AliAnalysisTaskNeutralMesonToPiPlPiMiPiZero *task=NULL;
  task= new AliAnalysisTaskNeutralMesonToPiPlPiMiPiZero(Form("GammaConvNeutralMesonPiPlPiMiPiZero_%i_%i",neutralPionMode, trainConfig));
  task->SetIsHeavyIon(2);
  task->SetIsMC(isMC);
  task->SetV0ReaderName(V0ReaderName);
  task->SetTolerance(tolerance);

  CutHandlerNeutralMixed cuts;
  
  Bool_t doEtaShiftIndCuts = kFALSE;
  TString stringShift = "";

  // Shifting in pPb direction

  doEtaShiftIndCuts = kTRUE;
  stringShift = "pPb";

  // EMCAL mode
  if( trainConfig == 1 ) {
    // everything open
    cuts.AddCut("80000113","00200009117000008260400000","1111111067032230000","000010400","0103503000000000","0103503000000000");
  } else if( trainConfig == 2 ) {
    // closing charged pion cuts, minimum TPC cluster = 80, TPC dEdx pi = \pm 3 sigma, min pt charged pi = 100 MeV
    cuts.AddCut("80000113","00200009117000008260400000","1111111067032230000","002010700","0103503000000000","0103503000000000");
  } else if( trainConfig == 3 ) {
    // closing charged pion cuts, minimum TPC cluster = 80, TPC dEdx pi = \pm 3 sigma, ITS dEdx = \pm 5 sigma, min pt charged pi = 100 MeV
    cuts.AddCut("80000113","00200009117000008260400000","1111111067032230000","002013700","0103503000000000","0103503000000000");
  } else if( trainConfig == 4 ) {
    // closing charged pion cuts, minimum TPC cluster = 80, TPC dEdx pi = \pm 3 sigma, ITS dEdx = \pm 4 sigma, min pt charged pi = 100 MeV
    cuts.AddCut("80000113","00200009117000008260400000","1111111067032230000","002016700","0103503000000000","0103503000000000");
  } else if( trainConfig == 5 ) {
    // closing charged pion cuts, minimum TPC cluster = 80, TPC dEdx pi = \pm 3 sigma, ITS dEdx = \pm 4 sigma, min pt charged pi = 100 MeV
    // closing neural pion cuts, 0.1 < M_gamma,gamma < 0.145
    cuts.AddCut("80000113","00200009117000008260400000","1111111067032230000","002016700","0103503100000000","0103503000000000");
  } else if( trainConfig == 6 ) {
    // closing charged pion cuts, minimum TPC cluster = 80, TPC dEdx pi = \pm 3 sigma, ITS dEdx = \pm 4 sigma, min pt charged pi = 100 MeV
    // closing neural pion cuts, 0.11 < M_gamma,gamma < 0.145
    cuts.AddCut("80000113","00200009117000008260400000","1111111067032230000","002016700","0103503200000000","0103503000000000");
  } else if( trainConfig == 7 ) {
    // closing charged pion cuts, minimum TPC cluster = 80, TPC dEdx pi = \pm 3 sigma, ITS dEdx = \pm 4 sigma, min pt charged pi = 100 MeV
    // closing neural pion cuts, 0.12 < M_gamma,gamma < 0.145
    cuts.AddCut("80000113","00200009117000008260400000","1111111067032230000","002016700","0103503300000000","0103503000000000");
  } else if( trainConfig == 8 ) {
    // closing charged pion cuts, minimum TPC cluster = 80, TPC dEdx pi = \pm 3 sigma, ITS dEdx = \pm 5 sigma, min pt charged pi = 100 MeV
    // closing neural pion cuts, 0.12 < M_gamma,gamma < 0.145
    cuts.AddCut("80000113","00200009117000008260400000","1111111067032230000","002013700","0103503300000000","0103503000000000");
  } else if( trainConfig == 9 ) {
    // closing charged pion cuts, minimum TPC cluster = 80, TPC dEdx pi = \pm 3 sigma, pi+pi- mass Cut at 0.75, min pt charged pi = 100 MeV
    // closing neural pion cuts, 0.1 < M_gamma,gamma < 0.145
    cuts.AddCut("80000113","00200009117000008260400000","1111111067032230000","002010702","0103503100000000","0103503000000000");
  } else if( trainConfig == 10 ){
    // eta < 0.9
    // closing charged pion cuts, minimum TPC cluster = 80, TPC dEdx pi = \pm 3 sigma, pi+pi- mass cut of 0.65, min pt charged pi = 100 MeV
    // closing neural pion cuts, 0.1 < M_gamma,gamma < 0.15
    // maxChi2 per cluster TPC <4, require TPC refit, DCA XY pT dependend 0.0182+0.0350/pt^1.01, DCA_Z = 3.0
    cuts.AddCut("80000113","00200009327000008250400000","1111111067032230000","302010706","0103503400000000","0153503000000000"); // only cuts above
    cuts.AddCut("80000113","00200009327000008250400000","1111111067032230000","30a010706","0103503400000000","0153503000000000"); // above + fMinClsTPC=80.+ fChi2PerClsTPC=4 + fRequireTPCRefit=kTRUE;
    cuts.AddCut("80000113","00200009327000008250400000","1111111067032230000","302310706","0103503400000000","0153503000000000"); // above + DCA pT dependent 0.0182+0.0350/pt^1.01 + DCA_Z < 3.0
    cuts.AddCut("80000113","00200009327000008250400000","1111111067032230000","302030706","0103503400000000","0153503000000000"); // above + pTmin=0.15
    cuts.AddCut("80000113","00200009327000008250400000","1111111067032230000","30a330706","0103503400000000","0153503000000000"); // all of the above
  } else {
    Error(Form("GammaConvNeutralMeson_MixedMode_%i",trainConfig), "wrong trainConfig variable no cuts have been specified for the configuration");
    return;
  }

  // PHOS mode
  if( trainConfig == 31 ) {
    // everything open
    cuts.AddCut("80000113","00200009117000008260400000","2444400030022000000","000010400","0103503000000000","0103503000000000");
  } else if( trainConfig == 32 ) {
    // closing charged pion cuts, minimum TPC cluster = 80, TPC dEdx pi = \pm 3 sigma, min pt charged pi = 100 MeV
    cuts.AddCut("80000113","00200009117000008260400000","2444400030022000000","002010700","0103503000000000","0103503000000000");
  } else if( trainConfig == 33 ) {
    // closing charged pion cuts, minimum TPC cluster = 80, TPC dEdx pi = \pm 3 sigma, ITS dEdx = \pm 5 sigma, min pt charged pi = 100 MeV
    cuts.AddCut("80000113","00200009117000008260400000","2444400030022000000","002013700","0103503000000000","0103503000000000");
  } else if( trainConfig == 34 ) {
    // closing charged pion cuts, minimum TPC cluster = 80, TPC dEdx pi = \pm 3 sigma, ITS dEdx = \pm 4 sigma, min pt charged pi = 100 MeV
    cuts.AddCut("80000113","00200009117000008260400000","2444400030022000000","002016700","0103503000000000","0103503000000000");
  } else if( trainConfig == 35 ) {
    // closing charged pion cuts, minimum TPC cluster = 80, TPC dEdx pi = \pm 3 sigma, ITS dEdx = \pm 4 sigma, min pt charged pi = 100 MeV
    // closing neural pion cuts, 0.1 < M_gamma,gamma < 0.145
    cuts.AddCut("80000113","00200009117000008260400000","2444400030022000000","002016700","0103503100000000","0103503000000000");
  } else if( trainConfig == 36 ) {
    // closing charged pion cuts, minimum TPC cluster = 80, TPC dEdx pi = \pm 3 sigma, ITS dEdx = \pm 4 sigma, min pt charged pi = 100 MeV
    // closing neural pion cuts, 0.11 < M_gamma,gamma < 0.145
    cuts.AddCut("80000113","00200009117000008260400000","2444400030022000000","002016700","0103503200000000","0103503000000000");
  } else if( trainConfig == 37 ) {
    // closing charged pion cuts, minimum TPC cluster = 80, TPC dEdx pi = \pm 3 sigma, ITS dEdx = \pm 4 sigma, min pt charged pi = 100 MeV
    // closing neural pion cuts, 0.12 < M_gamma,gamma < 0.145
    cuts.AddCut("80000113","00200009117000008260400000","2444400030022000000","002016700","0103503300000000","0103503000000000");
  } else if( trainConfig == 38 ) {
    // closing charged pion cuts, minimum TPC cluster = 80, TPC dEdx pi = \pm 3 sigma, ITS dEdx = \pm 5 sigma, min pt charged pi = 100 MeV
    // closing neural pion cuts, 0.12 < M_gamma,gamma < 0.145
    cuts.AddCut("80000113","00200009117000008260400000","2444400030022000000","002013700","0103503300000000","0103503000000000");
  } else if( trainConfig == 39 ) {
    // closing charged pion cuts, minimum TPC cluster = 80, TPC dEdx pi = \pm 3 sigma, pi+pi- mass Cut at 0.75, min pt charged pi = 100 MeV
    // closing neural pion cuts, 0.1 < M_gamma,gamma < 0.145
    cuts.AddCut("80000113","00200009117000008260400000","2444400030022000000","002010702","0103503100000000","0103503000000000");
  } else if( trainConfig == 49) {
    // eta < 0.9
    // closing charged pion cuts, minimum TPC cluster = 80, TPC dEdx pi = \pm 3 sigma, pi+pi- mass cut of 0.65, min pt charged pi = 100 MeV
    // closing neural pion cuts, 0.1 < M_gamma,gamma < 0.15
    // maxChi2 per cluster TPC <4, require TPC refit, DCA XY pT dependend 0.0182+0.0350/pt^1.01, DCA_Z = 3.0
    cuts.AddCut("80000113","00200009327000008250400000","2444400043013300000","302010706","0103503400000000","0153503000000000"); // only cuts above
    cuts.AddCut("80000113","00200009327000008250400000","2444400043013300000","30a010706","0103503400000000","0153503000000000"); // above + fMinClsTPC=80.+ fChi2PerClsTPC=4 + fRequireTPCRefit=kTRUE;
    cuts.AddCut("80000113","00200009327000008250400000","2444400043013300000","302310706","0103503400000000","0153503000000000"); // above + DCA pT dependent 0.0182+0.0350/pt^1.01 + DCA_Z < 3.0
    cuts.AddCut("80000113","00200009327000008250400000","2444400043013300000","302030706","0103503400000000","0153503000000000"); // above + pTmin=0.15
    cuts.AddCut("80000113","00200009327000008250400000","2444400043013300000","30a330706","0103503400000000","0153503000000000"); // all of the above
  } else {
    Error(Form("GammaConvNeutralMeson_MixedMode_%i",trainConfig), "wrong trainConfig variable no cuts have been specified for the configuration");
    return;
  }

  if(!cuts.AreValid()){
    cout << "\n\n****************************************************" << endl;
    cout << "ERROR: No valid cuts stored in CutHandlerNeutralMixed! Returning..." << endl;
    cout << "****************************************************\n\n" << endl;
    return;
  }

  Int_t numberOfCuts = cuts.GetNCuts();

  TList *EventCutList = new TList();
  TList *ConvCutList  = new TList();
  TList *ClusterCutList  = new TList();
  TList *NeutralPionCutList = new TList();
  TList *MesonCutList = new TList();
  TList *PionCutList  = new TList();

  TList *HeaderList = new TList();
  TObjString *Header1 = new TObjString("pi0_1");
  HeaderList->Add(Header1);
  TObjString *Header3 = new TObjString("eta_2");
  HeaderList->Add(Header3);
  
  EventCutList->SetOwner(kTRUE);
  AliConvEventCuts **analysisEventCuts = new AliConvEventCuts*[numberOfCuts];
  ConvCutList->SetOwner(kTRUE);
  AliConversionPhotonCuts **analysisCuts = new AliConversionPhotonCuts*[numberOfCuts];
  ClusterCutList->SetOwner(kTRUE);
  AliCaloPhotonCuts **analysisClusterCuts = new AliCaloPhotonCuts*[numberOfCuts];
  NeutralPionCutList->SetOwner(kTRUE);
  AliConversionMesonCuts **analysisNeutralPionCuts   = new AliConversionMesonCuts*[numberOfCuts];
  MesonCutList->SetOwner(kTRUE);
  AliConversionMesonCuts **analysisMesonCuts   = new AliConversionMesonCuts*[numberOfCuts];
  PionCutList->SetOwner(kTRUE);
  AliPrimaryPionCuts **analysisPionCuts     = new AliPrimaryPionCuts*[numberOfCuts];

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
    analysisEventCuts[i]->SetV0ReaderName(V0ReaderName);
    analysisEventCuts[i]->InitializeCutsFromCutString((cuts.GetEventCut(i)).Data());
    EventCutList->Add(analysisEventCuts[i]);
    analysisEventCuts[i]->SetFillCutHistograms("",kFALSE);

    analysisCuts[i] = new AliConversionPhotonCuts();
    analysisCuts[i]->SetV0ReaderName(V0ReaderName);
    if( ! analysisCuts[i]->InitializeCutsFromCutString((cuts.GetConversionCut(i)).Data()) ) {
      cout<<"ERROR: analysisCuts [" <<i<<"]"<<endl;
      return 0;
    } else {				
      ConvCutList->Add(analysisCuts[i]);
      analysisCuts[i]->SetFillCutHistograms("",kFALSE);	
    }

    analysisClusterCuts[i] = new AliCaloPhotonCuts();
    analysisClusterCuts[i]->SetV0ReaderName(V0ReaderName);
    analysisClusterCuts[i]->SetCaloTrackMatcherName(TrackMatcherName);
    if( ! analysisClusterCuts[i]->InitializeCutsFromCutString((cuts.GetClusterCut(i)).Data()) ) {
      cout<<"ERROR: analysisClusterCuts [" <<i<<"]"<<endl;
      return 0;
    } else {				
      ClusterCutList->Add(analysisClusterCuts[i]);
      analysisClusterCuts[i]->SetFillCutHistograms("");			
    }

    analysisNeutralPionCuts[i] = new AliConversionMesonCuts();
    if( ! analysisNeutralPionCuts[i]->InitializeCutsFromCutString((cuts.GetNeutralPionCut(i)).Data()) ) {
      cout<<"ERROR: analysisMesonCuts [ " <<i<<" ] "<<endl;
      return 0;
    } else {
      NeutralPionCutList->Add(analysisNeutralPionCuts[i]);
      analysisNeutralPionCuts[i]->SetFillCutHistograms("");
    }
  
    analysisMesonCuts[i] = new AliConversionMesonCuts();
    if( ! analysisMesonCuts[i]->InitializeCutsFromCutString((cuts.GetMesonCut(i)).Data()) ) {
      cout<<"ERROR: analysisMesonCuts [ " <<i<<" ] "<<endl;
      return 0;
    } else {
      MesonCutList->Add(analysisMesonCuts[i]);
      analysisMesonCuts[i]->SetFillCutHistograms("");
    }
    analysisEventCuts[i]->SetAcceptedHeader(HeaderList);
    
    TString cutName( Form("%s_%s_%s_%s_%s_%s",(cuts.GetEventCut(i)).Data(), (cuts.GetConversionCut(i)).Data(), (cuts.GetClusterCut(i)).Data(),(cuts.GetPionCut(i)).Data(),(cuts.GetNeutralPionCut(i)).Data(), (cuts.GetMesonCut(i)).Data() ) );
    analysisPionCuts[i] = new AliPrimaryPionCuts();
    if( !analysisPionCuts[i]->InitializeCutsFromCutString((cuts.GetPionCut(i)).Data())) {
      cout<< "ERROR:  analysisPionCuts [ " <<i<<" ] "<<endl;
      return 0;
    } else { 
      PionCutList->Add(analysisPionCuts[i]);
      analysisPionCuts[i]->SetFillCutHistograms("",kFALSE,cutName); 
    }
  }

  task->SetNeutralPionMode(neutralPionMode);
  task->SetEventCutList(numberOfCuts,EventCutList);
  task->SetConversionCutList(ConvCutList);
  task->SetClusterCutList(ClusterCutList);
  task->SetNeutralPionCutList(NeutralPionCutList);
  task->SetMesonCutList(MesonCutList);
  task->SetPionCutList(PionCutList);

  task->SetMoveParticleAccordingToVertex(kTRUE);

  if(enableQAMesonTask) task->SetDoMesonQA(kTRUE);

  //connect containers
  AliAnalysisDataContainer *coutput =
  mgr->CreateContainer(Form("GammaConvNeutralMesonPiPlPiMiPiZero_%i_%i",neutralPionMode, trainConfig), TList::Class(),
              AliAnalysisManager::kOutputContainer,Form("GammaConvNeutralMesonPiPlPiMiPiZero_%i_%i.root",neutralPionMode, trainConfig));

  mgr->AddTask(task);
  mgr->ConnectInput(task,0,cinput);
  mgr->ConnectOutput(task,1,coutput);

  return;

}
