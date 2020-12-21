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

    void AddCut(TString eventCut, TString conversionCut, TString clusterCut, TString pionCut, TString neutralPionCut, TString mesonCut){
      if(nCuts>=nMaxCuts) {cout << "ERROR in CutHandlerNeutralMixed: Exceeded maximum number of cuts!" << endl; validCuts = false; return;}
      if( eventCut.Length()!=8 || conversionCut.Length()!=26 || clusterCut.Length()!=19 || pionCut.Length()!=9 || neutralPionCut.Length()!=16 || mesonCut.Length()!=16 ) {cout << "ERROR in CutHandlerNeutralMixed: Incorrect length of cut string!" << endl; validCuts = false; return;}
      eventCutArray[nCuts]=eventCut; conversionCutArray[nCuts]=conversionCut; clusterCutArray[nCuts]=clusterCut; pionCutArray[nCuts]=pionCut; neutralPionCutArray[nCuts]=neutralPionCut; mesonCutArray[nCuts]=mesonCut;
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
void AddTask_GammaConvNeutralMesonPiPlPiMiPiZero_MixedMode_pp(
    Int_t trainConfig                 = 1,
    Bool_t isMC                       = kFALSE,                         //run MC
    Int_t enableQAMesonTask          = 1,                               //enable QA in AliAnalysisTaskNeutralMesonToPiPlPiMiPiZero
    TString fileNameInputForWeighting = "MCSpectraInput.root",          // path to file for weigting input
    Bool_t doWeighting                = kFALSE,                         //enable Weighting
    TString generatorName             = "HIJING",
    TString cutnumberAODBranch        = "000000006008400001001500000",
    Double_t tolerance                = -1,
    TString periodNameV0Reader        = "",                              // period Name for V0Reader
    Int_t   runLightOutput            = 0,                               // run light output option 0: no light output 1: most cut histos stiched off 2: unecessary omega hists turned off as well
    TString additionalTrainConfig     = "0"                              // additional counter for trainconfig, this has to be always the last parameter
  ) {

  Int_t trackMatcherRunningMode = 0; // CaloTrackMatcher running mode
  //parse additionalTrainConfig flag
  TObjArray *rAddConfigArr = additionalTrainConfig.Tokenize("_");
  if(rAddConfigArr->GetEntries()<1){cout << "ERROR during parsing of additionalTrainConfig String '" << additionalTrainConfig.Data() << "'" << endl; return;}
  TObjString* rAdditionalTrainConfig;
  for(Int_t i = 0; i<rAddConfigArr->GetEntries() ; i++){
    if(i==0) rAdditionalTrainConfig = (TObjString*)rAddConfigArr->At(i);
    else{
      TObjString* temp = (TObjString*) rAddConfigArr->At(i);
      TString tempStr = temp->GetString();
      if(tempStr.BeginsWith("TM")){
        TString tempType = tempStr;
        tempType.Replace(0,2,"");
        trackMatcherRunningMode = tempType.Atoi();
        cout << Form("INFO: AddTask_GammaConvNeutralMesonPiPlPiMiPiZero_MixedMode_pp will use running mode '%i' for the TrackMatcher!",trackMatcherRunningMode) << endl;
      }
    }
  }
  TString sAdditionalTrainConfig = rAdditionalTrainConfig->GetString();
  if (sAdditionalTrainConfig.Atoi() > 0){
    trainConfig = trainConfig + sAdditionalTrainConfig.Atoi();
    cout << "INFO: AddTask_GammaConvNeutralMesonPiPlPiMiPiZero_MixedMode_pp running additionalTrainConfig '" << sAdditionalTrainConfig.Atoi() << "', train config: '" << trainConfig << "'" << endl;
  }

  Int_t isHeavyIon = 0;
  Int_t neutralPionMode = 1;

  // ================== GetAnalysisManager ===============================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error(Form("AddTask_GammaConvNeutralMesonPiPlPiMiPiZero_MixedMode_pp_%i",trainConfig), "No analysis manager found.");
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
  if (  periodNameV0Reader.CompareTo("LHC16f") == 0 || periodNameV0Reader.CompareTo("LHC17g")==0 || periodNameV0Reader.CompareTo("LHC18c")==0 ||
        periodNameV0Reader.CompareTo("LHC17d1") == 0  || periodNameV0Reader.CompareTo("LHC17d12")==0 ||
        periodNameV0Reader.CompareTo("LHC17h3")==0 || periodNameV0Reader.CompareTo("LHC17k1")==0 ||
        periodNameV0Reader.CompareTo("LHC17f8b") == 0 ||
        periodNameV0Reader.CompareTo("LHC16P1JJLowB") == 0 || periodNameV0Reader.CompareTo("LHC16P1Pyt8LowB") == 0 )
    cutnumberPhoton         = "00000088400000000100000000";

  TString cutnumberEvent = "00000003";
  TString PionCuts      = "000000200";            //Electron Cuts


  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();

  //========= Add V0 Reader to  ANALYSIS manager if not yet existent =====
  TString V0ReaderName = Form("V0ReaderV1_%s_%s",cutnumberEvent.Data(),cutnumberPhoton.Data());
  if( !(AliV0ReaderV1*)mgr->GetTask(V0ReaderName.Data()) ){
    AliV0ReaderV1 *fV0ReaderV1 = new AliV0ReaderV1(V0ReaderName.Data());

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
      fEventCuts->SetV0ReaderName(V0ReaderName);
      if (periodNameV0Reader.CompareTo("") != 0) fEventCuts->SetPeriodEnum(periodNameV0Reader);
      if(runLightOutput>0) fEventCuts->SetLightOutput(kTRUE);
      if(fEventCuts->InitializeCutsFromCutString(cutnumberEvent.Data())){
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
      if(runLightOutput>0) fCuts->SetLightOutput(kTRUE);
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
  //========= Add Pion Selector ====================
  if( !(AliPrimaryPionSelector*)mgr->GetTask("PionSelector") ){
    AliPrimaryPionSelector *fPionSelector = new AliPrimaryPionSelector("PionSelector");
    AliPrimaryPionCuts *fPionCuts=0;
    if( PionCuts!=""){
      fPionCuts= new AliPrimaryPionCuts(PionCuts.Data(),PionCuts.Data());
      if(runLightOutput>0) fPionCuts->SetLightOutput(kTRUE);

      if(fPionCuts->InitializeCutsFromCutString(PionCuts.Data())){
        fPionSelector->SetPrimaryPionCuts(fPionCuts);
        fPionCuts->SetFillCutHistograms("",kTRUE);
      }
    }

    fPionSelector->Init();
    mgr->AddTask(fPionSelector);

    AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer();
    mgr->ConnectInput (fPionSelector,0,cinput1);
  }

  AliAnalysisTaskNeutralMesonToPiPlPiMiPiZero *task=NULL;
  task= new AliAnalysisTaskNeutralMesonToPiPlPiMiPiZero(Form("GammaConvNeutralMesonPiPlPiMiPiZero_%i_%i",neutralPionMode, trainConfig));
  task->SetIsHeavyIon(isHeavyIon);
  task->SetIsMC(isMC);
  task->SetV0ReaderName(V0ReaderName);
  if(runLightOutput>1) task->SetLightOutput(kTRUE);
  task->SetTolerance(tolerance);
  task->SetTrackMatcherRunningMode(trackMatcherRunningMode);

  CutHandlerNeutralMixed cuts;

  // 7 TeV
  // EMCAL mode
  if( trainConfig == 1 ) {
    // everything open
    cuts.AddCut("00000113","00200009327000008250400000","1111113047032230000","000010400","0103503a00000000","0103503000000000");
  } else if( trainConfig == 2 ) {
    // closing charged pion cuts, minimum TPC cluster = 80, TPC dEdx pi = \pm 3 sigma, min pt charged pi = 100 MeV
    cuts.AddCut("00000113","00200009327000008250400000","1111113047032230000","002010700","0103503a00000000","0103503000000000");
  } else if( trainConfig == 3 ) {
    // closing charged pion cuts, minimum TPC cluster = 80, TPC dEdx pi = \pm 3 sigma, ITS dEdx = \pm 5 sigma, min pt charged pi = 100 MeV
    cuts.AddCut("00000113","00200009327000008250400000","1111113047032230000","002013700","0103503a00000000","0103503000000000");
  } else if( trainConfig == 4 ) {
    // closing charged pion cuts, minimum TPC cluster = 80, TPC dEdx pi = \pm 3 sigma, ITS dEdx = \pm 4 sigma, min pt charged pi = 100 MeV
    cuts.AddCut("00000113","00200009327000008250400000","1111113047032230000","002016700","0103503a00000000","0103503000000000");
  } else if( trainConfig == 5 ) {
    // closing charged pion cuts, minimum TPC cluster = 80, TPC dEdx pi = \pm 3 sigma, ITS dEdx = \pm 4 sigma, min pt charged pi = 100 MeV
    // closing neural pion cuts, 0.1 < M_gamma,gamma < 0.15
    cuts.AddCut("00000113","00200009327000008250400000","1111113047032230000","002016700","0103503400000000","0103503000000000");
  } else if( trainConfig == 6 ) {
    // closing charged pion cuts, minimum TPC cluster = 80, TPC dEdx pi = \pm 3 sigma, ITS dEdx = \pm 4 sigma, min pt charged pi = 100 MeV
    // closing neural pion cuts, 0.11 < M_gamma,gamma < 0.145
    cuts.AddCut("00000113","00200009327000008250400000","1111113047032230000","002016700","0103503200000000","0103503000000000");
  } else if( trainConfig == 7 ) {
    // closing charged pion cuts, minimum TPC cluster = 80, TPC dEdx pi = \pm 3 sigma, ITS dEdx = \pm 4 sigma, min pt charged pi = 100 MeV
    // closing neural pion cuts, 0.12 < M_gamma,gamma < 0.145
    cuts.AddCut("00000113","00200009327000008250400000","1111113047032230000","002016700","0103503300000000","0103503000000000");
  } else if( trainConfig == 8 ) {
    // closing charged pion cuts, minimum TPC cluster = 80, TPC dEdx pi = \pm 3 sigma, ITS dEdx = \pm 5 sigma, min pt charged pi = 100 MeV
    // closing neural pion cuts, 0.12 < M_gamma,gamma < 0.145
    cuts.AddCut("00000113","00200009327000008250400000","1111113047032230000","002013700","0103503300000000","0103503000000000");
  } else if( trainConfig == 9 ) {
    // closing charged pion cuts, minimum TPC cluster = 80, TPC dEdx pi = \pm 3 sigma, pi+pi- mass Cut at 0.65, min pt charged pi = 100 MeV
    // closing neural pion cuts, 0.1 < M_gamma,gamma < 0.15
    cuts.AddCut("00000113","00200009327000008250400000","1111113047032230000","002010706","0103503400000000","0153503000000000");
  } else if( trainConfig == 10 ) {
    // closing charged pion cuts, minimum TPC cluster = 80, TPC dEdx pi = \pm 3 sigma, pi+pi- mass Cut at 0.6, min pt charged pi = 100 MeV
    // closing neural pion cuts, 0.1 < M_gamma,gamma < 0.15
    cuts.AddCut("00000113","00200009327000008250400000","1111113047032230000","002010703","0103503400000000","0103503000000000");
  } else if( trainConfig == 11 ) {
    // closing charged pion cuts, minimum TPC cluster = 80, TPC dEdx pi = \pm 3 sigma, pi+pi- mass Cut at 0.5, min pt charged pi = 100 MeV
    // closing neural pion cuts, 0.1 < M_gamma,gamma < 0.15
    cuts.AddCut("00000113","00200009327000008250400000","1111113047032230000","002010705","0103503400000000","0103503000000000");

  } else if( trainConfig == 20 ) {
    // closing charged pion cuts, minimum TPC cluster = 80, TPC dEdx pi = \pm 3 sigma, pi+pi- mass Cut at 0.65, min pt charged pi = 100 MeV
    // closing neural pion cuts, 0.1 < M_gamma,gamma < 0.15
    cuts.AddCut("00000113","00200009327000008250400000","1111113047032230000","002010706","0103503400000000","0103503000000000");
  } else if( trainConfig == 21 ) {
    // closing charged pion cuts, minimum TPC cluster = 80, TPC dEdx pi = \pm 3 sigma, pi+pi- mass Cut at 0.65, min pt charged pi = 100 MeV
    // closing neural pion cuts, 0.1 < M_gamma,gamma < 0.15
    cuts.AddCut("00000113","00200009327000008250400000","1111113047032230000","002010706","0103503400000000","0103503000000000");
  } else if( trainConfig == 22 ) {
    // closing charged pion cuts, minimum TPC cluster = 80, TPC dEdx pi = \pm 3 sigma, pi+pi- mass Cut at 0.65, min pt charged pi = 100 MeV
    // closing neural pion cuts, 0.1 < M_gamma,gamma < 0.15
    cuts.AddCut("00000113","00200009327000008250400000","1111113047032230000","002010706","0103503400000000","0103503000000000");
  } else if( trainConfig == 23 ) {
    // closing charged pion cuts, minimum TPC cluster = 80, TPC dEdx pi = \pm 3 sigma, pi+pi- mass Cut at 0.65, min pt charged pi = 100 MeV
    // closing neural pion cuts, 0.1 < M_gamma,gamma < 0.15
    cuts.AddCut("00000113","00200009327000008250400000","1111113047032230000","002010706","0103503400000000","0103503000000000");
  } else if( trainConfig == 24 ) {
    // closing charged pion cuts, minimum TPC cluster = 80, TPC dEdx pi = \pm 3 sigma, pi+pi- mass Cut at 0.65, min pt charged pi = 100 MeV
    // closing neural pion cuts, 0.1 < M_gamma,gamma < 0.15
    cuts.AddCut("00000113","00200009327000008250400000","1111113047032230000","002010706","0103503400000000","0103503000000000");
  } else if( trainConfig == 25 ) {
    // closing charged pion cuts, minimum TPC cluster = 80, TPC dEdx pi = \pm 3 sigma, pi+pi- mass Cut at 0.65, min pt charged pi = 100 MeV
    // closing neural pion cuts, 0.1 < M_gamma,gamma < 0.15
    cuts.AddCut("00000113","00200009327000008250400000","1111113047032230000","002010706","0103503400000000","0103503000000000");
  } else if( trainConfig == 26 ) {
    // closing charged pion cuts, minimum TPC cluster = 80, TPC dEdx pi = \pm 3 sigma, pi+pi- mass Cut at 0.65, min pt charged pi = 100 MeV
    // closing neural pion cuts, 0.1 < M_gamma,gamma < 0.15
    cuts.AddCut("00000113","00200009327000008250400000","1111113047032230000","002010706","0103503400000000","0103503000000000");
  } else if( trainConfig == 27 ) {
    // closing charged pion cuts, minimum TPC cluster = 80, TPC dEdx pi = \pm 3 sigma, pi+pi- mass Cut at 0.65, min pt charged pi = 100 MeV
    // closing neural pion cuts, 0.1 < M_gamma,gamma < 0.15
    cuts.AddCut("00000113","00200009327000008250400000","1111113047032230000","002010706","0103503400000000","0103503000000000");
  } else if( trainConfig == 28 ) {
    // eta < 0.9
    // closing charged pion cuts, minimum TPC cluster = 80, TPC dEdx pi = \pm 3 sigma, pi+pi- mass cut of 0.75, min pt charged pi = 100 MeV
    // closing neural pion cuts, 0.120 < M_gamma,gamma < 0.15
    // maxChi2 per cluster TPC <4, require TPC refit, DCA XY pT dependend 0.0182+0.0350/pt^1.01, DCA_Z = 3.0
    cuts.AddCut("00000113","00200009327000008250400000","1111113047032230000","302010702","0103503600000000","0153503000000000"); // only cuts above
    cuts.AddCut("00000113","00200009327000008250400000","1111113047032230000","302010702","0103503500000000","0153503000000000"); // 0.120 < M_gamma,gamma < 0.15
    cuts.AddCut("00000113","00200009327000008250400000","1111113047032230000","302010702","0103503400000000","0153503000000000"); // 0.100 < M_gamma,gamma < 0.15
    //cuts.AddCut("00000113","00200009327000008250400000","1111113047032230000","30a010706","0103503600000000","0153503000000000"); // above + fMinClsTPC=80.+ fChi2PerClsTPC=4 + fRequireTPCRefit=kTRUE;
    //cuts.AddCut("00000113","00200009327000008250400000","1111113047032230000","302310706","0103503600000000","0153503000000000"); // above + DCA pT dependent 0.0182+0.0350/pt^1.01 + DCA_Z < 3.0
    //cuts.AddCut("00000113","00200009327000008250400000","1111113047032230000","302030706","0103503600000000","0153503000000000"); // above + pTmin=0.15
    //cuts.AddCut("00000113","00200009327000008250400000","1111113047032230000","30a330706","0103503600000000","0153503000000000"); // all of the above
  } else if( trainConfig == 29 ) {
    // eta < 0.9
    // closing charged pion cuts, minimum TPC cluster = 80, TPC dEdx pi = \pm 3 sigma, pi+pi- mass cut of 0.85, min pt charged pi = 100 MeV
    // closing neural pion cuts, 0.110 < M_gamma,gamma < 0.155
    // maxChi2 per cluster TPC <4, require TPC refit, DCA XY pT dependend 0.0182+0.0350/pt^1.01, DCA_Z = 3.0
    cuts.AddCut("00000113","00200009327000008250400000","1111113047032230000","302010708","0103503900000000","0153503000000000"); // normal mixing
    cuts.AddCut("00000113","00200009327000008250400000","1111113047032230000","302010708","0103503900000000","0d53503000000000"); // pi0 sideband mixing both sides
//  cuts.AddCut("00000113","00200009327000008250400000","1111113047032230000","302040708","0103503900000000","0153503000000000"); // normal mixing + charged pi pt cut > 400 MeV
//  cuts.AddCut("00000113","00200009327000008250400000","1111113047032230000","302010708","0103503900000000","0a53503000000000"); // likesign mixing
//  cuts.AddCut("00000113","00200009327000008250400000","1111113047032230000","302010708","0103503900000000","0b53503000000000"); // pi0 sideband mixing right (0.180-0.220)
//  cuts.AddCut("00000113","00200009327000008250400000","1111113047032230000","302010708","0103503900000000","0c53503000000000"); // pi0 sideband mixing left  (0.01-0.05)
  } else if( trainConfig == 30 ) {
    // eta < 0.9
    // closing charged pion cuts, minimum TPC cluster = 80, TPC dEdx pi = \pm 3 sigma, pi+pi- mass cut of 0.85, min pt charged pi = 100 MeV
    // closing neural pion cuts, 0.110 < M_gamma,gamma < 0.155
    // maxChi2 per cluster TPC <4, require TPC refit, DCA XY pT dependend 0.0182+0.0350/pt^1.01, DCA_Z = 3.0
    // nmb bck events 80
    cuts.AddCut("00000113","00200009327000008250400000","1111113047032230000","302010708","0103603900000000","0153503000000000"); // normal mixing
    cuts.AddCut("00000113","00200009327000008250400000","1111113047032230000","302010708","0103613900000000","0153503000000000"); // pT min pi0 = 0.4
    cuts.AddCut("00000113","00200009327000008250400000","1111113047032230000","302010708","0103623900000000","0153503000000000"); // pT min pi0 = 0.7
    cuts.AddCut("00000113","00200009327000008250400000","1111113047032230000","302010708","0103653900000000","0153503000000000"); // pT min pi0 = 1.2
    cuts.AddCut("00000113","00200009327000008250400000","1111113047032230000","302010708","0103663900000000","0153503000000000"); // pT min pi0 = 1.5

  // PHOS mode
  } else if( trainConfig == 31 ) {
    // everything open
    cuts.AddCut("00000113","00200009327000008250400000","2444400043013300000","000010400","0103503a00000000","0103503000000000");
  } else if( trainConfig == 32 ) {
    // closing charged pion cuts, minimum TPC cluster = 80, TPC dEdx pi = \pm 3 sigma, min pt charged pi = 100 MeV
    cuts.AddCut("00000113","00200009327000008250400000","2444400043013300000","002010700","0103503a00000000","0103503000000000");
  } else if( trainConfig == 33 ) {
    // closing charged pion cuts, minimum TPC cluster = 80, TPC dEdx pi = \pm 3 sigma, ITS dEdx = \pm 5 sigma, min pt charged pi = 100 MeV
    cuts.AddCut("00000113","00200009327000008250400000","2444400043013300000","002013700","0103503a00000000","0103503000000000");
  } else if( trainConfig == 34 ) {
    // closing charged pion cuts, minimum TPC cluster = 80, TPC dEdx pi = \pm 3 sigma, ITS dEdx = \pm 4 sigma, min pt charged pi = 100 MeV
    cuts.AddCut("00000113","00200009327000008250400000","2444400043013300000","002016700","0103503a00000000","0103503000000000");
  } else if( trainConfig == 35 ) {
    // closing charged pion cuts, minimum TPC cluster = 80, TPC dEdx pi = \pm 3 sigma, ITS dEdx = \pm 4 sigma, min pt charged pi = 100 MeV
    // closing neural pion cuts, 0.1 < M_gamma,gamma < 0.15
    cuts.AddCut("00000113","00200009327000008250400000","2444400043013300000","002016700","0103503400000000","0103503000000000");
  } else if( trainConfig == 36 ) {
    // closing charged pion cuts, minimum TPC cluster = 80, TPC dEdx pi = \pm 3 sigma, ITS dEdx = \pm 4 sigma, min pt charged pi = 100 MeV
    // closing neural pion cuts, 0.11 < M_gamma,gamma < 0.145
    cuts.AddCut("00000113","00200009327000008250400000","2444400043013300000","002016700","0103503200000000","0103503000000000");
  } else if( trainConfig == 37 ) {
    // closing charged pion cuts, minimum TPC cluster = 80, TPC dEdx pi = \pm 3 sigma, ITS dEdx = \pm 4 sigma, min pt charged pi = 100 MeV
    // closing neural pion cuts, 0.12 < M_gamma,gamma < 0.145
    cuts.AddCut("00000113","00200009327000008250400000","2444400043013300000","002016700","0103503300000000","0103503000000000");
  } else if( trainConfig == 38 ) {
    // closing charged pion cuts, minimum TPC cluster = 80, TPC dEdx pi = \pm 3 sigma, ITS dEdx = \pm 5 sigma, min pt charged pi = 100 MeV
    // closing neural pion cuts, 0.12 < M_gamma,gamma < 0.145
    cuts.AddCut("00000113","00200009327000008250400000","2444400043013300000","002013700","0103503300000000","0103503000000000");
  } else if( trainConfig == 39 ) {
    // closing charged pion cuts, minimum TPC cluster = 80, TPC dEdx pi = \pm 3 sigma, pi+pi- mass Cut at 0.65, min pt charged pi = 100 MeV
    // closing neural pion cuts, 0.1 < M_gamma,gamma < 0.15
    cuts.AddCut("00000113","00200009327000008250400000","2444400043013300000","002010706","0103503400000000","0153503000000000");
  } else if( trainConfig == 40 ) {
    // eta < 0.9
    // closing charged pion cuts, minimum TPC cluster = 80, TPC dEdx pi = \pm 3 sigma, pi+pi- mass cut of 0.75, min pt charged pi = 100 MeV
    // closing neural pion cuts, 0.120 < M_gamma,gamma < 0.150
    // maxChi2 per cluster TPC <4, require TPC refit, DCA XY pT dependend 0.0182+0.0350/pt^1.01, DCA_Z = 3.0
    cuts.AddCut("00000113","00200009327000008250400000","2444400043013300000","302010702","0103503600000000","0153503000000000"); // only cuts above
    cuts.AddCut("00000113","00200009327000008250400000","2444400043013300000","302010702","0103503500000000","0153503000000000"); // 0.110 < M_gamma,gamma < 0.150
    cuts.AddCut("00000113","00200009327000008250400000","2444400043013300000","302010702","0103503400000000","0153503000000000"); // 0.100 < M_gamma,gamma < 0.150
    //cuts.AddCut("00000113","00200009327000008250400000","2444400043013300000","30a010706","0103503600000000","0153503000000000"); // above + fMinClsTPC=80.+ fChi2PerClsTPC=4 + fRequireTPCRefit=kTRUE;
    //cuts.AddCut("00000113","00200009327000008250400000","2444400043013300000","302310706","0103503600000000","0153503000000000"); // above + DCA pT dependent 0.0182+0.0350/pt^1.01 + DCA_Z < 3.0
    //cuts.AddCut("00000113","00200009327000008250400000","2444400043013300000","302030706","0103503600000000","0153503000000000"); // above + pTmin=0.15
    //cuts.AddCut("00000113","00200009327000008250400000","2444400043013300000","30a330706","0103503600000000","0153503000000000"); // all of the above
  } else if( trainConfig == 41 ) {
    // eta < 0.9
    // closing charged pion cuts, minimum TPC cluster = 80, TPC dEdx pi = \pm 3 sigma, pi+pi- mass cut of 0.85, min pt charged pi = 100 MeV
    // closing neural pion cuts, 0.120 < M_gamma,gamma < 0.150
    // maxChi2 per cluster TPC <4, require TPC refit, DCA XY pT dependend 0.0182+0.0350/pt^1.01, DCA_Z = 3.0
    cuts.AddCut("00000113","00200009327000008250400000","2444400043013300000","302010708","0103503600000000","0153503000000000"); // normal event mixing
    cuts.AddCut("00000113","00200009327000008250400000","2444400043013300000","302010708","0103503600000000","0d53503000000000"); // pi0 sideband mixing both sides
  //  cuts.AddCut("00000113","00200009327000008250400000","2444400043013300000","302040708","0103503600000000","0153503000000000"); // normal event mixing + charged pi pt cut > 400 MeV
  //  cuts.AddCut("00000113","00200009327000008250400000","2444400043013300000","302010708","0103503600000000","0a53503000000000"); // likesign event mixing
  //  cuts.AddCut("00000113","00200009327000008250400000","2444400043013300000","302010708","0103503600000000","0b53503000000000"); // pi0 sideband mixing right (0.180-0.220)
  //  cuts.AddCut("00000113","00200009327000008250400000","2444400043013300000","302010708","0103503600000000","0c53503000000000"); // pi0 sideband mixing left (0.01-0.05)
  } else if( trainConfig == 42 ) {
      // eta < 0.9
    // closing charged pion cuts, minimum TPC cluster = 80, TPC dEdx pi = \pm 3 sigma, pi+pi- mass cut of 0.85, min pt charged pi = 100 MeV
    // closing neural pion cuts, 0.120 < M_gamma,gamma < 0.150
    // maxChi2 per cluster TPC <4, require TPC refit, DCA XY pT dependend 0.0182+0.0350/pt^1.01, DCA_Z = 3.0
    // nmb bck events 80
    // phos nonlin gammaconv 12
    cuts.AddCut("00000113","00200009327000008250400000","2444411043013300000","302010708","0103603600000000","0153503000000000"); // normal event mixing
    cuts.AddCut("00000113","00200009327000008250400000","2444411043013300000","302010708","0103613600000000","0153503000000000"); // min pT pi0 = 0.4
    cuts.AddCut("00000113","00200009327000008250400000","2444411043013300000","302010708","0103623600000000","0153503000000000"); // min pT pi0 = 0.7
    cuts.AddCut("00000113","00200009327000008250400000","2444411043013300000","302010708","0103653600000000","0153503000000000"); // min pT pi0 = 1.2
    cuts.AddCut("00000113","00200009327000008250400000","2444411043013300000","302010708","0103663600000000","0153503000000000"); // min pT pi0 = 1.5
  } else if( trainConfig == 43 ) { // with PHOS settings from Public Note
    cuts.AddCut("00000113","00200009327000008250400000","2444400043012300000","40b440708","0103563600000000","0153503000000000"); // no NonLin
    cuts.AddCut("00000113","00200009327000008250400000","2444401043012300000","40b440708","0103563600000000","0153503000000000"); // ext PHOS NonLin
    cuts.AddCut("00000113","00200009327000008250400000","2444412043012300000","40b440708","0103563600000000","0153503000000000"); // const NonLin
  } else if( trainConfig == 44 ) { // thesis cuts PHOS
    cuts.AddCut("00000113","00200009227000008250400000","2444411043012300000","302010708","0103603600000000","0153503000000000"); //PCM-PHOS non lin
    cuts.AddCut("00000113","00200009227000008250400000","2444411043012300000","322010708","0103603600000000","0153503000000000"); //with ITS requirement
    // PCM-EMCal 7 TeV Sys
  } else if( trainConfig == 50)  { // Standard
    cuts.AddCut("00000113","00200009227000008250400000","1111111047032230000","302010708","0103603800000000","0153503000000000");
    cuts.AddCut("00000113","00200009227000008250400000","1111111047032230000","322010708","0103603800000000","0153503000000000"); // with ITS requirement

    // *************Variations in AliConvEventCuts**************************
  } else if( trainConfig == 51)  { // RemovePileUp
    cuts.AddCut("00000013","00200009227000008250400000","1111111047032230000","302010708","0103603800000000","0153503000000000");

    // *************Variations in AliConversionPhotonCuts*******************
  } else if ( trainConfig == 52) { // EtaCut
    cuts.AddCut("00000113","01200009227000008250400000","1111111047032230000","302010708","0103603800000000","0153503000000000"); // eta 0.6
    cuts.AddCut("00000113","04200009227000008250400000","1111111047032230000","302010708","0103603800000000","0153503000000000"); // eta 0.75
    cuts.AddCut("00000113","0d200009227000008250400000","1111111047032230000","302010708","0103603800000000","0153503000000000"); // eta 0.8
  } else if ( trainConfig == 53) { // SinglePtCut
    cuts.AddCut("00000113","00200019227000008250400000","1111111047032230000","302010708","0103603800000000","0153503000000000"); // 0.100 GeV
    cuts.AddCut("00000113","00200049227000008250400000","1111111047032230000","302010708","0103603800000000","0153503000000000"); // 0.075 GeV
    cuts.AddCut("00000113","00200069227000008250400000","1111111047032230000","302010708","0103603800000000","0153503000000000"); // 0.04 GeV
    cuts.AddCut("00000113","00200059227000008250400000","1111111047032230000","302010708","0103603800000000","0153503000000000"); // 0.125 GeV
  } else if ( trainConfig == 54) { // ClsTPCCut
    cuts.AddCut("00000113","00200008227000008250400000","1111111047032230000","302010708","0103603800000000","0153503000000000"); // fMinClsTPCToF= 0.35;fUseCorrectedTPCClsInfo=1;
    cuts.AddCut("00000113","00200006227000008250400000","1111111047032230000","302010708","0103603800000000","0153503000000000"); // fMinClsTPCToF= 0.70;fUseCorrectedTPCClsInfo=1;
    cuts.AddCut("00000113","00200001227000008250400000","1111111047032230000","302010708","0103603800000000","0153503000000000"); // fMinClsTPCToF= 60;fUseCorrectedTPCClsInfo=0;
    cuts.AddCut("00000113","00200002227000008250400000","1111111047032230000","302010708","0103603800000000","0153503000000000"); // fMinClsTPCToF= 80;fUseCorrectedTPCClsInfo=0;
    cuts.AddCut("00000113","00200003227000008250400000","1111111047032230000","302010708","0103603800000000","0153503000000000"); // fMinClsTPCToF= 100;fUseCorrectedTPCClsInfo=0;
  } else if ( trainConfig == 55) { // TPCdEdxCutElectron
    cuts.AddCut("00000113","00200009327000008250400000","1111111047032230000","302010708","0103603800000000","0153503000000000"); // -4,5
    cuts.AddCut("00000113","00200009627000008250400000","1111111047032230000","302010708","0103603800000000","0153503000000000"); // -2.5,4
    cuts.AddCut("00000113","00200009427000008250400000","1111111047032230000","302010708","0103603800000000","0153503000000000"); // -6,7
    cuts.AddCut("00000113","00200009527000008250400000","1111111047032230000","302010708","0103603800000000","0153503000000000"); // -4,4
    cuts.AddCut("00000113","00200009627000008250400000","1111111047032230000","302010708","0103603800000000","0153503000000000"); // -2.5,4
  } else if ( trainConfig == 56) { // TPCdEdxCutPion
    cuts.AddCut("00000113","00200009217000008250400000","1111111047032230000","302010708","0103603800000000","0153503000000000"); // fPIDnSigmaAbovePionLine=0; fPIDnSigmaAbovePionLineHighPt=-10;
    cuts.AddCut("00000113","00200009237000008250400000","1111111047032230000","302010708","0103603800000000","0153503000000000"); // fPIDnSigmaAbovePionLine=2.5; fPIDnSigmaAbovePionLineHighPt=-10;
    cuts.AddCut("00000113","00200009247000008250400000","1111111047032230000","302010708","0103603800000000","0153503000000000"); // fPIDnSigmaAbovePionLine=0.5; fPIDnSigmaAbovePionLineHighPt=-10;
    cuts.AddCut("00000113","00200009257000008250400000","1111111047032230000","302010708","0103603800000000","0153503000000000"); // fPIDnSigmaAbovePionLine=2; fPIDnSigmaAbovePionLineHighPt=-10;
  } else if ( trainConfig == 57) { // piMomdedxSigmaCut
    cuts.AddCut("00000113","00200009225000008250400000","1111111047032230000","302010708","0103603800000000","0153503000000000"); // 0.3 GeV
    cuts.AddCut("00000113","00200009220000008250400000","1111111047032230000","302010708","0103603800000000","0153503000000000"); // 0.5 GeV
    cuts.AddCut("00000113","00200009228000008250400000","1111111047032230000","302010708","0103603800000000","0153503000000000"); // 0.2 GeV
    cuts.AddCut("00000113","00200009226000008250400000","1111111047032230000","302010708","0103603800000000","0153503000000000"); // 0.25 GeV
  } else if ( trainConfig == 58) { // QtMaxCut
    cuts.AddCut("00000113","00200009227000003250400000","1111111047032230000","302010708","0103603800000000","0153503000000000"); // fQtMax=0.05; fDo2DQt=kFALSE;
    cuts.AddCut("00000113","00200009227000009250400000","1111111047032230000","302010708","0103603800000000","0153503000000000"); // fQtMax=0.03; fDo2DQt=kTRUE;
    cuts.AddCut("00000113","00200009227000002250400000","1111111047032230000","302010708","0103603800000000","0153503000000000"); // fQtMax=0.06; fDo2DQt=kTRUE;
    cuts.AddCut("00000113","00200009227000006250400000","1111111047032230000","302010708","0103603800000000","0153503000000000"); // fQtMax=0.02; fDo2DQt=kTRUE;
  } else if ( trainConfig == 59) { // Chi2GammaCut
    cuts.AddCut("00000113","00200009227000008150400000","1111111047032230000","302010708","0103603800000000","0153503000000000"); // fChi2CutConversion = 50.;
    cuts.AddCut("00000113","00200009227000008850400000","1111111047032230000","302010708","0103603800000000","0153503000000000"); // fChi2CutConversion = 20.;
    cuts.AddCut("00000113","00200009227000008a50400000","1111111047032230000","302010708","0103603800000000","0153503000000000"); // fChi2CutConversion = 25.;
    cuts.AddCut("00000113","00200009227000008950400000","1111111047032230000","302010708","0103603800000000","0153503000000000"); // fChi2CutConversion = 15.;
  } else if ( trainConfig == 60) { // PsiPair
    cuts.AddCut("00000113","00200009227000008260400000","1111111047032230000","302010708","0103603800000000","0153503000000000"); // 0.05 and 2D
    cuts.AddCut("00000113","00200009227000008280400000","1111111047032230000","302010708","0103603800000000","0153503000000000"); // 0.2 and 2D
    cuts.AddCut("00000113","00200009227000008210400000","1111111047032230000","302010708","0103603800000000","0153503000000000"); // 0.1 and 1D
    cuts.AddCut("00000113","00200009227000008220400000","1111111047032230000","302010708","0103603800000000","0153503000000000"); // 0.05 and 1D
  } else if ( trainConfig == 61) { // CosinePointingAngle
    cuts.AddCut("00000113","00200009227000008250700000","1111111047032230000","302010708","0103603800000000","0153503000000000"); // fCosPAngleCut = 0.95;
    cuts.AddCut("00000113","00200009227000008250300000","1111111047032230000","302010708","0103603800000000","0153503000000000"); // fCosPAngleCut = 0.75;
    cuts.AddCut("00000113","00200009227000008250600000","1111111047032230000","302010708","0103603800000000","0153503000000000"); // fCosPAngleCut = 0.9;

    // *************Variations in AliCaloPhotonsCut**************************
  } else if( trainConfig == 62)  { // NonLinearity
    cuts.AddCut("00000113","00200009227000008250400000","1111112047032230000","302010708","0103603800000000","0153503000000000"); // NonLinearity pp Calo - only shifting MC - no timing cut
    cuts.AddCut("00000113","00200009227000008250400000","1111113047032230000","302010708","0103603800000000","0153503000000000"); // NonLinearity ConvCalo - kTestBeamv3 + shifting MC
    cuts.AddCut("00000113","00200009227000008250400000","1111121047032230000","302010708","0103603800000000","0153503000000000"); // NonLinearity pp ConvCalo - only shifting MC - no timing cut (Fits)
    cuts.AddCut("00000113","00200009227000008250400000","1111122047032230000","302010708","0103603800000000","0153503000000000"); // NonLinearity pp ConvCalo - only shifting MC - no timing cut (Fits)
    cuts.AddCut("00000113","00200009227000008250400000","1111123047032230000","302010708","0103603800000000","0153503000000000"); // NonLinearity ConvCalo - kTestBeamv3 + shifting MC
  } else if( trainConfig == 63)  { // Timing diff (std is -100ns to 100ns)
    cuts.AddCut("00000113","00200009227000008250400000","1111111037032230000","302010708","0103603800000000","0153503000000000"); // 200ns
    cuts.AddCut("00000113","00200009227000008250400000","1111111057032230000","302010708","0103603800000000","0153503000000000"); // 50ns
    cuts.AddCut("00000113","00200009227000008250400000","1111111077032230000","302010708","0103603800000000","0153503000000000"); // 30ns
    cuts.AddCut("00000113","00200009227000008250400000","1111111097032230000","302010708","0103603800000000","0153503000000000"); // 20-25ns
  } else if( trainConfig == 64)  { // TrackMatching
    cuts.AddCut("00000113","00200009227000008250400000","1111111046032230000","302010708","0103603800000000","0153503000000000");
    cuts.AddCut("00000113","00200009227000008250400000","1111111048032230000","302010708","0103603800000000","0153503000000000");
    cuts.AddCut("00000113","00200009227000008250400000","1111111049032230000","302010708","0103603800000000","0153503000000000");
    cuts.AddCut("00000113","00200009227000008250400000","111111104a032230000","302010708","0103603800000000","0153503000000000");
    cuts.AddCut("00000113","00200009227000008250400000","111111104b032230000","302010708","0103603800000000","0153503000000000");
    cuts.AddCut("00000113","00200009227000008250400000","1111111043032230000","302010708","0103603800000000","0153503000000000"); //non pt dependent
  } else if( trainConfig == 65)  { // MinEnergy (of cluster)
    cuts.AddCut("00000113","00200009227000008250400000","1111111047022230000","302010708","0103603800000000","0153503000000000"); // 0.6
    cuts.AddCut("00000113","00200009227000008250400000","1111111047042230000","302010708","0103603800000000","0153503000000000"); // 0.8
    cuts.AddCut("00000113","00200009227000008250400000","1111111047052230000","302010708","0103603800000000","0153503000000000"); // 0.9
  } else if( trainConfig == 66)  { // MinNCells
    cuts.AddCut("00000113","00200009227000008250400000","1111111047031230000","302010708","0103603800000000","0153503000000000"); // 1
    cuts.AddCut("00000113","00200009227000008250400000","1111111047033230000","302010708","0103603800000000","0153503000000000"); // 3
  } else if( trainConfig == 67)  { // MinMaxM02
    cuts.AddCut("00000113","00200009227000008250400000","1111111047032330000","302010708","0103603800000000","0153503000000000"); // 0.2 - 0.5
    cuts.AddCut("00000113","00200009227000008250400000","1111111047032130000","302010708","0103603800000000","0153503000000000"); // 0.002 - 0.5
    cuts.AddCut("00000113","00200009227000008250400000","1111111047032240000","302010708","0103603800000000","0153503000000000"); // 0.1 - 0.4
    cuts.AddCut("00000113","00200009227000008250400000","1111111047032220000","302010708","0103603800000000","0153503000000000"); // 0.1 - 0.7

    // *************Variations in AliPrimaryPionCuts******************
  } else if( trainConfig == 68)  { // ClsTPCCut
    cuts.AddCut("00000113","00200009227000008250400000","1111111047032230000","301010708","0103603800000000","0153503000000000"); // 70
    cuts.AddCut("00000113","00200009227000008250400000","1111111047032230000","303010708","0103603800000000","0153503000000000"); // 100
    cuts.AddCut("00000113","00200009227000008250400000","1111111047032230000","305010708","0103603800000000","0153503000000000"); // 35% of findable clusters
    cuts.AddCut("00000113","00200009227000008250400000","1111111047032230000","306010708","0103603800000000","0153503000000000"); // 60% of findable clusters
    cuts.AddCut("00000113","00200009227000008250400000","1111111047032230000","30b010708","0103603800000000","0153503000000000"); // PHOS public note
  } else if ( trainConfig == 69) { // DCACut
    cuts.AddCut("00000113","00200009227000008250400000","1111111047032230000","302110708","0103603800000000","0153503000000000");  // XYPtDep("0.0182+0.0350/pt^1.01");
    cuts.AddCut("00000113","00200009227000008250400000","1111111047032230000","302210708","0103603800000000","0153503000000000");  // z=2cm xy=1cm
    cuts.AddCut("00000113","00200009227000008250400000","1111111047032230000","302310708","0103603800000000","0153503000000000");  // z=3cm XYPtDep("0.0182+0.0350/pt^1.01");
    cuts.AddCut("00000113","00200009227000008250400000","1111111047032230000","302410708","0103603800000000","0153503000000000");  // z=3cm xy=0.5
  } else if ( trainConfig == 70) { // pT Cut
    cuts.AddCut("00000113","00200009227000008250400000","1111111047032230000","302000708","0103603800000000","0153503000000000");  // pt>0.075
    cuts.AddCut("00000113","00200009227000008250400000","1111111047032230000","302020708","0103603800000000","0153503000000000");  // pt>0.125
    cuts.AddCut("00000113","00200009227000008250400000","1111111047032230000","302030708","0103603800000000","0153503000000000");  // pt>0.15
    cuts.AddCut("00000113","00200009227000008250400000","1111111047032230000","302040708","0103603800000000","0153503000000000");  // pt>0.4
  } else if ( trainConfig == 71) { // TPCdEdxCutPion
    cuts.AddCut("00000113","00200009227000008250400000","1111111047032230000","302010508","0103603800000000","0153503000000000"); // -4,4
    cuts.AddCut("00000113","00200009227000008250400000","1111111047032230000","302010808","0103603800000000","0153503000000000"); // -2,3
    cuts.AddCut("00000113","00200009227000008250400000","1111111047032230000","302010208","0103603800000000","0153503000000000"); // -6,7
    cuts.AddCut("00000113","00200009227000008250400000","1111111047032230000","302010308","0103603800000000","0153503000000000"); // -5,5
    cuts.AddCut("00000113","00200009227000008250400000","1111111047032230000","302010408","0103603800000000","0153503000000000"); // -4,5
  } else if ( trainConfig == 72) { // MassCut
    cuts.AddCut("00000113","00200009227000008250400000","1111111047032230000","302010707","0103603800000000","0153503000000000"); // 0.7 GeV
    cuts.AddCut("00000113","00200009227000008250400000","1111111047032230000","302010706","0103603800000000","0153503000000000"); // 0.65 GeV
    cuts.AddCut("00000113","00200009227000008250400000","1111111047032230000","302010701","0103603800000000","0153503000000000"); // 1 GeV
    cuts.AddCut("00000113","00200009227000008250400000","1111111047032230000","302010702","0103603800000000","0153503000000000"); // 0.75 GeV
    cuts.AddCut("00000113","00200009227000008250400000","1111111047032230000","302010704","0103603800000000","0153503000000000"); // 0.54 eta mass
    // *************Variations in AliConversionMesonCuts (NeutralPion) ******************
  } else if ( trainConfig == 73) { // RapidityMesonCut
    cuts.AddCut("00000113","00200009227000008250400000","1111111047032230000","302010708","0103503800000000","0153503000000000"); // 0.85
    cuts.AddCut("00000113","00200009227000008250400000","1111111047032230000","302010708","0103303800000000","0153503000000000"); // 0.6
    cuts.AddCut("00000113","00200009227000008250400000","1111111047032230000","302010708","0103203800000000","0153503000000000"); // 0.7
  } else if ( trainConfig == 74) { // pT
    cuts.AddCut("00000113","00200009227000008250400000","1111111047032230000","302010708","0103613800000000","0153503000000000"); // min pT pi0 = 0.4
    cuts.AddCut("00000113","00200009227000008250400000","1111111047032230000","302010708","0103623800000000","0153503000000000"); // min pT pi0 = 0.7
    cuts.AddCut("00000113","00200009227000008250400000","1111111047032230000","302010708","0103673800000000","0153503000000000"); // min pT pi0 = 0.5
  } else if ( trainConfig == 75) { // alphaMesonCut
    cuts.AddCut("00000113","00200009227000008250400000","1111111047032230000","302010708","0103607800000000","0153503000000000"); // 0-0.85
    cuts.AddCut("00000113","00200009227000008250400000","1111111047032230000","302010708","0103605800000000","0153503000000000"); // 0-0.75
    cuts.AddCut("00000113","00200009227000008250400000","1111111047032230000","302010708","0103600800000000","0153503000000000"); // 0-0.7
    cuts.AddCut("00000113","00200009227000008250400000","1111111047032230000","302010708","0103604800000000","0153503000000000"); // 0-0.65
  } else if ( trainConfig == 76) { // selectionWindow
    cuts.AddCut("00000113","00200009227000008250400000","1111111047032230000","302010708","0103603900000000","0153503000000000"); // 0.11 - 0.155
    cuts.AddCut("00000113","00200009227000008250400000","1111111047032230000","302010708","0103603300000000","0153503000000000"); // 0.12 - 0.145
    cuts.AddCut("00000113","00200009227000008250400000","1111111047032230000","302010708","0103603a00000000","0153503000000000"); // 0.08 - 0.145
    cuts.AddCut("00000113","00200009227000008250400000","1111111047032230000","302010708","0103603400000000","0153503000000000"); // 0.11 - 0.15

    // *************Variations in AliConversionMesonCuts (omega) ******************
  } else if ( trainConfig == 77) { // Background scheme
    cuts.AddCut("00000113","00200009227000008250400000","1111111047032230000","302010708","0103603800000000","0a53503000000000"); // likesign
    cuts.AddCut("00000113","00200009227000008250400000","1111111047032230000","302010708","0103603800000000","0b53503000000000"); // sideband right
    cuts.AddCut("00000113","00200009227000008250400000","1111111047032230000","302010708","0103603800000000","0c53503000000000"); // sideband left
    cuts.AddCut("00000113","00200009227000008250400000","1111111047032230000","302010708","0103603800000000","0d53503000000000"); // sideband both sides
  } else if ( trainConfig == 78) { // NumberBckEvents
    cuts.AddCut("00000113","00200009227000008250400000","1111111047032230000","302010708","0103603800000000","0133503000000000"); // 20
    cuts.AddCut("00000113","00200009227000008250400000","1111111047032230000","302010708","0103603800000000","0163503000000000"); // 80
  } else if ( trainConfig == 79) { // RapidityCut
    cuts.AddCut("00000113","00200009227000008250400000","1111111047032230000","302010708","0103603800000000","0153203000000000"); // 0.7
    cuts.AddCut("00000113","00200009227000008250400000","1111111047032230000","302010708","0103603800000000","0153003000000000"); // 1.35
    cuts.AddCut("00000113","00200009227000008250400000","1111111047032230000","302010708","0103603800000000","0153303000000000"); // 0.6
  } else if ( trainConfig == 80) { // AlphaCut
    cuts.AddCut("00000113","00200009227000008250400000","1111111047032230000","302010708","0103603800000000","0153507000000000"); // 0-0.85
    cuts.AddCut("00000113","00200009227000008250400000","1111111047032230000","302010708","0103603800000000","0153505000000000"); // 0-0.75
    cuts.AddCut("00000113","00200009227000008250400000","1111111047032230000","302010708","0103603800000000","0153500000000000"); // 0-0.7

    // 8TeV
  } else if( trainConfig == 101 ) {
    // closing charged pion cuts, minimum TPC cluster = 80, TPC dEdx pi = \pm 3 sigma, pi+pi- mass Cut at 0.65, min pt charged pi = 100 MeV
    // closing neural pion cuts, 0.1 < M_gamma,gamma < 0.15
    cuts.AddCut("00010113","00200009327000008250400000","1111111067032230000","002010706","0163103400000010","0163503000000000");
    cuts.AddCut("00052113","00200009327000008250400000","1111111067032230000","002010706","0163103400000010","0163503000000000");
    cuts.AddCut("00081113","00200009327000008250400000","1111111067032230000","002010706","0163103400000010","0163503000000000");
  } else if( trainConfig == 102 ) {
    // closing charged pion cuts, minimum TPC cluster = 80, TPC dEdx pi = \pm 3 sigma, pi+pi- mass Cut at 0.65/0.7/0.75, min pt charged pi = 100 MeV
    // closing neural pion cuts, 0.1 < M_gamma,gamma < 0.15
    cuts.AddCut("00010113","00200009327000008250400000","1111111067032230000","002010706","0163103400000010","0163503000000000");
    cuts.AddCut("00010113","00200009327000008250400000","1111111067032230000","002010707","0163103400000010","0163503000000000");
    cuts.AddCut("00010113","00200009327000008250400000","1111111067032230000","002010702","0163103400000010","0163503000000000");
  } else if( trainConfig == 103 ) {
    // closing charged pion cuts, minimum TPC cluster = 80, TPC dEdx pi = \pm 3 sigma, pi+pi- mass Cut at 0.65/0.7/0.75, min pt charged pi = 100 MeV
    // closing neural pion cuts, 0.1 < M_gamma,gamma < 0.15
    cuts.AddCut("00052113","00200009327000008250400000","1111111067032230000","002010706","0163103400000010","0163503000000000");
    cuts.AddCut("00052113","00200009327000008250400000","1111111067032230000","002010707","0163103400000010","0163503000000000");
    cuts.AddCut("00052113","00200009327000008250400000","1111111067032230000","002010702","0163103400000010","0163503000000000");
  } else if( trainConfig == 104 ) {
    // eta < 0.9
    // closing charged pion cuts, minimum TPC cluster = 80, TPC dEdx pi = \pm 3 sigma, pi+pi- mass cut of 0.65, min pt charged pi = 100 MeV
    // closing neural pion cuts, 0.1 < M_gamma,gamma < 0.15
    // maxChi2 per cluster TPC <4, require TPC refit, DCA XY pT dependend 0.0182+0.0350/pt^1.01, DCA_Z = 3.0
    cuts.AddCut("00010113","00200009327000008250400000","1111111067032230000","302010706","0103503400000000","0153503000000000"); // only cuts above
    cuts.AddCut("00010113","00200009327000008250400000","1111111067032230000","30a010706","0103503400000000","0153503000000000"); // above + fMinClsTPC=80.+ fChi2PerClsTPC=4 + fRequireTPCRefit=kTRUE;
    cuts.AddCut("00010113","00200009327000008250400000","1111111067032230000","302310706","0103503400000000","0153503000000000"); // above + DCA pT dependent 0.0182+0.0350/pt^1.01 + DCA_Z < 3.0
    cuts.AddCut("00010113","00200009327000008250400000","1111111067032230000","302030706","0103503400000000","0153503000000000"); // above + pTmin=0.15
    cuts.AddCut("00010113","00200009327000008250400000","1111111067032230000","30a330706","0103503400000000","0153503000000000"); // all of the above
  } else if( trainConfig == 105 ) {
    // eta < 0.9
    // closing charged pion cuts, minimum TPC cluster = 80, TPC dEdx pi = \pm 3 sigma, pi+pi- mass cut of 0.65, min pt charged pi = 100 MeV
    // closing neural pion cuts, 0.1 < M_gamma,gamma < 0.15
    // maxChi2 per cluster TPC <4, require TPC refit, DCA XY pT dependend 0.0182+0.0350/pt^1.01, DCA_Z = 3.0
    cuts.AddCut("00052113","00200009327000008250400000","1111111067032230000","302010706","0103503400000000","0153503000000000"); // only cuts above
    cuts.AddCut("00052113","00200009327000008250400000","1111111067032230000","30a010706","0103503400000000","0153503000000000"); // above + fMinClsTPC=80.+ fChi2PerClsTPC=4 + fRequireTPCRefit=kTRUE;
    cuts.AddCut("00052113","00200009327000008250400000","1111111067032230000","302310706","0103503400000000","0153503000000000"); // above + DCA pT dependent 0.0182+0.0350/pt^1.01 + DCA_Z < 3.0
    cuts.AddCut("00052113","00200009327000008250400000","1111111067032230000","302030706","0103503400000000","0153503000000000"); // above + pTmin=0.15
    cuts.AddCut("00052113","00200009327000008250400000","1111111067032230000","30a330706","0103503400000000","0153503000000000"); // all of the above
  } else if( trainConfig == 106 ) {
    // eta < 0.9
    // closing charged pion cuts, minimum TPC cluster = 80, TPC dEdx pi = \pm 3 sigma, pi+pi- mass cut of 0.65, min pt charged pi = 100 MeV
    // closing neural pion cuts, 0.1 < M_gamma,gamma < 0.15
    // maxChi2 per cluster TPC <4, require TPC refit, DCA XY pT dependend 0.0182+0.0350/pt^1.01, DCA_Z = 3.0
    cuts.AddCut("00081113","00200009327000008250400000","1111111067032230000","302010706","0103503400000000","0153503000000000"); // only cuts above
    cuts.AddCut("00081113","00200009327000008250400000","1111111067032230000","30a010706","0103503400000000","0153503000000000"); // above + fMinClsTPC=80.+ fChi2PerClsTPC=4 + fRequireTPCRefit=kTRUE;
    cuts.AddCut("00081113","00200009327000008250400000","1111111067032230000","302310706","0103503400000000","0153503000000000"); // above + DCA pT dependent 0.0182+0.0350/pt^1.01 + DCA_Z < 3.0
    cuts.AddCut("00081113","00200009327000008250400000","1111111067032230000","302030706","0103503400000000","0153503000000000"); // above + pTmin=0.15
    cuts.AddCut("00081113","00200009327000008250400000","1111111067032230000","30a330706","0103503400000000","0153503000000000"); // all of the above
  // 8 TeV PHOS
  } else if( trainConfig == 110) {
    // eta < 0.9
    // closing charged pion cuts, minimum TPC cluster = 80, TPC dEdx pi = \pm 3 sigma, pi+pi- mass cut of 0.65, min pt charged pi = 100 MeV
    // closing neural pion cuts, 0.1 < M_gamma,gamma < 0.15
    // maxChi2 per cluster TPC <4, require TPC refit, DCA XY pT dependend 0.0182+0.0350/pt^1.01, DCA_Z = 3.0
    cuts.AddCut("00010113","00200009327000008250400000","2444400043013300000","302010706","0103503400000000","0153503000000000"); // only cuts above
    cuts.AddCut("00010113","00200009327000008250400000","2444400043013300000","30a010706","0103503400000000","0153503000000000"); // above + fMinClsTPC=80.+ fChi2PerClsTPC=4 + fRequireTPCRefit=kTRUE;
    cuts.AddCut("00010113","00200009327000008250400000","2444400043013300000","302310706","0103503400000000","0153503000000000"); // above + DCA pT dependent 0.0182+0.0350/pt^1.01 + DCA_Z < 3.0
    cuts.AddCut("00010113","00200009327000008250400000","2444400043013300000","302030706","0103503400000000","0153503000000000"); // above + pTmin=0.15
    cuts.AddCut("00010113","00200009327000008250400000","2444400043013300000","30a330706","0103503400000000","0153503000000000"); // all of the above
  } else if( trainConfig == 111) {
    // eta < 0.9
    // closing charged pion cuts, minimum TPC cluster = 80, TPC dEdx pi = \pm 3 sigma, pi+pi- mass cut of 0.65, min pt charged pi = 100 MeV
    // closing neural pion cuts, 0.1 < M_gamma,gamma < 0.15
    // maxChi2 per cluster TPC <4, require TPC refit, DCA XY pT dependend 0.0182+0.0350/pt^1.01, DCA_Z = 3.0
    cuts.AddCut("00062113","00200009327000008250400000","2444400043013300000","302010706","0103503400000000","0153503000000000"); // only cuts above
    cuts.AddCut("00062113","00200009327000008250400000","2444400043013300000","30a010706","0103503400000000","0153503000000000"); // above + fMinClsTPC=80.+ fChi2PerClsTPC=4 + fRequireTPCRefit=kTRUE;
    cuts.AddCut("00062113","00200009327000008250400000","2444400043013300000","302310706","0103503400000000","0153503000000000"); // above + DCA pT dependent 0.0182+0.0350/pt^1.01 + DCA_Z < 3.0
    cuts.AddCut("00062113","00200009327000008250400000","2444400043013300000","302030706","0103503400000000","0153503000000000"); // above + pTmin=0.15
    cuts.AddCut("00062113","00200009327000008250400000","2444400043013300000","30a330706","0103503400000000","0153503000000000"); // all of the above

  // 13 TeV
  } else if( trainConfig == 201 ) {
    // eta < 0.9
    // closing charged pion cuts, minimum TPC cluster = 80, TPC dEdx pi = \pm 3 sigma, pi+pi- mass cut of 0.65, min pt charged pi = 100 MeV
    // closing neural pion cuts, 0.1 < M_gamma,gamma < 0.15
    // maxChi2 per cluster TPC <4, require TPC refit, DCA XY pT dependend 0.0182+0.0350/pt^1.01, DCA_Z = 3.0
    // timing cluster cut open
    cuts.AddCut("00010113","00200009327000008250400000","1111111017032230000","30a330706","0103503400000000","0153503000000000"); // all of the above
    cuts.AddCut("00052113","00200009327000008250400000","1111111017032230000","30a330706","0103503400000000","0153503000000000"); // all of the above
    cuts.AddCut("00083113","00200009327000008250400000","1111111017032230000","30a330706","0103503400000000","0153503000000000"); // all of the above
    cuts.AddCut("00085113","00200009327000008250400000","1111111017032230000","30a330706","0103503400000000","0153503000000000"); // all of the above
  } else if( trainConfig == 202 ) {
    // eta < 0.9
    // closing charged pion cuts, minimum TPC cluster = 80, TPC dEdx pi = \pm 3 sigma, pi+pi- mass cut of 0.85, min pt charged pi = 100 MeV
    // closing neural pion cuts, 0.110 < M_gamma,gamma < 0.155
    // maxChi2 per cluster TPC <4, require TPC refit, DCA XY pT dependend 0.0182+0.0350/pt^1.01, DCA_Z = 3.0
    cuts.AddCut("00010113","00200009327000008250400000","1111100047032230000","302010708","0103503900000000","0153503000000000"); // normal event mixing; Triggers: V0AND
    cuts.AddCut("00052113","00200009327000008250400000","1111100047032230000","302010708","0103503900000000","0153503000000000"); // normal event mixing; Triggers: CEMC7: V0AND and EMCAL fired
    cuts.AddCut("00083113","00200009327000008250400000","1111100047032230000","302010708","0103503900000000","0153503000000000"); // normal event mixing; Triggers: 7EG1 - CINT7 EG1
    cuts.AddCut("00085113","00200009327000008250400000","1111100047032230000","302010708","0103503900000000","0153503000000000"); // normal event mixing; Triggers: 7EG2 - CINT7 EG2
  // 13 TeV PHOS
  } else if( trainConfig == 203) {
    //As above (201) but only MB
    cuts.AddCut("00010113","00200009327000008250400000","1111111017032230000","30a330706","0103503400000000","0153503000000000");
  } else if( trainConfig == 204) {
    //As above (202) but only MB
    cuts.AddCut("00010113","00200009327000008250400000","1111100047032230000","302010708","0103503900000000","0153503000000000");
  } else if( trainConfig == 210) {
    // eta < 0.9
    // closing charged pion cuts, minimum TPC cluster = 80, TPC dEdx pi = \pm 3 sigma, pi+pi- mass cut of 0.65, min pt charged pi = 100 MeV
    // closing neural pion cuts, 0.1 < M_gamma,gamma < 0.15
    // maxChi2 per cluster TPC <4, require TPC refit, DCA XY pT dependend 0.0182+0.0350/pt^1.01, DCA_Z = 3.0
    // timing cluster cut open
    cuts.AddCut("00010113","00200009327000008250400000","2444400013013300000","30a330706","0103503400000000","0153503000000000"); // all of the above
    cuts.AddCut("00062113","00200009327000008250400000","2444400013013300000","30a330706","0103503400000000","0153503000000000"); // all of the above
  } else if( trainConfig == 211) {
    // eta < 0.9
    // closing charged pion cuts, minimum TPC cluster = 80, TPC dEdx pi = \pm 3 sigma, pi+pi- mass cut of 0.85, min pt charged pi = 100 MeV
    // closing neural pion cuts, 0.120 < M_gamma,gamma < 0.150
    // maxChi2 per cluster TPC <4, require TPC refit, DCA XY pT dependend 0.0182+0.0350/pt^1.01, DCA_Z = 3.0
    cuts.AddCut("00010113","00200009327000008250400000","2444400043013300000","302010708","0103503600000000","0153503000000000"); // normal event mixing; Triggers: V0AND
    cuts.AddCut("00062113","00200009327000008250400000","2444400043013300000","302010708","0103503600000000","0153503000000000"); // normal event mixing; Triggers: V0AND and EMCal fired
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
    TString TrackMatcherName = Form("CaloTrackMatcher_%s_%i",caloCutPos.Data(),trackMatcherRunningMode);
    if( !(AliCaloTrackMatcher*)mgr->GetTask(TrackMatcherName.Data()) ){
      AliCaloTrackMatcher* fTrackMatcher = new AliCaloTrackMatcher(TrackMatcherName.Data(),caloCutPos.Atoi(),trackMatcherRunningMode);
      fTrackMatcher->SetV0ReaderName(V0ReaderName);
      mgr->AddTask(fTrackMatcher);
      mgr->ConnectInput(fTrackMatcher,0,cinput);
    }

    analysisEventCuts[i] = new AliConvEventCuts();
    analysisEventCuts[i]->SetV0ReaderName(V0ReaderName);
    if(runLightOutput>0) analysisEventCuts[i]->SetLightOutput(kTRUE);
    analysisEventCuts[i]->InitializeCutsFromCutString((cuts.GetEventCut(i)).Data());
    if (periodNameV0Reader.CompareTo("") != 0) analysisEventCuts[i]->SetPeriodEnum(periodNameV0Reader);
    EventCutList->Add(analysisEventCuts[i]);
    analysisEventCuts[i]->SetFillCutHistograms("",kFALSE);

    analysisCuts[i] = new AliConversionPhotonCuts();
    analysisCuts[i]->SetV0ReaderName(V0ReaderName);
    if(runLightOutput>0) analysisCuts[i]->SetLightOutput(kTRUE);
    if( ! analysisCuts[i]->InitializeCutsFromCutString((cuts.GetConversionCut(i)).Data()) ) {
      cout<<"ERROR: analysisCuts [" <<i<<"]"<<endl;
      return 0;
    } else {
      ConvCutList->Add(analysisCuts[i]);
      analysisCuts[i]->SetFillCutHistograms("",kFALSE);
    }

    analysisClusterCuts[i] = new AliCaloPhotonCuts();
    analysisClusterCuts[i]->SetV0ReaderName(V0ReaderName);
    if(runLightOutput>0) analysisClusterCuts[i]->SetLightOutput(kTRUE);
    analysisClusterCuts[i]->SetCaloTrackMatcherName(TrackMatcherName);
    if( ! analysisClusterCuts[i]->InitializeCutsFromCutString((cuts.GetClusterCut(i)).Data()) ) {
      cout<<"ERROR: analysisClusterCuts [" <<i<<"]"<<endl;
      return 0;
    } else {
      ClusterCutList->Add(analysisClusterCuts[i]);
      analysisClusterCuts[i]->SetFillCutHistograms("");
    }

    analysisNeutralPionCuts[i] = new AliConversionMesonCuts();
    if(runLightOutput>0) analysisNeutralPionCuts[i]->SetLightOutput(kTRUE);
    if( ! analysisNeutralPionCuts[i]->InitializeCutsFromCutString((cuts.GetNeutralPionCut(i)).Data()) ) {
      cout<<"ERROR: analysisMesonCuts [ " <<i<<" ] "<<endl;
      return 0;
    } else {
      NeutralPionCutList->Add(analysisNeutralPionCuts[i]);
      analysisNeutralPionCuts[i]->SetFillCutHistograms("");
    }

    analysisMesonCuts[i] = new AliConversionMesonCuts();
    if(runLightOutput>0) analysisMesonCuts[i]->SetLightOutput(kTRUE);
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
    if(runLightOutput>0) analysisPionCuts[i]->SetLightOutput(kTRUE);

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

  task->SetDoMesonQA(enableQAMesonTask );

  //connect containers
  AliAnalysisDataContainer *coutput =
  mgr->CreateContainer(Form("GammaConvNeutralMesonPiPlPiMiPiZero_%i_%i",neutralPionMode, trainConfig), TList::Class(),
              AliAnalysisManager::kOutputContainer,Form("GammaConvNeutralMesonPiPlPiMiPiZero_%i_%i.root",neutralPionMode, trainConfig));

  mgr->AddTask(task);
  mgr->ConnectInput(task,0,cinput);
  mgr->ConnectOutput(task,1,coutput);

  return;

}
