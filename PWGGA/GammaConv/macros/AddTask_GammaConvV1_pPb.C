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

class CutHandlerConvpPb{
  public:
    CutHandlerConvpPb(Int_t nMax=10){
      nCuts=0; nMaxCuts=nMax; validCuts = true;
      eventCutArray = new TString[nMaxCuts]; photonCutArray = new TString[nMaxCuts]; mesonCutArray = new TString[nMaxCuts]; clusterCutArray = new TString[nMaxCuts];
      for(Int_t i=0; i<nMaxCuts; i++) {eventCutArray[i] = ""; photonCutArray[i] = ""; mesonCutArray[i] = ""; clusterCutArray[i] = "";}
    }

    void AddCut(TString eventCut, TString photonCut, TString mesonCut){
      if(nCuts>=nMaxCuts) {cout << "ERROR in CutHandlerConvpPb: Exceeded maximum number of cuts!" << endl; validCuts = false; return;}
      if( eventCut.Length()!=8 || photonCut.Length()!=26 || mesonCut.Length()!=16 ) {cout << "ERROR in CutHandlerConvpPb: Incorrect length of cut string!" << endl; validCuts = false; return;}
      eventCutArray[nCuts]=eventCut; photonCutArray[nCuts]=photonCut; mesonCutArray[nCuts]=mesonCut; clusterCutArray[nCuts]="";
      nCuts++;
      return;
    }
    void AddCut(TString eventCut, TString photonCut, TString mesonCut, TString clusterCut){
      if(nCuts>=nMaxCuts) {cout << "ERROR in CutHandlerConvpPb: Exceeded maximum number of cuts!" << endl; validCuts = false; return;}
      if( eventCut.Length()!=8 || photonCut.Length()!=26 || mesonCut.Length()!=16 || clusterCut.Length()!=19 ) {cout << "ERROR in CutHandlerConvpPb: Incorrect length of cut string!" << endl; validCuts = false; return;}
      eventCutArray[nCuts]=eventCut; photonCutArray[nCuts]=photonCut; mesonCutArray[nCuts]=mesonCut; clusterCutArray[nCuts]=clusterCut;
      nCuts++;
      return;
    }

    Bool_t AreValid(){return validCuts;}
    Int_t GetNCuts(){if(validCuts) return nCuts; else return 0;}
    TString GetEventCut(Int_t i){if(validCuts&&i<nMaxCuts&&i>=0) return eventCutArray[i]; else{cout << "ERROR in CutHandlerConvpPb: GetEventCut wrong index i" << endl;return "";}}
    TString GetPhotonCut(Int_t i){if(validCuts&&i<nMaxCuts&&i>=0) return photonCutArray[i]; else {cout << "ERROR in CutHandlerConvpPb: GetPhotonCut wrong index i" << endl;return "";}}
    TString GetClusterCut(Int_t i){if(validCuts&&i<nMaxCuts&&i>=0) return clusterCutArray[i]; else {cout << "ERROR in CutHandlerConvpPb: GetClusterCut wrong index i" << endl;return "";}}
    TString GetMesonCut(Int_t i){if(validCuts&&i<nMaxCuts&&i>=0) return mesonCutArray[i]; else {cout << "ERROR in CutHandlerConvpPb: GetMesonCut wrong index i" << endl;return "";}}
  private:
    Bool_t validCuts;
    Int_t nCuts; Int_t nMaxCuts;
    TString* eventCutArray;
    TString* photonCutArray;
    TString* mesonCutArray;
    TString* clusterCutArray;
};

void AddTask_GammaConvV1_pPb(   Int_t     trainConfig                   = 1,                              // change different set of cuts
                                Int_t     isMC                          = 0,                              // run MC
                                Int_t     enableQAMesonTask             = 0,                              // enable QA in AliAnalysisTaskGammaConvV1
                                Int_t     enableQAPhotonTask            = 0,                              // enable additional QA task
                                TString   fileNameInputForWeighting     = "MCSpectraInput.root",          // path to file for weigting input
                                Int_t     doWeightingPart               = 0,                              // enable Weighting
                                TString   generatorName                 = "DPMJET",                       // generator Name
                                TString   cutnumberAODBranch            = "800000006008400000150000000",  // cutnumber for AOD branch
                                Bool_t    enableV0findingEffi           = kFALSE,                         // enables V0finding efficiency histograms
                                Bool_t    enablePlotVsCentrality        = kFALSE,
                                Bool_t    enableTriggerMimicking        = kFALSE,                         // enable trigger mimicking
                                Bool_t    enableTriggerOverlapRej       = kFALSE,                         // enable trigger overlap rejection
                                Float_t   maxFacPtHard                  = 3.,                             // maximum factor between hardest jet and ptHard generated
                                TString   periodNameV0Reader            = "",
                                Bool_t    doMultiplicityWeighting       = kFALSE,                          //
                                TString   fileNameInputForMultWeighing  = "Multiplicity.root",             //
                                TString   periodNameAnchor              = "",                              //
                                Bool_t    runLightOutput                = kFALSE,                         // switch to run light output (only essential histograms for afterburner)
                                Bool_t    runTHnSparse                  = kTRUE,                          // switch on THNsparse
                                Int_t     enableMatBudWeightsPi0        = 0,                              // 1 = three radial bins, 2 = 10 radial bins
                                TString   filenameMatBudWeights         = "MCInputFileMaterialBudgetWeights.root",
                                Int_t     debugLevel                    = 0,                              // introducing debug levels for grid running
				Bool_t    enableElecDeDxPostCalibration = kFALSE,
				TString   fileNameElecDeDxPostCalibration = "dEdxCorrectionMap_Period_Pass.root",
                                TString   additionalTrainConfig         = "0"                             // additional counter for trainconfig, this has to be always the last parameter
                          ) {


  Int_t trackMatcherRunningMode = 0; // CaloTrackMatcher running mode
  TString corrTaskSetting = ""; // select which correction task setting to use
  //parse additionalTrainConfig flag
  TObjArray *rAddConfigArr = additionalTrainConfig.Tokenize("_");
  if(rAddConfigArr->GetEntries()<1){cout << "ERROR: AddTask_GammaConvV1_pPb during parsing of additionalTrainConfig String '" << additionalTrainConfig.Data() << "'" << endl; return;}
  TObjString* rAdditionalTrainConfig;
  for(Int_t i = 0; i<rAddConfigArr->GetEntries() ; i++){
    if(i==0) rAdditionalTrainConfig = (TObjString*)rAddConfigArr->At(i);
    else{
      TObjString* temp = (TObjString*) rAddConfigArr->At(i);
      TString tempStr = temp->GetString();
      cout<< tempStr.Data()<<endl;

      if(tempStr.Contains("MaterialBudgetWeights") && enableMatBudWeightsPi0 > 0){
        TObjArray *fileNameMatBudWeightsArr = filenameMatBudWeights.Tokenize("/");
        if(fileNameMatBudWeightsArr->GetEntries()<1 ){cout<<"ERROR: AddTask_GammaConvV1_pPb when reading material budget weights file name" << filenameMatBudWeights.Data()<< "'" << endl; return;}
        TObjString * oldMatObjStr = (TObjString*)fileNameMatBudWeightsArr->At( fileNameMatBudWeightsArr->GetEntries()-1);
        TString  oldfileName  = oldMatObjStr->GetString();
        TString  newFileName  = Form("MCInputFile%s.root",tempStr.Data());
        cout<<newFileName.Data()<<endl;
        if( oldfileName.EqualTo(newFileName.Data()) == 0 ){
          filenameMatBudWeights.ReplaceAll(oldfileName.Data(),newFileName.Data());
          cout << "INFO: AddTask_GammaConvV1_pPb the material budget weights file has been change to " <<filenameMatBudWeights.Data()<<"'"<< endl;
        }
      } else if(tempStr.BeginsWith("CF")){
        cout << "INFO: AddTask_GammaCalo_pPb will use custom branch from Correction Framework!" << endl;
        corrTaskSetting = tempStr;
        corrTaskSetting.Replace(0,2,"");
      } else if(tempStr.BeginsWith("TM")){
        TString tempType = tempStr;
        tempType.Replace(0,2,"");
        trackMatcherRunningMode = tempType.Atoi();
        cout << Form("INFO: AddTask_GammaConvV1_pp will use running mode '%i' for the TrackMatcher!",trackMatcherRunningMode) << endl;
      }
    }
  }


  TString sAdditionalTrainConfig = rAdditionalTrainConfig->GetString();
  if (sAdditionalTrainConfig.Atoi() > 0){
    trainConfig = trainConfig + sAdditionalTrainConfig.Atoi();
    cout << "INFO: AddTask_GammaConvV1_pPb running additionalTrainConfig '" << sAdditionalTrainConfig.Atoi() << "', train config: '" << trainConfig << "'" << endl;
  }
  if(corrTaskSetting.CompareTo(""))
    cout << "corrTaskSetting: " << corrTaskSetting.Data() << endl;

  cout << endl << endl;
  cout << "************************************************************************" << endl;
  cout << "************************************************************************" << endl;
  cout << "INFO: Initializing GammaConvV1 for pPb - config: " <<  trainConfig << endl;
  cout << "************************************************************************" << endl;
  cout << "************************************************************************" << endl;

  Int_t isHeavyIon = 2;

  if (debugLevel > 0){
    cout << "enabled debugging for trainconfig: " << debugLevel << endl;
  }
  // ================== GetAnalysisManager ===============================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error(Form("AddTask_GammaConvV1_%i",trainConfig), "No analysis manager found.");
    return ;
  }

  // ================== GetInputEventHandler =============================
  AliVEventHandler *inputHandler=mgr->GetInputEventHandler();
  Bool_t isAOD          = kFALSE;
  if(inputHandler->IsA()==AliAODInputHandler::Class()) isAOD  = kTRUE;


  //========= Check whether PID Reponse is there ====
  if(!(AliPIDResponse*)mgr->GetTask("PIDResponseTask")){
    Error(Form("AddTask_GammaConvV1_%i",trainConfig), "No PID response has been initialized aborting.");
    return;
  }

  //=========  Set Cutnumber for V0Reader ================================
  TString cutnumberPhoton     = "00000008400100001500000000";
  TString cutnumberEvent      = "80000003";

  Bool_t doEtaShift = kFALSE;
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  //========= Add V0 Reader to  ANALYSIS manager if not yet existent =====
  TString V0ReaderName        = Form("V0ReaderV1_%s_%s",cutnumberEvent.Data(),cutnumberPhoton.Data());
  AliV0ReaderV1 *fV0ReaderV1  =  NULL;
  if( !(AliV0ReaderV1*)mgr->GetTask(V0ReaderName.Data()) ){
    fV0ReaderV1               = new AliV0ReaderV1(V0ReaderName.Data());
    // AOD mode
    if (isAOD) fV0ReaderV1->AliV0ReaderV1::SetDeltaAODBranchName(Form("GammaConv_%s_gamma",cutnumberAODBranch.Data()));
    if (periodNameV0Reader.CompareTo("") != 0) fV0ReaderV1->SetPeriodName(periodNameV0Reader);
    fV0ReaderV1->SetUseOwnXYZCalculation(kTRUE);
    fV0ReaderV1->SetCreateAODs(kFALSE);// AOD Output
    fV0ReaderV1->SetUseAODConversionPhoton(kTRUE);
    fV0ReaderV1->SetProduceV0FindingEfficiency(enableV0findingEffi);

    AliConvEventCuts *fEventCuts=NULL;
    if(cutnumberEvent!=""){
      fEventCuts= new AliConvEventCuts(cutnumberEvent.Data(),cutnumberEvent.Data());
      fEventCuts->SetPreSelectionCutFlag(kTRUE);
      fEventCuts->SetV0ReaderName(V0ReaderName);
      if (debugLevel > 0) fEventCuts->SetDebugLevel(debugLevel);
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
      if(trainConfig==193 || trainConfig==194 || trainConfig==195 || trainConfig==196){
        fCuts->SetDodEdxSigmaCut(kFALSE);
      }
      if(fCuts->InitializeCutsFromCutString(cutnumberPhoton.Data())){
        fV0ReaderV1->SetConversionCuts(fCuts);
        fCuts->SetFillCutHistograms("",kTRUE);
      }
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

  CutHandlerConvpPb cuts;
  Bool_t doEtaShiftIndCuts = kFALSE;
  TString stringShift = "";

  // new standard configurations  MB
  if(trainConfig == 1){
    cuts.AddCut("80010113", "00200009327000008250404000", "0162103500900000"); // new default
    cuts.AddCut("80010113", "00200009327000008250400000", "0162103500900000"); // new default, no to close
  } else if (trainConfig == 2){
    cuts.AddCut("80010113", "00200009327000008250404000", "0162103500000000"); // new default
    cuts.AddCut("80010113", "00200009327000008250400000", "0162103500000000"); // new default, no to close
  } else if (trainConfig == 3){
    cuts.AddCut("80010123", "00200009327000008250404000", "0162103500900000"); // new default
    cuts.AddCut("80010123", "00200009327000008250400000", "0162103500900000"); // new default, no to close
  } else if (trainConfig == 4){
    cuts.AddCut("80010123", "00200009327000008250404000", "0162103500000000"); // new default
    cuts.AddCut("80010123", "00200009327000008250400000", "0162103500000000"); // new default, no to close
  } else if (trainConfig == 5){ // past-future protection
    cuts.AddCut("80010213", "00200009327000008250404000", "0162103500900000"); // new default, +-2.225\mus no other interaction
    cuts.AddCut("80010213", "00200009327000008250400000", "0162103500900000"); // new default, no to close, +-2.225\mus no other interaction
    cuts.AddCut("80010513", "00200009327000008250404000", "0162103500900000"); // new default, +-1.075\mus no other interaction
    cuts.AddCut("80010513", "00200009327000008250400000", "0162103500900000"); // new default, no to close, +-1.075\mus no other interaction
  } else if (trainConfig == 6){ // past-future protection
    cuts.AddCut("80010213", "00200009327000008250404000", "0162103500000000"); // new default,  +-2.225\mus no other interaction
    cuts.AddCut("80010213", "00200009327000008250400000", "0162103500000000"); // new default, no to close, +-2.225\mus no other interaction
    cuts.AddCut("80010513", "00200009327000008250404000", "0162103500000000"); // new default,  +-1.075\mus no other interaction
    cuts.AddCut("80010513", "00200009327000008250400000", "0162103500000000"); // new default, no to close,  +-1.075\mus no other interaction
  } else if (trainConfig == 7){
    cuts.AddCut("80010123", "00200009327000008250404000", "0162103500900000"); // new default
    cuts.AddCut("80010123", "00200009327000008250400000", "0162103500900000"); // new default, no to close
  } else if (trainConfig == 8){
    cuts.AddCut("80010123", "00200009327000008250404000", "0162103500000000"); // new default
    cuts.AddCut("80010123", "00200009327000008250400000", "0162103500000000"); // new default, no to close

  // default cut all cents without smearing and to close V0
  } else if (trainConfig == 10) {
    cuts.AddCut("80210113", "00200009327000008250400000", "0162103500000000"); // 0-20%
    cuts.AddCut("82410113", "00200009327000008250400000", "0162103500000000"); // 20-40%
    cuts.AddCut("84610113", "00200009327000008250400000", "0162103500000000"); // 40-60%
    cuts.AddCut("86010113", "00200009327000008250400000", "0162103500000000"); // 60-100%
  } else if (trainConfig == 11) { // past future protection: +-2.225\mus no other interaction
    cuts.AddCut("80210213", "00200009327000008250400000", "0162103500000000"); // 0-20%
    cuts.AddCut("82410213", "00200009327000008250400000", "0162103500000000"); // 20-40%
    cuts.AddCut("84610213", "00200009327000008250400000", "0162103500000000"); // 40-60%
    cuts.AddCut("86010213", "00200009327000008250400000", "0162103500000000"); // 60-100%
  } else if (trainConfig == 12) { // past future protection: +-1.075\mus no other interaction
    cuts.AddCut("80210513", "00200009327000008250400000", "0162103500000000"); // 0-20%
    cuts.AddCut("82410513", "00200009327000008250400000", "0162103500000000"); // 20-40%
    cuts.AddCut("84610513", "00200009327000008250400000", "0162103500000000"); // 40-60%
    cuts.AddCut("86010513", "00200009327000008250400000", "0162103500000000"); // 60-100%
  } else if (trainConfig == 13) {
    cuts.AddCut("90010113", "00200009327000008250400000", "0162103500900000"); // new default, no to close
    cuts.AddCut("90210113", "00200009327000008250400000", "0162103500900000"); // new default, no to close
    cuts.AddCut("92410113", "00200009327000008250400000", "0162103500900000"); // new default, no to close
    cuts.AddCut("94610113", "00200009327000008250400000", "0162103500900000"); // new default, no to close
    cuts.AddCut("96010113", "00200009327000008250400000", "0162103500900000"); // new default, no to close
  } else if (trainConfig == 14) {
    cuts.AddCut("e0010113", "00200009327000008250400000", "0162103500900000"); // new default, no to close
    cuts.AddCut("e0210113", "00200009327000008250400000", "0162103500900000"); // new default, no to close
    cuts.AddCut("e2410113", "00200009327000008250400000", "0162103500900000"); // new default, no to close
    cuts.AddCut("e4610113", "00200009327000008250400000", "0162103500900000"); // new default, no to close
    cuts.AddCut("e6010113", "00200009327000008250400000", "0162103500900000"); // new default, no to close

  // new standard configurations 0-20
  } else if (trainConfig == 20){
    cuts.AddCut("80210113", "00200009327000008250404000", "0162103500900000"); //new default
    cuts.AddCut("80210113", "00200009327000008250400000", "0162103500900000"); //new default, no to close
  } else if (trainConfig == 21){
    cuts.AddCut("80210123", "00200009327000008250404000", "0162103500900000"); //new default
    cuts.AddCut("80210123", "00200009327000008250400000", "0162103500900000"); //new default, no to close
  } else if (trainConfig == 22){
    cuts.AddCut("80210123", "00200009327000008250404000", "0162103500000000"); //new default
    cuts.AddCut("80210123", "00200009327000008250400000", "0162103500000000"); //new default, no to close

  // new standard configurations 20-40
  } else if (trainConfig == 30){
    cuts.AddCut("82410113", "00200009327000008250404000", "0162103500900000"); //new default
    cuts.AddCut("82410113", "00200009327000008250400000", "0162103500900000"); //new default, no to close
  } else if (trainConfig == 31){
    cuts.AddCut("82410123", "00200009327000008250404000", "0162103500900000"); //new default
    cuts.AddCut("82410123", "00200009327000008250400000", "0162103500900000"); //new default, no to close
  } else if (trainConfig == 32){
    cuts.AddCut("82410123", "00200009327000008250404000", "0162103500000000"); //new default
    cuts.AddCut("82410123", "00200009327000008250400000", "0162103500000000"); //new default, no to close

  // new standard configurations 40-60
  } else if (trainConfig == 40){
    cuts.AddCut("84610113", "00200009327000008250404000", "0162103500900000"); //new default
    cuts.AddCut("84610113", "00200009327000008250400000", "0162103500900000"); //new default, no to close
  } else if (trainConfig == 41){
    cuts.AddCut("84610123", "00200009327000008250404000", "0162103500900000"); //new default
    cuts.AddCut("84610123", "00200009327000008250400000", "0162103500900000"); //new default, no to close
  } else if (trainConfig == 42){
    cuts.AddCut("84610123", "00200009327000008250404000", "0162103500000000"); //new default
    cuts.AddCut("84610123", "00200009327000008250400000", "0162103500000000"); //new default, no to close

  // new standard configurations 60-80
  } else if (trainConfig == 50){
    cuts.AddCut("86010113", "00200009327000008250404000", "0162103500900000"); //new default
    cuts.AddCut("86010113", "00200009327000008250400000", "0162103500900000"); //new default, no to close
  } else if (trainConfig == 51){
    cuts.AddCut("86010123", "00200009327000008250404000", "0162103500900000"); //new default
    cuts.AddCut("86010123", "00200009327000008250400000", "0162103500900000"); //new default, no to close
  } else if (trainConfig == 52){
    cuts.AddCut("86010123", "00200009327000008250404000", "0162103500000000"); //new default
    cuts.AddCut("86010123", "00200009327000008250400000", "0162103500000000"); //new default, no to close

  //--------------------------------------------------------------------------
  // Systematics variations for standard ana w/o to close V0, wo smearing
  //--------------------------------------------------------------------------
  } else if (trainConfig == 100) { //default + dEdx Variation
    cuts.AddCut("80010113", "00200009327000008250400000", "0162103500000000"); // new default
    cuts.AddCut("80010113", "00200009217000008260400000", "0162103500000000"); // old standard cut
    cuts.AddCut("80010113", "00200009327000008250404000", "0162103500000000"); // new default w/ to clos v0
  } else if (trainConfig == 101) { // gamma eta & meson rap var
    cuts.AddCut("80010113", "03200009327000008250400000", "0162303500000000"); // |eta| < 0.65, |y| < 0.6
    cuts.AddCut("80010113", "04200009327000008250400000", "0162203500000000"); // |eta| < 0.75, |y| < 0.7
    cuts.AddCut("80010113", "01200009327000008250400000", "0162403500000000"); // |eta| < 0.6, |y| < 0.5
  } else if (trainConfig == 102) { // minR and single pt var
    cuts.AddCut("80010113", "00100009327000008250400000", "0162103500000000"); // minR 2.8
    cuts.AddCut("80010113", "00900009327000008250400000", "0162103500000000"); // minR 7.5
    cuts.AddCut("80010113", "00200079327000008250400000", "0162103500000000"); // single pT 0. GeV/c
    cuts.AddCut("80010113", "00200019327000008250400000", "0162103500000000"); // single pT 0.1 GeV/c
  } else if (trainConfig == 103) { // TPC cluster & edEdx var
    cuts.AddCut("80010113", "00200008327000008250400000", "0162103500000000"); // TPC Cluster 0.35
    cuts.AddCut("80010113", "00200006327000008250400000", "0162103500000000"); // TPC Cluster 0.7
    cuts.AddCut("80010113", "00200009227000008250400000", "0162103500000000"); // edEdx -4,5
    cuts.AddCut("80010113", "00200009627000008250400000", "0162103500000000"); // edEdx -2.5,4
    cuts.AddCut("80010113", "00200009127000008250400000", "0162103500000000"); // edEdx 5,5
  } else if (trainConfig == 104) { //PidEdx Variation
    cuts.AddCut("80010113", "00200009357000008250400000", "0162103500000000"); // PidEdx(2, >0.4GeV)
    cuts.AddCut("80010113", "00200009387300008250400000", "0162103500000000"); // PidEdx(2, >0.4GeV; 1>3.5GeV)
    cuts.AddCut("80010113", "00200009320000008250400000", "0162103500000000"); // PidEdx(1, >0.5GeV)
    cuts.AddCut("80010113", "00200009325000008250400000", "0162103500000000"); // PidEdx(1, >0.3GeV)
  } else if (trainConfig == 105) { //PidEdx Variation 2
    cuts.AddCut("80010113", "00200009327300008250400000", "0162103500000000"); // PidEdx(1, >0.4GeV; -10>3.5GeV)
    cuts.AddCut("80010113", "00200009326000008250400000", "0162103500000000"); // PidEdx(1, >0.25GeV)
    cuts.AddCut("80010113", "00200009326200008250400000", "0162103500000000"); // PidEdx(1, >0.25GeV; -10>4GeV)
    cuts.AddCut("80010113", "00200009327200008250400000", "0162103500000000"); // PidEdx(1, >0.4GeV; -10>4GeV)
  } else if (trainConfig == 106) { // qt & psipair variation
    cuts.AddCut("80010113", "00200009327000003250400000", "0162103500000000"); // qT 0.05 1D
    cuts.AddCut("80010113", "00200009327000009250400000", "0162103500000000"); // qT 0.03 2D
    cuts.AddCut("80010113", "00200009327000002250400000", "0162103500000000"); // qT 0.07 1D
    cuts.AddCut("80010113", "00200009327000008240400000", "0162103500000000"); // Psi Pair: 1D 0.2
  } else if (trainConfig == 107) { // PsiPair Variation
    cuts.AddCut("80010113", "00200009327000008210400000", "0162103500000000"); // Psi Pair: 1D 0.1
    cuts.AddCut("80010113", "00200009327000008150400000", "0162103500000000"); // chi2 50  2D
    cuts.AddCut("80010113", "00200009327000008850400000", "0162103500000000"); // chi2 20  2D
    cuts.AddCut("80010113", "00200009327000008260400000", "0162103500000000"); // psi pair 0.05
  } else if (trainConfig == 108) { // cos point & meson alpha variation
    cuts.AddCut("80010113", "00200009327000008250300000", "0162103500000000"); // cos pointing angle 0.75
    cuts.AddCut("80010113", "00200009327000008250600000", "0162103500000000"); // cos pointing angle 0.9
    cuts.AddCut("80010113", "00200009327000008250400000", "0162106500000000"); // alpha meson cut 0.8
  } else if (trainConfig == 109) { // BG mixing & smear variations
    cuts.AddCut("80010113", "00200009327000008250400000", "0262103500000000"); // BG track multiplicity
    cuts.AddCut("80010113", "00200009327000008250400000", "0162103500800000"); // fPSigSmearingCte=0.014;
    cuts.AddCut("80010113", "00200009327000008250400000", "0162103500900000"); // fPSigSmearingCte=0.014;

  //--------------------------------------------------------------------------
  // Systematics variations for standard ana w/o to close V0, wo smearing added signals
  //--------------------------------------------------------------------------
  } else if (trainConfig == 120) { //default + dEdx Variation added sig
    cuts.AddCut("80010123", "00200009327000008250400000", "0162103500000000"); // new default
    cuts.AddCut("80010123", "00200009217000008260400000", "0162103500000000"); // old standard cut
    cuts.AddCut("80010123", "00200009327000008250404000", "0162103500000000"); // new default w/ to close V0
  } else if (trainConfig == 121) {  // gamma eta & meson rap var added sig
    cuts.AddCut("80010123", "03200009327000008250400000", "0162303500000000"); // |eta| < 0.65, |y| < 0.6
    cuts.AddCut("80010123", "04200009327000008250400000", "0162203500000000"); // |eta| < 0.75, |y| < 0.7
    cuts.AddCut("80010123", "01200009327000008250400000", "0162403500000000"); // |eta| < 0.6, |y| < 0.5
  } else if (trainConfig == 122) { // minR and single pt var added sig
    cuts.AddCut("80010123", "00100009327000008250400000", "0162103500000000"); // minR 2.8
    cuts.AddCut("80010123", "00900009327000008250400000", "0162103500000000"); // minR 7.5
    cuts.AddCut("80010123", "00200079327000008250400000", "0162103500000000"); // single pT 0. GeV/c
    cuts.AddCut("80010123", "00200019327000008250400000", "0162103500000000"); // single pT 0.1 GeV/c
  } else if (trainConfig == 123) { // TPC cluster & edEdx var add sig
    cuts.AddCut("80010123", "00200008327000008250400000", "0162103500000000"); // TPC Cluster 0.35
    cuts.AddCut("80010123", "00200006327000008250400000", "0162103500000000"); // TPC Cluster 0.7
    cuts.AddCut("80010123", "00200009227000008250400000", "0162103500000000"); // edEdx -3,5
    cuts.AddCut("80010123", "00200009627000008250400000", "0162103500000000"); // edEdx -2.5,4
    cuts.AddCut("80010123", "00200009127000008250400000", "0162103500000000"); // edEdx -5,5
  } else if (trainConfig == 124) { //PidEdx Variation
    cuts.AddCut("80010123", "00200009357000008250400000", "0162103500000000"); // PidEdx(2, >0.4GeV)
    cuts.AddCut("80010123", "00200009387300008250400000", "0162103500000000"); // PidEdx(2, >0.4GeV; 1>3.5GeV)
    cuts.AddCut("80010123", "00200009320000008250400000", "0162103500000000"); // PidEdx(1, >0.5GeV)
    cuts.AddCut("80010123", "00200009325000008250400000", "0162103500000000"); // PidEdx(1, >0.3GeV)
  } else if (trainConfig == 125) { //PidEdx Variation 2
    cuts.AddCut("80010123", "00200009327300008250400000", "0162103500000000"); // PidEdx(1, >0.4GeV; -10>3.5GeV)
    cuts.AddCut("80010123", "00200009326000008250400000", "0162103500000000"); // PidEdx(1, >0.25GeV)
    cuts.AddCut("80010123", "00200009326200008250400000", "0162103500000000"); // PidEdx(1, >0.25GeV; -10>4GeV)
    cuts.AddCut("80010123", "00200009327200008250400000", "0162103500000000"); // PidEdx(1, >0.4GeV; -10>4GeV)
  } else if (trainConfig == 126) { // qt & psipair variation
    cuts.AddCut("80010123", "00200009327000003250400000", "0162103500000000"); // qT 0.05 1D
    cuts.AddCut("80010123", "00200009327000009250400000", "0162103500000000"); // qT 0.03 2D
    cuts.AddCut("80010123", "00200009327000002250400000", "0162103500000000"); // qT 0.07 1D
    cuts.AddCut("80010123", "00200009327000008240400000", "0162103500000000"); // Psi Pair: 1D 0.2
  } else if (trainConfig == 127) { // PsiPair Variation
    cuts.AddCut("80010123", "00200009327000008210400000", "0162103500000000"); // Psi Pair: 1D 0.1
    cuts.AddCut("80010123", "00200009327000008150400000", "0162103500000000"); // chi2 50  2D
    cuts.AddCut("80010123", "00200009327000008850400000", "0162103500000000"); // chi2 20  2D
    cuts.AddCut("80010123", "00200009327000008260400000", "0162103500000000"); // psi pair 0.05
  } else if (trainConfig == 128) { // cos point & meson alpha variation
    cuts.AddCut("80010123", "00200009327000008250300000", "0162103500000000"); // cos pointing angle 0.75
    cuts.AddCut("80010123", "00200009327000008250600000", "0162103500000000"); // cos pointing angle 0.9
    cuts.AddCut("80010123", "00200009327000008250400000", "0162106500000000"); // alpha meson cut 0.8
  } else if (trainConfig == 129) { // BG mixing & smear variations
    cuts.AddCut("80010123", "00200009327000008250400000", "0262103500000000"); // BG track multiplicity
    cuts.AddCut("80010123", "00200009327000008250400000", "0162103500800000"); // fPSigSmearingCte=0.014;
    cuts.AddCut("80010123", "00200009327000008250400000", "0162103500900000"); // fPSigSmearingCte=0.014;
    // pseudorapidity studies
    //-----------------------------------------------------------------------
  } else if (trainConfig == 130) {
    cuts.AddCut("80010113", "0a200009327000008250400000", "0162103500000000"); //Eta cut -0.9 - -0.2 and 0.2 - 0.9
    cuts.AddCut("80010113", "0b200009327000008250400000", "0162103500000000"); //Eta cut -0.9 - -0.2 and 0.2 - 0.9 with LineCut
  } else if (trainConfig == 131) {
    cuts.AddCut("80010123", "0a200009327000008250400000", "0162103500000000"); //Eta cut -0.9 - -0.2 and 0.2 - 0.9 Added Signals
    cuts.AddCut("80010123", "0b200009327000008250400000", "0162103500000000"); //Eta cut -0.9 - -0.2 and 0.2 - 0.9 Added Signals

  //--------------------------------------------------------------------------
  // purity studies (kappa cut)
  //--------------------------------------------------------------------------
  } else if (trainConfig == 150) {
    cuts.AddCut("80010113", "00200009300000008250400000", "0162103500000000"); // -3 < kappa <  5
    cuts.AddCut("80010113", "00200009500000008250400000", "0162103500000000"); // -5 < kappa < 10
    cuts.AddCut("80010113", "00200009600000008250400000", "0162103500000000"); // -3 < kappa < 10
    cuts.AddCut("80010113", "00200009700000008250400000", "0162103500000000"); //  0 < kappa < 10
  } else if (trainConfig == 152) {
    cuts.AddCut("80210113", "00200009300000008250400000", "0162103500000000"); // -3 < kappa <  5
    cuts.AddCut("80210113", "00200009500000008250400000", "0162103500000000"); // -5 < kappa < 10
    cuts.AddCut("80210113", "00200009600000008250400000", "0162103500000000"); // -3 < kappa < 10
    cuts.AddCut("80210113", "00200009700000008250400000", "0162103500000000"); //  0 < kappa < 10
  } else if (trainConfig == 153) {
    cuts.AddCut("82410113", "00200009300000008250400000", "0162103500000000"); // -3 < kappa <  5
    cuts.AddCut("82410113", "00200009500000008250400000", "0162103500000000"); // -5 < kappa < 10
    cuts.AddCut("82410113", "00200009600000008250400000", "0162103500000000"); // -3 < kappa < 10
    cuts.AddCut("82410113", "00200009700000008250400000", "0162103500000000"); //  0 < kappa < 10
  } else if (trainConfig == 154) {
    cuts.AddCut("84610113", "00200009300000008250400000", "0162103500000000"); // -3 < kappa <  5
    cuts.AddCut("84610113", "00200009500000008250400000", "0162103500000000"); // -5 < kappa < 10
    cuts.AddCut("84610113", "00200009600000008250400000", "0162103500000000"); // -3 < kappa < 10
    cuts.AddCut("84610113", "00200009700000008250400000", "0162103500000000"); //  0 < kappa < 10
  } else if (trainConfig == 155) {
    cuts.AddCut("86010113", "00200009300000008250400000", "0162103500000000"); // -3 < kappa <  5
    cuts.AddCut("86010113", "00200009500000008250400000", "0162103500000000"); // -5 < kappa < 10
    cuts.AddCut("86010113", "00200009600000008250400000", "0162103500000000"); // -3 < kappa < 10
    cuts.AddCut("86010113", "00200009700000008250400000", "0162103500000000"); //  0 < kappa < 10

  //--------------------------------------------------------------------------
  // Material weight studies
  //--------------------------------------------------------------------------
  // standard cuts
  } else if (trainConfig == 200) {
    cuts.AddCut("80010113", "00200009327000008250404000", "0162103500000000");
  } else if (trainConfig == 201) {
    cuts.AddCut("80010113", "00200009327000008250404000", "0162103500000000");
  } else if (trainConfig == 202) {
    cuts.AddCut("80010113", "00200009327000008250404000", "0162103500000000");

  // standard cuts added signals
  } else if (trainConfig == 210) {
    cuts.AddCut("80010123", "00200009327000008250404000", "0162103500000000");
  } else if (trainConfig == 211) {
    cuts.AddCut("80010123", "00200009327000008250404000", "0162103500000000");
  } else if (trainConfig == 212) {
    cuts.AddCut("80010123", "00200009327000008250404000", "0162103500000000");

  //--------------------------------------------------------------------------
  // 2016 pPb w/o past future protection
  //--------------------------------------------------------------------------
  // Min Bias
  } else if (trainConfig == 300) {
    cuts.AddCut("80010113", "00200009397302001280004000", "0162103500000000"); // Min Bias
    cuts.AddCut("80010113", "00200009397302001280004000", "0162101500000000"); // Min Bias , alpha pT dependent
  // TRD trigger
  } else if (trainConfig == 301) {
    cuts.AddCut("80047113", "00200009397302001280004000", "0162103500000000"); // TRD trigger HQU for 8TeV
    cuts.AddCut("80043113", "00200009397302001280004000", "0162103500000000"); // TRD trigger HSE for 8TeV
  // Calo triggers
  } else if (trainConfig == 302) {  // EMC triggers
    cuts.AddCut("80010113", "00200009327000008250404000", "0162103500000000", "1111100007032230000"); // Min Bias
    cuts.AddCut("80052113", "00200009327000008250404000", "0162103500000000", "1111100007032230000"); // EMC7
    cuts.AddCut("80085113", "00200009327000008250404000", "0162103500000000", "1111100007032230000"); // EG2
    cuts.AddCut("80083113", "00200009327000008250404000", "0162103500000000", "1111100007032230000"); // EG1
  } else if (trainConfig == 303) {  // PHOS triggers
    cuts.AddCut("80010113", "00200009327000008250404000", "0162103500000000", "2444400041013200000"); // MinBias
    cuts.AddCut("80052113", "00200009327000008250404000", "0162103500000000", "2444400041013200000"); // PHI7
  //
  } else if (trainConfig == 304) {
    cuts.AddCut("80110113", "00200009327000008250404000", "0162103500000000"); // 0-10
    cuts.AddCut("81210113", "00200009327000008250404000", "0162103500000000"); // 0-20
    cuts.AddCut("82410113", "00200009327000008250404000", "0162103500000000"); // 20-40
    cuts.AddCut("84610113", "00200009327000008250404000", "0162103500000000"); // 40-60
    cuts.AddCut("86810113", "00200009327000008250404000", "0162103500000000"); // 60-80
    cuts.AddCut("88010113", "00200009327000008250404000", "0162103500000000"); // 80-100
  } else if (trainConfig == 305) {
    cuts.AddCut("80010113", "00200009327000008250404000", "0162103500000000"); // new default for 5TeV
    cuts.AddCut("80210113", "00200009327000008250404000", "0162103500000000"); // 0-20
    cuts.AddCut("86010113", "00200009327000008250404000", "0162103500000000"); // 60-100
    cuts.AddCut("a0110113", "00200009327000008250404000", "0162103500000000"); // 0-5
    cuts.AddCut("a1210113", "00200009327000008250404000", "0162103500000000"); // 5-10

  } else if (trainConfig == 306){ // pPb 2013 TeV defaults
    cuts.AddCut("80110113", "00200009327000008250400000", "0162103500000000"); // 0-10
    cuts.AddCut("81210113", "00200009327000008250400000", "0162103500000000"); // 10-20
    cuts.AddCut("82410113", "00200009327000008250400000", "0162103500000000"); // 20-40
    cuts.AddCut("84610113", "00200009327000008250400000", "0162103500000000"); // 40-60
    cuts.AddCut("86810113", "00200009327000008250400000", "0162103500000000"); // 60-80
    cuts.AddCut("88010113", "00200009327000008250400000", "0162103500000000"); // 80-100
  } else if (trainConfig == 307){ // pPb 2013 TeV defaults
    cuts.AddCut("80210113", "00200009327000008250400000", "0162103500000000"); // 0-20
    cuts.AddCut("86010113", "00200009327000008250400000", "0162103500000000"); // 60-100
    cuts.AddCut("a0110113", "00200009327000008250400000", "0162103500000000"); // 0-5
    cuts.AddCut("a1210113", "00200009327000008250400000", "0162103500000000"); // 5-10

  } else if (trainConfig == 308) { // pPb 5 TeV Run2 narrow cent
    cuts.AddCut("c0110113", "00200009327000008250404000", "0162103500000000"); // 0-1
    cuts.AddCut("c0210113", "00200009327000008250404000", "0162103500000000"); // 0-2
  //--------------------------------------------------------------------------
  // 2016 pPb w/ past future protection 2.24 \mus protected
  //--------------------------------------------------------------------------
  // Min Bias
  } else if (trainConfig == 310) {
    cuts.AddCut("80010213", "00200009397302001280004000", "0162103500000000"); // Min Bias
    cuts.AddCut("80010213", "00200009397302001280004000", "0162101500000000"); // Min Bias , alpha pT dependent
  // TRD trigger
  } else if (trainConfig == 311) {
    cuts.AddCut("80047213", "00200009397302001280004000", "0162103500000000"); // TRD trigger HQU for 8TeV
    cuts.AddCut("80043213", "00200009397302001280004000", "0162103500000000"); // TRD trigger HSE for 8TeV
  // Calo triggers
  } else if (trainConfig == 312) {  // EMC triggers
    cuts.AddCut("80010213", "00200009327000008250404000", "0162103500000000", "1111100007032230000"); // Min Bias
    cuts.AddCut("80052213", "00200009327000008250404000", "0162103500000000", "1111100007032230000"); // EMC7
    cuts.AddCut("80085213", "00200009327000008250404000", "0162103500000000", "1111100007032230000"); // EG2
    cuts.AddCut("80083213", "00200009327000008250404000", "0162103500000000", "1111100007032230000"); // EG1
  } else if (trainConfig == 313) {  // PHOS triggers
    cuts.AddCut("80010213", "00200009327000008250404000", "0162103500000000", "2444400041013200000"); // MinBias
    cuts.AddCut("80052213", "00200009327000008250404000", "0162103500000000", "2444400041013200000"); // PHI7
  } else if (trainConfig == 314) {
    cuts.AddCut("80010213", "00200009327000008250404000", "0162103500000000"); // new default for 5TeV
    cuts.AddCut("80010213", "00200009327000008250404000", "0162101500000000"); // new default, alpha pT dependent for 5TeV
  //--------------------------------------------------------------------------
  // Systematics variations for direct photons pPb 5 TeV Run1 and Run2 - 0-100%
  //--------------------------------------------------------------------------
  } else if (trainConfig == 400) { //default + dEdx Variation
    cuts.AddCut("80010113", "00200009327000008250400000", "0162103500000000"); // new default
    cuts.AddCut("80010113", "00200009217000008260400000", "0162103500000000"); // old standard cut
    cuts.AddCut("80010113", "00200009327000008250404000", "0162103500000000"); // new default w/ to clos v0
  } else if (trainConfig == 401) { // gamma eta & meson rap var
    cuts.AddCut("80010113", "03200009327000008250400000", "0162303500000000"); // |eta| < 0.65, |y| < 0.6
    cuts.AddCut("80010113", "04200009327000008250400000", "0162203500000000"); // |eta| < 0.75, |y| < 0.7
    cuts.AddCut("80010113", "01200009327000008250400000", "0162403500000000"); // |eta| < 0.6, |y| < 0.5
  } else if (trainConfig == 402) { // minR and single pt var
    cuts.AddCut("80010113", "00100009327000008250400000", "0162103500000000"); // minR 2.8
    cuts.AddCut("80010113", "00900009327000008250400000", "0162103500000000"); // minR 7.5
    cuts.AddCut("80010113", "00200079327000008250400000", "0162103500000000"); // single pT 0. GeV/c
    cuts.AddCut("80010113", "00200019327000008250400000", "0162103500000000"); // single pT 0.1 GeV/c
  } else if (trainConfig == 403) { // TPC cluster & edEdx var
    cuts.AddCut("80010113", "00200008327000008250400000", "0162103500000000"); // TPC Cluster 0.35
    cuts.AddCut("80010113", "00200006327000008250400000", "0162103500000000"); // TPC Cluster 0.7
    cuts.AddCut("80010113", "00200009227000008250400000", "0162103500000000"); // edEdx -4,5
    cuts.AddCut("80010113", "00200009627000008250400000", "0162103500000000"); // edEdx -2.5,4
    cuts.AddCut("80010113", "00200009127000008250400000", "0162103500000000"); // edEdx 5,5
  } else if (trainConfig == 404) { //PidEdx Variation
    cuts.AddCut("80010113", "00200009357000008250400000", "0162103500000000"); // PidEdx(2, >0.4GeV)
    cuts.AddCut("80010113", "00200009387300008250400000", "0162103500000000"); // PidEdx(2, >0.4GeV; 1>3.5GeV)
    cuts.AddCut("80010113", "00200009320000008250400000", "0162103500000000"); // PidEdx(1, >0.5GeV)
    cuts.AddCut("80010113", "00200009325000008250400000", "0162103500000000"); // PidEdx(1, >0.3GeV)
  } else if (trainConfig == 405) { //PidEdx Variation 2
    cuts.AddCut("80010113", "00200009327300008250400000", "0162103500000000"); // PidEdx(1, >0.4GeV; -10>3.5GeV)
    cuts.AddCut("80010113", "00200009326000008250400000", "0162103500000000"); // PidEdx(1, >0.25GeV)
    cuts.AddCut("80010113", "00200009326200008250400000", "0162103500000000"); // PidEdx(1, >0.25GeV; -10>4GeV)
    cuts.AddCut("80010113", "00200009327200008250400000", "0162103500000000"); // PidEdx(1, >0.4GeV; -10>4GeV)
  } else if (trainConfig == 406) { // qt & psipair variation
    cuts.AddCut("80010113", "00200009327000003250400000", "0162103500000000"); // qT 0.05 1D
    cuts.AddCut("80010113", "00200009327000009250400000", "0162103500000000"); // qT 0.03 2D
    cuts.AddCut("80010113", "00200009327000002250400000", "0162103500000000"); // qT 0.07 1D
    cuts.AddCut("80010113", "00200009327000008240400000", "0162103500000000"); // Psi Pair: 1D 0.2
  } else if (trainConfig == 407) { // PsiPair Variation
    cuts.AddCut("80010113", "00200009327000008210400000", "0162103500000000"); // Psi Pair: 1D 0.1
    cuts.AddCut("80010113", "00200009327000008150400000", "0162103500000000"); // chi2 50  2D
    cuts.AddCut("80010113", "00200009327000008850400000", "0162103500000000"); // chi2 20  2D
    cuts.AddCut("80010113", "00200009327000008260400000", "0162103500000000"); // psi pair 0.05
  } else if (trainConfig == 408) { // cos point & meson alpha variation
    cuts.AddCut("80010113", "00200009327000008250300000", "0162103500000000"); // cos pointing angle 0.75
    cuts.AddCut("80010113", "00200009327000008250600000", "0162103500000000"); // cos pointing angle 0.9
    cuts.AddCut("80010113", "00200009327000008250400000", "0162106500000000"); // alpha meson cut 0.8
  } else if (trainConfig == 409) { // BG mixing & smear variations
    cuts.AddCut("80010113", "00200009327000008250400000", "0262103500000000"); // BG track multiplicity
    cuts.AddCut("80010113", "00200009327000008250400000", "0162103500800000"); // fPSigSmearingCte=0.014;
    cuts.AddCut("80010113", "00200009327000008250400000", "0162103500900000"); // fPSigSmearingCte=0.014;

  //--------------------------------------------------------------------------
  // Systematics variations for direct photons pPb 5 TeV Run1 and Run2 - 0-20%
  //--------------------------------------------------------------------------
  } else if (trainConfig == 410) { //default + dEdx Variation
    cuts.AddCut("80210113", "00200009327000008250400000", "0162103500000000"); // new default
    cuts.AddCut("80210113", "00200009217000008260400000", "0162103500000000"); // old standard cut
    cuts.AddCut("80210113", "00200009327000008250404000", "0162103500000000"); // new default w/ to clos v0
  } else if (trainConfig == 411) { // gamma eta & meson rap var
    cuts.AddCut("80210113", "03200009327000008250400000", "0162303500000000"); // |eta| < 0.65, |y| < 0.6
    cuts.AddCut("80210113", "04200009327000008250400000", "0162203500000000"); // |eta| < 0.75, |y| < 0.7
    cuts.AddCut("80210113", "01200009327000008250400000", "0162403500000000"); // |eta| < 0.6, |y| < 0.5
  } else if (trainConfig == 412) { // minR and single pt var
    cuts.AddCut("80210113", "00100009327000008250400000", "0162103500000000"); // minR 2.8
    cuts.AddCut("80210113", "00900009327000008250400000", "0162103500000000"); // minR 7.5
    cuts.AddCut("80210113", "00200079327000008250400000", "0162103500000000"); // single pT 0. GeV/c
    cuts.AddCut("80210113", "00200019327000008250400000", "0162103500000000"); // single pT 0.1 GeV/c
  } else if (trainConfig == 413) { // TPC cluster & edEdx var
    cuts.AddCut("80210113", "00200008327000008250400000", "0162103500000000"); // TPC Cluster 0.35
    cuts.AddCut("80210113", "00200006327000008250400000", "0162103500000000"); // TPC Cluster 0.7
    cuts.AddCut("80210113", "00200009227000008250400000", "0162103500000000"); // edEdx -4,5
    cuts.AddCut("80210113", "00200009627000008250400000", "0162103500000000"); // edEdx -2.5,4
    cuts.AddCut("80210113", "00200009127000008250400000", "0162103500000000"); // edEdx 5,5
  } else if (trainConfig == 414) { //PidEdx Variation
    cuts.AddCut("80210113", "00200009357000008250400000", "0162103500000000"); // PidEdx(2, >0.4GeV)
    cuts.AddCut("80210113", "00200009387300008250400000", "0162103500000000"); // PidEdx(2, >0.4GeV; 1>3.5GeV)
    cuts.AddCut("80210113", "00200009320000008250400000", "0162103500000000"); // PidEdx(1, >0.5GeV)
    cuts.AddCut("80210113", "00200009325000008250400000", "0162103500000000"); // PidEdx(1, >0.3GeV)
  } else if (trainConfig == 415) { //PidEdx Variation 2
    cuts.AddCut("80210113", "00200009327300008250400000", "0162103500000000"); // PidEdx(1, >0.4GeV; -10>3.5GeV)
    cuts.AddCut("80210113", "00200009326000008250400000", "0162103500000000"); // PidEdx(1, >0.25GeV)
    cuts.AddCut("80210113", "00200009326200008250400000", "0162103500000000"); // PidEdx(1, >0.25GeV; -10>4GeV)
    cuts.AddCut("80210113", "00200009327200008250400000", "0162103500000000"); // PidEdx(1, >0.4GeV; -10>4GeV)
  } else if (trainConfig == 416) { // qt & psipair variation
    cuts.AddCut("80210113", "00200009327000003250400000", "0162103500000000"); // qT 0.05 1D
    cuts.AddCut("80210113", "00200009327000009250400000", "0162103500000000"); // qT 0.03 2D
    cuts.AddCut("80210113", "00200009327000002250400000", "0162103500000000"); // qT 0.07 1D
    cuts.AddCut("80210113", "00200009327000008240400000", "0162103500000000"); // Psi Pair: 1D 0.2
  } else if (trainConfig == 417) { // PsiPair Variation
    cuts.AddCut("80210113", "00200009327000008210400000", "0162103500000000"); // Psi Pair: 1D 0.1
    cuts.AddCut("80210113", "00200009327000008150400000", "0162103500000000"); // chi2 50  2D
    cuts.AddCut("80210113", "00200009327000008850400000", "0162103500000000"); // chi2 20  2D
    cuts.AddCut("80210113", "00200009327000008260400000", "0162103500000000"); // psi pair 0.05
  } else if (trainConfig == 418) { // cos point & meson alpha variation
    cuts.AddCut("80210113", "00200009327000008250300000", "0162103500000000"); // cos pointing angle 0.75
    cuts.AddCut("80210113", "00200009327000008250600000", "0162103500000000"); // cos pointing angle 0.9
    cuts.AddCut("80210113", "00200009327000008250400000", "0162106500000000"); // alpha meson cut 0.8
  } else if (trainConfig == 419) { // BG mixing & smear variations
    cuts.AddCut("80210113", "00200009327000008250400000", "0262103500000000"); // BG track multiplicity
    cuts.AddCut("80210113", "00200009327000008250400000", "0162103500800000"); // fPSigSmearingCte=0.014;
    cuts.AddCut("80210113", "00200009327000008250400000", "0162103500900000"); // fPSigSmearingCte=0.014;

  //--------------------------------------------------------------------------
  // Systematics variations for direct photons pPb 5 TeV Run1 and Run2 - 20-40%
  //--------------------------------------------------------------------------
  } else if (trainConfig == 420) { //default + dEdx Variation
    cuts.AddCut("82410113", "00200009327000008250400000", "0162103500000000"); // new default
    cuts.AddCut("82410113", "00200009217000008260400000", "0162103500000000"); // old standard cut
    cuts.AddCut("82410113", "00200009327000008250404000", "0162103500000000"); // new default w/ to clos v0
  } else if (trainConfig == 421) { // gamma eta & meson rap var
    cuts.AddCut("82410113", "03200009327000008250400000", "0162303500000000"); // |eta| < 0.65, |y| < 0.6
    cuts.AddCut("82410113", "04200009327000008250400000", "0162203500000000"); // |eta| < 0.75, |y| < 0.7
    cuts.AddCut("82410113", "01200009327000008250400000", "0162403500000000"); // |eta| < 0.6, |y| < 0.5
  } else if (trainConfig == 422) { // minR and single pt var
    cuts.AddCut("82410113", "00100009327000008250400000", "0162103500000000"); // minR 2.8
    cuts.AddCut("82410113", "00900009327000008250400000", "0162103500000000"); // minR 7.5
    cuts.AddCut("82410113", "00200079327000008250400000", "0162103500000000"); // single pT 0. GeV/c
    cuts.AddCut("82410113", "00200019327000008250400000", "0162103500000000"); // single pT 0.1 GeV/c
  } else if (trainConfig == 423) { // TPC cluster & edEdx var
    cuts.AddCut("82410113", "00200008327000008250400000", "0162103500000000"); // TPC Cluster 0.35
    cuts.AddCut("82410113", "00200006327000008250400000", "0162103500000000"); // TPC Cluster 0.7
    cuts.AddCut("82410113", "00200009227000008250400000", "0162103500000000"); // edEdx -4,5
    cuts.AddCut("82410113", "00200009627000008250400000", "0162103500000000"); // edEdx -2.5,4
    cuts.AddCut("82410113", "00200009127000008250400000", "0162103500000000"); // edEdx 5,5
  } else if (trainConfig == 424) { //PidEdx Variation
    cuts.AddCut("82410113", "00200009357000008250400000", "0162103500000000"); // PidEdx(2, >0.4GeV)
    cuts.AddCut("82410113", "00200009387300008250400000", "0162103500000000"); // PidEdx(2, >0.4GeV; 1>3.5GeV)
    cuts.AddCut("82410113", "00200009320000008250400000", "0162103500000000"); // PidEdx(1, >0.5GeV)
    cuts.AddCut("82410113", "00200009325000008250400000", "0162103500000000"); // PidEdx(1, >0.3GeV)
  } else if (trainConfig == 425) { //PidEdx Variation 2
    cuts.AddCut("82410113", "00200009327300008250400000", "0162103500000000"); // PidEdx(1, >0.4GeV; -10>3.5GeV)
    cuts.AddCut("82410113", "00200009326000008250400000", "0162103500000000"); // PidEdx(1, >0.25GeV)
    cuts.AddCut("82410113", "00200009326200008250400000", "0162103500000000"); // PidEdx(1, >0.25GeV; -10>4GeV)
    cuts.AddCut("82410113", "00200009327200008250400000", "0162103500000000"); // PidEdx(1, >0.4GeV; -10>4GeV)
  } else if (trainConfig == 426) { // qt & psipair variation
    cuts.AddCut("82410113", "00200009327000003250400000", "0162103500000000"); // qT 0.05 1D
    cuts.AddCut("82410113", "00200009327000009250400000", "0162103500000000"); // qT 0.03 2D
    cuts.AddCut("82410113", "00200009327000002250400000", "0162103500000000"); // qT 0.07 1D
    cuts.AddCut("82410113", "00200009327000008240400000", "0162103500000000"); // Psi Pair: 1D 0.2
  } else if (trainConfig == 427) { // PsiPair Variation
    cuts.AddCut("82410113", "00200009327000008210400000", "0162103500000000"); // Psi Pair: 1D 0.1
    cuts.AddCut("82410113", "00200009327000008150400000", "0162103500000000"); // chi2 50  2D
    cuts.AddCut("82410113", "00200009327000008850400000", "0162103500000000"); // chi2 20  2D
    cuts.AddCut("82410113", "00200009327000008260400000", "0162103500000000"); // psi pair 0.05
  } else if (trainConfig == 428) { // cos point & meson alpha variation
    cuts.AddCut("82410113", "00200009327000008250300000", "0162103500000000"); // cos pointing angle 0.75
    cuts.AddCut("82410113", "00200009327000008250600000", "0162103500000000"); // cos pointing angle 0.9
    cuts.AddCut("82410113", "00200009327000008250400000", "0162106500000000"); // alpha meson cut 0.8
  } else if (trainConfig == 429) { // BG mixing & smear variations
    cuts.AddCut("82410113", "00200009327000008250400000", "0262103500000000"); // BG track multiplicity
    cuts.AddCut("82410113", "00200009327000008250400000", "0162103500800000"); // fPSigSmearingCte=0.014;
    cuts.AddCut("82410113", "00200009327000008250400000", "0162103500900000"); // fPSigSmearingCte=0.014;

  //--------------------------------------------------------------------------
  // Systematics variations for direct photons pPb 5 TeV Run1 and Run2 - 40-60%
  //--------------------------------------------------------------------------
  } else if (trainConfig == 430) { //default + dEdx Variation
    cuts.AddCut("84610113", "00200009327000008250400000", "0162103500000000"); // new default
    cuts.AddCut("84610113", "00200009217000008260400000", "0162103500000000"); // old standard cut
    cuts.AddCut("84610113", "00200009327000008250404000", "0162103500000000"); // new default w/ to clos v0
  } else if (trainConfig == 431) { // gamma eta & meson rap var
    cuts.AddCut("84610113", "03200009327000008250400000", "0162303500000000"); // |eta| < 0.65, |y| < 0.6
    cuts.AddCut("84610113", "04200009327000008250400000", "0162203500000000"); // |eta| < 0.75, |y| < 0.7
    cuts.AddCut("84610113", "01200009327000008250400000", "0162403500000000"); // |eta| < 0.6, |y| < 0.5
  } else if (trainConfig == 432) { // minR and single pt var
    cuts.AddCut("84610113", "00100009327000008250400000", "0162103500000000"); // minR 2.8
    cuts.AddCut("84610113", "00900009327000008250400000", "0162103500000000"); // minR 7.5
    cuts.AddCut("84610113", "00200079327000008250400000", "0162103500000000"); // single pT 0. GeV/c
    cuts.AddCut("84610113", "00200019327000008250400000", "0162103500000000"); // single pT 0.1 GeV/c
  } else if (trainConfig == 433) { // TPC cluster & edEdx var
    cuts.AddCut("84610113", "00200008327000008250400000", "0162103500000000"); // TPC Cluster 0.35
    cuts.AddCut("84610113", "00200006327000008250400000", "0162103500000000"); // TPC Cluster 0.7
    cuts.AddCut("84610113", "00200009227000008250400000", "0162103500000000"); // edEdx -4,5
    cuts.AddCut("84610113", "00200009627000008250400000", "0162103500000000"); // edEdx -2.5,4
    cuts.AddCut("84610113", "00200009127000008250400000", "0162103500000000"); // edEdx 5,5
  } else if (trainConfig == 434) { //PidEdx Variation
    cuts.AddCut("84610113", "00200009357000008250400000", "0162103500000000"); // PidEdx(2, >0.4GeV)
    cuts.AddCut("84610113", "00200009387300008250400000", "0162103500000000"); // PidEdx(2, >0.4GeV; 1>3.5GeV)
    cuts.AddCut("84610113", "00200009320000008250400000", "0162103500000000"); // PidEdx(1, >0.5GeV)
    cuts.AddCut("84610113", "00200009325000008250400000", "0162103500000000"); // PidEdx(1, >0.3GeV)
  } else if (trainConfig == 435) { //PidEdx Variation 2
    cuts.AddCut("84610113", "00200009327300008250400000", "0162103500000000"); // PidEdx(1, >0.4GeV; -10>3.5GeV)
    cuts.AddCut("84610113", "00200009326000008250400000", "0162103500000000"); // PidEdx(1, >0.25GeV)
    cuts.AddCut("84610113", "00200009326200008250400000", "0162103500000000"); // PidEdx(1, >0.25GeV; -10>4GeV)
    cuts.AddCut("84610113", "00200009327200008250400000", "0162103500000000"); // PidEdx(1, >0.4GeV; -10>4GeV)
  } else if (trainConfig == 436) { // qt & psipair variation
    cuts.AddCut("84610113", "00200009327000003250400000", "0162103500000000"); // qT 0.05 1D
    cuts.AddCut("84610113", "00200009327000009250400000", "0162103500000000"); // qT 0.03 2D
    cuts.AddCut("84610113", "00200009327000002250400000", "0162103500000000"); // qT 0.07 1D
    cuts.AddCut("84610113", "00200009327000008240400000", "0162103500000000"); // Psi Pair: 1D 0.2
  } else if (trainConfig == 437) { // PsiPair Variation
    cuts.AddCut("84610113", "00200009327000008210400000", "0162103500000000"); // Psi Pair: 1D 0.1
    cuts.AddCut("84610113", "00200009327000008150400000", "0162103500000000"); // chi2 50  2D
    cuts.AddCut("84610113", "00200009327000008850400000", "0162103500000000"); // chi2 20  2D
    cuts.AddCut("84610113", "00200009327000008260400000", "0162103500000000"); // psi pair 0.05
  } else if (trainConfig == 438) { // cos point & meson alpha variation
    cuts.AddCut("84610113", "00200009327000008250300000", "0162103500000000"); // cos pointing angle 0.75
    cuts.AddCut("84610113", "00200009327000008250600000", "0162103500000000"); // cos pointing angle 0.9
    cuts.AddCut("84610113", "00200009327000008250400000", "0162106500000000"); // alpha meson cut 0.8
  } else if (trainConfig == 439) { // BG mixing & smear variations
    cuts.AddCut("84610113", "00200009327000008250400000", "0262103500000000"); // BG track multiplicity
    cuts.AddCut("84610113", "00200009327000008250400000", "0162103500800000"); // fPSigSmearingCte=0.014;
    cuts.AddCut("84610113", "00200009327000008250400000", "0162103500900000"); // fPSigSmearingCte=0.014;

  //--------------------------------------------------------------------------
  // Systematics variations for direct photons pPb 5 TeV Run1 and Run2 - 60-100%
  //--------------------------------------------------------------------------
  } else if (trainConfig == 440) { //default + dEdx Variation
    cuts.AddCut("86010113", "00200009327000008250400000", "0162103500000000"); // new default
    cuts.AddCut("86010113", "00200009217000008260400000", "0162103500000000"); // old standard cut
    cuts.AddCut("86010113", "00200009327000008250404000", "0162103500000000"); // new default w/ to clos v0
  } else if (trainConfig == 441) { // gamma eta & meson rap var
    cuts.AddCut("86010113", "03200009327000008250400000", "0162303500000000"); // |eta| < 0.65, |y| < 0.6
    cuts.AddCut("86010113", "04200009327000008250400000", "0162203500000000"); // |eta| < 0.75, |y| < 0.7
    cuts.AddCut("86010113", "01200009327000008250400000", "0162403500000000"); // |eta| < 0.6, |y| < 0.5
  } else if (trainConfig == 442) { // minR and single pt var
    cuts.AddCut("86010113", "00100009327000008250400000", "0162103500000000"); // minR 2.8
    cuts.AddCut("86010113", "00900009327000008250400000", "0162103500000000"); // minR 7.5
    cuts.AddCut("86010113", "00200079327000008250400000", "0162103500000000"); // single pT 0. GeV/c
    cuts.AddCut("86010113", "00200019327000008250400000", "0162103500000000"); // single pT 0.1 GeV/c
  } else if (trainConfig == 443) { // TPC cluster & edEdx var
    cuts.AddCut("86010113", "00200008327000008250400000", "0162103500000000"); // TPC Cluster 0.35
    cuts.AddCut("86010113", "00200006327000008250400000", "0162103500000000"); // TPC Cluster 0.7
    cuts.AddCut("86010113", "00200009227000008250400000", "0162103500000000"); // edEdx -4,5
    cuts.AddCut("86010113", "00200009627000008250400000", "0162103500000000"); // edEdx -2.5,4
    cuts.AddCut("86010113", "00200009127000008250400000", "0162103500000000"); // edEdx 5,5
  } else if (trainConfig == 444) { //PidEdx Variation
    cuts.AddCut("86010113", "00200009357000008250400000", "0162103500000000"); // PidEdx(2, >0.4GeV)
    cuts.AddCut("86010113", "00200009387300008250400000", "0162103500000000"); // PidEdx(2, >0.4GeV; 1>3.5GeV)
    cuts.AddCut("86010113", "00200009320000008250400000", "0162103500000000"); // PidEdx(1, >0.5GeV)
    cuts.AddCut("86010113", "00200009325000008250400000", "0162103500000000"); // PidEdx(1, >0.3GeV)
  } else if (trainConfig == 445) { //PidEdx Variation 2
    cuts.AddCut("86010113", "00200009327300008250400000", "0162103500000000"); // PidEdx(1, >0.4GeV; -10>3.5GeV)
    cuts.AddCut("86010113", "00200009326000008250400000", "0162103500000000"); // PidEdx(1, >0.25GeV)
    cuts.AddCut("86010113", "00200009326200008250400000", "0162103500000000"); // PidEdx(1, >0.25GeV; -10>4GeV)
    cuts.AddCut("86010113", "00200009327200008250400000", "0162103500000000"); // PidEdx(1, >0.4GeV; -10>4GeV)
  } else if (trainConfig == 446) { // qt & psipair variation
    cuts.AddCut("86010113", "00200009327000003250400000", "0162103500000000"); // qT 0.05 1D
    cuts.AddCut("86010113", "00200009327000009250400000", "0162103500000000"); // qT 0.03 2D
    cuts.AddCut("86010113", "00200009327000002250400000", "0162103500000000"); // qT 0.07 1D
    cuts.AddCut("86010113", "00200009327000008240400000", "0162103500000000"); // Psi Pair: 1D 0.2
  } else if (trainConfig == 447) { // PsiPair Variation
    cuts.AddCut("86010113", "00200009327000008210400000", "0162103500000000"); // Psi Pair: 1D 0.1
    cuts.AddCut("86010113", "00200009327000008150400000", "0162103500000000"); // chi2 50  2D
    cuts.AddCut("86010113", "00200009327000008850400000", "0162103500000000"); // chi2 20  2D
    cuts.AddCut("86010113", "00200009327000008260400000", "0162103500000000"); // psi pair 0.05
  } else if (trainConfig == 448) { // cos point & meson alpha variation
    cuts.AddCut("86010113", "00200009327000008250300000", "0162103500000000"); // cos pointing angle 0.75
    cuts.AddCut("86010113", "00200009327000008250600000", "0162103500000000"); // cos pointing angle 0.9
    cuts.AddCut("86010113", "00200009327000008250400000", "0162106500000000"); // alpha meson cut 0.8
  } else if (trainConfig == 449) { // BG mixing & smear variations
    cuts.AddCut("86010113", "00200009327000008250400000", "0262103500000000"); // BG track multiplicity
    cuts.AddCut("86010113", "00200009327000008250400000", "0162103500800000"); // fPSigSmearingCte=0.014;
    cuts.AddCut("86010113", "00200009327000008250400000", "0162103500900000"); // fPSigSmearingCte=0.014;


  //Systematics HM
  } else if (trainConfig == 800) {//MB
    cuts.AddCut("80010113", "00200009a27000008250a04120", "0162103500000000"); // default
    cuts.AddCut("80010113", "00200079a27000008250a04120", "0162103500000000"); // min pT no cut
    cuts.AddCut("80010113", "00200069a27000008250a04120", "0162103500000000"); // min pT 40 cut
    cuts.AddCut("80010113", "00200049a27000008250a04120", "0162103500000000"); // min pT 75 cut
    cuts.AddCut("80010113", "00200008a27000008250a04120", "0162103500000000"); // TPC cluster 35%
    cuts.AddCut("80010113", "00200007a27000008250a04120", "0162103500000000"); // TPC cluster 70%
    cuts.AddCut("80010113", "00200009b27000008250a04120", "0162103500000000"); // edEdx -3.2,3.2
    cuts.AddCut("80010113", "00200009c27000008250a04120", "0162103500000000"); // edEdx -2.8,2.8
 } else if (trainConfig == 801) {//0-10
    cuts.AddCut("80110113", "00200009a27000008250a04120", "0162103500000000"); // default
    cuts.AddCut("80110113", "00200079a27000008250a04120", "0162103500000000"); // min pT no cut
    cuts.AddCut("80110113", "00200069a27000008250a04120", "0162103500000000"); // min pT 40 cut
    cuts.AddCut("80110113", "00200049a27000008250a04120", "0162103500000000"); // min pT 75 cut
    cuts.AddCut("80110113", "00200008a27000008250a04120", "0162103500000000"); // TPC cluster 35%
    cuts.AddCut("80110113", "00200006a27000008250a04120", "0162103500000000"); // TPC cluster 70%
    cuts.AddCut("80110113", "00200009b27000008250a04120", "0162103500000000"); // edEdx -3,2,3.2
    cuts.AddCut("80110113", "00200009c27000008250a04120", "0162103500000000"); // edEdx -2.8,2.8
  } else if (trainConfig == 802) {//0-20
    cuts.AddCut("80210113", "00200009a27000008250a04120", "0162103500000000"); // default
    cuts.AddCut("80210113", "00200079a27000008250a04120", "0162103500000000"); // min pT no cut
    cuts.AddCut("80210113", "00200069a27000008250a04120", "0162103500000000"); // min pT 40 cut
    cuts.AddCut("80210113", "00200049a27000008250a04120", "0162103500000000"); // min pT 75 cut
    cuts.AddCut("80210113", "00200008a27000008250a04120", "0162103500000000"); // TPC cluster 35%
    cuts.AddCut("80210113", "00200006a27000008250a04120", "0162103500000000"); // TPC cluster 70%
    cuts.AddCut("80210113", "00200009b27000008250a04120", "0162103500000000"); // edEdx -3,2,3.2
    cuts.AddCut("80210113", "00200009c27000008250a04120", "0162103500000000"); // edEdx -2.8,2.8
  } else if (trainConfig == 803) {//20-40
    cuts.AddCut("82410113", "00200009a27000008250a04120", "0162103500000000"); // default
    cuts.AddCut("82410113", "00200079a27000008250a04120", "0162103500000000"); // min pT no cut
    cuts.AddCut("82410113", "00200069a27000008250a04120", "0162103500000000"); // min pT 40 cut
    cuts.AddCut("82410113", "00200049a27000008250a04120", "0162103500000000"); // min pT 75 cut
    cuts.AddCut("82410113", "00200008a27000008250a04120", "0162103500000000"); // TPC cluster 35%
    cuts.AddCut("82410113", "00200006a27000008250a04120", "0162103500000000"); // TPC cluster 70%
    cuts.AddCut("82410113", "00200009b27000008250a04120", "0162103500000000"); // edEdx -3.2,3.2
    cuts.AddCut("82410113", "00200009c27000008250a04120", "0162103500000000"); // edEdx -2.8,2.8
  } else if (trainConfig == 804) {//40-60
    cuts.AddCut("84610113", "00200009a27000008250a04120", "0162103500000000"); // default
    cuts.AddCut("84610113", "00200079a27000008250a04120", "0162103500000000"); // min pT no cut
    cuts.AddCut("84610113", "00200069a27000008250a04120", "0162103500000000"); // min pT 40 cut
    cuts.AddCut("84610113", "00200049a27000008250a04120", "0162103500000000"); // min pT 75 cut
    cuts.AddCut("84610113", "00200008a27000008250a04120", "0162103500000000"); // TPC cluster 35%
    cuts.AddCut("84610113", "00200006a27000008250a04120", "0162103500000000"); // TPC cluster 70%
    cuts.AddCut("84610113", "00200009b27000008250a04120", "0162103500000000"); // edEdx -3.2,3.2
    cuts.AddCut("84610113", "00200009c27000008250a04120", "0162103500000000"); // edEdx -2.8,2.8
  } else if (trainConfig == 805) {//60-80
    cuts.AddCut("86810113", "00200009a27000008250a04120", "0162103500000000"); // default
    cuts.AddCut("86810113", "00200079a27000008250a04120", "0162103500000000"); // min pT no cut
    cuts.AddCut("86810113", "00200069a27000008250a04120", "0162103500000000"); // min pT 40 cut
    cuts.AddCut("86810113", "00200049a27000008250a04120", "0162103500000000"); // min pT 75 cut
    cuts.AddCut("86810113", "00200008a27000008250a04120", "0162103500000000"); // TPC cluster 35%
    cuts.AddCut("86810113", "00200006a27000008250a04120", "0162103500000000"); // TPC cluster 70%
    cuts.AddCut("86810113", "00200009b27000008250a04120", "0162103500000000"); // edEdx -3.2,3.2
    cuts.AddCut("86810113", "00200009c27000008250a04120", "0162103500000000"); // edEdx -2.8,2.8
  } else if (trainConfig == 806) {//80-100
    cuts.AddCut("88010113", "00200009a27000008250a04120", "0162103500000000"); // default
    cuts.AddCut("88010113", "00200079a27000008250a04120", "0162103500000000"); // min pT no cut
    cuts.AddCut("88010113", "00200069a27000008250a04120", "0162103500000000"); // min pT 40 cut
    cuts.AddCut("88010113", "00200049a27000008250a04120", "0162103500000000"); // min pT 75 cut
    cuts.AddCut("88010113", "00200008a27000008250a04120", "0162103500000000"); // TPC cluster 35%
    cuts.AddCut("88010113", "00200006a27000008250a04120", "0162103500000000"); // TPC cluster 70%
    cuts.AddCut("88010113", "00200009b27000008250a04120", "0162103500000000"); // edEdx -3.2,3.2
    cuts.AddCut("88010113", "00200009c27000008250a04120", "0162103500000000"); // edEdx -2.8,2.8
  } else if (trainConfig == 807) {//MB
    cuts.AddCut("80010113", "00200009a57000008250a04120", "0162103500000000"); // pidEdx 2,-10
    cuts.AddCut("80010113", "00200009a17000008250a04120", "0162103500000000"); // pidEdx 0,-10
    cuts.AddCut("80010113", "00200009a27000003250a04120", "0162103500000000"); // qT max 0.05 1D
    cuts.AddCut("80010113", "00200009a27000002250a04120", "0162103500000000"); // qT max 0.06 2D
    cuts.AddCut("80010113", "00200009a27000008250904120", "0162103500000000"); // CosPA 0.99
    cuts.AddCut("80010113", "00200009a27000008250b04120", "0162103500000000"); // CosPA 0.985
  } else if (trainConfig == 808) {//0-10
    cuts.AddCut("80110113", "00200009a57000008250a04120", "0162103500000000"); // pidEdx 2,-10
    cuts.AddCut("80110113", "00200009a17000008250a04120", "0162103500000000"); // pidEdx 0,-10
    cuts.AddCut("80110113", "00200009a27000003250a04120", "0162103500000000"); // qT max 0.05 1D
    cuts.AddCut("80110113", "00200009a27000002250a04120", "0162103500000000"); // qT max 0.06 2D
    cuts.AddCut("80110113", "00200009a27000008250904120", "0162103500000000"); // CosPA 0.99
    cuts.AddCut("80110113", "00200009a27000008250b04120", "0162103500000000"); // CosPA 0.985
  } else if (trainConfig == 809) {//0-20
    cuts.AddCut("80210113", "00200009a57000008250a04120", "0162103500000000"); // pidEdx 2,-10
    cuts.AddCut("80210113", "00200009a17000008250a04120", "0162103500000000"); // pidEdx 0,-10
    cuts.AddCut("80210113", "00200009a27000003250a04120", "0162103500000000"); // qT max 0.05 1D
    cuts.AddCut("80210113", "00200009a27000002250a04120", "0162103500000000"); // qT max 0.06 2D
    cuts.AddCut("80210113", "00200009a27000008250904120", "0162103500000000"); // CosPA 0.99
    cuts.AddCut("80210113", "00200009a27000008250b04120", "0162103500000000"); // CosPA 0.985
  } else if (trainConfig == 810) {//20-40
    cuts.AddCut("82410113", "00200009a57000008250a04120", "0162103500000000"); // pidEdx 2,-10
    cuts.AddCut("82410113", "00200009a17000008250a04120", "0162103500000000"); // pidEdx 0,-10
    cuts.AddCut("82410113", "00200009a27000003250a04120", "0162103500000000"); // qT max 0.05 1D
    cuts.AddCut("82410113", "00200009a27000002250a04120", "0162103500000000"); // qT max 0.06 2D
    cuts.AddCut("82410113", "00200009a27000008250904120", "0162103500000000"); // CosPA 0.99
    cuts.AddCut("82410113", "00200009a27000008250b04120", "0162103500000000"); // CosPA 0.985
  } else if (trainConfig == 811) {//40-60
    cuts.AddCut("84610113", "00200009a57000008250a04120", "0162103500000000"); // pidEdx 2,-10
    cuts.AddCut("84610113", "00200009a17000008250a04120", "0162103500000000"); // pidEdx 0,-10
    cuts.AddCut("84610113", "00200009a27000003250a04120", "0162103500000000"); // qT max 0.05 1D
    cuts.AddCut("84610113", "00200009a27000002250a04120", "0162103500000000"); // qT max 0.06 2D
    cuts.AddCut("84610113", "00200009a27000008250904120", "0162103500000000"); // CosPA 0.99
    cuts.AddCut("84610113", "00200009a27000008250b04120", "0162103500000000"); // CosPA 0.985
  } else if (trainConfig == 812) {//60-80
    cuts.AddCut("86810113", "00200009a27000008250a04120", "0162103500000000"); // pidEdx 2,-10
    cuts.AddCut("86810113", "00200009a27000008250a04120", "0162103500000000"); // pidEdx 0,-10
    cuts.AddCut("86810113", "00200009a27000003250a04120", "0162103500000000"); // qT max 0.05 1D
    cuts.AddCut("86810113", "00200009a27000002250a04120", "0162103500000000"); // qT max 0.06 2D
    cuts.AddCut("86810113", "00200009a27000008250904120", "0162103500000000"); // CosPA 0.99
    cuts.AddCut("86810113", "00200009a27000008250b04120", "0162103500000000"); // CosPA 0.985
  } else if (trainConfig == 813) {//80-100
    cuts.AddCut("88010113", "00200009a27000008250a04120", "0162103500000000"); // pidEdx 2,-10
    cuts.AddCut("88010113", "00200009a27000008250a04120", "0162103500000000"); // pidEdx 0,-10
    cuts.AddCut("88010113", "00200009a27000003250a04120", "0162103500000000"); // qT max 0.05 1D
    cuts.AddCut("88010113", "00200009a27000002250a04120", "0162103500000000"); // qT max 0.06 2D
    cuts.AddCut("88010113", "00200009a27000008250904120", "0162103500000000"); // CosPA 0.99
    cuts.AddCut("88010113", "00200009a27000008250b04120", "0162103500000000"); // CosPA 0.985
  } else if (trainConfig == 814) {//MB
    cuts.AddCut("80010113", "00200009a27000008250a04120", "0162103500000000"); // default
    cuts.AddCut("80010113", "00200009a27000008250c04120", "0162103500000000"); // cosPA 0.996
    cuts.AddCut("80010113", "00200009a27000008250d04120", "0162103500000000"); // cosPA 0.997
    cuts.AddCut("80010113", "00200009a27000008250e04120", "0162103500000000"); // cosPA 0.998
    cuts.AddCut("80010113", "00200009a27000008250f04120", "0162103500000000"); // cosPA 0.999
  } else if (trainConfig == 815) {//0-10
    cuts.AddCut("80110113", "00200009a27000008250a04120", "0162103500000000"); // default
    cuts.AddCut("80110113", "00200009a27000008250c04120", "0162103500000000"); // cosPA 0.996
    cuts.AddCut("80110113", "00200009a27000008250d04120", "0162103500000000"); // cosPA 0.997
    cuts.AddCut("80110113", "00200009a27000008250e04120", "0162103500000000"); // cosPA 0.998
    cuts.AddCut("80110113", "00200009a27000008250f04120", "0162103500000000"); // cosPA 0.999
  } else if (trainConfig == 816) {//0-20
    cuts.AddCut("80210113", "00200009a27000008250a04120", "0162103500000000"); // default
    cuts.AddCut("80210113", "00200009a27000008250c04120", "0162103500000000"); // cosPA 0.996
    cuts.AddCut("80210113", "00200009a27000008250d04120", "0162103500000000"); // cosPA 0.997
    cuts.AddCut("80210113", "00200009a27000008250e04120", "0162103500000000"); // cosPA 0.998
    cuts.AddCut("80210113", "00200009a27000008250f04120", "0162103500000000"); // cosPA 0.999
  } else if (trainConfig == 817) {//20-40
    cuts.AddCut("82410113", "00200009a27000008250a04120", "0162103500000000"); // default
    cuts.AddCut("82410113", "00200009a27000008250c04120", "0162103500000000"); // cosPA 0.996
    cuts.AddCut("82410113", "00200009a27000008250d04120", "0162103500000000"); // cosPA 0.997
    cuts.AddCut("82410113", "00200009a27000008250e04120", "0162103500000000"); // cosPA 0.998
    cuts.AddCut("82410113", "00200009a27000008250f04120", "0162103500000000"); // cosPA 0.999
  } else if (trainConfig == 818) {//40-60
    cuts.AddCut("84610113", "00200009a27000008250a04120", "0162103500000000"); // default
    cuts.AddCut("84610113", "00200009a27000008250c04120", "0162103500000000"); // cosPA 0.996
    cuts.AddCut("84610113", "00200009a27000008250d04120", "0162103500000000"); // cosPA 0.997
    cuts.AddCut("84610113", "00200009a27000008250e04120", "0162103500000000"); // cosPA 0.998
    cuts.AddCut("84610113", "00200009a27000008250f04120", "0162103500000000"); // cosPA 0.999
  } else if (trainConfig == 819) {//60-80
    cuts.AddCut("86810113", "00200009a27000008250a04120", "0162103500000000"); // default
    cuts.AddCut("86810113", "00200009a27000008250c04120", "0162103500000000"); // cosPA 0.996
    cuts.AddCut("86810113", "00200009a27000008250d04120", "0162103500000000"); // cosPA 0.997
    cuts.AddCut("86810113", "00200009a27000008250e04120", "0162103500000000"); // cosPA 0.998
    cuts.AddCut("86810113", "00200009a27000008250f04120", "0162103500000000"); // cosPA 0.999
  } else if (trainConfig == 820) {//80-100
    cuts.AddCut("88010113", "00200009a27000008250a04120", "0162103500000000"); // default
    cuts.AddCut("88010113", "00200009a27000008250c04120", "0162103500000000"); // cosPA 0.996
    cuts.AddCut("88010113", "00200009a27000008250d04120", "0162103500000000"); // cosPA 0.997
    cuts.AddCut("88010113", "00200009a27000008250e04120", "0162103500000000"); // cosPA 0.998
    cuts.AddCut("88010113", "00200009a27000008250f04120", "0162103500000000"); // cosPA 0.999







  } else {
    Error(Form("GammaConvV1_%i",trainConfig), "wrong trainConfig variable no cuts have been specified for the configuration");
    return;
  }

  if(!cuts.AreValid()){
    cout << "\n\n****************************************************" << endl;
    cout << "ERROR: No valid cuts stored in CutHandlerConvpPb! Returning..." << endl;
    cout << "****************************************************\n\n" << endl;
    return;
  }

  Int_t numberOfCuts = cuts.GetNCuts();

  TList *EventCutList = new TList();
  TList *ConvCutList = new TList();
  TList *MesonCutList = new TList();
  TList *ClusterCutList = new TList();

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
  if (periodNameV0Reader.Contains("LHC18b9")||periodNameV0Reader.Contains("LHC17g8")){
    TObjString *HeaderPMB = new TObjString("EPOSLHC_0");
    TObjString *HeaderP8J = new TObjString("Pythia8Jets_1");
    if (doWeightingPart==4) {
      HeaderList->Add(HeaderPMB);
      HeaderList->Add(HeaderP8J);
    } else {
      HeaderList->Add(HeaderP8J);
    }
  }

  Bool_t doWeighting = kFALSE;
  if (doWeightingPart == 1 || doWeightingPart == 2 || doWeightingPart == 3) doWeighting = kTRUE;

  EventCutList->SetOwner(kTRUE);
  AliConvEventCuts **analysisEventCuts = new AliConvEventCuts*[numberOfCuts];
  ConvCutList->SetOwner(kTRUE);
  AliConversionPhotonCuts **analysisCuts = new AliConversionPhotonCuts*[numberOfCuts];
  MesonCutList->SetOwner(kTRUE);
  AliConversionMesonCuts **analysisMesonCuts = new AliConversionMesonCuts*[numberOfCuts];
  ClusterCutList->SetOwner(kTRUE);
  AliCaloPhotonCuts **analysisClusterCuts     = new AliCaloPhotonCuts*[numberOfCuts];
  Bool_t enableClustersForTrigger             = kFALSE;
  Bool_t initializedMatBudWeigths_existing    = kFALSE;

  if (doWeighting) Printf("weighting has been switched on");

  for(Int_t i = 0; i<numberOfCuts; i++){
    cout << "initialization of cutnumber: " << i << endl;

    analysisEventCuts[i] = new AliConvEventCuts();
    if ( trainConfig == 1 || trainConfig == 2 ||  trainConfig == 5 || trainConfig == 6 ||
        ( trainConfig > 99 && trainConfig < 120 ) ||
        ( trainConfig > 199 && trainConfig < 210 ) ){
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
    if ( trainConfig == 3 || trainConfig == 4 || trainConfig == 7 || trainConfig == 8 ||
        ( trainConfig > 119 && trainConfig < 140 ) ||
        ( trainConfig > 209 && trainConfig < 220 ) ){
      if (doWeighting){
        analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13e7_addSig_pPb_5023GeV_MBV0A",
                                      "Eta_Hijing_LHC13e7_addSig_pPb_5023GeV_MBV0A", "","Pi0_Fit_Data_pPb_5023GeV_MBV0A","Eta_Fit_Data_pPb_5023GeV_MBV0A");
      }
    }

    TString dataInputMultHisto  = "";
    TString mcInputMultHisto    = "";
    TString triggerString       = (cuts.GetEventCut(i)).Data();
    triggerString               = triggerString(3,2);

    dataInputMultHisto          = Form("%s_%s", periodNameAnchor.Data(), triggerString.Data());
    mcInputMultHisto            = Form("%s_%s", periodNameV0Reader.Data(), triggerString.Data());

    if (doMultiplicityWeighting){
      cout << "enabling mult weighting" << endl;
      analysisEventCuts[i]->SetUseWeightMultiplicityFromFile( kTRUE, fileNameInputForMultWeighing, dataInputMultHisto, mcInputMultHisto );
    }

    if (debugLevel > 0) analysisEventCuts[i]->SetDebugLevel(debugLevel);
    analysisEventCuts[i]->SetTriggerMimicking(enableTriggerMimicking);
    analysisEventCuts[i]->SetTriggerOverlapRejecion(enableTriggerOverlapRej);
    analysisEventCuts[i]->SetMaxFacPtHard(maxFacPtHard);
    analysisEventCuts[i]->SetCorrectionTaskSetting(corrTaskSetting);
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
    cout << "initialized event cut: " << (cuts.GetEventCut(i)).Data() << endl;

    if ( trainConfig == 302 || trainConfig == 303  || trainConfig == 312 || trainConfig == 313 ){
        TString caloCutPos = cuts.GetClusterCut(i);
        caloCutPos.Resize(1);
        TString TrackMatcherName = Form("CaloTrackMatcher_%s_%i",caloCutPos.Data(),trackMatcherRunningMode);
        if( !(AliCaloTrackMatcher*)mgr->GetTask(TrackMatcherName.Data()) ){
          AliCaloTrackMatcher* fTrackMatcher = new AliCaloTrackMatcher(TrackMatcherName.Data(),caloCutPos.Atoi(),trackMatcherRunningMode);
          fTrackMatcher->SetV0ReaderName(V0ReaderName);
          fTrackMatcher->SetCorrectionTaskSetting(corrTaskSetting);
          mgr->AddTask(fTrackMatcher);
          mgr->ConnectInput(fTrackMatcher,0,cinput);
        }

        enableClustersForTrigger  = kTRUE;
        analysisClusterCuts[i]    = new AliCaloPhotonCuts();
        analysisClusterCuts[i]->SetV0ReaderName(V0ReaderName);
        analysisClusterCuts[i]->SetCorrectionTaskSetting(corrTaskSetting);
        analysisClusterCuts[i]->SetLightOutput(runLightOutput);
        analysisClusterCuts[i]->SetCaloTrackMatcherName(TrackMatcherName);
        analysisClusterCuts[i]->InitializeCutsFromCutString((cuts.GetClusterCut(i)).Data());
        ClusterCutList->Add(analysisClusterCuts[i]);
        analysisClusterCuts[i]->SetFillCutHistograms("");
    }


    analysisCuts[i] = new AliConversionPhotonCuts();
    if ( trainConfig > 149 && trainConfig < 156 ){
      analysisCuts[i]->SetSwitchToKappaInsteadOfNSigdEdxTPC(kTRUE);
    }

    if (enableMatBudWeightsPi0 > 0){
      if (isMC > 0){
        if (analysisCuts[i]->InitializeMaterialBudgetWeights(enableMatBudWeightsPi0,filenameMatBudWeights)){
          initializedMatBudWeigths_existing = kTRUE;}
        else {cout << "ERROR The initialization of the materialBudgetWeights did not work out." << endl;}
      } else {cout << "ERROR 'enableMatBudWeightsPi0'-flag was set > 0 even though this is not a MC task. It was automatically reset to 0." << endl;}
    }
    if (enableElecDeDxPostCalibration>0){
      if (isMC == 0){
	if( analysisCuts[i]->InitializeElecDeDxPostCalibration(fileNameElecDeDxPostCalibration)){
	  analysisCuts[i]->SetDoElecDeDxPostCalibration(enableElecDeDxPostCalibration);
	} else {
	  enableElecDeDxPostCalibration=kFALSE;
	  analysisCuts[i]->SetDoElecDeDxPostCalibration(enableElecDeDxPostCalibration);
	}
      } else{
	cout << "ERROR enableElecDeDxPostCalibration set to True even if MC file. Automatically reset to 0"<< endl;
	enableElecDeDxPostCalibration=kFALSE;
	analysisCuts[i]->SetDoElecDeDxPostCalibration(enableElecDeDxPostCalibration);
      }
    }
    analysisCuts[i]->SetV0ReaderName(V0ReaderName);
    analysisCuts[i]->SetLightOutput(runLightOutput);
    analysisCuts[i]->InitializeCutsFromCutString((cuts.GetPhotonCut(i)).Data());

    analysisCuts[i]->SetIsHeavyIon(isHeavyIon);
    ConvCutList->Add(analysisCuts[i]);
    analysisCuts[i]->SetFillCutHistograms("",kFALSE);
    cout << "initialized photon cut: " << (cuts.GetPhotonCut(i)).Data() << endl;

    analysisMesonCuts[i] = new AliConversionMesonCuts();
    analysisMesonCuts[i]->SetLightOutput(runLightOutput);
    analysisMesonCuts[i]->InitializeCutsFromCutString((cuts.GetMesonCut(i)).Data());
    MesonCutList->Add(analysisMesonCuts[i]);
    analysisMesonCuts[i]->SetFillCutHistograms("");
    analysisEventCuts[i]->SetAcceptedHeader(HeaderList);
    cout << "initialized meson cut: " << (cuts.GetMesonCut(i)).Data() << endl;
  }

  task->SetDoTHnSparse(runTHnSparse);
  task->SetEventCutList(numberOfCuts,EventCutList);
  task->SetConversionCutList(numberOfCuts,ConvCutList);
  task->SetMesonCutList(numberOfCuts,MesonCutList);
  task->SetMoveParticleAccordingToVertex(kTRUE);
  task->SetDoMesonAnalysis(kTRUE);
  task->SetDoMesonQA(enableQAMesonTask); //Attention new switch for Pi0 QA
  task->SetDoPhotonQA(enableQAPhotonTask); //Attention new switch small for Photon QA
  task->SetDoPlotVsCentrality(enablePlotVsCentrality);
  if (enableClustersForTrigger){
    task->SetDoClusterSelectionForTriggerNorm(enableClustersForTrigger);
    task->SetClusterCutList(numberOfCuts,ClusterCutList);
  }

  if (initializedMatBudWeigths_existing) {
      task->SetDoMaterialBudgetWeightingOfGammasForTrueMesons(kTRUE);
  }


  //connect containers
  AliAnalysisDataContainer *coutput =
      mgr->CreateContainer(!(corrTaskSetting.CompareTo("")) ? Form("GammaConvV1_%i",trainConfig) : Form("GammaConvV1_%i_%s",trainConfig,corrTaskSetting.Data()), TList::Class(),
                AliAnalysisManager::kOutputContainer, Form("GammaConvV1_%i.root",trainConfig) );

  mgr->AddTask(task);
  mgr->ConnectInput(task,0,cinput);
  mgr->ConnectOutput(task,1,coutput);
  Int_t nContainer = 2;
  for(Int_t i = 0; i<numberOfCuts; i++){
    if(enableQAPhotonTask>1){
      mgr->ConnectOutput(task,nContainer,mgr->CreateContainer(Form("%s_%s_%s Photon DCA tree",(cuts.GetEventCut(i)).Data(),(cuts.GetPhotonCut(i)).Data(),(cuts.GetMesonCut(i)).Data()), TTree::Class(), AliAnalysisManager::kOutputContainer, Form("GammaConvV1_%i.root",trainConfig)) );
      nContainer++;
    }
    if(enableQAMesonTask>1){
      mgr->ConnectOutput(task,nContainer,mgr->CreateContainer(Form("%s_%s_%s Meson DCA tree",(cuts.GetEventCut(i)).Data(),(cuts.GetPhotonCut(i)).Data(),(cuts.GetMesonCut(i)).Data()), TTree::Class(), AliAnalysisManager::kOutputContainer, Form("GammaConvV1_%i.root",trainConfig)) );
      nContainer++;
    }
  }

  return;

}
