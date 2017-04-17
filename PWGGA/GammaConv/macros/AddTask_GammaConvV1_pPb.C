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
                                TString   additionalTrainConfig         = "0"                             // additional counter for trainconfig, this has to be always the last parameter
                          ) {


  
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
      }
    }
  }
  
  
  TString sAdditionalTrainConfig = rAdditionalTrainConfig->GetString();
  if (sAdditionalTrainConfig.Atoi() > 0){
    trainConfig = trainConfig + sAdditionalTrainConfig.Atoi();
    cout << "INFO: AddTask_GammaConvV1_pPb running additionalTrainConfig '" << sAdditionalTrainConfig.Atoi() << "', train config: '" << trainConfig << "'" << endl;
  }

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

  Bool_t isMCForOtherSettings = 0;
  if (isMC > 0) isMCForOtherSettings = 1;
  //========= Add PID Reponse to ANALYSIS manager ====
  if(!(AliPIDResponse*)mgr->GetTask("PIDResponseTask")){
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
    AddTaskPIDResponse(isMCForOtherSettings);
  }

  //=========  Set Cutnumber for V0Reader ================================
  TString cutnumberPhoton = "06000008000100001500000000";
  TString cutnumberEvent = "80000103";
  if (trainConfig==189 || trainConfig==190 || trainConfig==191 || trainConfig==192) cutnumberEvent = "80000003";

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
    if(inputHandler->IsA()==AliAODInputHandler::Class()){
    // AOD mode
      cout << "AOD handler: adding " << cutnumberAODBranch.Data() << " as conversion branch" << endl;
      fV0ReaderV1->SetDeltaAODBranchName(Form("GammaConv_%s_gamma",cutnumberAODBranch.Data()));
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
  
  if(trainConfig == 1){
    cuts.AddCut("80000113", "03200009217000008260400000", "0162303500900000"); //New STANDARD CUT |eta| < 0.65, |y| < 0.6
    cuts.AddCut("80000113", "04200009217000008260400000", "0162203500900000"); //New STANDARD CUT |eta| < 0.75, |y| < 0.7
    cuts.AddCut("80000113", "01200009217000008260400000", "0162403500900000"); //New STANDARD CUT |eta| < 0.6, |y| < 0.5
    cuts.AddCut("80000113", "00200009237002003220000000", "0162103500900000");
  } else if (trainConfig == 2) {
    cuts.AddCut("80000123", "03200009217000008260400000", "0162303500900000"); //New STANDARD CUT |eta| < 0.65, |y| < 0.6
    cuts.AddCut("80000123", "04200009217000008260400000", "0162203500900000"); //New STANDARD CUT |eta| < 0.75, |y| < 0.7
    cuts.AddCut("80000123", "01200009217000008260400000", "0162403500900000"); //New STANDARD CUT |eta| < 0.6, |y| < 0.5
    cuts.AddCut("80000123", "00200009237002003220000000", "0162103500900000");
  } else if (trainConfig == 3) {
    cuts.AddCut("80000113", "00200009277002003220000000", "0162103500900000");
    cuts.AddCut("80000113", "00200009255102003220000000", "0162103500900000");
    cuts.AddCut("80000113", "00200009217000003220000000", "0162103500900000"); //just tighten Psi pair
    cuts.AddCut("80000113", "00200009217000003260000000", "0162103500900000"); //tighten Psi pair and chi2 in 2D
  } else if (trainConfig == 4) {
    cuts.AddCut("80000123", "00200009277002003220000000", "0162103500900000");
    cuts.AddCut("80000123", "00200009255102003220000000", "0162103500900000");
    cuts.AddCut("80000123", "00200009217000003220000000", "0162103500900000"); //just tighten Psi pair
    cuts.AddCut("80000123", "00200009217000003260000000", "0162103500900000"); //tighten Psi pair and chi2 in 2D
  } else if (trainConfig == 5) {
    cuts.AddCut("80000113", "00200009217000008220000000", "0162103500900000"); //tighten psi pair and qt in 2D
    cuts.AddCut("80000113", "00200009217000008260000000", "0162103500900000"); //tighten psi pair and chi2 in 2D and qt in 2D
    cuts.AddCut("80000113", "00200009217000008220400000", "0162103500900000"); //clean cuts
    cuts.AddCut("80000113", "00200009217000008260400000", "0162103500900000"); //clean cuts ///New STANDARD CUT
  } else if (trainConfig == 6) {
    cuts.AddCut("80000123", "00200009217000008220000000", "0162103500900000"); //tighten psi pair and qt in 2D
    cuts.AddCut("80000123", "00200009217000008260000000", "0162103500900000"); //tighten psi pair and chi2 in 2D and qt in 2D
    cuts.AddCut("80000123", "00200009217000008220400000", "0162103500900000"); //clean cuts
    cuts.AddCut("80000123", "00200009217000008260400000", "0162103500900000"); //clean cuts ///New STANDARD CUT
  } else if (trainConfig == 7) {
    cuts.AddCut("80000113", "00100009217000008260400000", "0162103500900000");  //new standard eta=0.9 y=0.8  // minR 2.8
    cuts.AddCut("80000113", "00900009217000008260400000", "0162103500900000");  //new standard eta=0.9 y=0.8  //minR 7.5
    cuts.AddCut("80000113", "00200079217000008260400000", "0162103500900000");  //new standard eta=0.9 y=0.8  //single pT 0. GeV/c
    cuts.AddCut("80000113", "00200019217000008260400000", "0162103500900000");  //new standard eta=0.9 y=0.8  //single pT 0.1 GeV/c
  } else if (trainConfig == 8) {
    cuts.AddCut("80000123", "00100009217000008260400000", "0162103500900000");  //new standard eta=0.9 y=0.8  // minR 2.8
    cuts.AddCut("80000123", "00900009217000008260400000", "0162103500900000");  //new standard eta=0.9 y=0.8  //minR 7.5
    cuts.AddCut("80000123", "00200079217000008260400000", "0162103500900000");  //new standard eta=0.9 y=0.8  //single pT 0. GeV/c
    cuts.AddCut("80000123", "00200019217000008260400000", "0162103500900000");  //new standard eta=0.9 y=0.8  //single pT 0.1 GeV/c
  } else if (trainConfig == 9) {
    cuts.AddCut("80000113", "00200008217000008260400000", "0162103500900000");  //new standard eta=0.9 y=0.8 //TPC Cluster 0.35
    cuts.AddCut("80000113", "00200006217000008260400000", "0162103500900000");  //new standard eta=0.9 y=0.8 //TPC Cluster 0.7
    cuts.AddCut("80000113", "00200009317000008260400000", "0162103500900000");  //new standard eta=0.9 y=0.8 //edEdx -4,5
    cuts.AddCut("80000113", "00200009617000008260400000", "0162103500900000");  //new standard eta=0.9 y=0.8 //edEdx -2.5,4
  } else if (trainConfig == 10) {
    cuts.AddCut("80000123", "00200008217000008260400000", "0162103500900000");  //new standard eta=0.9 y=0.8 //TPC Cluster 0.35
    cuts.AddCut("80000123", "00200006217000008260400000", "0162103500900000");  //new standard eta=0.9 y=0.8 //TPC Cluster 0.7
    cuts.AddCut("80000123", "00200009317000008260400000", "0162103500900000");  //new standard eta=0.9 y=0.8 //edEdx -4,5
    cuts.AddCut("80000123", "00200009617000008260400000", "0162103500900000");  //new standard eta=0.9 y=0.8 //edEdx -2.5,4
  } else if (trainConfig == 11) {
    cuts.AddCut("80000113", "00200009227000008260400000", "0162103500900000");  //new standard eta=0.9 y=0.8 // pi dEdx 1,-10
    cuts.AddCut("80000113", "00200009257000008260400000", "0162103500900000");  //new standard eta=0.9 y=0.8 // pi dEdx 2,-10
    cuts.AddCut("80000113", "00200009216000008260400000", "0162103500900000");  //new standard eta=0.9 y=0.8 //pi minMom 0.25
    cuts.AddCut("80000113", "00200009215000008260400000", "0162103500900000");  //new standard eta=0.9 y=0.8 //pi min Mom 0.3
  } else if (trainConfig == 12) {
    cuts.AddCut("80000123", "00200009227000008260400000", "0162103500900000");  //new standard eta=0.9 y=0.8 // pi dEdx 1,-10
    cuts.AddCut("80000123", "00200009257000008260400000", "0162103500900000");  //new standard eta=0.9 y=0.8 // pi dEdx 2,-10
    cuts.AddCut("80000123", "00200009216000008260400000", "0162103500900000");  //new standard eta=0.9 y=0.8 //pi minMom 0.25
    cuts.AddCut("80000123", "00200009215000008260400000", "0162103500900000");  //new standard eta=0.9 y=0.8 //pi min Mom 0.3
  } else if (trainConfig == 13) {
    cuts.AddCut("80000113", "00200009217200008260400000", "0162103500900000");  //new standard eta=0.9 y=0.8 //pi maxMom 4GeV
    cuts.AddCut("80000113", "00200009216200008260400000", "0162103500900000");  //new standard eta=0.9 y=0.8 // pi minMom 0.25GeV maxMom 4GeV
    cuts.AddCut("80000113", "00200009226000008260400000", "0162103500900000");  //new standard eta=0.9 y=0.8 //pi dEdx 1,-10  minMom 0.25GeV
    cuts.AddCut("80000113", "00200009226200008260400000", "0162103500900000");  //new standard eta=0.9 y=0.8 //pi dEdx 1,-10  minMom 0.25GeV maxMom 4GeV
  } else if (trainConfig == 14) {
    cuts.AddCut("80000123", "00200009217200008260400000", "0162103500900000");  //new standard eta=0.9 y=0.8 //pi maxMom 4GeV
    cuts.AddCut("80000123", "00200009216200008260400000", "0162103500900000");  //new standard eta=0.9 y=0.8 // pi minMom 0.25GeV maxMom 4GeV
    cuts.AddCut("80000123", "00200009226000008260400000", "0162103500900000");  //new standard eta=0.9 y=0.8 //pi dEdx 1,-10  minMom 0.25GeV
    cuts.AddCut("80000123", "00200009226200008260400000", "0162103500900000");  //new standard eta=0.9 y=0.8 //pi dEdx 1,-10  minMom 0.25GeV maxMom 4GeV
  } else if (trainConfig == 15) {
    cuts.AddCut("80000113", "00200009217000003260400000", "0162103500900000");  //new standard eta=0.9 y=0.8 //qT 0.05 1D
    cuts.AddCut("80000113", "00200009217000009260400000", "0162103500900000");  //new standard eta=0.9 y=0.8 //qT 0.03 2D
    cuts.AddCut("80000113", "00200009217000002000260400", "0162103500900000");  //new standard eta=0.9 y=0.8 //qT 0.07 1D
    cuts.AddCut("80000113", "00200009217000008220400000", "0162103500900000");  //new standard eta=0.9 y=0.8 //chi2 and psi Pair 1D
  } else if (trainConfig == 16) {
    cuts.AddCut("80000123", "00200009217000003260400000", "0162103500900000");  //new standard eta=0.9 y=0.8 //qT 0.05 1D
    cuts.AddCut("80000123", "00200009217000009260400000", "0162103500900000");  //new standard eta=0.9 y=0.8 //qT 0.03 2D
    cuts.AddCut("80000123", "00200009217000002000260400", "0162103500900000");  //new standard eta=0.9 y=0.8 //qT 0.07 1D
    cuts.AddCut("80000123", "00200009217000008220400000", "0162103500900000");  //new standard eta=0.9 y=0.8 //chi2 and psi Pair 1D
  } else if (trainConfig == 17) {
    cuts.AddCut("80000113", "00200009217000008160400000", "0162103500900000");  //new standard eta=0.9 y=0.8 //chi2 50  2D
    cuts.AddCut("80000113", "00200009217000008860400000", "0162103500900000");  //new standard eta=0.9 y=0.8 //chi2 20  2D
    cuts.AddCut("80000113", "00200009217000008250400000", "0162103500900000");  //new standard eta=0.9 y=0.8 //psi pair 0.1
    cuts.AddCut("80000113", "00200009217000008270400000", "0162103500900000");  //new standard eta=0.9 y=0.8 //psi pair 0.035
  } else if (trainConfig == 18) {
    cuts.AddCut("80000123", "00200009217000008160400000", "0162103500900000");  //new standard eta=0.9 y=0.8 //chi2 50  2D
    cuts.AddCut("80000123", "00200009217000008860400000", "0162103500900000");  //new standard eta=0.9 y=0.8 //chi2 20  2D
    cuts.AddCut("80000123", "00200009217000008250400000", "0162103500900000");  //new standard eta=0.9 y=0.8 //psi pair 0.1
    cuts.AddCut("80000123", "00200009217000008270400000", "0162103500900000");  //new standard eta=0.9 y=0.8 //psi pair 0.035
  } else if (trainConfig == 19) {
    cuts.AddCut("80000113", "00200009217000008260300000", "0162103500900000");  //new standard eta=0.9 y=0.8 // cos pointing angle 0.75
    cuts.AddCut("80000113", "00200009217000008260600000", "0162103500900000");  //new standard eta=0.9 y=0.8 //cos pointing angle 0.9
    cuts.AddCut("80000113", "00200009217000008260400000", "0162106500900000");  //new standard eta=0.9 y=0.8 // alpha meson cut 0.8
    cuts.AddCut("80000113", "00200009217000008260400000", "0162103400900000");  //new standard eta=0.9 y=0.8 //chi2 meson 500
  } else if (trainConfig == 20) {
    cuts.AddCut("80000123", "00200009217000008260300000", "0162103500900000");  //new standard eta=0.9 y=0.8 // cos pointing angle 0.75
    cuts.AddCut("80000123", "00200009217000008260600000", "0162103500900000");  //new standard eta=0.9 y=0.8 //cos pointing angle 0.9
    cuts.AddCut("80000123", "00200009217000008260400000", "0162106500900000");  //new standard eta=0.9 y=0.8 // alpha meson cut 0.8
    cuts.AddCut("80000123", "00200009217000008260400000", "0162103400900000");  //new standard eta=0.9 y=0.8 //chi2 meson 500
  } else if (trainConfig == 21) {
    cuts.AddCut("80000113", "00200009217000008260400000", "0262103500900000");  //new standard eta=0.9 y=0.8 // BG track multiplicity
    cuts.AddCut("80000113", "00200009217000008260400000", "0162103500000000");  //new standard eta=0.9 y=0.8 //no MCP smearing
    cuts.AddCut("80000113", "00200009317200003290000000", "0162103500900000");  //old standard eta=0.9 y=0.8
    cuts.AddCut("80000113", "00200009217000008260400000", "0162103500800000");  //new standard eta=0.9 y=0.8 // fPSigSmearingCte=0.014;
  } else if (trainConfig == 22) {
    cuts.AddCut("80000123", "00200009217000008260400000", "0262103500900000");  //new standard eta=0.9 y=0.8 // BG track multiplicity
    cuts.AddCut("80000123", "00200009217000008260400000", "0162103500000000");  //new standard eta=0.9 y=0.8 //no MCP smearing
    cuts.AddCut("80000123", "00200009317200003290000000", "0162103500900000");  //old standard eta=0.9 y=0.8
    cuts.AddCut("80000123", "00200009217000008260400000", "0162103500800000");  //new standard eta=0.9 y=0.8 // fPSigSmearingCte=0.014;
  } else if (trainConfig == 23){
    cuts.AddCut("80000113", "00200009217000008260400000", "0162103500900000"); // standard cut
    cuts.AddCut("80000123", "00200009217000008260400000", "0162103500900000"); // standard cut + added signals
    
  } else if (trainConfig == 120) {
    cuts.AddCut("80200113", "00200009300000008260400000", "0162103500000000"); // -3 < kappa <  5
    cuts.AddCut("80200113", "00200009500000008260400000", "0162103500000000"); // -5 < kappa < 10
    cuts.AddCut("80200113", "00200009600000008260400000", "0162103500000000"); // -3 < kappa < 10
    cuts.AddCut("80200113", "00200009700000008260400000", "0162103500000000"); //  0 < kappa < 10
  } else if (trainConfig == 121) {
    cuts.AddCut("82400113", "00200009300000008260400000", "0162103500000000"); // -3 < kappa <  5
    cuts.AddCut("82400113", "00200009500000008260400000", "0162103500000000"); // -5 < kappa < 10
    cuts.AddCut("82400113", "00200009600000008260400000", "0162103500000000"); // -3 < kappa < 10
    cuts.AddCut("82400113", "00200009700000008260400000", "0162103500000000"); //  0 < kappa < 10
  } else if (trainConfig == 122) {
    cuts.AddCut("84600113", "00200009300000008260400000", "0162103500000000"); // -3 < kappa <  5
    cuts.AddCut("84600113", "00200009500000008260400000", "0162103500000000"); // -5 < kappa < 10
    cuts.AddCut("84600113", "00200009600000008260400000", "0162103500000000"); // -3 < kappa < 10
    cuts.AddCut("84600113", "00200009700000008260400000", "0162103500000000"); //  0 < kappa < 10
  } else if (trainConfig == 123) {
    cuts.AddCut("86000113", "00200009300000008260400000", "0162103500000000"); // -3 < kappa <  5
    cuts.AddCut("86000113", "00200009500000008260400000", "0162103500000000"); // -5 < kappa < 10
    cuts.AddCut("86000113", "00200009600000008260400000", "0162103500000000"); // -3 < kappa < 10
    cuts.AddCut("86000113", "00200009700000008260400000", "0162103500000000"); //  0 < kappa < 10

  } else if (trainConfig == 130) { // default cut
    cuts.AddCut("80000113", "00200009217000008260400000", "0162103500000000"); //clean cuts
    cout << "setting cut 130" << endl;
  //--------- purity studies (kappa cut)    -------------------------//
  } else if (trainConfig == 131) {
    cuts.AddCut("80000113", "00200009300000008260400000", "0162103500000000"); // -3 < kappa <  5
    cuts.AddCut("80000113", "00200009500000008260400000", "0162103500000000"); // -5 < kappa < 10
    cuts.AddCut("80000113", "00200009600000008260400000", "0162103500000000"); // -3 < kappa < 10
    cuts.AddCut("80000113", "00200009700000008260400000", "0162103500000000"); //  0 < kappa < 10
  } else if (trainConfig == 132) {
    cuts.AddCut("80200113", "00200009217000008260400000", "0162103500000000"); // 0-20%
    cuts.AddCut("82400113", "00200009217000008260400000", "0162103500000000"); // 20-40%
    cuts.AddCut("84600113", "00200009217000008260400000", "0162103500000000"); // 40-60%
    cuts.AddCut("86000113", "00200009217000008260400000", "0162103500000000"); // 60-100%
  
  // default cuts MB without MC smearing    
  } else if (trainConfig == 133) {
    cuts.AddCut("80000113", "03200009217000008260400000", "0162303500000000"); //New STANDARD CUT |eta| < 0.65, |y| < 0.6
    cuts.AddCut("80000113", "04200009217000008260400000", "0162203500000000"); //New STANDARD CUT |eta| < 0.75, |y| < 0.7
    cuts.AddCut("80000113", "01200009217000008260400000", "0162403500000000"); //New STANDARD CUT |eta| < 0.6, |y| < 0.5
    cuts.AddCut("80000113", "00200009217000003260000000", "0162103500000000"); //tighten Psi pair and chi2 in 2D
  } else if (trainConfig == 134) {
    cuts.AddCut("80000123", "03200009217000008260400000", "0162303500000000"); //New STANDARD CUT |eta| < 0.65, |y| < 0.6
    cuts.AddCut("80000123", "04200009217000008260400000", "0162203500000000"); //New STANDARD CUT |eta| < 0.75, |y| < 0.7
    cuts.AddCut("80000123", "01200009217000008260400000", "0162403500000000"); //New STANDARD CUT |eta| < 0.6, |y| < 0.5
    cuts.AddCut("80000123", "00200009217000003260000000", "0162103500000000"); //tighten Psi pair and chi2 in 2D
  } else if (trainConfig == 135) {
    cuts.AddCut("80000113", "00200009217000008220000000", "0162103500000000"); //tighten psi pair and qt in 2D
    cuts.AddCut("80000113", "00200009217000008260000000", "0162103500000000"); //tighten psi pair and chi2 in 2D and qt in 2D
    cuts.AddCut("80000113", "00200009217000008220400000", "0162103500000000"); //clean cuts
    cuts.AddCut("80000113", "00200009217000008260400000", "0162103500000000"); //clean cuts
  } else if (trainConfig == 136) {
    cuts.AddCut("80000123", "00200009217000008220000000", "0162103500000000"); //tighten psi pair and qt in 2D
    cuts.AddCut("80000123", "00200009217000008260000000", "0162103500000000"); //tighten psi pair and chi2 in 2D and qt in 2D
    cuts.AddCut("80000123", "00200009217000008220400000", "0162103500000000"); //clean cuts
    cuts.AddCut("80000123", "00200009217000008260400000", "0162103500000000"); //clean cuts
  } else if (trainConfig == 137) {
    cuts.AddCut("80000113", "00100009217000008260400000", "0162103500000000");  //new standard eta=0.9 y=0.8  // minR 2.8
    cuts.AddCut("80000113", "00900009217000008260400000", "0162103500000000");  //new standard eta=0.9 y=0.8  //minR 7.5
    cuts.AddCut("80000113", "00200079217000008260400000", "0162103500000000");  //new standard eta=0.9 y=0.8  //single pT 0. GeV/c
    cuts.AddCut("80000113", "00200019217000008260400000", "0162103500000000");  //new standard eta=0.9 y=0.8  //single pT 0.1 GeV/c
  } else if (trainConfig == 138) {
    cuts.AddCut("80000123", "00100009217000008260400000", "0162103500000000");  //new standard eta=0.9 y=0.8  // minR 2.8
    cuts.AddCut("80000123", "00900009217000008260400000", "0162103500000000");  //new standard eta=0.9 y=0.8  //minR 7.5
    cuts.AddCut("80000123", "00200079217000008260400000", "0162103500000000");  //new standard eta=0.9 y=0.8  //single pT 0. GeV/c
    cuts.AddCut("80000123", "00200019217000008260400000", "0162103500000000");  //new standard eta=0.9 y=0.8  //single pT 0.1 GeV/c
  } else if (trainConfig == 139) {
    cuts.AddCut("80000113", "00200008217000008260400000", "0162103500000000");  //new standard eta=0.9 y=0.8 //TPC Cluster 0.35
    cuts.AddCut("80000113", "00200006217000008260400000", "0162103500000000");  //new standard eta=0.9 y=0.8 //TPC Cluster 0.7
    cuts.AddCut("80000113", "00200009317000008260400000", "0162103500000000");  //new standard eta=0.9 y=0.8 //edEdx -4,5
    cuts.AddCut("80000113", "00200009617000008260400000", "0162103500000000");  //new standard eta=0.9 y=0.8 //edEdx -2.5,4
  } else if (trainConfig == 140) {
    cuts.AddCut("80000123", "00200008217000008260400000", "0162103500000000");  //new standard eta=0.9 y=0.8 //TPC Cluster 0.35
    cuts.AddCut("80000123", "00200006217000008260400000", "0162103500000000");  //new standard eta=0.9 y=0.8 //TPC Cluster 0.7
    cuts.AddCut("80000123", "00200009317000008260400000", "0162103500000000");  //new standard eta=0.9 y=0.8 //edEdx -4,5
    cuts.AddCut("80000123", "00200009617000008260400000", "0162103500000000");  //new standard eta=0.9 y=0.8 //edEdx -2.5,4
  } else if (trainConfig == 141) {
    cuts.AddCut("80000113", "00200009227000008260400000", "0162103500000000");  //new standard eta=0.9 y=0.8 // pi dEdx 1,-10
    cuts.AddCut("80000113", "00200009257000008260400000", "0162103500000000");  //new standard eta=0.9 y=0.8 // pi dEdx 2,-10
    cuts.AddCut("80000113", "00200009216000008260400000", "0162103500000000");  //new standard eta=0.9 y=0.8 //pi minMom 0.25
    cuts.AddCut("80000113", "00200009215000008260400000", "0162103500000000");  //new standard eta=0.9 y=0.8 //pi min Mom 0.3
  } else if (trainConfig == 142) {
    cuts.AddCut("80000123", "00200009227000008260400000", "0162103500000000");  //new standard eta=0.9 y=0.8 // pi dEdx 1,-10
    cuts.AddCut("80000123", "00200009257000008260400000", "0162103500000000");  //new standard eta=0.9 y=0.8 // pi dEdx 2,-10
    cuts.AddCut("80000123", "00200009216000008260400000", "0162103500000000");  //new standard eta=0.9 y=0.8 //pi minMom 0.25
    cuts.AddCut("80000123", "00200009215000008260400000", "0162103500000000");  //new standard eta=0.9 y=0.8 //pi min Mom 0.3
  } else if (trainConfig == 143) {
    cuts.AddCut("80000113", "00200009217200008260400000", "0162103500000000");  //new standard eta=0.9 y=0.8 //pi maxMom 4GeV
    cuts.AddCut("80000113", "00200009216200008260400000", "0162103500000000");  //new standard eta=0.9 y=0.8 // pi minMom 0.25GeV maxMom 4GeV
    cuts.AddCut("80000113", "00200009226000008260400000", "0162103500000000");  //new standard eta=0.9 y=0.8 //pi dEdx 1,-10  minMom 0.25GeV
    cuts.AddCut("80000113", "00200009226200008260400000", "0162103500000000");  //new standard eta=0.9 y=0.8 //pi dEdx 1,-10  minMom 0.25GeV maxMom 4GeV
  } else if (trainConfig == 144) {
    cuts.AddCut("80000123", "00200009217200008260400000", "0162103500000000");  //new standard eta=0.9 y=0.8 //pi maxMom 4GeV
    cuts.AddCut("80000123", "00200009216200008260400000", "0162103500000000");  //new standard eta=0.9 y=0.8 // pi minMom 0.25GeV maxMom 4GeV
    cuts.AddCut("80000123", "00200009226000008260400000", "0162103500000000");  //new standard eta=0.9 y=0.8 //pi dEdx 1,-10  minMom 0.25GeV
    cuts.AddCut("80000123", "00200009226200008260400000", "0162103500000000");  //new standard eta=0.9 y=0.8 //pi dEdx 1,-10  minMom 0.25GeV maxMom 4GeV
  } else if (trainConfig == 145) {
    cuts.AddCut("80000113", "00200009217000003260400000", "0162103500000000");  //new standard eta=0.9 y=0.8 //qT 0.05 1D
    cuts.AddCut("80000113", "00200009217000009260400000", "0162103500000000");  //new standard eta=0.9 y=0.8 //qT 0.03 2D
    cuts.AddCut("80000113", "00200009217000002000260400", "0162103500000000");  //new standard eta=0.9 y=0.8 //qT 0.07 1D
    cuts.AddCut("80000113", "00200009217000008220400000", "0162103500000000");  //new standard eta=0.9 y=0.8 //chi2 and psi Pair 1D
  } else if (trainConfig == 146) {
    cuts.AddCut("80000123", "00200009217000003260400000", "0162103500000000");  //new standard eta=0.9 y=0.8 //qT 0.05 1D
    cuts.AddCut("80000123", "00200009217000009260400000", "0162103500000000");  //new standard eta=0.9 y=0.8 //qT 0.03 2D
    cuts.AddCut("80000123", "00200009217000002000260400", "0162103500000000");  //new standard eta=0.9 y=0.8 //qT 0.07 1D
    cuts.AddCut("80000123", "00200009217000008220400000", "0162103500000000");  //new standard eta=0.9 y=0.8 //chi2 and psi Pair 1D
  } else if (trainConfig == 147) {
    cuts.AddCut("80000113", "00200009217000008160400000", "0162103500000000");  //new standard eta=0.9 y=0.8 //chi2 50  2D
    cuts.AddCut("80000113", "00200009217000008860400000", "0162103500000000");  //new standard eta=0.9 y=0.8 //chi2 20  2D
    cuts.AddCut("80000113", "00200009217000008250400000", "0162103500000000");  //new standard eta=0.9 y=0.8 //psi pair 0.1
    cuts.AddCut("80000113", "00200009217000008270400000", "0162103500000000");  //new standard eta=0.9 y=0.8 //psi pair 0.035
  } else if (trainConfig == 148) {
    cuts.AddCut("80000123", "00200009217000008160400000", "0162103500000000");  //new standard eta=0.9 y=0.8 //chi2 50  2D
    cuts.AddCut("80000123", "00200009217000008860400000", "0162103500000000");  //new standard eta=0.9 y=0.8 //chi2 20  2D
    cuts.AddCut("80000123", "00200009217000008250400000", "0162103500000000");  //new standard eta=0.9 y=0.8 //psi pair 0.1
    cuts.AddCut("80000123", "00200009217000008270400000", "0162103500000000");  //new standard eta=0.9 y=0.8 //psi pair 0.035
  } else if (trainConfig == 149) {
    cuts.AddCut("80000113", "00200009217000008260300000", "0162103500000000");  //new standard eta=0.9 y=0.8 // cos pointing angle 0.75
    cuts.AddCut("80000113", "00200009217000008260600000", "0162103500000000");  //new standard eta=0.9 y=0.8 //cos pointing angle 0.9
    cuts.AddCut("80000113", "00200009217000008260400000", "0162106500000000");  //new standard eta=0.9 y=0.8 // alpha meson cut 0.8
    cuts.AddCut("80000113", "00200009217000008260400000", "0162103400000000");  //new standard eta=0.9 y=0.8 //chi2 meson 500
  } else if (trainConfig == 150) {
    cuts.AddCut("80000123", "00200009217000008260300000", "0162103500000000");  //new standard eta=0.9 y=0.8 // cos pointing angle 0.75
    cuts.AddCut("80000123", "00200009217000008260600000", "0162103500000000");  //new standard eta=0.9 y=0.8 //cos pointing angle 0.9
    cuts.AddCut("80000123", "00200009217000008260400000", "0162106500000000");  //new standard eta=0.9 y=0.8 // alpha meson cut 0.8
    cuts.AddCut("80000123", "00200009217000008260400000", "0162103400000000");  //new standard eta=0.9 y=0.8 //chi2 meson 500
  } else if (trainConfig == 151) {
    cuts.AddCut("80000113", "00200009217000008260400000", "0262103500000000");  //new standard eta=0.9 y=0.8 // BG track multiplicity
    cuts.AddCut("80000113", "00200009217000008260400000", "0162103500000000");  //new standard eta=0.9 y=0.8 //no MCP smearing
    cuts.AddCut("80000113", "00200009317200003290000000", "0162103500000000");  //old standard eta=0.9 y=0.8
    cuts.AddCut("80000113", "00200009217000008260400000", "0162103500800000");  //new standard eta=0.9 y=0.8 // fPSigSmearingCte=0.014;
  } else if (trainConfig == 152) {
    cuts.AddCut("80000123", "00200009217000008260400000", "0262103500000000");  //new standard eta=0.9 y=0.8 // BG track multiplicity
    cuts.AddCut("80000123", "00200009217000008260400000", "0162103500000000");  //new standard eta=0.9 y=0.8 //no MCP smearing
    cuts.AddCut("80000123", "00200009317200003290000000", "0162103500000000");  //old standard eta=0.9 y=0.8
    cuts.AddCut("80000123", "00200009217000008260400000", "0162103500800000");  //new standard eta=0.9 y=0.8 // fPSigSmearingCte=0.014;
  
      
  } else if (trainConfig == 173) {
    cuts.AddCut("80000113", "00200009217000008260400000", "0162103500900000");  // all Photon Qualities
    cuts.AddCut("80000113", "00200009217000008260420000", "0162103500900000");  // only photon quality 1
    cuts.AddCut("80000113", "00200009217000008260430000", "0162103500900000");  // only photon quality 2
    cuts.AddCut("80000113", "00200009217000008260440000", "0162103500900000");  // only photon quality 3
  } else if (trainConfig == 174) {
    cuts.AddCut("80000123", "00200009217000008260400000", "0162103500900000");  // all Photon Qualities
    cuts.AddCut("80000123", "00200009217000008260420000", "0162103500900000");  // only photon quality 1
    cuts.AddCut("80000123", "00200009217000008260430000", "0162103500900000");  // only photon quality 2
    cuts.AddCut("80000123", "00200009217000008260440000", "0162103500900000");  // only photon quality 3
  } else if (trainConfig == 175) {
    cuts.AddCut("80000113", "00200009217000008260410000", "0162103500900000");  //no shared electron cut on
    cuts.AddCut("80000113", "00200019217000008260400000", "0162103500900000");  //single pt cut 0.1
    cuts.AddCut("80000113", "00200029217000008260400000", "0162103500900000");  //single pt cut 0.15
    cuts.AddCut("80000113", "00200049217000008260400000", "0162103500900000");  //single pt cut 0.075
  } else if (trainConfig == 176) {
    cuts.AddCut("80000123", "00200009217000008260410000", "0162103500900000");  //no shared electron cut on
    cuts.AddCut("80000123", "00200019217000008260400000", "0162103500900000");  //single pt cut 0.1
    cuts.AddCut("80000123", "00200029217000008260400000", "0162103500900000");  //single pt cut 0.15
    cuts.AddCut("80000123", "00200049217000008260400000", "0162103500900000");  //single pt cut 0.075
  } else if (trainConfig == 177) {
    cuts.AddCut("80000113", "00700009217000008260400000", "0162103500900000");  // min R =35 cm
    cuts.AddCut("80000113", "00700009217000008260420000", "0162103500900000");  // min R =35 cm & only photon quality 1
    cuts.AddCut("80000113", "00700009217000008260430000", "0162103500900000");  // min R =35 cm & only photon quality 2
    cuts.AddCut("80000113", "00700009217000008260440000", "0162103500900000");  // min R =35 cm & only photon quality 3
  } else if (trainConfig == 178) {
    cuts.AddCut("80000123", "00700009217000008260400000", "0162103500900000");  // min R =35 cm
    cuts.AddCut("80000123", "00700009217000008260420000", "0162103500900000");  // min R =35 cm & only photon quality 1
    cuts.AddCut("80000123", "00700009217000008260430000", "0162103500900000");  // min R =35 cm & only photon quality 2
    cuts.AddCut("80000123", "00700009217000008260440000", "0162103500900000");  // min R =35 cm & only photon quality 3
  } else if (trainConfig == 179) {
    cuts.AddCut("80000113", "00200009217000008260400000", "0162803500900000"); //new standard rapidity 0-0.5 in cms
    cuts.AddCut("80000113", "00500009217000008260400000", "0162103500900000"); //new standard RCut=10cm
    cuts.AddCut("80000113", "00800009217000008260400000", "0162103500900000"); //new standard RCut=12.5cm
    cuts.AddCut("80000113", "00600009217000008260400000", "0162103500900000"); //new standard RCut=20cm
  } else if (trainConfig == 180) {
    cuts.AddCut("80000123", "00200009217000008260400000", "0162803500900000"); //new standard rapidity 0-0.5 in cms
    cuts.AddCut("80000123", "00500009217000008260400000", "0162103500900000"); //new standard RCut=10cm
    cuts.AddCut("80000123", "00800009217000008260400000", "0162103500900000"); //new standard RCut=12.5cm
    cuts.AddCut("80000123", "00600009217000008260400000", "0162103500900000"); //new standard RCut=20cm
  } else if (trainConfig == 181) {
    cuts.AddCut("80000113", "00200009297002008250400000", "0152506500000000");
    cuts.AddCut("80000113", "00200009297002008260400000", "0152506500000000");
    cuts.AddCut("80000113", "00200009297002008270400000", "0152506500000000");
    cuts.AddCut("80000113", "00200009297002008250000000", "0152506500000000");
  } else if (trainConfig == 182) {
    cuts.AddCut("80000123", "00200009297002008250400000", "0152506500000000");
    cuts.AddCut("80000123", "00200009297002008260400000", "0152506500000000");
    cuts.AddCut("80000123", "00200009297002008270400000", "0152506500000000");
    cuts.AddCut("80000123", "00200009297002008250000000", "0152506500000000");
  } else if (trainConfig == 183) {
    cuts.AddCut("80000113", "00000009217000008260440000", "0162103500900000");  // min R =0 cm & only photon quality 3
    cuts.AddCut("80000113", "00100009217000008260440000", "0162103500900000");  // min R =2.8 cm & only photon quality 3
    cuts.AddCut("80000113", "00200009217000008260440000", "0162103500900000");  // min R =5 cm & only photon quality 3
    cuts.AddCut("80000113", "00900009217000008260440000", "0162103500900000");  // min R =7.5 cm & only photon quality 3 
  } else if (trainConfig == 184) {
    cuts.AddCut("80000123", "00000009217000008260440000", "0162103500900000");  // min R =0 cm & only photon quality 3
    cuts.AddCut("80000123", "00100009217000008260440000", "0162103500900000");  // min R =2.8 cm & only photon quality 3
    cuts.AddCut("80000123", "00200009217000008260440000", "0162103500900000");  // min R =5 cm & only photon quality 3
    cuts.AddCut("80000123", "00900009217000008260440000", "0162103500900000");  // min R =7.5 cm & only photon quality 3
  } else if (trainConfig == 185) {
    cuts.AddCut("80000113", "00500009217000008260440000", "0162103500900000");  // min R =10 cm & only photon quality 3
    cuts.AddCut("80000113", "00800009217000008260440000", "0162103500900000");  // min R =12.5 cm & only photon quality 3
    cuts.AddCut("80000113", "00600009217000008260440000", "0162103500900000");  // min R =20 cm & only photon quality 3
    cuts.AddCut("80000113", "00700009217000008260440000", "0162103500900000");  // min R =35 cm & only photon quality 3 
  } else if (trainConfig == 186) {
    cuts.AddCut("80000123", "00500009217000008260440000", "0162103500900000");  // min R =10 cm & only photon quality 3
    cuts.AddCut("80000123", "00800009217000008260440000", "0162103500900000");  // min R =12.5 cm & only photon quality 3
    cuts.AddCut("80000123", "00600009217000008260440000", "0162103500900000");  // min R =20 cm & only photon quality 3
    cuts.AddCut("80000123", "00700009217000008260440000", "0162103500900000");  // min R =35 cm & only photon quality 3  
  } else if (trainConfig == 187) {
    cuts.AddCut("80000113", "00000009217000008260400000", "0162103500900000");  // min R =0 cm
    cuts.AddCut("80000113", "00100009217000008260400000", "0162103500900000");  // min R =2.8 cm
    cuts.AddCut("80000113", "00200009217000008260400000", "0162103500900000");  // min R =5 cm
    cuts.AddCut("80000113", "00900009217000008260400000", "0162103500900000");  // min R =7.5 cm
  } else if (trainConfig == 188) {
    cuts.AddCut("80000123", "00000009217000008260400000", "0162103500900000");  // min R =0 cm
    cuts.AddCut("80000123", "00100009217000008260400000", "0162103500900000");  // min R =2.8 cm
    cuts.AddCut("80000123", "00200009217000008260400000", "0162103500900000");  // min R =5 cm
    cuts.AddCut("80000123", "00900009217000008260400000", "0162103500900000");  // min R =7.5 cm 
  } else if (trainConfig == 189) {
    cuts.AddCut("80000113", "00200009217000008260404000", "0162103500900000");  // standard dc cut 4
    cuts.AddCut("80000113", "00200009217000008260404000", "0162101500900000");  //  standard dc cut 4 alpha cut 1
    cuts.AddCut("80000013", "00200009217000008260404000", "0162103500900000");  // standard dc cut 4 no Pileup cut
    cuts.AddCut("80000013", "00200009217000008260404000", "0162101500900000");  //  standard dc cut 4 alpha cut 1 no Pileup cut
  } else if (trainConfig == 190) {
    cuts.AddCut("80000123", "00200009217000008260404000", "0162103500900000");  // standard dc cut 4
    cuts.AddCut("80000123", "00200009217000008260404000", "0162101500900000");  // standard dc cut 4 alpha cut 1
    cuts.AddCut("80000023", "00200009217000008260404000", "0162103500900000");  // standard dc cut 4  no Pileup cut
    cuts.AddCut("80000023", "00200009217000008260404000", "0162101500900000");  // standard dc cut 4 alpha cut 1 no Pileup cut
  } else if (trainConfig == 191) {
    cuts.AddCut("80000113", "00200009217000008260404000", "0162103500000000");  // standard dc cut 4//eta meson
    cuts.AddCut("80000113", "00200009217000008260404000", "0162101500000000");  // standard dc cut 4 alpha cut 1
    cuts.AddCut("80000013", "00200009217000008260404000", "0162103500000000");  // standard dc cut 4  no Pileup cut
    cuts.AddCut("80000013", "00200009217000008260404000", "0162101500000000");  // standard dc cut 4 alpha cut 1 no Pileup cut
  } else if (trainConfig == 192) {
    cuts.AddCut("80000123", "00200009217000008260404000", "0162103500000000");  // standard dc cut 4//eta meson
    cuts.AddCut("80000123", "00200009217000008260404000", "0162101500000000");  // standard dc cut 4 alpha cut 1
    cuts.AddCut("80000023", "00200009217000008260404000", "0162103500000000");  // standard dc cut 4  no Pileup cut
    cuts.AddCut("80000023", "00200009217000008260404000", "0162101500000000");  // standard dc cut 4 alpha cut 1 no Pileup cut
  } else if (trainConfig == 193) {
    cuts.AddCut("80000113", "00200009217000008260404000", "0162103500900000");  // standard //pi0 dc cut4
    cuts.AddCut("80000113", "00200009000000008260404000", "0162103500900000");  // open dEdx Cut -10;10
    cuts.AddCut("80000113", "00200009400000008260404000", "0162103500900000");  // open dEdx cut; -6;7
    cuts.AddCut("80000113", "00200009000100008260404000", "0162103500900000");  // no dEdx Cut
  } else if (trainConfig == 194) {
    cuts.AddCut("80000123", "00200009217000008260404000", "0162103500900000");  // standard //pi0 dc cut4
    cuts.AddCut("80000123", "00200009000000008260404000", "0162103500900000");  //  open dEdx Cut -10;10
    cuts.AddCut("80000123", "00200009400000008260404000", "0162103500900000");  //  open dEdx cut; -6;7
    cuts.AddCut("80000123", "00200009000100008260404000", "0162103500900000");  //no dEdx Cut 
  } else if (trainConfig == 195) {
    cuts.AddCut("80000113", "00200009217000008260404000", "0162103500000000");  // standard //eta dc cut4
    cuts.AddCut("80000113", "00200009000000008260404000", "0162103500000000");  // open dEdx Cut -10;10
    cuts.AddCut("80000113", "00200009400000008260404000", "0162103500000000");  // open dEdx cut; -6;7
    cuts.AddCut("80000113", "00200009000100008260404000", "0162103500000000");  // no dEdx Cut
  } else if (trainConfig == 196) {
    cuts.AddCut("80000123", "00200009217000008260404000", "0162103500000000");  // standard //eta dc cut4
    cuts.AddCut("80000123", "00200009000000008260404000", "0162103500000000");  //  open dEdx Cut -10;10
    cuts.AddCut("80000123", "00200009400000008260404000", "0162103500000000");  //  open dEdx cut; -6;7
    cuts.AddCut("80000123", "00200009000100008260404000", "0162103500000000");  //no dEdx Cut
  } else if (trainConfig == 197) {
    cuts.AddCut("80000113", "00200009217000008260404000", "0162101500900000");  // standard // dc cut4 alpha1 
    cuts.AddCut("80052113", "00200009217000008260404000", "0162101500900000");  //  kEMC7
    cuts.AddCut("80083113", "00200009217000008260404000", "0162101500900000");  //  kEMCEG1 based on INT7
    cuts.AddCut("80085113", "00200009217000008260404000", "0162101500900000");  //  kEMCEG2 based on INT7 
  } else if (trainConfig == 198) {
    cuts.AddCut("80000113", "00200009217000008260404000", "0162101500000000");  // standard //eta dc cut4 alpha1 
    cuts.AddCut("80052113", "00200009217000008260404000", "0162101500000000");  //  kEMC7
    cuts.AddCut("80083113", "00200009217000008260404000", "0162101500000000");  //  kEMCEG1 based on INT7
    cuts.AddCut("80085113", "00200009217000008260404000", "0162101500000000");  //  kEMCEG2 based on INT7
  } else if (trainConfig == 199) {   //Stephan weights
    cuts.AddCut("80000113", "00200009217000008260404000", "0162103500900000");  // standard dc cut 4
  } else if (trainConfig == 200) {
    cuts.AddCut("80000123", "00200009217000008260404000", "0162103500900000");  // standard dc cut 4
  } else if (trainConfig == 201) {
    cuts.AddCut("80000113", "00200009217000008260404000", "0162103500900000");  // standard dc cut 4
  } else if (trainConfig == 202) {
    cuts.AddCut("80000123", "00200009217000008260404000", "0162103500900000");  // standard dc cut 4
  } else if (trainConfig == 203) {
    cuts.AddCut("80000113", "00200009217000008260404000", "0162103500900000");  // standard dc cut 4
  } else if (trainConfig == 204) {
    cuts.AddCut("80000123", "00200009217000008260404000", "0162103500900000");  // standard dc cut 4
   
 // -------------------------A. Marin,   2016 pPb, open cuts --------------------------------------

 // Min Bias
  } else if (trainConfig == 300) {  
    cuts.AddCut("80010113", "00200009297302001280004000", "0152103500000000"); // Min Bias
    cuts.AddCut("80010113", "00200009297302001280004000", "0152101500000000"); // Min Bias , alpha pT dependent

    // TRD trigger
  } else if (trainConfig == 301) {  
    cuts.AddCut("80047113", "00200009297302001280004000", "0152103500000000"); // TRD trigger HQU for 8TeV
    cuts.AddCut("80043113", "00200009297302001280004000", "0152103500000000"); // TRD trigger HSE for 8TeV

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
  
  
  Bool_t doWeighting = kFALSE;
  if (doWeightingPart == 1 || doWeightingPart == 2 || doWeightingPart == 3) doWeighting = kTRUE;
  
  EventCutList->SetOwner(kTRUE);
  AliConvEventCuts **analysisEventCuts = new AliConvEventCuts*[numberOfCuts];
  ConvCutList->SetOwner(kTRUE);
  AliConversionPhotonCuts **analysisCuts = new AliConversionPhotonCuts*[numberOfCuts];
  MesonCutList->SetOwner(kTRUE);
  AliConversionMesonCuts **analysisMesonCuts = new AliConversionMesonCuts*[numberOfCuts];
  Bool_t initializedMatBudWeigths_existing    = kFALSE;
  
  if (doWeighting) Printf("weighting has been switched on");
  
  for(Int_t i = 0; i<numberOfCuts; i++){
    cout << "initialization of cutnumber: " << i << endl;
    
    analysisEventCuts[i] = new AliConvEventCuts();
    if ( trainConfig == 1 || trainConfig == 3 || trainConfig == 5 || trainConfig == 7 || trainConfig == 9 || trainConfig == 11 || trainConfig == 13 || trainConfig == 15|| trainConfig == 17||
        trainConfig == 19 || trainConfig == 21 || trainConfig == 130 || trainConfig == 131|| trainConfig == 133 || trainConfig == 135 || trainConfig == 137 || trainConfig == 139 || trainConfig == 141 || trainConfig == 143 || trainConfig == 145 || trainConfig == 147 || trainConfig == 149 || trainConfig == 151 || trainConfig == 173 || trainConfig == 175 || trainConfig == 177 || trainConfig == 179 ||
        trainConfig == 181 || trainConfig == 183 || trainConfig == 185 || trainConfig == 187 || trainConfig == 189 || trainConfig == 191 || trainConfig == 193 || trainConfig == 195 || trainConfig == 197 || trainConfig == 198 || trainConfig == 199 || trainConfig == 201 || trainConfig == 203 ){
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
    if ( trainConfig == 2 || trainConfig == 4 || trainConfig == 6 || trainConfig == 8 || trainConfig == 10 || trainConfig == 12 || trainConfig == 14 || trainConfig == 16|| trainConfig == 18||
        trainConfig == 20|| trainConfig == 22 || trainConfig == 130 || trainConfig == 134 || trainConfig == 136 || trainConfig == 138 || trainConfig == 140 || trainConfig == 142 || trainConfig == 144 || 
        trainConfig == 146 || trainConfig == 148 || trainConfig == 150 || trainConfig == 152 || trainConfig == 174 || trainConfig == 176 || trainConfig == 178 || trainConfig == 180 ||
        trainConfig == 182 || trainConfig == 184 || trainConfig == 186 || trainConfig == 188 || trainConfig == 190 || trainConfig == 192 || trainConfig == 194 || trainConfig == 196 ||
        trainConfig == 200 || trainConfig == 202 || trainConfig == 204){
      if (doWeighting){
        analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13e7_addSig_pPb_5023GeV_MBV0A",
                                      "Eta_Hijing_LHC13e7_addSig_pPb_5023GeV_MBV0A", "","Pi0_Fit_Data_pPb_5023GeV_MBV0A","Eta_Fit_Data_pPb_5023GeV_MBV0A");
      }
    }  
    
    if (trainConfig == 23){
      if (doWeighting){
        if (i == 0){
          if (generatorName.CompareTo("DPMJET")==0){
          analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE, fileNameInputForWeighting, "Pi0_DPMJET_LHC13b2_efix_pPb_5023GeV_MBV0A",
                                        "Eta_DPMJET_LHC13b2_efix_pPb_5023GeV_MBV0A", "","Pi0_Fit_Data_pPb_5023GeV_MBV0A","Eta_Fit_Data_pPb_5023GeV_MBV0A");
          } else if (generatorName.CompareTo("HIJING")==0){
            analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13e7_pPb_5023GeV_MBV0A",
                                          "Eta_Hijing_LHC13e7_pPb_5023GeV_MBV0A", "","Pi0_Fit_Data_pPb_5023GeV_MBV0A","Eta_Fit_Data_pPb_5023GeV_MBV0A");
          }
        }  
        if (i == 1){
          analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13e7_addSig_pPb_5023GeV_MBV0A",
                                      "Eta_Hijing_LHC13e7_addSig_pPb_5023GeV_MBV0A", "","Pi0_Fit_Data_pPb_5023GeV_MBV0A","Eta_Fit_Data_pPb_5023GeV_MBV0A");
        }  
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
    analysisEventCuts[i]->SetV0ReaderName(V0ReaderName);
    if (periodNameV0Reader.CompareTo("") != 0) analysisEventCuts[i]->SetPeriodEnum(periodNameV0Reader);
    analysisEventCuts[i]->SetLightOutput(runLightOutput);
    analysisEventCuts[i]->InitializeCutsFromCutString((cuts.GetEventCut(i)).Data());    
    if (doEtaShiftIndCuts) {
      analysisEventCuts[i]->DoEtaShift(doEtaShiftIndCuts);
      analysisEventCuts[i]->SetEtaShift(stringShift);
    }
    if (trainConfig == 179 || trainConfig == 180) {
      analysisEventCuts[0]->DoEtaShift(kTRUE);
      analysisEventCuts[0]->SetEtaShift(-0.215);
    }
    
    EventCutList->Add(analysisEventCuts[i]);
    analysisEventCuts[i]->SetFillCutHistograms("",kFALSE);
    cout << "initialized event cut: " << (cuts.GetEventCut(i)).Data() << endl;
    
    analysisCuts[i] = new AliConversionPhotonCuts();
    if ( trainConfig == 131 || trainConfig == 120 || trainConfig == 121 || trainConfig == 122 || trainConfig == 123){
      analysisCuts[i]->SetSwitchToKappaInsteadOfNSigdEdxTPC(kTRUE);
    }

    if (enableMatBudWeightsPi0 > 0){
      if (isMC > 0){
        if (analysisCuts[i]->InitializeMaterialBudgetWeights(enableMatBudWeightsPi0,filenameMatBudWeights)){
          initializedMatBudWeigths_existing = kTRUE;}
        else {cout << "ERROR The initialization of the materialBudgetWeights did not work out." << endl;}
      } else {cout << "ERROR 'enableMatBudWeightsPi0'-flag was set > 0 even though this is not a MC task. It was automatically reset to 0." << endl;}
    }
    analysisCuts[i]->SetV0ReaderName(V0ReaderName);
    analysisCuts[i]->SetLightOutput(runLightOutput);
    analysisCuts[i]->InitializeCutsFromCutString((cuts.GetPhotonCut(i)).Data());
    if ((trainConfig == 193 || trainConfig == 194 || trainConfig == 195 || trainConfig == 196 )&& i==3 ) {
      analysisCuts[i]->SetDodEdxSigmaCut(kFALSE);
    }
    analysisCuts[i]->SetIsHeavyIon(isHeavyIon);
    ConvCutList->Add(analysisCuts[i]);
    analysisCuts[i]->SetFillCutHistograms("",kFALSE);
    cout << "initialized photon cut: " << (cuts.GetPhotonCut(i)).Data() << endl;
    
    analysisMesonCuts[i] = new AliConversionMesonCuts();
    if (trainConfig == 189 || trainConfig == 190 || trainConfig == 191 || trainConfig == 192) {
      analysisMesonCuts[i]->SetOpeningAngleCut(0.0);
    }
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
  task->SetDoPhotonQA(enableQAPhotonTask);  //Attention new switch small for Photon QA
  task->SetDoPlotVsCentrality(enablePlotVsCentrality);
  if (initializedMatBudWeigths_existing) {
      task->SetDoMaterialBudgetWeightingOfGammasForTrueMesons(kTRUE);
  }
  
  
  //connect containers
  AliAnalysisDataContainer *coutput =
    mgr->CreateContainer(Form("GammaConvV1_%i",trainConfig), TList::Class(),
              AliAnalysisManager::kOutputContainer,Form("GammaConvV1_%i.root",trainConfig));
  
  mgr->AddTask(task);
  mgr->ConnectInput(task,0,cinput);
  mgr->ConnectOutput(task,1,coutput);
  
  return;
  
}
