/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: Friederike Bock, Lucia Leardini                                *
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
//pp together with all supporting classes,
//this is mainly for running the dca tree or other memory extensive output versions
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
      eventCutArray = new TString[nMaxCuts]; photonCutArray = new TString[nMaxCuts]; mesonCutArray = new TString[nMaxCuts]; clusterCutArray = new TString[nMaxCuts];
      for(Int_t i=0; i<nMaxCuts; i++) {eventCutArray[i] = ""; photonCutArray[i] = ""; mesonCutArray[i] = ""; clusterCutArray[i] = "";}
    }

    void AddCut(TString eventCut, TString photonCut, TString mesonCut){
      if(nCuts>=nMaxCuts) {cout << "ERROR in CutHandlerConv: Exceeded maximum number of cuts!" << endl; validCuts = false; return;}
      if( eventCut.Length()!=8 || photonCut.Length()!=26 || mesonCut.Length()!=16 ) {cout << "ERROR in CutHandlerConv: Incorrect length of cut string!" << endl; validCuts = false; return;}
      eventCutArray[nCuts]=eventCut; photonCutArray[nCuts]=photonCut; mesonCutArray[nCuts]=mesonCut; clusterCutArray[nCuts]="";
      nCuts++;
      return;
    }
    void AddCut(TString eventCut, TString photonCut, TString mesonCut, TString clusterCut){
      if(nCuts>=nMaxCuts) {cout << "ERROR in CutHandlerConv: Exceeded maximum number of cuts!" << endl; validCuts = false; return;}
      if( eventCut.Length()!=8 || photonCut.Length()!=26 || mesonCut.Length()!=16 || clusterCut.Length()!=19 ) {cout << "ERROR in CutHandlerConv: Incorrect length of cut string!" << endl; validCuts = false; return;}
      eventCutArray[nCuts]=eventCut; photonCutArray[nCuts]=photonCut; mesonCutArray[nCuts]=mesonCut; clusterCutArray[nCuts]=clusterCut;
      nCuts++;
      return;
    }

    Bool_t AreValid(){return validCuts;}
    Int_t GetNCuts(){if(validCuts) return nCuts; else return 0;}
    TString GetEventCut(Int_t i){if(validCuts&&i<nMaxCuts&&i>=0) return eventCutArray[i]; else{cout << "ERROR in CutHandlerConv: GetEventCut wrong index i" << endl;return "";}}
    TString GetPhotonCut(Int_t i){if(validCuts&&i<nMaxCuts&&i>=0) return photonCutArray[i]; else {cout << "ERROR in CutHandlerConv: GetPhotonCut wrong index i" << endl;return "";}}
    TString GetMesonCut(Int_t i){if(validCuts&&i<nMaxCuts&&i>=0) return mesonCutArray[i]; else {cout << "ERROR in CutHandlerConv: GetMesonCut wrong index i" << endl;return "";}}
    TString GetClusterCut(Int_t i){if(validCuts&&i<nMaxCuts&&i>=0) return clusterCutArray[i]; else {cout << "ERROR in CutHandlerConv: GetClusterCut wrong index i" << endl;return "";}}
  private:
    Bool_t validCuts;
    Int_t nCuts; Int_t nMaxCuts;
    TString* eventCutArray;
    TString* photonCutArray;
    TString* mesonCutArray;
    TString* clusterCutArray;
};

void AddTask_GammaConvV1_pp2(   Int_t    trainConfig                 = 1,                               //change different set of cuts
                                Int_t    isMC                        = 0,                               //run MC
                                Int_t    enableQAMesonTask           = 0,                               //enable QA in AliAnalysisTaskGammaConvV1
                                Int_t    enableQAPhotonTask          = 0,                               // enable additional QA task
                                TString  fileNameInputForWeighting   = "MCSpectraInput.root",           // path to file for weigting input
                                TString  cutnumberAODBranch          = "000000006008400001001500000",   // cutnumber with which AODs have been filtered
                                Bool_t   enableV0findingEffi         = kFALSE,                          // enables V0finding efficiency histograms
                                Bool_t   enableTriggerMimicking      = kFALSE,                          // enable trigger mimicking
                                Bool_t   enableTriggerOverlapRej     = kFALSE,                          // enable trigger overlap rejection
                                Float_t  maxFacPtHard                = 3.,                              // maximum factor between hardest jet and ptHard generated
                                TString  periodNameV0Reader          = "",                              //
                                Bool_t  doMultiplicityWeighting      = kFALSE,                          //
                                TString fileNameInputForMultWeighing = "Multiplicity.root",             //
                                TString periodNameAnchor             = "",                              //
                                Bool_t   runLightOutput             = kFALSE,                          // switch to run light output (only essential histograms for afterburner)
				Bool_t   enableElecDeDxPostCalibration = kFALSE,
				TString  fileNameElecDeDxPostCalibration = "dEdxCorrectionMap_Period_Pass.root",
                                TString  additionalTrainConfig       = "0"                              // additional counter for trainconfig, this has to be always the last parameter
                           ) {

  Int_t isHeavyIon = 0;


  TString corrTaskSetting = ""; // select which correction task setting to use
  //parse additionalTrainConfig flag
  TObjArray *rAddConfigArr = additionalTrainConfig.Tokenize("_");
  if(rAddConfigArr->GetEntries()<1){cout << "ERROR: AddTask_GammaConvV1_pp2 during parsing of additionalTrainConfig String '" << additionalTrainConfig.Data() << "'" << endl; return;}
  TObjString* rAdditionalTrainConfig;
  for(Int_t i = 0; i<rAddConfigArr->GetEntries() ; i++){
    if(i==0){ rAdditionalTrainConfig = (TObjString*)rAddConfigArr->At(i);
    } else {
      TObjString* temp = (TObjString*) rAddConfigArr->At(i);
      TString tempStr = temp->GetString();
      cout<< tempStr.Data()<<endl;

      if(tempStr.Contains("MaterialBudgetWeights") && enableMatBudWeightsPi0 > 0){
         if(tempStr.Contains("MaterialBudgetWeightsNONE")){
            enableMatBudWeightsPi0 = 0;
            cout << "INFO:  AddTask_GammaConvV1_pp2 materialBudgetWeights switched off signaled by additionalTrainConfigFlag" << endl;
         } else {
            TObjArray *fileNameMatBudWeightsArr = filenameMatBudWeights.Tokenize("/");
            if(fileNameMatBudWeightsArr->GetEntries()<1 ){
                cout<<"ERROR: AddTask_GammaConvV1_pp2 when reading material budget weights file name" << filenameMatBudWeights.Data()<< "'" << endl;
                return;
            }
            TObjString * oldMatObjStr = (TObjString*)fileNameMatBudWeightsArr->At( fileNameMatBudWeightsArr->GetEntries()-1);
            TString  oldfileName  = oldMatObjStr->GetString();
            TString  newFileName  = Form("MCInputFile%s.root",tempStr.Data());
            cout<<newFileName.Data()<<endl;
            if( oldfileName.EqualTo(newFileName.Data()) == 0 ){
              filenameMatBudWeights.ReplaceAll(oldfileName.Data(),newFileName.Data());
              cout << "INFO: AddTask_GammaConvV1_pp2 the material budget weights file has been change to " <<filenameMatBudWeights.Data()<<"'"<< endl;
          }
        }
      } else if(tempStr.BeginsWith("CF")){
        cout << "INFO: AddTask_GammaCalo_pp2 will use custom branch from Correction Framework!" << endl;
        corrTaskSetting = tempStr;
        corrTaskSetting.Replace(0,2,"");
      }
    }
  }

  TString sAdditionalTrainConfig = rAdditionalTrainConfig->GetString();
  if (sAdditionalTrainConfig.Atoi() > 0){
    trainConfig = trainConfig + sAdditionalTrainConfig.Atoi();
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

  TString cutnumberPhoton = "00000008400000000100000000";
  if (  periodNameV0Reader.CompareTo("LHC16f") == 0 || periodNameV0Reader.CompareTo("LHC17g")==0 || periodNameV0Reader.CompareTo("LHC18c")==0 ||
        periodNameV0Reader.CompareTo("LHC17d1") == 0  || periodNameV0Reader.CompareTo("LHC17d12")==0 ||
        periodNameV0Reader.CompareTo("LHC17h3")==0 || periodNameV0Reader.CompareTo("LHC17k1")==0 ||
        periodNameV0Reader.CompareTo("LHC17f8b") == 0 ||
        periodNameV0Reader.CompareTo("LHC16P1JJLowB") == 0 || periodNameV0Reader.CompareTo("LHC16P1Pyt8LowB") == 0 )
    cutnumberPhoton         = "00000088400000000100000000";
  TString cutnumberEvent = "00000003";
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
        if (trainConfig == 21){
          fCuts->SetDodEdxSigmaCut(kFALSE);
        }
      }
    }
    if(inputHandler->IsA()==AliAODInputHandler::Class()){
    // AOD mode
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
  AliAnalysisTaskGammaConvV1 *task=NULL;
  task= new AliAnalysisTaskGammaConvV1(Form("GammaConvV1_%i",trainConfig));
  task->SetIsHeavyIon(isHeavyIon);
  task->SetIsMC(isMC);
  task->SetV0ReaderName(V0ReaderName);
  task->SetLightOutput(runLightOutput);
  // Cut Numbers to use in Analysis

  CutHandlerConv cuts;

  //----------------------------- configuration for 2.76TeV standard cuts ----------------------------------------------------
  if (trainConfig == 1){
    cuts.AddCut("00000113", "00200009397300008250400000", "0163103100900000"); // new pi0/eta cut 2.76TeV
  } else if (trainConfig == 2) {
    cuts.AddCut("00000113", "00200009397300008250400000", "0163103100000000"); // new pi0/eta cut 2.76TeV without MC smearing
  } else if (trainConfig == 3) {
    cuts.AddCut("00000113", "00200009366300003800000000", "0163103100900000"); // standard cut Pi0 pp 2.76TeV PbPb paper 2012
  } else if (trainConfig == 4) {
    cuts.AddCut("00000113", "00200009297002008250400000", "0163103100900000"); // standard cut LHC11h pp 2.76TeV
  } else if (trainConfig == 5) {
    cuts.AddCut("00000113", "00200009227302008250404000", "0163101500000000"); // Ana eta analysis prefered 2.76TeV
  } else if (trainConfig == 6) {
    cuts.AddCut("00000113", "00200009327000008250400000", "0163103100900000"); // go with 1sigma pi rejec to infty 2.76TeV
  } else if (trainConfig == 7) {
    cuts.AddCut("00000113", "00200009317000008250400000", "0163103100900000"); // go with 0sigma pi rejec to infty 2.76TeV
  } else if (trainConfig == 8) {
    cuts.AddCut("00000113", "00200009357000008250400000", "0163103100900000"); // go with 2sigma pi reject to infy 2.76TeV
  } else if (trainConfig == 9) {
    cuts.AddCut("00003113", "00200009397300008250400000", "0163103100900000"); // new pi0/eta cut 2.76TeV wSDD
  } else if (trainConfig == 10) {
    cuts.AddCut("00051013", "00200009397300008250400000", "0163103100900000"); // new pi0/eta cut 2.76TeV wSDD & EMC1

  //----------------------------- configuration for  8 TeV standard  --------------------------------------------------------
  } else if (trainConfig == 20) {
    cuts.AddCut("00010113", "00200009227300008250404000", "0152103500000000"); //standard cut pp 8 TeV
  } else if (trainConfig == 21) {
    cuts.AddCut("00052113", "00200009227302008250400000", "0152103500000000"); //standard cut pp 8 TeV EMC7
  } else if (trainConfig == 22) {
    cuts.AddCut("00081113", "00200009227302008250400000", "0152103500000000"); //standard cut pp 8 TeV EGA
  } else if (trainConfig == 23) {
    cuts.AddCut("00010213", "00200009227300008250404000", "0152103500000000"); //standard cut pp 8 TeV + past future max rejection
  } else if (trainConfig == 24) {
    cuts.AddCut("00010513", "00200009227300008250404000", "0152103500000000"); //standard cut pp 8 TeV + past future medium rejection
  } else if (trainConfig == 25) {
    cuts.AddCut("00010113", "0a200009227300008250404000", "0152103500000000"); //eta cut 0.2 < |eta| < 0.9
  } else if (trainConfig == 26) {
    cuts.AddCut("00010613", "00200009227300008250404000", "0152103500000000"); //V0M vs TPCout cut6
  } else if (trainConfig == 27) {
    cuts.AddCut("00010113", "002000j9227300008250404000", "0152103500000000"); //asym pT cut: 0.100 GeV and 0.075 GeV
  } else if (trainConfig == 28) {
    cuts.AddCut("00010113", "002000l9227300008250404000", "0152103500000000"); //asym pT cut: 0.200 GeV and 0.075 GeV
  } else if (trainConfig == 29) {
    cuts.AddCut("00010113", "002000f9227300008250404000", "0152103500000000"); //gamma pT > 0.2 GeV/c
  //----------------------------- configuration for  7 TeV standard cuts -----------------------------------------------------
  } else if (trainConfig == 30) {
    cuts.AddCut("00000113", "00200009227300008250404000", "0152103500000000"); //New standard cut pp 7 TeV direct photon
  } else if (trainConfig == 31) {
    cuts.AddCut("00000113", "00200009227302008250400000", "0152103500000000"); //standard cut pp 7 TeV
  } else if (trainConfig == 32) {
    cuts.AddCut("00000113", "00200008366300000200000000", "0163103100900000"); //old standard cut pp 7 TeV
  } else if (trainConfig == 33) {
    cuts.AddCut("00000113", "002000c9227300008250404000", "0152103500000000"); // 7 TeV std, but min electron pT > 0.6 for all configs
  } else if (trainConfig == 34) {
    cuts.AddCut("00000113", "002000b9227300008250404000", "0152103500000000"); // gamma pT > 0.1 GeV/c
  } else if (trainConfig == 35) {
    cuts.AddCut("00000113", "002000e9227300008250404000", "0152103500000000"); // gamma pT > 0.15 GeV/c
  } else if (trainConfig == 36) {
    cuts.AddCut("00000113", "002000f9227300008250404000", "0152103500000000"); // gamma pT > 0.2 GeV/c

  //----------------------------- configuration for run 2 analysis 13 TeV ----------------------------------------------------
  } else if (trainConfig == 40){
    cuts.AddCut("00010113", "00200009227302008254404000", "0152101500000000"); //standard cut Gamma pp 13TeV, V0AND
  } else if (trainConfig == 41){
    cuts.AddCut("00074113", "00200009227302008254404000", "0152101500000000"); //standard cut Gamma pp 13TeV, V0 HM
  } else if (trainConfig == 42){
    cuts.AddCut("00075113", "00200009227302008254404000", "0152101500000000"); //standard cut Gamma pp 13TeV, SPD HM
  } else if (trainConfig == 43){
    cuts.AddCut("00010113", "00200009227300008250404000", "0152103500000000"); //New standard cut Gamma Pi0 Eta pp 13TeV, V0AND
  } else if (trainConfig == 44){
    cuts.AddCut("00010113", "00200009266300008854404000", "0152101500000000"); // A. Marin alpha pT dependent and gamma asym cut
  } else if (trainConfig == 45){
    cuts.AddCut("00010113", "00200009267300008254404000", "0152103500000000"); // A. Marin alpha pT dependent and gamma asym cut
  } else if (trainConfig == 46){
    cuts.AddCut("00010113", "00a00009267300008254404000", "0152103500000000"); // A. Marin alpha pT dependent and gamma asym cut
  } else if (trainConfig == 47){
    cuts.AddCut("00010113", "00b00009267300008254404000", "0152103500000000"); // A. Marin alpha pT dependent and gamma asym cut
  } else if (trainConfig == 48){
    cuts.AddCut("00010113", "00c00009267300008254404000", "0152103500000000"); // A. Marin alpha pT dependent and gamma asym cut
  } else if (trainConfig == 49){
    cuts.AddCut("00010113", "00200009227300008250404000", "0163103100000000"); // J. Luehder AOD Compare

  //----------------------------- configuration for run 2 analysis 5 TeV ----------------------------------------------------
  } else if (trainConfig == 50){
    cuts.AddCut("00010113", "00200009227300008250404000", "0152101500000000"); //old standard cut pp 5 TeV VAND
  } else if (trainConfig == 51){
    cuts.AddCut("00010113", "00200009227300008250404000", "0152103500000000"); //new standard cut pp 5 TeV VAND
  } else if (trainConfig == 52){
    cuts.AddCut("00010113", "00a00009227300008250404000", "0152103500000000"); //new standard cut pp 5 TeV VAND
  } else if (trainConfig == 53){
    cuts.AddCut("00010113", "00b00009227300008250404000", "0152103500000000"); //new standard cut pp 5 TeV VAND
  } else if (trainConfig == 54){
    cuts.AddCut("00010113", "00c00009227300008250404000", "0152103500000000"); //new standard cut pp 5 TeV VAND
  } else if (trainConfig == 55){
    cuts.AddCut("00010113", "00200009227300008250704000", "0152103500000000"); //test cosPA scan
  } else if (trainConfig == 56){
    cuts.AddCut("00010113", "00200009227300008250804000", "0152103500000000"); //test cosPA scan
  } else if (trainConfig == 57){
    cuts.AddCut("00010113", "00200009227300008250904000", "0152103500000000"); //test cosPA scan
  } else if (trainConfig == 58){
    cuts.AddCut("00010113", "00200009227300008250a04000", "0152103500000000"); //test cosPA scan
  } else if (trainConfig == 59){
    cuts.AddCut("00010113", "00200009a27300008250904120", "0152103500000000"); //cosPA, 0.99 eta 0.9
  } else if (trainConfig == 60){
    cuts.AddCut("00010113", "0d200009a27300008250904120", "0152103500000000"); //cosPA, 0.99 eta 0.8
  } else if (trainConfig == 61){
    cuts.AddCut("00010113", "00200009a27300008250a04120", "0152103500000000"); //cosPA, 0.995 eta 0.9
  } else if (trainConfig == 62){
    cuts.AddCut("00010113", "0d200009a27300008250a04120", "0152103500000000"); //cosPA, 0.995 eta 0.8
  //----------------------------- configuration for run 2 analysis 13 TeV Triggers --------------------------------------------
  } else if (trainConfig == 70) { // EMC triggers -50, +30 ns
    cuts.AddCut("00010113", "00200009227300008250404000", "0163103100000000","1111100060032220000"); //INT7
  } else if (trainConfig == 71) { // EMC triggers -50, +30 ns
    cuts.AddCut("00085113", "00200009227300008250404000", "0163103100000000","1111100060032220000"); //EG2
  } else if (trainConfig == 72) { // EMC triggers -50, +30 ns
    cuts.AddCut("00083113", "00200009227300008250404000", "0163103100000000","1111100060032220000"); //EG1
  } else if (trainConfig == 73) { // DCAL triggers -50, +30 ns
    cuts.AddCut("00089113", "00200009227300008250404000", "0163103100000000","3885500060032220000"); //DG2
  } else if (trainConfig == 74) { // DCAL triggers -50, +30 ns
    cuts.AddCut("0008b113", "00200009227300008250404000", "0163103100000000","3885500060032220000"); //DG1

   //----------------------------- configuration for run 2 analysis 13 TeVLowB --------------------------------------------
  } else if (trainConfig == 80) { //
    cuts.AddCut("00010113", "00200089227300008280404000", "0152103500000000"); //standard cut

  } else if (trainConfig == 90) { //Standard cut for pp 5 TeV analysis VAND
    cuts.AddCut("00010113", "0d200009227300008250404000", "0152103500000000"); //

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
  TList *ClusterCutList = new TList();

  TList *HeaderList = new TList();
  TObjString *Header2 = new TObjString("BOX");
  HeaderList->Add(Header2);

  EventCutList->SetOwner(kTRUE);
  AliConvEventCuts **analysisEventCuts = new AliConvEventCuts*[numberOfCuts];
  ConvCutList->SetOwner(kTRUE);
  AliConversionPhotonCuts **analysisCuts = new AliConversionPhotonCuts*[numberOfCuts];
  MesonCutList->SetOwner(kTRUE);
  AliConversionMesonCuts **analysisMesonCuts = new AliConversionMesonCuts*[numberOfCuts];
  ClusterCutList->SetOwner(kTRUE);
  AliCaloPhotonCuts **analysisClusterCuts     = new AliCaloPhotonCuts*[numberOfCuts];
  Bool_t enableClustersForTrigger             = kFALSE;

  for(Int_t i = 0; i<numberOfCuts; i++){
    analysisEventCuts[i] = new AliConvEventCuts();

    TString dataInputMultHisto    = "";
    TString mcInputMultHisto      = "";
    TString triggerString         = (cuts.GetEventCut(i)).Data();
    triggerString                 = triggerString(3,2);
    if (triggerString.CompareTo("03")==0)
      triggerString               = "00";
    if (periodNameAnchor.CompareTo("LHC13g") == 0 && triggerString.CompareTo("10")== 0 )
      triggerString               = "00";

    dataInputMultHisto            = Form("%s_%s", periodNameAnchor.Data(), triggerString.Data());
    mcInputMultHisto              = Form("%s_%s", periodNameV0Reader.Data(), triggerString.Data());

    if (doMultiplicityWeighting){
      cout << "enabling mult weighting" << endl;
      analysisEventCuts[i]->SetUseWeightMultiplicityFromFile( kTRUE, fileNameInputForMultWeighing, dataInputMultHisto, mcInputMultHisto );
    }

    analysisEventCuts[i]->SetTriggerMimicking(enableTriggerMimicking);
    analysisEventCuts[i]->SetTriggerOverlapRejecion(enableTriggerOverlapRej);
    analysisEventCuts[i]->SetMaxFacPtHard(maxFacPtHard);
    analysisEventCuts[i]->SetV0ReaderName(V0ReaderName);
    if (periodNameV0Reader.CompareTo("") != 0) analysisEventCuts[i]->SetPeriodEnum(periodNameV0Reader);
    analysisEventCuts[i]->SetLightOutput(runLightOutput);
    analysisEventCuts[i]->InitializeCutsFromCutString((cuts.GetEventCut(i)).Data());
    EventCutList->Add(analysisEventCuts[i]);
    analysisEventCuts[i]->SetFillCutHistograms("",kFALSE);

    if ( (trainConfig >= 70 && trainConfig <= 79) ){
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

        enableClustersForTrigger  = kTRUE;
        analysisClusterCuts[i]    = new AliCaloPhotonCuts();
        analysisClusterCuts[i]->SetV0ReaderName(V0ReaderName);
        analysisClusterCuts[i]->SetCorrectionTaskSetting(corrTaskSetting);
        analysisClusterCuts[i]->SetLightOutput(runLightOutput);
        analysisClusterCuts[i]->InitializeCutsFromCutString((cuts.GetClusterCut(i)).Data());
        ClusterCutList->Add(analysisClusterCuts[i]);
        analysisClusterCuts[i]->SetFillCutHistograms("");
    }

    analysisCuts[i] = new AliConversionPhotonCuts();

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
    if (trainConfig == 21){
      analysisCuts[i]->SetDodEdxSigmaCut(kFALSE);
    }

    ConvCutList->Add(analysisCuts[i]);
    analysisCuts[i]->SetFillCutHistograms("",kFALSE);

    analysisMesonCuts[i] = new AliConversionMesonCuts();
    analysisMesonCuts[i]->SetLightOutput(runLightOutput);
    analysisMesonCuts[i]->InitializeCutsFromCutString((cuts.GetMesonCut(i)).Data());
    MesonCutList->Add(analysisMesonCuts[i]);
    analysisMesonCuts[i]->SetFillCutHistograms("");

    analysisEventCuts[i]->SetAcceptedHeader(HeaderList);
  }

  task->SetEventCutList(numberOfCuts,EventCutList);
  task->SetConversionCutList(numberOfCuts,ConvCutList);
  task->SetMesonCutList(numberOfCuts,MesonCutList);
  task->SetMoveParticleAccordingToVertex(kTRUE);
  task->SetDoMesonAnalysis(kTRUE);
  task->SetDoMesonQA(enableQAMesonTask); //Attention new switch for Pi0 QA
  task->SetDoPhotonQA(enableQAPhotonTask);  //Attention new switch small for Photon QA
  if (enableClustersForTrigger){
    task->SetDoClusterSelectionForTriggerNorm(enableClustersForTrigger);
    task->SetClusterCutList(numberOfCuts,ClusterCutList);
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
