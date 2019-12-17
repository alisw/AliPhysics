class CutHandlerConvDalitz{
  public:
    CutHandlerConvDalitz(Int_t nMax=10){
      nCuts=0; nMaxCuts=nMax; validCuts = true;
      eventCutArray = new TString[nMaxCuts]; photonCutArray = new TString[nMaxCuts]; mesonCutArray = new TString[nMaxCuts]; elecCutArray = new TString[nMaxCuts];
      for(Int_t i=0; i<nMaxCuts; i++) {eventCutArray[i] = ""; photonCutArray[i] = ""; mesonCutArray[i] = ""; elecCutArray[i] = "";}
    }

    void AddCut(TString eventCut, TString photonCut, TString elecCut, TString mesonCut){
      if(nCuts>=nMaxCuts) {cout << "ERROR in CutHandlerConvDalitz: Exceeded maximum number of cuts!" << endl; validCuts = false; return;}
      if( eventCut.Length()!=8 || photonCut.Length()!=26 || mesonCut.Length()!=16 || elecCut.Length()!=20 ) {cout << "ERROR in CutHandlerConvDalitz: Incorrect length of cut string!" << endl; validCuts = false; return;}
      eventCutArray[nCuts]=eventCut; photonCutArray[nCuts]=photonCut; mesonCutArray[nCuts]=mesonCut; elecCutArray[nCuts]=elecCut;
      nCuts++;
      return;
    }

    Bool_t AreValid(){return validCuts;}
    Int_t GetNCuts(){if(validCuts) return nCuts; else return 0;}
    TString GetEventCut(Int_t i){if(validCuts&&i<nMaxCuts&&i>=0) return eventCutArray[i]; else{cout << "ERROR in CutHandlerConvDalitz: GetEventCut wrong index i" << endl;return "";}}
    TString GetPhotonCut(Int_t i){if(validCuts&&i<nMaxCuts&&i>=0) return photonCutArray[i]; else {cout << "ERROR in CutHandlerConvDalitz: GetPhotonCut wrong index i" << endl;return "";}}
    TString GetMesonCut(Int_t i){if(validCuts&&i<nMaxCuts&&i>=0) return mesonCutArray[i]; else {cout << "ERROR in CutHandlerConvDalitz: GetMesonCut wrong index i" << endl;return "";}}
    TString GetElecCut(Int_t i){if(validCuts&&i<nMaxCuts&&i>=0) return elecCutArray[i]; else {cout << "ERROR in CutHandlerConvDalitz: GetElecCut wrong index i" << endl;return "";}}
  private:
    Bool_t validCuts;
    Int_t nCuts; Int_t nMaxCuts;
    TString* eventCutArray;
    TString* photonCutArray;
    TString* mesonCutArray;
    TString* elecCutArray;
};

void AddTask_GammaConvDalitzV1_PbPb(  Int_t   trainConfig               = 1,
                                      Bool_t  isMC                      = kFALSE, //run MC
                                      Bool_t  enableQAMesonTask         = kFALSE, //enable QA in AliAnalysisTaskGammaConvDalitzV1
                                      Bool_t  enableDoMesonChic         = kFALSE, // enable additional Chic analysis
                                      TString fileNameInputForWeighting = "MCSpectraInput.root", // path to file for weigting input
                                      Bool_t  doWeighting               = kFALSE,  //enable Weighting
                                      Bool_t  enableUseTHnSparse        = kTRUE,                          // enable THnSparse for mixed event BG
                                      TString cutnumberAODBranch        = "0000000060084001001500000",
                                      TString   periodName                      = "LHC13d2",
                                      Int_t     enableV0EffiStudies             = 0,
                                      Int_t     enableTriggerMimicking  = 0,
                                      Bool_t    enableTriggerOverlapRej = kFALSE,    // enable trigger overlap rejection
                                      TString   settingMaxFacPtHard     = "3.",       // maximum factor between hardest jet and ptHard generated
                                      TString   periodNameV0Reader      = "",
                                      Bool_t    runLightOutput          = kFALSE,
                                      Int_t     enableMatBudWeightsPi0          = 0,
                                      TString   additionalTrainConfig   = "0"
                                 )  {

  Int_t isHeavyIon = 1;

  if (additionalTrainConfig.Atoi() > 0){
    trainConfig = trainConfig + additionalTrainConfig.Atoi();
  }

  // ================== GetAnalysisManager ===============================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error(Form("AddTask_GammaConvDalitzV1_PbPb_%i",trainConfig), "No analysis manager found.");
    return ;
  }

  TObjArray *rmaxFacPtHardSetting = settingMaxFacPtHard.Tokenize("_");
  if(rmaxFacPtHardSetting->GetEntries()<1){cout << "ERROR: AddTask_GammaConvDalitzV1_PbPb during parsing of settingMaxFacPtHard String '" << settingMaxFacPtHard.Data() << "'" << endl; return;}
  Bool_t fMinPtHardSet        = kFALSE;
  Double_t minFacPtHard       = -1;
  Bool_t fMaxPtHardSet        = kFALSE;
  Double_t maxFacPtHard       = 100;
  Bool_t fSingleMaxPtHardSet  = kFALSE;
  Double_t maxFacPtHardSingle = 100;
  for(Int_t i = 0; i<rmaxFacPtHardSetting->GetEntries() ; i++){
    TObjString* tempObjStrPtHardSetting     = (TObjString*) rmaxFacPtHardSetting->At(i);
    TString strTempSetting                  = tempObjStrPtHardSetting->GetString();
    if(strTempSetting.BeginsWith("MINPTHFAC:")){
      strTempSetting.Replace(0,10,"");
      minFacPtHard               = strTempSetting.Atof();
      cout << "running with min pT hard jet fraction of: " << minFacPtHard << endl;
      fMinPtHardSet        = kTRUE;
    } else if(strTempSetting.BeginsWith("MAXPTHFAC:")){
      strTempSetting.Replace(0,10,"");
      maxFacPtHard               = strTempSetting.Atof();
      cout << "running with max pT hard jet fraction of: " << maxFacPtHard << endl;
      fMaxPtHardSet        = kTRUE;
    } else if(strTempSetting.BeginsWith("MAXPTHFACSINGLE:")){
      strTempSetting.Replace(0,16,"");
      maxFacPtHardSingle         = strTempSetting.Atof();
      cout << "running with max single particle pT hard fraction of: " << maxFacPtHardSingle << endl;
      fSingleMaxPtHardSet        = kTRUE;
    } else if(rmaxFacPtHardSetting->GetEntries()==1 && strTempSetting.Atof()>0){
      maxFacPtHard               = strTempSetting.Atof();
      cout << "running with max pT hard jet fraction of: " << maxFacPtHard << endl;
      fMaxPtHardSet        = kTRUE;
    }
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
  TString cutnumberPhoton = "00000070000000000500004000";
  if (periodNameV0Reader.CompareTo("LHC17n") == 0 || periodNameV0Reader.Contains("LHC17j7"))
    cutnumberPhoton         = "00000088400100001500000000";
  TString cutnumberEvent  = "10000003";
  TString ElecCuts        = "30105400000003300000";            //Electron Cuts

  Bool_t enableV0findingEffi = kFALSE;
  if(enableV0EffiStudies > 0){
    enableV0findingEffi = kTRUE;
    cutnumberPhoton = "00000070000000000500004000";
    if(enableV0EffiStudies == 1){
      cutnumberEvent = "50100013";
    }else if(enableV0EffiStudies == 2){
      cutnumberEvent = "52500013";
    }
  }
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();

  //========= Add V0 Reader to  ANALYSIS manager if not yet existent =====
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
      if (runLightOutput > 0) fEventCuts->SetLightOutput(kTRUE);
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
      fCuts->SetProcessAODCheck(kFALSE); // if processAODcheckForV0s is kTRUE, also check for V0s to be contained in AliAODs and AliAODGammaConversion.root
      if(trainConfig == 182 || trainConfig == 183 || trainConfig == 184 || trainConfig == 185){
        fCuts->SetDodEdxSigmaCut(kFALSE);
      }
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
    AliLog::SetGlobalLogLevel(AliLog::kInfo);

    //connect input V0Reader
    mgr->AddTask(fV0ReaderV1);
    mgr->ConnectInput(fV0ReaderV1,0,cinput);

  }
  //================================================
  //========= Add Electron Selector ================

  if( !(AliDalitzElectronSelector*)mgr->GetTask("ElectronSelector") ){

    AliDalitzElectronSelector *fElectronSelector = new AliDalitzElectronSelector("ElectronSelector");
    // Set AnalysisCut Number
    AliDalitzElectronCuts *fElecCuts=0;
    if( ElecCuts!=""){
      fElecCuts= new AliDalitzElectronCuts(ElecCuts.Data(),ElecCuts.Data());
      fElecCuts->SetUseCrossedRows(kTRUE);
      if(fElecCuts->InitializeCutsFromCutString(ElecCuts.Data())){
        fElectronSelector->SetDalitzElectronCuts(fElecCuts);
        fElecCuts->SetFillCutHistograms("",kTRUE);
      }
    }
    fElectronSelector->Init();
    mgr->AddTask(fElectronSelector);
    AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer();

    //connect input V0Reader
    mgr->ConnectInput (fElectronSelector,0,cinput1);
  }

  //================================================
  //========= Add task to the ANALYSIS manager =====
  //================================================
  //            find input container
   AliAnalysisTaskGammaConvDalitzV1 *task=NULL;

  task= new AliAnalysisTaskGammaConvDalitzV1(Form("GammaConvDalitzV1_%i",trainConfig));
  task->SetIsHeavyIon(isHeavyIon);
  task->SetIsMC(isMC);
  task->SetV0ReaderName(V0ReaderName);

  CutHandlerConvDalitz cuts; // object to add the cut

  if ( trainConfig == 1 ) { // LHC15o, kINT7, cent. from V0M, reject added particles
    cuts.AddCut("10110013", "00200009247602008250404000", "20405400233202223710", "0263103500900000"); //  0-10%
    cuts.AddCut("11210013", "00200009247602008250404000", "20405400233202223710", "0263103500900000"); // 10-20%
    cuts.AddCut("12510013", "00200009247602008250404000", "20405400233202223710", "0263103500900000"); // 20-50%
    cuts.AddCut("15910013", "00200009247602008250404000", "20405400233202223710", "0263103500900000"); // 50-90%
    cuts.AddCut("10910013", "00200009247602008250404000", "20405400233202223710", "0263103500900000"); //  0-90%
  } else if ( trainConfig == 2 ) { //XeXe configurations
    cuts.AddCut("10210113", "00200009327000008250400000", "10885400233102227610", "0263103500900000"); // 0-20 // NOTE Check electron
    cuts.AddCut("12410113", "00200009327000008250400000", "10885400233102227610", "0263103500900000"); // 20-40
    cuts.AddCut("10410113", "00200009327000008250400000", "10885400233102227610", "0263103500900000"); // 0-40
    cuts.AddCut("14810113", "00200009327000008250400000", "10885400233102227610", "0263103500900000"); // 40-80
  } else if ( trainConfig == 3 ) { //XeXe configurations, only 40-80 for testing purposes
    cuts.AddCut("14810113", "00200009327000008250400000", "10885400233102227610", "0263103500900000"); // 0-20 // NOTE Check electron
  }

  if(!cuts.AreValid()){
    cout << "\n\n****************************************************" << endl;
    cout << "ERROR: No valid cuts stored in CutHandlerConv! Returning..." << endl;
    cout << "****************************************************\n\n" << endl;
    return;
  }

  Int_t numberOfCuts = cuts.GetNCuts();



  TList *EventCutList = new TList();
  TList *ConvCutList  = new TList();
  TList *MesonCutList = new TList();
  TList *ElecCutList  = new TList();


  TList *HeaderList = new TList();
  if (periodName.CompareTo("LHC13d2")==0){
    TObjString *Header1 = new TObjString("pi0_1");
    HeaderList->Add(Header1);

  } else if (periodName.CompareTo("LHC12a17x_fix")==0){
    TObjString *Header1 = new TObjString("PARAM");
    HeaderList->Add(Header1);
  } else if (periodName.CompareTo("LHC14a1a")==0){
    if (headerSelectionInt == 1){
      TObjString *Header1 = new TObjString("pi0_1");
      HeaderList->Add(Header1);
    } else if (headerSelectionInt == 2){
      TObjString *Header1 = new TObjString("eta_2");
      HeaderList->Add(Header1);
    } else if (headerSelectionInt == 3){
      TString nameHeaders[2]    = { "pi0_1", "eta_2" };
      for (Int_t iHead = 0; iHead < 2; iHead++ ){
        TObjString *Header = new TObjString(nameHeaders[iHead]);
        HeaderList->Add(Header);
      }
    } else if (headerSelectionInt == 4){
      TObjString *Header1 = new TObjString("pi0EMC_3");
      HeaderList->Add(Header1);
    } else if (headerSelectionInt == 5){
      TObjString *Header1 = new TObjString("etaEMC_5");
      HeaderList->Add(Header1);
    } else if (headerSelectionInt == 6){
      TString nameHeaders[2]    = { "pi0EMC_3", "etaEMC_5" };
      for (Int_t iHead = 0; iHead < 2; iHead++ ){
        TObjString *Header = new TObjString(nameHeaders[iHead]);
        HeaderList->Add(Header);
      }
    } else if (headerSelectionInt == 7){
      TString nameHeaders[4]    = { "pi0_1", "eta_2", "pi0EMC_3", "etaEMC_5" };
      for (Int_t iHead = 0; iHead < 4; iHead++ ){
        TObjString *Header = new TObjString(nameHeaders[iHead]);
        HeaderList->Add(Header);
      }
    } else if (headerSelectionInt == 8){
      TObjString *Header1 = new TObjString("gEMCPhoton_7");
      HeaderList->Add(Header1);
    } else if (headerSelectionInt == 9){
      TString nameHeaders[10]   = { "Pythia_Jets_PtHard_1_10", "Pythia_Jets_PtHard_2_10", "Pythia_Jets_PtHard_3_10", "Pythia_Jets_PtHard_4_10", "Pythia_Jets_PtHard_5_10",
        "Pythia_Jets_PtHard_6_10", "Pythia_Jets_PtHard_7_10", "Pythia_Jets_PtHard_8_10", "Pythia_Jets_PtHard_9_10", "Pythia_Jets_PtHard_10_10"
      };
      for (Int_t iHead = 0; iHead < 10; iHead++ ){
        TObjString *Header = new TObjString(nameHeaders[iHead]);
        HeaderList->Add(Header);
      }
    } else if (headerSelectionInt == 10){
      TObjString *Header1 = new TObjString("pythia_bele_10_10");
      HeaderList->Add(Header1);
    } else if (headerSelectionInt == 11){
      TString nameHeaders[34]   = { "Pythia_Jets_PtHard_1_10", "Pythia_Jets_PtHard_2_10", "Pythia_Jets_PtHard_3_10", "Pythia_Jets_PtHard_4_10", "Pythia_Jets_PtHard_5_10",
        "Pythia_Jets_PtHard_6_10", "Pythia_Jets_PtHard_7_10", "Pythia_Jets_PtHard_8_10", "Pythia_Jets_PtHard_9_10", "Pythia_Jets_PtHard_10_10",
        "gEMCPhoton_7", "flat pt kstar_8", "flat pt kstarbar_9", "pythia_cele_10_10", "pythia_cele_18_10",
        "pythia_cele_30_10", "pythia_cele_50_10", "pythia_ccbar_10_10", "pythia_ccbar_18_10", "pythia_ccbar_30_10",
        "pythia_ccbar_50_10", "pythia_bele_10_10", "pythia_bele_18_10", "pythia_bele_30_10", "pythia_bele_50_10",
        "pythia_bbbar_10_10", "pythia_bbbar_18_10", "pythia_bbbar_30_10", "pythia_bbbar_50_10", "pi0_1",
        "eta_2", "pi0EMC_3", "etaEMC_5", "hijing_0"
      };
      for (Int_t iHead = 0; iHead < 34; iHead++ ){
        TObjString *Header = new TObjString(nameHeaders[iHead]);
        HeaderList->Add(Header);
      }
    } else if (headerSelectionInt == 12){
      TObjString *Header1 = new TObjString("pi0PHS_4");
      HeaderList->Add(Header1);
    } else if (headerSelectionInt == 13){
      TObjString *Header1 = new TObjString("etaPHS_6");
      HeaderList->Add(Header1);
    } else if (headerSelectionInt == 14){
      TString nameHeaders[2]    = { "pi0PHS_4", "etaPHS_6" };
      for (Int_t iHead = 0; iHead < 2; iHead++ ){
        TObjString *Header = new TObjString(nameHeaders[iHead]);
        HeaderList->Add(Header);
      }
    } else {
      TString nameHeaders[2]    = { "pi0_1", "eta_2" };
      for (Int_t iHead = 0; iHead < 2; iHead++ ){
        TObjString *Header = new TObjString(nameHeaders[iHead]);
        HeaderList->Add(Header);
      }
    }
  } else if (periodName.CompareTo("LHC14a1b")==0 || periodName.CompareTo("LHC14a1c")==0){
    if (headerSelectionInt == 1 || headerSelectionInt == 2 || headerSelectionInt == 3 ){
      TObjString *Header1 = new TObjString("BOX");
      HeaderList->Add(Header1);
    } if (headerSelectionInt == 4 || headerSelectionInt == 5 || headerSelectionInt == 6 ){
      TObjString *Header1 = new TObjString("PARAM_EMC");
      HeaderList->Add(Header1);
    } if (headerSelectionInt == 12 || headerSelectionInt == 13 || headerSelectionInt == 14 ){
      TObjString *Header1 = new TObjString("PARAM_PHOS");
      HeaderList->Add(Header1);
    }
  } else if (periodName.CompareTo("LHC16h4")==0){
    if (headerSelectionInt == 1){
      TObjString *Header1 = new TObjString("Injector (pi0)_1");
      HeaderList->Add(Header1);
    } else if (headerSelectionInt == 2){
      TObjString *Header1 = new TObjString("Injector (eta)_2");
      HeaderList->Add(Header1);
    } else {
      TObjString *Header1 = new TObjString("Injector (pi0)_1");
      HeaderList->Add(Header1);
      TObjString *Header2 = new TObjString("Injector (eta)_2");
      HeaderList->Add(Header2);
    }
  }

  EventCutList->SetOwner(kTRUE);
  AliConvEventCuts **analysisEventCuts          = new AliConvEventCuts*[numberOfCuts];

  ConvCutList->SetOwner(kTRUE);
  AliConversionPhotonCuts **analysisCuts       = new AliConversionPhotonCuts*[numberOfCuts];

  MesonCutList->SetOwner(kTRUE);
  AliConversionMesonCuts **analysisMesonCuts   = new AliConversionMesonCuts*[numberOfCuts];

  ElecCutList->SetOwner(kTRUE);
  AliDalitzElectronCuts **analysisElecCuts     = new AliDalitzElectronCuts*[numberOfCuts];

  for(Int_t i = 0; i<numberOfCuts; i++){


    analysisEventCuts[i] = new AliConvEventCuts();

    analysisEventCuts[i]->SetTriggerMimicking(enableTriggerMimicking);
    analysisEventCuts[i]->SetTriggerOverlapRejecion(enableTriggerOverlapRej);
    if(fMinPtHardSet)
      analysisEventCuts[i]->SetMinFacPtHard(minFacPtHard);
    if(fMaxPtHardSet)
      analysisEventCuts[i]->SetMaxFacPtHard(maxFacPtHard);
    if(fSingleMaxPtHardSet)
      analysisEventCuts[i]->SetMaxFacPtHardSingleParticle(maxFacPtHardSingle);
    analysisEventCuts[i]->SetV0ReaderName(V0ReaderName);
    if (periodNameV0Reader.CompareTo("") != 0) analysisEventCuts[i]->SetPeriodEnum(periodNameV0Reader);
    analysisEventCuts[i]->SetLightOutput(runLightOutput);

    if( ! analysisEventCuts[i]->InitializeCutsFromCutString((cuts.GetEventCut(i)).Data()) ){
      cout<<"ERROR:  analysisEventCuts ["<<i<<"]"<<endl;
      return 0;
    }
    if (periodName.CompareTo("LHC14a1b") ==0 || periodName.CompareTo("LHC14a1c") ==0 ){
      if (headerSelectionInt == 1 || headerSelectionInt == 4 || headerSelectionInt == 12 ) analysisEventCuts[i]->SetAddedSignalPDGCode(111);
      if (headerSelectionInt == 2 || headerSelectionInt == 5 || headerSelectionInt == 13 ) analysisEventCuts[i]->SetAddedSignalPDGCode(221);
    }
    EventCutList->Add(analysisEventCuts[i]);
    EventCutList->Add(analysisEventCuts[i]);
    analysisEventCuts[i]->SetFillCutHistograms("",kFALSE);

    analysisCuts[i] = new AliConversionPhotonCuts();
    if (enableMatBudWeightsPi0 > 0){
      if (isMC > 0){
        if (analysisCuts[i]->InitializeMaterialBudgetWeights(enableMatBudWeightsPi0,filenameMatBudWeights)){
          initializedMatBudWeigths_existing = kTRUE;}
        else {cout << "ERROR The initialization of the materialBudgetWeights did not work out." << endl;}
      }
      else {cout << "ERROR 'enableMatBudWeightsPi0'-flag was set > 0 even though this is not a MC task. It was automatically reset to 0." << endl;}
    }
    analysisCuts[i]->SetIsHeavyIon(isHeavyIon);
    analysisCuts[i]->SetV0ReaderName(V0ReaderName);
    analysisCuts[i]->SetLightOutput(runLightOutput);
    if( ! analysisCuts[i]->InitializeCutsFromCutString((cuts.GetPhotonCut(i)).Data() ) ){
      cout<<"ERROR: analysisCuts [" <<i<<"]"<<endl;
      return 0;
    }

    ConvCutList->Add(analysisCuts[i]);
    analysisCuts[i]->SetFillCutHistograms("",kFALSE);

    analysisMesonCuts[i] = new AliConversionMesonCuts();
    analysisMesonCuts[i]->SetLightOutput(runLightOutput);
    if ( ! analysisMesonCuts[i]->InitializeCutsFromCutString((cuts.GetMesonCut(i)).Data()) ) {
      cout<<"ERROR: analysisMesonCuts [ " <<i<<" ] "<<endl;
      return 0;
    }
    MesonCutList->Add(analysisMesonCuts[i]);
    analysisMesonCuts[i]->SetFillCutHistograms("");

    analysisEventCuts[i]->SetAcceptedHeader(HeaderList);
    analysisElecCuts[i] = new AliDalitzElectronCuts();
    analysisElecCuts[i]->SetUseCrossedRows(kTRUE);

    if( !analysisElecCuts[i]->InitializeCutsFromCutString((cuts.GetElecCut(i)).Data()) ){
      cout<< "ERROR:  analysisElecCuts [ " <<i<<" ] "<<endl;
      return 0;
    }
    ElecCutList->Add(analysisElecCuts[i]);
    analysisElecCuts[i]->SetFillCutHistograms("",kFALSE,(cuts.GetElecCut(i)).Data());
  }

  task->SetEventCutList(numberOfCuts,EventCutList);
  task->SetConversionCutList(numberOfCuts,ConvCutList);
  task->SetMesonCutList(MesonCutList);
  task->SetElectronCutList(ElecCutList);
  task->SetMoveParticleAccordingToVertex(kTRUE);
  if(enableQAMesonTask) task->SetDoMesonQA(kTRUE);
  if(enableDoMesonChic) task->SetDoChicAnalysis(kTRUE);
  task->SetDoTHnSparse(enableUseTHnSparse);

  //connect containers
  AliAnalysisDataContainer *coutput =
  mgr->CreateContainer( Form("GammaConvDalitzV1_%i",trainConfig), TList::Class(),
                        AliAnalysisManager::kOutputContainer,Form("GammaConvDalitzV1_%i.root",trainConfig));

  mgr->AddTask(task);
  mgr->ConnectInput(task,0,cinput);
  mgr->ConnectOutput(task,1,coutput);

  return;
}
