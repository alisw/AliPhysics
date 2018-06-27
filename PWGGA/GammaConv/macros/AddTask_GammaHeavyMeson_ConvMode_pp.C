/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: Friederike Bock, Sandro Wenzel                                 *
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
//($ALIPHYSICS/PWGGA/GammaConv/AliAnalysisTaskHeavyNeutralMesonToGG.cxx) for
//pPb together with all supporting classes
//***************************************************************************************

//***************************************************************************************
//CutHandler contains all cuts for a certain analysis and trainconfig,
//it automatically checks length of cutStrings and takes care of the number of added cuts,
// no specification of the variable 'numberOfCuts' needed anymore.
//***************************************************************************************
class CutHandlerHeavyMesonConv{
  public:
    CutHandlerHeavyMesonConv(Int_t nMax=10){
      nCuts=0; nMaxCuts=nMax; validCuts = true;
      eventCutArray = new TString[nMaxCuts]; photonCutArray = new TString[nMaxCuts]; mesonCutArray = new TString[nMaxCuts]; clusterCutArray = new TString[nMaxCuts];
      for(Int_t i=0; i<nMaxCuts; i++) {eventCutArray[i] = ""; photonCutArray[i] = ""; mesonCutArray[i] = ""; clusterCutArray[i] = "";}
    }

    void AddCut(TString eventCut, TString photonCut, TString mesonCut){
      if(nCuts>=nMaxCuts) {cout << "ERROR in CutHandlerHeavyMesonConv: Exceeded maximum number of cuts!" << endl; validCuts = false; return;}
      if( eventCut.Length()!=8 || photonCut.Length()!=26 || mesonCut.Length()!=16 ) {cout << "ERROR in CutHandlerHeavyMesonConv: Incorrect length of cut string!" << endl; validCuts = false; return;}
      eventCutArray[nCuts]=eventCut; photonCutArray[nCuts]=photonCut; mesonCutArray[nCuts]=mesonCut; clusterCutArray[nCuts]="";
      nCuts++;
      return;
    }
    void AddCut(TString eventCut, TString photonCut, TString mesonCut, TString clusterCut){
      if(nCuts>=nMaxCuts) {cout << "ERROR in CutHandlerHeavyMesonConv: Exceeded maximum number of cuts!" << endl; validCuts = false; return;}
      if( eventCut.Length()!=8 || photonCut.Length()!=26 || mesonCut.Length()!=16 || clusterCut.Length()!=19 ) {cout << "ERROR in CutHandlerHeavyMesonConv: Incorrect length of cut string!" << endl; validCuts = false; return;}
      eventCutArray[nCuts]=eventCut; photonCutArray[nCuts]=photonCut; mesonCutArray[nCuts]=mesonCut; clusterCutArray[nCuts]=clusterCut;
      nCuts++;
      return;
    }

    Bool_t AreValid(){return validCuts;}
    Int_t GetNCuts(){if(validCuts) return nCuts; else return 0;}
    TString GetEventCut(Int_t i){if(validCuts&&i<nMaxCuts&&i>=0) return eventCutArray[i]; else{cout << "ERROR in CutHandlerHeavyMesonConv: GetEventCut wrong index i" << endl;return "";}}
    TString GetPhotonCut(Int_t i){if(validCuts&&i<nMaxCuts&&i>=0) return photonCutArray[i]; else {cout << "ERROR in CutHandlerHeavyMesonConv: GetPhotonCut wrong index i" << endl;return "";}}
    TString GetMesonCut(Int_t i){if(validCuts&&i<nMaxCuts&&i>=0) return mesonCutArray[i]; else {cout << "ERROR in CutHandlerHeavyMesonConv: GetMesonCut wrong index i" << endl;return "";}}
    TString GetClusterCut(Int_t i){if(validCuts&&i<nMaxCuts&&i>=0) return clusterCutArray[i]; else {cout << "ERROR in CutHandlerHeavyMesonConv: GetClusterCut wrong index i" << endl;return "";}}
  private:
    Bool_t validCuts;
    Int_t nCuts; Int_t nMaxCuts;
    TString* eventCutArray;
    TString* photonCutArray;
    TString* mesonCutArray;
    TString* clusterCutArray;
};

//***************************************************************************************
//main function
//***************************************************************************************
void AddTask_GammaHeavyMeson_ConvMode_pp(   Int_t     selectedMeson                 = 0,                    // select the corresponding meson: 0 pi0, 1 eta, 2 eta'
                                            Int_t     trainConfig                   = 1,                    // change different set of cuts
                                            Int_t     isMC                          = 0,                    // run MC
                                            Int_t     enableQAMesonTask             = 0,                    // enable QA in AliAnalysisTaskGammaConvV1
                                            Int_t     enableQAPhotonTask            = 0,                    // enable additional QA task
                                            TString   fileNameInputForWeighting     = "MCSpectraInput.root",// path to file for weigting input / modified acceptance
                                            Int_t     doWeightingPart               = 0,                    // enable Weighting
                                            TString   periodName                    = "DPMJET",             // generator Name
                                            TString   cutnumberAODBranch            = "000000006008400000001500000",  // cutnumber for AOD branch
                                            Int_t     enableExtMatchAndQA           = 0,                    // disabled (0), extMatch (1), extQA_noCellQA (2), extMatch+extQA_noCellQA (3), extQA+cellQA (4), extMatch+extQA+cellQA (5)
                                            Bool_t    isUsingTHnSparse              = kTRUE,                // enable or disable usage of THnSparses for background estimation
                                            Bool_t    enableV0findingEffi           = kFALSE,               // enables V0finding efficiency histograms
                                            Bool_t    enableTriggerMimicking        = kFALSE,               // enable trigger mimicking
                                            Bool_t    enableTriggerOverlapRej       = kFALSE,               // enable trigger overlap rejection
                                            Float_t   maxFacPtHard                  = 3,                    // maximum factor between hardest jet and ptHard generated
                                            TString   periodNameV0Reader            = "",                   // period Name for V0Reader
                                            Bool_t    doMultiplicityWeighting       = kFALSE,               // enable multiplicity weights
                                            TString   fileNameInputForMultWeighing  = "Multiplicity.root",  // file for multiplicity weights
                                            TString   periodNameAnchor              = "",                   // anchor period name for mult weighting
                                            Int_t     runLightOutput                = 0,                    // switch to run light output 0 (disabled), 1 (for CutClasses), 2 (for cutClasses and task)
                                            TString   additionalTrainConfig         = "0"                   // additional counter for trainconfig, this has to be always the last parameter
                                          ) {

  TH1S* histoAcc = 0x0;                     // histo for modified acceptance
  TString corrTaskSetting = ""; // select which correction task setting to use
  //parse additionalTrainConfig flag
  TObjArray *rAddConfigArr = additionalTrainConfig.Tokenize("_");
  if(rAddConfigArr->GetEntries()<1){cout << "ERROR: AddTask_GammaHeavyMeson_ConvMode_pp during parsing of additionalTrainConfig String '" << additionalTrainConfig.Data() << "'" << endl; return;}
  TObjString* rAdditionalTrainConfig = (TObjString*)rAddConfigArr->At(0);
  TString sAdditionalTrainConfig = rAdditionalTrainConfig->GetString();
  if (sAdditionalTrainConfig.Atoi() > 0){
    trainConfig = trainConfig + sAdditionalTrainConfig.Atoi();
    cout << "INFO: AddTask_GammaHeavyMeson_ConvMode_pp running additionalTrainConfig '" << sAdditionalTrainConfig.Atoi() << "', train config: '" << trainConfig << "'" << endl;
  }
  cout << "corrTaskSetting: " << corrTaskSetting.Data() << endl;

  Int_t isHeavyIon    = 2;
  // meson reco mode: 0 - PCM-PCM, 1 - PCM-Calo, 2 - Calo-Calo
  Int_t mesonRecoMode = 0;

  // ================== GetAnalysisManager ===============================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error(Form("AddTask_GammaHeavyMeson_ConvMode_pp_%i",trainConfig), "No analysis manager found.");
    return ;
  }

  // ================== GetInputEventHandler =============================
  AliVEventHandler *inputHandler=mgr->GetInputEventHandler();

  Bool_t isMCForOtherTasks = kFALSE;
  if (isMC > 0) isMCForOtherTasks = kTRUE;


  //========= Add PID Reponse to ANALYSIS manager ====
  if(!(AliPIDResponse*)mgr->GetTask("PIDResponseTask")){
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
    AddTaskPIDResponse(isMCForOtherTasks);
  }

  Printf("here \n");

  //=========  Set Cutnumber for V0Reader ================================
  TString cutnumberPhoton   = "00000008400100001500000000";
  if (  periodNameV0Reader.CompareTo("LHC16f") == 0 || periodNameV0Reader.CompareTo("LHC17g")==0 || periodNameV0Reader.CompareTo("LHC18c")==0 ||
        periodNameV0Reader.CompareTo("LHC17d1") == 0  || periodNameV0Reader.CompareTo("LHC17d12")==0 ||
        periodNameV0Reader.CompareTo("LHC17h3")==0 || periodNameV0Reader.CompareTo("LHC17k1")==0 ||
        periodNameV0Reader.CompareTo("LHC17f8b") == 0 ||
        periodNameV0Reader.CompareTo("LHC16P1JJLowB") == 0 || periodNameV0Reader.CompareTo("LHC16P1Pyt8LowB") == 0 )
    cutnumberPhoton         = "00000088400000000100000000";

  TString cutnumberEvent    = "00000003";
  Bool_t doEtaShift         = kFALSE;
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

    AliConvEventCuts *fEventCuts=NULL;
    if(cutnumberEvent!=""){
      fEventCuts= new AliConvEventCuts(cutnumberEvent.Data(),cutnumberEvent.Data());
      fEventCuts->SetPreSelectionCutFlag(kTRUE);
      fEventCuts->SetV0ReaderName(V0ReaderName);
      if (periodNameV0Reader.CompareTo("") != 0) fEventCuts->SetPeriodEnum(periodNameV0Reader);
      if (runLightOutput > 0) fEventCuts->SetLightOutput(kTRUE);
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
      if (runLightOutput > 0) fCuts->SetLightOutput(kTRUE);
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
  AliAnalysisTaskHeavyNeutralMesonToGG *task=NULL;
  task= new AliAnalysisTaskHeavyNeutralMesonToGG(Form("HeavyNeutralMesonToGG_%i_%i_%i", mesonRecoMode, selectedMeson, trainConfig));
  task->SetIsHeavyIon(isHeavyIon);
  task->SetIsMC(isMC);
  task->SetV0ReaderName(V0ReaderName);
  task->SetCorrectionTaskSetting(corrTaskSetting);

  task->SetMesonRecoMode(mesonRecoMode); // meson reco mode: 0 - PCM-PCM, 1 - PCM-Calo, 2 - Calo-Calo
  if (runLightOutput > 1) task->SetLightOutput(kTRUE);

  //create cut handler
  CutHandlerHeavyMesonConv cuts;

  // *********************************************************************************************************
  // 2.76 TeV  pp Run1
  // *********************************************************************************************************
  if(trainConfig == 1){    // various standard cuts
    cuts.AddCut("00000113", "00200009397300008250400000", "0163103100900000"); // new pi0/eta cut 2.76TeV
    cuts.AddCut("00000113", "00200009397300008250400000", "0163103100000000"); // new pi0/eta cut 2.76TeV without MC smearing
    cuts.AddCut("00000113", "00200009366300003800000000", "0163103100900000"); // standard cut Pi0 pp 2.76TeV PbPb paper 2012
  } else if (trainConfig == 2) { // various standard cuts added signals
    cuts.AddCut("00000123", "00200009397300008250400000", "0163103100900000"); // new pi0/eta cut 2.76TeV
    cuts.AddCut("00000123", "00200009397300008250400000", "0163103100000000"); // new pi0/eta cut 2.76TeV without MC smearing
    cuts.AddCut("00000123", "00200009366300003800000000", "0163103100900000"); // standard cut Pi0 pp 2.76TeV PbPb paper 2012
  } else if (trainConfig == 3) { // additional standards
    cuts.AddCut("00000113", "00200009297002008250400000", "0163103100900000"); // standard cut LHC11h pp 2.76TeV
    cuts.AddCut("00000113", "00200009227302008250404000", "0163101500000000"); // Ana eta analysis prefered
    cuts.AddCut("00000113", "00200008366300000200000000", "0163103100900000"); // standard cut Pi0 pp 7TeV, all photon qualities
  } else if (trainConfig == 4) { // additional standards added signals
    cuts.AddCut("00000123", "00200009297002008250400000", "0163103100900000"); // standard cut LHC11h pp 2.76TeV
    cuts.AddCut("00000123", "00200009227302008250404000", "0163101500000000"); // Ana eta analysis prefered
    cuts.AddCut("00000123", "00200008366300000200000000", "0163103100900000"); // standard cut Pi0 pp 7TeV, all photon qualities

  // *********************************************************************************************************
  // 8 TeV  pp Run1
  // *********************************************************************************************************
  } else if (trainConfig == 100) {
    cuts.AddCut("00010113", "00200009227300008250404000", "0152103500000000"); // New standard cut for 8TeV analysis V0AND with double counting cut, TOF removed
    cuts.AddCut("00010013", "00200009227300008250404000", "0152103500000000"); // no SPD pileup cut

  // *********************************************************************************************************
  // 7 TeV  pp Run1
  // *********************************************************************************************************
  } else if (trainConfig == 200) {
    cuts.AddCut("00000113", "00200009227300008250404000", "0152103500000000"); //New standard cut for 7TeV analysis V0OR with double counting cut, TOF removed
    cuts.AddCut("00000013", "00200009227300008250404000", "0152103500000000"); // no SPD pileup cut
    cuts.AddCut("00000113", "00200009227302008250400000", "0152103500000000"); //New standard cut for 7TeV analysis V0OR
    cuts.AddCut("00000113", "00200008366300000200000000", "0163103100900000"); //Old standard cut for 7TeV analysis V0OR
    cuts.AddCut("00000113", "00200009360300007800004000", "0263103100900000"); //dalitz: New Standard Only MB, standard pp7Tev cut dalitz

  // *********************************************************************************************************
  // 5 TeV  pp Run2
  // *********************************************************************************************************
  } else if (trainConfig == 300) {
    cuts.AddCut("00010113", "00200009227302008250400000", "0163103100000000"); //
  } else if (trainConfig == 301) {
    cuts.AddCut("00010113", "00200009227300008250404000", "0163103100000000","1111100017032220000"); //INT7
    cuts.AddCut("00052113", "00200009227300008250404000", "0163103100000000","1111100017032220000"); //EMC7
    cuts.AddCut("00085113", "00200009227300008250404000", "0163103100000000","1111100017032220000"); //EG2
    cuts.AddCut("00083113", "00200009227300008250404000", "0163103100000000","1111100017032220000"); //EG1

  // *********************************************************************************************************
  // 13 TeV  pp Run2
  // *********************************************************************************************************
  } else if (trainConfig == 400) {
    cuts.AddCut("00010113", "00200009227300008250404000", "0163103100000000"); //INT7
  } else if (trainConfig == 401) {
    cuts.AddCut("00010113", "00200009227300008250404000", "0163103100000000","1111100017032220000"); //INT7
    cuts.AddCut("00052113", "00200009227300008250404000", "0163103100000000","1111100017032220000"); //EMC7
    cuts.AddCut("00085113", "00200009227300008250404000", "0163103100000000","1111100017032220000"); //EG2
    cuts.AddCut("00083113", "00200009227300008250404000", "0163103100000000","1111100017032220000"); //EG1
  } else if (trainConfig == 402) {
    cuts.AddCut("00010113", "00200009227300008250404000", "0163103100000000","1111100067032220000"); //INT7
    cuts.AddCut("00052113", "00200009227300008250404000", "0163103100000000","1111100067032220000"); //EMC7
    cuts.AddCut("00085113", "00200009227300008250404000", "0163103100000000","1111100067032220000"); //EG2
    cuts.AddCut("00083113", "00200009227300008250404000", "0163103100000000","1111100067032220000"); //EG1

  } else {
    Error(Form("HeavyNeutralMesonToGG_%i",trainConfig), "wrong trainConfig variable no cuts have been specified for the configuration");
    return;
  }

  if(!cuts.AreValid()){
    cout << "\n\n****************************************************" << endl;
    cout << "ERROR: No valid cuts stored in CutHandlerHeavyMesonConv! Returning..." << endl;
    cout << "****************************************************\n\n" << endl;
    return;
  }

  TList *EventCutList   = new TList();
  TList *ConvCutList    = new TList();
  TList *MesonCutList   = new TList();
  TList *ClusterCutList = new TList();

  Int_t numberOfCuts = cuts.GetNCuts();

  TList *HeaderList = new TList();
  if (periodName.Contains("LHC12i3")){
    TObjString *Header2 = new TObjString("BOX");
    HeaderList->Add(Header2);
  } else if (periodName.CompareTo("LHC14e2b")==0){
    TObjString *Header2 = new TObjString("pi0_1");
    HeaderList->Add(Header2);
    TObjString *Header3 = new TObjString("eta_2");
    HeaderList->Add(Header3);
  }

  TString energy = "";
  TString mcName = "";
  TString mcNameAdd = "";
  if (periodName.Contains("WOSDD")){
    mcNameAdd = "_WOSDD";
  } else if (periodName.Contains("WSDD")){
    mcNameAdd = "_WSDD";
  }
  if (periodName.Contains("LHC12i3")){
    energy = "2760GeV";
    mcName = "Pythia8_LHC12i3";
  } else if (periodName.Contains("LHC12f1a")){
    energy = "2760GeV";
    mcName = "Pythia8_LHC12f1a";
  } else if (periodName.Contains("LHC12f1b")){
    energy = "2760GeV";
    mcName = "Phojet_LHC12f1b";
  } else if (periodName.Contains("LHC14e2a")){
    energy = "8TeV";
    mcName = "Pythia8_LHC14e2a";
  } else if (periodName.Contains("LHC14e2b")){
    energy = "8TeV";
    mcName = "Pythia8_LHC14e2b";
  } else if (periodName.Contains("LHC14e2c")){
    energy = "8TeV";
    mcName = "Phojet_LHC14e2c";
  }

  EventCutList->SetOwner(kTRUE);
  AliConvEventCuts **analysisEventCuts        = new AliConvEventCuts*[numberOfCuts];
  ConvCutList->SetOwner(kTRUE);
  AliConversionPhotonCuts **analysisCuts      = new AliConversionPhotonCuts*[numberOfCuts];
  MesonCutList->SetOwner(kTRUE);
  AliConversionMesonCuts **analysisMesonCuts  = new AliConversionMesonCuts*[numberOfCuts];
  ClusterCutList->SetOwner(kTRUE);
  AliCaloPhotonCuts **analysisClusterCuts     = new AliCaloPhotonCuts*[numberOfCuts];
  Bool_t enableClustersForTrigger             = kFALSE;

  for(Int_t i = 0; i<numberOfCuts; i++){
    analysisEventCuts[i] = new AliConvEventCuts();
    TString fitNamePi0            = Form("Pi0_Fit_Data_%s",energy.Data());
    TString fitNameEta            = Form("Eta_Fit_Data_%s",energy.Data());
    TString fAddedSignalString    = (cuts.GetEventCut(i)).Data();
    fAddedSignalString            = fAddedSignalString(6,1);
    Bool_t fAddedSignal           = kFALSE;
    if (fAddedSignalString.CompareTo("2") == 0)
      fAddedSignal                = kTRUE;

    TString mcInputNamePi0        = "";
    TString mcInputNameEta        = "";
    if (fAddedSignal && (periodName.Contains("LHC12i3") || periodName.CompareTo("LHC14e2b")==0)){
      mcInputNamePi0              = Form("Pi0_%s%s_addSig_%s", mcName.Data(), mcNameAdd.Data(), energy.Data() );
      mcInputNameEta              = Form("Eta_%s%s_addSig_%s", mcName.Data(), mcNameAdd.Data(), energy.Data() );
    } else {
      mcInputNamePi0              = Form("Pi0_%s%s_%s", mcName.Data(), mcNameAdd.Data(), energy.Data() );
      mcInputNameEta              = Form("Eta_%s%s_%s", mcName.Data(), mcNameAdd.Data(), energy.Data() );
    }
    if (doWeightingPart) analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE, fileNameInputForWeighting, mcInputNamePi0, mcInputNameEta, "",fitNamePi0,fitNameEta);


    TString dataInputMultHisto  = "";
    TString mcInputMultHisto    = "";
    TString triggerString       = cuts.GetEventCut(i);
    triggerString               = triggerString(3,2);

    dataInputMultHisto          = Form("%s_%s", periodNameAnchor.Data(), triggerString.Data());
    mcInputMultHisto            = Form("%s_%s", periodNameV0Reader.Data(), triggerString.Data());

    if (doMultiplicityWeighting){
      cout << "enabling mult weighting" << endl;
      analysisEventCuts[i]->SetUseWeightMultiplicityFromFile( kTRUE, fileNameInputForMultWeighing, dataInputMultHisto, mcInputMultHisto );
    }

    analysisEventCuts[i]->SetTriggerMimicking(enableTriggerMimicking);
    analysisEventCuts[i]->SetTriggerOverlapRejecion(enableTriggerOverlapRej);
    analysisEventCuts[i]->SetMaxFacPtHard(maxFacPtHard);
    analysisEventCuts[i]->SetV0ReaderName(V0ReaderName);
    analysisEventCuts[i]->SetCorrectionTaskSetting(corrTaskSetting);
    if (periodNameV0Reader.CompareTo("") != 0) analysisEventCuts[i]->SetPeriodEnum(periodNameV0Reader);
    if (runLightOutput > 0) analysisEventCuts[i]->SetLightOutput(kTRUE);
    analysisEventCuts[i]->InitializeCutsFromCutString((cuts.GetEventCut(i)).Data());
    EventCutList->Add(analysisEventCuts[i]);


    analysisEventCuts[i]->SetFillCutHistograms("",kFALSE);

    if ( trainConfig == 301 || trainConfig == 401 || trainConfig == 402   ){
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
    analysisCuts[i]->SetV0ReaderName(V0ReaderName);
    if (runLightOutput > 0) analysisCuts[i]->SetLightOutput(kTRUE);
    analysisCuts[i]->InitializeCutsFromCutString((cuts.GetPhotonCut(i)).Data());
    analysisCuts[i]->SetIsHeavyIon(isHeavyIon);
    ConvCutList->Add(analysisCuts[i]);
    analysisCuts[i]->SetFillCutHistograms("",kFALSE);

    analysisMesonCuts[i] = new AliConversionMesonCuts();
    if (runLightOutput > 0) analysisMesonCuts[i]->SetLightOutput(kTRUE);
    analysisMesonCuts[i]->SetRunningMode(0);
    analysisMesonCuts[i]->InitializeCutsFromCutString((cuts.GetMesonCut(i)).Data());
    MesonCutList->Add(analysisMesonCuts[i]);
    analysisMesonCuts[i]->SetFillCutHistograms("");
    analysisEventCuts[i]->SetAcceptedHeader(HeaderList);
  }

  task->SetEventCutList(numberOfCuts,EventCutList);
  task->SetConversionCutList(numberOfCuts,ConvCutList);
  task->SetMesonCutList(numberOfCuts,MesonCutList);
  task->SetMoveParticleAccordingToVertex(kTRUE);
  task->SetCorrectionTaskSetting(corrTaskSetting);
  if (enableClustersForTrigger){
    task->SetDoClusterSelectionForTriggerNorm(enableClustersForTrigger);
    task->SetCaloCutList(numberOfCuts,ClusterCutList);
    task->SetDoClusterQA(1);  //Attention new switch small for Cluster QA
    if(enableExtMatchAndQA > 1){ task->SetPlotHistsExtQA(kTRUE);}
  }
  task->SetMesonType(selectedMeson);
  task->SetDoMesonQA(enableQAMesonTask); //Attention new switch for Pi0 QA
  task->SetDoPhotonQA(enableQAPhotonTask);  //Attention new switch small for Photon QA
  task->SetUseTHnSparse(isUsingTHnSparse);

  //connect containers
  AliAnalysisDataContainer *coutput =
  mgr->CreateContainer(!(corrTaskSetting.CompareTo("")) ? Form("HeavyNeutralMesonToGG_%i_%i_%i", mesonRecoMode, selectedMeson, trainConfig)
                                                        :  Form("HeavyNeutralMesonToGG_%i_%i_%i_%s", mesonRecoMode, selectedMeson, trainConfig, corrTaskSetting.Data()), TList::Class(),
                        AliAnalysisManager::kOutputContainer,Form("HeavyNeutralMesonToGG_%i_%i.root",mesonRecoMode,trainConfig) );

  mgr->AddTask(task);
  mgr->ConnectInput(task,0,cinput);
  mgr->ConnectOutput(task,1,coutput);

  return;

}
