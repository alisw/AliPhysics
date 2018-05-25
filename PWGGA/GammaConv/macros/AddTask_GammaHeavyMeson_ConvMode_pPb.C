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
      eventCutArray = new TString[nMaxCuts]; photonCutArray = new TString[nMaxCuts]; mesonCutArray = new TString[nMaxCuts];
      for(Int_t i=0; i<nMaxCuts; i++) {eventCutArray[i] = ""; photonCutArray[i] = ""; mesonCutArray[i] = "";}
    }

    void AddCut(TString eventCut, TString photonCut, TString mesonCut){
      if(nCuts>=nMaxCuts) {cout << "ERROR in CutHandlerHeavyMesonConv: Exceeded maximum number of cuts!" << endl; validCuts = false; return;}
      if( eventCut.Length()!=8 || photonCut.Length()!=26 || mesonCut.Length()!=16 ) {cout << "ERROR in CutHandlerHeavyMesonConv: Incorrect length of cut string!" << endl; validCuts = false; return;}
      eventCutArray[nCuts]=eventCut; photonCutArray[nCuts]=photonCut; mesonCutArray[nCuts]=mesonCut;
      nCuts++;
      return;
    }
    Bool_t AreValid(){return validCuts;}
    Int_t GetNCuts(){if(validCuts) return nCuts; else return 0;}
    TString GetEventCut(Int_t i){if(validCuts&&i<nMaxCuts&&i>=0) return eventCutArray[i]; else{cout << "ERROR in CutHandlerHeavyMesonConv: GetEventCut wrong index i" << endl;return "";}}
    TString GetPhotonCut(Int_t i){if(validCuts&&i<nMaxCuts&&i>=0) return photonCutArray[i]; else {cout << "ERROR in CutHandlerHeavyMesonConv: GetPhotonCut wrong index i" << endl;return "";}}
    TString GetMesonCut(Int_t i){if(validCuts&&i<nMaxCuts&&i>=0) return mesonCutArray[i]; else {cout << "ERROR in CutHandlerHeavyMesonConv: GetMesonCut wrong index i" << endl;return "";}}
  private:
    Bool_t validCuts;
    Int_t nCuts; Int_t nMaxCuts;
    TString* eventCutArray;
    TString* photonCutArray;
    TString* mesonCutArray;
};

//***************************************************************************************
//main function
//***************************************************************************************
void AddTask_GammaHeavyMeson_ConvMode_pPb(  Int_t     selectedMeson                 = 0,                    // select the corresponding meson: 0 pi0, 1 eta, 2 eta'
                                            Int_t     trainConfig                   = 1,                    // change different set of cuts
                                            Int_t     isMC                          = 0,                    // run MC
                                            Int_t     enableQAMesonTask             = 0,                    // enable QA in AliAnalysisTaskGammaConvV1
                                            Int_t     enableQAPhotonTask            = 0,                    // enable additional QA task
                                            TString   fileNameInputForWeighting     = "MCSpectraInput.root",// path to file for weigting input / modified acceptance
                                            Int_t     doWeightingPart               = 0,                    // enable Weighting
                                            TString   generatorName                 = "DPMJET",             // generator Name
                                            TString   cutnumberAODBranch            = "800000006008400000001500000",  // cutnumber for AOD branch
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
  if(rAddConfigArr->GetEntries()<1){cout << "ERROR: AddTask_GammaHeavyMeson_ConvMode_pPb during parsing of additionalTrainConfig String '" << additionalTrainConfig.Data() << "'" << endl; return;}
  TObjString* rAdditionalTrainConfig = (TObjString*)rAddConfigArr->At(0);
  TString sAdditionalTrainConfig = rAdditionalTrainConfig->GetString();
  if (sAdditionalTrainConfig.Atoi() > 0){
    trainConfig = trainConfig + sAdditionalTrainConfig.Atoi();
    cout << "INFO: AddTask_GammaHeavyMeson_ConvMode_pPb running additionalTrainConfig '" << sAdditionalTrainConfig.Atoi() << "', train config: '" << trainConfig << "'" << endl;
  }
  cout << "corrTaskSetting: " << corrTaskSetting.Data() << endl;

  Int_t isHeavyIon    = 2;
  // meson reco mode: 0 - PCM-PCM, 1 - PCM-Calo, 2 - Calo-Calo
  Int_t mesonRecoMode = 0;

  // ================== GetAnalysisManager ===============================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error(Form("AddTask_GammaHeavyMeson_ConvMode_pPb_%i",trainConfig), "No analysis manager found.");
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
  TString cutnumberEvent    = "80000003";
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
  task= new AliAnalysisTaskHeavyNeutralMesonToGG(Form("HeavyNeutralMesonToGG_%i",trainConfig));
  task->SetIsHeavyIon(isHeavyIon);
  task->SetIsMC(isMC);
  task->SetV0ReaderName(V0ReaderName);
  task->SetCorrectionTaskSetting(corrTaskSetting);

  task->SetMesonRecoMode(mesonRecoMode); // meson reco mode: 0 - PCM-PCM, 1 - PCM-Calo, 2 - Calo-Calo
  if (runLightOutput > 1) task->SetLightOutput(kTRUE);

  //create cut handler
  CutHandlerHeavyMesonConv cuts;

  // Run1 default
  if (trainConfig == 1){
    cuts.AddCut("80010113","00200009327000008250400000","0163103100000010");

  // Run2 default
  } else if (trainConfig == 100){
    cuts.AddCut("80010113","00200009327000008250400000","0163103100000010");
    cuts.AddCut("80210113","00200009327000008250400000","0163103100000010");
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

  Int_t numberOfCuts = cuts.GetNCuts();

  TList *HeaderList     = new TList();
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
  ConvCutList->SetOwner(kTRUE);
  AliConversionPhotonCuts **analysisCuts      = new AliConversionPhotonCuts*[numberOfCuts];
  MesonCutList->SetOwner(kTRUE);
  AliConversionMesonCuts **analysisMesonCuts  = new AliConversionMesonCuts*[numberOfCuts];

  for(Int_t i = 0; i<numberOfCuts; i++){
    analysisEventCuts[i] = new AliConvEventCuts();

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

    analysisCuts[i] = new AliConversionPhotonCuts();
    analysisCuts[i]->SetV0ReaderName(V0ReaderName);
    if (runLightOutput > 0) analysisCuts[i]->SetLightOutput(kTRUE);
    analysisCuts[i]->InitializeCutsFromCutString((cuts.GetPhotonCut(i)).Data());
    analysisCuts[i]->SetIsHeavyIon(isHeavyIon);
    ConvCutList->Add(analysisCuts[i]);
    analysisCuts[i]->SetFillCutHistograms("",kFALSE);

    analysisMesonCuts[i] = new AliConversionMesonCuts();
    if (runLightOutput > 0) analysisMesonCuts[i]->SetLightOutput(kTRUE);
    analysisMesonCuts[i]->SetRunningMode(2);
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
  task->SetMesonType(selectedMeson);
  task->SetDoMesonQA(enableQAMesonTask); //Attention new switch for Pi0 QA
  task->SetDoPhotonQA(enableQAPhotonTask);  //Attention new switch small for Photon QA
  task->SetDoClusterQA(1);  //Attention new switch small for Cluster QA
  task->SetUseTHnSparse(isUsingTHnSparse);
  if(enableExtMatchAndQA > 1){ task->SetPlotHistsExtQA(kTRUE);}

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
