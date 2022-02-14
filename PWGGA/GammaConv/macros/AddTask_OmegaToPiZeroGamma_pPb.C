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
//($ALIPHYSICS/PWGGA/GammaConv/AliAnalysisTaskOmegaToPiZeroGamma.cxx) for
//pPb together with all supporting classes
//***************************************************************************************

//***************************************************************************************
//CutHandler contains all cuts for a certain analysis and trainconfig,
//it automatically checks length of cutStrings and takes care of the number of added cuts,
//no specification of the variable 'numberOfCuts' needed anymore.
//***************************************************************************************
class CutHandlerPiZeroGamma{
  public:
    CutHandlerPiZeroGamma(Int_t nMax=10){
      nCuts=0; nMaxCuts=nMax; validCuts = true;
      eventCutArray = new TString[nMaxCuts]; photonCutArray = new TString[nMaxCuts]; clusterCutArray = new TString[nMaxCuts]; neutralPionCutArray = new TString[nMaxCuts]; mesonCutArray = new TString[nMaxCuts];
      for(Int_t i=0; i<nMaxCuts; i++) {eventCutArray[i] = ""; photonCutArray[i] = ""; clusterCutArray[i] = ""; neutralPionCutArray[i] = ""; mesonCutArray[i] = "";}
    }

    void AddCut(TString eventCut, TString photonCut, TString clusterCut, TString neutralPionCut, TString mesonCut){
      if(nCuts>=nMaxCuts) {cout << "ERROR in CutHandlerPiZeroGamma: Exceeded maximum number of cuts!" << endl; validCuts = false; return;}
      if( eventCut.Length()!=8 || photonCut.Length()!=26 || clusterCut.Length()!=19 || neutralPionCut.Length()!=16 || mesonCut.Length()!=16 ) {cout << "ERROR in CutHandlerPiZeroGamma: Incorrect length of cut string!" << endl; validCuts = false; return;}
      eventCutArray[nCuts]=eventCut; photonCutArray[nCuts]=photonCut; clusterCutArray[nCuts]=clusterCut; neutralPionCutArray[nCuts]=neutralPionCut; mesonCutArray[nCuts]=mesonCut;
      nCuts++;
      return;
    }
    Bool_t AreValid(){return validCuts;}
    Int_t GetNCuts(){if(validCuts) return nCuts; else return 0;}
    TString GetEventCut(Int_t i){if(validCuts&&i<nMaxCuts&&i>=0) return eventCutArray[i]; else{cout << "ERROR in CutHandlerPiZeroGamma: GetEventCut wrong index i" << endl;return "";}}
    TString GetPhotonCut(Int_t i){if(validCuts&&i<nMaxCuts&&i>=0) return photonCutArray[i]; else {cout << "ERROR in CutHandlerPiZeroGamma: GetPhotonCut wrong index i" << endl;return "";}}
    TString GetClusterCut(Int_t i){if(validCuts&&i<nMaxCuts&&i>=0) return clusterCutArray[i]; else {cout << "ERROR in CutHandlerPiZeroGamma: GetClusterCut wrong index i" << endl;return "";}}
    TString GetNeutralPionCut(Int_t i){if(validCuts&&i<nMaxCuts&&i>=0) return neutralPionCutArray[i]; else {cout << "ERROR in CutHandlerPiZeroGamma: GetNeutralPionCut wrong index i" << endl;return "";}}
    TString GetMesonCut(Int_t i){if(validCuts&&i<nMaxCuts&&i>=0) return mesonCutArray[i]; else {cout << "ERROR in CutHandlerPiZeroGamma: GetMesonCut wrong index i" << endl;return "";}}
  private:
    Bool_t validCuts;
    Int_t nCuts; Int_t nMaxCuts;
    TString* eventCutArray;
    TString* photonCutArray;
    TString* clusterCutArray;
    TString* neutralPionCutArray;
    TString* mesonCutArray;
};

//***************************************************************************************
//main function
//***************************************************************************************
void AddTask_OmegaToPiZeroGamma_pPb( Int_t     trainConfig                 = 1,                    // change different set of cuts
                                Int_t     isMC                        = 0,                    // run MC
                                Int_t     enableQAMesonTask           = 1,                    // enable QA in AliAnalysisTaskGammaConvV1
                                Int_t     enableQAPhotonTask          = 1,                    // enable additional QA task
                                Bool_t    DoPiZeroGammaAngleCut       = kFALSE,               // flag for enabling cut on pi0-gamma angle
                                Double_t  lowerFactor                 = 0.75,                 // scale factor for lower limit in pi0-gamma angle cut
                                Double_t  upperFactor                 = 2.5,                  // scale factor for upper limit in pi0-gamma angle cut
                                TString   fileNameInputForWeighting   = "MCSpectraInput.root",// path to file for weigting input
                                Int_t     doWeightingPart             = 0,                    // enable Weighting
                                TString   generatorName               = "DPMJET",             // generator Name
                                TString   cutnumberAODBranch          = "800000006008400000001500000",  // cutnumber for AOD branch
                                Int_t     enableExtMatchAndQA         = 0,                    // disabled (0), extMatch (1), extQA_noCellQA (2), extMatch+extQA_noCellQA (3), extQA+cellQA (4), extMatch+extQA+cellQA (5)
                                Bool_t    isUsingTHnSparse            = kTRUE,                // enable or disable usage of THnSparses for background estimation
                                Bool_t    enableV0findingEffi         = kFALSE,               // enables V0finding efficiency histograms
                                Int_t     enableTriggerMimicking      = 0,                    // enable trigger mimicking
                                Bool_t    enableTriggerOverlapRej     = kFALSE,               // enable trigger overlap rejection
                                TString   settingMaxFacPtHard         = "3.",                 // maximum factor between hardest jet and ptHard generated
                                TString   periodNameV0Reader          = "",                   // period Name for V0Reader
                                Bool_t    enableSortingMCLabels       = kTRUE,                // enable sorting for MC cluster labels
                                Bool_t    runLightOutput              = kFALSE,               // switch to run light output (only essential histograms for afterburner)
                                TString   additionalTrainConfig       = "0"                   // additional counter for trainconfig, this has to be always the last parameter
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
        cout << Form("INFO: AddTask_OmegaToPiZeroGamma_pPb will use running mode '%i' for the TrackMatcher!",trackMatcherRunningMode) << endl;
      }
    }
  }
  TString sAdditionalTrainConfig = rAdditionalTrainConfig->GetString();
  if (sAdditionalTrainConfig.Atoi() > 0){
    trainConfig = trainConfig + sAdditionalTrainConfig.Atoi();
    cout << "INFO: AddTask_OmegaToPiZeroGamma_pPb running additionalTrainConfig '" << sAdditionalTrainConfig.Atoi() << "', train config: '" << trainConfig << "'" << endl;
  }

  TObjArray *rmaxFacPtHardSetting = settingMaxFacPtHard.Tokenize("_");
  if(rmaxFacPtHardSetting->GetEntries()<1){cout << "ERROR: AddTask_OmegaToPiZeroGamma_pPb during parsing of settingMaxFacPtHard String '" << settingMaxFacPtHard.Data() << "'" << endl; return;}
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

  Int_t isHeavyIon = 2;

  // fReconMethod = first digit of trainconfig
  Int_t ReconMethod = trainConfig/100;

  // ================== GetAnalysisManager ===============================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error(Form("AddTask_OmegaToPiZeroGamma_pPb_%i_%i", trainConfig, isMC), "No analysis manager found.");
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
  AliAnalysisTaskOmegaToPiZeroGamma *task=NULL;
  task= new AliAnalysisTaskOmegaToPiZeroGamma(Form("OmegaToPiZeroGamma_pPb_%i_%i",trainConfig,isMC));
  task->SetIsHeavyIon(isHeavyIon);
  task->SetIsMC(isMC);
  task->SetV0ReaderName(V0ReaderName);
  task->SetReconMethod(ReconMethod);
  task->SetDoPiZeroGammaAngleCut(DoPiZeroGammaAngleCut);
  task->SetlowerFactor(lowerFactor);
  task->SetupperFactor(upperFactor);
  task->SetTrackMatcherRunningMode(trackMatcherRunningMode);

  //create cut handler
  CutHandlerPiZeroGamma cuts;

  // cluster cuts
  // 0 "ClusterType",  1 "EtaMin", 2 "EtaMax", 3 "PhiMin", 4 "PhiMax", 5 "DistanceToBadChannel", 6 "Timing", 7 "TrackMatching", 8 "ExoticCell",
  // 9 "MinEnergy", 10 "MinNCells", 11 "MinM02", 12 "MaxM02", 13 "MinM20", 14 "MaxM20", 15 "MaximumDispersion", 16 "NLM"

  //************************************************ EMCAL clusters **********************************************************
  if (trainConfig == 1){ // EMCAL clusters standard cuts
    cuts.AddCut("80000113","00200009327000008250400000","1111141053032230000","0163103100000010","0163103100000010");

  } else if (trainConfig == 2){ // same as std, for varying angle cut
    cuts.AddCut("80000113","00200009327000008250400000","1111141053032230000","0163103100000010","0163103100000010");

  } else if (trainConfig == 3){ // same as std, for varying angle cut
    cuts.AddCut("80000113","00200009327000008250400000","1111141053032230000","0163103100000010","0163103100000010");

  } else if (trainConfig == 4){ // same as std, for varying angle cut
    cuts.AddCut("80000113","00200009327000008250400000","1111141053032230000","0163103100000010","0163103100000010");

  } else if (trainConfig == 101){ // EMCAL clusters standard cuts
    cuts.AddCut("80000113","00200009327000008250400000","1111141053032230000","0163103100000010","0163103100000010");

  } else if (trainConfig == 102){ // same as std, for varying angle cut
    cuts.AddCut("80000113","00200009327000008250400000","1111141053032230000","0163103100000010","0163103100000010");

  } else if (trainConfig == 103){ // same as std, for varying angle cut
    cuts.AddCut("80000113","00200009327000008250400000","1111141053032230000","0163103100000010","0163103100000010");

  } else if (trainConfig == 104){ // same as std, for varying angle cut
    cuts.AddCut("80000113","00200009327000008250400000","1111141053032230000","0163103100000010","0163103100000010");

  } else if (trainConfig == 201){ // EMCAL clusters standard cuts
    cuts.AddCut("80000113","00200009327000008250400000","1111141053032230000","0163103100000010","0163103100000010");

  } else if (trainConfig == 202){ // same as std, for varying angle cut
    cuts.AddCut("80000113","00200009327000008250400000","1111141053032230000","0163103100000010","0163103100000010");

  } else if (trainConfig == 203){ // same as std, for varying angle cut
    cuts.AddCut("80000113","00200009327000008250400000","1111141053032230000","0163103100000010","0163103100000010");

  } else if (trainConfig == 204){ // same as std, for varying angle cut
    cuts.AddCut("80000113","00200009327000008250400000","1111141053032230000","0163103100000010","0163103100000010");

  } else if (trainConfig == 301){ // EMCAL clusters standard cuts
    cuts.AddCut("80000113","00200009327000008250400000","1111141053032230000","0163103100000010","0163103100000010");

  } else if (trainConfig == 302){ // same as std, for varying angle cut
    cuts.AddCut("80000113","00200009327000008250400000","1111141053032230000","0163103100000010","0163103100000010");

  } else if (trainConfig == 303){ // same as std, for varying angle cut
    cuts.AddCut("80000113","00200009327000008250400000","1111141053032230000","0163103100000010","0163103100000010");

  } else if (trainConfig == 304){ // same as std, for varying angle cut
    cuts.AddCut("80000113","00200009327000008250400000","1111141053032230000","0163103100000010","0163103100000010");

  } else if (trainConfig == 401){ // EMCAL clusters standard cuts
    cuts.AddCut("80000113","00200009327000008250400000","1111141053032230000","0163103100000010","0163103100000010");

  } else if (trainConfig == 402){ // same as std, for varying angle cut
    cuts.AddCut("80000113","00200009327000008250400000","1111141053032230000","0163103100000010","0163103100000010");

  } else if (trainConfig == 403){ // same as std, for varying angle cut
    cuts.AddCut("80000113","00200009327000008250400000","1111141053032230000","0163103100000010","0163103100000010");

  } else if (trainConfig == 404){ // same as std, for varying angle cut
    cuts.AddCut("80000113","00200009327000008250400000","1111141053032230000","0163103100000010","0163103100000010");

  } else if (trainConfig == 501){ // EMCAL clusters standard cuts
    cuts.AddCut("80000113","00200009327000008250400000","1111141053032230000","0163103100000010","0163103100000010");

  } else if (trainConfig == 502){ // same as std, for varying angle cut
    cuts.AddCut("80000113","00200009327000008250400000","1111141053032230000","0163103100000010","0163103100000010");

  } else if (trainConfig == 503){ // same as std, for varying angle cut
    cuts.AddCut("80000113","00200009327000008250400000","1111141053032230000","0163103100000010","0163103100000010");

  } else if (trainConfig == 504){ // same as std, for varying angle cut
    cuts.AddCut("80000113","00200009327000008250400000","1111141053032230000","0163103100000010","0163103100000010");

  } else {
    Error(Form("OmegaToPiZeroGamma_pPb_%i_%i", trainConfig, isMC), "wrong trainConfig variable no cuts have been specified for the configuration");
    return;
  }


  if(!cuts.AreValid()){
    cout << "\n\n****************************************************" << endl;
    cout << "ERROR: No valid cuts stored in CutHandlerPiZeroGamma! Returning..." << endl;
    cout << "****************************************************\n\n" << endl;
    return;
  }

  TList *EventCutList   = new TList();
  TList *ConvCutList    = new TList();
  TList *ClusterCutList = new TList();
  TList *neutralPionCutList = new TList();
  TList *MesonCutList   = new TList();

  Int_t numberOfCuts = cuts.GetNCuts();

  EventCutList->SetOwner(kTRUE);
  AliConvEventCuts **analysisEventCuts        = new AliConvEventCuts*[numberOfCuts];
  ConvCutList->SetOwner(kTRUE);
  AliConversionPhotonCuts **analysisCuts      = new AliConversionPhotonCuts*[numberOfCuts];
  ClusterCutList->SetOwner(kTRUE);
  AliCaloPhotonCuts **analysisClusterCuts     = new AliCaloPhotonCuts*[numberOfCuts];
  neutralPionCutList->SetOwner(kTRUE);
  AliConversionMesonCuts **analysisNeutralPionCuts = new AliConversionMesonCuts*[numberOfCuts];
  MesonCutList->SetOwner(kTRUE);
  AliConversionMesonCuts **analysisMesonCuts  = new AliConversionMesonCuts*[numberOfCuts];

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
    analysisEventCuts[i]->InitializeCutsFromCutString((cuts.GetEventCut(i)).Data());
    EventCutList->Add(analysisEventCuts[i]);
    analysisEventCuts[i]->SetFillCutHistograms("",kFALSE);

    analysisCuts[i] = new AliConversionPhotonCuts();
    analysisCuts[i]->SetV0ReaderName(V0ReaderName);
    analysisCuts[i]->SetLightOutput(runLightOutput);
    analysisCuts[i]->InitializeCutsFromCutString((cuts.GetPhotonCut(i)).Data());
    analysisCuts[i]->SetIsHeavyIon(isHeavyIon);
    ConvCutList->Add(analysisCuts[i]);
    analysisCuts[i]->SetFillCutHistograms("",kFALSE);

    analysisClusterCuts[i] = new AliCaloPhotonCuts((isMC==2));
    analysisClusterCuts[i]->SetV0ReaderName(V0ReaderName);
    analysisClusterCuts[i]->SetCaloTrackMatcherName(TrackMatcherName);
    analysisClusterCuts[i]->SetLightOutput(runLightOutput);
    if(ReconMethod == 2) analysisClusterCuts[i]->SetIsPureCaloCut(2);
    analysisClusterCuts[i]->InitializeCutsFromCutString((cuts.GetClusterCut(i)).Data());
    ClusterCutList->Add(analysisClusterCuts[i]);
    analysisClusterCuts[i]->SetExtendedMatchAndQA(enableExtMatchAndQA);
    analysisClusterCuts[i]->SetFillCutHistograms("");

    analysisNeutralPionCuts[i] = new AliConversionMesonCuts();
    analysisNeutralPionCuts[i]->SetLightOutput(runLightOutput);
    analysisNeutralPionCuts[i]->SetRunningMode(2);
    analysisNeutralPionCuts[i]->InitializeCutsFromCutString((cuts.GetNeutralPionCut(i)).Data());
    neutralPionCutList->Add(analysisNeutralPionCuts[i]);
    analysisNeutralPionCuts[i]->SetFillCutHistograms("NeutralPionCuts");

    analysisMesonCuts[i] = new AliConversionMesonCuts();
    analysisMesonCuts[i]->SetLightOutput(runLightOutput);
    analysisMesonCuts[i]->SetRunningMode(2);
    analysisMesonCuts[i]->InitializeCutsFromCutString((cuts.GetMesonCut(i)).Data());
    MesonCutList->Add(analysisMesonCuts[i]);
    analysisMesonCuts[i]->SetFillCutHistograms("MesonCuts");
  }

  task->SetEventCutList(numberOfCuts,EventCutList);
  task->SetConversionCutList(numberOfCuts,ConvCutList);
  task->SetCaloCutList(numberOfCuts,ClusterCutList);
  task->SetNeutralPionCutList(numberOfCuts,neutralPionCutList);
  task->SetMesonCutList(numberOfCuts,MesonCutList);
  task->SetMoveParticleAccordingToVertex(kTRUE);
  task->SetDoMesonQA(enableQAMesonTask); //Attention new switch for Pi0 QA
  task->SetDoPhotonQA(enableQAPhotonTask);  //Attention new switch small for Photon QA
  task->SetEnableSortingOfMCClusLabels(enableSortingMCLabels);
  if(enableExtMatchAndQA > 1){ task->SetPlotHistsExtQA(kTRUE);}

  //connect containers
  AliAnalysisDataContainer *coutput =
    mgr->CreateContainer(Form("OmegaToPiZeroGamma_pPb_%i_%i",trainConfig,isMC), TList::Class(),
              AliAnalysisManager::kOutputContainer,Form("OmegaToPiZeroGamma_pPb_%i_%i.root",trainConfig,isMC));

  mgr->AddTask(task);
  mgr->ConnectInput(task,0,cinput);
  mgr->ConnectOutput(task,1,coutput);

  return;

}
