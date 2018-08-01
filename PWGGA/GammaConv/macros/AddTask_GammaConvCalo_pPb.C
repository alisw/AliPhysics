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
//($ALIPHYSICS/PWGGA/GammaConv/AliAnalysisTaskGammaConvCalo.cxx) for
//pPb together with all supporting classes
//***************************************************************************************

//***************************************************************************************
//CutHandler contains all cuts for a certain analysis and trainconfig,
//it automatically checks length of cutStrings and takes care of the number of added cuts,
// no specification of the variable 'numberOfCuts' needed anymore.
//***************************************************************************************
class CutHandlerConvCalo{
  public:
    CutHandlerConvCalo(Int_t nMax=10){
      nCuts=0; nMaxCuts=nMax; validCuts = true;
      eventCutArray = new TString[nMaxCuts]; photonCutArray = new TString[nMaxCuts]; clusterCutArray = new TString[nMaxCuts]; mesonCutArray = new TString[nMaxCuts];
      for(Int_t i=0; i<nMaxCuts; i++) {eventCutArray[i] = ""; photonCutArray[i] = ""; clusterCutArray[i] = ""; mesonCutArray[i] = "";}
    }

    void AddCut(TString eventCut, TString photonCut, TString clusterCut, TString mesonCut){
      if(nCuts>=nMaxCuts) {cout << "ERROR in CutHandlerConvCalo: Exceeded maximum number of cuts!" << endl; validCuts = false; return;}
      if( eventCut.Length()!=8 || photonCut.Length()!=26 || clusterCut.Length()!=19 || mesonCut.Length()!=16 ) {cout << "ERROR in CutHandlerConvCalo: Incorrect length of cut string!" << endl; validCuts = false; return;}
      eventCutArray[nCuts]=eventCut; photonCutArray[nCuts]=photonCut; clusterCutArray[nCuts]=clusterCut; mesonCutArray[nCuts]=mesonCut;
      nCuts++;
      return;
    }
    Bool_t AreValid(){return validCuts;}
    Int_t GetNCuts(){if(validCuts) return nCuts; else return 0;}
    TString GetEventCut(Int_t i){if(validCuts&&i<nMaxCuts&&i>=0) return eventCutArray[i]; else{cout << "ERROR in CutHandlerConvCalo: GetEventCut wrong index i" << endl;return "";}}
    TString GetPhotonCut(Int_t i){if(validCuts&&i<nMaxCuts&&i>=0) return photonCutArray[i]; else {cout << "ERROR in CutHandlerConvCalo: GetPhotonCut wrong index i" << endl;return "";}}
    TString GetClusterCut(Int_t i){if(validCuts&&i<nMaxCuts&&i>=0) return clusterCutArray[i]; else {cout << "ERROR in CutHandlerConvCalo: GetClusterCut wrong index i" << endl;return "";}}
    TString GetMesonCut(Int_t i){if(validCuts&&i<nMaxCuts&&i>=0) return mesonCutArray[i]; else {cout << "ERROR in CutHandlerConvCalo: GetMesonCut wrong index i" << endl;return "";}}
  private:
    Bool_t validCuts;
    Int_t nCuts; Int_t nMaxCuts;
    TString* eventCutArray;
    TString* photonCutArray;
    TString* clusterCutArray;
    TString* mesonCutArray;
};

//***************************************************************************************
//main function
//***************************************************************************************
void AddTask_GammaConvCalo_pPb( Int_t     trainConfig                   = 1,                    // change different set of cuts
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
                                Bool_t    enableSortingMCLabels         = kTRUE,                // enable sorting for MC cluster labels
                                Int_t     runLightOutput                = 0,                    // switch to run light output 0 (disabled), 1 (for CutClasses), 2 (for cutClasses and task)
                                Bool_t    doPrimaryTrackMatching        = kTRUE,                // enable basic track matching for all primary tracks to cluster
                                TString   additionalTrainConfig         = "0"                   // additional counter for trainconfig, this has to be always the last parameter
) {

  Bool_t doTreeClusterShowerShape = kFALSE; // enable tree for meson cand EMCal shower shape studies
  TH1S* histoAcc = 0x0;                     // histo for modified acceptance
  TString corrTaskSetting = ""; // select which correction task setting to use
  //parse additionalTrainConfig flag
  TObjArray *rAddConfigArr = additionalTrainConfig.Tokenize("_");
  if(rAddConfigArr->GetEntries()<1){cout << "ERROR: AddTask_GammaConvCalo_pPb during parsing of additionalTrainConfig String '" << additionalTrainConfig.Data() << "'" << endl; return;}
  TObjString* rAdditionalTrainConfig;
  for(Int_t i = 0; i<rAddConfigArr->GetEntries() ; i++){
    if(i==0) rAdditionalTrainConfig = (TObjString*)rAddConfigArr->At(i);
    else{
      TObjString* temp = (TObjString*) rAddConfigArr->At(i);
      TString tempStr = temp->GetString();
      if(tempStr.CompareTo("INVMASSCLUSTree") == 0){
        cout << "INFO: AddTask_GammaConvCalo_pPb activating 'INVMASSCLUSTree'" << endl;
        doTreeClusterShowerShape = kTRUE;
      }else if(tempStr.BeginsWith("MODIFYACC")){
        cout << "INFO: AddTask_GammaConvCalo_pPb activating 'MODIFYACC'" << endl;
        TString tempType = tempStr;
        tempType.Replace(0,9,"");
        cout << "INFO: connecting to alien..." << endl;
        TGrid::Connect("alien://");
        cout << "done!" << endl;
        TFile *w = TFile::Open(fileNameInputForWeighting.Data());
        if(!w){cout << "ERROR: Could not open file: " << fileNameInputForWeighting.Data() << endl;return;}
        histoAcc = (TH1S*) w->Get(tempType.Data());
        if(!histoAcc) {cout << "ERROR: Could not find histo: " << tempType.Data() << endl;return;}
        cout << "found: " << histoAcc << endl;
      }else if(tempStr.BeginsWith("CF")){
        cout << "INFO: AddTask_GammaConvCalo_pPb will use custom branch from Correction Framework!" << endl;
        corrTaskSetting = tempStr;
        corrTaskSetting.Replace(0,2,"");
      }
    }
  }
  TString sAdditionalTrainConfig = rAdditionalTrainConfig->GetString();
  if (sAdditionalTrainConfig.Atoi() > 0){
    trainConfig = trainConfig + sAdditionalTrainConfig.Atoi();
    cout << "INFO: AddTask_GammaConvCalo_pPb running additionalTrainConfig '" << sAdditionalTrainConfig.Atoi() << "', train config: '" << trainConfig << "'" << endl;
  }
  cout << "corrTaskSetting: " << corrTaskSetting.Data() << endl;

  Int_t isHeavyIon = 2;

  // ================== GetAnalysisManager ===============================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error(Form("AddTask_GammaConvCalo_pPb_%i",trainConfig), "No analysis manager found.");
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
  AliAnalysisTaskGammaConvCalo *task=NULL;
  task= new AliAnalysisTaskGammaConvCalo(Form("GammaConvCalo_%i",trainConfig));
  task->SetIsHeavyIon(isHeavyIon);
  task->SetIsMC(isMC);
  task->SetV0ReaderName(V0ReaderName);
  task->SetCorrectionTaskSetting(corrTaskSetting);
  if (runLightOutput > 1) task->SetLightOutput(kTRUE);
  task->SetDoPrimaryTrackMatching(doPrimaryTrackMatching);

  //create cut handler
  CutHandlerConvCalo cuts;

  // cluster cuts
  // 0 "ClusterType",  1 "EtaMin", 2 "EtaMax", 3 "PhiMin", 4 "PhiMax", 5 "DistanceToBadChannel", 6 "Timing", 7 "TrackMatching", 8 "ExoticCell",
  // 9 "MinEnergy", 10 "MinNCells", 11 "MinM02", 12 "MaxM02", 13 "MinMaxM20", 14 "RecConv", 15 "MaximumDispersion", 16 "NLM"


  //************************************************ EMCAL clusters **********************************************************
  //---------------------------------------------------------------------------------------------
  // no non linearity cuts
  //---------------------------------------------------------------------------------------------
  if (trainConfig == 1){ // min energy = 0.3 GeV/c
    cuts.AddCut("80010113","00200009327000008250400000","1111100057022230000","0163103100000010"); //standart cut, kINT7
  } else if (trainConfig == 2){  // min energy = 0.3 GeV/c
    cuts.AddCut("80052113","00200009327000008250400000","1111100057022230000","0163103100000010"); //standard cut, kEMC7
    cuts.AddCut("80083113","00200009327000008250400000","1111100057022230000","0163103100000010"); //standard cut, kEMCEG1 based on INT7
    cuts.AddCut("80085113","00200009327000008250400000","1111100057022230000","0163103100000010"); //standard cut, kEMCEG2 based on INT7
  } else if (trainConfig == 3){ // min energy = 0.4 GeV/c
    cuts.AddCut("80010113","00200009327000008250400000","1111100057032230000","0163103100000010"); //standart cut, kINT7
  } else if (trainConfig == 4){ // min energy = 0.4 GeV/
    cuts.AddCut("80052113","00200009327000008250400000","1111100057032230000","0163103100000010"); //standard cut, kEMC7
    cuts.AddCut("80083113","00200009327000008250400000","1111100057032230000","0163103100000010"); //standard cut, kEMCEG1 based on INT7
    cuts.AddCut("80085113","00200009327000008250400000","1111100057032230000","0163103100000010"); //standard cut, kEMCEG2 based on INT7
  } else if (trainConfig == 5){ // min energy = 0.4 GeV/c w/o time cut
    cuts.AddCut("80010113","00200009327000008250400000","1111100007032230000","0163103100000010"); //standart cut, kINT7
  } else if (trainConfig == 6){ // min energy = 0.4 GeV/ w/o time cut
    cuts.AddCut("80052113","00200009327000008250400000","1111100007032230000","0163103100000010"); //standard cut, kEMC7
    cuts.AddCut("80083113","00200009327000008250400000","1111100007032230000","0163103100000010"); //standard cut, kEMCEG1 based on INT7
    cuts.AddCut("80085113","00200009327000008250400000","1111100007032230000","0163103100000010"); //standard cut, kEMCEG2 based on INT7
  } else if (trainConfig == 7){ // introduce fast future protection
    cuts.AddCut("80010113","00200009327000008250400000","1111100007032230000","0163103100000010"); // no PF
    cuts.AddCut("80010313","00200009327000008250400000","1111100007032230000","0163103100000010"); // 0.1 \mus protected
    cuts.AddCut("80010413","00200009327000008250400000","1111100007032230000","0163103100000010"); // 0.25 \mus protected
    cuts.AddCut("80010513","00200009327000008250400000","1111100007032230000","0163103100000010"); // 1.075 \mus protected
  } else if (trainConfig == 8){ // PF for cent var, 1.075 \mus protected
    cuts.AddCut("80210513","00200009327000008250400000","1111100007032230000","0163103100000010"); // 0-20
    cuts.AddCut("82410513","00200009327000008250400000","1111100007032230000","0163103100000010"); // 20-40
    cuts.AddCut("84610513","00200009327000008250400000","1111100007032230000","0163103100000010"); // 40-60
    cuts.AddCut("86010513","00200009327000008250400000","1111100007032230000","0163103100000010"); // 60-100
  //---------------------------------------------------------------------------------------------
  // standard cuts
  //---------------------------------------------------------------------------------------------
  } else if (trainConfig == 10){ // introduce fast future protection
    cuts.AddCut("80010113","00200009327000008250400000","1111141007032230000","0163103100000010"); // no PF
    cuts.AddCut("80010313","00200009327000008250400000","1111141007032230000","0163103100000010"); // 0.1 \mus protected
    cuts.AddCut("80010413","00200009327000008250400000","1111141007032230000","0163103100000010"); // 0.25 \mus protected
    cuts.AddCut("80010513","00200009327000008250400000","1111141007032230000","0163103100000010"); // 1.075 \mus protected
  } else if (trainConfig == 11) {
    cuts.AddCut("80010113","00200009327000008250400000","2444451044013200000","0163103100000010"); // standart cut, kINT7
    cuts.AddCut("80062113","00200009327000008250400000","2444451044013200000","0163103100000010"); // standard cut, kPHI7

  } else if (trainConfig == 12){ // EMCAL clusters standard cuts cent dependent
    cuts.AddCut("80210113","00200009327000008250400000","1111141057032230000","0163103100000010"); // 0-20
  } else if (trainConfig == 13){ // EMCAL clusters standard cuts cent dependent
    cuts.AddCut("82410113","00200009327000008250400000","1111141057032230000","0163103100000010"); // 20-40
  } else if (trainConfig == 14){ // EMCAL clusters standard cuts cent dependent
    cuts.AddCut("84610113","00200009327000008250400000","1111141057032230000","0163103100000010"); // 40-60
  } else if (trainConfig == 15){ // EMCAL clusters standard cuts cent dependent
    cuts.AddCut("86010113","00200009327000008250400000","1111141057032230000","0163103100000010"); // 60-100
  } else if (trainConfig == 16){ // EMCAL clusters standard cut MB
    cuts.AddCut("80010113","00200009327000008250400000","1111141057032230000","0163103100000010"); // INT7
  } else if (trainConfig == 17){ // EMCAL clusters standard cut MB + cent dependent
    cuts.AddCut("80010113","00200009327000008250400000","1111141057032230000","0163103100000010"); // 0-100
    cuts.AddCut("80210113","00200009327000008250400000","1111141057032230000","0163103100000010"); // 0-20
    cuts.AddCut("82410113","00200009327000008250400000","1111141057032230000","0163103100000010"); // 20-40
    cuts.AddCut("84610113","00200009327000008250400000","1111141057032230000","0163103100000010"); // 40-60
    cuts.AddCut("86010113","00200009327000008250400000","1111141057032230000","0163103100000010"); // 60-100
  } else if (trainConfig == 18){ // EMCAL clusters standard cut MB + cent dependent CL1 est
    cuts.AddCut("90010113","00200009327000008250400000","1111141057032230000","0163103100000010"); // 0-20
    cuts.AddCut("90210113","00200009327000008250400000","1111141057032230000","0163103100000010"); // 0-20
    cuts.AddCut("92410113","00200009327000008250400000","1111141057032230000","0163103100000010"); // 20-40
    cuts.AddCut("94610113","00200009327000008250400000","1111141057032230000","0163103100000010"); // 40-60
    cuts.AddCut("96010113","00200009327000008250400000","1111141057032230000","0163103100000010"); // 60-100
  } else if (trainConfig == 19){ // EMCAL clusters standard cut MB + cent dependent CL1 est
    cuts.AddCut("e0010113","00200009327000008250400000","1111141057032230000","0163103100000010"); // 0-20
    cuts.AddCut("e0210113","00200009327000008250400000","1111141057032230000","0163103100000010"); // 0-20
    cuts.AddCut("e2410113","00200009327000008250400000","1111141057032230000","0163103100000010"); // 20-40
    cuts.AddCut("e4610113","00200009327000008250400000","1111141057032230000","0163103100000010"); // 40-60
    cuts.AddCut("e6010113","00200009327000008250400000","1111141057032230000","0163103100000010"); // 60-100

  } else if (trainConfig == 20){ // EMCAL clusters standard cuts
    cuts.AddCut("80010113","00200009327000008250400000","1111141057032230000","0163103100000010"); // INT7
    cuts.AddCut("80052113","00200009327000008250400000","1111141057032230000","0163103100000010"); // EMC7
    cuts.AddCut("80085113","00200009327000008250400000","1111141057032230000","0163103100000010"); // EG2
    cuts.AddCut("80083113","00200009327000008250400000","1111141057032230000","0163103100000010"); // EG1

  //---------------------------------------------------------------------------------------------
  // minimum bias variations
  //---------------------------------------------------------------------------------------------
  } else if (trainConfig == 21){ //EMCAL minEnergy variation
    cuts.AddCut("80010113","00200009327000008250400000","1111141057032230000","0163103100000010"); // 0.7 GeV/c default
    cuts.AddCut("80010113","00200009327000008250400000","1111141057042230000","0163103100000010"); // 0.8 GeV/c
    cuts.AddCut("80010113","00200009327000008250400000","1111141057052230000","0163103100000010"); // 0.9 GeV/c
  } else if (trainConfig == 22){ //EMCAL minNCells variation
    cuts.AddCut("80010113","00200009327000008250400000","1111141057031230000","0163103100000010"); // n cells >= 1
    cuts.AddCut("80010113","00200009327000008250400000","1111141057033230000","0163103100000010"); // n cells >= 3
    cuts.AddCut("80010113","00200009327000008250400000","1111141057032200000","0163103100000010"); // no max M02 cut
    cuts.AddCut("80010113","00200009327000008250400000","1111141057032250000","0163103100000010"); // M02 < 0.3
    cuts.AddCut("80010113","00200009327000008250400000","1111141057032260000","0163103100000010"); // M02 < 0.27
    cuts.AddCut("80010113","00200009327000008250400000","1112141057032230000","0163103100000010"); // only modules with TRD infront
    cuts.AddCut("80010113","00200009327000008250400000","1111341057032230000","0163103100000010"); // no modules with TRD infront
  } else if (trainConfig == 23){ // EMCAL track matching variations
    cuts.AddCut("80010113","00200009327000008250400000","1111141053032230000","0163103100000010"); // fixed TM
    cuts.AddCut("80010113","00200009327000008250400000","1111141056032230000","0163103100000010"); // pt dep var 1
    cuts.AddCut("80010113","00200009327000008250400000","1111141058032230000","0163103100000010"); // pt dep var 1
    cuts.AddCut("80010113","00200009327000008250400000","1111141059032230000","0163103100000010"); // pt dep var 1
  } else if (trainConfig == 24){ // PCM variations
    cuts.AddCut("80010113","00200009227000008250400000","1111141057032230000","0163103100000010"); // dEdx e -3, 5
    cuts.AddCut("80010113","00200009127000008250400000","1111141057032230000","0163103100000010"); // dEdx e -5, 5
    cuts.AddCut("80010113","00200009357000008250400000","1111141057032230000","0163103100000010"); // dEdx pi 2
    cuts.AddCut("80010113","00200009317000008250400000","1111141057032230000","0163103100000010"); // dEdx pi 0
    cuts.AddCut("80010113","00200009387300008250400000","1111141057032230000","0163103100000010"); // dEdx pi 2 high 1 (> 3.5 GeV)
  } else if (trainConfig == 25){ // PCM variations
    cuts.AddCut("80010113","00200009327000009250400000","1111141057032230000","0163103100000010"); // qt 2D 0.03
    cuts.AddCut("80010113","00200009327000003250400000","1111141057032230000","0163103100000010"); // qt 1D 0.05
    cuts.AddCut("80010113","00200009327000002250400000","1111141057032230000","0163103100000010"); // qt 1D 0.07
    cuts.AddCut("80010113","00200049327000008250400000","1111141057032230000","0163103100000010"); // single pt > 0.075
    cuts.AddCut("80010113","00200019327000008250400000","1111141057032230000","0163103100000010"); // single pt > 0.1
  } else if (trainConfig == 26){ // PCM variations
    cuts.AddCut("80010113","00200009327000008850400000","1111141057032230000","0163103100000010"); // 2D psi pair chi2 var
    cuts.AddCut("80010113","00200009327000008260400000","1111141057032230000","0163103100000010"); // 2D psi pair chi2 var
    cuts.AddCut("80010113","00200009327000008860400000","1111141057032230000","0163103100000010"); // 2D psi pair chi2 var
    cuts.AddCut("80010113","00200009327000008280400000","1111141057032230000","0163103100000010"); // 2D psi pair chi2 var
    cuts.AddCut("80010113","00200009327000008880400000","1111141057032230000","0163103100000010"); // 2D psi pair chi2 var
  } else if (trainConfig == 27){ // PCM variations
    cuts.AddCut("80010113","00200006327000008250400000","1111141057032230000","0163103100000010"); // min TPC cl > 0.7
    cuts.AddCut("80010113","00200008327000008250400000","1111141057032230000","0163103100000010"); // min TPC cl > 0.35
    cuts.AddCut("80010113","00200009327000008250400000","1111141057032230000","0163106100000010"); // alpha < 0.8
    cuts.AddCut("80010113","00200009327000008250400000","1111141057032230000","0163105100000010"); // alpha < 0.75
  } else if (trainConfig == 28){ // PCM variations
    cuts.AddCut("80010113","00202209327000008250400000","1111141057032230000","0163103100000010"); // restrict acceptance to EMCAL loose
    cuts.AddCut("80010113","00204409327000008250400000","1111141057032230000","0163103100000010"); // restrict acceptance to EMCAL tight
  } else if (trainConfig == 29){ // PCM variations pi dEdx
    cuts.AddCut("80010113","00200009317300008250400000","1111141057032230000","0163103100000010"); // dEdx pi: 0: 0.4-3.5, -10: 3.5 ->
    cuts.AddCut("80010113","00200009327300008250400000","1111141057032230000","0163103100000010"); // dEdx pi: 1: 0.4-3.5, -10: 3.5 ->
    cuts.AddCut("80010113","00200009325000008250400000","1111141057032230000","0163103100000010"); // dEdx pi: 1: 0.3 ->
    cuts.AddCut("80010113","00200009320000008250400000","1111141057032230000","0163103100000010"); // dEdx pi: 1: 0.5 ->
  } else if (trainConfig == 30){ // PCM variations pi dEdx
    cuts.AddCut("80010113","00200009327600008250400000","1111141057032230000","0163103100000010"); // dEdx pi: 1: 0.4-2, -10: 2. ->
    cuts.AddCut("80010113","00200009327400008250400000","1111141057032230000","0163103100000010"); // dEdx pi: 1: 0.4-3, -10: 3. ->
    cuts.AddCut("80010113","00200009315600008250400000","1111141057032230000","0163103100000010"); // dEdx pi: 0: 0.3-2, -10: 2. ->
    cuts.AddCut("80010113","00200009367400008250400000","1111141057032230000","0163103100000010"); // dEdx pi: 2: 0.4-3, 0.5: 3. ->
    cuts.AddCut("80010113","00200009347400008250400000","1111141057032230000","0163103100000010"); // dEdx pi: 3: 0.4-3, 1: 3. ->
  } else if (trainConfig == 31){ // PCM variations to close V0s
    cuts.AddCut("80010113","00200009327000008250400000","1111141057032230000","0163103100000010");
    cuts.AddCut("80010113","00200009327000008250401000","1111141057032230000","0163103100000010");
    cuts.AddCut("80010113","00200009327000008250402000","1111141057032230000","0163103100000010");
    cuts.AddCut("80010113","00200009327000008250403000","1111141057032230000","0163103100000010");
  } else if (trainConfig == 32){ // EMCal cluster, non lin variations
    cuts.AddCut("80010113","00200009327000008250400000","1111141057032230000","0163103100000010");
    cuts.AddCut("80010113","00200009327000008250400000","1111142057032230000","0163103100000010");
    cuts.AddCut("80010113","00200009327000008250400000","1111151057032230000","0163103100000010");
    cuts.AddCut("80010113","00200009327000008250400000","1111152057032230000","0163103100000010");
  } else if (trainConfig == 33){  // EMCal cluster, non lin variations
    cuts.AddCut("80010113","00200009327000008250400000","1111100057032230000","0163103100000010");
    cuts.AddCut("80010113","00200009327000008250400000","1111101057032230000","0163103100000010");
    cuts.AddCut("80010113","00200009327000008250400000","1111102057032230000","0163103100000010");
  } else if (trainConfig == 34){ // EMCal cluster, non lin variations
    cuts.AddCut("80010113","00200009327000008250400000","1111143057032230000","0163103100000010");
    cuts.AddCut("80010113","00200009327000008250400000","1111144057032230000","0163103100000010");
    cuts.AddCut("80010113","00200009327000008250400000","1111153057032230000","0163103100000010");
    cuts.AddCut("80010113","00200009327000008250400000","1111154057032230000","0163103100000010");

  //---------------------------------------------------------------------------------------------
  // non lin variations triggers
  //---------------------------------------------------------------------------------------------
  } else if (trainConfig == 40){ // non linearity variations EMC7
    cuts.AddCut("80052113","00200009327000008250400000","1111100057032230000","0163103100000010"); // no non linearity
    cuts.AddCut("80052113","00200009327000008250400000","1111101057032230000","0163103100000010"); // kSDM
    cuts.AddCut("80052113","00200009327000008250400000","1111141057032230000","0163103100000010"); // conv calo
    cuts.AddCut("80052113","00200009327000008250400000","1111142057032230000","0163103100000010"); // calo
    cuts.AddCut("80052113","00200009327000008250400000","1111151057032230000","0163103100000010"); // conv calo
    cuts.AddCut("80052113","00200009327000008250400000","1111152057032230000","0163103100000010"); // calo
  } else if (trainConfig == 41){ // non linearity variations EG2
    cuts.AddCut("80085113","00200009327000008250400000","1111100057032230000","0163103100000010"); // no non linearity
    cuts.AddCut("80085113","00200009327000008250400000","1111101057032230000","0163103100000010"); // kSDM
    cuts.AddCut("80085113","00200009327000008250400000","1111141057032230000","0163103100000010"); // conv calo
    cuts.AddCut("80085113","00200009327000008250400000","1111142057032230000","0163103100000010"); // calo
    cuts.AddCut("80085113","00200009327000008250400000","1111151057032230000","0163103100000010"); // conv calo
    cuts.AddCut("80085113","00200009327000008250400000","1111152057032230000","0163103100000010"); // calo
  } else if (trainConfig == 42){ // non linearity variations EG1
    cuts.AddCut("80083113","00200009327000008250400000","1111100057032230000","0163103100000010"); // no non linearity
    cuts.AddCut("80083113","00200009327000008250400000","1111101057032230000","0163103100000010"); // kSDM
    cuts.AddCut("80083113","00200009327000008250400000","1111141057032230000","0163103100000010"); // conv calo
    cuts.AddCut("80083113","00200009327000008250400000","1111142057032230000","0163103100000010"); // calo
    cuts.AddCut("80083113","00200009327000008250400000","1111151057032230000","0163103100000010"); // conv calo
    cuts.AddCut("80083113","00200009327000008250400000","1111152057032230000","0163103100000010"); // calo

  //---------------------------------------------------------------------------------------------
  // 0-20 % variations
  //---------------------------------------------------------------------------------------------
  } else if (trainConfig == 60){ // EMCAL clusters standard cuts
    cuts.AddCut("80210113","00200009327000008250400000","1111141057032230000","0163103100000010"); // INT7
    cuts.AddCut("80252113","00200009327000008250400000","1111141057032230000","0163103100000010"); // EMC7
    cuts.AddCut("80285113","00200009327000008250400000","1111141057022230000","0163103100000010"); // EG2
    cuts.AddCut("80283113","00200009327000008250400000","1111141057022230000","0163103100000010"); // EG1
  } else if (trainConfig == 61){ //EMCAL minEnergy variation
    cuts.AddCut("80210113","00200009327000008250400000","1111141057032230000","0163103100000010"); // 0.6 GeV/c default
    cuts.AddCut("80210113","00200009327000008250400000","1111141057042230000","0163103100000010"); // 0.7 GeV/c
    cuts.AddCut("80210113","00200009327000008250400000","1111141057052230000","0163103100000010"); // 0.9 GeV/c
  } else if (trainConfig == 62){ //EMCAL minNCells variation
    cuts.AddCut("80210113","00200009327000008250400000","1111141057031230000","0163103100000010"); // n cells >= 1
    cuts.AddCut("80210113","00200009327000008250400000","1111141057033230000","0163103100000010"); // n cells >= 3
    cuts.AddCut("80210113","00200009327000008250400000","1111141057032200000","0163103100000010"); // no max M02 cut
    cuts.AddCut("80210113","00200009327000008250400000","1111141057032250000","0163103100000010"); // M02 < 0.3
    cuts.AddCut("80210113","00200009327000008250400000","1111141057032260000","0163103100000010"); // M02 < 0.27
    cuts.AddCut("80210113","00200009327000008250400000","1112141057032230000","0163103100000010"); // only modules with TRD infront
    cuts.AddCut("80210113","00200009327000008250400000","1111341057032230000","0163103100000010"); // no modules with TRD infront
  } else if (trainConfig == 63){ // EMCAL track matching variations
    cuts.AddCut("80210113","00200009327000008250400000","1111141053032230000","0163103100000010"); //
    cuts.AddCut("80210113","00200009327000008250400000","1111141056032230000","0163103100000010"); //
    cuts.AddCut("80210113","00200009327000008250400000","1111141058032230000","0163103100000010"); //
    cuts.AddCut("80210113","00200009327000008250400000","1111141059032230000","0163103100000010"); //
  } else if (trainConfig == 64){ // PCM variations
    cuts.AddCut("80210113","00200009227000008250400000","1111141057032230000","0163103100000010"); // dEdx e -3, 5
    cuts.AddCut("80210113","00200009127000008250400000","1111141057032230000","0163103100000010"); // dEdx e -5, 5
    cuts.AddCut("80210113","00200009357000008250400000","1111141057032230000","0163103100000010"); // dEdx pi 2
    cuts.AddCut("80210113","00200009317000008250400000","1111141057032230000","0163103100000010"); // dEdx pi 0
    cuts.AddCut("80210113","00200009387300008250400000","1111141057032230000","0163103100000010"); // dEdx pi 2 high 1 (> 3.5 GeV)
  } else if (trainConfig == 65){ // PCM variations
    cuts.AddCut("80210113","00200009327000009250400000","1111141057032230000","0163103100000010"); // qt 2D 0.03
    cuts.AddCut("80210113","00200009327000003250400000","1111141057032230000","0163103100000010"); // qt 1D 0.05
    cuts.AddCut("80210113","00200009327000002250400000","1111141057032230000","0163103100000010"); // qt 1D 0.07
    cuts.AddCut("80210113","00200049327000008250400000","1111141057032230000","0163103100000010"); // single pt > 0.075
    cuts.AddCut("80210113","00200019327000008250400000","1111141057032230000","0163103100000010"); // single pt > 0.1
  } else if (trainConfig == 66){ // PCM variations
    cuts.AddCut("80210113","00200009327000008850400000","1111141057032230000","0163103100000010"); // 2D psi pair chi2 var
    cuts.AddCut("80210113","00200009327000008260400000","1111141057032230000","0163103100000010"); // 2D psi pair chi2 var
    cuts.AddCut("80210113","00200009327000008860400000","1111141057032230000","0163103100000010"); // 2D psi pair chi2 var
    cuts.AddCut("80210113","00200009327000008280400000","1111141057032230000","0163103100000010"); // 2D psi pair chi2 var
    cuts.AddCut("80210113","00200009327000008880400000","1111141057032230000","0163103100000010"); // 2D psi pair chi2 var
  } else if (trainConfig == 67){ // PCM variations
    cuts.AddCut("80210113","00200006327000008250400000","1111141057032230000","0163103100000010"); // min TPC cl > 0.7
    cuts.AddCut("80210113","00200008327000008250400000","1111141057032230000","0163103100000010"); // min TPC cl > 0.35
    cuts.AddCut("80210113","00200009327000008250400000","1111141057032230000","0163106100000010"); // alpha < 0.8
    cuts.AddCut("80210113","00200009327000008250400000","1111141057032230000","0163105100000010"); // alpha < 0.75
  } else if (trainConfig == 68){ // PCM variations
    cuts.AddCut("80210113","00202209327000008250400000","1111141057032230000","0163103100000010"); // restrict acceptance to EMCAL loose
    cuts.AddCut("80210113","00204409327000008250400000","1111141057032230000","0163103100000010"); // restrict acceptance to EMCAL tight
  } else if (trainConfig == 69){ // PCM variations pi dEdx
    cuts.AddCut("80210113","00200009317300008250400000","1111141057032230000","0163103100000010"); // dEdx pi: 0: 0.4-3.5, -10: 3.5 ->
    cuts.AddCut("80210113","00200009327300008250400000","1111141057032230000","0163103100000010"); // dEdx pi: 1: 0.4-3.5, -10: 3.5 ->
    cuts.AddCut("80210113","00200009325000008250400000","1111141057032230000","0163103100000010"); // dEdx pi: 1: 0.3 ->
    cuts.AddCut("80210113","00200009320000008250400000","1111141057032230000","0163103100000010"); // dEdx pi: 1: 0.5 ->
  } else if (trainConfig == 70){ // PCM variations pi dEdx
    cuts.AddCut("80210113","00200009327600008250400000","1111141057032230000","0163103100000010"); // dEdx pi: 1: 0.4-2, -10: 2. ->
    cuts.AddCut("80210113","00200009327400008250400000","1111141057032230000","0163103100000010"); // dEdx pi: 1: 0.4-3, -10: 3. ->
    cuts.AddCut("80210113","00200009315600008250400000","1111141057032230000","0163103100000010"); // dEdx pi: 0: 0.3-2, -10: 2. ->
    cuts.AddCut("80210113","00200009367400008250400000","1111141057032230000","0163103100000010"); // dEdx pi: 2: 0.4-3, 0.5: 3. ->
    cuts.AddCut("80210113","00200009347400008250400000","1111141057032230000","0163103100000010"); // dEdx pi: 3: 0.4-3, 1: 3. ->
  } else if (trainConfig == 71){ // PCM variations to close V0s
    cuts.AddCut("80210113","00200009327000008250400000","1111141057032230000","0163103100000010");
    cuts.AddCut("80210113","00200009327000008250401000","1111141057032230000","0163103100000010");
    cuts.AddCut("80210113","00200009327000008250402000","1111141057032230000","0163103100000010");
    cuts.AddCut("80210113","00200009327000008250403000","1111141057032230000","0163103100000010");
  } else if (trainConfig == 72){ // EMCal cluster, non lin variations
    cuts.AddCut("80210113","00200009327000008250400000","1111142057032230000","0163103100000010");
    cuts.AddCut("80210113","00200009327000008250400000","1111143057032230000","0163103100000010");
    cuts.AddCut("80210113","00200009327000008250400000","1111151057032230000","0163103100000010");
    cuts.AddCut("80210113","00200009327000008250400000","1111152057032230000","0163103100000010");
  } else if (trainConfig == 73){  // EMCal cluster, non lin variations
    cuts.AddCut("80210113","00200009327000008250400000","1111100057032230000","0163103100000010");
    cuts.AddCut("80210113","00200009327000008250400000","1111101057032230000","0163103100000010");
    cuts.AddCut("80210113","00200009327000008250400000","1111102057032230000","0163103100000010");
    cuts.AddCut("80210113","00200009327000008250400000","1111144057032230000","0163103100000010");

  //---------------------------------------------------------------------------------------------
  // 20-40 % variations
  //---------------------------------------------------------------------------------------------
  } else if (trainConfig == 80){ // EMCAL clusters standard cuts
    cuts.AddCut("82410113","00200009327000008250400000","1111141057032230000","0163103100000010"); // INT7
    cuts.AddCut("82452113","00200009327000008250400000","1111141057032230000","0163103100000010"); // EMC7
    cuts.AddCut("82485113","00200009327000008250400000","1111141057022230000","0163103100000010"); // EG2
    cuts.AddCut("82483113","00200009327000008250400000","1111141057022230000","0163103100000010"); // EG1
  } else if (trainConfig == 81){ //EMCAL minEnergy variation
    cuts.AddCut("82410113","00200009327000008250400000","1111141057032230000","0163103100000010"); // 0.7 GeV/c default
    cuts.AddCut("82410113","00200009327000008250400000","1111141057042230000","0163103100000010"); // 0.8 GeV/c
    cuts.AddCut("82410113","00200009327000008250400000","1111141057052230000","0163103100000010"); // 0.9 GeV/c
  } else if (trainConfig == 82){ //EMCAL minNCells variation
    cuts.AddCut("82410113","00200009327000008250400000","1111141057031230000","0163103100000010"); // n cells >= 1
    cuts.AddCut("82410113","00200009327000008250400000","1111141057033230000","0163103100000010"); // n cells >= 3
    cuts.AddCut("82410113","00200009327000008250400000","1111141057032200000","0163103100000010"); // no max M02 cut
    cuts.AddCut("82410113","00200009327000008250400000","1111141057032250000","0163103100000010"); // M02 < 0.3
    cuts.AddCut("82410113","00200009327000008250400000","1111141057032260000","0163103100000010"); // M02 < 0.27
    cuts.AddCut("82410113","00200009327000008250400000","1112141057032230000","0163103100000010"); // only modules with TRD infront
    cuts.AddCut("82410113","00200009327000008250400000","1111341057032230000","0163103100000010"); // no modules with TRD infront
  } else if (trainConfig == 83){ // EMCAL track matching variations
    cuts.AddCut("82410113","00200009327000008250400000","1111141053032230000","0163103100000010"); //
    cuts.AddCut("82410113","00200009327000008250400000","1111141056032230000","0163103100000010"); //
    cuts.AddCut("82410113","00200009327000008250400000","1111141058032230000","0163103100000010"); //
    cuts.AddCut("82410113","00200009327000008250400000","1111141059032230000","0163103100000010"); //
  } else if (trainConfig == 84){ // PCM variations
    cuts.AddCut("82410113","00200009227000008250400000","1111141057032230000","0163103100000010"); // dEdx e -3, 5
    cuts.AddCut("82410113","00200009127000008250400000","1111141057032230000","0163103100000010"); // dEdx e -5, 5
    cuts.AddCut("82410113","00200009357000008250400000","1111141057032230000","0163103100000010"); // dEdx pi 2
    cuts.AddCut("82410113","00200009317000008250400000","1111141057032230000","0163103100000010"); // dEdx pi 0
    cuts.AddCut("82410113","00200009387300008250400000","1111141057032230000","0163103100000010"); // dEdx pi 2 high 1 (> 3.5 GeV)
  } else if (trainConfig == 85){ // PCM variations
    cuts.AddCut("82410113","00200009327000009250400000","1111141057032230000","0163103100000010"); // qt 2D 0.03
    cuts.AddCut("82410113","00200009327000003250400000","1111141057032230000","0163103100000010"); // qt 1D 0.05
    cuts.AddCut("82410113","00200009327000002250400000","1111141057032230000","0163103100000010"); // qt 1D 0.07
    cuts.AddCut("82410113","00200049327000008250400000","1111141057032230000","0163103100000010"); // single pt > 0.075
    cuts.AddCut("82410113","00200019327000008250400000","1111141057032230000","0163103100000010"); // single pt > 0.1
  } else if (trainConfig == 86){ // PCM variations
    cuts.AddCut("82410113","00200009327000008850400000","1111141057032230000","0163103100000010"); // 2D psi pair chi2 var
    cuts.AddCut("82410113","00200009327000008260400000","1111141057032230000","0163103100000010"); // 2D psi pair chi2 var
    cuts.AddCut("82410113","00200009327000008860400000","1111141057032230000","0163103100000010"); // 2D psi pair chi2 var
    cuts.AddCut("82410113","00200009327000008280400000","1111141057032230000","0163103100000010"); // 2D psi pair chi2 var
    cuts.AddCut("82410113","00200009327000008880400000","1111141057032230000","0163103100000010"); // 2D psi pair chi2 var
  } else if (trainConfig == 87){ // PCM variations
    cuts.AddCut("82410113","00200006327000008250400000","1111141057032230000","0163103100000010"); // min TPC cl > 0.7
    cuts.AddCut("82410113","00200008327000008250400000","1111141057032230000","0163103100000010"); // min TPC cl > 0.35
    cuts.AddCut("82410113","00200009327000008250400000","1111141057032230000","0163106100000010"); // alpha < 0.8
    cuts.AddCut("82410113","00200009327000008250400000","1111141057032230000","0163105100000010"); // alpha < 0.75
  } else if (trainConfig == 88){ // PCM variations
    cuts.AddCut("82410113","00202209327000008250400000","1111141057032230000","0163103100000010"); // restrict acceptance to EMCAL loose
    cuts.AddCut("82410113","00204409327000008250400000","1111141057032230000","0163103100000010"); // restrict acceptance to EMCAL tight
  } else if (trainConfig == 89){ // PCM variations pi dEdx
    cuts.AddCut("82410113","00200009317300008250400000","1111141057032230000","0163103100000010"); // dEdx pi: 0: 0.4-3.5, -10: 3.5 ->
    cuts.AddCut("82410113","00200009327300008250400000","1111141057032230000","0163103100000010"); // dEdx pi: 1: 0.4-3.5, -10: 3.5 ->
    cuts.AddCut("82410113","00200009325000008250400000","1111141057032230000","0163103100000010"); // dEdx pi: 1: 0.3 ->
    cuts.AddCut("82410113","00200009320000008250400000","1111141057032230000","0163103100000010"); // dEdx pi: 1: 0.5 ->
  } else if (trainConfig == 90){ // PCM variations pi dEdx
    cuts.AddCut("82410113","00200009327600008250400000","1111141057032230000","0163103100000010"); // dEdx pi: 1: 0.4-2, -10: 2. ->
    cuts.AddCut("82410113","00200009327400008250400000","1111141057032230000","0163103100000010"); // dEdx pi: 1: 0.4-3, -10: 3. ->
    cuts.AddCut("82410113","00200009315600008250400000","1111141057032230000","0163103100000010"); // dEdx pi: 0: 0.3-2, -10: 2. ->
    cuts.AddCut("82410113","00200009367400008250400000","1111141057032230000","0163103100000010"); // dEdx pi: 2: 0.4-3, 0.5: 3. ->
    cuts.AddCut("82410113","00200009347400008250400000","1111141057032230000","0163103100000010"); // dEdx pi: 3: 0.4-3, 1: 3. ->
  } else if (trainConfig == 91){ // PCM variations to close V0s
    cuts.AddCut("82410113","00200009327000008250400000","1111141057032230000","0163103100000010");
    cuts.AddCut("82410113","00200009327000008250401000","1111141057032230000","0163103100000010");
    cuts.AddCut("82410113","00200009327000008250402000","1111141057032230000","0163103100000010");
    cuts.AddCut("82410113","00200009327000008250403000","1111141057032230000","0163103100000010");
  } else if (trainConfig == 92){ // EMCal cluster, non lin variations
    cuts.AddCut("82410113","00200009327000008250400000","1111142057032230000","0163103100000010");
    cuts.AddCut("82410113","00200009327000008250400000","1111143057032230000","0163103100000010");
    cuts.AddCut("82410113","00200009327000008250400000","1111151057032230000","0163103100000010");
    cuts.AddCut("82410113","00200009327000008250400000","1111152057032230000","0163103100000010");
  } else if (trainConfig == 93){  // EMCal cluster, non lin variations
    cuts.AddCut("82410113","00200009327000008250400000","1111100057032230000","0163103100000010");
    cuts.AddCut("82410113","00200009327000008250400000","1111101057032230000","0163103100000010");
    cuts.AddCut("82410113","00200009327000008250400000","1111102057032230000","0163103100000010");
    cuts.AddCut("82410113","00200009327000008250400000","1111144057032230000","0163103100000010");

  //---------------------------------------------------------------------------------------------
  // 40-60 % variations
  //---------------------------------------------------------------------------------------------
  } else if (trainConfig == 100){ // EMCAL clusters standard cuts
    cuts.AddCut("84610113","00200009327000008250400000","1111141057032230000","0163103100000010"); // INT7
    cuts.AddCut("84652113","00200009327000008250400000","1111141057032230000","0163103100000010"); // EMC7
    cuts.AddCut("84685113","00200009327000008250400000","1111141057022230000","0163103100000010"); // EG2
    cuts.AddCut("84683113","00200009327000008250400000","1111141057022230000","0163103100000010"); // EG1
  } else if (trainConfig == 101){ //EMCAL minEnergy variation
    cuts.AddCut("84610113","00200009327000008250400000","1111141057032230000","0163103100000010"); // 0.7 GeV/c default
    cuts.AddCut("84610113","00200009327000008250400000","1111141057042230000","0163103100000010"); // 0.8 GeV/c
    cuts.AddCut("84610113","00200009327000008250400000","1111141057052230000","0163103100000010"); // 0.9 GeV/c
  } else if (trainConfig == 102){ //EMCAL minNCells variation
    cuts.AddCut("84610113","00200009327000008250400000","1111141057031230000","0163103100000010"); // n cells >= 1
    cuts.AddCut("84610113","00200009327000008250400000","1111141057033230000","0163103100000010"); // n cells >= 3
    cuts.AddCut("84610113","00200009327000008250400000","1111141057032200000","0163103100000010"); // no max M02 cut
    cuts.AddCut("84610113","00200009327000008250400000","1111141057032250000","0163103100000010"); // M02 < 0.3
    cuts.AddCut("84610113","00200009327000008250400000","1111141057032260000","0163103100000010"); // M02 < 0.27
    cuts.AddCut("84610113","00200009327000008250400000","1112141057032230000","0163103100000010"); // only modules with TRD infront
    cuts.AddCut("84610113","00200009327000008250400000","1111341057032230000","0163103100000010"); // no modules with TRD infront
  } else if (trainConfig == 103){ // EMCAL track matching variations
    cuts.AddCut("84610113","00200009327000008250400000","1111141053032230000","0163103100000010"); //
    cuts.AddCut("84610113","00200009327000008250400000","1111141056032230000","0163103100000010"); //
    cuts.AddCut("84610113","00200009327000008250400000","1111141058032230000","0163103100000010"); //
    cuts.AddCut("84610113","00200009327000008250400000","1111141059032230000","0163103100000010"); //
  } else if (trainConfig == 104){ // PCM variations
    cuts.AddCut("84610113","00200009227000008250400000","1111141057032230000","0163103100000010"); // dEdx e -3, 5
    cuts.AddCut("84610113","00200009127000008250400000","1111141057032230000","0163103100000010"); // dEdx e -5, 5
    cuts.AddCut("84610113","00200009357000008250400000","1111141057032230000","0163103100000010"); // dEdx pi 2
    cuts.AddCut("84610113","00200009317000008250400000","1111141057032230000","0163103100000010"); // dEdx pi 0
    cuts.AddCut("84610113","00200009387300008250400000","1111141057032230000","0163103100000010"); // dEdx pi 2 high 1 (> 3.5 GeV)
  } else if (trainConfig == 105){ // PCM variations
    cuts.AddCut("84610113","00200009327000009250400000","1111141057032230000","0163103100000010"); // qt 2D 0.03
    cuts.AddCut("84610113","00200009327000003250400000","1111141057032230000","0163103100000010"); // qt 1D 0.05
    cuts.AddCut("84610113","00200009327000002250400000","1111141057032230000","0163103100000010"); // qt 1D 0.07
    cuts.AddCut("84610113","00200049327000008250400000","1111141057032230000","0163103100000010"); // single pt > 0.075
    cuts.AddCut("84610113","00200019327000008250400000","1111141057032230000","0163103100000010"); // single pt > 0.1
  } else if (trainConfig == 106){ // PCM variations
    cuts.AddCut("84610113","00200009327000008850400000","1111141057032230000","0163103100000010"); // 2D psi pair chi2 var
    cuts.AddCut("84610113","00200009327000008260400000","1111141057032230000","0163103100000010"); // 2D psi pair chi2 var
    cuts.AddCut("84610113","00200009327000008860400000","1111141057032230000","0163103100000010"); // 2D psi pair chi2 var
    cuts.AddCut("84610113","00200009327000008280400000","1111141057032230000","0163103100000010"); // 2D psi pair chi2 var
    cuts.AddCut("84610113","00200009327000008880400000","1111141057032230000","0163103100000010"); // 2D psi pair chi2 var
  } else if (trainConfig == 107){ // PCM variations
    cuts.AddCut("84610113","00200006327000008250400000","1111141057032230000","0163103100000010"); // min TPC cl > 0.7
    cuts.AddCut("84610113","00200008327000008250400000","1111141057032230000","0163103100000010"); // min TPC cl > 0.35
    cuts.AddCut("84610113","00200009327000008250400000","1111141057032230000","0163106100000010"); // alpha < 0.8
    cuts.AddCut("84610113","00200009327000008250400000","1111141057032230000","0163105100000010"); // alpha < 0.75
  } else if (trainConfig == 108){ // PCM variations
    cuts.AddCut("84610113","00202209327000008250400000","1111141057032230000","0163103100000010"); // restrict acceptance to EMCAL loose
    cuts.AddCut("84610113","00204409327000008250400000","1111141057032230000","0163103100000010"); // restrict acceptance to EMCAL tight
  } else if (trainConfig == 109){ // PCM variations pi dEdx
    cuts.AddCut("84610113","00200009317300008250400000","1111141057032230000","0163103100000010"); // dEdx pi: 0: 0.4-3.5, -10: 3.5 ->
    cuts.AddCut("84610113","00200009327300008250400000","1111141057032230000","0163103100000010"); // dEdx pi: 1: 0.4-3.5, -10: 3.5 ->
    cuts.AddCut("84610113","00200009325000008250400000","1111141057032230000","0163103100000010"); // dEdx pi: 1: 0.3 ->
    cuts.AddCut("84610113","00200009320000008250400000","1111141057032230000","0163103100000010"); // dEdx pi: 1: 0.5 ->
  } else if (trainConfig == 110){ // PCM variations pi dEdx
    cuts.AddCut("84610113","00200009327600008250400000","1111141057032230000","0163103100000010"); // dEdx pi: 1: 0.4-2, -10: 2. ->
    cuts.AddCut("84610113","00200009327400008250400000","1111141057032230000","0163103100000010"); // dEdx pi: 1: 0.4-3, -10: 3. ->
    cuts.AddCut("84610113","00200009315600008250400000","1111141057032230000","0163103100000010"); // dEdx pi: 0: 0.3-2, -10: 2. ->
    cuts.AddCut("84610113","00200009367400008250400000","1111141057032230000","0163103100000010"); // dEdx pi: 2: 0.4-3, 0.5: 3. ->
    cuts.AddCut("84610113","00200009347400008250400000","1111141057032230000","0163103100000010"); // dEdx pi: 3: 0.4-3, 1: 3. ->
  } else if (trainConfig == 111){ // PCM variations to close V0s
    cuts.AddCut("84610113","00200009327000008250400000","1111141057032230000","0163103100000010");
    cuts.AddCut("84610113","00200009327000008250401000","1111141057032230000","0163103100000010");
    cuts.AddCut("84610113","00200009327000008250402000","1111141057032230000","0163103100000010");
    cuts.AddCut("84610113","00200009327000008250403000","1111141057032230000","0163103100000010");
  } else if (trainConfig == 112){ // EMCal cluster, non lin variations
    cuts.AddCut("84610113","00200009327000008250400000","1111142057032230000","0163103100000010");
    cuts.AddCut("84610113","00200009327000008250400000","1111143057032230000","0163103100000010");
    cuts.AddCut("84610113","00200009327000008250400000","1111151057032230000","0163103100000010");
    cuts.AddCut("84610113","00200009327000008250400000","1111152057032230000","0163103100000010");
  } else if (trainConfig == 113){  // EMCal cluster, non lin variations
    cuts.AddCut("84610113","00200009327000008250400000","1111100057032230000","0163103100000010");
    cuts.AddCut("84610113","00200009327000008250400000","1111101057032230000","0163103100000010");
    cuts.AddCut("84610113","00200009327000008250400000","1111102057032230000","0163103100000010");
    cuts.AddCut("84610113","00200009327000008250400000","1111144057032230000","0163103100000010");

  //---------------------------------------------------------------------------------------------
  // 60-100 % variations
  //---------------------------------------------------------------------------------------------
  } else if (trainConfig == 120){ // EMCAL clusters standard cuts
    cuts.AddCut("86010113","00200009327000008250400000","1111141057032230000","0163103100000010"); // INT7
    cuts.AddCut("86052113","00200009327000008250400000","1111141057032230000","0163103100000010"); // EMC7
    cuts.AddCut("86085113","00200009327000008250400000","1111141057022230000","0163103100000010"); // EG2
    cuts.AddCut("86083113","00200009327000008250400000","1111141057022230000","0163103100000010"); // EG1
  } else if (trainConfig == 121){ //EMCAL minEnergy variation
    cuts.AddCut("86010113","00200009327000008250400000","1111141057032230000","0163103100000010"); // 0.7 GeV/c default
    cuts.AddCut("86010113","00200009327000008250400000","1111141057042230000","0163103100000010"); // 0.8 GeV/c
    cuts.AddCut("86010113","00200009327000008250400000","1111141057052230000","0163103100000010"); // 0.9 GeV/c
  } else if (trainConfig == 122){ //EMCAL minNCells variation
    cuts.AddCut("86010113","00200009327000008250400000","1111141057031230000","0163103100000010"); // n cells >= 1
    cuts.AddCut("86010113","00200009327000008250400000","1111141057033230000","0163103100000010"); // n cells >= 3
    cuts.AddCut("86010113","00200009327000008250400000","1111141057032200000","0163103100000010"); // no max M02 cut
    cuts.AddCut("86010113","00200009327000008250400000","1111141057032250000","0163103100000010"); // M02 < 0.3
    cuts.AddCut("86010113","00200009327000008250400000","1111141057032260000","0163103100000010"); // M02 < 0.27
    cuts.AddCut("86010113","00200009327000008250400000","1112141057032230000","0163103100000010"); // only modules with TRD infront
    cuts.AddCut("86010113","00200009327000008250400000","1111341057032230000","0163103100000010"); // no modules with TRD infront
  } else if (trainConfig == 123){ // EMCAL track matching variations
    cuts.AddCut("86010113","00200009327000008250400000","1111141053032230000","0163103100000010"); //
    cuts.AddCut("86010113","00200009327000008250400000","1111141056032230000","0163103100000010"); //
    cuts.AddCut("86010113","00200009327000008250400000","1111141058032230000","0163103100000010"); //
    cuts.AddCut("86010113","00200009327000008250400000","1111141059032230000","0163103100000010"); //
  } else if (trainConfig == 124){ // PCM variations
    cuts.AddCut("86010113","00200009227000008250400000","1111141057032230000","0163103100000010"); // dEdx e -3, 5
    cuts.AddCut("86010113","00200009127000008250400000","1111141057032230000","0163103100000010"); // dEdx e -5, 5
    cuts.AddCut("86010113","00200009357000008250400000","1111141057032230000","0163103100000010"); // dEdx pi 2
    cuts.AddCut("86010113","00200009317000008250400000","1111141057032230000","0163103100000010"); // dEdx pi 0
    cuts.AddCut("86010113","00200009387300008250400000","1111141057032230000","0163103100000010"); // dEdx pi 2 high 1 (> 3.5 GeV)
  } else if (trainConfig == 125){ // PCM variations
    cuts.AddCut("86010113","00200009327000009250400000","1111141057032230000","0163103100000010"); // qt 2D 0.03
    cuts.AddCut("86010113","00200009327000003250400000","1111141057032230000","0163103100000010"); // qt 1D 0.05
    cuts.AddCut("86010113","00200009327000002250400000","1111141057032230000","0163103100000010"); // qt 1D 0.07
    cuts.AddCut("86010113","00200049327000008250400000","1111141057032230000","0163103100000010"); // single pt > 0.075
    cuts.AddCut("86010113","00200019327000008250400000","1111141057032230000","0163103100000010"); // single pt > 0.1
  } else if (trainConfig == 126){ // PCM variations
    cuts.AddCut("86010113","00200009327000008850400000","1111141057032230000","0163103100000010"); // 2D psi pair chi2 var
    cuts.AddCut("86010113","00200009327000008260400000","1111141057032230000","0163103100000010"); // 2D psi pair chi2 var
    cuts.AddCut("86010113","00200009327000008860400000","1111141057032230000","0163103100000010"); // 2D psi pair chi2 var
    cuts.AddCut("86010113","00200009327000008280400000","1111141057032230000","0163103100000010"); // 2D psi pair chi2 var
    cuts.AddCut("86010113","00200009327000008880400000","1111141057032230000","0163103100000010"); // 2D psi pair chi2 var
  } else if (trainConfig == 127){ // PCM variations
    cuts.AddCut("86010113","00200006327000008250400000","1111141057032230000","0163103100000010"); // min TPC cl > 0.7
    cuts.AddCut("86010113","00200008327000008250400000","1111141057032230000","0163103100000010"); // min TPC cl > 0.35
    cuts.AddCut("86010113","00200009327000008250400000","1111141057032230000","0163106100000010"); // alpha < 0.8
    cuts.AddCut("86010113","00200009327000008250400000","1111141057032230000","0163105100000010"); // alpha < 0.75
  } else if (trainConfig == 128){ // PCM variations
    cuts.AddCut("86010113","00202209327000008250400000","1111141057032230000","0163103100000010"); // restrict acceptance to EMCAL loose
    cuts.AddCut("86010113","00204409327000008250400000","1111141057032230000","0163103100000010"); // restrict acceptance to EMCAL tight
  } else if (trainConfig == 129){ // PCM variations pi dEdx
    cuts.AddCut("86010113","00200009317300008250400000","1111141057032230000","0163103100000010"); // dEdx pi: 0: 0.4-3.5, -10: 3.5 ->
    cuts.AddCut("86010113","00200009327300008250400000","1111141057032230000","0163103100000010"); // dEdx pi: 1: 0.4-3.5, -10: 3.5 ->
    cuts.AddCut("86010113","00200009325000008250400000","1111141057032230000","0163103100000010"); // dEdx pi: 1: 0.3 ->
    cuts.AddCut("86010113","00200009320000008250400000","1111141057032230000","0163103100000010"); // dEdx pi: 1: 0.5 ->
  } else if (trainConfig == 130){ // PCM variations pi dEdx
    cuts.AddCut("86010113","00200009327600008250400000","1111141057032230000","0163103100000010"); // dEdx pi: 1: 0.4-2, -10: 2. ->
    cuts.AddCut("86010113","00200009327400008250400000","1111141057032230000","0163103100000010"); // dEdx pi: 1: 0.4-3, -10: 3. ->
    cuts.AddCut("86010113","00200009315600008250400000","1111141057032230000","0163103100000010"); // dEdx pi: 0: 0.3-2, -10: 2. ->
    cuts.AddCut("86010113","00200009367400008250400000","1111141057032230000","0163103100000010"); // dEdx pi: 2: 0.4-3, 0.5: 3. ->
    cuts.AddCut("86010113","00200009347400008250400000","1111141057032230000","0163103100000010"); // dEdx pi: 3: 0.4-3, 1: 3. ->
  } else if (trainConfig == 131){ // PCM variations to close V0s
    cuts.AddCut("86010113","00200009327000008250400000","1111141057032230000","0163103100000010");
    cuts.AddCut("86010113","00200009327000008250401000","1111141057032230000","0163103100000010");
    cuts.AddCut("86010113","00200009327000008250402000","1111141057032230000","0163103100000010");
    cuts.AddCut("86010113","00200009327000008250403000","1111141057032230000","0163103100000010");
  } else if (trainConfig == 132){ // EMCal cluster, non lin variations
    cuts.AddCut("86010113","00200009327000008250400000","1111142057032230000","0163103100000010");
    cuts.AddCut("86010113","00200009327000008250400000","1111143057032230000","0163103100000010");
    cuts.AddCut("86010113","00200009327000008250400000","1111151057032230000","0163103100000010");
    cuts.AddCut("86010113","00200009327000008250400000","1111152057032230000","0163103100000010");
  } else if (trainConfig == 133){  // EMCal cluster, non lin variations
    cuts.AddCut("86010113","00200009327000008250400000","1111100057032230000","0163103100000010");
    cuts.AddCut("86010113","00200009327000008250400000","1111101057032230000","0163103100000010");
    cuts.AddCut("86010113","00200009327000008250400000","1111102057032230000","0163103100000010");
    cuts.AddCut("86010113","00200009327000008250400000","1111144057032230000","0163103100000010");

  // ===============================================================================================
  // Run 2 data EMC clusters pPb 5TeV
  // ===============================================================================================
  } else if (trainConfig == 200){ // EMCAL clusters standard cuts, triggers, no nonlin, open timing
    cuts.AddCut("80010113","00200009327000008250400000","1111100057032230000","0163103100000010"); // INT7
    cuts.AddCut("80010113","00200009327000008250400000","1111100017032230000","0163103100000010"); // INT7
  } else if (trainConfig == 201){
    cuts.AddCut("80052113","00200009327000008250400000","1111100057032230000","0163103100000010"); // EMC7
    cuts.AddCut("80085113","00200009327000008250400000","1111100057032230000","0163103100000010"); // EG2
    cuts.AddCut("80083113","00200009327000008250400000","1111100057032230000","0163103100000010"); // EG1
  } else if (trainConfig == 202){
    cuts.AddCut("80052113","00200009327000008250400000","1111100017032230000","0163103100000010"); // EMC7
    cuts.AddCut("80085113","00200009327000008250400000","1111100017032230000","0163103100000010"); // EG2
    cuts.AddCut("80083113","00200009327000008250400000","1111100017032230000","0163103100000010"); // EG1
  } else if (trainConfig == 203){ // EMCAL clusters standard cuts, triggers, no nonlin, open timing
    cuts.AddCut("80110113","00200009327000008250400000","1111100057032230000","0163103100000010"); // 0-10
    cuts.AddCut("81210113","00200009327000008250400000","1111100057032230000","0163103100000010"); // 10-20
    cuts.AddCut("82410113","00200009327000008250400000","1111100057032230000","0163103100000010"); // 20-40
    cuts.AddCut("84610113","00200009327000008250400000","1111100057032230000","0163103100000010"); // 40-60
    cuts.AddCut("86810113","00200009327000008250400000","1111100057032230000","0163103100000010"); // 60-80
    cuts.AddCut("88010113","00200009327000008250400000","1111100057032230000","0163103100000010"); // 80-100
  } else if (trainConfig == 204){ // EMCAL clusters standard cuts, triggers, no nonlin, open timing
    cuts.AddCut("80210113","00200009327000008250400000","1111100057032230000","0163103100000010"); // 0-20
    cuts.AddCut("86010113","00200009327000008250400000","1111100057032230000","0163103100000010"); // 60-100
    cuts.AddCut("a0110113","00200009327000008250400000","1111100057032230000","0163103100000010"); // 0-5
    cuts.AddCut("a1210113","00200009327000008250400000","1111100057032230000","0163103100000010"); // 5-10
  } else if (trainConfig == 205){ // EMCAL clusters standard cuts, triggers, no nonlin, open timing
    cuts.AddCut("80110113","00200009327000008250400000","1111100017032230000","0163103100000010"); // 0-10
    cuts.AddCut("81210113","00200009327000008250400000","1111100017032230000","0163103100000010"); // 10-20
    cuts.AddCut("82410113","00200009327000008250400000","1111100017032230000","0163103100000010"); // 20-40
    cuts.AddCut("84610113","00200009327000008250400000","1111100017032230000","0163103100000010"); // 40-60
    cuts.AddCut("86810113","00200009327000008250400000","1111100017032230000","0163103100000010"); // 60-80
    cuts.AddCut("88010113","00200009327000008250400000","1111100017032230000","0163103100000010"); // 80-100
  } else if (trainConfig == 206){ // EMCAL clusters standard cuts, triggers, no nonlin, open timing
    cuts.AddCut("80210113","00200009327000008250400000","1111100017032230000","0163103100000010"); // 0-20
    cuts.AddCut("86010113","00200009327000008250400000","1111100017032230000","0163103100000010"); // 60-100
    cuts.AddCut("a0110113","00200009327000008250400000","1111100017032230000","0163103100000010"); // 0-5
    cuts.AddCut("a1210113","00200009327000008250400000","1111100017032230000","0163103100000010"); // 5-10

  } else if (trainConfig == 207){ // EMCAL clusters standard cuts, nonlin variations
    cuts.AddCut("80010113","00200009327000008250400000","1111141057032230000","0163103100000010"); // 0-100
    cuts.AddCut("80010113","00200009327000008250400000","1111151057032230000","0163103100000010"); // 0-100
    cuts.AddCut("80010113","00200009327000008250400000","1111142057032230000","0163103100000010"); // 0-100
    cuts.AddCut("80010113","00200009327000008250400000","1111152057032230000","0163103100000010"); // 0-100
    cuts.AddCut("80010113","00200009327000008250400000","1111155057032230000","0163103100000010"); // 0-100
    cuts.AddCut("80010113","00200009327000008250400000","1111156057032230000","0163103100000010"); // 0-100
  } else if (trainConfig == 208){ // EMCAL clusters standard cuts, nonlin variations
    cuts.AddCut("80010113","00200009327000008250400000","1111100057032230000","0163103100000010"); // 0-100
    cuts.AddCut("80010113","00200009327000008250400000","1111102057032230000","0163103100000010"); // 0-100
  } else if (trainConfig == 209){ // EMCAL clusters standard cuts, nonlin variations
    cuts.AddCut("80010113","00200009327000008250400000","1111141017032230000","0163103100000010"); // 0-100
    cuts.AddCut("80010113","00200009327000008250400000","1111151017032230000","0163103100000010"); // 0-100
    cuts.AddCut("80010113","00200009327000008250400000","1111142017032230000","0163103100000010"); // 0-100
    cuts.AddCut("80010113","00200009327000008250400000","1111152017032230000","0163103100000010"); // 0-100
    cuts.AddCut("80010113","00200009327000008250400000","1111155017032230000","0163103100000010"); // 0-100
    cuts.AddCut("80010113","00200009327000008250400000","1111156017032230000","0163103100000010"); // 0-100
  } else if (trainConfig == 210){ // EMCAL clusters standard cuts, nonlin variations
    cuts.AddCut("80010113","00200009327000008250400000","1111100017032230000","0163103100000010"); // 0-100
    cuts.AddCut("80010113","00200009327000008250400000","1111102017032230000","0163103100000010"); // 0-100

  // non lin variations with cent
  } else if (trainConfig == 211){ // EMCAL clusters standard cuts, triggers, no nonlin, open timing
    cuts.AddCut("80110113","00200009327000008250400000","1111141057032230000","0163103100000010"); // 0-10
    cuts.AddCut("81210113","00200009327000008250400000","1111141057032230000","0163103100000010"); // 10-20
    cuts.AddCut("82410113","00200009327000008250400000","1111141057032230000","0163103100000010"); // 20-40
    cuts.AddCut("84610113","00200009327000008250400000","1111141057032230000","0163103100000010"); // 40-60
    cuts.AddCut("86810113","00200009327000008250400000","1111141057032230000","0163103100000010"); // 60-80
    cuts.AddCut("88010113","00200009327000008250400000","1111141057032230000","0163103100000010"); // 80-100
  } else if (trainConfig == 212){ // EMCAL clusters standard cuts, triggers, no nonlin, open timing
    cuts.AddCut("80210113","00200009327000008250400000","1111141057032230000","0163103100000010"); // 0-20
    cuts.AddCut("86010113","00200009327000008250400000","1111141057032230000","0163103100000010"); // 60-100
    cuts.AddCut("a0110113","00200009327000008250400000","1111141057032230000","0163103100000010"); // 0-5
    cuts.AddCut("a1210113","00200009327000008250400000","1111141057032230000","0163103100000010"); // 5-10
  } else if (trainConfig == 213){ // EMCAL clusters standard cuts, triggers, no nonlin, open timing
    cuts.AddCut("80110113","00200009327000008250400000","1111151057032230000","0163103100000010"); // 0-10
    cuts.AddCut("81210113","00200009327000008250400000","1111151057032230000","0163103100000010"); // 10-20
    cuts.AddCut("82410113","00200009327000008250400000","1111151057032230000","0163103100000010"); // 20-40
    cuts.AddCut("84610113","00200009327000008250400000","1111151057032230000","0163103100000010"); // 40-60
    cuts.AddCut("86810113","00200009327000008250400000","1111151057032230000","0163103100000010"); // 60-80
    cuts.AddCut("88010113","00200009327000008250400000","1111151057032230000","0163103100000010"); // 80-100
  } else if (trainConfig == 214){ // EMCAL clusters standard cuts, triggers, no nonlin, open timing
    cuts.AddCut("80210113","00200009327000008250400000","1111151057032230000","0163103100000010"); // 0-20
    cuts.AddCut("86010113","00200009327000008250400000","1111151057032230000","0163103100000010"); // 60-100
    cuts.AddCut("a0110113","00200009327000008250400000","1111151057032230000","0163103100000010"); // 0-5
    cuts.AddCut("a1210113","00200009327000008250400000","1111151057032230000","0163103100000010"); // 5-10

  // non lin variations with diff minE cut
  } else if (trainConfig == 215){ // EMCAL clusters nonlin variations, min E=0.6 GeV
    cuts.AddCut("80010113","00200009327000008250400000","1111100057022230000","0163103100000010"); // 0-100
    cuts.AddCut("80010113","00200009327000008250400000","1111141057022230000","0163103100000010"); // 0-100
    cuts.AddCut("80010113","00200009327000008250400000","1111151057022230000","0163103100000010"); // 0-100
    cuts.AddCut("80010113","00200009327000008250400000","1111142057022230000","0163103100000010"); // 0-100
    cuts.AddCut("80010113","00200009327000008250400000","1111152057022230000","0163103100000010"); // 0-100
  } else if (trainConfig == 216){ // EMCAL clusters nonlin variations, min E=0.65 GeV
    cuts.AddCut("80010113","00200009327000008250400000","11111000570b2230000","0163103100000010"); // 0-100
    cuts.AddCut("80010113","00200009327000008250400000","11111410570b2230000","0163103100000010"); // 0-100
    cuts.AddCut("80010113","00200009327000008250400000","11111510570b2230000","0163103100000010"); // 0-100
    cuts.AddCut("80010113","00200009327000008250400000","11111420570b2230000","0163103100000010"); // 0-100
    cuts.AddCut("80010113","00200009327000008250400000","11111520570b2230000","0163103100000010"); // 0-100
  } else if (trainConfig == 217){ // EMCAL clusters nonlin variations, min E=0.675 GeV
    cuts.AddCut("80010113","00200009327000008250400000","11111000570c2230000","0163103100000010"); // 0-100
    cuts.AddCut("80010113","00200009327000008250400000","11111410570c2230000","0163103100000010"); // 0-100
    cuts.AddCut("80010113","00200009327000008250400000","11111510570c2230000","0163103100000010"); // 0-100
    cuts.AddCut("80010113","00200009327000008250400000","11111420570c2230000","0163103100000010"); // 0-100
    cuts.AddCut("80010113","00200009327000008250400000","11111520570c2230000","0163103100000010"); // 0-100
  } else if (trainConfig == 218){ // EMCAL clusters nonlin variations, min E=0.625 GeV
    cuts.AddCut("80010113","00200009327000008250400000","11111000570d2230000","0163103100000010"); // 0-100
    cuts.AddCut("80010113","00200009327000008250400000","11111410570d2230000","0163103100000010"); // 0-100
    cuts.AddCut("80010113","00200009327000008250400000","11111510570d2230000","0163103100000010"); // 0-100
    cuts.AddCut("80010113","00200009327000008250400000","11111420570d2230000","0163103100000010"); // 0-100
    cuts.AddCut("80010113","00200009327000008250400000","11111520570d2230000","0163103100000010"); // 0-100

  // narrow cent variations
  } else if (trainConfig == 219){
    cuts.AddCut("c0110113","00200009327000008250400000","1111151057032230000","0163103100000010"); // 0-1
    cuts.AddCut("c0210113","00200009327000008250400000","1111151057032230000","0163103100000010"); // 0-2
  } else if (trainConfig == 220){
    cuts.AddCut("c0110113","00200009327000008250400000","1111141057032230000","0163103100000010"); // 0-1
    cuts.AddCut("c0210113","00200009327000008250400000","1111141057032230000","0163103100000010"); // 0-2
  // AOD validation
  } else if (trainConfig == 221){
    cuts.AddCut("80010113","00200009327000008250400000","1111151057032230000","0163103100000010"); // 0-100

  // ===============================================================================================
  // Run 1 data PHOS clusters
  // ===============================================================================================
  } else if (trainConfig == 301) {  // min energy = 0.3 GeV/c
    cuts.AddCut("80010113","00200009327000008250400000","2444400041013200000","0163103100000010"); // standart cut, kINT7
    cuts.AddCut("80062113","00200009327000008250400000","2444400041013200000","0163103100000010"); // standard cut, kPHI7
  } else if (trainConfig == 302) { //PHOS
    cuts.AddCut("80010113","00200009327000008250400000","2444400041013200000","0163103100000010");
    cuts.AddCut("80010113","00200009327000008250400000","2444400042013200000","0163103100000010");
    cuts.AddCut("80010113","00200009327000008250400000","2444400043013200000","0163103100000010");
  } else if (trainConfig == 303) { //PHOS
    cuts.AddCut("80010023","00200009327000008250400000","2444400041013200000","0163103100000010");
    cuts.AddCut("80010023","00200009327000008250400000","2444400042013200000","0163103100000010");
    cuts.AddCut("80010023","00200009327000008250400000","2444400043013200000","0163103100000010");
  } else if (trainConfig == 305) { //PHOS timing cut variations
    cuts.AddCut("80010113","00200009327000008250400000","2444400011013200000","0163103100000010"); // 1000ns
    cuts.AddCut("80010113","00200009327000008250400000","2444400031013200000","0163103100000010"); // 200ns
    cuts.AddCut("80010113","00200009327000008250400000","2444400051013200000","0163103100000010"); // 50ns
  } else if (trainConfig == 306) {  // Non lin variations PHOS INT7
    cuts.AddCut("80010113","00200009327000008250400000","2444401041013200000","0163103100000010"); // PHOS group default
    cuts.AddCut("80010113","00200009327000008250400000","2444451041013200000","0163103100000010"); // CCMF PHOS
    cuts.AddCut("80010113","00200009327000008250400000","2444452041013200000","0163103100000010"); // CMF PHOS
  } else if (trainConfig == 307) {  // Non lin variations PHOS PHI7
    cuts.AddCut("80062113","00200009327000008250400000","2444401041013200000","0163103100000010"); // PHOS group default
    cuts.AddCut("80062113","00200009327000008250400000","2444451041013200000","0163103100000010"); // CCMF PHOS
    cuts.AddCut("80062113","00200009327000008250400000","2444452041013200000","0163103100000010"); // CMF PHOS
  } else if (trainConfig == 308) {  // CCMF cent dependet
    cuts.AddCut("80210113","00200009327000008250400000","2444451041013200000","0163103100000010"); // 0-20
    cuts.AddCut("82410113","00200009327000008250400000","2444451041013200000","0163103100000010"); // 20-40
    cuts.AddCut("84610113","00200009327000008250400000","2444451041013200000","0163103100000010"); // 40-60
    cuts.AddCut("86010113","00200009327000008250400000","2444451041013200000","0163103100000010"); // 60-100
  } else if (trainConfig == 309) {  // PHOS default cent dependet
    cuts.AddCut("80210113","00200009327000008250400000","2444401041013200000","0163103100000010"); // 0-20
    cuts.AddCut("82410113","00200009327000008250400000","2444401041013200000","0163103100000010"); // 20-40
    cuts.AddCut("84610113","00200009327000008250400000","2444401041013200000","0163103100000010"); // 40-60
    cuts.AddCut("86010113","00200009327000008250400000","2444401041013200000","0163103100000010"); // 60-100

  // testing past future protection
  } else if (trainConfig == 310) {  // PHOS INT7 PF test
    cuts.AddCut("80010213","00200009327000008250400000","2444451041013200000","0163103100000010"); // PF 2.25 \mu
    cuts.AddCut("80010513","00200009327000008250400000","2444451041013200000","0163103100000010"); // PF 1.075 \mu
    cuts.AddCut("80010413","00200009327000008250400000","2444451041013200000","0163103100000010"); // PF 0.25 \mu
  } else if (trainConfig == 311) {  // PHOS INT7 PF test, open cluster time
    cuts.AddCut("80010213","00200009327000008250400000","2444451011013200000","0163103100000010"); // PF 2.25 \mu
    cuts.AddCut("80010513","00200009327000008250400000","2444451011013200000","0163103100000010"); // PF 1.075 \mu
    cuts.AddCut("80010413","00200009327000008250400000","2444451011013200000","0163103100000010"); // PF 0.25 \mu
  } else if (trainConfig == 312) {  // PHOS PHI7 PF test
    cuts.AddCut("80062213","00200009327000008250400000","2444451041013200000","0163103100000010"); // PF 2.25 \mu
    cuts.AddCut("80062513","00200009327000008250400000","2444451041013200000","0163103100000010"); // PF 1.075 \mu
    cuts.AddCut("80062413","00200009327000008250400000","2444451041013200000","0163103100000010"); // PF 0.25 \mu
  } else if (trainConfig == 313) {  // PHOS PHI7 PF test, open cluster time
    cuts.AddCut("80062213","00200009327000008250400000","2444451011013200000","0163103100000010"); // PF 2.25 \mu
    cuts.AddCut("80062513","00200009327000008250400000","2444451011013200000","0163103100000010"); // PF 1.075 \mu
    cuts.AddCut("80062413","00200009327000008250400000","2444451011013200000","0163103100000010"); // PF 0.25 \mu

    //PCM-PHOS systematics
  } else if(trainConfig == 320){//dEdx e-line variation and dE/dx pi-line variation
    cuts.AddCut("80010113","00200009127000008250400000","2444451041013200000","0163103100000010"); //-5 < sigma < 5
    cuts.AddCut("80010113","00200009227000008250400000","2444451041013200000","0163103100000010"); //-3 < sigma < 5
    cuts.AddCut("80010113","00200009327400008250400000","2444451041013200000","0163103100000010"); //1, -10, 0.4, 3
    cuts.AddCut("80010113","00200009367400008250400000","2444451041013200000","0163103100000010"); //2, 0.5, 0.4, 3
  } else if(trainConfig == 321){//dE/dx pi-line variation and single pT variation and 2D triangular chi2 and psi pair
    cuts.AddCut("80010113","00200009317400008250400000","2444451041013200000","0163103100000010"); //0, -10, 0.4, 3
    cuts.AddCut("80010113","00200049327000008250400000","2444451041013200000","0163103100000010"); //single pT 0.075 GeV/c
    cuts.AddCut("80010113","00200019327000008250400000","2444451041013200000","0163103100000010"); //single pT 0.1   GeV/c
    cuts.AddCut("80010113","00200009327000008850400000","2444451041013200000","0163103100000010"); //20 & 0.1
  } else if(trainConfig == 322){//2D triangular chi2 and psi pair
    cuts.AddCut("80010113","00200009327000008260400000","2444451041013200000","0163103100000010"); //30 & 0.05
    cuts.AddCut("80010113","00200009327000008860400000","2444451041013200000","0163103100000010"); //20 & 0.05
    cuts.AddCut("80010113","00200009327000008280400000","2444451041013200000","0163103100000010"); //30 & 0.2
    cuts.AddCut("80010113","00200009327000008880400000","2444451041013200000","0163103100000010"); //20 & 0.2
  } else if(trainConfig == 323){//min TPC clusters and //elliptic cut Qt and alpha
    cuts.AddCut("80010113","00200000327000008250400000","2444451041013200000","0163103100000010"); //0
    cuts.AddCut("80010113","00200008327000008250400000","2444451041013200000","0163103100000010"); //0.35
    cuts.AddCut("80010113","00200009327000009250400000","2444451041013200000","0163103100000010"); // qT 0.03 no quadratic
    cuts.AddCut("80010113","00200009327000003250400000","2444451041013200000","0163103100000010"); // qT 0.05 y  quadratic
  } else if(trainConfig == 324){// and min phi max phi and NL variations
    cuts.AddCut("80010113","00209909327000008250400000","2444451041013200000","0163103100000010"); //4.54 > phi > 5.59
    cuts.AddCut("80010113","00200009327000008250400000","2444452041013200000","0163103100000010"); // PHOS calo NL
    cuts.AddCut("80010113","00200009327000008250400000","2444401041013200000","0163103100000010"); // PHOS people NL

    //Cluster Cuts varation
  } else if(trainConfig == 325){ // first set of variations CLUSTER
    cuts.AddCut("80010113","00200009327000008250400000","2444451041073200000","0163103100000010"); // min energy 0.2 GeV
    cuts.AddCut("80010113","00200009327000008250400000","2444451041083200000","0163103100000010"); // min energy 0.4 GeV
    cuts.AddCut("80010113","00200009327000008250400000","2444451041023200000","0163103100000010"); // min energy 0.5 GeV
  } else if(trainConfig == 326){ // second set of variations CLUSTER
    cuts.AddCut("80010113","00200009327000008250400000","2444451041013270000","0163103100000010"); // min/max M02  0.1<M<1
    cuts.AddCut("80010113","00200009327000008250400000","2444451041013280000","0163103100000010"); // min/max M02  0.1<M<1
    cuts.AddCut("80010113","00200009327000008250400000","2444451041012200000","0163103100000010"); // min number 2 cells
    cuts.AddCut("80010113","00200009327000008250400000","2444451041014200000","0163103100000010"); // min number 4 cells
  } else if(trainConfig == 327){ // MESON
    cuts.AddCut("80010113","00200009327000008250400000","2444451041013200000","0163403100000010"); // rapidity variation  y<0.5
    cuts.AddCut("80010113","00200009327000008250400000","2444451041013200000","0163803100000010"); // rapidity variation  y<0.25
    cuts.AddCut("80010113","00200009327000008250400000","2444451041013200000","0163106100000010"); // alpha meson variation 1 0<alpha<0.8
    cuts.AddCut("80010113","00200009327000008250400000","2444451041013200000","0163105100000010"); // alpha meson variation 2 0<alpha<0.75
  } else if(trainConfig == 328){ // fourth set of variations
    cuts.AddCut("80010113","00200009327000008250400000","2444451044013200000","0163103100000010"); // tm variation
    cuts.AddCut("80010113","00200009327000008250400000","2444451045013200000","0163103100000010"); // tm variation
    cuts.AddCut("80010113","00200009327000008250400000","2444451046013200000","0163103100000010"); // tm variation
    cuts.AddCut("80010113","00200009327000008250400000","2444451041013200000","0163103100000000"); // min opening angle 0    -> open
    cuts.AddCut("80010113","00200009327000008250400000","2444451041013200000","0163103100000030"); // min opening angle 0.01 -> 2 cell diag
  } else if(trainConfig == 329){ // centrality dependent and with latest TM
    cuts.AddCut("80010113","00200009327000008250400000","2444451044013200000","0163103100000010"); // 0-100% with PCM NL
    cuts.AddCut("80210113","00200009327000008250400000","2444451044013200000","0163103100000010"); // 0-20% with PCM NL
    cuts.AddCut("82410113","00200009327000008250400000","2444451044013200000","0163103100000010"); // 20-40% with PCM NL
    cuts.AddCut("84610113","00200009327000008250400000","2444451044013200000","0163103100000010"); // 40-60% with PCM NL
    cuts.AddCut("86010113","00200009327000008250400000","2444451044013200000","0163103100000010"); // 60-100% with PCM NL

    //-----------------------------------------------------------------------------------------------
    // PCM-PHOS systematics run 1
    //-----------------------------------------------------------------------------------------------
    //0 - 20
  } else if(trainConfig == 330){//dEdx e-line variation and dE/dx pi-line variation
    cuts.AddCut("80210113","00200009127000008250400000","2444451041013200000","0163103100000010"); //-5 < sigma < 5
    cuts.AddCut("80210113","00200009227000008250400000","2444451041013200000","0163103100000010"); //-3 < sigma < 5
    cuts.AddCut("80210113","00200009327400008250400000","2444451041013200000","0163103100000010"); //1, -10, 0.4, 3
    cuts.AddCut("80210113","00200009367400008250400000","2444451041013200000","0163103100000010"); //2, 0.5, 0.4, 3
  } else if(trainConfig == 331){//dE/dx pi-line variation and single pT variation and 2D triangular chi2 and psi pair
    cuts.AddCut("80210113","00200009317400008250400000","2444451041013200000","0163103100000010"); //0, -10, 0.4, 3
    cuts.AddCut("80210113","00200049327000008250400000","2444451041013200000","0163103100000010"); //single pT 0.075 GeV/c
    cuts.AddCut("80210113","00200019327000008250400000","2444451041013200000","0163103100000010"); //single pT 0.1   GeV/c
    cuts.AddCut("80210113","00200009327000008850400000","2444451041013200000","0163103100000010"); //20 & 0.1
  } else if(trainConfig == 332){//2D triangular chi2 and psi pair
    cuts.AddCut("80210113","00200009327000008260400000","2444451041013200000","0163103100000010"); //30 & 0.05
    cuts.AddCut("80210113","00200009327000008860400000","2444451041013200000","0163103100000010"); //20 & 0.05
    cuts.AddCut("80210113","00200009327000008280400000","2444451041013200000","0163103100000010"); //30 & 0.2
    cuts.AddCut("80210113","00200009327000008880400000","2444451041013200000","0163103100000010"); //20 & 0.2
  } else if(trainConfig == 333){//min TPC clusters and //elliptic cut Qt and alpha
    cuts.AddCut("80210113","00200000327000008250400000","2444451041013200000","0163103100000010"); //0
    cuts.AddCut("80210113","00200008327000008250400000","2444451041013200000","0163103100000010"); //0.35
    cuts.AddCut("80210113","00200009327000009250400000","2444451041013200000","0163103100000010"); // qT 0.03 no quadratic
    cuts.AddCut("80210113","00200009327000003250400000","2444451041013200000","0163103100000010"); // qT 0.05 y  quadratic
  } else if(trainConfig == 333){// and min phi max phi and NL variations
    cuts.AddCut("80210113","00209909327000008250400000","2444451041013200000","0163103100000010"); //4.54 > phi > 5.59
    cuts.AddCut("80210113","00200009327000008250400000","2444452041013200000","0163103100000010"); // PHOS calo NL
    cuts.AddCut("80210113","00200009327000008250400000","2444401041013200000","0163103100000010"); // PHOS people NL

    //Cluster Cuts varation
  } else if(trainConfig == 334){ // first set of variations CLUSTER
    cuts.AddCut("80210113","00200009327000008250400000","2444451041073200000","0163103100000010"); // min energy 0.2 GeV
    cuts.AddCut("80210113","00200009327000008250400000","2444451041083200000","0163103100000010"); // min energy 0.4 GeV
    cuts.AddCut("80210113","00200009327000008250400000","2444451041023200000","0163103100000010"); // min energy 0.5 GeV
  } else if(trainConfig == 335){ // second set of variations CLUSTER
    cuts.AddCut("80210113","00200009327000008250400000","2444451041013270000","0163103100000010"); // min/max M02  0.1<M<1
    cuts.AddCut("80210113","00200009327000008250400000","2444451041013280000","0163103100000010"); // min/max M02  0.1<M<1
    cuts.AddCut("80210113","00200009327000008250400000","2444451041012200000","0163103100000010"); // min number 2 cells
    cuts.AddCut("80210113","00200009327000008250400000","2444451041014200000","0163103100000010"); // min number 4 cells
  } else if(trainConfig == 336){ // MESON
    cuts.AddCut("80210113","00200009327000008250400000","2444451041013200000","0163403100000010"); // rapidity variation  y<0.5
    cuts.AddCut("80210113","00200009327000008250400000","2444451041013200000","0163803100000010"); // rapidity variation  y<0.25
    cuts.AddCut("80210113","00200009327000008250400000","2444451041013200000","0163106100000010"); // alpha meson variation 1 0<alpha<0.8
    cuts.AddCut("80210113","00200009327000008250400000","2444451041013200000","0163105100000010"); // alpha meson variation 2 0<alpha<0.75
  } else if(trainConfig == 337){ // fourth set of variations
    cuts.AddCut("80210113","00200009327000008250400000","2444451044013200000","0163103100000010"); // tm variation
    cuts.AddCut("80210113","00200009327000008250400000","2444451045013200000","0163103100000010"); // tm variation
    cuts.AddCut("80210113","00200009327000008250400000","2444451046013200000","0163103100000010"); // tm variation
    cuts.AddCut("80210113","00200009327000008250400000","2444451041013200000","0163103100000000"); // min opening angle 0    -> open
    cuts.AddCut("80210113","00200009327000008250400000","2444451041013200000","0163103100000030"); // min opening angle 0.01 -> 2 cell diag
    //20 - 40
  } else if(trainConfig == 340){//dEdx e-line variation and dE/dx pi-line variation
    cuts.AddCut("82410113","00200009127000008250400000","2444451041013200000","0163103100000010"); //-5 < sigma < 5
    cuts.AddCut("82410113","00200009227000008250400000","2444451041013200000","0163103100000010"); //-3 < sigma < 5
    cuts.AddCut("82410113","00200009327400008250400000","2444451041013200000","0163103100000010"); //1, -10, 0.4, 3
    cuts.AddCut("82410113","00200009367400008250400000","2444451041013200000","0163103100000010"); //2, 0.5, 0.4, 3
  } else if(trainConfig == 341){//dE/dx pi-line variation and single pT variation and 2D triangular chi2 and psi pair
    cuts.AddCut("82410113","00200009317400008250400000","2444451041013200000","0163103100000010"); //0, -10, 0.4, 3
    cuts.AddCut("82410113","00200049327000008250400000","2444451041013200000","0163103100000010"); //single pT 0.075 GeV/c
    cuts.AddCut("82410113","00200019327000008250400000","2444451041013200000","0163103100000010"); //single pT 0.1   GeV/c
    cuts.AddCut("82410113","00200009327000008850400000","2444451041013200000","0163103100000010"); //20 & 0.1
  } else if(trainConfig == 342){//2D triangular chi2 and psi pair
    cuts.AddCut("82410113","00200009327000008260400000","2444451041013200000","0163103100000010"); //30 & 0.05
    cuts.AddCut("82410113","00200009327000008860400000","2444451041013200000","0163103100000010"); //20 & 0.05
    cuts.AddCut("82410113","00200009327000008280400000","2444451041013200000","0163103100000010"); //30 & 0.2
    cuts.AddCut("82410113","00200009327000008880400000","2444451041013200000","0163103100000010"); //20 & 0.2
  } else if(trainConfig == 343){//min TPC clusters and //elliptic cut Qt and alpha
    cuts.AddCut("82410113","00200000327000008250400000","2444451041013200000","0163103100000010"); //0
    cuts.AddCut("82410113","00200008327000008250400000","2444451041013200000","0163103100000010"); //0.35
    cuts.AddCut("82410113","00200009327000009250400000","2444451041013200000","0163103100000010"); // qT 0.03 no quadratic
    cuts.AddCut("82410113","00200009327000003250400000","2444451041013200000","0163103100000010"); // qT 0.05 y  quadratic
  } else if(trainConfig == 344){// and min phi max phi and NL variations
    cuts.AddCut("82410113","00209909327000008250400000","2444451041013200000","0163103100000010"); //4.54 > phi > 5.59
    cuts.AddCut("82410113","00200009327000008250400000","2444452041013200000","0163103100000010"); // PHOS calo NL
    cuts.AddCut("82410113","00200009327000008250400000","2444401041013200000","0163103100000010"); // PHOS people NL

    //Cluster Cuts varation
  } else if(trainConfig == 345){ // first set of variations CLUSTER
    cuts.AddCut("82410113","00200009327000008250400000","2444451041073200000","0163103100000010"); // min energy 0.2 GeV
    cuts.AddCut("82410113","00200009327000008250400000","2444451041083200000","0163103100000010"); // min energy 0.4 GeV
    cuts.AddCut("82410113","00200009327000008250400000","2444451041023200000","0163103100000010"); // min energy 0.5 GeV
  } else if(trainConfig == 346){ // second set of variations CLUSTER
    cuts.AddCut("82410113","00200009327000008250400000","2444451041013270000","0163103100000010"); // min/max M02  0.1<M<1
    cuts.AddCut("82410113","00200009327000008250400000","2444451041013280000","0163103100000010"); // min/max M02  0.1<M<1
    cuts.AddCut("82410113","00200009327000008250400000","2444451041012200000","0163103100000010"); // min number 2 cells
    cuts.AddCut("82410113","00200009327000008250400000","2444451041014200000","0163103100000010"); // min number 4 cells
  } else if(trainConfig == 347){ // MESON
    cuts.AddCut("82410113","00200009327000008250400000","2444451041013200000","0163403100000010"); // rapidity variation  y<0.5
    cuts.AddCut("82410113","00200009327000008250400000","2444451041013200000","0163803100000010"); // rapidity variation  y<0.25
    cuts.AddCut("82410113","00200009327000008250400000","2444451041013200000","0163106100000010"); // alpha meson variation 1 0<alpha<0.8
    cuts.AddCut("82410113","00200009327000008250400000","2444451041013200000","0163105100000010"); // alpha meson variation 2 0<alpha<0.75
  } else if(trainConfig == 348){ // fourth set of variations
    cuts.AddCut("82410113","00200009327000008250400000","2444451044013200000","0163103100000010"); // tm variation
    cuts.AddCut("82410113","00200009327000008250400000","2444451045013200000","0163103100000010"); // tm variation
    cuts.AddCut("82410113","00200009327000008250400000","2444451046013200000","0163103100000010"); // tm variation
    cuts.AddCut("82410113","00200009327000008250400000","2444451041013200000","0163103100000000"); // min opening angle 0    -> open
    cuts.AddCut("82410113","00200009327000008250400000","2444451041013200000","0163103100000030"); // min opening angle 0.01 -> 2 cell diag
    //40 - 60
  } else if(trainConfig == 350){//dEdx e-line variation and dE/dx pi-line variation
    cuts.AddCut("84610113","00200009127000008250400000","2444451041013200000","0163103100000010"); //-5 < sigma < 5
    cuts.AddCut("84610113","00200009227000008250400000","2444451041013200000","0163103100000010"); //-3 < sigma < 5
    cuts.AddCut("84610113","00200009327400008250400000","2444451041013200000","0163103100000010"); //1, -10, 0.4, 3
    cuts.AddCut("84610113","00200009367400008250400000","2444451041013200000","0163103100000010"); //2, 0.5, 0.4, 3
  } else if(trainConfig == 351){//dE/dx pi-line variation and single pT variation and 2D triangular chi2 and psi pair
    cuts.AddCut("84610113","00200009317400008250400000","2444451041013200000","0163103100000010"); //0, -10, 0.4, 3
    cuts.AddCut("84610113","00200049327000008250400000","2444451041013200000","0163103100000010"); //single pT 0.075 GeV/c
    cuts.AddCut("84610113","00200019327000008250400000","2444451041013200000","0163103100000010"); //single pT 0.1   GeV/c
    cuts.AddCut("84610113","00200009327000008850400000","2444451041013200000","0163103100000010"); //20 & 0.1
  } else if(trainConfig == 352){//2D triangular chi2 and psi pair
    cuts.AddCut("84610113","00200009327000008260400000","2444451041013200000","0163103100000010"); //30 & 0.05
    cuts.AddCut("84610113","00200009327000008860400000","2444451041013200000","0163103100000010"); //20 & 0.05
    cuts.AddCut("84610113","00200009327000008280400000","2444451041013200000","0163103100000010"); //30 & 0.2
    cuts.AddCut("84610113","00200009327000008880400000","2444451041013200000","0163103100000010"); //20 & 0.2
  } else if(trainConfig == 353){//min TPC clusters and //elliptic cut Qt and alpha
    cuts.AddCut("84610113","00200000327000008250400000","2444451041013200000","0163103100000010"); //0
    cuts.AddCut("84610113","00200008327000008250400000","2444451041013200000","0163103100000010"); //0.35
    cuts.AddCut("84610113","00200009327000009250400000","2444451041013200000","0163103100000010"); // qT 0.03 no quadratic
    cuts.AddCut("84610113","00200009327000003250400000","2444451041013200000","0163103100000010"); // qT 0.05 y  quadratic
  } else if(trainConfig == 354){// and min phi max phi and NL variations
    cuts.AddCut("84610113","00209909327000008250400000","2444451041013200000","0163103100000010"); //4.54 > phi > 5.59
    cuts.AddCut("84610113","00200009327000008250400000","2444452041013200000","0163103100000010"); // PHOS calo NL
    cuts.AddCut("84610113","00200009327000008250400000","2444401041013200000","0163103100000010"); // PHOS people NL

    //Cluster Cuts varation
  } else if(trainConfig == 355){ // first set of variations CLUSTER
    cuts.AddCut("84610113","00200009327000008250400000","2444451041073200000","0163103100000010"); // min energy 0.2 GeV
    cuts.AddCut("84610113","00200009327000008250400000","2444451041083200000","0163103100000010"); // min energy 0.4 GeV
    cuts.AddCut("84610113","00200009327000008250400000","2444451041023200000","0163103100000010"); // min energy 0.5 GeV
  } else if(trainConfig == 356){ // second set of variations CLUSTER
    cuts.AddCut("84610113","00200009327000008250400000","2444451041013270000","0163103100000010"); // min/max M02  0.1<M<1
    cuts.AddCut("84610113","00200009327000008250400000","2444451041013280000","0163103100000010"); // min/max M02  0.1<M<1
    cuts.AddCut("84610113","00200009327000008250400000","2444451041012200000","0163103100000010"); // min number 2 cells
    cuts.AddCut("84610113","00200009327000008250400000","2444451041014200000","0163103100000010"); // min number 4 cells
  } else if(trainConfig == 357){ // MESON
    cuts.AddCut("84610113","00200009327000008250400000","2444451041013200000","0163403100000010"); // rapidity variation  y<0.5
    cuts.AddCut("84610113","00200009327000008250400000","2444451041013200000","0163803100000010"); // rapidity variation  y<0.25
    cuts.AddCut("84610113","00200009327000008250400000","2444451041013200000","0163106100000010"); // alpha meson variation 1 0<alpha<0.8
    cuts.AddCut("84610113","00200009327000008250400000","2444451041013200000","0163105100000010"); // alpha meson variation 2 0<alpha<0.75
  } else if(trainConfig == 358){ // fourth set of variations
    cuts.AddCut("84610113","00200009327000008250400000","2444451044013200000","0163103100000010"); // tm variation
    cuts.AddCut("84610113","00200009327000008250400000","2444451045013200000","0163103100000010"); // tm variation
    cuts.AddCut("84610113","00200009327000008250400000","2444451046013200000","0163103100000010"); // tm variation
    cuts.AddCut("84610113","00200009327000008250400000","2444451041013200000","0163103100000000"); // min opening angle 0    -> open
    cuts.AddCut("84610113","00200009327000008250400000","2444451041013200000","0163103100000030"); // min opening angle 0.01 -> 2 cell diag
    //60 - 100
  } else if(trainConfig == 360){//dEdx e-line variation and dE/dx pi-line variation
    cuts.AddCut("86010113","00200009127000008250400000","2444451041013200000","0163103100000010"); //-5 < sigma < 5
    cuts.AddCut("86010113","00200009227000008250400000","2444451041013200000","0163103100000010"); //-3 < sigma < 5
    cuts.AddCut("86010113","00200009327400008250400000","2444451041013200000","0163103100000010"); //1, -10, 0.4, 3
    cuts.AddCut("86010113","00200009367400008250400000","2444451041013200000","0163103100000010"); //2, 0.5, 0.4, 3
  } else if(trainConfig == 361){//dE/dx pi-line variation and single pT variation and 2D triangular chi2 and psi pair
    cuts.AddCut("86010113","00200009317400008250400000","2444451041013200000","0163103100000010"); //0, -10, 0.4, 3
    cuts.AddCut("86010113","00200049327000008250400000","2444451041013200000","0163103100000010"); //single pT 0.075 GeV/c
    cuts.AddCut("86010113","00200019327000008250400000","2444451041013200000","0163103100000010"); //single pT 0.1   GeV/c
    cuts.AddCut("86010113","00200009327000008850400000","2444451041013200000","0163103100000010"); //20 & 0.1
  } else if(trainConfig == 362){//2D triangular chi2 and psi pair
    cuts.AddCut("86010113","00200009327000008260400000","2444451041013200000","0163103100000010"); //30 & 0.05
    cuts.AddCut("86010113","00200009327000008860400000","2444451041013200000","0163103100000010"); //20 & 0.05
    cuts.AddCut("86010113","00200009327000008280400000","2444451041013200000","0163103100000010"); //30 & 0.2
    cuts.AddCut("86010113","00200009327000008880400000","2444451041013200000","0163103100000010"); //20 & 0.2
  } else if(trainConfig == 363){//min TPC clusters and //elliptic cut Qt and alpha
    cuts.AddCut("86010113","00200000327000008250400000","2444451041013200000","0163103100000010"); //0
    cuts.AddCut("86010113","00200008327000008250400000","2444451041013200000","0163103100000010"); //0.35
    cuts.AddCut("86010113","00200009327000009250400000","2444451041013200000","0163103100000010"); // qT 0.03 no quadratic
    cuts.AddCut("86010113","00200009327000003250400000","2444451041013200000","0163103100000010"); // qT 0.05 y  quadratic
  } else if(trainConfig == 364){// and min phi max phi and NL variations
    cuts.AddCut("86010113","00209909327000008250400000","2444451041013200000","0163103100000010"); //4.54 > phi > 5.59
    cuts.AddCut("86010113","00200009327000008250400000","2444452041013200000","0163103100000010"); // PHOS calo NL
    cuts.AddCut("86010113","00200009327000008250400000","2444401041013200000","0163103100000010"); // PHOS people NL

    //Cluster Cuts varation
  } else if(trainConfig == 365){ // first set of variations CLUSTER
    cuts.AddCut("86010113","00200009327000008250400000","2444451041073200000","0163103100000010"); // min energy 0.2 GeV
    cuts.AddCut("86010113","00200009327000008250400000","2444451041083200000","0163103100000010"); // min energy 0.4 GeV
    cuts.AddCut("86010113","00200009327000008250400000","2444451041023200000","0163103100000010"); // min energy 0.5 GeV
  } else if(trainConfig == 366){ // second set of variations CLUSTER
    cuts.AddCut("86010113","00200009327000008250400000","2444451041013270000","0163103100000010"); // min/max M02  0.1<M<1
    cuts.AddCut("86010113","00200009327000008250400000","2444451041013280000","0163103100000010"); // min/max M02  0.1<M<1
    cuts.AddCut("86010113","00200009327000008250400000","2444451041012200000","0163103100000010"); // min number 2 cells
    cuts.AddCut("86010113","00200009327000008250400000","2444451041014200000","0163103100000010"); // min number 4 cells
  } else if(trainConfig == 367){ // MESON
    cuts.AddCut("86010113","00200009327000008250400000","2444451041013200000","0163403100000010"); // rapidity variation  y<0.5
    cuts.AddCut("86010113","00200009327000008250400000","2444451041013200000","0163803100000010"); // rapidity variation  y<0.25
    cuts.AddCut("86010113","00200009327000008250400000","2444451041013200000","0163106100000010"); // alpha meson variation 1 0<alpha<0.8
    cuts.AddCut("86010113","00200009327000008250400000","2444451041013200000","0163105100000010"); // alpha meson variation 2 0<alpha<0.75
  } else if(trainConfig == 368){ // fourth set of variations
    cuts.AddCut("86010113","00200009327000008250400000","2444451044013200000","0163103100000010"); // tm variation
    cuts.AddCut("86010113","00200009327000008250400000","2444451045013200000","0163103100000010"); // tm variation
    cuts.AddCut("86010113","00200009327000008250400000","2444451046013200000","0163103100000010"); // tm variation
    cuts.AddCut("86010113","00200009327000008250400000","2444451041013200000","0163103100000000"); // min opening angle 0    -> open
    cuts.AddCut("86010113","00200009327000008250400000","2444451041013200000","0163103100000030"); // min opening angle 0.01 -> 2 cell diag


  // ===============================================================================================
  // Run 2 data EMC clusters pPb 8TeV
  // ===============================================================================================
  } else if (trainConfig == 400){ // EMCAL clusters standard cuts, triggers, no nonlin
    cuts.AddCut("80010113","00200009327000008250400000","1111100017032230000","0163103100000010"); // INT7
    cuts.AddCut("80010113","00200009327000008250400000","1111100057032230000","0163103100000010"); // INT7
  } else if (trainConfig == 401){ // EMCAL clusters standard cuts, triggers, no nonlin
    cuts.AddCut("80052113","00200009327000008250400000","1111100057032230000","0163103100000010"); // EMC7
    cuts.AddCut("80085113","00200009327000008250400000","1111100057032230000","0163103100000010"); // EG2
    cuts.AddCut("80083113","00200009327000008250400000","1111100057032230000","0163103100000010"); // EG1
  } else if (trainConfig == 402){ // EMCAL clusters standard cuts, triggers, no nonlin, open timing
    cuts.AddCut("80052113","00200009327000008250400000","1111100017032230000","0163103100000010"); // EMC7
    cuts.AddCut("80085113","00200009327000008250400000","1111100017032230000","0163103100000010"); // EG2
    cuts.AddCut("80083113","00200009327000008250400000","1111100017032230000","0163103100000010"); // EG1
    
  } else if ( trainConfig == 480){ // EMCAL standard cut but CALO+CALOFAST - NO TM
    cuts.AddCut("800a0113","00200009327000008250400000","1111100050032230000","0163103100000010"); // std INT7
    cuts.AddCut("800a1113","00200009327000008250400000","1111100050032230000","0163103100000010"); // std EMC7
    cuts.AddCut("800a2113","00200009327000008250400000","1111100050032230000","0163103100000010"); // std EG2
    cuts.AddCut("800a3113","00200009327000008250400000","1111100050032230000","0163103100000010"); // std EG1
  } else if ( trainConfig == 481){ // EMCAL standard cut but standard readout - NO TM
    cuts.AddCut("80010113","00200009327000008250400000","1111100050032230000","0163103100000010"); // std INT7
    cuts.AddCut("80052113","00200009327000008250400000","1111100050032230000","0163103100000010"); // std EMC7
    cuts.AddCut("80085113","00200009327000008250400000","1111100050032230000","0163103100000010"); // std EG2
    cuts.AddCut("80083113","00200009327000008250400000","1111100050032230000","0163103100000010"); // std EG1
  
  // ===============================================================================================
  // Run 2 data PHOS clusters 5TeV
  // ===============================================================================================
  // INT7 triggers
  } else if (trainConfig == 500) {  // PHOS  INT7
    cuts.AddCut("80010113","00200009327000008250400000","2446600051013200000","0163103100000010"); // PHOS group standard +-50ns
  } else if (trainConfig == 501) {  // PHOS  INT7
    cuts.AddCut("80010113","00200009327000008250400000","2446600041013200000","0163103100000010"); // PHOS group standard
    cuts.AddCut("80010113","00200009327000008250400000","2446600011013200000","0163103100000010"); // PHOS group standard, 1000 \mus
    cuts.AddCut("80010113","00200009327000008250400000","2446600061013200000","0163103100000010"); // PHOS group standard, -30, 50ns
    cuts.AddCut("80010113","00200009327000008250400000","24466000a1013200000","0163103100000010"); // PHOS group standard, -12.5, 13ns

  } else if (trainConfig == 502) {  // PHOS  INT7 with cents
    cuts.AddCut("80110113","00200009327000008250400000","2446600051013200000","0163103100000010"); // non lin 0-10%
    cuts.AddCut("81210113","00200009327000008250400000","2446600051013200000","0163103100000010"); // non lin 10-20%
    cuts.AddCut("82410113","00200009327000008250400000","2446600051013200000","0163103100000010"); // non lin 20-40%
    cuts.AddCut("84610113","00200009327000008250400000","2446600051013200000","0163103100000010"); // non lin 40-60%
    cuts.AddCut("86810113","00200009327000008250400000","2446600051013200000","0163103100000010"); // non lin 60-80%
    cuts.AddCut("88010113","00200009327000008250400000","2446600051013200000","0163103100000010"); // non lin 80-100%
  } else if (trainConfig == 503) {  // PHOS  INT7 with cents
    cuts.AddCut("80210113","00200009327000008250400000","2446600051013200000","0163103100000010"); // non lin 0-10%
    cuts.AddCut("86010113","00200009327000008250400000","2446600051013200000","0163103100000010"); // non lin 10-20%
    cuts.AddCut("a0110113","00200009327000008250400000","2446600051013200000","0163103100000010"); // non lin 0-5%
    cuts.AddCut("a1210113","00200009327000008250400000","2446600051013200000","0163103100000010"); // non lin 5-10%
  } else if (trainConfig == 504) {  // PHOS  INT7
    cuts.AddCut("80010113","00200009327000008250400000","2446600051013200000","0163103100000010"); // non lin none
    cuts.AddCut("80010113","00200009327000008250400000","2446641051013200000","0163103100000010"); // non lin 0-100%
    cuts.AddCut("80010113","00200009327000008250400000","2446651051013200000","0163103100000010"); // non lin 0-100%
    cuts.AddCut("80010113","00200009327000008250400000","2446601051013200000","0163103100000010"); // non lin PHOS 0-100%
  } else if (trainConfig == 505) {  // PHOS  INT7 with cents
    cuts.AddCut("80110113","00200009327000008250400000","2446601051013200000","0163103100000010"); // non lin 0-10%
    cuts.AddCut("81210113","00200009327000008250400000","2446601051013200000","0163103100000010"); // non lin 10-20%
    cuts.AddCut("82410113","00200009327000008250400000","2446601051013200000","0163103100000010"); // non lin 20-40%
    cuts.AddCut("84610113","00200009327000008250400000","2446601051013200000","0163103100000010"); // non lin 40-60%
    cuts.AddCut("86810113","00200009327000008250400000","2446601051013200000","0163103100000010"); // non lin 60-80%
    cuts.AddCut("88010113","00200009327000008250400000","2446601051013200000","0163103100000010"); // non lin 80-100%
  } else if (trainConfig == 506) {  // PHOS  INT7 with cents
    cuts.AddCut("80210113","00200009327000008250400000","2446601051013200000","0163103100000010"); // non lin 0-10%
    cuts.AddCut("86010113","00200009327000008250400000","2446601051013200000","0163103100000010"); // non lin 10-20%
    cuts.AddCut("a0110113","00200009327000008250400000","2446601051013200000","0163103100000010"); // non lin 0-5%
    cuts.AddCut("a1210113","00200009327000008250400000","2446601051013200000","0163103100000010"); // non lin 5-10%
  } else if (trainConfig == 507) {  // PHOS  INT7 with narrow cents
    cuts.AddCut("c0110113","00200009327000008250400000","2446641051013200000","0163103100000010"); // non lin 0-1%
    cuts.AddCut("c0210113","00200009327000008250400000","2446641051013200000","0163103100000010"); // non lin 0-2%
  } else if (trainConfig == 508) {  // AOD validation
    cuts.AddCut("80010113","00200009327000008250400000","2446641051013200000","0163103100000010"); // non lin 0-100%
  } else if (trainConfig == 509) {  // AOD validation PHOS NL
    cuts.AddCut("80010113","00200009327000008250400000","2446601051013200000","0163103100000010"); // non lin 0-100%
  } else if (trainConfig == 510) {  // AOD validation no NL
    cuts.AddCut("80010113","00200009327000008250400000","2446600051013200000","0163103100000010"); // non lin 0-100%
  } else if (trainConfig == 511) {  // JJ AOD validation no NL
    cuts.AddCut("80010113","00200009327000008250400000","2446600051013200000","0163103100000010"); // no non lin 0-100%
    cuts.AddCut("80010123","00200009327000008250400000","2446600051013200000","0163103100000010"); // no non lin 0-100%

  // ===============================================================================================
  // Run 2 data PHOS clusters 8TeV
  // ===============================================================================================
  } else if (trainConfig == 600) {  // PHOS PHI7
    cuts.AddCut("80001113","00200009327000008250400000","2446600051013200000","0163103100000010");// PHI7 triggers
    cuts.AddCut("80001113","00200009327000008250400000","2446600011013200000","0163103100000010");// PHI7 triggers
  } else if (trainConfig == 601) {  // PHOS PHI7
    cuts.AddCut("80062113","00200009327000008250400000","2446600051013200000","0163103100000010");// PHI7 triggers
    cuts.AddCut("80062113","00200009327000008250400000","2446600011013200000","0163103100000010");// PHI7 triggers


  // ===============================================================================================
  // Run 2 data DMC clusters pPb 5TeV
  // ===============================================================================================
  } else if (trainConfig == 700){ // DCal clusters standard cuts, triggers, no nonlin, open timing
    cuts.AddCut("80010113","00200009327000008250400000","3885500057032230000","0163103100000010"); // INT7
    cuts.AddCut("80010113","00200009327000008250400000","3885500017032230000","0163103100000010"); // INT7
  } else if (trainConfig == 701){  // DCal clusters standard cuts, triggers, no nonlin
    cuts.AddCut("00055113","00200009327000008250400000","3885500057032230000","0163103100000010"); // EMC7
    cuts.AddCut("00089113","00200009327000008250400000","3885500057032230000","0163103100000010"); // EG2
    cuts.AddCut("0008b113","00200009327000008250400000","3885500057032230000","0163103100000010"); // EG1
  } else if (trainConfig == 702){  // DCal clusters standard cuts, triggers, no nonlin, open timing
    cuts.AddCut("00055113","00200009327000008250400000","3885500017032230000","0163103100000010"); // EMC7
    cuts.AddCut("00089113","00200009327000008250400000","3885500017032230000","0163103100000010"); // EG2
    cuts.AddCut("0008b113","00200009327000008250400000","3885500017032230000","0163103100000010"); // EG1
  } else if (trainConfig == 703){ // EMCAL clusters standard cuts, triggers, no nonlin, open timing
    cuts.AddCut("80110113","00200009327000008250400000","3885500057032230000","0163103100000010"); // 0-10
    cuts.AddCut("81210113","00200009327000008250400000","3885500057032230000","0163103100000010"); // 10-20
    cuts.AddCut("82410113","00200009327000008250400000","3885500057032230000","0163103100000010"); // 20-40
    cuts.AddCut("84610113","00200009327000008250400000","3885500057032230000","0163103100000010"); // 40-60
    cuts.AddCut("86810113","00200009327000008250400000","3885500057032230000","0163103100000010"); // 60-80
    cuts.AddCut("88010113","00200009327000008250400000","3885500057032230000","0163103100000010"); // 80-100
  } else if (trainConfig == 704){ // EMCAL clusters standard cuts, triggers, no nonlin, open timing
    cuts.AddCut("80210113","00200009327000008250400000","3885500057032230000","0163103100000010"); // 0-20
    cuts.AddCut("86010113","00200009327000008250400000","3885500057032230000","0163103100000010"); // 60-100
    cuts.AddCut("a0110113","00200009327000008250400000","3885500057032230000","0163103100000010"); // 0-5
    cuts.AddCut("a1210113","00200009327000008250400000","3885500057032230000","0163103100000010"); // 5-10
  } else if (trainConfig == 705){ // EMCAL clusters standard cuts, triggers, no nonlin, open timing
    cuts.AddCut("80110113","00200009327000008250400000","3885500017032230000","0163103100000010"); // 0-10
    cuts.AddCut("81210113","00200009327000008250400000","3885500017032230000","0163103100000010"); // 10-20
    cuts.AddCut("82410113","00200009327000008250400000","3885500017032230000","0163103100000010"); // 20-40
    cuts.AddCut("84610113","00200009327000008250400000","3885500017032230000","0163103100000010"); // 40-60
    cuts.AddCut("86810113","00200009327000008250400000","3885500017032230000","0163103100000010"); // 60-80
    cuts.AddCut("88010113","00200009327000008250400000","3885500017032230000","0163103100000010"); // 80-100
  } else if (trainConfig == 706){ // EMCAL clusters standard cuts, triggers, no nonlin, open timing
    cuts.AddCut("80210113","00200009327000008250400000","3885500017032230000","0163103100000010"); // 0-20
    cuts.AddCut("86010113","00200009327000008250400000","3885500017032230000","0163103100000010"); // 60-100
    cuts.AddCut("a0110113","00200009327000008250400000","3885500017032230000","0163103100000010"); // 0-5
    cuts.AddCut("a1210113","00200009327000008250400000","3885500017032230000","0163103100000010"); // 5-10
  } else if (trainConfig == 707){ // AOD validation
    cuts.AddCut("80010113","00200009327000008250400000","3885500057032230000","0163103100000010"); // INT7
    
  } else if (trainConfig == 780){ // DCal clusters standard cuts, CALO+CALOFAST triggers, no nonlin, open timing, no TM
    cuts.AddCut("800a0113","00200009327000008250400000","3885500050032230000","0163103100000010"); // INT7
    cuts.AddCut("000a6113","00200009327000008250400000","3885500050032230000","0163103100000010"); // EMC7
    cuts.AddCut("000a7113","00200009327000008250400000","3885500050032230000","0163103100000010"); // EG2
    cuts.AddCut("000a8113","00200009327000008250400000","3885500050032230000","0163103100000010"); // EG1
  } else if (trainConfig == 781){ // DCal clusters standard cuts, triggers, no nonlin, open timing, no TM
    cuts.AddCut("80010113","00200009327000008250400000","3885500050032230000","0163103100000010"); // INT7
    cuts.AddCut("00055113","00200009327000008250400000","3885500050032230000","0163103100000010"); // EMC7
    cuts.AddCut("00089113","00200009327000008250400000","3885500050032230000","0163103100000010"); // EG2
    cuts.AddCut("0008b113","00200009327000008250400000","3885500050032230000","0163103100000010"); // EG1
  // ===============================================================================================
  // Run 2 data DMC clusters pPb 8TeV
  // ===============================================================================================
  } else if (trainConfig == 800){  // DCal clusters standard cuts, triggers, no nonlin
    cuts.AddCut("80010113","00200009327000008250400000","3885500057032230000","0163103100000010"); // INT7
    cuts.AddCut("80010113","00200009327000008250400000","3885500017032230000","0163103100000010"); // INT7
  } else if (trainConfig == 801){  // DCal clusters standard cuts, triggers, no nonlin
    cuts.AddCut("00055113","00200009327000008250400000","3885500057032230000","0163103100000010"); // EMC7
    cuts.AddCut("00089113","00200009327000008250400000","3885500057032230000","0163103100000010"); // EG2
    cuts.AddCut("0008b113","00200009327000008250400000","3885500057032230000","0163103100000010"); // EG1
  } else if (trainConfig == 802){  // DCal clusters standard cuts, triggers, no nonlin, open timing
    cuts.AddCut("80052113","00200009327000008250400000","3885500017032230000","0163103100000010"); // EMC7
    cuts.AddCut("00089113","00200009327000008250400000","3885500017032230000","0163103100000010"); // EG2
    cuts.AddCut("0008b113","00200009327000008250400000","3885500017032230000","0163103100000010"); // EG1

  } else {
    Error(Form("GammaConvCalo_%i",trainConfig), "wrong trainConfig variable no cuts have been specified for the configuration");
    return;
  }

  if(!cuts.AreValid()){
    cout << "\n\n****************************************************" << endl;
    cout << "ERROR: No valid cuts stored in CutHandlerConvCalo! Returning..." << endl;
    cout << "****************************************************\n\n" << endl;
    return;
  }

  TList *EventCutList   = new TList();
  TList *ConvCutList    = new TList();
  TList *ClusterCutList = new TList();
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

  EventCutList->SetOwner(kTRUE);
  AliConvEventCuts **analysisEventCuts        = new AliConvEventCuts*[numberOfCuts];
  ConvCutList->SetOwner(kTRUE);
  AliConversionPhotonCuts **analysisCuts      = new AliConversionPhotonCuts*[numberOfCuts];
  ClusterCutList->SetOwner(kTRUE);
  AliCaloPhotonCuts **analysisClusterCuts     = new AliCaloPhotonCuts*[numberOfCuts];
  MesonCutList->SetOwner(kTRUE);
  AliConversionMesonCuts **analysisMesonCuts  = new AliConversionMesonCuts*[numberOfCuts];

  for(Int_t i = 0; i<numberOfCuts; i++){
    //create AliCaloTrackMatcher instance, if there is none present
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

    analysisClusterCuts[i] = new AliCaloPhotonCuts(isMC);
    analysisClusterCuts[i]->SetHistoToModifyAcceptance(histoAcc);
    analysisClusterCuts[i]->SetV0ReaderName(V0ReaderName);
    analysisClusterCuts[i]->SetCorrectionTaskSetting(corrTaskSetting);
    analysisClusterCuts[i]->SetCaloTrackMatcherName(TrackMatcherName);
    if (runLightOutput > 0) analysisClusterCuts[i]->SetLightOutput(kTRUE);
    analysisClusterCuts[i]->InitializeCutsFromCutString((cuts.GetClusterCut(i)).Data());
    ClusterCutList->Add(analysisClusterCuts[i]);
    analysisClusterCuts[i]->SetExtendedMatchAndQA(enableExtMatchAndQA);
    analysisClusterCuts[i]->SetFillCutHistograms("");

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
  task->SetCaloCutList(numberOfCuts,ClusterCutList);
  task->SetMesonCutList(numberOfCuts,MesonCutList);
  task->SetMoveParticleAccordingToVertex(kTRUE);
  task->SetCorrectionTaskSetting(corrTaskSetting);
  task->SetDoMesonAnalysis(kTRUE);
  task->SetDoMesonQA(enableQAMesonTask); //Attention new switch for Pi0 QA
  task->SetDoPhotonQA(enableQAPhotonTask);  //Attention new switch small for Photon QA
  task->SetDoClusterQA(1);  //Attention new switch small for Cluster QA
  task->SetUseTHnSparse(isUsingTHnSparse);
  task->SetDoTreeInvMassShowerShape(doTreeClusterShowerShape);
  task->SetEnableSortingOfMCClusLabels(enableSortingMCLabels);
  if(enableExtMatchAndQA > 1){ task->SetPlotHistsExtQA(kTRUE);}

  //connect containers
  AliAnalysisDataContainer *coutput =
    mgr->CreateContainer(!(corrTaskSetting.CompareTo("")) ? Form("GammaConvCalo_%i",trainConfig) : Form("GammaConvCalo_%i_%s",trainConfig,corrTaskSetting.Data()), TList::Class(),
              AliAnalysisManager::kOutputContainer,Form("GammaConvCalo_%i.root",trainConfig) );

  mgr->AddTask(task);
  mgr->ConnectInput(task,0,cinput);
  mgr->ConnectOutput(task,1,coutput);

  return;

}
