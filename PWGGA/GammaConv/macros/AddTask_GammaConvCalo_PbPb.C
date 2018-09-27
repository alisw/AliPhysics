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
//PbPb together with all supporting classes
//***************************************************************************************

//***************************************************************************************
//CutHandler contains all cuts for a certain analysis and trainconfig,
//it automatically checks length of cutStrings and takes care of the number of added cuts,
//no specification of the variable 'numberOfCuts' needed anymore.
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
void AddTask_GammaConvCalo_PbPb(  Int_t     trainConfig                     = 1,                      // change different set of cuts
                                  Int_t     isMC                            = 0,                      // run MC
                                  Int_t     enableQAMesonTask               = 0,                      // enable QA in AliAnalysisTaskGammaConvV1
                                  Int_t     enableQAPhotonTask              = 0,                      // enable additional QA task
                                  TString   fileNameInputForWeighting       = "MCSpectraInput.root",  // path to file for weigting input / modified acceptance
                                  Int_t     headerSelectionInt              = 0,                      // 1 pi0 header, 2 eta header, 3 both (only for "named" boxes)
                                  Bool_t    enableHeaderOverlap             = kTRUE,                // allow overlapping header for the clusters
                                  TString   cutnumberAODBranch              = "100000006008400000001500000",
                                  TString   periodName                      = "LHC13d2",              // name of the period for added signals and weighting
                                  Bool_t    doWeighting                     = kFALSE,                 // enable Weighting
                                  Int_t     enableExtMatchAndQA             = 0,                      // disabled (0), extMatch (1), extQA_noCellQA (2), extMatch+extQA_noCellQA (3), extQA+cellQA (4), extMatch+extQA+cellQA (5)
                                  Bool_t    isUsingTHnSparse                = kTRUE,                  // enable or disable usage of THnSparses for background estimation
                                  Bool_t    enableV0findingEffi             = kFALSE,                 // enables V0finding efficiency histograms
                                  TString   periodNameV0Reader              = "",                     // period Name for V0Reader
                                  Bool_t    enableSortingMCLabels           = kTRUE,                  // enable sorting for MC cluster labels
                                  Int_t     runLightOutput                  = 0,                      // switch to run light output 0 (disabled), 1 (for CutClasses), 2 (for cutClasses and task)
                                  Bool_t    doFlattening                    = kFALSE,                 // switch on centrality flattening for LHC11h
                                  TString   fileNameInputForCentFlattening  = "",                     // file name for centrality flattening
                                  Bool_t    doPrimaryTrackMatching          = kTRUE,                  // enable basic track matching for all primary tracks to cluster
                                  Bool_t    doMultiplicityWeighting         = kFALSE,                         //
                                  TString   fileNameInputForMultWeighing    = "Multiplicity.root",            //
                                  TString   periodNameAnchor                = "",
                                  TString   additionalTrainConfig           = "0"                     // additional counter for trainconfig
                                ) {

  Int_t trackMatcherRunningMode = 0; // CaloTrackMatcher running mode
  Bool_t doTreeClusterShowerShape = kFALSE; // enable tree for meson cand EMCal shower shape studies
  TH1S* histoAcc = 0x0;                     // histo for modified acceptance
  TString corrTaskSetting = ""; // select which correction task setting to use
  //parse additionalTrainConfig flag
  TObjArray *rAddConfigArr = additionalTrainConfig.Tokenize("_");
  if(rAddConfigArr->GetEntries()<1){cout << "ERROR: AddTask_GammaConvCalo_PbPb during parsing of additionalTrainConfig String '" << additionalTrainConfig.Data() << "'" << endl; return;}
  TObjString* rAdditionalTrainConfig;
  for(Int_t i = 0; i<rAddConfigArr->GetEntries() ; i++){
    if(i==0) rAdditionalTrainConfig = (TObjString*)rAddConfigArr->At(i);
    else{
      TObjString* temp = (TObjString*) rAddConfigArr->At(i);
      TString tempStr = temp->GetString();
      if(tempStr.CompareTo("INVMASSCLUSTree") == 0){
        cout << "INFO: AddTask_GammaConvCalo_PbPb activating 'INVMASSCLUSTree'" << endl;
        doTreeClusterShowerShape = kTRUE;
      }else if(tempStr.BeginsWith("MODIFYACC")){
        cout << "INFO: AddTask_GammaConvCalo_PbPb activating 'MODIFYACC'" << endl;
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
        cout << "INFO: AddTask_GammaConvCalo_PbPb will use custom branch from Correction Framework!" << endl;
        corrTaskSetting = tempStr;
        corrTaskSetting.Replace(0,2,"");
      }else if(tempStr.BeginsWith("TM")){
        TString tempType = tempStr;
        tempType.Replace(0,2,"");
        trackMatcherRunningMode = tempType.Atoi();
        cout << Form("INFO: AddTask_GammaConvCalo_PbPb will use running mode '%i' for the TrackMatcher!",trackMatcherRunningMode) << endl;
      }
    }
  }
  TString sAdditionalTrainConfig = rAdditionalTrainConfig->GetString();
  if (sAdditionalTrainConfig.Atoi() > 0){
    trainConfig = trainConfig + sAdditionalTrainConfig.Atoi();
    cout << "INFO: AddTask_GammaConvCalo_PbPb running additionalTrainConfig '" << sAdditionalTrainConfig.Atoi() << "', train config: '" << trainConfig << "'" << endl;
  }
  cout << "corrTaskSetting: " << corrTaskSetting.Data() << endl;

  Int_t isHeavyIon = 1;

  // ================== GetAnalysisManager ===============================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error(Form("AddTask_GammaConvCalo_PbPb_%i",trainConfig), "No analysis manager found.");
    return ;
  }

  // ================== GetInputEventHandler =============================
  AliVEventHandler *inputHandler = mgr->GetInputEventHandler();

  Bool_t isMCForOtherTasks = kFALSE;
  if (isMC > 0) isMCForOtherTasks = kTRUE;

  //========= Add PID Reponse to ANALYSIS manager ====
  if(!(AliPIDResponse*)mgr->GetTask("PIDResponseTask")){
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
    AddTaskPIDResponse(isMCForOtherTasks);
  }

  //=========  Set Cutnumber for V0Reader ================================
  TString cutnumberPhoton   = "00000008400100001500000000";
  if (periodNameV0Reader.CompareTo("LHC17n") == 0 || periodNameV0Reader.Contains("LHC17j7"))
    cutnumberPhoton         = "00000088400100001500000000";
  TString cutnumberEvent    = "10000003";
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
  task->SetTrackMatcherRunningMode(trackMatcherRunningMode);

  //create cut handler
  CutHandlerConvCalo cuts;

  // cluster cuts
  // 0 "ClusterType",  1 "EtaMin", 2 "EtaMax", 3 "PhiMin", 4 "PhiMax", 5 "DistanceToBadChannel", 6 "Timing", 7 "TrackMatching", 8 "ExoticCell",
  // 9 "MinEnergy", 10 "MinNCells", 11 "MinM02", 12 "MaxM02", 13 "MinMaxM20", 14 "RecConv", 15 "MaximumDispersion", 16 "NLM"

  //****************************************************************************************************
  // EMCal 2.76 TeV Pb-Pb  LHC11h & LHC10h MB
  //****************************************************************************************************
  if (trainConfig == 1){ // EMCAL clusters
    cuts.AddCut("60100013","00200009297002008250400000","1111100053032230000","0163103100000010"); // 0-5%
    cuts.AddCut("61200013","00200009297002008250400000","1111100053032230000","0163103100000010"); // 5-10%
    cuts.AddCut("50100013","00200009297002008250400000","1111100053032230000","0163103100000010"); // 0-10%
    cuts.AddCut("52400013","00200009297002008250400000","1111100053032230000","0163103100000010"); // 20-40%
    cuts.AddCut("52500013","00200009297002008250400000","1111100053032230000","0163103100000010"); // 20-50%
  } else if (trainConfig == 2){ // EMCAL clusters no timing
    cuts.AddCut("60100013","00200009297002008250400000","1111100003032230000","0163103100000010"); // 0-5%
    cuts.AddCut("61200013","00200009297002008250400000","1111100003032230000","0163103100000010"); // 5-10%
    cuts.AddCut("50100013","00200009297002008250400000","1111100003032230000","0163103100000010"); // 0-10%
    cuts.AddCut("52400013","00200009297002008250400000","1111100003032230000","0163103100000010"); // 20-40%
    cuts.AddCut("52500013","00200009297002008250400000","1111100003032230000","0163103100000010"); // 20-50%
  } else if (trainConfig == 3){ // EMCAL clusters
    cuts.AddCut("60100013","00200009297002008250400000","1111100053032230000","0163103100000010"); // 0-5%
    cuts.AddCut("61200013","00200009297002008250400000","1111100053032230000","0163103100000010"); // 5-10%
    cuts.AddCut("50100013","00200009297002008250400000","1111100053032230000","0163103100000010"); // 0-10%
    cuts.AddCut("51200013","00200009297002008250400000","1111100053032230000","0163103100000010"); // 10-20%
    cuts.AddCut("52400013","00200009297002008250400000","1111100053032230000","0163103100000010"); // 20-40%
  } else if (trainConfig == 4){ // EMCAL clusters no timing
    cuts.AddCut("60100013","00200009297002008250400000","1111100003032230000","0163103100000010"); // 0-5%
    cuts.AddCut("61200013","00200009297002008250400000","1111100003032230000","0163103100000010"); // 5-10%
    cuts.AddCut("50100013","00200009297002008250400000","1111100003032230000","0163103100000010"); // 0-10%
    cuts.AddCut("51200013","00200009297002008250400000","1111100003032230000","0163103100000010"); // 10-20%
    cuts.AddCut("52400013","00200009297002008250400000","1111100003032230000","0163103100000010"); // 20-40%
  } else if (trainConfig == 5){ // EMCAL clusters
    cuts.AddCut("54600013","00200009297002008250400000","1111100053032230000","0163103100000010"); // 40-60%
    cuts.AddCut("56800013","00200009297002008250400000","1111100053032230000","0163103100000010"); // 60-80%
    cuts.AddCut("52600013","00200009297002008250400000","1111100053032230000","0163103100000010"); // 20-60%
    cuts.AddCut("54800013","00200009297002008250400000","1111100053032230000","0163103100000010"); // 40-80%
    cuts.AddCut("52500013","00200009297002008250400000","1111100053032230000","0163103100000010"); // 20-50%
  } else if (trainConfig == 6){ // EMCAL clusters no timing
    cuts.AddCut("54600013","00200009297002008250400000","1111100003032230000","0163103100000010"); // 40-60%
    cuts.AddCut("56800013","00200009297002008250400000","1111100003032230000","0163103100000010"); // 60-80%
    cuts.AddCut("52600013","00200009297002008250400000","1111100003032230000","0163103100000010"); // 20-60%
    cuts.AddCut("54800013","00200009297002008250400000","1111100003032230000","0163103100000010"); // 40-80%
    cuts.AddCut("52500013","00200009297002008250400000","1111100003032230000","0163103100000010"); // 20-50%

  } else if (trainConfig == 7){ // EMCAL clusters
    cuts.AddCut("60100013","00200009297002008250400000","1111100053032230000","0163103100000010"); // 0-5%
    cuts.AddCut("61200013","00200009297002008250400000","1111100053032230000","0163103100000010"); // 5-10%
    cuts.AddCut("50100013","00200009297002008250400000","1111100053032230000","0163103100000010"); // 0-10%
  } else if (trainConfig == 8){ // EMCAL clusters no timing
    cuts.AddCut("60100013","00200009297002008250400000","1111100003032230000","0163103100000010"); // 0-5%
    cuts.AddCut("61200013","00200009297002008250400000","1111100003032230000","0163103100000010"); // 5-10%
    cuts.AddCut("50100013","00200009297002008250400000","1111100003032230000","0163103100000010"); // 0-10%
  } else if (trainConfig == 9){ // EMCAL clusters
    cuts.AddCut("51200013","00200009297002008250400000","1111100053032230000","0163103100000010"); // 10-20%
    cuts.AddCut("52400013","00200009297002008250400000","1111100053032230000","0163103100000010"); // 20-40%
    cuts.AddCut("52500013","00200009297002008250400000","1111100053032230000","0163103100000010"); // 20-50%
  } else if (trainConfig == 10){ // EMCAL clusters no timing
    cuts.AddCut("51200013","00200009297002008250400000","1111100003032230000","0163103100000010"); // 10-20%
    cuts.AddCut("52400013","00200009297002008250400000","1111100003032230000","0163103100000010"); // 20-40%
    cuts.AddCut("52500013","00200009297002008250400000","1111100003032230000","0163103100000010"); // 20-50%
  } else if (trainConfig == 11){ // EMCAL clusters
    cuts.AddCut("54600013","00200009297002008250400000","1111100053032230000","0163103100000010"); // 40-60%
    cuts.AddCut("56800013","00200009297002008250400000","1111100053032230000","0163103100000010"); // 60-80%
  } else if (trainConfig == 12){ // EMCAL clusters no timing
    cuts.AddCut("54600013","00200009297002008250400000","1111100003032230000","0163103100000010"); // 40-60%
    cuts.AddCut("56800013","00200009297002008250400000","1111100003032230000","0163103100000010"); // 60-80%


  } else if (trainConfig == 13){ // EMCAL clusters
    cuts.AddCut("50900013","00200009297002008250400000","1111100053032230000","0163103100000010"); // 0-90%
  } else if (trainConfig == 14){ // EMCAL clusters
    cuts.AddCut("50900013","00200009297002008250400000","1111100003032230000","0163103100000010"); // 0-90%

  //****************************************************************************************************
  // EMCal 2.76TeV Pb-Pb for LHC11h EMC trigger
  //****************************************************************************************************
  } else if (trainConfig == 30){ // EMCAL clusters
    cuts.AddCut("50980013","00200009297002008250400000","1111100053032230000","0163103100000010"); // 0-90%
    cuts.AddCut("50180013","00200009297002008250400000","1111100053032230000","0163103100000010"); // 0-10%
    cuts.AddCut("51280013","00200009297002008250400000","1111100053032230000","0163103100000010"); // 10-20%
    cuts.AddCut("52580013","00200009297002008250400000","1111100053032230000","0163103100000010"); // 20-50%
    cuts.AddCut("55880013","00200009297002008250400000","1111100053032230000","0163103100000010"); // 50-80%
  } else if (trainConfig == 31){ // EMCAL clusters no timing cut
    cuts.AddCut("50980013","00200009297002008250400000","1111100003032230000","0163103100000010"); // 0-90%
    cuts.AddCut("50180013","00200009297002008250400000","1111100003032230000","0163103100000010"); // 0-10%
    cuts.AddCut("51280013","00200009297002008250400000","1111100003032230000","0163103100000010"); // 10-20%
    cuts.AddCut("52580013","00200009297002008250400000","1111100003032230000","0163103100000010"); // 20-50%
    cuts.AddCut("55880013","00200009297002008250400000","1111100003032230000","0163103100000010"); // 50-80%
  } else if (trainConfig == 32){ // EMCAL clusters
    cuts.AddCut("50980013","00200009297002008250400000","1111100053032230000","0163103100000010"); // 0-90%
  } else if (trainConfig == 33){ // EMCAL clusters no timing cut
    cuts.AddCut("50980013","00200009297002008250400000","1111100003032230000","0163103100000010"); // 0-90%
  } else if (trainConfig == 34){ // EMCAL clusters
    cuts.AddCut("50100013","00200009297002008250400000","1111102053032230000","0163103100000010"); // 0-10
    cuts.AddCut("50100013","00200009297002008250400000","1111171053032230000","0163103100000010"); // 0-10
  } else if (trainConfig == 35){ // EMCAL clusters
    cuts.AddCut("50100013","00200009297002008250400000","11111020530a2230000","0163103100000010"); // 0-10 // using cluster cuts of Astrid
    cuts.AddCut("50100013","00200009297002008250400000","11111710530a2230000","0163103100000010"); // 0-10
    cuts.AddCut("50100013","00200009297002008250400000","1111171053032230000","0163103100000010");
  } else if (trainConfig == 36){ // EMCAL clusters
    cuts.AddCut("52500013","00200009297002008250400000","1111102053032230000","0163103100000010"); // 20-50
    cuts.AddCut("52500013","00200009297002008250400000","1111172053032230000","0163103100000010"); // 20-50
  } else if (trainConfig == 37){ // EMCAL clusters
    cuts.AddCut("52500013","00200009297002008250400000","11111020530a2230000","0163103100000010"); // 20-50
    cuts.AddCut("52500013","00200009297002008250400000","11111720530a2230000","0163103100000010"); // 20-50
  } else if (trainConfig == 38){ // EMCAL clusters - added signals
    cuts.AddCut("50100023","00200009297002008250400000","1111171053032230000","0163103100000010"); // 0-10
    cuts.AddCut("50100023","00200009297002008250400000","11111710530a2230000","0163103100000010"); // 0-10
    cuts.AddCut("50100023","00200009297002008250400000","11111710530b2230000","0163103100000010"); // 0-10
  } else if (trainConfig == 39){ // EMCAL clusters - added signals
    cuts.AddCut("50100023","00200009297002008250400000","1111171053032230000","0163103100000010"); // 0-10
    cuts.AddCut("50100023","00200009297002008250400000","11111710530a2230000","0163103100000010"); // 0-10
    cuts.AddCut("50100023","00200009297002008250400000","11111710530b2230000","0163103100000010"); // 0-10
  } else if (trainConfig == 40){ // EMCAL clusters - added signals
    cuts.AddCut("52500023","00200009297002008250400000","1111172053032230000","0163103100000010"); // 20-50
    cuts.AddCut("52500023","00200009297002008250400000","11111720530a2230000","0163103100000010"); // 20-50
    cuts.AddCut("52500023","00200009297002008250400000","11111720530b2230000","0163103100000010"); // 20-50
  } else if (trainConfig == 41){ // EMCAL clusters - added signals
    cuts.AddCut("52500023","00200009297002008250400000","1111172053032230000","0163103100000010"); // 20-50
    cuts.AddCut("52500023","00200009297002008250400000","11111720530a2230000","0163103100000010"); // 20-50
    cuts.AddCut("52500023","00200009297002008250400000","11111720530b2230000","0163103100000010"); // 20-50
  } else if (trainConfig == 42){ // EMCAL clusters all headers
    cuts.AddCut("50100003","00200009297002008250400000","11111020530a2230000","0163103100000010"); // 0-10
    cuts.AddCut("50100003","00200009297002008250400000","11111710530a2230000","0163103100000010"); // 0-10
    cuts.AddCut("50100003","00200009297002008250400000","1111171053032230000","0163103100000010");
  } else if (trainConfig == 43){ // EMCAL clusters added signals pi0 header forseen
    cuts.AddCut("50100023","00200009297002008250400000","11111020530a2230000","0163103100000010"); // 0-10 // using cluster cuts of Astrid
    cuts.AddCut("50100023","00200009297002008250400000","11111710530a2230000","0163103100000010"); // 0-10
    cuts.AddCut("50100023","00200009297002008250400000","1111171053032230000","0163103100000010");
  } else if (trainConfig == 44){ // EMCAL clusters added signals eta header forseen
    cuts.AddCut("50100023","00200009297002008250400000","11111020530a2230000","0163103100000010"); // 0-10 // using cluster cuts of Astrid
    cuts.AddCut("50100023","00200009297002008250400000","11111710530a2230000","0163103100000010"); // 0-10
    cuts.AddCut("50100023","00200009297002008250400000","1111171053032230000","0163103100000010");
  } else if (trainConfig == 45){ // EMCAL clusters added signals other header forseen
    cuts.AddCut("50100023","00200009297002008250400000","11111020530a2230000","0163103100000010"); // 0-10 // using cluster cuts of Astrid
    cuts.AddCut("50100023","00200009297002008250400000","11111710530a2230000","0163103100000010"); // 0-10
    cuts.AddCut("50100023","00200009297002008250400000","1111171053032230000","0163103100000010");
  } else if (trainConfig == 46){ // EMCAL clusters added signals other header forseen
    cuts.AddCut("50100023","00200009297002008250400000","11111020530a2230000","0163103100000010"); // 0-10 // using cluster cuts of Astrid
    cuts.AddCut("50100023","00200009297002008250400000","11111710530a2230000","0163103100000010"); // 0-10
    cuts.AddCut("50100023","00200009297002008250400000","1111171053032230000","0163103100000010");
  } else if (trainConfig == 47){ // EMCAL clusters V0 Cent
    cuts.AddCut("10100013","00200009297002008250400000","11111020530a2230000","0163103100000010"); // 0-10 // using cluster cuts of Astrid
    cuts.AddCut("10100013","00200009297002008250400000","11111710530a2230000","0163103100000010"); // 0-10
    cuts.AddCut("10100013","00200009297002008250400000","1111171053032230000","0163103100000010");

  //****************************************************************************************************
  // PHOS 2.76TeV Pb-Pb LHC10h & LHC11h
  //****************************************************************************************************
  } else if (trainConfig == 101){ // PHOS clusters
    cuts.AddCut("60100013","00200009297002008250400000","2444400042033200000","0163103100000010"); // 0-5%
    cuts.AddCut("61200013","00200009297002008250400000","2444400042033200000","0163103100000010"); // 5-10%
    cuts.AddCut("50100013","00200009297002008250400000","2444400042033200000","0163103100000010"); // 0-10%
    cuts.AddCut("52400013","00200009297002008250400000","2444400042033200000","0163103100000010"); // 20-40%
    cuts.AddCut("52500013","00200009297002008250400000","2444400042033200000","0163103100000010"); // 20-50%
  } else if (trainConfig == 102){ // PHOS clusters
    cuts.AddCut("60100013","00200009297002008250400000","2444400042033200000","0163103100000010"); // 0-5%
    cuts.AddCut("61200013","00200009297002008250400000","2444400042033200000","0163103100000010"); // 5-10%
    cuts.AddCut("50100013","00200009297002008250400000","2444400042033200000","0163103100000010"); // 0-10%
    cuts.AddCut("51200013","00200009297002008250400000","2444400042033200000","0163103100000010"); // 10-20%
    cuts.AddCut("52400013","00200009297002008250400000","2444400042033200000","0163103100000010"); // 20-40%
  } else if (trainConfig == 103){ // PHOS clusters
    cuts.AddCut("54600013","00200009297002008250400000","2444400042033200000","0163103100000010"); // 40-60%
    cuts.AddCut("56800013","00200009297002008250400000","2444400042033200000","0163103100000010"); // 60-80%
    cuts.AddCut("52600013","00200009297002008250400000","2444400042033200000","0163103100000010"); // 20-60%
    cuts.AddCut("54800013","00200009297002008250400000","2444400042033200000","0163103100000010"); // 40-80%
    cuts.AddCut("52500013","00200009297002008250400000","2444400042033200000","0163103100000010"); // 20-50%
  } else if (trainConfig == 104){ // PHOS clusters no timing
    cuts.AddCut("60100013","00200009297002008250400000","2444400002033200000","0163103100000010"); // 0-5%
    cuts.AddCut("61200013","00200009297002008250400000","2444400002033200000","0163103100000010"); // 5-10%
    cuts.AddCut("50100013","00200009297002008250400000","2444400002033200000","0163103100000010"); // 0-10%
    cuts.AddCut("52400013","00200009297002008250400000","2444400002033200000","0163103100000010"); // 20-40%
    cuts.AddCut("52500013","00200009297002008250400000","2444400002033200000","0163103100000010"); // 20-50%
  } else if (trainConfig == 105){ // PHOS clusters no timing
    cuts.AddCut("60100013","00200009297002008250400000","2444400002033200000","0163103100000010"); // 0-5%
    cuts.AddCut("61200013","00200009297002008250400000","2444400002033200000","0163103100000010"); // 5-10%
    cuts.AddCut("50100013","00200009297002008250400000","2444400002033200000","0163103100000010"); // 0-10%
    cuts.AddCut("51200013","00200009297002008250400000","2444400002033200000","0163103100000010"); // 10-20%
    cuts.AddCut("52400013","00200009297002008250400000","2444400002033200000","0163103100000010"); // 20-40%
  } else if (trainConfig == 106){ // PHOS clusters no timing
    cuts.AddCut("54600013","00200009297002008250400000","2444400002033200000","0163103100000010"); // 40-60%
    cuts.AddCut("56800013","00200009297002008250400000","2444400002033200000","0163103100000010"); // 60-80%
    cuts.AddCut("52600013","00200009297002008250400000","2444400002033200000","0163103100000010"); // 20-60%
    cuts.AddCut("54800013","00200009297002008250400000","2444400002033200000","0163103100000010"); // 40-80%
    cuts.AddCut("52500013","00200009297002008250400000","2444400002033200000","0163103100000010"); // 20-50%
  } else if (trainConfig == 107){ // PHOS clusters no timing
    cuts.AddCut("50900013","00200009297002008250400000","2444400002033200000","0163103100000010"); // 0-90%

  //****************************************************************************************************
  // EMCal 5TeV Pb-Pb LHC15o
  //****************************************************************************************************
  } else if (trainConfig == 201){ // EMCAL clusters central
    cuts.AddCut("10110013","00200009327000008250400000","1111100057032230000","0163103100000010"); // 0-10
    cuts.AddCut("30110013","00200009327000008250400000","1111100057032230000","0163103100000010"); // 0-5
    cuts.AddCut("31210013","00200009327000008250400000","1111100057032230000","0163103100000010"); // 5-10
  } else if (trainConfig == 202){ // EMCAL clusters semi-central
    cuts.AddCut("11210013","00200009327000008250400000","1111100057032230000","0163103100000010"); // 10-20
    cuts.AddCut("12310013","00200009327000008250400000","1111100057032230000","0163103100000010"); // 20-30
    cuts.AddCut("13410013","00200009327000008250400000","1111100057032230000","0163103100000010"); // 30-40
    cuts.AddCut("12410013","00200009327000008250400000","1111100057032230000","0163103100000010"); // 20-40
  } else if (trainConfig == 203){ // EMCAL clusters semi-central
    cuts.AddCut("14510013","00200009327000008250400000","1111100057032230000","0163103100000010"); // 40-50
    cuts.AddCut("14610013","00200009327000008250400000","1111100057032230000","0163103100000010"); // 40-60
    cuts.AddCut("15610013","00200009327000008250400000","1111100057032230000","0163103100000010"); // 40-60
  } else if (trainConfig == 204){ // EMCAL clusters peripheral
    cuts.AddCut("16710013","00200009327000008250400000","1111100057032230000","0163103100000010"); // 60-70
    cuts.AddCut("17810013","00200009327000008250400000","1111100057032230000","0163103100000010"); // 70-80
    cuts.AddCut("18910013","00200009327000008250400000","1111100057032230000","0163103100000010"); // 80-90
    cuts.AddCut("16810013","00200009327000008250400000","1111100057032230000","0163103100000010"); // 60-80

  } else if (trainConfig == 205){ // EMCAL clusters central add sig
    cuts.AddCut("10110023","00200009327000008250400000","1111100057032230000","0163103100000010"); // 0-10
    cuts.AddCut("30110023","00200009327000008250400000","1111100057032230000","0163103100000010"); // 0-5
    cuts.AddCut("31210023","00200009327000008250400000","1111100057032230000","0163103100000010"); // 5-10
  } else if (trainConfig == 206){ // EMCAL clusters semi-central add sig
    cuts.AddCut("11210023","00200009327000008250400000","1111100057032230000","0163103100000010"); // 10-20
    cuts.AddCut("12310023","00200009327000008250400000","1111100057032230000","0163103100000010"); // 20-30
    cuts.AddCut("13410023","00200009327000008250400000","1111100057032230000","0163103100000010"); // 30-40
    cuts.AddCut("12410023","00200009327000008250400000","1111100057032230000","0163103100000010"); // 20-40
  } else if (trainConfig == 207){ // EMCAL clusters semi-central add sig
    cuts.AddCut("14510023","00200009327000008250400000","1111100057032230000","0163103100000010"); // 40-50
    cuts.AddCut("14610023","00200009327000008250400000","1111100057032230000","0163103100000010"); // 40-60
    cuts.AddCut("15610023","00200009327000008250400000","1111100057032230000","0163103100000010"); // 40-60
  } else if (trainConfig == 208){ // EMCAL clusters peripheral add sig
    cuts.AddCut("16710023","00200009327000008250400000","1111100057032230000","0163103100000010"); // 60-70
    cuts.AddCut("17810023","00200009327000008250400000","1111100057032230000","0163103100000010"); // 70-80
    cuts.AddCut("18910023","00200009327000008250400000","1111100057032230000","0163103100000010"); // 80-90
    cuts.AddCut("16810023","00200009327000008250400000","1111100057032230000","0163103100000010"); // 60-80

  } else if (trainConfig == 209){ // EMCAL clusters central add sig
    cuts.AddCut("10110023","00200009327000008250400000","1111100057032230000","0163103100000010"); // 0-10
    cuts.AddCut("30110023","00200009327000008250400000","1111100057032230000","0163103100000010"); // 0-5
    cuts.AddCut("31210023","00200009327000008250400000","1111100057032230000","0163103100000010"); // 5-10
  } else if (trainConfig == 210){ // EMCAL clusters semi-central add sig
    cuts.AddCut("11210023","00200009327000008250400000","1111100057032230000","0163103100000010"); // 10-20
    cuts.AddCut("12310023","00200009327000008250400000","1111100057032230000","0163103100000010"); // 20-30
    cuts.AddCut("13410023","00200009327000008250400000","1111100057032230000","0163103100000010"); // 30-40
    cuts.AddCut("12410023","00200009327000008250400000","1111100057032230000","0163103100000010"); // 20-40
  } else if (trainConfig == 211){ // EMCAL clusters semi-central add sig
    cuts.AddCut("14510023","00200009327000008250400000","1111100057032230000","0163103100000010"); // 40-50
    cuts.AddCut("14610023","00200009327000008250400000","1111100057032230000","0163103100000010"); // 40-60
    cuts.AddCut("15610023","00200009327000008250400000","1111100057032230000","0163103100000010"); // 40-60
  } else if (trainConfig == 212){ // EMCAL clusters peripheral add sig
    cuts.AddCut("16710023","00200009327000008250400000","1111100057032230000","0163103100000010"); // 60-70
    cuts.AddCut("17810023","00200009327000008250400000","1111100057032230000","0163103100000010"); // 70-80
    cuts.AddCut("18910023","00200009327000008250400000","1111100057032230000","0163103100000010"); // 80-90
    cuts.AddCut("16810023","00200009327000008250400000","1111100057032230000","0163103100000010"); // 60-80

  } else if (trainConfig == 213){ // EMCAL clusters central TB
    cuts.AddCut("10110013","00200009327000008250400000","1111102057032230000","0163103100000010"); // 0-10
    cuts.AddCut("30110013","00200009327000008250400000","1111102057032230000","0163103100000010"); // 0-5
    cuts.AddCut("31210013","00200009327000008250400000","1111102057032230000","0163103100000010"); // 5-10
  } else if (trainConfig == 214){ // EMCAL clusters semi-central TB
    cuts.AddCut("11210013","00200009327000008250400000","1111102057032230000","0163103100000010"); // 10-20
    cuts.AddCut("12310013","00200009327000008250400000","1111102057032230000","0163103100000010"); // 20-30
    cuts.AddCut("13410013","00200009327000008250400000","1111102057032230000","0163103100000010"); // 30-40
    cuts.AddCut("12410013","00200009327000008250400000","1111102057032230000","0163103100000010"); // 20-40
  } else if (trainConfig == 215){ // EMCAL clusters semi-central TB
    cuts.AddCut("14510013","00200009327000008250400000","1111102057032230000","0163103100000010"); // 40-50
    cuts.AddCut("14610013","00200009327000008250400000","1111102057032230000","0163103100000010"); // 40-60
    cuts.AddCut("15610013","00200009327000008250400000","1111102057032230000","0163103100000010"); // 40-60
  } else if (trainConfig == 216){ // EMCAL clusters peripheral TB
    cuts.AddCut("16710013","00200009327000008250400000","1111102057032230000","0163103100000010"); // 60-70
    cuts.AddCut("17810013","00200009327000008250400000","1111102057032230000","0163103100000010"); // 70-80
    cuts.AddCut("18910013","00200009327000008250400000","1111102057032230000","0163103100000010"); // 80-90
    cuts.AddCut("16810013","00200009327000008250400000","1111102057032230000","0163103100000010"); // 60-80

  } else if (trainConfig == 226){ // EMCAL clusters - peripheral centrality selection for PbPb EMCal
    cuts.AddCut("15910013","00200009327000008250400000","11111020530a2230000","0163103100000010"); //
    cuts.AddCut("15910013","00200009327000008250400000","11111870530a2230000","0163103100000010"); //
    cuts.AddCut("15910013","00200009327000008250400000","11111020530b2230000","0163103100000010"); //
    cuts.AddCut("15910013","00200009327000008250400000","11111870530b2230000","0163103100000010"); //

  } else if (trainConfig == 240){ // EMCAL clusters - 0-90% centrality for PbPb EMCal cluster QA
    cuts.AddCut("50910013","00200009327000008250400000","1111100053032230000","0163103100000010"); //  0-90 calo correction cent dep
    cuts.AddCut("50910013","00200009327000008250400000","1111102053032230000","0163103100000010"); //  0-90 calo correction cent dep
    cuts.AddCut("50910013","00200009327000008250400000","1111183053032230000","0163103100000010"); //  0-90 calo correction cent dep
  } else if (trainConfig == 241){ // EMCAL clusters - 0-90% centrality for PbPb EMCal cluster QA
    cuts.AddCut("50910613","00200009327000008250400000","1111100053032230000","0163103100000010"); //  0-90 calo correction cent dep
    cuts.AddCut("50910613","00200009327000008250400000","1111102053032230000","0163103100000010"); //  0-90 calo correction cent dep
    cuts.AddCut("50910613","00200009327000008250400000","1111183053032230000","0163103100000010"); //  0-90 calo correction cent dep
  } else if (trainConfig == 242){ // EMCAL clusters - centrality selection for PbPb EMCal
    cuts.AddCut("50110013","00200009327000008250400000","1111100053032230000","0163103100000010"); //  0-10 calo correction cent dep
    cuts.AddCut("51210013","00200009327000008250400000","1111100053032230000","0163103100000010"); // 10-20 calo correction cent dep
    cuts.AddCut("52510013","00200009327000008250400000","1111100053032230000","0163103100000010"); // 20-50 calo correction cent dep
    cuts.AddCut("55910013","00200009327000008250400000","1111100053032230000","0163103100000010"); // 50-90 calo correction cent dep
  } else if (trainConfig == 243){ // EMCAL clusters - centrality selection for PbPb EMCal
    cuts.AddCut("50110613","00200009327000008250400000","1111100053032230000","0163103100000010"); //  0-10 calo correction cent dep
    cuts.AddCut("51210613","00200009327000008250400000","1111100053032230000","0163103100000010"); // 10-20 calo correction cent dep
    cuts.AddCut("52510613","00200009327000008250400000","1111100053032230000","0163103100000010"); // 20-50 calo correction cent dep
    cuts.AddCut("55910613","00200009327000008250400000","1111100053032230000","0163103100000010"); // 50-90 calo correction cent dep
  } else if (trainConfig == 244){ // EMCAL clusters - centrality selection for PbPb EMCal
    cuts.AddCut("50110613","00200009327000008250400000","1111102053032230000","0163103100000010"); //  0-10 calo correction cent dep
    cuts.AddCut("51210613","00200009327000008250400000","1111102053032230000","0163103100000010"); // 10-20 calo correction cent dep
    cuts.AddCut("52510613","00200009327000008250400000","1111102053032230000","0163103100000010"); // 20-50 calo correction cent dep
    cuts.AddCut("55910613","00200009327000008250400000","1111102053032230000","0163103100000010"); // 50-90 calo correction cent dep
  } else if (trainConfig == 245){ // EMCAL clusters - centrality selection for PbPb EMCal
    cuts.AddCut("50110613","00200009327000008250400000","1111184053032230000","0163103100000010"); //  0-10 calo correction cent dep
    cuts.AddCut("51210613","00200009327000008250400000","1111185053032230000","0163103100000010"); // 10-20 calo correction cent dep
    cuts.AddCut("52510613","00200009327000008250400000","1111186053032230000","0163103100000010"); // 20-50 calo correction cent dep
    cuts.AddCut("55910613","00200009327000008250400000","1111187053032230000","0163103100000010"); // 50-90 calo correction cent dep
  } else if (trainConfig == 246){ // EMCAL clusters - 0-90% centrality for PbPb EMCal cluster QA
    cuts.AddCut("50910613","00200009327000008250400000","1111183050032230000","0163103100000010"); //  0-90 calo correction cent dep
    cuts.AddCut("50910613","00200009327000008250400000","1111183051032230000","0163103100000010"); //  0-90 calo correction cent dep
    cuts.AddCut("50910613","00200009327000008250400000","1111183057032230000","0163103100000010"); //  0-90 calo correction cent dep
  } else if (trainConfig == 247){ // EMCAL clusters - centrality selection for PbPb EMCal
    cuts.AddCut("50110613","00200009327000008250400000","1111184050032230000","0163103100000010"); //  0-10 calo correction cent dep
    cuts.AddCut("51210613","00200009327000008250400000","1111185050032230000","0163103100000010"); // 10-20 calo correction cent dep
    cuts.AddCut("52510613","00200009327000008250400000","1111186050032230000","0163103100000010"); // 20-50 calo correction cent dep
    cuts.AddCut("55910613","00200009327000008250400000","1111187050032230000","0163103100000010"); // 50-90 calo correction cent dep
  } else if (trainConfig == 248){ // EMCAL clusters - centrality selection for PbPb EMCal
    cuts.AddCut("50110613","00200009327000008250400000","1111184051032230000","0163103100000010"); //  0-10 calo correction cent dep
    cuts.AddCut("51210613","00200009327000008250400000","1111185051032230000","0163103100000010"); // 10-20 calo correction cent dep
    cuts.AddCut("52510613","00200009327000008250400000","1111186051032230000","0163103100000010"); // 20-50 calo correction cent dep
    cuts.AddCut("55910613","00200009327000008250400000","1111187051032230000","0163103100000010"); // 50-90 calo correction cent dep
  } else if (trainConfig == 249){ // EMCAL clusters - centrality selection for PbPb EMCal
    cuts.AddCut("50110613","00200009327000008250400000","1111184057032230000","0163103100000010"); //  0-10 calo correction cent dep
    cuts.AddCut("51210613","00200009327000008250400000","1111185057032230000","0163103100000010"); // 10-20 calo correction cent dep
    cuts.AddCut("52510613","00200009327000008250400000","1111186057032230000","0163103100000010"); // 20-50 calo correction cent dep
    cuts.AddCut("55910613","00200009327000008250400000","1111187057032230000","0163103100000010"); // 50-90 calo correction cent dep
  } else if (trainConfig == 250){ // EMCAL clusters - 0-90% centrality for PbPb EMCal cluster QA
    cuts.AddCut("10910013","00200009327000008250400000","1111183051032230000","0163103100000010"); //  0-90 calo correction cent dep
    cuts.AddCut("10910a13","00200009327000008250400000","1111183051032230000","0163103100000010"); //  0-90 calo correction cent dep
  } else if (trainConfig == 251){ // EMCAL clusters - 0-90% centrality for PbPb EMCal cluster QA
    cuts.AddCut("10110013","00200009327000008250400000","1111183051032230000","0163103100000010"); //  0-90 calo correction cent dep
    cuts.AddCut("11210013","00200009327000008250400000","1111183051032230000","0163103100000010"); //  0-90 calo correction cent dep
    cuts.AddCut("12510013","00200009327000008250400000","1111183051032230000","0163103100000010"); //  0-90 calo correction cent dep
    cuts.AddCut("15910013","00200009327000008250400000","1111183051032230000","0163103100000010"); //  0-90 calo correction cent dep
  } else if (trainConfig == 252){ // EMCAL clusters - 0-90% centrality for PbPb EMCal cluster QA
    cuts.AddCut("10110a13","00200009327000008250400000","1111183051032230000","0163103100000010"); //  0-90 calo correction cent dep
    cuts.AddCut("11210a13","00200009327000008250400000","1111183051032230000","0163103100000010"); //  0-90 calo correction cent dep
    cuts.AddCut("12510a13","00200009327000008250400000","1111183051032230000","0163103100000010"); //  0-90 calo correction cent dep
    cuts.AddCut("15910a13","00200009327000008250400000","1111183051032230000","0163103100000010"); //  0-90 calo correction cent dep
  } else if (trainConfig == 253){ // EMCAL clusters - 20180718 - default without corrections
    cuts.AddCut("10110a13","00200009327000008250400000","1111100051032230000","0163103100000010"); //  0-90 calo correction cent dep
    cuts.AddCut("11210a13","00200009327000008250400000","1111100051032230000","0163103100000010"); //  0-90 calo correction cent dep
    cuts.AddCut("12410a13","00200009327000008250400000","1111100051032230000","0163103100000010"); //  0-90 calo correction cent dep
    cuts.AddCut("14610a13","00200009327000008250400000","1111100051032230000","0163103100000010"); //  0-90 calo correction cent dep
    cuts.AddCut("16810a13","00200009327000008250400000","1111100051032230000","0163103100000010"); //  0-90 calo correction cent dep
  } else if (trainConfig == 254){ // EMCAL clusters - 20180718 - default with peripheral corrections
    cuts.AddCut("10110a13","00200009327000008250400000","1111187051032230000","0163103100000010"); //
    cuts.AddCut("11210a13","00200009327000008250400000","1111187051032230000","0163103100000010"); //
    cuts.AddCut("12410a13","00200009327000008250400000","1111187051032230000","0163103100000010"); //
    cuts.AddCut("14610a13","00200009327000008250400000","1111187051032230000","0163103100000010"); //
    cuts.AddCut("16810a13","00200009327000008250400000","1111187051032230000","0163103100000010"); //
  } else if (trainConfig == 255){ // EMCAL clusters - 20180718 - default with peripheral corrections
    cuts.AddCut("10110a13","00200009327000008250400000","1111184051032230000","0163103100000010"); //
    cuts.AddCut("11210a13","00200009327000008250400000","1111185051032230000","0163103100000010"); //
    cuts.AddCut("12410a13","00200009327000008250400000","1111186051032230000","0163103100000010"); //
    cuts.AddCut("14610a13","00200009327000008250400000","1111187051032230000","0163103100000010"); //
    cuts.AddCut("16810a13","00200009327000008250400000","1111187051032230000","0163103100000010"); //


  } else if (trainConfig == 290){ // EMCAL clusters - correction convcalo f1
    cuts.AddCut("10110013","00200009327000008250400000","1111181003032230000","0163103100000010"); // 0-10
    cuts.AddCut("11210013","00200009327000008250400000","1111181003032230000","0163103100000010"); // 10-20
    cuts.AddCut("12510013","00200009327000008250400000","1111181003032230000","0163103100000010"); // 20-50
    cuts.AddCut("15910013","00200009327000008250400000","1111181003032230000","0163103100000010"); // 50-90
  } else if (trainConfig == 291){ // EMCAL clusters - correction calocalo f2
    cuts.AddCut("10110013","00200009327000008250400000","1111192003032230000","0163103100000010"); // 0-10
    cuts.AddCut("11210013","00200009327000008250400000","1111192003032230000","0163103100000010"); // 10-20
    cuts.AddCut("12510013","00200009327000008250400000","1111192003032230000","0163103100000010"); // 20-50
    cuts.AddCut("15910013","00200009327000008250400000","1111192003032230000","0163103100000010"); // 50-90
  } else if (trainConfig == 292){ // EMCAL clusters - 0-90% centrality for PbPb EMCal cluster QA
    cuts.AddCut("10910013","00200009327000008250400000","1111100003032230000","0163103100000010"); // 0-90
  } else if (trainConfig == 293){ // EMCAL clusters - centrality selection for PbPb EMCal
    cuts.AddCut("10110013","00200009327000008250400000","1111184053032230000","0163103100000010"); //  0-10 calo correction cent dep
    cuts.AddCut("11210013","00200009327000008250400000","1111185053032230000","0163103100000010"); // 10-20 calo correction cent dep
    cuts.AddCut("12510013","00200009327000008250400000","1111186053032230000","0163103100000010"); // 20-50 calo correction cent dep
    cuts.AddCut("15910013","00200009327000008250400000","1111187053032230000","0163103100000010"); // 50-90 calo correction cent dep


  //****************************************************************************************************
  // EMCal 5TeV Xe-Xe LHC17n
  //****************************************************************************************************
  } else if (trainConfig == 300){ // EMCAL clusters - 0-80% centrality
    cuts.AddCut("10810013","00200009327000008250400000","1111100017032230000","0163103100000010"); // 0-80
  } else if (trainConfig == 301){ // EMCAL clusters
    cuts.AddCut("10210013","00200009327000008250400000","1111100017032230000","0163103100000010"); // 0-20
    cuts.AddCut("12410013","00200009327000008250400000","1111100017032230000","0163103100000010"); // 20-40
    cuts.AddCut("10410013","00200009327000008250400000","1111100017032230000","0163103100000010"); // 0-40
    cuts.AddCut("14810013","00200009327000008250400000","1111100017032230000","0163103100000010"); // 40-80
  } else if (trainConfig == 302){ // EMCAL clusters - 0-80% centrality for EMCal cluster QA TB calib
    cuts.AddCut("10810013","00200009327000008250400000","1111102017032230000","0163103100000010"); // 0-80
  } else if (trainConfig == 303){ // EMCAL clusters TB calib
    cuts.AddCut("10210013","00200009327000008250400000","1111102017032230000","0163103100000010"); // 0-20
    cuts.AddCut("12410013","00200009327000008250400000","1111102017032230000","0163103100000010"); // 20-40
    cuts.AddCut("10410013","00200009327000008250400000","1111102017032230000","0163103100000010"); // 0-40
    cuts.AddCut("14810013","00200009327000008250400000","1111102017032230000","0163103100000010"); // 40-80
  } else if (trainConfig == 304){ // EMCAL clusters -  EMCal cluster QA TB calib, PCM min pt 0.02
    cuts.AddCut("10810013","00200089327000008250400000","1111102017032230000","0163103100000010"); // 0-80
    cuts.AddCut("10210013","00200089327000008250400000","1111102017032230000","0163103100000010"); // 0-20
    cuts.AddCut("12410013","00200089327000008250400000","1111102017032230000","0163103100000010"); // 20-40
    cuts.AddCut("10410013","00200089327000008250400000","1111102017032230000","0163103100000010"); // 0-40
    cuts.AddCut("14810013","00200089327000008250400000","1111102017032230000","0163103100000010"); // 40-80

  } else if (trainConfig == 310){ // EMCAL clusters  non lin vars 0-20
    cuts.AddCut("10210013","00200009327000008250400000","1111171017032230000","0163103100000010");
    cuts.AddCut("10210013","00200009327000008250400000","1111172017032230000","0163103100000010");
    cuts.AddCut("10210013","00200009327000008250400000","1111181017032230000","0163103100000010");
    cuts.AddCut("10210013","00200009327000008250400000","1111182017032230000","0163103100000010");
  } else if (trainConfig == 320){ // EMCAL clusters  non lin vars 20-40
    cuts.AddCut("12410013","00200009327000008250400000","1111171007032230000","0163103100000010");
    cuts.AddCut("12410013","00200009327000008250400000","1111172007032230000","0163103100000010");
    cuts.AddCut("12410013","00200009327000008250400000","1111181007032230000","0163103100000010");
    cuts.AddCut("12410013","00200009327000008250400000","1111182007032230000","0163103100000010");
  } else if (trainConfig == 330){ // EMCAL clusters  non lin vars 40-80
    cuts.AddCut("14810013","00200009327000008250400000","1111171007032230000","0163103100000010");
    cuts.AddCut("14810013","00200009327000008250400000","1111172007032230000","0163103100000010");
    cuts.AddCut("14810013","00200009327000008250400000","1111181007032230000","0163103100000010");
    cuts.AddCut("14810013","00200009327000008250400000","1111182007032230000","0163103100000010");
  } else if (trainConfig == 340){ // EMCAL clusters  non lin vars 0-40
    cuts.AddCut("10410013","00200009327000008250400000","1111171007032230000","0163103100000010");
    cuts.AddCut("10410013","00200009327000008250400000","1111172007032230000","0163103100000010");
    cuts.AddCut("10410013","00200009327000008250400000","1111181007032230000","0163103100000010");
    cuts.AddCut("10410013","00200009327000008250400000","1111182007032230000","0163103100000010");
  } else if (trainConfig == 350){ // EMCAL clusters  non lin vars 0-80
    cuts.AddCut("10810013","00200009327000008250400000","1111171007032230000","0163103100000010");
    cuts.AddCut("10810013","00200009327000008250400000","1111172007032230000","0163103100000010");
    cuts.AddCut("10810013","00200009327000008250400000","1111181007032230000","0163103100000010");
    cuts.AddCut("10810013","00200009327000008250400000","1111182007032230000","0163103100000010");

  //****************************************************************************************************
  // PHOS 5TeV Xe-Xe LHC17n
  //****************************************************************************************************
  } else if (trainConfig == 400){ // PHOS clusters - 0-90% centrality
    cuts.AddCut("10810013","00200009327000008250400000","2446600011013200000","0163103100000010"); // 0-80%
    cuts.AddCut("10810013","00200009327000008250400000","2446601011013200000","0163103100000010"); // 0-80%
  } else if (trainConfig == 401) {  // PHOS Cent dep Xe-Xe
    cuts.AddCut("10210013","00200009327000008250400000","2446600011013200000","0163103100000010"); // 0-20%
    cuts.AddCut("12410013","00200009327000008250400000","2446600011013200000","0163103100000010"); // 20-40%
    cuts.AddCut("10410013","00200009327000008250400000","2446600011013200000","0163103100000010"); // 0-40%
    cuts.AddCut("14810013","00200009327000008250400000","2446600011013200000","0163103100000010"); // 40-80%
  } else if (trainConfig == 402) {  // PHOS Cent dep diff timing
    cuts.AddCut("10210013","00200009327000008250400000","2446601011013200000","0163103100000010"); // 0-20%
    cuts.AddCut("12410013","00200009327000008250400000","2446601011013200000","0163103100000010"); // 20-40%
    cuts.AddCut("10410013","00200009327000008250400000","2446601011013200000","0163103100000010"); // 0-40%
    cuts.AddCut("14810013","00200009327000008250400000","2446601011013200000","0163103100000010"); // 40-80%
  } else if (trainConfig == 403) {  // PHOS Cent dep diff timing, min pT PCM 0.02
    cuts.AddCut("10210013","00200089327000008250400000","2446601011013200000","0163103100000010"); // 0-20%
    cuts.AddCut("12410013","00200089327000008250400000","2446601011013200000","0163103100000010"); // 20-40%
    cuts.AddCut("10410013","00200089327000008250400000","2446601011013200000","0163103100000010"); // 0-40%
    cuts.AddCut("14810013","00200089327000008250400000","2446601011013200000","0163103100000010"); // 40-80%
  } else if (trainConfig == 410) {  // PHOS Cluster non lin vars 0-20
    cuts.AddCut("10210013","00200009327000008250400000","2446601011013200000","0163103100000010");
    cuts.AddCut("10210013","00200009327000008250400000","2446671011013200000","0163103100000010");
    cuts.AddCut("10210013","00200009327000008250400000","2446672011013200000","0163103100000010");
    cuts.AddCut("10210013","00200009327000008250400000","2446681011013200000","0163103100000010");
    cuts.AddCut("10210013","00200009327000008250400000","2446682011013200000","0163103100000010");
  } else if (trainConfig == 420) {  // PHOS Cluster non lin vars 20-40
    cuts.AddCut("12410013","00200009327000008250400000","2446601011013200000","0163103100000010");
    cuts.AddCut("12410013","00200009327000008250400000","2446671011013200000","0163103100000010");
    cuts.AddCut("12410013","00200009327000008250400000","2446672011013200000","0163103100000010");
    cuts.AddCut("12410013","00200009327000008250400000","2446681011013200000","0163103100000010");
    cuts.AddCut("12410013","00200009327000008250400000","2446682011013200000","0163103100000010");
  } else if (trainConfig == 430) {  // PHOS Cluster non lin vars 40-80
    cuts.AddCut("14810013","00200009327000008250400000","2446601011013200000","0163103100000010");
    cuts.AddCut("14810013","00200009327000008250400000","2446671011013200000","0163103100000010");
    cuts.AddCut("14810013","00200009327000008250400000","2446672011013200000","0163103100000010");
    cuts.AddCut("14810013","00200009327000008250400000","2446681011013200000","0163103100000010");
    cuts.AddCut("14810013","00200009327000008250400000","2446682011013200000","0163103100000010");
  } else if (trainConfig == 440) {  // PHOS Cluster non lin vars 0-40
    cuts.AddCut("10410013","00200009327000008250400000","2446601011013200000","0163103100000010");
    cuts.AddCut("10410013","00200009327000008250400000","2446671011013200000","0163103100000010");
    cuts.AddCut("10410013","00200009327000008250400000","2446672011013200000","0163103100000010");
    cuts.AddCut("10410013","00200009327000008250400000","2446681011013200000","0163103100000010");
    cuts.AddCut("10410013","00200009327000008250400000","2446682011013200000","0163103100000010");
  } else if (trainConfig == 450) {  // PHOS Cluster non lin vars 0-80
    cuts.AddCut("10810013","00200009327000008250400000","2446601011013200000","0163103100000010");
    cuts.AddCut("10810013","00200009327000008250400000","2446671011013200000","0163103100000010");
    cuts.AddCut("10810013","00200009327000008250400000","2446672011013200000","0163103100000010");
    cuts.AddCut("10810013","00200009327000008250400000","2446681011013200000","0163103100000010");
    cuts.AddCut("10810013","00200009327000008250400000","2446682011013200000","0163103100000010");

  //****************************************************************************************************
  // PHOS 5TeV Pb-Pb LHC15o
  //****************************************************************************************************
  } else if (trainConfig == 601){ // EMCAL clusters central
    cuts.AddCut("10110013","00200009327000008250400000","2446600051013200000","0163103100000010"); // 0-10
    cuts.AddCut("30110013","00200009327000008250400000","2446600051013200000","0163103100000010"); // 0-5
    cuts.AddCut("31210013","00200009327000008250400000","2446600051013200000","0163103100000010"); // 5-10
  } else if (trainConfig == 602){ // EMCAL clusters semi-central
    cuts.AddCut("11210013","00200009327000008250400000","2446600051013200000","0163103100000010"); // 10-20
    cuts.AddCut("12310013","00200009327000008250400000","2446600051013200000","0163103100000010"); // 20-30
    cuts.AddCut("13410013","00200009327000008250400000","2446600051013200000","0163103100000010"); // 30-40
    cuts.AddCut("12410013","00200009327000008250400000","2446600051013200000","0163103100000010"); // 20-40
  } else if (trainConfig == 603){ // EMCAL clusters semi-central
    cuts.AddCut("14510013","00200009327000008250400000","2446600051013200000","0163103100000010"); // 40-50
    cuts.AddCut("14610013","00200009327000008250400000","2446600051013200000","0163103100000010"); // 40-60
    cuts.AddCut("15610013","00200009327000008250400000","2446600051013200000","0163103100000010"); // 40-60
  } else if (trainConfig == 604){ // EMCAL clusters peripheral
    cuts.AddCut("16710013","00200009327000008250400000","2446600051013200000","0163103100000010"); // 60-70
    cuts.AddCut("17810013","00200009327000008250400000","2446600051013200000","0163103100000010"); // 70-80
    cuts.AddCut("18910013","00200009327000008250400000","2446600051013200000","0163103100000010"); // 80-90
    cuts.AddCut("16810013","00200009327000008250400000","2446600051013200000","0163103100000010"); // 60-80
  } else if (trainConfig == 605){ // EMCAL clusters central
    cuts.AddCut("10110013","00200009327000008250400000","2446601051013200000","0163103100000010"); // 0-10
    cuts.AddCut("30110013","00200009327000008250400000","2446601051013200000","0163103100000010"); // 0-5
    cuts.AddCut("31210013","00200009327000008250400000","2446601051013200000","0163103100000010"); // 5-10
  } else if (trainConfig == 606){ // EMCAL clusters semi-central
    cuts.AddCut("11210013","00200009327000008250400000","2446601051013200000","0163103100000010"); // 10-20
    cuts.AddCut("12310013","00200009327000008250400000","2446601051013200000","0163103100000010"); // 20-30
    cuts.AddCut("13410013","00200009327000008250400000","2446601051013200000","0163103100000010"); // 30-40
    cuts.AddCut("12410013","00200009327000008250400000","2446601051013200000","0163103100000010"); // 20-40
  } else if (trainConfig == 607){ // EMCAL clusters semi-central
    cuts.AddCut("14510013","00200009327000008250400000","2446601051013200000","0163103100000010"); // 40-50
    cuts.AddCut("14610013","00200009327000008250400000","2446601051013200000","0163103100000010"); // 40-60
    cuts.AddCut("15610013","00200009327000008250400000","2446601051013200000","0163103100000010"); // 40-60
  } else if (trainConfig == 608){ // EMCAL clusters peripheral
    cuts.AddCut("16710013","00200009327000008250400000","2446601051013200000","0163103100000010"); // 60-70
    cuts.AddCut("17810013","00200009327000008250400000","2446601051013200000","0163103100000010"); // 70-80
    cuts.AddCut("18910013","00200009327000008250400000","2446601051013200000","0163103100000010"); // 80-90
    cuts.AddCut("16810013","00200009327000008250400000","2446601051013200000","0163103100000010"); // 60-80


  //****************************************************************************************************
  // EMCal 5TeV Pb-Pb LHC15o EMC triggers
  //****************************************************************************************************
  } else if (trainConfig == 201){ // EMCAL clusters central
    cuts.AddCut("10183013","00200009327000008250400000","1111100057032230000","0163103100000010"); // 0-10
    cuts.AddCut("30183013","00200009327000008250400000","1111100057032230000","0163103100000010"); // 0-5
    cuts.AddCut("31283013","00200009327000008250400000","1111100057032230000","0163103100000010"); // 5-10
  } else if (trainConfig == 202){ // EMCAL clusters semi-central
    cuts.AddCut("11283013","00200009327000008250400000","1111100057032230000","0163103100000010"); // 10-20
    cuts.AddCut("12383013","00200009327000008250400000","1111100057032230000","0163103100000010"); // 20-30
    cuts.AddCut("13483013","00200009327000008250400000","1111100057032230000","0163103100000010"); // 30-40
    cuts.AddCut("12483013","00200009327000008250400000","1111100057032230000","0163103100000010"); // 20-40
  } else if (trainConfig == 203){ // EMCAL clusters semi-central
    cuts.AddCut("14583013","00200009327000008250400000","1111100057032230000","0163103100000010"); // 40-50
    cuts.AddCut("14683013","00200009327000008250400000","1111100057032230000","0163103100000010"); // 40-60
    cuts.AddCut("15683013","00200009327000008250400000","1111100057032230000","0163103100000010"); // 40-60
  } else if (trainConfig == 204){ // EMCAL clusters peripheral
    cuts.AddCut("16783013","00200009327000008250400000","1111100057032230000","0163103100000010"); // 60-70
    cuts.AddCut("17883013","00200009327000008250400000","1111100057032230000","0163103100000010"); // 70-80
    cuts.AddCut("18983013","00200009327000008250400000","1111100057032230000","0163103100000010"); // 80-90
    cuts.AddCut("16883013","00200009327000008250400000","1111100057032230000","0163103100000010"); // 60-80


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

  Int_t numberOfCuts = cuts.GetNCuts();

  TList *EventCutList     = new TList();
  TList *ConvCutList      = new TList();
  TList *ClusterCutList   = new TList();
  TList *MesonCutList     = new TList();

  TList *HeaderList       = new TList();
  if (periodName.CompareTo("LHC13d2")==0){
    TObjString *Header1 = new TObjString("pi0_1");
    HeaderList->Add(Header1);
  //    TObjString *Header3 = new TObjString("eta_2");
  //    HeaderList->Add(Header3);

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
  } else if (periodName.CompareTo("LHC16h4b")==0 || periodName.CompareTo("LHC16h4b2")==0){
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
    TString TrackMatcherName = Form("CaloTrackMatcher_%s_%i",caloCutPos.Data(),trackMatcherRunningMode);
    if( !(AliCaloTrackMatcher*)mgr->GetTask(TrackMatcherName.Data()) ){
      AliCaloTrackMatcher* fTrackMatcher = new AliCaloTrackMatcher(TrackMatcherName.Data(),caloCutPos.Atoi(),trackMatcherRunningMode);
      fTrackMatcher->SetV0ReaderName(V0ReaderName);
      fTrackMatcher->SetCorrectionTaskSetting(corrTaskSetting);
      mgr->AddTask(fTrackMatcher);
      mgr->ConnectInput(fTrackMatcher,0,cinput);
    }

    analysisEventCuts[i] = new AliConvEventCuts();
//     if ( trainConfig == 1){
//       if (periodName.CompareTo("LHC14a1a") ==0 || periodName.CompareTo("LHC14a1b") ==0 || periodName.CompareTo("LHC14a1c") ==0 ){
//         if ( i == 0 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_0005TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_0005TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0005V0M","Eta_Fit_Data_PbPb_2760GeV_0005V0M");
//         if ( i == 1 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_0510TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_0510TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0510V0M","Eta_Fit_Data_PbPb_2760GeV_0510V0M");
//         if ( i == 2 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_0010TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_0010TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0010V0M","Eta_Fit_Data_PbPb_2760GeV_0010V0M");
//         if ( i == 3 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_2040TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_2040TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_2040V0M","Eta_Fit_Data_PbPb_2760GeV_2040V0M");
//         if ( i == 4 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_2050TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_2050TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_2050V0M","Eta_Fit_Data_PbPb_2760GeV_2050V0M");
//       }
//     }
    analysisEventCuts[i]->SetV0ReaderName(V0ReaderName);
    if (periodNameV0Reader.CompareTo("") != 0) analysisEventCuts[i]->SetPeriodEnum(periodNameV0Reader);

    if(periodName.Contains("LHC11h") && doFlattening){
      cout << "entering the cent. flattening loop -> searching for file: " << fileNameInputForCentFlattening.Data() << endl;

      if( fileNameInputForCentFlattening.Contains("FlatFile") ){
        analysisEventCuts[i]->SetUseWeightFlatCentralityFromFile(doFlattening, fileNameInputForCentFlattening, "Cent");
      } else if( fileNameInputForCentFlattening.Contains("Good") ){
        analysisEventCuts[i]->SetUseWeightFlatCentralityFromFile(doFlattening, fileNameInputForCentFlattening, "CentGoodRuns");
      }else if( fileNameInputForCentFlattening.Contains("SemiGood") ){
        analysisEventCuts[i]->SetUseWeightFlatCentralityFromFile(doFlattening, fileNameInputForCentFlattening, "CentSemiGoodRuns");
      }else {
        analysisEventCuts[i]->SetUseWeightFlatCentralityFromFile(doFlattening, fileNameInputForCentFlattening, "CentTotalRuns");
      }
    }

    TString dataInputMultHisto  = "";
    TString mcInputMultHisto    = "";
    if (doMultiplicityWeighting){
      cout << "INFO enableling mult weighting" << endl;
      if(periodNameAnchor.CompareTo("LHC15o")==0){
        TString cutNumber = cuts.GetEventCut(i);
        TString centCut = cutNumber(0,3);  // first three digits of event cut
        dataInputMultHisto = Form("%s_%s", periodNameAnchor.Data(), centCut.Data());
        mcInputMultHisto   = Form("%s_%s", periodName.Data(), centCut.Data());
        cout << "INFO read " << dataInputMultHisto.Data() << " and " <<  mcInputMultHisto.Data() << " from " << fileNameInputForMultWeighing.Data() << endl;
      }
      analysisEventCuts[i]->SetUseWeightMultiplicityFromFile(kTRUE, fileNameInputForMultWeighing, dataInputMultHisto, mcInputMultHisto );
    }

    if (trainConfig == 34 || trainConfig == 35){
      if (periodName.CompareTo("LHC14a1a") ==0 || periodName.CompareTo("LHC14a1b") ==0 || periodName.CompareTo("LHC14a1c") ==0 ){
        if ( doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_0010TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_0010TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0010V0M","Eta_Fit_Data_PbPb_2760GeV_0010V0M");
      }
    }
    if (trainConfig == 36 || trainConfig == 37){
      if (periodName.CompareTo("LHC14a1a") ==0 || periodName.CompareTo("LHC14a1b") ==0 || periodName.CompareTo("LHC14a1c") ==0 ){
        if ( doWeighting)   analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_2050TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_2050TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_2050V0M","Eta_Fit_Data_PbPb_2760GeV_2050V0M");
      }
    }

    if (trainConfig == 38 || trainConfig == 39 || trainConfig == 43 || trainConfig == 44){
      if (periodName.CompareTo("LHC14a1a") ==0 || periodName.CompareTo("LHC14a1b") ==0 || periodName.CompareTo("LHC14a1c") ==0 ){
        if ( doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_0010TPC",periodName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_0010TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0010V0M","Eta_Fit_Data_PbPb_2760GeV_0010V0M");
      }
    }

    if (trainConfig == 40 || trainConfig == 41){
      if (periodName.CompareTo("LHC14a1a") ==0 || periodName.CompareTo("LHC14a1b") ==0 || periodName.CompareTo("LHC14a1c") ==0 ){
        if ( doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_2050TPC",periodName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_2050TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_2050V0M","Eta_Fit_Data_PbPb_2760GeV_2050V0M");
      }
    }

    if (runLightOutput > 0) analysisEventCuts[i]->SetLightOutput(kTRUE);
    analysisEventCuts[i]->InitializeCutsFromCutString((cuts.GetEventCut(i)).Data());
    if (periodName.CompareTo("LHC14a1b") ==0 || periodName.CompareTo("LHC14a1c") ==0 ){
      if (headerSelectionInt == 1 || headerSelectionInt == 4 || headerSelectionInt == 12 ) analysisEventCuts[i]->SetAddedSignalPDGCode(111);
      if (headerSelectionInt == 2 || headerSelectionInt == 5 || headerSelectionInt == 13 ) analysisEventCuts[i]->SetAddedSignalPDGCode(221);
    }
    analysisEventCuts[i]->SetCorrectionTaskSetting(corrTaskSetting);
    EventCutList->Add(analysisEventCuts[i]);
    analysisEventCuts[i]->SetFillCutHistograms("",kFALSE);

    analysisCuts[i] = new AliConversionPhotonCuts();
    analysisCuts[i]->SetV0ReaderName(V0ReaderName);
    if (runLightOutput > 0) analysisCuts[i]->SetLightOutput(kTRUE);
    analysisCuts[i]->InitializeCutsFromCutString((cuts.GetPhotonCut(i)).Data());
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

  task->SetAllowOverlapHeaders(enableHeaderOverlap);
  task->SetEventCutList(numberOfCuts,EventCutList);
  task->SetConversionCutList(numberOfCuts,ConvCutList);
  task->SetCaloCutList(numberOfCuts,ClusterCutList);
  task->SetMesonCutList(numberOfCuts,MesonCutList);
  task->SetCorrectionTaskSetting(corrTaskSetting);
  task->SetMoveParticleAccordingToVertex(kTRUE);
  task->SetDoMesonAnalysis(kTRUE);
  task->SetDoMesonQA(enableQAMesonTask); //Attention new switch for Pi0 QA
  task->SetDoPhotonQA(enableQAPhotonTask);  //Attention new switch small for Photon QA
  task->SetDoClusterQA(1);  //Attention new switch small for Cluster QA
  task->SetEnableSortingOfMCClusLabels(enableSortingMCLabels);
  task->SetDoTreeInvMassShowerShape(doTreeClusterShowerShape);
  task->SetUseTHnSparse(isUsingTHnSparse);
  if(enableExtMatchAndQA > 1){ task->SetPlotHistsExtQA(kTRUE);}

  //connect containers
  AliAnalysisDataContainer *coutput =
    mgr->CreateContainer(!(corrTaskSetting.CompareTo("")) ? Form("GammaConvCalo_%i",trainConfig) : Form("GammaConvCalo_%i_%s",trainConfig,corrTaskSetting.Data()), TList::Class(),
              AliAnalysisManager::kOutputContainer,Form("GammaConvCalo_%i.root",trainConfig));

  mgr->AddTask(task);
  mgr->ConnectInput(task,0,cinput);
  mgr->ConnectOutput(task,1,coutput);
  if(enableQAPhotonTask>1){
    for(Int_t i = 0; i<numberOfCuts; i++){
      mgr->ConnectOutput(task,2+i,mgr->CreateContainer(Form("%s_%s_%s_%s Photon DCA tree",(cuts.GetEventCut(i)).Data(),(cuts.GetPhotonCut(i)).Data(),cuts.GetClusterCut(i).Data(),(cuts.GetMesonCut(i)).Data()), TTree::Class(), AliAnalysisManager::kOutputContainer, Form("GammaConvCalo_%i.root",trainConfig)) );
    }
  }

  return;

}
