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
//($ALIPHYSICS/PWGGA/GammaConv/AliAnalysisTaskGammaCalo.cxx) for
//PbPb together with all supporting classes
//***************************************************************************************

//***************************************************************************************
//CutHandler contains all cuts for a certain analysis and trainconfig,
//it automatically checks length of cutStrings and takes care of the number of added cuts,
//no specification of the variable 'numberOfCuts' needed anymore.
//***************************************************************************************
class CutHandlerCalo{
  public:
    CutHandlerCalo(Int_t nMax=10){
      nCuts=0; nMaxCuts=nMax; validCuts = true;
      eventCutArray = new TString[nMaxCuts]; clusterCutArray = new TString[nMaxCuts]; mesonCutArray = new TString[nMaxCuts];
      for(Int_t i=0; i<nMaxCuts; i++) {eventCutArray[i] = ""; clusterCutArray[i] = ""; mesonCutArray[i] = "";}
    }

    void AddCut(TString eventCut, TString clusterCut, TString mesonCut){
      if(nCuts>=nMaxCuts) {cout << "ERROR in CutHandlerCalo: Exceeded maximum number of cuts!" << endl; validCuts = false; return;}
      if( eventCut.Length()!=8 || clusterCut.Length()!=19 || mesonCut.Length()!=16 ) {cout << "ERROR in CutHandlerCalo: Incorrect length of cut string!" << endl; validCuts = false; return;}
      eventCutArray[nCuts]=eventCut; clusterCutArray[nCuts]=clusterCut; mesonCutArray[nCuts]=mesonCut;
      nCuts++;
      return;
    }
    Bool_t AreValid(){return validCuts;}
    Int_t GetNCuts(){if(validCuts) return nCuts; else return 0;}
    TString GetEventCut(Int_t i){if(validCuts&&i<nMaxCuts&&i>=0) return eventCutArray[i]; else{cout << "ERROR in CutHandlerCalo: GetEventCut wrong index i" << endl;return "";}}
    TString GetClusterCut(Int_t i){if(validCuts&&i<nMaxCuts&&i>=0) return clusterCutArray[i]; else {cout << "ERROR in CutHandlerCalo: GetClusterCut wrong index i" << endl;return "";}}
    TString GetMesonCut(Int_t i){if(validCuts&&i<nMaxCuts&&i>=0) return mesonCutArray[i]; else {cout << "ERROR in CutHandlerCalo: GetMesonCut wrong index i" << endl;return "";}}
  private:
    Bool_t validCuts;
    Int_t nCuts; Int_t nMaxCuts;
    TString* eventCutArray;
    TString* clusterCutArray;
    TString* mesonCutArray;
};

//***************************************************************************************
//main function
//***************************************************************************************
void AddTask_GammaCalo_PbPb(  Int_t     trainConfig                     = 1,                    // change different set of cuts
                              Int_t     isMC                            = 0,                    // run MC
                              Int_t     enableQAMesonTask               = 0,                    // enable QA in AliAnalysisTaskGammaConvV1
                              Int_t     enableQAClusterTask             = 0,                    // enable additional QA task
                              TString   fileNameInputForWeighting       = "MCSpectraInput.root",// path to file for weigting input / modified acceptance
                              Int_t     headerSelectionInt              = 0,                    // 1 pi0 header, 2 eta header, 3 both (only for "named" boxes)
                              Bool_t    enableHeaderOverlap             = kTRUE,                // allow overlapping header for the clusters
                              TString   cutnumberAODBranch              = "111110006008400000001500000",
                              TString   periodName                      = "LHC13d2",            // name of the period for added signals and weighting
                              Bool_t    doWeighting                     = kFALSE,               // enable Weighting
                              Bool_t    isUsingTHnSparse                = kTRUE,                // enable or disable usage of THnSparses for background estimation
                              Int_t     enableExtMatchAndQA             = 0,                    // disabled (0), extMatch (1), extQA_noCellQA (2), extMatch+extQA_noCellQA (3), extQA+cellQA (4), extMatch+extQA+cellQA (5)
                              TString   periodNameV0Reader              = "",                   // name of period for V0Reader
                              Bool_t    enableSortingMCLabels           = kTRUE,                // enable sorting for MC cluster labels
                              Bool_t    runLightOutput                  = kFALSE,               // switch to run light output (only essential histograms for afterburner)
                              Bool_t    doFlattening                    = kFALSE,               // switch on centrality flattening for LHC11h
                              TString   fileNameInputForCentFlattening  = "",                   // file name for centrality flattening
                              Bool_t    doMultiplicityWeighting         = kFALSE,                         //
                              TString   fileNameInputForMultWeighing    = "Multiplicity.root",            //
                              TString   periodNameAnchor                = "",
                              TString   additionalTrainConfig           = "0"                   // additional counter for trainconfig
                ) {

  Int_t trackMatcherRunningMode = 0; // CaloTrackMatcher running mode
  Bool_t doTreeEOverP = kFALSE; // switch to produce EOverP tree
  TH1S* histoAcc = 0x0;         // histo for modified acceptance
  TString corrTaskSetting = ""; // select which correction task setting to use
  //parse additionalTrainConfig flag
  TObjArray *rAddConfigArr = additionalTrainConfig.Tokenize("_");
  if(rAddConfigArr->GetEntries()<1){cout << "ERROR: AddTask_GammaCalo_PbPb during parsing of additionalTrainConfig String '" << additionalTrainConfig.Data() << "'" << endl; return;}
  TObjString* rAdditionalTrainConfig;
  for(Int_t i = 0; i<rAddConfigArr->GetEntries() ; i++){
    if(i==0) rAdditionalTrainConfig = (TObjString*)rAddConfigArr->At(i);
    else{
      TObjString* temp = (TObjString*) rAddConfigArr->At(i);
      TString tempStr = temp->GetString();
      if(tempStr.CompareTo("EPCLUSTree") == 0){
        cout << "INFO: AddTask_GammaCalo_PbPb activating 'EPCLUSTree'" << endl;
        doTreeEOverP = kTRUE;
      }else if(tempStr.BeginsWith("MODIFYACC")){
        cout << "INFO: AddTask_GammaCalo_PbPb activating 'MODIFYACC'" << endl;
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
        cout << "INFO: AddTask_GammaCalo_PbPb will use custom branch from Correction Framework!" << endl;
        corrTaskSetting = tempStr;
        corrTaskSetting.Replace(0,2,"");
      }else if(tempStr.BeginsWith("TM")){
        TString tempType = tempStr;
        tempType.Replace(0,2,"");
        trackMatcherRunningMode = tempType.Atoi();
        cout << Form("INFO: AddTask_GammaCalo_PbPb will use running mode '%i' for the TrackMatcher!",trackMatcherRunningMode) << endl;
      }
    }
  }
  TString sAdditionalTrainConfig = rAdditionalTrainConfig->GetString();
  if (sAdditionalTrainConfig.Atoi() > 0){
    trainConfig = trainConfig + sAdditionalTrainConfig.Atoi();
    cout << "INFO: AddTask_GammaCalo_PbPb running additionalTrainConfig '" << sAdditionalTrainConfig.Atoi() << "', train config: '" << trainConfig << "'" << endl;
  }
  cout << "corrTaskSetting: " << corrTaskSetting.Data() << endl;

  Int_t isHeavyIon = 1;

  // ================== GetAnalysisManager ===============================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error(Form("AddTask_GammaCalo_PbPb_%i",trainConfig), "No analysis manager found.");
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
  AliAnalysisTaskGammaCalo *task=NULL;
  task= new AliAnalysisTaskGammaCalo(Form("GammaCalo_%i",trainConfig));
  task->SetIsHeavyIon(isHeavyIon);
  task->SetIsMC(isMC);
  task->SetV0ReaderName(V0ReaderName);
  task->SetCorrectionTaskSetting(corrTaskSetting);
  task->SetLightOutput(runLightOutput);
  task->SetTrackMatcherRunningMode(trackMatcherRunningMode);

  //create cut handler
  CutHandlerCalo cuts;

  // meson cuts
  // meson type (Dalitz or not), BG scheme, pool depth, rotation degrees, rapidity cut, radius cut, alpha, chi2, shared electrons, reject to close v0, MC smearing, dca, dca, dca

  if (trainConfig == 1){ // EMCAL clusters
    cuts.AddCut("60100013","1111100053032230000","0163103100000050"); // 0-5%
    cuts.AddCut("61200013","1111100053032230000","0163103100000050"); // 5-10%
    cuts.AddCut("50100013","1111100053032230000","0163103100000050"); // 0-10%
    cuts.AddCut("52400013","1111100053032230000","0163103100000050"); // 20-40%
    cuts.AddCut("52500013","1111100053032230000","0163103100000050"); // 20-50%
  } else if (trainConfig == 2){ // EMCAL clusters no timing cut
    cuts.AddCut("60100013","1111100003032230000","0163103100000050"); // 0-5%
    cuts.AddCut("61200013","1111100003032230000","0163103100000050"); // 5-10%
    cuts.AddCut("50100013","1111100003032230000","0163103100000050"); // 0-10%
    cuts.AddCut("52400013","1111100003032230000","0163103100000050"); // 20-40%
    cuts.AddCut("52500013","1111100003032230000","0163103100000050"); // 20-50%
  } else if (trainConfig == 3){ // EMCAL clusters
    cuts.AddCut("60100013","1111100053032230000","0163103100000050"); // 0-5%
    cuts.AddCut("61200013","1111100053032230000","0163103100000050"); // 5-10%
    cuts.AddCut("50100013","1111100053032230000","0163103100000050"); // 0-10%
    cuts.AddCut("51200013","1111100053032230000","0163103100000050"); // 10-20%
    cuts.AddCut("52400013","1111100053032230000","0163103100000050"); // 20-40%
  } else if (trainConfig == 4){ // EMCAL clusters no timing cut
    cuts.AddCut("60100013","1111100003032230000","0163103100000050"); // 0-5%
    cuts.AddCut("61200013","1111100003032230000","0163103100000050"); // 5-10%
    cuts.AddCut("50100013","1111100003032230000","0163103100000050"); // 0-10%
    cuts.AddCut("51200013","1111100003032230000","0163103100000050"); // 10-20%
    cuts.AddCut("52400013","1111100003032230000","0163103100000050"); // 20-40%
  } else if (trainConfig == 5){ // EMCAL clusters
    cuts.AddCut("54600013","1111100053032230000","0163103100000050"); // 40-60%
    cuts.AddCut("56800013","1111100053032230000","0163103100000050"); // 60-80%
    cuts.AddCut("52600013","1111100053032230000","0163103100000050"); // 20-60%
    cuts.AddCut("54800013","1111100053032230000","0163103100000050"); // 40-80%
    cuts.AddCut("52500013","1111100053032230000","0163103100000050"); // 20-50%
  } else if (trainConfig == 6){ // EMCAL clusters  no timing cut
    cuts.AddCut("54600013","1111100003032230000","0163103100000050"); // 40-60%
    cuts.AddCut("56800013","1111100003032230000","0163103100000050"); // 60-80%
    cuts.AddCut("52600013","1111100003032230000","0163103100000050"); // 20-60%
    cuts.AddCut("54800013","1111100003032230000","0163103100000050"); // 40-80%
    cuts.AddCut("52500013","1111100003032230000","0163103100000050"); // 20-50%
  } else if (trainConfig == 7){ // EMCAL clusters
    cuts.AddCut("60100013","1111100053032230000","0163103100000050"); // 0-5%
    cuts.AddCut("61200013","1111100053032230000","0163103100000050"); // 5-10%
    cuts.AddCut("50100013","1111100053032230000","0163103100000050"); // 0-10%
  } else if (trainConfig == 8){ // EMCAL clusters no timing cut
    cuts.AddCut("60100013","1111100003032230000","0163103100000050"); // 0-5%
    cuts.AddCut("61200013","1111100003032230000","0163103100000050"); // 5-10%
    cuts.AddCut("50100013","1111100003032230000","0163103100000050"); // 0-10%
  } else if (trainConfig == 9){ // EMCAL clusters
    cuts.AddCut("51200013","1111100053032230000","0163103100000050"); // 10-20%
    cuts.AddCut("52400013","1111100053032230000","0163103100000050"); // 20-40%
    cuts.AddCut("52500013","1111100053032230000","0163103100000050"); // 20-50%
  } else if (trainConfig == 10){ // EMCAL clusters no timing cut
    cuts.AddCut("51200013","1111100003032230000","0163103100000050"); // 10-20%
    cuts.AddCut("52400013","1111100003032230000","0163103100000050"); // 20-40%
    cuts.AddCut("52500013","1111100003032230000","0163103100000050"); // 20-50%
  } else if (trainConfig == 11){ // EMCAL clusters
    cuts.AddCut("54600013","1111100053032230000","0163103100000050"); // 40-60%
    cuts.AddCut("56800013","1111100053032230000","0163103100000050"); // 60-80%
  } else if (trainConfig == 12){ // EMCAL clusters no timing cut
    cuts.AddCut("54600013","1111100003032230000","0163103100000050"); // 40-60%
    cuts.AddCut("56800013","1111100003032230000","0163103100000050"); // 60-80%

  } else if (trainConfig == 13){ // EMCAL clusters
    cuts.AddCut("50900013","1111100053032230000","0163103100000050"); // 0-90%
  } else if (trainConfig == 14){ // EMCAL clusters no timing cut
    cuts.AddCut("50900013","1111100003032230000","0163103100000050"); // 0-90%

  // EMCal trigger for LHC11h
  } else if (trainConfig == 30){ // EMCAL clusters
    cuts.AddCut("50980013","1111100053032230000","0163103100000050"); // 0-90%
    cuts.AddCut("50180013","1111100053032230000","0163103100000050"); // 0-10%
    cuts.AddCut("51280013","1111100053032230000","0163103100000050"); // 10-20%
    cuts.AddCut("52580013","1111100053032230000","0163103100000050"); // 20-50%
    cuts.AddCut("55880013","1111100053032230000","0163103100000050"); // 50-80%
  } else if (trainConfig == 31){ // EMCAL clusters no timing cut
    cuts.AddCut("50980013","1111100003032230000","0163103100000050"); // 0-90%
    cuts.AddCut("50180013","1111100003032230000","0163103100000050"); // 0-10%
    cuts.AddCut("51280013","1111100003032230000","0163103100000050"); // 10-20%
    cuts.AddCut("52580013","1111100003032230000","0163103100000050"); // 20-50%
    cuts.AddCut("55880013","1111100003032230000","0163103100000050"); // 50-80%
  } else if (trainConfig == 32){ // EMCAL clusters
    cuts.AddCut("50980013","1111100053032230000","0163103100000050"); // 0-90%
  } else if (trainConfig == 33){ // EMCAL clusters no timing cut
    cuts.AddCut("50980013","1111100003032230000","0163103100000050"); // 0-90%
  } else if (trainConfig == 34){ // EMCAL clusters
    cuts.AddCut("50100013","1111102053032230000","01631031000000d0"); // 0-10
    cuts.AddCut("50100013","1111171053032230000","01631031000000d0"); // 0-10
  } else if (trainConfig == 35){ // EMCAL clusters
    cuts.AddCut("50100013","11111020530a2230000","01631031000000d0"); // 0-10 // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("50100013","11111020530a2230000","0163103100000000"); // 0-10 // reproduce Astrids cuts without opening angle cut
    cuts.AddCut("50100013","11111710530a2230000","01631031000000d0"); // 0-10
  } else if (trainConfig == 36){ // EMCAL clusters
    cuts.AddCut("52500013","1111102053032230000","01631031000000d0"); // 20-50
    cuts.AddCut("52500013","1111172053032230000","01631031000000d0"); // 20-50
  } else if (trainConfig == 37){ // EMCAL clusters
    cuts.AddCut("52500013","11111020530a2230000","01631031000000d0"); // 20-50 // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("52500013","11111020530a2230000","0163103100000000"); // 20-50 // reproduce Astrids cuts without opening angle cut
    cuts.AddCut("52500013","11111720530a2230000","01631031000000d0"); // 20-50
  } else if (trainConfig == 38){ // EMCAL clusters - added signals
    cuts.AddCut("50100023","1111171053032230000","01631031000000d0"); // 0-10
    cuts.AddCut("50100023","11111710530a2230000","01631031000000d0"); // 0-10 // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("50100023","11111710530b2230000","01631031000000d0"); // 0-10
  } else if (trainConfig == 39){ // EMCAL clusters - added signals
    cuts.AddCut("50100023","1111171053032230000","01631031000000d0"); // 0-10
    cuts.AddCut("50100023","11111710530a2230000","01631031000000d0"); // 0-10
    cuts.AddCut("50100023","11111710530b2230000","01631031000000d0"); // 0-10
  } else if (trainConfig == 40){ // EMCAL clusters - added signals
    cuts.AddCut("52500023","1111172053032230000","01631031000000d0"); // 20-50
    cuts.AddCut("52500023","11111720530a2230000","01631031000000d0"); // 20-50
    cuts.AddCut("52500023","11111720530b2230000","01631031000000d0"); // 20-50
  } else if (trainConfig == 41){ // EMCAL clusters - added signals
    cuts.AddCut("52500023","1111172053032230000","01631031000000d0"); // 20-50
    cuts.AddCut("52500023","11111720530a2230000","01631031000000d0"); // 20-50
    cuts.AddCut("52500023","11111720530b2230000","01631031000000d0"); // 20-50
  } else if (trainConfig == 42){ // EMCAL clusters - all headers
    cuts.AddCut("50100003","11111020530a2230000","01631031000000d0"); // 0-10   // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("50100003","11111020530a2230000","0163103100000000"); // 0-10   // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("50100003","11111710530a2230000","01631031000000d0"); // 0-10   // new calib
  } else if (trainConfig == 43){ // EMCAL clusters - added signals pi0 forseen
    cuts.AddCut("50100023","11111020530a2230000","01631031000000d0"); // 0-10  // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("50100023","11111020530a2230000","0163103100000000"); // 0-10  // reproduce Astrids cuts without opening angle cut
    cuts.AddCut("50100023","11111710530a2230000","01631031000000d0");
  } else if (trainConfig == 44){ // EMCAL clusters - added signals eta forseen
    cuts.AddCut("50100023","11111020530a2230000","01631031000000d0"); // 0-10  // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("50100023","11111020530a2230000","0163103100000000"); // 0-10  // reproduce Astrids cuts without opening angle cut
    cuts.AddCut("50100023","11111710530a2230000","01631031000000d0");
  } else if (trainConfig == 45){ // EMCAL clusters - added signals other
    cuts.AddCut("50100023","11111020530a2230000","01631031000000d0"); // 0-10  // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("50100023","11111020530a2230000","0163103100000000"); // 0-10  // reproduce Astrids cuts without opening angle cut
    cuts.AddCut("50100023","11111710530a2230000","01631031000000d0");
  } else if (trainConfig == 46){ // EMCAL clusters - added signals other
    cuts.AddCut("50100023","11111020530a2230000","01631031000000d0"); // 0-10  // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("50100023","11111020530a2230000","0163103100000000"); // 0-10  // reproduce Astrids cuts without opening angle cut
    cuts.AddCut("50100023","11111710530a2230000","01631031000000d0");
  } else if (trainConfig == 47){ // EMCAL clusters V0 Cent
    cuts.AddCut("10100013","11111020530a2230000","01631031000000d0"); // 0-10 // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("10100013","11111020530a2230000","0163103100000000"); // 0-10 // reproduce Astrids cuts without opening angle cut
    cuts.AddCut("10100013","11111710530a2230000","01631031000000d0"); // 0-10

  // trainconfig for PbPb studies in 2.76 TeV with TB  nonlin
  } else if (trainConfig == 48){
    cuts.AddCut("10100013","11111020530a2230000","01631031000000d0"); // 0-10 V0M  // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("11200013","11111020530a2230000","01631031000000d0"); // 10-20 V0M // reproduce Astrids cuts with opening angle cut
  } else if (trainConfig == 49){
    cuts.AddCut("30100013","11111020530a2230000","01631031000000d0"); // 0-5 V0M   // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("31200013","11111020530a2230000","01631031000000d0"); // 5-10 V0M  // reproduce Astrids cuts with opening angle cut
  } else if (trainConfig == 50){
    cuts.AddCut("12300013","11111020530a2230000","01631031000000d0"); // 20-30 V0M // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("13400013","11111020530a2230000","01631031000000d0"); // 30-40 V0M // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("14500013","11111020530a2230000","01631031000000d0"); // 40-50 V0M // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("15600013","11111020530a2230000","01631031000000d0"); // 50-60 V0M // reproduce Astrids cuts with opening angle cut
  } else if (trainConfig == 51){
    cuts.AddCut("16700013","11111020530a2230000","01631031000000d0"); // 60-70 V0M // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("17800013","11111020530a2230000","01631031000000d0"); // 70-80 V0M // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("18900013","11111020530a2230000","01631031000000d0"); // 80-90 V0M // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("14600013","11111020530a2230000","01631031000000d0"); // 40-60 V0M // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("16800013","11111020530a2230000","01631031000000d0"); // 60-80 V0M // reproduce Astrids cuts with opening angle cut
    // trainconfig for PbPb studies in 2.76 TeV with no  nonlin
  } else if (trainConfig == 52){
    cuts.AddCut("10100013","11111000530a2230000","01631031000000d0"); // 0-10 V0M  // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("11200013","11111000530a2230000","01631031000000d0"); // 10-20 V0M // reproduce Astrids cuts with opening angle cut
  } else if (trainConfig == 53){
    cuts.AddCut("30100013","11111000530a2230000","01631031000000d0"); // 0-5 V0M   // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("31200013","11111000530a2230000","01631031000000d0"); // 5-10 V0M  // reproduce Astrids cuts with opening angle cut
  } else if (trainConfig == 54){
    cuts.AddCut("12300013","11111000530a2230000","01631031000000d0"); // 20-30 V0M // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("13400013","11111000530a2230000","01631031000000d0"); // 30-40 V0M // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("14500013","11111000530a2230000","01631031000000d0"); // 40-50 V0M // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("15600013","11111000530a2230000","01631031000000d0"); // 50-60 V0M // reproduce Astrids cuts with opening angle cut
  } else if (trainConfig == 55){
    cuts.AddCut("16700013","11111000530a2230000","01631031000000d0"); // 60-70 V0M // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("17800013","11111000530a2230000","01631031000000d0"); // 70-80 V0M // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("18900013","11111000530a2230000","01631031000000d0"); // 80-90 V0M // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("14600013","11111000530a2230000","01631031000000d0"); // 40-60 V0M // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("16800013","11111000530a2230000","01631031000000d0"); // 60-80 V0M // reproduce Astrids cuts with opening angle cut
    // trainconfig for PbPb studies in 2.76 TeV with TB nonlin
  } else if (trainConfig == 56){
    cuts.AddCut("50100013","11111020530a2230000","01631031000000d0"); // 0-10 TM  // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("51200013","11111020530a2230000","01631031000000d0"); // 10-20 TM // reproduce Astrids cuts with opening angle cut
  } else if (trainConfig == 57){
    cuts.AddCut("60100013","11111020530a2230000","01631031000000d0"); // 0-5 TM   // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("61100013","11111020530a2230000","01631031000000d0"); // 5-10 TM  // reproduce Astrids cuts with opening angle cut
  } else if (trainConfig == 58){
    cuts.AddCut("52300013","11111020530a2230000","01631031000000d0"); // 20-30 TM // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("53400013","11111020530a2230000","01631031000000d0"); // 30-40 TM // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("54500013","11111020530a2230000","01631031000000d0"); // 40-50 TM // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("55600013","11111020530a2230000","01631031000000d0"); // 50-60 TM // reproduce Astrids cuts with opening angle cut
  } else if (trainConfig == 59){
    cuts.AddCut("56700013","11111020530a2230000","01631031000000d0"); // 60-70 TM // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("57800013","11111020530a2230000","01631031000000d0"); // 70-80 TM // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("58900013","11111020530a2230000","01631031000000d0"); // 80-90 TM // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("54600013","11111020530a2230000","01631031000000d0"); // 40-60 TM // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("56800013","11111020530a2230000","01631031000000d0"); // 60-80 TM // reproduce Astrids cuts with opening angle cut

  // trainconfig for PbPb studies in 2.76 TeV with TB  nonlin
  } else if (trainConfig == 60){
    cuts.AddCut("10100023","11111020530a2230000","01631031000000d0"); // 0-10 V0M  // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("11200023","11111020530a2230000","01631031000000d0"); // 10-20 V0M // reproduce Astrids cuts with opening angle cut
  } else if (trainConfig == 61){
    cuts.AddCut("30100023","11111020530a2230000","01631031000000d0"); // 0-5 V0M   // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("31200023","11111020530a2230000","01631031000000d0"); // 5-10 V0M  // reproduce Astrids cuts with opening angle cut
  } else if (trainConfig == 62){
    cuts.AddCut("10100023","11111020530a2230000","01631031000000d0"); // 0-10 V0M  // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("11200023","11111020530a2230000","01631031000000d0"); // 10-20 V0M // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("30100023","11111020530a2230000","01631031000000d0"); // 0-5 V0M   // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("31200023","11111020530a2230000","01631031000000d0"); // 5-10 V0M  // reproduce Astrids cuts with opening angle cut
  } else if (trainConfig == 63){
    cuts.AddCut("12300023","11111020530a2230000","01631031000000d0"); // 20-30 V0M // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("13400023","11111020530a2230000","01631031000000d0"); // 30-40 V0M // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("14500023","11111020530a2230000","01631031000000d0"); // 40-50 V0M // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("15600023","11111020530a2230000","01631031000000d0"); // 50-60 V0M // reproduce Astrids cuts with opening angle cut
  } else if (trainConfig == 64){
    cuts.AddCut("12300023","11111020530a2230000","01631031000000d0"); // 20-30 V0M // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("13400023","11111020530a2230000","01631031000000d0"); // 30-40 V0M // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("14500023","11111020530a2230000","01631031000000d0"); // 40-50 V0M // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("15600023","11111020530a2230000","01631031000000d0"); // 50-60 V0M // reproduce Astrids cuts with opening angle cut
  } else if (trainConfig == 65){
    cuts.AddCut("16700023","11111020530a2230000","01631031000000d0"); // 60-70 V0M // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("17800023","11111020530a2230000","01631031000000d0"); // 70-80 V0M // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("18900023","11111020530a2230000","01631031000000d0"); // 80-90 V0M // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("14600023","11111020530a2230000","01631031000000d0"); // 40-60 V0M // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("16800023","11111020530a2230000","01631031000000d0"); // 60-80 V0M // reproduce Astrids cuts with opening angle cut
  } else if (trainConfig == 66){
    cuts.AddCut("16700023","11111020530a2230000","01631031000000d0"); // 60-70 V0M // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("17800023","11111020530a2230000","01631031000000d0"); // 70-80 V0M // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("18900023","11111020530a2230000","01631031000000d0"); // 80-90 V0M // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("14600023","11111020530a2230000","01631031000000d0"); // 40-60 V0M // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("16800023","11111020530a2230000","01631031000000d0"); // 60-80 V0M // reproduce Astrids cuts with opening angle cut

  } else if (trainConfig == 67){
    cuts.AddCut("40100013","11111020530a2230000","01631031000000d0"); // 0-10 new track array cuts in MC  // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("41200013","11111020530a2230000","01631031000000d0"); // 10-20 new track array cuts in MC // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("42300013","11111020530a2230000","01631031000000d0"); // 20-30 new track array cuts in MC // reproduce Astrids cuts with opening angle cut
  } else if (trainConfig == 68){
    cuts.AddCut("70100013","11111020530a2230000","01631031000000d0"); // 0-5 new track array cuts in MC   // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("71200013","11111020530a2230000","01631031000000d0"); // 5-10 new track array cuts in MC  // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("72300013","11111020530a2230000","01631031000000d0"); // 10-15 new track array cuts in MC  // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("73400013","11111020530a2230000","01631031000000d0"); // 15-20 new track array cuts in MC  // reproduce Astrids cuts with opening angle cut
  } else if (trainConfig == 69){
    cuts.AddCut("43400013","11111020530a2230000","01631031000000d0"); // 30-40 new track array cuts in MC  // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("44500013","11111020530a2230000","01631031000000d0"); // 40-50 new track array cuts in MC // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("45600013","11111020530a2230000","01631031000000d0"); // 50-60 new track array cuts in MC // reproduce Astrids cuts with opening angle cut
  } else if (trainConfig == 70){
    cuts.AddCut("46700013","11111020530a2230000","01631031000000d0"); // 60-70 new track array cuts in MC  // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("47800013","11111020530a2230000","01631031000000d0"); // 70-80 new track array cuts in MC // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("48900013","11111020530a2230000","01631031000000d0"); // 80-90 new track array cuts in MC // reproduce Astrids cuts with opening angle cut

  } else if (trainConfig == 71){ // added signals
    cuts.AddCut("40100023","11111020530a2230000","01631031000000d0"); // 0-10 new track array cuts in MC  // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("41200023","11111020530a2230000","01631031000000d0"); // 10-20 new track array cuts in MC // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("42300023","11111020530a2230000","01631031000000d0"); // 20-30 new track array cuts in MC // reproduce Astrids cuts with opening angle cut
  } else if (trainConfig == 72){ // added signals
    cuts.AddCut("70100023","11111020530a2230000","01631031000000d0"); // 0-5 new track array cuts in MC   // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("71200023","11111020530a2230000","01631031000000d0"); // 5-10 new track array cuts in MC  // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("72300023","11111020530a2230000","01631031000000d0"); // 10-15 new track array cuts in MC  // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("73400023","11111020530a2230000","01631031000000d0"); // 15-20 new track array cuts in MC  // reproduce Astrids cuts with opening angle cut
  } else if (trainConfig == 73){ // added signals
    cuts.AddCut("43400023","11111020530a2230000","01631031000000d0"); // 30-40 new track array cuts in MC  // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("44500023","11111020530a2230000","01631031000000d0"); // 40-50 new track array cuts in MC // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("45600023","11111020530a2230000","01631031000000d0"); // 60-70 new track array cuts in MC // reproduce Astrids cuts with opening angle cut
  } else if (trainConfig == 74){ // added signals
    cuts.AddCut("46700023","11111020530a2230000","01631031000000d0"); // 60-70 new track array cuts in MC  // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("47800023","11111020530a2230000","01631031000000d0"); // 70-80 new track array cuts in MC // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("48900023","11111020530a2230000","01631031000000d0"); // 80-90 new track array cuts in MC // reproduce Astrids cuts with opening angle cut
  } else if (trainConfig == 75){ // added signals duplicate
    cuts.AddCut("40100023","11111020530a2230000","01631031000000d0"); // 0-10 new track array cuts in MC  // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("41200023","11111020530a2230000","01631031000000d0"); // 10-20 new track array cuts in MC // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("42300023","11111020530a2230000","01631031000000d0"); // 20-30 new track array cuts in MC // reproduce Astrids cuts with opening angle cut
  } else if (trainConfig == 76){ // added signals duplicate
    cuts.AddCut("70100023","11111020530a2230000","01631031000000d0"); // 0-5 new track array cuts in MC   // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("71200023","11111020530a2230000","01631031000000d0"); // 5-10 new track array cuts in MC  // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("72300023","11111020530a2230000","01631031000000d0"); // 10-15 new track array cuts in MC  // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("73400023","11111020530a2230000","01631031000000d0"); // 15-20 new track array cuts in MC  // reproduce Astrids cuts with opening angle cut
  } else if (trainConfig == 78){ // added signals duplicate
    cuts.AddCut("43400023","11111020530a2230000","01631031000000d0"); // 30-40 new track array cuts in MC  // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("44500023","11111020530a2230000","01631031000000d0"); // 40-50 new track array cuts in MC // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("45600023","11111020530a2230000","01631031000000d0"); // 60-70 new track array cuts in MC // reproduce Astrids cuts with opening angle cut
  } else if (trainConfig == 79){ // added signals duplicate
    cuts.AddCut("46700023","11111020530a2230000","01631031000000d0"); // 60-70 new track array cuts in MC  // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("47800023","11111020530a2230000","01631031000000d0"); // 70-80 new track array cuts in MC // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("48900023","11111020530a2230000","01631031000000d0"); // 80-90 new track array cuts in MC // reproduce Astrids cuts with opening angle cut


  // **********************************************************************************************************
  // ***************************** PHOS configurations run 1 **************************************************
  // **********************************************************************************************************
  } else if (trainConfig == 101){ // PHOS clusters
    cuts.AddCut("60100013","2444400040033200000","0163103100000030"); // 0-5%
    cuts.AddCut("61200013","2444400040033200000","0163103100000030"); // 5-10%
    cuts.AddCut("50100013","2444400040033200000","0163103100000030"); // 0-10%
    cuts.AddCut("52400013","2444400040033200000","0163103100000030"); // 20-40%
    cuts.AddCut("52500013","2444400040033200000","0163103100000030"); // 20-50%
  } else if (trainConfig == 102){ // PHOS clusters
    cuts.AddCut("60100013","2444400040033200000","0163103100000030"); // 0-5%
    cuts.AddCut("61200013","2444400040033200000","0163103100000030"); // 5-10%
    cuts.AddCut("50100013","2444400040033200000","0163103100000030"); // 0-10%
    cuts.AddCut("51200013","2444400040033200000","0163103100000030"); // 10-20%
    cuts.AddCut("52400013","2444400040033200000","0163103100000030"); // 20-40%
  } else if (trainConfig == 103){ // PHOS clusters
    cuts.AddCut("54600013","2444400040033200000","0163103100000030"); // 40-60%
    cuts.AddCut("56800013","2444400040033200000","0163103100000030"); // 60-80%
    cuts.AddCut("52600013","2444400040033200000","0163103100000030"); // 20-60%
    cuts.AddCut("54800013","2444400040033200000","0163103100000030"); // 40-80%
    cuts.AddCut("52500013","2444400040033200000","0163103100000030"); // 20-50%
  } else if (trainConfig == 104){ // PHOS clusters with TM
    cuts.AddCut("60100013","2444400042033200000","0163103100000030"); // 0-5%
    cuts.AddCut("61200013","2444400042033200000","0163103100000030"); // 5-10%
    cuts.AddCut("50100013","2444400042033200000","0163103100000030"); // 0-10%
    cuts.AddCut("52400013","2444400042033200000","0163103100000030"); // 20-40%
    cuts.AddCut("52500013","2444400042033200000","0163103100000030"); // 20-50%
  } else if (trainConfig == 105){ // PHOS clusters with TM
    cuts.AddCut("60100013","2444400042033200000","0163103100000030"); // 0-5%
    cuts.AddCut("61200013","2444400042033200000","0163103100000030"); // 5-10%
    cuts.AddCut("50100013","2444400042033200000","0163103100000030"); // 0-10%
    cuts.AddCut("51200013","2444400042033200000","0163103100000030"); // 10-20%
    cuts.AddCut("52400013","2444400042033200000","0163103100000030"); // 20-40%
  } else if (trainConfig == 106){ // PHOS clusters with TM
    cuts.AddCut("54600013","2444400042033200000","0163103100000030"); // 40-60%
    cuts.AddCut("56800013","2444400042033200000","0163103100000030"); // 60-80%
    cuts.AddCut("52600013","2444400042033200000","0163103100000030"); // 20-60%
    cuts.AddCut("54800013","2444400042033200000","0163103100000030"); // 40-80%
    cuts.AddCut("52500013","2444400042033200000","0163103100000030"); // 20-50%
  } else if (trainConfig == 107){ // PHOS clusters with TM
    cuts.AddCut("50900013","2444400002033200000","0163103100000030"); // 0-90%


  // **********************************************************************************************************
  // ***************************** EMC configurations PbPb run 2 **********************************************
  // **********************************************************************************************************
  } else if (trainConfig == 201){ // EMCAL clusters central
    cuts.AddCut("10110013","1111100057032230000","01631031000000d0"); // 0-10
    cuts.AddCut("30110013","1111100057032230000","01631031000000d0"); // 0-5
    cuts.AddCut("31210013","1111100057032230000","01631031000000d0"); // 5-10
  } else if (trainConfig == 202){ // EMCAL clusters semi-central
    cuts.AddCut("11210013","1111100057032230000","01631031000000d0"); // 10-20
    cuts.AddCut("12310013","1111100057032230000","01631031000000d0"); // 20-30
    cuts.AddCut("13410013","1111100057032230000","01631031000000d0"); // 30-40
    cuts.AddCut("12410013","1111100057032230000","01631031000000d0"); // 20-40
  } else if (trainConfig == 203){ // EMCAL clusters semi-central
    cuts.AddCut("14510013","1111100057032230000","01631031000000d0"); // 40-50
    cuts.AddCut("14610013","1111100057032230000","01631031000000d0"); // 40-60
    cuts.AddCut("15610013","1111100057032230000","01631031000000d0"); // 50-60
  } else if (trainConfig == 204){ // EMCAL clusters peripheral
    cuts.AddCut("16710013","1111100057032230000","01631031000000d0"); // 60-70
    cuts.AddCut("17810013","1111100057032230000","01631031000000d0"); // 70-80
    cuts.AddCut("18910013","1111100057032230000","01631031000000d0"); // 80-90
    cuts.AddCut("16810013","1111100057032230000","01631031000000d0"); // 60-80

  } else if (trainConfig == 205){ // EMCAL clusters central add sig
    cuts.AddCut("10110023","1111100057032230000","01631031000000d0"); // 0-10
    cuts.AddCut("30110023","1111100057032230000","01631031000000d0"); // 0-5
    cuts.AddCut("31210023","1111100057032230000","01631031000000d0"); // 5-10
  } else if (trainConfig == 206){ // EMCAL clusters semi-central add sig
    cuts.AddCut("11210023","1111100057032230000","01631031000000d0"); // 10-20
    cuts.AddCut("12310023","1111100057032230000","01631031000000d0"); // 20-30
    cuts.AddCut("13410023","1111100057032230000","01631031000000d0"); // 30-40
    cuts.AddCut("12410023","1111100057032230000","01631031000000d0"); // 20-40
  } else if (trainConfig == 207){ // EMCAL clusters semi-central add sig
    cuts.AddCut("14510023","1111100057032230000","01631031000000d0"); // 40-50
    cuts.AddCut("14610023","1111100057032230000","01631031000000d0"); // 40-60
    cuts.AddCut("15610023","1111100057032230000","01631031000000d0"); // 50-60
  } else if (trainConfig == 208){ // EMCAL clusters peripheral add sig
    cuts.AddCut("16710023","1111100057032230000","01631031000000d0"); // 60-70
    cuts.AddCut("17810023","1111100057032230000","01631031000000d0"); // 70-80
    cuts.AddCut("18910023","1111100057032230000","01631031000000d0"); // 80-90
    cuts.AddCut("16810023","1111100057032230000","01631031000000d0"); // 60-80

  } else if (trainConfig == 209){ // EMCAL clusters - correction convcalo f1
    cuts.AddCut("10110013","1111181053032230000","0163103100000050"); // 0-10
    cuts.AddCut("11210013","1111181053032230000","0163103100000050"); // 10-20
    cuts.AddCut("12510013","1111181053032230000","0163103100000050"); // 20-50
    cuts.AddCut("15910013","1111181053032230000","0163103100000050"); // 50-90
  } else if (trainConfig == 210){ // EMCAL clusters - correction calocalo f2
    cuts.AddCut("10110013","1111192053032230000","0163103100000050"); // 0-10
    cuts.AddCut("11210013","1111192053032230000","0163103100000050"); // 10-20
    cuts.AddCut("12510013","1111192053032230000","0163103100000050"); // 20-50
    cuts.AddCut("15910013","1111192053032230000","0163103100000050"); // 50-90
  } else if (trainConfig == 211){ // EMCAL clusters - 0-90% centrality for PbPb EMCal cluster QA
    cuts.AddCut("10910013","1111100003032230000","0163103100000050"); // 0-90

  } else if (trainConfig == 212){ // EMCAL clusters - 0-90% centrality for PbPb EMCal
    cuts.AddCut("50910113","1111183053032230000","0163103100000050"); //  0-90 calo correction cent dep
  } else if (trainConfig == 213){ // EMCAL clusters - centrality selection for PbPb EMCal
    cuts.AddCut("50110113","1111184053032230000","0163103100000050"); //  0-10 calo correction cent dep
    cuts.AddCut("51210113","1111185053032230000","0163103100000050"); // 10-20 calo correction cent dep
    cuts.AddCut("52510113","1111186053032230000","0163103100000050"); // 20-50 calo correction cent dep
    cuts.AddCut("55910113","1111187053032230000","0163103100000050"); // 50-90 calo correction cent dep

  // trainconfig for PbPb studies in 5 TeV with TB nonlin
  } else if (trainConfig == 214){
    cuts.AddCut("10110113","11111020530a2230000","01631031000000d0"); // 0-10 V0M  // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("11210113","11111020530a2230000","01631031000000d0"); // 10-20 V0M // reproduce Astrids cuts with opening angle cut
  } else if (trainConfig == 215){
    cuts.AddCut("30110113","11111020530a2230000","01631031000000d0"); // 0-5 V0M   // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("31210113","11111020530a2230000","01631031000000d0"); // 5-10 V0M  // reproduce Astrids cuts with opening angle cut
  } else if (trainConfig == 216){
    cuts.AddCut("12310113","11111020530a2230000","01631031000000d0"); // 20-30 V0M // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("13410113","11111020530a2230000","01631031000000d0"); // 30-40 V0M // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("14510113","11111020530a2230000","01631031000000d0"); // 40-50 V0M // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("15610113","11111020530a2230000","01631031000000d0"); // 50-60 V0M // reproduce Astrids cuts with opening angle cut
  } else if (trainConfig == 217){
    cuts.AddCut("16710113","11111020530a2230000","01631031000000d0"); // 60-70 V0M // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("17810113","11111020530a2230000","01631031000000d0"); // 70-80 V0M // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("18910113","11111020530a2230000","01631031000000d0"); // 80-90 V0M // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("14610113","11111020530a2230000","01631031000000d0"); // 40-60 V0M // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("16810113","11111020530a2230000","01631031000000d0"); // 60-80 V0M // reproduce Astrids cuts with opening angle cut
  // trainconfig for PbPb studies in 5 TeV with no  nonlin
  } else if (trainConfig == 218){
    cuts.AddCut("10110113","11111000530a2230000","01631031000000d0"); // 0-10 V0M  // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("11210113","11111000530a2230000","01631031000000d0"); // 10-20 V0M // reproduce Astrids cuts with opening angle cut
  } else if (trainConfig == 219){
    cuts.AddCut("30110113","11111000530a2230000","01631031000000d0"); // 0-5 V0M   // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("31210113","11111000530a2230000","01631031000000d0"); // 5-10 V0M  // reproduce Astrids cuts with opening angle cut
  } else if (trainConfig == 220){
    cuts.AddCut("12310113","11111000530a2230000","01631031000000d0"); // 20-30 V0M // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("13410113","11111000530a2230000","01631031000000d0"); // 30-40 V0M // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("14510113","11111000530a2230000","01631031000000d0"); // 40-50 V0M // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("15610113","11111000530a2230000","01631031000000d0"); // 50-60 V0M // reproduce Astrids cuts with opening angle cut
  } else if (trainConfig == 221){
    cuts.AddCut("16710113","11111000530a2230000","01631031000000d0"); // 60-70 V0M // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("17810113","11111000530a2230000","01631031000000d0"); // 70-80 V0M // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("18910113","11111000530a2230000","01631031000000d0"); // 80-90 V0M // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("14610113","11111000530a2230000","01631031000000d0"); // 40-60 V0M // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("16810113","11111000530a2230000","01631031000000d0"); // 60-80 V0M // reproduce Astrids cuts with opening angle cut
    // trainconfig for PbPb studies in 5 TeV with TB nonlin
  } else if (trainConfig == 222){
    cuts.AddCut("50110113","11111020530a2230000","01631031000000d0"); // 0-10 TM  // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("51210113","11111020530a2230000","01631031000000d0"); // 10-20 TM // reproduce Astrids cuts with opening angle cut
  } else if (trainConfig == 223){
    cuts.AddCut("60110113","11111020530a2230000","01631031000000d0"); // 0-5 TM   // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("61110113","11111020530a2230000","01631031000000d0"); // 5-10 TM  // reproduce Astrids cuts with opening angle cut
  } else if (trainConfig == 224){
    cuts.AddCut("52310113","11111020530a2230000","01631031000000d0"); // 20-30 TM // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("53410113","11111020530a2230000","01631031000000d0"); // 30-40 TM // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("54510113","11111020530a2230000","01631031000000d0"); // 40-50 TM // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("55610113","11111020530a2230000","01631031000000d0"); // 50-60 TM // reproduce Astrids cuts with opening angle cut
  } else if (trainConfig == 225){
    cuts.AddCut("56710113","11111020530a2230000","01631031000000d0"); // 60-70 TM // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("57810113","11111020530a2230000","01631031000000d0"); // 70-80 TM // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("58910113","11111020530a2230000","01631031000000d0"); // 80-90 TM // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("54610113","11111020530a2230000","01631031000000d0"); // 40-60 TM // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("56810113","11111020530a2230000","01631031000000d0"); // 60-80 TM // reproduce Astrids cuts with opening angle cut

  } else if (trainConfig == 226){ // EMCAL clusters - peripheral centrality selection for PbPb EMCal
    cuts.AddCut("15910113","11111020530a2230000","01631031000000d0"); //
    cuts.AddCut("15910113","11111870530a2230000","01631031000000d0"); //
    cuts.AddCut("15910113","11111020530b2230000","01631031000000d0"); //
    cuts.AddCut("15910113","11111870530b2230000","01631031000000d0"); //

  } else if (trainConfig == 230){ // EMCAL clusters central add sig
    cuts.AddCut("10110023","1111100057032230000","01631031000000d0"); // 0-10
    cuts.AddCut("30110023","1111100057032230000","01631031000000d0"); // 0-5
    cuts.AddCut("31210023","1111100057032230000","01631031000000d0"); // 5-10
  } else if (trainConfig == 231){ // EMCAL clusters semi-central add sig
    cuts.AddCut("11210023","1111100057032230000","01631031000000d0"); // 10-20
    cuts.AddCut("12310023","1111100057032230000","01631031000000d0"); // 20-30
    cuts.AddCut("13410023","1111100057032230000","01631031000000d0"); // 30-40
    cuts.AddCut("12410023","1111100057032230000","01631031000000d0"); // 20-40
  } else if (trainConfig == 232){ // EMCAL clusters semi-central add sig
    cuts.AddCut("14510023","1111100057032230000","01631031000000d0"); // 40-50
    cuts.AddCut("14610023","1111100057032230000","01631031000000d0"); // 40-60
    cuts.AddCut("15610023","1111100057032230000","01631031000000d0"); // 50-60
  } else if (trainConfig == 233){ // EMCAL clusters peripheral add sig
    cuts.AddCut("16710023","1111100057032230000","01631031000000d0"); // 60-70
    cuts.AddCut("17810023","1111100057032230000","01631031000000d0"); // 70-80
    cuts.AddCut("18910023","1111100057032230000","01631031000000d0"); // 80-90
    cuts.AddCut("16810023","1111100057032230000","01631031000000d0"); // 60-80


  } else if (trainConfig == 234){ // EMCAL clusters central, TB
    cuts.AddCut("10110113","1111102057032230000","01631031000000d0"); // 0-10
    cuts.AddCut("30110113","1111102057032230000","01631031000000d0"); // 0-5
    cuts.AddCut("31210113","1111102057032230000","01631031000000d0"); // 5-10
  } else if (trainConfig == 235){ // EMCAL clusters semi-central, TB
    cuts.AddCut("11210113","1111102057032230000","01631031000000d0"); // 10-20
    cuts.AddCut("12310113","1111102057032230000","01631031000000d0"); // 20-30
    cuts.AddCut("13410113","1111102057032230000","01631031000000d0"); // 30-40
    cuts.AddCut("12410113","1111102057032230000","01631031000000d0"); // 20-40
  } else if (trainConfig == 236){ // EMCAL clusters semi-central, TB
    cuts.AddCut("14510113","1111102057032230000","01631031000000d0"); // 40-50
    cuts.AddCut("14610113","1111102057032230000","01631031000000d0"); // 40-60
    cuts.AddCut("15610113","1111102057032230000","01631031000000d0"); // 50-60
  } else if (trainConfig == 237){ // EMCAL clusters peripheral, TB
    cuts.AddCut("16710113","1111102057032230000","01631031000000d0"); // 60-70
    cuts.AddCut("17810113","1111102057032230000","01631031000000d0"); // 70-80
    cuts.AddCut("18910113","1111102057032230000","01631031000000d0"); // 80-90
    cuts.AddCut("16810113","1111102057032230000","01631031000000d0"); // 60-80


  } else if (trainConfig == 240){ // EMCAL clusters - 0-90% centrality for PbPb EMCal
    cuts.AddCut("50910013","1111100053032230000","0163103100000050"); //  0-90 calo correction cent dep
    cuts.AddCut("50910013","1111102053032230000","0163103100000050"); //  0-90 calo correction cent dep
    cuts.AddCut("50910013","1111183053032230000","0163103100000050"); //  0-90 calo correction cent dep
  } else if (trainConfig == 241){ // EMCAL clusters - 0-90% centrality for PbPb EMCal
    cuts.AddCut("50910613","1111100053032230000","0163103100000050"); //  0-90 calo correction cent dep
    cuts.AddCut("50910613","1111102053032230000","0163103100000050"); //  0-90 calo correction cent dep
    cuts.AddCut("50910613","1111183053032230000","0163103100000050"); //  0-90 calo correction cent dep
  } else if (trainConfig == 242){ // EMCAL clusters - centrality selection for PbPb EMCal
    cuts.AddCut("50110013","1111100053032230000","0163103100000050"); //  0-10 calo correction cent dep
    cuts.AddCut("51210013","1111100053032230000","0163103100000050"); // 10-20 calo correction cent dep
    cuts.AddCut("52510013","1111100053032230000","0163103100000050"); // 20-50 calo correction cent dep
    cuts.AddCut("55910013","1111100053032230000","0163103100000050"); // 50-90 calo correction cent dep
  } else if (trainConfig == 243){ // EMCAL clusters - centrality selection for PbPb EMCal
    cuts.AddCut("50110613","1111100053032230000","0163103100000050"); //  0-10 calo correction cent dep
    cuts.AddCut("51210613","1111100053032230000","0163103100000050"); // 10-20 calo correction cent dep
    cuts.AddCut("52510613","1111100053032230000","0163103100000050"); // 20-50 calo correction cent dep
    cuts.AddCut("55910613","1111100053032230000","0163103100000050"); // 50-90 calo correction cent dep
  } else if (trainConfig == 244){ // EMCAL clusters - centrality selection for PbPb EMCal
    cuts.AddCut("50110613","1111102053032230000","0163103100000050"); //  0-10 calo correction cent dep
    cuts.AddCut("51210613","1111102053032230000","0163103100000050"); // 10-20 calo correction cent dep
    cuts.AddCut("52510613","1111102053032230000","0163103100000050"); // 20-50 calo correction cent dep
    cuts.AddCut("55910613","1111102053032230000","0163103100000050"); // 50-90 calo correction cent dep
  } else if (trainConfig == 245){ // EMCAL clusters - centrality selection for PbPb EMCal
    cuts.AddCut("50110613","1111184053032230000","0163103100000050"); //  0-10 calo correction cent dep
    cuts.AddCut("51210613","1111185053032230000","0163103100000050"); // 10-20 calo correction cent dep
    cuts.AddCut("52510613","1111186053032230000","0163103100000050"); // 20-50 calo correction cent dep
    cuts.AddCut("55910613","1111187053032230000","0163103100000050"); // 50-90 calo correction cent dep
  } else if (trainConfig == 246){ // EMCAL clusters - 0-90% centrality for PbPb EMCal
    cuts.AddCut("50910613","1111183050032230000","0163103100000050"); //  0-90 calo correction cent dep
    cuts.AddCut("50910613","1111183051032230000","0163103100000050"); //  0-90 calo correction cent dep
    cuts.AddCut("50910613","1111183057032230000","0163103100000050"); //  0-90 calo correction cent dep
  } else if (trainConfig == 247){ // EMCAL clusters - centrality selection for PbPb EMCal
    cuts.AddCut("50110613","1111184050032230000","0163103100000050"); //  0-10 calo correction cent dep
    cuts.AddCut("51210613","1111185050032230000","0163103100000050"); // 10-20 calo correction cent dep
    cuts.AddCut("52510613","1111186050032230000","0163103100000050"); // 20-50 calo correction cent dep
    cuts.AddCut("55910613","1111187050032230000","0163103100000050"); // 50-90 calo correction cent dep
  } else if (trainConfig == 248){ // EMCAL clusters - centrality selection for PbPb EMCal
    cuts.AddCut("50110613","1111184051032230000","0163103100000050"); //  0-10 calo correction cent dep
    cuts.AddCut("51210613","1111185051032230000","0163103100000050"); // 10-20 calo correction cent dep
    cuts.AddCut("52510613","1111186051032230000","0163103100000050"); // 20-50 calo correction cent dep
    cuts.AddCut("55910613","1111187051032230000","0163103100000050"); // 50-90 calo correction cent dep
  } else if (trainConfig == 249){ // EMCAL clusters - centrality selection for PbPb EMCal
    cuts.AddCut("50110613","1111184057032230000","0163103100000050"); //  0-10 calo correction cent dep
    cuts.AddCut("51210613","1111185057032230000","0163103100000050"); // 10-20 calo correction cent dep
    cuts.AddCut("52510613","1111186057032230000","0163103100000050"); // 20-50 calo correction cent dep
    cuts.AddCut("55910613","1111187057032230000","0163103100000050"); // 50-90 calo correction cent dep
  } else if (trainConfig == 250){ // EMCAL clusters - centrality selection for PbPb EMCal
    cuts.AddCut("10910013","1111183051032230000","0163103100000050"); //  0-90 calo correction
    cuts.AddCut("10910a13","1111183051032230000","0163103100000050"); //  0-90 calo correction
  } else if (trainConfig == 251){ // EMCAL clusters - centrality selection for PbPb EMCal
    cuts.AddCut("10110013","1111183051032230000","0163103100000050"); //  0-90 calo correction
    cuts.AddCut("11210013","1111183051032230000","0163103100000050"); //  0-90 calo correction
    cuts.AddCut("12510013","1111183051032230000","0163103100000050"); //  0-90 calo correction
    cuts.AddCut("15910013","1111183051032230000","0163103100000050"); //  0-90 calo correction
  } else if (trainConfig == 252){ // EMCAL clusters - centrality selection for PbPb EMCal
    cuts.AddCut("10110a13","1111183051032230000","0163103100000050"); //  0-90 calo correction
    cuts.AddCut("11210a13","1111183051032230000","0163103100000050"); //  0-90 calo correction
    cuts.AddCut("12510a13","1111183051032230000","0163103100000050"); //  0-90 calo correction
    cuts.AddCut("15910a13","1111183051032230000","0163103100000050"); //  0-90 calo correction
  } else if (trainConfig == 253){ // EMCAL clusters - 20180718 - default without corrections
    cuts.AddCut("10110a13","1111100051032230000","0163103100000050"); //  0-90 calo correction
    cuts.AddCut("11210a13","1111100051032230000","0163103100000050"); //  0-90 calo correction
    cuts.AddCut("12410a13","1111100051032230000","0163103100000050"); //  0-90 calo correction
    cuts.AddCut("14610a13","1111100051032230000","0163103100000050"); //  0-90 calo correction
    cuts.AddCut("16810a13","1111100051032230000","0163103100000050"); //  0-90 calo correction
  } else if (trainConfig == 254){ // EMCAL clusters - 20180718 - default with peripheral corrections
    cuts.AddCut("10110a13","1111187051032230000","0163103100000050"); //
    cuts.AddCut("11210a13","1111187051032230000","0163103100000050"); //
    cuts.AddCut("12410a13","1111187051032230000","0163103100000050"); //
    cuts.AddCut("14610a13","1111187051032230000","0163103100000050"); //
    cuts.AddCut("16810a13","1111187051032230000","0163103100000050"); //



  // **********************************************************************************************************
  // ***************************** EMC configurations XeXe run 2 **********************************************
  // **********************************************************************************************************
  } else if (trainConfig == 300){ // EMCAL clusters - 0-80% centrality for XeXe EMCal cluster QA
    cuts.AddCut("10810013","1111100007032230000","01631031000000d0"); // 0-80
  } else if (trainConfig == 301){ // EMCAL clusters - 0-80% centrality for XeXe EMCal cluster QA
    cuts.AddCut("10210013","1111100007032230000","01631031000000d0"); // 0-20
    cuts.AddCut("12410013","1111100007032230000","01631031000000d0"); // 20-40
    cuts.AddCut("10410013","1111100007032230000","01631031000000d0"); // 0-40
    cuts.AddCut("14810013","1111100007032230000","01631031000000d0"); // 40-80
  } else if (trainConfig == 302){ // EMCAL clusters - 0-80% centrality for XeXe EMCal cluster QA - TB nl
    cuts.AddCut("10810013","1111102007032230000","01631031000000d0"); // 0-80
  } else if (trainConfig == 303){ // EMCAL clusters - diff centralities for XeXe EMCal cluster QA - TB nl
    cuts.AddCut("10210013","1111102007032230000","01631031000000d0"); // 0-20
    cuts.AddCut("12410013","1111102007032230000","01631031000000d0"); // 20-40
    cuts.AddCut("10410013","1111102007032230000","01631031000000d0"); // 0-40
    cuts.AddCut("14810013","1111102007032230000","01631031000000d0"); // 40-80

  //-----------------------------------------------------------------------------------------------
  // Systematics variations 0-20 EMC
  //-----------------------------------------------------------------------------------------------
  } else if (trainConfig == 310){ // EMCal clusters non lin variations 0-20 %
    cuts.AddCut("10210013","1111171007032230000","01631031000000d0");
    cuts.AddCut("10210013","1111172007032230000","01631031000000d0");
    cuts.AddCut("10210013","1111181007032230000","01631031000000d0");
    cuts.AddCut("10210013","1111182007032230000","01631031000000d0");
  } else if (trainConfig == 311){ // second set of variations CLUSTER
    cuts.AddCut("10210013","1111102007022230000","01631031000000d0"); // min energy cluster variation 1  600 MeV
    cuts.AddCut("10210013","1111102007042230000","01631031000000d0"); // min energy cluster variation 2  800 MeV
    cuts.AddCut("10210013","1111102007052230000","01631031000000d0"); // min energy cluster variation 3  900 MeV
    cuts.AddCut("10210013","1111102007032230000","01631031000000d0"); // min/max M02  0.1<M<0.5
    cuts.AddCut("10210013","1111102007032200000","01631031000000d0"); // min/max M02  0.1<M<100
  } else if (trainConfig == 312){ // third set of variations CLUSTER
    cuts.AddCut("10210013","1111102007032250000","01631031000000d0"); // min/max M02  0.1<M<0.3
    cuts.AddCut("10210013","1111102007032260000","01631031000000d0"); // min/max M02  0.1<M<0.27
    cuts.AddCut("10210013","1111102007031230000","01631031000000d0"); // min number of cells variation 1  1 cell
    cuts.AddCut("10210013","1112102007032230000","01631031000000d0"); // only modules with TRD infront
    cuts.AddCut("10210013","1111302007032230000","01631031000000d0"); // no modules with TRD infront
  } else if (trainConfig == 313){ // third set of variations MESON
    cuts.AddCut("10210013","1111102007032230000","01633031000000d0"); // rapidity variation  y<0.6
    cuts.AddCut("10210013","1111102007032230000","01634031000000d0"); // rapidity variation  y<0.5
    cuts.AddCut("10210013","1111102007032230000","01631061000000d0"); // alpha meson variation 1   0<alpha<0.8
    cuts.AddCut("10210013","1111102007032230000","01631051000000d0"); // alpha meson variation 2  0<alpha<0.75
  } else if (trainConfig == 314){ // opening angle variations
    cuts.AddCut("10210013","1111102007032230000","0163103100000040"); // min opening angle 0.0152
    cuts.AddCut("10210013","1111102007032230000","0163103100000050"); // min opening angle 0.0202
    cuts.AddCut("10210013","1111102007032230000","0163103100000070"); // min opening angle 0.016
    cuts.AddCut("10210013","1111102007032230000","0163103100000080"); // min opening angle 0.018
    cuts.AddCut("10210013","1111102007032230000","0163103100000090"); // min opening angle 0.018
  } else if (trainConfig == 315){ // TM variations
    cuts.AddCut("10210013","1111102003032230000","01631031000000d0"); // fixed window
    cuts.AddCut("10210013","1111102006032230000","01631031000000d0"); // tm pt dependent var 1
    cuts.AddCut("10210013","1111102008032230000","01631031000000d0"); // tm pt dependent var 2
    cuts.AddCut("10210013","1111102009032230000","01631031000000d0"); // tm pt dependent var 3

  //-----------------------------------------------------------------------------------------------
  // Systematics variations 20-40 EMC
  //-----------------------------------------------------------------------------------------------
  } else if (trainConfig == 320){ // EMCal clusters non lin variations 20-40 %
    cuts.AddCut("12410013","1111171007032230000","01631031000000d0");
    cuts.AddCut("12410013","1111172007032230000","01631031000000d0");
    cuts.AddCut("12410013","1111181007032230000","01631031000000d0");
    cuts.AddCut("12410013","1111182007032230000","01631031000000d0");
  } else if (trainConfig == 321){ // second set of variations CLUSTER
    cuts.AddCut("12410013","1111102007022230000","01631031000000d0"); // min energy cluster variation 1  600 MeV
    cuts.AddCut("12410013","1111102007042230000","01631031000000d0"); // min energy cluster variation 2  800 MeV
    cuts.AddCut("12410013","1111102007052230000","01631031000000d0"); // min energy cluster variation 3  900 MeV
    cuts.AddCut("12410013","1111102007032230000","01631031000000d0"); // min/max M02  0.1<M<0.5
    cuts.AddCut("12410013","1111102007032200000","01631031000000d0"); // min/max M02  0.1<M<100
  } else if (trainConfig == 322){ // third set of variations CLUSTER
    cuts.AddCut("12410013","1111102007032250000","01631031000000d0"); // min/max M02  0.1<M<0.3
    cuts.AddCut("12410013","1111102007032260000","01631031000000d0"); // min/max M02  0.1<M<0.27
    cuts.AddCut("12410013","1111102007031230000","01631031000000d0"); // min number of cells variation 1  1 cell
    cuts.AddCut("12410013","1112102007032230000","01631031000000d0"); // only modules with TRD infront
    cuts.AddCut("12410013","1111302007032230000","01631031000000d0"); // no modules with TRD infront
  } else if (trainConfig == 323){ // third set of variations MESON
    cuts.AddCut("12410013","1111102007032230000","01633031000000d0"); // rapidity variation  y<0.6
    cuts.AddCut("12410013","1111102007032230000","01634031000000d0"); // rapidity variation  y<0.5
    cuts.AddCut("12410013","1111102007032230000","01631061000000d0"); // alpha meson variation 1   0<alpha<0.8
    cuts.AddCut("12410013","1111102007032230000","01631051000000d0"); // alpha meson variation 2  0<alpha<0.75
  } else if (trainConfig == 324){ // opening angle variations
    cuts.AddCut("12410013","1111102007032230000","0163103100000040"); // min opening angle 0.0152
    cuts.AddCut("12410013","1111102007032230000","0163103100000050"); // min opening angle 0.0202
    cuts.AddCut("12410013","1111102007032230000","0163103100000070"); // min opening angle 0.016
    cuts.AddCut("12410013","1111102007032230000","0163103100000080"); // min opening angle 0.018
    cuts.AddCut("12410013","1111102007032230000","0163103100000090"); // min opening angle 0.018
  } else if (trainConfig == 325){ // TM variations
    cuts.AddCut("12410013","1111102003032230000","01631031000000d0"); // fixed window
    cuts.AddCut("12410013","1111102006032230000","01631031000000d0"); // tm pt dependent var 1
    cuts.AddCut("12410013","1111102008032230000","01631031000000d0"); // tm pt dependent var 2
    cuts.AddCut("12410013","1111102009032230000","01631031000000d0"); // tm pt dependent var 3

  //-----------------------------------------------------------------------------------------------
  // Systematics variations 40-80 EMC
  //-----------------------------------------------------------------------------------------------
  } else if (trainConfig == 330){ // EMCal clusters non lin variations 40-80 %
    cuts.AddCut("14810013","1111171007032230000","01631031000000d0");
    cuts.AddCut("14810013","1111172007032230000","01631031000000d0");
    cuts.AddCut("14810013","1111181007032230000","01631031000000d0");
    cuts.AddCut("14810013","1111182007032230000","01631031000000d0");
  } else if (trainConfig == 331){ // second set of variations CLUSTER
    cuts.AddCut("14810013","1111102007022230000","01631031000000d0"); // min energy cluster variation 1  600 MeV
    cuts.AddCut("14810013","1111102007042230000","01631031000000d0"); // min energy cluster variation 2  800 MeV
    cuts.AddCut("14810013","1111102007052230000","01631031000000d0"); // min energy cluster variation 3  900 MeV
    cuts.AddCut("14810013","1111102007032230000","01631031000000d0"); // min/max M02  0.1<M<0.5
    cuts.AddCut("14810013","1111102007032200000","01631031000000d0"); // min/max M02  0.1<M<100
  } else if (trainConfig == 332){ // third set of variations CLUSTER
    cuts.AddCut("14810013","1111102007032250000","01631031000000d0"); // min/max M02  0.1<M<0.3
    cuts.AddCut("14810013","1111102007032260000","01631031000000d0"); // min/max M02  0.1<M<0.27
    cuts.AddCut("14810013","1111102007031230000","01631031000000d0"); // min number of cells variation 1  1 cell
    cuts.AddCut("14810013","1112102007032230000","01631031000000d0"); // only modules with TRD infront
    cuts.AddCut("14810013","1111302007032230000","01631031000000d0"); // no modules with TRD infront
  } else if (trainConfig == 333){ // third set of variations MESON
    cuts.AddCut("14810013","1111102007032230000","01633031000000d0"); // rapidity variation  y<0.6
    cuts.AddCut("14810013","1111102007032230000","01634031000000d0"); // rapidity variation  y<0.5
    cuts.AddCut("14810013","1111102007032230000","01631061000000d0"); // alpha meson variation 1   0<alpha<0.8
    cuts.AddCut("14810013","1111102007032230000","01631051000000d0"); // alpha meson variation 2  0<alpha<0.75
  } else if (trainConfig == 334){ // opening angle variations
    cuts.AddCut("14810013","1111102007032230000","0163103100000040"); // min opening angle 0.0152
    cuts.AddCut("14810013","1111102007032230000","0163103100000050"); // min opening angle 0.0202
    cuts.AddCut("14810013","1111102007032230000","0163103100000070"); // min opening angle 0.016
    cuts.AddCut("14810013","1111102007032230000","0163103100000080"); // min opening angle 0.018
    cuts.AddCut("14810013","1111102007032230000","0163103100000090"); // min opening angle 0.018
  } else if (trainConfig == 335){ // TM variations
    cuts.AddCut("14810013","1111102003032230000","01631031000000d0"); // fixed window
    cuts.AddCut("14810013","1111102006032230000","01631031000000d0"); // tm pt dependent var 1
    cuts.AddCut("14810013","1111102008032230000","01631031000000d0"); // tm pt dependent var 2
    cuts.AddCut("14810013","1111102009032230000","01631031000000d0"); // tm pt dependent var 3

  //-----------------------------------------------------------------------------------------------
  // Systematics variations 0-40 EMC
  //-----------------------------------------------------------------------------------------------
  } else if (trainConfig == 340){ // EMCal clusters non lin variations 0-40 %
    cuts.AddCut("10410013","1111171007032230000","01631031000000d0");
    cuts.AddCut("10410013","1111172007032230000","01631031000000d0");
    cuts.AddCut("10410013","1111181007032230000","01631031000000d0");
    cuts.AddCut("10410013","1111182007032230000","01631031000000d0");
  } else if (trainConfig == 341){ // second set of variations CLUSTER
    cuts.AddCut("10410013","1111102007022230000","01631031000000d0"); // min energy cluster variation 1  600 MeV
    cuts.AddCut("10410013","1111102007042230000","01631031000000d0"); // min energy cluster variation 2  800 MeV
    cuts.AddCut("10410013","1111102007052230000","01631031000000d0"); // min energy cluster variation 3  900 MeV
    cuts.AddCut("10410013","1111102007032230000","01631031000000d0"); // min/max M02  0.1<M<0.5
    cuts.AddCut("10410013","1111102007032200000","01631031000000d0"); // min/max M02  0.1<M<100
  } else if (trainConfig == 342){ // third set of variations CLUSTER
    cuts.AddCut("10410013","1111102007032250000","01631031000000d0"); // min/max M02  0.1<M<0.3
    cuts.AddCut("10410013","1111102007032260000","01631031000000d0"); // min/max M02  0.1<M<0.27
    cuts.AddCut("10410013","1111102007031230000","01631031000000d0"); // min number of cells variation 1  1 cell
    cuts.AddCut("10410013","1112102007032230000","01631031000000d0"); // only modules with TRD infront
    cuts.AddCut("10410013","1111302007032230000","01631031000000d0"); // no modules with TRD infront
  } else if (trainConfig == 343){ // third set of variations MESON
    cuts.AddCut("10410013","1111102007032230000","01633031000000d0"); // rapidity variation  y<0.6
    cuts.AddCut("10410013","1111102007032230000","01634031000000d0"); // rapidity variation  y<0.5
    cuts.AddCut("10410013","1111102007032230000","01631061000000d0"); // alpha meson variation 1   0<alpha<0.8
    cuts.AddCut("10410013","1111102007032230000","01631051000000d0"); // alpha meson variation 2  0<alpha<0.75
  } else if (trainConfig == 344){ // opening angle variations
    cuts.AddCut("10410013","1111102007032230000","0163103100000040"); // min opening angle 0.0152
    cuts.AddCut("10410013","1111102007032230000","0163103100000050"); // min opening angle 0.0202
    cuts.AddCut("10410013","1111102007032230000","0163103100000070"); // min opening angle 0.016
    cuts.AddCut("10410013","1111102007032230000","0163103100000080"); // min opening angle 0.018
    cuts.AddCut("10410013","1111102007032230000","0163103100000090"); // min opening angle 0.018
  } else if (trainConfig == 345){ // TM variations
    cuts.AddCut("10410013","1111102003032230000","01631031000000d0"); // fixed window
    cuts.AddCut("10410013","1111102006032230000","01631031000000d0"); // tm pt dependent var 1
    cuts.AddCut("10410013","1111102008032230000","01631031000000d0"); // tm pt dependent var 2
    cuts.AddCut("10410013","1111102009032230000","01631031000000d0"); // tm pt dependent var 3

  //-----------------------------------------------------------------------------------------------
  // Systematics variations 0-80 EMC
  //-----------------------------------------------------------------------------------------------
  } else if (trainConfig == 350){ // EMCal clusters non lin variations 0-80 %
    cuts.AddCut("10810013","1111171007032230000","01631031000000d0");
    cuts.AddCut("10810013","1111172007032230000","01631031000000d0");
    cuts.AddCut("10810013","1111181007032230000","01631031000000d0");
    cuts.AddCut("10810013","1111182007032230000","01631031000000d0");
  } else if (trainConfig == 351){ // second set of variations CLUSTER
    cuts.AddCut("10810013","1111102007022230000","01631031000000d0"); // min energy cluster variation 1  600 MeV
    cuts.AddCut("10810013","1111102007042230000","01631031000000d0"); // min energy cluster variation 2  800 MeV
    cuts.AddCut("10810013","1111102007052230000","01631031000000d0"); // min energy cluster variation 3  900 MeV
    cuts.AddCut("10810013","1111102007032230000","01631031000000d0"); // min/max M02  0.1<M<0.5
    cuts.AddCut("10810013","1111102007032200000","01631031000000d0"); // min/max M02  0.1<M<100
  } else if (trainConfig == 352){ // third set of variations CLUSTER
    cuts.AddCut("10810013","1111102007032250000","01631031000000d0"); // min/max M02  0.1<M<0.3
    cuts.AddCut("10810013","1111102007032260000","01631031000000d0"); // min/max M02  0.1<M<0.27
    cuts.AddCut("10810013","1111102007031230000","01631031000000d0"); // min number of cells variation 1  1 cell
    cuts.AddCut("10810013","1112102007032230000","01631031000000d0"); // only modules with TRD infront
    cuts.AddCut("10810013","1111302007032230000","01631031000000d0"); // no modules with TRD infront
  } else if (trainConfig == 353){ // third set of variations MESON
    cuts.AddCut("10810013","1111102007032230000","01633031000000d0"); // rapidity variation  y<0.6
    cuts.AddCut("10810013","1111102007032230000","01634031000000d0"); // rapidity variation  y<0.5
    cuts.AddCut("10810013","1111102007032230000","01631061000000d0"); // alpha meson variation 1   0<alpha<0.8
    cuts.AddCut("10810013","1111102007032230000","01631051000000d0"); // alpha meson variation 2  0<alpha<0.75
  } else if (trainConfig == 354){ // opening angle variations
    cuts.AddCut("10810013","1111102007032230000","0163103100000040"); // min opening angle 0.0152
    cuts.AddCut("10810013","1111102007032230000","0163103100000050"); // min opening angle 0.0202
    cuts.AddCut("10810013","1111102007032230000","0163103100000070"); // min opening angle 0.016
    cuts.AddCut("10810013","1111102007032230000","0163103100000080"); // min opening angle 0.018
    cuts.AddCut("10810013","1111102007032230000","0163103100000090"); // min opening angle 0.018
  } else if (trainConfig == 355){ // TM variations
    cuts.AddCut("10810013","1111102003032230000","01631031000000d0"); // fixed window
    cuts.AddCut("10810013","1111102006032230000","01631031000000d0"); // tm pt dependent var 1
    cuts.AddCut("10810013","1111102008032230000","01631031000000d0"); // tm pt dependent var 2
    cuts.AddCut("10810013","1111102009032230000","01631031000000d0"); // tm pt dependent var 3

  // **********************************************************************************************************
  // ***************************** PHOS configurations XeXe run 2 *********************************************
  // **********************************************************************************************************
  } else if (trainConfig == 400){ // PHOS clusters - 0-80% centrality for XeXe PHOS cluster QA
    cuts.AddCut("10810013","2446600040013300000","0163103100000010"); // 0-80
    cuts.AddCut("10810013","2446601040013300000","0163103100000010"); // 0-80
  } else if (trainConfig == 401){ // PHOS clusters - centrality for XeXe PHOS cluster QA
    cuts.AddCut("10210013","2446600040013300000","0163103100000010"); // 0-20
    cuts.AddCut("12410013","2446600040013300000","0163103100000010"); // 20-40
    cuts.AddCut("10410013","2446600040013300000","0163103100000010"); // 0-40
    cuts.AddCut("14810013","2446600040013300000","0163103100000010"); // 40-80
  } else if (trainConfig == 402){ // PHOS clusters - centrality for XeXe PHOS cluster QA
    cuts.AddCut("10210013","2446601040013300000","0163103100000010"); // 0-20
    cuts.AddCut("12410013","2446601040013300000","0163103100000010"); // 20-40
    cuts.AddCut("10410013","2446601040013300000","0163103100000010"); // 0-40
    cuts.AddCut("14810013","2446601040013300000","0163103100000010"); // 40-80

  //-----------------------------------------------------------------------------------------------
  // Systematics variations 0-20 PHOS
  //-----------------------------------------------------------------------------------------------
  } else if (trainConfig == 410){ // PHOS clusters - non lin variations 0-20%
    cuts.AddCut("10210013","2446601040013300000","0163103100000010");
    cuts.AddCut("10210013","2446671040013300000","0163103100000010");
    cuts.AddCut("10210013","2446672040013300000","0163103100000010");
    cuts.AddCut("10210013","2446681040013300000","0163103100000010");
    cuts.AddCut("10210013","2446682040013300000","0163103100000010");

  //-----------------------------------------------------------------------------------------------
  // Systematics variations 20-40 PHOS
  //-----------------------------------------------------------------------------------------------
  } else if (trainConfig == 420){ // PHOS clusters - non lin variations 20-40%
    cuts.AddCut("12410013","2446601040013300000","0163103100000010");
    cuts.AddCut("12410013","2446671040013300000","0163103100000010");
    cuts.AddCut("12410013","2446672040013300000","0163103100000010");
    cuts.AddCut("12410013","2446681040013300000","0163103100000010");
    cuts.AddCut("12410013","2446682040013300000","0163103100000010");

  //-----------------------------------------------------------------------------------------------
  // Systematics variations 40-80 PHOS
  //-----------------------------------------------------------------------------------------------
  } else if (trainConfig == 430){ // PHOS clusters - non lin variations 40-80%
    cuts.AddCut("14810013","2446601040013300000","0163103100000010");
    cuts.AddCut("14810013","2446671040013300000","0163103100000010");
    cuts.AddCut("14810013","2446672040013300000","0163103100000010");
    cuts.AddCut("14810013","2446681040013300000","0163103100000010");
    cuts.AddCut("14810013","2446682040013300000","0163103100000010");

  //-----------------------------------------------------------------------------------------------
  // Systematics variations 0-40 PHOS
  //-----------------------------------------------------------------------------------------------
  } else if (trainConfig == 440){ // PHOS clusters - non lin variations 0-40%
    cuts.AddCut("10410013","2446601040013300000","0163103100000010");
    cuts.AddCut("10410013","2446671040013300000","0163103100000010");
    cuts.AddCut("10410013","2446672040013300000","0163103100000010");
    cuts.AddCut("10410013","2446681040013300000","0163103100000010");
    cuts.AddCut("10410013","2446682040013300000","0163103100000010");

  //-----------------------------------------------------------------------------------------------
  // Systematics variations 0-80 PHOS
  //-----------------------------------------------------------------------------------------------
  } else if (trainConfig == 450){ // PHOS clusters - non lin variations 0-80%
    cuts.AddCut("10810013","2446601040013300000","0163103100000010");
    cuts.AddCut("10810013","2446671040013300000","0163103100000010");
    cuts.AddCut("10810013","2446672040013300000","0163103100000010");
    cuts.AddCut("10810013","2446681040013300000","0163103100000010");
    cuts.AddCut("10810013","2446682040013300000","0163103100000010");

  // **********************************************************************************************************
  // ***************************** EMC configurations PbPb run 2 EMC trigger **********************************************
  // **********************************************************************************************************
  } else if (trainConfig == 701){ // EMCAL clusters central
    cuts.AddCut("10183013","1111100057032230000","01631031000000d0"); // 0-10
    cuts.AddCut("30183013","1111100057032230000","01631031000000d0"); // 0-5
    cuts.AddCut("31283013","1111100057032230000","01631031000000d0"); // 5-10
  } else if (trainConfig == 702){ // EMCAL clusters semi-central
    cuts.AddCut("11283013","1111100057032230000","01631031000000d0"); // 10-20
    cuts.AddCut("12383013","1111100057032230000","01631031000000d0"); // 20-30
    cuts.AddCut("13483013","1111100057032230000","01631031000000d0"); // 30-40
    cuts.AddCut("12483013","1111100057032230000","01631031000000d0"); // 20-40
  } else if (trainConfig == 703){ // EMCAL clusters semi-central
    cuts.AddCut("14583013","1111100057032230000","01631031000000d0"); // 40-50
    cuts.AddCut("14683013","1111100057032230000","01631031000000d0"); // 40-60
    cuts.AddCut("15683013","1111100057032230000","01631031000000d0"); // 50-60
  } else if (trainConfig == 704){ // EMCAL clusters peripheral
    cuts.AddCut("16783013","1111100057032230000","01631031000000d0"); // 60-70
    cuts.AddCut("17883013","1111100057032230000","01631031000000d0"); // 70-80
    cuts.AddCut("18983013","1111100057032230000","01631031000000d0"); // 80-90
    cuts.AddCut("16883013","1111100057032230000","01631031000000d0"); // 60-80


  } else {
    Error(Form("GammaConvCalo_%i",trainConfig), "wrong trainConfig variable no cuts have been specified for the configuration");
    return;
  }

  if(!cuts.AreValid()){
    cout << "\n\n****************************************************" << endl;
    cout << "ERROR: No valid cuts stored in CutHandlerCalo! Returning..." << endl;
    cout << "****************************************************\n\n" << endl;
    return;
  }

  Int_t numberOfCuts = cuts.GetNCuts();

  TList *HeaderList = new TList();
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

  TList *EventCutList   = new TList();
  TList *ClusterCutList = new TList();
  TList *MesonCutList   = new TList();


  EventCutList->SetOwner(kTRUE);
  AliConvEventCuts **analysisEventCuts        = new AliConvEventCuts*[numberOfCuts];
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

    analysisEventCuts[i]    = new AliConvEventCuts();
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

    if (trainConfig == 34 || trainConfig == 35 ){
      if (periodName.CompareTo("LHC14a1a") ==0 || periodName.CompareTo("LHC14a1b") ==0 || periodName.CompareTo("LHC14a1c") ==0 ){
        if ( i == 0 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_0010TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_0010TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0010V0M","Eta_Fit_Data_PbPb_2760GeV_0010V0M");
        if ( i == 1 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_0010TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_0010TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0010V0M","Eta_Fit_Data_PbPb_2760GeV_0010V0M");
      }
    }
    if (trainConfig == 36 || trainConfig == 37){
      if (periodName.CompareTo("LHC14a1a") ==0 || periodName.CompareTo("LHC14a1b") ==0 || periodName.CompareTo("LHC14a1c") ==0 ){
        if ( i == 0 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_2050TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_2050TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_2050V0M","Eta_Fit_Data_PbPb_2760GeV_2050V0M");
        if ( i == 1 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_2050TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_2050TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_2050V0M","Eta_Fit_Data_PbPb_2760GeV_2050V0M");
      }
    }

    if (trainConfig == 38 || trainConfig == 39 || trainConfig == 43 || trainConfig == 44 ){
      if (periodName.CompareTo("LHC14a1a") ==0 || periodName.CompareTo("LHC14a1b") ==0 || periodName.CompareTo("LHC14a1c") ==0 ){
        if ( i == 0 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_0010TPC",periodName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_0010TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0010V0M","Eta_Fit_Data_PbPb_2760GeV_0010V0M");
        if ( i == 1 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_0010TPC",periodName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_0010TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0010V0M","Eta_Fit_Data_PbPb_2760GeV_0010V0M");
        if ( i == 2 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_0010TPC",periodName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_0010TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0010V0M","Eta_Fit_Data_PbPb_2760GeV_0010V0M");
      }
    }

    if (trainConfig == 40 || trainConfig == 41){
      if (periodName.CompareTo("LHC14a1a") ==0 || periodName.CompareTo("LHC14a1b") ==0 || periodName.CompareTo("LHC14a1c") ==0 ){
        if ( i == 0 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_2050TPC",periodName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_2050TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_2050V0M","Eta_Fit_Data_PbPb_2760GeV_2050V0M");
        if ( i == 1 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_2050TPC",periodName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_2050TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_2050V0M","Eta_Fit_Data_PbPb_2760GeV_2050V0M");
        if ( i == 2 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_2050TPC",periodName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_2050TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_2050V0M","Eta_Fit_Data_PbPb_2760GeV_2050V0M");
      }
    }

    analysisEventCuts[i]->SetDebugLevel(0);
    analysisEventCuts[i]->SetLightOutput(runLightOutput);
    analysisEventCuts[i]->InitializeCutsFromCutString((cuts.GetEventCut(i)).Data());
    EventCutList->Add(analysisEventCuts[i]);
    analysisEventCuts[i]->SetCorrectionTaskSetting(corrTaskSetting);
    analysisEventCuts[i]->SetFillCutHistograms("",kFALSE);
    if (periodName.CompareTo("LHC14a1b") ==0 || periodName.CompareTo("LHC14a1c") ==0 ){
      if (headerSelectionInt == 1 || headerSelectionInt == 4 || headerSelectionInt == 12 ) analysisEventCuts[i]->SetAddedSignalPDGCode(111);
      if (headerSelectionInt == 2 || headerSelectionInt == 5 || headerSelectionInt == 13 ) analysisEventCuts[i]->SetAddedSignalPDGCode(221);
    }

    analysisClusterCuts[i]  = new AliCaloPhotonCuts(isMC);
    analysisClusterCuts[i]->SetHistoToModifyAcceptance(histoAcc);
    analysisClusterCuts[i]->SetV0ReaderName(V0ReaderName);
    analysisClusterCuts[i]->SetCorrectionTaskSetting(corrTaskSetting);
    analysisClusterCuts[i]->SetCaloTrackMatcherName(TrackMatcherName);
    analysisClusterCuts[i]->SetLightOutput(runLightOutput);
    analysisClusterCuts[i]->InitializeCutsFromCutString((cuts.GetClusterCut(i)).Data());
    ClusterCutList->Add(analysisClusterCuts[i]);
    analysisClusterCuts[i]->SetExtendedMatchAndQA(enableExtMatchAndQA);
    analysisClusterCuts[i]->SetFillCutHistograms("");

    analysisMesonCuts[i]    = new AliConversionMesonCuts();
    analysisMesonCuts[i]->SetLightOutput(runLightOutput);
    analysisMesonCuts[i]->InitializeCutsFromCutString((cuts.GetMesonCut(i)).Data());
    analysisMesonCuts[i]->SetIsMergedClusterCut(2);
    analysisMesonCuts[i]->SetCaloMesonCutsObject(analysisClusterCuts[i]);
    MesonCutList->Add(analysisMesonCuts[i]);
    analysisMesonCuts[i]->SetFillCutHistograms("");
    analysisEventCuts[i]->SetAcceptedHeader(HeaderList);
  }
  task->SetAllowOverlapHeaders(enableHeaderOverlap);
  task->SetEventCutList(numberOfCuts,EventCutList);
  task->SetCaloCutList(numberOfCuts,ClusterCutList);
  task->SetMesonCutList(numberOfCuts,MesonCutList);
  task->SetDoMesonAnalysis(kTRUE);
  task->SetCorrectionTaskSetting(corrTaskSetting);
  task->SetDoMesonQA(enableQAMesonTask);        //Attention new switch for Pi0 QA
  task->SetDoClusterQA(enableQAClusterTask);    //Attention new switch small for Cluster QA
  task->SetDoTHnSparse(isUsingTHnSparse);
  task->SetProduceTreeEOverP(doTreeEOverP);
  task->SetEnableSortingOfMCClusLabels(enableSortingMCLabels);
  if(enableExtMatchAndQA > 1){ task->SetPlotHistsExtQA(kTRUE);}

  //connect containers
  AliAnalysisDataContainer *coutput =
    mgr->CreateContainer(!(corrTaskSetting.CompareTo("")) ? Form("GammaCalo_%i",trainConfig) : Form("GammaCalo_%i_%s",trainConfig,corrTaskSetting.Data()), TList::Class(),
              AliAnalysisManager::kOutputContainer,Form("GammaCalo_%i.root",trainConfig) );

  mgr->AddTask(task);
  mgr->ConnectInput(task,0,cinput);
  mgr->ConnectOutput(task,1,coutput);
  return;
}
