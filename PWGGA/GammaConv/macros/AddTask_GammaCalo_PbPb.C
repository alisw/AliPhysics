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
                              TString   additionalTrainConfig           = "0"                   // additional counter for trainconfig
                ) {

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
      fV0ReaderV1->SetDeltaAODBranchName(Form("GammaConv_%s_gamma",cutnumberAODBranch.Data()));
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
    cuts.AddCut("20100013","11111020530a2230000","01631031000000d0"); // 0-5 V0M   // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("21100013","11111020530a2230000","01631031000000d0"); // 5-10 V0M  // reproduce Astrids cuts with opening angle cut
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
    cuts.AddCut("20100013","11111000530a2230000","01631031000000d0"); // 0-5 V0M   // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("21100013","11111000530a2230000","01631031000000d0"); // 5-10 V0M  // reproduce Astrids cuts with opening angle cut
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
    cuts.AddCut("20100023","11111020530a2230000","01631031000000d0"); // 0-5 V0M   // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("21100023","11111020530a2230000","01631031000000d0"); // 5-10 V0M  // reproduce Astrids cuts with opening angle cut
  } else if (trainConfig == 62){
    cuts.AddCut("10100023","11111020530a2230000","01631031000000d0"); // 0-10 V0M  // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("11200023","11111020530a2230000","01631031000000d0"); // 10-20 V0M // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("20100023","11111020530a2230000","01631031000000d0"); // 0-5 V0M   // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("21100023","11111020530a2230000","01631031000000d0"); // 5-10 V0M  // reproduce Astrids cuts with opening angle cut
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


  // run 2 configs
  // first look
  } else if (trainConfig == 201){ // EMCAL clusters
    cuts.AddCut("20110113","1111100003032230000","0163103100000050"); // 0-10
    cuts.AddCut("21210113","1111100003032230000","0163103100000050"); // 10-20
    cuts.AddCut("22510113","1111100003032230000","0163103100000050"); // 20-50
    cuts.AddCut("25910113","1111100003032230000","0163103100000050"); // 50-90
    cuts.AddCut("20010113","1111100003032230000","0163103100000050"); // 0-100
  } else if (trainConfig == 202){ // EMCAL clusters
    cuts.AddCut("10110113","1111100053032230000","0163103100000050"); // 0-10
    cuts.AddCut("11210113","1111100053032230000","0163103100000050"); // 10-20
    cuts.AddCut("12510113","1111100053032230000","0163103100000050"); // 20-50
    cuts.AddCut("15910113","1111100053032230000","0163103100000050"); // 50-90
    cuts.AddCut("10010113","1111100053032230000","0163103100000050"); // 0-100

  } else if (trainConfig == 203){ // EMCAL clusters - change opening angle
    cuts.AddCut("10110113","1111100003032230000","0163103100000040"); // 0-10
    cuts.AddCut("11210113","1111100003032230000","0163103100000040"); // 10-20
    cuts.AddCut("12510113","1111100003032230000","0163103100000040"); // 20-50
    cuts.AddCut("15910113","1111100003032230000","0163103100000040"); // 50-90
    cuts.AddCut("10010113","1111100003032230000","0163103100000040"); // 0-100
  } else if (trainConfig == 204){ // EMCAL clusters - timing cut +- 50ns
    cuts.AddCut("10110113","1111100053032230000","0163103100000050"); // 0-10
    cuts.AddCut("11210113","1111100053032230000","0163103100000050"); // 10-20
    cuts.AddCut("12510113","1111100053032230000","0163103100000050"); // 20-50
    cuts.AddCut("15910113","1111100053032230000","0163103100000050"); // 50-90
    cuts.AddCut("10010113","1111100053032230000","0163103100000050"); // 0-100

  } else if (trainConfig == 205){ // EMCAL clusters - user defined header!
    cuts.AddCut("10110123","1111100003032230000","0163103100000050"); // 0-10
    cuts.AddCut("11210123","1111100003032230000","0163103100000050"); // 10-20
    cuts.AddCut("12510123","1111100003032230000","0163103100000050"); // 20-50
    cuts.AddCut("15910123","1111100003032230000","0163103100000050"); // 50-90
    cuts.AddCut("10010123","1111100003032230000","0163103100000050"); // 0-100

  } else if (trainConfig == 206){ // EMCAL clusters - 205 duplicant for header setting
    cuts.AddCut("10110123","1111100003032230000","0163103100000050"); // 0-10
    cuts.AddCut("11210123","1111100003032230000","0163103100000050"); // 10-20
    cuts.AddCut("12510123","1111100003032230000","0163103100000050"); // 20-50
    cuts.AddCut("15910123","1111100003032230000","0163103100000050"); // 50-90
    cuts.AddCut("10010123","1111100003032230000","0163103100000050"); // 0-100

  } else if (trainConfig == 207){ // EMCAL clusters - timing cut +- 50ns
    cuts.AddCut("10110113","1111100053032230000","0163103100000050"); // 0-10
    cuts.AddCut("11210113","1111100053032230000","0163103100000050"); // 10-20
    cuts.AddCut("12510113","1111100053032230000","0163103100000050"); // 20-50
    cuts.AddCut("15910113","1111100053032230000","0163103100000050"); // 50-90
  } else if (trainConfig == 208){ // EMCAL clusters - correction convcalo f1
    cuts.AddCut("10110113","1111181053032230000","0163103100000050"); // 0-10
    cuts.AddCut("11210113","1111181053032230000","0163103100000050"); // 10-20
    cuts.AddCut("12510113","1111181053032230000","0163103100000050"); // 20-50
    cuts.AddCut("15910113","1111181053032230000","0163103100000050"); // 50-90
  } else if (trainConfig == 209){ // EMCAL clusters - correction calocalo f2
    cuts.AddCut("10110113","1111192053032230000","0163103100000050"); // 0-10
    cuts.AddCut("11210113","1111192053032230000","0163103100000050"); // 10-20
    cuts.AddCut("12510113","1111192053032230000","0163103100000050"); // 20-50
    cuts.AddCut("15910113","1111192053032230000","0163103100000050"); // 50-90
  } else if (trainConfig == 210){ // EMCAL clusters - 0-90% centrality for PbPb EMCal cluster QA
    cuts.AddCut("10910113","1111100003032230000","0163103100000050"); // 0-90
  } else if (trainConfig == 211){ // EMCAL clusters - 0-90% centrality for PbPb EMCal cluster QA
    cuts.AddCut("10910113","1111181003032230000","0163103100000050"); // 0-90 convcalo correction f1
    cuts.AddCut("10910113","1111182003032230000","0163103100000050"); // 0-90 calocalo correction f1
    cuts.AddCut("10910113","1111191003032230000","0163103100000050"); // 0-90 convcalo correction f2
    cuts.AddCut("10910113","1111192003032230000","0163103100000050"); // 0-90 calocalo correction f2

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
    cuts.AddCut("20110113","11111020530a2230000","01631031000000d0"); // 0-5 V0M   // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("21110113","11111020530a2230000","01631031000000d0"); // 5-10 V0M  // reproduce Astrids cuts with opening angle cut
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
    cuts.AddCut("20110113","11111000530a2230000","01631031000000d0"); // 0-5 V0M   // reproduce Astrids cuts with opening angle cut
    cuts.AddCut("21110113","11111000530a2230000","01631031000000d0"); // 5-10 V0M  // reproduce Astrids cuts with opening angle cut
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

  // Xe-xe configurations
  } else if (trainConfig == 300){ // EMCAL clusters - 0-80% centrality for XeXe EMCal cluster QA
    cuts.AddCut("10810113","1111100007032230000","01631031000000d0"); // 0-80
  } else if (trainConfig == 301){ // EMCAL clusters - 0-80% centrality for XeXe EMCal cluster QA
    cuts.AddCut("10210113","1111100007032230000","01631031000000d0"); // 0-20
    cuts.AddCut("12410113","1111100007032230000","01631031000000d0"); // 20-40
    cuts.AddCut("10410113","1111100007032230000","01631031000000d0"); // 0-40
    cuts.AddCut("14810113","1111100007032230000","01631031000000d0"); // 40-80
  } else if (trainConfig == 302){ // EMCAL clusters - 0-80% centrality for XeXe EMCal cluster QA - TB nl
    cuts.AddCut("10810113","1111102007032230000","01631031000000d0"); // 0-80
  } else if (trainConfig == 303){ // EMCAL clusters - diff centralities for XeXe EMCal cluster QA - TB nl
    cuts.AddCut("10210113","1111102007032230000","01631031000000d0"); // 0-20
    cuts.AddCut("12410113","1111102007032230000","01631031000000d0"); // 20-40
    cuts.AddCut("10410113","1111102007032230000","01631031000000d0"); // 0-40
    cuts.AddCut("14810113","1111102007032230000","01631031000000d0"); // 40-80


  // Xe-Xe configurations PHOS
  } else if (trainConfig == 400){ // PHOS clusters - 0-80% centrality for XeXe PHOS cluster QA
    cuts.AddCut("10810113","2446600040013300000","0163103100000010"); // 0-80
    cuts.AddCut("10810113","2446601040013300000","0163103100000010"); // 0-80
  } else if (trainConfig == 401){ // PHOS clusters - centrality for XeXe PHOS cluster QA
    cuts.AddCut("10210113","2446600040013300000","0163103100000010"); // 0-20
    cuts.AddCut("12410113","2446600040013300000","0163103100000010"); // 20-40
    cuts.AddCut("10410113","2446600040013300000","0163103100000010"); // 0-40
    cuts.AddCut("14810113","2446600040013300000","0163103100000010"); // 40-80
  } else if (trainConfig == 402){ // PHOS clusters - centrality for XeXe PHOS cluster QA
    cuts.AddCut("10210113","2446601040013300000","0163103100000010"); // 0-20
    cuts.AddCut("12410113","2446601040013300000","0163103100000010"); // 20-40
    cuts.AddCut("10410113","2446601040013300000","0163103100000010"); // 0-40
    cuts.AddCut("14810113","2446601040013300000","0163103100000010"); // 40-80

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
    TObjString *Header1 = new TObjString("BOX");
    HeaderList->Add(Header1);
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
    TString TrackMatcherName = Form("CaloTrackMatcher_%s",caloCutPos.Data());
    if( !(AliCaloTrackMatcher*)mgr->GetTask(TrackMatcherName.Data()) ){
      AliCaloTrackMatcher* fTrackMatcher = new AliCaloTrackMatcher(TrackMatcherName.Data(),caloCutPos.Atoi());
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

//     analysisEventCuts[i]->SetDebugLevel(3);
    analysisEventCuts[i]->SetLightOutput(runLightOutput);
    analysisEventCuts[i]->InitializeCutsFromCutString((cuts.GetEventCut(i)).Data());
    EventCutList->Add(analysisEventCuts[i]);
    analysisEventCuts[i]->SetCorrectionTaskSetting(corrTaskSetting);
    analysisEventCuts[i]->SetFillCutHistograms("",kFALSE);

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
              AliAnalysisManager::kOutputContainer,!(corrTaskSetting.CompareTo("")) ? Form("GammaCalo_%i.root",trainConfig) : Form("GammaCalo_%i_%s.root",trainConfig,corrTaskSetting.Data()));

  mgr->AddTask(task);
  mgr->ConnectInput(task,0,cinput);
  mgr->ConnectOutput(task,1,coutput);
  return;
}
