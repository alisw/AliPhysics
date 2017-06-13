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
      }
    }
  }
  TString sAdditionalTrainConfig = rAdditionalTrainConfig->GetString();
  if (sAdditionalTrainConfig.Atoi() > 0){
    trainConfig = trainConfig + sAdditionalTrainConfig.Atoi();
    cout << "INFO: AddTask_GammaCalo_PbPb running additionalTrainConfig '" << sAdditionalTrainConfig.Atoi() << "', train config: '" << trainConfig << "'" << endl;
  }

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
    cuts.AddCut("50180013","1111100053032230000","0163103100000050"); // 0-10%
    cuts.AddCut("51280013","1111100053032230000","0163103100000050"); // 10-20%
    cuts.AddCut("52580013","1111100053032230000","0163103100000050"); // 20-50%
    cuts.AddCut("55880013","1111100053032230000","0163103100000050"); // 50-80%
  } else if (trainConfig == 35){ // EMCAL clusters no timing cut
    cuts.AddCut("50180013","1111100003032230000","0163103100000050"); // 0-10%
    cuts.AddCut("51280013","1111100003032230000","0163103100000050"); // 10-20%
    cuts.AddCut("52580013","1111100003032230000","0163103100000050"); // 20-50%
    cuts.AddCut("55880013","1111100003032230000","0163103100000050"); // 50-80%
    
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
    cuts.AddCut("50110113","1111100053032230000","0163103100000050"); // 0-10
    cuts.AddCut("51210113","1111100053032230000","0163103100000050"); // 10-20
    cuts.AddCut("52510113","1111100053032230000","0163103100000050"); // 20-50
    cuts.AddCut("55910113","1111100053032230000","0163103100000050"); // 50-90
    cuts.AddCut("50010113","1111100053032230000","0163103100000050"); // 0-100
    
  } else if (trainConfig == 203){ // EMCAL clusters - change opening angle
    cuts.AddCut("50110113","1111100003032230000","0163103100000040"); // 0-10
    cuts.AddCut("51210113","1111100003032230000","0163103100000040"); // 10-20
    cuts.AddCut("52510113","1111100003032230000","0163103100000040"); // 20-50
    cuts.AddCut("55910113","1111100003032230000","0163103100000040"); // 50-90
    cuts.AddCut("50010113","1111100003032230000","0163103100000040"); // 0-100
  } else if (trainConfig == 204){ // EMCAL clusters - timing cut +- 50ns
    cuts.AddCut("50110113","1111100053032230000","0163103100000050"); // 0-10
    cuts.AddCut("51210113","1111100053032230000","0163103100000050"); // 10-20
    cuts.AddCut("52510113","1111100053032230000","0163103100000050"); // 20-50
    cuts.AddCut("55910113","1111100053032230000","0163103100000050"); // 50-90
    cuts.AddCut("50010113","1111100053032230000","0163103100000050"); // 0-100
    
  } else if (trainConfig == 205){ // EMCAL clusters - user defined header!
    cuts.AddCut("50110123","1111100003032230000","0163103100000050"); // 0-10
    cuts.AddCut("51210123","1111100003032230000","0163103100000050"); // 10-20
    cuts.AddCut("52510123","1111100003032230000","0163103100000050"); // 20-50
    cuts.AddCut("55910123","1111100003032230000","0163103100000050"); // 50-90
    cuts.AddCut("50010123","1111100003032230000","0163103100000050"); // 0-100
    
  } else if (trainConfig == 206){ // EMCAL clusters - 205 duplicant for header setting
    cuts.AddCut("50110123","1111100003032230000","0163103100000050"); // 0-10
    cuts.AddCut("51210123","1111100003032230000","0163103100000050"); // 10-20
    cuts.AddCut("52510123","1111100003032230000","0163103100000050"); // 20-50
    cuts.AddCut("55910123","1111100003032230000","0163103100000050"); // 50-90
    cuts.AddCut("50010123","1111100003032230000","0163103100000050"); // 0-100
    
  } else if (trainConfig == 210){ // EMCAL clusters - 0-90% centrality for PbPb EMCal cluster QA
    cuts.AddCut("50910113","1111100003032230000","0163103100000050"); // 0-90
    

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
    } else {
      TObjString *Header1 = new TObjString("pi0_1");
      HeaderList->Add(Header1);
      TObjString *Header2 = new TObjString("eta_2");
      HeaderList->Add(Header2);
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
      mgr->AddTask(fTrackMatcher);
      mgr->ConnectInput(fTrackMatcher,0,cinput);
    }

    analysisEventCuts[i]    = new AliConvEventCuts();   
    analysisEventCuts[i]->SetV0ReaderName(V0ReaderName);
    if (periodNameV0Reader.CompareTo("") != 0) analysisEventCuts[i]->SetPeriodEnum(periodNameV0Reader);
    
    if(periodName.CompareTo("LHC11h") && (doFlattening > 0)){
      cout << "entering the flattening loop -> searching for file: " << fileNameInputForCentFlattening.Data() << endl;  
      if( fileNameInputForCentFlattening.Contains("Low") ){
        analysisEventCuts[i]->SetUseWeightFlatCentralityFromFile(doFlattening, fileNameInputForCentFlattening, "CentLowRange");
      }else if( fileNameInputForCentFlattening.Contains("Middle") ){
        analysisEventCuts[i]->SetUseWeightFlatCentralityFromFile(doFlattening, fileNameInputForCentFlattening, "CentMiddleRange");
      }else if( fileNameInputForCentFlattening.Contains("High") ){
        analysisEventCuts[i]->SetUseWeightFlatCentralityFromFile(doFlattening, fileNameInputForCentFlattening, "CentHighRange");
      }else {
        analysisEventCuts[i]->SetUseWeightFlatCentralityFromFile(doFlattening, fileNameInputForCentFlattening, "Cent");
      }
    }

    analysisEventCuts[i]->SetLightOutput(runLightOutput);
    analysisEventCuts[i]->InitializeCutsFromCutString((cuts.GetEventCut(i)).Data());
    EventCutList->Add(analysisEventCuts[i]);
    analysisEventCuts[i]->SetFillCutHistograms("",kFALSE);
      
    analysisClusterCuts[i]  = new AliCaloPhotonCuts(isMC);
    analysisClusterCuts[i]->SetHistoToModifyAcceptance(histoAcc);
    analysisClusterCuts[i]->SetV0ReaderName(V0ReaderName);
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
  task->SetEventCutList(numberOfCuts,EventCutList);
  task->SetCaloCutList(numberOfCuts,ClusterCutList);
  task->SetMesonCutList(numberOfCuts,MesonCutList);
  task->SetDoMesonAnalysis(kTRUE);
  task->SetDoMesonQA(enableQAMesonTask);        //Attention new switch for Pi0 QA
  task->SetDoClusterQA(enableQAClusterTask);    //Attention new switch small for Cluster QA
  task->SetDoTHnSparse(isUsingTHnSparse);
  task->SetProduceTreeEOverP(doTreeEOverP);
  task->SetEnableSortingOfMCClusLabels(enableSortingMCLabels);
  if(enableExtMatchAndQA > 1){ task->SetPlotHistsExtQA(kTRUE);}

  //connect containers
  AliAnalysisDataContainer *coutput =
    mgr->CreateContainer(Form("GammaCalo_%i",trainConfig), TList::Class(),
              AliAnalysisManager::kOutputContainer,Form("GammaCalo_%i.root",trainConfig));

  mgr->AddTask(task);
  mgr->ConnectInput(task,0,cinput);
  mgr->ConnectOutput(task,1,coutput);

  return;

}
