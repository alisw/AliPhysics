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
void AddTask_GammaConvCalo_PbPb(  Int_t     trainConfig                 = 1,                      // change different set of cuts
                                  Int_t     isMC                        = 0,                      // run MC 
                                  Int_t     enableQAMesonTask           = 0,                      // enable QA in AliAnalysisTaskGammaConvV1
                                  Int_t     enableQAPhotonTask          = 0,                      // enable additional QA task
                                  TString   fileNameInputForWeighting   = "MCSpectraInput.root",  // path to file for weigting input
                                  Int_t     headerSelectionInt          = 0,                      // 1 pi0 header, 2 eta header, 3 both (only for "named" boxes)
                                  TString   cutnumberAODBranch          = "100000006008400000001500000",
                                  TString   periodName                  = "LHC13d2",              // name of the period for added signals and weighting
                                  Bool_t    doWeighting                 = kFALSE,                 // enable Weighting
                                  Int_t     enableExtMatchAndQA         = 0,                      // enable matching histograms (1) and extended QA (2), only QA(3), all disabled (0)
                                  Bool_t    isUsingTHnSparse            = kTRUE,                  // enable or disable usage of THnSparses for background estimation
                                  Bool_t    enableV0findingEffi         = kFALSE,                 // enables V0finding efficiency histograms
                                  TString   periodNameV0Reader          = "",                     // period Name for V0Reader
                                  Bool_t    enableSortingMCLabels       = kTRUE,                  // enable sorting for MC cluster labels
                                  Bool_t    runLightOutput              = kFALSE                  // switch to run light output (only essential histograms for afterburner)
                                ) {

  // ================= Load Librariers =================================
  gSystem->Load("libCore");  
  gSystem->Load("libTree");
  gSystem->Load("libGeom");
  gSystem->Load("libVMC");
  gSystem->Load("libPhysics");
  gSystem->Load("libMinuit");
  gSystem->Load("libSTEERBase");
  gSystem->Load("libESD");
  gSystem->Load("libAOD");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");  
  gSystem->Load("libCDB");
  gSystem->Load("libSTEER");
  gSystem->Load("libSTEERBase");
  gSystem->Load("libTender");
  gSystem->Load("libTenderSupplies");
  gSystem->Load("libPWGflowBase");
  gSystem->Load("libPWGflowTasks");
  gSystem->Load("libPWGGAGammaConv");

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
    if (!mgr) {
      Error("AddTask_V0ReaderV1", "No analysis manager found.");
      return;
    }

    AliConvEventCuts *fEventCuts=NULL;
    if(cutnumberEvent!=""){
      fEventCuts= new AliConvEventCuts(cutnumberEvent.Data(),cutnumberEvent.Data());
      fEventCuts->SetPreSelectionCutFlag(kTRUE);
      fEventCuts->SetV0ReaderName(V0ReaderName);
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
  AliAnalysisTaskGammaConvCalo *task=NULL;
  task= new AliAnalysisTaskGammaConvCalo(Form("GammaConvCalo_%i",trainConfig));
  task->SetIsHeavyIon(isHeavyIon);
  task->SetIsMC(isMC);
  task->SetV0ReaderName(V0ReaderName);
  task->SetLightOutput(runLightOutput);

  //create cut handler
  CutHandlerConvCalo cuts;

  // cluster cuts
  // 0 "ClusterType",  1 "EtaMin", 2 "EtaMax", 3 "PhiMin", 4 "PhiMax", 5 "DistanceToBadChannel", 6 "Timing", 7 "TrackMatching", 8 "ExoticCell",
  // 9 "MinEnergy", 10 "MinNCells", 11 "MinM02", 12 "MaxM02", 13 "MinM20", 14 "MaxM20", 15 "MaximumDispersion", 16 "NLM"
  
  if (trainConfig == 1){ // EMCAL clusters
    cuts.AddCut("60100013","00200009297002008250400000","1111100050032230000","0163103100000010"); // 0-5%
    cuts.AddCut("61200013","00200009297002008250400000","1111100050032230000","0163103100000010"); // 5-10%
    cuts.AddCut("50100013","00200009297002008250400000","1111100050032230000","0163103100000010"); // 0-10%
    cuts.AddCut("52400013","00200009297002008250400000","1111100050032230000","0163103100000010"); // 20-40%
    cuts.AddCut("52500013","00200009297002008250400000","1111100050032230000","0163103100000010"); // 20-50%
  } else if (trainConfig == 2){ // EMCAL clusters
    cuts.AddCut("60100013","00200009297002008250400000","1111100050032230000","0163103100000010"); // 0-5%
    cuts.AddCut("61200013","00200009297002008250400000","1111100050032230000","0163103100000010"); // 5-10%
    cuts.AddCut("50100013","00200009297002008250400000","1111100050032230000","0163103100000010"); // 0-10%
    cuts.AddCut("51200013","00200009297002008250400000","1111100050032230000","0163103100000010"); // 10-20%
    cuts.AddCut("52400013","00200009297002008250400000","1111100050032230000","0163103100000010"); // 20-40%
  } else if (trainConfig == 3){ // EMCAL clusters
    cuts.AddCut("54600013","00200009297002008250400000","1111100050032230000","0163103100000010"); // 40-60%
    cuts.AddCut("56800013","00200009297002008250400000","1111100050032230000","0163103100000010"); // 60-80%
    cuts.AddCut("52600013","00200009297002008250400000","1111100050032230000","0163103100000010"); // 20-60%
    cuts.AddCut("54800013","00200009297002008250400000","1111100050032230000","0163103100000010"); // 40-80%
    cuts.AddCut("52500013","00200009297002008250400000","1111100050032230000","0163103100000010"); // 20-50%
  } else if (trainConfig == 4){ // EMCAL clusters
    cuts.AddCut("60100013","00200009297002008250400000","1111100053032230000","0163103100000010"); // 0-5%
  } else if (trainConfig == 5){ // EMCAL clusters
    cuts.AddCut("61200013","00200009297002008250400000","1111100053032230000","0163103100000010"); // 5-10%
  } else if (trainConfig == 6){ // EMCAL clusters
    cuts.AddCut("50100013","00200009297002008250400000","1111100053032230000","0163103100000010"); // 0-10%
  } else if (trainConfig == 7){ // EMCAL clusters  
    cuts.AddCut("51200013","00200009297002008250400000","1111100053032230000","0163103100000010"); // 10-20%
  } else if (trainConfig == 8){ // EMCAL clusters
    cuts.AddCut("52400013","00200009297002008250400000","1111100053032230000","0163103100000010"); // 20-40%
  } else if (trainConfig == 9){ // EMCAL clusters
    cuts.AddCut("52500013","00200009297002008250400000","1111100053032230000","0163103100000010"); // 20-50%
  } else if (trainConfig == 10){ // EMCAL clusters
    cuts.AddCut("54600013","00200009297002008250400000","1111100053032230000","0163103100000010"); // 40-60%
  } else if (trainConfig == 11){ // EMCAL clusters  
    cuts.AddCut("56800013","00200009297002008250400000","1111100053032230000","0163103100000010"); // 60-80%
  } else if (trainConfig == 31){ // PHOS clusters
    cuts.AddCut("60100013","00200009297002008250400000","2444400048033200000","0163103100000010"); // 0-5%
    cuts.AddCut("61200013","00200009297002008250400000","2444400048033200000","0163103100000010"); // 5-10%
    cuts.AddCut("50100013","00200009297002008250400000","2444400048033200000","0163103100000010"); // 0-10%
    cuts.AddCut("52400013","00200009297002008250400000","2444400048033200000","0163103100000010"); // 20-40%
    cuts.AddCut("52500013","00200009297002008250400000","2444400048033200000","0163103100000010"); // 20-50%
  } else if (trainConfig == 32){ // PHOS clusters
    cuts.AddCut("60100013","00200009297002008250400000","2444400048033200000","0163103100000010"); // 0-5%
    cuts.AddCut("61200013","00200009297002008250400000","2444400048033200000","0163103100000010"); // 5-10%
    cuts.AddCut("50100013","00200009297002008250400000","2444400048033200000","0163103100000010"); // 0-10%
    cuts.AddCut("51200013","00200009297002008250400000","2444400048033200000","0163103100000010"); // 10-20%
    cuts.AddCut("52400013","00200009297002008250400000","2444400048033200000","0163103100000010"); // 20-40%
  } else if (trainConfig == 33){ // PHOS clusters
    cuts.AddCut("54600013","00200009297002008250400000","2444400048033200000","0163103100000010"); // 40-60%
    cuts.AddCut("56800013","00200009297002008250400000","2444400048033200000","0163103100000010"); // 60-80%
    cuts.AddCut("52600013","00200009297002008250400000","2444400048033200000","0163103100000010"); // 20-60%
    cuts.AddCut("54800013","00200009297002008250400000","2444400048033200000","0163103100000010"); // 40-80%
    cuts.AddCut("52500013","00200009297002008250400000","2444400048033200000","0163103100000010"); // 20-50%
  } else if (trainConfig == 34){ // PHOS clusters
    cuts.AddCut("60100013","00200009297002008250400000","2444400048033200000","0163103100000010"); // 0-5%
    cuts.AddCut("61200013","00200009297002008250400000","2444400048033200000","0163103100000010"); // 5-10%
    cuts.AddCut("50100013","00200009297002008250400000","2444400048033200000","0163103100000010"); // 0-10%
    cuts.AddCut("52400013","00200009297002008250400000","2444400048033200000","0163103100000010"); // 20-40%
    cuts.AddCut("52500013","00200009297002008250400000","2444400048033200000","0163103100000010"); // 20-50%
  } else if (trainConfig == 35){ // PHOS clusters
    cuts.AddCut("60100013","00200009297002008250400000","2444400048033200000","0163103100000010"); // 0-5%
    cuts.AddCut("61200013","00200009297002008250400000","2444400048033200000","0163103100000010"); // 5-10%
    cuts.AddCut("50100013","00200009297002008250400000","2444400048033200000","0163103100000010"); // 0-10%
    cuts.AddCut("51200013","00200009297002008250400000","2444400048033200000","0163103100000010"); // 10-20%
    cuts.AddCut("52400013","00200009297002008250400000","2444400048033200000","0163103100000010"); // 20-40%
  } else if (trainConfig == 36){ // PHOS clusters
    cuts.AddCut("54600013","00200009297002008250400000","2444400048033200000","0163103100000010"); // 40-60%
    cuts.AddCut("56800013","00200009297002008250400000","2444400048033200000","0163103100000010"); // 60-80%
    cuts.AddCut("52600013","00200009297002008250400000","2444400048033200000","0163103100000010"); // 20-60%
    cuts.AddCut("54800013","00200009297002008250400000","2444400048033200000","0163103100000010"); // 40-80%
    cuts.AddCut("52500013","00200009297002008250400000","2444400048033200000","0163103100000010"); // 20-50%
    
   } else if (trainConfig == 201){ // EMCAL clusters
    cuts.AddCut("20110113","00200009327000008250400000","1111100000022220000","0163103100000050"); // 0-10
    cuts.AddCut("21210113","00200009327000008250400000","1111100000022220000","0163103100000050"); // 10-20
    cuts.AddCut("22510113","00200009327000008250400000","1111100000022220000","0163103100000050"); // 20-50
    cuts.AddCut("25910113","00200009327000008250400000","1111100000022220000","0163103100000050"); // 50-90
    cuts.AddCut("20010113","00200009327000008250400000","1111100000022220000","0163103100000050"); // 0-100
  } else if (trainConfig == 202){ // EMCAL clusters
    cuts.AddCut("50110113","00200009327000008250400000","1111100000022220000","0163103100000050"); // 0-10
    cuts.AddCut("51210113","00200009327000008250400000","1111100000022220000","0163103100000050"); // 10-20
    cuts.AddCut("52510113","00200009327000008250400000","1111100000022220000","0163103100000050"); // 20-50
    cuts.AddCut("55910113","00200009327000008250400000","1111100000022220000","0163103100000050"); // 50-90
    cuts.AddCut("50010113","00200009327000008250400000","1111100000022220000","0163103100000050"); // 0-100
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
    } else {
      TObjString *Header1 = new TObjString("pi0_1");
      HeaderList->Add(Header1);
      TObjString *Header2 = new TObjString("eta_2");
      HeaderList->Add(Header2);
    }  
  } else if (periodName.CompareTo("LHC14a1b")==0 || periodName.CompareTo("LHC14a1c")==0){
    TObjString *Header1 = new TObjString("BOX");
    HeaderList->Add(Header1);
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
    analysisEventCuts[i]->SetLightOutput(runLightOutput);
    analysisEventCuts[i]->InitializeCutsFromCutString((cuts.GetEventCut(i)).Data());
    if (periodName.CompareTo("LHC14a1b") ==0 || periodName.CompareTo("LHC14a1c") ==0 ){
      if (headerSelectionInt == 1) analysisEventCuts[i]->SetAddedSignalPDGCode(111);
      if (headerSelectionInt == 2) analysisEventCuts[i]->SetAddedSignalPDGCode(221);
    }
    EventCutList->Add(analysisEventCuts[i]);
    analysisEventCuts[i]->SetFillCutHistograms("",kFALSE);

    analysisCuts[i] = new AliConversionPhotonCuts();
    analysisCuts[i]->SetV0ReaderName(V0ReaderName);
    analysisCuts[i]->SetLightOutput(runLightOutput);
    analysisCuts[i]->InitializeCutsFromCutString((cuts.GetPhotonCut(i)).Data());
    ConvCutList->Add(analysisCuts[i]);
    analysisCuts[i]->SetFillCutHistograms("",kFALSE);
        
    analysisClusterCuts[i] = new AliCaloPhotonCuts((isMC==2));
    analysisClusterCuts[i]->SetV0ReaderName(V0ReaderName);
    analysisClusterCuts[i]->SetLightOutput(runLightOutput);
    analysisClusterCuts[i]->InitializeCutsFromCutString((cuts.GetClusterCut(i)).Data());
    ClusterCutList->Add(analysisClusterCuts[i]);
    analysisClusterCuts[i]->SetExtendedMatchAndQA(enableExtMatchAndQA);
    analysisClusterCuts[i]->SetFillCutHistograms("");
    
    analysisMesonCuts[i] = new AliConversionMesonCuts();
    analysisMesonCuts[i]->SetLightOutput(runLightOutput);
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
  task->SetDoMesonAnalysis(kTRUE);
  task->SetDoMesonQA(enableQAMesonTask); //Attention new switch for Pi0 QA
  task->SetDoPhotonQA(enableQAPhotonTask);  //Attention new switch small for Photon QA
  task->SetDoClusterQA(1);  //Attention new switch small for Cluster QA
  task->SetEnableSortingOfMCClusLabels(enableSortingMCLabels);
    task->SetUseTHnSparse(isUsingTHnSparse);
  if(enableExtMatchAndQA == 2 || enableExtMatchAndQA == 3){ task->SetPlotHistsExtQA(kTRUE);}
  
  //connect containers
  AliAnalysisDataContainer *coutput =
    mgr->CreateContainer(Form("GammaConvCalo_%i",trainConfig), TList::Class(),
              AliAnalysisManager::kOutputContainer,Form("GammaConvCalo_%i.root",trainConfig));

  mgr->AddTask(task);
  mgr->ConnectInput(task,0,cinput);
  mgr->ConnectOutput(task,1,coutput);

  return;

}
