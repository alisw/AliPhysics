/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: Friederike Bock, Lucia Leardini                                *
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
//($ALIPHYSICS/PWGGA/GammaConv/AliAnalysisTaskGammaConvV1.cxx) for
//pp together with all supporting classes
//***************************************************************************************

//***************************************************************************************
//CutHandler contains all cuts for a certain analysis and trainconfig,
//it automatically checks length of cutStrings and takes care of the number of added cuts,
//no specification of the variable 'numberOfCuts' needed anymore.
//***************************************************************************************

class CutHandlerConv{
  public:
    CutHandlerConv(Int_t nMax=10){
      nCuts=0; nMaxCuts=nMax; validCuts = true;
      eventCutArray = new TString[nMaxCuts]; photonCutArray = new TString[nMaxCuts]; mesonCutArray = new TString[nMaxCuts]; clusterCutArray = new TString[nMaxCuts];
      for(Int_t i=0; i<nMaxCuts; i++) {eventCutArray[i] = ""; photonCutArray[i] = ""; mesonCutArray[i] = ""; clusterCutArray[i] = "";}
    }

    void AddCut(TString eventCut, TString photonCut, TString mesonCut){
      if(nCuts>=nMaxCuts) {cout << "ERROR in CutHandlerConv: Exceeded maximum number of cuts!" << endl; validCuts = false; return;}
      if( eventCut.Length()!=8 || photonCut.Length()!=26 || mesonCut.Length()!=16 ) {cout << "ERROR in CutHandlerConv: Incorrect length of cut string!" << endl; validCuts = false; return;}
      eventCutArray[nCuts]=eventCut; photonCutArray[nCuts]=photonCut; mesonCutArray[nCuts]=mesonCut; clusterCutArray[nCuts]="";
      nCuts++;
      return;
    }
    void AddCut(TString eventCut, TString photonCut, TString mesonCut, TString clusterCut){
      if(nCuts>=nMaxCuts) {cout << "ERROR in CutHandlerConv: Exceeded maximum number of cuts!" << endl; validCuts = false; return;}
      if( eventCut.Length()!=8 || photonCut.Length()!=26 || mesonCut.Length()!=16 || clusterCut.Length()!=19 ) {cout << "ERROR in CutHandlerConv: Incorrect length of cut string!" << endl; validCuts = false; return;}
      eventCutArray[nCuts]=eventCut; photonCutArray[nCuts]=photonCut; mesonCutArray[nCuts]=mesonCut; clusterCutArray[nCuts]=clusterCut;
      nCuts++;
      return;
    }

    Bool_t AreValid(){return validCuts;}
    Int_t GetNCuts(){if(validCuts) return nCuts; else return 0;}
    TString GetEventCut(Int_t i){if(validCuts&&i<nMaxCuts&&i>=0) return eventCutArray[i]; else{cout << "ERROR in CutHandlerConv: GetEventCut wrong index i" << endl;return "";}}
    TString GetPhotonCut(Int_t i){if(validCuts&&i<nMaxCuts&&i>=0) return photonCutArray[i]; else {cout << "ERROR in CutHandlerConv: GetPhotonCut wrong index i" << endl;return "";}}
    TString GetMesonCut(Int_t i){if(validCuts&&i<nMaxCuts&&i>=0) return mesonCutArray[i]; else {cout << "ERROR in CutHandlerConv: GetMesonCut wrong index i" << endl;return "";}}
    TString GetClusterCut(Int_t i){if(validCuts&&i<nMaxCuts&&i>=0) return clusterCutArray[i]; else {cout << "ERROR in CutHandlerConv: GetClusterCut wrong index i" << endl;return "";}}
  private:
    Bool_t validCuts;
    Int_t nCuts; Int_t nMaxCuts;
    TString* eventCutArray;
    TString* photonCutArray;
    TString* mesonCutArray;
    TString* clusterCutArray;
};

void AddTask_GammaConvV1_pp(  Int_t   trainConfig                     = 1,                    // change different set of cuts
                              Int_t   isMC                            = 0,                    // run MC
                              Int_t   enableQAMesonTask               = 0,                    // enable meson QA in AliAnalysisTaskGammaConvV1
                              Int_t   enableQAPhotonTask              = 0,                    // enable photon QA in AliAnalysisTaskGammaConvV1
                              TString fileNameInputForPartWeighting   = "MCSpectraInput.root",// path to file for weigting input
                              TString cutnumberAODBranch              = "000000006008400001001500000",  // cutnumber for AOD branch
                              TString periodname                      = "LHC12f1x",           // period name
                              Bool_t  doParticleWeighting             = kFALSE,               // enables weighting
                              Bool_t  enableV0findingEffi             = kFALSE,               // enables V0finding efficiency histograms
                              Bool_t  enableTriggerMimicking          = kFALSE,               // enable trigger mimicking
                              Bool_t  enableTriggerOverlapRej         = kFALSE,               // enable trigger overlap rejection
                              Float_t maxFacPtHard                    = 3.,                   // maximum factor between hardest jet and ptHard generated
                              TString periodNameV0Reader              = "",
                              Bool_t  doMultiplicityWeighting         = kFALSE,                  //
                              TString fileNameInputForMultWeighing    = "Multiplicity.root",    //
                              TString periodNameAnchor                = ""
                              
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

  Int_t isHeavyIon = 0;
  // ================== GetAnalysisManager ===============================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error(Form("AddTask_GammaConvV1_%i",trainConfig), "No analysis manager found.");
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
  TString cutnumberPhoton = "00200008400000002200000000";
  TString cutnumberEvent = "00000003"; 
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
    if(trainConfig == 101){
      fV0ReaderV1->SetProduceImpactParamHistograms(kTRUE);
      cutnumberPhoton="00200009227302008250400400";
      cutnumberEvent ="00010113";
    }
    if (!mgr) {
      Error("AddTask_V0ReaderV1", "No analysis manager found.");
      return;
    }

    AliConvEventCuts *fEventCuts=NULL;
    if(cutnumberEvent!=""){
      fEventCuts= new AliConvEventCuts(cutnumberEvent.Data(),cutnumberEvent.Data());
      fEventCuts->SetPreSelectionCutFlag(kTRUE);
      fEventCuts->SetV0ReaderName(V0ReaderName);
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
      if(fCuts->InitializeCutsFromCutString(cutnumberPhoton.Data())){
        fV0ReaderV1->SetConversionCuts(fCuts);
        fCuts->SetFillCutHistograms("",kTRUE);
      }
      if( trainConfig  == 73 || trainConfig  == 74 || (trainConfig  >= 80 && trainConfig  <= 87)   ){
        fCuts->SetDodEdxSigmaCut(kFALSE);
      }
    }

    
    if(inputHandler->IsA()==AliAODInputHandler::Class()){
    // AOD mode
      fV0ReaderV1->SetDeltaAODBranchName(Form("GammaConv_%s_gamma",cutnumberAODBranch.Data()));
    }
    fV0ReaderV1->Init();

    AliLog::SetGlobalLogLevel(AliLog::kInfo);

    //connect input V0Reader
    mgr->AddTask(fV0ReaderV1);
    mgr->ConnectInput(fV0ReaderV1,0,cinput);

  }

  //================================================
  //========= Add task to the ANALYSIS manager =====
  //            find input container
  AliAnalysisTaskGammaConvV1 *task=NULL;
  task= new AliAnalysisTaskGammaConvV1(Form("GammaConvV1_%i",trainConfig));
  task->SetIsHeavyIon(isHeavyIon);
  task->SetIsMC(isMC);
  task->SetV0ReaderName(V0ReaderName);
  // Cut Numbers to use in Analysis
  
  CutHandlerConv cuts;
  
  if(trainConfig == 1){
    cuts.AddCut("00000113", "00200009366300003800000000", "0163103100900000"); //standard cut Pi0 pp 2.76TeV without SDD , only Minbias MC
    cuts.AddCut("00012113", "00200009366300003800000000", "0163103100900000"); //standard cut Pi0 pp 2.76TeV without SDD, V0AND
    cuts.AddCut("00000113", "00200009326000003800000000", "0163103100900000"); //standard cut Gamma pp 2-76TeV
    cuts.AddCut("00000113", "00200009366000003800000000", "0163103100900000"); //standard cut Gamma pp 2-76TeV
  } else if (trainConfig == 2) {
    cuts.AddCut("00000123", "00200009366300003800000000", "0163103100900000"); //standard cut Pi0 pp 2.76TeV without SDD , only boxes
    cuts.AddCut("00012123", "00200009366300003800000000", "0163103100900000"); //standard cut Pi0 pp 2.76TeV without SDD, V0AND , only boxes
    cuts.AddCut("00000123", "00200009326000003800000000", "0163103100900000"); //standard cut Gamma pp 2-76TeV , only boxes
    cuts.AddCut("00000123", "00200009366000003800000000", "0163103100900000"); //standard cut Gamma pp 2-76TeV 
  } else if (trainConfig == 3) {
    cuts.AddCut("00003113", "00200009366300003800000000", "0163103100900000"); //standard cut Pi0 pp 2.76TeV with SDD , only Minbias MC
    cuts.AddCut("00013113", "00200009366300003800000000", "0163103100900000"); //standard cut Pi0 pp 2.76TeV with SDD, V0AND , only Minbias MC
    cuts.AddCut("00003123", "00200009366300003800000000", "0163103100900000"); //standard cut Pi0 pp 2.76TeV with SDD , only Boxes MC
    cuts.AddCut("00013123", "00200009366300003800000000", "0163103100900000"); //standard cut Pi0 pp 2.76TeV with SDD, V0AND, only Boxes MC
  } else if (trainConfig == 4) {
    cuts.AddCut("00000113", "00200009366300003800000000", "0163103100900000"); //standard cut Pi0 pp 2.76TeV without SDD , all photon qualities
    cuts.AddCut("00000113", "00200009366300003800020000", "0163103100900000"); //standard cut Pi0 pp 2.76TeV without SDD , photon quality 1
    cuts.AddCut("00000113", "00200009366300003800030000", "0163103100900000"); //standard cut Pi0 pp 2.76TeV without SDD , photon quality 2
    cuts.AddCut("00000113", "00200009366300003800040000", "0163103100900000"); //standard cut Pi0 pp 2.76TeV without SDD , photon quality 3
  } else if (trainConfig == 5) {
    cuts.AddCut("00000113", "00700009366300003800000000", "0163103100900000"); //standard cut Pi0 pp 2.76TeV without SDD , all photon qualities, min R = 35 cm
    cuts.AddCut("00000113", "00700009366300003800020000", "0163103100900000"); //standard cut Pi0 pp 2.76TeV without SDD , photon quality 1, min R = 35 cm
    cuts.AddCut("00000113", "00700009366300003800030000", "0163103100900000"); //standard cut Pi0 pp 2.76TeV without SDD , photon quality 2, min R = 35 cm
    cuts.AddCut("00000113", "00700009366300003800040000", "0163103100900000"); //standard cut Pi0 pp 2.76TeV without SDD , photon quality 3, min R = 35 cm
  } else if (trainConfig == 6) {
    cuts.AddCut("00000113", "00200008366300003200000000", "0163103100900000"); //standard cut Pi0 pp 7TeV, with qt cut 0.05, all photon qualities
    cuts.AddCut("00000113", "00200008366300003200020000", "0163103100900000"); //standard cut Pi0 pp 7TeV, with qt cut 0.05, photon quality 1
    cuts.AddCut("00000113", "00200008366300003200030000", "0163103100900000"); //standard cut Pi0 pp 7TeV, with qt cut 0.05, photon quality 2
    cuts.AddCut("00000113", "00200008366300003200040000", "0163103100900000"); //standard cut Pi0 pp 7TeV, with qt cut 0.05, photon quality 3   
  } else if (trainConfig == 7) {
    cuts.AddCut("00000113", "00700008366300003200000000", "0163103100900000"); //standard cut Pi0 pp 7TeV, with qt cut 0.05, all photon qualities, min R = 35 cm
    cuts.AddCut("00000113", "00700008366300003200020000", "0163103100900000"); //standard cut Pi0 pp 7TeV, with qt cut 0.05, photon quality 1, min R = 35 cm
    cuts.AddCut("00000113", "00700008366300003200030000", "0163103100900000"); //standard cut Pi0 pp 7TeV, with qt cut 0.05, photon quality 2, min R = 35 cm
    cuts.AddCut("00000113", "00700008366300003200040000", "0163103100900000"); //standard cut Pi0 pp 7TeV, with qt cut 0.05, photon quality 3, min R = 35 cm
  } else if (trainConfig == 8) {
    cuts.AddCut("00000113", "00200008366300000200000000", "0163103100900000"); //standard cut Pi0 pp 7TeV, all photon qualities
    cuts.AddCut("00000113", "00200008366300000200020000", "0163103100900000"); //standard cut Pi0 pp 7TeV, photon quality 1
    cuts.AddCut("00000113", "00200008366300000200030000", "0163103100900000"); //standard cut Pi0 pp 7TeV, photon quality 2
    cuts.AddCut("00000113", "00200008366300000200040000", "0163103100900000"); //standard cut Pi0 pp 7TeV, photon quality 3   
  } else if (trainConfig == 9) {
    cuts.AddCut("00000113", "00700008366300000200000000", "0163103100900000"); //standard cut Pi0 pp 7TeV, all photon qualities, min R = 35 cm
    cuts.AddCut("00000113", "00700008366300000200020000", "0163103100900000"); //standard cut Pi0 pp 7TeV, photon quality 1, min R = 35 cm
    cuts.AddCut("00000113", "00700008366300000200030000", "0163103100900000"); //standard cut Pi0 pp 7TeV, photon quality 2, min R = 35 cm
    cuts.AddCut("00000113", "00700008366300000200040000", "0163103100900000"); //standard cut Pi0 pp 7TeV, photon quality 3, min R = 35 cm	   
  } else if (trainConfig == 10) {
    cuts.AddCut("00003113", "00200009366300003800000000", "0163103100900000"); //standard cut Pi0 pp 2.76TeV with SDD , all photon qualities
    cuts.AddCut("00003113", "00200009366300003800020000", "0163103100900000"); //standard cut Pi0 pp 2.76TeV with SDD , photon quality 1
    cuts.AddCut("00003113", "00200009366300003800030000", "0163103100900000"); //standard cut Pi0 pp 2.76TeV with SDD , photon quality 2
    cuts.AddCut("00003113", "00200009366300003800040000", "0163103100900000"); //standard cut Pi0 pp 2.76TeV with SDD , photon quality 3
  } else if (trainConfig == 11) {
    cuts.AddCut("00003113", "00700009366300003800000000", "0163103100900000"); //standard cut Pi0 pp 2.76TeV with SDD , all photon qualities, min R = 35 cm
    cuts.AddCut("00003113", "00700009366300003800020000", "0163103100900000"); //standard cut Pi0 pp 2.76TeV with SDD , photon quality 1, min R = 35 cm
    cuts.AddCut("00003113", "00700009366300003800030000", "0163103100900000"); //standard cut Pi0 pp 2.76TeV with SDD , photon quality 2, min R = 35 cm
    cuts.AddCut("00003113", "00700009366300003800040000", "0163103100900000"); //standard cut Pi0 pp 2.76TeV with SDD , photon quality 3, min R = 35 cm
  } else if (trainConfig == 12) {
    cuts.AddCut("00000113", "00200009297002008250400000", "0152506500000000"); //standard cut LHC11h pp 2.76TeV 
    cuts.AddCut("00000113", "03200009297002008250400000", "0152506500000000"); //variation eta 0.65
    cuts.AddCut("00000113", "04200009297002008250400000", "0152506500000000"); //variation eta 0.75
    cuts.AddCut("00000113", "00200009295002008250400000", "0152506500000000"); //variation pion p dEdx 0.3-5.
  } else if (trainConfig == 13) { //added signals
    cuts.AddCut("00000123", "00200009297002008250400000", "0152506500000000"); //standard cut LHC11h pp 2.76TeV 
    cuts.AddCut("00000123", "03200009297002008250400000", "0152506500000000"); //variation eta 0.65
    cuts.AddCut("00000123", "04200009297002008250400000", "0152506500000000"); //variation eta 0.75
    cuts.AddCut("00000123", "00200009295002008250400000", "0152506500000000"); //variation pion p dEdx 0.3-5.
  } else if (trainConfig == 14) {
    cuts.AddCut("00000113", "00200049297002008250400000", "0152506500000000"); //variation pt 0.075 
    cuts.AddCut("00000113", "00200019297002008250400000", "0152506500000000"); //variation pt 0.1
    cuts.AddCut("00000113", "00200006297002008250400000", "0152506500000000"); //variation TPC cls 0.7
    cuts.AddCut("00000113", "00200008297002008250400000", "0152506500000000"); //variation TPC cls 0.35 
  } else if (trainConfig == 15) { //added signals
    cuts.AddCut("00000123", "00200049297002008250400000", "0152506500000000"); //variation pt 0.075 
    cuts.AddCut("00000123", "00200019297002008250400000", "0152506500000000"); //variation pt 0.1
    cuts.AddCut("00000123", "00200006297002008250400000", "0152506500000000"); //variation TPC cls 0.7
    cuts.AddCut("00000123", "00200008297002008250400000", "0152506500000000"); //variation TPC cls 0.35 
  } else if (trainConfig == 16) {
    cuts.AddCut("00000113", "00200009397002008250400000", "0152506500000000"); //variation edEdx -4,5
    cuts.AddCut("00000113", "00200009697002008250400000", "0152506500000000"); //variation edEdx -2.5,4
    cuts.AddCut("00000113", "00200009297003008250400000", "0152506500000000"); //variation TOF el. PID -3,5
    cuts.AddCut("00000113", "00200009297004008250400000", "0152506500000000"); //variation TOF el. PID -2,3
  } else if (trainConfig == 17) { //added signals
    cuts.AddCut("00000123", "00200009397002008250400000", "0152506500000000"); //variation edEdx -4,5
    cuts.AddCut("00000123", "00200009697002008250400000", "0152506500000000"); //variation edEdx -2.5,4
    cuts.AddCut("00000123", "00200009297003008250400000", "0152506500000000"); //variation TOF el. PID -3,5
    cuts.AddCut("00000123", "00200009297004008250400000", "0152506500000000"); //variation TOF el. PID -2,3
  } else if (trainConfig == 18) {
    cuts.AddCut("00000113", "00200009297002009250400000", "0152506500000000"); //variation qt 0.03
    cuts.AddCut("00000113", "00200009297002002250400000", "0152506500000000"); //variation qt 0.07 no2D
    cuts.AddCut("00000113", "00200009297002008150400000", "0152506500000000"); //variation chi2 50.
    cuts.AddCut("00000113", "00200009297002008850400000", "0152506500000000"); //variation chi2 20.
  } else if (trainConfig == 19) { //added signals
    cuts.AddCut("00000123", "00200009297002009250400000", "0152506500000000"); //variation qt 0.03
    cuts.AddCut("00000123", "00200009297002002250400000", "0152506500000000"); //variation qt 0.07 no2D
    cuts.AddCut("00000123", "00200009297002008150400000", "0152506500000000"); //variation chi2 50.
    cuts.AddCut("00000123", "00200009297002008850400000", "0152506500000000"); //variation chi2 20.
  } else if (trainConfig == 20) {
    cuts.AddCut("00000113", "00200009297002008260400000", "0152506500000000"); //variation psi pair 0.05
    cuts.AddCut("00000113", "00200009297002008280400000", "0152506500000000"); //variation psi pair 0.2
    cuts.AddCut("00000113", "00200009297002008250000000", "0152506500000000"); //variation cosPA -1
    cuts.AddCut("00000113", "00200009297002008250400000", "0152505500000000"); //variation alpha 0.75
  } else if (trainConfig == 21) { //added signals
    cuts.AddCut("00000123", "00200009297002008260400000", "0152506500000000"); //variation psi pair 0.05
    cuts.AddCut("00000123", "00200009297002008280400000", "0152506500000000"); //variation psi pair 0.2
    cuts.AddCut("00000123", "00200009297002008250000000", "0152506500000000"); //variation cosPA -1
    cuts.AddCut("00000123", "00200009297002008250400000", "0152505500000000"); //variation alpha 0.75
  } else if (trainConfig == 22) {
    cuts.AddCut("00040113", "00200009297002008250400000", "0152506500000000"); // trigger kTRD
    cuts.AddCut("00050113", "00200009297002008250400000", "0152506500000000"); // trigger kEMC
    cuts.AddCut("00060113", "00200009297002008250400000", "0152506500000000"); // trigger kPHI
    cuts.AddCut("00070113", "00200009297002008250400000", "0152506500000000"); // trigger kHighMult
  } else if (trainConfig == 23) {
    cuts.AddCut("00080113", "00200009297002008250400000", "0152506500000000"); // trigger kEMCEGA
    cuts.AddCut("00090113", "00200009297002008250400000", "0152506500000000"); // trigger kEMCEJE
    cuts.AddCut("00000113", "00200009297002008250400000", "0152506500000000"); // minimum bias
    cuts.AddCut("00011113", "00200009297002008250400000", "0152506500000000"); // trigger kINT8
  } else if (trainConfig == 24) {
    cuts.AddCut("00042113", "00200009297002008250400000", "0152506500000000"); // trigger kTRD CINT8 HEE
    cuts.AddCut("00044113", "00200009297002008250400000", "0152506500000000"); // trigger kTRD CINT8 HSE
    cuts.AddCut("00046113", "00200009297002008250400000", "0152506500000000"); // trigger kTRD CINT8 HJE
    cuts.AddCut("00048113", "00200009297002008250400000", "0152506500000000"); // trigger kTRD CINT8 HQU
  } else if (trainConfig == 25) {
    cuts.AddCut("00041113", "00200009297002008250400000", "0152506500000000"); // trigger kTRD CINT7 HEE
    cuts.AddCut("00043113", "00200009297002008250400000", "0152506500000000"); // trigger kTRD CINT7 HSE
    cuts.AddCut("00045113", "00200009297002008250400000", "0152506500000000"); // trigger kTRD CINT7 HJE
    cuts.AddCut("00047113", "00200009297002008250400000", "0152506500000000"); // trigger kTRD CINT7 HQU
  } else if (trainConfig == 26) {
    cuts.AddCut("00052113", "00200009297002008250400000", "0152506500000000"); // trigger kEMC7
    cuts.AddCut("00053113", "00200009297002008250400000", "0152506500000000"); // trigger kEMC8
    cuts.AddCut("00062113", "00200009297002008250400000", "0152506500000000"); // trigger kPHI7
    cuts.AddCut("00063113", "00200009297002008250400000", "0152506500000000"); // trigger kPHI8
  } else if (trainConfig == 27) {
    cuts.AddCut("00051113", "00200009297002008250400000", "0152506500000000"); // trigger kEMC1
    cuts.AddCut("00071113", "00200009297002008250400000", "0152506500000000"); // trigger kSHM1
    cuts.AddCut("00072113", "00200009297002008250400000", "0152506500000000"); // trigger kSHM7
    cuts.AddCut("00073113", "00200009297002008250400000", "0152506500000000"); // trigger kSHM8
  } else if (trainConfig == 28) {
    cuts.AddCut("00081113", "00200009297002008250400000", "0152506500000000"); // trigger kEMCEGA + CINT7
    cuts.AddCut("00082113", "00200009297002008250400000", "0152506500000000"); // trigger kEMCEGA + CINT8
    cuts.AddCut("00083113", "00200009297002008250400000", "0152506500000000"); // trigger kEMCEG1 + CINT7
    cuts.AddCut("00084113", "00200009297002008250400000", "0152506500000000"); // trigger kEMCEG1 + CINT8
  } else if (trainConfig == 29) {
    cuts.AddCut("00085113", "00200009297002008250400000", "0152506500000000"); // trigger kEMCEG2 + CINT7
    cuts.AddCut("00086113", "00200009297002008250400000", "0152506500000000"); // trigger kEMCEG2 + CINT8
    cuts.AddCut("00091113", "00200009297002008250400000", "0152506500000000"); // trigger kEMCEJE + CINT7
    cuts.AddCut("00092113", "00200009297002008250400000", "0152506500000000"); // trigger kEMCEJE + CINT8
  } else if (trainConfig == 30) {
    cuts.AddCut("00093113", "00200009297002008250400000", "0152506500000000"); // trigger kEMCEJ1 + CINT7
    cuts.AddCut("00094113", "00200009297002008250400000", "0152506500000000"); // trigger kEMCEJ1 + CINT8
    cuts.AddCut("00095113", "00200009297002008250400000", "0152506500000000"); // trigger kEMCEJ2 + CINT7
    cuts.AddCut("00096113", "00200009297002008250400000", "0152506500000000"); // trigger kEMCEJ2 + CINT8		
  } else if (trainConfig == 31) {
    cuts.AddCut("00000113", "00200009257002008250400000", "0152106500000000"); //new standard cut for pp 8 TeV
    cuts.AddCut("00000113", "00200009357002008250400000", "0152106500000000"); //variation edEdx -4,5
    cuts.AddCut("00000113", "00200009657002008250400000", "0152106500000000"); //variation edEdx -2.5,4
    cuts.AddCut("00000113", "00200009255002008250400000", "0152106500000000"); //variation pion p dEdx 0.3-5.
  } else if (trainConfig == 32) { //added signals
    cuts.AddCut("00000123", "00200009257002008250400000", "0152106500000000"); //new standard cut for pp 8 TeV
    cuts.AddCut("00000123", "00200009357002008250400000", "0152106500000000"); //variation edEdx -4,5
    cuts.AddCut("00000123", "00200009657002008250400000", "0152106500000000"); //variation edEdx -2.5,4
    cuts.AddCut("00000123", "00200009255002008250400000", "0152106500000000"); //variation pion p dEdx 0.3-5.
  } else if (trainConfig == 33) {
    cuts.AddCut("00000113", "00200049257002008250400000", "0152106500000000"); //variation pt 0.075 
    cuts.AddCut("00000113", "00200019257002008250400000", "0152106500000000"); //variation pt 0.1
    cuts.AddCut("00000113", "00200006257002008250400000", "0152106500000000"); //variation TPC cls 0.7
    cuts.AddCut("00000113", "00200008257002008250400000", "0152106500000000"); //variation TPC cls 0.35 
  } else if (trainConfig == 34) { //added signals
    cuts.AddCut("00000123", "00200049257002008250400000", "0152106500000000"); //variation pt 0.075 
    cuts.AddCut("00000123", "00200019257002008250400000", "0152106500000000"); //variation pt 0.1
    cuts.AddCut("00000123", "00200006257002008250400000", "0152106500000000"); //variation TPC cls 0.7
    cuts.AddCut("00000123", "00200008257002008250400000", "0152106500000000"); //variation TPC cls 0.35 
  } else if (trainConfig == 35) {
    cuts.AddCut("00000113", "00200009227002008250400000", "0152106500000000"); //variation pidEdx 1,-10
    cuts.AddCut("00000113", "00200009237002008250400000", "0152106500000000"); //variation pidEdx 2.5,-10
    cuts.AddCut("00000113", "00200009297002008250400000", "0152106500000000"); //variation pidEdx 3,-10
    cuts.AddCut("00000113", "00200009250002008250400000", "0152106500000000"); //variation pion p dEdx 0.5-5
  } else if (trainConfig == 36) { //added signals
    cuts.AddCut("00000123", "00200009227002008250400000", "0152106500000000"); //variation pidEdx 1,-10
    cuts.AddCut("00000123", "00200009237002008250400000", "0152106500000000"); //variation pidEdx 2.5,-10
    cuts.AddCut("00000123", "00200009297002008250400000", "0152106500000000"); //variation pidEdx 3,-10
    cuts.AddCut("00000123", "00200009250002008250400000", "0152106500000000"); //variation pion p dEdx 0.5-5
  } else if (trainConfig == 37) {
    cuts.AddCut("00000113", "00200009257002009250400000", "0152106500000000"); //variation qt 0.03
    cuts.AddCut("00000113", "00200009257002002250400000", "0152106500000000"); //variation qt 0.07 no2D
    cuts.AddCut("00000113", "00200009257002008150400000", "0152106500000000"); //variation chi2 50.
    cuts.AddCut("00000113", "00200009257002008850400000", "0152106500000000"); //variation chi2 20.
  } else if (trainConfig == 38) { //added signals
    cuts.AddCut("00000123", "00200009257002009250400000", "0152106500000000"); //variation qt 0.03
    cuts.AddCut("00000123", "00200009257002002250400000", "0152106500000000"); //variation qt 0.07 no2D
    cuts.AddCut("00000123", "00200009257002008150400000", "0152106500000000"); //variation chi2 50.
    cuts.AddCut("00000123", "00200009257002008850400000", "0152106500000000"); //variation chi2 20.
  } else if (trainConfig == 39) {
    cuts.AddCut("00000113", "00200009257002008260400000", "0152106500000000"); //variation psi pair 0.05
    cuts.AddCut("00000113", "00200009257002008280400000", "0152106500000000"); //variation psi pair 0.2
    cuts.AddCut("00000113", "00200009257002008250000000", "0152106500000000"); //variation cosPA -1
    cuts.AddCut("00000113", "00200009257002008250600000", "0152106500000000"); //variation cosPA 0.9
  } else if (trainConfig == 40) { //added signals
    cuts.AddCut("00000123", "00200009257002008260400000", "0152106500000000"); //variation psi pair 0.05
    cuts.AddCut("00000123", "00200009257002008280400000", "0152106500000000"); //variation psi pair 0.2
    cuts.AddCut("00000123", "00200009257002008250000000", "0152106500000000"); //variation cosPA -1
    cuts.AddCut("00000123", "00200009257002008250600000", "0152106500000000"); //variation cosPA 0.9
  } else if (trainConfig == 41) {
    cuts.AddCut("00000113", "00200009257002008950400000", "0152106500000000"); //variation chi2 15
    cuts.AddCut("00000113", "00200009257002008230400000", "0152106500000000"); //variation psi pair 0.035
    cuts.AddCut("00000113", "00200009257002008250400000", "0152105500000000"); //variation alpha 0.75
    cuts.AddCut("00000113", "00200009257002008250400000", "0152107500000000"); //variation alpha 0.85
  } else if (trainConfig == 42) { //added signals
    cuts.AddCut("00000123", "00200009257002008950400000", "0152106500000000"); //variation chi2 15
    cuts.AddCut("00000123", "00200009257002008230400000", "0152106500000000"); //variation psi pair 0.035
    cuts.AddCut("00000123", "00200009257002008250400000", "0152105500000000"); //variation alpha 0.75
    cuts.AddCut("00000123", "00200009257002008250400000", "0152107500000000"); //variation alpha 0.85
  } else if (trainConfig == 43) {
    cuts.AddCut("00040113", "00200009257002008250400000", "0152106500000000"); // trigger kTRD with y 0.8
    cuts.AddCut("00050113", "00200009257002008250400000", "0152106500000000"); // trigger kEMC with y 0.8
    cuts.AddCut("00060113", "00200009257002008250400000", "0152106500000000"); // trigger kPHI with y 0.8
    cuts.AddCut("00070113", "00200009257002008250400000", "0152106500000000"); // trigger kHighMult with y 0.8
  } else if (trainConfig == 44) {
    cuts.AddCut("00080113", "00200009257002008250400000", "0152106500000000"); // trigger kEMCEGA with y 0.8
    cuts.AddCut("00090113", "00200009257002008250400000", "0152106500000000"); // trigger kEMCEJE with y 0.8
    cuts.AddCut("00010113", "00200009257002008250400000", "0152106500000000"); // minimum bias with y 0.8
    cuts.AddCut("00011113", "00200009257002008250400000", "0152106500000000"); // trigger kINT8 with y 0.8
  } else if (trainConfig == 45) {
    cuts.AddCut("00042113", "00200009257002008250400000", "0152106500000000"); // trigger kTRD CINT8 HEE
    cuts.AddCut("00044113", "00200009257002008250400000", "0152106500000000"); // trigger kTRD CINT8 HSE
    cuts.AddCut("00046113", "00200009257002008250400000", "0152106500000000"); // trigger kTRD CINT8 HJE
    cuts.AddCut("00048113", "00200009257002008250400000", "0152106500000000"); // trigger kTRD CINT8 HQU
  } else if (trainConfig == 46) {
    cuts.AddCut("00041113", "00200009257002008250400000", "0152106500000000"); // trigger kTRD CINT7 HEE
    cuts.AddCut("00043113", "00200009257002008250400000", "0152106500000000"); // trigger kTRD CINT7 HSE
    cuts.AddCut("00045113", "00200009257002008250400000", "0152106500000000"); // trigger kTRD CINT7 HJE
    cuts.AddCut("00047113", "00200009257002008250400000", "0152106500000000"); // trigger kTRD CINT7 HQU
  } else if (trainConfig == 47) {
    cuts.AddCut("00052113", "00200009257002008250400000", "0152106500000000"); // trigger kEMC7
    cuts.AddCut("00053113", "00200009257002008250400000", "0152106500000000"); // trigger kEMC8
    cuts.AddCut("00062113", "00200009257002008250400000", "0152106500000000"); // trigger kPHI7
    cuts.AddCut("00063113", "00200009257002008250400000", "0152106500000000"); // trigger kPHI8
  } else if (trainConfig == 48) {
    cuts.AddCut("00051113", "00200009257002008250400000", "0152106500000000"); // trigger kEMC1
    cuts.AddCut("00071113", "00200009257002008250400000", "0152106500000000"); // trigger kSHM1
    cuts.AddCut("00072113", "00200009257002008250400000", "0152106500000000"); // trigger kSHM7
    cuts.AddCut("00073113", "00200009257002008250400000", "0152106500000000"); // trigger kSHM8
  } else if (trainConfig == 49) {
    cuts.AddCut("00081113", "00200009257002008250400000", "0152106500000000"); // trigger kEMCEGA + CINT7
    cuts.AddCut("00082113", "00200009257002008250400000", "0152106500000000"); // trigger kEMCEGA + CINT8
    cuts.AddCut("00083113", "00200009257002008250400000", "0152106500000000"); // trigger kEMCEG1 + CINT7
    cuts.AddCut("00084113", "00200009257002008250400000", "0152106500000000"); // trigger kEMCEG1 + CINT8
  } else if (trainConfig == 50) {
    cuts.AddCut("00085113", "00200009257002008250400000", "0152106500000000"); // trigger kEMCEG2 + CINT7
    cuts.AddCut("00086113", "00200009257002008250400000", "0152106500000000"); // trigger kEMCEG2 + CINT8
    cuts.AddCut("00091113", "00200009257002008250400000", "0152106500000000"); // trigger kEMCEJE + CINT7
    cuts.AddCut("00092113", "00200009257002008250400000", "0152106500000000"); // trigger kEMCEJE + CINT8
  } else if (trainConfig == 51) {
    cuts.AddCut("00093113", "00200009257002008250400000", "0152106500000000"); // trigger kEMCEJ1 + CINT7
    cuts.AddCut("00094113", "00200009257002008250400000", "0152106500000000"); // trigger kEMCEJ1 + CINT8
    cuts.AddCut("00095113", "00200009257002008250400000", "0152106500000000"); // trigger kEMCEJ2 + CINT7
    cuts.AddCut("00096113", "00200009257002008250400000", "0152106500000000"); // trigger kEMCEJ2 + CINT8		
  } else if (trainConfig == 52) { //pp 2.76TeV cuts
    cuts.AddCut("00042113", "00200009366300003800000000", "0163103100900000"); // trigger kTRD CINT8 HEE
    cuts.AddCut("00044113", "00200009366300003800000000", "0163103100900000"); // trigger kTRD CINT8 HSE
    cuts.AddCut("00046113", "00200009366300003800000000", "0163103100900000"); // trigger kTRD CINT8 HJE
    cuts.AddCut("00048113", "00200009366300003800000000", "0163103100900000"); // trigger kTRD CINT8 HQU
  } else if (trainConfig == 53) { //pp 2.76TeV cuts
    cuts.AddCut("00041113", "00200009366300003800000000", "0163103100900000"); // trigger kTRD CINT7 HEE
    cuts.AddCut("00043113", "00200009366300003800000000", "0163103100900000"); // trigger kTRD CINT7 HSE
    cuts.AddCut("00045113", "00200009366300003800000000", "0163103100900000"); // trigger kTRD CINT7 HJE
    cuts.AddCut("00047113", "00200009366300003800000000", "0163103100900000"); // trigger kTRD CINT7 HQU
  } else if (trainConfig == 54) { //pp 2.76TeV cuts
    cuts.AddCut("00052113", "00200009366300003800000000", "0163103100900000"); // trigger kEMC7
    cuts.AddCut("00053113", "00200009366300003800000000", "0163103100900000"); // trigger kEMC8
    cuts.AddCut("00062113", "00200009366300003800000000", "0163103100900000"); // trigger kPHI7
    cuts.AddCut("00063113", "00200009366300003800000000", "0163103100900000"); // trigger kPHI8
  } else if (trainConfig == 55) { //pp 2.76TeV cuts
    cuts.AddCut("00051113", "00200009366300003800000000", "0163103100900000"); // trigger kEMC1
    cuts.AddCut("00071113", "00200009366300003800000000", "0163103100900000"); // trigger kSHM1
    cuts.AddCut("00072113", "00200009366300003800000000", "0163103100900000"); // trigger kSHM7
    cuts.AddCut("00073113", "00200009366300003800000000", "0163103100900000"); // trigger kSHM8
  } else if (trainConfig == 56) { //pp 2.76TeV cuts
    cuts.AddCut("00081113", "00200009366300003800000000", "0163103100900000"); // trigger kEMCEGA + CINT7
    cuts.AddCut("00082113", "00200009366300003800000000", "0163103100900000"); // trigger kEMCEGA + CINT8
    cuts.AddCut("00083113", "00200009366300003800000000", "0163103100900000"); // trigger kEMCEG1 + CINT7
    cuts.AddCut("00084113", "00200009366300003800000000", "0163103100900000"); // trigger kEMCEG1 + CINT8
  } else if (trainConfig == 57) { //pp 2.76TeV cuts
    cuts.AddCut("00085113", "00200009366300003800000000", "0163103100900000"); // trigger kEMCEG2 + CINT7
    cuts.AddCut("00086113", "00200009366300003800000000", "0163103100900000"); // trigger kEMCEG2 + CINT8
    cuts.AddCut("00091113", "00200009366300003800000000", "0163103100900000"); // trigger kEMCEJE + CINT7
    cuts.AddCut("00092113", "00200009366300003800000000", "0163103100900000"); // trigger kEMCEJE + CINT8
  } else if (trainConfig == 58) { //pp 2.76TeV cuts
    cuts.AddCut("00093113", "00200009366300003800000000", "0163103100900000"); // trigger kEMCEJ1 + CINT7
    cuts.AddCut("00094113", "00200009366300003800000000", "0163103100900000"); // trigger kEMCEJ1 + CINT8
    cuts.AddCut("00095113", "00200009366300003800000000", "0163103100900000"); // trigger kEMCEJ2 + CINT7
    cuts.AddCut("00096113", "00200009366300003800000000", "0163103100900000"); // trigger kEMCEJ2 + CINT8		
  } else if (trainConfig == 59) {
    cuts.AddCut("00000113", "00200009257002008250400000", "0152103500000000"); //alpha meson 1
    cuts.AddCut("00000113", "00200009217302008250400000", "0152106500000000"); //pion 0-sigma cut for 0.4GeV<p<3.5GeV above -10-sigma
    cuts.AddCut("00000113", "00200009227302008250400000", "0152106500000000"); //pion 1-sigma cut for 0.4GeV<p<3.5GeV above -10-sigma
    cuts.AddCut("00000113", "00200009287302008250400000", "0152106500000000"); //pion 2-sigma cut for 0.4GeV<p<3.5GeV above   1-sigma
  //---------configs for V0OR 7TeV--------------------------//
  } else if (trainConfig == 60) {
    cuts.AddCut("00000113", "00200009227302008250400000", "0152103500000000"); //New standard cut for eta analysis
    cuts.AddCut("00000113", "00200009227302008250400000", "0152107500000000"); //variation alpha 0.85
    cuts.AddCut("00000113", "00200049227302008250400000", "0152103500000000"); //variation single pt 0.075
    cuts.AddCut("00000113", "00200019227302008250400000", "0152103500000000"); //variation single pt 0.1
  } else if (trainConfig == 61) {
	cuts.AddCut("00000113", "00200006227302008250400000", "0152103500000000"); //variation 70% TPC Cluster
    cuts.AddCut("00000113", "00200008227302008250400000", "0152103500000000"); //variation 35% TPC Cluster
    cuts.AddCut("00000113", "00200009327302008250400000", "0152103500000000"); //variation dE/dx e -4<sigma<5
    cuts.AddCut("00000113", "00200009627302008250400000", "0152103500000000"); //variation dE/dx e -2.5<sigma<4 
  } else if (trainConfig == 62) {
    cuts.AddCut("00000113", "00200009225302008250400000", "0152103500000000"); //variation min p 0.3GeV pion line
    cuts.AddCut("00000113", "00200009220302008250400000", "0152103500000000"); //variation min p 0.5GeV pion line
    cuts.AddCut("00000113", "00200009227102008250400000", "0152103500000000"); //variation max p 5GeV pion line
    cuts.AddCut("00000113", "00200009217302008250400000", "0152103500000000"); //variation nsigma>0 pion line
  } else if (trainConfig == 63) {
    cuts.AddCut("00000113", "00200009257302008250400000", "0152103500000000"); //variation nsigma>2 pion line
    cuts.AddCut("00000113", "00200009227302002250400000", "0152103500000000"); //variation max qt 0.07 1D
    cuts.AddCut("00000113", "00200009227302009250400000", "0152103500000000"); //variation max qt 0.03 2D
    cuts.AddCut("00000113", "00200009227302003250400000", "0152103500000000"); //variation max qt 0.05 1D
  } else if (trainConfig == 64) {
    cuts.AddCut("00000113", "00200009227302008180400000", "0152103500000000"); //variation chi2 50 psi pair 0.2 2D
    cuts.AddCut("00000113", "00200009227302008210400000", "0152103500000000"); //variation chi2 30 psi pair 0.1 1D
    cuts.AddCut("00000113", "00200009227302008860400000", "0152103500000000"); //variation chi2 20 psi pair 0.05 2D
    cuts.AddCut("00000113", "00200009227302008250400000", "0252103500000000"); //variation BG scheme track mult
    //---------configs for V0AND 8TeV--------------------------//
  } else if (trainConfig == 65) {
    cuts.AddCut("00010113", "00200009227302008250400000", "0152103500000000"); //New standard cut for eta analysis
    cuts.AddCut("00010113", "00200009227302008250400000", "0152107500000000"); //variation alpha 0.85
    cuts.AddCut("00010113", "00200049227302008250400000", "0152103500000000"); //variation single pt 0.075
    cuts.AddCut("00010113", "00200019227302008250400000", "0152103500000000"); //variation single pt 0.1
  } else if (trainConfig == 66) {
    cuts.AddCut("00010113", "00200006227302008250400000", "0152103500000000"); //variation 70% TPC Cluster
    cuts.AddCut("00010113", "00200008227302008250400000", "0152103500000000"); //variation 35% TPC Cluster
    cuts.AddCut("00010113", "00200009327302008250400000", "0152103500000000"); //variation dE/dx e -4<sigma<5
    cuts.AddCut("00010113", "00200009627302008250400000", "0152103500000000"); //variation dE/dx e -2.5<sigma<4
  } else if (trainConfig == 67) {
    cuts.AddCut("00010113", "00200009225302008250400000", "0152103500000000"); //variation min p 0.3GeV pion line
    cuts.AddCut("00010113", "00200009220302008250400000", "0152103500000000"); //variation min p 0.5GeV pion line
    cuts.AddCut("00010113", "00200009227102008250400000", "0152103500000000"); //variation max p 5GeV pion line
    cuts.AddCut("00010113", "00200009217302008250400000", "0152103500000000"); //variation nsigma>0 pion line
  } else if (trainConfig == 68) {
    cuts.AddCut("00010113", "00200009257302008250400000", "0152103500000000"); //variation nsigma>2 pion line
    cuts.AddCut("00010113", "00200009227302002250400000", "0152103500000000"); //variation max qt 0.07 1D
    cuts.AddCut("00010113", "00200009227302009250400000", "0152103500000000"); //variation max qt 0.03 2D
    cuts.AddCut("00010113", "00200009227302003250400000", "0152103500000000"); //variation max qt 0.05 1D
  } else if (trainConfig == 69) {
    cuts.AddCut("00010113", "00200009227302008180400000", "0152103500000000"); //variation chi2 50 psi pair 0.2 2D
    cuts.AddCut("00010113", "00200009227302008210400000", "0152103500000000"); //variation chi2 30 psi pair 0.1 1D
    cuts.AddCut("00010113", "00200009227302008860400000", "0152103500000000"); //variation chi2 20 psi pair 0.05 2D
    cuts.AddCut("00010113", "00200009227302008250400000", "0252103500000000"); //variation BG scheme track mult
  } else if (trainConfig == 70) {
    cuts.AddCut("00000113", "00200009227302008250400000", "0152103500000000"); //New standard cut for 7TeV analysis V0OR
    cuts.AddCut("00010113", "00200009227302008250400000", "0152103500000000"); //New standard cut for 8TeV analysis V0AND
    cuts.AddCut("00000113", "00200008386302001200400000", "0152103500000000"); //Old standard cut for 7TeV analysis V0OR
  } else if (trainConfig == 71) {
    cuts.AddCut("00000113", "00200009227302008250400000", "0152103500000000"); //New standard cut for eta analysis
    cuts.AddCut("00000113", "00200009227302008250400000", "0152101500000000"); //variation alpha pT dependent
    cuts.AddCut("00000113", "00200009227302008250400000", "0152109500000000"); //variation alpha
    cuts.AddCut("00000113", "00200009227302008250400000", "0152101500000002"); //variation alpha/opan Max
  } else if (trainConfig == 72) {//added signals
    cuts.AddCut("00000123", "00200009227302008250400000", "0152103500000000"); //New standard cut for eta analysis
    cuts.AddCut("00000123", "00200009227302008250400000", "0152101500000000"); //variation alpha 0.85
    cuts.AddCut("00000123", "00200009227302008250400000", "0152109500000000"); //variation alpha
    cuts.AddCut("00000123", "00200009227302008250400000", "0152101500000002"); //variation alpha opan max				
  } else if (trainConfig == 73) {
    cuts.AddCut("00000113", "00200009227302008250400000", "0152103500000000"); //New standard cut for eta analysis
    cuts.AddCut("00000113", "00200009227302008250404000", "0152101500000000"); //variation alpha pT dependent
    cuts.AddCut("00000113", "00200009227302008250404000", "0152109500000000"); //variation alpha
    cuts.AddCut("00000113", "00200009227302008250400000", "0152101500000000"); //double counting
  } else if (trainConfig == 74) {//added signals
    cuts.AddCut("00000123", "00200009227302008250400000", "0152103500000000"); //New standard cut for eta analysis
    cuts.AddCut("00000123", "00200009227302008250404000", "0152101500000000"); //variation alpha 0.85
    cuts.AddCut("00000123", "00200009227302008250404000", "0152109500000000"); //variation alpha
    cuts.AddCut("00000123", "00200009227302008250400000", "0152101500000000"); //double counting	
  } else if (trainConfig == 75) { //pp 8TeV cuts triggers
    cuts.AddCut("00010113", "00200009227302008250400000", "0152103500000000"); //New standard cut for eta/pi0 analysis
    cuts.AddCut("00052013", "00200009227302008250400000", "0152103500000000"); // trigger kEMC7
    cuts.AddCut("00062013", "00200009227302008250400000", "0152103500000000"); // trigger kPHI7
    cuts.AddCut("00081013", "00200009227302008250400000", "0152103500000000"); // trigger kEMCEGA + CINT7
  } else if (trainConfig == 76) { //pp 8TeV cuts with kappa
    cuts.AddCut("00000113", "00200009227302008250400000", "0152103500000000"); //New standard cut for eta/pi0 analysis kappa 0-2.5
    cuts.AddCut("00000113", "00200009327302008250400000", "0152103500000000"); // kappa 0-2.0
    cuts.AddCut("00000113", "00200009427302008250400000", "0152103500000000"); // kappa 0-1.5
    cuts.AddCut("00000113", "00200009527302008250400000", "0152103500000000"); // kappa 0-1.0
  } else if (trainConfig == 77) { //pp 8TeV cuts with pT dependend alpha cut
    cuts.AddCut("00000113", "00200009227302008250400000", "0152101500000000"); // 0.65, 1.8
    cuts.AddCut("00000113", "00200009227302008250400000", "0152102500000000"); // 0.80, 1.2
    cuts.AddCut("00000113", "00200009227302008250400000", "0152109500000000"); // 0.65, 1.2
  } else if (trainConfig == 78) { //pp 8TeV cuts with EMC triggers
    cuts.AddCut("00010113", "00200009227302008250400000", "0152103500000000","1111111063032230000"); // standard cut 8tev
    cuts.AddCut("00052113", "00200009227302008250400000", "0152103500000000","1111111063032230000"); // trigger kEMC7
    cuts.AddCut("00081113", "00200009227302008250400000", "0152103500000000","1111111063032230000"); // trigger kEGA7
    //---------systematic studies July 2015--------------------------//
  } else if (trainConfig == 80) {
    cuts.AddCut("00000113", "00200009227302008250404000", "0152101500000000"); //New standard cut for eta: alpha pT dependent
    cuts.AddCut("00000113", "00200009227302008250404000", "0152106500000000"); // alpha variation 0.8		
    cuts.AddCut("00000113", "00200009227302008250404000", "0152103500000000"); // alpha variation  1
    cuts.AddCut("00000113", "00200049227302008250404000", "0152101500000000"); // single pT cut 0.075
  } else if (trainConfig == 81) {
    cuts.AddCut("00000123", "00200009227302008250404000", "0152101500000000"); //New standard cut for eta: alpha pT dependent
    cuts.AddCut("00000123", "00200009227302008250404000", "0152106500000000"); // alpha variation  0.8
    cuts.AddCut("00000123", "00200009227302008250404000", "0152103500000000"); // alpha variation  1
    cuts.AddCut("00000123", "00200049227302008250404000", "0152101500000000"); // single pT cut 0.075
  } else if (trainConfig == 82) {
    cuts.AddCut("00000113", "00200019227302008250404000", "0152101500000000"); // single pT cut 0.1
    cuts.AddCut("00000113", "00200008227302008250404000", "0152101500000000"); // TPC cls 0.35	
    cuts.AddCut("00000113", "00200006227302008250404000", "0152101500000000"); // TPC cls 0.7
    cuts.AddCut("00000113", "00200009227302009250404000", "0152101500000000"); // qT cut 0.03 2D
  } else if (trainConfig == 83) {
    cuts.AddCut("00000123", "00200019227302008250404000", "0152101500000000"); // single pT cut 0.1
    cuts.AddCut("00000123", "00200008227302008250404000", "0152101500000000"); // TPC cls 0.35
    cuts.AddCut("00000123", "00200006227302008250404000", "0152101500000000"); // TPC cls 0.7
    cuts.AddCut("00000123", "00200009227302009250404000", "0152101500000000"); // qT cut 0.03 2D
  } else if (trainConfig == 84) {
    cuts.AddCut("00000113", "00200009227302008210404000", "0152101500000000"); // variation chi2 30 psi pair 0.1 1D
    cuts.AddCut("00000113", "00200009227302008180404000", "0152101500000000"); // variation chi2 50 psi pair 0.2 2D	
    cuts.AddCut("00000113", "00200009227302008860404000", "0152101500000000"); // variation chi2 20 psi pair 0.05 2D
    cuts.AddCut("00000113", "00200009227302009250404000", "0252101500000000"); // variation BG scheme track mult
  } else if (trainConfig == 85) {
    cuts.AddCut("00000123", "00200009227302008210404000", "0152101500000000"); // variation chi2 30 psi pair 0.1 1D
    cuts.AddCut("00000123", "00200008227302008180404000", "0152101500000000"); // variation chi2 50 psi pair 0.2 2D
    cuts.AddCut("00000123", "00200006227302008860404000", "0152101500000000"); // variation chi2 20 psi pair 0.05 2D
    cuts.AddCut("00000123", "00200009227302009250404000", "0252101500000000"); // variation BG scheme track mult
  } else if (trainConfig == 86) {
    cuts.AddCut("00000113", "00200009227302008250404000", "0152101500000000"); //New standard cut for eta: alpha pT dependent
    cuts.AddCut("00000113", "00200009227302002250404000", "0152101500000000"); // qT		
    cuts.AddCut("00000113", "00200009227302009250404000", "0152101500000000"); // qT
    cuts.AddCut("00000113", "00200009227302008250004000", "0152101500000000"); // cosPA
  } else if (trainConfig == 87) {
    cuts.AddCut("00000123", "00200009227302008250404000", "0152101500000000"); //New standard cut for eta: alpha pT dependent
    cuts.AddCut("00000123", "00200009227302002250404000", "0152101500000000"); // qT
    cuts.AddCut("00000123", "00200009227302009250404000", "0152101500000000"); // qT
    cuts.AddCut("00000123", "00200009227302008250004000", "0152101500000000"); // cosPA
    //---------systematic studies August 2015 with dEdx--------------------------//
  } else if (trainConfig == 90) {
    cuts.AddCut("00000113", "00200009217302008250404000", "0152101500000000"); //New standard cut for eta: alpha pT dependent
    cuts.AddCut("00000113", "00200009217302008250404000", "0152106500000000"); // alpha variation 0.8		
    cuts.AddCut("00000113", "00200009217302008250404000", "0152103500000000"); // alpha variation  1
    cuts.AddCut("00000113", "00200049217302008250404000", "0152101500000000"); // single pT cut 0.075
  } else if (trainConfig == 91) {
    cuts.AddCut("00000123", "00200009217302008250404000", "0152101500000000"); //New standard cut for eta: alpha pT dependent
    cuts.AddCut("00000123", "00200009217302008250404000", "0152106500000000"); // alpha variation  0.8
    cuts.AddCut("00000123", "00200009217302008250404000", "0152103500000000"); // alpha variation  1
    cuts.AddCut("00000123", "00200049217302008250404000", "0152101500000000"); // single pT cut 0.075
  } else if (trainConfig == 92) {
    cuts.AddCut("00000113", "00200019217302008250404000", "0152101500000000"); // single pT cut 0.1
    cuts.AddCut("00000113", "00200008217302008250404000", "0152101500000000"); // TPC cls 0.35	
    cuts.AddCut("00000113", "00200006217302008250404000", "0152101500000000"); // TPC cls 0.7
    cuts.AddCut("00000113", "00200009217302009250404000", "0152101500000000"); // qT cut 0.03 2D
  } else if (trainConfig == 93) {
    cuts.AddCut("00000123", "00200019217302008250404000", "0152101500000000"); // single pT cut 0.1
    cuts.AddCut("00000123", "00200008217302008250404000", "0152101500000000"); // TPC cls 0.35
    cuts.AddCut("00000123", "00200006217302008250404000", "0152101500000000"); // TPC cls 0.7
    cuts.AddCut("00000123", "00200009217302009250404000", "0152101500000000"); // qT cut 0.03 2D
  } else if (trainConfig == 94) {
    cuts.AddCut("00000113", "00200009217302008210404000", "0152101500000000"); // variation chi2 30 psi pair 0.1 1D
    cuts.AddCut("00000113", "00200009217302008180404000", "0152101500000000"); // variation chi2 50 psi pair 0.2 2D	
    cuts.AddCut("00000113", "00200009217302008860404000", "0152101500000000"); // variation chi2 20 psi pair 0.05 2D
    cuts.AddCut("00000113", "00200009217302009250404000", "0252101500000000"); // variation BG scheme track mult
  } else if (trainConfig == 95) {
    cuts.AddCut("00000123", "00200009217302008210404000", "0152101500000000"); // variation chi2 30 psi pair 0.1 1D
    cuts.AddCut("00000123", "00200009217302008180404000", "0152101500000000"); // variation chi2 50 psi pair 0.2 2D
    cuts.AddCut("00000123", "00200009217302008860404000", "0152101500000000"); // variation chi2 20 psi pair 0.05 2D
    cuts.AddCut("00000123", "00200009217302009250404000", "0252101500000000"); // variation BG scheme track mult
  } else if (trainConfig == 96) {
    cuts.AddCut("00000113", "00200009227302008250404000", "0152101500000000"); //New standard cut for eta: alpha pT dependent
    cuts.AddCut("00000113", "00200009217302002250404000", "0152101500000000"); // qT		
    cuts.AddCut("00000113", "00200009217302009250404000", "0152101500000000"); // qT
    cuts.AddCut("00000113", "00200009217302008250004000", "0152101500000000"); // cosPA
  } else if (trainConfig == 97) {
    cuts.AddCut("00000123", "00200009227302008250404000", "0152101500000000"); //New standard cut for eta: alpha pT dependent
    cuts.AddCut("00000123", "00200009217302002250404000", "0152101500000000"); // qT
    cuts.AddCut("00000123", "00200009217302009250404000", "0152101500000000"); // qT
    cuts.AddCut("00000123", "00200009217302008250004000", "0152101500000000"); // cosPA
  } else if (trainConfig == 98) {
    cuts.AddCut("00000113", "00200009317302008250404000", "0152101500000000"); //dEdx variation
    cuts.AddCut("00000113", "00200009617302008250404000", "0152101500000000"); //dEdx variation
    cuts.AddCut("00000113", "00200009215302008250404000", "0152101500000000"); //dEdx variation
    cuts.AddCut("00000113", "00200009210302008250404000", "0152101500000000"); //dEdx variation
  } else if (trainConfig == 99) {
    cuts.AddCut("00000123", "00200009317302008250404000", "0152101500000000"); //dEdx variation
    cuts.AddCut("00000123", "00200009617302008250404000", "0152101500000000"); //dEdx variation
    cuts.AddCut("00000123", "00200009215302008250404000", "0152101500000000"); //dEdx variation
    cuts.AddCut("00000123", "00200009210302008250404000", "0152101500000000"); //dEdx variation
  } else if (trainConfig == 100) { // kMB
    cuts.AddCut("00100113", "00200009227302008250400000", "0152103500000000"); // 0 -2
    cuts.AddCut("01200113", "00200009227302008250400000", "0152103500000000"); // 2 -5
    cuts.AddCut("02300113", "00200009227302008250400000", "0152103500000000"); // 5 -10
    cuts.AddCut("03400113", "00200009227302008250400000", "0152103500000000"); // 10-30
    cuts.AddCut("04500113", "00200009227302008250400000", "0152103500000000"); // 30-100
  } else if (trainConfig == 101) {  // like trainConfig 71, just with special trigger kINT7
    cuts.AddCut("00010113", "00200009227302008250400000", "0152103500000000"); //New standard cut for eta analysis
    cuts.AddCut("00010113", "00200009227302008250400000", "0152101500000000"); //variation alpha pT dependent
    cuts.AddCut("00010113", "00200009227302008250400000", "0152109500000000"); //variation alpha
    cuts.AddCut("00010113", "00200009227302008250400000", "0152101500000002"); //variation alpha opan max
  } else if (trainConfig == 102) {  // like trainConfig 71, except with smearing added to meson cut
    cuts.AddCut("00010113", "00200009227302008250400000", "0152103500900000"); //New standard cut for eta analysis
    cuts.AddCut("00010113", "00200009227302008250400000", "0152101500900000"); //variation alpha pT dependent
    cuts.AddCut("00010113", "00200009227302008250400000", "0152109500900000"); //variation alpha
    cuts.AddCut("00010113", "00200009227302008250400000", "0152101500900002"); //variation alpha opan max
  } else if (trainConfig == 103) {  // for V0 High-Mult trigger by HM
    cuts.AddCut("00074113", "00200009227302008250400000", "0152103500000000"); //standard cut
    cuts.AddCut("00074113", "00200009227302008250400000", "0152101500000000"); //variation alpha pT dependent
    cuts.AddCut("00074113", "00200009227302008250400000", "0152109500000000"); //variation alpha
    cuts.AddCut("00074113", "00200009227302008250400000", "0152101500000002"); //variation alpha opan max
  } else if (trainConfig == 104) {  // for SPD High-Mult trigger by HM
    cuts.AddCut("00075113", "00200009227302008250400000", "0152103500000000"); //standard cut
    cuts.AddCut("00075113", "00200009227302008250400000", "0152101500000000"); //variation alpha pT dependent
    cuts.AddCut("00075113", "00200009227302008250400000", "0152109500000000"); //variation alpha
    cuts.AddCut("00075113", "00200009227302008250400000", "0152101500000002"); //variation alpha opan max
  } else if (trainConfig == 105) {  // check # of entries w/ pileup rejection cut for V0HM
    cuts.AddCut("00074013", "00200009227302008250400000", "0152103500000000"); //standard cut
    cuts.AddCut("00074013", "00200009227302008250400000", "0152101500000000"); //variation alpha pT dependent
    cuts.AddCut("00074013", "00200009227302008250400000", "0152109500000000"); //variation alpha
    cuts.AddCut("00074013", "00200009227302008250400000", "0152101500000002"); //variation alpha opan max
  } else if (trainConfig == 106) {  // check # of entries w/o pileup rejection cut for SPHM
    cuts.AddCut("00075013", "00200009227302008250400000", "0152103500000000"); //standard cut
    cuts.AddCut("00075013", "00200009227302008250400000", "0152101500000000"); //variation alpha pT dependent
    cuts.AddCut("00075013", "00200009227302008250400000", "0152109500000000"); //variation alpha
    cuts.AddCut("00075013", "00200009227302008250400000", "0152101500000002"); //variation alpha opan max
  } else if (trainConfig == 107) { //LHC11a
    cuts.AddCut("00003113", "00200009227302008250400000", "0152103500000000","1111121050032220000"); //WSDD test with triggers INT1
    cuts.AddCut("00051113", "00200009227302008250400000", "0152103500000000","1111121050032220000"); //WSDD test with triggers EMC1
    cuts.AddCut("00051113", "00202209227302008250400000", "0152103500000000","1111121050032220000"); //WSDD test with triggers EMC1 restricted wide EMC range
    cuts.AddCut("00051113", "00204409227302008250400000", "0152103500000000","1111121050032220000"); //WSDD test with triggers EMC1 restricted tight EMC range
  } else if (trainConfig == 108) { //LHC13g full acceptance
    cuts.AddCut("00010113", "00200009227302008250400000", "0152103500000000","1111121060032220000"); //INT7
    cuts.AddCut("00052113", "00200009227302008250400000", "0152103500000000","1111121060032220000"); //EMC7
    cuts.AddCut("00085113", "00200009227302008250400000", "0152103500000000","1111121060032220000"); //EG2
    cuts.AddCut("00083113", "00200009227302008250400000", "0152103500000000","1111121060032220000"); //EG1
  } else if (trainConfig == 109) { //LHC13g loose EMC acceptance for PCM photons
    cuts.AddCut("00010113", "00202209227302008250400000", "0152103500000000","1111121060032220000"); //INT7
    cuts.AddCut("00052113", "00202209227302008250400000", "0152103500000000","1111121060032220000"); //EMC7
    cuts.AddCut("00085113", "00202209227302008250400000", "0152103500000000","1111121060032220000"); //EG2
    cuts.AddCut("00083113", "00202209227302008250400000", "0152103500000000","1111121060032220000"); //EG1
  } else if (trainConfig == 110) { //LHC13g tight EMC acceptance for PCM photons
    cuts.AddCut("00010113", "00204409227302008250400000", "0152103500000000","1111121060032220000"); //INT7
    cuts.AddCut("00052113", "00204409227302008250400000", "0152103500000000","1111121060032220000"); //EMC7
    cuts.AddCut("00085113", "00204409227302008250400000", "0152103500000000","1111121060032220000"); //EG2
    cuts.AddCut("00083113", "00204409227302008250400000", "0152103500000000","1111121060032220000"); //EG1
  } else if (trainConfig == 111) { //LHC11a
    cuts.AddCut("00003113", "00200009227302008250400000", "0152103500000000","1111121050032220000"); //WSDD test with triggers INT1
    cuts.AddCut("00051013", "00200009227302008250400000", "0152103500000000","1111121050032220000"); //WSDD test with triggers EMC1
    cuts.AddCut("00051013", "00202209227302008250400000", "0152103500000000","1111121050032220000"); //WSDD test with triggers EMC1 restricted wide EMC range
    cuts.AddCut("00051013", "00204409227302008250400000", "0152103500000000","1111121050032220000"); //WSDD test with triggers EMC1 restricted tight EMC range
  } else if (trainConfig == 112) { //LHC13g full acceptance
    cuts.AddCut("00010113", "00200009227302008250400000", "0152103500000000","1111121060032220000"); //INT7
    cuts.AddCut("00052013", "00200009227302008250400000", "0152103500000000","1111121060032220000"); //EMC7
    cuts.AddCut("00085013", "00200009227302008250400000", "0152103500000000","1111121060032220000"); //EG2
    cuts.AddCut("00083013", "00200009227302008250400000", "0152103500000000","1111121060032220000"); //EG1
  } else if (trainConfig == 113) { //LHC13g loose EMC acceptance for PCM photons
    cuts.AddCut("00010113", "00202209227302008250400000", "0152103500000000","1111121060032220000"); //INT7
    cuts.AddCut("00052013", "00202209227302008250400000", "0152103500000000","1111121060032220000"); //EMC7
    cuts.AddCut("00085013", "00202209227302008250400000", "0152103500000000","1111121060032220000"); //EG2
    cuts.AddCut("00083013", "00202209227302008250400000", "0152103500000000","1111121060032220000"); //EG1
  } else if (trainConfig == 114) { //LHC13g tight EMC acceptance for PCM photons
    cuts.AddCut("00010113", "00204409227302008250400000", "0152103500000000","1111121060032220000"); //INT7
    cuts.AddCut("00052013", "00204409227302008250400000", "0152103500000000","1111121060032220000"); //EMC7
    cuts.AddCut("00085013", "00204409227302008250400000", "0152103500000000","1111121060032220000"); //EG2
    cuts.AddCut("00083013", "00204409227302008250400000", "0152103500000000","1111121060032220000"); //EG1
  } else if (trainConfig == 115) { //kINT7
    cuts.AddCut("00010113", "00200009227302008250400000", "0152101500000000"); 
    cuts.AddCut("00010113", "00200009227302008250404000", "0152101500000000"); // double counting cut
    cuts.AddCut("00010113", "00200009327302008250400000", "0152101500000000"); // dEdx 4 sigma below e
    cuts.AddCut("00010113", "00200009327302008250404000", "0152101500000000");
  }  else if (trainConfig == 116) {  // like trainconfig 101, with changes in conversion cuts
    cuts.AddCut("00010113", "00200009227302008250400000", "0152101500000000");
    cuts.AddCut("00010113", "00200009227302008250404000", "0152101500000000"); //double counting
    cuts.AddCut("00010113", "00200009327302008250400000", "0152101500000000"); //variation dE/dx e -2.5<sigma<4
    cuts.AddCut("00010113", "00200009227312008250400000", "0152101500000000"); //variation nSigma rejection cut
  } else if (trainConfig == 117) { //kINT7
    cuts.AddCut("00010113", "00200009227302008250404000", "0152101500000000");
    cuts.AddCut("00010113", "00200009327302008250404000", "0152101500000000"); // dEdx 4 sigma below e
    cuts.AddCut("00010113", "00200079227302008250404000", "0152101500000000"); // pT cut at 0
    cuts.AddCut("00010113", "00200079327302008250404000", "0152101500000000");
  } else {
    Error(Form("GammaConvV1_%i",trainConfig), "wrong trainConfig variable no cuts have been specified for the configuration");
    return;
  }
  
	if(!cuts.AreValid()){
    cout << "\n\n****************************************************" << endl;
    cout << "ERROR: No valid cuts stored in CutHandlerConv! Returning..." << endl;
    cout << "****************************************************\n\n" << endl;
    return;
  }

  Int_t numberOfCuts = cuts.GetNCuts(); 
   
  TList *EventCutList = new TList();
  TList *ConvCutList = new TList();
  TList *MesonCutList = new TList();
  TList *ClusterCutList = new TList();

  TList *HeaderList = new TList();
  if (periodname.Contains("LHC12i3")){	
    TObjString *Header2 = new TObjString("BOX");
    HeaderList->Add(Header2);
  } else if (periodname.CompareTo("LHC14e2b")==0){
    TObjString *Header2 = new TObjString("pi0_1");
    HeaderList->Add(Header2);
    TObjString *Header3 = new TObjString("eta_2");
    HeaderList->Add(Header3);
  }	
  
  TString energy = "";
  TString mcName = "";
  TString mcNameAdd = "";
  if (periodname.Contains("WOSDD")){
    mcNameAdd = "_WOSDD";
  } else if (periodname.Contains("WSDD")){
    mcNameAdd = "_WSDD";
  } 	
  if (periodname.Contains("LHC12i3")){
    energy = "2760GeV";
    mcName = "Pythia8_LHC12i3";
  } else if (periodname.Contains("LHC12f1a")){	
    energy = "2760GeV";
    mcName = "Pythia8_LHC12f1a";	
  } else if (periodname.Contains("LHC12f1b")){	
    energy = "2760GeV";
    mcName = "Phojet_LHC12f1b";			
  } else if (periodname.Contains("LHC14e2a")){	
    energy = "8TeV";
    mcName = "Pythia8_LHC14e2a";			
  } else if (periodname.Contains("LHC14e2b")){	
    energy = "8TeV";
    mcName = "Pythia8_LHC14e2b";				
  } else if (periodname.Contains("LHC14e2c")){		
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
    analysisEventCuts[i]          = new AliConvEventCuts();
    TString fitNamePi0            = Form("Pi0_Fit_Data_%s",energy.Data());
    TString fitNameEta            = Form("Eta_Fit_Data_%s",energy.Data());
    TString fAddedSignalString    = (cuts.GetEventCut(i)).Data();
    fAddedSignalString            = fAddedSignalString(6,1);
    Bool_t fAddedSignal           = kFALSE;
    if (fAddedSignalString.CompareTo("2") == 0) 
      fAddedSignal                = kTRUE;

    TString mcInputNamePi0        = "";
    TString mcInputNameEta        = "";
    if (fAddedSignal && (periodname.Contains("LHC12i3") || periodname.CompareTo("LHC14e2b")==0)){
      mcInputNamePi0              = Form("Pi0_%s%s_addSig_%s", mcName.Data(), mcNameAdd.Data(), energy.Data() );
      mcInputNameEta              = Form("Eta_%s%s_addSig_%s", mcName.Data(), mcNameAdd.Data(), energy.Data() );
    } else {
      mcInputNamePi0              = Form("Pi0_%s%s_%s", mcName.Data(), mcNameAdd.Data(), energy.Data() );
      mcInputNameEta              = Form("Eta_%s%s_%s", mcName.Data(), mcNameAdd.Data(), energy.Data() );
    }	
    //		if (doParticleWeighting) analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kFALSE, kFALSE, kFALSE, fileNameInputForPartWeighting, mcInputNamePi0, mcInputNameEta, "",fitNamePi0,fitNameEta);
    if (doParticleWeighting) analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE, fileNameInputForPartWeighting, mcInputNamePi0, mcInputNameEta, "",fitNamePi0,fitNameEta);

    TString dataInputMultHisto    = "";
    TString mcInputMultHisto      = "";
    TString triggerString         = (cuts.GetEventCut(i)).Data();
    triggerString                 = triggerString(3,2);
    if (triggerString.CompareTo("03")==0) 
      triggerString               = "00";
    if (periodNameAnchor.CompareTo("LHC13g") == 0 && triggerString.CompareTo("10")== 0 )
      triggerString               = "00";

    dataInputMultHisto            = Form("%s_%s", periodNameAnchor.Data(), triggerString.Data());
    mcInputMultHisto              = Form("%s_%s", periodNameV0Reader.Data(), triggerString.Data());
   
    if (doMultiplicityWeighting){
      cout << "enableling mult weighting" << endl;
      analysisEventCuts[i]->SetUseWeightMultiplicityFromFile( kTRUE, fileNameInputForMultWeighing, dataInputMultHisto, mcInputMultHisto );
    }

    
    analysisEventCuts[i]->SetTriggerMimicking(enableTriggerMimicking);
    analysisEventCuts[i]->SetTriggerOverlapRejecion(enableTriggerOverlapRej);
    analysisEventCuts[i]->SetMaxFacPtHard(maxFacPtHard);
    analysisEventCuts[i]->SetV0ReaderName(V0ReaderName);
    analysisEventCuts[i]->InitializeCutsFromCutString((cuts.GetEventCut(i)).Data());

    EventCutList->Add(analysisEventCuts[i]);
    analysisEventCuts[i]->SetFillCutHistograms("",kFALSE);

    if (trainConfig > 106 && trainConfig < 115 || trainConfig == 78){
        enableClustersForTrigger  = kTRUE;
        analysisClusterCuts[i]    = new AliCaloPhotonCuts();
        analysisClusterCuts[i]->SetV0ReaderName(V0ReaderName);
        analysisClusterCuts[i]->InitializeCutsFromCutString((cuts.GetClusterCut(i)).Data());
        ClusterCutList->Add(analysisClusterCuts[i]);
        analysisClusterCuts[i]->SetFillCutHistograms("");   
    }  
    
    analysisCuts[i]               = new AliConversionPhotonCuts();
    if (trainConfig == 76 ){
      analysisCuts[i]->SetSwitchToKappaInsteadOfNSigdEdxTPC(kTRUE);
    }
    analysisCuts[i]->SetV0ReaderName(V0ReaderName);
    analysisCuts[i]->InitializeCutsFromCutString((cuts.GetPhotonCut(i)).Data());
    
    if( trainConfig  == 73 || trainConfig  == 74 || (trainConfig  >= 80 && trainConfig  <= 87) ){
      analysisCuts[i]->SetDodEdxSigmaCut(kFALSE);
    }
    ConvCutList->Add(analysisCuts[i]);
    analysisCuts[i]->SetFillCutHistograms("",kFALSE);

    analysisMesonCuts[i]          = new AliConversionMesonCuts();
    analysisMesonCuts[i]->InitializeCutsFromCutString((cuts.GetMesonCut(i)).Data());
    MesonCutList->Add(analysisMesonCuts[i]);
    analysisMesonCuts[i]->SetFillCutHistograms("");		
    analysisEventCuts[i]->SetAcceptedHeader(HeaderList);
    
  
  }

  task->SetEventCutList(numberOfCuts,EventCutList);
  task->SetConversionCutList(numberOfCuts,ConvCutList);
  task->SetMesonCutList(numberOfCuts,MesonCutList);
  task->SetMoveParticleAccordingToVertex(kTRUE);
  task->SetDoMesonAnalysis(kTRUE);
  task->SetDoMesonQA(enableQAMesonTask); //Attention new switch for Pi0 QA
  task->SetDoPhotonQA(enableQAPhotonTask);  //Attention new switch small for Photon QA
  if (enableClustersForTrigger){
    task->SetDoClusterSelectionForTriggerNorm(enableClustersForTrigger);
    task->SetClusterCutList(numberOfCuts,ClusterCutList);
  }  
  
  //connect containers
  AliAnalysisDataContainer *coutput =
    mgr->CreateContainer(Form("GammaConvV1_%i",trainConfig), TList::Class(),
              AliAnalysisManager::kOutputContainer,Form("GammaConvV1_%i.root",trainConfig));

  mgr->AddTask(task);
  mgr->ConnectInput(task,0,cinput);
  mgr->ConnectOutput(task,1,coutput);

  return;

}
