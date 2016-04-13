/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: Friederike Bock, Annika Passfeld                               *
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
      eventCutArray = new TString[nMaxCuts]; photonCutArray = new TString[nMaxCuts]; mesonCutArray = new TString[nMaxCuts];
      for(Int_t i=0; i<nMaxCuts; i++) {eventCutArray[i] = ""; photonCutArray[i] = ""; mesonCutArray[i] = "";}
    }

    void AddCut(TString eventCut, TString photonCut, TString mesonCut){
      if(nCuts>=nMaxCuts) {cout << "ERROR in CutHandlerConv: Exceeded maximum number of cuts!" << endl; validCuts = false; return;}
      if( eventCut.Length()!=8 || photonCut.Length()!=26 || mesonCut.Length()!=16 ) {cout << "ERROR in CutHandlerConv: Incorrect length of cut string!" << endl; validCuts = false; return;}
      eventCutArray[nCuts]=eventCut; photonCutArray[nCuts]=photonCut; mesonCutArray[nCuts]=mesonCut;
      nCuts++;
      return;
    }
    Bool_t AreValid(){return validCuts;}
    Int_t GetNCuts(){if(validCuts) return nCuts; else return 0;}
    TString GetEventCut(Int_t i){if(validCuts&&i<nMaxCuts&&i>=0) return eventCutArray[i]; else{cout << "ERROR in CutHandlerConv: GetEventCut wrong index i" << endl;return "";}}
    TString GetPhotonCut(Int_t i){if(validCuts&&i<nMaxCuts&&i>=0) return photonCutArray[i]; else {cout << "ERROR in CutHandlerConv: GetPhotonCut wrong index i" << endl;return "";}}
    TString GetMesonCut(Int_t i){if(validCuts&&i<nMaxCuts&&i>=0) return mesonCutArray[i]; else {cout << "ERROR in CutHandlerConv: GetMesonCut wrong index i" << endl;return "";}}
  private:
    Bool_t validCuts;
    Int_t nCuts; Int_t nMaxCuts;
    TString* eventCutArray;
    TString* photonCutArray;
    TString* mesonCutArray;
};

void AddTask_GammaConvV1_pPb3(  Int_t         trainConfig                   = 1,                                // change different set of cuts
                                Int_t         isMC                          = 0,                                // run MC
                                Int_t         enableQAMesonTask             = 0,                                // enable QA in AliAnalysisTaskGammaConvV1
                                Int_t         enableQAPhotonTask            = 0,                                // enable additional QA task
                                TString       fileNameInputForWeighting     = "MCSpectraInput.root",            // path to file for weigting input
                                Bool_t        doWeightingPart               = kFALSE,                           // enable Weighting
                                TString       generatorName                 = "DPMJET",                         // generator Name    
                                TString       cutnumberAODBranch            = "800000006008400000001500000",    // cutnumber for AOD branch
                                Bool_t        enableV0findingEffi           = kFALSE,                           // enables V0finding efficiency histograms
                                Bool_t        enablePlotVsCentrality        = kFALSE,
                                Bool_t        enableTriggerMimicking        = kFALSE,                           // enable trigger mimicking
                                Bool_t        enableTriggerOverlapRej       = kFALSE,                           // enable trigger overlap rejection
                                Float_t       maxFacPtHard                  = 3.,                               // maximum factor between hardest jet and ptHard generated
                                TString       periodNameV0Reader            = "",
                                Bool_t        doMultiplicityWeighting       = kFALSE,                           //
                                TString       fileNameInputForMultWeighing  = "Multiplicity.root",              //
                                TString       periodNameAnchor              = ""
                          ) {

  Int_t isHeavyIon = 2;
  
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
  //offline V0Finder
  TString cutnumberPhoton = "16000008400100001500000000";
  TString cutnumberEvent = "80000103";
  Bool_t doEtaShift = kFALSE;
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
      if(fCuts->InitializeCutsFromCutString(cutnumberPhoton.Data())){
        fV0ReaderV1->SetConversionCuts(fCuts);
        fCuts->SetFillCutHistograms("",kTRUE);
      }
    }
    if(inputHandler->IsA()==AliAODInputHandler::Class()){
    // AOD mode
      cout << "AOD handler: adding " << cutnumberAODBranch.Data() << " as conversion branch" << endl;
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
  //================================================
  //            find input container
  AliAnalysisTaskGammaConvV1 *task=NULL;
  task= new AliAnalysisTaskGammaConvV1(Form("GammaConvV1_%i",trainConfig));
  task->SetIsHeavyIon(isHeavyIon);
  task->SetIsMC(isMC);
    task->SetV0ReaderName(V0ReaderName);
  // Cut Numbers to use in Analysis
  
  CutHandlerConv cuts;
  
  Bool_t doEtaShiftIndCuts = kFALSE;
  TString stringShift = "";
  
  if (trainConfig == 1) {
    cuts.AddCut("80000113", "10200009217000008260400000", "0162103500900000");  // all Photon Qualities  //offline V0Finder
    cuts.AddCut("80000113", "10200009217000008260420000", "0162103500900000");  // only photon quality 1  //offline V0Finder
    cuts.AddCut("80000113", "10200009217000008260430000", "0162103500900000");  // only photon quality 2  //offline V0Finder
    cuts.AddCut("80000113", "10200009217000008260440000", "0162103500900000");  // only photon quality 3  //offline V0Finder
  } else if (trainConfig == 2) {
    cuts.AddCut("80000123", "10200009217000008260400000", "0162103500900000");  // all Photon Qualities  //offline V0Finder
    cuts.AddCut("80000123", "10200009217000008260420000", "0162103500900000");  // only photon quality 1  //offline V0Finder
    cuts.AddCut("80000123", "10200009217000008260430000", "0162103500900000");  // only photon quality 2  //offline V0Finder
    cuts.AddCut("80000123", "10200009217000008260440000", "0162103500900000");  // only photon quality 3  //offline V0Finder
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
  
  TList *HeaderList = new TList();
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
  
  
  Bool_t doWeighting = kFALSE;
  if (doWeightingPart == 1 || doWeightingPart == 2 || doWeightingPart == 3) doWeighting = kTRUE;
  
  EventCutList->SetOwner(kTRUE);
  AliConvEventCuts **analysisEventCuts = new AliConvEventCuts*[numberOfCuts];
  ConvCutList->SetOwner(kTRUE);
  AliConversionPhotonCuts **analysisCuts = new AliConversionPhotonCuts*[numberOfCuts];
  MesonCutList->SetOwner(kTRUE);
  AliConversionMesonCuts **analysisMesonCuts = new AliConversionMesonCuts*[numberOfCuts];
  
  
  for(Int_t i = 0; i<numberOfCuts; i++){
    
    analysisEventCuts[i] = new AliConvEventCuts();
    if ( trainConfig == 1  ){
      if (doWeighting){
        if (generatorName.CompareTo("DPMJET")==0){
          analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE, fileNameInputForWeighting, "Pi0_DPMJET_LHC13b2_efix_pPb_5023GeV_MBV0A", "Eta_DPMJET_LHC13b2_efix_pPb_5023GeV_MBV0A", "","Pi0_Fit_Data_pPb_5023GeV_MBV0A","Eta_Fit_Data_pPb_5023GeV_MBV0A");
        } else if (generatorName.CompareTo("HIJING")==0){
          analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13e7_pPb_5023GeV_MBV0A", "Eta_Hijing_LHC13e7_pPb_5023GeV_MBV0A", "","Pi0_Fit_Data_pPb_5023GeV_MBV0A","Eta_Fit_Data_pPb_5023GeV_MBV0A");
        }
      }
    }   
    if ( trainConfig == 2 ){
      if (doWeighting){
        analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13e7_addSig_pPb_5023GeV_MBV0A", "Eta_Hijing_LHC13e7_addSig_pPb_5023GeV_MBV0A", "","Pi0_Fit_Data_pPb_5023GeV_MBV0A","Eta_Fit_Data_pPb_5023GeV_MBV0A");
      }
      
    }   
    analysisEventCuts[i]->SetTriggerMimicking(enableTriggerMimicking);
    analysisEventCuts[i]->SetTriggerOverlapRejecion(enableTriggerOverlapRej);
    analysisEventCuts[i]->SetMaxFacPtHard(maxFacPtHard);
    analysisEventCuts[i]->SetV0ReaderName(V0ReaderName);
    analysisEventCuts[i]->InitializeCutsFromCutString((cuts.GetEventCut(i)).Data());
        
    if (doEtaShiftIndCuts) {
      analysisEventCuts[i]->DoEtaShift(doEtaShiftIndCuts);
      analysisEventCuts[i]->SetEtaShift(stringShift);
    }
    if (trainConfig == 179 || trainConfig == 180) {
      analysisEventCuts[0]->DoEtaShift(kTRUE);
      analysisEventCuts[0]->SetEtaShift(-0.215);
    }
    
    EventCutList->Add(analysisEventCuts[i]);
    analysisEventCuts[i]->SetFillCutHistograms("",kFALSE);
    
    analysisCuts[i] = new AliConversionPhotonCuts();
    analysisCuts[i]->SetV0ReaderName(V0ReaderName);
    analysisCuts[i]->InitializeCutsFromCutString((cuts.GetPhotonCut(i)).Data());
    analysisCuts[i]->SetIsHeavyIon(isHeavyIon);
    ConvCutList->Add(analysisCuts[i]);
    analysisCuts[i]->SetFillCutHistograms("",kFALSE);
    
    analysisMesonCuts[i] = new AliConversionMesonCuts();
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
  task->SetDoPlotVsCentrality(enablePlotVsCentrality);
  
  //connect containers
  AliAnalysisDataContainer *coutput =
    mgr->CreateContainer(Form("GammaConvV1_%i",trainConfig), TList::Class(),
              AliAnalysisManager::kOutputContainer,Form("GammaConvV1_%i.root",trainConfig));
  
  mgr->AddTask(task);
  mgr->ConnectInput(task,0,cinput);
  mgr->ConnectOutput(task,1,coutput);
  
  return;
  
}
