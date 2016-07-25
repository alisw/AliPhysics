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
//pp together with all supporting classes
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
void AddTask_OmegaToPiZeroGamma(  Int_t     trainConfig                   = 1,                      // change different set of cuts
                                Int_t     isMC                          = 0,                      // run MC
                                Int_t     enableQAMesonTask             = 1,                      // enable QA in AliAnalysisTaskGammaConvV1
                                Int_t     enableQAPhotonTask            = 1,                      // enable additional QA task
                                TString   cutnumberAODBranch            = "000000006008400001001500000",
                                Int_t     enableExtMatchAndQA           = 0,                      // enable matching histograms (1) and extended QA (2), only QA(3), all disabled (0)
                                Bool_t    enableV0findingEffi           = kFALSE,                 // enables V0finding efficiency histograms
                                Bool_t    enableTriggerMimicking        = kFALSE,                 // enable trigger mimicking
                                Bool_t    enableTriggerOverlapRej       = kFALSE,                 // enable trigger overlap rejection
                                Float_t   maxFacPtHard                  = 3.,                     // maximum factor between hardest jet and ptHard generated
                                TString   periodNameV0Reader            = "",                     // period Name for V0 Reader 
                                Bool_t    doMultiplicityWeighting       = kFALSE,                 // enable multiplicity weights
                                TString   fileNameInputForMultWeighing  = "Multiplicity.root",    // file for multiplicity weights
                                TString   periodNameAnchor              = "",                     // anchor period name for mult weighting
                                Bool_t    enableSortingMCLabels         = kTRUE,                  // enable sorting for MC cluster labels
                                Bool_t    runLightOutput                = kFALSE,                 // switch to run light output (only essential histograms for afterburner)
                                Bool_t    doSmear                       = kFALSE,                 // switches to run user defined smearing
                                Double_t  bremSmear                     = 1.,
                                Double_t  smearPar                      = 0.,                     // conv photon smearing params
                                Double_t  smearParConst                 = 0.                      // conv photon smearing params
              ) {
  
  Int_t isHeavyIon = 0;

  // fReconMethod = first digit of trainconfig
  Int_t ReconMethod = trainConfig/100;
  
  // ================== GetAnalysisManager ===============================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error(Form("AddTask_OmegaToPiZeroGamma_%i_%i", trainConfig, isMC), "No analysis manager found.");
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
  
  Printf("here \n");
  
  //=========  Set Cutnumber for V0Reader ================================
  TString cutnumberPhoton = "00000008400100001500000000";
  TString cutnumberEvent = "00000003";
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
  AliAnalysisTaskOmegaToPiZeroGamma *task=NULL;
  task= new AliAnalysisTaskOmegaToPiZeroGamma(Form("OmegaToPiZeroGamma_%i_%i", trainConfig, isMC));
  task->SetIsHeavyIon(isHeavyIon);
  task->SetIsMC(isMC);
  task->SetV0ReaderName(V0ReaderName);
  task->SetReconMethod(ReconMethod);

  //create cut handler
  CutHandlerPiZeroGamma cuts;

  // cluster cuts
  // 0 "ClusterType",  1 "EtaMin", 2 "EtaMax", 3 "PhiMin", 4 "PhiMax", 5 "DistanceToBadChannel", 6 "Timing", 7 "TrackMatching", 8 "ExoticCell",
  // 9 "MinEnergy", 10 "MinNCells", 11 "MinM02", 12 "MaxM02", 13 "MinM20", 14 "MaxM20", 15 "MaximumDispersion", 16 "NLM"

  // ************************************* EMCAL cuts ****************************************************

  // 7 TeV
  // cuts for ReconMethod==0
  if(trainConfig == 1){ // EMCAL clusters pp 7 TeV
    cuts.AddCut("00000113","00200009327000008250400000","1111111013032230000","0163103000000010","0163103000000010"); // pion mass (0.08,0.145), last string is not used as of now

  } else if(trainConfig == 2){ // std 8TeV
    cuts.AddCut("00010113","00200009327000008250400000","1111111063032230000","0163103100000010","0163103000000010");
    cuts.AddCut("00052113","00200009327000008250400000","1111111063032230000","0163103100000010","0163103000000010");
    cuts.AddCut("00081113","00200009327000008250400000","1111111063032230000","0163103100000010","0163103000000010");

  // cuts for ReconMethod==1
  } else if(trainConfig == 101){ // EMCAL clusters pp 7 TeV
    cuts.AddCut("00000113","00200009327000008250400000","1111111013032230000","0163103000000010","0163103000000010"); // pion mass (0.08,0.145)

  } else if(trainConfig == 102){ // std 8TeV
    cuts.AddCut("00010113","00200009327000008250400000","1111111063032230000","0163103100000010","0163103000000010");
    cuts.AddCut("00052113","00200009327000008250400000","1111111063032230000","0163103100000010","0163103000000010");
    cuts.AddCut("00081113","00200009327000008250400000","1111111063032230000","0163103100000010","0163103000000010");

  // cuts for ReconMethod==2
  } else if(trainConfig == 201){ // EMCAL clusters pp 7 TeV
    cuts.AddCut("00000113","00200009327000008250400000","1111111013032230000","0163103000000010","0163103000000010"); // pion mass (0.08,0.145)

  } else if(trainConfig == 202){ // std 8TeV
    cuts.AddCut("00010113","00200009327000008250400000","1111111063032230000","0163103100000010","0163103000000010");
    cuts.AddCut("00052113","00200009327000008250400000","1111111063032230000","0163103100000010","0163103000000010");
    cuts.AddCut("00081113","00200009327000008250400000","1111111063032230000","0163103100000010","0163103000000010");

  // cuts for ReconMethod==3
  } else if(trainConfig == 301){ // EMCAL clusters pp 7 TeV
    cuts.AddCut("00000113","00200009327000008250400000","1111111013032230000","0163103000000010","0163103000000010"); // pion mass (0.08,0.145)

  } else if(trainConfig == 302){ // std 8TeV
    cuts.AddCut("00010113","00200009327000008250400000","1111111063032230000","0163103100000010","0163103000000010");
    cuts.AddCut("00052113","00200009327000008250400000","1111111063032230000","0163103100000010","0163103000000010");
    cuts.AddCut("00081113","00200009327000008250400000","1111111063032230000","0163103100000010","0163103000000010");

  // cuts for ReconMethod==4
  } else if(trainConfig == 401){ // EMCAL clusters pp 7 TeV
    cuts.AddCut("00000113","00200009327000008250400000","1111111013032230000","0163103000000010","0163103000000010"); // pion mass (0.08,0.145)

  } else if(trainConfig == 402){ // std 8TeV
    cuts.AddCut("00010113","00200009327000008250400000","1111111063032230000","0163103100000010","0163103000000010");
    cuts.AddCut("00052113","00200009327000008250400000","1111111063032230000","0163103100000010","0163103000000010");
    cuts.AddCut("00081113","00200009327000008250400000","1111111063032230000","0163103100000010","0163103000000010");

  // cuts for ReconMethod==5
  } else if(trainConfig == 501){ // EMCAL clusters pp 7 TeV
    cuts.AddCut("00000113","00200009327000008250400000","1111111013032230000","0163103000000010","0163103000000010"); // pion mass (0.08,0.145)

  } else if(trainConfig == 502){ // std 8TeV
    cuts.AddCut("00010113","00200009327000008250400000","1111111063032230000","0163103100000010","0163103000000010");
    cuts.AddCut("00052113","00200009327000008250400000","1111111063032230000","0163103100000010","0163103000000010");
    cuts.AddCut("00081113","00200009327000008250400000","1111111063032230000","0163103100000010","0163103000000010");

  } else {
    Error(Form("OmegaToPiZeroGamma_%i_%i", trainConfig, isMC), "wrong trainConfig variable no cuts have been specified for the configuration");
    return;
  }

  if(!cuts.AreValid()){
    cout << "\n\n****************************************************" << endl;
    cout << "ERROR: No valid cuts stored in CutHandlerPiZeroGamma! Returning..." << endl;
    cout << "****************************************************\n\n" << endl;
    return;
  }

  Int_t numberOfCuts = cuts.GetNCuts();

  TList *EventCutList = new TList();
  TList *ConvCutList = new TList();
  TList *ClusterCutList = new TList();
  TList *neutralPionCutList = new TList();
  TList *MesonCutList = new TList();

  EventCutList->SetOwner(kTRUE);
  AliConvEventCuts **analysisEventCuts = new AliConvEventCuts*[numberOfCuts];
  ConvCutList->SetOwner(kTRUE);
  AliConversionPhotonCuts **analysisCuts = new AliConversionPhotonCuts*[numberOfCuts];
  ClusterCutList->SetOwner(kTRUE);
  AliCaloPhotonCuts **analysisClusterCuts = new AliCaloPhotonCuts*[numberOfCuts];
  neutralPionCutList->SetOwner(kTRUE);
  AliConversionMesonCuts **analysisNeutralPionCuts = new AliConversionMesonCuts*[numberOfCuts];
  MesonCutList->SetOwner(kTRUE);
  AliConversionMesonCuts **analysisMesonCuts = new AliConversionMesonCuts*[numberOfCuts];

  for(Int_t i = 0; i<numberOfCuts; i++){
    analysisEventCuts[i] = new AliConvEventCuts();
    
    TString dataInputMultHisto  = "";
    TString mcInputMultHisto    = "";
    TString triggerString   = cuts.GetEventCut(i);
    triggerString           = triggerString(3,2);
    if (triggerString.CompareTo("03")==0) 
      triggerString         = "00";
    if (periodNameAnchor.CompareTo("LHC13g") == 0 && triggerString.CompareTo("10")== 0 )
      triggerString         = "00";

    dataInputMultHisto      = Form("%s_%s", periodNameAnchor.Data(), triggerString.Data());
    mcInputMultHisto        = Form("%s_%s", periodNameV0Reader.Data(), triggerString.Data());
   
    if (doMultiplicityWeighting){
      cout << "enableling mult weighting" << endl;
      analysisEventCuts[i]->SetUseWeightMultiplicityFromFile( kTRUE, fileNameInputForMultWeighing, dataInputMultHisto, mcInputMultHisto );
    }
              
    
    analysisEventCuts[i]->SetTriggerMimicking(enableTriggerMimicking);
    analysisEventCuts[i]->SetTriggerOverlapRejecion(enableTriggerOverlapRej);
    analysisEventCuts[i]->SetMaxFacPtHard(maxFacPtHard);
    analysisEventCuts[i]->SetV0ReaderName(V0ReaderName);
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
    analysisClusterCuts[i]->SetLightOutput(runLightOutput);
    analysisClusterCuts[i]->InitializeCutsFromCutString((cuts.GetClusterCut(i)).Data());
    ClusterCutList->Add(analysisClusterCuts[i]);
    analysisClusterCuts[i]->SetExtendedMatchAndQA(enableExtMatchAndQA);
    analysisClusterCuts[i]->SetFillCutHistograms("");
    
    analysisNeutralPionCuts[i] = new AliConversionMesonCuts();
    analysisNeutralPionCuts[i]->SetLightOutput(runLightOutput);
    analysisNeutralPionCuts[i]->SetRunningMode(2);
    analysisNeutralPionCuts[i]->InitializeCutsFromCutString((cuts.GetNeutralPionCut(i)).Data());
    neutralPionCutList->Add(analysisNeutralPionCuts[i]);
    analysisNeutralPionCuts[i]->SetFillCutHistograms("");
    if(doSmear) analysisNeutralPionCuts[i]->SetDefaultSmearing(bremSmear,smearPar,smearParConst);

    analysisMesonCuts[i] = new AliConversionMesonCuts();
    analysisMesonCuts[i]->SetLightOutput(runLightOutput);
    analysisMesonCuts[i]->SetRunningMode(2);
    analysisMesonCuts[i]->InitializeCutsFromCutString((cuts.GetMesonCut(i)).Data());
    MesonCutList->Add(analysisMesonCuts[i]);
    analysisMesonCuts[i]->SetFillCutHistograms("");
    if(doSmear) analysisMesonCuts[i]->SetDefaultSmearing(bremSmear,smearPar,smearParConst);
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
  if(enableExtMatchAndQA == 2 || enableExtMatchAndQA == 3){ task->SetPlotHistsExtQA(kTRUE);}

  //connect containers
  AliAnalysisDataContainer *coutput =
    mgr->CreateContainer(Form("OmegaToPiZeroGamma_%i_%i", trainConfig, isMC), TList::Class(),
              AliAnalysisManager::kOutputContainer,Form("OmegaToPiZeroGamma_%i_%i.root", trainConfig, isMC));

  mgr->AddTask(task);
  mgr->ConnectInput(task,0,cinput);
  mgr->ConnectOutput(task,1,coutput);

  return;

}
