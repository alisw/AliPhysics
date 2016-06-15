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

void AddTask_GammaConvV1_PbPb2( Int_t         trainConfig                   = 1,                                // change different set of cuts
                                Int_t         isMC                          = 0,                                // run MC
                                Int_t         enableQAMesonTask             = 0,                                // enable QA in AliAnalysisTaskGammaConvV1
                                Int_t         enableQAPhotonTask            = 0,                                // enable additional QA task
                                TString     fileNameInputForWeighting       = "MCSpectraInput.root",            // path to file for weigting input
                                Bool_t         doWeighting                  = kFALSE,                           // enable Weighting
                                TString     cutnumberAODBranch              = "100000006008400000001500000",    // cutnumber with which AODs have been filtered
                                Bool_t         enableV0findingEffi          = kFALSE,                           // enables V0finding efficiency histograms
                                Bool_t         enableTriggerMimicking       = kFALSE,                           // enable trigger mimicking
                                Bool_t         enableTriggerOverlapRej      = kFALSE,                           // enable trigger overlap rejection
                                Float_t        maxFacPtHard                 = 3.,                               // maximum factor between hardest jet and ptHard generated
                                TString        periodNameV0Reader           = "",
                                Bool_t         runLightOutput               = kFALSE                            // switch to run light output (only essential histograms for afterburner)
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
  TString cutnumberPhoton = "00000008400100001500000000";
  TString cutnumberEvent = "10000003";
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

    AliLog::SetGlobalLogLevel(AliLog::kInfo);

    //connect input V0Reader
    mgr->AddTask(fV0ReaderV1);
    mgr->ConnectInput(fV0ReaderV1,0,cinput);

  }

  //================================================
  //========= Add task to the ANALYSIS manager =====
  //================================================
  AliAnalysisTaskGammaConvV1 *task=NULL;
  task= new AliAnalysisTaskGammaConvV1(Form("GammaConvV1_%i",trainConfig));
  task->SetIsHeavyIon(isHeavyIon);
  task->SetIsMC(isMC);
  task->SetV0ReaderName(V0ReaderName);
  task->SetLightOutput(runLightOutput);
  // Cut Numbers to use in Analysis
  
  CutHandlerConv cuts;

  if (trainConfig == 1){ 
    cuts.AddCut("60100013", "04200009297002003220000000", "0152204500900000"); 
  } else if (trainConfig == 2) { 
    cuts.AddCut("61200013", "04200009297002003220000000", "0152204500900000"); 
  } else if (trainConfig == 3) { 
    cuts.AddCut("50100013", "04200009297002003220000000", "0152204500900000"); 
  } else if (trainConfig == 4) { 
    cuts.AddCut("50200013", "04200009297002003220000000", "0152204500900000");    
  } else if (trainConfig == 5) { 
    cuts.AddCut("51200013", "04200009297002003220000000", "0152204500900000");    
  } else if (trainConfig == 6) { 
    cuts.AddCut("52400013", "04200009297002003220000000", "0152204500900000");       
  } else if (trainConfig == 7) {    
    cuts.AddCut("54600013", "04200009297002003220000000", "0152206500900000"); 
  } else if (trainConfig == 8) {    
    cuts.AddCut("54800013", "04200009297002003220000000", "0152206500900000");    
  } else if (trainConfig == 9) {    
    cuts.AddCut("54500013", "04200009297002003220000000", "0152206500900000"); 
  } else if (trainConfig == 10) { 
    cuts.AddCut("55600013", "04200009297002003220000000", "0152206500900000");
  } else if (trainConfig == 11) { 
    cuts.AddCut("56800013", "04200009297002003220000000", "0152206500900000");    
  } else if (trainConfig == 12) { 
    cuts.AddCut("56700013", "04200009297002003220000000", "0152206500900000"); 
  } else if (trainConfig == 13) { 
    cuts.AddCut("57800013", "04200009297002003220000000", "0152206500900000"); 
  } else if (trainConfig == 14) { 
    cuts.AddCut("46900013", "04200009297002003220000000", "0152206500900000");
  } else if (trainConfig == 15) { 
    cuts.AddCut("58900013", "04200009297002003220000000", "0152206500900000");    
  } else  if (trainConfig == 16){ 
    cuts.AddCut("60100013", "00200009247602008250400000", "0152501500000000"); 
  } else if (trainConfig == 17) { 
    cuts.AddCut("61200013", "00200009247602008250400000", "0152501500000000"); 
  } else if (trainConfig == 18) { 
    cuts.AddCut("50100013", "00200009247602008250400000", "0152501500000000"); 
  } else if (trainConfig == 19) { 
    cuts.AddCut("50200013", "00200009247602008250400000", "0152501500000000");    
  } else if (trainConfig == 20) { 
    cuts.AddCut("51200013", "00200009247602008250400000", "0152501500000000");    
  } else if (trainConfig == 21) { 
    cuts.AddCut("52400013", "00200009247602008250400000", "0152501500000000");       
  } else if (trainConfig == 22) {    
    cuts.AddCut("54600013", "00200009247602008250400000", "0152501500000000"); 
  } else if (trainConfig == 23) {    
    cuts.AddCut("54800013", "00200009247602008250400000", "0152501500000000");    
  } else if (trainConfig == 24) {    
    cuts.AddCut("54500013", "00200009247602008250400000", "0152501500000000"); 
  } else if (trainConfig == 25) { 
    cuts.AddCut("55600013", "00200009247602008250400000", "0152501500000000");
  } else if (trainConfig == 26) { 
    cuts.AddCut("56800013", "00200009247602008250400000", "0152501500000000");    
  } else if (trainConfig == 27) { 
    cuts.AddCut("56700013", "00200009247602008250400000", "0152501500000000"); 
  } else if (trainConfig == 28) { 
    cuts.AddCut("57800013", "00200009247602008250400000", "0152501500000000"); 
  } else if (trainConfig == 29) { 
    cuts.AddCut("46900013", "00200009247602008250400000", "0152501500000000");
  } else if (trainConfig == 30) { 
    cuts.AddCut("58900013", "00200009247602008250400000", "0152501500000000");    
  } else if (trainConfig == 31) { 
    cuts.AddCut("50800013", "00200009247602008250400000", "0152501500000000");    
  } else if (trainConfig == 32) {
    cuts.AddCut("52500013", "00200009247602008250400000", "0152501500000000");
  } else if (trainConfig == 33) { 
    cuts.AddCut("53500013", "00200009247602008250400000", "0152501500000000");       
  } else if (trainConfig == 34) { 
    cuts.AddCut("54500013", "00200009247602008250400000", "0152501500000000");       
  } else if (trainConfig == 35) { 
    cuts.AddCut("53400013", "00200009247602008250400000", "0152501500000000");       
  } else if (trainConfig == 36) { 
    cuts.AddCut("52300013", "00200009247602008250400000", "0152501500000000");     
  } else if (trainConfig == 37) { 
    cuts.AddCut("52500013", "00700009247602008250400000", "0152501500000000");               
  } else if (trainConfig == 38) { 
    cuts.AddCut("52400013", "00700009247602008250400000", "0152501500000000");       
  } else if (trainConfig == 39) { //latest std 11h cut
    cuts.AddCut("52400013", "00200009247602008250404000", "0152501500000000");
  } else if (trainConfig == 40) {
    cuts.AddCut("52500013", "00200009247602008250404000", "0152501500000000");

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
  TObjString *Header1 = new TObjString("BOX");
  HeaderList->Add(Header1);
  //    TObjString *Header3 = new TObjString("eta_2");
  //    HeaderList->Add(Header3);
  
  EventCutList->SetOwner(kTRUE);
  AliConvEventCuts **analysisEventCuts = new AliConvEventCuts*[numberOfCuts];
  ConvCutList->SetOwner(kTRUE);
  AliConversionPhotonCuts **analysisCuts = new AliConversionPhotonCuts*[numberOfCuts];
  MesonCutList->SetOwner(kTRUE);
  AliConversionMesonCuts **analysisMesonCuts = new AliConversionMesonCuts*[numberOfCuts];

  for(Int_t i = 0; i<numberOfCuts; i++){
    analysisEventCuts[i] = new AliConvEventCuts();
    analysisEventCuts[i]->SetTriggerMimicking(enableTriggerMimicking);
    analysisEventCuts[i]->SetTriggerOverlapRejecion(enableTriggerOverlapRej);
    analysisEventCuts[i]->SetMaxFacPtHard(maxFacPtHard);
    analysisEventCuts[i]->SetV0ReaderName(V0ReaderName);
    analysisEventCuts[i]->SetLightOutput(runLightOutput);
    analysisEventCuts[i]->InitializeCutsFromCutString((cuts.GetEventCut(i)).Data());
        
    EventCutList->Add(analysisEventCuts[i]);
    analysisEventCuts[i]->SetFillCutHistograms("",kFALSE);
    
    analysisCuts[i] = new AliConversionPhotonCuts();
    analysisCuts[i]->SetIsHeavyIon(isHeavyIon);
    analysisCuts[i]->SetV0ReaderName(V0ReaderName);
    analysisCuts[i]->SetLightOutput(runLightOutput);
    analysisCuts[i]->InitializeCutsFromCutString((cuts.GetPhotonCut(i)).Data());
    ConvCutList->Add(analysisCuts[i]);
    analysisCuts[i]->SetFillCutHistograms("",kFALSE);

    analysisMesonCuts[i] = new AliConversionMesonCuts();
    analysisMesonCuts[i]->SetLightOutput(runLightOutput);
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
  task->SetDoPlotVsCentrality(kTRUE);

  //connect containers
  AliAnalysisDataContainer *coutput =
    mgr->CreateContainer(Form("GammaConvV1_%i",trainConfig), TList::Class(),
              AliAnalysisManager::kOutputContainer,Form("GammaConvV1_%i.root",trainConfig));

  mgr->AddTask(task);
  mgr->ConnectInput(task,0,cinput);
  mgr->ConnectOutput(task,1,coutput);

  return;

}
