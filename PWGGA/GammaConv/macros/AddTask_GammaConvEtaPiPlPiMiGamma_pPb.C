void AddTask_GammaConvEtaPiPlPiMiGamma_pPb(
                                            Int_t   trainConfig               = 1,
                                            Bool_t  isMC                      = kFALSE, //run MC
                                            Bool_t  enableQAMesonTask         = kTRUE, //enable QA in AliAnalysisTaskEtaToPiPlPiMiGamma
                                            TString fileNameInputForWeighting = "MCSpectraInput.root", // path to file for weigting input
                                            Bool_t  doWeighting               = kFALSE,  //enable Weighting
                                            TString generatorName             = "HIJING",
                                            TString cutnumberAODBranch        = "000000006008400001001500000"
                                          ) {

  Int_t isHeavyIon = 2;

  // ================== GetAnalysisManager ===============================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error(Form("AddTask_GammaConvEtaPiPlPiMiGamma_pPb_%i",trainConfig), "No analysis manager found.");
    return ;
  }

  // ================== GetInputEventHandler =============================
  AliVEventHandler *inputHandler=mgr->GetInputEventHandler();

  //========= Add PID Reponse to ANALYSIS manager ====
  if(!(AliPIDResponse*)mgr->GetTask("PIDResponseTask")){
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
    AddTaskPIDResponse(isMC);
  }

  //=========  Set Cutnumber for V0Reader ================================
  TString cutnumberPhoton = "06000008400100001500000000";
  TString cutnumberEvent  = "80000003";
  TString PionCuts        = "000000200";            //Electron Cuts

  Bool_t doEtaShift = kFALSE;

  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();

  //========= Add V0 Reader to  ANALYSIS manager if not yet existent =====
  TString V0ReaderName = Form("V0ReaderV1_%s_%s",cutnumberEvent.Data(),cutnumberPhoton.Data());
  if( !(AliV0ReaderV1*)mgr->GetTask(V0ReaderName.Data()) ){
    AliV0ReaderV1 *fV0ReaderV1 = new AliV0ReaderV1(V0ReaderName.Data());

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
      fV0ReaderV1->AliV0ReaderV1::SetDeltaAODBranchName(Form("GammaConv_%s_gamma",cutnumberAODBranch.Data()));
    }
    fV0ReaderV1->Init();

    AliLog::SetGlobalLogLevel(AliLog::kInfo);

    //connect input V0Reader
    mgr->AddTask(fV0ReaderV1);
    mgr->ConnectInput(fV0ReaderV1,0,cinput);
  }

  //================================================
  //========= Add Electron Selector ================


  if( !(AliPrimaryPionSelector*)mgr->GetTask("PionSelector") ){

    AliPrimaryPionSelector *fPionSelector = new AliPrimaryPionSelector("PionSelector");
    // Set AnalysisCut Number

    AliPrimaryPionCuts *fPionCuts=0;
    if( PionCuts!=""){
      fPionCuts= new AliPrimaryPionCuts(PionCuts.Data(),PionCuts.Data());
      if(fPionCuts->InitializeCutsFromCutString(PionCuts.Data())){
        fPionSelector->SetPrimaryPionCuts(fPionCuts);
        fPionCuts->SetFillCutHistograms("",kTRUE);

      }
    }

    fPionSelector->Init();
    mgr->AddTask(fPionSelector);

    AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer();

    //connect input V0Reader
    mgr->ConnectInput (fPionSelector,0,cinput1);

  }



  AliAnalysisTaskEtaToPiPlPiMiGamma *task=NULL;

  task= new AliAnalysisTaskEtaToPiPlPiMiGamma(Form("GammaConvEtaPiPlPiMiGamma_%i",trainConfig));

  task->SetIsHeavyIon(2);
  task->SetIsMC(isMC);
  task->SetV0ReaderName(V0ReaderName);

  // Cut Numbers to use in Analysis
  Int_t numberOfCuts = 1;

  TString *eventCutArray = new TString[numberOfCuts];
  TString *ConvCutarray    = new TString[numberOfCuts];
  TString *PionCutarray    = new TString[numberOfCuts];
  TString *MesonCutarray   = new TString[numberOfCuts];

  Bool_t doEtaShiftIndCuts = kFALSE;
  TString stringShift = "";

  // Shifting in pPb direction

  doEtaShiftIndCuts = kTRUE;
  stringShift = "pPb";

  if( trainConfig == 1 ) {
    eventCutArray[ 0] = "80000113"; ConvCutarray[0] = "00200009117000008260400000"; PionCutarray[0] = "000000400"; MesonCutarray[0] = "0103503000000000"; //standard cut Pi0 PbPb 00-100
  } else if( trainConfig == 1 ) {
    eventCutArray[ 0] = "80000113"; ConvCutarray[0] = "00200009117000008260400000"; PionCutarray[0] = "000000403"; MesonCutarray[0] = "0103503000000000"; //standard cut Pi0 PbPb 00-100
  } else if( trainConfig == 1 ) {
    eventCutArray[ 0] = "80000113"; ConvCutarray[0] = "00200009117000008260400000"; PionCutarray[0] = "000000404"; MesonCutarray[0] = "0103503000000000"; //standard cut Pi0 PbPb 00-100
  } else if( trainConfig == 1 ) {
    eventCutArray[ 0] = "80000113"; ConvCutarray[0] = "00200009117000008260400000"; PionCutarray[0] = "000000405"; MesonCutarray[0] = "0103503000000000"; //standard cut Pi0 PbPb 00-100
  }

  TList *EventCutList = new TList();
  TList *ConvCutList  = new TList();
  TList *MesonCutList = new TList();
  TList *PionCutList  = new TList();

  TList *HeaderList = new TList();
  TObjString *Header1 = new TObjString("pi0_1");
  HeaderList->Add(Header1);
  TObjString *Header3 = new TObjString("eta_2");
  HeaderList->Add(Header3);

  if (periodNameV0Reader.Contains("LHC17g6a2") || periodNameV0Reader.Contains("LHC17g6a3") ){
    TObjString *HeaderPMB = new TObjString("Dpmjet_0");
    TObjString *HeaderP8J = new TObjString("Pythia8JetsGammaTrg_1");
    if (doWeightingPart==4) { // all headers
      HeaderList->Add(HeaderPMB);
      HeaderList->Add(HeaderP8J);
    } else if (doWeightingPart==5) { // only MB header
      HeaderList->Add(HeaderPMB);
    } else { // only JJ header
      HeaderList->Add(HeaderP8J);
    }
  }
  
  EventCutList->SetOwner(kTRUE);
  AliConvEventCuts **analysisEventCuts = new AliConvEventCuts*[numberOfCuts];
  ConvCutList->SetOwner(kTRUE);
  AliConversionPhotonCuts **analysisCuts = new AliConversionPhotonCuts*[numberOfCuts];
  MesonCutList->SetOwner(kTRUE);
  AliConversionMesonCuts **analysisMesonCuts   = new AliConversionMesonCuts*[numberOfCuts];
  PionCutList->SetOwner(kTRUE);
  AliPrimaryPionCuts **analysisPionCuts     = new AliPrimaryPionCuts*[numberOfCuts];

  for(Int_t i = 0; i<numberOfCuts; i++){
    analysisEventCuts[i] = new AliConvEventCuts();
    analysisEventCuts[i]->SetV0ReaderName(V0ReaderName);
    analysisEventCuts[i]->InitializeCutsFromCutString(eventCutArray[i].Data());
    EventCutList->Add(analysisEventCuts[i]);
    analysisEventCuts[i]->SetFillCutHistograms("",kFALSE);

    analysisCuts[i] = new AliConversionPhotonCuts();
    analysisCuts[i]->SetV0ReaderName(V0ReaderName);
    if( ! analysisCuts[i]->InitializeCutsFromCutString(ConvCutarray[i].Data()) ) {
        cout<<"ERROR: analysisCuts [" <<i<<"]"<<endl;
        return 0;
    } else {
      ConvCutList->Add(analysisCuts[i]);
      analysisCuts[i]->SetFillCutHistograms("",kFALSE);

    }

    analysisMesonCuts[i] = new AliConversionMesonCuts();

    if( ! analysisMesonCuts[i]->InitializeCutsFromCutString(MesonCutarray[i].Data()) ) {
      cout<<"ERROR: analysisMesonCuts [ " <<i<<" ] "<<endl;
      return 0;
    } else {
      MesonCutList->Add(analysisMesonCuts[i]);
      analysisMesonCuts[i]->SetFillCutHistograms("");
    }
    analysisEventCuts[i]->SetAcceptedHeader(HeaderList);

    TString cutName( Form("%s_%s_%s_%s",eventCutArray[i].Data(), ConvCutarray[i].Data(),PionCutarray[i].Data(),MesonCutarray[i].Data() ) );
    analysisPionCuts[i] = new AliPrimaryPionCuts();
    if( !analysisPionCuts[i]->InitializeCutsFromCutString(PionCutarray[i].Data())) {
      cout<< "ERROR:  analysisPionCuts [ " <<i<<" ] "<<endl;
      return 0;
    } else {
      PionCutList->Add(analysisPionCuts[i]);
      analysisPionCuts[i]->SetFillCutHistograms("",kFALSE,cutName);
    }


  }

  task->SetEventCutList(numberOfCuts,EventCutList);
  task->SetConversionCutList(numberOfCuts,ConvCutList);
  task->SetMesonCutList(MesonCutList);
  task->SetPionCutList(PionCutList);

  task->SetMoveParticleAccordingToVertex(kTRUE);

  if(enableQAMesonTask) task->SetDoMesonQA(kTRUE);

  //connect containers
  AliAnalysisDataContainer *coutput =
  mgr->CreateContainer(Form("GammaConvEtaPiPlPiMiGamma_%i",trainConfig), TList::Class(),
              AliAnalysisManager::kOutputContainer,Form("GammaConvEtaPiPlPiMiGamma_%i.root",trainConfig));

  mgr->AddTask(task);
  mgr->ConnectInput(task,0,cinput);
  mgr->ConnectOutput(task,1,coutput);

  return;

}
