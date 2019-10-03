/////////////////////////////////////////////////////////////////////////////////////////////
//
// AddTaskGammaConvFlow_pp macro
// Author: Friederike Bock, Lucas Altenkaemper, Mike Sas
//         Commented where necessary
/////////////////////////////////////////////////////////////////////////////////////////////

class AliAnalysisDataContainer;
class AliFlowTrackCuts;
class AliFlowTrackSimpleCuts;
class AliFlowEventCuts;
class AliFlowEventSimpleCuts;

class CutHandlerConvFlow{
  public:
    CutHandlerConvFlow(Int_t nMax=10){
      nCuts=0; nMaxCuts=nMax; validCuts = true;
      eventCutArray = new TString[nMaxCuts]; photonCutArray = new TString[nMaxCuts];
      for(Int_t i=0; i<nMaxCuts; i++) {eventCutArray[i] = ""; photonCutArray[i] = "";}
    }

    void AddCut(TString eventCut, TString photonCut){
      if(nCuts>=nMaxCuts) {cout << "ERROR in CutHandlerConvFlow: Exceeded maximum number of cuts!" << endl; validCuts = false; return;}
      if( eventCut.Length()!=8 || photonCut.Length()!=26) {cout << "ERROR in CutHandlerConvFlow: Incorrect length of cut string!" << endl; validCuts = false; return;}
      eventCutArray[nCuts]=eventCut; photonCutArray[nCuts]=photonCut;
      nCuts++;
      return;
    }
    Bool_t AreValid(){return validCuts;}
    Int_t GetNCuts(){if(validCuts) return nCuts; else return 0;}
    TString GetEventCut(Int_t i){if(validCuts&&i<nMaxCuts&&i>=0) return eventCutArray[i]; else{cout << "ERROR in CutHandlerConvFlow: GetEventCut wrong index i" << endl;return "";}}
    TString GetPhotonCut(Int_t i){if(validCuts&&i<nMaxCuts&&i>=0) return photonCutArray[i]; else {cout << "ERROR in CutHandlerConvFlow: GetPhotonCut wrong index i" << endl;return "";}}
  private:
    Bool_t validCuts;
    Int_t nCuts; Int_t nMaxCuts;
    TString* eventCutArray;
    TString* photonCutArray;
};

void AddTask_GammaConvFlow_pPb(
                                TString uniqueID              = "",
                                Int_t trainConfig             = 1,                            //change different set of cuts
                                Int_t enableQAMesonTask       = 0,                            //enable QA in AddTask_GammaConvFlow_PbPb2
                                Int_t enableQAPhotonTask      = 0,                            // enable additional QA task
                                TString cutnumberAODBranch    = "100000006008400000001500000",
                                Bool_t debug                  = kFALSE,
                                Bool_t UseMassSel             = kFALSE,
                                Float_t MinMass               = 0,
                                Float_t MaxMass               = 0.2,
                                Bool_t UseKappaSel            = kFALSE,
                                Float_t MinKappa              = 0,
                                Float_t MaxKappa              = 10,
                                Bool_t isMC                   = kFALSE,
                                TString additionalTrainConfig = "0"                           // additional counter for trainconfig, always has to be last parameter
                              ) {

  Int_t isHeavyIon = 2;

  // make use of train subwagon feature
  if (additionalTrainConfig.Atoi() > 0){
    trainConfig = trainConfig + additionalTrainConfig.Atoi();
  }



  // ================== GetAnalysisManager ===============================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
      Error(Form("AddTask_GammaConvV1_%i",trainConfig), "No analysis manager found.");
      return ;
  }

  // ================== GetInputEventHandler =============================
  AliVEventHandler *inputHandler=mgr->GetInputEventHandler();

  //=========  Set Cutnumber for V0Reader ================================
  TString cutnumberPhoton = "06000008000100001500000000";
  TString cutnumberEvent = "80000103";
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  //========= Add V0 Reader to  ANALYSIS manager if not yet existent =====
  TString V0ReaderName = Form("V0ReaderV1_%s_%s",cutnumberEvent.Data(),cutnumberPhoton.Data());
  if( !(AliV0ReaderV1*)mgr->GetTask(V0ReaderName.Data()) ){
      AliV0ReaderV1 *fV0ReaderV1 = new AliV0ReaderV1(V0ReaderName.Data());

      fV0ReaderV1->SetUseOwnXYZCalculation(kTRUE);
      fV0ReaderV1->SetCreateAODs(kFALSE);// AOD Output
      fV0ReaderV1->SetUseAODConversionPhoton(kTRUE);
      fV0ReaderV1->SetUseMassToZero(kFALSE);
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


  CutHandlerConvFlow cuts;

  if (trainConfig == 1){
    cuts.AddCut("80000013", "00200009000000008260400000");
  } else if (trainConfig == 2){
    cuts.AddCut("80200013", "00200009000000008260400000");
  } else if (trainConfig == 3){
    cuts.AddCut("82400013", "00200009000000008260400000");
  } else if (trainConfig == 4){
    cuts.AddCut("84600013", "00200009000000008260400000");
  } else if (trainConfig == 5){
    cuts.AddCut("86000013", "00200009000000008260400000");
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

  const int numberOfCuts = cuts.GetNCuts();

  //================================================
  //========= Add task to the ANALYSIS manager =====
  //================================================
  AliAnalysisTaskGammaConvFlow *task=NULL;
  task= new AliAnalysisTaskGammaConvFlow(Form("GammaConvV1_%i",trainConfig),numberOfCuts);
  task->SetIsHeavyIon(isHeavyIon);
  task->AliAnalysisTaskGammaConvFlow::SetIsMC(isMC);
  task->SetV0ReaderName(V0ReaderName);

  AliFlowTrackCuts* cutsRP = new AliFlowTrackCuts(Form("RFPcuts%s",uniqueID));
  if(!cutsRP) {
      if(debug) cout << " Fatal error: no RP cuts found, could be a library problem! " << endl;
      return 0x0;
  }
  cutsRP = cutsRP->GetStandardVZEROOnlyTrackCuts(); // select vzero tracks
  cutsRP->SetVZEROgainEqualizationPerRing(kFALSE);
  cutsRP->SetApplyRecentering(kTRUE);

  task->SetRPCuts(cutsRP);

  if(UseMassSel==kTRUE)  task->SetMassWindow(MinMass,MaxMass);
  if(UseKappaSel==kTRUE) task->SetKappaWindow(MinKappa,MaxKappa);

  //======================================================================
  TList *EventCutList = new TList();
  TList *ConvCutList = new TList();

  TList *HeaderList = new TList();
  TObjString *Header1 = new TObjString("BOX");
  HeaderList->Add(Header1);
  //    TObjString *Header3 = new TObjString("eta_2");
  //    HeaderList->Add(Header3);

  EventCutList->SetOwner(kTRUE);
  AliConvEventCuts **analysisEventCuts = new AliConvEventCuts*[numberOfCuts];
  ConvCutList->SetOwner(kTRUE);
  AliConversionPhotonCuts **analysisCuts = new AliConversionPhotonCuts*[numberOfCuts];

  for(Int_t i = 0; i<numberOfCuts; i++){
    analysisEventCuts[i] = new AliConvEventCuts();
    analysisEventCuts[i]->InitializeCutsFromCutString((cuts.GetEventCut(i)).Data());
    analysisEventCuts[i]->SetV0ReaderName(V0ReaderName);
    EventCutList->Add(analysisEventCuts[i]);
    analysisEventCuts[i]->SetFillCutHistograms("",kFALSE);

    analysisCuts[i] = new AliConversionPhotonCuts();
    analysisCuts[i]->InitializeCutsFromCutString((cuts.GetPhotonCut(i)).Data());
    analysisCuts[i]->SetV0ReaderName(V0ReaderName);
    ConvCutList->Add(analysisCuts[i]);
    analysisCuts[i]->SetFillCutHistograms("",kFALSE);

    analysisEventCuts[i]->SetAcceptedHeader(HeaderList);
  }

  task->SetEventCutList(numberOfCuts,EventCutList);
  task->SetConversionCutList(numberOfCuts,ConvCutList);
  task->SetDoMesonAnalysis(kFALSE);
  task->SetDoMesonQA(enableQAMesonTask); //Attention new switch for Pi0 QA
  task->SetDoPhotonQA(enableQAPhotonTask);  //Attention new switch small for Photon QA

  //connect containers
  AliAnalysisDataContainer *coutput =
  mgr->CreateContainer(Form("GammaConvV1_%i",trainConfig), TList::Class(),
                      AliAnalysisManager::kOutputContainer,Form("GammaConvFlow_%i.root",trainConfig));

  mgr->AddTask(task);
  mgr->ConnectInput(task,0,cinput);
  mgr->ConnectOutput(task,1,coutput);

  return;

}
