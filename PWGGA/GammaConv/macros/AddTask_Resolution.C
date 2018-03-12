void AddTask_Resolution(  TString   V0ReaderEventCutNumber  = "80000003",
                          TString   V0ReaderPhotonCutNumber = "060000084001001500000000",
                          TString   TaskEventCutnumber      = "80000113",
                          TString   TaskPhotonCutnumber     = "092000092170008260400000",
                          Bool_t    isMC                    = kTRUE,
                          Int_t     IsHeavyIon              = 0,
                          TString   cutnumberAODBranch      = "0000000060084001001500000",
                          Bool_t    doEtaShiftV0Reader      = kFALSE,
                          Bool_t    enableV0findingEffi     = kFALSE                        // enables V0finding efficiency histograms
                      ){

  //get the current analysis manager
  // ================== GetAnalysisManager ===============================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error(Form("AddTask_GammaConvV1_%i",trainConfig), "No analysis manager found.");
    return ;
  }

  // ================== GetInputEventHandler =============================
  AliVEventHandler *inputHandler=mgr->GetInputEventHandler();

  //========= Add PID Reponse to ANALYSIS manager ====
  if(!(AliPIDResponse*)mgr->GetTask("PIDResponseTask")){
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
    AddTaskPIDResponse(isMC);
  }

  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  TString V0ReaderName = Form("V0ReaderV1_%s_%s",cutnumberEvent.Data(),cutnumberPhoton.Data());
  if( !(AliV0ReaderV1*)mgr->GetTask(V0ReaderName.Data()) ){
    AliV0ReaderV1 *fV0ReaderV1 = new AliV0ReaderV1(V0ReaderName.Data());

    fV0ReaderV1->SetUseOwnXYZCalculation(kTRUE);
    fV0ReaderV1->SetCreateAODs(kFALSE);// AOD Output
    fV0ReaderV1->SetUseAODConversionPhoton(kTRUE);
    fV0ReaderV1->SetProduceV0FindingEfficiency(enableV0findingEffi);

    AliConvEventCuts *fEventCuts=NULL;
    if(V0ReaderEventCutNumber!=""){
      fEventCuts= new AliConvEventCuts(V0ReaderEventCutNumber.Data(),V0ReaderEventCutNumber.Data());
      fEventCuts->SetPreSelectionCutFlag(kTRUE);
            fEventCuts->SetV0ReaderName(V0ReaderName);
      if(fEventCuts->InitializeCutsFromCutString(V0ReaderEventCutNumber.Data())){
        fV0ReaderV1->SetEventCuts(fEventCuts);
        fEventCuts->SetFillCutHistograms("",kTRUE);
        if (IsHeavyIon==2){
          fEventCuts->SelectCollisionCandidates(AliVEvent::kINT7);
          fEventCuts->DoEtaShift(doEtaShiftV0Reader);
        }
      }
    }

    // Set AnalysisCut Number
    AliConversionPhotonCuts *fCuts=NULL;
    if(V0ReaderPhotonCutNumber!=""){
      fCuts= new AliConversionPhotonCuts(V0ReaderPhotonCutNumber.Data(),V0ReaderPhotonCutNumber.Data());
      fCuts->SetPreSelectionCutFlag(kTRUE);
      fCuts->SetIsHeavyIon(IsHeavyIon);
      fCuts->SetV0ReaderName(V0ReaderName);
      if(fCuts->InitializeCutsFromCutString(V0ReaderPhotonCutNumber.Data())){
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

  } else {
    Error("AddTask_V0ReaderV1", "Cannot execute AddTask, V0ReaderV1 already exists.");
  }

  AliConvEventCuts *analysisEventCuts = new AliConvEventCuts();
  analysisEventCuts->SetV0ReaderName(V0ReaderName);
  analysisEventCuts->InitializeCutsFromCutString(TaskEventCutnumber.Data());
  analysisEventCuts->SetFillCutHistograms("",kFALSE);

  AliConversionPhotonCuts *analysisCuts = new AliConversionPhotonCuts();
  analysisCuts->SetV0ReaderName(V0ReaderName);
  analysisCuts->InitializeCutsFromCutString(TaskPhotonCutnumber.Data());
  analysisCuts->SetFillCutHistograms("",kFALSE);

  AliAnalysisTaskResolution *fResolution= new AliAnalysisTaskResolution(Form("%s_%s_Resolution",(analysisEventCuts->GetCutNumber()).Data(), (analysisCuts->GetCutNumber()).Data()));
  fResolution->SetEventCuts(analysisEventCuts,IsHeavyIon);
  fResolution->SetConversionCuts(analysisCuts,IsHeavyIon);
  fResolution->SetIsMC(isMC);
  fResolution->SetV0ReaderName(V0ReaderName);
  mgr->AddTask(fResolution);

  AliAnalysisDataContainer *coutput1 =
  mgr->CreateContainer(Form("GammaConvResolution_%s_%s", TaskEventCutnumber.Data(), TaskPhotonCutnumber.Data()), TList::Class(),
            AliAnalysisManager::kOutputContainer, Form("GammaConv_Resolution_%s_%s.root", TaskEventCutnumber.Data(), TaskPhotonCutnumber.Data()));

  AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer();
  mgr->ConnectInput(fResolution,  0, cinput1 );
  mgr->ConnectOutput (fResolution,  1, coutput1);
  //connect containers
  return;
}

