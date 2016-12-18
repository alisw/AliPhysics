AliAnalysisTaskSE* AddTaskSigma0Run2(Bool_t isAOD=kFALSE, Bool_t isMC=kFALSE, Bool_t isHeavyIon=kFALSE, Bool_t isQA=kTRUE)
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskSigma0Run2()", "No analysis manager found.");
    return 0x0;
  }  
  
  // ================== GetInputEventHandler =============================
  AliVEventHandler *inputHandler=mgr->GetInputEventHandler();
  
  //========= Add PID Reponse to ANALYSIS manager ====
  if(!(AliPIDResponse*)mgr->GetTask("PIDResponseTask")){
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
    AddTaskPIDResponse(isMC);
  }
  
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  
  
  
 //=========  Set Cutnumber for V0Reader ================================
  TString cutnumberPhoton   = "06000008400100001500000000";
  TString cutnumberEvent    = "00000003";
  TString cutnumberV0       = "1600";
  TString PionCuts          = "000000200";            //Electron Cuts
  
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  
  //========= Add V0 Reader to analysis manager =========================
  TString V0ReaderName = Form("V0ReaderV1_%s_%s",cutnumberEvent.Data(),cutnumberPhoton.Data());
  TString V0ReaderStrangeName = Form("V0ReaderStrange_%s_%s",cutnumberEvent.Data(),cutnumberV0.Data());
  
  if( !(AliV0ReaderV1*)mgr->GetTask(V0ReaderName.Data()) && !(AliV0ReaderStrange*)mgr->GetTask(V0ReaderStrangeName.Data())){

    AliV0ReaderV1 *fV0ReaderV1 = new AliV0ReaderV1(V0ReaderName.Data());
    fV0ReaderV1->SetUseOwnXYZCalculation(kTRUE);
    fV0ReaderV1->SetCreateAODs(kFALSE);// AOD Output
    fV0ReaderV1->SetUseAODConversionPhoton(kTRUE);
    
    AliV0ReaderStrange *fV0ReaderStrange = new AliV0ReaderStrange(V0ReaderStrangeName.Data());
    
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
        fV0ReaderStrange->SetEventCuts(fEventCuts);
        fEventCuts->SetFillCutHistograms("",kTRUE);
      }
    }

    // Set AnalysisCut Number
    AliConversionPhotonCuts *fCuts=NULL;
    if(cutnumberPhoton!=""){
      fCuts= new AliConversionPhotonCuts(cutnumberPhoton.Data(),cutnumberPhoton.Data());
      fCuts->SetPreSelectionCutFlag(kTRUE);
//       fCuts->SetIsHeavyIon(isHeavyIon);
      fCuts->SetIsHeavyIon(kFALSE);
      fCuts->SetV0ReaderName(V0ReaderName);
      if(fCuts->InitializeCutsFromCutString(cutnumberPhoton.Data())){
        fV0ReaderV1->SetConversionCuts(fCuts);
        fCuts->SetFillCutHistograms("",kTRUE);
      }
    }
    
    AliV0CutsStrange *fV0CutsStrange = NULL;
    if(cutnumberV0!=""){
      fV0CutsStrange = new AliV0CutsStrange(cutnumberV0.Data(), cutnumberV0.Data());
      fV0CutsStrange->SetIsQA(isQA);
      if(fV0CutsStrange->InitializeCutsFromCutString(cutnumberV0.Data())){
        fV0ReaderStrange->SetV0Cuts(fV0CutsStrange);
        fV0CutsStrange->SetFillCutHistograms("",kTRUE);
      }
    }

//     if(inputHandler->IsA()==AliAODInputHandler::Class()){
    // AOD mode
//       fV0ReaderV1->SetDeltaAODBranchName(Form("GammaConv_%s_gamma",cutnumberAODBranch.Data()));
//     }


    fV0ReaderV1->Init();
    fV0ReaderStrange->Init();

    AliLog::SetGlobalLogLevel(AliLog::kFatal);

    //connect input V0Reader
    mgr->AddTask(fV0ReaderV1);
    mgr->ConnectInput(fV0ReaderV1,0,cinput);
    mgr->AddTask(fV0ReaderStrange);
    mgr->ConnectInput(fV0ReaderStrange,0,cinput);
  }
    
  AliAnalysisTaskSigma0Run2 *task = new AliAnalysisTaskSigma0Run2("TaskSigma0Run2");
  task->SetV0ReaderName(V0ReaderName);
  task->SetV0ReaderStrangeName(V0ReaderStrangeName);
  task->SetIsHeavyIon(isHeavyIon);
  task->SetIsMC(isMC);
  task->SetIsQA(isQA);
  mgr->AddTask(task);
  
  TString containerName = mgr->GetCommonFileName();
  containerName += ":Sigma0Analysis_Run2";
  
  AliAnalysisDataContainer *cOutputList = mgr->CreateContainer("histos", TList::Class(), AliAnalysisManager::kOutputContainer, containerName.Data() );
  
  mgr->ConnectInput(task,  0, cinput);
  mgr->ConnectOutput(task, 1, cOutputList);
  
  return task;
}

