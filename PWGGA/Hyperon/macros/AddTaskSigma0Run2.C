AliAnalysisTaskSE *AddTaskSigma0Run2(bool isAOD = false, bool isMC = false,
                                     bool isHeavyIon = false, bool isQA = true,
                                     bool isRun1 = false) {
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskSigma0Run2()", "No analysis manager found.");
    return 0x0;
  }

  // ================== GetInputEventHandler =============================
  AliVEventHandler *inputHandler = mgr->GetInputEventHandler();

  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();

  //========= Add PID Reponse to ANALYSIS manager ========================
  if (isMC) {
    // IMPORTANT - SET WHEN USING DIFFERENT PASS -  1 here for pass1, 2 for
    // pass2
    if (isRun1)
      AliAnalysisTaskPIDResponse *pidResponse =
          reinterpret_cast<AliAnalysisTaskPIDResponse *>(
              gInterpreter->ExecuteMacro("$ALICE_ROOT/ANALYSIS/macros/"
                                         "AddTaskPIDResponse.C (kTRUE, kTRUE, "
                                         "kTRUE, \"4\")"));
    else
      AliAnalysisTaskPIDResponse *pidResponse =
          reinterpret_cast<AliAnalysisTaskPIDResponse *>(
              gInterpreter->ExecuteMacro("$ALICE_ROOT/ANALYSIS/macros/"
                                         "AddTaskPIDResponse.C (kTRUE, kTRUE, "
                                         "kTRUE, \"1\")"));
  } else {
    AliAnalysisTaskPIDResponse *pidResponse =
        reinterpret_cast<AliAnalysisTaskPIDResponse *>(
            gInterpreter->ExecuteMacro(
                "$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C)"));
  }

  //=========  Set Cutnumber for V0Reader ================================
  TString cutnumberPhoton = "00200008400020002282000000";
  //      00000008800020002282000000 tighter TPC dEdx cut
  //      00000088400020002282000000 pt,ele > 0.02
  //      00000008400020002282000000 pt,ele > 0.05
  TString cutnumberEvent = "00000000";
  TString periodNameV0Reader = "";
  Bool_t enableV0findingEffi = kFALSE;
  Bool_t runLightOutput = kFALSE;
  TString cutnumberAODBranch = "00000003_06000008400100001000000000";

  //========= Add V0 Reader to  ANALYSIS manager if not yet existent =====
  TString V0ReaderName =
      Form("V0ReaderV1_%s_%s", cutnumberEvent.Data(), cutnumberPhoton.Data());
  AliConvEventCuts *fEventCuts = NULL;

  if (!(AliV0ReaderV1 *)mgr->GetTask(V0ReaderName.Data())) {
    AliV0ReaderV1 *fV0ReaderV1 = new AliV0ReaderV1(V0ReaderName.Data());
    if (periodNameV0Reader.CompareTo("") != 0)
      fV0ReaderV1->SetPeriodName(periodNameV0Reader);
    fV0ReaderV1->SetUseOwnXYZCalculation(kTRUE);
    fV0ReaderV1->SetCreateAODs(kFALSE);  // AOD Output
    fV0ReaderV1->SetUseAODConversionPhoton(kTRUE);
    fV0ReaderV1->SetProduceV0FindingEfficiency(enableV0findingEffi);

    if (!mgr) {
      Error("AddTask_V0ReaderV1", "No analysis manager found.");
      return nullptr;
    }

    if (cutnumberEvent != "") {
      fEventCuts =
          new AliConvEventCuts(cutnumberEvent.Data(), cutnumberEvent.Data());
      fEventCuts->SetPreSelectionCutFlag(kTRUE);
      fEventCuts->SetV0ReaderName(V0ReaderName);
      if (periodNameV0Reader.CompareTo("") != 0)
        fEventCuts->SetPeriodEnum(periodNameV0Reader);
      fV0ReaderV1->SetEventCuts(fEventCuts);
    }

    // Set AnalysisCut Number
    AliConversionPhotonCuts *fCuts = NULL;
    if (cutnumberPhoton != "") {
      fCuts = new AliConversionPhotonCuts(cutnumberPhoton.Data(),
                                          cutnumberPhoton.Data());
      fCuts->SetPreSelectionCutFlag(kTRUE);
      fCuts->SetIsHeavyIon(isHeavyIon);
      fCuts->SetV0ReaderName(V0ReaderName);
      fCuts->SetLightOutput(runLightOutput);
      if (fCuts->InitializeCutsFromCutString(cutnumberPhoton.Data())) {
        fV0ReaderV1->SetConversionCuts(fCuts);
        fCuts->SetFillCutHistograms("", kTRUE);
      }
    }

    if (inputHandler->IsA() == AliAODInputHandler::Class()) {
      fV0ReaderV1->SetDeltaAODBranchName(
          Form("GammaConv_%s_gamma", cutnumberAODBranch.Data()));
    }

    fV0ReaderV1->Init();

    AliLog::SetGlobalLogLevel(AliLog::kFatal);

    // connect input V0Reader
    mgr->AddTask(fV0ReaderV1);
    mgr->ConnectInput(fV0ReaderV1, 0, cinput);
  }

  //========= Init subtasks and start analyis ============================

  AliSigma0EventCuts *evCuts = new AliSigma0EventCuts;
  evCuts->SetV0ReaderName(V0ReaderName.Data());
  if (isRun1) {
    evCuts->SetTrigger(AliVEvent::kMB);
    evCuts->SetUseAliEventCuts(false);
    evCuts->SetUseConversionCuts(true);
    evCuts->SetZVertexCut(10.f);
    evCuts->SetNVertexContributors(1);
  } else {
    evCuts->SetTrigger(AliVEvent::kINT7);
    evCuts->SetUseAliEventCuts(true);
    // if the V0 High Multiplicity trigger is selected, an additional cut on the
    // V0 percentile is necessary
    if (evCuts->GetTrigger() == AliVEvent::kHighMultV0)
      evCuts->SetV0Percentile(0.1);
  }
  evCuts->InitCutHistograms();

  AliSigma0SingleParticleCuts *spCuts =
      AliSigma0SingleParticleCuts::ElectronCuts();
  spCuts->SetIsMC(isMC);
  spCuts->SetExtendedQA(isQA);
  spCuts->InitCutHistograms();

  AliSigma0SingleParticleCuts *spCutsProton =
      AliSigma0SingleParticleCuts::DefaultCuts();
  spCutsProton->SetIsMC(isMC);
  spCutsProton->SetExtendedQA(isQA);
  spCutsProton->InitCutHistograms();

  AliSigma0V0Cuts *v0Cuts = AliSigma0V0Cuts::Sigma0Cuts();
  v0Cuts->SetIsMC(isMC);
  v0Cuts->SetPileUpRejection(!isRun1);
  v0Cuts->SetSingleParticleCuts(spCuts);
  v0Cuts->SetExtendedQA(isQA);
  v0Cuts->InitCutHistograms("Sigma0");

  AliSigma0V0Cuts *v0LambdaCuts = AliSigma0V0Cuts::DefaultCuts();
  v0LambdaCuts->SetIsMC(isMC);
  v0LambdaCuts->SetPileUpRejection(!isRun1);
  v0LambdaCuts->SetSingleParticleCuts(spCuts);
  v0LambdaCuts->SetExtendedQA(isQA);
  v0LambdaCuts->InitCutHistograms("Lambda");

  AliSigma0PhotonCuts *photonCuts = AliSigma0PhotonCuts::DefaultCuts();
  photonCuts->SetIsMC(isMC);
  photonCuts->SetV0ReaderName(V0ReaderName.Data());
  photonCuts->InitCutHistograms();

  AliSigma0PhotonMotherCuts *photonMotherCuts =
      AliSigma0PhotonMotherCuts::DefaultCuts();
  photonMotherCuts->SetIsMC(isMC);
  photonMotherCuts->SetSigmaMass(1.192642);
  photonMotherCuts->SetSigmaMassCut(0.005);
  photonMotherCuts->SetSigmaSideband(0.015, 0.05);
  photonMotherCuts->InitCutHistograms();

  AliSigma0EventContainer *evCont = new AliSigma0EventContainer();
  evCont->SetIsMC(isMC);
  evCont->SetExtendedQA(isQA);
  evCont->SetSigmaMass(1.192642);
  evCont->SetSigmaMassCut(0.005);
  evCont->SetSigmaSideband(0.015, 0.05);
  evCont->SetProtonMixingDepth(10);
  evCont->SetLambdaMixingDepth(10);
  evCont->SetPhotonMixingDepth(10);
  evCont->SetSigmaMixingDepth(25);
  //  evCont->SetZvertexBins(10, -10, 2);
  evCont->InitCutHistograms();

  AliAnalysisTaskSigma0Run2 *task =
      new AliAnalysisTaskSigma0Run2("AnalysisTaskSigma0");
  task->SelectCollisionCandidates(evCuts->GetTrigger());
  task->SetV0ReaderName(V0ReaderName.Data());
  task->SetEventCuts(evCuts);
  task->SetSingleParticleCuts(spCuts);
  task->SetProtonCuts(spCutsProton);
  task->SetV0Cuts(v0Cuts);
  task->SetV0LambdaCuts(v0LambdaCuts);
  task->SetPhotonCuts(photonCuts);
  task->SetPhotonMotherCuts(photonMotherCuts);
  task->SetEventContainer(evCont);
  task->SetIsHeavyIon(isHeavyIon);
  task->SetIsMC(isMC);
  mgr->AddTask(task);

  TString containerName = mgr->GetCommonFileName();
  containerName += ":AnalysisTaskSigma0";

  AliAnalysisDataContainer *cOutputList = mgr->CreateContainer(
      "histos", TList::Class(), AliAnalysisManager::kOutputContainer,
      containerName.Data());

  mgr->ConnectInput(task, 0, cinput);
  mgr->ConnectOutput(task, 1, cOutputList);

  return task;
}
