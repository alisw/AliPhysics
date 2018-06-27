AliAnalysisTaskSE *AddTaskSigma0Run2(bool isMC = false, bool isHeavyIon = false,
                                     TString trigger = "kINT7",
                                     const char *cutVariation = "0") {
  TString suffix;
  suffix.Form("%s", cutVariation);

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskSigma0Run2()", "No analysis manager found.");
    return 0x0;
  }

  // ================== GetInputEventHandler =============================
  AliVEventHandler *inputHandler = mgr->GetInputEventHandler();

  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();

  //=========  Set Cutnumber for V0Reader ================================
  TString cutnumberPhoton;
  cutnumberPhoton = "00200008400020002282000000";
  if (suffix == "10") cutnumberPhoton = "00000008800020002282000000";
  if (suffix == "11") cutnumberPhoton = "00000088400020002282000000";
  if (suffix == "12") cutnumberPhoton = "00000008400020002282000000";
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

    fV0ReaderV1->Init();

    AliLog::SetGlobalLogLevel(AliLog::kFatal);

    // connect input V0Reader
    mgr->AddTask(fV0ReaderV1);
    mgr->ConnectInput(fV0ReaderV1, 0, cinput);
  }

  //========= Init subtasks and start analyis ============================
  AliSigma0V0Cuts *v0Cuts = AliSigma0V0Cuts::LambdaCuts();
  v0Cuts->SetIsMC(isMC);
  v0Cuts->SetPID(3122);
  v0Cuts->SetPosPID(AliPID::kProton);
  v0Cuts->SetNegPID(AliPID::kPion);
  AliSigma0V0Cuts *antiv0Cuts = AliSigma0V0Cuts::LambdaCuts();
  antiv0Cuts->SetIsMC(isMC);
  antiv0Cuts->SetPID(-3122);
  antiv0Cuts->SetPosPID(AliPID::kPion);
  antiv0Cuts->SetNegPID(AliPID::kProton);

  if (suffix != "0") {
    v0Cuts->SetLightweight(true);
    antiv0Cuts->SetLightweight(true);
  }
  if (suffix == "1") {
    v0Cuts->SetV0PtMin(0);
    antiv0Cuts->SetV0PtMin(0);
  }
  if (suffix == "2") {
    v0Cuts->SetV0CosPAMin(0.97);
    antiv0Cuts->SetV0CosPAMin(0.97);
  }
  if (suffix == "3") {
    v0Cuts->SetV0CosPAMin(0.95);
    antiv0Cuts->SetV0CosPAMin(0.95);
  }
  if (suffix == "4") {
    v0Cuts->SetTPCclusterMin(0);
    v0Cuts->SetTPCRatioFindable(0.8);
    antiv0Cuts->SetTPCclusterMin(0);
    antiv0Cuts->SetTPCRatioFindable(0.8);
  }
  if (suffix == "5") {
    v0Cuts->SetTPCclusterMin(0);
    v0Cuts->SetTPCRatioFindable(0.6);
    antiv0Cuts->SetTPCclusterMin(0);
    antiv0Cuts->SetTPCRatioFindable(0.6);
  }
  if (suffix == "6") {
    v0Cuts->SetTPCclusterMin(0);
    v0Cuts->SetTPCRatioFindable(0.6);
    antiv0Cuts->SetTPCclusterMin(0);
    antiv0Cuts->SetTPCRatioFindable(0.6);
  }
  if (suffix == "7") {
    v0Cuts->SetV0RadiusMax(200.f);
    v0Cuts->SetV0RadiusMin(0.);
    v0Cuts->SetV0DecayVertexMax(200.f);
    antiv0Cuts->SetV0RadiusMax(200.f);
    antiv0Cuts->SetV0RadiusMin(0.);
    antiv0Cuts->SetV0DecayVertexMax(200.f);
  }
  if (suffix == "8") {
    v0Cuts->SetK0Rejection(0., 0.);
    antiv0Cuts->SetK0Rejection(0., 0.);
  }
  if (suffix == "9") {
    v0Cuts->SetLambdaSelection(1.115683 - 0.008, 1.115683 + 0.008);
    antiv0Cuts->SetLambdaSelection(1.115683 - 0.008, 1.115683 + 0.008);
  }

  v0Cuts->InitCutHistograms(TString("Lambda"));
  antiv0Cuts->InitCutHistograms(TString("AntiLambda"));

  AliSigma0PhotonMotherCuts *sigmaCuts =
      AliSigma0PhotonMotherCuts::DefaultCuts();
  sigmaCuts->SetIsMC(isMC);
  sigmaCuts->SetSigmaMass(1.192642);
  sigmaCuts->SetSigmaMassCut(0.005);
  sigmaCuts->SetSigmaSideband(0.015, 0.05);
  if (suffix != "0") sigmaCuts->SetLightweight(true);
  sigmaCuts->InitCutHistograms(TString("Sigma0"));

  AliSigma0PhotonMotherCuts *antiSigmaCuts =
      AliSigma0PhotonMotherCuts::DefaultCuts();
  antiSigmaCuts->SetIsMC(isMC);
  antiSigmaCuts->SetSigmaMass(1.192642);
  antiSigmaCuts->SetSigmaMassCut(0.005);
  antiSigmaCuts->SetSigmaSideband(0.015, 0.05);
  if (suffix != "0") antiSigmaCuts->SetLightweight(true);
  antiSigmaCuts->InitCutHistograms(TString("AntiSigma0"));

  AliAnalysisTaskSigma0Run2 *task =
      new AliAnalysisTaskSigma0Run2("AnalysisTaskSigma0");
  if (trigger == "kINT7") {
    task->SetTrigger(AliVEvent::kINT7);
  } else if (trigger == "kHighMultV0") {
    task->SetTrigger(AliVEvent::kHighMultV0);
    task->SetV0Percentile(0.1);
  }
  task->SetV0ReaderName(V0ReaderName.Data());
  task->SetIsHeavyIon(isHeavyIon);
  task->SetIsMC(isMC);
  task->SetV0Cuts(v0Cuts);
  task->SetAntiV0Cuts(antiv0Cuts);
  task->SetSigmaCuts(sigmaCuts);
  task->SetAntiSigmaCuts(antiSigmaCuts);

  if (suffix != "0") task->SetLightweight(true);

  mgr->AddTask(task);

  TString containerName = mgr->GetCommonFileName();
  containerName += ":AnalysisTaskSigma0_";
  containerName += suffix;

  TString name = "histo_" + suffix;
  AliAnalysisDataContainer *cOutputList = mgr->CreateContainer(
      name, TList::Class(), AliAnalysisManager::kOutputContainer,
      containerName.Data());

  mgr->ConnectInput(task, 0, cinput);
  mgr->ConnectOutput(task, 1, cOutputList);

  return task;
}
