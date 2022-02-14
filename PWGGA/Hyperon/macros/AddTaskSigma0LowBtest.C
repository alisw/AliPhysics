AliAnalysisTaskSE *AddTaskSigma0LowBtest(bool isMC = false, TString trigger =
                                             "kINT7",
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
  cutnumberPhoton = "00200008400000002280f20090";
  if (suffix == "1") {
    cutnumberPhoton = "00200088400000002280f20090";
  }
  TString cutnumberEvent = "00000000";
  TString periodNameV0Reader = "";
  Bool_t enableV0findingEffi = kFALSE;
  Bool_t fillHistos = kTRUE;
  Bool_t runLightOutput = kFALSE;
  if (suffix != "0" && suffix != "999" && suffix != "3" && suffix != "4") {
    runLightOutput = kTRUE;
    fillHistos = kFALSE;
  }

  //========= Add V0 Reader to  ANALYSIS manager if not yet existent =====
  TString V0ReaderName = Form("V0ReaderV1_%s_%s", cutnumberEvent.Data(),
                              cutnumberPhoton.Data());
  AliConvEventCuts *fEventCuts = NULL;

  if (!(AliV0ReaderV1 *) mgr->GetTask(V0ReaderName.Data())) {
    AliV0ReaderV1 *fV0ReaderV1 = new AliV0ReaderV1(V0ReaderName.Data());
    if (periodNameV0Reader.CompareTo("") != 0)
      fV0ReaderV1->SetPeriodName(periodNameV0Reader);
    fV0ReaderV1->SetUseOwnXYZCalculation(kTRUE);
    fV0ReaderV1->SetCreateAODs(kFALSE);  // AOD Output
    fV0ReaderV1->SetUseAODConversionPhoton(kTRUE);
    fV0ReaderV1->SetProduceV0FindingEfficiency(enableV0findingEffi);

    if (!mgr) {
      Error("AddTask_V0ReaderV1", "No analysis manager found.");
      return NULL;
    }

    if (cutnumberEvent != "") {
      fEventCuts = new AliConvEventCuts(cutnumberEvent.Data(),
                                        cutnumberEvent.Data());
      fEventCuts->SetPreSelectionCutFlag(kTRUE);
      fEventCuts->SetV0ReaderName(V0ReaderName);
      fEventCuts->SetLightOutput(runLightOutput);
      fEventCuts->SetFillCutHistograms("", fillHistos);
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
      fCuts->SetIsHeavyIon(false);
      fCuts->SetV0ReaderName(V0ReaderName);
      fCuts->SetLightOutput(runLightOutput);
      fCuts->SetFillCutHistograms("", fillHistos);
      if (fCuts->InitializeCutsFromCutString(cutnumberPhoton.Data())) {
        fV0ReaderV1->SetConversionCuts(fCuts);
      }
    }
    fV0ReaderV1->Init();
    AliLog::SetGlobalLogLevel(AliLog::kFatal);

    // connect input V0Reader
    mgr->AddTask(fV0ReaderV1);
    mgr->ConnectInput(fV0ReaderV1, 0, cinput);
  }

  //========= Init subtasks and start analyis ============================
  // Track Cuts

  AliSigma0V0Cuts *v0Cuts = AliSigma0V0Cuts::LambdaCuts();
  v0Cuts->SetIsMC(isMC);
  v0Cuts->SetPID(3122);
  v0Cuts->SetPosPID(AliPID::kProton, 2212);
  v0Cuts->SetNegPID(AliPID::kPion, -211);

  AliSigma0V0Cuts *antiv0Cuts = AliSigma0V0Cuts::LambdaCuts();
  antiv0Cuts->SetIsMC(isMC);
  antiv0Cuts->SetPID(-3122);
  antiv0Cuts->SetPosPID(AliPID::kPion, 211);
  antiv0Cuts->SetNegPID(AliPID::kProton, -2212);

  if (suffix != "0" || suffix != "1") {
    v0Cuts->SetLightweight(true);
    antiv0Cuts->SetLightweight(true);
  }
  if (suffix == "999") {
    v0Cuts->SetCheckCutsMC(true);
    antiv0Cuts->SetCheckCutsMC(true);
    v0Cuts->SetLightweight(false);
    antiv0Cuts->SetLightweight(false);
  }

  AliSigma0PhotonMotherCuts *sigmaCuts =
      AliSigma0PhotonMotherCuts::DefaultCuts();
  sigmaCuts->SetIsMC(isMC);
  sigmaCuts->SetPDG(3212, 3122, 22);
  sigmaCuts->SetLambdaCuts(v0Cuts);
  sigmaCuts->SetV0ReaderName(V0ReaderName.Data());
  if (suffix != "0" && suffix != "1" && suffix != "999") {
    sigmaCuts->SetLightweight(true);
  }

  AliSigma0PhotonMotherCuts *antiSigmaCuts =
      AliSigma0PhotonMotherCuts::DefaultCuts();
  antiSigmaCuts->SetIsMC(isMC);
  antiSigmaCuts->SetPDG(-3212, -3122, 22);
  antiSigmaCuts->SetLambdaCuts(antiv0Cuts);
  antiSigmaCuts->SetV0ReaderName(V0ReaderName.Data());
  if (suffix != "0" && suffix != "1" && suffix != "999") {
    antiSigmaCuts->SetLightweight(true);
  }

  if (trigger == "kINT7") {
    sigmaCuts->SetMultiplicityMode(AliVEvent::kINT7);
    antiSigmaCuts->SetMultiplicityMode(AliVEvent::kINT7);
  } else if (trigger == "kHighMultV0") {
    sigmaCuts->SetMultiplicityMode(AliVEvent::kHighMultV0);
    antiSigmaCuts->SetMultiplicityMode(AliVEvent::kHighMultV0);
  } else if (trigger == "AliVEvent::kMB") {
    sigmaCuts->SetMultiplicityMode(AliVEvent::kINT7);
    antiSigmaCuts->SetMultiplicityMode(AliVEvent::kINT7);
  }

  AliAnalysisTaskSigma0Run2 *task = new AliAnalysisTaskSigma0Run2(
      "AnalysisTaskSigma0Run2");
  if (trigger == "kINT7") {
    task->SetTrigger(AliVEvent::kINT7);
    task->SetMultiplicityMode(AliVEvent::kINT7);
    task->SelectCollisionCandidates(AliVEvent::kINT7);
  } else if (trigger == "kHighMultV0") {
    if (isMC) {
      task->SetTrigger(AliVEvent::kINT7);
      task->SelectCollisionCandidates(AliVEvent::kINT7);
      task->SetMultiplicityMode(AliVEvent::kHighMultV0);
    } else {
      task->SetTrigger(AliVEvent::kHighMultV0);
      task->SelectCollisionCandidates(AliVEvent::kHighMultV0);
      task->SetMultiplicityMode(AliVEvent::kHighMultV0);
    }
  } else if (trigger == "AliVEvent::kMB") {
    task->SetTrigger(AliVEvent::kMB);
    task->SelectCollisionCandidates(AliVEvent::kMB);
    task->SetMultiplicityMode(AliVEvent::kINT7);
  }
  task->SetV0ReaderName(V0ReaderName.Data());
  task->SetIsRun1(false);
  task->SetIsHeavyIon(false);
  task->SetIsMC(isMC);
  task->SetPhotonDCArCut(true);
  task->SetV0Cuts(v0Cuts);
  task->SetAntiV0Cuts(antiv0Cuts);
  task->SetSigmaCuts(sigmaCuts);
  task->SetAntiSigmaCuts(antiSigmaCuts);

  if (suffix != "0" && suffix != "1" && suffix != "999") {
    task->SetLightweight(true);
  }

  mgr->AddTask(task);

  TString containerName = mgr->GetCommonFileName();
  containerName += ":Sigma0_";
  if (trigger == "kHighMultV0")
    containerName += "HighMultV0_";
  containerName += suffix;

  TString name = "histo_";
  if (trigger == "kHighMultV0")
    name += "HighMultV0_";
  name += suffix;
  AliAnalysisDataContainer *cOutputList = mgr->CreateContainer(
      name, TList::Class(), AliAnalysisManager::kOutputContainer,
      containerName.Data());

  mgr->ConnectInput(task, 0, cinput);
  mgr->ConnectOutput(task, 1, cOutputList);

  return task;
}
