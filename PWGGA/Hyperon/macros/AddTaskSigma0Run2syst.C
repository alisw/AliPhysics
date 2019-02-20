AliAnalysisTaskSE *AddTaskSigma0Run2syst(bool isMC = false,
                                         bool isHeavyIon = false,
                                         TString trigger = "kINT7",
                                         const char *cutVariation = "19") {
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
  cutnumberPhoton = "10200008400000002282000000";

  // Standard
  // if (suffix == "00") cutnumberPhoton = "00200008400000002280920000";
  // Lower pT requirement for e (0|40) and gamma (0|10)
  if (suffix == "77") cutnumberPhoton = "00200078400000002280920000";
  if (suffix == "76") cutnumberPhoton = "00200068400000002280920000";
  // Lower Cos PA 0.98 and -1
  if (suffix == "75") cutnumberPhoton = "00200008400000002280820000";
  if (suffix == "74") cutnumberPhoton = "00200008400000002280020000";
  // Smaller Conversion-Radius 0 and 2.8
  if (suffix == "73") cutnumberPhoton = "00000008400000002280920000";
  if (suffix == "72") cutnumberPhoton = "00100008400000002280920000";

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
  v0Cuts->SetPosPID(AliPID::kProton, 2212);
  v0Cuts->SetNegPID(AliPID::kPion, -211);
  AliSigma0V0Cuts *antiv0Cuts = AliSigma0V0Cuts::LambdaCuts();
  antiv0Cuts->SetIsMC(isMC);
  antiv0Cuts->SetPID(-3122);
  antiv0Cuts->SetPosPID(AliPID::kPion, 211);
  antiv0Cuts->SetNegPID(AliPID::kProton, -2212);

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
    v0Cuts->SetTPCRatioFindable(0.35);
    antiv0Cuts->SetTPCclusterMin(0);
    antiv0Cuts->SetTPCRatioFindable(0.35);
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
  if (suffix == "10") {
    v0Cuts->SetPileUpRejectionMode(AliSigma0V0Cuts::OneDaughterCombined);
    antiv0Cuts->SetPileUpRejectionMode(AliSigma0V0Cuts::OneDaughterCombined);
  }
  if (suffix == "11") {
    v0Cuts->SetPileUpRejectionMode(AliSigma0V0Cuts::None);
    antiv0Cuts->SetPileUpRejectionMode(AliSigma0V0Cuts::None);
  }
  if (suffix == "999") {
    v0Cuts->SetCheckCutsMC(true);
    antiv0Cuts->SetCheckCutsMC(true);
    v0Cuts->SetLightweight(false);
    antiv0Cuts->SetLightweight(false);
  }
  if (suffix == "77") {
    v0Cuts->SetK0Rejection(0., 0.);
    antiv0Cuts->SetK0Rejection(0., 0.);
    v0Cuts->SetLambdaSelection(1.115683 - 0.008, 1.115683 + 0.008);
    antiv0Cuts->SetLambdaSelection(1.115683 - 0.008, 1.115683 + 0.008);
    v0Cuts->SetPileUpRejectionMode(AliSigma0V0Cuts::OneDaughterCombined);
    antiv0Cuts->SetPileUpRejectionMode(AliSigma0V0Cuts::OneDaughterCombined);
    v0Cuts->SetArmenterosCut(0.01, 0.12, 0.3, 0.95);
    antiv0Cuts->SetArmenterosCut(0.01, 0.12, 0.3, 0.95);
  }
  if (suffix == "76") {
    v0Cuts->SetK0Rejection(0., 0.);
    antiv0Cuts->SetK0Rejection(0., 0.);
    v0Cuts->SetLambdaSelection(1.115683 - 0.008, 1.115683 + 0.008);
    antiv0Cuts->SetLambdaSelection(1.115683 - 0.008, 1.115683 + 0.008);
    v0Cuts->SetPileUpRejectionMode(AliSigma0V0Cuts::OneDaughterCombined);
    antiv0Cuts->SetPileUpRejectionMode(AliSigma0V0Cuts::OneDaughterCombined);
    v0Cuts->SetArmenterosCut(0.01, 0.12, 0.3, 0.95);
    antiv0Cuts->SetArmenterosCut(0.01, 0.12, 0.3, 0.95);
  }
  if (suffix == "75") {
    v0Cuts->SetK0Rejection(0., 0.);
    antiv0Cuts->SetK0Rejection(0., 0.);
    v0Cuts->SetLambdaSelection(1.115683 - 0.008, 1.115683 + 0.008);
    antiv0Cuts->SetLambdaSelection(1.115683 - 0.008, 1.115683 + 0.008);
    v0Cuts->SetPileUpRejectionMode(AliSigma0V0Cuts::OneDaughterCombined);
    antiv0Cuts->SetPileUpRejectionMode(AliSigma0V0Cuts::OneDaughterCombined);
    v0Cuts->SetArmenterosCut(0.01, 0.12, 0.3, 0.95);
    antiv0Cuts->SetArmenterosCut(0.01, 0.12, 0.3, 0.95);
  }
  if (suffix == "74") {
    v0Cuts->SetK0Rejection(0., 0.);
    antiv0Cuts->SetK0Rejection(0., 0.);
    v0Cuts->SetLambdaSelection(1.115683 - 0.008, 1.115683 + 0.008);
    antiv0Cuts->SetLambdaSelection(1.115683 - 0.008, 1.115683 + 0.008);
    v0Cuts->SetPileUpRejectionMode(AliSigma0V0Cuts::OneDaughterCombined);
    antiv0Cuts->SetPileUpRejectionMode(AliSigma0V0Cuts::OneDaughterCombined);
    v0Cuts->SetArmenterosCut(0.01, 0.12, 0.3, 0.95);
    antiv0Cuts->SetArmenterosCut(0.01, 0.12, 0.3, 0.95);
  }
  if (suffix == "73") {
    v0Cuts->SetK0Rejection(0., 0.);
    antiv0Cuts->SetK0Rejection(0., 0.);
    v0Cuts->SetLambdaSelection(1.115683 - 0.008, 1.115683 + 0.008);
    antiv0Cuts->SetLambdaSelection(1.115683 - 0.008, 1.115683 + 0.008);
    v0Cuts->SetPileUpRejectionMode(AliSigma0V0Cuts::OneDaughterCombined);
    antiv0Cuts->SetPileUpRejectionMode(AliSigma0V0Cuts::OneDaughterCombined);
    v0Cuts->SetArmenterosCut(0.01, 0.12, 0.3, 0.95);
    antiv0Cuts->SetArmenterosCut(0.01, 0.12, 0.3, 0.95);
  }
  if (suffix == "72") {
    v0Cuts->SetK0Rejection(0., 0.);
    antiv0Cuts->SetK0Rejection(0., 0.);
    v0Cuts->SetLambdaSelection(1.115683 - 0.008, 1.115683 + 0.008);
    antiv0Cuts->SetLambdaSelection(1.115683 - 0.008, 1.115683 + 0.008);
    v0Cuts->SetPileUpRejectionMode(AliSigma0V0Cuts::OneDaughterCombined);
    antiv0Cuts->SetPileUpRejectionMode(AliSigma0V0Cuts::OneDaughterCombined);
    v0Cuts->SetArmenterosCut(0.01, 0.12, 0.3, 0.95);
    antiv0Cuts->SetArmenterosCut(0.01, 0.12, 0.3, 0.95);
  }
  AliSigma0V0Cuts *photonV0Cuts = AliSigma0V0Cuts::PhotonCuts();
  photonV0Cuts->SetIsMC(isMC);
  photonV0Cuts->SetPID(22);
  photonV0Cuts->SetPosPID(AliPID::kElectron, 11);
  photonV0Cuts->SetNegPID(AliPID::kElectron, -11);
  if (suffix != "0") photonV0Cuts->SetLightweight(true);
  if (suffix == "12") {
    v0Cuts->SetPileUpRejectionMode(AliSigma0V0Cuts::OneDaughterCombined);
  }
  if (suffix == "13") {
    v0Cuts->SetPileUpRejectionMode(AliSigma0V0Cuts::BothDaughtersCombined);
  }
  if (suffix == "999") {
    photonV0Cuts->SetCheckCutsMC(true);
    photonV0Cuts->SetLightweight(false);
  }
  if (suffix == "77") {
    photonV0Cuts->SetV0OnFlyStatus(true);
    photonV0Cuts->SetArmenterosCut(0, 0.06, -1, 1);
    photonV0Cuts->SetV0RadiusMax(180.f);
    photonV0Cuts->SetV0RadiusMin(5.f);
    photonV0Cuts->SetTPCRatioFindable(0.35f);
    photonV0Cuts->SetV0CosPAMin(0.99);
    photonV0Cuts->SetChi2Max(5);
  }
  if (suffix == "76") {
    photonV0Cuts->SetV0OnFlyStatus(true);
    photonV0Cuts->SetArmenterosCut(0, 0.06, -1, 1);
    photonV0Cuts->SetV0RadiusMax(180.f);
    photonV0Cuts->SetV0RadiusMin(5.f);
    photonV0Cuts->SetTPCRatioFindable(0.35f);
    photonV0Cuts->SetV0CosPAMin(0.99);
    photonV0Cuts->SetChi2Max(5);
  }
  if (suffix == "75") {
    photonV0Cuts->SetV0OnFlyStatus(true);
    photonV0Cuts->SetArmenterosCut(0, 0.06, -1, 1);
    photonV0Cuts->SetV0RadiusMax(180.f);
    photonV0Cuts->SetV0RadiusMin(5.f);
    photonV0Cuts->SetTPCRatioFindable(0.35f);
    photonV0Cuts->SetV0CosPAMin(0.98);
    photonV0Cuts->SetChi2Max(5);
  }
  if (suffix == "74") {
    photonV0Cuts->SetV0OnFlyStatus(true);
    photonV0Cuts->SetArmenterosCut(0, 0.06, -1, 1);
    photonV0Cuts->SetV0RadiusMax(180.f);
    photonV0Cuts->SetV0RadiusMin(5.f);
    photonV0Cuts->SetTPCRatioFindable(0.35f);
    photonV0Cuts->SetV0CosPAMin(-1);
    photonV0Cuts->SetChi2Max(5);
  }
  if (suffix == "73") {
    photonV0Cuts->SetV0OnFlyStatus(true);
    photonV0Cuts->SetArmenterosCut(0, 0.06, -1, 1);
    photonV0Cuts->SetV0RadiusMax(180.f);
    photonV0Cuts->SetV0RadiusMin(0.f);
    photonV0Cuts->SetTPCRatioFindable(0.35f);
    photonV0Cuts->SetV0CosPAMin(0.99);
    photonV0Cuts->SetChi2Max(5);
  }
  if (suffix == "72") {
    photonV0Cuts->SetV0OnFlyStatus(true);
    photonV0Cuts->SetArmenterosCut(0, 0.06, -1, 1);
    photonV0Cuts->SetV0RadiusMax(180.f);
    photonV0Cuts->SetV0RadiusMin(2.8f);
    photonV0Cuts->SetTPCRatioFindable(0.35f);
    photonV0Cuts->SetV0CosPAMin(0.99);
    photonV0Cuts->SetChi2Max(5);
  }
  AliSigma0PhotonMotherCuts *sigmaCuts =
      AliSigma0PhotonMotherCuts::DefaultCuts();
  sigmaCuts->SetIsMC(isMC);
  sigmaCuts->SetPDG(3212, 3122, 22);
  sigmaCuts->SetSigmaMassCut(0.2);
  sigmaCuts->SetLambdaCuts(v0Cuts);
  sigmaCuts->SetV0ReaderName(V0ReaderName.Data());
  if (suffix != "0") {
    sigmaCuts->SetLightweight(true);
  }
  if (suffix == "19") sigmaCuts->SetPhotonMaxPt(1.5);
  if (suffix == "20") sigmaCuts->SetPhotonMaxPt(1.);

  AliSigma0PhotonMotherCuts *antiSigmaCuts =
      AliSigma0PhotonMotherCuts::DefaultCuts();
  antiSigmaCuts->SetIsMC(isMC);
  // Fixed sign of PDG Codes.
  antiSigmaCuts->SetPDG(-3212, -3122, 22);
  antiSigmaCuts->SetSigmaMassCut(0.2);
  antiSigmaCuts->SetLambdaCuts(antiv0Cuts);
  antiSigmaCuts->SetV0ReaderName(V0ReaderName.Data());
  if (suffix != "0") {
    antiSigmaCuts->SetLightweight(true);
    antiSigmaCuts->SetTreeOutput(false);
  }
  if (suffix == "19") antiSigmaCuts->SetPhotonMaxPt(1.5);
  if (suffix == "20") antiSigmaCuts->SetPhotonMaxPt(1.);

  AliSigma0PhotonMotherCuts *sigmaPhotonCuts =
      AliSigma0PhotonMotherCuts::DefaultCuts();
  sigmaPhotonCuts->SetIsMC(isMC);
  sigmaPhotonCuts->SetPDG(3212, 3122, 22);
  sigmaPhotonCuts->SetSigmaMassCut(0.2);
  sigmaPhotonCuts->SetPhotonCuts(photonV0Cuts);
  sigmaPhotonCuts->SetLambdaCuts(v0Cuts);
  if (suffix != "0") {
    sigmaPhotonCuts->SetLightweight(true);
  }
  if (suffix == "19") sigmaPhotonCuts->SetPhotonMaxPt(1.5);
  if (suffix == "20") sigmaPhotonCuts->SetPhotonMaxPt(1.);

  AliSigma0PhotonMotherCuts *antiSigmaPhotonCuts =
      AliSigma0PhotonMotherCuts::DefaultCuts();
  antiSigmaPhotonCuts->SetIsMC(isMC);
  antiSigmaPhotonCuts->SetPDG(-3212, -3122, 22);
  antiSigmaPhotonCuts->SetSigmaMassCut(0.2);
  antiSigmaPhotonCuts->SetPhotonCuts(photonV0Cuts);
  antiSigmaPhotonCuts->SetLambdaCuts(antiv0Cuts);
  if (suffix != "0") {
    antiSigmaPhotonCuts->SetLightweight(true);
  }
  if (suffix == "19") antiSigmaPhotonCuts->SetPhotonMaxPt(1.5);
  if (suffix == "20") antiSigmaPhotonCuts->SetPhotonMaxPt(1.);

  AliAnalysisTaskSigma0Run2 *task =
      new AliAnalysisTaskSigma0Run2("AnalysisTaskSigma0");
  if (trigger == "kINT7") {
    task->SetTrigger(AliVEvent::kINT7);
    task->SelectCollisionCandidates(AliVEvent::kINT7);
  } else if (trigger == "kHighMultV0") {
    task->SetTrigger(AliVEvent::kHighMultV0);
    task->SelectCollisionCandidates(AliVEvent::kHighMultV0);
  }
  task->SetV0ReaderName(V0ReaderName.Data());
  task->SetIsHeavyIon(isHeavyIon);
  task->SetIsMC(isMC);
  task->SetV0Cuts(v0Cuts);
  task->SetAntiV0Cuts(antiv0Cuts);
  task->SetPhotonV0Cuts(photonV0Cuts);
  task->SetSigmaCuts(sigmaCuts);
  task->SetAntiSigmaCuts(antiSigmaCuts);
  task->SetSigmaPhotonCuts(sigmaPhotonCuts);
  task->SetAntiSigmaPhotonCuts(antiSigmaPhotonCuts);

  if (suffix != "0") task->SetLightweight(true);

  mgr->AddTask(task);

  TString containerName = mgr->GetCommonFileName();
  containerName += ":AnalysisTaskSigma0_";
  if (trigger == "kHighMultV0") containerName += "HighMultV0_";
  containerName += suffix;

  TString name = "histo_";
  if (trigger == "kHighMultV0") name += "HighMultV0_";
  name += suffix;
  AliAnalysisDataContainer *cOutputList = mgr->CreateContainer(
      name, TList::Class(), AliAnalysisManager::kOutputContainer,
      containerName.Data());

  mgr->ConnectInput(task, 0, cinput);
  mgr->ConnectOutput(task, 1, cOutputList);

  return task;
}
