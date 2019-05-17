AliAnalysisTaskSE *AddTaskSigma0FemtoAOD(bool isMC = false, bool MomRes = false,
                                         bool fullBlastQA = false,
                                         TString trigger = "kINT7",
                                         const char *cutVariation = "0") {
  TString suffix = TString::Format("%s", cutVariation);

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
  cutnumberPhoton = "00200008400000002280920000";
  TString cutnumberEvent = "00000000";
  TString periodNameV0Reader = "";
  Bool_t enableV0findingEffi = kFALSE;
  Bool_t fillHistos = kTRUE;
  Bool_t runLightOutput = kFALSE;
  if (suffix != "0" && suffix != "999") {
    runLightOutput = kTRUE;
    fillHistos = kFALSE;
  }
  if (suffix == "21") {
    // eta < 0.8
    cutnumberPhoton = "0d200008400000002280920000";
  } else if (suffix == "22") {
    // single pT > 0, gammapT > 0
    cutnumberPhoton = "00200078400000002280920000";
  } else if (suffix == "23") {
    // single pT > 0.050, gammapT > 0.150
    cutnumberPhoton = "002000a8400000002280920000";
  } else if (suffix == "24") {
    // TPC cluster, findable > 0.6
    cutnumberPhoton = "00200009400000002280920000";
  } else if (suffix == "25") {
    // TPC PID -10,10
    cutnumberPhoton = "00200008000000002280920000";
  } else if (suffix == "26") {
    // TPC PID -3,3
    cutnumberPhoton = "00200008a00000002280920000";
  } else if (suffix == "27") {
    // 1-D Qt cut, qt < 0.1
    cutnumberPhoton = "00200008400000001280920000";
  } else if (suffix == "28") {
    // 2-D Qt cut, qt < 0.05
    cutnumberPhoton = "00200008400000003280920000";
  } else if (suffix == "29") {
    // psiPair < 0.2, 1-D
    cutnumberPhoton = "00200008400000002240920000";
  } else if (suffix == "30") {
    // psiPair < 0.1, 2-D
    cutnumberPhoton = "00200008400000002250920000";
  } else if (suffix == "31") {
    // cosPA < 0.98
    cutnumberPhoton = "00200008400000002280820000";
  } else if (suffix == "32") {
    // cosPA < 0.995
    cutnumberPhoton = "00200008400000002280a20000";
  }

  //========= Add V0 Reader to  ANALYSIS manager if not yet existent =====
  TString V0ReaderName = TString::Format(
      "V0ReaderV1_%s_%s", cutnumberEvent.Data(), cutnumberPhoton.Data());
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
      return NULL;
    }

    if (cutnumberEvent != "") {
      fEventCuts =
          new AliConvEventCuts(cutnumberEvent.Data(), cutnumberEvent.Data());
      fEventCuts->SetPreSelectionCutFlag(kTRUE);
      fEventCuts->SetV0ReaderName(V0ReaderName);
      fEventCuts->SetLightOutput(runLightOutput);
      if (periodNameV0Reader.CompareTo("") != 0)
        fEventCuts->SetPeriodEnum(periodNameV0Reader);
      fV0ReaderV1->SetEventCuts(fEventCuts);
      fEventCuts->SetFillCutHistograms("", kFALSE);
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
  // Event Cuts
  AliFemtoDreamEventCuts *evtCuts = AliFemtoDreamEventCuts::StandardCutsRun2();
  evtCuts->CleanUpMult(false, false, false, true);

  // Track Cuts
  AliFemtoDreamTrackCuts *TrackCuts =
      AliFemtoDreamTrackCuts::PrimProtonCuts(isMC, true, false, false);
  TrackCuts->SetFilterBit(128);
  TrackCuts->SetCutCharge(1);

  AliFemtoDreamTrackCuts *AntiTrackCuts =
      AliFemtoDreamTrackCuts::PrimProtonCuts(isMC, true, false, false);
  AntiTrackCuts->SetFilterBit(128);
  AntiTrackCuts->SetCutCharge(-1);

  if (suffix != "0" && suffix != "999") {
    TrackCuts->SetMinimalBooking(true);
    AntiTrackCuts->SetMinimalBooking(true);
  }
  if (suffix == "1") {
    TrackCuts->SetPtRange(0.4, 4.05);
    AntiTrackCuts->SetPtRange(0.4, 4.05);
  } else if (suffix == "2") {
    TrackCuts->SetPtRange(0.6, 4.05);
    AntiTrackCuts->SetPtRange(0.6, 4.05);
  } else if (suffix == "3") {
    TrackCuts->SetEtaRange(-0.7, 0.7);
    AntiTrackCuts->SetEtaRange(-0.7, 0.7);
  } else if (suffix == "4") {
    TrackCuts->SetEtaRange(-0.9, 0.9);
    AntiTrackCuts->SetEtaRange(-0.9, 0.9);
  } else if (suffix == "5") {
    TrackCuts->SetPID(AliPID::kProton, 0.75, 2);
    AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 2);
  } else if (suffix == "6") {
    TrackCuts->SetPID(AliPID::kProton, 0.75, 5);
    AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 5);
  } else if (suffix == "7") {
    TrackCuts->SetNClsTPC(70);
    AntiTrackCuts->SetNClsTPC(70);
  } else if (suffix == "8") {
    TrackCuts->SetNClsTPC(90);
    AntiTrackCuts->SetNClsTPC(90);
  } else if (suffix == "9") {
    TrackCuts->SetFilterBit(96);
    AntiTrackCuts->SetFilterBit(96);
  }

  // Lambda Cuts
  AliFemtoDreamv0Cuts *v0Cuts =
      AliFemtoDreamv0Cuts::LambdaSigma0Cuts(isMC, false, false);
  AliFemtoDreamTrackCuts *Posv0Daug =
      AliFemtoDreamTrackCuts::DecayProtonCuts(isMC, false, false);
  Posv0Daug->SetEtaRange(-0.9, 0.9);
  AliFemtoDreamTrackCuts *Negv0Daug =
      AliFemtoDreamTrackCuts::DecayPionCuts(isMC, false, false);
  Negv0Daug->SetEtaRange(-0.9, 0.9);

  AliFemtoDreamv0Cuts *antiv0Cuts =
      AliFemtoDreamv0Cuts::LambdaSigma0Cuts(isMC, false, false);
  AliFemtoDreamTrackCuts *PosAntiv0Daug =
      AliFemtoDreamTrackCuts::DecayPionCuts(isMC, false, false);
  PosAntiv0Daug->SetCutCharge(1);
  PosAntiv0Daug->SetEtaRange(-0.9, 0.9);
  AliFemtoDreamTrackCuts *NegAntiv0Daug =
      AliFemtoDreamTrackCuts::DecayProtonCuts(isMC, false, false);
  NegAntiv0Daug->SetCutCharge(-1);
  NegAntiv0Daug->SetEtaRange(-0.9, 0.9);

  if (suffix != "0") {
    v0Cuts->SetMinimalBooking(true);
    antiv0Cuts->SetMinimalBooking(true);
  }
  if (suffix == "10") {
    v0Cuts->SetPtRange(0.24, 999.9);
    antiv0Cuts->SetPtRange(0.24, 999.9);
  } else if (suffix == "11") {
    v0Cuts->SetPtRange(0.36, 999.9);
    antiv0Cuts->SetPtRange(0.36, 999.9);
  } else if (suffix == "12") {
    v0Cuts->SetCutCPA(0.995);
    antiv0Cuts->SetCutCPA(0.995);
  } else if (suffix == "13") {
    v0Cuts->SetCutCPA(0.98);
    antiv0Cuts->SetCutCPA(0.998);
  } else if (suffix == "14") {
    Posv0Daug->SetPID(AliPID::kProton, 999.9, 3);
    Negv0Daug->SetPID(AliPID::kPion, 999.9, 3);
    PosAntiv0Daug->SetPID(AliPID::kPion, 999.9, 3);
    NegAntiv0Daug->SetPID(AliPID::kProton, 999.9, 3);
  } else if (suffix == "15") {
    Posv0Daug->SetPID(AliPID::kProton, 999.9, 6);
    Negv0Daug->SetPID(AliPID::kPion, 999.9, 6);
    PosAntiv0Daug->SetPID(AliPID::kPion, 999.9, 6);
    NegAntiv0Daug->SetPID(AliPID::kProton, 999.9, 6);
  } else if (suffix == "16") {
    v0Cuts->SetArmenterosCut(0., 0.15, 0.2, 1);
    antiv0Cuts->SetArmenterosCut(0, 0.15, 0.2, 1);
  } else if (suffix == "17") {
    Posv0Daug->SetNClsTPC(80);
    Negv0Daug->SetNClsTPC(80);
    PosAntiv0Daug->SetNClsTPC(80);
    NegAntiv0Daug->SetNClsTPC(80);
  } else if (suffix == "18") {
    Posv0Daug->SetNClsTPC(60);
    Negv0Daug->SetNClsTPC(60);
    PosAntiv0Daug->SetNClsTPC(60);
    NegAntiv0Daug->SetNClsTPC(60);
  } else if (suffix == "19") {
    Posv0Daug->SetEtaRange(-0.8, 0.8);
    Negv0Daug->SetEtaRange(-0.8, 0.8);
    PosAntiv0Daug->SetEtaRange(-0.8, 0.8);
    NegAntiv0Daug->SetEtaRange(-0.8, 0.8);
  } else if (suffix == "20") {
    v0Cuts->SetCutInvMass(0.008);
    antiv0Cuts->SetCutInvMass(0.008);
  }

  v0Cuts->SetPosDaugterTrackCuts(Posv0Daug);
  v0Cuts->SetNegDaugterTrackCuts(Negv0Daug);
  v0Cuts->SetPDGCodePosDaug(2212);  // Proton
  v0Cuts->SetPDGCodeNegDaug(211);   // Pion
  v0Cuts->SetPDGCodev0(3122);       // Lambda
  antiv0Cuts->SetPosDaugterTrackCuts(PosAntiv0Daug);
  antiv0Cuts->SetNegDaugterTrackCuts(NegAntiv0Daug);
  antiv0Cuts->SetPDGCodePosDaug(211);   // Pion
  antiv0Cuts->SetPDGCodeNegDaug(2212);  // Proton
  antiv0Cuts->SetPDGCodev0(-3122);      // Lambda

  AliSigma0PhotonCuts *photon = AliSigma0PhotonCuts::PhotonCuts();
  photon->SetV0ReaderName(V0ReaderName.Data());

  AliSigma0AODPhotonMotherCuts *sigmaCuts =
      AliSigma0AODPhotonMotherCuts::DefaultCuts();
  sigmaCuts->SetIsMC(isMC);
  sigmaCuts->SetPDG(3212, 3122, 22);
  if (suffix != "0" && suffix != "999") {
    sigmaCuts->SetLightweight(true);
  }

  AliSigma0AODPhotonMotherCuts *antiSigmaCuts =
      AliSigma0AODPhotonMotherCuts::DefaultCuts();
  antiSigmaCuts->SetIsMC(isMC);
  antiSigmaCuts->SetPDG(-3212, -3122, 22);
  if (suffix != "0" && suffix != "999") {
    antiSigmaCuts->SetLightweight(true);
  }

  // vary the sidebands
  if (suffix == "33") {
    sigmaCuts->SetSigmaSideband(0.01, 0.075);
    antiSigmaCuts->SetSigmaSideband(0.01, 0.075);
  } else if (suffix == "34") {
    sigmaCuts->SetSigmaSideband(0.01, 0.1);
    antiSigmaCuts->SetSigmaSideband(0.01, 0.1);
  } else if (suffix == "35") {
    sigmaCuts->SetSigmaSideband(0.025, 0.05);
    antiSigmaCuts->SetSigmaSideband(0.025, 0.05);
  } else if (suffix == "36") {
    sigmaCuts->SetSigmaSideband(0.025, 0.075);
    antiSigmaCuts->SetSigmaSideband(0.025, 0.075);
  } else if (suffix == "37") {
    sigmaCuts->SetSigmaSideband(0.025, 0.1);
    antiSigmaCuts->SetSigmaSideband(0.025, 0.1);
  } else if (suffix == "38") {
    sigmaCuts->SetSigmaSideband(0.005, 0.025);
    antiSigmaCuts->SetSigmaSideband(0.005, 0.025);
  } else if (suffix == "39") {
    sigmaCuts->SetSigmaSideband(0.005, 0.05);
    antiSigmaCuts->SetSigmaSideband(0.005, 0.05);
  } else if (suffix == "40") {
    sigmaCuts->SetSigmaSideband(0.005, 0.075);
    antiSigmaCuts->SetSigmaSideband(0.005, 0.075);
  } else if (suffix == "41") {
    sigmaCuts->SetSigmaSideband(0.005, 0.1);
    antiSigmaCuts->SetSigmaSideband(0.005, 0.1);
  }

  // Femto Collection
  std::vector<int> PDGParticles;
  PDGParticles.push_back(2212);
  PDGParticles.push_back(2212);
  PDGParticles.push_back(3212);
  PDGParticles.push_back(3212);
  PDGParticles.push_back(3212);
  PDGParticles.push_back(3212);
  PDGParticles.push_back(3212);
  PDGParticles.push_back(3212);

  std::vector<float> ZVtxBins;
  ZVtxBins.push_back(-10);
  ZVtxBins.push_back(-8);
  ZVtxBins.push_back(-6);
  ZVtxBins.push_back(-4);
  ZVtxBins.push_back(-2);
  ZVtxBins.push_back(0);
  ZVtxBins.push_back(2);
  ZVtxBins.push_back(4);
  ZVtxBins.push_back(6);
  ZVtxBins.push_back(8);
  ZVtxBins.push_back(10);

  std::vector<int> NBins;
  std::vector<float> kMin;
  std::vector<float> kMax;
  std::vector<int> pairQA;
  std::vector<bool> closeRejection;
  const int nPairs = 36;
  for (int i = 0; i < nPairs; ++i) {
    pairQA.push_back(0);
    closeRejection.push_back(false);
    if (suffix == "0") {
      NBins.push_back(600);
      kMin.push_back(0.);
      kMax.push_back(3.);
    } else {
      NBins.push_back(200);
      kMin.push_back(0.);
      kMax.push_back(1.);
    }
  }

  NBins[0] = 250;  // pp
  NBins[8] = 250;  // barp barp

  closeRejection[0] = true;  // pp
  closeRejection[8] = true;  // barp barp

  // do extended QA for the pairs in default mode
  if (suffix == "0") {
    NBins[0] = 750;  // pp
    NBins[8] = 750;  // barp barp

    pairQA[0] = 11;   // pp
    pairQA[2] = 14;   // pSigma
    pairQA[8] = 11;   // barp barp
    pairQA[10] = 14;  // barp bSigma
  }

  AliFemtoDreamCollConfig *config =
      new AliFemtoDreamCollConfig("Femto", "Femto");
  std::vector<int> MultBins;
  if (trigger == "kHighMultV0") {
    std::vector<int> MultBins;
    MultBins.push_back(0);
    MultBins.push_back(4);
    MultBins.push_back(8);
    MultBins.push_back(12);
    MultBins.push_back(16);
    MultBins.push_back(20);
    MultBins.push_back(24);
    MultBins.push_back(28);
    MultBins.push_back(32);
    MultBins.push_back(36);
    MultBins.push_back(40);
    MultBins.push_back(44);
    MultBins.push_back(48);
    MultBins.push_back(52);
    MultBins.push_back(56);
    MultBins.push_back(60);
    MultBins.push_back(64);
    MultBins.push_back(68);
    MultBins.push_back(72);
    MultBins.push_back(76);
    MultBins.push_back(80);
    MultBins.push_back(84);
    MultBins.push_back(88);
    MultBins.push_back(92);
    MultBins.push_back(96);
    MultBins.push_back(100);
    config->SetMultBins(MultBins);
  } else {
    std::vector<int> MultBins;
    MultBins.push_back(0);
    MultBins.push_back(4);
    MultBins.push_back(8);
    MultBins.push_back(12);
    MultBins.push_back(16);
    MultBins.push_back(20);
    MultBins.push_back(24);
    MultBins.push_back(28);
    MultBins.push_back(32);
    MultBins.push_back(36);
    MultBins.push_back(40);
    MultBins.push_back(60);
    MultBins.push_back(80);
    config->SetMultBins(MultBins);
  }
  config->SetMultBinning(true);

  config->SetExtendedQAPairs(pairQA);
  config->SetdPhidEtaPlotsSmallK(false);
  config->SetZBins(ZVtxBins);
  if (MomRes && isMC) {
    config->SetMomentumResolution(true);
  }

  if (trigger == "kHighMultV0") {
    config->SetDeltaEtaMax(0.01);
    config->SetDeltaPhiMax(0.01);
    config->SetClosePairRejection(closeRejection);
  }

  if (suffix == "0" && fullBlastQA) {
    config->SetPhiEtaBinnign(true);
    config->SetkTBinning(true);
    config->SetmTBinning(true);
    config->SetPtQA(true);
  }
  config->SetdPhidEtaPlots(false);
  config->SetPDGCodes(PDGParticles);
  config->SetNBinsHist(NBins);
  config->SetMinKRel(kMin);
  config->SetMaxKRel(kMax);
  config->SetMixingDepth(10);
  config->SetUseEventMixing(true);
  config->SetMultiplicityEstimator(AliFemtoDreamEvent::kRef08);
  if (suffix != "0") {
    config->SetMinimalBookingME(true);
  }

  AliAnalysisTaskAODSigma0Femto *task =
      new AliAnalysisTaskAODSigma0Femto("AliAnalysisTaskAODSigma0Femto", isMC);
  if (trigger == "kINT7") {
    task->SetTrigger(AliVEvent::kINT7);
    task->SelectCollisionCandidates(AliVEvent::kINT7);
  } else if (trigger == "kHighMultV0") {
    if (isMC) {
      task->SetTrigger(AliVEvent::kINT7);
      task->SelectCollisionCandidates(AliVEvent::kINT7);
    } else {
      task->SetTrigger(AliVEvent::kHighMultV0);
      task->SelectCollisionCandidates(AliVEvent::kHighMultV0);
    }
  } else if (trigger == "AliVEvent::kMB") {
    task->SetTrigger(AliVEvent::kMB);
    task->SelectCollisionCandidates(AliVEvent::kMB);
  }
  task->SetEventCuts(evtCuts);
  task->SetV0ReaderName(V0ReaderName.Data());
  task->SetProtonCuts(TrackCuts);
  task->SetAntiProtonCuts(AntiTrackCuts);
  task->SetV0Cuts(v0Cuts);
  task->SetAntiV0Cuts(antiv0Cuts);
  task->SetPhotonCuts(photon);
  task->SetSigmaCuts(sigmaCuts);
  task->SetAntiSigmaCuts(antiSigmaCuts);
  task->SetCollectionConfig(config);

  if (suffix != "0" && suffix != "999") {
    task->SetLightweight(true);
  }

  mgr->AddTask(task);

  TString addon = "";
  if (trigger == "kINT7") {
    addon += "MB";
  } else if (trigger == "kHighMultV0") {
    addon += "HM";
  }

  TString file = AliAnalysisManager::GetCommonFileName();

  mgr->ConnectInput(task, 0, cinput);

  TString QAName = Form("%sQA%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *coutputQA = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      QAName.Data(), TList::Class(), AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), QAName.Data()));
  mgr->ConnectOutput(task, 1, coutputQA);

  TString EvtCutsName = Form("%sEvtCuts%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *coutputEvtCuts = mgr->CreateContainer(
      EvtCutsName.Data(), TList::Class(), AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), EvtCutsName.Data()));
  mgr->ConnectOutput(task, 2, coutputEvtCuts);

  TString TrackCutsName = Form("%sTrackCuts%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *couputTrkCuts =
      mgr->CreateContainer(TrackCutsName.Data(), TList::Class(),
                           AliAnalysisManager::kOutputContainer,
                           Form("%s:%s", file.Data(), TrackCutsName.Data()));
  mgr->ConnectOutput(task, 3, couputTrkCuts);

  TString AntiTrackCutsName =
      Form("%sAntiTrackCuts%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *coutputAntiTrkCuts = mgr->CreateContainer(
      AntiTrackCutsName.Data(), TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), AntiTrackCutsName.Data()));
  mgr->ConnectOutput(task, 4, coutputAntiTrkCuts);

  TString v0CutsName = Form("%sv0Cuts%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *coutputv0Cuts = mgr->CreateContainer(
      v0CutsName.Data(), TList::Class(), AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), v0CutsName.Data()));
  mgr->ConnectOutput(task, 5, coutputv0Cuts);

  TString Antiv0CutsName = Form("%sAntiv0Cuts%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *coutputAntiv0Cuts =
      mgr->CreateContainer(Antiv0CutsName.Data(), TList::Class(),
                           AliAnalysisManager::kOutputContainer,
                           Form("%s:%s", file.Data(), Antiv0CutsName.Data()));
  mgr->ConnectOutput(task, 6, coutputAntiv0Cuts);

  TString Photonv0CutsName =
      Form("%sPhotonCuts%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *coutputPhotonv0Cuts =
      mgr->CreateContainer(Photonv0CutsName.Data(), TList::Class(),
                           AliAnalysisManager::kOutputContainer,
                           Form("%s:%s", file.Data(), Photonv0CutsName.Data()));
  mgr->ConnectOutput(task, 7, coutputPhotonv0Cuts);

  TString SigmaCutsName = Form("%sSigma0Cuts%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *coutputSigmaCuts =
      mgr->CreateContainer(SigmaCutsName.Data(), TList::Class(),
                           AliAnalysisManager::kOutputContainer,
                           Form("%s:%s", file.Data(), SigmaCutsName.Data()));
  mgr->ConnectOutput(task, 8, coutputSigmaCuts);

  TString AntiSigmaCutsName =
      Form("%sAntiSigma0Cuts%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *coutputAntiSigmaCuts = mgr->CreateContainer(
      AntiSigmaCutsName.Data(), TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), AntiSigmaCutsName.Data()));
  mgr->ConnectOutput(task, 9, coutputAntiSigmaCuts);

  AliAnalysisDataContainer *coutputResults;
  TString ResultsName = Form("%sResults%s", addon.Data(), suffix.Data());
  coutputResults = mgr->CreateContainer(
      ResultsName.Data(), TList::Class(), AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), ResultsName.Data()));
  mgr->ConnectOutput(task, 10, coutputResults);

  AliAnalysisDataContainer *coutputResultQA;
  TString ResultQAName = Form("%sResultQA%s", addon.Data(), suffix.Data());
  coutputResultQA = mgr->CreateContainer(
      ResultQAName.Data(), TList::Class(), AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), ResultQAName.Data()));
  mgr->ConnectOutput(task, 11, coutputResultQA);

  if (isMC) {
    TString TrkCutsMCName =
        Form("%sTrackCutsMC%s", addon.Data(), suffix.Data());
    AliAnalysisDataContainer *coutputTrkCutsMC =
        mgr->CreateContainer(TrkCutsMCName.Data(), TList::Class(),
                             AliAnalysisManager::kOutputContainer,
                             Form("%s:%s", file.Data(), TrkCutsMCName.Data()));
    mgr->ConnectOutput(task, 12, coutputTrkCutsMC);

    TString AntiTrkCutsMCName =
        Form("%sAntiTrackCutsMC%s", addon.Data(), suffix.Data());
    AliAnalysisDataContainer *coutputAntiTrkCutsMC = mgr->CreateContainer(
        AntiTrkCutsMCName.Data(), TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), AntiTrkCutsMCName.Data()));
    mgr->ConnectOutput(task, 13, coutputAntiTrkCutsMC);

    TString V0CutsMCName = Form("%sv0CutsMC%s", addon.Data(), suffix.Data());
    AliAnalysisDataContainer *coutputV0CutsMC =
        mgr->CreateContainer(V0CutsMCName.Data(), TList::Class(),
                             AliAnalysisManager::kOutputContainer,
                             Form("%s:%s", file.Data(), V0CutsMCName.Data()));
    mgr->ConnectOutput(task, 14, coutputV0CutsMC);

    TString AntiV0CutsMCName =
        Form("%sAntiv0CutsMC%s", addon.Data(), suffix.Data());
    AliAnalysisDataContainer *coutputAntiV0CutsMC = mgr->CreateContainer(
        AntiTrkCutsMCName.Data(), TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), AntiV0CutsMCName.Data()));
    mgr->ConnectOutput(task, 15, coutputAntiV0CutsMC);
  }

  return task;
}
