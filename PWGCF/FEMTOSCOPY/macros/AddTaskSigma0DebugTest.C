#include <vector>

AliAnalysisTaskSE *AddTaskSigma0DebugTest(bool isMC = false,
                                          bool MomRes = false,
                                          bool fullBlastQA = false,
                                          TString trigger = "kINT7",
                                          const char *cutVariation = "0") {
  TString suffix = TString::Format("%s", cutVariation);

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskSigma0DebugTest()", "No analysis manager found.");
    return 0x0;
  }

  // ================== GetInputEventHandler =============================
  AliVEventHandler *inputHandler = mgr->GetInputEventHandler();
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();

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

  if (suffix == "12") {
    v0Cuts->SetCutInvMass(0.0055);
    antiv0Cuts->SetCutInvMass(0.0055);
  } else if (suffix == "13") {
    v0Cuts->SetCutInvMass(0.005);
    antiv0Cuts->SetCutInvMass(0.005);
  } else if (suffix == "14") {
    v0Cuts->SetCutInvMass(0.0045);
    antiv0Cuts->SetCutInvMass(0.0045);
  } else if (suffix == "15") {
    v0Cuts->SetCutInvMass(0.004);
    antiv0Cuts->SetCutInvMass(0.004);
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
  if (suffix == "1") {
    photon->SetDCArMax(1);
  } else if (suffix == "2") {
    photon->SetDCArMax(0.95);
  } else if (suffix == "3") {
    photon->SetDCArMax(0.9);
  } else if (suffix == "4") {
    photon->SetDCArMax(0.85);
  } else if (suffix == "5") {
    photon->SetDCArMax(0.8);
  } else if (suffix == "6") {
    photon->SetDCArMax(0.75);
  } else if (suffix == "7") {
    photon->SetDCArMax(0.7);
  } else if (suffix == "8") {
    photon->SetDCArMax(0.65);
  } else if (suffix == "9") {
    photon->SetDCArMax(0.6);
  } else if (suffix == "10") {
    photon->SetDCArMax(0.55);
  } else if (suffix == "11") {
    photon->SetDCArMax(0.5);
  }

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

  if (suffix == "16") {
    sigmaCuts->SetSigmaMassCut(0.00225);
    antiSigmaCuts->SetSigmaMassCut(0.00225);
  } else if (suffix == "17") {
    sigmaCuts->SetSigmaMassCut(0.0025);
    antiSigmaCuts->SetSigmaMassCut(0.0025);
  } else if (suffix == "18") {
    sigmaCuts->SetSigmaMassCut(0.00275);
    antiSigmaCuts->SetSigmaMassCut(0.00275);
  } else if (suffix == "19") {
    sigmaCuts->SetSigmaMassCut(0.00325);
    antiSigmaCuts->SetSigmaMassCut(0.00325);
  } else if (suffix == "20") {
    sigmaCuts->SetSigmaMassCut(0.0035);
    antiSigmaCuts->SetSigmaMassCut(0.0035);
  } else if (suffix == "21") {
    sigmaCuts->SetSigmaMassCut(0.00375);
    antiSigmaCuts->SetSigmaMassCut(0.00375);
  } else if (suffix == "22") {
    sigmaCuts->SetSigmaMassCut(0.004);
    antiSigmaCuts->SetSigmaMassCut(0.004);
  } else if (suffix == "23") {
    sigmaCuts->SetSigmaMassCut(0.00425);
    antiSigmaCuts->SetSigmaMassCut(0.00425);
  } else if (suffix == "24") {
    sigmaCuts->SetSigmaMassCut(0.0045);
    antiSigmaCuts->SetSigmaMassCut(0.0045);
  } else if (suffix == "25") {
    sigmaCuts->SetSigmaMassCut(0.00475);
    antiSigmaCuts->SetSigmaMassCut(0.00475);
  } else if (suffix == "26") {
    sigmaCuts->SetSigmaMassCut(0.005);
    antiSigmaCuts->SetSigmaMassCut(0.005);
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

  // do extended QA for the pairs in default mode
  if (suffix == "0") {
    NBins[0] = 750;  // pp
    NBins[8] = 750;  // barp barp
  }
  pairQA[0] = 11;   // pp
  pairQA[2] = 14;   // pSigma
  pairQA[8] = 11;   // barp barp
  pairQA[10] = 14;  // barp bSigma

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
    // no close pair rejection since we don't care about pp
    config->SetDeltaEtaMax(0.);
    config->SetDeltaPhiMax(0.);
    config->SetClosePairRejection(closeRejection);
  }

  if (suffix == "0" && fullBlastQA) {
    config->SetPhiEtaBinnign(true);
    config->SetkTBinning(true);
    config->SetmTBinning(true);
  }
  config->SetPtQA(true);
  config->SetdPhidEtaPlots(false);
  config->SetPDGCodes(PDGParticles);
  config->SetNBinsHist(NBins);
  config->SetMinKRel(kMin);
  config->SetMaxKRel(kMax);
  config->SetMixingDepth(10);
  config->SetUseEventMixing(true);
  config->SetMultiplicityEstimator(AliFemtoDreamEvent::kRef08);
  config->SetMinimalBookingME(false);

  AliAnalysisTaskNanoAODSigma0Femto *task =
      new AliAnalysisTaskNanoAODSigma0Femto("AliAnalysisTaskNanoAODSigma0Femto", isMC);

  task->SetEventCuts(evtCuts);
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
    addon += "MBSigma0";
  } else if (trigger == "kHighMultV0") {
    addon += "HMSigma0";
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
