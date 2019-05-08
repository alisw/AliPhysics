#include "TROOT.h"
#include "TSystem.h"
#include <stdio.h>

AliAnalysisTaskSE* AddTaskFemtoDreamSysVar(bool isMC = false,
                                           bool isESD = false, TString CentEst =
                                               "kInt7",
                                           bool notpp = true,  //1
                                           bool fineBinning = true,  //2
                                           bool kTBinning = false,  //3
                                           bool kTCentBinning = false,  //4
                                           bool mTBinning = false,  //5
                                           bool minimalBooking = true,  // 6
                                           const char *swuffix = "")  //7
                                           {
  TString suffix = Form("%s", swuffix);
  bool DCAPlots = false;
  bool CPAPlots = false;
  bool CombSigma = false;
  bool ContributionSplitting = false;
  bool ContributionSplittingDaug = false;
  bool PileUpRej = true;
  bool eventMixing = true;
  bool phiSpin = false;
  // the manager is static, so get the existing manager via the static method
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  if (!mgr) {
    printf("No analysis manager to connect to!\n");
    return nullptr;
  }
  // just to see if all went well, check if the input event handler has been
  // connected
  if (!mgr->GetInputEventHandler()) {
    printf("This task requires an input event handler!\n");
    return nullptr;
  }

  if (!(AliPIDResponse*) mgr->GetTask("PIDResponseTask")) {
    if (isMC) {
      // IMPORTANT - SET WHEN USING DIFFERENT PASS
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
  }
  AliFemtoDreamEventCuts *evtCuts = AliFemtoDreamEventCuts::StandardCutsRun2();
  evtCuts->SetMinimalBooking(true);

  //Track Cuts
  AliFemtoDreamTrackCuts *TrackCuts = AliFemtoDreamTrackCuts::PrimProtonCuts(
      isMC, DCAPlots, CombSigma, ContributionSplitting);
  TrackCuts->SetCutCharge(1);
  TrackCuts->SetMinimalBooking(minimalBooking);
  AliFemtoDreamTrackCuts *AntiTrackCuts =
      AliFemtoDreamTrackCuts::PrimProtonCuts(isMC, DCAPlots, CombSigma,
                                             ContributionSplitting);
  AntiTrackCuts->SetCutCharge(-1);
  AntiTrackCuts->SetMinimalBooking(minimalBooking);

  AliFemtoDreamv0Cuts *v0Cuts;
  AliFemtoDreamv0Cuts *Antiv0Cuts;
  AliFemtoDreamCascadeCuts *CascadeCuts;
  AliFemtoDreamCascadeCuts *AntiCascadeCuts;

  //Lambda Cuts
  v0Cuts = AliFemtoDreamv0Cuts::LambdaCuts(isMC, CPAPlots,
                                           ContributionSplitting);
  AliFemtoDreamTrackCuts *Posv0Daug = AliFemtoDreamTrackCuts::DecayProtonCuts(
      isMC, PileUpRej, false);
  AliFemtoDreamTrackCuts *Negv0Daug = AliFemtoDreamTrackCuts::DecayPionCuts(
      isMC, PileUpRej, false);
  v0Cuts->SetMinimalBooking(minimalBooking);
  Antiv0Cuts = AliFemtoDreamv0Cuts::LambdaCuts(isMC, CPAPlots,
                                               ContributionSplitting);
  Antiv0Cuts->SetMinimalBooking(minimalBooking);
  AliFemtoDreamTrackCuts *PosAntiv0Daug = AliFemtoDreamTrackCuts::DecayPionCuts(
      isMC, PileUpRej, false);
  PosAntiv0Daug->SetCutCharge(1);
  AliFemtoDreamTrackCuts *NegAntiv0Daug =
      AliFemtoDreamTrackCuts::DecayProtonCuts(isMC, PileUpRej, false);
  NegAntiv0Daug->SetCutCharge(-1);

  v0Cuts->SetPosDaugterTrackCuts(Posv0Daug);
  v0Cuts->SetNegDaugterTrackCuts(Negv0Daug);
  v0Cuts->SetPDGCodePosDaug(2212);  //Proton
  v0Cuts->SetPDGCodeNegDaug(211);  //Pion
  v0Cuts->SetPDGCodev0(3122);  //Lambda

  Antiv0Cuts->SetPosDaugterTrackCuts(PosAntiv0Daug);
  Antiv0Cuts->SetNegDaugterTrackCuts(NegAntiv0Daug);
  Antiv0Cuts->SetPDGCodePosDaug(211);  //Pion
  Antiv0Cuts->SetPDGCodeNegDaug(2212);  //Proton
  Antiv0Cuts->SetPDGCodev0(-3122);  //Lambda

  //Cascade Cuts
  CascadeCuts = AliFemtoDreamCascadeCuts::XiCuts(isMC, ContributionSplitting);
  CascadeCuts->SetMinimalBooking(minimalBooking);
  CascadeCuts->SetXiCharge(-1);
  AliFemtoDreamTrackCuts *XiNegCuts = AliFemtoDreamTrackCuts::Xiv0PionCuts(
      isMC, PileUpRej, false);
  AliFemtoDreamTrackCuts *XiPosCuts = AliFemtoDreamTrackCuts::Xiv0ProtonCuts(
      isMC, PileUpRej, false);
  AliFemtoDreamTrackCuts *XiBachCuts = AliFemtoDreamTrackCuts::XiBachPionCuts(
      isMC, PileUpRej, false);

  AntiCascadeCuts = AliFemtoDreamCascadeCuts::XiCuts(isMC,
                                                     ContributionSplitting);
  AntiCascadeCuts->SetXiCharge(1);
  AntiCascadeCuts->SetMinimalBooking(minimalBooking);
  AliFemtoDreamTrackCuts *AntiXiNegCuts =
      AliFemtoDreamTrackCuts::Xiv0ProtonCuts(isMC, PileUpRej, false);
  AntiXiNegCuts->SetCutCharge(-1);
  AliFemtoDreamTrackCuts *AntiXiPosCuts = AliFemtoDreamTrackCuts::Xiv0PionCuts(
      isMC, PileUpRej, false);
  AntiXiPosCuts->SetCutCharge(1);
  AliFemtoDreamTrackCuts *AntiXiBachCuts =
      AliFemtoDreamTrackCuts::XiBachPionCuts(isMC, PileUpRej, false);
  AntiXiBachCuts->SetCutCharge(1);

  //wide wariation
  if (suffix == "101") {
    CascadeCuts->SetXiMassRange(1.360, 0.028);
    AntiCascadeCuts->SetXiMassRange(1.360, 0.028);
  } else if (suffix == "102") {
    CascadeCuts->SetXiMassRange(1.284, 0.028);
    AntiCascadeCuts->SetXiMassRange(1.284, 0.028);
    //thinner wariation
  } else if (suffix == "103") {
    CascadeCuts->SetXiMassRange(1.346, 0.014);
    AntiCascadeCuts->SetXiMassRange(1.346, 0.014);
  } else if (suffix == "104") {
    CascadeCuts->SetXiMassRange(1.298, 0.014);
    AntiCascadeCuts->SetXiMassRange(1.298, 0.014);
    //Further away from the peak
  } else if (suffix == "105") {
    CascadeCuts->SetXiMassRange(1.360, 0.014);
    AntiCascadeCuts->SetXiMassRange(1.360, 0.014);
  } else if (suffix == "106") {
    CascadeCuts->SetXiMassRange(1.284, 0.014);
    AntiCascadeCuts->SetXiMassRange(1.284, 0.014);
    //furthest away from the peak
  } else if (suffix == "107") {
    CascadeCuts->SetXiMassRange(1.374, 0.014);
    AntiCascadeCuts->SetXiMassRange(1.374, 0.014);
  } else if (suffix == "108") {
    CascadeCuts->SetXiMassRange(1.270, 0.014);
    AntiCascadeCuts->SetXiMassRange(1.270, 0.014);
  }

  CascadeCuts->Setv0Negcuts(XiNegCuts);
  CascadeCuts->Setv0PosCuts(XiPosCuts);
  CascadeCuts->SetBachCuts(XiBachCuts);
  CascadeCuts->SetPDGCodeCasc(3312);
  CascadeCuts->SetPDGCodev0(3122);
  CascadeCuts->SetPDGCodePosDaug(2212);
  CascadeCuts->SetPDGCodeNegDaug(-211);
  CascadeCuts->SetPDGCodeBach(-211);

  AntiCascadeCuts->Setv0Negcuts(AntiXiNegCuts);
  AntiCascadeCuts->Setv0PosCuts(AntiXiPosCuts);
  AntiCascadeCuts->SetBachCuts(AntiXiBachCuts);
  AntiCascadeCuts->SetPDGCodeCasc(-3312);
  AntiCascadeCuts->SetPDGCodev0(-3122);
  AntiCascadeCuts->SetPDGCodePosDaug(211);
  AntiCascadeCuts->SetPDGCodeNegDaug(-2212);
  AntiCascadeCuts->SetPDGCodeBach(-211);

  //Thanks, CINT - will not compile due to an illegal constructor
  std::vector<int> PDGParticles;
  PDGParticles.push_back(2212);
  PDGParticles.push_back(2212);
  PDGParticles.push_back(3122);
  PDGParticles.push_back(3122);
  PDGParticles.push_back(3312);
  PDGParticles.push_back(3312);
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
  if (fineBinning) {
    NBins.push_back(250);  // p p
    NBins.push_back(250);  // p barp
    NBins.push_back(250);  // p Lambda
    NBins.push_back(250);  // p barLambda
    NBins.push_back(250);  // p Xi
    NBins.push_back(250);  // p barXi
    NBins.push_back(250);  // barp barp
    NBins.push_back(250);  // barp Lambda
    NBins.push_back(250);  // barp barLambda
    NBins.push_back(250);  // barp Xi
    NBins.push_back(250);  // barp barXi
    NBins.push_back(250);  // Lambda Lambda
    NBins.push_back(250);  // Lambda barLambda
    NBins.push_back(250);  // Lambda Xi
    NBins.push_back(250);  // Lambda barXi
    NBins.push_back(250);  // barLambda barLambda
    NBins.push_back(250);  // barLambda Xi
    NBins.push_back(250);  // barLambda barXi
    NBins.push_back(250);  // Xi Xi
    NBins.push_back(250);  // Xi barXi
    NBins.push_back(250);  // barXi barXi
  } else {  //standard binning Run1
    NBins.push_back(750);  // p p
    NBins.push_back(750);  // p barp
    NBins.push_back(150);  // p Lambda
    NBins.push_back(150);  // p barLambda
    NBins.push_back(150);  // p Xi
    NBins.push_back(150);  // p barXi
    NBins.push_back(750);  // barp barp
    NBins.push_back(150);  // barp Lambda
    NBins.push_back(150);  // barp barLambda
    NBins.push_back(150);  // barp Xi
    NBins.push_back(150);  // barp barXi
    NBins.push_back(150);  // Lambda Lambda
    NBins.push_back(150);  // Lambda barLambda
    NBins.push_back(150);  // Lambda Xi
    NBins.push_back(150);  // Lambda barXi
    NBins.push_back(150);  // barLambda barLambda
    NBins.push_back(150);  // barLambda Xi
    NBins.push_back(150);  // barLambda barXi
    NBins.push_back(150);  // Xi Xi
    NBins.push_back(150);  // Xi barXi
    NBins.push_back(150);  // barXi barXi
  }
  std::vector<float> kMin;
  kMin.push_back(0.);
  kMin.push_back(0.);
  kMin.push_back(0.);
  kMin.push_back(0.);
  kMin.push_back(0.);
  kMin.push_back(0.);
  kMin.push_back(0.);
  kMin.push_back(0.);
  kMin.push_back(0.);
  kMin.push_back(0.);
  kMin.push_back(0.);
  kMin.push_back(0.);
  kMin.push_back(0.);
  kMin.push_back(0.);
  kMin.push_back(0.);
  kMin.push_back(0.);
  kMin.push_back(0.);
  kMin.push_back(0.);
  kMin.push_back(0.);
  kMin.push_back(0.);
  kMin.push_back(0.);
  std::vector<float> kMax;
  kMax.push_back(1.);
  kMax.push_back(1.);
  kMax.push_back(1.);
  kMax.push_back(1.);
  kMax.push_back(1.);
  kMax.push_back(1.);
  kMax.push_back(1.);
  kMax.push_back(1.);
  kMax.push_back(1.);
  kMax.push_back(1.);
  kMax.push_back(1.);
  kMax.push_back(1.);
  kMax.push_back(1.);
  kMax.push_back(1.);
  kMax.push_back(1.);
  kMax.push_back(1.);
  kMax.push_back(1.);
  kMax.push_back(1.);
  kMax.push_back(1.);
  kMax.push_back(1.);
  kMax.push_back(1.);
  AliFemtoDreamCollConfig *config = new AliFemtoDreamCollConfig("Femto",
                                                                "Femto");
  if (notpp) {
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
  config->SetZBins(ZVtxBins);
  config->SetPDGCodes(PDGParticles);
  config->SetNBinsHist(NBins);
  config->SetMinKRel(kMin);
  config->SetMaxKRel(kMax);
  config->SetMixingDepth(10);
  config->SetSpinningDepth(10);
  config->SetkTBinning(kTBinning);
  config->SetkTCentralityBinning(kTCentBinning);
  config->SetmTBinning(true);
  config->SetExtendedQAPairs(config->GetStandardPairs());
  config->SetUseEventMixing(eventMixing);
  config->SetUsePhiSpinning(phiSpin);

  config->SetDeltaEtaMax(0.010);
  config->SetDeltaPhiMax(0.010);

  config->SetClosePairRejection(config->GetStandardPairRejection());
  if (minimalBooking) {
    config->SetMinimalBookingME(true);
    config->SetMinimalBookingSample(true);
  } else {
    config->SetMinimalBookingME(false);
    config->SetMinimalBookingSample(false);
  }
  if (!notpp) {
    config->SetMultiplicityEstimator(AliFemtoDreamEvent::kSPD);
  } else {
    config->SetMultiplicityEstimator(AliFemtoDreamEvent::kRef08);
  }

  TString TaskName = Form("FemtoDream_%s", suffix.Data());
  AliAnalysisTaskFemtoDream *task = new AliAnalysisTaskFemtoDream(
      TaskName.Data(), isESD, isMC);
  if (CentEst == "kInt7") {
    task->SelectCollisionCandidates(AliVEvent::kINT7);
    task->SetMVPileUp(kTRUE);
  } else if (CentEst == "kMB") {
    task->SelectCollisionCandidates(AliVEvent::kMB);
    task->SetMVPileUp(kFALSE);
  } else if (CentEst == "kHM") {
    task->SelectCollisionCandidates(AliVEvent::kHighMultV0);
    task->SetMVPileUp(kFALSE);
  } else {
    std::cout
        << "====================================================================="
        << std::endl;
    std::cout
        << "====================================================================="
        << std::endl;
    std::cout
        << "Centrality Estimator not set, fix it else your Results will be empty!"
        << std::endl;
    std::cout
        << "====================================================================="
        << std::endl;
    std::cout
        << "====================================================================="
        << std::endl;
  }
  task->SetDebugLevel(0);
  task->SetEvtCutQA(false);
  task->SetTrackBufferSize(2000);
  task->SetEventCuts(evtCuts);

  if (suffix == "1") {
    TrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);
    AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);

    TrackCuts->SetNClsTPC(90);
    AntiTrackCuts->SetNClsTPC(90);

    v0Cuts->SetCutCPA(0.995);
    Antiv0Cuts->SetCutCPA(0.995);

    Posv0Daug->SetNClsTPC(80);
    Negv0Daug->SetNClsTPC(80);
    PosAntiv0Daug->SetNClsTPC(80);
    NegAntiv0Daug->SetNClsTPC(80);

    v0Cuts->SetCutDCADaugToPrimVtx(0.06);
    Antiv0Cuts->SetCutDCADaugToPrimVtx(0.06);

    CascadeCuts->SetCutXiDaughterDCA(1.9);
    AntiCascadeCuts->SetCutXiDaughterDCA(1.9);

    CascadeCuts->SetCutXiMinDistBachToPrimVtx(0.06);
    AntiCascadeCuts->SetCutXiMinDistBachToPrimVtx(0.06);

    CascadeCuts->SetCutXiTransverseRadius(0.6, 200);
    AntiCascadeCuts->SetCutXiTransverseRadius(0.6, 200);

    CascadeCuts->SetCutv0MaxDaughterDCA(1.4);
    AntiCascadeCuts->SetCutv0MaxDaughterDCA(1.4);

    CascadeCuts->SetCutv0CPA(0.96);
    AntiCascadeCuts->SetCutv0CPA(0.96);

    CascadeCuts->SetCutv0TransverseRadius(1.7, 200);
    AntiCascadeCuts->SetCutv0TransverseRadius(1.7, 200);

    CascadeCuts->SetCutv0MinDistToPrimVtx(0.08);
    AntiCascadeCuts->SetCutv0MinDistToPrimVtx(0.08);

    XiNegCuts->SetEtaRange(-0.83, 0.83);
    XiPosCuts->SetEtaRange(-0.83, 0.83);
    XiBachCuts->SetEtaRange(-0.83, 0.83);
    AntiXiNegCuts->SetEtaRange(-0.83, 0.83);
    AntiXiPosCuts->SetEtaRange(-0.83, 0.83);
    AntiXiBachCuts->SetEtaRange(-0.83, 0.83);

    XiNegCuts->SetPID(AliPID::kPion, 999, 3);
    XiPosCuts->SetPID(AliPID::kProton, 999, 3);
    XiBachCuts->SetPID(AliPID::kPion, 999, 3);
    AntiXiNegCuts->SetPID(AliPID::kProton, 999, 3);
    AntiXiPosCuts->SetPID(AliPID::kPion, 999, 3);
    AntiXiBachCuts->SetPID(AliPID::kPion, 999, 3);

    CascadeCuts->SetPtRangeXi(0.8, 999.5);
    AntiCascadeCuts->SetPtRangeXi(0.8, 999.5);

  } else if (suffix == "2") {
    TrackCuts->SetPtRange(0.6, 4.05);
    AntiTrackCuts->SetPtRange(0.6, 4.05);

    TrackCuts->SetEtaRange(-0.77, 0.77);
    AntiTrackCuts->SetEtaRange(-0.77, 0.77);

    config->SetDeltaEtaMax(0.012);
    config->SetDeltaPhiMax(0.012);

    Posv0Daug->SetPID(AliPID::kProton, 999.9, 4);
    Negv0Daug->SetPID(AliPID::kPion, 999.9, 4);
    PosAntiv0Daug->SetPID(AliPID::kPion, 999.9, 4);
    NegAntiv0Daug->SetPID(AliPID::kProton, 999.9, 4);

    CascadeCuts->SetCutXiDaughterDCA(1.9);
    AntiCascadeCuts->SetCutXiDaughterDCA(1.9);

    CascadeCuts->SetCutXiMinDistBachToPrimVtx(0.04);
    AntiCascadeCuts->SetCutXiMinDistBachToPrimVtx(0.04);

    CascadeCuts->SetCutv0CPA(0.99);
    AntiCascadeCuts->SetCutv0CPA(0.99);

    CascadeCuts->SetCutv0MinDistToPrimVtx(0.06);
    AntiCascadeCuts->SetCutv0MinDistToPrimVtx(0.06);

    CascadeCuts->SetCutv0MinDaugDistToPrimVtx(0.06);
    AntiCascadeCuts->SetCutv0MinDaugDistToPrimVtx(0.06);

    XiNegCuts->SetPID(AliPID::kPion, 999, 3);
    XiPosCuts->SetPID(AliPID::kProton, 999, 3);
    XiBachCuts->SetPID(AliPID::kPion, 999, 3);
    AntiXiNegCuts->SetPID(AliPID::kProton, 999, 3);
    AntiXiPosCuts->SetPID(AliPID::kPion, 999, 3);
    AntiXiBachCuts->SetPID(AliPID::kPion, 999, 3);
  } else if (suffix == "3") {
    v0Cuts->SetCutCPA(0.995);
    Antiv0Cuts->SetCutCPA(0.995);

    Posv0Daug->SetEtaRange(-0.77, 0.77);
    Negv0Daug->SetEtaRange(-0.77, 0.77);
    PosAntiv0Daug->SetEtaRange(-0.77, 0.77);
    NegAntiv0Daug->SetEtaRange(-0.77, 0.77);

    v0Cuts->SetCutDCADaugTov0Vtx(1.2);
    Antiv0Cuts->SetCutDCADaugTov0Vtx(1.2);

    CascadeCuts->SetCutXiDaughterDCA(1.3);
    AntiCascadeCuts->SetCutXiDaughterDCA(1.3);

    CascadeCuts->SetCutXiMinDistBachToPrimVtx(0.06);
    AntiCascadeCuts->SetCutXiMinDistBachToPrimVtx(0.06);

    CascadeCuts->SetCutXiCPA(0.99);
    AntiCascadeCuts->SetCutXiCPA(0.99);

    CascadeCuts->SetCutv0MaxDaughterDCA(1.3);
    AntiCascadeCuts->SetCutv0MaxDaughterDCA(1.3);

    CascadeCuts->SetCutv0CPA(0.99);
    AntiCascadeCuts->SetCutv0CPA(0.99);

    CascadeCuts->SetCutv0TransverseRadius(1.7, 200);
    AntiCascadeCuts->SetCutv0TransverseRadius(1.7, 200);

    CascadeCuts->SetCutv0MinDistToPrimVtx(0.08);
    AntiCascadeCuts->SetCutv0MinDistToPrimVtx(0.08);

    CascadeCuts->SetCutv0MinDaugDistToPrimVtx(0.06);
    AntiCascadeCuts->SetCutv0MinDaugDistToPrimVtx(0.06);

    CascadeCuts->SetPtRangeXi(1.2, 999.5);
    AntiCascadeCuts->SetPtRangeXi(1.2, 999.5);

  } else if (suffix == "4") {
    TrackCuts->SetEtaRange(-0.77, 0.77);
    AntiTrackCuts->SetEtaRange(-0.77, 0.77);

    TrackCuts->SetNClsTPC(70);
    AntiTrackCuts->SetNClsTPC(70);

    Posv0Daug->SetNClsTPC(80);
    Negv0Daug->SetNClsTPC(80);
    PosAntiv0Daug->SetNClsTPC(80);
    NegAntiv0Daug->SetNClsTPC(80);

    v0Cuts->SetCutDCADaugToPrimVtx(0.06);
    Antiv0Cuts->SetCutDCADaugToPrimVtx(0.06);

    CascadeCuts->SetCutXiDaughterDCA(1.9);
    AntiCascadeCuts->SetCutXiDaughterDCA(1.9);

    CascadeCuts->SetCutXiTransverseRadius(1.0, 200);
    AntiCascadeCuts->SetCutXiTransverseRadius(1.0, 200);

    CascadeCuts->SetCutv0MaxDaughterDCA(1.3);
    AntiCascadeCuts->SetCutv0MaxDaughterDCA(1.3);

    CascadeCuts->SetCutv0CPA(0.99);
    AntiCascadeCuts->SetCutv0CPA(0.99);

    CascadeCuts->SetCutv0MinDaugDistToPrimVtx(0.06);
    AntiCascadeCuts->SetCutv0MinDaugDistToPrimVtx(0.06);

    XiNegCuts->SetPID(AliPID::kPion, 999, 4.5);
    XiPosCuts->SetPID(AliPID::kProton, 999, 4.5);
    XiBachCuts->SetPID(AliPID::kPion, 999, 4.5);
    AntiXiNegCuts->SetPID(AliPID::kProton, 999, 4.5);
    AntiXiPosCuts->SetPID(AliPID::kPion, 999, 4.5);
    AntiXiBachCuts->SetPID(AliPID::kPion, 999, 4.5);

  } else if (suffix == "5") {
    TrackCuts->SetPtRange(0.4, 4.05);
    AntiTrackCuts->SetPtRange(0.4, 4.05);

    TrackCuts->SetEtaRange(-0.77, 0.77);
    AntiTrackCuts->SetEtaRange(-0.77, 0.77);

    TrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);
    AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);

    config->SetDeltaEtaMax(0.012);
    config->SetDeltaPhiMax(0.012);

    v0Cuts->SetCutCPA(0.995);
    Antiv0Cuts->SetCutCPA(0.995);

    Posv0Daug->SetEtaRange(-0.83, 0.83);
    Negv0Daug->SetEtaRange(-0.83, 0.83);
    PosAntiv0Daug->SetEtaRange(-0.83, 0.83);
    NegAntiv0Daug->SetEtaRange(-0.83, 0.83);

    v0Cuts->SetCutDCADaugTov0Vtx(1.2);
    Antiv0Cuts->SetCutDCADaugTov0Vtx(1.2);

    CascadeCuts->SetCutXiMinDistBachToPrimVtx(0.06);
    AntiCascadeCuts->SetCutXiMinDistBachToPrimVtx(0.06);

    CascadeCuts->SetCutXiCPA(0.99);
    AntiCascadeCuts->SetCutXiCPA(0.99);

    CascadeCuts->SetCutXiTransverseRadius(1.0, 200);
    AntiCascadeCuts->SetCutXiTransverseRadius(1.0, 200);

    CascadeCuts->SetCutv0MaxDaughterDCA(1.4);
    AntiCascadeCuts->SetCutv0MaxDaughterDCA(1.4);

    CascadeCuts->SetCutv0CPA(0.96);
    AntiCascadeCuts->SetCutv0CPA(0.96);

    CascadeCuts->SetCutv0TransverseRadius(1.7, 200);
    AntiCascadeCuts->SetCutv0TransverseRadius(1.7, 200);

    CascadeCuts->SetCutv0MinDaugDistToPrimVtx(0.06);
    AntiCascadeCuts->SetCutv0MinDaugDistToPrimVtx(0.06);

    XiNegCuts->SetEtaRange(-0.77, 0.77);
    XiPosCuts->SetEtaRange(-0.77, 0.77);
    XiBachCuts->SetEtaRange(-0.77, 0.77);
    AntiXiNegCuts->SetEtaRange(-0.77, 0.77);
    AntiXiPosCuts->SetEtaRange(-0.77, 0.77);
    AntiXiBachCuts->SetEtaRange(-0.77, 0.77);

    XiNegCuts->SetPID(AliPID::kPion, 999, 4.5);
    XiPosCuts->SetPID(AliPID::kProton, 999, 4.5);
    XiBachCuts->SetPID(AliPID::kPion, 999, 4.5);
    AntiXiNegCuts->SetPID(AliPID::kProton, 999, 4.5);
    AntiXiPosCuts->SetPID(AliPID::kPion, 999, 4.5);
    AntiXiBachCuts->SetPID(AliPID::kPion, 999, 4.5);

    CascadeCuts->SetPtRangeXi(0.8, 999.5);
    AntiCascadeCuts->SetPtRangeXi(0.8, 999.5);
  } else if (suffix == "6") {
    TrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);
    AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);

    TrackCuts->SetNClsTPC(90);
    AntiTrackCuts->SetNClsTPC(90);

    v0Cuts->SetCutCPA(0.995);
    Antiv0Cuts->SetCutCPA(0.995);

    Posv0Daug->SetPID(AliPID::kProton, 999.9, 4);
    Negv0Daug->SetPID(AliPID::kPion, 999.9, 4);
    PosAntiv0Daug->SetPID(AliPID::kPion, 999.9, 4);
    NegAntiv0Daug->SetPID(AliPID::kProton, 999.9, 4);

    Posv0Daug->SetEtaRange(-0.83, 0.83);
    Negv0Daug->SetEtaRange(-0.83, 0.83);
    PosAntiv0Daug->SetEtaRange(-0.83, 0.83);
    NegAntiv0Daug->SetEtaRange(-0.83, 0.83);

    CascadeCuts->SetCutXiDaughterDCA(1.9);
    AntiCascadeCuts->SetCutXiDaughterDCA(1.9);

    CascadeCuts->SetCutXiMinDistBachToPrimVtx(0.06);
    AntiCascadeCuts->SetCutXiMinDistBachToPrimVtx(0.06);

    CascadeCuts->SetCutXiTransverseRadius(1.0, 200);
    AntiCascadeCuts->SetCutXiTransverseRadius(1.0, 200);

    CascadeCuts->SetCutv0MaxDaughterDCA(1.4);
    AntiCascadeCuts->SetCutv0MaxDaughterDCA(1.4);

    CascadeCuts->SetCutv0TransverseRadius(1.7, 200);
    AntiCascadeCuts->SetCutv0TransverseRadius(1.7, 200);

    CascadeCuts->SetCutv0MinDistToPrimVtx(0.08);
    AntiCascadeCuts->SetCutv0MinDistToPrimVtx(0.08);

    CascadeCuts->SetCutv0MinDaugDistToPrimVtx(0.06);
    AntiCascadeCuts->SetCutv0MinDaugDistToPrimVtx(0.06);

    XiNegCuts->SetEtaRange(-0.77, 0.77);
    XiPosCuts->SetEtaRange(-0.77, 0.77);
    XiBachCuts->SetEtaRange(-0.77, 0.77);
    AntiXiNegCuts->SetEtaRange(-0.77, 0.77);
    AntiXiPosCuts->SetEtaRange(-0.77, 0.77);
    AntiXiBachCuts->SetEtaRange(-0.77, 0.77);

    CascadeCuts->SetPtRangeXi(1.2, 999.5);
    AntiCascadeCuts->SetPtRangeXi(1.2, 999.5);
  } else if (suffix == "7") {
    TrackCuts->SetEtaRange(-0.85, 0.85);
    AntiTrackCuts->SetEtaRange(-0.85, 0.85);

    TrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);
    AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);

    config->SetDeltaEtaMax(0.012);
    config->SetDeltaPhiMax(0.012);

    v0Cuts->SetCutCPA(0.995);
    Antiv0Cuts->SetCutCPA(0.995);

    Posv0Daug->SetEtaRange(-0.83, 0.83);
    Negv0Daug->SetEtaRange(-0.83, 0.83);
    PosAntiv0Daug->SetEtaRange(-0.83, 0.83);
    NegAntiv0Daug->SetEtaRange(-0.83, 0.83);

    v0Cuts->SetCutDCADaugTov0Vtx(1.2);
    Antiv0Cuts->SetCutDCADaugTov0Vtx(1.2);

    v0Cuts->SetCutDCADaugToPrimVtx(0.06);
    Antiv0Cuts->SetCutDCADaugToPrimVtx(0.06);

    CascadeCuts->SetCutv0MinDistToPrimVtx(0.08);
    AntiCascadeCuts->SetCutv0MinDistToPrimVtx(0.08);

    XiNegCuts->SetEtaRange(-0.77, 0.77);
    XiPosCuts->SetEtaRange(-0.77, 0.77);
    XiBachCuts->SetEtaRange(-0.77, 0.77);
    AntiXiNegCuts->SetEtaRange(-0.77, 0.77);
    AntiXiPosCuts->SetEtaRange(-0.77, 0.77);
    AntiXiBachCuts->SetEtaRange(-0.77, 0.77);

    XiNegCuts->SetPID(AliPID::kPion, 999, 3.5);
    XiPosCuts->SetPID(AliPID::kProton, 999, 3.5);
    XiBachCuts->SetPID(AliPID::kPion, 999, 3.5);
    AntiXiNegCuts->SetPID(AliPID::kProton, 999, 3.5);
    AntiXiPosCuts->SetPID(AliPID::kPion, 999, 3.5);
    AntiXiBachCuts->SetPID(AliPID::kPion, 999, 3.5);

    CascadeCuts->SetPtRangeXi(1.2, 999.5);
    AntiCascadeCuts->SetPtRangeXi(1.2, 999.5);

  } else if (suffix == "8") {
    TrackCuts->SetPtRange(0.4, 4.05);
    AntiTrackCuts->SetPtRange(0.4, 4.05);

    TrackCuts->SetEtaRange(-0.77, 0.77);
    AntiTrackCuts->SetEtaRange(-0.77, 0.77);

    TrackCuts->SetNClsTPC(90);
    AntiTrackCuts->SetNClsTPC(90);

    config->SetDeltaEtaMax(0.012);
    config->SetDeltaPhiMax(0.012);

    Posv0Daug->SetEtaRange(-0.83, 0.83);
    Negv0Daug->SetEtaRange(-0.83, 0.83);
    PosAntiv0Daug->SetEtaRange(-0.83, 0.83);
    NegAntiv0Daug->SetEtaRange(-0.83, 0.83);

    CascadeCuts->SetCutXiDaughterDCA(1.3);
    AntiCascadeCuts->SetCutXiDaughterDCA(1.3);

    CascadeCuts->SetCutXiMinDistBachToPrimVtx(0.06);
    AntiCascadeCuts->SetCutXiMinDistBachToPrimVtx(0.06);

    CascadeCuts->SetCutXiTransverseRadius(1.0, 200);
    AntiCascadeCuts->SetCutXiTransverseRadius(1.0, 200);

    CascadeCuts->SetCutv0CPA(0.96);
    AntiCascadeCuts->SetCutv0CPA(0.96);

    CascadeCuts->SetCutv0MinDaugDistToPrimVtx(0.06);
    AntiCascadeCuts->SetCutv0MinDaugDistToPrimVtx(0.06);

    XiNegCuts->SetEtaRange(-0.77, 0.77);
    XiPosCuts->SetEtaRange(-0.77, 0.77);
    XiBachCuts->SetEtaRange(-0.77, 0.77);
    AntiXiNegCuts->SetEtaRange(-0.77, 0.77);
    AntiXiPosCuts->SetEtaRange(-0.77, 0.77);
    AntiXiBachCuts->SetEtaRange(-0.77, 0.77);

    XiNegCuts->SetPID(AliPID::kPion, 999, 4.5);
    XiPosCuts->SetPID(AliPID::kProton, 999, 4.5);
    XiBachCuts->SetPID(AliPID::kPion, 999, 4.5);
    AntiXiNegCuts->SetPID(AliPID::kProton, 999, 4.5);
    AntiXiPosCuts->SetPID(AliPID::kPion, 999, 4.5);
    AntiXiBachCuts->SetPID(AliPID::kPion, 999, 4.5);

    CascadeCuts->SetPtRangeXi(1.2, 999.5);
    AntiCascadeCuts->SetPtRangeXi(1.2, 999.5);
  } else if (suffix == "9") {
    Posv0Daug->SetNClsTPC(80);
    Negv0Daug->SetNClsTPC(80);
    PosAntiv0Daug->SetNClsTPC(80);
    NegAntiv0Daug->SetNClsTPC(80);

    Posv0Daug->SetEtaRange(-0.77, 0.77);
    Negv0Daug->SetEtaRange(-0.77, 0.77);
    PosAntiv0Daug->SetEtaRange(-0.77, 0.77);
    NegAntiv0Daug->SetEtaRange(-0.77, 0.77);

    v0Cuts->SetCutDCADaugTov0Vtx(1.2);
    Antiv0Cuts->SetCutDCADaugTov0Vtx(1.2);

    CascadeCuts->SetCutXiDaughterDCA(1.3);
    AntiCascadeCuts->SetCutXiDaughterDCA(1.3);

    CascadeCuts->SetCutXiCPA(0.99);
    AntiCascadeCuts->SetCutXiCPA(0.99);

    CascadeCuts->SetCutXiTransverseRadius(1.0, 200);
    AntiCascadeCuts->SetCutXiTransverseRadius(1.0, 200);

    CascadeCuts->SetCutv0MaxDaughterDCA(1.4);
    AntiCascadeCuts->SetCutv0MaxDaughterDCA(1.4);

    CascadeCuts->SetCutv0TransverseRadius(1.7, 200);
    AntiCascadeCuts->SetCutv0TransverseRadius(1.7, 200);

    CascadeCuts->SetCutv0MinDistToPrimVtx(0.06);
    AntiCascadeCuts->SetCutv0MinDistToPrimVtx(0.06);

    CascadeCuts->SetCutv0MinDaugDistToPrimVtx(0.06);
    AntiCascadeCuts->SetCutv0MinDaugDistToPrimVtx(0.06);

    XiNegCuts->SetPID(AliPID::kPion, 999, 3);
    XiPosCuts->SetPID(AliPID::kProton, 999, 3);
    XiBachCuts->SetPID(AliPID::kPion, 999, 3);
    AntiXiNegCuts->SetPID(AliPID::kProton, 999, 3);
    AntiXiPosCuts->SetPID(AliPID::kPion, 999, 3);
    AntiXiBachCuts->SetPID(AliPID::kPion, 999, 3);

    CascadeCuts->SetPtRangeXi(0.8, 999.5);
    AntiCascadeCuts->SetPtRangeXi(0.8, 999.5);

  } else if (suffix == "10") {
    TrackCuts->SetPtRange(0.4, 4.05);
    AntiTrackCuts->SetPtRange(0.4, 4.05);

    TrackCuts->SetEtaRange(-0.85, 0.85);
    AntiTrackCuts->SetEtaRange(-0.85, 0.85);

    Posv0Daug->SetPID(AliPID::kProton, 999.9, 4);
    Negv0Daug->SetPID(AliPID::kPion, 999.9, 4);
    PosAntiv0Daug->SetPID(AliPID::kPion, 999.9, 4);
    NegAntiv0Daug->SetPID(AliPID::kProton, 999.9, 4);

    Posv0Daug->SetEtaRange(-0.77, 0.77);
    Negv0Daug->SetEtaRange(-0.77, 0.77);
    PosAntiv0Daug->SetEtaRange(-0.77, 0.77);
    NegAntiv0Daug->SetEtaRange(-0.77, 0.77);

    v0Cuts->SetCutDCADaugToPrimVtx(0.06);
    Antiv0Cuts->SetCutDCADaugToPrimVtx(0.06);

    CascadeCuts->SetCutXiDaughterDCA(1.3);
    AntiCascadeCuts->SetCutXiDaughterDCA(1.3);

    CascadeCuts->SetCutXiCPA(0.99);
    AntiCascadeCuts->SetCutXiCPA(0.99);

    CascadeCuts->SetCutXiTransverseRadius(0.6, 200);
    AntiCascadeCuts->SetCutXiTransverseRadius(0.6, 200);

    CascadeCuts->SetCutv0MinDistToPrimVtx(0.08);
    AntiCascadeCuts->SetCutv0MinDistToPrimVtx(0.08);

    CascadeCuts->SetCutv0MinDaugDistToPrimVtx(0.06);
    AntiCascadeCuts->SetCutv0MinDaugDistToPrimVtx(0.06);

    XiNegCuts->SetEtaRange(-0.87, 0.87);
    XiPosCuts->SetEtaRange(-0.87, 0.87);
    XiBachCuts->SetEtaRange(-0.87, 0.87);
    AntiXiNegCuts->SetEtaRange(-0.87, 0.87);
    AntiXiPosCuts->SetEtaRange(-0.87, 0.87);
    AntiXiBachCuts->SetEtaRange(-0.87, 0.87);

    XiNegCuts->SetPID(AliPID::kPion, 999, 3.5);
    XiPosCuts->SetPID(AliPID::kProton, 999, 3.5);
    XiBachCuts->SetPID(AliPID::kPion, 999, 3.5);
    AntiXiNegCuts->SetPID(AliPID::kProton, 999, 3.5);
    AntiXiPosCuts->SetPID(AliPID::kPion, 999, 3.5);
    AntiXiBachCuts->SetPID(AliPID::kPion, 999, 3.5);

    CascadeCuts->SetPtRangeXi(1.2, 999.5);
    AntiCascadeCuts->SetPtRangeXi(1.2, 999.5);
  } else if (suffix == "11") {
    TrackCuts->SetEtaRange(-0.77, 0.77);
    AntiTrackCuts->SetEtaRange(-0.77, 0.77);

    TrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);
    AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);

    Posv0Daug->SetNClsTPC(80);
    Negv0Daug->SetNClsTPC(80);
    PosAntiv0Daug->SetNClsTPC(80);
    NegAntiv0Daug->SetNClsTPC(80);

    Posv0Daug->SetEtaRange(-0.77, 0.77);
    Negv0Daug->SetEtaRange(-0.77, 0.77);
    PosAntiv0Daug->SetEtaRange(-0.77, 0.77);
    NegAntiv0Daug->SetEtaRange(-0.77, 0.77);

    v0Cuts->SetCutDCADaugTov0Vtx(1.2);
    Antiv0Cuts->SetCutDCADaugTov0Vtx(1.2);

    CascadeCuts->SetCutXiDaughterDCA(1.9);
    AntiCascadeCuts->SetCutXiDaughterDCA(1.9);

    CascadeCuts->SetCutXiMinDistBachToPrimVtx(0.06);
    AntiCascadeCuts->SetCutXiMinDistBachToPrimVtx(0.06);

    CascadeCuts->SetCutXiCPA(0.99);
    AntiCascadeCuts->SetCutXiCPA(0.99);

    CascadeCuts->SetCutXiTransverseRadius(0.6, 200);
    AntiCascadeCuts->SetCutXiTransverseRadius(0.6, 200);

    CascadeCuts->SetCutv0MaxDaughterDCA(1.3);
    AntiCascadeCuts->SetCutv0MaxDaughterDCA(1.3);

    CascadeCuts->SetCutv0CPA(0.99);
    AntiCascadeCuts->SetCutv0CPA(0.99);

    CascadeCuts->SetCutv0TransverseRadius(1.7, 200);
    AntiCascadeCuts->SetCutv0TransverseRadius(1.7, 200);

    CascadeCuts->SetCutv0MinDistToPrimVtx(0.06);
    AntiCascadeCuts->SetCutv0MinDistToPrimVtx(0.06);

    XiNegCuts->SetEtaRange(-0.83, 0.87);
    XiPosCuts->SetEtaRange(-0.83, 0.87);
    XiBachCuts->SetEtaRange(-0.83, 0.87);
    AntiXiNegCuts->SetEtaRange(-0.83, 0.87);
    AntiXiPosCuts->SetEtaRange(-0.83, 0.87);
    AntiXiBachCuts->SetEtaRange(-0.83, 0.87);

    XiNegCuts->SetPID(AliPID::kPion, 999, 4.5);
    XiPosCuts->SetPID(AliPID::kProton, 999, 4.5);
    XiBachCuts->SetPID(AliPID::kPion, 999, 4.5);
    AntiXiNegCuts->SetPID(AliPID::kProton, 999, 4.5);
    AntiXiPosCuts->SetPID(AliPID::kPion, 999, 4.5);
    AntiXiBachCuts->SetPID(AliPID::kPion, 999, 4.5);
  } else if (suffix == "12") {
    TrackCuts->SetPtRange(0.4, 4.05);
    AntiTrackCuts->SetPtRange(0.4, 4.05);

    TrackCuts->SetEtaRange(-0.85, 0.85);
    AntiTrackCuts->SetEtaRange(-0.85, 0.85);

    TrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);
    AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);

    config->SetDeltaEtaMax(0.012);
    config->SetDeltaPhiMax(0.012);

    Posv0Daug->SetPID(AliPID::kProton, 999.9, 4);
    Negv0Daug->SetPID(AliPID::kPion, 999.9, 4);
    PosAntiv0Daug->SetPID(AliPID::kPion, 999.9, 4);
    NegAntiv0Daug->SetPID(AliPID::kProton, 999.9, 4);

    Posv0Daug->SetEtaRange(-0.77, 0.77);
    Negv0Daug->SetEtaRange(-0.77, 0.77);
    PosAntiv0Daug->SetEtaRange(-0.77, 0.77);
    NegAntiv0Daug->SetEtaRange(-0.77, 0.77);

    CascadeCuts->SetCutXiDaughterDCA(1.3);
    AntiCascadeCuts->SetCutXiDaughterDCA(1.3);

    CascadeCuts->SetCutXiMinDistBachToPrimVtx(0.06);
    AntiCascadeCuts->SetCutXiMinDistBachToPrimVtx(0.06);

    CascadeCuts->SetCutv0MaxDaughterDCA(1.3);
    AntiCascadeCuts->SetCutv0MaxDaughterDCA(1.3);

    CascadeCuts->SetCutv0CPA(0.99);
    AntiCascadeCuts->SetCutv0CPA(0.99);

    CascadeCuts->SetCutv0MinDistToPrimVtx(0.08);
    AntiCascadeCuts->SetCutv0MinDistToPrimVtx(0.08);

    CascadeCuts->SetCutv0MinDaugDistToPrimVtx(0.06);
    AntiCascadeCuts->SetCutv0MinDaugDistToPrimVtx(0.06);

    XiNegCuts->SetPID(AliPID::kPion, 999, 3);
    XiPosCuts->SetPID(AliPID::kProton, 999, 3);
    XiBachCuts->SetPID(AliPID::kPion, 999, 3);
    AntiXiNegCuts->SetPID(AliPID::kProton, 999, 3);
    AntiXiPosCuts->SetPID(AliPID::kPion, 999, 3);
    AntiXiBachCuts->SetPID(AliPID::kPion, 999, 3);

  } else if (suffix == "13") {
    TrackCuts->SetEtaRange(-0.85, 0.85);
    AntiTrackCuts->SetEtaRange(-0.85, 0.85);

    TrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);
    AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);

    TrackCuts->SetNClsTPC(90);
    AntiTrackCuts->SetNClsTPC(90);

    Posv0Daug->SetEtaRange(-0.83, 0.83);
    Negv0Daug->SetEtaRange(-0.83, 0.83);
    PosAntiv0Daug->SetEtaRange(-0.83, 0.83);
    NegAntiv0Daug->SetEtaRange(-0.83, 0.83);

    v0Cuts->SetCutDCADaugTov0Vtx(1.2);
    Antiv0Cuts->SetCutDCADaugTov0Vtx(1.2);

    v0Cuts->SetCutDCADaugToPrimVtx(0.06);
    Antiv0Cuts->SetCutDCADaugToPrimVtx(0.06);

    CascadeCuts->SetCutXiDaughterDCA(1.9);
    AntiCascadeCuts->SetCutXiDaughterDCA(1.9);

    CascadeCuts->SetCutXiMinDistBachToPrimVtx(0.06);
    AntiCascadeCuts->SetCutXiMinDistBachToPrimVtx(0.06);

    CascadeCuts->SetCutXiCPA(0.985);
    AntiCascadeCuts->SetCutXiCPA(0.985);

    CascadeCuts->SetCutXiTransverseRadius(1.0, 200);
    AntiCascadeCuts->SetCutXiTransverseRadius(1.0, 200);

    CascadeCuts->SetCutv0CPA(0.96);
    AntiCascadeCuts->SetCutv0CPA(0.96);

    CascadeCuts->SetCutv0TransverseRadius(1.7, 200);
    AntiCascadeCuts->SetCutv0TransverseRadius(1.7, 200);

    CascadeCuts->SetCutv0MinDistToPrimVtx(0.06);
    AntiCascadeCuts->SetCutv0MinDistToPrimVtx(0.06);

    CascadeCuts->SetCutv0MinDaugDistToPrimVtx(0.06);
    AntiCascadeCuts->SetCutv0MinDaugDistToPrimVtx(0.06);

    XiNegCuts->SetEtaRange(-0.77, 0.77);
    XiPosCuts->SetEtaRange(-0.77, 0.77);
    XiBachCuts->SetEtaRange(-0.77, 0.77);
    AntiXiNegCuts->SetEtaRange(-0.77, 0.77);
    AntiXiPosCuts->SetEtaRange(-0.77, 0.77);
    AntiXiBachCuts->SetEtaRange(-0.77, 0.77);

    XiNegCuts->SetPID(AliPID::kPion, 999, 3.5);
    XiPosCuts->SetPID(AliPID::kProton, 999, 3.5);
    XiBachCuts->SetPID(AliPID::kPion, 999, 3.5);
    AntiXiNegCuts->SetPID(AliPID::kProton, 999, 3.5);
    AntiXiPosCuts->SetPID(AliPID::kPion, 999, 3.5);
    AntiXiBachCuts->SetPID(AliPID::kPion, 999, 3.5);

    CascadeCuts->SetPtRangeXi(1.2, 999.5);
    AntiCascadeCuts->SetPtRangeXi(1.2, 999.5);

  } else if (suffix == "14") {
    TrackCuts->SetPtRange(0.4, 4.05);
    AntiTrackCuts->SetPtRange(0.4, 4.05);

    TrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);
    AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);

    config->SetDeltaEtaMax(0.012);
    config->SetDeltaPhiMax(0.012);

    Posv0Daug->SetNClsTPC(80);
    Negv0Daug->SetNClsTPC(80);
    PosAntiv0Daug->SetNClsTPC(80);
    NegAntiv0Daug->SetNClsTPC(80);

    Posv0Daug->SetEtaRange(-0.83, 0.83);
    Negv0Daug->SetEtaRange(-0.83, 0.83);
    PosAntiv0Daug->SetEtaRange(-0.83, 0.83);
    NegAntiv0Daug->SetEtaRange(-0.83, 0.83);

    CascadeCuts->SetCutv0MaxDaughterDCA(1.3);
    AntiCascadeCuts->SetCutv0MaxDaughterDCA(1.3);

    CascadeCuts->SetCutv0TransverseRadius(1.1, 200);
    AntiCascadeCuts->SetCutv0TransverseRadius(1.1, 200);

    CascadeCuts->SetCutv0MinDistToPrimVtx(0.08);
    AntiCascadeCuts->SetCutv0MinDistToPrimVtx(0.08);

    XiNegCuts->SetEtaRange(-0.83, 0.87);
    XiPosCuts->SetEtaRange(-0.83, 0.87);
    XiBachCuts->SetEtaRange(-0.83, 0.87);
    AntiXiNegCuts->SetEtaRange(-0.83, 0.87);
    AntiXiPosCuts->SetEtaRange(-0.83, 0.87);
    AntiXiBachCuts->SetEtaRange(-0.83, 0.87);

  } else if (suffix == "15") {
    TrackCuts->SetPtRange(0.6, 4.05);
    AntiTrackCuts->SetPtRange(0.6, 4.05);

    TrackCuts->SetEtaRange(-0.85, 0.85);
    AntiTrackCuts->SetEtaRange(-0.85, 0.85);

    TrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);
    AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);

    TrackCuts->SetNClsTPC(90);
    AntiTrackCuts->SetNClsTPC(90);

    v0Cuts->SetCutCPA(0.995);
    Antiv0Cuts->SetCutCPA(0.995);

    Posv0Daug->SetNClsTPC(80);
    Negv0Daug->SetNClsTPC(80);
    PosAntiv0Daug->SetNClsTPC(80);
    NegAntiv0Daug->SetNClsTPC(80);

    v0Cuts->SetCutDCADaugTov0Vtx(1.2);
    Antiv0Cuts->SetCutDCADaugTov0Vtx(1.2);

    v0Cuts->SetCutDCADaugToPrimVtx(0.06);
    Antiv0Cuts->SetCutDCADaugToPrimVtx(0.06);

    CascadeCuts->SetCutXiDaughterDCA(1.3);
    AntiCascadeCuts->SetCutXiDaughterDCA(1.3);

    CascadeCuts->SetCutXiMinDistBachToPrimVtx(0.04);
    AntiCascadeCuts->SetCutXiMinDistBachToPrimVtx(0.04);

    CascadeCuts->SetCutXiCPA(0.985);
    AntiCascadeCuts->SetCutXiCPA(0.985);

    CascadeCuts->SetCutXiTransverseRadius(0.6, 200);
    AntiCascadeCuts->SetCutXiTransverseRadius(0.6, 200);

    CascadeCuts->SetCutv0MaxDaughterDCA(1.4);
    AntiCascadeCuts->SetCutv0MaxDaughterDCA(1.4);

    CascadeCuts->SetCutv0CPA(0.99);
    AntiCascadeCuts->SetCutv0CPA(0.99);

    CascadeCuts->SetCutv0TransverseRadius(1.1, 200);
    AntiCascadeCuts->SetCutv0TransverseRadius(1.1, 200);

    CascadeCuts->SetCutv0MinDistToPrimVtx(0.06);
    AntiCascadeCuts->SetCutv0MinDistToPrimVtx(0.06);

    XiNegCuts->SetEtaRange(-0.83, 0.87);
    XiPosCuts->SetEtaRange(-0.83, 0.87);
    XiBachCuts->SetEtaRange(-0.83, 0.87);
    AntiXiNegCuts->SetEtaRange(-0.83, 0.87);
    AntiXiPosCuts->SetEtaRange(-0.83, 0.87);
    AntiXiBachCuts->SetEtaRange(-0.83, 0.87);

    CascadeCuts->SetPtRangeXi(0.8, 999.5);
    AntiCascadeCuts->SetPtRangeXi(0.8, 999.5);

  } else if (suffix == "16") {
    TrackCuts->SetPtRange(0.4, 4.05);
    AntiTrackCuts->SetPtRange(0.4, 4.05);

    TrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);
    AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);

    TrackCuts->SetNClsTPC(70);
    AntiTrackCuts->SetNClsTPC(70);

    config->SetDeltaEtaMax(0.012);
    config->SetDeltaPhiMax(0.012);

    Posv0Daug->SetNClsTPC(80);
    Negv0Daug->SetNClsTPC(80);
    PosAntiv0Daug->SetNClsTPC(80);
    NegAntiv0Daug->SetNClsTPC(80);

    Posv0Daug->SetEtaRange(-0.77, 0.77);
    Negv0Daug->SetEtaRange(-0.77, 0.77);
    PosAntiv0Daug->SetEtaRange(-0.77, 0.77);
    NegAntiv0Daug->SetEtaRange(-0.77, 0.77);

    v0Cuts->SetCutDCADaugToPrimVtx(0.06);
    Antiv0Cuts->SetCutDCADaugToPrimVtx(0.06);

    CascadeCuts->SetCutXiDaughterDCA(1.9);
    AntiCascadeCuts->SetCutXiDaughterDCA(1.9);

    CascadeCuts->SetCutXiTransverseRadius(1.0, 200);
    AntiCascadeCuts->SetCutXiTransverseRadius(1.0, 200);

    CascadeCuts->SetCutv0TransverseRadius(1.1, 200);
    AntiCascadeCuts->SetCutv0TransverseRadius(1.1, 200);

    CascadeCuts->SetPtRangeXi(1.2, 999.5);
    AntiCascadeCuts->SetPtRangeXi(1.2, 999.5);

  } else if (suffix == "17") {
    TrackCuts->SetPtRange(0.4, 4.05);
    AntiTrackCuts->SetPtRange(0.4, 4.05);

    TrackCuts->SetEtaRange(-0.85, 0.85);
    AntiTrackCuts->SetEtaRange(-0.85, 0.85);

    TrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);
    AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);

    config->SetDeltaEtaMax(0.012);
    config->SetDeltaPhiMax(0.012);

    v0Cuts->SetCutCPA(0.995);
    Antiv0Cuts->SetCutCPA(0.995);

    Posv0Daug->SetEtaRange(-0.77, 0.77);
    Negv0Daug->SetEtaRange(-0.77, 0.77);
    PosAntiv0Daug->SetEtaRange(-0.77, 0.77);
    NegAntiv0Daug->SetEtaRange(-0.77, 0.77);

    v0Cuts->SetCutDCADaugToPrimVtx(0.06);
    Antiv0Cuts->SetCutDCADaugToPrimVtx(0.06);

    CascadeCuts->SetCutXiDaughterDCA(1.3);
    AntiCascadeCuts->SetCutXiDaughterDCA(1.3);

    CascadeCuts->SetCutXiMinDistBachToPrimVtx(0.04);
    AntiCascadeCuts->SetCutXiMinDistBachToPrimVtx(0.04);

    CascadeCuts->SetCutXiTransverseRadius(1.0, 200);
    AntiCascadeCuts->SetCutXiTransverseRadius(1.0, 200);

    CascadeCuts->SetCutv0MaxDaughterDCA(1.4);
    AntiCascadeCuts->SetCutv0MaxDaughterDCA(1.4);

    CascadeCuts->SetCutv0CPA(0.99);
    AntiCascadeCuts->SetCutv0CPA(0.99);

    CascadeCuts->SetCutv0TransverseRadius(1.7, 200);
    AntiCascadeCuts->SetCutv0TransverseRadius(1.7, 200);

    CascadeCuts->SetCutv0MinDistToPrimVtx(0.06);
    AntiCascadeCuts->SetCutv0MinDistToPrimVtx(0.06);

    XiNegCuts->SetEtaRange(-0.77, 0.77);
    XiPosCuts->SetEtaRange(-0.77, 0.77);
    XiBachCuts->SetEtaRange(-0.77, 0.77);
    AntiXiNegCuts->SetEtaRange(-0.77, 0.77);
    AntiXiPosCuts->SetEtaRange(-0.77, 0.77);
    AntiXiBachCuts->SetEtaRange(-0.77, 0.77);
  } else if (suffix == "18") {
    TrackCuts->SetPtRange(0.4, 4.05);
    AntiTrackCuts->SetPtRange(0.4, 4.05);

    TrackCuts->SetEtaRange(-0.85, 0.85);
    AntiTrackCuts->SetEtaRange(-0.85, 0.85);

    TrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);
    AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);

    TrackCuts->SetNClsTPC(90);
    AntiTrackCuts->SetNClsTPC(90);

    Posv0Daug->SetNClsTPC(80);
    Negv0Daug->SetNClsTPC(80);
    PosAntiv0Daug->SetNClsTPC(80);
    NegAntiv0Daug->SetNClsTPC(80);

    Posv0Daug->SetEtaRange(-0.77, 0.77);
    Negv0Daug->SetEtaRange(-0.77, 0.77);
    PosAntiv0Daug->SetEtaRange(-0.77, 0.77);
    NegAntiv0Daug->SetEtaRange(-0.77, 0.77);

    CascadeCuts->SetCutXiDaughterDCA(1.9);
    AntiCascadeCuts->SetCutXiDaughterDCA(1.9);

    CascadeCuts->SetCutXiCPA(0.99);
    AntiCascadeCuts->SetCutXiCPA(0.99);

    CascadeCuts->SetCutv0MaxDaughterDCA(1.4);
    AntiCascadeCuts->SetCutv0MaxDaughterDCA(1.4);

    CascadeCuts->SetCutv0CPA(0.96);
    AntiCascadeCuts->SetCutv0CPA(0.96);

    XiNegCuts->SetPID(AliPID::kPion, 999, 3);
    XiPosCuts->SetPID(AliPID::kProton, 999, 3);
    XiBachCuts->SetPID(AliPID::kPion, 999, 3);
    AntiXiNegCuts->SetPID(AliPID::kProton, 999, 3);
    AntiXiPosCuts->SetPID(AliPID::kPion, 999, 3);
    AntiXiBachCuts->SetPID(AliPID::kPion, 999, 3);

    CascadeCuts->SetPtRangeXi(0.8, 999.5);
    AntiCascadeCuts->SetPtRangeXi(0.8, 999.5);

  } else if (suffix == "19") {
    TrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);
    AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);

    TrackCuts->SetNClsTPC(70);
    AntiTrackCuts->SetNClsTPC(70);

    config->SetDeltaEtaMax(0.012);
    config->SetDeltaPhiMax(0.012);

    v0Cuts->SetCutCPA(0.995);
    Antiv0Cuts->SetCutCPA(0.995);

    Posv0Daug->SetNClsTPC(80);
    Negv0Daug->SetNClsTPC(80);
    PosAntiv0Daug->SetNClsTPC(80);
    NegAntiv0Daug->SetNClsTPC(80);

    v0Cuts->SetCutDCADaugTov0Vtx(1.2);
    Antiv0Cuts->SetCutDCADaugTov0Vtx(1.2);

    CascadeCuts->SetCutv0MaxDaughterDCA(1.4);
    AntiCascadeCuts->SetCutv0MaxDaughterDCA(1.4);

    CascadeCuts->SetCutv0TransverseRadius(1.7, 200);
    AntiCascadeCuts->SetCutv0TransverseRadius(1.7, 200);

    CascadeCuts->SetCutv0MinDaugDistToPrimVtx(0.06);
    AntiCascadeCuts->SetCutv0MinDaugDistToPrimVtx(0.06);

    XiNegCuts->SetEtaRange(-0.83, 0.87);
    XiPosCuts->SetEtaRange(-0.83, 0.87);
    XiBachCuts->SetEtaRange(-0.83, 0.87);
    AntiXiNegCuts->SetEtaRange(-0.83, 0.87);
    AntiXiPosCuts->SetEtaRange(-0.83, 0.87);
    AntiXiBachCuts->SetEtaRange(-0.83, 0.87);

    XiNegCuts->SetPID(AliPID::kPion, 999, 3);
    XiPosCuts->SetPID(AliPID::kProton, 999, 3);
    XiBachCuts->SetPID(AliPID::kPion, 999, 3);
    AntiXiNegCuts->SetPID(AliPID::kProton, 999, 3);
    AntiXiPosCuts->SetPID(AliPID::kPion, 999, 3);
    AntiXiBachCuts->SetPID(AliPID::kPion, 999, 3);

    CascadeCuts->SetPtRangeXi(0.8, 999.5);
    AntiCascadeCuts->SetPtRangeXi(0.8, 999.5);

  } else if (suffix == "20") {
    TrackCuts->SetPtRange(0.4, 4.05);
    AntiTrackCuts->SetPtRange(0.4, 4.05);

    TrackCuts->SetNClsTPC(90);
    AntiTrackCuts->SetNClsTPC(90);

    config->SetDeltaEtaMax(0.012);
    config->SetDeltaPhiMax(0.012);

    Posv0Daug->SetNClsTPC(80);
    Negv0Daug->SetNClsTPC(80);
    PosAntiv0Daug->SetNClsTPC(80);
    NegAntiv0Daug->SetNClsTPC(80);

    Posv0Daug->SetEtaRange(-0.83, 0.83);
    Negv0Daug->SetEtaRange(-0.83, 0.83);
    PosAntiv0Daug->SetEtaRange(-0.83, 0.83);
    NegAntiv0Daug->SetEtaRange(-0.83, 0.83);

    CascadeCuts->SetCutXiMinDistBachToPrimVtx(0.06);
    AntiCascadeCuts->SetCutXiMinDistBachToPrimVtx(0.06);

    CascadeCuts->SetCutXiCPA(0.985);
    AntiCascadeCuts->SetCutXiCPA(0.985);

    CascadeCuts->SetCutv0MaxDaughterDCA(1.3);
    AntiCascadeCuts->SetCutv0MaxDaughterDCA(1.3);

    CascadeCuts->SetCutv0CPA(0.99);
    AntiCascadeCuts->SetCutv0CPA(0.99);

    CascadeCuts->SetCutv0MinDistToPrimVtx(0.06);
    AntiCascadeCuts->SetCutv0MinDistToPrimVtx(0.06);

    XiNegCuts->SetEtaRange(-0.87, 0.87);
    XiPosCuts->SetEtaRange(-0.87, 0.87);
    XiBachCuts->SetEtaRange(-0.87, 0.87);
    AntiXiNegCuts->SetEtaRange(-0.87, 0.87);
    AntiXiPosCuts->SetEtaRange(-0.87, 0.87);
    AntiXiBachCuts->SetEtaRange(-0.87, 0.87);

    CascadeCuts->SetPtRangeXi(1.2, 999.5);
    AntiCascadeCuts->SetPtRangeXi(1.2, 999.5);
  } else if (suffix == "21") {

    TrackCuts->SetEtaRange(-0.85, 0.85);
    AntiTrackCuts->SetEtaRange(-0.85, 0.85);

    TrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);
    AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);

    TrackCuts->SetNClsTPC(70);
    AntiTrackCuts->SetNClsTPC(70);

    config->SetDeltaEtaMax(0.012);
    config->SetDeltaPhiMax(0.012);

    Posv0Daug->SetNClsTPC(80);
    Negv0Daug->SetNClsTPC(80);
    PosAntiv0Daug->SetNClsTPC(80);
    NegAntiv0Daug->SetNClsTPC(80);

    Posv0Daug->SetEtaRange(-0.77, 0.77);
    Negv0Daug->SetEtaRange(-0.77, 0.77);
    PosAntiv0Daug->SetEtaRange(-0.77, 0.77);
    NegAntiv0Daug->SetEtaRange(-0.77, 0.77);

    CascadeCuts->SetCutXiMinDistBachToPrimVtx(0.04);
    AntiCascadeCuts->SetCutXiMinDistBachToPrimVtx(0.04);

    CascadeCuts->SetCutXiTransverseRadius(1.0, 200);
    AntiCascadeCuts->SetCutXiTransverseRadius(1.0, 200);

    CascadeCuts->SetCutv0MaxDaughterDCA(1.4);
    AntiCascadeCuts->SetCutv0MaxDaughterDCA(1.4);

    CascadeCuts->SetCutv0TransverseRadius(1.7, 200);
    AntiCascadeCuts->SetCutv0TransverseRadius(1.7, 200);

    CascadeCuts->SetCutv0MinDaugDistToPrimVtx(0.06);
    AntiCascadeCuts->SetCutv0MinDaugDistToPrimVtx(0.06);

    XiNegCuts->SetEtaRange(-0.77, 0.77);
    XiPosCuts->SetEtaRange(-0.77, 0.77);
    XiBachCuts->SetEtaRange(-0.77, 0.77);
    AntiXiNegCuts->SetEtaRange(-0.77, 0.77);
    AntiXiPosCuts->SetEtaRange(-0.77, 0.77);
    AntiXiBachCuts->SetEtaRange(-0.77, 0.77);

    XiNegCuts->SetPID(AliPID::kPion, 999, 3);
    XiPosCuts->SetPID(AliPID::kProton, 999, 3);
    XiBachCuts->SetPID(AliPID::kPion, 999, 3);
    AntiXiNegCuts->SetPID(AliPID::kProton, 999, 3);
    AntiXiPosCuts->SetPID(AliPID::kPion, 999, 3);
    AntiXiBachCuts->SetPID(AliPID::kPion, 999, 3);
  } else if (suffix == "22") {
    TrackCuts->SetPtRange(0.4, 4.05);
    AntiTrackCuts->SetPtRange(0.4, 4.05);

    TrackCuts->SetEtaRange(-0.77, 0.77);
    AntiTrackCuts->SetEtaRange(-0.77, 0.77);

    TrackCuts->SetNClsTPC(70);
    AntiTrackCuts->SetNClsTPC(70);

    config->SetDeltaEtaMax(0.012);
    config->SetDeltaPhiMax(0.012);

    v0Cuts->SetCutCPA(0.995);
    Antiv0Cuts->SetCutCPA(0.995);

    Posv0Daug->SetPID(AliPID::kProton, 999.9, 4);
    Negv0Daug->SetPID(AliPID::kPion, 999.9, 4);
    PosAntiv0Daug->SetPID(AliPID::kPion, 999.9, 4);
    NegAntiv0Daug->SetPID(AliPID::kProton, 999.9, 4);

    CascadeCuts->SetCutXiDaughterDCA(1.9);
    AntiCascadeCuts->SetCutXiDaughterDCA(1.9);

    CascadeCuts->SetCutXiMinDistBachToPrimVtx(0.04);
    AntiCascadeCuts->SetCutXiMinDistBachToPrimVtx(0.04);

    CascadeCuts->SetCutXiTransverseRadius(0.6, 200);
    AntiCascadeCuts->SetCutXiTransverseRadius(0.6, 200);

    CascadeCuts->SetCutv0MaxDaughterDCA(1.3);
    AntiCascadeCuts->SetCutv0MaxDaughterDCA(1.3);

    CascadeCuts->SetCutv0CPA(0.99);
    AntiCascadeCuts->SetCutv0CPA(0.99);

    CascadeCuts->SetCutv0MinDistToPrimVtx(0.06);
    AntiCascadeCuts->SetCutv0MinDistToPrimVtx(0.06);

    XiNegCuts->SetEtaRange(-0.83, 0.87);
    XiPosCuts->SetEtaRange(-0.83, 0.87);
    XiBachCuts->SetEtaRange(-0.83, 0.87);
    AntiXiNegCuts->SetEtaRange(-0.83, 0.87);
    AntiXiPosCuts->SetEtaRange(-0.83, 0.87);
    AntiXiBachCuts->SetEtaRange(-0.83, 0.87);

    XiNegCuts->SetPID(AliPID::kPion, 999, 3);
    XiPosCuts->SetPID(AliPID::kProton, 999, 3);
    XiBachCuts->SetPID(AliPID::kPion, 999, 3);
    AntiXiNegCuts->SetPID(AliPID::kProton, 999, 3);
    AntiXiPosCuts->SetPID(AliPID::kPion, 999, 3);
    AntiXiBachCuts->SetPID(AliPID::kPion, 999, 3);

    CascadeCuts->SetPtRangeXi(0.8, 999.5);
    AntiCascadeCuts->SetPtRangeXi(0.8, 999.5);
  } else if (suffix == "23") {
    TrackCuts->SetEtaRange(-0.77, 0.77);
    AntiTrackCuts->SetEtaRange(-0.77, 0.77);

    TrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);
    AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);

    Posv0Daug->SetPID(AliPID::kProton, 999.9, 4);
    Negv0Daug->SetPID(AliPID::kPion, 999.9, 4);
    PosAntiv0Daug->SetPID(AliPID::kPion, 999.9, 4);
    NegAntiv0Daug->SetPID(AliPID::kProton, 999.9, 4);

    Posv0Daug->SetEtaRange(-0.83, 0.83);
    Negv0Daug->SetEtaRange(-0.83, 0.83);
    PosAntiv0Daug->SetEtaRange(-0.83, 0.83);
    NegAntiv0Daug->SetEtaRange(-0.83, 0.83);

    v0Cuts->SetCutDCADaugTov0Vtx(1.2);
    Antiv0Cuts->SetCutDCADaugTov0Vtx(1.2);

    CascadeCuts->SetCutXiDaughterDCA(1.3);
    AntiCascadeCuts->SetCutXiDaughterDCA(1.3);

    CascadeCuts->SetCutXiCPA(0.985);
    AntiCascadeCuts->SetCutXiCPA(0.985);

    CascadeCuts->SetCutXiTransverseRadius(0.6, 200);
    AntiCascadeCuts->SetCutXiTransverseRadius(0.6, 200);

    XiNegCuts->SetEtaRange(-0.83, 0.87);
    XiPosCuts->SetEtaRange(-0.83, 0.87);
    XiBachCuts->SetEtaRange(-0.83, 0.87);
    AntiXiNegCuts->SetEtaRange(-0.83, 0.87);
    AntiXiPosCuts->SetEtaRange(-0.83, 0.87);
    AntiXiBachCuts->SetEtaRange(-0.83, 0.87);

    XiNegCuts->SetPID(AliPID::kPion, 999, 4.5);
    XiPosCuts->SetPID(AliPID::kProton, 999, 4.5);
    XiBachCuts->SetPID(AliPID::kPion, 999, 4.5);
    AntiXiNegCuts->SetPID(AliPID::kProton, 999, 4.5);
    AntiXiPosCuts->SetPID(AliPID::kPion, 999, 4.5);
    AntiXiBachCuts->SetPID(AliPID::kPion, 999, 4.5);

    CascadeCuts->SetPtRangeXi(0.8, 999.5);
    AntiCascadeCuts->SetPtRangeXi(0.8, 999.5);

  } else if (suffix == "24") {
    TrackCuts->SetEtaRange(-0.77, 0.77);
    AntiTrackCuts->SetEtaRange(-0.77, 0.77);

    TrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);
    AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);

    TrackCuts->SetNClsTPC(90);
    AntiTrackCuts->SetNClsTPC(90);

    config->SetDeltaEtaMax(0.012);
    config->SetDeltaPhiMax(0.012);

    v0Cuts->SetCutCPA(0.995);
    Antiv0Cuts->SetCutCPA(0.995);

    Posv0Daug->SetPID(AliPID::kProton, 999.9, 4);
    Negv0Daug->SetPID(AliPID::kPion, 999.9, 4);
    PosAntiv0Daug->SetPID(AliPID::kPion, 999.9, 4);
    NegAntiv0Daug->SetPID(AliPID::kProton, 999.9, 4);

    Posv0Daug->SetEtaRange(-0.83, 0.83);
    Negv0Daug->SetEtaRange(-0.83, 0.83);
    PosAntiv0Daug->SetEtaRange(-0.83, 0.83);
    NegAntiv0Daug->SetEtaRange(-0.83, 0.83);

    //XI
    CascadeCuts->SetCutXiDaughterDCA(1.3);
    AntiCascadeCuts->SetCutXiDaughterDCA(1.3);

    CascadeCuts->SetCutXiMinDistBachToPrimVtx(0.06);
    AntiCascadeCuts->SetCutXiMinDistBachToPrimVtx(0.06);

    CascadeCuts->SetCutXiTransverseRadius(1.0, 200);
    AntiCascadeCuts->SetCutXiTransverseRadius(1.0, 200);

    CascadeCuts->SetCutv0CPA(0.99);
    AntiCascadeCuts->SetCutv0CPA(0.99);

    CascadeCuts->SetCutv0TransverseRadius(1.7, 200);
    AntiCascadeCuts->SetCutv0TransverseRadius(1.7, 200);

    CascadeCuts->SetCutv0MinDistToPrimVtx(0.08);
    AntiCascadeCuts->SetCutv0MinDistToPrimVtx(0.08);

    XiNegCuts->SetPID(AliPID::kPion, 999, 4.5);
    XiPosCuts->SetPID(AliPID::kProton, 999, 4.5);
    XiBachCuts->SetPID(AliPID::kPion, 999, 4.5);
    AntiXiNegCuts->SetPID(AliPID::kProton, 999, 4.5);
    AntiXiPosCuts->SetPID(AliPID::kPion, 999, 4.5);
    AntiXiBachCuts->SetPID(AliPID::kPion, 999, 4.5);

  } else if (suffix == "25") {
    TrackCuts->SetPtRange(0.4, 4.05);
    AntiTrackCuts->SetPtRange(0.4, 4.05);

    TrackCuts->SetEtaRange(-0.83, 0.83);
    AntiTrackCuts->SetEtaRange(-0.83, 0.83);

    config->SetDeltaEtaMax(0.012);
    config->SetDeltaPhiMax(0.012);

    Posv0Daug->SetPID(AliPID::kProton, 999.9, 4);
    Negv0Daug->SetPID(AliPID::kPion, 999.9, 4);
    PosAntiv0Daug->SetPID(AliPID::kPion, 999.9, 4);
    NegAntiv0Daug->SetPID(AliPID::kProton, 999.9, 4);

    Posv0Daug->SetEtaRange(-0.77, 0.77);
    Negv0Daug->SetEtaRange(-0.77, 0.77);
    PosAntiv0Daug->SetEtaRange(-0.77, 0.77);
    NegAntiv0Daug->SetEtaRange(-0.77, 0.77);

    v0Cuts->SetCutDCADaugToPrimVtx(0.06);
    Antiv0Cuts->SetCutDCADaugToPrimVtx(0.06);

    CascadeCuts->SetCutXiMinDistBachToPrimVtx(0.06);
    AntiCascadeCuts->SetCutXiMinDistBachToPrimVtx(0.06);

    CascadeCuts->SetCutXiCPA(0.99);
    AntiCascadeCuts->SetCutXiCPA(0.99);

    CascadeCuts->SetCutXiTransverseRadius(0.6, 200);
    AntiCascadeCuts->SetCutXiTransverseRadius(0.6, 200);

    CascadeCuts->SetCutv0MaxDaughterDCA(1.3);
    AntiCascadeCuts->SetCutv0MaxDaughterDCA(1.3);

    CascadeCuts->SetCutv0TransverseRadius(1.1, 200);
    AntiCascadeCuts->SetCutv0TransverseRadius(1.1, 200);

    CascadeCuts->SetCutv0MinDaugDistToPrimVtx(0.06);
    AntiCascadeCuts->SetCutv0MinDaugDistToPrimVtx(0.06);

    XiNegCuts->SetEtaRange(-0.83, 0.87);
    XiPosCuts->SetEtaRange(-0.83, 0.87);
    XiBachCuts->SetEtaRange(-0.83, 0.87);
    AntiXiNegCuts->SetEtaRange(-0.83, 0.87);
    AntiXiPosCuts->SetEtaRange(-0.83, 0.87);
    AntiXiBachCuts->SetEtaRange(-0.83, 0.87);

    XiNegCuts->SetPID(AliPID::kPion, 999, 4.5);
    XiPosCuts->SetPID(AliPID::kProton, 999, 4.5);
    XiBachCuts->SetPID(AliPID::kPion, 999, 4.5);
    AntiXiNegCuts->SetPID(AliPID::kProton, 999, 4.5);
    AntiXiPosCuts->SetPID(AliPID::kPion, 999, 4.5);
    AntiXiBachCuts->SetPID(AliPID::kPion, 999, 4.5);

  } else if (suffix == "26") {
    TrackCuts->SetPtRange(0.6, 4.05);
    AntiTrackCuts->SetPtRange(0.6, 4.05);

    TrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);
    AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);

    config->SetDeltaEtaMax(0.012);
    config->SetDeltaPhiMax(0.012);

    Posv0Daug->SetEtaRange(-0.77, 0.77);
    Negv0Daug->SetEtaRange(-0.77, 0.77);
    PosAntiv0Daug->SetEtaRange(-0.77, 0.77);
    NegAntiv0Daug->SetEtaRange(-0.77, 0.77);

    CascadeCuts->SetCutXiDaughterDCA(1.3);
    AntiCascadeCuts->SetCutXiDaughterDCA(1.3);

    CascadeCuts->SetCutXiMinDistBachToPrimVtx(0.04);
    AntiCascadeCuts->SetCutXiMinDistBachToPrimVtx(0.04);

    CascadeCuts->SetCutXiCPA(0.99);
    AntiCascadeCuts->SetCutXiCPA(0.99);

    CascadeCuts->SetCutXiTransverseRadius(0.6, 200);
    AntiCascadeCuts->SetCutXiTransverseRadius(0.6, 200);

    CascadeCuts->SetCutv0MaxDaughterDCA(1.3);
    AntiCascadeCuts->SetCutv0MaxDaughterDCA(1.3);

    CascadeCuts->SetCutv0CPA(0.96);
    AntiCascadeCuts->SetCutv0CPA(0.96);

    CascadeCuts->SetCutv0TransverseRadius(1.7, 200);
    AntiCascadeCuts->SetCutv0TransverseRadius(1.7, 200);

    CascadeCuts->SetCutv0MinDistToPrimVtx(0.06);
    AntiCascadeCuts->SetCutv0MinDistToPrimVtx(0.06);

    XiNegCuts->SetEtaRange(-0.77, 0.77);
    XiPosCuts->SetEtaRange(-0.77, 0.77);
    XiBachCuts->SetEtaRange(-0.77, 0.77);
    AntiXiNegCuts->SetEtaRange(-0.77, 0.77);
    AntiXiPosCuts->SetEtaRange(-0.77, 0.77);
    AntiXiBachCuts->SetEtaRange(-0.77, 0.77);

    XiNegCuts->SetPID(AliPID::kPion, 999, 3.5);
    XiPosCuts->SetPID(AliPID::kProton, 999, 3.5);
    XiBachCuts->SetPID(AliPID::kPion, 999, 3.5);
    AntiXiNegCuts->SetPID(AliPID::kProton, 999, 3.5);
    AntiXiPosCuts->SetPID(AliPID::kPion, 999, 3.5);
    AntiXiBachCuts->SetPID(AliPID::kPion, 999, 3.5);

  } else if (suffix == "27") {
    TrackCuts->SetEtaRange(-0.77, 0.77);
    AntiTrackCuts->SetEtaRange(-0.77, 0.77);

    TrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);
    AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);

    TrackCuts->SetNClsTPC(90);
    AntiTrackCuts->SetNClsTPC(90);

    Posv0Daug->SetPID(AliPID::kProton, 999.9, 4);
    Negv0Daug->SetPID(AliPID::kPion, 999.9, 4);
    PosAntiv0Daug->SetPID(AliPID::kPion, 999.9, 4);
    NegAntiv0Daug->SetPID(AliPID::kProton, 999.9, 4);

    Posv0Daug->SetEtaRange(-0.77, 0.77);
    Negv0Daug->SetEtaRange(-0.77, 0.77);
    PosAntiv0Daug->SetEtaRange(-0.77, 0.77);
    NegAntiv0Daug->SetEtaRange(-0.77, 0.77);

    v0Cuts->SetCutDCADaugTov0Vtx(1.2);
    Antiv0Cuts->SetCutDCADaugTov0Vtx(1.2);

    CascadeCuts->SetCutXiDaughterDCA(1.9);
    AntiCascadeCuts->SetCutXiDaughterDCA(1.9);

    CascadeCuts->SetCutXiMinDistBachToPrimVtx(0.06);
    AntiCascadeCuts->SetCutXiMinDistBachToPrimVtx(0.06);

    CascadeCuts->SetCutv0TransverseRadius(1.1, 200);
    AntiCascadeCuts->SetCutv0TransverseRadius(1.1, 200);

    CascadeCuts->SetCutv0MinDistToPrimVtx(0.06);
    AntiCascadeCuts->SetCutv0MinDistToPrimVtx(0.06);

    CascadeCuts->SetCutv0MinDaugDistToPrimVtx(0.06);
    AntiCascadeCuts->SetCutv0MinDaugDistToPrimVtx(0.06);

    XiNegCuts->SetEtaRange(-0.83, 0.87);
    XiPosCuts->SetEtaRange(-0.83, 0.87);
    XiBachCuts->SetEtaRange(-0.83, 0.87);
    AntiXiNegCuts->SetEtaRange(-0.83, 0.87);
    AntiXiPosCuts->SetEtaRange(-0.83, 0.87);
    AntiXiBachCuts->SetEtaRange(-0.83, 0.87);

    XiNegCuts->SetPID(AliPID::kPion, 999, 4.5);
    XiPosCuts->SetPID(AliPID::kProton, 999, 4.5);
    XiBachCuts->SetPID(AliPID::kPion, 999, 4.5);
    AntiXiNegCuts->SetPID(AliPID::kProton, 999, 4.5);
    AntiXiPosCuts->SetPID(AliPID::kPion, 999, 4.5);
    AntiXiBachCuts->SetPID(AliPID::kPion, 999, 4.5);

    CascadeCuts->SetPtRangeXi(0.8, 999.5);
    AntiCascadeCuts->SetPtRangeXi(0.8, 999.5);

  } else if (suffix == "28") {
    TrackCuts->SetNClsTPC(90);
    AntiTrackCuts->SetNClsTPC(90);

    v0Cuts->SetCutCPA(0.995);
    Antiv0Cuts->SetCutCPA(0.995);

    v0Cuts->SetCutDCADaugToPrimVtx(0.06);
    Antiv0Cuts->SetCutDCADaugToPrimVtx(0.06);

    CascadeCuts->SetCutXiDaughterDCA(1.9);
    AntiCascadeCuts->SetCutXiDaughterDCA(1.9);

    CascadeCuts->SetCutXiMinDistBachToPrimVtx(0.04);
    AntiCascadeCuts->SetCutXiMinDistBachToPrimVtx(0.04);

    CascadeCuts->SetCutv0CPA(0.99);
    AntiCascadeCuts->SetCutv0CPA(0.99);

    CascadeCuts->SetCutv0TransverseRadius(1.1, 200);
    AntiCascadeCuts->SetCutv0TransverseRadius(1.1, 200);

    CascadeCuts->SetCutv0MinDistToPrimVtx(0.08);
    AntiCascadeCuts->SetCutv0MinDistToPrimVtx(0.08);

    CascadeCuts->SetPtRangeXi(1.2, 999.5);
    AntiCascadeCuts->SetPtRangeXi(1.2, 999.5);

  } else if (suffix == "29") {
    TrackCuts->SetPtRange(0.4, 4.05);
    AntiTrackCuts->SetPtRange(0.4, 4.05);

    TrackCuts->SetEtaRange(-0.77, 0.77);
    AntiTrackCuts->SetEtaRange(-0.77, 0.77);

    TrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);
    AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);

    config->SetDeltaEtaMax(0.012);
    config->SetDeltaPhiMax(0.012);

    Posv0Daug->SetNClsTPC(80);
    Negv0Daug->SetNClsTPC(80);
    PosAntiv0Daug->SetNClsTPC(80);
    NegAntiv0Daug->SetNClsTPC(80);

    v0Cuts->SetCutDCADaugToPrimVtx(0.06);
    Antiv0Cuts->SetCutDCADaugToPrimVtx(0.06);

    CascadeCuts->SetCutXiDaughterDCA(1.3);
    AntiCascadeCuts->SetCutXiDaughterDCA(1.3);

    CascadeCuts->SetCutXiMinDistBachToPrimVtx(0.04);
    AntiCascadeCuts->SetCutXiMinDistBachToPrimVtx(0.04);

    CascadeCuts->SetCutXiCPA(0.99);
    AntiCascadeCuts->SetCutXiCPA(0.99);

    CascadeCuts->SetCutXiTransverseRadius(1.0, 200);
    AntiCascadeCuts->SetCutXiTransverseRadius(1.0, 200);

    CascadeCuts->SetCutv0MaxDaughterDCA(1.4);
    AntiCascadeCuts->SetCutv0MaxDaughterDCA(1.4);

    CascadeCuts->SetCutv0CPA(0.96);
    AntiCascadeCuts->SetCutv0CPA(0.96);

    CascadeCuts->SetCutv0TransverseRadius(1.7, 200);
    AntiCascadeCuts->SetCutv0TransverseRadius(1.7, 200);

    CascadeCuts->SetCutv0MinDistToPrimVtx(0.06);
    AntiCascadeCuts->SetCutv0MinDistToPrimVtx(0.06);

    XiNegCuts->SetEtaRange(-0.83, 0.87);
    XiPosCuts->SetEtaRange(-0.83, 0.87);
    XiBachCuts->SetEtaRange(-0.83, 0.87);
    AntiXiNegCuts->SetEtaRange(-0.83, 0.87);
    AntiXiPosCuts->SetEtaRange(-0.83, 0.87);
    AntiXiBachCuts->SetEtaRange(-0.83, 0.87);

    XiNegCuts->SetPID(AliPID::kPion, 999, 3);
    XiPosCuts->SetPID(AliPID::kProton, 999, 3);
    XiBachCuts->SetPID(AliPID::kPion, 999, 3);
    AntiXiNegCuts->SetPID(AliPID::kProton, 999, 3);
    AntiXiPosCuts->SetPID(AliPID::kPion, 999, 3);
    AntiXiBachCuts->SetPID(AliPID::kPion, 999, 3);

    CascadeCuts->SetPtRangeXi(0.8, 999.5);
    AntiCascadeCuts->SetPtRangeXi(0.8, 999.5);
  } else if (suffix == "30") {
    TrackCuts->SetEtaRange(-0.85, 0.85);
    AntiTrackCuts->SetEtaRange(-0.85, 0.85);

    config->SetDeltaEtaMax(0.012);
    config->SetDeltaPhiMax(0.012);

    v0Cuts->SetCutCPA(0.995);
    Antiv0Cuts->SetCutCPA(0.995);

    Posv0Daug->SetPID(AliPID::kProton, 999.9, 4);
    Negv0Daug->SetPID(AliPID::kPion, 999.9, 4);
    PosAntiv0Daug->SetPID(AliPID::kPion, 999.9, 4);
    NegAntiv0Daug->SetPID(AliPID::kProton, 999.9, 4);

    Posv0Daug->SetNClsTPC(80);
    Negv0Daug->SetNClsTPC(80);
    PosAntiv0Daug->SetNClsTPC(80);
    NegAntiv0Daug->SetNClsTPC(80);

    Posv0Daug->SetEtaRange(-0.83, 0.83);
    Negv0Daug->SetEtaRange(-0.83, 0.83);
    PosAntiv0Daug->SetEtaRange(-0.83, 0.83);
    NegAntiv0Daug->SetEtaRange(-0.83, 0.83);

    v0Cuts->SetCutDCADaugTov0Vtx(1.2);
    Antiv0Cuts->SetCutDCADaugTov0Vtx(1.2);

    CascadeCuts->SetCutXiDaughterDCA(1.3);
    AntiCascadeCuts->SetCutXiDaughterDCA(1.3);

    CascadeCuts->SetCutXiMinDistBachToPrimVtx(0.04);
    AntiCascadeCuts->SetCutXiMinDistBachToPrimVtx(0.04);

    CascadeCuts->SetCutXiCPA(0.99);
    AntiCascadeCuts->SetCutXiCPA(0.99);

    CascadeCuts->SetCutv0MaxDaughterDCA(1.3);
    AntiCascadeCuts->SetCutv0MaxDaughterDCA(1.3);

    CascadeCuts->SetCutv0CPA(0.96);
    AntiCascadeCuts->SetCutv0CPA(0.96);

    CascadeCuts->SetCutv0TransverseRadius(1.7, 200);
    AntiCascadeCuts->SetCutv0TransverseRadius(1.7, 200);

    CascadeCuts->SetCutv0MinDistToPrimVtx(0.06);
    AntiCascadeCuts->SetCutv0MinDistToPrimVtx(0.06);

    CascadeCuts->SetCutv0MinDaugDistToPrimVtx(0.06);
    AntiCascadeCuts->SetCutv0MinDaugDistToPrimVtx(0.06);

    XiNegCuts->SetEtaRange(-0.83, 0.83);
    XiPosCuts->SetEtaRange(-0.83, 0.83);
    XiBachCuts->SetEtaRange(-0.83, 0.83);
    AntiXiNegCuts->SetEtaRange(-0.83, 0.83);
    AntiXiPosCuts->SetEtaRange(-0.83, 0.83);
    AntiXiBachCuts->SetEtaRange(-0.83, 0.83);

    XiNegCuts->SetPID(AliPID::kPion, 999, 3.5);
    XiPosCuts->SetPID(AliPID::kProton, 999, 3.5);
    XiBachCuts->SetPID(AliPID::kPion, 999, 3.5);
    AntiXiNegCuts->SetPID(AliPID::kProton, 999, 3.5);
    AntiXiPosCuts->SetPID(AliPID::kPion, 999, 3.5);
    AntiXiBachCuts->SetPID(AliPID::kPion, 999, 3.5);

    CascadeCuts->SetPtRangeXi(1.2, 999.5);
    AntiCascadeCuts->SetPtRangeXi(1.2, 999.5);

  } else if (suffix == "31") {
    TrackCuts->SetPtRange(0.4, 4.05);
    AntiTrackCuts->SetPtRange(0.4, 4.05);

    TrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);
    AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);

    TrackCuts->SetNClsTPC(70);
    AntiTrackCuts->SetNClsTPC(70);

    config->SetDeltaEtaMax(0.012);
    config->SetDeltaPhiMax(0.012);

    v0Cuts->SetCutCPA(0.995);
    Antiv0Cuts->SetCutCPA(0.995);

    Posv0Daug->SetEtaRange(-0.83, 0.83);
    Negv0Daug->SetEtaRange(-0.83, 0.83);
    PosAntiv0Daug->SetEtaRange(-0.83, 0.83);
    NegAntiv0Daug->SetEtaRange(-0.83, 0.83);

    v0Cuts->SetCutDCADaugTov0Vtx(1.2);
    Antiv0Cuts->SetCutDCADaugTov0Vtx(1.2);

    CascadeCuts->SetCutXiDaughterDCA(1.9);
    AntiCascadeCuts->SetCutXiDaughterDCA(1.9);

    CascadeCuts->SetCutXiCPA(0.985);
    AntiCascadeCuts->SetCutXiCPA(0.985);

    CascadeCuts->SetCutv0MaxDaughterDCA(1.3);
    AntiCascadeCuts->SetCutv0MaxDaughterDCA(1.3);

    CascadeCuts->SetCutv0CPA(0.99);
    AntiCascadeCuts->SetCutv0CPA(0.99);

    CascadeCuts->SetCutv0TransverseRadius(1.1, 200);
    AntiCascadeCuts->SetCutv0TransverseRadius(1.1, 200);

    CascadeCuts->SetCutv0MinDaugDistToPrimVtx(0.06);
    AntiCascadeCuts->SetCutv0MinDaugDistToPrimVtx(0.06);

    XiNegCuts->SetPID(AliPID::kPion, 999, 3);
    XiPosCuts->SetPID(AliPID::kProton, 999, 3);
    XiBachCuts->SetPID(AliPID::kPion, 999, 3);
    AntiXiNegCuts->SetPID(AliPID::kProton, 999, 3);
    AntiXiPosCuts->SetPID(AliPID::kPion, 999, 3);
    AntiXiBachCuts->SetPID(AliPID::kPion, 999, 3);

  } else if (suffix == "32") {

    TrackCuts->SetEtaRange(-0.77, 0.77);
    AntiTrackCuts->SetEtaRange(-0.77, 0.77);

    TrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);
    AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);

    TrackCuts->SetNClsTPC(90);
    AntiTrackCuts->SetNClsTPC(90);

    v0Cuts->SetCutCPA(0.995);
    Antiv0Cuts->SetCutCPA(0.995);

    Posv0Daug->SetNClsTPC(80);
    Negv0Daug->SetNClsTPC(80);
    PosAntiv0Daug->SetNClsTPC(80);
    NegAntiv0Daug->SetNClsTPC(80);

    v0Cuts->SetCutDCADaugToPrimVtx(0.06);
    Antiv0Cuts->SetCutDCADaugToPrimVtx(0.06);

    CascadeCuts->SetCutv0MaxDaughterDCA(1.3);
    AntiCascadeCuts->SetCutv0MaxDaughterDCA(1.3);

    CascadeCuts->SetCutv0CPA(0.99);
    AntiCascadeCuts->SetCutv0CPA(0.99);

    CascadeCuts->SetCutv0TransverseRadius(1.1, 200);
    AntiCascadeCuts->SetCutv0TransverseRadius(1.1, 200);

    CascadeCuts->SetCutv0MinDistToPrimVtx(0.08);
    AntiCascadeCuts->SetCutv0MinDistToPrimVtx(0.08);

    CascadeCuts->SetCutv0MinDaugDistToPrimVtx(0.06);
    AntiCascadeCuts->SetCutv0MinDaugDistToPrimVtx(0.06);

    XiNegCuts->SetPID(AliPID::kPion, 999, 4.5);
    XiPosCuts->SetPID(AliPID::kProton, 999, 4.5);
    XiBachCuts->SetPID(AliPID::kPion, 999, 4.5);
    AntiXiNegCuts->SetPID(AliPID::kProton, 999, 4.5);
    AntiXiPosCuts->SetPID(AliPID::kPion, 999, 4.5);
    AntiXiBachCuts->SetPID(AliPID::kPion, 999, 4.5);

    CascadeCuts->SetPtRangeXi(0.8, 999.5);
    AntiCascadeCuts->SetPtRangeXi(0.8, 999.5);
  } else if (suffix == "33") {
    TrackCuts->SetEtaRange(-0.77, 0.77);
    AntiTrackCuts->SetEtaRange(-0.77, 0.77);

    config->SetDeltaEtaMax(0.012);
    config->SetDeltaPhiMax(0.012);

    v0Cuts->SetCutCPA(0.995);
    Antiv0Cuts->SetCutCPA(0.995);

    Posv0Daug->SetNClsTPC(80);
    Negv0Daug->SetNClsTPC(80);
    PosAntiv0Daug->SetNClsTPC(80);
    NegAntiv0Daug->SetNClsTPC(80);

    Posv0Daug->SetEtaRange(-0.77, 0.77);
    Negv0Daug->SetEtaRange(-0.77, 0.77);
    PosAntiv0Daug->SetEtaRange(-0.77, 0.77);
    NegAntiv0Daug->SetEtaRange(-0.77, 0.77);

    CascadeCuts->SetCutXiCPA(0.985);
    AntiCascadeCuts->SetCutXiCPA(0.985);

    CascadeCuts->SetCutXiTransverseRadius(1.0, 200);
    AntiCascadeCuts->SetCutXiTransverseRadius(1.0, 200);

    CascadeCuts->SetCutv0CPA(0.96);
    AntiCascadeCuts->SetCutv0CPA(0.96);

    CascadeCuts->SetCutv0MinDaugDistToPrimVtx(0.06);
    AntiCascadeCuts->SetCutv0MinDaugDistToPrimVtx(0.06);

    XiNegCuts->SetEtaRange(-0.77, 0.77);
    XiPosCuts->SetEtaRange(-0.77, 0.77);
    XiBachCuts->SetEtaRange(-0.77, 0.77);
    AntiXiNegCuts->SetEtaRange(-0.77, 0.77);
    AntiXiPosCuts->SetEtaRange(-0.77, 0.77);
    AntiXiBachCuts->SetEtaRange(-0.77, 0.77);

    XiNegCuts->SetPID(AliPID::kPion, 999, 3.5);
    XiPosCuts->SetPID(AliPID::kProton, 999, 3.5);
    XiBachCuts->SetPID(AliPID::kPion, 999, 3.5);
    AntiXiNegCuts->SetPID(AliPID::kProton, 999, 3.5);
    AntiXiPosCuts->SetPID(AliPID::kPion, 999, 3.5);
    AntiXiBachCuts->SetPID(AliPID::kPion, 999, 3.5);
  } else if (suffix == "34") {
    TrackCuts->SetEtaRange(-0.77, 0.77);
    AntiTrackCuts->SetEtaRange(-0.77, 0.77);

    TrackCuts->SetNClsTPC(70);
    AntiTrackCuts->SetNClsTPC(70);

    config->SetDeltaEtaMax(0.012);
    config->SetDeltaPhiMax(0.012);

    Posv0Daug->SetNClsTPC(80);
    Negv0Daug->SetNClsTPC(80);
    PosAntiv0Daug->SetNClsTPC(80);
    NegAntiv0Daug->SetNClsTPC(80);

    Posv0Daug->SetEtaRange(-0.77, 0.77);
    Negv0Daug->SetEtaRange(-0.77, 0.77);
    PosAntiv0Daug->SetEtaRange(-0.77, 0.77);
    NegAntiv0Daug->SetEtaRange(-0.77, 0.77);

    CascadeCuts->SetCutXiDaughterDCA(1.3);
    AntiCascadeCuts->SetCutXiDaughterDCA(1.3);

    CascadeCuts->SetCutXiCPA(0.985);
    AntiCascadeCuts->SetCutXiCPA(0.985);

    CascadeCuts->SetCutXiTransverseRadius(0.6, 200);
    AntiCascadeCuts->SetCutXiTransverseRadius(0.6, 200);

    CascadeCuts->SetCutv0MaxDaughterDCA(1.4);
    AntiCascadeCuts->SetCutv0MaxDaughterDCA(1.4);

    CascadeCuts->SetCutv0TransverseRadius(1.1, 200);
    AntiCascadeCuts->SetCutv0TransverseRadius(1.1, 200);

    XiNegCuts->SetEtaRange(-0.83, 0.83);
    XiPosCuts->SetEtaRange(-0.83, 0.83);
    XiBachCuts->SetEtaRange(-0.83, 0.83);
    AntiXiNegCuts->SetEtaRange(-0.83, 0.83);
    AntiXiPosCuts->SetEtaRange(-0.83, 0.83);
    AntiXiBachCuts->SetEtaRange(-0.83, 0.83);

    CascadeCuts->SetPtRangeXi(0.8, 999.5);
    AntiCascadeCuts->SetPtRangeXi(0.8, 999.5);

  } else if (suffix == "35") {
    TrackCuts->SetPtRange(0.4, 4.05);
    AntiTrackCuts->SetPtRange(0.4, 4.05);

    TrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);
    AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);

    TrackCuts->SetNClsTPC(70);
    AntiTrackCuts->SetNClsTPC(70);

    config->SetDeltaEtaMax(0.012);
    config->SetDeltaPhiMax(0.012);

    v0Cuts->SetCutCPA(0.995);
    Antiv0Cuts->SetCutCPA(0.995);

    Posv0Daug->SetNClsTPC(80);
    Negv0Daug->SetNClsTPC(80);
    PosAntiv0Daug->SetNClsTPC(80);
    NegAntiv0Daug->SetNClsTPC(80);

    v0Cuts->SetCutDCADaugTov0Vtx(1.2);
    Antiv0Cuts->SetCutDCADaugTov0Vtx(1.2);

    CascadeCuts->SetCutXiDaughterDCA(1.9);
    AntiCascadeCuts->SetCutXiDaughterDCA(1.9);

    CascadeCuts->SetCutXiTransverseRadius(0.6, 200);
    AntiCascadeCuts->SetCutXiTransverseRadius(0.6, 200);

    CascadeCuts->SetCutv0MaxDaughterDCA(1.4);
    AntiCascadeCuts->SetCutv0MaxDaughterDCA(1.4);

    CascadeCuts->SetCutv0MinDistToPrimVtx(0.06);
    AntiCascadeCuts->SetCutv0MinDistToPrimVtx(0.06);

    XiNegCuts->SetEtaRange(-0.83, 0.83);
    XiPosCuts->SetEtaRange(-0.83, 0.83);
    XiBachCuts->SetEtaRange(-0.83, 0.83);
    AntiXiNegCuts->SetEtaRange(-0.83, 0.83);
    AntiXiPosCuts->SetEtaRange(-0.83, 0.83);
    AntiXiBachCuts->SetEtaRange(-0.83, 0.83);
  } else if (suffix == "36") {
    TrackCuts->SetPtRange(0.6, 4.05);
    AntiTrackCuts->SetPtRange(0.6, 4.05);

    TrackCuts->SetEtaRange(-0.85, 0.85);
    AntiTrackCuts->SetEtaRange(-0.85, 0.85);

    TrackCuts->SetNClsTPC(70);
    AntiTrackCuts->SetNClsTPC(70);

    config->SetDeltaEtaMax(0.012);
    config->SetDeltaPhiMax(0.012);

    Posv0Daug->SetPID(AliPID::kProton, 999.9, 4);
    Negv0Daug->SetPID(AliPID::kPion, 999.9, 4);
    PosAntiv0Daug->SetPID(AliPID::kPion, 999.9, 4);
    NegAntiv0Daug->SetPID(AliPID::kProton, 999.9, 4);

    Posv0Daug->SetNClsTPC(80);
    Negv0Daug->SetNClsTPC(80);
    PosAntiv0Daug->SetNClsTPC(80);
    NegAntiv0Daug->SetNClsTPC(80);

    Posv0Daug->SetEtaRange(-0.83, 0.83);
    Negv0Daug->SetEtaRange(-0.83, 0.83);
    PosAntiv0Daug->SetEtaRange(-0.83, 0.83);
    NegAntiv0Daug->SetEtaRange(-0.83, 0.83);

    v0Cuts->SetCutDCADaugTov0Vtx(1.2);
    Antiv0Cuts->SetCutDCADaugTov0Vtx(1.2);

    CascadeCuts->SetCutXiMinDistBachToPrimVtx(0.04);
    AntiCascadeCuts->SetCutXiMinDistBachToPrimVtx(0.04);

    CascadeCuts->SetCutXiTransverseRadius(1.0, 200);
    AntiCascadeCuts->SetCutXiTransverseRadius(1.0, 200);

    CascadeCuts->SetCutv0CPA(0.96);
    AntiCascadeCuts->SetCutv0CPA(0.96);

    CascadeCuts->SetCutv0TransverseRadius(1.7, 200);
    AntiCascadeCuts->SetCutv0TransverseRadius(1.7, 200);

    CascadeCuts->SetCutv0MinDistToPrimVtx(0.06);
    AntiCascadeCuts->SetCutv0MinDistToPrimVtx(0.06);

    XiNegCuts->SetEtaRange(-0.77, 0.77);
    XiPosCuts->SetEtaRange(-0.77, 0.77);
    XiBachCuts->SetEtaRange(-0.77, 0.77);
    AntiXiNegCuts->SetEtaRange(-0.77, 0.77);
    AntiXiPosCuts->SetEtaRange(-0.77, 0.77);
    AntiXiBachCuts->SetEtaRange(-0.77, 0.77);

    XiNegCuts->SetPID(AliPID::kPion, 999, 4.5);
    XiPosCuts->SetPID(AliPID::kProton, 999, 4.5);
    XiBachCuts->SetPID(AliPID::kPion, 999, 4.5);
    AntiXiNegCuts->SetPID(AliPID::kProton, 999, 4.5);
    AntiXiPosCuts->SetPID(AliPID::kPion, 999, 4.5);
    AntiXiBachCuts->SetPID(AliPID::kPion, 999, 4.5);

    CascadeCuts->SetPtRangeXi(1.2, 999.5);
    AntiCascadeCuts->SetPtRangeXi(1.2, 999.5);

  } else if (suffix == "37") {
    TrackCuts->SetEtaRange(-0.85, 0.85);
    AntiTrackCuts->SetEtaRange(-0.85, 0.85);

    TrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);
    AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);

    TrackCuts->SetNClsTPC(90);
    AntiTrackCuts->SetNClsTPC(90);

    config->SetDeltaEtaMax(0.012);
    config->SetDeltaPhiMax(0.012);

    v0Cuts->SetCutCPA(0.995);
    Antiv0Cuts->SetCutCPA(0.995);

    Posv0Daug->SetEtaRange(-0.83, 0.83);
    Negv0Daug->SetEtaRange(-0.83, 0.83);
    PosAntiv0Daug->SetEtaRange(-0.83, 0.83);
    NegAntiv0Daug->SetEtaRange(-0.83, 0.83);

    v0Cuts->SetCutDCADaugTov0Vtx(1.2);
    Antiv0Cuts->SetCutDCADaugTov0Vtx(1.2);

    CascadeCuts->SetCutXiMinDistBachToPrimVtx(0.06);
    AntiCascadeCuts->SetCutXiMinDistBachToPrimVtx(0.06);

    CascadeCuts->SetCutXiCPA(0.99);
    AntiCascadeCuts->SetCutXiCPA(0.99);

    CascadeCuts->SetCutv0CPA(0.96);
    AntiCascadeCuts->SetCutv0CPA(0.96);

    CascadeCuts->SetCutv0TransverseRadius(1.7, 200);
    AntiCascadeCuts->SetCutv0TransverseRadius(1.7, 200);

    CascadeCuts->SetCutv0MinDaugDistToPrimVtx(0.06);
    AntiCascadeCuts->SetCutv0MinDaugDistToPrimVtx(0.06);

    XiNegCuts->SetPID(AliPID::kPion, 999, 3);
    XiPosCuts->SetPID(AliPID::kProton, 999, 3);
    XiBachCuts->SetPID(AliPID::kPion, 999, 3);
    AntiXiNegCuts->SetPID(AliPID::kProton, 999, 3);
    AntiXiPosCuts->SetPID(AliPID::kPion, 999, 3);
    AntiXiBachCuts->SetPID(AliPID::kPion, 999, 3);

  } else if (suffix == "38") {
    TrackCuts->SetEtaRange(-0.77, 0.77);
    AntiTrackCuts->SetEtaRange(-0.77, 0.77);

    TrackCuts->SetNClsTPC(90);
    AntiTrackCuts->SetNClsTPC(90);

    config->SetDeltaEtaMax(0.012);
    config->SetDeltaPhiMax(0.012);

    v0Cuts->SetCutCPA(0.995);
    Antiv0Cuts->SetCutCPA(0.995);

    Posv0Daug->SetNClsTPC(80);
    Negv0Daug->SetNClsTPC(80);
    PosAntiv0Daug->SetNClsTPC(80);
    NegAntiv0Daug->SetNClsTPC(80);

    Posv0Daug->SetEtaRange(-0.77, 0.77);
    Negv0Daug->SetEtaRange(-0.77, 0.77);
    PosAntiv0Daug->SetEtaRange(-0.77, 0.77);
    NegAntiv0Daug->SetEtaRange(-0.77, 0.77);

    CascadeCuts->SetCutXiDaughterDCA(1.9);
    AntiCascadeCuts->SetCutXiDaughterDCA(1.9);

    CascadeCuts->SetCutXiTransverseRadius(0.6, 200);
    AntiCascadeCuts->SetCutXiTransverseRadius(0.6, 200);

    CascadeCuts->SetCutv0MaxDaughterDCA(1.3);
    AntiCascadeCuts->SetCutv0MaxDaughterDCA(1.3);

    CascadeCuts->SetCutv0CPA(0.96);
    AntiCascadeCuts->SetCutv0CPA(0.96);

    CascadeCuts->SetCutv0MinDaugDistToPrimVtx(0.06);
    AntiCascadeCuts->SetCutv0MinDaugDistToPrimVtx(0.06);

    XiNegCuts->SetEtaRange(-0.77, 0.77);
    XiPosCuts->SetEtaRange(-0.77, 0.77);
    XiBachCuts->SetEtaRange(-0.77, 0.77);
    AntiXiNegCuts->SetEtaRange(-0.77, 0.77);
    AntiXiPosCuts->SetEtaRange(-0.77, 0.77);
    AntiXiBachCuts->SetEtaRange(-0.77, 0.77);

    CascadeCuts->SetPtRangeXi(1.2, 999.5);
    AntiCascadeCuts->SetPtRangeXi(1.2, 999.5);

  } else if (suffix == "39") {
    TrackCuts->SetPtRange(0.4, 4.05);
    AntiTrackCuts->SetPtRange(0.4, 4.05);

    TrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);
    AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);

    TrackCuts->SetNClsTPC(90);
    AntiTrackCuts->SetNClsTPC(90);

    config->SetDeltaEtaMax(0.012);
    config->SetDeltaPhiMax(0.012);

    Posv0Daug->SetPID(AliPID::kProton, 999.9, 4);
    Negv0Daug->SetPID(AliPID::kPion, 999.9, 4);
    PosAntiv0Daug->SetPID(AliPID::kPion, 999.9, 4);
    NegAntiv0Daug->SetPID(AliPID::kProton, 999.9, 4);

    Posv0Daug->SetNClsTPC(80);
    Negv0Daug->SetNClsTPC(80);
    PosAntiv0Daug->SetNClsTPC(80);
    NegAntiv0Daug->SetNClsTPC(80);

    Posv0Daug->SetEtaRange(-0.77, 0.77);
    Negv0Daug->SetEtaRange(-0.77, 0.77);
    PosAntiv0Daug->SetEtaRange(-0.77, 0.77);
    NegAntiv0Daug->SetEtaRange(-0.77, 0.77);

    v0Cuts->SetCutDCADaugTov0Vtx(1.2);
    Antiv0Cuts->SetCutDCADaugTov0Vtx(1.2);

    CascadeCuts->SetCutXiMinDistBachToPrimVtx(0.04);
    AntiCascadeCuts->SetCutXiMinDistBachToPrimVtx(0.04);

    CascadeCuts->SetCutXiCPA(0.985);
    AntiCascadeCuts->SetCutXiCPA(0.985);

    CascadeCuts->SetCutXiTransverseRadius(1.0, 200);
    AntiCascadeCuts->SetCutXiTransverseRadius(1.0, 200);

    XiNegCuts->SetEtaRange(-0.77, 0.77);
    XiPosCuts->SetEtaRange(-0.77, 0.77);
    XiBachCuts->SetEtaRange(-0.77, 0.77);
    AntiXiNegCuts->SetEtaRange(-0.77, 0.77);
    AntiXiPosCuts->SetEtaRange(-0.77, 0.77);
    AntiXiBachCuts->SetEtaRange(-0.77, 0.77);

    XiNegCuts->SetPID(AliPID::kPion, 999, 4.5);
    XiPosCuts->SetPID(AliPID::kProton, 999, 4.5);
    XiBachCuts->SetPID(AliPID::kPion, 999, 4.5);
    AntiXiNegCuts->SetPID(AliPID::kProton, 999, 4.5);
    AntiXiPosCuts->SetPID(AliPID::kPion, 999, 4.5);
    AntiXiBachCuts->SetPID(AliPID::kPion, 999, 4.5);

    CascadeCuts->SetPtRangeXi(0.8, 999.5);
    AntiCascadeCuts->SetPtRangeXi(0.8, 999.5);
  } else if (suffix == "40") {
    TrackCuts->SetPtRange(0.6, 4.05);
    AntiTrackCuts->SetPtRange(0.6, 4.05);

    TrackCuts->SetNClsTPC(70);
    AntiTrackCuts->SetNClsTPC(70);

    Posv0Daug->SetPID(AliPID::kProton, 999.9, 4);
    Negv0Daug->SetPID(AliPID::kPion, 999.9, 4);
    PosAntiv0Daug->SetPID(AliPID::kPion, 999.9, 4);
    NegAntiv0Daug->SetPID(AliPID::kProton, 999.9, 4);

    v0Cuts->SetCutDCADaugTov0Vtx(1.2);
    Antiv0Cuts->SetCutDCADaugTov0Vtx(1.2);

    CascadeCuts->SetCutXiMinDistBachToPrimVtx(0.06);
    AntiCascadeCuts->SetCutXiMinDistBachToPrimVtx(0.06);

    CascadeCuts->SetCutXiTransverseRadius(1.0, 200);
    AntiCascadeCuts->SetCutXiTransverseRadius(1.0, 200);

    CascadeCuts->SetCutv0MaxDaughterDCA(1.3);
    AntiCascadeCuts->SetCutv0MaxDaughterDCA(1.3);

    CascadeCuts->SetCutv0CPA(0.96);
    AntiCascadeCuts->SetCutv0CPA(0.96);

    CascadeCuts->SetCutv0TransverseRadius(1.1, 200);
    AntiCascadeCuts->SetCutv0TransverseRadius(1.1, 200);

    XiNegCuts->SetEtaRange(-0.83, 0.83);
    XiPosCuts->SetEtaRange(-0.83, 0.83);
    XiBachCuts->SetEtaRange(-0.83, 0.83);
    AntiXiNegCuts->SetEtaRange(-0.83, 0.83);
    AntiXiPosCuts->SetEtaRange(-0.83, 0.83);
    AntiXiBachCuts->SetEtaRange(-0.83, 0.9);

    XiNegCuts->SetPID(AliPID::kPion, 999, 4.5);
    XiPosCuts->SetPID(AliPID::kProton, 999, 4.5);
    XiBachCuts->SetPID(AliPID::kPion, 999, 4.5);
    AntiXiNegCuts->SetPID(AliPID::kProton, 999, 4.5);
    AntiXiPosCuts->SetPID(AliPID::kPion, 999, 4.5);
    AntiXiBachCuts->SetPID(AliPID::kPion, 999, 4.5);

    CascadeCuts->SetPtRangeXi(0.8, 999.5);
    AntiCascadeCuts->SetPtRangeXi(0.8, 999.5);
  } else if (suffix == "41") {
    TrackCuts->SetPtRange(0.4, 4.05);
    AntiTrackCuts->SetPtRange(0.4, 4.05);

    TrackCuts->SetEtaRange(-0.77, 0.77);
    AntiTrackCuts->SetEtaRange(-0.77, 0.77);

    TrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);
    AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);

    config->SetDeltaEtaMax(0.012);
    config->SetDeltaPhiMax(0.012);

    v0Cuts->SetCutCPA(0.995);
    Antiv0Cuts->SetCutCPA(0.995);

    Posv0Daug->SetEtaRange(-0.77, 0.77);
    Negv0Daug->SetEtaRange(-0.77, 0.77);
    PosAntiv0Daug->SetEtaRange(-0.77, 0.77);
    NegAntiv0Daug->SetEtaRange(-0.77, 0.77);

    v0Cuts->SetCutDCADaugTov0Vtx(1.2);
    Antiv0Cuts->SetCutDCADaugTov0Vtx(1.2);

    CascadeCuts->SetCutXiDaughterDCA(1.3);
    AntiCascadeCuts->SetCutXiDaughterDCA(1.3);

    CascadeCuts->SetCutXiTransverseRadius(0.6, 200);
    AntiCascadeCuts->SetCutXiTransverseRadius(0.6, 200);

    CascadeCuts->SetCutv0CPA(0.96);
    AntiCascadeCuts->SetCutv0CPA(0.96);

    CascadeCuts->SetCutv0TransverseRadius(1.7, 200);
    AntiCascadeCuts->SetCutv0TransverseRadius(1.7, 200);

    XiNegCuts->SetEtaRange(-0.77, 0.77);
    XiPosCuts->SetEtaRange(-0.77, 0.77);
    XiBachCuts->SetEtaRange(-0.77, 0.77);
    AntiXiNegCuts->SetEtaRange(-0.77, 0.77);
    AntiXiPosCuts->SetEtaRange(-0.77, 0.77);
    AntiXiBachCuts->SetEtaRange(-0.77, 0.77);

    XiNegCuts->SetPID(AliPID::kPion, 999, 4.5);
    XiPosCuts->SetPID(AliPID::kProton, 999, 4.5);
    XiBachCuts->SetPID(AliPID::kPion, 999, 4.5);
    AntiXiNegCuts->SetPID(AliPID::kProton, 999, 4.5);
    AntiXiPosCuts->SetPID(AliPID::kPion, 999, 4.5);
    AntiXiBachCuts->SetPID(AliPID::kPion, 999, 4.5);

    CascadeCuts->SetPtRangeXi(0.8, 999.5);
    AntiCascadeCuts->SetPtRangeXi(0.8, 999.5);

  } else if (suffix == "42") {
    TrackCuts->SetPtRange(0.6, 4.05);
    AntiTrackCuts->SetPtRange(0.6, 4.05);

    TrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);
    AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);

    v0Cuts->SetCutCPA(0.995);
    Antiv0Cuts->SetCutCPA(0.995);

    Posv0Daug->SetPID(AliPID::kProton, 999.9, 4);
    Negv0Daug->SetPID(AliPID::kPion, 999.9, 4);
    PosAntiv0Daug->SetPID(AliPID::kPion, 999.9, 4);
    NegAntiv0Daug->SetPID(AliPID::kProton, 999.9, 4);

    Posv0Daug->SetNClsTPC(80);
    Negv0Daug->SetNClsTPC(80);
    PosAntiv0Daug->SetNClsTPC(80);
    NegAntiv0Daug->SetNClsTPC(80);

    v0Cuts->SetCutDCADaugTov0Vtx(1.2);
    Antiv0Cuts->SetCutDCADaugTov0Vtx(1.2);

    CascadeCuts->SetCutXiDaughterDCA(1.3);
    AntiCascadeCuts->SetCutXiDaughterDCA(1.3);

    CascadeCuts->SetCutXiCPA(0.99);
    AntiCascadeCuts->SetCutXiCPA(0.99);

    CascadeCuts->SetCutv0MaxDaughterDCA(1.3);
    AntiCascadeCuts->SetCutv0MaxDaughterDCA(1.3);

    CascadeCuts->SetCutv0CPA(0.96);
    AntiCascadeCuts->SetCutv0CPA(0.96);

    CascadeCuts->SetCutv0TransverseRadius(1.7, 200);
    AntiCascadeCuts->SetCutv0TransverseRadius(1.7, 200);

    CascadeCuts->SetCutv0MinDistToPrimVtx(0.06);
    AntiCascadeCuts->SetCutv0MinDistToPrimVtx(0.06);

    XiNegCuts->SetPID(AliPID::kPion, 999, 4.5);
    XiPosCuts->SetPID(AliPID::kProton, 999, 4.5);
    XiBachCuts->SetPID(AliPID::kPion, 999, 4.5);
    AntiXiNegCuts->SetPID(AliPID::kProton, 999, 4.5);
    AntiXiPosCuts->SetPID(AliPID::kPion, 999, 4.5);
    AntiXiBachCuts->SetPID(AliPID::kPion, 999, 4.5);
  } else if (suffix == "43") {
    TrackCuts->SetPtRange(0.6, 4.05);
    AntiTrackCuts->SetPtRange(0.6, 4.05);

    TrackCuts->SetEtaRange(-0.85, 0.85);
    AntiTrackCuts->SetEtaRange(-0.85, 0.85);

    TrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);
    AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);

    config->SetDeltaEtaMax(0.012);
    config->SetDeltaPhiMax(0.012);

    Posv0Daug->SetEtaRange(-0.83, 0.83);
    Negv0Daug->SetEtaRange(-0.83, 0.83);
    PosAntiv0Daug->SetEtaRange(-0.83, 0.83);
    NegAntiv0Daug->SetEtaRange(-0.83, 0.83);

    CascadeCuts->SetCutXiDaughterDCA(1.3);
    AntiCascadeCuts->SetCutXiDaughterDCA(1.3);

    CascadeCuts->SetCutXiCPA(0.985);
    AntiCascadeCuts->SetCutXiCPA(0.985);

    CascadeCuts->SetCutXiTransverseRadius(1.0, 200);
    AntiCascadeCuts->SetCutXiTransverseRadius(1.0, 200);

    CascadeCuts->SetCutv0MaxDaughterDCA(1.3);
    AntiCascadeCuts->SetCutv0MaxDaughterDCA(1.3);

    CascadeCuts->SetCutv0CPA(0.96);
    AntiCascadeCuts->SetCutv0CPA(0.96);

    CascadeCuts->SetCutv0MinDistToPrimVtx(0.06);
    AntiCascadeCuts->SetCutv0MinDistToPrimVtx(0.06);

    CascadeCuts->SetCutv0MinDaugDistToPrimVtx(0.06);
    AntiCascadeCuts->SetCutv0MinDaugDistToPrimVtx(0.06);
  } else if (suffix == "44") {
    TrackCuts->SetEtaRange(-0.85, 0.85);
    AntiTrackCuts->SetEtaRange(-0.85, 0.85);

    TrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);
    AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);

    config->SetDeltaEtaMax(0.012);
    config->SetDeltaPhiMax(0.012);

    v0Cuts->SetCutCPA(0.995);
    Antiv0Cuts->SetCutCPA(0.995);

    Posv0Daug->SetPID(AliPID::kProton, 999.9, 4);
    Negv0Daug->SetPID(AliPID::kPion, 999.9, 4);
    PosAntiv0Daug->SetPID(AliPID::kPion, 999.9, 4);
    NegAntiv0Daug->SetPID(AliPID::kProton, 999.9, 4);

    Posv0Daug->SetNClsTPC(80);
    Negv0Daug->SetNClsTPC(80);
    PosAntiv0Daug->SetNClsTPC(80);
    NegAntiv0Daug->SetNClsTPC(80);

    Posv0Daug->SetEtaRange(-0.83, 0.83);
    Negv0Daug->SetEtaRange(-0.83, 0.83);
    PosAntiv0Daug->SetEtaRange(-0.83, 0.83);
    NegAntiv0Daug->SetEtaRange(-0.83, 0.83);

    v0Cuts->SetCutDCADaugToPrimVtx(0.06);
    Antiv0Cuts->SetCutDCADaugToPrimVtx(0.06);

    CascadeCuts->SetCutXiDaughterDCA(1.9);
    AntiCascadeCuts->SetCutXiDaughterDCA(1.9);

    CascadeCuts->SetCutv0MaxDaughterDCA(1.4);
    AntiCascadeCuts->SetCutv0MaxDaughterDCA(1.4);

    CascadeCuts->SetCutv0CPA(0.96);
    AntiCascadeCuts->SetCutv0CPA(0.96);

    XiNegCuts->SetPID(AliPID::kPion, 999, 4.5);
    XiPosCuts->SetPID(AliPID::kProton, 999, 4.5);
    XiBachCuts->SetPID(AliPID::kPion, 999, 4.5);
    AntiXiNegCuts->SetPID(AliPID::kProton, 999, 4.5);
    AntiXiPosCuts->SetPID(AliPID::kPion, 999, 4.5);
    AntiXiBachCuts->SetPID(AliPID::kPion, 999, 4.5);

    CascadeCuts->SetPtRangeXi(0.8, 999.5);
    AntiCascadeCuts->SetPtRangeXi(0.8, 999.5);

  }

  task->SetTrackCuts(TrackCuts);
  task->SetAntiTrackCuts(AntiTrackCuts);
  task->Setv0Cuts(v0Cuts);
  task->SetAntiv0Cuts(Antiv0Cuts);
  task->SetCascadeCuts(CascadeCuts);
  task->SetAntiCascadeCuts(AntiCascadeCuts);
  task->SetCollectionConfig(config);
  mgr->AddTask(task);

  TString file = AliAnalysisManager::GetCommonFileName();
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();

  mgr->ConnectInput(task, 0, cinput);
  TString addon = "";
  if (CentEst == "kInt7") {
    addon += "MB";
  } else if (CentEst == "kHM") {
    addon += "HM";
  }
  AliAnalysisDataContainer *coutputQA;

  std::cout << "CONTAINTER NAME: " << addon.Data() << std::endl;
  TString QAName = Form("%sQA%s", addon.Data(), suffix.Data());
  coutputQA = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      QAName.Data(),
      TList::Class(), AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), QAName.Data()));
  mgr->ConnectOutput(task, 1, coutputQA);

  AliAnalysisDataContainer *coutputEvtCuts;
  TString EvtCutsName = Form("%sEvtCuts%s", addon.Data(), suffix.Data());
  coutputEvtCuts = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      EvtCutsName.Data(),
      TList::Class(), AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), EvtCutsName.Data()));
  mgr->ConnectOutput(task, 2, coutputEvtCuts);

  AliAnalysisDataContainer *couputTrkCuts;
  TString TrackCutsName = Form("%sTrackCuts%s", addon.Data(), suffix.Data());
  couputTrkCuts = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      TrackCutsName.Data(),
      TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), TrackCutsName.Data()));
  mgr->ConnectOutput(task, 3, couputTrkCuts);

  AliAnalysisDataContainer *coutputAntiTrkCuts;
  TString AntiTrackCutsName = Form("%sAntiTrackCuts%s", addon.Data(),
                                   suffix.Data());
  coutputAntiTrkCuts = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      AntiTrackCutsName.Data(),
      TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), AntiTrackCutsName.Data()));
  mgr->ConnectOutput(task, 4, coutputAntiTrkCuts);

  AliAnalysisDataContainer *coutputv0Cuts;
  TString v0CutsName = Form("%sv0Cuts%s", addon.Data(), suffix.Data());
  coutputv0Cuts = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      v0CutsName.Data(),
      TList::Class(), AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), v0CutsName.Data()));
  mgr->ConnectOutput(task, 5, coutputv0Cuts);

  AliAnalysisDataContainer *coutputAntiv0Cuts;
  TString Antiv0CutsName = Form("%sAntiv0Cuts%s", addon.Data(), suffix.Data());
  coutputAntiv0Cuts = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      Antiv0CutsName.Data(),
      TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), Antiv0CutsName.Data()));
  mgr->ConnectOutput(task, 6, coutputAntiv0Cuts);

  AliAnalysisDataContainer *coutputCascadeCuts;
  TString CascadeCutsName = Form("%sCascadeCuts%s", addon.Data(),
                                 suffix.Data());
  coutputCascadeCuts = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      CascadeCutsName.Data(),
      TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), CascadeCutsName.Data()));
  mgr->ConnectOutput(task, 7, coutputCascadeCuts);

  AliAnalysisDataContainer *coutputAntiCascadeCuts;
  TString AntiCascadeCutsName = Form("%sAntiCascadeCuts%s", addon.Data(),
                                     suffix.Data());
  coutputAntiCascadeCuts = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      AntiCascadeCutsName.Data(),
      TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), AntiCascadeCutsName.Data()));
  mgr->ConnectOutput(task, 8, coutputAntiCascadeCuts);

  AliAnalysisDataContainer *coutputResults;
  TString ResultsName = Form("%sResults%s", addon.Data(), suffix.Data());
  coutputResults = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      ResultsName.Data(),
      TList::Class(), AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), ResultsName.Data()));
  mgr->ConnectOutput(task, 9, coutputResults);

  AliAnalysisDataContainer *coutputResultQA;
  TString ResultQAName = Form("%sResultQA%s", addon.Data(), suffix.Data());
  coutputResultQA = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      ResultQAName.Data(),
      TList::Class(), AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), ResultQAName.Data()));
  mgr->ConnectOutput(task, 10, coutputResultQA);

  AliAnalysisDataContainer *coutputResultsSample;
  TString ResultsSampleName = Form("%sResultsSample%s", addon.Data(),
                                   suffix.Data());
  coutputResultsSample = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      ResultsSampleName.Data(),
      TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), ResultsSampleName.Data()));
  mgr->ConnectOutput(task, 11, coutputResultsSample);

  AliAnalysisDataContainer *coutputResultQASample;
  TString ResultQASampleName = Form("%sResultQASample%s", addon.Data(),
                                    suffix.Data());
  coutputResultQASample = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      ResultQASampleName.Data(),
      TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), ResultQASampleName.Data()));
  mgr->ConnectOutput(task, 12, coutputResultQASample);

  if (isMC) {
    AliAnalysisDataContainer *coutputTrkCutsMC;
    TString TrkCutsMCName = Form("%sTrkCutsMC%s", addon.Data(), suffix.Data());
    coutputTrkCutsMC = mgr->CreateContainer(
        //@suppress("Invalid arguments") it works ffs
        TrkCutsMCName.Data(),
        TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), TrkCutsMCName.Data()));
    mgr->ConnectOutput(task, 13, coutputTrkCutsMC);

    AliAnalysisDataContainer *coutputAntiTrkCutsMC;
    TString AntiTrkCutsMCName = Form("%sAntiTrkCutsMC%s", addon.Data(),
                                     suffix.Data());
    coutputAntiTrkCutsMC = mgr->CreateContainer(
        //@suppress("Invalid arguments") it works ffs
        AntiTrkCutsMCName.Data(),
        TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), AntiTrkCutsMCName.Data()));
    mgr->ConnectOutput(task, 14, coutputAntiTrkCutsMC);

    AliAnalysisDataContainer *coutputv0CutsMC;
    TString v0CutsMCName = Form("%sv0CutsMC%s", addon.Data(), suffix.Data());
    coutputv0CutsMC = mgr->CreateContainer(
        //@suppress("Invalid arguments") it works ffs
        v0CutsMCName.Data(),
        TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), v0CutsMCName.Data()));
    mgr->ConnectOutput(task, 15, coutputv0CutsMC);

    AliAnalysisDataContainer *coutputAntiv0CutsMC;
    TString Antiv0CutsMCName = Form("%sAntiv0CutsMC%s", addon.Data(),
                                    suffix.Data());
    coutputAntiv0CutsMC = mgr->CreateContainer(
        //@suppress("Invalid arguments") it works ffs
        Antiv0CutsMCName.Data(),
        TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), Antiv0CutsMCName.Data()));
    mgr->ConnectOutput(task, 16, coutputAntiv0CutsMC);

    AliAnalysisDataContainer *coutputXiCutsMC;
    TString XiCutsMCName = Form("%sXiCutsMC%s", addon.Data(), suffix.Data());
    coutputXiCutsMC = mgr->CreateContainer(
        //@suppress("Invalid arguments") it works ffs
        XiCutsMCName.Data(),
        TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), XiCutsMCName.Data()));
    mgr->ConnectOutput(task, 17, coutputXiCutsMC);

    AliAnalysisDataContainer *coutputAntiXiCutsMC;
    TString AntiXiCutsMCName = Form("%sAntiXiCutsMC%s", addon.Data(),
                                    suffix.Data());
    coutputAntiXiCutsMC = mgr->CreateContainer(
        //@suppress("Invalid arguments") it works ffs
        AntiXiCutsMCName.Data(),
        TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), AntiXiCutsMCName.Data()));
    mgr->ConnectOutput(task, 18, coutputAntiXiCutsMC);
  }
  return task;
}

//
//if (suffix == "1") {
//    TrackCuts->SetPtRange(0.4, 4.05);
//    AntiTrackCuts->SetPtRange(0.4, 4.05);
//  } else if (suffix == "2") {
//    TrackCuts->SetPtRange(0.6, 4.05);
//    AntiTrackCuts->SetPtRange(0.6, 4.05);
//  } else if (suffix == "3") {
//    TrackCuts->SetEtaRange(-0.7, 0.7);
//    AntiTrackCuts->SetEtaRange(-0.7, 0.7);
//  } else if (suffix == "4") {
//    TrackCuts->SetEtaRange(-0.9, 0.9);
//    AntiTrackCuts->SetEtaRange(-0.9, 0.9);
//  } else if (suffix == "5") {
//    TrackCuts->SetPID(AliPID::kProton, 0.75, 2);
//    AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 2);
//  } else if (suffix == "6") {
//    TrackCuts->SetPID(AliPID::kProton, 0.75, 5);
//    AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 5);
//  } else if (suffix == "7") {
//    TrackCuts->SetFilterBit(96);
//    AntiTrackCuts->SetFilterBit(96);
//  } else if (suffix == "8") {
//    TrackCuts->SetNClsTPC(70);
//    AntiTrackCuts->SetNClsTPC(70);
//  } else if (suffix == "9") {
//    TrackCuts->SetNClsTPC(90);
//    AntiTrackCuts->SetNClsTPC(90);
//  }
// if (suffix == "10") {
//   v0Cuts->SetPtRange(0.24, 999.9);
//   Antiv0Cuts->SetPtRange(0.24, 999.9);
// } else if (suffix == "11") {
//   v0Cuts->SetPtRange(0.36, 999.9);
//   Antiv0Cuts->SetPtRange(0.36, 999.9);
// } else if (suffix == "12") {
//   v0Cuts->SetCutCPA(0.998);
//   Antiv0Cuts->SetCutCPA(0.998);
// } else if (suffix == "13") {
//   Posv0Daug->SetPID(AliPID::kProton, 999.9, 4);
//   Negv0Daug->SetPID(AliPID::kPion, 999.9, 4);
//   PosAntiv0Daug->SetPID(AliPID::kPion, 999.9, 4);
//   NegAntiv0Daug->SetPID(AliPID::kProton, 999.9, 4);
// } else if (suffix == "14") {
//   Posv0Daug->SetNClsTPC(80);
//   Negv0Daug->SetNClsTPC(80);
//   PosAntiv0Daug->SetNClsTPC(80);
//   NegAntiv0Daug->SetNClsTPC(80);
// } else if (suffix == "15") {
//   Posv0Daug->SetEtaRange(-0.7, 0.7);
//   Negv0Daug->SetEtaRange(-0.7, 0.7);
//   PosAntiv0Daug->SetEtaRange(-0.7, 0.7);
//   NegAntiv0Daug->SetEtaRange(-0.7, 0.7);
// } else if (suffix == "16") {
//   Posv0Daug->SetEtaRange(-0.9, 0.90.83)/   Negv0Daug->SetEtaRange(-0.9, 0.9);
//   PosAntiv0Daug->SetEtaRange(-0.9, 0.9);
//   NegAntiv0Daug->SetEtaRange(-0.9, 0.9);
// } else if (suffix == "17") {
//   v0Cuts->SetCutDCADaugTov0Vtx(1.2);
//   Antiv0Cuts->SetCutDCADaugTov0Vtx(1.2);
// } else if (suffix == "18") {
//   v0Cuts->SetCutDCADaugToPrimVtx(0.06);
//   Antiv0Cuts->SetCutDCADaugToPrimVtx(0.06);
// }
//if (suffix == "19") {
//  CascadeCuts->SetCutXiDaughterDCA(1.3);
//  AntiCascadeCuts->SetCutXiDaughterDCA(1.3);
//} else if (suffix == "20") {
//  CascadeCuts->SetCutXiDaughterDCA(1.9);
//  AntiCascadeCuts->SetCutXiDaughterDCA(1.9);
//} else if (suffix == "21") {
//  CascadeCuts->SetCutXiMinDistBachToPrimVtx(0.06);
//  AntiCascadeCuts->SetCutXiMinDistBachToPrimVtx(0.06);
//} else if (suffix == "22") {
//  CascadeCuts->SetCutXiMinDistBachToPrimVtx(0.04);
//  AntiCascadeCuts->SetCutXiMinDistBachToPrimVtx(0.04);
//} else if (suffix == "23") {
//  CascadeCuts->SetCutXiCPA(0.99);
//  AntiCascadeCuts->SetCutXiCPA(0.99);
//} else if (suffix == "24") {
//  CascadeCuts->SetCutXiCPA(0.985);
//  AntiCascadeCuts->SetCutXiCPA(0.985);
//} else if (suffix == "25") {
//  CascadeCuts->SetCutXiTransverseRadius(0.6, 200);
//  AntiCascadeCuts->SetCutXiTransverseRadius(0.6, 200);
//} else if (suffix == "26") {
//  CascadeCuts->SetCutXiTransverseRadius(1.0, 200);
//  AntiCascadeCuts->SetCutXiTransverseRadius(1.0, 200);
//} else if (suffix == "27") {
//  CascadeCuts->SetCutv0MaxDaughterDCA(1.3);
//  AntiCascadeCuts->SetCutv0MaxDaughterDCA(1.3);
//} else if (suffix == "28") {
//  CascadeCuts->SetCutv0MaxDaughterDCA(1.4);
//  AntiCascadeCuts->SetCutv0MaxDaughterDCA(1.4);
//} else if (suffix == "29") {
//  CascadeCuts->SetCutv0CPA(0.96);
//  AntiCascadeCuts->SetCutv0CPA(0.96);
//} else if (suffix == "30") {
//  CascadeCuts->SetCutv0CPA(0.99);
//  AntiCascadeCuts->SetCutv0CPA(0.99);
//} else if (suffix == "31") {
//  CascadeCuts->SetCutv0TransverseRadius(1.1, 200);
//  AntiCascadeCuts->SetCutv0TransverseRadius(1.1, 200);
//} else if (suffix == "32") {
//  CascadeCuts->SetCutv0TransverseRadius(1.7, 200);
//  AntiCascadeCuts->SetCutv0TransverseRadius(1.7, 200);
//} else if (suffix == "33") {
//  CascadeCuts->SetCutv0MinDistToPrimVtx(0.08);
//  AntiCascadeCuts->SetCutv0MinDistToPrimVtx(0.08);
//} else if (suffix == "34") {
//  CascadeCuts->SetCutv0MinDistToPrimVtx(0.06);
//  AntiCascadeCuts->SetCutv0MinDistToPrimVtx(0.06);
//} else if (suffix == "35") {
//  CascadeCuts->SetCutv0MinDaugDistToPrimVtx(0.05);
//  AntiCascadeCuts->SetCutv0MinDaugDistToPrimVtx(0.05);
//} else if (suffix == "36") {
//  CascadeCuts->SetCutv0MinDaugDistToPrimVtx(0.06);
//  AntiCascadeCuts->SetCutv0MinDaugDistToPrimVtx(0.06);
//} else if (suffix == "37") {
//  XiNegCuts->SetEtaRange(-0.7, 0.7);
//  XiPosCuts->SetEtaRange(-0.7, 0.7);
//  XiBachCuts->SetEtaRange(-0.7, 0.7);
//  AntiXiNegCuts->SetEtaRange(-0.7, 0.7);
//  AntiXiPosCuts->SetEtaRange(-0.7, 0.7);
//  AntiXiBachCuts->SetEtaRange(-0.7, 0.7);
//} else if (suffix == "38") {
//  XiNegCuts->SetEtaRange(-0.9, 0.9);
//  XiPosCuts->SetEtaRange(-0.9, 0.9);
//  XiBachCuts->SetEtaRange(-0.9, 0.9);
//  AntiXiNegCuts->SetEtaRange(-0.9, 0.9);
//  AntiXiPosCuts->SetEtaRange(-0.9, 0.9);
//  AntiXiBachCuts->SetEtaRange(-0.9, 0.9);
//} else if (suffix == "39") {
//  XiNegCuts->SetPID(AliPID::kPion, 999, 3);
//  XiPosCuts->SetPID(AliPID::kProton, 999, 3);
//  XiBachCuts->SetPID(AliPID::kPion, 999, 3);
//  AntiXiNegCuts->SetPID(AliPID::kProton, 999, 3);
//  AntiXiPosCuts->SetPID(AliPID::kPion, 999, 3);
//  AntiXiBachCuts->SetPID(AliPID::kPion, 999, 3);
//} else if (suffix == "40") {
//  XiNegCuts->SetPID(AliPID::kPion, 999, 4.5);
//  XiPosCuts->SetPID(AliPID::kProton, 999, 4.5);
//  XiBachCuts->SetPID(AliPID::kPion, 999, 4.5);
//  AntiXiNegCuts->SetPID(AliPID::kProton, 999, 4.5);
//  AntiXiPosCuts->SetPID(AliPID::kPion, 999, 4.5);
//  AntiXiBachCuts->SetPID(AliPID::kPion, 999, 4.5);
//} else if (suffix == "41") {
//  CascadeCuts->SetPtRangeXi(0.8, 999.5);
//  AntiCascadeCuts->SetPtRangeXi(0.8, 999.5);
//} else if (suffix == "42") {
//  CascadeCuts->SetPtRangeXi(1.2, 999.5);
//  AntiCascadeCuts->SetPtRangeXi(1.2, 999.5);
//}
