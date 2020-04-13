#include <vector>
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskNanoXioton.h"
#include "AliAnalysisTaskAODXioton.h"
#include "AliFemtoDreamEventCuts.h"
#include "AliFemtoDreamTrackCuts.h"
#include "AliFemtoDreamCascadeCuts.h"
#include "AliFemtoDreamCollConfig.h"

AliAnalysisTaskSE *AddTaskFemtoXoton(bool fullBlastQA = false,
                                     bool isMC = false, bool isNano = true,
                                     int iDepth = 10,
                                     bool Systematic = false,
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

  //========= Init subtasks and start analyis ============================
  // Event Cuts
  AliFemtoDreamEventCuts *evtCuts = AliFemtoDreamEventCuts::StandardCutsRun2();
  evtCuts->CleanUpMult(false, false, false, true);
  evtCuts->SetMultVsCentPlots(true);
  // Track Cuts
  AliFemtoDreamTrackCuts *TrackCuts = AliFemtoDreamTrackCuts::PrimProtonCuts(
      isMC, true, false, false);
  TrackCuts->SetFilterBit(128);
  TrackCuts->SetCutCharge(1);

  AliFemtoDreamTrackCuts *AntiTrackCuts =
      AliFemtoDreamTrackCuts::PrimProtonCuts(isMC, true, false, false);
  AntiTrackCuts->SetFilterBit(128);
  AntiTrackCuts->SetCutCharge(-1);

  //Cascade Cuts
  AliFemtoDreamCascadeCuts* CascadeCuts = AliFemtoDreamCascadeCuts::XiCuts(
      isMC, false);
  CascadeCuts->SetXiCharge(-1);
  AliFemtoDreamTrackCuts *XiNegCuts = AliFemtoDreamTrackCuts::Xiv0PionCuts(
      isMC, true, false);
  XiNegCuts->SetCheckTPCRefit(false);  //for nanos this is already done while prefiltering
  AliFemtoDreamTrackCuts *XiPosCuts = AliFemtoDreamTrackCuts::Xiv0ProtonCuts(
      isMC, true, false);
  XiPosCuts->SetCheckTPCRefit(false);  //for nanos this is already done while prefiltering
  AliFemtoDreamTrackCuts *XiBachCuts = AliFemtoDreamTrackCuts::XiBachPionCuts(
      isMC, true, false);
  XiBachCuts->SetCheckTPCRefit(false);  //for nanos this is already done while prefiltering

  CascadeCuts->Setv0Negcuts(XiNegCuts);
  CascadeCuts->Setv0PosCuts(XiPosCuts);
  CascadeCuts->SetBachCuts(XiBachCuts);
  CascadeCuts->SetPDGCodeCasc(3312);
  CascadeCuts->SetPDGCodev0(3122);
  CascadeCuts->SetPDGCodePosDaug(2212);
  CascadeCuts->SetPDGCodeNegDaug(-211);
  CascadeCuts->SetPDGCodeBach(-211);

  AliFemtoDreamCascadeCuts* AntiCascadeCuts = AliFemtoDreamCascadeCuts::XiCuts(
      isMC, false);
  AntiCascadeCuts->SetXiCharge(1);
  AliFemtoDreamTrackCuts *AntiXiNegCuts =
      AliFemtoDreamTrackCuts::Xiv0ProtonCuts(isMC, true, false);
  AntiXiNegCuts->SetCutCharge(-1);
  AntiXiNegCuts->SetCheckTPCRefit(false);  //for nanos this is already done while prefiltering
  AliFemtoDreamTrackCuts *AntiXiPosCuts = AliFemtoDreamTrackCuts::Xiv0PionCuts(
      isMC, true, false);
  AntiXiPosCuts->SetCutCharge(1);
  AntiXiPosCuts->SetCheckTPCRefit(false);  //for nanos this is already done while prefiltering
  AliFemtoDreamTrackCuts *AntiXiBachCuts =
      AliFemtoDreamTrackCuts::XiBachPionCuts(isMC, true, false);
  AntiXiBachCuts->SetCutCharge(1);
  AntiXiBachCuts->SetCheckTPCRefit(false);  //for nanos this is already done while prefiltering

  AntiCascadeCuts->Setv0Negcuts(AntiXiNegCuts);
  AntiCascadeCuts->Setv0PosCuts(AntiXiPosCuts);
  AntiCascadeCuts->SetBachCuts(AntiXiBachCuts);
  AntiCascadeCuts->SetPDGCodeCasc(-3312);
  AntiCascadeCuts->SetPDGCodev0(-3122);
  AntiCascadeCuts->SetPDGCodePosDaug(211);
  AntiCascadeCuts->SetPDGCodeNegDaug(-2212);
  AntiCascadeCuts->SetPDGCodeBach(211);

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

  if (Systematic) {
    evtCuts->SetMinimalBooking(true);
    TrackCuts->SetMinimalBooking(true);
    AntiTrackCuts->SetMinimalBooking(true);
    CascadeCuts->SetMinimalBooking(true);
    AntiCascadeCuts->SetMinimalBooking(true);
  }

  AliFemtoDreamCollConfig *config = new AliFemtoDreamCollConfig("Femto",
                                                                "Femto");

  // Femto Collection
  std::vector<int> PDGParticles;
  PDGParticles.push_back(2212);
  PDGParticles.push_back(2212);
  PDGParticles.push_back(3312);
  PDGParticles.push_back(3312);

  std::vector<int> NBins;
  std::vector<float> kMin;
  std::vector<float> kMax;
  std::vector<int> pairQA;
  std::vector<bool> closeRejection;
  //pairs:
  //pp                0
  //p bar p           1
  //p Xi              2
  //p bar Xi          3
  //bar p bar p       4
  //bar p Xi          5
  //bar p bar Xi      6
  //Xi Xi             7
  //Xi bar Xi         8
  //bar Xi bar Xi     9
  const int nPairs = 10;
  for (int i = 0; i < nPairs; ++i) {
    pairQA.push_back(0);
    closeRejection.push_back(false);
    NBins.push_back(1000);
    kMin.push_back(0.);
    kMax.push_back(1.);
  }
  pairQA[0] = 11;
  pairQA[4] = 11;
  pairQA[2] = 13;
  pairQA[6] = 13;

  closeRejection[0] = true;  // pp
  closeRejection[4] = true;  // barp barp

  config->SetPDGCodes(PDGParticles);
  config->SetNBinsHist(NBins);
  config->SetMinKRel(kMin);
  config->SetMaxKRel(kMax);
  config->SetClosePairRejection(closeRejection);
  config->SetDeltaEtaMax(0.02);
  config->SetDeltaPhiMax(0.020);
  config->SetExtendedQAPairs(pairQA);

  config->SetMixingDepth(iDepth);
  config->SetUseEventMixing(true);
  config->SetMultiplicityEstimator(AliFemtoDreamEvent::kRef08);
  if (isMC) {
    config->SetMomentumResolution(true);
  }
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

  config->SetZBins(ZVtxBins);

  config->SetMultBinning(true);
  config->SetmTBinning(true);

  config->SetdPhidEtaPlotsSmallK(true);
  config->SetdPhidEtaPlots(true);
  config->SetPhiEtaBinnign(true);

  if (fullBlastQA) {
    config->SetkTBinning(true);
    config->SetPtQA(true);
    config->SetMassQA(true);
  }

  if (!fullBlastQA || Systematic) {
    config->SetMinimalBookingME(true);
    config->SetMinimalBookingSample(true);
  }

  if (Systematic) {
    if (suffix == "1") {
      TrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);
      AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);

      TrackCuts->SetNClsTPC(90);
      AntiTrackCuts->SetNClsTPC(90);

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

      config->SetDeltaEtaMax(0.022);
      config->SetDeltaPhiMax(0.022);

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

      config->SetDeltaEtaMax(0.022);
      config->SetDeltaPhiMax(0.022);

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

      config->SetDeltaEtaMax(0.022);
      config->SetDeltaPhiMax(0.022);

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

      config->SetDeltaEtaMax(0.022);
      config->SetDeltaPhiMax(0.022);

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

      config->SetDeltaEtaMax(0.022);
      config->SetDeltaPhiMax(0.022);

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

      config->SetDeltaEtaMax(0.022);
      config->SetDeltaPhiMax(0.022);

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

      config->SetDeltaEtaMax(0.022);
      config->SetDeltaPhiMax(0.022);

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

      config->SetDeltaEtaMax(0.022);
      config->SetDeltaPhiMax(0.022);

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

      config->SetDeltaEtaMax(0.022);
      config->SetDeltaPhiMax(0.022);

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

      config->SetDeltaEtaMax(0.022);
      config->SetDeltaPhiMax(0.022);

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

      config->SetDeltaEtaMax(0.022);
      config->SetDeltaPhiMax(0.022);

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

      config->SetDeltaEtaMax(0.022);
      config->SetDeltaPhiMax(0.022);

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

      config->SetDeltaEtaMax(0.022);
      config->SetDeltaPhiMax(0.022);

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

      config->SetDeltaEtaMax(0.022);
      config->SetDeltaPhiMax(0.022);

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

      config->SetDeltaEtaMax(0.022);
      config->SetDeltaPhiMax(0.022);

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

      config->SetDeltaEtaMax(0.022);
      config->SetDeltaPhiMax(0.022);

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

      config->SetDeltaEtaMax(0.022);
      config->SetDeltaPhiMax(0.022);

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

      config->SetDeltaEtaMax(0.022);
      config->SetDeltaPhiMax(0.022);

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

      config->SetDeltaEtaMax(0.022);
      config->SetDeltaPhiMax(0.022);

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

      config->SetDeltaEtaMax(0.022);
      config->SetDeltaPhiMax(0.022);

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

      config->SetDeltaEtaMax(0.022);
      config->SetDeltaPhiMax(0.022);

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

      config->SetDeltaEtaMax(0.022);
      config->SetDeltaPhiMax(0.022);

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

      config->SetDeltaEtaMax(0.022);
      config->SetDeltaPhiMax(0.022);

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

      config->SetDeltaEtaMax(0.022);
      config->SetDeltaPhiMax(0.022);

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

      config->SetDeltaEtaMax(0.022);
      config->SetDeltaPhiMax(0.022);

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

      config->SetDeltaEtaMax(0.022);
      config->SetDeltaPhiMax(0.022);

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

      config->SetDeltaEtaMax(0.022);
      config->SetDeltaPhiMax(0.022);

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

      config->SetDeltaEtaMax(0.022);
      config->SetDeltaPhiMax(0.022);

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
  }

  TString addon = "PXi";

  TString file = AliAnalysisManager::GetCommonFileName();

  TString EvtCutsName = Form("%sEvtCuts%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *coutputEvtCuts = mgr->CreateContainer(
      EvtCutsName.Data(), TList::Class(), AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), EvtCutsName.Data()));

  TString TrackCutsName = Form("%sTrackCuts%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *couputTrkCuts = mgr->CreateContainer(
      TrackCutsName.Data(), TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), TrackCutsName.Data()));

  TString AntiTrackCutsName = Form("%sAntiTrackCuts%s", addon.Data(),
                                   suffix.Data());
  AliAnalysisDataContainer *coutputAntiTrkCuts = mgr->CreateContainer(
      AntiTrackCutsName.Data(), TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), AntiTrackCutsName.Data()));

  AliAnalysisDataContainer *coutputCascadeCuts;
  TString CascadeCutsName = Form("%sCascadeCuts%s", addon.Data(),
                                 suffix.Data());
  coutputCascadeCuts = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      CascadeCutsName.Data(),
      TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), CascadeCutsName.Data()));

  AliAnalysisDataContainer *coutputAntiCascadeCuts;
  TString AntiCascadeCutsName = Form("%sAntiCascadeCuts%s", addon.Data(),
                                     suffix.Data());
  coutputAntiCascadeCuts = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      AntiCascadeCutsName.Data(),
      TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), AntiCascadeCutsName.Data()));

  AliAnalysisDataContainer *coutputResults;
  TString ResultsName = Form("%sResults%s", addon.Data(), suffix.Data());
  coutputResults = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      ResultsName.Data(),
      TList::Class(), AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), ResultsName.Data()));

  AliAnalysisDataContainer *coutputResultsQA;
  TString ResultsQAName = Form("%sResultsQA%s", addon.Data(), suffix.Data());
  coutputResultsQA = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      ResultsQAName.Data(),
      TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), ResultsQAName.Data()));

  AliAnalysisDataContainer *coutputDumpster;
  TString DumpsterName = Form("%sDumpster%s", addon.Data(), suffix.Data());
  coutputDumpster = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      DumpsterName.Data(),
      TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), DumpsterName.Data()));

  AliAnalysisDataContainer *coutputTrkCutsMC;
  AliAnalysisDataContainer *coutputAntiTrkCutsMC;
  AliAnalysisDataContainer *coutputCascCutsMC;
  AliAnalysisDataContainer *coutputAntiCascCutsMC;
  if (isMC) {
    TString TrkCutsMCName = Form("%sTrkCutsMC%s", addon.Data(), suffix.Data());
    coutputTrkCutsMC = mgr->CreateContainer(
        //@suppress("Invalid arguments") it works ffs
        TrkCutsMCName.Data(),
        TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), TrkCutsMCName.Data()));

    TString AntiTrkCutsMCName = Form("%sAntiTrkCutsMC%s", addon.Data(),
                                     suffix.Data());
    coutputAntiTrkCutsMC = mgr->CreateContainer(
        //@suppress("Invalid arguments") it works ffs
        AntiTrkCutsMCName.Data(),
        TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), AntiTrkCutsMCName.Data()));

    TString CascCutsMCName = Form("%sCascCutsMC%s", addon.Data(),
                                  suffix.Data());
    coutputCascCutsMC = mgr->CreateContainer(
        //@suppress("Invalid arguments") it works ffs
        CascCutsMCName.Data(),
        TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), CascCutsMCName.Data()));

    TString AntiCascCutsMCName = Form("%sAntiCascCutsMC%s", addon.Data(),
                                      suffix.Data());
    coutputAntiCascCutsMC = mgr->CreateContainer(
        //@suppress("Invalid arguments") it works ffs
        AntiCascCutsMCName.Data(),
        TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), AntiCascCutsMCName.Data()));
  }

  AliAnalysisTaskNanoXioton* taskNano;
  AliAnalysisTaskAODXioton* taskAOD;

  if (isNano) {
    taskNano = new AliAnalysisTaskNanoXioton("femtoNanoXoton", isMC);
    if (!fullBlastQA) {
      taskNano->SetRunTaskLightWeight(true);
    }
    taskNano->SelectCollisionCandidates(AliVEvent::kHighMultV0);
    taskNano->SetEventCuts(evtCuts);
    taskNano->SetProtonCuts(TrackCuts);
    taskNano->SetAntiProtonCuts(AntiTrackCuts);
    taskNano->SetXiCuts(CascadeCuts);
    taskNano->SetAntiXiCuts(AntiCascadeCuts);
    taskNano->SetCorrelationConfig(config);
    mgr->AddTask(taskNano);

    mgr->ConnectInput(taskNano, 0, cinput);
    mgr->ConnectOutput(taskNano, 1, coutputEvtCuts);



    mgr->ConnectOutput(taskNano, 2, couputTrkCuts);
    mgr->ConnectOutput(taskNano, 3, coutputAntiTrkCuts);
    mgr->ConnectOutput(taskNano, 4, coutputCascadeCuts);
    mgr->ConnectOutput(taskNano, 5, coutputAntiCascadeCuts);
    mgr->ConnectOutput(taskNano, 6, coutputResults);
    mgr->ConnectOutput(taskNano, 7, coutputResultsQA);
    mgr->ConnectOutput(taskNano, 8, coutputDumpster);
    if (isMC) {
      mgr->ConnectOutput(taskNano, 9, coutputTrkCutsMC);
      mgr->ConnectOutput(taskNano, 10, coutputAntiTrkCutsMC);
      mgr->ConnectOutput(taskNano, 11, coutputCascCutsMC);
      mgr->ConnectOutput(taskNano, 12, coutputAntiCascCutsMC);
    }
  } else {
    taskAOD = new AliAnalysisTaskAODXioton("femtoAODXoton", isMC);
    if (!fullBlastQA) {
      taskAOD->SetRunTaskLightWeight(true);
    }
    taskAOD->SelectCollisionCandidates(AliVEvent::kHighMultV0);
    taskAOD->SetEventCuts(evtCuts);
    taskAOD->SetProtonCuts(TrackCuts);
    taskAOD->SetAntiProtonCuts(AntiTrackCuts);
    taskAOD->SetXiCuts(CascadeCuts);
    taskAOD->SetAntiXiCuts(AntiCascadeCuts);
    taskAOD->SetCorrelationConfig(config);
    mgr->AddTask(taskAOD);
    mgr->ConnectInput(taskAOD, 0, cinput);
    mgr->ConnectOutput(taskAOD, 1, coutputEvtCuts);
    mgr->ConnectOutput(taskAOD, 2, couputTrkCuts);
    mgr->ConnectOutput(taskAOD, 3, coutputAntiTrkCuts);
    mgr->ConnectOutput(taskAOD, 4, coutputCascadeCuts);
    mgr->ConnectOutput(taskAOD, 5, coutputAntiCascadeCuts);
    mgr->ConnectOutput(taskAOD, 6, coutputResults);
    mgr->ConnectOutput(taskAOD, 7, coutputResultsQA);
    if (isMC) {
      mgr->ConnectOutput(taskAOD, 8, coutputTrkCutsMC);
      mgr->ConnectOutput(taskAOD, 9, coutputAntiTrkCutsMC);
      mgr->ConnectOutput(taskAOD, 10, coutputCascCutsMC);
      mgr->ConnectOutput(taskAOD, 11, coutputAntiCascCutsMC);
    }
  }
  if (isNano) {
    return taskNano;
  } else {
    return taskAOD;
  }

}
