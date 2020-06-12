#include <vector>
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliFemtoDreamEventCuts.h"
#include "AliFemtoDreamTrackCuts.h"
#include "AliFemtoDreamCollConfig.h"
#include "AliAnalysisTaskNanoPt.h"
AliAnalysisTaskSE* AddTaskNanoPt( bool isMC = true, bool fIsMCTruth = true,
                                  TString trigger = "kINT7", //2
                                  bool DCAPlots = false, //3
                                  bool CombSigma = false, //4
                                  bool ContributionSplitting = false, //5,
                                  bool DumpPdApAd = true, //6
                                  bool fullBlastQA = true, bool RefMult08 = true, bool Systematic = false,
                                  bool SystematicpTCutVariation = true, const char *cutVariation = "0") {

  TString suffix = TString::Format("%s", cutVariation);

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  if (!mgr) {
    printf("No analysis manager to connect to!\n");
    return nullptr;
  }

  if (!mgr->GetInputEventHandler()) {
    printf("This task requires an input event handler!\n");
    return nullptr;
  }
  //========= Init subtasks and start analyis ============================
  // Event Cuts
  AliFemtoDreamEventCuts *evtCuts = AliFemtoDreamEventCuts::StandardCutsRun2();
  evtCuts->CleanUpMult(false, false, false, true);

  // =====================================================================
  //Proton track Cuts----------------------------------------------------------------------------
  AliFemtoDreamTrackCuts *TrackCuts = AliFemtoDreamTrackCuts::PrimProtonCuts(
                                        isMC, true, CombSigma, ContributionSplitting);
  TrackCuts->SetMinimalBooking(SystematicpTCutVariation);
  TrackCuts->SetCutCharge(1);
  //Antiproton track Cuts-------------------------------------------------------------------------
  AliFemtoDreamTrackCuts *AntiTrackCuts = AliFemtoDreamTrackCuts::PrimProtonCuts(
      isMC, true, CombSigma, ContributionSplitting);
  AntiTrackCuts->SetMinimalBooking(SystematicpTCutVariation);
  AntiTrackCuts->SetCutCharge(-1);
  //deuteron track cuts----------------------------------------------------------------------------
  AliFemtoDreamTrackCuts *TrackCutsDeuteron =
    AliFemtoDreamTrackCuts::PrimDeuteronCuts(isMC, true, CombSigma,
        ContributionSplitting);
  TrackCutsDeuteron->SetMinimalBooking(SystematicpTCutVariation);
  TrackCutsDeuteron->SetCutCharge(1);
  //Antideuteron track cuts----------------------------------------------------------------------------
  AliFemtoDreamTrackCuts *AntiTrackCutsDeuteron =
    AliFemtoDreamTrackCuts::PrimDeuteronCuts(isMC, true, CombSigma,
        ContributionSplitting);
  AntiTrackCutsDeuteron->SetMinimalBooking(SystematicpTCutVariation);
  AntiTrackCutsDeuteron->SetCutCharge(-1);
/////////////////////For no NSigmaTOF information///
// =====================================================================
  //deuteron track cuts----------------------------------------------------------------------------
  AliFemtoDreamTrackCuts *TrackCutsDeuteronNoTOF =
    AliFemtoDreamTrackCuts::PrimDeuteronCuts(isMC, true, CombSigma,
        ContributionSplitting);
  TrackCutsDeuteronNoTOF->SetMinimalBooking(SystematicpTCutVariation);
  TrackCutsDeuteronNoTOF->SetCutCharge(1);
  TrackCutsDeuteronNoTOF->SetPID(AliPID::kDeuteron, 999.);
  //Antideuteron track cuts----------------------------------------------------------------------------
  AliFemtoDreamTrackCuts *AntiTrackCutsDeuteronNoTOF =
    AliFemtoDreamTrackCuts::PrimDeuteronCuts(isMC, true, CombSigma,
        ContributionSplitting);
  AntiTrackCutsDeuteronNoTOF->SetMinimalBooking(SystematicpTCutVariation);
  AntiTrackCutsDeuteronNoTOF->SetCutCharge(-1);
  AntiTrackCutsDeuteronNoTOF->SetPID(AliPID::kDeuteron, 999.);
//====================================================================================================================================
  std::vector<int> PDGParticles;
  PDGParticles.push_back(2212);
  PDGParticles.push_back(2212);
  PDGParticles.push_back(1000010020);
  PDGParticles.push_back(1000010020);

  std::vector<int> NBins;
  std::vector<float> kMin;
  std::vector<float> kMax;
  std::vector<bool> closeRejection;
  std::vector<float> mTBins = { 1.14, 1.26, 999. };
  std::vector<int> pairQA;
  //pairs:
  // pp             0
  // p bar p        1
  // p d            2
  // p bar d        3
  // bar p bar p    4
  // bar p d        5
  // bar p bar d    6
  // d d            7
  // d bar d        8
  // bar d bar d    9
  const int nPairs = 10;
  for (int i = 0; i < nPairs; ++i) {
    closeRejection.push_back(false);
    pairQA.push_back(0);
    NBins.push_back(1000);
    kMin.push_back(0.);
    kMax.push_back(1.);
  }

  closeRejection[0] = true; // pp
  closeRejection[2] = true; // pd
  closeRejection[4] = true; // barp barp
  closeRejection[6] = true; // barp bard
  closeRejection[7] = true; // dd
  closeRejection[9] = true; // bard bard
  pairQA[0] = 11;   // pp
  pairQA[2] = 11;   // pd
  pairQA[4] = 11;   // barp barp
  pairQA[6] = 11;   // barp bard
  pairQA[7] = 11;   // dd
  pairQA[9] = 11;   // bard bard

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

  AliFemtoDreamCollConfig *config = new AliFemtoDreamCollConfig("Femto", "Femto",
      false);
  config->SetZBins(ZVtxBins);
  config->SetMultBins(MultBins);
  config->SetMultBinning(true);
  config->SetPDGCodes(PDGParticles);
  config->SetNBinsHist(NBins);
  config->SetMinKRel(kMin);
  config->SetMaxKRel(kMax);
  config->SetClosePairRejection(closeRejection);
  config->SetPtQA(true);
  config->SetExtendedQAPairs(pairQA);
  config->SetDeltaEtaMax(0.017);  // and here you set the actual values
  config->SetDeltaPhiMax(0.017);  // and here you set the actual values
  //Here we set the mixing depth.
  config->SetMixingDepth(10);

  config->SetmTBins(mTBins);
  config->SetDomTMultBinning(true);
  config->SetmTBinning(true);

  if (isMC) {
    config->SetMomentumResolution(true);
  }
  if (fullBlastQA) {
    config->SetkTBinning(true);
    config->SetPtQA(true);
  }

  if (!fullBlastQA) {
    config->SetMinimalBookingME(true);
    config->SetMinimalBookingSample(true);
  }

  if (RefMult08) {
    config->SetMultiplicityEstimator(AliFemtoDreamEvent::kRef08);
  }

  if (Systematic) {
    if (suffix == "1") {
      TrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);
      AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);

      TrackCuts->SetNClsTPC(90);
      AntiTrackCuts->SetNClsTPC(90);

    } else if (suffix == "2") {
      TrackCuts->SetPtRange(0.6, 4.05);
      AntiTrackCuts->SetPtRange(0.6, 4.05);

      TrackCuts->SetEtaRange(-0.77, 0.77);
      AntiTrackCuts->SetEtaRange(-0.77, 0.77);

      config->SetDeltaEtaMax(0.019);
      config->SetDeltaPhiMax(0.019);

    } else if (suffix == "3") {
      TrackCuts->SetEtaRange(-0.77, 0.77);
      AntiTrackCuts->SetEtaRange(-0.77, 0.77);

      TrackCuts->SetNClsTPC(70);
      AntiTrackCuts->SetNClsTPC(70);

    } else if (suffix == "4") {
      TrackCuts->SetPtRange(0.4, 4.05);
      AntiTrackCuts->SetPtRange(0.4, 4.05);

      TrackCuts->SetEtaRange(-0.77, 0.77);
      AntiTrackCuts->SetEtaRange(-0.77, 0.77);

      TrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);
      AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);

      config->SetDeltaEtaMax(0.019);
      config->SetDeltaPhiMax(0.019);

    } else if (suffix == "5") {
      TrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);
      AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);

      TrackCuts->SetNClsTPC(90);
      AntiTrackCuts->SetNClsTPC(90);

    } else if (suffix == "6") {
      TrackCuts->SetEtaRange(-0.85, 0.85);
      AntiTrackCuts->SetEtaRange(-0.85, 0.85);

      TrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);
      AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);

      config->SetDeltaEtaMax(0.019);
      config->SetDeltaPhiMax(0.019);

    } else if (suffix == "7") {
      TrackCuts->SetPtRange(0.4, 4.05);
      AntiTrackCuts->SetPtRange(0.4, 4.05);

      TrackCuts->SetEtaRange(-0.77, 0.77);
      AntiTrackCuts->SetEtaRange(-0.77, 0.77);

      TrackCuts->SetNClsTPC(90);
      AntiTrackCuts->SetNClsTPC(90);

      config->SetDeltaEtaMax(0.019);
      config->SetDeltaPhiMax(0.019);

    } else if (suffix == "8") {

    } else if (suffix == "9") {
      TrackCuts->SetPtRange(0.4, 4.05);
      AntiTrackCuts->SetPtRange(0.4, 4.05);

      TrackCuts->SetEtaRange(-0.85, 0.85);
      AntiTrackCuts->SetEtaRange(-0.85, 0.85);

    } else if (suffix == "10") {
      TrackCuts->SetEtaRange(-0.77, 0.77);
      AntiTrackCuts->SetEtaRange(-0.77, 0.77);

      TrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);
      AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);

    } else if (suffix == "11") {
      TrackCuts->SetPtRange(0.4, 4.05);
      AntiTrackCuts->SetPtRange(0.4, 4.05);

      TrackCuts->SetEtaRange(-0.85, 0.85);
      AntiTrackCuts->SetEtaRange(-0.85, 0.85);

      TrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);
      AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);

      config->SetDeltaEtaMax(0.019);
      config->SetDeltaPhiMax(0.019);

    } else if (suffix == "12") {
      TrackCuts->SetEtaRange(-0.85, 0.85);
      AntiTrackCuts->SetEtaRange(-0.85, 0.85);

      TrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);
      AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);

      TrackCuts->SetNClsTPC(90);
      AntiTrackCuts->SetNClsTPC(90);

    } else if (suffix == "13") {
      TrackCuts->SetPtRange(0.4, 4.05);
      AntiTrackCuts->SetPtRange(0.4, 4.05);

      TrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);
      AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);

      config->SetDeltaEtaMax(0.019);
      config->SetDeltaPhiMax(0.019);

    } else if (suffix == "14") {
      TrackCuts->SetPtRange(0.6, 4.05);
      AntiTrackCuts->SetPtRange(0.6, 4.05);

      TrackCuts->SetEtaRange(-0.85, 0.85);
      AntiTrackCuts->SetEtaRange(-0.85, 0.85);

      TrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);
      AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);

      TrackCuts->SetNClsTPC(90);
      AntiTrackCuts->SetNClsTPC(90);

    } else if (suffix == "15") {
      TrackCuts->SetPtRange(0.4, 4.05);
      AntiTrackCuts->SetPtRange(0.4, 4.05);

      TrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);
      AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);

      TrackCuts->SetNClsTPC(70);
      AntiTrackCuts->SetNClsTPC(70);

      config->SetDeltaEtaMax(0.019);
      config->SetDeltaPhiMax(0.019);

    } else if (suffix == "16") {
      TrackCuts->SetPtRange(0.4, 4.05);
      AntiTrackCuts->SetPtRange(0.4, 4.05);

      TrackCuts->SetEtaRange(-0.85, 0.85);
      AntiTrackCuts->SetEtaRange(-0.85, 0.85);

      TrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);
      AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);

      config->SetDeltaEtaMax(0.019);
      config->SetDeltaPhiMax(0.019);

    } else if (suffix == "17") {
      TrackCuts->SetPtRange(0.4, 4.05);
      AntiTrackCuts->SetPtRange(0.4, 4.05);

      TrackCuts->SetEtaRange(-0.85, 0.85);
      AntiTrackCuts->SetEtaRange(-0.85, 0.85);

      TrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);
      AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);

      TrackCuts->SetNClsTPC(90);
      AntiTrackCuts->SetNClsTPC(90);

    } else if (suffix == "18") {
      TrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);
      AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);

      TrackCuts->SetNClsTPC(70);
      AntiTrackCuts->SetNClsTPC(70);

      config->SetDeltaEtaMax(0.019);
      config->SetDeltaPhiMax(0.019);

    } else if (suffix == "19") {
      TrackCuts->SetPtRange(0.4, 4.05);
      AntiTrackCuts->SetPtRange(0.4, 4.05);

      TrackCuts->SetNClsTPC(90);
      AntiTrackCuts->SetNClsTPC(90);

      config->SetDeltaEtaMax(0.019);
      config->SetDeltaPhiMax(0.019);

    } else if (suffix == "20") {

      TrackCuts->SetEtaRange(-0.85, 0.85);
      AntiTrackCuts->SetEtaRange(-0.85, 0.85);

      TrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);
      AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);

      TrackCuts->SetNClsTPC(70);
      AntiTrackCuts->SetNClsTPC(70);

      config->SetDeltaEtaMax(0.019);
      config->SetDeltaPhiMax(0.019);

    } else if (suffix == "21") {
      TrackCuts->SetPtRange(0.4, 4.05);
      AntiTrackCuts->SetPtRange(0.4, 4.05);

      TrackCuts->SetEtaRange(-0.77, 0.77);
      AntiTrackCuts->SetEtaRange(-0.77, 0.77);

      TrackCuts->SetNClsTPC(70);
      AntiTrackCuts->SetNClsTPC(70);

      config->SetDeltaEtaMax(0.019);
      config->SetDeltaPhiMax(0.019);

    } else if (suffix == "22") {
      TrackCuts->SetEtaRange(-0.77, 0.77);
      AntiTrackCuts->SetEtaRange(-0.77, 0.77);

      TrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);
      AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);

    } else if (suffix == "23") {
      TrackCuts->SetEtaRange(-0.77, 0.77);
      AntiTrackCuts->SetEtaRange(-0.77, 0.77);

      TrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);
      AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);

      TrackCuts->SetNClsTPC(90);
      AntiTrackCuts->SetNClsTPC(90);

      config->SetDeltaEtaMax(0.019);
      config->SetDeltaPhiMax(0.019);

    } else if (suffix == "24") {
      TrackCuts->SetPtRange(0.4, 4.05);
      AntiTrackCuts->SetPtRange(0.4, 4.05);

      TrackCuts->SetEtaRange(-0.83, 0.83);
      AntiTrackCuts->SetEtaRange(-0.83, 0.83);

      config->SetDeltaEtaMax(0.019);
      config->SetDeltaPhiMax(0.019);

    } else if (suffix == "25") {
      TrackCuts->SetPtRange(0.6, 4.05);
      AntiTrackCuts->SetPtRange(0.6, 4.05);

      TrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);
      AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);

      config->SetDeltaEtaMax(0.019);
      config->SetDeltaPhiMax(0.019);

    } else if (suffix == "26") {
      TrackCuts->SetEtaRange(-0.77, 0.77);
      AntiTrackCuts->SetEtaRange(-0.77, 0.77);

      TrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);
      AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);

      TrackCuts->SetNClsTPC(90);
      AntiTrackCuts->SetNClsTPC(90);

    } else if (suffix == "27") {
      TrackCuts->SetNClsTPC(90);
      AntiTrackCuts->SetNClsTPC(90);

    } else if (suffix == "28") {
      TrackCuts->SetPtRange(0.4, 4.05);
      AntiTrackCuts->SetPtRange(0.4, 4.05);

      TrackCuts->SetEtaRange(-0.77, 0.77);
      AntiTrackCuts->SetEtaRange(-0.77, 0.77);

      TrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);
      AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);

      config->SetDeltaEtaMax(0.019);
      config->SetDeltaPhiMax(0.019);

    } else if (suffix == "29") {
      TrackCuts->SetEtaRange(-0.85, 0.85);
      AntiTrackCuts->SetEtaRange(-0.85, 0.85);

      config->SetDeltaEtaMax(0.019);
      config->SetDeltaPhiMax(0.019);

    } else if (suffix == "30") {
      TrackCuts->SetPtRange(0.4, 4.05);
      AntiTrackCuts->SetPtRange(0.4, 4.05);

      TrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);
      AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);

      TrackCuts->SetNClsTPC(70);
      AntiTrackCuts->SetNClsTPC(70);

      config->SetDeltaEtaMax(0.019);
      config->SetDeltaPhiMax(0.019);

    } else if (suffix == "31") {

      TrackCuts->SetEtaRange(-0.77, 0.77);
      AntiTrackCuts->SetEtaRange(-0.77, 0.77);

      TrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);
      AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);

      TrackCuts->SetNClsTPC(90);
      AntiTrackCuts->SetNClsTPC(90);

    } else if (suffix == "32") {
      TrackCuts->SetEtaRange(-0.77, 0.77);
      AntiTrackCuts->SetEtaRange(-0.77, 0.77);

      config->SetDeltaEtaMax(0.019);
      config->SetDeltaPhiMax(0.019);

    } else if (suffix == "33") {
      TrackCuts->SetEtaRange(-0.77, 0.77);
      AntiTrackCuts->SetEtaRange(-0.77, 0.77);

      TrackCuts->SetNClsTPC(70);
      AntiTrackCuts->SetNClsTPC(70);

      config->SetDeltaEtaMax(0.019);
      config->SetDeltaPhiMax(0.019);

    } else if (suffix == "34") {
      TrackCuts->SetPtRange(0.4, 4.05);
      AntiTrackCuts->SetPtRange(0.4, 4.05);

      TrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);
      AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);

      TrackCuts->SetNClsTPC(70);
      AntiTrackCuts->SetNClsTPC(70);

      config->SetDeltaEtaMax(0.019);
      config->SetDeltaPhiMax(0.019);

    } else if (suffix == "35") {
      TrackCuts->SetPtRange(0.6, 4.05);
      AntiTrackCuts->SetPtRange(0.6, 4.05);

      TrackCuts->SetEtaRange(-0.85, 0.85);
      AntiTrackCuts->SetEtaRange(-0.85, 0.85);

      TrackCuts->SetNClsTPC(70);
      AntiTrackCuts->SetNClsTPC(70);

      config->SetDeltaEtaMax(0.019);
      config->SetDeltaPhiMax(0.019);

    } else if (suffix == "36") {
      TrackCuts->SetEtaRange(-0.85, 0.85);
      AntiTrackCuts->SetEtaRange(-0.85, 0.85);

      TrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);
      AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);

      TrackCuts->SetNClsTPC(90);
      AntiTrackCuts->SetNClsTPC(90);

      config->SetDeltaEtaMax(0.019);
      config->SetDeltaPhiMax(0.019);

    } else if (suffix == "37") {
      TrackCuts->SetEtaRange(-0.77, 0.77);
      AntiTrackCuts->SetEtaRange(-0.77, 0.77);

      TrackCuts->SetNClsTPC(90);
      AntiTrackCuts->SetNClsTPC(90);

      config->SetDeltaEtaMax(0.019);
      config->SetDeltaPhiMax(0.019);

    } else if (suffix == "38") {
      TrackCuts->SetPtRange(0.4, 4.05);
      AntiTrackCuts->SetPtRange(0.4, 4.05);

      TrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);
      AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);

      TrackCuts->SetNClsTPC(90);
      AntiTrackCuts->SetNClsTPC(90);

      config->SetDeltaEtaMax(0.019);
      config->SetDeltaPhiMax(0.019);

    } else if (suffix == "39") {
      TrackCuts->SetPtRange(0.6, 4.05);
      AntiTrackCuts->SetPtRange(0.6, 4.05);

      TrackCuts->SetNClsTPC(70);
      AntiTrackCuts->SetNClsTPC(70);

    } else if (suffix == "40") {
      TrackCuts->SetPtRange(0.4, 4.05);
      AntiTrackCuts->SetPtRange(0.4, 4.05);

      TrackCuts->SetEtaRange(-0.77, 0.77);
      AntiTrackCuts->SetEtaRange(-0.77, 0.77);

      TrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);
      AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);

      config->SetDeltaEtaMax(0.019);
      config->SetDeltaPhiMax(0.019);

    } else if (suffix == "41") {
      TrackCuts->SetPtRange(0.6, 4.05);
      AntiTrackCuts->SetPtRange(0.6, 4.05);

      TrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);
      AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);

    } else if (suffix == "42") {
      TrackCuts->SetPtRange(0.6, 4.05);
      AntiTrackCuts->SetPtRange(0.6, 4.05);

      TrackCuts->SetEtaRange(-0.85, 0.85);
      AntiTrackCuts->SetEtaRange(-0.85, 0.85);

      TrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);
      AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);

      config->SetDeltaEtaMax(0.019);
      config->SetDeltaPhiMax(0.019);

    } else if (suffix == "43") {
      TrackCuts->SetEtaRange(-0.85, 0.85);
      AntiTrackCuts->SetEtaRange(-0.85, 0.85);

      TrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);
      AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);

      config->SetDeltaEtaMax(0.019);
      config->SetDeltaPhiMax(0.019);

    }
  }

//Only pT cut variatons
  if (SystematicpTCutVariation) {
    if (suffix == "1") {
      TrackCutsDeuteron->SetPtRange(0.5, 1.4);
      AntiTrackCutsDeuteron->SetPtRange(0.5, 1.4);
      //TrackCutsDeuteronNoTOF->SetPtRange(0.5, 1.4);
      //AntiTrackCutsDeuteronNoTOF->SetPtRange(0.5, 1.4);
     // config->SetDeltaEtaMax(0.005);
      //config->SetDeltaPhiMax(0.005);
    } else if (suffix == "2") {
      TrackCutsDeuteron->SetPtRange(0.5, 1.5);
      AntiTrackCutsDeuteron->SetPtRange(0.5, 1.5);

      //config->SetDeltaEtaMax(0.007);
     // config->SetDeltaPhiMax(0.007);

    } else if (suffix == "3") {
      TrackCutsDeuteron->SetPtRange(0.5, 1.7);
      AntiTrackCutsDeuteron->SetPtRange(0.5, 1.7);
      //config->SetDeltaEtaMax(0.009);
      //config->SetDeltaPhiMax(0.009);
    } else if (suffix == "4") {

      TrackCutsDeuteron->SetPtRange(0.5, 1.9);
      AntiTrackCutsDeuteron->SetPtRange(0.5, 1.9);
      //config->SetDeltaEtaMax(0.011);
      //config->SetDeltaPhiMax(0.011);

    } else if (suffix == "5") {
      TrackCutsDeuteron->SetPtRange(0.5, 2.1);
      AntiTrackCutsDeuteron->SetPtRange(0.5, 2.1);
     // config->SetDeltaEtaMax(0.013);
      //config->SetDeltaPhiMax(0.013);

    } else if (suffix == "6") {


      TrackCutsDeuteron->SetPtRange(0.5, 2.3);
      AntiTrackCutsDeuteron->SetPtRange(0.5, 2.3);
     // config->SetDeltaEtaMax(0.015);
      //config->SetDeltaPhiMax(0.015);

    } else if (suffix == "7") {
      TrackCutsDeuteron->SetPtRange(0.5, 2.7);
      AntiTrackCutsDeuteron->SetPtRange(0.5, 2.7);
      //config->SetDeltaEtaMax(0.019);
     // config->SetDeltaPhiMax(0.019);

    } else if (suffix == "8") {
      TrackCutsDeuteron->SetPtRange(0.5, 3.0);
      AntiTrackCutsDeuteron->SetPtRange(0.5, 3.0);

      //config->SetDeltaEtaMax(0.021);
      ///config->SetDeltaPhiMax(0.021);
    } else if (suffix == "9") {

      TrackCutsDeuteron->SetPtRange(0.8, 2.5);
      AntiTrackCutsDeuteron->SetPtRange(0.8, 2.5);

      //config->SetDeltaEtaMax(0.023);
     // config->SetDeltaPhiMax(0.023);

    } else if (suffix == "10") {//A check with Micheals Jung task
      TrackCutsDeuteron->SetFilterBit(128);
      AntiTrackCutsDeuteron->SetFilterBit(128);
      TrackCutsDeuteron->SetPtRange(0.8, 2.5);
      AntiTrackCutsDeuteron->SetPtRange(0.8, 2.5);
      //A check for Michaels task
      //config->SetDeltaEtaMax(0.023);
     // config->SetDeltaPhiMax(0.023);

    }else if (suffix == "11") {
      TrackCutsDeuteron->SetFilterBit(128);
      AntiTrackCutsDeuteron->SetFilterBit(128);
      TrackCutsDeuteron->SetPtRange(0.4, 2.5);
      AntiTrackCutsDeuteron->SetPtRange(0.4, 2.5);
     // config->SetDeltaEtaMax(0.025);
     // config->SetDeltaPhiMax(0.025);

    } else if (suffix == "12") {
      TrackCutsDeuteron->SetFilterBit(128);
      AntiTrackCutsDeuteron->SetFilterBit(128);
      TrackCutsDeuteron->SetPtRange(0.5, 1.4);
      AntiTrackCutsDeuteron->SetPtRange(0.5, 1.4);
    } else if (suffix == "13") {
      TrackCutsDeuteron->SetFilterBit(128);
      AntiTrackCutsDeuteron->SetFilterBit(128);
      TrackCutsDeuteron->SetPtRange(0.5, 1.7);
      AntiTrackCutsDeuteron->SetPtRange(0.5, 1.7);
    } else if (suffix == "14") {
      TrackCutsDeuteron->SetFilterBit(128);
      AntiTrackCutsDeuteron->SetFilterBit(128);
      TrackCutsDeuteron->SetPtRange(0.5, 1.7);
      AntiTrackCutsDeuteron->SetPtRange(0.5, 1.7);
    } else if (suffix == "15") {
      TrackCutsDeuteron->SetFilterBit(128);
      AntiTrackCutsDeuteron->SetFilterBit(128);
      TrackCutsDeuteron->SetPtRange(0.5, 2.1);
      AntiTrackCutsDeuteron->SetPtRange(0.5, 2.1);
    }
  }

  AliAnalysisTaskNanoPt *task = new AliAnalysisTaskNanoPt(
    "AliAnalysisTaskNanoPt", isMC);

  if (trigger == "kINT7") {
    task->SelectCollisionCandidates(AliVEvent::kINT7);
    std::cout << "Added kINT7 Trigger \n";
  } else if (trigger == "kHM") {
    task->SelectCollisionCandidates(AliVEvent::kHighMultV0);
    std::cout << "Added kHighMult Trigger \n";
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

  task->SetEventCuts(evtCuts);
  task->SetProtonCuts(TrackCuts);
  task->SetAntiProtonCuts(AntiTrackCuts);
  task->SetDeuteronCuts(TrackCutsDeuteron);
  task->SetAntiDeuteronCuts(AntiTrackCutsDeuteron);
  task->SetDeuteronCutsNoTOF(TrackCutsDeuteronNoTOF);
  task->SetAntiDeuteronCutsNoTOF(AntiTrackCutsDeuteronNoTOF);
  task->SetCollectionConfig(config);
  task->SetUseDumpster(DumpPdApAd);
  task->SetMCTruth(fIsMCTruth);
  mgr->AddTask(task);

  TString file = AliAnalysisManager::GetCommonFileName();
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  mgr->ConnectInput(task, 0, cinput);

  TString addon = "";

  if (trigger == "kINT7") {
    addon += "MB";
  } else if (trigger == "kHM") {
    addon += "HM";
  }

  TString EvtCutsName = Form("%sEvtCuts%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *coutputEvtCuts = mgr->CreateContainer(
        EvtCutsName.Data(), TList::Class(), AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), EvtCutsName.Data()));
  mgr->ConnectOutput(task, 1, coutputEvtCuts);

  TString TrackCutsName = Form("%sProton%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *coutputTrkCuts = mgr->CreateContainer(
        TrackCutsName.Data(), TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), TrackCutsName.Data()));
  mgr->ConnectOutput(task, 2, coutputTrkCuts);

  TString AntiTrackCutsName = Form("%sAntiProton%s", addon.Data(),
                                   suffix.Data());
  AliAnalysisDataContainer *coutputAntiTrkCuts = mgr->CreateContainer(
        AntiTrackCutsName.Data(), TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), AntiTrackCutsName.Data()));
  mgr->ConnectOutput(task, 3, coutputAntiTrkCuts);

  TString TrackCutsDeuteronName = Form("%sDeuteron%s", addon.Data(),
                                       suffix.Data());
  AliAnalysisDataContainer *coutputTrkCutsDeuteron = mgr->CreateContainer(
        TrackCutsDeuteronName.Data(), TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), TrackCutsDeuteronName.Data()));
  mgr->ConnectOutput(task, 4, coutputTrkCutsDeuteron);

  TString AntiTrackCutsDeuteronName = Form("%sAntiDeuteron%s", addon.Data(),
                                      suffix.Data());
  AliAnalysisDataContainer *coutputAntiTrkCutsDeuteron = mgr->CreateContainer(
        AntiTrackCutsDeuteronName.Data(), TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), AntiTrackCutsDeuteronName.Data()));
  mgr->ConnectOutput(task, 5, coutputAntiTrkCutsDeuteron);
//============== NoTOF STUFF========================================
  TString TrackCutsDeuteronNoTOFName = Form("%sDeuteronNoTOF%s", addon.Data(),
                                       suffix.Data());
  AliAnalysisDataContainer *coutputTrkCutsDeuteronNoTOF = mgr->CreateContainer(
        TrackCutsDeuteronNoTOFName.Data(), TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), TrackCutsDeuteronNoTOFName.Data()));
  mgr->ConnectOutput(task, 6, coutputTrkCutsDeuteronNoTOF);

  TString AntiTrackCutsDeuteronNoTOFName = Form("%sAntiDeuteronNoTOF%s",
      addon.Data(), suffix.Data());
  AliAnalysisDataContainer *coutputAntiTrkCutsDeuteronNoTOF = mgr
      ->CreateContainer(
        AntiTrackCutsDeuteronNoTOFName.Data(), TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), AntiTrackCutsDeuteronNoTOFName.Data()));
  mgr->ConnectOutput(task, 7, coutputAntiTrkCutsDeuteronNoTOF);

  AliAnalysisDataContainer *coutputResults;
  TString ResultsName = Form("%sResults%s", addon.Data(), suffix.Data());
  coutputResults = mgr->CreateContainer(
                     //@suppress("Invalid arguments") it works ffs
                     ResultsName.Data(),
                     TList::Class(), AliAnalysisManager::kOutputContainer,
                     Form("%s:%s", file.Data(), ResultsName.Data()));
  mgr->ConnectOutput(task, 8, coutputResults);

  AliAnalysisDataContainer *coutputResultsQA;
  TString ResultsQAName = Form("%sResultsQA%s", addon.Data(), suffix.Data());
  coutputResultsQA = mgr->CreateContainer(
                       //@suppress("Invalid arguments") it works ffs
                       ResultsQAName.Data(),
                       TList::Class(),
                       AliAnalysisManager::kOutputContainer,
                       Form("%s:%s", file.Data(), ResultsQAName.Data()));
  mgr->ConnectOutput(task, 9, coutputResultsQA);

  AliAnalysisDataContainer *coutputDumpster;
  TString DumpsterName = Form("%sDumpster%s", addon.Data(), suffix.Data());
  coutputDumpster = mgr->CreateContainer(
                      //@suppress("Invalid arguments") it works ffs
                      DumpsterName.Data(),
                      TList::Class(), AliAnalysisManager::kOutputContainer,
                      Form("%s:%s", file.Data(), DumpsterName.Data()));
  mgr->ConnectOutput(task, 10, coutputDumpster);

  if (isMC) {
    AliAnalysisDataContainer *coutputTrkCutsMC;
    TString TrkCutsMCName = Form("%sProtonMC%s", addon.Data(), suffix.Data());
    coutputTrkCutsMC = mgr->CreateContainer(
                         TrkCutsMCName.Data(), TList::Class(),
                         AliAnalysisManager::kOutputContainer,
                         Form("%s:%s", file.Data(), TrkCutsMCName.Data()));
    mgr->ConnectOutput(task, 11, coutputTrkCutsMC);

    AliAnalysisDataContainer *coutputAntiTrkCutsMC;
    TString AntiTrkCutsMCName = Form("%sAntiProtonMC%s", addon.Data(),
                                     suffix.Data());
    coutputAntiTrkCutsMC = mgr->CreateContainer(
                             //@suppress("Invalid arguments") it works ffs
                             AntiTrkCutsMCName.Data(),
                             TList::Class(),
                             AliAnalysisManager::kOutputContainer,
                             Form("%s:%s", file.Data(), AntiTrkCutsMCName.Data()));
    mgr->ConnectOutput(task, 12, coutputAntiTrkCutsMC);

    AliAnalysisDataContainer *coutputv0CutsMC;
    TString v0CutsMCName = Form("%sDeuteronMC%s", addon.Data(), suffix.Data());
    coutputv0CutsMC = mgr->CreateContainer(
                        //@suppress("Invalid arguments") it works ffs
                        v0CutsMCName.Data(),
                        TList::Class(),
                        AliAnalysisManager::kOutputContainer,
                        Form("%s:%s", file.Data(), v0CutsMCName.Data()));
    mgr->ConnectOutput(task, 13, coutputv0CutsMC);

    AliAnalysisDataContainer *coutputAntiv0CutsMC;
    TString Antiv0CutsMCName = Form("%sAntiDeuteronMC%s", addon.Data(),
                                    suffix.Data());
    coutputAntiv0CutsMC = mgr->CreateContainer(
                            //@suppress("Invalid arguments") it works ffs
                            Antiv0CutsMCName.Data(),
                            TList::Class(),
                            AliAnalysisManager::kOutputContainer,
                            Form("%s:%s", file.Data(), Antiv0CutsMCName.Data()));
    mgr->ConnectOutput(task, 14, coutputAntiv0CutsMC);

  }

  return task;
}

