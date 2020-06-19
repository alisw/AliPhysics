#include <vector>
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliFemtoDreamEventCuts.h"
#include "AliFemtoDreamTrackCuts.h"
#include "AliFemtoDreamCollConfig.h"
#include "AliAnalysisTaskNanoPPCoalescence.h"
AliAnalysisTaskSE* AddTaskNanoPPCoalescence( bool isMC = true, bool fIsMCTruth = true,
                                  TString trigger = "kINT7", //2
                                  bool DCAPlots = false, //3
                                  bool CombSigma = false, //4
                                  bool ContributionSplitting = false, //5,
                                  bool fullBlastQA = true, bool RefMult08 = true,
                                  bool Systematic = false,
                                  const char *cutVariation = "0") {

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

  // Event Cuts
  AliFemtoDreamEventCuts *evtCuts = AliFemtoDreamEventCuts::StandardCutsRun2();
  evtCuts->CleanUpMult(false, false, false, true);


  //Track Cuts for p/barp
  AliFemtoDreamTrackCuts *TrackCuts = AliFemtoDreamTrackCuts::PrimProtonCuts(
                                        isMC, true, CombSigma, ContributionSplitting);
  TrackCuts->SetMinimalBooking(Systematic);
  TrackCuts->SetCutCharge(1);

  AliFemtoDreamTrackCuts *AntiTrackCuts = AliFemtoDreamTrackCuts::PrimProtonCuts(
      isMC, true, CombSigma, ContributionSplitting);
  AntiTrackCuts->SetMinimalBooking(Systematic);
  AntiTrackCuts->SetCutCharge(-1);

  std::vector<int> PDGParticles;
  PDGParticles.push_back(2212);
  PDGParticles.push_back(2212);

  std::vector<int> NBins;
  std::vector<float> kMin;
  std::vector<float> kMax;
  std::vector<bool> closeRejection;
  std::vector<float> mTBins = { 1.14, 1.26, 999.};
  std::vector<int> pairQA;
  //pairs:
  // pp             0
  // p bar p        1
  // bar p p        2
  // bar p bar p    3

  const int nPairs = 3;
  for (int i = 0; i < nPairs; ++i) {
    closeRejection.push_back(false);
    pairQA.push_back(0);
    NBins.push_back(1000);
    kMin.push_back(0.);
    kMax.push_back(1.);
  }

  closeRejection[0] = true; // pp
  closeRejection[3] = true; // barp barp
  pairQA[0] = 11;   // pp
  pairQA[3] = 11;   // barp barp

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

  AliFemtoDreamCollConfig *config = new AliFemtoDreamCollConfig("Femto","Femto",false);
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


  AliAnalysisTaskNanoPPCoalescence *task = new AliAnalysisTaskNanoPPCoalescence(
    "AliAnalysisTaskNanoPPCoalescence", isMC);

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
  task->SetCollectionConfig(config);
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

  AliAnalysisDataContainer *coutputResults;
  TString ResultsName = Form("%sResults%s", addon.Data(), suffix.Data());
  coutputResults = mgr->CreateContainer(
                     //@suppress("Invalid arguments") it works ffs
                     ResultsName.Data(),
                     TList::Class(), AliAnalysisManager::kOutputContainer,
                     Form("%s:%s", file.Data(), ResultsName.Data()));
  mgr->ConnectOutput(task, 4, coutputResults);

  AliAnalysisDataContainer *coutputResultsQA;
  TString ResultsQAName = Form("%sResultsQA%s", addon.Data(), suffix.Data());
  coutputResultsQA = mgr->CreateContainer(
                       //@suppress("Invalid arguments") it works ffs
                       ResultsQAName.Data(),
                       TList::Class(),
                       AliAnalysisManager::kOutputContainer,
                       Form("%s:%s", file.Data(), ResultsQAName.Data()));
  mgr->ConnectOutput(task, 5, coutputResultsQA);

  if (isMC) {
    AliAnalysisDataContainer *coutputTrkCutsMC;
    TString TrkCutsMCName = Form("%sProtonMC%s", addon.Data(), suffix.Data());
    coutputTrkCutsMC = mgr->CreateContainer(
                         TrkCutsMCName.Data(), TList::Class(),
                         AliAnalysisManager::kOutputContainer,
                         Form("%s:%s", file.Data(), TrkCutsMCName.Data()));
    mgr->ConnectOutput(task, 6, coutputTrkCutsMC);

    AliAnalysisDataContainer *coutputAntiTrkCutsMC;
    TString AntiTrkCutsMCName = Form("%sAntiProtonMC%s", addon.Data(),
                                     suffix.Data());
    coutputAntiTrkCutsMC = mgr->CreateContainer(
                             //@suppress("Invalid arguments") it works ffs
                             AntiTrkCutsMCName.Data(),
                             TList::Class(),
                             AliAnalysisManager::kOutputContainer,
                             Form("%s:%s", file.Data(), AntiTrkCutsMCName.Data()));
    mgr->ConnectOutput(task, 7, coutputAntiTrkCutsMC);

  }

  return task;
}

