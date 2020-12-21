#include <vector>
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliFemtoDreamEventCuts.h"
#include "AliFemtoDreamTrackCuts.h"
#include "AliFemtoDreamCollConfig.h"
#include "AliAnalysisTaskNanoPt.h"
AliAnalysisTaskSE* AddTaskNanoPt(bool isMC = true, //1
                                 bool IsMCTruth = true, //2
                                 TString trigger = "kINT7", //3
                                 bool DCAPlots = false, //4
                                 bool CombSigma = false, //5
                                 bool ContributionSplitting = false, //6,
                                 bool DumpPdApAd = true, //7
                                 bool fullBlastQA = true, //8,
                                 bool MinBook = true, //8,
                                 bool RefMult08 = true, //9
                                 bool SidebandStudy = true, //10
                                 bool Systematic = false, //11
                                 bool SystematicpTCutVariation = false, //12
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
  //========= Init subtasks and start analyis ============================
  // Event Cuts
  AliFemtoDreamEventCuts *evtCuts = AliFemtoDreamEventCuts::StandardCutsRun2();
  evtCuts->CleanUpMult(false, false, false, true);

  // =====================================================================
  //Proton track Cuts----------------------------------------------------------------------------
  AliFemtoDreamTrackCuts *TrackCuts = AliFemtoDreamTrackCuts::PrimProtonCuts(
                                        isMC, true, CombSigma, ContributionSplitting);
  TrackCuts->SetMinimalBooking(MinBook);
  TrackCuts->SetCutCharge(1);
  //Antiproton track Cuts-------------------------------------------------------------------------
  AliFemtoDreamTrackCuts *AntiTrackCuts = AliFemtoDreamTrackCuts::PrimProtonCuts(
      isMC, true, CombSigma, ContributionSplitting);
  AntiTrackCuts->SetMinimalBooking(MinBook);
  AntiTrackCuts->SetCutCharge(-1);
  //deuteron track cuts----------------------------------------------------------------------------
  AliFemtoDreamTrackCuts *TrackCutsDeuteron =
    AliFemtoDreamTrackCuts::PrimDeuteronCuts(isMC, true, CombSigma,
        ContributionSplitting);
  TrackCutsDeuteron->SetMinimalBooking(MinBook);
  TrackCutsDeuteron->SetCutCharge(1);
  //Antideuteron track cuts----------------------------------------------------------------------------
  AliFemtoDreamTrackCuts *AntiTrackCutsDeuteron =
    AliFemtoDreamTrackCuts::PrimDeuteronCuts(isMC, true, CombSigma,
        ContributionSplitting);
  AntiTrackCutsDeuteron->SetMinimalBooking(MinBook);
  AntiTrackCutsDeuteron->SetCutCharge(-1);
/////////////////////For no NSigmaTOF information///
// =====================================================================
  //deuteron track cuts----------------------------------------------------------------------------
  AliFemtoDreamTrackCuts *TrackCutsDeuteronNoTOF =
    AliFemtoDreamTrackCuts::PrimDeuteronCuts(isMC, true, CombSigma,
        ContributionSplitting);
  TrackCutsDeuteronNoTOF->SetMinimalBooking(MinBook);
  TrackCutsDeuteronNoTOF->SetCutCharge(1);
  TrackCutsDeuteronNoTOF->SetPID(AliPID::kDeuteron, 999.);
  //Antideuteron track cuts----------------------------------------------------------------------------
  AliFemtoDreamTrackCuts *AntiTrackCutsDeuteronNoTOF =
    AliFemtoDreamTrackCuts::PrimDeuteronCuts(isMC, true, CombSigma,
        ContributionSplitting);
  AntiTrackCutsDeuteronNoTOF->SetMinimalBooking(MinBook);
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
  pairQA[0] = 11; // pp
  pairQA[2] = 11; // pd
  pairQA[4] = 11; // barp barp
  pairQA[6] = 11; // barp bard
  pairQA[7] = 11; // dd
  pairQA[9] = 11; // bard bard

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
  //Cut on TOF mass of Deuteron and Anti-deuteron
  if (SidebandStudy) {

    TrackCutsDeuteron->SetPtRange(1.4, 4.05);
    AntiTrackCutsDeuteron->SetPtRange(1.4, 4.05);
    TrackCutsDeuteron->SetPID(AliPID::kDeuteron, 1.4);
    AntiTrackCutsDeuteron->SetPID(AliPID::kDeuteron, 1.4);
    TrackCuts->SetPtRange(0.5, 4.05);
    AntiTrackCuts->SetPtRange(0.5, 4.05);

    TrackCutsDeuteron->SetPlotTOFMassSq(true);
    AntiTrackCutsDeuteron->SetPlotTOFMassSq(true);
    TrackCutsDeuteron->SetCutTOFInvMass(true);
    AntiTrackCutsDeuteron->SetCutTOFInvMass(true);
    //3.5179128721 Nominal mass^2 peak
    if (suffix == "0") {
         TrackCutsDeuteron->SetPtRange(0.5, 1.4);
         AntiTrackCutsDeuteron->SetPtRange(0.5,1.4);
         TrackCutsDeuteron->SetPID(AliPID::kDeuteron, 1.4);
         AntiTrackCutsDeuteron->SetPID(AliPID::kDeuteron, 1.4);
         TrackCuts->SetPtRange(0.5, 4.05);
         AntiTrackCuts->SetPtRange(0.5, 4.05);
         TrackCutsDeuteron->SetCutTOFInvMass(false);
         AntiTrackCutsDeuteron->SetCutTOFInvMass(false);
      }else if (suffix == "1") {
        TrackCutsDeuteron->SetCutPeakTOFInvMass(0.155);
        AntiTrackCutsDeuteron->SetCutPeakTOFInvMass(0.155);
      } else if (suffix == "2") {
        TrackCutsDeuteron->SetCutPeakTOFInvMass(0.310);
        AntiTrackCutsDeuteron->SetCutPeakTOFInvMass(0.310);
      } else if (suffix == "3") {
        TrackCutsDeuteron->SetCutPeakTOFInvMass(0.21791);
        AntiTrackCutsDeuteron->SetCutPeakTOFInvMass(0.21791);
      } else if (suffix == "4") {
        TrackCutsDeuteron->SetCutPeakTOFInvMass(0.482);
        AntiTrackCutsDeuteron->SetCutPeakTOFInvMass(0.482);
      } else if (suffix == "5") {//LeftSideBand lower edge = 2.0GeV/c^2
        TrackCutsDeuteron->SetCutTOFMassForSB(2.0, 3.5179128721-0.310);//two sigma
        AntiTrackCutsDeuteron->SetCutTOFMassForSB(2.0, 3.5179128721-0.310);//two sigma
      } else if (suffix == "6") {
        TrackCutsDeuteron->SetCutTOFMassForSB(2.0, 3.5179128721-0.465);//three sigma
        AntiTrackCutsDeuteron->SetCutTOFMassForSB(2.0, 3.5179128721-0.465);//three sigma
      } else if (suffix == "7") {
        TrackCutsDeuteron->SetCutTOFMassForSB(2.8, 3.2);//Out of three sigma
        AntiTrackCutsDeuteron->SetCutTOFMassForSB(2.8, 3.2);//Out of three sigma
      }else if (suffix == "8") {//LeftSideBand lower edge = 2.5GeV/c^2
        TrackCutsDeuteron->SetCutTOFMassForSB(2.5, 3.5179128721-0.310);//two sigma
        AntiTrackCutsDeuteron->SetCutTOFMassForSB(2.5, 3.5179128721-0.310);//two sigma
      } else if (suffix == "9") {
        TrackCutsDeuteron->SetCutTOFMassForSB(2.5, 3.5179128721-0.465);//three sigma
        AntiTrackCutsDeuteron->SetCutTOFMassForSB(2.5, 3.5179128721-0.465);//three sigma
      } else if (suffix == "10") {
        TrackCutsDeuteron->SetCutTOFMassForSB(2.5, 2.9);//Out of three sigma
        AntiTrackCutsDeuteron->SetCutTOFMassForSB(2.5, 2.9);//Out of three sigma
      }else if (suffix == "11") {//RightSideBand upper edge = 5.5GeV/c^2
        TrackCutsDeuteron->SetCutTOFMassForSB(3.5179128721+0.310,5.5);//two sigma
        AntiTrackCutsDeuteron->SetCutTOFMassForSB(3.5179128721+0.310,5.5);//two sigma
      } else if (suffix == "12") {
        TrackCutsDeuteron->SetCutTOFMassForSB(3.5179128721+0.465,5.5);//three sigma
        AntiTrackCutsDeuteron->SetCutTOFMassForSB(3.5179128721+0.465,5.5);//three sigma
      } else if (suffix == "13") {
        TrackCutsDeuteron->SetCutTOFMassForSB(4.01,5.5);//Out of three sigma
        AntiTrackCutsDeuteron->SetCutTOFMassForSB(4.01,5.5);//Out of three sigma
      }else if (suffix == "14") {//RightSideBand upper edge = 5.0GeV/c^2
        TrackCutsDeuteron->SetCutTOFMassForSB(3.5179128721+0.310,5.0);//two sigma
        AntiTrackCutsDeuteron->SetCutTOFMassForSB(3.5179128721+0.310,5.0);//two sigma
      } else if (suffix == "15") {
        TrackCutsDeuteron->SetCutTOFMassForSB(3.5179128721+0.465,5.0);//three sigma
        AntiTrackCutsDeuteron->SetCutTOFMassForSB(3.5179128721+0.465,5.0);//three sigma
      } else if (suffix == "16") {
        TrackCutsDeuteron->SetCutTOFMassForSB(3.85,4.3);//Out of three sigma
        AntiTrackCutsDeuteron->SetCutTOFMassForSB(3.85,4.3);//Out of three sigma
      }
  }
  if (Systematic) {
    if (suffix == "1") {
      TrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);
      AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);
      TrackCuts->SetNClsTPC(90);
      AntiTrackCuts->SetNClsTPC(90);

      TrackCutsDeuteron->SetPID(AliPID::kDeuteron, 1.4, 2.7);
      AntiTrackCutsDeuteron->SetPID(AliPID::kDeuteron, 1.4, 2.7);
      TrackCutsDeuteron->SetNClsTPC(90);
      AntiTrackCutsDeuteron->SetNClsTPC(90);

    } else if (suffix == "2") {
      TrackCuts->SetPtRange(0.6, 2.5);
      AntiTrackCuts->SetPtRange(0.6, 2.5);
      TrackCuts->SetEtaRange(-0.77, 0.77);
      AntiTrackCuts->SetEtaRange(-0.77, 0.77);

      TrackCutsDeuteron->SetPtRange(0.6, 2.5);
      AntiTrackCutsDeuteron->SetPtRange(0.6, 2.5);
      TrackCutsDeuteron->SetEtaRange(-0.77, 0.77);
      AntiTrackCutsDeuteron->SetEtaRange(-0.77, 0.77);

      config->SetDeltaEtaMax(0.019);
      config->SetDeltaPhiMax(0.019);
    } else if (suffix == "3") {
      TrackCuts->SetEtaRange(-0.77, 0.77);
      AntiTrackCuts->SetEtaRange(-0.77, 0.77);
      TrackCuts->SetNClsTPC(70);
      AntiTrackCuts->SetNClsTPC(70);

      TrackCutsDeuteron->SetEtaRange(-0.77, 0.77);
      AntiTrackCutsDeuteron->SetEtaRange(-0.77, 0.77);
      TrackCutsDeuteron->SetNClsTPC(70);
      AntiTrackCutsDeuteron->SetNClsTPC(70);
    } else if (suffix == "4") {

      TrackCuts->SetPtRange(0.4, 2.5);
      AntiTrackCuts->SetPtRange(0.4, 2.5);
      TrackCuts->SetEtaRange(-0.77, 0.77);
      AntiTrackCuts->SetEtaRange(-0.77, 0.77);
      TrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);
      AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);

      TrackCutsDeuteron->SetPtRange(0.4, 2.0);
      AntiTrackCutsDeuteron->SetPtRange(0.4, 2.0);
      TrackCutsDeuteron->SetEtaRange(-0.77, 0.77);
      AntiTrackCutsDeuteron->SetEtaRange(-0.77, 0.77);
      TrackCutsDeuteron->SetPID(AliPID::kDeuteron, 1.4, 3.3);
      AntiTrackCutsDeuteron->SetPID(AliPID::kDeuteron, 1.4, 3.3);

      config->SetDeltaEtaMax(0.019);
      config->SetDeltaPhiMax(0.019);

    } else if (suffix == "5") {
      TrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);
      AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);
      TrackCuts->SetNClsTPC(90);
      AntiTrackCuts->SetNClsTPC(90);

      TrackCutsDeuteron->SetPID(AliPID::kDeuteron, 1.4, 3.3);
      AntiTrackCutsDeuteron->SetPID(AliPID::kDeuteron, 1.4, 3.3);
      TrackCutsDeuteron->SetNClsTPC(90);
      AntiTrackCutsDeuteron->SetNClsTPC(90);

    } else if (suffix == "6") {
      TrackCuts->SetEtaRange(-0.83, 0.83);
      AntiTrackCuts->SetEtaRange(-0.83, 0.83);
      TrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);
      AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);

      TrackCutsDeuteron->SetEtaRange(-0.83, 0.83);
      AntiTrackCutsDeuteron->SetEtaRange(-0.83, 0.83);
      TrackCutsDeuteron->SetPID(AliPID::kDeuteron, 1.4, 3.3);
      AntiTrackCutsDeuteron->SetPID(AliPID::kDeuteron, 1.4, 3.3);

      config->SetDeltaEtaMax(0.019);
      config->SetDeltaPhiMax(0.019);

    } else if (suffix == "7") {
      TrackCuts->SetPtRange(0.4, 2.5);
      AntiTrackCuts->SetPtRange(0.4, 2.5);
      TrackCuts->SetEtaRange(-0.77, 0.77);
      AntiTrackCuts->SetEtaRange(-0.77, 0.77);
      TrackCuts->SetNClsTPC(90);
      AntiTrackCuts->SetNClsTPC(90);

      TrackCutsDeuteron->SetPtRange(0.4, 2.0);
      AntiTrackCutsDeuteron->SetPtRange(0.4, 2.0);
      TrackCutsDeuteron->SetEtaRange(-0.77, 0.77);
      AntiTrackCutsDeuteron->SetEtaRange(-0.77, 0.77);
      TrackCutsDeuteron->SetNClsTPC(90);
      AntiTrackCutsDeuteron->SetNClsTPC(90);

      config->SetDeltaEtaMax(0.019);
      config->SetDeltaPhiMax(0.019);

    } else if (suffix == "8") {
      TrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);
      AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);
      TrackCuts->SetNClsTPC(90);
      AntiTrackCuts->SetNClsTPC(90);

      TrackCutsDeuteron->SetPID(AliPID::kDeuteron, 1.4, 2);
      AntiTrackCutsDeuteron->SetPID(AliPID::kDeuteron, 1.4, 2);
      TrackCutsDeuteron->SetNClsTPC(90);
      AntiTrackCutsDeuteron->SetNClsTPC(90);
    } else if (suffix == "9") {

      TrackCuts->SetPtRange(0.4, 2.5);
      AntiTrackCuts->SetPtRange(0.4, 2.5);
      TrackCuts->SetEtaRange(-0.83, 0.83);
      AntiTrackCuts->SetEtaRange(-0.83, 0.83);

      TrackCutsDeuteron->SetPtRange(0.4, 2.0);
      AntiTrackCutsDeuteron->SetPtRange(0.4, 2.0);
      TrackCutsDeuteron->SetEtaRange(-0.83, 0.83);
      AntiTrackCutsDeuteron->SetEtaRange(-0.83, 0.83);

    } else if (suffix == "10") {

      TrackCuts->SetEtaRange(-0.77, 0.77);
      AntiTrackCuts->SetEtaRange(-0.77, 0.77);
      TrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);
      AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);

      TrackCutsDeuteron->SetEtaRange(-0.77, 0.77);
      AntiTrackCutsDeuteron->SetEtaRange(-0.77, 0.77);
      TrackCutsDeuteron->SetPID(AliPID::kDeuteron, 1.4, 2.7);
      AntiTrackCutsDeuteron->SetPID(AliPID::kDeuteron, 1.4, 2.7);

    } else if (suffix == "11") {

      TrackCuts->SetPtRange(0.4, 2.5);
      AntiTrackCuts->SetPtRange(0.4, 2.5);
      TrackCuts->SetEtaRange(-0.83, 0.83);
      AntiTrackCuts->SetEtaRange(-0.83, 0.83);
      TrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);
      AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);

      TrackCutsDeuteron->SetPtRange(0.4, 2.0);
      AntiTrackCutsDeuteron->SetPtRange(0.4, 2.0);
      TrackCutsDeuteron->SetEtaRange(-0.83, 0.83);
      AntiTrackCutsDeuteron->SetEtaRange(-0.83, 0.83);
      TrackCutsDeuteron->SetPID(AliPID::kDeuteron, 1.4, 2.7);
      AntiTrackCutsDeuteron->SetPID(AliPID::kDeuteron, 1.4, 2.7);

      config->SetDeltaEtaMax(0.019);
      config->SetDeltaPhiMax(0.019);

    } else if (suffix == "12") {
      TrackCuts->SetEtaRange(-0.83, 0.83);
      AntiTrackCuts->SetEtaRange(-0.83, 0.83);
      TrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);
      AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);
      TrackCuts->SetNClsTPC(90);
      AntiTrackCuts->SetNClsTPC(90);

      TrackCutsDeuteron->SetEtaRange(-0.83, 0.83);
      AntiTrackCutsDeuteron->SetEtaRange(-0.83, 0.83);
      TrackCutsDeuteron->SetPID(AliPID::kDeuteron, 1.4, 2.7);
      AntiTrackCutsDeuteron->SetPID(AliPID::kDeuteron, 1.4, 2.7);
      TrackCutsDeuteron->SetNClsTPC(90);
      AntiTrackCutsDeuteron->SetNClsTPC(90);

    } else if (suffix == "13") {
      TrackCuts->SetPtRange(0.4, 2.5);
      AntiTrackCuts->SetPtRange(0.4, 2.5);
      TrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);
      AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);

      TrackCutsDeuteron->SetPtRange(0.4, 2.0);
      AntiTrackCutsDeuteron->SetPtRange(0.4, 2.0);
      TrackCutsDeuteron->SetPID(AliPID::kDeuteron, 1.4, 2.7);
      AntiTrackCutsDeuteron->SetPID(AliPID::kDeuteron, 1.4, 2.7);

      config->SetDeltaEtaMax(0.019);
      config->SetDeltaPhiMax(0.019);

    } else if (suffix == "14") {
      TrackCuts->SetPtRange(0.6, 2.5);
      AntiTrackCuts->SetPtRange(0.6, 2.5);
      TrackCuts->SetEtaRange(-0.83, 0.83);
      AntiTrackCuts->SetEtaRange(-0.83, 0.83);
      TrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);
      AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);
      TrackCuts->SetNClsTPC(90);
      AntiTrackCuts->SetNClsTPC(90);

      TrackCutsDeuteron->SetPtRange(0.6, 2.0);
      AntiTrackCutsDeuteron->SetPtRange(0.6, 2.0);
      TrackCutsDeuteron->SetEtaRange(-0.83, 0.83);
      AntiTrackCutsDeuteron->SetEtaRange(-0.83, 0.83);
      TrackCutsDeuteron->SetPID(AliPID::kDeuteron, 1.4, 2.7);
      AntiTrackCutsDeuteron->SetPID(AliPID::kDeuteron, 1.4, 2.7);
      TrackCutsDeuteron->SetNClsTPC(90);
      AntiTrackCutsDeuteron->SetNClsTPC(90);

    } else if (suffix == "15") {
      TrackCuts->SetPtRange(0.4, 2.5);
      AntiTrackCuts->SetPtRange(0.4, 2.5);
      TrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);
      AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);
      TrackCuts->SetNClsTPC(70);
      AntiTrackCuts->SetNClsTPC(70);

      TrackCutsDeuteron->SetPtRange(0.4, 2.0);
      AntiTrackCutsDeuteron->SetPtRange(0.4, 2.0);
      TrackCutsDeuteron->SetPID(AliPID::kDeuteron, 1.4, 3.3);
      AntiTrackCutsDeuteron->SetPID(AliPID::kDeuteron, 1.4, 3.3);
      TrackCutsDeuteron->SetNClsTPC(90);
      AntiTrackCutsDeuteron->SetNClsTPC(90);

      config->SetDeltaEtaMax(0.019);
      config->SetDeltaPhiMax(0.019);

    } else if (suffix == "16") {
      TrackCuts->SetPtRange(0.4, 2.5);
      AntiTrackCuts->SetPtRange(0.4,2.5);
      TrackCuts->SetEtaRange(-0.83, 0.83);
      AntiTrackCuts->SetEtaRange(-0.83, 0.83);
      TrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);
      AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);

      TrackCutsDeuteron->SetPtRange(0.4, 2.0);
      AntiTrackCutsDeuteron->SetPtRange(0.4, 2.0);
      TrackCutsDeuteron->SetEtaRange(-0.83, 0.83);
      AntiTrackCutsDeuteron->SetEtaRange(-0.83, 0.83);
      TrackCutsDeuteron->SetPID(AliPID::kDeuteron, 1.4, 2.7);
      AntiTrackCutsDeuteron->SetPID(AliPID::kDeuteron, 1.4, 2.7);

      config->SetDeltaEtaMax(0.019);
      config->SetDeltaPhiMax(0.019);

    } else if (suffix == "17") {
      TrackCuts->SetPtRange(0.4, 2.5);
      AntiTrackCuts->SetPtRange(0.4, 2.5);
      TrackCuts->SetEtaRange(-0.83, 0.83);
      AntiTrackCuts->SetEtaRange(-0.83, 0.83);
      TrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);
      AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);
      TrackCuts->SetNClsTPC(90);
      AntiTrackCuts->SetNClsTPC(90);

      TrackCutsDeuteron->SetPtRange(0.4, 2.0);
      AntiTrackCutsDeuteron->SetPtRange(0.4, 2.0);
      TrackCutsDeuteron->SetEtaRange(-0.83, 0.83);
      AntiTrackCutsDeuteron->SetEtaRange(-0.83, 0.83);
      TrackCutsDeuteron->SetPID(AliPID::kDeuteron, 1.4, 2.7);
      AntiTrackCutsDeuteron->SetPID(AliPID::kDeuteron, 1.4, 2.7);
      TrackCutsDeuteron->SetNClsTPC(90);
      AntiTrackCutsDeuteron->SetNClsTPC(90);

    } else if (suffix == "18") {
      TrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);
      AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);
      TrackCuts->SetNClsTPC(70);
      AntiTrackCuts->SetNClsTPC(70);

      TrackCutsDeuteron->SetPID(AliPID::kDeuteron, 1.4, 2.7);
      AntiTrackCutsDeuteron->SetPID(AliPID::kDeuteron, 1.4, 2.7);
      TrackCutsDeuteron->SetNClsTPC(70);
      AntiTrackCutsDeuteron->SetNClsTPC(70);

      config->SetDeltaEtaMax(0.019);
      config->SetDeltaPhiMax(0.019);

    } else if (suffix == "19") {
      TrackCuts->SetPtRange(0.4, 2.5);
      AntiTrackCuts->SetPtRange(0.4, 2.5);
      TrackCuts->SetNClsTPC(90);
      AntiTrackCuts->SetNClsTPC(90);

      TrackCutsDeuteron->SetPtRange(0.4, 2.0);
      AntiTrackCutsDeuteron->SetPtRange(0.4, 2.0);
      TrackCutsDeuteron->SetNClsTPC(90);
      AntiTrackCutsDeuteron->SetNClsTPC(90);
      config->SetDeltaEtaMax(0.019);
      config->SetDeltaPhiMax(0.019);

    } else if (suffix == "20") {

      TrackCuts->SetEtaRange(-0.83, 0.83);
      AntiTrackCuts->SetEtaRange(-0.83, 0.83);
      TrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);
      AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);
      TrackCuts->SetNClsTPC(70);
      AntiTrackCuts->SetNClsTPC(70);

      TrackCutsDeuteron->SetEtaRange(-0.83, 0.83);
      AntiTrackCutsDeuteron->SetEtaRange(-0.83, 0.83);
      TrackCutsDeuteron->SetPID(AliPID::kDeuteron, 1.4, 3.3);
      AntiTrackCutsDeuteron->SetPID(AliPID::kDeuteron, 1.4, 3.3);
      TrackCutsDeuteron->SetNClsTPC(70);
      AntiTrackCutsDeuteron->SetNClsTPC(70);

      config->SetDeltaEtaMax(0.019);
      config->SetDeltaPhiMax(0.019);

    } else if (suffix == "21") {
      TrackCuts->SetPtRange(0.4, 2.5);
      AntiTrackCuts->SetPtRange(0.4, 2.5);
      TrackCuts->SetEtaRange(-0.77, 0.77);
      AntiTrackCuts->SetEtaRange(-0.77, 0.77);
      TrackCuts->SetNClsTPC(70);
      AntiTrackCuts->SetNClsTPC(70);

      TrackCutsDeuteron->SetPtRange(0.4, 2.0);
      AntiTrackCutsDeuteron->SetPtRange(0.4, 2.0);
      TrackCutsDeuteron->SetEtaRange(-0.77, 0.77);
      AntiTrackCutsDeuteron->SetEtaRange(-0.77, 0.77);
      TrackCutsDeuteron->SetNClsTPC(70);
      AntiTrackCutsDeuteron->SetNClsTPC(70);

      config->SetDeltaEtaMax(0.019);
      config->SetDeltaPhiMax(0.019);

    } else if (suffix == "22") {
      TrackCuts->SetEtaRange(-0.77, 0.77);
      AntiTrackCuts->SetEtaRange(-0.77, 0.77);
      TrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);
      AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);

      TrackCutsDeuteron->SetEtaRange(-0.77, 0.77);
      AntiTrackCutsDeuteron->SetEtaRange(-0.77, 0.77);
      TrackCutsDeuteron->SetPID(AliPID::kDeuteron, 1.4, 2.7);
      AntiTrackCutsDeuteron->SetPID(AliPID::kDeuteron, 1.4, 2.7);

    } else if (suffix == "23") {
      TrackCuts->SetEtaRange(-0.77, 0.77);
      AntiTrackCuts->SetEtaRange(-0.77, 0.77);
      TrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);
      AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);
      TrackCuts->SetNClsTPC(90);
      AntiTrackCuts->SetNClsTPC(90);

      TrackCutsDeuteron->SetEtaRange(-0.77, 0.77);
      AntiTrackCutsDeuteron->SetEtaRange(-0.77, 0.77);
      TrackCutsDeuteron->SetPID(AliPID::kDeuteron, 1.4, 3.3);
      AntiTrackCutsDeuteron->SetPID(AliPID::kDeuteron, 1.4, 3.3);
      TrackCutsDeuteron->SetNClsTPC(90);
      AntiTrackCutsDeuteron->SetNClsTPC(90);

      config->SetDeltaEtaMax(0.019);
      config->SetDeltaPhiMax(0.019);

    } else if (suffix == "24") {
      TrackCuts->SetPtRange(0.4, 2.5);
      AntiTrackCuts->SetPtRange(0.4, 2.5);
      TrackCuts->SetEtaRange(-0.83, 0.83);
      AntiTrackCuts->SetEtaRange(-0.83, 0.83);

      TrackCutsDeuteron->SetPtRange(0.4, 2.0);
      AntiTrackCutsDeuteron->SetPtRange(0.4, 2.0);
      TrackCutsDeuteron->SetEtaRange(-0.83, 0.83);
      AntiTrackCutsDeuteron->SetEtaRange(-0.83, 0.83);

      config->SetDeltaEtaMax(0.019);
      config->SetDeltaPhiMax(0.019);

    } else if (suffix == "25") {
      TrackCuts->SetPtRange(0.6, 2.5);
      AntiTrackCuts->SetPtRange(0.6, 2.5);
      TrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);
      AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);

      TrackCutsDeuteron->SetPtRange(0.6, 2.0);
      AntiTrackCutsDeuteron->SetPtRange(0.6, 2.0);
      TrackCutsDeuteron->SetPID(AliPID::kDeuteron, 1.4, 2.7);
      AntiTrackCutsDeuteron->SetPID(AliPID::kDeuteron, 1.4, 2.7);

      config->SetDeltaEtaMax(0.019);
      config->SetDeltaPhiMax(0.019);

    } else if (suffix == "26") {

      TrackCuts->SetEtaRange(-0.77, 0.77);
      AntiTrackCuts->SetEtaRange(-0.77, 0.77);
      TrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);
      AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);
      TrackCuts->SetNClsTPC(90);
      AntiTrackCuts->SetNClsTPC(90);

      TrackCutsDeuteron->SetEtaRange(-0.77, 0.77);
      AntiTrackCutsDeuteron->SetEtaRange(-0.77, 0.77);
      TrackCutsDeuteron->SetPID(AliPID::kDeuteron, 1.4, 3.3);
      AntiTrackCutsDeuteron->SetPID(AliPID::kDeuteron, 1.4, 3.3);
      TrackCutsDeuteron->SetNClsTPC(90);
      AntiTrackCutsDeuteron->SetNClsTPC(90);

    } else if (suffix == "27") {
      TrackCuts->SetNClsTPC(90);
      AntiTrackCuts->SetNClsTPC(90);
      TrackCutsDeuteron->SetNClsTPC(90);
      AntiTrackCutsDeuteron->SetNClsTPC(90);

    } else if (suffix == "28") {
      TrackCuts->SetPtRange(0.4, 2.5);
      AntiTrackCuts->SetPtRange(0.4, 2.5);
      TrackCuts->SetEtaRange(-0.77, 0.77);
      AntiTrackCuts->SetEtaRange(-0.77, 0.77);
      TrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);
      AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);

      TrackCutsDeuteron->SetPtRange(0.4, 2.0);
      AntiTrackCutsDeuteron->SetPtRange(0.4, 2.0);
      TrackCutsDeuteron->SetEtaRange(-0.77, 0.77);
      AntiTrackCutsDeuteron->SetEtaRange(-0.77, 0.77);
      TrackCutsDeuteron->SetPID(AliPID::kDeuteron, 1.4, 2.7);
      AntiTrackCutsDeuteron->SetPID(AliPID::kDeuteron, 1.4, 2.7);

      config->SetDeltaEtaMax(0.019);
      config->SetDeltaPhiMax(0.019);

    } else if (suffix == "29") {
      TrackCuts->SetEtaRange(-0.83, 0.83);
      AntiTrackCuts->SetEtaRange(-0.83, 0.83);

      TrackCutsDeuteron->SetEtaRange(-0.83, 0.83);
      AntiTrackCutsDeuteron->SetEtaRange(-0.83, 0.83);

      config->SetDeltaEtaMax(0.019);
      config->SetDeltaPhiMax(0.019);

    } else if (suffix == "30") {
      TrackCuts->SetPtRange(0.4, 2.5);
      AntiTrackCuts->SetPtRange(0.4, 2.5);
      TrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);
      AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);
      TrackCuts->SetNClsTPC(70);
      AntiTrackCuts->SetNClsTPC(70);

      TrackCutsDeuteron->SetPtRange(0.4, 2.0);
      AntiTrackCutsDeuteron->SetPtRange(0.4, 2.0);
      TrackCutsDeuteron->SetPID(AliPID::kDeuteron, 1.4, 3.3);
      AntiTrackCutsDeuteron->SetPID(AliPID::kDeuteron, 1.4, 3.3);
      TrackCutsDeuteron->SetNClsTPC(70);
      AntiTrackCutsDeuteron->SetNClsTPC(70);

      config->SetDeltaEtaMax(0.019);
      config->SetDeltaPhiMax(0.019);

    } else if (suffix == "31") {

      TrackCuts->SetEtaRange(-0.77, 0.77);
      AntiTrackCuts->SetEtaRange(-0.77, 0.77);
      TrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);
      AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);
      TrackCuts->SetNClsTPC(90);
      AntiTrackCuts->SetNClsTPC(90);

      TrackCutsDeuteron->SetEtaRange(-0.77, 0.77);
      AntiTrackCutsDeuteron->SetEtaRange(-0.77, 0.77);
      TrackCutsDeuteron->SetPID(AliPID::kDeuteron, 1.4, 2.7);
      AntiTrackCutsDeuteron->SetPID(AliPID::kDeuteron, 1.4, 2.7);
      TrackCutsDeuteron->SetNClsTPC(90);
      AntiTrackCutsDeuteron->SetNClsTPC(90);

    } else if (suffix == "32") {
      TrackCuts->SetEtaRange(-0.77, 0.77);
      AntiTrackCuts->SetEtaRange(-0.77, 0.77);

      TrackCutsDeuteron->SetEtaRange(-0.77, 0.77);
      AntiTrackCutsDeuteron->SetEtaRange(-0.77, 0.77);

      config->SetDeltaEtaMax(0.019);
      config->SetDeltaPhiMax(0.019);

    } else if (suffix == "33") {
      TrackCuts->SetEtaRange(-0.77, 0.77);
      AntiTrackCuts->SetEtaRange(-0.77, 0.77);
      TrackCuts->SetNClsTPC(70);
      AntiTrackCuts->SetNClsTPC(70);

      TrackCutsDeuteron->SetEtaRange(-0.77, 0.77);
      AntiTrackCutsDeuteron->SetEtaRange(-0.77, 0.77);
      TrackCutsDeuteron->SetNClsTPC(70);
      AntiTrackCutsDeuteron->SetNClsTPC(70);

      config->SetDeltaEtaMax(0.019);
      config->SetDeltaPhiMax(0.019);

    } else if (suffix == "34") {
      TrackCuts->SetPtRange(0.4, 2.5);
      AntiTrackCuts->SetPtRange(0.4, 2.5);
      TrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);
      AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);
      TrackCuts->SetNClsTPC(70);
      AntiTrackCuts->SetNClsTPC(70);

      TrackCutsDeuteron->SetPtRange(0.4, 2.0);
      AntiTrackCutsDeuteron->SetPtRange(0.4, 2.0);
      TrackCutsDeuteron->SetPID(AliPID::kDeuteron, 1.4, 2.7);
      AntiTrackCutsDeuteron->SetPID(AliPID::kDeuteron, 1.4, 2.7);
      TrackCutsDeuteron->SetNClsTPC(70);
      AntiTrackCutsDeuteron->SetNClsTPC(70);

      config->SetDeltaEtaMax(0.019);
      config->SetDeltaPhiMax(0.019);

    } else if (suffix == "35") {
      TrackCuts->SetPtRange(0.6, 2.5);
      AntiTrackCuts->SetPtRange(0.6, 2.5);
      TrackCuts->SetEtaRange(-0.83, 0.83);
      AntiTrackCuts->SetEtaRange(-0.83, 0.83);
      TrackCuts->SetNClsTPC(70);
      AntiTrackCuts->SetNClsTPC(70);

      TrackCutsDeuteron->SetPtRange(0.6, 2.0);
      AntiTrackCutsDeuteron->SetPtRange(0.6, 2.0);
      TrackCutsDeuteron->SetEtaRange(-0.83, 0.83);
      AntiTrackCutsDeuteron->SetEtaRange(-0.83, 0.83);
      TrackCutsDeuteron->SetNClsTPC(70);
      AntiTrackCutsDeuteron->SetNClsTPC(70);

      config->SetDeltaEtaMax(0.019);
      config->SetDeltaPhiMax(0.019);

    } else if (suffix == "36") {
      TrackCuts->SetEtaRange(-0.83, 0.83);
      AntiTrackCuts->SetEtaRange(-0.83, 0.83);
      TrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);
      AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);
      TrackCuts->SetNClsTPC(90);
      AntiTrackCuts->SetNClsTPC(90);

      TrackCutsDeuteron->SetEtaRange(-0.83, 0.83);
      AntiTrackCutsDeuteron->SetEtaRange(-0.83, 0.83);
      TrackCutsDeuteron->SetPID(AliPID::kDeuteron, 1.4, 3.3);
      AntiTrackCutsDeuteron->SetPID(AliPID::kDeuteron, 1.4, 3.3);
      TrackCutsDeuteron->SetNClsTPC(90);
      AntiTrackCutsDeuteron->SetNClsTPC(90);

      config->SetDeltaEtaMax(0.019);
      config->SetDeltaPhiMax(0.019);

    } else if (suffix == "37") {
      TrackCuts->SetEtaRange(-0.77, 0.77);
      AntiTrackCuts->SetEtaRange(-0.77, 0.77);
      TrackCuts->SetNClsTPC(90);
      AntiTrackCuts->SetNClsTPC(90);

      TrackCutsDeuteron->SetEtaRange(-0.77, 0.77);
      AntiTrackCutsDeuteron->SetEtaRange(-0.77, 0.77);
      TrackCutsDeuteron->SetNClsTPC(90);
      AntiTrackCutsDeuteron->SetNClsTPC(90);

      config->SetDeltaEtaMax(0.019);
      config->SetDeltaPhiMax(0.019);

    } else if (suffix == "38") {

      TrackCuts->SetPtRange(0.4, 2.5);
      AntiTrackCuts->SetPtRange(0.4, 2.5);
      TrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);
      AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);
      TrackCuts->SetNClsTPC(90);
      AntiTrackCuts->SetNClsTPC(90);

      TrackCutsDeuteron->SetPtRange(0.4, 2.0);
      AntiTrackCutsDeuteron->SetPtRange(0.4, 2.0);
      TrackCutsDeuteron->SetPID(AliPID::kDeuteron, 1.4, 3.3);
      AntiTrackCutsDeuteron->SetPID(AliPID::kDeuteron, 1.4, 3.3);
      TrackCutsDeuteron->SetNClsTPC(90);
      AntiTrackCutsDeuteron->SetNClsTPC(90);

      config->SetDeltaEtaMax(0.019);
      config->SetDeltaPhiMax(0.019);

    } else if (suffix == "39") {
      TrackCuts->SetPtRange(0.6, 2.5);
      AntiTrackCuts->SetPtRange(0.6, 2.5);
      TrackCuts->SetNClsTPC(70);
      AntiTrackCuts->SetNClsTPC(70);

      TrackCutsDeuteron->SetPtRange(0.6, 2.0);
      AntiTrackCutsDeuteron->SetPtRange(0.6, 2.0);
      TrackCutsDeuteron->SetNClsTPC(90);
      AntiTrackCutsDeuteron->SetNClsTPC(90);

    } else if (suffix == "40") {
      TrackCuts->SetPtRange(0.4, 2.5);
      AntiTrackCuts->SetPtRange(0.4, 2.5);
      TrackCuts->SetEtaRange(-0.77, 0.77);
      AntiTrackCuts->SetEtaRange(-0.77, 0.77);
      TrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);
      AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);

      TrackCutsDeuteron->SetPtRange(0.4, 2.0);
      AntiTrackCutsDeuteron->SetPtRange(0.4, 2.0);
      TrackCutsDeuteron->SetEtaRange(-0.77, 0.77);
      AntiTrackCutsDeuteron->SetEtaRange(-0.77, 0.77);
      TrackCutsDeuteron->SetPID(AliPID::kDeuteron, 1.4, 3.3);
      AntiTrackCutsDeuteron->SetPID(AliPID::kDeuteron, 1.4, 3.3);

      config->SetDeltaEtaMax(0.019);
      config->SetDeltaPhiMax(0.019);

    } else if (suffix == "41") {
      TrackCuts->SetPtRange(0.6, 2.5);
      AntiTrackCuts->SetPtRange(0.6, 2.5);
      TrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);
      AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);

      TrackCutsDeuteron->SetPtRange(0.6, 2.0);
      AntiTrackCutsDeuteron->SetPtRange(0.6, 2.0);
      TrackCutsDeuteron->SetPID(AliPID::kDeuteron, 1.4, 3.3);
      AntiTrackCutsDeuteron->SetPID(AliPID::kDeuteron, 1.4, 3.3);

    } else if (suffix == "42") {
      TrackCuts->SetPtRange(0.6, 2.5);
      AntiTrackCuts->SetPtRange(0.6, 2.5);
      TrackCuts->SetEtaRange(-0.83, 0.83);
      AntiTrackCuts->SetEtaRange(-0.83, 0.83);
      TrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);
      AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);

      TrackCutsDeuteron->SetPtRange(0.6, 2.0);
      AntiTrackCutsDeuteron->SetPtRange(0.6, 2.0);
      TrackCutsDeuteron->SetEtaRange(-0.83, 0.83);
      AntiTrackCutsDeuteron->SetEtaRange(-0.83, 0.83);
      TrackCutsDeuteron->SetPID(AliPID::kDeuteron, 1.4, 3.3);
      AntiTrackCutsDeuteron->SetPID(AliPID::kDeuteron, 1.4, 3.3);

      config->SetDeltaEtaMax(0.019);
      config->SetDeltaPhiMax(0.019);

    } else if (suffix == "43") {
      TrackCuts->SetEtaRange(-0.83, 0.83);
      AntiTrackCuts->SetEtaRange(-0.83, 0.83);
      TrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);
      AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);

      TrackCutsDeuteron->SetEtaRange(-0.83, 0.83);
      AntiTrackCutsDeuteron->SetEtaRange(-0.83, 0.83);
      TrackCutsDeuteron->SetPID(AliPID::kDeuteron, 1.4, 3.3);
      AntiTrackCutsDeuteron->SetPID(AliPID::kDeuteron, 1.4, 3.3);

      config->SetDeltaEtaMax(0.019);
      config->SetDeltaPhiMax(0.019);

    } else if (suffix == "44") {

      TrackCuts->SetPtRange(0.4, 2.5);
      AntiTrackCuts->SetPtRange(0.4, 2.5);
      TrackCuts->SetEtaRange(-0.77, 0.77);
      AntiTrackCuts->SetEtaRange(-0.77, 0.77);
      TrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);
      AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);

      TrackCutsDeuteron->SetPtRange(0.4, 2.0);
      AntiTrackCutsDeuteron->SetPtRange(0.4, 2.0);
      TrackCutsDeuteron->SetEtaRange(-0.77, 0.77);
      AntiTrackCutsDeuteron->SetEtaRange(-0.77, 0.77);
      TrackCutsDeuteron->SetPID(AliPID::kDeuteron, 1.4, 3.3);
      AntiTrackCutsDeuteron->SetPID(AliPID::kDeuteron, 1.4, 3.3);

      config->SetDeltaEtaMax(0.015);
      config->SetDeltaPhiMax(0.015);

    } else if (suffix == "45") {
      TrackCuts->SetEtaRange(-0.83, 0.83);
      AntiTrackCuts->SetEtaRange(-0.83, 0.83);
      TrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);
      AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);

      TrackCutsDeuteron->SetEtaRange(-0.83, 0.83);
      AntiTrackCutsDeuteron->SetEtaRange(-0.83, 0.83);
      TrackCutsDeuteron->SetPID(AliPID::kDeuteron, 1.4, 3.3);
      AntiTrackCutsDeuteron->SetPID(AliPID::kDeuteron, 1.4, 3.3);

      config->SetDeltaEtaMax(0.015);
      config->SetDeltaPhiMax(0.015);

    } else if (suffix == "46") {
      TrackCuts->SetPtRange(0.4, 2.5);
      AntiTrackCuts->SetPtRange(0.4, 2.5);
      TrackCuts->SetEtaRange(-0.77, 0.77);
      AntiTrackCuts->SetEtaRange(-0.77, 0.77);
      TrackCuts->SetNClsTPC(90);
      AntiTrackCuts->SetNClsTPC(90);

      TrackCutsDeuteron->SetPtRange(0.4, 2.0);
      AntiTrackCutsDeuteron->SetPtRange(0.4, 2.0);
      TrackCutsDeuteron->SetEtaRange(-0.77, 0.77);
      AntiTrackCutsDeuteron->SetEtaRange(-0.77, 0.77);
      TrackCutsDeuteron->SetNClsTPC(90);
      AntiTrackCutsDeuteron->SetNClsTPC(90);

      config->SetDeltaEtaMax(0.015);
      config->SetDeltaPhiMax(0.015);

    } else if (suffix == "47") {

      TrackCuts->SetPtRange(0.4, 2.5);
      AntiTrackCuts->SetPtRange(0.4, 2.5);
      TrackCuts->SetEtaRange(-0.83, 0.83);
      AntiTrackCuts->SetEtaRange(-0.83, 0.83);
      TrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);
      AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);

      TrackCutsDeuteron->SetPtRange(0.4, 2.0);
      AntiTrackCutsDeuteron->SetPtRange(0.4, 2.0);
      TrackCutsDeuteron->SetEtaRange(-0.83, 0.83);
      AntiTrackCutsDeuteron->SetEtaRange(-0.83, 0.83);
      TrackCutsDeuteron->SetPID(AliPID::kDeuteron, 1.4, 2.7);
      AntiTrackCutsDeuteron->SetPID(AliPID::kDeuteron, 1.4, 2.7);

      config->SetDeltaEtaMax(0.015);
      config->SetDeltaPhiMax(0.015);

    } else if (suffix == "48") {
      TrackCuts->SetPtRange(0.4, 2.5);
      AntiTrackCuts->SetPtRange(0.4, 2.5);
      TrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);
      AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);

      TrackCutsDeuteron->SetPtRange(0.4, 2.0);
      AntiTrackCutsDeuteron->SetPtRange(0.4, 2.0);
      TrackCutsDeuteron->SetPID(AliPID::kDeuteron, 1.4, 2.7);
      AntiTrackCutsDeuteron->SetPID(AliPID::kDeuteron, 1.4, 2.7);

      config->SetDeltaEtaMax(0.015);
      config->SetDeltaPhiMax(0.015);

    } else if (suffix == "49") {
      TrackCuts->SetPtRange(0.4, 2.5);
      AntiTrackCuts->SetPtRange(0.4, 2.5);
      TrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);
      AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);
      TrackCuts->SetNClsTPC(70);
      AntiTrackCuts->SetNClsTPC(70);

      TrackCutsDeuteron->SetPtRange(0.4, 2.0);
      AntiTrackCutsDeuteron->SetPtRange(0.4, 2.0);
      TrackCutsDeuteron->SetPID(AliPID::kDeuteron, 1.4, 3.3);
      AntiTrackCutsDeuteron->SetPID(AliPID::kDeuteron, 1.4, 3.3);
      TrackCutsDeuteron->SetNClsTPC(90);
      AntiTrackCutsDeuteron->SetNClsTPC(90);

      config->SetDeltaEtaMax(0.015);
      config->SetDeltaPhiMax(0.015);

    } else if (suffix == "50") {
      TrackCuts->SetPtRange(0.4, 2.5);
      AntiTrackCuts->SetPtRange(0.4, 2.5);
      TrackCuts->SetEtaRange(-0.83, 0.83);
      AntiTrackCuts->SetEtaRange(-0.83, 0.83);
      TrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);
      AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);

      TrackCutsDeuteron->SetPtRange(0.4, 2.0);
      AntiTrackCutsDeuteron->SetPtRange(0.4, 2.0);
      TrackCutsDeuteron->SetEtaRange(-0.83, 0.83);
      AntiTrackCutsDeuteron->SetEtaRange(-0.83, 0.83);
      TrackCutsDeuteron->SetPID(AliPID::kDeuteron, 1.4, 2.7);
      AntiTrackCutsDeuteron->SetPID(AliPID::kDeuteron, 1.4, 2.7);

      config->SetDeltaEtaMax(0.015);
      config->SetDeltaPhiMax(0.015);

    } else if (suffix == "51") {
      TrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);
      AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);
      TrackCuts->SetNClsTPC(70);
      AntiTrackCuts->SetNClsTPC(70);

      TrackCutsDeuteron->SetPID(AliPID::kDeuteron, 1.4, 2.7);
      AntiTrackCutsDeuteron->SetPID(AliPID::kDeuteron, 1.4, 2.7);
      TrackCutsDeuteron->SetNClsTPC(70);
      AntiTrackCutsDeuteron->SetNClsTPC(70);

      config->SetDeltaEtaMax(0.015);
      config->SetDeltaPhiMax(0.015);

    } else if (suffix == "52") {
      TrackCuts->SetPtRange(0.4, 2.5);
      AntiTrackCuts->SetPtRange(0.4, 2.5);
      TrackCuts->SetNClsTPC(90);
      AntiTrackCuts->SetNClsTPC(90);

      TrackCutsDeuteron->SetPtRange(0.4, 2.0);
      AntiTrackCutsDeuteron->SetPtRange(0.4, 2.0);
      TrackCutsDeuteron->SetNClsTPC(90);
      AntiTrackCutsDeuteron->SetNClsTPC(90);
      config->SetDeltaEtaMax(0.015);
      config->SetDeltaPhiMax(0.015);

    } else if (suffix == "53") {

      TrackCuts->SetEtaRange(-0.83, 0.83);
      AntiTrackCuts->SetEtaRange(-0.83, 0.83);
      TrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);
      AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);
      TrackCuts->SetNClsTPC(70);
      AntiTrackCuts->SetNClsTPC(70);

      TrackCutsDeuteron->SetEtaRange(-0.83, 0.83);
      AntiTrackCutsDeuteron->SetEtaRange(-0.83, 0.83);
      TrackCutsDeuteron->SetPID(AliPID::kDeuteron, 1.4, 3.3);
      AntiTrackCutsDeuteron->SetPID(AliPID::kDeuteron, 1.4, 3.3);
      TrackCutsDeuteron->SetNClsTPC(70);
      AntiTrackCutsDeuteron->SetNClsTPC(70);

      config->SetDeltaEtaMax(0.015);
      config->SetDeltaPhiMax(0.015);

    } else if (suffix == "54") {
      TrackCuts->SetPtRange(0.4, 2.5);
      AntiTrackCuts->SetPtRange(0.4, 2.5);
      TrackCuts->SetEtaRange(-0.77, 0.77);
      AntiTrackCuts->SetEtaRange(-0.77, 0.77);
      TrackCuts->SetNClsTPC(70);
      AntiTrackCuts->SetNClsTPC(70);

      TrackCutsDeuteron->SetPtRange(0.4, 2.0);
      AntiTrackCutsDeuteron->SetPtRange(0.4, 2.0);
      TrackCutsDeuteron->SetEtaRange(-0.77, 0.77);
      AntiTrackCutsDeuteron->SetEtaRange(-0.77, 0.77);
      TrackCutsDeuteron->SetNClsTPC(70);
      AntiTrackCutsDeuteron->SetNClsTPC(70);

      config->SetDeltaEtaMax(0.015);
      config->SetDeltaPhiMax(0.015);

    } else if (suffix == "55") {
      TrackCuts->SetEtaRange(-0.77, 0.77);
      AntiTrackCuts->SetEtaRange(-0.77, 0.77);
      TrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);
      AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);
      TrackCuts->SetNClsTPC(90);
      AntiTrackCuts->SetNClsTPC(90);

      TrackCutsDeuteron->SetEtaRange(-0.77, 0.77);
      AntiTrackCutsDeuteron->SetEtaRange(-0.77, 0.77);
      TrackCutsDeuteron->SetPID(AliPID::kDeuteron, 1.4, 3.3);
      AntiTrackCutsDeuteron->SetPID(AliPID::kDeuteron, 1.4, 3.3);
      TrackCutsDeuteron->SetNClsTPC(90);
      AntiTrackCutsDeuteron->SetNClsTPC(90);

      config->SetDeltaEtaMax(0.015);
      config->SetDeltaPhiMax(0.015);

    } else if (suffix == "56") {
      TrackCuts->SetPtRange(0.4, 2.5);
      AntiTrackCuts->SetPtRange(0.4, 2.5);
      TrackCuts->SetEtaRange(-0.83, 0.83);
      AntiTrackCuts->SetEtaRange(-0.83, 0.83);

      TrackCutsDeuteron->SetPtRange(0.4, 2.0);
      AntiTrackCutsDeuteron->SetPtRange(0.4, 2.0);
      TrackCutsDeuteron->SetEtaRange(-0.83, 0.83);
      AntiTrackCutsDeuteron->SetEtaRange(-0.83, 0.83);

      config->SetDeltaEtaMax(0.015);
      config->SetDeltaPhiMax(0.015);

    } else if (suffix == "57") {
      TrackCuts->SetPtRange(0.4, 2.5);
      AntiTrackCuts->SetPtRange(0.4, 2.5);
      TrackCuts->SetEtaRange(-0.83, 0.83);
      AntiTrackCuts->SetEtaRange(-0.83, 0.83);

      TrackCutsDeuteron->SetPtRange(0.4, 2.0);
      AntiTrackCutsDeuteron->SetPtRange(0.4, 2.0);
      TrackCutsDeuteron->SetEtaRange(-0.83, 0.83);
      AntiTrackCutsDeuteron->SetEtaRange(-0.83, 0.83);

      config->SetDeltaEtaMax(0.019);
      config->SetDeltaPhiMax(0.019);

    } else if (suffix == "58") {
      TrackCuts->SetPtRange(0.6, 2.5);
      AntiTrackCuts->SetPtRange(0.6, 2.5);
      TrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);
      AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);

      TrackCutsDeuteron->SetPtRange(0.6, 2.0);
      AntiTrackCutsDeuteron->SetPtRange(0.6, 2.0);
      TrackCutsDeuteron->SetPID(AliPID::kDeuteron, 1.4, 2.7);
      AntiTrackCutsDeuteron->SetPID(AliPID::kDeuteron, 1.4, 2.7);

      config->SetDeltaEtaMax(0.015);
      config->SetDeltaPhiMax(0.015);

    } else if (suffix == "59") {
      TrackCuts->SetPtRange(0.4, 2.5);
      AntiTrackCuts->SetPtRange(0.4, 2.5);
      TrackCuts->SetEtaRange(-0.77, 0.77);
      AntiTrackCuts->SetEtaRange(-0.77, 0.77);
      TrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);
      AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);

      TrackCutsDeuteron->SetPtRange(0.4, 2.0);
      AntiTrackCutsDeuteron->SetPtRange(0.4, 2.0);
      TrackCutsDeuteron->SetEtaRange(-0.77, 0.77);
      AntiTrackCutsDeuteron->SetEtaRange(-0.77, 0.77);
      TrackCutsDeuteron->SetPID(AliPID::kDeuteron, 1.4, 2.7);
      AntiTrackCutsDeuteron->SetPID(AliPID::kDeuteron, 1.4, 2.7);

      config->SetDeltaEtaMax(0.015);
      config->SetDeltaPhiMax(0.015);

    } else if (suffix == "60") {
      TrackCuts->SetPtRange(0.4, 2.5);
      AntiTrackCuts->SetPtRange(0.4, 2.5);
      TrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);
      AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);
      TrackCuts->SetNClsTPC(70);
      AntiTrackCuts->SetNClsTPC(70);

      TrackCutsDeuteron->SetPtRange(0.4, 2.0);
      AntiTrackCutsDeuteron->SetPtRange(0.4, 2.0);
      TrackCutsDeuteron->SetPID(AliPID::kDeuteron, 1.4, 3.3);
      AntiTrackCutsDeuteron->SetPID(AliPID::kDeuteron, 1.4, 3.3);
      TrackCutsDeuteron->SetNClsTPC(70);
      AntiTrackCutsDeuteron->SetNClsTPC(70);

      config->SetDeltaEtaMax(0.015);
      config->SetDeltaPhiMax(0.015);

    } else if (suffix == "61") {
      TrackCuts->SetEtaRange(-0.77, 0.77);
      AntiTrackCuts->SetEtaRange(-0.77, 0.77);

      TrackCutsDeuteron->SetEtaRange(-0.77, 0.77);
      AntiTrackCutsDeuteron->SetEtaRange(-0.77, 0.77);

      config->SetDeltaEtaMax(0.015);
      config->SetDeltaPhiMax(0.015);

    } else if (suffix == "62") {
      TrackCuts->SetEtaRange(-0.77, 0.77);
      AntiTrackCuts->SetEtaRange(-0.77, 0.77);
      TrackCuts->SetNClsTPC(70);
      AntiTrackCuts->SetNClsTPC(70);

      TrackCutsDeuteron->SetEtaRange(-0.77, 0.77);
      AntiTrackCutsDeuteron->SetEtaRange(-0.77, 0.77);
      TrackCutsDeuteron->SetNClsTPC(70);
      AntiTrackCutsDeuteron->SetNClsTPC(70);

      config->SetDeltaEtaMax(0.015);
      config->SetDeltaPhiMax(0.015);

    } else if (suffix == "63") {
      TrackCuts->SetPtRange(0.4, 2.5);
      AntiTrackCuts->SetPtRange(0.4, 2.5);
      TrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);
      AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);
      TrackCuts->SetNClsTPC(70);
      AntiTrackCuts->SetNClsTPC(70);

      TrackCutsDeuteron->SetPtRange(0.4, 2.0);
      AntiTrackCutsDeuteron->SetPtRange(0.4, 2.0);
      TrackCutsDeuteron->SetPID(AliPID::kDeuteron, 1.4, 2.7);
      AntiTrackCutsDeuteron->SetPID(AliPID::kDeuteron, 1.4, 2.7);
      TrackCutsDeuteron->SetNClsTPC(70);
      AntiTrackCutsDeuteron->SetNClsTPC(70);

      config->SetDeltaEtaMax(0.019);
      config->SetDeltaPhiMax(0.019);

    } else if (suffix == "64") {
      TrackCuts->SetPtRange(0.6, 2.5);
      AntiTrackCuts->SetPtRange(0.6, 2.5);
      TrackCuts->SetEtaRange(-0.83, 0.83);
      AntiTrackCuts->SetEtaRange(-0.83, 0.83);
      TrackCuts->SetNClsTPC(70);
      AntiTrackCuts->SetNClsTPC(70);

      TrackCutsDeuteron->SetPtRange(0.6, 2.0);
      AntiTrackCutsDeuteron->SetPtRange(0.6, 2.0);
      TrackCutsDeuteron->SetEtaRange(-0.83, 0.83);
      AntiTrackCutsDeuteron->SetEtaRange(-0.83, 0.83);
      TrackCutsDeuteron->SetNClsTPC(70);
      AntiTrackCutsDeuteron->SetNClsTPC(70);

      config->SetDeltaEtaMax(0.015);
      config->SetDeltaPhiMax(0.015);

    } else if (suffix == "65") {
      TrackCuts->SetEtaRange(-0.83, 0.83);
      AntiTrackCuts->SetEtaRange(-0.83, 0.83);
      TrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);
      AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);
      TrackCuts->SetNClsTPC(90);
      AntiTrackCuts->SetNClsTPC(90);

      TrackCutsDeuteron->SetEtaRange(-0.83, 0.83);
      AntiTrackCutsDeuteron->SetEtaRange(-0.83, 0.83);
      TrackCutsDeuteron->SetPID(AliPID::kDeuteron, 1.4, 3.3);
      AntiTrackCutsDeuteron->SetPID(AliPID::kDeuteron, 1.4, 3.3);
      TrackCutsDeuteron->SetNClsTPC(90);
      AntiTrackCutsDeuteron->SetNClsTPC(90);

      config->SetDeltaEtaMax(0.019);
      config->SetDeltaPhiMax(0.019);

    } else if (suffix == "66") {
      TrackCuts->SetEtaRange(-0.77, 0.77);
      AntiTrackCuts->SetEtaRange(-0.77, 0.77);
      TrackCuts->SetNClsTPC(90);
      AntiTrackCuts->SetNClsTPC(90);

      TrackCutsDeuteron->SetEtaRange(-0.77, 0.77);
      AntiTrackCutsDeuteron->SetEtaRange(-0.77, 0.77);
      TrackCutsDeuteron->SetNClsTPC(90);
      AntiTrackCutsDeuteron->SetNClsTPC(90);

      config->SetDeltaEtaMax(0.015);
      config->SetDeltaPhiMax(0.015);

    } else if (suffix == "67") {

      TrackCuts->SetPtRange(0.4, 2.5);
      AntiTrackCuts->SetPtRange(0.4, 2.5);
      TrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);
      AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);
      TrackCuts->SetNClsTPC(90);
      AntiTrackCuts->SetNClsTPC(90);

      TrackCutsDeuteron->SetPtRange(0.4, 2.0);
      AntiTrackCutsDeuteron->SetPtRange(0.4, 2.0);
      TrackCutsDeuteron->SetPID(AliPID::kDeuteron, 1.4, 3.3);
      AntiTrackCutsDeuteron->SetPID(AliPID::kDeuteron, 1.4, 3.3);
      TrackCutsDeuteron->SetNClsTPC(90);
      AntiTrackCutsDeuteron->SetNClsTPC(90);

      config->SetDeltaEtaMax(0.015);
      config->SetDeltaPhiMax(0.015);

    } else if (suffix == "68") {
      TrackCuts->SetPtRange(0.4, 2.5);
      AntiTrackCuts->SetPtRange(0.4, 2.5);
      TrackCuts->SetEtaRange(-0.77, 0.77);
      AntiTrackCuts->SetEtaRange(-0.77, 0.77);
      TrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);
      AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);

      TrackCutsDeuteron->SetPtRange(0.4, 2.0);
      AntiTrackCutsDeuteron->SetPtRange(0.4, 2.0);
      TrackCutsDeuteron->SetEtaRange(-0.77, 0.77);
      AntiTrackCutsDeuteron->SetEtaRange(-0.77, 0.77);
      TrackCutsDeuteron->SetPID(AliPID::kDeuteron, 1.4, 3.3);
      AntiTrackCutsDeuteron->SetPID(AliPID::kDeuteron, 1.4, 3.3);

      config->SetDeltaEtaMax(0.015);
      config->SetDeltaPhiMax(0.015);

    } else if (suffix == "69") {
      TrackCuts->SetPtRange(0.6, 2.5);
      AntiTrackCuts->SetPtRange(0.6, 2.5);
      TrackCuts->SetEtaRange(-0.83, 0.83);
      AntiTrackCuts->SetEtaRange(-0.83, 0.83);
      TrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);
      AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);

      TrackCutsDeuteron->SetPtRange(0.6, 2.0);
      AntiTrackCutsDeuteron->SetPtRange(0.6, 2.0);
      TrackCutsDeuteron->SetEtaRange(-0.83, 0.83);
      AntiTrackCutsDeuteron->SetEtaRange(-0.83, 0.83);
      TrackCutsDeuteron->SetPID(AliPID::kDeuteron, 1.4, 3.3);
      AntiTrackCutsDeuteron->SetPID(AliPID::kDeuteron, 1.4, 3.3);

      config->SetDeltaEtaMax(0.015);
      config->SetDeltaPhiMax(0.015);

    } else if (suffix == "70") {
      TrackCuts->SetEtaRange(-0.83, 0.83);
      AntiTrackCuts->SetEtaRange(-0.83, 0.83);
      TrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);
      AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);

      TrackCutsDeuteron->SetEtaRange(-0.83, 0.83);
      AntiTrackCutsDeuteron->SetEtaRange(-0.83, 0.83);
      TrackCutsDeuteron->SetPID(AliPID::kDeuteron, 1.4, 3.3);
      AntiTrackCutsDeuteron->SetPID(AliPID::kDeuteron, 1.4, 3.3);

      config->SetDeltaEtaMax(0.015);
      config->SetDeltaPhiMax(0.015);

    }

  }

//Only pT cut variatons
  if (SystematicpTCutVariation) {
    //Splitting into pure TPC region
    if (suffix == "1") {
      TrackCutsDeuteron->SetPtRange(0.5, 1.4);
      AntiTrackCutsDeuteron->SetPtRange(0.5, 1.4);
      TrackCutsDeuteron->SetPID(AliPID::kDeuteron, 1.4);
      AntiTrackCutsDeuteron->SetPID(AliPID::kDeuteron, 1.4);
      TrackCuts->SetPtRange(0.5, 2.5);
      AntiTrackCuts->SetPtRange(0.5, 2.5);
      //Splitting into pure TPC region
    } else if (suffix == "2") {
      TrackCutsDeuteron->SetPtRange(1.5, 2.5);
      AntiTrackCutsDeuteron->SetPtRange(1.5, 2.5);
      TrackCutsDeuteron->SetPID(AliPID::kDeuteron, 1.4);
      AntiTrackCutsDeuteron->SetPID(AliPID::kDeuteron, 1.4);
      TrackCuts->SetPtRange(0.5, 2.5);
      AntiTrackCuts->SetPtRange(0.5, 2.5);
    } else if (suffix == "3") {

      TrackCutsDeuteron->SetPtRange(1.5, 4.05);
      AntiTrackCutsDeuteron->SetPtRange(1.5, 4.05);
      TrackCutsDeuteron->SetPID(AliPID::kDeuteron, 1.4);
      AntiTrackCutsDeuteron->SetPID(AliPID::kDeuteron, 1.4);
      TrackCuts->SetPtRange(0.5, 2.5);
      AntiTrackCuts->SetPtRange(0.5, 2.5);

    }else if(suffix == "4"){
      TrackCutsDeuteron->SetPtRange(1.5, 4.05);
      AntiTrackCutsDeuteron->SetPtRange(1.5, 4.05);
      TrackCutsDeuteron->SetPID(AliPID::kDeuteron, 1.4);
      AntiTrackCutsDeuteron->SetPID(AliPID::kDeuteron, 1.4);
      TrackCuts->SetPtRange(0.5, 4.05);
      AntiTrackCuts->SetPtRange(0.5,4.05);
    } else if (suffix == "5") {
      TrackCutsDeuteron->SetPtRange(2.0, 4.05);
      AntiTrackCutsDeuteron->SetPtRange(2.0, 4.05);
      TrackCutsDeuteron->SetPID(AliPID::kDeuteron, 1.4);
      AntiTrackCutsDeuteron->SetPID(AliPID::kDeuteron, 1.4);
      TrackCuts->SetPtRange(0.5, 4.05);
      AntiTrackCuts->SetPtRange(0.5, 4.05);

    }else if (suffix == "6") {
      TrackCutsDeuteron->SetPtRange(2.0, 4.05);
      AntiTrackCutsDeuteron->SetPtRange(2.0, 4.05);
      TrackCutsDeuteron->SetPID(AliPID::kDeuteron, 1.4);
      AntiTrackCutsDeuteron->SetPID(AliPID::kDeuteron, 1.4);
      TrackCuts->SetPtRange(0.5, 2.5);
      AntiTrackCuts->SetPtRange(0.5, 2.5);

    } else if (suffix == "7") {
      TrackCutsDeuteron->SetPID(AliPID::kDeuteron, 1.4);
      AntiTrackCutsDeuteron->SetPID(AliPID::kDeuteron, 1.4);

      TrackCutsDeuteron->SetPtRange(0.5, 2.5);
      AntiTrackCutsDeuteron->SetPtRange(0.5, 2.5);
      TrackCuts->SetPtRange(0.5, 2.5);
      AntiTrackCuts->SetPtRange(0.5, 2.5);

    } else if (suffix == "8") {
      TrackCutsDeuteron->SetPID(AliPID::kDeuteron, 1.4);
      AntiTrackCutsDeuteron->SetPID(AliPID::kDeuteron, 1.4);

      TrackCutsDeuteron->SetPtRange(0.5, 4.05);
      AntiTrackCutsDeuteron->SetPtRange(0.5, 4.05);
      TrackCuts->SetPtRange(0.5, 2.5);
      AntiTrackCuts->SetPtRange(0.5, 2.5);
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
  task->SetMCTruth(IsMCTruth);
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

