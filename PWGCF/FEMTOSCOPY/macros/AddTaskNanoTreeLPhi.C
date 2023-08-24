#ifndef __CINT__
#include <vector>
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskNanoTreeLPhi.h"
#include "AliFemtoDreamEventCuts.h"
#include "AliFemtoDreamTrackCuts.h"
#include "AliFemtoDreamv0Cuts.h"
#include "AliFemtoDreamCollConfig.h"
#endif

AliAnalysisTaskSE *AddTaskNanoTreeLPhi(bool isMC = false,
                                       TString CentEst = "kInt7",
                                       const char *cutVariation = "0")
{
  TString suffix = TString::Format("%s", cutVariation);

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  if (!mgr)
  {
    printf("No analysis manager to connect to!\n");
    return nullptr;
  }
  if (!mgr->GetInputEventHandler())
  {
    printf("This task requires an input event handler!\n");
    return nullptr;
  }

  AliFemtoDreamEventCuts *evtCuts = AliFemtoDreamEventCuts::StandardCutsRun2();
  evtCuts->CleanUpMult(false, false, false, true);

  //Lambda Cuts
  AliFemtoDreamv0Cuts *v0Cuts = AliFemtoDreamv0Cuts::LambdaCuts(isMC, true, true);
  AliFemtoDreamTrackCuts *Posv0Daug = AliFemtoDreamTrackCuts::DecayProtonCuts(isMC, true, false); //PileUpRej, false
  AliFemtoDreamTrackCuts *Negv0Daug = AliFemtoDreamTrackCuts::DecayPionCuts(isMC, true, false);
  Posv0Daug->SetEtaRange(-0.9, 0.9);
  Posv0Daug->SetNClsTPC(60);
  Posv0Daug->SetPID(AliPID::kProton, 999., 6);
  Negv0Daug->SetEtaRange(-0.9, 0.9);
  Negv0Daug->SetNClsTPC(60);
  Negv0Daug->SetPID(AliPID::kPion, 999., 6);
  v0Cuts->SetPosDaugterTrackCuts(Posv0Daug);
  v0Cuts->SetNegDaugterTrackCuts(Negv0Daug);
  v0Cuts->SetPDGCodePosDaug(2212); //Proton
  v0Cuts->SetPDGCodeNegDaug(211);  //Pion
  v0Cuts->SetPDGCodev0(3122);      //Lambda

  v0Cuts->SetCutDCADaugToPrimVtx(0.03);
  v0Cuts->SetPtRange(0.2, 999.);
  v0Cuts->SetCutMaxDecayVtx(110);
  v0Cuts->SetCutTransverseRadius(0.15, 110);
  v0Cuts->SetCutDCADaugTov0Vtx(1.7);
  v0Cuts->SetCutCPA(0.97);

  AliFemtoDreamv0Cuts *Antiv0Cuts = AliFemtoDreamv0Cuts::LambdaCuts(isMC, true, true);
  AliFemtoDreamTrackCuts *PosAntiv0Daug = AliFemtoDreamTrackCuts::DecayPionCuts(isMC, true, false);
  PosAntiv0Daug->SetCutCharge(1);
  PosAntiv0Daug->SetEtaRange(-0.9, 0.9);
  PosAntiv0Daug->SetNClsTPC(60);
  PosAntiv0Daug->SetPID(AliPID::kPion, 999., 6);
  AliFemtoDreamTrackCuts *NegAntiv0Daug = AliFemtoDreamTrackCuts::DecayProtonCuts(isMC, true, false);
  NegAntiv0Daug->SetCutCharge(-1);
  NegAntiv0Daug->SetEtaRange(-0.9, 0.9);
  NegAntiv0Daug->SetNClsTPC(60);
  NegAntiv0Daug->SetPID(AliPID::kProton, 999., 6);

  Antiv0Cuts->SetPosDaugterTrackCuts(PosAntiv0Daug);
  Antiv0Cuts->SetNegDaugterTrackCuts(NegAntiv0Daug);
  Antiv0Cuts->SetPDGCodePosDaug(211);  //Pion
  Antiv0Cuts->SetPDGCodeNegDaug(2212); //Proton
  Antiv0Cuts->SetPDGCodev0(-3122);     //Lambda

  Antiv0Cuts->SetPtRange(0.2, 999.);
  Antiv0Cuts->SetCutMaxDecayVtx(110);
  Antiv0Cuts->SetCutTransverseRadius(0.15, 110);
  Antiv0Cuts->SetCutDCADaugTov0Vtx(1.7);
  Antiv0Cuts->SetCutCPA(0.97);

  AliFemtoDreamTrackCuts *TrackPosKaonCuts =
      AliFemtoDreamTrackCuts::PrimKaonCuts(isMC, true, false, false);
  TrackPosKaonCuts->SetCutCharge(1);
  TrackPosKaonCuts->SetFilterBit(128);

  AliFemtoDreamTrackCuts *TrackNegKaonCuts =
      AliFemtoDreamTrackCuts::PrimKaonCuts(isMC, true, false, false);
  TrackNegKaonCuts->SetCutCharge(-1);
  TrackNegKaonCuts->SetFilterBit(128);

  TrackPosKaonCuts->SetDCAVtxZ(0.4);
  TrackNegKaonCuts->SetDCAVtxZ(0.4);
  TrackPosKaonCuts->SetDCAVtxXY(0.8);
  TrackNegKaonCuts->SetDCAVtxXY(0.8);

  if (suffix != "0" && suffix != "999")
  {
    TrackPosKaonCuts->SetMinimalBooking(true);
    TrackNegKaonCuts->SetMinimalBooking(true);
  }

  AliFemtoDreamv0Cuts *TrackCutsPhi = new AliFemtoDreamv0Cuts();
  TrackCutsPhi->SetIsMonteCarlo(isMC);
  TrackCutsPhi->SetAxisInvMassPlots(400, 0.95, 2);
  TrackCutsPhi->SetCutInvMass(0.008);
  AliFemtoDreamTrackCuts *dummyCutsPos = new AliFemtoDreamTrackCuts();
  dummyCutsPos->SetIsMonteCarlo(isMC);
  AliFemtoDreamTrackCuts *dummyCutsNeg = new AliFemtoDreamTrackCuts();
  dummyCutsNeg->SetIsMonteCarlo(isMC);
  TrackCutsPhi->SetPosDaugterTrackCuts(dummyCutsPos);
  TrackCutsPhi->SetNegDaugterTrackCuts(dummyCutsNeg);
  TrackCutsPhi->SetPDGCodePosDaug(321);
  TrackCutsPhi->SetPDGCodeNegDaug(321);
  TrackCutsPhi->SetPDGCodev0(333);

  double Phimass = TDatabasePDG::Instance()->GetParticle(333)->Mass();

  if (suffix != "0")
  {
    evtCuts->SetMinimalBooking(true);
    TrackPosKaonCuts->SetMinimalBooking(true);
    TrackNegKaonCuts->SetMinimalBooking(true);
    v0Cuts->SetMinimalBooking(true);
    Antiv0Cuts->SetMinimalBooking(true);
    TrackCutsPhi->SetMinimalBooking(true);
  }

  // Now we define stuff we want for our Particle collection
  // Thanks, CINT - will not compile due to an illegal constructor
  // std::vector<int> PDGParticles ={2212,2212,3122,3122,3312,3312};
  // First we need to tell him about the particles we mix, from the
  // PDG code the mass is obtained.
  std::vector<int> PDGParticles;
  PDGParticles.push_back(3122); // 0 Lambda
  PDGParticles.push_back(3122); // 1 antiLambda
  PDGParticles.push_back(333);  // 2 Phi particle
  PDGParticles.push_back(321);  // 3 Kaon Plus
  PDGParticles.push_back(321);  // 4 Kaon Minus

  // We need to set the ZVtx bins
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
  // The Multiplicity bins are set here
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

  // Number of bins
  std::vector<int> NBins;
  // minimum k* value
  std::vector<float> kMin;
  // maximum k* value
  std::vector<float> kMax;
  // pair QA extended
  std::vector<int> pairQA;
  std::vector<bool> closeRejection;
  /// Pairs
  /// ΛΛ          0
  /// ΛantiΛ      1
  /// Λφ          2
  /// ΛK+         3
  /// ΛK-         4
  /// antiΛantiΛ  5
  /// antiΛφ      6
  /// antiΛK+     7
  /// antiΛK-     8
  /// φφ          9
  /// φK+         10
  /// φΚ-         11
  /// K+K+        12
  /// K+K-        13
  /// K-K-        14

  const int npairs = 15;
  for (int i = 0; i < npairs; i++)
  {
    NBins.push_back(750);
    closeRejection.push_back(false);
    kMin.push_back(0.);
    kMax.push_back(3.);
    pairQA.push_back(0);
  }

  if (suffix != "0")
  {
    pairQA[0] = 22;
    pairQA[1] = 22;
    pairQA[2] = 22;
    pairQA[3] = 21;
    pairQA[4] = 21;
    pairQA[5] = 22;
    pairQA[6] = 22;
    pairQA[7] = 21;
    pairQA[8] = 21;
  }
  else
  {
    pairQA[0] = 22;
    pairQA[1] = 22;
    pairQA[2] = 22;
    pairQA[3] = 21;
    pairQA[4] = 21;
    pairQA[5] = 22;
    pairQA[6] = 22;
    pairQA[7] = 21;
    pairQA[8] = 21;
    pairQA[12] = 11;
    pairQA[13] = 11;
    pairQA[14] = 11;
  }

  AliFemtoDreamCollConfig *config =
      new AliFemtoDreamCollConfig("Femto", "Femto");
  config->SetPtQA(true);
  config->SetMassQA(true);
  config->SetExtendedQAPairs(pairQA);
  config->SetZBins(ZVtxBins);
  config->SetMultBins(MultBins);
  config->SetMultBinning(true);
  config->SetPDGCodes(PDGParticles);
  config->SetNBinsHist(NBins);
  config->SetMinKRel(kMin);
  config->SetMaxKRel(kMax);
  config->SetUseEventMixing(true);
  config->SetMixingDepth(30);
  config->SetMultiplicityEstimator(AliFemtoDreamEvent::kRef08);

  if (suffix == "0")
  {
    config->SetkTBinning(true);
    config->SetMultBinning(true);
    config->SetmTBinning(true);
  }

  if (isMC)
  {
    config->SetMomentumResolution(true);
  }
  else
  {
    std::cout << "You are trying to request the Momentum Resolution without MC Info; fix it wont work! \n";
  }

  // Variation cuts
  const float KaonPtlow = 0.075;
  const float KaonPtup = 0.225;
  const float KaonEtaLow = 0.75;
  const float KaonEtaUp = 0.85;
  const float KaonNsigmaLow = 4.25;
  const float KaonNsigmaUp = 5.75;
  const float KaonNClsLow = 70;
  const float KaonNClsUp = 90;
  const float KaonPtMax = 999;

  /// Systematic variations (taken from pφ and ΛΚ)
  if (suffix != "0")
  {
    if (suffix == "1")
    {
      TrackPosKaonCuts->SetNClsTPC(KaonNClsUp);
      TrackNegKaonCuts->SetNClsTPC(KaonNClsUp);

      v0Cuts->SetCutCPA(0.995);
      Antiv0Cuts->SetCutCPA(0.995);

      Posv0Daug->SetNClsTPC(80);
      Negv0Daug->SetNClsTPC(80);
      PosAntiv0Daug->SetNClsTPC(80);
      NegAntiv0Daug->SetNClsTPC(80);

      v0Cuts->SetCutDCADaugToPrimVtx(0.06);
      Antiv0Cuts->SetCutDCADaugToPrimVtx(0.06);
    }
    else if (suffix == "2")
    {
      TrackPosKaonCuts->SetPtRange(KaonPtlow, KaonPtMax);
      TrackNegKaonCuts->SetPtRange(KaonPtlow, KaonPtMax);

      TrackPosKaonCuts->SetEtaRange(-KaonEtaLow, KaonEtaLow);
      TrackNegKaonCuts->SetEtaRange(-KaonEtaLow, KaonEtaLow);

      Posv0Daug->SetPID(AliPID::kProton, 999.9, 4);
      Negv0Daug->SetPID(AliPID::kPion, 999.9, 4);
      PosAntiv0Daug->SetPID(AliPID::kPion, 999.9, 4);
      NegAntiv0Daug->SetPID(AliPID::kProton, 999.9, 4);
    }
    else if (suffix == "3")
    {
      v0Cuts->SetCutCPA(0.995);
      Antiv0Cuts->SetCutCPA(0.995);

      Posv0Daug->SetEtaRange(-0.77, 0.77);
      Negv0Daug->SetEtaRange(-0.77, 0.77);
      PosAntiv0Daug->SetEtaRange(-0.77, 0.77);
      NegAntiv0Daug->SetEtaRange(-0.77, 0.77);
      v0Cuts->SetCutDCADaugTov0Vtx(1.2);
      Antiv0Cuts->SetCutDCADaugTov0Vtx(1.2);
    }
    else if (suffix == "4")
    {
      TrackPosKaonCuts->SetEtaRange(-KaonEtaLow, KaonEtaLow);
      TrackNegKaonCuts->SetEtaRange(-KaonEtaLow, KaonEtaLow);

      TrackPosKaonCuts->SetNClsTPC(KaonNClsLow);
      TrackNegKaonCuts->SetNClsTPC(KaonNClsLow);

      Posv0Daug->SetNClsTPC(80);
      Negv0Daug->SetNClsTPC(80);
      PosAntiv0Daug->SetNClsTPC(80);
      NegAntiv0Daug->SetNClsTPC(80);
      v0Cuts->SetCutDCADaugToPrimVtx(0.06);
      Antiv0Cuts->SetCutDCADaugToPrimVtx(0.06);
    }
    else if (suffix == "5")
    {
      TrackPosKaonCuts->SetPtRange(KaonPtlow, KaonPtMax);
      TrackNegKaonCuts->SetPtRange(KaonPtlow, KaonPtMax);

      TrackPosKaonCuts->SetEtaRange(-KaonEtaLow, KaonEtaLow);
      TrackNegKaonCuts->SetEtaRange(-KaonEtaLow, KaonEtaLow);

      TrackPosKaonCuts->SetPID(AliPID::kKaon, 0.4, KaonNsigmaUp);
      TrackNegKaonCuts->SetPID(AliPID::kKaon, 0.4, KaonNsigmaUp);

      v0Cuts->SetCutCPA(0.995);
      Antiv0Cuts->SetCutCPA(0.995);

      Posv0Daug->SetEtaRange(-0.83, 0.83);
      Negv0Daug->SetEtaRange(-0.83, 0.83);
      PosAntiv0Daug->SetEtaRange(-0.83, 0.83);
      NegAntiv0Daug->SetEtaRange(-0.83, 0.83);

      v0Cuts->SetCutDCADaugTov0Vtx(1.2);
      Antiv0Cuts->SetCutDCADaugTov0Vtx(1.2);
    }
    else if (suffix == "6")
    {
      TrackPosKaonCuts->SetPID(AliPID::kKaon, 0.4, KaonNsigmaLow);
      TrackNegKaonCuts->SetPID(AliPID::kKaon, 0.4, KaonNsigmaLow);

      TrackPosKaonCuts->SetNClsTPC(KaonNClsUp);
      TrackNegKaonCuts->SetNClsTPC(KaonNClsUp);

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
    }
    else if (suffix == "7")
    {
      TrackPosKaonCuts->SetEtaRange(-KaonEtaUp, KaonEtaUp);
      TrackNegKaonCuts->SetEtaRange(-KaonEtaUp, KaonEtaUp);

      TrackPosKaonCuts->SetPID(AliPID::kKaon, 0.4, KaonNsigmaUp);
      TrackNegKaonCuts->SetPID(AliPID::kKaon, 0.4, KaonNsigmaUp);

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
    }
    else if (suffix == "8")
    {
      TrackPosKaonCuts->SetPtRange(KaonPtlow, KaonPtMax);
      TrackNegKaonCuts->SetPtRange(KaonPtlow, KaonPtMax);

      TrackPosKaonCuts->SetEtaRange(-KaonEtaLow, KaonEtaLow);
      TrackNegKaonCuts->SetEtaRange(-KaonEtaLow, KaonEtaLow);

      TrackPosKaonCuts->SetNClsTPC(KaonNClsUp);
      TrackNegKaonCuts->SetNClsTPC(KaonNClsUp);

      Posv0Daug->SetEtaRange(-0.83, 0.83);
      Negv0Daug->SetEtaRange(-0.83, 0.83);
      PosAntiv0Daug->SetEtaRange(-0.83, 0.83);
      NegAntiv0Daug->SetEtaRange(-0.83, 0.83);
    }
    else if (suffix == "9")
    {
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
    }
    else if (suffix == "10")
    {
      TrackPosKaonCuts->SetPtRange(KaonPtlow, KaonPtMax);
      TrackNegKaonCuts->SetPtRange(KaonPtlow, KaonPtMax);

      TrackPosKaonCuts->SetEtaRange(-KaonEtaUp, KaonEtaUp);
      TrackNegKaonCuts->SetEtaRange(-KaonEtaUp, KaonEtaUp);

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
    }
    else if (suffix == "11")
    {
      TrackPosKaonCuts->SetEtaRange(-KaonEtaLow, KaonEtaLow);
      TrackNegKaonCuts->SetEtaRange(-KaonEtaLow, KaonEtaLow);

      TrackPosKaonCuts->SetPID(AliPID::kKaon, 0.4, KaonNsigmaLow);
      TrackNegKaonCuts->SetPID(AliPID::kKaon, 0.4, KaonNsigmaLow);

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
    }
    else if (suffix == "12")
    {
      TrackPosKaonCuts->SetPtRange(KaonPtlow, KaonPtMax);
      TrackNegKaonCuts->SetPtRange(KaonPtlow, KaonPtMax);

      TrackPosKaonCuts->SetEtaRange(-KaonEtaUp, KaonEtaUp);
      TrackNegKaonCuts->SetEtaRange(-KaonEtaUp, KaonEtaUp);

      TrackPosKaonCuts->SetPID(AliPID::kKaon, 0.4, KaonNsigmaLow);
      TrackNegKaonCuts->SetPID(AliPID::kKaon, 0.4, KaonNsigmaLow);

      Posv0Daug->SetPID(AliPID::kProton, 999.9, 4);
      Negv0Daug->SetPID(AliPID::kPion, 999.9, 4);
      PosAntiv0Daug->SetPID(AliPID::kPion, 999.9, 4);
      NegAntiv0Daug->SetPID(AliPID::kProton, 999.9, 4);

      Posv0Daug->SetEtaRange(-0.77, 0.77);
      Negv0Daug->SetEtaRange(-0.77, 0.77);
      PosAntiv0Daug->SetEtaRange(-0.77, 0.77);
      NegAntiv0Daug->SetEtaRange(-0.77, 0.77);
    }
    else if (suffix == "13")
    {
      TrackPosKaonCuts->SetEtaRange(-KaonEtaUp, KaonEtaUp);
      TrackNegKaonCuts->SetEtaRange(-KaonEtaUp, KaonEtaUp);

      TrackPosKaonCuts->SetPID(AliPID::kKaon, 0.4, KaonNsigmaLow);
      TrackNegKaonCuts->SetPID(AliPID::kKaon, 0.4, KaonNsigmaLow);

      TrackPosKaonCuts->SetNClsTPC(KaonNClsUp);
      TrackNegKaonCuts->SetNClsTPC(KaonNClsUp);

      Posv0Daug->SetEtaRange(-0.83, 0.83);
      Negv0Daug->SetEtaRange(-0.83, 0.83);
      PosAntiv0Daug->SetEtaRange(-0.83, 0.83);
      NegAntiv0Daug->SetEtaRange(-0.83, 0.83);

      v0Cuts->SetCutDCADaugTov0Vtx(1.2);
      Antiv0Cuts->SetCutDCADaugTov0Vtx(1.2);

      v0Cuts->SetCutDCADaugToPrimVtx(0.06);
      Antiv0Cuts->SetCutDCADaugToPrimVtx(0.06);
    }
    else if (suffix == "14")
    {
      TrackPosKaonCuts->SetPtRange(KaonPtlow, KaonPtMax);
      TrackNegKaonCuts->SetPtRange(KaonPtlow, KaonPtMax);

      TrackPosKaonCuts->SetPID(AliPID::kKaon, 0.4, KaonNsigmaLow);
      TrackNegKaonCuts->SetPID(AliPID::kKaon, 0.4, KaonNsigmaLow);

      Posv0Daug->SetNClsTPC(80);
      Negv0Daug->SetNClsTPC(80);
      PosAntiv0Daug->SetNClsTPC(80);
      NegAntiv0Daug->SetNClsTPC(80);

      Posv0Daug->SetEtaRange(-0.83, 0.83);
      Negv0Daug->SetEtaRange(-0.83, 0.83);
      PosAntiv0Daug->SetEtaRange(-0.83, 0.83);
      NegAntiv0Daug->SetEtaRange(-0.83, 0.83);
    }
    else if (suffix == "15")
    {
      TrackPosKaonCuts->SetPtRange(KaonPtup, KaonPtMax);
      TrackNegKaonCuts->SetPtRange(KaonPtup, KaonPtMax);

      TrackPosKaonCuts->SetEtaRange(-KaonEtaUp, KaonEtaUp);
      TrackNegKaonCuts->SetEtaRange(-KaonEtaUp, KaonEtaUp);

      TrackPosKaonCuts->SetPID(AliPID::kKaon, 0.4, KaonNsigmaLow);
      TrackNegKaonCuts->SetPID(AliPID::kKaon, 0.4, KaonNsigmaLow);

      TrackPosKaonCuts->SetNClsTPC(KaonNClsUp);
      TrackNegKaonCuts->SetNClsTPC(KaonNClsUp);

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
    }
    else if (suffix == "16")
    {
      TrackPosKaonCuts->SetPtRange(KaonPtlow, KaonPtMax);
      TrackNegKaonCuts->SetPtRange(KaonPtlow, KaonPtMax);

      TrackPosKaonCuts->SetPID(AliPID::kKaon, 0.4, KaonNsigmaUp);
      TrackNegKaonCuts->SetPID(AliPID::kKaon, 0.4, KaonNsigmaUp);

      TrackPosKaonCuts->SetNClsTPC(KaonNClsLow);
      TrackNegKaonCuts->SetNClsTPC(KaonNClsLow);

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
    }
    else if (suffix == "17")
    {
      TrackPosKaonCuts->SetPtRange(KaonPtlow, KaonPtMax);
      TrackNegKaonCuts->SetPtRange(KaonPtlow, KaonPtMax);

      TrackPosKaonCuts->SetEtaRange(-KaonEtaUp, KaonEtaUp);
      TrackNegKaonCuts->SetEtaRange(-KaonEtaUp, KaonEtaUp);

      TrackPosKaonCuts->SetPID(AliPID::kKaon, 0.4, KaonNsigmaLow);
      TrackNegKaonCuts->SetPID(AliPID::kKaon, 0.4, KaonNsigmaLow);

      v0Cuts->SetCutCPA(0.995);
      Antiv0Cuts->SetCutCPA(0.995);

      Posv0Daug->SetEtaRange(-0.77, 0.77);
      Negv0Daug->SetEtaRange(-0.77, 0.77);
      PosAntiv0Daug->SetEtaRange(-0.77, 0.77);
      NegAntiv0Daug->SetEtaRange(-0.77, 0.77);

      v0Cuts->SetCutDCADaugToPrimVtx(0.06);
      Antiv0Cuts->SetCutDCADaugToPrimVtx(0.06);
    }
    else if (suffix == "18")
    {
      TrackPosKaonCuts->SetPtRange(KaonPtlow, KaonPtMax);
      TrackNegKaonCuts->SetPtRange(KaonPtlow, KaonPtMax);

      TrackPosKaonCuts->SetEtaRange(-KaonEtaUp, KaonEtaUp);
      TrackNegKaonCuts->SetEtaRange(-KaonEtaUp, KaonEtaUp);

      TrackPosKaonCuts->SetPID(AliPID::kKaon, 0.4, KaonNsigmaLow);
      TrackNegKaonCuts->SetPID(AliPID::kKaon, 0.4, KaonNsigmaLow);

      TrackPosKaonCuts->SetNClsTPC(KaonNClsUp);
      TrackNegKaonCuts->SetNClsTPC(KaonNClsUp);

      Posv0Daug->SetNClsTPC(80);
      Negv0Daug->SetNClsTPC(80);
      PosAntiv0Daug->SetNClsTPC(80);
      NegAntiv0Daug->SetNClsTPC(80);

      Posv0Daug->SetEtaRange(-0.77, 0.77);
      Negv0Daug->SetEtaRange(-0.77, 0.77);
      PosAntiv0Daug->SetEtaRange(-0.77, 0.77);
      NegAntiv0Daug->SetEtaRange(-0.77, 0.77);
    }
    else if (suffix == "19")
    {
      TrackPosKaonCuts->SetPID(AliPID::kKaon, 0.4, KaonNsigmaLow);
      TrackNegKaonCuts->SetPID(AliPID::kKaon, 0.4, KaonNsigmaLow);

      TrackPosKaonCuts->SetNClsTPC(KaonNClsLow);
      TrackNegKaonCuts->SetNClsTPC(KaonNClsLow);

      v0Cuts->SetCutCPA(0.995);
      Antiv0Cuts->SetCutCPA(0.995);

      Posv0Daug->SetNClsTPC(80);
      Negv0Daug->SetNClsTPC(80);
      PosAntiv0Daug->SetNClsTPC(80);
      NegAntiv0Daug->SetNClsTPC(80);

      v0Cuts->SetCutDCADaugTov0Vtx(1.2);
      Antiv0Cuts->SetCutDCADaugTov0Vtx(1.2);
    }
    else if (suffix == "20")
    {
      TrackPosKaonCuts->SetPtRange(KaonPtlow, KaonPtMax);
      TrackNegKaonCuts->SetPtRange(KaonPtlow, KaonPtMax);

      TrackPosKaonCuts->SetNClsTPC(KaonNClsUp);
      TrackNegKaonCuts->SetNClsTPC(KaonNClsUp);

      Posv0Daug->SetNClsTPC(80);
      Negv0Daug->SetNClsTPC(80);
      PosAntiv0Daug->SetNClsTPC(80);
      NegAntiv0Daug->SetNClsTPC(80);

      Posv0Daug->SetEtaRange(-0.83, 0.83);
      Negv0Daug->SetEtaRange(-0.83, 0.83);
      PosAntiv0Daug->SetEtaRange(-0.83, 0.83);
      NegAntiv0Daug->SetEtaRange(-0.83, 0.83);
    }
    else if (suffix == "21")
    {

      TrackPosKaonCuts->SetEtaRange(-KaonEtaUp, KaonEtaUp);
      TrackNegKaonCuts->SetEtaRange(-KaonEtaUp, KaonEtaUp);

      TrackPosKaonCuts->SetPID(AliPID::kKaon, 0.4, KaonNsigmaUp);
      TrackNegKaonCuts->SetPID(AliPID::kKaon, 0.4, KaonNsigmaUp);

      TrackPosKaonCuts->SetNClsTPC(KaonNClsLow);
      TrackNegKaonCuts->SetNClsTPC(KaonNClsLow);

      Posv0Daug->SetNClsTPC(80);
      Negv0Daug->SetNClsTPC(80);
      PosAntiv0Daug->SetNClsTPC(80);
      NegAntiv0Daug->SetNClsTPC(80);

      Posv0Daug->SetEtaRange(-0.77, 0.77);
      Negv0Daug->SetEtaRange(-0.77, 0.77);
      PosAntiv0Daug->SetEtaRange(-0.77, 0.77);
      NegAntiv0Daug->SetEtaRange(-0.77, 0.77);
    }
    else if (suffix == "22")
    {
      TrackPosKaonCuts->SetPtRange(KaonPtlow, KaonPtMax);
      TrackNegKaonCuts->SetPtRange(KaonPtlow, KaonPtMax);

      TrackPosKaonCuts->SetEtaRange(-KaonEtaLow, KaonEtaLow);
      TrackNegKaonCuts->SetEtaRange(-KaonEtaLow, KaonEtaLow);

      TrackPosKaonCuts->SetNClsTPC(KaonNClsLow);
      TrackNegKaonCuts->SetNClsTPC(KaonNClsLow);

      v0Cuts->SetCutCPA(0.995);
      Antiv0Cuts->SetCutCPA(0.995);

      Posv0Daug->SetPID(AliPID::kProton, 999.9, 4);
      Negv0Daug->SetPID(AliPID::kPion, 999.9, 4);
      PosAntiv0Daug->SetPID(AliPID::kPion, 999.9, 4);
      NegAntiv0Daug->SetPID(AliPID::kProton, 999.9, 4);
    }
    else if (suffix == "23")
    {
      TrackPosKaonCuts->SetEtaRange(-KaonEtaLow, KaonEtaLow);
      TrackNegKaonCuts->SetEtaRange(-KaonEtaLow, KaonEtaLow);

      TrackPosKaonCuts->SetPID(AliPID::kKaon, 0.4, KaonNsigmaLow);
      TrackNegKaonCuts->SetPID(AliPID::kKaon, 0.4, KaonNsigmaLow);

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
    }
    else if (suffix == "24")
    {
      TrackPosKaonCuts->SetEtaRange(-KaonEtaLow, KaonEtaLow);
      TrackNegKaonCuts->SetEtaRange(-KaonEtaLow, KaonEtaLow);

      TrackPosKaonCuts->SetPID(AliPID::kKaon, 0.4, KaonNsigmaUp);
      TrackNegKaonCuts->SetPID(AliPID::kKaon, 0.4, KaonNsigmaUp);

      TrackPosKaonCuts->SetNClsTPC(KaonNClsUp);
      TrackNegKaonCuts->SetNClsTPC(KaonNClsUp);

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
    }
    else if (suffix == "25")
    {
      TrackPosKaonCuts->SetPtRange(KaonPtlow, KaonPtMax);
      TrackNegKaonCuts->SetPtRange(KaonPtlow, KaonPtMax);

      TrackPosKaonCuts->SetEtaRange(-0.83, 0.83);
      TrackNegKaonCuts->SetEtaRange(-0.83, 0.83);

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
    }
    else if (suffix == "26")
    {
      TrackPosKaonCuts->SetPtRange(KaonPtup, KaonPtMax);
      TrackNegKaonCuts->SetPtRange(KaonPtup, KaonPtMax);

      TrackPosKaonCuts->SetPID(AliPID::kKaon, 0.4, KaonNsigmaLow);
      TrackNegKaonCuts->SetPID(AliPID::kKaon, 0.4, KaonNsigmaLow);

      Posv0Daug->SetEtaRange(-0.77, 0.77);
      Negv0Daug->SetEtaRange(-0.77, 0.77);
      PosAntiv0Daug->SetEtaRange(-0.77, 0.77);
      NegAntiv0Daug->SetEtaRange(-0.77, 0.77);
    }
    else if (suffix == "27")
    {
      TrackPosKaonCuts->SetEtaRange(-KaonEtaLow, KaonEtaLow);
      TrackNegKaonCuts->SetEtaRange(-KaonEtaLow, KaonEtaLow);

      TrackPosKaonCuts->SetPID(AliPID::kKaon, 0.4, KaonNsigmaUp);
      TrackNegKaonCuts->SetPID(AliPID::kKaon, 0.4, KaonNsigmaUp);

      TrackPosKaonCuts->SetNClsTPC(KaonNClsUp);
      TrackNegKaonCuts->SetNClsTPC(KaonNClsUp);

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
    }
    else if (suffix == "28")
    {
      TrackPosKaonCuts->SetNClsTPC(KaonNClsUp);
      TrackNegKaonCuts->SetNClsTPC(KaonNClsUp);

      v0Cuts->SetCutCPA(0.995);
      Antiv0Cuts->SetCutCPA(0.995);

      v0Cuts->SetCutDCADaugToPrimVtx(0.06);
      Antiv0Cuts->SetCutDCADaugToPrimVtx(0.06);
    }
    else if (suffix == "29")
    {
      TrackPosKaonCuts->SetPtRange(KaonPtlow, KaonPtMax);
      TrackNegKaonCuts->SetPtRange(KaonPtlow, KaonPtMax);

      TrackPosKaonCuts->SetEtaRange(-KaonEtaLow, KaonEtaLow);
      TrackNegKaonCuts->SetEtaRange(-KaonEtaLow, KaonEtaLow);

      TrackPosKaonCuts->SetPID(AliPID::kKaon, 0.4, KaonNsigmaLow);
      TrackNegKaonCuts->SetPID(AliPID::kKaon, 0.4, KaonNsigmaLow);

      Posv0Daug->SetNClsTPC(80);
      Negv0Daug->SetNClsTPC(80);
      PosAntiv0Daug->SetNClsTPC(80);
      NegAntiv0Daug->SetNClsTPC(80);

      v0Cuts->SetCutDCADaugToPrimVtx(0.06);
      Antiv0Cuts->SetCutDCADaugToPrimVtx(0.06);
    }
    else if (suffix == "30")
    {
      TrackPosKaonCuts->SetEtaRange(-KaonEtaUp, KaonEtaUp);
      TrackNegKaonCuts->SetEtaRange(-KaonEtaUp, KaonEtaUp);

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
    }
    else if (suffix == "31")
    {
      TrackPosKaonCuts->SetPtRange(KaonPtlow, KaonPtMax);
      TrackNegKaonCuts->SetPtRange(KaonPtlow, KaonPtMax);

      TrackPosKaonCuts->SetPID(AliPID::kKaon, 0.4, KaonNsigmaUp);
      TrackNegKaonCuts->SetPID(AliPID::kKaon, 0.4, KaonNsigmaUp);

      TrackPosKaonCuts->SetNClsTPC(KaonNClsLow);
      TrackNegKaonCuts->SetNClsTPC(KaonNClsLow);

      v0Cuts->SetCutCPA(0.995);
      Antiv0Cuts->SetCutCPA(0.995);

      Posv0Daug->SetEtaRange(-0.83, 0.83);
      Negv0Daug->SetEtaRange(-0.83, 0.83);
      PosAntiv0Daug->SetEtaRange(-0.83, 0.83);
      NegAntiv0Daug->SetEtaRange(-0.83, 0.83);

      v0Cuts->SetCutDCADaugTov0Vtx(1.2);
      Antiv0Cuts->SetCutDCADaugTov0Vtx(1.2);
    }
    else if (suffix == "32")
    {

      TrackPosKaonCuts->SetEtaRange(-KaonEtaLow, KaonEtaLow);
      TrackNegKaonCuts->SetEtaRange(-KaonEtaLow, KaonEtaLow);

      TrackPosKaonCuts->SetPID(AliPID::kKaon, 0.4, KaonNsigmaLow);
      TrackNegKaonCuts->SetPID(AliPID::kKaon, 0.4, KaonNsigmaLow);

      TrackPosKaonCuts->SetNClsTPC(KaonNClsUp);
      TrackNegKaonCuts->SetNClsTPC(KaonNClsUp);

      v0Cuts->SetCutCPA(0.995);
      Antiv0Cuts->SetCutCPA(0.995);

      Posv0Daug->SetNClsTPC(80);
      Negv0Daug->SetNClsTPC(80);
      PosAntiv0Daug->SetNClsTPC(80);
      NegAntiv0Daug->SetNClsTPC(80);

      v0Cuts->SetCutDCADaugToPrimVtx(0.06);
      Antiv0Cuts->SetCutDCADaugToPrimVtx(0.06);
    }
    else if (suffix == "33")
    {
      TrackPosKaonCuts->SetEtaRange(-KaonEtaLow, KaonEtaLow);
      TrackNegKaonCuts->SetEtaRange(-KaonEtaLow, KaonEtaLow);

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
    }
    else if (suffix == "34")
    {
      TrackPosKaonCuts->SetEtaRange(-KaonEtaLow, KaonEtaLow);
      TrackNegKaonCuts->SetEtaRange(-KaonEtaLow, KaonEtaLow);

      TrackPosKaonCuts->SetNClsTPC(KaonNClsLow);
      TrackNegKaonCuts->SetNClsTPC(KaonNClsLow);

      Posv0Daug->SetNClsTPC(80);
      Negv0Daug->SetNClsTPC(80);
      PosAntiv0Daug->SetNClsTPC(80);
      NegAntiv0Daug->SetNClsTPC(80);

      Posv0Daug->SetEtaRange(-0.77, 0.77);
      Negv0Daug->SetEtaRange(-0.77, 0.77);
      PosAntiv0Daug->SetEtaRange(-0.77, 0.77);
      NegAntiv0Daug->SetEtaRange(-0.77, 0.77);
    }
    else if (suffix == "35")
    {
      TrackPosKaonCuts->SetPtRange(KaonPtlow, KaonPtMax);
      TrackNegKaonCuts->SetPtRange(KaonPtlow, KaonPtMax);

      TrackPosKaonCuts->SetPID(AliPID::kKaon, 0.4, KaonNsigmaLow);
      TrackNegKaonCuts->SetPID(AliPID::kKaon, 0.4, KaonNsigmaLow);

      TrackPosKaonCuts->SetNClsTPC(KaonNClsLow);
      TrackNegKaonCuts->SetNClsTPC(KaonNClsLow);

      v0Cuts->SetCutCPA(0.995);
      Antiv0Cuts->SetCutCPA(0.995);

      Posv0Daug->SetNClsTPC(80);
      Negv0Daug->SetNClsTPC(80);
      PosAntiv0Daug->SetNClsTPC(80);
      NegAntiv0Daug->SetNClsTPC(80);

      v0Cuts->SetCutDCADaugTov0Vtx(1.2);
      Antiv0Cuts->SetCutDCADaugTov0Vtx(1.2);
    }
    else if (suffix == "36")
    {
      TrackPosKaonCuts->SetPtRange(KaonPtup, KaonPtMax);
      TrackNegKaonCuts->SetPtRange(KaonPtup, KaonPtMax);

      TrackPosKaonCuts->SetEtaRange(-KaonEtaUp, KaonEtaUp);
      TrackNegKaonCuts->SetEtaRange(-KaonEtaUp, KaonEtaUp);

      TrackPosKaonCuts->SetNClsTPC(KaonNClsLow);
      TrackNegKaonCuts->SetNClsTPC(KaonNClsLow);

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
    }
    else if (suffix == "37")
    {
      TrackPosKaonCuts->SetEtaRange(-KaonEtaUp, KaonEtaUp);
      TrackNegKaonCuts->SetEtaRange(-KaonEtaUp, KaonEtaUp);

      TrackPosKaonCuts->SetPID(AliPID::kKaon, 0.4, KaonNsigmaUp);
      TrackNegKaonCuts->SetPID(AliPID::kKaon, 0.4, KaonNsigmaUp);

      TrackPosKaonCuts->SetNClsTPC(KaonNClsUp);
      TrackNegKaonCuts->SetNClsTPC(KaonNClsUp);

      v0Cuts->SetCutCPA(0.995);
      Antiv0Cuts->SetCutCPA(0.995);

      Posv0Daug->SetEtaRange(-0.83, 0.83);
      Negv0Daug->SetEtaRange(-0.83, 0.83);
      PosAntiv0Daug->SetEtaRange(-0.83, 0.83);
      NegAntiv0Daug->SetEtaRange(-0.83, 0.83);

      v0Cuts->SetCutDCADaugTov0Vtx(1.2);
      Antiv0Cuts->SetCutDCADaugTov0Vtx(1.2);
    }
    else if (suffix == "38")
    {
      TrackPosKaonCuts->SetEtaRange(-KaonEtaLow, KaonEtaLow);
      TrackNegKaonCuts->SetEtaRange(-KaonEtaLow, KaonEtaLow);

      TrackPosKaonCuts->SetNClsTPC(KaonNClsUp);
      TrackNegKaonCuts->SetNClsTPC(KaonNClsUp);

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
    }
    else if (suffix == "39")
    {
      TrackPosKaonCuts->SetPtRange(KaonPtlow, KaonPtMax);
      TrackNegKaonCuts->SetPtRange(KaonPtlow, KaonPtMax);

      TrackPosKaonCuts->SetPID(AliPID::kKaon, 0.4, KaonNsigmaUp);
      TrackNegKaonCuts->SetPID(AliPID::kKaon, 0.4, KaonNsigmaUp);

      TrackPosKaonCuts->SetNClsTPC(KaonNClsUp);
      TrackNegKaonCuts->SetNClsTPC(KaonNClsUp);

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
    }
    else if (suffix == "40")
    {
      TrackPosKaonCuts->SetPtRange(KaonPtup, KaonPtMax);
      TrackNegKaonCuts->SetPtRange(KaonPtup, KaonPtMax);

      TrackPosKaonCuts->SetNClsTPC(KaonNClsLow);
      TrackNegKaonCuts->SetNClsTPC(KaonNClsLow);

      Posv0Daug->SetPID(AliPID::kProton, 999.9, 4);
      Negv0Daug->SetPID(AliPID::kPion, 999.9, 4);
      PosAntiv0Daug->SetPID(AliPID::kPion, 999.9, 4);
      NegAntiv0Daug->SetPID(AliPID::kProton, 999.9, 4);

      v0Cuts->SetCutDCADaugTov0Vtx(1.2);
      Antiv0Cuts->SetCutDCADaugTov0Vtx(1.2);
    }
    else if (suffix == "41")
    {
      TrackPosKaonCuts->SetPtRange(KaonPtlow, KaonPtMax);
      TrackNegKaonCuts->SetPtRange(KaonPtlow, KaonPtMax);

      TrackPosKaonCuts->SetEtaRange(-KaonEtaLow, KaonEtaLow);
      TrackNegKaonCuts->SetEtaRange(-KaonEtaLow, KaonEtaLow);

      TrackPosKaonCuts->SetPID(AliPID::kKaon, 0.4, KaonNsigmaUp);
      TrackNegKaonCuts->SetPID(AliPID::kKaon, 0.4, KaonNsigmaUp);

      v0Cuts->SetCutCPA(0.995);
      Antiv0Cuts->SetCutCPA(0.995);

      Posv0Daug->SetEtaRange(-0.77, 0.77);
      Negv0Daug->SetEtaRange(-0.77, 0.77);
      PosAntiv0Daug->SetEtaRange(-0.77, 0.77);
      NegAntiv0Daug->SetEtaRange(-0.77, 0.77);

      v0Cuts->SetCutDCADaugTov0Vtx(1.2);
      Antiv0Cuts->SetCutDCADaugTov0Vtx(1.2);
    }
    else if (suffix == "42")
    {
      TrackPosKaonCuts->SetPtRange(KaonPtup, KaonPtMax);
      TrackNegKaonCuts->SetPtRange(KaonPtup, KaonPtMax);

      TrackPosKaonCuts->SetPID(AliPID::kKaon, 0.4, KaonNsigmaLow);
      TrackNegKaonCuts->SetPID(AliPID::kKaon, 0.4, KaonNsigmaLow);

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
    }
    else if (suffix == "43")
    {
      TrackPosKaonCuts->SetPtRange(KaonPtup, KaonPtMax);
      TrackNegKaonCuts->SetPtRange(KaonPtup, KaonPtMax);

      TrackPosKaonCuts->SetEtaRange(-KaonEtaUp, KaonEtaUp);
      TrackNegKaonCuts->SetEtaRange(-KaonEtaUp, KaonEtaUp);

      TrackPosKaonCuts->SetPID(AliPID::kKaon, 0.4, KaonNsigmaUp);
      TrackNegKaonCuts->SetPID(AliPID::kKaon, 0.4, KaonNsigmaUp);

      Posv0Daug->SetEtaRange(-0.83, 0.83);
      Negv0Daug->SetEtaRange(-0.83, 0.83);
      PosAntiv0Daug->SetEtaRange(-0.83, 0.83);
      NegAntiv0Daug->SetEtaRange(-0.83, 0.83);
    }
    else if (suffix == "44")
    {
      TrackPosKaonCuts->SetEtaRange(-KaonEtaUp, KaonEtaUp);
      TrackNegKaonCuts->SetEtaRange(-KaonEtaUp, KaonEtaUp);

      TrackPosKaonCuts->SetPID(AliPID::kKaon, 0.4, KaonNsigmaUp);
      TrackNegKaonCuts->SetPID(AliPID::kKaon, 0.4, KaonNsigmaUp);

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
    }
  }

  // now we create the task
  AliAnalysisTaskNanoTreeLPhi *task =
      new AliAnalysisTaskNanoTreeLPhi(
          "AliAnalysisTaskNanoTreeLPhi", isMC);
  // THIS IS VERY IMPORTANT ELSE YOU DONT PROCESS ANY EVENTS
  // kINT7 == Minimum bias
  // kHighMultV0 high multiplicity triggered by the V0 detector
  if (CentEst == "kInt7")
  {
    task->SetTrigger(AliVEvent::kINT7);
    task->SelectCollisionCandidates(AliVEvent::kINT7);
    std::cout << "Added kINT7 Trigger \n";
  }
  else if (CentEst == "kHM")
  {
    task->SetTrigger(AliVEvent::kHighMultV0);
    task->SelectCollisionCandidates(AliVEvent::kHighMultV0);
    std::cout << "Added kHighMultV0 Trigger \n";
  }
  else
  {
    std::cout << "============================================================="
                 "========"
              << std::endl;
    std::cout << "============================================================="
                 "========"
              << std::endl;
    std::cout << "Centrality Estimator not set, fix it else your Results will "
                 "be empty!"
              << std::endl;
    std::cout << "============================================================="
                 "========"
              << std::endl;
    std::cout << "============================================================="
                 "========"
              << std::endl;
  }

  // Throw all our settings to the task
  task->SetEventCuts(evtCuts);
  task->SetLambdaCuts(v0Cuts);
  task->SetAntiLambdaCuts(Antiv0Cuts);
  task->SetPhiCuts(TrackCutsPhi);
  task->SetPosKaonCuts(TrackPosKaonCuts);
  task->SetNegKaonCuts(TrackNegKaonCuts);
  task->SetCollectionConfig(config);

  mgr->AddTask(task);

  TString file = AliAnalysisManager::GetCommonFileName();

  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();

  mgr->ConnectInput(task, 0, cinput);

  TString addon = "";
  if (CentEst == "kInt7")
  {
    addon += "MB";
  }
  else if (CentEst == "kHM")
  {
    addon += "HM";
  }

  TString QAName = Form("%sQA%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *coutputQA = mgr->CreateContainer(
      QAName.Data(), TList::Class(), AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), QAName.Data()));
  mgr->ConnectOutput(task, 1, coutputQA);

  TString EvtCutsName = Form("%sEvtCuts%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *coutputEvtCuts = mgr->CreateContainer(
      EvtCutsName.Data(), TList::Class(), AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), EvtCutsName.Data()));
  mgr->ConnectOutput(task, 2, coutputEvtCuts);

  AliAnalysisDataContainer *coutputv0Cuts;
  TString v0CutsName = Form("%sLambdaCuts%s", addon.Data(), suffix.Data());
  coutputv0Cuts = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      v0CutsName.Data(),
      TList::Class(), AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), v0CutsName.Data()));
  mgr->ConnectOutput(task, 3, coutputv0Cuts);

  AliAnalysisDataContainer *coutputAntiv0Cuts;
  TString Antiv0CutsName = Form("%sAntiLambdaCuts%s", addon.Data(), suffix.Data());
  coutputAntiv0Cuts = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      Antiv0CutsName.Data(),
      TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), Antiv0CutsName.Data()));
  mgr->ConnectOutput(task, 4, coutputAntiv0Cuts);

  TString TrackCutsName = Form("%sKaonPlusCuts%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *couputTrkCuts = mgr->CreateContainer(
      TrackCutsName.Data(), TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), TrackCutsName.Data()));
  mgr->ConnectOutput(task, 5, couputTrkCuts);

  TString AntiTrackCutsName = Form("%sKaonMinusCuts%s", addon.Data(),
                                   suffix.Data());
  AliAnalysisDataContainer *coutputAntiTrkCuts = mgr->CreateContainer(
      AntiTrackCutsName.Data(), TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), AntiTrackCutsName.Data()));
  mgr->ConnectOutput(task, 6, coutputAntiTrkCuts);

  AliAnalysisDataContainer *coutputPhiCuts;
  TString PhiCutsName = Form("%sPhiCuts%s", addon.Data(), suffix.Data());
  coutputPhiCuts = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      PhiCutsName.Data(),
      TList::Class(), AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), PhiCutsName.Data()));
  mgr->ConnectOutput(task, 7, coutputPhiCuts);

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

  AliAnalysisDataContainer *coutputResultsTree;
  TString ResultsTreeName = Form("%sResultsTree%s", addon.Data(), suffix.Data());
  coutputResultsTree = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      ResultsTreeName.Data(),
      TTree::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), ResultsTreeName.Data()));
  mgr->ConnectOutput(task, 10, coutputResultsTree);

  return task;
}
