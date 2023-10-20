#ifndef __CINT__
#include <vector>
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskLambdaPion.h"
#include "AliFemtoDreamEventCuts.h"
#include "AliFemtoDreamTrackCuts.h"
#include "AliFemtoDreamv0Cuts.h"
#include "AliFemtoDreamCollConfig.h"
#endif

AliAnalysisTaskSE *AddTaskFemtoLambdaPion(bool isMC = false,
                                              TString CentEst = "kInt7",
                                              int filterBit = 128,
                                              int WhichPionCut = 0,
                                              const char *sTcut = "0",
                                              bool DoAncestors = false,
                                              bool IsSystematics = false,
                                              bool isNewPC = false,
                                              const char *cutVariation = "0")
{
  TString suffix = TString::Format("%s", cutVariation);
  TString sTsuffix = TString::Format("%s", sTcut);

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

  //========= Init subtasks and start analysis ============================
  // Event Cuts 
  AliFemtoDreamEventCuts *evtCuts = AliFemtoDreamEventCuts::StandardCutsRun2();
  evtCuts->CleanUpMult(false, false, false, true);

  if (sTsuffix == "1")
  {
    evtCuts->SetSphericityCuts(0.7,1);
  }

  // Lambda Cuts
  AliFemtoDreamv0Cuts *v0Cuts = AliFemtoDreamv0Cuts::LambdaCuts(isMC, true, true);
  AliFemtoDreamTrackCuts *Posv0Daug = AliFemtoDreamTrackCuts::DecayProtonCuts(isMC, true, false); // PileUpRej, false
  AliFemtoDreamTrackCuts *Negv0Daug = AliFemtoDreamTrackCuts::DecayPionCuts(isMC, true, false);

  v0Cuts->SetPosDaugterTrackCuts(Posv0Daug);
  v0Cuts->SetNegDaugterTrackCuts(Negv0Daug);
  v0Cuts->SetPDGCodePosDaug(2212); // Proton
  v0Cuts->SetPDGCodeNegDaug(211);  // Pion
  v0Cuts->SetPDGCodev0(3122);      // Lambda

  AliFemtoDreamv0Cuts *Antiv0Cuts = AliFemtoDreamv0Cuts::LambdaCuts(isMC, true, true);
  AliFemtoDreamTrackCuts *PosAntiv0Daug = AliFemtoDreamTrackCuts::DecayPionCuts(isMC, true, false);
  PosAntiv0Daug->SetCutCharge(1);
  AliFemtoDreamTrackCuts *NegAntiv0Daug = AliFemtoDreamTrackCuts::DecayProtonCuts(isMC, true, false);
  NegAntiv0Daug->SetCutCharge(-1);

  Antiv0Cuts->SetPosDaugterTrackCuts(PosAntiv0Daug);
  Antiv0Cuts->SetNegDaugterTrackCuts(NegAntiv0Daug);
  Antiv0Cuts->SetPDGCodePosDaug(211);  // Pion
  Antiv0Cuts->SetPDGCodeNegDaug(2212); // Proton
  Antiv0Cuts->SetPDGCodev0(-3122);     // Lambda

  // MARCEL

  AliFemtoDreamTrackCuts *TrackPosPionCuts = NULL;
  AliFemtoDreamTrackCuts *TrackCutsAntiPion = NULL;

  TrackPosPionCuts = AliFemtoDreamTrackCuts::PrimPionCuts(isMC, true, false, false);
  TrackPosPionCuts->SetFilterBit(96);
  TrackPosPionCuts->SetCutCharge(1);
  TrackCutsAntiPion = AliFemtoDreamTrackCuts::PrimPionCuts(isMC, true, false, false);
  TrackCutsAntiPion->SetFilterBit(96);
  TrackCutsAntiPion->SetCutCharge(-1);

  // ORIGINAL

  // AliFemtoDreamTrackCuts *TrackPosPionCuts =
  //     AliFemtoDreamTrackCuts::PrimPionCuts(isMC, true, false, false);
  // TrackPosPionCuts->SetCutCharge(1);  // Positive Pion charge
  // TrackPosPionCuts->SetFilterBit(filterBit);
  // if (WhichPionCut == 0)
  // {
  //   TrackPosPionCuts->SetPIDkd(); // Oton
  // }
  // else if (WhichPionCut == 1)
  // {
  //   TrackPosPionCuts->SetPIDkd(true, true); // Ramona
  // }
// 
  // AliFemtoDreamTrackCuts *TrackCutsAntiPion =
  //     AliFemtoDreamTrackCuts::PrimPionCuts(isMC, true, false, false);
  // TrackCutsAntiPion->SetCutCharge(-1);  // Negative Pion charge
  // TrackCutsAntiPion->SetFilterBit(filterBit);
  // if (WhichPionCut == 0)
  // {
  //   TrackCutsAntiPion->SetPIDkd(); // Oton
  // }
  // else if (WhichPionCut == 1)
  // {
  //   TrackCutsAntiPion->SetPIDkd(true, true); // Ramona
  // }

  if (IsSystematics)
  {
    evtCuts->SetMinimalBooking(true);
    TrackPosPionCuts->SetMinimalBooking(true);
    TrackCutsAntiPion->SetMinimalBooking(true);
    v0Cuts->SetMinimalBooking(true);
    Antiv0Cuts->SetMinimalBooking(true);
  }

  // Now we define stuff we want for our Particle collection
  // Thanks, CINT - will not compile due to an illegal constructor
  // std::vector<int> PDGParticles ={2212,2212,3122,3122,3312,3312};
  // First we need to tell him about the particles we mix, from the
  // PDG code the mass is obtained.
  std::vector<int> PDGParticles;
  PDGParticles.push_back(211);  // 0 Pion Plus
  PDGParticles.push_back(211);  // 1 Pion Minus
  PDGParticles.push_back(3122); // 2 Lambda
  PDGParticles.push_back(3122); // 3 antiLambda

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
  std::vector<float> mTBins = {0.9, 1.15, 1.4, 4.5, 999.};

  // pairs:
  // K+K+               0
  // K+K-               1
  // K+ La              2
  // K+ bar La          3
  // K-K-               4
  // K- La              5
  // K- bar La          6
  // La La              7
  // La La bar          8
  // La bar La bar      9

  const int npairs = 10;
  for (int i = 0; i < npairs; i++)
  {
    NBins.push_back(1500);
    closeRejection.push_back(false);
    kMin.push_back(0.);
    kMax.push_back(6.);
    pairQA.push_back(0);
  }

  if (IsSystematics)
  {
    pairQA[2] = 12;
    pairQA[3] = 12;
    pairQA[5] = 12;
    pairQA[6] = 12;
  }
  else
  {
    pairQA[0] = 11;
    pairQA[1] = 11;
    pairQA[2] = 12;
    pairQA[3] = 12;
    pairQA[4] = 11;
    pairQA[5] = 12;
    pairQA[6] = 12;
    pairQA[7] = 22;
    pairQA[8] = 22;
    pairQA[9] = 22;
    closeRejection[0] = true;
    closeRejection[1] = true;
    closeRejection[4] = true;
  }

  AliFemtoDreamCollConfig *config = new AliFemtoDreamCollConfig("Femto", "Femto");

  config->SetPDGCodes(PDGParticles);
  config->SetZBins(ZVtxBins);
  config->SetMultBins(MultBins);
  config->SetMultBinning(true);
  config->SetClosePairRejection(closeRejection);
  config->SetDeltaEtaMax(0.012);
  config->SetDeltaPhiMax(0.012);
  config->SetExtendedQAPairs(pairQA);
  config->SetNBinsHist(NBins);
  config->SetMinKRel(kMin);
  config->SetMaxKRel(kMax);
  config->SetUseEventMixing(true);
  config->SetMixingDepth(30);
  config->SetMultiplicityEstimator(AliFemtoDreamEvent::kRef08);

  config->SetMultBinning(true);
  config->SetmTBins(mTBins);
  config->SetDomTMultBinning(true);
  config->SetmTBinning(true);

  if (IsSystematics)
  {
    config->SetMinimalBookingME(true);
  }
  else if (!IsSystematics)
  {
    config->SetPtQA(true);
    config->SetMassQA(true);
    config->SetkTBinning(true);
  }

  if (isMC)
  {
    config->SetMomentumResolution(true);
    if (DoAncestors)
    {
      config->SetAncestors(true);
      config->GetDoAncestorsPlots();
    }
  }
  else
  {
    std::cout << "You are trying to request the Momentum Resolution without MC Info; fix it wont work! \n";
  }

  // Variation cuts
  const float PionPtlow = 0.09;
  const float PionPtup = 0.19;
  const float PionEtaLow = 0.75;
  const float PionEtaUp = 0.85;
  const float PionNClsLow = 70;
  const float PionNClsUp = 90;
  const float PionPtMax = 999;

  AliPID::EParticleType aliPIDParticle;
  aliPIDParticle = AliPID::kPion;
  std::map<std::string, float> PionPIDTight;
  std::map<std::string, float> PionPIDLoose;

  PionPIDTight = {
      {"COMB", 2.7},
      {"TPC", 2.7},
      {"EXCLUSION", 3.3},
  }; // for SetPIDkd() when using oton's K selection
  PionPIDLoose = {
      {"COMB", 3.3},
      {"TPC", 3.3},
      {"EXCLUSION", 2.7},
  };

  /// Systematic variations (taken from pφ and ΛΚ)
  if (IsSystematics)
  {
    if (suffix == "1")
    {
      TrackPosPionCuts->SetNClsTPC(PionNClsUp);
      TrackCutsAntiPion->SetNClsTPC(PionNClsUp);

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
      TrackPosPionCuts->SetPtRange(PionPtlow, PionPtMax);
      TrackCutsAntiPion->SetPtRange(PionPtlow, PionPtMax);

      TrackPosPionCuts->SetEtaRange(-PionEtaLow, PionEtaLow);
      TrackCutsAntiPion->SetEtaRange(-PionEtaLow, PionEtaLow);

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
      TrackPosPionCuts->SetEtaRange(-PionEtaLow, PionEtaLow);
      TrackCutsAntiPion->SetEtaRange(-PionEtaLow, PionEtaLow);

      TrackPosPionCuts->SetNClsTPC(PionNClsLow);
      TrackCutsAntiPion->SetNClsTPC(PionNClsLow);

      Posv0Daug->SetNClsTPC(80);
      Negv0Daug->SetNClsTPC(80);
      PosAntiv0Daug->SetNClsTPC(80);
      NegAntiv0Daug->SetNClsTPC(80);
      v0Cuts->SetCutDCADaugToPrimVtx(0.06);
      Antiv0Cuts->SetCutDCADaugToPrimVtx(0.06);
    }
    else if (suffix == "5")
    {
      TrackPosPionCuts->SetPtRange(PionPtlow, PionPtMax);
      TrackCutsAntiPion->SetPtRange(PionPtlow, PionPtMax);

      TrackPosPionCuts->SetEtaRange(-PionEtaLow, PionEtaLow);
      TrackCutsAntiPion->SetEtaRange(-PionEtaLow, PionEtaLow);

      TrackPosPionCuts->SetPIDkd(true, false, PionPIDLoose["COMB"], PionPIDTight["TPC"], PionPIDTight["EXCLUSION"]);
      TrackCutsAntiPion->SetPIDkd(true, false, PionPIDLoose["COMB"], PionPIDTight["TPC"], PionPIDTight["EXCLUSION"]);

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
      TrackPosPionCuts->SetPIDkd(true, false, PionPIDLoose["COMB"], PionPIDLoose["TPC"], PionPIDTight["EXCLUSION"]);
      TrackCutsAntiPion->SetPIDkd(true, false, PionPIDLoose["COMB"], PionPIDLoose["TPC"], PionPIDTight["EXCLUSION"]);

      TrackPosPionCuts->SetNClsTPC(PionNClsUp);
      TrackCutsAntiPion->SetNClsTPC(PionNClsUp);

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
      TrackPosPionCuts->SetEtaRange(-PionEtaUp, PionEtaUp);
      TrackCutsAntiPion->SetEtaRange(-PionEtaUp, PionEtaUp);

      TrackPosPionCuts->SetPIDkd(true, false, PionPIDLoose["COMB"], PionPIDLoose["TPC"], PionPIDLoose["EXCLUSION"]);
      TrackCutsAntiPion->SetPIDkd(true, false, PionPIDLoose["COMB"], PionPIDLoose["TPC"], PionPIDLoose["EXCLUSION"]);

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
      TrackPosPionCuts->SetPtRange(PionPtlow, PionPtMax);
      TrackCutsAntiPion->SetPtRange(PionPtlow, PionPtMax);

      TrackPosPionCuts->SetEtaRange(-PionEtaLow, PionEtaLow);
      TrackCutsAntiPion->SetEtaRange(-PionEtaLow, PionEtaLow);

      TrackPosPionCuts->SetNClsTPC(PionNClsUp);
      TrackCutsAntiPion->SetNClsTPC(PionNClsUp);

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
      TrackPosPionCuts->SetPtRange(PionPtlow, PionPtMax);
      TrackCutsAntiPion->SetPtRange(PionPtlow, PionPtMax);

      TrackPosPionCuts->SetEtaRange(-PionEtaUp, PionEtaUp);
      TrackCutsAntiPion->SetEtaRange(-PionEtaUp, PionEtaUp);

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
      TrackPosPionCuts->SetEtaRange(-PionEtaLow, PionEtaLow);
      TrackCutsAntiPion->SetEtaRange(-PionEtaLow, PionEtaLow);

      TrackPosPionCuts->SetPIDkd(true, false, PionPIDLoose["COMB"], PionPIDTight["TPC"], PionPIDTight["EXCLUSION"]);
      TrackCutsAntiPion->SetPIDkd(true, false, PionPIDLoose["COMB"], PionPIDTight["TPC"], PionPIDTight["EXCLUSION"]);

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
      TrackPosPionCuts->SetPtRange(PionPtlow, PionPtMax);
      TrackCutsAntiPion->SetPtRange(PionPtlow, PionPtMax);

      TrackPosPionCuts->SetEtaRange(-PionEtaUp, PionEtaUp);
      TrackCutsAntiPion->SetEtaRange(-PionEtaUp, PionEtaUp);

      TrackPosPionCuts->SetPIDkd(true, false, PionPIDTight["COMB"], PionPIDTight["TPC"], PionPIDTight["EXCLUSION"]);
      TrackCutsAntiPion->SetPIDkd(true, false, PionPIDTight["COMB"], PionPIDTight["TPC"], PionPIDTight["EXCLUSION"]);

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
      TrackPosPionCuts->SetEtaRange(-PionEtaUp, PionEtaUp);
      TrackCutsAntiPion->SetEtaRange(-PionEtaUp, PionEtaUp);

      TrackPosPionCuts->SetPIDkd(true, false, PionPIDTight["COMB"], PionPIDLoose["TPC"], PionPIDTight["EXCLUSION"]);
      TrackCutsAntiPion->SetPIDkd(true, false, PionPIDTight["COMB"], PionPIDLoose["TPC"], PionPIDTight["EXCLUSION"]);

      TrackPosPionCuts->SetNClsTPC(PionNClsUp);
      TrackCutsAntiPion->SetNClsTPC(PionNClsUp);

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
      TrackPosPionCuts->SetPtRange(PionPtlow, PionPtMax);
      TrackCutsAntiPion->SetPtRange(PionPtlow, PionPtMax);

      TrackPosPionCuts->SetPIDkd(true, false, PionPIDTight["COMB"], PionPIDLoose["TPC"], PionPIDLoose["EXCLUSION"]);
      TrackCutsAntiPion->SetPIDkd(true, false, PionPIDTight["COMB"], PionPIDLoose["TPC"], PionPIDLoose["EXCLUSION"]);

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
      TrackPosPionCuts->SetPtRange(PionPtup, PionPtMax);
      TrackCutsAntiPion->SetPtRange(PionPtup, PionPtMax);

      TrackPosPionCuts->SetEtaRange(-PionEtaUp, PionEtaUp);
      TrackCutsAntiPion->SetEtaRange(-PionEtaUp, PionEtaUp);

      TrackPosPionCuts->SetPIDkd(true, false, PionPIDTight["COMB"], PionPIDTight["TPC"], PionPIDLoose["EXCLUSION"]);
      TrackCutsAntiPion->SetPIDkd(true, false, PionPIDTight["COMB"], PionPIDTight["TPC"], PionPIDLoose["EXCLUSION"]);

      TrackPosPionCuts->SetNClsTPC(PionNClsUp);
      TrackCutsAntiPion->SetNClsTPC(PionNClsUp);

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
      TrackPosPionCuts->SetPtRange(PionPtlow, PionPtMax);
      TrackCutsAntiPion->SetPtRange(PionPtlow, PionPtMax);

      TrackPosPionCuts->SetPIDkd(true, false, PionPIDLoose["COMB"], PionPIDTight["TPC"], PionPIDLoose["EXCLUSION"]);
      TrackCutsAntiPion->SetPIDkd(true, false, PionPIDLoose["COMB"], PionPIDTight["TPC"], PionPIDLoose["EXCLUSION"]);

      TrackPosPionCuts->SetNClsTPC(PionNClsLow);
      TrackCutsAntiPion->SetNClsTPC(PionNClsLow);

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
      TrackPosPionCuts->SetPtRange(PionPtlow, PionPtMax);
      TrackCutsAntiPion->SetPtRange(PionPtlow, PionPtMax);

      TrackPosPionCuts->SetEtaRange(-PionEtaUp, PionEtaUp);
      TrackCutsAntiPion->SetEtaRange(-PionEtaUp, PionEtaUp);

      TrackPosPionCuts->SetPIDkd(true, false, PionPIDLoose["COMB"], PionPIDTight["TPC"], PionPIDTight["EXCLUSION"]);
      TrackCutsAntiPion->SetPIDkd(true, false, PionPIDLoose["COMB"], PionPIDTight["TPC"], PionPIDTight["EXCLUSION"]);

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
      TrackPosPionCuts->SetPtRange(PionPtlow, PionPtMax);
      TrackCutsAntiPion->SetPtRange(PionPtlow, PionPtMax);

      TrackPosPionCuts->SetEtaRange(-PionEtaUp, PionEtaUp);
      TrackCutsAntiPion->SetEtaRange(-PionEtaUp, PionEtaUp);

      TrackPosPionCuts->SetPIDkd(true, false, PionPIDTight["COMB"], PionPIDTight["TPC"], PionPIDLoose["EXCLUSION"]);
      TrackCutsAntiPion->SetPIDkd(true, false, PionPIDTight["COMB"], PionPIDTight["TPC"], PionPIDLoose["EXCLUSION"]);

      TrackPosPionCuts->SetNClsTPC(PionNClsUp);
      TrackCutsAntiPion->SetNClsTPC(PionNClsUp);

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
      TrackPosPionCuts->SetPIDkd(true, false, PionPIDTight["COMB"], PionPIDTight["TPC"], PionPIDTight["EXCLUSION"]);
      TrackCutsAntiPion->SetPIDkd(true, false, PionPIDTight["COMB"], PionPIDTight["TPC"], PionPIDTight["EXCLUSION"]);

      TrackPosPionCuts->SetNClsTPC(PionNClsLow);
      TrackCutsAntiPion->SetNClsTPC(PionNClsLow);

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
      TrackPosPionCuts->SetPtRange(PionPtlow, PionPtMax);
      TrackCutsAntiPion->SetPtRange(PionPtlow, PionPtMax);

      TrackPosPionCuts->SetNClsTPC(PionNClsUp);
      TrackCutsAntiPion->SetNClsTPC(PionNClsUp);

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

      TrackPosPionCuts->SetEtaRange(-PionEtaUp, PionEtaUp);
      TrackCutsAntiPion->SetEtaRange(-PionEtaUp, PionEtaUp);

      TrackPosPionCuts->SetPIDkd(true, false, PionPIDLoose["COMB"], PionPIDTight["TPC"], PionPIDTight["EXCLUSION"]);
      TrackCutsAntiPion->SetPIDkd(true, false, PionPIDLoose["COMB"], PionPIDTight["TPC"], PionPIDTight["EXCLUSION"]);

      TrackPosPionCuts->SetNClsTPC(PionNClsLow);
      TrackCutsAntiPion->SetNClsTPC(PionNClsLow);

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
      TrackPosPionCuts->SetPtRange(PionPtlow, PionPtMax);
      TrackCutsAntiPion->SetPtRange(PionPtlow, PionPtMax);

      TrackPosPionCuts->SetEtaRange(-PionEtaLow, PionEtaLow);
      TrackCutsAntiPion->SetEtaRange(-PionEtaLow, PionEtaLow);

      TrackPosPionCuts->SetNClsTPC(PionNClsLow);
      TrackCutsAntiPion->SetNClsTPC(PionNClsLow);

      v0Cuts->SetCutCPA(0.995);
      Antiv0Cuts->SetCutCPA(0.995);

      Posv0Daug->SetPID(AliPID::kProton, 999.9, 4);
      Negv0Daug->SetPID(AliPID::kPion, 999.9, 4);
      PosAntiv0Daug->SetPID(AliPID::kPion, 999.9, 4);
      NegAntiv0Daug->SetPID(AliPID::kProton, 999.9, 4);
    }
    else if (suffix == "23")
    {
      TrackPosPionCuts->SetEtaRange(-PionEtaLow, PionEtaLow);
      TrackCutsAntiPion->SetEtaRange(-PionEtaLow, PionEtaLow);

      TrackPosPionCuts->SetPIDkd(true, false, PionPIDLoose["COMB"], PionPIDLoose["TPC"], PionPIDTight["EXCLUSION"]);
      TrackCutsAntiPion->SetPIDkd(true, false, PionPIDLoose["COMB"], PionPIDLoose["TPC"], PionPIDTight["EXCLUSION"]);

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
      TrackPosPionCuts->SetEtaRange(-PionEtaLow, PionEtaLow);
      TrackCutsAntiPion->SetEtaRange(-PionEtaLow, PionEtaLow);

      TrackPosPionCuts->SetPIDkd(true, false, PionPIDLoose["COMB"], PionPIDLoose["TPC"], PionPIDLoose["EXCLUSION"]);
      TrackCutsAntiPion->SetPIDkd(true, false, PionPIDLoose["COMB"], PionPIDLoose["TPC"], PionPIDLoose["EXCLUSION"]);

      TrackPosPionCuts->SetNClsTPC(PionNClsUp);
      TrackCutsAntiPion->SetNClsTPC(PionNClsUp);

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

      // XI
    }
    else if (suffix == "25")
    {
      TrackPosPionCuts->SetPtRange(PionPtlow, PionPtMax);
      TrackCutsAntiPion->SetPtRange(PionPtlow, PionPtMax);

      TrackPosPionCuts->SetEtaRange(-0.83, 0.83);
      TrackCutsAntiPion->SetEtaRange(-0.83, 0.83);

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
      TrackPosPionCuts->SetPtRange(PionPtup, PionPtMax);
      TrackCutsAntiPion->SetPtRange(PionPtup, PionPtMax);

      TrackPosPionCuts->SetPIDkd(true, false, PionPIDTight["COMB"], PionPIDLoose["TPC"], PionPIDLoose["EXCLUSION"]);
      TrackCutsAntiPion->SetPIDkd(true, false, PionPIDTight["COMB"], PionPIDLoose["TPC"], PionPIDLoose["EXCLUSION"]);

      Posv0Daug->SetEtaRange(-0.77, 0.77);
      Negv0Daug->SetEtaRange(-0.77, 0.77);
      PosAntiv0Daug->SetEtaRange(-0.77, 0.77);
      NegAntiv0Daug->SetEtaRange(-0.77, 0.77);
    }
    else if (suffix == "27")
    {
      TrackPosPionCuts->SetEtaRange(-PionEtaLow, PionEtaLow);
      TrackCutsAntiPion->SetEtaRange(-PionEtaLow, PionEtaLow);

      TrackPosPionCuts->SetPIDkd(true, false, PionPIDTight["COMB"], PionPIDTight["TPC"], PionPIDLoose["EXCLUSION"]);
      TrackCutsAntiPion->SetPIDkd(true, false, PionPIDTight["COMB"], PionPIDTight["TPC"], PionPIDLoose["EXCLUSION"]);

      TrackPosPionCuts->SetNClsTPC(PionNClsUp);
      TrackCutsAntiPion->SetNClsTPC(PionNClsUp);

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
      TrackPosPionCuts->SetNClsTPC(PionNClsUp);
      TrackCutsAntiPion->SetNClsTPC(PionNClsUp);

      v0Cuts->SetCutCPA(0.995);
      Antiv0Cuts->SetCutCPA(0.995);

      v0Cuts->SetCutDCADaugToPrimVtx(0.06);
      Antiv0Cuts->SetCutDCADaugToPrimVtx(0.06);
    }
    else if (suffix == "29")
    {
      TrackPosPionCuts->SetPtRange(PionPtlow, PionPtMax);
      TrackCutsAntiPion->SetPtRange(PionPtlow, PionPtMax);

      TrackPosPionCuts->SetEtaRange(-PionEtaLow, PionEtaLow);
      TrackCutsAntiPion->SetEtaRange(-PionEtaLow, PionEtaLow);

      TrackPosPionCuts->SetPIDkd(true, false, PionPIDTight["COMB"], PionPIDLoose["TPC"], PionPIDTight["EXCLUSION"]);
      TrackCutsAntiPion->SetPIDkd(true, false, PionPIDTight["COMB"], PionPIDLoose["TPC"], PionPIDTight["EXCLUSION"]);

      Posv0Daug->SetNClsTPC(80);
      Negv0Daug->SetNClsTPC(80);
      PosAntiv0Daug->SetNClsTPC(80);
      NegAntiv0Daug->SetNClsTPC(80);

      v0Cuts->SetCutDCADaugToPrimVtx(0.06);
      Antiv0Cuts->SetCutDCADaugToPrimVtx(0.06);
    }
    else if (suffix == "30")
    {
      TrackPosPionCuts->SetEtaRange(-PionEtaUp, PionEtaUp);
      TrackCutsAntiPion->SetEtaRange(-PionEtaUp, PionEtaUp);

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
      TrackPosPionCuts->SetPtRange(PionPtlow, PionPtMax);
      TrackCutsAntiPion->SetPtRange(PionPtlow, PionPtMax);

      TrackPosPionCuts->SetPIDkd(true, false, PionPIDTight["COMB"], PionPIDTight["TPC"], PionPIDLoose["EXCLUSION"]);
      TrackCutsAntiPion->SetPIDkd(true, false, PionPIDTight["COMB"], PionPIDTight["TPC"], PionPIDLoose["EXCLUSION"]);

      TrackPosPionCuts->SetNClsTPC(PionNClsLow);
      TrackCutsAntiPion->SetNClsTPC(PionNClsLow);

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

      TrackPosPionCuts->SetEtaRange(-PionEtaLow, PionEtaLow);
      TrackCutsAntiPion->SetEtaRange(-PionEtaLow, PionEtaLow);

      TrackPosPionCuts->SetPIDkd(true, false, PionPIDLoose["COMB"], PionPIDTight["TPC"], PionPIDLoose["EXCLUSION"]);
      TrackCutsAntiPion->SetPIDkd(true, false, PionPIDLoose["COMB"], PionPIDTight["TPC"], PionPIDLoose["EXCLUSION"]);

      TrackPosPionCuts->SetNClsTPC(PionNClsUp);
      TrackCutsAntiPion->SetNClsTPC(PionNClsUp);

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
      TrackPosPionCuts->SetEtaRange(-PionEtaLow, PionEtaLow);
      TrackCutsAntiPion->SetEtaRange(-PionEtaLow, PionEtaLow);

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
      TrackPosPionCuts->SetEtaRange(-PionEtaLow, PionEtaLow);
      TrackCutsAntiPion->SetEtaRange(-PionEtaLow, PionEtaLow);

      TrackPosPionCuts->SetNClsTPC(PionNClsLow);
      TrackCutsAntiPion->SetNClsTPC(PionNClsLow);

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
      TrackPosPionCuts->SetPtRange(PionPtlow, PionPtMax);
      TrackCutsAntiPion->SetPtRange(PionPtlow, PionPtMax);

      TrackPosPionCuts->SetPIDkd(true, false, PionPIDLoose["COMB"], PionPIDTight["TPC"], PionPIDTight["EXCLUSION"]);
      TrackCutsAntiPion->SetPIDkd(true, false, PionPIDLoose["COMB"], PionPIDTight["TPC"], PionPIDTight["EXCLUSION"]);

      TrackPosPionCuts->SetNClsTPC(PionNClsLow);
      TrackCutsAntiPion->SetNClsTPC(PionNClsLow);

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
      TrackPosPionCuts->SetPtRange(PionPtup, PionPtMax);
      TrackCutsAntiPion->SetPtRange(PionPtup, PionPtMax);

      TrackPosPionCuts->SetEtaRange(-PionEtaUp, PionEtaUp);
      TrackCutsAntiPion->SetEtaRange(-PionEtaUp, PionEtaUp);

      TrackPosPionCuts->SetNClsTPC(PionNClsLow);
      TrackCutsAntiPion->SetNClsTPC(PionNClsLow);

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
      TrackPosPionCuts->SetEtaRange(-PionEtaUp, PionEtaUp);
      TrackCutsAntiPion->SetEtaRange(-PionEtaUp, PionEtaUp);

      TrackPosPionCuts->SetPIDkd(true, false, PionPIDTight["COMB"], PionPIDTight["TPC"], PionPIDTight["EXCLUSION"]);
      TrackCutsAntiPion->SetPIDkd(true, false, PionPIDTight["COMB"], PionPIDTight["TPC"], PionPIDTight["EXCLUSION"]);

      TrackPosPionCuts->SetNClsTPC(PionNClsUp);
      TrackCutsAntiPion->SetNClsTPC(PionNClsUp);

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
      TrackPosPionCuts->SetEtaRange(-PionEtaLow, PionEtaLow);
      TrackCutsAntiPion->SetEtaRange(-PionEtaLow, PionEtaLow);

      TrackPosPionCuts->SetNClsTPC(PionNClsUp);
      TrackCutsAntiPion->SetNClsTPC(PionNClsUp);

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
      TrackPosPionCuts->SetPtRange(PionPtlow, PionPtMax);
      TrackCutsAntiPion->SetPtRange(PionPtlow, PionPtMax);

      TrackPosPionCuts->SetPIDkd(true, false, PionPIDTight["COMB"], PionPIDLoose["TPC"], PionPIDTight["EXCLUSION"]);
      TrackCutsAntiPion->SetPIDkd(true, false, PionPIDTight["COMB"], PionPIDLoose["TPC"], PionPIDTight["EXCLUSION"]);

      TrackPosPionCuts->SetNClsTPC(PionNClsUp);
      TrackCutsAntiPion->SetNClsTPC(PionNClsUp);

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
      TrackPosPionCuts->SetPtRange(PionPtup, PionPtMax);
      TrackCutsAntiPion->SetPtRange(PionPtup, PionPtMax);

      TrackPosPionCuts->SetNClsTPC(PionNClsLow);
      TrackCutsAntiPion->SetNClsTPC(PionNClsLow);

      Posv0Daug->SetPID(AliPID::kProton, 999.9, 4);
      Negv0Daug->SetPID(AliPID::kPion, 999.9, 4);
      PosAntiv0Daug->SetPID(AliPID::kPion, 999.9, 4);
      NegAntiv0Daug->SetPID(AliPID::kProton, 999.9, 4);

      v0Cuts->SetCutDCADaugTov0Vtx(1.2);
      Antiv0Cuts->SetCutDCADaugTov0Vtx(1.2);
    }
    else if (suffix == "41")
    {
      TrackPosPionCuts->SetPtRange(PionPtlow, PionPtMax);
      TrackCutsAntiPion->SetPtRange(PionPtlow, PionPtMax);

      TrackPosPionCuts->SetEtaRange(-PionEtaLow, PionEtaLow);
      TrackCutsAntiPion->SetEtaRange(-PionEtaLow, PionEtaLow);

      TrackPosPionCuts->SetPIDkd(true, false, PionPIDLoose["COMB"], PionPIDTight["TPC"], PionPIDTight["EXCLUSION"]);
      TrackCutsAntiPion->SetPIDkd(true, false, PionPIDLoose["COMB"], PionPIDTight["TPC"], PionPIDTight["EXCLUSION"]);

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
      TrackPosPionCuts->SetPtRange(PionPtup, PionPtMax);
      TrackCutsAntiPion->SetPtRange(PionPtup, PionPtMax);

      TrackPosPionCuts->SetPIDkd(true, false, PionPIDLoose["COMB"], PionPIDLoose["TPC"], PionPIDTight["EXCLUSION"]);
      TrackCutsAntiPion->SetPIDkd(true, false, PionPIDLoose["COMB"], PionPIDLoose["TPC"], PionPIDTight["EXCLUSION"]);

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
      TrackPosPionCuts->SetPtRange(PionPtup, PionPtMax);
      TrackCutsAntiPion->SetPtRange(PionPtup, PionPtMax);

      TrackPosPionCuts->SetEtaRange(-PionEtaUp, PionEtaUp);
      TrackCutsAntiPion->SetEtaRange(-PionEtaUp, PionEtaUp);

      TrackPosPionCuts->SetPIDkd(true, false, PionPIDLoose["COMB"], PionPIDLoose["TPC"], PionPIDLoose["EXCLUSION"]);
      TrackCutsAntiPion->SetPIDkd(true, false, PionPIDLoose["COMB"], PionPIDLoose["TPC"], PionPIDLoose["EXCLUSION"]);

      Posv0Daug->SetEtaRange(-0.83, 0.83);
      Negv0Daug->SetEtaRange(-0.83, 0.83);
      PosAntiv0Daug->SetEtaRange(-0.83, 0.83);
      NegAntiv0Daug->SetEtaRange(-0.83, 0.83);
    }
    else if (suffix == "44")
    {
      TrackPosPionCuts->SetEtaRange(-PionEtaUp, PionEtaUp);
      TrackCutsAntiPion->SetEtaRange(-PionEtaUp, PionEtaUp);

      TrackPosPionCuts->SetPIDkd(true, false, PionPIDLoose["COMB"], PionPIDTight["TPC"], PionPIDLoose["EXCLUSION"]);
      TrackCutsAntiPion->SetPIDkd(true, false, PionPIDLoose["COMB"], PionPIDTight["TPC"], PionPIDLoose["EXCLUSION"]);

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
  AliAnalysisTaskLambdaPion *task =
      new AliAnalysisTaskLambdaPion(
          "AliAnalysisTaskLambdaPion", isMC, isNewPC);

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
  task->SetPosPionCuts(TrackPosPionCuts);
  task->SetNegPionCuts(TrackCutsAntiPion);
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

  TString TrackCutsName = Form("%sPionPlusCuts%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *couputTrkCuts = mgr->CreateContainer(
      TrackCutsName.Data(), TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), TrackCutsName.Data()));
  mgr->ConnectOutput(task, 5, couputTrkCuts);

  TString AntiTrackCutsName = Form("%sPionMinusCuts%s", addon.Data(),
                                   suffix.Data());
  AliAnalysisDataContainer *coutputAntiTrkCuts = mgr->CreateContainer(
      AntiTrackCutsName.Data(), TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), AntiTrackCutsName.Data()));
  mgr->ConnectOutput(task, 6, coutputAntiTrkCuts);

  AliAnalysisDataContainer *coutputResults;
  TString ResultsName = Form("%sResults%s", addon.Data(), suffix.Data());
  coutputResults = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      ResultsName.Data(),
      TList::Class(), AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), ResultsName.Data()));
  mgr->ConnectOutput(task, 7, coutputResults);

  AliAnalysisDataContainer *coutputResultsQA;
  TString ResultsQAName = Form("%sResultsQA%s", addon.Data(), suffix.Data());
  coutputResultsQA = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      ResultsQAName.Data(),
      TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), ResultsQAName.Data()));
  mgr->ConnectOutput(task, 8, coutputResultsQA);

  if (isMC)
  {
    AliAnalysisDataContainer *coutputv0CutsMC;
    TString v0CutsMCName = Form("%sv0CutsMC%s", addon.Data(), suffix.Data());
    coutputv0CutsMC = mgr->CreateContainer(
        //@suppress("Invalid arguments") it works ffs
        v0CutsMCName.Data(),
        TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), v0CutsMCName.Data()));
    mgr->ConnectOutput(task, 9, coutputv0CutsMC);

    AliAnalysisDataContainer *coutputAntiv0CutsMC;
    TString Antiv0CutsMCName = Form("%sAntiv0CutsMC%s", addon.Data(), suffix.Data());
    coutputAntiv0CutsMC = mgr->CreateContainer(
        //@suppress("Invalid arguments") it works ffs
        Antiv0CutsMCName.Data(),
        TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), Antiv0CutsMCName.Data()));
    mgr->ConnectOutput(task, 10, coutputAntiv0CutsMC);

    AliAnalysisDataContainer *coutputTrkCutsMC;
    TString TrkCutsMCName = Form("%sTrkCutsMC%s", addon.Data(), suffix.Data());
    coutputTrkCutsMC = mgr->CreateContainer(
        //@suppress("Invalid arguments") it works ffs
        TrkCutsMCName.Data(),
        TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), TrkCutsMCName.Data()));
    mgr->ConnectOutput(task, 11, coutputTrkCutsMC);

    AliAnalysisDataContainer *coutputAntiTrkCutsMC;
    TString AntiTrkCutsMCName = Form("%sAntiTrkCutsMC%s", addon.Data(), suffix.Data());
    coutputAntiTrkCutsMC = mgr->CreateContainer(
        //@suppress("Invalid arguments") it works ffs
        AntiTrkCutsMCName.Data(),
        TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), AntiTrkCutsMCName.Data()));
    mgr->ConnectOutput(task, 12, coutputAntiTrkCutsMC);
  }

  return task;
}
