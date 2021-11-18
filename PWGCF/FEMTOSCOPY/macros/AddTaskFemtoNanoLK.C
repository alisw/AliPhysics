#ifndef __CINT__
#include <vector>
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskNanoLK.h"
#include "AliFemtoDreamEventCuts.h"
#include "AliFemtoDreamTrackCuts.h"
#include "AliFemtoDreamCascadeCuts.h"
#include "AliFemtoDreamCollConfig.h"
#endif

AliAnalysisTaskSE *AddTaskFemtoNanoLK(bool fullBlastQA = false,           //1
                                           bool isMC = false,             //2
                                           int fFilterBit = 128,          //3
                                           TString triggerData = "kInt7", //4
                                           bool DodPhidEtaPlots = false,  //5
                                           int WhichKaonCut = 0,          //6
                                           const char *cutVariation = "0")
{

  TString suffix = TString::Format("%s", cutVariation);

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    Error("AddTaskFemtoNano_pLambda()", "No analysis manager found.");
    return 0x0;
  }

  // ================== GetInputEventHandler =============================
  AliVEventHandler *inputHandler = mgr->GetInputEventHandler();
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();

  //========= Init subtasks and start analysis ============================
  // Event Cuts
  AliFemtoDreamEventCuts *evtCuts = AliFemtoDreamEventCuts::StandardCutsRun2();
  evtCuts->CleanUpMult(false, false, false, true);

  /// Kaon cuts
  AliFemtoDreamTrackCuts *TrackPosKaonCuts =
      AliFemtoDreamTrackCuts::PrimKaonCuts(isMC, true, false, false);
  TrackPosKaonCuts->SetCutCharge(1);
  TrackPosKaonCuts->SetFilterBit(fFilterBit);
  if(WhichKaonCut==0){
    TrackPosKaonCuts->SetPIDkd();// Oton
  } else if (WhichKaonCut==1){
    TrackPosKaonCuts->SetPIDkd(true,true); //Ramona
  }

  AliFemtoDreamTrackCuts *TrackNegKaonCuts =
      AliFemtoDreamTrackCuts::PrimKaonCuts(isMC, true, false, false);
  TrackNegKaonCuts->SetCutCharge(-1);
  TrackNegKaonCuts->SetFilterBit(fFilterBit);
  if (WhichKaonCut == 0)
  {
    TrackNegKaonCuts->SetPIDkd(); // Oton
  }
  else if (WhichKaonCut == 1)
  {
    TrackNegKaonCuts->SetPIDkd(true, true); //Ramona
  }

  // TrackPosKaonCuts->SetDCAVtxZ(0.4); ///Emma´s cut optimized for φ reco
  // TrackNegKaonCuts->SetDCAVtxZ(0.4);
  // TrackPosKaonCuts->SetDCAVtxXY(0.8);
  // TrackNegKaonCuts->SetDCAVtxXY(0.8);

  /// Lambda cuts
  AliFemtoDreamv0Cuts *v0Cuts = AliFemtoDreamv0Cuts::LambdaCuts(isMC, true, true);
  AliFemtoDreamTrackCuts *Posv0Daug = AliFemtoDreamTrackCuts::DecayProtonCuts(isMC, true, false); //PileUpRej, false
  AliFemtoDreamTrackCuts *Negv0Daug = AliFemtoDreamTrackCuts::DecayPionCuts(isMC, true, false);
  v0Cuts->SetPosDaugterTrackCuts(Posv0Daug);
  v0Cuts->SetNegDaugterTrackCuts(Negv0Daug);
  Posv0Daug->SetCutCharge(1);
  Negv0Daug->SetCutCharge(-1);
  v0Cuts->SetPDGCodePosDaug(2212); //Proton
  v0Cuts->SetPDGCodeNegDaug(211);  //Pion
  v0Cuts->SetPDGCodev0(3122);      //Lambda

  AliFemtoDreamv0Cuts *Antiv0Cuts = AliFemtoDreamv0Cuts::LambdaCuts(isMC, true, true);
  AliFemtoDreamTrackCuts *PosAntiv0Daug = AliFemtoDreamTrackCuts::DecayPionCuts(isMC, true, false);
  AliFemtoDreamTrackCuts *NegAntiv0Daug = AliFemtoDreamTrackCuts::DecayProtonCuts(isMC, true, false);
  Antiv0Cuts->SetPosDaugterTrackCuts(PosAntiv0Daug);
  Antiv0Cuts->SetNegDaugterTrackCuts(NegAntiv0Daug);
  PosAntiv0Daug->SetCutCharge(1);
  NegAntiv0Daug->SetCutCharge(-1);
  Antiv0Cuts->SetPDGCodePosDaug(211);  //Pion
  Antiv0Cuts->SetPDGCodeNegDaug(2212); //Proton
  Antiv0Cuts->SetPDGCodev0(-3122);     //Lambda

  if (!fullBlastQA || suffix != "0")
  {
    evtCuts->SetMinimalBooking(true);
    TrackPosKaonCuts->SetMinimalBooking(true);
    TrackNegKaonCuts->SetMinimalBooking(true);
    v0Cuts->SetMinimalBooking(true);
    Antiv0Cuts->SetMinimalBooking(true);
  }

  // First we need to tell him about the particles we mix, from the
  // PDG code the mass is obtained.
  AliFemtoDreamCollConfig *config = new AliFemtoDreamCollConfig("Femto",
                                                                "Femto", false);
  // Femto Collection
  std::vector<int> PDGParticles;
  PDGParticles.push_back(321); //K+
  PDGParticles.push_back(321); //K-
  PDGParticles.push_back(3122); //Lambda
  PDGParticles.push_back(3122); //AntiLambda


  std::vector<int> NBins;
  std::vector<float> kMin;
  std::vector<float> kMax;
  std::vector<int> pairQA;
  std::vector<int> pairQASyst;
  std::vector<bool> closeRejection;
  //pairs:
  //K+K+               0
  //K+K-               1
  //K+ La              2
  //K+ bar La          3
  //K-K-               4
  //K- La              5
  //K- bar La          6
  //La La              7
  //La La bar          8
  //La bar La bar      9

  const int nPairs = 10;
  for (int i = 0; i < nPairs; ++i)
  {
    pairQA.push_back(0);
    closeRejection.push_back(false);
    NBins.push_back(1500);
    kMin.push_back(0.);
    kMax.push_back(6.);
  }
  if (suffix != "0")
  {
    pairQA[2] = 12;
    pairQA[3] = 12;
    pairQA[5] = 12;
    pairQA[6] = 12;
    pairQA[7] = 22;
    pairQA[8] = 22;
    pairQA[9] = 22;
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
    closeRejection[0] = true; // K+K+
    closeRejection[4] = true; // K-K-
  }

  config->SetPDGCodes(PDGParticles);
  config->SetNBinsHist(NBins);
  config->SetMinKRel(kMin);
  config->SetMaxKRel(kMax);
  config->SetClosePairRejection(closeRejection);
  config->SetDeltaEtaMax(0.012);
  config->SetDeltaPhiMax(0.012);
  config->SetExtendedQAPairs(pairQA);
  config->SetMixingDepth(10);
  config->SetUseEventMixing(true);
  config->SetMultiplicityEstimator(AliFemtoDreamEvent::kRef08);

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

  config->SetdPhidEtaPlots(DodPhidEtaPlots);

  if (fullBlastQA)
  {
    config->SetkTBinning(true);
    config->SetPtQA(true);
    config->SetMultBinning(true);
    config->SetmTBinning(true);
  }

  if (!fullBlastQA || suffix != "0")
  {
    config->SetMinimalBookingME(true);
    config->SetMinimalBookingSample(true);
    config->SetMultBinning(true);
    config->SetmTBinning(true);
  }

  if (isMC)
  {
    config->SetMomentumResolution(true); //kstar true vs. kstar reco
  }
  else
  {
    std::cout
        << "You are trying to request the Momentum Resolution without MC Info; fix it wont work! \n";
  }

///Setting all cuts for systematics

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

  AliAnalysisTaskNanoLK *task = new AliAnalysisTaskNanoLK("femtoLK", isMC);

  if (!fullBlastQA)
  {
    task->SetRunTaskLightWeight(true);
  }
  if (triggerData == "kINT7")
  {
    task->SelectCollisionCandidates(AliVEvent::kINT7);
  }
  else if (triggerData == "kHM")
  {
    task->SelectCollisionCandidates(AliVEvent::kHighMultV0);
  }
  task->SetEventCuts(evtCuts);
  task->SetPosKaonCuts(TrackPosKaonCuts);
  task->SetNegKaonCuts(TrackNegKaonCuts);
  task->Setv0Cuts(v0Cuts);
  task->SetAntiv0Cuts(Antiv0Cuts);
  task->SetCorrelationConfig(config);
  mgr->AddTask(task);

  TString addon = "";
  if (triggerData == "kINT7")
  {
    addon += "MBLK";
  }
  else if (triggerData == "kHM")
  {
    addon += "HMLK";
  }

  TString file = AliAnalysisManager::GetCommonFileName();

  mgr->ConnectInput(task, 0, cinput);

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

  TString TrackCutsName = Form("%sTrackCuts%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *couputTrkCuts = mgr->CreateContainer(
      TrackCutsName.Data(), TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), TrackCutsName.Data()));
  mgr->ConnectOutput(task, 3, couputTrkCuts);

  TString AntiTrackCutsName = Form("%sAntiTrackCuts%s", addon.Data(),
                                   suffix.Data());
  AliAnalysisDataContainer *coutputAntiTrkCuts = mgr->CreateContainer(
      AntiTrackCutsName.Data(), TList::Class(),
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

  AliAnalysisDataContainer *coutputResultsSample;
  TString ResultsSampleName = Form("%sResultsSample%s", addon.Data(),
                                   suffix.Data());
  coutputResultsSample = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      ResultsSampleName.Data(),
      TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), ResultsSampleName.Data()));
  mgr->ConnectOutput(task, 9, coutputResultsSample);

  AliAnalysisDataContainer *coutputResultsSampleQA;
  TString ResultsSampleQAName = Form("%sResultsSampleQA%s", addon.Data(),
                                     suffix.Data());
  coutputResultsSampleQA = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      ResultsSampleQAName.Data(),
      TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), ResultsSampleQAName.Data()));
  mgr->ConnectOutput(task, 10, coutputResultsSampleQA);

  if (isMC)
  {
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

    AliAnalysisDataContainer *coutputv0CutsMC;
    TString v0CutsMCName = Form("%sv0CutsMC%s", addon.Data(), suffix.Data());
    coutputv0CutsMC = mgr->CreateContainer(
        //@suppress("Invalid arguments") it works ffs
        v0CutsMCName.Data(),
        TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), v0CutsMCName.Data()));
    mgr->ConnectOutput(task, 13, coutputv0CutsMC);

    AliAnalysisDataContainer *coutputAntiv0CutsMC;
    TString Antiv0CutsMCName = Form("%sAntiv0CutsMC%s", addon.Data(), suffix.Data());
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