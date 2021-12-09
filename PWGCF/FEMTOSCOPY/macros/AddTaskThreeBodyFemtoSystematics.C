#if !defined(__CINT__) || defined(__CLING__)
#include <vector>
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskThreeBodyFemto.h"
#include "AliAnalysisTaskThreeBodyFemtoAOD.h"
#include "AliFemtoDreamEventCuts.h"
#include "AliFemtoDreamTrackCuts.h"
#include "AliFemtoDreamCascadeCuts.h"
#include "AliFemtoDreamCollConfig.h"
#endif

AliAnalysisTaskSE *AddTaskThreeBodyFemtoSystematics(int trigger = 0, bool fullBlastQA = true,
                                     bool isMC = false, bool isNano = true,
                                     bool plotAdditionalPlots=true,  
                                     int mixingDepthFromTask = 20,
                                     bool ClosePairRejectionForAll = "false", int mixinfChoice = 0,
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

  // Track Cuts
  AliFemtoDreamTrackCuts *TrackCuts = AliFemtoDreamTrackCuts::PrimProtonCuts(
      isMC, true, false, false);
  TrackCuts->SetFilterBit(128);
  TrackCuts->SetCutCharge(1);

  AliFemtoDreamTrackCuts *AntiTrackCuts =
      AliFemtoDreamTrackCuts::PrimProtonCuts(isMC, true, false, false);
  AntiTrackCuts->SetFilterBit(128);
  AntiTrackCuts->SetCutCharge(-1);
  //Lambda Cuts
  AliFemtoDreamv0Cuts *v0Cuts = AliFemtoDreamv0Cuts::LambdaCuts(isMC, true,
                                                                false);
  AliFemtoDreamTrackCuts *Posv0Daug = AliFemtoDreamTrackCuts::DecayProtonCuts(
      isMC, true, false);

  AliFemtoDreamTrackCuts *Negv0Daug = AliFemtoDreamTrackCuts::DecayPionCuts(
      isMC, true, false);
 

  v0Cuts->SetPosDaugterTrackCuts(Posv0Daug);
  v0Cuts->SetNegDaugterTrackCuts(Negv0Daug);
  v0Cuts->SetPDGCodePosDaug(2212);  //Proton
  v0Cuts->SetPDGCodeNegDaug(211);  //Pion
  v0Cuts->SetPDGCodev0(3122);  //Lambda

  AliFemtoDreamv0Cuts *Antiv0Cuts = AliFemtoDreamv0Cuts::LambdaCuts(isMC, true,
                                                                    false);
  AliFemtoDreamTrackCuts *PosAntiv0Daug = AliFemtoDreamTrackCuts::DecayPionCuts(
      isMC, true, false);
  PosAntiv0Daug->SetCutCharge(1);
  AliFemtoDreamTrackCuts *NegAntiv0Daug =
      AliFemtoDreamTrackCuts::DecayProtonCuts(isMC, true, false);
  NegAntiv0Daug->SetCutCharge(-1);

  Antiv0Cuts->SetPosDaugterTrackCuts(PosAntiv0Daug);
  Antiv0Cuts->SetNegDaugterTrackCuts(NegAntiv0Daug);
  Antiv0Cuts->SetPDGCodePosDaug(211);  //Pion
  Antiv0Cuts->SetPDGCodeNegDaug(2212);  //Proton
  Antiv0Cuts->SetPDGCodev0(-3122);  //Lambda

  if (!fullBlastQA) {
    evtCuts->SetMinimalBooking(true);
    TrackCuts->SetMinimalBooking(true);
    AntiTrackCuts->SetMinimalBooking(true);
    v0Cuts->SetMinimalBooking(true);
    Antiv0Cuts->SetMinimalBooking(true);
  }


 


  AliFemtoDreamCollConfig *config = new AliFemtoDreamCollConfig("Femto",
                                                                "Femto", false);
  // Femto Collection
  std::vector<int> PDGParticles;
  PDGParticles.push_back(2212);
  PDGParticles.push_back(2212);
  PDGParticles.push_back(3122);
  PDGParticles.push_back(3122);

  std::vector<int> NBins;
  std::vector<float> kMin;
  std::vector<float> kMax;
  std::vector<int> pairQA;
  std::vector<bool> closeRejection;

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
  pairQA[2] = 12;
  pairQA[6] = 12;

  closeRejection[0] = true;  // pp
  closeRejection[4] = true;  // barp barp

  config->SetPDGCodes(PDGParticles);
  config->SetNBinsHist(NBins);
  config->SetMinKRel(kMin);
  config->SetMaxKRel(kMax);
  config->SetClosePairRejection(closeRejection);
  config->SetDeltaEtaMax(0.017);
  config->SetDeltaPhiMax(0.017);
  config->SetExtendedQAPairs(pairQA);
  config->SetMixingDepth(mixingDepthFromTask);
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

  config->SetMultBinning(true);


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

  } else if (suffix == "2") {
    TrackCuts->SetPtRange(0.6, 4.05);
    AntiTrackCuts->SetPtRange(0.6, 4.05);

    TrackCuts->SetEtaRange(-0.77, 0.77);
    AntiTrackCuts->SetEtaRange(-0.77, 0.77);

    config->SetDeltaEtaMax(0.019);
    config->SetDeltaPhiMax(0.019);

    Posv0Daug->SetPID(AliPID::kProton, 999.9, 4);
    Negv0Daug->SetPID(AliPID::kPion, 999.9, 4);
    PosAntiv0Daug->SetPID(AliPID::kPion, 999.9, 4);
    NegAntiv0Daug->SetPID(AliPID::kProton, 999.9, 4);

  } else if (suffix == "3") {
    v0Cuts->SetCutCPA(0.995);
    Antiv0Cuts->SetCutCPA(0.995);

    Posv0Daug->SetEtaRange(-0.77, 0.77);
    Negv0Daug->SetEtaRange(-0.77, 0.77);
    PosAntiv0Daug->SetEtaRange(-0.77, 0.77);
    NegAntiv0Daug->SetEtaRange(-0.77, 0.77);

    v0Cuts->SetCutDCADaugTov0Vtx(1.2);
    Antiv0Cuts->SetCutDCADaugTov0Vtx(1.2);

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

  } else if (suffix == "5") {
    TrackCuts->SetPtRange(0.4, 4.05);
    AntiTrackCuts->SetPtRange(0.4, 4.05);

    TrackCuts->SetEtaRange(-0.77, 0.77);
    AntiTrackCuts->SetEtaRange(-0.77, 0.77);

    TrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);
    AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);

    config->SetDeltaEtaMax(0.019);
    config->SetDeltaPhiMax(0.019);

    v0Cuts->SetCutCPA(0.995);
    Antiv0Cuts->SetCutCPA(0.995);

    Posv0Daug->SetEtaRange(-0.83, 0.83);
    Negv0Daug->SetEtaRange(-0.83, 0.83);
    PosAntiv0Daug->SetEtaRange(-0.83, 0.83);
    NegAntiv0Daug->SetEtaRange(-0.83, 0.83);

    v0Cuts->SetCutDCADaugTov0Vtx(1.2);
    Antiv0Cuts->SetCutDCADaugTov0Vtx(1.2);

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

  } else if (suffix == "7") {
    TrackCuts->SetEtaRange(-0.85, 0.85);
    AntiTrackCuts->SetEtaRange(-0.85, 0.85);

    TrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);
    AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);

    config->SetDeltaEtaMax(0.019);
    config->SetDeltaPhiMax(0.019);

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

  } else if (suffix == "8") {
    TrackCuts->SetPtRange(0.4, 4.05);
    AntiTrackCuts->SetPtRange(0.4, 4.05);

    TrackCuts->SetEtaRange(-0.77, 0.77);
    AntiTrackCuts->SetEtaRange(-0.77, 0.77);

    TrackCuts->SetNClsTPC(90);
    AntiTrackCuts->SetNClsTPC(90);

    config->SetDeltaEtaMax(0.019);
    config->SetDeltaPhiMax(0.019);

    Posv0Daug->SetEtaRange(-0.83, 0.83);
    Negv0Daug->SetEtaRange(-0.83, 0.83);
    PosAntiv0Daug->SetEtaRange(-0.83, 0.83);
    NegAntiv0Daug->SetEtaRange(-0.83, 0.83);

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

  } else if (suffix == "12") {
    TrackCuts->SetPtRange(0.4, 4.05);
    AntiTrackCuts->SetPtRange(0.4, 4.05);

    TrackCuts->SetEtaRange(-0.85, 0.85);
    AntiTrackCuts->SetEtaRange(-0.85, 0.85);

    TrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);
    AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);

    config->SetDeltaEtaMax(0.019);
    config->SetDeltaPhiMax(0.019);

    Posv0Daug->SetPID(AliPID::kProton, 999.9, 4);
    Negv0Daug->SetPID(AliPID::kPion, 999.9, 4);
    PosAntiv0Daug->SetPID(AliPID::kPion, 999.9, 4);
    NegAntiv0Daug->SetPID(AliPID::kProton, 999.9, 4);

    Posv0Daug->SetEtaRange(-0.77, 0.77);
    Negv0Daug->SetEtaRange(-0.77, 0.77);
    PosAntiv0Daug->SetEtaRange(-0.77, 0.77);
    NegAntiv0Daug->SetEtaRange(-0.77, 0.77);

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

  } else if (suffix == "14") {
    TrackCuts->SetPtRange(0.4, 4.05);
    AntiTrackCuts->SetPtRange(0.4, 4.05);

    TrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);
    AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);

    config->SetDeltaEtaMax(0.019);
    config->SetDeltaPhiMax(0.019);

    Posv0Daug->SetNClsTPC(80);
    Negv0Daug->SetNClsTPC(80);
    PosAntiv0Daug->SetNClsTPC(80);
    NegAntiv0Daug->SetNClsTPC(80);

    Posv0Daug->SetEtaRange(-0.83, 0.83);
    Negv0Daug->SetEtaRange(-0.83, 0.83);
    PosAntiv0Daug->SetEtaRange(-0.83, 0.83);
    NegAntiv0Daug->SetEtaRange(-0.83, 0.83);

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

  } else if (suffix == "16") {
    TrackCuts->SetPtRange(0.4, 4.05);
    AntiTrackCuts->SetPtRange(0.4, 4.05);

    TrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);
    AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);

    TrackCuts->SetNClsTPC(70);
    AntiTrackCuts->SetNClsTPC(70);

    config->SetDeltaEtaMax(0.019);
    config->SetDeltaPhiMax(0.019);

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

  } else if (suffix == "17") {
    TrackCuts->SetPtRange(0.4, 4.05);
    AntiTrackCuts->SetPtRange(0.4, 4.05);

    TrackCuts->SetEtaRange(-0.85, 0.85);
    AntiTrackCuts->SetEtaRange(-0.85, 0.85);

    TrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);
    AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);

    config->SetDeltaEtaMax(0.019);
    config->SetDeltaPhiMax(0.019);

    v0Cuts->SetCutCPA(0.995);
    Antiv0Cuts->SetCutCPA(0.995);

    Posv0Daug->SetEtaRange(-0.77, 0.77);
    Negv0Daug->SetEtaRange(-0.77, 0.77);
    PosAntiv0Daug->SetEtaRange(-0.77, 0.77);
    NegAntiv0Daug->SetEtaRange(-0.77, 0.77);

    v0Cuts->SetCutDCADaugToPrimVtx(0.06);
    Antiv0Cuts->SetCutDCADaugToPrimVtx(0.06);

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

  } else if (suffix == "19") {
    TrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);
    AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);

    TrackCuts->SetNClsTPC(70);
    AntiTrackCuts->SetNClsTPC(70);

    config->SetDeltaEtaMax(0.019);
    config->SetDeltaPhiMax(0.019);

    v0Cuts->SetCutCPA(0.995);
    Antiv0Cuts->SetCutCPA(0.995);

    Posv0Daug->SetNClsTPC(80);
    Negv0Daug->SetNClsTPC(80);
    PosAntiv0Daug->SetNClsTPC(80);
    NegAntiv0Daug->SetNClsTPC(80);

    v0Cuts->SetCutDCADaugTov0Vtx(1.2);
    Antiv0Cuts->SetCutDCADaugTov0Vtx(1.2);

  } else if (suffix == "20") {
    TrackCuts->SetPtRange(0.4, 4.05);
    AntiTrackCuts->SetPtRange(0.4, 4.05);

    TrackCuts->SetNClsTPC(90);
    AntiTrackCuts->SetNClsTPC(90);

    config->SetDeltaEtaMax(0.019);
    config->SetDeltaPhiMax(0.019);

    Posv0Daug->SetNClsTPC(80);
    Negv0Daug->SetNClsTPC(80);
    PosAntiv0Daug->SetNClsTPC(80);
    NegAntiv0Daug->SetNClsTPC(80);

    Posv0Daug->SetEtaRange(-0.83, 0.83);
    Negv0Daug->SetEtaRange(-0.83, 0.83);
    PosAntiv0Daug->SetEtaRange(-0.83, 0.83);
    NegAntiv0Daug->SetEtaRange(-0.83, 0.83);

  } else if (suffix == "21") {

    TrackCuts->SetEtaRange(-0.85, 0.85);
    AntiTrackCuts->SetEtaRange(-0.85, 0.85);

    TrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);
    AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);

    TrackCuts->SetNClsTPC(70);
    AntiTrackCuts->SetNClsTPC(70);

    config->SetDeltaEtaMax(0.019);
    config->SetDeltaPhiMax(0.019);

    Posv0Daug->SetNClsTPC(80);
    Negv0Daug->SetNClsTPC(80);
    PosAntiv0Daug->SetNClsTPC(80);
    NegAntiv0Daug->SetNClsTPC(80);

    Posv0Daug->SetEtaRange(-0.77, 0.77);
    Negv0Daug->SetEtaRange(-0.77, 0.77);
    PosAntiv0Daug->SetEtaRange(-0.77, 0.77);
    NegAntiv0Daug->SetEtaRange(-0.77, 0.77);

  } else if (suffix == "22") {
    TrackCuts->SetPtRange(0.4, 4.05);
    AntiTrackCuts->SetPtRange(0.4, 4.05);

    TrackCuts->SetEtaRange(-0.77, 0.77);
    AntiTrackCuts->SetEtaRange(-0.77, 0.77);

    TrackCuts->SetNClsTPC(70);
    AntiTrackCuts->SetNClsTPC(70);

    config->SetDeltaEtaMax(0.019);
    config->SetDeltaPhiMax(0.019);

    v0Cuts->SetCutCPA(0.995);
    Antiv0Cuts->SetCutCPA(0.995);

    Posv0Daug->SetPID(AliPID::kProton, 999.9, 4);
    Negv0Daug->SetPID(AliPID::kPion, 999.9, 4);
    PosAntiv0Daug->SetPID(AliPID::kPion, 999.9, 4);
    NegAntiv0Daug->SetPID(AliPID::kProton, 999.9, 4);

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

  } else if (suffix == "24") {
    TrackCuts->SetEtaRange(-0.77, 0.77);
    AntiTrackCuts->SetEtaRange(-0.77, 0.77);

    TrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);
    AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);

    TrackCuts->SetNClsTPC(90);
    AntiTrackCuts->SetNClsTPC(90);

    config->SetDeltaEtaMax(0.019);
    config->SetDeltaPhiMax(0.019);

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

  } else if (suffix == "25") {
    TrackCuts->SetPtRange(0.4, 4.05);
    AntiTrackCuts->SetPtRange(0.4, 4.05);

    TrackCuts->SetEtaRange(-0.83, 0.83);
    AntiTrackCuts->SetEtaRange(-0.83, 0.83);

    config->SetDeltaEtaMax(0.019);
    config->SetDeltaPhiMax(0.019);

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

  } else if (suffix == "26") {
    TrackCuts->SetPtRange(0.6, 4.05);
    AntiTrackCuts->SetPtRange(0.6, 4.05);

    TrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);
    AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);

    config->SetDeltaEtaMax(0.019);
    config->SetDeltaPhiMax(0.019);

    Posv0Daug->SetEtaRange(-0.77, 0.77);
    Negv0Daug->SetEtaRange(-0.77, 0.77);
    PosAntiv0Daug->SetEtaRange(-0.77, 0.77);
    NegAntiv0Daug->SetEtaRange(-0.77, 0.77);

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

  } else if (suffix == "28") {
    TrackCuts->SetNClsTPC(90);
    AntiTrackCuts->SetNClsTPC(90);

    v0Cuts->SetCutCPA(0.995);
    Antiv0Cuts->SetCutCPA(0.995);

    v0Cuts->SetCutDCADaugToPrimVtx(0.06);
    Antiv0Cuts->SetCutDCADaugToPrimVtx(0.06);

  } else if (suffix == "29") {
    TrackCuts->SetPtRange(0.4, 4.05);
    AntiTrackCuts->SetPtRange(0.4, 4.05);

    TrackCuts->SetEtaRange(-0.77, 0.77);
    AntiTrackCuts->SetEtaRange(-0.77, 0.77);

    TrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);
    AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);

    config->SetDeltaEtaMax(0.019);
    config->SetDeltaPhiMax(0.019);

    Posv0Daug->SetNClsTPC(80);
    Negv0Daug->SetNClsTPC(80);
    PosAntiv0Daug->SetNClsTPC(80);
    NegAntiv0Daug->SetNClsTPC(80);

    v0Cuts->SetCutDCADaugToPrimVtx(0.06);
    Antiv0Cuts->SetCutDCADaugToPrimVtx(0.06);
  } else if (suffix == "30") {
    TrackCuts->SetEtaRange(-0.85, 0.85);
    AntiTrackCuts->SetEtaRange(-0.85, 0.85);

    config->SetDeltaEtaMax(0.019);
    config->SetDeltaPhiMax(0.019);

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
  } else if (suffix == "31") {
    TrackCuts->SetPtRange(0.4, 4.05);
    AntiTrackCuts->SetPtRange(0.4, 4.05);

    TrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);
    AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);

    TrackCuts->SetNClsTPC(70);
    AntiTrackCuts->SetNClsTPC(70);

    config->SetDeltaEtaMax(0.019);
    config->SetDeltaPhiMax(0.019);

    v0Cuts->SetCutCPA(0.995);
    Antiv0Cuts->SetCutCPA(0.995);

    Posv0Daug->SetEtaRange(-0.83, 0.83);
    Negv0Daug->SetEtaRange(-0.83, 0.83);
    PosAntiv0Daug->SetEtaRange(-0.83, 0.83);
    NegAntiv0Daug->SetEtaRange(-0.83, 0.83);

    v0Cuts->SetCutDCADaugTov0Vtx(1.2);
    Antiv0Cuts->SetCutDCADaugTov0Vtx(1.2);
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
  } else if (suffix == "33") {
    TrackCuts->SetEtaRange(-0.77, 0.77);
    AntiTrackCuts->SetEtaRange(-0.77, 0.77);

    config->SetDeltaEtaMax(0.019);
    config->SetDeltaPhiMax(0.019);

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
  } else if (suffix == "34") {
    TrackCuts->SetEtaRange(-0.77, 0.77);
    AntiTrackCuts->SetEtaRange(-0.77, 0.77);

    TrackCuts->SetNClsTPC(70);
    AntiTrackCuts->SetNClsTPC(70);

    config->SetDeltaEtaMax(0.019);
    config->SetDeltaPhiMax(0.019);

    Posv0Daug->SetNClsTPC(80);
    Negv0Daug->SetNClsTPC(80);
    PosAntiv0Daug->SetNClsTPC(80);
    NegAntiv0Daug->SetNClsTPC(80);

    Posv0Daug->SetEtaRange(-0.77, 0.77);
    Negv0Daug->SetEtaRange(-0.77, 0.77);
    PosAntiv0Daug->SetEtaRange(-0.77, 0.77);
    NegAntiv0Daug->SetEtaRange(-0.77, 0.77);

  } else if (suffix == "35") {
    TrackCuts->SetPtRange(0.4, 4.05);
    AntiTrackCuts->SetPtRange(0.4, 4.05);

    TrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);
    AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);

    TrackCuts->SetNClsTPC(70);
    AntiTrackCuts->SetNClsTPC(70);

    config->SetDeltaEtaMax(0.019);
    config->SetDeltaPhiMax(0.019);

    v0Cuts->SetCutCPA(0.995);
    Antiv0Cuts->SetCutCPA(0.995);

    Posv0Daug->SetNClsTPC(80);
    Negv0Daug->SetNClsTPC(80);
    PosAntiv0Daug->SetNClsTPC(80);
    NegAntiv0Daug->SetNClsTPC(80);

    v0Cuts->SetCutDCADaugTov0Vtx(1.2);
    Antiv0Cuts->SetCutDCADaugTov0Vtx(1.2);
  } else if (suffix == "36") {
    TrackCuts->SetPtRange(0.6, 4.05);
    AntiTrackCuts->SetPtRange(0.6, 4.05);

    TrackCuts->SetEtaRange(-0.85, 0.85);
    AntiTrackCuts->SetEtaRange(-0.85, 0.85);

    TrackCuts->SetNClsTPC(70);
    AntiTrackCuts->SetNClsTPC(70);

    config->SetDeltaEtaMax(0.019);
    config->SetDeltaPhiMax(0.019);

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
  } else if (suffix == "37") {
    TrackCuts->SetEtaRange(-0.85, 0.85);
    AntiTrackCuts->SetEtaRange(-0.85, 0.85);

    TrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);
    AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);

    TrackCuts->SetNClsTPC(90);
    AntiTrackCuts->SetNClsTPC(90);

    config->SetDeltaEtaMax(0.019);
    config->SetDeltaPhiMax(0.019);

    v0Cuts->SetCutCPA(0.995);
    Antiv0Cuts->SetCutCPA(0.995);

    Posv0Daug->SetEtaRange(-0.83, 0.83);
    Negv0Daug->SetEtaRange(-0.83, 0.83);
    PosAntiv0Daug->SetEtaRange(-0.83, 0.83);
    NegAntiv0Daug->SetEtaRange(-0.83, 0.83);

    v0Cuts->SetCutDCADaugTov0Vtx(1.2);
    Antiv0Cuts->SetCutDCADaugTov0Vtx(1.2);

  } else if (suffix == "38") {
    TrackCuts->SetEtaRange(-0.77, 0.77);
    AntiTrackCuts->SetEtaRange(-0.77, 0.77);

    TrackCuts->SetNClsTPC(90);
    AntiTrackCuts->SetNClsTPC(90);

    config->SetDeltaEtaMax(0.019);
    config->SetDeltaPhiMax(0.019);

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
  } else if (suffix == "39") {
    TrackCuts->SetPtRange(0.4, 4.05);
    AntiTrackCuts->SetPtRange(0.4, 4.05);

    TrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);
    AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);

    TrackCuts->SetNClsTPC(90);
    AntiTrackCuts->SetNClsTPC(90);

    config->SetDeltaEtaMax(0.019);
    config->SetDeltaPhiMax(0.019);

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
  } else if (suffix == "41") {
    TrackCuts->SetPtRange(0.4, 4.05);
    AntiTrackCuts->SetPtRange(0.4, 4.05);

    TrackCuts->SetEtaRange(-0.77, 0.77);
    AntiTrackCuts->SetEtaRange(-0.77, 0.77);

    TrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);
    AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);

    config->SetDeltaEtaMax(0.019);
    config->SetDeltaPhiMax(0.019);

    v0Cuts->SetCutCPA(0.995);
    Antiv0Cuts->SetCutCPA(0.995);

    Posv0Daug->SetEtaRange(-0.77, 0.77);
    Negv0Daug->SetEtaRange(-0.77, 0.77);
    PosAntiv0Daug->SetEtaRange(-0.77, 0.77);
    NegAntiv0Daug->SetEtaRange(-0.77, 0.77);

    v0Cuts->SetCutDCADaugTov0Vtx(1.2);
    Antiv0Cuts->SetCutDCADaugTov0Vtx(1.2);
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
  } else if (suffix == "43") {
    TrackCuts->SetPtRange(0.6, 4.05);
    AntiTrackCuts->SetPtRange(0.6, 4.05);

    TrackCuts->SetEtaRange(-0.85, 0.85);
    AntiTrackCuts->SetEtaRange(-0.85, 0.85);

    TrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);
    AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);

    config->SetDeltaEtaMax(0.019);
    config->SetDeltaPhiMax(0.019);

    Posv0Daug->SetEtaRange(-0.83, 0.83);
    Negv0Daug->SetEtaRange(-0.83, 0.83);
    PosAntiv0Daug->SetEtaRange(-0.83, 0.83);
    NegAntiv0Daug->SetEtaRange(-0.83, 0.83);
  } else if (suffix == "44") {
    TrackCuts->SetEtaRange(-0.85, 0.85);
    AntiTrackCuts->SetEtaRange(-0.85, 0.85);

    TrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);
    AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);

    config->SetDeltaEtaMax(0.019);
    config->SetDeltaPhiMax(0.019);

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
  else if (suffix == "100") {
    TrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);
    AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);
  }
  else if (suffix == "101") {
    TrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);
    AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);
  }
  else if (suffix == "102") {
    Posv0Daug->SetPID(AliPID::kProton, 999.9, 4);
    Negv0Daug->SetPID(AliPID::kPion, 999.9, 4);
    PosAntiv0Daug->SetPID(AliPID::kPion, 999.9, 4);
    NegAntiv0Daug->SetPID(AliPID::kProton, 999.9, 4);
  }



  if (isMC) {
    config->SetMomentumResolution(true);
  }

  if (fullBlastQA) {
    config->SetkTBinning(true);
    config->SetPtQA(true);
    config->SetMassQA(true);
  }

  if (!fullBlastQA) {
    config->SetMinimalBookingME(true);
    config->SetMinimalBookingSample(true);
  }

  TString addon = "PL";
  TString file = AliAnalysisManager::GetCommonFileName();

  TString EvtCutsName = Form("%sEvtCuts_%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *coutputEvtCuts = mgr->CreateContainer(
      EvtCutsName.Data(), TList::Class(), AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), EvtCutsName.Data()));

  TString TrackCutsName = Form("%sTrackCuts_%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *couputTrkCuts = mgr->CreateContainer(
      TrackCutsName.Data(), TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), TrackCutsName.Data()));

  TString AntiTrackCutsName = Form("%sAntiTrackCuts_%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *coutputAntiTrkCuts = mgr->CreateContainer(
      AntiTrackCutsName.Data(), TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), AntiTrackCutsName.Data()));

  AliAnalysisDataContainer *coutputv0Cuts;
  TString v0CutsName = Form("%sv0Cuts_%s", addon.Data(), suffix.Data());
  coutputv0Cuts = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      v0CutsName.Data(),
      TList::Class(), AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), v0CutsName.Data()));

  AliAnalysisDataContainer *coutputAntiv0Cuts;
  TString Antiv0CutsName = Form("%sAntiv0Cuts_%s", addon.Data(), suffix.Data());
  coutputAntiv0Cuts = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      Antiv0CutsName.Data(),
      TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), Antiv0CutsName.Data()));

  AliAnalysisDataContainer *coutputResults;
  TString ResultsName = Form("%sResults_%s", addon.Data(), suffix.Data());
  coutputResults = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      ResultsName.Data(),
      TList::Class(), AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), ResultsName.Data()));

  AliAnalysisDataContainer *coutputResultsQA;
  TString ResultsQAName = Form("%sResultsQA_%s", addon.Data(), suffix.Data());
  coutputResultsQA = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      ResultsQAName.Data(),
      TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), ResultsQAName.Data()));

  AliAnalysisDataContainer *coutputResultsSample;
  TString ResultsSampleName = Form("%sResultsSample_%s", addon.Data(), suffix.Data());
  coutputResultsSample = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      ResultsSampleName.Data(),
      TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), ResultsSampleName.Data()));

  AliAnalysisDataContainer *coutputResultsSampleQA;
  TString ResultsSampleQAName = Form("%sResultsSampleQA_%s", addon.Data(), suffix.Data());
  coutputResultsSampleQA = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      ResultsSampleQAName.Data(),
      TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), ResultsSampleQAName.Data()));

  AliAnalysisDataContainer *coutputTrkCutsMC;
  AliAnalysisDataContainer *coutputAntiTrkCutsMC;
  AliAnalysisDataContainer *coutputv0CutsMC;
  AliAnalysisDataContainer *coutputAntiv0CutsMC;
  if (isMC) {
    TString TrkCutsMCName = Form("%sTrkCutsMC_%s", addon.Data(), suffix.Data());
    coutputTrkCutsMC = mgr->CreateContainer(
        //@suppress("Invalid arguments") it works ffs
        TrkCutsMCName.Data(),
        TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), TrkCutsMCName.Data()));

    TString AntiTrkCutsMCName = Form("%sAntiTrkCutsMC_%s", addon.Data(), suffix.Data());
    coutputAntiTrkCutsMC = mgr->CreateContainer(
        //@suppress("Invalid arguments") it works ffs
        AntiTrkCutsMCName.Data(),
        TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), AntiTrkCutsMCName.Data()));

    TString v0CutsMCName = Form("%sv0CutsMC_%s", addon.Data(), suffix.Data());
    coutputv0CutsMC = mgr->CreateContainer(
        //@suppress("Invalid arguments") it works ffs
        v0CutsMCName.Data(),
        TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), v0CutsMCName.Data()));

    TString Antiv0CutsMCName = Form("%sAntiv0CutsMC_%s", addon.Data(), suffix.Data());
    coutputAntiv0CutsMC = mgr->CreateContainer(
        //@suppress("Invalid arguments") it works ffs
        Antiv0CutsMCName.Data(),
        TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), Antiv0CutsMCName.Data()));

  }

  TString ThreeBodyName = Form("%sThreeBody_%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *coutputThreeBody = mgr->CreateContainer(
    //@suppress("Invalid arguments") it works ffs
    ThreeBodyName.Data(),
    TList::Class(),
    AliAnalysisManager::kOutputContainer,
    Form("%s:%s", file.Data(), ThreeBodyName.Data()));

  AliAnalysisDataContainer *coutputEvtCutsTrigger;
  AliAnalysisDataContainer *couputTrkCutsTrigger;
  AliAnalysisDataContainer *coutputAntiTrkCutsTrigger;
  AliAnalysisDataContainer *coutputv0CutsTrigger;
  AliAnalysisDataContainer *coutputAntiv0CutsTrigger;
 

 
  bool turnoffClosePairRejectionCompletely = false;
  bool ClosePairRejectionPPPorPPL = true;

  AliAnalysisTaskThreeBodyFemto* taskNano;
  AliAnalysisTaskThreeBodyFemtoAOD* taskAOD;
  if(isNano){
    taskNano= new AliAnalysisTaskThreeBodyFemto("femtoNanoThreeBody", isMC);
    if (!fullBlastQA)
    { 
      taskNano->SetRunTaskLightWeight(true);
    }

    if (trigger == 0) { 
        taskNano->SelectCollisionCandidates(AliVEvent::kHighMultV0);  
      } else if (trigger == 1){     
        taskNano->SelectCollisionCandidates(AliVEvent::kINT7);  
      } 
      //Should I plot all the little things?
    if(plotAdditionalPlots == false){
      taskNano->SetRunPlotInvMassTriplet(false);
      taskNano->SetRunPlotQ3Vsq(false);
      taskNano->SetRunPlotPhiTheta(false);
      taskNano->SetRunPlotOtherHistos(false);
      taskNano->SetRunPlotMult(false);
    }
    if(plotAdditionalPlots == true){
      taskNano->SetRunPlotInvMassTriplet(true);
      taskNano->SetRunPlotQ3Vsq(true);
      taskNano->SetRunPlotPhiTheta(true);
      taskNano->SetRunPlotOtherHistos(true);
      taskNano->SetRunPlotMult(true);
    }
    taskNano->SetEventCuts(evtCuts);  
    taskNano->SetProtonCuts(TrackCuts); 
    taskNano->SetAntiProtonCuts(AntiTrackCuts); 
    taskNano->Setv0Cuts(v0Cuts);  
    taskNano->SetAntiv0Cuts(Antiv0Cuts);  
    taskNano->SetCorrelationConfig(config); 
    taskNano->SetRunThreeBodyHistograms(true);
    taskNano->SetClosePairRejectionForAll(ClosePairRejectionForAll);
    taskNano->SetturnoffClosePairRejectionCompletely(turnoffClosePairRejectionCompletely);
    taskNano->SetClosePairRejectionPPPorPPL(ClosePairRejectionPPPorPPL);


    taskNano->SetMixingChoice(mixinfChoice);

    mgr->AddTask(taskNano); 
    
    mgr->ConnectInput(taskNano, 0, cinput); 
    mgr->ConnectOutput(taskNano, 1, coutputEvtCuts);  
    mgr->ConnectOutput(taskNano, 2, couputTrkCuts); 
    mgr->ConnectOutput(taskNano, 3, coutputAntiTrkCuts);  
    mgr->ConnectOutput(taskNano, 4, coutputv0Cuts); 
    mgr->ConnectOutput(taskNano, 5, coutputAntiv0Cuts); 
    mgr->ConnectOutput(taskNano, 6, coutputResults);  
    mgr->ConnectOutput(taskNano, 7, coutputResultsQA);  
    mgr->ConnectOutput(taskNano, 8, coutputResultsSample);  
    mgr->ConnectOutput(taskNano, 9, coutputResultsSampleQA);  
    mgr->ConnectOutput(taskNano, 10, coutputThreeBody);  
    if (isMC) { 
      mgr->ConnectOutput(taskNano, 11, coutputTrkCutsMC); 
      mgr->ConnectOutput(taskNano, 12, coutputAntiTrkCutsMC); 
      mgr->ConnectOutput(taskNano, 13, coutputv0CutsMC);  
      mgr->ConnectOutput(taskNano, 14, coutputAntiv0CutsMC);  
    } 
  }
  else{
    taskAOD= new AliAnalysisTaskThreeBodyFemtoAOD("femtoAODThreeBody", isMC, false);
    if (!fullBlastQA)
    { 
      taskAOD->SetRunTaskLightWeight(true);
    }

    if (trigger == 0) { 
        taskAOD->SelectCollisionCandidates(AliVEvent::kHighMultV0);  
      } else if (trigger == 1){     
        taskAOD->SelectCollisionCandidates(AliVEvent::kINT7);  
      } 
    taskAOD->SetEventCuts(evtCuts);  
    taskAOD->SetProtonCuts(TrackCuts); 
    taskAOD->SetAntiProtonCuts(AntiTrackCuts); 
    taskAOD->Setv0Cuts(v0Cuts);  
    taskAOD->SetAntiv0Cuts(Antiv0Cuts);   
    taskAOD->SetCorrelationConfig(config); 
    taskAOD->SetRunThreeBodyHistograms(true);
    taskAOD->SetTriggerOn(false);
    taskAOD->SetIsMC(isMC);
    mgr->AddTask(taskAOD); 
    
    mgr->ConnectInput(taskAOD, 0, cinput); 
    mgr->ConnectOutput(taskAOD, 1, coutputEvtCuts);  
    mgr->ConnectOutput(taskAOD, 2, couputTrkCuts); 
    mgr->ConnectOutput(taskAOD, 3, coutputAntiTrkCuts);  
    mgr->ConnectOutput(taskAOD, 4, coutputv0Cuts); 
    mgr->ConnectOutput(taskAOD, 5, coutputAntiv0Cuts); 
    mgr->ConnectOutput(taskAOD, 6, coutputResults);  
    mgr->ConnectOutput(taskAOD, 7, coutputResultsQA);  
    mgr->ConnectOutput(taskAOD, 8, coutputResultsSample);  
    mgr->ConnectOutput(taskAOD, 9, coutputResultsSampleQA); 
    mgr->ConnectOutput(taskAOD, 10, coutputThreeBody); 
    if (isMC) { 
      mgr->ConnectOutput(taskAOD, 16, coutputTrkCutsMC); 
      mgr->ConnectOutput(taskAOD, 17, coutputAntiTrkCutsMC); 
      mgr->ConnectOutput(taskAOD, 18, coutputv0CutsMC);  
      mgr->ConnectOutput(taskAOD, 19, coutputAntiv0CutsMC);  
    } 
  }

      
  if (isNano) {
    return taskNano;
  } else {
    return taskAOD;
  }

  
}
