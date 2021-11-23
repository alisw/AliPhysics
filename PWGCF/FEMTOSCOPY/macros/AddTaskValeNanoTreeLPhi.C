#ifndef __CINT__
#include <vector>
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskValeNanoTreeLPhi.h"
#include "AliFemtoDreamEventCuts.h"
#include "AliFemtoDreamTrackCuts.h"
#include "AliFemtoDreamv0Cuts.h"
#include "AliFemtoDreamCollConfig.h"
#endif

AliAnalysisTaskSE *AddTaskValeNanoTreeLPhi(bool isMC = false,
                                           bool fullblastQA = false,
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
  v0Cuts->SetCutDCADaugToPrimVtx(0.03);
  Negv0Daug->SetEtaRange(-0.9, 0.9);
  Negv0Daug->SetNClsTPC(60);
  Negv0Daug->SetPID(AliPID::kPion, 999., 6);

  //Loosening cuts for Lambdas
  //1. DCA to V0 vtx (def. < 1.5 cm)
  //2. DCA daugh to PV (def > 0.05 cm)
  v0Cuts->SetPtRange(0.2, 999.);
  v0Cuts->SetCutMaxDecayVtx(110);
  v0Cuts->SetCutTransverseRadius(0.15, 110);
  v0Cuts->SetCutDCADaugTov0Vtx(1.7);
  v0Cuts->SetCutCPA(0.97);

  v0Cuts->SetPosDaugterTrackCuts(Posv0Daug);
  v0Cuts->SetNegDaugterTrackCuts(Negv0Daug);
  v0Cuts->SetPDGCodePosDaug(2212); //Proton
  v0Cuts->SetPDGCodeNegDaug(211);  //Pion
  v0Cuts->SetPDGCodev0(3122);      //Lambda

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
  Antiv0Cuts->SetCutDCADaugToPrimVtx(0.03);
  //Loosening cuts for AntiLambdas
  //1. DCA to V0 vtx (def. < 1.5 cm)
  //2. DCA daugh to PV (def > 0.05 cm)
  Antiv0Cuts->SetPtRange(0.2, 999.);
  Antiv0Cuts->SetCutMaxDecayVtx(110);
  Antiv0Cuts->SetCutTransverseRadius(0.15, 110);
  Antiv0Cuts->SetCutDCADaugTov0Vtx(1.7);
  Antiv0Cuts->SetCutCPA(0.97);

  Antiv0Cuts->SetPosDaugterTrackCuts(PosAntiv0Daug);
  Antiv0Cuts->SetNegDaugterTrackCuts(NegAntiv0Daug);
  Antiv0Cuts->SetPDGCodePosDaug(211);  //Pion
  Antiv0Cuts->SetPDGCodeNegDaug(2212); //Proton
  Antiv0Cuts->SetPDGCodev0(-3122);     //Lambda

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
    TrackCutsPhi->SetMinimalBooking(true);
  }

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

  std::vector<int> NBins;
  std::vector<float> kMin;
  std::vector<float> kMax;
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

  if (fullblastQA) {
    config->SetPtQA(true);
    config->SetMassQA(true);
  }
  if(!fullblastQA){
    evtCuts->SetMinimalBooking(true);
    v0Cuts->SetMinimalBooking(true);
    Antiv0Cuts->SetMinimalBooking(true);
    TrackPosKaonCuts->SetMinimalBooking(true);
    TrackNegKaonCuts->SetMinimalBooking(true);
    TrackCutsPhi->SetMinimalBooking(true);
  }

  AliAnalysisTaskValeNanoTreeLPhi *task =
      new AliAnalysisTaskValeNanoTreeLPhi("FemtoDreamLPhiTree", isMC);
  task->SetTrigger(AliVEvent::kHighMultV0);
  task->SelectCollisionCandidates(AliVEvent::kHighMultV0);

  if (!fullblastQA)
  {
    task->SetRunTaskLightWeight(true);
  }
  task->SetEventCuts(evtCuts);
  task->SetLambdaCuts(v0Cuts);
  task->SetAntiLambdaCuts(Antiv0Cuts);
  task->SetPhiCuts(TrackCutsPhi);
  task->SetPosKaonCuts(TrackPosKaonCuts);
  task->SetNegKaonCuts(TrackNegKaonCuts);
  task->SetCollectionConfig(config);

  mgr->AddTask(task);
  TString addon = "LPhi";

  TString file = AliAnalysisManager::GetCommonFileName();
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
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
  TString ResultsTreeName = Form("%sTree%s", addon.Data(), suffix.Data());
  coutputResultsTree = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      ResultsTreeName.Data(),
      TTree::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), ResultsTreeName.Data()));
  mgr->ConnectOutput(task, 10, coutputResultsTree);

  return task;
}