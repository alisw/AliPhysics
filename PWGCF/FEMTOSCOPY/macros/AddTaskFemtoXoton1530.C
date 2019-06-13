#include <AliAnalysisTaskNanoXioton1530.h>
#include <vector>
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliFemtoDreamEventCuts.h"
#include "AliFemtoDreamTrackCuts.h"
#include "AliFemtoDreamCascadeCuts.h"
#include "AliFemtoDreamCollConfig.h"

AliAnalysisTaskSE *AddTaskFemtoXoton1530(bool fullBlastQA = false,
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
      false, true, false, false);
  TrackCuts->SetFilterBit(128);
  TrackCuts->SetCutCharge(1);

  AliFemtoDreamTrackCuts *AntiTrackCuts =
      AliFemtoDreamTrackCuts::PrimProtonCuts(false, true, false, false);
  AntiTrackCuts->SetFilterBit(128);
  AntiTrackCuts->SetCutCharge(-1);

  //Track Cuts
  AliFemtoDreamTrackCuts *PiPCuts = new AliFemtoDreamTrackCuts();
  PiPCuts->SetPtRange(0.15, 4.0);
  PiPCuts->SetEtaRange(-0.8, 0.8);
  PiPCuts->SetNClsTPC(70);
  PiPCuts->SetDCAReCalculation(true);
  PiPCuts->SetDCAVtxXY(0.3);
  PiPCuts->SetDCAVtxZ(0.3);
  PiPCuts->SetFilterBit(96);
  PiPCuts->SetCutSharedCls(true);
  PiPCuts->SetCutTPCCrossedRows(true, 70, 0.5);
  PiPCuts->SetPID(AliPID::kPion, 99.5, 3);
//  PiPCuts->SetCutSmallestSig(true);
  PiPCuts->SetCutCharge(1);

  AliFemtoDreamTrackCuts *PiMCuts = new AliFemtoDreamTrackCuts();
  PiMCuts->SetPtRange(0.15, 4.0);
  PiMCuts->SetEtaRange(-0.8, 0.8);
  PiMCuts->SetNClsTPC(70);
  PiMCuts->SetDCAReCalculation(true);
  PiMCuts->SetDCAVtxXY(0.3);
  PiMCuts->SetDCAVtxZ(0.3);
  PiMCuts->SetFilterBit(96);
  PiMCuts->SetCutSharedCls(true);
  PiMCuts->SetCutTPCCrossedRows(true, 70, 0.5);
  PiMCuts->SetPID(AliPID::kPion, 99.5, 3);
//  PiMCuts->SetCutSmallestSig(true);
  PiMCuts->SetCutCharge(-1);

  //Cascade Cuts
  AliFemtoDreamCascadeCuts* CascadeCuts = AliFemtoDreamCascadeCuts::XiFor1530Cuts(
      false, false);
  CascadeCuts->SetXiCharge(-1);
  AliFemtoDreamTrackCuts *XiNegCuts = AliFemtoDreamTrackCuts::Xiv0PionCuts(
      false, true, false);
  XiNegCuts->SetCheckTPCRefit(false);  //for nanos this is already done while prefiltering
  AliFemtoDreamTrackCuts *XiPosCuts = AliFemtoDreamTrackCuts::Xiv0ProtonCuts(
      false, true, false);
  XiPosCuts->SetCheckTPCRefit(false);  //for nanos this is already done while prefiltering
  AliFemtoDreamTrackCuts *XiBachCuts = AliFemtoDreamTrackCuts::XiBachPionCuts(
      false, true, false);
  XiBachCuts->SetCheckTPCRefit(false);  //for nanos this is already done while prefiltering

  CascadeCuts->Setv0Negcuts(XiNegCuts);
  CascadeCuts->Setv0PosCuts(XiPosCuts);
  CascadeCuts->SetBachCuts(XiBachCuts);
  CascadeCuts->SetPDGCodeCasc(3312);
  CascadeCuts->SetPDGCodev0(3122);
  CascadeCuts->SetPDGCodePosDaug(2212);
  CascadeCuts->SetPDGCodeNegDaug(-211);
  CascadeCuts->SetPDGCodeBach(-211);

  AliFemtoDreamCascadeCuts* AntiCascadeCuts = AliFemtoDreamCascadeCuts::XiFor1530Cuts(
      false, false);
  AntiCascadeCuts->SetXiCharge(1);
  AliFemtoDreamTrackCuts *AntiXiNegCuts =
      AliFemtoDreamTrackCuts::Xiv0ProtonCuts(false, true, false);
  AntiXiNegCuts->SetCutCharge(-1);
  AntiXiNegCuts->SetCheckTPCRefit(false);  //for nanos this is already done while prefiltering
  AliFemtoDreamTrackCuts *AntiXiPosCuts = AliFemtoDreamTrackCuts::Xiv0PionCuts(
      false, true, false);
  AntiXiPosCuts->SetCutCharge(1);
  AntiXiPosCuts->SetCheckTPCRefit(false);  //for nanos this is already done while prefiltering
  AliFemtoDreamTrackCuts *AntiXiBachCuts =
      AliFemtoDreamTrackCuts::XiBachPionCuts(false, true, false);
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

  AliFemtoDreamv0Cuts *Xi1530Cuts = new AliFemtoDreamv0Cuts();
  Xi1530Cuts->SetIsMonteCarlo(false);
  Xi1530Cuts->SetAxisInvMassPlots(400, 1.46, 1.6);
  Xi1530Cuts->SetCutInvMass(0.007);
  Xi1530Cuts->SetPDGCodePosDaug(211);
  Xi1530Cuts->SetPDGCodeNegDaug(3312);
  Xi1530Cuts->SetPDGCodev0(3324);
  Xi1530Cuts->SetCutDaughters(false);

  AliFemtoDreamv0Cuts *AntiXi1530Cuts = new AliFemtoDreamv0Cuts();
  AntiXi1530Cuts->SetIsMonteCarlo(false);
  AntiXi1530Cuts->SetAxisInvMassPlots(400, 1.46, 1.6);
  AntiXi1530Cuts->SetCutInvMass(0.007);
  AntiXi1530Cuts->SetPDGCodePosDaug(3312);
  AntiXi1530Cuts->SetPDGCodeNegDaug(211);
  AntiXi1530Cuts->SetPDGCodev0(3324);
  AntiXi1530Cuts->SetCutDaughters(false);

  if (suffix != "0" && suffix != "999") {
    evtCuts->SetMinimalBooking(true);
    TrackCuts->SetMinimalBooking(true);
    AntiTrackCuts->SetMinimalBooking(true);
    PiPCuts->SetMinimalBooking(true);
    PiMCuts->SetMinimalBooking(true);
    CascadeCuts->SetMinimalBooking(true);
    AntiCascadeCuts->SetMinimalBooking(true);
    Xi1530Cuts->SetMinimalBooking(true);
    AntiXi1530Cuts->SetMinimalBooking(true);
  }

  AliFemtoDreamCollConfig *config = new AliFemtoDreamCollConfig("Femto",
                                                                "Femto");

  // Femto Collection
  std::vector<int> PDGParticles;
  PDGParticles.push_back(2212);
  PDGParticles.push_back(2212);
  PDGParticles.push_back(3324);
  PDGParticles.push_back(3324);

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
    if (suffix == "0") {
      NBins.push_back(1500);
      kMin.push_back(0.);
      kMax.push_back(6.);
    } else {
      NBins.push_back(250);
      kMin.push_back(0.);
      kMax.push_back(1.);
    }
  }
  pairQA[0] = 11;
  pairQA[4] = 11;
  pairQA[2] = 14;
  pairQA[6] = 14;

  closeRejection[0] = true;  // pp
  closeRejection[4] = true;  // barp barp

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

  config->SetMultBinning(true);
  config->SetdPhidEtaPlotsSmallK(false);
  config->SetdPhidEtaPlots(false);

  config->SetPhiEtaBinnign(false);

  if (suffix == "0" && fullBlastQA) {
    config->SetkTBinning(true);
    config->SetmTBinning(true);
    config->SetPtQA(true);
  }

  if (suffix != "0") {
    config->SetMinimalBookingME(true);
  }
  AliAnalysisTaskNanoXioton1530* task = new AliAnalysisTaskNanoXioton1530(
      "femtoXoton");
  if (suffix != "0" && suffix != "999") {
    task->SetRunTaskLightWeight(true);
  }
  task->SelectCollisionCandidates(AliVEvent::kHighMultV0);
  task->SetEventCuts(evtCuts);
  task->SetProtonCuts(TrackCuts);
  task->SetAntiProtonCuts(AntiTrackCuts);
  task->SetPionCuts(PiPCuts);
  task->SetAntiPionCuts(PiMCuts);
  task->SetXiCuts(CascadeCuts);
  task->SetAntiXiCuts(AntiCascadeCuts);
  task->SetXi1530Cuts(Xi1530Cuts);
  task->SetAntiXi1530Cuts(AntiXi1530Cuts);
  task->SetCorrelationConfig(config);
  mgr->AddTask(task);

  TString addon = "PXi";

  TString file = AliAnalysisManager::GetCommonFileName();

  mgr->ConnectInput(task, 0, cinput);

  TString EvtCutsName = Form("%sEvtCuts%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *coutputEvtCuts = mgr->CreateContainer(
      EvtCutsName.Data(), TList::Class(), AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), EvtCutsName.Data()));
  mgr->ConnectOutput(task, 1, coutputEvtCuts);

  TString TrackCutsName = Form("%sTrackCuts%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *couputTrkCuts = mgr->CreateContainer(
      TrackCutsName.Data(), TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), TrackCutsName.Data()));
  mgr->ConnectOutput(task, 2, couputTrkCuts);

  TString AntiTrackCutsName = Form("%sAntiTrackCuts%s", addon.Data(),
                                   suffix.Data());
  AliAnalysisDataContainer *coutputAntiTrkCuts = mgr->CreateContainer(
      AntiTrackCutsName.Data(), TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), AntiTrackCutsName.Data()));
  mgr->ConnectOutput(task, 3, coutputAntiTrkCuts);

  TString PionCutsName = Form("%sPionCuts%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *couputPionCuts = mgr->CreateContainer(
      PionCutsName.Data(), TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), PionCutsName.Data()));
  mgr->ConnectOutput(task, 4, couputPionCuts);

  TString AntiPionCutsName = Form("%sAntiPionCuts%s", addon.Data(),
                                   suffix.Data());
  AliAnalysisDataContainer *coutputAntiPionCuts = mgr->CreateContainer(
      AntiPionCutsName.Data(), TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), AntiPionCutsName.Data()));
  mgr->ConnectOutput(task, 5, coutputAntiPionCuts);

  AliAnalysisDataContainer *coutputCascadeCuts;
  TString CascadeCutsName = Form("%sCascadeCuts%s", addon.Data(),
                                 suffix.Data());
  coutputCascadeCuts = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      CascadeCutsName.Data(),
      TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), CascadeCutsName.Data()));
  mgr->ConnectOutput(task, 6, coutputCascadeCuts);

  AliAnalysisDataContainer *coutputAntiCascadeCuts;
  TString AntiCascadeCutsName = Form("%sAntiCascadeCuts%s", addon.Data(),
                                     suffix.Data());
  coutputAntiCascadeCuts = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      AntiCascadeCutsName.Data(),
      TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), AntiCascadeCutsName.Data()));
  mgr->ConnectOutput(task, 7, coutputAntiCascadeCuts);

  AliAnalysisDataContainer *coutputXi1530Cuts;
  TString Xi1530CutsName = Form("%sXi1530Cuts%s", addon.Data(),
                                 suffix.Data());
  coutputXi1530Cuts = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      Xi1530CutsName.Data(),
      TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), Xi1530CutsName.Data()));
  mgr->ConnectOutput(task, 8, coutputXi1530Cuts);

  AliAnalysisDataContainer *coutputAntiXi1530Cuts;
  TString AntiXi1530CutsName = Form("%sAntiXi1530Cuts%s", addon.Data(),
                                 suffix.Data());
  coutputAntiXi1530Cuts = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      AntiXi1530CutsName.Data(),
      TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), AntiXi1530CutsName.Data()));
  mgr->ConnectOutput(task, 9, coutputAntiXi1530Cuts);

  AliAnalysisDataContainer *coutputResults;
  TString ResultsName = Form("%sResults%s", addon.Data(), suffix.Data());
  coutputResults = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      ResultsName.Data(),
      TList::Class(), AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), ResultsName.Data()));
  mgr->ConnectOutput(task, 10, coutputResults);

  AliAnalysisDataContainer *coutputResultsQA;
  TString ResultsQAName = Form("%sResultsQA%s", addon.Data(), suffix.Data());
  coutputResultsQA = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      ResultsQAName.Data(),
      TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), ResultsQAName.Data()));
  mgr->ConnectOutput(task, 11, coutputResultsQA);

  return task;
}
