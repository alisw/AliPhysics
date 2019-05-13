#include <vector>
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskNanoXioton.h"
#include "AliFemtoDreamEventCuts.h"
#include "AliFemtoDreamTrackCuts.h"
#include "AliFemtoDreamCascadeCuts.h"
#include "AliFemtoDreamCollConfig.h"

AliAnalysisTaskSE *AddTaskFemtoXoton(bool fullBlastQA = false,
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
      false, true, false, false);
  TrackCuts->SetFilterBit(128);
  TrackCuts->SetCutCharge(1);

  AliFemtoDreamTrackCuts *AntiTrackCuts =
      AliFemtoDreamTrackCuts::PrimProtonCuts(false, true, false, false);
  AntiTrackCuts->SetFilterBit(128);
  AntiTrackCuts->SetCutCharge(-1);

  if (suffix != "0" && suffix != "999") {
    evtCuts->SetMinimalBooking(true);
    TrackCuts->SetMinimalBooking(true);
    AntiTrackCuts->SetMinimalBooking(true);
  }

  AliFemtoDreamCollConfig *config = new AliFemtoDreamCollConfig("Femto",
                                                                "Femto");

  // Femto Collection
  std::vector<int> PDGParticles;
  PDGParticles.push_back(2212);
  PDGParticles.push_back(2212);
  //  PDGParticles.push_back(3212);
  //  PDGParticles.push_back(3212);

  std::vector<int> NBins;
  std::vector<float> kMin;
  std::vector<float> kMax;
  std::vector<int> pairQA;
  std::vector<bool> closeRejection;
  //pairs:
  //pp                1
  //p bar p           2
  //bar p bar p       3
  const int nPairs = 3;
  for (int i = 0; i < nPairs; ++i) {
    pairQA.push_back(0);
    closeRejection.push_back(false);
    if (suffix == "0") {
      NBins.push_back(750);
      kMin.push_back(0.);
      kMax.push_back(3.);
    } else {
      NBins.push_back(250);
      kMin.push_back(0.);
      kMax.push_back(1.);
    }
  }
  pairQA[0] = 11;
  pairQA[2] = 11;

  closeRejection[0] = true;  // pp
  closeRejection[2] = true;  // barp barp

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

  if (suffix == "0" && fullBlastQA) {
    config->SetPhiEtaBinnign(true);
    config->SetkTBinning(true);
    config->SetmTBinning(true);
    config->SetPtQA(true);
  }
  config->SetdPhidEtaPlots(false);
  config->SetInvMassPairs(false);

  if (suffix != "0") {
    config->SetMinimalBookingME(true);
  }
  AliAnalysisTaskNanoXioton* task = new AliAnalysisTaskNanoXioton("femtoXoton");
  if (suffix != "0" && suffix != "999") {
    task->SetRunTaskLightWeight(true);
  }
  task->SetEventCuts(evtCuts);
  task->SetProtonCuts(TrackCuts);
  task->SetAntiProtonCuts(AntiTrackCuts);

  mgr->AddTask(task);

  TString addon = "HM";

  TString file = AliAnalysisManager::GetCommonFileName();

  mgr->ConnectInput(task, 0, cinput);

  TString EvtCutsName = Form("%sEvtCuts%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *coutputEvtCuts = mgr->CreateContainer(
      EvtCutsName.Data(), TList::Class(), AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), EvtCutsName.Data()));
  mgr->ConnectOutput(task, 1, coutputEvtCuts);

  TString TrackCutsName = Form("%sTrackCuts%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *couputTrkCuts =
      mgr->CreateContainer(TrackCutsName.Data(), TList::Class(),
                           AliAnalysisManager::kOutputContainer,
                           Form("%s:%s", file.Data(), TrackCutsName.Data()));
  mgr->ConnectOutput(task, 2, couputTrkCuts);

  TString AntiTrackCutsName =
      Form("%sAntiTrackCuts%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *coutputAntiTrkCuts = mgr->CreateContainer(
      AntiTrackCutsName.Data(), TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), AntiTrackCutsName.Data()));
  mgr->ConnectOutput(task, 3, coutputAntiTrkCuts);

  return task;
}
