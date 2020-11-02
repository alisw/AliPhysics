#ifndef __CINT__
#include <vector>
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliFemtoDreamEventCuts.h"
#include "AliFemtoDreamTrackCuts.h"
#include "AliFemtoDreamCollConfig.h"
#include "AliAnalysisTaskCharmingFemto.h"
#endif

AliAnalysisTaskSE *AddTaskCharmingFemto(bool isMC = false, bool fullBlastQA =
                                            false,
                                        TString trigger = "kINT7",
                                        const char *cutVariation = "0") {
  TString suffix = TString::Format("%s", cutVariation);

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskCharmingFemto()", "No analysis manager found.");
    return nullptr;
  }

  // ================== GetInputEventHandler =============================
  AliVEventHandler *inputHandler = mgr->GetInputEventHandler();
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();

  //========= Init subtasks and start analyis ============================
  // Event Cuts
  AliFemtoDreamEventCuts *evtCuts = AliFemtoDreamEventCuts::StandardCutsRun2();
  evtCuts->CleanUpMult(false, false, false, true);

  // =====================================================================
  // Proton cut variations
  const float ProtonPtlow = 0.4;
  const float ProtonPtup = 0.6;
  const float ProtonEtaLow = 0.75;
  const float ProtonEtaUp = 0.85;
  const float ProtonNsigmaLow = 2.5;
  const float ProtonNsigmaUp = 3.5;
  const float ProtonNClsLow = 70;
  const float ProtonNClsUp = 90;

  // Track Cuts
  AliFemtoDreamTrackCuts *TrackCuts = AliFemtoDreamTrackCuts::PrimProtonCuts(
      isMC, true, false, false);
  TrackCuts->SetFilterBit(128);
  TrackCuts->SetCutCharge(1);

  AliFemtoDreamTrackCuts *AntiTrackCuts =
      AliFemtoDreamTrackCuts::PrimProtonCuts(isMC, true, false, false);
  AntiTrackCuts->SetFilterBit(128);
  AntiTrackCuts->SetCutCharge(-1);

  TrackCuts->SetMinimalBooking(suffix != "0");
  AntiTrackCuts->SetMinimalBooking(suffix != "0");

  // =====================================================================
  // D mesons - to be done


  // =====================================================================
  // Femto Collection
  std::vector<int> PDGParticles;
  PDGParticles.push_back(2212);  // 0 - proton
  PDGParticles.push_back(2212);  // 1 - antiproton
  PDGParticles.push_back(411);   // 2 - dplus
  PDGParticles.push_back(411);   // 3 - dminus

  std::vector<float> ZVtxBins = AliFemtoDreamCollConfig::GetDefaultZbins();

  std::vector<int> NBins;
  std::vector<float> kMin;
  std::vector<float> kMax;
  std::vector<int> pairQA;
  std::vector<bool> closeRejection;
  const int nPairs = 10;
  for (int i = 0; i < nPairs; ++i) {
    pairQA.push_back(0);
    closeRejection.push_back(false);
    NBins.push_back((suffix == "0") ? 600 : 200);
    kMin.push_back(0.);
    kMax.push_back((suffix == "0") ? 3. : 1.);
  }

  NBins[0] = (suffix == "0") ? 3000 : 1000;  // pp
  NBins[4] = (suffix == "0") ? 3000 : 1000;  // barp barp

  closeRejection[0] = true;  // pp
  closeRejection[4] = true;  // barp barp

  pairQA[0] = 11;   // pp
  pairQA[4] = 11;   // pbarpbar
  pairQA[2] = 10;   // pDplus
  pairQA[6] = 10;   // barp dminus

  AliFemtoDreamCollConfig *config = new AliFemtoDreamCollConfig("Femto",
                                                                "Femto");
  std::vector<int> MultBins;
  if (trigger == "kHighMultV0") {
    std::vector<int> MultBins = AliFemtoDreamCollConfig::GetHMMultBins();
    config->SetMultBins(MultBins);
  } else {
    std::vector<int> MultBins = AliFemtoDreamCollConfig::GetMBMultBins();
    config->SetMultBins(MultBins);
  }
  config->SetMultBinning(true);

  config->SetExtendedQAPairs(pairQA);
  config->SetdPhidEtaPlotsSmallK(false);
  config->SetZBins(ZVtxBins);
  config->SetMomentumResolution(isMC);

  config->SetDeltaEtaMax(0.012);
  config->SetDeltaPhiMax(0.012);
  config->SetClosePairRejection(closeRejection);

  config->SetPhiEtaBinnign((suffix == "0" && fullBlastQA));
  config->SetkTBinning((suffix == "0" && fullBlastQA));
  config->SetmTBinning((suffix == "0" && fullBlastQA));

  config->SetPtQA(true);
  config->SetMassQA(true);
  config->SetdPhidEtaPlots(false);
  config->SetPDGCodes(PDGParticles);
  config->SetNBinsHist(NBins);
  config->SetMinKRel(kMin);
  config->SetMaxKRel(kMax);
  config->SetMixingDepth(10);
  config->SetUseEventMixing(true);
  config->SetMultiplicityEstimator(AliFemtoDreamEvent::kRef08);
  config->SetMinimalBookingME(false);

  AliAnalysisTaskCharmingFemto *task = new AliAnalysisTaskCharmingFemto(
      "AliAnalysisTaskCharmingFemto", isMC);

  task->SetEventCuts(evtCuts);
  task->SetProtonCuts(TrackCuts);
  task->SetAntiProtonCuts(AntiTrackCuts);
  task->SetCollectionConfig(config);

  task->SetLightweight(suffix != "0");

  mgr->AddTask(task);

  TString addon = "";
  if (trigger == "kINT7") {
    addon += "MB_CharmFemto_";
  } else if (trigger == "kHighMultV0") {
    addon += "HM_CharmFemto_";
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

  AliAnalysisDataContainer *coutputResults;
  TString ResultsName = Form("%sResults%s", addon.Data(), suffix.Data());
  coutputResults = mgr->CreateContainer(
      ResultsName.Data(), TList::Class(), AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), ResultsName.Data()));
  mgr->ConnectOutput(task, 5, coutputResults);

  AliAnalysisDataContainer *coutputResultQA;
  TString ResultQAName = Form("%sResultQA%s", addon.Data(), suffix.Data());
  coutputResultQA = mgr->CreateContainer(
      ResultQAName.Data(), TList::Class(), AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), ResultQAName.Data()));
  mgr->ConnectOutput(task, 6, coutputResultQA);

  if (isMC) {
    TString TrkCutsMCName = Form("%sTrackCutsMC%s", addon.Data(),
                                 suffix.Data());
    AliAnalysisDataContainer *coutputTrkCutsMC = mgr->CreateContainer(
        TrkCutsMCName.Data(), TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), TrkCutsMCName.Data()));
    mgr->ConnectOutput(task, 7, coutputTrkCutsMC);

    TString AntiTrkCutsMCName = Form("%sAntiTrackCutsMC%s", addon.Data(),
                                     suffix.Data());
    AliAnalysisDataContainer *coutputAntiTrkCutsMC = mgr->CreateContainer(
        AntiTrkCutsMCName.Data(), TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), AntiTrkCutsMCName.Data()));
    mgr->ConnectOutput(task, 8, coutputAntiTrkCutsMC);
  }

  return task;
}
