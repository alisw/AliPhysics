//PT for beginning of ProtonDeuteron
#include <vector>
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliFemtoDreamEventCuts.h"
#include "AliFemtoDreamTrackCuts.h"
#include "AliFemtoDreamCollConfig.h"
#include "AliAnalysisTaskNanoPt.h"

AliAnalysisTaskSE *AddTaskNanoPt(  bool isMC = true,
                                   TString trigger = "kINT7",
                                   bool DCAPlots = false,//3
                                   bool CombSigma = false,//4
                                   bool ContributionSplitting = false,
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
  //========= Init subtasks and start analyis ============================
  // Event Cuts
  AliFemtoDreamEventCuts *evtCuts = AliFemtoDreamEventCuts::StandardCutsRun2();
  evtCuts->CleanUpMult(false, false, false, true);

  // =====================================================================
  //Proton track Cuts----------------------------------------------------------------------------
  AliFemtoDreamTrackCuts *TrackCuts =
    AliFemtoDreamTrackCuts::PrimProtonCuts(isMC, true, CombSigma, ContributionSplitting);
  TrackCuts->SetMinimalBooking(false);
  TrackCuts->SetCutCharge(1);
  //Antiproton track Cuts-------------------------------------------------------------------------
  AliFemtoDreamTrackCuts *AntiTrackCuts =
    AliFemtoDreamTrackCuts::PrimProtonCuts(isMC, true, CombSigma, ContributionSplitting);
  AntiTrackCuts->SetMinimalBooking(false);
  AntiTrackCuts->SetCutCharge(-1);
  //deuteron track cuts----------------------------------------------------------------------------
  AliFemtoDreamTrackCuts *TrackCutsDeuteron = AliFemtoDreamTrackCuts::PrimDeuteronCuts(isMC, true,
      CombSigma, ContributionSplitting);
  TrackCutsDeuteron->SetMinimalBooking(false);
  TrackCutsDeuteron->SetCutCharge(1);
  //Antideuteron track cuts----------------------------------------------------------------------------
  AliFemtoDreamTrackCuts *AntiTrackCutsDeuteron = AliFemtoDreamTrackCuts::PrimDeuteronCuts( isMC, true,
      CombSigma, ContributionSplitting);
  AntiTrackCutsDeuteron->SetMinimalBooking(false);
  AntiTrackCutsDeuteron->SetCutCharge(-1);
/////////////////////For no NSigmaTOF information///
// =====================================================================
  //Proton track Cuts----------------------------------------------------------------------------
  AliFemtoDreamTrackCuts *TrackCutsNoTOF =
    AliFemtoDreamTrackCuts::PrimProtonCuts(isMC, true, CombSigma, ContributionSplitting);
  TrackCutsNoTOF->SetMinimalBooking(false);
  TrackCutsNoTOF->SetCutCharge(1);
  TrackCutsNoTOF->SetPID(AliPID::kProton, 999.);
  //Antiproton track Cuts-------------------------------------------------------------------------
  AliFemtoDreamTrackCuts *AntiTrackCutsNoTOF =
    AliFemtoDreamTrackCuts::PrimProtonCuts(isMC, true, CombSigma, ContributionSplitting);
  AntiTrackCutsNoTOF->SetMinimalBooking(false);
  AntiTrackCutsNoTOF->SetCutCharge(-1);
  AntiTrackCutsNoTOF-> SetPID(AliPID::kProton, 999.);
  //deuteron track cuts----------------------------------------------------------------------------
  AliFemtoDreamTrackCuts *TrackCutsDeuteronNoTOF = AliFemtoDreamTrackCuts::PrimDeuteronCuts(isMC, true,
      CombSigma, ContributionSplitting);
  TrackCutsDeuteronNoTOF->SetMinimalBooking(false);
  TrackCutsDeuteronNoTOF->SetCutCharge(1);
  TrackCutsDeuteronNoTOF->SetPID(AliPID::kDeuteron, 999.,3.,false,3.,true);
  //Antideuteron track cuts----------------------------------------------------------------------------
  AliFemtoDreamTrackCuts *AntiTrackCutsDeuteronNoTOF = AliFemtoDreamTrackCuts::PrimDeuteronCuts( isMC, true,
      CombSigma, ContributionSplitting);
  AntiTrackCutsDeuteronNoTOF->SetMinimalBooking(false);
  AntiTrackCutsDeuteronNoTOF->SetCutCharge(-1);
  AntiTrackCutsDeuteronNoTOF->SetPID(AliPID::kDeuteron, 999.,3.,false,3.,true);
//====================================================================================================================================
  std::vector<int> PDGParticles;
  PDGParticles.push_back(2212);
  PDGParticles.push_back(2212);
  PDGParticles.push_back(1000010020);
  PDGParticles.push_back(1000010020);

  std::vector<bool> closeRejection;
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
  }

  closeRejection[0] = true;  // pp
  closeRejection[2] = true;  // pd
  closeRejection[4] = true;  // barp barp
  closeRejection[6] = true;  // barp bard
  closeRejection[7] = true;  // dd
  closeRejection[9] = true;  // bard bard
  pairQA[0] = 11;    // pp
  pairQA[2] = 11;    // pd
  pairQA[4] = 11;    // barp barp
  pairQA[6] = 11;    // barp bard
  pairQA[7] = 11;    // dd
  pairQA[9] = 11;    // bard bard

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

  std::vector<int> NBins;
  NBins.push_back(750);
  NBins.push_back(750);
  NBins.push_back(750);
  NBins.push_back(750);
  NBins.push_back(750);
  NBins.push_back(750);
  NBins.push_back(750);
  NBins.push_back(750);
  NBins.push_back(750);
  NBins.push_back(750);

  std::vector<float> kMin;
  //minimum k* value
  kMin.push_back(0.);
  kMin.push_back(0.);
  kMin.push_back(0.);
  kMin.push_back(0.);
  kMin.push_back(0.);
  kMin.push_back(0.);
  kMin.push_back(0.);
  kMin.push_back(0.);
  kMin.push_back(0.);
  kMin.push_back(0.);
  //maximum k* value

  std::vector<float> kMax;
  kMax.push_back(3.);
  kMax.push_back(3.);
  kMax.push_back(3.);
  kMax.push_back(3.);
  kMax.push_back(3.);
  kMax.push_back(3.);
  kMax.push_back(3.);
  kMax.push_back(3.);
  kMax.push_back(3.);
  kMax.push_back(3.);

  AliFemtoDreamCollConfig *config = new AliFemtoDreamCollConfig("Femto", "Femto");
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
  config->SetDeltaEtaMax(0.012); // and here you set the actual values
  config->SetDeltaPhiMax(0.012); // and here you set the actual values
  //Here we set the mixing depth.
  config->SetMixingDepth(10);

  AliAnalysisTaskNanoPt *task =
    new AliAnalysisTaskNanoPt("AliAnalysisTaskNanoPt", isMC);

  if (trigger == "kINT7") {
    task->SelectCollisionCandidates(AliVEvent::kINT7);
    std::cout << "Added kINT7 Trigger \n";
  } else if (trigger == "kHM") {
    task->SelectCollisionCandidates(AliVEvent::kHighMultV0);
    std::cout << "Added kHighMult Trigger \n";
  } else {
    std::cout << "=====================================================================" << std::endl;
    std::cout << "=====================================================================" << std::endl;
    std::cout << "Centrality Estimator not set, fix it else your Results will be empty!" << std::endl;
    std::cout << "=====================================================================" << std::endl;
    std::cout << "=====================================================================" << std::endl;
  }

  task->SetEventCuts(evtCuts);
  task->SetProtonCuts(TrackCuts);
  task->SetAntiProtonCuts(AntiTrackCuts);
  task->SetDeuteronCuts(TrackCutsDeuteron);
  task->SetAntiDeuteronCuts(AntiTrackCutsDeuteron);
  task->SetProtonCutsNoTOF(TrackCutsNoTOF);
  task->SetAntiProtonCutsNoTOF(AntiTrackCutsNoTOF);
  task->SetDeuteronCutsNoTOF(TrackCutsDeuteronNoTOF);
  task->SetAntiDeuteronCutsNoTOF(AntiTrackCutsDeuteronNoTOF);
  task->SetCollectionConfig(config);
  mgr->AddTask(task);

  TString file = AliAnalysisManager::GetCommonFileName();
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  mgr->ConnectInput(task, 0, cinput);

  TString addon = "";

  if (trigger == "kInt7") {
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

  TString AntiTrackCutsName = Form("%sAntiProton%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *coutputAntiTrkCuts = mgr->CreateContainer(
        AntiTrackCutsName.Data(), TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), AntiTrackCutsName.Data()));
  mgr->ConnectOutput(task, 3, coutputAntiTrkCuts);

  TString TrackCutsDeuteronName = Form("%sDeuteron%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *coutputTrkCutsDeuteron = mgr->CreateContainer(
        TrackCutsDeuteronName.Data(), TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), TrackCutsDeuteronName.Data()));
  mgr->ConnectOutput(task, 4, coutputTrkCutsDeuteron);

  TString AntiTrackCutsDeuteronName = Form("%sAntiDeuteron%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *coutputAntiTrkCutsDeuteron = mgr->CreateContainer(
        AntiTrackCutsDeuteronName.Data(), TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), AntiTrackCutsDeuteronName.Data()));
  mgr->ConnectOutput(task, 5, coutputAntiTrkCutsDeuteron);
//============== NoTOF STUFF========================================
  TString TrackCutsNoTOFName = Form("%sProtonNoTOF%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *coutputTrkCutsNoTOF = mgr->CreateContainer(
        TrackCutsNoTOFName.Data(), TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), TrackCutsNoTOFName.Data()));
  mgr->ConnectOutput(task, 6, coutputTrkCutsNoTOF);

  TString AntiTrackCutsNoTOFName = Form("%sAntiProtonNoTOF%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *coutputAntiTrkCutsNoTOF = mgr->CreateContainer(
        AntiTrackCutsNoTOFName.Data(), TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), AntiTrackCutsNoTOFName.Data()));
  mgr->ConnectOutput(task, 7, coutputAntiTrkCutsNoTOF);

  TString TrackCutsDeuteronNoTOFName = Form("%sDeuteronNoTOF%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *coutputTrkCutsDeuteronNoTOF = mgr->CreateContainer(
        TrackCutsDeuteronNoTOFName.Data(), TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), TrackCutsDeuteronNoTOFName.Data()));
  mgr->ConnectOutput(task, 8, coutputTrkCutsDeuteronNoTOF);

  TString AntiTrackCutsDeuteronNoTOFName = Form("%sAntiDeuteronNoTOF%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *coutputAntiTrkCutsDeuteronNoTOF = mgr->CreateContainer(
        AntiTrackCutsDeuteronNoTOFName.Data(), TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), AntiTrackCutsDeuteronNoTOFName.Data()));
  mgr->ConnectOutput(task, 9, coutputAntiTrkCutsDeuteronNoTOF);


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

  if (isMC) {
    AliAnalysisDataContainer *coutputTrkCutsMC;
    TString TrkCutsMCName = Form("%sProtonMC%s", addon.Data(), suffix.Data());
    coutputTrkCutsMC = mgr->CreateContainer(
                         TrkCutsMCName.Data(),
                         TList::Class(),
                         AliAnalysisManager::kOutputContainer,
                         Form("%s:%s", file.Data(), TrkCutsMCName.Data()));
    mgr->ConnectOutput(task, 12, coutputTrkCutsMC);

    AliAnalysisDataContainer *coutputAntiTrkCutsMC;
    TString AntiTrkCutsMCName = Form("%sAntiProtonMC%s", addon.Data(),
                                     suffix.Data());
    coutputAntiTrkCutsMC = mgr->CreateContainer(
                             //@suppress("Invalid arguments") it works ffs
                             AntiTrkCutsMCName.Data(),
                             TList::Class(),
                             AliAnalysisManager::kOutputContainer,
                             Form("%s:%s", file.Data(), AntiTrkCutsMCName.Data()));
    mgr->ConnectOutput(task, 13, coutputAntiTrkCutsMC);

    AliAnalysisDataContainer *coutputv0CutsMC;
    TString v0CutsMCName = Form("%sDeuteronMC%s", addon.Data(), suffix.Data());
    coutputv0CutsMC = mgr->CreateContainer(
                        //@suppress("Invalid arguments") it works ffs
                        v0CutsMCName.Data(),
                        TList::Class(),
                        AliAnalysisManager::kOutputContainer,
                        Form("%s:%s", file.Data(), v0CutsMCName.Data()));
    mgr->ConnectOutput(task, 14, coutputv0CutsMC);

    AliAnalysisDataContainer *coutputAntiv0CutsMC;
    TString Antiv0CutsMCName = Form("%sAntiDeuteronMC%s", addon.Data(),
                                    suffix.Data());
    coutputAntiv0CutsMC = mgr->CreateContainer(
                            //@suppress("Invalid arguments") it works ffs
                            Antiv0CutsMCName.Data(),
                            TList::Class(),
                            AliAnalysisManager::kOutputContainer,
                            Form("%s:%s", file.Data(), Antiv0CutsMCName.Data()));
    mgr->ConnectOutput(task, 15, coutputAntiv0CutsMC);

  }

  return task;
}

