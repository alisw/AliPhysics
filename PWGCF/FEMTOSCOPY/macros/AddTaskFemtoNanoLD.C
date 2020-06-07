#include <vector>
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskNanoLD.h"
#include "AliFemtoDreamEventCuts.h"
#include "AliFemtoDreamTrackCuts.h"
#include "AliFemtoDreamCascadeCuts.h"
#include "AliFemtoDreamCollConfig.h"

AliAnalysisTaskSE *AddTaskFemtoNanoLD(bool fullBlastQA = false,
                                      bool isMC = false,
                                      bool Systematic = false,
                                      const char *cutVariation = "0") {

  TString suffix = TString::Format("%s", cutVariation);

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskFemtoNanoLD()", "No analysis manager found.");
    return 0x0;
  }

  // ================== GetInputEventHandler =============================
  AliVEventHandler *inputHandler = mgr->GetInputEventHandler();
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();

  //========= Init subtasks and start analyis ============================
  // Event Cuts
  AliFemtoDreamEventCuts *evtCuts = AliFemtoDreamEventCuts::StandardCutsRun2();
  evtCuts->CleanUpMult(false, false, false, true);

  // Track Cuts for Deuterons
  AliFemtoDreamTrackCuts *DeuteronCuts = new AliFemtoDreamTrackCuts();
  
  DeuteronCuts->SetIsMonteCarlo(isMC);
  DeuteronCuts->SetPlotDCADist(true);
  DeuteronCuts->SetPlotCombSigma(false);
  DeuteronCuts->SetPlotContrib(false);

  DeuteronCuts->SetFilterBit(128);
  DeuteronCuts->SetCutCharge(1);
  DeuteronCuts->SetPtRange(0.4, 4.);
  DeuteronCuts->SetEtaRange(-0.8, 0.8);
  DeuteronCuts->SetNClsTPC(80);
  DeuteronCuts->SetDCAReCalculation(true);  // Get the dca from PropagateToVertex
  DeuteronCuts->SetDCAVtxZ(0.2);
  DeuteronCuts->SetDCAVtxXY(0.1);
  DeuteronCuts->SetCutSharedCls(true);
  DeuteronCuts->SetCutTPCCrossedRows(true, 70, 0.83);
  DeuteronCuts->SetPID(AliPID::kDeuteron, 1.4);
  DeuteronCuts->SetRejLowPtPionsTOF(true);
  DeuteronCuts->SetCutSmallestSig(true);

  /*
  if (suffix == "1") {
    DeuteronCuts->SetFilterBit(256);
  }
  else if (suffix == "2") {
    DeuteronCuts->SetPID(AliPID::kDeuteron, 999.);
  }
  else if (suffix == "3") {
    DeuteronCuts->SetFilterBit(256);
    DeuteronCuts->SetPID(AliPID::kDeuteron, 999.);
    //DeuteronCuts->SetCutITSPID(1.4, -2., 1e30);
  }
  else if (suffix == "4") {
    DeuteronCuts->SetFilterBit(256);
    DeuteronCuts->SetDCAVtxZ(5.);
    DeuteronCuts->SetDCAVtxXY(5.);
  }
  */

  // Track Cuts for Anti-Deuterons
  AliFemtoDreamTrackCuts *AntiDeuteronCuts = new AliFemtoDreamTrackCuts();

  AntiDeuteronCuts->SetIsMonteCarlo(isMC);
  AntiDeuteronCuts->SetPlotDCADist(true);
  AntiDeuteronCuts->SetPlotCombSigma(false);
  AntiDeuteronCuts->SetPlotContrib(false);
  
  AntiDeuteronCuts->SetFilterBit(128);
  AntiDeuteronCuts->SetCutCharge(-1);
  AntiDeuteronCuts->SetPtRange(0.4, 4.);
  AntiDeuteronCuts->SetEtaRange(-0.8, 0.8);
  AntiDeuteronCuts->SetNClsTPC(80);
  AntiDeuteronCuts->SetDCAReCalculation(true);  // Get the dca from PropagateToVertex
  AntiDeuteronCuts->SetDCAVtxZ(0.2);
  AntiDeuteronCuts->SetDCAVtxXY(0.1);
  AntiDeuteronCuts->SetCutSharedCls(true);
  AntiDeuteronCuts->SetCutTPCCrossedRows(true, 70, 0.83);
  AntiDeuteronCuts->SetPID(AliPID::kDeuteron, 1.4);
  AntiDeuteronCuts->SetRejLowPtPionsTOF(true);
  AntiDeuteronCuts->SetCutSmallestSig(true);

  /*
  if (suffix == "1") {
    AntiDeuteronCuts->SetFilterBit(256);
  }
  else if (suffix == "2") {
    AntiDeuteronCuts->SetPID(AliPID::kDeuteron, 999.);
  }
  else if (suffix == "3") {
    AntiDeuteronCuts->SetFilterBit(256);
    AntiDeuteronCuts->SetPID(AliPID::kDeuteron, 999.);
    //AntiDeuteronCuts->SetCutITSPID(1.4, -2., 1e30);
  }
  else if (suffix == "4") {
    AntiDeuteronCuts->SetFilterBit(256);
    AntiDeuteronCuts->SetDCAVtxZ(5.);
    AntiDeuteronCuts->SetDCAVtxXY(5.);
  }
  */

  // Lambda Cuts
  AliFemtoDreamv0Cuts *v0Cuts = AliFemtoDreamv0Cuts::LambdaCuts(isMC, true,
                                                                false);
  v0Cuts->SetCutWindow(1.0775, 1.10768);
  if (suffix == "1") {
    v0Cuts->SetCutWindow(1.0775, 1.10568);
  }
  else if (suffix == "2") {
    v0Cuts->SetCutWindow(1.12368, 1.2);
  }
  else if (suffix == "3") {
    v0Cuts->SetCutWindow(1.12568, 1.2);
  }
  else if (suffix == "4") {
    v0Cuts->SetCutWindow(1.12368, 1.15386);
  }
  else if (suffix == "5") {
    v0Cuts->SetCutWindow(1.12568, 1.15386);
  }

  AliFemtoDreamTrackCuts *Posv0Daug = AliFemtoDreamTrackCuts::DecayProtonCuts(
      isMC, true, false);
  Posv0Daug->SetCutCharge(1);

  AliFemtoDreamTrackCuts *Negv0Daug = AliFemtoDreamTrackCuts::DecayPionCuts(
      isMC, true, false);
  Negv0Daug->SetCutCharge(-1);

  v0Cuts->SetPosDaugterTrackCuts(Posv0Daug);
  v0Cuts->SetNegDaugterTrackCuts(Negv0Daug);
  v0Cuts->SetPDGCodePosDaug(2212);  //Proton
  v0Cuts->SetPDGCodeNegDaug(211);  //Pion
  v0Cuts->SetPDGCodev0(3122);  //Lambda

  // Anti-Lambda Cuts
  AliFemtoDreamv0Cuts *Antiv0Cuts = AliFemtoDreamv0Cuts::LambdaCuts(isMC, true,
                                                                    false);
  Antiv0Cuts->SetCutWindow(1.0775, 1.10768);
  if (suffix == "1") {
    Antiv0Cuts->SetCutWindow(1.0775, 1.10568);
  }
  else if (suffix == "2") {
    Antiv0Cuts->SetCutWindow(1.12368, 1.2);
  }
  else if (suffix == "3") {
    Antiv0Cuts->SetCutWindow(1.12568, 1.2);
  }
  else if (suffix == "4") {
    Antiv0Cuts->SetCutWindow(1.12368, 1.15386);
  }
  else if (suffix == "5") {
    Antiv0Cuts->SetCutWindow(1.12568, 1.15386);
  }

  AliFemtoDreamTrackCuts *PosAntiv0Daug = AliFemtoDreamTrackCuts::DecayPionCuts(
      isMC, true, false);
  PosAntiv0Daug->SetCutCharge(1);

  AliFemtoDreamTrackCuts *NegAntiv0Daug = AliFemtoDreamTrackCuts::DecayProtonCuts(
        isMC, true, false);
  NegAntiv0Daug->SetCutCharge(-1);

  Antiv0Cuts->SetPosDaugterTrackCuts(PosAntiv0Daug);
  Antiv0Cuts->SetNegDaugterTrackCuts(NegAntiv0Daug);
  Antiv0Cuts->SetPDGCodePosDaug(211);  //Pion
  Antiv0Cuts->SetPDGCodeNegDaug(2212);  //Proton
  Antiv0Cuts->SetPDGCodev0(-3122);  //Lambda

  if (!fullBlastQA) {
    evtCuts->SetMinimalBooking(true);
    DeuteronCuts->SetMinimalBooking(true);
    AntiDeuteronCuts->SetMinimalBooking(true);
    v0Cuts->SetMinimalBooking(true);
    Antiv0Cuts->SetMinimalBooking(true);
  }

  AliFemtoDreamCollConfig *config = new AliFemtoDreamCollConfig("Femto",
                                                                "Femto", false);
  // Femto Collection
  std::vector<int> PDGParticles;
  PDGParticles.push_back(1000010020);
  PDGParticles.push_back(1000010020);
  PDGParticles.push_back(3122);
  PDGParticles.push_back(3122);

  std::vector<int> NBins;
  std::vector<float> kMin;
  std::vector<float> kMax;
  std::vector<int> pairQA;
  std::vector<bool> closeRejection;

  // Matrix of particles:
  //
  //       | d bar(d) L bar(L)
  //--------------------------
  //     d | 0   1    2   3
  // bar(d)|     4    5   6 
  //     L |          7   8
  // bar(L)|              9
  //
  // Corresponding pairs:
  // d d              0
  // d bar(d)         1
  // d L              2
  // d bar(L)         3
  // bar(d) bar(d)    4
  // bar(d) L         5
  // bar(d) bar(L)    6
  // L L              7
  // L bar(L)         8
  // bar(L) bar(L)    9

  const int nPairs = 10;
  for (int i = 0; i < nPairs; ++i) {
    pairQA.push_back(0);
    closeRejection.push_back(false);
    NBins.push_back(750);
    kMin.push_back(0.);
    kMax.push_back(3.);
  }
  pairQA[2] = 12;
  pairQA[6] = 12;

  config->SetPDGCodes(PDGParticles);
  config->SetNBinsHist(NBins);
  config->SetMinKRel(kMin);
  config->SetMaxKRel(kMax);
  config->SetClosePairRejection(closeRejection);
  //config->SetDeltaEtaMax(0.012);
  //config->SetDeltaPhiMax(0.012);
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
  config->SetmTBinning(true);

  config->SetdPhidEtaPlotsSmallK(false);
  config->SetdPhidEtaPlots(false);
  config->SetPhiEtaBinnign(false);

  if (fullBlastQA) {
    config->SetdPhidEtaPlotsSmallK(true);
    config->SetdPhidEtaPlots(true);
    config->SetPhiEtaBinnign(true);
    config->SetkTBinning(true);
    config->SetPtQA(true);
    //config->SetMassQA(true);
  }

  if (!fullBlastQA) {
    config->SetMinimalBookingME(true);
  }


  // Create the task
  AliAnalysisTaskNanoLD* task = new AliAnalysisTaskNanoLD("femtoLD");
  if (!fullBlastQA) {
    task->SetRunTaskLightWeight(true);
  }
  task->SelectCollisionCandidates(AliVEvent::kHighMultV0);
  task->SetEventCuts(evtCuts);
  task->SetDeuteronCuts(DeuteronCuts);
  task->SetAntiDeuteronCuts(AntiDeuteronCuts);
  task->Setv0Cuts(v0Cuts);
  task->SetAntiv0Cuts(Antiv0Cuts);
  task->SetCorrelationConfig(config);
  mgr->AddTask(task);

  TString addon = "LD";

  TString file = AliAnalysisManager::GetCommonFileName();

  mgr->ConnectInput(task, 0, cinput);

  // Output container for event cuts
  TString EvtCutsName = Form("%sEvtCuts%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *coutputEvtCuts = mgr->CreateContainer(
      EvtCutsName.Data(), TList::Class(), AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), EvtCutsName.Data()));
  mgr->ConnectOutput(task, 1, coutputEvtCuts);

  // Output container for Deuteron cuts
  TString DeuteronCutsName = Form("%sDeuteronCuts%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *couputTrkCuts = mgr->CreateContainer(
      DeuteronCutsName.Data(), TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), DeuteronCutsName.Data()));
  mgr->ConnectOutput(task, 2, couputTrkCuts);

  // Output container for Anti-Deuteron cuts
  TString AntiDeuteronCutsName = Form("%sAntiDeuteronCuts%s", addon.Data(),
                                   suffix.Data());
  AliAnalysisDataContainer *coutputAntiTrkCuts = mgr->CreateContainer(
      AntiDeuteronCutsName.Data(), TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), AntiDeuteronCutsName.Data()));
  mgr->ConnectOutput(task, 3, coutputAntiTrkCuts);

  // Output container for Lambda cuts
  TString v0CutsName = Form("%sv0Cuts%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *coutputv0Cuts = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      v0CutsName.Data(),
      TList::Class(), AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), v0CutsName.Data()));
   mgr->ConnectOutput(task, 4, coutputv0Cuts);
 
  // Output container for Anti-Lambda cuts
  TString Antiv0CutsName = Form("%sAntiv0Cuts%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *coutputAntiv0Cuts = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      Antiv0CutsName.Data(),
      TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), Antiv0CutsName.Data()));
   mgr->ConnectOutput(task, 5, coutputAntiv0Cuts);

  // Output container for results
  AliAnalysisDataContainer *coutputResults;
  TString ResultsName = Form("%sResults%s", addon.Data(), suffix.Data());
  coutputResults = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      ResultsName.Data(),
      TList::Class(), AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), ResultsName.Data()));
  mgr->ConnectOutput(task, 6, coutputResults);

  // Output container for results QA
  AliAnalysisDataContainer *coutputResultsQA;
  TString ResultsQAName = Form("%sResultsQA%s", addon.Data(), suffix.Data());
  coutputResultsQA = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      ResultsQAName.Data(),
      TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), ResultsQAName.Data()));
  mgr->ConnectOutput(task, 7, coutputResultsQA);

  return task;
}
