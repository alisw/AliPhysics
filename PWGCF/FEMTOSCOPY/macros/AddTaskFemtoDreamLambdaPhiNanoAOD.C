#ifndef __CINT__
#include <vector>
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskNanoAODFemtoDreamLambdaPhi.h"
#include "AliFemtoDreamEventCuts.h"
#include "AliFemtoDreamTrackCuts.h"
#include "AliFemtoDreamv0Cuts.h"
#include "AliFemtoDreamCollConfig.h"
#endif

AliAnalysisTaskSE *AddTaskFemtoDreamLambdaPhiNanoAOD(bool isMC = false,
                                               TString CentEst = "kInt7",
                                               const char* sTcut = "0",
                                               double maxDCAv0vtx = 1.5,
                                               double minDCADaughPV = 0.05,
                                               const char *cutVariation = "0") {
  TString suffix = TString::Format("%s", cutVariation);
  TString sTsuffix = TString::Format("%s", sTcut);

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  if (!mgr) {
    printf("No analysis manager to connect to!\n");
    return nullptr;
  }
  if (!mgr->GetInputEventHandler()) {
    printf("This task requires an input event handler!\n");
    return nullptr;
  }

  AliFemtoDreamEventCuts *evtCuts = AliFemtoDreamEventCuts::StandardCutsRun2();
  evtCuts->CleanUpMult(false, false, false, true);
  if(sTsuffix=="1") {
    evtCuts->SetSphericityCuts(0.7, 1);
  }

  //Lambda Cuts
  AliFemtoDreamv0Cuts *v0Cuts = AliFemtoDreamv0Cuts::LambdaCuts(isMC, true, true);
  AliFemtoDreamTrackCuts *Posv0Daug = AliFemtoDreamTrackCuts::DecayProtonCuts(isMC, true, false); //PileUpRej, false
  AliFemtoDreamTrackCuts *Negv0Daug = AliFemtoDreamTrackCuts::DecayPionCuts(isMC, true, false);
  //Loosening cuts for Lambdas
  //1. DCA to V0 vtx (def. < 1.5 cm)
  //2. DCA daugh to PV (def > 0.05 cm)
  v0Cuts->SetCutDCADaugTov0Vtx(maxDCAv0vtx);
  v0Cuts->SetCutDCADaugToPrimVtx(minDCADaughPV)

  v0Cuts->SetPosDaugterTrackCuts(Posv0Daug);
  v0Cuts->SetNegDaugterTrackCuts(Negv0Daug);
  v0Cuts->SetPDGCodePosDaug(2212); //Proton
  v0Cuts->SetPDGCodeNegDaug(211);  //Pion
  v0Cuts->SetPDGCodev0(3122);      //Lambda

  AliFemtoDreamv0Cuts *Antiv0Cuts = AliFemtoDreamv0Cuts::LambdaCuts(isMC, true, true);
  AliFemtoDreamTrackCuts *PosAntiv0Daug = AliFemtoDreamTrackCuts::DecayPionCuts(isMC, true, false);
  PosAntiv0Daug->SetCutCharge(1);
  AliFemtoDreamTrackCuts *NegAntiv0Daug = AliFemtoDreamTrackCuts::DecayProtonCuts(isMC, true, false);
  NegAntiv0Daug->SetCutCharge(-1);
  //Loosening cuts for AntiLambdas
  //1. DCA to V0 vtx (def. < 1.5 cm)
  //2. DCA daugh to PV (def > 0.05 cm)
  Antiv0Cuts->SetCutDCADaugTov0Vtx(maxDCAv0vtx);
  Antiv0Cuts->SetCutDCADaugToPrimVtx(minDCADaughPV);

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

  if (suffix != "0" && suffix != "999") {
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

  if (suffix != "0") {
    TrackCutsPhi->SetMinimalBooking(true);
  }


  // Now we define stuff we want for our Particle collection
  // Thanks, CINT - will not compile due to an illegal constructor
  // std::vector<int> PDGParticles ={2212,2212,3122,3122,3312,3312};
  // First we need to tell him about the particles we mix, from the
  // PDG code the mass is obtained.
  std::vector<int> PDGParticles;
  PDGParticles.push_back(3122);  // 0 Lambda
  PDGParticles.push_back(3122);  // 1 antiLambda
  PDGParticles.push_back(333);   // 2 Phi particle
  PDGParticles.push_back(321);   // 3 Kaon Plus
  PDGParticles.push_back(321);   // 4 Kaon Minus


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
  //  NBins.push_back(750);
  //  NBins.push_back(750);
  //  NBins.push_back(750);
  //  NBins.push_back(750);
  //  NBins.push_back(750);
  //  NBins.push_back(750);
  std::vector<float> kMin;
  // minimum k* value
  //  kMin.push_back(0.);
  //  kMin.push_back(0.);
  //  kMin.push_back(0.);
  //  kMin.push_back(0.);
  //  kMin.push_back(0.);
  //  kMin.push_back(0.);
  // maximum k* value
  std::vector<float> kMax;
  //  kMax.push_back(3.);
  //  kMax.push_back(3.);
  //  kMax.push_back(3.);
  //  kMax.push_back(3.);
  //  kMax.push_back(3.);
  //  kMax.push_back(3.);
  // pair QA extended
  std::vector<int> pairQA;
  //  pairQA.push_back(11); //pp
  //  pairQA.push_back(0); //pap
  //  pairQA.push_back(12); //pphi
  //  pairQA.push_back(11); //apap
  //  pairQA.push_back(12); //apphi
  //  pairQA.push_back(22); //phiphi

  for (int i = 0; i < (1 + 2 + 3 + 4 + 5); i++) {
    NBins.push_back(750);
    kMin.push_back(0.);
    kMax.push_back(3.);
    pairQA.push_back(0);
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

  if (isMC) {
    config->SetMomentumResolution(true);
  } else {
    std::cout << "You are trying to request the Momentum Resolution without MC Info; fix it wont work! \n";
  }

  // now we create the task
  AliAnalysisTaskNanoAODFemtoDreamLambdaPhi *task =
      new AliAnalysisTaskNanoAODFemtoDreamLambdaPhi(
          "AliAnalysisTaskNanoAODFemtoDreamLambdaPhi", isMC);
  // THIS IS VERY IMPORTANT ELSE YOU DONT PROCESS ANY EVENTS
  // kINT7 == Minimum bias
  // kHighMultV0 high multiplicity triggered by the V0 detector
  if (CentEst == "kInt7") {
    task->SetTrigger(AliVEvent::kINT7);
    task->SelectCollisionCandidates(AliVEvent::kINT7);
    std::cout << "Added kINT7 Trigger \n";
  } else if (CentEst == "kHM") {
    task->SetTrigger(AliVEvent::kHighMultV0);
    task->SelectCollisionCandidates(AliVEvent::kHighMultV0);
    std::cout << "Added kHighMultV0 Trigger \n";
  } else {
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
  if (CentEst == "kInt7") {
    addon += "MB";
  } else if (CentEst == "kHM") {
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

  // TString QAName = Form("%sLambdaPhiResults%s", addon.Data(), suffix.Data());
  // coutputQA = mgr->CreateContainer(QAName.Data(), TList::Class(),
  //                                  AliAnalysisManager::kOutputContainer,
  //                                  Form("%s:%s", file.Data(), QAName.Data()));
  // mgr->ConnectOutput(task, 1, coutputQA);

  return task;
}
