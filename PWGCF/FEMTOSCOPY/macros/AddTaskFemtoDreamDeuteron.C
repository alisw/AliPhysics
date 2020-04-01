#include "TROOT.h"
#include "TSystem.h"
AliAnalysisTaskSE* AddTaskFemtoDreamDeuteron(bool isMC = false,//1
    TString CentEst = "kInt7",//2
    bool DCAPlots = false,//3
    bool CombSigma = false,//4
    bool ContributionSplitting = false,//5,
    bool DumpPdApAd = true,//6
    bool DumpRest = false,//7
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

  AliFemtoDreamEventCuts *evtCuts =
    AliFemtoDreamEventCuts::StandardCutsRun2();
  evtCuts->CleanUpMult(false, false, false, true);

  AliFemtoDreamTrackCuts *TrackCutsDeuteronDCA = AliFemtoDreamTrackCuts::PrimDeuteronCuts(
        isMC, true, CombSigma, ContributionSplitting);
  TrackCutsDeuteronDCA->SetCutCharge(1);
  TrackCutsDeuteronDCA->SetMinimalBooking(false);
  AliFemtoDreamTrackCuts *TrackCutsDeuteronMass =  AliFemtoDreamTrackCuts::PrimDeuteronCuts(
        isMC, true, CombSigma, ContributionSplitting);
  TrackCutsDeuteronMass->SetCutCharge(1);
  TrackCutsDeuteronMass->SetMinimalBooking(false);
  TrackCutsDeuteronMass->SetPID(AliPID::kDeuteron, 999.,60.);
  //anti deuterons
  AliFemtoDreamTrackCuts *TrackCutsAntiDeuteronDCA = AliFemtoDreamTrackCuts::PrimDeuteronCuts(
        isMC, true, CombSigma, ContributionSplitting);
  TrackCutsAntiDeuteronDCA->SetCutCharge(-1);
  TrackCutsAntiDeuteronDCA->SetMinimalBooking(false);
  AliFemtoDreamTrackCuts *TrackCutsAntiDeuteronMass =  AliFemtoDreamTrackCuts::PrimDeuteronCuts(
        isMC, true, CombSigma, ContributionSplitting);
  TrackCutsAntiDeuteronMass->SetCutCharge(-1);
  TrackCutsAntiDeuteronMass->SetMinimalBooking(false);
  TrackCutsAntiDeuteronMass->SetPID(AliPID::kDeuteron, 999.,60.);
  //proton
  AliFemtoDreamTrackCuts *TrackCutsProtonDCA = AliFemtoDreamTrackCuts::PrimProtonCuts(
        isMC, true, CombSigma, ContributionSplitting);
  TrackCutsProtonDCA->SetCutCharge(1);
  TrackCutsProtonDCA->SetMinimalBooking(false);
  AliFemtoDreamTrackCuts *TrackCutsProtonMass =  AliFemtoDreamTrackCuts::PrimProtonCuts(
        isMC, true, CombSigma, ContributionSplitting);
  TrackCutsProtonMass->SetCutCharge(1);
  TrackCutsProtonMass->SetMinimalBooking(false);
  TrackCutsProtonMass->SetPID(AliPID::kProton, 999.);
  //antiproton
  AliFemtoDreamTrackCuts *TrackCutsAntiProtonDCA = AliFemtoDreamTrackCuts::PrimProtonCuts(
        isMC, true, CombSigma, ContributionSplitting);
  TrackCutsAntiProtonDCA->SetCutCharge(-1);
  TrackCutsAntiProtonDCA->SetMinimalBooking(false);
  AliFemtoDreamTrackCuts *TrackCutsAntiProtonMass =  AliFemtoDreamTrackCuts::PrimProtonCuts(
        isMC, true, CombSigma, ContributionSplitting);
  TrackCutsAntiProtonMass->SetCutCharge(-1);
  TrackCutsAntiProtonMass->SetMinimalBooking(false);
  TrackCutsAntiProtonMass->SetPID(AliPID::kProton, 999.);

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
  //We need to set the ZVtx bins
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
  //The Multiplicity bins are set here
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
  config->SetMixingDepth(10);

  AliAnalysisTaskFemtoDreamDeuteron *task =
    new AliAnalysisTaskFemtoDreamDeuteron("FemtoDreamDefault", isMC);
  if (CentEst == "kInt7") {
    task->SelectCollisionCandidates(AliVEvent::kINT7);
    std::cout << "Added kINT7 Trigger \n";
  } else if (CentEst == "kHM") {
    task->SelectCollisionCandidates(AliVEvent::kHighMultV0);
    std::cout << "Added kHighMultV0 Trigger \n";
  } else {
    std::cout << "=====================================================================" << std::endl;
    std::cout << "=====================================================================" << std::endl;
    std::cout << "Centrality Estimator not set, fix it else your Results will be empty!" << std::endl;
    std::cout << "=====================================================================" << std::endl;
    std::cout << "=====================================================================" << std::endl;
  }

  //Throw all our settings to the task
  task->SetEventCuts(evtCuts);
  task->SetTrackCutsDeuteronDCA(TrackCutsDeuteronDCA);
  task->SetTrackCutsDeuteronMass(TrackCutsDeuteronMass);
  task->SetTrackCutsAntiDeuteronDCA(TrackCutsAntiDeuteronDCA);
  task->SetTrackCutsAntiDeuteronMass(TrackCutsAntiDeuteronMass);
  task->SetTrackCutsProtonDCA(TrackCutsProtonDCA);
  task->SetTrackCutsProtonMass(TrackCutsProtonMass);
  task->SetTrackCutsAntiProtonDCA(TrackCutsAntiProtonDCA);
  task->SetTrackCutsAntiProtonMass(TrackCutsAntiProtonMass);
  task->SetCollectionConfig(config);
  task->SetUseDumpster(DumpPdApAd);
  task->SetUseDumpsterRestPairs(DumpRest);
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

  TString EvtCutsName = Form("%sEvtCuts%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *coutputEvtCuts = mgr->CreateContainer(
        EvtCutsName.Data(), TList::Class(), AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), EvtCutsName.Data()));
  mgr->ConnectOutput(task, 1, coutputEvtCuts);

  TString TrackCutsName = Form("%sProtonDCA%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *coutputTrkCuts = mgr->CreateContainer(
        TrackCutsName.Data(), TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), TrackCutsName.Data()));
  mgr->ConnectOutput(task, 2, coutputTrkCuts);

  TString AntiTrackCutsName = Form("%sAntiProtonDCA%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *coutputAntiTrkCuts = mgr->CreateContainer(
        AntiTrackCutsName.Data(), TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), AntiTrackCutsName.Data()));
  mgr->ConnectOutput(task, 3, coutputAntiTrkCuts);

  TString TrackCutsDeuteronName = Form("%sDeuteronDCA%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *coutputTrkCutsDeuteron = mgr->CreateContainer(
        TrackCutsDeuteronName.Data(), TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), TrackCutsDeuteronName.Data()));
  mgr->ConnectOutput(task, 4, coutputTrkCutsDeuteron);

  TString AntiTrackCutsDeuteronName = Form("%sAntiDeuteronDCA%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *coutputAntiTrkCutsDeuteron = mgr->CreateContainer(
        AntiTrackCutsDeuteronName.Data(), TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), AntiTrackCutsDeuteronName.Data()));
  mgr->ConnectOutput(task, 5, coutputAntiTrkCutsDeuteron);
//============== NoTOF STUFF========================================
  TString TrackCutsNoTOFName = Form("%sProtonMass%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *coutputTrkCutsNoTOF = mgr->CreateContainer(
        TrackCutsNoTOFName.Data(), TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), TrackCutsNoTOFName.Data()));
  mgr->ConnectOutput(task, 6, coutputTrkCutsNoTOF);

  TString AntiTrackCutsNoTOFName = Form("%sAntiProtonMass%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *coutputAntiTrkCutsNoTOF = mgr->CreateContainer(
        AntiTrackCutsNoTOFName.Data(), TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), AntiTrackCutsNoTOFName.Data()));
  mgr->ConnectOutput(task, 7, coutputAntiTrkCutsNoTOF);

  TString TrackCutsDeuteronNoTOFName = Form("%sDeuteronMass%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *coutputTrkCutsDeuteronNoTOF = mgr->CreateContainer(
        TrackCutsDeuteronNoTOFName.Data(), TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), TrackCutsDeuteronNoTOFName.Data()));
  mgr->ConnectOutput(task, 8, coutputTrkCutsDeuteronNoTOF);

  TString AntiTrackCutsDeuteronNoTOFName = Form("%sAntiDeuteronMass%s", addon.Data(), suffix.Data());
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

  AliAnalysisDataContainer *coutputDumpster;
  TString DumpsterName = Form("%sDumpster%s", addon.Data(), suffix.Data());
  coutputDumpster = mgr->CreateContainer(
                      //@suppress("Invalid arguments") it works ffs
                      DumpsterName.Data(),
                      TList::Class(),
                      AliAnalysisManager::kOutputContainer,
                      Form("%s:%s", file.Data(), DumpsterName.Data()));
  mgr->ConnectOutput(task, 12, coutputDumpster);

  if (isMC) {
    AliAnalysisDataContainer *coutputTrkCutsMC;
    TString TrkCutsMCName = Form("%sProtonMC%s", addon.Data(), suffix.Data());
    coutputTrkCutsMC = mgr->CreateContainer(
                         TrkCutsMCName.Data(),
                         TList::Class(),
                         AliAnalysisManager::kOutputContainer,
                         Form("%s:%s", file.Data(), TrkCutsMCName.Data()));
    mgr->ConnectOutput(task, 13, coutputTrkCutsMC);

    AliAnalysisDataContainer *coutputAntiTrkCutsMC;
    TString AntiTrkCutsMCName = Form("%sAntiProtonMC%s", addon.Data(),
                                     suffix.Data());
    coutputAntiTrkCutsMC = mgr->CreateContainer(
                             //@suppress("Invalid arguments") it works ffs
                             AntiTrkCutsMCName.Data(),
                             TList::Class(),
                             AliAnalysisManager::kOutputContainer,
                             Form("%s:%s", file.Data(), AntiTrkCutsMCName.Data()));
    mgr->ConnectOutput(task, 14, coutputAntiTrkCutsMC);

    AliAnalysisDataContainer *coutputv0CutsMC;
    TString v0CutsMCName = Form("%sDeuteronMC%s", addon.Data(), suffix.Data());
    coutputv0CutsMC = mgr->CreateContainer(
                        //@suppress("Invalid arguments") it works ffs
                        v0CutsMCName.Data(),
                        TList::Class(),
                        AliAnalysisManager::kOutputContainer,
                        Form("%s:%s", file.Data(), v0CutsMCName.Data()));
    mgr->ConnectOutput(task, 15, coutputv0CutsMC);

    AliAnalysisDataContainer *coutputAntiv0CutsMC;
    TString Antiv0CutsMCName = Form("%sAntiDeuteronMC%s", addon.Data(),
                                    suffix.Data());
    coutputAntiv0CutsMC = mgr->CreateContainer(
                            //@suppress("Invalid arguments") it works ffs
                            Antiv0CutsMCName.Data(),
                            TList::Class(),
                            AliAnalysisManager::kOutputContainer,
                            Form("%s:%s", file.Data(), Antiv0CutsMCName.Data()));
    mgr->ConnectOutput(task, 16, coutputAntiv0CutsMC);
  }
  return task;
}

