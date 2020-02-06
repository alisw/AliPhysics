#include "TROOT.h"
#include "TSystem.h"
AliAnalysisTaskSE* AddTaskFemtoDreamDeuteron(
  bool isMC = false,//1
  TString CentEst = "kInt7",//2
  bool DCAPlots = false,//3
  bool CombSigma = false,//4
  bool ContributionSplitting = false//5
) {

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
  AliFemtoDreamTrackCuts *TrackCutsDeuteronMass =  AliFemtoDreamTrackCuts::PrimDeuteronCuts(
        isMC, true, CombSigma, ContributionSplitting);
  TrackCutsDeuteronMass->SetCutCharge(1);
  TrackCutsDeuteronMass->SetPID(AliPID::kDeuteron, 999.);
  //anti deuterons
  AliFemtoDreamTrackCuts *TrackCutsAntiDeuteronDCA = AliFemtoDreamTrackCuts::PrimDeuteronCuts(
        isMC, true, CombSigma, ContributionSplitting);
  TrackCutsAntiDeuteronDCA->SetCutCharge(-1);
  AliFemtoDreamTrackCuts *TrackCutsAntiDeuteronMass =  AliFemtoDreamTrackCuts::PrimDeuteronCuts(
        isMC, true, CombSigma, ContributionSplitting);
  TrackCutsAntiDeuteronMass->SetCutCharge(-1);
  TrackCutsAntiDeuteronMass->SetPID(AliPID::kDeuteron, 999.);
//proton
  AliFemtoDreamTrackCuts *TrackCutsProtonDCA = AliFemtoDreamTrackCuts::PrimProtonCuts(
        isMC, true, CombSigma, ContributionSplitting);
  TrackCutsProtonDCA->SetCutCharge(1);
  AliFemtoDreamTrackCuts *TrackCutsProtonMass =  AliFemtoDreamTrackCuts::PrimProtonCuts(
        isMC, true, CombSigma, ContributionSplitting);
  TrackCutsProtonMass->SetCutCharge(1);
  TrackCutsProtonMass->SetPID(AliPID::kProton, 999.);
  //antiproton
  AliFemtoDreamTrackCuts *TrackCutsAntiProtonDCA = AliFemtoDreamTrackCuts::PrimProtonCuts(
        isMC, true, CombSigma, ContributionSplitting);
  TrackCutsAntiProtonDCA->SetCutCharge(-1);
  AliFemtoDreamTrackCuts *TrackCutsAntiProtonMass =  AliFemtoDreamTrackCuts::PrimProtonCuts(
        isMC, true, CombSigma, ContributionSplitting);
  TrackCutsAntiProtonMass->SetCutCharge(-1);
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

  AliAnalysisDataContainer *coutputQA;
  TString QAName = Form("%sQA", addon.Data());
  coutputQA = mgr->CreateContainer(
                QAName.Data(), TList::Class(),
                AliAnalysisManager::kOutputContainer,
                Form("%s:%s", file.Data(), QAName.Data()));
  mgr->ConnectOutput(task, 1, coutputQA);

  return task;
}

