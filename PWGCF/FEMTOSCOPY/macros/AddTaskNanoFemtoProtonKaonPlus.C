#include <vector>
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliFemtoDreamEventCuts.h"
#include "AliFemtoDreamTrackCuts.h"
#include "AliFemtoDreamCollConfig.h"
#include "AliAnalysisTaskNanoFemtoProtonKaonPlus.h"

AliAnalysisTaskSE* AddTaskNanoFemtoProtonKaonPlus(
    bool isMC = false, //1
    TString trigger = "kHM", //2
    bool fullBlastQA = true,//3
    bool UseSphericityCut = true,//4
    float SphericityMinPt = 0.5, //5
    bool UseRamonaCut = false, //6
    int filterBit = 128, //7
    bool DoPairCleaning = false, //8
    const char *cutVariation = "0", //9
    bool DoAncestors = false //10
    ) {

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

  //Event cuts ------------------------------------------------------------------------
  AliFemtoDreamEventCuts *evtCuts = AliFemtoDreamEventCuts::StandardCutsRun2();
  evtCuts->CleanUpMult(false, false, false, true);

  if (UseSphericityCut){
    float SpherDown = 0.7;
    evtCuts->SetSphericityCuts(SpherDown, 1.0, SphericityMinPt); 
  }


  //Track cuts ------------------------------------------------------------------------
  
  AliFemtoDreamTrackCuts *TrackPosKaonCuts =
      AliFemtoDreamTrackCuts::PrimKaonCuts(isMC, true, false, false);
  TrackPosKaonCuts->SetCutCharge(1);
  TrackPosKaonCuts->SetFilterBit(filterBit);
  if (!UseRamonaCut)
  {
    TrackPosKaonCuts->SetPIDkd(); // Oton
  }
  else
  {
    TrackPosKaonCuts->SetPIDkd(true, true); // Ramona
  }

  AliFemtoDreamTrackCuts *TrackNegKaonCuts =
      AliFemtoDreamTrackCuts::PrimKaonCuts(isMC, true, false, false);
  TrackNegKaonCuts->SetCutCharge(-1);
  TrackNegKaonCuts->SetFilterBit(filterBit);
  if (!UseRamonaCut)
  {
    TrackNegKaonCuts->SetPIDkd(); // Oton
  }
  else
  {
    TrackNegKaonCuts->SetPIDkd(true, true); // Ramona
  }

  //Proton and AntiProton cuts
  AliFemtoDreamTrackCuts *TrackCutsProton = AliFemtoDreamTrackCuts::PrimProtonCuts(
        isMC, true, false, false);
  TrackCutsProton->SetFilterBit(128);
  TrackCutsProton->SetCutCharge(1);

  AliFemtoDreamTrackCuts *TrackCutsAntiProton = AliFemtoDreamTrackCuts::PrimProtonCuts(
        isMC, true, false, false);
  TrackCutsAntiProton->SetFilterBit(128);
  TrackCutsAntiProton->SetCutCharge(-1);


  //Set-up output ------------------------------------------------------------------------
  std::vector<int> PDGParticles;
  PDGParticles.push_back(2212); 
  PDGParticles.push_back(-2212); 
  PDGParticles.push_back(321); 
  PDGParticles.push_back(-321);


  //BRELOOM Will have to enter mT bins
  std::vector<bool> closeRejection;
  std::vector<float> mTBins;
  mTBins.push_back(0.71); 
  mTBins.push_back(1.); 
  mTBins.push_back(1.2); 
  mTBins.push_back(1.4); 
  mTBins.push_back(1.7); 
  mTBins.push_back(1.9); 
  mTBins.push_back(10.0); 
  std::vector<int> pairQA;
  //pairs: 
  // pp             0
  // p bar p        1
  // p K+          2
  // p K-          3
  // bar p bar p    4
  // bar p K+      5
  // bar p K-      6
  // K+ K+        7
  // K+ K-        8
  // K- K-        9
  const int nPairs = 10;
  std::vector<int> NBins;
  std::vector<float> kMin;   //minimum k* value
  std::vector<float> kMax; //maximum k* value

  for (int i = 0; i < nPairs; ++i) {
    closeRejection.push_back(false);
    pairQA.push_back(0);
    NBins.push_back(750);
    kMin.push_back(0.);
    kMax.push_back(3.);
  }

  closeRejection[0] = true;  // pp
  closeRejection[2] = true;  // pK+
  closeRejection[3] = false;  // pK-
  closeRejection[4] = true;  // barp barp
  closeRejection[5] = false;  // barp K+
  closeRejection[6] = true;  // barp K-
  closeRejection[7] = true;  // K+K+
  closeRejection[9] = true;  // K-K-

  pairQA[0] = 11;  // pp
  pairQA[2] = 11;  // pK+
  pairQA[3] = 11;  // pK-
  pairQA[4] = 11;  // barp barp
  pairQA[5] = 11;  // barp K+
  pairQA[6] = 11;  // barp K-
  pairQA[7] = 11;  // K+K+
  pairQA[9] = 11;  // K-K-

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
  MultBins.push_back(84);
  MultBins.push_back(88);
  MultBins.push_back(92);
  MultBins.push_back(96);
  MultBins.push_back(100);

  AliFemtoDreamCollConfig *config = new AliFemtoDreamCollConfig("Femto", "Femto", false);
  config->SetZBins(ZVtxBins);
  config->SetMultBins(MultBins);
  config->SetMultBinning(true);
  config->SetPDGCodes(PDGParticles);
  config->SetNBinsHist(NBins);
  config->SetMinKRel(kMin);
  config->SetMaxKRel(kMax);
  config->SetClosePairRejection(closeRejection);
  config->SetExtendedQAPairs(pairQA);
  //BRELOOM Will have to set values here
  config->SetDeltaEtaMax(0.017); // and here you set the actual values 
  config->SetDeltaPhiMax(0.04); // and here you set the actual values 
  config->SetMixingDepth(10);
  config->SetmTBins(mTBins);
  config->SetDomTMultBinning(true);
  config->SetmTBinning(true);
  config->SetMultiplicityEstimator(AliFemtoDreamEvent::kRef08);

  if (isMC) {
    config->SetMomentumResolution(true);
    if(DoAncestors){
      config->SetAncestors(true);
      config->GetDoAncestorsPlots();
    }
  }
  if (fullBlastQA) {
    // config->SetkTBinning(true);
    config->SetPtQA(true);
    config->SetdPhidEtaPlots(true);
    config->SetdPhidEtaPlotsSmallK(true);
    config->SetPhiEtaBinnign(true);
  } else {
    evtCuts->SetMinimalBooking(true);
    TrackCutsProton->SetMinimalBooking(true);
    TrackCutsAntiProton->SetMinimalBooking(true);
    TrackPosKaonCuts->SetMinimalBooking(true);
    TrackNegKaonCuts->SetMinimalBooking(true);
    config->SetMinimalBookingME(true);
    config->SetMinimalBookingSample(true);
  }

  AliAnalysisTaskNanoFemtoProtonKaonPlus *task =
  new AliAnalysisTaskNanoFemtoProtonKaonPlus("FemtoDreamDefault", isMC);
  if (trigger == "kINT7") {
    task->SelectCollisionCandidates(AliVEvent::kINT7);
    std::cout << "Added kINT7 Trigger \n";
  } else if (trigger == "kHM") {
    task->SelectCollisionCandidates(AliVEvent::kHighMultV0);
    std::cout << "Added kHighMult Trigger \n";
  } else {
    std::cout
        << "====================================================================="
        << std::endl;
    std::cout
        << "====================================================================="
        << std::endl;
    std::cout
        << "Centrality Estimator not set, fix it else your Results will be empty!"
        << std::endl;
    std::cout
        << "====================================================================="
        << std::endl;
    std::cout
        << "====================================================================="
        << std::endl;
  }
  if (!fullBlastQA) {
    task->SetRunTaskLightWeight(true);
  }
  task->SetEventCuts(evtCuts);
  task->SetTrackCutsKaon(TrackPosKaonCuts);
  task->SetTrackCutsAntiKaon(TrackNegKaonCuts);
  task->SetTrackCutsProton(TrackCutsProton);
  task->SetTrackCutsAntiProton(TrackCutsAntiProton);
  task->SetCollectionConfig(config);
  task->SetDoPairCleaning(DoPairCleaning);

  mgr->AddTask(task);

  TString addon = "";

  if (trigger == "kINT7") {
    addon += "kINT7";
  } else if (trigger == "kHM") {
    addon += "HM";
  }

  suffix = ""; 
  TString file = AliAnalysisManager::GetCommonFileName();
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  mgr->ConnectInput(task, 0, cinput);

  TString EvtCutsName = Form("%s_EvtCuts_%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *coutputEvtCuts = mgr->CreateContainer(
        EvtCutsName.Data(), TList::Class(), AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), EvtCutsName.Data()));
  mgr->ConnectOutput(task, 1, coutputEvtCuts);

  TString TrackCutsName = Form("%s_TrackProton_%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *coutputTrkCuts = mgr->CreateContainer(
        TrackCutsName.Data(), TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), TrackCutsName.Data()));
  mgr->ConnectOutput(task, 2, coutputTrkCuts);

  TString AntiTrackCutsName = Form("%s_AntiTrackProton_%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *coutputAntiTrkCuts = mgr->CreateContainer(
        AntiTrackCutsName.Data(), TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), AntiTrackCutsName.Data()));
  mgr->ConnectOutput(task, 3, coutputAntiTrkCuts);

//BRELOOM change to kaon
  TString TrackCutsKaonName = Form("%s_TrackKaon_%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *coutputTrkCutsKaon = mgr->CreateContainer(
        TrackCutsKaonName.Data(), TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), TrackCutsKaonName.Data()));
  mgr->ConnectOutput(task, 4, coutputTrkCutsKaon);

  TString AntiTrackCutsKaonName = Form("%s_AntiTrackKaon_%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *coutputAntiTrkCutsKaon = mgr->CreateContainer(
        AntiTrackCutsKaonName.Data(), TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), AntiTrackCutsKaonName.Data()));
  mgr->ConnectOutput(task, 5, coutputAntiTrkCutsKaon);

  AliAnalysisDataContainer *coutputResults;
  TString ResultsName = Form("%s_Results_%s", addon.Data(), suffix.Data());
  coutputResults = mgr->CreateContainer(ResultsName.Data(),
                     TList::Class(), AliAnalysisManager::kOutputContainer,
                     Form("%s:%s", file.Data(), ResultsName.Data()));
  mgr->ConnectOutput(task, 6, coutputResults);

  AliAnalysisDataContainer *coutputResultsQA;
  TString ResultsQAName = Form("%s_ResultsQA_%s", addon.Data(), suffix.Data());
  coutputResultsQA = mgr->CreateContainer(
                       //@suppress("Invalid arguments") it works ffs
                       ResultsQAName.Data(),
                       TList::Class(),
                       AliAnalysisManager::kOutputContainer,
                       Form("%s:%s", file.Data(), ResultsQAName.Data()));
  mgr->ConnectOutput(task, 7, coutputResultsQA);

  if (isMC) {
    AliAnalysisDataContainer *coutputTrkCutsMC;
    TString TrkCutsMCName = Form("%s_ProtonMC_%s", addon.Data(), suffix.Data());
    coutputTrkCutsMC = mgr->CreateContainer(
                         TrkCutsMCName.Data(),
                         TList::Class(),
                         AliAnalysisManager::kOutputContainer,
                         Form("%s:%s", file.Data(), TrkCutsMCName.Data()));
    mgr->ConnectOutput(task, 8, coutputTrkCutsMC);

    AliAnalysisDataContainer *coutputAntiTrkCutsMC;
    TString AntiTrkCutsMCName = Form("%s_AntiProtonMC_%s", addon.Data(),
                                     suffix.Data());
    coutputAntiTrkCutsMC = mgr->CreateContainer(
                             //@suppress("Invalid arguments") it works ffs
                             AntiTrkCutsMCName.Data(),
                             TList::Class(),
                             AliAnalysisManager::kOutputContainer,
                             Form("%s:%s", file.Data(), AntiTrkCutsMCName.Data()));
    mgr->ConnectOutput(task, 9, coutputAntiTrkCutsMC);

    AliAnalysisDataContainer *coutputv0CutsMC;
    TString v0CutsMCName = Form("%s_KaonMC_%s", addon.Data(), suffix.Data());
    coutputv0CutsMC = mgr->CreateContainer(
                        //@suppress("Invalid arguments") it works ffs
                        v0CutsMCName.Data(),
                        TList::Class(),
                        AliAnalysisManager::kOutputContainer,
                        Form("%s:%s", file.Data(), v0CutsMCName.Data()));
    mgr->ConnectOutput(task, 10, coutputv0CutsMC);

    AliAnalysisDataContainer *coutputAntiv0CutsMC;
    TString Antiv0CutsMCName = Form("%s_AntiKaonMC_%s", addon.Data(),
                                    suffix.Data());
    coutputAntiv0CutsMC = mgr->CreateContainer(
                            //@suppress("Invalid arguments") it works ffs
                            Antiv0CutsMCName.Data(),
                            TList::Class(),
                            AliAnalysisManager::kOutputContainer,
                            Form("%s:%s", file.Data(), Antiv0CutsMCName.Data()));
    mgr->ConnectOutput(task, 11, coutputAntiv0CutsMC);
  }
  return task;
}

