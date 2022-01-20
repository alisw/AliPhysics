#if !defined(__CINT__) || defined(__CLING__)
#include <vector>
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskThreeBodyFemtoMixedCharge.h"
#include "AliAnalysisTaskThreeBodyFemtoAOD.h"
#include "AliFemtoDreamEventCuts.h"
#include "AliFemtoDreamTrackCuts.h"
#include "AliFemtoDreamCascadeCuts.h"
#include "AliFemtoDreamCollConfig.h"
#endif

AliAnalysisTaskSE *AddTaskThreeBodyFemtoMixedCharge(int trigger = 0, bool fullBlastQA = true,
                                     bool isMC = false,
                                     int mixingDepthFromTask = 20,
                                     const char *cutVariation = "0", bool ClosePairRejectionForAll = false, 
                                     bool run2Body = false, int mixinfChoice = 0, bool mix21 = false,
                                     int whichTripletsToRun = 11,
                                     const char *triggerVariation = "0") {



  TString suffix = TString::Format("%s", cutVariation);
  TString suffixTrigger = TString::Format("%s", triggerVariation);
  bool isNano = true;

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
      isMC, true, false, false);
  TrackCuts->SetFilterBit(128);
  TrackCuts->SetCutCharge(1);

  AliFemtoDreamTrackCuts *AntiTrackCuts =
      AliFemtoDreamTrackCuts::PrimProtonCuts(isMC, true, false, false);
  AntiTrackCuts->SetFilterBit(128);
  AntiTrackCuts->SetCutCharge(-1);
  // If  cut variation needed
  if(suffix=="1" || suffix=="8"){
    TrackCuts->SetEtaRange(-0.9, 0.9);
    AntiTrackCuts->SetEtaRange(-0.9, 0.9);
  }
  //Lambda Cuts
  AliFemtoDreamv0Cuts *v0Cuts = AliFemtoDreamv0Cuts::LambdaCuts(isMC, true,
                                                                false);
  AliFemtoDreamTrackCuts *Posv0Daug = AliFemtoDreamTrackCuts::DecayProtonCuts(
      isMC, true, false);

  AliFemtoDreamTrackCuts *Negv0Daug = AliFemtoDreamTrackCuts::DecayPionCuts(
      isMC, true, false);
  // If  cut variation needed
  if(suffix=="2" || suffix=="7" || suffix=="8"){
    Posv0Daug->SetEtaRange(-0.9, 0.9);
    Negv0Daug->SetEtaRange(-0.9, 0.9);
  }
  v0Cuts->SetPosDaugterTrackCuts(Posv0Daug);
  v0Cuts->SetNegDaugterTrackCuts(Negv0Daug);
  v0Cuts->SetPDGCodePosDaug(2212);  //Proton
  v0Cuts->SetPDGCodeNegDaug(211);  //Pion
  v0Cuts->SetPDGCodev0(3122);  //Lambda

  AliFemtoDreamv0Cuts *Antiv0Cuts = AliFemtoDreamv0Cuts::LambdaCuts(isMC, true,
                                                                    false);
  AliFemtoDreamTrackCuts *PosAntiv0Daug = AliFemtoDreamTrackCuts::DecayPionCuts(
      isMC, true, false);
  PosAntiv0Daug->SetCutCharge(1);
  AliFemtoDreamTrackCuts *NegAntiv0Daug =
      AliFemtoDreamTrackCuts::DecayProtonCuts(isMC, true, false);
  NegAntiv0Daug->SetCutCharge(-1);
  // If  cut variation needed
  if(suffix=="2" || suffix=="7" || suffix=="8"){
    PosAntiv0Daug->SetEtaRange(-0.9, 0.9);
    NegAntiv0Daug->SetEtaRange(-0.9, 0.9);
  }
  Antiv0Cuts->SetPosDaugterTrackCuts(PosAntiv0Daug);
  Antiv0Cuts->SetNegDaugterTrackCuts(NegAntiv0Daug);
  Antiv0Cuts->SetPDGCodePosDaug(211);  //Pion
  Antiv0Cuts->SetPDGCodeNegDaug(2212);  //Proton
  Antiv0Cuts->SetPDGCodev0(-3122);  //Lambda

  // If  cut variations eta
  if(suffix=="3" || suffix=="6" || suffix=="7" || suffix=="8"){
    v0Cuts->SetCutInvMass(0.006);
    Antiv0Cuts->SetCutInvMass(0.006);
  }
  if(suffix=="4"|| suffix=="6" || suffix=="7" || suffix=="8"){
    v0Cuts->SetCutCPA(0.999);
    Antiv0Cuts->SetCutCPA(0.999);
  }
  if(suffix=="5"|| suffix=="6" || suffix=="7" || suffix=="8"){
    v0Cuts->SetDaughterTimingCut(AliFemtoDreamv0Cuts::OneDaughterCombined);
    Antiv0Cuts->SetDaughterTimingCut(AliFemtoDreamv0Cuts::OneDaughterCombined);
  }


  if (!fullBlastQA) {
    evtCuts->SetMinimalBooking(true);
    TrackCuts->SetMinimalBooking(true);
    AntiTrackCuts->SetMinimalBooking(true);
    v0Cuts->SetMinimalBooking(true);
    Antiv0Cuts->SetMinimalBooking(true);
  }


  AliFemtoDreamEventCuts *evtCutsTrigger;
  AliFemtoDreamTrackCuts *TrackCutsTrigger;
  AliFemtoDreamTrackCuts *AntiTrackCutsTrigger;
  AliFemtoDreamv0Cuts *v0CutsTrigger;
  AliFemtoDreamTrackCuts *Posv0DaugTrigger;
  AliFemtoDreamTrackCuts *Negv0DaugTrigger;
  AliFemtoDreamv0Cuts *Antiv0CutsTrigger;
  AliFemtoDreamTrackCuts *PosAntiv0DaugTrigger;
  AliFemtoDreamTrackCuts *NegAntiv0DaugTrigger;

  bool TriggerOnSample = false;

  AliFemtoDreamCollConfig *config = new AliFemtoDreamCollConfig("Femto",
                                                                "Femto", false);
  // Femto Collection
  std::vector<int> PDGParticles;
  PDGParticles.push_back(2212);
  PDGParticles.push_back(2212);
  PDGParticles.push_back(3122);
  PDGParticles.push_back(3122);

  std::vector<int> NBins;
  std::vector<float> kMin;
  std::vector<float> kMax;
  std::vector<int> pairQA;
  std::vector<bool> closeRejection;

  const int nPairs = 10;
  for (int i = 0; i < nPairs; ++i) {
    pairQA.push_back(0);
    closeRejection.push_back(false);
    NBins.push_back(1000);
    kMin.push_back(0.);
    kMax.push_back(1.);
  }
  pairQA[0] = 11;
  pairQA[4] = 11;
  pairQA[2] = 12;
  pairQA[6] = 12;

  closeRejection[0] = true;  // pp
  closeRejection[4] = true;  // barp barp

  config->SetPDGCodes(PDGParticles);
  config->SetNBinsHist(NBins);
  config->SetMinKRel(kMin);
  config->SetMaxKRel(kMax);
  config->SetClosePairRejection(closeRejection);
  config->SetDeltaEtaMax(0.017);
  config->SetDeltaPhiMax(0.017);

  if(suffixTrigger=="699"){
    config->SetDeltaEtaMax(0.02);
    config->SetDeltaPhiMax(0.02);
  }  
  if(suffixTrigger=="669"){
    config->SetDeltaEtaMax(0.03);
    config->SetDeltaPhiMax(0.03);
  }
  config->SetExtendedQAPairs(pairQA);
  config->SetMixingDepth(mixingDepthFromTask);
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

  if (isMC) {
    config->SetMomentumResolution(true);
  }

  if (fullBlastQA) {
    config->SetkTBinning(true);
    config->SetPtQA(true);
    config->SetMassQA(true);
  }

  if (!fullBlastQA) {
    config->SetMinimalBookingME(true);
    config->SetMinimalBookingSample(true);
  }

  TString addon = "PL";
  TString file = AliAnalysisManager::GetCommonFileName();

  TString EvtCutsName = Form("%sEvtCuts_%s_%s", addon.Data(), suffixTrigger.Data(), suffix.Data());
  AliAnalysisDataContainer *coutputEvtCuts = mgr->CreateContainer(
      EvtCutsName.Data(), TList::Class(), AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), EvtCutsName.Data()));

  TString TrackCutsName = Form("%sTrackCuts_%s_%s", addon.Data(), suffixTrigger.Data(), suffix.Data());
  AliAnalysisDataContainer *couputTrkCuts = mgr->CreateContainer(
      TrackCutsName.Data(), TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), TrackCutsName.Data()));

  TString AntiTrackCutsName = Form("%sAntiTrackCuts_%s_%s", addon.Data(), suffixTrigger.Data(), suffix.Data());
  AliAnalysisDataContainer *coutputAntiTrkCuts = mgr->CreateContainer(
      AntiTrackCutsName.Data(), TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), AntiTrackCutsName.Data()));

  AliAnalysisDataContainer *coutputv0Cuts;
  TString v0CutsName = Form("%sv0Cuts_%s_%s", addon.Data(), suffixTrigger.Data(), suffix.Data());
  coutputv0Cuts = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      v0CutsName.Data(),
      TList::Class(), AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), v0CutsName.Data()));

  AliAnalysisDataContainer *coutputAntiv0Cuts;
  TString Antiv0CutsName = Form("%sAntiv0Cuts_%s_%s", addon.Data(), suffixTrigger.Data(), suffix.Data());
  coutputAntiv0Cuts = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      Antiv0CutsName.Data(),
      TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), Antiv0CutsName.Data()));

  AliAnalysisDataContainer *coutputResults;
  TString ResultsName = Form("%sResults_%s_%s", addon.Data(), suffixTrigger.Data(), suffix.Data());
  coutputResults = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      ResultsName.Data(),
      TList::Class(), AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), ResultsName.Data()));

  AliAnalysisDataContainer *coutputResultsQA;
  TString ResultsQAName = Form("%sResultsQA_%s_%s", addon.Data(), suffixTrigger.Data(), suffix.Data());
  coutputResultsQA = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      ResultsQAName.Data(),
      TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), ResultsQAName.Data()));

  AliAnalysisDataContainer *coutputResultsSample;
  TString ResultsSampleName = Form("%sResultsSample_%s_%s", addon.Data(), suffixTrigger.Data(), suffix.Data());
  coutputResultsSample = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      ResultsSampleName.Data(),
      TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), ResultsSampleName.Data()));

  AliAnalysisDataContainer *coutputResultsSampleQA;
  TString ResultsSampleQAName = Form("%sResultsSampleQA_%s_%s", addon.Data(), suffixTrigger.Data(), suffix.Data());
  coutputResultsSampleQA = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      ResultsSampleQAName.Data(),
      TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), ResultsSampleQAName.Data()));

  AliAnalysisDataContainer *coutputTrkCutsMC;
  AliAnalysisDataContainer *coutputAntiTrkCutsMC;
  AliAnalysisDataContainer *coutputv0CutsMC;
  AliAnalysisDataContainer *coutputAntiv0CutsMC;
  if (isMC) {
    TString TrkCutsMCName = Form("%sTrkCutsMC_%s_%s", addon.Data(), suffixTrigger.Data(), suffix.Data());
    coutputTrkCutsMC = mgr->CreateContainer(
        //@suppress("Invalid arguments") it works ffs
        TrkCutsMCName.Data(),
        TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), TrkCutsMCName.Data()));

    TString AntiTrkCutsMCName = Form("%sAntiTrkCutsMC_%s_%s", addon.Data(), suffixTrigger.Data(), suffix.Data());
    coutputAntiTrkCutsMC = mgr->CreateContainer(
        //@suppress("Invalid arguments") it works ffs
        AntiTrkCutsMCName.Data(),
        TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), AntiTrkCutsMCName.Data()));

    TString v0CutsMCName = Form("%sv0CutsMC_%s_%s", addon.Data(), suffixTrigger.Data(), suffix.Data());
    coutputv0CutsMC = mgr->CreateContainer(
        //@suppress("Invalid arguments") it works ffs
        v0CutsMCName.Data(),
        TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), v0CutsMCName.Data()));

    TString Antiv0CutsMCName = Form("%sAntiv0CutsMC_%s_%s", addon.Data(), suffixTrigger.Data(), suffix.Data());
    coutputAntiv0CutsMC = mgr->CreateContainer(
        //@suppress("Invalid arguments") it works ffs
        Antiv0CutsMCName.Data(),
        TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), Antiv0CutsMCName.Data()));

  }

  TString ThreeBodyName = Form("%sThreeBody_%s_%s", addon.Data(), suffixTrigger.Data(), suffix.Data());
  AliAnalysisDataContainer *coutputThreeBody = mgr->CreateContainer(
    //@suppress("Invalid arguments") it works ffs
    ThreeBodyName.Data(),
    TList::Class(),
    AliAnalysisManager::kOutputContainer,
    Form("%s:%s", file.Data(), ThreeBodyName.Data()));

  AliAnalysisDataContainer *coutputEvtCutsTrigger;
  AliAnalysisDataContainer *couputTrkCutsTrigger;
  AliAnalysisDataContainer *coutputAntiTrkCutsTrigger;
  AliAnalysisDataContainer *coutputv0CutsTrigger;
  AliAnalysisDataContainer *coutputAntiv0CutsTrigger;
 

  AliAnalysisTaskThreeBodyFemtoMixedCharge* taskNano;

  taskNano= new AliAnalysisTaskThreeBodyFemtoMixedCharge("femtoNanoThreeBody", isMC);
  if (!fullBlastQA)
  { 
    taskNano->SetRunTaskLightWeight(true);
  }

  if (trigger == 0) { 
      taskNano->SelectCollisionCandidates(AliVEvent::kHighMultV0);  
    } else if (trigger == 1){     
      taskNano->SelectCollisionCandidates(AliVEvent::kINT7);  
    } 
  taskNano->SetEventCuts(evtCuts);  
  taskNano->SetProtonCuts(TrackCuts); 
  taskNano->SetAntiProtonCuts(AntiTrackCuts); 
  taskNano->Setv0Cuts(v0Cuts);  
  taskNano->SetAntiv0Cuts(Antiv0Cuts);  
  taskNano->SetCorrelationConfig(config); 
  taskNano->SetRunThreeBodyHistograms(true);
  taskNano->SetClosePairRejectionForAll(ClosePairRejectionForAll);
  taskNano->SetRun2Body(run2Body);
  taskNano->SetMixingChoice(mixinfChoice);
  taskNano->SetSame2Mixed1Choice(mix21);
  taskNano->SetWhichTripletsToRun(whichTripletsToRun);
  
  

  mgr->AddTask(taskNano); 
  
  mgr->ConnectInput(taskNano, 0, cinput); 
  mgr->ConnectOutput(taskNano, 1, coutputEvtCuts);  
  mgr->ConnectOutput(taskNano, 2, couputTrkCuts); 
  mgr->ConnectOutput(taskNano, 3, coutputAntiTrkCuts);  
  mgr->ConnectOutput(taskNano, 4, coutputv0Cuts); 
  mgr->ConnectOutput(taskNano, 5, coutputAntiv0Cuts); 
  mgr->ConnectOutput(taskNano, 6, coutputResults);  
  mgr->ConnectOutput(taskNano, 7, coutputResultsQA);  
  mgr->ConnectOutput(taskNano, 8, coutputResultsSample);  
  mgr->ConnectOutput(taskNano, 9, coutputResultsSampleQA);  
  mgr->ConnectOutput(taskNano, 10, coutputThreeBody);  
  if (isMC) { 
    mgr->ConnectOutput(taskNano, 11, coutputTrkCutsMC); 
    mgr->ConnectOutput(taskNano, 12, coutputAntiTrkCutsMC); 
    mgr->ConnectOutput(taskNano, 13, coutputv0CutsMC);  
    mgr->ConnectOutput(taskNano, 14, coutputAntiv0CutsMC);  
  } 

    
  return taskNano;
  

  
}
