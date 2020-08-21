#if !defined(__CINT__) || defined(__CLING__)
#include <vector>
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskThreeBodyFemto.h"
#include "AliAnalysisTaskThreeBodyFemtoAOD.h"
#include "AliFemtoDreamEventCuts.h"
#include "AliFemtoDreamTrackCuts.h"
#include "AliFemtoDreamCascadeCuts.h"
#include "AliFemtoDreamCollConfig.h"
#endif

AliAnalysisTaskSE *AddTaskThreeBodyFemto(int trigger = 0, bool fullBlastQA = true,
                                     bool isMC = false, bool isNano = true, bool triggerOn = false,
                                     bool triggerCutVariation = false, const char *cutVariation = "0") {



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
  if(triggerCutVariation){
    TrackCuts->SetPtRange(0.3, 0.7);
    TrackCuts->SetEtaRange(-0.7, 0.7);
    TrackCuts->SetNClsTPC(65);
    TrackCuts->SetPID(AliPID::kProton, 0.75,5.); 
    TrackCuts->SetRejLowPtPionsTOF(true);
    TrackCuts->SetCutSmallestSig(true);

    AntiTrackCuts->SetPtRange(0.3, 0.7);
    AntiTrackCuts->SetEtaRange(-0.7, 0.7);
    AntiTrackCuts->SetNClsTPC(65);
    AntiTrackCuts->SetPID(AliPID::kProton, 0.75,5.); 
    AntiTrackCuts->SetRejLowPtPionsTOF(true);
    AntiTrackCuts->SetCutSmallestSig(true);

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
  if(triggerCutVariation){ 
    v0Cuts->SetPtRange(0.2, 999.);
    v0Cuts->SetCutCPA(0.98);
    v0Cuts->SetCutDCADaugToPrimVtx(0.07);
    v0Cuts->SetCutDCADaugTov0Vtx(1.3);
    Posv0Daug->SetNClsTPC(65);
    Negv0Daug->SetNClsTPC(65);
    Posv0Daug->SetEtaRange(-0.7, 0.7);
    Negv0Daug->SetEtaRange(-0.7, 0.7);
    Posv0Daug->SetPID(AliPID::kProton, 999., 6.5);
    Negv0Daug->SetPID(AliPID::kPion, 999., 6.5);
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
  if(triggerCutVariation){ 
    Antiv0Cuts->SetPtRange(0.2, 999.);
    Antiv0Cuts->SetCutCPA(0.98);
    Antiv0Cuts->SetCutDCADaugToPrimVtx(0.07);
    Antiv0Cuts->SetCutDCADaugTov0Vtx(1.3);
    PosAntiv0Daug->SetNClsTPC(65);
    NegAntiv0Daug->SetNClsTPC(65);
    PosAntiv0Daug->SetEtaRange(-0.7, 0.7);
    NegAntiv0Daug->SetEtaRange(-0.7, 0.7);
    PosAntiv0Daug->SetPID(AliPID::kPion, 999., 6.5);
    NegAntiv0Daug->SetPID(AliPID::kProton, 999., 6.5);
  }
  Antiv0Cuts->SetPosDaugterTrackCuts(PosAntiv0Daug);
  Antiv0Cuts->SetNegDaugterTrackCuts(NegAntiv0Daug);
  Antiv0Cuts->SetPDGCodePosDaug(211);  //Pion
  Antiv0Cuts->SetPDGCodeNegDaug(2212);  //Proton
  Antiv0Cuts->SetPDGCodev0(-3122);  //Lambda

  // If  cut variations needed
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

  TString EvtCutsName = Form("%sEvtCuts%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *coutputEvtCuts = mgr->CreateContainer(
      EvtCutsName.Data(), TList::Class(), AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), EvtCutsName.Data()));

  TString TrackCutsName = Form("%sTrackCuts%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *couputTrkCuts = mgr->CreateContainer(
      TrackCutsName.Data(), TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), TrackCutsName.Data()));

  TString AntiTrackCutsName = Form("%sAntiTrackCuts%s", addon.Data(),
                                   suffix.Data());
  AliAnalysisDataContainer *coutputAntiTrkCuts = mgr->CreateContainer(
      AntiTrackCutsName.Data(), TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), AntiTrackCutsName.Data()));

  AliAnalysisDataContainer *coutputv0Cuts;
  TString v0CutsName = Form("%sv0Cuts%s", addon.Data(), suffix.Data());
  coutputv0Cuts = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      v0CutsName.Data(),
      TList::Class(), AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), v0CutsName.Data()));

  AliAnalysisDataContainer *coutputAntiv0Cuts;
  TString Antiv0CutsName = Form("%sAntiv0Cuts%s", addon.Data(), suffix.Data());
  coutputAntiv0Cuts = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      Antiv0CutsName.Data(),
      TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), Antiv0CutsName.Data()));

  AliAnalysisDataContainer *coutputResults;
  TString ResultsName = Form("%sResults%s", addon.Data(), suffix.Data());
  coutputResults = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      ResultsName.Data(),
      TList::Class(), AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), ResultsName.Data()));

  AliAnalysisDataContainer *coutputResultsQA;
  TString ResultsQAName = Form("%sResultsQA%s", addon.Data(), suffix.Data());
  coutputResultsQA = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      ResultsQAName.Data(),
      TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), ResultsQAName.Data()));

  AliAnalysisDataContainer *coutputResultsSample;
  TString ResultsSampleName = Form("%sResultsSample%s", addon.Data(),
                                   suffix.Data());
  coutputResultsSample = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      ResultsSampleName.Data(),
      TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), ResultsSampleName.Data()));

  AliAnalysisDataContainer *coutputResultsSampleQA;
  TString ResultsSampleQAName = Form("%sResultsSampleQA%s", addon.Data(),
                                     suffix.Data());
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
    TString TrkCutsMCName = Form("%sTrkCutsMC%s", addon.Data(), suffix.Data());
    coutputTrkCutsMC = mgr->CreateContainer(
        //@suppress("Invalid arguments") it works ffs
        TrkCutsMCName.Data(),
        TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), TrkCutsMCName.Data()));

    TString AntiTrkCutsMCName = Form("%sAntiTrkCutsMC%s", addon.Data(),
                                     suffix.Data());
    coutputAntiTrkCutsMC = mgr->CreateContainer(
        //@suppress("Invalid arguments") it works ffs
        AntiTrkCutsMCName.Data(),
        TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), AntiTrkCutsMCName.Data()));

    TString v0CutsMCName = Form("%sv0CutsMC%s", addon.Data(), suffix.Data());
    coutputv0CutsMC = mgr->CreateContainer(
        //@suppress("Invalid arguments") it works ffs
        v0CutsMCName.Data(),
        TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), v0CutsMCName.Data()));

    TString Antiv0CutsMCName = Form("%sAntiv0CutsMC%s", addon.Data(),
                                    suffix.Data());
    coutputAntiv0CutsMC = mgr->CreateContainer(
        //@suppress("Invalid arguments") it works ffs
        Antiv0CutsMCName.Data(),
        TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), Antiv0CutsMCName.Data()));

  }

  AliAnalysisDataContainer *coutputThreeBody;
  TString ThreeBodyName = Form("%sThreeBody%s", addon.Data(), suffix.Data());
  coutputThreeBody = mgr->CreateContainer(
    //@suppress("Invalid arguments") it works ffs
    ThreeBodyName.Data(),
    TList::Class(),
    AliAnalysisManager::kOutputContainer,
    Form("%s:%s", file.Data(), ThreeBodyName.Data()));

  AliAnalysisTaskThreeBodyFemto* taskNano;
  AliAnalysisTaskThreeBodyFemtoAOD* taskAOD;
  if(isNano){
    taskNano= new AliAnalysisTaskThreeBodyFemto("femtoNanoThreeBody", isMC);
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
  }
  else{
    taskAOD= new AliAnalysisTaskThreeBodyFemtoAOD("femtoAODThreeBody", isMC, triggerOn);
    if (!fullBlastQA)
    { 
      taskAOD->SetRunTaskLightWeight(true);
    }

    if (trigger == 0) { 
        taskAOD->SelectCollisionCandidates(AliVEvent::kHighMultV0);  
      } else if (trigger == 1){     
        taskAOD->SelectCollisionCandidates(AliVEvent::kINT7);  
      } 
    taskAOD->SetEventCuts(evtCuts);  
    taskAOD->SetProtonCuts(TrackCuts); 
    taskAOD->SetAntiProtonCuts(AntiTrackCuts); 
    taskAOD->Setv0Cuts(v0Cuts);  
    taskAOD->SetAntiv0Cuts(Antiv0Cuts);  
    taskAOD->SetCorrelationConfig(config); 
    taskAOD->SetRunThreeBodyHistograms(true);
    mgr->AddTask(taskAOD); 
    
    mgr->ConnectInput(taskAOD, 0, cinput); 
    mgr->ConnectOutput(taskAOD, 1, coutputEvtCuts);  
    mgr->ConnectOutput(taskAOD, 2, couputTrkCuts); 
    mgr->ConnectOutput(taskAOD, 3, coutputAntiTrkCuts);  
    mgr->ConnectOutput(taskAOD, 4, coutputv0Cuts); 
    mgr->ConnectOutput(taskAOD, 5, coutputAntiv0Cuts); 
    mgr->ConnectOutput(taskAOD, 6, coutputResults);  
    mgr->ConnectOutput(taskAOD, 7, coutputResultsQA);  
    mgr->ConnectOutput(taskAOD, 8, coutputResultsSample);  
    mgr->ConnectOutput(taskAOD, 9, coutputResultsSampleQA);  
    mgr->ConnectOutput(taskAOD, 10, coutputThreeBody);  
    if (isMC) { 
      mgr->ConnectOutput(taskAOD, 11, coutputTrkCutsMC); 
      mgr->ConnectOutput(taskAOD, 12, coutputAntiTrkCutsMC); 
      mgr->ConnectOutput(taskAOD, 13, coutputv0CutsMC);  
      mgr->ConnectOutput(taskAOD, 14, coutputAntiv0CutsMC);  
    } 
  }

      
  if (isNano) {
    return taskNano;
  } else {
    return taskAOD;
  }

  
}
