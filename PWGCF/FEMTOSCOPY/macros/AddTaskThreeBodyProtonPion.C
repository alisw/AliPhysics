#if !defined(__CINT__) || defined(__CLING__)
#include <vector>
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
//#include "AliAnalysisTaskProtonProtonKaons.h"
#include "AliAnalysisTaskThreeBodyProtonPrimary.h"
//#include "AliAnalysisTaskThreeBodyFemtoAOD.h"
#include "AliFemtoDreamEventCuts.h"
#include "AliFemtoDreamTrackCuts.h"
#include "AliFemtoDreamCascadeCuts.h"
#include "AliFemtoDreamCollConfig.h"
#endif

AliAnalysisTaskSE *AddTaskThreeBodyProtonPion(int trigger = 0, bool fullBlastQA = true, bool isMC = false, bool isNano = true, float Q3Limit = 0.6, float Q3LimitSample = 3.0,float Q3LimitSample2 = 3.0, float Q3LimitFraction = 0.5, float Q3LimitSampleFraction = 0.01, float Q3LimitSampleFraction2 = 0.01, const char *cutVariation = "0", bool ClosePairRejectionForAll = false, bool UseSphericityCut = false){


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

  if (UseSphericityCut){
    float SpherDown = 0.7;
    evtCuts->SetSphericityCuts(SpherDown, 1.0, 0.5); // THINK IF NEEDED FOR THREE BODY
  }

// Track Cuts for Protons and Antiprotons ================================================================================================
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

  // Track Cuts for Pions  =================================================================================================================

  AliFemtoDreamTrackCuts *TrackCutsPion=new AliFemtoDreamTrackCuts();
  TrackCutsPion->SetIsMonteCarlo(isMC);
  TrackCutsPion->SetCutCharge(1);
  TrackCutsPion->SetPtRange(0.14, 4.0);
  TrackCutsPion->SetEtaRange(-0.8, 0.8);
  TrackCutsPion->SetNClsTPC(80);
  // Not mention in AN oder Indico
  TrackCutsPion->SetDCAReCalculation(true);//Get the dca from the PropagateToVetex
  TrackCutsPion->SetFilterBit(128);//96); // Filterbit 5+6
  TrackCutsPion->SetDCAVtxZ(0.3);
  TrackCutsPion->SetDCAVtxXY(0.3);
  // Cut on avrg. separation in TPC: <Dr> < 12 cm (10 cm, 3 cm); Share quality < 1.0; share fraction < 0.05

  // FOR NOW OF BEACUSE OF MAX's INFORMATION
  //TrackCutsPion->SetCutSharedCls(true);

  TrackCutsPion->SetNClsTPC(80); // In Indico + additional Chi²/NDF <4
  TrackCutsPion->SetPID(AliPID::kPion, 0.5);
  TrackCutsPion->SetRejLowPtPionsTOF(false);
  TrackCutsPion->SetMinimalBooking(false);
  //this checks if the sigma of the wanted hypothesis is the smallest, and if
  //another particle has a smaller sigma, the track is rejected.
  // Not mention in AN oder Indico
  //TrackCutsPion->SetCutSmallestSig(true);
  TrackCutsPion->SetPlotDCADist(true);

 //The same things for negative pions
  AliFemtoDreamTrackCuts *TrackCutsAntiPion=new AliFemtoDreamTrackCuts();
  TrackCutsAntiPion->SetIsMonteCarlo(isMC);
  TrackCutsAntiPion->SetCutCharge(-1);
  TrackCutsAntiPion->SetPtRange(0.14, 4.0);
  TrackCutsAntiPion->SetEtaRange(-0.8, 0.8);
  TrackCutsAntiPion->SetNClsTPC(80);
  TrackCutsAntiPion->SetDCAReCalculation(true);

  // FOR NOW OF BEACUSE OF MAX's INFORMATION
  //TrackCutsAntiPion->SetCutSharedCls(true);}

  TrackCutsAntiPion->SetNClsTPC(80);
  TrackCutsAntiPion->SetPID(AliPID::kPion, 0.5);
  TrackCutsAntiPion->SetRejLowPtPionsTOF(false);
  TrackCutsAntiPion->SetMinimalBooking(false);
  //TrackCutsAntiPion->SetCutSmallestSig(true);
  TrackCutsAntiPion->SetPlotDCADist(true);

  TrackCutsAntiPion->SetFilterBit(128);//96);
  TrackCutsAntiPion->SetDCAVtxZ(0.3);
  TrackCutsAntiPion->SetDCAVtxXY(0.3);

  //TrackCutsAntiPion->SetCutSharedCls(true);


  if (!fullBlastQA) {
    evtCuts->SetMinimalBooking(true);
    TrackCuts->SetMinimalBooking(true);
    AntiTrackCuts->SetMinimalBooking(true);
    TrackCutsPion->SetMinimalBooking(true);
    TrackCutsAntiPion->SetMinimalBooking(true);
    //GANESHA add here stuff for v0
  }



  AliFemtoDreamCollConfig *config = new AliFemtoDreamCollConfig("Femto",
                                                                "Femto", false);
  // Femto Collection
  std::vector<int> PDGParticles;
  PDGParticles.push_back(2212); //Proton
  PDGParticles.push_back(2212);
  PDGParticles.push_back(211); //Pion
  PDGParticles.push_back(211);

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

  TString addon = "PPion";
  TString file = AliAnalysisManager::GetCommonFileName();

  TString EvtCutsName = Form("%sEvtCuts_%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *coutputEvtCuts = mgr->CreateContainer(
      EvtCutsName.Data(), TList::Class(), AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), EvtCutsName.Data()));

  TString TrackCutsName = Form("%sTrackCuts_%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *couputTrkCuts = mgr->CreateContainer(
      TrackCutsName.Data(), TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), TrackCutsName.Data()));

  TString AntiTrackCutsName = Form("%sAntiTrackCuts_%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *coutputAntiTrkCuts = mgr->CreateContainer(
      AntiTrackCutsName.Data(), TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), AntiTrackCutsName.Data()));

  AliAnalysisDataContainer *coutputPionCuts;
  TString PionCutsName = Form("%sPionCuts_%s", addon.Data(), suffix.Data());
  coutputPionCuts = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      PionCutsName.Data(),
      TList::Class(), AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), PionCutsName.Data()));

  AliAnalysisDataContainer *coutputAntiPionCuts;
  TString AntiPionCutsName = Form("%sAntiPionCuts_%s", addon.Data(), suffix.Data());
  coutputAntiPionCuts = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      AntiPionCutsName.Data(),
      TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), AntiPionCutsName.Data()));

  AliAnalysisDataContainer *coutputResults;
  TString ResultsName = Form("%sResults_%s", addon.Data(), suffix.Data());
  coutputResults = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      ResultsName.Data(),
      TList::Class(), AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), ResultsName.Data()));

  AliAnalysisDataContainer *coutputResultsQA;
  TString ResultsQAName = Form("%sResultsQA_%s", addon.Data(), suffix.Data());
  coutputResultsQA = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      ResultsQAName.Data(),
      TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), ResultsQAName.Data()));

  AliAnalysisDataContainer *coutputResultsSample;
  TString ResultsSampleName = Form("%sResultsSample_%s", addon.Data(), suffix.Data());
  coutputResultsSample = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      ResultsSampleName.Data(),
      TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), ResultsSampleName.Data()));

  AliAnalysisDataContainer *coutputResultsSampleQA;
  TString ResultsSampleQAName = Form("%sResultsSampleQA_%s", addon.Data(), suffix.Data());
  coutputResultsSampleQA = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      ResultsSampleQAName.Data(),
      TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), ResultsSampleQAName.Data()));

  AliAnalysisDataContainer *coutputTrkCutsMC;
  AliAnalysisDataContainer *coutputAntiTrkCutsMC;
  AliAnalysisDataContainer *coutputPionCutsMC;
  AliAnalysisDataContainer *coutputAntiPionCutsMC;
  if (isMC) {
    TString TrkCutsMCName = Form("%sTrkCutsMC_%s", addon.Data(), suffix.Data());
    coutputTrkCutsMC = mgr->CreateContainer(
        //@suppress("Invalid arguments") it works ffs
        TrkCutsMCName.Data(),
        TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), TrkCutsMCName.Data()));

    TString AntiTrkCutsMCName = Form("%sAntiTrkCutsMC_%s", addon.Data(), suffix.Data());
    coutputAntiTrkCutsMC = mgr->CreateContainer(
        //@suppress("Invalid arguments") it works ffs
        AntiTrkCutsMCName.Data(),
        TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), AntiTrkCutsMCName.Data()));

    TString PionCutsMCName = Form("%sPionCutsMC_%s", addon.Data(), suffix.Data());
    coutputPionCutsMC = mgr->CreateContainer(
        //@suppress("Invalid arguments") it works ffs
        PionCutsMCName.Data(),
        TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), PionCutsMCName.Data()));

    TString AntiPionCutsMCName = Form("%sAntiPionCutsMC_%s", addon.Data(), suffix.Data());
    coutputAntiPionCutsMC = mgr->CreateContainer(
        //@suppress("Invalid arguments") it works ffs
        AntiPionCutsMCName.Data(),
        TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), AntiPionCutsMCName.Data()));

  }

  TString ThreeBodyName = Form("%sThreeBodyProtonPion_%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *coutputThreeBody = mgr->CreateContainer(
    //@suppress("Invalid arguments") it works ffs
    ThreeBodyName.Data(),
    TList::Class(),
    AliAnalysisManager::kOutputContainer,
    Form("%s:%s", file.Data(), ThreeBodyName.Data()));



  AliAnalysisTaskThreeBodyProtonPrimary* taskNano;
//  AliAnalysisTaskThreeBodyProtonPrimaryAOD* taskAOD;
  if(isNano){
    taskNano= new AliAnalysisTaskThreeBodyProtonPrimary("femtoNanoThreeBodyProtonPion", isMC);
    if (!fullBlastQA)
    {
      taskNano->SetRunTaskLightWeight(true);
    }

    taskNano->SetEventCuts(evtCuts);
    taskNano->SetProtonCuts(TrackCuts);
    taskNano->SetAntiProtonCuts(AntiTrackCuts);
    taskNano->SetPrimaryCuts(TrackCutsPion);
    taskNano->SetAntiPrimaryCuts(TrackCutsAntiPion);
    taskNano->SetCorrelationConfig(config);
    taskNano->SetRunThreeBodyHistograms(true);
    taskNano->SetClosePairRejectionForAll(ClosePairRejectionForAll);
    taskNano->SetCleanWithLambdas(false);
    mgr->AddTask(taskNano);

    mgr->ConnectInput(taskNano, 0, cinput);
    mgr->ConnectOutput(taskNano, 1, coutputEvtCuts);
    mgr->ConnectOutput(taskNano, 2, couputTrkCuts);
    mgr->ConnectOutput(taskNano, 3, coutputAntiTrkCuts);
    mgr->ConnectOutput(taskNano, 4, coutputPionCuts);
    mgr->ConnectOutput(taskNano, 5, coutputAntiPionCuts);
    mgr->ConnectOutput(taskNano, 8, coutputResults);
    mgr->ConnectOutput(taskNano, 9, coutputResultsQA);
    mgr->ConnectOutput(taskNano, 10, coutputResultsSample);
    mgr->ConnectOutput(taskNano, 11, coutputResultsSampleQA);
    mgr->ConnectOutput(taskNano, 12, coutputThreeBody);
    if (isMC) {
      mgr->ConnectOutput(taskNano, 13, coutputTrkCutsMC);
      mgr->ConnectOutput(taskNano, 14, coutputAntiTrkCutsMC);
      mgr->ConnectOutput(taskNano, 15, coutputPionCutsMC);
      mgr->ConnectOutput(taskNano, 16, coutputAntiPionCutsMC);
    }
  }
  else{
/*
    taskAOD= new AliAnalysisTaskThreeBodyProtonPrimaryAOD("femtoAODProtonProtonPions", isMC, triggerOn);
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
    taskAOD->Setv0Cuts(TrackCutsPion);
    taskAOD->SetAntiv0Cuts(TrackCutsAntiPion);
    taskAOD->SetCorrelationConfig(config);
    taskAOD->SetRunThreeBodyHistograms(true);
    //taskAOD->SetTriggerOn(triggerOn);
    taskAOD->SetIsMC(isMC);




    taskAOD->SetQ3Limit(Q3Limit);
    taskAOD->SetQ3LimitSample(Q3LimitSample) ;
    taskAOD->SetQ3LimitSample2(Q3LimitSample2) ;
    taskAOD->SetQ3LimitSampleFraction( Q3LimitSampleFraction) ;
    taskAOD->SetQ3LimitSampleFraction2( Q3LimitSampleFraction2) ;
    taskAOD->SetQ3LimitFraction( Q3LimitFraction) ;


    mgr->AddTask(taskAOD);

    mgr->ConnectInput(taskAOD, 0, cinput);
    mgr->ConnectOutput(taskAOD, 1, coutputEvtCuts);
    mgr->ConnectOutput(taskAOD, 2, couputTrkCuts);
    mgr->ConnectOutput(taskAOD, 3, coutputAntiTrkCuts);
    mgr->ConnectOutput(taskAOD, 4, coutputPionCuts);
    mgr->ConnectOutput(taskAOD, 5, coutputAntiPionCuts);
    mgr->ConnectOutput(taskAOD, 6, coutputResults);
    mgr->ConnectOutput(taskAOD, 7, coutputResultsQA);
    mgr->ConnectOutput(taskAOD, 8, coutputResultsSample);
    mgr->ConnectOutput(taskAOD, 9, coutputResultsSampleQA);
    mgr->ConnectOutput(taskAOD, 10, coutputThreeBody);
    if (isMC) {
      mgr->ConnectOutput(taskAOD, 16, coutputTrkCutsMC);
      mgr->ConnectOutput(taskAOD, 17, coutputAntiTrkCutsMC);
      mgr->ConnectOutput(taskAOD, 18, coutputPionCutsMC);
      mgr->ConnectOutput(taskAOD, 19, coutputAntiPionCutsMC);
    }
*/
  }


  if (isNano) {
    return taskNano;
  } else {
  //  return taskAOD;
  }


}
