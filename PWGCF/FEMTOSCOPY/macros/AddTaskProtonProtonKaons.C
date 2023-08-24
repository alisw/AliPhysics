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

AliAnalysisTaskSE *AddTaskProtonProtonKaons(int trigger = 0, bool fullBlastQA = true,
                                     bool isMC = false, bool isNano = true, bool triggerOn = false, int MixingDepth = 30,
                                     float Q3Limit = 0.6, float Q3LimitSample = 3.0,float Q3LimitSample2 = 3.0, float Q3LimitFraction = 0.5, float Q3LimitSampleFraction = 0.01, float Q3LimitSampleFraction2 = 0.01,
                                     const char *cutVariation = "0", bool turnoffClosePairRejectionCompletely = false, bool ClosePairRejectionForAll = "false",
                                     const char *triggerVariation = "0", bool RunPlotPt = true, bool RunPlotQ3Vsq = false, bool UseSphericityCut = false, bool DoOnlyThreeBody = false, bool RunOfficialTwoBody=false, int KaonCut = 1, bool DoTwoPrimary = false, bool StandardMixing = false, bool DoKinematicPlots = false, bool RunPlotMult = true, bool RunCorrDeltaPhi = true, bool RunPlotPhiTheta = false, bool RunPlotP1 = false, bool CutElectrons = false) {


  TString suffix = TString::Format("%s", cutVariation);
  TString suffixTrigger = TString::Format("%s", triggerVariation);

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
  // Track Cuts
  AliFemtoDreamTrackCuts *TrackCuts = AliFemtoDreamTrackCuts::PrimProtonCuts(
      isMC, true, false, false);
  if(trigger==0){
    TrackCuts->SetFilterBit(128);
  }else{
    TrackCuts->SetFilterBit(256);
  }
  TrackCuts->SetCutCharge(1);

  AliFemtoDreamTrackCuts *AntiTrackCuts =
      AliFemtoDreamTrackCuts::PrimProtonCuts(isMC, true, false, false);
  if(trigger==0){
    AntiTrackCuts->SetFilterBit(128);
  }else{
    AntiTrackCuts->SetFilterBit(256);
  }
  AntiTrackCuts->SetCutCharge(-1);
  // If  cut variation needed
  if(suffix=="1" || suffix=="8"){
    TrackCuts->SetEtaRange(-0.9, 0.9);
    AntiTrackCuts->SetEtaRange(-0.9, 0.9);
  }

  if(suffix=="69"){
    TrackCuts->SetRapidityRange(-0.5, 0.5, 0.9382720881);
    AntiTrackCuts->SetRapidityRange(-0.5, 0.5, 0.9382720881);
  }

  //Kaon Cuts
  AliFemtoDreamTrackCuts *KaonCuts = AliFemtoDreamTrackCuts::PrimKaonCuts(
    isMC, true, false, false);
  KaonCuts->SetPtRange(0.15, 10);
  if(trigger==0){
    KaonCuts->SetFilterBit(128);
  }else{
    KaonCuts->SetFilterBit(256);
  }
  KaonCuts->SetCutCharge(1);
  if(KaonCut==0){ // cuts by Oton
   KaonCuts->SetPIDkd(true,false,3,3,3);
  }else if(KaonCut==1){ // cuts by Ramona
   KaonCuts->SetPIDkd(true,true);
  }

  //AntiKaon Cuts
  AliFemtoDreamTrackCuts *AntiKaonCuts = AliFemtoDreamTrackCuts::PrimKaonCuts(
    isMC, true, false, false);
  AntiKaonCuts->SetPtRange(0.15, 10);
  if(trigger==0){  
    AntiKaonCuts->SetFilterBit(128);
  }else{
    AntiKaonCuts->SetFilterBit(256);
  }
  AntiKaonCuts->SetCutCharge(-1);
  if(KaonCut==0){ // cuts by Oton
   AntiKaonCuts->SetPIDkd(true,false,3,3,3);
  }else if(KaonCut==1){ // cuts by Ramona
   AntiKaonCuts->SetPIDkd(true,true);
  }
  if(suffix=="69"){
    KaonCuts->SetRapidityRange(-0.5, 0.5, 0.493677);
    AntiKaonCuts->SetRapidityRange(-0.5, 0.5, 0.493677);
  }

  if (!fullBlastQA) {
    evtCuts->SetMinimalBooking(true);
    TrackCuts->SetMinimalBooking(true);
    AntiTrackCuts->SetMinimalBooking(true);
    KaonCuts->SetMinimalBooking(true);
    AntiKaonCuts->SetMinimalBooking(true);
  }

  double DeltaPhiMaxpp = 0.017;      
  double DeltaEtaMaxpp = 0.017;
  double DeltaPhiMaxpKplus = 0.04;
  double DeltaEtaMaxpKplus = 0.012;

  AliFemtoDreamCollConfig *config = new AliFemtoDreamCollConfig("Femto",
                                                                "Femto", false);
  // Femto Collection
  std::vector<int> PDGParticles;
  PDGParticles.push_back(2212);
  PDGParticles.push_back(2212);
  PDGParticles.push_back(321);
  PDGParticles.push_back(321);

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
  config->SetDeltaEtaMax(0.);
  config->SetDeltaPhiMax(0.);
//  config->SetDeltaEtaMax(0.017);
//  config->SetDeltaPhiMax(0.017);
  config->SetExtendedQAPairs(pairQA);
  config->SetMixingDepth(MixingDepth);
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

  TString addon = "PK";
  TString file = AliAnalysisManager::GetCommonFileName();

  TString EvtCutsName = Form("%sEvtCuts_%s_%s", addon.Data(), suffixTrigger.Data(), suffix.Data());
  AliAnalysisDataContainer *coutputEvtCuts = mgr->CreateContainer(
      EvtCutsName.Data(), TList::Class(), AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), EvtCutsName.Data()));

  TString TrackCutsName = Form("%sProtonCuts_%s_%s", addon.Data(), suffixTrigger.Data(), suffix.Data());
  AliAnalysisDataContainer *couputTrkCuts = mgr->CreateContainer(
      TrackCutsName.Data(), TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), TrackCutsName.Data()));

  TString AntiTrackCutsName = Form("%sAntiProtonCuts_%s_%s", addon.Data(), suffixTrigger.Data(), suffix.Data());
  AliAnalysisDataContainer *coutputAntiTrkCuts = mgr->CreateContainer(
      AntiTrackCutsName.Data(), TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), AntiTrackCutsName.Data()));

  AliAnalysisDataContainer *coutputKaonCuts;
  TString KaonCutsName = Form("%sKaonCuts_%s_%s", addon.Data(), suffixTrigger.Data(), suffix.Data());
  coutputKaonCuts = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      KaonCutsName.Data(),
      TList::Class(), AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), KaonCutsName.Data()));

  AliAnalysisDataContainer *coutputAntiKaonCuts;
  TString AntiKaonCutsName = Form("%sAntiKaonCuts_%s_%s", addon.Data(), suffixTrigger.Data(), suffix.Data());
  coutputAntiKaonCuts = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      AntiKaonCutsName.Data(),
      TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), AntiKaonCutsName.Data()));

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
  AliAnalysisDataContainer *coutputKaonCutsMC;
  AliAnalysisDataContainer *coutputAntiKaonCutsMC;
  if (isMC) {
    TString TrkCutsMCName = Form("%sProtonCutsMC_%s_%s", addon.Data(), suffixTrigger.Data(), suffix.Data());
    coutputTrkCutsMC = mgr->CreateContainer(
        //@suppress("Invalid arguments") it works ffs
        TrkCutsMCName.Data(),
        TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), TrkCutsMCName.Data()));

    TString AntiTrkCutsMCName = Form("%sAntiProtonCutsMC_%s_%s", addon.Data(), suffixTrigger.Data(), suffix.Data());
    coutputAntiTrkCutsMC = mgr->CreateContainer(
        //@suppress("Invalid arguments") it works ffs
        AntiTrkCutsMCName.Data(),
        TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), AntiTrkCutsMCName.Data()));

    TString KaonCutsMCName = Form("%sKaonCutsMC_%s_%s", addon.Data(), suffixTrigger.Data(), suffix.Data());
    coutputKaonCutsMC = mgr->CreateContainer(
        //@suppress("Invalid arguments") it works ffs
        KaonCutsMCName.Data(),
        TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), KaonCutsMCName.Data()));

    TString AntiKaonCutsMCName = Form("%sAntiKaonCutsMC_%s_%s", addon.Data(), suffixTrigger.Data(), suffix.Data());
    coutputAntiKaonCutsMC = mgr->CreateContainer(
        //@suppress("Invalid arguments") it works ffs
        AntiKaonCutsMCName.Data(),
        TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), AntiKaonCutsMCName.Data()));

  }

  TString ThreeBodyName = Form("%sProtonProtonKaons_%s_%s", addon.Data(), suffixTrigger.Data(), suffix.Data());
  AliAnalysisDataContainer *coutputThreeBody = mgr->CreateContainer(
    //@suppress("Invalid arguments") it works ffs
    ThreeBodyName.Data(),
    TList::Class(),
    AliAnalysisManager::kOutputContainer,
    Form("%s:%s", file.Data(), ThreeBodyName.Data()));



//    AliAnalysisTaskProtonProtonKaons* taskNano;
  AliAnalysisTaskThreeBodyProtonPrimary* taskNano;
//  AliAnalysisTaskProtonProtonKaonsAOD* taskAOD;
  if(isNano){
//    taskNano= new AliAnalysisTaskProtonProtonKaons("femtoNanoProtonProtonKaons", isMC);
    taskNano= new AliAnalysisTaskThreeBodyProtonPrimary("femtoNanoProtonProtonKaons", isMC);
    if (!fullBlastQA)
    {
      taskNano->SetRunTaskLightWeight(true);
    }

    taskNano->SetEventCuts(evtCuts);
    taskNano->SetProtonCuts(TrackCuts);
    taskNano->SetAntiProtonCuts(AntiTrackCuts);
    taskNano->SetPrimaryCuts(KaonCuts);
    taskNano->SetAntiPrimaryCuts(AntiKaonCuts);
    taskNano->SetCorrelationConfig(config);
    taskNano->SetRunThreeBodyHistograms(true);
    taskNano->SetClosePairRejectionForAll(ClosePairRejectionForAll);
    taskNano->SetturnoffClosePairRejectionCompletely(turnoffClosePairRejectionCompletely);
    taskNano->SetCleanWithLambdas(false);
    taskNano->SetDoOnlyThreeBody(DoOnlyThreeBody);
    taskNano->SetRunOfficialTwoBody(RunOfficialTwoBody);
    taskNano->SetRunPlotOtherHistos(false);
    //taskNano->SetDoTwoPrimary(DoTwoPrimary);
    taskNano->SetStandardMixing(StandardMixing);
    taskNano->SetRunPlotQ3Vsq(RunPlotQ3Vsq);
    taskNano->SetRunPlotPt(RunPlotPt);
    taskNano->SetDeltaPhiMaxPP(DeltaPhiMaxpp);      
    taskNano->SetDeltaEtaMaxPP(DeltaEtaMaxpp);
    taskNano->SetDeltaPhiMaxPPrim(DeltaPhiMaxpKplus);
    taskNano->SetDeltaEtaMaxPPrim(DeltaEtaMaxpKplus);
    taskNano->SetDeltaPhiMaxPAPrim(0.);
    taskNano->SetDeltaEtaMaxPAPrim(0.);
    taskNano->SetDoKinematicsPlots(DoKinematicPlots);
    taskNano->SetRunPlotInvMass(DoKinematicPlots);
    taskNano->SetRunPlotMult(RunPlotMult);
    taskNano->SetRunCorrDeltaPhi(RunCorrDeltaPhi);
    taskNano->SetRunPlotPhiTheta(RunPlotPhiTheta);
    taskNano->SetPlotP1(RunPlotP1);
    taskNano->SetCutElectrons(CutElectrons);
    if (isMC) taskNano->SetPlotsMC(true);

    mgr->AddTask(taskNano);

    mgr->ConnectInput(taskNano, 0, cinput);
    mgr->ConnectOutput(taskNano, 1, coutputEvtCuts);
    mgr->ConnectOutput(taskNano, 2, couputTrkCuts);
    mgr->ConnectOutput(taskNano, 3, coutputAntiTrkCuts);
    mgr->ConnectOutput(taskNano, 4, coutputKaonCuts);
    mgr->ConnectOutput(taskNano, 5, coutputAntiKaonCuts);
    mgr->ConnectOutput(taskNano, 8, coutputResults);
    mgr->ConnectOutput(taskNano, 9, coutputResultsQA);
    mgr->ConnectOutput(taskNano, 10, coutputResultsSample);
    mgr->ConnectOutput(taskNano, 11, coutputResultsSampleQA);
    mgr->ConnectOutput(taskNano, 12, coutputThreeBody);
    if (isMC) {
      mgr->ConnectOutput(taskNano, 13, coutputTrkCutsMC);
      mgr->ConnectOutput(taskNano, 14, coutputAntiTrkCutsMC);
      mgr->ConnectOutput(taskNano, 15, coutputKaonCutsMC);
      mgr->ConnectOutput(taskNano, 16, coutputAntiKaonCutsMC);
    }
  }

  else{
/*
    taskAOD= new AliAnalysisTaskProtonProtonKaonsAOD("femtoAODProtonProtonKaons", isMC, triggerOn);
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
    taskAOD->Setv0Cuts(KaonCuts);
    taskAOD->SetAntiv0Cuts(AntiKaonCuts);
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
    mgr->ConnectOutput(taskAOD, 4, coutputKaonCuts);
    mgr->ConnectOutput(taskAOD, 5, coutputAntiKaonCuts);
    mgr->ConnectOutput(taskAOD, 6, coutputResults);
    mgr->ConnectOutput(taskAOD, 7, coutputResultsQA);
    mgr->ConnectOutput(taskAOD, 8, coutputResultsSample);
    mgr->ConnectOutput(taskAOD, 9, coutputResultsSampleQA);
    mgr->ConnectOutput(taskAOD, 10, coutputThreeBody);
    if (isMC) {
      mgr->ConnectOutput(taskAOD, 16, coutputTrkCutsMC);
      mgr->ConnectOutput(taskAOD, 17, coutputAntiTrkCutsMC);
      mgr->ConnectOutput(taskAOD, 18, coutputKaonCutsMC);
      mgr->ConnectOutput(taskAOD, 19, coutputAntiKaonCutsMC);
    }
*/
  }


  if (isNano) {
    return taskNano;
  } else {
  //  return taskAOD;
  }


}
