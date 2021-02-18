#include <vector>
//#include "AliAnalysisTaskSE.h"
//#include "AliAnalysisManager.h"
//#include "AliAnalysisTaskGeorgiosNTuple.h"
//#include "AliFemtoDreamEventCuts.h"
//#include "AliFemtoDreamTrackCuts.h"
//#include "AliFemtoDreamCascadeCuts.h"
//#include "AliFemtoDreamCollConfig.h"

#define MONTECARLO

AliAnalysisTaskSE *AddTaskGeorgiosNTuple(bool fullBlastQA = true,
		                         bool isMC=true,
					 const char *cutVariation = "0") {
  //set fullBlastQA and suffix (cut variation)
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
  evtCuts->SetMultVsCentPlots(true);

//v0 Cuts (Georgios)
  AliFemtoDreamv0Cuts *v0Cuts = AliFemtoDreamv0Cuts::LambdaCuts(isMC, true, false);
  AliFemtoDreamTrackCuts *Posv0Daug = AliFemtoDreamTrackCuts::DecayProtonCuts(isMC, true, false);//PileUpRej, false
  AliFemtoDreamTrackCuts *Negv0Daug = AliFemtoDreamTrackCuts::DecayPionCuts(isMC, true, false);
  v0Cuts->SetPosDaugterTrackCuts(Posv0Daug);
  v0Cuts->SetNegDaugterTrackCuts(Negv0Daug);
  v0Cuts->SetPDGCodePosDaug(2212);  //Proton
  v0Cuts->SetPDGCodeNegDaug(-211);  //Pion
  v0Cuts->SetPDGCodev0(3122);  //Lambda
  v0Cuts->SetCutInvMass(0.03); //include Background
//Anti v0 Cuts
  AliFemtoDreamv0Cuts *Antiv0Cuts = AliFemtoDreamv0Cuts::LambdaCuts(isMC, true, false);
  AliFemtoDreamTrackCuts *PosAntiv0Daug = AliFemtoDreamTrackCuts::DecayPionCuts(isMC, true, false);
  PosAntiv0Daug->SetCutCharge(1);
  AliFemtoDreamTrackCuts *NegAntiv0Daug = AliFemtoDreamTrackCuts::DecayProtonCuts(isMC, true, false);
  NegAntiv0Daug->SetCutCharge(-1);
  Antiv0Cuts->SetPosDaugterTrackCuts(PosAntiv0Daug);
  Antiv0Cuts->SetNegDaugterTrackCuts(NegAntiv0Daug);
  Antiv0Cuts->SetPDGCodePosDaug(211);  //Pion
  Antiv0Cuts->SetPDGCodeNegDaug(-2212);  //Proton
  Antiv0Cuts->SetPDGCodev0(-3122);  //Lambda
  Antiv0Cuts->SetCutInvMass(0.03); //include Background

  //Systematic variations for v0
  //v0Cuts->SetCutDCADaugTov0Vtx(1.9);
  //Antiv0Cuts->SetCutDCADaugTov0Vtx(1.9);
  v0Cuts->SetCutDCADaugToPrimVtx(0.05);
  Antiv0Cuts->SetCutDCADaugToPrimVtx(0.05);
  v0Cuts->SetPtRange(0.22,999.);
  Antiv0Cuts->SetPtRange(0.22,999.);

  //Cascade Cuts (Background)
  AliFemtoDreamCascadeCuts* CascadeXiBGRCuts = AliFemtoDreamCascadeCuts::XiCuts(isMC, false);
  CascadeXiBGRCuts->SetXiCharge(-1);
  AliFemtoDreamTrackCuts *XiBGRNegCuts = AliFemtoDreamTrackCuts::Xiv0PionCuts(isMC, true, false);
  XiBGRNegCuts->SetCheckTPCRefit(false);//for nanos this is already done while prefiltering
  AliFemtoDreamTrackCuts *XiBGRPosCuts = AliFemtoDreamTrackCuts::Xiv0ProtonCuts(isMC, true, false);
  XiBGRPosCuts->SetCheckTPCRefit(false);
  AliFemtoDreamTrackCuts *XiBGRBachCuts = AliFemtoDreamTrackCuts::XiBachPionCuts(isMC, true, false);
  XiBGRBachCuts->SetCheckTPCRefit(false);
  CascadeXiBGRCuts->Setv0Negcuts(XiBGRNegCuts);
  CascadeXiBGRCuts->Setv0PosCuts(XiBGRPosCuts);
  CascadeXiBGRCuts->SetBachCuts(XiBGRBachCuts);
  CascadeXiBGRCuts->SetPDGCodeCasc(3312);                
  CascadeXiBGRCuts->SetPDGCodev0(3122);
  CascadeXiBGRCuts->SetPDGCodePosDaug(2212);
  CascadeXiBGRCuts->SetPDGCodeNegDaug(-211);
  CascadeXiBGRCuts->SetPDGCodeBach(-211);
  //AntiCascade cuts (Background)
  AliFemtoDreamCascadeCuts* AntiCascadeXiBGRCuts = AliFemtoDreamCascadeCuts::XiCuts(isMC, false);
  AntiCascadeXiBGRCuts->SetXiCharge(1);
  AliFemtoDreamTrackCuts *AntiXiBGRNegCuts = AliFemtoDreamTrackCuts::Xiv0ProtonCuts(isMC, true, false);
  AntiXiBGRNegCuts->SetCutCharge(-1);
  AntiXiBGRNegCuts->SetCheckTPCRefit(false);
  AliFemtoDreamTrackCuts *AntiXiBGRPosCuts = AliFemtoDreamTrackCuts::Xiv0PionCuts(isMC, true, false);
  AntiXiBGRPosCuts->SetCutCharge(1);
  AntiXiBGRPosCuts->SetCheckTPCRefit(false);
  AliFemtoDreamTrackCuts *AntiXiBGRBachCuts =AliFemtoDreamTrackCuts::XiBachPionCuts(isMC, true, false);
  AntiXiBGRBachCuts->SetCutCharge(1);
  AntiXiBGRBachCuts->SetCheckTPCRefit(false);
  AntiCascadeXiBGRCuts->Setv0Negcuts(AntiXiBGRNegCuts);
  AntiCascadeXiBGRCuts->Setv0PosCuts(AntiXiBGRPosCuts);
  AntiCascadeXiBGRCuts->SetBachCuts(AntiXiBGRBachCuts);
  AntiCascadeXiBGRCuts->SetPDGCodeCasc(-3312);
  AntiCascadeXiBGRCuts->SetPDGCodev0(-3122);
  AntiCascadeXiBGRCuts->SetPDGCodePosDaug(211);
  AntiCascadeXiBGRCuts->SetPDGCodeNegDaug(-2212);
  AntiCascadeXiBGRCuts->SetPDGCodeBach(211); 

  //Cascade Cuts 
  AliFemtoDreamCascadeCuts* CascadeXiCuts = AliFemtoDreamCascadeCuts::XiCuts(isMC, false);
  CascadeXiCuts->SetXiCharge(-1);
  AliFemtoDreamTrackCuts *XiNegCuts = AliFemtoDreamTrackCuts::Xiv0PionCuts(isMC, true, false);
  XiNegCuts->SetCheckTPCRefit(false);
  AliFemtoDreamTrackCuts *XiPosCuts = AliFemtoDreamTrackCuts::Xiv0ProtonCuts(isMC, true, false);
  XiPosCuts->SetCheckTPCRefit(false);
  AliFemtoDreamTrackCuts *XiBachCuts = AliFemtoDreamTrackCuts::XiBachPionCuts(isMC, true, false);
  XiBachCuts->SetCheckTPCRefit(false);
  CascadeXiCuts->Setv0Negcuts(XiNegCuts);
  CascadeXiCuts->Setv0PosCuts(XiPosCuts);
  CascadeXiCuts->SetBachCuts(XiBachCuts);
  CascadeXiCuts->SetPDGCodeCasc(3312);
  CascadeXiCuts->SetPDGCodev0(3122);
  CascadeXiCuts->SetPDGCodePosDaug(2212);
  CascadeXiCuts->SetPDGCodeNegDaug(-211);
  CascadeXiCuts->SetPDGCodeBach(-211);
  CascadeXiCuts->SetXiMassRange(1.322, 0.06); //include Background
  

  //AntiCascade cuts 
  AliFemtoDreamCascadeCuts* AntiCascadeXiCuts = AliFemtoDreamCascadeCuts::XiCuts(isMC, false);
  AntiCascadeXiCuts->SetXiCharge(1);
  AliFemtoDreamTrackCuts *AntiXiNegCuts = AliFemtoDreamTrackCuts::Xiv0ProtonCuts(isMC, true, false);
  AntiXiNegCuts->SetCutCharge(-1);
  AntiXiNegCuts->SetCheckTPCRefit(false);
  AliFemtoDreamTrackCuts *AntiXiPosCuts = AliFemtoDreamTrackCuts::Xiv0PionCuts(isMC, true, false);
  AntiXiPosCuts->SetCutCharge(1);
  AntiXiPosCuts->SetCheckTPCRefit(false);
  AliFemtoDreamTrackCuts *AntiXiBachCuts =AliFemtoDreamTrackCuts::XiBachPionCuts(isMC, true, false);
  AntiXiBachCuts->SetCutCharge(1);
  AntiXiBachCuts->SetCheckTPCRefit(false);
  AntiCascadeXiCuts->Setv0Negcuts(AntiXiNegCuts);
  AntiCascadeXiCuts->Setv0PosCuts(AntiXiPosCuts);
  AntiCascadeXiCuts->SetBachCuts(AntiXiBachCuts);
  AntiCascadeXiCuts->SetPDGCodeCasc(-3312);
  AntiCascadeXiCuts->SetPDGCodev0(-3122);
  AntiCascadeXiCuts->SetPDGCodePosDaug(211);
  AntiCascadeXiCuts->SetPDGCodeNegDaug(-2212);
  AntiCascadeXiCuts->SetPDGCodeBach(211);
  AntiCascadeXiCuts->SetXiMassRange(1.322, 0.06); //include Background

  //systematics for Xis
  CascadeXiCuts->SetCutXiDaughterDCA(1.9);
  AntiCascadeXiCuts->SetCutXiDaughterDCA(1.9);
  CascadeXiCuts->SetCutXiMinDistBachToPrimVtx(0.04);
  AntiCascadeXiCuts->SetCutXiMinDistBachToPrimVtx(0.04);
  CascadeXiCuts->SetCutv0MinDistToPrimVtx(0.06);
  AntiCascadeXiCuts->SetCutv0MinDistToPrimVtx(0.06);
  CascadeXiCuts->SetCutv0MinDaugDistToPrimVtx(0.04);
  AntiCascadeXiCuts->SetCutv0MinDaugDistToPrimVtx(0.04);
  CascadeXiCuts->SetCutXiTransverseRadius(0.6, 200);
  AntiCascadeXiCuts->SetCutXiTransverseRadius(0.6, 200);
  CascadeXiCuts->SetCutv0TransverseRadius(1.1, 200);
  AntiCascadeXiCuts->SetCutv0TransverseRadius(1.1, 200);
  XiNegCuts->SetEtaRange(-0.91, 0.91);
  XiPosCuts->SetEtaRange(-0.91, 0.91);
  XiBachCuts->SetEtaRange(-0.91, 0.91);
  AntiXiNegCuts->SetEtaRange(-0.91, 0.91);
  AntiXiPosCuts->SetEtaRange(-0.91, 0.91);
  AntiXiBachCuts->SetEtaRange(-0.91, 0.91);

  //BGR 
  CascadeXiBGRCuts->SetCutXiDaughterDCA(1.9);
  AntiCascadeXiBGRCuts->SetCutXiDaughterDCA(1.9);
  CascadeXiBGRCuts->SetCutXiMinDistBachToPrimVtx(0.04);
  AntiCascadeXiBGRCuts->SetCutXiMinDistBachToPrimVtx(0.04);
  CascadeXiBGRCuts->SetCutv0MinDistToPrimVtx(0.06);
  AntiCascadeXiBGRCuts->SetCutv0MinDistToPrimVtx(0.06);   
  CascadeXiBGRCuts->SetCutv0MinDaugDistToPrimVtx(0.04);
  AntiCascadeXiBGRCuts->SetCutv0MinDaugDistToPrimVtx(0.04);
  CascadeXiBGRCuts->SetCutXiTransverseRadius(0.6, 200);
  AntiCascadeXiBGRCuts->SetCutXiTransverseRadius(0.6, 200);
  CascadeXiBGRCuts->SetCutv0TransverseRadius(1.1, 200);
  AntiCascadeXiBGRCuts->SetCutv0TransverseRadius(1.1, 200);
  XiBGRNegCuts->SetEtaRange(-0.91, 0.91);
  XiBGRPosCuts->SetEtaRange(-0.91, 0.91);
  XiBGRBachCuts->SetEtaRange(-0.91, 0.91);
  AntiXiBGRNegCuts->SetEtaRange(-0.91, 0.91);
  AntiXiBGRPosCuts->SetEtaRange(-0.91, 0.91);
  AntiXiBGRBachCuts->SetEtaRange(-0.91, 0.91);


  if (suffix != "0" && suffix != "999") {
    evtCuts->SetMinimalBooking(true);
    v0Cuts->SetMinimalBooking(true);
    Antiv0Cuts->SetMinimalBooking(true);
    CascadeXiBGRCuts->SetMinimalBooking(true);
    AntiCascadeXiBGRCuts->SetMinimalBooking(true);
    CascadeXiCuts->SetMinimalBooking(true);
    AntiCascadeXiCuts->SetMinimalBooking(true);
  }

  AliFemtoDreamCollConfig *config = new AliFemtoDreamCollConfig("Femto","Femto");


  // Femto Collection
  std::vector<int> PDGParticles;   
  PDGParticles.push_back(3122);
  PDGParticles.push_back(3122);
  PDGParticles.push_back(3312);
  PDGParticles.push_back(3312);
  PDGParticles.push_back(3312);
  PDGParticles.push_back(3312);



  std::vector<int> NBins;
  std::vector<float> kMin;
  std::vector<float> kMax;
  std::vector<int> pairQA;
  std::vector<bool> closeRejection;
  //pairs:
  //Lambda Lambda               0
  //Lambda barLambda            1
  //Lambda XiBGR                2
  //Lambda barXiBGR             3
  //Lambda Xi                   4
  //Lambda barXi                5
  //barLambda barLambda         6
  //barLambda XiBGR             7
  //barLambda barXiBGR          8
  //barp Xi                     9
  //barp barXi                  10
  //XiBGR XiBGR                 11
  //XiBGR barXiBGR              12
  //XiBGR Xi                    13
  //XiBGR barXi                 14
  //barXiBGR barXiBGR           15
  //barXiBGR Xi                 16
  //barXiBGR barXi              17
  //Xi Xi                       18
  //Xi barXi                    19
  //barXi barXi                 20


  const int nPairs = 21;
  for (int i = 0; i < nPairs; ++i) {
    pairQA.push_back(0);
    closeRejection.push_back(false);
    if (suffix == "0") {
      NBins.push_back(1500);
      kMin.push_back(0.);
      kMax.push_back(6.);
    } else {
      NBins.push_back(250);
      kMin.push_back(0.);
      kMax.push_back(1.);
    }
  }

  pairQA[4] = 13;//p-Omega  //to be changed 
  pairQA[10] = 13;//pbar-antiOmega  //to be changed

  closeRejection[0] = true;  // pp             //to be changed
  closeRejection[6] = true;  // barp barp      //to be changed

  config->SetPDGCodes(PDGParticles);
  config->SetNBinsHist(NBins);
  config->SetMinKRel(kMin);
  config->SetMaxKRel(kMax);
  config->SetClosePairRejection(closeRejection);
  config->SetDeltaEtaMax(0.012);
  config->SetDeltaPhiMax(0.012);
  config->SetExtendedQAPairs(pairQA);

  config->SetMixingDepth(10);
  config->SetUseEventMixing(true);
  config->SetMultiplicityEstimator(AliFemtoDreamEvent::kRef08);

  if (isMC) {
   config->SetMomentumResolution(true);
  }

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
  config->SetdPhidEtaPlotsSmallK(false);
  config->SetdPhidEtaPlots(false);

  config->SetPhiEtaBinnign(false);

  if (suffix == "0" && fullBlastQA) {
    config->SetkTBinning(true);
    config->SetmTBinning(true);
    config->SetPtQA(true);
  }

  if (suffix != "0") {
    config->SetMinimalBookingME(true);
  }
  AliAnalysisTaskGeorgiosNTuple* task = new AliAnalysisTaskGeorgiosNTuple("GeorgiosNTuple",true);
  if (suffix != "0" && suffix != "999") {
    task->SetRunTaskLightWeight(true);
  }
  task->SelectCollisionCandidates(AliVEvent::kHighMultV0);
  task->SetEventCuts(evtCuts);
  task->SetLambdaCuts(v0Cuts);
  task->SetAntiLambdaCuts(Antiv0Cuts);
  task->SetXiBGRCuts(CascadeXiBGRCuts);
  task->SetAntiXiBGRCuts(AntiCascadeXiBGRCuts);
  task->SetXiCuts(CascadeXiCuts);
  task->SetAntiXiCuts(AntiCascadeXiCuts);
  task->SetCorrelationConfig(config);
  mgr->AddTask(task);

  TString addon = "LambdaXi";

  TString file = AliAnalysisManager::GetCommonFileName();

  mgr->ConnectInput(task, 0, cinput);


  TString EvtCutsName = Form("%sEvtCuts%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *coutputEvtCuts = mgr->CreateContainer(
      EvtCutsName.Data(), TList::Class(), AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), EvtCutsName.Data()));
  mgr->ConnectOutput(task, 1, coutputEvtCuts);

  TString v0CutsName = Form("%sv0Cuts%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *couputv0Cuts = mgr->CreateContainer(
      v0CutsName.Data(), TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), v0CutsName.Data()));
  mgr->ConnectOutput(task, 2, couputv0Cuts);

  TString Antiv0CutsName = Form("%sAntiv0Cuts%s", addon.Data(),
                                   suffix.Data());
  AliAnalysisDataContainer *coutputAntiv0Cuts = mgr->CreateContainer(
      Antiv0CutsName.Data(), TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), Antiv0CutsName.Data()));
  mgr->ConnectOutput(task, 3, coutputAntiv0Cuts);

  AliAnalysisDataContainer *coutputCascadeXiBGRCuts;
  TString CascadeXiBGRCutsName = Form("%sCascadeXiBGRCuts%s", addon.Data(), suffix.Data());
  coutputCascadeXiBGRCuts = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      CascadeXiBGRCutsName.Data(),
      TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), CascadeXiBGRCutsName.Data()));
  mgr->ConnectOutput(task, 4, coutputCascadeXiBGRCuts);


  AliAnalysisDataContainer *coutputAntiCascadeXiBGRCuts;
  TString AntiCascadeXiBGRCutsName = Form("%sAntiCascadeXiBGRCuts%s", addon.Data(), suffix.Data());
  coutputAntiCascadeXiBGRCuts = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      AntiCascadeXiBGRCutsName.Data(),
      TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), AntiCascadeXiBGRCutsName.Data()));
  mgr->ConnectOutput(task, 5, coutputAntiCascadeXiBGRCuts);


  AliAnalysisDataContainer *coutputCascadeXiCuts;
  TString CascadeXiCutsName = Form("%sCascadeXiCuts%s", addon.Data(), suffix.Data());
  coutputCascadeXiCuts = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      CascadeXiCutsName.Data(),
      TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), CascadeXiCutsName.Data()));
  mgr->ConnectOutput(task, 6, coutputCascadeXiCuts);


  AliAnalysisDataContainer *coutputAntiCascadeXiCuts;
  TString AntiCascadeXiCutsName = Form("%sAntiCascadeXiCuts%s", addon.Data(), suffix.Data());
  coutputAntiCascadeXiCuts = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      AntiCascadeXiCutsName.Data(),
      TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), AntiCascadeXiCutsName.Data()));
  mgr->ConnectOutput(task, 7, coutputAntiCascadeXiCuts);









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


  //Georgios tree:
  AliAnalysisDataContainer *coutputTreeGeorgos;
  TString TreeGeorgiosName = Form("%sTreeGeorgios",addon.Data());
  coutputTreeGeorgios = mgr->CreateContainer(
    //@suppress("Invalid arguments") it works ffs
    TreeGeorgiosName.Data(),
    TTree::Class(),
    AliAnalysisManager::kOutputContainer,
    Form("%s:%s", file.Data(), TreeGeorgiosName.Data()));
  mgr->ConnectOutput(task, 10, coutputTreeGeorgios);


#ifdef MONTECARLO

  AliAnalysisDataContainer *coutputv0CutsMC;
  AliAnalysisDataContainer *coutputAntiv0CutsMC;
  AliAnalysisDataContainer *coutputCascadeXiBGRCutsMC;
  AliAnalysisDataContainer *coutputAntiCascadeXiBGRCutsMC;
  AliAnalysisDataContainer *coutputCascadeXiCutsMC;
  AliAnalysisDataContainer *coutputAntiCascadeXiCutsMC;
  
    TString v0CutsMCName = Form("%sv0CutsMC%s", addon.Data(), suffix.Data());
    coutputv0CutsMC = mgr->CreateContainer(
        //@suppress("Invalid arguments") it works ffs
        v0CutsMCName.Data(),
        TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), v0CutsMCName.Data()));
    mgr->ConnectOutput(task, 11, coutputv0CutsMC);

    TString Antiv0CutsMCName = Form("%sAntiv0CutsMC%s", addon.Data(), suffix.Data());
    coutputAntiv0CutsMC = mgr->CreateContainer(
        //@suppress("Invalid arguments") it works ffs
        Antiv0CutsMCName.Data(),
        TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), Antiv0CutsMCName.Data()));
    mgr->ConnectOutput(task, 12, coutputAntiv0CutsMC);

   TString CascadeXiBGRCutsMCName = Form("%sCascadeXiBGRCutsMC%s", addon.Data(), suffix.Data());
    coutputCascadeXiBGRCutsMC = mgr->CreateContainer(
        //@suppress("Invalid arguments") it works ffs
        CascadeXiBGRCutsMCName.Data(),
        TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), CascadeXiBGRCutsMCName.Data()));
    mgr->ConnectOutput(task, 13, coutputCascadeXiBGRCutsMC);

    TString AntiCascadeXiBGRCutsMCName = Form("%sAntiCascadeXiBGRCutsMC%s", addon.Data(), suffix.Data());
    coutputAntiCascadeXiBGRCutsMC = mgr->CreateContainer(
        //@suppress("Invalid arguments") it works ffs
        AntiCascadeXiBGRCutsMCName.Data(),
        TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), AntiCascadeXiBGRCutsMCName.Data()));
    mgr->ConnectOutput(task, 14, coutputAntiCascadeXiBGRCutsMC);

    TString CascadeXiCutsMCName = Form("%sCascadeXiCutsMC%s", addon.Data(), suffix.Data());
    coutputCascadeXiCutsMC = mgr->CreateContainer(
        //@suppress("Invalid arguments") it works ffs
        CascadeXiCutsMCName.Data(),
        TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), CascadeXiCutsMCName.Data()));
    mgr->ConnectOutput(task, 15, coutputCascadeXiCutsMC);

    TString AntiCascadeXiCutsMCName = Form("%sAntiCascadeXiCutsMC%s", addon.Data(), suffix.Data());
    coutputAntiCascadeXiCutsMC = mgr->CreateContainer(
        //@suppress("Invalid arguments") it works ffs
        AntiCascadeXiCutsMCName.Data(),
        TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), AntiCascadeXiCutsMCName.Data()));
    mgr->ConnectOutput(task, 16, coutputAntiCascadeXiCutsMC);
#endif

    return task;
}
