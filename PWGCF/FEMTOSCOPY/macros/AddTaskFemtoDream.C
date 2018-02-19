#ifndef ADDTASKFEMTODREAM_C
#define ADDTASKFEMTODREAM_C
#include "AliFemtoDreamEventCuts.h"
#include "AliFemtoDreamTrackCuts.h"
#include "AliAnalysisTaskFemtoDream.h"
#include "TROOT.h"

void AddTaskFemtoDream(bool isMC, TString CentEst) {
  gROOT->ProcessLine(".include $ALICE_ROOT/include");
  gROOT->ProcessLine(".include $ALICE_PHYSICS/include");
  gROOT->ProcessLine(".include $ROOTSYS/include");

  // the manager is static, so get the existing manager via the static method
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  if (!mgr)
  {
    printf("No analysis manager to connect to!\n");
    return;
  }
  // just to see if all went well, check if the input event handler has been
  // connected
  if (!mgr->GetInputEventHandler()) {
    printf("This task requires an input event handler!\n");
    return;
  }
  AliFemtoDreamEventCuts *evtCuts=
      AliFemtoDreamEventCuts::StandardCutsRun1();
  bool DCAPlots=false;
  bool CPAPlots=false;
  bool CombSigma=false;
  bool ContributionSplitting=false;
  bool ContributionSplittingDaug=false;

  //Track Cuts
  AliFemtoDreamTrackCuts *TrackCuts=
      AliFemtoDreamTrackCuts::PrimProtonCuts(
          isMC,DCAPlots,CombSigma,ContributionSplitting);
  TrackCuts->SetCutCharge(1);
  AliFemtoDreamTrackCuts *AntiTrackCuts=
      AliFemtoDreamTrackCuts::PrimProtonCuts(
          isMC,DCAPlots,CombSigma,ContributionSplitting);
  AntiTrackCuts->SetCutCharge(-1);

  //Lambda Cuts
  AliFemtoDreamv0Cuts *v0Cuts=
      AliFemtoDreamv0Cuts::LambdaCuts(isMC,CPAPlots,ContributionSplitting);
  AliFemtoDreamTrackCuts *Posv0Daug=AliFemtoDreamTrackCuts::DecayProtonCuts(
      isMC,ContributionSplittingDaug);
  AliFemtoDreamTrackCuts *Negv0Daug=AliFemtoDreamTrackCuts::DecayPionCuts(
      isMC,ContributionSplittingDaug);
  v0Cuts->SetPosDaugterTrackCuts(Posv0Daug);
  v0Cuts->SetNegDaugterTrackCuts(Negv0Daug);
  v0Cuts->SetPDGCodePosDaug(2212);//Proton
  v0Cuts->SetPDGCodeNegDaug(211);//Pion
  v0Cuts->SetPDGCodev0(3122);//Lambda
  AliFemtoDreamv0Cuts *Antiv0Cuts=
      AliFemtoDreamv0Cuts::LambdaCuts(isMC,CPAPlots,ContributionSplitting);
  AliFemtoDreamTrackCuts *PosAntiv0Daug=AliFemtoDreamTrackCuts::DecayPionCuts(
      isMC,ContributionSplittingDaug);
  PosAntiv0Daug->SetCutCharge(1);
  AliFemtoDreamTrackCuts *NegAntiv0Daug=AliFemtoDreamTrackCuts::DecayProtonCuts(
      isMC,ContributionSplittingDaug);
  NegAntiv0Daug->SetCutCharge(-1);
  Antiv0Cuts->SetPosDaugterTrackCuts(PosAntiv0Daug);
  Antiv0Cuts->SetNegDaugterTrackCuts(NegAntiv0Daug);
  Antiv0Cuts->SetPDGCodePosDaug(211);//Pion
  Antiv0Cuts->SetPDGCodeNegDaug(2212);//Proton
  Antiv0Cuts->SetPDGCodev0(3122);//Lambda

  //Cascade Cuts
  AliFemtoDreamCascadeCuts *CascadeCuts=
      AliFemtoDreamCascadeCuts::XiCuts(false);
  CascadeCuts->SetXiCharge(-1);
  AliFemtoDreamTrackCuts *XiNegCuts=
      AliFemtoDreamTrackCuts::Xiv0PionCuts(false,false);
  AliFemtoDreamTrackCuts *XiPosCuts=
      AliFemtoDreamTrackCuts::Xiv0ProtonCuts(false,false);
  AliFemtoDreamTrackCuts *XiBachCuts=
      AliFemtoDreamTrackCuts::XiBachPionCuts(false,false);
  AliFemtoDreamCascadeCuts *AntiCascadeCuts=
      AliFemtoDreamCascadeCuts::XiCuts(false);
  CascadeCuts->Setv0Negcuts(XiNegCuts);
  CascadeCuts->Setv0PosCuts(XiPosCuts);
  CascadeCuts->SetBachCuts(XiBachCuts);
  AntiCascadeCuts->SetXiCharge(1);

  AliFemtoDreamTrackCuts *AntiXiNegCuts=
      AliFemtoDreamTrackCuts::Xiv0ProtonCuts(false,false);
  AntiXiNegCuts->SetCutCharge(-1);
  AliFemtoDreamTrackCuts *AntiXiPosCuts=
      AliFemtoDreamTrackCuts::Xiv0PionCuts(false,false);
  AntiXiPosCuts->SetCutCharge(1);
  AliFemtoDreamTrackCuts *AntiXiBachCuts=
      AliFemtoDreamTrackCuts::XiBachPionCuts(false,false);
  AntiCascadeCuts->Setv0Negcuts(AntiXiNegCuts);
  AntiCascadeCuts->Setv0PosCuts(AntiXiPosCuts);
  AntiCascadeCuts->SetBachCuts(AntiXiBachCuts);
  AntiXiBachCuts->SetCutCharge(1);

  std::vector<int> PDGParticles ={2212,2212,3122,3122};
  std::vector<double> ZVtxBins = {-10,-8,-6,-4,-2,0,2,4,6,8,10};
  std::vector<int> MultBins = {0,4,8,12,16,20,24,28,32,36,40,44,48,52,56,60,64,80};
  std::vector<int> NBins={750,750,150,150,750,150,150,150,150,150};
  std::vector<double> kMin={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  std::vector<double> kMax={3.,3.,3.,3.,3.,3.,3.,3.,3.,3.};
  AliFemtoDreamCollConfig *config=new AliFemtoDreamCollConfig("Femto","Femto");
  config->SetZBins(ZVtxBins);
  config->SetMultBins(MultBins);
  config->SetPDGCodes(PDGParticles);
  config->SetNBinsHist(NBins);
  config->SetMinKRel(kMin);
  config->SetMaxKRel(kMax);
  config->SetMixingDepth(10);

  AliAnalysisTaskFemtoDream *task=new AliAnalysisTaskFemtoDream("FemtoDream",isMC);
  if(CentEst == "kInt7"){
    task->SelectCollisionCandidates(AliVEvent::kINT7);
    task->SetMVPileUp(kTRUE);
  }else if(CentEst == "kMB"){
    task->SelectCollisionCandidates(AliVEvent::kMB);
    task->SetMVPileUp(kFALSE);
  }else{
    std::cout << "=====================================================================" << std::endl;
    std::cout << "=====================================================================" << std::endl;
    std::cout << "Centrality Estimator not set, fix it else your Results will be empty!" << std::endl;
    std::cout << "=====================================================================" << std::endl;
    std::cout << "=====================================================================" << std::endl;
  }
  task->SetDebugLevel(0);
  task->SetEvtCutQA(true);
  task->SetTrackBufferSize(2500);
  task->SetEventCuts(evtCuts);
  task->SetTrackCuts(TrackCuts);
  task->SetAntiTrackCuts(AntiTrackCuts);
  task->Setv0Cuts(v0Cuts);
  task->SetAntiv0Cuts(Antiv0Cuts);
  task->SetCascadeCuts(CascadeCuts);
  task->SetAntiCascadeCuts(AntiCascadeCuts);
  task->SetCollectionConfig(config);
  mgr->AddTask(task);
//  TString QAOutput="AnalysisQA.root";
//  TString TrackOutput="AnalysisQATracks.root";
//  TString v0Output="AnalysisQAv0.root";
//  TString CascadeOutput="AnalysisQACascade.root";

  TString file = AliAnalysisManager::GetCommonFileName();

  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();

  mgr->ConnectInput(task, 0, cinput);

  AliAnalysisDataContainer *coutputQA;
  TString QAName = Form("QA");
  coutputQA = mgr->CreateContainer(//@suppress("Invalid arguments") it works ffs
      QAName.Data(), TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), QAName.Data()));
  mgr->ConnectOutput(task, 1, coutputQA);

  AliAnalysisDataContainer *coutputEvtCuts;
  TString EvtCutsName = Form("EvtCuts");
  coutputEvtCuts = mgr->CreateContainer(//@suppress("Invalid arguments") it works ffs
      EvtCutsName.Data(), TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), EvtCutsName.Data()));
  mgr->ConnectOutput(task, 2, coutputEvtCuts);

  AliAnalysisDataContainer *couputTrkCuts;
  TString TrackCutsName = Form("TrackCuts");
  couputTrkCuts = mgr->CreateContainer(//@suppress("Invalid arguments") it works ffs
      TrackCutsName.Data(), TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), TrackCutsName.Data()));
  mgr->ConnectOutput(task, 3, couputTrkCuts);

  AliAnalysisDataContainer *coutputAntiTrkCuts;
  TString AntiTrackCutsName = Form("AntiTrackCuts");
  coutputAntiTrkCuts = mgr->CreateContainer(//@suppress("Invalid arguments") it works ffs
      AntiTrackCutsName.Data(), TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), AntiTrackCutsName.Data()));
  mgr->ConnectOutput(task, 4, coutputAntiTrkCuts);

  AliAnalysisDataContainer *coutputv0Cuts;
  TString v0CutsName = Form("v0Cuts");
  coutputv0Cuts = mgr->CreateContainer(//@suppress("Invalid arguments") it works ffs
      v0CutsName.Data(), TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), v0CutsName.Data()));
  mgr->ConnectOutput(task, 5, coutputv0Cuts);

  AliAnalysisDataContainer *coutputAntiv0Cuts;
  TString Antiv0CutsName = Form("Antiv0Cuts");
  coutputAntiv0Cuts = mgr->CreateContainer(//@suppress("Invalid arguments") it works ffs
      Antiv0CutsName.Data(), TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), Antiv0CutsName.Data()));
  mgr->ConnectOutput(task, 6, coutputAntiv0Cuts);

  AliAnalysisDataContainer *coutputCascadeCuts;
  TString CascadeCutsName = Form("CascadeCuts");
  coutputCascadeCuts = mgr->CreateContainer(//@suppress("Invalid arguments") it works ffs
      CascadeCutsName.Data(), TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), CascadeCutsName.Data()));
  mgr->ConnectOutput(task, 7, coutputCascadeCuts);

  AliAnalysisDataContainer *coutputAntiCascadeCuts;
  TString AntiCascadeCutsName = Form("AntiCascadeCuts");
  coutputAntiCascadeCuts = mgr->CreateContainer(//@suppress("Invalid arguments") it works ffs
      AntiCascadeCutsName.Data(), TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), AntiCascadeCutsName.Data()));
  mgr->ConnectOutput(task, 8, coutputAntiCascadeCuts);

  AliAnalysisDataContainer *coutputResults;
  TString ResultsName = Form("Results");
  coutputResults = mgr->CreateContainer(//@suppress("Invalid arguments") it works ffs
      ResultsName.Data(), TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), ResultsName.Data()));
  mgr->ConnectOutput(task, 9, coutputResults);

  AliAnalysisDataContainer *coutputResultQA;
  TString ResultQAName = Form("ResultQA");
  coutputResultQA = mgr->CreateContainer(//@suppress("Invalid arguments") it works ffs
      ResultQAName.Data(), TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), ResultQAName.Data()));
  mgr->ConnectOutput(task, 10, coutputResultQA);
  if (isMC) {
    AliAnalysisDataContainer *coutputTrkCutsMC;
    TString TrkCutsMCName = Form("TrkCutsMC");
    coutputTrkCutsMC = mgr->CreateContainer(//@suppress("Invalid arguments") it works ffs
        TrkCutsMCName.Data(), TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), TrkCutsMCName.Data()));
    mgr->ConnectOutput(task, 11, coutputTrkCutsMC);

    AliAnalysisDataContainer *coutputAntiTrkCutsMC;
    TString AntiTrkCutsMCName = Form("AntiTrkCutsMC");
    coutputAntiTrkCutsMC = mgr->CreateContainer(//@suppress("Invalid arguments") it works ffs
        AntiTrkCutsMCName.Data(), TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), AntiTrkCutsMCName.Data()));
    mgr->ConnectOutput(task, 12, coutputAntiTrkCutsMC);

    AliAnalysisDataContainer *coutputv0CutsMC;
    TString v0CutsMCName = Form("v0CutsMC");
    coutputv0CutsMC = mgr->CreateContainer(//@suppress("Invalid arguments") it works ffs
        v0CutsMCName.Data(), TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), v0CutsMCName.Data()));
    mgr->ConnectOutput(task, 13, coutputv0CutsMC);

    AliAnalysisDataContainer *coutputAntiv0CutsMC;
    TString Antiv0CutsMCName = Form("Antiv0CutsMC");
    coutputAntiv0CutsMC = mgr->CreateContainer(//@suppress("Invalid arguments") it works ffs
        Antiv0CutsMCName.Data(), TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), Antiv0CutsMCName.Data()));
    mgr->ConnectOutput(task, 14, coutputAntiv0CutsMC);
  }
  if (!mgr->InitAnalysis()) {
    return;
  }
  return;
}
#endif
