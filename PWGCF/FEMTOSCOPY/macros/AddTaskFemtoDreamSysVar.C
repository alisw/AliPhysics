#ifndef ADDTASKFEMTODREAM_C
#define ADDTASKFEMTODREAM_C
#include "AliFemtoDreamEventCuts.h"
#include "AliFemtoDreamTrackCuts.h"
#include "AliAnalysisTaskFemtoDream.h"
#include "TROOT.h"

AliAnalysisTask* AddTaskFemtoDreamSysVar(
    bool isMC=false, TString CentEst="kInt7",int sysVar=0,
    bool notpp=true, bool DCAPlots=false,bool CPAPlots=false,
    bool CombSigma=false,bool ContributionSplitting=false,
    bool ContributionSplittingDaug=false)
{
  gROOT->ProcessLine(".include $ALICE_ROOT/include");
  gROOT->ProcessLine(".include $ALICE_PHYSICS/include");
  gROOT->ProcessLine(".include $ROOTSYS/include");

  // the manager is static, so get the existing manager via the static method
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  if (!mgr)
  {
    printf("No analysis manager to connect to!\n");
    return nullptr;
  }
  // just to see if all went well, check if the input event handler has been
  // connected
  if (!mgr->GetInputEventHandler()) {
    printf("This task requires an input event handler!\n");
    return nullptr;
  }
  AliFemtoDreamEventCuts *evtCuts=
      AliFemtoDreamEventCuts::StandardCutsRun2();

  //Track Cuts
  AliFemtoDreamTrackCuts *TrackCuts=
      AliFemtoDreamTrackCuts::PrimProtonCuts(
          isMC,DCAPlots,CombSigma,ContributionSplitting);
  TrackCuts->SetCutCharge(1);
  AliFemtoDreamTrackCuts *AntiTrackCuts=
      AliFemtoDreamTrackCuts::PrimProtonCuts(
          isMC,DCAPlots,CombSigma,ContributionSplitting);
  AntiTrackCuts->SetCutCharge(-1);

  switch(sysVar){
    case 1:
      TrackCuts->SetPtRange(0.4,4.05);
      AntiTrackCuts->SetPtRange(0.4,4.05);
      break;
    case 2:
      TrackCuts->SetPtRange(0.5,4.05);
      AntiTrackCuts->SetPtRange(0.5,4.05);
      break;
    case 3:
      TrackCuts->SetEtaRange(-0.7,0.7);
      AntiTrackCuts->SetEtaRange(-0.7,0.7);
      break;
    case 3:
      TrackCuts->SetEtaRange(-0.9,0.9);
      AntiTrackCuts->SetEtaRange(-0.7,0.7);
      break;
    case 5:
      TrackCuts->SetPID(AliPID::kProton,0.75,2);
      AntiTrackCuts->SetPID(AliPID::kProton,0.75,2);
      break;
    case 6:
      TrackCuts->SetPID(AliPID::kProton,0.75,5);
      AntiTrackCuts->SetPID(AliPID::kProton,0.75,5);
      break;
    case 7:
      TrackCuts->SetFilterBit(96);
      AntiTrackCuts->SetFilterBit(96);
      break;
    case 8:
      TrackCuts->SetNClsTPC(70);
      AntiTrackCuts->SetNClsTPC(70);
      break;
    case 9:
      TrackCuts->SetNClsTPC(90);
      AntiTrackCuts->SetNClsTPC(90);
      break;
    default:
      break;
  }

  AliFemtoDreamv0Cuts *v0Cuts;
  AliFemtoDreamv0Cuts *Antiv0Cuts;
  AliFemtoDreamCascadeCuts *CascadeCuts;
  AliFemtoDreamCascadeCuts *AntiCascadeCuts;

  //Lambda Cuts
  v0Cuts=
      AliFemtoDreamv0Cuts::LambdaCuts(isMC,CPAPlots,ContributionSplitting);
  AliFemtoDreamTrackCuts *Posv0Daug=AliFemtoDreamTrackCuts::DecayProtonCuts(
      isMC,false);
  AliFemtoDreamTrackCuts *Negv0Daug=AliFemtoDreamTrackCuts::DecayPionCuts(
      isMC,false);

  Antiv0Cuts=
      AliFemtoDreamv0Cuts::LambdaCuts(isMC,CPAPlots,ContributionSplitting);
  AliFemtoDreamTrackCuts *PosAntiv0Daug=AliFemtoDreamTrackCuts::DecayPionCuts(
      isMC,false);
  PosAntiv0Daug->SetCutCharge(1);
  AliFemtoDreamTrackCuts *NegAntiv0Daug=AliFemtoDreamTrackCuts::DecayProtonCuts(
      isMC,false);
  NegAntiv0Daug->SetCutCharge(-1);


  switch(sysVar){
    case 10:
      v0Cuts->SetPtRange(0.24,999.9);
      Antiv0Cuts->SetPtRange(0.24,999.9);
      break;
    case 11:
      v0Cuts->SetPtRange(0.36,999.9);
      Antiv0Cuts->SetPtRange(0.36,999.9);
      break;
    case 12:
      v0Cuts->SetCutCPA(0.998);
      Antiv0Cuts->SetCutCPA(0.998);
      break;
    case 13:
      Posv0Daug->SetPID(AliPID::kProton,999.9,4);
      Negv0Daug->SetPID(AliPID::kPion,999.9,4);
      PosAntiv0Daug->SetPID(AliPID::kPion,999.9,4);
      NegAntiv0Daug->SetPID(AliPID::kProton,999.9,4);
      break;
    case 14:
      Posv0Daug->SetNClsTPC(80);
      Negv0Daug->SetNClsTPC(80);
      PosAntiv0Daug->SetNClsTPC(80);
      NegAntiv0Daug->SetNClsTPC(80);
      break;
    case 15:
      Posv0Daug->SetEtaRange(-0.7,0.7);
      Negv0Daug->SetEtaRange(-0.7,0.7);
      PosAntiv0Daug->SetEtaRange(-0.7,0.7);
      NegAntiv0Daug->SetEtaRange(-0.7,0.7);
      break;
    case 16:
      Posv0Daug->SetEtaRange(-0.9,0.9);
      Negv0Daug->SetEtaRange(-0.9,0.9);
      PosAntiv0Daug->SetEtaRange(-0.9,0.9);
      NegAntiv0Daug->SetEtaRange(-0.9,0.9);
      break;
    case 17:
      v0Cuts->SetCutDCADaugTov0Vtx(1.2);
      Antiv0Cuts->SetCutDCADaugTov0Vtx(1.2);
      break;
    case 18:
      v0Cuts->SetCutDCADaugToPrimVtx(0.06);
      Antiv0Cuts->SetCutDCADaugToPrimVtx(0.06);
      break;
    default:
      break;
  }

  v0Cuts->SetPosDaugterTrackCuts(Posv0Daug);
  v0Cuts->SetNegDaugterTrackCuts(Negv0Daug);
  v0Cuts->SetPDGCodePosDaug(2212);//Proton
  v0Cuts->SetPDGCodeNegDaug(211);//Pion
  v0Cuts->SetPDGCodev0(3122);//Lambda

  Antiv0Cuts->SetPosDaugterTrackCuts(PosAntiv0Daug);
  Antiv0Cuts->SetNegDaugterTrackCuts(NegAntiv0Daug);
  Antiv0Cuts->SetPDGCodePosDaug(211);//Pion
  Antiv0Cuts->SetPDGCodeNegDaug(2212);//Proton
  Antiv0Cuts->SetPDGCodev0(-3122);//Lambda

  //Cascade Cuts
  CascadeCuts=
      AliFemtoDreamCascadeCuts::XiCuts(isMC,ContributionSplitting);
  CascadeCuts->SetXiCharge(-1);
  AliFemtoDreamTrackCuts *XiNegCuts=
      AliFemtoDreamTrackCuts::Xiv0PionCuts(isMC,false);
  AliFemtoDreamTrackCuts *XiPosCuts=
      AliFemtoDreamTrackCuts::Xiv0ProtonCuts(isMC,false);
  AliFemtoDreamTrackCuts *XiBachCuts=
      AliFemtoDreamTrackCuts::XiBachPionCuts(isMC,false);

  AntiCascadeCuts=
      AliFemtoDreamCascadeCuts::XiCuts(isMC,ContributionSplitting);
  AntiCascadeCuts->SetXiCharge(1);
  AliFemtoDreamTrackCuts *AntiXiNegCuts=
      AliFemtoDreamTrackCuts::Xiv0ProtonCuts(isMC,false);
  AntiXiNegCuts->SetCutCharge(-1);
  AliFemtoDreamTrackCuts *AntiXiPosCuts=
      AliFemtoDreamTrackCuts::Xiv0PionCuts(isMC,false);
  AntiXiPosCuts->SetCutCharge(1);
  AliFemtoDreamTrackCuts *AntiXiBachCuts=
      AliFemtoDreamTrackCuts::XiBachPionCuts(isMC,false);
  AntiXiBachCuts->SetCutCharge(1);

  switch(sysVar){
    case 19:
      CascadeCuts->SetCutXiDaughterDCA(1.3);
      AntiCascadeCuts->SetCutXiDaughterDCA(1.3);
      break;
    case 20:
      CascadeCuts->SetCutXiMinDistBachToPrimVtx(0.06);
      AntiCascadeCuts->SetCutXiMinDistBachToPrimVtx(0.06);
      break;
    case 21:
      CascadeCuts->SetCutXiCPA(0.99);
      AntiCascadeCuts->SetCutXiCPA(0.99);
      break;
    case 22:
      CascadeCuts->SetCutXiTransverseRadius(1.0,200);
      AntiCascadeCuts->SetCutXiTransverseRadius(1.0,200);
      break;
    case 23:
      CascadeCuts->SetCutv0MaxDaughterDCA(1.3);
      AntiCascadeCuts->SetCutv0MaxDaughterDCA(1.3);
      break;
    case 24:
      CascadeCuts->SetCutv0CPA(0.99);
      AntiCascadeCuts->SetCutv0CPA(0.99);
      break;
    case 25:
      CascadeCuts->SetCutv0TransverseRadius(1.7,200);
      AntiCascadeCuts->SetCutv0TransverseRadius(1.7,200);
      break;
    case 26:
      CascadeCuts->SetCutv0MinDistToPrimVtx(0.08);
      AntiCascadeCuts->SetCutv0MinDistToPrimVtx(0.08);
      break;
    case 27:
      CascadeCuts->SetCutv0MinDaugDistToPrimVtx(0.05);
      AntiCascadeCuts->SetCutv0MinDaugDistToPrimVtx(0.05);
      break;
    case 28:
      XiNegCuts->SetEtaRange(-0.7,0.7);
      XiPosCuts->SetEtaRange(-0.7,0.7);
      XiBachCuts->SetEtaRange(-0.7,0.7);
      AntiXiNegCuts->SetEtaRange(-0.7,0.7);
      AntiXiPosCuts->SetEtaRange(-0.7,0.7);
      AntiXiBachCuts->SetEtaRange(-0.7,0.7);
      break;
    case 29:
      XiNegCuts->SetEtaRange(-0.9,0.9);
      XiPosCuts->SetEtaRange(-0.9,0.9);
      XiBachCuts->SetEtaRange(-0.9,0.9);
      AntiXiNegCuts->SetEtaRange(-0.9,0.9);
      AntiXiPosCuts->SetEtaRange(-0.9,0.9);
      AntiXiBachCuts->SetEtaRange(-0.9,0.9);
      break;
    case 30:
      XiNegCuts->SetPID(AliPID::kPion,999,3);
      XiPosCuts->SetPID(AliPID::kProton,999,3);
      XiBachCuts->SetPID(AliPID::kPion,999,3);
      AntiXiNegCuts->SetPID(AliPID::kProton,999,3);
      AntiXiPosCuts->SetPID(AliPID::kPion,999,3);
      AntiXiBachCuts->SetPID(AliPID::kPion,999,3);
      break;
    default:
      break;
  }
  CascadeCuts->Setv0Negcuts(XiNegCuts);
  CascadeCuts->Setv0PosCuts(XiPosCuts);
  CascadeCuts->SetBachCuts(XiBachCuts);
  CascadeCuts->SetPDGCodeCasc(-3312);
  CascadeCuts->SetPDGCodev0(3122);
  CascadeCuts->SetPDGCodePosDaug(2212);
  CascadeCuts->SetPDGCodeNegDaug(-211);
  CascadeCuts->SetPDGCodeBach(-211);

  AntiCascadeCuts->Setv0Negcuts(AntiXiNegCuts);
  AntiCascadeCuts->Setv0PosCuts(AntiXiPosCuts);
  AntiCascadeCuts->SetBachCuts(AntiXiBachCuts);
  AntiCascadeCuts->SetPDGCodeCasc(3312);
  AntiCascadeCuts->SetPDGCodev0(-3122);
  AntiCascadeCuts->SetPDGCodePosDaug(211);
  AntiCascadeCuts->SetPDGCodeNegDaug(-2212);
  AntiCascadeCuts->SetPDGCodeBach(-211);

  std::vector<int> PDGParticles ={2212,2212,3122,3122,3312,3312};
  std::vector<double> ZVtxBins = {-10,-8,-6,-4,-2,0,2,4,6,8,10};
  std::vector<int> NBins= {750,750,150,150,150,150,750,150,150,150,150,150,150,150,150,150,150,150,150,150,150};
  std::vector<double> kMin= {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  std::vector<double> kMax= {3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3.};
  AliFemtoDreamCollConfig *config=new AliFemtoDreamCollConfig("Femto","Femto");
  if (notpp) {
    std::vector<int> MultBins = {0,4,8,12,16,20,24,28,32,36,40,44,48,52,56,60,64,68,72,80};
    config->SetMultBins(MultBins);
  } else {
    std::vector<int> MultBins = {0,4,8,12,16,20,24,28,32,36,40,60,80};
    config->SetMultBins(MultBins);
  }
  config->SetMultBinning(true);
  config->SetZBins(ZVtxBins);
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

    AliAnalysisDataContainer *coutputXiCutsMC;
    TString XiCutsMCName = Form("XiCutsMC");
    coutputXiCutsMC = mgr->CreateContainer(//@suppress("Invalid arguments") it works ffs
        XiCutsMCName.Data(), TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), XiCutsMCName.Data()));
    mgr->ConnectOutput(task, 15, coutputXiCutsMC);

    AliAnalysisDataContainer *coutputAntiXiCutsMC;
    TString AntiXiCutsMCName = Form("AntiXiCutsMC");
    coutputAntiXiCutsMC = mgr->CreateContainer(//@suppress("Invalid arguments") it works ffs
        AntiXiCutsMCName.Data(), TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), AntiXiCutsMCName.Data()));
    mgr->ConnectOutput(task, 16, coutputAntiXiCutsMC);
  }
  if (!mgr->InitAnalysis()) {
    return nullptr;
  }
  return task;
}
#endif

