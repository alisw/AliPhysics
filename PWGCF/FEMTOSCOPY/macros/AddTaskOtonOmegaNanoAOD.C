#include <vector>
//#include "AliAnalysisTaskSE.h"
//#include "AliAnalysisManager.h"
//#include "AliAnalysisTaskOtonOmegaNanoAOD.h"
//#include "AliFemtoDreamEventCuts.h"
//#include "AliFemtoDreamTrackCuts.h"
//#include "AliFemtoDreamCascadeCuts.h"
//#include "AliFemtoDreamCollConfig.h"

#ifdef __CLING__
R__ADD_INCLUDE_PATH($ALICE_PHYSICS)
#include <PWGCF/FEMTOSCOPY/macros/ConfigOtonOmega.C>
#endif

AliAnalysisTaskSE *AddTaskOtonOmegaNanoAOD(
                                        bool GetConfigFromAlien = false,
                                        TString cFileName = "ConfigOtonOmega.C",
                                        Int_t iConfigCuts = 0
) {

  //set fullBlastQA and suffix (cut variation)
  bool fullBlastQA = true;
  const char *cutVariation = "0";
  TString suffix = TString::Format("%s", cutVariation);

  //Get Config File:
  TString configBasePath= "$ALICE_PHYSICS/PWGCF/FEMTOSCOPY/macros/";
  if(GetConfigFromAlien && (!gSystem->Exec(Form("alien_cp alien:///alice/cern.ch/user/o/ovazquez/%s .",cFileName.Data()))) ){
      configBasePath=Form("%s/",gSystem->pwd());
  } else {
      cFileName = "ConfigOtonOmega.C"; //if not getting it from alien, take the standard one
  }
  TString configFilePath(configBasePath+cFileName);
  std::cout << "Configpath: " << configFilePath << std::endl;
  TString configFunction(cFileName(0,cFileName.Length() - 2));
  if (!gROOT->GetListOfGlobalFunctions()->FindObject(configFunction.Data())) gROOT->LoadMacro(configFilePath.Data());

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
  // Track Cuts
  AliFemtoDreamTrackCuts *TrackCuts = AliFemtoDreamTrackCuts::PrimProtonCuts(false, true, false, false);
  TrackCuts->SetFilterBit(128);
  TrackCuts->SetCutCharge(1);
  // AntiTrack Cuts
  AliFemtoDreamTrackCuts *AntiTrackCuts =AliFemtoDreamTrackCuts::PrimProtonCuts(false, true, false, false);
  AntiTrackCuts->SetFilterBit(128);
  AntiTrackCuts->SetCutCharge(-1);
  //Cascade Cuts (bkg)
  AliFemtoDreamCascadeCuts* CascadeCuts = AliFemtoDreamCascadeCuts::OmegaCuts(false, false);
  //AliOtonOmegaCascadeCuts* CascadeCuts = AliOtonOmegaCascadeCuts::OmegaCuts(false, false);
  CascadeCuts->SetXiCharge(-1);
  AliFemtoDreamTrackCuts *XiNegCuts = AliFemtoDreamTrackCuts::Xiv0PionCuts(false, true, false);
  XiNegCuts->SetCheckTPCRefit(false);//for nanos this is already done while prefiltering
  AliFemtoDreamTrackCuts *XiPosCuts = AliFemtoDreamTrackCuts::Xiv0ProtonCuts(false, true, false);
  XiPosCuts->SetCheckTPCRefit(false);
  AliFemtoDreamTrackCuts *XiBachCuts = AliFemtoDreamTrackCuts::OmegaBachKaonCuts(false, true, false);
  XiBachCuts->SetCheckTPCRefit(false);
  CascadeCuts->Setv0Negcuts(XiNegCuts);
  CascadeCuts->Setv0PosCuts(XiPosCuts);
  CascadeCuts->SetBachCuts(XiBachCuts);
  CascadeCuts->SetPDGCodeCasc(3334);
  CascadeCuts->SetPDGCodev0(3122);
  CascadeCuts->SetPDGCodePosDaug(2212);
  CascadeCuts->SetPDGCodeNegDaug(-211);
  CascadeCuts->SetPDGCodeBach(-321);
  //AntiCascade cuts (bkg)
  AliFemtoDreamCascadeCuts* AntiCascadeCuts = AliFemtoDreamCascadeCuts::OmegaCuts(false, false);
  AntiCascadeCuts->SetXiCharge(1);
  AliFemtoDreamTrackCuts *AntiXiNegCuts = AliFemtoDreamTrackCuts::Xiv0ProtonCuts(false, true, false);
  AntiXiNegCuts->SetCutCharge(-1);
  AntiXiNegCuts->SetCheckTPCRefit(false);
  AliFemtoDreamTrackCuts *AntiXiPosCuts = AliFemtoDreamTrackCuts::Xiv0PionCuts(false, true, false);
  AntiXiPosCuts->SetCutCharge(1);
  AntiXiPosCuts->SetCheckTPCRefit(false);
  AliFemtoDreamTrackCuts *AntiXiBachCuts =AliFemtoDreamTrackCuts::OmegaBachKaonCuts(false, true, false);
  AntiXiBachCuts->SetCutCharge(1);
  AntiXiBachCuts->SetCheckTPCRefit(false);
  AntiCascadeCuts->Setv0Negcuts(AntiXiNegCuts);
  AntiCascadeCuts->Setv0PosCuts(AntiXiPosCuts);
  AntiCascadeCuts->SetBachCuts(AntiXiBachCuts);
  AntiCascadeCuts->SetPDGCodeCasc(-3334);
  AntiCascadeCuts->SetPDGCodev0(-3122);
  AntiCascadeCuts->SetPDGCodePosDaug(211);
  AntiCascadeCuts->SetPDGCodeNegDaug(-2212);
  AntiCascadeCuts->SetPDGCodeBach(321);
  //Omega Cuts 
  AliFemtoDreamCascadeCuts* CascadeOmegaCuts = AliFemtoDreamCascadeCuts::OmegaCuts(false, false);
  CascadeOmegaCuts->SetXiCharge(-1);
  AliFemtoDreamTrackCuts *OmegaNegCuts = AliFemtoDreamTrackCuts::Xiv0PionCuts(false, true, false);
  OmegaNegCuts->SetCheckTPCRefit(false);
  AliFemtoDreamTrackCuts *OmegaPosCuts = AliFemtoDreamTrackCuts::Xiv0ProtonCuts(false, true, false);
  OmegaPosCuts->SetCheckTPCRefit(false);
  AliFemtoDreamTrackCuts *OmegaBachCuts = AliFemtoDreamTrackCuts::OmegaBachKaonCuts(false, true, false);
  OmegaBachCuts->SetCheckTPCRefit(false);
  CascadeOmegaCuts->Setv0Negcuts(OmegaNegCuts);
  CascadeOmegaCuts->Setv0PosCuts(OmegaPosCuts);
  CascadeOmegaCuts->SetBachCuts(OmegaBachCuts);
  CascadeOmegaCuts->SetPDGCodeCasc(3334);
  CascadeOmegaCuts->SetPDGCodev0(3122);
  CascadeOmegaCuts->SetPDGCodePosDaug(2212);
  CascadeOmegaCuts->SetPDGCodeNegDaug(-211);
  CascadeOmegaCuts->SetPDGCodeBach(-321);
  //AntiCascade cuts (bkg)
  AliFemtoDreamCascadeCuts* AntiCascadeOmegaCuts = AliFemtoDreamCascadeCuts::OmegaCuts(false, false);
  AntiCascadeOmegaCuts->SetXiCharge(1);
  AliFemtoDreamTrackCuts *AntiOmegaNegCuts = AliFemtoDreamTrackCuts::Xiv0ProtonCuts(false, true, false);
  AntiOmegaNegCuts->SetCutCharge(-1);
  AntiOmegaNegCuts->SetCheckTPCRefit(false);
  AliFemtoDreamTrackCuts *AntiOmegaPosCuts = AliFemtoDreamTrackCuts::Xiv0PionCuts(false, true, false);
  AntiOmegaPosCuts->SetCutCharge(1);
  AntiOmegaPosCuts->SetCheckTPCRefit(false);
  AliFemtoDreamTrackCuts *AntiOmegaBachCuts =AliFemtoDreamTrackCuts::OmegaBachKaonCuts(false, true, false);
  AntiOmegaBachCuts->SetCutCharge(1);
  AntiOmegaBachCuts->SetCheckTPCRefit(false);
  AntiCascadeOmegaCuts->Setv0Negcuts(AntiOmegaNegCuts);
  AntiCascadeOmegaCuts->Setv0PosCuts(AntiOmegaPosCuts);
  AntiCascadeOmegaCuts->SetBachCuts(AntiOmegaBachCuts);
  AntiCascadeOmegaCuts->SetPDGCodeCasc(-3334);
  AntiCascadeOmegaCuts->SetPDGCodev0(-3122);
  AntiCascadeOmegaCuts->SetPDGCodePosDaug(211);
  AntiCascadeOmegaCuts->SetPDGCodeNegDaug(-2212);
  AntiCascadeOmegaCuts->SetPDGCodeBach(321);

  //Pass cuts to config file (to be further setup of modified there)_
  ConfigOtonOmega(CascadeCuts,XiPosCuts,XiNegCuts,XiBachCuts,AntiCascadeCuts,AntiXiPosCuts,AntiXiNegCuts,AntiXiBachCuts,CascadeOmegaCuts,OmegaPosCuts,OmegaNegCuts,OmegaBachCuts,AntiCascadeOmegaCuts,AntiOmegaPosCuts,AntiOmegaNegCuts,AntiOmegaBachCuts,TrackCuts,AntiTrackCuts,iConfigCuts);

  //???????????????????
//  if (suffix != "0" && suffix != "999") {
    evtCuts->SetMinimalBooking(true);
    TrackCuts->SetMinimalBooking(true);
    AntiTrackCuts->SetMinimalBooking(true);
    CascadeCuts->SetMinimalBooking(true);
    AntiCascadeCuts->SetMinimalBooking(true);
    CascadeOmegaCuts->SetMinimalBooking(true);
    AntiCascadeOmegaCuts->SetMinimalBooking(true);
//  }

  AliFemtoDreamCollConfig *config = new AliFemtoDreamCollConfig("Femto","Femto");

  // Femto Collection
  std::vector<int> PDGParticles;
  PDGParticles.push_back(2212);
  PDGParticles.push_back(2212);
  PDGParticles.push_back(3334);
  PDGParticles.push_back(3334);
  PDGParticles.push_back(3334);
  PDGParticles.push_back(3334);


  std::vector<int> NBins;
  std::vector<float> kMin;
  std::vector<float> kMax;
  std::vector<int> pairQA;
  std::vector<bool> closeRejection;
  //pairs:
  //p p               0
  //p barp            1
  //p Xi              2  // here and in the following Xi stands for Omega background
  //p barXi           3
  //p Omega           4
  //p barOmega        5
  //barp barp         6
  //barp Xi           7
  //barp barXi        8
  //barp Omega        9
  //barp barOmega     10
  //Xi Xi             11
  //Xi barXi          12
  //Xi Omega          13
  //Xi barOmega       14
  //barXi barXi       15
  //barXi Omega       16
  //barXi barOmega    17
  //Omega Omega       18
  //Omega barOmega    19
  //barOmega barOmega 20


//this shit is all by hand and I want to die:
/////////////////////////////////////////////
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

  pairQA[4] = 13;//p-Omega
  pairQA[10] = 13;//pbar-antiOmega

  closeRejection[0] = true;  // pp              ????
  closeRejection[6] = true;  // barp barp       ????

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

//  if (suffix != "0") {
    config->SetMinimalBookingME(true);
config->SetMinimalBookingSample(true);
//  }
  AliAnalysisTaskOtonOmegaNanoAOD* task = new AliAnalysisTaskOtonOmegaNanoAOD("OtonOmegaNanoAOD",true);
  if (suffix != "0" && suffix != "999") {
    task->SetRunTaskLightWeight(true);
  }
  task->SelectCollisionCandidates(AliVEvent::kHighMultV0);
  task->SetEventCuts(evtCuts);
  task->SetProtonCuts(TrackCuts);
  task->SetAntiProtonCuts(AntiTrackCuts);
  task->SetXiCuts(CascadeCuts);
  task->SetAntiXiCuts(AntiCascadeCuts);
  task->SetOmegaCuts(CascadeOmegaCuts);
  task->SetAntiOmegaCuts(AntiCascadeOmegaCuts);
  task->SetCorrelationConfig(config);
  mgr->AddTask(task);

  TString addon = "POmega";

  TString file = AliAnalysisManager::GetCommonFileName();

  mgr->ConnectInput(task, 0, cinput);

  if (cFileName == "ConfigOtonOmega.C"){
    std::cout << "No CONTAINERaddon from configfile " << std::endl;
  }else{
    addon += "_";
    for(int ii=15;ii<=22;ii++){ addon += cFileName(ii);};
    std::cout << "CONTAINERaddon from configfile " << cFileName(15,22).Data() << std::endl;
  }

  TString EvtCutsName = Form("%sEvtCuts%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *coutputEvtCuts = mgr->CreateContainer(
      EvtCutsName.Data(), TList::Class(), AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), EvtCutsName.Data()));
  mgr->ConnectOutput(task, 1, coutputEvtCuts);

  TString TrackCutsName = Form("%sTrackCuts%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *couputTrkCuts = mgr->CreateContainer(
      TrackCutsName.Data(), TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), TrackCutsName.Data()));
  mgr->ConnectOutput(task, 2, couputTrkCuts);

  TString AntiTrackCutsName = Form("%sAntiTrackCuts%s", addon.Data(),
                                   suffix.Data());
  AliAnalysisDataContainer *coutputAntiTrkCuts = mgr->CreateContainer(
      AntiTrackCutsName.Data(), TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), AntiTrackCutsName.Data()));
  mgr->ConnectOutput(task, 3, coutputAntiTrkCuts);

  AliAnalysisDataContainer *coutputCascadeCuts;
  TString CascadeCutsName = Form("%sCascadeCuts%s", addon.Data(), suffix.Data());
  coutputCascadeCuts = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      CascadeCutsName.Data(),
      TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), CascadeCutsName.Data()));
  mgr->ConnectOutput(task, 4, coutputCascadeCuts);


  AliAnalysisDataContainer *coutputAntiCascadeCuts;
  TString AntiCascadeCutsName = Form("%sAntiCascadeCuts%s", addon.Data(), suffix.Data());
  coutputAntiCascadeCuts = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      AntiCascadeCutsName.Data(),
      TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), AntiCascadeCutsName.Data()));
  mgr->ConnectOutput(task, 5, coutputAntiCascadeCuts);


  AliAnalysisDataContainer *coutputCascadeOmegaCuts;
  TString CascadeOmegaCutsName = Form("%sCascadeOmegaCuts%s", addon.Data(), suffix.Data());
  coutputCascadeOmegaCuts = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      CascadeOmegaCutsName.Data(),
      TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), CascadeOmegaCutsName.Data()));
  mgr->ConnectOutput(task, 6, coutputCascadeOmegaCuts);


  AliAnalysisDataContainer *coutputAntiCascadeOmegaCuts;
  TString AntiCascadeOmegaCutsName = Form("%sAntiCascadeOmegaCuts%s", addon.Data(), suffix.Data());
  coutputAntiCascadeOmegaCuts = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      AntiCascadeOmegaCutsName.Data(),
      TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), AntiCascadeOmegaCutsName.Data()));
  mgr->ConnectOutput(task, 7, coutputAntiCascadeOmegaCuts);









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


  //omega tree:
  AliAnalysisDataContainer *coutputTreeOmega;
  TString TreeOmegaName = Form("%sTreeOmega",addon.Data());
  coutputTreeOmega = mgr->CreateContainer(
    //@suppress("Invalid arguments") it works ffs
    TreeOmegaName.Data(),
    TTree::Class(),
    AliAnalysisManager::kOutputContainer,
    Form("%s:%s", file.Data(), TreeOmegaName.Data()));
  mgr->ConnectOutput(task, 10, coutputTreeOmega);



  return task;
}
