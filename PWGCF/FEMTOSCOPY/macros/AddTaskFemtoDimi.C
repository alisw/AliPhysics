#include "TROOT.h"
#include "TSystem.h"

AliAnalysisTaskSE* AddTaskFemtoDimi(
    bool isMC = false,                //1
    TString CentEst = "kInt7",        //2
    bool MomReso = false,             //3 to set to true only when running on MC
    bool etaPhiPlotsAtTPCRadii=false,  //4 to set to true only when running on MC but very Mem. Consuming
    bool dPhidEtaPlots=true,          //5
    const char *mixmethod = "0") {    //6

  TString suffix=Form("%s",mixmethod);

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  if (!mgr) {
    printf("No analysis manager to connect to!\n");
    return nullptr;
  }
  // just to see if all went well, check if the input event handler has been
  // connected
  if (!mgr->GetInputEventHandler()) {
    printf("This task requires an input event handler!\n");
    return nullptr;
  }

  if(!(AliPIDResponse*)mgr->GetTask("PIDResponseTask")){
    if (isMC) {
      // IMPORTANT - SET WHEN USING DIFFERENT PASS
      AliAnalysisTaskPIDResponse *pidResponse =
          reinterpret_cast<AliAnalysisTaskPIDResponse *>(
              gInterpreter->ExecuteMacro("$AL5ICE_ROOT/ANALYSIS/macros/"
                  "AddTaskPIDResponse.C (kTRUE, kTRUE, "
                  "kTRUE, \"1\")"));
    } else {
      AliAnalysisTaskPIDResponse *pidResponse =
          reinterpret_cast<AliAnalysisTaskPIDResponse *>(
              gInterpreter->ExecuteMacro(
                  "$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C)"));
    }
  }

  AliFemtoDreamEventCuts *evtCuts = AliFemtoDreamEventCuts::StandardCutsRun2();
  evtCuts->CleanUpMult(false, false, false, true);
  evtCuts->SetMultVsCentPlots(true);
  evtCuts->SetDoSphericityCuts(false);
  evtCuts->SetDoSpherocityCuts(false);


  AliAnalysisTaskGrandma *task = new AliAnalysisTaskGrandma("myFirstTask",
                                                            isMC);

//Track cuts
  AliFemtoDreamTrackCuts *TrackCuts = AliFemtoDreamTrackCuts::PrimProtonCuts(
      isMC, true, false, true);
  TrackCuts->SetCutCharge(1);
  AliFemtoDreamTrackCuts *AntiTrackCuts =
      AliFemtoDreamTrackCuts::PrimProtonCuts(isMC, true, false, true);
  AntiTrackCuts->SetCutCharge(-1);

  //V0 cuts
  AliFemtoDreamv0Cuts *v0Cuts = AliFemtoDreamv0Cuts::LambdaCuts(
      isMC,true,true);
  AliFemtoDreamTrackCuts *Posv0Daug=AliFemtoDreamTrackCuts::DecayProtonCuts(
          isMC,true,false);
  AliFemtoDreamTrackCuts *Negv0Daug=AliFemtoDreamTrackCuts::DecayPionCuts(
          isMC,true,false);
  v0Cuts->SetPosDaugterTrackCuts(Posv0Daug);
  v0Cuts->SetNegDaugterTrackCuts(Negv0Daug);
  v0Cuts->SetPDGCodePosDaug(2212);//Proton
  v0Cuts->SetPDGCodeNegDaug(211);//Pion
  v0Cuts->SetPDGCodev0(3122);//Lambda

  AliFemtoDreamv0Cuts *Antiv0Cuts = AliFemtoDreamv0Cuts::LambdaCuts(
      isMC, true,true);
  AliFemtoDreamTrackCuts *PosAntiv0Daug=AliFemtoDreamTrackCuts::DecayPionCuts(
          isMC,true,false);
  PosAntiv0Daug->SetCutCharge(1);
  AliFemtoDreamTrackCuts *NegAntiv0Daug=AliFemtoDreamTrackCuts::DecayProtonCuts(
          isMC,true,false);
  NegAntiv0Daug->SetCutCharge(-1);
  Antiv0Cuts->SetPosDaugterTrackCuts(PosAntiv0Daug);
  Antiv0Cuts->SetNegDaugterTrackCuts(NegAntiv0Daug);
  Antiv0Cuts->SetPDGCodePosDaug(211);//Pion
  Antiv0Cuts->SetPDGCodeNegDaug(2212);//Proton
  Antiv0Cuts->SetPDGCodev0(-3122);//Lambda

      //Cascade Cuts
  AliFemtoDreamCascadeCuts* CascadeCuts = AliFemtoDreamCascadeCuts::XiCuts(
      isMC, false);
  CascadeCuts->SetXiCharge(-1);
    AliFemtoDreamTrackCuts *XiNegCuts = AliFemtoDreamTrackCuts::Xiv0PionCuts(
      isMC, true, false);
  XiNegCuts->SetCheckTPCRefit(false);  //for nanos this is already done while prefiltering
  AliFemtoDreamTrackCuts *XiPosCuts = AliFemtoDreamTrackCuts::Xiv0ProtonCuts(
      isMC, true, false);
  XiPosCuts->SetCheckTPCRefit(false);  //for nanos this is already done while prefiltering
  AliFemtoDreamTrackCuts *XiBachCuts = AliFemtoDreamTrackCuts::XiBachPionCuts(
      isMC, true, false);
  XiBachCuts->SetCheckTPCRefit(false);  //for nanos this is already done while prefiltering

  CascadeCuts->Setv0Negcuts(XiNegCuts);
  CascadeCuts->Setv0PosCuts(XiPosCuts);
  CascadeCuts->SetBachCuts(XiBachCuts);
  CascadeCuts->SetPDGCodeCasc(3312);
  CascadeCuts->SetPDGCodev0(3122);
  CascadeCuts->SetPDGCodePosDaug(2212);
  CascadeCuts->SetPDGCodeNegDaug(-211);
  CascadeCuts->SetPDGCodeBach(-211);

  AliFemtoDreamCascadeCuts* AntiCascadeCuts = AliFemtoDreamCascadeCuts::XiCuts(
      isMC, false);
  AntiCascadeCuts->SetXiCharge(1);
  AliFemtoDreamTrackCuts *AntiXiNegCuts =
      AliFemtoDreamTrackCuts::Xiv0ProtonCuts(isMC, true, false);
  AntiXiNegCuts->SetCutCharge(-1);
  AntiXiNegCuts->SetCheckTPCRefit(true);  //for nanos this is already done while prefiltering
  AliFemtoDreamTrackCuts *AntiXiPosCuts = AliFemtoDreamTrackCuts::Xiv0PionCuts(
      isMC, true, false);
  AntiXiPosCuts->SetCutCharge(1);
  AntiXiPosCuts->SetCheckTPCRefit(true);  //for nanos this is already done while prefiltering
  AliFemtoDreamTrackCuts *AntiXiBachCuts =
      AliFemtoDreamTrackCuts::XiBachPionCuts(isMC, true, false);
  AntiXiBachCuts->SetCutCharge(1);
  AntiXiBachCuts->SetCheckTPCRefit(true);  //for nanos this is already done while prefiltering

  AntiCascadeCuts->Setv0Negcuts(AntiXiNegCuts);
  AntiCascadeCuts->Setv0PosCuts(AntiXiPosCuts);
  AntiCascadeCuts->SetBachCuts(AntiXiBachCuts);
  AntiCascadeCuts->SetPDGCodeCasc(-3312);
  AntiCascadeCuts->SetPDGCodev0(-3122);
  AntiCascadeCuts->SetPDGCodePosDaug(211);
  AntiCascadeCuts->SetPDGCodeNegDaug(-2212);
  AntiCascadeCuts->SetPDGCodeBach(211);


  AliFemtoDreamCollConfig *config = new AliFemtoDreamCollConfig("Femto",
                                                                "Femto");

  std::vector<int> pairQA;
  std::vector<int> NBins;
  std::vector<float> kMin;
  std::vector<float> kMax;
  std::vector<bool> closeRejection;

  std::vector<int> PDGParticles;
  PDGParticles.push_back(2212);//proton
  PDGParticles.push_back(2212);
  PDGParticles.push_back(3122);//Lambda
  PDGParticles.push_back(3122);
  PDGParticles.push_back(3312);//Cascade
  PDGParticles.push_back(3312);

  //pairs:
  //pp                0
  //p bar p           1
  //p La              2
  //p bar La          3
  //bar p bar p       4
  //bar p La          5
  //bar p bar La      6
  //p Xi              7
  //p bar Xi          8
  //bar p Xi          9
  //bar p bar Xi      10
  //La La             11
  //La bar La         12
  //bar La bar La     13
  //La Xi             14
  //bar La Xi         15
  //La bar Xi         16
  //bar La bar Xi     17
  //Xi Xi             18
  //Xi bar Xi         19
  //Xi bar Xi bar     20

  const int nPairs = 21;
  for (int i = 0; i < nPairs; ++i) {
    pairQA.push_back(0);
    closeRejection.push_back(false);
    NBins.push_back(1500);
    kMin.push_back(0.);
    kMax.push_back(6.);
  }
  pairQA[0] = 11;
  pairQA[1] = 11;
  pairQA[2] = 12;
  pairQA[3] = 12;
  pairQA[4] = 11;
  pairQA[5] = 12;
  pairQA[6] = 12;
  pairQA[7] = 13;
  pairQA[8] = 13;
  pairQA[9] = 13;
  pairQA[10] = 13;
  pairQA[11] = 22;
  pairQA[12] = 22;
  pairQA[13] = 22;
  pairQA[14] = 23;
  pairQA[15] = 23;
  pairQA[16] = 23;
  pairQA[17] = 23;
  pairQA[18] = 33;
  pairQA[19] = 33;
  pairQA[20] = 33;

  closeRejection[0] = true;  // pp
  closeRejection[4] = true;  // barp barp

  closeRejection[18] = true;  // Xi Xi
  closeRejection[20] = true;  // barXi barXi

  config->SetPDGCodes(PDGParticles);
  config->SetNBinsHist(NBins);
  config->SetMinKRel(kMin);
  config->SetMaxKRel(kMax);
  config->SetClosePairRejection(closeRejection);
  config->SetDeltaEtaMax(0.012);
  config->SetDeltaPhiMax(0.012);
  config->SetExtendedQAPairs(pairQA);

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

  std::vector<int> centBins;
  centBins.push_back(20);
  centBins.push_back(40);
  centBins.push_back(90);
  config->SetCentBins(centBins);
  config->SetkTCentralityBinning(false);


if(isMC)
  {
  config->SetdPhidEtaPlots(dPhidEtaPlots);  // warsaw like plots
} else {
  config->SetdPhidEtaPlots(dPhidEtaPlots);  // warsaw like plots
}

if (MomReso) {
  if (isMC) {
    config->SetMomentumResolution(true);//kstar true vs. kstar reco
  } else {
    std::cout
        << "You are trying to request the Momentum Resolution without MC Info; fix it wont work! \n";
  }
}


  config->SetNBinsHist(NBins);
  config->SetMinKRel(kMin);
  config->SetMaxKRel(kMax);
  config->SetMultBinning(true);
  config->SetCentBinning(false);
  config->SetkTBinning(true);
  config->SetmTBinning(true);

  config->SetMinimalBookingME(false);
  config->SetMinimalBookingSample(true);

  //Setting the configurations of the mixing methods to run on trains:

  if(suffix == "0"){
    config->SetMixingDepth(10);
    config->SetUseEventMixing(true);
  } else if (suffix == "1a"){
    config->SetUseEventMixing(false);
    config->SetUsePhiSpinning(true);
    config->SetControlMethod(AliFemtoDreamCollConfig::kCorrelatedPhi);
    config->SetCorrelationRange(0.1);
    config->SetSpinningDepth(1);
  } else if (suffix == "1b"){
    config->SetUseEventMixing(false);
    config->SetUsePhiSpinning(true);
    config->SetControlMethod(AliFemtoDreamCollConfig::kCorrelatedPhi);
    config->SetCorrelationRange(0.2);
    config->SetSpinningDepth(1);
  } else if (suffix == "1c"){
    config->SetUseEventMixing(false);
    config->SetUsePhiSpinning(true);
    config->SetControlMethod(AliFemtoDreamCollConfig::kCorrelatedPhi);
    config->SetCorrelationRange(0.3);
    config->SetSpinningDepth(1);
  } else if (suffix == "2"){
    config->SetUseEventMixing(false);
    config->SetUsePhiSpinning(true);
    config->SetControlMethod(AliFemtoDreamCollConfig::kStravinsky);
    config->SetSpinningDepth(1);
  } else if (suffix == "3"){
    config->SetUseEventMixing(false);
    config->SetUsePhiSpinning(true);
    config->SetControlMethod(AliFemtoDreamCollConfig::kPhiSpin);
    config->SetSpinningDepth(1);
  }

  config->SetMultiplicityEstimator(AliFemtoDreamEvent::kRef08);


  if (CentEst == "kInt7") {
    task->SelectCollisionCandidates(AliVEvent::kINT7);
    std::cout << "Added kINT7 Trigger \n";
  } else if (CentEst == "kMB") {
    task->SelectCollisionCandidates(AliVEvent::kMB);
    std::cout << "Added kMB Trigger \n";
  } else if (CentEst == "kHM") {
    task->SelectCollisionCandidates(AliVEvent::kHighMultV0);
    std::cout << "Added kHighMultV0 Trigger \n";
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




  task->SetEvtCutQA(true);
  task->SetTrackBufferSize(2000);
  task->SetEventCuts(evtCuts);
  task->SetTrackCuts(TrackCuts);
  task->SetAntiTrackCuts(AntiTrackCuts);
  task->Setv0Cuts(v0Cuts);
  task->SetAntiv0Cuts(Antiv0Cuts);
  task->SetXiCuts(CascadeCuts);
  task->SetAntiXiCuts(AntiCascadeCuts);
  task->SetCollectionConfig(config);
  mgr->AddTask(task);

  TString file = AliAnalysisManager::GetCommonFileName();

  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();

  mgr->ConnectInput(task, 0, cinput);
  TString addon="";
  if (CentEst=="kInt7") {
    addon+="MBDimi";
  } else if (CentEst=="kHM") {
    addon+="HMDimi";
  }

  AliAnalysisDataContainer *coutputQA;

  std::cout << "CONTAINTER NAME: " << addon.Data() << " " << suffix.Data() << std::endl;
  TString QAName;
   QAName = Form("%sQA%s",addon.Data(),suffix.Data());
  coutputQA = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      QAName.Data(),
      TList::Class(), AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), QAName.Data()));
  mgr->ConnectOutput(task, 1, coutputQA);

  AliAnalysisDataContainer *coutputEvtCuts;
  TString EvtCutsName;
   EvtCutsName = Form("%sEvtCuts%s",addon.Data(),suffix.Data());
  coutputEvtCuts = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      EvtCutsName.Data(),
      TList::Class(), AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), EvtCutsName.Data()));
  mgr->ConnectOutput(task, 2, coutputEvtCuts);

  AliAnalysisDataContainer *couputTrkCuts;
  TString TrackCutsName;
 TrackCutsName = Form("%sTrackCuts%s",addon.Data(),suffix.Data());
  couputTrkCuts = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      TrackCutsName.Data(),
      TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), TrackCutsName.Data()));
  mgr->ConnectOutput(task, 3, couputTrkCuts);

  AliAnalysisDataContainer *coutputAntiTrkCuts;
  TString AntiTrackCutsName;
  AntiTrackCutsName = Form("%sAntiTrackCuts%s",addon.Data(),suffix.Data());
  coutputAntiTrkCuts = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      AntiTrackCutsName.Data(),
      TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), AntiTrackCutsName.Data()));
  mgr->ConnectOutput(task, 4, coutputAntiTrkCuts);

  AliAnalysisDataContainer *couputv0Cuts;
  TString v0CutsName;
   v0CutsName = Form("%sv0Cuts%s",addon.Data(),suffix.Data());
  couputv0Cuts = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      v0CutsName.Data(),
      TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), v0CutsName.Data()));
  mgr->ConnectOutput(task, 5, couputv0Cuts);

  AliAnalysisDataContainer *coutputAntiv0Cuts;
  TString Antiv0CutsName;
   Antiv0CutsName = Form("%sAntiv0Cuts%s",addon.Data(),suffix.Data());
  coutputAntiv0Cuts = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      Antiv0CutsName.Data(),
      TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), Antiv0CutsName.Data()));
  mgr->ConnectOutput(task, 6, coutputAntiv0Cuts);

    AliAnalysisDataContainer *couputXiCuts;
  TString XiCutsName;
   XiCutsName = Form("%sXiCuts%s",addon.Data(),suffix.Data());
  couputXiCuts = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      XiCutsName.Data(),
      TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), XiCutsName.Data()));
  mgr->ConnectOutput(task, 7, couputXiCuts);

  AliAnalysisDataContainer *coutputAntiXiCuts;
  TString AntiXiCutsName;
   AntiXiCutsName = Form("%sAntiXiCuts%s",addon.Data(),suffix.Data());
  coutputAntiXiCuts = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      AntiXiCutsName.Data(),
      TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), AntiXiCutsName.Data()));
  mgr->ConnectOutput(task, 8, coutputAntiXiCuts);

  AliAnalysisDataContainer *coutputResults;
  TString ResultsName;
   ResultsName = Form("%sResults%s",addon.Data(),suffix.Data());
  coutputResults = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      ResultsName.Data(),
      TList::Class(), AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), ResultsName.Data()));
  mgr->ConnectOutput(task, 9, coutputResults);

  AliAnalysisDataContainer *coutputResultQA;
  TString ResultQAName;
  ResultQAName = Form("%sResultQA%s",addon.Data(),suffix.Data());
  coutputResultQA = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      ResultQAName.Data(),
      TList::Class(), AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), ResultQAName.Data()));
  mgr->ConnectOutput(task, 10, coutputResultQA);

    AliAnalysisDataContainer *coutputResultsSample;
  TString ResultsSampleName = Form("%sResultsSample%s", addon.Data(),
                                   suffix.Data());
  coutputResultsSample = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      ResultsSampleName.Data(),
      TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), ResultsSampleName.Data()));
  mgr->ConnectOutput(task, 11, coutputResultsSample);

  AliAnalysisDataContainer *coutputResultsSampleQA;
  TString ResultsSampleQAName = Form("%sResultsSampleQA%s", addon.Data(),
                                     suffix.Data());
  coutputResultsSampleQA = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      ResultsSampleQAName.Data(),
      TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), ResultsSampleQAName.Data()));
  mgr->ConnectOutput(task, 12, coutputResultsSampleQA);

  if (isMC) {
    AliAnalysisDataContainer *coutputTrkCutsMC;
    TString TrkCutsMCName = Form("%sTrkCutsMC%s",addon.Data(),suffix.Data());
    coutputTrkCutsMC = mgr->CreateContainer(
        //@suppress("Invalid arguments") it works ffs
        TrkCutsMCName.Data(),
        TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), TrkCutsMCName.Data()));
    mgr->ConnectOutput(task, 13, coutputTrkCutsMC);

    AliAnalysisDataContainer *coutputv0CutsMC;
    TString v0CutsMCName = Form("%sv0CutsMC%s",addon.Data(),suffix.Data());
    coutputv0CutsMC = mgr->CreateContainer(
        //@suppress("Invalid arguments") it works ffs
        v0CutsMCName.Data(),
        TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), v0CutsMCName.Data()));
    mgr->ConnectOutput(task, 14, coutputv0CutsMC);

    AliAnalysisDataContainer *coutputAntiTrkCutsMC;
    TString AntiTrkCutsMCName = Form("%sAntiTrkCutsMC%s",addon.Data(),suffix.Data());
    coutputAntiTrkCutsMC = mgr->CreateContainer(
        //@suppress("Invalid arguments") it works ffs
        AntiTrkCutsMCName.Data(),
        TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), AntiTrkCutsMCName.Data()));
    mgr->ConnectOutput(task, 15, coutputAntiTrkCutsMC);

    AliAnalysisDataContainer *coutputAntiv0CutsMC;
    TString Antiv0CutsMCName = Form("%sAntiv0CutsMC%s",addon.Data(),suffix.Data());
    coutputAntiv0CutsMC = mgr->CreateContainer(
        //@suppress("Invalid arguments") it works ffs
        Antiv0CutsMCName.Data(),
        TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), Antiv0CutsMCName.Data()));
    mgr->ConnectOutput(task, 16, coutputAntiv0CutsMC);

    AliAnalysisDataContainer *coutputXiCutsMC;
    TString XiCutsMCName = Form("%sXiCutsMC%s",addon.Data(),suffix.Data());
    coutputXiCutsMC = mgr->CreateContainer(
        //@suppress("Invalid arguments") it works ffs
        XiCutsMCName.Data(),
        TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), XiCutsMCName.Data()));
    mgr->ConnectOutput(task, 17, coutputXiCutsMC);

     AliAnalysisDataContainer *coutputAntiXiCutsMC;
    TString AntiXiCutsMCName = Form("%sAntiXiCutsMC%s",addon.Data(),suffix.Data());
    coutputAntiXiCutsMC = mgr->CreateContainer(
        //@suppress("Invalid arguments") it works ffs
        AntiXiCutsMCName.Data(),
        TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), AntiXiCutsMCName.Data()));
    mgr->ConnectOutput(task, 18, coutputAntiXiCutsMC);
  }

  return task;

}
