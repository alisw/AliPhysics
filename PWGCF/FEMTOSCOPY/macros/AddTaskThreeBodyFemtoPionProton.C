#if !defined(__CINT__) || defined(__CLING__)
#include <vector>
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskThreeBodyFemtoAODPionProton.h"
#include "AliAnalysisTaskThreeBodyFemtoPionProton.h"
#include "AliFemtoDreamEventCuts.h"
#include "AliFemtoDreamTrackCuts.h"
#include "AliFemtoDreamCascadeCuts.h"
#include "AliFemtoDreamCollConfig.h"
#endif

AliAnalysisTaskSE *AddTaskThreeBodyFemtoPionProton(int trigger = 0, bool fullBlastQA = true,
                                    bool isMC = false, bool isNano = true, bool triggerOn = false, bool sphericityON = false,
                                    const char *cutVariation = "0", const char *triggerVariation = "0") {



  TString suffix = TString::Format("%s", cutVariation);
  TString suffixTrigger = TString::Format("%s", triggerVariation);


  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskThreeBodyFemtoPionProton()", "No analysis manager found.");
    return 0x0;
  }

  // ================== GetInputEventHandler =============================
  AliVEventHandler *inputHandler = mgr->GetInputEventHandler();
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();

  //========= Init subtasks and start analyis ============================
  // Event Cuts
  AliFemtoDreamEventCuts *evtCuts = AliFemtoDreamEventCuts::StandardCutsRun2();
  evtCuts->CleanUpMult(false, false, false, true);


  if (sphericityON==true){
    float SpherDown = 0.7;
    evtCuts->SetSphericityCuts(SpherDown, 1.0, 0.5); // THINK IF NEEDED FOR THREE BODY
  }
  
  


  // Track Cuts
  AliFemtoDreamTrackCuts *TrackCuts = AliFemtoDreamTrackCuts::PrimProtonCuts(
      isMC, true, false, false);
  TrackCuts->SetFilterBit(128);
  TrackCuts->SetCutCharge(1);

  AliFemtoDreamTrackCuts *AntiTrackCuts =
      AliFemtoDreamTrackCuts::PrimProtonCuts(isMC, true, false, false);
  AntiTrackCuts->SetFilterBit(128);
  AntiTrackCuts->SetCutCharge(-1);

  // Track Cuts for pions ____________________________________________________________________________________

  AliFemtoDreamTrackCuts *fTrackCutsPosPion=new AliFemtoDreamTrackCuts();
  fTrackCutsPosPion->SetIsMonteCarlo(isMC);
  fTrackCutsPosPion->SetCutCharge(1);
  fTrackCutsPosPion->SetPtRange(0.14, 4.0);
  fTrackCutsPosPion->SetEtaRange(-0.8, 0.8);
  fTrackCutsPosPion->SetNClsTPC(80);
  // Not mention in AN oder Indico
  fTrackCutsPosPion->SetDCAReCalculation(true);//Get the dca from the PropagateToVetex
  fTrackCutsPosPion->SetFilterBit(128);//96); // Filterbit 5+6
  fTrackCutsPosPion->SetDCAVtxZ(0.3);
  fTrackCutsPosPion->SetDCAVtxXY(0.3);
  // Cut on avrg. separation in TPC: <Dr> < 12 cm (10 cm, 3 cm); Share quality < 1.0; share fraction < 0.05

  // FOR NOW OF BEACUSE OF MAX's INFORMATION
  //fTrackCutsPosPion->SetCutSharedCls(true);

  fTrackCutsPosPion->SetNClsTPC(80); // In Indico + additional Chi²/NDF <4
  fTrackCutsPosPion->SetPID(AliPID::kPion, 0.5);
  fTrackCutsPosPion->SetRejLowPtPionsTOF(false);
  fTrackCutsPosPion->SetMinimalBooking(false);
  //this checks if the sigma of the wanted hypothesis is the smallest, and if
  //another particle has a smaller sigma, the track is rejected.
  // Not mention in AN oder Indico
  //fTrackCutsPosPion->SetCutSmallestSig(true);
  fTrackCutsPosPion->SetPlotDCADist(true);

 //The same things for negative pions
  AliFemtoDreamTrackCuts *fTrackCutsNegPion=new AliFemtoDreamTrackCuts();
  fTrackCutsNegPion->SetIsMonteCarlo(isMC);
  fTrackCutsNegPion->SetCutCharge(-1);
  fTrackCutsNegPion->SetPtRange(0.14, 4.0);
  fTrackCutsNegPion->SetEtaRange(-0.8, 0.8);
  fTrackCutsNegPion->SetNClsTPC(80);
  fTrackCutsNegPion->SetDCAReCalculation(true);

  // FOR NOW OF BEACUSE OF MAX's INFORMATION
  //fTrackCutsNegPion->SetCutSharedCls(true);}

  fTrackCutsNegPion->SetNClsTPC(80);
  fTrackCutsNegPion->SetPID(AliPID::kPion, 0.5);
  fTrackCutsNegPion->SetRejLowPtPionsTOF(false);
  fTrackCutsNegPion->SetMinimalBooking(false);
  //fTrackCutsNegPion->SetCutSmallestSig(true);
  fTrackCutsNegPion->SetPlotDCADist(true);

  fTrackCutsNegPion->SetFilterBit(128);//96);
  fTrackCutsNegPion->SetDCAVtxZ(0.3);
  fTrackCutsNegPion->SetDCAVtxXY(0.3);

  //fTrackCutsNegPion->SetCutSharedCls(true);


  if (!fullBlastQA) {
    evtCuts->SetMinimalBooking(true);
    TrackCuts->SetMinimalBooking(true);
    AntiTrackCuts->SetMinimalBooking(true);
    fTrackCutsPosPion->SetMinimalBooking(true);
    fTrackCutsNegPion->SetMinimalBooking(true);
  }


  float Q3Limit=0.;

  if(triggerOn){
    if(suffixTrigger=="0"){
      Q3Limit = 1.0;
    }
  }
  AliFemtoDreamCollConfig *config = new AliFemtoDreamCollConfig("Femto",
                                                                "Femto", false);
  // Femto Collection
  std::vector<int> PDGParticles;
  PDGParticles.push_back(2212); //proton
  PDGParticles.push_back(2212); //antiproton
  PDGParticles.push_back(211); // pi+
  PDGParticles.push_back(211); // pi-

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

  pairQA[7] = 11;
  pairQA[8] = 11;
  pairQA[9] = 11;

  closeRejection[7] = true;  
  closeRejection[8] = true;  
  closeRejection[9] = true;  

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

  TString addon = "PPION";
  if(sphericityON) addon="PPIONSphericity";
  if(!sphericityON) addon="PPIONnoSphericity";
  TString file = AliAnalysisManager::GetCommonFileName();

  TString EvtCutsName = Form("%sEvtCuts_%s_%s", addon.Data(), suffixTrigger.Data(), suffix.Data());
  AliAnalysisDataContainer *coutputEvtCuts = mgr->CreateContainer(
      EvtCutsName.Data(), TList::Class(), AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), EvtCutsName.Data()));

  TString TrackCutsName = Form("%sTrackCutsProton_%s_%s", addon.Data(), suffixTrigger.Data(), suffix.Data());
  AliAnalysisDataContainer *couputTrkCuts = mgr->CreateContainer(
      TrackCutsName.Data(), TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), TrackCutsName.Data()));

  TString AntiTrackCutsName = Form("%sAntiTrackCutsProton_%s_%s", addon.Data(), suffixTrigger.Data(), suffix.Data());
  AliAnalysisDataContainer *coutputAntiTrkCuts = mgr->CreateContainer(
      AntiTrackCutsName.Data(), TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), AntiTrackCutsName.Data()));

  AliAnalysisDataContainer *couputTrkCutsPion;
  TString TrackCutsNamePion = Form("%sTrackCutsPion_%s_%s", addon.Data(), suffixTrigger.Data(), suffix.Data());
  couputTrkCutsPion = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      TrackCutsNamePion.Data(),
      TList::Class(), AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), TrackCutsNamePion.Data()));

  AliAnalysisDataContainer *couputAntiTrkCutsPion;
  TString AntiTrackCutsNamePion = Form("%sAntiTrackCutsPion_%s_%s", addon.Data(), suffixTrigger.Data(), suffix.Data());
  couputAntiTrkCutsPion = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      AntiTrackCutsNamePion.Data(),
      TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), AntiTrackCutsNamePion.Data()));

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
  AliAnalysisDataContainer *couputTrkCutsPionMC;
  AliAnalysisDataContainer *couputAntiTrkCutsPionMC;
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
    couputTrkCutsPionMC = mgr->CreateContainer(
        //@suppress("Invalid arguments") it works ffs
        v0CutsMCName.Data(),
        TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), v0CutsMCName.Data()));

    TString Antiv0CutsMCName = Form("%sAntiv0CutsMC_%s_%s", addon.Data(), suffixTrigger.Data(), suffix.Data());
    couputAntiTrkCutsPionMC = mgr->CreateContainer(
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




  AliAnalysisTaskThreeBodyFemtoAODPionProton* taskAOD;
  AliAnalysisTaskThreeBodyFemtoPionProton* taskNano;

  if (isNano) {
    taskNano= new AliAnalysisTaskThreeBodyFemtoPionProton("femtoThreeBodyPionProton", isMC, triggerOn);
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
    taskNano->SetPionCuts(fTrackCutsPosPion);  
    taskNano->SetAntiPionCuts(fTrackCutsNegPion);  
    if(triggerOn){
      taskNano->SetQ3Limit(Q3Limit);
    }
    taskNano->SetCorrelationConfig(config); 
    taskNano->SetRunThreeBodyHistograms(true);
    taskNano->SetTriggerOn(triggerOn);
    taskNano->SetIsMC(isMC);
    mgr->AddTask(taskNano); 

    mgr->ConnectInput(taskNano, 0, cinput); 
    mgr->ConnectOutput(taskNano, 1, coutputEvtCuts);  
    mgr->ConnectOutput(taskNano, 2, couputTrkCuts); 
    mgr->ConnectOutput(taskNano, 3, coutputAntiTrkCuts);  
    mgr->ConnectOutput(taskNano, 4, couputTrkCutsPion); 
    mgr->ConnectOutput(taskNano, 5, couputAntiTrkCutsPion); 
    mgr->ConnectOutput(taskNano, 6, coutputResults);  
    mgr->ConnectOutput(taskNano, 7, coutputResultsQA);  
    mgr->ConnectOutput(taskNano, 8, coutputResultsSample);  
    mgr->ConnectOutput(taskNano, 9, coutputResultsSampleQA); 
    mgr->ConnectOutput(taskNano, 10, coutputThreeBody); 
    if (isMC) { 
      mgr->ConnectOutput(taskNano, 11, coutputTrkCutsMC); 
      mgr->ConnectOutput(taskNano, 12, coutputAntiTrkCutsMC); 
      mgr->ConnectOutput(taskNano, 13, couputTrkCutsPionMC);  
      mgr->ConnectOutput(taskNano, 14, couputAntiTrkCutsPionMC);  
    } 

  }
  else{
    taskAOD= new AliAnalysisTaskThreeBodyFemtoAODPionProton("femtoAODThreeBodyPionProton", isMC, triggerOn);
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
    taskAOD->SetPionCuts(fTrackCutsPosPion);  
    taskAOD->SetAntiPionCuts(fTrackCutsNegPion);  
    if(triggerOn){
      taskAOD->SetQ3Limit(Q3Limit);
    }
    taskAOD->SetCorrelationConfig(config); 
    taskAOD->SetRunThreeBodyHistograms(true);
    taskAOD->SetTriggerOn(triggerOn);
    taskAOD->SetIsMC(isMC);
    mgr->AddTask(taskAOD); 

    mgr->ConnectInput(taskAOD, 0, cinput); 
    mgr->ConnectOutput(taskAOD, 1, coutputEvtCuts);  
    mgr->ConnectOutput(taskAOD, 2, couputTrkCuts); 
    mgr->ConnectOutput(taskAOD, 3, coutputAntiTrkCuts);  
    mgr->ConnectOutput(taskAOD, 4, couputTrkCutsPion); 
    mgr->ConnectOutput(taskAOD, 5, couputAntiTrkCutsPion); 
    mgr->ConnectOutput(taskAOD, 6, coutputResults);  
    mgr->ConnectOutput(taskAOD, 7, coutputResultsQA);  
    mgr->ConnectOutput(taskAOD, 8, coutputResultsSample);  
    mgr->ConnectOutput(taskAOD, 9, coutputResultsSampleQA); 
    mgr->ConnectOutput(taskAOD, 10, coutputThreeBody); 
    if (isMC) { 
      mgr->ConnectOutput(taskAOD, 11, coutputTrkCutsMC); 
      mgr->ConnectOutput(taskAOD, 12, coutputAntiTrkCutsMC); 
      mgr->ConnectOutput(taskAOD, 13, couputTrkCutsPionMC);  
      mgr->ConnectOutput(taskAOD, 14, couputAntiTrkCutsPionMC);  
    } 
  }

      
  if (isNano) {
    return taskNano;
  } else {
    return taskAOD;
  }

  
}
