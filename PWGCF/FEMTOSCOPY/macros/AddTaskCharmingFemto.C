#ifndef __CINT__
#include <vector>
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliFemtoDreamEventCuts.h"
#include "AliFemtoDreamTrackCuts.h"
#include "AliFemtoDreamCollConfig.h"
#include "AliAnalysisTaskCharmingFemto.h"
#include "AliRDHFCuts.h"
#include "AliRDHFCutsDplustoKpipi.h"
#endif

AliAnalysisTaskSE *AddTaskCharmingFemto(
    bool isMC = false, bool fullBlastQA = true, TString trigger = "kINT7",
    int channelHF = AliAnalysisTaskCharmingFemto::kDplustoKpipi,
    TString fileCutObjHF = "HFCuts.root", TString cutObjHFName = "AnalysisCuts",
    bool applyML = false, TString configML = "config_ML.yml",
    int useAODProtection = 0, const char *cutVariation = "0") {
  TString suffix = TString::Format("%s", cutVariation);

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskCharmingFemto()", "No analysis manager found.");
    return nullptr;
  }

  // ================== GetInputEventHandler =============================
  AliVEventHandler *inputHandler = mgr->GetInputEventHandler();
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();

  //========= Init subtasks and start analyis ============================
  // Event Cuts
  AliFemtoDreamEventCuts *evtCuts = AliFemtoDreamEventCuts::StandardCutsRun2();
  evtCuts->CleanUpMult(false, false, false, true);
  evtCuts->SetSphericityCuts(0.7, 1.0);

  if (suffix == "1") {
    evtCuts->SetSphericityCuts(0., 1.);
  }

  // =====================================================================
  // Proton cut variations
  const float ProtonPtlow = 0.4;
  const float ProtonPtup = 0.6;
  const float ProtonEtaLow = 0.75;
  const float ProtonEtaUp = 0.85;
  const float ProtonNsigmaLow = 2.5;
  const float ProtonNsigmaUp = 3.5;
  const float ProtonNClsLow = 70;
  const float ProtonNClsUp = 90;

  // Track Cuts
  AliFemtoDreamTrackCuts *TrackCuts = AliFemtoDreamTrackCuts::PrimProtonCuts(
      isMC, true, false, false);
  TrackCuts->SetFilterBit(128);
  TrackCuts->SetCutCharge(1);

  AliFemtoDreamTrackCuts *AntiTrackCuts =
      AliFemtoDreamTrackCuts::PrimProtonCuts(isMC, true, false, false);
  AntiTrackCuts->SetFilterBit(128);
  AntiTrackCuts->SetCutCharge(-1);

  TrackCuts->SetMinimalBooking(suffix != "0");
  AntiTrackCuts->SetMinimalBooking(suffix != "0");

  // =====================================================================
  // D mesons
  TFile* fileCuts = TFile::Open(fileCutObjHF.Data());
  if(!fileCuts ||(fileCuts&& !fileCuts->IsOpen())){
      Error("AddTaskCharmingFemto()", "Input HF cut object file not found.");
      return nullptr;
  }

  AliRDHFCuts *analysisCutsHF = nullptr;
  TString HFPartName = "";
  switch(channelHF) {
    case AliAnalysisTaskCharmingFemto::kDplustoKpipi:
      HFPartName = "Dplus";
      analysisCutsHF = (AliRDHFCutsDplustoKpipi*)fileCuts->Get(cutObjHFName);
    break;
    default:
      Error("AddTaskCharmingFemto()", "Wrong HF hadron setting, particle not implemented.");
      return nullptr;
    break;
  }


  // =====================================================================
  // Femto Collection
  std::vector<int> PDGParticles;
  PDGParticles.push_back(2212);  // 0 - proton
  PDGParticles.push_back(2212);  // 1 - antiproton
  PDGParticles.push_back(411);   // 2 - dplus
  PDGParticles.push_back(411);   // 3 - dminus

  std::vector<float> ZVtxBins = AliFemtoDreamCollConfig::GetDefaultZbins();

  std::vector<int> NBins;
  std::vector<float> kMin;
  std::vector<float> kMax;
  std::vector<int> pairQA;
  std::vector<bool> closeRejection;
  const int nPairs = 10;
  for (int i = 0; i < nPairs; ++i) {
    pairQA.push_back(0);
    closeRejection.push_back(false);
    NBins.push_back(600);
    kMin.push_back(0.);
    kMax.push_back(3.);
  }

  closeRejection[0] = true;  // pp
  closeRejection[4] = true;  // barp barp

  pairQA[0] = 11;   // pp
  pairQA[4] = 11;   // pbarpbar
  pairQA[2] = 13;   // pDplus
  pairQA[3] = 13;   // pDminus
  pairQA[5] = 13;   // barp Dplus
  pairQA[6] = 13;   // barp Dminus

  AliFemtoDreamCollConfig *config = new AliFemtoDreamCollConfig("Femto",
                                                                "Femto");
  if (trigger == "kHighMultV0") {
    std::vector<int> MultBins = AliFemtoDreamCollConfig::GetHMMultBins();
    config->SetMultBins(MultBins);
  } else if (isMC) {  // Min bias trigger in MC, but we of course want the HM binning for the CF
    std::vector<int> MultBins = AliFemtoDreamCollConfig::GetHMMultBins();
    config->SetMultBins(MultBins);
  } else {
    std::vector<int> MultBins = AliFemtoDreamCollConfig::GetMBMultBins();
    config->SetMultBins(MultBins);
  }
  config->SetMultBinning(true);

  config->SetExtendedQAPairs(pairQA);
  config->SetZBins(ZVtxBins);
  config->SetMomentumResolution(isMC);

  config->SetDeltaEtaMax(0.012);
  config->SetDeltaPhiMax(0.012);
  config->SetClosePairRejection(closeRejection);

  config->SetPhiEtaBinnign((suffix == "0" && fullBlastQA));
  config->SetmTBinning((suffix == "0" && fullBlastQA));
  config->SetdPhidEtaPlots((suffix == "0" && fullBlastQA));

  config->SetPtQA((suffix == "0" && fullBlastQA));
  config->SetMassQA((suffix == "0" && fullBlastQA));
  config->SetPDGCodes(PDGParticles);
  config->SetNBinsHist(NBins);
  config->SetMinKRel(kMin);
  config->SetMaxKRel(kMax);
  config->SetMixingDepth(10);
  config->SetUseEventMixing(true);
  config->SetMultiplicityEstimator(AliFemtoDreamEvent::kRef08);
  config->SetMinimalBookingME(false);

  if (suffix == "2") {
    config->SetPtQA(true);
    config->SetMassQA(true);
  }

  if (isMC && (suffix == "11" || suffix == "12" || suffix == "13") ) {
    config->SetdPhidEtaPlots(true);
  }

  AliAnalysisTaskCharmingFemto *task = new AliAnalysisTaskCharmingFemto(
      "AliAnalysisTaskCharmingFemto", isMC);
  task->SetLightweight(suffix != "0");
  task->SetEventCuts(evtCuts);
  task->SetProtonCuts(TrackCuts);
  task->SetAntiProtonCuts(AntiTrackCuts);
  task->SetCollectionConfig(config);
  task->SetDecayChannel(channelHF);
  task->SetHFCuts(analysisCutsHF);
  task->SetAODMismatchProtection(useAODProtection);
  if(applyML) {
    task->SetDoMLApplication(applyML);
    task->SetMLConfigFile(configML);
  }

  if (trigger == "kINT7") {
    task->SelectCollisionCandidates(AliVEvent::kINT7);
    task->SetTrigger(AliVEvent::kINT7);
  } else if (trigger == "kHighMultV0") {
    task->SelectCollisionCandidates(AliVEvent::kHighMultV0);
    task->SetTrigger(AliVEvent::kHighMultV0);
  }

  if (suffix == "2") {
    task->SetNSigmaSelection(3);
  } else if (suffix == "3") {
    task->SetMassWindow(1.9, 1.95);  // upper sideband, 5 sigma away from the peak - far off the D*
  } else if (suffix == "4") {
    task->SetMassWindow(1.95, 1.98);  // upper sideband
  } else if (suffix == "5") {
    task->SetMassWindow(1.98, 2.05);  // around D*
  } else if (suffix == "6") {
    task->SetMassWindow(2.05, 2.1);  // upper sideband
  } else if (suffix == "7") {
    task->SetMassWindow(1.79, 1.84);  // lower sideband, 5 sigma away from the peak
  } else if (suffix == "8") {
    task->SetMassWindow(1.74, 1.79);  // lower sideband
  } else if (suffix == "9") {
    task->SetMassWindow(1.69, 1.74);  // lower sideband
  } else if (suffix == "10") {
    task->SetMassWindow(1.64, 1.69);  // lower sideband
  }

  if (isMC) {
    task->ScaleMCBeautyFraction(0.5, 0.05);
    if (suffix == "11") {
      task->ScaleMCBeautyFraction(0.5, 0.);
    } else if (suffix == "12") {
      task->ScaleMCBeautyFraction(0.5, 0.5);
    } else if (suffix == "13") {
      task->UseTrueDOnly();
    }
  }

  mgr->AddTask(task);

  TString addon = "";
  if (trigger == "kINT7") {
    addon += "MB_CharmFemto_";
  } else if (trigger == "kHighMultV0") {
    addon += "HM_CharmFemto_";
  }

  TString file = AliAnalysisManager::GetCommonFileName();

  mgr->ConnectInput(task, 0, cinput);

  TString QAName = Form("%sQA%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *coutputQA = mgr->CreateContainer(
      QAName.Data(), TList::Class(), AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), QAName.Data()));
  mgr->ConnectOutput(task, 1, coutputQA);

  TString EvtCutsName = Form("%sEvtCuts%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *coutputEvtCuts = mgr->CreateContainer(
      EvtCutsName.Data(), TList::Class(), AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), EvtCutsName.Data()));
  mgr->ConnectOutput(task, 2, coutputEvtCuts);

  TString TrackCutsName = Form("%sTrackCuts%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *couputTrkCuts = mgr->CreateContainer(
      TrackCutsName.Data(), TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), TrackCutsName.Data()));
  mgr->ConnectOutput(task, 3, couputTrkCuts);

  TString AntiTrackCutsName = Form("%sAntiTrackCuts%s", addon.Data(),
                                   suffix.Data());
  AliAnalysisDataContainer *coutputAntiTrkCuts = mgr->CreateContainer(
      AntiTrackCutsName.Data(), TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), AntiTrackCutsName.Data()));
  mgr->ConnectOutput(task, 4, coutputAntiTrkCuts);

  AliAnalysisDataContainer *coutputDplus;
  TString DplusName = Form("%sDChargedQA%s", addon.Data(), suffix.Data());
  coutputDplus = mgr->CreateContainer(
		  DplusName.Data(), TList::Class(), AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), DplusName.Data()));
  mgr->ConnectOutput(task, 5, coutputDplus);

  AliAnalysisDataContainer *coutputResults;
  TString ResultsName = Form("%sResults%s", addon.Data(), suffix.Data());
  coutputResults = mgr->CreateContainer(
      ResultsName.Data(), TList::Class(), AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), ResultsName.Data()));
  mgr->ConnectOutput(task, 6, coutputResults);

  AliAnalysisDataContainer *coutputResultQA;
  TString ResultQAName = Form("%sResultQA%s", addon.Data(), suffix.Data());
  coutputResultQA = mgr->CreateContainer(
      ResultQAName.Data(), TList::Class(), AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), ResultQAName.Data()));
  mgr->ConnectOutput(task, 7, coutputResultQA);

  AliAnalysisDataContainer *coutputCutObjHF;
  TString CutObjHFName = Form("%sCutObject%s%s", addon.Data(), HFPartName.Data(), suffix.Data());
  coutputCutObjHF = mgr->CreateContainer(
      CutObjHFName.Data(), AliRDHFCuts::Class(), AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), CutObjHFName.Data()));
  mgr->ConnectOutput(task, 8, coutputCutObjHF);

  if (isMC) {
    TString TrkCutsMCName = Form("%sTrackCutsMC%s", addon.Data(),
                                 suffix.Data());
    AliAnalysisDataContainer *coutputTrkCutsMC = mgr->CreateContainer(
        TrkCutsMCName.Data(), TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), TrkCutsMCName.Data()));
    mgr->ConnectOutput(task, 9, coutputTrkCutsMC);

    TString AntiTrkCutsMCName = Form("%sAntiTrackCutsMC%s", addon.Data(),
                                     suffix.Data());
    AliAnalysisDataContainer *coutputAntiTrkCutsMC = mgr->CreateContainer(
        AntiTrkCutsMCName.Data(), TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), AntiTrkCutsMCName.Data()));
    mgr->ConnectOutput(task, 10, coutputAntiTrkCutsMC);
  }

  return task;
}
