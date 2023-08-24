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

AliAnalysisTaskSE *AddTaskCharmingFemtoSyst(
    bool isMC = false,
    bool fullBlastQA = true,
    TString trigger = "kINT7",
    int channelHF = AliAnalysisTaskCharmingFemto::kDplustoKpipi,
    TString fileCutObjHF = "HFCuts.root",
    TString cutObjHFName = "AnalysisCuts",
    TString cutHFsuffix = "",
    bool applyML = false,
    TString configML = "config_ML.yml",
    int useAODProtection = 0,
    int massSelection = AliAnalysisTaskCharmingFemto::kSignal,
    unsigned int pdgCodeBuddy = 211,
    const char *cutVariation = "0"
  ) {
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

  // =====================================================================
  // Cut variations
  float buddyPtlow;
  float buddyPtup;
  float buddyEtaLow = 0.75;
  float buddyEtaUp = 0.85;
  float buddyNsigmaLow = 2.7;
  float buddyNsigmaUp = 3.3;
  float buddyNClsLow = 70;
  float buddyNClsUp = 90;
  float buddyMaxPt;

  std::map<std::string, float> kaonPIDTight;
  std::map<std::string, float> kaonPIDLoose;
  
  AliPID::EParticleType aliPIDParticle;
  
  if (pdgCodeBuddy == 2212){ // proton
    buddyPtlow = 0.45;
    buddyPtup = 0.55;
    buddyNsigmaLow = 2.7;
    buddyNsigmaUp = 3.3;
    buddyMaxPt = 4.05
    aliPIDParticle = AliPID::kProton;
  } else if (pdgCodeBuddy == 211){ // pion+
    buddyPtlow = 0.1;
    buddyPtup = 0.2;
    buddyNsigmaLow = 2.7;
    buddyNsigmaUp = 3.3;
    buddyMaxPt = 4.0;
    aliPIDParticle = AliPID::kPion;
  } else if (pdgCodeBuddy == 321){ // kaon+
    buddyPtlow = 0.09;
    buddyPtup = 0.19;
    buddyNsigmaLow = 3;
    buddyNsigmaUp = 3;
    buddyMaxPt = 4.0;
    aliPIDParticle = AliPID::kKaon;
    kaonPIDTight = { {"COMB", 3.7}, {"TPC", 2.7}, {"EXCLUSION", 3.3}, }; // for SetPIDkd() when using oton's K selection
    kaonPIDLoose = { {"COMB", 4.3}, {"TPC", 3.3}, {"EXCLUSION", 2.7}, };
  } 

  // Track Cuts
  AliFemtoDreamTrackCuts *TrackCuts = nullptr;
  AliFemtoDreamTrackCuts *AntiTrackCuts = nullptr;

  if(std::abs(pdgCodeBuddy) == 2212) {
    TrackCuts = AliFemtoDreamTrackCuts::PrimProtonCuts(isMC, true, false, false);
    TrackCuts->SetFilterBit(128);
    TrackCuts->SetCutCharge(1);
    AntiTrackCuts = AliFemtoDreamTrackCuts::PrimProtonCuts(isMC, true, false, false);
    AntiTrackCuts->SetFilterBit(128);
    AntiTrackCuts->SetCutCharge(-1);
  }
  else if(std::abs(pdgCodeBuddy) == 211) {
    TrackCuts = AliFemtoDreamTrackCuts::PrimPionCuts(isMC, true, false, false);
    TrackCuts->SetFilterBit(96);
    TrackCuts->SetCutCharge(1);
    AntiTrackCuts = AliFemtoDreamTrackCuts::PrimPionCuts(isMC, true, false, false);
    AntiTrackCuts->SetFilterBit(96);
    AntiTrackCuts->SetCutCharge(-1);
  }
  else if(std::abs(pdgCodeBuddy) == 321) {
    TrackCuts = AliFemtoDreamTrackCuts::PrimKaonCuts(isMC, true, false, false);
    TrackCuts->SetFilterBit(128);
    TrackCuts->SetPtRange(0.15, 4.0);
    TrackCuts->SetCutCharge(1);
    if(!useMCTruthReco) {
      TrackCuts->SetPIDkd();
    }
    AntiTrackCuts = AliFemtoDreamTrackCuts::PrimKaonCuts(isMC, true, false, false);
    AntiTrackCuts->SetFilterBit(128);
    AntiTrackCuts->SetPtRange(0.15, 4.0);
    AntiTrackCuts->SetCutCharge(-1);
    if(!useMCTruthReco) {
      AntiTrackCuts->SetPIDkd();
    }
  }
  else {
    Error("AddTaskAnyCharmingFemto()", "Particle not implemented.");
    return nullptr;
  }

  TrackCuts->SetMinimalBooking(suffix != "0");
  AntiTrackCuts->SetMinimalBooking(suffix != "0");

  if (suffix == "1") {
    TrackCuts->SetPtRange(buddyPtlow, buddyMaxPt);
    AntiTrackCuts->SetPtRange(buddyPtlow, buddyMaxPt);
    TrackCuts->SetEtaRange(-buddyEtaLow, buddyEtaLow);
    AntiTrackCuts->SetEtaRange(-buddyEtaLow, buddyEtaLow);
  } else if (suffix == "2") {
    TrackCuts->SetPtRange(buddyPtup, buddyMaxPt);
    AntiTrackCuts->SetPtRange(buddyPtup, buddyMaxPt);
    
    // Set PID variations
    if (aliPIDparticle == AliPID::kProton){
      TrackCuts->SetPID(aliPIDParticle, 0.75, buddyNsigmaLow);
      AntiTrackCuts->SetPID(aliPIDParticle, 0.75, buddyNsigmaLow);
    } else if (aliPIDparticle == AliPID::kPion){
      TrackCuts->SetPID(aliPIDParticle, 0.5, buddyNsigmaLow);
      AntiTrackCuts->SetPID(aliPIDParticle, 0.5, buddyNsigmaLow);
    } else if (aliPIDparticle == AliPID::kKaon){
      TrackCuts->SetPIDkd(true, false, kaonPIDLoose["COMB"], kaonPIDTight["TPC"], kaonPIDTight["EXCLUSION"]);
      AntiTrackCuts->SetPIDkd(true, false, kaonPIDLoose["COMB"], kaonPIDTight["TPC"], kaonPIDTight["EXCLUSION"]);
    }
  } else if (suffix == "3") {
    TrackCuts->SetEtaRange(-buddyEtaUp, buddyEtaUp);
    AntiTrackCuts->SetEtaRange(-buddyEtaUp, buddyEtaUp);
    TrackCuts->SetNClsTPC(buddyNClsUp);
    AntiTrackCuts->SetNClsTPC(buddyNClsUp);
  } else if (suffix == "4") {
    TrackCuts->SetPtRange(buddyPtup, buddyMaxPt);
    AntiTrackCuts->SetPtRange(buddyPtup, buddyMaxPt);
    
    // Set PID variations
    if (aliPIDparticle == AliPID::kProton){
      TrackCuts->SetPID(aliPIDParticle, 0.75, buddyNsigmaUp);
      AntiTrackCuts->SetPID(aliPIDParticle, 0.75, buddyNsigmaUp);
    } else if (aliPIDparticle == AliPID::kPion){
      TrackCuts->SetPID(aliPIDParticle, 0.5, buddyNsigmaUp);
      AntiTrackCuts->SetPID(aliPIDParticle, 0.5, buddyNsigmaUp);
    } else if (aliPIDparticle == AliPID::kKaon){
      TrackCuts->SetPIDkd(true, false, kaonPIDTight["COMB"], kaonPIDTight["TPC"], kaonPIDTight["EXCLUSION"]);
      AntiTrackCuts->SetPIDkd(true, false, kaonPIDTight["COMB"], kaonPIDTight["TPC"], kaonPIDTight["EXCLUSION"]);
    }
  } else if (suffix == "5") {
    TrackCuts->SetEtaRange(-buddyEtaLow, buddyEtaLow);
    AntiTrackCuts->SetEtaRange(-buddyEtaLow, buddyEtaLow);
    TrackCuts->SetNClsTPC(buddyNClsLow);
    AntiTrackCuts->SetNClsTPC(buddyNClsLow);
  } else if (suffix == "6") {
    TrackCuts->SetPtRange(buddyPtup, buddyMaxPt);
    AntiTrackCuts->SetPtRange(buddyPtup, buddyMaxPt);
    
    // Set PID variations
    if (aliPIDparticle == AliPID::kProton){
      TrackCuts->SetPID(aliPIDParticle, 0.75, buddyNsigmaLow);
      AntiTrackCuts->SetPID(aliPIDParticle, 0.75, buddyNsigmaLow);
    } else if (aliPIDparticle == AliPID::kPion){
      TrackCuts->SetPID(aliPIDParticle, 0.5, buddyNsigmaLow);
      AntiTrackCuts->SetPID(aliPIDParticle, 0.5, buddyNsigmaLow);
    } else if (aliPIDparticle == AliPID::kKaon){
      TrackCuts->SetPIDkd(true, false, kaonPIDLoose["COMB"], kaonPIDTight["TPC"], kaonPIDTight["EXCLUSION"]);
      AntiTrackCuts->SetPIDkd(true, false, kaonPIDLoose["COMB"], kaonPIDTight["TPC"], kaonPIDTight["EXCLUSION"]);
    }
  } else if (suffix == "7") {
    TrackCuts->SetEtaRange(-buddyEtaLow, buddyEtaLow);
    AntiTrackCuts->SetEtaRange(-buddyEtaLow, buddyEtaLow);
    TrackCuts->SetNClsTPC(buddyNClsUp);
    AntiTrackCuts->SetNClsTPC(buddyNClsUp);
  } else if (suffix == "8") {
    TrackCuts->SetPtRange(buddyPtup, buddyMaxPt);
    AntiTrackCuts->SetPtRange(buddyPtup, buddyMaxPt);
    
    // Set PID variations
    if (aliPIDparticle == AliPID::kProton){
      TrackCuts->SetPID(aliPIDParticle, 0.75, buddyNsigmaUp);
      AntiTrackCuts->SetPID(aliPIDParticle, 0.75, buddyNsigmaUp);
    } else if (aliPIDparticle == AliPID::kPion){
      TrackCuts->SetPID(aliPIDParticle, 0.5, buddyNsigmaUp);
      AntiTrackCuts->SetPID(aliPIDParticle, 0.5, buddyNsigmaUp);
    } else if (aliPIDparticle == AliPID::kKaon){
      TrackCuts->SetPIDkd(true, false, kaonPIDTight["COMB"], kaonPIDTight["TPC"], kaonPIDTight["EXCLUSION"]);
      AntiTrackCuts->SetPIDkd(true, false, kaonPIDTight["COMB"], kaonPIDTight["TPC"], kaonPIDTight["EXCLUSION"]);
    }
  } else if (suffix == "9") {
    TrackCuts->SetPtRange(buddyPtup, buddyMaxPt);
    AntiTrackCuts->SetPtRange(buddyPtup, buddyMaxPt);
    TrackCuts->SetNClsTPC(buddyNClsUp);
    AntiTrackCuts->SetNClsTPC(buddyNClsUp);
  } else if (suffix == "10") {
    
    // Set PID variations
    if (aliPIDparticle == AliPID::kProton){
      TrackCuts->SetPID(aliPIDParticle, 0.75, buddyNsigmaLow);
      AntiTrackCuts->SetPID(aliPIDParticle, 0.75, buddyNsigmaLow);
    } else if (aliPIDparticle == AliPID::kPion){
      TrackCuts->SetPID(aliPIDParticle, 0.5, buddyNsigmaLow);
      AntiTrackCuts->SetPID(aliPIDParticle, 0.5, buddyNsigmaLow);
    } else if (aliPIDparticle == AliPID::kKaon){
      TrackCuts->SetPIDkd(true, false, kaonPIDLoose["COMB"], kaonPIDTight["TPC"], kaonPIDTight["EXCLUSION"]);
      AntiTrackCuts->SetPIDkd(true, false, kaonPIDLoose["COMB"], kaonPIDTight["TPC"], kaonPIDTight["EXCLUSION"]);
    }
    TrackCuts->SetNClsTPC(buddyNClsLow);
    AntiTrackCuts->SetNClsTPC(buddyNClsLow);
  } else if (suffix == "11") {
    TrackCuts->SetPtRange(buddyPtup, buddyMaxPt);
    AntiTrackCuts->SetPtRange(buddyPtup, buddyMaxPt);
    
    // Set PID variations
    if (aliPIDparticle == AliPID::kProton){
      TrackCuts->SetPID(aliPIDParticle, 0.75, buddyNsigmaLow);
      AntiTrackCuts->SetPID(aliPIDParticle, 0.75, buddyNsigmaLow);
    } else if (aliPIDparticle == AliPID::kPion){
      TrackCuts->SetPID(aliPIDParticle, 0.5, buddyNsigmaLow);
      AntiTrackCuts->SetPID(aliPIDParticle, 0.5, buddyNsigmaLow);
    } else if (aliPIDparticle == AliPID::kKaon){
      TrackCuts->SetPIDkd(true, false, kaonPIDLoose["COMB"], kaonPIDTight["TPC"], kaonPIDTight["EXCLUSION"]);
      AntiTrackCuts->SetPIDkd(true, false, kaonPIDLoose["COMB"], kaonPIDTight["TPC"], kaonPIDTight["EXCLUSION"]);
    }
  } else if (suffix == "12") {
    TrackCuts->SetPtRange(buddyPtup, buddyMaxPt);
    AntiTrackCuts->SetPtRange(buddyPtup, buddyMaxPt);
    TrackCuts->SetNClsTPC(buddyNClsUp);
    AntiTrackCuts->SetNClsTPC(buddyNClsUp);
  } else if (suffix == "13") {
    TrackCuts->SetEtaRange(-buddyEtaUp, buddyEtaUp);
    AntiTrackCuts->SetEtaRange(-buddyEtaUp, buddyEtaUp);
    TrackCuts->SetPtRange(buddyPtlow, buddyMaxPt);
    AntiTrackCuts->SetPtRange(buddyPtlow, buddyMaxPt);
  } else if (suffix == "14") {
    
    // Set PID variations
    if (aliPIDparticle == AliPID::kProton){
      TrackCuts->SetPID(aliPIDParticle, 0.75, buddyNsigmaLow);
      AntiTrackCuts->SetPID(aliPIDParticle, 0.75, buddyNsigmaLow);
    } else if (aliPIDparticle == AliPID::kPion){
      TrackCuts->SetPID(aliPIDParticle, 0.5, buddyNsigmaLow);
      AntiTrackCuts->SetPID(aliPIDParticle, 0.5, buddyNsigmaLow);
    } else if (aliPIDparticle == AliPID::kKaon){
      TrackCuts->SetPIDkd(true, false, kaonPIDLoose["COMB"], kaonPIDTight["TPC"], kaonPIDTight["EXCLUSION"]);
      AntiTrackCuts->SetPIDkd(true, false, kaonPIDLoose["COMB"], kaonPIDTight["TPC"], kaonPIDTight["EXCLUSION"]);
    }
    TrackCuts->SetEtaRange(-buddyEtaLow, buddyEtaLow);
    AntiTrackCuts->SetEtaRange(-buddyEtaLow, buddyEtaLow);
  } else if (suffix == "15") {
    TrackCuts->SetNClsTPC(buddyNClsUp);
    AntiTrackCuts->SetNClsTPC(buddyNClsUp);
    TrackCuts->SetPtRange(buddyPtup, buddyMaxPt);
    AntiTrackCuts->SetPtRange(buddyPtup, buddyMaxPt);
  } else if (suffix == "16") {
    TrackCuts->SetEtaRange(-buddyEtaUp, buddyEtaUp);
    AntiTrackCuts->SetEtaRange(-buddyEtaUp, buddyEtaUp);
    
    // Set PID variations
    if (aliPIDparticle == AliPID::kProton){
      TrackCuts->SetPID(aliPIDParticle, 0.75, buddyNsigmaLow);
      AntiTrackCuts->SetPID(aliPIDParticle, 0.75, buddyNsigmaLow);
    } else if (aliPIDparticle == AliPID::kPion){
      TrackCuts->SetPID(aliPIDParticle, 0.5, buddyNsigmaLow);
      AntiTrackCuts->SetPID(aliPIDParticle, 0.5, buddyNsigmaLow);
    } else if (aliPIDparticle == AliPID::kKaon){
      TrackCuts->SetPIDkd(true, false, kaonPIDLoose["COMB"], kaonPIDTight["TPC"], kaonPIDTight["EXCLUSION"]);
      AntiTrackCuts->SetPIDkd(true, false, kaonPIDLoose["COMB"], kaonPIDTight["TPC"], kaonPIDTight["EXCLUSION"]);
    }
  } else if (suffix == "17") {
    TrackCuts->SetPtRange(buddyPtup, buddyMaxPt);
    AntiTrackCuts->SetPtRange(buddyPtup, buddyMaxPt);
    TrackCuts->SetNClsTPC(buddyNClsLow);
    AntiTrackCuts->SetNClsTPC(buddyNClsLow);
  } else if (suffix == "18") {
    TrackCuts->SetPtRange(buddyPtup, buddyMaxPt);
    AntiTrackCuts->SetPtRange(buddyPtup, buddyMaxPt);
    
    // Set PID variations
    if (aliPIDparticle == AliPID::kProton){
      TrackCuts->SetPID(aliPIDParticle, 0.75, buddyNsigmaUp);
      AntiTrackCuts->SetPID(aliPIDParticle, 0.75, buddyNsigmaUp);
    } else if (aliPIDparticle == AliPID::kPion){
      TrackCuts->SetPID(aliPIDParticle, 0.5, buddyNsigmaUp);
      AntiTrackCuts->SetPID(aliPIDParticle, 0.5, buddyNsigmaUp);
    } else if (aliPIDparticle == AliPID::kKaon){
      TrackCuts->SetPIDkd(true, false, kaonPIDTight["COMB"], kaonPIDTight["TPC"], kaonPIDTight["EXCLUSION"]);
      AntiTrackCuts->SetPIDkd(true, false, kaonPIDTight["COMB"], kaonPIDTight["TPC"], kaonPIDTight["EXCLUSION"]);
    }
  } else if (suffix == "19") {
    TrackCuts->SetPtRange(buddyPtlow, buddyMaxPt);
    AntiTrackCuts->SetPtRange(buddyPtlow, buddyMaxPt);
    TrackCuts->SetEtaRange(-buddyEtaUp, buddyEtaUp);
    AntiTrackCuts->SetEtaRange(-buddyEtaUp, buddyEtaUp);
  } else if (suffix == "20") {
    TrackCuts->SetNClsTPC(buddyNClsLow);
    AntiTrackCuts->SetNClsTPC(buddyNClsLow);
    
    // Set PID variations
    if (aliPIDparticle == AliPID::kProton){
      TrackCuts->SetPID(aliPIDParticle, 0.75, buddyNsigmaUp);
      AntiTrackCuts->SetPID(aliPIDParticle, 0.75, buddyNsigmaUp);
    } else if (aliPIDparticle == AliPID::kPion){
      TrackCuts->SetPID(aliPIDParticle, 0.5, buddyNsigmaUp);
      AntiTrackCuts->SetPID(aliPIDParticle, 0.5, buddyNsigmaUp);
    } else if (aliPIDparticle == AliPID::kKaon){
      TrackCuts->SetPIDkd(true, false, kaonPIDTight["COMB"], kaonPIDTight["TPC"], kaonPIDTight["EXCLUSION"]);
      AntiTrackCuts->SetPIDkd(true, false, kaonPIDTight["COMB"], kaonPIDTight["TPC"], kaonPIDTight["EXCLUSION"]);
    }
  }

  // =====================================================================
  // D mesons
  TFile* fileCuts = TFile::Open(fileCutObjHF.Data());
  if(!fileCuts ||(fileCuts&& !fileCuts->IsOpen())){
      Error("AddTaskCharmingFemto()", "Input HF cut object file not found.");
      return nullptr;
  }

  AliRDHFCuts *analysisCutsHF = nullptr;
  TString HFPartName = "";
  Int_t pdgDmeson;
  switch(channelHF) {
    case AliAnalysisTaskCharmingFemto::kDplustoKpipi:
      HFPartName = "Dplus";
      analysisCutsHF = (AliRDHFCutsDplustoKpipi*)fileCuts->Get(cutObjHFName);
      pdgDmeson = 411;
      break;
    case AliAnalysisTaskCharmingFemto::kDstartoKpipi:
      HFPartName = "Dstar";
      analysisCutsHF = (AliRDHFCutsDStartoKpipi*)fileCuts->Get(cutObjHFName);
      pdgDmeson = 413;
      break;
    default:
      Error("AddTaskCharmingFemto()", "Wrong HF hadron setting, particle not implemented.");
      return nullptr;
    break;
  }

  // =====================================================================
  // Femto Collection
  std::vector<int> PDGParticles;
  PDGParticles.push_back(pdgCodeBuddy);  //
  PDGParticles.push_back(pdgCodeBuddy);  //
  PDGParticles.push_back(pdgDmeson);   // 2 - dplus or dstar+
  PDGParticles.push_back(pdgDmeson);   // 3 - dminus or dstar-

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

  config->SetmTBinning((suffix == "0" && fullBlastQA));
  config->SetPtQA(true);
  config->SetMassQA((suffix == "0" && fullBlastQA));
  config->SetPDGCodes(PDGParticles);
  config->SetNBinsHist(NBins);
  config->SetMinKRel(kMin);
  config->SetMaxKRel(kMax);
  config->SetMixingDepth(10000);
  config->SetUseEventMixing(true);
  config->SetMultiplicityEstimator(AliFemtoDreamEvent::kRef08);
  config->SetMinimalBookingME(false);

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
  task->SetMassSelection(massSelection);
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

  if (isMC) {
    task->ScaleMCBeautyFraction(0.5, 0.05);
  }

  mgr->AddTask(task);

  TString addon = "";
  if (trigger == "kINT7") {
    addon += "MB_CharmFemto_";
  } else if (trigger == "kHighMultV0") {
    addon += "HM_CharmFemto_";
  }
  if (massSelection == AliAnalysisTaskCharmingFemto::kSidebandRight) {
    addon += "SBRight_";
  } else if (massSelection == AliAnalysisTaskCharmingFemto::kSidebandLeft) {
    addon += "SBLeft_";
  }
  if(!cutHFsuffix.EqualTo("")) {
    addon += Form("%s_", cutHFsuffix.Data());
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
