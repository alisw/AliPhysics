#ifndef __CINT__
#include <vector>
#include <map>
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliFemtoDreamEventCuts.h"
#include "AliFemtoDreamTrackCuts.h"
#include "AliFemtoDreamCollConfig.h"
#include "AliAnalysisTaskCharmingFemto.h"
#include "AliRDHFCuts.h"
#include "AliRDHFCutsDplustoKpipi.h"
#include "AliRDHFCutsDStartoKpipi.h"
#include "TTree.h"
#endif

AliAnalysisTaskSE *AddTaskAnyCharmingFemto(
    bool isMC = false,
    bool useMCTruthReco = false,
    bool isMCtruth = false,
    bool useTree = false,
    bool useLooseSelections = false,
    bool fullBlastQA = true,
    TString trigger = "kINT7",
    int channelHF = AliAnalysisTaskCharmingFemto::kDplustoKpipi,
    TString fileCutObjHF = "HFCuts.root",
    TString cutObjHFName = "AnalysisCuts",
    TString cutHFsuffix = "",
    bool applyML = false, TString configML = "config_ML.yml",
    int useAODProtection = 0,
    AliAnalysisTaskCharmingFemto::MassSelectionType massSelection = AliAnalysisTaskCharmingFemto::kSignal,
    int pdgBuddy = 2212,
    int mixingDepth = 10,
    std::vector<std::string> colsToSave = {
      "mult",
      "kStar",
      "is_oldpcrm",
      "is_newpcrm",
      "heavy_mult",
      "heavy_invmass",
      "heavy_pt",
      "heavy_origin",
      "light_mult",
      "light_px",
      "light_py",
      "light_eta",
      "light_nsigtpc",
      "light_nsigtof",
      "light_dcaxy",
      "light_dcaz"},
    const char *cutVariation = "0"
  ) {

  TString suffix = TString::Format("%s", cutVariation);

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskAnyCharmingFemto()", "No analysis manager found.");
    return nullptr;
  }

  // ================== GetInputEventHandler =============================
  AliVEventHandler *inputHandler = mgr->GetInputEventHandler();
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();

  //========= Init subtasks and start analyis ============================
  // Event Cuts
  AliFemtoDreamEventCuts *evtCuts = AliFemtoDreamEventCuts::StandardCutsRun2();
  evtCuts->CleanUpMult(false, false, false, true);

  // Track Cuts
  AliFemtoDreamTrackCuts *TrackCuts = nullptr;
  AliFemtoDreamTrackCuts *AntiTrackCuts = nullptr;

  // overwritten later if useTree == true
  if(std::abs(pdgBuddy) == 2212) {
    TrackCuts = AliFemtoDreamTrackCuts::PrimProtonCuts(isMC, true, false, false);
    TrackCuts->SetFilterBit(128);
    TrackCuts->SetCutCharge(1);
    AntiTrackCuts = AliFemtoDreamTrackCuts::PrimProtonCuts(isMC, true, false, false);
    AntiTrackCuts->SetFilterBit(128);
    AntiTrackCuts->SetCutCharge(-1);
  }
  else if(std::abs(pdgBuddy) == 211) {
    TrackCuts = AliFemtoDreamTrackCuts::PrimPionCuts(isMC, true, false, false);
    TrackCuts->SetFilterBit(96);
    TrackCuts->SetCutCharge(1);
    TrackCuts->SetPtRange(0.14, 4.0);

    AntiTrackCuts = AliFemtoDreamTrackCuts::PrimPionCuts(isMC, true, false, false);
    AntiTrackCuts->SetFilterBit(96);
    AntiTrackCuts->SetCutCharge(-1);
    AntiTrackCuts->SetPtRange(0.14, 4.0);
  }
  else if(std::abs(pdgBuddy) == 321) {
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
  
  if (pdgBuddy == 2212){ // proton
    buddyPtlow = 0.45;
    buddyPtup = 0.55;
    buddyNsigmaLow = 2.7;
    buddyNsigmaUp = 3.3;
    buddyMaxPt = 4.05;
    aliPIDParticle = AliPID::kProton;
  } else if (pdgBuddy == 211){ // pion+
    buddyPtlow = 0.11;
    buddyPtup = 0.17;
    buddyNsigmaLow = 2.85;
    buddyNsigmaUp = 3.15;
    buddyMaxPt = 4.0;
    // buddyNClsLow = 75;
    // buddyNClsUp = 85;
    aliPIDParticle = AliPID::kPion;
  } else if (pdgBuddy == 321){ // kaon+
    buddyPtlow = 0.1;
    buddyPtup = 0.2;
    buddyNsigmaLow = 3;
    buddyNsigmaUp = 3;
    buddyMaxPt = 4.0;
    aliPIDParticle = AliPID::kKaon;
    kaonPIDTight = { {"COMB", 2.7}, {"TPC", 2.7}, {"EXCLUSION", 3.3}, }; // for SetPIDkd() when using oton's K selection
    kaonPIDLoose = { {"COMB", 3.3}, {"TPC", 3.3}, {"EXCLUSION", 2.7}, };
  } 

if (!isMC) {
  TrackCuts->SetMinimalBooking(suffix != "0");
  AntiTrackCuts->SetMinimalBooking(suffix != "0");
}

  if (suffix == "1") {
    TrackCuts->SetPtRange(buddyPtlow, buddyMaxPt);
    AntiTrackCuts->SetPtRange(buddyPtlow, buddyMaxPt);
    TrackCuts->SetEtaRange(-buddyEtaLow, buddyEtaLow);
    AntiTrackCuts->SetEtaRange(-buddyEtaLow, buddyEtaLow);
  } else if (suffix == "2") {
    TrackCuts->SetPtRange(buddyPtup, buddyMaxPt);
    AntiTrackCuts->SetPtRange(buddyPtup, buddyMaxPt);
    
    // Set PID variations
    if (aliPIDParticle == AliPID::kProton){
      TrackCuts->SetPID(aliPIDParticle, 0.75, buddyNsigmaLow);
      AntiTrackCuts->SetPID(aliPIDParticle, 0.75, buddyNsigmaLow);
    } else if (aliPIDParticle == AliPID::kPion){
      TrackCuts->SetPID(aliPIDParticle, 0.5, buddyNsigmaLow);
      AntiTrackCuts->SetPID(aliPIDParticle, 0.5, buddyNsigmaLow);
    } else if (aliPIDParticle == AliPID::kKaon){
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
    if (aliPIDParticle == AliPID::kProton){
      TrackCuts->SetPID(aliPIDParticle, 0.75, buddyNsigmaUp);
      AntiTrackCuts->SetPID(aliPIDParticle, 0.75, buddyNsigmaUp);
    } else if (aliPIDParticle == AliPID::kPion){
      TrackCuts->SetPID(aliPIDParticle, 0.5, buddyNsigmaUp);
      AntiTrackCuts->SetPID(aliPIDParticle, 0.5, buddyNsigmaUp);
    } else if (aliPIDParticle == AliPID::kKaon){
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
    if (aliPIDParticle == AliPID::kProton){
      TrackCuts->SetPID(aliPIDParticle, 0.75, buddyNsigmaLow);
      AntiTrackCuts->SetPID(aliPIDParticle, 0.75, buddyNsigmaLow);
    } else if (aliPIDParticle == AliPID::kPion){
      TrackCuts->SetPID(aliPIDParticle, 0.5, buddyNsigmaLow);
      AntiTrackCuts->SetPID(aliPIDParticle, 0.5, buddyNsigmaLow);
    } else if (aliPIDParticle == AliPID::kKaon){
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
    if (aliPIDParticle == AliPID::kProton){
      TrackCuts->SetPID(aliPIDParticle, 0.75, buddyNsigmaUp);
      AntiTrackCuts->SetPID(aliPIDParticle, 0.75, buddyNsigmaUp);
    } else if (aliPIDParticle == AliPID::kPion){
      TrackCuts->SetPID(aliPIDParticle, 0.5, buddyNsigmaUp);
      AntiTrackCuts->SetPID(aliPIDParticle, 0.5, buddyNsigmaUp);
    } else if (aliPIDParticle == AliPID::kKaon){
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
    if (aliPIDParticle == AliPID::kProton){
      TrackCuts->SetPID(aliPIDParticle, 0.75, buddyNsigmaLow);
      AntiTrackCuts->SetPID(aliPIDParticle, 0.75, buddyNsigmaLow);
    } else if (aliPIDParticle == AliPID::kPion){
      TrackCuts->SetPID(aliPIDParticle, 0.5, buddyNsigmaLow);
      AntiTrackCuts->SetPID(aliPIDParticle, 0.5, buddyNsigmaLow);
    } else if (aliPIDParticle == AliPID::kKaon){
      TrackCuts->SetPIDkd(true, false, kaonPIDLoose["COMB"], kaonPIDTight["TPC"], kaonPIDTight["EXCLUSION"]);
      AntiTrackCuts->SetPIDkd(true, false, kaonPIDLoose["COMB"], kaonPIDTight["TPC"], kaonPIDTight["EXCLUSION"]);
    }
    TrackCuts->SetNClsTPC(buddyNClsLow);
    AntiTrackCuts->SetNClsTPC(buddyNClsLow);
  } else if (suffix == "11") {
    TrackCuts->SetPtRange(buddyPtup, buddyMaxPt);
    AntiTrackCuts->SetPtRange(buddyPtup, buddyMaxPt);
    
    // Set PID variations
    if (aliPIDParticle == AliPID::kProton){
      TrackCuts->SetPID(aliPIDParticle, 0.75, buddyNsigmaLow);
      AntiTrackCuts->SetPID(aliPIDParticle, 0.75, buddyNsigmaLow);
    } else if (aliPIDParticle == AliPID::kPion){
      TrackCuts->SetPID(aliPIDParticle, 0.5, buddyNsigmaLow);
      AntiTrackCuts->SetPID(aliPIDParticle, 0.5, buddyNsigmaLow);
    } else if (aliPIDParticle == AliPID::kKaon){
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
    if (aliPIDParticle == AliPID::kProton){
      TrackCuts->SetPID(aliPIDParticle, 0.75, buddyNsigmaLow);
      AntiTrackCuts->SetPID(aliPIDParticle, 0.75, buddyNsigmaLow);
    } else if (aliPIDParticle == AliPID::kPion){
      TrackCuts->SetPID(aliPIDParticle, 0.5, buddyNsigmaLow);
      AntiTrackCuts->SetPID(aliPIDParticle, 0.5, buddyNsigmaLow);
    } else if (aliPIDParticle == AliPID::kKaon){
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
    if (aliPIDParticle == AliPID::kProton){
      TrackCuts->SetPID(aliPIDParticle, 0.75, buddyNsigmaLow);
      AntiTrackCuts->SetPID(aliPIDParticle, 0.75, buddyNsigmaLow);
    } else if (aliPIDParticle == AliPID::kPion){
      TrackCuts->SetPID(aliPIDParticle, 0.5, buddyNsigmaLow);
      AntiTrackCuts->SetPID(aliPIDParticle, 0.5, buddyNsigmaLow);
    } else if (aliPIDParticle == AliPID::kKaon){
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
    if (aliPIDParticle == AliPID::kProton){
      TrackCuts->SetPID(aliPIDParticle, 0.75, buddyNsigmaUp);
      AntiTrackCuts->SetPID(aliPIDParticle, 0.75, buddyNsigmaUp);
    } else if (aliPIDParticle == AliPID::kPion){
      TrackCuts->SetPID(aliPIDParticle, 0.5, buddyNsigmaUp);
      AntiTrackCuts->SetPID(aliPIDParticle, 0.5, buddyNsigmaUp);
    } else if (aliPIDParticle == AliPID::kKaon){
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
    if (aliPIDParticle == AliPID::kProton){
      TrackCuts->SetPID(aliPIDParticle, 0.75, buddyNsigmaUp);
      AntiTrackCuts->SetPID(aliPIDParticle, 0.75, buddyNsigmaUp);
    } else if (aliPIDParticle == AliPID::kPion){
      TrackCuts->SetPID(aliPIDParticle, 0.5, buddyNsigmaUp);
      AntiTrackCuts->SetPID(aliPIDParticle, 0.5, buddyNsigmaUp);
    } else if (aliPIDParticle == AliPID::kKaon){
      TrackCuts->SetPIDkd(true, false, kaonPIDTight["COMB"], kaonPIDTight["TPC"], kaonPIDTight["EXCLUSION"]);
      AntiTrackCuts->SetPIDkd(true, false, kaonPIDTight["COMB"], kaonPIDTight["TPC"], kaonPIDTight["EXCLUSION"]);
    }
  }

  // overwrite previous settings and use loose selections if trees are used
  if (useTree && useLooseSelections) {
    // eta
    TrackCuts->SetEtaRange(-0.9, 0.9);
    AntiTrackCuts->SetEtaRange(-0.9, 0.9);
    
    // ncls
    TrackCuts->SetNClsTPC(70);
    AntiTrackCuts->SetNClsTPC(70);
    
    if (aliPIDParticle == AliPID::kPion){
      TrackCuts->SetPtRange(0.11, 4);
      AntiTrackCuts->SetPtRange(0.11, 4);
      TrackCuts->SetPID(aliPIDParticle, 0.5, 3.5);
      AntiTrackCuts->SetPID(aliPIDParticle, 0.5, 3.5);
    } else if (aliPIDParticle == AliPID::kKaon){
      TrackCuts->SetPtRange(0.1, 4);
      AntiTrackCuts->SetPtRange(0.1, 4);
      TrackCuts->SetPIDkd(true, false, 3.5, 3.5, -999); // last parameter is dummy
      AntiTrackCuts->SetPIDkd(true, false, 3.5, 3.5, -999); // last parameter is dummy
    }
  }
  
  // =====================================================================
  // D mesons
  TFile* fileCuts = TFile::Open(fileCutObjHF.Data());
  if(!fileCuts ||(fileCuts&& !fileCuts->IsOpen())){
      Error("AddTaskAnyCharmingFemto()", "Input HF cut object file not found.");
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
      Error("AddTaskAnyCharmingFemto()", "Wrong HF hadron setting, particle not implemented.");
      return nullptr;
    break;
  }

  // =====================================================================
  // Femto Collection
  std::vector<int> PDGParticles;
  PDGParticles.push_back(pdgBuddy);  //
  PDGParticles.push_back(pdgBuddy);  //
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

  if (!isMCtruth) {
    closeRejection[0] = true;   // light-light
    closeRejection[4] = true;   // antilight-antilight
    closeRejection[2] = true;   // light-charm
    closeRejection[3] = true;   // light-anticharm
    closeRejection[5] = true;   // antilight charm
    closeRejection[6] = true;   // antilight anticharm
  }
  else {
    closeRejection[0] = false;  // pp
    closeRejection[4] = false;  // barp barp
  }

  pairQA[0] = 11;   // light-light
  pairQA[4] = 11;   // antilight-antilight
  pairQA[2] = 13;   // light-charm
  pairQA[3] = 13;   // light-anticharm
  pairQA[5] = 13;   // antilight charm
  pairQA[6] = 13;   // antilight anticharm

  AliFemtoDreamCollConfig *config = new AliFemtoDreamCollConfig("Femto", "Femto");
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
  config->SetPtQA((suffix == "0" && fullBlastQA));
  if (!isMCtruth) {
    config->SetPhiEtaBinnign(suffix == "0" && fullBlastQA);
  }
  config->SetMassQA((suffix == "0" && fullBlastQA));
  config->SetPDGCodes(PDGParticles);
  config->SetNBinsHist(NBins);
  config->SetMinKRel(kMin);
  config->SetMaxKRel(kMax);
  config->SetMixingDepth(mixingDepth);
  config->SetUseEventMixing(true);
  config->SetMultiplicityEstimator(AliFemtoDreamEvent::kRef08);
  config->SetMinimalBookingME(suffix != "0");

  AliAnalysisTaskCharmingFemto *task = new AliAnalysisTaskCharmingFemto(
      "AliAnalysisTaskCharmingFemto", isMC, isMCtruth, useTree);
  task->SetLightweight(suffix != "0");
  task->SetEventCuts(evtCuts);
  task->SetProtonCuts(TrackCuts);
  task->SetAntiProtonCuts(AntiTrackCuts);
  task->SetCollectionConfig(config);
  task->SetDecayChannel(channelHF);
  task->SetHFCuts(analysisCutsHF);
  task->SetAODMismatchProtection(useAODProtection);
  task->SetMassSelection(massSelection);
  task->SetUseMCTruthReco(useMCTruthReco);
  task->SetColsToSave(colsToSave);

  if(applyML) {
    task->SetDoMLApplication(applyML);
    task->SetMLConfigFile(configML);
  }

  task->SetLightPDG(pdgBuddy);
  
  if (trigger == "kINT7") {
    task->SelectCollisionCandidates(AliVEvent::kINT7);
    task->SetTrigger(AliVEvent::kINT7);
  } else if (trigger == "kHighMultV0") {
    task->SelectCollisionCandidates(AliVEvent::kHighMultV0);
    task->SetTrigger(AliVEvent::kHighMultV0);
  }

  if (isMCtruth){
    if (std::abs(pdgBuddy) == 321) {
      task->SetBuddypTLowMCTRUTH(0.15);
      task->SetBuddypTHighMCTRUTH(4.0);
    }
    if (std::abs(pdgBuddy) == 211) {
      task->SetBuddypTLowMCTRUTH(0.14);
      task->SetBuddypTHighMCTRUTH(4.0);
    }

    task->SetBuddyEtaMCTRUTH(0.8);
    task->SetBuddyOriginMCTRUTH(0);
    task->SetDmesonOriginMCTRUTH(0);
  }
  // Cutvariation currently not implemented for mctruth gen/reco

  // if(useMCTruthReco){
  //   task->SetBuddyOriginMCTRUTH(0);
  //   task->SetDmesonOriginMCTRUTH(0);
  //   if (suffix == "1") {
  //     task->SetBuddyOriginMCTRUTH(0);
  //     task->SetDmesonOriginMCTRUTH(1);
  //   }
  //   if (suffix == "2") {
  //     task->SetBuddyOriginMCTRUTH(0);
  //     task->SetDmesonOriginMCTRUTH(2);
  //   }
  //   if (suffix == "3") {
  //     task->SetBuddyOriginMCTRUTH(1);
  //     task->SetDmesonOriginMCTRUTH(1);
  //   }
  //   if (suffix == "4") {
  //     task->SetBuddyOriginMCTRUTH(1);
  //     task->SetDmesonOriginMCTRUTH(2);
  //   }
  //   if (suffix == "5") {
  //     task->SetBuddyOriginMCTRUTH(2);
  //     task->SetDmesonOriginMCTRUTH(1);
  //   }
  //   if (suffix == "6") {
  //     task->SetBuddyOriginMCTRUTH(2);
  //     task->SetDmesonOriginMCTRUTH(2);
  //   }
  // }

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

  int nOutput = 9;
  if (isMC) {
    TString TrkCutsMCName = Form("%sTrackCutsMC%s", addon.Data(),
                                 suffix.Data());
    AliAnalysisDataContainer *coutputTrkCutsMC = mgr->CreateContainer(
        TrkCutsMCName.Data(), TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), TrkCutsMCName.Data()));
    mgr->ConnectOutput(task, nOutput++, coutputTrkCutsMC);

    TString AntiTrkCutsMCName = Form("%sAntiTrackCutsMC%s", addon.Data(),
                                     suffix.Data());
    AliAnalysisDataContainer *coutputAntiTrkCutsMC = mgr->CreateContainer(
        AntiTrkCutsMCName.Data(), TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), AntiTrkCutsMCName.Data()));
    mgr->ConnectOutput(task, nOutput++, coutputAntiTrkCutsMC);
  }

  if(useTree) {
    TString treeDir = Form("%s:%sTrees%s", file.Data(), addon.Data(), suffix.Data());

    AliAnalysisDataContainer *coutputTreeSE_pp = mgr->CreateContainer(
        Form("%stSE_pp%s", addon.Data(), suffix.Data()), TTree::Class(),
        AliAnalysisManager::kOutputContainer, treeDir);
    coutputTreeSE_pp->SetSpecialOutput();
    mgr->ConnectOutput(task, nOutput++, coutputTreeSE_pp);

    AliAnalysisDataContainer *coutputTreeSE_mm = mgr->CreateContainer(
        Form("%stSE_mm%s", addon.Data(), suffix.Data()), TTree::Class(),
        AliAnalysisManager::kOutputContainer, treeDir);
    coutputTreeSE_mm->SetSpecialOutput();
    mgr->ConnectOutput(task, nOutput++, coutputTreeSE_mm);

    AliAnalysisDataContainer *coutputTreeSE_mp = mgr->CreateContainer(
        Form("%stSE_mp%s", addon.Data(), suffix.Data()), TTree::Class(),
        AliAnalysisManager::kOutputContainer, treeDir);
    coutputTreeSE_mp->SetSpecialOutput();
    mgr->ConnectOutput(task, nOutput++, coutputTreeSE_mp);

    AliAnalysisDataContainer *coutputTreeSE_pm = mgr->CreateContainer(
        Form("%stSE_pm%s", addon.Data(), suffix.Data()), TTree::Class(),
        AliAnalysisManager::kOutputContainer, treeDir);
    coutputTreeSE_pm->SetSpecialOutput();
    mgr->ConnectOutput(task, nOutput++, coutputTreeSE_pm);

    AliAnalysisDataContainer *coutputTreeME_pp = mgr->CreateContainer(
        Form("%stME_pp%s", addon.Data(), suffix.Data()), TTree::Class(),
        AliAnalysisManager::kOutputContainer, treeDir);
    coutputTreeME_pp->SetSpecialOutput();
    mgr->ConnectOutput(task, nOutput++, coutputTreeME_pp);

    AliAnalysisDataContainer *coutputTreeME_mm = mgr->CreateContainer(
        Form("%stME_mm%s", addon.Data(), suffix.Data()), TTree::Class(),
        AliAnalysisManager::kOutputContainer, treeDir);
    coutputTreeME_mm->SetSpecialOutput();
    mgr->ConnectOutput(task, nOutput++, coutputTreeME_mm);

    AliAnalysisDataContainer *coutputTreeME_mp = mgr->CreateContainer(
        Form("%stME_mp%s", addon.Data(), suffix.Data()), TTree::Class(),
        AliAnalysisManager::kOutputContainer, treeDir);
    coutputTreeME_mp->SetSpecialOutput();
    mgr->ConnectOutput(task, nOutput++, coutputTreeME_mp);

    AliAnalysisDataContainer *coutputTreeME_pm = mgr->CreateContainer(
        Form("%stME_pm%s", addon.Data(), suffix.Data()), TTree::Class(),
        AliAnalysisManager::kOutputContainer, treeDir);
    coutputTreeME_pm->SetSpecialOutput();
    mgr->ConnectOutput(task, nOutput++, coutputTreeME_pm);

  }
  return task;
}
