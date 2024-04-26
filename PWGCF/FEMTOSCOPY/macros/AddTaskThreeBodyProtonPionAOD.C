#if !defined(__CINT__) || defined(__CLING__)
#include <vector>
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskAODThreeBodyProtonPrimary.h"
#include "AliFemtoDreamEventCuts.h"
#include "AliFemtoDreamTrackCuts.h"
#include "AliFemtoDreamCascadeCuts.h"
#include "AliFemtoDreamCollConfig.h"
#endif

AliAnalysisTaskSE *AddTaskThreeBodyProtonPionAOD(
  bool fullBlastQA = true, //1
  bool isMC = false, //2
  TString trigger = "kHM", //3
  bool UseSphericityCut = true, //4
  bool DoThreeBody = true, //5
  bool turnoffClosePairRejectionCompletely = false, //6
  bool ClosePairRejectionForAll = false, //7
  bool RunPlotPhiTheta = false, //8
  double Q3LimitForDeltaPhiDeltaEta = 0.4, //9
  bool StandardMixing = true, //10
  int MixingDepth = 10, //11
  bool RunPlotMult = true, //12
  bool RunPlotInvMass = false, //13
  double Q3MinValue = 0., //14
  double Q3cutValue = 1.0,//15
  bool RunmTPlots = false, //16,
  bool RunPairMultThreeBody = false, //17
  float PionMaxPt = 4.0, //18
  bool RemoveResonances = false, //19
  bool GetProjector = false, //20
  bool GetMomentumResolution = false, //21
  int FilterBitProton = 128, //22
  const char *cutVariation = "0" //23
  ){

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

  if (UseSphericityCut){
    float SpherDown = 0.7;
    evtCuts->SetSphericityCuts(SpherDown, 1.0, 0.5); // THINK IF NEEDED FOR THREE BODY
  }

// Track Cuts for Protons and Antiprotons ================================================================================================
  AliFemtoDreamTrackCuts *TrackCutsProton = AliFemtoDreamTrackCuts::PrimProtonCuts(
      isMC, true, false, false);
  TrackCutsProton->SetFilterBit(FilterBitProton);
  TrackCutsProton->SetCutCharge(1);

  AliFemtoDreamTrackCuts *TrackCutsAntiProton =
      AliFemtoDreamTrackCuts::PrimProtonCuts(isMC, true, false, false);
  TrackCutsAntiProton->SetFilterBit(FilterBitProton);
  TrackCutsAntiProton->SetCutCharge(-1);

  // Track Cuts for Pions  =================================================================================================================

  AliFemtoDreamTrackCuts *TrackCutsPion = NULL;
  AliFemtoDreamTrackCuts *TrackCutsAntiPion = NULL;

  TrackCutsPion = AliFemtoDreamTrackCuts::PrimPionCuts(isMC, true, false, false);
  TrackCutsPion->SetPtRange(0.14, PionMaxPt);
  TrackCutsPion->SetFilterBit(96);
  TrackCutsPion->SetCutCharge(1);

  TrackCutsAntiPion = AliFemtoDreamTrackCuts::PrimPionCuts(isMC, true, false, false);
  TrackCutsAntiPion->SetPtRange(0.14, PionMaxPt);
  TrackCutsAntiPion->SetFilterBit(96);
  TrackCutsAntiPion->SetCutCharge(-1);

  if (!fullBlastQA) {
    evtCuts->SetMinimalBooking(true);
    TrackCutsProton->SetMinimalBooking(true);
    TrackCutsAntiProton->SetMinimalBooking(true);
    TrackCutsPion->SetMinimalBooking(true);
    TrackCutsAntiPion->SetMinimalBooking(true);
  }

  AliFemtoDreamCollConfig *config = new AliFemtoDreamCollConfig("Femto",
                                                                "Femto", false);
  // Femto Collection
  std::vector<int> PDGParticles;
  PDGParticles.push_back(2212); //Proton
  PDGParticles.push_back(2212);
  PDGParticles.push_back(211); //Pion
  PDGParticles.push_back(211);

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
  config->SetDeltaEtaMax(0.0);
  config->SetDeltaPhiMax(0.0);
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

  double DeltaPhiPP = 0.017; //default 
  double DeltaEtaPP = 0.017; 

  double DeltaPhiPPrim = 0.04; 
  double DeltaEtaPPrim = 0.017; 

 Float_t Proton_pT_VarLow = 0.4;
 Float_t Proton_pT_VarHigh = 0.6;
 Float_t Proton_Eta_VarLow = 0.77;
 Float_t Proton_Eta_VarHigh = 0.85;
 Float_t Proton_Clusters_VarLow = 70.;
 Float_t Proton_Clusters_VarHigh = 90.;
 Float_t Proton_Sigma_VarLow = 2.5;
 Float_t Proton_Sigma_VarHigh = 3.5;

 Float_t Pion_pT_VarLow = 0.12;
 Float_t Pion_pT_VarHigh = 0.15;
 Float_t Pion_Eta_VarLow = 0.77;
 Float_t Pion_Eta_VarHigh = 0.85;
 //Float_t Pion_Eta_VarLow = 0.7; //old 
 //Float_t Pion_Eta_VarHigh = 0.9; //old
 Float_t Pion_Clusters_VarLow = 70;
 Float_t Pion_Clusters_VarHigh = 90;
 Float_t Pion_Sigma_VarLow = 2.5;
 Float_t Pion_Sigma_VarHigh = 3.5;

 //Float_t Pion_Sigma_VarLow = 2.7; //old
 //Float_t Pion_Sigma_VarHigh = 3.3; //old

 Float_t DPhi_VarLow = 0.035;
 Float_t DPhi_VarHigh = 0.045;
 Float_t DEta_VarLow = 0.015;
 Float_t DEta_VarHigh = 0.019;

 Float_t DPhi_PP_Default = 0.017;
 Float_t DPhi_PP_VarHigh = 0.019;
 Float_t DEta_PP_Default = 0.017;
 Float_t DEta_PP_VarHigh = 0.019;

 if (suffix == "1"){
 
    TrackCutsProton->SetNClsTPC(Proton_Clusters_VarLow);
    TrackCutsAntiProton->SetNClsTPC(Proton_Clusters_VarLow);
 
    TrackCutsPion->SetPID(AliPID::kPion, 0.5,Pion_Sigma_VarHigh);
    TrackCutsAntiPion->SetPID(AliPID::kPion, 0.5,Pion_Sigma_VarHigh);
 
    DeltaEtaPPrim = DEta_VarLow;
 
    DeltaPhiPP = DEta_VarLow;
    DeltaEtaPP = DPhi_PP_Default;
 
  } else if (suffix == "2"){
 
    TrackCutsProton->SetEtaRange(-Proton_Eta_VarHigh,Proton_Eta_VarHigh);
    TrackCutsAntiProton->SetEtaRange(-Proton_Eta_VarHigh,Proton_Eta_VarHigh);
    TrackCutsProton->SetPtRange(Proton_pT_VarHigh,4.05);
    TrackCutsAntiProton->SetPtRange(Proton_pT_VarHigh,4.05);
    TrackCutsProton->SetNClsTPC(Proton_Clusters_VarHigh);
    TrackCutsAntiProton->SetNClsTPC(Proton_Clusters_VarHigh);
 
    TrackCutsPion->SetPID(AliPID::kPion, 0.5,Pion_Sigma_VarHigh);
    TrackCutsAntiPion->SetPID(AliPID::kPion, 0.5,Pion_Sigma_VarHigh);
 
    DeltaPhiPPrim = DPhi_VarLow;
    DeltaEtaPPrim = DEta_VarHigh;
 
    DeltaPhiPP = DPhi_PP_VarHigh;
 
  } else if (suffix == "3"){
 
    TrackCutsProton->SetPtRange(Proton_pT_VarHigh,4.05);
    TrackCutsAntiProton->SetPtRange(Proton_pT_VarHigh,4.05);
    TrackCutsProton->SetNClsTPC(Proton_Clusters_VarLow);
    TrackCutsAntiProton->SetNClsTPC(Proton_Clusters_VarLow);
 
    TrackCutsPion->SetPID(AliPID::kPion, 0.5,Pion_Sigma_VarHigh);
    TrackCutsAntiPion->SetPID(AliPID::kPion, 0.5,Pion_Sigma_VarHigh);
    TrackCutsPion->SetEtaRange(-Pion_Eta_VarLow,Pion_Eta_VarLow);
    TrackCutsAntiPion->SetEtaRange(-Pion_Eta_VarLow,Pion_Eta_VarLow);
    TrackCutsPion->SetPtRange(Pion_pT_VarHigh,4.0);
    TrackCutsAntiPion->SetPtRange(Pion_pT_VarHigh,4.0);
 
    DeltaPhiPPrim = DPhi_VarHigh;
    DeltaEtaPPrim = DEta_VarHigh;
 
    DeltaPhiPP = DEta_VarLow;
    DeltaEtaPP = DPhi_PP_Default;
 
  } else if (suffix == "4"){
 
    TrackCutsProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarHigh);
    TrackCutsAntiProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarHigh);
    TrackCutsProton->SetEtaRange(-Proton_Eta_VarHigh,Proton_Eta_VarHigh);
    TrackCutsAntiProton->SetEtaRange(-Proton_Eta_VarHigh,Proton_Eta_VarHigh);
    TrackCutsProton->SetPtRange(Proton_pT_VarHigh,4.05);
    TrackCutsAntiProton->SetPtRange(Proton_pT_VarHigh,4.05);
    TrackCutsProton->SetNClsTPC(Proton_Clusters_VarLow);
    TrackCutsAntiProton->SetNClsTPC(Proton_Clusters_VarLow);
 
    TrackCutsPion->SetPID(AliPID::kPion, 0.5,Pion_Sigma_VarHigh);
    TrackCutsAntiPion->SetPID(AliPID::kPion, 0.5,Pion_Sigma_VarHigh);
    TrackCutsPion->SetEtaRange(-Pion_Eta_VarHigh,Pion_Eta_VarHigh);
    TrackCutsAntiPion->SetEtaRange(-Pion_Eta_VarHigh,Pion_Eta_VarHigh);
    TrackCutsPion->SetPtRange(Pion_pT_VarLow,4.0);
    TrackCutsAntiPion->SetPtRange(Pion_pT_VarLow,4.0);
    TrackCutsPion->SetNClsTPC(Pion_Clusters_VarHigh);
    TrackCutsAntiPion->SetNClsTPC(Pion_Clusters_VarHigh);
 
 
    DeltaPhiPP = DPhi_PP_VarHigh;
    DeltaEtaPP = DPhi_PP_Default;
 
  } else if (suffix == "5"){
 
 
    TrackCutsPion->SetPID(AliPID::kPion, 0.5,Pion_Sigma_VarLow);
    TrackCutsAntiPion->SetPID(AliPID::kPion, 0.5,Pion_Sigma_VarLow);
    TrackCutsPion->SetPtRange(Pion_pT_VarLow,4.0);
    TrackCutsAntiPion->SetPtRange(Pion_pT_VarLow,4.0);
 
    DeltaPhiPPrim = DPhi_VarHigh;
    DeltaEtaPPrim = DEta_VarHigh;
 
    DeltaPhiPP = DPhi_PP_VarHigh;
    DeltaEtaPP = DEta_PP_VarHigh;
 
  } else if (suffix == "6"){
 
    TrackCutsProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarHigh);
    TrackCutsAntiProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarHigh);
    TrackCutsProton->SetEtaRange(-Proton_Eta_VarLow,Proton_Eta_VarLow);
    TrackCutsAntiProton->SetEtaRange(-Proton_Eta_VarLow,Proton_Eta_VarLow);
    TrackCutsProton->SetPtRange(Proton_pT_VarLow,4.05);
    TrackCutsAntiProton->SetPtRange(Proton_pT_VarLow,4.05);
 
    TrackCutsPion->SetPID(AliPID::kPion, 0.5,Pion_Sigma_VarLow);
    TrackCutsAntiPion->SetPID(AliPID::kPion, 0.5,Pion_Sigma_VarLow);
    TrackCutsPion->SetPtRange(Pion_pT_VarLow,4.0);
    TrackCutsAntiPion->SetPtRange(Pion_pT_VarLow,4.0);
    TrackCutsPion->SetNClsTPC(Pion_Clusters_VarLow);
    TrackCutsAntiPion->SetNClsTPC(Pion_Clusters_VarLow);
 
    DeltaPhiPPrim = DPhi_VarLow;
    DeltaEtaPPrim = DEta_VarLow;
 
    DeltaPhiPP = DEta_VarLow;
 
  } else if (suffix == "7"){
 
    TrackCutsProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarHigh);
    TrackCutsAntiProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarHigh);
    TrackCutsProton->SetEtaRange(-Proton_Eta_VarHigh,Proton_Eta_VarHigh);
    TrackCutsAntiProton->SetEtaRange(-Proton_Eta_VarHigh,Proton_Eta_VarHigh);
    TrackCutsProton->SetPtRange(Proton_pT_VarLow,4.05);
    TrackCutsAntiProton->SetPtRange(Proton_pT_VarLow,4.05);
    TrackCutsProton->SetNClsTPC(Proton_Clusters_VarHigh);
    TrackCutsAntiProton->SetNClsTPC(Proton_Clusters_VarHigh);
 
    TrackCutsPion->SetPID(AliPID::kPion, 0.5,Pion_Sigma_VarLow);
    TrackCutsAntiPion->SetPID(AliPID::kPion, 0.5,Pion_Sigma_VarLow);
    TrackCutsPion->SetEtaRange(-Pion_Eta_VarHigh,Pion_Eta_VarHigh);
    TrackCutsAntiPion->SetEtaRange(-Pion_Eta_VarHigh,Pion_Eta_VarHigh);
    TrackCutsPion->SetNClsTPC(Pion_Clusters_VarHigh);
    TrackCutsAntiPion->SetNClsTPC(Pion_Clusters_VarHigh);
 
    DeltaEtaPPrim = DEta_VarHigh;
 
 
  } else if (suffix == "8"){
 
    TrackCutsProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarLow);
    TrackCutsAntiProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarLow);
    TrackCutsProton->SetPtRange(Proton_pT_VarLow,4.05);
    TrackCutsAntiProton->SetPtRange(Proton_pT_VarLow,4.05);
    TrackCutsProton->SetNClsTPC(Proton_Clusters_VarLow);
    TrackCutsAntiProton->SetNClsTPC(Proton_Clusters_VarLow);
 
    TrackCutsPion->SetPID(AliPID::kPion, 0.5,Pion_Sigma_VarHigh);
    TrackCutsAntiPion->SetPID(AliPID::kPion, 0.5,Pion_Sigma_VarHigh);
    TrackCutsPion->SetEtaRange(-Pion_Eta_VarHigh,Pion_Eta_VarHigh);
    TrackCutsAntiPion->SetEtaRange(-Pion_Eta_VarHigh,Pion_Eta_VarHigh);
 
    DeltaPhiPPrim = DPhi_VarLow;
 
    DeltaEtaPP = DPhi_PP_Default;
 
  } else if (suffix == "9"){
 
    TrackCutsProton->SetEtaRange(-Proton_Eta_VarHigh,Proton_Eta_VarHigh);
    TrackCutsAntiProton->SetEtaRange(-Proton_Eta_VarHigh,Proton_Eta_VarHigh);
    TrackCutsProton->SetNClsTPC(Proton_Clusters_VarLow);
    TrackCutsAntiProton->SetNClsTPC(Proton_Clusters_VarLow);
 
    TrackCutsPion->SetEtaRange(-Pion_Eta_VarHigh,Pion_Eta_VarHigh);
    TrackCutsAntiPion->SetEtaRange(-Pion_Eta_VarHigh,Pion_Eta_VarHigh);
    TrackCutsPion->SetNClsTPC(Pion_Clusters_VarLow);
    TrackCutsAntiPion->SetNClsTPC(Pion_Clusters_VarLow);
 
    DeltaPhiPPrim = DPhi_VarHigh;
    DeltaEtaPPrim = DEta_VarHigh;
 
    DeltaPhiPP = DEta_VarLow;
 
  } else if (suffix == "10"){
 
    TrackCutsProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarHigh);
    TrackCutsAntiProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarHigh);
    TrackCutsProton->SetEtaRange(-Proton_Eta_VarLow,Proton_Eta_VarLow);
    TrackCutsAntiProton->SetEtaRange(-Proton_Eta_VarLow,Proton_Eta_VarLow);
    TrackCutsProton->SetNClsTPC(Proton_Clusters_VarHigh);
    TrackCutsAntiProton->SetNClsTPC(Proton_Clusters_VarHigh);
 
    TrackCutsPion->SetPID(AliPID::kPion, 0.5,Pion_Sigma_VarHigh);
    TrackCutsAntiPion->SetPID(AliPID::kPion, 0.5,Pion_Sigma_VarHigh);
    TrackCutsPion->SetEtaRange(-Pion_Eta_VarHigh,Pion_Eta_VarHigh);
    TrackCutsAntiPion->SetEtaRange(-Pion_Eta_VarHigh,Pion_Eta_VarHigh);
    TrackCutsPion->SetNClsTPC(Pion_Clusters_VarHigh);
    TrackCutsAntiPion->SetNClsTPC(Pion_Clusters_VarHigh);
 
    DeltaEtaPPrim = DEta_VarHigh;
 
    DeltaPhiPP = DEta_VarLow;
 
  } else if (suffix == "11"){
 
    TrackCutsProton->SetEtaRange(-Proton_Eta_VarLow,Proton_Eta_VarLow);
    TrackCutsAntiProton->SetEtaRange(-Proton_Eta_VarLow,Proton_Eta_VarLow);
    TrackCutsProton->SetPtRange(Proton_pT_VarHigh,4.05);
    TrackCutsAntiProton->SetPtRange(Proton_pT_VarHigh,4.05);
 
    TrackCutsPion->SetPtRange(Pion_pT_VarLow,4.0);
    TrackCutsAntiPion->SetPtRange(Pion_pT_VarLow,4.0);
    TrackCutsPion->SetNClsTPC(Pion_Clusters_VarHigh);
    TrackCutsAntiPion->SetNClsTPC(Pion_Clusters_VarHigh);
 
    DeltaPhiPPrim = DPhi_VarLow;
 
    DeltaPhiPP = DPhi_PP_VarHigh;
    DeltaEtaPP = DPhi_PP_Default;
 
  } else if (suffix == "12"){
 
    TrackCutsProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarHigh);
    TrackCutsAntiProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarHigh);
    TrackCutsProton->SetEtaRange(-Proton_Eta_VarLow,Proton_Eta_VarLow);
    TrackCutsAntiProton->SetEtaRange(-Proton_Eta_VarLow,Proton_Eta_VarLow);
    TrackCutsProton->SetPtRange(Proton_pT_VarLow,4.05);
    TrackCutsAntiProton->SetPtRange(Proton_pT_VarLow,4.05);
 
    TrackCutsPion->SetPID(AliPID::kPion, 0.5,Pion_Sigma_VarLow);
    TrackCutsAntiPion->SetPID(AliPID::kPion, 0.5,Pion_Sigma_VarLow);
    TrackCutsPion->SetEtaRange(-Pion_Eta_VarLow,Pion_Eta_VarLow);
    TrackCutsAntiPion->SetEtaRange(-Pion_Eta_VarLow,Pion_Eta_VarLow);
    TrackCutsPion->SetPtRange(Pion_pT_VarHigh,4.0);
    TrackCutsAntiPion->SetPtRange(Pion_pT_VarHigh,4.0);
    TrackCutsPion->SetNClsTPC(Pion_Clusters_VarHigh);
    TrackCutsAntiPion->SetNClsTPC(Pion_Clusters_VarHigh);
 
    DeltaEtaPPrim = DEta_VarLow;
 
    DeltaPhiPP = DEta_VarLow;
 
  } else if (suffix == "13"){
 
    TrackCutsProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarLow);
    TrackCutsAntiProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarLow);
    TrackCutsProton->SetEtaRange(-Proton_Eta_VarHigh,Proton_Eta_VarHigh);
    TrackCutsAntiProton->SetEtaRange(-Proton_Eta_VarHigh,Proton_Eta_VarHigh);
    TrackCutsProton->SetNClsTPC(Proton_Clusters_VarHigh);
    TrackCutsAntiProton->SetNClsTPC(Proton_Clusters_VarHigh);
 
    TrackCutsPion->SetPID(AliPID::kPion, 0.5,Pion_Sigma_VarHigh);
    TrackCutsAntiPion->SetPID(AliPID::kPion, 0.5,Pion_Sigma_VarHigh);
    TrackCutsPion->SetNClsTPC(Pion_Clusters_VarLow);
    TrackCutsAntiPion->SetNClsTPC(Pion_Clusters_VarLow);
 
    DeltaPhiPPrim = DPhi_VarLow;
 
    DeltaPhiPP = DPhi_PP_VarHigh;
 
  } else if (suffix == "14"){
 
    TrackCutsProton->SetEtaRange(-Proton_Eta_VarHigh,Proton_Eta_VarHigh);
    TrackCutsAntiProton->SetEtaRange(-Proton_Eta_VarHigh,Proton_Eta_VarHigh);
    TrackCutsProton->SetPtRange(Proton_pT_VarHigh,4.05);
    TrackCutsAntiProton->SetPtRange(Proton_pT_VarHigh,4.05);
    TrackCutsProton->SetNClsTPC(Proton_Clusters_VarLow);
    TrackCutsAntiProton->SetNClsTPC(Proton_Clusters_VarLow);
 
    TrackCutsPion->SetEtaRange(-Pion_Eta_VarHigh,Pion_Eta_VarHigh);
    TrackCutsAntiPion->SetEtaRange(-Pion_Eta_VarHigh,Pion_Eta_VarHigh);
    TrackCutsPion->SetPtRange(Pion_pT_VarHigh,4.0);
    TrackCutsAntiPion->SetPtRange(Pion_pT_VarHigh,4.0);
    TrackCutsPion->SetNClsTPC(Pion_Clusters_VarLow);
    TrackCutsAntiPion->SetNClsTPC(Pion_Clusters_VarLow);
 
    DeltaPhiPPrim = DPhi_VarLow;
    DeltaEtaPPrim = DEta_VarHigh;
 
 
  } else if (suffix == "15"){
 
    TrackCutsProton->SetEtaRange(-Proton_Eta_VarHigh,Proton_Eta_VarHigh);
    TrackCutsAntiProton->SetEtaRange(-Proton_Eta_VarHigh,Proton_Eta_VarHigh);
    TrackCutsProton->SetPtRange(Proton_pT_VarHigh,4.05);
    TrackCutsAntiProton->SetPtRange(Proton_pT_VarHigh,4.05);
    TrackCutsProton->SetNClsTPC(Proton_Clusters_VarHigh);
    TrackCutsAntiProton->SetNClsTPC(Proton_Clusters_VarHigh);
 
    TrackCutsPion->SetPID(AliPID::kPion, 0.5,Pion_Sigma_VarLow);
    TrackCutsAntiPion->SetPID(AliPID::kPion, 0.5,Pion_Sigma_VarLow);
    TrackCutsPion->SetEtaRange(-Pion_Eta_VarHigh,Pion_Eta_VarHigh);
    TrackCutsAntiPion->SetEtaRange(-Pion_Eta_VarHigh,Pion_Eta_VarHigh);
    TrackCutsPion->SetPtRange(Pion_pT_VarLow,4.0);
    TrackCutsAntiPion->SetPtRange(Pion_pT_VarLow,4.0);
 
    DeltaEtaPPrim = DEta_VarHigh;
 
    DeltaPhiPP = DEta_VarLow;
    DeltaEtaPP = DPhi_PP_Default;
 
  } else if (suffix == "16"){
 
    TrackCutsProton->SetEtaRange(-Proton_Eta_VarHigh,Proton_Eta_VarHigh);
    TrackCutsAntiProton->SetEtaRange(-Proton_Eta_VarHigh,Proton_Eta_VarHigh);
    TrackCutsProton->SetPtRange(Proton_pT_VarHigh,4.05);
    TrackCutsAntiProton->SetPtRange(Proton_pT_VarHigh,4.05);
    TrackCutsProton->SetNClsTPC(Proton_Clusters_VarHigh);
    TrackCutsAntiProton->SetNClsTPC(Proton_Clusters_VarHigh);
 
    TrackCutsPion->SetPtRange(Pion_pT_VarLow,4.0);
    TrackCutsAntiPion->SetPtRange(Pion_pT_VarLow,4.0);
    TrackCutsPion->SetNClsTPC(Pion_Clusters_VarHigh);
    TrackCutsAntiPion->SetNClsTPC(Pion_Clusters_VarHigh);
 
 
    DeltaPhiPP = DEta_VarLow;
    DeltaEtaPP = DPhi_PP_Default;
 
  } else if (suffix == "17"){
 
    TrackCutsProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarHigh);
    TrackCutsAntiProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarHigh);
 
    TrackCutsPion->SetPID(AliPID::kPion, 0.5,Pion_Sigma_VarHigh);
    TrackCutsAntiPion->SetPID(AliPID::kPion, 0.5,Pion_Sigma_VarHigh);
    TrackCutsPion->SetPtRange(Pion_pT_VarLow,4.0);
    TrackCutsAntiPion->SetPtRange(Pion_pT_VarLow,4.0);
    TrackCutsPion->SetNClsTPC(Pion_Clusters_VarLow);
    TrackCutsAntiPion->SetNClsTPC(Pion_Clusters_VarLow);
 
    DeltaEtaPPrim = DEta_VarHigh;
 
    DeltaEtaPP = DPhi_PP_Default;
 
  } else if (suffix == "18"){
 
    TrackCutsProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarHigh);
    TrackCutsAntiProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarHigh);
    TrackCutsProton->SetEtaRange(-Proton_Eta_VarLow,Proton_Eta_VarLow);
    TrackCutsAntiProton->SetEtaRange(-Proton_Eta_VarLow,Proton_Eta_VarLow);
    TrackCutsProton->SetNClsTPC(Proton_Clusters_VarHigh);
    TrackCutsAntiProton->SetNClsTPC(Proton_Clusters_VarHigh);
 
    TrackCutsPion->SetEtaRange(-Pion_Eta_VarHigh,Pion_Eta_VarHigh);
    TrackCutsAntiPion->SetEtaRange(-Pion_Eta_VarHigh,Pion_Eta_VarHigh);
    TrackCutsPion->SetPtRange(Pion_pT_VarLow,4.0);
    TrackCutsAntiPion->SetPtRange(Pion_pT_VarLow,4.0);
 
    DeltaEtaPPrim = DEta_VarHigh;
 
    DeltaPhiPP = DPhi_PP_VarHigh;
 
  } else if (suffix == "19"){
 
    TrackCutsProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarLow);
    TrackCutsAntiProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarLow);
    TrackCutsProton->SetNClsTPC(Proton_Clusters_VarHigh);
    TrackCutsAntiProton->SetNClsTPC(Proton_Clusters_VarHigh);
 
    TrackCutsPion->SetPID(AliPID::kPion, 0.5,Pion_Sigma_VarHigh);
    TrackCutsAntiPion->SetPID(AliPID::kPion, 0.5,Pion_Sigma_VarHigh);
    TrackCutsPion->SetEtaRange(-Pion_Eta_VarLow,Pion_Eta_VarLow);
    TrackCutsAntiPion->SetEtaRange(-Pion_Eta_VarLow,Pion_Eta_VarLow);
    TrackCutsPion->SetNClsTPC(Pion_Clusters_VarHigh);
    TrackCutsAntiPion->SetNClsTPC(Pion_Clusters_VarHigh);
 
    DeltaPhiPPrim = DPhi_VarHigh;
    DeltaEtaPPrim = DEta_VarHigh;
 
    DeltaPhiPP = DPhi_PP_VarHigh;
    DeltaEtaPP = DEta_PP_VarHigh;
 
  } else if (suffix == "20"){
 
    TrackCutsProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarHigh);
    TrackCutsAntiProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarHigh);
    TrackCutsProton->SetEtaRange(-Proton_Eta_VarHigh,Proton_Eta_VarHigh);
    TrackCutsAntiProton->SetEtaRange(-Proton_Eta_VarHigh,Proton_Eta_VarHigh);
    TrackCutsProton->SetPtRange(Proton_pT_VarHigh,4.05);
    TrackCutsAntiProton->SetPtRange(Proton_pT_VarHigh,4.05);
 
    TrackCutsPion->SetPID(AliPID::kPion, 0.5,Pion_Sigma_VarLow);
    TrackCutsAntiPion->SetPID(AliPID::kPion, 0.5,Pion_Sigma_VarLow);
    TrackCutsPion->SetPtRange(Pion_pT_VarHigh,4.0);
    TrackCutsAntiPion->SetPtRange(Pion_pT_VarHigh,4.0);
    TrackCutsPion->SetNClsTPC(Pion_Clusters_VarLow);
    TrackCutsAntiPion->SetNClsTPC(Pion_Clusters_VarLow);
 
    DeltaPhiPPrim = DPhi_VarHigh;
    DeltaEtaPPrim = DEta_VarHigh;
 
    DeltaPhiPP = DPhi_PP_VarHigh;
    DeltaEtaPP = DPhi_PP_Default;
 
  } else if (suffix == "21"){
 
    TrackCutsProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarLow);
    TrackCutsAntiProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarLow);
    TrackCutsProton->SetNClsTPC(Proton_Clusters_VarHigh);
    TrackCutsAntiProton->SetNClsTPC(Proton_Clusters_VarHigh);
 
    TrackCutsPion->SetPID(AliPID::kPion, 0.5,Pion_Sigma_VarLow);
    TrackCutsAntiPion->SetPID(AliPID::kPion, 0.5,Pion_Sigma_VarLow);
    TrackCutsPion->SetPtRange(Pion_pT_VarHigh,4.0);
    TrackCutsAntiPion->SetPtRange(Pion_pT_VarHigh,4.0);
 
    DeltaEtaPPrim = DEta_VarLow;
 
    DeltaPhiPP = DPhi_PP_VarHigh;
    DeltaEtaPP = DEta_PP_VarHigh;
 
  } else if (suffix == "22"){
 
    TrackCutsProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarHigh);
    TrackCutsAntiProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarHigh);
    TrackCutsProton->SetEtaRange(-Proton_Eta_VarHigh,Proton_Eta_VarHigh);
    TrackCutsAntiProton->SetEtaRange(-Proton_Eta_VarHigh,Proton_Eta_VarHigh);
    TrackCutsProton->SetPtRange(Proton_pT_VarHigh,4.05);
    TrackCutsAntiProton->SetPtRange(Proton_pT_VarHigh,4.05);
 
    TrackCutsPion->SetPID(AliPID::kPion, 0.5,Pion_Sigma_VarLow);
    TrackCutsAntiPion->SetPID(AliPID::kPion, 0.5,Pion_Sigma_VarLow);
    TrackCutsPion->SetPtRange(Pion_pT_VarLow,4.0);
    TrackCutsAntiPion->SetPtRange(Pion_pT_VarLow,4.0);
 
    DeltaEtaPPrim = DEta_VarLow;
 
 
  } else if (suffix == "23"){
 
    TrackCutsProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarLow);
    TrackCutsAntiProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarLow);
    TrackCutsProton->SetEtaRange(-Proton_Eta_VarHigh,Proton_Eta_VarHigh);
    TrackCutsAntiProton->SetEtaRange(-Proton_Eta_VarHigh,Proton_Eta_VarHigh);
    TrackCutsProton->SetPtRange(Proton_pT_VarHigh,4.05);
    TrackCutsAntiProton->SetPtRange(Proton_pT_VarHigh,4.05);
    TrackCutsProton->SetNClsTPC(Proton_Clusters_VarHigh);
    TrackCutsAntiProton->SetNClsTPC(Proton_Clusters_VarHigh);
 
    TrackCutsPion->SetPID(AliPID::kPion, 0.5,Pion_Sigma_VarLow);
    TrackCutsAntiPion->SetPID(AliPID::kPion, 0.5,Pion_Sigma_VarLow);
    TrackCutsPion->SetEtaRange(-Pion_Eta_VarLow,Pion_Eta_VarLow);
    TrackCutsAntiPion->SetEtaRange(-Pion_Eta_VarLow,Pion_Eta_VarLow);
    TrackCutsPion->SetPtRange(Pion_pT_VarHigh,4.0);
    TrackCutsAntiPion->SetPtRange(Pion_pT_VarHigh,4.0);
    TrackCutsPion->SetNClsTPC(Pion_Clusters_VarLow);
    TrackCutsAntiPion->SetNClsTPC(Pion_Clusters_VarLow);
 
    DeltaPhiPPrim = DPhi_VarLow;
    DeltaEtaPPrim = DEta_VarHigh;
 
    DeltaPhiPP = DEta_VarLow;
    DeltaEtaPP = DEta_PP_VarHigh;
 
  } else if (suffix == "24"){
 
    TrackCutsProton->SetNClsTPC(Proton_Clusters_VarLow);
    TrackCutsAntiProton->SetNClsTPC(Proton_Clusters_VarLow);
 
    TrackCutsPion->SetPID(AliPID::kPion, 0.5,Pion_Sigma_VarLow);
    TrackCutsAntiPion->SetPID(AliPID::kPion, 0.5,Pion_Sigma_VarLow);
    TrackCutsPion->SetEtaRange(-Pion_Eta_VarHigh,Pion_Eta_VarHigh);
    TrackCutsAntiPion->SetEtaRange(-Pion_Eta_VarHigh,Pion_Eta_VarHigh);
    TrackCutsPion->SetPtRange(Pion_pT_VarLow,4.0);
    TrackCutsAntiPion->SetPtRange(Pion_pT_VarLow,4.0);
    TrackCutsPion->SetNClsTPC(Pion_Clusters_VarHigh);
    TrackCutsAntiPion->SetNClsTPC(Pion_Clusters_VarHigh);
 
    DeltaPhiPPrim = DPhi_VarHigh;
    DeltaEtaPPrim = DEta_VarHigh;
 
    DeltaPhiPP = DPhi_PP_VarHigh;
    DeltaEtaPP = DPhi_PP_Default;
 
  } else if (suffix == "25"){
 
    TrackCutsProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarLow);
    TrackCutsAntiProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarLow);
    TrackCutsProton->SetEtaRange(-Proton_Eta_VarHigh,Proton_Eta_VarHigh);
    TrackCutsAntiProton->SetEtaRange(-Proton_Eta_VarHigh,Proton_Eta_VarHigh);
    TrackCutsProton->SetPtRange(Proton_pT_VarHigh,4.05);
    TrackCutsAntiProton->SetPtRange(Proton_pT_VarHigh,4.05);
    TrackCutsProton->SetNClsTPC(Proton_Clusters_VarLow);
    TrackCutsAntiProton->SetNClsTPC(Proton_Clusters_VarLow);
 
    TrackCutsPion->SetPID(AliPID::kPion, 0.5,Pion_Sigma_VarHigh);
    TrackCutsAntiPion->SetPID(AliPID::kPion, 0.5,Pion_Sigma_VarHigh);
    TrackCutsPion->SetPtRange(Pion_pT_VarLow,4.0);
    TrackCutsAntiPion->SetPtRange(Pion_pT_VarLow,4.0);
 
    DeltaEtaPPrim = DEta_VarLow;
 
    DeltaPhiPP = DEta_VarLow;
 
  } else if (suffix == "26"){
 
    TrackCutsProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarHigh);
    TrackCutsAntiProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarHigh);
    TrackCutsProton->SetEtaRange(-Proton_Eta_VarHigh,Proton_Eta_VarHigh);
    TrackCutsAntiProton->SetEtaRange(-Proton_Eta_VarHigh,Proton_Eta_VarHigh);
    TrackCutsProton->SetPtRange(Proton_pT_VarLow,4.05);
    TrackCutsAntiProton->SetPtRange(Proton_pT_VarLow,4.05);
    TrackCutsProton->SetNClsTPC(Proton_Clusters_VarHigh);
    TrackCutsAntiProton->SetNClsTPC(Proton_Clusters_VarHigh);
 
    TrackCutsPion->SetNClsTPC(Pion_Clusters_VarLow);
    TrackCutsAntiPion->SetNClsTPC(Pion_Clusters_VarLow);
 
    DeltaEtaPPrim = DEta_VarLow;
 
    DeltaPhiPP = DPhi_PP_VarHigh;
    DeltaEtaPP = DEta_PP_VarHigh;
 
  } else if (suffix == "27"){
 
    TrackCutsProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarLow);
    TrackCutsAntiProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarLow);
    TrackCutsProton->SetNClsTPC(Proton_Clusters_VarLow);
    TrackCutsAntiProton->SetNClsTPC(Proton_Clusters_VarLow);
 
    TrackCutsPion->SetPID(AliPID::kPion, 0.5,Pion_Sigma_VarHigh);
    TrackCutsAntiPion->SetPID(AliPID::kPion, 0.5,Pion_Sigma_VarHigh);
    TrackCutsPion->SetEtaRange(-Pion_Eta_VarLow,Pion_Eta_VarLow);
    TrackCutsAntiPion->SetEtaRange(-Pion_Eta_VarLow,Pion_Eta_VarLow);
 
    DeltaPhiPPrim = DPhi_VarHigh;
    DeltaEtaPPrim = DEta_VarLow;
 
    DeltaPhiPP = DEta_VarLow;
    DeltaEtaPP = DPhi_PP_Default;
 
  } else if (suffix == "28"){
 
    TrackCutsProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarLow);
    TrackCutsAntiProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarLow);
    TrackCutsProton->SetPtRange(Proton_pT_VarLow,4.05);
    TrackCutsAntiProton->SetPtRange(Proton_pT_VarLow,4.05);
    TrackCutsProton->SetNClsTPC(Proton_Clusters_VarLow);
    TrackCutsAntiProton->SetNClsTPC(Proton_Clusters_VarLow);
 
    TrackCutsPion->SetEtaRange(-Pion_Eta_VarLow,Pion_Eta_VarLow);
    TrackCutsAntiPion->SetEtaRange(-Pion_Eta_VarLow,Pion_Eta_VarLow);
    TrackCutsPion->SetNClsTPC(Pion_Clusters_VarHigh);
    TrackCutsAntiPion->SetNClsTPC(Pion_Clusters_VarHigh);
 
    DeltaPhiPPrim = DPhi_VarLow;
    DeltaEtaPPrim = DEta_VarLow;
 
    DeltaPhiPP = DPhi_PP_VarHigh;
    DeltaEtaPP = DEta_PP_VarHigh;
 
  } else if (suffix == "29"){
 
    TrackCutsProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarHigh);
    TrackCutsAntiProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarHigh);
    TrackCutsProton->SetPtRange(Proton_pT_VarLow,4.05);
    TrackCutsAntiProton->SetPtRange(Proton_pT_VarLow,4.05);
    TrackCutsProton->SetNClsTPC(Proton_Clusters_VarHigh);
    TrackCutsAntiProton->SetNClsTPC(Proton_Clusters_VarHigh);
 
    TrackCutsPion->SetPtRange(Pion_pT_VarHigh,4.0);
    TrackCutsAntiPion->SetPtRange(Pion_pT_VarHigh,4.0);
    TrackCutsPion->SetNClsTPC(Pion_Clusters_VarHigh);
    TrackCutsAntiPion->SetNClsTPC(Pion_Clusters_VarHigh);
 
    DeltaPhiPPrim = DPhi_VarLow;
    DeltaEtaPPrim = DEta_VarHigh;
 
    DeltaPhiPP = DEta_VarLow;
 
  } else if (suffix == "30"){
 
    TrackCutsProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarHigh);
    TrackCutsAntiProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarHigh);
    TrackCutsProton->SetEtaRange(-Proton_Eta_VarLow,Proton_Eta_VarLow);
    TrackCutsAntiProton->SetEtaRange(-Proton_Eta_VarLow,Proton_Eta_VarLow);
    TrackCutsProton->SetNClsTPC(Proton_Clusters_VarLow);
    TrackCutsAntiProton->SetNClsTPC(Proton_Clusters_VarLow);
 
    TrackCutsPion->SetNClsTPC(Pion_Clusters_VarLow);
    TrackCutsAntiPion->SetNClsTPC(Pion_Clusters_VarLow);
 
    DeltaPhiPPrim = DPhi_VarLow;
    DeltaEtaPPrim = DEta_VarHigh;
 
    DeltaPhiPP = DPhi_PP_VarHigh;
    DeltaEtaPP = DEta_PP_VarHigh;
 
  } else if (suffix == "31"){
 
    TrackCutsProton->SetPtRange(Proton_pT_VarLow,4.05);
    TrackCutsAntiProton->SetPtRange(Proton_pT_VarLow,4.05);
    TrackCutsProton->SetNClsTPC(Proton_Clusters_VarLow);
    TrackCutsAntiProton->SetNClsTPC(Proton_Clusters_VarLow);
 
    TrackCutsPion->SetPID(AliPID::kPion, 0.5,Pion_Sigma_VarLow);
    TrackCutsAntiPion->SetPID(AliPID::kPion, 0.5,Pion_Sigma_VarLow);
    TrackCutsPion->SetPtRange(Pion_pT_VarHigh,4.0);
    TrackCutsAntiPion->SetPtRange(Pion_pT_VarHigh,4.0);
 
 
    DeltaEtaPP = DPhi_PP_Default;
 
  } else if (suffix == "32"){
 
    TrackCutsProton->SetEtaRange(-Proton_Eta_VarHigh,Proton_Eta_VarHigh);
    TrackCutsAntiProton->SetEtaRange(-Proton_Eta_VarHigh,Proton_Eta_VarHigh);
    TrackCutsProton->SetPtRange(Proton_pT_VarLow,4.05);
    TrackCutsAntiProton->SetPtRange(Proton_pT_VarLow,4.05);
 
    TrackCutsPion->SetPID(AliPID::kPion, 0.5,Pion_Sigma_VarHigh);
    TrackCutsAntiPion->SetPID(AliPID::kPion, 0.5,Pion_Sigma_VarHigh);
    TrackCutsPion->SetEtaRange(-Pion_Eta_VarLow,Pion_Eta_VarLow);
    TrackCutsAntiPion->SetEtaRange(-Pion_Eta_VarLow,Pion_Eta_VarLow);
    TrackCutsPion->SetPtRange(Pion_pT_VarHigh,4.0);
    TrackCutsAntiPion->SetPtRange(Pion_pT_VarHigh,4.0);
    TrackCutsPion->SetNClsTPC(Pion_Clusters_VarHigh);
    TrackCutsAntiPion->SetNClsTPC(Pion_Clusters_VarHigh);
 
    DeltaEtaPPrim = DEta_VarLow;
 
    DeltaPhiPP = DEta_VarLow;
    DeltaEtaPP = DEta_PP_VarHigh;
 
  } else if (suffix == "33"){
 
    TrackCutsProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarLow);
    TrackCutsAntiProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarLow);
    TrackCutsProton->SetEtaRange(-Proton_Eta_VarHigh,Proton_Eta_VarHigh);
    TrackCutsAntiProton->SetEtaRange(-Proton_Eta_VarHigh,Proton_Eta_VarHigh);
    TrackCutsProton->SetPtRange(Proton_pT_VarHigh,4.05);
    TrackCutsAntiProton->SetPtRange(Proton_pT_VarHigh,4.05);
 
    TrackCutsPion->SetPID(AliPID::kPion, 0.5,Pion_Sigma_VarHigh);
    TrackCutsAntiPion->SetPID(AliPID::kPion, 0.5,Pion_Sigma_VarHigh);
    TrackCutsPion->SetNClsTPC(Pion_Clusters_VarHigh);
    TrackCutsAntiPion->SetNClsTPC(Pion_Clusters_VarHigh);
 
    DeltaPhiPPrim = DPhi_VarHigh;
    DeltaEtaPPrim = DEta_VarLow;
 
    DeltaPhiPP = DPhi_PP_VarHigh;
    DeltaEtaPP = DPhi_PP_Default;
 
  } else if (suffix == "34"){
 
    TrackCutsProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarLow);
    TrackCutsAntiProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarLow);
    TrackCutsProton->SetEtaRange(-Proton_Eta_VarLow,Proton_Eta_VarLow);
    TrackCutsAntiProton->SetEtaRange(-Proton_Eta_VarLow,Proton_Eta_VarLow);
    TrackCutsProton->SetPtRange(Proton_pT_VarHigh,4.05);
    TrackCutsAntiProton->SetPtRange(Proton_pT_VarHigh,4.05);
    TrackCutsProton->SetNClsTPC(Proton_Clusters_VarHigh);
    TrackCutsAntiProton->SetNClsTPC(Proton_Clusters_VarHigh);
 
    TrackCutsPion->SetPID(AliPID::kPion, 0.5,Pion_Sigma_VarLow);
    TrackCutsAntiPion->SetPID(AliPID::kPion, 0.5,Pion_Sigma_VarLow);
    TrackCutsPion->SetEtaRange(-Pion_Eta_VarLow,Pion_Eta_VarLow);
    TrackCutsAntiPion->SetEtaRange(-Pion_Eta_VarLow,Pion_Eta_VarLow);
    TrackCutsPion->SetNClsTPC(Pion_Clusters_VarHigh);
    TrackCutsAntiPion->SetNClsTPC(Pion_Clusters_VarHigh);
 
    DeltaPhiPPrim = DPhi_VarHigh;
    DeltaEtaPPrim = DEta_VarHigh;
 
    DeltaPhiPP = DPhi_PP_VarHigh;
    DeltaEtaPP = DEta_PP_VarHigh;
 
  } else if (suffix == "35"){
 
    TrackCutsProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarLow);
    TrackCutsAntiProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarLow);
    TrackCutsProton->SetPtRange(Proton_pT_VarHigh,4.05);
    TrackCutsAntiProton->SetPtRange(Proton_pT_VarHigh,4.05);
 
    TrackCutsPion->SetPID(AliPID::kPion, 0.5,Pion_Sigma_VarLow);
    TrackCutsAntiPion->SetPID(AliPID::kPion, 0.5,Pion_Sigma_VarLow);
    TrackCutsPion->SetEtaRange(-Pion_Eta_VarHigh,Pion_Eta_VarHigh);
    TrackCutsAntiPion->SetEtaRange(-Pion_Eta_VarHigh,Pion_Eta_VarHigh);
    TrackCutsPion->SetPtRange(Pion_pT_VarHigh,4.0);
    TrackCutsAntiPion->SetPtRange(Pion_pT_VarHigh,4.0);
    TrackCutsPion->SetNClsTPC(Pion_Clusters_VarHigh);
    TrackCutsAntiPion->SetNClsTPC(Pion_Clusters_VarHigh);
 
 
    DeltaPhiPP = DEta_VarLow;
 
  } else if (suffix == "36"){
 
 
    TrackCutsPion->SetPID(AliPID::kPion, 0.5,Pion_Sigma_VarLow);
    TrackCutsAntiPion->SetPID(AliPID::kPion, 0.5,Pion_Sigma_VarLow);
    TrackCutsPion->SetEtaRange(-Pion_Eta_VarLow,Pion_Eta_VarLow);
    TrackCutsAntiPion->SetEtaRange(-Pion_Eta_VarLow,Pion_Eta_VarLow);
    TrackCutsPion->SetPtRange(Pion_pT_VarHigh,4.0);
    TrackCutsAntiPion->SetPtRange(Pion_pT_VarHigh,4.0);
    TrackCutsPion->SetNClsTPC(Pion_Clusters_VarHigh);
    TrackCutsAntiPion->SetNClsTPC(Pion_Clusters_VarHigh);
 
    DeltaPhiPPrim = DPhi_VarHigh;
    DeltaEtaPPrim = DEta_VarLow;
 
    DeltaEtaPP = DPhi_PP_Default;
 
  } else if (suffix == "37"){
 
    TrackCutsProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarLow);
    TrackCutsAntiProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarLow);
    TrackCutsProton->SetEtaRange(-Proton_Eta_VarHigh,Proton_Eta_VarHigh);
    TrackCutsAntiProton->SetEtaRange(-Proton_Eta_VarHigh,Proton_Eta_VarHigh);
    TrackCutsProton->SetPtRange(Proton_pT_VarHigh,4.05);
    TrackCutsAntiProton->SetPtRange(Proton_pT_VarHigh,4.05);
 
    TrackCutsPion->SetEtaRange(-Pion_Eta_VarLow,Pion_Eta_VarLow);
    TrackCutsAntiPion->SetEtaRange(-Pion_Eta_VarLow,Pion_Eta_VarLow);
    TrackCutsPion->SetPtRange(Pion_pT_VarHigh,4.0);
    TrackCutsAntiPion->SetPtRange(Pion_pT_VarHigh,4.0);
    TrackCutsPion->SetNClsTPC(Pion_Clusters_VarLow);
    TrackCutsAntiPion->SetNClsTPC(Pion_Clusters_VarLow);
 
    DeltaPhiPPrim = DPhi_VarHigh;
 
    DeltaPhiPP = DEta_VarLow;
    DeltaEtaPP = DPhi_PP_Default;
 
  } else if (suffix == "38"){
 
    TrackCutsProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarLow);
    TrackCutsAntiProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarLow);
    TrackCutsProton->SetEtaRange(-Proton_Eta_VarLow,Proton_Eta_VarLow);
    TrackCutsAntiProton->SetEtaRange(-Proton_Eta_VarLow,Proton_Eta_VarLow);
 
    TrackCutsPion->SetNClsTPC(Pion_Clusters_VarLow);
    TrackCutsAntiPion->SetNClsTPC(Pion_Clusters_VarLow);
 
 
    DeltaPhiPP = DPhi_PP_VarHigh;
 
  } else if (suffix == "39"){
 
    TrackCutsProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarLow);
    TrackCutsAntiProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarLow);
    TrackCutsProton->SetEtaRange(-Proton_Eta_VarLow,Proton_Eta_VarLow);
    TrackCutsAntiProton->SetEtaRange(-Proton_Eta_VarLow,Proton_Eta_VarLow);
    TrackCutsProton->SetNClsTPC(Proton_Clusters_VarLow);
    TrackCutsAntiProton->SetNClsTPC(Proton_Clusters_VarLow);
 
    TrackCutsPion->SetEtaRange(-Pion_Eta_VarHigh,Pion_Eta_VarHigh);
    TrackCutsAntiPion->SetEtaRange(-Pion_Eta_VarHigh,Pion_Eta_VarHigh);
    TrackCutsPion->SetPtRange(Pion_pT_VarHigh,4.0);
    TrackCutsAntiPion->SetPtRange(Pion_pT_VarHigh,4.0);
    TrackCutsPion->SetNClsTPC(Pion_Clusters_VarLow);
    TrackCutsAntiPion->SetNClsTPC(Pion_Clusters_VarLow);
 
    DeltaPhiPPrim = DPhi_VarLow;
 
    DeltaPhiPP = DPhi_PP_VarHigh;
    DeltaEtaPP = DEta_PP_VarHigh;
 
  } else if (suffix == "40"){
 
    TrackCutsProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarLow);
    TrackCutsAntiProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarLow);
    TrackCutsProton->SetEtaRange(-Proton_Eta_VarLow,Proton_Eta_VarLow);
    TrackCutsAntiProton->SetEtaRange(-Proton_Eta_VarLow,Proton_Eta_VarLow);
    TrackCutsProton->SetPtRange(Proton_pT_VarLow,4.05);
    TrackCutsAntiProton->SetPtRange(Proton_pT_VarLow,4.05);
    TrackCutsProton->SetNClsTPC(Proton_Clusters_VarLow);
    TrackCutsAntiProton->SetNClsTPC(Proton_Clusters_VarLow);
 
    TrackCutsPion->SetPID(AliPID::kPion, 0.5,Pion_Sigma_VarHigh);
    TrackCutsAntiPion->SetPID(AliPID::kPion, 0.5,Pion_Sigma_VarHigh);
    TrackCutsPion->SetEtaRange(-Pion_Eta_VarHigh,Pion_Eta_VarHigh);
    TrackCutsAntiPion->SetEtaRange(-Pion_Eta_VarHigh,Pion_Eta_VarHigh);
    TrackCutsPion->SetNClsTPC(Pion_Clusters_VarHigh);
    TrackCutsAntiPion->SetNClsTPC(Pion_Clusters_VarHigh);
 
    DeltaPhiPPrim = DPhi_VarHigh;
 
    DeltaPhiPP = DPhi_PP_VarHigh;
 
  } else if (suffix == "41"){
 
    TrackCutsProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarHigh);
    TrackCutsAntiProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarHigh);
    TrackCutsProton->SetEtaRange(-Proton_Eta_VarLow,Proton_Eta_VarLow);
    TrackCutsAntiProton->SetEtaRange(-Proton_Eta_VarLow,Proton_Eta_VarLow);
    TrackCutsProton->SetPtRange(Proton_pT_VarHigh,4.05);
    TrackCutsAntiProton->SetPtRange(Proton_pT_VarHigh,4.05);
    TrackCutsProton->SetNClsTPC(Proton_Clusters_VarHigh);
    TrackCutsAntiProton->SetNClsTPC(Proton_Clusters_VarHigh);
 
    TrackCutsPion->SetPID(AliPID::kPion, 0.5,Pion_Sigma_VarLow);
    TrackCutsAntiPion->SetPID(AliPID::kPion, 0.5,Pion_Sigma_VarLow);
 
    DeltaPhiPPrim = DPhi_VarHigh;
 
    DeltaPhiPP = DEta_VarLow;
 
  } else if (suffix == "42"){
 
    TrackCutsProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarLow);
    TrackCutsAntiProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarLow);
    TrackCutsProton->SetEtaRange(-Proton_Eta_VarHigh,Proton_Eta_VarHigh);
    TrackCutsAntiProton->SetEtaRange(-Proton_Eta_VarHigh,Proton_Eta_VarHigh);
    TrackCutsProton->SetPtRange(Proton_pT_VarHigh,4.05);
    TrackCutsAntiProton->SetPtRange(Proton_pT_VarHigh,4.05);
    TrackCutsProton->SetNClsTPC(Proton_Clusters_VarHigh);
    TrackCutsAntiProton->SetNClsTPC(Proton_Clusters_VarHigh);
 
    TrackCutsPion->SetPtRange(Pion_pT_VarHigh,4.0);
    TrackCutsAntiPion->SetPtRange(Pion_pT_VarHigh,4.0);
    TrackCutsPion->SetNClsTPC(Pion_Clusters_VarHigh);
    TrackCutsAntiPion->SetNClsTPC(Pion_Clusters_VarHigh);
 
    DeltaEtaPPrim = DEta_VarLow;
 
    DeltaPhiPP = DPhi_PP_VarHigh;
    DeltaEtaPP = DPhi_PP_Default;
 
  } else if (suffix == "43"){
 
    TrackCutsProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarLow);
    TrackCutsAntiProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarLow);
    TrackCutsProton->SetEtaRange(-Proton_Eta_VarHigh,Proton_Eta_VarHigh);
    TrackCutsAntiProton->SetEtaRange(-Proton_Eta_VarHigh,Proton_Eta_VarHigh);
    TrackCutsProton->SetPtRange(Proton_pT_VarLow,4.05);
    TrackCutsAntiProton->SetPtRange(Proton_pT_VarLow,4.05);
 
    TrackCutsPion->SetPID(AliPID::kPion, 0.5,Pion_Sigma_VarLow);
    TrackCutsAntiPion->SetPID(AliPID::kPion, 0.5,Pion_Sigma_VarLow);
    TrackCutsPion->SetPtRange(Pion_pT_VarLow,4.0);
    TrackCutsAntiPion->SetPtRange(Pion_pT_VarLow,4.0);
    TrackCutsPion->SetNClsTPC(Pion_Clusters_VarLow);
    TrackCutsAntiPion->SetNClsTPC(Pion_Clusters_VarLow);
 
    DeltaPhiPPrim = DPhi_VarLow;
 
    DeltaPhiPP = DPhi_PP_VarHigh;
    DeltaEtaPP = DEta_PP_VarHigh;
 
  } else if (suffix == "44"){
 
    TrackCutsProton->SetEtaRange(-Proton_Eta_VarLow,Proton_Eta_VarLow);
    TrackCutsAntiProton->SetEtaRange(-Proton_Eta_VarLow,Proton_Eta_VarLow);
    TrackCutsProton->SetPtRange(Proton_pT_VarLow,4.05);
    TrackCutsAntiProton->SetPtRange(Proton_pT_VarLow,4.05);
    TrackCutsProton->SetNClsTPC(Proton_Clusters_VarLow);
    TrackCutsAntiProton->SetNClsTPC(Proton_Clusters_VarLow);
 
    TrackCutsPion->SetPID(AliPID::kPion, 0.5,Pion_Sigma_VarLow);
    TrackCutsAntiPion->SetPID(AliPID::kPion, 0.5,Pion_Sigma_VarLow);
    TrackCutsPion->SetNClsTPC(Pion_Clusters_VarLow);
    TrackCutsAntiPion->SetNClsTPC(Pion_Clusters_VarLow);
 
 
    DeltaPhiPP = DEta_VarLow;
    DeltaEtaPP = DEta_PP_VarHigh;
 
  }


  TString addon = "";
  if (trigger == "kINT7") {
    addon += "kINT7";
  } else if (trigger == "kHM") {
    addon += "HM";
  }

  TString file = AliAnalysisManager::GetCommonFileName();

  TString EvtCutsName = Form("%s_EvtCuts_%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *coutputEvtCuts = mgr->CreateContainer(
      EvtCutsName.Data(), TList::Class(), AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), EvtCutsName.Data()));

  TString TrackCutsName = Form("%s_TrackCuts_%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *couputTrkCuts = mgr->CreateContainer(
      TrackCutsName.Data(), TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), TrackCutsName.Data()));

  TString AntiTrackCutsName = Form("%s_AntiTrackCuts_%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *coutputAntiTrkCuts = mgr->CreateContainer(
      AntiTrackCutsName.Data(), TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), AntiTrackCutsName.Data()));

  AliAnalysisDataContainer *coutputPionCuts;
  TString PionCutsName = Form("%s_PionCuts_%s", addon.Data(), suffix.Data());
  coutputPionCuts = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      PionCutsName.Data(),
      TList::Class(), AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), PionCutsName.Data()));

  AliAnalysisDataContainer *coutputAntiPionCuts;
  TString AntiPionCutsName = Form("%s_AntiPionCuts_%s", addon.Data(), suffix.Data());
  coutputAntiPionCuts = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      AntiPionCutsName.Data(),
      TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), AntiPionCutsName.Data()));

  AliAnalysisDataContainer *coutputResults;
  TString ResultsName = Form("%s_Results_%s", addon.Data(), suffix.Data());
  coutputResults = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      ResultsName.Data(),
      TList::Class(), AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), ResultsName.Data()));

  AliAnalysisDataContainer *coutputResultsQA;
  TString ResultsQAName = Form("%s_ResultsQA_%s", addon.Data(), suffix.Data());
  coutputResultsQA = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      ResultsQAName.Data(),
      TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), ResultsQAName.Data()));

  AliAnalysisDataContainer *coutputResultsSample;
  TString ResultsSampleName = Form("%s_ResultsSample_%s", addon.Data(), suffix.Data());
  coutputResultsSample = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      ResultsSampleName.Data(),
      TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), ResultsSampleName.Data()));

  AliAnalysisDataContainer *coutputResultsSampleQA;
  TString ResultsSampleQAName = Form("%s_ResultsSampleQA_%s", addon.Data(), suffix.Data());
  coutputResultsSampleQA = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      ResultsSampleQAName.Data(),
      TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), ResultsSampleQAName.Data()));

  AliAnalysisDataContainer *coutputTrkCutsMC;
  AliAnalysisDataContainer *coutputAntiTrkCutsMC;
  AliAnalysisDataContainer *coutputPionCutsMC;
  AliAnalysisDataContainer *coutputAntiPionCutsMC;
  if (isMC) {
    TString TrkCutsMCName = Form("%s_TrkCutsMC_%s", addon.Data(), suffix.Data());
    coutputTrkCutsMC = mgr->CreateContainer(
        //@suppress("Invalid arguments") it works ffs
        TrkCutsMCName.Data(),
        TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), TrkCutsMCName.Data()));

    TString AntiTrkCutsMCName = Form("%s_AntiTrkCutsMC_%s", addon.Data(), suffix.Data());
    coutputAntiTrkCutsMC = mgr->CreateContainer(
        //@suppress("Invalid arguments") it works ffs
        AntiTrkCutsMCName.Data(),
        TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), AntiTrkCutsMCName.Data()));

    TString PionCutsMCName = Form("%s_PionCutsMC_%s", addon.Data(), suffix.Data());
    coutputPionCutsMC = mgr->CreateContainer(
        //@suppress("Invalid arguments") it works ffs
        PionCutsMCName.Data(),
        TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), PionCutsMCName.Data()));

    TString AntiPionCutsMCName = Form("%s_AntiPionCutsMC_%s", addon.Data(), suffix.Data());
    coutputAntiPionCutsMC = mgr->CreateContainer(
        //@suppress("Invalid arguments") it works ffs
        AntiPionCutsMCName.Data(),
        TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), AntiPionCutsMCName.Data()));

  }

  TString ThreeBodyName = Form("%s_ThreeBodyProtonPion_%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *coutputThreeBody = mgr->CreateContainer(
    //@suppress("Invalid arguments") it works ffs
    ThreeBodyName.Data(),
    TList::Class(),
    AliAnalysisManager::kOutputContainer,
    Form("%s:%s", file.Data(), ThreeBodyName.Data()));

    AliAnalysisTaskAODThreeBodyProtonPrimary* taskAOD;
    TString taskName = "ThreeBodyProtonPion";
    taskAOD= new AliAnalysisTaskAODThreeBodyProtonPrimary(Form("femtoNanoThreeBodyProtonPion_%s", taskName.Data()), isMC);

    if (trigger == "kINT7") {
      taskAOD->SelectCollisionCandidates(AliVEvent::kINT7);
      std::cout << "Added kINT7 Trigger \n";
    } else if (trigger == "kHM") {
      taskAOD->SelectCollisionCandidates(AliVEvent::kHighMultV0);
      std::cout << "Added kHighMult Trigger \n";
    } else {
      std::cout << "=====================================================================" << std::endl;
      std::cout << "=====================================================================" << std::endl;
      std::cout << "Centrality Estimator not set, fix it else your Results will be empty!" << std::endl;
      std::cout << "=====================================================================" << std::endl;
      std::cout << "=====================================================================" << std::endl;
    }

    if (!fullBlastQA)
    {
      taskAOD->SetRunTaskLightWeight(true);
    }

  

    taskAOD->SetEventCuts(evtCuts);
    taskAOD->SetProtonCuts(TrackCutsProton);
    taskAOD->SetAntiProtonCuts(TrackCutsAntiProton);
    taskAOD->SetPrimaryCuts(TrackCutsPion);
    taskAOD->SetAntiPrimaryCuts(TrackCutsAntiPion);
    taskAOD->SetCorrelationConfig(config);
  
    taskAOD->SetDoOnlyThreeBody(DoThreeBody);
    taskAOD->SetStandardMixing(StandardMixing);

    taskAOD->SetturnoffClosePairRejectionCompletely(turnoffClosePairRejectionCompletely);
    taskAOD->SetClosePairRejectionForAll(ClosePairRejectionForAll);
    taskAOD->SetDeltaPhiMaxPP(DeltaPhiPP);    
    taskAOD->SetDeltaEtaMaxPP(DeltaEtaPP);
    taskAOD->SetDeltaPhiMaxPPrim(DeltaPhiPPrim);
    taskAOD->SetDeltaEtaMaxPPrim(DeltaEtaPPrim);
    taskAOD->SetDeltaPhiMaxPAPrim(0.0); 
    taskAOD->SetDeltaEtaMaxPAPrim(0.0); 

    taskAOD->SetQ3LimitForDeltaPhiDeltaEta(Q3LimitForDeltaPhiDeltaEta);
    taskAOD->SetQ3cutValue(Q3cutValue); 
    taskAOD->SetQ3MinValue(Q3MinValue);

    taskAOD->SetRunThreeBodyHistograms(true);
    taskAOD->SetRunPlotInvMass(RunPlotInvMass);
    taskAOD->SetRunPlotPhiTheta(RunPlotPhiTheta);
    taskAOD->SetRunPlotMult(RunPlotMult);
    taskAOD->SetRunPairMultThreeBody(RunPairMultThreeBody);
    taskAOD->SetRunmTPlots(RunmTPlots);

    taskAOD->SetMCAndReso(isMC,RemoveResonances); 
    taskAOD->SetRunProjector(GetProjector);
    taskAOD->SetGetMomentumResolution(GetMomentumResolution);
    
    mgr->AddTask(taskAOD);

    mgr->ConnectInput(taskAOD, 0, cinput);
    mgr->ConnectOutput(taskAOD, 1, coutputEvtCuts);
    mgr->ConnectOutput(taskAOD, 2, couputTrkCuts);
    mgr->ConnectOutput(taskAOD, 3, coutputAntiTrkCuts);
    mgr->ConnectOutput(taskAOD, 4, coutputPionCuts);
    mgr->ConnectOutput(taskAOD, 5, coutputAntiPionCuts);
    mgr->ConnectOutput(taskAOD, 6, coutputResults);
    mgr->ConnectOutput(taskAOD, 7, coutputResultsQA);
    mgr->ConnectOutput(taskAOD, 8, coutputResultsSample);
    mgr->ConnectOutput(taskAOD, 9, coutputResultsSampleQA);
    mgr->ConnectOutput(taskAOD, 10, coutputThreeBody);
    if (isMC) {
      mgr->ConnectOutput(taskAOD, 11, coutputTrkCutsMC);
      mgr->ConnectOutput(taskAOD, 12, coutputAntiTrkCutsMC);
      mgr->ConnectOutput(taskAOD, 13, coutputPionCutsMC);
      mgr->ConnectOutput(taskAOD, 14, coutputAntiPionCutsMC);
    }

    return taskAOD;
}
