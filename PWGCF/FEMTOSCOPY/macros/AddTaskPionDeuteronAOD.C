#include "TROOT.h"
#include "TSystem.h"

AliAnalysisTaskSE *AddTaskPionDeuteronAOD(bool isMC = false,         // 1
                                          TString trigger = "kINT7", // 2
                                          bool fullBlastQA = true,   // 3
                                          bool MoreChecks = false,   // 4
                                          bool systmatics = false,   // 5
                                          bool DoFunWithPhaseSpace = false,//6
                                          bool DopTPionNucleonHisto = false,
                                          float pTOnepTTwokStarCutOff = 3., //7
                                          int mTBinningChoice = 0, //8
                                          const char *cutVariation = "0") {
  TString suffix = TString::Format("%s", cutVariation);
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    printf("No analysis manager to connect to!\n");
    return nullptr;
  }
  if (!mgr->GetInputEventHandler()) {
    printf("This task requires an input event handler!\n");
    return nullptr;
  }
  // Oton's selections for deuterons since its large range more statistics at
  // the end
  AliFemtoDreamEventCuts *evtCuts = AliFemtoDreamEventCuts::StandardCutsRun2();
  evtCuts->CleanUpMult(false, false, false, true);
  AliFemtoDreamTrackCuts *TrackCutsDeuteronDCA =
      AliFemtoDreamTrackCuts::PrimDeuteronCuts(isMC, true, false, false);
  TrackCutsDeuteronDCA->SetCutCharge(1);
  TrackCutsDeuteronDCA->SetPtRange(0.8, 2.4);
  TrackCutsDeuteronDCA->SetFilterBit(128);

  AliFemtoDreamTrackCuts *TrackCutsDeuteronMass =
      AliFemtoDreamTrackCuts::PrimDeuteronCuts(isMC, true, false, false);
  TrackCutsDeuteronMass->SetCutCharge(1);
  TrackCutsDeuteronMass->SetPID(AliPID::kDeuteron, 999.0);

  AliFemtoDreamTrackCuts *TrackCutsAntiDeuteronDCA =
      AliFemtoDreamTrackCuts::PrimDeuteronCuts(isMC, true, false, false);
  TrackCutsAntiDeuteronDCA->SetCutCharge(-1);
  TrackCutsAntiDeuteronDCA->SetPtRange(0.8, 2.4);
  TrackCutsAntiDeuteronDCA->SetFilterBit(128);

  AliFemtoDreamTrackCuts *TrackCutsAntiDeuteronMass =
      AliFemtoDreamTrackCuts::PrimDeuteronCuts(isMC, true, false, false);
  TrackCutsAntiDeuteronMass->SetCutCharge(-1);
  TrackCutsAntiDeuteronMass->SetPID(AliPID::kDeuteron, 999.0);


  //---------- posPions--------------------
  AliFemtoDreamTrackCuts *posPions =
      AliFemtoDreamTrackCuts::PrimPionCuts(isMC, true, false, false);
  posPions->SetCutCharge(1);
  posPions->SetFilterBit(128);
  posPions->SetPtRange(0.14, 4.0);

  //---------- NegPions--------------------
  AliFemtoDreamTrackCuts *NegPions =
      AliFemtoDreamTrackCuts::PrimPionCuts(isMC, true, false, false);
  NegPions->SetCutCharge(-1);
  NegPions->SetFilterBit(128);
  NegPions->SetPtRange(0.14, 4.0);

  // a few more checks for deuterons
  if (MoreChecks) {
    // Oton's selections for deuterons
    if (suffix == "1") {
      TrackCutsDeuteronDCA->SetPtRange(0.5, 1.4);
      TrackCutsAntiDeuteronDCA->SetPtRange(0.5, 1.4);
    } else if (suffix == "2") {
      TrackCutsDeuteronDCA->SetPtRange(0.8, 1.4);
      TrackCutsAntiDeuteronDCA->SetPtRange(0.8, 1.4);
      TrackCutsDeuteronDCA->SetPIDkd(false, false, 3, 3);
      TrackCutsAntiDeuteronDCA->SetPIDkd(false, false, 3, 3);
    } else if (suffix == "3") {
      TrackCutsDeuteronDCA->SetPtRange(0.5, 2.4);
      TrackCutsAntiDeuteronDCA->SetPtRange(0.5, 2.4);
      TrackCutsDeuteronDCA->SetPIDkd(false, false, 3, 3);
      TrackCutsAntiDeuteronDCA->SetPIDkd(false, false, 3, 3);
    } else if (suffix == "4"){
      posPions->SetPtRange(0.16, 4.0);
      NegPions->SetPtRange(0.16, 4.0);
    }
  } else {
    TrackCutsDeuteronDCA->SetPtRange(0.8, 2.4);
    TrackCutsAntiDeuteronDCA->SetPtRange(0.8, 2.4);
    TrackCutsDeuteronDCA->SetPIDkd(false, false, 3, 3);
    TrackCutsAntiDeuteronDCA->SetPIDkd(false, false, 3, 3);
  }




  std::vector<int> PDGParticles;
  PDGParticles.push_back(211);        // pi+
  PDGParticles.push_back(211);        // pi-
  PDGParticles.push_back(1000010020); // d
  PDGParticles.push_back(1000010020); // dbar

  std::vector<bool> closeRejection;
  std::vector<float> mTBins;
  if(mTBinningChoice == 0){
  mTBins.push_back(1.06);
  mTBins.push_back(1.19);
  mTBins.push_back(1.3);
  mTBins.push_back(999.);
  }else if(mTBinningChoice == 1){
    mTBins.push_back(1.03);
    mTBins.push_back(1.18);
    mTBins.push_back(1.28);
    mTBins.push_back(1.4);
    mTBins.push_back(1.55);
    mTBins.push_back(999.);
  }else if(mTBinningChoice == 2){
    mTBins.push_back(1.03);
    mTBins.push_back(1.16);
    mTBins.push_back(1.24);
    mTBins.push_back(1.38);
    mTBins.push_back(1.52);
    mTBins.push_back(999.);
  }
  else{
    mTBins.push_back(1.03);
    mTBins.push_back(1.20);
    mTBins.push_back(1.38);
    mTBins.push_back(1.50);
    mTBins.push_back(999.);
  }

  std::vector<int> pairQA;
  // pairs:
  // pipi             0
  // pi bar p        1
  // pi d            2
  // pi bar d        3
  // bar pi bar pi    4
  // bar pi d        5
  // bar pi bar d    6
  // d d            7
  // d bar d        8
  // bar d bar d    9
  const int nPairs = 10;
  for (int i = 0; i < nPairs; ++i) {
    closeRejection.push_back(false);
    pairQA.push_back(0);
  }
  closeRejection[2] = true; // pid
  closeRejection[3] = true; // pi bard
  closeRejection[5] = true; // barpi barpi
  closeRejection[6] = true; // barpi d
 // closeRejection[7] = true; // dd
 // closeRejection[9] = true; // bard bar

  pairQA[2] = 11;           // pid
  pairQA[3] = 11;           // pi bard
  pairQA[5] = 11;           // barpi barpi
  pairQA[6] = 11;           // barpi bard
 // pairQA[7] = 11;           // dd
 // pairQA[9] = 11;           // bard bard
  // We need to set the ZVtx bins
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
  // The Multiplicity bins are set here
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

  std::vector<int> NBins;
  NBins.push_back(750);
  NBins.push_back(750);
  NBins.push_back(750);
  NBins.push_back(750);
  NBins.push_back(750);
  NBins.push_back(750);
  NBins.push_back(750);
  NBins.push_back(750);
  NBins.push_back(750);
  NBins.push_back(750);
  std::vector<float> kMin;
  // minimum k* value
  kMin.push_back(0.);
  kMin.push_back(0.);
  kMin.push_back(0.);
  kMin.push_back(0.);
  kMin.push_back(0.);
  kMin.push_back(0.);
  kMin.push_back(0.);
  kMin.push_back(0.);
  kMin.push_back(0.);
  kMin.push_back(0.);
  // maximum k* value
  std::vector<float> kMax;
  kMax.push_back(3.);
  kMax.push_back(3.);
  kMax.push_back(3.);
  kMax.push_back(3.);
  kMax.push_back(3.);
  kMax.push_back(3.);
  kMax.push_back(3.);
  kMax.push_back(3.);
  kMax.push_back(3.);
  kMax.push_back(3.);

  AliFemtoDreamCollConfig *config =
      new AliFemtoDreamCollConfig("Femto", "Femto", false);
  config->SetZBins(ZVtxBins);
  config->SetMultBins(MultBins);
  config->SetMultBinning(true);
  config->SetPDGCodes(PDGParticles);
  config->SetNBinsHist(NBins);
  config->SetMinKRel(kMin);
  config->SetMaxKRel(kMax);
  config->SetClosePairRejection(closeRejection);
  config->SetExtendedQAPairs(pairQA);
  config->SetDeltaEtaMax(0.017); // and here you set the actual values
  config->SetDeltaPhiMax(0.04);  // and here you set the actual values
  config->SetMixingDepth(10);
  config->SetmTBins(mTBins);
  config->SetDomTMultBinning(true);
  config->SetmTBinning(true);
  config->SetMultiplicityEstimator(AliFemtoDreamEvent::kRef08);
  
  if(DoFunWithPhaseSpace){
    config->SetpTOnepTTwokStarPlotsmT(true, pTOnepTTwokStarCutOff);
  } else {
    config->SetpTOnepTTwokStarPlotsmT(false, pTOnepTTwokStarCutOff);
  }
  if(DopTPionNucleonHisto){
    config->SetpTPionNucleonkStarPlotsmT(true);
  }
  if (isMC) {
    config->SetMomentumResolution(true);
  }
  if (fullBlastQA) {
    config->SetPtQA(true);
    config->SetdPhidEtaPlots(
        true); /// to check delta eta deltaphi plot carefully
    config->SetdPhidEtaPlotsSmallK(true);
    config->SetPhiEtaBinnign(true);
  } else {
    evtCuts->SetMinimalBooking(true);
    NegPions->SetMinimalBooking(true);
    posPions->SetMinimalBooking(true);
    TrackCutsDeuteronDCA->SetMinimalBooking(true);
    TrackCutsDeuteronMass->SetMinimalBooking(true);
    TrackCutsAntiDeuteronDCA->SetMinimalBooking(true);
    TrackCutsAntiDeuteronMass->SetMinimalBooking(true);
    config->SetMinimalBookingME(true);
    config->SetMinimalBookingSample(true);
  }
  // here comes systematics for beatiful pi-d anaylsis: doors to break the
  // riddle of thermal vs Coalescence:D

  if (systmatics) {
    //=============================================================================================================
    // Systematics
    Float_t Deuteron_pT_VarLow = 0.75;
    Float_t Deuteron_pT_VarHigh = 0.85;
    Float_t Deuteron_Eta_VarLow = 0.76;
    Float_t Deuteron_Eta_VarHigh = 0.84;
    Float_t Deuteron_Clusters_VarLow = 70.;
    Float_t Deuteron_Clusters_VarHigh = 90.;
    Float_t Deuteron_Sigma_VarLow = 2.5;
    Float_t Deuteron_Sigma_VarHigh = 3.5;

    Float_t Pion_pT_VarLow = 0.12;
    Float_t Pion_pT_VarHigh = 0.15;
    Float_t Pion_Eta_VarLow = 0.77;  // old 0.7
    Float_t Pion_Eta_VarHigh = 0.85; // old 0.9
    Float_t Pion_Clusters_VarLow = 70;
    Float_t Pion_Clusters_VarHigh = 90;
    Float_t Pion_Sigma_VarLow = 2.5;  // old 2.7;
    Float_t Pion_Sigma_VarHigh = 3.5; // old 3.5

    Float_t DPhi_VarLow = 0.035;
    Float_t DPhi_VarHigh = 0.045;
    Float_t DEta_VarLow = 0.015;
    Float_t DEta_VarHigh = 0.019;

    /// look carefully at the variations
    if (suffix == "1") {

      TrackCutsDeuteronDCA->SetPIDkd(false, false, Deuteron_Sigma_VarLow,
                                     Deuteron_Sigma_VarLow);
      TrackCutsAntiDeuteronDCA->SetPIDkd(false, false, Deuteron_Sigma_VarLow,
                                         Deuteron_Sigma_VarLow);
      TrackCutsDeuteronDCA->SetEtaRange(-Deuteron_Eta_VarLow,
                                        Deuteron_Eta_VarLow);
      TrackCutsAntiDeuteronDCA->SetEtaRange(-Deuteron_Eta_VarLow,
                                            Deuteron_Eta_VarLow);

      posPions->SetPID(AliPID::kPion, 0.5, Pion_Sigma_VarHigh);
      NegPions->SetPID(AliPID::kPion, 0.5, Pion_Sigma_VarHigh);
      posPions->SetEtaRange(-Pion_Eta_VarHigh, Pion_Eta_VarHigh);
      NegPions->SetEtaRange(-Pion_Eta_VarHigh, Pion_Eta_VarHigh);
      posPions->SetPtRange(Pion_pT_VarLow, 4.0);
      NegPions->SetPtRange(Pion_pT_VarLow, 4.0);

      config->SetDeltaPhiMax(DPhi_VarHigh);

    } else if (suffix == "2") {

      TrackCutsDeuteronDCA->SetPtRange(Deuteron_pT_VarLow, 2.4);
      TrackCutsAntiDeuteronDCA->SetPtRange(Deuteron_pT_VarLow, 2.4);
      TrackCutsDeuteronDCA->SetNClsTPC(Deuteron_Clusters_VarLow);
      TrackCutsAntiDeuteronDCA->SetNClsTPC(Deuteron_Clusters_VarLow);

      posPions->SetPID(AliPID::kPion, 0.5, Pion_Sigma_VarLow);
      NegPions->SetPID(AliPID::kPion, 0.5, Pion_Sigma_VarLow);
      posPions->SetEtaRange(-Pion_Eta_VarHigh, Pion_Eta_VarHigh);
      NegPions->SetEtaRange(-Pion_Eta_VarHigh, Pion_Eta_VarHigh);
      posPions->SetPtRange(Pion_pT_VarHigh, 4.0);
      NegPions->SetPtRange(Pion_pT_VarHigh, 4.0);
      posPions->SetNClsTPC(Pion_Clusters_VarLow);
      NegPions->SetNClsTPC(Pion_Clusters_VarLow);

      config->SetDeltaPhiMax(DPhi_VarHigh);

    } else if (suffix == "3") {

      TrackCutsDeuteronDCA->SetPtRange(Deuteron_pT_VarHigh, 2.4);
      TrackCutsAntiDeuteronDCA->SetPtRange(Deuteron_pT_VarHigh, 2.4);
      TrackCutsDeuteronDCA->SetNClsTPC(Deuteron_Clusters_VarLow);
      TrackCutsAntiDeuteronDCA->SetNClsTPC(Deuteron_Clusters_VarLow);

      posPions->SetPID(AliPID::kPion, 0.5, Pion_Sigma_VarHigh);
      NegPions->SetPID(AliPID::kPion, 0.5, Pion_Sigma_VarHigh);
      posPions->SetEtaRange(-Pion_Eta_VarHigh, Pion_Eta_VarHigh);
      NegPions->SetEtaRange(-Pion_Eta_VarHigh, Pion_Eta_VarHigh);
      posPions->SetPtRange(Pion_pT_VarLow, 4.0);
      NegPions->SetPtRange(Pion_pT_VarLow, 4.0);

    } else if (suffix == "4") {

      TrackCutsDeuteronDCA->SetPIDkd(false, false, Deuteron_Sigma_VarLow,
                                     Deuteron_Sigma_VarLow);
      TrackCutsAntiDeuteronDCA->SetPIDkd(false, false, Deuteron_Sigma_VarLow,
                                         Deuteron_Sigma_VarLow);
      TrackCutsDeuteronDCA->SetPtRange(Deuteron_pT_VarHigh, 2.4);
      TrackCutsAntiDeuteronDCA->SetPtRange(Deuteron_pT_VarHigh, 2.4);
      TrackCutsDeuteronDCA->SetNClsTPC(Deuteron_Clusters_VarLow);
      TrackCutsAntiDeuteronDCA->SetNClsTPC(Deuteron_Clusters_VarLow);

    } else if (suffix == "5") {

      TrackCutsDeuteronDCA->SetEtaRange(-Deuteron_Eta_VarHigh,
                                        Deuteron_Eta_VarHigh);
      TrackCutsAntiDeuteronDCA->SetEtaRange(-Deuteron_Eta_VarHigh,
                                            Deuteron_Eta_VarHigh);

      posPions->SetPID(AliPID::kPion, 0.5, Pion_Sigma_VarHigh);
      NegPions->SetPID(AliPID::kPion, 0.5, Pion_Sigma_VarHigh);
      posPions->SetEtaRange(-Pion_Eta_VarHigh, Pion_Eta_VarHigh);
      NegPions->SetEtaRange(-Pion_Eta_VarHigh, Pion_Eta_VarHigh);
      posPions->SetNClsTPC(Pion_Clusters_VarHigh);
      NegPions->SetNClsTPC(Pion_Clusters_VarHigh);

      config->SetDeltaPhiMax(DPhi_VarHigh);
      config->SetDeltaEtaMax(DEta_VarHigh);

    } else if (suffix == "6") {

      TrackCutsDeuteronDCA->SetPIDkd(false, false, Deuteron_Sigma_VarLow,
                                     Deuteron_Sigma_VarLow);
      TrackCutsAntiDeuteronDCA->SetPIDkd(false, false, Deuteron_Sigma_VarLow,
                                         Deuteron_Sigma_VarLow);
      TrackCutsDeuteronDCA->SetEtaRange(-Deuteron_Eta_VarHigh,
                                        Deuteron_Eta_VarHigh);
      TrackCutsAntiDeuteronDCA->SetEtaRange(-Deuteron_Eta_VarHigh,
                                            Deuteron_Eta_VarHigh);
      TrackCutsDeuteronDCA->SetPtRange(Deuteron_pT_VarHigh, 2.4);
      TrackCutsAntiDeuteronDCA->SetPtRange(Deuteron_pT_VarHigh, 2.4);
      TrackCutsDeuteronDCA->SetNClsTPC(Deuteron_Clusters_VarHigh);
      TrackCutsAntiDeuteronDCA->SetNClsTPC(Deuteron_Clusters_VarHigh);

      posPions->SetPID(AliPID::kPion, 0.5, Pion_Sigma_VarHigh);
      NegPions->SetPID(AliPID::kPion, 0.5, Pion_Sigma_VarHigh);
      posPions->SetEtaRange(-Pion_Eta_VarHigh, Pion_Eta_VarHigh);
      NegPions->SetEtaRange(-Pion_Eta_VarHigh, Pion_Eta_VarHigh);
      posPions->SetPtRange(Pion_pT_VarHigh, 4.0);
      NegPions->SetPtRange(Pion_pT_VarHigh, 4.0);
      posPions->SetNClsTPC(Pion_Clusters_VarHigh);
      NegPions->SetNClsTPC(Pion_Clusters_VarHigh);

      config->SetDeltaPhiMax(DPhi_VarHigh);
      config->SetDeltaEtaMax(DEta_VarLow);

    } else if (suffix == "7") {

      TrackCutsDeuteronDCA->SetPIDkd(false, false, Deuteron_Sigma_VarLow,
                                     Deuteron_Sigma_VarLow);
      TrackCutsAntiDeuteronDCA->SetPIDkd(false, false, Deuteron_Sigma_VarLow,
                                         Deuteron_Sigma_VarLow);
      TrackCutsDeuteronDCA->SetPtRange(Deuteron_pT_VarHigh, 2.4);
      TrackCutsAntiDeuteronDCA->SetPtRange(Deuteron_pT_VarHigh, 2.4);
      TrackCutsDeuteronDCA->SetNClsTPC(Deuteron_Clusters_VarHigh);
      TrackCutsAntiDeuteronDCA->SetNClsTPC(Deuteron_Clusters_VarHigh);

      posPions->SetPID(AliPID::kPion, 0.5, Pion_Sigma_VarLow);
      NegPions->SetPID(AliPID::kPion, 0.5, Pion_Sigma_VarLow);
      posPions->SetPtRange(Pion_pT_VarLow, 4.0);
      NegPions->SetPtRange(Pion_pT_VarLow, 4.0);

      config->SetDeltaEtaMax(DEta_VarHigh);

    } else if (suffix == "8") {

      TrackCutsDeuteronDCA->SetPIDkd(false, false, Deuteron_Sigma_VarHigh,
                                     Deuteron_Sigma_VarHigh);
      TrackCutsAntiDeuteronDCA->SetPIDkd(false, false, Deuteron_Sigma_VarHigh,
                                         Deuteron_Sigma_VarHigh);
      TrackCutsDeuteronDCA->SetEtaRange(-Deuteron_Eta_VarLow,
                                        Deuteron_Eta_VarLow);
      TrackCutsAntiDeuteronDCA->SetEtaRange(-Deuteron_Eta_VarLow,
                                            Deuteron_Eta_VarLow);
      TrackCutsDeuteronDCA->SetPtRange(Deuteron_pT_VarHigh, 2.4);
      TrackCutsAntiDeuteronDCA->SetPtRange(Deuteron_pT_VarHigh, 2.4);
      TrackCutsDeuteronDCA->SetNClsTPC(Deuteron_Clusters_VarLow);
      TrackCutsAntiDeuteronDCA->SetNClsTPC(Deuteron_Clusters_VarLow);

      posPions->SetPID(AliPID::kPion, 0.5, Pion_Sigma_VarHigh);
      NegPions->SetPID(AliPID::kPion, 0.5, Pion_Sigma_VarHigh);
      posPions->SetEtaRange(-Pion_Eta_VarLow, Pion_Eta_VarLow);
      NegPions->SetEtaRange(-Pion_Eta_VarLow, Pion_Eta_VarLow);
      posPions->SetPtRange(Pion_pT_VarHigh, 4.0);
      NegPions->SetPtRange(Pion_pT_VarHigh, 4.0);
      posPions->SetNClsTPC(Pion_Clusters_VarHigh);
      NegPions->SetNClsTPC(Pion_Clusters_VarHigh);

      config->SetDeltaEtaMax(DEta_VarHigh);

    } else if (suffix == "9") {

      TrackCutsDeuteronDCA->SetPIDkd(false, false, Deuteron_Sigma_VarHigh,
                                     Deuteron_Sigma_VarHigh);
      TrackCutsAntiDeuteronDCA->SetPIDkd(false, false, Deuteron_Sigma_VarHigh,
                                         Deuteron_Sigma_VarHigh);
      TrackCutsDeuteronDCA->SetEtaRange(-Deuteron_Eta_VarLow,
                                        Deuteron_Eta_VarLow);
      TrackCutsAntiDeuteronDCA->SetEtaRange(-Deuteron_Eta_VarLow,
                                            Deuteron_Eta_VarLow);
      TrackCutsDeuteronDCA->SetPtRange(Deuteron_pT_VarLow, 2.4);
      TrackCutsAntiDeuteronDCA->SetPtRange(Deuteron_pT_VarLow, 2.4);
      TrackCutsDeuteronDCA->SetNClsTPC(Deuteron_Clusters_VarLow);
      TrackCutsAntiDeuteronDCA->SetNClsTPC(Deuteron_Clusters_VarLow);

      posPions->SetEtaRange(-Pion_Eta_VarLow, Pion_Eta_VarLow);
      NegPions->SetEtaRange(-Pion_Eta_VarLow, Pion_Eta_VarLow);
      posPions->SetPtRange(Pion_pT_VarLow, 4.0);
      NegPions->SetPtRange(Pion_pT_VarLow, 4.0);

      config->SetDeltaPhiMax(DPhi_VarLow);

    } else if (suffix == "10") {

      TrackCutsDeuteronDCA->SetPIDkd(false, false, Deuteron_Sigma_VarLow,
                                     Deuteron_Sigma_VarLow);
      TrackCutsAntiDeuteronDCA->SetPIDkd(false, false, Deuteron_Sigma_VarLow,
                                         Deuteron_Sigma_VarLow);
      TrackCutsDeuteronDCA->SetEtaRange(-Deuteron_Eta_VarHigh,
                                        Deuteron_Eta_VarHigh);
      TrackCutsAntiDeuteronDCA->SetEtaRange(-Deuteron_Eta_VarHigh,
                                            Deuteron_Eta_VarHigh);
      TrackCutsDeuteronDCA->SetPtRange(Deuteron_pT_VarLow, 2.4);
      TrackCutsAntiDeuteronDCA->SetPtRange(Deuteron_pT_VarLow, 2.4);

      posPions->SetEtaRange(-Pion_Eta_VarLow, Pion_Eta_VarLow);
      NegPions->SetEtaRange(-Pion_Eta_VarLow, Pion_Eta_VarLow);
      posPions->SetNClsTPC(Pion_Clusters_VarHigh);
      NegPions->SetNClsTPC(Pion_Clusters_VarHigh);

      config->SetDeltaPhiMax(DPhi_VarLow);
      config->SetDeltaEtaMax(DEta_VarHigh);

    } else if (suffix == "11") {

      TrackCutsDeuteronDCA->SetPIDkd(false, false, Deuteron_Sigma_VarLow,
                                     Deuteron_Sigma_VarLow);
      TrackCutsAntiDeuteronDCA->SetPIDkd(false, false, Deuteron_Sigma_VarLow,
                                         Deuteron_Sigma_VarLow);
      TrackCutsDeuteronDCA->SetEtaRange(-Deuteron_Eta_VarLow,
                                        Deuteron_Eta_VarLow);
      TrackCutsAntiDeuteronDCA->SetEtaRange(-Deuteron_Eta_VarLow,
                                            Deuteron_Eta_VarLow);
      TrackCutsDeuteronDCA->SetPtRange(Deuteron_pT_VarHigh, 2.4);
      TrackCutsAntiDeuteronDCA->SetPtRange(Deuteron_pT_VarHigh, 2.4);
      TrackCutsDeuteronDCA->SetNClsTPC(Deuteron_Clusters_VarLow);
      TrackCutsAntiDeuteronDCA->SetNClsTPC(Deuteron_Clusters_VarLow);

      posPions->SetPID(AliPID::kPion, 0.5, Pion_Sigma_VarLow);
      NegPions->SetPID(AliPID::kPion, 0.5, Pion_Sigma_VarLow);
      posPions->SetEtaRange(-Pion_Eta_VarHigh, Pion_Eta_VarHigh);
      NegPions->SetEtaRange(-Pion_Eta_VarHigh, Pion_Eta_VarHigh);
      posPions->SetNClsTPC(Pion_Clusters_VarLow);
      NegPions->SetNClsTPC(Pion_Clusters_VarLow);

      config->SetDeltaEtaMax(DEta_VarLow);

    } else if (suffix == "12") {

      TrackCutsDeuteronDCA->SetEtaRange(-Deuteron_Eta_VarLow,
                                        Deuteron_Eta_VarLow);
      TrackCutsAntiDeuteronDCA->SetEtaRange(-Deuteron_Eta_VarLow,
                                            Deuteron_Eta_VarLow);
      TrackCutsDeuteronDCA->SetNClsTPC(Deuteron_Clusters_VarHigh);
      TrackCutsAntiDeuteronDCA->SetNClsTPC(Deuteron_Clusters_VarHigh);

      posPions->SetPID(AliPID::kPion, 0.5, Pion_Sigma_VarHigh);
      NegPions->SetPID(AliPID::kPion, 0.5, Pion_Sigma_VarHigh);
      posPions->SetEtaRange(-Pion_Eta_VarLow, Pion_Eta_VarLow);
      NegPions->SetEtaRange(-Pion_Eta_VarLow, Pion_Eta_VarLow);
      posPions->SetPtRange(Pion_pT_VarHigh, 4.0);
      NegPions->SetPtRange(Pion_pT_VarHigh, 4.0);
      posPions->SetNClsTPC(Pion_Clusters_VarLow);
      NegPions->SetNClsTPC(Pion_Clusters_VarLow);

      config->SetDeltaPhiMax(DPhi_VarLow);

    } else if (suffix == "13") {

      TrackCutsDeuteronDCA->SetPtRange(Deuteron_pT_VarHigh, 2.4);
      TrackCutsAntiDeuteronDCA->SetPtRange(Deuteron_pT_VarHigh, 2.4);
      TrackCutsDeuteronDCA->SetNClsTPC(Deuteron_Clusters_VarHigh);
      TrackCutsAntiDeuteronDCA->SetNClsTPC(Deuteron_Clusters_VarHigh);

      posPions->SetEtaRange(-Pion_Eta_VarHigh, Pion_Eta_VarHigh);
      NegPions->SetEtaRange(-Pion_Eta_VarHigh, Pion_Eta_VarHigh);
      posPions->SetPtRange(Pion_pT_VarLow, 4.0);
      NegPions->SetPtRange(Pion_pT_VarLow, 4.0);
      posPions->SetNClsTPC(Pion_Clusters_VarLow);
      NegPions->SetNClsTPC(Pion_Clusters_VarLow);

      config->SetDeltaEtaMax(DEta_VarLow);

    } else if (suffix == "14") {

      TrackCutsDeuteronDCA->SetPIDkd(false, false, Deuteron_Sigma_VarHigh,
                                     Deuteron_Sigma_VarHigh);
      TrackCutsAntiDeuteronDCA->SetPIDkd(false, false, Deuteron_Sigma_VarHigh,
                                         Deuteron_Sigma_VarHigh);
      TrackCutsDeuteronDCA->SetEtaRange(-Deuteron_Eta_VarLow,
                                        Deuteron_Eta_VarLow);
      TrackCutsAntiDeuteronDCA->SetEtaRange(-Deuteron_Eta_VarLow,
                                            Deuteron_Eta_VarLow);
      TrackCutsDeuteronDCA->SetNClsTPC(Deuteron_Clusters_VarLow);
      TrackCutsAntiDeuteronDCA->SetNClsTPC(Deuteron_Clusters_VarLow);

      posPions->SetPID(AliPID::kPion, 0.5, Pion_Sigma_VarHigh);
      NegPions->SetPID(AliPID::kPion, 0.5, Pion_Sigma_VarHigh);
      posPions->SetPtRange(Pion_pT_VarHigh, 4.0);
      NegPions->SetPtRange(Pion_pT_VarHigh, 4.0);

      config->SetDeltaEtaMax(DEta_VarLow);

    } else if (suffix == "15") {

      TrackCutsDeuteronDCA->SetEtaRange(-Deuteron_Eta_VarHigh,
                                        Deuteron_Eta_VarHigh);
      TrackCutsAntiDeuteronDCA->SetEtaRange(-Deuteron_Eta_VarHigh,
                                            Deuteron_Eta_VarHigh);
      TrackCutsDeuteronDCA->SetPtRange(Deuteron_pT_VarHigh, 2.4);
      TrackCutsAntiDeuteronDCA->SetPtRange(Deuteron_pT_VarHigh, 2.4);
      TrackCutsDeuteronDCA->SetNClsTPC(Deuteron_Clusters_VarHigh);
      TrackCutsAntiDeuteronDCA->SetNClsTPC(Deuteron_Clusters_VarHigh);

      posPions->SetPID(AliPID::kPion, 0.5, Pion_Sigma_VarLow);
      NegPions->SetPID(AliPID::kPion, 0.5, Pion_Sigma_VarLow);
      posPions->SetPtRange(Pion_pT_VarLow, 4.0);
      NegPions->SetPtRange(Pion_pT_VarLow, 4.0);
      posPions->SetNClsTPC(Pion_Clusters_VarHigh);
      NegPions->SetNClsTPC(Pion_Clusters_VarHigh);

      config->SetDeltaEtaMax(DEta_VarLow);

    } else if (suffix == "16") {

      TrackCutsDeuteronDCA->SetPIDkd(false, false, Deuteron_Sigma_VarLow,
                                     Deuteron_Sigma_VarLow);
      TrackCutsAntiDeuteronDCA->SetPIDkd(false, false, Deuteron_Sigma_VarLow,
                                         Deuteron_Sigma_VarLow);
      TrackCutsDeuteronDCA->SetPtRange(Deuteron_pT_VarLow, 2.4);
      TrackCutsAntiDeuteronDCA->SetPtRange(Deuteron_pT_VarLow, 2.4);
      TrackCutsDeuteronDCA->SetNClsTPC(Deuteron_Clusters_VarLow);
      TrackCutsAntiDeuteronDCA->SetNClsTPC(Deuteron_Clusters_VarLow);

      posPions->SetPID(AliPID::kPion, 0.5, Pion_Sigma_VarLow);
      NegPions->SetPID(AliPID::kPion, 0.5, Pion_Sigma_VarLow);
      posPions->SetEtaRange(-Pion_Eta_VarLow, Pion_Eta_VarLow);
      NegPions->SetEtaRange(-Pion_Eta_VarLow, Pion_Eta_VarLow);
      posPions->SetPtRange(Pion_pT_VarLow, 4.0);
      NegPions->SetPtRange(Pion_pT_VarLow, 4.0);

      config->SetDeltaPhiMax(DPhi_VarHigh);
      config->SetDeltaEtaMax(DEta_VarLow);

    } else if (suffix == "17") {

      TrackCutsDeuteronDCA->SetPIDkd(false, false, Deuteron_Sigma_VarHigh,
                                     Deuteron_Sigma_VarHigh);
      TrackCutsAntiDeuteronDCA->SetPIDkd(false, false, Deuteron_Sigma_VarHigh,
                                         Deuteron_Sigma_VarHigh);
      TrackCutsDeuteronDCA->SetPtRange(Deuteron_pT_VarLow, 2.4);
      TrackCutsAntiDeuteronDCA->SetPtRange(Deuteron_pT_VarLow, 2.4);

      posPions->SetPID(AliPID::kPion, 0.5, Pion_Sigma_VarLow);
      NegPions->SetPID(AliPID::kPion, 0.5, Pion_Sigma_VarLow);
      posPions->SetEtaRange(-Pion_Eta_VarHigh, Pion_Eta_VarHigh);
      NegPions->SetEtaRange(-Pion_Eta_VarHigh, Pion_Eta_VarHigh);
      posPions->SetPtRange(Pion_pT_VarLow, 4.0);
      NegPions->SetPtRange(Pion_pT_VarLow, 4.0);
      posPions->SetNClsTPC(Pion_Clusters_VarHigh);
      NegPions->SetNClsTPC(Pion_Clusters_VarHigh);

      config->SetDeltaPhiMax(DPhi_VarLow);

    } else if (suffix == "18") {

      TrackCutsDeuteronDCA->SetEtaRange(-Deuteron_Eta_VarHigh,
                                        Deuteron_Eta_VarHigh);
      TrackCutsAntiDeuteronDCA->SetEtaRange(-Deuteron_Eta_VarHigh,
                                            Deuteron_Eta_VarHigh);
      TrackCutsDeuteronDCA->SetPtRange(Deuteron_pT_VarLow, 2.4);
      TrackCutsAntiDeuteronDCA->SetPtRange(Deuteron_pT_VarLow, 2.4);

      posPions->SetPID(AliPID::kPion, 0.5, Pion_Sigma_VarHigh);
      NegPions->SetPID(AliPID::kPion, 0.5, Pion_Sigma_VarHigh);
      posPions->SetEtaRange(-Pion_Eta_VarHigh, Pion_Eta_VarHigh);
      NegPions->SetEtaRange(-Pion_Eta_VarHigh, Pion_Eta_VarHigh);
      posPions->SetPtRange(Pion_pT_VarLow, 4.0);
      NegPions->SetPtRange(Pion_pT_VarLow, 4.0);

      config->SetDeltaPhiMax(DPhi_VarHigh);

    } else if (suffix == "19") {

      TrackCutsDeuteronDCA->SetEtaRange(-Deuteron_Eta_VarLow,
                                        Deuteron_Eta_VarLow);
      TrackCutsAntiDeuteronDCA->SetEtaRange(-Deuteron_Eta_VarLow,
                                            Deuteron_Eta_VarLow);
      TrackCutsDeuteronDCA->SetNClsTPC(Deuteron_Clusters_VarHigh);
      TrackCutsAntiDeuteronDCA->SetNClsTPC(Deuteron_Clusters_VarHigh);

      posPions->SetEtaRange(-Pion_Eta_VarLow, Pion_Eta_VarLow);
      NegPions->SetEtaRange(-Pion_Eta_VarLow, Pion_Eta_VarLow);
      posPions->SetPtRange(Pion_pT_VarLow, 4.0);
      NegPions->SetPtRange(Pion_pT_VarLow, 4.0);
      posPions->SetNClsTPC(Pion_Clusters_VarLow);
      NegPions->SetNClsTPC(Pion_Clusters_VarLow);

      config->SetDeltaPhiMax(DPhi_VarHigh);
      config->SetDeltaEtaMax(DEta_VarHigh);

    } else if (suffix == "20") {

      TrackCutsDeuteronDCA->SetPIDkd(false, false, Deuteron_Sigma_VarLow,
                                     Deuteron_Sigma_VarLow);
      TrackCutsAntiDeuteronDCA->SetPIDkd(false, false, Deuteron_Sigma_VarLow,
                                         Deuteron_Sigma_VarLow);
      TrackCutsDeuteronDCA->SetEtaRange(-Deuteron_Eta_VarLow,
                                        Deuteron_Eta_VarLow);
      TrackCutsAntiDeuteronDCA->SetEtaRange(-Deuteron_Eta_VarLow,
                                            Deuteron_Eta_VarLow);
      TrackCutsDeuteronDCA->SetPtRange(Deuteron_pT_VarLow, 2.4);
      TrackCutsAntiDeuteronDCA->SetPtRange(Deuteron_pT_VarLow, 2.4);

      posPions->SetPID(AliPID::kPion, 0.5, Pion_Sigma_VarLow);
      NegPions->SetPID(AliPID::kPion, 0.5, Pion_Sigma_VarLow);
      posPions->SetEtaRange(-Pion_Eta_VarLow, Pion_Eta_VarLow);
      NegPions->SetEtaRange(-Pion_Eta_VarLow, Pion_Eta_VarLow);
      posPions->SetPtRange(Pion_pT_VarHigh, 4.0);
      NegPions->SetPtRange(Pion_pT_VarHigh, 4.0);
      posPions->SetNClsTPC(Pion_Clusters_VarLow);
      NegPions->SetNClsTPC(Pion_Clusters_VarLow);

      config->SetDeltaEtaMax(DEta_VarLow);

    } else if (suffix == "21") {

      TrackCutsDeuteronDCA->SetPIDkd(false, false, Deuteron_Sigma_VarLow,
                                     Deuteron_Sigma_VarLow);
      TrackCutsAntiDeuteronDCA->SetPIDkd(false, false, Deuteron_Sigma_VarLow,
                                         Deuteron_Sigma_VarLow);
      TrackCutsDeuteronDCA->SetEtaRange(-Deuteron_Eta_VarLow,
                                        Deuteron_Eta_VarLow);
      TrackCutsAntiDeuteronDCA->SetEtaRange(-Deuteron_Eta_VarLow,
                                            Deuteron_Eta_VarLow);
      TrackCutsDeuteronDCA->SetNClsTPC(Deuteron_Clusters_VarHigh);
      TrackCutsAntiDeuteronDCA->SetNClsTPC(Deuteron_Clusters_VarHigh);

      posPions->SetPID(AliPID::kPion, 0.5, Pion_Sigma_VarHigh);
      NegPions->SetPID(AliPID::kPion, 0.5, Pion_Sigma_VarHigh);
      posPions->SetPtRange(Pion_pT_VarLow, 4.0);
      NegPions->SetPtRange(Pion_pT_VarLow, 4.0);

      config->SetDeltaPhiMax(DPhi_VarHigh);
      config->SetDeltaEtaMax(DEta_VarLow);

    } else if (suffix == "22") {

      TrackCutsDeuteronDCA->SetPIDkd(false, false, Deuteron_Sigma_VarLow,
                                     Deuteron_Sigma_VarLow);
      TrackCutsAntiDeuteronDCA->SetPIDkd(false, false, Deuteron_Sigma_VarLow,
                                         Deuteron_Sigma_VarLow);
      TrackCutsDeuteronDCA->SetEtaRange(-Deuteron_Eta_VarHigh,
                                        Deuteron_Eta_VarHigh);
      TrackCutsAntiDeuteronDCA->SetEtaRange(-Deuteron_Eta_VarHigh,
                                            Deuteron_Eta_VarHigh);
      TrackCutsDeuteronDCA->SetPtRange(Deuteron_pT_VarLow, 2.4);
      TrackCutsAntiDeuteronDCA->SetPtRange(Deuteron_pT_VarLow, 2.4);
      TrackCutsDeuteronDCA->SetNClsTPC(Deuteron_Clusters_VarHigh);
      TrackCutsAntiDeuteronDCA->SetNClsTPC(Deuteron_Clusters_VarHigh);

      posPions->SetPID(AliPID::kPion, 0.5, Pion_Sigma_VarHigh);
      NegPions->SetPID(AliPID::kPion, 0.5, Pion_Sigma_VarHigh);
      posPions->SetEtaRange(-Pion_Eta_VarHigh, Pion_Eta_VarHigh);
      NegPions->SetEtaRange(-Pion_Eta_VarHigh, Pion_Eta_VarHigh);
      posPions->SetPtRange(Pion_pT_VarLow, 4.0);
      NegPions->SetPtRange(Pion_pT_VarLow, 4.0);

      config->SetDeltaPhiMax(DPhi_VarLow);

    } else if (suffix == "23") {

      TrackCutsDeuteronDCA->SetPIDkd(false, false, Deuteron_Sigma_VarLow,
                                     Deuteron_Sigma_VarLow);
      TrackCutsAntiDeuteronDCA->SetPIDkd(false, false, Deuteron_Sigma_VarLow,
                                         Deuteron_Sigma_VarLow);
      TrackCutsDeuteronDCA->SetEtaRange(-Deuteron_Eta_VarHigh,
                                        Deuteron_Eta_VarHigh);
      TrackCutsAntiDeuteronDCA->SetEtaRange(-Deuteron_Eta_VarHigh,
                                            Deuteron_Eta_VarHigh);
      TrackCutsDeuteronDCA->SetPtRange(Deuteron_pT_VarHigh, 2.4);
      TrackCutsAntiDeuteronDCA->SetPtRange(Deuteron_pT_VarHigh, 2.4);
      TrackCutsDeuteronDCA->SetNClsTPC(Deuteron_Clusters_VarLow);
      TrackCutsAntiDeuteronDCA->SetNClsTPC(Deuteron_Clusters_VarLow);

      posPions->SetEtaRange(-Pion_Eta_VarLow, Pion_Eta_VarLow);
      NegPions->SetEtaRange(-Pion_Eta_VarLow, Pion_Eta_VarLow);
      posPions->SetNClsTPC(Pion_Clusters_VarHigh);
      NegPions->SetNClsTPC(Pion_Clusters_VarHigh);

      config->SetDeltaEtaMax(DEta_VarHigh);

    } else if (suffix == "24") {

      TrackCutsDeuteronDCA->SetPIDkd(false, false, Deuteron_Sigma_VarLow,
                                     Deuteron_Sigma_VarLow);
      TrackCutsAntiDeuteronDCA->SetPIDkd(false, false, Deuteron_Sigma_VarLow,
                                         Deuteron_Sigma_VarLow);
      TrackCutsDeuteronDCA->SetEtaRange(-Deuteron_Eta_VarLow,
                                        Deuteron_Eta_VarLow);
      TrackCutsAntiDeuteronDCA->SetEtaRange(-Deuteron_Eta_VarLow,
                                            Deuteron_Eta_VarLow);
      TrackCutsDeuteronDCA->SetPtRange(Deuteron_pT_VarLow, 2.4);
      TrackCutsAntiDeuteronDCA->SetPtRange(Deuteron_pT_VarLow, 2.4);
      TrackCutsDeuteronDCA->SetNClsTPC(Deuteron_Clusters_VarLow);
      TrackCutsAntiDeuteronDCA->SetNClsTPC(Deuteron_Clusters_VarLow);

      posPions->SetEtaRange(-Pion_Eta_VarHigh, Pion_Eta_VarHigh);
      NegPions->SetEtaRange(-Pion_Eta_VarHigh, Pion_Eta_VarHigh);
      posPions->SetPtRange(Pion_pT_VarHigh, 4.0);
      NegPions->SetPtRange(Pion_pT_VarHigh, 4.0);
      posPions->SetNClsTPC(Pion_Clusters_VarHigh);
      NegPions->SetNClsTPC(Pion_Clusters_VarHigh);

      config->SetDeltaPhiMax(DPhi_VarLow);
      config->SetDeltaEtaMax(DEta_VarHigh);

    } else if (suffix == "25") {

      TrackCutsDeuteronDCA->SetPtRange(Deuteron_pT_VarLow, 2.4);
      TrackCutsAntiDeuteronDCA->SetPtRange(Deuteron_pT_VarLow, 2.4);
      TrackCutsDeuteronDCA->SetNClsTPC(Deuteron_Clusters_VarHigh);
      TrackCutsAntiDeuteronDCA->SetNClsTPC(Deuteron_Clusters_VarHigh);

      posPions->SetEtaRange(-Pion_Eta_VarLow, Pion_Eta_VarLow);
      NegPions->SetEtaRange(-Pion_Eta_VarLow, Pion_Eta_VarLow);
      posPions->SetNClsTPC(Pion_Clusters_VarHigh);
      NegPions->SetNClsTPC(Pion_Clusters_VarHigh);

    } else if (suffix == "26") {

      TrackCutsDeuteronDCA->SetEtaRange(-Deuteron_Eta_VarLow,
                                        Deuteron_Eta_VarLow);
      TrackCutsAntiDeuteronDCA->SetEtaRange(-Deuteron_Eta_VarLow,
                                            Deuteron_Eta_VarLow);

      posPions->SetPID(AliPID::kPion, 0.5, Pion_Sigma_VarLow);
      NegPions->SetPID(AliPID::kPion, 0.5, Pion_Sigma_VarLow);
      posPions->SetEtaRange(-Pion_Eta_VarLow, Pion_Eta_VarLow);
      NegPions->SetEtaRange(-Pion_Eta_VarLow, Pion_Eta_VarLow);
      posPions->SetPtRange(Pion_pT_VarLow, 4.0);
      NegPions->SetPtRange(Pion_pT_VarLow, 4.0);
      posPions->SetNClsTPC(Pion_Clusters_VarHigh);
      NegPions->SetNClsTPC(Pion_Clusters_VarHigh);

      config->SetDeltaPhiMax(DPhi_VarHigh);
      config->SetDeltaEtaMax(DEta_VarLow);

    } else if (suffix == "27") {

      TrackCutsDeuteronDCA->SetEtaRange(-Deuteron_Eta_VarHigh,
                                        Deuteron_Eta_VarHigh);
      TrackCutsAntiDeuteronDCA->SetEtaRange(-Deuteron_Eta_VarHigh,
                                            Deuteron_Eta_VarHigh);
      TrackCutsDeuteronDCA->SetNClsTPC(Deuteron_Clusters_VarHigh);
      TrackCutsAntiDeuteronDCA->SetNClsTPC(Deuteron_Clusters_VarHigh);

      posPions->SetEtaRange(-Pion_Eta_VarHigh, Pion_Eta_VarHigh);
      NegPions->SetEtaRange(-Pion_Eta_VarHigh, Pion_Eta_VarHigh);
      posPions->SetPtRange(Pion_pT_VarHigh, 4.0);
      NegPions->SetPtRange(Pion_pT_VarHigh, 4.0);
      posPions->SetNClsTPC(Pion_Clusters_VarLow);
      NegPions->SetNClsTPC(Pion_Clusters_VarLow);

      config->SetDeltaPhiMax(DPhi_VarHigh);

    } else if (suffix == "28") {

      TrackCutsDeuteronDCA->SetPIDkd(false, false, Deuteron_Sigma_VarLow,
                                     Deuteron_Sigma_VarLow);
      TrackCutsAntiDeuteronDCA->SetPIDkd(false, false, Deuteron_Sigma_VarLow,
                                         Deuteron_Sigma_VarLow);
      TrackCutsDeuteronDCA->SetEtaRange(-Deuteron_Eta_VarHigh,
                                        Deuteron_Eta_VarHigh);
      TrackCutsAntiDeuteronDCA->SetEtaRange(-Deuteron_Eta_VarHigh,
                                            Deuteron_Eta_VarHigh);
      TrackCutsDeuteronDCA->SetNClsTPC(Deuteron_Clusters_VarHigh);
      TrackCutsAntiDeuteronDCA->SetNClsTPC(Deuteron_Clusters_VarHigh);

      posPions->SetPID(AliPID::kPion, 0.5, Pion_Sigma_VarLow);
      NegPions->SetPID(AliPID::kPion, 0.5, Pion_Sigma_VarLow);
      posPions->SetEtaRange(-Pion_Eta_VarLow, Pion_Eta_VarLow);
      NegPions->SetEtaRange(-Pion_Eta_VarLow, Pion_Eta_VarLow);
      posPions->SetPtRange(Pion_pT_VarLow, 4.0);
      NegPions->SetPtRange(Pion_pT_VarLow, 4.0);
      posPions->SetNClsTPC(Pion_Clusters_VarLow);
      NegPions->SetNClsTPC(Pion_Clusters_VarLow);

      config->SetDeltaEtaMax(DEta_VarLow);

    } else if (suffix == "29") {

      TrackCutsDeuteronDCA->SetPtRange(Deuteron_pT_VarLow, 2.4);
      TrackCutsAntiDeuteronDCA->SetPtRange(Deuteron_pT_VarLow, 2.4);
      TrackCutsDeuteronDCA->SetNClsTPC(Deuteron_Clusters_VarHigh);
      TrackCutsAntiDeuteronDCA->SetNClsTPC(Deuteron_Clusters_VarHigh);

      posPions->SetPID(AliPID::kPion, 0.5, Pion_Sigma_VarHigh);
      NegPions->SetPID(AliPID::kPion, 0.5, Pion_Sigma_VarHigh);
      posPions->SetEtaRange(-Pion_Eta_VarHigh, Pion_Eta_VarHigh);
      NegPions->SetEtaRange(-Pion_Eta_VarHigh, Pion_Eta_VarHigh);

      config->SetDeltaPhiMax(DPhi_VarHigh);
      config->SetDeltaEtaMax(DEta_VarLow);

    } else if (suffix == "30") {

      TrackCutsDeuteronDCA->SetPIDkd(false, false, Deuteron_Sigma_VarHigh,
                                     Deuteron_Sigma_VarHigh);
      TrackCutsAntiDeuteronDCA->SetPIDkd(false, false, Deuteron_Sigma_VarHigh,
                                         Deuteron_Sigma_VarHigh);
      TrackCutsDeuteronDCA->SetEtaRange(-Deuteron_Eta_VarHigh,
                                        Deuteron_Eta_VarHigh);
      TrackCutsAntiDeuteronDCA->SetEtaRange(-Deuteron_Eta_VarHigh,
                                            Deuteron_Eta_VarHigh);
      TrackCutsDeuteronDCA->SetPtRange(Deuteron_pT_VarHigh, 2.4);
      TrackCutsAntiDeuteronDCA->SetPtRange(Deuteron_pT_VarHigh, 2.4);
      TrackCutsDeuteronDCA->SetNClsTPC(Deuteron_Clusters_VarLow);
      TrackCutsAntiDeuteronDCA->SetNClsTPC(Deuteron_Clusters_VarLow);

      posPions->SetPID(AliPID::kPion, 0.5, Pion_Sigma_VarLow);
      NegPions->SetPID(AliPID::kPion, 0.5, Pion_Sigma_VarLow);
      posPions->SetPtRange(Pion_pT_VarHigh, 4.0);
      NegPions->SetPtRange(Pion_pT_VarHigh, 4.0);
      posPions->SetNClsTPC(Pion_Clusters_VarLow);
      NegPions->SetNClsTPC(Pion_Clusters_VarLow);

      config->SetDeltaEtaMax(DEta_VarHigh);

    } else if (suffix == "31") {

      TrackCutsDeuteronDCA->SetPIDkd(false, false, Deuteron_Sigma_VarLow,
                                     Deuteron_Sigma_VarLow);
      TrackCutsAntiDeuteronDCA->SetPIDkd(false, false, Deuteron_Sigma_VarLow,
                                         Deuteron_Sigma_VarLow);
      TrackCutsDeuteronDCA->SetEtaRange(-Deuteron_Eta_VarHigh,
                                        Deuteron_Eta_VarHigh);
      TrackCutsAntiDeuteronDCA->SetEtaRange(-Deuteron_Eta_VarHigh,
                                            Deuteron_Eta_VarHigh);
      TrackCutsDeuteronDCA->SetNClsTPC(Deuteron_Clusters_VarHigh);
      TrackCutsAntiDeuteronDCA->SetNClsTPC(Deuteron_Clusters_VarHigh);

      posPions->SetPID(AliPID::kPion, 0.5, Pion_Sigma_VarHigh);
      NegPions->SetPID(AliPID::kPion, 0.5, Pion_Sigma_VarHigh);
      posPions->SetEtaRange(-Pion_Eta_VarLow, Pion_Eta_VarLow);
      NegPions->SetEtaRange(-Pion_Eta_VarLow, Pion_Eta_VarLow);
      posPions->SetPtRange(Pion_pT_VarHigh, 4.0);
      NegPions->SetPtRange(Pion_pT_VarHigh, 4.0);
      posPions->SetNClsTPC(Pion_Clusters_VarLow);
      NegPions->SetNClsTPC(Pion_Clusters_VarLow);

      config->SetDeltaPhiMax(DPhi_VarLow);

    } else if (suffix == "32") {

      TrackCutsDeuteronDCA->SetPIDkd(false, false, Deuteron_Sigma_VarHigh,
                                     Deuteron_Sigma_VarHigh);
      TrackCutsAntiDeuteronDCA->SetPIDkd(false, false, Deuteron_Sigma_VarHigh,
                                         Deuteron_Sigma_VarHigh);
      TrackCutsDeuteronDCA->SetNClsTPC(Deuteron_Clusters_VarLow);
      TrackCutsAntiDeuteronDCA->SetNClsTPC(Deuteron_Clusters_VarLow);

      posPions->SetPtRange(Pion_pT_VarHigh, 4.0);
      NegPions->SetPtRange(Pion_pT_VarHigh, 4.0);

      config->SetDeltaPhiMax(DPhi_VarLow);
      config->SetDeltaEtaMax(DEta_VarHigh);

    } else if (suffix == "33") {

      TrackCutsDeuteronDCA->SetPIDkd(false, false, Deuteron_Sigma_VarLow,
                                     Deuteron_Sigma_VarLow);
      TrackCutsAntiDeuteronDCA->SetPIDkd(false, false, Deuteron_Sigma_VarLow,
                                         Deuteron_Sigma_VarLow);
      TrackCutsDeuteronDCA->SetEtaRange(-Deuteron_Eta_VarLow,
                                        Deuteron_Eta_VarLow);
      TrackCutsAntiDeuteronDCA->SetEtaRange(-Deuteron_Eta_VarLow,
                                            Deuteron_Eta_VarLow);
      TrackCutsDeuteronDCA->SetNClsTPC(Deuteron_Clusters_VarLow);
      TrackCutsAntiDeuteronDCA->SetNClsTPC(Deuteron_Clusters_VarLow);

      posPions->SetEtaRange(-Pion_Eta_VarHigh, Pion_Eta_VarHigh);
      NegPions->SetEtaRange(-Pion_Eta_VarHigh, Pion_Eta_VarHigh);
      posPions->SetPtRange(Pion_pT_VarLow, 4.0);
      NegPions->SetPtRange(Pion_pT_VarLow, 4.0);

      config->SetDeltaEtaMax(DEta_VarHigh);

    } else if (suffix == "34") {

      TrackCutsDeuteronDCA->SetPIDkd(false, false, Deuteron_Sigma_VarLow,
                                     Deuteron_Sigma_VarLow);
      TrackCutsAntiDeuteronDCA->SetPIDkd(false, false, Deuteron_Sigma_VarLow,
                                         Deuteron_Sigma_VarLow);

      posPions->SetPID(AliPID::kPion, 0.5, Pion_Sigma_VarHigh);
      NegPions->SetPID(AliPID::kPion, 0.5, Pion_Sigma_VarHigh);
      posPions->SetEtaRange(-Pion_Eta_VarLow, Pion_Eta_VarLow);
      NegPions->SetEtaRange(-Pion_Eta_VarLow, Pion_Eta_VarLow);
      posPions->SetNClsTPC(Pion_Clusters_VarLow);
      NegPions->SetNClsTPC(Pion_Clusters_VarLow);

      config->SetDeltaEtaMax(DEta_VarHigh);

    } else if (suffix == "35") {

      TrackCutsDeuteronDCA->SetPIDkd(false, false, Deuteron_Sigma_VarLow,
                                     Deuteron_Sigma_VarLow);
      TrackCutsAntiDeuteronDCA->SetPIDkd(false, false, Deuteron_Sigma_VarLow,
                                         Deuteron_Sigma_VarLow);
      TrackCutsDeuteronDCA->SetPtRange(Deuteron_pT_VarHigh, 2.4);
      TrackCutsAntiDeuteronDCA->SetPtRange(Deuteron_pT_VarHigh, 2.4);
      TrackCutsDeuteronDCA->SetNClsTPC(Deuteron_Clusters_VarHigh);
      TrackCutsAntiDeuteronDCA->SetNClsTPC(Deuteron_Clusters_VarHigh);

      posPions->SetPID(AliPID::kPion, 0.5, Pion_Sigma_VarLow);
      NegPions->SetPID(AliPID::kPion, 0.5, Pion_Sigma_VarLow);
      posPions->SetEtaRange(-Pion_Eta_VarLow, Pion_Eta_VarLow);
      NegPions->SetEtaRange(-Pion_Eta_VarLow, Pion_Eta_VarLow);
      posPions->SetPtRange(Pion_pT_VarHigh, 4.0);
      NegPions->SetPtRange(Pion_pT_VarHigh, 4.0);
      posPions->SetNClsTPC(Pion_Clusters_VarLow);
      NegPions->SetNClsTPC(Pion_Clusters_VarLow);

      config->SetDeltaPhiMax(DPhi_VarLow);
      config->SetDeltaEtaMax(DEta_VarHigh);

    } else if (suffix == "36") {

      TrackCutsDeuteronDCA->SetPIDkd(false, false, Deuteron_Sigma_VarHigh,
                                     Deuteron_Sigma_VarHigh);
      TrackCutsAntiDeuteronDCA->SetPIDkd(false, false, Deuteron_Sigma_VarHigh,
                                         Deuteron_Sigma_VarHigh);
      TrackCutsDeuteronDCA->SetEtaRange(-Deuteron_Eta_VarHigh,
                                        Deuteron_Eta_VarHigh);
      TrackCutsAntiDeuteronDCA->SetEtaRange(-Deuteron_Eta_VarHigh,
                                            Deuteron_Eta_VarHigh);
      TrackCutsDeuteronDCA->SetPtRange(Deuteron_pT_VarLow, 2.4);
      TrackCutsAntiDeuteronDCA->SetPtRange(Deuteron_pT_VarLow, 2.4);
      TrackCutsDeuteronDCA->SetNClsTPC(Deuteron_Clusters_VarHigh);
      TrackCutsAntiDeuteronDCA->SetNClsTPC(Deuteron_Clusters_VarHigh);

      posPions->SetPID(AliPID::kPion, 0.5, Pion_Sigma_VarHigh);
      NegPions->SetPID(AliPID::kPion, 0.5, Pion_Sigma_VarHigh);
      posPions->SetEtaRange(-Pion_Eta_VarLow, Pion_Eta_VarLow);
      NegPions->SetEtaRange(-Pion_Eta_VarLow, Pion_Eta_VarLow);
      posPions->SetPtRange(Pion_pT_VarHigh, 4.0);
      NegPions->SetPtRange(Pion_pT_VarHigh, 4.0);
      posPions->SetNClsTPC(Pion_Clusters_VarLow);
      NegPions->SetNClsTPC(Pion_Clusters_VarLow);

      config->SetDeltaPhiMax(DPhi_VarLow);
      config->SetDeltaEtaMax(DEta_VarHigh);

    } else if (suffix == "37") {

      TrackCutsDeuteronDCA->SetPIDkd(false, false, Deuteron_Sigma_VarLow,
                                     Deuteron_Sigma_VarLow);
      TrackCutsAntiDeuteronDCA->SetPIDkd(false, false, Deuteron_Sigma_VarLow,
                                         Deuteron_Sigma_VarLow);
      TrackCutsDeuteronDCA->SetPtRange(Deuteron_pT_VarHigh, 2.4);
      TrackCutsAntiDeuteronDCA->SetPtRange(Deuteron_pT_VarHigh, 2.4);
      TrackCutsDeuteronDCA->SetNClsTPC(Deuteron_Clusters_VarHigh);
      TrackCutsAntiDeuteronDCA->SetNClsTPC(Deuteron_Clusters_VarHigh);

      posPions->SetEtaRange(-Pion_Eta_VarLow, Pion_Eta_VarLow);
      NegPions->SetEtaRange(-Pion_Eta_VarLow, Pion_Eta_VarLow);
      posPions->SetPtRange(Pion_pT_VarLow, 4.0);
      NegPions->SetPtRange(Pion_pT_VarLow, 4.0);

      config->SetDeltaPhiMax(DPhi_VarHigh);

    } else if (suffix == "38") {

      TrackCutsDeuteronDCA->SetPIDkd(false, false, Deuteron_Sigma_VarLow,
                                     Deuteron_Sigma_VarLow);
      TrackCutsAntiDeuteronDCA->SetPIDkd(false, false, Deuteron_Sigma_VarLow,
                                         Deuteron_Sigma_VarLow);
      TrackCutsDeuteronDCA->SetPtRange(Deuteron_pT_VarLow, 2.4);
      TrackCutsAntiDeuteronDCA->SetPtRange(Deuteron_pT_VarLow, 2.4);
      TrackCutsDeuteronDCA->SetNClsTPC(Deuteron_Clusters_VarHigh);
      TrackCutsAntiDeuteronDCA->SetNClsTPC(Deuteron_Clusters_VarHigh);

      posPions->SetEtaRange(-Pion_Eta_VarHigh, Pion_Eta_VarHigh);
      NegPions->SetEtaRange(-Pion_Eta_VarHigh, Pion_Eta_VarHigh);
      posPions->SetPtRange(Pion_pT_VarHigh, 4.0);
      NegPions->SetPtRange(Pion_pT_VarHigh, 4.0);
      posPions->SetNClsTPC(Pion_Clusters_VarHigh);
      NegPions->SetNClsTPC(Pion_Clusters_VarHigh);

      config->SetDeltaPhiMax(DPhi_VarHigh);
      config->SetDeltaEtaMax(DEta_VarLow);

    } else if (suffix == "39") {

      TrackCutsDeuteronDCA->SetPIDkd(false, false, Deuteron_Sigma_VarLow,
                                     Deuteron_Sigma_VarLow);
      TrackCutsAntiDeuteronDCA->SetPIDkd(false, false, Deuteron_Sigma_VarLow,
                                         Deuteron_Sigma_VarLow);
      TrackCutsDeuteronDCA->SetEtaRange(-Deuteron_Eta_VarLow,
                                        Deuteron_Eta_VarLow);
      TrackCutsAntiDeuteronDCA->SetEtaRange(-Deuteron_Eta_VarLow,
                                            Deuteron_Eta_VarLow);
      TrackCutsDeuteronDCA->SetPtRange(Deuteron_pT_VarHigh, 2.4);
      TrackCutsAntiDeuteronDCA->SetPtRange(Deuteron_pT_VarHigh, 2.4);
      TrackCutsDeuteronDCA->SetNClsTPC(Deuteron_Clusters_VarHigh);
      TrackCutsAntiDeuteronDCA->SetNClsTPC(Deuteron_Clusters_VarHigh);

      posPions->SetPID(AliPID::kPion, 0.5, Pion_Sigma_VarHigh);
      NegPions->SetPID(AliPID::kPion, 0.5, Pion_Sigma_VarHigh);
      posPions->SetPtRange(Pion_pT_VarHigh, 4.0);
      NegPions->SetPtRange(Pion_pT_VarHigh, 4.0);
      posPions->SetNClsTPC(Pion_Clusters_VarHigh);
      NegPions->SetNClsTPC(Pion_Clusters_VarHigh);

      config->SetDeltaPhiMax(DPhi_VarLow);
      config->SetDeltaEtaMax(DEta_VarLow);

    } else if (suffix == "40") {

      TrackCutsDeuteronDCA->SetPIDkd(false, false, Deuteron_Sigma_VarLow,
                                     Deuteron_Sigma_VarLow);
      TrackCutsAntiDeuteronDCA->SetPIDkd(false, false, Deuteron_Sigma_VarLow,
                                         Deuteron_Sigma_VarLow);
      TrackCutsDeuteronDCA->SetEtaRange(-Deuteron_Eta_VarHigh,
                                        Deuteron_Eta_VarHigh);
      TrackCutsAntiDeuteronDCA->SetEtaRange(-Deuteron_Eta_VarHigh,
                                            Deuteron_Eta_VarHigh);
      TrackCutsDeuteronDCA->SetPtRange(Deuteron_pT_VarHigh, 2.4);
      TrackCutsAntiDeuteronDCA->SetPtRange(Deuteron_pT_VarHigh, 2.4);
      TrackCutsDeuteronDCA->SetNClsTPC(Deuteron_Clusters_VarLow);
      TrackCutsAntiDeuteronDCA->SetNClsTPC(Deuteron_Clusters_VarLow);

      posPions->SetPID(AliPID::kPion, 0.5, Pion_Sigma_VarHigh);
      NegPions->SetPID(AliPID::kPion, 0.5, Pion_Sigma_VarHigh);
      posPions->SetEtaRange(-Pion_Eta_VarHigh, Pion_Eta_VarHigh);
      NegPions->SetEtaRange(-Pion_Eta_VarHigh, Pion_Eta_VarHigh);
      posPions->SetNClsTPC(Pion_Clusters_VarHigh);
      NegPions->SetNClsTPC(Pion_Clusters_VarHigh);

      config->SetDeltaPhiMax(DPhi_VarHigh);
      config->SetDeltaEtaMax(DEta_VarLow);

    } else if (suffix == "41") {

      TrackCutsDeuteronDCA->SetEtaRange(-Deuteron_Eta_VarLow,
                                        Deuteron_Eta_VarLow);
      TrackCutsAntiDeuteronDCA->SetEtaRange(-Deuteron_Eta_VarLow,
                                            Deuteron_Eta_VarLow);
      TrackCutsDeuteronDCA->SetNClsTPC(Deuteron_Clusters_VarHigh);
      TrackCutsAntiDeuteronDCA->SetNClsTPC(Deuteron_Clusters_VarHigh);

      posPions->SetPID(AliPID::kPion, 0.5, Pion_Sigma_VarLow);
      NegPions->SetPID(AliPID::kPion, 0.5, Pion_Sigma_VarLow);
      posPions->SetEtaRange(-Pion_Eta_VarLow, Pion_Eta_VarLow);
      NegPions->SetEtaRange(-Pion_Eta_VarLow, Pion_Eta_VarLow);

      config->SetDeltaPhiMax(DPhi_VarLow);
      config->SetDeltaEtaMax(DEta_VarHigh);

    } else if (suffix == "42") {

      TrackCutsDeuteronDCA->SetEtaRange(-Deuteron_Eta_VarLow,
                                        Deuteron_Eta_VarLow);
      TrackCutsAntiDeuteronDCA->SetEtaRange(-Deuteron_Eta_VarLow,
                                            Deuteron_Eta_VarLow);
      TrackCutsDeuteronDCA->SetPtRange(Deuteron_pT_VarLow, 2.4);
      TrackCutsAntiDeuteronDCA->SetPtRange(Deuteron_pT_VarLow, 2.4);

      posPions->SetEtaRange(-Pion_Eta_VarLow, Pion_Eta_VarLow);
      NegPions->SetEtaRange(-Pion_Eta_VarLow, Pion_Eta_VarLow);

      config->SetDeltaPhiMax(DPhi_VarHigh);
      config->SetDeltaEtaMax(DEta_VarLow);

    } else if (suffix == "43") {

      TrackCutsDeuteronDCA->SetPIDkd(false, false, Deuteron_Sigma_VarHigh,
                                     Deuteron_Sigma_VarHigh);
      TrackCutsAntiDeuteronDCA->SetPIDkd(false, false, Deuteron_Sigma_VarHigh,
                                         Deuteron_Sigma_VarHigh);
      TrackCutsDeuteronDCA->SetPtRange(Deuteron_pT_VarLow, 2.4);
      TrackCutsAntiDeuteronDCA->SetPtRange(Deuteron_pT_VarLow, 2.4);
      TrackCutsDeuteronDCA->SetNClsTPC(Deuteron_Clusters_VarHigh);
      TrackCutsAntiDeuteronDCA->SetNClsTPC(Deuteron_Clusters_VarHigh);

      posPions->SetPID(AliPID::kPion, 0.5, Pion_Sigma_VarHigh);
      NegPions->SetPID(AliPID::kPion, 0.5, Pion_Sigma_VarHigh);
      posPions->SetPtRange(Pion_pT_VarLow, 4.0);
      NegPions->SetPtRange(Pion_pT_VarLow, 4.0);
      posPions->SetNClsTPC(Pion_Clusters_VarHigh);
      NegPions->SetNClsTPC(Pion_Clusters_VarHigh);

      config->SetDeltaEtaMax(DEta_VarLow);

    } else if (suffix == "44") {

      TrackCutsDeuteronDCA->SetPIDkd(false, false, Deuteron_Sigma_VarLow,
                                     Deuteron_Sigma_VarLow);
      TrackCutsAntiDeuteronDCA->SetPIDkd(false, false, Deuteron_Sigma_VarLow,
                                         Deuteron_Sigma_VarLow);
      TrackCutsDeuteronDCA->SetPtRange(Deuteron_pT_VarLow, 2.4);
      TrackCutsAntiDeuteronDCA->SetPtRange(Deuteron_pT_VarLow, 2.4);
      TrackCutsDeuteronDCA->SetNClsTPC(Deuteron_Clusters_VarLow);
      TrackCutsAntiDeuteronDCA->SetNClsTPC(Deuteron_Clusters_VarLow);

      posPions->SetEtaRange(-Pion_Eta_VarLow, Pion_Eta_VarLow);
      NegPions->SetEtaRange(-Pion_Eta_VarLow, Pion_Eta_VarLow);
      posPions->SetPtRange(Pion_pT_VarLow, 4.0);
      NegPions->SetPtRange(Pion_pT_VarLow, 4.0);
      posPions->SetNClsTPC(Pion_Clusters_VarLow);
      NegPions->SetNClsTPC(Pion_Clusters_VarLow);

      config->SetDeltaEtaMax(DEta_VarHigh);
    }
  }

  AliAnalysisTaskFemtoDreamDeuteron *task =
      new AliAnalysisTaskFemtoDreamDeuteron("FemtoDreamDefault", isMC);
  if (trigger == "kINT7") {
    task->SelectCollisionCandidates(AliVEvent::kINT7);
    std::cout << "Added kINT7 Trigger \n";
  } else if (trigger == "kHM") {
    task->SelectCollisionCandidates(AliVEvent::kHighMultV0);
    std::cout << "Added kHighMult Trigger \n";
  } else {
    std::cout << "============================================================="
                 "========"
              << std::endl;
    std::cout << "============================================================="
                 "========"
              << std::endl;
    std::cout << "Centrality Estimator not set, fix it else your Results will "
                 "be empty!"
              << std::endl;
    std::cout << "============================================================="
                 "========"
              << std::endl;
    std::cout << "============================================================="
                 "========"
              << std::endl;
  }
  if (!fullBlastQA) {
    task->SetRunTaskLightWeight(true);
  }
  task->SetEventCuts(evtCuts);
  task->SetTrackCutsProtonDCA(posPions);
  task->SetTrackCutsAntiProtonDCA(NegPions);
  task->SetTrackCutsDeuteronDCA(TrackCutsDeuteronDCA);
  task->SetTrackCutsAntiDeuteronDCA(TrackCutsAntiDeuteronDCA);
  task->SetTrackCutsDeuteronMass(TrackCutsDeuteronMass);
  task->SetTrackCutsAntiDeuteronMass(TrackCutsAntiDeuteronMass);

  task->SetCollectionConfig(config);

  mgr->AddTask(task);

  TString addon = "";

  if (trigger == "kINT7") {
    addon += "MB";
  } else if (trigger == "kHM") {
    addon += "HM";
  }

  TString file = AliAnalysisManager::GetCommonFileName();
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  mgr->ConnectInput(task, 0, cinput);

  TString EvtCutsName = Form("%sEvtCuts%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *coutputEvtCuts = mgr->CreateContainer(
      EvtCutsName.Data(), TList::Class(), AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), EvtCutsName.Data()));
  mgr->ConnectOutput(task, 1, coutputEvtCuts);

  TString TrackCutsName = Form("%sPosPions%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *coutputTrkCuts =
      mgr->CreateContainer(TrackCutsName.Data(), TList::Class(),
                           AliAnalysisManager::kOutputContainer,
                           Form("%s:%s", file.Data(), TrackCutsName.Data()));
  mgr->ConnectOutput(task, 2, coutputTrkCuts);

  TString AntiTrackCutsName = Form("%sNegPions%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *coutputAntiTrkCuts = mgr->CreateContainer(
      AntiTrackCutsName.Data(), TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), AntiTrackCutsName.Data()));
  mgr->ConnectOutput(task, 3, coutputAntiTrkCuts);

  TString TrackCutsDeuteronName =
      Form("%sDeuteronDCA%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *coutputTrkCutsDeuteron = mgr->CreateContainer(
      TrackCutsDeuteronName.Data(), TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), TrackCutsDeuteronName.Data()));
  mgr->ConnectOutput(task, 4, coutputTrkCutsDeuteron);

  TString AntiTrackCutsDeuteronName =
      Form("%sAntiDeuteronDCA%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *coutputAntiTrkCutsDeuteron = mgr->CreateContainer(
      AntiTrackCutsDeuteronName.Data(), TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), AntiTrackCutsDeuteronName.Data()));
  mgr->ConnectOutput(task, 5, coutputAntiTrkCutsDeuteron);

  TString TrackCutsDeuteronNoTOFName =
      Form("%sDeuteronMass%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *coutputTrkCutsDeuteronNoTOF = mgr->CreateContainer(
      TrackCutsDeuteronNoTOFName.Data(), TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), TrackCutsDeuteronNoTOFName.Data()));
  mgr->ConnectOutput(task, 6, coutputTrkCutsDeuteronNoTOF);

  TString AntiTrackCutsDeuteronNoTOFName =
      Form("%sAntiDeuteronMass%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *coutputAntiTrkCutsDeuteronNoTOF =
      mgr->CreateContainer(
          AntiTrackCutsDeuteronNoTOFName.Data(), TList::Class(),
          AliAnalysisManager::kOutputContainer,
          Form("%s:%s", file.Data(), AntiTrackCutsDeuteronNoTOFName.Data()));
  mgr->ConnectOutput(task, 7, coutputAntiTrkCutsDeuteronNoTOF);

  AliAnalysisDataContainer *coutputResults;
  TString ResultsName = Form("%sResults%s", addon.Data(), suffix.Data());
  coutputResults = mgr->CreateContainer(
      ResultsName.Data(), TList::Class(), AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), ResultsName.Data()));
  mgr->ConnectOutput(task, 8, coutputResults);

  AliAnalysisDataContainer *coutputResultsQA;
  TString ResultsQAName = Form("%sResultsQA%s", addon.Data(), suffix.Data());
  coutputResultsQA = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      ResultsQAName.Data(), TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), ResultsQAName.Data()));
  mgr->ConnectOutput(task, 9, coutputResultsQA);

  if (isMC) {
    AliAnalysisDataContainer *coutputTrkCutsMC;
    TString TrkCutsMCName = Form("%sPosPionsMC%s", addon.Data(), suffix.Data());
    coutputTrkCutsMC =
        mgr->CreateContainer(TrkCutsMCName.Data(), TList::Class(),
                             AliAnalysisManager::kOutputContainer,
                             Form("%s:%s", file.Data(), TrkCutsMCName.Data()));
    mgr->ConnectOutput(task, 10, coutputTrkCutsMC);

    AliAnalysisDataContainer *coutputAntiTrkCutsMC;
    TString AntiTrkCutsMCName =
        Form("%sNegPionsMC%s", addon.Data(), suffix.Data());
    coutputAntiTrkCutsMC = mgr->CreateContainer(
        //@suppress("Invalid arguments") it works ffs
        AntiTrkCutsMCName.Data(), TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), AntiTrkCutsMCName.Data()));
    mgr->ConnectOutput(task, 11, coutputAntiTrkCutsMC);

    AliAnalysisDataContainer *coutputv0CutsMC;
    TString v0CutsMCName = Form("%sDeuteronMC%s", addon.Data(), suffix.Data());
    coutputv0CutsMC = mgr->CreateContainer(
        //@suppress("Invalid arguments") it works ffs
        v0CutsMCName.Data(), TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), v0CutsMCName.Data()));
    mgr->ConnectOutput(task, 12, coutputv0CutsMC);

    AliAnalysisDataContainer *coutputAntiv0CutsMC;
    TString Antiv0CutsMCName =
        Form("%sAntiDeuteronMC%s", addon.Data(), suffix.Data());
    coutputAntiv0CutsMC = mgr->CreateContainer(
        //@suppress("Invalid arguments") it works ffs
        Antiv0CutsMCName.Data(), TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), Antiv0CutsMCName.Data()));
    mgr->ConnectOutput(task, 13, coutputAntiv0CutsMC);
  }
  return task;
}
