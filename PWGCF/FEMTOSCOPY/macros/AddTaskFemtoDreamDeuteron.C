#include "TROOT.h"
#include "TSystem.h"

AliAnalysisTaskSE* AddTaskFemtoDreamDeuteron(bool isMC = false,//1
    bool fullBlastQA = true,//2
    bool SystematicLowpT = false,//3
    bool SidebandStudy = false, //4
    bool SystematicFullpT = false,//5
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

  AliFemtoDreamEventCuts *evtCuts = AliFemtoDreamEventCuts::StandardCutsRun2();
  evtCuts->CleanUpMult(false, false, false, true);
  AliFemtoDreamTrackCuts *TrackCutsDeuteronDCA = AliFemtoDreamTrackCuts::PrimDeuteronCuts(
        isMC, true, false, false);
  TrackCutsDeuteronDCA->SetCutCharge(1);

  AliFemtoDreamTrackCuts *TrackCutsDeuteronMass =  AliFemtoDreamTrackCuts::PrimDeuteronCuts(
        isMC, true, false, false);
  TrackCutsDeuteronMass->SetCutCharge(1);
  TrackCutsDeuteronMass->SetPID(AliPID::kDeuteron, 999.0);

  AliFemtoDreamTrackCuts *TrackCutsAntiDeuteronDCA = AliFemtoDreamTrackCuts::PrimDeuteronCuts(
        isMC, true, false, false);
  TrackCutsAntiDeuteronDCA->SetCutCharge(-1);

  AliFemtoDreamTrackCuts *TrackCutsAntiDeuteronMass =  AliFemtoDreamTrackCuts::PrimDeuteronCuts(
        isMC, true, false, false);
  TrackCutsAntiDeuteronMass->SetCutCharge(-1);
  TrackCutsAntiDeuteronMass->SetPID(AliPID::kDeuteron, 999.0);

  AliFemtoDreamTrackCuts *TrackCutsProtonDCA = AliFemtoDreamTrackCuts::PrimProtonCuts(
        isMC, true, false, false);
  TrackCutsProtonDCA->SetCutCharge(1);

  AliFemtoDreamTrackCuts *TrackCutsAntiProtonDCA = AliFemtoDreamTrackCuts::PrimProtonCuts(
        isMC, true, false, false);
  TrackCutsAntiProtonDCA->SetCutCharge(-1);

  std::vector<int> PDGParticles;
  PDGParticles.push_back(2212);
  PDGParticles.push_back(2212);
  PDGParticles.push_back(1000010020);
  PDGParticles.push_back(1000010020);
  std::vector<bool> closeRejection;
 // std::vector<float> mTBins;
//  mTBins.push_back(1.14);
//  mTBins.push_back(1.26);
//  mTBins.push_back(999.);
  std::vector<int> pairQA;
  //pairs:
  // pp             0
  // p bar p        1
  // p d            2
  // p bar d        3
  // bar p bar p    4
  // bar p d        5
  // bar p bar d    6
  // d d            7
  // d bar d        8
  // bar d bar d    9
  const int nPairs = 10;
  for (int i = 0; i < nPairs; ++i) {
    closeRejection.push_back(false);
    pairQA.push_back(0);
  }
  closeRejection[0] = true;  // pp
  closeRejection[2] = true;  // pd
  closeRejection[4] = true;  // barp barp
  closeRejection[6] = true;  // barp bard
  closeRejection[7] = true;  // dd
  closeRejection[9] = true;  // bard bar
  pairQA[0] = 11;    // pp
  pairQA[2] = 11;    // pd
  pairQA[4] = 11;    // barp barp
  pairQA[6] = 11;    // barp bard
  pairQA[7] = 11;    // dd
  pairQA[9] = 11;    // bard bard
  //We need to set the ZVtx bins
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
  //The Multiplicity bins are set here
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
  //minimum k* value
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
  //maximum k* value
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

  AliFemtoDreamCollConfig *config = new AliFemtoDreamCollConfig("Femto", "Femto", false);
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
  config->SetDeltaPhiMax(0.017); // and here you set the actual values
  config->SetMixingDepth(10);
  //config->SetmTBins(mTBins);
  //config->SetDomTMultBinning(true);
  //config->SetmTBinning(true);
  config->SetMultiplicityEstimator(AliFemtoDreamEvent::kRef08);

  if (isMC) {
    config->SetMomentumResolution(true);
  }
  if (fullBlastQA) {
   // config->SetkTBinning(true);
    config->SetPtQA(true);
  }

  if (!fullBlastQA) {
    evtCuts->SetMinimalBooking(true);
    TrackCutsProtonDCA->SetMinimalBooking(true);
    TrackCutsAntiProtonDCA->SetMinimalBooking(true);
    TrackCutsDeuteronDCA->SetMinimalBooking(true);
    TrackCutsDeuteronMass->SetMinimalBooking(true);
    TrackCutsAntiDeuteronDCA->SetMinimalBooking(true);
    TrackCutsAntiDeuteronMass->SetMinimalBooking(true);
    config->SetMinimalBookingME(true);
    config->SetMinimalBookingSample(true);
  }


 if (SystematicLowpT) {
   TrackCutsDeuteronDCA->SetPtRange(0.5, 1.4);
   TrackCutsAntiDeuteronDCA->SetPtRange(0.5, 1.4);
   if (suffix == "0"){
     TrackCutsDeuteronDCA->SetPtRange(0.5, 1.4);
     TrackCutsAntiDeuteronDCA->SetPtRange(0.5, 1.4);
   }else if (suffix == "1") {
      TrackCutsProtonDCA->SetPID(AliPID::kProton, 0.75, 2.5);
      TrackCutsAntiProtonDCA->SetPID(AliPID::kProton, 0.75, 2.5);
      TrackCutsProtonDCA->SetNClsTPC(90);
      TrackCutsAntiProtonDCA->SetNClsTPC(90);
      TrackCutsDeuteronDCA->SetPID(AliPID::kDeuteron,1.4, 2.7);
      TrackCutsAntiDeuteronDCA->SetPID(AliPID::kDeuteron, 1.4, 2.7);
      TrackCutsDeuteronDCA->SetNClsTPC(90);
      TrackCutsAntiDeuteronDCA->SetNClsTPC(90);
    } else if (suffix == "2") {
      TrackCutsProtonDCA->SetPtRange(0.6, 4.05);
      TrackCutsAntiProtonDCA->SetPtRange(0.6, 4.05);
      TrackCutsProtonDCA->SetEtaRange(-0.77, 0.77);
      TrackCutsAntiProtonDCA->SetEtaRange(-0.77, 0.77);
      TrackCutsDeuteronDCA->SetEtaRange(-0.77, 0.77);
      TrackCutsAntiDeuteronDCA->SetEtaRange(-0.77, 0.77);
      config->SetDeltaEtaMax(0.019);
      config->SetDeltaPhiMax(0.019);
    }  else if (suffix == "3") {
      TrackCutsProtonDCA->SetEtaRange(-0.77, 0.77);
      TrackCutsAntiProtonDCA->SetEtaRange(-0.77, 0.77);
      TrackCutsProtonDCA->SetNClsTPC(70);
      TrackCutsAntiProtonDCA->SetNClsTPC(70);
      TrackCutsDeuteronDCA->SetEtaRange(-0.77, 0.77);
      TrackCutsAntiDeuteronDCA->SetEtaRange(-0.77, 0.77);
      TrackCutsDeuteronDCA->SetNClsTPC(70);
      TrackCutsAntiDeuteronDCA->SetNClsTPC(70);
    } else if (suffix == "4") {
      TrackCutsProtonDCA->SetPtRange(0.4, 4.05);
      TrackCutsAntiProtonDCA->SetPtRange(0.4, 4.05);
      TrackCutsProtonDCA->SetEtaRange(-0.77, 0.77);
      TrackCutsAntiProtonDCA->SetEtaRange(-0.77, 0.77);
      TrackCutsProtonDCA->SetPID(AliPID::kProton, 0.75, 3.5);
      TrackCutsAntiProtonDCA->SetPID(AliPID::kProton, 0.75, 3.5);
      TrackCutsDeuteronDCA->SetPtRange(0.4, 1.4);
      TrackCutsAntiProtonDCA->SetPtRange(0.4, 1.4);
      TrackCutsDeuteronDCA->SetEtaRange(-0.77, 0.77);
      TrackCutsAntiProtonDCA->SetEtaRange(-0.77, 0.77);
      TrackCutsDeuteronDCA->SetPID(AliPID::kDeuteron,1.4, 3.3);
      TrackCutsAntiDeuteronDCA->SetPID(AliPID::kDeuteron, 1.4, 3.3);
      config->SetDeltaEtaMax(0.019);
      config->SetDeltaPhiMax(0.019);
    } else if (suffix == "5") {
      TrackCutsProtonDCA->SetPID(AliPID::kProton, 0.75, 2.5);
      TrackCutsAntiProtonDCA->SetPID(AliPID::kProton, 0.75, 2.5);
      TrackCutsProtonDCA->SetNClsTPC(90);
      TrackCutsAntiProtonDCA->SetNClsTPC(90);
      TrackCutsDeuteronDCA->SetPID(AliPID::kDeuteron,1.4, 2.7);
      TrackCutsAntiDeuteronDCA->SetPID(AliPID::kDeuteron, 1.4, 2.7);
      TrackCutsDeuteronDCA->SetNClsTPC(90);
      TrackCutsAntiDeuteronDCA->SetNClsTPC(90);
    } else if (suffix == "6") {
      TrackCutsProtonDCA->SetEtaRange(-0.83, 0.83);
      TrackCutsAntiProtonDCA->SetEtaRange(-0.83, 0.83);
      TrackCutsProtonDCA->SetPID(AliPID::kProton, 0.75, 3.5);
      TrackCutsAntiProtonDCA->SetPID(AliPID::kProton, 0.75, 3.5);
      TrackCutsDeuteronDCA->SetEtaRange(-0.83, 0.83);
      TrackCutsAntiDeuteronDCA->SetEtaRange(-0.83, 0.83);
      TrackCutsDeuteronDCA->SetPID(AliPID::kDeuteron,1.4, 3.3);
      TrackCutsAntiDeuteronDCA->SetPID(AliPID::kDeuteron, 1.4, 3.3);
      config->SetDeltaEtaMax(0.015);
      config->SetDeltaPhiMax(0.015);
    } else if (suffix == "7") {
      TrackCutsProtonDCA->SetPtRange(0.4, 4.05);
      TrackCutsAntiProtonDCA->SetPtRange(0.4, 4.05);
      TrackCutsProtonDCA->SetEtaRange(-0.77, 0.77);
      TrackCutsAntiProtonDCA->SetEtaRange(-0.77, 0.77);
      TrackCutsProtonDCA->SetNClsTPC(90);
      TrackCutsAntiProtonDCA->SetNClsTPC(90);
      TrackCutsDeuteronDCA->SetPtRange(0.4, 1.4);
      TrackCutsAntiDeuteronDCA->SetPtRange(0.4, 1.4);
      TrackCutsDeuteronDCA->SetEtaRange(-0.77, 0.77);
      TrackCutsAntiDeuteronDCA->SetEtaRange(-0.77, 0.77);
      TrackCutsDeuteronDCA->SetNClsTPC(90);
      TrackCutsAntiDeuteronDCA->SetNClsTPC(90);
      config->SetDeltaEtaMax(0.019);
      config->SetDeltaPhiMax(0.019);
    } else if (suffix == "8") {
      TrackCutsProtonDCA->SetPtRange(0.4, 4.05);
      TrackCutsAntiProtonDCA->SetPtRange(0.4, 4.05);
      TrackCutsProtonDCA->SetEtaRange(-0.83, 0.83);
      TrackCutsAntiProtonDCA->SetEtaRange(-0.83, 0.83);
      TrackCutsDeuteronDCA->SetPtRange(0.4, 1.4);
      TrackCutsAntiDeuteronDCA->SetPtRange(0.4, 1.4);
      TrackCutsDeuteronDCA->SetEtaRange(-0.83, 0.83);
      TrackCutsAntiDeuteronDCA->SetEtaRange(-0.83, 0.83);
    } else if (suffix == "9") {
      TrackCutsProtonDCA->SetEtaRange(-0.77, 0.77);
      TrackCutsAntiProtonDCA->SetEtaRange(-0.77, 0.77);
      TrackCutsProtonDCA->SetPID(AliPID::kProton, 0.75, 2.5);
      TrackCutsAntiProtonDCA->SetPID(AliPID::kProton, 0.75, 2.5);
      TrackCutsDeuteronDCA->SetEtaRange(-0.77, 0.77);
      TrackCutsAntiDeuteronDCA->SetEtaRange(-0.77, 0.77);
      TrackCutsDeuteronDCA->SetPID(AliPID::kDeuteron,1.4, 2.7);
      TrackCutsAntiDeuteronDCA->SetPID(AliPID::kDeuteron, 1.4, 2.7);
    } else if (suffix == "10") {
      TrackCutsProtonDCA->SetPtRange(0.4, 4.05);
      TrackCutsAntiProtonDCA->SetPtRange(0.4, 4.05);
      TrackCutsProtonDCA->SetEtaRange(-0.83, 0.83);
      TrackCutsAntiProtonDCA->SetEtaRange(-0.83, 0.83);
      TrackCutsProtonDCA->SetPID(AliPID::kProton, 0.75, 2.5);
      TrackCutsAntiProtonDCA->SetPID(AliPID::kProton, 0.75, 2.5);
      TrackCutsDeuteronDCA->SetPtRange(0.4, 1.4);
      TrackCutsAntiDeuteronDCA->SetPtRange(0.4, 1.4);
      TrackCutsDeuteronDCA->SetEtaRange(-0.83, 0.83);
      TrackCutsAntiDeuteronDCA->SetEtaRange(-0.83, 0.83);
      TrackCutsDeuteronDCA->SetPID(AliPID::kDeuteron,1.4, 2.7);
      TrackCutsAntiDeuteronDCA->SetPID(AliPID::kDeuteron, 1.4, 2.7);
      config->SetDeltaEtaMax(0.019);
      config->SetDeltaPhiMax(0.019);
    } else if (suffix == "11") {
      TrackCutsProtonDCA->SetEtaRange(-0.83, 0.83);
      TrackCutsAntiProtonDCA->SetEtaRange(-0.83, 0.83);
      TrackCutsProtonDCA->SetPID(AliPID::kProton, 0.75, 2.5);
      TrackCutsAntiProtonDCA->SetPID(AliPID::kProton, 0.75, 2.5);
      TrackCutsProtonDCA->SetNClsTPC(90);
      TrackCutsAntiProtonDCA->SetNClsTPC(90);
      TrackCutsDeuteronDCA->SetEtaRange(-0.83, 0.83);
      TrackCutsAntiDeuteronDCA->SetEtaRange(-0.83, 0.83);
      TrackCutsDeuteronDCA->SetPID(AliPID::kDeuteron,1.4, 2.7);
      TrackCutsAntiDeuteronDCA->SetPID(AliPID::kDeuteron, 1.4, 2.7);
      TrackCutsDeuteronDCA->SetNClsTPC(90);
      TrackCutsAntiDeuteronDCA->SetNClsTPC(90);
    } else if (suffix == "12") {
      TrackCutsProtonDCA->SetPtRange(0.4, 4.05);
      TrackCutsAntiProtonDCA->SetPtRange(0.4, 4.05);
      TrackCutsProtonDCA->SetPID(AliPID::kProton, 0.75, 2.5);
      TrackCutsAntiProtonDCA->SetPID(AliPID::kProton, 0.75, 2.5);
      TrackCutsDeuteronDCA->SetPtRange(0.4, 1.4);
      TrackCutsAntiDeuteronDCA->SetPtRange(0.4, 1.4);
      TrackCutsDeuteronDCA->SetPID(AliPID::kDeuteron,1.4, 2.7);
      TrackCutsAntiDeuteronDCA->SetPID(AliPID::kDeuteron, 1.4, 2.7);
      config->SetDeltaEtaMax(0.019);
      config->SetDeltaPhiMax(0.019);
    } else if (suffix == "13") {
      TrackCutsProtonDCA->SetPtRange(0.6, 4.05);
      TrackCutsAntiProtonDCA->SetPtRange(0.6, 4.05);
      TrackCutsProtonDCA->SetEtaRange(-0.83, 0.83);
      TrackCutsAntiProtonDCA->SetEtaRange(-0.83, 0.83);
      TrackCutsProtonDCA->SetPID(AliPID::kProton, 0.75, 2.5);
      TrackCutsAntiProtonDCA->SetPID(AliPID::kProton, 0.75, 2.5);
      TrackCutsProtonDCA->SetNClsTPC(90);
      TrackCutsAntiProtonDCA->SetNClsTPC(90);
      TrackCutsDeuteronDCA->SetPtRange(0.6, 1.4);
      TrackCutsAntiDeuteronDCA->SetPtRange(0.6, 1.4);
      TrackCutsDeuteronDCA->SetEtaRange(-0.83, 0.83);
      TrackCutsAntiDeuteronDCA->SetEtaRange(-0.83, 0.83);
      TrackCutsDeuteronDCA->SetPID(AliPID::kDeuteron,1.4, 2.7);
      TrackCutsAntiDeuteronDCA->SetPID(AliPID::kDeuteron, 1.4, 2.7);
      TrackCutsDeuteronDCA->SetNClsTPC(90);
      TrackCutsAntiDeuteronDCA->SetNClsTPC(90);
    } else if (suffix == "14") {
      TrackCutsProtonDCA->SetPtRange(0.4, 4.05);
      TrackCutsAntiProtonDCA->SetPtRange(0.4, 4.05);
      TrackCutsProtonDCA->SetPID(AliPID::kProton, 0.75, 3.5);
      TrackCutsAntiProtonDCA->SetPID(AliPID::kProton, 0.75, 3.5);
      TrackCutsProtonDCA->SetNClsTPC(70);
      TrackCutsAntiProtonDCA->SetNClsTPC(70);
      TrackCutsDeuteronDCA->SetPtRange(0.4, 1.4);
      TrackCutsAntiDeuteronDCA->SetPtRange(0.4, 1.4);
      TrackCutsDeuteronDCA->SetPID(AliPID::kDeuteron,1.4, 3.3);
      TrackCutsAntiDeuteronDCA->SetPID(AliPID::kDeuteron, 1.4, 3.3);
      TrackCutsDeuteronDCA->SetNClsTPC(70);
      TrackCutsAntiDeuteronDCA->SetNClsTPC(70);
      config->SetDeltaEtaMax(0.019);
      config->SetDeltaPhiMax(0.019);
    } else if (suffix == "15") {
      TrackCutsProtonDCA->SetPtRange(0.4, 4.05);
      TrackCutsAntiProtonDCA->SetPtRange(0.4, 4.05);
      TrackCutsProtonDCA->SetEtaRange(-0.83, 0.83);
      TrackCutsAntiProtonDCA->SetEtaRange(-0.83, 0.83);
      TrackCutsProtonDCA->SetPID(AliPID::kProton, 0.75, 2.5);
      TrackCutsAntiProtonDCA->SetPID(AliPID::kProton, 0.75, 2.5);
      TrackCutsDeuteronDCA->SetPtRange(0.4, 1.4);
      TrackCutsAntiDeuteronDCA->SetPtRange(0.4, 1.4);
      TrackCutsDeuteronDCA->SetEtaRange(-0.83, 0.83);
      TrackCutsAntiDeuteronDCA->SetEtaRange(-0.83, 0.83);
      TrackCutsDeuteronDCA->SetPID(AliPID::kDeuteron,1.4, 2.7);
      TrackCutsAntiDeuteronDCA->SetPID(AliPID::kDeuteron, 1.4, 2.7);
      config->SetDeltaEtaMax(0.015);
      config->SetDeltaPhiMax(0.015);
    } else if (suffix == "16") {
      TrackCutsProtonDCA->SetPtRange(0.4, 4.05);
      TrackCutsAntiProtonDCA->SetPtRange(0.4, 4.05);
      TrackCutsProtonDCA->SetEtaRange(-0.83, 0.83);
      TrackCutsAntiProtonDCA->SetEtaRange(-0.83, 0.83);
      TrackCutsProtonDCA->SetPID(AliPID::kProton, 0.75, 2.5);
      TrackCutsAntiProtonDCA->SetPID(AliPID::kProton, 0.75, 2.5);
      TrackCutsProtonDCA->SetNClsTPC(90);
      TrackCutsAntiProtonDCA->SetNClsTPC(90);
      TrackCutsDeuteronDCA->SetPtRange(0.4, 1.4);
      TrackCutsAntiDeuteronDCA->SetPtRange(0.4, 1.4);
      TrackCutsDeuteronDCA->SetEtaRange(-0.83, 0.83);
      TrackCutsAntiDeuteronDCA->SetEtaRange(-0.83, 0.83);
      TrackCutsDeuteronDCA->SetPID(AliPID::kDeuteron,1.4, 2.7);
      TrackCutsAntiDeuteronDCA->SetPID(AliPID::kDeuteron, 1.4, 2.7);
      TrackCutsDeuteronDCA->SetNClsTPC(90);
      TrackCutsAntiDeuteronDCA->SetNClsTPC(90);
    } else if (suffix == "17") {
      TrackCutsProtonDCA->SetPID(AliPID::kProton, 0.75, 2.5);
      TrackCutsAntiProtonDCA->SetPID(AliPID::kProton, 0.75, 2.5);
      TrackCutsProtonDCA->SetNClsTPC(70);
      TrackCutsAntiProtonDCA->SetNClsTPC(70);
      TrackCutsDeuteronDCA->SetPID(AliPID::kDeuteron,1.4, 2.7);
      TrackCutsAntiDeuteronDCA->SetPID(AliPID::kDeuteron, 1.4, 2.7);
      TrackCutsDeuteronDCA->SetNClsTPC(70);
      TrackCutsAntiDeuteronDCA->SetNClsTPC(70);
      config->SetDeltaEtaMax(0.019);
      config->SetDeltaPhiMax(0.019);
    } else if (suffix == "18") {
      TrackCutsProtonDCA->SetPtRange(0.4, 4.05);
      TrackCutsAntiProtonDCA->SetPtRange(0.4, 4.05);
      TrackCutsProtonDCA->SetNClsTPC(90);
      TrackCutsAntiProtonDCA->SetNClsTPC(90);
      TrackCutsDeuteronDCA->SetPtRange(0.4, 1.4);
      TrackCutsAntiDeuteronDCA->SetPtRange(0.4, 1.4);
      TrackCutsDeuteronDCA->SetNClsTPC(90);
      TrackCutsAntiDeuteronDCA->SetNClsTPC(90);
      config->SetDeltaEtaMax(0.019);
      config->SetDeltaPhiMax(0.019);
    } else if (suffix == "19") {
      TrackCutsProtonDCA->SetEtaRange(-0.83, 0.83);
      TrackCutsAntiProtonDCA->SetEtaRange(-0.83, 0.83);
      TrackCutsProtonDCA->SetPID(AliPID::kProton, 0.75, 3.5);
      TrackCutsAntiProtonDCA->SetPID(AliPID::kProton, 0.75, 3.5);
      TrackCutsProtonDCA->SetNClsTPC(70);
      TrackCutsAntiProtonDCA->SetNClsTPC(70);
      TrackCutsDeuteronDCA->SetEtaRange(-0.83, 0.83);
      TrackCutsAntiDeuteronDCA->SetEtaRange(-0.83, 0.83);
      TrackCutsDeuteronDCA->SetPID(AliPID::kDeuteron,1.4, 3.3);
      TrackCutsAntiDeuteronDCA->SetPID(AliPID::kDeuteron, 1.4, 3.3);
      TrackCutsDeuteronDCA->SetNClsTPC(70);
      TrackCutsAntiDeuteronDCA->SetNClsTPC(70);
      config->SetDeltaEtaMax(0.019);
      config->SetDeltaPhiMax(0.019);
    } else if (suffix == "20") {
      TrackCutsProtonDCA->SetPtRange(0.4, 4.05);
      TrackCutsAntiProtonDCA->SetPtRange(0.4, 4.05);
      TrackCutsProtonDCA->SetEtaRange(-0.77, 0.77);
      TrackCutsAntiProtonDCA->SetEtaRange(-0.77, 0.77);
      TrackCutsProtonDCA->SetNClsTPC(70);
      TrackCutsAntiProtonDCA->SetNClsTPC(70);
      TrackCutsDeuteronDCA->SetPtRange(0.4, 1.4);
      TrackCutsAntiDeuteronDCA->SetPtRange(0.4, 1.4);
      TrackCutsDeuteronDCA->SetEtaRange(-0.77, 0.77);
      TrackCutsAntiDeuteronDCA->SetEtaRange(-0.77, 0.77);
      TrackCutsDeuteronDCA->SetNClsTPC(70);
      TrackCutsAntiDeuteronDCA->SetNClsTPC(70);
      config->SetDeltaEtaMax(0.019);
      config->SetDeltaPhiMax(0.019);
    } else if (suffix == "21") {
      TrackCutsProtonDCA->SetEtaRange(-0.77, 0.77);
      TrackCutsAntiProtonDCA->SetEtaRange(-0.77, 0.77);
      TrackCutsProtonDCA->SetPID(AliPID::kProton, 0.75, 2.5);
      TrackCutsAntiProtonDCA->SetPID(AliPID::kProton, 0.75, 2.5);
      TrackCutsDeuteronDCA->SetEtaRange(-0.77, 0.77);
      TrackCutsAntiDeuteronDCA->SetEtaRange(-0.77, 0.77);
      TrackCutsDeuteronDCA->SetPID(AliPID::kDeuteron,1.4, 2.7);
      TrackCutsAntiDeuteronDCA->SetPID(AliPID::kDeuteron, 1.4, 2.7);
    } else if (suffix == "22") {
      TrackCutsProtonDCA->SetEtaRange(-0.77, 0.77);
      TrackCutsAntiProtonDCA->SetEtaRange(-0.77, 0.77);
      TrackCutsProtonDCA->SetPID(AliPID::kProton, 0.75, 3.5);
      TrackCutsAntiProtonDCA->SetPID(AliPID::kProton, 0.75, 3.5);
      TrackCutsProtonDCA->SetNClsTPC(90);
      TrackCutsAntiProtonDCA->SetNClsTPC(90);
      TrackCutsDeuteronDCA->SetEtaRange(-0.77, 0.77);
      TrackCutsAntiDeuteronDCA->SetEtaRange(-0.77, 0.77);
      TrackCutsDeuteronDCA->SetPID(AliPID::kDeuteron,1.4, 3.3);
      TrackCutsAntiDeuteronDCA->SetPID(AliPID::kDeuteron, 1.4, 3.3);
      TrackCutsDeuteronDCA->SetNClsTPC(90);
      TrackCutsAntiDeuteronDCA->SetNClsTPC(90);
      config->SetDeltaEtaMax(0.019);
      config->SetDeltaPhiMax(0.019);
    } else if (suffix == "23") {
      TrackCutsProtonDCA->SetPtRange(0.4, 4.05);
      TrackCutsAntiProtonDCA->SetPtRange(0.4, 4.05);
      TrackCutsProtonDCA->SetEtaRange(-0.83, 0.83);
      TrackCutsAntiProtonDCA->SetEtaRange(-0.83, 0.83);
      TrackCutsDeuteronDCA->SetPtRange(0.4, 1.4);
      TrackCutsAntiDeuteronDCA->SetPtRange(0.4, 1.4);
      TrackCutsDeuteronDCA->SetEtaRange(-0.83, 0.83);
      TrackCutsAntiDeuteronDCA->SetEtaRange(-0.83, 0.83);
      config->SetDeltaEtaMax(0.015);
      config->SetDeltaPhiMax(0.015);
    } else if (suffix == "24") {
      TrackCutsProtonDCA->SetPtRange(0.6, 4.05);
      TrackCutsAntiProtonDCA->SetPtRange(0.6, 4.05);
      TrackCutsProtonDCA->SetPID(AliPID::kProton, 0.75, 2.5);
      TrackCutsAntiProtonDCA->SetPID(AliPID::kProton, 0.75, 2.5);
      TrackCutsDeuteronDCA->SetPtRange(0.6, 1.4);
      TrackCutsAntiDeuteronDCA->SetPtRange(0.6, 1.4);
      TrackCutsDeuteronDCA->SetPID(AliPID::kDeuteron,1.4, 2.7);
      TrackCutsAntiDeuteronDCA->SetPID(AliPID::kDeuteron, 1.4, 2.7);
      config->SetDeltaEtaMax(0.019);
      config->SetDeltaPhiMax(0.019);
    } else if (suffix == "25") {
      TrackCutsProtonDCA->SetEtaRange(-0.77, 0.77);
      TrackCutsAntiProtonDCA->SetEtaRange(-0.77, 0.77);
      TrackCutsProtonDCA->SetPID(AliPID::kProton, 0.75, 3.5);
      TrackCutsAntiProtonDCA->SetPID(AliPID::kProton, 0.75, 3.5);
      TrackCutsProtonDCA->SetNClsTPC(90);
      TrackCutsAntiProtonDCA->SetNClsTPC(90);
      TrackCutsDeuteronDCA->SetEtaRange(-0.77, 0.77);
      TrackCutsAntiDeuteronDCA->SetEtaRange(-0.77, 0.77);
      TrackCutsDeuteronDCA->SetPID(AliPID::kDeuteron,1.4, 3.3);
      TrackCutsAntiDeuteronDCA->SetPID(AliPID::kDeuteron, 1.4, 3.3);
      TrackCutsDeuteronDCA->SetNClsTPC(90);
      TrackCutsAntiDeuteronDCA->SetNClsTPC(90);
    } else if (suffix == "26") {
      TrackCutsProtonDCA->SetNClsTPC(90);
      TrackCutsAntiProtonDCA->SetNClsTPC(90);
      TrackCutsDeuteronDCA->SetNClsTPC(90);
      TrackCutsAntiDeuteronDCA->SetNClsTPC(90);
    } else if (suffix == "27") {
      TrackCutsProtonDCA->SetPtRange(0.4, 4.05);
      TrackCutsAntiProtonDCA->SetPtRange(0.4, 4.05);
      TrackCutsProtonDCA->SetEtaRange(-0.77, 0.77);
      TrackCutsAntiProtonDCA->SetEtaRange(-0.77, 0.77);
      TrackCutsProtonDCA->SetPID(AliPID::kProton, 0.75, 2.5);
      TrackCutsAntiProtonDCA->SetPID(AliPID::kProton, 0.75, 2.5);
      TrackCutsDeuteronDCA->SetPtRange(0.4, 1.4);
      TrackCutsAntiDeuteronDCA->SetPtRange(0.4, 1.4);
      TrackCutsDeuteronDCA->SetEtaRange(-0.77, 0.77);
      TrackCutsAntiDeuteronDCA->SetEtaRange(-0.77, 0.77);
      TrackCutsDeuteronDCA->SetPID(AliPID::kDeuteron,1.4, 2.7);
      TrackCutsAntiDeuteronDCA->SetPID(AliPID::kDeuteron, 1.4, 2.7);
      config->SetDeltaEtaMax(0.019);
      config->SetDeltaPhiMax(0.019);
    } else if (suffix == "28") {
      TrackCutsProtonDCA->SetEtaRange(-0.83, 0.83);
      TrackCutsAntiProtonDCA->SetEtaRange(-0.83, 0.83);
      TrackCutsDeuteronDCA->SetEtaRange(-0.83, 0.83);
      TrackCutsAntiDeuteronDCA->SetEtaRange(-0.83, 0.83);
      config->SetDeltaEtaMax(0.019);
      config->SetDeltaPhiMax(0.019);
    } else if (suffix == "29") {
      TrackCutsProtonDCA->SetPtRange(0.4, 4.05);
      TrackCutsAntiProtonDCA->SetPtRange(0.4, 4.05);
      TrackCutsProtonDCA->SetPID(AliPID::kProton, 0.75, 3.5);
      TrackCutsAntiProtonDCA->SetPID(AliPID::kProton, 0.75, 3.5);
      TrackCutsProtonDCA->SetNClsTPC(70);
      TrackCutsAntiProtonDCA->SetNClsTPC(70);
      TrackCutsDeuteronDCA->SetPtRange(0.4, 1.4);
      TrackCutsAntiDeuteronDCA->SetPtRange(0.4, 1.4);
      TrackCutsDeuteronDCA->SetPID(AliPID::kDeuteron,1.4, 3.3);
      TrackCutsAntiDeuteronDCA->SetPID(AliPID::kDeuteron, 1.4, 3.3);
      TrackCutsDeuteronDCA->SetNClsTPC(70);
      TrackCutsAntiDeuteronDCA->SetNClsTPC(70);
      config->SetDeltaEtaMax(0.015);
      config->SetDeltaPhiMax(0.015);
    } else if (suffix == "30") {
      TrackCutsProtonDCA->SetEtaRange(-0.77, 0.77);
      TrackCutsAntiProtonDCA->SetEtaRange(-0.77, 0.77);
      TrackCutsProtonDCA->SetPID(AliPID::kProton, 0.75, 2.5);
      TrackCutsAntiProtonDCA->SetPID(AliPID::kProton, 0.75, 2.5);
      TrackCutsProtonDCA->SetNClsTPC(90);
      TrackCutsAntiProtonDCA->SetNClsTPC(90);
      TrackCutsDeuteronDCA->SetEtaRange(-0.77, 0.77);
      TrackCutsAntiDeuteronDCA->SetEtaRange(-0.77, 0.77);
      TrackCutsDeuteronDCA->SetPID(AliPID::kDeuteron,1.4, 2.3);
      TrackCutsAntiDeuteronDCA->SetPID(AliPID::kDeuteron, 1.4, 2.3);
      TrackCutsDeuteronDCA->SetNClsTPC(90);
      TrackCutsAntiDeuteronDCA->SetNClsTPC(90);
    } else if (suffix == "31") {
      TrackCutsProtonDCA->SetEtaRange(-0.77, 0.77);
      TrackCutsAntiProtonDCA->SetEtaRange(-0.77, 0.77);
      TrackCutsDeuteronDCA->SetEtaRange(-0.77, 0.77);
      TrackCutsAntiDeuteronDCA->SetEtaRange(-0.77, 0.77);
      config->SetDeltaEtaMax(0.015);
      config->SetDeltaPhiMax(0.015);
    } else if (suffix == "32") {
      TrackCutsProtonDCA->SetEtaRange(-0.77, 0.77);
      TrackCutsAntiProtonDCA->SetEtaRange(-0.77, 0.77);
      TrackCutsProtonDCA->SetNClsTPC(70);
      TrackCutsAntiProtonDCA->SetNClsTPC(70);
      TrackCutsDeuteronDCA->SetEtaRange(-0.77, 0.77);
      TrackCutsAntiDeuteronDCA->SetEtaRange(-0.77, 0.77);
      TrackCutsDeuteronDCA->SetNClsTPC(70);
      TrackCutsAntiDeuteronDCA->SetNClsTPC(70);
      config->SetDeltaEtaMax(0.019);
      config->SetDeltaPhiMax(0.019);
    } else if (suffix == "33") {
      TrackCutsProtonDCA->SetPtRange(0.4, 4.05);
      TrackCutsAntiProtonDCA->SetPtRange(0.4, 4.05);
      TrackCutsProtonDCA->SetPID(AliPID::kProton, 0.75, 2.5);
      TrackCutsAntiProtonDCA->SetPID(AliPID::kProton, 0.75, 2.5);
      TrackCutsProtonDCA->SetNClsTPC(70);
      TrackCutsAntiProtonDCA->SetNClsTPC(70);
      TrackCutsDeuteronDCA->SetPtRange(0.4, 1.4);
      TrackCutsAntiDeuteronDCA->SetPtRange(0.4, 1.4);
      TrackCutsDeuteronDCA->SetPID(AliPID::kDeuteron,1.4, 2.3);
      TrackCutsAntiDeuteronDCA->SetPID(AliPID::kDeuteron, 1.4, 2.3);
      TrackCutsDeuteronDCA->SetNClsTPC(70);
      TrackCutsAntiDeuteronDCA->SetNClsTPC(70);
      config->SetDeltaEtaMax(0.019);
      config->SetDeltaPhiMax(0.019);
    } else if (suffix == "34") {
      TrackCutsProtonDCA->SetPtRange(0.6, 4.05);
      TrackCutsAntiProtonDCA->SetPtRange(0.6, 4.05);
      TrackCutsProtonDCA->SetEtaRange(-0.83, 0.83);
      TrackCutsAntiProtonDCA->SetEtaRange(-0.83, 0.83);
      TrackCutsProtonDCA->SetNClsTPC(70);
      TrackCutsAntiProtonDCA->SetNClsTPC(70);
      TrackCutsDeuteronDCA->SetPtRange(0.6, 1.4);
      TrackCutsAntiDeuteronDCA->SetPtRange(0.6, 1.4);
      TrackCutsDeuteronDCA->SetEtaRange(-0.83, 0.83);
      TrackCutsAntiDeuteronDCA->SetEtaRange(-0.83, 0.83);
      TrackCutsDeuteronDCA->SetNClsTPC(70);
      TrackCutsAntiDeuteronDCA->SetNClsTPC(70);
      config->SetDeltaEtaMax(0.019);
      config->SetDeltaPhiMax(0.019);
    } else if (suffix == "35") {
      TrackCutsProtonDCA->SetEtaRange(-0.83, 0.83);
      TrackCutsAntiProtonDCA->SetEtaRange(-0.83, 0.83);
      TrackCutsProtonDCA->SetPID(AliPID::kProton, 0.75, 3.5);
      TrackCutsAntiProtonDCA->SetPID(AliPID::kProton, 0.75, 3.5);
      TrackCutsProtonDCA->SetNClsTPC(90);
      TrackCutsAntiProtonDCA->SetNClsTPC(90);
      TrackCutsDeuteronDCA->SetEtaRange(-0.83, 0.83);
      TrackCutsAntiDeuteronDCA->SetEtaRange(-0.83, 0.83);
      TrackCutsDeuteronDCA->SetPID(AliPID::kDeuteron,1.4, 3.3);
      TrackCutsAntiDeuteronDCA->SetPID(AliPID::kDeuteron, 1.4, 3.3);
      TrackCutsDeuteronDCA->SetNClsTPC(90);
      TrackCutsAntiDeuteronDCA->SetNClsTPC(90);
      config->SetDeltaEtaMax(0.019);
      config->SetDeltaPhiMax(0.019);
    } else if (suffix == "36") {
      TrackCutsProtonDCA->SetEtaRange(-0.77, 0.77);
      TrackCutsAntiProtonDCA->SetEtaRange(-0.77, 0.77);
      TrackCutsProtonDCA->SetNClsTPC(90);
      TrackCutsAntiProtonDCA->SetNClsTPC(90);
      TrackCutsDeuteronDCA->SetEtaRange(-0.77, 0.77);
      TrackCutsAntiDeuteronDCA->SetEtaRange(-0.77, 0.77);
      TrackCutsDeuteronDCA->SetNClsTPC(90);
      TrackCutsAntiDeuteronDCA->SetNClsTPC(90);
      config->SetDeltaEtaMax(0.019);
      config->SetDeltaPhiMax(0.019);
    } else if (suffix == "37") {
      TrackCutsProtonDCA->SetPtRange(0.4, 4.05);
      TrackCutsAntiProtonDCA->SetPtRange(0.4, 4.05);
      TrackCutsProtonDCA->SetPID(AliPID::kProton, 0.75, 3.5);
      TrackCutsAntiProtonDCA->SetPID(AliPID::kProton, 0.75, 3.5);
      TrackCutsProtonDCA->SetNClsTPC(90);
      TrackCutsAntiProtonDCA->SetNClsTPC(90);
      TrackCutsDeuteronDCA->SetPtRange(0.4, 1.4);
      TrackCutsAntiDeuteronDCA->SetPtRange(0.4, 1.4);
      TrackCutsDeuteronDCA->SetPID(AliPID::kDeuteron,1.4, 3.3);
      TrackCutsAntiDeuteronDCA->SetPID(AliPID::kDeuteron, 1.4, 3.3);
      TrackCutsDeuteronDCA->SetNClsTPC(90);
      TrackCutsAntiDeuteronDCA->SetNClsTPC(90);
      config->SetDeltaEtaMax(0.019);
      config->SetDeltaPhiMax(0.019);
    } else if (suffix == "38") {
      TrackCutsProtonDCA->SetPtRange(0.6, 4.05);
      TrackCutsAntiProtonDCA->SetPtRange(0.6, 4.05);
      TrackCutsProtonDCA->SetNClsTPC(70);
      TrackCutsAntiProtonDCA->SetNClsTPC(70);
      TrackCutsDeuteronDCA->SetPtRange(0.6, 1.4);
      TrackCutsAntiDeuteronDCA->SetPtRange(0.6, 1.4);
      TrackCutsDeuteronDCA->SetNClsTPC(70);
      TrackCutsAntiDeuteronDCA->SetNClsTPC(70);
    } else if (suffix == "39") {
      TrackCutsProtonDCA->SetPtRange(0.4, 4.05);
      TrackCutsAntiProtonDCA->SetPtRange(0.4, 4.05);
      TrackCutsProtonDCA->SetEtaRange(-0.77, 0.77);
      TrackCutsAntiProtonDCA->SetEtaRange(-0.77, 0.77);
      TrackCutsProtonDCA->SetPID(AliPID::kProton, 0.75, 3.5);
      TrackCutsAntiProtonDCA->SetPID(AliPID::kProton, 0.75, 3.5);
      TrackCutsDeuteronDCA->SetPtRange(0.4, 1.4);
      TrackCutsAntiDeuteronDCA->SetPtRange(0.4, 1.4);
      TrackCutsDeuteronDCA->SetEtaRange(-0.77, 0.77);
      TrackCutsAntiDeuteronDCA->SetEtaRange(-0.77, 0.77);
      TrackCutsDeuteronDCA->SetPID(AliPID::kDeuteron,1.4, 3.3);
      TrackCutsAntiDeuteronDCA->SetPID(AliPID::kDeuteron, 1.4, 3.3);
      config->SetDeltaEtaMax(0.015);
      config->SetDeltaPhiMax(0.015);
    } else if (suffix == "40") {
      TrackCutsProtonDCA->SetPtRange(0.6, 4.05);
      TrackCutsAntiProtonDCA->SetPtRange(0.6, 4.05);
      TrackCutsProtonDCA->SetPID(AliPID::kProton, 0.75, 2.5);
      TrackCutsAntiProtonDCA->SetPID(AliPID::kProton, 0.75, 2.5);
      TrackCutsDeuteronDCA->SetPtRange(0.6, 1.4);
      TrackCutsAntiDeuteronDCA->SetPtRange(0.6, 1.4);
      TrackCutsDeuteronDCA->SetPID(AliPID::kDeuteron,1.4, 2.3);
      TrackCutsAntiDeuteronDCA->SetPID(AliPID::kDeuteron, 1.4, 2.3);
    } else if (suffix == "41") {
      TrackCutsProtonDCA->SetPtRange(0.6, 4.05);
      TrackCutsAntiProtonDCA->SetPtRange(0.6, 4.05);
      TrackCutsProtonDCA->SetEtaRange(-0.83, 0.83);
      TrackCutsAntiProtonDCA->SetEtaRange(-0.83, 0.83);
      TrackCutsProtonDCA->SetPID(AliPID::kProton, 0.75, 3.5);
      TrackCutsAntiProtonDCA->SetPID(AliPID::kProton, 0.75, 3.5);
      TrackCutsDeuteronDCA->SetPtRange(0.6, 1.4);
      TrackCutsAntiDeuteronDCA->SetPtRange(0.6, 1.4);
      TrackCutsDeuteronDCA->SetEtaRange(-0.83, 0.83);
      TrackCutsAntiDeuteronDCA->SetEtaRange(-0.83, 0.83);
      TrackCutsDeuteronDCA->SetPID(AliPID::kDeuteron,1.4, 3.3);
      TrackCutsAntiDeuteronDCA->SetPID(AliPID::kDeuteron, 1.4, 3.3);
      config->SetDeltaEtaMax(0.019);
      config->SetDeltaPhiMax(0.019);
    } else if (suffix == "42") {
      TrackCutsProtonDCA->SetEtaRange(-0.83, 0.83);
      TrackCutsAntiProtonDCA->SetEtaRange(-0.83, 0.83);
      TrackCutsProtonDCA->SetPID(AliPID::kProton, 0.75, 3.5);
      TrackCutsAntiProtonDCA->SetPID(AliPID::kProton, 0.75, 3.5);
      TrackCutsDeuteronDCA->SetEtaRange(-0.83, 0.83);
      TrackCutsAntiDeuteronDCA->SetEtaRange(-0.83, 0.83);
      TrackCutsDeuteronDCA->SetPID(AliPID::kDeuteron,1.4, 3.3);
      TrackCutsAntiDeuteronDCA->SetPID(AliPID::kDeuteron, 1.4, 3.3);
      config->SetDeltaEtaMax(0.019);
      config->SetDeltaPhiMax(0.019);
    }else if (suffix == "43") {
      TrackCutsProtonDCA->SetDCAVtxZ(0.16);
      TrackCutsProtonDCA->SetDCAVtxXY(0.08);
      TrackCutsProtonDCA->SetDCAVtxZ(0.16);
      TrackCutsProtonDCA->SetDCAVtxXY(0.08);
      TrackCutsDeuteronDCA->SetDCAVtxZ(0.16);
      TrackCutsDeuteronDCA->SetDCAVtxXY(0.08);
      TrackCutsDeuteronDCA->SetDCAVtxZ(0.16);
      TrackCutsAntiDeuteronDCA->SetDCAVtxXY(0.08);
    }else if (suffix == "44") {
      TrackCutsProtonDCA->SetDCAVtxZ(0.24);
      TrackCutsProtonDCA->SetDCAVtxXY(0.12);
      TrackCutsProtonDCA->SetDCAVtxZ(0.24);
      TrackCutsProtonDCA->SetDCAVtxXY(0.12);
      TrackCutsDeuteronDCA->SetDCAVtxZ(0.24);
      TrackCutsDeuteronDCA->SetDCAVtxXY(0.12);
      TrackCutsDeuteronDCA->SetDCAVtxZ(0.24);
      TrackCutsAntiDeuteronDCA->SetDCAVtxXY(0.12);
    }
 }

 if (SystematicFullpT) {

    if (suffix == "0"){
      TrackCutsDeuteronDCA->SetPtRange(0.5, 4.05);
      TrackCutsAntiDeuteronDCA->SetPtRange(0.5, 4.05);
    }else if (suffix == "1") {
       TrackCutsProtonDCA->SetPID(AliPID::kProton, 0.75, 2.5);
       TrackCutsAntiProtonDCA->SetPID(AliPID::kProton, 0.75, 2.5);
       TrackCutsProtonDCA->SetNClsTPC(90);
       TrackCutsAntiProtonDCA->SetNClsTPC(90);
       TrackCutsDeuteronDCA->SetPID(AliPID::kDeuteron,1.4, 2.7);
       TrackCutsAntiDeuteronDCA->SetPID(AliPID::kDeuteron, 1.4, 2.7);
       TrackCutsDeuteronDCA->SetNClsTPC(90);
       TrackCutsAntiDeuteronDCA->SetNClsTPC(90);
     } else if (suffix == "2") {
       TrackCutsProtonDCA->SetPtRange(0.6, 4.05);
       TrackCutsAntiProtonDCA->SetPtRange(0.6, 4.05);
       TrackCutsProtonDCA->SetEtaRange(-0.77, 0.77);
       TrackCutsAntiProtonDCA->SetEtaRange(-0.77, 0.77);
       TrackCutsDeuteronDCA->SetEtaRange(-0.77, 0.77);
       TrackCutsAntiDeuteronDCA->SetEtaRange(-0.77, 0.77);
       config->SetDeltaEtaMax(0.019);
       config->SetDeltaPhiMax(0.019);
     }  else if (suffix == "3") {
       TrackCutsProtonDCA->SetEtaRange(-0.77, 0.77);
       TrackCutsAntiProtonDCA->SetEtaRange(-0.77, 0.77);
       TrackCutsProtonDCA->SetNClsTPC(70);
       TrackCutsAntiProtonDCA->SetNClsTPC(70);
       TrackCutsDeuteronDCA->SetEtaRange(-0.77, 0.77);
       TrackCutsAntiDeuteronDCA->SetEtaRange(-0.77, 0.77);
       TrackCutsDeuteronDCA->SetNClsTPC(70);
       TrackCutsAntiDeuteronDCA->SetNClsTPC(70);
     } else if (suffix == "4") {
       TrackCutsProtonDCA->SetPtRange(0.4, 4.05);
       TrackCutsAntiProtonDCA->SetPtRange(0.4, 4.05);
       TrackCutsProtonDCA->SetEtaRange(-0.77, 0.77);
       TrackCutsAntiProtonDCA->SetEtaRange(-0.77, 0.77);
       TrackCutsProtonDCA->SetPID(AliPID::kProton, 0.75, 3.5);
       TrackCutsAntiProtonDCA->SetPID(AliPID::kProton, 0.75, 3.5);
       TrackCutsDeuteronDCA->SetPtRange(0.4, 4.05);
       TrackCutsAntiProtonDCA->SetPtRange(0.4, 4.05);
       TrackCutsDeuteronDCA->SetEtaRange(-0.77, 0.77);
       TrackCutsAntiProtonDCA->SetEtaRange(-0.77, 0.77);
       TrackCutsDeuteronDCA->SetPID(AliPID::kDeuteron,1.4, 3.3);
       TrackCutsAntiDeuteronDCA->SetPID(AliPID::kDeuteron, 1.4, 3.3);
       config->SetDeltaEtaMax(0.019);
       config->SetDeltaPhiMax(0.019);
     } else if (suffix == "5") {
       TrackCutsProtonDCA->SetPID(AliPID::kProton, 0.75, 2.5);
       TrackCutsAntiProtonDCA->SetPID(AliPID::kProton, 0.75, 2.5);
       TrackCutsProtonDCA->SetNClsTPC(90);
       TrackCutsAntiProtonDCA->SetNClsTPC(90);
       TrackCutsDeuteronDCA->SetPID(AliPID::kDeuteron,1.4, 2.7);
       TrackCutsAntiDeuteronDCA->SetPID(AliPID::kDeuteron, 1.4, 2.7);
       TrackCutsDeuteronDCA->SetNClsTPC(90);
       TrackCutsAntiDeuteronDCA->SetNClsTPC(90);
     } else if (suffix == "6") {
       TrackCutsProtonDCA->SetEtaRange(-0.83, 0.83);
       TrackCutsAntiProtonDCA->SetEtaRange(-0.83, 0.83);
       TrackCutsProtonDCA->SetPID(AliPID::kProton, 0.75, 3.5);
       TrackCutsAntiProtonDCA->SetPID(AliPID::kProton, 0.75, 3.5);
       TrackCutsDeuteronDCA->SetEtaRange(-0.83, 0.83);
       TrackCutsAntiDeuteronDCA->SetEtaRange(-0.83, 0.83);
       TrackCutsDeuteronDCA->SetPID(AliPID::kDeuteron,1.4, 3.3);
       TrackCutsAntiDeuteronDCA->SetPID(AliPID::kDeuteron, 1.4, 3.3);
       config->SetDeltaEtaMax(0.015);
       config->SetDeltaPhiMax(0.015);
     } else if (suffix == "7") {
       TrackCutsProtonDCA->SetPtRange(0.4, 4.05);
       TrackCutsAntiProtonDCA->SetPtRange(0.4, 4.05);
       TrackCutsProtonDCA->SetEtaRange(-0.77, 0.77);
       TrackCutsAntiProtonDCA->SetEtaRange(-0.77, 0.77);
       TrackCutsProtonDCA->SetNClsTPC(90);
       TrackCutsAntiProtonDCA->SetNClsTPC(90);
       TrackCutsDeuteronDCA->SetPtRange(0.4, 4.05);
       TrackCutsAntiDeuteronDCA->SetPtRange(0.4, 4.05);
       TrackCutsDeuteronDCA->SetEtaRange(-0.77, 0.77);
       TrackCutsAntiDeuteronDCA->SetEtaRange(-0.77, 0.77);
       TrackCutsDeuteronDCA->SetNClsTPC(90);
       TrackCutsAntiDeuteronDCA->SetNClsTPC(90);
       config->SetDeltaEtaMax(0.019);
       config->SetDeltaPhiMax(0.019);
     } else if (suffix == "8") {
       TrackCutsProtonDCA->SetPtRange(0.4, 4.05);
       TrackCutsAntiProtonDCA->SetPtRange(0.4, 4.05);
       TrackCutsProtonDCA->SetEtaRange(-0.83, 0.83);
       TrackCutsAntiProtonDCA->SetEtaRange(-0.83, 0.83);
       TrackCutsDeuteronDCA->SetPtRange(0.4, 4.05);
       TrackCutsAntiDeuteronDCA->SetPtRange(0.4, 4.05);
       TrackCutsDeuteronDCA->SetEtaRange(-0.83, 0.83);
       TrackCutsAntiDeuteronDCA->SetEtaRange(-0.83, 0.83);
     } else if (suffix == "9") {
       TrackCutsProtonDCA->SetEtaRange(-0.77, 0.77);
       TrackCutsAntiProtonDCA->SetEtaRange(-0.77, 0.77);
       TrackCutsProtonDCA->SetPID(AliPID::kProton, 0.75, 2.5);
       TrackCutsAntiProtonDCA->SetPID(AliPID::kProton, 0.75, 2.5);
       TrackCutsDeuteronDCA->SetEtaRange(-0.77, 0.77);
       TrackCutsAntiDeuteronDCA->SetEtaRange(-0.77, 0.77);
       TrackCutsDeuteronDCA->SetPID(AliPID::kDeuteron,1.4, 2.7);
       TrackCutsAntiDeuteronDCA->SetPID(AliPID::kDeuteron, 1.4, 2.7);
     } else if (suffix == "10") {
       TrackCutsProtonDCA->SetPtRange(0.4, 4.05);
       TrackCutsAntiProtonDCA->SetPtRange(0.4, 4.05);
       TrackCutsProtonDCA->SetEtaRange(-0.83, 0.83);
       TrackCutsAntiProtonDCA->SetEtaRange(-0.83, 0.83);
       TrackCutsProtonDCA->SetPID(AliPID::kProton, 0.75, 2.5);
       TrackCutsAntiProtonDCA->SetPID(AliPID::kProton, 0.75, 2.5);
       TrackCutsDeuteronDCA->SetPtRange(0.4, 4.05);
       TrackCutsAntiDeuteronDCA->SetPtRange(0.4, 4.05);
       TrackCutsDeuteronDCA->SetEtaRange(-0.83, 0.83);
       TrackCutsAntiDeuteronDCA->SetEtaRange(-0.83, 0.83);
       TrackCutsDeuteronDCA->SetPID(AliPID::kDeuteron,1.4, 2.7);
       TrackCutsAntiDeuteronDCA->SetPID(AliPID::kDeuteron, 1.4, 2.7);
       config->SetDeltaEtaMax(0.019);
       config->SetDeltaPhiMax(0.019);
     } else if (suffix == "11") {
       TrackCutsProtonDCA->SetEtaRange(-0.83, 0.83);
       TrackCutsAntiProtonDCA->SetEtaRange(-0.83, 0.83);
       TrackCutsProtonDCA->SetPID(AliPID::kProton, 0.75, 2.5);
       TrackCutsAntiProtonDCA->SetPID(AliPID::kProton, 0.75, 2.5);
       TrackCutsProtonDCA->SetNClsTPC(90);
       TrackCutsAntiProtonDCA->SetNClsTPC(90);
       TrackCutsDeuteronDCA->SetEtaRange(-0.83, 0.83);
       TrackCutsAntiDeuteronDCA->SetEtaRange(-0.83, 0.83);
       TrackCutsDeuteronDCA->SetPID(AliPID::kDeuteron,1.4, 2.7);
       TrackCutsAntiDeuteronDCA->SetPID(AliPID::kDeuteron, 1.4, 2.7);
       TrackCutsDeuteronDCA->SetNClsTPC(90);
       TrackCutsAntiDeuteronDCA->SetNClsTPC(90);
     } else if (suffix == "12") {
       TrackCutsProtonDCA->SetPtRange(0.4, 4.05);
       TrackCutsAntiProtonDCA->SetPtRange(0.4, 4.05);
       TrackCutsProtonDCA->SetPID(AliPID::kProton, 0.75, 2.5);
       TrackCutsAntiProtonDCA->SetPID(AliPID::kProton, 0.75, 2.5);
       TrackCutsDeuteronDCA->SetPtRange(0.4, 4.05);
       TrackCutsAntiDeuteronDCA->SetPtRange(0.4, 4.05);
       TrackCutsDeuteronDCA->SetPID(AliPID::kDeuteron,1.4, 2.7);
       TrackCutsAntiDeuteronDCA->SetPID(AliPID::kDeuteron, 1.4, 2.7);
       config->SetDeltaEtaMax(0.019);
       config->SetDeltaPhiMax(0.019);
     } else if (suffix == "13") {
       TrackCutsProtonDCA->SetPtRange(0.6, 4.05);
       TrackCutsAntiProtonDCA->SetPtRange(0.6, 4.05);
       TrackCutsProtonDCA->SetEtaRange(-0.83, 0.83);
       TrackCutsAntiProtonDCA->SetEtaRange(-0.83, 0.83);
       TrackCutsProtonDCA->SetPID(AliPID::kProton, 0.75, 2.5);
       TrackCutsAntiProtonDCA->SetPID(AliPID::kProton, 0.75, 2.5);
       TrackCutsProtonDCA->SetNClsTPC(90);
       TrackCutsAntiProtonDCA->SetNClsTPC(90);
       TrackCutsDeuteronDCA->SetPtRange(0.6, 4.05);
       TrackCutsAntiDeuteronDCA->SetPtRange(0.6, 4.05);
       TrackCutsDeuteronDCA->SetEtaRange(-0.83, 0.83);
       TrackCutsAntiDeuteronDCA->SetEtaRange(-0.83, 0.83);
       TrackCutsDeuteronDCA->SetPID(AliPID::kDeuteron,1.4, 2.7);
       TrackCutsAntiDeuteronDCA->SetPID(AliPID::kDeuteron, 1.4, 2.7);
       TrackCutsDeuteronDCA->SetNClsTPC(90);
       TrackCutsAntiDeuteronDCA->SetNClsTPC(90);
     } else if (suffix == "14") {
       TrackCutsProtonDCA->SetPtRange(0.4, 4.05);
       TrackCutsAntiProtonDCA->SetPtRange(0.4, 4.05);
       TrackCutsProtonDCA->SetPID(AliPID::kProton, 0.75, 3.5);
       TrackCutsAntiProtonDCA->SetPID(AliPID::kProton, 0.75, 3.5);
       TrackCutsProtonDCA->SetNClsTPC(70);
       TrackCutsAntiProtonDCA->SetNClsTPC(70);
       TrackCutsDeuteronDCA->SetPtRange(0.4, 4.05);
       TrackCutsAntiDeuteronDCA->SetPtRange(0.4, 4.05);
       TrackCutsDeuteronDCA->SetPID(AliPID::kDeuteron,1.4, 3.3);
       TrackCutsAntiDeuteronDCA->SetPID(AliPID::kDeuteron, 1.4, 3.3);
       TrackCutsDeuteronDCA->SetNClsTPC(70);
       TrackCutsAntiDeuteronDCA->SetNClsTPC(70);
       config->SetDeltaEtaMax(0.019);
       config->SetDeltaPhiMax(0.019);
     } else if (suffix == "15") {
       TrackCutsProtonDCA->SetPtRange(0.4, 4.05);
       TrackCutsAntiProtonDCA->SetPtRange(0.4, 4.05);
       TrackCutsProtonDCA->SetEtaRange(-0.83, 0.83);
       TrackCutsAntiProtonDCA->SetEtaRange(-0.83, 0.83);
       TrackCutsProtonDCA->SetPID(AliPID::kProton, 0.75, 2.5);
       TrackCutsAntiProtonDCA->SetPID(AliPID::kProton, 0.75, 2.5);
       TrackCutsDeuteronDCA->SetPtRange(0.4, 4.05);
       TrackCutsAntiDeuteronDCA->SetPtRange(0.4, 4.05);
       TrackCutsDeuteronDCA->SetEtaRange(-0.83, 0.83);
       TrackCutsAntiDeuteronDCA->SetEtaRange(-0.83, 0.83);
       TrackCutsDeuteronDCA->SetPID(AliPID::kDeuteron,1.4, 2.7);
       TrackCutsAntiDeuteronDCA->SetPID(AliPID::kDeuteron, 1.4, 2.7);
       config->SetDeltaEtaMax(0.015);
       config->SetDeltaPhiMax(0.015);
     } else if (suffix == "16") {
       TrackCutsProtonDCA->SetPtRange(0.4, 4.05);
       TrackCutsAntiProtonDCA->SetPtRange(0.4, 4.05);
       TrackCutsProtonDCA->SetEtaRange(-0.83, 0.83);
       TrackCutsAntiProtonDCA->SetEtaRange(-0.83, 0.83);
       TrackCutsProtonDCA->SetPID(AliPID::kProton, 0.75, 2.5);
       TrackCutsAntiProtonDCA->SetPID(AliPID::kProton, 0.75, 2.5);
       TrackCutsProtonDCA->SetNClsTPC(90);
       TrackCutsAntiProtonDCA->SetNClsTPC(90);
       TrackCutsDeuteronDCA->SetPtRange(0.4, 4.05);
       TrackCutsAntiDeuteronDCA->SetPtRange(0.4, 4.05);
       TrackCutsDeuteronDCA->SetEtaRange(-0.83, 0.83);
       TrackCutsAntiDeuteronDCA->SetEtaRange(-0.83, 0.83);
       TrackCutsDeuteronDCA->SetPID(AliPID::kDeuteron,1.4, 2.7);
       TrackCutsAntiDeuteronDCA->SetPID(AliPID::kDeuteron, 1.4, 2.7);
       TrackCutsDeuteronDCA->SetNClsTPC(90);
       TrackCutsAntiDeuteronDCA->SetNClsTPC(90);
     } else if (suffix == "17") {
       TrackCutsProtonDCA->SetPID(AliPID::kProton, 0.75, 2.5);
       TrackCutsAntiProtonDCA->SetPID(AliPID::kProton, 0.75, 2.5);
       TrackCutsProtonDCA->SetNClsTPC(70);
       TrackCutsAntiProtonDCA->SetNClsTPC(70);
       TrackCutsDeuteronDCA->SetPID(AliPID::kDeuteron,1.4, 2.7);
       TrackCutsAntiDeuteronDCA->SetPID(AliPID::kDeuteron, 1.4, 2.7);
       TrackCutsDeuteronDCA->SetNClsTPC(70);
       TrackCutsAntiDeuteronDCA->SetNClsTPC(70);
       config->SetDeltaEtaMax(0.019);
       config->SetDeltaPhiMax(0.019);
     } else if (suffix == "18") {
       TrackCutsProtonDCA->SetPtRange(0.4, 4.05);
       TrackCutsAntiProtonDCA->SetPtRange(0.4, 4.05);
       TrackCutsProtonDCA->SetNClsTPC(90);
       TrackCutsAntiProtonDCA->SetNClsTPC(90);
       TrackCutsDeuteronDCA->SetPtRange(0.4, 4.05);
       TrackCutsAntiDeuteronDCA->SetPtRange(0.4, 4.05);
       TrackCutsDeuteronDCA->SetNClsTPC(90);
       TrackCutsAntiDeuteronDCA->SetNClsTPC(90);
       config->SetDeltaEtaMax(0.019);
       config->SetDeltaPhiMax(0.019);
     } else if (suffix == "19") {
       TrackCutsProtonDCA->SetEtaRange(-0.83, 0.83);
       TrackCutsAntiProtonDCA->SetEtaRange(-0.83, 0.83);
       TrackCutsProtonDCA->SetPID(AliPID::kProton, 0.75, 3.5);
       TrackCutsAntiProtonDCA->SetPID(AliPID::kProton, 0.75, 3.5);
       TrackCutsProtonDCA->SetNClsTPC(70);
       TrackCutsAntiProtonDCA->SetNClsTPC(70);
       TrackCutsDeuteronDCA->SetEtaRange(-0.83, 0.83);
       TrackCutsAntiDeuteronDCA->SetEtaRange(-0.83, 0.83);
       TrackCutsDeuteronDCA->SetPID(AliPID::kDeuteron,1.4, 3.3);
       TrackCutsAntiDeuteronDCA->SetPID(AliPID::kDeuteron, 1.4, 3.3);
       TrackCutsDeuteronDCA->SetNClsTPC(70);
       TrackCutsAntiDeuteronDCA->SetNClsTPC(70);
       config->SetDeltaEtaMax(0.019);
       config->SetDeltaPhiMax(0.019);
     } else if (suffix == "20") {
       TrackCutsProtonDCA->SetPtRange(0.4, 4.05);
       TrackCutsAntiProtonDCA->SetPtRange(0.4, 4.05);
       TrackCutsProtonDCA->SetEtaRange(-0.77, 0.77);
       TrackCutsAntiProtonDCA->SetEtaRange(-0.77, 0.77);
       TrackCutsProtonDCA->SetNClsTPC(70);
       TrackCutsAntiProtonDCA->SetNClsTPC(70);
       TrackCutsDeuteronDCA->SetPtRange(0.4, 4.05);
       TrackCutsAntiDeuteronDCA->SetPtRange(0.4, 4.05);
       TrackCutsDeuteronDCA->SetEtaRange(-0.77, 0.77);
       TrackCutsAntiDeuteronDCA->SetEtaRange(-0.77, 0.77);
       TrackCutsDeuteronDCA->SetNClsTPC(70);
       TrackCutsAntiDeuteronDCA->SetNClsTPC(70);
       config->SetDeltaEtaMax(0.019);
       config->SetDeltaPhiMax(0.019);
     } else if (suffix == "21") {
       TrackCutsProtonDCA->SetEtaRange(-0.77, 0.77);
       TrackCutsAntiProtonDCA->SetEtaRange(-0.77, 0.77);
       TrackCutsProtonDCA->SetPID(AliPID::kProton, 0.75, 2.5);
       TrackCutsAntiProtonDCA->SetPID(AliPID::kProton, 0.75, 2.5);
       TrackCutsDeuteronDCA->SetEtaRange(-0.77, 0.77);
       TrackCutsAntiDeuteronDCA->SetEtaRange(-0.77, 0.77);
       TrackCutsDeuteronDCA->SetPID(AliPID::kDeuteron,1.4, 2.7);
       TrackCutsAntiDeuteronDCA->SetPID(AliPID::kDeuteron, 1.4, 2.7);
     } else if (suffix == "22") {
       TrackCutsProtonDCA->SetEtaRange(-0.77, 0.77);
       TrackCutsAntiProtonDCA->SetEtaRange(-0.77, 0.77);
       TrackCutsProtonDCA->SetPID(AliPID::kProton, 0.75, 3.5);
       TrackCutsAntiProtonDCA->SetPID(AliPID::kProton, 0.75, 3.5);
       TrackCutsProtonDCA->SetNClsTPC(90);
       TrackCutsAntiProtonDCA->SetNClsTPC(90);
       TrackCutsDeuteronDCA->SetEtaRange(-0.77, 0.77);
       TrackCutsAntiDeuteronDCA->SetEtaRange(-0.77, 0.77);
       TrackCutsDeuteronDCA->SetPID(AliPID::kDeuteron,1.4, 3.3);
       TrackCutsAntiDeuteronDCA->SetPID(AliPID::kDeuteron, 1.4, 3.3);
       TrackCutsDeuteronDCA->SetNClsTPC(90);
       TrackCutsAntiDeuteronDCA->SetNClsTPC(90);
       config->SetDeltaEtaMax(0.019);
       config->SetDeltaPhiMax(0.019);
     } else if (suffix == "23") {
       TrackCutsProtonDCA->SetPtRange(0.4, 4.05);
       TrackCutsAntiProtonDCA->SetPtRange(0.4, 4.05);
       TrackCutsProtonDCA->SetEtaRange(-0.83, 0.83);
       TrackCutsAntiProtonDCA->SetEtaRange(-0.83, 0.83);
       TrackCutsDeuteronDCA->SetPtRange(0.4, 4.05);
       TrackCutsAntiDeuteronDCA->SetPtRange(0.4, 4.05);
       TrackCutsDeuteronDCA->SetEtaRange(-0.83, 0.83);
       TrackCutsAntiDeuteronDCA->SetEtaRange(-0.83, 0.83);
       config->SetDeltaEtaMax(0.015);
       config->SetDeltaPhiMax(0.015);
     } else if (suffix == "24") {
       TrackCutsProtonDCA->SetPtRange(0.6, 4.05);
       TrackCutsAntiProtonDCA->SetPtRange(0.6, 4.05);
       TrackCutsProtonDCA->SetPID(AliPID::kProton, 0.75, 2.5);
       TrackCutsAntiProtonDCA->SetPID(AliPID::kProton, 0.75, 2.5);
       TrackCutsDeuteronDCA->SetPtRange(0.6, 4.05);
       TrackCutsAntiDeuteronDCA->SetPtRange(0.6, 4.05);
       TrackCutsDeuteronDCA->SetPID(AliPID::kDeuteron,1.4, 2.7);
       TrackCutsAntiDeuteronDCA->SetPID(AliPID::kDeuteron, 1.4, 2.7);
       config->SetDeltaEtaMax(0.019);
       config->SetDeltaPhiMax(0.019);
     } else if (suffix == "25") {
       TrackCutsProtonDCA->SetEtaRange(-0.77, 0.77);
       TrackCutsAntiProtonDCA->SetEtaRange(-0.77, 0.77);
       TrackCutsProtonDCA->SetPID(AliPID::kProton, 0.75, 3.5);
       TrackCutsAntiProtonDCA->SetPID(AliPID::kProton, 0.75, 3.5);
       TrackCutsProtonDCA->SetNClsTPC(90);
       TrackCutsAntiProtonDCA->SetNClsTPC(90);
       TrackCutsDeuteronDCA->SetEtaRange(-0.77, 0.77);
       TrackCutsAntiDeuteronDCA->SetEtaRange(-0.77, 0.77);
       TrackCutsDeuteronDCA->SetPID(AliPID::kDeuteron,1.4, 3.3);
       TrackCutsAntiDeuteronDCA->SetPID(AliPID::kDeuteron, 1.4, 3.3);
       TrackCutsDeuteronDCA->SetNClsTPC(90);
       TrackCutsAntiDeuteronDCA->SetNClsTPC(90);
     } else if (suffix == "26") {
       TrackCutsProtonDCA->SetNClsTPC(90);
       TrackCutsAntiProtonDCA->SetNClsTPC(90);
       TrackCutsDeuteronDCA->SetNClsTPC(90);
       TrackCutsAntiDeuteronDCA->SetNClsTPC(90);
     } else if (suffix == "27") {
       TrackCutsProtonDCA->SetPtRange(0.4, 4.05);
       TrackCutsAntiProtonDCA->SetPtRange(0.4, 4.05);
       TrackCutsProtonDCA->SetEtaRange(-0.77, 0.77);
       TrackCutsAntiProtonDCA->SetEtaRange(-0.77, 0.77);
       TrackCutsProtonDCA->SetPID(AliPID::kProton, 0.75, 2.5);
       TrackCutsAntiProtonDCA->SetPID(AliPID::kProton, 0.75, 2.5);
       TrackCutsDeuteronDCA->SetPtRange(0.4, 4.05);
       TrackCutsAntiDeuteronDCA->SetPtRange(0.4, 4.05);
       TrackCutsDeuteronDCA->SetEtaRange(-0.77, 0.77);
       TrackCutsAntiDeuteronDCA->SetEtaRange(-0.77, 0.77);
       TrackCutsDeuteronDCA->SetPID(AliPID::kDeuteron,1.4, 2.7);
       TrackCutsAntiDeuteronDCA->SetPID(AliPID::kDeuteron, 1.4, 2.7);
       config->SetDeltaEtaMax(0.019);
       config->SetDeltaPhiMax(0.019);
     } else if (suffix == "28") {
       TrackCutsProtonDCA->SetEtaRange(-0.83, 0.83);
       TrackCutsAntiProtonDCA->SetEtaRange(-0.83, 0.83);
       TrackCutsDeuteronDCA->SetEtaRange(-0.83, 0.83);
       TrackCutsAntiDeuteronDCA->SetEtaRange(-0.83, 0.83);
       config->SetDeltaEtaMax(0.019);
       config->SetDeltaPhiMax(0.019);
     } else if (suffix == "29") {
       TrackCutsProtonDCA->SetPtRange(0.4, 4.05);
       TrackCutsAntiProtonDCA->SetPtRange(0.4, 4.05);
       TrackCutsProtonDCA->SetPID(AliPID::kProton, 0.75, 3.5);
       TrackCutsAntiProtonDCA->SetPID(AliPID::kProton, 0.75, 3.5);
       TrackCutsProtonDCA->SetNClsTPC(70);
       TrackCutsAntiProtonDCA->SetNClsTPC(70);
       TrackCutsDeuteronDCA->SetPtRange(0.4, 4.05);
       TrackCutsAntiDeuteronDCA->SetPtRange(0.4, 4.05);
       TrackCutsDeuteronDCA->SetPID(AliPID::kDeuteron,1.4, 3.3);
       TrackCutsAntiDeuteronDCA->SetPID(AliPID::kDeuteron, 1.4, 3.3);
       TrackCutsDeuteronDCA->SetNClsTPC(70);
       TrackCutsAntiDeuteronDCA->SetNClsTPC(70);
       config->SetDeltaEtaMax(0.015);
       config->SetDeltaPhiMax(0.015);
     } else if (suffix == "30") {
       TrackCutsProtonDCA->SetEtaRange(-0.77, 0.77);
       TrackCutsAntiProtonDCA->SetEtaRange(-0.77, 0.77);
       TrackCutsProtonDCA->SetPID(AliPID::kProton, 0.75, 2.5);
       TrackCutsAntiProtonDCA->SetPID(AliPID::kProton, 0.75, 2.5);
       TrackCutsProtonDCA->SetNClsTPC(90);
       TrackCutsAntiProtonDCA->SetNClsTPC(90);
       TrackCutsDeuteronDCA->SetEtaRange(-0.77, 0.77);
       TrackCutsAntiDeuteronDCA->SetEtaRange(-0.77, 0.77);
       TrackCutsDeuteronDCA->SetPID(AliPID::kDeuteron,1.4, 2.3);
       TrackCutsAntiDeuteronDCA->SetPID(AliPID::kDeuteron, 1.4, 2.3);
       TrackCutsDeuteronDCA->SetNClsTPC(90);
       TrackCutsAntiDeuteronDCA->SetNClsTPC(90);
     } else if (suffix == "31") {
       TrackCutsProtonDCA->SetEtaRange(-0.77, 0.77);
       TrackCutsAntiProtonDCA->SetEtaRange(-0.77, 0.77);
       TrackCutsDeuteronDCA->SetEtaRange(-0.77, 0.77);
       TrackCutsAntiDeuteronDCA->SetEtaRange(-0.77, 0.77);
       config->SetDeltaEtaMax(0.015);
       config->SetDeltaPhiMax(0.015);
     } else if (suffix == "32") {
       TrackCutsProtonDCA->SetEtaRange(-0.77, 0.77);
       TrackCutsAntiProtonDCA->SetEtaRange(-0.77, 0.77);
       TrackCutsProtonDCA->SetNClsTPC(70);
       TrackCutsAntiProtonDCA->SetNClsTPC(70);
       TrackCutsDeuteronDCA->SetEtaRange(-0.77, 0.77);
       TrackCutsAntiDeuteronDCA->SetEtaRange(-0.77, 0.77);
       TrackCutsDeuteronDCA->SetNClsTPC(70);
       TrackCutsAntiDeuteronDCA->SetNClsTPC(70);
       config->SetDeltaEtaMax(0.019);
       config->SetDeltaPhiMax(0.019);
     } else if (suffix == "33") {
       TrackCutsProtonDCA->SetPtRange(0.4, 4.05);
       TrackCutsAntiProtonDCA->SetPtRange(0.4, 4.05);
       TrackCutsProtonDCA->SetPID(AliPID::kProton, 0.75, 2.5);
       TrackCutsAntiProtonDCA->SetPID(AliPID::kProton, 0.75, 2.5);
       TrackCutsProtonDCA->SetNClsTPC(70);
       TrackCutsAntiProtonDCA->SetNClsTPC(70);
       TrackCutsDeuteronDCA->SetPtRange(0.4, 4.05);
       TrackCutsAntiDeuteronDCA->SetPtRange(0.4, 4.05);
       TrackCutsDeuteronDCA->SetPID(AliPID::kDeuteron,1.4, 2.3);
       TrackCutsAntiDeuteronDCA->SetPID(AliPID::kDeuteron, 1.4, 2.3);
       TrackCutsDeuteronDCA->SetNClsTPC(70);
       TrackCutsAntiDeuteronDCA->SetNClsTPC(70);
       config->SetDeltaEtaMax(0.019);
       config->SetDeltaPhiMax(0.019);
     } else if (suffix == "34") {
       TrackCutsProtonDCA->SetPtRange(0.6, 4.05);
       TrackCutsAntiProtonDCA->SetPtRange(0.6, 4.05);
       TrackCutsProtonDCA->SetEtaRange(-0.83, 0.83);
       TrackCutsAntiProtonDCA->SetEtaRange(-0.83, 0.83);
       TrackCutsProtonDCA->SetNClsTPC(70);
       TrackCutsAntiProtonDCA->SetNClsTPC(70);
       TrackCutsDeuteronDCA->SetPtRange(0.6, 4.05);
       TrackCutsAntiDeuteronDCA->SetPtRange(0.6, 4.05);
       TrackCutsDeuteronDCA->SetEtaRange(-0.83, 0.83);
       TrackCutsAntiDeuteronDCA->SetEtaRange(-0.83, 0.83);
       TrackCutsDeuteronDCA->SetNClsTPC(70);
       TrackCutsAntiDeuteronDCA->SetNClsTPC(70);
       config->SetDeltaEtaMax(0.019);
       config->SetDeltaPhiMax(0.019);
     } else if (suffix == "35") {
       TrackCutsProtonDCA->SetEtaRange(-0.83, 0.83);
       TrackCutsAntiProtonDCA->SetEtaRange(-0.83, 0.83);
       TrackCutsProtonDCA->SetPID(AliPID::kProton, 0.75, 3.5);
       TrackCutsAntiProtonDCA->SetPID(AliPID::kProton, 0.75, 3.5);
       TrackCutsProtonDCA->SetNClsTPC(90);
       TrackCutsAntiProtonDCA->SetNClsTPC(90);
       TrackCutsDeuteronDCA->SetEtaRange(-0.83, 0.83);
       TrackCutsAntiDeuteronDCA->SetEtaRange(-0.83, 0.83);
       TrackCutsDeuteronDCA->SetPID(AliPID::kDeuteron,1.4, 3.3);
       TrackCutsAntiDeuteronDCA->SetPID(AliPID::kDeuteron, 1.4, 3.3);
       TrackCutsDeuteronDCA->SetNClsTPC(90);
       TrackCutsAntiDeuteronDCA->SetNClsTPC(90);
       config->SetDeltaEtaMax(0.019);
       config->SetDeltaPhiMax(0.019);
     } else if (suffix == "36") {
       TrackCutsProtonDCA->SetEtaRange(-0.77, 0.77);
       TrackCutsAntiProtonDCA->SetEtaRange(-0.77, 0.77);
       TrackCutsProtonDCA->SetNClsTPC(90);
       TrackCutsAntiProtonDCA->SetNClsTPC(90);
       TrackCutsDeuteronDCA->SetEtaRange(-0.77, 0.77);
       TrackCutsAntiDeuteronDCA->SetEtaRange(-0.77, 0.77);
       TrackCutsDeuteronDCA->SetNClsTPC(90);
       TrackCutsAntiDeuteronDCA->SetNClsTPC(90);
       config->SetDeltaEtaMax(0.019);
       config->SetDeltaPhiMax(0.019);
     } else if (suffix == "37") {
       TrackCutsProtonDCA->SetPtRange(0.4, 4.05);
       TrackCutsAntiProtonDCA->SetPtRange(0.4, 4.05);
       TrackCutsProtonDCA->SetPID(AliPID::kProton, 0.75, 3.5);
       TrackCutsAntiProtonDCA->SetPID(AliPID::kProton, 0.75, 3.5);
       TrackCutsProtonDCA->SetNClsTPC(90);
       TrackCutsAntiProtonDCA->SetNClsTPC(90);
       TrackCutsDeuteronDCA->SetPtRange(0.4, 4.05);
       TrackCutsAntiDeuteronDCA->SetPtRange(0.4, 4.05);
       TrackCutsDeuteronDCA->SetPID(AliPID::kDeuteron,1.4, 3.3);
       TrackCutsAntiDeuteronDCA->SetPID(AliPID::kDeuteron, 1.4, 3.3);
       TrackCutsDeuteronDCA->SetNClsTPC(90);
       TrackCutsAntiDeuteronDCA->SetNClsTPC(90);
       config->SetDeltaEtaMax(0.019);
       config->SetDeltaPhiMax(0.019);
     } else if (suffix == "38") {
       TrackCutsProtonDCA->SetPtRange(0.6, 4.05);
       TrackCutsAntiProtonDCA->SetPtRange(0.6, 4.05);
       TrackCutsProtonDCA->SetNClsTPC(70);
       TrackCutsAntiProtonDCA->SetNClsTPC(70);
       TrackCutsDeuteronDCA->SetPtRange(0.6, 4.05);
       TrackCutsAntiDeuteronDCA->SetPtRange(0.6, 4.05);
       TrackCutsDeuteronDCA->SetNClsTPC(70);
       TrackCutsAntiDeuteronDCA->SetNClsTPC(70);
     } else if (suffix == "39") {
       TrackCutsProtonDCA->SetPtRange(0.4, 4.05);
       TrackCutsAntiProtonDCA->SetPtRange(0.4, 4.05);
       TrackCutsProtonDCA->SetEtaRange(-0.77, 0.77);
       TrackCutsAntiProtonDCA->SetEtaRange(-0.77, 0.77);
       TrackCutsProtonDCA->SetPID(AliPID::kProton, 0.75, 3.5);
       TrackCutsAntiProtonDCA->SetPID(AliPID::kProton, 0.75, 3.5);
       TrackCutsDeuteronDCA->SetPtRange(0.4, 4.05);
       TrackCutsAntiDeuteronDCA->SetPtRange(0.4, 4.05);
       TrackCutsDeuteronDCA->SetEtaRange(-0.77, 0.77);
       TrackCutsAntiDeuteronDCA->SetEtaRange(-0.77, 0.77);
       TrackCutsDeuteronDCA->SetPID(AliPID::kDeuteron,1.4, 3.3);
       TrackCutsAntiDeuteronDCA->SetPID(AliPID::kDeuteron, 1.4, 3.3);
       config->SetDeltaEtaMax(0.015);
       config->SetDeltaPhiMax(0.015);
     } else if (suffix == "40") {
       TrackCutsProtonDCA->SetPtRange(0.6, 4.05);
       TrackCutsAntiProtonDCA->SetPtRange(0.6, 4.05);
       TrackCutsProtonDCA->SetPID(AliPID::kProton, 0.75, 2.5);
       TrackCutsAntiProtonDCA->SetPID(AliPID::kProton, 0.75, 2.5);
       TrackCutsDeuteronDCA->SetPtRange(0.6, 4.05);
       TrackCutsAntiDeuteronDCA->SetPtRange(0.6, 4.05);
       TrackCutsDeuteronDCA->SetPID(AliPID::kDeuteron,1.4, 2.3);
       TrackCutsAntiDeuteronDCA->SetPID(AliPID::kDeuteron, 1.4, 2.3);
     } else if (suffix == "41") {
       TrackCutsProtonDCA->SetPtRange(0.6, 4.05);
       TrackCutsAntiProtonDCA->SetPtRange(0.6, 4.05);
       TrackCutsProtonDCA->SetEtaRange(-0.83, 0.83);
       TrackCutsAntiProtonDCA->SetEtaRange(-0.83, 0.83);
       TrackCutsProtonDCA->SetPID(AliPID::kProton, 0.75, 3.5);
       TrackCutsAntiProtonDCA->SetPID(AliPID::kProton, 0.75, 3.5);
       TrackCutsDeuteronDCA->SetPtRange(0.6, 4.05);
       TrackCutsAntiDeuteronDCA->SetPtRange(0.6, 4.05);
       TrackCutsDeuteronDCA->SetEtaRange(-0.83, 0.83);
       TrackCutsAntiDeuteronDCA->SetEtaRange(-0.83, 0.83);
       TrackCutsDeuteronDCA->SetPID(AliPID::kDeuteron,1.4, 3.3);
       TrackCutsAntiDeuteronDCA->SetPID(AliPID::kDeuteron, 1.4, 3.3);
       config->SetDeltaEtaMax(0.019);
       config->SetDeltaPhiMax(0.019);
     } else if (suffix == "42") {
       TrackCutsProtonDCA->SetEtaRange(-0.83, 0.83);
       TrackCutsAntiProtonDCA->SetEtaRange(-0.83, 0.83);
       TrackCutsProtonDCA->SetPID(AliPID::kProton, 0.75, 3.5);
       TrackCutsAntiProtonDCA->SetPID(AliPID::kProton, 0.75, 3.5);
       TrackCutsDeuteronDCA->SetEtaRange(-0.83, 0.83);
       TrackCutsAntiDeuteronDCA->SetEtaRange(-0.83, 0.83);
       TrackCutsDeuteronDCA->SetPID(AliPID::kDeuteron,1.4, 3.3);
       TrackCutsAntiDeuteronDCA->SetPID(AliPID::kDeuteron, 1.4, 3.3);
       config->SetDeltaEtaMax(0.019);
       config->SetDeltaPhiMax(0.019);
     }else if (suffix == "43") {
       TrackCutsProtonDCA->SetDCAVtxZ(0.16);
       TrackCutsProtonDCA->SetDCAVtxXY(0.08);
       TrackCutsProtonDCA->SetDCAVtxZ(0.16);
       TrackCutsProtonDCA->SetDCAVtxXY(0.08);
       TrackCutsDeuteronDCA->SetDCAVtxZ(0.16);
       TrackCutsDeuteronDCA->SetDCAVtxXY(0.08);
       TrackCutsDeuteronDCA->SetDCAVtxZ(0.16);
       TrackCutsAntiDeuteronDCA->SetDCAVtxXY(0.08);
     }else if (suffix == "44") {
       TrackCutsProtonDCA->SetDCAVtxZ(0.24);
       TrackCutsProtonDCA->SetDCAVtxXY(0.12);
       TrackCutsProtonDCA->SetDCAVtxZ(0.24);
       TrackCutsProtonDCA->SetDCAVtxXY(0.12);
       TrackCutsDeuteronDCA->SetDCAVtxZ(0.24);
       TrackCutsDeuteronDCA->SetDCAVtxXY(0.12);
       TrackCutsDeuteronDCA->SetDCAVtxZ(0.24);
       TrackCutsAntiDeuteronDCA->SetDCAVtxXY(0.12);
     }
   }

  AliAnalysisTaskFemtoDreamDeuteron *task =
    new AliAnalysisTaskFemtoDreamDeuteron("FemtoDreamDefault", isMC);
  task->SelectCollisionCandidates(AliVEvent::kHighMultV0);
  if (!fullBlastQA) {
    task->SetRunTaskLightWeight(true);
  }
  if (SidebandStudy) {
    TrackCutsDeuteronDCA->SetPtRange(1.4, 4.05);
    TrackCutsAntiDeuteronDCA->SetPtRange(1.4, 4.05);
    TrackCutsDeuteronDCA->SetPID(AliPID::kDeuteron, 999.0);
    TrackCutsAntiDeuteronDCA->SetPID(AliPID::kDeuteron, 999.0);
    float up  = 0.0;
    float low = 0.0;
    if (suffix == "0"){//signal we keep fixed signal with tune from m^2 TPC+TOF
       up  = 3.0;
       low = -3.5;
    }else if (suffix == "1") {//lower sideband defaul
      up  = -4.1;
      low = -10.6;
    }else if (suffix == "2") {//lower sideband
      up  = -4.07;
      low = -10.6;
    }else if (suffix == "3") {//lower sideband
      up  = -4.13;
      low = -10.6;
    }else if (suffix == "4") {//lower sideband
      up  = -4.05;
      low = -10.6;
    }else if (suffix == "5") {//lower sideband
      up  = -4.15;
      low = -10.6;
    }else if (suffix == "6") {//lower sideband
      up  = -4.1;
      low = -10.3;
    }else if (suffix == "7") {//lower sideband
      up  = -4.1;
      low = -10.9;
    }else if (suffix == "8") {//lower sideband
      up  = -4.1;
      low = -10.0;
    }else if (suffix == "9") {//lower sideband
      up  = -4.1;
      low = -11.2;
    }else if (suffix == "10") {//upper sideband default
       up  = 10.1;
       low = 3.6;
    }else if (suffix == "11") {//upper sideband
      up  = 10.1;
      low = 3.57;
    }else if (suffix == "12") {//upper sideband
      up  = 10.1;
      low = 3.63;
    }else if (suffix == "13") {//upper sideband
      up  = 10.1;
      low = 3.55;
    }else if (suffix == "14") {//upper sideband
      up  = 10.1;
      low = 3.65;
    }else if (suffix == "15") {//upper sideband
      up  = 9.8;
      low = 3.6;
    }else if (suffix == "16") {//upper sideband
      up  = 10.4;
      low = 3.6;
    }else if (suffix == "17") {//upper sideband
      up  = 11.0;
      low = 3.6;
    }else if (suffix == "18") {//upper sideband
      up  = 9.8;
      low = 3.6;
    }
    task->SetSideband(up,low);
    task->SetSideband(up,low);
  }
  //Throw all our settings to the task
  task->SetEventCuts(evtCuts);
  task->SetTrackCutsDeuteronDCA(TrackCutsDeuteronDCA);
  task->SetTrackCutsDeuteronMass(TrackCutsDeuteronMass);
  task->SetTrackCutsAntiDeuteronDCA(TrackCutsAntiDeuteronDCA);
  task->SetTrackCutsAntiDeuteronMass(TrackCutsAntiDeuteronMass);
  task->SetTrackCutsProtonDCA(TrackCutsProtonDCA);
  task->SetTrackCutsAntiProtonDCA(TrackCutsAntiProtonDCA);
  task->SetCollectionConfig(config);

  mgr->AddTask(task);

  TString addon = "HM";

  TString file = AliAnalysisManager::GetCommonFileName();
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  mgr->ConnectInput(task, 0, cinput);

  TString EvtCutsName = Form("%sEvtCuts%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *coutputEvtCuts = mgr->CreateContainer(
        EvtCutsName.Data(), TList::Class(), AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), EvtCutsName.Data()));
  mgr->ConnectOutput(task, 1, coutputEvtCuts);

  TString TrackCutsName = Form("%sProtonDCA%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *coutputTrkCuts = mgr->CreateContainer(
        TrackCutsName.Data(), TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), TrackCutsName.Data()));
  mgr->ConnectOutput(task, 2, coutputTrkCuts);

  TString AntiTrackCutsName = Form("%sAntiProtonDCA%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *coutputAntiTrkCuts = mgr->CreateContainer(
        AntiTrackCutsName.Data(), TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), AntiTrackCutsName.Data()));
  mgr->ConnectOutput(task, 3, coutputAntiTrkCuts);

  TString TrackCutsDeuteronName = Form("%sDeuteronDCA%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *coutputTrkCutsDeuteron = mgr->CreateContainer(
        TrackCutsDeuteronName.Data(), TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), TrackCutsDeuteronName.Data()));
  mgr->ConnectOutput(task, 4, coutputTrkCutsDeuteron);

  TString AntiTrackCutsDeuteronName = Form("%sAntiDeuteronDCA%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *coutputAntiTrkCutsDeuteron = mgr->CreateContainer(
        AntiTrackCutsDeuteronName.Data(), TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), AntiTrackCutsDeuteronName.Data()));
  mgr->ConnectOutput(task, 5, coutputAntiTrkCutsDeuteron);

  TString TrackCutsDeuteronNoTOFName = Form("%sDeuteronMass%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *coutputTrkCutsDeuteronNoTOF = mgr->CreateContainer(
        TrackCutsDeuteronNoTOFName.Data(), TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), TrackCutsDeuteronNoTOFName.Data()));
  mgr->ConnectOutput(task, 6, coutputTrkCutsDeuteronNoTOF);

  TString AntiTrackCutsDeuteronNoTOFName = Form("%sAntiDeuteronMass%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *coutputAntiTrkCutsDeuteronNoTOF = mgr->CreateContainer(
        AntiTrackCutsDeuteronNoTOFName.Data(), TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), AntiTrackCutsDeuteronNoTOFName.Data()));
  mgr->ConnectOutput(task, 7, coutputAntiTrkCutsDeuteronNoTOF);

  AliAnalysisDataContainer *coutputResults;
  TString ResultsName = Form("%sResults%s", addon.Data(), suffix.Data());
  coutputResults = mgr->CreateContainer(ResultsName.Data(),
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

  if (isMC) {
    AliAnalysisDataContainer *coutputTrkCutsMC;
    TString TrkCutsMCName = Form("%sProtonMC%s", addon.Data(), suffix.Data());
    coutputTrkCutsMC = mgr->CreateContainer(
                         TrkCutsMCName.Data(),
                         TList::Class(),
                         AliAnalysisManager::kOutputContainer,
                         Form("%s:%s", file.Data(), TrkCutsMCName.Data()));
    mgr->ConnectOutput(task, 10, coutputTrkCutsMC);

    AliAnalysisDataContainer *coutputAntiTrkCutsMC;
    TString AntiTrkCutsMCName = Form("%sAntiProtonMC%s", addon.Data(),
                                     suffix.Data());
    coutputAntiTrkCutsMC = mgr->CreateContainer(
                             //@suppress("Invalid arguments") it works ffs
                             AntiTrkCutsMCName.Data(),
                             TList::Class(),
                             AliAnalysisManager::kOutputContainer,
                             Form("%s:%s", file.Data(), AntiTrkCutsMCName.Data()));
    mgr->ConnectOutput(task, 11, coutputAntiTrkCutsMC);

    AliAnalysisDataContainer *coutputv0CutsMC;
    TString v0CutsMCName = Form("%sDeuteronMC%s", addon.Data(), suffix.Data());
    coutputv0CutsMC = mgr->CreateContainer(
                        //@suppress("Invalid arguments") it works ffs
                        v0CutsMCName.Data(),
                        TList::Class(),
                        AliAnalysisManager::kOutputContainer,
                        Form("%s:%s", file.Data(), v0CutsMCName.Data()));
    mgr->ConnectOutput(task, 12, coutputv0CutsMC);

    AliAnalysisDataContainer *coutputAntiv0CutsMC;
    TString Antiv0CutsMCName = Form("%sAntiDeuteronMC%s", addon.Data(),
                                    suffix.Data());
    coutputAntiv0CutsMC = mgr->CreateContainer(
                            //@suppress("Invalid arguments") it works ffs
                            Antiv0CutsMCName.Data(),
                            TList::Class(),
                            AliAnalysisManager::kOutputContainer,
                            Form("%s:%s", file.Data(), Antiv0CutsMCName.Data()));
    mgr->ConnectOutput(task, 13, coutputAntiv0CutsMC);
  }
  return task;
}

