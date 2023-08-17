#include <vector>
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliFemtoDreamEventCuts.h"
#include "AliFemtoDreamTrackCuts.h"
#include "AliFemtoDreamCollConfig.h"
#include "AliAnalysisTaskNanoFemtoProtonPion.h"
#include "AliAnalysisTaskFemtoProtonPion.h"

AliAnalysisTaskSE* AddTaskNanoFemtoProtonPion(
    bool isNano = true, //0
    bool isMC = false, //1
    bool doOfficialFemto = true, //2
    TString trigger = "kHM", //3
    bool fullBlastQA = true,//4
    bool UseSphericityCut = true,//5
    float SphericityMinPt = 0.5, //6
    int PionFilterbit = 96, //7
    bool DoPairCleaning = false, //8
    bool DoAncestors = false, //9
    bool RemoveMCResonances = true, //10
    bool RemoveMCResonanceDaughters = true, //11
    bool DoInvMass = false, //12
    bool DoResonanceLorentzFactor = false, //13
    bool DoFinemTBinning = false, //14
    int mTBinningChoice = 1, //15
    bool DoFunWithPhaseSpace = false, //16
    float pTOnepTTwokStarCutOff = 3., //17
    bool DoKine = false, //18
    bool DoReco = false, //19
    const char *cutVariation = "0" //20
    ) {

  TString suffix = TString::Format("%s", cutVariation);
  bool isSystematic = false; 
  if(suffix != "0"){
    isSystematic = true; 
  }
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    printf("No analysis manager to connect to!\n");
    return nullptr;
  }
  if (!mgr->GetInputEventHandler()) {
    printf("This taskNano requires an input event handler!\n");
    return nullptr;
  }

  //Event cuts ------------------------------------------------------------------------
  AliFemtoDreamEventCuts *evtCuts = AliFemtoDreamEventCuts::StandardCutsRun2();
  evtCuts->CleanUpMult(false, false, false, true);

  if (UseSphericityCut){
    float SpherDown = 0.7;
    evtCuts->SetSphericityCuts(SpherDown, 1.0, SphericityMinPt); 
  }

  //Proton & AntiProton cuts --------------------------------------------------------
  AliFemtoDreamTrackCuts *TrackCutsProton = AliFemtoDreamTrackCuts::PrimProtonCuts(
        isMC, true, false, false);
  TrackCutsProton->SetFilterBit(128);
  TrackCutsProton->SetCutCharge(1);

  AliFemtoDreamTrackCuts *TrackCutsAntiProton = AliFemtoDreamTrackCuts::PrimProtonCuts(
        isMC, true, false, false);
  TrackCutsAntiProton->SetFilterBit(128);
  TrackCutsAntiProton->SetCutCharge(-1);

  //Pion & AntiPion Track cuts -------------------------------------------------------
  AliFemtoDreamTrackCuts *TrackCutsPion = NULL;
  AliFemtoDreamTrackCuts *TrackCutsAntiPion = NULL;
 
  TrackCutsPion = AliFemtoDreamTrackCuts::PrimPionCuts(isMC, true, false, false);
  TrackCutsPion->SetFilterBit(PionFilterbit);
  TrackCutsPion->SetCutCharge(1);

  if(isMC && PionFilterbit == 128){ //for MC template fits
    TrackCutsPion->CheckParticleMothers(true);
    TrackCutsPion->SetPlotDCADist(true);
    TrackCutsPion->SetFillQALater(false);
  }

  TrackCutsAntiPion = AliFemtoDreamTrackCuts::PrimPionCuts(isMC, true, false, false);
  TrackCutsAntiPion->SetFilterBit(PionFilterbit);
  TrackCutsAntiPion->SetCutCharge(-1);
  
  if(isMC && PionFilterbit == 128){ //for MC template fits
    TrackCutsAntiPion->CheckParticleMothers(true);
    TrackCutsAntiPion->SetPlotDCADist(true);
    TrackCutsAntiPion->SetFillQALater(false);
  }

  //Set-up output ------------------------------------------------------------------------
  std::vector<int> PDGParticles;
  PDGParticles.push_back(2212); 
  PDGParticles.push_back(-2212); 
  PDGParticles.push_back(211); 
  PDGParticles.push_back(-211); 

  std::vector<bool> closeRejection;
  std::vector<float> mTBins;

  if(mTBinningChoice == 0){ //old Binning
    mTBins.push_back(0.53); 
    mTBins.push_back(0.7); 
    mTBins.push_back(0.8); 
    mTBins.push_back(1.0); 
    mTBins.push_back(1.2); 
    mTBins.push_back(1.5); 
    mTBins.push_back(2.0);
    if(DoFinemTBinning){
     mTBins.push_back(3.0); 
    } 
    mTBins.push_back(4.0); 
  } else if (mTBinningChoice == 1) { //new Binning
    //{0.53, 0.75, 0.95, 1.2, 1.5, 2.0, 2.5, 4.0}
    mTBins.push_back(0.53); 
    mTBins.push_back(0.75); 
    mTBins.push_back(0.95); 
    mTBins.push_back(1.2); 
    mTBins.push_back(1.5); 
    mTBins.push_back(2.0); 
    mTBins.push_back(2.5);
    mTBins.push_back(4.0); 
  } else {
    mTBins.push_back(0.53); 
    mTBins.push_back(0.75); 
    mTBins.push_back(0.85); 
    mTBins.push_back(0.95); 
    mTBins.push_back(1.05); 
    mTBins.push_back(1.2); 
    mTBins.push_back(1.35); 
    mTBins.push_back(1.5); 
    mTBins.push_back(1.75); 
    mTBins.push_back(2.0); 
    mTBins.push_back(2.25);
    mTBins.push_back(2.5);
    mTBins.push_back(4.0); 
  }
  
  std::vector<int> pairQA;
  //pairs: 
  // pp             0
  // p bar p        1
  // p pi+          2
  // p pi-          3
  // bar p bar p    4
  // bar p pi+      5
  // bar p pi-      6
  // pi+ pi+        7
  // pi+ pi-        8
  // pi- pi-        9
  const int nPairs = 10;
  std::vector<int> NBins;
  std::vector<float> kMin;   //minimum k* value
  std::vector<float> kMax; //maximum k* value

  for (int i = 0; i < nPairs; ++i) {
    closeRejection.push_back(false);
    pairQA.push_back(0);
    NBins.push_back(750);
    kMin.push_back(0.);
    kMax.push_back(3.);
  }

  closeRejection[0] = true;  // pp
  closeRejection[2] = true;  // ppi+
  closeRejection[3] = false;  // ppi-
  closeRejection[4] = true;  // barp barp
  closeRejection[5] = false;  // barp pi+
  closeRejection[6] = true;  // barp pi-
  closeRejection[7] = true;  // pi+pi+
  closeRejection[9] = true;  // pi-pi-

  /*if(DoKine){
    closeRejection[0] = false;  // pp
    closeRejection[2] = false;  // ppi+
    closeRejection[3] = false;  // ppi-
    closeRejection[4] = false;  // barp barp
    closeRejection[5] = false;  // barp pi+
    closeRejection[6] = false;  // barp pi-
    closeRejection[7] = false;  // pi+pi+
    closeRejection[9] = false;  // pi-pi-
  }*/

  pairQA[0] = 11;  // pp
  pairQA[2] = 11;  // ppi+
  pairQA[3] = 11;  // ppi-
  pairQA[4] = 11;  // barp barp
  pairQA[5] = 11;  // barp pi+
  pairQA[6] = 11;  // barp pi-

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
  MultBins.push_back(84);
  MultBins.push_back(88);
  MultBins.push_back(92);
  MultBins.push_back(96);
  MultBins.push_back(100);

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
  config->SetDeltaPhiMax(0.04); // and here you set the actual values 
  config->SetMixingDepth(10);
  config->SetmTBins(mTBins);
  config->SetDomTMultBinning(true);
  config->SetmTBinning(true);
  config->SetMultiplicityEstimator(AliFemtoDreamEvent::kRef08);

  if (isMC) {
    config->SetMomentumResolution(true);
    if(DoAncestors){
      config->SetAncestors(true);
      config->GetDoAncestorsPlots();
    }
  }
  if (fullBlastQA) {
    // config->SetkTBinning(true);
    config->SetPtQA(true);
    config->SetdPhidEtaPlots(true);
    config->SetdPhidEtaPlotsSmallK(true);
    config->SetPhiEtaBinnign(true);
  } else {
    evtCuts->SetMinimalBooking(true);
    TrackCutsProton->SetMinimalBooking(true);
    TrackCutsAntiProton->SetMinimalBooking(true);
    TrackCutsPion->SetMinimalBooking(true);
    TrackCutsAntiPion->SetMinimalBooking(true);
    config->SetMinimalBookingME(true);
    config->SetMinimalBookingSample(true);
  }

  /*if(DoKine){
    config->SetdPhidEtaPlots(false);
    config->SetdPhidEtaPlotsSmallK(false);
    config->SetPhiEtaBinnign(false);
  }*/

  if(DoFunWithPhaseSpace){
    config->SetpTOnepTTwokStarPlotsmT(true, pTOnepTTwokStarCutOff);
  } else {
    config->SetpTOnepTTwokStarPlotsmT(false, pTOnepTTwokStarCutOff);
  }

  //============================================================================================================= Systematics
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
   Float_t Pion_Eta_VarLow = 0.7;
   Float_t Pion_Eta_VarHigh = 0.9;
   Float_t Pion_Clusters_VarLow = 70;
   Float_t Pion_Clusters_VarHigh = 90;
   Float_t Pion_Sigma_VarLow = 2.7;
   Float_t Pion_Sigma_VarHigh = 3.3;

   Float_t DPhi_VarLow = 0.035;
   Float_t DPhi_VarHigh = 0.045;
   Float_t DEta_VarLow = 0.015;
   Float_t DEta_VarHigh = 0.019;

  if (suffix == "1"){
 
    TrackCutsProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarLow);
    TrackCutsAntiProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarLow);
    TrackCutsProton->SetEtaRange(-Proton_Eta_VarLow,Proton_Eta_VarLow);
    TrackCutsAntiProton->SetEtaRange(-Proton_Eta_VarLow,Proton_Eta_VarLow);
 
    TrackCutsPion->SetPID(AliPID::kPion, 0.5,Pion_Sigma_VarHigh);
    TrackCutsAntiPion->SetPID(AliPID::kPion, 0.5,Pion_Sigma_VarHigh);
    TrackCutsPion->SetEtaRange(-Pion_Eta_VarHigh,Pion_Eta_VarHigh);
    TrackCutsAntiPion->SetEtaRange(-Pion_Eta_VarHigh,Pion_Eta_VarHigh);
    TrackCutsPion->SetPtRange(Pion_pT_VarLow,4.0);
    TrackCutsAntiPion->SetPtRange(Pion_pT_VarLow,4.0);
 
    config->SetDeltaPhiMax(DPhi_VarHigh);
 
  } else if (suffix == "2"){
 
    TrackCutsProton->SetPtRange(Proton_pT_VarLow,4.05);
    TrackCutsAntiProton->SetPtRange(Proton_pT_VarLow,4.05);
    TrackCutsProton->SetNClsTPC(Proton_Clusters_VarLow);
    TrackCutsAntiProton->SetNClsTPC(Proton_Clusters_VarLow);
 
    TrackCutsPion->SetPID(AliPID::kPion, 0.5,Pion_Sigma_VarLow);
    TrackCutsAntiPion->SetPID(AliPID::kPion, 0.5,Pion_Sigma_VarLow);
    TrackCutsPion->SetEtaRange(-Pion_Eta_VarHigh,Pion_Eta_VarHigh);
    TrackCutsAntiPion->SetEtaRange(-Pion_Eta_VarHigh,Pion_Eta_VarHigh);
    TrackCutsPion->SetPtRange(Pion_pT_VarHigh,4.0);
    TrackCutsAntiPion->SetPtRange(Pion_pT_VarHigh,4.0);
    TrackCutsPion->SetNClsTPC(Pion_Clusters_VarLow);
    TrackCutsAntiPion->SetNClsTPC(Pion_Clusters_VarLow);
 
    config->SetDeltaPhiMax(DPhi_VarHigh);
 
  } else if (suffix == "3"){
 
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
 
  } else if (suffix == "4"){
 
    TrackCutsProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarLow);
    TrackCutsAntiProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarLow);
    TrackCutsProton->SetPtRange(Proton_pT_VarHigh,4.05);
    TrackCutsAntiProton->SetPtRange(Proton_pT_VarHigh,4.05);
    TrackCutsProton->SetNClsTPC(Proton_Clusters_VarLow);
    TrackCutsAntiProton->SetNClsTPC(Proton_Clusters_VarLow);
 
  } else if (suffix == "5"){
 
    TrackCutsProton->SetEtaRange(-Proton_Eta_VarHigh,Proton_Eta_VarHigh);
    TrackCutsAntiProton->SetEtaRange(-Proton_Eta_VarHigh,Proton_Eta_VarHigh);
 
    TrackCutsPion->SetPID(AliPID::kPion, 0.5,Pion_Sigma_VarHigh);
    TrackCutsAntiPion->SetPID(AliPID::kPion, 0.5,Pion_Sigma_VarHigh);
    TrackCutsPion->SetEtaRange(-Pion_Eta_VarHigh,Pion_Eta_VarHigh);
    TrackCutsAntiPion->SetEtaRange(-Pion_Eta_VarHigh,Pion_Eta_VarHigh);
    TrackCutsPion->SetNClsTPC(Pion_Clusters_VarHigh);
    TrackCutsAntiPion->SetNClsTPC(Pion_Clusters_VarHigh);
 
    config->SetDeltaPhiMax(DPhi_VarHigh);
    config->SetDeltaEtaMax(DEta_VarHigh);
 
  } else if (suffix == "6"){
 
    TrackCutsProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarLow);
    TrackCutsAntiProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarLow);
    TrackCutsProton->SetEtaRange(-Proton_Eta_VarHigh,Proton_Eta_VarHigh);
    TrackCutsAntiProton->SetEtaRange(-Proton_Eta_VarHigh,Proton_Eta_VarHigh);
    TrackCutsProton->SetPtRange(Proton_pT_VarHigh,4.05);
    TrackCutsAntiProton->SetPtRange(Proton_pT_VarHigh,4.05);
    TrackCutsProton->SetNClsTPC(Proton_Clusters_VarHigh);
    TrackCutsAntiProton->SetNClsTPC(Proton_Clusters_VarHigh);
 
    TrackCutsPion->SetPID(AliPID::kPion, 0.5,Pion_Sigma_VarHigh);
    TrackCutsAntiPion->SetPID(AliPID::kPion, 0.5,Pion_Sigma_VarHigh);
    TrackCutsPion->SetEtaRange(-Pion_Eta_VarHigh,Pion_Eta_VarHigh);
    TrackCutsAntiPion->SetEtaRange(-Pion_Eta_VarHigh,Pion_Eta_VarHigh);
    TrackCutsPion->SetPtRange(Pion_pT_VarHigh,4.0);
    TrackCutsAntiPion->SetPtRange(Pion_pT_VarHigh,4.0);
    TrackCutsPion->SetNClsTPC(Pion_Clusters_VarHigh);
    TrackCutsAntiPion->SetNClsTPC(Pion_Clusters_VarHigh);
 
    config->SetDeltaPhiMax(DPhi_VarHigh);
    config->SetDeltaEtaMax(DEta_VarLow);
 
  } else if (suffix == "7"){
 
    TrackCutsProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarLow);
    TrackCutsAntiProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarLow);
    TrackCutsProton->SetPtRange(Proton_pT_VarHigh,4.05);
    TrackCutsAntiProton->SetPtRange(Proton_pT_VarHigh,4.05);
    TrackCutsProton->SetNClsTPC(Proton_Clusters_VarHigh);
    TrackCutsAntiProton->SetNClsTPC(Proton_Clusters_VarHigh);
 
    TrackCutsPion->SetPID(AliPID::kPion, 0.5,Pion_Sigma_VarLow);
    TrackCutsAntiPion->SetPID(AliPID::kPion, 0.5,Pion_Sigma_VarLow);
    TrackCutsPion->SetPtRange(Pion_pT_VarLow,4.0);
    TrackCutsAntiPion->SetPtRange(Pion_pT_VarLow,4.0);
 
    config->SetDeltaEtaMax(DEta_VarHigh);
 
  } else if (suffix == "8"){
 
    TrackCutsProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarHigh);
    TrackCutsAntiProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarHigh);
    TrackCutsProton->SetEtaRange(-Proton_Eta_VarLow,Proton_Eta_VarLow);
    TrackCutsAntiProton->SetEtaRange(-Proton_Eta_VarLow,Proton_Eta_VarLow);
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
    TrackCutsPion->SetNClsTPC(Pion_Clusters_VarHigh);
    TrackCutsAntiPion->SetNClsTPC(Pion_Clusters_VarHigh);
 
    config->SetDeltaEtaMax(DEta_VarHigh);
 
  } else if (suffix == "9"){
 
    TrackCutsProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarHigh);
    TrackCutsAntiProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarHigh);
    TrackCutsProton->SetEtaRange(-Proton_Eta_VarLow,Proton_Eta_VarLow);
    TrackCutsAntiProton->SetEtaRange(-Proton_Eta_VarLow,Proton_Eta_VarLow);
    TrackCutsProton->SetPtRange(Proton_pT_VarLow,4.05);
    TrackCutsAntiProton->SetPtRange(Proton_pT_VarLow,4.05);
    TrackCutsProton->SetNClsTPC(Proton_Clusters_VarLow);
    TrackCutsAntiProton->SetNClsTPC(Proton_Clusters_VarLow);
 
    TrackCutsPion->SetEtaRange(-Pion_Eta_VarLow,Pion_Eta_VarLow);
    TrackCutsAntiPion->SetEtaRange(-Pion_Eta_VarLow,Pion_Eta_VarLow);
    TrackCutsPion->SetPtRange(Pion_pT_VarLow,4.0);
    TrackCutsAntiPion->SetPtRange(Pion_pT_VarLow,4.0);
 
    config->SetDeltaPhiMax(DPhi_VarLow);
 
  } else if (suffix == "10"){
 
    TrackCutsProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarLow);
    TrackCutsAntiProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarLow);
    TrackCutsProton->SetEtaRange(-Proton_Eta_VarHigh,Proton_Eta_VarHigh);
    TrackCutsAntiProton->SetEtaRange(-Proton_Eta_VarHigh,Proton_Eta_VarHigh);
    TrackCutsProton->SetPtRange(Proton_pT_VarLow,4.05);
    TrackCutsAntiProton->SetPtRange(Proton_pT_VarLow,4.05);
 
    TrackCutsPion->SetEtaRange(-Pion_Eta_VarLow,Pion_Eta_VarLow);
    TrackCutsAntiPion->SetEtaRange(-Pion_Eta_VarLow,Pion_Eta_VarLow);
    TrackCutsPion->SetNClsTPC(Pion_Clusters_VarHigh);
    TrackCutsAntiPion->SetNClsTPC(Pion_Clusters_VarHigh);
 
    config->SetDeltaPhiMax(DPhi_VarLow);
    config->SetDeltaEtaMax(DEta_VarHigh);
 
  } else if (suffix == "11"){
 
    TrackCutsProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarLow);
    TrackCutsAntiProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarLow);
    TrackCutsProton->SetEtaRange(-Proton_Eta_VarLow,Proton_Eta_VarLow);
    TrackCutsAntiProton->SetEtaRange(-Proton_Eta_VarLow,Proton_Eta_VarLow);
    TrackCutsProton->SetPtRange(Proton_pT_VarHigh,4.05);
    TrackCutsAntiProton->SetPtRange(Proton_pT_VarHigh,4.05);
    TrackCutsProton->SetNClsTPC(Proton_Clusters_VarLow);
    TrackCutsAntiProton->SetNClsTPC(Proton_Clusters_VarLow);
 
    TrackCutsPion->SetPID(AliPID::kPion, 0.5,Pion_Sigma_VarLow);
    TrackCutsAntiPion->SetPID(AliPID::kPion, 0.5,Pion_Sigma_VarLow);
    TrackCutsPion->SetEtaRange(-Pion_Eta_VarHigh,Pion_Eta_VarHigh);
    TrackCutsAntiPion->SetEtaRange(-Pion_Eta_VarHigh,Pion_Eta_VarHigh);
    TrackCutsPion->SetNClsTPC(Pion_Clusters_VarLow);
    TrackCutsAntiPion->SetNClsTPC(Pion_Clusters_VarLow);
 
    config->SetDeltaEtaMax(DEta_VarLow);
 
  } else if (suffix == "12"){
 
    TrackCutsProton->SetEtaRange(-Proton_Eta_VarLow,Proton_Eta_VarLow);
    TrackCutsAntiProton->SetEtaRange(-Proton_Eta_VarLow,Proton_Eta_VarLow);
    TrackCutsProton->SetNClsTPC(Proton_Clusters_VarHigh);
    TrackCutsAntiProton->SetNClsTPC(Proton_Clusters_VarHigh);
 
    TrackCutsPion->SetPID(AliPID::kPion, 0.5,Pion_Sigma_VarHigh);
    TrackCutsAntiPion->SetPID(AliPID::kPion, 0.5,Pion_Sigma_VarHigh);
    TrackCutsPion->SetEtaRange(-Pion_Eta_VarLow,Pion_Eta_VarLow);
    TrackCutsAntiPion->SetEtaRange(-Pion_Eta_VarLow,Pion_Eta_VarLow);
    TrackCutsPion->SetPtRange(Pion_pT_VarHigh,4.0);
    TrackCutsAntiPion->SetPtRange(Pion_pT_VarHigh,4.0);
    TrackCutsPion->SetNClsTPC(Pion_Clusters_VarLow);
    TrackCutsAntiPion->SetNClsTPC(Pion_Clusters_VarLow);
 
    config->SetDeltaPhiMax(DPhi_VarLow);
 
  } else if (suffix == "13"){
 
    TrackCutsProton->SetPtRange(Proton_pT_VarHigh,4.05);
    TrackCutsAntiProton->SetPtRange(Proton_pT_VarHigh,4.05);
    TrackCutsProton->SetNClsTPC(Proton_Clusters_VarHigh);
    TrackCutsAntiProton->SetNClsTPC(Proton_Clusters_VarHigh);
 
    TrackCutsPion->SetEtaRange(-Pion_Eta_VarHigh,Pion_Eta_VarHigh);
    TrackCutsAntiPion->SetEtaRange(-Pion_Eta_VarHigh,Pion_Eta_VarHigh);
    TrackCutsPion->SetPtRange(Pion_pT_VarLow,4.0);
    TrackCutsAntiPion->SetPtRange(Pion_pT_VarLow,4.0);
    TrackCutsPion->SetNClsTPC(Pion_Clusters_VarLow);
    TrackCutsAntiPion->SetNClsTPC(Pion_Clusters_VarLow);
 
    config->SetDeltaEtaMax(DEta_VarLow);
 
  } else if (suffix == "14"){
 
    TrackCutsProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarHigh);
    TrackCutsAntiProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarHigh);
    TrackCutsProton->SetEtaRange(-Proton_Eta_VarLow,Proton_Eta_VarLow);
    TrackCutsAntiProton->SetEtaRange(-Proton_Eta_VarLow,Proton_Eta_VarLow);
    TrackCutsProton->SetNClsTPC(Proton_Clusters_VarLow);
    TrackCutsAntiProton->SetNClsTPC(Proton_Clusters_VarLow);
 
    TrackCutsPion->SetPID(AliPID::kPion, 0.5,Pion_Sigma_VarHigh);
    TrackCutsAntiPion->SetPID(AliPID::kPion, 0.5,Pion_Sigma_VarHigh);
    TrackCutsPion->SetPtRange(Pion_pT_VarHigh,4.0);
    TrackCutsAntiPion->SetPtRange(Pion_pT_VarHigh,4.0);
 
    config->SetDeltaEtaMax(DEta_VarLow);
 
  } else if (suffix == "15"){
 
    TrackCutsProton->SetEtaRange(-Proton_Eta_VarHigh,Proton_Eta_VarHigh);
    TrackCutsAntiProton->SetEtaRange(-Proton_Eta_VarHigh,Proton_Eta_VarHigh);
    TrackCutsProton->SetPtRange(Proton_pT_VarHigh,4.05);
    TrackCutsAntiProton->SetPtRange(Proton_pT_VarHigh,4.05);
    TrackCutsProton->SetNClsTPC(Proton_Clusters_VarHigh);
    TrackCutsAntiProton->SetNClsTPC(Proton_Clusters_VarHigh);
 
    TrackCutsPion->SetPID(AliPID::kPion, 0.5,Pion_Sigma_VarLow);
    TrackCutsAntiPion->SetPID(AliPID::kPion, 0.5,Pion_Sigma_VarLow);
    TrackCutsPion->SetPtRange(Pion_pT_VarLow,4.0);
    TrackCutsAntiPion->SetPtRange(Pion_pT_VarLow,4.0);
    TrackCutsPion->SetNClsTPC(Pion_Clusters_VarHigh);
    TrackCutsAntiPion->SetNClsTPC(Pion_Clusters_VarHigh);
 
    config->SetDeltaEtaMax(DEta_VarLow);
 
  } else if (suffix == "16"){
 
    TrackCutsProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarLow);
    TrackCutsAntiProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarLow);
    TrackCutsProton->SetPtRange(Proton_pT_VarLow,4.05);
    TrackCutsAntiProton->SetPtRange(Proton_pT_VarLow,4.05);
    TrackCutsProton->SetNClsTPC(Proton_Clusters_VarLow);
    TrackCutsAntiProton->SetNClsTPC(Proton_Clusters_VarLow);
 
    TrackCutsPion->SetPID(AliPID::kPion, 0.5,Pion_Sigma_VarLow);
    TrackCutsAntiPion->SetPID(AliPID::kPion, 0.5,Pion_Sigma_VarLow);
    TrackCutsPion->SetEtaRange(-Pion_Eta_VarLow,Pion_Eta_VarLow);
    TrackCutsAntiPion->SetEtaRange(-Pion_Eta_VarLow,Pion_Eta_VarLow);
    TrackCutsPion->SetPtRange(Pion_pT_VarLow,4.0);
    TrackCutsAntiPion->SetPtRange(Pion_pT_VarLow,4.0);
 
    config->SetDeltaPhiMax(DPhi_VarHigh);
    config->SetDeltaEtaMax(DEta_VarLow);
 
  } else if (suffix == "17"){
 
    TrackCutsProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarHigh);
    TrackCutsAntiProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarHigh);
    TrackCutsProton->SetPtRange(Proton_pT_VarLow,4.05);
    TrackCutsAntiProton->SetPtRange(Proton_pT_VarLow,4.05);
 
    TrackCutsPion->SetPID(AliPID::kPion, 0.5,Pion_Sigma_VarLow);
    TrackCutsAntiPion->SetPID(AliPID::kPion, 0.5,Pion_Sigma_VarLow);
    TrackCutsPion->SetEtaRange(-Pion_Eta_VarHigh,Pion_Eta_VarHigh);
    TrackCutsAntiPion->SetEtaRange(-Pion_Eta_VarHigh,Pion_Eta_VarHigh);
    TrackCutsPion->SetPtRange(Pion_pT_VarLow,4.0);
    TrackCutsAntiPion->SetPtRange(Pion_pT_VarLow,4.0);
    TrackCutsPion->SetNClsTPC(Pion_Clusters_VarHigh);
    TrackCutsAntiPion->SetNClsTPC(Pion_Clusters_VarHigh);
 
    config->SetDeltaPhiMax(DPhi_VarLow);
 
  } else if (suffix == "18"){
 
    TrackCutsProton->SetEtaRange(-Proton_Eta_VarHigh,Proton_Eta_VarHigh);
    TrackCutsAntiProton->SetEtaRange(-Proton_Eta_VarHigh,Proton_Eta_VarHigh);
    TrackCutsProton->SetPtRange(Proton_pT_VarLow,4.05);
    TrackCutsAntiProton->SetPtRange(Proton_pT_VarLow,4.05);
 
    TrackCutsPion->SetPID(AliPID::kPion, 0.5,Pion_Sigma_VarHigh);
    TrackCutsAntiPion->SetPID(AliPID::kPion, 0.5,Pion_Sigma_VarHigh);
    TrackCutsPion->SetEtaRange(-Pion_Eta_VarHigh,Pion_Eta_VarHigh);
    TrackCutsAntiPion->SetEtaRange(-Pion_Eta_VarHigh,Pion_Eta_VarHigh);
    TrackCutsPion->SetPtRange(Pion_pT_VarLow,4.0);
    TrackCutsAntiPion->SetPtRange(Pion_pT_VarLow,4.0);
 
    config->SetDeltaPhiMax(DPhi_VarHigh);
 
  } else if (suffix == "19"){
 
    TrackCutsProton->SetEtaRange(-Proton_Eta_VarLow,Proton_Eta_VarLow);
    TrackCutsAntiProton->SetEtaRange(-Proton_Eta_VarLow,Proton_Eta_VarLow);
    TrackCutsProton->SetNClsTPC(Proton_Clusters_VarHigh);
    TrackCutsAntiProton->SetNClsTPC(Proton_Clusters_VarHigh);
 
    TrackCutsPion->SetEtaRange(-Pion_Eta_VarLow,Pion_Eta_VarLow);
    TrackCutsAntiPion->SetEtaRange(-Pion_Eta_VarLow,Pion_Eta_VarLow);
    TrackCutsPion->SetPtRange(Pion_pT_VarLow,4.0);
    TrackCutsAntiPion->SetPtRange(Pion_pT_VarLow,4.0);
    TrackCutsPion->SetNClsTPC(Pion_Clusters_VarLow);
    TrackCutsAntiPion->SetNClsTPC(Pion_Clusters_VarLow);
 
    config->SetDeltaPhiMax(DPhi_VarHigh);
    config->SetDeltaEtaMax(DEta_VarHigh);
 
  } else if (suffix == "20"){
 
    TrackCutsProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarLow);
    TrackCutsAntiProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarLow);
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
    TrackCutsPion->SetNClsTPC(Pion_Clusters_VarLow);
    TrackCutsAntiPion->SetNClsTPC(Pion_Clusters_VarLow);
 
    config->SetDeltaEtaMax(DEta_VarLow);
 
  } else if (suffix == "21"){
 
    TrackCutsProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarLow);
    TrackCutsAntiProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarLow);
    TrackCutsProton->SetEtaRange(-Proton_Eta_VarLow,Proton_Eta_VarLow);
    TrackCutsAntiProton->SetEtaRange(-Proton_Eta_VarLow,Proton_Eta_VarLow);
    TrackCutsProton->SetNClsTPC(Proton_Clusters_VarHigh);
    TrackCutsAntiProton->SetNClsTPC(Proton_Clusters_VarHigh);
 
    TrackCutsPion->SetPID(AliPID::kPion, 0.5,Pion_Sigma_VarHigh);
    TrackCutsAntiPion->SetPID(AliPID::kPion, 0.5,Pion_Sigma_VarHigh);
    TrackCutsPion->SetPtRange(Pion_pT_VarLow,4.0);
    TrackCutsAntiPion->SetPtRange(Pion_pT_VarLow,4.0);
 
    config->SetDeltaPhiMax(DPhi_VarHigh);
    config->SetDeltaEtaMax(DEta_VarLow);
 
  } else if (suffix == "22"){
 
    TrackCutsProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarLow);
    TrackCutsAntiProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarLow);
    TrackCutsProton->SetEtaRange(-Proton_Eta_VarHigh,Proton_Eta_VarHigh);
    TrackCutsAntiProton->SetEtaRange(-Proton_Eta_VarHigh,Proton_Eta_VarHigh);
    TrackCutsProton->SetPtRange(Proton_pT_VarLow,4.05);
    TrackCutsAntiProton->SetPtRange(Proton_pT_VarLow,4.05);
    TrackCutsProton->SetNClsTPC(Proton_Clusters_VarHigh);
    TrackCutsAntiProton->SetNClsTPC(Proton_Clusters_VarHigh);
 
    TrackCutsPion->SetPID(AliPID::kPion, 0.5,Pion_Sigma_VarHigh);
    TrackCutsAntiPion->SetPID(AliPID::kPion, 0.5,Pion_Sigma_VarHigh);
    TrackCutsPion->SetEtaRange(-Pion_Eta_VarHigh,Pion_Eta_VarHigh);
    TrackCutsAntiPion->SetEtaRange(-Pion_Eta_VarHigh,Pion_Eta_VarHigh);
    TrackCutsPion->SetPtRange(Pion_pT_VarLow,4.0);
    TrackCutsAntiPion->SetPtRange(Pion_pT_VarLow,4.0);
 
    config->SetDeltaPhiMax(DPhi_VarLow);
 
  } else if (suffix == "23"){
 
    TrackCutsProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarLow);
    TrackCutsAntiProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarLow);
    TrackCutsProton->SetEtaRange(-Proton_Eta_VarHigh,Proton_Eta_VarHigh);
    TrackCutsAntiProton->SetEtaRange(-Proton_Eta_VarHigh,Proton_Eta_VarHigh);
    TrackCutsProton->SetPtRange(Proton_pT_VarHigh,4.05);
    TrackCutsAntiProton->SetPtRange(Proton_pT_VarHigh,4.05);
    TrackCutsProton->SetNClsTPC(Proton_Clusters_VarLow);
    TrackCutsAntiProton->SetNClsTPC(Proton_Clusters_VarLow);
 
    TrackCutsPion->SetEtaRange(-Pion_Eta_VarLow,Pion_Eta_VarLow);
    TrackCutsAntiPion->SetEtaRange(-Pion_Eta_VarLow,Pion_Eta_VarLow);
    TrackCutsPion->SetNClsTPC(Pion_Clusters_VarHigh);
    TrackCutsAntiPion->SetNClsTPC(Pion_Clusters_VarHigh);
 
    config->SetDeltaEtaMax(DEta_VarHigh);
 
  } else if (suffix == "24"){
 
    TrackCutsProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarLow);
    TrackCutsAntiProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarLow);
    TrackCutsProton->SetEtaRange(-Proton_Eta_VarLow,Proton_Eta_VarLow);
    TrackCutsAntiProton->SetEtaRange(-Proton_Eta_VarLow,Proton_Eta_VarLow);
    TrackCutsProton->SetPtRange(Proton_pT_VarLow,4.05);
    TrackCutsAntiProton->SetPtRange(Proton_pT_VarLow,4.05);
    TrackCutsProton->SetNClsTPC(Proton_Clusters_VarLow);
    TrackCutsAntiProton->SetNClsTPC(Proton_Clusters_VarLow);
 
    TrackCutsPion->SetEtaRange(-Pion_Eta_VarHigh,Pion_Eta_VarHigh);
    TrackCutsAntiPion->SetEtaRange(-Pion_Eta_VarHigh,Pion_Eta_VarHigh);
    TrackCutsPion->SetPtRange(Pion_pT_VarHigh,4.0);
    TrackCutsAntiPion->SetPtRange(Pion_pT_VarHigh,4.0);
    TrackCutsPion->SetNClsTPC(Pion_Clusters_VarHigh);
    TrackCutsAntiPion->SetNClsTPC(Pion_Clusters_VarHigh);
 
    config->SetDeltaPhiMax(DPhi_VarLow);
    config->SetDeltaEtaMax(DEta_VarHigh);
 
  } else if (suffix == "25"){
 
    TrackCutsProton->SetPtRange(Proton_pT_VarLow,4.05);
    TrackCutsAntiProton->SetPtRange(Proton_pT_VarLow,4.05);
    TrackCutsProton->SetNClsTPC(Proton_Clusters_VarHigh);
    TrackCutsAntiProton->SetNClsTPC(Proton_Clusters_VarHigh);
 
    TrackCutsPion->SetEtaRange(-Pion_Eta_VarLow,Pion_Eta_VarLow);
    TrackCutsAntiPion->SetEtaRange(-Pion_Eta_VarLow,Pion_Eta_VarLow);
    TrackCutsPion->SetNClsTPC(Pion_Clusters_VarHigh);
    TrackCutsAntiPion->SetNClsTPC(Pion_Clusters_VarHigh);
 
  } else if (suffix == "26"){
 
    TrackCutsProton->SetEtaRange(-Proton_Eta_VarLow,Proton_Eta_VarLow);
    TrackCutsAntiProton->SetEtaRange(-Proton_Eta_VarLow,Proton_Eta_VarLow);
 
    TrackCutsPion->SetPID(AliPID::kPion, 0.5,Pion_Sigma_VarLow);
    TrackCutsAntiPion->SetPID(AliPID::kPion, 0.5,Pion_Sigma_VarLow);
    TrackCutsPion->SetEtaRange(-Pion_Eta_VarLow,Pion_Eta_VarLow);
    TrackCutsAntiPion->SetEtaRange(-Pion_Eta_VarLow,Pion_Eta_VarLow);
    TrackCutsPion->SetPtRange(Pion_pT_VarLow,4.0);
    TrackCutsAntiPion->SetPtRange(Pion_pT_VarLow,4.0);
    TrackCutsPion->SetNClsTPC(Pion_Clusters_VarHigh);
    TrackCutsAntiPion->SetNClsTPC(Pion_Clusters_VarHigh);
 
    config->SetDeltaPhiMax(DPhi_VarHigh);
    config->SetDeltaEtaMax(DEta_VarLow);
 
  } else if (suffix == "27"){
 
    TrackCutsProton->SetEtaRange(-Proton_Eta_VarHigh,Proton_Eta_VarHigh);
    TrackCutsAntiProton->SetEtaRange(-Proton_Eta_VarHigh,Proton_Eta_VarHigh);
    TrackCutsProton->SetNClsTPC(Proton_Clusters_VarHigh);
    TrackCutsAntiProton->SetNClsTPC(Proton_Clusters_VarHigh);
 
    TrackCutsPion->SetEtaRange(-Pion_Eta_VarHigh,Pion_Eta_VarHigh);
    TrackCutsAntiPion->SetEtaRange(-Pion_Eta_VarHigh,Pion_Eta_VarHigh);
    TrackCutsPion->SetPtRange(Pion_pT_VarHigh,4.0);
    TrackCutsAntiPion->SetPtRange(Pion_pT_VarHigh,4.0);
    TrackCutsPion->SetNClsTPC(Pion_Clusters_VarLow);
    TrackCutsAntiPion->SetNClsTPC(Pion_Clusters_VarLow);
 
    config->SetDeltaPhiMax(DPhi_VarHigh);
 
  } else if (suffix == "28"){
 
    TrackCutsProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarLow);
    TrackCutsAntiProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarLow);
    TrackCutsProton->SetEtaRange(-Proton_Eta_VarHigh,Proton_Eta_VarHigh);
    TrackCutsAntiProton->SetEtaRange(-Proton_Eta_VarHigh,Proton_Eta_VarHigh);
    TrackCutsProton->SetNClsTPC(Proton_Clusters_VarHigh);
    TrackCutsAntiProton->SetNClsTPC(Proton_Clusters_VarHigh);
 
    TrackCutsPion->SetPID(AliPID::kPion, 0.5,Pion_Sigma_VarLow);
    TrackCutsAntiPion->SetPID(AliPID::kPion, 0.5,Pion_Sigma_VarLow);
    TrackCutsPion->SetEtaRange(-Pion_Eta_VarLow,Pion_Eta_VarLow);
    TrackCutsAntiPion->SetEtaRange(-Pion_Eta_VarLow,Pion_Eta_VarLow);
    TrackCutsPion->SetPtRange(Pion_pT_VarLow,4.0);
    TrackCutsAntiPion->SetPtRange(Pion_pT_VarLow,4.0);
    TrackCutsPion->SetNClsTPC(Pion_Clusters_VarLow);
    TrackCutsAntiPion->SetNClsTPC(Pion_Clusters_VarLow);
 
    config->SetDeltaEtaMax(DEta_VarLow);
 
  } else if (suffix == "29"){
 
    TrackCutsProton->SetPtRange(Proton_pT_VarLow,4.05);
    TrackCutsAntiProton->SetPtRange(Proton_pT_VarLow,4.05);
    TrackCutsProton->SetNClsTPC(Proton_Clusters_VarHigh);
    TrackCutsAntiProton->SetNClsTPC(Proton_Clusters_VarHigh);
 
    TrackCutsPion->SetPID(AliPID::kPion, 0.5,Pion_Sigma_VarHigh);
    TrackCutsAntiPion->SetPID(AliPID::kPion, 0.5,Pion_Sigma_VarHigh);
    TrackCutsPion->SetEtaRange(-Pion_Eta_VarHigh,Pion_Eta_VarHigh);
    TrackCutsAntiPion->SetEtaRange(-Pion_Eta_VarHigh,Pion_Eta_VarHigh);
 
    config->SetDeltaPhiMax(DPhi_VarHigh);
    config->SetDeltaEtaMax(DEta_VarLow);
 
  } else if (suffix == "30"){
 
    TrackCutsProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarHigh);
    TrackCutsAntiProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarHigh);
    TrackCutsProton->SetEtaRange(-Proton_Eta_VarHigh,Proton_Eta_VarHigh);
    TrackCutsAntiProton->SetEtaRange(-Proton_Eta_VarHigh,Proton_Eta_VarHigh);
    TrackCutsProton->SetPtRange(Proton_pT_VarHigh,4.05);
    TrackCutsAntiProton->SetPtRange(Proton_pT_VarHigh,4.05);
    TrackCutsProton->SetNClsTPC(Proton_Clusters_VarLow);
    TrackCutsAntiProton->SetNClsTPC(Proton_Clusters_VarLow);
 
    TrackCutsPion->SetPID(AliPID::kPion, 0.5,Pion_Sigma_VarLow);
    TrackCutsAntiPion->SetPID(AliPID::kPion, 0.5,Pion_Sigma_VarLow);
    TrackCutsPion->SetPtRange(Pion_pT_VarHigh,4.0);
    TrackCutsAntiPion->SetPtRange(Pion_pT_VarHigh,4.0);
    TrackCutsPion->SetNClsTPC(Pion_Clusters_VarLow);
    TrackCutsAntiPion->SetNClsTPC(Pion_Clusters_VarLow);
 
    config->SetDeltaEtaMax(DEta_VarHigh);
 
  } else if (suffix == "31"){
 
    TrackCutsProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarLow);
    TrackCutsAntiProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarLow);
    TrackCutsProton->SetEtaRange(-Proton_Eta_VarHigh,Proton_Eta_VarHigh);
    TrackCutsAntiProton->SetEtaRange(-Proton_Eta_VarHigh,Proton_Eta_VarHigh);
    TrackCutsProton->SetNClsTPC(Proton_Clusters_VarHigh);
    TrackCutsAntiProton->SetNClsTPC(Proton_Clusters_VarHigh);
 
    TrackCutsPion->SetPID(AliPID::kPion, 0.5,Pion_Sigma_VarHigh);
    TrackCutsAntiPion->SetPID(AliPID::kPion, 0.5,Pion_Sigma_VarHigh);
    TrackCutsPion->SetEtaRange(-Pion_Eta_VarLow,Pion_Eta_VarLow);
    TrackCutsAntiPion->SetEtaRange(-Pion_Eta_VarLow,Pion_Eta_VarLow);
    TrackCutsPion->SetPtRange(Pion_pT_VarHigh,4.0);
    TrackCutsAntiPion->SetPtRange(Pion_pT_VarHigh,4.0);
    TrackCutsPion->SetNClsTPC(Pion_Clusters_VarLow);
    TrackCutsAntiPion->SetNClsTPC(Pion_Clusters_VarLow);
 
    config->SetDeltaPhiMax(DPhi_VarLow);
 
  } else if (suffix == "32"){
 
    TrackCutsProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarHigh);
    TrackCutsAntiProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarHigh);
    TrackCutsProton->SetNClsTPC(Proton_Clusters_VarLow);
    TrackCutsAntiProton->SetNClsTPC(Proton_Clusters_VarLow);
 
    TrackCutsPion->SetPtRange(Pion_pT_VarHigh,4.0);
    TrackCutsAntiPion->SetPtRange(Pion_pT_VarHigh,4.0);
 
    config->SetDeltaPhiMax(DPhi_VarLow);
    config->SetDeltaEtaMax(DEta_VarHigh);
 
  } else if (suffix == "33"){
 
    TrackCutsProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarLow);
    TrackCutsAntiProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarLow);
    TrackCutsProton->SetEtaRange(-Proton_Eta_VarLow,Proton_Eta_VarLow);
    TrackCutsAntiProton->SetEtaRange(-Proton_Eta_VarLow,Proton_Eta_VarLow);
    TrackCutsProton->SetNClsTPC(Proton_Clusters_VarLow);
    TrackCutsAntiProton->SetNClsTPC(Proton_Clusters_VarLow);
 
    TrackCutsPion->SetEtaRange(-Pion_Eta_VarHigh,Pion_Eta_VarHigh);
    TrackCutsAntiPion->SetEtaRange(-Pion_Eta_VarHigh,Pion_Eta_VarHigh);
    TrackCutsPion->SetPtRange(Pion_pT_VarLow,4.0);
    TrackCutsAntiPion->SetPtRange(Pion_pT_VarLow,4.0);
 
    config->SetDeltaEtaMax(DEta_VarHigh);
 
  } else if (suffix == "34"){
 
    TrackCutsProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarLow);
    TrackCutsAntiProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarLow);
 
    TrackCutsPion->SetPID(AliPID::kPion, 0.5,Pion_Sigma_VarHigh);
    TrackCutsAntiPion->SetPID(AliPID::kPion, 0.5,Pion_Sigma_VarHigh);
    TrackCutsPion->SetEtaRange(-Pion_Eta_VarLow,Pion_Eta_VarLow);
    TrackCutsAntiPion->SetEtaRange(-Pion_Eta_VarLow,Pion_Eta_VarLow);
    TrackCutsPion->SetNClsTPC(Pion_Clusters_VarLow);
    TrackCutsAntiPion->SetNClsTPC(Pion_Clusters_VarLow);
 
    config->SetDeltaEtaMax(DEta_VarHigh);
 
  } else if (suffix == "35"){
 
    TrackCutsProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarLow);
    TrackCutsAntiProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarLow);
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
 
    config->SetDeltaPhiMax(DPhi_VarLow);
    config->SetDeltaEtaMax(DEta_VarHigh);
 
  } else if (suffix == "36"){
 
    TrackCutsProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarHigh);
    TrackCutsAntiProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarHigh);
    TrackCutsProton->SetEtaRange(-Proton_Eta_VarHigh,Proton_Eta_VarHigh);
    TrackCutsAntiProton->SetEtaRange(-Proton_Eta_VarHigh,Proton_Eta_VarHigh);
    TrackCutsProton->SetPtRange(Proton_pT_VarLow,4.05);
    TrackCutsAntiProton->SetPtRange(Proton_pT_VarLow,4.05);
    TrackCutsProton->SetNClsTPC(Proton_Clusters_VarHigh);
    TrackCutsAntiProton->SetNClsTPC(Proton_Clusters_VarHigh);
 
    TrackCutsPion->SetPID(AliPID::kPion, 0.5,Pion_Sigma_VarHigh);
    TrackCutsAntiPion->SetPID(AliPID::kPion, 0.5,Pion_Sigma_VarHigh);
    TrackCutsPion->SetEtaRange(-Pion_Eta_VarLow,Pion_Eta_VarLow);
    TrackCutsAntiPion->SetEtaRange(-Pion_Eta_VarLow,Pion_Eta_VarLow);
    TrackCutsPion->SetPtRange(Pion_pT_VarHigh,4.0);
    TrackCutsAntiPion->SetPtRange(Pion_pT_VarHigh,4.0);
    TrackCutsPion->SetNClsTPC(Pion_Clusters_VarLow);
    TrackCutsAntiPion->SetNClsTPC(Pion_Clusters_VarLow);
 
    config->SetDeltaPhiMax(DPhi_VarLow);
    config->SetDeltaEtaMax(DEta_VarHigh);
 
  } else if (suffix == "37"){
 
    TrackCutsProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarLow);
    TrackCutsAntiProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarLow);
    TrackCutsProton->SetPtRange(Proton_pT_VarHigh,4.05);
    TrackCutsAntiProton->SetPtRange(Proton_pT_VarHigh,4.05);
    TrackCutsProton->SetNClsTPC(Proton_Clusters_VarHigh);
    TrackCutsAntiProton->SetNClsTPC(Proton_Clusters_VarHigh);
 
    TrackCutsPion->SetEtaRange(-Pion_Eta_VarLow,Pion_Eta_VarLow);
    TrackCutsAntiPion->SetEtaRange(-Pion_Eta_VarLow,Pion_Eta_VarLow);
    TrackCutsPion->SetPtRange(Pion_pT_VarLow,4.0);
    TrackCutsAntiPion->SetPtRange(Pion_pT_VarLow,4.0);
 
    config->SetDeltaPhiMax(DPhi_VarHigh);
 
  } else if (suffix == "38"){
 
    TrackCutsProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarLow);
    TrackCutsAntiProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarLow);
    TrackCutsProton->SetPtRange(Proton_pT_VarLow,4.05);
    TrackCutsAntiProton->SetPtRange(Proton_pT_VarLow,4.05);
    TrackCutsProton->SetNClsTPC(Proton_Clusters_VarHigh);
    TrackCutsAntiProton->SetNClsTPC(Proton_Clusters_VarHigh);
 
    TrackCutsPion->SetEtaRange(-Pion_Eta_VarHigh,Pion_Eta_VarHigh);
    TrackCutsAntiPion->SetEtaRange(-Pion_Eta_VarHigh,Pion_Eta_VarHigh);
    TrackCutsPion->SetPtRange(Pion_pT_VarHigh,4.0);
    TrackCutsAntiPion->SetPtRange(Pion_pT_VarHigh,4.0);
    TrackCutsPion->SetNClsTPC(Pion_Clusters_VarHigh);
    TrackCutsAntiPion->SetNClsTPC(Pion_Clusters_VarHigh);
 
    config->SetDeltaPhiMax(DPhi_VarHigh);
    config->SetDeltaEtaMax(DEta_VarLow);
 
  } else if (suffix == "39"){
 
    TrackCutsProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarLow);
    TrackCutsAntiProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarLow);
    TrackCutsProton->SetEtaRange(-Proton_Eta_VarLow,Proton_Eta_VarLow);
    TrackCutsAntiProton->SetEtaRange(-Proton_Eta_VarLow,Proton_Eta_VarLow);
    TrackCutsProton->SetPtRange(Proton_pT_VarHigh,4.05);
    TrackCutsAntiProton->SetPtRange(Proton_pT_VarHigh,4.05);
    TrackCutsProton->SetNClsTPC(Proton_Clusters_VarHigh);
    TrackCutsAntiProton->SetNClsTPC(Proton_Clusters_VarHigh);
 
    TrackCutsPion->SetPID(AliPID::kPion, 0.5,Pion_Sigma_VarHigh);
    TrackCutsAntiPion->SetPID(AliPID::kPion, 0.5,Pion_Sigma_VarHigh);
    TrackCutsPion->SetPtRange(Pion_pT_VarHigh,4.0);
    TrackCutsAntiPion->SetPtRange(Pion_pT_VarHigh,4.0);
    TrackCutsPion->SetNClsTPC(Pion_Clusters_VarHigh);
    TrackCutsAntiPion->SetNClsTPC(Pion_Clusters_VarHigh);
 
    config->SetDeltaPhiMax(DPhi_VarLow);
    config->SetDeltaEtaMax(DEta_VarLow);
 
  } else if (suffix == "40"){
 
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
    TrackCutsPion->SetEtaRange(-Pion_Eta_VarHigh,Pion_Eta_VarHigh);
    TrackCutsAntiPion->SetEtaRange(-Pion_Eta_VarHigh,Pion_Eta_VarHigh);
    TrackCutsPion->SetNClsTPC(Pion_Clusters_VarHigh);
    TrackCutsAntiPion->SetNClsTPC(Pion_Clusters_VarHigh);
 
    config->SetDeltaPhiMax(DPhi_VarHigh);
    config->SetDeltaEtaMax(DEta_VarLow);
 
  } else if (suffix == "41"){
 
    TrackCutsProton->SetEtaRange(-Proton_Eta_VarLow,Proton_Eta_VarLow);
    TrackCutsAntiProton->SetEtaRange(-Proton_Eta_VarLow,Proton_Eta_VarLow);
    TrackCutsProton->SetNClsTPC(Proton_Clusters_VarHigh);
    TrackCutsAntiProton->SetNClsTPC(Proton_Clusters_VarHigh);
 
    TrackCutsPion->SetPID(AliPID::kPion, 0.5,Pion_Sigma_VarLow);
    TrackCutsAntiPion->SetPID(AliPID::kPion, 0.5,Pion_Sigma_VarLow);
    TrackCutsPion->SetEtaRange(-Pion_Eta_VarLow,Pion_Eta_VarLow);
    TrackCutsAntiPion->SetEtaRange(-Pion_Eta_VarLow,Pion_Eta_VarLow);
 
    config->SetDeltaPhiMax(DPhi_VarLow);
    config->SetDeltaEtaMax(DEta_VarHigh);
 
  } else if (suffix == "42"){
 
    TrackCutsProton->SetEtaRange(-Proton_Eta_VarLow,Proton_Eta_VarLow);
    TrackCutsAntiProton->SetEtaRange(-Proton_Eta_VarLow,Proton_Eta_VarLow);
    TrackCutsProton->SetPtRange(Proton_pT_VarLow,4.05);
    TrackCutsAntiProton->SetPtRange(Proton_pT_VarLow,4.05);
 
    TrackCutsPion->SetEtaRange(-Pion_Eta_VarLow,Pion_Eta_VarLow);
    TrackCutsAntiPion->SetEtaRange(-Pion_Eta_VarLow,Pion_Eta_VarLow);
 
    config->SetDeltaPhiMax(DPhi_VarHigh);
    config->SetDeltaEtaMax(DEta_VarLow);
 
  } else if (suffix == "43"){
 
    TrackCutsProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarHigh);
    TrackCutsAntiProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarHigh);
    TrackCutsProton->SetPtRange(Proton_pT_VarLow,4.05);
    TrackCutsAntiProton->SetPtRange(Proton_pT_VarLow,4.05);
    TrackCutsProton->SetNClsTPC(Proton_Clusters_VarHigh);
    TrackCutsAntiProton->SetNClsTPC(Proton_Clusters_VarHigh);
 
    TrackCutsPion->SetPID(AliPID::kPion, 0.5,Pion_Sigma_VarHigh);
    TrackCutsAntiPion->SetPID(AliPID::kPion, 0.5,Pion_Sigma_VarHigh);
    TrackCutsPion->SetPtRange(Pion_pT_VarLow,4.0);
    TrackCutsAntiPion->SetPtRange(Pion_pT_VarLow,4.0);
    TrackCutsPion->SetNClsTPC(Pion_Clusters_VarHigh);
    TrackCutsAntiPion->SetNClsTPC(Pion_Clusters_VarHigh);
 
    config->SetDeltaEtaMax(DEta_VarLow);
 
  } else if (suffix == "44"){
 
    TrackCutsProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarLow);
    TrackCutsAntiProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarLow);
    TrackCutsProton->SetPtRange(Proton_pT_VarLow,4.05);
    TrackCutsAntiProton->SetPtRange(Proton_pT_VarLow,4.05);
    TrackCutsProton->SetNClsTPC(Proton_Clusters_VarLow);
    TrackCutsAntiProton->SetNClsTPC(Proton_Clusters_VarLow);
 
    TrackCutsPion->SetEtaRange(-Pion_Eta_VarLow,Pion_Eta_VarLow);
    TrackCutsAntiPion->SetEtaRange(-Pion_Eta_VarLow,Pion_Eta_VarLow);
    TrackCutsPion->SetPtRange(Pion_pT_VarLow,4.0);
    TrackCutsAntiPion->SetPtRange(Pion_pT_VarLow,4.0);
    TrackCutsPion->SetNClsTPC(Pion_Clusters_VarLow);
    TrackCutsAntiPion->SetNClsTPC(Pion_Clusters_VarLow);
 
    config->SetDeltaEtaMax(DEta_VarHigh);
 
  }

  //============================================================================================================= Output Management
  TString addon = "";
  if (trigger == "kINT7") {
    addon += "kINT7";
  } else if (trigger == "kHM") {
    addon += "HM";
  }

  TString file = AliAnalysisManager::GetCommonFileName();
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer(); 

  TString EvtCutsName = Form("%s_EvtCuts_%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *coutputEvtCuts = mgr->CreateContainer(
        EvtCutsName.Data(), TList::Class(), AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), EvtCutsName.Data()));
  

  TString TrackCutsName = Form("%s_TrackProton_%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *coutputTrkCuts = mgr->CreateContainer(
        TrackCutsName.Data(), TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), TrackCutsName.Data()));


  TString AntiTrackCutsName = Form("%s_AntiTrackProton_%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *coutputAntiTrkCuts = mgr->CreateContainer(
        AntiTrackCutsName.Data(), TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), AntiTrackCutsName.Data()));


  TString TrackCutsPionName = Form("%s_TrackPion_%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *coutputTrkCutsPion = mgr->CreateContainer(
        TrackCutsPionName.Data(), TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), TrackCutsPionName.Data()));
 

  TString AntiTrackCutsPionName = Form("%s_AntiTrackPion_%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *coutputAntiTrkCutsPion = mgr->CreateContainer(
        AntiTrackCutsPionName.Data(), TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), AntiTrackCutsPionName.Data()));


  AliAnalysisDataContainer *coutputResults;
  TString ResultsName = Form("%s_Results_%s", addon.Data(), suffix.Data());
  coutputResults = mgr->CreateContainer(ResultsName.Data(),
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


  AliAnalysisDataContainer *coutputResultsThreeD;
  TString ResultsThreeDName = Form("%s_ResultsThreeD_%s", addon.Data(), suffix.Data());
  coutputResultsThreeD = mgr->CreateContainer(ResultsThreeDName.Data(),
                     TList::Class(), AliAnalysisManager::kOutputContainer,
                     Form("%s:%s", file.Data(), ResultsThreeDName.Data()));
 
  AliAnalysisDataContainer *coutputTrkCutsMC;
  AliAnalysisDataContainer *coutputAntiTrkCutsMC;
  AliAnalysisDataContainer *coutputTrkCutsPionMC;
  AliAnalysisDataContainer *coutputTrkCutsAntiPionMC;

  if (isMC) {
    TString TrkCutsMCName = Form("%s_ProtonMC_%s", addon.Data(), suffix.Data());
    coutputTrkCutsMC = mgr->CreateContainer(
                         TrkCutsMCName.Data(),
                         TList::Class(),
                         AliAnalysisManager::kOutputContainer,
                         Form("%s:%s", file.Data(), TrkCutsMCName.Data()));

    TString AntiTrkCutsMCName = Form("%s_AntiProtonMC_%s", addon.Data(), suffix.Data());
    coutputAntiTrkCutsMC = mgr->CreateContainer(
                             //@suppress("Invalid arguments") it works ffs
                             AntiTrkCutsMCName.Data(),
                             TList::Class(),
                             AliAnalysisManager::kOutputContainer,
                             Form("%s:%s", file.Data(), AntiTrkCutsMCName.Data()));

    TString PionCutsMCName = Form("%s_PionMC_%s", addon.Data(), suffix.Data());
    coutputTrkCutsPionMC = mgr->CreateContainer(
                        //@suppress("Invalid arguments") it works ffs
                        PionCutsMCName.Data(),
                        TList::Class(),
                        AliAnalysisManager::kOutputContainer,
                        Form("%s:%s", file.Data(), PionCutsMCName.Data()));

    TString AntiPionCutsMCName = Form("%s_AntiPionMC_%s", addon.Data(), suffix.Data());
    coutputTrkCutsAntiPionMC = mgr->CreateContainer(
                            //@suppress("Invalid arguments") it works ffs
                            AntiPionCutsMCName.Data(),
                            TList::Class(),
                            AliAnalysisManager::kOutputContainer,
                            Form("%s:%s", file.Data(), AntiPionCutsMCName.Data()));
  }


  AliAnalysisTaskNanoFemtoProtonPion *taskNano;
  AliAnalysisTaskFemtoProtonPion *taskAOD;

  if(isNano){
    taskNano = new AliAnalysisTaskNanoFemtoProtonPion("FemtoDreamDefault", isMC);
    
    if (trigger == "kINT7") {
      taskNano->SelectCollisionCandidates(AliVEvent::kINT7);
      std::cout << "Added kINT7 Trigger \n";
    } else if (trigger == "kHM") {
      taskNano->SelectCollisionCandidates(AliVEvent::kHighMultV0);
      std::cout << "Added kHighMult Trigger \n";
    } else {
      std::cout << "=====================================================================" << std::endl;
      std::cout << "=====================================================================" << std::endl;
      std::cout << "Centrality Estimator not set, fix it else your Results will be empty!" << std::endl;
      std::cout << "=====================================================================" << std::endl;
      std::cout << "=====================================================================" << std::endl;
    }
    if (!fullBlastQA) {
      taskNano->SetRunTaskLightWeight(true);
    }
    taskNano->SetEventCuts(evtCuts);
    taskNano->SetTrackCutsPion(TrackCutsPion);
    taskNano->SetTrackCutsAntiPion(TrackCutsAntiPion);
    taskNano->SetTrackCutsProton(TrackCutsProton);
    taskNano->SetTrackCutsAntiProton(TrackCutsAntiProton);
    taskNano->SetCollectionConfig(config);
    taskNano->SetDoPairCleaning(DoPairCleaning);
    taskNano->SetDoOfficialFemto(doOfficialFemto); 
    taskNano->SetDoResonanceLorentzFactor(DoResonanceLorentzFactor);
    
    //Set-up for own looping & calculus -> needed for 3D studies
    //IMPORTANT: 0, 1, 2, 3 and the names has to correspond to the order given to the offical femto framework!!!!
    taskNano->SetCombinationInput("00 11 02 13 03 12"); //p-p barp-barp p-pion barp-barpion p-barpion barp-pion
    taskNano->SetClosePairRejectionInput("true true true true false false");
    taskNano->SetNameTagInput("Proton AntiProton Pion AntiPion");
    taskNano->SetDoOwnFemto(!doOfficialFemto); //Do own looping and calculus 
    taskNano->SetDoThreeDFemto(false); //No 3D femto for now
    taskNano->SetRunPlotMult(true);
    taskNano->SetRunPlotPhiTheta(fullBlastQA); 
    taskNano->SetDoAncestors(DoAncestors); //Does not affect official femto part
    taskNano->SetRemoveMCResonances(RemoveMCResonances, RemoveMCResonanceDaughters);
    taskNano->SetDoInvMassPlot(DoInvMass);
    
    mgr->AddTask(taskNano);
    mgr->ConnectInput(taskNano, 0, cinput);
    mgr->ConnectOutput(taskNano, 1, coutputEvtCuts);
    mgr->ConnectOutput(taskNano, 2, coutputTrkCuts);      
    mgr->ConnectOutput(taskNano, 3, coutputAntiTrkCuts);
    mgr->ConnectOutput(taskNano, 4, coutputTrkCutsPion);
    mgr->ConnectOutput(taskNano, 5, coutputAntiTrkCutsPion);
    mgr->ConnectOutput(taskNano, 6, coutputResults);
    mgr->ConnectOutput(taskNano, 7, coutputResultsQA);
    mgr->ConnectOutput(taskNano, 8, coutputResultsThreeD);
    if(isMC){
      mgr->ConnectOutput(taskNano, 9, coutputTrkCutsMC);
      mgr->ConnectOutput(taskNano, 10, coutputAntiTrkCutsMC);
      mgr->ConnectOutput(taskNano, 11, coutputTrkCutsPionMC);
      mgr->ConnectOutput(taskNano, 12, coutputTrkCutsAntiPionMC);
    }
  } else {
    taskAOD = new AliAnalysisTaskFemtoProtonPion("FemtoDreamDefault", isMC);
    
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
    if (!fullBlastQA) {
      taskAOD->SetRunTaskLightWeight(true);
    }
    taskAOD->SetEventCuts(evtCuts);
    taskAOD->SetTrackCutsPion(TrackCutsPion);
    taskAOD->SetTrackCutsAntiPion(TrackCutsAntiPion);
    taskAOD->SetTrackCutsProton(TrackCutsProton);
    taskAOD->SetTrackCutsAntiProton(TrackCutsAntiProton);
    taskAOD->SetCollectionConfig(config);
    taskAOD->SetDoPairCleaning(DoPairCleaning);
    taskAOD->SetDoOfficialFemto(doOfficialFemto); 
    taskAOD->SetDoResonanceLorentzFactor(DoResonanceLorentzFactor);

    taskAOD->SetDoKine(DoKine); 
    taskAOD->SetDoReco(DoReco);

    //Set-up for own looping & calculus -> needed for 3D studies
    //IMPORTANT: 0, 1, 2, 3 and the names has to correspond to the order given to the offical femto framework!!!!
    taskAOD->SetCombinationInput("00 11 02 13 03 12"); //p-p barp-barp p-pion barp-barpion p-barpion barp-pion
    taskAOD->SetClosePairRejectionInput("true true true true false false");
    taskAOD->SetNameTagInput("Proton AntiProton Pion AntiPion");
    taskAOD->SetDoOwnFemto(!doOfficialFemto); //Do own looping and calculus 
    taskAOD->SetDoThreeDFemto(false); //No 3D femto for now
    taskAOD->SetRunPlotMult(true);
    taskAOD->SetRunPlotPhiTheta(fullBlastQA); 
    taskAOD->SetDoAncestors(DoAncestors); //Does not affect official femto part
    taskAOD->SetRemoveMCResonances(RemoveMCResonances, RemoveMCResonanceDaughters);
    taskAOD->SetDoInvMassPlot(DoInvMass);
    
    mgr->AddTask(taskAOD);
    mgr->ConnectInput(taskAOD, 0, cinput);
    mgr->ConnectOutput(taskAOD, 1, coutputEvtCuts);
    mgr->ConnectOutput(taskAOD, 2, coutputTrkCuts);      
    mgr->ConnectOutput(taskAOD, 3, coutputAntiTrkCuts);
    mgr->ConnectOutput(taskAOD, 4, coutputTrkCutsPion);
    mgr->ConnectOutput(taskAOD, 5, coutputAntiTrkCutsPion);
    mgr->ConnectOutput(taskAOD, 6, coutputResults);
    mgr->ConnectOutput(taskAOD, 7, coutputResultsQA);
    mgr->ConnectOutput(taskAOD, 8, coutputResultsThreeD);
    if(isMC){
      mgr->ConnectOutput(taskAOD, 9, coutputTrkCutsMC);
      mgr->ConnectOutput(taskAOD, 10, coutputAntiTrkCutsMC);
      mgr->ConnectOutput(taskAOD, 11, coutputTrkCutsPionMC);
      mgr->ConnectOutput(taskAOD, 12, coutputTrkCutsAntiPionMC);
    }
  }
   
  if (isNano) {
    return taskNano;
  } else {
    return taskAOD;
  }
  
}

