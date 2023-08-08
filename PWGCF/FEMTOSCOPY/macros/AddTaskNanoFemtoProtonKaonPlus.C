#include <vector>
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliFemtoDreamEventCuts.h"
#include "AliFemtoDreamTrackCuts.h"
#include "AliFemtoDreamCollConfig.h"
#include "AliAnalysisTaskNanoFemtoProtonKaonPlus.h"

AliAnalysisTaskSE* AddTaskNanoFemtoProtonKaonPlus(
    bool isNano = true, //1
    bool isMC = false, //2
    TString trigger = "kHM", //3
    bool fullBlastQA = true,//4
    bool UseSphericityCut = true,//5
    float SphericityMinPt = 0.5, //6
    bool UseRamonaCut = false, //7
    int filterBit = 128, //8
    bool DoPairCleaning = false, //9
    bool DoAncestors = false, //10
    float DCAxy = 0.1, //11
    float DCAz = 0.2, //12
    bool doDCA = false, //13
    bool DopTOnepTTwo = false, //14
    float pTOnepTTwokStarCutOff = 3.0, //15
    float pThighK = 999.0, //16
    float pTlowp = 0.4, //17
    int whichmTbinning = 1, //18
    const char *cutVariation = "0" //19
    ) {

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

  //Event cuts ------------------------------------------------------------------------
  AliFemtoDreamEventCuts *evtCuts = AliFemtoDreamEventCuts::StandardCutsRun2();
  evtCuts->CleanUpMult(false, false, false, true);

  if (UseSphericityCut){
    float SpherDown = 0.7;
    evtCuts->SetSphericityCuts(SpherDown, 1.0, SphericityMinPt); 
  }


  //Track cuts ------------------------------------------------------------------------
  
  AliFemtoDreamTrackCuts *TrackPosKaonCuts =
      AliFemtoDreamTrackCuts::PrimKaonCuts(isMC, true, false, false);
  TrackPosKaonCuts->SetCutCharge(1);
  if(isMC){
    TrackPosKaonCuts->SetPlotDCADist(true);
  }
  if(doDCA){
    TrackPosKaonCuts->SetPtRange(0.15, 5.0);
  }
  if (!UseRamonaCut)
  {
    TrackPosKaonCuts->SetPIDkd(); // Oton
    TrackPosKaonCuts->SetFilterBit(filterBit);
    TrackPosKaonCuts->SetPtRange(0.15, pThighK);
    //Default: 0.15 to 999 for Oton and Ramona
  }
  else
  {
    TrackPosKaonCuts->SetFilterBit(filterBit);
    TrackPosKaonCuts->SetPIDkd(true, true); // Ramona
    TrackPosKaonCuts->SetDCAVtxZ(DCAz);
    TrackPosKaonCuts->SetDCAVtxXY(DCAxy);
    TrackPosKaonCuts->SetCutTPCCrossedRows(false, 0, 0);
    TrackPosKaonCuts->SetCutSharedCls(false);
    TrackPosKaonCuts->SetCutSmallestSig(false);
    TrackPosKaonCuts->SetRejLowPtPionsTOF(false);
    TrackPosKaonCuts->SetPtRange(0.15, pThighK);
  }

  AliFemtoDreamTrackCuts *TrackNegKaonCuts =
      AliFemtoDreamTrackCuts::PrimKaonCuts(isMC, true, false, false);
  TrackNegKaonCuts->SetCutCharge(-1);
  if(isMC){
    TrackNegKaonCuts->SetPlotDCADist(true);
  }
  if(doDCA){
    TrackNegKaonCuts->SetPtRange(0.15, 5.0);
  }
  if (!UseRamonaCut)
  {
    TrackNegKaonCuts->SetPIDkd(); // Oton
    TrackNegKaonCuts->SetFilterBit(filterBit);
    TrackNegKaonCuts->SetPtRange(0.15, pThighK);
  }
  else
  {
    TrackNegKaonCuts->SetFilterBit(filterBit);
    TrackNegKaonCuts->SetPIDkd(true, true); // Ramona
    TrackNegKaonCuts->SetDCAVtxZ(DCAz);
    TrackNegKaonCuts->SetDCAVtxXY(DCAxy);
    TrackNegKaonCuts->SetCutTPCCrossedRows(false, 0, 0);
    TrackNegKaonCuts->SetCutSharedCls(false);
    TrackNegKaonCuts->SetCutSmallestSig(false);
    TrackNegKaonCuts->SetRejLowPtPionsTOF(false);
    TrackNegKaonCuts->SetPtRange(0.15, pThighK);
  }

  //Proton and AntiProton cuts
  AliFemtoDreamTrackCuts *TrackCutsProton = AliFemtoDreamTrackCuts::PrimProtonCuts(
        isMC, true, false, false);
  TrackCutsProton->SetCutCharge(1);
  TrackCutsProton->SetPtRange(pTlowp, 4.05);
  if(UseRamonaCut){
    TrackCutsProton->SetDCAVtxZ(DCAz);
    TrackCutsProton->SetDCAVtxXY(DCAxy);
    TrackCutsProton->SetCutTPCCrossedRows(false, 0, 0);
    TrackCutsProton->SetCutSharedCls(false);
    TrackCutsProton->SetCutSmallestSig(false);
    TrackCutsProton->SetRejLowPtPionsTOF(false);
    TrackCutsProton->SetPtRange(pTlowp, 3.0);
    //Default: 0.4 to 3.0 for Ramona, 0.5 to 4.05 for Oton
  }

  AliFemtoDreamTrackCuts *TrackCutsAntiProton = AliFemtoDreamTrackCuts::PrimProtonCuts(
        isMC, true, false, false);
  TrackCutsAntiProton->SetCutCharge(-1);
  TrackCutsAntiProton->SetPtRange(pTlowp, 4.05);
  if(UseRamonaCut){
    TrackCutsAntiProton->SetDCAVtxZ(DCAz);
    TrackCutsAntiProton->SetDCAVtxXY(DCAxy);
    TrackCutsAntiProton->SetCutTPCCrossedRows(false, 0, 0);
    TrackCutsAntiProton->SetCutSharedCls(false);
    TrackCutsAntiProton->SetCutSmallestSig(false);
    TrackCutsAntiProton->SetRejLowPtPionsTOF(false);
    TrackCutsAntiProton->SetPtRange(pTlowp, 3.0);
  }

  //Set-up output ------------------------------------------------------------------------
  std::vector<int> PDGParticles;
  PDGParticles.push_back(2212); 
  PDGParticles.push_back(-2212); 
  PDGParticles.push_back(321); 
  PDGParticles.push_back(-321);


  //BRELOOM Will have to enter mT bins
  std::vector<bool> closeRejection;
  std::vector<float> mTBins;
  //BRELOOM correct mT bins
  if(whichmTbinning == 1){//new binning
    mTBins.push_back(0.7);
    mTBins.push_back(1.0);
    mTBins.push_back(1.2);
    mTBins.push_back(1.4);
    mTBins.push_back(1.5);
    mTBins.push_back(1.8);
    mTBins.push_back(2.0);
    mTBins.push_back(100.0);
  } else if(whichmTbinning == 2){//old binning
    mTBins.push_back(0.71); 
    mTBins.push_back(1.); 
    mTBins.push_back(1.2); 
    mTBins.push_back(1.4); 
    mTBins.push_back(1.7); 
    mTBins.push_back(1.9); 
    mTBins.push_back(10.0); 
  }

  std::vector<int> pairQA;
  //pairs: 
  // pp             0
  // p bar p        1
  // p K+          2
  // p K-          3
  // bar p bar p    4
  // bar p K+      5
  // bar p K-      6
  // K+ K+        7
  // K+ K-        8
  // K- K-        9
  const int nPairs = 10;
  std::vector<int> NBins;
  std::vector<float> kMin;   //minimum k* value
  std::vector<float> kMax; //maximum k* value

  for (int i = 0; i < nPairs; ++i) {
    closeRejection.push_back(false);
    pairQA.push_back(0);
    NBins.push_back(1500);
    kMin.push_back(0.);
    kMax.push_back(3.);
  }

  closeRejection[0] = true;  // pp
  closeRejection[2] = true;  // pK+
  closeRejection[3] = false;  // pK-
  closeRejection[4] = true;  // barp barp
  closeRejection[5] = false;  // barp K+
  closeRejection[6] = true;  // barp K-
  closeRejection[7] = true;  // K+K+
  closeRejection[9] = true;  // K-K-

  pairQA[0] = 11;  // pp
  pairQA[2] = 11;  // pK+
  pairQA[3] = 11;  // pK-
  pairQA[4] = 11;  // barp barp
  pairQA[5] = 11;  // barp K+
  pairQA[6] = 11;  // barp K-
  pairQA[7] = 11;  // K+K+
  pairQA[9] = 11;  // K-K-

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
  //BRELOOM Will have to set values here
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
    TrackPosKaonCuts->SetMinimalBooking(true);
    TrackNegKaonCuts->SetMinimalBooking(true);
    config->SetMinimalBookingME(true);
    config->SetMinimalBookingSample(true);
  }

  if(DopTOnepTTwo){
    config->SetpTOnepTTwokStarPlotsmT(true, pTOnepTTwokStarCutOff);
  } else {
    config->SetpTOnepTTwokStarPlotsmT(false, pTOnepTTwokStarCutOff);
  }

  //---------------------------------------------------- Systematics
  Float_t Proton_pT_VarLow = 0.4;
  Float_t Proton_pT_VarHigh = 0.6;
  Float_t Proton_Eta_VarLow = 0.77;
  Float_t Proton_Eta_VarHigh = 0.83;
  Float_t Proton_Clusters_VarLow = 70.;
  Float_t Proton_Clusters_VarHigh = 80.;
  Float_t Proton_Sigma_VarLow = 2.5;
  Float_t Proton_Sigma_VarHigh = 3.5;

  Float_t Kaon_pT_VarLow = 0.1;
  Float_t Kaon_pT_VarHigh = 0.2;
  Float_t Kaon_Eta_VarLow = 0.75;
  Float_t Kaon_Eta_VarHigh = 0.85;
  Float_t Kaon_Clusters_VarLow = 70;
  Float_t Kaon_Clusters_VarHigh = 90;

  Float_t Kaon_SigmaTPC_Default = 3.; 
  Float_t Kaon_SigmaTPC_VarLow = 2.7;
  Float_t Kaon_SigmaTPC_VarHigh = 3.3;

  Float_t Kaon_SigmaCombined_Default = 3.; 
  Float_t Kaon_SigmaCombined_VarLow = 2.7;
  Float_t Kaon_SigmaCombined_VarHigh = 3.3;

  Float_t Kaon_SigmaExclusion_Default = 3.; 
  Float_t Kaon_SigmaExclusion_VarLow = 3.3;
  Float_t Kaon_SigmaExclusion_VarHigh = 2.7;
  
  Float_t DPhi_VarLow = 0.035; 
  Float_t DPhi_VarHigh = 0.045;
  Float_t DEta_VarLow = 0.015;
  Float_t DEta_VarHigh = 0.019;

 if (suffix == "1"){
 
    TrackCutsProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarHigh);
    TrackCutsAntiProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarHigh);
    TrackCutsProton->SetEtaRange(-Proton_Eta_VarHigh,Proton_Eta_VarHigh);
    TrackCutsAntiProton->SetEtaRange(-Proton_Eta_VarHigh,Proton_Eta_VarHigh);
    TrackCutsProton->SetPtRange(Proton_pT_VarHigh,4.05);
    TrackCutsAntiProton->SetPtRange(Proton_pT_VarHigh,4.05);
 
    TrackPosKaonCuts->SetPIDkd(true, UseRamonaCut,Kaon_SigmaCombined_Default, Kaon_SigmaTPC_VarLow, Kaon_SigmaExclusion_Default);
    TrackNegKaonCuts->SetPIDkd(true, UseRamonaCut,Kaon_SigmaCombined_Default, Kaon_SigmaTPC_VarLow, Kaon_SigmaExclusion_Default);
    TrackPosKaonCuts->SetEtaRange(-Kaon_Eta_VarLow,Kaon_Eta_VarLow);
    TrackNegKaonCuts->SetEtaRange(-Kaon_Eta_VarLow,Kaon_Eta_VarLow);
 
 
  } else if (suffix == "2"){
 
    TrackCutsProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarHigh);
    TrackCutsAntiProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarHigh);
    TrackCutsProton->SetPtRange(Proton_pT_VarLow,4.05);
    TrackCutsAntiProton->SetPtRange(Proton_pT_VarLow,4.05);
 
    TrackPosKaonCuts->SetPIDkd(true, UseRamonaCut,Kaon_SigmaCombined_VarLow, Kaon_SigmaTPC_Default, Kaon_SigmaExclusion_VarLow);
    TrackNegKaonCuts->SetPIDkd(true, UseRamonaCut,Kaon_SigmaCombined_VarLow, Kaon_SigmaTPC_Default, Kaon_SigmaExclusion_VarLow);
    TrackPosKaonCuts->SetEtaRange(-Kaon_Eta_VarHigh,Kaon_Eta_VarHigh);
    TrackNegKaonCuts->SetEtaRange(-Kaon_Eta_VarHigh,Kaon_Eta_VarHigh);
    TrackPosKaonCuts->SetPtRange(Kaon_pT_VarHigh,pThighK);
    TrackNegKaonCuts->SetPtRange(Kaon_pT_VarHigh,pThighK);
    TrackPosKaonCuts->SetNClsTPC(Kaon_Clusters_VarLow);
    TrackNegKaonCuts->SetNClsTPC(Kaon_Clusters_VarLow);
 
    config->SetDeltaEtaMax(DEta_VarHigh);
 
  } else if (suffix == "3"){
 
    TrackCutsProton->SetPtRange(Proton_pT_VarLow,4.05);
    TrackCutsAntiProton->SetPtRange(Proton_pT_VarLow,4.05);
    TrackCutsProton->SetNClsTPC(Proton_Clusters_VarHigh);
    TrackCutsAntiProton->SetNClsTPC(Proton_Clusters_VarHigh);
 
    TrackPosKaonCuts->SetPIDkd(true, UseRamonaCut,Kaon_SigmaCombined_VarLow, Kaon_SigmaTPC_Default, Kaon_SigmaExclusion_VarHigh);
    TrackNegKaonCuts->SetPIDkd(true, UseRamonaCut,Kaon_SigmaCombined_VarLow, Kaon_SigmaTPC_Default, Kaon_SigmaExclusion_VarHigh);
    TrackPosKaonCuts->SetEtaRange(-Kaon_Eta_VarLow,Kaon_Eta_VarLow);
    TrackNegKaonCuts->SetEtaRange(-Kaon_Eta_VarLow,Kaon_Eta_VarLow);
    TrackPosKaonCuts->SetPtRange(Kaon_pT_VarLow,pThighK);
    TrackNegKaonCuts->SetPtRange(Kaon_pT_VarLow,pThighK);
    TrackPosKaonCuts->SetNClsTPC(Kaon_Clusters_VarLow);
    TrackNegKaonCuts->SetNClsTPC(Kaon_Clusters_VarLow);
 
    config->SetDeltaEtaMax(DEta_VarLow);
 
  } else if (suffix == "4"){
 
    TrackCutsProton->SetEtaRange(-Proton_Eta_VarHigh,Proton_Eta_VarHigh);
    TrackCutsAntiProton->SetEtaRange(-Proton_Eta_VarHigh,Proton_Eta_VarHigh);
    TrackCutsProton->SetPtRange(Proton_pT_VarLow,4.05);
    TrackCutsAntiProton->SetPtRange(Proton_pT_VarLow,4.05);
    TrackCutsProton->SetNClsTPC(Proton_Clusters_VarHigh);
    TrackCutsAntiProton->SetNClsTPC(Proton_Clusters_VarHigh);
 
    TrackPosKaonCuts->SetPIDkd(true, UseRamonaCut,Kaon_SigmaCombined_VarHigh, Kaon_SigmaTPC_Default, Kaon_SigmaExclusion_VarHigh);
    TrackNegKaonCuts->SetPIDkd(true, UseRamonaCut,Kaon_SigmaCombined_VarHigh, Kaon_SigmaTPC_Default, Kaon_SigmaExclusion_VarHigh);
    TrackPosKaonCuts->SetEtaRange(-Kaon_Eta_VarHigh,Kaon_Eta_VarHigh);
    TrackNegKaonCuts->SetEtaRange(-Kaon_Eta_VarHigh,Kaon_Eta_VarHigh);
    TrackPosKaonCuts->SetPtRange(Kaon_pT_VarHigh,pThighK);
    TrackNegKaonCuts->SetPtRange(Kaon_pT_VarHigh,pThighK);
 
    config->SetDeltaEtaMax(DEta_VarLow);
 
  } else if (suffix == "5"){
 
    TrackCutsProton->SetEtaRange(-Proton_Eta_VarLow,Proton_Eta_VarLow);
    TrackCutsAntiProton->SetEtaRange(-Proton_Eta_VarLow,Proton_Eta_VarLow);
    TrackCutsProton->SetPtRange(Proton_pT_VarHigh,4.05);
    TrackCutsAntiProton->SetPtRange(Proton_pT_VarHigh,4.05);
    TrackCutsProton->SetNClsTPC(Proton_Clusters_VarHigh);
    TrackCutsAntiProton->SetNClsTPC(Proton_Clusters_VarHigh);
 
    TrackPosKaonCuts->SetPIDkd(true, UseRamonaCut,Kaon_SigmaCombined_VarLow, Kaon_SigmaTPC_VarHigh, Kaon_SigmaExclusion_VarLow);
    TrackNegKaonCuts->SetPIDkd(true, UseRamonaCut,Kaon_SigmaCombined_VarLow, Kaon_SigmaTPC_VarHigh, Kaon_SigmaExclusion_VarLow);
    TrackPosKaonCuts->SetPtRange(Kaon_pT_VarLow,pThighK);
    TrackNegKaonCuts->SetPtRange(Kaon_pT_VarLow,pThighK);
    TrackPosKaonCuts->SetNClsTPC(Kaon_Clusters_VarLow);
    TrackNegKaonCuts->SetNClsTPC(Kaon_Clusters_VarLow);
 
    config->SetDeltaPhiMax(DPhi_VarLow);
    config->SetDeltaEtaMax(DEta_VarLow);
 
  } else if (suffix == "6"){
 
    TrackCutsProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarLow);
    TrackCutsAntiProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarLow);
    TrackCutsProton->SetPtRange(Proton_pT_VarHigh,4.05);
    TrackCutsAntiProton->SetPtRange(Proton_pT_VarHigh,4.05);
    TrackCutsProton->SetNClsTPC(Proton_Clusters_VarLow);
    TrackCutsAntiProton->SetNClsTPC(Proton_Clusters_VarLow);
 
    TrackPosKaonCuts->SetPIDkd(true, UseRamonaCut,Kaon_SigmaCombined_VarLow, Kaon_SigmaTPC_VarHigh, Kaon_SigmaExclusion_VarLow);
    TrackNegKaonCuts->SetPIDkd(true, UseRamonaCut,Kaon_SigmaCombined_VarLow, Kaon_SigmaTPC_VarHigh, Kaon_SigmaExclusion_VarLow);
    TrackPosKaonCuts->SetPtRange(Kaon_pT_VarLow,pThighK);
    TrackNegKaonCuts->SetPtRange(Kaon_pT_VarLow,pThighK);
    TrackPosKaonCuts->SetNClsTPC(Kaon_Clusters_VarLow);
    TrackNegKaonCuts->SetNClsTPC(Kaon_Clusters_VarLow);
 
    config->SetDeltaEtaMax(DEta_VarLow);
 
  } else if (suffix == "7"){
 
    TrackCutsProton->SetPtRange(Proton_pT_VarHigh,4.05);
    TrackCutsAntiProton->SetPtRange(Proton_pT_VarHigh,4.05);
    TrackCutsProton->SetNClsTPC(Proton_Clusters_VarHigh);
    TrackCutsAntiProton->SetNClsTPC(Proton_Clusters_VarHigh);
 
    TrackPosKaonCuts->SetPIDkd(true, UseRamonaCut,Kaon_SigmaCombined_VarHigh, Kaon_SigmaTPC_Default, Kaon_SigmaExclusion_VarLow);
    TrackNegKaonCuts->SetPIDkd(true, UseRamonaCut,Kaon_SigmaCombined_VarHigh, Kaon_SigmaTPC_Default, Kaon_SigmaExclusion_VarLow);
    TrackPosKaonCuts->SetPtRange(Kaon_pT_VarLow,pThighK);
    TrackNegKaonCuts->SetPtRange(Kaon_pT_VarLow,pThighK);
    TrackPosKaonCuts->SetNClsTPC(Kaon_Clusters_VarHigh);
    TrackNegKaonCuts->SetNClsTPC(Kaon_Clusters_VarHigh);
 
    config->SetDeltaEtaMax(DEta_VarHigh);
 
  } else if (suffix == "8"){
 
    TrackCutsProton->SetEtaRange(-Proton_Eta_VarHigh,Proton_Eta_VarHigh);
    TrackCutsAntiProton->SetEtaRange(-Proton_Eta_VarHigh,Proton_Eta_VarHigh);
    TrackCutsProton->SetPtRange(Proton_pT_VarHigh,4.05);
    TrackCutsAntiProton->SetPtRange(Proton_pT_VarHigh,4.05);
    TrackCutsProton->SetNClsTPC(Proton_Clusters_VarLow);
    TrackCutsAntiProton->SetNClsTPC(Proton_Clusters_VarLow);
 
    TrackPosKaonCuts->SetPIDkd(true, UseRamonaCut,Kaon_SigmaCombined_Default, Kaon_SigmaTPC_VarLow, Kaon_SigmaExclusion_Default);
    TrackNegKaonCuts->SetPIDkd(true, UseRamonaCut,Kaon_SigmaCombined_Default, Kaon_SigmaTPC_VarLow, Kaon_SigmaExclusion_Default);
    TrackPosKaonCuts->SetEtaRange(-Kaon_Eta_VarHigh,Kaon_Eta_VarHigh);
    TrackNegKaonCuts->SetEtaRange(-Kaon_Eta_VarHigh,Kaon_Eta_VarHigh);
    TrackPosKaonCuts->SetPtRange(Kaon_pT_VarLow,pThighK);
    TrackNegKaonCuts->SetPtRange(Kaon_pT_VarLow,pThighK);
 
    config->SetDeltaPhiMax(DPhi_VarLow);
 
  } else if (suffix == "9"){
 
    TrackCutsProton->SetEtaRange(-Proton_Eta_VarHigh,Proton_Eta_VarHigh);
    TrackCutsAntiProton->SetEtaRange(-Proton_Eta_VarHigh,Proton_Eta_VarHigh);
    TrackCutsProton->SetPtRange(Proton_pT_VarHigh,4.05);
    TrackCutsAntiProton->SetPtRange(Proton_pT_VarHigh,4.05);
    TrackCutsProton->SetNClsTPC(Proton_Clusters_VarLow);
    TrackCutsAntiProton->SetNClsTPC(Proton_Clusters_VarLow);
 
    TrackPosKaonCuts->SetPIDkd(true, UseRamonaCut,Kaon_SigmaCombined_Default, Kaon_SigmaTPC_VarLow, Kaon_SigmaExclusion_Default);
    TrackNegKaonCuts->SetPIDkd(true, UseRamonaCut,Kaon_SigmaCombined_Default, Kaon_SigmaTPC_VarLow, Kaon_SigmaExclusion_Default);
    TrackPosKaonCuts->SetPtRange(Kaon_pT_VarHigh,pThighK);
    TrackNegKaonCuts->SetPtRange(Kaon_pT_VarHigh,pThighK);
    TrackPosKaonCuts->SetNClsTPC(Kaon_Clusters_VarHigh);
    TrackNegKaonCuts->SetNClsTPC(Kaon_Clusters_VarHigh);
 
    config->SetDeltaPhiMax(DPhi_VarLow);
    config->SetDeltaEtaMax(DEta_VarLow);
 
  } else if (suffix == "10"){
 
    TrackCutsProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarHigh);
    TrackCutsAntiProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarHigh);
    TrackCutsProton->SetEtaRange(-Proton_Eta_VarHigh,Proton_Eta_VarHigh);
    TrackCutsAntiProton->SetEtaRange(-Proton_Eta_VarHigh,Proton_Eta_VarHigh);
    TrackCutsProton->SetPtRange(Proton_pT_VarHigh,4.05);
    TrackCutsAntiProton->SetPtRange(Proton_pT_VarHigh,4.05);
 
    TrackPosKaonCuts->SetPIDkd(true, UseRamonaCut,Kaon_SigmaCombined_VarHigh, Kaon_SigmaTPC_VarHigh, Kaon_SigmaExclusion_VarLow);
    TrackNegKaonCuts->SetPIDkd(true, UseRamonaCut,Kaon_SigmaCombined_VarHigh, Kaon_SigmaTPC_VarHigh, Kaon_SigmaExclusion_VarLow);
    TrackPosKaonCuts->SetEtaRange(-Kaon_Eta_VarHigh,Kaon_Eta_VarHigh);
    TrackNegKaonCuts->SetEtaRange(-Kaon_Eta_VarHigh,Kaon_Eta_VarHigh);
 
    config->SetDeltaPhiMax(DPhi_VarLow);
    config->SetDeltaEtaMax(DEta_VarHigh);
 
  } else if (suffix == "11"){
 
    TrackCutsProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarHigh);
    TrackCutsAntiProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarHigh);
    TrackCutsProton->SetEtaRange(-Proton_Eta_VarHigh,Proton_Eta_VarHigh);
    TrackCutsAntiProton->SetEtaRange(-Proton_Eta_VarHigh,Proton_Eta_VarHigh);
    TrackCutsProton->SetPtRange(Proton_pT_VarLow,4.05);
    TrackCutsAntiProton->SetPtRange(Proton_pT_VarLow,4.05);
    TrackCutsProton->SetNClsTPC(Proton_Clusters_VarLow);
    TrackCutsAntiProton->SetNClsTPC(Proton_Clusters_VarLow);
 
    TrackPosKaonCuts->SetPIDkd(true, UseRamonaCut,Kaon_SigmaCombined_VarLow, Kaon_SigmaTPC_VarHigh, Kaon_SigmaExclusion_VarHigh);
    TrackNegKaonCuts->SetPIDkd(true, UseRamonaCut,Kaon_SigmaCombined_VarLow, Kaon_SigmaTPC_VarHigh, Kaon_SigmaExclusion_VarHigh);
 
 
  } else if (suffix == "12"){
 
    TrackCutsProton->SetEtaRange(-Proton_Eta_VarHigh,Proton_Eta_VarHigh);
    TrackCutsAntiProton->SetEtaRange(-Proton_Eta_VarHigh,Proton_Eta_VarHigh);
    TrackCutsProton->SetNClsTPC(Proton_Clusters_VarHigh);
    TrackCutsAntiProton->SetNClsTPC(Proton_Clusters_VarHigh);
 
    TrackPosKaonCuts->SetPIDkd(true, UseRamonaCut,Kaon_SigmaCombined_VarLow, Kaon_SigmaTPC_VarHigh, Kaon_SigmaExclusion_VarLow);
    TrackNegKaonCuts->SetPIDkd(true, UseRamonaCut,Kaon_SigmaCombined_VarLow, Kaon_SigmaTPC_VarHigh, Kaon_SigmaExclusion_VarLow);
    TrackPosKaonCuts->SetEtaRange(-Kaon_Eta_VarHigh,Kaon_Eta_VarHigh);
    TrackNegKaonCuts->SetEtaRange(-Kaon_Eta_VarHigh,Kaon_Eta_VarHigh);
    TrackPosKaonCuts->SetPtRange(Kaon_pT_VarHigh,pThighK);
    TrackNegKaonCuts->SetPtRange(Kaon_pT_VarHigh,pThighK);
    TrackPosKaonCuts->SetNClsTPC(Kaon_Clusters_VarHigh);
    TrackNegKaonCuts->SetNClsTPC(Kaon_Clusters_VarHigh);
 
    config->SetDeltaPhiMax(DPhi_VarHigh);
 
  } else if (suffix == "13"){
 
    TrackCutsProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarHigh);
    TrackCutsAntiProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarHigh);
    TrackCutsProton->SetPtRange(Proton_pT_VarHigh,4.05);
    TrackCutsAntiProton->SetPtRange(Proton_pT_VarHigh,4.05);
    TrackCutsProton->SetNClsTPC(Proton_Clusters_VarHigh);
    TrackCutsAntiProton->SetNClsTPC(Proton_Clusters_VarHigh);
 
    TrackPosKaonCuts->SetPIDkd(true, UseRamonaCut,Kaon_SigmaCombined_VarLow, Kaon_SigmaTPC_VarLow, Kaon_SigmaExclusion_VarHigh);
    TrackNegKaonCuts->SetPIDkd(true, UseRamonaCut,Kaon_SigmaCombined_VarLow, Kaon_SigmaTPC_VarLow, Kaon_SigmaExclusion_VarHigh);
    TrackPosKaonCuts->SetEtaRange(-Kaon_Eta_VarHigh,Kaon_Eta_VarHigh);
    TrackNegKaonCuts->SetEtaRange(-Kaon_Eta_VarHigh,Kaon_Eta_VarHigh);
    TrackPosKaonCuts->SetNClsTPC(Kaon_Clusters_VarLow);
    TrackNegKaonCuts->SetNClsTPC(Kaon_Clusters_VarLow);
 
    config->SetDeltaPhiMax(DPhi_VarHigh);
    config->SetDeltaEtaMax(DEta_VarLow);
 
  } else if (suffix == "14"){
 
    TrackCutsProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarLow);
    TrackCutsAntiProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarLow);
    TrackCutsProton->SetPtRange(Proton_pT_VarLow,4.05);
    TrackCutsAntiProton->SetPtRange(Proton_pT_VarLow,4.05);
    TrackCutsProton->SetNClsTPC(Proton_Clusters_VarHigh);
    TrackCutsAntiProton->SetNClsTPC(Proton_Clusters_VarHigh);
 
    TrackPosKaonCuts->SetPIDkd(true, UseRamonaCut,Kaon_SigmaCombined_VarLow, Kaon_SigmaTPC_VarLow, Kaon_SigmaExclusion_VarHigh);
    TrackNegKaonCuts->SetPIDkd(true, UseRamonaCut,Kaon_SigmaCombined_VarLow, Kaon_SigmaTPC_VarLow, Kaon_SigmaExclusion_VarHigh);
    TrackPosKaonCuts->SetEtaRange(-Kaon_Eta_VarHigh,Kaon_Eta_VarHigh);
    TrackNegKaonCuts->SetEtaRange(-Kaon_Eta_VarHigh,Kaon_Eta_VarHigh);
    TrackPosKaonCuts->SetNClsTPC(Kaon_Clusters_VarHigh);
    TrackNegKaonCuts->SetNClsTPC(Kaon_Clusters_VarHigh);
 
 
  } else if (suffix == "15"){
 
    TrackCutsProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarLow);
    TrackCutsAntiProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarLow);
    TrackCutsProton->SetEtaRange(-Proton_Eta_VarLow,Proton_Eta_VarLow);
    TrackCutsAntiProton->SetEtaRange(-Proton_Eta_VarLow,Proton_Eta_VarLow);
    TrackCutsProton->SetPtRange(Proton_pT_VarHigh,4.05);
    TrackCutsAntiProton->SetPtRange(Proton_pT_VarHigh,4.05);
    TrackCutsProton->SetNClsTPC(Proton_Clusters_VarHigh);
    TrackCutsAntiProton->SetNClsTPC(Proton_Clusters_VarHigh);
 
    TrackPosKaonCuts->SetPIDkd(true, UseRamonaCut,Kaon_SigmaCombined_VarLow, Kaon_SigmaTPC_Default, Kaon_SigmaExclusion_VarHigh);
    TrackNegKaonCuts->SetPIDkd(true, UseRamonaCut,Kaon_SigmaCombined_VarLow, Kaon_SigmaTPC_Default, Kaon_SigmaExclusion_VarHigh);
    TrackPosKaonCuts->SetEtaRange(-Kaon_Eta_VarLow,Kaon_Eta_VarLow);
    TrackNegKaonCuts->SetEtaRange(-Kaon_Eta_VarLow,Kaon_Eta_VarLow);
    TrackPosKaonCuts->SetPtRange(Kaon_pT_VarLow,pThighK);
    TrackNegKaonCuts->SetPtRange(Kaon_pT_VarLow,pThighK);
    TrackPosKaonCuts->SetNClsTPC(Kaon_Clusters_VarLow);
    TrackNegKaonCuts->SetNClsTPC(Kaon_Clusters_VarLow);
 
    config->SetDeltaPhiMax(DPhi_VarHigh);
 
  } else if (suffix == "16"){
 
    TrackCutsProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarLow);
    TrackCutsAntiProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarLow);
    TrackCutsProton->SetPtRange(Proton_pT_VarHigh,4.05);
    TrackCutsAntiProton->SetPtRange(Proton_pT_VarHigh,4.05);
    TrackCutsProton->SetNClsTPC(Proton_Clusters_VarHigh);
    TrackCutsAntiProton->SetNClsTPC(Proton_Clusters_VarHigh);
 
    TrackPosKaonCuts->SetPIDkd(true, UseRamonaCut,Kaon_SigmaCombined_VarHigh, Kaon_SigmaTPC_VarLow, Kaon_SigmaExclusion_VarLow);
    TrackNegKaonCuts->SetPIDkd(true, UseRamonaCut,Kaon_SigmaCombined_VarHigh, Kaon_SigmaTPC_VarLow, Kaon_SigmaExclusion_VarLow);
    TrackPosKaonCuts->SetPtRange(Kaon_pT_VarHigh,pThighK);
    TrackNegKaonCuts->SetPtRange(Kaon_pT_VarHigh,pThighK);
    TrackPosKaonCuts->SetNClsTPC(Kaon_Clusters_VarLow);
    TrackNegKaonCuts->SetNClsTPC(Kaon_Clusters_VarLow);
 
    config->SetDeltaPhiMax(DPhi_VarHigh);
    config->SetDeltaEtaMax(DEta_VarLow);
 
  } else if (suffix == "17"){
 
    TrackCutsProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarHigh);
    TrackCutsAntiProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarHigh);
    TrackCutsProton->SetPtRange(Proton_pT_VarLow,4.05);
    TrackCutsAntiProton->SetPtRange(Proton_pT_VarLow,4.05);
 
    TrackPosKaonCuts->SetPIDkd(true, UseRamonaCut,Kaon_SigmaCombined_VarHigh, Kaon_SigmaTPC_Default, Kaon_SigmaExclusion_VarLow);
    TrackNegKaonCuts->SetPIDkd(true, UseRamonaCut,Kaon_SigmaCombined_VarHigh, Kaon_SigmaTPC_Default, Kaon_SigmaExclusion_VarLow);
    TrackPosKaonCuts->SetPtRange(Kaon_pT_VarHigh,pThighK);
    TrackNegKaonCuts->SetPtRange(Kaon_pT_VarHigh,pThighK);
 
    config->SetDeltaEtaMax(DEta_VarLow);
 
  } else if (suffix == "18"){
 
    TrackCutsProton->SetEtaRange(-Proton_Eta_VarHigh,Proton_Eta_VarHigh);
    TrackCutsAntiProton->SetEtaRange(-Proton_Eta_VarHigh,Proton_Eta_VarHigh);
    TrackCutsProton->SetPtRange(Proton_pT_VarLow,4.05);
    TrackCutsAntiProton->SetPtRange(Proton_pT_VarLow,4.05);
    TrackCutsProton->SetNClsTPC(Proton_Clusters_VarHigh);
    TrackCutsAntiProton->SetNClsTPC(Proton_Clusters_VarHigh);
 
    TrackPosKaonCuts->SetPIDkd(true, UseRamonaCut,Kaon_SigmaCombined_VarLow, Kaon_SigmaTPC_Default, Kaon_SigmaExclusion_VarHigh);
    TrackNegKaonCuts->SetPIDkd(true, UseRamonaCut,Kaon_SigmaCombined_VarLow, Kaon_SigmaTPC_Default, Kaon_SigmaExclusion_VarHigh);
    TrackPosKaonCuts->SetNClsTPC(Kaon_Clusters_VarHigh);
    TrackNegKaonCuts->SetNClsTPC(Kaon_Clusters_VarHigh);
 
    config->SetDeltaPhiMax(DPhi_VarLow);
    config->SetDeltaEtaMax(DEta_VarHigh);
 
  } else if (suffix == "19"){
 
    TrackCutsProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarLow);
    TrackCutsAntiProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarLow);
    TrackCutsProton->SetEtaRange(-Proton_Eta_VarLow,Proton_Eta_VarLow);
    TrackCutsAntiProton->SetEtaRange(-Proton_Eta_VarLow,Proton_Eta_VarLow);
    TrackCutsProton->SetPtRange(Proton_pT_VarLow,4.05);
    TrackCutsAntiProton->SetPtRange(Proton_pT_VarLow,4.05);
    TrackCutsProton->SetNClsTPC(Proton_Clusters_VarLow);
    TrackCutsAntiProton->SetNClsTPC(Proton_Clusters_VarLow);
 
    TrackPosKaonCuts->SetPIDkd(true, UseRamonaCut,Kaon_SigmaCombined_VarHigh, Kaon_SigmaTPC_VarLow, Kaon_SigmaExclusion_Default);
    TrackNegKaonCuts->SetPIDkd(true, UseRamonaCut,Kaon_SigmaCombined_VarHigh, Kaon_SigmaTPC_VarLow, Kaon_SigmaExclusion_Default);
    TrackPosKaonCuts->SetEtaRange(-Kaon_Eta_VarLow,Kaon_Eta_VarLow);
    TrackNegKaonCuts->SetEtaRange(-Kaon_Eta_VarLow,Kaon_Eta_VarLow);
    TrackPosKaonCuts->SetPtRange(Kaon_pT_VarHigh,pThighK);
    TrackNegKaonCuts->SetPtRange(Kaon_pT_VarHigh,pThighK);
 
    config->SetDeltaEtaMax(DEta_VarLow);
 
  } else if (suffix == "20"){
 
    TrackCutsProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarLow);
    TrackCutsAntiProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarLow);
    TrackCutsProton->SetEtaRange(-Proton_Eta_VarHigh,Proton_Eta_VarHigh);
    TrackCutsAntiProton->SetEtaRange(-Proton_Eta_VarHigh,Proton_Eta_VarHigh);
    TrackCutsProton->SetPtRange(Proton_pT_VarLow,4.05);
    TrackCutsAntiProton->SetPtRange(Proton_pT_VarLow,4.05);
 
    TrackPosKaonCuts->SetPIDkd(true, UseRamonaCut,Kaon_SigmaCombined_VarHigh, Kaon_SigmaTPC_VarLow, Kaon_SigmaExclusion_VarHigh);
    TrackNegKaonCuts->SetPIDkd(true, UseRamonaCut,Kaon_SigmaCombined_VarHigh, Kaon_SigmaTPC_VarLow, Kaon_SigmaExclusion_VarHigh);
    TrackPosKaonCuts->SetPtRange(Kaon_pT_VarLow,pThighK);
    TrackNegKaonCuts->SetPtRange(Kaon_pT_VarLow,pThighK);
 
    config->SetDeltaEtaMax(DEta_VarLow);
 
  } else if (suffix == "21"){
 
    TrackCutsProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarLow);
    TrackCutsAntiProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarLow);
    TrackCutsProton->SetEtaRange(-Proton_Eta_VarHigh,Proton_Eta_VarHigh);
    TrackCutsAntiProton->SetEtaRange(-Proton_Eta_VarHigh,Proton_Eta_VarHigh);
    TrackCutsProton->SetNClsTPC(Proton_Clusters_VarLow);
    TrackCutsAntiProton->SetNClsTPC(Proton_Clusters_VarLow);
 
    TrackPosKaonCuts->SetPIDkd(true, UseRamonaCut,Kaon_SigmaCombined_VarHigh, Kaon_SigmaTPC_VarHigh, Kaon_SigmaExclusion_VarLow);
    TrackNegKaonCuts->SetPIDkd(true, UseRamonaCut,Kaon_SigmaCombined_VarHigh, Kaon_SigmaTPC_VarHigh, Kaon_SigmaExclusion_VarLow);
    TrackPosKaonCuts->SetEtaRange(-Kaon_Eta_VarHigh,Kaon_Eta_VarHigh);
    TrackNegKaonCuts->SetEtaRange(-Kaon_Eta_VarHigh,Kaon_Eta_VarHigh);
    TrackPosKaonCuts->SetPtRange(Kaon_pT_VarHigh,pThighK);
    TrackNegKaonCuts->SetPtRange(Kaon_pT_VarHigh,pThighK);
 
    config->SetDeltaPhiMax(DPhi_VarLow);
 
  } else if (suffix == "22"){
 
    TrackCutsProton->SetEtaRange(-Proton_Eta_VarHigh,Proton_Eta_VarHigh);
    TrackCutsAntiProton->SetEtaRange(-Proton_Eta_VarHigh,Proton_Eta_VarHigh);
    TrackCutsProton->SetPtRange(Proton_pT_VarLow,4.05);
    TrackCutsAntiProton->SetPtRange(Proton_pT_VarLow,4.05);
    TrackCutsProton->SetNClsTPC(Proton_Clusters_VarLow);
    TrackCutsAntiProton->SetNClsTPC(Proton_Clusters_VarLow);
 
    TrackPosKaonCuts->SetPIDkd(true, UseRamonaCut,Kaon_SigmaCombined_VarLow, Kaon_SigmaTPC_Default, Kaon_SigmaExclusion_VarLow);
    TrackNegKaonCuts->SetPIDkd(true, UseRamonaCut,Kaon_SigmaCombined_VarLow, Kaon_SigmaTPC_Default, Kaon_SigmaExclusion_VarLow);
    TrackPosKaonCuts->SetEtaRange(-Kaon_Eta_VarHigh,Kaon_Eta_VarHigh);
    TrackNegKaonCuts->SetEtaRange(-Kaon_Eta_VarHigh,Kaon_Eta_VarHigh);
    TrackPosKaonCuts->SetPtRange(Kaon_pT_VarLow,pThighK);
    TrackNegKaonCuts->SetPtRange(Kaon_pT_VarLow,pThighK);
    TrackPosKaonCuts->SetNClsTPC(Kaon_Clusters_VarLow);
    TrackNegKaonCuts->SetNClsTPC(Kaon_Clusters_VarLow);
 
    config->SetDeltaPhiMax(DPhi_VarHigh);
    config->SetDeltaEtaMax(DEta_VarLow);
 
  } else if (suffix == "23"){
 
    TrackCutsProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarLow);
    TrackCutsAntiProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarLow);
    TrackCutsProton->SetEtaRange(-Proton_Eta_VarLow,Proton_Eta_VarLow);
    TrackCutsAntiProton->SetEtaRange(-Proton_Eta_VarLow,Proton_Eta_VarLow);
    TrackCutsProton->SetPtRange(Proton_pT_VarHigh,4.05);
    TrackCutsAntiProton->SetPtRange(Proton_pT_VarHigh,4.05);
 
    TrackPosKaonCuts->SetPIDkd(true, UseRamonaCut,Kaon_SigmaCombined_VarLow, Kaon_SigmaTPC_VarLow, Kaon_SigmaExclusion_VarLow);
    TrackNegKaonCuts->SetPIDkd(true, UseRamonaCut,Kaon_SigmaCombined_VarLow, Kaon_SigmaTPC_VarLow, Kaon_SigmaExclusion_VarLow);
    TrackPosKaonCuts->SetEtaRange(-Kaon_Eta_VarHigh,Kaon_Eta_VarHigh);
    TrackNegKaonCuts->SetEtaRange(-Kaon_Eta_VarHigh,Kaon_Eta_VarHigh);
    TrackPosKaonCuts->SetPtRange(Kaon_pT_VarHigh,pThighK);
    TrackNegKaonCuts->SetPtRange(Kaon_pT_VarHigh,pThighK);
    TrackPosKaonCuts->SetNClsTPC(Kaon_Clusters_VarHigh);
    TrackNegKaonCuts->SetNClsTPC(Kaon_Clusters_VarHigh);
 
    config->SetDeltaEtaMax(DEta_VarHigh);
 
  } else if (suffix == "24"){
 
    TrackCutsProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarLow);
    TrackCutsAntiProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarLow);
    TrackCutsProton->SetEtaRange(-Proton_Eta_VarHigh,Proton_Eta_VarHigh);
    TrackCutsAntiProton->SetEtaRange(-Proton_Eta_VarHigh,Proton_Eta_VarHigh);
    TrackCutsProton->SetPtRange(Proton_pT_VarLow,4.05);
    TrackCutsAntiProton->SetPtRange(Proton_pT_VarLow,4.05);
 
    TrackPosKaonCuts->SetPIDkd(true, UseRamonaCut,Kaon_SigmaCombined_VarHigh, Kaon_SigmaTPC_VarLow, Kaon_SigmaExclusion_VarHigh);
    TrackNegKaonCuts->SetPIDkd(true, UseRamonaCut,Kaon_SigmaCombined_VarHigh, Kaon_SigmaTPC_VarLow, Kaon_SigmaExclusion_VarHigh);
    TrackPosKaonCuts->SetPtRange(Kaon_pT_VarLow,pThighK);
    TrackNegKaonCuts->SetPtRange(Kaon_pT_VarLow,pThighK);
 
    config->SetDeltaPhiMax(DPhi_VarLow);
    config->SetDeltaEtaMax(DEta_VarLow);
 
  } else if (suffix == "25"){
 
    TrackCutsProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarHigh);
    TrackCutsAntiProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarHigh);
    TrackCutsProton->SetEtaRange(-Proton_Eta_VarHigh,Proton_Eta_VarHigh);
    TrackCutsAntiProton->SetEtaRange(-Proton_Eta_VarHigh,Proton_Eta_VarHigh);
    TrackCutsProton->SetPtRange(Proton_pT_VarHigh,4.05);
    TrackCutsAntiProton->SetPtRange(Proton_pT_VarHigh,4.05);
    TrackCutsProton->SetNClsTPC(Proton_Clusters_VarHigh);
    TrackCutsAntiProton->SetNClsTPC(Proton_Clusters_VarHigh);
 
    TrackPosKaonCuts->SetPIDkd(true, UseRamonaCut,Kaon_SigmaCombined_Default, Kaon_SigmaTPC_Default, Kaon_SigmaExclusion_VarLow);
    TrackNegKaonCuts->SetPIDkd(true, UseRamonaCut,Kaon_SigmaCombined_Default, Kaon_SigmaTPC_Default, Kaon_SigmaExclusion_VarLow);
    TrackPosKaonCuts->SetEtaRange(-Kaon_Eta_VarLow,Kaon_Eta_VarLow);
    TrackNegKaonCuts->SetEtaRange(-Kaon_Eta_VarLow,Kaon_Eta_VarLow);
    TrackPosKaonCuts->SetPtRange(Kaon_pT_VarHigh,pThighK);
    TrackNegKaonCuts->SetPtRange(Kaon_pT_VarHigh,pThighK);
    TrackPosKaonCuts->SetNClsTPC(Kaon_Clusters_VarHigh);
    TrackNegKaonCuts->SetNClsTPC(Kaon_Clusters_VarHigh);
 
 
  } else if (suffix == "26"){
 
    TrackCutsProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarHigh);
    TrackCutsAntiProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarHigh);
    TrackCutsProton->SetPtRange(Proton_pT_VarHigh,4.05);
    TrackCutsAntiProton->SetPtRange(Proton_pT_VarHigh,4.05);
 
    TrackPosKaonCuts->SetPIDkd(true, UseRamonaCut,Kaon_SigmaCombined_Default, Kaon_SigmaTPC_Default, Kaon_SigmaExclusion_VarHigh);
    TrackNegKaonCuts->SetPIDkd(true, UseRamonaCut,Kaon_SigmaCombined_Default, Kaon_SigmaTPC_Default, Kaon_SigmaExclusion_VarHigh);
    TrackPosKaonCuts->SetEtaRange(-Kaon_Eta_VarHigh,Kaon_Eta_VarHigh);
    TrackNegKaonCuts->SetEtaRange(-Kaon_Eta_VarHigh,Kaon_Eta_VarHigh);
    TrackPosKaonCuts->SetPtRange(Kaon_pT_VarHigh,pThighK);
    TrackNegKaonCuts->SetPtRange(Kaon_pT_VarHigh,pThighK);
    TrackPosKaonCuts->SetNClsTPC(Kaon_Clusters_VarHigh);
    TrackNegKaonCuts->SetNClsTPC(Kaon_Clusters_VarHigh);
 
    config->SetDeltaPhiMax(DPhi_VarLow);
    config->SetDeltaEtaMax(DEta_VarLow);
 
  } else if (suffix == "27"){
 
    TrackCutsProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarLow);
    TrackCutsAntiProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarLow);
    TrackCutsProton->SetEtaRange(-Proton_Eta_VarLow,Proton_Eta_VarLow);
    TrackCutsAntiProton->SetEtaRange(-Proton_Eta_VarLow,Proton_Eta_VarLow);
    TrackCutsProton->SetPtRange(Proton_pT_VarLow,4.05);
    TrackCutsAntiProton->SetPtRange(Proton_pT_VarLow,4.05);
    TrackCutsProton->SetNClsTPC(Proton_Clusters_VarHigh);
    TrackCutsAntiProton->SetNClsTPC(Proton_Clusters_VarHigh);
 
    TrackPosKaonCuts->SetPIDkd(true, UseRamonaCut,Kaon_SigmaCombined_Default, Kaon_SigmaTPC_VarHigh, Kaon_SigmaExclusion_VarLow);
    TrackNegKaonCuts->SetPIDkd(true, UseRamonaCut,Kaon_SigmaCombined_Default, Kaon_SigmaTPC_VarHigh, Kaon_SigmaExclusion_VarLow);
    TrackPosKaonCuts->SetEtaRange(-Kaon_Eta_VarHigh,Kaon_Eta_VarHigh);
    TrackNegKaonCuts->SetEtaRange(-Kaon_Eta_VarHigh,Kaon_Eta_VarHigh);
    TrackPosKaonCuts->SetPtRange(Kaon_pT_VarLow,pThighK);
    TrackNegKaonCuts->SetPtRange(Kaon_pT_VarLow,pThighK);
    TrackPosKaonCuts->SetNClsTPC(Kaon_Clusters_VarHigh);
    TrackNegKaonCuts->SetNClsTPC(Kaon_Clusters_VarHigh);
 
    config->SetDeltaPhiMax(DPhi_VarLow);
 
  } else if (suffix == "28"){
 
    TrackCutsProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarHigh);
    TrackCutsAntiProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarHigh);
    TrackCutsProton->SetEtaRange(-Proton_Eta_VarHigh,Proton_Eta_VarHigh);
    TrackCutsAntiProton->SetEtaRange(-Proton_Eta_VarHigh,Proton_Eta_VarHigh);
    TrackCutsProton->SetPtRange(Proton_pT_VarHigh,4.05);
    TrackCutsAntiProton->SetPtRange(Proton_pT_VarHigh,4.05);
    TrackCutsProton->SetNClsTPC(Proton_Clusters_VarHigh);
    TrackCutsAntiProton->SetNClsTPC(Proton_Clusters_VarHigh);
 
    TrackPosKaonCuts->SetPIDkd(true, UseRamonaCut,Kaon_SigmaCombined_Default, Kaon_SigmaTPC_Default, Kaon_SigmaExclusion_VarHigh);
    TrackNegKaonCuts->SetPIDkd(true, UseRamonaCut,Kaon_SigmaCombined_Default, Kaon_SigmaTPC_Default, Kaon_SigmaExclusion_VarHigh);
    TrackPosKaonCuts->SetEtaRange(-Kaon_Eta_VarHigh,Kaon_Eta_VarHigh);
    TrackNegKaonCuts->SetEtaRange(-Kaon_Eta_VarHigh,Kaon_Eta_VarHigh);
    TrackPosKaonCuts->SetPtRange(Kaon_pT_VarHigh,pThighK);
    TrackNegKaonCuts->SetPtRange(Kaon_pT_VarHigh,pThighK);
    TrackPosKaonCuts->SetNClsTPC(Kaon_Clusters_VarLow);
    TrackNegKaonCuts->SetNClsTPC(Kaon_Clusters_VarLow);
 
    config->SetDeltaPhiMax(DPhi_VarHigh);
    config->SetDeltaEtaMax(DEta_VarLow);
 
  } else if (suffix == "29"){
 
    TrackCutsProton->SetEtaRange(-Proton_Eta_VarHigh,Proton_Eta_VarHigh);
    TrackCutsAntiProton->SetEtaRange(-Proton_Eta_VarHigh,Proton_Eta_VarHigh);
    TrackCutsProton->SetPtRange(Proton_pT_VarLow,4.05);
    TrackCutsAntiProton->SetPtRange(Proton_pT_VarLow,4.05);
    TrackCutsProton->SetNClsTPC(Proton_Clusters_VarHigh);
    TrackCutsAntiProton->SetNClsTPC(Proton_Clusters_VarHigh);
 
    TrackPosKaonCuts->SetPIDkd(true, UseRamonaCut,Kaon_SigmaCombined_Default, Kaon_SigmaTPC_VarLow, Kaon_SigmaExclusion_Default);
    TrackNegKaonCuts->SetPIDkd(true, UseRamonaCut,Kaon_SigmaCombined_Default, Kaon_SigmaTPC_VarLow, Kaon_SigmaExclusion_Default);
    TrackPosKaonCuts->SetEtaRange(-Kaon_Eta_VarLow,Kaon_Eta_VarLow);
    TrackNegKaonCuts->SetEtaRange(-Kaon_Eta_VarLow,Kaon_Eta_VarLow);
    TrackPosKaonCuts->SetPtRange(Kaon_pT_VarLow,pThighK);
    TrackNegKaonCuts->SetPtRange(Kaon_pT_VarLow,pThighK);
    TrackPosKaonCuts->SetNClsTPC(Kaon_Clusters_VarLow);
    TrackNegKaonCuts->SetNClsTPC(Kaon_Clusters_VarLow);
 
    config->SetDeltaPhiMax(DPhi_VarLow);
    config->SetDeltaEtaMax(DEta_VarLow);
 
  } else if (suffix == "30"){
 
    TrackCutsProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarLow);
    TrackCutsAntiProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarLow);
    TrackCutsProton->SetPtRange(Proton_pT_VarLow,4.05);
    TrackCutsAntiProton->SetPtRange(Proton_pT_VarLow,4.05);
    TrackCutsProton->SetNClsTPC(Proton_Clusters_VarLow);
    TrackCutsAntiProton->SetNClsTPC(Proton_Clusters_VarLow);
 
    TrackPosKaonCuts->SetPIDkd(true, UseRamonaCut,Kaon_SigmaCombined_VarHigh, Kaon_SigmaTPC_Default, Kaon_SigmaExclusion_VarHigh);
    TrackNegKaonCuts->SetPIDkd(true, UseRamonaCut,Kaon_SigmaCombined_VarHigh, Kaon_SigmaTPC_Default, Kaon_SigmaExclusion_VarHigh);
    TrackPosKaonCuts->SetEtaRange(-Kaon_Eta_VarHigh,Kaon_Eta_VarHigh);
    TrackNegKaonCuts->SetEtaRange(-Kaon_Eta_VarHigh,Kaon_Eta_VarHigh);
    TrackPosKaonCuts->SetPtRange(Kaon_pT_VarHigh,pThighK);
    TrackNegKaonCuts->SetPtRange(Kaon_pT_VarHigh,pThighK);
    TrackPosKaonCuts->SetNClsTPC(Kaon_Clusters_VarHigh);
    TrackNegKaonCuts->SetNClsTPC(Kaon_Clusters_VarHigh);
 
    config->SetDeltaPhiMax(DPhi_VarLow);
    config->SetDeltaEtaMax(DEta_VarHigh);
 
  } else if (suffix == "31"){
 
    TrackCutsProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarHigh);
    TrackCutsAntiProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarHigh);
 
    TrackPosKaonCuts->SetPIDkd(true, UseRamonaCut,Kaon_SigmaCombined_VarLow, Kaon_SigmaTPC_Default, Kaon_SigmaExclusion_Default);
    TrackNegKaonCuts->SetPIDkd(true, UseRamonaCut,Kaon_SigmaCombined_VarLow, Kaon_SigmaTPC_Default, Kaon_SigmaExclusion_Default);
    TrackPosKaonCuts->SetEtaRange(-Kaon_Eta_VarHigh,Kaon_Eta_VarHigh);
    TrackNegKaonCuts->SetEtaRange(-Kaon_Eta_VarHigh,Kaon_Eta_VarHigh);
    TrackPosKaonCuts->SetNClsTPC(Kaon_Clusters_VarHigh);
    TrackNegKaonCuts->SetNClsTPC(Kaon_Clusters_VarHigh);
 
    config->SetDeltaPhiMax(DPhi_VarLow);
    config->SetDeltaEtaMax(DEta_VarHigh);
 
  } else if (suffix == "32"){
 
    TrackCutsProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarLow);
    TrackCutsAntiProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarLow);
    TrackCutsProton->SetPtRange(Proton_pT_VarLow,4.05);
    TrackCutsAntiProton->SetPtRange(Proton_pT_VarLow,4.05);
    TrackCutsProton->SetNClsTPC(Proton_Clusters_VarLow);
    TrackCutsAntiProton->SetNClsTPC(Proton_Clusters_VarLow);
 
    TrackPosKaonCuts->SetPIDkd(true, UseRamonaCut,Kaon_SigmaCombined_VarHigh, Kaon_SigmaTPC_Default, Kaon_SigmaExclusion_VarHigh);
    TrackNegKaonCuts->SetPIDkd(true, UseRamonaCut,Kaon_SigmaCombined_VarHigh, Kaon_SigmaTPC_Default, Kaon_SigmaExclusion_VarHigh);
 
    config->SetDeltaPhiMax(DPhi_VarHigh);
    config->SetDeltaEtaMax(DEta_VarHigh);
 
  } else if (suffix == "33"){
 
    TrackCutsProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarHigh);
    TrackCutsAntiProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarHigh);
    TrackCutsProton->SetPtRange(Proton_pT_VarHigh,4.05);
    TrackCutsAntiProton->SetPtRange(Proton_pT_VarHigh,4.05);
    TrackCutsProton->SetNClsTPC(Proton_Clusters_VarHigh);
    TrackCutsAntiProton->SetNClsTPC(Proton_Clusters_VarHigh);
 
    TrackPosKaonCuts->SetPIDkd(true, UseRamonaCut,Kaon_SigmaCombined_VarLow, Kaon_SigmaTPC_VarHigh, Kaon_SigmaExclusion_VarLow);
    TrackNegKaonCuts->SetPIDkd(true, UseRamonaCut,Kaon_SigmaCombined_VarLow, Kaon_SigmaTPC_VarHigh, Kaon_SigmaExclusion_VarLow);
    TrackPosKaonCuts->SetEtaRange(-Kaon_Eta_VarLow,Kaon_Eta_VarLow);
    TrackNegKaonCuts->SetEtaRange(-Kaon_Eta_VarLow,Kaon_Eta_VarLow);
    TrackPosKaonCuts->SetNClsTPC(Kaon_Clusters_VarHigh);
    TrackNegKaonCuts->SetNClsTPC(Kaon_Clusters_VarHigh);
 
    config->SetDeltaPhiMax(DPhi_VarHigh);
 
  } else if (suffix == "34"){
 
    TrackCutsProton->SetEtaRange(-Proton_Eta_VarHigh,Proton_Eta_VarHigh);
    TrackCutsAntiProton->SetEtaRange(-Proton_Eta_VarHigh,Proton_Eta_VarHigh);
    TrackCutsProton->SetNClsTPC(Proton_Clusters_VarHigh);
    TrackCutsAntiProton->SetNClsTPC(Proton_Clusters_VarHigh);
 
    TrackPosKaonCuts->SetPIDkd(true, UseRamonaCut,Kaon_SigmaCombined_Default, Kaon_SigmaTPC_VarLow, Kaon_SigmaExclusion_VarLow);
    TrackNegKaonCuts->SetPIDkd(true, UseRamonaCut,Kaon_SigmaCombined_Default, Kaon_SigmaTPC_VarLow, Kaon_SigmaExclusion_VarLow);
    TrackPosKaonCuts->SetEtaRange(-Kaon_Eta_VarHigh,Kaon_Eta_VarHigh);
    TrackNegKaonCuts->SetEtaRange(-Kaon_Eta_VarHigh,Kaon_Eta_VarHigh);
    TrackPosKaonCuts->SetNClsTPC(Kaon_Clusters_VarLow);
    TrackNegKaonCuts->SetNClsTPC(Kaon_Clusters_VarLow);
 
    config->SetDeltaEtaMax(DEta_VarLow);
 
  } else if (suffix == "35"){
 
    TrackCutsProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarHigh);
    TrackCutsAntiProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarHigh);
    TrackCutsProton->SetEtaRange(-Proton_Eta_VarLow,Proton_Eta_VarLow);
    TrackCutsAntiProton->SetEtaRange(-Proton_Eta_VarLow,Proton_Eta_VarLow);
    TrackCutsProton->SetPtRange(Proton_pT_VarLow,4.05);
    TrackCutsAntiProton->SetPtRange(Proton_pT_VarLow,4.05);
    TrackCutsProton->SetNClsTPC(Proton_Clusters_VarHigh);
    TrackCutsAntiProton->SetNClsTPC(Proton_Clusters_VarHigh);
 
    TrackPosKaonCuts->SetEtaRange(-Kaon_Eta_VarLow,Kaon_Eta_VarLow);
    TrackNegKaonCuts->SetEtaRange(-Kaon_Eta_VarLow,Kaon_Eta_VarLow);
    TrackPosKaonCuts->SetPtRange(Kaon_pT_VarHigh,pThighK);
    TrackNegKaonCuts->SetPtRange(Kaon_pT_VarHigh,pThighK);
    TrackPosKaonCuts->SetNClsTPC(Kaon_Clusters_VarHigh);
    TrackNegKaonCuts->SetNClsTPC(Kaon_Clusters_VarHigh);
 
    config->SetDeltaEtaMax(DEta_VarLow);
 
  } else if (suffix == "36"){
 
    TrackCutsProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarLow);
    TrackCutsAntiProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarLow);
    TrackCutsProton->SetEtaRange(-Proton_Eta_VarLow,Proton_Eta_VarLow);
    TrackCutsAntiProton->SetEtaRange(-Proton_Eta_VarLow,Proton_Eta_VarLow);
    TrackCutsProton->SetPtRange(Proton_pT_VarLow,4.05);
    TrackCutsAntiProton->SetPtRange(Proton_pT_VarLow,4.05);
    TrackCutsProton->SetNClsTPC(Proton_Clusters_VarHigh);
    TrackCutsAntiProton->SetNClsTPC(Proton_Clusters_VarHigh);
 
    TrackPosKaonCuts->SetPIDkd(true, UseRamonaCut,Kaon_SigmaCombined_VarLow, Kaon_SigmaTPC_VarLow, Kaon_SigmaExclusion_VarHigh);
    TrackNegKaonCuts->SetPIDkd(true, UseRamonaCut,Kaon_SigmaCombined_VarLow, Kaon_SigmaTPC_VarLow, Kaon_SigmaExclusion_VarHigh);
    TrackPosKaonCuts->SetEtaRange(-Kaon_Eta_VarHigh,Kaon_Eta_VarHigh);
    TrackNegKaonCuts->SetEtaRange(-Kaon_Eta_VarHigh,Kaon_Eta_VarHigh);
    TrackPosKaonCuts->SetPtRange(Kaon_pT_VarLow,pThighK);
    TrackNegKaonCuts->SetPtRange(Kaon_pT_VarLow,pThighK);
    TrackPosKaonCuts->SetNClsTPC(Kaon_Clusters_VarLow);
    TrackNegKaonCuts->SetNClsTPC(Kaon_Clusters_VarLow);
 
    config->SetDeltaEtaMax(DEta_VarLow);
 
  } else if (suffix == "37"){
 
    TrackCutsProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarLow);
    TrackCutsAntiProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarLow);
    TrackCutsProton->SetEtaRange(-Proton_Eta_VarHigh,Proton_Eta_VarHigh);
    TrackCutsAntiProton->SetEtaRange(-Proton_Eta_VarHigh,Proton_Eta_VarHigh);
    TrackCutsProton->SetPtRange(Proton_pT_VarLow,4.05);
    TrackCutsAntiProton->SetPtRange(Proton_pT_VarLow,4.05);
    TrackCutsProton->SetNClsTPC(Proton_Clusters_VarHigh);
    TrackCutsAntiProton->SetNClsTPC(Proton_Clusters_VarHigh);
 
    TrackPosKaonCuts->SetPIDkd(true, UseRamonaCut,Kaon_SigmaCombined_VarHigh, Kaon_SigmaTPC_VarHigh, Kaon_SigmaExclusion_VarHigh);
    TrackNegKaonCuts->SetPIDkd(true, UseRamonaCut,Kaon_SigmaCombined_VarHigh, Kaon_SigmaTPC_VarHigh, Kaon_SigmaExclusion_VarHigh);
    TrackPosKaonCuts->SetEtaRange(-Kaon_Eta_VarLow,Kaon_Eta_VarLow);
    TrackNegKaonCuts->SetEtaRange(-Kaon_Eta_VarLow,Kaon_Eta_VarLow);
    TrackPosKaonCuts->SetNClsTPC(Kaon_Clusters_VarLow);
    TrackNegKaonCuts->SetNClsTPC(Kaon_Clusters_VarLow);
 
    config->SetDeltaPhiMax(DPhi_VarLow);
 
  } else if (suffix == "38"){
 
    TrackCutsProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarHigh);
    TrackCutsAntiProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarHigh);
    TrackCutsProton->SetEtaRange(-Proton_Eta_VarHigh,Proton_Eta_VarHigh);
    TrackCutsAntiProton->SetEtaRange(-Proton_Eta_VarHigh,Proton_Eta_VarHigh);
 
    TrackPosKaonCuts->SetPIDkd(true, UseRamonaCut,Kaon_SigmaCombined_Default, Kaon_SigmaTPC_VarHigh, Kaon_SigmaExclusion_VarHigh);
    TrackNegKaonCuts->SetPIDkd(true, UseRamonaCut,Kaon_SigmaCombined_Default, Kaon_SigmaTPC_VarHigh, Kaon_SigmaExclusion_VarHigh);
    TrackPosKaonCuts->SetEtaRange(-Kaon_Eta_VarLow,Kaon_Eta_VarLow);
    TrackNegKaonCuts->SetEtaRange(-Kaon_Eta_VarLow,Kaon_Eta_VarLow);
    TrackPosKaonCuts->SetPtRange(Kaon_pT_VarHigh,pThighK);
    TrackNegKaonCuts->SetPtRange(Kaon_pT_VarHigh,pThighK);
    TrackPosKaonCuts->SetNClsTPC(Kaon_Clusters_VarHigh);
    TrackNegKaonCuts->SetNClsTPC(Kaon_Clusters_VarHigh);
 
    config->SetDeltaEtaMax(DEta_VarLow);
 
  } else if (suffix == "39"){
 
    TrackCutsProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarLow);
    TrackCutsAntiProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarLow);
    TrackCutsProton->SetEtaRange(-Proton_Eta_VarHigh,Proton_Eta_VarHigh);
    TrackCutsAntiProton->SetEtaRange(-Proton_Eta_VarHigh,Proton_Eta_VarHigh);
    TrackCutsProton->SetPtRange(Proton_pT_VarLow,4.05);
    TrackCutsAntiProton->SetPtRange(Proton_pT_VarLow,4.05);
 
    TrackPosKaonCuts->SetPIDkd(true, UseRamonaCut,Kaon_SigmaCombined_VarLow, Kaon_SigmaTPC_VarHigh, Kaon_SigmaExclusion_VarLow);
    TrackNegKaonCuts->SetPIDkd(true, UseRamonaCut,Kaon_SigmaCombined_VarLow, Kaon_SigmaTPC_VarHigh, Kaon_SigmaExclusion_VarLow);
    TrackPosKaonCuts->SetNClsTPC(Kaon_Clusters_VarLow);
    TrackNegKaonCuts->SetNClsTPC(Kaon_Clusters_VarLow);
 
    config->SetDeltaPhiMax(DPhi_VarLow);
    config->SetDeltaEtaMax(DEta_VarHigh);
 
  } else if (suffix == "40"){
 
    TrackCutsProton->SetEtaRange(-Proton_Eta_VarLow,Proton_Eta_VarLow);
    TrackCutsAntiProton->SetEtaRange(-Proton_Eta_VarLow,Proton_Eta_VarLow);
    TrackCutsProton->SetPtRange(Proton_pT_VarLow,4.05);
    TrackCutsAntiProton->SetPtRange(Proton_pT_VarLow,4.05);
 
    TrackPosKaonCuts->SetPIDkd(true, UseRamonaCut,Kaon_SigmaCombined_VarLow, Kaon_SigmaTPC_VarHigh, Kaon_SigmaExclusion_Default);
    TrackNegKaonCuts->SetPIDkd(true, UseRamonaCut,Kaon_SigmaCombined_VarLow, Kaon_SigmaTPC_VarHigh, Kaon_SigmaExclusion_Default);
    TrackPosKaonCuts->SetNClsTPC(Kaon_Clusters_VarHigh);
    TrackNegKaonCuts->SetNClsTPC(Kaon_Clusters_VarHigh);
 
    config->SetDeltaPhiMax(DPhi_VarLow);
    config->SetDeltaEtaMax(DEta_VarLow);
 
  } else if (suffix == "41"){
 
    TrackCutsProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarHigh);
    TrackCutsAntiProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarHigh);
    TrackCutsProton->SetEtaRange(-Proton_Eta_VarHigh,Proton_Eta_VarHigh);
    TrackCutsAntiProton->SetEtaRange(-Proton_Eta_VarHigh,Proton_Eta_VarHigh);
    TrackCutsProton->SetPtRange(Proton_pT_VarHigh,4.05);
    TrackCutsAntiProton->SetPtRange(Proton_pT_VarHigh,4.05);
 
    TrackPosKaonCuts->SetPIDkd(true, UseRamonaCut,Kaon_SigmaCombined_VarHigh, Kaon_SigmaTPC_VarLow, Kaon_SigmaExclusion_VarHigh);
    TrackNegKaonCuts->SetPIDkd(true, UseRamonaCut,Kaon_SigmaCombined_VarHigh, Kaon_SigmaTPC_VarLow, Kaon_SigmaExclusion_VarHigh);
    TrackPosKaonCuts->SetPtRange(Kaon_pT_VarHigh,pThighK);
    TrackNegKaonCuts->SetPtRange(Kaon_pT_VarHigh,pThighK);
 
    config->SetDeltaPhiMax(DPhi_VarHigh);
    config->SetDeltaEtaMax(DEta_VarLow);
 
  } else if (suffix == "42"){
 
    TrackCutsProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarHigh);
    TrackCutsAntiProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarHigh);
    TrackCutsProton->SetPtRange(Proton_pT_VarHigh,4.05);
    TrackCutsAntiProton->SetPtRange(Proton_pT_VarHigh,4.05);
    TrackCutsProton->SetNClsTPC(Proton_Clusters_VarHigh);
    TrackCutsAntiProton->SetNClsTPC(Proton_Clusters_VarHigh);
 
    TrackPosKaonCuts->SetPIDkd(true, UseRamonaCut,Kaon_SigmaCombined_Default, Kaon_SigmaTPC_VarHigh, Kaon_SigmaExclusion_VarLow);
    TrackNegKaonCuts->SetPIDkd(true, UseRamonaCut,Kaon_SigmaCombined_Default, Kaon_SigmaTPC_VarHigh, Kaon_SigmaExclusion_VarLow);
    TrackPosKaonCuts->SetEtaRange(-Kaon_Eta_VarLow,Kaon_Eta_VarLow);
    TrackNegKaonCuts->SetEtaRange(-Kaon_Eta_VarLow,Kaon_Eta_VarLow);
    TrackPosKaonCuts->SetPtRange(Kaon_pT_VarHigh,pThighK);
    TrackNegKaonCuts->SetPtRange(Kaon_pT_VarHigh,pThighK);
    TrackPosKaonCuts->SetNClsTPC(Kaon_Clusters_VarHigh);
    TrackNegKaonCuts->SetNClsTPC(Kaon_Clusters_VarHigh);
 
    config->SetDeltaPhiMax(DPhi_VarLow);
 
  } else if (suffix == "43"){
 
    TrackCutsProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarHigh);
    TrackCutsAntiProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarHigh);
    TrackCutsProton->SetEtaRange(-Proton_Eta_VarLow,Proton_Eta_VarLow);
    TrackCutsAntiProton->SetEtaRange(-Proton_Eta_VarLow,Proton_Eta_VarLow);
    TrackCutsProton->SetNClsTPC(Proton_Clusters_VarLow);
    TrackCutsAntiProton->SetNClsTPC(Proton_Clusters_VarLow);
 
    TrackPosKaonCuts->SetPIDkd(true, UseRamonaCut,Kaon_SigmaCombined_VarHigh, Kaon_SigmaTPC_VarHigh, Kaon_SigmaExclusion_Default);
    TrackNegKaonCuts->SetPIDkd(true, UseRamonaCut,Kaon_SigmaCombined_VarHigh, Kaon_SigmaTPC_VarHigh, Kaon_SigmaExclusion_Default);
    TrackPosKaonCuts->SetPtRange(Kaon_pT_VarLow,pThighK);
    TrackNegKaonCuts->SetPtRange(Kaon_pT_VarLow,pThighK);
    TrackPosKaonCuts->SetNClsTPC(Kaon_Clusters_VarHigh);
    TrackNegKaonCuts->SetNClsTPC(Kaon_Clusters_VarHigh);
 
    config->SetDeltaPhiMax(DPhi_VarLow);
    config->SetDeltaEtaMax(DEta_VarHigh);
 
  } else if (suffix == "44"){
 
    TrackCutsProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarHigh);
    TrackCutsAntiProton->SetPID(AliPID::kProton, 0.75,Proton_Sigma_VarHigh);
    TrackCutsProton->SetPtRange(Proton_pT_VarHigh,4.05);
    TrackCutsAntiProton->SetPtRange(Proton_pT_VarHigh,4.05);
 
    TrackPosKaonCuts->SetPIDkd(true, UseRamonaCut,Kaon_SigmaCombined_VarLow, Kaon_SigmaTPC_VarHigh, Kaon_SigmaExclusion_Default);
    TrackNegKaonCuts->SetPIDkd(true, UseRamonaCut,Kaon_SigmaCombined_VarLow, Kaon_SigmaTPC_VarHigh, Kaon_SigmaExclusion_Default);
    TrackPosKaonCuts->SetEtaRange(-Kaon_Eta_VarHigh,Kaon_Eta_VarHigh);
    TrackNegKaonCuts->SetEtaRange(-Kaon_Eta_VarHigh,Kaon_Eta_VarHigh);
    TrackPosKaonCuts->SetPtRange(Kaon_pT_VarHigh,pThighK);
    TrackNegKaonCuts->SetPtRange(Kaon_pT_VarHigh,pThighK);
 
    config->SetDeltaPhiMax(DPhi_VarLow);
 
  }

  //---------------------------------------------------- Config of Output Containers

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

  TString TrackCutsKaonName = Form("%s_TrackKaon_%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *coutputTrkCutsKaon = mgr->CreateContainer(
        TrackCutsKaonName.Data(), TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), TrackCutsKaonName.Data()));

  TString AntiTrackCutsKaonName = Form("%s_AntiTrackKaon_%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *coutputAntiTrkCutsKaon = mgr->CreateContainer(
        AntiTrackCutsKaonName.Data(), TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), AntiTrackCutsKaonName.Data()));

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

  AliAnalysisDataContainer *coutputTrkCutsMC;
  AliAnalysisDataContainer *coutputAntiTrkCutsMC;
  AliAnalysisDataContainer *coutputv0CutsMC;
  AliAnalysisDataContainer *coutputAntiv0CutsMC;

  if (isMC) {
    TString TrkCutsMCName = Form("%s_ProtonMC_%s", addon.Data(), suffix.Data());
    coutputTrkCutsMC = mgr->CreateContainer(
                         TrkCutsMCName.Data(),
                         TList::Class(),
                         AliAnalysisManager::kOutputContainer,
                         Form("%s:%s", file.Data(), TrkCutsMCName.Data()));

    TString AntiTrkCutsMCName = Form("%s_AntiProtonMC_%s", addon.Data(),
                                     suffix.Data());
    coutputAntiTrkCutsMC = mgr->CreateContainer(
                             //@suppress("Invalid arguments") it works ffs
                             AntiTrkCutsMCName.Data(),
                             TList::Class(),
                             AliAnalysisManager::kOutputContainer,
                             Form("%s:%s", file.Data(), AntiTrkCutsMCName.Data()));

    TString v0CutsMCName = Form("%s_KaonMC_%s", addon.Data(), suffix.Data());
    coutputv0CutsMC = mgr->CreateContainer(
                        //@suppress("Invalid arguments") it works ffs
                        v0CutsMCName.Data(),
                        TList::Class(),
                        AliAnalysisManager::kOutputContainer,
                        Form("%s:%s", file.Data(), v0CutsMCName.Data()));

    TString Antiv0CutsMCName = Form("%s_AntiKaonMC_%s", addon.Data(),
                                    suffix.Data());
    coutputAntiv0CutsMC = mgr->CreateContainer(
                            //@suppress("Invalid arguments") it works ffs
                            Antiv0CutsMCName.Data(),
                            TList::Class(),
                            AliAnalysisManager::kOutputContainer,
                            Form("%s:%s", file.Data(), Antiv0CutsMCName.Data()));
  }

 //-----------------------------------------------------------------------------

 AliAnalysisTaskNanoFemtoProtonKaonPlus *taskNano;
 AliAnalysisTaskFemtoProtonKaonPlus *taskAOD;

  if(isNano){
    taskNano = new AliAnalysisTaskNanoFemtoProtonKaonPlus("FemtoDreamDefault", isMC);
    
    if (trigger == "kINT7") {
      taskNano->SelectCollisionCandidates(AliVEvent::kINT7);
      std::cout << "Added kINT7 Trigger \n";
    } else if (trigger == "kHM") {
      taskNano->SelectCollisionCandidates(AliVEvent::kHighMultV0);
      std::cout << "Added kHighMult Trigger \n";
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
    if (!fullBlastQA) {
      taskNano->SetRunTaskLightWeight(true);
    }
    taskNano->SetEventCuts(evtCuts);
    taskNano->SetTrackCutsKaon(TrackPosKaonCuts);
    taskNano->SetTrackCutsAntiKaon(TrackNegKaonCuts);
    taskNano->SetTrackCutsProton(TrackCutsProton);
    taskNano->SetTrackCutsAntiProton(TrackCutsAntiProton);
    taskNano->SetCollectionConfig(config);
    taskNano->SetDoPairCleaning(DoPairCleaning);

    mgr->AddTask(taskNano);

    mgr->ConnectInput(taskNano, 0, cinput);
    mgr->ConnectOutput(taskNano, 1, coutputEvtCuts);
    mgr->ConnectOutput(taskNano, 2, coutputTrkCuts);
    mgr->ConnectOutput(taskNano, 3, coutputAntiTrkCuts);
    mgr->ConnectOutput(taskNano, 4, coutputTrkCutsKaon);
    mgr->ConnectOutput(taskNano, 5, coutputAntiTrkCutsKaon);
    mgr->ConnectOutput(taskNano, 6, coutputResults);
    mgr->ConnectOutput(taskNano, 7, coutputResultsQA);
    if(isMC){
      mgr->ConnectOutput(taskNano, 8, coutputTrkCutsMC);
      mgr->ConnectOutput(taskNano, 9, coutputAntiTrkCutsMC);
      mgr->ConnectOutput(taskNano, 10, coutputv0CutsMC);
      mgr->ConnectOutput(taskNano, 11, coutputAntiv0CutsMC);
    }
  } else {
    taskAOD = new AliAnalysisTaskFemtoProtonKaonPlus("FemtoDreamDefault", isMC);
    
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
    taskAOD->SetTrackCutsKaon(TrackPosKaonCuts);
    taskAOD->SetTrackCutsAntiKaon(TrackNegKaonCuts);
    taskAOD->SetTrackCutsProton(TrackCutsProton);
    taskAOD->SetTrackCutsAntiProton(TrackCutsAntiProton);
    taskAOD->SetCollectionConfig(config);
    taskAOD->SetDoPairCleaning(DoPairCleaning);

    mgr->AddTask(taskAOD);

    mgr->ConnectInput(taskAOD, 0, cinput);
    mgr->ConnectOutput(taskAOD, 1, coutputEvtCuts);
    mgr->ConnectOutput(taskAOD, 2, coutputTrkCuts);
    mgr->ConnectOutput(taskAOD, 3, coutputAntiTrkCuts);
    mgr->ConnectOutput(taskAOD, 4, coutputTrkCutsKaon);
    mgr->ConnectOutput(taskAOD, 5, coutputAntiTrkCutsKaon);
    mgr->ConnectOutput(taskAOD, 6, coutputResults);
    mgr->ConnectOutput(taskAOD, 7, coutputResultsQA);
    if(isMC){
      mgr->ConnectOutput(taskAOD, 8, coutputTrkCutsMC);
      mgr->ConnectOutput(taskAOD, 9, coutputAntiTrkCutsMC);
      mgr->ConnectOutput(taskAOD, 10, coutputv0CutsMC);
      mgr->ConnectOutput(taskAOD, 11, coutputAntiv0CutsMC);
    }
  }

  if (isNano) {
    return taskNano;
  } else {
    return taskAOD;
  }
}

