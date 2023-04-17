#include "TROOT.h"
#include "TSystem.h"

AliAnalysisTaskSE* AddTaskOtonXx(int isMCint = 0,
    int KaonCut = 0,
    int XiCut = 0,
    int ProtonCut = 0,
    int DoFDpairing = 0
    ) {


bool isMC = false;

//do pions:
  // KaonCut = 0 // Kaons
  // KaonCut = 1 // pions
  bool isPi = false;
  if(KaonCut==1) isPi = true;

//do omegas
  // XiCut = 0 // Xi
  // XiCut = 1 // Omega
  bool isOmega = false;
  if(XiCut==1) isOmega = true;


  const char fullBlastQA = true; //moved from arguments
  const char *cutVariation = "0"; //moved from arguments, for the moment I don't use it

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

//DO I NEED to incLUDE A EVENTHANDLER????????
//  AliVEventHandler *inputHandler = mgr->GetInputEventHandler();
  AliFemtoDreamEventCuts *evtCuts = AliFemtoDreamEventCuts::StandardCutsRun2();
  evtCuts->CleanUpMult(false, false, false, true);


  //kaons:
  AliFemtoDreamTrackCuts *TrackCutsKaon;
  if(KaonCut!=1){ // std, use kaons
   TrackCutsKaon = AliFemtoDreamTrackCuts::PrimKaonCuts(isMC, true, false, false);
   TrackCutsKaon->SetPIDkd();
  }else{ //use pions
/*
   TrackCutsKaon = AliFemtoDreamTrackCuts::PrimPionCuts(isMC, true, false, false);
*/
   TrackCutsKaon = new AliFemtoDreamTrackCuts();
   TrackCutsKaon->SetIsMonteCarlo(isMC);
   TrackCutsKaon->SetPtRange(0.14, 4.0);
   TrackCutsKaon->SetEtaRange(-0.8, 0.8);
   TrackCutsKaon->SetNClsTPC(80);
   TrackCutsKaon->SetDCAReCalculation(true);//Get the dca from the PropagateToVetex
   TrackCutsKaon->SetDCAVtxZ(0.3);
   TrackCutsKaon->SetDCAVtxXY(0.3);
   TrackCutsKaon->SetNClsTPC(80); // In Indico + additrotonal Chi²/NDF <4
   TrackCutsKaon->SetPID(AliPID::kPion, 0.5);
   TrackCutsKaon->SetRejLowPtPionsTOF(false);
   TrackCutsKaon->SetMinimalBooking(false);
   TrackCutsKaon->SetPlotDCADist(true);
   TrackCutsKaon->SetCheckPileUpSPDTOF(true);
  }
  TrackCutsKaon->SetCutCharge(1);
  TrackCutsKaon->SetFilterBit(128);

  //antikaons
  AliFemtoDreamTrackCuts *TrackCutsAntiKaon;
  if(KaonCut!=1){ // std, use kaons
   TrackCutsAntiKaon = AliFemtoDreamTrackCuts::PrimKaonCuts(isMC, true, false, false);
   TrackCutsAntiKaon->SetPIDkd();
  }else{ //use pions
/*
   TrackCutsAntiKaon = AliFemtoDreamTrackCuts::PrimPionCuts(isMC, true, false, false); 
*/
   TrackCutsAntiKaon = new AliFemtoDreamTrackCuts();
   TrackCutsAntiKaon->SetIsMonteCarlo(isMC);
   TrackCutsAntiKaon->SetPtRange(0.14, 4.0);
   TrackCutsAntiKaon->SetEtaRange(-0.8, 0.8);
   TrackCutsAntiKaon->SetNClsTPC(80);
   TrackCutsAntiKaon->SetDCAReCalculation(true);
   TrackCutsAntiKaon->SetDCAVtxZ(0.3);
   TrackCutsAntiKaon->SetDCAVtxXY(0.3);
   TrackCutsAntiKaon->SetNClsTPC(80);
   TrackCutsAntiKaon->SetPID(AliPID::kPion, 0.5);
   TrackCutsAntiKaon->SetRejLowPtPionsTOF(false);
   TrackCutsAntiKaon->SetMinimalBooking(false);
   TrackCutsAntiKaon->SetPlotDCADist(true);
   TrackCutsAntiKaon->SetCheckPileUpSPDTOF(true);
  }
  TrackCutsAntiKaon->SetCutCharge(-1);
  TrackCutsAntiKaon->SetFilterBit(128);


  //Cascade Cuts 
  AliFemtoDreamCascadeCuts* CascadeXiCuts;
  if(!isOmega){
   CascadeXiCuts = AliFemtoDreamCascadeCuts::XiCuts(isMC, false);
  }else{
   CascadeXiCuts = AliFemtoDreamCascadeCuts::OmegaCuts(isMC, false);
  }
  CascadeXiCuts->SetXiCharge(-1);
  AliFemtoDreamTrackCuts *XiNegCuts = AliFemtoDreamTrackCuts::Xiv0PionCuts(isMC, true, false);
  XiNegCuts->SetCheckTPCRefit(false);
  XiNegCuts->SetCheckPileUp(false);  //release timing information 
  AliFemtoDreamTrackCuts *XiPosCuts = AliFemtoDreamTrackCuts::Xiv0ProtonCuts(isMC, true, false);
  XiPosCuts->SetCheckTPCRefit(false);
  AliFemtoDreamTrackCuts *XiBachCuts;
  if(!isOmega){
   XiBachCuts = AliFemtoDreamTrackCuts::XiBachPionCuts(isMC, true, false);
  }else{
   XiBachCuts = AliFemtoDreamTrackCuts::OmegaBachKaonCuts(isMC, true, false);
  }
  XiBachCuts->SetCheckTPCRefit(false);
  XiNegCuts->SetEtaRange(-0.8, 0.8);
  XiPosCuts->SetEtaRange(-0.8, 0.8);
  XiBachCuts->SetEtaRange(-0.8, 0.8);
  CascadeXiCuts->Setv0Negcuts(XiNegCuts);
  CascadeXiCuts->Setv0PosCuts(XiPosCuts);
  CascadeXiCuts->SetBachCuts(XiBachCuts);
  if(!isOmega){
   CascadeXiCuts->SetPDGCodeCasc(3312);
  }else{
   CascadeXiCuts->SetPDGCodeCasc(3334);
  }
  CascadeXiCuts->SetPDGCodev0(3122);
  CascadeXiCuts->SetPDGCodePosDaug(2212);
  CascadeXiCuts->SetPDGCodeNegDaug(-211);
  if(!isOmega){
   CascadeXiCuts->SetPDGCodeBach(-211);
  }else{
   CascadeXiCuts->SetPDGCodeBach(-321);
  }


  //AntiCascade cuts 
  AliFemtoDreamCascadeCuts* AntiCascadeXiCuts;
  if(!isOmega){
   AntiCascadeXiCuts = AliFemtoDreamCascadeCuts::XiCuts(isMC, false);
  }else{
   AntiCascadeXiCuts = AliFemtoDreamCascadeCuts::OmegaCuts(isMC, false);
  }
  AntiCascadeXiCuts->SetXiCharge(1);
  AliFemtoDreamTrackCuts *AntiXiNegCuts = AliFemtoDreamTrackCuts::Xiv0ProtonCuts(isMC, true, false);
  AntiXiNegCuts->SetCutCharge(-1);
  AntiXiNegCuts->SetCheckTPCRefit(false);
  AliFemtoDreamTrackCuts *AntiXiPosCuts = AliFemtoDreamTrackCuts::Xiv0PionCuts(isMC, true, false);
  AntiXiPosCuts->SetCutCharge(1);
  AntiXiPosCuts->SetCheckTPCRefit(false);
  AntiXiPosCuts->SetCheckPileUp(false);  //release timing information 
  AliFemtoDreamTrackCuts *AntiXiBachCuts;
  if(!isOmega){
   AntiXiBachCuts =AliFemtoDreamTrackCuts::XiBachPionCuts(isMC, true, false);
  }else{
   AntiXiBachCuts =AliFemtoDreamTrackCuts::OmegaBachKaonCuts(isMC, true, false);
  }
  AntiXiBachCuts->SetCutCharge(1);
  AntiXiBachCuts->SetCheckTPCRefit(false); // THREE TIMES ??? why ???
  AntiXiNegCuts->SetEtaRange(-0.8, 0.8);
  AntiXiPosCuts->SetEtaRange(-0.8, 0.8);
  AntiXiBachCuts->SetEtaRange(-0.8, 0.8);
  AntiCascadeXiCuts->Setv0Negcuts(AntiXiNegCuts);
  AntiCascadeXiCuts->Setv0PosCuts(AntiXiPosCuts);
  AntiCascadeXiCuts->SetBachCuts(AntiXiBachCuts);
  if(!isOmega){
   AntiCascadeXiCuts->SetPDGCodeCasc(-3312);
  }else{
   AntiCascadeXiCuts->SetPDGCodeCasc(-3334);
  }
  AntiCascadeXiCuts->SetPDGCodev0(-3122);
  AntiCascadeXiCuts->SetPDGCodePosDaug(211);
  AntiCascadeXiCuts->SetPDGCodeNegDaug(-2212);
  if(!isOmega){
   AntiCascadeXiCuts->SetPDGCodeBach(211);
  }else{
   AntiCascadeXiCuts->SetPDGCodeBach(321);
  }


  if(!isOmega){
   //Set the Xi cuts: // THIS ARE SUPPOSED TO BE GEORGIOS L-Xi CUTS
   CascadeXiCuts->SetCutXiDaughterDCA(1.5);
   AntiCascadeXiCuts->SetCutXiDaughterDCA(1.5);
   CascadeXiCuts->SetCutXiMinDistBachToPrimVtx(0.05);
   AntiCascadeXiCuts->SetCutXiMinDistBachToPrimVtx(0.05);
   CascadeXiCuts->SetCutv0MinDistToPrimVtx(0.07);
   AntiCascadeXiCuts->SetCutv0MinDistToPrimVtx(0.07);
   CascadeXiCuts->SetCutv0MinDaugDistToPrimVtx(0.05);
   AntiCascadeXiCuts->SetCutv0MinDaugDistToPrimVtx(0.05);
   CascadeXiCuts->SetCutXiCPA(0.97);
   AntiCascadeXiCuts->SetCutXiCPA(0.97);
   CascadeXiCuts->SetCutXiTransverseRadius(0.8, 200);
   AntiCascadeXiCuts->SetCutXiTransverseRadius(0.8, 200);
   CascadeXiCuts->SetCutv0TransverseRadius(1.4, 200);
   AntiCascadeXiCuts->SetCutv0TransverseRadius(1.4, 200);
   CascadeXiCuts->SetXiMassRange(1.322, 0.015);      // mass is open!!!
   AntiCascadeXiCuts->SetXiMassRange(1.322, 0.015);  // mass is open!!!
  }else{
Float_t XiDaughterDCA = .8; //std
Float_t v0MaxDaughterDCA = 1.2; //std
Float_t v0MinDistToPrimVtx = 0.06; //std
Float_t XiCPA = 0.995; //std
Float_t v0CPA = 0.97; //std
Float_t XiMinDistBachToPrimVtx = 0.04; //std
Float_t v0MinDaugDistToPrimVtx = 0.04; //std
Float_t Eta = 0.8; //std
Float_t nSigma = 4.; //std
Int_t TPCcrR = 60; //std
Float_t CrF = 0.75; //std
Float_t v0TransverseRadius = .000001; //std, no cut (cut in Vertexer?)
Float_t XiTransverseRadius = .000001; //std, no cut (cut in Vertexer?)
Float_t v0MassRange = 0.006; //std
Float_t RejectXis = 0.008; //std
Float_t PtRangeXi = 0.000001; //std, no cut (cut in Vertexer?)
bool PionCheckPileUp = false; //std

  CascadeXiCuts->SetXiMassRange(1.67245, 0.015);                      // mass is open !!!
  CascadeXiCuts->SetCutXiDaughterDCA(XiDaughterDCA);
  CascadeXiCuts->SetCutv0MaxDaughterDCA(v0MaxDaughterDCA);
  CascadeXiCuts->SetCutv0MinDistToPrimVtx(v0MinDistToPrimVtx);
  CascadeXiCuts->SetCutXiCPA(XiCPA);
  CascadeXiCuts->SetCutv0CPA(v0CPA);
  CascadeXiCuts->SetCutXiMinDistBachToPrimVtx(XiMinDistBachToPrimVtx);
  CascadeXiCuts->SetCutv0MinDaugDistToPrimVtx(v0MinDaugDistToPrimVtx);
  XiPosCuts->SetPID(AliPID::kProton, 999., nSigma);
  XiNegCuts->SetPID(AliPID::kPion, 999., nSigma);
  //XiBachCuts->SetPID(AliPID::kKaon, 999., 4,true, 3);
  XiBachCuts->SetPID(AliPID::kKaon, 999., nSigma);
  XiPosCuts->SetEtaRange(-1.*Eta,Eta);
  XiNegCuts->SetEtaRange(-1.*Eta,Eta);
  XiBachCuts->SetEtaRange(-1.*Eta,Eta);
  XiPosCuts->SetCutTPCCrossedRows(true, TPCcrR, CrF);
  XiNegCuts->SetCutTPCCrossedRows(true, TPCcrR, CrF);
  XiBachCuts->SetCutTPCCrossedRows(true, TPCcrR, CrF);
  CascadeXiCuts->SetCutv0TransverseRadius(v0TransverseRadius, 200);
  CascadeXiCuts->SetCutXiTransverseRadius(XiTransverseRadius, 200);
  CascadeXiCuts->Setv0MassRange(1.116, v0MassRange);
  CascadeXiCuts->SetRejectXis(1.322, RejectXis);
  CascadeXiCuts->SetPtRangeXi(PtRangeXi, 999.9);
  XiNegCuts->SetCheckPileUp(PionCheckPileUp);

  AntiCascadeXiCuts->SetXiMassRange(1.67245, 0.015);                      // mass is open !!!
  AntiCascadeXiCuts->SetCutXiDaughterDCA(XiDaughterDCA);
  AntiCascadeXiCuts->SetCutv0MaxDaughterDCA(v0MaxDaughterDCA);
  AntiCascadeXiCuts->SetCutv0MinDistToPrimVtx(v0MinDistToPrimVtx);
  AntiCascadeXiCuts->SetCutXiCPA(XiCPA);
  AntiCascadeXiCuts->SetCutv0CPA(v0CPA);
  AntiCascadeXiCuts->SetCutXiMinDistBachToPrimVtx(XiMinDistBachToPrimVtx);
  AntiCascadeXiCuts->SetCutv0MinDaugDistToPrimVtx(v0MinDaugDistToPrimVtx);
  AntiXiPosCuts->SetEtaRange(-1.*Eta,Eta);
  AntiXiNegCuts->SetEtaRange(-1.*Eta,Eta);
  AntiXiBachCuts->SetEtaRange(-1.*Eta,Eta);
  AntiXiPosCuts->SetPID(AliPID::kPion, 999., nSigma);
  AntiXiNegCuts->SetPID(AliPID::kProton, 999., nSigma);
  //AntiXiBachCuts->SetPID(AliPID::kKaon, 999., 4,true,3);
  AntiXiBachCuts->SetPID(AliPID::kKaon, 999., nSigma);
  AntiXiPosCuts->SetCutTPCCrossedRows(true, TPCcrR, CrF);
  AntiXiNegCuts->SetCutTPCCrossedRows(true, TPCcrR, CrF);
  AntiXiBachCuts->SetCutTPCCrossedRows(true, TPCcrR, CrF);
  AntiCascadeXiCuts->SetCutv0TransverseRadius(v0TransverseRadius, 200);
  AntiCascadeXiCuts->SetCutXiTransverseRadius(XiTransverseRadius, 200);
  AntiCascadeXiCuts->Setv0MassRange(1.116, v0MassRange);     
  AntiCascadeXiCuts->SetRejectXis(1.322, RejectXis);
  AntiCascadeXiCuts->SetPtRangeXi(PtRangeXi, 999.9);
  AntiXiPosCuts->SetCheckPileUp(PionCheckPileUp);
}


  // Femto Collection
//TO WHICH ORDER CORRESPONDS THIS???????????
  std::vector<int> PDGParticles;
  if(KaonCut!=1){ //use Kaons
   PDGParticles.push_back(321);
   PDGParticles.push_back(321);
  }else { // use pions
   PDGParticles.push_back(211);
   PDGParticles.push_back(211);
  }
  if(!isOmega){
   PDGParticles.push_back(3312);
   PDGParticles.push_back(3312);
  }else{
   PDGParticles.push_back(3334);
   PDGParticles.push_back(3334);
  }

  std::vector<float> kMin;
  std::vector<float> kMax;
  std::vector<int> NBins;
  std::vector<bool> closeRejection;
  std::vector<int> pairQA;
  //pairs:
//k k   0
//k bark 1
//k Xi 2
//k barXi 3
//bark bark 4
//bark Xi 5
//bark barXi 6
  // Xi Xi           7
  // Xi barXi        8
  // barXi barXi    9

//this shit is still all by hand and I want to die again:
/////////////////////////////////////////////////////////
  const int nPairs = 10;//
  for (int i = 0; i < nPairs; ++i) {
    closeRejection.push_back(false);
    pairQA.push_back(0);
    if (suffix == "0") {
      NBins.push_back(1500);
      kMin.push_back(0.);
      kMax.push_back(6.);
    } else {
      NBins.push_back(250);
      kMin.push_back(0.);
      kMax.push_back(1.);
    }
  }

//  closeRejection[0] = true;  // k k

//  pairQA[0] = 13;  // k k

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
    TrackCutsKaon->SetMinimalBooking(true);
    TrackCutsAntiKaon->SetMinimalBooking(true);
    CascadeXiCuts->SetMinimalBooking(true);
    AntiCascadeXiCuts->SetMinimalBooking(true);
    config->SetMinimalBookingME(true);
    config->SetMinimalBookingSample(true);
  }

  //Define here the analysis task
  AliAnalysisTaskOtonXx *task =
   new AliAnalysisTaskOtonXx("ThisNameApparentlyStillUseless", false,false,false,isOmega,isPi);
  task->SelectCollisionCandidates(AliVEvent::kHighMultV0);
  if (!fullBlastQA) {
    task->SetRunTaskLightWeight(true);
  }

  //Throw all our settings to the task
  task->SetEventCuts(evtCuts);
  task->SetTrackCutsKaon(TrackCutsKaon);
  task->SetTrackCutsAntiKaon(TrackCutsAntiKaon);
  task->SetXiCuts(CascadeXiCuts);
  task->SetAntiXiCuts(AntiCascadeXiCuts);

  task->SetCollectionConfig(config);

  mgr->AddTask(task);

  TString addon = "XX";

  TString file = AliAnalysisManager::GetCommonFileName();
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  mgr->ConnectInput(task, 0, cinput);

  TString EvtCutsName = Form("%sEvtCuts%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *coutputEvtCuts = mgr->CreateContainer(
        EvtCutsName.Data(), TList::Class(), AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), EvtCutsName.Data()));
  mgr->ConnectOutput(task, 1, coutputEvtCuts);

  TString TrackCutsKaonName = Form("%sKaon%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *coutputTrkCutsKaon = mgr->CreateContainer(
        TrackCutsKaonName.Data(), TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), TrackCutsKaonName.Data()));
  mgr->ConnectOutput(task, 2, coutputTrkCutsKaon);

  TString AntiTrackCutsKaonName = Form("%sAntiKaon%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *coutputAntiTrkCutsKaon = mgr->CreateContainer(
        AntiTrackCutsKaonName.Data(), TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), AntiTrackCutsKaonName.Data()));
  mgr->ConnectOutput(task, 3, coutputAntiTrkCutsKaon);

  TString XiCutsName = Form("%sXi%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *coutputXiCuts = mgr->CreateContainer(
        XiCutsName.Data(), TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), XiCutsName.Data()));
  mgr->ConnectOutput(task, 4, coutputXiCuts);
 
  TString AntiXiCutsName = Form("%sAntiXi%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *coutputAntiXiCuts = mgr->CreateContainer(
        AntiXiCutsName.Data(), TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), AntiXiCutsName.Data()));
  mgr->ConnectOutput(task, 5, coutputAntiXiCuts);

  AliAnalysisDataContainer *coutputResults;
  TString ResultsName = Form("%sResults%s", addon.Data(), suffix.Data());
  coutputResults = mgr->CreateContainer(ResultsName.Data(),
                     TList::Class(), AliAnalysisManager::kOutputContainer,
                     Form("%s:%s", file.Data(), ResultsName.Data()));
  mgr->ConnectOutput(task, 6, coutputResults);

  AliAnalysisDataContainer *coutputResultsQA;
  TString ResultsQAName = Form("%sResultsQA%s", addon.Data(), suffix.Data());
  coutputResultsQA = mgr->CreateContainer(
                       //@suppress("Invalid arguments") it works ffs
                       ResultsQAName.Data(),
                       TList::Class(),
                       AliAnalysisManager::kOutputContainer,
                       Form("%s:%s", file.Data(), ResultsQAName.Data()));
  mgr->ConnectOutput(task, 7, coutputResultsQA);


  //The tree:
  AliAnalysisDataContainer *coutputTree;
  TString TreeName = Form("%sTree",addon.Data());
  coutputTree = mgr->CreateContainer(
    //@suppress("Invalid arguments") it works ffs
    TreeName.Data(),
    TTree::Class(),
    AliAnalysisManager::kOutputContainer,
    Form("%s:%s", file.Data(), TreeName.Data()));
  mgr->ConnectOutput(task, 8, coutputTree);

  return task;
}

