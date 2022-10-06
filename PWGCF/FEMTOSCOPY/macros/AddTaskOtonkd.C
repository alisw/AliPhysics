#include "TROOT.h"
#include "TSystem.h"

AliAnalysisTaskSE* AddTaskOtonkd(int isMCint = 0,
    int KaonCut = 0,
    int DeuteronCut = 0,
    int ProtonCut = 0,
    int DoFDpairing = 0
    ) {

  // isMCint = 0 //  DATA: isMC=false 
  // isMCint = 1 //  MC STD: isMC=true,  isMCtruth=false;
  // isMCint = 2 //  MC TRUTH: isMC=true,  isMCtruth=true;
  // isMCint = -1 //  DATA AGAIN WITH SOME PROTONS AS DEUTERONS: isMC=false,  isMCtruth=false, isIncludeSomeProtons=true

  // KaonCut = 0 // Oton Std
  // KaonCut = 1 // Ramona
  // KaonCut = 2 // Oton open

  // DeuteronCut = 0 // Oton 
  // DeuteronCut = 1 // FD Std
  // DeuteronCut = 2 // FD open

  // ProtonCut = 0 // protons
  // ProtonCut = 1 // pions
  // ProtonCut = 2 // perform p-d analysis (substitute kaons in tree for protons)

  // DoFDpairing = 0 // no FD pairing
  // DoFDpairing = 1 // YES do FD pairing


  bool isMC = false;
  if(isMCint>0) isMC = true;
  bool isMCtruth=false;
  if(isMCint==2)  isMCtruth=true;
  bool isIncludeSomeProtons = false;
  if(isMCint==-1) isIncludeSomeProtons = true;
  bool doFDpairing = false;
  if(DoFDpairing==1)  doFDpairing=true;
  //add a boolean when running with pions since we don't need all of them:
  bool isPions = false;
  if(ProtonCut==1) isPions = true;
  bool DOpd = false;
  if(ProtonCut==2) DOpd = true;


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
  AliFemtoDreamTrackCuts *TrackCutsKaon = AliFemtoDreamTrackCuts::PrimKaonCuts(
        isMC, true, false, false);
  TrackCutsKaon->SetCutCharge(1);
  TrackCutsKaon->SetFilterBit(128);
  if(KaonCut==0){// cuts by Oton
   TrackCutsKaon->SetPIDkd();
  }else if(KaonCut==1){ // cuts by Ramona
   TrackCutsKaon->SetPIDkd(true,true);
   //TrackCutsKaon->SetNClsTPC(70);
   //TrackCutsKaon->SetDCAVtxZ(0.1);
  }else if(KaonCut==2){ // Open cuts for syst by Oton
   //open cuts for syst:
   TrackCutsKaon->SetPIDkd(true,false,3.5,3.5);
   TrackCutsKaon->SetEtaRange(-0.84, 0.84);
   TrackCutsKaon->SetNClsTPC(70);
   TrackCutsKaon->SetDCAVtxZ(0.24);
  }else if(KaonCut==7){ // Kaons by Oton with checkpileupSPDTOF
   TrackCutsKaon->SetPIDkd();
   TrackCutsKaon->SetCheckPileUpSPDTOF(true);
  }

  AliFemtoDreamTrackCuts *TrackCutsAntiKaon = AliFemtoDreamTrackCuts::PrimKaonCuts(
        isMC, true, false, false);
  TrackCutsAntiKaon->SetCutCharge(-1);
  TrackCutsAntiKaon->SetFilterBit(128);
  if(KaonCut==0){ // cuts by Oton
   TrackCutsAntiKaon->SetPIDkd();
  }else if(KaonCut==1){ // cuts by Ramona
   TrackCutsAntiKaon->SetPIDkd(true,true);
   //TrackCutsAntiKaon->SetNClsTPC(70);
   //TrackCutsAntiKaon->SetDCAVtxZ(0.1);
  }else if(KaonCut==2){ // Open cuts for syst by Oton
   // Oton open cuts for syst:
   TrackCutsAntiKaon->SetPIDkd(true,false,3.5,3.5);
   TrackCutsAntiKaon->SetEtaRange(-0.84, 0.84);
   TrackCutsAntiKaon->SetNClsTPC(70);
   TrackCutsAntiKaon->SetDCAVtxZ(0.24);
  }else if(KaonCut==7){ // Kaons by Oton with checkpileupSPDTOF
   TrackCutsAntiKaon->SetPIDkd();
   TrackCutsAntiKaon->SetCheckPileUpSPDTOF(true);
  }

  //deuterons
  AliFemtoDreamTrackCuts *TrackCutsDeuteron = AliFemtoDreamTrackCuts::PrimDeuteronCuts(
        isMC, true, false, false);
  TrackCutsDeuteron->SetCutCharge(1);
  TrackCutsDeuteron->SetFilterBit(128);
  if(DeuteronCut==0){ // cuts by Oton
   TrackCutsDeuteron->SetPIDkd(false,false,3,3);
  }else if(DeuteronCut==1){ // cuts by FemtoDream
   TrackCutsDeuteron->SetPtRange(0.5,1.4);
  }else if(DeuteronCut==2){ // Open cuts for syst by FemtoDream
   TrackCutsDeuteron->SetPtRange(0.4,1.5);
   TrackCutsDeuteron->SetEtaRange(-0.84, 0.84);
   TrackCutsDeuteron->SetNClsTPC(70);
   TrackCutsDeuteron->SetDCAVtxZ(0.25);
   TrackCutsDeuteron->SetPID(AliPID::kDeuteron,1.4, 3.5);
  }else if(DeuteronCut==3){ // Deuterons by Oton with open cuts for PID/DCA fits
   TrackCutsDeuteron->SetPIDkd(false,false,3,5);
   TrackCutsDeuteron->SetDCAVtxXY(0.5);
  }else if(DeuteronCut==4){ // cuts by FemtoDream with tight PID req
   TrackCutsDeuteron->SetPtRange(0.5,1.4);
   TrackCutsDeuteron->SetPID(AliPID::kDeuteron,1.4, 1.4);//1.4 sigmas as well!
  }else if(DeuteronCut==5){ // Deuterons by Oton with open cuts SYST
   TrackCutsDeuteron->SetPIDkd(false,false,3.5,3.5);
   TrackCutsDeuteron->SetEtaRange(-0.84, 0.84);
   TrackCutsDeuteron->SetNClsTPC(70);
   TrackCutsDeuteron->SetDCAVtxZ(0.24);
  }else if(DeuteronCut==6){ // Deuterons by Oton with cut in min pt
   TrackCutsDeuteron->SetPIDkd(false,false,3,3);
   TrackCutsDeuteron->SetPtRange(0.8,2.9);
  }else if(DeuteronCut==7){ // Deuterons by Oton with checkpileupSPDTOF
   TrackCutsDeuteron->SetPIDkd(false,false,3,3);
   TrackCutsDeuteron->SetCheckPileUpSPDTOF(true);
  }

  AliFemtoDreamTrackCuts *TrackCutsAntiDeuteron = AliFemtoDreamTrackCuts::PrimDeuteronCuts(
        isMC, true, false, false);
  TrackCutsAntiDeuteron->SetCutCharge(-1);
  TrackCutsAntiDeuteron->SetFilterBit(128);
  if(DeuteronCut==0){ // cuts by Oton
   TrackCutsAntiDeuteron->SetPIDkd(false,false,3,3);
  }else if(DeuteronCut==1){ // cuts by FemtoDream
   TrackCutsAntiDeuteron->SetPtRange(0.5,1.4);
  }else if(DeuteronCut==2){ // Open cuts for syst by FemtoDream
   TrackCutsAntiDeuteron->SetPtRange(0.4,1.5);
   TrackCutsAntiDeuteron->SetEtaRange(-0.84, 0.84);
   TrackCutsAntiDeuteron->SetNClsTPC(70);
   TrackCutsAntiDeuteron->SetDCAVtxZ(0.25);
   TrackCutsAntiDeuteron->SetPID(AliPID::kDeuteron,1.4, 3.5);
  }else if(DeuteronCut==3){ // Deuterons by Oton with open cuts for PID/DCA fits
   TrackCutsAntiDeuteron->SetPIDkd(false,false,3,5);
   TrackCutsAntiDeuteron->SetDCAVtxXY(0.5);
  }else if(DeuteronCut==4){ // cuts by FemtoDream with tight PID req
   TrackCutsAntiDeuteron->SetPtRange(0.5,1.4);
   TrackCutsAntiDeuteron->SetPID(AliPID::kDeuteron,1.4,1.4);//1.4 sigmas as well!
  }else if(DeuteronCut==5){ // Deuterons by Oton with open cuts SYST
   TrackCutsAntiDeuteron->SetPIDkd(false,false,3.5,3.5);
   TrackCutsAntiDeuteron->SetEtaRange(-0.84, 0.84);
   TrackCutsAntiDeuteron->SetNClsTPC(70);
   TrackCutsAntiDeuteron->SetDCAVtxZ(0.24);
  }else if(DeuteronCut==6){ // Deuterons by Oton with cut in min pt
   TrackCutsAntiDeuteron->SetPIDkd(false,false,3,3);
   TrackCutsAntiDeuteron->SetPtRange(0.8,2.9);
  }else if(DeuteronCut==7){ // Deuterons by Oton with checkpileupSPDTOF
   TrackCutsAntiDeuteron->SetPIDkd(false,false,3,3);
   TrackCutsAntiDeuteron->SetCheckPileUpSPDTOF(true);
  }

  //protons
  AliFemtoDreamTrackCuts *TrackCutsProton;
  if(ProtonCut==0||ProtonCut==2){ //use protons
   TrackCutsProton = AliFemtoDreamTrackCuts::PrimProtonCuts(isMC, true, false, false);
  }
  if(ProtonCut==1){ //use pions
   TrackCutsProton = new AliFemtoDreamTrackCuts();
   TrackCutsProton->SetIsMonteCarlo(isMC);
   TrackCutsProton->SetPtRange(0.14, 4.0);
   TrackCutsProton->SetEtaRange(-0.8, 0.8);
   TrackCutsProton->SetNClsTPC(80);
   TrackCutsProton->SetDCAReCalculation(true);//Get the dca from the PropagateToVetex
   TrackCutsProton->SetFilterBit(128);//96); // Filterbit 5+6
   TrackCutsProton->SetDCAVtxZ(0.3);
   TrackCutsProton->SetDCAVtxXY(0.3);
   TrackCutsProton->SetNClsTPC(80); // In Indico + additrotonal Chi²/NDF <4
   TrackCutsProton->SetPID(AliPID::kPion, 0.5);
   TrackCutsProton->SetRejLowPtPionsTOF(false);
   TrackCutsProton->SetMinimalBooking(false);
   TrackCutsProton->SetPlotDCADist(true);
   TrackCutsProton->SetCheckPileUpSPDTOF(true);
  }
  TrackCutsProton->SetCutCharge(1);

  AliFemtoDreamTrackCuts *TrackCutsAntiProton;
  if(ProtonCut==0||ProtonCut==2){ //use protons
   TrackCutsAntiProton = AliFemtoDreamTrackCuts::PrimProtonCuts(isMC, true, false, false);
  }
  if(ProtonCut==1){ //use pions
   TrackCutsAntiProton = new AliFemtoDreamTrackCuts();
   TrackCutsAntiProton->SetIsMonteCarlo(isMC);
   TrackCutsAntiProton->SetPtRange(0.14, 4.0);
   TrackCutsAntiProton->SetEtaRange(-0.8, 0.8);
   TrackCutsAntiProton->SetNClsTPC(80);
   TrackCutsAntiProton->SetDCAReCalculation(true);
   TrackCutsAntiProton->SetFilterBit(128);//96);
   TrackCutsAntiProton->SetDCAVtxZ(0.3);
   TrackCutsAntiProton->SetDCAVtxXY(0.3);
   TrackCutsAntiProton->SetNClsTPC(80);
   TrackCutsAntiProton->SetPID(AliPID::kPion, 0.5);
   TrackCutsAntiProton->SetRejLowPtPionsTOF(false);
   TrackCutsAntiProton->SetMinimalBooking(false);
   TrackCutsAntiProton->SetPlotDCADist(true);
   TrackCutsAntiProton->SetCheckPileUpSPDTOF(true);
  }
  TrackCutsAntiProton->SetCutCharge(-1);


//TO WHICH ORDER CORRESPONDS THIS???????????
  std::vector<int> PDGParticles;
  PDGParticles.push_back(321);//k
  PDGParticles.push_back(321);//k
  if(ProtonCut==0||ProtonCut==2){ //use protons
   PDGParticles.push_back(2212);
   PDGParticles.push_back(2212);
  }
  if(ProtonCut==1){ //use pions
   PDGParticles.push_back(211);
   PDGParticles.push_back(211);
  }
  PDGParticles.push_back(1000010020);//d
  PDGParticles.push_back(1000010020);//d

/*
// if we need ITS//  WOW!!!!
  TrackCutsDeuteron->SetITSnSigmaCut(true);
  TrackCutsAntiDeuteron->SetITSnSigmaCut(true);
  TrackCutsDeuteron->SetCutITSPID(1.4, -2.0, 1e30);
  TrackCutsAntiDeuteron->SetCutITSPID(1.4, -2.0, 1e30);
*/

  std::vector<float> kMin;
  std::vector<float> kMax;
  std::vector<int> NBins;
  std::vector<bool> closeRejection;
 // std::vector<float> mTBins;
//  mTBins.push_back(1.14);
//  mTBins.push_back(1.26);
//  mTBins.push_back(999.);
  std::vector<int> pairQA;
  //pairs:
//k k   0
//k bark 1
//k p 2
//k barp 3
//k d 4
//k bard 5
//bark bark 6
//bark p 7
//bark barp 8
//bark d 9
//bark bard 10
  // p p           11
  // p barp        12
  // p d            13
  // p bard        14
  // barp barp    15
  // barp d        16
  // barp bard    17
  // d d            18
  // d bard        19
  // bard bard    20

//this shit is still all by hand and I want to die again:
/////////////////////////////////////////////////////////
  const int nPairs = 21;//
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

  closeRejection[0] = true;  // k k
  closeRejection[2] = true;  // k p
  closeRejection[4] = true;  // k d
  closeRejection[5] = true;  // k bard
  closeRejection[6] = true;  // bark bark
  closeRejection[8] = true;  // bark barp
  closeRejection[9] = true;  // bark d
  closeRejection[10] = true;  // bark bard
  closeRejection[11] = true;  // pp
  closeRejection[13] = true;  // pd
  closeRejection[14] = true;  // pbard
  closeRejection[15] = true;  // barp barp
  closeRejection[16] = true;  // barp d
  closeRejection[17] = true;  // barp bard
  closeRejection[18] = true;  // dd
  closeRejection[20] = true;  // bard bar

  pairQA[0] = 11;  // k k
  pairQA[2] = 11;  // k p
  pairQA[4] = 11;  // k d
  pairQA[6] = 11;  // bark bark
  pairQA[8] = 11;  // bark barp
  pairQA[10] = 11;  // bark bard
  pairQA[11] = 11;  // pp
  pairQA[13] = 11;  // pd
  pairQA[15] = 11;  // barp barp
  pairQA[17] = 11;  // barp bard
  pairQA[18] = 11;  // dd
  pairQA[20] = 11;  // bard bar

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
    TrackCutsProton->SetMinimalBooking(true);
    TrackCutsAntiProton->SetMinimalBooking(true);
    TrackCutsDeuteron->SetMinimalBooking(true);
    //TrackCutsDeuteronMass->SetMinimalBooking(true);
    TrackCutsAntiDeuteron->SetMinimalBooking(true);
    //TrackCutsAntiDeuteronMass->SetMinimalBooking(true);
    config->SetMinimalBookingME(true);
    config->SetMinimalBookingSample(true);
  }

  //Define here the analysis task
  AliAnalysisTaskOtonkd *task =
   new AliAnalysisTaskOtonkd("ThisNameApparentlyStillUseless", isMC, isMCtruth, isIncludeSomeProtons, isPions, doFDpairing, DOpd);
  task->SelectCollisionCandidates(AliVEvent::kHighMultV0);
  if (!fullBlastQA) {
    task->SetRunTaskLightWeight(true);
  }

  //Throw all our settings to the task
  task->SetEventCuts(evtCuts);
  task->SetTrackCutsDeuteron(TrackCutsDeuteron);
  //task->SetTrackCutsDeuteronMass(TrackCutsDeuteronMass);
  task->SetTrackCutsAntiDeuteron(TrackCutsAntiDeuteron);
  //task->SetTrackCutsAntiDeuteronMass(TrackCutsAntiDeuteronMass);
  task->SetTrackCutsProton(TrackCutsProton);
  task->SetTrackCutsAntiProton(TrackCutsAntiProton);
  task->SetTrackCutsKaon(TrackCutsKaon);
  task->SetTrackCutsAntiKaon(TrackCutsAntiKaon);
  task->SetCollectionConfig(config);

  mgr->AddTask(task);

  TString addon = "KD";

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

  TString TrackCutsName = Form("%sProton%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *coutputTrkCuts = mgr->CreateContainer(
        TrackCutsName.Data(), TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), TrackCutsName.Data()));
  mgr->ConnectOutput(task, 4, coutputTrkCuts);

  TString AntiTrackCutsName = Form("%sAntiProton%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *coutputAntiTrkCuts = mgr->CreateContainer(
        AntiTrackCutsName.Data(), TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), AntiTrackCutsName.Data()));
  mgr->ConnectOutput(task, 5, coutputAntiTrkCuts);

  TString TrackCutsDeuteronName = Form("%sDeuteron%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *coutputTrkCutsDeuteron = mgr->CreateContainer(
        TrackCutsDeuteronName.Data(), TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), TrackCutsDeuteronName.Data()));
  mgr->ConnectOutput(task, 6, coutputTrkCutsDeuteron);

  TString AntiTrackCutsDeuteronName = Form("%sAntiDeuteron%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *coutputAntiTrkCutsDeuteron = mgr->CreateContainer(
        AntiTrackCutsDeuteronName.Data(), TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), AntiTrackCutsDeuteronName.Data()));
  mgr->ConnectOutput(task, 7, coutputAntiTrkCutsDeuteron);

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


  //The tree:
  AliAnalysisDataContainer *coutputTree;
  TString TreeName = Form("%sTree",addon.Data());
  coutputTree = mgr->CreateContainer(
    //@suppress("Invalid arguments") it works ffs
    TreeName.Data(),
    TTree::Class(),
    AliAnalysisManager::kOutputContainer,
    Form("%s:%s", file.Data(), TreeName.Data()));
  mgr->ConnectOutput(task, 10, coutputTree);

  return task;
}

