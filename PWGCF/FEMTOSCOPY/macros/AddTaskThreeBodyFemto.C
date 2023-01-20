#if !defined(__CINT__) || defined(__CLING__)
#include <vector>
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskThreeBodyFemto.h"
#include "AliAnalysisTaskThreeBodyFemtoAOD.h"
#include "AliFemtoDreamEventCuts.h"
#include "AliFemtoDreamTrackCuts.h"
#include "AliFemtoDreamCascadeCuts.h"
#include "AliFemtoDreamCollConfig.h"
#endif

AliAnalysisTaskSE *AddTaskThreeBodyFemto(int trigger = 0, bool fullBlastQA = true,
                                     bool isMC = false, bool isNano = true, bool triggerOn = false, 
                                     int mixingDepthFromTask = 20,
                                     float Q3Limit = 0.6, float Q3LimitSample = 3.0,float Q3LimitSample2 = 3.0, float Q3LimitFraction = 0.5, float Q3LimitSampleFraction = 0.01, float Q3LimitSampleFraction2 = 0.01,
                                     const char *cutVariation = "0", bool ClosePairRejectionForAll = false, 
                                     bool run2Body = false, int mixinfChoice = 0, bool mix21 = false,
                                     int whichTripletsToRun = 11, bool Runppp0ppL1 = false,
                                     const char *triggerVariation = "0") {



  TString suffix = TString::Format("%s", cutVariation);
  TString suffixTrigger = TString::Format("%s", triggerVariation);

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

  // Track Cuts
  AliFemtoDreamTrackCuts *TrackCuts = AliFemtoDreamTrackCuts::PrimProtonCuts(
      isMC, true, false, false);
  TrackCuts->SetFilterBit(128);
  TrackCuts->SetCutCharge(1);

  AliFemtoDreamTrackCuts *AntiTrackCuts =
      AliFemtoDreamTrackCuts::PrimProtonCuts(isMC, true, false, false);
  AntiTrackCuts->SetFilterBit(128);
  AntiTrackCuts->SetCutCharge(-1);
  // If  cut variation needed
  if(suffix=="1" || suffix=="8"){
    TrackCuts->SetEtaRange(-0.9, 0.9);
    AntiTrackCuts->SetEtaRange(-0.9, 0.9);
  }
  //Lambda Cuts
  AliFemtoDreamv0Cuts *v0Cuts = AliFemtoDreamv0Cuts::LambdaCuts(isMC, true,
                                                                false);
  AliFemtoDreamTrackCuts *Posv0Daug = AliFemtoDreamTrackCuts::DecayProtonCuts(
      isMC, true, false);

  AliFemtoDreamTrackCuts *Negv0Daug = AliFemtoDreamTrackCuts::DecayPionCuts(
      isMC, true, false);
  // If  cut variation needed
  if(suffix=="2" || suffix=="7" || suffix=="8"){
    Posv0Daug->SetEtaRange(-0.9, 0.9);
    Negv0Daug->SetEtaRange(-0.9, 0.9);
  }
  v0Cuts->SetPosDaugterTrackCuts(Posv0Daug);
  v0Cuts->SetNegDaugterTrackCuts(Negv0Daug);
  v0Cuts->SetPDGCodePosDaug(2212);  //Proton
  v0Cuts->SetPDGCodeNegDaug(211);  //Pion
  v0Cuts->SetPDGCodev0(3122);  //Lambda

  AliFemtoDreamv0Cuts *Antiv0Cuts = AliFemtoDreamv0Cuts::LambdaCuts(isMC, true,
                                                                    false);
  AliFemtoDreamTrackCuts *PosAntiv0Daug = AliFemtoDreamTrackCuts::DecayPionCuts(
      isMC, true, false);
  PosAntiv0Daug->SetCutCharge(1);
  AliFemtoDreamTrackCuts *NegAntiv0Daug =
      AliFemtoDreamTrackCuts::DecayProtonCuts(isMC, true, false);
  NegAntiv0Daug->SetCutCharge(-1);
  // If  cut variation needed
  if(suffix=="2" || suffix=="7" || suffix=="8"){
    PosAntiv0Daug->SetEtaRange(-0.9, 0.9);
    NegAntiv0Daug->SetEtaRange(-0.9, 0.9);
  }
  Antiv0Cuts->SetPosDaugterTrackCuts(PosAntiv0Daug);
  Antiv0Cuts->SetNegDaugterTrackCuts(NegAntiv0Daug);
  Antiv0Cuts->SetPDGCodePosDaug(211);  //Pion
  Antiv0Cuts->SetPDGCodeNegDaug(2212);  //Proton
  Antiv0Cuts->SetPDGCodev0(-3122);  //Lambda

  // If  cut variations eta
  if(suffix=="3" || suffix=="6" || suffix=="7" || suffix=="8"){
    v0Cuts->SetCutInvMass(0.006);
    Antiv0Cuts->SetCutInvMass(0.006);
  }
  if(suffix=="4"|| suffix=="6" || suffix=="7" || suffix=="8"){
    v0Cuts->SetCutCPA(0.999);
    Antiv0Cuts->SetCutCPA(0.999);
  }
  if(suffix=="5"|| suffix=="6" || suffix=="7" || suffix=="8"){
    v0Cuts->SetDaughterTimingCut(AliFemtoDreamv0Cuts::OneDaughterCombined);
    Antiv0Cuts->SetDaughterTimingCut(AliFemtoDreamv0Cuts::OneDaughterCombined);
  }


  if (!fullBlastQA) {
    evtCuts->SetMinimalBooking(true);
    TrackCuts->SetMinimalBooking(true);
    AntiTrackCuts->SetMinimalBooking(true);
    v0Cuts->SetMinimalBooking(true);
    Antiv0Cuts->SetMinimalBooking(true);
  }


  AliFemtoDreamEventCuts *evtCutsTrigger;
  AliFemtoDreamTrackCuts *TrackCutsTrigger;
  AliFemtoDreamTrackCuts *AntiTrackCutsTrigger;
  AliFemtoDreamv0Cuts *v0CutsTrigger;
  AliFemtoDreamTrackCuts *Posv0DaugTrigger;
  AliFemtoDreamTrackCuts *Negv0DaugTrigger;
  AliFemtoDreamv0Cuts *Antiv0CutsTrigger;
  AliFemtoDreamTrackCuts *PosAntiv0DaugTrigger;
  AliFemtoDreamTrackCuts *NegAntiv0DaugTrigger;

  bool TriggerOnSample = false;

  if(triggerOn){
    if(suffixTrigger=="0"){
      Q3Limit = 3.0;
      evtCutsTrigger = new AliFemtoDreamEventCuts();
      evtCutsTrigger->UseDontWorryEvtCuts(false);
      evtCutsTrigger->SetZVtxPosition(-15., 15.);
      //evtCutsTrigger->CleanUpMult(false, false, false, true); will I have this in Run3 for trigger???

      TrackCutsTrigger =  new AliFemtoDreamTrackCuts();
      TrackCutsTrigger->SetPlotDCADist(true);
      TrackCutsTrigger->SetIsMonteCarlo(isMC);
      
      TrackCutsTrigger->SetFilterBit(128); //will I have this in Run3 for trigger???
      TrackCutsTrigger->SetCutCharge(1);
      TrackCutsTrigger->SetPtRange(0.3, 100000.);
      TrackCutsTrigger->SetEtaRange(-0.9, 0.9);
      TrackCutsTrigger->SetCutTPCCrossedRows(true, 60., 0.);
      //TrackCutsTrigger->SetNClsTPC(60);// for now use only crossed pad rows
      TrackCutsTrigger->SetPID(AliPID::kProton, 0.75,5., false, 5.);  
      TrackCutsTrigger->SetDCAReCalculation(true);  //will I have this in Run3 for trigger???
      TrackCutsTrigger->SetDCAVtxZ(0.8);
      TrackCutsTrigger->SetDCAVtxXY(0.4);

      AntiTrackCutsTrigger = new AliFemtoDreamTrackCuts();
      AntiTrackCutsTrigger->SetPlotDCADist(true);
      AntiTrackCutsTrigger->SetIsMonteCarlo(isMC);

      AntiTrackCutsTrigger->SetFilterBit(128); //will I have this in Run3 for trigger???
      AntiTrackCutsTrigger->SetCutCharge(-1);
      AntiTrackCutsTrigger->SetPtRange(0.3, 100000.);
      AntiTrackCutsTrigger->SetEtaRange(-0.9, 0.9);
      AntiTrackCutsTrigger->SetCutTPCCrossedRows(true, 60., 0.);
      //AntiTrackCutsTrigger->SetNClsTPC(60); // for now use only crossed pad rows
      AntiTrackCutsTrigger->SetPID(AliPID::kProton, 0.75,5., false, 5.); 
      AntiTrackCutsTrigger->SetDCAReCalculation(true);  //will I have this in Run3 for trigger???
      AntiTrackCutsTrigger->SetDCAVtxZ(0.8);
      AntiTrackCutsTrigger->SetDCAVtxXY(0.4);

      v0CutsTrigger = new AliFemtoDreamv0Cuts();
      v0CutsTrigger->SetIsMonteCarlo(isMC);
      v0CutsTrigger->SetPlotCPADist(true);

      Posv0DaugTrigger = new AliFemtoDreamTrackCuts(); // proton
      Posv0DaugTrigger->SetIsMonteCarlo(isMC);
      Posv0DaugTrigger->SetFillQALater(true);
      Posv0DaugTrigger->SetCutCharge(1);

      Negv0DaugTrigger = new AliFemtoDreamTrackCuts(); //pion
      Negv0DaugTrigger->SetIsMonteCarlo(isMC);
      Negv0DaugTrigger->SetFillQALater(true);
      Negv0DaugTrigger->SetCutCharge(-1);

      v0CutsTrigger->SetCutCharge(0);
      v0CutsTrigger->SetPtRange(0., 999.);
      v0CutsTrigger->SetCutCPA(0.9);
      v0CutsTrigger->SetCutDCADaugToPrimVtx(0.02);
      v0CutsTrigger->SetCutDCADaugTov0Vtx(3.);
      v0CutsTrigger->SetKaonRejection(0.49, 0.505);
      v0CutsTrigger->SetCutMaxDecayVtx(150);
      v0CutsTrigger->SetCutTransverseRadius(0.1, 150);
      v0CutsTrigger->SetCutInvMass(0.010);
      v0CutsTrigger->SetAxisInvMassPlots(400, 1.0, 1.2);

      Posv0DaugTrigger->SetCutTPCCrossedRows(true, 60, 0.);
      Negv0DaugTrigger->SetCutTPCCrossedRows(true, 60, 0.);
      Posv0DaugTrigger->SetEtaRange(-0.9, 0.9);
      Negv0DaugTrigger->SetEtaRange(-0.9, 0.9);
      Posv0DaugTrigger->SetPID(AliPID::kProton, 999., 7.);
      Negv0DaugTrigger->SetPID(AliPID::kPion, 999., 7.);
      Negv0DaugTrigger->SetDCAReCalculation(true);
      Posv0DaugTrigger->SetDCAReCalculation(true);

      v0CutsTrigger->SetPosDaugterTrackCuts(Posv0DaugTrigger);
      v0CutsTrigger->SetNegDaugterTrackCuts(Negv0DaugTrigger);
      v0CutsTrigger->SetPDGCodePosDaug(2212);  //Proton
      v0CutsTrigger->SetPDGCodeNegDaug(211);  //Pion
      v0CutsTrigger->SetPDGCodev0(3122);  //Lambda

      Antiv0CutsTrigger = new AliFemtoDreamv0Cuts();
      Antiv0CutsTrigger->SetIsMonteCarlo(isMC);
      Antiv0CutsTrigger->SetPlotCPADist(true);
      Antiv0CutsTrigger->SetPlotContrib(false);

      PosAntiv0DaugTrigger = new AliFemtoDreamTrackCuts();
      PosAntiv0DaugTrigger->SetIsMonteCarlo(isMC);
      PosAntiv0DaugTrigger->SetFillQALater(true);
      PosAntiv0DaugTrigger->SetCutCharge(1);
      NegAntiv0DaugTrigger =new AliFemtoDreamTrackCuts();
      NegAntiv0DaugTrigger->SetIsMonteCarlo(isMC);
      NegAntiv0DaugTrigger->SetFillQALater(true);
      NegAntiv0DaugTrigger->SetCutCharge(-1);

      Antiv0CutsTrigger->SetCutCharge(0);
      Antiv0CutsTrigger->SetPtRange(0., 999.);
      Antiv0CutsTrigger->SetCutCPA(0.9);
      Antiv0CutsTrigger->SetCutDCADaugToPrimVtx(0.02);
      Antiv0CutsTrigger->SetCutDCADaugTov0Vtx(3.);
      Antiv0CutsTrigger->SetKaonRejection(0.49, 0.505);
      Antiv0CutsTrigger->SetCutMaxDecayVtx(150);
      Antiv0CutsTrigger->SetCutTransverseRadius(0.1, 150);
      Antiv0CutsTrigger->SetCutInvMass(0.010);
      Antiv0CutsTrigger->SetAxisInvMassPlots(400, 1.0, 1.2);

      PosAntiv0DaugTrigger->SetCutTPCCrossedRows(true, 60, 0.);
      NegAntiv0DaugTrigger->SetCutTPCCrossedRows(true, 60, 0.);
      PosAntiv0DaugTrigger->SetEtaRange(-0.9, 0.9);
      NegAntiv0DaugTrigger->SetEtaRange(-0.9, 0.9);
      PosAntiv0DaugTrigger->SetPID(AliPID::kPion, 999., 7.);
      NegAntiv0DaugTrigger->SetPID(AliPID::kProton, 999., 7.);
      PosAntiv0DaugTrigger->SetDCAReCalculation(true);
      NegAntiv0DaugTrigger->SetDCAReCalculation(true);

      Antiv0CutsTrigger->SetPosDaugterTrackCuts(PosAntiv0DaugTrigger);
      Antiv0CutsTrigger->SetNegDaugterTrackCuts(NegAntiv0DaugTrigger);
      Antiv0CutsTrigger->SetPDGCodePosDaug(211);  //Pion
      Antiv0CutsTrigger->SetPDGCodeNegDaug(2212);  //Proton
      Antiv0CutsTrigger->SetPDGCodev0(-3122);  //Lambda
      if (!fullBlastQA) {
        evtCutsTrigger->SetMinimalBooking(true);
        TrackCutsTrigger->SetMinimalBooking(true);
        AntiTrackCutsTrigger->SetMinimalBooking(true);
        v0CutsTrigger->SetMinimalBooking(true);
        Antiv0CutsTrigger->SetMinimalBooking(true);
      }
    }

    if(suffixTrigger=="1"){
      Q3Limit = 1.5;
      evtCutsTrigger = new AliFemtoDreamEventCuts();
      evtCutsTrigger->UseDontWorryEvtCuts(false);
      evtCutsTrigger->SetZVtxPosition(-15., 15.);
      //evtCutsTrigger->CleanUpMult(false, false, false, true); will I have this in Run3 for trigger???

      TrackCutsTrigger =  new AliFemtoDreamTrackCuts();
      TrackCutsTrigger->SetPlotDCADist(true);
      TrackCutsTrigger->SetIsMonteCarlo(isMC);
      
      TrackCutsTrigger->SetFilterBit(128); //will I have this in Run3 for trigger???
      TrackCutsTrigger->SetCutCharge(1);
      TrackCutsTrigger->SetPtRange(0.3, 100000.);
      TrackCutsTrigger->SetEtaRange(-0.9, 0.9);
      TrackCutsTrigger->SetCutTPCCrossedRows(true, 60., 0.);
      //TrackCutsTrigger->SetNClsTPC(60);// for now use only crossed pad rows
      TrackCutsTrigger->SetPID(AliPID::kProton, 0.75,5., false, 5.);  
      TrackCutsTrigger->SetDCAReCalculation(true);  //will I have this in Run3 for trigger???
      TrackCutsTrigger->SetDCAVtxZ(0.8);
      TrackCutsTrigger->SetDCAVtxXY(0.4);

      AntiTrackCutsTrigger = new AliFemtoDreamTrackCuts();
      AntiTrackCutsTrigger->SetPlotDCADist(true);
      AntiTrackCutsTrigger->SetIsMonteCarlo(isMC);

      AntiTrackCutsTrigger->SetFilterBit(128); //will I have this in Run3 for trigger???
      AntiTrackCutsTrigger->SetCutCharge(-1);
      AntiTrackCutsTrigger->SetPtRange(0.3, 100000.);
      AntiTrackCutsTrigger->SetEtaRange(-0.9, 0.9);
      AntiTrackCutsTrigger->SetCutTPCCrossedRows(true, 60., 0.);
      //AntiTrackCutsTrigger->SetNClsTPC(60); // for now use only crossed pad rows
      AntiTrackCutsTrigger->SetPID(AliPID::kProton, 0.75,5., false, 5.); 
      AntiTrackCutsTrigger->SetDCAReCalculation(true);  //will I have this in Run3 for trigger???
      AntiTrackCutsTrigger->SetDCAVtxZ(0.8);
      AntiTrackCutsTrigger->SetDCAVtxXY(0.4);

      v0CutsTrigger = new AliFemtoDreamv0Cuts();
      v0CutsTrigger->SetIsMonteCarlo(isMC);
      v0CutsTrigger->SetPlotCPADist(true);

      Posv0DaugTrigger = new AliFemtoDreamTrackCuts(); // proton
      Posv0DaugTrigger->SetIsMonteCarlo(isMC);
      Posv0DaugTrigger->SetFillQALater(true);
      Posv0DaugTrigger->SetCutCharge(1);

      Negv0DaugTrigger = new AliFemtoDreamTrackCuts(); //pion
      Negv0DaugTrigger->SetIsMonteCarlo(isMC);
      Negv0DaugTrigger->SetFillQALater(true);
      Negv0DaugTrigger->SetCutCharge(-1);

      v0CutsTrigger->SetCutCharge(0);
      v0CutsTrigger->SetPtRange(0., 999.);
      v0CutsTrigger->SetCutCPA(0.9);
      v0CutsTrigger->SetCutDCADaugToPrimVtx(0.02);
      v0CutsTrigger->SetCutDCADaugTov0Vtx(3.);
      v0CutsTrigger->SetKaonRejection(0.49, 0.505);
      v0CutsTrigger->SetCutMaxDecayVtx(150);
      v0CutsTrigger->SetCutTransverseRadius(0.1, 150);
      v0CutsTrigger->SetCutInvMass(0.010);
      v0CutsTrigger->SetAxisInvMassPlots(400, 1.0, 1.2);

      Posv0DaugTrigger->SetCutTPCCrossedRows(true, 60, 0.);
      Negv0DaugTrigger->SetCutTPCCrossedRows(true, 60, 0.);
      Posv0DaugTrigger->SetEtaRange(-0.9, 0.9);
      Negv0DaugTrigger->SetEtaRange(-0.9, 0.9);
      Posv0DaugTrigger->SetPID(AliPID::kProton, 999., 7.);
      Negv0DaugTrigger->SetPID(AliPID::kPion, 999., 7.);
      Negv0DaugTrigger->SetDCAReCalculation(true);
      Posv0DaugTrigger->SetDCAReCalculation(true);

      v0CutsTrigger->SetPosDaugterTrackCuts(Posv0DaugTrigger);
      v0CutsTrigger->SetNegDaugterTrackCuts(Negv0DaugTrigger);
      v0CutsTrigger->SetPDGCodePosDaug(2212);  //Proton
      v0CutsTrigger->SetPDGCodeNegDaug(211);  //Pion
      v0CutsTrigger->SetPDGCodev0(3122);  //Lambda

      Antiv0CutsTrigger = new AliFemtoDreamv0Cuts();
      Antiv0CutsTrigger->SetIsMonteCarlo(isMC);
      Antiv0CutsTrigger->SetPlotCPADist(true);
      Antiv0CutsTrigger->SetPlotContrib(false);

      PosAntiv0DaugTrigger = new AliFemtoDreamTrackCuts();
      PosAntiv0DaugTrigger->SetIsMonteCarlo(isMC);
      PosAntiv0DaugTrigger->SetFillQALater(true);
      PosAntiv0DaugTrigger->SetCutCharge(1);
      NegAntiv0DaugTrigger =new AliFemtoDreamTrackCuts();
      NegAntiv0DaugTrigger->SetIsMonteCarlo(isMC);
      NegAntiv0DaugTrigger->SetFillQALater(true);
      NegAntiv0DaugTrigger->SetCutCharge(-1);

      Antiv0CutsTrigger->SetCutCharge(0);
      Antiv0CutsTrigger->SetPtRange(0., 999.);
      Antiv0CutsTrigger->SetCutCPA(0.9);
      Antiv0CutsTrigger->SetCutDCADaugToPrimVtx(0.02);
      Antiv0CutsTrigger->SetCutDCADaugTov0Vtx(3.);
      Antiv0CutsTrigger->SetKaonRejection(0.49, 0.505);
      Antiv0CutsTrigger->SetCutMaxDecayVtx(150);
      Antiv0CutsTrigger->SetCutTransverseRadius(0.1, 150);
      Antiv0CutsTrigger->SetCutInvMass(0.010);
      Antiv0CutsTrigger->SetAxisInvMassPlots(400, 1.0, 1.2);

      PosAntiv0DaugTrigger->SetCutTPCCrossedRows(true, 60, 0.);
      NegAntiv0DaugTrigger->SetCutTPCCrossedRows(true, 60, 0.);
      PosAntiv0DaugTrigger->SetEtaRange(-0.9, 0.9);
      NegAntiv0DaugTrigger->SetEtaRange(-0.9, 0.9);
      PosAntiv0DaugTrigger->SetPID(AliPID::kPion, 999., 7.);
      NegAntiv0DaugTrigger->SetPID(AliPID::kProton, 999., 7.);
      PosAntiv0DaugTrigger->SetDCAReCalculation(true);
      NegAntiv0DaugTrigger->SetDCAReCalculation(true);

      Antiv0CutsTrigger->SetPosDaugterTrackCuts(PosAntiv0DaugTrigger);
      Antiv0CutsTrigger->SetNegDaugterTrackCuts(NegAntiv0DaugTrigger);
      Antiv0CutsTrigger->SetPDGCodePosDaug(211);  //Pion
      Antiv0CutsTrigger->SetPDGCodeNegDaug(2212);  //Proton
      Antiv0CutsTrigger->SetPDGCodev0(-3122);  //Lambda
      if (!fullBlastQA) {
        evtCutsTrigger->SetMinimalBooking(true);
        TrackCutsTrigger->SetMinimalBooking(true);
        AntiTrackCutsTrigger->SetMinimalBooking(true);
        v0CutsTrigger->SetMinimalBooking(true);
        Antiv0CutsTrigger->SetMinimalBooking(true);
      }
    }


    if(suffixTrigger=="2"){
      Q3Limit = 1.5;
      evtCutsTrigger = new AliFemtoDreamEventCuts();
      evtCutsTrigger->UseDontWorryEvtCuts(false);
      evtCutsTrigger->SetZVtxPosition(-15., 15.);
      //evtCutsTrigger->CleanUpMult(false, false, false, true); will I have this in Run3 for trigger???

      TrackCutsTrigger =  new AliFemtoDreamTrackCuts();
      TrackCutsTrigger->SetPlotDCADist(true);
      TrackCutsTrigger->SetIsMonteCarlo(isMC);
      
      TrackCutsTrigger->SetFilterBit(128); //will I have this in Run3 for trigger???
      TrackCutsTrigger->SetCutCharge(1);
      TrackCutsTrigger->SetPtRange(0.3, 5.);
      TrackCutsTrigger->SetEtaRange(-0.9, 0.9);
      TrackCutsTrigger->SetCutTPCCrossedRows(true, 60., 0.);
      //TrackCutsTrigger->SetNClsTPC(60);// for now use only crossed pad rows
      TrackCutsTrigger->SetPID(AliPID::kProton, 0.75,4., false, 4.);  
      TrackCutsTrigger->SetDCAReCalculation(true);  //will I have this in Run3 for trigger???
      TrackCutsTrigger->SetDCAVtxZ(0.8);
      TrackCutsTrigger->SetDCAVtxXY(0.4);

      AntiTrackCutsTrigger = new AliFemtoDreamTrackCuts();
      AntiTrackCutsTrigger->SetPlotDCADist(true);
      AntiTrackCutsTrigger->SetIsMonteCarlo(isMC);

      AntiTrackCutsTrigger->SetFilterBit(128); //will I have this in Run3 for trigger???
      AntiTrackCutsTrigger->SetCutCharge(-1);
      AntiTrackCutsTrigger->SetPtRange(0.3, 5.);
      AntiTrackCutsTrigger->SetEtaRange(-0.9, 0.9);
      AntiTrackCutsTrigger->SetCutTPCCrossedRows(true, 60., 0.);
      //AntiTrackCutsTrigger->SetNClsTPC(60); // for now use only crossed pad rows
      AntiTrackCutsTrigger->SetPID(AliPID::kProton, 0.75,4., false, 4.); 
      AntiTrackCutsTrigger->SetDCAReCalculation(true);  //will I have this in Run3 for trigger???
      AntiTrackCutsTrigger->SetDCAVtxZ(0.8);
      AntiTrackCutsTrigger->SetDCAVtxXY(0.4);

      v0CutsTrigger = new AliFemtoDreamv0Cuts();
      v0CutsTrigger->SetIsMonteCarlo(isMC);
      v0CutsTrigger->SetPlotCPADist(true);

      Posv0DaugTrigger = new AliFemtoDreamTrackCuts(); // proton
      Posv0DaugTrigger->SetIsMonteCarlo(isMC);
      Posv0DaugTrigger->SetFillQALater(true);
      Posv0DaugTrigger->SetCutCharge(1);

      Negv0DaugTrigger = new AliFemtoDreamTrackCuts(); //pion
      Negv0DaugTrigger->SetIsMonteCarlo(isMC);
      Negv0DaugTrigger->SetFillQALater(true);
      Negv0DaugTrigger->SetCutCharge(-1);

      v0CutsTrigger->SetCutCharge(0);
      v0CutsTrigger->SetPtRange(0., 999.);
      v0CutsTrigger->SetCutCPA(0.9);
      v0CutsTrigger->SetCutDCADaugToPrimVtx(0.02);
      v0CutsTrigger->SetCutDCADaugTov0Vtx(3.);
      v0CutsTrigger->SetKaonRejection(0.49, 0.505);
      v0CutsTrigger->SetCutMaxDecayVtx(120);
      v0CutsTrigger->SetCutTransverseRadius(0.1, 120);
      v0CutsTrigger->SetCutInvMass(0.006);
      v0CutsTrigger->SetAxisInvMassPlots(400, 1.0, 1.2);

      Posv0DaugTrigger->SetCutTPCCrossedRows(true, 60, 0.);
      Negv0DaugTrigger->SetCutTPCCrossedRows(true, 60, 0.);
      Posv0DaugTrigger->SetEtaRange(-0.9, 0.9);
      Negv0DaugTrigger->SetEtaRange(-0.9, 0.9);
      Posv0DaugTrigger->SetPID(AliPID::kProton, 999., 7.);
      Negv0DaugTrigger->SetPID(AliPID::kPion, 999., 7.);
      Negv0DaugTrigger->SetDCAReCalculation(true);
      Posv0DaugTrigger->SetDCAReCalculation(true);

      v0CutsTrigger->SetPosDaugterTrackCuts(Posv0DaugTrigger);
      v0CutsTrigger->SetNegDaugterTrackCuts(Negv0DaugTrigger);
      v0CutsTrigger->SetPDGCodePosDaug(2212);  //Proton
      v0CutsTrigger->SetPDGCodeNegDaug(211);  //Pion
      v0CutsTrigger->SetPDGCodev0(3122);  //Lambda

      Antiv0CutsTrigger = new AliFemtoDreamv0Cuts();
      Antiv0CutsTrigger->SetIsMonteCarlo(isMC);
      Antiv0CutsTrigger->SetPlotCPADist(true);
      Antiv0CutsTrigger->SetPlotContrib(false);

      PosAntiv0DaugTrigger = new AliFemtoDreamTrackCuts();
      PosAntiv0DaugTrigger->SetIsMonteCarlo(isMC);
      PosAntiv0DaugTrigger->SetFillQALater(true);
      PosAntiv0DaugTrigger->SetCutCharge(1);
      NegAntiv0DaugTrigger =new AliFemtoDreamTrackCuts();
      NegAntiv0DaugTrigger->SetIsMonteCarlo(isMC);
      NegAntiv0DaugTrigger->SetFillQALater(true);
      NegAntiv0DaugTrigger->SetCutCharge(-1);

      Antiv0CutsTrigger->SetCutCharge(0);
      Antiv0CutsTrigger->SetPtRange(0., 999.);
      Antiv0CutsTrigger->SetCutCPA(0.9);
      Antiv0CutsTrigger->SetCutDCADaugToPrimVtx(0.02);
      Antiv0CutsTrigger->SetCutDCADaugTov0Vtx(3.);
      Antiv0CutsTrigger->SetKaonRejection(0.49, 0.505);
      Antiv0CutsTrigger->SetCutMaxDecayVtx(120);
      Antiv0CutsTrigger->SetCutTransverseRadius(0.1, 120);
      Antiv0CutsTrigger->SetCutInvMass(0.006);
      Antiv0CutsTrigger->SetAxisInvMassPlots(400, 1.0, 1.2);

      PosAntiv0DaugTrigger->SetCutTPCCrossedRows(true, 60, 0.);
      NegAntiv0DaugTrigger->SetCutTPCCrossedRows(true, 60, 0.);
      PosAntiv0DaugTrigger->SetEtaRange(-0.9, 0.9);
      NegAntiv0DaugTrigger->SetEtaRange(-0.9, 0.9);
      PosAntiv0DaugTrigger->SetPID(AliPID::kPion, 999., 7.);
      NegAntiv0DaugTrigger->SetPID(AliPID::kProton, 999., 7.);
      PosAntiv0DaugTrigger->SetDCAReCalculation(true);
      NegAntiv0DaugTrigger->SetDCAReCalculation(true);

      Antiv0CutsTrigger->SetPosDaugterTrackCuts(PosAntiv0DaugTrigger);
      Antiv0CutsTrigger->SetNegDaugterTrackCuts(NegAntiv0DaugTrigger);
      Antiv0CutsTrigger->SetPDGCodePosDaug(211);  //Pion
      Antiv0CutsTrigger->SetPDGCodeNegDaug(2212);  //Proton
      Antiv0CutsTrigger->SetPDGCodev0(-3122);  //Lambda
      if (!fullBlastQA) {
        evtCutsTrigger->SetMinimalBooking(true);
        TrackCutsTrigger->SetMinimalBooking(true);
        AntiTrackCutsTrigger->SetMinimalBooking(true);
        v0CutsTrigger->SetMinimalBooking(true);
        Antiv0CutsTrigger->SetMinimalBooking(true);
      }
    }


    if(suffixTrigger=="3"){
      Q3Limit = 1.5;
      evtCutsTrigger = new AliFemtoDreamEventCuts();
      evtCutsTrigger->UseDontWorryEvtCuts(false);
      evtCutsTrigger->SetZVtxPosition(-15., 15.);
      //evtCutsTrigger->CleanUpMult(false, false, false, true); will I have this in Run3 for trigger???

      TrackCutsTrigger =  new AliFemtoDreamTrackCuts();
      TrackCutsTrigger->SetPlotDCADist(true);
      TrackCutsTrigger->SetIsMonteCarlo(isMC);
      
      TrackCutsTrigger->SetFilterBit(128); //will I have this in Run3 for trigger???
      TrackCutsTrigger->SetCutCharge(1);
      TrackCutsTrigger->SetPtRange(0.3, 5.);
      TrackCutsTrigger->SetEtaRange(-0.9, 0.9);
      TrackCutsTrigger->SetCutTPCCrossedRows(true, 60., 0.);
      //TrackCutsTrigger->SetNClsTPC(60);// for now use only crossed pad rows
      TrackCutsTrigger->SetPID(AliPID::kProton, 0.75,4., false, 4.);  
      TrackCutsTrigger->SetDCAReCalculation(true);  //will I have this in Run3 for trigger???
      TrackCutsTrigger->SetDCAVtxZ(0.8);
      TrackCutsTrigger->SetDCAVtxXY(0.4);

      AntiTrackCutsTrigger = new AliFemtoDreamTrackCuts();
      AntiTrackCutsTrigger->SetPlotDCADist(true);
      AntiTrackCutsTrigger->SetIsMonteCarlo(isMC);

      AntiTrackCutsTrigger->SetFilterBit(128); //will I have this in Run3 for trigger???
      AntiTrackCutsTrigger->SetCutCharge(-1);
      AntiTrackCutsTrigger->SetPtRange(0.3, 5.);
      AntiTrackCutsTrigger->SetEtaRange(-0.9, 0.9);
      AntiTrackCutsTrigger->SetCutTPCCrossedRows(true, 60., 0.);
      //AntiTrackCutsTrigger->SetNClsTPC(60); // for now use only crossed pad rows
      AntiTrackCutsTrigger->SetPID(AliPID::kProton, 0.75,4., false, 4.); 
      AntiTrackCutsTrigger->SetDCAReCalculation(true);  //will I have this in Run3 for trigger???
      AntiTrackCutsTrigger->SetDCAVtxZ(0.8);
      AntiTrackCutsTrigger->SetDCAVtxXY(0.4);

      v0CutsTrigger = new AliFemtoDreamv0Cuts();
      v0CutsTrigger->SetIsMonteCarlo(isMC);
      v0CutsTrigger->SetPlotCPADist(true);

      Posv0DaugTrigger = new AliFemtoDreamTrackCuts(); // proton
      Posv0DaugTrigger->SetIsMonteCarlo(isMC);
      Posv0DaugTrigger->SetFillQALater(true);
      Posv0DaugTrigger->SetCutCharge(1);
      Posv0DaugTrigger->SetCheckPileUp(true);

      Negv0DaugTrigger = new AliFemtoDreamTrackCuts(); //pion
      Negv0DaugTrigger->SetIsMonteCarlo(isMC);
      Negv0DaugTrigger->SetFillQALater(true);
      Negv0DaugTrigger->SetCutCharge(-1);
      Negv0DaugTrigger->SetCheckPileUp(true);

      v0CutsTrigger->SetCutCharge(0);
      v0CutsTrigger->SetPtRange(0., 999.);
      v0CutsTrigger->SetCutCPA(0.9);
      v0CutsTrigger->SetCutDCADaugToPrimVtx(0.02);
      v0CutsTrigger->SetCutDCADaugTov0Vtx(3.);
      v0CutsTrigger->SetKaonRejection(0.49, 0.505);
      v0CutsTrigger->SetCutMaxDecayVtx(120);
      v0CutsTrigger->SetCutTransverseRadius(0.1, 120);
      v0CutsTrigger->SetCutInvMass(0.006);
      v0CutsTrigger->SetAxisInvMassPlots(400, 1.0, 1.2);

      Posv0DaugTrigger->SetCutTPCCrossedRows(true, 60, 0.);
      Negv0DaugTrigger->SetCutTPCCrossedRows(true, 60, 0.);
      Posv0DaugTrigger->SetEtaRange(-0.9, 0.9);
      Negv0DaugTrigger->SetEtaRange(-0.9, 0.9);
      Posv0DaugTrigger->SetPID(AliPID::kProton, 999., 7.);
      Negv0DaugTrigger->SetPID(AliPID::kPion, 999., 7.);
      Negv0DaugTrigger->SetDCAReCalculation(true);
      Posv0DaugTrigger->SetDCAReCalculation(true);

      v0CutsTrigger->SetPosDaugterTrackCuts(Posv0DaugTrigger);
      v0CutsTrigger->SetNegDaugterTrackCuts(Negv0DaugTrigger);
      v0CutsTrigger->SetPDGCodePosDaug(2212);  //Proton
      v0CutsTrigger->SetPDGCodeNegDaug(211);  //Pion
      v0CutsTrigger->SetPDGCodev0(3122);  //Lambda

      Antiv0CutsTrigger = new AliFemtoDreamv0Cuts();
      Antiv0CutsTrigger->SetIsMonteCarlo(isMC);
      Antiv0CutsTrigger->SetPlotCPADist(true);
      Antiv0CutsTrigger->SetPlotContrib(false);

      PosAntiv0DaugTrigger = new AliFemtoDreamTrackCuts();
      PosAntiv0DaugTrigger->SetIsMonteCarlo(isMC);
      PosAntiv0DaugTrigger->SetFillQALater(true);
      PosAntiv0DaugTrigger->SetCutCharge(1);
      PosAntiv0DaugTrigger->SetCheckPileUp(true);

      NegAntiv0DaugTrigger =new AliFemtoDreamTrackCuts();
      NegAntiv0DaugTrigger->SetIsMonteCarlo(isMC);
      NegAntiv0DaugTrigger->SetFillQALater(true);
      NegAntiv0DaugTrigger->SetCutCharge(-1);
      NegAntiv0DaugTrigger->SetCheckPileUp(true);

      Antiv0CutsTrigger->SetCutCharge(0);
      Antiv0CutsTrigger->SetPtRange(0., 999.);
      Antiv0CutsTrigger->SetCutCPA(0.9);
      Antiv0CutsTrigger->SetCutDCADaugToPrimVtx(0.02);
      Antiv0CutsTrigger->SetCutDCADaugTov0Vtx(3.);
      Antiv0CutsTrigger->SetKaonRejection(0.49, 0.505);
      Antiv0CutsTrigger->SetCutMaxDecayVtx(120);
      Antiv0CutsTrigger->SetCutTransverseRadius(0.1, 120);
      Antiv0CutsTrigger->SetCutInvMass(0.006);
      Antiv0CutsTrigger->SetAxisInvMassPlots(400, 1.0, 1.2);

      PosAntiv0DaugTrigger->SetCutTPCCrossedRows(true, 60, 0.);
      NegAntiv0DaugTrigger->SetCutTPCCrossedRows(true, 60, 0.);
      PosAntiv0DaugTrigger->SetEtaRange(-0.9, 0.9);
      NegAntiv0DaugTrigger->SetEtaRange(-0.9, 0.9);
      PosAntiv0DaugTrigger->SetPID(AliPID::kPion, 999., 7.);
      NegAntiv0DaugTrigger->SetPID(AliPID::kProton, 999., 7.);
      PosAntiv0DaugTrigger->SetDCAReCalculation(true);
      NegAntiv0DaugTrigger->SetDCAReCalculation(true);

      Antiv0CutsTrigger->SetPosDaugterTrackCuts(PosAntiv0DaugTrigger);
      Antiv0CutsTrigger->SetNegDaugterTrackCuts(NegAntiv0DaugTrigger);
      Antiv0CutsTrigger->SetPDGCodePosDaug(211);  //Pion
      Antiv0CutsTrigger->SetPDGCodeNegDaug(2212);  //Proton
      Antiv0CutsTrigger->SetPDGCodev0(-3122);  //Lambda
      if (!fullBlastQA) {
        evtCutsTrigger->SetMinimalBooking(true);
        TrackCutsTrigger->SetMinimalBooking(true);
        AntiTrackCutsTrigger->SetMinimalBooking(true);
        v0CutsTrigger->SetMinimalBooking(true);
        Antiv0CutsTrigger->SetMinimalBooking(true);
      }
    }

    if(suffixTrigger=="4"){
      Q3Limit = 3.;
      evtCutsTrigger = new AliFemtoDreamEventCuts();
      evtCutsTrigger->UseDontWorryEvtCuts(false);
      evtCutsTrigger->SetZVtxPosition(-15., 15.);
      //evtCutsTrigger->CleanUpMult(false, false, false, true); will I have this in Run3 for trigger???

      TrackCutsTrigger =  new AliFemtoDreamTrackCuts();
      TrackCutsTrigger->SetPlotDCADist(true);
      TrackCutsTrigger->SetIsMonteCarlo(isMC);
      
      TrackCutsTrigger->SetFilterBit(128); //will I have this in Run3 for trigger???
      TrackCutsTrigger->SetCutCharge(1);
      TrackCutsTrigger->SetPtRange(0.3, 5.);
      TrackCutsTrigger->SetEtaRange(-0.9, 0.9);
      TrackCutsTrigger->SetCutTPCCrossedRows(true, 60., 0.);
      //TrackCutsTrigger->SetNClsTPC(60);// for now use only crossed pad rows
      TrackCutsTrigger->SetPID(AliPID::kProton, 0.75,4., false, 4.);  
      TrackCutsTrigger->SetDCAReCalculation(true);  //will I have this in Run3 for trigger???
      TrackCutsTrigger->SetDCAVtxZ(0.8);
      TrackCutsTrigger->SetDCAVtxXY(0.4);

      AntiTrackCutsTrigger = new AliFemtoDreamTrackCuts();
      AntiTrackCutsTrigger->SetPlotDCADist(true);
      AntiTrackCutsTrigger->SetIsMonteCarlo(isMC);

      AntiTrackCutsTrigger->SetFilterBit(128); //will I have this in Run3 for trigger???
      AntiTrackCutsTrigger->SetCutCharge(-1);
      AntiTrackCutsTrigger->SetPtRange(0.3, 5.);
      AntiTrackCutsTrigger->SetEtaRange(-0.9, 0.9);
      AntiTrackCutsTrigger->SetCutTPCCrossedRows(true, 60., 0.);
      //AntiTrackCutsTrigger->SetNClsTPC(60); // for now use only crossed pad rows
      AntiTrackCutsTrigger->SetPID(AliPID::kProton, 0.75,4., false, 4.); 
      AntiTrackCutsTrigger->SetDCAReCalculation(true);  //will I have this in Run3 for trigger???
      AntiTrackCutsTrigger->SetDCAVtxZ(0.8);
      AntiTrackCutsTrigger->SetDCAVtxXY(0.4);

      v0CutsTrigger = new AliFemtoDreamv0Cuts();
      v0CutsTrigger->SetIsMonteCarlo(isMC);
      v0CutsTrigger->SetPlotCPADist(true);

      Posv0DaugTrigger = new AliFemtoDreamTrackCuts(); // proton
      Posv0DaugTrigger->SetIsMonteCarlo(isMC);
      Posv0DaugTrigger->SetFillQALater(true);
      Posv0DaugTrigger->SetCutCharge(1);
      Posv0DaugTrigger->SetCheckPileUp(true);

      Negv0DaugTrigger = new AliFemtoDreamTrackCuts(); //pion
      Negv0DaugTrigger->SetIsMonteCarlo(isMC);
      Negv0DaugTrigger->SetFillQALater(true);
      Negv0DaugTrigger->SetCutCharge(-1);
      Negv0DaugTrigger->SetCheckPileUp(true);

      v0CutsTrigger->SetCutCharge(0);
      v0CutsTrigger->SetPtRange(0., 999.);
      v0CutsTrigger->SetCutCPA(0.9);
      v0CutsTrigger->SetCutDCADaugToPrimVtx(0.02);
      v0CutsTrigger->SetCutDCADaugTov0Vtx(3.);
      v0CutsTrigger->SetKaonRejection(0.49, 0.505);
      v0CutsTrigger->SetCutMaxDecayVtx(120);
      v0CutsTrigger->SetCutTransverseRadius(0.1, 120);
      v0CutsTrigger->SetCutInvMass(0.006);
      v0CutsTrigger->SetAxisInvMassPlots(400, 1.0, 1.2);

      Posv0DaugTrigger->SetCutTPCCrossedRows(true, 60, 0.);
      Negv0DaugTrigger->SetCutTPCCrossedRows(true, 60, 0.);
      Posv0DaugTrigger->SetEtaRange(-0.9, 0.9);
      Negv0DaugTrigger->SetEtaRange(-0.9, 0.9);
      Posv0DaugTrigger->SetPID(AliPID::kProton, 999., 7.);
      Negv0DaugTrigger->SetPID(AliPID::kPion, 999., 7.);
      Negv0DaugTrigger->SetDCAReCalculation(true);
      Posv0DaugTrigger->SetDCAReCalculation(true);

      v0CutsTrigger->SetPosDaugterTrackCuts(Posv0DaugTrigger);
      v0CutsTrigger->SetNegDaugterTrackCuts(Negv0DaugTrigger);
      v0CutsTrigger->SetPDGCodePosDaug(2212);  //Proton
      v0CutsTrigger->SetPDGCodeNegDaug(211);  //Pion
      v0CutsTrigger->SetPDGCodev0(3122);  //Lambda

      Antiv0CutsTrigger = new AliFemtoDreamv0Cuts();
      Antiv0CutsTrigger->SetIsMonteCarlo(isMC);
      Antiv0CutsTrigger->SetPlotCPADist(true);
      Antiv0CutsTrigger->SetPlotContrib(false);

      PosAntiv0DaugTrigger = new AliFemtoDreamTrackCuts();
      PosAntiv0DaugTrigger->SetIsMonteCarlo(isMC);
      PosAntiv0DaugTrigger->SetFillQALater(true);
      PosAntiv0DaugTrigger->SetCutCharge(1);
      PosAntiv0DaugTrigger->SetCheckPileUp(true);

      NegAntiv0DaugTrigger =new AliFemtoDreamTrackCuts();
      NegAntiv0DaugTrigger->SetIsMonteCarlo(isMC);
      NegAntiv0DaugTrigger->SetFillQALater(true);
      NegAntiv0DaugTrigger->SetCutCharge(-1);
      NegAntiv0DaugTrigger->SetCheckPileUp(true);

      Antiv0CutsTrigger->SetCutCharge(0);
      Antiv0CutsTrigger->SetPtRange(0., 999.);
      Antiv0CutsTrigger->SetCutCPA(0.9);
      Antiv0CutsTrigger->SetCutDCADaugToPrimVtx(0.02);
      Antiv0CutsTrigger->SetCutDCADaugTov0Vtx(3.);
      Antiv0CutsTrigger->SetKaonRejection(0.49, 0.505);
      Antiv0CutsTrigger->SetCutMaxDecayVtx(120);
      Antiv0CutsTrigger->SetCutTransverseRadius(0.1, 120);
      Antiv0CutsTrigger->SetCutInvMass(0.006);
      Antiv0CutsTrigger->SetAxisInvMassPlots(400, 1.0, 1.2);

      PosAntiv0DaugTrigger->SetCutTPCCrossedRows(true, 60, 0.);
      NegAntiv0DaugTrigger->SetCutTPCCrossedRows(true, 60, 0.);
      PosAntiv0DaugTrigger->SetEtaRange(-0.9, 0.9);
      NegAntiv0DaugTrigger->SetEtaRange(-0.9, 0.9);
      PosAntiv0DaugTrigger->SetPID(AliPID::kPion, 999., 7.);
      NegAntiv0DaugTrigger->SetPID(AliPID::kProton, 999., 7.);
      PosAntiv0DaugTrigger->SetDCAReCalculation(true);
      NegAntiv0DaugTrigger->SetDCAReCalculation(true);

      Antiv0CutsTrigger->SetPosDaugterTrackCuts(PosAntiv0DaugTrigger);
      Antiv0CutsTrigger->SetNegDaugterTrackCuts(NegAntiv0DaugTrigger);
      Antiv0CutsTrigger->SetPDGCodePosDaug(211);  //Pion
      Antiv0CutsTrigger->SetPDGCodeNegDaug(2212);  //Proton
      Antiv0CutsTrigger->SetPDGCodev0(-3122);  //Lambda
      if (!fullBlastQA) {
        evtCutsTrigger->SetMinimalBooking(true);
        TrackCutsTrigger->SetMinimalBooking(true);
        AntiTrackCutsTrigger->SetMinimalBooking(true);
        v0CutsTrigger->SetMinimalBooking(true);
        Antiv0CutsTrigger->SetMinimalBooking(true);
      }
    }

if(suffixTrigger=="5"){
      Q3Limit = 3.;
      evtCutsTrigger = new AliFemtoDreamEventCuts();
      evtCutsTrigger->UseDontWorryEvtCuts(false);
      evtCutsTrigger->SetZVtxPosition(-15., 15.);
      //evtCutsTrigger->CleanUpMult(false, false, false, true); will I have this in Run3 for trigger???

      TrackCutsTrigger =  new AliFemtoDreamTrackCuts();
      TrackCutsTrigger->SetPlotDCADist(true);
      TrackCutsTrigger->SetIsMonteCarlo(isMC);
      
      TrackCutsTrigger->SetFilterBit(128); //will I have this in Run3 for trigger???
      TrackCutsTrigger->SetCutCharge(1);
      TrackCutsTrigger->SetPtRange(0.3, 5.);
      TrackCutsTrigger->SetEtaRange(-0.9, 0.9);
      TrackCutsTrigger->SetCutTPCCrossedRows(true, 60., 0.);
      //TrackCutsTrigger->SetNClsTPC(60);// for now use only crossed pad rows
      TrackCutsTrigger->SetPID(AliPID::kProton, 0.75,4., false, 4.);  
      TrackCutsTrigger->SetDCAReCalculation(true);  //will I have this in Run3 for trigger???
      TrackCutsTrigger->SetDCAVtxZ(0.8);
      TrackCutsTrigger->SetDCAVtxXY(0.4);

      AntiTrackCutsTrigger = new AliFemtoDreamTrackCuts();
      AntiTrackCutsTrigger->SetPlotDCADist(true);
      AntiTrackCutsTrigger->SetIsMonteCarlo(isMC);

      AntiTrackCutsTrigger->SetFilterBit(128); //will I have this in Run3 for trigger???
      AntiTrackCutsTrigger->SetCutCharge(-1);
      AntiTrackCutsTrigger->SetPtRange(0.3, 5.);
      AntiTrackCutsTrigger->SetEtaRange(-0.9, 0.9);
      AntiTrackCutsTrigger->SetCutTPCCrossedRows(true, 60., 0.);
      //AntiTrackCutsTrigger->SetNClsTPC(60); // for now use only crossed pad rows
      AntiTrackCutsTrigger->SetPID(AliPID::kProton, 0.75,4., false, 4.); 
      AntiTrackCutsTrigger->SetDCAReCalculation(true);  //will I have this in Run3 for trigger???
      AntiTrackCutsTrigger->SetDCAVtxZ(0.8);
      AntiTrackCutsTrigger->SetDCAVtxXY(0.4);

      v0CutsTrigger = new AliFemtoDreamv0Cuts();
      v0CutsTrigger->SetIsMonteCarlo(isMC);
      v0CutsTrigger->SetPlotCPADist(true);

      Posv0DaugTrigger = new AliFemtoDreamTrackCuts(); // proton
      Posv0DaugTrigger->SetIsMonteCarlo(isMC);
      Posv0DaugTrigger->SetFillQALater(true);
      Posv0DaugTrigger->SetCutCharge(1);

      Negv0DaugTrigger = new AliFemtoDreamTrackCuts(); //pion
      Negv0DaugTrigger->SetIsMonteCarlo(isMC);
      Negv0DaugTrigger->SetFillQALater(true);
      Negv0DaugTrigger->SetCutCharge(-1);

      v0CutsTrigger->SetCutCharge(0);
      v0CutsTrigger->SetPtRange(0., 999.);
      v0CutsTrigger->SetCutCPA(0.9);
      v0CutsTrigger->SetCutDCADaugToPrimVtx(0.02);
      v0CutsTrigger->SetCutDCADaugTov0Vtx(3.);
      v0CutsTrigger->SetKaonRejection(0.49, 0.505);
      v0CutsTrigger->SetCutMaxDecayVtx(120);
      v0CutsTrigger->SetCutTransverseRadius(0.1, 120);
      v0CutsTrigger->SetCutInvMass(0.006);
      v0CutsTrigger->SetAxisInvMassPlots(400, 1.0, 1.2);

      Posv0DaugTrigger->SetCutTPCCrossedRows(true, 60, 0.);
      Negv0DaugTrigger->SetCutTPCCrossedRows(true, 60, 0.);
      Posv0DaugTrigger->SetEtaRange(-0.9, 0.9);
      Negv0DaugTrigger->SetEtaRange(-0.9, 0.9);
      Posv0DaugTrigger->SetPID(AliPID::kProton, 999., 7.);
      Negv0DaugTrigger->SetPID(AliPID::kPion, 999., 7.);
      Negv0DaugTrigger->SetDCAReCalculation(true);
      Posv0DaugTrigger->SetDCAReCalculation(true);

      v0CutsTrigger->SetPosDaugterTrackCuts(Posv0DaugTrigger);
      v0CutsTrigger->SetNegDaugterTrackCuts(Negv0DaugTrigger);
      v0CutsTrigger->SetPDGCodePosDaug(2212);  //Proton
      v0CutsTrigger->SetPDGCodeNegDaug(211);  //Pion
      v0CutsTrigger->SetPDGCodev0(3122);  //Lambda

      Antiv0CutsTrigger = new AliFemtoDreamv0Cuts();
      Antiv0CutsTrigger->SetIsMonteCarlo(isMC);
      Antiv0CutsTrigger->SetPlotCPADist(true);
      Antiv0CutsTrigger->SetPlotContrib(false);

      PosAntiv0DaugTrigger = new AliFemtoDreamTrackCuts();
      PosAntiv0DaugTrigger->SetIsMonteCarlo(isMC);
      PosAntiv0DaugTrigger->SetFillQALater(true);
      PosAntiv0DaugTrigger->SetCutCharge(1);

      NegAntiv0DaugTrigger =new AliFemtoDreamTrackCuts();
      NegAntiv0DaugTrigger->SetIsMonteCarlo(isMC);
      NegAntiv0DaugTrigger->SetFillQALater(true);
      NegAntiv0DaugTrigger->SetCutCharge(-1);

      Antiv0CutsTrigger->SetCutCharge(0);
      Antiv0CutsTrigger->SetPtRange(0., 999.);
      Antiv0CutsTrigger->SetCutCPA(0.9);
      Antiv0CutsTrigger->SetCutDCADaugToPrimVtx(0.02);
      Antiv0CutsTrigger->SetCutDCADaugTov0Vtx(3.);
      Antiv0CutsTrigger->SetKaonRejection(0.49, 0.505);
      Antiv0CutsTrigger->SetCutMaxDecayVtx(120);
      Antiv0CutsTrigger->SetCutTransverseRadius(0.1, 120);
      Antiv0CutsTrigger->SetCutInvMass(0.006);
      Antiv0CutsTrigger->SetAxisInvMassPlots(400, 1.0, 1.2);

      PosAntiv0DaugTrigger->SetCutTPCCrossedRows(true, 60, 0.);
      NegAntiv0DaugTrigger->SetCutTPCCrossedRows(true, 60, 0.);
      PosAntiv0DaugTrigger->SetEtaRange(-0.9, 0.9);
      NegAntiv0DaugTrigger->SetEtaRange(-0.9, 0.9);
      PosAntiv0DaugTrigger->SetPID(AliPID::kPion, 999., 7.);
      NegAntiv0DaugTrigger->SetPID(AliPID::kProton, 999., 7.);
      PosAntiv0DaugTrigger->SetDCAReCalculation(true);
      NegAntiv0DaugTrigger->SetDCAReCalculation(true);

      Antiv0CutsTrigger->SetPosDaugterTrackCuts(PosAntiv0DaugTrigger);
      Antiv0CutsTrigger->SetNegDaugterTrackCuts(NegAntiv0DaugTrigger);
      Antiv0CutsTrigger->SetPDGCodePosDaug(211);  //Pion
      Antiv0CutsTrigger->SetPDGCodeNegDaug(2212);  //Proton
      Antiv0CutsTrigger->SetPDGCodev0(-3122);  //Lambda
      if (!fullBlastQA) {
        evtCutsTrigger->SetMinimalBooking(true);
        TrackCutsTrigger->SetMinimalBooking(true);
        AntiTrackCutsTrigger->SetMinimalBooking(true);
        v0CutsTrigger->SetMinimalBooking(true);
        Antiv0CutsTrigger->SetMinimalBooking(true);
      }
    }

    if(suffixTrigger=="6"){
      Q3Limit = 1.5;
      evtCutsTrigger = new AliFemtoDreamEventCuts();
      evtCutsTrigger->UseDontWorryEvtCuts(false);
      evtCutsTrigger->SetZVtxPosition(-12., 12.);
      //evtCutsTrigger->CleanUpMult(false, false, false, true); will I have this in Run3 for trigger???

      TrackCutsTrigger =  new AliFemtoDreamTrackCuts();
      TrackCutsTrigger->SetPlotDCADist(true);
      TrackCutsTrigger->SetIsMonteCarlo(isMC);
      
      TrackCutsTrigger->SetFilterBit(128); //will I have this in Run3 for trigger???
      TrackCutsTrigger->SetCutCharge(1);
      TrackCutsTrigger->SetPtRange(0.3, 5.);
      TrackCutsTrigger->SetEtaRange(-0.9, 0.9);
      TrackCutsTrigger->SetCutTPCCrossedRows(true, 60., 0.);
      //TrackCutsTrigger->SetNClsTPC(60);// for now use only crossed pad rows
      TrackCutsTrigger->SetPID(AliPID::kProton, 0.75,4., false, 4.);  
      TrackCutsTrigger->SetDCAReCalculation(true);  //will I have this in Run3 for trigger???
      TrackCutsTrigger->SetDCAVtxZ(0.8);
      TrackCutsTrigger->SetDCAVtxXY(0.4);

      AntiTrackCutsTrigger = new AliFemtoDreamTrackCuts();
      AntiTrackCutsTrigger->SetPlotDCADist(true);
      AntiTrackCutsTrigger->SetIsMonteCarlo(isMC);

      AntiTrackCutsTrigger->SetFilterBit(128); //will I have this in Run3 for trigger???
      AntiTrackCutsTrigger->SetCutCharge(-1);
      AntiTrackCutsTrigger->SetPtRange(0.3, 5.);
      AntiTrackCutsTrigger->SetEtaRange(-0.9, 0.9);
      AntiTrackCutsTrigger->SetCutTPCCrossedRows(true, 60., 0.);
      //AntiTrackCutsTrigger->SetNClsTPC(60); // for now use only crossed pad rows
      AntiTrackCutsTrigger->SetPID(AliPID::kProton, 0.75,4., false, 4.); 
      AntiTrackCutsTrigger->SetDCAReCalculation(true);  //will I have this in Run3 for trigger???
      AntiTrackCutsTrigger->SetDCAVtxZ(0.8);
      AntiTrackCutsTrigger->SetDCAVtxXY(0.4);

      v0CutsTrigger = new AliFemtoDreamv0Cuts();
      v0CutsTrigger->SetIsMonteCarlo(isMC);
      v0CutsTrigger->SetPlotCPADist(true);

      Posv0DaugTrigger = new AliFemtoDreamTrackCuts(); // proton
      Posv0DaugTrigger->SetIsMonteCarlo(isMC);
      Posv0DaugTrigger->SetFillQALater(true);
      Posv0DaugTrigger->SetCutCharge(1);
      Posv0DaugTrigger->SetCheckPileUp(true);

      Negv0DaugTrigger = new AliFemtoDreamTrackCuts(); //pion
      Negv0DaugTrigger->SetIsMonteCarlo(isMC);
      Negv0DaugTrigger->SetFillQALater(true);
      Negv0DaugTrigger->SetCutCharge(-1);
      Negv0DaugTrigger->SetCheckPileUp(true);

      v0CutsTrigger->SetCutCharge(0);
      v0CutsTrigger->SetPtRange(0., 999.);
      v0CutsTrigger->SetCutCPA(0.9);
      v0CutsTrigger->SetCutDCADaugToPrimVtx(0.02);
      v0CutsTrigger->SetCutDCADaugTov0Vtx(3.);
      v0CutsTrigger->SetKaonRejection(0.49, 0.505);
      v0CutsTrigger->SetCutMaxDecayVtx(120);
      v0CutsTrigger->SetCutTransverseRadius(0.1, 120);
      v0CutsTrigger->SetCutInvMass(0.006);
      v0CutsTrigger->SetAxisInvMassPlots(400, 1.0, 1.2);

      Posv0DaugTrigger->SetCutTPCCrossedRows(true, 60, 0.);
      Negv0DaugTrigger->SetCutTPCCrossedRows(true, 60, 0.);
      Posv0DaugTrigger->SetEtaRange(-0.9, 0.9);
      Negv0DaugTrigger->SetEtaRange(-0.9, 0.9);
      Posv0DaugTrigger->SetPID(AliPID::kProton, 999., 7.);
      Negv0DaugTrigger->SetPID(AliPID::kPion, 999., 7.);
      Negv0DaugTrigger->SetDCAReCalculation(true);
      Posv0DaugTrigger->SetDCAReCalculation(true);

      v0CutsTrigger->SetPosDaugterTrackCuts(Posv0DaugTrigger);
      v0CutsTrigger->SetNegDaugterTrackCuts(Negv0DaugTrigger);
      v0CutsTrigger->SetPDGCodePosDaug(2212);  //Proton
      v0CutsTrigger->SetPDGCodeNegDaug(211);  //Pion
      v0CutsTrigger->SetPDGCodev0(3122);  //Lambda

      Antiv0CutsTrigger = new AliFemtoDreamv0Cuts();
      Antiv0CutsTrigger->SetIsMonteCarlo(isMC);
      Antiv0CutsTrigger->SetPlotCPADist(true);
      Antiv0CutsTrigger->SetPlotContrib(false);

      PosAntiv0DaugTrigger = new AliFemtoDreamTrackCuts();
      PosAntiv0DaugTrigger->SetIsMonteCarlo(isMC);
      PosAntiv0DaugTrigger->SetFillQALater(true);
      PosAntiv0DaugTrigger->SetCutCharge(1);
      PosAntiv0DaugTrigger->SetCheckPileUp(true);

      NegAntiv0DaugTrigger =new AliFemtoDreamTrackCuts();
      NegAntiv0DaugTrigger->SetIsMonteCarlo(isMC);
      NegAntiv0DaugTrigger->SetFillQALater(true);
      NegAntiv0DaugTrigger->SetCutCharge(-1);
      NegAntiv0DaugTrigger->SetCheckPileUp(true);

      Antiv0CutsTrigger->SetCutCharge(0);
      Antiv0CutsTrigger->SetPtRange(0., 999.);
      Antiv0CutsTrigger->SetCutCPA(0.9);
      Antiv0CutsTrigger->SetCutDCADaugToPrimVtx(0.02);
      Antiv0CutsTrigger->SetCutDCADaugTov0Vtx(3.);
      Antiv0CutsTrigger->SetKaonRejection(0.49, 0.505);
      Antiv0CutsTrigger->SetCutMaxDecayVtx(120);
      Antiv0CutsTrigger->SetCutTransverseRadius(0.1, 120);
      Antiv0CutsTrigger->SetCutInvMass(0.006);
      Antiv0CutsTrigger->SetAxisInvMassPlots(400, 1.0, 1.2);

      PosAntiv0DaugTrigger->SetCutTPCCrossedRows(true, 60, 0.);
      NegAntiv0DaugTrigger->SetCutTPCCrossedRows(true, 60, 0.);
      PosAntiv0DaugTrigger->SetEtaRange(-0.9, 0.9);
      NegAntiv0DaugTrigger->SetEtaRange(-0.9, 0.9);
      PosAntiv0DaugTrigger->SetPID(AliPID::kPion, 999., 7.);
      NegAntiv0DaugTrigger->SetPID(AliPID::kProton, 999., 7.);
      PosAntiv0DaugTrigger->SetDCAReCalculation(true);
      NegAntiv0DaugTrigger->SetDCAReCalculation(true);

      Antiv0CutsTrigger->SetPosDaugterTrackCuts(PosAntiv0DaugTrigger);
      Antiv0CutsTrigger->SetNegDaugterTrackCuts(NegAntiv0DaugTrigger);
      Antiv0CutsTrigger->SetPDGCodePosDaug(211);  //Pion
      Antiv0CutsTrigger->SetPDGCodeNegDaug(2212);  //Proton
      Antiv0CutsTrigger->SetPDGCodev0(-3122);  //Lambda
      if (!fullBlastQA) {
        evtCutsTrigger->SetMinimalBooking(true);
        TrackCutsTrigger->SetMinimalBooking(true);
        AntiTrackCutsTrigger->SetMinimalBooking(true);
        v0CutsTrigger->SetMinimalBooking(true);
        Antiv0CutsTrigger->SetMinimalBooking(true);
      }
    }

    if(suffixTrigger=="7"){
      Q3Limit = 1.5;
      evtCutsTrigger = new AliFemtoDreamEventCuts();
      evtCutsTrigger->UseDontWorryEvtCuts(false);
      evtCutsTrigger->SetZVtxPosition(-12., 12.);
      //evtCutsTrigger->CleanUpMult(false, false, false, true); will I have this in Run3 for trigger???

      TrackCutsTrigger =  new AliFemtoDreamTrackCuts();
      TrackCutsTrigger->SetPlotDCADist(true);
      TrackCutsTrigger->SetIsMonteCarlo(isMC);
      
      TrackCutsTrigger->SetFilterBit(128); //will I have this in Run3 for trigger???
      TrackCutsTrigger->SetCutCharge(1);
      TrackCutsTrigger->SetPtRange(0.3, 5.);
      TrackCutsTrigger->SetEtaRange(-0.9, 0.9);
      TrackCutsTrigger->SetCutTPCCrossedRows(true, 60., 0.);
      //TrackCutsTrigger->SetNClsTPC(60);// for now use only crossed pad rows
      TrackCutsTrigger->SetPID(AliPID::kProton, 0.75,4., false, 4.);  
      TrackCutsTrigger->SetDCAReCalculation(true);  //will I have this in Run3 for trigger???
      TrackCutsTrigger->SetDCAVtxZ(0.4);
      TrackCutsTrigger->SetDCAVtxXY(0.2);

      AntiTrackCutsTrigger = new AliFemtoDreamTrackCuts();
      AntiTrackCutsTrigger->SetPlotDCADist(true);
      AntiTrackCutsTrigger->SetIsMonteCarlo(isMC);

      AntiTrackCutsTrigger->SetFilterBit(128); //will I have this in Run3 for trigger???
      AntiTrackCutsTrigger->SetCutCharge(-1);
      AntiTrackCutsTrigger->SetPtRange(0.3, 5.);
      AntiTrackCutsTrigger->SetEtaRange(-0.9, 0.9);
      AntiTrackCutsTrigger->SetCutTPCCrossedRows(true, 60., 0.);
      //AntiTrackCutsTrigger->SetNClsTPC(60); // for now use only crossed pad rows
      AntiTrackCutsTrigger->SetPID(AliPID::kProton, 0.75,4., false, 4.); 
      AntiTrackCutsTrigger->SetDCAReCalculation(true);  //will I have this in Run3 for trigger???
      AntiTrackCutsTrigger->SetDCAVtxZ(0.4);
      AntiTrackCutsTrigger->SetDCAVtxXY(0.2);

      v0CutsTrigger = new AliFemtoDreamv0Cuts();
      v0CutsTrigger->SetIsMonteCarlo(isMC);
      v0CutsTrigger->SetPlotCPADist(true);

      Posv0DaugTrigger = new AliFemtoDreamTrackCuts(); // proton
      Posv0DaugTrigger->SetIsMonteCarlo(isMC);
      Posv0DaugTrigger->SetFillQALater(true);
      Posv0DaugTrigger->SetCutCharge(1);
      Posv0DaugTrigger->SetCheckPileUp(true);

      Negv0DaugTrigger = new AliFemtoDreamTrackCuts(); //pion
      Negv0DaugTrigger->SetIsMonteCarlo(isMC);
      Negv0DaugTrigger->SetFillQALater(true);
      Negv0DaugTrigger->SetCutCharge(-1);
      Negv0DaugTrigger->SetCheckPileUp(true);

      v0CutsTrigger->SetCutCharge(0);
      v0CutsTrigger->SetPtRange(0., 999.);
      v0CutsTrigger->SetCutCPA(0.95);
      v0CutsTrigger->SetCutDCADaugToPrimVtx(0.02);
      v0CutsTrigger->SetCutDCADaugTov0Vtx(2.5);
      v0CutsTrigger->SetKaonRejection(0.49, 0.505);
      v0CutsTrigger->SetCutMaxDecayVtx(120);
      v0CutsTrigger->SetCutTransverseRadius(0.1, 120);
      v0CutsTrigger->SetCutInvMass(0.006);
      v0CutsTrigger->SetAxisInvMassPlots(400, 1.0, 1.2);

      Posv0DaugTrigger->SetCutTPCCrossedRows(true, 60, 0.);
      Negv0DaugTrigger->SetCutTPCCrossedRows(true, 60, 0.);
      Posv0DaugTrigger->SetEtaRange(-0.9, 0.9);
      Negv0DaugTrigger->SetEtaRange(-0.9, 0.9);
      Posv0DaugTrigger->SetPID(AliPID::kProton, 999., 7.);
      Negv0DaugTrigger->SetPID(AliPID::kPion, 999., 7.);
      Negv0DaugTrigger->SetDCAReCalculation(true);
      Posv0DaugTrigger->SetDCAReCalculation(true);

      v0CutsTrigger->SetPosDaugterTrackCuts(Posv0DaugTrigger);
      v0CutsTrigger->SetNegDaugterTrackCuts(Negv0DaugTrigger);
      v0CutsTrigger->SetPDGCodePosDaug(2212);  //Proton
      v0CutsTrigger->SetPDGCodeNegDaug(211);  //Pion
      v0CutsTrigger->SetPDGCodev0(3122);  //Lambda

      Antiv0CutsTrigger = new AliFemtoDreamv0Cuts();
      Antiv0CutsTrigger->SetIsMonteCarlo(isMC);
      Antiv0CutsTrigger->SetPlotCPADist(true);
      Antiv0CutsTrigger->SetPlotContrib(false);

      PosAntiv0DaugTrigger = new AliFemtoDreamTrackCuts();
      PosAntiv0DaugTrigger->SetIsMonteCarlo(isMC);
      PosAntiv0DaugTrigger->SetFillQALater(true);
      PosAntiv0DaugTrigger->SetCutCharge(1);
      PosAntiv0DaugTrigger->SetCheckPileUp(true);

      NegAntiv0DaugTrigger =new AliFemtoDreamTrackCuts();
      NegAntiv0DaugTrigger->SetIsMonteCarlo(isMC);
      NegAntiv0DaugTrigger->SetFillQALater(true);
      NegAntiv0DaugTrigger->SetCutCharge(-1);
      NegAntiv0DaugTrigger->SetCheckPileUp(true);

      Antiv0CutsTrigger->SetCutCharge(0);
      Antiv0CutsTrigger->SetPtRange(0., 999.);
      Antiv0CutsTrigger->SetCutCPA(0.95);
      Antiv0CutsTrigger->SetCutDCADaugToPrimVtx(0.02);
      Antiv0CutsTrigger->SetCutDCADaugTov0Vtx(2.5);
      Antiv0CutsTrigger->SetKaonRejection(0.49, 0.505);
      Antiv0CutsTrigger->SetCutMaxDecayVtx(120);
      Antiv0CutsTrigger->SetCutTransverseRadius(0.1, 120);
      Antiv0CutsTrigger->SetCutInvMass(0.006);
      Antiv0CutsTrigger->SetAxisInvMassPlots(400, 1.0, 1.2);

      PosAntiv0DaugTrigger->SetCutTPCCrossedRows(true, 60, 0.);
      NegAntiv0DaugTrigger->SetCutTPCCrossedRows(true, 60, 0.);
      PosAntiv0DaugTrigger->SetEtaRange(-0.9, 0.9);
      NegAntiv0DaugTrigger->SetEtaRange(-0.9, 0.9);
      PosAntiv0DaugTrigger->SetPID(AliPID::kPion, 999., 7.);
      NegAntiv0DaugTrigger->SetPID(AliPID::kProton, 999., 7.);
      PosAntiv0DaugTrigger->SetDCAReCalculation(true);
      NegAntiv0DaugTrigger->SetDCAReCalculation(true);

      Antiv0CutsTrigger->SetPosDaugterTrackCuts(PosAntiv0DaugTrigger);
      Antiv0CutsTrigger->SetNegDaugterTrackCuts(NegAntiv0DaugTrigger);
      Antiv0CutsTrigger->SetPDGCodePosDaug(211);  //Pion
      Antiv0CutsTrigger->SetPDGCodeNegDaug(2212);  //Proton
      Antiv0CutsTrigger->SetPDGCodev0(-3122);  //Lambda

      if (!fullBlastQA) {
        evtCutsTrigger->SetMinimalBooking(true);
        TrackCutsTrigger->SetMinimalBooking(true);
        AntiTrackCutsTrigger->SetMinimalBooking(true);
        v0CutsTrigger->SetMinimalBooking(true);
        Antiv0CutsTrigger->SetMinimalBooking(true);
      }
    }

    if(suffixTrigger=="8"){
      Q3Limit = 1.;
      TriggerOnSample = true;
      Q3LimitSample = 1.5;
            evtCutsTrigger = new AliFemtoDreamEventCuts();
      evtCutsTrigger->UseDontWorryEvtCuts(false);
      evtCutsTrigger->SetZVtxPosition(-12., 12.);
      //evtCutsTrigger->CleanUpMult(false, false, false, true); will I have this in Run3 for trigger???

      TrackCutsTrigger =  new AliFemtoDreamTrackCuts();
      TrackCutsTrigger->SetPlotDCADist(true);
      TrackCutsTrigger->SetIsMonteCarlo(isMC);
      
      TrackCutsTrigger->SetFilterBit(128); //will I have this in Run3 for trigger???
      TrackCutsTrigger->SetCutCharge(1);
      TrackCutsTrigger->SetPtRange(0.3, 5.);
      TrackCutsTrigger->SetEtaRange(-0.9, 0.9);
      TrackCutsTrigger->SetCutTPCCrossedRows(true, 60., 0.);
      //TrackCutsTrigger->SetNClsTPC(60);// for now use only crossed pad rows
      TrackCutsTrigger->SetPID(AliPID::kProton, 0.75,4., false, 4.);  
      TrackCutsTrigger->SetDCAReCalculation(true);  //will I have this in Run3 for trigger???
      TrackCutsTrigger->SetDCAVtxZ(0.4);
      TrackCutsTrigger->SetDCAVtxXY(0.2);

      AntiTrackCutsTrigger = new AliFemtoDreamTrackCuts();
      AntiTrackCutsTrigger->SetPlotDCADist(true);
      AntiTrackCutsTrigger->SetIsMonteCarlo(isMC);

      AntiTrackCutsTrigger->SetFilterBit(128); //will I have this in Run3 for trigger???
      AntiTrackCutsTrigger->SetCutCharge(-1);
      AntiTrackCutsTrigger->SetPtRange(0.3, 5.);
      AntiTrackCutsTrigger->SetEtaRange(-0.9, 0.9);
      AntiTrackCutsTrigger->SetCutTPCCrossedRows(true, 60., 0.);
      //AntiTrackCutsTrigger->SetNClsTPC(60); // for now use only crossed pad rows
      AntiTrackCutsTrigger->SetPID(AliPID::kProton, 0.75,4., false, 4.); 
      AntiTrackCutsTrigger->SetDCAReCalculation(true);  //will I have this in Run3 for trigger???
      AntiTrackCutsTrigger->SetDCAVtxZ(0.4);
      AntiTrackCutsTrigger->SetDCAVtxXY(0.2);

      v0CutsTrigger = new AliFemtoDreamv0Cuts();
      v0CutsTrigger->SetIsMonteCarlo(isMC);
      v0CutsTrigger->SetPlotCPADist(true);

      Posv0DaugTrigger = new AliFemtoDreamTrackCuts(); // proton
      Posv0DaugTrigger->SetIsMonteCarlo(isMC);
      Posv0DaugTrigger->SetFillQALater(true);
      Posv0DaugTrigger->SetCutCharge(1);
      Posv0DaugTrigger->SetCheckPileUp(true);

      Negv0DaugTrigger = new AliFemtoDreamTrackCuts(); //pion
      Negv0DaugTrigger->SetIsMonteCarlo(isMC);
      Negv0DaugTrigger->SetFillQALater(true);
      Negv0DaugTrigger->SetCutCharge(-1);
      Negv0DaugTrigger->SetCheckPileUp(true);

      v0CutsTrigger->SetCutCharge(0);
      v0CutsTrigger->SetPtRange(0., 999.);
      v0CutsTrigger->SetCutCPA(0.95);
      v0CutsTrigger->SetCutDCADaugToPrimVtx(0.02);
      v0CutsTrigger->SetCutDCADaugTov0Vtx(2.5);
      v0CutsTrigger->SetKaonRejection(0.49, 0.505);
      v0CutsTrigger->SetCutMaxDecayVtx(120);
      v0CutsTrigger->SetCutTransverseRadius(0.1, 120);
      v0CutsTrigger->SetCutInvMass(0.006);
      v0CutsTrigger->SetAxisInvMassPlots(400, 1.0, 1.2);

      Posv0DaugTrigger->SetCutTPCCrossedRows(true, 60, 0.);
      Negv0DaugTrigger->SetCutTPCCrossedRows(true, 60, 0.);
      Posv0DaugTrigger->SetEtaRange(-0.9, 0.9);
      Negv0DaugTrigger->SetEtaRange(-0.9, 0.9);
      Posv0DaugTrigger->SetPID(AliPID::kProton, 999., 7.);
      Negv0DaugTrigger->SetPID(AliPID::kPion, 999., 7.);
      Negv0DaugTrigger->SetDCAReCalculation(true);
      Posv0DaugTrigger->SetDCAReCalculation(true);

      v0CutsTrigger->SetPosDaugterTrackCuts(Posv0DaugTrigger);
      v0CutsTrigger->SetNegDaugterTrackCuts(Negv0DaugTrigger);
      v0CutsTrigger->SetPDGCodePosDaug(2212);  //Proton
      v0CutsTrigger->SetPDGCodeNegDaug(211);  //Pion
      v0CutsTrigger->SetPDGCodev0(3122);  //Lambda

      Antiv0CutsTrigger = new AliFemtoDreamv0Cuts();
      Antiv0CutsTrigger->SetIsMonteCarlo(isMC);
      Antiv0CutsTrigger->SetPlotCPADist(true);
      Antiv0CutsTrigger->SetPlotContrib(false);

      PosAntiv0DaugTrigger = new AliFemtoDreamTrackCuts();
      PosAntiv0DaugTrigger->SetIsMonteCarlo(isMC);
      PosAntiv0DaugTrigger->SetFillQALater(true);
      PosAntiv0DaugTrigger->SetCutCharge(1);
      PosAntiv0DaugTrigger->SetCheckPileUp(true);

      NegAntiv0DaugTrigger =new AliFemtoDreamTrackCuts();
      NegAntiv0DaugTrigger->SetIsMonteCarlo(isMC);
      NegAntiv0DaugTrigger->SetFillQALater(true);
      NegAntiv0DaugTrigger->SetCutCharge(-1);
      NegAntiv0DaugTrigger->SetCheckPileUp(true);

      Antiv0CutsTrigger->SetCutCharge(0);
      Antiv0CutsTrigger->SetPtRange(0., 999.);
      Antiv0CutsTrigger->SetCutCPA(0.95);
      Antiv0CutsTrigger->SetCutDCADaugToPrimVtx(0.02);
      Antiv0CutsTrigger->SetCutDCADaugTov0Vtx(2.5);
      Antiv0CutsTrigger->SetKaonRejection(0.49, 0.505);
      Antiv0CutsTrigger->SetCutMaxDecayVtx(120);
      Antiv0CutsTrigger->SetCutTransverseRadius(0.1, 120);
      Antiv0CutsTrigger->SetCutInvMass(0.006);
      Antiv0CutsTrigger->SetAxisInvMassPlots(400, 1.0, 1.2);

      PosAntiv0DaugTrigger->SetCutTPCCrossedRows(true, 60, 0.);
      NegAntiv0DaugTrigger->SetCutTPCCrossedRows(true, 60, 0.);
      PosAntiv0DaugTrigger->SetEtaRange(-0.9, 0.9);
      NegAntiv0DaugTrigger->SetEtaRange(-0.9, 0.9);
      PosAntiv0DaugTrigger->SetPID(AliPID::kPion, 999., 7.);
      NegAntiv0DaugTrigger->SetPID(AliPID::kProton, 999., 7.);
      PosAntiv0DaugTrigger->SetDCAReCalculation(true);
      NegAntiv0DaugTrigger->SetDCAReCalculation(true);

      Antiv0CutsTrigger->SetPosDaugterTrackCuts(PosAntiv0DaugTrigger);
      Antiv0CutsTrigger->SetNegDaugterTrackCuts(NegAntiv0DaugTrigger);
      Antiv0CutsTrigger->SetPDGCodePosDaug(211);  //Pion
      Antiv0CutsTrigger->SetPDGCodeNegDaug(2212);  //Proton
      Antiv0CutsTrigger->SetPDGCodev0(-3122);  //Lambda

      if (!fullBlastQA) {
        evtCutsTrigger->SetMinimalBooking(true);
        TrackCutsTrigger->SetMinimalBooking(true);
        AntiTrackCutsTrigger->SetMinimalBooking(true);
        v0CutsTrigger->SetMinimalBooking(true);
        Antiv0CutsTrigger->SetMinimalBooking(true);
      }
    }


     if(suffixTrigger=="9"){
      Q3Limit = 1.;
      TriggerOnSample = false;
      Q3LimitSample = 1.5;
            evtCutsTrigger = new AliFemtoDreamEventCuts();
      evtCutsTrigger->UseDontWorryEvtCuts(false);
      evtCutsTrigger->SetZVtxPosition(-12., 12.);
      //evtCutsTrigger->CleanUpMult(false, false, false, true); will I have this in Run3 for trigger???

      TrackCutsTrigger =  new AliFemtoDreamTrackCuts();
      TrackCutsTrigger->SetPlotDCADist(true);
      TrackCutsTrigger->SetIsMonteCarlo(isMC);
      
      TrackCutsTrigger->SetFilterBit(128); //will I have this in Run3 for trigger???
      TrackCutsTrigger->SetCutCharge(1);
      TrackCutsTrigger->SetPtRange(0.3, 5.);
      TrackCutsTrigger->SetEtaRange(-0.9, 0.9);
      TrackCutsTrigger->SetCutTPCCrossedRows(true, 60., 0.);
      //TrackCutsTrigger->SetNClsTPC(60);// for now use only crossed pad rows
      TrackCutsTrigger->SetPID(AliPID::kProton, 0.75,4., false, 4.);  
      TrackCutsTrigger->SetDCAReCalculation(true);  //will I have this in Run3 for trigger???
      TrackCutsTrigger->SetDCAVtxZ(0.4);
      TrackCutsTrigger->SetDCAVtxXY(0.2);

      AntiTrackCutsTrigger = new AliFemtoDreamTrackCuts();
      AntiTrackCutsTrigger->SetPlotDCADist(true);
      AntiTrackCutsTrigger->SetIsMonteCarlo(isMC);

      AntiTrackCutsTrigger->SetFilterBit(128); //will I have this in Run3 for trigger???
      AntiTrackCutsTrigger->SetCutCharge(-1);
      AntiTrackCutsTrigger->SetPtRange(0.3, 5.);
      AntiTrackCutsTrigger->SetEtaRange(-0.9, 0.9);
      AntiTrackCutsTrigger->SetCutTPCCrossedRows(true, 60., 0.);
      //AntiTrackCutsTrigger->SetNClsTPC(60); // for now use only crossed pad rows
      AntiTrackCutsTrigger->SetPID(AliPID::kProton, 0.75,4., false, 4.); 
      AntiTrackCutsTrigger->SetDCAReCalculation(true);  //will I have this in Run3 for trigger???
      AntiTrackCutsTrigger->SetDCAVtxZ(0.4);
      AntiTrackCutsTrigger->SetDCAVtxXY(0.2);

      v0CutsTrigger = new AliFemtoDreamv0Cuts();
      v0CutsTrigger->SetIsMonteCarlo(isMC);
      v0CutsTrigger->SetPlotCPADist(true);

      Posv0DaugTrigger = new AliFemtoDreamTrackCuts(); // proton
      Posv0DaugTrigger->SetIsMonteCarlo(isMC);
      Posv0DaugTrigger->SetFillQALater(true);
      Posv0DaugTrigger->SetCutCharge(1);
      Posv0DaugTrigger->SetCheckPileUp(true);

      Negv0DaugTrigger = new AliFemtoDreamTrackCuts(); //pion
      Negv0DaugTrigger->SetIsMonteCarlo(isMC);
      Negv0DaugTrigger->SetFillQALater(true);
      Negv0DaugTrigger->SetCutCharge(-1);
      Negv0DaugTrigger->SetCheckPileUp(true);

      v0CutsTrigger->SetCutCharge(0);
      v0CutsTrigger->SetPtRange(0., 999.);
      v0CutsTrigger->SetCutCPA(0.95);
      v0CutsTrigger->SetCutDCADaugToPrimVtx(0.02);
      v0CutsTrigger->SetCutDCADaugTov0Vtx(2.5);
      v0CutsTrigger->SetKaonRejection(0.49, 0.505);
      v0CutsTrigger->SetCutMaxDecayVtx(120);
      v0CutsTrigger->SetCutTransverseRadius(0.1, 120);
      v0CutsTrigger->SetCutInvMass(0.006);
      v0CutsTrigger->SetAxisInvMassPlots(400, 1.0, 1.2);

      Posv0DaugTrigger->SetCutTPCCrossedRows(true, 60, 0.);
      Negv0DaugTrigger->SetCutTPCCrossedRows(true, 60, 0.);
      Posv0DaugTrigger->SetEtaRange(-0.9, 0.9);
      Negv0DaugTrigger->SetEtaRange(-0.9, 0.9);
      Posv0DaugTrigger->SetPID(AliPID::kProton, 999., 7.);
      Negv0DaugTrigger->SetPID(AliPID::kPion, 999., 7.);
      Negv0DaugTrigger->SetDCAReCalculation(true);
      Posv0DaugTrigger->SetDCAReCalculation(true);

      v0CutsTrigger->SetPosDaugterTrackCuts(Posv0DaugTrigger);
      v0CutsTrigger->SetNegDaugterTrackCuts(Negv0DaugTrigger);
      v0CutsTrigger->SetPDGCodePosDaug(2212);  //Proton
      v0CutsTrigger->SetPDGCodeNegDaug(211);  //Pion
      v0CutsTrigger->SetPDGCodev0(3122);  //Lambda

      Antiv0CutsTrigger = new AliFemtoDreamv0Cuts();
      Antiv0CutsTrigger->SetIsMonteCarlo(isMC);
      Antiv0CutsTrigger->SetPlotCPADist(true);
      Antiv0CutsTrigger->SetPlotContrib(false);

      PosAntiv0DaugTrigger = new AliFemtoDreamTrackCuts();
      PosAntiv0DaugTrigger->SetIsMonteCarlo(isMC);
      PosAntiv0DaugTrigger->SetFillQALater(true);
      PosAntiv0DaugTrigger->SetCutCharge(1);
      PosAntiv0DaugTrigger->SetCheckPileUp(true);

      NegAntiv0DaugTrigger =new AliFemtoDreamTrackCuts();
      NegAntiv0DaugTrigger->SetIsMonteCarlo(isMC);
      NegAntiv0DaugTrigger->SetFillQALater(true);
      NegAntiv0DaugTrigger->SetCutCharge(-1);
      NegAntiv0DaugTrigger->SetCheckPileUp(true);

      Antiv0CutsTrigger->SetCutCharge(0);
      Antiv0CutsTrigger->SetPtRange(0., 999.);
      Antiv0CutsTrigger->SetCutCPA(0.95);
      Antiv0CutsTrigger->SetCutDCADaugToPrimVtx(0.02);
      Antiv0CutsTrigger->SetCutDCADaugTov0Vtx(2.5);
      Antiv0CutsTrigger->SetKaonRejection(0.49, 0.505);
      Antiv0CutsTrigger->SetCutMaxDecayVtx(120);
      Antiv0CutsTrigger->SetCutTransverseRadius(0.1, 120);
      Antiv0CutsTrigger->SetCutInvMass(0.006);
      Antiv0CutsTrigger->SetAxisInvMassPlots(400, 1.0, 1.2);

      PosAntiv0DaugTrigger->SetCutTPCCrossedRows(true, 60, 0.);
      NegAntiv0DaugTrigger->SetCutTPCCrossedRows(true, 60, 0.);
      PosAntiv0DaugTrigger->SetEtaRange(-0.9, 0.9);
      NegAntiv0DaugTrigger->SetEtaRange(-0.9, 0.9);
      PosAntiv0DaugTrigger->SetPID(AliPID::kPion, 999., 7.);
      NegAntiv0DaugTrigger->SetPID(AliPID::kProton, 999., 7.);
      PosAntiv0DaugTrigger->SetDCAReCalculation(true);
      NegAntiv0DaugTrigger->SetDCAReCalculation(true);

      Antiv0CutsTrigger->SetPosDaugterTrackCuts(PosAntiv0DaugTrigger);
      Antiv0CutsTrigger->SetNegDaugterTrackCuts(NegAntiv0DaugTrigger);
      Antiv0CutsTrigger->SetPDGCodePosDaug(211);  //Pion
      Antiv0CutsTrigger->SetPDGCodeNegDaug(2212);  //Proton
      Antiv0CutsTrigger->SetPDGCodev0(-3122);  //Lambda

      if (!fullBlastQA) {
        evtCutsTrigger->SetMinimalBooking(true);
        TrackCutsTrigger->SetMinimalBooking(true);
        AntiTrackCutsTrigger->SetMinimalBooking(true);
        v0CutsTrigger->SetMinimalBooking(true);
        Antiv0CutsTrigger->SetMinimalBooking(true);
      }
    }

    if(suffixTrigger=="10"){
      Q3Limit = 0.;
      TriggerOnSample = true;
      Q3LimitSample = 1.;
            evtCutsTrigger = new AliFemtoDreamEventCuts();
      evtCutsTrigger->UseDontWorryEvtCuts(false);
      evtCutsTrigger->SetZVtxPosition(-12., 12.);
      //evtCutsTrigger->CleanUpMult(false, false, false, true); will I have this in Run3 for trigger???

      TrackCutsTrigger =  new AliFemtoDreamTrackCuts();
      TrackCutsTrigger->SetPlotDCADist(true);
      TrackCutsTrigger->SetIsMonteCarlo(isMC);
      
      TrackCutsTrigger->SetFilterBit(128); //will I have this in Run3 for trigger???
      TrackCutsTrigger->SetCutCharge(1);
      TrackCutsTrigger->SetPtRange(0.3, 5.);
      TrackCutsTrigger->SetEtaRange(-0.9, 0.9);
      TrackCutsTrigger->SetCutTPCCrossedRows(true, 60., 0.);
      //TrackCutsTrigger->SetNClsTPC(60);// for now use only crossed pad rows
      TrackCutsTrigger->SetPID(AliPID::kProton, 0.75,4., false, 4.);  
      TrackCutsTrigger->SetDCAReCalculation(true);  //will I have this in Run3 for trigger???
      TrackCutsTrigger->SetDCAVtxZ(0.4);
      TrackCutsTrigger->SetDCAVtxXY(0.2);

      AntiTrackCutsTrigger = new AliFemtoDreamTrackCuts();
      AntiTrackCutsTrigger->SetPlotDCADist(true);
      AntiTrackCutsTrigger->SetIsMonteCarlo(isMC);

      AntiTrackCutsTrigger->SetFilterBit(128); //will I have this in Run3 for trigger???
      AntiTrackCutsTrigger->SetCutCharge(-1);
      AntiTrackCutsTrigger->SetPtRange(0.3, 5.);
      AntiTrackCutsTrigger->SetEtaRange(-0.9, 0.9);
      AntiTrackCutsTrigger->SetCutTPCCrossedRows(true, 60., 0.);
      //AntiTrackCutsTrigger->SetNClsTPC(60); // for now use only crossed pad rows
      AntiTrackCutsTrigger->SetPID(AliPID::kProton, 0.75,4., false, 4.); 
      AntiTrackCutsTrigger->SetDCAReCalculation(true);  //will I have this in Run3 for trigger???
      AntiTrackCutsTrigger->SetDCAVtxZ(0.4);
      AntiTrackCutsTrigger->SetDCAVtxXY(0.2);

      v0CutsTrigger = new AliFemtoDreamv0Cuts();
      v0CutsTrigger->SetIsMonteCarlo(isMC);
      v0CutsTrigger->SetPlotCPADist(true);

      Posv0DaugTrigger = new AliFemtoDreamTrackCuts(); // proton
      Posv0DaugTrigger->SetIsMonteCarlo(isMC);
      Posv0DaugTrigger->SetFillQALater(true);
      Posv0DaugTrigger->SetCutCharge(1);
      Posv0DaugTrigger->SetCheckPileUp(true);

      Negv0DaugTrigger = new AliFemtoDreamTrackCuts(); //pion
      Negv0DaugTrigger->SetIsMonteCarlo(isMC);
      Negv0DaugTrigger->SetFillQALater(true);
      Negv0DaugTrigger->SetCutCharge(-1);
      Negv0DaugTrigger->SetCheckPileUp(true);

      v0CutsTrigger->SetCutCharge(0);
      v0CutsTrigger->SetPtRange(0., 999.);
      v0CutsTrigger->SetCutCPA(0.95);
      v0CutsTrigger->SetCutDCADaugToPrimVtx(0.02);
      v0CutsTrigger->SetCutDCADaugTov0Vtx(2.5);
      v0CutsTrigger->SetKaonRejection(0.49, 0.505);
      v0CutsTrigger->SetCutMaxDecayVtx(120);
      v0CutsTrigger->SetCutTransverseRadius(0.1, 120);
      v0CutsTrigger->SetCutInvMass(0.006);
      v0CutsTrigger->SetAxisInvMassPlots(400, 1.0, 1.2);

      Posv0DaugTrigger->SetCutTPCCrossedRows(true, 60, 0.);
      Negv0DaugTrigger->SetCutTPCCrossedRows(true, 60, 0.);
      Posv0DaugTrigger->SetEtaRange(-0.9, 0.9);
      Negv0DaugTrigger->SetEtaRange(-0.9, 0.9);
      Posv0DaugTrigger->SetPID(AliPID::kProton, 999., 7.);
      Negv0DaugTrigger->SetPID(AliPID::kPion, 999., 7.);
      Negv0DaugTrigger->SetDCAReCalculation(true);
      Posv0DaugTrigger->SetDCAReCalculation(true);

      v0CutsTrigger->SetPosDaugterTrackCuts(Posv0DaugTrigger);
      v0CutsTrigger->SetNegDaugterTrackCuts(Negv0DaugTrigger);
      v0CutsTrigger->SetPDGCodePosDaug(2212);  //Proton
      v0CutsTrigger->SetPDGCodeNegDaug(211);  //Pion
      v0CutsTrigger->SetPDGCodev0(3122);  //Lambda

      Antiv0CutsTrigger = new AliFemtoDreamv0Cuts();
      Antiv0CutsTrigger->SetIsMonteCarlo(isMC);
      Antiv0CutsTrigger->SetPlotCPADist(true);
      Antiv0CutsTrigger->SetPlotContrib(false);

      PosAntiv0DaugTrigger = new AliFemtoDreamTrackCuts();
      PosAntiv0DaugTrigger->SetIsMonteCarlo(isMC);
      PosAntiv0DaugTrigger->SetFillQALater(true);
      PosAntiv0DaugTrigger->SetCutCharge(1);
      PosAntiv0DaugTrigger->SetCheckPileUp(true);

      NegAntiv0DaugTrigger =new AliFemtoDreamTrackCuts();
      NegAntiv0DaugTrigger->SetIsMonteCarlo(isMC);
      NegAntiv0DaugTrigger->SetFillQALater(true);
      NegAntiv0DaugTrigger->SetCutCharge(-1);
      NegAntiv0DaugTrigger->SetCheckPileUp(true);

      Antiv0CutsTrigger->SetCutCharge(0);
      Antiv0CutsTrigger->SetPtRange(0., 999.);
      Antiv0CutsTrigger->SetCutCPA(0.95);
      Antiv0CutsTrigger->SetCutDCADaugToPrimVtx(0.02);
      Antiv0CutsTrigger->SetCutDCADaugTov0Vtx(2.5);
      Antiv0CutsTrigger->SetKaonRejection(0.49, 0.505);
      Antiv0CutsTrigger->SetCutMaxDecayVtx(120);
      Antiv0CutsTrigger->SetCutTransverseRadius(0.1, 120);
      Antiv0CutsTrigger->SetCutInvMass(0.006);
      Antiv0CutsTrigger->SetAxisInvMassPlots(400, 1.0, 1.2);

      PosAntiv0DaugTrigger->SetCutTPCCrossedRows(true, 60, 0.);
      NegAntiv0DaugTrigger->SetCutTPCCrossedRows(true, 60, 0.);
      PosAntiv0DaugTrigger->SetEtaRange(-0.9, 0.9);
      NegAntiv0DaugTrigger->SetEtaRange(-0.9, 0.9);
      PosAntiv0DaugTrigger->SetPID(AliPID::kPion, 999., 7.);
      NegAntiv0DaugTrigger->SetPID(AliPID::kProton, 999., 7.);
      PosAntiv0DaugTrigger->SetDCAReCalculation(true);
      NegAntiv0DaugTrigger->SetDCAReCalculation(true);

      Antiv0CutsTrigger->SetPosDaugterTrackCuts(PosAntiv0DaugTrigger);
      Antiv0CutsTrigger->SetNegDaugterTrackCuts(NegAntiv0DaugTrigger);
      Antiv0CutsTrigger->SetPDGCodePosDaug(211);  //Pion
      Antiv0CutsTrigger->SetPDGCodeNegDaug(2212);  //Proton
      Antiv0CutsTrigger->SetPDGCodev0(-3122);  //Lambda

      if (!fullBlastQA) {
        evtCutsTrigger->SetMinimalBooking(true);
        TrackCutsTrigger->SetMinimalBooking(true);
        AntiTrackCutsTrigger->SetMinimalBooking(true);
        v0CutsTrigger->SetMinimalBooking(true);
        Antiv0CutsTrigger->SetMinimalBooking(true);
      }
    }

    if(suffixTrigger=="11"){
      Q3Limit = 1.5;
      TriggerOnSample = false;
      //Q3LimitSample = 1.;
            evtCutsTrigger = new AliFemtoDreamEventCuts();
      evtCutsTrigger->UseDontWorryEvtCuts(false);
      evtCutsTrigger->SetZVtxPosition(-10., 10.);
      //evtCutsTrigger->CleanUpMult(false, false, false, true); will I have this in Run3 for trigger???

      TrackCutsTrigger =  new AliFemtoDreamTrackCuts();
      TrackCutsTrigger->SetPlotDCADist(true);
      TrackCutsTrigger->SetIsMonteCarlo(isMC);
      
      TrackCutsTrigger->SetFilterBit(128); //will I have this in Run3 for trigger???
      TrackCutsTrigger->SetCutCharge(1);
      TrackCutsTrigger->SetPtRange(0.5, 4.05);
      TrackCutsTrigger->SetEtaRange(-0.8, 0.8);
      TrackCutsTrigger->SetCutTPCCrossedRows(true, 70., 0.83);
      //TrackCutsTrigger->SetNClsTPC(60);// for now use only crossed pad rows
      TrackCutsTrigger->SetPID(AliPID::kProton, 0.75,3., false, 3.);  
      TrackCutsTrigger->SetDCAReCalculation(true);  //will I have this in Run3 for trigger???
      TrackCutsTrigger->SetDCAVtxZ(0.2);
      TrackCutsTrigger->SetDCAVtxXY(0.1);

      AntiTrackCutsTrigger = new AliFemtoDreamTrackCuts();
      AntiTrackCutsTrigger->SetPlotDCADist(true);
      AntiTrackCutsTrigger->SetIsMonteCarlo(isMC);

      AntiTrackCutsTrigger->SetFilterBit(128); //will I have this in Run3 for trigger???
      AntiTrackCutsTrigger->SetCutCharge(-1);
      AntiTrackCutsTrigger->SetPtRange(0.5, 4.05);
      AntiTrackCutsTrigger->SetEtaRange(-0.8, 0.8);
      AntiTrackCutsTrigger->SetCutTPCCrossedRows(true, 70., 0.83);
      //AntiTrackCutsTrigger->SetNClsTPC(60); // for now use only crossed pad rows
      AntiTrackCutsTrigger->SetPID(AliPID::kProton, 0.75,3., false, 3.); 
      AntiTrackCutsTrigger->SetDCAReCalculation(true);  //will I have this in Run3 for trigger???
      AntiTrackCutsTrigger->SetDCAVtxZ(0.2);
      AntiTrackCutsTrigger->SetDCAVtxXY(0.1);

      v0CutsTrigger = new AliFemtoDreamv0Cuts();
      v0CutsTrigger->SetIsMonteCarlo(isMC);
      v0CutsTrigger->SetPlotCPADist(true);

      Posv0DaugTrigger = new AliFemtoDreamTrackCuts(); // proton
      Posv0DaugTrigger->SetIsMonteCarlo(isMC);
      Posv0DaugTrigger->SetFillQALater(true);
      Posv0DaugTrigger->SetCutCharge(1);
      Posv0DaugTrigger->SetCheckPileUp(true);

      Negv0DaugTrigger = new AliFemtoDreamTrackCuts(); //pion
      Negv0DaugTrigger->SetIsMonteCarlo(isMC);
      Negv0DaugTrigger->SetFillQALater(true);
      Negv0DaugTrigger->SetCutCharge(-1);
      Negv0DaugTrigger->SetCheckPileUp(true);

      v0CutsTrigger->SetCutCharge(0);
      v0CutsTrigger->SetPtRange(0.3, 999.);
      v0CutsTrigger->SetCutCPA(0.99);
      v0CutsTrigger->SetCutDCADaugToPrimVtx(0.05);
      v0CutsTrigger->SetCutDCADaugTov0Vtx(1.5);
      v0CutsTrigger->SetKaonRejection(0.48, 0.515);
      v0CutsTrigger->SetCutMaxDecayVtx(100);
      v0CutsTrigger->SetCutTransverseRadius(0.2, 100);
      v0CutsTrigger->SetCutInvMass(0.004);
      v0CutsTrigger->SetAxisInvMassPlots(400, 1.0, 1.2);

      Posv0DaugTrigger->SetCutTPCCrossedRows(true, 60, 0.);
      Negv0DaugTrigger->SetCutTPCCrossedRows(true, 60, 0.);
      Posv0DaugTrigger->SetEtaRange(-0.8, 0.8);
      Negv0DaugTrigger->SetEtaRange(-0.8, 0.8);
      Posv0DaugTrigger->SetPID(AliPID::kProton, 999., 5.);
      Negv0DaugTrigger->SetPID(AliPID::kPion, 999., 5.);
      Negv0DaugTrigger->SetDCAReCalculation(true);
      Posv0DaugTrigger->SetDCAReCalculation(true);

      v0CutsTrigger->SetPosDaugterTrackCuts(Posv0DaugTrigger);
      v0CutsTrigger->SetNegDaugterTrackCuts(Negv0DaugTrigger);
      v0CutsTrigger->SetPDGCodePosDaug(2212);  //Proton
      v0CutsTrigger->SetPDGCodeNegDaug(211);  //Pion
      v0CutsTrigger->SetPDGCodev0(3122);  //Lambda

      Antiv0CutsTrigger = new AliFemtoDreamv0Cuts();
      Antiv0CutsTrigger->SetIsMonteCarlo(isMC);
      Antiv0CutsTrigger->SetPlotCPADist(true);
      Antiv0CutsTrigger->SetPlotContrib(false);

      PosAntiv0DaugTrigger = new AliFemtoDreamTrackCuts();
      PosAntiv0DaugTrigger->SetIsMonteCarlo(isMC);
      PosAntiv0DaugTrigger->SetFillQALater(true);
      PosAntiv0DaugTrigger->SetCutCharge(1);
      PosAntiv0DaugTrigger->SetCheckPileUp(true);

      NegAntiv0DaugTrigger =new AliFemtoDreamTrackCuts();
      NegAntiv0DaugTrigger->SetIsMonteCarlo(isMC);
      NegAntiv0DaugTrigger->SetFillQALater(true);
      NegAntiv0DaugTrigger->SetCutCharge(-1);
      NegAntiv0DaugTrigger->SetCheckPileUp(true);

      Antiv0CutsTrigger->SetCutCharge(0);
      Antiv0CutsTrigger->SetPtRange(0.3, 999.);
      Antiv0CutsTrigger->SetCutCPA(0.99);
      Antiv0CutsTrigger->SetCutDCADaugToPrimVtx(0.05);
      Antiv0CutsTrigger->SetCutDCADaugTov0Vtx(1.5);
      Antiv0CutsTrigger->SetKaonRejection(0.48, 0.515);
      Antiv0CutsTrigger->SetCutMaxDecayVtx(100);
      Antiv0CutsTrigger->SetCutTransverseRadius(0.2, 100);
      Antiv0CutsTrigger->SetCutInvMass(0.004);
      Antiv0CutsTrigger->SetAxisInvMassPlots(400, 1.0, 1.2);

      PosAntiv0DaugTrigger->SetCutTPCCrossedRows(true, 60, 0.);
      NegAntiv0DaugTrigger->SetCutTPCCrossedRows(true, 60, 0.);
      PosAntiv0DaugTrigger->SetEtaRange(-0.8, 0.8);
      NegAntiv0DaugTrigger->SetEtaRange(-0.8, 0.8);
      PosAntiv0DaugTrigger->SetPID(AliPID::kPion, 999., 5.);
      NegAntiv0DaugTrigger->SetPID(AliPID::kProton, 999., 5.);
      PosAntiv0DaugTrigger->SetDCAReCalculation(true);
      NegAntiv0DaugTrigger->SetDCAReCalculation(true);

      Antiv0CutsTrigger->SetPosDaugterTrackCuts(PosAntiv0DaugTrigger);
      Antiv0CutsTrigger->SetNegDaugterTrackCuts(NegAntiv0DaugTrigger);
      Antiv0CutsTrigger->SetPDGCodePosDaug(211);  //Pion
      Antiv0CutsTrigger->SetPDGCodeNegDaug(2212);  //Proton
      Antiv0CutsTrigger->SetPDGCodev0(-3122);  //Lambda

      if (!fullBlastQA) {
        evtCutsTrigger->SetMinimalBooking(true);
        TrackCutsTrigger->SetMinimalBooking(true);
        AntiTrackCutsTrigger->SetMinimalBooking(true);
        v0CutsTrigger->SetMinimalBooking(true);
        Antiv0CutsTrigger->SetMinimalBooking(true);
      }
    }
    if(suffixTrigger=="12"){
      Q3Limit = 0.6;
      TriggerOnSample = false;
      //Q3LimitSample = 1.;
            evtCutsTrigger = new AliFemtoDreamEventCuts();
      evtCutsTrigger->UseDontWorryEvtCuts(false);
      evtCutsTrigger->SetZVtxPosition(-10., 10.);
      //evtCutsTrigger->CleanUpMult(false, false, false, true); will I have this in Run3 for trigger???

      TrackCutsTrigger =  new AliFemtoDreamTrackCuts();
      TrackCutsTrigger->SetPlotDCADist(true);
      TrackCutsTrigger->SetIsMonteCarlo(isMC);
      
      TrackCutsTrigger->SetFilterBit(128); //will I have this in Run3 for trigger???
      TrackCutsTrigger->SetCutCharge(1);
      TrackCutsTrigger->SetPtRange(0.5, 4.05);
      TrackCutsTrigger->SetEtaRange(-0.8, 0.8);
      TrackCutsTrigger->SetCutTPCCrossedRows(true, 70., 0.83);
      //TrackCutsTrigger->SetNClsTPC(60);// for now use only crossed pad rows
      TrackCutsTrigger->SetPID(AliPID::kProton, 0.75,3., false, 3.);  
      TrackCutsTrigger->SetDCAReCalculation(true);  //will I have this in Run3 for trigger???
      TrackCutsTrigger->SetDCAVtxZ(0.2);
      TrackCutsTrigger->SetDCAVtxXY(0.1);

      AntiTrackCutsTrigger = new AliFemtoDreamTrackCuts();
      AntiTrackCutsTrigger->SetPlotDCADist(true);
      AntiTrackCutsTrigger->SetIsMonteCarlo(isMC);

      AntiTrackCutsTrigger->SetFilterBit(128); //will I have this in Run3 for trigger???
      AntiTrackCutsTrigger->SetCutCharge(-1);
      AntiTrackCutsTrigger->SetPtRange(0.5, 4.05);
      AntiTrackCutsTrigger->SetEtaRange(-0.8, 0.8);
      AntiTrackCutsTrigger->SetCutTPCCrossedRows(true, 70., 0.83);
      //AntiTrackCutsTrigger->SetNClsTPC(60); // for now use only crossed pad rows
      AntiTrackCutsTrigger->SetPID(AliPID::kProton, 0.75,3., false, 3.); 
      AntiTrackCutsTrigger->SetDCAReCalculation(true);  //will I have this in Run3 for trigger???
      AntiTrackCutsTrigger->SetDCAVtxZ(0.2);
      AntiTrackCutsTrigger->SetDCAVtxXY(0.1);

      v0CutsTrigger = new AliFemtoDreamv0Cuts();
      v0CutsTrigger->SetIsMonteCarlo(isMC);
      v0CutsTrigger->SetPlotCPADist(true);

      Posv0DaugTrigger = new AliFemtoDreamTrackCuts(); // proton
      Posv0DaugTrigger->SetIsMonteCarlo(isMC);
      Posv0DaugTrigger->SetFillQALater(true);
      Posv0DaugTrigger->SetCutCharge(1);
      Posv0DaugTrigger->SetCheckPileUp(true);

      Negv0DaugTrigger = new AliFemtoDreamTrackCuts(); //pion
      Negv0DaugTrigger->SetIsMonteCarlo(isMC);
      Negv0DaugTrigger->SetFillQALater(true);
      Negv0DaugTrigger->SetCutCharge(-1);
      Negv0DaugTrigger->SetCheckPileUp(true);

      v0CutsTrigger->SetCutCharge(0);
      v0CutsTrigger->SetPtRange(0.3, 999.);
      v0CutsTrigger->SetCutCPA(0.99);
      v0CutsTrigger->SetCutDCADaugToPrimVtx(0.05);
      v0CutsTrigger->SetCutDCADaugTov0Vtx(1.5);
      v0CutsTrigger->SetKaonRejection(0.48, 0.515);
      v0CutsTrigger->SetCutMaxDecayVtx(100);
      v0CutsTrigger->SetCutTransverseRadius(0.2, 100);
      v0CutsTrigger->SetCutInvMass(0.004);
      v0CutsTrigger->SetAxisInvMassPlots(400, 1.0, 1.2);

      Posv0DaugTrigger->SetCutTPCCrossedRows(true, 60, 0.);
      Negv0DaugTrigger->SetCutTPCCrossedRows(true, 60, 0.);
      Posv0DaugTrigger->SetEtaRange(-0.8, 0.8);
      Negv0DaugTrigger->SetEtaRange(-0.8, 0.8);
      Posv0DaugTrigger->SetPID(AliPID::kProton, 999., 5.);
      Negv0DaugTrigger->SetPID(AliPID::kPion, 999., 5.);
      Negv0DaugTrigger->SetDCAReCalculation(true);
      Posv0DaugTrigger->SetDCAReCalculation(true);

      v0CutsTrigger->SetPosDaugterTrackCuts(Posv0DaugTrigger);
      v0CutsTrigger->SetNegDaugterTrackCuts(Negv0DaugTrigger);
      v0CutsTrigger->SetPDGCodePosDaug(2212);  //Proton
      v0CutsTrigger->SetPDGCodeNegDaug(211);  //Pion
      v0CutsTrigger->SetPDGCodev0(3122);  //Lambda

      Antiv0CutsTrigger = new AliFemtoDreamv0Cuts();
      Antiv0CutsTrigger->SetIsMonteCarlo(isMC);
      Antiv0CutsTrigger->SetPlotCPADist(true);
      Antiv0CutsTrigger->SetPlotContrib(false);

      PosAntiv0DaugTrigger = new AliFemtoDreamTrackCuts();
      PosAntiv0DaugTrigger->SetIsMonteCarlo(isMC);
      PosAntiv0DaugTrigger->SetFillQALater(true);
      PosAntiv0DaugTrigger->SetCutCharge(1);
      PosAntiv0DaugTrigger->SetCheckPileUp(true);

      NegAntiv0DaugTrigger =new AliFemtoDreamTrackCuts();
      NegAntiv0DaugTrigger->SetIsMonteCarlo(isMC);
      NegAntiv0DaugTrigger->SetFillQALater(true);
      NegAntiv0DaugTrigger->SetCutCharge(-1);
      NegAntiv0DaugTrigger->SetCheckPileUp(true);

      Antiv0CutsTrigger->SetCutCharge(0);
      Antiv0CutsTrigger->SetPtRange(0.3, 999.);
      Antiv0CutsTrigger->SetCutCPA(0.99);
      Antiv0CutsTrigger->SetCutDCADaugToPrimVtx(0.05);
      Antiv0CutsTrigger->SetCutDCADaugTov0Vtx(1.5);
      Antiv0CutsTrigger->SetKaonRejection(0.48, 0.515);
      Antiv0CutsTrigger->SetCutMaxDecayVtx(100);
      Antiv0CutsTrigger->SetCutTransverseRadius(0.2, 100);
      Antiv0CutsTrigger->SetCutInvMass(0.004);
      Antiv0CutsTrigger->SetAxisInvMassPlots(400, 1.0, 1.2);

      PosAntiv0DaugTrigger->SetCutTPCCrossedRows(true, 60, 0.);
      NegAntiv0DaugTrigger->SetCutTPCCrossedRows(true, 60, 0.);
      PosAntiv0DaugTrigger->SetEtaRange(-0.8, 0.8);
      NegAntiv0DaugTrigger->SetEtaRange(-0.8, 0.8);
      PosAntiv0DaugTrigger->SetPID(AliPID::kPion, 999., 5.);
      NegAntiv0DaugTrigger->SetPID(AliPID::kProton, 999., 5.);
      PosAntiv0DaugTrigger->SetDCAReCalculation(true);
      NegAntiv0DaugTrigger->SetDCAReCalculation(true);

      Antiv0CutsTrigger->SetPosDaugterTrackCuts(PosAntiv0DaugTrigger);
      Antiv0CutsTrigger->SetNegDaugterTrackCuts(NegAntiv0DaugTrigger);
      Antiv0CutsTrigger->SetPDGCodePosDaug(211);  //Pion
      Antiv0CutsTrigger->SetPDGCodeNegDaug(2212);  //Proton
      Antiv0CutsTrigger->SetPDGCodev0(-3122);  //Lambda

      if (!fullBlastQA) {
        evtCutsTrigger->SetMinimalBooking(true);
        TrackCutsTrigger->SetMinimalBooking(true);
        AntiTrackCutsTrigger->SetMinimalBooking(true);
        v0CutsTrigger->SetMinimalBooking(true);
        Antiv0CutsTrigger->SetMinimalBooking(true);
      }
    }

    if(suffixTrigger=="13"){


      TriggerOnSample = true;

      //Q3LimitSample = 1.;
      evtCutsTrigger = new AliFemtoDreamEventCuts();
      evtCutsTrigger->UseDontWorryEvtCuts(false);
      evtCutsTrigger->SetZVtxPosition(-10., 10.);
      //evtCutsTrigger->CleanUpMult(false, false, false, true); will I have this in Run3 for trigger???

      TrackCutsTrigger =  new AliFemtoDreamTrackCuts();
      TrackCutsTrigger->SetPlotDCADist(true);
      TrackCutsTrigger->SetIsMonteCarlo(isMC);
      
      TrackCutsTrigger->SetFilterBit(128); //will I have this in Run3 for trigger???
      TrackCutsTrigger->SetCutCharge(1);
      TrackCutsTrigger->SetPtRange(0.5, 4.05);
      TrackCutsTrigger->SetEtaRange(-0.8, 0.8);
      TrackCutsTrigger->SetCutTPCCrossedRows(true, 70., 0.83);
      //TrackCutsTrigger->SetNClsTPC(60);// for now use only crossed pad rows
      TrackCutsTrigger->SetPID(AliPID::kProton, 0.75,3., false, 3.);  
      TrackCutsTrigger->SetDCAReCalculation(true);  //will I have this in Run3 for trigger???
      TrackCutsTrigger->SetDCAVtxZ(0.2);
      TrackCutsTrigger->SetDCAVtxXY(0.1);

      AntiTrackCutsTrigger = new AliFemtoDreamTrackCuts();
      AntiTrackCutsTrigger->SetPlotDCADist(true);
      AntiTrackCutsTrigger->SetIsMonteCarlo(isMC);

      AntiTrackCutsTrigger->SetFilterBit(128); //will I have this in Run3 for trigger???
      AntiTrackCutsTrigger->SetCutCharge(-1);
      AntiTrackCutsTrigger->SetPtRange(0.5, 4.05);
      AntiTrackCutsTrigger->SetEtaRange(-0.8, 0.8);
      AntiTrackCutsTrigger->SetCutTPCCrossedRows(true, 70., 0.83);
      //AntiTrackCutsTrigger->SetNClsTPC(60); // for now use only crossed pad rows
      AntiTrackCutsTrigger->SetPID(AliPID::kProton, 0.75,3., false, 3.); 
      AntiTrackCutsTrigger->SetDCAReCalculation(true);  //will I have this in Run3 for trigger???
      AntiTrackCutsTrigger->SetDCAVtxZ(0.2);
      AntiTrackCutsTrigger->SetDCAVtxXY(0.1);

      v0CutsTrigger = new AliFemtoDreamv0Cuts();
      v0CutsTrigger->SetIsMonteCarlo(isMC);
      v0CutsTrigger->SetPlotCPADist(true);

      Posv0DaugTrigger = new AliFemtoDreamTrackCuts(); // proton
      Posv0DaugTrigger->SetIsMonteCarlo(isMC);
      Posv0DaugTrigger->SetFillQALater(true);
      Posv0DaugTrigger->SetCutCharge(1);
      Posv0DaugTrigger->SetCheckPileUp(true);

      Negv0DaugTrigger = new AliFemtoDreamTrackCuts(); //pion
      Negv0DaugTrigger->SetIsMonteCarlo(isMC);
      Negv0DaugTrigger->SetFillQALater(true);
      Negv0DaugTrigger->SetCutCharge(-1);
      Negv0DaugTrigger->SetCheckPileUp(true);

      v0CutsTrigger->SetCutCharge(0);
      v0CutsTrigger->SetPtRange(0.3, 999.);
      v0CutsTrigger->SetCutCPA(0.99);
      v0CutsTrigger->SetCutDCADaugToPrimVtx(0.05);
      v0CutsTrigger->SetCutDCADaugTov0Vtx(1.5);
      v0CutsTrigger->SetKaonRejection(0.48, 0.515);
      v0CutsTrigger->SetCutMaxDecayVtx(100);
      v0CutsTrigger->SetCutTransverseRadius(0.2, 100);
      v0CutsTrigger->SetCutInvMass(0.004);
      v0CutsTrigger->SetAxisInvMassPlots(400, 1.0, 1.2);

      Posv0DaugTrigger->SetCutTPCCrossedRows(true, 60, 0.);
      Negv0DaugTrigger->SetCutTPCCrossedRows(true, 60, 0.);
      Posv0DaugTrigger->SetEtaRange(-0.8, 0.8);
      Negv0DaugTrigger->SetEtaRange(-0.8, 0.8);
      Posv0DaugTrigger->SetPID(AliPID::kProton, 999., 5.);
      Negv0DaugTrigger->SetPID(AliPID::kPion, 999., 5.);
      Negv0DaugTrigger->SetDCAReCalculation(true);
      Posv0DaugTrigger->SetDCAReCalculation(true);

      v0CutsTrigger->SetPosDaugterTrackCuts(Posv0DaugTrigger);
      v0CutsTrigger->SetNegDaugterTrackCuts(Negv0DaugTrigger);
      v0CutsTrigger->SetPDGCodePosDaug(2212);  //Proton
      v0CutsTrigger->SetPDGCodeNegDaug(211);  //Pion
      v0CutsTrigger->SetPDGCodev0(3122);  //Lambda

      Antiv0CutsTrigger = new AliFemtoDreamv0Cuts();
      Antiv0CutsTrigger->SetIsMonteCarlo(isMC);
      Antiv0CutsTrigger->SetPlotCPADist(true);
      Antiv0CutsTrigger->SetPlotContrib(false);

      PosAntiv0DaugTrigger = new AliFemtoDreamTrackCuts();
      PosAntiv0DaugTrigger->SetIsMonteCarlo(isMC);
      PosAntiv0DaugTrigger->SetFillQALater(true);
      PosAntiv0DaugTrigger->SetCutCharge(1);
      PosAntiv0DaugTrigger->SetCheckPileUp(true);

      NegAntiv0DaugTrigger =new AliFemtoDreamTrackCuts();
      NegAntiv0DaugTrigger->SetIsMonteCarlo(isMC);
      NegAntiv0DaugTrigger->SetFillQALater(true);
      NegAntiv0DaugTrigger->SetCutCharge(-1);
      NegAntiv0DaugTrigger->SetCheckPileUp(true);

      Antiv0CutsTrigger->SetCutCharge(0);
      Antiv0CutsTrigger->SetPtRange(0.3, 999.);
      Antiv0CutsTrigger->SetCutCPA(0.99);
      Antiv0CutsTrigger->SetCutDCADaugToPrimVtx(0.05);
      Antiv0CutsTrigger->SetCutDCADaugTov0Vtx(1.5);
      Antiv0CutsTrigger->SetKaonRejection(0.48, 0.515);
      Antiv0CutsTrigger->SetCutMaxDecayVtx(100);
      Antiv0CutsTrigger->SetCutTransverseRadius(0.2, 100);
      Antiv0CutsTrigger->SetCutInvMass(0.004);
      Antiv0CutsTrigger->SetAxisInvMassPlots(400, 1.0, 1.2);

      PosAntiv0DaugTrigger->SetCutTPCCrossedRows(true, 60, 0.);
      NegAntiv0DaugTrigger->SetCutTPCCrossedRows(true, 60, 0.);
      PosAntiv0DaugTrigger->SetEtaRange(-0.8, 0.8);
      NegAntiv0DaugTrigger->SetEtaRange(-0.8, 0.8);
      PosAntiv0DaugTrigger->SetPID(AliPID::kPion, 999., 5.);
      NegAntiv0DaugTrigger->SetPID(AliPID::kProton, 999., 5.);
      PosAntiv0DaugTrigger->SetDCAReCalculation(true);
      NegAntiv0DaugTrigger->SetDCAReCalculation(true);

      Antiv0CutsTrigger->SetPosDaugterTrackCuts(PosAntiv0DaugTrigger);
      Antiv0CutsTrigger->SetNegDaugterTrackCuts(NegAntiv0DaugTrigger);
      Antiv0CutsTrigger->SetPDGCodePosDaug(211);  //Pion
      Antiv0CutsTrigger->SetPDGCodeNegDaug(2212);  //Proton
      Antiv0CutsTrigger->SetPDGCodev0(-3122);  //Lambda

      if (!fullBlastQA) {
        evtCutsTrigger->SetMinimalBooking(true);
        TrackCutsTrigger->SetMinimalBooking(true);
        AntiTrackCutsTrigger->SetMinimalBooking(true);
        v0CutsTrigger->SetMinimalBooking(true);
        Antiv0CutsTrigger->SetMinimalBooking(true);
      }
    }


    if(suffixTrigger=="14"){
      evtCutsTrigger = new AliFemtoDreamEventCuts();
      evtCutsTrigger->UseDontWorryEvtCuts(false);
      evtCutsTrigger->SetZVtxPosition(-12., 12.);
      //evtCutsTrigger->CleanUpMult(false, false, false, true); will I have this in Run3 for trigger???

      TrackCutsTrigger =  new AliFemtoDreamTrackCuts();
      TrackCutsTrigger->SetPlotDCADist(true);
      TrackCutsTrigger->SetIsMonteCarlo(isMC);
      
      TrackCutsTrigger->SetFilterBit(128); //will I have this in Run3 for trigger???
      TrackCutsTrigger->SetCutCharge(1);
      TrackCutsTrigger->SetPtRange(0.35, 5.);
      TrackCutsTrigger->SetEtaRange(-0.85, 0.85);
      TrackCutsTrigger->SetCutTPCCrossedRows(true, 70., 0.83);
      TrackCutsTrigger->SetNClsTPC(65);
      //TrackCutsTrigger->SetNClsTPC(60);// for now use only crossed pad rows
      TrackCutsTrigger->SetPID(AliPID::kProton, 0.75,4., false, 4.);  
      TrackCutsTrigger->SetDCAReCalculation(true);  //will I have this in Run3 for trigger???
      TrackCutsTrigger->SetDCAVtxZ(0.3);
      TrackCutsTrigger->SetDCAVtxXY(0.15);
      TrackCutsTrigger->SetRejLowPtPionsTOF(true);
      TrackCutsTrigger->SetCutSmallestSig(true);

      AntiTrackCutsTrigger = new AliFemtoDreamTrackCuts();
      AntiTrackCutsTrigger->SetPlotDCADist(true);
      AntiTrackCutsTrigger->SetIsMonteCarlo(isMC);

      AntiTrackCutsTrigger->SetFilterBit(128); //will I have this in Run3 for trigger???
      AntiTrackCutsTrigger->SetCutCharge(-1);
      AntiTrackCutsTrigger->SetPtRange(0.35, 5.);
      AntiTrackCutsTrigger->SetEtaRange(-0.85, 0.85);
      AntiTrackCutsTrigger->SetCutTPCCrossedRows(true, 70., 0.83);
      AntiTrackCutsTrigger->SetNClsTPC(65);
      //AntiTrackCutsTrigger->SetNClsTPC(60); // for now use only crossed pad rows
      AntiTrackCutsTrigger->SetPID(AliPID::kProton, 0.75,4., false, 4.); 
      AntiTrackCutsTrigger->SetDCAReCalculation(true);  //will I have this in Run3 for trigger???
      AntiTrackCutsTrigger->SetDCAVtxZ(0.3);
      AntiTrackCutsTrigger->SetDCAVtxXY(0.15);
      AntiTrackCutsTrigger->SetRejLowPtPionsTOF(true);
      AntiTrackCutsTrigger->SetCutSmallestSig(true);

      v0CutsTrigger = new AliFemtoDreamv0Cuts();
      v0CutsTrigger->SetIsMonteCarlo(isMC);
      v0CutsTrigger->SetPlotCPADist(true);

      Posv0DaugTrigger = new AliFemtoDreamTrackCuts(); // proton
      Posv0DaugTrigger->SetIsMonteCarlo(isMC);
      Posv0DaugTrigger->SetFillQALater(true);
      Posv0DaugTrigger->SetCutCharge(1);
      Posv0DaugTrigger->SetCheckPileUp(true);

      Negv0DaugTrigger = new AliFemtoDreamTrackCuts(); //pion
      Negv0DaugTrigger->SetIsMonteCarlo(isMC);
      Negv0DaugTrigger->SetFillQALater(true);
      Negv0DaugTrigger->SetCutCharge(-1);
      Negv0DaugTrigger->SetCheckPileUp(true);

      v0CutsTrigger->SetCutCharge(0);
      v0CutsTrigger->SetPtRange(0., 999.);
      v0CutsTrigger->SetCutCPA(0.985);
      v0CutsTrigger->SetCutDCADaugToPrimVtx(0.04);
      v0CutsTrigger->SetCutDCADaugTov0Vtx(1.8);
      v0CutsTrigger->SetKaonRejection(0.49, 0.505);
      v0CutsTrigger->SetCutMaxDecayVtx(100);
      v0CutsTrigger->SetCutTransverseRadius(0.2, 100);
      v0CutsTrigger->SetCutInvMass(0.006);
      v0CutsTrigger->SetAxisInvMassPlots(400, 1.0, 1.2);
      Posv0DaugTrigger->SetNClsTPC(60);
      Negv0DaugTrigger->SetNClsTPC(60);
      Posv0DaugTrigger->SetEtaRange(-0.85, 0.85);
      Negv0DaugTrigger->SetEtaRange(-0.85, 0.85);
      Posv0DaugTrigger->SetPID(AliPID::kProton, 999., 6.);
      Negv0DaugTrigger->SetPID(AliPID::kPion, 999., 6.);
      Negv0DaugTrigger->SetDCAReCalculation(true);
      Posv0DaugTrigger->SetDCAReCalculation(true);

      v0CutsTrigger->SetPosDaugterTrackCuts(Posv0DaugTrigger);
      v0CutsTrigger->SetNegDaugterTrackCuts(Negv0DaugTrigger);
      v0CutsTrigger->SetPDGCodePosDaug(2212);  //Proton
      v0CutsTrigger->SetPDGCodeNegDaug(211);  //Pion
      v0CutsTrigger->SetPDGCodev0(3122);  //Lambda

      Antiv0CutsTrigger = new AliFemtoDreamv0Cuts();
      Antiv0CutsTrigger->SetIsMonteCarlo(isMC);
      Antiv0CutsTrigger->SetPlotCPADist(true);
      Antiv0CutsTrigger->SetPlotContrib(false);

      PosAntiv0DaugTrigger = new AliFemtoDreamTrackCuts();
      PosAntiv0DaugTrigger->SetIsMonteCarlo(isMC);
      PosAntiv0DaugTrigger->SetFillQALater(true);
      PosAntiv0DaugTrigger->SetCutCharge(1);
      PosAntiv0DaugTrigger->SetCheckPileUp(true);

      NegAntiv0DaugTrigger =new AliFemtoDreamTrackCuts();
      NegAntiv0DaugTrigger->SetIsMonteCarlo(isMC);
      NegAntiv0DaugTrigger->SetFillQALater(true);
      NegAntiv0DaugTrigger->SetCutCharge(-1);
      NegAntiv0DaugTrigger->SetCheckPileUp(true);

      Antiv0CutsTrigger->SetCutCharge(0);
      Antiv0CutsTrigger->SetPtRange(0., 999.);
      Antiv0CutsTrigger->SetCutCPA(0.985);
      Antiv0CutsTrigger->SetCutDCADaugToPrimVtx(0.04);
      Antiv0CutsTrigger->SetCutDCADaugTov0Vtx(1.8);
      Antiv0CutsTrigger->SetKaonRejection(0.49, 0.505);
      Antiv0CutsTrigger->SetCutMaxDecayVtx(100);
      Antiv0CutsTrigger->SetCutTransverseRadius(0.2, 100);
      Antiv0CutsTrigger->SetCutInvMass(0.006);
      Antiv0CutsTrigger->SetAxisInvMassPlots(400, 1.0, 1.2);
      PosAntiv0DaugTrigger->SetNClsTPC(60);
      NegAntiv0DaugTrigger->SetNClsTPC(60);
      PosAntiv0DaugTrigger->SetEtaRange(-0.85, 0.85);
      NegAntiv0DaugTrigger->SetEtaRange(-0.85, 0.85);
      PosAntiv0DaugTrigger->SetPID(AliPID::kPion, 999., 6.);
      NegAntiv0DaugTrigger->SetPID(AliPID::kProton, 999., 6.);
      PosAntiv0DaugTrigger->SetDCAReCalculation(true);
      NegAntiv0DaugTrigger->SetDCAReCalculation(true);

      Antiv0CutsTrigger->SetPosDaugterTrackCuts(PosAntiv0DaugTrigger);
      Antiv0CutsTrigger->SetNegDaugterTrackCuts(NegAntiv0DaugTrigger);
      Antiv0CutsTrigger->SetPDGCodePosDaug(211);  //Pion
      Antiv0CutsTrigger->SetPDGCodeNegDaug(2212);  //Proton
      Antiv0CutsTrigger->SetPDGCodev0(-3122);  //Lambda

      if (!fullBlastQA) {
        evtCutsTrigger->SetMinimalBooking(true);
        TrackCutsTrigger->SetMinimalBooking(true);
        AntiTrackCutsTrigger->SetMinimalBooking(true);
        v0CutsTrigger->SetMinimalBooking(true);
        Antiv0CutsTrigger->SetMinimalBooking(true);
      }
    }




  }
  AliFemtoDreamCollConfig *config = new AliFemtoDreamCollConfig("Femto",
                                                                "Femto", false);
  // Femto Collection
  std::vector<int> PDGParticles;
  PDGParticles.push_back(2212);
  PDGParticles.push_back(2212);
  PDGParticles.push_back(3122);
  PDGParticles.push_back(3122);

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
  config->SetDeltaEtaMax(0.012);
  config->SetDeltaPhiMax(0.012);

  if(suffixTrigger=="699"){
    config->SetDeltaEtaMax(0.02);
    config->SetDeltaPhiMax(0.02);
  }  
  if(suffixTrigger=="669"){
    config->SetDeltaEtaMax(0.03);
    config->SetDeltaPhiMax(0.03);
  }
  config->SetExtendedQAPairs(pairQA);
  config->SetMixingDepth(mixingDepthFromTask);
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

  TString addon = "PL";
  TString file = AliAnalysisManager::GetCommonFileName();

  TString EvtCutsName = Form("%sEvtCuts_%s_%s", addon.Data(), suffixTrigger.Data(), suffix.Data());
  AliAnalysisDataContainer *coutputEvtCuts = mgr->CreateContainer(
      EvtCutsName.Data(), TList::Class(), AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), EvtCutsName.Data()));

  TString TrackCutsName = Form("%sTrackCuts_%s_%s", addon.Data(), suffixTrigger.Data(), suffix.Data());
  AliAnalysisDataContainer *couputTrkCuts = mgr->CreateContainer(
      TrackCutsName.Data(), TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), TrackCutsName.Data()));

  TString AntiTrackCutsName = Form("%sAntiTrackCuts_%s_%s", addon.Data(), suffixTrigger.Data(), suffix.Data());
  AliAnalysisDataContainer *coutputAntiTrkCuts = mgr->CreateContainer(
      AntiTrackCutsName.Data(), TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), AntiTrackCutsName.Data()));

  AliAnalysisDataContainer *coutputv0Cuts;
  TString v0CutsName = Form("%sv0Cuts_%s_%s", addon.Data(), suffixTrigger.Data(), suffix.Data());
  coutputv0Cuts = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      v0CutsName.Data(),
      TList::Class(), AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), v0CutsName.Data()));

  AliAnalysisDataContainer *coutputAntiv0Cuts;
  TString Antiv0CutsName = Form("%sAntiv0Cuts_%s_%s", addon.Data(), suffixTrigger.Data(), suffix.Data());
  coutputAntiv0Cuts = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      Antiv0CutsName.Data(),
      TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), Antiv0CutsName.Data()));

  AliAnalysisDataContainer *coutputResults;
  TString ResultsName = Form("%sResults_%s_%s", addon.Data(), suffixTrigger.Data(), suffix.Data());
  coutputResults = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      ResultsName.Data(),
      TList::Class(), AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), ResultsName.Data()));

  AliAnalysisDataContainer *coutputResultsQA;
  TString ResultsQAName = Form("%sResultsQA_%s_%s", addon.Data(), suffixTrigger.Data(), suffix.Data());
  coutputResultsQA = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      ResultsQAName.Data(),
      TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), ResultsQAName.Data()));

  AliAnalysisDataContainer *coutputResultsSample;
  TString ResultsSampleName = Form("%sResultsSample_%s_%s", addon.Data(), suffixTrigger.Data(), suffix.Data());
  coutputResultsSample = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      ResultsSampleName.Data(),
      TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), ResultsSampleName.Data()));

  AliAnalysisDataContainer *coutputResultsSampleQA;
  TString ResultsSampleQAName = Form("%sResultsSampleQA_%s_%s", addon.Data(), suffixTrigger.Data(), suffix.Data());
  coutputResultsSampleQA = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      ResultsSampleQAName.Data(),
      TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), ResultsSampleQAName.Data()));

  AliAnalysisDataContainer *coutputTrkCutsMC;
  AliAnalysisDataContainer *coutputAntiTrkCutsMC;
  AliAnalysisDataContainer *coutputv0CutsMC;
  AliAnalysisDataContainer *coutputAntiv0CutsMC;
  if (isMC) {
    TString TrkCutsMCName = Form("%sTrkCutsMC_%s_%s", addon.Data(), suffixTrigger.Data(), suffix.Data());
    coutputTrkCutsMC = mgr->CreateContainer(
        //@suppress("Invalid arguments") it works ffs
        TrkCutsMCName.Data(),
        TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), TrkCutsMCName.Data()));

    TString AntiTrkCutsMCName = Form("%sAntiTrkCutsMC_%s_%s", addon.Data(), suffixTrigger.Data(), suffix.Data());
    coutputAntiTrkCutsMC = mgr->CreateContainer(
        //@suppress("Invalid arguments") it works ffs
        AntiTrkCutsMCName.Data(),
        TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), AntiTrkCutsMCName.Data()));

    TString v0CutsMCName = Form("%sv0CutsMC_%s_%s", addon.Data(), suffixTrigger.Data(), suffix.Data());
    coutputv0CutsMC = mgr->CreateContainer(
        //@suppress("Invalid arguments") it works ffs
        v0CutsMCName.Data(),
        TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), v0CutsMCName.Data()));

    TString Antiv0CutsMCName = Form("%sAntiv0CutsMC_%s_%s", addon.Data(), suffixTrigger.Data(), suffix.Data());
    coutputAntiv0CutsMC = mgr->CreateContainer(
        //@suppress("Invalid arguments") it works ffs
        Antiv0CutsMCName.Data(),
        TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), Antiv0CutsMCName.Data()));

  }

  TString ThreeBodyName = Form("%sThreeBody_%s_%s", addon.Data(), suffixTrigger.Data(), suffix.Data());
  AliAnalysisDataContainer *coutputThreeBody = mgr->CreateContainer(
    //@suppress("Invalid arguments") it works ffs
    ThreeBodyName.Data(),
    TList::Class(),
    AliAnalysisManager::kOutputContainer,
    Form("%s:%s", file.Data(), ThreeBodyName.Data()));

  AliAnalysisDataContainer *coutputEvtCutsTrigger;
  AliAnalysisDataContainer *couputTrkCutsTrigger;
  AliAnalysisDataContainer *coutputAntiTrkCutsTrigger;
  AliAnalysisDataContainer *coutputv0CutsTrigger;
  AliAnalysisDataContainer *coutputAntiv0CutsTrigger;
  if(triggerOn){
    TString EvtCutsTriggerName = Form("%sEvtCutsTrigger_%s_%s", addon.Data(), suffixTrigger.Data(), suffix.Data());
    coutputEvtCutsTrigger = mgr->CreateContainer(
        EvtCutsTriggerName.Data(), TList::Class(), AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), EvtCutsTriggerName.Data()));

    TString TrackCutsTriggerName = Form("%sTrackCutsTrigger_%s_%s", addon.Data(), suffixTrigger.Data(), suffix.Data());
    couputTrkCutsTrigger = mgr->CreateContainer(
        TrackCutsTriggerName.Data(), TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), TrackCutsTriggerName.Data()));

    TString AntiTrackCutsTriggerName = Form("%sAntiTrackCutsTrigger_%s_%s", addon.Data(), suffixTrigger.Data(), suffix.Data());
    coutputAntiTrkCutsTrigger = mgr->CreateContainer(
        AntiTrackCutsTriggerName.Data(), TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), AntiTrackCutsTriggerName.Data()));

    
    TString v0CutsTriggerName = Form("%sv0CutsTrigger_%s_%s", addon.Data(), suffixTrigger.Data(), suffix.Data());
    coutputv0CutsTrigger = mgr->CreateContainer(
        //@suppress("Invalid arguments") it works ffs
        v0CutsTriggerName.Data(),
        TList::Class(), AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), v0CutsTriggerName.Data()));

    
    TString Antiv0CutsTriggerName = Form("%sAntiv0CutsTrigger_%s_%s", addon.Data(), suffixTrigger.Data(), suffix.Data());
    coutputAntiv0CutsTrigger = mgr->CreateContainer(
        //@suppress("Invalid arguments") it works ffs
        Antiv0CutsTriggerName.Data(),
        TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), Antiv0CutsTriggerName.Data()));
  }


  AliAnalysisTaskThreeBodyFemto* taskNano;
  AliAnalysisTaskThreeBodyFemtoAOD* taskAOD;
  if(isNano){
    taskNano= new AliAnalysisTaskThreeBodyFemto("femtoNanoThreeBody", isMC);
    if (!fullBlastQA)
    { 
      taskNano->SetRunTaskLightWeight(true);
    }

    if (trigger == 0) { 
        taskNano->SelectCollisionCandidates(AliVEvent::kHighMultV0);  
      } else if (trigger == 1){     
        taskNano->SelectCollisionCandidates(AliVEvent::kINT7);  
      } 
    taskNano->SetEventCuts(evtCuts);  
    taskNano->SetProtonCuts(TrackCuts); 
    taskNano->SetAntiProtonCuts(AntiTrackCuts); 
    taskNano->Setv0Cuts(v0Cuts);  
    taskNano->SetAntiv0Cuts(Antiv0Cuts);  
    taskNano->SetCorrelationConfig(config); 
    taskNano->SetRunThreeBodyHistograms(true);
    taskNano->SetClosePairRejectionForAll(ClosePairRejectionForAll);
    taskNano->SetRun2Body(run2Body);
    taskNano->SetMixingChoice(mixinfChoice);
    taskNano->SetSame2Mixed1Choice(mix21);
    taskNano->SetWhichTripletsToRun(whichTripletsToRun);
    
    

    mgr->AddTask(taskNano); 
    
    mgr->ConnectInput(taskNano, 0, cinput); 
    mgr->ConnectOutput(taskNano, 1, coutputEvtCuts);  
    mgr->ConnectOutput(taskNano, 2, couputTrkCuts); 
    mgr->ConnectOutput(taskNano, 3, coutputAntiTrkCuts);  
    mgr->ConnectOutput(taskNano, 4, coutputv0Cuts); 
    mgr->ConnectOutput(taskNano, 5, coutputAntiv0Cuts); 
    mgr->ConnectOutput(taskNano, 6, coutputResults);  
    mgr->ConnectOutput(taskNano, 7, coutputResultsQA);  
    mgr->ConnectOutput(taskNano, 8, coutputResultsSample);  
    mgr->ConnectOutput(taskNano, 9, coutputResultsSampleQA);  
    mgr->ConnectOutput(taskNano, 10, coutputThreeBody);  
    if (isMC) { 
      mgr->ConnectOutput(taskNano, 11, coutputTrkCutsMC); 
      mgr->ConnectOutput(taskNano, 12, coutputAntiTrkCutsMC); 
      mgr->ConnectOutput(taskNano, 13, coutputv0CutsMC);  
      mgr->ConnectOutput(taskNano, 14, coutputAntiv0CutsMC);  
    } 
  }
  else{
    taskAOD= new AliAnalysisTaskThreeBodyFemtoAOD("femtoAODThreeBody", isMC, triggerOn);
    if (!fullBlastQA)
    { 
      taskAOD->SetRunTaskLightWeight(true);
    }

    if (trigger == 0) { 
        taskAOD->SelectCollisionCandidates(AliVEvent::kHighMultV0);  
      } else if (trigger == 1){     
        taskAOD->SelectCollisionCandidates(AliVEvent::kINT7);  
      } 
    taskAOD->SetEventCuts(evtCuts);  
    taskAOD->SetProtonCuts(TrackCuts); 
    taskAOD->SetAntiProtonCuts(AntiTrackCuts); 
    taskAOD->Setv0Cuts(v0Cuts);  
    taskAOD->SetAntiv0Cuts(Antiv0Cuts);   
    if(triggerOn){
      taskAOD->SetEventCutsTrigger(evtCutsTrigger);  
      taskAOD->SetProtonCutsTrigger(TrackCutsTrigger); 
      taskAOD->SetAntiProtonCutsTrigger(AntiTrackCutsTrigger); 
      taskAOD->Setv0CutsTrigger(v0CutsTrigger);  
      taskAOD->SetAntiv0CutsTrigger(Antiv0CutsTrigger);  
      taskAOD->SetQ3Limit(Q3Limit);
      taskAOD->SetRunppp0ppL1(Runppp0ppL1);
      if(TriggerOnSample){
        taskAOD->SetQ3LimitSample(Q3LimitSample);
        taskAOD->SetTriggerOnSample(TriggerOnSample);
      }
    }
    taskAOD->SetCorrelationConfig(config); 
    taskAOD->SetRunThreeBodyHistograms(true);
    taskAOD->SetTriggerOn(triggerOn);
    taskAOD->SetIsMC(isMC);




    taskAOD->SetQ3Limit(Q3Limit);
    taskAOD->SetQ3LimitSample(Q3LimitSample) ;
    taskAOD->SetQ3LimitSample2(Q3LimitSample2) ;
    taskAOD->SetQ3LimitSampleFraction( Q3LimitSampleFraction) ;
    taskAOD->SetQ3LimitSampleFraction2( Q3LimitSampleFraction2) ;
    taskAOD->SetQ3LimitFraction( Q3LimitFraction) ;


    mgr->AddTask(taskAOD); 
    
    mgr->ConnectInput(taskAOD, 0, cinput); 
    mgr->ConnectOutput(taskAOD, 1, coutputEvtCuts);  
    mgr->ConnectOutput(taskAOD, 2, couputTrkCuts); 
    mgr->ConnectOutput(taskAOD, 3, coutputAntiTrkCuts);  
    mgr->ConnectOutput(taskAOD, 4, coutputv0Cuts); 
    mgr->ConnectOutput(taskAOD, 5, coutputAntiv0Cuts); 
    mgr->ConnectOutput(taskAOD, 6, coutputResults);  
    mgr->ConnectOutput(taskAOD, 7, coutputResultsQA);  
    mgr->ConnectOutput(taskAOD, 8, coutputResultsSample);  
    mgr->ConnectOutput(taskAOD, 9, coutputResultsSampleQA); 
    mgr->ConnectOutput(taskAOD, 10, coutputThreeBody); 
    if(triggerOn){
      mgr->ConnectOutput(taskAOD, 11, coutputEvtCutsTrigger);  
      mgr->ConnectOutput(taskAOD, 12, couputTrkCutsTrigger); 
      mgr->ConnectOutput(taskAOD, 13, coutputAntiTrkCutsTrigger);  
      mgr->ConnectOutput(taskAOD, 14, coutputv0CutsTrigger); 
      mgr->ConnectOutput(taskAOD, 15, coutputAntiv0CutsTrigger); 
    }
    if (isMC) { 
      mgr->ConnectOutput(taskAOD, 16, coutputTrkCutsMC); 
      mgr->ConnectOutput(taskAOD, 17, coutputAntiTrkCutsMC); 
      mgr->ConnectOutput(taskAOD, 18, coutputv0CutsMC);  
      mgr->ConnectOutput(taskAOD, 19, coutputAntiv0CutsMC);  
    } 
  }

      
  if (isNano) {
    return taskNano;
  } else {
    return taskAOD;
  }

  
}
