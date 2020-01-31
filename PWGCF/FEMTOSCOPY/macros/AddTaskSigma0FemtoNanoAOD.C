#ifndef __CINT__
#include <vector>
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliFemtoDreamEventCuts.h"
#include "AliFemtoDreamTrackCuts.h"
#include "AliFemtoDreamv0Cuts.h"
#include "AliSigma0PhotonCuts.h"
#include "AliSigma0AODPhotonMotherCuts.h"
#include "AliFemtoDreamCollConfig.h"
#include "AliAnalysisTaskNanoAODSigma0Femto.h"
#endif

AliAnalysisTaskSE *AddTaskSigma0FemtoNanoAOD(bool isMC = false, bool MomRes =
                                                 false,
                                             bool fullBlastQA = false,
                                             TString trigger = "kINT7",
                                             const char *cutVariation = "0") {
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

  // =====================================================================
  // Proton cuts
  const float ProtonPtlow = 0.4;
  const float ProtonPtup = 0.6;
  const float ProtonEtaLow = 0.75;
  const float ProtonEtaUp = 0.85;
  const float ProtonNsigmaLow = 2.5;
  const float ProtonNsigmaUp = 3.5;
  const float ProtonNClsLow = 70;
  const float ProtonNClsUp = 90;

  // Track Cuts
  AliFemtoDreamTrackCuts *TrackCuts = AliFemtoDreamTrackCuts::PrimProtonCuts(
      isMC, true, false, false);
  TrackCuts->SetFilterBit(128);
  TrackCuts->SetCutCharge(1);

  AliFemtoDreamTrackCuts *AntiTrackCuts =
      AliFemtoDreamTrackCuts::PrimProtonCuts(isMC, true, false, false);
  AntiTrackCuts->SetFilterBit(128);
  AntiTrackCuts->SetCutCharge(-1);

  if (suffix != "0" && suffix != "999") {
    TrackCuts->SetMinimalBooking(true);
    AntiTrackCuts->SetMinimalBooking(true);
  }
  if (suffix == "1") {
    TrackCuts->SetPtRange(ProtonPtlow, 4.05);
    AntiTrackCuts->SetPtRange(ProtonPtlow, 4.05);
    TrackCuts->SetEtaRange(-ProtonEtaUp, ProtonEtaUp);
    AntiTrackCuts->SetEtaRange(-ProtonEtaUp, ProtonEtaUp);
  } else if (suffix == "2") {
    TrackCuts->SetPtRange(ProtonPtup, 4.05);
    AntiTrackCuts->SetPtRange(ProtonPtup, 4.05);
    TrackCuts->SetPID(AliPID::kProton, 0.75, ProtonNsigmaLow);
    AntiTrackCuts->SetPID(AliPID::kProton, 0.75, ProtonNsigmaLow);
  } else if (suffix == "3") {
    TrackCuts->SetEtaRange(-ProtonEtaUp, ProtonEtaUp);
    AntiTrackCuts->SetEtaRange(-ProtonEtaUp, ProtonEtaUp);
    TrackCuts->SetNClsTPC(ProtonNClsUp);
    AntiTrackCuts->SetNClsTPC(ProtonNClsUp);
  } else if (suffix == "4") {
    TrackCuts->SetPtRange(ProtonPtup, 4.05);
    AntiTrackCuts->SetPtRange(ProtonPtup, 4.05);
    TrackCuts->SetPID(AliPID::kProton, 0.75, ProtonNsigmaUp);
    AntiTrackCuts->SetPID(AliPID::kProton, 0.75, ProtonNsigmaUp);
  } else if (suffix == "5") {
    TrackCuts->SetEtaRange(-ProtonEtaLow, ProtonEtaLow);
    AntiTrackCuts->SetEtaRange(-ProtonEtaLow, ProtonEtaLow);
    TrackCuts->SetNClsTPC(ProtonNClsLow);
    AntiTrackCuts->SetNClsTPC(ProtonNClsLow);
  } else if (suffix == "6") {
    TrackCuts->SetPtRange(ProtonPtup, 4.05);
    AntiTrackCuts->SetPtRange(ProtonPtup, 4.05);
    TrackCuts->SetPID(AliPID::kProton, 0.75, ProtonNsigmaUp);
    AntiTrackCuts->SetPID(AliPID::kProton, 0.75, ProtonNsigmaUp);
  } else if (suffix == "7") {
    TrackCuts->SetEtaRange(-ProtonEtaUp, ProtonEtaUp);
    AntiTrackCuts->SetEtaRange(-ProtonEtaUp, ProtonEtaUp);
    TrackCuts->SetNClsTPC(ProtonNClsLow);
    AntiTrackCuts->SetNClsTPC(ProtonNClsLow);
  } else if (suffix == "8") {
    TrackCuts->SetPtRange(ProtonPtup, 4.05);
    AntiTrackCuts->SetPtRange(ProtonPtup, 4.05);
    TrackCuts->SetPID(AliPID::kProton, 0.75, ProtonNsigmaLow);
    AntiTrackCuts->SetPID(AliPID::kProton, 0.75, ProtonNsigmaLow);
  } else if (suffix == "9") {
    TrackCuts->SetPtRange(ProtonPtup, 4.05);
    AntiTrackCuts->SetPtRange(ProtonPtup, 4.05);
    TrackCuts->SetNClsTPC(ProtonNClsUp);
    AntiTrackCuts->SetNClsTPC(ProtonNClsUp);
  } else if (suffix == "10") {
    TrackCuts->SetPID(AliPID::kProton, 0.75, ProtonNsigmaLow);
    AntiTrackCuts->SetPID(AliPID::kProton, 0.75, ProtonNsigmaLow);
    TrackCuts->SetNClsTPC(ProtonNClsLow);
    AntiTrackCuts->SetNClsTPC(ProtonNClsLow);
  } else if (suffix == "11") {
    TrackCuts->SetPtRange(ProtonPtup, 4.05);
    AntiTrackCuts->SetPtRange(ProtonPtup, 4.05);
    TrackCuts->SetPID(AliPID::kProton, 0.75, ProtonNsigmaLow);
    AntiTrackCuts->SetPID(AliPID::kProton, 0.75, ProtonNsigmaLow);
  } else if (suffix == "12") {
    TrackCuts->SetPtRange(ProtonPtup, 4.05);
    AntiTrackCuts->SetPtRange(ProtonPtup, 4.05);
    TrackCuts->SetNClsTPC(ProtonNClsUp);
    AntiTrackCuts->SetNClsTPC(ProtonNClsUp);
  } else if (suffix == "13") {
    TrackCuts->SetEtaRange(-ProtonEtaUp, ProtonEtaUp);
    AntiTrackCuts->SetEtaRange(-ProtonEtaUp, ProtonEtaUp);
    TrackCuts->SetPtRange(ProtonPtlow, 4.05);
    AntiTrackCuts->SetPtRange(ProtonPtlow, 4.05);
  } else if (suffix == "14") {
    TrackCuts->SetPID(AliPID::kProton, 0.75, ProtonNsigmaLow);
    AntiTrackCuts->SetPID(AliPID::kProton, 0.75, ProtonNsigmaLow);
    TrackCuts->SetEtaRange(-ProtonEtaLow, ProtonEtaLow);
    AntiTrackCuts->SetEtaRange(-ProtonEtaLow, ProtonEtaLow);
  } else if (suffix == "15") {
    TrackCuts->SetNClsTPC(ProtonNClsUp);
    AntiTrackCuts->SetNClsTPC(ProtonNClsUp);
    TrackCuts->SetPtRange(ProtonPtup, 4.05);
    AntiTrackCuts->SetPtRange(ProtonPtup, 4.05);
  } else if (suffix == "16") {
    TrackCuts->SetEtaRange(-ProtonEtaLow, ProtonEtaLow);
    AntiTrackCuts->SetEtaRange(-ProtonEtaLow, ProtonEtaLow);
    TrackCuts->SetPID(AliPID::kProton, 0.75, ProtonNsigmaLow);
    AntiTrackCuts->SetPID(AliPID::kProton, 0.75, ProtonNsigmaLow);
  } else if (suffix == "17") {
    TrackCuts->SetPtRange(ProtonPtup, 4.05);
    AntiTrackCuts->SetPtRange(ProtonPtup, 4.05);
    TrackCuts->SetNClsTPC(ProtonNClsLow);
    AntiTrackCuts->SetNClsTPC(ProtonNClsLow);
  } else if (suffix == "18") {
    TrackCuts->SetPtRange(ProtonPtup, 4.05);
    AntiTrackCuts->SetPtRange(ProtonPtup, 4.05);
    TrackCuts->SetPID(AliPID::kProton, 0.75, ProtonNsigmaUp);
    AntiTrackCuts->SetPID(AliPID::kProton, 0.75, ProtonNsigmaUp);
  } else if (suffix == "19") {
    TrackCuts->SetPtRange(ProtonPtup, 4.05);
    AntiTrackCuts->SetPtRange(ProtonPtup, 4.05);
    TrackCuts->SetPID(AliPID::kProton, 0.75, ProtonNsigmaLow);
    AntiTrackCuts->SetPID(AliPID::kProton, 0.75, ProtonNsigmaLow);
  } else if (suffix == "20") {
    TrackCuts->SetPtRange(ProtonPtup, 4.05);
    AntiTrackCuts->SetPtRange(ProtonPtup, 4.05);
    TrackCuts->SetPID(AliPID::kProton, 0.75, ProtonNsigmaLow);
    AntiTrackCuts->SetPID(AliPID::kProton, 0.75, ProtonNsigmaLow);
  } else if (suffix == "21") {
    TrackCuts->SetPtRange(ProtonPtup, 4.05);
    AntiTrackCuts->SetPtRange(ProtonPtup, 4.05);
    TrackCuts->SetNClsTPC(ProtonNClsLow);
    AntiTrackCuts->SetNClsTPC(ProtonNClsLow);
  } else if (suffix == "22") {
    TrackCuts->SetPtRange(ProtonPtup, 4.05);
    AntiTrackCuts->SetPtRange(ProtonPtup, 4.05);
    TrackCuts->SetPID(AliPID::kProton, 0.75, ProtonNsigmaLow);
    AntiTrackCuts->SetPID(AliPID::kProton, 0.75, ProtonNsigmaLow);
  } else if (suffix == "23") {
    TrackCuts->SetPtRange(ProtonPtlow, 4.05);
    AntiTrackCuts->SetPtRange(ProtonPtlow, 4.05);
    TrackCuts->SetNClsTPC(ProtonNClsUp);
    AntiTrackCuts->SetNClsTPC(ProtonNClsUp);
  } else if (suffix == "24") {
    TrackCuts->SetEtaRange(-ProtonEtaLow, ProtonEtaLow);
    AntiTrackCuts->SetEtaRange(-ProtonEtaLow, ProtonEtaLow);
    TrackCuts->SetPID(AliPID::kProton, 0.75, ProtonNsigmaUp);
    AntiTrackCuts->SetPID(AliPID::kProton, 0.75, ProtonNsigmaUp);
  } else if (suffix == "25") {
    TrackCuts->SetEtaRange(-ProtonEtaUp, ProtonEtaUp);
    AntiTrackCuts->SetEtaRange(-ProtonEtaUp, ProtonEtaUp);
    TrackCuts->SetNClsTPC(ProtonNClsLow);
    AntiTrackCuts->SetNClsTPC(ProtonNClsLow);
  }

  // =====================================================================
  // Lambda Cuts
  const float LambdaPtLow = 0.24;
  const float LambdaPtUp = 0.36;
  const float LambdaCPALow = 0.995;
  const float LambdaCPAUp = 0.99925;
  const float LambdaEtaLow = 0.85;
  const float LambdaNsigmaLow = 4;
  const float LambdaNsigmaUp = 6;
  const float LambdaNClsLow = 60;
  const float LambdaNClsUp = 80;
  const float LambdaMassUp = 0.008;

  AliFemtoDreamv0Cuts *v0Cuts = AliFemtoDreamv0Cuts::LambdaSigma0Cuts(isMC,
                                                                      false,
                                                                      false);
  AliFemtoDreamTrackCuts *Posv0Daug = AliFemtoDreamTrackCuts::DecayProtonCuts(
      isMC, false, false);
  Posv0Daug->SetEtaRange(-0.9, 0.9);
  AliFemtoDreamTrackCuts *Negv0Daug = AliFemtoDreamTrackCuts::DecayPionCuts(
      isMC, false, false);
  Negv0Daug->SetEtaRange(-0.9, 0.9);

  AliFemtoDreamv0Cuts *antiv0Cuts = AliFemtoDreamv0Cuts::LambdaSigma0Cuts(
      isMC, false, false);
  AliFemtoDreamTrackCuts *PosAntiv0Daug = AliFemtoDreamTrackCuts::DecayPionCuts(
      isMC, false, false);
  PosAntiv0Daug->SetCutCharge(1);
  PosAntiv0Daug->SetEtaRange(-0.9, 0.9);
  AliFemtoDreamTrackCuts *NegAntiv0Daug =
      AliFemtoDreamTrackCuts::DecayProtonCuts(isMC, false, false);
  NegAntiv0Daug->SetCutCharge(-1);
  NegAntiv0Daug->SetEtaRange(-0.9, 0.9);

  if (suffix != "0") {
    v0Cuts->SetMinimalBooking(true);
    antiv0Cuts->SetMinimalBooking(true);
  }
  if (suffix == "1") {
    v0Cuts->SetCutCPA(LambdaCPALow);
    antiv0Cuts->SetCutCPA(LambdaCPALow);
    Posv0Daug->SetEtaRange(-LambdaEtaLow, LambdaEtaLow);
    Negv0Daug->SetEtaRange(-LambdaEtaLow, LambdaEtaLow);
    PosAntiv0Daug->SetEtaRange(-LambdaEtaLow, LambdaEtaLow);
    NegAntiv0Daug->SetEtaRange(-LambdaEtaLow, LambdaEtaLow);
  } else if (suffix == "2") {
    v0Cuts->SetPtRange(LambdaPtLow, 999.f);
    antiv0Cuts->SetPtRange(LambdaPtLow, 999.f);
    Posv0Daug->SetEtaRange(-LambdaEtaLow, LambdaEtaLow);
    Negv0Daug->SetEtaRange(-LambdaEtaLow, LambdaEtaLow);
    PosAntiv0Daug->SetEtaRange(-LambdaEtaLow, LambdaEtaLow);
    NegAntiv0Daug->SetEtaRange(-LambdaEtaLow, LambdaEtaLow);
    Posv0Daug->SetPID(AliPID::kProton, 999.9, LambdaNsigmaUp);
    Negv0Daug->SetPID(AliPID::kPion, 999.9, LambdaNsigmaUp);
    PosAntiv0Daug->SetPID(AliPID::kPion, 999.9, LambdaNsigmaUp);
    NegAntiv0Daug->SetPID(AliPID::kProton, 999.9, LambdaNsigmaUp);
    Posv0Daug->SetNClsTPC(LambdaNClsLow);
    Negv0Daug->SetNClsTPC(LambdaNClsLow);
    PosAntiv0Daug->SetNClsTPC(LambdaNClsLow);
    NegAntiv0Daug->SetNClsTPC(LambdaNClsLow);
  } else if (suffix == "3") {
    v0Cuts->SetPtRange(LambdaPtUp, 999.f);
    antiv0Cuts->SetPtRange(LambdaPtUp, 999.f);
    Posv0Daug->SetPID(AliPID::kProton, 999.9, LambdaNsigmaUp);
    Negv0Daug->SetPID(AliPID::kPion, 999.9, LambdaNsigmaUp);
    PosAntiv0Daug->SetPID(AliPID::kPion, 999.9, LambdaNsigmaUp);
    NegAntiv0Daug->SetPID(AliPID::kProton, 999.9, LambdaNsigmaUp);
  } else if (suffix == "4") {
    Posv0Daug->SetPID(AliPID::kProton, 999.9, LambdaNsigmaUp);
    Negv0Daug->SetPID(AliPID::kPion, 999.9, LambdaNsigmaUp);
    PosAntiv0Daug->SetPID(AliPID::kPion, 999.9, LambdaNsigmaUp);
    NegAntiv0Daug->SetPID(AliPID::kProton, 999.9, LambdaNsigmaUp);
    Posv0Daug->SetEtaRange(-LambdaEtaLow, LambdaEtaLow);
    Negv0Daug->SetEtaRange(-LambdaEtaLow, LambdaEtaLow);
    PosAntiv0Daug->SetEtaRange(-LambdaEtaLow, LambdaEtaLow);
    NegAntiv0Daug->SetEtaRange(-LambdaEtaLow, LambdaEtaLow);
    v0Cuts->SetCutInvMass(LambdaMassUp);
    antiv0Cuts->SetCutInvMass(LambdaMassUp);
  } else if (suffix == "5") {
    v0Cuts->SetCutCPA(LambdaCPALow);
    antiv0Cuts->SetCutCPA(LambdaCPALow);
  } else if (suffix == "6") {
    v0Cuts->SetPtRange(LambdaPtUp, 999.f);
    antiv0Cuts->SetPtRange(LambdaPtUp, 999.f);
    Posv0Daug->SetPID(AliPID::kProton, 999.9, LambdaNsigmaLow);
    Negv0Daug->SetPID(AliPID::kPion, 999.9, LambdaNsigmaLow);
    PosAntiv0Daug->SetPID(AliPID::kPion, 999.9, LambdaNsigmaLow);
    NegAntiv0Daug->SetPID(AliPID::kProton, 999.9, LambdaNsigmaLow);
    Posv0Daug->SetEtaRange(-LambdaEtaLow, LambdaEtaLow);
    Negv0Daug->SetEtaRange(-LambdaEtaLow, LambdaEtaLow);
    PosAntiv0Daug->SetEtaRange(-LambdaEtaLow, LambdaEtaLow);
    NegAntiv0Daug->SetEtaRange(-LambdaEtaLow, LambdaEtaLow);
  } else if (suffix == "7") {
    v0Cuts->SetPtRange(LambdaPtLow, 999.f);
    antiv0Cuts->SetPtRange(LambdaPtLow, 999.f);
    Posv0Daug->SetNClsTPC(LambdaNClsUp);
    Negv0Daug->SetNClsTPC(LambdaNClsUp);
    PosAntiv0Daug->SetNClsTPC(LambdaNClsUp);
    NegAntiv0Daug->SetNClsTPC(LambdaNClsUp);
  } else if (suffix == "8") {
    v0Cuts->SetCutCPA(LambdaCPAUp);
    antiv0Cuts->SetCutCPA(LambdaCPAUp);
  } else if (suffix == "9") {
    v0Cuts->SetCutCPA(LambdaCPAUp);
    antiv0Cuts->SetCutCPA(LambdaCPAUp);
    Posv0Daug->SetPID(AliPID::kProton, 999.9, LambdaNsigmaUp);
    Negv0Daug->SetPID(AliPID::kPion, 999.9, LambdaNsigmaUp);
    PosAntiv0Daug->SetPID(AliPID::kPion, 999.9, LambdaNsigmaUp);
    NegAntiv0Daug->SetPID(AliPID::kProton, 999.9, LambdaNsigmaUp);
    Posv0Daug->SetNClsTPC(LambdaNClsLow);
    Negv0Daug->SetNClsTPC(LambdaNClsLow);
    PosAntiv0Daug->SetNClsTPC(LambdaNClsLow);
    NegAntiv0Daug->SetNClsTPC(LambdaNClsLow);
  } else if (suffix == "10") {
    v0Cuts->SetPtRange(LambdaPtUp, 999.f);
    antiv0Cuts->SetPtRange(LambdaPtUp, 999.f);
    Posv0Daug->SetNClsTPC(LambdaNClsLow);
    Negv0Daug->SetNClsTPC(LambdaNClsLow);
    PosAntiv0Daug->SetNClsTPC(LambdaNClsLow);
    NegAntiv0Daug->SetNClsTPC(LambdaNClsLow);
    Posv0Daug->SetEtaRange(-LambdaEtaLow, LambdaEtaLow);
    Negv0Daug->SetEtaRange(-LambdaEtaLow, LambdaEtaLow);
    PosAntiv0Daug->SetEtaRange(-LambdaEtaLow, LambdaEtaLow);
    NegAntiv0Daug->SetEtaRange(-LambdaEtaLow, LambdaEtaLow);
  } else if (suffix == "11") {
    v0Cuts->SetPtRange(LambdaPtUp, 999.f);
    antiv0Cuts->SetPtRange(LambdaPtUp, 999.f);
    v0Cuts->SetCutCPA(LambdaCPALow);
    antiv0Cuts->SetCutCPA(LambdaCPALow);
  } else if (suffix == "12") {
    v0Cuts->SetCutCPA(LambdaCPAUp);
    antiv0Cuts->SetCutCPA(LambdaCPAUp);
    Posv0Daug->SetPID(AliPID::kProton, 999.9, LambdaNsigmaLow);
    Negv0Daug->SetPID(AliPID::kPion, 999.9, LambdaNsigmaLow);
    PosAntiv0Daug->SetPID(AliPID::kPion, 999.9, LambdaNsigmaLow);
    NegAntiv0Daug->SetPID(AliPID::kProton, 999.9, LambdaNsigmaLow);
    Posv0Daug->SetNClsTPC(LambdaNClsUp);
    Negv0Daug->SetNClsTPC(LambdaNClsUp);
    PosAntiv0Daug->SetNClsTPC(LambdaNClsUp);
    NegAntiv0Daug->SetNClsTPC(LambdaNClsUp);
  } else if (suffix == "13") {
    v0Cuts->SetCutCPA(LambdaCPAUp);
    antiv0Cuts->SetCutCPA(LambdaCPAUp);
    Posv0Daug->SetNClsTPC(LambdaNClsLow);
    Negv0Daug->SetNClsTPC(LambdaNClsLow);
    PosAntiv0Daug->SetNClsTPC(LambdaNClsLow);
    NegAntiv0Daug->SetNClsTPC(LambdaNClsLow);
  } else if (suffix == "14") {
    Posv0Daug->SetPID(AliPID::kProton, 999.9, LambdaNsigmaUp);
    Negv0Daug->SetPID(AliPID::kPion, 999.9, LambdaNsigmaUp);
    PosAntiv0Daug->SetPID(AliPID::kPion, 999.9, LambdaNsigmaUp);
    NegAntiv0Daug->SetPID(AliPID::kProton, 999.9, LambdaNsigmaUp);
    Posv0Daug->SetEtaRange(-LambdaEtaLow, LambdaEtaLow);
    Negv0Daug->SetEtaRange(-LambdaEtaLow, LambdaEtaLow);
    PosAntiv0Daug->SetEtaRange(-LambdaEtaLow, LambdaEtaLow);
    NegAntiv0Daug->SetEtaRange(-LambdaEtaLow, LambdaEtaLow);
    Posv0Daug->SetNClsTPC(LambdaNClsLow);
    Negv0Daug->SetNClsTPC(LambdaNClsLow);
    PosAntiv0Daug->SetNClsTPC(LambdaNClsLow);
    NegAntiv0Daug->SetNClsTPC(LambdaNClsLow);
  } else if (suffix == "15") {
    v0Cuts->SetCutCPA(LambdaCPAUp);
    antiv0Cuts->SetCutCPA(LambdaCPAUp);
    v0Cuts->SetCutInvMass(LambdaMassUp);
    antiv0Cuts->SetCutInvMass(LambdaMassUp);
    Posv0Daug->SetNClsTPC(LambdaNClsLow);
    Negv0Daug->SetNClsTPC(LambdaNClsLow);
    PosAntiv0Daug->SetNClsTPC(LambdaNClsLow);
    NegAntiv0Daug->SetNClsTPC(LambdaNClsLow);
  } else if (suffix == "16") {
    v0Cuts->SetPtRange(LambdaPtUp, 999.f);
    antiv0Cuts->SetPtRange(LambdaPtUp, 999.f);
    Posv0Daug->SetNClsTPC(LambdaNClsLow);
    Negv0Daug->SetNClsTPC(LambdaNClsLow);
    PosAntiv0Daug->SetNClsTPC(LambdaNClsLow);
    NegAntiv0Daug->SetNClsTPC(LambdaNClsLow);
  } else if (suffix == "17") {
    v0Cuts->SetPtRange(LambdaPtUp, 999.f);
    antiv0Cuts->SetPtRange(LambdaPtUp, 999.f);
    v0Cuts->SetCutCPA(LambdaCPALow);
    antiv0Cuts->SetCutCPA(LambdaCPALow);
    Posv0Daug->SetNClsTPC(LambdaNClsLow);
    Negv0Daug->SetNClsTPC(LambdaNClsLow);
    PosAntiv0Daug->SetNClsTPC(LambdaNClsLow);
    NegAntiv0Daug->SetNClsTPC(LambdaNClsLow);
  } else if (suffix == "18") {
    v0Cuts->SetCutCPA(LambdaCPAUp);
    antiv0Cuts->SetCutCPA(LambdaCPAUp);
    Posv0Daug->SetNClsTPC(LambdaNClsLow);
    Negv0Daug->SetNClsTPC(LambdaNClsLow);
    PosAntiv0Daug->SetNClsTPC(LambdaNClsLow);
    NegAntiv0Daug->SetNClsTPC(LambdaNClsLow);
  } else if (suffix == "19") {
    v0Cuts->SetCutCPA(LambdaCPALow);
    antiv0Cuts->SetCutCPA(LambdaCPALow);
    Posv0Daug->SetPID(AliPID::kProton, 999.9, LambdaNsigmaLow);
    Negv0Daug->SetPID(AliPID::kPion, 999.9, LambdaNsigmaLow);
    PosAntiv0Daug->SetPID(AliPID::kPion, 999.9, LambdaNsigmaLow);
    NegAntiv0Daug->SetPID(AliPID::kProton, 999.9, LambdaNsigmaLow);
  } else if (suffix == "20") {
    v0Cuts->SetCutCPA(LambdaCPAUp);
    antiv0Cuts->SetCutCPA(LambdaCPAUp);
    v0Cuts->SetCutInvMass(LambdaMassUp);
    antiv0Cuts->SetCutInvMass(LambdaMassUp);
    v0Cuts->SetPtRange(LambdaPtUp, 999.f);
    antiv0Cuts->SetPtRange(LambdaPtUp, 999.f);
  } else if (suffix == "21") {
    v0Cuts->SetCutCPA(LambdaCPAUp);
    antiv0Cuts->SetCutCPA(LambdaCPAUp);
    Posv0Daug->SetNClsTPC(LambdaNClsLow);
    Negv0Daug->SetNClsTPC(LambdaNClsLow);
    PosAntiv0Daug->SetNClsTPC(LambdaNClsLow);
    NegAntiv0Daug->SetNClsTPC(LambdaNClsLow);
    Posv0Daug->SetEtaRange(-LambdaEtaLow, LambdaEtaLow);
    Negv0Daug->SetEtaRange(-LambdaEtaLow, LambdaEtaLow);
    PosAntiv0Daug->SetEtaRange(-LambdaEtaLow, LambdaEtaLow);
    NegAntiv0Daug->SetEtaRange(-LambdaEtaLow, LambdaEtaLow);
  } else if (suffix == "22") {
    v0Cuts->SetPtRange(LambdaPtUp, 999.f);
    antiv0Cuts->SetPtRange(LambdaPtUp, 999.f);
    Posv0Daug->SetPID(AliPID::kProton, 999.9, LambdaNsigmaLow);
    Negv0Daug->SetPID(AliPID::kPion, 999.9, LambdaNsigmaLow);
    PosAntiv0Daug->SetPID(AliPID::kPion, 999.9, LambdaNsigmaLow);
    NegAntiv0Daug->SetPID(AliPID::kProton, 999.9, LambdaNsigmaLow);
  } else if (suffix == "23") {
    v0Cuts->SetPtRange(LambdaPtLow, 999.f);
    antiv0Cuts->SetPtRange(LambdaPtLow, 999.f);
    Posv0Daug->SetEtaRange(-LambdaEtaLow, LambdaEtaLow);
    Negv0Daug->SetEtaRange(-LambdaEtaLow, LambdaEtaLow);
    PosAntiv0Daug->SetEtaRange(-LambdaEtaLow, LambdaEtaLow);
    NegAntiv0Daug->SetEtaRange(-LambdaEtaLow, LambdaEtaLow);
  } else if (suffix == "24") {
    v0Cuts->SetPtRange(LambdaPtUp, 999.f);
    antiv0Cuts->SetPtRange(LambdaPtUp, 999.f);
    Posv0Daug->SetNClsTPC(LambdaNClsLow);
    Negv0Daug->SetNClsTPC(LambdaNClsLow);
    PosAntiv0Daug->SetNClsTPC(LambdaNClsLow);
    NegAntiv0Daug->SetNClsTPC(LambdaNClsLow);
  } else if (suffix == "25") {
    v0Cuts->SetPtRange(LambdaPtLow, 999.f);
    antiv0Cuts->SetPtRange(LambdaPtLow, 999.f);
    v0Cuts->SetCutCPA(LambdaCPAUp);
    antiv0Cuts->SetCutCPA(LambdaCPAUp);
    Posv0Daug->SetPID(AliPID::kProton, 999.9, LambdaNsigmaUp);
    Negv0Daug->SetPID(AliPID::kPion, 999.9, LambdaNsigmaUp);
    PosAntiv0Daug->SetPID(AliPID::kPion, 999.9, LambdaNsigmaUp);
    NegAntiv0Daug->SetPID(AliPID::kProton, 999.9, LambdaNsigmaUp);
  }

  v0Cuts->SetPosDaugterTrackCuts(Posv0Daug);
  v0Cuts->SetNegDaugterTrackCuts(Negv0Daug);
  v0Cuts->SetPDGCodePosDaug(2212);  // Proton
  v0Cuts->SetPDGCodeNegDaug(211);   // Pion
  v0Cuts->SetPDGCodev0(3122);       // Lambda
  antiv0Cuts->SetPosDaugterTrackCuts(PosAntiv0Daug);
  antiv0Cuts->SetNegDaugterTrackCuts(NegAntiv0Daug);
  antiv0Cuts->SetPDGCodePosDaug(211);   // Pion
  antiv0Cuts->SetPDGCodeNegDaug(2212);  // Proton
  antiv0Cuts->SetPDGCodev0(-3122);      // Lambda

  // =====================================================================
  // Photon Cuts
  const float PhotonPtLow = 0.;
  const float PhotonPtUp = 0.15;
  const float PhotonDaughterPtLow = 0.0;
  const float PhotonDaughterPtUp = 0.075;
  const float PhotonCPALow = 0.995;
  const float PhotonCPAUp = 0.99925;
  const float PhotonDaughterEtaLow = 0.85;
  const float PhotonNsigmaDownLow = -7;
  const float PhotonNsigmaDownUp = -5;
  const float PhotonNsigmaHighLow = 6;
  const float PhotonNsigmaHighUp = 8;
  const float PhotonFClsLow = 0.3;
  const float PhotonFClsUp = 0.5;
  const float PhotonPsiPairLow = 0.15;
  const float PhotonPsiPairUp = 0.3;
  const float PhotonArmenterosLow = 0.05;
  const float PhotonArmenterosUp = 0.1;

  AliSigma0PhotonCuts *photon = AliSigma0PhotonCuts::PhotonCuts();
  photon->SetIsMC(isMC);
  if (suffix != "0") {
    photon->SetLightweight(true);
  }
  if (suffix == "1") {
    photon->SetElectronRatioFindable(PhotonFClsLow);
    photon->SetElectronEtaMax(PhotonDaughterEtaLow);
  } else if (suffix == "2") {
    photon->SetArmenterosQtMax(PhotonArmenterosLow);
    photon->SetPhotonPtMin(PhotonPtUp);
    photon->SetElectronEtaMax(PhotonDaughterEtaLow);
  } else if (suffix == "3") {
    photon->SetElectronPtMin(PhotonDaughterPtUp);
    photon->SetPsiPairMax(PhotonPsiPairUp);
  } else if (suffix == "4") {
    photon->SetElectronNSigmaTPCMax(PhotonNsigmaHighLow);
    photon->SetElectronNSigmaTPCMin(PhotonNsigmaDownLow);
  } else if (suffix == "5") {
    photon->SetElectronEtaMax(PhotonDaughterEtaLow);
    photon->SetCPAMin(PhotonCPAUp);
  } else if (suffix == "6") {
    photon->SetPhotonPtMin(PhotonPtLow);
    photon->SetPsiPairMax(0.2);
    photon->SetCPAMin(PhotonCPAUp);
  } else if (suffix == "7") {
    photon->SetArmenterosQtMax(PhotonArmenterosLow);
    photon->SetElectronNSigmaTPCMax(PhotonNsigmaHighLow);
    photon->SetElectronNSigmaTPCMin(PhotonNsigmaDownUp);
  } else if (suffix == "8") {
    photon->SetElectronRatioFindable(PhotonFClsLow);
    photon->SetCPAMin(PhotonCPALow);
    photon->SetPsiPairMax(PhotonPsiPairUp);
  } else if (suffix == "9") {
    photon->SetElectronEtaMax(PhotonDaughterEtaLow);
    photon->SetElectronRatioFindable(PhotonFClsUp);
  } else if (suffix == "10") {
    photon->SetPhotonPtMin(PhotonPtUp);
    photon->SetArmenterosQtMax(PhotonArmenterosUp);
  } else if (suffix == "11") {
    photon->SetElectronPtMin(PhotonDaughterPtLow);
    photon->SetPsiPairMax(PhotonPsiPairLow);
  } else if (suffix == "12") {
    photon->SetElectronEtaMax(PhotonDaughterEtaLow);
  } else if (suffix == "13") {
    photon->SetPhotonPtMin(PhotonPtLow);
    photon->SetArmenterosQtMax(PhotonArmenterosLow);
    photon->SetElectronNSigmaTPCMax(PhotonNsigmaHighUp);
    photon->SetElectronNSigmaTPCMin(PhotonNsigmaDownLow);
    photon->SetCPAMin(PhotonCPAUp);
  } else if (suffix == "14") {
    photon->SetPsiPairMax(PhotonPsiPairLow);
  } else if (suffix == "15") {
    photon->SetElectronEtaMax(PhotonDaughterEtaLow);
    photon->SetPsiPairMax(PhotonPsiPairUp);
  } else if (suffix == "16") {
    photon->SetPhotonPtMin(PhotonPtLow);
    photon->SetElectronRatioFindable(PhotonFClsUp);
  } else if (suffix == "17") {
    photon->SetPhotonPtMin(PhotonPtUp);
    photon->SetArmenterosQtMax(PhotonArmenterosUp);
    photon->SetCPAMin(PhotonCPALow);
  } else if (suffix == "18") {
    photon->SetElectronPtMin(PhotonDaughterPtLow);
    photon->SetPsiPairMax(PhotonPsiPairUp);
  } else if (suffix == "19") {
    photon->SetElectronNSigmaTPCMax(PhotonNsigmaHighUp);
    photon->SetElectronNSigmaTPCMin(PhotonNsigmaDownUp);
  } else if (suffix == "20") {
    photon->SetArmenterosQtMax(PhotonArmenterosUp);
    photon->SetElectronEtaMax(PhotonDaughterEtaLow);
    photon->SetPsiPairMax(PhotonPsiPairLow);
  } else if (suffix == "21") {
    photon->SetElectronRatioFindable(PhotonFClsLow);
  } else if (suffix == "22") {
    photon->SetPhotonPtMin(PhotonPtUp);
    photon->SetCPAMin(PhotonCPALow);
    photon->SetElectronEtaMax(PhotonDaughterEtaLow);
  } else if (suffix == "23") {
    photon->SetArmenterosQtMax(PhotonArmenterosLow);
    photon->SetElectronPtMin(PhotonDaughterPtLow);
    photon->SetCPAMin(PhotonCPALow);
  } else if (suffix == "24") {
    photon->SetElectronNSigmaTPCMax(PhotonNsigmaHighUp);
    photon->SetElectronNSigmaTPCMin(PhotonNsigmaDownLow);
    photon->SetPsiPairMax(PhotonPsiPairUp);
  } else if (suffix == "25") {
    photon->SetElectronPtMin(PhotonDaughterPtUp);
    photon->SetElectronRatioFindable(PhotonFClsLow);
  }

  // =====================================================================
  // Sigma Cuts
  const float sidebandDownLow = 0.0035;
  const float sidebandDownDefault = 0.005;
  const float sidebandDownUp = 0.0075;
  const float sidebandHighUp = 0.075;
  const float sidebandHighDefault = 0.05;
  const float sidebandHighLow = 0.035;

  AliSigma0AODPhotonMotherCuts *sigmaCuts =
      AliSigma0AODPhotonMotherCuts::DefaultCuts();
  sigmaCuts->SetIsMC(isMC);
  sigmaCuts->SetPDG(3212, 3122, 22);
  if (suffix != "0" && suffix != "999") {
    sigmaCuts->SetLightweight(true);
  }

  AliSigma0AODPhotonMotherCuts *antiSigmaCuts =
      AliSigma0AODPhotonMotherCuts::DefaultCuts();
  antiSigmaCuts->SetIsMC(isMC);
  antiSigmaCuts->SetPDG(-3212, -3122, 22);
  if (suffix != "0" && suffix != "999") {
    antiSigmaCuts->SetLightweight(true);
  }

  // vary the sidebands
  if (suffix == "1") {
    sigmaCuts->SetSigmaSideband(sidebandDownLow, sidebandHighLow);
    antiSigmaCuts->SetSigmaSideband(sidebandDownLow, sidebandHighLow);
  } else if (suffix == "2") {
    sigmaCuts->SetSigmaSideband(sidebandDownLow, sidebandHighDefault);
    antiSigmaCuts->SetSigmaSideband(sidebandDownLow, sidebandHighDefault);
  } else if (suffix == "3") {
    sigmaCuts->SetSigmaSideband(sidebandDownUp, sidebandHighUp);
    antiSigmaCuts->SetSigmaSideband(sidebandDownUp, sidebandHighUp);
  } else if (suffix == "4") {
    sigmaCuts->SetSigmaSideband(sidebandDownUp, sidebandHighUp);
    antiSigmaCuts->SetSigmaSideband(sidebandDownUp, sidebandHighUp);
  } else if (suffix == "5") {
    sigmaCuts->SetSigmaSideband(sidebandDownDefault, sidebandHighUp);
    antiSigmaCuts->SetSigmaSideband(sidebandDownDefault, sidebandHighUp);
  } else if (suffix == "6") {
    sigmaCuts->SetSigmaSideband(sidebandDownUp, sidebandHighDefault);
    antiSigmaCuts->SetSigmaSideband(sidebandDownUp, sidebandHighDefault);
  } else if (suffix == "7") {
    sigmaCuts->SetSigmaSideband(sidebandDownDefault, sidebandHighLow);
    antiSigmaCuts->SetSigmaSideband(sidebandDownDefault, sidebandHighLow);
  } else if (suffix == "8") {
    sigmaCuts->SetSigmaSideband(sidebandDownUp, sidebandHighDefault);
    antiSigmaCuts->SetSigmaSideband(sidebandDownUp, sidebandHighDefault);
  } else if (suffix == "9") {
    sigmaCuts->SetSigmaSideband(sidebandDownDefault, sidebandHighLow);
    antiSigmaCuts->SetSigmaSideband(sidebandDownDefault, sidebandHighLow);
  } else if (suffix == "10") {
    sigmaCuts->SetSigmaSideband(sidebandDownUp, sidebandHighUp);
    antiSigmaCuts->SetSigmaSideband(sidebandDownUp, sidebandHighUp);
  } else if (suffix == "11") {
    sigmaCuts->SetSigmaSideband(sidebandDownUp, sidebandHighDefault);
    antiSigmaCuts->SetSigmaSideband(sidebandDownUp, sidebandHighDefault);
  } else if (suffix == "12") {
    sigmaCuts->SetSigmaSideband(sidebandDownUp, sidebandHighUp);
    antiSigmaCuts->SetSigmaSideband(sidebandDownUp, sidebandHighUp);
  } else if (suffix == "13") {
    sigmaCuts->SetSigmaSideband(sidebandDownLow, sidebandHighUp);
    antiSigmaCuts->SetSigmaSideband(sidebandDownLow, sidebandHighUp);
  } else if (suffix == "14") {
    sigmaCuts->SetSigmaSideband(sidebandDownUp, sidebandHighDefault);
    antiSigmaCuts->SetSigmaSideband(sidebandDownUp, sidebandHighDefault);
  } else if (suffix == "15") {
    sigmaCuts->SetSigmaSideband(sidebandDownLow, sidebandHighDefault);
    antiSigmaCuts->SetSigmaSideband(sidebandDownLow, sidebandHighDefault);
  } else if (suffix == "16") {
    sigmaCuts->SetSigmaSideband(sidebandDownUp, sidebandHighUp);
    antiSigmaCuts->SetSigmaSideband(sidebandDownUp, sidebandHighUp);
  } else if (suffix == "17") {
    sigmaCuts->SetSigmaSideband(sidebandDownUp, sidebandHighLow);
    antiSigmaCuts->SetSigmaSideband(sidebandDownUp, sidebandHighLow);
  } else if (suffix == "18") {
    sigmaCuts->SetSigmaSideband(sidebandDownUp, sidebandHighDefault);
    antiSigmaCuts->SetSigmaSideband(sidebandDownUp, sidebandHighDefault);
  } else if (suffix == "19") {
    sigmaCuts->SetSigmaSideband(sidebandDownUp, sidebandHighUp);
    antiSigmaCuts->SetSigmaSideband(sidebandDownUp, sidebandHighUp);
  } else if (suffix == "20") {
    sigmaCuts->SetSigmaSideband(sidebandDownDefault, sidebandHighLow);
    antiSigmaCuts->SetSigmaSideband(sidebandDownDefault, sidebandHighLow);
  } else if (suffix == "21") {
    sigmaCuts->SetSigmaSideband(sidebandDownLow, sidebandHighDefault);
    antiSigmaCuts->SetSigmaSideband(sidebandDownLow, sidebandHighDefault);
  } else if (suffix == "22") {
    sigmaCuts->SetSigmaSideband(sidebandDownLow, sidebandHighDefault);
    antiSigmaCuts->SetSigmaSideband(sidebandDownLow, sidebandHighDefault);
  } else if (suffix == "23") {
    sigmaCuts->SetSigmaSideband(sidebandDownUp, sidebandHighUp);
    antiSigmaCuts->SetSigmaSideband(sidebandDownUp, sidebandHighUp);
  } else if (suffix == "24") {
    sigmaCuts->SetSigmaSideband(sidebandDownUp, sidebandHighLow);
    antiSigmaCuts->SetSigmaSideband(sidebandDownUp, sidebandHighLow);
  } else if (suffix == "25") {
    sigmaCuts->SetSigmaSideband(sidebandDownUp, sidebandHighDefault);
    antiSigmaCuts->SetSigmaSideband(sidebandDownUp, sidebandHighDefault);
  }

  // Femto Collection
  std::vector<int> PDGParticles;
  PDGParticles.push_back(2212);  // 0
  PDGParticles.push_back(2212);  // 1
  PDGParticles.push_back(3212);  // 2
  PDGParticles.push_back(3212);  // 3
  PDGParticles.push_back(3212);  // 4
  PDGParticles.push_back(3212);  // 5
  PDGParticles.push_back(3212);  // 6
  PDGParticles.push_back(3212);  // 7

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

  std::vector<int> NBins;
  std::vector<float> kMin;
  std::vector<float> kMax;
  std::vector<int> pairQA;
  std::vector<bool> closeRejection;
  const int nPairs = 36;
  for (int i = 0; i < nPairs; ++i) {
    pairQA.push_back(0);
    closeRejection.push_back(false);
    if (suffix == "0") {
      NBins.push_back(600);
      kMin.push_back(0.);
      kMax.push_back(3.);
    } else {
      NBins.push_back(200);
      kMin.push_back(0.);
      kMax.push_back(1.);
    }
  }

  if (suffix == "0") {
    NBins[0] = 3000;  // pp
    NBins[8] = 3000;  // barp barp
  } else {
    NBins[0] = 1000;  // pp
    NBins[8] = 1000;  // barp barp
  }

  closeRejection[0] = true;  // pp
  closeRejection[8] = true;  // barp barp

  pairQA[0] = 11;   // pp
  pairQA[8] = 11;   // pbarpbar
  pairQA[2] = 14;  // pSigma
  pairQA[10] = 14;  // barp bSigma

  AliFemtoDreamCollConfig *config = new AliFemtoDreamCollConfig("Femto",
                                                                "Femto");
  std::vector<int> MultBins;
  if (trigger == "kHighMultV0") {
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
  } else {
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
    MultBins.push_back(60);
    MultBins.push_back(80);
    config->SetMultBins(MultBins);
  }
  config->SetMultBinning(true);

  config->SetExtendedQAPairs(pairQA);
  config->SetdPhidEtaPlotsSmallK(false);
  config->SetZBins(ZVtxBins);
  if (MomRes && isMC) {
    config->SetMomentumResolution(true);
  }

  config->SetDeltaEtaMax(0.012);
  config->SetDeltaPhiMax(0.012);
  config->SetClosePairRejection(closeRejection);

  if (suffix == "0" && fullBlastQA) {
    config->SetPhiEtaBinnign(true);
    config->SetkTBinning(true);
    config->SetmTBinning(true);
  }

  config->SetUsePhiSpinning(false);
  config->SetControlMethod(AliFemtoDreamCollConfig::kCorrelatedPhi);
  config->SetCorrelationRange(0.1);  // to be validated
  config->SetSpinningDepth(1);      // to be validated

  config->SetPtQA(true);
  config->SetMassQA(true);
  config->SetdPhidEtaPlots(false);
  config->SetPDGCodes(PDGParticles);
  config->SetNBinsHist(NBins);
  config->SetMinKRel(kMin);
  config->SetMaxKRel(kMax);
  config->SetMixingDepth(10);
  config->SetUseEventMixing(true);
  config->SetMultiplicityEstimator(AliFemtoDreamEvent::kRef08);
  config->SetMinimalBookingME(false);
  if (suffix != "0") {
    config->SetMinimalBookingSample(true);
  }

  AliAnalysisTaskNanoAODSigma0Femto *task =
      new AliAnalysisTaskNanoAODSigma0Femto("AliAnalysisTaskNanoAODSigma0Femto",
                                            isMC);

  task->SetEventCuts(evtCuts);
  task->SetProtonCuts(TrackCuts);
  task->SetAntiProtonCuts(AntiTrackCuts);
  task->SetV0Cuts(v0Cuts);
  task->SetAntiV0Cuts(antiv0Cuts);
  task->SetPhotonCuts(photon);
  task->SetSigmaCuts(sigmaCuts);
  task->SetAntiSigmaCuts(antiSigmaCuts);
  task->SetCollectionConfig(config);
  task->SetUseDumpster(true);

  if (suffix != "0" && suffix != "999") {
    task->SetLightweight(true);
  }

  mgr->AddTask(task);

  TString addon = "";
  if (trigger == "kINT7") {
    addon += "MBSigma0";
  } else if (trigger == "kHighMultV0") {
    addon += "HMSigma0";
  }

  TString file = AliAnalysisManager::GetCommonFileName();

  mgr->ConnectInput(task, 0, cinput);

  TString QAName = Form("%sQA%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *coutputQA = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      QAName.Data(),
      TList::Class(), AliAnalysisManager::kOutputContainer,
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

  TString v0CutsName = Form("%sv0Cuts%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *coutputv0Cuts = mgr->CreateContainer(
      v0CutsName.Data(), TList::Class(), AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), v0CutsName.Data()));
  mgr->ConnectOutput(task, 5, coutputv0Cuts);

  TString Antiv0CutsName = Form("%sAntiv0Cuts%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *coutputAntiv0Cuts = mgr->CreateContainer(
      Antiv0CutsName.Data(), TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), Antiv0CutsName.Data()));
  mgr->ConnectOutput(task, 6, coutputAntiv0Cuts);

  TString Photonv0CutsName = Form("%sPhotonCuts%s", addon.Data(),
                                  suffix.Data());
  AliAnalysisDataContainer *coutputPhotonv0Cuts = mgr->CreateContainer(
      Photonv0CutsName.Data(), TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), Photonv0CutsName.Data()));
  mgr->ConnectOutput(task, 7, coutputPhotonv0Cuts);

  TString SigmaCutsName = Form("%sSigma0Cuts%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *coutputSigmaCuts = mgr->CreateContainer(
      SigmaCutsName.Data(), TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), SigmaCutsName.Data()));
  mgr->ConnectOutput(task, 8, coutputSigmaCuts);

  TString AntiSigmaCutsName = Form("%sAntiSigma0Cuts%s", addon.Data(),
                                   suffix.Data());
  AliAnalysisDataContainer *coutputAntiSigmaCuts = mgr->CreateContainer(
      AntiSigmaCutsName.Data(), TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), AntiSigmaCutsName.Data()));
  mgr->ConnectOutput(task, 9, coutputAntiSigmaCuts);

  AliAnalysisDataContainer *coutputResults;
  TString ResultsName = Form("%sResults%s", addon.Data(), suffix.Data());
  coutputResults = mgr->CreateContainer(
      ResultsName.Data(), TList::Class(), AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), ResultsName.Data()));
  mgr->ConnectOutput(task, 10, coutputResults);

  AliAnalysisDataContainer *coutputResultQA;
  TString ResultQAName = Form("%sResultQA%s", addon.Data(), suffix.Data());
  coutputResultQA = mgr->CreateContainer(
      ResultQAName.Data(), TList::Class(), AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), ResultQAName.Data()));
  mgr->ConnectOutput(task, 11, coutputResultQA);

  AliAnalysisDataContainer *coutputDumpsterQA;
  TString DumpsterName = Form("%sDumpster%s", addon.Data(), suffix.Data());
  coutputDumpsterQA = mgr->CreateContainer(
      DumpsterName.Data(), TList::Class(), AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), DumpsterName.Data()));
  mgr->ConnectOutput(task, 12, coutputDumpsterQA);

  if (isMC) {
    TString TrkCutsMCName = Form("%sTrackCutsMC%s", addon.Data(),
                                 suffix.Data());
    AliAnalysisDataContainer *coutputTrkCutsMC = mgr->CreateContainer(
        TrkCutsMCName.Data(), TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), TrkCutsMCName.Data()));
    mgr->ConnectOutput(task, 13, coutputTrkCutsMC);

    TString AntiTrkCutsMCName = Form("%sAntiTrackCutsMC%s", addon.Data(),
                                     suffix.Data());
    AliAnalysisDataContainer *coutputAntiTrkCutsMC = mgr->CreateContainer(
        AntiTrkCutsMCName.Data(), TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), AntiTrkCutsMCName.Data()));
    mgr->ConnectOutput(task, 14, coutputAntiTrkCutsMC);

    TString V0CutsMCName = Form("%sv0CutsMC%s", addon.Data(), suffix.Data());
    AliAnalysisDataContainer *coutputV0CutsMC = mgr->CreateContainer(
        V0CutsMCName.Data(), TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), V0CutsMCName.Data()));
    mgr->ConnectOutput(task, 15, coutputV0CutsMC);

    TString AntiV0CutsMCName = Form("%sAntiv0CutsMC%s", addon.Data(),
                                    suffix.Data());
    AliAnalysisDataContainer *coutputAntiV0CutsMC = mgr->CreateContainer(
        AntiV0CutsMCName.Data(), TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), AntiV0CutsMCName.Data()));
    mgr->ConnectOutput(task, 16, coutputAntiV0CutsMC);
  }

  return task;
}
