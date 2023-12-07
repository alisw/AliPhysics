/* Copyright(c) 1998-1999, ALICfalseriment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* Add a description of your MPI analysis */

#ifndef AliAnalysisTaskDataSpeedOfSoundSim_H
#define AliAnalysisTaskDataSpeedOfSoundSim_H

class AliESDtrackCuts;
class AliESDEvent;
class TList;
class TH1F;
class TH2F;
class TH3F;
class TH3D;
class TProfile;

#include "AliAnalysisTaskSE.h"
#include "AliEventCuts.h"
#include "AliGenEventHeader.h"
#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliMCParticle.h"
#include "AliMultSelection.h"
#include "AliStack.h"
#include "AliVEvent.h"
#include "TParticle.h"

class AliAnalysisTaskDataSpeedOfSoundSim : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskDataSpeedOfSoundSim();
  AliAnalysisTaskDataSpeedOfSoundSim(const char* name);
  virtual ~AliAnalysisTaskDataSpeedOfSoundSim();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t* option);
  virtual void Terminate(Option_t* option);
  void DCAxyDistributions();
  void TrackingEfficiency();
  void DetectorResponse();
  void ReadMCEvent();
  void TrueMultiplicityDistributions();
  void MultiplicityDistributions();
  void GetCalibratedV0Amplitude();
  void SetV0Mmin(double V0Mmin) { fV0Mmin = V0Mmin; }  // Set V0M min value
  void SetV0Mmax(double V0Mmax) { fV0Mmax = V0Mmax; }  // Set V0M max value
  void SetHMCut(double HMcut) { fHMCut = HMcut; }      // Set V0M max value
  void SetRandomNumberCut(double rdcut) { fRandomNumberCut = rdcut; }
  void SetUseMC(bool mc = false) { fUseMC = mc; }  // use to analyse MC data
  void SetEtaCut(const double& etacut) { fEtaCut = etacut; }
  void SetEtaMinCut(const double& etamin) { fEtaMin = etamin; }
  void SetEtaMaxCut(const double& etamax) { fEtaMax = etamax; }
  void SetPtMin(const double& ptmin) { fPtMin = ptmin; }
  void SetTrackCuts(bool TPConly = true) { fIsTPConly = TPConly; }
  void SetTrigger(UInt_t offlineTriggerMask = AliVEvent::kINT7) {
    fTrigger = offlineTriggerMask;
  }
  bool HasRecVertex();
  void GetSPDMultiplicity();

 protected:
 private:
  AliESDEvent* fESD;
  AliEventCuts fEventCuts;
  AliStack* fMCStack;
  AliMCEvent* fMC;
  bool fUseMC;
  bool fIsTPConly;
  UInt_t fTrigger;
  AliAnalysisFilter* fTrackFilter;
  AliAnalysisFilter* fTrackFilterwoDCA;
  TList* fOutputList;
  double fEtaCut;
  double fEtaMin;
  double fEtaMax;
  double fPtMin;
  double fV0Mmin;
  double fV0Mmax;
  double fHMCut;
  double fRandomNumberCut;
  double fv0mpercentile;
  float fv0mamplitude;
  int fRecNch;
  int fTrueNch;
  int fTrueNch14;
  int fTrueNch10;
  int fTrueNchEtaPos;
  int fTrueNchEtaNeg;
  int fTrueV0;
  int fTracklets14;
  int fTracklets10;
  float fdcaxy;
  float fdcaz;
  AliMultSelection* fMultSelection;
  TH2D* hNchvsV0M;
  TH2D* hNchvsV0MAmp;
  TH2D* hV0MvsV0MAmp;
  TProfile* pV0MAmpChannel;
  TH1D* hV0MAmplitude;
  TH1F* hCounter;
  TH1F* hV0Mmult;
  TH1F* hSPDmult14;
  TH1F* hSPDmult10;
  TProfile* pPtvsNch;
  TProfile* pPtEtaNegvsNchEtaPos;
  TProfile* pPtEtaPosvsNchEtaNeg;
  TProfile* pPtvsV0MAmp;
  TH2D* hPtvsV0MAmp;
  TH2D* hNchEtaPosvsNchEtaNeg;
  TH2D* hPtEtaNegvsNchEtaPos;
  TH2D* hPtEtaPosvsNchEtaNeg;
  TH1F* hTrueVtxZ;
  TH1F* hTrueNch;
  TH1F* hTrueV0MAmp;
  TH2D* hNchResponse;
  TH2D* hTruePtvsTrueNch;
  TH2D* hTrueNchEtaPosvsTrueNchEtaNeg;
  TH2D* hTruePtEtaNegvsTrueNchEtaPos;
  TH2D* hTruePtEtaPosvsTrueNchEtaNeg;
  TH1F* hTrueNch14;
  TH1F* hTrueNch10;
  TH2D* hTruePtvsTrueNch14;
  TH2D* hTruePtvsTrueNch10;
  TH2F* hDCAxyPri[1];
  TH2F* hDCAxyWeDe[1];
  TH2F* hDCAxyMaIn[1];
  TH1F* hPtInPrim_ch;
  TH1F* hPtInPrim_pion;
  TH1F* hPtInPrim_kaon;
  TH1F* hPtInPrim_proton;
  TH1F* hPtInPrim_sigmap;
  TH1F* hPtInPrim_sigmam;
  TH1F* hPtInPrim_omega;
  TH1F* hPtInPrim_xi;
  TH1F* hPtInPrim_rest;
  TH1F* hPtOutAll_ch;
  TH1F* hPtOutPrim_ch;
  TH1F* hPtOutPrim_pion;
  TH1F* hPtOutPrim_kaon;
  TH1F* hPtOutPrim_proton;
  TH1F* hPtOutPrim_sigmap;
  TH1F* hPtOutPrim_sigmam;
  TH1F* hPtOutPrim_omega;
  TH1F* hPtOutPrim_xi;
  TH1F* hPtOutPrim_rest;
  TH2D* hPhiEtaSPD;
  TH2D* hVtxZvsTracklets;
  TH2D* hTrackletsvsV0MAmp;
  TH2D* hPtvsTracklets14;
  TH2D* hPtvsTracklets10;

  AliAnalysisTaskDataSpeedOfSoundSim(
      const AliAnalysisTaskDataSpeedOfSoundSim&);  // not implemented
  AliAnalysisTaskDataSpeedOfSoundSim& operator=(
      const AliAnalysisTaskDataSpeedOfSoundSim&);  // not implemented

  ClassDef(AliAnalysisTaskDataSpeedOfSoundSim, 3);
};
#endif
