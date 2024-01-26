/* Copyright(c) 1998-1999, ALICfalseriment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* Add a description of your MPI analysis */

#ifndef AliAnalysisTaskDataSpeedOfSound_H
#define AliAnalysisTaskDataSpeedOfSound_H

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

class AliAnalysisTaskDataSpeedOfSound : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskDataSpeedOfSound();
  AliAnalysisTaskDataSpeedOfSound(const char* name);
  virtual ~AliAnalysisTaskDataSpeedOfSound();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t* option);
  virtual void Terminate(Option_t* option);
  void DCAxyDistributions() const;
  void GetSPDMultiplicity();
  void MultiplicityDistributions();
  void GetCalibratedV0Amplitude();
  void GetZDCCentrality();
  void SetV0Mmin(double V0Mmin) { fV0Mmin = V0Mmin; }  // Set V0M min value
  void SetV0Mmax(double V0Mmax) { fV0Mmax = V0Mmax; }  // Set V0M max value
  void SetHMCut(double HMcut) { fHMCut = HMcut; }      // Set V0M max value
  void SetUseMC(bool mc = false) { fUseMC = mc; }      // use to analyse MC data
  void SetEtaCut(const double& etacut) { fEtaCut = etacut; }
  void SetEtaMinCut(const double& etamin) { fEtaMin = etamin; }
  void SetEtaMaxCut(const double& etamax) { fEtaMax = etamax; }
  void SetEtaGappT(const double& eta) { fEtaGappT = eta; }
  void SetEtaGapNch(const double& etamin, const double& etamax) {
    fEtaGapNchMin = etamin;
    fEtaGapNchMax = etamax;
  }
  void SetPtMin(const double& ptmin) { fPtMin = ptmin; }
  void SetTrigger(UInt_t trigger = AliVEvent::kINT7) { fTrigger = trigger; }
  bool HasRecVertex();
  void SetSystematics(bool issystematics = true, int systematic = 1) {
    fIsSystematics = issystematics;
    fSystematic = systematic;
  }
  void ChangeCut(AliESDtrackCuts* fCuts);

 protected:
 private:
  AliESDEvent* fESD;
  AliEventCuts fEventCuts;
  AliStack* fMCStack;
  AliMCEvent* fMC;
  bool fUseMC;
  bool fIsSystematics;
  int fSystematic;
  UInt_t fTrigger;
  AliAnalysisFilter* fTrackFilter;
  AliAnalysisFilter* fTrackFilterwoDCA;
  TList* fOutputList;
  double fEtaCut;
  double fEtaMin;
  double fEtaMax;
  double fEtaGappT;
  double fEtaGapNchMin;
  double fEtaGapNchMax;
  double fPtMin;
  double fV0Mmin;
  double fV0Mmax;
  double fHMCut;
  double ftrackmult08;
  double fv0mpercentile;
  float fv0mamplitude;
  int fTracklets14;
  int fTracklets10;
  int fTrackletsEtaGap;
  double fza;
  double fzc;
  double fzn;
  float fdcaxy;
  float fdcaz;
  AliMultSelection* fMultSelection;
  TH2F* hNchvsV0M;
  TH2F* hNchvsV0MAmp;
  TH2F* hV0MvsV0MAmp;
  TProfile* pV0MAmpChannel;
  TH1F* hV0MAmplitude;
  TH1F* hV0Mmult;
  TProfile* pPtvsNch;
  TProfile* pPtEtaNegvsNchEtaPos;
  TProfile* pPtEtaPosvsNchEtaNeg;
  TProfile* pPtvsV0MAmp;
  TH2D* hPtvsV0MAmp;
  TH2D* hNchEtaPosvsNchEtaNeg;
  TH2D* hPtEtaNegvsNchEtaPos;
  TH2D* hPtEtaPosvsNchEtaNeg;
  TH2F* hDCAxyData[1];
  TH2F* hZAvsNchHM;
  TH2F* hZCvsNchHM;
  TH2F* hZNvsNchHM;
  TH2F* hZNvsV0MAmpHM;
  TH2D* hPtvsZAHM;
  TH2D* hPtvsZCHM;
  TH2D* hPtvsZNHM;
  TH2F* hPhiEtaSPD;
  TH1F* hEtaGapSPD;
  TH2F* hVtxZvsTracklets;
  TH2F* hTrackletsvsV0MAmp14;
  TH2F* hTrackletsvsV0MAmp10;
  TH1F* hTrackletsEtaGap;
  TH2D* hPtvsTracklets14;
  TH2D* hPtvsTracklets10;
  TH2D* hPtvsTrackletsEtaGap;
  TProfile* pPtvsTracklets14;
  TProfile* pPtvsTracklets10;
  TProfile* pPtvsTrackletsEtaGap;
  TProfile* pPtvsZA;
  TProfile* pPtvsZC;
  TProfile* pPtvsZN;

  AliAnalysisTaskDataSpeedOfSound(
      const AliAnalysisTaskDataSpeedOfSound&);  // not implemented
  AliAnalysisTaskDataSpeedOfSound& operator=(
      const AliAnalysisTaskDataSpeedOfSound&);  // not implemented

  ClassDef(AliAnalysisTaskDataSpeedOfSound, 3);
};
#endif
