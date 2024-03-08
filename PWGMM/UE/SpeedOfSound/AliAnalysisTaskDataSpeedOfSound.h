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
  void VertexPosition();
  void GetZDCCentrality();
  void SetV0Mmin(double V0Mmin) { fV0Mmin = V0Mmin; }  // Set V0M min value
  void SetV0Mmax(double V0Mmax) { fV0Mmax = V0Mmax; }  // Set V0M max value
  void SetHMCut(double HMcut) { fHMCut = HMcut; }      // Set V0M max value
  void SetUseMC(bool mc = false) { fUseMC = mc; }      // use to analyse MC data
  void SetEtaCut(const double& etacut, const double& etamin_halftpc,
                 const double& etamax_halftpc, const double& etamin_spd,
                 const double& etamax_spd, const double& etamin_tpc,
                 const double& etamax_tpc, const double& etacut_4ptwspd,
                 const double& etacut_4ptwtpc) {
    fEtaCut = etacut;
    fEtaCutHalfTPCMin = etamin_halftpc;
    fEtaCutHalfTPCMax = etamax_halftpc;
    fEtaCutSPDGapMin = etamin_spd;
    fEtaCutSPDGapMax = etamax_spd;
    fEtaCutTPCGapMin = etamin_tpc;
    fEtaCutTPCGapMax = etamax_tpc;
    fEtaCutForpTwSPDGap = etacut_4ptwspd;
    fEtaCutForpTwTPCGap = etacut_4ptwtpc;
  }
  void SetPtCut(const double& ptmin, const double& ptmin_cent,
                const double& ptmax_cent) {
    fPtMin = ptmin;
    fPtMinCent = ptmin_cent;
    fPtMaxCent = ptmax_cent;
  }
  void SetTrigger(UInt_t trigger = AliVEvent::kINT7) { fTrigger = trigger; }
  bool HasRecVertex();
  void SetSystematics(bool issystematics = true, int systematic = 1) {
    fIsSystematics = issystematics;
    fSystematic = systematic;
  }
  void SetSystematicsVtxZ(bool varyZpos = false, const float& minz = -5.0,
                          const float& maxz = 5.0) {
    fVaryVtxZPos = varyZpos;
    fMinVtxZPos = minz;
    fMaxVtxZPos = maxz;
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
  bool fVaryVtxZPos;
  float fMinVtxZPos;
  float fMaxVtxZPos;
  int fSystematic;
  UInt_t fTrigger;
  AliAnalysisFilter* fTrackFilter;
  AliAnalysisFilter* fTrackFilterwoDCA;
  TList* fOutputList;
  double fEtaCut;
  double fEtaCutHalfTPCMin;
  double fEtaCutHalfTPCMax;
  double fEtaCutForpTwSPDGap;
  double fEtaCutSPDGapMin;
  double fEtaCutSPDGapMax;
  double fEtaCutForpTwTPCGap;
  double fEtaCutTPCGapMin;
  double fEtaCutTPCGapMax;
  double fPtMin;
  double fPtMinCent;
  double fPtMaxCent;
  double fV0Mmin;
  double fV0Mmax;
  double fHMCut;
  double ftrackmult08;
  double fv0mpercentile;
  float fv0mamplitude;
  int fTracklets14;
  int fTracklets10;
  int fTrackletsEtaGap;
  int fTracksEtaGapTPC;
  double fza;
  double fzc;
  double fzn;
  float fdcaxy;
  float fdcaz;
  AliMultSelection* fMultSelection;
  TH1F* hNch;
  TH2F* hNchvsV0MAmp;
  TH2F* hV0MvsV0MAmp;
  TProfile* pV0MAmpChannel;
  TH1F* hV0MAmplitude;
  TH1F* hV0Percentile;
  TH2D* hPtvsNch;
  TProfile* pPtvsNch;
  TProfile* pPtEtaNegvsNchEtaPos;
  TProfile* pPtEtaPosvsNchEtaNeg;
  TProfile* pPtvsV0MAmp;
  TH2D* hPtvsV0MAmp;
  TH1F* hNchEtaPos;
  TH1F* hNchEtaNeg;
  TH2D* hPtEtaNegvsNchEtaPos;
  TH2D* hPtEtaPosvsNchEtaNeg;
  TH2F* hDCAxyData[1];
  TH2F* hPhiEtaSPD;
  TH2F* hPhiEtaGapTPC;
  TH1F* hPtWithCutForCent;
  TH2F* hPhiEtaPosHalfTPC;
  TH2F* hPhiEtaNegHalfTPC;
  TH2F* hPhiEtaGapSPD;
  TH1F* hBestVtxZ;
  TH1F* hTracklets14;
  TH1F* hTracklets10;
  TH1F* hTrackletsEtaGap;
  TH1F* hTracksEtaGapTPC;
  TH2D* hPtvsTracklets14;
  TH2D* hPtvsTracklets10;
  TH2D* hPtvsTrackletsEtaGap;
  TH2D* hPtvsTracksEtaGapTPC;
  TProfile* pPtvsTracklets14;
  TProfile* pPtvsTracklets10;
  TProfile* pPtvsTrackletsEtaGap;
  TProfile* pPtvsTracksEtaGapTPC;

  AliAnalysisTaskDataSpeedOfSound(
      const AliAnalysisTaskDataSpeedOfSound&);  // not implemented
  AliAnalysisTaskDataSpeedOfSound& operator=(
      const AliAnalysisTaskDataSpeedOfSound&);  // not implemented

  ClassDef(AliAnalysisTaskDataSpeedOfSound, 3);
};
#endif
