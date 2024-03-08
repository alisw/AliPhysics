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
  void VertexPosition();
  void DCAxyDistributions();
  void TrackingEfficiency();
  void DCA();
  int GetPidCode(int pdgCode) const;
  void ReadMCEvent();
  void TrueMultiplicityDistributions();
  void MultiplicityDistributions();
  void GetTPCMultiplicity();
  void GetCalibratedV0Amplitude();
  void SetV0Mmin(double V0Mmin) { fV0Mmin = V0Mmin; }  // Set V0M min value
  void SetV0Mmax(double V0Mmax) { fV0Mmax = V0Mmax; }  // Set V0M max value
  void SetHMCut(double HMcut) { fHMCut = HMcut; }      // Set V0M max value
  void SetRandomNumberCut(double rdcut) { fRandomNumberCut = rdcut; }
  void SetUseMC(bool mc = false) { fUseMC = mc; }  // use to analyse MC data
  void SetEtaCut(const double& etacut, const double& etamin,
                 const double& etamax, const double& min_spd,
                 const double& max_spd, const double& min_tpc,
                 const double& max_tpc, const double& eta_cut4ptwspd,
                 const double& eta_cut4ptwtpc) {
    fEtaCut = etacut;
    fEtaMin = etamin;
    fEtaMax = etamax;
    fEtaGapSPDNchMin = min_spd;
    fEtaGapSPDNchMax = max_spd;
    fEtaGapTPCNchMin = min_tpc;
    fEtaGapTPCNchMax = max_tpc;
    fEtaGapSPDpT = eta_cut4ptwspd;
    fEtaGapTPCpT = eta_cut4ptwtpc;
  }
  void SetPtCut(const double& ptmin, const double& ptmin_cent,
                const double& ptmax_cent) {
    fPtMin = ptmin;
    fPtMinCent = ptmin_cent;
    fPtMaxCent = ptmax_cent;
  }
  void SetTrigger(UInt_t offlineTriggerMask = AliVEvent::kINT7) {
    fTrigger = offlineTriggerMask;
  }
  bool HasRecVertex();
  void GetSPDMultiplicity();
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
  double fEtaMin;
  double fEtaMax;
  double fEtaGapSPDpT;
  double fEtaGapSPDNchMin;
  double fEtaGapSPDNchMax;
  double fEtaGapTPCpT;
  double fPtMinCent;
  double fPtMaxCent;
  double fEtaGapTPCNchMin;
  double fEtaGapTPCNchMax;
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
  int fTrueNchEtaGap;
  int fTrueNchEtaGapTPC;
  int fTrueNchEtaPos;
  int fTrueNchEtaNeg;
  int fTrueV0;
  int fTracklets14;
  int fTracklets10;
  int fTrackletsEtaGap;
  int fTracksEtaGapTPC;
  int fTracksEtaNeg;
  int fTracksEtaPos;
  int fTracksEta08;
  float fdcaxy;
  float fdcaz;
  AliMultSelection* fMultSelection;
  TH1F* hBestVtxZ;
  TH2D* hNchvsV0M;
  TH2D* hNchvsV0MAmp;
  TH2D* hV0MvsV0MAmp;
  TProfile* pV0MAmpChannel;
  TH1D* hV0MAmplitude;
  TH1F* hCounter;
  TH1F* hV0Mmult;
  TH2D* hPtvsV0MAmp;
  TH2D* hNchEtaPosvsNchEtaNeg;
  TH2D* hPtEtaNegvsNchEtaPos;
  TH2D* hPtEtaPosvsNchEtaNeg;
  TH1F* hTrueVtxZ;
  TH1F* hTrueNch;
  TH1F* hTrueV0MAmp;
  TH2D* hTruePtvsTrueNch;
  TH2D* hTrueNchEtaPosvsTrueNchEtaNeg;
  TH2D* hTruePtEtaNegvsTrueNchEtaPos;
  TH2D* hTruePtEtaPosvsTrueNchEtaNeg;
  TH1F* hTrueNchEtaGap;
  TH1F* hTrueNchEtaGapTPC;
  TH2D* hTruePtvsTrueNchEtaGap;
  TH2D* hTruePtvsTrueNchEtaGapTPC;
  TH2D* hTrueNchvsMeasNchEtaGapSPD;
  TH2D* hTrueNchvsMeasNchEtaGapTPC;
  TH2D* hTrueNchvsMeasNchHalfTPC;
  TH1F* hPtOutAll_ch;
  TH2F* hPhiEta;
  TH2F* hPhiEtaGap_SPD;
  TH2F* hPhiEtaGap_TPC;
  TH1F* hTrackletsEtaGap;
  TH1F* hTracksEtaGapTPC;
  TH2D* hPtvsTrackletsEtaGap;
  TH2D* hPtvsTracksEtaGapTPC;
  TH1F* hPtInPrim[8];
  TH2D* hPtInPrim_EtaGap_SPD[8];
  TH2D* hPtInPrim_EtaGap_TPC[8];
  TH2D* hPtInPrim_HalfEta[8];
  TH1F* hPtOutPrim[8];
  TH2D* hPtOutPrim_EtaGap_SPD[8];
  TH2D* hPtOutPrim_EtaGap_TPC[8];
  TH2D* hPtOutPrim_HalfEta[8];
  TH2F* hDCAxy[3];

  AliAnalysisTaskDataSpeedOfSoundSim(
      const AliAnalysisTaskDataSpeedOfSoundSim&);  // not implemented
  AliAnalysisTaskDataSpeedOfSoundSim& operator=(
      const AliAnalysisTaskDataSpeedOfSoundSim&);  // not implemented

  ClassDef(AliAnalysisTaskDataSpeedOfSoundSim, 3);
};
#endif
