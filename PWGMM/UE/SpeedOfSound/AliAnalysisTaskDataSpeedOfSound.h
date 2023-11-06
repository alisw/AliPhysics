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
  void AnalyzeRecEvent(int& rec_nch, int& rec_nch_neg_eta, int& rec_nch_pos_eta,
                       std::vector<float>& vec_rec_pt,
                       std::vector<float>& vec_rec_pt_neg_eta,
                       std::vector<float>& vec_rec_pt_pos_eta) const;
  void DCAxyDistributions() const;
  void TrackingEfficiency() const;
  void AnalyzeMCevent(int& true_nch, int& true_nch_v0, int& true_nch_neg_eta,
                      int& true_nch_pos_eta, std::vector<float>& vec_true_pt,
                      std::vector<float>& vec_true_pt_neg_eta,
                      std::vector<float>& vec_true_pt_pos_eta) const;
  void DetectorResponse(const int& true_nch, const int& rec_nch) const;
  void RecMultiplicityDistributions(
      const int& rec_nch, const int& rec_nch_neg_eta,
      const int& rec_nch_pos_eta, const std::vector<float>& vec_rec_pt,
      const std::vector<float>& vec_rec_pt_neg_eta,
      const std::vector<float>& vec_rec_pt_pos_eta) const;
  void TrueMultiplicityDistributions(
      const int& true_nch, const int& true_nch_v0,
      const std::vector<float>& vec_true_pt) const;
  void TrueMultiplicityDistributions(
      const int& true_nch, const int& true_nch_v0, const int& true_nch_neg_eta,
      const int& true_nch_pos_eta, const std::vector<float>& vec_true_pt,
      const std::vector<float>& vec_true_pt_neg_eta,
      const std::vector<float>& vec_true_pt_pos_eta) const;
  void MultiplicityDistributions(
      const int& rec_nch, const int& rec_nch_neg_eta,
      const int& rec_nch_pos_eta, const std::vector<float>& vec_rec_pt,
      const std::vector<float>& vec_rec_pt_neg_eta,
      const std::vector<float>& vec_rec_pt_pos_eta) const;
  void GetCalibratedV0Amplitude();
  void SetV0Mmin(double V0Mmin) { fV0Mmin = V0Mmin; }  // Set V0M min value
  void SetV0Mmax(double V0Mmax) { fV0Mmax = V0Mmax; }  // Set V0M max value
  void SetHMCut(double HMcut) { fHMCut = HMcut; }      // Set V0M max value
  void SetUseMC(bool mc = false) { fUseMC = mc; }      // use to analyse MC data
  void SetUseZDC(bool zdc = false) { fUseZDC = zdc; }  // use ZDC selection
  void SetEtaCut(const double& etacut) { fEtaCut = etacut; }
  void SetEtaMinCut(const double& etamin) { fEtaMin = etamin; }
  void SetEtaMaxCut(const double& etamax) { fEtaMax = etamax; }
  void SetPtMin(const double& ptmin) { fPtMin = ptmin; }
  void SetTrackCuts(bool TPConly = true) { fIsTPConly = TPConly; }
  void SetTrigger(UInt_t offlineTriggerMask = AliVEvent::kINT7) {
    fTrigger = offlineTriggerMask;
  }
  bool HasRecVertex();
  void ZDC(const int& rec_nch, const std::vector<float>& vec_rec_pt) const;

 protected:
 private:
  AliESDEvent* fESD;
  AliEventCuts fEventCuts;
  AliStack* fMCStack;
  AliMCEvent* fMC;
  bool fUseZDC;
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
  double ftrackmult08;
  double fv0mpercentile;
  float fv0mamplitude;
  float fdcaxy;
  float fdcaz;
  AliMultSelection* fMultSelection;
  TH2D* hNchvsV0M;
  TH2D* hNchvsV0MAmp;
  TH2D* hV0MvsV0MAmp;
  TProfile* pV0MAmpChannel;
  TH1D* hV0MAmplitude;
  TH1F* hV0Mmult;
  TProfile* pPtvsNch;
  TProfile* pPtEtaNegvsNchEtaPos;
  TProfile* pPtEtaPosvsNchEtaNeg;
  TProfile* pPtvsV0MAmp;
  TH3D* hPtvsNchvsV0MAmp;
  TH2D* hNchEtaPosvsNchEtaNeg;
  TH3D* hPtEtaNegvsNchEtaPosvsV0MAmp;
  TH3D* hPtEtaPosvsNchEtaNegvsV0MAmp;
  TH1F* hTrueVtxZ;
  TH3F* hTrueNchvsTrueV0MAmp;
  TH3F* hRecNchvsRecV0MAmp;
  TH2D* hNchResponse;
  TH2D* hPtTruePrivsV0M;
  TH2D* hPtRecPrivsV0M;
  TH3D* hRecPtvsRecNchvsRecV0MAmp;
  TH2D* hNchEtaPosvsNchEtaNeg_MCRec;
  TH3D* hPtEtaNegvsNchEtaPosvsV0MAmp_MCRec;
  TH3D* hPtEtaPosvsNchEtaNegvsV0MAmp_MCRec;
  TH3D* hTruePtvsTrueNchvsTrueV0MAmp;
  TH2D* hNchEtaPosvsNchEtaNeg_MCTrue;
  TH3D* hPtEtaNegvsNchEtaPosvsV0MAmp_MCTrue;
  TH3D* hPtEtaPosvsNchEtaNegvsV0MAmp_MCTrue;
  TH2F* hDCAxyPri[1];
  TH2F* hDCAxyWeDe[1];
  TH2F* hDCAxyMaIn[1];
  TH2F* hDCAxyData[1];
  TH2F* hAllpTRec;
  TH2F* hAllpTTrue;
  TH2F* hPripTRec;
  TH2F* hPripTTrue;
  TH2F* hTrueNchHM;
  TH2F* hTrueNchHMWithTrigger;
  TH2F* hTrueNchHMWithEventCuts;
  TH2F* hTrueNchHMWithVtxSel;
  TH2D* hZNCvsZNA;
  TH2D* hZPCvsZPA;
  TH2D* hZEM;
  TH2D* hZNvsNch;
  TH2D* hZAvsNchHM;
  TH2D* hZCvsNchHM;
  TH2D* hZNvsNchHM;
  TH2D* hZNvsV0MAmp;
  TH2D* hZNvsV0MAmpHM;
  TH3D* hPtvsNchvsZAHM;
  TH3D* hPtvsNchvsZCHM;
  TH3D* hPtvsNchvsZNHM;
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
