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
  void GetZDC();
  void SetV0Mmin(double V0Mmin) { fV0Mmin = V0Mmin; }  // Set V0M min value
  void SetV0Mmax(double V0Mmax) { fV0Mmax = V0Mmax; }  // Set V0M max value
  void SetHMCut(double HMcut) { fHMCut = HMcut; }      // Set V0M max value
  void SetUseMC(bool mc = false) { fUseMC = mc; }      // use to analyse MC data
  void SetEtaCut(const double& etacut, const double& etamin_spd,
                 const double& etamax_spd, const double& etamin_tpc,
                 const double& etamax_tpc, const double& etacut_4ptwspd,
                 const double& etacut_4ptwtpc) {
    fEtaCut = etacut;
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
  void SetSPDVtxZ(double cut) { fSPDVtxCut = cut; }

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
  int fTrackletsEtaGap;
  int fTracksEtaGapTPC;
  AliMultSelection* fMultSelection;
  TH1F* hNch;
  TH2F* hNchvsV0MAmp;
  TH2F* hV0MvsV0MAmp;
  TProfile* pV0MAmpChannel;
  TH1F* hV0MAmplitude;
  TH1F* hV0Percentile;
  TH2D* hPtvsNch;
  TProfile* pPtvsNch;
  TProfile* pPtvsV0MAmp;
  TH2D* hPtvsV0MAmp;
  TH2F* hDCAxyData[1];
  TH2F* hPhiEtaGapTPC;
  TH1F* hPtWithCutForCent;
  TH2F* hPhiEtaGapSPD;
  TH1F* hBestVtxZ;
  TH1F* hTrackletsEtaGap;
  TH1F* hTracksEtaGapTPC;
  TH2D* hPtvsTrackletsEtaGap;
  TH2D* hPtvsTracksEtaGapTPC;
  TProfile* pPtvsTrackletsEtaGap;
  TProfile* pPtvsTracksEtaGapTPC;
  TH1F* hSPDFull;
  TH1F* hSPDEtaAdj;
  TH1F* hSPDEtaGapW;
  TH1F* hSPDEtaGapWW;
  TH1F* hEtFull;
  TH1F* hEtEtaGap;
  TH2D* hPtvsSPDFull;
  TH2D* hPtvsSPDEtaAdj;
  TH2D* hPtvsSPDEtaGapW;
  TH2D* hPtvsSPDEtaGapWW;
  TH2D* hPtvsEtFull;
  TH2D* hPtvsEtEtaGap;
  TH2F* hPtvsTPCEtaGapWidepT;
  TH2F* hPtvsEtEtaGapWidepT;
  TProfile* pZVtxvsSPDClus;
  TProfile* pSPDClusvsEta;
  double fSPDVtxCut;
  int fSPDFull;
  int fSPDEtaAdj;
  int fSPDEtaGapW;
  int fSPDEtaGapWW;
  double fZN;
  double fZP;
  double fZDC;
  TProfile* pZDCvsV0MAmp;
  TProfile* pZDCvsTPCFull;
  TProfile* pZDCvsTPCEtaGap;
  TProfile* pZDCvsSPDFull;
  TProfile* pZDCvsSPDEtaGap;
  TProfile* pZDCvsSPDEtaAdj;
  TProfile* pZDCvsSPDEtaGapW;
  TProfile* pZDCvsEtFull;
  TProfile* pZDCvsEtEtaGap;
  TH2F* hPtvsEtFullWidepT;
  TH2F* hPtvsSPDFullWidepT;
  TH2F* hPtvsTPCFullWidepT;
  TH2F* hPtvsSPDEtaGapWidepT;
  TH2F* hPtvsSPDEtaGapWWidepT;
  TProfile* pZNvsV0MAmp;
  TProfile* pZNvsTPCFull;
  TProfile* pZNvsTPCEtaGap;
  TProfile* pZNvsSPDFull;
  TProfile* pZNvsSPDEtaGap;
  TProfile* pZNvsSPDEtaAdj;
  TProfile* pZNvsSPDEtaGapW;
  TProfile* pZNvsEtFull;
  TProfile* pZNvsEtEtaGap;
  TProfile* pZPvsV0MAmp;
  TProfile* pZPvsTPCFull;
  TProfile* pZPvsTPCEtaGap;
  TProfile* pZPvsSPDFull;
  TProfile* pZPvsSPDEtaGap;
  TProfile* pZPvsSPDEtaAdj;
  TProfile* pZPvsSPDEtaGapW;
  TProfile* pZPvsEtFull;
  TProfile* pZPvsEtEtaGap;

  AliAnalysisTaskDataSpeedOfSound(
      const AliAnalysisTaskDataSpeedOfSound&);  // not implemented
  AliAnalysisTaskDataSpeedOfSound& operator=(
      const AliAnalysisTaskDataSpeedOfSound&);  // not implemented

  ClassDef(AliAnalysisTaskDataSpeedOfSound, 5);
};
#endif
