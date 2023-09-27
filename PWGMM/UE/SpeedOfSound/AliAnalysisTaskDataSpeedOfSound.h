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
class TProfile;
// class THnSparse;

#include "AliAnalysisTaskSE.h"
#include "AliEventCuts.h"
#include "AliGenEventHeader.h"
#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliMCParticle.h"
#include "AliMultSelection.h"
#include "AliStack.h"
#include "TParticle.h"

class AliAnalysisTaskDataSpeedOfSound : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskDataSpeedOfSound();
  AliAnalysisTaskDataSpeedOfSound(const char* name);
  virtual ~AliAnalysisTaskDataSpeedOfSound();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t* option);
  virtual void Terminate(Option_t* option);
  void AnalyzeRecEvent(int& rec_nch, std::vector<float>& vec_rec_pt) const;
  void DCAxyDistributions() const;
  void TrackingEfficiency() const;
  void AnalyzeMCevent(int& true_nch, std::vector<float>& vec_true_pt) const;
  void DetectorResponse(const int& true_nch, const int& rec_nch) const;
  void RecMultiplicityDistributions(const int& rec_nch,
                                    const std::vector<float>& vec_rec_pt) const;
  void TrueMultiplicityDistributions(
      const int& true_nch, const std::vector<float>& vec_true_pt) const;
  void MultiplicityDistributions(const int& rec_nch,
                                 const std::vector<float>& vec_rec_pt) const;
  void SetV0Mmin(double V0Mmin) { fV0Mmin = V0Mmin; }  // Set V0M min value
  void SetV0Mmax(double V0Mmax) { fV0Mmax = V0Mmax; }  // Set V0M max value
  void SetUseMC(bool mc = false) { fUseMC = mc; }      // use to analyse MC data
  void SetTrackCuts(bool TPConly = true) { fIsTPConly = TPConly; }
  bool HasRecVertex();

 protected:
 private:
  AliESDEvent* fESD;
  AliEventCuts fEventCuts;
  AliStack* fMCStack;
  AliMCEvent* fMC;
  bool fUseMC;
  bool fIsTPConly;
  AliAnalysisFilter* fTrackFilter;
  AliAnalysisFilter* fTrackFilterwoDCA;
  TList* fOutputList;
  double fEtaCut;
  double fPtMin;
  double fV0Mmin;
  double fV0Mmax;
  double ftrackmult08;
  double fv0mpercentile;
  // double fv0mpercentilebefvtx;
  float fdcaxy;
  float fdcaz;
  AliMultSelection* fMultSelection;
  // AliMultSelection* fMultSelectionbefvtx;
  TH1F* hRefMult;
  TH2F* hNchUCvsV0M;
  TH1F* hV0Mmult;
  TH2F* hTrackletvsV0M;
  TH2F* hPtvsNch[5];
  TH1F* hTrueVtxZ;
  TH2F* hTrueNchvsV0M_UC;
  TH2F* hRecNchvsV0M_UC;
  // TH3F* hTrueNchvsTruePt;
  TH2F* hFullNchResponse;
  TH2F* hNchResponse;
  TH2F* hPtTruePrivsV0M;
  // TH2F* hPtTrueSecvsV0M;
  // TH2F* hPtTrueAllvsV0M;
  TH2F* hPtRecPrivsV0M;
  TH2F* hRecNchvsRecPt[5];
  TH2F* hTrueNchvsTruePt[5];
  TH2F* hDCAxyPri[6];
  TH2F* hDCAxyWeDe[6];
  TH2F* hDCAxyMaIn[6];
  TH2F* hDCAxyData[6];
  TH1F* hTrueNch;
  TH2F* hTrueNchHM;
  TH1F* hTrueNchWithTrigger;
  TH2F* hTrueNchHMWithTrigger;
  TH1F* hTrueNchWithEventCuts;
  TH2F* hTrueNchHMWithEventCuts;
  TH1F* hTrueNchWithVtxSel;
  TH2F* hTrueNchHMWithVtxSel;

  AliAnalysisTaskDataSpeedOfSound(
      const AliAnalysisTaskDataSpeedOfSound&);  // not implemented
  AliAnalysisTaskDataSpeedOfSound& operator=(
      const AliAnalysisTaskDataSpeedOfSound&);  // not implemented

  ClassDef(AliAnalysisTaskDataSpeedOfSound, 3);
};
#endif
