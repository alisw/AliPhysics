/* Copyright(c) 1998-1999, ALICfalseriment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* Add a description of your MPI analysis */

#ifndef AliAnalysisTaskZNZP_H
#define AliAnalysisTaskZNZP_H

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
#include <string>

class AliAnalysisTaskZNZP : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskZNZP();
  AliAnalysisTaskZNZP(const char* name);
  virtual ~AliAnalysisTaskZNZP();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t* option);
  virtual void Terminate(Option_t* option);
  void DCAxyDistributions() const;
  void GetSPDMultiplicity();
  void GetCalibratedV0Amplitude();
  void VertexPosition();
  void GetZDC();
  void SetPeriod(std::string period) { fPeriod = period; }
  void SetV0Mmin(double V0Mmin) { fV0Mmin = V0Mmin; }  // Set V0M min value
  void SetV0Mmax(double V0Mmax) { fV0Mmax = V0Mmax; }  // Set V0M max value
  void SetUseMC(bool mc = false) { fUseMC = mc; }      // use to analyse MC data
  void SetEtaCut(const double& etacut) { fEtaCut = etacut; }
  void SetPtCut(const double& ptmin) { fPtMin = ptmin; }
  void SetTrigger(UInt_t trigger = AliVEvent::kINT7) { fTrigger = trigger; }
  bool HasRecVertex();
  void SPDActivity();
  void TPCActivity();
  void IsTowerEnergy(bool sv) { fTowerEnergy = sv; };
  /*bool MeanSigmaZN(double& mean, double& sigma, const std::string& ZN);*/
  /*int CentBin();*/

 protected:
 private:
  AliESDEvent* fESD;
  AliEventCuts fEventCuts;
  AliStack* fMCStack;
  AliMCEvent* fMC;
  bool fUseMC;
  std::string fPeriod;
  bool fTowerEnergy;
  UInt_t fTrigger;
  AliMultSelection* fMultSelection;
  AliAnalysisFilter* fTrackFilter;
  AliAnalysisFilter* fTrackFilterwoDCA;
  TList* fOutputList;
  double fEtaCut;
  double fPtMin;
  double fV0Mmin;
  double fV0Mmax;
  double ftrackmult08;
  double fv0mpercentile;
  float fv0mamplitude;
  int fSPD;
  int fNchTPC;
  double fET;
  TF1* fZNCvsV0M;
  TF1* fZNAvsV0M;

  TProfile* pV0MAmpChannel;
  TH1F* hV0Percentile;
  TH1F* hBestVtxZ;
  TH2F* hZNAvsV0M;
  TH2F* hZNCvsV0M;
  TH2F* hAsyN;
  TH2F* hZPAvsV0M;
  TH2F* hZPCvsV0M;
  TH2F* hZNCNorm;
  TH2F* hZNANorm;
  TProfile* pZNChannel;
  TProfile* pZPChannel;
  TProfile* pZNCvsV0Amp;
  TProfile* pZNAvsV0Amp;
  TProfile* pZPCvsV0Amp;
  TProfile* pZPAvsV0Amp;
  TProfile* pZNCvsNch;
  TProfile* pZNAvsNch;
  TProfile* pZPCvsNch;
  TProfile* pZPAvsNch;
  TProfile* pZNCvsEt;
  TProfile* pZNAvsEt;
  TProfile* pZPCvsEt;
  TProfile* pZPAvsEt;
  TProfile* pZNCvsSPD;
  TProfile* pZNAvsSPD;
  TProfile* pZPCvsSPD;
  TProfile* pZPAvsSPD;
  TH2F* hZNCTowvsNCEn;
  TH2F* hZPATowvsPAEn;

  AliAnalysisTaskZNZP(const AliAnalysisTaskZNZP&);  // not implemented
  AliAnalysisTaskZNZP& operator=(
      const AliAnalysisTaskZNZP&);  // not implemented

  ClassDef(AliAnalysisTaskZNZP, 5);
};
#endif
