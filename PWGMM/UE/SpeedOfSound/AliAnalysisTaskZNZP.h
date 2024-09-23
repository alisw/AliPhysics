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
  void SetV0Mmin(double V0Mmin) { fV0Mmin = V0Mmin; }  // Set V0M min value
  void SetV0Mmax(double V0Mmax) { fV0Mmax = V0Mmax; }  // Set V0M max value
  void SetUseMC(bool mc = false) { fUseMC = mc; }      // use to analyse MC data
  void SetEtaCut(const double& etacut) { fEtaCut = etacut; }
  void SetPtCut(const double& ptmin) { fPtMin = ptmin; }
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
  void SaveAsymmetry(bool sv) { fSaveAsy = sv; };
  bool MeanSigmaZN(double& mean, double& sigma, const std::string& ZN);
  int CentBin();

 protected:
 private:
  AliESDEvent* fESD;
  AliEventCuts fEventCuts;
  AliStack* fMCStack;
  AliMCEvent* fMC;
  bool fUseMC;
  bool fIsSystematics;
  bool fVaryVtxZPos;
  bool fSaveAsy;
  float fMinVtxZPos;
  float fMaxVtxZPos;
  int fSystematic;
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

  TProfile* pV0MAmpChannel;
  TH1F* hV0Percentile;
  TH1F* hBestVtxZ;
  TH2F* hZNvsV0MPer;
  TH2F* hZNvsV0M;
  TH2F* hZNAvsV0M;
  TH2F* hZNCvsV0M;
  TH2F* hAsyN;
  TH2F* hZPvsV0MPer;
  TH2F* hZPvsV0M;
  TH2F* hZPAvsV0M;
  TH2F* hZPCvsV0M;
  /*TH2F* hAsyP;*/
  TH1F* hZNCpmc;
  TH1F* hZNApmc;
  TH1F* hZPCpmc;
  TH1F* hZPApmc;
  TH2F* hZNCNorm;
  TH2F* hZNANorm;
  /*TH2F* hZNCNormSca;*/
  /*TH2F* hZNANormSca;*/

  AliAnalysisTaskZNZP(const AliAnalysisTaskZNZP&);  // not implemented
  AliAnalysisTaskZNZP& operator=(
      const AliAnalysisTaskZNZP&);  // not implemented

  ClassDef(AliAnalysisTaskZNZP, 5);
};
#endif
