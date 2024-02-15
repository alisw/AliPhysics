#ifndef AliAnalysisTaskThreePartCorr_h
#define AliAnalysisTaskThreePartCorr_h

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// Analysis task to produce trees of lightweight events
// evgeny.kryshen@cern.ch

#include "AliAnalysisTaskSE.h"
#include <iostream>
#include <fstream>
using namespace std;

const Double_t massLambda = 1.115683;
const Double_t DGaussSigma = 0.004759;

class TList;
class AliESDEvent;
class AliAODEvent;
class AliESDtrack;
class AliAODTrack;
class AliAODv0;
class AliAODcascade;
class AliAnalysisUtils;
class AliPIDResponse;
class TH1I;
class TH1D;
class TH2D;
class AliVVertex;
class AliEventPoolManager;

class AliAnalysisTaskThreePartCorr : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskThreePartCorr();
  AliAnalysisTaskThreePartCorr(const char *name);
  virtual ~AliAnalysisTaskThreePartCorr();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);

  Bool_t V0TriggerCuts(AliAODv0 *V0, Double_t PrimVertex[3]);
  Bool_t XiTriggerCuts(AliAODcascade *Xi);
  Bool_t AssociateCuts(AliAODTrack *Track);
  Double_t PhiCorrection(Double_t TriggerPhi, Double_t AssociatePhi);

  AliPID::EParticleType TrackPID(AliAODTrack *Track, Bool_t FillHist = kFALSE);
  TObjArray *CloneAndReduceTrackList(TObjArray *Tracks);
  
 protected:
  AliAnalysisTaskThreePartCorr(const AliAnalysisTaskThreePartCorr &task);
  AliAnalysisTaskThreePartCorr& operator=(const AliAnalysisTaskThreePartCorr &task);

  AliESDEvent *fESD;
  AliAODEvent *fAOD;
  AliAnalysisUtils *fUtils;      //! analysis utils to detect pileup
  AliPIDResponse *fPIDResponse;  //! points to class for PID
  TList *fListOfHistos;          //! list of output histograms
  TList *fQAList, *fCorrList;    //! lists of correlation and QA histograms
  
  TH1I *hEventStatistics;        //! cut-by-cut counter of events

  TAxis *Centaxis, *Zvtxaxis, *Triggptaxis, *Ptaxis, *LambdaInvMassaxis, *XiInvMassaxis;
  AliEventPoolManager *fPoolMgr; // event pool manager
  
  TH1D *hTrackPt;
  TH1D *hTrackEta;
  TH1D *hTrackPhi;
  TH1D *hEventCent;
  TH1D *hEventZvtx;

  TH2D *hTrackdEdx;
  TH2D *hNSigmaPion;
  TH2D *hNSigmaKaon;
  TH2D *hNSigmaProton;

  TH3D *hInvMassLambda, *hInvMassXi;
  TH3D *****hSameLambda_SGNL, *****hSameLambda_SB;
  TH3D *****hSameXi;
  TH3D *****hMixLambda_SGNL, *****hMixLambda_SB;
  TH3D *****hMixXi;

  Float_t fCentV0M;
  
  Bool_t fIsAOD;
  Int_t Nmaxmixevents;
  Int_t Nmaxmixtracks;
 
  ClassDef(AliAnalysisTaskThreePartCorr, 1);
};
#endif
