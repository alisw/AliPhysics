/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* Add a description of your MPI analysis */

#ifndef AliAnalysisTaskFlatenicity_H
#define AliAnalysisTaskFlatenicity_H

class AliESDtrackCuts;
class AliESDEvent;
class TList;
class TH1D;
class TH2D;
class TH1I;
class TProfile;
class THnSparse;

#include "AliAnalysisTaskSE.h"

#include "AliGenEventHeader.h"
#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliMCParticle.h"
#include "AliStack.h"
#include "TParticle.h"

class AliAnalysisTaskFlatenicity : public AliAnalysisTaskSE {
public:
  AliAnalysisTaskFlatenicity();
  AliAnalysisTaskFlatenicity(const char *name);
  virtual ~AliAnalysisTaskFlatenicity();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);
  Double_t GetFlatenicity();
  Double_t GetFlatenicityMC();
  void MakeMCanalysis();
  void MakeDataanalysis();

  void SetPtMin(Double_t val) {
    fPtMin = val;
  } // Set pT cut for associated particles
  void SetUseMC(Bool_t flat_flag = kFALSE) {
    fUseMC = flat_flag;
  } // use to analyse MC data
  void SetMCclosureTest(Bool_t flat_flag = kFALSE) { fIsMCclosure = flat_flag; }
  void SetRemoveTrivialScaling(Bool_t flat_flag = kFALSE) {
    fRemoveTrivialScaling = flat_flag;
  }
  bool HasRecVertex();

protected:
private:
  AliESDEvent *fESD; //! input ESD event
  AliEventCuts fEventCuts;
  AliStack *fMCStack; //! MC stack
  AliMCEvent *fMC;    //! MC Event
  Bool_t fUseMC;      // analyze MC events
  Int_t fV0Mindex;
  Bool_t fIsMCclosure;
  Bool_t fRemoveTrivialScaling;
  Int_t fnGen;
  AliPIDResponse *fPIDResponse;
  AliAnalysisFilter *fTrackFilter;
  TList *fOutputList; //! output list in the root file
  Double_t fEtaCut;
  Double_t fPtMin;
  Double_t ftrackmult08;
  Double_t fv0mpercentile;
  Float_t fFlat;
  Float_t fFlatMC;
  AliMultSelection *fMultSelection;
  TH1D *hPtPrimIn;
  TH1D *hPtPrimOut;
  TH1D *hPtSecOut;
  TH1D *hPtOut;
  TH1D *hFlatenicity;
  TH1D *hFlatenicityMC;
  TH2D *hFlatResponse;
  TH2D *hFlatVsPt;
  TH2D *hFlatVsPtMC;
  TProfile *hActivityV0DataSect;
  TProfile *hActivityV0McSect;
  TH2D *hFlatVsNchMC;
  TH2D *hFlatVsV0M;
  TH1D *hCounter;
  TH2D *hFlatVsPtV0M[9];

  AliAnalysisTaskFlatenicity(
      const AliAnalysisTaskFlatenicity &); // not implemented
  AliAnalysisTaskFlatenicity &
  operator=(const AliAnalysisTaskFlatenicity &); // not implemented

  ClassDef(AliAnalysisTaskFlatenicity, 3);
};

#endif
