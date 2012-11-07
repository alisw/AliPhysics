// -*- C++ -*-
// $Id: AliAnalysisTaskLongRangeCorrelations.h 215 2012-10-31 16:57:09Z cmayer $
#ifndef _AliAnalysisTaskLongRangeCorrelations_H_
#define _AliAnalysisTaskLongRangeCorrelations_H_

/* Copyright(c) 2012, ALICE Experiment at CERN, All rights reserved.      *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////////////////////////////
//
// Analysis class for Long-range Correlations
//
// Look for correlations in eta (and in phi)
//
// This class needs input AODs.
// The output is a list of analysis-specific containers.
//
//    Author:
//    Christoph Mayer
// 
////////////////////////////////////////////////////////////////////////


class TList;
class TObjArray;

class AliAODEvent;
class AliAODHeader;
class AliEventPoolManager;

#include <TObject.h>
#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskLongRangeCorrelations : public AliAnalysisTaskSE { 
public: 
  AliAnalysisTaskLongRangeCorrelations(const char *name="AliAnalysisTaskLongRangeCorrelations");
  virtual ~AliAnalysisTaskLongRangeCorrelations();

  /*   virtual void NotifyRun(); */
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *);
  virtual void Terminate(Option_t *);

  void SetRunMixing(Bool_t runMixing)      { fRunMixing    = runMixing; }
  void SetMixingTracks(Int_t mixingTracks) { fMixingTracks = mixingTracks; }
  void SetTrackFilter(Int_t trackFilter)   { fTrackFilter  = trackFilter; }

  void SetCentralityRange(Double_t centMin, Double_t centMax) {
    fCentMin = centMin; fCentMax = centMax;
  }
  void SetPtRange(Double_t ptMin, Double_t ptMax) {
    fPtMin = ptMin; fPtMax = ptMax;
  }
  void SetPhiRange(Double_t phiMin, Double_t phiMax) {
    fPhiMin = phiMin; fPhiMax = phiMax;
  }

  TString GetOutputListName() const;

protected:
  // for now up to second moments:
  // <n_1>(phi_1,eta_1)
  // <n_2>(phi_2,eta_2) 
  // <n_1, n_2>(phi_1,eta_1;phi_2,eta_2)
  void CalculateMoments(TObjArray*, TObjArray*, Double_t weight=1.); 
  void       ComputeN2ForThisEvent(THnSparse*, THnSparse*, const char*, Double_t weight=1);
  THnSparse* ComputeNForThisEvent(TObjArray*, const char*) const;

  TObjArray* GetAcceptedTracks(AliAODEvent* , AliAODHeader*, Double_t);

  // filling histograms by name
  void Fill(const char* histName, Double_t x);                           // TH1 weight=1
  void Fill(const char* histName, Double_t x, Double_t y);               // TH2 weight=1
  void Fill(const char* histName, const Double_t* x, Double_t weight=1); // THnSparse

  void SetupForMixing();

  THnSparse* MakeHistSparsePhiEta(const char* name) const;
  THnSparse* MakeHistSparsePhiEtaPhiEta(const char* name) const;

private:
  TList*               fOutputList;         //! Output list
  Bool_t               fRunMixing;          //
  AliEventPoolManager* fPoolMgr;            //! event pool manager  
  Int_t                fMixingTracks;       //
  Int_t                fTrackFilter;        //
  Double_t             fCentMin, fCentMax;  // centrality range
  Double_t             fPtMin, fPtMax;      // P_{T} range
  Double_t             fPhiMin, fPhiMax;    // #phi range

  // histogram data
  Int_t    fnBinsCent, fnBinsPt, fnBinsPhi, fnBinsEta;
  Double_t fxMinCent,  fxMinPt,  fxMinPhi,  fxMinEta;
  Double_t fxMaxCent,  fxMaxPt,  fxMaxPhi,  fxMaxEta;

  AliAnalysisTaskLongRangeCorrelations(const AliAnalysisTaskLongRangeCorrelations&); // not implemented
  AliAnalysisTaskLongRangeCorrelations& operator=(const AliAnalysisTaskLongRangeCorrelations&); // not implemented
  
  ClassDef(AliAnalysisTaskLongRangeCorrelations, 1);
} ; 

class LRCParticle : public TObject {
public:
  LRCParticle(Double_t eta=0, Double_t phi=0)
    : fEta(eta), fPhi(phi) {}
  virtual ~LRCParticle() {}

  Double_t Eta() const { return fEta; }
  Double_t Phi() const { return fPhi; }

protected:
private:
  LRCParticle(const LRCParticle&);
  LRCParticle& operator=(const LRCParticle&);

  Double_t fEta;
  Double_t fPhi;
  ClassDef(LRCParticle, 1);
} ;
#endif // _AliAnalysisTaskLongRangeCorrelations_H_

