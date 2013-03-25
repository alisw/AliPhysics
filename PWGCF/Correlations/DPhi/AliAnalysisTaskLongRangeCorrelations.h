// -*- C++ -*-
// $Id: AliAnalysisTaskLongRangeCorrelations.h 233 2012-12-02 14:46:41Z cmayer $
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


class TAxis;
class TClonesArray;
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
  void SetMaxAbsVertexZ(Double_t maxAbsVertexZ) { fMaxAbsVertexZ = maxAbsVertexZ; }

  TString GetOutputListName() const;

protected:
  // for now up to second moments:
  // <n_1>(phi_1,eta_1)
  // <n_2>(phi_2,eta_2) 
  // <n_1, n_2>(phi_1,eta_1;phi_2,eta_2)
  void       CalculateMoments(TString, TObjArray*, TObjArray*, Double_t, Double_t); 
  void       ComputeNXForThisEvent(TObjArray*, const char*, Double_t, Double_t);
  THnSparse* ComputeNForThisEvent(TObjArray*, const char*, Double_t) const;
  void       FillNEtaHist(TString, THnSparse*, Double_t);

  TObjArray* GetAcceptedTracks(AliAODEvent*, Double_t);
  TObjArray* GetAcceptedTracks(TClonesArray*, Double_t);

  // filling histograms by name
  void Fill(const char*, Double_t);                           // TH1 weight=1
  void Fill(const char*, Double_t, Double_t );                // TH2 weight=1
  void Fill(const char*, Double_t, Double_t, Double_t);       // TH3 weight=1
  void Fill(const char*, const Double_t*, Double_t weight=1); // THnSparse

  void SetupForMixing();

  THnSparse* MakeHistSparsePhiEta(const char* name) const;
  THnSparse* MakeHistSparsePhiEtaPhiEta(const char* name) const;

private:
  TList*               fOutputList;         //! Output list
  TAxis*               fVertexZ;            //! vertex z bins
  Bool_t               fRunMixing;          //
  AliEventPoolManager* fPoolMgr;            //! event pool manager  
  Int_t                fMixingTracks;       //
  Int_t                fTrackFilter;        //
  Double_t             fCentMin, fCentMax;  // centrality range
  Double_t             fPtMin, fPtMax;      // P_{T} range
  Double_t             fPhiMin, fPhiMax;    // #phi range
  Double_t             fMaxAbsVertexZ;      // max abs(zvertex)
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

