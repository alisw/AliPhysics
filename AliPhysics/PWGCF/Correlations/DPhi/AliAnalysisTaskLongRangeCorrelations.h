// -*- C++ -*-
// $Id: AliAnalysisTaskLongRangeCorrelations.h 405 2014-03-21 11:49:16Z cmayer $
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
#include  "AliTHn.h" // cannot forward declare 

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

  void SetSelectPrimaryMCParticles(Int_t flagMC, Int_t flagMCData) {
    fSelectPrimaryMCParticles     = flagMC;
    fSelectPrimaryMCDataParticles = flagMCData;
  }

  void SetRangeN(Int_t nMin, Int_t nMax, Double_t deltaEta) {
    fNMin = nMin;
    fNMax = nMax;
    fDeltaEta = deltaEta;
  }

  Double_t GetDeltaEta() const { return fDeltaEta; }
  Int_t GetNMin() const { return fNMin; }
  Int_t GetNMax() const { return fNMax; }

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

  TObjArray* GetAcceptedTracks(AliAODEvent*, TClonesArray*, Double_t);
  TObjArray* GetAcceptedTracks(TClonesArray*, Double_t);

  // filling histograms by name
  void Fill(const char*, Double_t);                           // TH1 weight=1
  void Fill(const char*, Double_t, Double_t );                // TH2 weight=1
  void Fill(const char*, Double_t, Double_t, Double_t);       // TH3 weight=1
  void Fill(const char*, const Double_t*, Double_t weight=1); // THnSparse

  void SetupForMixing();

  THnSparse* MakeHistSparsePhiEta(const char* name) const;
  AliTHn* MakeHistPhiEta(const char* name) const;
  AliTHn* MakeHistPhiEtaPhiEta(const char* name) const; 

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
  Int_t                fSelectPrimaryMCParticles;     // 0: no, 1: yes, -1: only non-primary particles
  Int_t                fSelectPrimaryMCDataParticles; // 0: no, 1: yes, -1: only non-primary particles
  Int_t                fNMin;
  Int_t                fNMax;
  Double_t             fDeltaEta;
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

