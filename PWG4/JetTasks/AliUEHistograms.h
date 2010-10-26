#ifndef AliUEHistograms_H
#define AliUEHistograms_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliUEHistograms.h 20164 2007-08-14 15:31:50Z morsch $ */

// encapsulates several AliUEHist objects for a full UE analysis plus additional control histograms

#include "TObject.h"
#include "AliUEHist.h"

class AliVParticle;

class TList;
class TH1F;
class TH2F;
class TH3F;

class AliUEHistograms : public TObject
{
 public:
  AliUEHistograms();
  virtual ~AliUEHistograms();
  
  void Fill(Int_t eventType, AliUEHist::CFStep step, AliVParticle* leading, TList* toward, TList* away, TList* min, TList* max);
  void Fill(AliVParticle* leadingMC, AliVParticle* leadingReco);
  void FillEvent(Int_t eventType, Int_t step);
  void FillTrackingEfficiency(TObjArray* mc, TObjArray* recoPrim, TObjArray* recoAll, Int_t particleType);
  
  void CopyReconstructedData(AliUEHistograms* from);
  
  AliUEHist* GetUEHist(Int_t id);
  
  AliUEHist* GetNumberDensitypT() { return fNumberDensitypT; }
  AliUEHist* GetSumpT() { return fSumpT; }
  AliUEHist* GetNumberDensityPhi() { return fNumberDensityPhi; }
  
  TH2F* GetCorrelationpT()  { return fCorrelationpT; }
  TH2F* GetCorrelationEta() { return fCorrelationEta; }
  TH2F* GetCorrelationPhi() { return fCorrelationPhi; }
  TH2F* GetCorrelationR()   { return fCorrelationR; }
  TH2F* GetCorrelationLeading2Phi() { return fCorrelationLeading2Phi; }
  TH2F* GetCorrelationMultiplicity() { return fCorrelationMultiplicity; }
  
  TH2F* GetEventCount()     { return fEventCount; }
  TH3F* GetEventCountDifferential() { return fEventCountDifferential; }
  TH1F* GetVertexContributors() { return fVertexContributors; }
  
  void Correct(AliUEHistograms* corrections);
  
  void SetEtaRange(Float_t etaMin, Float_t etaMax);
  void SetPtRange(Float_t ptMin, Float_t ptMax);
  void SetContaminationEnhancement(TH1F* hist);
  void SetCombineMinMax(Bool_t flag);
  
  AliUEHistograms(const AliUEHistograms &c);
  AliUEHistograms& operator=(const AliUEHistograms& c);
  virtual void Copy(TObject& c) const;

  virtual Long64_t Merge(TCollection* list);
  
protected:
  void FillRegion(AliUEHist::Region region, AliUEHist::CFStep step, AliVParticle* leading, TList* list, Int_t multiplicity);
  Int_t CountParticles(TList* list, Float_t ptMin);

  AliUEHist* fNumberDensitypT;   // d^2N/dphideta vs pT,lead
  AliUEHist* fSumpT;             // d^2 sum(pT)/dphideta vs pT,lead
  AliUEHist* fNumberDensityPhi;  // d^2N/dphideta vs delta phi,lead (in pT,lead bins)
  
  TH2F* fCorrelationpT;         // pT,lead: true vs reco
  TH2F* fCorrelationEta;        // #eta,lead; true vs reco
  TH2F* fCorrelationPhi;        // #phi,lead; true vs reco
  TH2F* fCorrelationR;          // R = sqrt(delta eta^2 + delta phi^2) (true vs reco) vs pT,lead,MC
  TH2F* fCorrelationLeading2Phi;// delta phi (true vs reco) vs pT,lead,MC
  TH2F* fCorrelationMultiplicity; // number of mc particls vs reco particles (for pT > 0.5 GeV/c)
  
  TH2F* fEventCount;            // event count as function of step, event type (plus additional step -1 for all events without vertex range even in MC)
  TH3F* fEventCountDifferential;// event count as function of leading pT, step, event type
  
  TH1F* fVertexContributors;    // number of contributors to the vertex
  
  ClassDef(AliUEHistograms, 1)  // underlying event histogram container
};

#endif
