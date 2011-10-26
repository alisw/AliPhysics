#ifndef AliUEHistograms_H
#define AliUEHistograms_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliUEHistograms.h 20164 2007-08-14 15:31:50Z morsch $ */

// encapsulates several AliUEHist objects for a full UE analysis plus additional control histograms

#include "TNamed.h"
#include "AliUEHist.h"

class AliVParticle;

class TList;
class TSeqCollection;
class TObjArray;
class TH1F;
class TH2F;
class TH3F;

class AliUEHistograms : public TNamed
{
 public:
  AliUEHistograms(const char* name = "AliUEHistograms", const char* histograms = "123");
  virtual ~AliUEHistograms();
  
  void Fill(Int_t eventType, Float_t zVtx, AliUEHist::CFStep step, AliVParticle* leading, TList* toward, TList* away, TList* min, TList* max);
  void FillCorrelations(Double_t centrality, Float_t zVtx, AliUEHist::CFStep step, TObjArray* particles, TObjArray* mixed = 0, Float_t weight = 1, Bool_t firstTime = kTRUE);
  void Fill(AliVParticle* leadingMC, AliVParticle* leadingReco);
  void FillEvent(Int_t eventType, Int_t step);
  void FillEvent(Double_t centrality, Int_t step);
  void FillTrackingEfficiency(TObjArray* mc, TObjArray* recoPrim, TObjArray* recoAll, Int_t particleType, Double_t centrality = 0);
  
  void TwoTrackEfficiency(TObjArray* tracks, TObjArray* mixed, Float_t bSign);
  
  void CopyReconstructedData(AliUEHistograms* from);
  
  AliUEHist* GetUEHist(Int_t id);
  
  AliUEHist* GetNumberDensitypT() { return fNumberDensitypT; }
  AliUEHist* GetSumpT() { return fSumpT; }
  AliUEHist* GetNumberDensityPhi() { return fNumberDensityPhi; }
  
  void SetNumberDensitypT(AliUEHist* obj) { fNumberDensitypT = obj; }
  void SetSumpT(AliUEHist* obj) { fSumpT = obj; }
  void SetNumberDensityPhi(AliUEHist* obj) { fNumberDensityPhi = obj; }
  
  void SetRunNumber(Long64_t runNumber) { fRunNumber = runNumber; }
  
  TH2F* GetCorrelationpT()  { return fCorrelationpT; }
  TH2F* GetCorrelationEta() { return fCorrelationEta; }
  TH2F* GetCorrelationPhi() { return fCorrelationPhi; }
  TH2F* GetCorrelationR()   { return fCorrelationR; }
  TH2F* GetCorrelationLeading2Phi() { return fCorrelationLeading2Phi; }
  TH2F* GetCorrelationMultiplicity() { return fCorrelationMultiplicity; }
  
  TH2F* GetEventCount()     { return fEventCount; }
  TH3F* GetEventCountDifferential() { return fEventCountDifferential; }
  TH1F* GetVertexContributors() { return fVertexContributors; }
  TH1F* GetCentralityDistribution() { return fCentralityDistribution; }
  Long64_t GetRunNumber() { return fRunNumber; }
  TH3F* GetTwoTrackDistance(Int_t i) { return fTwoTrackDistancePt[i]; }
  
  void Correct(AliUEHistograms* corrections);
  
  void SetEtaRange(Float_t etaMin, Float_t etaMax);
  void SetPtRange(Float_t ptMin, Float_t ptMax);
  void SetZVtxRange(Float_t min, Float_t max);
  void SetContaminationEnhancement(TH1F* hist);
  void SetCombineMinMax(Bool_t flag);
  void SetSelectCharge(Int_t selectCharge) { fSelectCharge = selectCharge; }
  
  void ExtendTrackingEfficiency(Bool_t verbose = kFALSE);
  void Reset();

  AliUEHistograms(const AliUEHistograms &c);
  AliUEHistograms& operator=(const AliUEHistograms& c);
  virtual void Copy(TObject& c) const;

  virtual Long64_t Merge(TCollection* list);
  void Scale(Double_t factor);
  
protected:
  void FillRegion(AliUEHist::Region region, Float_t zVtx, AliUEHist::CFStep step, AliVParticle* leading, TList* list, Int_t multiplicity);
  Int_t CountParticles(TList* list, Float_t ptMin);
  
  static const Int_t fgkUEHists; // number of histograms

  AliUEHist* fNumberDensitypT;   // d^2N/dphideta vs pT,lead
  AliUEHist* fSumpT;             // d^2 sum(pT)/dphideta vs pT,lead
  AliUEHist* fNumberDensityPhi;  // d^2N/dphideta vs delta phi,lead (in pT,lead bins)
  
  TH2F* fCorrelationpT;         // pT,lead: true vs reco
  TH2F* fCorrelationEta;        // #eta,lead; true vs reco
  TH2F* fCorrelationPhi;        // #phi,lead; true vs reco
  TH2F* fCorrelationR;          // R = sqrt(delta eta^2 + delta phi^2) (true vs reco) vs pT,lead,MC
  TH2F* fCorrelationLeading2Phi;// delta phi (true vs reco) vs pT,lead,MC
  TH2F* fCorrelationMultiplicity; // number of mc particls vs reco particles (for pT > 0.5 GeV/c)
  
  TH2F* fEventCount;            // event count as function of step, (for pp: event type (plus additional step -1 for all events without vertex range even in MC)) (for PbPb: centrality)
  TH3F* fEventCountDifferential;// event count as function of leading pT, step, event type
  
  TH1F* fVertexContributors;    // number of contributors to the vertex
  TH1F* fCentralityDistribution; // distribution of the variable used for centrality selection
  TH2F* fCentralityCorrelation;  // centrality vs multiplicity
  
  TH3F* fITSClusterMap;          // its cluster map vs centrality vs pT
  
  TH3F* fTwoTrackDistancePt[2];    // control histograms for two-track efficiency study: dphi*_min vs deta (0 = before cut, 1 = after cut)
  
  Int_t fSelectCharge;           // (un)like sign selection when building correlations: 0: no selection; 1: unlike sign; 2: like sign
  
  Long64_t fRunNumber;           // run number that has been processed
  
  ClassDef(AliUEHistograms, 10)  // underlying event histogram container
};

#endif
