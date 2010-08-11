#ifndef AliUEHist_H
#define AliUEHist_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliUEHist.h 20164 2007-08-14 15:31:50Z morsch $ */

// encapsulate histogram and corrections for one underlying event histogram

#include "TObject.h"

class AliCFContainer;
class TH1;
class TH3;
class TH1D;
class TH2D;
class TCollection;
class AliCFGridSparse;

class AliUEHist : public TObject
{
 public:
  AliUEHist(const char* reqHist = "");
  virtual ~AliUEHist();
  
  const Int_t fkRegions;
  enum Region { kToward = 0, kAway, kMin, kMax };
  
  static const Int_t fgkCFSteps;
  enum CFStep { kCFStepAll = 0, kCFStepTriggered, kCFStepVertex, kCFStepAnaTopology, kCFStepTrackedOnlyPrim, kCFStepTracked, kCFStepReconstructed, kCFStepRealLeading, kCFStepBiasStudy, kCFStepBiasStudy2 };
  
  const char* GetRegionTitle(Region region);
  const char* GetStepTitle(CFStep step);
  
  AliCFContainer* GetTrackHist(Region region) { return fTrackHist[region]; }
  AliCFContainer* GetEventHist() { return fEventHist; }
  
  void SetTrackHist(Region region, AliCFContainer* hist) { fTrackHist[region] = hist; }
  void SetEventHist(AliCFContainer* hist) { fEventHist = hist; }
  
  TH1D* GetUEHist(CFStep step, Region region, Float_t ptLeadMin = -1, Float_t ptLeadMax = -1);
  
  TH1* GetTrackEfficiency(CFStep step1, CFStep step2, Int_t axis1, Int_t axis2 = -1);
  TH1* GetEventEfficiency(CFStep step1, CFStep step2, Int_t axis1, Int_t axis2 = -1, Float_t ptLeadMin = -1, Float_t ptLeadMax = -1);
  TH1* GetBias(CFStep step1, CFStep step2, const char* axis, Float_t leadPtMin = 0, Float_t leadPtMax = 0);
  
  TH1D* GetTrackingEfficiency(Int_t axis);
  TH2D* GetTrackingEfficiency();
  
  TH1D* GetTrackingContamination(Int_t axis);
  TH2D* GetTrackingContamination();
  
  TH1D* GetTrackingCorrection(Int_t axis);
  TH2D* GetTrackingCorrection();
  
  void Correct(AliUEHist* corrections);
  void CorrectTracks(CFStep step1, CFStep step2, TH1* trackCorrection, Int_t var1, Int_t var2 = -1);
  void CorrectEvents(CFStep step1, CFStep step2, TH1D* eventCorrection, Int_t var);
  
  void SetCombineMinMax(Bool_t flag) { fCombineMinMax = flag; }
  
  void SetEtaRange(Float_t etaMin, Float_t etaMax) { fEtaMin = etaMin; fEtaMax = etaMax; }
  void SetPtRange(Float_t ptMin, Float_t ptMax)    { fPtMin = ptMin; fPtMax = ptMax; }
  
  AliUEHist(const AliUEHist &c);
  AliUEHist& operator=(const AliUEHist& corr);
  virtual void Copy(TObject& c) const;

  virtual Long64_t Merge(TCollection* list);
  
protected:
  void SetStepNames(AliCFContainer* container);
  void WeightHistogram(TH3* hist1, TH1* hist2);
  void SetBinLimits(AliCFGridSparse* grid);
  void ResetBinLimits(AliCFGridSparse* grid);

  AliCFContainer* fTrackHist[4];      // container for track level distributions in four regions (toward, away, min, max) and at four analysis steps
  AliCFContainer* fEventHist;         // container for event level distribution at four analysis steps
  
  Float_t fEtaMin;                    // eta min for projections
  Float_t fEtaMax;                    // eta max for projections
  Float_t fPtMin;                     // pT min for projections (for track pT, not pT,lead)
  Float_t fPtMax;                     // pT max for projections (for track pT, not pT,lead)
  
  Bool_t fCombineMinMax;              // flag to combine min and max to a general towards region
  
  AliCFContainer* fCache;             //! cache variable for GetTrackEfficiency
  
  ClassDef(AliUEHist, 1) // underlying event histogram container
};

#endif
