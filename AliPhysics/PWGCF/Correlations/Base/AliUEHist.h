#ifndef AliUEHist_H
#define AliUEHist_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliUEHist.h 20164 2007-08-14 15:31:50Z morsch $ */

// encapsulate histogram and corrections for one underlying event histogram

#include "TObject.h"
#include "TString.h"

class AliCFContainer;
class TH1;
class TH1F;
class TH3;
class TH3F;
class TH1D;
class TH2;
class TH2D;
class TCollection;
class AliCFGridSparse;
class THnSparse;
class THnBase;

class AliUEHist : public TObject
{
 public:
  AliUEHist(const char* reqHist = "", const char* binning = 0);
  virtual ~AliUEHist();
  
  const UInt_t fkRegions;
  enum Region { kToward = 0, kAway, kMin, kMax };
  
  static const Int_t fgkCFSteps;
  enum CFStep { kCFStepAll = 0, kCFStepTriggered, kCFStepVertex, kCFStepAnaTopology, kCFStepTrackedOnlyPrim, kCFStepTracked, kCFStepReconstructed, kCFStepRealLeading, kCFStepBiasStudy, kCFStepBiasStudy2, kCFStepCorrected };
  
  const char* GetRegionTitle(Region region);
  const char* GetStepTitle(CFStep step);
  
  AliCFContainer* GetTrackHist(Region region) { return fTrackHist[region]; }
  AliCFContainer* GetEventHist() { return fEventHist; }
  AliCFContainer* GetTrackHistEfficiency()     { return fTrackHistEfficiency; }
  TH3F* GetMCRecoPtCorrelation() { return fFakePt; } 
 
  void SetTrackHist(Region region, AliCFContainer* hist) { fTrackHist[region] = hist; }
  void SetEventHist(AliCFContainer* hist) { fEventHist = hist; }
  void SetTrackHistEfficiency(AliCFContainer* hist) { fTrackHistEfficiency = hist; }
  
  void CopyReconstructedData(AliUEHist* from);
  void DeepCopy(AliUEHist* from);
  
  TH1* GetUEHist(CFStep step, Region region, Float_t ptLeadMin = -1, Float_t ptLeadMax = -1, Int_t multBinBegin = 0, Int_t multBinEnd = -1, Int_t twoD = 0, Bool_t etaNorm = kTRUE, Long64_t* normEvents = 0);
  TH1* GetPtHist(CFStep step, Region region, Float_t ptLeadMin, Float_t ptLeadMax, Int_t multBinBegin, Int_t multBinEnd, Float_t phiMin, Float_t phiMax, Float_t etaMin, Float_t etaMax, Bool_t skipPhiNormalization = kFALSE);
  TH2* GetSumOfRatios(AliUEHist* mixed, CFStep step, Region region, Float_t ptLeadMin, Float_t ptLeadMax, Int_t multBinBegin, Int_t multBinEnd, Bool_t etaNorm = kTRUE, Bool_t useVertexBins = kFALSE);
  
  void GetHistsZVtx(AliUEHist::CFStep step, AliUEHist::Region region, Float_t ptLeadMin, Float_t ptLeadMax, Int_t multBinBegin, Int_t multBinEnd, TH3** trackHist, TH1** eventHist);
  void GetHistsZVtxMult(AliUEHist::CFStep step, AliUEHist::Region region, Float_t ptLeadMin, Float_t ptLeadMax, THnBase** trackHist, TH2** eventHist);
  
  TH2* GetSumOfRatios2(AliUEHist* mixed, AliUEHist::CFStep step, AliUEHist::Region region, Float_t ptLeadMin, Float_t ptLeadMax, Int_t multBinBegin, Int_t multBinEnd, Bool_t normalizePerTrigger = kTRUE, Int_t stepForMixed = -1, Int_t *trigger = NULL);
  
  TH1* GetTriggersAsFunctionOfMultiplicity(AliUEHist::CFStep step, Float_t ptLeadMin, Float_t ptLeadMax);

  TH1* GetTrackEfficiency(CFStep step1, CFStep step2, Int_t axis1, Int_t axis2 = -1, Int_t source = 1, Int_t axis3 = -1);
  THnBase* GetTrackEfficiencyND(CFStep step1, CFStep step2);
  TH1* GetEventEfficiency(CFStep step1, CFStep step2, Int_t axis1, Int_t axis2 = -1, Float_t ptLeadMin = -1, Float_t ptLeadMax = -1);
  TH1* GetBias(CFStep step1, CFStep step2, Int_t region, const char* axis, Float_t leadPtMin = 0, Float_t leadPtMax = -1, Int_t weighting = 0);
  
  TH1D* GetTrackingEfficiency(Int_t axis);
  TH2D* GetTrackingEfficiency();
  TH2D* GetTrackingEfficiencyCentrality();
  
  TH2D* GetFakeRate();
  TH1D* GetFakeRate(Int_t axis);

  TH1D* GetTrackingContamination(Int_t axis);
  TH2D* GetTrackingContamination();
  TH2D* GetTrackingContaminationCentrality();
  
  TH1D* GetTrackingCorrection(Int_t axis);
  TH2D* GetTrackingCorrection();
  
  TH1D* GetTrackingEfficiencyCorrection(Int_t axis);
  TH2D* GetTrackingEfficiencyCorrection();
  TH2D* GetTrackingEfficiencyCorrectionCentrality();
  
  TH2* GetCorrelatedContamination();

  void ExtendTrackingEfficiency(Bool_t verbose = kFALSE);
  
  void Correct(AliUEHist* corrections);
  void CorrectTracks(CFStep step1, CFStep step2, TH1* trackCorrection, Int_t var1, Int_t var2 = -1);
  void CorrectTracks(CFStep step1, CFStep step2, Int_t region, TH1* trackCorrection, Int_t var1, Int_t var2 = -1);
  void CorrectEvents(CFStep step1, CFStep step2, TH1* eventCorrection, Int_t var1, Int_t var2 = -1);
  void CorrectCorrelatedContamination(CFStep step1, Int_t region, TH1* trackCorrection);
  
  void CondenseBin(THnSparse* grid, THnSparse* target, Int_t axis, Float_t targetValue, Float_t from, Float_t to);
  void CondenseBin(CFStep step, Int_t trackAxis, Int_t eventAxis, Float_t targetValue, Float_t from = 0, Float_t to = -1, CFStep tmpStep = AliUEHist::kCFStepBiasStudy2);
  void SymmetrizepTBins();
  
  void SetCombineMinMax(Bool_t flag) { fCombineMinMax = flag; }
  
  void SetEtaRange(Float_t etaMin, Float_t etaMax) { fEtaMin = etaMin; fEtaMax = etaMax; }
  void SetPtRange(Float_t ptMin, Float_t ptMax)    { fPtMin = ptMin; fPtMax = ptMax; }
  void SetPartSpecies(Int_t species)    { fPartSpecies = species;}
  void SetCentralityRange(Float_t min, Float_t max)    { fCentralityMin = min; fCentralityMax = max; }
  void SetZVtxRange(Float_t min, Float_t max)          { fZVtxMin = min; fZVtxMax = max; }
  void SetPt2Min(Float_t ptMin)			       { fPt2Min = ptMin; }
  void SetPt2Max(Float_t ptMin)			       { fPt2Max = ptMin; }
  
  Float_t GetTrackEtaCut() { return fTrackEtaCut; }
  void SetTrackEtaCut(Float_t value) { fTrackEtaCut = value; }
  void SetWeightPerEvent(Bool_t flag)   { fWeightPerEvent = flag; }
  void SetSkipScaleMixedEvent(Bool_t flag)  { fSkipScaleMixedEvent = flag; }
  
  void SetContaminationEnhancement(TH1F* hist)    { fContaminationEnhancement = hist; }
  
  void SetHistogramType(const char* histogramType)  { fHistogramType = histogramType; }
  
  void CountEmptyBins(AliUEHist::CFStep step, Float_t ptLeadMin, Float_t ptLeadMax);
  
  void AdditionalDPhiCorrection(Int_t step);
  
  void SetBinLimits(AliCFGridSparse* grid);
  void SetBinLimits(THnBase* grid);

  void ResetBinLimits(AliCFGridSparse* grid);
  void ResetBinLimits(THnBase* grid);
  
  void SetGetMultCache(Bool_t flag = kTRUE) { fGetMultCacheOn = flag; }
  
  AliUEHist(const AliUEHist &c);
  AliUEHist& operator=(const AliUEHist& corr);
  virtual void Copy(TObject& c) const;

  virtual Long64_t Merge(TCollection* list);
  void Scale(Double_t factor);
  void Reset();
  THnBase* ChangeToThn(THnBase* sparse);
  
  static TString CombineBinning(TString defaultBinning, TString customBinning);
  
protected:
  Double_t* GetBinning(const char* configuration, const char* tag, Int_t& nBins);
  void SetStepNames(AliCFContainer* container);
  void WeightHistogram(TH3* hist1, TH1* hist2);
  void MultiplyHistograms(THnSparse* grid, THnSparse* target, TH1* histogram, Int_t var1, Int_t var2);

  AliCFContainer* fTrackHist[4];      // container for track level distributions in four regions (toward, away, min, max) and at all analysis steps
  AliCFContainer* fEventHist;         // container for event level distribution at all analysis steps
  AliCFContainer* fTrackHistEfficiency; // container for tracking efficiency and contamination (all particles filled including leading one): axes: eta, pT, particle species
  TH3F* fFakePt;
 
  Float_t fEtaMin;                    // eta min for projections
  Float_t fEtaMax;                    // eta max for projections
  Float_t fPtMin;                     // pT min for projections (for track pT, not pT,lead)
  Float_t fPtMax;                     // pT max for projections (for track pT, not pT,lead)
  Int_t fPartSpecies;                   // Particle species for projections 
  Float_t fCentralityMin;             // centrality min for projections
  Float_t fCentralityMax;             // centrality max for projections
  Float_t fZVtxMin;                   // z vtx min for projections
  Float_t fZVtxMax;                   // z vtx max for projections
  Float_t fPt2Min;		      // pT min for projections (for pT,2 (only 2+1 corr case))
  Float_t fPt2Max;		      // pT max for projections (for pT,2 (only 2+1 corr case))
  
  TH1F* fContaminationEnhancement;    // histogram that contains the underestimation of secondaries in the MC as function of pT
  
  Bool_t fCombineMinMax;              // flag to combine min and max to a general towards region
  Float_t fTrackEtaCut;               // cut used during production of histograms (needed for finite bin correction in GetSumOfRatios)
  Bool_t fWeightPerEvent;	      // weight with the number of trigger particles per event
  Bool_t fSkipScaleMixedEvent;        // scale the mixed event with (0, 0) plus finite bin correction (default: kTRUE)
  
  AliCFContainer* fCache;             //! cache variable for GetTrackEfficiency
  
  Bool_t fGetMultCacheOn;             //! cache for GetHistsZVtxMult function active
  THnBase* fGetMultCache;             //! cache for GetHistsZVtxMult function
  
  TString fHistogramType;             // what is stored in this histogram
  
  ClassDef(AliUEHist, 16) // underlying event histogram container
};

#endif
