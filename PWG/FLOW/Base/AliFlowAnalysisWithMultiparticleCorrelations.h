/* 
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. 
 * See cxx source for full Copyright notice 
 * $Id$ 
 */

 /************************************ 
 * flow analysis with multi-particle *
 *           correlations            * 
 *                                   * 
 * author: Ante Bilandzic            * 
 *        (abilandzic@gmail.com)     *
 ************************************/ 

#ifndef ALIFLOWANALYSISWITHMULTIPARTICLECORRELATIONS_H
#define ALIFLOWANALYSISWITHMULTIPARTICLECORRELATIONS_H

#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TFile.h"
#include "TComplex.h"
#include "TDirectoryFile.h"
#include "Riostream.h"
#include "AliFlowEventSimple.h"
#include "AliFlowTrackSimple.h"

class AliFlowAnalysisWithMultiparticleCorrelations{
 public:
  AliFlowAnalysisWithMultiparticleCorrelations();
  virtual ~AliFlowAnalysisWithMultiparticleCorrelations(); 
  // Member functions are grouped as:
  // 0.) Methods called in the constructor;
  // 1.) Method Init() and methods called in it (!);
  // 2.) Method Make() and methods called in it;
  // 3.) Method Finish() and methods called in it;
  // 4.) Method GetOutputHistograms() and methods called in it;
  // 5.) Setters and getters;
  // 6.) The rest.

  // 0.) Methods called in the constructor:
  virtual void InitializeArraysForControlHistograms(); 
  virtual void InitializeArraysForQvector(); 
  virtual void InitializeArraysForCorrelations(); 

  // 1.) Method Init() and methods called in it (!):
  virtual void Init();
   virtual void CrossCheckSettings();
   virtual void BookAndNestAllLists(); 
   virtual void BookEverythingForBase();
   virtual void BookEverythingForControlHistograms();
   virtual void BookEverythingForWeights();
   virtual void BookEverythingForCorrelations();
   virtual void BookEverythingForCumulants();
   virtual void BookEverythingForNestedLoops();
   virtual void BookEverythingForStandardCandles();
   
  // 2.) Method Make() and methods called in it:
  virtual void Make(AliFlowEventSimple *anEvent);
   virtual Bool_t CrossCheckInternalFlags(AliFlowEventSimple *anEvent);
   virtual void FillControlHistograms(AliFlowEventSimple *anEvent);
   virtual void FillQvector(AliFlowEventSimple *anEvent);
   virtual void CalculateCorrelations();
   virtual void CalculateCumulants();
   virtual void ResetQvector();
   virtual void CrossCheckWithNestedLoops(AliFlowEventSimple *anEvent);

  // 3.) Method Finish() and methods called in it:
  virtual void Finish();
   //virtual void CrossCheck... TBI 
   virtual void CalculateStandardCandles();

  // 4.) Method GetOutputHistograms() and methods called in it: 
  virtual void GetOutputHistograms(TList *histList);
   virtual void GetPointersForControlHistograms(); 
   virtual void GetPointersForCorrelations(); 
   virtual void GetPointersForStandardCandles(); 

  // 5.) Setters and getters:
  //  5.0.) Base list and internal flags:
  void SetHistList(TList* const hlist) {this->fHistList = hlist;} 
  TList* GetHistList() const {return this->fHistList;} 
  void SetInternalFlagsPro(TProfile* const ifp) {this->fInternalFlagsPro = ifp;};
  TProfile* GetInternalFlagsPro() const {return this->fInternalFlagsPro;}; 
  void SetMinNoRPs(Int_t const min) {fUseInternalFlags = kTRUE; this->fMinNoRPs = min;};
  Int_t GetMinNoRPs() const {return this->fMinNoRPs;};
  void SetMaxNoRPs(Int_t const max) {fUseInternalFlags = kTRUE; this->fMaxNoRPs = max;};
  Int_t GetMaxNoRPs() const {return this->fMaxNoRPs;};
  void SetExactNoRPs(Int_t const exact) {fUseInternalFlags = kTRUE; this->fExactNoRPs = exact;};
  Int_t GetExactNoRPs() const {return this->fExactNoRPs;};

  //  5.1.) Control histograms:  
  void SetControlHistogramsList(TList* const chl) {this->fControlHistogramsList = chl;};
  TList* GetControlHistogramsList() const {return this->fControlHistogramsList;} 
  void SetControlHistogramsFlagsPro(TProfile* const chfp) {this->fControlHistogramsFlagsPro = chfp;};
  TProfile* GetControlHistogramsFlagsPro() const {return this->fControlHistogramsFlagsPro;}; 
  void SetFillControlHistograms(Bool_t const fch) {this->fFillControlHistograms = fch;};
  Bool_t GetFillControlHistograms() const {return this->fFillControlHistograms;};
  //  5.2.) Q-vector:
  // ...

  //  5.3.) Correlations:
  void SetCorrelationsList(TList* const cl) {this->fCorrelationsList = cl;};
  TList* GetCorrelationsList() const {return this->fCorrelationsList;} 
  void SetCorrelationsFlagsPro(TProfile* const cfp) {this->fCorrelationsFlagsPro = cfp;};
  TProfile* GetCorrelationsFlagsPro() const {return this->fCorrelationsFlagsPro;}; 
  void SetCalculateCorrelations(Bool_t const cc) {this->fCalculateCorrelations = cc;};
  Bool_t GetCalculateCorrelations() const {return this->fCalculateCorrelations;};

  //  5.4.) Cumulants:
  void SetCumulantsList(TList* const cl) {this->fCumulantsList = cl;};
  TList* GetCumulantsList() const {return this->fCumulantsList;} 
  void SetCumulantsFlagsPro(TProfile* const cfp) {this->fCumulantsFlagsPro = cfp;};
  TProfile* GetCumulantsFlagsPro() const {return this->fCumulantsFlagsPro;}; 
  void SetCalculateCumulants(Bool_t const cc) {this->fCalculateCumulants = cc;};
  Bool_t GetCalculateCumulants() const {return this->fCalculateCumulants;};

  //  5.5.) Weights: 
  void SetWeightsList(TList* const wl) {this->fWeightsList = (TList*)wl->Clone();};
  TList* GetWeightsList() const {return this->fWeightsList;} 
  void SetWeightsFlagsPro(TProfile* const wfp) {this->fWeightsFlagsPro = wfp;};
  TProfile* GetWeightsFlagsPro() const {return this->fWeightsFlagsPro;}; 
  void SetPhiWeightsHist(TH1D* const phwh) 
  { 
   fUsePhiWeights = kTRUE;
   phwh->SetDirectory(0); // TBI
   this->fPhiWeightsHist = (TH1D*)phwh->Clone(); // TBI 
  } 
  TH1D* GetPhiWeightsHist() const {return this->fPhiWeightsHist;};
  void SetPtWeightsHist(TH1D* const ptwh) 
  {
   fUsePtWeights = kTRUE;
   ptwh->SetDirectory(0); // TBI
   this->fPtWeightsHist = (TH1D*)ptwh->Clone(); // TBI
  } 
  TH1D* GetPtWeightsHist() const {return this->fPtWeightsHist;};
  void SetEtaWeightsHist(TH1D* const ewh)
  {
   fUseEtaWeights = kTRUE;
   ewh->SetDirectory(0); // TBI
   this->fEtaWeightsHist = (TH1D*)ewh->Clone(); // TBI
  } 
  TH1D* GetEtaWeightsHist() const {return this->fEtaWeightsHist;};

  //  5.6.) Nested loops:
  void SetNestedLoopsList(TList* const nll) {this->fNestedLoopsList = nll;} 
  TList* GetNestedLoopsList() const {return this->fNestedLoopsList;} 
  void SetNestedLoopsFlagsPro(TProfile* const nlfp) {this->fNestedLoopsFlagsPro = nlfp;};
  TProfile* GetNestedLoopsFlagsPro() const {return this->fNestedLoopsFlagsPro;}; 
  void SetCrossCheckWithNestedLoops(Bool_t const ccwnl) {this->fCrossCheckWithNestedLoops = ccwnl;};
  Bool_t GetCrossCheckWithNestedLoops() const {return this->fCrossCheckWithNestedLoops;};
 
  // 5.7.) 'Standard candles':
  void SetStandardCandlesList(TList* const scl) {this->fStandardCandlesList = scl;} 
  TList* GetStandardCandlesList() const {return this->fStandardCandlesList;} 
  void SetStandardCandlesFlagsPro(TProfile* const scfp) {this->fStandardCandlesFlagsPro = scfp;};
  TProfile* GetStandardCandlesFlagsPro() const {return this->fStandardCandlesFlagsPro;}; 
  void SetCalculateStandardCandles(Bool_t const csc) {this->fCalculateStandardCandles = csc;};
  Bool_t GetCalculateStandardCandles() const {return this->fCalculateStandardCandles;};
  void SetStandardCandlesHist(TH1D* const sch) {this->fStandardCandlesHist = sch;};
  TH1D* GetStandardCandlesHist() const {return this->fStandardCandlesHist;}; 

  // 6.) The rest:
  virtual void WriteHistograms(TString outputFileName);
  virtual void WriteHistograms(TDirectoryFile *outputFileName);
  virtual TComplex Q(Int_t n, Int_t p);
  virtual TComplex One(Int_t n1);
  virtual TComplex Two(Int_t n1, Int_t n2);
  virtual TComplex Three(Int_t n1, Int_t n2, Int_t n3);
  virtual TComplex Four(Int_t n1, Int_t n2, Int_t n3, Int_t n4);
  virtual TComplex Five(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Int_t n5);
  virtual TComplex Six(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Int_t n5, Int_t n6);
  virtual TComplex Seven(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Int_t n5, Int_t n6, Int_t n7);
  virtual TComplex Eight(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Int_t n5, Int_t n6, Int_t n7, Int_t n8);
  virtual Double_t PhiWeight(const Double_t &dPhi);
  virtual Double_t PtWeight(const Double_t &dPt); 
  virtual Double_t EtaWeight(const Double_t &dEta); 
 
 private:
  AliFlowAnalysisWithMultiparticleCorrelations(const AliFlowAnalysisWithMultiparticleCorrelations& afawQc);
  AliFlowAnalysisWithMultiparticleCorrelations& operator=(const AliFlowAnalysisWithMultiparticleCorrelations& afawQc); 
  // Data members are grouped as:
  // 0.) Base list and internal flags;
  // 1.) Control histograms;  
  // 2.) Q-vector;
  // 3.) Correlations;
  // 4.) Cumulants;
  // 5.) Weights;
  // 6.) Nested loops;
  // 7.) 'Standard candles'.

  // 0.) Base list and internal flags:
  TList* fHistList; // base list to hold all output object (a.k.a. grandmother of all lists)
  TProfile *fInternalFlagsPro; // profile to hold all internal flags and settings
  Bool_t fUseInternalFlags; // use internal flags (automatically set if some internal flag is used)
  Int_t fMinNoRPs; // minimum number of RPs required for the analysis 
  Int_t fMaxNoRPs; // maximum number of RPs allowed for the analysis 
  Int_t fExactNoRPs; // exact (randomly shuffled) number of RPs selected for the analysis 

  // 1.) Control histograms:  
  TList *fControlHistogramsList;        // list to hold all control histograms
  TProfile *fControlHistogramsFlagsPro; // profile to hold all flags for control histograms
  Bool_t fFillControlHistograms;        // fill or not control histograms (by default they are filled)
  TH1D *fKinematicsHist[2][3];          // [RP,POI][phi,pt,eta] distributions
  TH1D *fMultDistributionsHist[3];      // multiplicity distribution [RP,POI,reference multiplicity]
  TH2D *fMultCorrelationsHist[3];       // [RP vs. POI, RP vs. refMult, POI vs. refMult]  
  
  // 2.) Q-vector:
  TComplex fQvector[49][9]; // Q-vector components [fMaxHarmonic*fMaxCorrelator+1][fMaxCorrelator+1] = [6*8+1][8+1]  

  // 3.) Correlations:
  TList *fCorrelationsList;        // list to hold all correlations objects
  TProfile *fCorrelationsFlagsPro; // profile to hold all flags for correlations
  TProfile *fCorrelationsPro[8];   // multiparticle correlations 
  Bool_t fCalculateCorrelations;   // calculate and store correlations, or perhaps not, if the weather is bad...
  Int_t fMaxHarmonic;              // 6 (not going beyond v6, if you change this value, change also fQvector[49][9]) 
  Int_t fMaxCorrelator;            // 8 (not going beyond 8-p correlations, if you change this value, change also fQvector[49][9]) 

  // 4.) Cumulants:
  TList *fCumulantsList;        // list to hold all cumulants objects
  TProfile *fCumulantsFlagsPro; // profile to hold all flags for cumulants
  TProfile *f2pCumulantsPro;    // two-particle cumulants 
  Bool_t fCalculateCumulants;   // calculate and store cumulants, or perhaps not, if the weather is bad...
  
  // 5.) Weights:
  TList *fWeightsList;        // list to hold all weights objects
  TProfile *fWeightsFlagsPro; // profile to hold all flags for weights
  Bool_t fUsePhiWeights;      // use phi weights
  Bool_t fUsePtWeights;       // use pt weights
  Bool_t fUseEtaWeights;      // use eta weights
  TH1D *fPhiWeightsHist;      // histogram holding phi weights
  TH1D *fPtWeightsHist;       // histogram holding pt weights
  TH1D *fEtaWeightsHist;      // histogram holding eta weights 
  
  // 6.) Nested loops:
  TList *fNestedLoopsList;             // list to hold all nested loops objects
  TProfile *fNestedLoopsFlagsPro;      // profile to hold all flags for nested loops
  Bool_t fCrossCheckWithNestedLoops;   // cross-check results with nested loops
  TProfile *fNestedLoopsResultsCosPro; // profile to hold nested loops results (cosine)
  TProfile *fNestedLoopsResultsSinPro; // profile to hold nested loops results (sinus)

  // 7.) 'Standard candles':
  TList *fStandardCandlesList;        // list to hold all 'standard candles' objects
  TProfile *fStandardCandlesFlagsPro; // profile to hold all flags fo 'standard candles'
  Bool_t fCalculateStandardCandles;   // calculate and store 'standard candles'
  TH1D *fStandardCandlesHist;         // histogeam to hold results for 'standard candles'

  ClassDef(AliFlowAnalysisWithMultiparticleCorrelations,0);

};

//================================================================================================================

#endif





