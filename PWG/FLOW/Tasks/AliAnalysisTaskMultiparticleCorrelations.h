/* 
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. 
 * See cxx source for full Copyright notice 
 * $Id$ 
 */

/****************************************
 * analysis task for flow analysis with *
 *     multi-particle correlations      * 
 *                                      * 
 * author: Ante Bilandzic               *
 *         (abilandzic@gmail.com)       * 
 ***************************************/

#ifndef ALIANALYSISTASKMULTIPARTICLECORRELATIONS_H
#define ALIANALYSISTASKMULTIPARTICLECORRELATIONS_H

#include "AliAnalysisTaskSE.h"
#include "AliFlowAnalysisWithMultiparticleCorrelations.h"
#include "AliFlowEventSimple.h"

//================================================================================================================

class AliAnalysisTaskMultiparticleCorrelations : public AliAnalysisTaskSE{
 public:
  AliAnalysisTaskMultiparticleCorrelations();
  AliAnalysisTaskMultiparticleCorrelations(const char *name, Bool_t useParticleWeights=kFALSE);
  virtual ~AliAnalysisTaskMultiparticleCorrelations(){}; 
  
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *);
  
  // Internal flags:
  void SetMinNoRPs(Int_t const min) {fUseInternalFlags = kTRUE; this->fMinNoRPs = min;};
  Int_t GetMinNoRPs() const {return this->fMinNoRPs;};
  void SetMaxNoRPs(Int_t const max) {fUseInternalFlags = kTRUE; this->fMaxNoRPs = max;};
  Int_t GetMaxNoRPs() const {return this->fMaxNoRPs;};
  void SetExactNoRPs(Int_t const exact) {fUseInternalFlags = kTRUE; this->fExactNoRPs = exact;};
  Int_t GetExactNoRPs() const {return this->fExactNoRPs;};

  // Control histograms:
  void SetFillControlHistograms(Bool_t const fch) {this->fFillControlHistograms = fch;};
  Bool_t GetFillControlHistograms() const {return this->fFillControlHistograms;};
  void SetFillKinematicsHist(Bool_t const fkh) {this->fFillKinematicsHist = fkh;};
  Bool_t GetFillKinematicsHist() const {return this->fFillKinematicsHist;};
  void SetFillMultDistributionsHist(Bool_t const mdh) {this->fFillMultDistributionsHist = mdh;};
  Bool_t GetFillMultDistributionsHist() const {return this->fFillMultDistributionsHist;};
  void SetFillMultCorrelationsHist(Bool_t const mch) {this->fFillMultCorrelationsHist = mch;};
  Bool_t GetFillMultCorrelationsHist() const {return this->fFillMultCorrelationsHist;};

  // Q-vector:
  void SetCalculateQvector(Bool_t const cqv) {this->fCalculateQvector = cqv;};
  Bool_t GetCalculateQvector() const {return this->fCalculateQvector;};

  // Weights:               // TBI KG
  void SetPhiWeightsHist(TH1D* const phwh) {phwh->SetDirectory(0);this->fPhiWeightsHist = (TH1D*)phwh->Clone();}; // TBI
  TH1D* GetPhiWeightsHist() const {return this->fPhiWeightsHist;};
  void SetPtWeightsHist(TH1D* const ptwh) {ptwh->SetDirectory(0);this->fPtWeightsHist = (TH1D*)ptwh->Clone();}; // TBI
  TH1D* GetPtWeightsHist() const {return this->fPtWeightsHist;};
  void SetEtaWeightsHist(TH1D* const ewh) {ewh->SetDirectory(0);this->fEtaWeightsHist = (TH1D*)ewh->Clone();}; // TBI
  TH1D* GetEtaWeightsHist() const {return this->fEtaWeightsHist;};

  // Correlations:
  void SetCalculateCorrelations(Bool_t const cc) {this->fCalculateCorrelations = cc;};
  Bool_t GetCalculateCorrelations() const {return this->fCalculateCorrelations;};
  void SetCalculateIsotropic(Bool_t const ci) {this->fCalculateIsotropic = ci;};
  Bool_t GetCalculateIsotropic() const {return this->fCalculateIsotropic;};
  void SetCalculateSame(Bool_t const csh) {this->fCalculateSame = csh;};
  Bool_t GetCalculateSame() const {return this->fCalculateSame;};
  void SetSkipZeroHarmonics(Bool_t const szh) {this->fSkipZeroHarmonics = szh;};
  Bool_t GetSkipZeroHarmonics() const {return this->fSkipZeroHarmonics;};
  void SetCalculateSameIsotropic(Bool_t const csi) {this->fCalculateSameIsotropic = csi;};
  Bool_t GetCalculateSameIsotropic() const {return this->fCalculateSameIsotropic;};
  void SetCalculateAll(Bool_t const ca) {this->fCalculateAll = ca;};
  Bool_t GetCalculateAll() const {return this->fCalculateAll;};
  void SetDontGoBeyond(Int_t const dgb) {this->fDontGoBeyond = dgb;};
  Int_t GetDontGoBeyond() const {return this->fDontGoBeyond;};

  // Cumulants:
  void SetCalculateCumulants(Bool_t const cc) {this->fCalculateCumulants = cc;};
  Bool_t GetCalculateCumulants() const {return this->fCalculateCumulants;};

  // Nested loops:
  void SetCrossCheckWithNestedLoops(Bool_t const ccwnl) {this->fCrossCheckWithNestedLoops = ccwnl;};
  Bool_t GetCrossCheckWithNestedLoops() const {return this->fCrossCheckWithNestedLoops;};
  
  // 'Standard candles':
  void SetCalculateStandardCandles(Bool_t const csc) {this->fCalculateStandardCandles = csc;};
  Bool_t GetCalculateStandardCandles() const {return this->fCalculateStandardCandles;};

 private:
  AliAnalysisTaskMultiparticleCorrelations(const AliAnalysisTaskMultiparticleCorrelations& aatqc);
  AliAnalysisTaskMultiparticleCorrelations& operator=(const AliAnalysisTaskMultiparticleCorrelations& aatqc);
  
  AliFlowEventSimple *fEvent; // the input event
  AliFlowAnalysisWithMultiparticleCorrelations *fMPC; // "multi-particle correlations" object
  TList *fHistList; // base list to hold all output object (a.k.a. grandmother of all lists)

  // Internal flags:
  Bool_t fUseInternalFlags; // use internal flags (automatically set if some internal flag is used)
  Int_t fMinNoRPs; // minimum number of RPs required for the analysis 
  Int_t fMaxNoRPs; // maximum number of RPs allowed for the analysis 
  Int_t fExactNoRPs; // exact (randomly shuffled) number of RPs selected for the analysis 

  // Control histograms:
  Bool_t fFillControlHistograms;     // fill or not control histograms (by default they are filled)
  Bool_t fFillKinematicsHist;        // fill or not fKinematicsHist[2][3]
  Bool_t fFillMultDistributionsHist; // fill or not TH1D *fMultDistributionsHist[3]    
  Bool_t fFillMultCorrelationsHist;  // fill or not TH2D *fMultCorrelationsHist[3] 

  // Q-vector:
  Bool_t fCalculateQvector; // to calculate or not to calculate Q-vector components, that's a Boolean...

  // Weights:
  TH1D *fPhiWeightsHist; // histogram holding phi weights
  TH1D *fPtWeightsHist;  // histogram holding pt weights
  TH1D *fEtaWeightsHist; // histogram holding eta weights 

  // Correlations:
  Bool_t fCalculateCorrelations;  // calculate and store correlations, or perhaps not, if the weather is bad...
  Bool_t fCalculateIsotropic;     // calculate only isotropic correlations
  Bool_t fCalculateSame;          // calculate only 'same abs harmonics' correlations TBI 
  Bool_t fSkipZeroHarmonics;      // skip correlations which have some of the harmonicc equal to zero
  Bool_t fCalculateSameIsotropic; // calculate all isotropic correlations in 'same abs harmonic' TBI this can be implemented better
  Bool_t fCalculateAll;           // calculate all possible correlations 
  Int_t  fDontGoBeyond;           // do not go beyond fDontGoBeyond-p correlators

  // Cumulants:
  Bool_t fCalculateCumulants; // calculate and store cumulants, or perhaps not, if the weather is bad...

  // Nested loops:
  Bool_t fCrossCheckWithNestedLoops; // cross-check results with nested loops

  // 'Standard candles':
  Bool_t fCalculateStandardCandles; // calculate and store 'standard candles'
  
  ClassDef(AliAnalysisTaskMultiparticleCorrelations,1); 

};

//================================================================================================================

#endif











