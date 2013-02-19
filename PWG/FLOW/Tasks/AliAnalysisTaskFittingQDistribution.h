/* 
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. 
 * See cxx source for full Copyright notice 
 * $Id$ 
 */

/**************************************
 *    analysis task for fitting       * 
 *         q-distribution             *
 *                                    * 
 * authors: Naomi van der Kolk        *
 *           (kolk@nikhef.nl)         *  
 *          Raimond Snellings         *
 *           (snelling@nikhef.nl)     * 
 *          Ante Bilandzic            *
 *           (abilandzic@gmail.com)   * 
 * ***********************************/

#ifndef ALIANALYSISTASKFITTINGQDISTRIBUTION_H
#define ALIANALYSISTASKFITTINGQDISTRIBUTION_H

#include "AliAnalysisTaskSE.h"
#include "AliFlowCommonConstants.h"

class TString;
class TList;
class AliFlowEventSimple;
class AliFlowAnalysisWithFittingQDistribution;
class TH2D;

//================================================================================================================

class AliAnalysisTaskFittingQDistribution : public AliAnalysisTaskSE{
 public:
  AliAnalysisTaskFittingQDistribution();
  AliAnalysisTaskFittingQDistribution(const char *name, Bool_t useWeights=kFALSE);
  virtual ~AliAnalysisTaskFittingQDistribution(){}; 

  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);

  void SetBookOnlyBasicCCH(Bool_t const bobcch) {this->fBookOnlyBasicCCH = bobcch;};
  Bool_t GetBookOnlyBasicCCH() const {return this->fBookOnlyBasicCCH;};   
  void SetUsePhiWeights(Bool_t const uPhiW) {this->fUsePhiWeights = uPhiW;};
  Bool_t GetUsePhiWeights() const {return this->fUsePhiWeights;};
  void SetHarmonic(Int_t const harmonic) {this->fHarmonic = harmonic;};
  Int_t GetHarmonic() const {return this->fHarmonic;};  
  // q-distribution:
  void SetqMin(Double_t const qmin) {this->fqMin = qmin;};
  Double_t GetqMin() const {return this->fqMin;};
  void SetqMax(Double_t const qmax) {this->fqMax = qmax;};
  Double_t GetqMax() const {return this->fqMax;};
  void SetqNbins(Int_t const qNbins) {this->fqNbins = qNbins;};
  Int_t GetqNbins() const {return this->fqNbins;};  
  void SetStoreqDistributionVsMult(Bool_t const sqdvm) {this->fStoreqDistributionVsMult = sqdvm;};
  Bool_t GetStoreqDistributionVsMult() const {return this->fStoreqDistributionVsMult;};  
  void SetqDistributionVsMult(TH2D* const qdvm) {this->fqDistributionVsMult = qdvm;};
  TH2D* GetqDistributionVsMult() const {return this->fqDistributionVsMult;};
  void SetMinMult(Double_t const minm) {this->fMinMult = minm;};
  Double_t GetMinMult() const {return this->fMinMult;};
  void SetMaxMult(Double_t const maxm) {this->fMaxMult = maxm;};
  Double_t GetMaxMult() const {return this->fMaxMult;};
  void SetnBinsMult(Int_t const nbm) {this->fnBinsMult = nbm;};
  Int_t GetnBinsMult() const {return this->fnBinsMult;};  
  void SetDoFit(Bool_t df) {this->fDoFit = df;};
  Bool_t GetDoFit() const {return this->fDoFit;}; 
  void SetMultiplicityIs(AliFlowCommonConstants::ERefMultSource mi) {this->fMultiplicityIs = mi;};
 
 private:
  AliAnalysisTaskFittingQDistribution(const AliAnalysisTaskFittingQDistribution& aatfqd);
  AliAnalysisTaskFittingQDistribution& operator=(const AliAnalysisTaskFittingQDistribution& aatfqd);

  AliFlowEventSimple* fEvent;                    // the input event
  AliFlowAnalysisWithFittingQDistribution* fFQD; // Fitting q-distribution (FQD) object
  TList *fListHistos;                            // collection of output 
     
  Bool_t fBookOnlyBasicCCH;              // book only basis common control histrograms (by default book them all) 
  Bool_t fUseWeights;                    // use any weights
  Bool_t fUsePhiWeights;                 // phi weights
  TList* fListWeights;                   // list with weights  
  Int_t fHarmonic;                       // harmonic   
  Double_t fqMin;                        // lower boundary of TH1D *fqDistribution
  Double_t fqMax;                        // upper boundary of TH1D *fqDistribution
  Int_t fqNbins;                         // number of bins of TH1D *fqDistribution                                             
  Bool_t fStoreqDistributionVsMult;      // store q-distributions vs M
  TH2D *fqDistributionVsMult;            // distribution of Q/sqrt{M} vs multiplicity
  Double_t fMinMult;                     // minimum multiplicity
  Double_t fMaxMult;                     // maximum multiplicity
  Int_t fnBinsMult;                      // number of multiplicity bins
  Bool_t fDoFit;                         // do the final fit
  AliFlowCommonConstants::ERefMultSource fMultiplicityIs;              // by default kRP
                                                           
  ClassDef(AliAnalysisTaskFittingQDistribution, 1); 
};

//================================================================================================================

#endif











