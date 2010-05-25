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
 *           (anteb@nikhef.nl)        * 
 * ***********************************/

#ifndef ALIANALYSISTASKFITTINGQDISTRIBUTION_H
#define ALIANALYSISTASKFITTINGQDISTRIBUTION_H

#include "TString.h"
#include "AliAnalysisTaskSE.h"

class TList;
class AliFlowEventSimple;
class AliFlowAnalysisWithFittingQDistribution;

//================================================================================================================

class AliAnalysisTaskFittingQDistribution : public AliAnalysisTaskSE{
 public:
  AliAnalysisTaskFittingQDistribution();
  AliAnalysisTaskFittingQDistribution(const char *name, Bool_t useWeights=kFALSE);
  virtual ~AliAnalysisTaskFittingQDistribution(){}; 

  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  
  void SetUsePhiWeights(Bool_t const uPhiW) {this->fUsePhiWeights = uPhiW;};
  Bool_t GetUsePhiWeights() const {return this->fUsePhiWeights;};
  // q-distribution:
  void SetqMin(Double_t const qmin) {this->fqMin = qmin;};
  Double_t GetqMin() const {return this->fqMin;};
  void SetqMax(Double_t const qmax) {this->fqMax = qmax;};
  Double_t GetqMax() const {return this->fqMax;};
  void SetqNbins(Int_t const qNbins) {this->fqNbins = qNbins;};
  Int_t GetqNbins() const {return this->fqNbins;};  
 
 private:
  AliAnalysisTaskFittingQDistribution(const AliAnalysisTaskFittingQDistribution& aatfqd);
  AliAnalysisTaskFittingQDistribution& operator=(const AliAnalysisTaskFittingQDistribution& aatfqd);

  AliFlowEventSimple* fEvent;                    // the input event
  AliFlowAnalysisWithFittingQDistribution* fFQD; // Fitting q-distribution (FQD) object
  TList *fListHistos;                            // collection of output 
     
  Bool_t fUseWeights;    // use any weights
  Bool_t fUsePhiWeights; // phi weights
  TList* fListWeights;   // list with weights  
  Double_t fqMin;        // lower boundary of TH1D *fqDistribution
  Double_t fqMax;        // upper boundary of TH1D *fqDistribution
  Int_t fqNbins;         // number of bins of TH1D *fqDistribution                                             
                                                           
  ClassDef(AliAnalysisTaskFittingQDistribution, 1); 
};

//================================================================================================================

#endif











