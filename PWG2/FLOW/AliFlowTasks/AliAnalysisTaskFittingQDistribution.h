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
 
 private:
  AliAnalysisTaskFittingQDistribution(const AliAnalysisTaskFittingQDistribution& aatfqd);
  AliAnalysisTaskFittingQDistribution& operator=(const AliAnalysisTaskFittingQDistribution& aatfqd);

  AliFlowEventSimple* fEvent;                    // the input event
  AliFlowAnalysisWithFittingQDistribution* fFQD; // Fitting q-distribution (FQD) object
  TList *fListHistos;                            // collection of output 
     
  Bool_t fUseWeights;                            // use any weights
  Bool_t fUsePhiWeights;                         // phi weights
  TList* fListWeights;                           // list with weights                                               
                                                           
  ClassDef(AliAnalysisTaskFittingQDistribution, 0); 
};

//================================================================================================================

#endif











