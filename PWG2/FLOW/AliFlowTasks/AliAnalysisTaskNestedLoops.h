/* 
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. 
 * See cxx source for full Copyright notice 
 * $Id$ 
 */

/**********************************
 * analysis task for nested loops * 
 *                                * 
 * authors: Naomi van der Kolk    *
 *           (kolk@nikhef.nl)     *  
 *          Raimond Snellings     *
 *           (snelling@nikhef.nl) * 
 *          Ante Bilandzic        *
 *           (anteb@nikhef.nl)    * 
 * *******************************/

#ifndef ALIANALYSISTASKNESTEDLOOPS_H
#define ALIANALYSISTASKNESTEDLOOPS_H

#include "TString.h"
#include "AliAnalysisTaskSE.h"

class TList;
class AliFlowEventSimple;
class AliFlowAnalysisWithNestedLoops;

//================================================================================================================

class AliAnalysisTaskNestedLoops : public AliAnalysisTaskSE{
 public:
  AliAnalysisTaskNestedLoops();
  AliAnalysisTaskNestedLoops(const char *name, Bool_t useParticleWeights=kFALSE);
  virtual ~AliAnalysisTaskNestedLoops(){}; 
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  
  // particle weights:
  void SetUsePhiWeights(Bool_t const uPhiW) {this->fUsePhiWeights = uPhiW;};
  Bool_t GetUsePhiWeights() const {return this->fUsePhiWeights;};
  void SetUsePtWeights(Bool_t const uPtW) {this->fUsePtWeights = uPtW;};
  Bool_t GetUsePtWeights() const {return this->fUsePtWeights;};
  void SetUseEtaWeights(Bool_t const uEtaW) {this->fUseEtaWeights = uEtaW;};
  Bool_t GetUseEtaWeights() const {return this->fUseEtaWeights;};
 
 private:
  AliAnalysisTaskNestedLoops(const AliAnalysisTaskNestedLoops& aatmh);
  AliAnalysisTaskNestedLoops& operator=(const AliAnalysisTaskNestedLoops& aatmh);
  
  AliFlowEventSimple *fEvent; // the input event
  AliFlowAnalysisWithNestedLoops *fNL; // nested loops object
  TList *fListHistos; // collection of output 
  // particle weights:
  Bool_t fUseParticleWeights; // use any particle weights
  Bool_t fUsePhiWeights; // use phi weights
  Bool_t fUsePtWeights; // use pt weights
  Bool_t fUseEtaWeights; // use eta weights  
  TList *fWeightsList; // list with weights
  
  ClassDef(AliAnalysisTaskNestedLoops, 1); 
};

//================================================================================================================

#endif











