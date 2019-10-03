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

class TString;
class TList;
class AliFlowEventSimple;
class AliFlowAnalysisWithNestedLoops;

//================================================================================================================

class AliAnalysisTaskNestedLoops : public AliAnalysisTaskSE{
 public:
  AliAnalysisTaskNestedLoops();
  AliAnalysisTaskNestedLoops(const char *name, Bool_t useParticleWeights=kFALSE);
  virtual ~AliAnalysisTaskNestedLoops(){}; 
  
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *);
  
  // Setters and getters:
  // 0.) Common:
  void SetHarmonic(Int_t const h) {this->fHarmonic = h;};
  Int_t GetHarmonic() const {return this->fHarmonic;}; 
  void SetOppositeChargesPOI(Bool_t const ocp) {this->fOppositeChargesPOI = ocp;};
  Bool_t GetOppositeChargesPOI() const {return this->fOppositeChargesPOI;}; 
  void SetEvaluateDifferential3pCorrelator(Bool_t const ed3pc) {this->fEvaluateDifferential3pCorrelator = ed3pc;};
  Bool_t GetEvaluateDifferential3pCorrelator() const {return this->fEvaluateDifferential3pCorrelator;};       
  // 1.) Particle weights:
  void SetUsePhiWeights(Bool_t const uPhiW) {this->fUsePhiWeights = uPhiW;};
  Bool_t GetUsePhiWeights() const {return this->fUsePhiWeights;};
  void SetUsePtWeights(Bool_t const uPtW) {this->fUsePtWeights = uPtW;};
  Bool_t GetUsePtWeights() const {return this->fUsePtWeights;};
  void SetUseEtaWeights(Bool_t const uEtaW) {this->fUseEtaWeights = uEtaW;};
  Bool_t GetUseEtaWeights() const {return this->fUseEtaWeights;};
  // 2.) Nested loops for relative angle distribution (RAD): 
  void SetEvaluateNestedLoopsForRAD(Bool_t const enlfRAD) {this->fEvaluateNestedLoopsForRAD = enlfRAD;};
  Bool_t GetEvaluateNestedLoopsForRAD() const {return this->fEvaluateNestedLoopsForRAD;};
  // 3.) Debugging and cross-checking Q-cumulants:
  void SetEvaluateNestedLoopsForQC(Bool_t const enlfQC) {this->fEvaluateNestedLoopsForQC = enlfQC;};
  Bool_t GetEvaluateNestedLoopsForQC() const {return this->fEvaluateNestedLoopsForQC;};
  // 4.) Debugging and cross-checking mixed harmonics:
  void SetEvaluateNestedLoopsForMH(Bool_t const enlfMH) {this->fEvaluateNestedLoopsForMH = enlfMH;};
  Bool_t GetEvaluateNestedLoopsForMH() const {return this->fEvaluateNestedLoopsForMH;};
 
 private:
  AliAnalysisTaskNestedLoops(const AliAnalysisTaskNestedLoops& aatmh);
  AliAnalysisTaskNestedLoops& operator=(const AliAnalysisTaskNestedLoops& aatmh);
  
  AliFlowEventSimple *fEvent; // the input event
  AliFlowAnalysisWithNestedLoops *fNL; // nested loops object
  TList *fListHistos; // collection of output 
  Int_t fHarmonic; // integer n in correlators  
  Bool_t fOppositeChargesPOI; // two POIs, psi1 and psi2, in correlator <<cos[psi1+psi2-2phi3)]>> will be taken with opposite charges 
  Bool_t fEvaluateDifferential3pCorrelator; // evaluate <<cos[psi1+psi2-2phi3)]>>, where psi1 and psi2 are two POIs      
  // Particle weights:
  Bool_t fUseParticleWeights; // use any particle weights
  Bool_t fUsePhiWeights; // use phi weights
  Bool_t fUsePtWeights; // use pt weights
  Bool_t fUseEtaWeights; // use eta weights  
  TList *fWeightsList; // list with weights
  // Nested loops for relative angle distribution (RAD):
  Bool_t fEvaluateNestedLoopsForRAD; // evaluate nested loops for relative angle distribution (RAD)
  // Debugging and cross-checking Q-cumulants:
  Bool_t fEvaluateNestedLoopsForQC; // evaluate nested loops for Q-cumulants
  // Debugging and cross-checking mixed harmonics:
  Bool_t fEvaluateNestedLoopsForMH; // evaluate nested loops for mixed harmonics
  
  ClassDef(AliAnalysisTaskNestedLoops, 1); 
};

//================================================================================================================

#endif











