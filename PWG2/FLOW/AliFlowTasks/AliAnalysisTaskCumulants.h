/* 
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. 
 * See cxx source for full Copyright notice 
 * $Id$ 
 */

/**************************************
 * analysis task for cumulant method  * 
 *                                    * 
 * authors: Naomi van der Kolk        *
 *           (kolk@nikhef.nl)         *  
 *          Raimond Snellings         *
 *           (snelling@nikhef.nl)     * 
 *          Ante Bilandzic            *
 *           (anteb@nikhef.nl)        * 
 * ***********************************/

#ifndef ALIANALYSISTASKCUMULANTS_H
#define ALIANALYSISTASKCUMULANTS_H

#include "TString.h"
#include "AliAnalysisTaskSE.h"

class TList;
class AliFlowEventSimple;
class AliFlowAnalysisWithCumulants;

//================================================================================================================

class AliAnalysisTaskCumulants : public AliAnalysisTaskSE{
 public:
  AliAnalysisTaskCumulants();
  AliAnalysisTaskCumulants(const char *name, Bool_t useWeights=kFALSE);
  virtual ~AliAnalysisTaskCumulants(){}; 
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *);
  // setters and getters:
  void SetHarmonic(Int_t const harmonic) {this->fHarmonic = harmonic;};
  Int_t GetHarmonic() const {return this->fHarmonic;};   
  void SetMultiple(Int_t const multiple) {this->fMultiple = multiple;};
  Int_t GetMultiple() const {return this->fMultiple;};    
  void SetCalculateVsMultiplicity(Bool_t const ecvm) {this->fCalculateVsMultiplicity = ecvm;};
  Bool_t GetCalculateVsMultiplicity() const {return this->fCalculateVsMultiplicity;};  
  void SetnBinsMult(Int_t const nbm) {this->fnBinsMult = nbm;};
  Int_t GetnBinsMult() const {return this->fnBinsMult;};  
  void SetMinMult(Double_t const minm) {this->fMinMult = minm;};
  Double_t GetMinMult() const {return this->fMinMult;};
  void SetMaxMult(Double_t const maxm) {this->fMaxMult = maxm;};
  Double_t GetMaxMult() const {return this->fMaxMult;};
  void SetUsePhiWeights(Bool_t const uPhiW) {this->fUsePhiWeights = uPhiW;};
  Bool_t GetUsePhiWeights() const {return this->fUsePhiWeights;};
  void SetUsePtWeights(Bool_t const uPtW) {this->fUsePtWeights = uPtW;};
  Bool_t GetUsePtWeights() const {return this->fUsePtWeights;};
  void SetUseEtaWeights(Bool_t const uEtaW) {this->fUseEtaWeights = uEtaW;};
  Bool_t GetUseEtaWeights() const {return this->fUseEtaWeights;}; 
  void SetTuneParameters(Bool_t const tp) {this->fTuneParameters = tp;};
  Bool_t GetTuneParameters() const {return this->fTuneParameters;};  
  void SetTuningR0(Double_t const tr0, Int_t const r) {this->fTuningR0[r] = tr0;};
  Double_t GetTuningR0(Int_t const r) const {return this->fTuningR0[r];};

 private:
  AliAnalysisTaskCumulants(const AliAnalysisTaskCumulants& aatc);
  AliAnalysisTaskCumulants& operator=(const AliAnalysisTaskCumulants& aatc);
  AliFlowEventSimple *fEvent;         // the input event
  AliFlowAnalysisWithCumulants *fGFC; // Generating Function Cumulant object
  TList *fListHistos;                 // collection of output 
  // common:
  Int_t fHarmonic;                    // harmonic  
  Int_t fMultiple;                    // the multiple m in p=m*n, where n is harmonic (relevant for differential flow)   
  // cumulants vs multiplicity:
  Bool_t fCalculateVsMultiplicity;    // perform flow analysis independently for each multiplicity bin 
  Int_t fnBinsMult;                   // number of multiplicity bins for flow analysis versus multiplicity  
  Double_t fMinMult;                  // minimal multiplicity for flow analysis versus multiplicity  
  Double_t fMaxMult;                  // maximal multiplicity for flow analysis versus multiplicity    
  // particle weights:  
  Bool_t fUseWeights;                 // use any weights
  Bool_t fUsePhiWeights;              // use phi weights
  Bool_t fUsePtWeights;               // use pt weights
  Bool_t fUseEtaWeights;              // use eta weights  
  TList *fWeightsList;                // list with weights
  // tuning:
  Bool_t fTuneParameters;             // tune r0 and cut series at different order
  Double_t fTuningR0[10];             // different r0 values (at maximum 10 different values allowed)
           
  ClassDef(AliAnalysisTaskCumulants, 1); 
};

//================================================================================================================

#endif











