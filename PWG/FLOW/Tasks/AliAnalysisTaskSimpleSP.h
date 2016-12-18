/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice */
/* $Id: $ */

#ifndef AliAnalysisTaskSimpleSP_H
#define AliAnalysisTaskSimpleSP_H

/////////////////////////////////////////////////
// AliAnalysisTaskSimpleSP:
// analysis task for Scalar Product method
// Author: Naomi van der Kolk (kolk@nikhef.nl)
/////////////////////////////////////////////////

class AliFlowEventSimple;
class AliFlowAnalysisWithSimpleSP;
class TList;

#include "TString.h"
#include "AliAnalysisTaskSE.h"

//===============================================================

class AliAnalysisTaskSimpleSP : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskSimpleSP();
  AliAnalysisTaskSimpleSP(const char *name, Bool_t usePhiWeights);
  virtual ~AliAnalysisTaskSimpleSP();
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);

  void   SetUsePhiWeights(Bool_t const aPhiW) {this->fUsePhiWeights = aPhiW;}
  Bool_t GetUsePhiWeights() const             {return this->fUsePhiWeights;}

  void     SetRelDiffMsub(Double_t diff) { this->fRelDiffMsub = diff; }
  Double_t GetRelDiffMsub() const        { return this->fRelDiffMsub; }
  
  void SetApplyCorrectionForNUA(Bool_t const applyCorrectionForNUA) {this->fApplyCorrectionForNUA = applyCorrectionForNUA;};
  Bool_t GetApplyCorrectionForNUA() const {return this->fApplyCorrectionForNUA;};
  
  void SetHarmonic(Int_t const harmonic) {this->fHarmonic = harmonic;};
  void SetUseWeights(Bool_t weights) {fWeights = weights; }
  void SetUseScaling(Bool_t scaling) {fScaling = scaling; }
  Int_t GetHarmonic() const {return this->fHarmonic;};   

  void SetBehaveAsEP() { fNormalizationType = 0; }
  void SetV0SanityCheck(Bool_t t){ fV0SanityCheck = t;}
  
  void SetTotalQvector(const char *tqv) {*this->fTotalQvector = tqv;}; 

  void SetBookOnlyBasicCCH(Bool_t const aMB) {this->fMinimalBook = aMB;}

 private:

  AliAnalysisTaskSimpleSP(const AliAnalysisTaskSimpleSP& aAnalysisTask);
  AliAnalysisTaskSimpleSP& operator=(const AliAnalysisTaskSimpleSP& aAnalysisTask); 

  AliFlowEventSimple*               fEvent;         //input event
  AliFlowAnalysisWithSimpleSP* fSP;            // analysis object
  TList*                            fListHistos;    // collection of output

  Bool_t    fMinimalBook;   // flag to turn off QA and minimize FlowCommonHist                                                                                                                                                
  Bool_t    fUsePhiWeights; // use phi weights
  TList*    fListWeights;   // list with weights

  Double_t  fRelDiffMsub;   // the relative difference the two subevent multiplicities can have
  
  Bool_t fApplyCorrectionForNUA; // apply automatic correction for non-uniform acceptance 
  
  Int_t fHarmonic;               // harmonic
  Bool_t fWeights;               // use event weights
  Bool_t fScaling;               // use q-vector scaling
  Int_t fNormalizationType;      // 0: EP mode || 1: SP mode (default)
  Bool_t fV0SanityCheck;          // 14102016 test flag

  TString   *fTotalQvector;      // total Q-vector is: "QaQb" (means Qa+Qb), "Qa"  or "Qb"  
  
  ClassDef(AliAnalysisTaskSimpleSP, 1); // example of analysis
};

//==================================================================

#endif
