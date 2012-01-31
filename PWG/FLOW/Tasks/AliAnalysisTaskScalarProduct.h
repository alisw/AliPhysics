/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice */
/* $Id: $ */

#ifndef AliAnalysisTaskScalarProduct_H
#define AliAnalysisTaskScalarProduct_H

/////////////////////////////////////////////////
// AliAnalysisTaskScalarProduct:
// analysis task for Scalar Product method
// Author: Naomi van der Kolk (kolk@nikhef.nl)
/////////////////////////////////////////////////

class AliFlowEventSimple;
class AliFlowAnalysisWithScalarProduct;
class TList;

#include "TString.h"
#include "AliAnalysisTaskSE.h"

//===============================================================

class AliAnalysisTaskScalarProduct : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskScalarProduct();
  AliAnalysisTaskScalarProduct(const char *name, Bool_t usePhiWeights);
  virtual ~AliAnalysisTaskScalarProduct();
  
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
  Int_t GetHarmonic() const {return this->fHarmonic;};   

  void SetBehaveAsEP() { fNormalizationType = 0; }
  
  void SetTotalQvector(const char *tqv) {*this->fTotalQvector = tqv;}; 

 private:

  AliAnalysisTaskScalarProduct(const AliAnalysisTaskScalarProduct& aAnalysisTask);
  AliAnalysisTaskScalarProduct& operator=(const AliAnalysisTaskScalarProduct& aAnalysisTask); 

  AliFlowEventSimple*               fEvent;         //input event
  AliFlowAnalysisWithScalarProduct* fSP;            // analysis object
  TList*                            fListHistos;    // collection of output

  Bool_t    fUsePhiWeights; // use phi weights
  TList*    fListWeights;   // list with weights

  Double_t  fRelDiffMsub;   // the relative difference the two subevent multiplicities can have
  
  Bool_t fApplyCorrectionForNUA; // apply automatic correction for non-uniform acceptance 
  
  Int_t fHarmonic;               // harmonic
  Int_t fNormalizationType;      // 0: EP mode || 1: SP mode (default)

  TString   *fTotalQvector;      // total Q-vector is: "QaQb" (means Qa+Qb), "Qa"  or "Qb"  
  
  ClassDef(AliAnalysisTaskScalarProduct, 1); // example of analysis
};

//==================================================================

#endif
