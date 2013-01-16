/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice */
/* $Id: $ */

#ifndef AliAnalysisTaskMSP_H
#define AliAnalysisTaskMSP_H

/////////////////////////////////////////////////
// AliAnalysisTaskMSP:
// analysis task for Scalar Product method
// Author: Paul Kuijer (Paul.Kuijer@nikhef.nl)
/////////////////////////////////////////////////

class AliFlowEventSimple;
class AliFlowAnalysisWithMSP;
class TList;

#include "TString.h"
#include "AliAnalysisTaskSE.h"

//===============================================================

class AliAnalysisTaskMSP : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskMSP();
  AliAnalysisTaskMSP(const char *name, Bool_t usePhiWeights);
  virtual ~AliAnalysisTaskMSP();
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);

  void   SetUsePhiWeights(Bool_t const aPhiW) {this->fUsePhiWeights = aPhiW;}
  Bool_t GetUsePhiWeights() const             {return this->fUsePhiWeights;}
  
  void SetApplyCorrectionForNUA(Bool_t const applyCorrectionForNUA) {this->fApplyCorrectionForNUA = applyCorrectionForNUA;};
  Bool_t GetApplyCorrectionForNUA() const {return this->fApplyCorrectionForNUA;};
  
  void SetHarmonic(Int_t const harmonic) {this->fHarmonic = harmonic;};
  Int_t GetHarmonic() const {return this->fHarmonic;};   

  void SetBookOnlyBasicCCH(Bool_t const aMB) {this->fMinimalBook = aMB;}

 private:

  AliAnalysisTaskMSP(const AliAnalysisTaskMSP& aAnalysisTask);
  AliAnalysisTaskMSP& operator=(const AliAnalysisTaskMSP& aAnalysisTask); 

  AliFlowEventSimple*     fEvent;         //input event
  AliFlowAnalysisWithMSP* fSP;            // analysis object
  TList*                  fListHistos;    // collection of output

  Bool_t    fMinimalBook;   // flag to turn off QA and minimize FlowCommonHist                                                                                                                                                
  Bool_t    fUsePhiWeights; // use phi weights
  TList*    fListWeights;   // list with weights

  Bool_t fApplyCorrectionForNUA; // apply automatic correction for non-uniform acceptance 
  
  Int_t fHarmonic;               // harmonic
 
  ClassDef(AliAnalysisTaskMSP, 2); // example of analysis
};

//==================================================================

#endif
