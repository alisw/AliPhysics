/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef AliAnalysisTaskScalarProduct_H
#define AliAnalysisTaskScalarProduct_H

// AliAnalysisTaskScalarProduct:
// analysis task for Scalar Product method
// Author: Naomi van der Kolk (kolk@nikhef.nl)

class AliESDEvent;
class AliAODEvent;
class AliFlowAnalysisWithScalarProduct;
class AliFlowEventSimpleMaker;

#include "TString.h"
#include "AliAnalysisTask.h"

class AliAnalysisTaskScalarProduct : public AliAnalysisTask {
 public:
  AliAnalysisTaskScalarProduct(const char *name = "AliAnalysisTaskScalarProduct");
  virtual ~AliAnalysisTaskScalarProduct() {}
  
  virtual void   ConnectInputData(Option_t *);
  virtual void   CreateOutputObjects();
  virtual void   Exec(Option_t *option);
  virtual void   Terminate(Option_t *);

  void SetAnalysisType(TString type) { this->fAnalysisType = type; }
  TString GetAnalysisType() const    { return this->fAnalysisType; }
  
 private:
  AliESDEvent *fESD;                      // ESD object
  AliAODEvent *fAOD;                      // AOD object
  AliFlowAnalysisWithScalarProduct* fSP;  // analysis object
  AliFlowEventSimpleMaker* fEventMaker;   // FlowEventSimple maker object
  TString fAnalysisType;                  // can be MC, ESD or AOD

  ClassDef(AliAnalysisTaskScalarProduct, 1); // example of analysis
};

#endif