/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice */
/* $Id: $ */

#ifndef AliAnalysisTaskScalarProduct_H
#define AliAnalysisTaskScalarProduct_H

// AliAnalysisTaskScalarProduct:
// analysis task for Scalar Product method
// Author: Naomi van der Kolk (kolk@nikhef.nl)

class AliFlowEventSimple;
class AliFlowAnalysisWithScalarProduct;
class TList;

#include "TString.h"
#include "AliAnalysisTask.h"

class AliAnalysisTaskScalarProduct : public AliAnalysisTask {
 public:
  AliAnalysisTaskScalarProduct();
  AliAnalysisTaskScalarProduct(const char *name);
  virtual ~AliAnalysisTaskScalarProduct();
  
  virtual void   ConnectInputData(Option_t *);
  virtual void   CreateOutputObjects();
  virtual void   Exec(Option_t *option);
  virtual void   Terminate(Option_t *);

 private:

  AliAnalysisTaskScalarProduct(const AliAnalysisTaskScalarProduct& aAnalysisTask);
  AliAnalysisTaskScalarProduct& operator=(const AliAnalysisTaskScalarProduct& aAnalysisTask); 

  AliFlowEventSimple*               fEvent;      //input event
  AliFlowAnalysisWithScalarProduct* fSP;         // analysis object
  TList*                            fListHistos; // collection of output

  
  ClassDef(AliAnalysisTaskScalarProduct, 0); // example of analysis
};

#endif
