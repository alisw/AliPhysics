/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice */
/* $Id: $ */

#ifndef AliAnalysisTaskScalarProduct_H
#define AliAnalysisTaskScalarProduct_H

// AliAnalysisTaskScalarProduct:
// analysis task for Scalar Product method
// Author: Naomi van der Kolk (kolk@nikhef.nl)

class AliESDEvent;
class AliAODEvent;
class AliCFManager;
class AliFlowAnalysisWithScalarProduct;
class AliFlowEventSimpleMaker;
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

  void SetAnalysisType(TString type) { this->fAnalysisType = type; }
  TString GetAnalysisType() const    { return this->fAnalysisType; }

  void SetCFManager1(AliCFManager* cfmgr) {this->fCFManager1 = cfmgr; } 
  AliCFManager* GetCFManager1()           {return this->fCFManager1; }
  void SetCFManager2(AliCFManager* cfmgr) {this->fCFManager2 = cfmgr; } 
  AliCFManager* GetCFManager2()           {return this->fCFManager2; }

 private:

  AliAnalysisTaskScalarProduct(const AliAnalysisTaskScalarProduct& aAnalysisTask);
  AliAnalysisTaskScalarProduct& operator=(const AliAnalysisTaskScalarProduct& aAnalysisTask); 

  AliESDEvent *fESD;                      // ESD object
  AliAODEvent *fAOD;                      // AOD object
  AliFlowAnalysisWithScalarProduct* fSP;  // analysis object
  AliFlowEventSimpleMaker* fEventMaker;   // FlowEventSimple maker object
  TString fAnalysisType;                  // can be MC, ESD or AOD
  AliCFManager* fCFManager1;              // correction framework manager
  AliCFManager* fCFManager2;              // correction framework manager
  TList  *fListHistos;                    // collection of output

  ClassDef(AliAnalysisTaskScalarProduct, 1); // example of analysis
};

#endif
