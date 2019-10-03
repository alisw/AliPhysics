#ifndef ALIANALYSISTASKUNICOR_H
#define ALIANALYSISTASKUNICOR_H

/* Copyright(c) 1998-2048, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice */
/* $Id$ */

// Author: Dariusz Miskowiec 2007

//=============================================================================
// unicor analysis task
//=============================================================================
#include "AliAnalysisTaskSE.h"
class AliUnicorEventAliceESD;

//=============================================================================
class AliAnalysisTaskUnicor : public AliAnalysisTaskSE {
   
 public:                                        
  AliAnalysisTaskUnicor(const char *name="dali");           // default constructor
  AliAnalysisTaskUnicor(const AliAnalysisTaskUnicor &ta): AliAnalysisTaskSE(ta), fEv0(ta.fEv0), fOutputList(ta.fOutputList) {}
  AliAnalysisTaskUnicor &operator=(const AliAnalysisTaskUnicor &) {return (*this);} 
  virtual ~AliAnalysisTaskUnicor() {}                       // destructor
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);

 protected:
  AliUnicorEventAliceESD *fEv0;                       //! data/analysis interface
  TList          *fOutputList;                //! list of analysis objects

  ClassDef(AliAnalysisTaskUnicor,1)
};
//=============================================================================

#endif
