#ifndef ALIANALYSISTASKUNICORME_H
#define ALIANALYSISTASKUNICORME_H

/* Copyright(c) 1998-2048, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice */
/* $Id$ */

// Author: Dariusz Miskowiec 2007

//=============================================================================
// unicor analysis task
//=============================================================================
#include "AliAnalysisTaskME.h"
class AliUnicorEventAliceESD;

//=============================================================================
class AliAnalysisTaskUnicorME : public AliAnalysisTaskME {
   
 public:                                        
  AliAnalysisTaskUnicorME(const char *name="dali");           // default constructor
  AliAnalysisTaskUnicorME(const AliAnalysisTaskUnicorME &ta): AliAnalysisTaskME(ta), fEv0(ta.fEv0), fEv1(ta.fEv1), fOutputList(ta.fOutputList) {}
  AliAnalysisTaskUnicorME &operator=(const AliAnalysisTaskUnicorME &) {return (*this);} 
  virtual ~AliAnalysisTaskUnicorME() {}                       // destructor
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);

 protected:
  AliUnicorEventAliceESD *fEv0;                       //! data/analysis interface
  AliUnicorEventAliceESD *fEv1;                       //! data/analysis interface
  TList          *fOutputList;                //! list of analysis objects

  ClassDef(AliAnalysisTaskUnicorME,1)
};
//=============================================================================

#endif
