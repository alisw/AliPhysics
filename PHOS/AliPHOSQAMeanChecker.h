#ifndef ALIPHOSQAMEANCHECKER_H
#define ALIPHOSQAMEANCHECKER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//_________________________________________________________________________
// QA checker that compares a number with an average value plus or minus
// a width     
//                  
//*-- Author: Yves Schutz (SUBATECH)


// --- ROOT system ---

#include "TTask.h"

// --- Standard library ---

#include <assert.h>

// --- AliRoot header files ---

#include "AliPHOSQAChecker.h"  

class AliPHOSQAMeanChecker : public AliPHOSQAChecker {

 public:

  AliPHOSQAMeanChecker(){} ;          // default ctor (not to be used)
  AliPHOSQAMeanChecker(const char * name) ; // ctor
  AliPHOSQAMeanChecker(const char * name, Float_t mean, Float_t rms) ; //ctor
  virtual ~AliPHOSQAMeanChecker() ; // dtor
  virtual TString CheckingOperation() ; // where the checking operation is implemented
  //  virtual void  Exec(Option_t *option);  
  virtual void Print() ; 

  void SetMean(Float_t value) { fMean = value ; } 
  void SetRms(Float_t value) { fRms = value ; } 
  void Set(Float_t mean, Float_t rms) { fMean = mean ; fRms = rms ; } 

 private:

  Float_t fMean ; // the value that the checkable will be compared to
  Float_t fRms ;  // the range around the mean in which the test is OK
  
  ClassDef(AliPHOSQAMeanChecker,1)  // description 

};

#endif // ALIPHOSQAMeanChecker_H
