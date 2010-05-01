#ifndef ALIT0QACHECKER_H
#define ALIT0QACHECKER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


/* $Id$ */

//
//  Checks the quality assurance. 
//  By comparing with reference data
//  Skeleton for T0
//


// --- ROOT system ---
class TFile ; 
class TH1F ; 

// --- Standard library ---

// --- AliRoot header files ---
#include "AliQACheckerBase.h"

class AliT0QAChecker: public AliQACheckerBase {

public:
  AliT0QAChecker() : AliQACheckerBase("T0","T0 Quality Assurance Data Checker") {;}          // ctor
  AliT0QAChecker(const AliT0QAChecker& qac) : AliQACheckerBase(qac.GetName(), qac.GetTitle()) {;} // cpy ctor   
 // dtor
 virtual ~AliT0QAChecker() {;}
 Double_t CheckRaw(TObjArray *listrec , TObjArray *listref) const ;
private:
  virtual void Check(Double_t * test, AliQAv1::ALITASK_t, TObjArray ** list, const AliDetectorRecoParam * recoParam) ;
  
  ClassDef(AliT0QAChecker,1)  // description 

};

#endif // AliT0QAChecker_H
