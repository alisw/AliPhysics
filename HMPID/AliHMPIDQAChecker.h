#ifndef ALIHMPIDQACHECKER_H
#define ALIHMPIDQACHECKER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


/* $Id$ */

//
//  Checks the quality assurance. 
//  By comparing with reference data
//  Skeleton for HMPID
//


// --- ROOT system ---
class TFile ; 
class TH1F ; 
class TObjArray ; 

// --- Standard library ---

// --- AliRoot header files ---
#include "AliQACheckerBase.h"

class AliHMPIDQAChecker: public AliQACheckerBase {

public:
  AliHMPIDQAChecker() : AliQACheckerBase("HMPID","HMPID Quality Assurance Data Checker") {;}          // ctor
  AliHMPIDQAChecker(const AliHMPIDQAChecker& qac) : AliQACheckerBase(qac.GetName(), qac.GetTitle()) {;} // cpy ctor   
  virtual ~AliHMPIDQAChecker() {;} // dtor

  virtual Double_t * Check(AliQA::ALITASK_t /*index*/) ;
  virtual Double_t * Check(AliQA::ALITASK_t index, TObjArray ** list) ;
  
  Double_t CheckEntries(TObjArray * list) const ;
  Double_t CheckRecPoints(TObjArray *listrec, TObjArray *listref) const ;

private:
  
  ClassDef(AliHMPIDQAChecker,1)  // description 

};

#endif // AliHMPIDQAChecker_H
