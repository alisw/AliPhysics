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

// --- Standard library ---

// --- AliRoot header files ---
#include "AliQACheckerBase.h"

class AliHMPIDQAChecker: public AliQACheckerBase {

public:
  AliHMPIDQAChecker() : AliQACheckerBase("HMPID","HMPID Quality Assurance Data Checker") {;}          // ctor
  AliHMPIDQAChecker(const AliHMPIDQAChecker& qac) : AliQACheckerBase(qac.GetName(), qac.GetTitle()) {;} // cpy ctor   
  AliHMPIDQAChecker& operator = (const AliHMPIDQAChecker& qac) ;
  virtual ~AliHMPIDQAChecker() {;} // dtor

private:
  
  ClassDef(AliHMPIDQAChecker,1)  // description 

};

#endif // AliHMPIDQAChecker_H
