#ifndef ALIHMPIDQUALASSCHECKER_H
#define ALIHMPIDQUALASSCHECKER_H
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
#include "AliQualAssCheckerBase.h"

class AliHMPIDQualAssChecker: public AliQualAssCheckerBase {

public:
  AliHMPIDQualAssChecker() : AliQualAssCheckerBase("HMPID","HMPID Quality Assurance Data Checker") {;}          // ctor
  AliHMPIDQualAssChecker(const AliHMPIDQualAssChecker& qac) : AliQualAssCheckerBase(qac.GetName(), qac.GetTitle()) {;} // cpy ctor   
  AliHMPIDQualAssChecker& operator = (const AliHMPIDQualAssChecker& qac) ;
  virtual ~AliHMPIDQualAssChecker() {;} // dtor

private:
  
  ClassDef(AliHMPIDQualAssChecker,1)  // description 

};

#endif // AliHMPIDQualAssChecker_H
