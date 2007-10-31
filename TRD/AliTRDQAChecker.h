#ifndef ALITRDQUALASSCHECKER_H
#define ALITRDQUALASSCHECKER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*
  Checks the quality assurance. 
  By comparing with reference data
  S. Radomski Uni-Heidelberg October 2007
*/


// --- ROOT system ---
class TFile ; 
class TH1F ; 
class TH1I ; 

// --- Standard library ---

// --- AliRoot header files ---
#include "AliQACheckerBase.h"
class AliTRDLoader ; 

class AliTRDQAChecker: public AliQACheckerBase {

public:
  AliTRDQAChecker() : AliQACheckerBase("TRD","TRD Quality Assurance Data Maker") {;}          // ctor
  AliTRDQAChecker(const AliTRDQAChecker& qac) : AliQACheckerBase(qac.GetName(), qac.GetTitle()) {;} // cpy ctor   
  AliTRDQAChecker& operator = (const AliTRDQAChecker& qac) ;
  virtual ~AliTRDQAChecker() {;} // dtor

private:
  
  ClassDef(AliTRDQAChecker,1)  // description 

};

#endif // AliTRDQAChecker_H
