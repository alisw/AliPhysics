#ifndef ALIVZEROQACHECKER_H
#define ALIVZEROQACHECKER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


/*
  Checks the quality assurance. 
  By comparing with reference data
*/


// --- ROOT system ---
class TFile ; 
class TH1F ; 
class TH1I ; 

// --- Standard library ---

// --- AliRoot header files ---
#include "AliQACheckerBase.h"

class AliVZEROLoader ; 

class AliVZEROQAChecker: public AliQACheckerBase {

public:
  AliVZEROQAChecker() : AliQACheckerBase("VZERO","VZERO Quality Assurance Data Maker") {;}          // ctor
  AliVZEROQAChecker(const AliVZEROQAChecker& qac) : AliQACheckerBase(qac.GetName(), qac.GetTitle()) {;} // cpy ctor   
  AliVZEROQAChecker& operator = (const AliVZEROQAChecker& qac) ;
  virtual ~AliVZEROQAChecker() {;} // dtor

private:
  
  ClassDef(AliVZEROQAChecker,1)  // description 

};

#endif // AliVZEROQAChecker_H
