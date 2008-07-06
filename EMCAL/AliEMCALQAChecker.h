#ifndef ALIEMCALQACHECKER_H
#define ALIEMCALQACHECKER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


/*
  Checks the quality assurance. 
  By comparing with reference data

  Based on PHOS code written by
  Y. Schutz CERN July 2007
*/


// --- ROOT system ---
class TFile ; 
class TH1F ; 
class TH1I ; 

// --- Standard library ---

// --- AliRoot header files ---
#include "AliQACheckerBase.h"
class AliEMCALLoader ; 

class AliEMCALQAChecker: public AliQACheckerBase {

public:
  AliEMCALQAChecker() : AliQACheckerBase("EMCAL","EMCAL Quality Assurance Data Maker") {;}          // ctor
  virtual ~AliEMCALQAChecker() {;} // dtor

private:
  AliEMCALQAChecker(const AliEMCALQAChecker& qac);
  AliEMCALQAChecker& operator = (const AliEMCALQAChecker& qac) ;
  
  ClassDef(AliEMCALQAChecker,1)  // description 

};

#endif // AliEMCALQAChecker_H
