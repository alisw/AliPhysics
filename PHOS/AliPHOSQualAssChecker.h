#ifndef ALIPHOSQUALASSCHECKER_H
#define ALIPHOSQUALASSCHECKER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


/* $Id$ */

/*
  Checks the quality assurance. 
  By comparing with reference data
  Y. Schutz CERN July 2007
*/


// --- ROOT system ---
class TFile ; 
class TH1F ; 
class TH1I ; 

// --- Standard library ---

// --- AliRoot header files ---
#include "AliQualAssCheckerBase.h"
class AliPHOSLoader ; 

class AliPHOSQualAssChecker: public AliQualAssCheckerBase {

public:
  AliPHOSQualAssChecker() : AliQualAssCheckerBase("PHOS","PHOS Quality Assurance Data Maker") {;}          // ctor
  AliPHOSQualAssChecker(const AliPHOSQualAssChecker& qac) : AliQualAssCheckerBase(qac.GetName(), qac.GetTitle()) {;} // cpy ctor   
  AliPHOSQualAssChecker& operator = (const AliPHOSQualAssChecker& qac) ;
  virtual ~AliPHOSQualAssChecker() {;} // dtor

private:
  
  ClassDef(AliPHOSQualAssChecker,1)  // description 

};

#endif // AliPHOSQualAssChecker_H
