#ifndef ALIPHOSQACHECKER_H
#define ALIPHOSQACHECKER_H
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
#include "AliQACheckerBase.h"
class AliPHOSLoader ; 

class AliPHOSQAChecker: public AliQACheckerBase {

public:
  AliPHOSQAChecker() : AliQACheckerBase("PHOS","PHOS Quality Assurance Data Maker") {;}          // ctor
  AliPHOSQAChecker(const AliPHOSQAChecker& qac) : AliQACheckerBase(qac.GetName(), qac.GetTitle()) {;} // cpy ctor   
  virtual ~AliPHOSQAChecker() {;} // dtor

private:
  AliPHOSQAChecker & operator = (const AliPHOSQAChecker & /*qac*/);
  
  ClassDef(AliPHOSQAChecker,1)  // description 

};

#endif // AliPHOSQAChecker_H
