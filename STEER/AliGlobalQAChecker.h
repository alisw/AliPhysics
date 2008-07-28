#ifndef ALIGLOBALQACHECKER_H
#define ALIGLOBALQACHECKER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


/* $Id: AliGlobalQAChecker.h 27115 2008-07-04 15:12:14Z hristov $ */

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
class AliGlobalLoader ; 

class AliGlobalQAChecker: public AliQACheckerBase {

public:
  AliGlobalQAChecker() : AliQACheckerBase("Global","Global Quality Assurance Data Maker") {;}          // ctor
  AliGlobalQAChecker(const AliGlobalQAChecker& qac) : AliQACheckerBase(qac.GetName(), qac.GetTitle()) {;} // cpy ctor   
  virtual ~AliGlobalQAChecker() {;} // dtor

private:
  
  ClassDef(AliGlobalQAChecker,1)  // description 

};

#endif // AliGLOBALQAChecker_H
