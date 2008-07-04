#ifndef ALITPCQACHECKER_H
#define ALITPCQACHECKER_H
/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


/* $Id: $ */

/*
  Based on AliPHOSQAChecker.
  Checks the quality assurance by comparing with reference data.
  P. Christiansen, Lund, January 2008.
*/

// --- AliRoot header files ---
#include "AliQACheckerBase.h"

class AliTPCQAChecker: public AliQACheckerBase {

public:
  AliTPCQAChecker() : AliQACheckerBase("TPC","TPC Quality Assurance Checker") {;}          // ctor
  AliTPCQAChecker(const AliTPCQAChecker& qac) : AliQACheckerBase(qac.GetName(), qac.GetTitle()) {;} // cpy ctor   
  virtual ~AliTPCQAChecker() {;} // dtor

private:
  
  ClassDef(AliTPCQAChecker,1)  // TPC Quality Assurance Checker

};

#endif // AliTPCQAChecker_H
