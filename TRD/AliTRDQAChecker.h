#ifndef ALITRDQACHECKER_H
#define ALITRDQACHECKER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  Checks the quality assurance by comparing with reference data         //
//                                                                        //
//  Author:                                                               //
//    Sylwester Radomski (radomski@physi.uni-heidelberg.de)               //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

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

  AliTRDQAChecker() : AliQACheckerBase("TRD","TRD Quality Assurance Data Maker") {;} 
  AliTRDQAChecker(const AliTRDQAChecker& qac) : AliQACheckerBase(qac.GetName(), qac.GetTitle()) {;} 
  AliTRDQAChecker& operator = (const AliTRDQAChecker& qac) ;
  virtual ~AliTRDQAChecker() {;} 

private:
  
  ClassDef(AliTRDQAChecker,1)  // TRD QA checker

};
#endif
