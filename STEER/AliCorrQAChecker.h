#ifndef ALICORRQACHECKER_H
#define ALICORRQACHECKER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


/* $Id: AliCorrQAChecker.h 27115 2008-07-04 15:12:14Z hristov $ */

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
class AliCorrLoader ; 

class AliCorrQAChecker: public AliQACheckerBase {

public:
  AliCorrQAChecker() : AliQACheckerBase("Corr","Corr Quality Assurance Data Maker") {;}          // ctor
  AliCorrQAChecker(const AliCorrQAChecker& qac) : AliQACheckerBase(qac.GetName(), qac.GetTitle()) {;} // cpy ctor   
  virtual ~AliCorrQAChecker() {;} // dtor

  virtual Double_t *  Check(AliQAv1::ALITASK_t index, AliDetectorRecoParam * recoParam) { return AliQACheckerBase::Check(index, recoParam) ;}
  virtual Double_t * Check(AliQAv1::ALITASK_t index, TObjArray ** obj, AliDetectorRecoParam * recoParam) { return AliQACheckerBase::Check(index, obj, recoParam) ;}
  Double_t * Check(AliQAv1::ALITASK_t index, TNtupleD ** nData) ; 

private:
  
  ClassDef(AliCorrQAChecker,1)  // description 

};

#endif // AliCORRQAChecker_H
