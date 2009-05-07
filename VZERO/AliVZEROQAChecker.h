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
class TObjArray ; 

// --- Standard library ---

// --- AliRoot header files ---
#include "AliQACheckerBase.h"

class AliVZEROLoader ; 

class AliVZEROQAChecker: public AliQACheckerBase {

public:
  AliVZEROQAChecker() : AliQACheckerBase("VZERO","VZERO Quality Assurance Data Checker") {;}          // ctor
  AliVZEROQAChecker(const AliVZEROQAChecker& qac) : AliQACheckerBase(qac.GetName(), qac.GetTitle()) {;} // cpy ctor   
  virtual ~AliVZEROQAChecker() {;} // destructor
  
  virtual void   Init(const AliQAv1::DETECTORINDEX_t det) ; 

protected:  
  virtual  Double_t * Check(AliQAv1::ALITASK_t index, TObjArray ** list);
  virtual  Double_t * Check(AliQAv1::ALITASK_t ) ; 
  Double_t CheckEntries(TObjArray * list) const ;
  Double_t CheckEsds(TObjArray * list) const;
  
  virtual void SetQA(AliQAv1::ALITASK_t index, Double_t * value) const ;
  
private:
  
  ClassDef(AliVZEROQAChecker,1)  // description 

};

#endif // AliVZEROQAChecker_H
