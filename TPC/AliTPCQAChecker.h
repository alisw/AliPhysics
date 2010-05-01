#ifndef ALITPCQACHECKER_H
#define ALITPCQACHECKER_H
/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


/* $Id: $ */

/*
  Checks implemented a la AliMUONQAChecker.
  Checks the quality assurance by very simple checks on histogram content.
  P. Christiansen, Lund, September 2009.
*/

// --- ROOT header files ---
#include <TObjArray.h>

// --- AliRoot header files ---
#include "AliQACheckerBase.h"
#include "AliDetectorRecoParam.h"

class AliTPCQAChecker: public AliQACheckerBase {
  
 public:
 AliTPCQAChecker() : AliQACheckerBase("TPC","TPC Quality Assurance Checker"), fDebug(0) {;}          // ctor
 AliTPCQAChecker(const AliTPCQAChecker& qac) : AliQACheckerBase(qac.GetName(), qac.GetTitle()), fDebug(qac.GetDebugLevel()) {;} // cpy ctor   
  virtual ~AliTPCQAChecker() {;} // dtor

  virtual void Check(Double_t *  test, AliQAv1::ALITASK_t, TObjArray **, const AliDetectorRecoParam * recoParam); 
  void Init(const AliQAv1::DETECTORINDEX_t det); 
  void SetQA(AliQAv1::ALITASK_t index, Double_t * value) const;

  Int_t GetDebugLevel() const {return fDebug;}
  void  SetDebugLevel(Int_t value) {fDebug = value;}
  
private:
  
  Double_t CheckRAW(Int_t specie, TObjArray* list);
  Double_t CheckSIM(Int_t specie, TObjArray* list);
  Double_t CheckREC(Int_t specie, TObjArray* list);
  Double_t CheckESD(Int_t specie, TObjArray* list);

  Int_t fDebug;
  
  ClassDef(AliTPCQAChecker,2)  // TPC Quality Assurance Checker

};

#endif // AliTPCQAChecker_H
