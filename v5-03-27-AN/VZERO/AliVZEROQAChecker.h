#ifndef ALIVZEROQACHECKER_H
#define ALIVZEROQACHECKER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


/*
  Checks the quality of the data
  by comparing with reference data
  which should be loaded from QA ref DB
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
  AliVZEROQAChecker();
  virtual ~AliVZEROQAChecker() {;} // destructor
  
  virtual void   Init(const AliQAv1::DETECTORINDEX_t det);

  void SetLowEventCut(Int_t nEvents) {fLowEventCut = nEvents;}
  void SetORvsANDCut(Double_t cut) {fORvsANDCut = cut;}
  void SetBGvsBBCut(Double_t cut) {fBGvsBBCut = cut;}

protected:  
  virtual void Check( Double_t * test, AliQAv1::ALITASK_t index, TObjArray ** list, const AliDetectorRecoParam * recoParam);
  Double_t CheckRaws(TObjArray * list) const ;
  Double_t CheckEsds(TObjArray * list) const;
  
  virtual void SetQA(AliQAv1::ALITASK_t index, Double_t * value) const ;
  
private:
  AliVZEROQAChecker(const AliVZEROQAChecker& qac); // cpy ctor   
  AliVZEROQAChecker &operator=(const AliVZEROQAChecker& qac); // assignment operator

  Int_t    fLowEventCut; // Minimum number of events required by the QA checker
  Double_t fORvsANDCut; // VZERO OR vs VZERO AND counters cut
  Double_t fBGvsBBCut; // VZERO beam-gas vs beam-beam counters cut
  
  ClassDef(AliVZEROQAChecker,1)  // description 

};

#endif // AliVZEROQAChecker_H
