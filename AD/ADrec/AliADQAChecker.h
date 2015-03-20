#ifndef ALIADQACHECKER_H
#define ALIADQACHECKER_H
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

class AliADLoader ; 

class AliADQAChecker: public AliQACheckerBase {

public:
  AliADQAChecker();
  virtual ~AliADQAChecker() {;} // destructor
  
  virtual void   Init(const AliQAv1::DETECTORINDEX_t det);

  void SetLowEventCut(Int_t nEvents) {fLowEventCut = nEvents;}
  void SetORvsANDCut(Double_t cut) {fORvsANDCut = cut;}
  void SetBGvsBBCut(Double_t cut) {fBGvsBBCut = cut;}

protected:  
  virtual void Check( Double_t * test, AliQAv1::ALITASK_t index, TObjArray ** list, const AliDetectorRecoParam * recoParam);
  Double_t CheckRaws(TObjArray * list) const ;
  Double_t CheckEsds(TObjArray * list) const;
  
  virtual void   MakeImage( TObjArray ** list, AliQAv1::TASKINDEX_t task, AliQAv1::MODE_t mode) ;  
  virtual void SetQA(AliQAv1::ALITASK_t index, Double_t * value) const ;
  
private:
  AliADQAChecker(const AliADQAChecker& qac); // cpy ctor   
  AliADQAChecker &operator=(const AliADQAChecker& qac); // assignment operator

  Int_t    fLowEventCut; // Minimum number of events required by the QA checker
  Double_t fORvsANDCut; // AD OR vs AD AND counters cut
  Double_t fBGvsBBCut; // AD beam-gas vs beam-beam counters cut
  
  ClassDef(AliADQAChecker,1)  // description 

};

#endif // AliADQAChecker_H
