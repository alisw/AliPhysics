#ifndef ALIMUONQACHECKER_H
#define ALIMUONQACHECKER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup rec 
/// \class AliMUONQAChecker
/// \brief Implementation of AliQACheckerBase for MCH and MTR
///
//  Author: Laurent Aphecetche

// --- AliRoot header files ---
#include "AliQACheckerBase.h"

class TH1;
class TObjArray;

class AliMUONQAChecker: public AliQACheckerBase {

public:
  AliMUONQAChecker();
  AliMUONQAChecker(const AliMUONQAChecker& qac);
  virtual ~AliMUONQAChecker();

  virtual void   Init(const AliQAv1::DETECTORINDEX_t det) ; 

protected:

  using AliQACheckerBase::Check;
  
  virtual Double_t * Check(AliQAv1::ALITASK_t index) ;
  virtual Double_t * Check(AliQAv1::ALITASK_t index, TObjArray ** list) ;
  virtual void SetQA(AliQAv1::ALITASK_t index, Double_t * value) const ;	
	
  Double_t * CheckRaws(TObjArray** list);
  Double_t * CheckRecPoints(TObjArray** list);
  Double_t * CheckESD(TObjArray** list);
  TH1* GetHisto(TObjArray* list, const char* hname, Int_t specie) const;
  Double_t MarkHisto(TH1& histo, Double_t value) const;
  
private:
  
  ClassDef(AliMUONQAChecker,1)  // MUON quality assurance checker

};
#endif 
