//-*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTQACHECKER_H
#define ALIHLTQACHECKER_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/** @file   AliHLTQAChecker.h
    @author Matthias Richter
    @date   2009-11-24
    @brief  HLT QA checker instance
*/

#include "AliQACheckerBase.h"
#include <TObjArray.h>

class AliHLTQAChecker: public AliQACheckerBase {
  
 public:
  AliHLTQAChecker();
  virtual ~AliHLTQAChecker();

  virtual Double_t * Check(AliQAv1::ALITASK_t, TObjArray **, const AliDetectorRecoParam * recoParam); 
  void Init(const AliQAv1::DETECTORINDEX_t det); 
  void SetQA(AliQAv1::ALITASK_t index, Double_t * value) const;

private:
  AliHLTQAChecker(const AliHLTQAChecker& src);
  AliHLTQAChecker& operator=(const AliHLTQAChecker& src);
  
  Double_t CheckRAW(Int_t specie, TObjArray* list);
  Double_t CheckREC(Int_t specie, TObjArray* list);
  Double_t CheckESD(Int_t specie, TObjArray* list);

  ClassDef(AliHLTQAChecker,1)  // HLT Quality Assurance Checker

};

#endif // ALIHLTQACHECKER_H
