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
class AliMUONRecoParam;
class AliMUONVQAChecker;

class AliMUONQAChecker: public AliQACheckerBase {

public:
  AliMUONQAChecker();
  virtual ~AliMUONQAChecker();

  virtual void Init(const AliQAv1::DETECTORINDEX_t det); 

protected:

  virtual void Check(Double_t* test, AliQAv1::ALITASK_t index, TObjArray ** list, const AliDetectorRecoParam * recoParam); 

  virtual void SetQA(AliQAv1::ALITASK_t index, Double_t * value) const;	

private:
  /// Not implemented
  AliMUONQAChecker(const AliMUONQAChecker& qac);
  /// Not implemented
  AliMUONQAChecker& operator=(const AliMUONQAChecker& qac);
  
  TObjArray* fCheckers; ///< internal checkers

  ClassDef(AliMUONQAChecker,1)  // MUON quality assurance checker

};
#endif 
