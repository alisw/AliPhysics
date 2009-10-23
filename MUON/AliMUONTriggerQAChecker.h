#ifndef ALIMUONTRIGGERQACHECKER_H
#define ALIMUONTRIGGERQACHECKER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id: AliMUONTriggerQAChecker.h 34140 2009-08-06 13:02:16Z hristov $

/// \ingroup rec 
/// \class AliMUONTriggerQAChecker
/// \brief Implementation of QAChecker for MTR
///
//  Author: Laurent Aphecetche

#include "AliMUONVQAChecker.h"

class TObjArray;
class TH1;

class AliMUONTriggerQAChecker: public AliMUONVQAChecker {
public:
  AliMUONTriggerQAChecker();
  virtual ~AliMUONTriggerQAChecker();

  virtual ECheckCode * CheckRaws(TObjArray** list, AliMUONRecoParam* recoParam);
  virtual ECheckCode * CheckRecPoints(TObjArray** list, AliMUONRecoParam* recoParam);
  virtual ECheckCode * CheckESD(TObjArray** list, AliMUONRecoParam* recoParam);
  
private:

  AliMUONVQAChecker::ECheckCode MarkHisto(TH1& histo, AliMUONVQAChecker::ECheckCode value) const;

  ClassDef(AliMUONTriggerQAChecker,1)  // MUON quality assurance checker

};
#endif 
