#ifndef ALIMUONQACHECKER_H
#define ALIMUONQACHECKER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup rec 
/// \class AliMUONQAChecker
/// \brief MUON quality assurance checker
///
//  Author: Christian Finck



// --- AliRoot header files ---
#include "AliQACheckerBase.h"

class AliMUONQAChecker: public AliQACheckerBase {

public:
  AliMUONQAChecker();
  AliMUONQAChecker(const AliMUONQAChecker& qac);
  AliMUONQAChecker& operator=(const AliMUONQAChecker& qac);
  virtual ~AliMUONQAChecker();

private:
  
  ClassDef(AliMUONQAChecker,1)  // MUON quality assurance checker

};
#endif 
