#ifndef ALIITSQASPDCHECKER_H
#define ALIITSQASPDCHECKER_H
/* Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


/* $Id$ */

//
//  Checks the quality assurance. 
//  By comparing with reference data
//  INFN Torino
//  P. Cerello - apr 2008
//


// --- ROOT system ---
class TFile ; 
class TH2F ;  

// --- AliRoot header files ---
#include "AliQA.h"
#include "AliQACheckerBase.h"
#include "AliITSQAChecker.h"
class AliITSLoader ; 

class AliITSQASPDChecker: public TObject {

public:
  AliITSQASPDChecker() {;}          // ctor
  AliITSQASPDChecker& operator = (const AliITSQASPDChecker& qac) ; //operator =
  virtual ~AliITSQASPDChecker() {;} // dtor
  const Double_t Check(AliQA::ALITASK_t index);
private:
  
  AliITSQASPDChecker(const AliITSQASPDChecker& qac){;}  // cpy ctor   
  ClassDef(AliITSQASPDChecker,1)  // description 

};

#endif // AliITSQASPDChecker_H
