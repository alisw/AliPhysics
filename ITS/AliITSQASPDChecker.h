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
#include "AliQAv1.h"
#include "AliQACheckerBase.h"
#include "AliITSQAChecker.h"
class AliITSLoader ; 

class AliITSQASPDChecker: public TObject {

public:
  AliITSQASPDChecker():fSubDetOffset(0) {;}          // ctor
  AliITSQASPDChecker& operator = (const AliITSQASPDChecker& qac) ; //operator =
  virtual ~AliITSQASPDChecker() {;} // dtor
  Double_t Check(AliQAv1::ALITASK_t index, TObjArray * list);
  void SetTaskOffset(Int_t TaskOffset);
private:
  
  AliITSQASPDChecker(const AliITSQASPDChecker& /*qac*/):TObject(),fSubDetOffset(0){;}  // cpy ctor   
  Int_t fSubDetOffset;            // checking operation starting point
  ClassDef(AliITSQASPDChecker,1)  // description 

};

#endif // AliITSQASPDChecker_H
