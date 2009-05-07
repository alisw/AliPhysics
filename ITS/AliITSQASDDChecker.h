#ifndef ALIITSQASDDCHECKER_H
#define ALIITSQASDDCHECKER_H
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

class AliITSQASDDChecker: public TObject {

public:
  AliITSQASDDChecker():fSubDetOffset(0) {;}          // ctor
  AliITSQASDDChecker& operator = (const AliITSQASDDChecker& qac) ; //operator =
  virtual ~AliITSQASDDChecker() {;} // dtor
  Double_t Check(AliQAv1::ALITASK_t index, TObjArray * list);
  void SetTaskOffset(Int_t TaskOffset);

private:
  AliITSQASDDChecker(const AliITSQASDDChecker& /*qac*/):TObject(),fSubDetOffset(0) {;} // cpy ctor   
  Int_t fSubDetOffset;            // checking operation starting point
  ClassDef(AliITSQASDDChecker,1)  // description 

};

#endif // AliITSQASDDChecker_H
