#ifndef ALIITSQASSDCHECKER_H
#define ALIITSQASSDCHECKER_H
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

class AliITSQASSDChecker: public TObject {

public:
  AliITSQASSDChecker():fSubDetOffset(0) {;}          // ctor
  AliITSQASSDChecker& operator = (const AliITSQASSDChecker& qac) ; //operator =
  virtual ~AliITSQASSDChecker() {;} // dtor
  Double_t Check(AliQAv1::ALITASK_t /*index*/, TObjArray * /*list*/);
  void SetTaskOffset(Int_t TaskOffset);


private:
  
  AliITSQASSDChecker(const AliITSQASSDChecker& /*qac*/):TObject(),fSubDetOffset(0) {;} // cpy ctor   
  Int_t fSubDetOffset;            // checking operation starting point
  ClassDef(AliITSQASSDChecker,1)  // description 

};

#endif // AliITSQASSDChecker_H
