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
#include "AliQA.h"
#include "AliQACheckerBase.h"
#include "AliITSQAChecker.h"
class AliITSLoader ; 

class AliITSQASSDChecker: public TObject {

public:
  AliITSQASSDChecker() {;}          // ctor
  AliITSQASSDChecker& operator = (const AliITSQASSDChecker& qac) ; //operator =
  virtual ~AliITSQASSDChecker() {;} // dtor
  const Double_t Check(AliQA::ALITASK_t /*index*/, TObjArray * /*list*/, Int_t SubDetOffset);

private:
  
  AliITSQASSDChecker(const AliITSQASSDChecker& /*qac*/):TObject(){;} // cpy ctor   
  ClassDef(AliITSQASSDChecker,1)  // description 

};

#endif // AliITSQASSDChecker_H
