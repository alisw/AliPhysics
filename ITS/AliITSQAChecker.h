#ifndef ALIITSQACHECKER_H
#define ALIITSQACHECKER_H
/* Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


/* $Id$ */

//
//  Checks the quality assurance. 
//  By comparing with reference data
//  INFN Torino
//  W. Ferrarese Oct 2007
//


// --- ROOT system ---
class TFile ; 
class TH2F ;  

// --- AliRoot header files ---
#include "AliQACheckerBase.h"
class AliITSLoader ; 

class AliITSQAChecker: public AliQACheckerBase {

public:
  AliITSQAChecker() : AliQACheckerBase("ITS","SDD Quality Assurance Data Maker") {;}          // ctor
  AliITSQAChecker(const AliITSQAChecker& qac) : AliQACheckerBase(qac.GetName(), qac.GetTitle()) {;} // cpy ctor   
  AliITSQAChecker& operator = (const AliITSQAChecker& qac) ; //operator =
  virtual ~AliITSQAChecker() {;} // dtor

private:
  
  ClassDef(AliITSQAChecker,1)  // description 

};

#endif // AliITSQAChecker_H
