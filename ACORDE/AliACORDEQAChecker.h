#ifndef ALIACORDEQACHECKER_H
#define ALIACORDEQACHECKER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//
//  Checks the quality assurance for ACORDE. 
//  Default implementation
//
//  Authors:
//      Mario Rodriguez Cahuantzi <mrodrigu@mail.cern.ch> (FCFM-BUAP)
//      Luciano Diaz Gonzalez <luciano.diaz@nucleares.unam.mx> (ICN-UNAM)
//      Arturo Fernandez <afernan@mail.cern.ch> (FCFM-BUAP)


// --- ROOT system ---
class TFile ; 
class TH1F ; 
class TObjArray ;

// --- Standard library ---

// --- AliRoot header files ---
#include "AliQACheckerBase.h"

class AliACORDEQAChecker: public AliQACheckerBase {

public:
  AliACORDEQAChecker() : AliQACheckerBase("ACORDE","ACORDE Quality Assurance Data Checker") {;}          // constructor
//  AliACORDEQAChecker(const AliACORDEQAChecker& qac) : AliQACheckerBase(qac.GetName(), qac.GetTitle()) {;} // constructor   
  virtual ~AliACORDEQAChecker() {;} // destructor

  virtual Double_t * Check(AliQAv1::ALITASK_t index) ;
  virtual Double_t * Check(AliQAv1::ALITASK_t index, TObjArray ** list) ;

  Double_t CheckAcordeRefHits(TObjArray *AcordeList, TObjArray *AcordeRef) const;

private:

  
  ClassDef(AliACORDEQAChecker,1)  // description 

};

#endif // AliACORDEQAChecker_H
