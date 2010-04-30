#ifndef ALIACORDEQACHECKER_H
#define ALIACORDEQACHECKER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//
//  Checks the quality assurance for ACORDE. 
//  Default implementation from Yves skeleton
//
//  Authors:
//      Mario Rodriguez Cahuantzi <mrodrigu@mail.cern.ch> (FCFM-BUAP)
//      Luciano Diaz Gonzalez <luciano.diaz@nucleares.unam.mx> (ICN-UNAM)
//      Arturo Fernandez <afernan@mail.cern.ch> (FCFM-BUAP)
//  Last update: Nov. 14t 2009 --> MRC <mrodrigu@mail.cern.ch> (FCFM-BUAP) 
//...


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
  virtual ~AliACORDEQAChecker() {;} // destructor
  virtual void Check(Double_t *, AliQAv1::ALITASK_t /*index*/) ;
  virtual void Check(Double_t *, AliQAv1::ALITASK_t /*index*/, TObjArray ** list, const AliDetectorRecoParam * /* recoParam*/) ;

  Double_t CheckAcordeRefHits(const TH1 * href, const TH1 * hdata) const;

private:

  
  ClassDef(AliACORDEQAChecker,1)  // description 

};

#endif // AliACORDEQAChecker_H
