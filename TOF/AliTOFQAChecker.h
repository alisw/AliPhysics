#ifndef ALITOFQACHECKER_H
#define ALITOFQACHECKER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/////////////////////////////////////////////////////////////////////
//                                                                 // 
//  Checks the quality assurance.                                  //
//  By analysis of the histograms & comparing with reference data  //
//  Author S.Arcelli                                               //
//                                                                 // 
/////////////////////////////////////////////////////////////////////

#include "AliQA.h"
#include "AliQACheckerBase.h"

//class TFile ; 
//class TH1F ; 
//class TH1I ; 

class AliTOFQAChecker: public AliQACheckerBase {

public:
  AliTOFQAChecker() : AliQACheckerBase("TOF","TOF Quality Assurance Data Maker") {;}          // ctor
  AliTOFQAChecker(const AliTOFQAChecker& qac) : AliQACheckerBase(qac.GetName(), qac.GetTitle()) {;} // cpy ctor   
  virtual ~AliTOFQAChecker() {;} // dtor

 protected:

  virtual Double_t * Check(AliQA::ALITASK_t /*index*/, TObjArray ** list) ;
  virtual Double_t * Check(AliQA::ALITASK_t /*index*/) {return NULL;} 
  
  ClassDef(AliTOFQAChecker,2)  // description 

};

#endif // AliTOFQAChecker_H
