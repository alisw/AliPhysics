#ifndef ALITRDQUALASSCHECKER_H
#define ALITRDQUALASSCHECKER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//////////////////////////////////////////////////////
//
//  Checks the quality assurance. 
//  By comparing with reference data
//  S. Radomski Uni-Heidelberg October 2007
//
///////////////////////////////////////////////////////


// --- ROOT system ---
class TFile ; 
class TH1F ; 
class TH1I ; 
class TList ;
class TObjArray;

// --- Standard library ---

// --- AliRoot header files ---
#include "AliQACheckerBase.h"
class AliTRDLoader ; 

class AliTRDQAChecker: public AliQACheckerBase {

public:
  AliTRDQAChecker() : AliQACheckerBase("TRD","TRD Quality Assurance Data Maker") {;}          // ctor
  AliTRDQAChecker(const AliTRDQAChecker& qac) : AliQACheckerBase(qac.GetName(), qac.GetTitle()) {;} // cpy ctor   
  AliTRDQAChecker& operator = (const AliTRDQAChecker& qac) ;
  virtual ~AliTRDQAChecker() {;} // dtor

  virtual const Double_t Check() {return 1.0;}
  virtual const Double_t Check(TList * /*list*/) {return 1.0;}
  virtual const Double_t Check(TObjArray * /*list*/) {return 1.0;}

private:
  
  ClassDef(AliTRDQAChecker,1)  // description 

};

#endif // AliTRDQAChecker_H
