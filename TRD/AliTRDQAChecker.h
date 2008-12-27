#ifndef ALITRDQACHECKER_H
#define ALITRDQACHECKER_H
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
#include "AliQA.h"
#include "AliQACheckerBase.h"
class AliTRDLoader ; 

class AliTRDQAChecker: public AliQACheckerBase {

public:
  AliTRDQAChecker() : AliQACheckerBase("TRD","TRD Quality Assurance Data Maker") {;}          // ctor
  AliTRDQAChecker(const AliTRDQAChecker& qac) : AliQACheckerBase(qac.GetName(), qac.GetTitle()) {;} // cpy ctor   
  virtual ~AliTRDQAChecker() {;} // dtor

  virtual Double_t * Check(AliQA::ALITASK_t /*index*/) {return NULL;}
  virtual Double_t * Check(TList * /*list*/) {return NULL;}
  virtual Double_t * Check(AliQA::ALITASK_t /*index*/, TObjArray ** /*list*/) {return NULL;}
  virtual Double_t * Check(AliQA::ALITASK_t /*index*/, TNtupleD** /*nt*/)     {return NULL;}

private:
  
  ClassDef(AliTRDQAChecker,1)  // description 

};

#endif // AliTRDQACHECKER_H
