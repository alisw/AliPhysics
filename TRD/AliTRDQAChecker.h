#ifndef ALITRDQACHECKER_H
#define ALITRDQACHECKER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  Checks the quality assurance.                                         //
//  By comparing with reference data                                      //
//  S.Radomski Uni-Heidelberg October 2007                                //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

// --- ROOT system ---
class TFile ; 
class TH1F ; 
class TH1I ; 
class TList ;
class TObjArray;

// --- Standard library ---

// --- AliRoot header files ---
#include "AliQAv1.h"
#include "AliQACheckerBase.h"

class AliDetectorRecoParam;
class AliTRDLoader ; 

class AliTRDQAChecker: public AliQACheckerBase {

public:
  AliTRDQAChecker() : AliQACheckerBase("TRD","TRD Quality Assurance Data Maker") {;}          // ctor
  AliTRDQAChecker(const AliTRDQAChecker& qac) : AliQACheckerBase(qac.GetName(), qac.GetTitle()) {;} // cpy ctor   
  virtual ~AliTRDQAChecker() {;} // dtor

  virtual void Check(Double_t * test, AliQAv1::ALITASK_t /*index*/, TObjArray** /*list*/, const AliDetectorRecoParam* /*param*/) ;

private:
  
  ClassDef(AliTRDQAChecker,1)  // description 

};

#endif // AliTRDQACHECKER_H
