#ifndef ALIPHOSV1_H
#define ALIPHOSV1_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//_________________________________________________________________________
// Implementation version v1 of PHOS Manager class 
// Layout EMC + PPSD has name GPS2  
// Layout EMC + CPV  has name IHEP
//*--                  
//*-- Author: Yves Schutz (SUBATECH)

// --- ROOT system ---
class TClonesArray ;
class TLorentzVector ;
class TFile;

// --- AliRoot header files ---
#include "AliPHOSv0.h"

class AliPHOSv1 : public AliPHOSv0 {

public:

  AliPHOSv1(void) ;
  AliPHOSv1(const char *name, const char *title="") ;
  AliPHOSv1(const AliPHOSv1 & phos) {
    // cpy ctor: no implementation yet
    // requested by the Coding Convention
    abort() ; 
  }
  virtual ~AliPHOSv1(void) ;

  virtual void   AddHit( Int_t shunt, Int_t primary, Int_t track, Int_t id, Float_t *hits) ; 
 virtual Int_t  IsVersion(void) const {
    // Gives the version number 
    return 1 ; 
  }

  virtual void   StepManager(void) ;                              
  virtual TString Version(void){ 
    // returns the version number 
    return TString("v1") ; 
  }

  AliPHOSv1 & operator = (const AliPHOSv1 & rvalue)  {
    // assignement operator requested by coding convention but not needed
    abort() ;
    return *this ; 
  }

  void       CPVDigitize (TLorentzVector p, Float_t *xy, Int_t moduleNumber, TClonesArray *digits) ;
  Float_t    CPVPadResponseFunction(Float_t qhit, Float_t zg, Float_t xg) ;
  Double_t   CPVCumulPadResponse(Double_t x, Double_t y) ;

protected:


  ClassDef(AliPHOSv1,1)  // Implementation of PHOS manager class for layout EMC+PPSD

};

#endif // AliPHOSV1_H
