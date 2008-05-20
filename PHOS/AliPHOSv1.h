#ifndef ALIPHOSV1_H
#define ALIPHOSV1_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/* History of cvs commits:
 *
 * $Log$
 * Revision 1.40  2006/09/13 07:31:01  kharlov
 * Effective C++ corrections (T.Pocheptsov)
 *
 * Revision 1.39  2005/05/28 14:19:05  schutz
 * Compilation warnings fixed by T.P.
 *
 */

//_________________________________________________________________________
// Implementation version v1 of PHOS Manager class 
// Layout EMC + CPV  has name IHEP
//--                  
//-- Author: Yves Schutz (SUBATECH)

// --- ROOT system ---
class TClonesArray ;
class TFile;
#include <TLorentzVector.h>

// --- AliRoot header files ---
#include "AliPHOSv0.h"

class AliPHOSv1 : public AliPHOSv0 {

public:

  AliPHOSv1(void) ;
  AliPHOSv1(const char *name, const char *title="") ;
  virtual ~AliPHOSv1(void) ;

  using AliPHOSv0::AddHit;
  virtual void   AddHit( Int_t shunt, Int_t primary, Int_t id, Float_t *hits) ; 
  virtual void   FinishEvent() ;
  virtual void   FinishPrimary() ;
  virtual Int_t  IsVersion(void) const {
    // Gives the version number 
    return 1 ; 
  }

  virtual void   StepManager(void) ;                              
  virtual const TString Version(void)const { return TString("v1") ;  }

  void       CPVDigitize (TLorentzVector p, Float_t *xy, TClonesArray *digits) ;
  Float_t    CPVPadResponseFunction(Float_t qhit, Float_t zg, Float_t xg) ;
  Double_t   CPVCumulPadResponse(Double_t x, Double_t y) ;

 private:
  AliPHOSv1(AliPHOSv1 & phos);
  AliPHOSv1 & operator = (const AliPHOSv1 & /*rvalue*/);
  TClonesArray fCPVDigits; //! Array of CPV digits per current CPV hit

  ClassDef(AliPHOSv1,5)  // Implementation of PHOS manager class for layout EMC+PPSD

};

#endif // AliPHOSV1_H
