#ifndef ALIPHOSHIT_H
#define ALIPHOSHIT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//_________________________________________________________________________
//  Hits class for PHOS     
//  A hit in PHOS is the sum of all hits in a single crystal
//               
//*-- Author: Maxime Volkov (RRC KI) & Yves Schutz (SUBATECH)

// --- ROOT system ---

// --- AliRoot header files ---
#include "AliHit.h"

// --- Standard library ---

#include <iostream.h>

class AliPHOSHit : public AliHit {

public:

  AliPHOSHit() {}
  AliPHOSHit(const AliPHOSHit & hit) ; 
  AliPHOSHit(Int_t primary, Int_t id, Float_t *hits) ;
  AliPHOSHit(Int_t primary, Int_t tracknumber, Int_t id, Float_t *hits) ;
  virtual ~AliPHOSHit(void) {}  
  
  Float_t GetEnergy(void)   const { return fELOS ; }
  Int_t   GetId(void)       const { return fId ; }
  Int_t   GetPrimary(void)  const { return fPrimary ; }

  Bool_t operator == (AliPHOSHit const &rValue) const ;
  AliPHOSHit operator + (const AliPHOSHit& rValue) const ;

  friend ostream& operator << (ostream&, const AliPHOSHit&) ;

private:

  Int_t     fId ;        // Absolute Id number of PHOS Xtal or PPSD pad
  Float_t   fELOS ;      // Energy deposited
  Int_t     fPrimary ;   // Primary particles at the origine of the hit

  ClassDef(AliPHOSHit,1)  // Hit for PHOS

} ;

//////////////////////////////////////////////////////////////////////////////

#endif // ALIPHOSHIT_H
