#ifndef ALIPHOSHIT_H
#define ALIPHOSHIT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////
//  Hits class for PHOS                     //
//  Version SUBATECH                          //
//  Author M. Volkov, RRC KI                  //
//  october 1999:                             // 
//            Modified by Y. Schutz SUBATECH  //
////////////////////////////////////////////////

// --- ROOT system ---

// --- AliRoot header files ---
#include "AliHit.h"
#include <iostream.h>

class AliPHOSHit : public AliHit {

protected:

  Int_t     fId ;        // Absolute Id number of PHOS Xtal or PPSD pad
  Float_t   fELOS ;      // Energy deposited
  
public:

  AliPHOSHit() {}
  AliPHOSHit(Int_t shunt, Int_t track, Int_t id, Float_t *hits) ;
  virtual ~AliPHOSHit(void) {}
  
  Float_t GetEnergy(void)   const { return fELOS ; }
  Int_t   GetId(void)       const { return fId ; }
  
  Bool_t operator == (AliPHOSHit const &rValue) const ;
  AliPHOSHit operator + (const AliPHOSHit& rValue) const ;

  friend ostream& operator << (ostream&, const AliPHOSHit&) ;

  ClassDef(AliPHOSHit,1)  // Hits object for PHOS

} ;

//////////////////////////////////////////////////////////////////////////////

#endif // ALIPHOSHIT_H
