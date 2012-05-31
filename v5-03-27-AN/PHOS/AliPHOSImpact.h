#ifndef ALIPHOSIMPACT_H
#define ALIPHOSIMPACT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/* History of cvs commits:
 *
 * $Log$
 * Revision 1.2  2005/05/28 14:19:04  schutz
 * Compilation warnings fixed by T.P.
 *
 */

//_________________________________________________________________________
//  Hit impact class for PHOS     
//  A hit impact in PHOS is a set of parameters of a track which
//  enters the detector for the first time.
//  Track parameters are:
//  - track number
//  - primary particle number
//  - type of a particle
//  - impact coordinate
//  - impact 4-momentum
//               
//*-- Author: Yuri Kharlov (IHEP, Protvino/SUBATECH, Nantes)

// --- ROOT system ---
#include <TLorentzVector.h>

// --- AliRoot header files ---
#include "AliHit.h"

// --- Standard library ---

class AliPHOSImpact : public AliHit {

public:
  
  AliPHOSImpact();              // default ctor 
  AliPHOSImpact(const AliPHOSImpact & hit) ; 
  AliPHOSImpact(Int_t shunt, Int_t primary, Int_t track, Int_t pid, TLorentzVector p, Float_t *xyz);
  virtual ~AliPHOSImpact(void) { } // dtor 
  
  Int_t   GetPid(void)         const { 
    // returns the particle PDG code which initiates this hit
    return fPid ; 
  }
  Int_t   GetPrimary(void)     const { 
    // returns the primary particle id at the origin of this hit 
    return fPrimary ; 
  }
  TLorentzVector GetMomentum() const {
    // returns momentum of the particle which initiated this hit
    return  fMomentum;
  }
  void Print(const Option_t * = "")const;

private:
  AliPHOSImpact & operator = (const AliPHOSImpact & /*impact*/);

private:

  Int_t          fPid ;       // type of the particle that initiates that hit 
  Int_t          fPrimary ;   // Primary particles at the origine of the hit
  TLorentzVector fMomentum;   // 4-momentum of the particle

  ClassDef(AliPHOSImpact,1)  // Hit impact for PHOS

} ;

//////////////////////////////////////////////////////////////////////////////

#endif // ALIPHOSIMPACT_H
