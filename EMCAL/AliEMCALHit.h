#ifndef ALIEMCALHIT_H
#define ALIEMCALHIT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: */

//_________________________________________________________________________
//  Hits class for EMCAL    
//  A hit in EMCAL is the sum of all hits from the same primary
//  in the same segment of scintillator. 
//               
//*-- Author: Sahal Yacoob (LBL /UCT) 
// Based on AliPHOSHit
#include <TLorentzVector.h>
// --- AliRoot header files ---
#include "AliHit.h"

#include <iostream.h>

class AliEMCALHit : public AliHit {
    friend ostream& operator << (ostream&,AliEMCALHit&);
 public:
    AliEMCALHit(); // default ctor
    AliEMCALHit(const AliEMCALHit & hit);
    AliEMCALHit(Int_t shunt, Int_t primary, Int_t tracknumber, Int_t id,
		Float_t *hits,TLorentzVector *p);
    virtual ~AliEMCALHit(void) {}// dtor
    //returns the energy loss for this hit
    Float_t GetEnergy(void) const{return fELOS;}
    // return the identificator of this his
    Int_t   GetId(void) const { return fId;}
    // returns the primary particle id at the origine of this hit 
    Int_t   GetPrimary(void) const{return fPrimary;}
    // returns the energy/momentum LorentzVector of the enetering particle.
    TLorentzVector& GetP(void) {return fP;}
    Bool_t operator == (AliEMCALHit const &rValue) const;
    AliEMCALHit operator + (const AliEMCALHit& rValue);

 private:
    Int_t          fId;        // Absolute Id number EMCAL segment
    Float_t        fELOS;      // Energy deposited
    Int_t          fPrimary;   // Primary particles at the origine of the hit
    TLorentzVector fP;         // Primary partical enetrence momentum/energy

    ClassDef(AliEMCALHit,1)  // Hit for EMCAL

};
#endif // ALIEMCALHIT_H
