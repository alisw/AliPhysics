#ifndef ALIEMCALHIT_H
#define ALIEMCALHIT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: */

//_________________________________________________________________________
//  Hits class for EMCAL    
//  A hit in EMCAL is the sum of all hits from the same Entering particle 
//  in the same segment of scintillator. 
//               
//*-- Author: Sahal Yacoob (LBL /UCT) 
// Based on AliPHOSHit

// --- AliRoot header files ---
#include "AliHit.h"


class AliEMCALHit : public AliHit {
    
  friend ostream& operator << (ostream&,AliEMCALHit&);
 
public:
  AliEMCALHit(); // default ctor
  AliEMCALHit(const AliEMCALHit & hit);
  AliEMCALHit(Int_t shunt, Int_t primary, Int_t tracknumber, Int_t iparent, Float_t ienergy, Int_t id, Float_t *hits,Float_t *p);
  virtual ~AliEMCALHit(void) {}// dtor
  //returns the energy loss for this hit
  Float_t GetEnergy(void) const{return fELOS;}
  // return the identificator of this his
  Int_t   GetId(void) const { return fId;}
  // returns the primary particle id at the origine of this hit 
  Int_t   GetIparent(void) const{return fIparent;}
  Float_t GetIenergy(void) const{return fIenergy;}
  Int_t   GetPrimary(void) const{return fPrimary;}
  // returns the energy/momentum LorentzVector of the enetering particle.
  Float_t GetTime(void)     const {
    // returns the time of the first energy deposition
    return fTime ;}

  Float_t GetPx(void) const{return fPx;}
  Float_t GetPy(void) const{return fPy;}
  Float_t GetPz(void) const{return fPz;}
  Float_t GetPe(void) const{return fPe;}

  void   SetIparent(Int_t i_parent) {fIparent=i_parent;}
  void   SetPrimary(Int_t primary)  {fPrimary=primary;}

  Bool_t operator == (AliEMCALHit const &rValue) const;
  AliEMCALHit operator + (const AliEMCALHit& rValue);
  
 private:

  Int_t          fId;        // Absolute Id number EMCAL segment
  Float_t        fELOS;      // Energy deposited
  Int_t          fPrimary;   // Primary particles at the origin of the hit
  Float_t        fPx;      // Primary particle entrance momentum/energy
  Float_t        fPy;      // Primary particle entrance momentum/energy
  Float_t        fPz;      // Primary particle entrance momentum/energy
  Float_t        fPe;      // Primary particle entrance momentum/energy
  Int_t          fIparent;   // Parent particle that entered emcal
  Float_t        fIenergy;   // Initial energy of parent particle that enterred the emcal
  Float_t        fTime ;      // Time of the energy deposition
  
  ClassDef(AliEMCALHit,3)  // Hit for EMCAL
    
    };

#endif // ALIEMCALHIT_H
