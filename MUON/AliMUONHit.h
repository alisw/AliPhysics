#ifndef ALIMUONHIT_H
#define ALIMUONHIT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "AliHit.h"

class AliMUONHit : public AliHit {
 public:
    Int_t     fChamber;       // Chamber number
    Float_t   fParticle;      // Geant3 particle type
    Float_t   fTheta ;        // Incident theta angle in degrees      
    Float_t   fPhi   ;        // Incident phi angle in degrees
    Float_t   fTlength;       // Track length inside the chamber
    Float_t   fEloss;         // ionisation energy loss in gas
    Float_t   fAge;           // Particle Age
    Int_t     fPHfirst;       // first padhit
    Int_t     fPHlast;        // last padhit

    Float_t   fPTot;          // hit local momentum P
    Float_t   fCxHit;         // Px/P
    Float_t   fCyHit;         // Py/P
    Float_t   fCzHit;         // Pz/P

 public:
    AliMUONHit() {}
    AliMUONHit(Int_t fIshunt, Int_t track, Int_t *vol, Float_t *hits);
    virtual ~AliMUONHit() {}
    
    ClassDef(AliMUONHit,1)  //Hit object for set:MUON
};
#endif
