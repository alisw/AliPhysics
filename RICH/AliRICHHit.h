#ifndef ALIRICHHIT_H
#define ALIRICHHIT_H


/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */


#include "AliHit.h"

class AliRICHHit : public AliHit {
 public:
    Int_t     fChamber;       // Chamber number
    Float_t   fParticle;      // Geant3 particle type
    Float_t   fTheta ;        // Incident theta angle in degrees      
    Float_t   fPhi   ;        // Incident phi angle in degrees
    Float_t   fTlength;       // Track length inside the chamber
    Float_t   fEloss;         // ionisation energy loss in gas   
    Float_t   fPHfirst;       // first padhit
    Float_t   fPHlast;        // last padhit
    Float_t   fLoss;          // did it hit the freon?
    Float_t   fMomX;            // Local Momentum
    Float_t   fMomY;            // Local Momentum
    Float_t   fMomZ;            // Local Momentum
    Float_t   fNPads;           // Pads hit
 public:
    AliRICHHit() {}
    AliRICHHit(Int_t fIshunt, Int_t track, Int_t *vol, Float_t *hits);
    virtual ~AliRICHHit() {}
    
    ClassDef(AliRICHHit,1)  //Hits object for set:RICH
};
#endif
