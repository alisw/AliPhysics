#ifndef ALIRICHCERENKOV_H
#define ALIRICHCERENKOV_H


/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */


#include "AliHit.h"
//------------------------------------------------
// Cerenkov photon  object
//------------------------------------------------

class AliRICHCerenkov: public AliHit {
 public:
    Int_t     fChamber;         // Chamber number
    Float_t   fTheta ;          // Incident theta angle in degrees      
    Float_t   fPhi   ;          // Incident phi angle in degrees
    Float_t   fTlength;         // Track length inside the chamber
    Float_t   fEloss;           // ionisation energy loss in gas
    Int_t     fPHfirst;         // first padhit
    Int_t     fPHlast;          // last padhit
    Int_t     fCMother;         // index of mother particle
    Float_t   fLoss;            // nature of particle loss
    Float_t   fIndex;           // Index of photon
    Float_t   fProduction;      // Point of production
    Float_t   fMomX;            // Local Momentum
    Float_t   fMomY;            // Local Momentum
    Float_t   fMomZ;            // Local Momentum
    Float_t   fNPads;           // Pads hit
    Float_t   fCerenkovAngle;   // Cerenkov Angle
 public:
    AliRICHCerenkov() {}
    AliRICHCerenkov(Int_t fIshunt, Int_t track, Int_t *vol, Float_t *Cerenkovs);
    virtual ~AliRICHCerenkov() {}
    
    ClassDef(AliRICHCerenkov,1)  //Cerenkovs object for set:RICH
};
#endif
