#ifndef ALIMUONDIGIT_H
#define ALIMUONDIGIT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include <TObject.h>

class AliMUONDigit : public TObject {
 public:
    Int_t     fPadX;          // Pad number along x
    Int_t     fPadY ;         // Pad number along y
    Int_t     fSignal;        // Signal amplitude
    Int_t     fTcharges[10];  // charge per track making this digit (up to 10)
    Int_t     fTracks[10];    // primary tracks making this digit (up to 10)
    Int_t     fPhysics;       // physics contribution to signal 
    Int_t     fHit;           // hit number - temporary solution
 public:
    AliMUONDigit() {}
    AliMUONDigit(Int_t *digits);
    AliMUONDigit(Int_t *tracks, Int_t *charges, Int_t *digits);
    virtual ~AliMUONDigit();
    ClassDef(AliMUONDigit,1)  //Digits for set:MUON
};
#endif
