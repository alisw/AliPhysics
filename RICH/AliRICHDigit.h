#ifndef ALIRICHDIGIT_H
#define ALIRICHDIGIT_H


/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "AliDigit.h"
class AliRICHDigit : public TObject {
 public:
    Int_t     fPadX;        // Pad number along x
    Int_t     fPadY ;       // Pad number along y
    Int_t     fSignal;      // Signal amplitude
    
    
    Int_t     fTcharges[100];  // charge per track making this digit (up to 10)
    Int_t     fTracks[100];    // tracks making this digit (up to 10)
    Int_t     fPhysics;        // physics contribution to signal 
    Int_t     fHit;            // hit number - temporary solution
    
 public:
    AliRICHDigit() {}
    AliRICHDigit(Int_t *digits);
    AliRICHDigit(Int_t *tracks, Int_t *charges, Int_t *digits);
    virtual ~AliRICHDigit() {}
    ClassDef(AliRICHDigit,1)  //Digits for set:RICH
};
#endif







