#ifndef ALIMUONDIGIT_H
#define ALIMUONDIGIT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include <TObject.h>

class AliMUONDigit : public TObject {

 public:
    AliMUONDigit() {}
    AliMUONDigit(Int_t *digits);
    AliMUONDigit(Int_t *tracks, Int_t *charges, Int_t *digits);
    virtual ~AliMUONDigit();
    
    virtual Int_t    PadX()               {return fPadX;}
    virtual Int_t    PadY()               {return fPadY;}
    virtual Int_t    Signal()             {return fSignal;}
    virtual Int_t    Physics()            {return fPhysics;}
    virtual Int_t    Hit()                {return fHit;}    
    virtual Int_t    Cathode()            {return fCathode;}
    virtual Int_t    Track(Int_t i)       {return fTracks[i];}
    virtual Int_t    TrackCharge(Int_t i) {return fTcharges[i];}    
    virtual void     AddSignal(Int_t q)   {fSignal += q;}
    virtual void     AddPhysicsSignal(Int_t q)   {fPhysics += q;}	    
 private:
    Int_t     fPadX;          // Pad number along x
    Int_t     fPadY;          // Pad number along y
    Int_t     fCathode;       // Cathode number
    
    Int_t     fSignal;        // Signal amplitude
    Int_t     fTcharges[10];  // charge per track making this digit (up to 10)
    Int_t     fTracks[10];    // primary tracks making this digit (up to 10)
    Int_t     fPhysics;       // physics contribution to signal 
    Int_t     fHit;           // hit number - temporary solution

    ClassDef(AliMUONDigit,1)  //Digits for MUON
};
#endif
