#ifndef ALIRICHDIGIT_H
#define ALIRICHDIGIT_H


/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "AliDigit.h"

static const Int_t kMAXTRACKSPERRICHDIGIT = 10;

class AliRICHDigit : public TObject {
 protected:
    Int_t     fPadX;        // Pad number along x
    Int_t     fPadY ;       // Pad number along y
    Int_t     fSignal;      // Signal amplitude
    
    
    Int_t     fTcharges[kMAXTRACKSPERRICHDIGIT];  // charge per track making this digit (up to 10)
    Int_t     fTracks[kMAXTRACKSPERRICHDIGIT];    // tracks making this digit (up to 10)
    Int_t     fPhysics;        // physics contribution to signal 
    Int_t     fHit;            // hit number - temporary solution
    
 public:
    AliRICHDigit() {}
    AliRICHDigit(Int_t *digits);
    AliRICHDigit(Int_t *tracks, Int_t *charges, Int_t *digits);
    virtual ~AliRICHDigit() {}
    
    virtual Int_t    PadX()               {return fPadX;}
    virtual Int_t    PadY()               {return fPadY;}
    virtual Int_t    Signal()             {return fSignal;}
    virtual Int_t    Physics()            {return fPhysics;}
    virtual Int_t    Hit()                {return fHit;}    
    virtual Int_t    Track(Int_t i)       {return fTracks[i];}
    virtual Int_t    TrackCharge(Int_t i) {return fTcharges[i];}    
    virtual void     AddSignal(Int_t q)   {fSignal += q;}
    virtual void     AddPhysicsSignal(Int_t q)   {fPhysics += q;}	
            void     Print(Option_t *option)const;      //virtual
    
    ClassDef(AliRICHDigit,1)  //Digits for set:RICH
};
#endif







