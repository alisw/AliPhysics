#ifndef ALIRICHSDIGIT_H
#define ALIRICHSDIGIT_H


/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include <TObject.h>

class AliRICHSDigit : public TObject {
 public:
    
    Int_t     fHitNumber;    // Hit number
    Int_t     fQpad  ;          // Total charge      
    Int_t     fPadX  ;       // Pad number along X
    Int_t     fPadY  ;       // Pad number along Y
    Int_t     fRSec  ;       // R -sector of pad
    
 public:
    AliRICHSDigit() {
	fHitNumber=fPadX=fPadY=fQpad=fRSec=0;   
    }
    AliRICHSDigit(Int_t *clhits);
    Int_t   HitNumber()  {return fHitNumber;}
    Int_t   PadX()       {return fPadX;}   
    Int_t   PadY()       {return fPadY;}
    Int_t   QPad()       {return fQpad;}
    Int_t   RSec()       {return fRSec;}
    
    virtual ~AliRICHSDigit() {}
       
    ClassDef(AliRICHSDigit,1)  //Cluster object for set:RICH
};
#endif
