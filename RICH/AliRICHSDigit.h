#ifndef ALIRICHSDIGIT_H
#define ALIRICHSDIGIT_H


/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include <TObject.h>

class AliRICHSDigit : public TObject {
 protected:
    
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
    Int_t   HitNumber() const {return fHitNumber;}
    Int_t   PadX() const      {return fPadX;}   
    Int_t   PadY() const      {return fPadY;}
    Int_t   QPad() const      {return fQpad;}
    Int_t   RSec() const      {return fRSec;}
    void    Print(Option_t *option)const;      //virtual
    
    virtual ~AliRICHSDigit() {}
       
    ClassDef(AliRICHSDigit,1)  //Cluster object for set:RICH
};
#endif
