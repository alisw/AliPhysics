#ifndef ALIRICHPADHIT_H
#define ALIRICHPADHIT_H


/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include <TObject.h>

class AliRICHPadHit : public TObject {
 public:
    
    Int_t     fHitNumber;    // Hit number
    Int_t     fQpad  ;          // Total charge      
    Int_t     fPadX  ;       // Pad number along X
    Int_t     fPadY  ;       // Pad number along Y
    Int_t     fRSec  ;       // R -sector of pad
    
 public:
    AliRICHPadHit() {
	fHitNumber=fPadX=fPadY=fQpad=fRSec=0;   
    }
    AliRICHPadHit(Int_t *clhits);
    virtual ~AliRICHPadHit() {}
       
    ClassDef(AliRICHPadHit,1)  //Cluster object for set:RICH
};
#endif
