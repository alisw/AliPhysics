#ifndef ALIVZERODIGIT_H
#define ALIVZERODIGIT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "AliDigit.h"

//_____________________________________________________________________________
class AliVZEROdigit: public AliDigit  {

 public:
    AliVZEROdigit();
    AliVZEROdigit(Int_t* tracks, Int_t* digits);
    AliVZEROdigit(Int_t /* PMnumber */, Int_t /* ADC */, Int_t /* Time */);
    virtual ~AliVZEROdigit() {};
    Int_t   PMNumber()  const {return fPMNumber;}    
    Int_t   ADC()  const {return fADC;}
    Int_t   Time() const {return fTime;}
     
  private:
    Int_t  fTrack;         // Track number
    
  protected:
    Int_t  fEvent;         // Event number  
    Int_t  fPMNumber;      // Photomultiplier number (0 to 63)
    Int_t  fADC;           // ADC response
    Int_t  fTime;          // Time of Flight
    
    ClassDef(AliVZEROdigit,1)  //Digit (Header) object for set : VZERO
};

#endif
