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
    AliVZEROdigit(Int_t /* cellnumber */, Int_t /* adc */);
    virtual ~AliVZEROdigit() {};
    Int_t   CellNumber()  const {return fCellNumber;}    
    Int_t   ADC() const {return fADC;}
     
  private:
    Int_t  fTrack;         // Track number
    
  protected:
    Int_t  fEvent;         // Event number  
    Int_t  fCellNumber;    // Scintillator cell number
    Int_t  fADC;           // ADC response
    
    ClassDef(AliVZEROdigit,1)  //Digit (Header) object for set : VZERO
};

#endif
