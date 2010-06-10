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
    AliVZEROdigit(Int_t   PMnumber, Float_t  adc, Float_t time);
    AliVZEROdigit(Int_t   PMnumber, Float_t  adc, Float_t time, 
                  Float_t TimeWidth,
		  Bool_t  Integrator,
		  Short_t *chargeADC = 0,
		  Int_t *labels = 0);
    virtual ~AliVZEROdigit() {};
    virtual void Print(const Option_t* option="") const;

    enum {kNClocks = 21};

    Int_t   PMNumber()   const {return fPMNumber;}    
    Float_t ADC()        const {return fADC;}
    Float_t Time()       const {return fTime;}
    Float_t Width()      const {return fWidth;} 
    Bool_t  Integrator() const {return fIntegrator;}
    Short_t ChargeADC(Int_t clock) const {return (clock >= 0 && clock < kNClocks) ? fChargeADC[clock] : 0;}
    
  protected:
    Int_t   fPMNumber;      // PhotoMultiplier number (0 to 63)
    Float_t fADC;           // ADC response
    Float_t fTime;          // Time of Flight
    Float_t fWidth;         // Width of the time distribution
    Bool_t  fIntegrator;    // Integrator used
    Short_t fChargeADC[kNClocks]; // ADC samples as present in raw data

    ClassDef(AliVZEROdigit,5)  // VZERO Digit class
};

#endif
