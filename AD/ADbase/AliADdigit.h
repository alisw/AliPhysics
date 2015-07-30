#ifndef ALIADDIGIT_H
#define ALIADDIGIT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliADdigit.h  $ */

#include "AliDigit.h"
#include "AliADConst.h"

//_____________________________________________________________________________
class AliADdigit: public AliDigit  {

public:
    AliADdigit();
    AliADdigit(Int_t   PMnumber, Float_t time, 
                  Float_t TimeWidth,
		  Bool_t  Integrator,
		  Short_t *chargeADC = 0,
		  Bool_t  BBflag = kFALSE,
		  Bool_t  BGflag = kFALSE,
		  Int_t *labels = 0);
    virtual ~AliADdigit() {};
    virtual void Print(const Option_t* option="") const;

    Int_t   PMNumber()   const {return fPMNumber;}    
    Short_t ADC()        const {return fChargeADC[kADNClocks/2];}
    Float_t Time()       const {return fTime;}
    Float_t Width()      const {return fWidth;} 
    Bool_t  Integrator() const {return fIntegrator;}
    Short_t ChargeADC(Int_t clock) const {return (clock >= 0 && clock < kADNClocks) ? fChargeADC[clock] : 0;}
    Bool_t  GetIntegratorFlag(Int_t clock);
    Bool_t  GetBBflag()  const {return fBBflag;}
    Bool_t  GetBGflag()  const {return fBGflag;}
    void    SetBBflag(Bool_t bbFlag) {fBBflag = bbFlag;}
    void    SetBGflag(Bool_t bgFlag) {fBGflag = bgFlag;}
    
  protected:
    Int_t   fPMNumber;      // PhotoMultiplier number (0 to 16)
    Float_t fTime;          // Time of Flight
    Float_t fWidth;         // Width of the time distribution
    Bool_t  fIntegrator;    // Integrator used in central clock
    Short_t fChargeADC[kADNClocks]; // ADC samples as present in raw data
    Bool_t  fBBflag;	    // BB flag in central clock
    Bool_t  fBGflag;	    // BB flag in central clock
    

  ClassDef(AliADdigit,1)  // AD Digit class
};
//typedef AliADdigit AliADdigit;
#endif
