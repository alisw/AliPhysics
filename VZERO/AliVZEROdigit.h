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
    AliVZEROdigit(Int_t   /* PMnumber */, Float_t  /* ADC */, Float_t /* Time */);
    AliVZEROdigit(Int_t   /* PMnumber */, Float_t  /* ADC */, Float_t /* Time */, 
                  Float_t /* TimeWidth*/, Bool_t /* BBFlag */, Bool_t /* BGFlag */);
    AliVZEROdigit(Int_t   /* PMnumber */, Float_t  /* ADC */, Float_t /* Time */, 
                  Float_t /* TimeWidth*/, Bool_t /* BBFlag */, Bool_t /* BGFlag */, Bool_t /* Integrator */);
    virtual ~AliVZEROdigit() {};
    virtual void Print(const Option_t* option="") const;
    
    Int_t   PMNumber()   const {return fPMNumber;}    
    Float_t ADC()        const {return fADC;}
    Float_t Time()       const {return fTime;}
    Float_t Width()      const {return fWidth;} 
    Bool_t  BBFlag()     const {return fBBFlag;} 
    Bool_t  BGFlag()     const {return fBGFlag;}
    Bool_t  Integrator() const {return fIntegrator;}
       
  private:
    Int_t  fTrack;         // Track number
    
  protected:
    Int_t   fEvent;         // Event number  
    Int_t   fPMNumber;      // PhotoMultiplier number (0 to 63)
    Float_t fADC;           // ADC response
    Float_t fTime;          // Time of Flight
    Float_t fWidth;         // Width of the time distribution
    Bool_t  fBBFlag;        // Beam-Beam Flag given by Yannick in Raw Data only
    Bool_t  fBGFlag;        // Beam-Gas  Flag given by Yannick in Raw Data only
    Bool_t  fIntegrator;    // Integrator used
    
    ClassDef(AliVZEROdigit,3)  //Digit (Header) object for set : VZERO
};

#endif
