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
    AliVZEROdigit(Int_t /* PMnumber */, Int_t  /* ADC */, Int_t /* Time */);
    AliVZEROdigit(Int_t /* PMnumber */, Int_t  /* ADC */, Int_t /* Time */, 
                  Int_t /* TimeWidth*/, Bool_t /* BBFlag */, Bool_t /* BGFlag */);
    virtual ~AliVZEROdigit() {};
    virtual void Print(const Option_t* option="") const;
    
    Int_t   PMNumber() const {return fPMNumber;}    
    Int_t   ADC()      const {return fADC;}
    Int_t   Time()     const {return fTime;}
    Int_t   Width()    const {return fWidth;} 
    Bool_t  BBFlag()   const {return fBBFlag;} 
    Bool_t  BGFlag()   const {return fBGFlag;}
       
  private:
    Int_t  fTrack;         // Track number
    
  protected:
    Int_t  fEvent;         // Event number  
    Int_t  fPMNumber;      // PhotoMultiplier number (0 to 63)
    Int_t  fADC;           // ADC response
    Int_t  fTime;          // Time of Flight
    Int_t  fWidth;         // Width of the time distribution
    Bool_t fBBFlag;        // Beam-Beam Flag given by Yannick in Raw Data only
    Bool_t fBGFlag;        // Beam-Gas  Flag given by Yannick in Raw Data only
    
    ClassDef(AliVZEROdigit,2)  //Digit (Header) object for set : VZERO
};

#endif
