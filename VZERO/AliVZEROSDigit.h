#ifndef ALIVZEROSDIGIT_H
#define ALIVZEROSDIGIT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliDigit.h"

//_____________________________________________________________________________
class AliVZEROSDigit: public AliDigit  {

 public:
    AliVZEROSDigit();
    AliVZEROSDigit(Int_t pmnumber,
		   Int_t nbins, 
		   Float_t *charges,
		   Int_t *labels = 0);
    virtual ~AliVZEROSDigit() {};
    virtual void Print(const Option_t* option="") const;

    Int_t   PMNumber()   const {return fPMNumber;}

  private:
    Int_t   fPMNumber;      // PhotoMultiplier number (0 to 63)
    Int_t   fNBins;         // Number of charge bins
    Float_t*fCharges;       //[fNBins] Array with charges

    ClassDef(AliVZEROSDigit,1)  // VZERO SDigit class
};

#endif
