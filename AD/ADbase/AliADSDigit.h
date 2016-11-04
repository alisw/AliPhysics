#ifndef ALIADSDIGIT_H
#define ALIADSDIGIT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliDigit.h"

//_____________________________________________________________________________
class AliADSDigit: public AliDigit  {

public:
  AliADSDigit();
  AliADSDigit(Int_t pmnumber,
	      Int_t nbins,
	      Float_t *charges,
	      Int_t *labels = 0);
  virtual            ~AliADSDigit();

  virtual void        Print(const Option_t* option="") const;

  Int_t               PMNumber()   const {return fPMNumber;}
  Int_t               GetNBins()   const {return fNBins;}
  Float_t*            GetCharges() const {return fCharges;}

  virtual void        Clear(Option_t*);

private:
  AliADSDigit(const AliADSDigit& /*sdigit*/);
  AliADSDigit& operator = (const AliADSDigit& /*sdigit*/);

  Int_t               fPMNumber;      // PhotoMultiplier number (0 to 16)
  Int_t               fNBins;         // Number of charge bins
  Float_t*            fCharges;       //[fNBins] Array with charges

  ClassDef(AliADSDigit,1);  // AD SDigit class
};

#endif
