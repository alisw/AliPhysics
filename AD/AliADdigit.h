#ifndef ALIADDIGIT_H
#define ALIADDIGIT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliADdigit.h  $ */

#include "AliDigit.h"

//_____________________________________________________________________________
class AliADdigit: public AliDigit  {

public:
                   AliADdigit();
		   AliADdigit(Int_t* tracks, Int_t module, Float_t cell);
		   AliADdigit(Int_t* module, Float_t cell);
		   AliADdigit(Int_t module, Float_t cell);
  virtual          ~AliADdigit();
  virtual void      Print(const Option_t* option="") const;

  Int_t GetModule() const {return fModule;}
  Float_t GetCell() const {return fCell;}

    
private:
  Int_t fModule; //! module producing the digit 
  Float_t           fCell;                  // Time of Flight

  ClassDef(AliADdigit,1)  // AD Digit class
};
//typedef AliADdigit AliADdigit;
#endif
