#ifndef ALIVZERODIGIT_H
#define ALIVZERODIGIT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "AliDigit.h"

//_____________________________________________________________________________
class AliVZEROdigit: public AliDigit  {

 public:
    AliVZEROdigit() {};
    AliVZEROdigit(Int_t* tracks, Int_t* digits);
    virtual ~AliVZEROdigit() {};
    
  private:
    Int_t fEvent;         // Event number
    Int_t fTrack;         // Track number
   
    ClassDef(AliVZEROdigit,1)  //Digit (Header) object for set : VZERO
};

#endif
