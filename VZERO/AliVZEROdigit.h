#ifndef VZERODIGIT_H
#define VZERODIGIT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "AliDigit.h"
#include "AliVZERO.h"

//___________________________________________
class AliVZEROdigit: public AliDigit  {

 public:
    Int_t fEvent;         // Event number
    Int_t fMulti;         // Multiplicity of charged particles
    Int_t fTrack;
    Int_t fNCerenkovs;
 
 public:
    AliVZEROdigit() {}
    AliVZEROdigit(Int_t* tracks, Int_t* digits);
    virtual ~AliVZEROdigit() {}

    ClassDef(AliVZEROdigit,1)  //Digit (Header) object for set : VZERO
};
#endif
