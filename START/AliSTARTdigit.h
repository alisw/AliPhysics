#ifndef STARTDIGIT_H
#define STARTDIGIT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "AliDigit.h"
#include "AliSTART.h"

//___________________________________________
class AliSTARTdigit: public AliDigit  {
////////////////////////////////////////////////////////////////////////
 public:
    Int_t fEvent;            // Event number
    Int_t fTime_average;     // Average time
    Int_t fTime_diff;  // Time difference

 public:
    AliSTARTdigit() {}
    AliSTARTdigit(Int_t *tracks, Int_t *digits);
    virtual ~AliSTARTdigit() {}

    ClassDef(AliSTARTdigit,1)  //Digit (Header) object for set:ITS
};
#endif
