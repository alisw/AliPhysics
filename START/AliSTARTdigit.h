#ifndef ALISTARTDIGIT_H
#define ALISTARTDIGIT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
#include <TObject.h>
#include "AliSTART.h"
//___________________________________________
class AliSTARTdigit: public TObject  {
////////////////////////////////////////////////////////////////////////
 protected:
    Int_t fTimeAverage;     // Average time
    Int_t fTimeDiff;  // Time difference


 public:
    AliSTARTdigit();
    AliSTARTdigit(Int_t , Int_t );
    virtual ~AliSTARTdigit() {}
    void Set(Int_t, Int_t);
    Int_t GetTime();
    void MyDump(); 


    ClassDef(AliSTARTdigit,1)  //Digit (Header) object for set:START
};

inline AliSTARTdigit::AliSTARTdigit(){fTimeAverage=999999;fTimeDiff=999999;}
inline Int_t AliSTARTdigit::GetTime(){return fTimeDiff;}
inline void AliSTARTdigit::Set(Int_t Timeav, Int_t Timediff)
  {fTimeAverage=Timeav; fTimeDiff=Timediff;}

inline void AliSTARTdigit::MyDump(){
  printf("AliSTARTdigit: fTimeAverage=%d, fTimeDiff=%d\n",
	 fTimeAverage, fTimeDiff);
}

#endif



