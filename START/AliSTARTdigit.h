#ifndef STARTDIGIT_H
#define STARTDIGIT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
#include <TObject.h>
#include "AliSTART.h"
//___________________________________________
class AliSTARTdigit: public TObject  {
////////////////////////////////////////////////////////////////////////
 public:
    Int_t fTime_average;     // Average time
    Int_t fTime_diff;  // Time difference


 public:
    AliSTARTdigit();
    AliSTARTdigit(Int_t , Int_t );
    virtual ~AliSTARTdigit() {}
    void Set(Int_t, Int_t);
    void MyDump(); 


    ClassDef(AliSTARTdigit,1)  //Digit (Header) object for set:START
};

inline AliSTARTdigit::AliSTARTdigit(){fTime_average=99999.;fTime_diff=99999.;}
inline void AliSTARTdigit::Set(Int_t Timeav, Int_t Timediff)
  {fTime_average=Timeav; fTime_diff=Timediff;}

inline void AliSTARTdigit::MyDump(){
  printf("AliSTARTdigit: fTime_average=%d, fTime_diff=%d\n",
	 fTime_average, fTime_diff);
}

#endif



