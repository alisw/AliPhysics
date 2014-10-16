#ifndef ALIACORDEDIGIT_H
#define ALIACORDEDIGIT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////////////////////////////////
//  ACORDE digit: Id
//
// The digits are made in FinishEvent() by summing all the hits in a 
// counter.
////////////////////////////////////////////////////////////////////////////

#include "AliDigit.h"

class AliACORDEdigit: public AliDigit  {

 public:
  AliACORDEdigit();
  AliACORDEdigit(Int_t* tracks, Int_t module, Float_t pulse_time);
  AliACORDEdigit(Int_t* modules,Float_t pulse_time);
  AliACORDEdigit(Int_t module, Float_t pulse_time);
  virtual ~AliACORDEdigit();
  virtual void Print(const Option_t* option="") const;

  Int_t GetModule() const { return fModule;}
  Float_t GetTime() const { return fTime;}

  
private:
  Int_t fModule; // module producing the digit (1-60)
  Float_t fTime; //  time of the start of the square pulse
  
  ClassDef(AliACORDEdigit,1)  //Digit (Header) object for set : ACORDE (ACORDE)

};

typedef AliACORDEdigit AliCRTdigit; // for backward compatibility

#endif // ALIACORDEDIGIT_H
