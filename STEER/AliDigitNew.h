#ifndef ALIDIGITNEW_H
#define ALIDIGITNEW_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//  Base class for Alice Digits               //
////////////////////////////////////////////////

#include "TObject.h"

class AliDigitNew : public TObject {

 public: 
  AliDigitNew() ;   
  ~AliDigitNew() {;}
  Int_t   GetAmp() const { return fAmp  ; } 
  Int_t   GetId() const { return fId ; }     
   
 protected:
 
  Int_t fId ;                // absolute id
  Int_t fAmp ;               // digitalized energy
     
  ClassDef(AliDigitNew,1)  //Base class for all Alice digits

} ;
#endif // ALIDIGITNEW_H
