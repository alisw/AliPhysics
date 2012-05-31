#ifndef ALIDIGITNEW_H
#define ALIDIGITNEW_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//  Base class for Alice Digits               //
////////////////////////////////////////////////

#include "TObject.h"
class TClonesArray;

typedef TClonesArray DigitsList ;       

class AliDigitNew : public TObject {

 public: 
  AliDigitNew();   
  ~AliDigitNew() {}
  Int_t   GetAmp() const { return fAmp  ; } 
  Int_t   GetId() const  { return fId ; }      
  Int_t   GetIndexInList() const { return fIndexInList ; } 
  void    SetIndexInList(Int_t val) { fIndexInList = val ; } 

 protected:
 
  Int_t fAmp ;               // digitalized energy
  Int_t fId ;                // absolute id
  Int_t fIndexInList ;       // the index of this digit in the list stored in TreeD
   
  ClassDef(AliDigitNew,1)  //Base class for all Alice digits

} ;
#endif // ALIDIGITNEW_H
