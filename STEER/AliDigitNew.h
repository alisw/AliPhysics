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
  Int_t     fTracks[3];   //tracks number making this digit (up to 3)
  
 public: 
  AliDigitNew() ;   
  AliDigitNew(Int_t *track);
  ~AliDigitNew() {;}
  inline virtual int *GetTracks() {return &fTracks[0];}
  inline virtual Int_t GetAmp() = 0 ;
  
 private:
  

 public:  
  ClassDef(AliDigitNew,1)  //Base class for all Alice digits

} ;
#endif // ALIDIGITNEW_H
