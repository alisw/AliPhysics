#ifndef AliDigit_H
#define AliDigit_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//  Base class for Alice Digits               //
////////////////////////////////////////////////

#include "TObject.h"

class AliDigit : public TObject {
public:
  Int_t     fTracks[3];   //tracks number making this digit (up to 3)

public:
  AliDigit();
  AliDigit(Int_t *track);
  ~AliDigit() {;}
  inline virtual int *GetTracks() {return &fTracks[0];}
  
  ClassDef(AliDigit,1)  //Base class for all Alice digits
};
#endif
