#ifndef ALIDIGIT_H
#define ALIDIGIT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//  Base class for Alice Digits               //
////////////////////////////////////////////////

#include <assert.h>
#include "TObject.h"

class AliDigit : public TObject {
public:
  AliDigit();
  AliDigit(Int_t *track);
  virtual ~AliDigit() {}
  virtual Int_t *GetTracks() {return fTracks;}
  virtual Int_t GetTrack(Int_t i) const 
  {assert(0<=i&&i<=2); return fTracks[i];}
  
protected:
  Int_t     fTracks[3];   //tracks number making this digit (up to 3)

  ClassDef(AliDigit,1)  //Base class for all Alice digits
};
#endif
