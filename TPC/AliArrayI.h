#ifndef ALIARRAYI_H
#define ALIARRAYI_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//  Manager class generaol Alice segment 
//  segment is for example one pad row in TPC //
////////////////////////////////////////////////

#include "TObject.h"
#include "TArrayI.h"

class AliArrayI: public TObject ,public TArrayI {
public:
  void Expand(Int_t n);
  ClassDef(AliArrayI,1) 
};

#endif //ALIARRAY_I

