#ifndef ALIARRAYS_H
#define ALIARRAYS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//  Manager class generaol Alice segment digits
//  segment is for example one pad row in TPC //
////////////////////////////////////////////////

#include "TObject.h"
#include "TArrayS.h"


class AliArrayS:  public TObject,public TArrayS {
public:
  void Expand(Int_t n);
  ClassDef(AliArrayS,1) // Array handling
};
#endif //ALIARRAYS_H
