#ifndef ALITRDARRAYI_H
#define ALITRDARRAYI_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTRDarrayI.h,v */

#include <TObject.h>
#include <TArrayI.h>

class AliTRDarrayI: public TObject, public TArrayI {

public:

  AliTRDarrayI();
  virtual ~AliTRDarrayI();
  void Copy(TObject &a);
  void Expand(Int_t n);  

  ClassDef(AliTRDarrayI,1)  // An array of integers

};

#endif 

