#ifndef AliTRDArrayF_H
#define AliTRDArrayF_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTRDarrayF.h,v */

#include "TObject.h"
#include "TArrayF.h"

class AliTRDarrayF: public TObject ,public TArrayF {

public:

  ~AliTRDarrayF();
  void Expand(Int_t n);  

  ClassDef(AliTRDarrayF,1)  

};

#endif 

