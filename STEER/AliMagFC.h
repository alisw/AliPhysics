#ifndef ALIMAGFC_H
#define ALIMAGFC_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "AliMagF.h"

class AliMagFC  : public AliMagF
{
  //Alice Constant Magnetic Field

public:
  AliMagFC(){}
  AliMagFC(const char *name, const char *title, const Int_t integ, 
	   const Int_t map, const Float_t factor, const Float_t fmax);
  virtual ~AliMagFC() {}
  virtual void Field(Float_t *x, Float_t *b);
  virtual void ReadField() {}
  
  ClassDef(AliMagFC,1)  //Class for all Alice Constant MagField 
};

#endif
