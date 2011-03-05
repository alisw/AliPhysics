#ifndef ALISHILV3_H
#define ALISHILV3_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id $ */

////////////////////////////////////////////////
//  Manager class for Module: SHIL          //
////////////////////////////////////////////////
 
#include "AliSHIL.h"
class TGeoPcon;

 
class AliSHILv3 : public AliSHIL {
  
public:
  AliSHILv3();
  AliSHILv3(const char *name, const char *title);
  virtual      ~AliSHILv3() {}
  virtual void  CreateGeometry();
  virtual void  Init();
 private:
  virtual void InvertPcon(TGeoPcon* pcon);
  virtual TGeoPcon* MakeShapeFromTemplate(const TGeoPcon* pcon, Float_t drIn, Float_t drOut);  
 protected:
  ClassDef(AliSHILv3,1)  // Small angle absorber as built
};

#endif
