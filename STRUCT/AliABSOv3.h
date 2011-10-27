#ifndef ALIABSOV3_H
#define ALIABSOV3_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id $ */

////////////////////////////////////////////////
//  Manager class for Module: ABSO          //
////////////////////////////////////////////////
 
#include "AliABSO.h"
class TGeoPcon;

 
class AliABSOv3: public AliABSO {
  
public:
  AliABSOv3();
  AliABSOv3(const char *name, const char *title);
  virtual      ~AliABSOv3() {}
  virtual void  Init(){;}
  virtual void  CreateGeometry();
 private:
  virtual TGeoPcon* MakeShapeFromTemplate(const TGeoPcon* pcon, Float_t drIn, Float_t drOut);  
  ClassDef(AliABSOv3,1)  // Front Absorber as built
};

#endif
