#ifndef ALISHIL_H
#define ALISHIL_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//  Manager class for Module: SHIL          //
////////////////////////////////////////////////
 
#include "AliModule.h"
 
 
class AliSHIL : public AliModule {
  
public:
  AliSHIL();
  AliSHIL(const char *name, const char *title);
  virtual      ~AliSHIL() {}
  virtual void  CreateGeometry();
  virtual void  CreateMaterials();
  virtual void  Init();
  virtual Int_t IsVersion() const {return 0;}
 
  ClassDef(AliSHIL,1)  // Muon Shield Class
};

#endif
