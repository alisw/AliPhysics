#ifndef SHIL_Hv0
#define SHIL_Hv0
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//  Manager class for Module: SHIL          //
////////////////////////////////////////////////
 
#include "AliSHIL.h"
 
 
class AliSHILv0 : public AliSHIL {
  
public:
  AliSHILv0();
  AliSHILv0(const char *name, const char *title);
  virtual      ~AliSHILv0() {}
  virtual void  CreateGeometry();
  virtual void  Init();
  virtual Int_t IsVersion() const {return 0;}
  ClassDef(AliSHILv0,1)  // Muon Shield Class (Open Geometry)
};

#endif
