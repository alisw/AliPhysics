#ifndef ABSO_H
#define ABSO_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//  Manager class for detector: ABSO          //
////////////////////////////////////////////////
 
#include "AliModule.h"
 
 
class AliABSO : public AliModule {
 
public:
  AliABSO();
  AliABSO(const char *name, const char *title);
  virtual      ~AliABSO() {}
  virtual void  CreateGeometry();
  virtual void  CreateMaterials();
  virtual void  Init();
  virtual Int_t IsVersion() const {return 0;}
  virtual void  DrawModule();
  
  
  ClassDef(AliABSO,1)  // Muon Absorber Class
};

#endif
