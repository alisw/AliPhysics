#ifndef ABSOv0_H
#define ABSOv0_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//  Manager class for detector: ABSO          //
////////////////////////////////////////////////
 
#include "AliABSO.h"
 
 
class AliABSOv0 : public AliABSO {
 
public:
  AliABSOv0();
  AliABSOv0(const char *name, const char *title);
  virtual      ~AliABSOv0() {}
  virtual void  CreateGeometry();
  virtual Int_t IsVersion() const {return 0;}
  virtual void  Init();
  
  ClassDef(AliABSOv0,1)  // Muon Absorber Class
};

#endif
