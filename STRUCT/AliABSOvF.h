#ifndef ALIABSOVF_H
#define ALIABSOVF_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//  Manager class for detector: ABSO          //
////////////////////////////////////////////////
 
#include "AliABSO.h"
 
 
class AliABSOvF : public AliABSO {
 
public:
  AliABSOvF();
  AliABSOvF(const char *name, const char *title);
  virtual      ~AliABSOvF() {}
  virtual void  CreateGeometry();
  virtual Int_t IsVersion() const {return 0;}
  virtual void  Init();
  
  ClassDef(AliABSOvF,1)  // Muon Absorber Class
};

#endif
