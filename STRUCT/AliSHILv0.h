#ifndef ALISHILV0_H
#define ALISHILV0_H
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
  virtual void  SetPbCone(Bool_t flag=kTRUE) {fPbCone=flag;}
	  
 protected:
  Bool_t fPbCone;
  
  ClassDef(AliSHILv0,1)  // Muon Shield Class (Open Geometry)
      
};

#endif
