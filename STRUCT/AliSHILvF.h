#ifndef ALISHILVF_H
#define ALISHILVF_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//  Manager class for Module: SHIL          //
////////////////////////////////////////////////
 
#include "AliSHIL.h"
 
 
class AliSHILvF : public AliSHIL {
  
public:
  AliSHILvF();
  AliSHILvF(const char *name, const char *title);
  virtual      ~AliSHILvF() {}
  virtual void  CreateGeometry();
  virtual void  Init();
  virtual Int_t IsVersion() const {return 0;}
  virtual void  SetPbCone(Bool_t flag=kTRUE) {fPbCone=flag;}
	  
 protected:
  Bool_t fPbCone;         // outer Pb cone option flag 
  
  ClassDef(AliSHILvF,1)  // Muon Shield Class (Open Geometry)
      
};

#endif
