#ifndef ALIDIPOV2_H
#define ALIDIPOV2_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/////////////////////////////////////////////////
//  Manager class for detector: DIPO version 2 //
/////////////////////////////////////////////////
 
#include "AliDIPO.h"
  
class AliDIPOv2 : public AliDIPO {
  
public:
  AliDIPOv2();
  AliDIPOv2(const char *name, const char *title);
  virtual      ~AliDIPOv2() {}
  virtual void  CreateGeometry();
  virtual void  CreateMaterials();
  virtual Int_t IsVersion() const {return 2;}
 private:
  virtual void  CreateSpectrometerDipole();
  virtual void  CreateCompensatorDipole();
  
  ClassDef(AliDIPOv2,1)  //Class manager for magnetic dipole version 2
};

#endif
