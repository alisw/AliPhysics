#ifndef DIPOv1_H
#define DIPOv1_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/////////////////////////////////////////////////
//  Manager class for detector: DIPO version 1 //
/////////////////////////////////////////////////
 
#include "AliDIPO.h"
 
class AliDIPOv1 : public AliDIPO {
 
public:
  AliDIPOv1();
  AliDIPOv1(const char *name, const char *title);
  virtual      ~AliDIPOv1() {}
  virtual void  CreateGeometry();
  virtual void  CreateMaterials();
  virtual Int_t IsVersion() const {return 1;}
  virtual void  DrawModule();
  
  ClassDef(AliDIPOv1,1)  //Class for the Magnetic Dipole version 1
};

#endif
