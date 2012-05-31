#ifndef ALIDIPOV3_H
#define ALIDIPOV3_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/////////////////////////////////////////////////
//  Manager class for detector: DIPO version 3 //
/////////////////////////////////////////////////
 
#include "AliDIPOv2.h"
  
class AliDIPOv3 : public AliDIPOv2 {
  
public:
  AliDIPOv3();
  AliDIPOv3(const char *name, const char *title);
  virtual      ~AliDIPOv3() {}
  virtual Int_t IsVersion() const {return 3;}
 private:
  virtual void  CreateSpectrometerDipole();
  ClassDef(AliDIPOv3,1)  //Class manager for magnetic dipole version 2
};

#endif
