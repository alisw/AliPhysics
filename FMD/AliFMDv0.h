#ifndef FMDV0_H
#define FMDV0_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//  Manager and hits classes for set:FMD     //
////////////////////////////////////////////////
 
#include "AliFMD.h"
 
class AliFMDv0 : public AliFMD {
  
public:
  AliFMDv0();
  AliFMDv0(const char *name, const char *title);
  virtual       ~AliFMDv0() {}
  virtual void   CreateGeometry();
  virtual void   CreateMaterials();
  virtual void   DrawModule();
  virtual Int_t  IsVersion() const {return 0;}
 
  ClassDef(AliFMDv0,1)  //Class for FMD version 0
};

#endif
