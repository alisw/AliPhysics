#ifndef STARTV0_H
#define STARTV0_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//  Manager and hits classes for set:START     //
////////////////////////////////////////////////
 
#include "AliSTART.h"
 
class AliSTARTv0 : public AliSTART {
  
public:
  AliSTARTv0() {};
  AliSTARTv0(const char *name, const char *title);
  virtual       ~AliSTARTv0() {}
  virtual void   CreateGeometry();
  virtual void   CreateMaterials();
  virtual void   DrawModule();
  virtual void   Init();
  virtual Int_t  IsVersion() const {return 0;}
  virtual void   StepManager();
  
protected:
   Int_t fIdSens1; // Sensetive volume  in START
 
  ClassDef(AliSTARTv0,1)  //Class for START version 0
};

#endif


