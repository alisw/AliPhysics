#ifndef STARTV1_H
#define STARTV1_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
////////////////////////////////////////////////
//  Manager and hits classes for set:START     //
////////////////////////////////////////////////
 
#include "AliSTART.h"
 
class AliSTARTv1 : public AliSTART {
  
public:
  AliSTARTv1() {};
  AliSTARTv1(const char *name, const char *title);
  virtual       ~AliSTARTv1() {}
  virtual void   CreateGeometry();
  virtual void   CreateMaterials();
  virtual void   DrawDetector();
  virtual void   Init();
  virtual Int_t  IsVersion() const {return 0;}
  virtual void   StepManager();
  
protected:
   Int_t fIdSens1; // Sensetive volume  in START
 
  ClassDef(AliSTARTv1,1)  //Class for START version 0
};

#endif


