#ifndef T0V0_H
#define T0V0_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//  Manager and hits classes for set:T0     //
////////////////////////////////////////////////
 
#include "AliT0.h"
 
class AliT0v0 : public AliT0 {
 
public:


  enum constants {kAir=1, kSc=2, kVac=3, kCer=4, kGlass=6, kSteel=8, kRibber=9, kBrass=11, kLucite=12, kC=13, kPP=14, kAl=15, kOpGlass=16, kOpAir=17};

  AliT0v0() {};
  AliT0v0(const char *name, const char *title);
  virtual       ~AliT0v0() {}
  virtual void   CreateGeometry();
  virtual void   CreateMaterials();
  virtual void   DrawModule() const;
  virtual void   Init();
  virtual Int_t  IsVersion() const {return 0;}
  virtual void   StepManager();
  
protected:
   Int_t fIdSens1; // Sensetive volume  in T0
 
  ClassDef(AliT0v0,1)  //Class for T0 version 0
};

typedef AliT0v0 AliSTARTv0; // for backward compatibility

#endif


