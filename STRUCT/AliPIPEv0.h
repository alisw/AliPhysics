#ifndef ALIPIPEV0_H
#define ALIPIPEV0_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/////////////////////////////////////////////////////////
//  Manager and class for detector: PIPE  version 0    //
/////////////////////////////////////////////////////////
 
#include "AliPIPE.h"

class AliPIPEv0 : public AliPIPE {
  
 public:
    enum constants {kC=6, kAlu=9, kInox=19, kGetter=20, kBe=5, kVac=16,
	  kAir=15, kAlBe=21, kPA = 22};
	
  AliPIPEv0();
  AliPIPEv0(const char *name, const char *title);
  virtual       ~AliPIPEv0() {}
  virtual void   CreateGeometry();
  virtual void   CreateMaterials();
  virtual void   SetPipeMaterial(Int_t mat = kBe) {fPipeMaterial = mat;}
  virtual Int_t  IsVersion() const {return 0;}
 protected:
  Int_t   fPipeMaterial; // Pipe material (Al, Be, or Inox)
  
  ClassDef(AliPIPEv0,2)  //Class for PIPE version 0
};
 
#endif
