#ifndef ALIPIPEVGEO_H
#define ALIPIPEVGEO_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/////////////////////////////////////////////////////////
// ALICE beam pipe geometry                            //
// This version uses TGeo.                             //
// Author:                                             //
// Andreas Morsch                                      //
// e-mail: andreas.morsch@cern.ch                      // 
/////////////////////////////////////////////////////////
 
#include "AliPIPE.h"

class AliPIPEvGEO : public AliPIPE {
  
 public:
    enum constants {kC=6, kAlu=9, kInox=19, kGetter=20, kBe=5, kVac=16,
	  kAir=15, kAlBe=21, kPA = 22};
	
  AliPIPEvGEO();
  AliPIPEvGEO(const char *name, const char *title);
  virtual       ~AliPIPEvGEO() {}
  virtual void   CreateGeometry();
  virtual void   CreateMaterials();
  virtual Int_t  IsVersion() const {return 0;}
 protected:
  ClassDef(AliPIPEvGEO,1)  //Class for PIPE version using TGeo
};
 
#endif
