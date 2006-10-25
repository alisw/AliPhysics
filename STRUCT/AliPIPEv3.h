#ifndef ALIPIPEVGEO_H
#define ALIPIPEVGEO_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$*/

/////////////////////////////////////////////////////////
// ALICE beam pipe geometry                            //
// This version uses TGeo.                             //
// Author:                                             //
// Andreas Morsch                                      //
// e-mail: andreas.morsch@cern.ch                      // 
/////////////////////////////////////////////////////////
 
#include "AliPIPE.h"
class TGeoPcon;
class TGeoVolume;


class AliPIPEv3 : public AliPIPE {
    
 public:
    enum constants {kC=6, kAlu=9, kInox=19, kGetter=20, kBe=5, kVac=16,
	  kAir=15, kAlBe=21, kPA = 22};
	
  AliPIPEv3();
  AliPIPEv3(const char *name, const char *title);
  virtual       ~AliPIPEv3() {}
  virtual void   CreateGeometry();
  virtual void   CreateMaterials();
  virtual Int_t  IsVersion() const {return 0;}
 private:
  virtual TGeoPcon*   MakeMotherFromTemplate(TGeoPcon* shape, Int_t imin = -1, Int_t imax = -1, Float_t r0 = 0., Int_t nz =-1);
  virtual TGeoPcon*   MakeInsulationFromTemplate(TGeoPcon* shape);
  virtual TGeoVolume* MakeBellow(char* ext, Int_t nc, Float_t rMin, Float_t rMax, Float_t dU, Float_t rPlie, Float_t dPlie);
 protected:
  ClassDef(AliPIPEv3,1)  //Class for PIPE version using TGeo
};
 
#endif
