#ifndef ALIPIPEV4_H
#define ALIPIPEV4_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//-------------------------------------------------------------------------
//  Beam pipe class for ALICE MFT upgrade
//  This version uses TGeo
//  Authors:
//  F. Manso 
//  A. Morsch
//  R. Tieulent
//-------------------------------------------------------------------------

 
#include "AliPIPE.h"
class TGeoPcon;
class TGeoVolume;


class AliPIPEv4 : public AliPIPE {
    
 public:
    enum constants {kC=6, kAlu=9, kInox=19, kGetter=20, kBe=5, kVac=16, kAir=15, kAlBe=21, kPA = 22};
	
  AliPIPEv4();
  AliPIPEv4(const char *name, const char *title);
  virtual       ~AliPIPEv4() {};
  virtual void   CreateGeometry();
  virtual void   CreateMaterials();
  virtual Int_t  IsVersion() const {return 0;}
 private:
  virtual TGeoPcon*   MakeMotherFromTemplate(TGeoPcon* shape, Int_t imin = -1, Int_t imax = -1, Float_t r0 = 0., Int_t nz =-1);
  virtual TGeoPcon*   MakeInsulationFromTemplate(TGeoPcon* shape);
  virtual TGeoVolume* MakeBellow(const char* ext, Int_t nc, Float_t rMin, Float_t rMax, Float_t dU, Float_t rPlie, Float_t dPlie);
  virtual TGeoVolume* MakeBellowCside(const char* ext, Int_t nc, Float_t rMin, Float_t rMax, Float_t rPlie, Float_t dPlie);

 protected:
  ClassDef(AliPIPEv4,2)  //Class for PIPE version using TGeo
};
 
#endif
