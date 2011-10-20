#ifndef ALIPIPEVGEO4_H
#define ALIPIPEVGEO4_H
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


class AliPIPEv4 : public AliPIPE {
    
 public:
    enum constants {kC=6, kAlu=9, kInox=19, kGetter=20, kBe=5, kVac=16,
	  kAir=15, kAlBe=21, kPA = 22};
	
  AliPIPEv4();
  AliPIPEv4(const char *name, const char *title);
  AliPIPEv4(const char *name, const char *title, const Float_t theta_cone,  const Float_t rmin1, 
	    const Float_t epaisseur, const Float_t sigmaz, const Float_t z_chambre);
  virtual       ~AliPIPEv4() {};
  virtual void   CreateGeometry();
  virtual void   CreateMaterials();
  virtual Int_t  IsVersion() const {return 0;}
 private:
  virtual TGeoPcon*   MakeMotherFromTemplate(TGeoPcon* shape, Int_t imin = -1, Int_t imax = -1, Float_t r0 = 0., Int_t nz =-1);
  virtual TGeoPcon*   MakeInsulationFromTemplate(TGeoPcon* shape);
  virtual TGeoVolume* MakeBellow(const char* ext, Int_t nc, Float_t rMin, Float_t rMax, Float_t dU, Float_t rPlie, Float_t dPlie);

  Float_t   ftheta_cone; // angle of conical beam pipe, if angle < 3 --> cylindrical beam pipe
  Float_t   frmin1;      // internal radius of Be beam pipe
  Float_t   fepaisseur;  // width of Be beam pipe
  Float_t   fsigmaz;     // dispersion of z location (1 sigma) of beam impact position
  Float_t   fz_chambre;  // first pixel chamber location, closest to the IP
  Float_t   fzdebut1;    // beginning of beam pipe z location (A side)
  Float_t   fzfin4;      // end of beamp pipe z location (C side)


 protected:
  ClassDef(AliPIPEv4,1)  //Class for PIPE version using TGeo
};
 
#endif
