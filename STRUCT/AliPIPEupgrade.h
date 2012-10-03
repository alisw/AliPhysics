#ifndef ALIPIPEUPGRADE_H
#define ALIPIPEUPGRADE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliPIPEupgrade.h 55648 2012-04-10 10:59:00Z cvetan $*/

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


class AliPIPEupgrade : public AliPIPE {
    
 public:
  enum constants {kC=6, kAlu=9, kInox=19, kGetter=20, kBe=5, kVac=16,
		  kAir=15, kAlBe=21, kPA = 22};
  
  AliPIPEupgrade(Bool_t coneIsBe=0, Float_t ro=1.8, Float_t width=0.08, Float_t hlength=61.52);
  AliPIPEupgrade(const char *name, const char *title, 
		 Bool_t coneIsBe=0, Float_t ro=1.8, Float_t width=0.08, Float_t hlength=61.52);
  AliPIPEupgrade(Option_t *opt);
  // virtual       ~AliPIPEupgrade() {}
  virtual void   CreateGeometry();
  virtual void   CreateMaterials();
  virtual Int_t  IsVersion() const {return 0;}
  virtual void   SetBeamBackgroundSimulation() {fBeamBackground = kTRUE;}
  virtual void   AddAlignableVolumes() const;
  
  Float_t GetRmin() {return fIpPipeRo-fIpPipeWidth; }
  Float_t GetRmax() {return fIpPipeRo; }
  Float_t GetWidth(){return fIpPipeWidth; }
  Float_t GetDz()   {return fIpHLength; } 
	  
 protected:
  virtual TGeoPcon*   MakeMotherFromTemplate(const TGeoPcon* shape, Int_t imin = -1, Int_t imax = -1, Float_t r0 = 0., Int_t nz =-1);
  virtual TGeoPcon*   MakeInsulationFromTemplate(TGeoPcon* shape);
  virtual TGeoVolume* MakeBellow(const char* ext, Int_t nc, Float_t rMin, Float_t rMax, Float_t dU, Float_t rPlie, Float_t dPlie);
  Bool_t  fBeamBackground; // Flag for beam background simulations
  
  Bool_t  fConeIsBe;         // Flag for Material of the "long cone", can be set to Be (std is Alu)
  Float_t fIpPipeRo;         // outer diameter of the beampipe around the IP
  Float_t fIpPipeWidth;      // width of the beampipe around the IP
  Float_t fIpHLength;        // half length of the beampipe around the IP // FixMe: up to now, hardcoded to 57.25cm
 
  ClassDef(AliPIPEupgrade, 1)  // Class for PIPE version using TGeo
};
 
#endif
 
