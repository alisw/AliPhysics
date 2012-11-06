#ifndef ALIITSVUPGRADE_H
#define ALIITSVUPGRADE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//========================================================================
//
//        Geometry for the Upgrade of the Inner Tracking System
//
// Mario Sitta (sitta@to.infn.it)
//
//========================================================================


// $Id: AliITSUv11.h 

#include "AliITSU.h"

class  AliITSUv11Layer;
class  AliITSv11GeomBeamPipe;
class  TGeoVolume;
class  TGeoVolumeAssembly;

class AliITSUv11 : public AliITSU {

 public:
  AliITSUv11();
  AliITSUv11(const char *title, const Int_t nlay);
  virtual       ~AliITSUv11() ;
  
  virtual void   AddAlignableVolumes() const;
  virtual void   CreateGeometry();
  virtual void   CreateMaterials();
  virtual void   DefineLayer(Int_t nlay,Double_t phi0,Double_t r,Double_t zlen,Int_t nladd,
			     Int_t nmod, Double_t lthick=0.,Double_t dthick=0.,UInt_t detType=0);
  virtual void   DefineLayerTurbo(Int_t nlay,Double_t phi0,Double_t r,Double_t zlen,Int_t nladd,
				  Int_t nmod,Double_t width,Double_t tilt,
				  Double_t lthick = 0.,Double_t dthick = 0.,UInt_t detType=0);
  virtual void   GetLayerParameters(Int_t nlay, Double_t &phi0,Double_t &r, Double_t &zlen,
				    Int_t &nladd, Int_t &nmod,
				    Double_t &width, Double_t &tilt,
				    Double_t &lthick, Double_t &mthick) const;
  virtual void   Init(); 
  virtual Bool_t IsLayerTurbo(Int_t nlay);
  virtual Int_t  IsVersion()                 const { return 20;}  // vUpgrade ? do we need this
  virtual void   SetDefaults();
  virtual void   StepManager();
  virtual void   SetLayerDetTypeID(Int_t lr, UInt_t id);
  virtual Int_t  GetLayerDetTypeID(Int_t lr);
  //
 private:
  AliITSUv11(const AliITSUv11 &source); // copy constructor
  AliITSUv11& operator=(const AliITSUv11 &source); // assignment operator
  //
  Bool_t   *fLayTurbo;       // True for "turbo" layers
  Double_t *fLayPhi0;        // Vector of layer's 1st ladder phi in lab
  Double_t *fLayRadii;       // Vector of layer radii
  Double_t *fLayZLength;     // Vector of layer length along Z
  Int_t    *fLaddPerLay;     // Vector of number of ladders per layer
  Int_t    *fModPerLadd;     // Vector of number of modules per ladder
  Double_t *fLadThick;       // Vector of ladder thicknesses
  Double_t *fLadWidth;       // Vector of ladder width (only used for turbo)
  Double_t *fLadTilt;        // Vector of ladder tilt (only used for turbo)
  Double_t *fDetThick;       // Vector of detector thicknesses
  UInt_t   *fDetTypeID;      // Vector of detector type id
  //  
  AliITSUv11Layer **fUpGeom; //! Geometry
  
  // Parameters for the Upgrade geometry
  
  ClassDef(AliITSUv11,0)                          
};
 
#endif
