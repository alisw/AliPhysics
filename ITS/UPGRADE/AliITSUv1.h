#ifndef ALIITSUV1_H
#define ALIITSUV1_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//========================================================================
//
//        Geometry for the Upgrade of the Inner Tracking System
//
// Mario Sitta (sitta@to.infn.it)
//
//========================================================================


// $Id: AliITSUv1.h 

#include "AliITSU.h"

class  AliITSUv1Layer;
//class  AliITSv11GeomBeamPipe;
class  TGeoVolume;
class  TGeoVolumeAssembly;

class AliITSUv1 : public AliITSU {

 public:


  typedef enum {
    kIBModelDummy=0,
    kIBModel0=1,
    kIBModel1=2, 
    kIBModel21=3,
    kIBModel22=4,
    kIBModel3=5,
    kOBModelDummy=6,
    kOBModel0=7,
    kOBModel1=8 
  } AliITSUModel_t;
  

  AliITSUv1();
  AliITSUv1(const char *title, const Int_t nlay);
  virtual       ~AliITSUv1() ;
  virtual void   SetNWrapVolumes(Int_t n);
  virtual void   AddAlignableVolumes() const;
  virtual void   CreateGeometry();
  virtual void   CreateMaterials();
  virtual void   DefineLayer(Int_t nlay,Double_t phi0,Double_t r,Double_t zlen,Int_t nladd,
			     Int_t nmod, Double_t lthick=0.,Double_t dthick=0.,UInt_t detType=0, Int_t buildFlag=0);
  virtual void   DefineLayerTurbo(Int_t nlay,Double_t phi0,Double_t r,Double_t zlen,Int_t nladd,
				  Int_t nmod,Double_t width,Double_t tilt,
				  Double_t lthick = 0.,Double_t dthick = 0.,UInt_t detType=0, Int_t buildFlag=0);
  virtual void   GetLayerParameters(Int_t nlay, Double_t &phi0,Double_t &r, Double_t &zlen,
				    Int_t &nladd, Int_t &nmod,
				    Double_t &width, Double_t &tilt,
				    Double_t &lthick, Double_t &mthick,
				    UInt_t &dettype) const;
  virtual void   DefineWrapVolume(Int_t id, Double_t rmin,Double_t rmax, Double_t zspan);
  virtual void   Init(); 
  virtual Bool_t IsLayerTurbo(Int_t nlay);
  virtual Int_t  IsVersion()                 const { return 20;}  // vUpgrade ? do we need this
  virtual void   SetDefaults();
  virtual void   StepManager();
  virtual void   SetLayerDetTypeID(Int_t lr, UInt_t id);
  virtual Int_t  GetLayerDetTypeID(Int_t lr);
  virtual void   SetStaveModelIB(AliITSUModel_t model) {fStaveModelIB=model;}
  virtual void   SetStaveModelOB(AliITSUModel_t model) {fStaveModelOB=model;}
  virtual AliITSUModel_t GetStaveModelIB() const {return fStaveModelIB;}
  virtual AliITSUModel_t GetStaveModelOB() const {return fStaveModelOB;}
  //
 private:
  AliITSUv1(const AliITSUv1 &source); // copy constructor
  AliITSUv1& operator=(const AliITSUv1 &source); // assignment operator

  TGeoVolume* CreateWrapperVolume(const Int_t nLay);

  //
  Int_t     fNWrapVol;       // number of wrapper volumes
  Double_t* fWrapRMin;       // min radius of wrapper volume
  Double_t* fWrapRMax;       // max radius of wrapper volume
  Double_t* fWrapZSpan;      // Z span of wrapper volume
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
  Int_t    *fBuildLevel;     // Vector of Material Budget Studies
  //  
  AliITSUv1Layer **fUpGeom; //! Geometry
  AliITSUModel_t fStaveModelIB; // The stave model for the Inner Barrel
  AliITSUModel_t fStaveModelOB; // The stave model for the Outer Barrel
  
  // Parameters for the Upgrade geometry
  
  ClassDef(AliITSUv1,0)                          
};
 
#endif
