#ifndef ALIITSUV2_H
#define ALIITSUV2_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//========================================================================
//
//        Geometry for the Upgrade of the Inner Tracking System
//
// Mario Sitta (sitta@to.infn.it)
//
//========================================================================


// $Id: AliITSUv2.h 

#include "AliITSU.h"

class  AliITSUv2Layer;
//class  AliITSv11GeomBeamPipe;
class  TGeoVolume;
class  TGeoVolumeAssembly;

class AliITSUv2 : public AliITSU {

 public:


  typedef enum {
    kIBModelDummy=0,
    kIBModel0=1,
    kIBModel1=2, 
    kIBModel21=3,
    kIBModel22=4,
    kIBModel3=5,
    kIBModel4=10,
    kOBModelDummy=6,
    kOBModel0=7,
    kOBModel1=8, 
    kOBModel2=9 
  } AliITSUModel_t;
  

  AliITSUv2();
  AliITSUv2(const char *title, Int_t nlay);
  virtual       ~AliITSUv2() ;
  virtual void   SetNWrapVolumes(Int_t n);
  virtual void   AddAlignableVolumes() const;
  void           AddAlignableVolumesLayer(int lr, TString& parent,Int_t &lastUID) const;
  void           AddAlignableVolumesStave(int lr, int st, TString& parent,Int_t &lastUID) const;
  void           AddAlignableVolumesHalfStave(int lr, int st, int sst, TString& parent,Int_t &lastUID) const;
  void           AddAlignableVolumesModule(int lr, int st, int sst, int md, TString& parent,Int_t &lastUID) const;
  void           AddAlignableVolumesChip(int lr, int st, int sst, int md, int ch, TString& parent,Int_t &lastUID) const;

          void   AddGammaConversionRods(const Double_t diam, const Int_t num);
  virtual void   CreateGeometry();
  	  void   CreateSuppCyl(const Bool_t innerBarrel,TGeoVolume *dest,const TGeoManager *mgr=gGeoManager);
  virtual void   CreateMaterials();
  virtual void   DefineLayer(Int_t nlay,Double_t phi0,Double_t r,Double_t zlen,Int_t nstav,
			     Int_t nunit, Double_t lthick=0.,Double_t dthick=0.,UInt_t detType=0, Int_t buildFlag=0);
  virtual void   DefineLayerTurbo(Int_t nlay,Double_t phi0,Double_t r,Double_t zlen,Int_t nstav,
				  Int_t nunit,Double_t width,Double_t tilt,
				  Double_t lthick = 0.,Double_t dthick = 0.,UInt_t detType=0, Int_t buildFlag=0);
  virtual void   GetLayerParameters(Int_t nlay, Double_t &phi0,Double_t &r, Double_t &zlen,
				    Int_t &nstav, Int_t &nmod,
				    Double_t &width, Double_t &tilt,
				    Double_t &lthick, Double_t &mthick,
				    UInt_t &dettype) const;
  virtual void   DefineWrapVolume(Int_t id, Double_t rmin,Double_t rmax, Double_t zspan);
  virtual void   Init(); 
  virtual Bool_t IsLayerTurbo(Int_t nlay);
  virtual Int_t  IsVersion()                 const { return 20;}  // vUpgrade ? do we need this
  virtual void   SetDefaults();
  virtual void   StepManager();
  virtual void   SetLayerChipTypeID(Int_t lr, UInt_t id);
  virtual Int_t  GetLayerChipTypeID(Int_t lr);
  virtual void   SetStaveModelIB(AliITSUModel_t model) {fStaveModelIB=model;}
  virtual void   SetStaveModelOB(AliITSUModel_t model) {fStaveModelOB=model;}
  virtual AliITSUModel_t GetStaveModelIB() const {return fStaveModelIB;}
  virtual AliITSUModel_t GetStaveModelOB() const {return fStaveModelOB;}
  //
 private:
  AliITSUv2(const AliITSUv2 &source); // copy constructor
  AliITSUv2& operator=(const AliITSUv2 &source); // assignment operator

  void        CreateGammaConversionRods(const Double_t diam, const Int_t num,
					const TGeoManager *mgr=gGeoManager);
  TGeoVolume* CreateWrapperVolume(Int_t nLay);

  //
  Int_t     fNWrapVol;       // number of wrapper volumes
  Double_t* fWrapRMin;       // min radius of wrapper volume
  Double_t* fWrapRMax;       // max radius of wrapper volume
  Double_t* fWrapZSpan;      // Z span of wrapper volume
  Int_t*    fLay2WrapV;      // id of wrapper layer to which layer belongs (-1 if not wrapped)
  Bool_t   *fLayTurbo;       // True for "turbo" layers
  Double_t *fLayPhi0;        // Vector of layer's 1st stave phi in lab
  Double_t *fLayRadii;       // Vector of layer radii
  Double_t *fLayZLength;     // Vector of layer length along Z
  Int_t    *fStavPerLay;     // Vector of number of staves per layer
  Int_t    *fUnitPerStave;   // Vector of number of "units" per stave
  Double_t *fChipThick;      // Vector of chip thicknesses
  Double_t *fStaveWidth;     // Vector of stave width (only used for turbo)
  Double_t *fStaveTilt;      // Vector of stave tilt (only used for turbo)
  Double_t *fSensThick;      // Vector of sensor thicknesses
  UInt_t   *fChipTypeID;     // Vector of detector type id
  Int_t    *fBuildLevel;     // Vector of Material Budget Studies
  //  
  AliITSUv2Layer **fUpGeom; //! Geometry
  AliITSUModel_t fStaveModelIB; // The stave model for the Inner Barrel
  AliITSUModel_t fStaveModelOB; // The stave model for the Outer Barrel
  //
  Bool_t    fAddGammaConv;   // True to add Gamma Conversion Rods
  Double_t  fGammaConvDiam;  // Gamma Conversion Rod Diameter
  Int_t     fNRodGammaConv;  // Number of Gamma Conversion Rods

  // Parameters for the Upgrade geometry
  
  ClassDef(AliITSUv2,0)
};
 
#endif
