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


// $Id: AliITSvUpgrade.h 

#include "AliITSUpg.h"
#include "AliITSInitGeometryUpg.h"

class  AliITSv11GeometryUpgrade;
class  TGeoVolume;
class  TGeoVolumeAssembly;

class AliITSvUpgrade : public AliITSUpg {

 public:
    AliITSvUpgrade();
    AliITSvUpgrade(const char *title);
    AliITSvUpgrade(const char *name, const char *title, const Int_t nlay);
    virtual       ~AliITSvUpgrade() ;

    virtual void   AddAlignableVolumes() const;
    virtual void   CreateGeometry();
    virtual void   CreateMaterials();
    virtual void   DefineLayer(const Int_t nlay, const Double_t r,
			       const Double_t zlen, const Int_t nladd,
			       const Int_t nmod, const Double_t lthick = 0.,
			       const Double_t dthick = 0.);
    virtual void   DefineLayerTurbo(const Int_t nlay, const Double_t r,
				    const Double_t zlen, const Int_t nladd,
				    const Int_t nmod, const Double_t width,
				    const Double_t tilt,
				    const Double_t lthick = 0.,
				    const Double_t dthick = 0.);
    virtual void   GetLayerParameters(const Int_t nlay,
				      Double_t &r, Double_t &zlen,
				      Int_t &nladd, Int_t &nmod,
				      Double_t &width, Double_t &tilt,
				      Double_t &lthick, Double_t &mthick);
    virtual Int_t  GetMajorVersion() const {return fMajorVersion;}
    virtual Int_t  GetMinorVersion() const {return fMinorVersion;}
    virtual Int_t  GetNumberOfLayers() const {return fNumberOfLayers;}
    virtual void   Init(); 
    virtual Bool_t IsLayerTurbo(const Int_t nlay);
    virtual Int_t  IsVersion() const {// returns the ITS version number 
                                      return 20;}  // vUpgrade
    virtual void   SetDefaults();
    virtual void   SetMinorVersion(Int_t v=2) {fMinorVersion = v;}
    virtual void   SetNumberOfLayers(Int_t n) {fNumberOfLayers = n;}
    virtual void   StepManager();

 protected:
    void SetT2Lmatrix(Int_t uid, Double_t yShift,
		      Bool_t yFlip, Bool_t yRot180=kFALSE) const; // Set T2L matrix in TGeoPNEntries

 private:
    AliITSvUpgrade(const AliITSvUpgrade &source); // copy constructor
    AliITSvUpgrade& operator=(const AliITSvUpgrade &source); // assignment operator

    Int_t   fMajorVersion;     // Major version number == IsVersion
    Int_t   fMinorVersion;     // Minor version number

    Int_t   fNumberOfLayers;   // Number of layers
    Bool_t   *fLayTurbo;       // True for "turbo" layers
    Double_t *fLayRadii;       // Vector of layer radii
    Double_t *fLayZLength;     // Vector of layer length along Z
    Int_t    *fLaddPerLay;     // Vector of number of ladders per layer
    Int_t    *fModPerLadd;     // Vector of number of modules per ladder
    Double_t *fLadThick;       // Vector of ladder thicknesses
    Double_t *fLadWidth;       // Vector of ladder width (only used for turbo)
    Double_t *fLadTilt;        // Vector of ladder tilt (only used for turbo)
    Double_t *fDetThick;       // Vector of detector thicknesses

    AliITSInitGeometryUpg fInitGeom;   //! Get access to decoding and AliITSgeom init functins

    AliITSv11GeometryUpgrade **fUpGeom; //! Geometry

  // Parameters for the Upgrade geometry

    ClassDef(AliITSvUpgrade,0)                          
};
 
#endif
