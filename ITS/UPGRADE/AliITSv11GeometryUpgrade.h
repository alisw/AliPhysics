#ifndef ALIITSV11GEOMETRYUPGRADE_H
#define ALIITSV11GEOMETRYUPGRADE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


//*************************************************************************
// This class Defines the Geometry for the ITS Upgrade using TGeo
// This is a work class used to study different configurations
// during the development of the new ITS structure.
//
//  Mario Sitta <sitta@to.infn.it>
//*************************************************************************


/*
  $Id: AliITSv11GeometryUpgrade.h
 */

#include "AliITSv11Geometry.h"
#include <TGeoManager.h>
#include <TGeoCompositeShape.h>
#include <TGeoXtru.h>

class TGeoVolume;

class AliITSv11GeometryUpgrade : public AliITSv11Geometry {
  public:
    AliITSv11GeometryUpgrade();
    AliITSv11GeometryUpgrade(Int_t debug);
    AliITSv11GeometryUpgrade(Int_t lay, Int_t debug);
    AliITSv11GeometryUpgrade(Int_t lay, Bool_t turbo, Int_t debug);
    AliITSv11GeometryUpgrade(const AliITSv11GeometryUpgrade &source);
    AliITSv11GeometryUpgrade& operator=(const AliITSv11GeometryUpgrade &source);
    virtual ~AliITSv11GeometryUpgrade();
    //
    Bool_t    IsTurbo() {return fIsTurbo;};

    Double_t  GetLadderThick() {return fLadderThick;};
    Double_t  GetLadderTilt()  {return fLadderTilt;};
    Double_t  GetLadderWidth() {return fLadderWidth;};
    Double_t  GetSensorThick() {return fSensorThick;};
    Double_t  GetNLadders()    {return fNLadders;};
    Double_t  GetNModules()    {return fNModules;};
    Double_t  GetRadius()      {return fLayRadius;};
    Double_t  GetZLength()     {return fZLength;};

    void      SetLadderThick(const Double_t t) {fLadderThick = t;};
    void      SetLadderTilt(const Double_t t);
    void      SetLadderWidth(const Double_t w);
    void      SetSensorThick(const Double_t t) {fSensorThick = t;};
    void      SetNLadders(const Int_t n) {fNLadders = n;};
    void      SetNModules(const Int_t m) {fNModules = m;};
    void      SetRadius(const Double_t r) {fLayRadius = r;};
    void      SetZLength(const Double_t z) {fZLength   = z;};

    virtual void CreateLayer(TGeoVolume *moth,
		       const TGeoManager *mgr=gGeoManager);

  private:
    void CreateLayerTurbo(TGeoVolume *moth,
		    const TGeoManager *mgr=gGeoManager);

    TGeoVolume* CreateLadder(const TGeoManager *mgr=gGeoManager);
    TGeoVolume* CreateModule(const Double_t x, const Double_t y,
			     const Double_t z,
			     const TGeoManager *mgr=gGeoManager);

    Int_t     fLayerNumber; // Current layer number
    Double_t  fLayRadius;   // Inner radius of this layer
    Double_t  fZLength;     // Z length of this layer
    Double_t  fSensorThick; // Sensor thickness
    Double_t  fLadderThick; // Ladder thickness
    Double_t  fLadderWidth; // Ladder width (for turbo layers only)
    Double_t  fLadderTilt;  // Ladder tilt angle (for turbo layers only)
    Int_t     fNLadders;    // Number of ladders in this layer
    Int_t     fNModules;    // Number of modules per ladder in this layer
    Bool_t    fIsTurbo;     // True if this layer is a "turbo" layer

  // Parameters for the Upgrade geometry

    static const Double_t fgkDefaultSensorThick; // Default sensor thickness
    static const Double_t fgkDefaultLadderThick; // Default ladder thickness

  ClassDef(AliITSv11GeometryUpgrade,0) // ITS v11 Upgrade geometry
};

#endif
