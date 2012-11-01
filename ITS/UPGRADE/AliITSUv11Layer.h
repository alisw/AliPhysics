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
  $Id: AliITSUv11Layer.h
 */

#include "AliITSv11Geometry.h"
#include <TGeoManager.h>
#include <TGeoCompositeShape.h>
#include <TGeoXtru.h>

class TGeoVolume;

class AliITSUv11Layer : public AliITSv11Geometry {
  public:
    AliITSUv11Layer();
    AliITSUv11Layer(Int_t debug);
    AliITSUv11Layer(Int_t lay, Int_t debug);
    AliITSUv11Layer(Int_t lay, Bool_t turbo, Int_t debug);
    AliITSUv11Layer(const AliITSUv11Layer &source);
    AliITSUv11Layer& operator=(const AliITSUv11Layer &source);
    virtual ~AliITSUv11Layer();
    //
    Bool_t    IsTurbo() {return fIsTurbo;};

    Double_t  GetLadderThick() const {return fLadderThick;};
    Double_t  GetLadderTilt()  const {return fLadderTilt;};
    Double_t  GetLadderWidth() const {return fLadderWidth;};
    Double_t  GetSensorThick() const {return fSensorThick;};
    Double_t  GetNLadders()    const {return fNLadders;};
    Double_t  GetNModules()    const {return fNModules;};
    Double_t  GetRadius()      const {return fLayRadius;};
    Double_t  GetPhi0()        const {return fPhi0;};
    Double_t  GetZLength()     const {return fZLength;};
    Int_t     GetDetType()     const {return fDetTypeID;}
    //
    void      SetLadderThick(Double_t t)    {fLadderThick = t;};
    void      SetLadderTilt(Double_t t);
    void      SetLadderWidth(Double_t w);
    void      SetSensorThick(Double_t t)    {fSensorThick = t;};
    void      SetNLadders(Int_t n)          {fNLadders = n;};
    void      SetNModules(Int_t m)          {fNModules = m;};
    void      SetRadius(Double_t r)         {fLayRadius = r;};
    void      SetPhi0(Double_t phi)         {fPhi0 = phi;}
    void      SetZLength(Double_t z)        {fZLength   = z;};
    void      SetDetType(Int_t tp)          {fDetTypeID = tp;}
    virtual void CreateLayer(TGeoVolume *moth,const TGeoManager *mgr=gGeoManager);

  private:
    void CreateLayerTurbo(TGeoVolume *moth,const TGeoManager *mgr=gGeoManager);

    Double_t RadiusOfTurboContainer();

    TGeoVolume* CreateLadder(const TGeoManager *mgr=gGeoManager);
    TGeoVolume* CreateModule(Double_t x,Double_t y, Double_t z, const TGeoManager *mgr=gGeoManager);

    Int_t     fLayerNumber; // Current layer number
    Double_t  fPhi0;        // lab phi of 1st ladder, in degrees!!!
    Double_t  fLayRadius;   // Inner radius of this layer
    Double_t  fZLength;     // Z length of this layer
    Double_t  fSensorThick; // Sensor thickness
    Double_t  fLadderThick; // Ladder thickness
    Double_t  fLadderWidth; // Ladder width (for turbo layers only)
    Double_t  fLadderTilt;  // Ladder tilt angle (for turbo layers only) in degrees
    Int_t     fNLadders;    // Number of ladders in this layer
    Int_t     fNModules;    // Number of modules per ladder in this layer
    UInt_t    fDetTypeID;   // detector type id
    Bool_t    fIsTurbo;     // True if this layer is a "turbo" layer
  // Parameters for the Upgrade geometry

    static const Double_t fgkDefaultSensorThick; // Default sensor thickness
    static const Double_t fgkDefaultLadderThick; // Default ladder thickness

  ClassDef(AliITSUv11Layer,0) // ITS v11 Upgrade geometry
};

#endif
