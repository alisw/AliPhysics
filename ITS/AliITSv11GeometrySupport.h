#ifndef ALIITSV11GEOMETRYSUPPORT_H
#define ALIITSV11GEOMETRYSUPPORT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*
  $Id$
 */
#include <AliITSv11Geometry.h>
class TGeoVolume;

class AliITSv11GeometrySupport : public AliITSv11Geometry {
  public:
    AliITSv11GeometrySupport(){};
    AliITSv11GeometrySupport(Int_t debug):AliITSv11Geometry(debug){};
    virtual ~AliITSv11GeometrySupport(){};
    //
    virtual void SPDCone(TGeoVolume *moth,TGeoManager *mgr=gGeoManager);
    virtual void SPDThermalSheald(TGeoVolume *moth,
                          TGeoManager *mgr=gGeoManager); // called by SPDCone.
    virtual void SDDCone(TGeoVolume *moth,TGeoManager *mgr=gGeoManager);
    virtual void SSDCone(TGeoVolume *moth,TGeoManager *mgr=gGeoManager);
    virtual void ServicesCableSupport(TGeoVolume *moth,
                                      TGeoManager *mgr=gGeoManager);

  private:
    ClassDef(AliITSv11GeometrySupport,1) // ITS v11 Support geometry
};

#endif
