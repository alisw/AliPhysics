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
    AliITSv11GeometrySupport(Bool_t debug):AliITSv11Geometry(debug){};
    virtual ~AliITSv11GeometrySupport(){};
    //
    virtual void SPDCone(TGeoVolume *Moth);
    virtual void SPDThermalSheald(TGeoVolume *Moth);
    virtual void SDDCone(TGeoVolume *Moth);
    virtual void SSDCone(TGeoVolume *Moth);
    virtual void ServicesCableSupport(TGeoVolume *Moth);

  private:
    ClassDef(AliITSv11GeometrySupport,1) // ITS v11 Support geometry
};

#endif
