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
    // some local useful functions
    static Double_t RmaxFromZpCone(const Double_t *ar,const Double_t *az,
                                   Double_t tc,Double_t z,
                                   Double_t th=0.0){
        return AliITSv11Geometry::RFromZpCone(ar,az,4,tc,z,th);};
    static Double_t RminFromZpCone(const Double_t *ar,const Double_t *az,
                                   Double_t tc,Double_t z,Double_t th=0.0){
        return AliITSv11Geometry::RFromZpCone(ar,az,3,tc,z,th);};
    static Double_t ZFromRmaxpCone(const Double_t *ar,const Double_t *az,
                                   Double_t tc,Double_t r,
                                   Double_t th=0.0){
        return AliITSv11Geometry::ZFromRmaxpCone(ar,az,4,tc,r,th);};
    //
    static Double_t RmaxFromZpCone(const TGeoPcon *p,Double_t tc,Double_t z,
                                   Double_t th=0.0){
        return AliITSv11Geometry::RmaxFromZpCone(p,4,tc,z,th);};
    static Double_t RminFromZpCone(const TGeoPcon *p,Double_t tc,Double_t z,
                                   Double_t th=0.0){
        return AliITSv11Geometry::RminFromZpCone(p,3,tc,z,th);};
    static Double_t ZFromRmaxpCone(const TGeoPcon *p,Double_t tc,Double_t r,
                                   Double_t th=0.0)
        {return AliITSv11Geometry::ZFromRmaxpCone(p,4,tc,r,th);};
    static Double_t ZFromRminpCone(const TGeoPcon *p,Double_t tc,Double_t r,
                                   Double_t th=0.0)
        {return AliITSv11Geometry::ZFromRminpCone(p,3,tc,r,th);};

  private:
    ClassDef(AliITSv11GeometrySupport,1) // ITS v11 Support geometry
};

#endif
